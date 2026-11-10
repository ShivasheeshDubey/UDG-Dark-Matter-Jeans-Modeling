"""Newtonian spherical Jeans models for DF2, DF4 and Dragonfly 44.

Physical inputs and model choices follow the accompanying B.Tech report.  This
module is an independent, minimalist implementation of the *Newtonian* models,
not a reproduction of the source studies' full likelihood machinery.

Units throughout: distance kpc; mass solar masses; velocity km/s.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import least_squares
from scipy.special import gamma

G_KPC_KMS2_MSUN = 4.30091e-6
H0_KMS_PER_KPC = 70.0 / 1000.0
RHO_CRITICAL = 3 * H0_KMS_PER_KPC**2 / (8 * np.pi * G_KPC_KMS2_MSUN)


@dataclass(frozen=True)
class JeansGrid:
    """Numerical grid for enclosing mass, Jeans pressure and LOS projection."""
    r: np.ndarray
    z: np.ndarray

    @classmethod
    def standard(cls):
        return cls(
            r=np.geomspace(1e-4, 200.0, 8000),
            z=np.r_[0.0, np.geomspace(1e-4, 180.0, 2600)],
        )


@dataclass
class TracerModel:
    tracer: np.ndarray
    stellar_mass: np.ndarray
    grid: JeansGrid


def enclosed_mass(r, density):
    """M(<r)=4 pi integral_0^r density(s) s^2 ds (lower cutoff ~=0)."""
    return 4 * np.pi * cumulative_trapezoid(density * r**2, r, initial=0)


def df2_df4_density(r, rho0, a_star):
    """Compact stellar density model from Bhatia et al. 2023.

    rho0: Msun/kpc^2; a_star: kpc; returns Msun/kpc^3.
    """
    return (40 * rho0 / (63 * a_star)) * np.exp(-((11 * r / (10 * a_star)) ** (4 / 3)))


def df44_stellar_tracer(r, *, n=0.94, re=4.7, luminosity=2.33e8):
    """Prugniel-Simien-type approximation to deprojected Sersic luminosity."""
    p = 1 - 0.6097 / n + 0.05463 / n**2
    b = 2 * n - 1 / 3 + 0.009876 / n
    a_s = re / b**n
    j0 = luminosity / (4 * np.pi * n * a_s**3 * gamma(n * (3 - p)))
    luminosity_density = j0 * (r / a_s) ** (-p) * np.exp(-(r / a_s) ** (1 / n))
    return luminosity_density


def make_tracer(grid, density, mass_to_light=1.0):
    density = np.asarray(density, dtype=float)
    return TracerModel(density, enclosed_mass(grid.r, density * mass_to_light), grid)


def nfw_mass(r, m200, c):
    """NFW *halo* mass within r, normalized to M200 within r200."""
    r200 = (3 * m200 / (4 * np.pi * 200 * RHO_CRITICAL)) ** (1 / 3)
    r_s = r200 / c
    f = lambda x: np.log1p(x) - x / (1 + x)
    return m200 * f(r / r_s) / f(c)


def nfw_scale_radius(m200, c):
    r200 = (3 * m200 / (4 * np.pi * 200 * RHO_CRITICAL)) ** (1 / 3)
    return r200, r200 / c


def radial_variance(tracer, total_enclosed_mass, beta):
    """Constant-anisotropy Jeans solution, with zero outer tracer pressure."""
    r = tracer.grid.r
    nu = tracer.tracer
    weight = r ** (2 * beta)
    integrand = nu * G_KPC_KMS2_MSUN * total_enclosed_mass / r**2 * weight
    # Backward integration from r_max toward r; reverse-axis trapezoid is negative.
    integrated = -cumulative_trapezoid(integrand[::-1], r[::-1], initial=0)[::-1]
    denom = nu * weight
    return np.divide(integrated, denom, out=np.zeros_like(integrated), where=denom > 1e-270)


def los_dispersion(projected_radii, tracer, total_enclosed_mass, beta):
    """LOS dispersion in km/s, integrating over the sight-line z coordinate."""
    radii = np.atleast_1d(np.asarray(projected_radii, dtype=float))
    if np.any(radii <= 0):
        raise ValueError('Projected radii must be >0 in this implementation.')
    r = tracer.grid.r
    z = tracer.grid.z
    radius = np.sqrt(radii[:, None] ** 2 + z[None, :] ** 2)
    nu_z = np.interp(radius.ravel(), r, tracer.tracer, left=tracer.tracer[0], right=0).reshape(radius.shape)
    vr2 = radial_variance(tracer, total_enclosed_mass, beta)
    vr2_z = np.interp(radius.ravel(), r, vr2, left=vr2[0], right=0).reshape(radius.shape)
    weight = 1 - beta * (radii[:, None] ** 2) / (radius**2)
    numerator = trapezoid(weight * nu_z * vr2_z, z, axis=1)
    denominator = trapezoid(nu_z, z, axis=1)
    return np.sqrt(np.maximum(numerator / denominator, 0))


def df44_fit(data, tracer, *, free_concentration=False):
    """Fit M200, beta and optionally c to 9 independent DF44 dispersions.

    Asymmetric errors: upper error when prediction lies above observation,
    lower error otherwise. Boundary c in [2,40] is a chosen search constraint.
    """
    R = data['r_kpc'].to_numpy(dtype=float)
    obs = data['sigma_kms'].to_numpy(dtype=float)
    upper = data['sigma_err_high'].to_numpy(dtype=float)
    lower = data['sigma_err_low'].to_numpy(dtype=float)

    def residuals(p):
        m200 = 10 ** p[0]
        beta = p[1]
        c = p[2] if free_concentration else 10.0
        total_mass = tracer.stellar_mass + nfw_mass(tracer.grid.r, m200, c)
        pred = los_dispersion(R, tracer, total_mass, beta)
        errors = np.where(pred > obs, upper, lower)
        return (pred - obs) / errors

    if free_concentration:
        start = [np.log10(1.2e12), -0.3, 2.5]
        lower_bounds = [9.0, -2.0, 2.0]
        upper_bounds = [14.0, 0.2, 40.0]
    else:
        start = [np.log10(3.5e10), -0.6]
        lower_bounds = [9.0, -2.0]
        upper_bounds = [14.0, 0.2]

    fit = least_squares(
        residuals, start, bounds=(lower_bounds, upper_bounds),
        ftol=1e-11, xtol=1e-11, gtol=1e-11, max_nfev=200,
    )
    m200 = float(10 ** fit.x[0])
    beta = float(fit.x[1])
    c = float(fit.x[2]) if free_concentration else 10.0
    chi2 = float(np.sum(fit.fun**2))
    r200, rs = nfw_scale_radius(m200, c)
    mass_half = float(nfw_mass(np.array([6.3]), m200, c)[0])
    return {
        'fit': 'free_c' if free_concentration else 'fixed_c10',
        'm200_msun': m200,
        'c': c,
        'beta': beta,
        'chi2': chi2,
        'dof_nominal': int(len(R) - len(fit.x)),
        'chi2_reduced_nominal': chi2 / (len(R) - len(fit.x)),
        'r200_kpc': float(r200),
        'rs_kpc': float(rs),
        'mhalo_at_6p3_kpc_msun': mass_half,
        'on_concentration_bound': bool(free_concentration and abs(c - 2) < 1e-5),
    }


def df44_model_dispersion(R, tracer, *, m200, c, beta):
    total_mass = tracer.stellar_mass + nfw_mass(tracer.grid.r, m200, c)
    return los_dispersion(R, tracer, total_mass, beta)


def load_data(data_dir):
    data_dir = Path(data_dir)
    return (
        pd.read_csv(data_dir / 'DF2_original_GC_2018.csv'),
        pd.read_csv(data_dir / 'DF4_original_GC_2019.csv'),
        pd.read_csv(data_dir / 'DF44_original_dispersion_2019.csv'),
    )
