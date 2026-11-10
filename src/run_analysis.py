"""Generate the numerical results and figures used in the UDG B.Tech project.

Run from repository root:
    python -m src.run_analysis

The script uses the included observational CSVs and does not download data.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .udg_jeans import (
    JeansGrid, make_tracer, df2_df4_density, df44_stellar_tracer,
    los_dispersion, nfw_mass, df44_fit, df44_model_dispersion, load_data,
    enclosed_mass,
)

BASE = Path(__file__).resolve().parents[1]
DATA_DIR = BASE / 'data'
FIG_DIR = BASE / 'figures'
RESULT_DIR = BASE / 'results'


def configure_plots():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'savefig.dpi': 200,
        'figure.dpi': 140,
    })


def plot_df2_df4(R, df2, df4, c2, c4):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), sharey=True)
    for ax, title, obs, curve, color in zip(
        axes, ['DF2', 'DF4'], [df2, df4], [c2, c4], ['#0378a6', '#09866f']
    ):
        ax.plot(R, curve, color=color, lw=2.3, label='Stellar-only Jeans prediction')
        err = (obs.v_err_low.to_numpy() + obs.v_err_high.to_numpy()) / 2
        ax.errorbar(
            obs.r_kpc, obs.delta_v_kms, yerr=err, fmt='o',
            markersize=5, color='#bf3c53', ecolor='#9c8b8e',
            capsize=2, label='Individual GC |velocity offset|', zorder=4,
        )
        ax.set(xlim=(0, 8), xlabel='Projected radius R (kpc)', title=title)
        ax.grid(alpha=0.15)
        ax.legend(fontsize=8, loc='upper right')
    axes[0].set_ylabel('Speed / LOS dispersion (km/s)')
    axes[0].set_ylim(0, 35)
    fig.suptitle('DF2 and DF4: qualitative comparison of stellar-only predictions with GC offsets', fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'df2_df4_newtonian.png', bbox_inches='tight')
    plt.close(fig)


def plot_df44_stars(R, obs, curves):
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for name, y, color, ls in [
        ('Isotropic beta=0', curves['isotropic'], '#167eae', '-'),
        ('Tangential beta=-0.6', curves['tangential'], '#0b8c74', '--'),
        ('Radial beta=+0.6 (formal only)', curves['radial_formal'], '#ee983f', ':'),
    ]:
        ax.plot(R, y, color=color, linestyle=ls, linewidth=2.0, label=name)
    ax.errorbar(obs.r_kpc, obs.sigma_kms,
                yerr=[obs.sigma_err_low, obs.sigma_err_high],
                fmt='o', color='#b9304c', ecolor='#b9304c', capsize=2,
                label='DF44 observed stellar dispersions')
    ax.set(xlabel='Projected radius R (kpc)', ylabel='LOS velocity dispersion (km/s)',
           xlim=(0, 6), ylim=(0, 50), title='DF44: stellar-only Newtonian models')
    ax.grid(alpha=0.15)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'df44_stellar_only.png', bbox_inches='tight')
    plt.close(fig)


def plot_df44_nfw(R, obs, curves):
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for name, y, color, ls in [
        ('Initial illustrative NFW', curves['initial'], '#e3a145', ':'),
        ('Best fit: c=10 fixed', curves['fixed_c10'], '#098d85', '-'),
        ('Best fit: c free, at bound 2', curves['free_c'], '#226bbd', '--'),
    ]:
        ax.plot(R, y, color=color, linestyle=ls, linewidth=2.0, label=name)
    ax.errorbar(obs.r_kpc, obs.sigma_kms,
                yerr=[obs.sigma_err_low, obs.sigma_err_high],
                fmt='o', color='#b9304c', ecolor='#b9304c', capsize=2,
                label='Observed stellar dispersions')
    ax.set(xlabel='Projected radius R (kpc)', ylabel='LOS velocity dispersion (km/s)',
           xlim=(0, 6), ylim=(0, 50), title='DF44: Newtonian stellar + NFW halo fits')
    ax.grid(alpha=0.15)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'df44_nfw_bestfit.png', bbox_inches='tight')
    plt.close(fig)


def plot_fit_residuals(R, obs, curves, fits):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6.4), sharex=True,
                                  gridspec_kw={'height_ratios': [1.8, 1.0]})
    for name, y, col in [
        ('Fixed c=10', curves['fixed_c10'], '#098d85'),
        ('Free c, at bound 2', curves['free_c'], '#226bbd'),
    ]:
        ax1.plot(R, y, color=col, lw=2, label=name)
    ax1.errorbar(obs.r_kpc, obs.sigma_kms,
                 yerr=[obs.sigma_err_low, obs.sigma_err_high],
                 fmt='o', color='#b9304c', capsize=2, label='Observed')
    ax1.set(ylabel='LOS velocity dispersion (km/s)', title='DF44 NFW fit residuals')
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.15)
    y_obs = obs.sigma_kms.to_numpy(dtype=float)
    for name, values, marker, col in [
        ('Fixed c=10', fits['pred_fixed'], 'o', '#098d85'),
        ('Free c=2', fits['pred_free'], '^', '#226bbd'),
    ]:
        errs = np.where(values > y_obs, obs.sigma_err_high.to_numpy(),
                        obs.sigma_err_low.to_numpy())
        ax2.plot(obs.r_kpc, (values - y_obs) / errs, marker=marker,
                 ls='none', color=col, label=name)
    ax2.axhline(0, color='#444444', lw=1)
    ax2.set(xlabel='Projected radius R (kpc)',
            ylabel='Standardized residual')
    ax2.grid(alpha=0.15)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'df44_nfw_fit_residuals.png', bbox_inches='tight')
    plt.close(fig)


def plot_halo_mass(fixed, free):
    r = np.geomspace(0.05, 250, 500)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(r, nfw_mass(r, fixed['m200_msun'], fixed['c']), lw=2, color='#098d85', label='Fixed c=10')
    ax.plot(r, nfw_mass(r, free['m200_msun'], free['c']), lw=2, color='#226bbd', ls='--', label='Free c=2 boundary')
    ax.axvline(6.3, color='#777777', ls=':', label='Approx. half-light radius 6.3 kpc')
    ax.set(xscale='log', yscale='log', xlabel='Three-dimensional radius r (kpc)',
           ylabel='Enclosed halo mass M(<r) (solar masses)',
           title='DF44: similar inner mass, very different extrapolated M200')
    ax.grid(alpha=0.18, which='both')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'df44_enclosed_halo_mass.png', bbox_inches='tight')
    plt.close(fig)


def main():
    configure_plots()
    FIG_DIR.mkdir(exist_ok=True)
    RESULT_DIR.mkdir(exist_ok=True)
    df2, df4, df44 = load_data(DATA_DIR)
    grid = JeansGrid.standard()

    # DF2 and DF4 globular-cluster positions and systemic velocity conventions.
    ra0, dec0, distance_mpc = 40.44542, -8.40333, 20.0
    deg_to_kpc = distance_mpc * 1000 * np.pi / 180
    dx = (df2.ra_deg - ra0) * np.cos(np.deg2rad(dec0))
    dy = df2.dec_deg - dec0
    df2['r_kpc'] = np.sqrt(dx**2 + dy**2) * deg_to_kpc
    df2['delta_v_kms'] = np.abs(df2.v_kms - 1803.0)
    df4['delta_v_kms'] = np.abs(df4.v_kms - 1445.0)

    rho_df2 = df2_df4_density(grid.r, rho0=1.25e7, a_star=2.0)
    rho_df4 = df2_df4_density(grid.r, rho0=1.15e7, a_star=1.0)
    t2, t4 = make_tracer(grid, rho_df2), make_tracer(grid, rho_df4)
    stellar_rows = []

    R24 = np.linspace(.1, 8.0, 140)
    curve2 = los_dispersion(R24, t2, t2.stellar_mass, 0)
    curve4 = los_dispersion(R24, t4, t4.stellar_mass, 0)
    plot_df2_df4(R24, df2, df4, curve2, curve4)
    for name, tracer in [('DF2', t2), ('DF4', t4)]:
        disp_1 = los_dispersion([1], tracer, tracer.stellar_mass, 0)[0]
        stellar_rows.append({'galaxy': name, 'beta': 0, 'stellar_mass_msun': tracer.stellar_mass[-1],
                             'sigma_los_at_1_kpc_kms': disp_1,
                             'note': 'stellar-only, isotropic'})

    # DF44: fixed stellar mass-to-light ratio and adopted deprojected Sersic tracer.
    nu_df44 = df44_stellar_tracer(grid.r)
    t44 = make_tracer(grid, nu_df44, mass_to_light=1.0)
    R44 = np.linspace(.15, 6.0, 125)
    stars = {
        'isotropic': los_dispersion(R44, t44, t44.stellar_mass, 0),
        'tangential': los_dispersion(R44, t44, t44.stellar_mass, -0.6),
        'radial_formal': los_dispersion(R44, t44, t44.stellar_mass, 0.6),
    }
    plot_df44_stars(R44, df44, stars)
    for beta, name in [(0.0, 'isotropic'), (-0.6, 'tangential'), (.6, 'radial_formal')]:
        stellar_rows.append({'galaxy': 'DF44', 'beta': beta,
                             'stellar_mass_msun': t44.stellar_mass[-1],
                             'sigma_los_at_1_kpc_kms': los_dispersion([1], t44, t44.stellar_mass, beta)[0],
                             'note': ('formal Jeans-only, fails central cusp-anisotropy condition' if beta == .6
                                      else 'stellar-only, Upsilon*=1')})
    pd.DataFrame(stellar_rows).to_csv(RESULT_DIR / 'stellar_only_summary.csv', index=False)

    fixed = df44_fit(df44, t44, free_concentration=False)
    free = df44_fit(df44, t44, free_concentration=True)
    pd.DataFrame([fixed, free]).to_csv(RESULT_DIR / 'df44_fit_summary.csv', index=False)
    curves = {
        'initial': df44_model_dispersion(R44, t44, m200=7e10, c=10, beta=-0.8),
        'fixed_c10': df44_model_dispersion(R44, t44, m200=fixed['m200_msun'], c=10, beta=fixed['beta']),
        'free_c': df44_model_dispersion(R44, t44, m200=free['m200_msun'], c=free['c'], beta=free['beta']),
    }
    plot_df44_nfw(R44, df44, curves)
    pred_fixed = df44_model_dispersion(df44.r_kpc, t44, m200=fixed['m200_msun'], c=10, beta=fixed['beta'])
    pred_free = df44_model_dispersion(df44.r_kpc, t44, m200=free['m200_msun'], c=free['c'], beta=free['beta'])
    pred_initial = df44_model_dispersion(df44.r_kpc, t44, m200=7e10, c=10, beta=-.8)
    df_pred = df44.copy()
    df_pred['sigma_initial_kms'] = pred_initial
    df_pred['sigma_fixed_c10_kms'] = pred_fixed
    df_pred['sigma_free_c_kms'] = pred_free
    df_pred.to_csv(RESULT_DIR / 'df44_fit_point_predictions.csv', index=False)
    plot_fit_residuals(R44, df44, curves, {'pred_fixed': pred_fixed, 'pred_free': pred_free})
    plot_halo_mass(fixed, free)

    print('Stellar-only models:')
    print(pd.DataFrame(stellar_rows).to_string(index=False, formatters={
        'stellar_mass_msun': '{:.5e}'.format,
        'sigma_los_at_1_kpc_kms': '{:.3f}'.format,
    }))
    print('\nDF44 NFW fits:')
    for f in [fixed, free]:
        print(f"  {f['fit']}: M200={f['m200_msun']:.5e} Msun, "
              f"c={f['c']:.6f}, beta={f['beta']:.6f}, "
              f"chi2={f['chi2']:.4f}, nominal dof={f['dof_nominal']}, "
              f"Mhalo(<6.3 kpc)={f['mhalo_at_6p3_kpc_msun']:.5e} Msun")
    print('\nSaved plots to figures/ and tables to results/.')


if __name__ == '__main__':
    main()
