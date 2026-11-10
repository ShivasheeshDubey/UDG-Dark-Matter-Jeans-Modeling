# Observational data provenance

This folder holds the compact CSV datasets used in the project's computational comparison.

| File | Columns | Role |
|---|---|---|
| `DF2_original_GC_2018.csv` | GC ID, RA/Dec, heliocentric LOS velocity, asymmetric errors | Individual DF2 globular-cluster velocity offsets; projected radii computed in `src/run_analysis.py` using the adopted 20 Mpc distance and centre. |
| `DF4_original_GC_2019.csv` | GC ID, projected radius, LOS velocity, asymmetric errors | Individual DF4 GC velocity offsets. |
| `DF44_original_dispersion_2019.csv` | Projected radius, stellar LOS dispersion, asymmetric errors | Nine stellar-dispersion bins used for the DF44 nonlinear fits. |

The source observations originate from the studies cited in the main README, particularly van Dokkum et al. (2018, 2019). These are **compiled/transcribed research inputs**, not observations made by this project. Check the original publication and any associated machine-readable tables for provenance and reuse rights before public redistribution.

**Important:** The GC data represent *individual tracer velocities*, whereas DF44 data represent *binned stellar velocity dispersions*. The DF2/DF4 figures overlay absolute offsets on a model dispersion curve for qualitative context; they should not be interpreted as a statistically fitted dispersion profile.
