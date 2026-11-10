# Notebooks

## Main research notebook

**[UDG_Combined_Notebook.ipynb](UDG_Combined_Notebook.ipynb)** — supplied as `btp-1-ultra-diffuse-galaxies.ipynb` and renamed for repository clarity. Includes the DF2/DF4 and DF44 stellar-only Jeans calculations, NFW mass modelling, asymmetric-error fitting with fixed/free concentration and free stellar mass-to-light ratio, diagnostics, saved plots and LOS convergence checks.

Only portability changes have been made: it can find the three CSVs from within `notebooks/` using `../data/`, and writes output plots into `../figures/` and output files into `../results/`. The underlying analysis equations and fitting choices were not changed.

## Guided walkthrough

**[UDG_Analysis_Walkthrough.ipynb](UDG_Analysis_Walkthrough.ipynb)** — shorter, reproducible walkthrough using `src/udg_jeans.py` and `src/run_analysis.py`.

## Running

From the repository root:

```bash
python -m pip install -r requirements.txt
jupyter lab
```

Open `notebooks/UDG_Combined_Notebook.ipynb` and run all cells in order. Some cells optimize multiple scenarios and will take longer than a simple demonstration notebook.
