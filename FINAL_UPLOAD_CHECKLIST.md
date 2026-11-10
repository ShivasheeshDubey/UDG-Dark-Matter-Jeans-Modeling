# Final upload checklist

- [x] Final user-supplied presentation retained in `presentation/BTP_1_Presentation.pptx`.
- [x] Final LaTeX source retained in `report/main.tex`.
- [x] Public PDF derived from LaTeX retained without faculty signatures and the similarity appendix.
- [x] Reusable numerical solver and executable analysis script retained in `src/`.
- [x] Short walkthrough notebook retained in `notebooks/`.
- [x] User-supplied notebook `btp-1-ultra-diffuse-galaxies.ipynb` included as `notebooks/UDG_Combined_Notebook.ipynb`, with only portable path adaptations.
- [x] No PowerPoint creation/source scripts included in `presentation/`.
- [x] Observation CSVs, derived plots, and fit summaries included.
- [ ] Review third-party dataset and image redistribution rights before making the repository public.
- [ ] Choose a code license deliberately (none assumed by default).
- [ ] Double-check personal information in report/presentation against desired public visibility.

## Publish the folder to GitHub

After creating an empty GitHub repository named `udg-dark-matter-jeans-modeling`, extract the ZIP, `cd` into `UDG_Dark_Matter_GitHub_UPLOAD_READY`, and run:

```bash
git init
git branch -M main
git add .
git commit -m "Publish B.Tech UDG Jeans modelling project"
git remote add origin https://github.com/YOUR_USERNAME/udg-dark-matter-jeans-modeling.git
git push -u origin main
```

Replace YOUR_USERNAME with your GitHub username. This ZIP is not already hosted on GitHub; these commands upload it.
