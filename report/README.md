# Project report

- [BTP_Report.pdf](BTP_Report.pdf) — public copy of the final B.Tech report, rebuilt from the supplied `main.tex`.
- [main.tex](main.tex) — editable LaTeX source supplied with the final report.

**Public-release choice:** The uploaded signed submission PDF was not copied directly into this GitHub package. It has supervisor signatures on its certificate page and three appended similarity/originality-report pages. The PDF here was compiled from the supplied LaTeX source, so it contains the 26-page academic report **without the handwritten signatures or the appended similarity report**. Please share the signed submission privately rather than publishing faculty signatures without permission.

From the repository root, to rebuild the report (requires a working TeX Live installation):

```bash
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=report report/main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=report report/main.tex
```

The source expects the plots in the top-level `figures/` folder. The optional institute logo may be absent; the LaTeX source handles that case.
