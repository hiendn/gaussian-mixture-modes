# Likelihood estimation of Gaussian mixtures with one mode per component

Working manuscript with native R computations, updated 17 November 2026.

The main document is `manuscript.tex`. Compile with XeLaTeX. The bibliography is included directly from `references.tex`, so BibTeX is unnecessary. `R_Updated_Manuscript.pdf` is the locally compiled R update; `First_Draft.pdf` retains the initial draft.

The public manuscript and computation repository is:

https://github.com/hiendn/gaussian-mixture-modes

The R implementation is in `R/`; `scripts/run_all.R` generates the numerical records, figures and generated LaTeX table and values. `scripts/bootstrap.R` installs missing R dependencies into a project-local library. From a downloaded copy of this project, run:

    Rscript scripts/bootstrap.R
    Rscript scripts/run_all.R
    Rscript scripts/verify.R

The scripts support this Overleaf folder layout and the GitHub layout. R computations must be run locally; Overleaf compiles LaTeX and does not execute the R scripts. Recompile the LaTeX document after replacing generated values or figures.

The R computations reproduce the three geometric examples, triangle contour, three fixed-budget steps, and one joint-budget step. The numerical section reports the native R results. `results/` contains all anchor, trial and accepted witnesses, exact and inner residuals, software versions and solver statuses. The joint-budget status remains `optimal_inaccurate`; no numerical primal-dual gap was exposed by the R interface.

The original `code/` directory is preserved. A second archival copy under `legacy/python/` supports cross-language comparisons. `Author_Note.md` and `Author_Note.pdf` describe the first draft. `R_Implementation_Update.md` records the subsequent computational changes. The comparative simulation study and substantive application remain outstanding.
