# Numerical evidence and provenance

The manuscript's original package was read from its existing Overleaf-linked Dropbox folder on 17 November 2026. `verification/original_sha256.json` records the source-file hashes. The original computational scripts, data and saved outputs are preserved under `legacy/python/`; the local `original/` backup also retains the initial manuscript and author notes and is excluded from Git.

The public manuscript reports a newly executed native R implementation. Geometry settings, observations, parameter envelopes, likelihood scaling, proximal norms and inner supports were retained. The fixed-budget and variable-budget separation supports remain distinct. The proximal norm includes all state coordinates, with both off-diagonal entries of symmetric matrices counted in the Frobenius norm.

The original contour coordinates and fitting initialisation represent different rotations and translations of the same triangle. The graphics use the displayed triangle with vertices `(0,0)`, `(2.8,0)` and `(1.4,sqrt(3)*1.4)`. The fit starts from the centred triangle in the original algorithm.

`results/exact_geometry.json` retains exact rational interval endpoints for all 31 strict geometric sign checks. Decimal midpoints are display values. The local-budget residuals near `1e-72` remain positive and are retained in `results/geometry.json`.

`results/fixed_budget.json` and `results/joint_budget.json` retain every full anchor, trial and accepted witness state, exact residuals, inner checks, acceptance fractions, solver status, objective changes and settings. Their RDS companions preserve R matrix and array dimensions directly. The joint record closes the original saved Python record's missing-witness limitation. The Python joint run cannot be fully reconstructed from its saved JSON alone.

`verification/independent_fit_review.json` and its Markdown companion record independent calculations for the saved accepted states. `scripts/verify.R` rechecks the geometric signs, physical-coordinate Gaussian likelihoods, full-state penalties, reported ascent inequalities, agreement with the archived examples and generated manuscript table. It also runs `verification/independent_fit_check.R`, which independently recomputes all nonlinear and inner constraints without sourcing production fitting helpers. These scripts use the saved states; they do not run an optimiser.

To generate a new trajectory, run `Rscript scripts/run_all.R`. This replaces the derived results, figures and generated LaTeX values. Preserve a commit or copy before replacing a result set. Solver versions and canonicalisation can change the retained witness trajectory, particularly when a solution is marked `optimal_inaccurate`. Review new residuals and software statuses before using new results in a paper.

The statistical convergence theorem and the algorithmic convergence theorem are mathematical results under their stated assumptions. The numerical illustrations do not independently establish those assumptions or diagnose a global empirical objective gap.
