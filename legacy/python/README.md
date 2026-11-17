# Numerical material for the first manuscript

This directory preserves the existing numerical evidence and adds a figure-generation script. The saved numerical records are not outputs of a new simulation study. The manuscript reports each check at its actual scope.

## Contents

- make_figures.py: evaluates the stated triangle density and generates figures/triangle_modes.pdf, a PNG preview, and the plotted mode coordinates.
- geometry/check_geometry.py: high-precision evaluations of the local-radius and local-budget examples.
- geometry/rigorous_geometry_checks.py: exact rational sign decisions with exponential enclosures for the strict geometric comparisons.
- operator/check_operator.py: interval triangle witness checks and the earlier bounded operator calculations.
- operator/full_inner_conic.py: the three-step, 96-observation fixed-budget fitting check. It fixes \(M=I\) and \(\nu=q\), uses the fixed-budget separation support, and frees all density coordinates.
- operator/full_inner_data.csv: the 96 observations used in that saved check.
- operator/full_inner_results.json: its full recorded results and final witnesses.
- joint_budget/check_joint_budget_step.py: the separate twelve-observation joint-budget one-step check. The deterministic data are constructed in the script.
- joint_budget/joint_budget_step_check.json: its saved output. Solver status is “optimal_inaccurate”; this record does not retain every final witness.
- validate_final_fit.py and final_fit_validation.json: direct physical-coordinate likelihood and final-point certificate validation for the saved 96-observation fit.
- source_manifest.json: provenance and hashes of copied evidence files.

The original numerical scripts are copied without mathematical changes. They write their results beside themselves. Running them will replace their local JSON outputs, so keep the supplied archive or copy the code to another directory before a rerun if the original recorded results are needed.

## Dependencies and commands

The saved conic computations used Python 3.12, NumPy 2.5.5, SciPy 1.18.1, CVXPY 1.9.3, and Clarabel 0.11.1. The geometric checks also use mpmath. These are recorded environment versions, not a claim that other versions reproduce every solver status or final digit.

From the manuscript directory, the main commands are:

    python code/geometry/rigorous_geometry_checks.py
    python code/geometry/check_geometry.py
    python code/operator/check_operator.py
    python code/operator/full_inner_conic.py
    python code/joint_budget/check_joint_budget_step.py
    python code/validate_final_fit.py
    python code/make_figures.py

The figure was newly generated using NumPy 2.3.5 and Matplotlib 3.10.8. Its posterior-mean residuals were below \(2\times10^{-15}\). The numerical root search is used only to place the displayed mode markers; the theorem and witnesses establish that the mode count is exhaustive.

Exact-arithmetic geometry checks do not certify floating-point fitted iterates. The small fitting checks do not establish comparative accuracy, high-dimensional scaling, the asymptotic inner-gap condition, or global optimality.
