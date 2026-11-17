# Independent review of the native R fitting checks

All four accepted iterates pass independently recomputed physical-coordinate likelihood, exact nonlinear feasibility, convex-inner feasibility, acceptance blending, and feasible-ascent checks. No production fitting functions were sourced for the execution checks.

| Run | Physical mean log likelihood | Log-likelihood gain | Auxiliary gain | Whole-state penalty | Minimum exact slack | Minimum inner slack |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| fixed 1 | -3.523072257001 | 0.0532493553 | 0.0525815983 | 0.0153370391 | 1.02e-05 | 1.964e-07 |
| fixed 2 | -3.505922962289 | 0.0171492947 | 0.0165042722 | 0.0056807279 | 4.971e-08 | 2.09e-06 |
| fixed 3 | -3.497818663365 | 0.0081042989 | 0.0076395681 | 0.0027068760 | 3.041e-08 | 1.126e-06 |
| joint 1 | -2.791449254922 | 0.1262244495 | 0.1232831668 | 0.0047911796 | 4.781e-06 | 6.362e-05 |

The original 96-observation data file is identical to the R input. The twelve-observation example has 12 rows and two coordinates. Covariances were recomputed as the inverse of B-transpose times B; component means were independently recovered by solving B mu = v. The maximum discrepancy between independently recomputed and stored log likelihoods, gains, and penalties is below 1e-12.

The inspected R implementation preserves the original fixed-budget support, the variable-budget Young support, matrix-product majorant, local-budget curvature bound, ordered-pair indexing, full Frobenius proximal norm, data, envelopes, and retreat fractions. Both the trial and accepted full states are saved. The joint run now retains all density coordinates and all witnesses, permitting direct audits that were impossible for the old partial JSON record.

The R joint likelihood gain differs from the saved Python gain by -1.714e-06; the whole-state proximal penalty differs by 0.0001197. Both runs report optimal_inaccurate. These are separately computed feasible ascent steps, and should not be described as bitwise or iterate-identical reproductions.

The results are floating-point implementation checks. They do not certify global optimality, an asymptotic inner-gap condition, or convergence of a longer trajectory. A numerical primal-dual gap was not independently extracted.

Geometry and figure checks completed separately before reading any exact-geometry worker output: reported values agree with the manuscript, a direct nonidentity-covariance Gaussian check differed by at most 4.44e-16, and independent contour-location evaluations differed by at most 2.78e-17. Triangle means, mode markers, certified circles, and axes were inspected visually.
