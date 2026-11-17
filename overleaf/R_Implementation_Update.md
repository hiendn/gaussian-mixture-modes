# R computation update, 17 November 2026

The numerical examples and graphics are now implemented in native R using CVXR and Clarabel. The current manuscript reports the executed R fitting results and links to the public repository at https://github.com/hiendn/gaussian-mixture-modes.

The R port retains the original observations, parameter envelopes, certificate settings, proximal norms and distinct fixed-budget and variable-budget supports. The three fixed-budget likelihood gains agree with the archived Python results within 7e-8. The joint-budget run differs by about 1.7e-6 in likelihood gain; its witness trajectory and proximal penalty also differ. The current manuscript uses the actual R numbers rather than implying identical iterates across solver interfaces.

All four accepted R states passed independent physical-coordinate likelihood, nonlinear feasibility, convex-inner feasibility, whole-state penalty and ascent checks. Every joint-budget witness is now saved, removing the original saved-record limitation described in the first-draft author note. The new records do not retroactively supply the missing witnesses from the old Python run.

The exact geometry module performs 31 strict rational sign checks. The fitted-state checks remain floating-point calculations. The joint-budget solve reports optimal_inaccurate, and the R interface did not expose a primal-dual gap. These bounded checks do not establish an asymptotic inner-gap condition, convergence of a longer run, or global likelihood optimality.

The geometry, conditional algorithmic convergence, statistical target and scientific qualifications in the theory sections are unchanged. Comparative experiments and a substantive application remain future work. The first-draft author note remains a historical record of the initial version.
