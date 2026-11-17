"""Independent physical-coordinate check of the reported conic example.

Reads the saved fit and observations. It does not run a new optimiser or
diagnose numerical global optimality. Output is a reproducible final-point
feasibility and likelihood check.
"""
from pathlib import Path
import json
import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal

ROOT = Path(__file__).resolve().parent
run = json.loads((ROOT / 'operator/full_inner_results.json').read_text())
X = np.loadtxt(ROOT / 'operator/full_inner_data.csv', delimiter=',', skiprows=1)
z = {key: np.asarray(value, dtype=float) for key, value in run['final_lift'].items()}
B = z['B']; v = z['v']; alpha = z['alpha']
mu = np.linalg.solve(B, v.T).T
Bi = np.linalg.inv(B)
C = Bi @ Bi.T
logpi = alpha - logsumexp(alpha)
physical_logpdf = np.column_stack([
    multivariate_normal.logpdf(X, mean=m, cov=C) + logpi[j]
    for j, m in enumerate(mu)
])
physical_ll = float(logsumexp(physical_logpdf, axis=1).mean())
reported_ll = run['steps'][-1]['log_likelihood_per_observation']
assert abs(physical_ll - reported_ll) < 1e-12
assert np.max(np.abs(C - np.array(run['final_shared_covariance']))) < 1e-12

K, p = v.shape
pairs = [(j, h) for j in range(K) for h in range(K) if j != h]
beta = run['settings']['beta']; tau = run['settings']['tau']
nu = run['settings']['nu']
first, second, multipliers, sep, caps = [], [], [], [], []
actualA = {}
for k, (j, h) in enumerate(pairs):
    delta = v[h] - v[j]
    d2 = float(delta @ delta)
    A = np.exp(alpha[h]-alpha[j]-d2/2+z['r'][j]*np.sqrt(d2))
    actualA[j, h] = A
    r = z['r'][j]
    sep.append(d2-(nu+1+tau+r*r+nu/(r*r)))
    caps.append(z['a'][k]-A)
    multipliers.append(np.linalg.eigvalsh(z['lam'][k])[0])
    dc = delta-z['c'][j]
    second.append(np.linalg.eigvalsh(beta*np.eye(p)+z['lam'][k]-np.outer(dc,dc))[0])
for j in range(K):
    weighted = sum(actualA[a,b]*z['lam'][k] for k,(a,b) in enumerate(pairs) if a==j)
    first.append(np.linalg.eigvalsh(beta*np.eye(p)-np.outer(z['c'][j],z['c'][j])-weighted)[0])
checks = {
    'physical_log_likelihood_per_observation': physical_ll,
    'reported_log_likelihood_per_observation': reported_ll,
    'absolute_log_likelihood_discrepancy': abs(physical_ll-reported_ll),
    'minimum_first_matrix_slack': float(min(first)),
    'minimum_second_matrix_slack': float(min(second)),
    'minimum_multiplier_eigenvalue': float(min(multipliers)),
    'minimum_squared_separation_slack': float(min(sep)),
    'minimum_odds_epigraph_slack': float(min(caps)),
    'method': 'Independent physical-coordinate Gaussian log densities and direct final matrix inequalities',
    'new_optimisation_runs': 0,
}
assert min(first+second+multipliers+sep+caps)>0
(ROOT / 'final_fit_validation.json').write_text(json.dumps(checks, indent=2)+'\n')
print(json.dumps(checks, indent=2))
