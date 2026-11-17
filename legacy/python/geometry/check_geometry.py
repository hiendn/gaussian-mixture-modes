"""Deterministic checks for the local-coverage research lane.

These are numerical evaluations supporting separately written proofs, not
interval certificates. No external data or prior review reports are used.
"""

import itertools
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 90
OUT = Path(__file__).resolve().parent


def rho_sq(d, c, tau=mp.mpf("0")):
    a = d * d - c - 1 - tau
    return 2 * c / (a + mp.sqrt(a * a - 4 * c))


def mstr(x):
    return mp.nstr(x, 24)


def saturated(v, pi, j, r):
    delta = [v[h] - v[j] for h in range(len(v)) if h != j]
    others = [h for h in range(len(v)) if h != j]
    ds = [mp.sqrt(mp.fdot(z, z)) for z in delta]
    a = [pi[h] / pi[j] * mp.exp(-d * d / 2 + r * d) for h, d in zip(others, ds)]
    s = mp.fsum(a)
    mean = sum((a_h * z for a_h, z in zip(a, delta)), mp.matrix(len(v[0]), 1)) / (1 + s)
    variance = mp.fsum(a_h * d * d for a_h, d in zip(a, ds)) / (1 + s) - mp.fdot(mean, mean)
    return {
        "S": mstr(s),
        "saturated_trace": mstr(variance),
        "saturated_mean_norm": mstr(mp.sqrt(mp.fdot(mean, mean))),
        "triangle_mean_bound": mstr(mp.fsum(a_h * d for a_h, d in zip(a, ds)) / (1 + s)),
        "odds": [mstr(x) for x in a],
    }


def local_radius_example():
    v = [mp.matrix([0]), mp.matrix([3]), mp.matrix([7])]
    pi = [mp.mpf(50) / 101, mp.mpf(50) / 101, mp.mpf(1) / 101]
    tau = mp.mpf("0.01")
    local_d = [mp.mpf(3), mp.mpf(3), mp.mpf(4)]
    radii = [mp.sqrt(rho_sq(d, mp.mpf(1), tau)) for d in local_d]
    return {
        "means": [0, 3, 7],
        "weights": ["50/101", "50/101", "1/101"],
        "beta": "0.9",
        "tau": "0.01",
        "local_radii": [mstr(x) for x in radii],
        "local": [saturated(v, pi, j, radii[j]) for j in range(3)],
        "common": [saturated(v, pi, j, radii[0]) for j in range(3)],
    }


def local_dimension_example():
    v = [mp.matrix([0, 0]), mp.matrix([mp.mpf("2.3"), 0]), mp.matrix([0, 20])]
    pi = [mp.mpf(1) / 3] * 3
    beta, tau = mp.mpf("0.95"), mp.mpf("0.01")
    c = mp.mpf("1.001")
    radii = [mp.mpf("0.6"), mp.mpf("0.6"), mp.mpf("0.06")]
    local_d = [mp.mpf("2.3"), mp.mpf("2.3"), mp.mpf(20)]
    residual_indices = [1, 1, 0]
    eps = []
    for j in range(3):
        value = mp.mpf(0)
        for h in range(3):
            if h == j:
                continue
            z = v[h] - v[j]
            d = mp.sqrt(mp.fdot(z, z))
            value += mp.exp(-d*d/2 + mp.sqrt(2)*d) * z[residual_indices[j]]**2
        eps.append(value)
    return {
        "means": [[0, 0], ["2.3", 0], [0, 20]],
        "weights": ["1/3"] * 3,
        "beta": mstr(beta),
        "tau": mstr(tau),
        "q": 2,
        "c": mstr(c),
        "local_radii": [mstr(x) for x in radii],
        "epsilon": [mstr(x) for x in eps],
        "trace_budget_slack": [mstr(c - 1 - x / beta) for x in eps],
        "local_separation_slack": [mstr(d*d - (c + 1 + tau + r*r + c/(r*r))) for d, r in zip(local_d, radii)],
        "global_q_floor_margin": mstr(local_d[0]**2 - 2*(2+1) - tau),
        "local": [saturated(v, pi, j, radii[j]) for j in range(3)],
        "operator_invariance_bounds": [mstr(mp.sqrt(beta * rho_sq(d / mp.sqrt(beta), c))) for d in local_d],
    }


def equilateral_example():
    d = mp.mpf("2.8")
    c = mp.mpf(2)
    beta, tau = mp.mpf("0.99"), mp.mpf("0.01")
    r = mp.sqrt(rho_sq(d, c, tau))
    a = mp.exp(-d*d/2 + r*d)
    delta = np.array([[float(d), 0.0], [float(d)/2, float(d)*math.sqrt(3)/2]])
    aa = float(a)
    centre = aa*delta.sum(axis=0)/(1+2*aa)
    vertices = []
    for active in itertools.product([0, 1], repeat=2):
        odds = aa*np.array(active)
        weights = np.r_[1.0, odds]/(1+odds.sum())
        residual = np.vstack([np.zeros(2), delta])-centre
        moment = np.einsum("i,ij,ik->jk", weights, residual, residual)
        vertices.append({"active": active, "lambda_max": float(np.linalg.eigvalsh(moment)[-1])})
    v = [mp.matrix([0, 0]), mp.matrix([d, 0]), mp.matrix([d/2, d*mp.sqrt(3)/2])]
    return {
        "d": mstr(d), "beta": mstr(beta), "tau": mstr(tau), "r": mstr(r),
        "saturated": saturated(v, [mp.mpf(1)/3]*3, 0, r),
        "fixed_centre_vertex_operator_values": vertices,
        "operator_invariance_bound": mstr(mp.sqrt(beta*rho_sq(d/mp.sqrt(beta), c))),
        "mass_bound": mstr(rho_sq(d/mp.sqrt(beta), c)),
    }


def generic_contraction_counterexample():
    # Equal covariance two-component mixture. Even with the separation-based
    # all-mode coverage radius, actual contraction does not imply invariance.
    d = mp.mpf(3)
    ratio = mp.mpf(10000)
    r = mp.mpf("0.4")
    outputs = []
    for y in [-r, mp.mpf(0), r]:
        odds = ratio*mp.exp(-d*d/2+d*y)
        w = odds/(1+odds)
        outputs.append({"y": mstr(y), "T": mstr(d*w), "DT": mstr(d*d*w*(1-w))})
    lower_odds = ratio*mp.exp(-d*d/2-d*r)
    upper_odds = ratio*mp.exp(-d*d/2+d*r)
    return {
        "means": [0, 3], "prior_odds": 10000, "radius": "0.4",
        "q": 1, "tau": "0.01",
        "coverage_separation_slack": mstr(d*d - (2+mp.mpf("0.01")+r*r+1/(r*r))),
        "posterior_odds_range": [mstr(lower_odds), mstr(upper_odds)],
        "samples": outputs,
        "sup_DT_on_ball": mstr(d*d*lower_odds/(1+lower_odds)**2),
        "inf_T_on_ball": mstr(d*lower_odds/(1+lower_odds)),
        "robust_polytope_contains_half_half": bool(upper_odds >= 1),
    }


if __name__ == "__main__":
    result = {
        "local_radius_strict_improvement": local_radius_example(),
        "local_dimension_strict_improvement": local_dimension_example(),
        "equilateral_operator_application": equilateral_example(),
        "actual_contraction_is_insufficient": generic_contraction_counterexample(),
    }
    destination = OUT / "geometry_checks.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
