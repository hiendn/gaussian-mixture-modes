"""Exact rational sign certificates for two local-geometry examples.

No floating-point arithmetic is used to decide a sign. Exponential brackets
use a rational Taylor sum and a geometric bound on its remaining positive
terms. Floating-point midpoints appear only as clearly marked display values.
"""

from fractions import Fraction as F
import json
from pathlib import Path


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def wrap(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, other):
        other = self.wrap(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.wrap(other)

    def __rsub__(self, other):
        return self.wrap(other) + -self

    def __mul__(self, other):
        other = self.wrap(other)
        products = [self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self):
        assert not self.lo <= 0 <= self.hi
        return Interval(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * self.wrap(other).reciprocal()

    def __rtruediv__(self, other):
        return self.wrap(other) * self.reciprocal()

    def square(self):
        if self.lo <= 0 <= self.hi:
            return Interval(0, max(self.lo**2, self.hi**2))
        return Interval(min(self.lo**2, self.hi**2), max(self.lo**2, self.hi**2))

    def display_midpoint(self):
        return float((self.lo + self.hi) / 2)


def exp_negative(x, n=100):
    """A rational enclosure for exp(-x), x >= 0.

    S_n < exp(x) < S_n + t_(n+1)/(1-x/(n+2)) for x>0,
    because successive tail ratios are at most x/(n+2).
    """
    x = F(x)
    assert 0 <= x < n + 2
    term = F(1)
    total = F(1)
    for k in range(1, n + 1):
        term *= x / k
        total += term
    next_term = term * x / (n + 1)
    upper = total + next_term / (1 - x / (n + 2))
    return Interval(1 / upper, 1 / total)


def separation_margin(d_squared, radius, nu, tau):
    return d_squared - (nu + 1 + tau + radius**2 + nu / radius**2)


def strict_bound(name, quantity, direction, threshold):
    quantity = Interval.wrap(quantity)
    threshold = F(threshold)
    if direction == "less_than":
        assert quantity.hi < threshold, name
    elif direction == "greater_than":
        assert quantity.lo > threshold, name
    else:
        raise ValueError(direction)
    return {
        "quantity": name,
        "certified_relation": direction,
        "rational_threshold": str(threshold),
        "threshold_decimal_for_display": float(threshold),
        "midpoint_for_display_only": quantity.display_midpoint(),
        "sign_decision_arithmetic": "exact rational",
    }


def one_dimensional_example():
    means = [F(0), F(3), F(7)]
    weights = [F(50, 101), F(50, 101), F(1, 101)]
    radii = [F("0.383"), F("0.383"), F("0.269")]
    local_distances = [F(3), F(3), F(4)]
    results = []
    for j in range(3):
        odds, deltas = [], []
        for h in range(3):
            if h == j:
                continue
            delta = means[h] - means[j]
            distance = abs(delta)
            odds.append(weights[h] / weights[j] *
                        exp_negative(distance**2 / 2 - radii[j]*distance))
            deltas.append(delta)
        s = sum(odds, Interval(0))
        denominator = 1 + s
        mean = sum((a*d for a, d in zip(odds, deltas)), Interval(0)) / denominator
        second = sum((a*d*d for a, d in zip(odds, deltas)), Interval(0)) / denominator
        variance = second - mean.square()
        margin = separation_margin(local_distances[j]**2, radii[j], F(1), F("0.01"))
        results.extend([
            strict_bound(f"j{j+1}_separation_slack", margin, "greater_than", ["0.026", "0.026", "0.098"][j]),
            strict_bound(f"j{j+1}_S", s, "less_than", "0.05"),
            strict_bound(f"j{j+1}_saturated_trace", variance, "less_than", ["0.30", "0.30", "0.72"][j]),
        ])

    # H_{1,.01}(r)^2 strictly decreases on (0,1). Its value at .382 is
    # greater than 9, so every admissible common radius is greater than .382.
    common_lower = F("0.382")
    common_exclusion = -separation_margin(F(9), common_lower, F(1), F("0.01"))
    results.append(strict_bound("common_radius_.382_separation_failure", common_exclusion,
                                "greater_than", "0.0088"))

    # At any admissible common radius, P_3 contains the fixed two-point
    # distribution with odds a*=50 exp(-8+4*.382) on means 7 and 3.
    a_star = F(50) * exp_negative(F(8) - 4*common_lower)
    two_point_variance = 16*a_star/(1+a_star).square()
    results.append(strict_bound("fixed_two_point_robust_variance_in_common_ball",
                                two_point_variance, "greater_than", "1.06"))
    return {
        "means": [0, 3, 7],
        "weights": ["50/101", "50/101", "1/101"],
        "beta": "0.9", "tau": "0.01",
        "radii": ["0.383", "0.383", "0.269"],
        "exponential_enclosure": "Rational Taylor sum through degree 100 plus geometric positive-tail bound",
        "checks": results,
    }


def remote_component_example():
    beta, tau, nu = F("0.95"), F("0.01"), F("1.001")
    radii = [F("0.6"), F("0.6"), F("0.06")]
    distances_squared = [F("5.29"), F("5.29"), F(400)]
    results = []

    # Every remote distance is between 20 and 21. With sqrt(2)<3/2,
    # the coarse exponent is less than -170. The residual squared norm
    # is at most 400. Since e>2, epsilon_j < 400*2^-170 for every j.
    eps_bound = F(400, 2**170)
    results.append(strict_bound("coarse_residual_epsilon_upper", eps_bound,
                                "less_than", F(1, 10**45)))
    results.append(strict_bound("effective_budget_slack_lower", nu-1-eps_bound/beta,
                                "greater_than", "0.0009"))

    for j in range(3):
        margin = separation_margin(distances_squared[j], radii[j], nu, tau)
        results.append(strict_bound(f"j{j+1}_separation_slack", margin,
                                    "greater_than", ["0.138", "0.138", "119.9"][j]))

    # For the close pair, a=exp(-1.265)<.283. This inequality alone
    # can already be certified with a Taylor lower sum through degree eight.
    x = F("1.265")
    term, taylor_lower = F(1), F(1)
    for k in range(1, 9):
        term *= x/k
        taylor_lower += term
    assert taylor_lower > F(1000, 283)
    pair_odds_bound = F("0.283")
    pair_variance_bound = F("5.29")*pair_odds_bound/(1+pair_odds_bound)**2
    results.append(strict_bound("close_pair_variance_upper", pair_variance_bound,
                                "less_than", "0.91"))

    # A remote odds cap at radius .6 is below exp(-188)<2^-188.
    # Adding a mass t at a point less than 21 away changes the base
    # two-point variance by at most 24^2*t, and t is below that odds.
    remote_odds_bound = F(1, 2**188)
    full_variance_bound = pair_variance_bound + 576*remote_odds_bound
    results.append(strict_bound("j1_j2_saturated_trace_upper", full_variance_bound,
                                "less_than", "0.911"))
    results.append(strict_bound("j1_j2_S_upper", pair_odds_bound+remote_odds_bound,
                                "less_than", "0.284"))

    # At radius .06, each remote exponent is below -198. Both squared
    # distances are less than 441. Discarding mean-centering is enough.
    isolated_odds_bound = F(1, 2**198)
    results.append(strict_bound("j3_saturated_trace_upper", 2*441*isolated_odds_bound,
                                "less_than", F(1, 10**50)))
    results.append(strict_bound("j3_S_upper", 2*isolated_odds_bound,
                                "less_than", F(1, 10**50)))
    results.append(strict_bound("old_global_q2_floor_failure", F("6.01")-F("5.29"),
                                "greater_than", "0.7"))
    return {
        "means": [[0, 0], ["2.3", 0], [0, 20]],
        "weights": ["1/3"]*3,
        "beta": "0.95", "tau": "0.01", "nu": "1.001",
        "radii": ["0.6", "0.6", "0.06"],
        "tail_proof": "sqrt(2)<3/2, monotonicity of -d^2/2+s*d for d>=20, and e>2",
        "checks": results,
    }


if __name__ == "__main__":
    result = {
        "arithmetic": "Python Fraction throughout all sign decisions",
        "one_dimensional_local_radius_example": one_dimensional_example(),
        "two_dimensional_local_budget_example": remote_component_example(),
    }
    output = Path(__file__).resolve().parent / "rigorous_geometry_checks.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
