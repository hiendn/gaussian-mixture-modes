"""Reproduce the density illustration in the manuscript.

This script evaluates a specified density and locates the three plotted modes
by posterior-mean iteration. Exhaustiveness is supplied by the theorem and
certificate, not by this numerical search.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
d, radius = 2.8, 0.68
means = np.array([[0.0, 0.0], [d, 0.0], [d / 2, np.sqrt(3) * d / 2]])

def posterior_mean(y):
    logw = -0.5 * ((means - y) ** 2).sum(axis=1)
    w = np.exp(logw - logw.max())
    w /= w.sum()
    return w @ means

modes, residuals = [], []
for mean in means:
    y = mean.copy()
    for iteration in range(10000):
        yn = posterior_mean(y)
        if np.linalg.norm(yn - y) < 1e-14:
            y = yn
            break
        y = yn
    modes.append(y)
    residuals.append(float(np.linalg.norm(posterior_mean(y) - y)))
modes = np.asarray(modes)
assert max(residuals) < 1e-12
assert np.all(np.linalg.norm(modes - means, axis=1) < radius)

xx = np.linspace(-1.45, 4.25, 520)
yy = np.linspace(-1.4, 3.9, 490)
x, y = np.meshgrid(xx, yy)
grid = np.stack([x, y], axis=-1)
density = np.exp(-0.5 * ((grid[..., None, :] - means) ** 2).sum(axis=-1)).mean(axis=-1) / (2 * np.pi)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})
fig, ax = plt.subplots(figsize=(6.4, 5.0), layout="constrained")
levels = np.linspace(0.004, 0.058, 10)
cs = ax.contour(x, y, density, levels=levels, cmap="Blues", linewidths=1.15)
ax.clabel(cs, cs.levels[1::2], inline=True, fontsize=8, fmt="%.3f")
for mean in means:
    ax.add_patch(Circle(mean, radius, fill=False, ec="#67717a", ls="--", lw=1.3))
ax.scatter(means[:, 0], means[:, 1], marker="x", s=65, lw=1.8,
           color="#bc570e", label="Component means", zorder=4)
ax.scatter(modes[:, 0], modes[:, 1], s=30, color="#173956",
           label="Modes", zorder=5)
ax.set(xlabel="Whitened coordinate 1", ylabel="Whitened coordinate 2",
       xlim=(xx.min(), xx.max()), ylim=(yy.min(), yy.max()), aspect="equal")
ax.legend(loc="upper right", frameon=False, fontsize=9)
fig.savefig(OUT / "triangle_modes.pdf")
fig.savefig(OUT / "triangle_modes.png", dpi=180)
plt.close(fig)
(OUT / "triangle_mode_coordinates.json").write_text(json.dumps({
    "means": means.tolist(), "modes": modes.tolist(),
    "posterior_mean_residuals": residuals, "certified_radius": radius,
    "covariance": [[1, 0], [0, 1]], "weights": [1/3]*3,
    "purpose": "Figure coordinates only; not an exhaustive numerical mode count."
}, indent=2) + "\n")
print(json.dumps({"max_fixed_point_residual": max(residuals),
                  "figure": str(OUT / "triangle_modes.pdf")}))
