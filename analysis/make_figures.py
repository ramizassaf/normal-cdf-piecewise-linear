"""Generate all paper figures from results/results.json and the model definitions.
Run after compute_results.py:  python analysis/make_figures.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "paper", "figures")
R = json.load(open(os.path.join(ROOT, "results", "results.json")))

plt.rcParams.update({
    "font.family": ["TeX Gyre Pagella", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.3,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

INK = "#1f2a36"
COL = {"A": "#1b5e9a", "B": "#c0392b", "C": "#d4892a", "D": "#6c4a9e",
       "E": "#1f8a70", "F": "#7f8c8d", "P": "#444444"}

z = np.linspace(0, 3, 6001)
phi = ndtr(z)


def pwl(x, k, y):
    return np.interp(x, k, y)


A = pwl(z, [0, 1, 2, 3], ndtr(np.array([0, 1, 2, 3.0])))
B = np.where(z <= 1, 0.5 + 0.34 * z, np.where(z <= 2, 0.841 + 0.14 * (z - 1), 0.977 + 0.02 * (z - 2)))
C = pwl(z, [0, 1, 2, 3], [0.50, 0.84, 0.98, 1.00])
D = pwl(z, R["model_D"]["knots"], R["model_D"]["values"])
E = pwl(z, R["model_E"]["knots"], R["model_E"]["values"])
F = pwl(z, R["model_F"]["knots"], R["model_F"]["values"])
P = 0.5 + 0.5 * np.sqrt(1 - np.exp(-2 * z ** 2 / np.pi))


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    fig.savefig(os.path.join(FIG, name + ".png"))
    plt.close(fig)


# Figure 1: CDF and approximations
fig, ax = plt.subplots(figsize=(5.6, 3.0))
ax.plot(z, phi, color=INK, lw=2.2, label=r"Exact $\Phi(z)$")
ax.plot(z, A, color=COL["A"], ls="--", label="Model A (exact secant)")
ax.plot(z, E, color=COL["E"], ls="-.", label="Model E (minimax, free breakpoints)")
ax.scatter([0, 1, 2, 3], ndtr(np.array([0, 1, 2, 3.0])), s=14, color=COL["A"], zorder=5)
ax.set_xlabel(r"$z$")
ax.set_ylabel("Cumulative probability")
ax.set_xlim(0, 3)
ax.set_ylim(0.48, 1.01)
ax.legend(frameon=False, loc="lower right")
ax.grid(alpha=0.25, lw=0.4)
save(fig, "fig1_cdf")

# Figure 2: signed error (two panels)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.3), sharey=True)
ax = axs[0]
ax.axhline(0, color="#999", lw=0.6)
ax.plot(z, 100 * (A - phi), color=COL["A"], label="A (exact secant)")
for seg in [(z <= 1), (z > 1) & (z <= 2), (z > 2)]:
    ax.plot(z[seg], 100 * (B - phi)[seg], color=COL["B"], ls="--",
            label="B (rounded, discontinuous)" if seg[0] else None)
ax.plot(z, 100 * (C - phi), color=COL["C"], ls=":", lw=1.6, label="C (50-84-98-100)")
ax.set_title("(a) Interpolation-based models")
ax.set_xlabel(r"$z$")
ax.set_ylabel(r"Signed error $e(z)$ (percentage points)")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))
ax.grid(alpha=0.25, lw=0.4)
ax = axs[1]
ax.axhline(0, color="#999", lw=0.6)
ax.plot(z, 100 * (D - phi), color=COL["D"], label="D (minimax, integer breakpoints)")
ax.plot(z, 100 * (E - phi), color=COL["E"], ls="-.", label="E (minimax, free breakpoints)")
ax.plot(z, 100 * (F - phi), color=COL["F"], ls="--", label="F (two-decimal constrained)")
ax.set_title("(b) Optimised models")
ax.set_xlabel(r"$z$")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))
ax.grid(alpha=0.25, lw=0.4)
for a in axs:
    a.set_xlim(0, 3)
save(fig, "fig2_error")

# Figure 3: accuracy-complexity frontier
fr = R["frontier"]
fig, ax = plt.subplots(figsize=(5.0, 3.0))
for key, lab, col, mk in [("equal_width_interp", "Equal-width interpolation", COL["A"], "o"),
                          ("equal_width_minimax", "Equal-width, minimax values", COL["D"], "s"),
                          ("free_minimax", "Free breakpoints, minimax values", COL["E"], "^")]:
    xs = [f[key]["constants"] for f in fr]
    ys = [100 * f[key]["max_abs_error"] for f in fr]
    ax.plot(xs, ys, marker=mk, ms=4, color=col, label=lab)
ax.set_yscale("log")
ax.set_xlabel("Number of stored constants")
ax.set_ylabel(r"Max $|e(z)|$ on $[0,3]$ (pp)")
ax.legend(frameon=False)
ax.grid(alpha=0.25, lw=0.4, which="both")
save(fig, "fig3_frontier")

# Figure 4: relative error in the upper tail probability Q(z) = 1 - Phi(z)
Q = 1 - phi
fig, ax = plt.subplots(figsize=(5.6, 3.0))
m = z >= 0.5
for arr, lab, col, ls in [(A, "Model A", COL["A"], "-"), (C, "Model C", COL["C"], ":"),
                          (E, "Model E", COL["E"], "-."), (P, "P\u00f3lya (1949)", COL["P"], "--")]:
    rel = 100 * ((1 - arr) - Q) / Q
    ax.plot(z[m], rel[m], color=col, ls=ls, label=lab)
ax.axhline(0, color="#999", lw=0.6)
ax.set_ylim(-110, 160)
ax.set_xlim(0.5, 3)
ax.set_xlabel(r"$z$")
ax.set_ylabel(r"Relative error in $Q(z)=1-\Phi(z)$ (%)")
ax.legend(frameon=False, loc="upper left", ncol=2)
ax.grid(alpha=0.25, lw=0.4)
save(fig, "fig4_tail")
print("figures written to", FIG)
