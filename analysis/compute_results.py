"""
Reproducible analysis for:
"Low-Complexity Piecewise-Linear Approximations of the Standard Normal CDF
 for Engineering Calculations: An Accuracy-Complexity Analysis"

Every number reported in the paper is produced by this script.
Run:  python analysis/compute_results.py
Outputs: results/results.json, results/error_grid.csv, paper/figures/*.pdf|png
"""
import json
import os
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import linprog, minimize
from scipy.special import ndtr
import mpmath

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "paper", "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

Z_MAX = 3.0
STEP = 1e-4
Z = np.round(np.arange(0.0, Z_MAX + STEP / 2, STEP), 10)  # 30,001 points
PHI = ndtr(Z)


# ---------------------------------------------------------------------------
# 0. Reference check: scipy ndtr against 50-digit mpmath
# ---------------------------------------------------------------------------
mpmath.mp.dps = 50
zc = Z[::10]  # 3,001 points
ref = np.array([float(mpmath.ncdf(mpmath.mpf(str(v)))) for v in zc])
ref_check = float(np.max(np.abs(ref - ndtr(zc))))


# ---------------------------------------------------------------------------
# 1. Model definitions
# ---------------------------------------------------------------------------
def pwl(z, knots, values):
    """Continuous piecewise-linear interpolation through (knots, values)."""
    return np.interp(z, knots, values)


def model_A(z):
    """Exact-secant interpolation through (k, Phi(k)), k = 0..3. Unrounded."""
    k = np.array([0.0, 1.0, 2.0, 3.0])
    return pwl(z, k, ndtr(k))


def model_B(z):
    """Rounded-coefficient model as originally published (discontinuous).
    Segments are left-open/right-closed except the first: [0,1], (1,2], (2,3]."""
    z = np.asarray(z, dtype=float)
    return np.where(
        z <= 1, 0.5 + 0.34 * z,
        np.where(z <= 2, 0.841 + 0.14 * (z - 1), 0.977 + 0.02 * (z - 2)))


def model_C(z):
    """Two-digit node model '50-84-98-100': nodes 0.50, 0.84, 0.98, 1.00 at 0,1,2,3.
    Continuous; slopes 0.34, 0.14, 0.02."""
    return pwl(z, [0, 1, 2, 3], [0.50, 0.84, 0.98, 1.00])


def stats(approx, z=Z, phi=PHI):
    e = approx - phi
    ae = np.abs(e)
    i = int(np.argmax(ae))
    return {
        "max_abs_error": float(ae[i]),
        "max_abs_error_pp": float(100 * ae[i]),
        "argmax_z": float(z[i]),
        "sign_at_max": float(np.sign(e[i])),
        "mean_signed_error": float(np.mean(e)),
        "mean_abs_error": float(np.mean(ae)),
        "rmse": float(np.sqrt(np.mean(e ** 2))),
        "p95_abs_error": float(np.percentile(ae, 95)),
        "max_positive_error": float(np.max(e)),
        "max_negative_error": float(np.min(e)),
        "max_rel_error_pct": float(np.max(ae / phi) * 100),
        "argmax_rel_z": float(z[int(np.argmax(ae / phi))]),
        "fraction_underestimating": float(np.mean(e < -1e-12)),
    }


# ---------------------------------------------------------------------------
# 2. Minimax optimisation (linear programme for fixed breakpoints)
# ---------------------------------------------------------------------------
def minimax_values(knots, z=Z, phi=PHI, fix_zero=True):
    """Given breakpoints, find node values minimising max|Lhat - Phi| on the grid.
    Variables: node values y_0..y_k and t. Constraint y_0 = 0.5 keeps
    symmetry Phi(-z) = 1 - Phi(z) continuous at z = 0."""
    knots = np.asarray(knots, float)
    n = len(knots)
    # Basis matrix: Lhat(z) = sum_j y_j * hat_j(z)
    Bm = np.zeros((len(z), n))
    for j in range(n):
        v = np.zeros(n); v[j] = 1.0
        Bm[:, j] = np.interp(z, knots, v)
    # minimise t   s.t.  B y - t <= phi,  -B y - t <= -phi
    c = np.zeros(n + 1); c[-1] = 1.0
    A_ub = np.vstack([np.hstack([Bm, -np.ones((len(z), 1))]),
                      np.hstack([-Bm, -np.ones((len(z), 1))])])
    b_ub = np.concatenate([phi, -phi])
    bounds = [(0.0, 1.0)] * n + [(0, None)]   # a CDF approximation must stay in [0, 1]
    if fix_zero:
        bounds[0] = (0.5, 0.5)
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    return res.x[:n], res.x[-1]


# Coarser grid for the outer breakpoint search (speed); final values re-checked on full grid
ZC = Z[::5]
PHIC = PHI[::5]


def optimal_breakpoints(k, z_max=Z_MAX):
    """Free-breakpoint minimax PWL with k segments on [0, z_max]."""
    if k == 1:
        knots = np.array([0.0, z_max])
        y, _ = minimax_values(knots)
        return knots, y
    best = None
    # multistart from equal-width and a few perturbations
    starts = [np.linspace(0, z_max, k + 1)[1:-1]]
    rng = np.random.default_rng(0)
    for _ in range(6):
        starts.append(np.sort(rng.uniform(0.2, z_max - 0.2, k - 1)))

    def obj(inner):
        inner = np.sort(np.clip(inner, 1e-3, z_max - 1e-3))
        if np.any(np.diff(np.concatenate([[0], inner, [z_max]])) < 0.02):
            return 1.0
        knots = np.concatenate([[0], inner, [z_max]])
        _, t = minimax_values(knots, ZC, PHIC)
        return t

    for s in starts:
        r = minimize(obj, s, method="Nelder-Mead",
                     options={"xatol": 1e-4, "fatol": 1e-7, "maxiter": 600})
        if best is None or r.fun < best.fun:
            best = r
    inner = np.sort(best.x)
    knots = np.concatenate([[0], inner, [z_max]])
    y, _ = minimax_values(knots)
    return knots, y


results = {
    "environment": {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
        "platform": platform.platform(),
        "reference_cdf": "scipy.special.ndtr (Cephes), cross-checked with mpmath.ncdf at 50 digits",
        "grid": f"z = 0 to {Z_MAX} in steps of {STEP} ({len(Z)} points)",
        "max_abs_diff_scipy_vs_mpmath": ref_check,
    }
}

# ---------------------------------------------------------------------------
# 3. Core models A, B, C
# ---------------------------------------------------------------------------
results["breakpoint_cdf"] = {str(k): float(ndtr(k)) for k in [0, 1, 2, 3]}
results["model_A_slopes"] = [float(ndtr(1) - ndtr(0)), float(ndtr(2) - ndtr(1)),
                             float(ndtr(3) - ndtr(2))]

models = {"A": model_A, "B": model_B, "C": model_C}
results["models"] = {}
for name, f in models.items():
    results["models"][name] = stats(f(Z))

# Boundary behaviour for model B (left value, right limit, jump)
eps = 1e-12
results["model_B_boundaries"] = {}
for b in [1.0, 2.0]:
    left = float(model_B(b))
    right = float(model_B(b + eps))
    results["model_B_boundaries"][str(b)] = {
        "value_at_breakpoint": left, "right_limit": round(right, 6),
        "jump": round(right - left, 6), "exact_phi": float(ndtr(b)),
        "error_left": left - float(ndtr(b)), "error_right": round(right, 6) - float(ndtr(b))}
results["model_B_value_at_3"] = float(model_B(3.0))

# Point table
pts = [0, 0.5, 1, 1.5, 2, 2.5, 3]
results["point_table"] = []
for p in pts:
    row = {"z": p, "phi": float(ndtr(p))}
    for name, f in models.items():
        v = float(f(np.array([p]))[0])
        row[name] = v
        row[f"e_{name}"] = v - float(ndtr(p))
    results["point_table"].append(row)

# ---------------------------------------------------------------------------
# 4. Optimised three-segment models
# ---------------------------------------------------------------------------
# D: integer breakpoints, minimax node values
kD = np.array([0.0, 1.0, 2.0, 3.0])
yD, _ = minimax_values(kD)
fD = lambda z: pwl(z, kD, yD)
results["models"]["D"] = stats(fD(Z))
results["model_D"] = {"knots": kD.tolist(), "values": yD.tolist()}

# E: free breakpoints, minimax node values
kE, yE = optimal_breakpoints(3)
fE = lambda z: pwl(z, kE, yE)
results["models"]["E"] = stats(fE(Z))
results["model_E"] = {"knots": kE.tolist(), "values": yE.tolist(),
                      "slopes": (np.diff(yE) / np.diff(kE)).tolist()}

# F: simplicity-constrained optimum.
#    Breakpoints restricted to multiples of 0.5, node values restricted to
#    two decimals, continuity kept, y0 = 0.50. Exhaustive search.
best = None
cands_b = np.arange(0.5, 3.0, 0.5)
for i, b1 in enumerate(cands_b):
    for b2 in cands_b[i + 1:]:
        knots = np.array([0, b1, b2, 3.0])
        yopt, _ = minimax_values(knots, ZC, PHIC)
        grids = [np.round(np.arange(np.round(v, 2) - 0.02, np.round(v, 2) + 0.0201, 0.01), 2)
                 for v in yopt[1:]]
        for y1 in grids[0]:
            for y2 in grids[1]:
                for y3 in grids[2]:
                    y = np.array([0.5, y1, y2, min(y3, 1.0)])
                    if np.any(np.diff(y) < 0):
                        continue
                    m = np.max(np.abs(pwl(ZC, knots, y) - PHIC))
                    if best is None or m < best[0]:
                        best = (m, knots.copy(), y.copy())
kF, yF = best[1], best[2]
fF = lambda z: pwl(z, kF, yF)
results["models"]["F"] = stats(fF(Z))
results["model_F"] = {"knots": kF.tolist(), "values": yF.tolist(),
                      "slopes": (np.diff(yF) / np.diff(kF)).tolist()}

# ---------------------------------------------------------------------------
# 5. Accuracy-complexity frontier
#    complexity = number of stored constants
#    equal-width interpolation: k+1 node values (breakpoints implied)
#    free-breakpoint minimax: (k+1 node values) + (k-1 interior breakpoints) = 2k
#    (y0 = 0.5 is counted although it is fixed)
# ---------------------------------------------------------------------------
frontier = []
for k in range(1, 7):
    kn = np.linspace(0, Z_MAX, k + 1)
    m_eq = float(np.max(np.abs(pwl(Z, kn, ndtr(kn)) - PHI)))
    y_eqmm, _ = minimax_values(kn)
    m_eqmm = float(np.max(np.abs(pwl(Z, kn, y_eqmm) - PHI)))
    ko, yo = optimal_breakpoints(k)
    m_opt = float(np.max(np.abs(pwl(Z, ko, yo) - PHI)))
    frontier.append({"segments": k,
                     "equal_width_interp": {"constants": k + 1, "max_abs_error": m_eq},
                     "equal_width_minimax": {"constants": k + 1, "max_abs_error": m_eqmm},
                     "free_minimax": {"constants": 2 * k, "max_abs_error": m_opt,
                                      "knots": ko.tolist()}})
results["frontier"] = frontier

# ---------------------------------------------------------------------------
# 6. Reference non-linear approximations (for context)
# ---------------------------------------------------------------------------
def polya(z):  # Polya (1949)
    return 0.5 + 0.5 * np.sqrt(1 - np.exp(-2 * z ** 2 / np.pi))


def as_26_2_17(z):  # Abramowitz & Stegun 26.2.17, |error| < 7.5e-8
    p = 0.2316419
    b = [0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429]
    t = 1 / (1 + p * z)
    poly = t * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * b[4]))))
    return 1 - np.exp(-z * z / 2) / np.sqrt(2 * np.pi) * poly


def bowling(z):  # Bowling et al. (2009) logistic approximation
    return 1 / (1 + np.exp(-(0.07056 * z ** 3 + 1.5976 * z)))


results["reference_approximations"] = {
    "Polya_1949": stats(polya(Z)),
    "AS_26_2_17": stats(as_26_2_17(Z)),
    "Bowling_2009": stats(bowling(Z)),
}

# ---------------------------------------------------------------------------
# 7. Tail probability Q(z) = 1 - Phi(z) and two-sided defect rates
# ---------------------------------------------------------------------------
tail = []
for p in [1.0, 1.5, 2.0, 2.5, 3.0]:
    Q = 1 - float(ndtr(p))
    row = {"z": p, "Q_exact": Q, "two_sided_ppm_exact": 2 * Q * 1e6,
           "two_sided_yield_exact": 1 - 2 * Q}
    for name, f in [("A", model_A), ("B", model_B), ("C", model_C), ("E", fE), ("F", fF)]:
        Qh = 1 - float(f(np.array([p]))[0])
        row[f"Q_{name}"] = Qh
        row[f"ppm_{name}"] = 2 * Qh * 1e6
        row[f"relQ_{name}_pct"] = 100 * (Qh - Q) / Q
    tail.append(row)
results["tail_table"] = tail

# Defects table (two-sided, centred, no 1.5-sigma shift), z up to 3.5 and C_pk 1.33
defects = []
for p in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
    Q = 1 - float(ndtr(p))
    defects.append({"z": p, "phi": float(ndtr(p)), "two_sided_coverage": 1 - 2 * Q,
                    "two_sided_ppm": 2 * Q * 1e6, "one_sided_ppm": Q * 1e6})
results["defects_table"] = defects

# ---------------------------------------------------------------------------
# 8. Worked examples
# ---------------------------------------------------------------------------
ex = {}
z = 2.5
for name, f in [("A", model_A), ("B", model_B), ("C", model_C), ("F", fF)]:
    ph = float(f(np.array([z]))[0])
    ex[f"bearing_{name}"] = {"phi_hat": ph, "yield": 2 * ph - 1}
ex["bearing_exact"] = {"phi": float(ndtr(z)), "yield": 2 * float(ndtr(z)) - 1}
for zz in [0.5, 1.2, 1.5, 1.8, 2.5]:
    ex[f"z_{zz}"] = {"phi": float(ndtr(zz)),
                     **{n: float(f(np.array([zz]))[0]) for n, f in
                        [("A", model_A), ("B", model_B), ("C", model_C), ("F", fF)]}}
ex["cpk_1.33_yield_exact"] = 2 * float(ndtr(4.0)) - 1
results["examples"] = ex

with open(os.path.join(RES, "results.json"), "w") as fh:
    json.dump(results, fh, indent=2)

# Error grid CSV (every 0.001)
with open(os.path.join(RES, "error_grid.csv"), "w") as fh:
    fh.write("z,phi,A,B,C,D,E,F,e_A,e_B,e_C,e_D,e_E,e_F\n")
    zz = Z[::10]
    ph = PHI[::10]
    cols = [model_A(zz), model_B(zz), model_C(zz), fD(zz), fE(zz), fF(zz)]
    for i in range(len(zz)):
        vals = [c[i] for c in cols]
        fh.write(",".join([f"{zz[i]:.3f}", f"{ph[i]:.10f}"] +
                          [f"{v:.10f}" for v in vals] +
                          [f"{v - ph[i]:.10f}" for v in vals]) + "\n")

print(json.dumps({k: results[k] for k in ["environment", "models", "model_B_boundaries",
                                          "model_D", "model_E", "model_F"]}, indent=1))
print(json.dumps(results["frontier"], indent=1))
print(json.dumps(results["reference_approximations"], indent=1))
