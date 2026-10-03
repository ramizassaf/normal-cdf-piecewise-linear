"""Write LaTeX tables for the paper directly from results/results.json,
so no number is typed by hand. Also prints supplementary checks.
Run after compute_results.py:  python analysis/make_tables.py
"""
import json
import os

import numpy as np
from scipy.special import ndtr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, "results", "results.json")))
OUT = os.path.join(ROOT, "paper", "tables")
os.makedirs(OUT, exist_ok=True)


def f4(x):
    return f"{x:.4f}".replace("-", "$-$")


def s4(x):  # signed, 4 decimals
    s = f"{x:+.4f}"
    return s.replace("-", "$-$").replace("+", "$+$") if abs(x) > 5e-5 else "0.0000"


def w(name, body):
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write(body)


# Table: breakpoint values and secant slopes
bp = R["breakpoint_cdf"]
sl = R["model_A_slopes"]
rows = []
for i, k in enumerate(["0", "1", "2", "3"]):
    rows.append(f"{k} & {bp[k]:.6f} & " + (f"[{k},{int(k)+1}] & {sl[i]:.6f}" if i < 3 else "-- & --") + r" \\")
w("tab_breakpoints.tex", "\n".join(rows))

# Table: point values A, B, C
rows = []
for r in R["point_table"]:
    rows.append(f"{r['z']:.1f} & {r['phi']:.4f} & {r['A']:.4f} & {s4(r['e_A'])} & "
                f"{r['B']:.4f} & {s4(r['e_B'])} & {r['C']:.4f} & {s4(r['e_C'])}" + r" \\")
w("tab_points.tex", "\n".join(rows))

# Table: summary statistics for all models
names = {"A": "A: exact secant", "B": "B: rounded (published)", "C": "C: 50--84--98--100",
         "D": "D: minimax, integer $b$", "E": "E: minimax, free $b$", "F": "F: two-decimal optimum"}
rows = []
for k in ["A", "B", "C", "D", "E", "F"]:
    m = R["models"][k]
    rows.append(f"{names[k]} & {m['max_abs_error']:.4f} & {m['argmax_z']:.2f} & "
                f"{s4(m['mean_signed_error'])} & {m['mean_abs_error']:.4f} & {m['rmse']:.4f} & "
                f"{m['p95_abs_error']:.4f} & {s4(m['max_positive_error'])} & "
                f"{m['max_rel_error_pct']:.2f}" + r" \\")
w("tab_summary.tex", "\n".join(rows))

# Table: optimised model constants
rows = []
for k in ["D", "E", "F"]:
    d = R[f"model_{k}"]
    kn = ", ".join(f"{v:.4f}".rstrip("0").rstrip(".") if v in (0, 1, 2, 3) else f"{v:.4f}" for v in d["knots"])
    vals = ", ".join(f"{v:.4f}" for v in d["values"])
    rows.append(f"{k} & ({kn}) & ({vals}) & {R['models'][k]['max_abs_error']:.4f}" + r" \\")
w("tab_optimised.tex", "\n".join(rows))

# Table: frontier
rows = []
for fr in R["frontier"]:
    a = fr["equal_width_interp"]; b = fr["equal_width_minimax"]; c = fr["free_minimax"]
    rows.append(f"{fr['segments']} & {a['constants']} & {a['max_abs_error']:.4f} & "
                f"{b['max_abs_error']:.4f} & {c['constants']} & {c['max_abs_error']:.4f}" + r" \\")
w("tab_frontier.tex", "\n".join(rows))

# Table: reference approximations
ref = R["reference_approximations"]
lab = {"Polya_1949": "P\\'olya (1949) [3]", "Bowling_2009": "Bowling et al. (2009) [4]",
       "AS_26_2_17": "Abramowitz--Stegun 26.2.17 [1]"}
rows = []
for k in ["Polya_1949", "Bowling_2009", "AS_26_2_17"]:
    m = ref[k]
    rows.append(f"{lab[k]} & {m['max_abs_error']:.2e} & {m['argmax_z']:.2f}".replace("e-0", r"\times 10^{-").replace("e-", r"\times 10^{-") + "" )
rows = []
for k, ops in [("Polya_1949", "1 exp, 1 sqrt"), ("Bowling_2009", "1 exp, 1 division"),
               ("AS_26_2_17", "1 exp, 6 constants")]:
    m = ref[k]
    mant, ex = f"{m['max_abs_error']:.2e}".split("e")
    rows.append(f"{lab[k]} & ${mant}\\times10^{{{int(ex)}}}$ & {m['argmax_z']:.2f} & {ops}" + r" \\")
w("tab_reference.tex", "\n".join(rows))

# Table: tail / defect rates
rows = []
for t in R["tail_table"]:
    rows.append(f"{t['z']:.1f} & {t['two_sided_ppm_exact']:,.0f} & {t['ppm_A']:,.0f} & "
                f"{t['relQ_A_pct']:+.1f} & {t['ppm_C']:,.0f} & {t['relQ_C_pct']:+.1f} & "
                f"{t['ppm_E']:,.0f} & {t['relQ_E_pct']:+.1f}".replace("-", "$-$").replace("+", "$+$") + r" \\")
w("tab_tail.tex", "\n".join(rows))

# Table: corrected defects table
rows = []
for d in R["defects_table"]:
    rows.append(f"{d['z']:.1f} & {d['phi']:.5f} & {100*d['two_sided_coverage']:.3f} & "
                f"{d['two_sided_ppm']:,.0f} & {d['one_sided_ppm']:,.0f}" + r" \\")
w("tab_defects.tex", "\n".join(rows))

# ---------------- supplementary checks ----------------
z = np.linspace(0, 3, 30001); phi = ndtr(z)
E = R["model_E"]
kE4 = np.round(E["knots"], 4); yE4 = np.round(E["values"], 4)
print("E rounded to 4 dp max|e|:", np.max(np.abs(np.interp(z, kE4, yE4) - phi)))
for k in ["D", "E", "F"]:
    print(k, "monotone:", bool(np.all(np.diff(R[f'model_{k}']['values']) >= 0)))
# classical bound for linear interpolation: h^2/8 * max|f''|, f'' = -z phi(z), max at z = 1
from scipy.stats import norm
print("h^2/8 max|Phi''| (h=1):", norm.pdf(1) / 8)
# spreadsheet formula check (symmetric Model C)
def sheet(x):
    a = abs(x)
    if a > 3: return float("nan")
    h = 0.34 * a if a <= 1 else (0.34 + 0.14 * (a - 1) if a <= 2 else 0.48 + 0.02 * (a - 2))
    return 0.5 + (1 if x >= 0 else -1) * h
C = lambda x: np.interp(abs(x), [0, 1, 2, 3], [0.5, 0.84, 0.98, 1.0])
xs = np.linspace(-3, 3, 601)
print("sheet vs model C (max diff incl. symmetry):",
      max(abs(sheet(x) - (C(x) if x >= 0 else 1 - C(x))) for x in xs))
print("tables written to", OUT)
