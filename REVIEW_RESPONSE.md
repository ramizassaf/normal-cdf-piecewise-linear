# Response to Reviewer

Paper: *Low-Complexity Piecewise-Linear Approximations of the Standard Normal CDF for Engineering Calculations*, version 2.0.

Thank you for the careful review. Every point below was addressed in code first, then in the text. All numbers come from `analysis/compute_results.py`; tables are generated from `results/results.json`.

## Final technical pass (8 points)

**1. Probability units vs percentage points.** Section 3 now defines both. Errors are reported in probability units throughout (0.0240) and converted to percentage points only where marked (2.40 pp). The v1 labels such as "0.0208 pp" were wrong by a factor of 100 and are gone. Table captions state the unit.

**2. Signed error.** `E_abs` is replaced by the signed error e(z) = Φ̂(z) − Φ(z). Absolute error |e(z)|, relative CDF error, and relative tail error r_Q are defined separately in Section 3. Tables report signed and absolute statistics in separate columns, including the largest positive error.

**3. 50-84-98 examples.** Recomputed: z = 0.5 gives 0.670 vs 0.6915 (−2.15 pp, v1 said 1.6); z = 1.5 gives 0.910 vs 0.9332 (−2.32 pp); z = 2.5 gives 0.990 vs 0.9938 (−0.38 pp, v1 said 0.6). The "z/3", "(z−1)/7" and "10% rule" shortcuts were removed: the 10% rule gave Φ(1.5) ≈ 0.95, which is wrong. The v1 two-sided example at z = 1.2 used the first-segment formula on the second segment; it is corrected in Section 7.2.

**4. Defects table.** Rebuilt from exact values for a centred process with no 1.5σ shift. v1 errors: z = 2.0 is 45,500 ppm (v1: 4,550); z = 2.5 is 12,419 ppm (v1: 124); z = 3.0 is 2,700 ppm (v1: 2.7); z = 3.5 is 465 ppm (v1: 0.032). The table now also gives one-sided ppm and z = 4.0.

**5. Optimisation.** The model is now solved. For fixed breakpoints the minimax problem is a linear programme (HiGHS via SciPy) with y₀ = 0.5 and 0 ≤ yⱼ ≤ 1. Free breakpoints use a multistart Nelder-Mead outer search, reported as "best found", not as a proven global optimum. A simplicity-constrained version (two-decimal node values, half-step breakpoints) is solved by exhaustive search. "Pareto frontier" now refers only to Table 6 and Figure 3, with an explicit complexity measure (number of stored constants), and the text states that the measure ignores arithmetic difficulty. A notable result: per stored constant, equal-width breakpoints with minimax values dominate free breakpoints.

**6. Engineering suitability.** The claim "suitable for quality control" is removed. Section 8 states the bound (|e| ≤ 0.0235 for Model C, so ±4.7 pp in two-sided coverage), lists suitable uses (rough yield screening, software sanity checks, teaching), and lists unsuitable ones (ppm, capability claims near limits, safety-critical reliability, z > 3). Section 6.4 and Figure 4 quantify the tail effect: at z = 2.5, Model A overestimates the defect rate by 94%.

**7. Model A vs Model B boundaries.** Model A is continuous by construction. Table 2 gives Model B's value, right limit, jump, and error at z = 1 (+0.001 jump) and z = 2 (−0.004 jump, overestimate of +0.0038 on the left). v1 stated jumps of 0.0003 and 0.0002, which was wrong. A continuous two-digit Model C replaces Model B as the recommended field rule.

**8. Validation statistics.** All recomputed on a 30,001-point grid. The reference CDF (`scipy.special.ndtr`) is checked against 50-digit `mpmath` (max difference 1.1 × 10⁻¹⁶). Environment versions are stored in `results.json`. The v1 statistics (attributed to "R 4.1.0") had not been computed and are withdrawn. Key corrections: Model A's maximum error is 0.0240 at z = 1.47, not 0.0215 at z = 0.5; its mean signed error is −0.0112, not "≈ 0". `tests/test_models.py` checks the main claims, and CI reruns the full analysis on every push.

## Earlier major comments

| Comment | Change |
|---|---|
| Model not continuous as written | Models A, B, C defined separately; continuity and jumps stated (Sections 4.1 to 4.3) |
| Exact secant vs rounded model | Separate definitions, separate statistics (Tables 3 and 5) |
| "Less than 2%" claim | Removed; metrics defined in Section 3 |
| Concavity argument wrong | Rewritten: Φ''(z) = −zφ(z) < 0 for z > 0, so Model A underestimates inside every segment; confirmed numerically |
| One-sided vs two-sided examples | All examples use 2Φ(z) − 1 for centred two-sided coverage (Section 7) |
| "Immaterial beyond z = 3" | Removed; replaced with the suggested restriction statement |
| Novelty claim | Reframed around the minimax formulation, frontier, and tail analysis; interpolation itself is described as classical |
| Six-segment 0.0005 claim | Computed value is 0.0071; claim withdrawn (Section 6.2) |
| Citations | Each reference checked against the claim it supports; Chase et al. (2011) replaced by Chase and Parkinson (1991); Bowling et al. (2009), Pólya (1949), de Boor, HiGHS and Nelder-Mead added |
