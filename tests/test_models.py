"""Fast checks of the key claims in the paper. Run: pytest -q"""
import json, os
import numpy as np
from scipy.special import ndtr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, "results", "results.json")))
Z = np.linspace(0, 3, 30001)
PHI = ndtr(Z)
A = np.interp(Z, [0, 1, 2, 3], ndtr(np.array([0, 1, 2, 3.0])))
C = np.interp(Z, [0, 1, 2, 3], [0.5, 0.84, 0.98, 1.0])


def B(z):
    z = np.asarray(z, float)
    return np.where(z <= 1, 0.5 + 0.34 * z, np.where(z <= 2, 0.841 + 0.14 * (z - 1), 0.977 + 0.02 * (z - 2)))


def test_model_A_never_overestimates():  # concavity of Phi on z > 0
    assert np.all(A - PHI <= 1e-15)


def test_model_A_max_error_location():
    e = np.abs(A - PHI)
    assert abs(e.max() - 0.0240) < 5e-5
    assert 1.44 < Z[e.argmax()] < 1.50


def test_model_A_within_interpolation_bound():
    assert np.abs(A - PHI).max() <= 0.2419707245 / 8  # h^2/8 * phi(1), h = 1


def test_model_B_jumps():
    assert abs(float(B(1.0 + 1e-12)) - float(B(1.0)) - 0.001) < 1e-9
    assert abs(float(B(2.0 + 1e-12)) - float(B(2.0)) + 0.004) < 1e-9
    assert float(B(2.0)) - ndtr(2.0) > 0.0037  # overestimates left of z = 2


def test_model_C_continuous_and_max():
    assert abs(np.abs(C - PHI).max() - 0.0235) < 5e-5


def test_two_sided_ppm_table():
    expect = {1.0: 317311, 1.5: 133614, 2.0: 45500, 2.5: 12419, 3.0: 2700, 3.5: 465, 4.0: 63}
    for z, ppm in expect.items():
        assert round(2 * (1 - ndtr(z)) * 1e6) == ppm


def test_results_json_matches_recomputation():
    assert abs(R["models"]["A"]["max_abs_error"] - np.abs(A - PHI).max()) < 1e-6
    assert abs(R["models"]["C"]["max_abs_error"] - np.abs(C - PHI).max()) < 1e-6
    assert R["models"]["A"]["max_positive_error"] == 0.0


def test_optimised_models_are_valid_cdfs():
    for k in "DEF":
        y = np.array(R[f"model_{k}"]["values"])
        assert y[0] == 0.5 and np.all(np.diff(y) >= 0) and np.all(y <= 1)
    assert R["models"]["D"]["max_abs_error"] < R["models"]["A"]["max_abs_error"]
    assert R["models"]["E"]["max_abs_error"] < R["models"]["D"]["max_abs_error"]


def test_polya_max_error():
    P = 0.5 + 0.5 * np.sqrt(1 - np.exp(-2 * Z ** 2 / np.pi))
    assert abs(np.abs(P - PHI).max() - 0.00315) < 5e-5
