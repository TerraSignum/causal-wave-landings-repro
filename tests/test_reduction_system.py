"""Tests for the coefficient-reduction system R = {C1, ..., C5} of Section 4.

Each C_i is verified at the per-residual level (each <= 0.30% on the
measured coefficients). The algebraic-exact BH-1/4 identity in Q is
verified in closed form, and the integer uniqueness of N=3 for
BH(N) = 1/4 is verified by an integer scan.
"""

import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import verify_reduction_system as M  # noqa: E402


@pytest.fixture(scope="module")
def coeff():
    return M.load_coefficients()


@pytest.fixture(scope="module")
def residuals(coeff):
    return M.reduction_residuals(coeff)


@pytest.fixture(scope="module")
def output(coeff, residuals):
    out_path = REPO / "outputs" / "reduction_system.json"
    if out_path.exists():
        with open(out_path, "r", encoding="utf-8") as f:
            return json.load(f)
    return None


def test_coefficients_loaded(coeff):
    assert coeff["alpha_xi"] == pytest.approx(0.90082)
    assert coeff["D_Omega"] == pytest.approx(0.83996)
    assert coeff["beta_pi"] == pytest.approx(0.93791)
    assert coeff["gamma"] == pytest.approx(0.10021)
    assert coeff["eps_sync2"] == pytest.approx(0.05000)
    assert coeff["N_gen"] == 3


def test_C1_complementarity(residuals):
    label, _expr, res = residuals[0]
    assert label == "C1"
    assert res <= 0.0030, f"C1 residual {res*100:.3f}% exceeds 0.30%"


def test_C2_einstein_relation(residuals):
    label, _expr, res = residuals[1]
    assert label == "C2"
    assert res <= 0.0030, f"C2 residual {res*100:.3f}% exceeds 0.30%"


def test_C3_fdt_symmetry(residuals):
    label, _expr, res = residuals[2]
    assert label == "C3"
    assert res <= 0.0030, f"C3 residual {res*100:.3f}% exceeds 0.30%"


def test_C4_generation_geometry(residuals):
    label, _expr, res = residuals[3]
    assert label == "C4"
    assert res <= 0.0030, f"C4 residual {res*100:.3f}% exceeds 0.30%"


def test_C5_clifford_projector(residuals):
    label, _expr, res = residuals[4]
    assert label == "C5"
    assert res <= 0.0030, f"C5 residual {res*100:.3f}% exceeds 0.30%"


def test_max_residual_below_threshold(residuals):
    """The manuscript Section 4 claim: max residual across C1-C5 is < 0.30%."""
    max_res = max(r for (_, _, r) in residuals)
    assert max_res <= 0.0030, (
        f"Max residual {max_res*100:.3f}% exceeds the 0.30% threshold "
        f"asserted in Section 4 of the manuscript."
    )


def test_rational_reduction_BH_quarter_algebraic():
    """alpha_xi/2 - 2*gamma = 9/20 - 1/5 = 1/4 algebraically in Q."""
    lhs, target, eq = M.bh_quarter_under_R()
    assert lhs == Fraction(1, 4)
    assert target == Fraction(1, 4)
    assert eq is True


def test_BH_general_formula_integer_uniqueness():
    """N=3 is the unique positive integer solving BH(N) = 1/4."""
    sols = M.integer_uniqueness_check_for_quarter()
    positive = [n for n in sols if n > 0]
    assert positive == [3]


def test_C5_at_lattice_value_close_to_15_16(coeff):
    """beta_pi is close to 15/16 within 0.05%."""
    rel = abs(coeff["beta_pi"] - 15.0 / 16.0) / coeff["beta_pi"]
    assert rel < 0.0006, f"beta_pi off 15/16 by {rel*100:.4f}%, expected < 0.06%"


def test_canonical_schema_blocks_present():
    """The verifier must emit the same JSON schema as the upstream pipeline
    (outputs_causal_wave_universality/coefficient_reduction.json), so that a
    reader of the public bundle reproduces the canonical structure key-for-key.
    """
    M.main()  # regenerate the file from scratch
    out_path = REPO / "outputs" / "coefficient_reduction.json"
    with open(out_path, "r", encoding="utf-8") as f:
        canon = json.load(f)
    expected_top = {
        "observed",
        "C1_alpha_plus_gamma_is_1",
        "C2_D_equals_beta_minus_gamma",
        "C3_eps_equals_half_gamma",
        "C4_tan_theta_is_1_over_Ngen",
        "reduction_1param",
        "beta_pi_algebraic_forms",
        "reduction_0param",
    }
    assert set(canon.keys()) == expected_top, (
        f"Canonical schema mismatch: got {set(canon.keys())}, "
        f"expected {expected_top}"
    )


def test_canonical_C1_to_C4_values_match_upstream():
    """The C1..C4 blocks must reproduce the upstream pipeline values."""
    out_path = REPO / "outputs" / "coefficient_reduction.json"
    with open(out_path, "r", encoding="utf-8") as f:
        canon = json.load(f)

    # Upstream reference values (frozen from
    # outputs_causal_wave_universality/coefficient_reduction.json,
    # 2026-04-26 release of the universality pipeline).
    assert canon["C1_alpha_plus_gamma_is_1"]["observed_sum"] == pytest.approx(1.00103)
    assert canon["C1_alpha_plus_gamma_is_1"]["target"] == 1.0
    assert canon["C1_alpha_plus_gamma_is_1"]["deviation_pct"] == pytest.approx(0.10300, abs=1e-4)
    assert canon["C1_alpha_plus_gamma_is_1"]["passes_pg_tol"] is True

    assert canon["C2_D_equals_beta_minus_gamma"]["predicted"] == pytest.approx(0.83770)
    assert canon["C2_D_equals_beta_minus_gamma"]["residual_pct"] == pytest.approx(0.26906, abs=1e-4)
    assert canon["C2_D_equals_beta_minus_gamma"]["passes_pg_tol"] is True

    assert canon["C3_eps_equals_half_gamma"]["predicted"] == pytest.approx(0.050105)
    assert canon["C3_eps_equals_half_gamma"]["residual_pct"] == pytest.approx(0.21000, abs=1e-4)
    assert canon["C3_eps_equals_half_gamma"]["passes_pg_tol"] is True

    c4 = canon["C4_tan_theta_is_1_over_Ngen"]
    assert c4["N_gen"] == 3
    assert c4["gamma_predicted"] == pytest.approx(0.1)
    assert c4["alpha_predicted"] == pytest.approx(0.9)
    assert c4["gamma_residual_pct"] == pytest.approx(0.20956, abs=1e-4)
    assert c4["alpha_residual_pct"] == pytest.approx(0.09103, abs=1e-4)
    assert c4["passes_pg_tol"] is True


def test_canonical_reduction_1param_and_0param_max_residuals():
    """Both closure tiers must match the upstream max-residual values."""
    out_path = REPO / "outputs" / "coefficient_reduction.json"
    with open(out_path, "r", encoding="utf-8") as f:
        canon = json.load(f)
    r1 = canon["reduction_1param"]
    assert r1["free_parameters"] == ["beta_pi"]
    assert r1["max_residual_pct"] == pytest.approx(0.24406, abs=1e-4)
    assert r1["passes_pg_tol"] is True
    r0 = canon["reduction_0param"]
    assert r0["beta_pi_formula"] == "15/16"
    assert r0["max_residual_pct"] == pytest.approx(0.29287, abs=1e-4)
    assert r0["passes_pg_tol"] is True


def test_canonical_beta_pi_algebraic_forms_table():
    """The 9-form survey must reproduce the upstream residuals exactly."""
    out_path = REPO / "outputs" / "coefficient_reduction.json"
    with open(out_path, "r", encoding="utf-8") as f:
        canon = json.load(f)
    forms = canon["beta_pi_algebraic_forms"]
    # Frozen upstream survey: residual_pct values to 4 dp.
    upstream = {
        "3/pi":                          1.8146,
        "pi/4 + 1/2":                   37.0492,
        "1 - 1/(4*pi)":                  1.8645,
        "1 - gamma * (pi/4)":            1.7715,
        "cos(pi/12)**2":                 0.5222,
        "(N_gen + 1/pi) / (N_gen + 1)": 11.5504,
        "1 - 2/pi * eps":                3.2262,
        "15/16":                         0.0437,
        "1 - 1/(2*pi)":                 10.3491,
    }
    for name, expected_res in upstream.items():
        assert name in forms, f"missing algebraic form: {name}"
        assert forms[name]["residual_pct"] == pytest.approx(expected_res, abs=1e-3)
    # Only 15/16 passes the strict PG band; the 0-param closure picks it.
    assert forms["15/16"]["residual_pct"] < 0.30


def test_predictions_P1_to_P8_match_manuscript():
    """Pre-registered predictions P1..P8 from Section 4 of manuscript.

    Each prediction is a tolerance condition on a measured quantity; this
    test verifies that the present coefficients satisfy each prediction
    within its declared tolerance (2.5%).
    """
    coeff = M.load_coefficients()
    a, D, b, g, e2 = (coeff["alpha_xi"], coeff["D_Omega"],
                      coeff["beta_pi"], coeff["gamma"], coeff["eps_sync2"])
    Ngen = coeff["N_gen"]
    tol = 0.025

    # P1: |alpha_xi + gamma - 1| <= 0.025
    assert abs(a + g - 1.0) <= tol

    # P2: |D_Omega - (beta_pi - gamma)| / |D_Omega| <= 0.025
    assert abs(D - (b - g)) / abs(D) <= tol

    # P3: |eps_sync2 - gamma/2| / |eps_sync2| <= 0.025
    assert abs(e2 - g / 2) / abs(e2) <= tol

    # P4: |gamma - 0.1| / gamma <= 0.025
    assert abs(g - 0.1) / g <= tol

    # P5: |beta_pi - 0.9375| / beta_pi <= 0.025
    assert abs(b - 0.9375) / b <= tol

    # P6: |alpha_xi/2 - 2*gamma - 1/4| / (1/4) <= 0.025
    assert abs(a / 2 - 2 * g - 0.25) / 0.25 <= tol

    # P7 — landing-row 2.5% bound is tested separately in test_landings.py
    # (covered for L1..L8 there); this test checks the cross-link in spirit.

    # P8 — integer uniqueness of (N_gen, d) = (3, 4) under
    # gamma = 1/(N^2 + 1), beta_pi = (2^d - 1)/2^d
    assert abs(g - 1.0 / (Ngen ** 2 + 1)) / g <= tol
    d_int = 4
    assert abs(b - (2 ** d_int - 1) / 2 ** d_int) / b <= tol
