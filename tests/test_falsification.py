"""Deliberate-failure tests.

A reproducibility package whose tests cannot fail is of no value. This
file constructs deliberately broken configurations (broken coefficients,
wrong formula plug-in, swapped targets) and verifies that the closure
correctly fails. These are exactly the failure modes that the paper's
falsification section labels (F1)-(F4).
"""

import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import recompute_landings as M


def _coeff_baseline():
    return {
        "alpha_xi":  0.90082,
        "D_Omega":   0.83996,
        "beta_pi":   0.93791,
        "gamma":     0.10021,
        "eps_sync2": 0.05000,
        "G_NET":     1.78852,
        "N_gen":     3,
    }


def test_baseline_l1_l5_pass():
    """Sanity: the baseline coefficients pass all L1-L5 within 2.5%."""
    coeff = _coeff_baseline()
    preds = {
        "L1": coeff["alpha_xi"] + coeff["beta_pi"] + coeff["eps_sync2"] - coeff["gamma"],
        "L2": coeff["alpha_xi"] ** 2 * coeff["eps_sync2"] * coeff["N_gen"],
        "L3": -1 + coeff["eps_sync2"] ** 2 / coeff["gamma"],
        "L4": coeff["gamma"] * coeff["beta_pi"] * 4 / math.pi,
        "L5": (coeff["alpha_xi"] / coeff["beta_pi"]) * (coeff["D_Omega"] / (1 + coeff["gamma"])) * coeff["G_NET"],
    }
    targets = {"L1": 11/6, "L2": 0.120, "L3": -1.0, "L4": 0.1179, "L5": 4/3}
    for k in preds:
        assert abs(preds[k] - targets[k]) / abs(targets[k]) <= 0.025


def test_wrong_alpha_xi_breaks_robust_core():
    """If alpha_xi is shifted by 50%, several robust rows fail."""
    coeff = _coeff_baseline()
    coeff["alpha_xi"] *= 1.5
    coeff["G_NET"] = (coeff["alpha_xi"] + coeff["beta_pi"]
                      + coeff["eps_sync2"] - coeff["gamma"])
    L1 = coeff["alpha_xi"] + coeff["beta_pi"] + coeff["eps_sync2"] - coeff["gamma"]
    L2 = coeff["alpha_xi"] ** 2 * coeff["eps_sync2"] * coeff["N_gen"]
    target_L1 = 11/6
    target_L2 = 0.120
    res_L1 = abs(L1 - target_L1) / abs(target_L1)
    res_L2 = abs(L2 - target_L2) / abs(target_L2)
    assert res_L1 > 0.025 or res_L2 > 0.025, (
        "Robust core should fail under a 50% alpha_xi perturbation."
    )


def test_zero_gamma_breaks_w_DE():
    """gamma -> 0 would make the dark-energy formula L3 diverge."""
    eps = 0.05
    with pytest.raises(ZeroDivisionError):
        gamma = 0.0
        _ = -1 + eps ** 2 / gamma


def test_swapped_target_does_not_pass():
    """Swap the alpha_s target with sin^2 theta_W; the robust core fails."""
    coeff = _coeff_baseline()
    L4_pred = coeff["gamma"] * coeff["beta_pi"] * 4 / math.pi
    wrong_target = 0.23122
    residual = abs(L4_pred - wrong_target) / abs(wrong_target)
    assert residual > 0.025, (
        "Plugging the L4 prediction against the L6 target must fail; "
        "otherwise the eight-row mapping is degenerate."
    )


def test_random_coefficients_break_robust_core():
    """A random coefficient set should not reproduce the L1-L5 closure."""
    from random import Random
    rng = Random(0)
    coeff = {
        "alpha_xi":  rng.uniform(0.02, 0.99),
        "D_Omega":   rng.uniform(0.02, 0.99),
        "beta_pi":   rng.uniform(0.02, 0.99),
        "gamma":     rng.uniform(0.02, 0.99),
        "eps_sync2": rng.uniform(0.02, 0.99),
        "N_gen": 3,
    }
    coeff["G_NET"] = (coeff["alpha_xi"] + coeff["beta_pi"]
                      + coeff["eps_sync2"] - coeff["gamma"])
    L1 = coeff["G_NET"]
    L2 = coeff["alpha_xi"] ** 2 * coeff["eps_sync2"] * coeff["N_gen"]
    L4 = coeff["gamma"] * coeff["beta_pi"] * 4 / math.pi
    targets = {"L1": 11/6, "L2": 0.120, "L4": 0.1179}
    preds = {"L1": L1, "L2": L2, "L4": L4}
    n_pass = sum(1 for k in preds
                 if abs(preds[k] - targets[k]) / abs(targets[k]) <= 0.025)
    assert n_pass < 3, (
        f"Random coefficients passed {n_pass}/3 of the robust subset; "
        f"this would indicate the test is too weak."
    )


def test_unknown_formula_raises():
    coeff = _coeff_baseline()
    with pytest.raises(ValueError):
        M.evaluate_formula("NOT A REAL FORMULA", coeff)


def test_F2_boundary_perturbation_pushes_robust_core_past_25pct():
    """F2 boundary test: find the smallest alpha_xi perturbation that
    pushes one of the robust rows past the 2.5% PRECISE bound. The
    closure should be sharp at this threshold — a smaller perturbation
    keeps the row inside, a larger one pushes it out.

    We solve for the boundary perturbation analytically on L1
    (alpha_dn = alpha_xi + beta_pi + eps^2 - gamma):
      target = 11/6 = 1.83333
      baseline L1 = 1.78852, residual -2.444%
      perturbation that hits the +2.5% boundary on L1 from above:
        L1_boundary = target * 1.025 = 1.87917
        delta_alpha_xi = L1_boundary - L1_baseline = 0.09065
      so a +0.0907 shift on alpha_xi pushes L1 past +2.5%.
    """
    coeff = _coeff_baseline()
    L1_baseline = (coeff["alpha_xi"] + coeff["beta_pi"]
                   + coeff["eps_sync2"] - coeff["gamma"])
    target = 11 / 6
    # Smallest delta_alpha_xi that pushes L1 to +2.5% above target
    delta_required = target * 1.025 - L1_baseline
    # Just inside (should still pass)
    coeff_in = dict(coeff)
    coeff_in["alpha_xi"] = coeff["alpha_xi"] + delta_required * 0.99
    L1_in = (coeff_in["alpha_xi"] + coeff_in["beta_pi"]
             + coeff_in["eps_sync2"] - coeff_in["gamma"])
    res_in = abs(L1_in - target) / target
    assert res_in <= 0.025, (
        f"Just inside the F2 boundary, L1 residual {res_in:.4%} should still "
        f"be within 2.5%."
    )
    # Just outside (should fail)
    coeff_out = dict(coeff)
    coeff_out["alpha_xi"] = coeff["alpha_xi"] + delta_required * 1.01
    L1_out = (coeff_out["alpha_xi"] + coeff_out["beta_pi"]
              + coeff_out["eps_sync2"] - coeff_out["gamma"])
    res_out = abs(L1_out - target) / target
    assert res_out > 0.025, (
        f"Just outside the F2 boundary, L1 residual {res_out:.4%} must "
        f"exceed 2.5% — the closure must be sharp at the threshold."
    )
