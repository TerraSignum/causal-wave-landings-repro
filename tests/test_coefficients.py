"""Verify that the five measured causal-wave coefficients load correctly
and have the expected values from the frozen JSON input.
"""

import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import recompute_landings as M


def test_five_coefficients_present():
    coeff = M.load_coefficients()
    for k in ("alpha_xi", "D_Omega", "beta_pi", "gamma", "eps_sync2"):
        assert k in coeff
        assert isinstance(coeff[k], float)


def test_alpha_xi_value():
    coeff = M.load_coefficients()
    assert coeff["alpha_xi"] == pytest.approx(0.90082, abs=1e-5)


def test_D_Omega_value():
    coeff = M.load_coefficients()
    assert coeff["D_Omega"] == pytest.approx(0.83996, abs=1e-5)


def test_beta_pi_value():
    coeff = M.load_coefficients()
    assert coeff["beta_pi"] == pytest.approx(0.93791, abs=1e-5)


def test_gamma_value():
    coeff = M.load_coefficients()
    assert coeff["gamma"] == pytest.approx(0.10021, abs=1e-5)


def test_eps_sync2_value():
    coeff = M.load_coefficients()
    assert coeff["eps_sync2"] == pytest.approx(0.05000, abs=1e-5)


def test_G_NET_consistency():
    """G_NET = alpha_xi + beta_pi + eps_sync2 - gamma."""
    coeff = M.load_coefficients()
    g_recompute = (coeff["alpha_xi"] + coeff["beta_pi"]
                   + coeff["eps_sync2"] - coeff["gamma"])
    assert g_recompute == pytest.approx(coeff["G_NET"], abs=1e-5)
    assert coeff["G_NET"] == pytest.approx(1.78852, abs=1e-5)


def test_N_gen_is_three():
    coeff = M.load_coefficients()
    assert coeff["N_gen"] == 3


def test_coefficients_are_in_unit_interval():
    """All five coefficients lie in (0, 1) by physical interpretation."""
    coeff = M.load_coefficients()
    for k in ("alpha_xi", "D_Omega", "beta_pi", "gamma", "eps_sync2"):
        assert 0.0 < coeff[k] < 1.0
