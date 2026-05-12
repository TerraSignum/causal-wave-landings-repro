"""Tests for the IR-hyperbolicity certificate of the Xi-d'Alembertian.

The IR-hyperbolicity criterion (Theorem 7.9 of the program's proof
collection) requires four moment conditions on the retarded-shell
stencil: vanishing linear spatial and mixed time-space moments, and
strictly positive time-normalisation Z_tau and spatial variance Z_x.
Together they yield the IR symbol expansion
    sigma_Xi = -Z_tau omega^2 + Z_x |k|^2 + m_eff^2
              + O(|omega|^3 + |k|^3),
emergent signal speed c_Xi^2 = Z_x / Z_tau, and Lorentz-deviation
Delta_Lor(k) = o(k^2) at small k. This test asserts each of the four
conditions holds on the bundled 1+1d symmetric retarded-shell stencil
and that the o(k^2) statement is reproduced numerically.
"""

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import verify_ir_hyperbolicity as M  # noqa: E402


@pytest.fixture(scope="module")
def output():
    M.main()
    out_path = REPO / "outputs" / "ir_hyperbolicity_certificate.json"
    with open(out_path, "r", encoding="utf-8") as f:
        return json.load(f)


def test_M1_linear_spatial_moment_vanishes():
    """M1: sum w * dx = 0 on a symmetric stencil."""
    assert M.linear_spatial_moment(M.STENCIL) == pytest.approx(0.0, abs=1e-12)


def test_M2_mixed_time_space_moment_vanishes():
    """M2: sum w * dt * dx = 0 on a symmetric stencil."""
    assert M.mixed_time_space_moment(M.STENCIL) == pytest.approx(0.0, abs=1e-12)


def test_M3_time_normalisation_strictly_positive():
    """M3: Z_tau > 0 (required for hyperbolicity)."""
    assert M.time_normalisation(M.STENCIL) > 0.0


def test_M4_spatial_variance_strictly_positive():
    """M4: Z_x > 0 (required for emergent positive c_Xi^2)."""
    assert M.spatial_variance(M.STENCIL) > 0.0


def test_emergent_signal_speed_strictly_positive():
    """c_Xi^2 = Z_x / Z_tau > 0 -> hyperbolic IR sector."""
    Z_tau = M.time_normalisation(M.STENCIL)
    Z_x = M.spatial_variance(M.STENCIL)
    c_Xi_sq = Z_x / Z_tau
    assert c_Xi_sq > 0.0
    assert c_Xi_sq == pytest.approx(1.0, abs=1e-12)


def test_lorentz_deviation_is_o_k_squared():
    """Delta_Lor(k) / k^2 -> 0 as k -> 0 (the o(k^2) statement).

    With the bundled stencil, Delta_Lor(k) / k^2 = c_Xi^2 k^2 / 2 to
    leading order, so the ratio of consecutive small-k samples should
    scale as k^2. Equivalently, the ratio at k=0.01 must be at least
    100x smaller than at k=0.1.
    """
    deviations, c_Xi_sq = M.lorentz_deviation(M.STENCIL,
                                               [0.01, 0.1])
    ratio_small_k = deviations[0][1]  # k=0.01
    ratio_large_k = deviations[1][1]  # k=0.1
    assert ratio_small_k > 0
    assert ratio_large_k > 0
    # The o(k^2) statement: ratio at k=0.01 must be << ratio at k=0.1
    assert ratio_small_k < ratio_large_k / 50.0, (
        f"Lorentz-deviation does not vanish faster than k^2: "
        f"{ratio_small_k} vs {ratio_large_k}"
    )
    # Specifically, leading order is c_Xi^2 k^2 / 2 -> ratio scales as k^2:
    # ratio(0.01)/ratio(0.1) should be 1/100 to leading order.
    scale_factor = ratio_small_k / ratio_large_k
    assert scale_factor == pytest.approx(0.01, rel=0.05)


def test_certificate_output_has_expected_keys(output):
    expected = {
        "criterion", "theorem", "stencil", "moments",
        "effective_mass_squared", "emergent_signal_speed_squared",
        "hyperbolicity_passes", "lorentz_deviation_sweep",
    }
    assert expected.issubset(set(output.keys()))
    assert output["hyperbolicity_passes"] is True
    for m in ("M1_linear_spatial", "M2_mixed_time_space",
              "M3_time_normalisation", "M4_spatial_variance"):
        assert output["moments"][m]["passes"] is True


def test_emergent_signal_speed_in_certificate(output):
    assert output["emergent_signal_speed_squared"] == pytest.approx(1.0, abs=1e-12)
