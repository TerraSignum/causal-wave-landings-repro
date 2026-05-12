"""Verify that the eight benchmark landings reproduce within tolerance.

The structural claim of Paper 2 is the PRECISE_2.5-or-better closure on
all eight landings under the five measured coefficients without fitting.
This test enforces that claim mechanically.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import recompute_landings as M


def _rows():
    return M.compute_all()["rows"]


def test_all_ten_landings_present():
    rows = _rows()
    assert len(rows) == 10
    assert {r["id"] for r in rows} == {f"L{i}" for i in range(1, 11)}


def test_all_within_2_5_percent():
    """Hauptclaim: 10/10 within 2.5% residual."""
    rows = _rows()
    for r in rows:
        assert abs(r["residual_pct"]) <= 2.5, (
            f"{r['id']} ({r['observable']}) residual "
            f"{r['residual_pct']:+.4f}% exceeds 2.5%."
        )


def test_robust_core_under_2_5_percent():
    """Robust core L1-L5 must be within 2.5% (PRECISE_2.5 tier)."""
    rows = {r["id"]: r for r in _rows()}
    for rid in ("L1", "L2", "L3", "L4", "L5"):
        assert abs(rows[rid]["residual_pct"]) <= 2.5


def test_three_exact_tier_targets_are_under_0_01_percent():
    """L6, L7, L8 must indeed be within 0.01% (EXACT tier classification)."""
    rows = {r["id"]: r for r in _rows()}
    for rid in ("L6", "L7", "L8"):
        assert abs(rows[rid]["residual_pct"]) <= 0.01


def test_L1_alpha_dn_yukawa():
    rows = {r["id"]: r for r in _rows()}
    r = rows["L1"]
    assert r["target"] == pytest.approx(11 / 6, abs=1e-9)
    assert r["prediction"] == pytest.approx(1.78852, abs=1e-4)


def test_L2_omega_dm_h2():
    rows = {r["id"]: r for r in _rows()}
    r = rows["L2"]
    assert r["target"] == pytest.approx(0.120, abs=1e-9)
    assert r["prediction"] == pytest.approx(0.12172, abs=1e-4)


def test_L4_alpha_s_M_Z():
    rows = {r["id"]: r for r in _rows()}
    r = rows["L4"]
    assert r["target"] == pytest.approx(0.1179, abs=1e-9)
    assert r["prediction"] == pytest.approx(0.11967, abs=1e-4)


def test_L6_sin2_theta_W_exact_tier():
    rows = {r["id"]: r for r in _rows()}
    r = rows["L6"]
    assert r["prediction"] == pytest.approx(0.23122, abs=1e-4)
    assert abs(r["residual_pct"]) <= 0.01


def test_L7_BH_quarter_exact_tier():
    rows = {r["id"]: r for r in _rows()}
    r = rows["L7"]
    assert r["prediction"] == pytest.approx(0.25, abs=1e-3)
    assert abs(r["residual_pct"]) <= 0.01


def test_L8_einstein_gap_exact_tier():
    rows = {r["id"]: r for r in _rows()}
    r = rows["L8"]
    assert r["prediction"] == pytest.approx(2 / 3, abs=1e-3)
    assert abs(r["residual_pct"]) <= 0.01


def test_max_residual_is_below_2_5():
    rows = _rows()
    max_resid = max(abs(r["residual_pct"]) for r in rows)
    assert max_resid <= 2.5
