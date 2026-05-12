"""Tests for the bundled Yukawa-cluster Fisher's combined-p statistic."""

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def fisher():
    with open(REPO / "data" / "yukawa_cluster_fisher.json", "r", encoding="utf-8") as f:
        return json.load(f)


def test_three_observables_present(fisher):
    obs = fisher["observables"]
    assert len(obs) == 3
    names = {o["name"] for o in obs}
    assert names == {"alpha_dn", "w_DE", "H_0"}


def test_two_distinct_sectors(fisher):
    obs = fisher["observables"]
    sectors = {o["sector"].split()[0] for o in obs}
    assert "QFT" in sectors
    assert "Einstein" in sectors


def test_anchor_residuals_below_zero_point_three_pct(fisher):
    """All three Yukawa-cluster observables sub-0.3% on the anchor target."""
    for o in fisher["observables"]:
        assert o["anchor_residual_pct"] < 0.3


def test_fishers_combined_p_recomputes(fisher):
    """Fisher's T = -2 * sum(ln p_i) on the three bundled per-observable
    p-values must match the bundled T to four significant figures."""
    import math
    pvals = [o["per_observable_p_value"] for o in fisher["observables"]]
    T_recompute = -2.0 * sum(math.log(p) for p in pvals)
    T_bundle = fisher["fishers_combined_statistic"]["T"]
    assert abs(T_recompute - T_bundle) / T_bundle < 1e-3, (
        f"Fisher's T mismatch: bundle={T_bundle}, recompute={T_recompute:.4f}"
    )


def test_combined_p_below_three_sigma(fisher):
    """The Fisher's-combined p must be below 2.7e-3 (3 sigma)."""
    p = fisher["fishers_combined_statistic"]["p_combined"]
    assert p < 2.7e-3


def test_significance_above_four_sigma(fisher):
    sigma = fisher["fishers_combined_statistic"]["approximate_significance_sigma"]
    assert sigma >= 4.0


def test_loop_class_factor(fisher):
    assert fisher["loop_class"] == "1+gamma/4"
    assert fisher["lemma"] == 1
