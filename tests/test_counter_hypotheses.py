"""Tests for the CH-1..CH-5 counter-hypothesis matrix (Section 5).

Each counter-hypothesis is paired with a specific bundled numerical
quantity that operationally distinguishes it. These tests assert that
each evidence quantity is present, has a sensible numerical value, and
carries the recorded verdict.
"""

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import verify_counter_hypotheses as M  # noqa: E402


@pytest.fixture(scope="module")
def evidence():
    return M.load_evidence()


@pytest.fixture(scope="module")
def output(evidence):
    M.main()
    out_path = REPO / "outputs" / "counter_hypotheses_summary.json"
    with open(out_path, "r", encoding="utf-8") as f:
        return json.load(f)


def test_all_five_counter_hypotheses_present(evidence):
    chs = evidence["counter_hypotheses"]
    assert len(chs) == 5
    ids = [c["id"] for c in chs]
    assert ids == ["CH-1", "CH-2", "CH-3", "CH-4", "CH-5"]


def test_carrier_to_arrow_S_bounce_strong(evidence):
    arrow = evidence["carrier_to_arrow"]
    assert arrow["S_bounce_canonical"] > 30.0
    assert arrow["time_asymmetry_ratio_canonical"] > 1e30
    assert arrow["phase_is_bounce_carrier"] is True
    assert arrow["strength"] == "STRONG"


def test_carrier_to_world_admissibility_six_of_seven(evidence):
    w = evidence["carrier_to_world_admissibility"]
    assert w["omega_world_axis_coverage_passed"] == 6
    assert w["omega_world_axis_coverage_total"] == 7
    assert w["carrier_threshold_world_guard_canonical"] == "WORLD_OK"


def test_carrier_to_irreversibility_gamow_dominant(evidence):
    i = evidence["carrier_to_irreversibility"]
    assert i["gamow_channel_share_canonical_pct"] >= 90.0


def test_carrier_to_geometry_quantitative(evidence):
    g = evidence["carrier_to_geometry"]
    assert 0.0 < g["metric_quality_score_canonical"] < 1.0
    assert abs(g["sin_theta_W_squared_residual_pct"]) < 0.01
    assert g["G_N_ratio"] == pytest.approx(1.0, abs=0.001)


def test_carrier_to_directed_continuation_temporal_admissibility(evidence):
    dc = evidence["carrier_to_directed_continuation"]
    assert dc["temporal_admissibility_axes_passed"] == \
           dc["temporal_admissibility_axes_total"]
    assert dc["S_bounce_canonical"] > 30.0


def test_all_counter_hypotheses_weakened(evidence):
    """Every CH-i records a weakened or strongly-weakened verdict
    against the bundled evidence."""
    for ch in evidence["counter_hypotheses"]:
        assert ch["verdict"] in ("weakened", "strongly_weakened"), (
            f"{ch['id']} has unexpected verdict {ch['verdict']!r}"
        )


def test_CH_4_is_strongly_weakened(evidence):
    """CH-4 (carrier is auxiliary) is the strongest counter-hypothesis
    a priori and is strongly weakened by the simultaneous derivation
    of the emergent metric, the Newtonian Schwarzschild far field, and
    the sin^2 theta_W reproduction."""
    ch4 = next(c for c in evidence["counter_hypotheses"]
               if c["id"] == "CH-4")
    assert ch4["verdict"] == "strongly_weakened"


def test_evidence_keys_resolve_for_each_CH(evidence):
    """Each CH lists evidence keys; each must resolve to a bundled value."""
    for ch in evidence["counter_hypotheses"]:
        for k in ch["evidence_against_keys"]:
            v = M.resolve_evidence_value(evidence, k)
            assert v is not None


def test_summary_output_aggregate(output):
    assert output["n_total"] == 5
    assert output["all_weakened"] is True
    assert output["n_strongly_weakened"] == 1
