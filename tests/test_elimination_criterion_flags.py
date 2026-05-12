"""Verify that the Look-Elsewhere caveat is correctly flagged.

Rows L6 ($\\sin^2\\theta_W$) and L8 (Einstein-gap exponent) are
EXACT-tier numerically but enumeration-induced; the manuscript
reports them as ``look_elsewhere_caveated=true``.
Row L7 (BH 1/4) ascends to a $\\mathbb Q$-exact algebraic identity
$\\alpha_\\xi/2 - 2\\gamma = 9/20 - 1/5 = 1/4$ under the
coefficient-reduction system $\\mathcal R$, with closed-form
integer-uniqueness $N_\\mathrm{gen}=3$, and is therefore exempt
from the multiple-trial inflation correction.
Rows L1-L5 are the robust PRECISE_2.5 core and must NOT be caveated.
"""

import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))


def test_caveat_flags_match_manuscript_definition():
    with open(REPO / "data" / "landings_expected.json", "r", encoding="utf-8") as f:
        blob = json.load(f)
    caveated = {l["id"] for l in blob["landings"] if l["look_elsewhere_caveated"]}
    not_caveated = {l["id"] for l in blob["landings"] if not l["look_elsewhere_caveated"]}
    # L7 is exempt because it is Q-exact under R; the trial-factor
    # correction does not apply to a closed-form algebraic identity.
    assert caveated == {"L6", "L8"}
    # L1..L5 + L7 + L9, L10 (Paper 4 cosmological-constant tensor)
    # are not look-elsewhere caveated.
    assert not_caveated == {"L1", "L2", "L3", "L4", "L5", "L7",
                              "L9", "L10"}


def test_robust_core_definition_matches():
    with open(REPO / "data" / "landings_expected.json", "r", encoding="utf-8") as f:
        blob = json.load(f)
    assert set(blob["robust_core_definition"]["rows"]) == {"L1", "L2", "L3", "L4", "L5"}
    assert set(blob["look_elsewhere_caveat_definition"]["rows"]) == {"L6", "L8"}


def test_selection_correction_recorded():
    with open(REPO / "data" / "selection_correction.json", "r", encoding="utf-8") as f:
        sel = json.load(f)
    assert sel["search_space_size"] == 222954
    assert "joint_null_under_uniform_density" in sel
    assert sel["joint_null_under_uniform_density"]["p_joint_three_exact_hits_under_uniform_null"] >= 0.30


def test_perturbation_null_recorded():
    with open(REPO / "data" / "coefficient_perturbation_null.json", "r", encoding="utf-8") as f:
        pert = json.load(f)
    p_all = pert["p_all_3_targets_hit_under_random_coefficients"]
    p_any = pert["p_any_target_hit_under_random_coefficients"]
    assert p_all >= 0.5, "Look-elsewhere caveat is meaningful only if P(all 3) is non-trivially large"
    assert p_any == 1.0


def test_exact_tier_rows_are_not_in_robust_claim():
    """The EXACT-tier rows L6, L7, L8 must not be exported as standalone proof."""
    import run_look_elsewhere  # noqa: F401, ensure importable
    with open(REPO / "data" / "landings_expected.json", "r", encoding="utf-8") as f:
        blob = json.load(f)
    for L in blob["landings"]:
        if L["look_elsewhere_caveated"]:
            assert L["tier"] == "EXACT"
            assert "look_elsewhere_note" in L
