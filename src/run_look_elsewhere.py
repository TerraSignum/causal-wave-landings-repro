r"""
Look-elsewhere control for the three EXACT-tier benchmark landings.

This script does NOT re-run the full enumeration of 222954 algebraic
compositions (that lives in the upstream pipeline
`outputs_causal_wave_universality/selection_correction.json`).
Instead, it loads the cached enumeration result, reports the
Bonferroni-corrected joint p-value across the three EXACT-tier targets,
and prints the look-elsewhere caveat in a form that is suitable for
inclusion in the manuscript and for downstream tests.

The structural claim of Paper 2 is the PRECISE-or-better closure of the
five robust rows L1-L5; the EXACT classification of L6, L7, L8 is
look-elsewhere-caveated and is reported as such, not as a standalone
proof.

Usage:
    python ./src/run_look_elsewhere.py
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUTPUTS = REPO / "outputs"
OUTPUTS.mkdir(parents=True, exist_ok=True)


def main():
    with open(DATA / "selection_correction.json", "r", encoding="utf-8") as f:
        sel = json.load(f)
    with open(DATA / "coefficient_perturbation_null.json", "r", encoding="utf-8") as f:
        pert = json.load(f)

    print("=" * 78)
    print("Look-elsewhere control for EXACT-tier landings (L6, L7, L8)")
    print("=" * 78)
    print()

    print("--- Enumeration null ---")
    print(f"  Search-space size:     {sel['search_space_size']}")
    print(f"  In-domain compositions: {sel['in_domain_size']}")
    print(f"  Targets considered:    {sel['n_targets']}")
    print()

    print("--- Per-target hit rates at residual cut 1e-4 ---")
    for tname, t in sel["targets"].items():
        hits = t.get("exact_hits_at_1e-4_cut", "n/a")
        nullhits = t.get("expected_uniform_null", "n/a")
        ratio = t.get("observed_over_expected", "n/a")
        print(f"  {tname:<28} hits={hits}  null={nullhits}  ratio={ratio}")
    print()

    print("--- Joint null (uniform combinatorial) ---")
    p_joint = sel["joint_null_under_uniform_density"]["p_joint_three_exact_hits_under_uniform_null"]
    print(f"  p_joint(3 EXACT hits | uniform null) = {p_joint:.3f}")
    print(f"  Interpretation:")
    print(f"    {sel['joint_null_under_uniform_density']['interpretation']}")
    print()

    print("--- Coefficient-perturbation null ---")
    print(f"  N random trials: {pert['n_random_trials']}")
    print(f"  Random range:    {pert['random_coefficient_range']}")
    print(f"  P(all 3 EXACT | random coefficients): "
          f"{pert['p_all_3_targets_hit_under_random_coefficients']:.3f}")
    print(f"  P(any EXACT  | random coefficients): "
          f"{pert['p_any_target_hit_under_random_coefficients']:.3f}")
    print()
    for tname, st in pert["null_statistics"].items():
        p_better = st["p_random_residual_le_baseline"]
        print(f"  {tname:<28} P(random residual <= baseline) = {p_better:.3f}")
    print()

    print("--- Look-elsewhere verdict ---")
    print("  L6 (sin^2 theta_W):              EXACT but look-elsewhere-caveated.")
    print("  L7 (BH entropy 1/4):             EXACT and Q-exact under R (LE-exempt:")
    print("    alpha_xi/2 - 2 gamma = 9/20 - 1/5 = 1/4 in Q with N_gen=3 unique).")
    print("  L8 (Einstein gap 2/3):           EXACT but look-elsewhere-caveated;")
    print("    the 2/3 exponent is independently supported by an analytic")
    print("    Delta_E(N) <= C_0 N^(-2/3) bound in the gravitational-closure paper.")
    print("  L1-L5 (robust PRECISE_2.5 core): NOT look-elsewhere-caveated;")
    print("    these are the robust closure claim of the paper.")
    print()

    out = {
        "look_elsewhere_caveated_rows": ["L6", "L8"],
        "Q_exact_under_R_rows": ["L7"],
        "robust_core_rows": ["L1", "L2", "L3", "L4", "L5"],
        "p_joint_uniform_null": p_joint,
        "p_all_3_random_coefficients": pert["p_all_3_targets_hit_under_random_coefficients"],
        "p_any_random_coefficients": pert["p_any_target_hit_under_random_coefficients"],
        "structural_claim": "PRECISE_2.5-or-better closure on rows L1-L5 with no fits to targets",
        "non_claim": "Standalone EXACT-tier proof of L6, L7, L8 as structural identities",
    }
    out_path = OUTPUTS / "look_elsewhere_report.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
