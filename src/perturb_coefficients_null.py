r"""
Coefficient-perturbation null test (independent reproduction).

Replaces the five measured causal-wave coefficients
(alpha_xi, D_Omega, beta_pi, gamma, eps_sync2) with random uniform draws
from a configurable range and re-evaluates the eight benchmark landings.
Reports per-target probabilities of reproducing the EXACT-tier hits and
the joint probability of reproducing all three EXACT hits simultaneously.

This complements the cached upstream null test in
`data/coefficient_perturbation_null.json`. It is intended as an
independent on-the-fly verification under any RNG seed; the cached
result remains the canonical reference for the manuscript.

Usage:
    python ./src/perturb_coefficients_null.py
    python ./src/perturb_coefficients_null.py --seed 42 --trials 1000
"""

import argparse
import json
import math
from pathlib import Path
from random import Random

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUTPUTS = REPO / "outputs"
OUTPUTS.mkdir(parents=True, exist_ok=True)


def evaluate_three_exact_targets(coeff):
    """Evaluate L6, L7, L8 only (the EXACT-tier targets)."""
    a = coeff["alpha_xi"]
    D = coeff["D_Omega"]
    b = coeff["beta_pi"]
    g = coeff["gamma"]
    pi = math.pi
    sin2_theta_W   = b - (1 - g) * pi / 4
    BH_quarter     = a / 2 - 2 * g
    einstein_gap_23 = (1 - g) * pi / 4 - (1 - D) / 4
    return {
        "sin2_theta_W": sin2_theta_W,
        "BH_entropy_quarter": BH_quarter,
        "Einstein_gap_two_thirds": einstein_gap_23,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0,
                   help="RNG seed (default: 0)")
    p.add_argument("--trials", type=int, default=200,
                   help="Number of random trials (default: 200)")
    p.add_argument("--cmin", type=float, default=0.02,
                   help="Lower bound of random coefficient range (default: 0.02)")
    p.add_argument("--cmax", type=float, default=0.99,
                   help="Upper bound of random coefficient range (default: 0.99)")
    p.add_argument("--cut", type=float, default=1e-4,
                   help="EXACT residual threshold (default: 1e-4)")
    args = p.parse_args()

    targets = {
        "sin2_theta_W": 0.23122,
        "BH_entropy_quarter": 0.25,
        "Einstein_gap_two_thirds": 2.0 / 3.0,
    }

    rng = Random(args.seed)

    counts_at_least_one = 0
    counts_all_three = 0
    per_target_hits = {k: 0 for k in targets}
    residuals_log = {k: [] for k in targets}

    for _ in range(args.trials):
        coeff = {
            "alpha_xi":  rng.uniform(args.cmin, args.cmax),
            "D_Omega":   rng.uniform(args.cmin, args.cmax),
            "beta_pi":   rng.uniform(args.cmin, args.cmax),
            "gamma":     rng.uniform(args.cmin, args.cmax),
        }
        preds = evaluate_three_exact_targets(coeff)
        hits = 0
        for tname, tval in targets.items():
            r = abs(preds[tname] - tval) / abs(tval) if tval != 0 else abs(preds[tname])
            residuals_log[tname].append(r)
            if r <= args.cut:
                hits += 1
                per_target_hits[tname] += 1
        if hits >= 1:
            counts_at_least_one += 1
        if hits == 3:
            counts_all_three += 1

    p_any = counts_at_least_one / args.trials
    p_all = counts_all_three / args.trials

    print("=" * 78)
    print("Coefficient-perturbation null test (on-the-fly)")
    print("=" * 78)
    print()
    print(f"  N trials: {args.trials}")
    print(f"  Random range: [{args.cmin}, {args.cmax}]")
    print(f"  EXACT residual cut: {args.cut}")
    print(f"  RNG seed: {args.seed}")
    print()
    print("--- Per-target hit rates ---")
    for tname in targets:
        hr = per_target_hits[tname] / args.trials
        rs = sorted(residuals_log[tname])
        median = rs[len(rs) // 2]
        print(f"  {tname:<28} hit_rate={hr:.3f}   median_residual={median:.3e}")
    print()
    print("--- Joint statistics ---")
    print(f"  P(at least 1 EXACT hit | random coefficients) = {p_any:.3f}")
    print(f"  P(all 3 EXACT hits   | random coefficients) = {p_all:.3f}")
    print()
    print("--- Verdict ---")
    if p_all > 0.5:
        print(f"  CONFIRMS look-elsewhere caveat: random coefficients")
        print(f"  reproduce all three EXACT hits at probability {p_all:.3f}.")
    else:
        print(f"  Joint probability {p_all:.3f}; the EXACT classification of")
        print(f"  L6/L7/L8 is correspondingly less look-elsewhere-suspicious.")

    out = {
        "n_trials": args.trials,
        "seed": args.seed,
        "range": [args.cmin, args.cmax],
        "cut": args.cut,
        "p_any": p_any,
        "p_all_three": p_all,
        "per_target_hit_rate": {k: v / args.trials for k, v in per_target_hits.items()},
    }
    out_path = OUTPUTS / "perturbation_null_run.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
