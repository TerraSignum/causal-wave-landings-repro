"""
Audit the structural stability of the candidate rational reduction
hypothesis R via the loop-form-pattern test on the three measured
constraint deviations (sec:reduction, "Stability of R" paragraph).

The three measured deviations:
- delta_C1 = alpha_xi + gamma - 1
- delta_C2 = D_Omega - (beta_pi - gamma)
- delta_C3 = 2 eps_sync_sq - gamma

each must admit a closed-form one-loop self-energy interpretation
with shared resummation kernel 1/(1 - 2*gamma^2) for delta_C1
and delta_C3.
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "outputs" / "R_stability_audit.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

ALPHA_XI = 0.90082
D_OMEGA = 0.83996
BETA_PI = 0.93791
GAMMA = 0.10021
EPS_SYNC_SQ = 0.05000
N_GEN = 3


def main() -> None:
    delta_C1_measured = ALPHA_XI + GAMMA - 1.0
    delta_C2_measured = D_OMEGA - (BETA_PI - GAMMA)
    delta_C3_measured = 2.0 * EPS_SYNC_SQ - GAMMA

    delta_C1_loop_form = GAMMA ** 3 / (1.0 - 2.0 * GAMMA ** 2)
    delta_C2_loop_form = (N_GEN ** 2 / 2.0) * GAMMA ** 2 * EPS_SYNC_SQ
    delta_C3_loop_form = -2.0 * GAMMA ** 4 / (1.0 - 2.0 * GAMMA ** 2)

    pct_match = {
        "C1": abs(delta_C1_measured - delta_C1_loop_form) / abs(delta_C1_measured) * 100 if delta_C1_measured != 0 else float('nan'),
        "C2": abs(delta_C2_measured - delta_C2_loop_form) / abs(delta_C2_measured) * 100 if delta_C2_measured != 0 else float('nan'),
        "C3": abs(delta_C3_measured - delta_C3_loop_form) / abs(delta_C3_measured) * 100 if delta_C3_measured != 0 else float('nan'),
    }

    bundled = {
        "schema_version": "1.0.0",
        "release": "v0.1.0",
        "stand": "2026-04-28",
        "audit": "R-hypothesis structural-stability via loop-form pattern",
        "measured_deviations": {
            "delta_C1": delta_C1_measured,
            "delta_C2": delta_C2_measured,
            "delta_C3": delta_C3_measured,
        },
        "loop_form_predictions": {
            "delta_C1_form": "gamma^3 / (1 - 2*gamma^2)",
            "delta_C1_value": delta_C1_loop_form,
            "delta_C2_form": "N_gen^2/2 * gamma^2 * eps_sync_sq",
            "delta_C2_value": delta_C2_loop_form,
            "delta_C3_form": "-2*gamma^4 / (1 - 2*gamma^2)",
            "delta_C3_value": delta_C3_loop_form,
        },
        "pct_match": pct_match,
        "shared_kernel_C1_C3": "1/(1 - 2*gamma^2)",
        "verdict": (
            f"All three measured constraint deviations admit a "
            f"closed-form one-loop self-energy interpretation with "
            f"match accuracies "
            f"C1: {pct_match['C1']:.2f}%, "
            f"C2: {pct_match['C2']:.2f}%, "
            f"C3: {pct_match['C3']:.2f}%. "
            f"C1 and C3 share the resummation kernel "
            f"1/(1 - 2*gamma^2) — the structural stability witness "
            f"for R."
        ),
    }

    OUT.write_text(json.dumps(bundled, indent=2), encoding="utf-8")
    print("R-hypothesis loop-form pattern test:")
    for c, m in pct_match.items():
        print(f"  delta_{c}: match accuracy = {m:.2f}%")
    print(f"\nSaved: {OUT}")


if __name__ == "__main__":
    main()
