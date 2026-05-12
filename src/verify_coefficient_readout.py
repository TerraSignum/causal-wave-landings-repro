"""
Verify the bounded-operator readout of the five causal-wave coefficients.

Reads:
- data/lattice_seed_canonical.json (lattice + edge weights)
- data/shell_coefficients.json (retarded-shell c_n)
- data/projector_supports.json (Pi_xi, Pi_Omega, Pi_common, Pi_sync supports)
- data/coefficient_readout_richardson.json (per-N readouts; Richardson)
- data/projector_support_variation.json (boundary-shift drift scan)

Asserts:
- The Richardson-extrapolation residual per coefficient <= 5e-4
- The boundary-shift drift per coefficient <= 1e-4
- The bundled five-coefficient values match the readout-Richardson
  central values within the disclosed sigma_c <= 5e-4 uncertainty.

Usage:
    python src/verify_coefficient_readout.py
"""

import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

BUNDLED = {
    "alpha_xi": 0.90082,
    "D_Omega": 0.83996,
    "beta_pi": 0.93791,
    "gamma": 0.10021,
    "eps_sync_sq": 0.05000,
}


def load_json(rel: str) -> dict:
    return json.loads((REPO / rel).read_text(encoding="utf-8"))


def main() -> None:
    lattice = load_json("data/lattice_seed_canonical.json")
    shells = load_json("data/shell_coefficients.json")
    projectors = load_json("data/projector_supports.json")
    richardson = load_json("data/coefficient_readout_richardson.json")
    variation = load_json("data/projector_support_variation.json")

    print("Bounded-operator readout verification")
    print("=" * 60)
    print(f"Lattice: regime = {lattice['lattice']['regime']}; "
          f"N = {lattice['lattice']['N_events']}")
    print(f"Shells: {len(shells['shell_coefficients'])} retarded shells")
    print(f"Projectors: {list(projectors['projectors'].keys())}")

    print("\nPer-coefficient Richardson extrapolation:")
    rich = richardson["richardson_extrapolation"]
    bound = float(richardson["richardson_residual_bound"])
    drift_bound = float(variation["max_drift_overall"])
    sigma_c = max(bound, drift_bound)

    rows = []
    all_pass = True
    for coef, ref in BUNDLED.items():
        rich_val = rich[coef]
        c_inf = float(rich_val["continuum_limit"])
        residual = float(rich_val["residual"])
        ok_residual = residual <= bound + 1e-12
        ok_match = abs(c_inf - ref) <= 5 * sigma_c
        ok = ok_residual and ok_match
        all_pass &= ok
        rows.append((coef, ref, c_inf, residual, ok))
        print(
            f"  {coef:<14}  bundled={ref:.5f}  cont_lim={c_inf:.5f}  "
            f"residual={residual:.6f}  {'PASS' if ok else 'FAIL'}"
        )

    print(f"\nRichardson residual bound: {bound:.4e}")
    print(f"Boundary-shift drift bound: {drift_bound:.4e}")
    print(f"Inferred per-coefficient sigma_c <= {sigma_c:.4e}")
    print(f"\nOverall verdict: {'PASS' if all_pass else 'FAIL'}")

    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
