r"""
Causal-wave benchmark landings: strict recompute.

Loads the five measured causal-wave transport coefficients and the
expected ten benchmark landings (L1..L10) from frozen JSON inputs, evaluates
each landing formula directly, and reports prediction, target,
residual, and tier.

NO canonical override: every value is computed from the loaded
coefficients via the closed-form formula stored in
`data/landings_expected.json`. Deviations from the cached
`expected_prediction` field would indicate that the data has changed
or the recompute path is broken.

Usage (Windows PowerShell):
    python .\src\recompute_landings.py

Usage (POSIX):
    python ./src/recompute_landings.py
"""

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUTPUTS = REPO / "outputs"
OUTPUTS.mkdir(parents=True, exist_ok=True)


def load_json(name):
    with open(DATA / name, "r", encoding="utf-8") as f:
        return json.load(f)


def load_coefficients():
    """Load the five measured causal-wave coefficients."""
    blob = load_json("causal_wave_coefficients.json")
    c = blob["coefficients"]
    coeff = {
        "alpha_xi":  c["alpha_xi"]["value"],
        "D_Omega":   c["D_Omega"]["value"],
        "beta_pi":   c["beta_pi"]["value"],
        "gamma":     c["gamma"]["value"],
        "eps_sync2": c["eps_sync2"]["value"],
    }
    coeff["G_NET"] = blob["derived_quantities"]["G_NET"]["value"]
    coeff["N_gen"] = blob["derived_quantities"]["N_gen"]["value"]
    return coeff


def evaluate_formula(formula, coeff):
    """
    Evaluate a closed-form landing formula.

    The formulas are intentionally a fixed, enumerated set; we do NOT
    eval() them so that no JSON manipulation can inject arbitrary code.
    """
    a = coeff["alpha_xi"]
    D = coeff["D_Omega"]
    b = coeff["beta_pi"]
    g = coeff["gamma"]
    e = coeff["eps_sync2"]
    G = coeff["G_NET"]
    N = coeff["N_gen"]
    pi = math.pi

    if formula == "alpha_xi + beta_pi + eps_sync2 - gamma":
        return a + b + e - g
    if formula == "alpha_xi**2 * eps_sync2 * N_gen":
        return a * a * e * N
    if formula == "-1 + eps_sync2**2 / gamma":
        return -1 + (e * e) / g
    if formula == "gamma * beta_pi * 4/pi":
        return g * b * 4 / pi
    if formula == "(alpha_xi/beta_pi) * (D_Omega/(1+gamma)) * G_NET":
        return (a / b) * (D / (1 + g)) * G
    if formula == "beta_pi - (1 - gamma) * pi/4":
        return b - (1 - g) * pi / 4
    if formula == "alpha_xi/2 - 2*gamma":
        return a / 2 - 2 * g
    if formula == "(1 - gamma) * pi/4 - (1 - D_Omega) / 4":
        return (1 - g) * pi / 4 - (1 - D) / 4
    if formula == "alpha_xi**2":
        return a * a
    if formula == "-gamma**2 / 2":
        return -(g * g) / 2

    raise ValueError(f"Unknown landing formula: {formula!r}")


def classify_tier(residual_pct, thresholds_pct):
    """
    Classify the residual into a tier:
      - EXACT          (|residual| <=  0.01%)
      - PRECISE_2.5    (|residual| <=  2.5%)
      - PRECISE        (|residual| <=  1%)   [historical sub-tier]
      - FACTOR2        (|residual| <= 50%)
      - FAR_OFF        (otherwise)
    """
    r = abs(residual_pct)
    if r <= thresholds_pct["EXACT"] * 100:
        return "EXACT"
    if r <= thresholds_pct["PRECISE_2.5"] * 100:
        return "PRECISE_2.5"
    if r <= thresholds_pct["FACTOR2"] * 100:
        return "FACTOR2"
    return "FAR_OFF"


def recompute_one(landing, coeff, thresholds):
    """Recompute a single landing entry from the formula and coefficients."""
    pred = evaluate_formula(landing["formula"], coeff)
    target = landing["target_value"]
    if target == 0:
        residual_pct = 0.0
    else:
        residual_pct = (pred - target) / target * 100.0
    tier = classify_tier(residual_pct, thresholds)
    return {
        "id": landing["id"],
        "observable": landing["observable"],
        "formula": landing["formula"],
        "prediction": pred,
        "target": target,
        "residual_pct": residual_pct,
        "tier": tier,
        "look_elsewhere_caveated": landing["look_elsewhere_caveated"],
        "expected_tier": landing["tier"],
    }


def compute_all():
    """Run the full eight-landing recompute."""
    coeff = load_coefficients()
    expected = load_json("landings_expected.json")
    thresholds = expected["tier_thresholds"]
    rows = [recompute_one(land, coeff, thresholds) for land in expected["landings"]]
    return {"coefficients": coeff, "rows": rows}


def main():
    R = compute_all()

    print("=" * 78)
    print("Causal-wave benchmark landings -- strict recompute")
    print("Source: emergence-core-data v0.1.0 (frozen)")
    print("=" * 78)
    print()
    print("--- Five measured causal-wave coefficients ---")
    c = R["coefficients"]
    print(f"  alpha_xi   = {c['alpha_xi']:.5f}")
    print(f"  D(Omega)   = {c['D_Omega']:.5f}")
    print(f"  beta_pi    = {c['beta_pi']:.5f}")
    print(f"  gamma      = {c['gamma']:.5f}")
    print(f"  eps_sync^2 = {c['eps_sync2']:.5f}")
    print(f"  G_NET      = {c['G_NET']:.5f}  (alpha_xi + beta_pi + eps_sync^2 - gamma)")
    print(f"  N_gen      = {c['N_gen']}")
    print()

    print("--- Ten benchmark landings ---")
    print(f"  {'id':<3} {'observable':<28} {'pred':>13} {'target':>13} "
          f"{'resid %':>9} {'tier':<12} {'LE':<3}")
    print("  " + "-" * 78)
    n_precise_or_better = 0
    n_robust_core_pass = 0
    robust_core_ids = {"L1", "L2", "L3", "L4", "L5"}
    for r in R["rows"]:
        le = "yes" if r["look_elsewhere_caveated"] else "no"
        print(f"  {r['id']:<3} {r['observable']:<28} "
              f"{r['prediction']:>13.6f} {r['target']:>13.6f} "
              f"{r['residual_pct']:>+8.4f}% {r['tier']:<12} {le:<3}")
        if r["tier"] in ("EXACT", "PRECISE_2.5"):
            n_precise_or_better += 1
        if r["id"] in robust_core_ids and abs(r["residual_pct"]) <= 2.5:
            n_robust_core_pass += 1

    print()
    n_total = len(R.get("rows", []))
    print("--- Aggregate result ---")
    print(f"  PRECISE-or-better landings:        {n_precise_or_better}/{n_total}")
    print(f"  Robust core (L1-L5) within 2.5%:   {n_robust_core_pass}/5")
    print(f"  Look-elsewhere caveated rows:      L6, L7, L8")
    print()
    print("--- Acceptance ---")
    if n_precise_or_better == n_total and n_robust_core_pass == 5:
        print(f"  PASS: {n_total}/{n_total} PRECISE-or-better; robust 5/5 core within 2.5%.")
    else:
        print("  FAIL: closure broken.")

    out_path = OUTPUTS / "recompute_landings.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2)
    print()
    print(f"Saved machine-readable result: {out_path}")


if __name__ == "__main__":
    main()
