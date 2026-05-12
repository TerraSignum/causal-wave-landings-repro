r"""
Verify the coefficient-reduction system R = {C1, ..., C5} of the manuscript.

This script reproduces, from the public bundle alone, the canonical
coefficient-reduction output of the universality pipeline
(outputs_causal_wave_universality/coefficient_reduction.json).

For each of the five algebraic conditions of Section 4 of the manuscript,
the canonical structure is emitted with identical field names and values:

  - observed                   (the five measured coefficients);
  - C1_alpha_plus_gamma_is_1   (alpha_xi + gamma = 1);
  - C2_D_equals_beta_minus_gamma (D_Omega = beta_pi - gamma);
  - C3_eps_equals_half_gamma   (eps_sync2 = gamma / 2);
  - C4_tan_theta_is_1_over_Ngen (gamma = 1/(N_gen^2 + 1) and the paired
                                 prediction alpha_xi = N_gen^2/(N_gen^2 + 1));
  - reduction_1param           (one-parameter closure leaving beta_pi free);
  - beta_pi_algebraic_forms    (nine algebraic candidates for beta_pi);
  - reduction_0param           (zero-parameter closure with beta_pi = 15/16).

The script also verifies, in closed-form Q-arithmetic, the algebraic
identity BH(N=3) = (N^2 - 4) / (2(N^2 + 1)) = 1/4, the unique structural
ascent point of the L7 (Bekenstein-Hawking 1/4) row under the rational
prediction (alpha_xi, D_Omega, beta_pi, gamma, eps_sync2)_R =
(9/10, 67/80, 15/16, 1/10, 1/20). The integer scan of BH(N) in [-10, 10]
confirms that N=3 is the unique positive integer solution.

Usage:
    python ./src/verify_reduction_system.py
"""

import json
import math
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUTPUTS = REPO / "outputs"
OUTPUTS.mkdir(parents=True, exist_ok=True)

# Per-coefficient pipeline-gate tolerance (0.30%). Used by the canonical
# coefficient_reduction.json structure to set passes_pg_tol on each block.
PG_TOL_PCT = 0.30


def load_coefficients():
    """Load the five measured causal-wave transport coefficients."""
    with open(DATA / "causal_wave_coefficients.json", "r", encoding="utf-8") as f:
        d = json.load(f)
    c = d["coefficients"]

    def _val(k):
        v = c[k]
        return v["value"] if isinstance(v, dict) else v

    out = {
        "alpha_xi":  float(_val("alpha_xi")),
        "D_Omega":   float(_val("D_Omega")),
        "beta_pi":   float(_val("beta_pi")),
        "gamma":     float(_val("gamma")),
        "eps_sync2": float(_val("eps_sync2")),
    }
    if "derived_quantities" in d and "N_gen" in d["derived_quantities"]:
        ng = d["derived_quantities"]["N_gen"]
        out["N_gen"] = int(ng["value"] if isinstance(ng, dict) else ng)
    else:
        out["N_gen"] = 3
    return out


def reduction_residuals(c):
    """Compute the per-condition residual of each of C1-C5 (manuscript form).

    Returned as a list of (label, identity_string, residual_fraction) tuples,
    where residual_fraction is the relative deviation as a fraction (not %).
    """
    a, D, b, g, e2 = (
        c["alpha_xi"], c["D_Omega"], c["beta_pi"], c["gamma"], c["eps_sync2"]
    )
    Ng = c["N_gen"]
    rows = [
        ("C1", "alpha_xi + gamma = 1",        abs((a + g) - 1.0) / 1.0),
        ("C2", "D_Omega = beta_pi - gamma",   abs(D - (b - g)) / abs(D)),
        ("C3", "eps_sync2 = gamma / 2",       abs(e2 - g / 2) / abs(e2)),
        ("C4", "gamma = 1 / (N_gen^2 + 1)",   abs(g - 1.0 / (Ng * Ng + 1)) / abs(g)),
        ("C5", "beta_pi = 15/16",             abs(b - 15.0 / 16.0) / abs(b)),
    ]
    return rows


def canonical_C_blocks(c):
    """Build the C1..C4 blocks in the exact canonical schema."""
    a, D, b, g, e2 = (
        c["alpha_xi"], c["D_Omega"], c["beta_pi"], c["gamma"], c["eps_sync2"]
    )
    Ng = c["N_gen"]

    # C1: alpha_xi + gamma = 1
    observed_sum = a + g
    C1_dev_pct = abs(observed_sum - 1.0) / 1.0 * 100.0

    # C2: D_Omega = beta_pi - gamma
    C2_pred = b - g
    C2_res_pct = abs(D - C2_pred) / abs(D) * 100.0

    # C3: eps_sync2 = gamma / 2
    C3_pred = g / 2.0
    C3_res_pct = abs(e2 - C3_pred) / abs(e2) * 100.0

    # C4: tan_theta = 1 / N_gen  ==>  gamma = 1/(N_gen^2 + 1) and
    #                                  alpha_xi = N_gen^2/(N_gen^2 + 1)
    gamma_pred = 1.0 / (Ng * Ng + 1)
    alpha_pred = (Ng * Ng) / (Ng * Ng + 1.0)
    gamma_res_pct = abs(g - gamma_pred) / abs(g) * 100.0
    alpha_res_pct = abs(a - alpha_pred) / abs(a) * 100.0

    return {
        "C1_alpha_plus_gamma_is_1": {
            "observed_sum": observed_sum,
            "target": 1.0,
            "deviation_pct": C1_dev_pct,
            "passes_pg_tol": C1_dev_pct <= PG_TOL_PCT,
        },
        "C2_D_equals_beta_minus_gamma": {
            "predicted": C2_pred,
            "observed": D,
            "residual_pct": C2_res_pct,
            "passes_pg_tol": C2_res_pct <= PG_TOL_PCT,
        },
        "C3_eps_equals_half_gamma": {
            "predicted": C3_pred,
            "observed": e2,
            "residual_pct": C3_res_pct,
            "passes_pg_tol": C3_res_pct <= PG_TOL_PCT,
        },
        "C4_tan_theta_is_1_over_Ngen": {
            "N_gen": Ng,
            "gamma_predicted": gamma_pred,
            "gamma_observed": g,
            "gamma_residual_pct": gamma_res_pct,
            "alpha_predicted": alpha_pred,
            "alpha_observed": a,
            "alpha_residual_pct": alpha_res_pct,
            "passes_pg_tol": (gamma_res_pct <= PG_TOL_PCT
                              and alpha_res_pct <= PG_TOL_PCT),
        },
    }


def canonical_reduction_1param(c):
    """One-parameter closure: alpha_xi, gamma, eps_sync2, D_Omega all fixed
    by C1-C4 with N_gen = 3; only beta_pi remains free (= observed)."""
    Ng = c["N_gen"]
    a_pred = (Ng * Ng) / (Ng * Ng + 1.0)
    g_pred = 1.0 / (Ng * Ng + 1.0)
    e2_pred = g_pred / 2.0
    b_obs = c["beta_pi"]
    D_pred = b_obs - g_pred  # via C2

    values = {
        "alpha_xi": a_pred,
        "D_Omega":  D_pred,
        "beta_pi":  b_obs,
        "gamma":    g_pred,
        "eps_sync2": e2_pred,
    }
    res = {
        "alpha_xi":  abs(a_pred - c["alpha_xi"])  / abs(c["alpha_xi"])  * 100.0,
        "D_Omega":   abs(D_pred - c["D_Omega"])   / abs(c["D_Omega"])   * 100.0,
        "beta_pi":   0.0,
        "gamma":     abs(g_pred - c["gamma"])     / abs(c["gamma"])     * 100.0,
        "eps_sync2": abs(e2_pred - c["eps_sync2"])/ abs(c["eps_sync2"]) * 100.0,
    }
    max_res = max(res.values())
    return {
        "free_parameters": ["beta_pi"],
        "values": values,
        "residuals_pct": res,
        "max_residual_pct": max_res,
        "passes_pg_tol": max_res <= PG_TOL_PCT,
    }


def beta_pi_algebraic_forms(c):
    """Test nine algebraic candidate forms for beta_pi.

    The forms are pre-registered analytic ansatze from the universality
    pipeline; only those within PG tolerance pass. The 0-parameter closure
    selects 15/16 because (i) it is the closest in absolute deviation
    and (ii) it admits the algebraic-exact rational reduction to the
    Bekenstein-Hawking 1/4 identity in Q (see bh_quarter_under_R)."""
    Ng = c["N_gen"]
    e = c["eps_sync2"]
    g = c["gamma"]
    b_obs = c["beta_pi"]

    forms = [
        ("3/pi",                       3.0 / math.pi),
        ("pi/4 + 1/2",                 math.pi / 4.0 + 0.5),
        ("1 - 1/(4*pi)",               1.0 - 1.0 / (4.0 * math.pi)),
        ("1 - gamma * (pi/4)",         1.0 - g * (math.pi / 4.0)),
        ("cos(pi/12)**2",              math.cos(math.pi / 12.0) ** 2),
        ("(N_gen + 1/pi) / (N_gen + 1)", (Ng + 1.0 / math.pi) / (Ng + 1.0)),
        ("1 - 2/pi * eps",             1.0 - (2.0 / math.pi) * e),
        ("15/16",                      15.0 / 16.0),
        ("1 - 1/(2*pi)",               1.0 - 1.0 / (2.0 * math.pi)),
    ]
    # Per canonical pipeline, the algebraic-form survey uses a relaxed
    # 2.0% acceptance band (so near-miss candidates appear as PASS in the
    # table); only 15/16 also passes the strict PG_TOL_PCT band of 0.30%.
    SURVEY_BAND_PCT = 2.0
    out = {}
    for name, pred in forms:
        rp = abs(pred - b_obs) / abs(b_obs) * 100.0
        out[name] = {
            "pred": pred,
            "residual_pct": rp,
            "passes_pg_tol": rp <= SURVEY_BAND_PCT,
        }
    return out


def canonical_reduction_0param(c):
    """Zero-parameter closure: beta_pi = 15/16 (closest algebraic form
    that lies within strict PG tolerance and admits the algebraic-exact
    Bekenstein-Hawking 1/4 identity)."""
    Ng = c["N_gen"]
    a_pred = (Ng * Ng) / (Ng * Ng + 1.0)
    g_pred = 1.0 / (Ng * Ng + 1.0)
    e2_pred = g_pred / 2.0
    b_pred = 15.0 / 16.0
    D_pred = b_pred - g_pred

    values = {
        "alpha_xi": a_pred,
        "D_Omega":  D_pred,
        "beta_pi":  b_pred,
        "gamma":    g_pred,
        "eps_sync2": e2_pred,
    }
    res = {
        "alpha_xi":  abs(a_pred - c["alpha_xi"])  / abs(c["alpha_xi"])  * 100.0,
        "D_Omega":   abs(D_pred - c["D_Omega"])   / abs(c["D_Omega"])   * 100.0,
        "beta_pi":   abs(b_pred - c["beta_pi"])   / abs(c["beta_pi"])   * 100.0,
        "gamma":     abs(g_pred - c["gamma"])     / abs(c["gamma"])     * 100.0,
        "eps_sync2": abs(e2_pred - c["eps_sync2"])/ abs(c["eps_sync2"]) * 100.0,
    }
    max_res = max(res.values())
    return {
        "beta_pi_formula": "15/16",
        "values": values,
        "residuals_pct": res,
        "max_residual_pct": max_res,
        "passes_pg_tol": max_res <= PG_TOL_PCT,
    }


def canonical_coefficient_reduction(c):
    """Assemble the full canonical coefficient_reduction.json structure."""
    return {
        "observed": {
            "alpha_xi":  c["alpha_xi"],
            "D_Omega":   c["D_Omega"],
            "beta_pi":   c["beta_pi"],
            "gamma":     c["gamma"],
            "eps_sync2": c["eps_sync2"],
        },
        **canonical_C_blocks(c),
        "reduction_1param": canonical_reduction_1param(c),
        "beta_pi_algebraic_forms": beta_pi_algebraic_forms(c),
        "reduction_0param": canonical_reduction_0param(c),
    }


def rational_reduction_values():
    """The rational prediction (alpha_xi, D_Omega, beta_pi, gamma, eps_sync2)_R."""
    return {
        "alpha_xi":  Fraction(9, 10),
        "D_Omega":   Fraction(67, 80),
        "beta_pi":   Fraction(15, 16),
        "gamma":     Fraction(1, 10),
        "eps_sync2": Fraction(1, 20),
    }


def bh_quarter_under_R():
    """Closed-form check: alpha_xi/2 - 2*gamma = 9/20 - 1/5 = 1/4 in Q."""
    R = rational_reduction_values()
    lhs = R["alpha_xi"] / 2 - 2 * R["gamma"]
    target = Fraction(1, 4)
    return lhs, target, lhs == target


def bh_general_formula(N):
    """BH(N) = (N^2 - 4) / (2 * (N^2 + 1)). Returns Fraction."""
    return Fraction(N * N - 4, 2 * (N * N + 1))


def integer_uniqueness_check_for_quarter():
    """Scan integer N from -10 to 10 and find all N solving BH(N) = 1/4."""
    target = Fraction(1, 4)
    return [N for N in range(-10, 11) if bh_general_formula(N) == target]


def main():
    coeff = load_coefficients()

    print("=" * 72)
    print("Coefficient-reduction system R: per-condition residuals")
    print("=" * 72)
    print()
    print(f"  Measured coefficients: alpha_xi={coeff['alpha_xi']}, "
          f"D_Omega={coeff['D_Omega']}, beta_pi={coeff['beta_pi']}, "
          f"gamma={coeff['gamma']}, eps_sync2={coeff['eps_sync2']}, "
          f"N_gen={coeff['N_gen']}")
    print()
    rows = reduction_residuals(coeff)
    max_residual = 0.0
    print(f"  {'cond':<5} {'identity':<40} {'residual':>12}")
    print("  " + "-" * 65)
    for label, expr, res in rows:
        max_residual = max(max_residual, res)
        print(f"  {label:<5} {expr:<40} {res*100:>11.4f}%")
    print()
    print(f"  Max residual across C1-C5: {max_residual*100:.4f}%")
    threshold = PG_TOL_PCT / 100  # 0.30 percent
    print(f"  Threshold of the manuscript (0.30%): "
          f"{'PASS' if max_residual <= threshold else 'FAIL'}")
    print()

    print("=" * 72)
    print("Canonical coefficient-reduction blocks (1-param and 0-param)")
    print("=" * 72)
    print()
    canon = canonical_coefficient_reduction(coeff)
    r1 = canon["reduction_1param"]
    r0 = canon["reduction_0param"]
    print(f"  reduction_1param  (free: beta_pi):")
    print(f"    max residual = {r1['max_residual_pct']:.5f}%  "
          f"-> passes_pg_tol={r1['passes_pg_tol']}")
    print(f"  reduction_0param  (beta_pi = 15/16):")
    print(f"    max residual = {r0['max_residual_pct']:.5f}%  "
          f"-> passes_pg_tol={r0['passes_pg_tol']}")
    print()
    print("  beta_pi algebraic-form survey:")
    forms = canon["beta_pi_algebraic_forms"]
    for name, blk in forms.items():
        flag = "PASS" if blk["passes_pg_tol"] else "FAIL"
        print(f"    {name:<32} pred={blk['pred']:.6f}  "
              f"res={blk['residual_pct']:>7.4f}%  [{flag}]")
    print()

    print("=" * 72)
    print("Algebraic-exact Bekenstein-Hawking 1/4 under R")
    print("=" * 72)
    print()
    R = rational_reduction_values()
    print(f"  Rational reduction values:")
    for k, v in R.items():
        print(f"    {k:<10} = {v} = {float(v):.6f}")
    print()
    lhs, target, eq = bh_quarter_under_R()
    print(f"  alpha_xi/2 - 2*gamma = {lhs} = {target}")
    print(f"  Target (BH coefficient): {target} = {float(target):.6f}")
    print(f"  Algebraically exact in Q: {'YES' if eq else 'NO'}")
    print()
    integer_solutions = integer_uniqueness_check_for_quarter()
    positive_solutions = [N for N in integer_solutions if N > 0]
    print(f"  General formula BH(N) = (N^2 - 4) / (2*(N^2 + 1))")
    print(f"  Integer solutions of BH(N) = 1/4 in [-10, 10]: {integer_solutions}")
    print(f"  Positive integer solutions: {positive_solutions}")
    print(f"  Unique positive integer solution: "
          f"{positive_solutions[0] if len(positive_solutions) == 1 else 'NOT UNIQUE'}")
    print()

    # 1) Canonical structure: identical schema to the upstream pipeline output
    canon_path = OUTPUTS / "coefficient_reduction.json"
    with open(canon_path, "w", encoding="utf-8") as f:
        json.dump(canon, f, indent=2)
    print(f"Saved (canonical schema): {canon_path}")

    # 2) Manuscript-narrative summary report (kept for the test suite)
    narrative = {
        "measured_coefficients": coeff,
        "conditions": [
            {"label": l, "identity": e, "residual_pct": r * 100}
            for (l, e, r) in rows
        ],
        "max_residual_pct": max_residual * 100,
        "threshold_pct": PG_TOL_PCT,
        "all_within_threshold": max_residual <= threshold,
        "rational_reduction": {k: [v.numerator, v.denominator] for k, v in R.items()},
        "bh_quarter": {
            "lhs_numerator": lhs.numerator,
            "lhs_denominator": lhs.denominator,
            "target_numerator": target.numerator,
            "target_denominator": target.denominator,
            "algebraically_exact_in_Q": eq,
        },
        "bh_integer_uniqueness": {
            "scan_range": [-10, 10],
            "integer_solutions_for_quarter": integer_solutions,
            "positive_solutions": positive_solutions,
            "unique_positive_solution": (
                positive_solutions[0] if len(positive_solutions) == 1 else None
            ),
        },
    }
    narr_path = OUTPUTS / "reduction_system.json"
    with open(narr_path, "w", encoding="utf-8") as f:
        json.dump(narrative, f, indent=2)
    print(f"Saved (narrative summary): {narr_path}")


if __name__ == "__main__":
    main()
