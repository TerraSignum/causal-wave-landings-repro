"""verify_carrier_operator_linear_stability.py

Linear-stability analysis of the carrier-transport operator
around the trivial vacuum. Two parts:

  (A) Algebraic, self-contained (no external data needed):
      - Rational reduction of the five carrier coefficients via R
      - Closed-form growth rate of common-phase mode (G_NET)
      - Closed-form critical eigenvalue lambda_crit
      Both quantities are returned as exact rationals.

  (B) Empirical cross-check (uses companion lattice data):
      - Builds the graph-Laplacian L = D - W from the lattice
        adjacency matrix and counts how many eigenvalues lie below
        lambda_crit. Predicted: exactly 1 (the constant mode at
        lambda = 0).
      - Reports the count for each lattice size in a within-canonical
        N-ladder. Empirical check uses companion lattice snapshot
        files which are not bundled in this repository to keep
        reproducibility lean; if the files are not visible the empirical
        check is skipped with a note.

Outputs to data/carrier_operator_linear_stability.json with full
algebraic and empirical results.
"""
from __future__ import annotations

import json
import os
import glob
from fractions import Fraction
from pathlib import Path

import numpy as np

# ============================================================
# PART A: Algebraic, self-contained
# ============================================================

# Five carrier coefficients (System R, Section "Coefficient-reduction")
ALPHA_XI  = Fraction(9, 10)     # Xi reaction rate
GAMMA     = Fraction(1, 10)     # damping rate
EPS_SYNC2 = Fraction(1, 20)     # sync persistence amplitude
BETA_PI   = Fraction(15, 16)    # common-phase / holonomy mode strength
D_OMEGA   = BETA_PI - GAMMA     # 67/80 (Einstein relation under R)

NUMERATOR    = ALPHA_XI - GAMMA + EPS_SYNC2     # 17/20
LAMBDA_CRIT  = NUMERATOR / D_OMEGA              # 68/67
G_NET        = ALPHA_XI + BETA_PI + EPS_SYNC2 - GAMMA  # 143/80


def algebraic_summary() -> dict:
    return {
        "alpha_xi": str(ALPHA_XI),
        "gamma": str(GAMMA),
        "eps_sync2": str(EPS_SYNC2),
        "beta_pi": str(BETA_PI),
        "D_Omega": str(D_OMEGA),
        "non_common_zero_numerator": str(NUMERATOR),
        "lambda_crit_rational": str(LAMBDA_CRIT),
        "lambda_crit_float": float(LAMBDA_CRIT),
        "G_NET_rational": str(G_NET),
        "G_NET_float": float(G_NET),
        "linear_growth_rate_non_common":
            "g(lambda) = -D(Omega) * lambda + (alpha_xi - gamma + eps_sync2) "
            "= (68 - 67*lambda) / 80",
        "linear_growth_rate_common_mode":
            "g_common = G_NET = alpha_xi + beta_pi + eps_sync2 - gamma = 143/80",
        "critical_eigenvalue_definition":
            "lambda_crit = (alpha_xi - gamma + eps_sync2) / D(Omega) = 68/67",
    }


# ============================================================
# PART B: Empirical cross-check on companion lattice ladder
# ============================================================

REPO_ROOT  = Path(__file__).resolve().parents[1]
WORKSPACE  = REPO_ROOT.parent
EMERGENCE_RESULTS = Path(os.environ.get("EMERGENCE_RESULTS_ROOT",
                                        str(WORKSPACE)))

# Within-canonical-regime ladder N -> [candidate snapshot dirs]
N_LADDER = {
    64:  ["results_d1_p5n64_24seeds",  "results_d1_p5n64"],
    72:  ["results_d1_p5n72_24seeds",  "results_d1_p5n72"],
    84:  ["results_d1_p5n84_24seeds",  "results_d1_p5n84"],
    100: ["results_d1_p5n100_24seeds", "results_d1_p5n100"],
    128: ["results_d1_p5n128_kq_fixed", "results_d1_p5n128"],
    200: ["results_d1_p5n200_8seeds",  "results_d1_p5n200"],
    256: ["results_d1_p5n256_12seeds", "results_d1_p5n256"],
    300: ["results_d1_p5n300_12seeds", "results_d1_p5n300"],
    512: ["results_d1_p5n512_12seeds", "results_d1_p5n512"],
}


def find_snapshot_npz(N: int) -> str | None:
    for sub in N_LADDER[N]:
        cand = sorted(glob.glob(str(EMERGENCE_RESULTS / sub / "*.snapshots.npz")))
        if cand:
            return cand[0]
    return None


def graph_laplacian_eigvals(W: np.ndarray) -> np.ndarray:
    Ws = (W + W.T) / 2.0
    deg = np.sum(Ws, axis=1)
    L = np.diag(deg) - Ws
    L = (L + L.T) / 2.0
    return np.linalg.eigvalsh(L)


def empirical_check() -> dict:
    lambda_crit_f = float(LAMBDA_CRIT)
    per_N: dict = {}
    for N in sorted(N_LADDER):
        path = find_snapshot_npz(N)
        if path is None:
            per_N[str(N)] = {"path": None, "status": "snapshot_not_visible"}
            continue
        try:
            d = np.load(path, allow_pickle=True)
            edge = d["edge_xi_snapshots"][:, -1, :, :]   # (S, N, N), final
            S = edge.shape[0]
            counts = []
            for s in range(S):
                ev = graph_laplacian_eigvals(edge[s])
                counts.append(int(np.sum(ev < lambda_crit_f)))
            counts_arr = np.array(counts)
            per_N[str(N)] = {
                "path": path,
                "n_seeds": int(S),
                "counts_below_crit": counts_arr.tolist(),
                "min_count": int(counts_arr.min()),
                "max_count": int(counts_arr.max()),
                "mean_count": float(counts_arr.mean()),
                "var_count": float(counts_arr.var()),
                "all_seeds_equal_one":
                    bool(np.all(counts_arr == 1)),
            }
        except Exception as exc:
            per_N[str(N)] = {"path": path, "status": "load_error",
                             "error": repr(exc)}

    # Symanzik 1/N on (mean count) / N — predicted: exactly 1/N
    Ns_visible = [int(k) for k, v in per_N.items()
                  if "mean_count" in v]
    if len(Ns_visible) >= 2:
        N_arr = np.array(sorted(Ns_visible), dtype=float)
        mean_arr = np.array([per_N[str(int(N))]["mean_count"]
                             for N in N_arr])
        frac = mean_arr / N_arr
        A = np.column_stack([np.ones_like(N_arr), 1.0 / N_arr])
        coef, *_ = np.linalg.lstsq(A, frac, rcond=None)
        residuals = frac - A @ coef
        rmse = float(np.sqrt(np.mean(residuals**2)))
        symanzik = {
            "Ns": [int(x) for x in N_arr],
            "growing_fraction": frac.tolist(),
            "fit_form": "f(N) = a + b/N",
            "a_inf": float(coef[0]),
            "b": float(coef[1]),
            "rmse": rmse,
            "expected_exact": "a=0, b=1",
        }
    else:
        symanzik = {"status": "insufficient_visible_lattice_sizes"}

    return {"per_N": per_N, "symanzik_one_over_N": symanzik}


# ============================================================
# Driver
# ============================================================

def main() -> int:
    out = {
        "audit": "carrier-operator-linear-stability",
        "stand": "2026-05-04",
        "algebraic": algebraic_summary(),
        "empirical": empirical_check(),
        "predicted_outcome":
            "for every lattice size in the within-canonical N-ladder "
            "the lattice graph-Laplacian has exactly one eigenvalue "
            "below lambda_crit = 68/67 (the constant mode), and the "
            "growing-mode fraction f(N) = (count below crit) / N "
            "satisfies f(N) = 1/N exactly to machine precision.",
    }
    out_dir = REPO_ROOT / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "carrier_operator_linear_stability.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    # Console summary
    print("=" * 70)
    print("Algebraic part:")
    print(f"  lambda_crit = {LAMBDA_CRIT} = {float(LAMBDA_CRIT):.10f}")
    print(f"  G_NET       = {G_NET} = {float(G_NET):.10f}")
    print()
    print("Empirical part:")
    emp = out["empirical"]
    sym = emp["symanzik_one_over_N"]
    if "a_inf" in sym:
        print(f"  Lattices visible: {sym['Ns']}")
        print(f"  Growing-mode fraction: {sym['growing_fraction']}")
        print(f"  Symanzik a + b/N: a_inf = {sym['a_inf']:+.6e}, "
              f"b = {sym['b']:+.6e}, rmse = {sym['rmse']:.2e}")
        print(f"  Expected exact: a=0, b=1")
    else:
        print(f"  Snapshot data not visible at "
              f"{EMERGENCE_RESULTS}; algebraic part remains valid.")
    print(f"\nWritten: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
