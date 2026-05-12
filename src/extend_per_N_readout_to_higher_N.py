"""Extend causal_wave_per_N_readout.json with the post-flip
regimes P5N256 (12-seed) and P5N512 (12-seed) that were not in
the original readout bundle.

The readout pipeline mirrors the per-seed coefficient readout of
peer_reviews/causal_wave_obs5_internal_external_calibration.py:
  alpha_raw = mean(|xi off-diagonal|) / (mean + 0.05)
  beta_raw  = 1 - 1/(1 + top1_eig * sqrt(N))
  diff_raw  = laplacian_diffusion_strength(xi)
followed by calibration with the bundled
(alpha_scale, beta_scale, diff_scale) factors, and the C1/C3
identities gamma = 1 - alpha, eps^2_sync = gamma/2, D_omega_C2
= beta - gamma.

The bundled JSON is written back in place (sorted by n_lat).

Usage:
  python src/extend_per_N_readout_to_higher_N.py
"""
from __future__ import annotations
import io
import json
import math
import os
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
EMERGENCE = os.path.dirname(REPO)

READOUT_PATH = os.path.join(REPO, "data", "causal_wave_per_N_readout.json")
# Two other repos bundle the same file; we mirror to all of them
MIRROR_PATHS = [
    os.path.join(EMERGENCE, "emergent-gr-anisotropic-source-dm-de-repro",
                 "data", "causal_wave_per_N_readout.json"),
    os.path.join(EMERGENCE, "emergent-gr-closure-repro",
                 "outputs", "causal_wave_per_N_readout.json"),
]

# Per-regime NPZ snapshot lookup (latest preferred run per regime)
NPZ_LOOKUP = {
    "P5N256": os.path.join(EMERGENCE, "results_d1_p5n256_12seeds",
                           "P5N256.snapshots.npz"),
    "P5N512": os.path.join(EMERGENCE, "results_d1_p5n512_12seeds",
                           "P5N512.snapshots.npz"),
}


def normalised_gram_eigvals(x: np.ndarray) -> np.ndarray:
    g = x @ x.T
    eig = np.linalg.eigvalsh(g)
    eig = eig[eig > 0]
    total = float(np.sum(eig))
    if total <= 1e-12:
        return np.zeros(x.shape[0])
    return np.sort(eig / total)[::-1]


def laplacian_diffusion_strength(xi: np.ndarray,
                                 xi_thresh: float = 0.5) -> float:
    x = np.array(np.where(np.isfinite(xi), xi, 0.0), dtype=float, copy=True)
    n = x.shape[0]
    np.fill_diagonal(x, 0.0)
    adj = (x > xi_thresh).astype(float)
    weight = x * adj
    deg = weight.sum(axis=1) + 1e-12
    deg_inv_sqrt = 1.0 / np.sqrt(deg)
    L = np.eye(n) - deg_inv_sqrt[:, None] * weight * deg_inv_sqrt[None, :]
    eig_L = np.linalg.eigvalsh(L)
    eig_L_pos = eig_L[eig_L > 1e-9]
    if eig_L_pos.size < 2:
        return float("nan")
    return float(np.sum(1.0 / eig_L_pos) / n)


def per_seed_readouts(xis: list[np.ndarray],
                      calib: dict) -> list[dict]:
    rows = []
    for xi in xis:
        x = np.array(np.where(np.isfinite(xi), xi, 0.0), dtype=float, copy=True)
        np.fill_diagonal(x, 0.0)
        vals = x[x != 0]
        mean_off = float(np.mean(np.abs(vals))) if vals.size else 0.0
        alpha_raw = mean_off / (mean_off + 0.05)

        x_gram = x.copy()
        np.fill_diagonal(x_gram, 1.0)
        eig = normalised_gram_eigvals(x_gram)
        top1 = float(eig[0]) if eig.size else 0.0
        top2 = float(eig[1]) if eig.size > 1 else 0.0
        beta_raw = 1.0 - 1.0 / (1.0 + top1 * (xi.shape[0] ** 0.5))
        diff_raw = laplacian_diffusion_strength(x)

        alpha = alpha_raw * calib["alpha_scale"]
        beta = beta_raw * calib["beta_scale"]
        diff = diff_raw * calib["diff_scale"] if math.isfinite(diff_raw) else float("nan")
        gamma = 1.0 - alpha
        rows.append({
            "alpha_xi_raw": alpha_raw,
            "beta_pi_raw": beta_raw,
            "D_omega_raw": diff_raw,
            "spectral_gap": top1 - top2,
            "alpha_xi": alpha,
            "gamma_C1": gamma,
            "eps_sync2_C3": gamma / 2.0,
            "beta_pi": beta,
            "D_omega_C2": beta - gamma,
            "D_omega_lattice": diff,
            "C1_residual": 0.0,
            "C2_residual": diff - (beta - gamma) if math.isfinite(diff) else float("nan"),
            "C3_residual": 0.0,
        })
    return rows


def load_snapshot_xis(path: str) -> list[np.ndarray]:
    with np.load(path, allow_pickle=False) as z:
        snaps = z["edge_xi_snapshots"]
        last = snaps.shape[1] - 1
        return [np.asarray(snaps[s, last], dtype=float)
                for s in range(snaps.shape[0])]


def aggregate(rows: list[dict], regime: str, n_lat: int) -> dict:
    keys = [k for k in rows[0] if not k.endswith("_raw") or True]
    out = {"regime": regime, "n_lat": n_lat}
    for k in rows[0]:
        vals = [r[k] for r in rows
                if isinstance(r[k], (int, float)) and math.isfinite(r[k])]
        out[k] = float(np.median(vals)) if vals else float("nan")
    return out


def main() -> int:
    bundle = json.loads(open(READOUT_PATH, encoding="utf-8").read())
    calib = bundle["calibration"]

    existing_regimes = {r["regime"] for r in bundle["p5_ladder_per_N_readout"]}
    print("=== Extending per-N readout ===")
    print(f"Existing regimes in bundle: {sorted(existing_regimes)}")
    print()

    new_rows = []
    for regime, npz_path in NPZ_LOOKUP.items():
        if regime in existing_regimes:
            print(f"  {regime}: already present, skipping")
            continue
        if not os.path.exists(npz_path):
            print(f"  {regime}: NPZ not found at {npz_path}, skipping")
            continue
        n_lat = int(regime.replace("P5N", ""))
        xis = load_snapshot_xis(npz_path)
        per_seed = per_seed_readouts(xis, calib)
        row = aggregate(per_seed, regime, n_lat)
        new_rows.append(row)
        print(f"  {regime} (N={n_lat}, {len(xis)} seeds): "
              f"alpha={row['alpha_xi']:.4f}, beta={row['beta_pi']:.4f}, "
              f"D_omega={row['D_omega_lattice']:.4f}, "
              f"gamma={row['gamma_C1']:.4f}")

    if not new_rows:
        print("\nNo new rows to add.")
        return 0

    # Merge and sort by n_lat
    merged = bundle["p5_ladder_per_N_readout"] + new_rows
    merged.sort(key=lambda r: r["n_lat"])
    bundle["p5_ladder_per_N_readout"] = merged

    # Add note in method if not present
    note = (" Extended 2026-05-11 with P5N256 (12-seed) and P5N512 "
            "(12-seed) by running the same per-seed coefficient "
            "readout against the bundled snapshot NPZs.")
    if "Extended 2026-05-11" not in bundle.get("method", ""):
        bundle["method"] = bundle.get("method", "") + note

    # Write back to primary path and mirrors
    out_text = json.dumps(bundle, indent=2, sort_keys=False, ensure_ascii=False)
    for p in [READOUT_PATH] + MIRROR_PATHS:
        if os.path.isdir(os.path.dirname(p)):
            with open(p, "w", encoding="utf-8") as f:
                f.write(out_text)
            print(f"\nWrote {p}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
