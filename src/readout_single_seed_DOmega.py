"""One-shot D_Omega readout from a single-seed slim-runner snapshot NPZ.

Reads the existing calibration constants from causal_wave_per_N_readout.json,
runs the per-seed readout (alpha/beta/D_omega) on the snapshot's final
xi-matrix, and prints the result. Does NOT modify any bundled JSON.

Usage:
    python src/readout_single_seed_DOmega.py \\
        --snapshot ../results_d1_p5n600_trial1seed/P5N600.snapshots.npz
"""
from __future__ import annotations
import argparse
import json
import math
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
EMERGENCE = os.path.dirname(REPO)
READOUT_JSON = os.path.join(REPO, "data", "causal_wave_per_N_readout.json")


def laplacian_diffusion_strength(xi, xi_thresh=0.5):
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


def normalised_gram_eigvals(x):
    g = x @ x.T
    eig = np.linalg.eigvalsh(g)
    eig = eig[eig > 0]
    total = float(np.sum(eig))
    if total <= 1e-12:
        return np.zeros(x.shape[0])
    return np.sort(eig / total)[::-1]


def readout_single_xi(xi, calib):
    x = np.array(np.where(np.isfinite(xi), xi, 0.0), dtype=float, copy=True)
    np.fill_diagonal(x, 0.0)
    vals = x[x != 0]
    mean_off = float(np.mean(np.abs(vals))) if vals.size else 0.0
    alpha_raw = mean_off / (mean_off + 0.05)

    x_gram = x.copy()
    np.fill_diagonal(x_gram, 1.0)
    eig = normalised_gram_eigvals(x_gram)
    top1 = float(eig[0]) if eig.size else 0.0
    beta_raw = 1.0 - 1.0 / (1.0 + top1 * (xi.shape[0] ** 0.5))
    diff_raw = laplacian_diffusion_strength(x)

    alpha = alpha_raw * calib["alpha_scale"]
    beta = beta_raw * calib["beta_scale"]
    diff = diff_raw * calib["diff_scale"] if math.isfinite(diff_raw) else float("nan")
    gamma = 1.0 - alpha
    return {
        "alpha_xi": alpha, "beta_pi": beta, "gamma_C1": gamma,
        "eps_sync2_C3": gamma / 2.0,
        "D_omega_C2": beta - gamma,
        "D_omega_lattice": diff,
        "alpha_xi_raw": alpha_raw, "beta_pi_raw": beta_raw,
        "D_omega_raw": diff_raw,
    }


def v2_of(n):
    k, n = 0, int(n)
    while n > 0 and n % 2 == 0:
        n //= 2
        k += 1
    return k


def theta_chir(n, n_star=50, n_gen=3, d=4):
    if n <= 0:
        return 0.0
    x = math.log(n / n_star) / math.log(d * n_gen)
    return math.atan(n_gen ** (2 * x - 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", required=True,
                    help="Path to *.snapshots.npz")
    ap.add_argument("--snapshot-index", type=int, default=-1,
                    help="Snapshot index to use (-1 = last)")
    args = ap.parse_args()

    cfg = json.loads(open(READOUT_JSON, encoding="utf-8").read())
    calib = cfg["calibration"]

    z = np.load(args.snapshot, allow_pickle=False)
    snaps = z["edge_xi_snapshots"]
    n_seeds, n_snaps_per_seed = snaps.shape[0], snaps.shape[1]
    n_lat = int(snaps.shape[2])
    idx = args.snapshot_index if args.snapshot_index >= 0 \
            else n_snaps_per_seed - 1
    print(f"=== Single-seed D_Omega readout ===")
    print(f"NPZ: {args.snapshot}")
    print(f"n_lat={n_lat}, n_seeds={n_seeds}, n_snaps/seed={n_snaps_per_seed}")
    print(f"using snapshot index {idx} (of 0..{n_snaps_per_seed-1})")
    print()

    rows = []
    for s in range(n_seeds):
        xi = np.asarray(snaps[s, idx], dtype=float)
        r = readout_single_xi(xi, calib)
        r["seed"] = s
        rows.append(r)
        print(f"  seed {s}: alpha={r['alpha_xi']:.4f} beta={r['beta_pi']:.4f}"
              f" gamma={r['gamma_C1']:.4f} D_Omega={r['D_omega_lattice']:.4f}")

    if n_seeds > 1:
        keys = ["alpha_xi", "beta_pi", "gamma_C1", "D_omega_lattice"]
        medians = {k: float(np.median([r[k] for r in rows])) for k in keys}
        print()
        print("  median over seeds:")
        for k, v in medians.items():
            print(f"    {k:<18} = {v:.4f}")
    else:
        medians = {k: rows[0][k] for k in
                    ["alpha_xi", "beta_pi", "gamma_C1", "D_omega_lattice"]}

    # Compare against v_2-interaction model prediction
    v2 = v2_of(n_lat)
    th = theta_chir(n_lat)
    sin2 = math.sin(th) ** 2
    d_slow = (67/80) * math.cos(th)**2 + (math.pi/4) * sin2
    b0, b1 = -0.1613, -0.0522
    r_pred = b0 + b1 * sin2 * v2
    d_pred = d_slow + r_pred
    d_obs = medians["D_omega_lattice"]
    print()
    print(f"=== v_2-interaction model comparison (N={n_lat}) ===")
    print(f"  v_2({n_lat}) = {v2}")
    print(f"  theta_chir   = {math.degrees(th):.2f} deg")
    print(f"  D_slow       = {d_slow:.4f}")
    print(f"  r_pred       = {r_pred:+.4f}")
    print(f"  D_Omega_pred = {d_pred:.4f}")
    print(f"  D_Omega_obs  = {d_obs:.4f}")
    print(f"  difference   = {d_obs - d_pred:+.4f}")
    print(f"  rel_err      = {((d_obs - d_pred) / d_pred * 100):+.2f}%"
          if abs(d_pred) > 1e-6 else "  rel_err      = (D_pred ~ 0)")


if __name__ == "__main__":
    main()
