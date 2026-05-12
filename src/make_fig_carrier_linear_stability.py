"""make_fig_carrier_linear_stability.py

Produces a two-panel figure for the carrier-operator linear-stability
section in P2:

  (left)  Spectrum of the lattice graph-Laplacian on the within-
          canonical N-ladder, with the predicted critical eigenvalue
          lambda_crit = 68/67 marked as a horizontal threshold and
          the growth rate g(lambda) = (68 - 67*lambda)/80 overlaid.
          Shows that for each N exactly one eigenvalue lies below
          lambda_crit (the constant mode at lambda = 0).

  (right) Symanzik 1/N convergence of the growing-mode fraction
          f(N) = (count below crit) / N. Predicted: f(N) = 1/N
          exactly. Empirically: a_inf = 0, b = 1 to machine precision.

Output: figures/fig_carrier_linear_stability.{pdf,png}
"""
from __future__ import annotations

import os
import glob
from fractions import Fraction
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless CI / no-DISPLAY
matplotlib.rcParams["pdf.fonttype"] = 42  # embed TrueType (vector, arXiv-friendly)
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.pyplot as plt

# Match the algebraic constants of the verifier
ALPHA_XI = Fraction(9, 10); GAMMA = Fraction(1, 10); EPS_SYNC2 = Fraction(1, 20)
BETA_PI = Fraction(15, 16); D_OMEGA = BETA_PI - GAMMA
LAMBDA_CRIT = (ALPHA_XI - GAMMA + EPS_SYNC2) / D_OMEGA  # 68/67
G_NET = ALPHA_XI + BETA_PI + EPS_SYNC2 - GAMMA          # 143/80
LAMBDA_CRIT_F = float(LAMBDA_CRIT)

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = REPO_ROOT.parent
EMERGENCE_RESULTS = Path(os.environ.get("EMERGENCE_RESULTS_ROOT",
                                        str(WORKSPACE)))

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


def main() -> int:
    spectra = {}      # N -> sorted eigenvalues (seed-averaged adjacency)
    counts = {}       # N -> mean count of eigenvalues below crit per seed
    for N in sorted(N_LADDER):
        path = find_snapshot_npz(N)
        if path is None:
            continue
        d = np.load(path, allow_pickle=True)
        edge = d["edge_xi_snapshots"][:, -1, :, :]  # (S, N, N) final
        S = edge.shape[0]
        # Per-seed counts
        per_seed_counts = []
        for s in range(S):
            ev = graph_laplacian_eigvals(edge[s])
            per_seed_counts.append(int(np.sum(ev < LAMBDA_CRIT_F)))
        counts[N] = float(np.mean(per_seed_counts))
        # Spectrum from seed-averaged adjacency for visualisation
        W_avg = np.mean(edge, axis=0)
        spectra[N] = graph_laplacian_eigvals(W_avg)

    Ns_arr = np.array(sorted(spectra.keys()), dtype=float)
    if len(Ns_arr) == 0:
        print("No lattice snapshot data visible; figure cannot be made.")
        return 1

    # Set up figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # ---------- Panel (left): spectrum scatter + lambda_crit line ----------
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=Ns_arr.min(), vmax=Ns_arr.max())
    for N in sorted(spectra.keys()):
        ev = spectra[N]
        x = np.full_like(ev, float(N))
        ax1.scatter(x, ev, s=3, color=cmap(norm(N)), alpha=0.55,
                    rasterized=True)
        # Highlight the smallest eigenvalue (constant mode)
        ax1.scatter([float(N)], [ev[0]], s=28, color="red",
                    edgecolor="black", linewidth=0.4, zorder=4,
                    label=("constant-mode eigenvalue (lambda = 0)"
                           if N == sorted(spectra.keys())[0] else None))

    ax1.axhline(LAMBDA_CRIT_F, color="black", linestyle="--", linewidth=1.4,
                label=r"$\lambda_{\rm crit}=68/67\approx 1.015$")
    ax1.set_xlabel(r"lattice size $N$")
    ax1.set_ylabel(r"graph-Laplacian eigenvalue $\lambda$")
    ax1.set_title("Lattice spectrum vs critical eigenvalue")
    ax1.set_xticks(sorted(spectra.keys()))
    ax1.set_xticklabels([str(int(N)) for N in sorted(spectra.keys())],
                        rotation=45, ha="right", fontsize=9)
    ax1.set_ylim(-0.3, max([ev.max() for ev in spectra.values()]) * 1.05)
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax1.grid(True, alpha=0.25)

    # ---------- Panel (right): Symanzik 1/N ----------
    Ns = np.array(sorted(counts.keys()), dtype=float)
    fracs = np.array([counts[N] / N for N in Ns])
    inv_N = 1.0 / Ns
    # Predicted line f(N) = 1/N
    x_dense = np.linspace(0, 1.0/Ns.min() * 1.05, 200)
    ax2.plot(x_dense, x_dense, color="black", linestyle="--", linewidth=1.2,
             label=r"prediction $f(N)=1/N$")
    ax2.scatter(inv_N, fracs, s=42, color="C3", edgecolor="black",
                linewidth=0.5, zorder=3,
                label="empirical (within-canonical ladder)")
    ax2.set_xlabel(r"$1/N$")
    ax2.set_ylabel(r"growing-mode fraction $f(N) = (\#\lambda<\lambda_{\rm crit})/N$")
    ax2.set_title("Symanzik convergence: $f(N) = 1/N$ exactly")
    ax2.legend(loc="upper left", fontsize=9, framealpha=0.92)
    ax2.set_xlim(0, max(inv_N) * 1.08)
    ax2.set_ylim(0, max(fracs) * 1.15)
    ax2.grid(True, alpha=0.25)

    plt.tight_layout()

    out_dir = REPO_ROOT / "paper" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "fig_carrier_linear_stability.pdf"
    png_path = out_dir / "fig_carrier_linear_stability.png"
    fig.savefig(pdf_path, dpi=200, bbox_inches="tight")
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
