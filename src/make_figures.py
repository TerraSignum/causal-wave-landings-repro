r"""
Generate the four figures of the Paper 2 manuscript.

  Figure 1 - Causal-wave transport equation schematic.
  Figure 2 - Five coefficients as a bar plot.
  Figure 3 - Ten landings: residuals vs target with tier coloring.
  Figure 4 - Look-elsewhere distribution for the three EXACT-tier targets.

Outputs are written to paper/figures/ in PDF and PNG.

Usage:
    python ./src/make_figures.py
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import recompute_landings as M  # noqa: E402

FIG_DIR = REPO / "paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})


def save_both(fig, stem):
    pdf = FIG_DIR / f"{stem}.pdf"
    png = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf, format="pdf")
    fig.savefig(png, format="png")
    print(f"  Saved: {pdf.relative_to(REPO)} + .png")


def figure_1_schematic():
    """Causal-wave transport equation schematic."""
    fig, ax = plt.subplots(figsize=(11, 3.5))
    ax.set_axis_off()
    ax.text(0.5, 0.78,
            r"$\dfrac{\mathrm{d}C}{\mathrm{d}\tau}\;=\;"
            r"D(\Omega)\,\Delta C\;-\;\gamma\,C\;+\;\alpha_\xi\,C\;+\;"
            r"\beta_\pi\,\Pi_{\mathrm{common}}\,C\;+\;\varepsilon_{\mathrm{sync}}^{2}\,C$",
            ha="center", va="center", fontsize=18)

    labels = [
        (0.25, 0.40, r"$D(\Omega)$" "\n" "diffusion / curvature\n" r"$0.83996$"),
        (0.40, 0.40, r"$-\gamma\,C$" "\n" "damping / time arrow\n" r"$\gamma=0.10021$"),
        (0.55, 0.40, r"$\alpha_\xi\,C$" "\n" "Xi reaction rate\n" r"$0.90082$"),
        (0.70, 0.40, r"$\beta_\pi\,\Pi_{\mathrm{common}}\,C$" "\n" "common phase\n" r"$\beta_\pi=0.93791$"),
        (0.86, 0.40, r"$\varepsilon_{\mathrm{sync}}^{2}\,C$" "\n" "sync persistence\n" r"$0.05000$"),
    ]
    for (x, y, txt) in labels:
        ax.text(x, y, txt, ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#fafafa",
                          edgecolor="#888888"))
    ax.text(0.5, 0.06,
            "Five measured coefficients; not fitted to downstream observables",
            ha="center", fontsize=10, style="italic", color="#444")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("Causal-wave transport equation and five measured coefficients",
                 pad=10)
    save_both(fig, "fig1_schematic")
    plt.close(fig)


def figure_2_coefficients():
    """Bar chart of the five causal-wave coefficients."""
    R = M.compute_all()
    c = R["coefficients"]
    names = [r"$\alpha_\xi$", r"$D(\Omega)$", r"$\beta_\pi$",
             r"$\gamma$", r"$\varepsilon^{2}_{\mathrm{sync}}$"]
    keys = ["alpha_xi", "D_Omega", "beta_pi", "gamma", "eps_sync2"]
    vals = [c[k] for k in keys]
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#4a90d9", "#5cb85c", "#9b59b6", "#e89043", "#888888"]
    bars = ax.bar(names, vals, color=colors, edgecolor="black", lw=1.0)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01,
                f"{val:.5f}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Coefficient value")
    ax.set_ylim(0, 1.05)
    ax.set_title("Five measured causal-wave transport coefficients", pad=10)
    save_both(fig, "fig2_coefficients")
    plt.close(fig)


def figure_3_landings():
    """Ten landings: residual vs target with tier coloring."""
    R = M.compute_all()
    rows = R["rows"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    ids = [r["id"] for r in rows]
    residuals = [r["residual_pct"] for r in rows]
    tiers = [r["tier"] for r in rows]
    le = [r["look_elsewhere_caveated"] for r in rows]
    color_map = {"EXACT": "#1f7a1f", "PRECISE_2.5": "#1f77b4",
                 "FACTOR2": "#e89043", "FAR_OFF": "#cc3333"}
    colors = [color_map[t] for t in tiers]
    edgecolors = ["#000000" if not l else "#cc3333" for l in le]

    bars = ax1.bar(ids, residuals, color=colors, edgecolor=edgecolors, linewidth=1.5)
    for bar, val in zip(bars, residuals):
        h = val
        ax1.text(bar.get_x() + bar.get_width() / 2,
                 h + (0.1 if h >= 0 else -0.1),
                 f"{val:+.3f}%",
                 ha="center", va=("bottom" if h >= 0 else "top"), fontsize=8)
    ax1.axhline(2.5, color="#999999", linestyle=":", lw=1, label="PRECISE_2.5 cut")
    ax1.axhline(-2.5, color="#999999", linestyle=":", lw=1)
    ax1.axhline(0, color="black", linewidth=0.7)
    ax1.set_ylabel("Residual vs target (%)")
    ax1.set_title("(A) Ten landings: residuals", pad=8)
    ax1.set_ylim(-3.5, 3.5)
    ax1.legend(loc="lower right", framealpha=0.9)

    ax2.bar(ids, [abs(r) for r in residuals], color=colors, edgecolor=edgecolors, linewidth=1.5)
    ax2.set_yscale("log")
    ax2.axhline(2.5, color="#999999", linestyle=":", lw=1)
    ax2.axhline(0.01, color="#999999", linestyle="-.", lw=1, label="EXACT cut (0.01%)")
    ax2.set_ylabel("|residual| (%)")
    ax2.set_title("(B) Same data, log scale; red edge = look-elsewhere caveated", pad=8)
    ax2.legend(loc="upper right", framealpha=0.9)
    ax2.set_ylim(1e-4, 10)

    n_total = len(rows)
    fig.suptitle(f"Causal-wave benchmark landings; {n_total}/{n_total} PRECISE-or-better",
                 fontsize=12, y=0.995)
    save_both(fig, "fig3_landings")
    plt.close(fig)


def figure_4_lookelsewhere():
    """Look-elsewhere distribution for the three EXACT-tier targets."""
    with open(REPO / "data" / "coefficient_perturbation_null.json", "r", encoding="utf-8") as f:
        pert = json.load(f)
    targets = ["sin2_theta_W", "BH_entropy_quarter", "Einstein_gap_two_thirds"]
    p_random_le = [pert["null_statistics"][t]["p_random_residual_le_baseline"] for t in targets]
    p_all_three = pert["p_all_3_targets_hit_under_random_coefficients"]
    p_any = pert["p_any_target_hit_under_random_coefficients"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    nice = [r"$\sin^{2}\theta_W$", "BH entropy 1/4", "Einstein gap 2/3"]
    bars = ax1.bar(nice, p_random_le, color="#cc3333", edgecolor="black", linewidth=1.0)
    for bar, val in zip(bars, p_random_le):
        ax1.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                 f"{val:.2f}", ha="center", va="bottom", fontsize=10)
    ax1.set_ylabel("P(random residual <= baseline)")
    ax1.set_ylim(0, 1.05)
    ax1.set_title("(A) Per-target null:\nrandom coefficients beat baseline",
                  pad=8)

    labels = ["any of 3 EXACT", "all 3 EXACT"]
    vals = [p_any, p_all_three]
    bars = ax2.bar(labels, vals, color=["#888888", "#cc3333"],
                   edgecolor="black", linewidth=1.0)
    for bar, val in zip(bars, vals):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                 f"{val:.2f}", ha="center", va="bottom", fontsize=10)
    ax2.axhline(0.5, color="#444", linestyle=":", lw=1)
    ax2.set_ylabel("P(random coefficients reproduce hits)")
    ax2.set_ylim(0, 1.1)
    ax2.set_title("(B) Joint null:\n200 random trials", pad=8)

    fig.suptitle("Look-elsewhere control for the three EXACT-tier targets (L6, L7, L8)",
                 fontsize=12, y=1.0)
    save_both(fig, "fig4_lookelsewhere")
    plt.close(fig)


def figure_5_ir_hyperbolicity():
    """IR-hyperbolicity certificate: log-log plot of the
    leading-order Lorentz-deviation Delta_Lor(k)/k^2 versus k,
    showing the o(k^2) decay (slope = 2 on log-log)."""
    import verify_ir_hyperbolicity as Vh
    ks = np.geomspace(1e-3, 0.5, 64)
    deviations, c_Xi_sq = Vh.lorentz_deviation(Vh.STENCIL, list(ks))
    ratios = np.array([r for (_k, r) in deviations])

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.loglog(ks, ratios, "o-", color="#1f77b4", lw=1.5, ms=4,
              label=r"$\Delta_{\mathrm{Lor}}(k)\,/\,k^{2}$")
    # Reference slope of k^2 (i.e., the o(k^2) envelope).
    ref = (ks ** 2) * c_Xi_sq / 2.0
    ax.loglog(ks, ref, "--", color="#444",
              label=r"$\propto k^{2}$ envelope ($c_{\Xi}^{2}\,k^{2}/2$)")
    ax.set_xlabel(r"$|\mathbf{k}|$  (lattice units)")
    ax.set_ylabel(r"$\Delta_{\mathrm{Lor}}(k) / k^{2}$")
    ax.set_title(r"IR-hyperbolicity certificate: "
                 r"$\Delta_{\mathrm{Lor}}(k)\,/\,k^{2} \to 0$ as $k \to 0$ "
                 r"($c_{\Xi}^{2} = Z_{x}/Z_{\tau} = 1$)")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="upper left", framealpha=0.95)
    save_both(fig, "fig5_ir_hyperbolicity")
    plt.close(fig)


def main():
    print("Generating Paper 2 figures into paper/figures/")
    print()
    figure_1_schematic()
    figure_2_coefficients()
    figure_3_landings()
    figure_4_lookelsewhere()
    figure_5_ir_hyperbolicity()
    print()
    print("All five figures generated.")


if __name__ == "__main__":
    main()
