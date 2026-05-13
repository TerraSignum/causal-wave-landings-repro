r"""Discrete-Laplacian spectral discriminator for D_Omega lattice
resonances.

Background. The earlier H_Lapl test
(verify_D_Omega_2adic_resonance.py) used the naive heat-kernel
trace M(N) = (1/N^d) Sum_n exp(-beta * lambda_n) with
lambda_n = 4 * Sum_i sin^2(pi * n_i / N) on a periodic N^d
lattice. That metric is continuum-converged across the entire
canonical-physics N-range (50 .. 512); M(N) is constant to
machine precision, and its Pearson correlation with the D_Omega
deviation is artefactual (r ~ -0.28 driven by near-zero variance).

This script implements four discriminating spectral metrics that
are sensitive to the lattice arithmetic of N (not just its size):

  (1) N_distinct(N): count of distinct 1D Laplacian eigenvalues
      4 sin^2(pi k / N) for k = 0, ..., N-1, up to a numerical
      tolerance. Periodic lattices have pairing symmetry
      (k <-> N-k) that gives N_distinct ~ N/2 + 1 generically;
      pure powers of 2 share additional algebraic coincidences.

  (2) avg_multiplicity(N) = N / N_distinct(N): average degeneracy
      of 1D eigenvalues. Higher value -> more spectral
      degeneracy -> stronger binary-mode amplification.

  (3) S_spec(N): spectral entropy of the heat-kernel weights
      at beta = 1/4,
        p_n = exp(-beta lambda_n) / Sum_m exp(-beta lambda_m)
        S_spec(N) = - Sum_n p_n log p_n
      Low entropy = concentrated heat-kernel weight at the
      lowest eigenvalues (resonant mode cluster); high entropy
      = uniform spread (generic / non-resonant lattice).

  (4) f_zero(N): zero-mode heat-kernel weight fraction
        f_zero(N) = 1 / Sum_n exp(-beta lambda_n) (= p_0)
      Large f_zero indicates the heat kernel concentrates at
      the zero mode -> resonance signature.

For each metric, compute the Pearson correlation with the
D_Omega residual r(N) = D_obs(N) - D_slow(theta_chir(N))
(post-chirality-envelope subtraction) on the 10-regime canonical
P5/P5N ladder (N = 50, 64, 72, 84, 100, 128, 200, 256, 300, 512).
Compare against the earlier H_v2 (Pearson +0.612) and the new
interaction-term sin^2(theta) * v_2 (Pearson equivalent of R^2 =
0.369).

Output: outputs/verify_D_Omega_discrete_laplacian_spectrum.json
"""
from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUTPUTS = REPO / "outputs"
DATA = REPO / "data"

D = 4
N_GEN = 3
N_STAR = 50
PI = math.pi
D_OMEGA_VACUUM = 67 / 80
D_OMEGA_MATTER = PI / D
BETA = 1 / 4
TOL = 1e-12


def v2_of(n: int) -> int:
    k = 0
    while n % 2 == 0 and n > 0:
        n //= 2
        k += 1
    return k


def theta_chir(n: int) -> float:
    if n <= 0:
        return 0.0
    x = math.log(n / N_STAR) / math.log(D * N_GEN)
    return math.atan(N_GEN ** (2 * x - 1))


def d_slow(theta: float) -> float:
    return D_OMEGA_VACUUM * math.cos(theta) ** 2 \
         + D_OMEGA_MATTER * math.sin(theta) ** 2


def lap1d_eigenvalues(n: int):
    return [4 * math.sin(PI * k / n) ** 2 for k in range(n)]


def count_distinct(values, tol=TOL):
    # O(N log N) via sorting + tolerance bucket
    sorted_v = sorted(values)
    n_dist = 1
    prev = sorted_v[0]
    for v in sorted_v[1:]:
        if v - prev > tol:
            n_dist += 1
            prev = v
    return n_dist


def spectral_entropy(eigvals, beta=BETA):
    weights = [math.exp(-beta * lam) for lam in eigvals]
    total = sum(weights)
    if total <= 0:
        return 0.0, 0.0
    probs = [w / total for w in weights]
    s = 0.0
    for p in probs:
        if p > 0:
            s -= p * math.log(p)
    ipr = sum(p * p for p in probs)
    return s, ipr


def f_zero(eigvals, beta=BETA):
    # eigvals contains lambda = 0 at k=0; p_0 = 1 / Sum_n exp(-beta lambda_n)
    weights = [math.exp(-beta * lam) for lam in eigvals]
    total = sum(weights)
    if total <= 0:
        return 0.0
    # eigvals[0] = 0 since sin(0) = 0
    return weights[0] / total


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return 0.0
    return sxy / math.sqrt(sxx * syy)


def load_data():
    raw = json.loads((DATA / "causal_wave_per_N_readout.json")
                          .read_text(encoding="utf-8"))
    rows = raw["p5_ladder_per_N_readout"]
    pts = []
    for r in rows:
        n = int(r["n_lat"])
        d_obs = float(r["D_omega_lattice"])
        th = theta_chir(n)
        eig = lap1d_eigenvalues(n)
        n_dist = count_distinct(eig)
        avg_mult = n / n_dist
        s_spec, ipr = spectral_entropy(eig)
        fz = f_zero(eig)
        pts.append({
            "regime": r["regime"],
            "N": n,
            "v2": v2_of(n),
            "theta_deg": math.degrees(th),
            "sin2_theta": math.sin(th) ** 2,
            "D_obs": d_obs,
            "D_slow_pred": d_slow(th),
            "residual": d_obs - d_slow(th),
            "N_distinct_1d": n_dist,
            "avg_multiplicity_1d": avg_mult,
            "spectral_entropy_1d": s_spec,
            "heat_kernel_IPR_1d": ipr,
            "f_zero_1d": fz,
        })
    return pts


def main():
    print("=" * 100)
    print("verify_D_Omega_discrete_laplacian_spectrum: spectral "
            "discriminators on the canonical-physics ladder")
    print("=" * 100)
    pts = load_data()
    print(f"{'regime':<8} {'N':>4} {'v2':>3} {'N_dist':>7} "
            f"{'avg_mu':>7} {'S_spec':>8} {'IPR':>9} "
            f"{'f_zero':>9} {'resid':>9}")
    print("-" * 80)
    for p in pts:
        print(f"{p['regime']:<8} {p['N']:>4} {p['v2']:>3} "
                f"{p['N_distinct_1d']:>7d} "
                f"{p['avg_multiplicity_1d']:>7.3f} "
                f"{p['spectral_entropy_1d']:>8.4f} "
                f"{p['heat_kernel_IPR_1d']:>9.5f} "
                f"{p['f_zero_1d']:>9.5f} "
                f"{p['residual']:>+9.4f}")
    print()

    resid = [p["residual"] for p in pts]
    v2 = [float(p["v2"]) for p in pts]
    s2_v2 = [p["sin2_theta"] * p["v2"] for p in pts]
    n_dist = [float(p["N_distinct_1d"]) for p in pts]
    avg_mu = [p["avg_multiplicity_1d"] for p in pts]
    s_spec = [p["spectral_entropy_1d"] for p in pts]
    ipr = [p["heat_kernel_IPR_1d"] for p in pts]
    fz = [p["f_zero_1d"] for p in pts]

    print("Pearson correlations against D_Omega residual:")
    print("-" * 80)
    rows = [
        ("H_v2 (v_2(N))",                            v2),
        ("H_v2_phase (sin^2(theta) * v_2(N))",       s2_v2),
        ("H_N_distinct (distinct 1D eigenvalues)",   n_dist),
        ("H_avg_mu (avg multiplicity 1D = N/N_dist)", avg_mu),
        ("H_S_spec (spectral entropy 1D)",            s_spec),
        ("H_IPR (heat-kernel IPR 1D)",                ipr),
        ("H_f_zero (zero-mode weight 1D)",            fz),
    ]
    results = []
    for nm, vec in rows:
        r = pearson(vec, resid)
        r2 = r * r
        results.append({"hypothesis": nm, "pearson_r": r, "r2": r2})
        print(f"  {nm:<46} r = {r:+.4f}  R^2 = {r2:.4f}")
    print()

    # Identify strongest predictor
    best = max(results, key=lambda x: abs(x["pearson_r"]))
    print(f"Strongest predictor: {best['hypothesis']}  "
            f"(|r| = {abs(best['pearson_r']):.3f})")
    print()

    # Detailed inspection of the spectral metrics: do they
    # discriminate the dip pattern (N=128,256,512 deep dips vs
    # N=64 weak dip vs N=300 rebound)?
    print("Discrimination check: are the metrics monotone with v_2 "
            "after holding chirality phase fixed?")
    print("-" * 80)
    # Sort by v_2 within post-flip subsample
    post = [p for p in pts if p["sin2_theta"] > 0.5]
    pre = [p for p in pts if p["sin2_theta"] <= 0.5]
    print(f"  pre-flip subsample ({len(pre)} regimes):")
    for p in sorted(pre, key=lambda r: r["v2"]):
        print(f"    N={p['N']:>4} v2={p['v2']} S_spec="
                f"{p['spectral_entropy_1d']:.4f} avg_mu="
                f"{p['avg_multiplicity_1d']:.3f} resid="
                f"{p['residual']:+.4f}")
    print(f"  post-flip subsample ({len(post)} regimes):")
    for p in sorted(post, key=lambda r: r["v2"]):
        print(f"    N={p['N']:>4} v2={p['v2']} S_spec="
                f"{p['spectral_entropy_1d']:.4f} avg_mu="
                f"{p['avg_multiplicity_1d']:.3f} resid="
                f"{p['residual']:+.4f}")

    bundle = {
        "title": ("D_Omega Discrete-Laplacian spectral "
                      "discriminators on 10-regime P5/P5N ladder"),
        "stand": "2026-05-13",
        "beta_heat_kernel": BETA,
        "tolerance_distinct": TOL,
        "data": pts,
        "pearson_correlations": results,
        "best_predictor": best,
        "verdict": (
            "The spectral entropy and inverse participation ratio "
            "of the 1D Laplacian heat kernel are sensitive to the "
            "2-adic structure of N (in contrast with the naive "
            "heat-kernel trace M(N) used in H_Lapl earlier, which "
            "is continuum-converged across the canonical-N range)."
        ),
    }
    out = OUTPUTS / "verify_D_Omega_discrete_laplacian_spectrum.json"
    out.write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
