r"""Full chirality-harmonic-basis closure of D_Omega(N),
extending the slow-envelope ansatz with the sin(2 theta) and
sin(4 theta) endpoint-preserving harmonics that already close
the matter-localised K, Q factor fields (cf. K+Q full-closure
audit 2026-05-06; iter-37 final; eight clean System-R rationals
at <1 sigma).

The K, Q harmonic closure is:

    F(N) = F_pre cos^2(theta) + F_post sin^2(theta)
             + a sin(2 theta) + b sin(4 theta)

with eight System-R rationals (Memory entry K + Q full closure
2026-05-06):

  K_pre  = 4/3 - gamma^2/3      = 133/100   (0.05 sigma)
  K_post = 4/3 + gamma^3                     (0.40 sigma)
  a_K    = -gamma^2/d           = -1/400    (0.10 sigma)
  b_K    = +gamma/16            = +1/160    (0.12 sigma)
  Q_pre  = 1/4 + gamma^2 (N_gen + d) = 8/25  (0.29 sigma)
  Q_post = 1/4 + gamma / (2 N_gen d) = 61/240 (0.09 sigma)
  a_Q    = -2 gamma^2           = -1/50     (0.30 sigma)
  b_Q    = -9 gamma^2/8         = -9/800    (0.02 sigma)

The D_Omega(N) model used so far is just the slow envelope:

    D_slow(theta) = (67/80) cos^2(theta) + (pi/d) sin^2(theta)
    D_Omega(N) ~ D_slow + lattice-resonance(v_2(N), sin^2 theta)

This script tests whether D_Omega also closes in the natural
harmonic basis, i.e. whether the residual

    r(N) = D_obs(N) - D_slow(theta_chir(N))

is well-fitted by the harmonic ansatz

    r(N) = a_D sin(2 theta_chir(N)) + b_D sin(4 theta_chir(N))
             + c_D sin^2(theta_chir(N)) v_2(N)

If the two harmonic coefficients (a_D, b_D) capture the N=84 /
N=128 excess dips that remain after the v_2 * sin^2 interaction
fit alone, then the "missed correlation" is just the K, Q-style
harmonic basis applied to D_Omega.

If they recover clean System-R rationals (e.g. a_D, b_D in
{gamma, gamma^2/d, gamma/16, ...}), then the closure extends
the established 8-rational K, Q pattern with two additional
rationals for D_Omega.

Output: outputs/verify_D_Omega_full_harmonic_closure.json
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
GAMMA = 1 / 10
PI = math.pi
D_OMEGA_VACUUM = 67 / 80
D_OMEGA_MATTER = PI / D


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


def ols(X, y):
    n = len(y)
    p = len(X[0])
    A = [[sum(X[r][i] * X[r][j] for r in range(n))
              for j in range(p)] + [sum(X[r][i] * y[r]
                                          for r in range(n))]
            for i in range(p)]
    for i in range(p):
        piv = A[i][i]
        if abs(piv) < 1e-14:
            for k in range(i + 1, p):
                if abs(A[k][i]) > 1e-14:
                    A[i], A[k] = A[k], A[i]
                    piv = A[i][i]
                    break
        for j in range(p + 1):
            A[i][j] /= piv
        for k in range(p):
            if k != i:
                f = A[k][i]
                for j in range(p + 1):
                    A[k][j] -= f * A[i][j]
    beta = [A[i][p] for i in range(p)]
    y_pred = [sum(X[k][j] * beta[j] for j in range(p))
                for k in range(n)]
    rss = sum((y[k] - y_pred[k]) ** 2 for k in range(n))
    ybar = sum(y) / n
    tss = sum((y[k] - ybar) ** 2 for k in range(n))
    r2 = 1 - rss / tss if tss > 0 else 0.0
    aic = n * math.log(rss / n) + 2 * p if rss > 0 else float("-inf")
    return {"beta": beta, "rss": rss, "r2": r2, "aic": aic,
              "y_pred": y_pred, "n": n, "p": p}


def loo_cv(y_full, build_X_fn):
    n = len(y_full)
    sq_errs = []
    for k in range(n):
        X_tr = [build_X_fn(i) for i in range(n) if i != k]
        y_tr = [y_full[i] for i in range(n) if i != k]
        res = ols(X_tr, y_tr)
        x_te = build_X_fn(k)
        y_hat = sum(x_te[j] * res["beta"][j]
                      for j in range(len(x_te)))
        sq_errs.append((y_tr[0] if len(y_full) == 1
                           else (y_full[k] - y_hat)) ** 2)
    rmse = math.sqrt(sum(sq_errs) / n)
    return {"rmse_cv": rmse}


def closest_rational(value, max_q=2000):
    best = (1, 1, abs(value - 1.0))
    for q in range(1, max_q + 1):
        p = round(value * q)
        if p == 0:
            continue
        err = abs(value - p / q)
        if err < best[2]:
            best = (p, q, err)
    return best


def system_r_label(value):
    """Try to identify the value as a simple System-R rational
    expression in (gamma, N_gen, d)."""
    candidates = [
        (GAMMA, "gamma"),
        (-GAMMA, "-gamma"),
        (GAMMA / 2, "gamma/2"),
        (-GAMMA / 2, "-gamma/2"),
        (GAMMA / 4, "gamma/4"),
        (-GAMMA / 4, "-gamma/4"),
        (GAMMA / 8, "gamma/8"),
        (GAMMA / 16, "gamma/16"),
        (-GAMMA / 16, "-gamma/16"),
        (GAMMA ** 2, "gamma^2"),
        (-GAMMA ** 2, "-gamma^2"),
        (GAMMA ** 2 / 2, "gamma^2/2"),
        (-GAMMA ** 2 / 2, "-gamma^2/2"),
        (GAMMA ** 2 / D, "gamma^2/d"),
        (-GAMMA ** 2 / D, "-gamma^2/d"),
        (2 * GAMMA ** 2, "2 gamma^2"),
        (-2 * GAMMA ** 2, "-2 gamma^2"),
        (GAMMA ** 2 * (N_GEN + D), "gamma^2 (N_gen+d)"),
        (PI / D, "pi/d"),
        (-PI / D, "-pi/d"),
        (1 / D, "1/d = 1/4"),
        (1 / (2 * D), "1/(2d)"),
        (1 / 8, "1/8"),
        (-1 / 8, "-1/8"),
        (1 / 16, "1/16"),
        (-1 / 16, "-1/16"),
    ]
    best = None
    for v, name in candidates:
        err_rel = abs((value - v) / v) if v != 0 else abs(value)
        if best is None or err_rel < best[2]:
            best = (v, name, err_rel)
    return best


def load_data():
    raw = json.loads((DATA / "causal_wave_per_N_readout.json")
                          .read_text(encoding="utf-8"))
    rows = raw["p5_ladder_per_N_readout"]
    pts = []
    for r in rows:
        n = int(r["n_lat"])
        d_obs = float(r["D_omega_lattice"])
        th = theta_chir(n)
        pts.append({
            "regime": r["regime"], "N": n,
            "v2": v2_of(n),
            "theta": th,
            "theta_deg": math.degrees(th),
            "cos2_theta": math.cos(th) ** 2,
            "sin2_theta": math.sin(th) ** 2,
            "sin_2theta": math.sin(2 * th),
            "sin_4theta": math.sin(4 * th),
            "D_obs": d_obs,
            "D_slow_pred": d_slow(th),
            "residual": d_obs - d_slow(th),
        })
    return pts


def main():
    print("=" * 100)
    print("verify_D_Omega_full_harmonic_closure: K+Q-style harmonic basis "
            "for D_Omega(N)")
    print("=" * 100)
    pts = load_data()
    print(f"{'regime':<8} {'N':>4} {'v2':>3} {'theta':>7} "
            f"{'sin(2t)':>9} {'sin(4t)':>9} {'D_obs':>7} "
            f"{'D_slow':>7} {'resid':>9}")
    print("-" * 80)
    for p in pts:
        print(f"{p['regime']:<8} {p['N']:>4} {p['v2']:>3} "
                f"{p['theta_deg']:>7.2f} "
                f"{p['sin_2theta']:>+9.4f} "
                f"{p['sin_4theta']:>+9.4f} "
                f"{p['D_obs']:>7.4f} "
                f"{p['D_slow_pred']:>7.4f} "
                f"{p['residual']:>+9.4f}")
    print()

    y = [p["residual"] for p in pts]

    # Model A: v_2 * sin^2 interaction only (previous winner, commit
    # bde4c71)
    X_A = [[1.0, p["sin2_theta"] * p["v2"]] for p in pts]
    res_A = ols(X_A, y)

    # Model B: harmonic-basis-only (no v_2)
    X_B = [[1.0, p["sin_2theta"], p["sin_4theta"]] for p in pts]
    res_B = ols(X_B, y)

    # Model C: full = interaction + harmonics
    X_C = [[1.0, p["sin2_theta"] * p["v2"],
              p["sin_2theta"], p["sin_4theta"]] for p in pts]
    res_C = ols(X_C, y)

    # Model D: harmonics only (no constant intercept; pure K+Q
    # closure shape on the residual)
    X_D = [[p["sin_2theta"], p["sin_4theta"]] for p in pts]
    res_D = ols(X_D, y)

    print("Model comparison (residual = D_obs - D_slow):")
    print("-" * 80)
    print(f"{'model':<48} {'p':>2} {'R^2':>7} {'AIC':>9} "
            f"{'beta':>8}")

    def loo(build_fn):
        return loo_cv(y, build_fn)["rmse_cv"]

    def b_A(i):
        p = pts[i]
        return [1.0, p["sin2_theta"] * p["v2"]]

    def b_B(i):
        p = pts[i]
        return [1.0, p["sin_2theta"], p["sin_4theta"]]

    def b_C(i):
        p = pts[i]
        return [1.0, p["sin2_theta"] * p["v2"],
                  p["sin_2theta"], p["sin_4theta"]]

    def b_D(i):
        p = pts[i]
        return [p["sin_2theta"], p["sin_4theta"]]

    loo_A = loo(b_A)
    loo_B = loo(b_B)
    loo_C = loo(b_C)
    loo_D = loo(b_D)

    print(f"  A) v_2 * sin^2 only (prior winner)            "
            f"{res_A['p']:>2} {res_A['r2']:>+.3f} {res_A['aic']:>+.3f} "
            f"LOO={loo_A:.4f}")
    print(f"  B) harmonics-only b0 + a sin(2t) + b sin(4t)  "
            f"{res_B['p']:>2} {res_B['r2']:>+.3f} {res_B['aic']:>+.3f} "
            f"LOO={loo_B:.4f}")
    print(f"  C) full = interaction + harmonics              "
            f"{res_C['p']:>2} {res_C['r2']:>+.3f} {res_C['aic']:>+.3f} "
            f"LOO={loo_C:.4f}")
    print(f"  D) pure harmonics (no intercept)               "
            f"{res_D['p']:>2} {res_D['r2']:>+.3f} {res_D['aic']:>+.3f} "
            f"LOO={loo_D:.4f}")
    print()

    # Report fitted coefficients with rational identification
    print("Fitted coefficients with closest System-R rational match:")
    print("-" * 80)
    print(f"Model B (harmonics-only with intercept):")
    names_B = ["intercept", "a_D (sin 2 theta)", "b_D (sin 4 theta)"]
    for nm, b in zip(names_B, res_B["beta"]):
        v_r, name_r, err = system_r_label(b)
        p_r, q_r, err_r = closest_rational(b, 2000)
        print(f"  {nm:<22} = {b:+.5f}  closest: {p_r}/{q_r} "
                f"(rel err {err_r/abs(b)*100:.1f}%); "
                f"System-R: {name_r} (rel err {err*100:.1f}%)")
    print()
    print(f"Model D (harmonics-only without intercept):")
    names_D = ["a_D (sin 2 theta)", "b_D (sin 4 theta)"]
    for nm, b in zip(names_D, res_D["beta"]):
        v_r, name_r, err = system_r_label(b)
        p_r, q_r, err_r = closest_rational(b, 2000)
        print(f"  {nm:<22} = {b:+.5f}  closest: {p_r}/{q_r} "
                f"(rel err {err_r/abs(b)*100:.1f}%); "
                f"System-R: {name_r} (rel err {err*100:.1f}%)")
    print()
    print(f"Model C (full = interaction + harmonics):")
    names_C = ["intercept", "c_D (sin^2 v_2)",
                "a_D (sin 2 theta)", "b_D (sin 4 theta)"]
    for nm, b in zip(names_C, res_C["beta"]):
        v_r, name_r, err = system_r_label(b)
        p_r, q_r, err_r = closest_rational(b, 2000)
        print(f"  {nm:<22} = {b:+.5f}  closest: {p_r}/{q_r} "
                f"(rel err {err_r/abs(b)*100:.1f}%); "
                f"System-R: {name_r} (rel err {err*100:.1f}%)")
    print()

    # Residuals per regime for the best model
    aics = {"A_interaction": res_A["aic"], "B_harmonic_full": res_B["aic"],
              "C_combined": res_C["aic"], "D_harmonic_pure": res_D["aic"]}
    loos = {"A_interaction": loo_A, "B_harmonic_full": loo_B,
              "C_combined": loo_C, "D_harmonic_pure": loo_D}
    winner_aic = min(aics, key=lambda k: aics[k])
    winner_loo = min(loos, key=lambda k: loos[k])
    print(f"Winners: AIC-min = {winner_aic}, LOO-RMSE-min = {winner_loo}")
    print()

    # Report per-regime fits
    print(f"Per-regime fits (winner = {winner_aic}):")
    print(f"{'N':>4} {'resid':>9} {'fit':>9} {'sub-resid':>11}")
    winner_res = {"A_interaction": res_A, "B_harmonic_full": res_B,
                     "C_combined": res_C,
                     "D_harmonic_pure": res_D}[winner_aic]
    for p, yhat in zip(pts, winner_res["y_pred"]):
        sr = p["residual"] - yhat
        print(f"{p['N']:>4} {p['residual']:>+9.4f} "
                f"{yhat:>+9.4f} {sr:>+11.4f}")

    bundle = {
        "title": ("D_Omega harmonic-basis closure (K+Q-style) "
                      "with sin(2 theta), sin(4 theta)"),
        "stand": "2026-05-13",
        "data": pts,
        "model_comparison": {
            "A_interaction_only": {
                "r2": res_A["r2"], "aic": res_A["aic"],
                "beta": res_A["beta"], "loo_rmse": loo_A,
            },
            "B_harmonic_with_intercept": {
                "r2": res_B["r2"], "aic": res_B["aic"],
                "beta": res_B["beta"], "loo_rmse": loo_B,
            },
            "C_full": {
                "r2": res_C["r2"], "aic": res_C["aic"],
                "beta": res_C["beta"], "loo_rmse": loo_C,
            },
            "D_harmonic_pure": {
                "r2": res_D["r2"], "aic": res_D["aic"],
                "beta": res_D["beta"], "loo_rmse": loo_D,
            },
        },
        "winners": {"aic": winner_aic, "loo": winner_loo},
    }
    out = OUTPUTS / "verify_D_Omega_full_harmonic_closure.json"
    out.write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
