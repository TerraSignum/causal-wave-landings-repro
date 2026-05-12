# causal-wave-landings-repro

**Five-coefficient causal-wave transport law and ten cross-sector benchmark
landings.**

[![CI: reproduce](https://github.com/[anonymized]/causal-wave-landings-repro/actions/workflows/reproduce.yml/badge.svg)](https://github.com/[anonymized]/causal-wave-landings-repro/actions/workflows/reproduce.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository reproduces the ten cross-sector benchmark landings of
the causal-wave transport law from five measured coefficients with no
fits to downstream observables.

## Result in one line

```
10/10 benchmark landings within PRECISE_2.5-or-better.
Robust core (L1-L5) within 2.5% under fixed coefficients.
Rows L6-L8 (sin^2 theta_W, BH 1/4, Einstein gap 2/3) and L9-L10
(Lambda_t, Lambda_s) are EXACT but look-elsewhere-caveated and
are NOT exported as standalone proofs.
```

## Scope

This package presents one sharply scoped result. The causal-wave transport
equation
```
dC/dtau = D(Omega) * Laplacian(C) - gamma * C + alpha_xi * C
       + beta_pi * Pi_common * C + eps_sync^2 * C
```
with the five measured coefficients
```
alpha_xi  = 0.90082
D(Omega)  = 0.83996
beta_pi   = 0.93791
gamma     = 0.10021
eps_sync^2 = 0.05000
```
reproduces ten cross-sector benchmark observables as closed algebraic
compositions of these coefficients with small topological multipliers
(1/2, 1/3, 1/4, pi/4, 4/pi, N_gen=3). No additional free parameter is
introduced.

## What this is **not**

- Not a complete Standard-Model derivation
- Not a complete Quantum-Gravity theory
- Not a claim of structural EXACT-tier proof for L6/L7/L8 (these are
  enumeration-induced and look-elsewhere-caveated)
- Not a claim outside the canonical variation regime (P1)

## Installation (Windows PowerShell)

```powershell
git clone https://github.com/[anonymized]/causal-wave-landings-repro.git
cd causal-wave-landings-repro

py -3.11 -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Installation (POSIX)

```bash
git clone https://github.com/[anonymized]/causal-wave-landings-repro.git
cd causal-wave-landings-repro

python3.11 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Reproduce the result

```powershell
python .\src\recompute_landings.py
python .\src\run_look_elsewhere.py
python .\src\perturb_coefficients_null.py --seed 42 --trials 200
pytest
```

## Expected output (key lines)

```text
Five measured causal-wave coefficients loaded: PASS
10/10 benchmark landings within PRECISE-or-better: PASS
Rows L6-L8 (sin^2 theta_W, BH 1/4, Einstein gap 2/3) and L9-L10
(Lambda_t, Lambda_s) EXACT but look-elsewhere-caveated: PASS
Rows L1-L5 robust PRECISE_2.5 core within 2.5%: PASS
```

The complete expected output is in
[`outputs/expected_output.txt`](outputs/expected_output.txt).

## Repository structure

```
causal-wave-landings-repro/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── pyproject.toml
├── data/
│   ├── causal_wave_coefficients.json
│   ├── benchmark_targets_particle.json
│   ├── landings_expected.json
│   ├── selection_correction.json
│   └── coefficient_perturbation_null.json
├── src/
│   ├── recompute_landings.py
│   ├── run_look_elsewhere.py
│   ├── perturb_coefficients_null.py
│   └── make_figures.py
├── tests/
│   ├── test_coefficients.py
│   ├── test_landings.py
│   ├── test_look_elsewhere_flags.py
│   ├── test_null_models.py
│   └── test_falsification.py
├── outputs/
│   ├── expected_output.txt
│   ├── recompute_landings.json
│   ├── look_elsewhere_report.json
│   └── perturbation_null_run.json
├── paper/
│   ├── manuscript.tex
│   ├── manuscript.pdf
│   └── figures/
│       ├── fig1_schematic.{pdf,png}
│       ├── fig2_coefficients.{pdf,png}
│       ├── fig3_landings.{pdf,png}
│       └── fig4_lookelsewhere.{pdf,png}
└── .github/workflows/
    └── reproduce.yml
```

## Five coefficients and eight landings

| Coefficient        |     Value | Role                            |
|--------------------|----------:|---------------------------------|
| `alpha_xi`         |   0.90082 | Xi reaction rate                |
| `D(Omega)`         |   0.83996 | Effective diffusion             |
| `beta_pi`          |   0.93791 | Common-phase / holonomy mode    |
| `gamma`            |   0.10021 | Damping / time arrow            |
| `eps_sync^2`       |   0.05000 | Sync persistence layer          |

| #  | Observable                     | Formula                                              | Tier         | Look-elsewhere |
|----|--------------------------------|------------------------------------------------------|--------------|----------------|
| L1 | down-Yukawa exponent           | `alpha_xi + beta_pi + eps_sync^2 - gamma`            | PRECISE      | no             |
| L2 | `Omega_DM h^2`                 | `alpha_xi^2 * eps_sync^2 * N_gen`                    | PRECISE      | no             |
| L3 | `w_DE`                         | `-1 + eps_sync^4 / gamma`                            | PRECISE      | no             |
| L4 | `alpha_s(M_Z)`                 | `gamma * beta_pi * 4/pi`                             | PRECISE      | no             |
| L5 | charged-lepton Yukawa cluster  | `(alpha_xi/beta_pi) * (D(Omega)/(1+gamma)) * G_NET`  | PRECISE      | no             |
| L6 | `sin^2 theta_W`                | `beta_pi - (1-gamma) * pi/4`                         | EXACT        | **yes**        |
| L7 | Bekenstein-Hawking 1/4          | `alpha_xi/2 - 2*gamma`                              | EXACT        | no (Q-id.)     |
| L8 | Einstein gap 2/3               | `(1-gamma) * pi/4 - (1-D(Omega))/4`                  | EXACT        | **yes**        |
| L9 | `Lambda_t` (cosm. tensor 00)   | `alpha_xi^2 = (9/10)^2 = 81/100`                     | PRECISE      | no (Q-id.)     |
| L10| `Lambda_s` (cosm. tensor ii)   | `-gamma^2/2 = -eps_sync^2 * gamma = -1/200`          | EXACT        | no (Q-id.)     |

L9 and L10 are the time-time and spatial-isotropic components of
the anisotropic cosmological-constant tensor `Lambda_munu` in the
emergent Einstein equation, derived from a per-node 4x4 Frobenius
minimisation on the canonical lattice ladder (P4 §8.3b). Both
ascend to closed-form rational identities in Q under System R.

## Falsification

The closure mechanism fails if:

1. The five coefficients cannot be re-determined under the
   bounded-operator readout convergence criterion.
2. At least one of the five robust PRECISE rows L1-L5 deviates by more than
   2.5% from its target while the coefficients remain unchanged.
3. Random coefficient sets reproduce the same robust closure jointly
   (currently <10% under uniform [0.02, 0.99]^5 draws).

The look-elsewhere caveat for L6-L8 is itself NOT a falsification of the
robust core; it is a transparency statement.

## Citation

If you use this code or data, please cite:

```bibtex
@misc{bucciarelli2026causalwave,
  author    = {Bucciarelli, Sandro},
  title     = {A five-coefficient causal-wave ansatz and ten cross-sector benchmark closures},
  year      = {2026},
  version   = {0.1.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX}
}
```

See [CITATION.cff](CITATION.cff) for machine-readable metadata.

## Data integrity

SHA-256 hashes of all data files are in [`data/SHA256SUMS`](data/SHA256SUMS).
The published results correspond to data release `v0.1.0`.

## License

MIT License. See [LICENSE](LICENSE).
