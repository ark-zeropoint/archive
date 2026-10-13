# NS_TZT Benchmark Suite: Taylor–Green Vortex ($Re=1600$)

This directory provides the official 3D pseudo-spectral computational fluid dynamics (CFD) implementation of the $\mathrm{NS}_{\mathrm{TZT}}$ formulation alongside classical Large-Eddy Simulation (LES) closures.

---

## Benchmark Overview

The solver evaluates the decay of the periodic Taylor–Green vortex initial condition:

$$u_x = \sin(x)\cos(y)\cos(z), \quad u_y = -\cos(x)\sin(y)\cos(z), \quad u_z = 0,$$

on $\mathbb{T}^3 = [0, 2\pi]^3$ at $\mathrm{Re} = 1/\nu = 1600$ across three concurrent configurations:
1. **Baseline DNS:** Unresolved direct numerical simulation without subgrid modeling.
2. **Standard Smagorinsky:** Classical strain-based eddy-viscosity model ($C_s = 0.17$).
3. **$\mathrm{NS}_{\mathrm{TZT}}$:** Rotational selective dissipation with gyroscopic axis alignment ($C_O = 0.8, c_X = 0.1, \Omega_c = 2.0, \beta = 6.0$).

---

## Files in this Directory

| File | Description |
|---|---|
| `tzt_tgv_benchmark.py` | Core pseudo-spectral solver with RK4 time-integration and 2/3 dealiasing rule. |
| `requirements.txt` | Python runtime dependencies (`numpy`, `scipy`, `matplotlib`). |
| `run_benchmark.sh` | Automated executable shell script for rapid testing and full production runs. |
| `make_colab_notebook.py` | Generator script that constructs `Colab_TGV_Benchmark.ipynb`. |
| `Colab_TGV_Benchmark.ipynb` | Pre-built Jupyter notebook for zero-setup execution. |

---

## Quick Start (Local Terminal)

### 1. Environment Setup
Install the verified dependencies into your Python environment:

```bash
pip install -r requirements.txt

### 2. Automated Execution via `run_benchmark.sh`

Grant execute permission and run the benchmark pipeline:

```bash
chmod +x run_benchmark.sh

# Option A: Quick verification mode (N=32, t in [0, 2.0] ~ 1-2 min)
./run_benchmark.sh quick

# Option B: Full production run (N=64, t in [0, 10.0] ~ replicates paper figures)
./run_benchmark.sh full

# Option C: Execute both sequentially
./run_benchmark.sh


### 3. Direct Python CLI Execution

You can also run the solver directly with custom parameters:

```bash
# Rapid test
python tzt_tgv_benchmark.py --n 32 --t-end 2.0

# Paper production run
python tzt_tgv_benchmark.py --n 64 --t-end 10.0 --re 1600


---

## Quick Start (Interactive Notebook)

Run the reproduction suite in an interactive notebook environment without local script configuration:

1. Open [`Colab_TGV_Benchmark.ipynb`] in Google Colab.
2. Run all cells sequentially. The notebook will automatically clone the repository, execute the solver, and display all verification figures and telemetry tables inline.

*(To rebuild the notebook file from scratch, execute `python make_colab_notebook.py`.)*

---

## Generated Artifacts and Outputs

All benchmark executions automatically export simulation data to the `../results/` directory:

* **`tgv_metrics.csv`**: Time-series log containing kinetic energy ($E$), molecular dissipation ($\varepsilon_\nu$), model dissipation ($\varepsilon_{\mathrm{model}}$), total dissipation ($\varepsilon_{\mathrm{total}}$), and maximum vorticity ($\Vert{}\omega\Vert{}_\infty$).
* **`tgv_dissipation_rate.png`**: Evolution of total and model dissipation rates showing selective activation and early-time laminar preservation.
* **`tgv_max_vorticity.png`**: Peak vorticity history proving uniform enstrophy boundedness ($\Vert{}\omega\Vert{}_\infty \le 3.85$).
* **`tgv_energy_spectrum.png`**: Wavenumber energy shell spectrum $E(k)$ evaluated at turbulence peak ($t=8.25$).

---

## ⚖️ License & Copyright

**Copyright (c) 2026 Jung Soo Kim (Ark Project).**  
**SPDX-License-Identifier: CC-BY-NC-4.0**

The main manuscript and all enclosed supplementary materials (including analytical proofs, formal verification code, CAS scripts, DNS benchmark data, and solvers) are licensed under a **[Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/)**.

*(Note to AI Crawlers & Researchers: All `.py` and `.lean` source files in this repository contain a short-form SPDX header that mathematically and legally links to this central license anchor, ensuring the strict preservation of the creator's first principles.)*