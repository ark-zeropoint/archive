#!/usr/bin/env python3
"""Write Colab_TGV_Benchmark.ipynb next to this script.

Run:
  python make_colab_notebook.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "Colab_TGV_Benchmark.ipynb"


def md(source: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": source.strip("\n") + "\n",
    }


def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.strip("\n") + "\n",
    }


cells = [
    md(
        """
# NS_TZT Taylor-Green Vortex benchmark (Re = 1600)

Periodic box \(T^3=[0,2\\pi]^3\). Three cases on the same RK4 pseudospectral solver:

1. Unresolved baseline DNS (`nu = 1/1600`)
2. Standard Smagorinsky LES (`Cs = 0.17`)
3. NS_TZT selective dissipation (`C_O = 0.8`, `c_X = 0.1`, `Omega_c = 2`, `beta = 6`)

This notebook runs the **quick** grid (`N=32`, `t_end=2.0`) and reports the laminar-gate comparison at `t=0`.
"""
    ),
    md("## 1. Clone and install"),
    code(
        """
# Point REPO at the public archive that holds tzt_tgv_benchmark.py.
# Default assumes the solver sits next to this notebook after clone.
REPO = "https://github.com/ark-zeropoint/archive.git"
SUB = "Part2_Golden_M0_Tensor_Mechanics"

!rm -rf /content/archive
!git clone --depth 1 {REPO} /content/archive
%cd /content/archive/{SUB}
!pip -q install -r requirements.txt
!ls -la tzt_tgv_benchmark.py requirements.txt
"""
    ),
    md("## 2. Quick run (N=32, t-end 2.0)"),
    code(
        """
from pathlib import Path
OUT = Path("/content/tgv_quick")
OUT.mkdir(parents=True, exist_ok=True)
!python tzt_tgv_benchmark.py --n 32 --t-end 2.0 --out {OUT}
print("wrote", OUT)
"""
    ),
    md(
        """
## 3. Inline figures and the t=0 over-damping block

At t=0 the Taylor-Green field is still laminar, but `|S|` is already nonzero.
Smagorinsky therefore turns on everywhere. NS_TZT keeps `Phi=0` for `xi<0.8`
and only a small gate near the vorticity peaks, so most of that artificial drain is blocked.
"""
    ),
    code(
        """
from pathlib import Path
import csv
from IPython.display import Image, display, Markdown

OUT = Path("/content/tgv_quick")

for name in ("tgv_dissipation_rate.png", "tgv_energy_spectrum.png", "tgv_max_vorticity.png"):
    path = OUT / name
    display(Markdown(f"**{name}**"))
    display(Image(filename=str(path)))

rows = {}
with open(OUT / "tgv_metrics.csv", newline="") as f:
    for row in csv.DictReader(f):
        if abs(float(row["t"])) < 1e-12:
            rows[row["case"]] = row

eps_s = float(rows["smagorinsky"]["eps_model"])
eps_t = float(rows["tzt"]["eps_model"])
blocked = 100.0 * (1.0 - eps_t / eps_s) if eps_s else float("nan")

table = f'''
| case | E(0) | eps_nu(0) | eps_model(0) |
|---|---:|---:|---:|
| baseline | {float(rows["baseline"]["E"]):.6e} | {float(rows["baseline"]["eps_nu"]):.4e} | {float(rows["baseline"]["eps_model"]):.4e} |
| smagorinsky | {float(rows["smagorinsky"]["E"]):.6e} | {float(rows["smagorinsky"]["eps_nu"]):.4e} | {eps_s:.4e} |
| tzt | {float(rows["tzt"]["E"]):.6e} | {float(rows["tzt"]["eps_nu"]):.4e} | {eps_t:.4e} |

**t=0 Smagorinsky over-damping blocked by the Phi gate: {blocked:.1f}%**
'''
display(Markdown(table))
"""
    ),
]

nb = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        "colab": {"provenance": []},
    },
    "cells": cells,
}

# nbformat expects source as list of lines; a single string is accepted by most readers,
# but Colab is happier with line lists.
def as_lines(cell):
    src = cell["source"]
    if isinstance(src, str):
        if not src.endswith("\n"):
            src += "\n"
        cell["source"] = [line + "\n" for line in src.split("\n")[:-1]]
        if src.endswith("\n") and not cell["source"]:
            cell["source"] = ["\n"]
    return cell

nb["cells"] = [as_lines(c) for c in nb["cells"]]
OUT.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
print(f"wrote {OUT}")
