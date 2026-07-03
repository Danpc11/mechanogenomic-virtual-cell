"""
================================================================================
 make_state_grid.py  —  Precompute virtual-cell states for visualization
================================================================================

Runs the VirtualCell over a grid of stiffness x time and writes a JSON of full
states (mechanical + mechanogenomic) that the renderers and the web demo read.
Grid: E = 1, 5, 13, 23 kPa   x   t = 2, 36, 72, 120 h.

    python visualization/make_state_grid.py   ->   docs/states.json
================================================================================
"""
from __future__ import annotations
import sys
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from virtual_cell import VirtualCell
from paths import REPO_ROOT

E_GRID = [1, 5, 13, 23]          # kPa (1 & 23 have full experimental support)
T_GRID = [2, 36, 72, 120]        # h


def build(phenotype="hepatocyte"):
    cell = VirtualCell(phenotype)
    grid = {}
    for E in E_GRID:
        for t in T_GRID:
            s = cell.simulate(E, t=t)
            grid[f"E{E}_t{t}"] = dict(
                E=s.E_kPa, t=s.t_h, stage=s.fibrosis_stage,
                traction=round(s.traction, 1), drive=round(s.nuclear_drive, 1),
                area=round(s.nuclear_area, 1), yap=round(s.yap_nc, 2),
                lamin=round(s.laminAC, 2), nc=round(s.nc_eff, 0),
                tau=round(s.tau_h, 0), func=round(s.function_index, 2),
                genes={k: v for k, v in s.gene_scores.items()})
    return dict(phenotype=phenotype, E_grid=E_GRID, t_grid=T_GRID, states=grid)


if __name__ == "__main__":
    out = build()
    dest = REPO_ROOT / "docs" / "states.json"
    dest.parent.mkdir(exist_ok=True)
    json.dump(out, open(dest, "w"), indent=1)
    print(f"Wrote {len(out['states'])} states to {dest}")
