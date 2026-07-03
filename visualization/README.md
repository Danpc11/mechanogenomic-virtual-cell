# Visualization

Model-generated visualizations of the mechanogenomic virtual cell — rendered
from the model's own state outputs, not hand-drawn cartoons.

## Variable → visual mapping

| Model output | Visual encoding |
|---|---|
| `E` (stiffness) | substrate bar thickness/color |
| `traction` T(E) | traction arrows at adhesions |
| `nuclear_drive` σ | nucleus fill warmth |
| `nuclear_area` A(t) | nucleus size |
| flattening (drive) | nucleus aspect ratio |
| `yap_nc` | dots inside the nucleus |
| `laminAC` | nuclear-envelope thickness |
| `nc_eff` | number of adhesions |
| gene scores | activation bars |

## Files

- `make_state_grid.py` — precompute states over E×t → `docs/states.json`
- `render_plotly_virtual_cell.py` — 3D Plotly HTML (level-1 embedded demo)
- `render_pyvista_virtual_cell.py` — PyVista+Trame app (level-2 scientific)

## Web demo

`docs/virtual_cell_demo.html` is a self-contained interactive demo (sliders for
stiffness and time) suitable for GitHub Pages. Enable Pages on the `docs/`
folder to publish it.
