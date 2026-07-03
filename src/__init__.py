"""
mechanogenomic-virtual-cell — a phenotype-aware physical-computational virtual
cell linking tissue stiffness to nuclear form and mechanogenomic state.

Dual-use import model:
  * In-place (cloned repo): add `src/` to sys.path and import modules by name,
    e.g. `import mvirtual_cell` (the modules import each other this way).
  * Installed package (`pip install .`): import under the `mvcell` namespace,
    e.g. `from mvcell.virtual_cell import VirtualCell`. To keep the internal
    bare-name imports working, this __init__ adds the package directory to
    sys.path on import.
"""
import os
import sys

__version__ = "0.2.0"

# make internal bare-name imports (e.g. `import mvirtual_cell`) resolve whether
# used in-place or as the installed `mvcell` package
_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)

try:
    from virtual_cell import VirtualCell, CellState            # noqa: F401
    from mvirtual_cell import PHENOTYPES                        # noqa: F401
except Exception:
    pass
