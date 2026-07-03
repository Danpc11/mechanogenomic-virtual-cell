"""
================================================================================
 paths.py  —  Centralized path resolution for the repository
================================================================================

Resolves the repository directories from this file's location, so every module
finds data/ and results/ regardless of the current working directory. Import
the constants below instead of hard-coding relative paths.

    from paths import DATA_DIR, RESULTS_DIR, FIGURES_DIR
    df = pd.read_csv(DATA_DIR / "hepatocyte_two_populations.csv")
================================================================================
"""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent   # repository root
SRC_DIR = REPO_ROOT / "src"
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
