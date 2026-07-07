"""
================================================================================
 hepatic_panel.py — extended hepatic mechanotransduction gene panel
================================================================================

Extends the core mechanosensitive panel (in gene_module.py) with the broader set
of mechanotransduction genes measured in the liver RNA-seq validation cohorts,
organized so predictions stay falsifiable.

Design principle (anti-circularity):
  The response SHAPE for each gene is NOT hand-assigned gene-by-gene. Each gene
  is placed in a mechanistic CATEGORY, and a fixed RULE maps category -> shape:

      YAP/TEAD effectors, contractile activation  -> sigmoid    (threshold switch)
      adhesion / mechanosensor / cytoskeleton     -> weak_power (saturating)
      nuclear envelope-LINC, epigenetic, graded TF-> linear     (dose-graded)

  DIRECTION (up / down with stiffness) is assigned from each gene's role:
  identity/function markers and stiffness-repressed regulators fall (down);
  effectors and structural reinforcers rise (up).

  Both assignments are made from prior biology, independent of the observed
  fibrosis-stage data, so the transcriptomic comparison remains a genuine test.

This module feeds the same Gene objects used by gene_module, so the validation
script and scoring work unchanged.
================================================================================
"""

from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gene_module import Gene


# ---------------------------------------------------------------------------
# The RULE: functional category -> predicted response shape
# ---------------------------------------------------------------------------
CATEGORY_SHAPE = {
    "yap_tead_effector":   "sigmoid",     # threshold switch on nuclear YAP/TAZ
    "fibrotic_effector":   "sigmoid",     # de-novo myofibroblast/matrix program
    "mechanosensor":       "weak_power",  # adhesion sensing, saturating
    "adhesion":            "weak_power",
    "cytoskeleton":        "weak_power",  # contractile machinery, saturating
    "nuclear_envelope":    "linear",      # LINC / lamina, stiffness-graded
    "epigenetic":          "linear",      # chromatin regulators, graded
    "graded_tf":           "linear",      # graded transcription factors
    "identity":            "linear",      # identity marker (direction = down)
}


def _g(symbol, category, direction, note="", actionable=False):
    """Build a Gene with shape set by the category RULE (not by hand)."""
    shape = CATEGORY_SHAPE[category]
    return dict(gene=Gene(symbol, shape, category.replace("_", " "),
                          actionable=actionable, note=note),
                category=category, direction=direction)


# ---------------------------------------------------------------------------
# Extended hepatic panel. Each entry: symbol -> (category, direction).
# Shapes are derived from the category rule above; nothing is hand-tuned.
# ---------------------------------------------------------------------------
HEPATIC_PANEL_SPEC = {
    # --- YAP/TAZ/TEAD effectors (sigmoid, up) ---
    "YAP1":   ("yap_tead_effector", "up"),
    "WWTR1":  ("yap_tead_effector", "up"),   # TAZ
    "TEAD2":  ("yap_tead_effector", "up"),
    "TEAD4":  ("yap_tead_effector", "up"),
    "CCN2":   ("yap_tead_effector", "up"),   # CTGF, canonical target
    "MKL1":   ("yap_tead_effector", "up"),   # MRTF-A
    "SMAD2":  ("yap_tead_effector", "up"),   # TGF-beta/YAP crosstalk
    # --- fibrotic effectors (sigmoid, up) ---
    "ACTA2":  ("fibrotic_effector", "up"),
    "COL1A1": ("fibrotic_effector", "up"),
    "COL1A2": ("fibrotic_effector", "up"),
    "LOX":    ("fibrotic_effector", "up"),
    # --- mechanosensors / adhesion (weak_power) ---
    "PIEZO1": ("mechanosensor", "up"),
    "PTK2":   ("mechanosensor", "up"),       # FAK
    "ILK":    ("adhesion", "up"),
    "VCL":    ("adhesion", "up"),
    "SRC":    ("adhesion", "up"),
    "FERMT2": ("adhesion", "down"),          # kindlin-2, context-repressed
    # --- contractile cytoskeleton (weak_power, up) ---
    "MYH9":   ("cytoskeleton", "up"),
    "MYH10":  ("cytoskeleton", "up"),
    "MYL9":   ("cytoskeleton", "up"),
    "ACTB":   ("cytoskeleton", "up"),
    "CFL1":   ("cytoskeleton", "up"),
    "FLNA":   ("cytoskeleton", "up"),
    # --- nuclear envelope / LINC (linear) ---
    "LMNA":   ("nuclear_envelope", "up"),
    "LMNB2":  ("nuclear_envelope", "up"),
    "SYNE2":  ("nuclear_envelope", "down"),  # nesprin-2
    "SYNE3":  ("nuclear_envelope", "down"),  # nesprin-3
    "NUP93":  ("nuclear_envelope", "up"),
    "TPR":    ("nuclear_envelope", "up"),
    "TMPO":   ("nuclear_envelope", "up"),    # LAP2
    # --- epigenetic regulators (linear) ---
    "HDAC1":  ("epigenetic", "up"),
    "HDAC2":  ("epigenetic", "up"),
    "HDAC3":  ("epigenetic", "down"),
    "DNMT3A": ("epigenetic", "up"),
    "EHMT1":  ("epigenetic", "down"),
    "KDM4A":  ("epigenetic", "down"),
    "KDM4B":  ("epigenetic", "down"),
    "SUV39H1":("epigenetic", "down"),
    "SUZ12":  ("epigenetic", "down"),
    "CREBBP": ("epigenetic", "down"),
    # --- graded TF (linear, up) ---
    "KLF2":   ("graded_tf", "up"),
    # --- identity marker (linear, down) ---
    "HNF4A":  ("identity", "down"),
}


def hepatic_panel():
    """Return {symbol: Gene} for the extended hepatic panel, shapes set by rule."""
    panel = {}
    for sym, (cat, _dir) in HEPATIC_PANEL_SPEC.items():
        panel[sym] = Gene(sym, CATEGORY_SHAPE[cat], cat.replace("_", " "))
    return panel


def hepatic_directions():
    """Return {symbol: 'up'|'down'} predicted direction for each panel gene."""
    return {sym: d for sym, (_c, d) in HEPATIC_PANEL_SPEC.items()}


def hepatic_categories():
    """Return {symbol: category} for reporting/grouping."""
    return {sym: c for sym, (c, _d) in HEPATIC_PANEL_SPEC.items()}


if __name__ == "__main__":
    panel = hepatic_panel()
    dirs = hepatic_directions()
    cats = hepatic_categories()
    print(f"Extended hepatic panel: {len(panel)} genes\n")
    print(f"  {'gene':10s}{'category':22s}{'shape':12s}{'direction'}")
    for sym in panel:
        print(f"  {sym:10s}{cats[sym]:22s}{panel[sym].shape:12s}{dirs[sym]}")
