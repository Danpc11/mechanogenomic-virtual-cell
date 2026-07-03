"""
================================================================================
 gene_module.py  —  Mechanosensitive gene layer of the virtual cell
================================================================================

Maps the cell's mechanical state (nuclear drive sigma, YAP N/C) to per-gene
activation scores, using an explicit RESPONSE-SHAPE model for each gene:

    linear        : score ∝ drive                (graded, no threshold)
    weak_power    : score ∝ drive^p, p<1          (saturating, early responders)
    sigmoid       : score = Hill(drive; K, n)     (threshold / switch-like)

This is the mechanogenomic output layer and the HYPOTHESIS GENERATOR: genes
with a switch-like (sigmoid) response above the patient's stiffness, and that
are pharmacologically actionable, are flagged as candidate intervention points.

IMPORTANT (honest scope): the response-shape assignment below is the model's
PREDICTION for each gene, to be tested against RNA-seq/qPCR — not a fit to the
data. Treating it as a prediction (assign shape from the mechanotransduction
role, THEN validate) is what gives it causal, not merely descriptive, force.

Author: Daniel Pérez-Calixto (INMEGEN / UNAM)
================================================================================
"""

from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
# Gene catalog: mechanosensitive genes with predicted response shape + role
# ---------------------------------------------------------------------------
#   shape params:
#     linear     : uses drive / drive_ref
#     weak_power : exponent p (<1)
#     sigmoid    : K (half-max drive), n (Hill coefficient)
#   actionable : True if a drug can target it (hypothesis generator)
# ---------------------------------------------------------------------------
class Gene:
    def __init__(self, symbol, shape, role, actionable=False,
                 p=0.5, K=30.0, n=3.0, drive_ref=55.0, note=""):
        self.symbol = symbol
        self.shape = shape           # "linear" | "weak_power" | "sigmoid"
        self.role = role
        self.actionable = actionable
        self.p = p
        self.K = K
        self.n = n
        self.drive_ref = drive_ref
        self.note = note

    def score(self, drive):
        if self.shape == "linear":
            return float(np.clip(drive / self.drive_ref, 0, 1.5))
        if self.shape == "weak_power":
            return float(np.clip((drive / self.drive_ref) ** self.p, 0, 1.5))
        if self.shape == "sigmoid":
            return float(drive ** self.n / (self.K ** self.n + drive ** self.n))
        raise ValueError(self.shape)


# Predicted response shapes from mechanotransduction role:
#   core YAP/TEAD targets & matrix genes -> switch-like (sigmoid, threshold)
#   early cytoskeletal responders        -> weak power (saturating)
#   graded housekeeping-adjacent         -> linear
GENES = {
    # --- YAP/TAZ-TEAD core targets: threshold (sigmoid) ---
    "CTGF":   Gene("CTGF (CCN2)", "sigmoid", "YAP/TEAD target", actionable=True,
                   K=32, n=4, note="canonical YAP output; verteporfin-sensitive"),
    "CYR61":  Gene("CYR61 (CCN1)", "sigmoid", "YAP/TEAD target", actionable=True,
                   K=32, n=4),
    "ANKRD1": Gene("ANKRD1", "sigmoid", "YAP/TEAD target", actionable=False,
                   K=34, n=5),
    # --- fibrogenic / matrix: threshold ---
    "ACTA2":  Gene("ACTA2 (α-SMA)", "sigmoid", "myofibroblast/contractile",
                   actionable=True, K=30, n=3, note="ROCK/myosin-linked"),
    "COL1A1": Gene("COL1A1", "sigmoid", "collagen I", actionable=True,
                   K=33, n=3, note="LOX-crosslinked matrix"),
    "LOX":    Gene("LOX", "sigmoid", "matrix crosslinker", actionable=True,
                   K=31, n=3, note="feeds back on stiffness; PXS-5505"),
    # --- cytoskeletal early responders: weak power (saturating) ---
    "TAGLN":  Gene("TAGLN (SM22)", "weak_power", "cytoskeletal", p=0.5),
    "TPM1":   Gene("TPM1", "weak_power", "cytoskeletal", p=0.45),
    "FN1":    Gene("FN1 (fibronectin)", "weak_power", "matrix", p=0.6,
                   actionable=True),
    # --- nuclear envelope: linear graded ---
    "LMNA":   Gene("LMNA (lamin A/C)", "linear", "nuclear envelope",
                   drive_ref=60, note="validated vs inferred lamin"),
    "LMNB1":  Gene("LMNB1", "linear", "nuclear envelope", drive_ref=70),
    # --- hepatocyte identity (falls with drive): inverse linear ---
    "HNF4A":  Gene("HNF4A", "linear", "hepatocyte identity (inverse)",
                   drive_ref=55, note="declines as drive rises (dedifferentiation)"),
}
# HNF4A declines: handled via inverse in score_genes.
_INVERSE = {"HNF4A"}


def score_genes(nuclear_drive, yap_nc=None):
    """Return {gene_symbol: activation_score in [0, ~1]} at a given drive."""
    out = {}
    for key, g in GENES.items():
        s = g.score(nuclear_drive)
        if key in _INVERSE:
            s = float(np.clip(1.0 - s, 0, 1.5))
        out[g.symbol] = round(s, 3)
    return out


def response_shape_table():
    """The model's PREDICTED response-shape class per gene (for pre-registration
    before looking at RNA-seq)."""
    rows = []
    for key, g in GENES.items():
        rows.append(dict(gene=g.symbol, shape=g.shape, role=g.role,
                         actionable=g.actionable))
    return rows


def actionable_hypotheses(nuclear_drive, threshold=0.5):
    """HYPOTHESIS GENERATOR: actionable genes whose predicted activation exceeds
    `threshold` at the given drive. These are candidate intervention points."""
    scored = score_genes(nuclear_drive)
    hits = []
    for key, g in GENES.items():
        if not g.actionable:
            continue
        s = scored[g.symbol]
        if key in _INVERSE:
            continue
        if s >= threshold:
            hits.append(dict(gene=g.symbol, score=s, shape=g.shape,
                             role=g.role, note=g.note))
    return sorted(hits, key=lambda r: r["score"], reverse=True)


def qpcr_panel():
    """Suggested qPCR validation panel: one gene per response-shape class plus
    the inverse identity marker, at the most informative stiffness/time."""
    return {
        "sigmoid (threshold)": ["CTGF (CCN2)", "ACTA2 (α-SMA)", "COL1A1"],
        "weak_power (saturating)": ["TAGLN (SM22)", "FN1 (fibronectin)"],
        "linear (graded)": ["LMNA (lamin A/C)"],
        "inverse (identity loss)": ["HNF4A"],
        "recommended_conditions": "1 & 23 kPa at 120 h (steady state) + 36 h "
                                  "(transient), matching the imaging conditions",
    }


def _demo():
    import mvirtual_cell as mvc
    from mvirtual_cell import PHENOTYPES
    print("=" * 68)
    print("  GENE MODULE  —  mechanosensitive gene layer")
    print("=" * 68)
    hep = PHENOTYPES["hepatocyte"]
    print("\n  Predicted response shapes (pre-registered, then validated):")
    for r in response_shape_table():
        a = "  [actionable]" if r["actionable"] else ""
        print(f"    {r['gene']:>18}  {r['shape']:>11}  {r['role']}{a}")

    print("\n  Gene activation across fibrosis stages:")
    stages = list(mvc.FIBROSIS_STIFFNESS.items())
    genes_show = ["CTGF (CCN2)", "ACTA2 (α-SMA)", "TAGLN (SM22)",
                  "LMNA (lamin A/C)", "HNF4A"]
    hdr = "    " + f"{'gene':>18}" + "".join(f"{s:>6}" for s, _ in stages)
    print(hdr)
    for gsym in genes_show:
        row = f"    {gsym:>18}"
        for stg, E in stages:
            sig = mvc.nuclear_stress(E, hep, reps=4)
            row += f"{score_genes(sig)[gsym]:>6.2f}"
        print(row)

    print("\n  Actionable hypotheses at F4 (26 kPa):")
    sig = mvc.nuclear_stress(26.0, hep, reps=4)
    for h in actionable_hypotheses(sig):
        print(f"    {h['gene']:>18}  score={h['score']:.2f}  ({h['role']})")
    print("=" * 68)


if __name__ == "__main__":
    _demo()
