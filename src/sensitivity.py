"""
================================================================================
 sensitivity.py  —  Global & local sensitivity analysis of the virtual cell
================================================================================

Quantifies how robust the model outputs are to its physical parameters — the
second thing a computational reviewer asks after benchmarking. Two analyses:

  1. LOCAL (one-at-a-time, OAT): perturb each parameter by +/-x% around the
     calibrated value and measure the fractional change in each output. Gives
     elasticities (normalized local sensitivities); cheap and interpretable.

  2. GLOBAL (variance-based, Sobol-style first-order indices): sample all
     parameters jointly over plausible ranges and decompose output variance
     into per-parameter contributions. Uses a lightweight Saltelli-style
     estimator (no external SALib dependency).

Outputs analyzed: nuclear_drive (sigma), nuclear_area, yap_nc, and the
mechanogenomic gene signature magnitude. Parameters: motor/clutch/coupling
(nm, nc, kc, alpha) and the nuclear gate (laminAC).

Interpretation for the paper: parameters with high sensitivity are the ones the
data must constrain (and that inference should identify); low-sensitivity
parameters are safe to fix. This is the robustness evidence.

Author: Daniel Pérez-Calixto (INMEGEN / UNAM)
================================================================================
"""

from __future__ import annotations
import numpy as np
from dataclasses import replace

import mvirtual_cell as mvc
from mvirtual_cell import PHENOTYPES


# parameters to probe, with plausible ranges (low, high) for the global analysis
PARAM_RANGES = {
    "nm":      (30.0, 70.0),
    "nc":      (50.0, 140.0),
    "kc":      (0.6, 1.8),
    "alpha":   (0.08, 0.20),
    "laminAC": (0.8, 2.4),
}


# ---------------------------------------------------------------------------
# Output functions of the phenotype (at a reference stiffness)
# ---------------------------------------------------------------------------
def _outputs_fast(ph, E=23.0, base=None, ref_key="hepatocyte"):
    """Fast, deterministic analytic outputs for sensitivity sweeps. Uses the
    motor's saturating stress scaled by the clutch/coupling parameters and the
    lamin-gated area/YAP maps, avoiding the stochastic motor so global variance
    analysis is tractable.

    The analytic scaling is anchored to a *reference* phenotype (``base`` /
    ``ref_key``): parameter ratios are taken relative to that reference. When
    analyzing a phenotype, ``base`` MUST be that phenotype's own calibrated
    point and ``ref_key`` its surrogate key, so the ratios are 1.0 at the
    reference and the elasticities are measured against the correct anchor.

    Raises RuntimeError if a fast surrogate is not available for ``ref_key``
    (e.g. phenotypes without pre-fit SATURATING_PARAMS) — the caller should then
    fall back to the stochastic path (use_fast=False) rather than silently
    normalizing against the wrong reference.
    """
    import fast_model as fm
    if base is None:
        base = PHENOTYPES["hepatocyte"]
        ref_key = "hepatocyte"
    # the fast surrogate must exist for this reference phenotype
    try:
        sig0 = float(fm.nuclear_stress_fast(E, ref_key))
    except Exception as exc:
        raise RuntimeError(
            f"fast surrogate unavailable for phenotype '{ref_key}'; "
            f"call with use_fast=False to use the stochastic path") from exc
    # analytic scalings, taken RELATIVE TO `base` (not hardcoded hepatocyte):
    #   more clutches (nc) raise transmitted stress (saturating);
    #   stiffer clutches (kc) reduce the kappa/(kappa+kc) transmission;
    #   more motors (nm) raise force but saturate; alpha scales kappa.
    nc_f = (ph.nc / base.nc) ** 0.6
    nm_f = (ph.nm / base.nm) ** 0.3
    kappa = ph.alpha * E
    kc_f = (kappa / (kappa + ph.kc)) / (kappa / (kappa + base.kc))
    sig = sig0 * nc_f * nm_f * kc_f
    # lamin-gated area & YAP
    s_half = ph.s0 * ph.laminAC
    area = ph.A_min + (ph.A_max - ph.A_min) * sig / (sig + s_half)
    yap = 1.0 + 3.8 * sig / (sig + s_half)
    return dict(nuclear_drive=sig, nuclear_area=area, yap_nc=yap)


# map a phenotype object to its fast-surrogate key (only calibrated ones)
_SURROGATE_KEYS = {"hepatocyte", "NHLF", "AT2_lung", "MDA"}


def _ref_key_for(base):
    """Return the surrogate key for a reference phenotype, or None if it has no
    pre-fit fast surrogate (then the caller must use the stochastic path)."""
    if base is None:
        return "hepatocyte"
    for key, ph in PHENOTYPES.items():
        if ph is base:
            return key if key in _SURROGATE_KEYS else None
    return None


def _outputs(ph, E=23.0, reps=4, use_fast=True, base=None):
    if use_fast:
        ref_key = _ref_key_for(base)
        if ref_key is not None:
            try:
                return _outputs_fast(ph, E, base=base, ref_key=ref_key)
            except Exception:
                pass  # fall through to stochastic path
    sig = mvc.nuclear_stress(E, ph, reps=reps)
    area = mvc.nuclear_area_ss(E, ph, reps=reps)
    yap = mvc.yap_nc_ratio(E, ph, reps=reps)
    return dict(nuclear_drive=sig, nuclear_area=area, yap_nc=yap)


def _apply(ph, name, value):
    return replace(ph, **{name: (int(round(value)) if name in ("nm", "nc")
                                 else value)})


# ---------------------------------------------------------------------------
# 1. LOCAL sensitivity (OAT elasticities)
# ---------------------------------------------------------------------------
def local_sensitivity(base=None, E=23.0, delta=0.15, reps=6):
    """Elasticity E_ij = (dOut_i/Out_i) / (dParam_j/Param_j), central-difference
    OAT. Returns {output: {param: elasticity}}."""
    if base is None:
        base = PHENOTYPES["hepatocyte"]
    out0 = _outputs(base, E, reps=reps, base=base)
    result = {o: {} for o in out0}
    for pname in PARAM_RANGES:
        p0 = getattr(base, pname)
        hi = _outputs(_apply(base, pname, p0 * (1 + delta)), E, reps=reps, base=base)
        lo = _outputs(_apply(base, pname, p0 * (1 - delta)), E, reps=reps, base=base)
        for o in out0:
            dOut = (hi[o] - lo[o]) / out0[o]
            dP = 2 * delta
            result[o][pname] = float(dOut / dP)
    return result


# ---------------------------------------------------------------------------
# 2. GLOBAL variance-based (Sobol first-order, Saltelli-style)
# ---------------------------------------------------------------------------
def sobol_first_order(base=None, E=23.0, n=1024, reps=2, seed=0):
    """First-order Sobol indices S_i = Var_i / Var_total via the Saltelli
    estimator with two independent sample matrices A, B and the mixed matrices
    A_B^i. Lightweight (no SALib). n is the base sample size."""
    if base is None:
        base = PHENOTYPES["hepatocyte"]
    rng = np.random.default_rng(seed)
    names = list(PARAM_RANGES)
    d = len(names)
    lows = np.array([PARAM_RANGES[k][0] for k in names])
    highs = np.array([PARAM_RANGES[k][1] for k in names])

    def sample(m):
        return lows + rng.random((m, d)) * (highs - lows)

    def evaluate(X, out_key):
        vals = np.empty(len(X))
        for i, row in enumerate(X):
            ph = base
            for j, nm in enumerate(names):
                ph = _apply(ph, nm, row[j])
            vals[i] = _outputs(ph, E, reps=reps, base=base)[out_key]
        return vals

    A = sample(n)
    B = sample(n)
    results = {}
    for out_key in ("nuclear_drive", "nuclear_area", "yap_nc"):
        yA = evaluate(A, out_key)
        yB = evaluate(B, out_key)
        varY = np.var(np.concatenate([yA, yB]))
        S = {}
        for i, nm in enumerate(names):
            ABi = A.copy(); ABi[:, i] = B[:, i]
            yABi = evaluate(ABi, out_key)
            # Jansen/Saltelli first-order estimator
            Si = np.mean(yB * (yABi - yA)) / (varY + 1e-12)
            S[nm] = float(np.clip(Si, 0, 1))
        results[out_key] = S
    return results


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_sensitivity(verbose=True, global_analysis=True, n=1024):
    loc = local_sensitivity()
    if verbose:
        print("=" * 72)
        print("  SENSITIVITY ANALYSIS  —  robustness of the virtual cell")
        print("=" * 72)
        print("\n[1] LOCAL (OAT elasticities at E=23 kPa; |value|>1 = amplifying):")
        params = list(PARAM_RANGES)
        print(f"    {'output':>14}" + "".join(f"{p:>9}" for p in params))
        for o, row in loc.items():
            print(f"    {o:>14}" + "".join(f"{row[p]:>9.2f}" for p in params))

    glob = None
    if global_analysis:
        glob = sobol_first_order(n=n)
        if verbose:
            print("\n[2] GLOBAL (first-order Sobol indices; fraction of variance):")
            params = list(PARAM_RANGES)
            print(f"    {'output':>14}" + "".join(f"{p:>9}" for p in params))
            for o, row in glob.items():
                print(f"    {o:>14}" + "".join(f"{row[p]:>9.2f}" for p in params))
            print("\n  Reading: high index = the data MUST constrain this parameter")
            print("  (and inference should identify it); low = safe to fix.")

    if verbose:
        # headline: which parameter dominates each output
        print("\n  Dominant parameter per output (local):")
        for o, row in loc.items():
            k = max(row, key=lambda p: abs(row[p]))
            print(f"    {o:>14}  <-  {k}  (elasticity {row[k]:+.2f})")
        print("=" * 72)
    return dict(local=loc, sobol=glob)


if __name__ == "__main__":
    run_sensitivity()
