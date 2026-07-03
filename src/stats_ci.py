"""
================================================================================
 stats_ci.py  —  Confidence intervals & bootstrap for the virtual cell
================================================================================

Adds statistical rigor: nonparametric bootstrap confidence intervals for model
outputs and calibrated parameters, and simple resampling-based comparisons. A
computational reviewer expects uncertainty bands, not point estimates.

Provides:
  * bootstrap_ci        — percentile & BCa CIs for any statistic of a sample
  * bootstrap_model_output — CI on a model output propagated from data resampling
  * bootstrap_parameter — CI on a fitted parameter via case resampling
  * permutation_test    — resampling test for a difference between conditions
  * fold_change_ci      — CI on the stiffness fold-change (the primary result)

No external dependencies beyond numpy/scipy.

Author: Daniel Pérez-Calixto (INMEGEN / UNAM)
================================================================================
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy import stats as _sps


# ---------------------------------------------------------------------------
# Core bootstrap
# ---------------------------------------------------------------------------
def bootstrap_ci(sample, statistic=np.mean, n_boot=5000, alpha=0.05,
                 method="bca", seed=0):
    """Bootstrap confidence interval for `statistic` of a 1-D sample.

    method : "percentile" or "bca" (bias-corrected and accelerated).
    Returns dict(estimate, ci_low, ci_high, se, n_boot, method).
    """
    x = np.asarray(sample, float)
    n = len(x)
    rng = np.random.default_rng(seed)
    theta_hat = float(statistic(x))
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = np.array([statistic(x[i]) for i in idx])
    se = float(boot.std(ddof=1))

    if method == "percentile":
        lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    elif method == "bca":
        # bias correction
        z0 = _sps.norm.ppf((np.sum(boot < theta_hat) + 0.5) / (n_boot + 1))
        # acceleration via jackknife
        jack = np.array([statistic(np.delete(x, i)) for i in range(n)])
        jbar = jack.mean()
        num = np.sum((jbar - jack) ** 3)
        den = 6.0 * (np.sum((jbar - jack) ** 2) ** 1.5) + 1e-12
        a = num / den
        z = _sps.norm.ppf
        zl, zu = z(alpha / 2), z(1 - alpha / 2)
        a1 = _sps.norm.cdf(z0 + (z0 + zl) / (1 - a * (z0 + zl)))
        a2 = _sps.norm.cdf(z0 + (z0 + zu) / (1 - a * (z0 + zu)))
        lo, hi = np.percentile(boot, [100 * a1, 100 * a2])
    else:
        raise ValueError(method)
    return dict(estimate=theta_hat, ci_low=float(lo), ci_high=float(hi),
                se=se, n_boot=n_boot, method=method)


# ---------------------------------------------------------------------------
# Fold-change CI (the primary mechanosensitivity result)
# ---------------------------------------------------------------------------
def fold_change_ci(soft_values, stiff_values, n_boot=5000, alpha=0.05, seed=0):
    """CI on the stiffness fold-change (stiff/soft) via independent bootstrap of
    each condition. This is the headline mechanosensitivity statistic (~2.2x)."""
    s = np.asarray(soft_values, float)
    t = np.asarray(stiff_values, float)
    rng = np.random.default_rng(seed)
    fc = np.array([t[rng.integers(0, len(t), len(t))].mean()
                   / s[rng.integers(0, len(s), len(s))].mean()
                   for _ in range(n_boot)])
    lo, hi = np.percentile(fc, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return dict(fold_change=float(t.mean() / s.mean()),
                ci_low=float(lo), ci_high=float(hi), n_boot=n_boot)


# ---------------------------------------------------------------------------
# Parameter CI via case resampling
# ---------------------------------------------------------------------------
def bootstrap_parameter(fit_fn, data_rows, n_boot=1000, alpha=0.05, seed=0):
    """CI on a fitted parameter by resampling data cases.

    fit_fn      : callable(rows) -> scalar parameter estimate
    data_rows   : list/array of data cases (rows) to resample
    """
    rows = list(data_rows)
    n = len(rows)
    rng = np.random.default_rng(seed)
    theta_hat = float(fit_fn(rows))
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            boot.append(float(fit_fn([rows[i] for i in idx])))
        except Exception:
            continue
    boot = np.array(boot)
    lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return dict(estimate=theta_hat, ci_low=float(lo), ci_high=float(hi),
                se=float(boot.std(ddof=1)), n_boot=len(boot))


# ---------------------------------------------------------------------------
# Permutation test
# ---------------------------------------------------------------------------
def permutation_test(a, b, statistic=lambda x, y: np.mean(y) - np.mean(x),
                     n_perm=10000, seed=0):
    """Two-sample permutation test for a difference between conditions.
    Returns dict(observed, p_value, n_perm)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    obs = statistic(a, b)
    pool = np.concatenate([a, b])
    na = len(a)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        s = statistic(pool[:na], pool[na:])
        if abs(s) >= abs(obs):
            count += 1
    return dict(observed=float(obs), p_value=(count + 1) / (n_perm + 1),
                n_perm=n_perm)


# ---------------------------------------------------------------------------
# Demo on the real hepatocyte data
# ---------------------------------------------------------------------------
def _load_area_by_condition():
    try:
        from paths import DATA_DIR
        p = DATA_DIR / "hepatocyte_complete_data.json"
    except Exception:
        p = Path("hepatocyte_complete_data.json")
    d = json.load(open(p))
    soft = [a for a in d["complete_timecourse"]["1_kPa"]["pop_high"] if a]
    stiff = [a for a in d["complete_timecourse"]["23_kPa"]["pop_high"] if a]
    return soft, stiff


def _demo():
    print("=" * 70)
    print("  STATISTICS & BOOTSTRAP  —  uncertainty for the virtual cell")
    print("=" * 70)
    soft, stiff = _load_area_by_condition()

    print("\n[1] Bootstrap CI on mean nuclear area per stiffness (BCa):")
    for label, vals in [("1 kPa (soft)", soft), ("23 kPa (stiff)", stiff)]:
        r = bootstrap_ci(vals, np.mean, n_boot=4000)
        print(f"    {label:>16}: {r['estimate']:6.1f} um^2  "
              f"95% CI [{r['ci_low']:.1f}, {r['ci_high']:.1f}]  (SE {r['se']:.1f})")

    print("\n[2] Fold-change CI (the primary mechanosensitivity result):")
    fc = fold_change_ci(soft, stiff, n_boot=5000)
    print(f"    stiffness response (23/1 kPa) = {fc['fold_change']:.2f}x  "
          f"95% CI [{fc['ci_low']:.2f}, {fc['ci_high']:.2f}]")

    print("\n[3] Permutation test: is stiff area > soft area?")
    pt = permutation_test(soft, stiff)
    print(f"    observed difference = {pt['observed']:.1f} um^2, "
          f"p = {pt['p_value']:.4f}")
    print("=" * 70)


if __name__ == "__main__":
    _demo()
