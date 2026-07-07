"""
================================================================================
 benchmark.py  —  Benchmarking the virtual cell against simple baselines
================================================================================

A reviewer's first question for any mechanistic model is: does the physics buy
you anything over a naive fit? This module answers that quantitatively by
comparing, on the SAME data, three models of increasing physical content:

  (a) LINEAR        area = a + b * E             (or log E)   — no dynamics
  (b) STATIC-DRIVE  area = A_min + (A_max-A_min)*sigma/(sigma+K)
                    uses the motor-clutch stress but NO temporal/nuclear layer
  (c) FULL          the mechanistic virtual cell: motor-clutch stress + lamin-
                    gated nucleus + stiffness-dependent relaxation tau(E)

Fairness measures:
  * same observable (mechanosensitive nuclear area over stiffness AND time),
  * leave-one-condition-out cross-validation (predict a held-out (E,t)),
  * AIC / BIC to penalize the extra parameters of the richer models,
  * report where each model fails (which conditions), not just a global number.

The point is NOT to show the full model always wins on raw R^2 (a flexible fit
can do that) but that it (i) generalizes to held-out conditions and (ii)
predicts the TEMPORAL structure the simple models structurally cannot.
================================================================================
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

import mvirtual_cell as mvc
from mvirtual_cell import PHENOTYPES
import fast_model as fm


# ---------------------------------------------------------------------------
# Data assembly: mechanosensitive nuclear area over (E, t)
# ---------------------------------------------------------------------------
def load_area_data():
    """Return arrays E, t, area for the mechanosensitive population.
    Uses the complete timecourse (1 & 23 kPa, 36-120 h) — the conditions with
    steady-state-relevant dynamics."""
    try:
        from paths import DATA_DIR
        p = DATA_DIR / "hepatocyte_complete_data.json"
    except Exception:
        p = Path("hepatocyte_complete_data.json")
    d = json.load(open(p))
    E, t, A = [], [], []
    for key, Eval in [("1_kPa", 1.0), ("23_kPa", 23.0)]:
        c = d["complete_timecourse"][key]
        for th, a in zip(c["t_h"], c["pop_high"]):
            if a is not None:
                E.append(Eval); t.append(float(th)); A.append(float(a))
    return np.array(E), np.array(t), np.array(A)


# ---------------------------------------------------------------------------
# The three models (each: fit params on train, predict on any (E, t))
# ---------------------------------------------------------------------------
class LinearModel:
    """(a) area = a + b * log10(E). No time dependence. 2 params."""
    name = "linear (E->area)"
    n_params = 2

    def fit(self, E, t, A):
        x = np.log10(E)
        self.b, self.a = np.polyfit(x, A, 1)
        return self

    def predict(self, E, t):
        return self.a + self.b * np.log10(np.atleast_1d(E))


class StaticDriveModel:
    """(b) motor-clutch stress -> area, but STATIC (no tau, no time).
    area = A_min + (A_max-A_min)*sigma(E)/(sigma+K). 3 params."""
    name = "static-drive (no nucleus dynamics)"
    n_params = 3

    def fit(self, E, t, A):
        sig = fm.nuclear_stress_fast(E, "hepatocyte")

        def f(sig, Amin, Amax, K):
            return Amin + (Amax - Amin) * sig / (sig + K)
        p0 = [60, 250, 20]
        self.p, _ = curve_fit(f, sig, A, p0=p0,
                              bounds=([30, 100, 1], [90, 400, 200]), maxfev=20000)
        return self

    def predict(self, E, t):
        sig = fm.nuclear_stress_fast(np.atleast_1d(E), "hepatocyte")
        Amin, Amax, K = self.p
        return Amin + (Amax - Amin) * sig / (sig + K)


class FullModel:
    """(c) full mechanistic virtual cell: motor-clutch stress + lamin-gated
    nucleus + stiffness-dependent relaxation tau(E). Fits A_min, A_max, lamin,
    A0 (4 params); tau(E) is set by the calibrated anchors (not refit here)."""
    name = "full mechanistic (nucleus + tau(E))"
    n_params = 4

    def __init__(self):
        self.ph = PHENOTYPES["hepatocyte"]

    def fit(self, E, t, A):
        ph = self.ph

        def f(X, Amin, Amax, lam, A0):
            Ev, tv = X
            sig = fm.nuclear_stress_fast(Ev, "hepatocyte")
            Ass = Amin + (Amax - Amin) * sig / (sig + ph.s0 * lam)
            tau = np.array([mvc.tau_of_E(e, ph) for e in np.atleast_1d(Ev)])
            return Ass + (A0 - Ass) * np.exp(-tv / tau)
        p0 = [55, 260, 1.2, 65]
        self.p, _ = curve_fit(f, (E, t), A, p0=p0,
                              bounds=([30, 150, 0.5, 40], [80, 400, 3.0, 90]),
                              maxfev=30000)
        return self

    def predict(self, E, t):
        ph = self.ph
        Amin, Amax, lam, A0 = self.p
        E = np.atleast_1d(E); t = np.atleast_1d(t)
        sig = fm.nuclear_stress_fast(E, "hepatocyte")
        Ass = Amin + (Amax - Amin) * sig / (sig + ph.s0 * lam)
        tau = np.array([mvc.tau_of_E(e, ph) for e in E])
        return Ass + (A0 - Ass) * np.exp(-t / tau)


MODELS = [LinearModel, StaticDriveModel, FullModel]


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def _r2(y, yhat):
    y = np.asarray(y); yhat = np.asarray(yhat)
    ss = np.sum((y - yhat) ** 2)
    return 1 - ss / np.sum((y - y.mean()) ** 2)


def _aic(y, yhat, k):
    n = len(y)
    rss = np.sum((np.asarray(y) - np.asarray(yhat)) ** 2)
    return n * np.log(rss / n) + 2 * k


def _bic(y, yhat, k):
    n = len(y)
    rss = np.sum((np.asarray(y) - np.asarray(yhat)) ** 2)
    return n * np.log(rss / n) + k * np.log(n)


def fit_and_score(model_cls, E, t, A):
    """In-sample fit + AIC/BIC (penalizes complexity)."""
    m = model_cls().fit(E, t, A)
    pred = np.array([m.predict(e, tt)[0] for e, tt in zip(E, t)])
    return dict(name=model_cls.name, k=model_cls.n_params,
                R2=_r2(A, pred), AIC=_aic(A, pred, model_cls.n_params),
                BIC=_bic(A, pred, model_cls.n_params), pred=pred)


def cross_validate(model_cls, E, t, A):
    """Leave-one-condition-out CV: fit on n-1 points, predict the held-out one.
    This is the honest generalization test (a flexible fit can overfit in-sample)."""
    n = len(A)
    preds = np.zeros(n)
    for i in range(n):
        mask = np.arange(n) != i
        m = model_cls().fit(E[mask], t[mask], A[mask])
        preds[i] = m.predict(E[i], t[i])[0]
    return dict(name=model_cls.name, cv_R2=_r2(A, preds),
                cv_RMSE=float(np.sqrt(np.mean((A - preds) ** 2))), cv_pred=preds)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_benchmark(verbose=True):
    E, t, A = load_area_data()
    insample = [fit_and_score(m, E, t, A) for m in MODELS]
    cv = [cross_validate(m, E, t, A) for m in MODELS]

    if verbose:
        print("=" * 74)
        print("  BENCHMARK  —  full mechanistic model vs simple baselines")
        print("=" * 74)
        print(f"  Observable: mechanosensitive nuclear area over (E, t), "
              f"n={len(A)} conditions\n")
        print(f"  {'model':>36}{'k':>3}{'R2':>7}{'AIC':>8}{'BIC':>8}{'CV-R2':>8}")
        for ins, c in zip(insample, cv):
            print(f"  {ins['name']:>36}{ins['k']:>3}{ins['R2']:>7.2f}"
                  f"{ins['AIC']:>8.1f}{ins['BIC']:>8.1f}{c['cv_R2']:>8.2f}")
        print("\n  Lower AIC/BIC = better (complexity-penalized). CV-R2 = "
              "held-out generalization.")

        # where the simple models fail: the temporal structure
        print("\n  Temporal test — predicted area at 23 kPa across time:")
        print(f"    {'t(h)':>6}{'obs':>7}", end="")
        for m in MODELS:
            print(f"{m.name.split()[0]:>10}", end="")
        print()
        fits = [m().fit(E, t, A) for m in MODELS]
        for th in [36, 72, 120]:
            obs = A[(E == 23) & (t == th)]
            row = f"    {th:>6}{obs[0] if len(obs) else float('nan'):>7.0f}"
            for fit in fits:
                row += f"{fit.predict(23.0, th)[0]:>10.0f}"
            print(row)
        print("\n  Note: the linear and static-drive models predict the SAME area")
        print("  at every time (no temporal term); only the full model captures")
        print("  the slow rise on stiff substrate (the tau(E) dynamical law).")
        print("=" * 74)

    return dict(insample=insample, cv=cv, E=E, t=t, A=A)


if __name__ == "__main__":
    run_benchmark()
