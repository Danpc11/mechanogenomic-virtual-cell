# Code review notes

Findings from a review pass over `src/`, `test/`, and `results/make_figures.py`.
Nothing here affects the model's physics or calibrated numbers — all four
items are bookkeeping/consistency issues found while reading the code.

Three are fixed below (dead code removal + a duplication fix, all verified
against the test suite). One is left open because fixing it properly requires
a design decision only you can make.

---

## Fixed

### 1. Dead `HNF4A` entry in `gene_module.py`

**Was:** the top-level `GENES` dict defined `"HNF4A"` (linear, inverse), but
`CORE_GENES` never included it, so it was unreachable — every phenotype
actually scored the *other* `HNF4A` defined in `PHENOTYPE_GENES["hepatocyte"]`
(different `note`). `_INVERSE = {"HNF4A"}` referenced the dead key too.

**Fix:** removed the unreachable `GENES["HNF4A"]` entry and its comment;
`_INVERSE` is now an explicit empty set with a comment pointing at
`_PHENO_INVERSE` as the real mechanism for phenotype-specific inverse genes.
Verified `score_genes()` still returns `HNF4A` for hepatocyte (via the
phenotype-specific panel) and gene counts per phenotype are unchanged.

### 2. Dead variable in `virtual_cell.py::VirtualCell._drive()`

**Was:** `base = PHENOTYPES.get("hepatocyte")` was computed and never used.

**Fix:** removed the line. No behavior change.

### 3. Fibrosis-stage cutoffs were duplicated (and could have silently diverged)

**Was:** the same E→stage boundaries (6, 8, 10, 14 kPa) were hardcoded
independently in `virtual_cell.py::_stage_of_E()` (short codes like `"F0"`)
and `pharmacology.py::map_patient()` (descriptive labels like `"F0 (no/minimal
fibrosis)"`). They agreed today, but nothing enforced that — updating one
without the other would have made `VirtualCell.simulate()` and
`pharmacology.map_patient()` report different stages for the same stiffness.

**Fix:** added a single source of truth in `mvirtual_cell.py`:

```python
STAGE_CUTOFFS_KPA = [6.0, 8.0, 10.0, 14.0]   # upper bound of F0, F1, F2, F3
STAGE_LABELS = {"F0": "F0 (no/minimal fibrosis)", ...}

def stage_of_stiffness(E):
    ...  # returns the short code, e.g. "F0"
```

`virtual_cell.py` now calls `mvc.stage_of_stiffness(E)` directly (short code);
`pharmacology.py::map_patient()` calls `mvc.STAGE_LABELS[mvc.stage_of_stiffness(E)]`
for the descriptive label. Both call sites now derive from the same cutoffs.
Verified they agree at all five representative stiffnesses (4, 7, 9.5, 13, 26
kPa) and the full test suite still passes (11/11).

---

## Still open — needs your call

### 4. `results/make_figures.py::figure2()` hardcodes R² instead of computing it

**Where:** the line `r2 = [0.70, 0.90, 0.96, 0.982]` in `figure2()`.

This panel compares linear/power/log/saturating functional forms against the
motor-clutch simulator, but the R² values plotted are typed-in constants
rather than the output of `symbolic.compare_physical_forms()` — which already
exists and computes exactly this comparison from real `(E, sigma)` data. If
the calibration changes, this figure keeps showing the old numbers with no
warning.

**What's needed to resolve it:** a call on cost vs. reproducibility. Wiring in
the real values means generating `(E, sigma)` via
`symbolic.generate_stress_data()`, which runs the stochastic motor-clutch
simulator dozens of times (seconds, not instant) every time the figure script
runs. That's a one-off cost when regenerating figures, but it adds a new
runtime dependency from `make_figures.py` on `symbolic.py`, and changes the
script from "instant" to "takes a few seconds." Two ways to go:

- compute it for real (accurate, slower, one more import), or
- keep the constants but label them explicitly as illustrative in the plot
  (fast, but the figure stays disconnected from the actual model).

I didn't want to make that tradeoff for you, so this one's parked.

---

## Also noted, not touched: test-count / README mismatch

`README.md` says the suite has **17 validations** (lines ~270 and ~844), but
`test/test_virtual_cell.py::ALL_TESTS` has **11** actual `assert`-based tests.
The other 6 claimed validations (gene response-shape predictions,
phenotype-specific panels, benchmark generalization, sensitivity dominance,
bootstrap CI) only exist as printed output in `_demo()` functions across
`gene_module.py`, `benchmark.py`, `sensitivity.py`, `stats_ci.py` — none of
them fail CI if the underlying behavior breaks.

Two ways to close this, and they're not equivalent in effort:
- write ~6 new `assert`-based tests promoting those `_demo()` checks into
  `ALL_TESTS` (real regression protection, more work), or
- correct the count in the README to 11 (one line, but leaves those modules
  unprotected by CI).

Left open on purpose — didn't want to guess which of those you'd want.

---

## Also noted, not touched: `sensitivity.py` ignores its own `base` argument

**Where:** `src/sensitivity.py::_outputs_fast()`, the line
`base = PHENOTYPES["hepatocyte"]`.

`local_sensitivity(base=...)` and `sobol_first_order(base=...)` both accept an
arbitrary phenotype to analyze, but the fast analytic path always normalizes
against hepatocyte's `nc`/`nm`/`kc`, never against whatever `base` was passed
in. Calling `local_sensitivity(base=PHENOTYPES["MDA"])` today would silently
return elasticities computed against the wrong reference point — no
exception, just a quiet wrong number. It's invisible now because
`run_sensitivity()`, the only caller, never overrides `base`.

Fixing it well means picking one of:
- re-derive the scaling from `base` itself, which needs a fallback for
  phenotypes without pre-fit `SATURATING_PARAMS` (A549, MCF10A, fibroblast —
  e.g. calibrate on the fly via `fast_model.calibrate(base)`, which isn't
  free), or
- restrict the function's contract to hepatocyte only and raise explicitly
  for any other `base`, pushing other phenotypes onto the slower
  `use_fast=False` path.

Left open — this changes what sensitivity numbers mean for non-hepatocyte
phenotypes, which is your call.
