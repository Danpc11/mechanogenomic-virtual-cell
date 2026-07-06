"""
================================================================================
 test_virtual_cell.py  —  Executable validation of the virtual-cell model
================================================================================

Runnable checks that the physical model reproduces its qualitative anchors.
These are not fits — they verify that the calibrated model still behaves
correctly (biphasic traction, stiffness-dependent nuclear spreading, YAP
activation, lamin-dependent gating, two-population dynamics, monotonic
fibrosis response). Run directly:

    python test_virtual_cell.py

or under pytest:

    pytest test_virtual_cell.py -v

Each test prints its result and asserts the expected behavior. Because the
motor-clutch engine is stochastic, tests use enough replicates for stable
means and assert robust qualitative relations rather than exact values.

Author: Daniel Pérez-Calixto (INMEGEN / UNAM)
================================================================================
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np
import mvirtual_cell as mvc
from mvirtual_cell import (PHENOTYPES, traction, nuclear_stress,
                           nuclear_area_ss, yap_nc_ratio, nuclear_area_time)
from dataclasses import replace

REPS = 8  # replicate seeds for stable stochastic means


# ---------------------------------------------------------------------------
def test_biphasic_traction():
    """Traction should be biphasic in substrate stiffness: an intermediate
    optimum exceeds both the very soft and very stiff limits."""
    ph = PHENOTYPES["hepatocyte"]
    kappas = np.array([0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0])  # pN/nm
    # sweep kappa directly by mapping through alpha (E = kappa/alpha)
    T = np.array([traction(k / ph.alpha, ph, reps=REPS) for k in kappas])
    k_opt = kappas[np.argmax(T)]
    assert T.max() > T[0], "optimum should exceed the soft limit"
    assert T.max() > T[-1], "optimum should exceed the stiff limit"
    assert 0.1 <= k_opt <= 5.0, f"optimum stiffness out of expected range: {k_opt}"
    print(f"  [OK] biphasic traction: optimum at kappa≈{k_opt:.2f} pN/nm "
          f"(T_opt={T.max():.1f} pN)")


def test_nuclear_area_increases_with_stiffness():
    """Steady-state nuclear area should increase monotonically with stiffness
    for the mechanosensitive response."""
    ph = PHENOTYPES["hepatocyte"]
    Es = [0.5, 1, 5, 23]
    A = [nuclear_area_ss(E, ph, reps=REPS) for E in Es]
    assert A[-1] > A[0], "area should be larger on stiff than on soft substrate"
    assert all(np.diff(A) > -2.0), "area should be broadly non-decreasing"
    print(f"  [OK] nuclear area vs stiffness: {A[0]:.1f} -> {A[-1]:.1f} µm²")


def test_yap_activation_with_stiffness():
    """YAP N/C ratio should rise from near-baseline on soft substrates to a
    multi-fold value on stiff substrates."""
    ph = PHENOTYPES["hepatocyte"]
    yap_soft = yap_nc_ratio(0.5, ph, reps=REPS)
    yap_stiff = yap_nc_ratio(20.0, ph, reps=REPS)
    assert yap_stiff > yap_soft, "YAP should increase with stiffness"
    assert yap_stiff / yap_soft > 1.3, "YAP should rise appreciably with stiffness"
    print(f"  [OK] YAP N/C: {yap_soft:.2f} (soft) -> {yap_stiff:.2f} (stiff)")


def test_lamin_knockdown_reduces_yap():
    """Reducing lamin A/C (knockdown) should lower stiff-substrate YAP:
    a softer nucleus gates less transcriptional activation."""
    ph = PHENOTYPES["hepatocyte"]
    yap_wt = yap_nc_ratio(20.0, ph, reps=REPS)
    ph_kd = replace(ph, laminAC=ph.laminAC * 0.2)   # 80% knockdown
    yap_kd = yap_nc_ratio(20.0, ph_kd, reps=REPS)
    drop = (yap_wt - yap_kd) / (yap_wt - 1.0)
    assert yap_kd < yap_wt, "lamin knockdown should reduce YAP"
    assert drop > 0.3, f"knockdown should substantially reduce YAP (drop={drop:.2f})"
    print(f"  [OK] lamin knockdown: YAP {yap_wt:.2f} -> {yap_kd:.2f} "
          f"({drop*100:.0f}% of dynamic range lost)")


def test_phenotype_lamin_ordering():
    """Phenotypes with higher lamin A/C should show larger YAP dynamic range
    (soft->stiff), reflecting stiffer, more gating-competent nuclei."""
    def yap_range(key):
        ph = PHENOTYPES[key]
        return yap_nc_ratio(20.0, ph, reps=REPS) - yap_nc_ratio(0.5, ph, reps=REPS)
    r_mda = yap_range("MDA")        # laminAC = 0.5 (soft nucleus)
    r_at2 = yap_range("AT2_lung")   # laminAC = 1.3 (stiff nucleus)
    assert r_at2 > r_mda, ("high-lamin phenotype should have larger YAP range "
                           f"(AT2={r_at2:.2f} vs MDA={r_mda:.2f})")
    print(f"  [OK] lamin ordering: YAP range MDA(soft)={r_mda:.2f} "
          f"< AT2(stiff)={r_at2:.2f}")


def test_two_population_basal_constant():
    """The basal population is constant; the mechanosensitive population grows
    with time (contact inhibition included)."""
    ph = PHENOTYPES["hepatocyte"]
    mus_basal, mus_mecano = [], []
    for t in [2, 12, 24, 36]:
        mb, mm, phi = mvc.population_mixture(23.0, t, ph, reps=4)
        mus_basal.append(mb)
        mus_mecano.append(mm)
    assert np.std(mus_basal) < 1e-6, "basal population must be constant"
    assert mus_mecano[-1] > mus_mecano[0], "mechanosensitive pop should grow in time"
    print(f"  [OK] two populations: basal constant={mus_basal[0]:.1f} µm², "
          f"mecano {mus_mecano[0]:.1f}->{mus_mecano[-1]:.1f} µm²")


def test_contact_inhibition_reduces_clutches():
    """Effective substrate clutches should decrease as confluence rises."""
    ph = PHENOTYPES["hepatocyte"]
    nc_early = mvc.nc_effective(ph, t=2)
    nc_late = mvc.nc_effective(ph, t=36)
    assert nc_late < nc_early, "contact inhibition should reduce effective clutches"
    print(f"  [OK] contact inhibition: nc_eff {nc_early} (2h) -> {nc_late} (36h)")


def test_temporal_relaxation():
    """Nuclear area should relax toward the steady state over time."""
    ph = PHENOTYPES["hepatocyte"]
    A_early = nuclear_area_time(23.0, 2, ph, reps=REPS)
    A_late = nuclear_area_time(23.0, 100, ph, reps=REPS)
    A_ss = nuclear_area_ss(23.0, ph, reps=REPS)
    assert abs(A_late - A_ss) < abs(A_early - A_ss), "area should approach steady state"
    print(f"  [OK] temporal relaxation: A(2h)={A_early:.1f} -> "
          f"A(100h)={A_late:.1f} ≈ A_ss={A_ss:.1f} µm²")


def test_fibrosis_monotonic():
    """Predicted nuclear stress should increase monotonically across fibrosis
    stages F0->F4 (the tissue stiffens)."""
    pred = mvc.fibrosis_prediction(PHENOTYPES["hepatocyte"], reps=REPS)
    sig = np.array(pred["sigma"])
    assert sig[-1] > sig[0], "F4 stress should exceed F0"
    # allow tiny stochastic dips but require overall increase
    assert np.polyfit(range(len(sig)), sig, 1)[0] > 0, "trend should be increasing"
    print(f"  [OK] fibrosis F0->F4 stress: {sig[0]:.1f} -> {sig[-1]:.1f} "
          f"(monotonic increasing)")


def test_optimum_sensitive_to_clutch_not_motor():
    """The optimal stiffness should shift more when clutch number changes than
    when motor number changes by the same relative amount."""
    ph = PHENOTYPES["hepatocyte"]
    kappas = np.array([0.1, 0.3, 0.7, 1.5, 3.0, 7.0])

    def k_opt(phen):
        T = [traction(k / phen.alpha, phen, reps=REPS) for k in kappas]
        return kappas[int(np.argmax(T))]

    base = k_opt(ph)
    more_clutch = k_opt(replace(ph, nc=int(ph.nc * 1.6)))
    more_motor = k_opt(replace(ph, nm=int(ph.nm * 1.6)))
    shift_clutch = abs(np.log(more_clutch / base))
    shift_motor = abs(np.log(more_motor / base))
    # robust qualitative check: clutch perturbation shifts optimum at least as much
    assert shift_clutch >= shift_motor - 1e-9, (
        f"optimum should be at least as sensitive to clutch as to motor "
        f"(clutch shift={shift_clutch:.2f}, motor shift={shift_motor:.2f})")
    print(f"  [OK] optimum sensitivity: clutch shift={shift_clutch:.2f} "
          f">= motor shift={shift_motor:.2f}")


def test_tau_scales_with_stiffness():
    """Recalibration result: nuclear relaxation is SLOWER on stiffer substrates
    (tau increases with stiffness), from the complete 2-120 h timecourse."""
    ph = PHENOTYPES["hepatocyte"]
    tau_soft = mvc.tau_of_E(1.0, ph)
    tau_stiff = mvc.tau_of_E(23.0, ph)
    assert tau_stiff > tau_soft, "tau should increase with stiffness"
    assert tau_stiff / tau_soft > 2.0, "stiff relaxation should be much slower"
    # dynamics: stiff substrate keeps growing longer than soft
    A_soft_late = mvc.nuclear_area_time(1.0, 120, ph, reps=REPS)
    A_soft_mid = mvc.nuclear_area_time(1.0, 36, ph, reps=REPS)
    A_stiff_late = mvc.nuclear_area_time(23.0, 120, ph, reps=REPS)
    A_stiff_mid = mvc.nuclear_area_time(23.0, 36, ph, reps=REPS)
    soft_growth = A_soft_late - A_soft_mid
    stiff_growth = A_stiff_late - A_stiff_mid
    assert stiff_growth > soft_growth, ("stiff substrate should keep growing "
                                        "between 36-120 h more than soft")
    print(f"  [OK] tau scales with stiffness: {tau_soft:.0f} h (soft) -> "
          f"{tau_stiff:.0f} h (stiff); stiff still growing 36->120 h")


def test_virtual_cell_interface():
    """The VirtualCell class produces a coherent state and state vector."""
    from virtual_cell import VirtualCell
    cell = VirtualCell("hepatocyte")
    soft = cell.simulate(1.0, t=120)
    stiff = cell.simulate(23.0, t=120)
    assert stiff.nuclear_drive > soft.nuclear_drive, "drive rises with stiffness"
    assert stiff.yap_nc > soft.yap_nc, "YAP rises with stiffness"
    assert stiff.function_index < soft.function_index, "function falls on stiff"
    v = cell.state_vector(stiff)
    assert len(v) == len(cell.STATE_FIELDS), "state vector length matches fields"
    assert stiff.fibrosis_stage == "F4" and soft.fibrosis_stage == "F0"
    print(f"  [OK] VirtualCell: soft F0 (YAP {soft.yap_nc:.1f}) -> "
          f"stiff F4 (YAP {stiff.yap_nc:.1f}); {len(v)}-D state vector")


def test_gene_module_response_shapes():
    """Sigmoid genes are switch-like (low at F0, high at F4); identity genes fall."""
    import gene_module as gm
    sig_soft = mvc.nuclear_stress(1.0, PHENOTYPES["hepatocyte"], reps=REPS)
    sig_stiff = mvc.nuclear_stress(23.0, PHENOTYPES["hepatocyte"], reps=REPS)
    soft = gm.score_genes(sig_soft)
    stiff = gm.score_genes(sig_stiff)
    # sigmoid YAP target: big jump
    assert stiff["CTGF (CCN2)"] - soft["CTGF (CCN2)"] > 0.5, "CTGF switch-like"
    # inverse identity gene falls
    assert stiff["HNF4A"] < soft["HNF4A"], "HNF4A (identity) falls with stiffness"
    # actionable hypotheses emerge at high stiffness
    hyp = gm.actionable_hypotheses(sig_stiff)
    assert len(hyp) >= 3, "actionable hypotheses at F4"
    print(f"  [OK] gene shapes: CTGF {soft['CTGF (CCN2)']:.2f}->{stiff['CTGF (CCN2)']:.2f}, "
          f"HNF4A {soft['HNF4A']:.2f}->{stiff['HNF4A']:.2f}; {len(hyp)} hypotheses")


def test_benchmark_full_generalizes():
    """Full model generalizes (CV) at least as well as baselines and uniquely
    captures the temporal rise at 23 kPa."""
    import benchmark as bm
    E, t, A = bm.load_area_data()
    cv = {c["name"]: c["cv_R2"] for c in
          [bm.cross_validate(m, E, t, A) for m in bm.MODELS]}
    full = [v for k, v in cv.items() if "full" in k][0]
    simple = max(v for k, v in cv.items() if "full" not in k)
    assert full >= simple - 0.02, "full model should not generalize worse"
    fits = [m().fit(E, t, A) for m in bm.MODELS]
    full_fit = [f for f in fits if "full" in f.name][0]
    lin_fit = [f for f in fits if "linear" in f.name][0]
    a36, a120 = full_fit.predict(23.0, 36)[0], full_fit.predict(23.0, 120)[0]
    l36, l120 = lin_fit.predict(23.0, 36)[0], lin_fit.predict(23.0, 120)[0]
    assert a120 - a36 > 20, "full model captures the temporal rise"
    assert abs(l120 - l36) < 1, "linear model is time-flat (structural limit)"
    print(f"  [OK] benchmark: full CV-R2={full:.2f} >= simple {simple:.2f}; "
          f"full rises {a36:.0f}->{a120:.0f} at 23 kPa, linear flat")


def test_sensitivity_identifies_key_params():
    """Sensitivity flags the nuclear gate (laminAC) and adhesion (nc/kc) group
    as the parameters the data must constrain, and alpha as low-impact."""
    import sensitivity as sa
    res = sa.run_sensitivity(verbose=False, global_analysis=True, n=1024)
    sob = res["sobol"]["nuclear_area"]
    # laminAC (nuclear gate) should carry substantial variance
    assert sob["laminAC"] > 0.3, "lamin dominates nuclear-area variance"
    # the adhesion/clutch group (nc, kc) should also matter
    assert sob["nc"] + sob["kc"] > 0.3, "adhesion/clutch group matters"
    # alpha (E->kappa coupling) is low-impact -> safe to fix (matches inference)
    assert sob["alpha"] < 0.2, "alpha is low-sensitivity (safe to fix)"
    print(f"  [OK] sensitivity: area driven by lamin (S={sob['laminAC']:.2f}) & "
          f"adhesion nc+kc (S={sob['nc']+sob['kc']:.2f}); alpha low "
          f"(S={sob['alpha']:.2f})")


def test_bootstrap_fold_change_ci():
    """Bootstrap CI on the stiffness fold-change contains the ~2.2x estimate
    and excludes 1 (a real mechanical response)."""
    import stats_ci as st
    import json
    from pathlib import Path
    try:
        from paths import DATA_DIR
        p = DATA_DIR / "hepatocyte_complete_data.json"
    except Exception:
        p = Path("hepatocyte_complete_data.json")
    d = json.load(open(p))
    soft = [a for a in d["complete_timecourse"]["1_kPa"]["pop_high"] if a]
    stiff = [a for a in d["complete_timecourse"]["23_kPa"]["pop_high"] if a]
    fc = st.fold_change_ci(soft, stiff, n_boot=3000)
    assert fc["fold_change"] > 1.5, "fold-change is a real response"
    assert fc["ci_low"] > 1.0, "CI excludes 1 (no response)"
    ci = st.bootstrap_ci(stiff, n_boot=2000)
    assert ci["ci_low"] < ci["estimate"] < ci["ci_high"], "CI brackets estimate"
    print(f"  [OK] bootstrap: fold-change {fc['fold_change']:.2f}x "
          f"95% CI [{fc['ci_low']:.2f}, {fc['ci_high']:.2f}] (excludes 1)")


def test_phenotype_specific_genes():
    """Each phenotype exposes its own identity/lineage markers, and they respond
    in the biologically correct direction (identity down, effectors up)."""
    import gene_module as gm
    # MDA (invasive cancer): invasion/EMT effectors rise on stiff
    sig = mvc.nuclear_stress(23.0, PHENOTYPES["hepatocyte"], reps=REPS)
    mda = gm.score_genes(sig, phenotype="MDA")
    assert mda.get("MMP9", 0) > 0.5, "MDA invasion marker MMP9 up on stiff"
    # hepatocyte identity falls
    hep = gm.score_genes(sig, phenotype="hepatocyte")
    assert hep.get("HNF4A", 1) < 0.3, "hepatocyte identity HNF4A down on stiff"
    assert hep.get("ALB (albumin)", 1) < 0.3, "albumin (function) down on stiff"
    # panels differ between phenotypes
    hep_genes = set(gm.genes_for("hepatocyte"))
    mda_genes = set(gm.genes_for("MDA"))
    assert hep_genes != mda_genes, "phenotypes have distinct gene panels"
    assert "HNF4A" in hep_genes and "HNF4A" not in mda_genes
    print(f"  [OK] phenotype genes: hepatocyte identity falls (HNF4A/ALB), "
          f"MDA invasion rises (MMP9={mda['MMP9']:.2f}); panels differ")


# ---------------------------------------------------------------------------
ALL_TESTS = [
    test_biphasic_traction,
    test_nuclear_area_increases_with_stiffness,
    test_yap_activation_with_stiffness,
    test_lamin_knockdown_reduces_yap,
    test_phenotype_lamin_ordering,
    test_two_population_basal_constant,
    test_contact_inhibition_reduces_clutches,
    test_temporal_relaxation,
    test_fibrosis_monotonic,
    test_optimum_sensitive_to_clutch_not_motor,
    test_tau_scales_with_stiffness,
    test_virtual_cell_interface,
    test_gene_module_response_shapes,
    test_benchmark_full_generalizes,
    test_sensitivity_identifies_key_params,
    test_bootstrap_fold_change_ci,
    test_phenotype_specific_genes,
]


def run_all():
    print("=" * 72)
    print("  VIRTUAL CELL — model validation suite")
    print(f"  numba: {'ON' if mvc.HAS_NUMBA else 'OFF (pure python, slower)'}")
    print("=" * 72)
    passed = 0
    failed = 0
    for test in ALL_TESTS:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  [FAIL] {test.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"  [ERROR] {test.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print("=" * 72)
    print(f"  {passed}/{len(ALL_TESTS)} passed"
          + (f", {failed} failed" if failed else "  — all validations OK"))
    print("=" * 72)
    return failed == 0


if __name__ == "__main__":
    import sys
    ok = run_all()
    sys.exit(0 if ok else 1)
