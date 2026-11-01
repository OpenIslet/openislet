"""Independent three-pass numerical and publication controls; fixtures are synthetic."""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import ingest_wetlab_results as w
from scripts import run_layer4_tissue_and_metabolic_twin as l4
from scripts import run_layer5_personalized_patient_compiler as l5
from scripts import wetlab_contract as contract


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def patient(**updates):
    value = dict(patient_id="SYNTHETIC-ASTRA", ancestry="EUR", nkg2a_neg_kir_educated_nk_fraction=.4,
                 pre_existing_gad65_ia2_znt8_t_cell_zscore=1., ifn1_ish_genetic_risk_percentile=40., sirpa_genotype="v1/v1")
    value.update(updates)
    return value


def observation(batch="SYNTHETIC-A", value=.91):
    return dict(batch_id=batch, obs_stage7_yield_ratio=value, obs_stage7_yield_se=.05)


@pytest.mark.parametrize("duration", [5e-13, 1e-15, 1e-100])
def test_rk4_integrates_every_representable_positive_interval(duration):
    actual = l4.rk4_integrate(lambda t, y: np.array([1./duration]), (0., duration), np.array([0.]), [duration], dt_max=duration)
    require(abs(actual[0, 0]-1.) < 1e-14, "absolute clock epsilon silently discarded a valid interval")


@pytest.mark.parametrize("value", [np.complex64(.4+2j), np.complex128(.4+2j), np.complex128(.4), np.bool_(True)])
def test_numpy_nonreal_or_boolean_patient_inputs_are_not_coerced(value):
    with pytest.raises(ValueError):
        l5.compile_personalized_patient_blueprint(patient(nkg2a_neg_kir_educated_nk_fraction=value))


@pytest.mark.parametrize("value", [np.float64(.4), np.float32(.4), np.int64(0)])
def test_numpy_real_patient_inputs_still_work(value):
    result = l5.compile_personalized_patient_blueprint(patient(nkg2a_neg_kir_educated_nk_fraction=value))
    require(result["patient_id"] == "SYNTHETIC-ASTRA", "valid numpy scalar rejected")


def test_hwe_degenerate_allele_frequencies_include_closed_quantile_boundary():
    for q in (0., .3, 1.):
        require(l5._hwe_alternate_alleles(q, 0.) == 0, "p=0 must have zero alternate alleles")
        require(l5._hwe_alternate_alleles(q, 1.) == 2, "p=1 must have two alternate alleles")


def test_copula_uses_same_cumulative_genotype_distribution_as_main_population():
    seen = {}
    def capture(record):
        # Copy only genotype receipts, avoiding dependence on assigned survival scores.
        if record["patient_id"].startswith("COP-"):
            _, ancestry, index = record["patient_id"].split("-")
            seen[(ancestry, int(index))] = record["tyk2_rs34536443_genotype"]
        return l5.compile_personalized_patient_blueprint(record)
    l5.simulate_10000_patient_global_population_trial(capture)
    frequencies = {"EUR": .042, "AFR": .008, "EAS": .002, "SAS": .015, "AMR": .024}
    for (ancestry, index), genotype in seen.items():
        p = frequencies[ancestry]
        q = ((index*997+len(ancestry)) % 10000)/10000.
        expected = "C/C" if q < p*p else "G/C" if q < 1.-(1.-p)**2 else "G/G"
        require(genotype == expected, f"Copula HWE mass differs at {ancestry}/{index}")
    require(len(seen) == 10000, "did not observe every copula patient")


@pytest.mark.parametrize("height", [.5, 25000., 1e6])
@pytest.mark.parametrize("n", [7, 15, 31])
def test_nodewise_catalase_matches_independent_manufactured_solution(n, height):
    half = .0175
    z = np.linspace(-half, half, n+2)[1:-1]
    truth = height*(1.-(z/half)**2)
    kcat, km, perfusion = 18., 25000., .01
    reaction = kcat*truth/(1.+truth/km)
    # The exact second derivative of the parabola is -2*height/half**2.
    source = 2*l4.D_H2O2_HYDROGEL_CM2_S*height/half**2 + reaction + perfusion*truth
    result = l4._solve_peroxide_slab_bvp(source, 2*half/(n+1), kcat, perfusion, km_uM=km)
    require(np.max(abs(result["concentration_uM"]-truth)) < max(1e-7, height*1e-8), "MM solve does not satisfy the declared nodewise equation")
    require(np.max(abs(result["decomposition_uM_s"]-reaction)) < max(1e-6, float(max(reaction))*1e-8), "oxygen source is not the same local decomposition reaction")
    require(result["converged"] is True, "valid manufactured solution did not converge")


def test_zero_oxygen_regime_is_a_defined_failed_control_without_artificial_floor():
    result = l4._run_single_pouch_pde_regime("SYNTHETIC-ZERO-OXYGEN", 0., "none", 68., True, 0.,
                                           boundary_po2_mmhg=0., nz_override=5, _assess_numerical_gate=False)
    require(result["passes_all_oxygen_peroxide_and_ph_gates"] is False, "zero oxygen passed")
    require(all(row["center_z0_po2_mmhg"] == 0. for row in result["trajectory_0h_to_336h"]), "unphysical positive oxygen floor survives")
    require(result["slab_oxygen_solver_diagnostics"]["max_residual_inf_mmhg_s"] == 0., "zero solution is not exact")
    json.dumps(result, allow_nan=False)


def test_nonconverged_slab_cannot_pass_with_favorable_spherical_solver(monkeypatch):
    original = np.linalg.solve
    def zero_oxygen_step(matrix, rhs):
        # Oxygen Newton systems have a diagonal MM derivative and tiny rhs;
        # identify them from source context without replacing peroxide helper.
        import inspect
        frame = inspect.currentframe().f_back
        if frame and frame.f_code.co_name == "_run_single_pouch_pde_regime" and frame.f_lineno:
            if "mat" in frame.f_locals and matrix is frame.f_locals["mat"]:
                return np.zeros_like(rhs)
        return original(matrix, rhs)
    monkeypatch.setattr(np.linalg, "solve", zero_oxygen_step)
    result = l4._run_single_pouch_pde_regime("SYNTHETIC-NONCONVERGED-SLAB", .45, "central", 68., True, .97,
                                           nz_override=5, _assess_numerical_gate=False)
    require(result["min_slab_po2_across_14_days_mmhg"] == 40., "fixture did not retain favorable initial guess")
    require(result["slab_oxygen_solver_diagnostics"]["converged"] is False, "nonzero slab residual hidden")
    require(result["passes_all_oxygen_peroxide_and_ph_gates"] is False, "unresolved slab passed qualification")


def coupling_row(**updates):
    value = dict(regime_id="SYNTHETIC-PDE", min_120um_cluster_core_po2_mmhg=4.,
                 min_slab_po2_across_14_days_mmhg=20., max_interstitial_h2o2_uM=1., max_interstitial_ph=7.4,
                 coupled_dynamic_14d_viable_endocrine_fraction_v80_066=1., passes_all_oxygen_peroxide_and_ph_gates=True,
                 passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate=True)
    value.update(updates)
    return value


@pytest.mark.parametrize("field,value", [
    ("max_interstitial_h2o2_uM", float("nan")), ("max_interstitial_h2o2_uM", -1.),
    ("max_interstitial_ph", float("nan")), ("max_interstitial_ph", 15.),
    ("coupled_dynamic_14d_viable_endocrine_fraction_v80_066", float("nan")),
    ("coupled_dynamic_14d_viable_endocrine_fraction_v80_066", 1.01),
    ("passes_all_oxygen_peroxide_and_ph_gates", "False"),
    ("passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate", 1),
    ("min_slab_po2_unrounded_mmhg", float("nan")), ("min_120um_cluster_core_po2_mmhg", True),
])
def test_imported_pde_state_fails_before_metabolic_execution(monkeypatch, field, value):
    def forbidden(*args, **kwargs):
        raise AssertionError("invalid import reached numerical integration")
    monkeypatch.setattr(l4, "rk4_integrate", forbidden)
    with pytest.raises(ValueError):
        l4.solve_24h_paracrine_and_whole_body_metabolic_ode([coupling_row(**{field: value})])


def test_pde_regime_identity_is_unique():
    with pytest.raises(ValueError, match="unique"):
        l4._validate_pde_coupling_regimes([coupling_row(), coupling_row()])


@pytest.mark.parametrize("negative", [-0., "-0", "-0.000e10"])
def test_signed_zero_cannot_replay_a_measurement_under_new_batch_id(negative):
    first = w.run_calibration(observation(value=0.), False)
    with pytest.raises(ValueError, match="replay"):
        w.run_calibration(observation("SYNTHETIC-B", negative), False, first)


def test_hash_helper_itself_canonicalizes_signed_zero():
    require(w.compute_measurement_payload_sha256(observation(value=-0.)) == w.compute_measurement_payload_sha256(observation(value=0.)), "standalone replay hash distinguishes signed zero")


@pytest.mark.parametrize("kind", ["fifo", "hardlink", "directory", "symlink"])
def test_cli_rejects_special_and_aliased_input_files(tmp_path, kind):
    target = tmp_path/"unsafe"
    original = tmp_path/"original.json"
    original.write_text("{}")
    if kind == "fifo": os.mkfifo(target)
    elif kind == "hardlink": os.link(original, target)
    elif kind == "directory": target.mkdir()
    else: target.symlink_to(original)
    with pytest.raises(ValueError):
        w._validate_safe_cli_path(str(target), "input")


def test_safe_read_rejects_duplicate_json_keys(tmp_path):
    target = tmp_path/"duplicates.json"
    target.write_text('{"value":1,"value":2}')
    with pytest.raises(ValueError, match="Duplicate"):
        w._read_safe_json_file(target)


def test_output_publication_rejects_stale_snapshots_without_changing_newer_state(tmp_path):
    target = tmp_path/"state.json"
    first = w.run_calibration(observation(), False)
    next_state = w.run_calibration(observation("SYNTHETIC-B", .93), False, first)
    w._atomic_write_json(target, first)
    snapshot = target.read_bytes()
    w._atomic_write_json(target, next_state)
    newer = target.read_bytes()
    with pytest.raises(ValueError, match="changed"):
        w._publish_calibration_report(target, next_state, {target: snapshot})
    require(target.read_bytes() == newer, "stale writer destroyed newer posterior")


def test_empirical_overwrite_requires_current_parent_report(tmp_path):
    target = tmp_path/"state.json"
    first = w.run_calibration(observation(), False)
    unrelated = w.run_calibration(observation("SYNTHETIC-B", .93), False)
    w._atomic_write_json(target, first)
    snapshot = target.read_bytes()
    with pytest.raises(ValueError, match="extend"):
        w._publish_calibration_report(target, unrelated, {target: snapshot})
    require(target.read_bytes() == snapshot, "unrelated initial prior silently replaced empirical history")


def test_correct_sequential_publication_still_succeeds(tmp_path):
    target = tmp_path/"state.json"
    first = w.run_calibration(observation(), False)
    w._publish_calibration_report(target, first, {target: None})
    second = w.run_calibration(observation("SYNTHETIC-B", .93), False, first)
    w._publish_calibration_report(target, second, {target: target.read_bytes()})
    require(json.loads(target.read_text())["ingested_batch_ids"] == ["SYNTHETIC-A", "SYNTHETIC-B"], "valid sequential update failed")


def test_alias_failure_rolls_back_canonical_report(tmp_path, monkeypatch):
    target, alias = tmp_path/"canonical.json", tmp_path/"legacy.json"
    monkeypatch.setattr(w, "DEFAULT_OUT", target)
    monkeypatch.setattr(w, "LEGACY_ALIAS_PATH", alias)
    first = w.run_calibration(observation(), False)
    for p in (target, alias): w._atomic_write_json(p, first)
    snapshots = {p: p.read_bytes() for p in (target, alias)}
    second = w.run_calibration(observation("SYNTHETIC-B", .93), False, first)
    original = w._atomic_write_json
    def fail_alias(path, data):
        if path == alias: raise OSError("synthetic second-target publication failure")
        original(path, data)
    monkeypatch.setattr(w, "_atomic_write_json", fail_alias)
    with pytest.raises(OSError, match="synthetic"):
        w._publish_calibration_report(target, second, snapshots)
    require(all(p.read_bytes() == snapshots[p] for p in snapshots), "alias failure left divergent posterior versions")


def test_concurrent_writer_lock_fails_closed(tmp_path):
    target = tmp_path/"state.json"
    with w._calibration_output_lock(target):
        with pytest.raises(ValueError, match="writer"):
            with w._calibration_output_lock(target):
                raise AssertionError("concurrent lock was acquired")


def test_cli_no_write_creates_neither_output_nor_lock(tmp_path):
    before = sorted(p.name for p in tmp_path.iterdir())
    result = subprocess.run([sys.executable, str(ROOT/"scripts/ingest_wetlab_results.py"), "--dry-run-calibration",
                             "--no-write", "--output", str(tmp_path/"must-not-exist.json")], capture_output=True, text=True)
    require(result.returncode == 0, result.stderr)
    require(sorted(p.name for p in tmp_path.iterdir()) == before, "read-only execution created publication state")


def test_layer5_duplicate_fasta_header_cannot_be_silently_overwritten(tmp_path, monkeypatch):
    source = (ROOT/"data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta").read_text()
    first = source.split(">", 2)[1]
    mutated = source+">"+first
    original = Path.read_text
    def fake_read(path, *args, **kwargs):
        if path == ROOT/"data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta": return mutated
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", fake_read)
    with pytest.raises((ValueError, RuntimeError), match="[Dd]uplicate"):
        l5.build_part_a_crispr_grna_and_polymorphism_atlas()


def test_parallel_spherical_diagnostics_remain_bound_to_their_own_solution():
    from concurrent.futures import ThreadPoolExecutor
    inputs = [0., 5., 25.] * 16
    with ThreadPoolExecutor(6) as pool:
        results = list(pool.map(lambda value: l4._solve_fixed_surface_spherical_bvp_with_diagnostics(value, .006), inputs))
    for (core, consumption, diagnostics), surface in zip(results, inputs):
        require(diagnostics["raw_core_po2_mmhg"] == core, "another concurrent solve overwrote diagnostics")
        require(0. <= core <= surface and consumption >= 0., "spherical solution violated nonnegative/max-principle bound")


@pytest.mark.parametrize("converged", [True, False])
def test_sensitivity_tables_do_not_promote_rounding_or_failed_solves(monkeypatch, converged):
    fixture = l4._run_single_pouch_pde_regime("SYNTHETIC-TABLE", .45, "central", 68., True, .97,
                                            _assess_numerical_gate=False)
    fixture.update(min_slab_po2_across_14_days_mmhg=15.0, min_slab_po2_unrounded_mmhg=14.999,
                   min_120um_cluster_core_po2_mmhg=0.5, min_120um_cluster_core_po2_unrounded_mmhg=.4999,
                   required_120um_oxygen_solvers_converged=converged)
    if not converged:
        # Even comfortably passing numbers cannot certify an unresolved solve.
        fixture.update(min_slab_po2_unrounded_mmhg=25., min_120um_cluster_core_po2_unrounded_mmhg=5.)
    monkeypatch.setattr(l4, "_run_single_pouch_pde_regime", lambda *a, **k: copy.deepcopy(fixture))
    result = l4.solve_336h_spatial_pouch_pde()
    rows = result["lead_regime_4_scenario_sensitivity_table_sol_v80_12"]
    for row in rows:
        require(row["passes_15mmHg_gsis_gate_at_2500_ieq_cm2"] is False, "rounded/unconverged slab promoted to pass")
        require(row["passes_0p5mmHg_anoxic_necrosis_gate"] is False, "rounded/unconverged core promoted to pass")
        require(row["passes_2p0mmHg_combined_core_safety_gate"] is False, "unconverged combined core promoted to pass")
        shrinking = row["shrinking_core_coupled_bvp_ultracode_v80_09"]
        require(shrinking["passes_15mmHg_gsis_gate_at_2500_ieq_cm2"] is False, "shrinking-core table promoted to pass")


def test_publication_rolls_back_failure_after_successful_rename(tmp_path, monkeypatch):
    import stat
    target = tmp_path / "posterior.json"
    first = w.run_calibration(observation(), is_dry_run=False)
    second = w.run_calibration(observation("SYNTHETIC-B", .92), is_dry_run=False, prior_posterior_report=first)
    target.write_text(json.dumps(first)); before = target.read_bytes()
    original, raised = w.os.fsync, [False]
    def post_rename_failure(fd):
        if stat.S_ISDIR(os.fstat(fd).st_mode) and not raised[0]:
            raised[0] = True
            raise OSError("SYNTHETIC post-rename parent fsync failure")
        return original(fd)
    monkeypatch.setattr(w.os, "fsync", post_rename_failure)
    with pytest.raises(OSError, match="post-rename"):
        w._publish_calibration_report(target, second, {target: before})
    require(target.read_bytes() == before, "failed publication left a new posterior after rename")


@pytest.mark.parametrize("solver,cutoff", [("fixed", 1e-6), ("two_zone", 1e-5), ("aggregate", 1e-4)])
@pytest.mark.parametrize("fraction", [0., .1, 1.])
def test_small_positive_oxygen_shortcut_is_not_exact_convergence(solver, cutoff, fraction):
    value = cutoff * fraction
    if solver == "fixed":
        core, uptake, diagnostics = l4._solve_fixed_surface_spherical_bvp_with_diagnostics(value, .006)
        require(core == 0. and uptake == 0., "shortcut approximation unexpectedly redesigned")
    elif solver == "two_zone":
        diagnostics = l4._solve_two_zone_fv_wigner_seitz_bvp(value, .006)
        require(diagnostics["mean_constraint_error_mmhg"] == -value, "prescribed mean mismatch hidden")
    else:
        result = l4.solve_spherical_cluster_po2_bvp(value, .006)
        diagnostics = result["spherical_bvp_solver_diagnostics"]
        require(l4._required_oxygen_solvers_converged(result) is (value == 0.), "positive shortcut qualified as solved")
    require(diagnostics["converged"] is (value == 0.), "positive approximation falsely claimed converged")
    require(diagnostics["constraint_error_mmhg"] == -value, "actual boundary or mean mismatch not reported")
    require(diagnostics["short_circuit_cutoff_mmhg"] == cutoff, "existing numerical cutoff changed")
    if value:
        require(diagnostics["residual_inf_mmhg_s"] is None, "unevaluated residual fabricated as exact zero")
    else:
        require(diagnostics["residual_inf_mmhg_s"] == 0., "exact zero control lost its exact residual")
