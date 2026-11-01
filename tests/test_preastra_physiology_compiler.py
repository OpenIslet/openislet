"""Independent negative controls for PREASTRA L4/L5 Category-A representations.

These tests execute the models and source receipts. They reject false completion of
unknown genome scans/yield and do not qualify any wet-lab or clinical hypothesis.
"""
from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import run_layer4_tissue_and_metabolic_twin as l4
from scripts import run_layer5_personalized_patient_compiler as l5
from scripts import rank_islet1_stacks as ranker
from scripts.engineering_priors import engineering_complexity_discount


def require(condition, message):
    if not condition:
        raise AssertionError(message)


@pytest.mark.parametrize("override", [
    {"phi_islet_in_slab": value} for value in (-.1, 0., 1.2, float("nan"), float("inf"))
] + [{"radius_cm": 0.}, {"radius_cm": -.1}, {"mean_po2_mmhg": -1.},
     {"mean_po2_mmhg": float("nan")}, {"mean_po2_mmhg": float("inf")},
     {"n_per_zone": 1}, {"n_per_zone": 2.5}, {"n_per_zone": True},
     {"supply_closure": "unrecognized"}, {"d_hydrogel_override": 0.}])
def test_l4_invalid_domains_rejected_before_zero_or_clamp_branch(override):
    args = dict(mean_po2_mmhg=19.93, radius_cm=.006)
    args.update(override)
    with pytest.raises(ValueError):
        l4._solve_two_zone_fv_wigner_seitz_bvp(**args)


def test_l4_phi_one_uniform_preserves_homogeneous_limit():
    result = l4._solve_two_zone_fv_wigner_seitz_bvp(19.93, .006, phi_islet_in_slab=1.)
    require(result["converged"], "homogeneous BVP did not converge")
    require(result["core_po2_mmhg"] == result["surface_po2_mmhg"] == 19.93, "phi1 homogeneous limit drift")
    require(abs(result["mean_constraint_error_mmhg"]) < 1e-10, "mean constraint failed")


@pytest.mark.parametrize("nz", (51, 101, 201, 401))
def test_l4_one_physical_lamina_conserves_profiles_on_every_mesh(nz):
    geometry = l4._slab_control_volume_profiles(nz)
    weights = geometry["widths"] / .035
    require(abs(float(np.dot(geometry["central_mask"], weights))*7.-2.38) < 1e-12, "mesh changes physical lamina volume")
    for key in ("central_source", "outer_demand"):
        require(abs(float(np.dot(geometry[key], weights))-1.) < 1e-10, "assigned continuum profile not normalized")
    alg = l4.verify_exact_slab_source_distribution_algebra()
    cad = alg["non_overlapping_phase_volume_and_cad_ledger_7cm3_ultracode_v80_18"]["layout_option_b_31wt_pct_cao2_confined_central_lamina_microbeads"]
    require(cad["central_lamina_total_envelope_volume_cm3"] == 2.38, "CAD and PDE geometry differ")
    require(alg["selected_source_lamina_half_width_fraction"] == l4.CENTRAL_LAMINA_HALF_WIDTH_FRACTION, "algebra uses different width")


@pytest.fixture(scope="module")
def near_threshold():
    return l4._run_single_pouch_pde_regime("PREASTRA_CRITICAL_F0447", .447, "central", 68., True, .97)


def test_l4_critical_gate_uses_unrounded_refinement_allowance(near_threshold):
    assessment = near_threshold["strict_gel_only_numerical_gate_assessment"]
    require(assessment["refinement_applied"], "critical run was not refined")
    require(assessment["reference_mesh_nodes"] == [201,401], "critical reference grids missing")
    lo, hi = assessment["reference_core_interval_mmhg"]
    passed = near_threshold["passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate"]
    require(not passed or lo >= 2., "strict PASS ignores refinement uncertainty")
    require(assessment["reference_401_core_unrounded_mmhg"] >= 2. or not passed, "coarse PASS masks resolved reference failure")
    require(assessment["estimated_core_discretization_error_mmhg"] > 0, "missing numerical allowance")
    require(lo <= assessment["reference_401_core_unrounded_mmhg"] <= hi, "reference outside allowance")


def test_l4_fast_catalase_loss_remains_real_failed_chemistry_control():
    result = l4._run_single_pouch_pde_regime("PREASTRA_CATALASE_1H", .45, "central", 68., True, .97, catalase_half_life_h=1.)
    require(result["max_interstitial_h2o2_uM"] > 60., "enzyme-loss control lost peroxide coupling")
    require(result["min_slab_po2_across_14_days_mmhg"] < 1., "enzyme-loss control lost oxygen coupling")
    require(not result["passes_all_oxygen_peroxide_and_ph_gates"], "failed chemistry passed product gate")


@pytest.fixture(scope="module")
def failed_pde_ode():
    return l4.solve_24h_paracrine_and_whole_body_metabolic_ode([{
        "regime_id": "PREASTRA_EXPLICIT_FAILED_PDE", "min_120um_cluster_core_po2_mmhg": .1,
        "min_slab_po2_across_14_days_mmhg": 1., "max_interstitial_h2o2_uM": 1.,
        "max_interstitial_ph": 7.4, "coupled_dynamic_14d_viable_endocrine_fraction_v80_066": .02,
        "passes_all_oxygen_peroxide_and_ph_gates": False,
    }])


def test_l4_reference_curves_cannot_impersonate_matched_pde(failed_pde_ode):
    for curve in failed_pde_ode["graft_architectures_compared"]:
        require(curve["pde_coupling_applied"] is False, "independent reference mislabeled PDE coupled")
        require(curve["evaluation_scope"] == "REFERENCE_FIXED_FUNCTIONAL_BETA_IEQ_NOT_PDE_DERIVED", "reference scope missing")
        require(curve["source_pde_regime_id"] is None, "reference inherited unrelated PDE")
        require("viable_endocrine_fraction_coupled" not in curve, "reference fraction mislabeled coupled")
        require(curve["assumed_reference_functional_fraction"] == 1., "fixed-dose reference changed")
    table = failed_pde_ode["matched_pde_ode_qualification_scenarios"]
    require(len(table) == 1, "coupled table must preserve named PDE rows")
    row = table[0]
    require(row["source_pde_regime_id"] == "PREASTRA_EXPLICIT_FAILED_PDE", "matched PDE provenance missing")
    require(row["coupled_viable_endocrine_fraction"] <= .02, "failed PDE viability did not reach ODE")
    require(row["coupled_24h_time_in_range_pct"] == 0., "failed matched ODE unexpectedly passes")
    require(not row["coupled_complete_product_pass"], "failed PDE produced complete-product PASS")


@pytest.fixture(scope="module")
def guide_atlas():
    return l5.build_part_a_crispr_grna_and_polymorphism_atlas()


def test_l5_unperformed_genome_scan_is_unknown_and_zero_needs_receipt(guide_atlas):
    for guide in guide_atlas["guides"]:
        require(guide["grch38_0_1_2_mismatch_coding_off_targets"] is None, "unenumerated count is a false zero")
        require(guide["genomewide_off_target_enumeration_status"] == "NOT_ENUMERATED", "scan status false completion")
        require(guide["off_target_specificity_breakdown_claude_v80_006"]["coding_exome_0_1_2_mismatch_hits"] is None, "nested unknown lost")
    with pytest.raises(ValueError):
        l5.compute_panel_and_reference_off_target_specificity("probe", "Cas9", 0, 1, [.05], 65.)


def test_l5_retains_all_source_identities_and_scopes_frequency_priors(guide_atlas):
    raw = json.loads(l5.RAW_LAYER5_DBSNP_PATH.read_text())
    audit = guide_atlas["gnomad_v4_1_and_dbsnp_guide_overlap_audit_claude_v80_045"]
    require("highest_frequency_overlapping_variants" not in json.dumps(audit), "unscoped frequency ranking survived")
    # The emitted source rows must still be read directly from the pinned 56-row receipt.
    source = raw["ucsc_dbsnp155_primary_interval_responses"]
    emitted = audit["per_guide_interval_polymorphism_summary"]
    total = 0
    for guide_id, entry in emitted.items():
        rows = entry["ucsc_dbsnp155_variant_rows"]
        expected = [{key: item[key] for key in ("rsid","chrom","chromStart","chromEnd","ref","alts","func","varType")}
                    for item in source[guide_id]["items"]]
        require(rows == expected, "primary variant identity/coordinate altered")
        total += len(rows)
        require(entry["assumed_screening_frequency_priors"]["source_class"] == "ASSUMED_NOT_DATABASE_MEASUREMENT", "assumed AF lacks class")
        for key, value in entry.items():
            if key.startswith("gnomad_v4_1_"):
                require(value is None, "measured frequency alias has no AC/AN receipt")
    require(total == 56, "source-derived variant count drift")


def patient(**overrides):
    obj = {"patient_id":"PREASTRA", "ancestry":"EUR", "nkg2a_neg_kir_educated_nk_fraction":.36,
           "pre_existing_gad65_ia2_znt8_t_cell_zscore":1.1,"ifn1_ish_genetic_risk_percentile":45.,"sirpa_genotype":"v1/v1"}
    obj.update(overrides)
    return obj


def test_l5_shared_registry_discount_separate_from_physical_yield(monkeypatch):
    monkeypatch.setitem(ranker.ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["locus_dsb_penalty_rate_lambda_K"], "point_estimate", .21)
    monkeypatch.setitem(ranker.ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["base_edit_nickase_penalty_factor_gamma_B"], "point_estimate", .93)
    for nk in (.36,.91):
        for ifn in (45.,95.):
            compiled = l5.compile_personalized_patient_blueprint(patient(nkg2a_neg_kir_educated_nk_fraction=nk, ifn1_ish_genetic_risk_percentile=ifn))
            counts = compiled["editing_complexity_and_manufacturing_accounting_v80_009"]
            k,b = counts["k_dsb_or_indel_chromosomal_loci"], counts["b_prime_or_base_edit_loci"]
            expected = ranker.engineering_complexity_discount(k,b,.21,.93)
            require(abs(counts["ranker_engineering_complexity_discount_factor"] - expected) < 1e-12, "L1/L5 discount registry drift")
            require(counts["uncloned_simultaneous_pass_yield_prior"] is None, "utility was reused as physical yield")
            require(counts["physical_manufacturing_yield_evidence_status"] == "NOT_MEASURED", "physical yield status false measured")
            product = compiled["product_definition"]
            require(product["K"] == k and product["B"] == b, "product identity/edit count divergence")
    require(engineering_complexity_discount(4,0) == 1., "K4 B0 anchor broken")
    require(abs(engineering_complexity_discount(5,0) - math.exp(-.15)) < 1e-12, "K5 anchor broken")


@pytest.fixture(scope="module")
def population():
    return l5.simulate_10000_patient_global_population_trial()


def test_l5_dynamic_scope_matches_actual_values_and_analytic_primary(population):
    account = population["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    text = account["scope_reconciliation_note"]
    for key in ("conditional_on_release_qc_mean_personalized_survival_pct", "conditional_on_release_qc_mean_standardized_plan_a_survival_pct", "ranker_complexity_discounted_personalized_utility_pct"):
        require(f"{account[key]:.4f}%" in text, "scope text drifted from computed scalar")
    require(account["physical_uncloned_yield_prior"] is None, "population utility mistaken for yield")
    require(len(account["engineering_utility_prior_sensitivity"]) == 4, "range endpoint sensitivity missing")
    require(population["analytic_max_ancestry_disparity_gap_personalized_pct_v80_065"] == .4, "analytic primary drift")
    require("analytic primary" in population["ancestry_disparity_gap_methodology_note_v80_065"], "analytic/sample gap ambiguity")
    require("sampled secondary" in population["ancestry_disparity_gap_methodology_note_v80_065"], "sample qualification missing")


def test_l5_only_prime_editor_is_executable_optional_tyk2(guide_atlas):
    guide = next(g for g in guide_atlas["guides"] if "TYK2" in g["guide_id"])
    require("PEMAX" in guide["guide_id"], "optional executable editor changed")
    require(guide["cgbe_alternative_implementation_status"] == "NOT_IMPLEMENTED_REQUIRES_ENGINEERED_PAM_WINDOW_AND_BYSTANDER_LEDGER", "CGBE contingency falsely completed")
    require(len(guide["epegRNA_extension_specification"]["primer_binding_site_pbs_13nt_5to3"]) == 13 and len(guide["epegRNA_extension_specification"]["rt_template_rtt_12nt_5to3"]) == 12, "PEmax sequence arithmetic regressed")


def test_l5_report_summary_uses_analytic_primary_and_preserves_unknown_count(monkeypatch, population, guide_atlas):
    monkeypatch.setattr(l5, "build_part_a_crispr_grna_and_polymorphism_atlas", lambda: guide_atlas)
    monkeypatch.setattr(l5, "simulate_10000_patient_global_population_trial", lambda: population)
    summary = l5.build_layer5_atlas()["summary"]
    require(summary["coding_off_targets_0_to_2_mismatches"] is None, "unknown genome counts summed as zero")
    require(summary["analytic_max_ancestry_disparity_gap_personalized_pct"] == population["analytic_max_ancestry_disparity_gap_personalized_pct_v80_065"], "summary selects sampled headline")
    require(summary["sampled_max_ancestry_disparity_gap_personalized_pct"] == population["sampled_max_ancestry_disparity_gap_personalized_pct"], "sampled field drift")
    require(summary["sampled_ancestry_gap_method"]["seed"] == 20261023 and summary["sampled_ancestry_gap_method"]["n_per_cohort"] == 2000, "sampling qualification missing")
