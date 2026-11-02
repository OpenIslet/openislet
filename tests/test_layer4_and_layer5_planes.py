#!/usr/bin/env python3
"""Executable consistency checks for conditional Layer 4/5 model artifacts.

Layer 4 checks conservative oxygen geometry, nonlinear catalase/peroxide/pH chemistry,
reference paracrine curves, and explicitly matched PDE-to-ODE viability states. These
calculations are not wet-lab validation. Layer 5 checks sequence/reference arithmetic,
unknown genome-wide off-target enumeration, and illustrative assigned-rule population
summaries; panel scores and utility discounts are not clinical or manufacturing estimates.
Assertions remain active under Python optimization through the explicit _check helper.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import generate_release_manifest as manifest_mod
from scripts import ingest_wetlab_results as wetlab
from scripts import run_layer4_tissue_and_metabolic_twin as l4_mod
from scripts import run_layer5_personalized_patient_compiler as l5_mod

L4_ATLAS = REPO_ROOT / "data" / "extracted" / "layer4_tissue_and_metabolic_twin_atlas.json"
L5_ATLAS = REPO_ROOT / "data" / "extracted" / "layer5_personalized_patient_and_crispr_atlas.json"
WETLAB_CALIB = REPO_ROOT / "data" / "extracted" / "wetlab_bayesian_posterior_calibration_v0.8.0.json"


def _check(cond: bool, msg: str = "Verification check failed") -> None:
    if not cond:
        raise AssertionError(msg)


def main() -> None:
    _check(L4_ATLAS.exists(), f"Missing {L4_ATLAS}")
    _check(L5_ATLAS.exists(), f"Missing {L5_ATLAS}")
    _check(WETLAB_CALIB.exists(), f"Missing {WETLAB_CALIB}")

    l4 = json.loads(L4_ATLAS.read_text(encoding="utf-8"))
    l5 = json.loads(L5_ATLAS.read_text(encoding="utf-8"))
    calib = json.loads(WETLAB_CALIB.read_text(encoding="utf-8"))

    # 1. Verify Layer 4 Part A exact algebra & 336h PDE (including AUDIT-REF & AUDIT-REF)
    _check(l4["summary"]["active_350um_pouch_volume_at_200cm2_cm3"] == 7.0, "Active pouch volume != 7.0 cm3")
    alg = l4["part_a_exact_oxygen_boundary_value_and_volume_algebra"]
    _check(1195.0 <= alg["passive_capacity_islet_only_ieq_per_cm2"] <= 1202.0, "Passive islet-only capacity out of range")
    _check(1040.0 <= alg["passive_capacity_with_15pct_iec_stroma_ieq_per_cm2"] <= 1045.0, "Passive +15% iEC capacity out of range")
    row_2500 = alg["exact_boundary_value_loading_table"][2]
    _check(row_2500["areal_loading_ieq_per_cm2"] == 2500.0, "Row 2 areal loading != 2500.0")
    _check(row_2500["active_350um_hydrogel_slab_volume_cm3"] == 7.0, "Row 2 active volume != 7.0")
    _check(abs(row_2500["f_chem_min_uniform_dispersion_islet_only"] - 0.5207) < 1e-3, "f_chem_min_uniform_dispersion_islet_only mismatch")
    _check(abs(row_2500["f_chem_min_central_core_sheet_islet_only"] - 0.3077) < 1e-3, "f_chem_min_central_core_sheet_islet_only mismatch")
    _check(
        abs(row_2500["f_chem_min_finite_width_lamina_a0p34_zero_order_including_15pct_iec_demand"] - 0.3544) < 1e-3,
        "Piecewise finite-width lamina threshold at 2500+iEC must equal 0.3544 (AUDIT-REF)",
    )
    row_3000 = alg["exact_boundary_value_loading_table"][3]
    _check(
        abs(row_3000["f_chem_min_finite_width_lamina_a0p34_zero_order_including_15pct_iec_demand"] - 0.4107) < 1e-3,
        "Piecewise finite-width lamina threshold at 3000+iEC must equal 0.4107 (AUDIT-REF)",
    )
    # AUDIT-REF: Verify phi->1 uniform supply limit in _solve_two_zone_fv_wigner_seitz_bvp returns mean_po2_mmhg
    ws_phi1 = l4_mod._solve_two_zone_fv_wigner_seitz_bvp(19.93, 0.006, phi_islet_in_slab=0.999999, supply_closure="uniform")
    _check(
        abs(ws_phi1["core_po2_mmhg"] - 19.93) < 1e-3 and abs(ws_phi1["mean_constraint_error_mmhg"]) < 1e-7,
        f"phi->1 uniform Wigner-Seitz limit failed: got {ws_phi1}",
    )

    pde_part = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    regimes = {r["regime_id"]: r for r in pde_part["regimes_simulated"]}
    _check(regimes["PASSIVE_NO_CAO2_CONTROL"]["passes_all_oxygen_peroxide_and_ph_gates"] is False, "PASSIVE_NO_CAO2_CONTROL should fail")
    _check(regimes["UNDERDOSED_CENTRAL_CAO2_F020_V70_12_TEST"]["passes_all_oxygen_peroxide_and_ph_gates"] is False, "UNDERDOSED_CENTRAL_CAO2_F020_V70_12_TEST should fail")
    _check(regimes["RAW_UNCOATED_CAO2_NO_CATALASE_CHEMICAL_TOXICITY"]["passes_all_oxygen_peroxide_and_ph_gates"] is False, "RAW_UNCOATED_CAO2_NO_CATALASE_CHEMICAL_TOXICITY should fail")
    _check(regimes["RAW_UNCOATED_CAO2_NO_CATALASE_CHEMICAL_TOXICITY"]["max_interstitial_h2o2_uM"] > 500.0, "Raw CaO2 H2O2 <= 500 uM")
    _check(7.44 <= regimes["RAW_UNCOATED_CAO2_NO_CATALASE_CHEMICAL_TOXICITY"]["max_interstitial_ph"] <= 7.60, "Raw CaO2 pH out of range")
    _check(regimes["FULL_LENGTH_CD47_WT_IEC_ANGIOGENESIS_IMPAIRED"]["passes_all_oxygen_peroxide_and_ph_gates"] is False, "FULL_LENGTH_CD47_WT_IEC_ANGIOGENESIS_IMPAIRED should fail")

    # AUDIT-REF: Isolated peroxide control at identical 28 mM/pH buffer capacity fails H2O2 gate while pH stays physiological
    no_cat_ctrl = regimes["PDMS_CAO2_F045_NO_CATALASE_ISOLATED_PEROXIDE_CONTROL"]
    _check(no_cat_ctrl["passes_all_oxygen_peroxide_and_ph_gates"] is False, "Isolated peroxide control should fail")
    _check(no_cat_ctrl["max_interstitial_h2o2_uM"] > 50.0, "Isolated peroxide control H2O2 <= 50 uM")
    _check(no_cat_ctrl["max_interstitial_ph"] <= 7.50, "Isolated peroxide control pH > 7.50")

    lead_pde = regimes["LEAD_CENTRAL_LAMINA_PDMS_CAO2_F045_PLUS_CATALASE"]
    _check(lead_pde["passes_all_oxygen_peroxide_and_ph_gates"] is True, "Lead PDE regime should pass all gates")
    _check(abs(lead_pde["min_slab_po2_across_14_days_mmhg"] - 19.93) < 0.05, f"Lead min_slab_po2 != 19.93 mmHg (got {lead_pde['min_slab_po2_across_14_days_mmhg']} — catches diffusivity mutant!)")
    _check(abs(lead_pde["min_120um_cluster_core_po2_mmhg"] - 8.65) < 0.05, f"Lead min_120um_cluster_core_po2 != 8.65 mmHg (got {lead_pde['min_120um_cluster_core_po2_mmhg']})")
    _check(lead_pde["min_300um_cluster_core_po2_mmhg"] == 0.0, "300 um cluster core pO2 should be 0.0 mmHg")
    _check(lead_pde["max_interstitial_h2o2_uM"] <= 5.0, "Lead max_interstitial_h2o2 > 5.0 uM")
    _check(lead_pde["max_interstitial_ph"] <= 7.60, "Lead max_interstitial_ph > 7.60")
    _check(abs(alg["chemical_o2_fraction_definitions_reconciled_sol_v80_12"]["f_chem_72h_time_averaged_integrated_fraction"] - 0.318821) < 1e-4, "f_chem_72h mismatch")
    # AUDIT-REF & AUDIT-REF: Exact CaO2 MW (72.078 g/mol -> 791.5 mg), H2O2->O2 coupling, 337-pt hourly grid, and shrinking-core comparator
    _check(abs(lead_pde["cao2_initial_mass_mg_in_7cm3_slab"] - 791.5) < 0.5, "CaO2 initial mass != 791.5 mg")
    _check(lead_pde["hour0_h2o2_to_o2_disproportionation_yield_fraction"] > 0.99, "Lead disproportionation yield <= 0.99")
    _check(no_cat_ctrl["hour0_h2o2_to_o2_disproportionation_yield_fraction"] < 0.10, "No-catalase disproportionation yield >= 0.10")
    _check(lead_pde["hourly_grid_extremum_evaluation_points"] == 337, "Hourly grid points != 337")
    _check("shrinking_core_2_3_remaining_fraction_at_36h_and_72h" in lead_pde, "Missing shrinking_core_2_3_remaining_fraction_at_36h_and_72h")
    # AUDIT-REF: Both fixed-sigma (139.8 um) and fixed-phi (221.2 um) critical passive slab thicknesses reported
    crit_slab = alg["critical_passive_slab_thickness_reconciliation_at_3000_ieq_cm2_sol_v80_r2_10"]
    _check(abs(crit_slab["fixed_areal_loading_sigma_3000_ieq_cm2_linear_scaling"]["h_crit_islet_only_um"] - 139.80) < 0.1, "Fixed-sigma h_crit != 139.80 um")
    _check(abs(crit_slab["fixed_volume_fraction_phi_0p1515_quadratic_scaling"]["h_crit_islet_only_um"] - 221.20) < 0.1, "Fixed-phi h_crit != 221.20 um")

    # AUDIT-REF: 4-scenario sensitivity table verifies conditional nature of oxygen pass
    sens_table = pde_part["lead_regime_4_scenario_sensitivity_table_sol_v80_12"]
    _check(len(sens_table) == 4, "Sensitivity table len != 4")
    _check(sens_table[0]["passes_15mmHg_gsis_gate_at_2500_ieq_cm2"] is True, "Sensitivity row 0 should pass")
    _check(sens_table[1]["passes_15mmHg_gsis_gate_at_2500_ieq_cm2"] is False, "Sensitivity row 1 should fail")
    _check(sens_table[2]["passes_15mmHg_gsis_gate_at_2500_ieq_cm2"] is False, "Sensitivity row 2 should fail")
    _check(sens_table[3]["passes_15mmHg_gsis_gate_at_2500_ieq_cm2"] is False, "Sensitivity row 3 should fail")

    # Fixed geometry and raw/refined conservative gate: one physical source zone on every mesh.
    grid_rows = pde_part["spatial_grid_convergence_51_101_201_401_preastra_a094"]
    _check(len(grid_rows) == 24, "Six scenarios must each retain all four mesh receipts")
    for row in grid_rows:
        ledger = row["spatial_control_volume_ledger"]
        _check(abs(ledger["lamina_volume_cm3"] - 2.38) < 1e-12, "Mesh changed physical lamina volume")
        _check(abs(ledger["assigned_chemical_source_relative_integral"] - 1.) < 1e-10, "Chemical source integral drift")
        _check(abs(ledger["assigned_total_demand_relative_integral"] - 1.) < 1e-10, "Demand integral drift")
    assessment = lead_pde["strict_gel_only_numerical_gate_assessment"]
    _check(assessment["refinement_applied"] is True, "Nominal conservative gate requires its refinement receipt")
    _check(assessment["reference_core_interval_mmhg"][0] >= 2., "Nominal conservative pass has no numerical margin")

    # DEC-080 / IIP-0002: Verify joint vascular timing-by-exchange sensitivity grid across both release laws
    iip2 = pde_part["joint_vascular_timing_and_exchange_sensitivity_iip_0002"]
    _check(iip2["proposal_id"] == "IIP-0002" and iip2["open_question_id"] == "OQ-005", "IIP-0002 metadata mismatch")
    _check("Marissa Famularo (ORCID 0009-0001-8032-2215, GitHub @MarissaFamularo)" in iip2["credit"], "IIP-0002 credit mismatch")
    _check(
        "Neither cited paper calibrates the chosen timing values, oxygen-transfer coefficients, or their joint distribution."
        in iip2["exploratory_engineering_scope_note"],
        "IIP-0002 missing required calibration disclaimer",
    )
    _check(
        abs(iip2["timing_axis_specification"]["realized_minus_target_midpoint_lag_hours"] - 2.9426) < 0.01,
        "IIP-0002 realized-minus-target lag != ~2.9426 h",
    )
    iip2_rows = iip2["joint_grid_rows"]
    _check(len(iip2_rows) == 74, f"Expected 74 IIP-0002 rows (2 laws x 9 timings x 4 kp + 2 zero-vessel bounds), got {len(iip2_rows)}")
    flips = iip2["release_law_direction_flip_and_gate_brackets_at_kp_0p14"]
    _check(
        flips["strict_2p0mmHg_gel_only_core_gate_flip_at_kp_0p14"]["first_order_t50_36h_gel_only_core_unrounded_mmhg"] > 2.0
        and flips["strict_2p0mmHg_gel_only_core_gate_flip_at_kp_0p14"]["first_order_t50_40h_gel_only_core_unrounded_mmhg"] < 2.0,
        "First-order strict 2.0 mmHg core gate must flip between 36 h and 40 h at k_p = 0.14",
    )
    _check(
        flips["slab_15mmHg_gate_flip_at_kp_0p14"]["first_order_t50_56h_min_slab_unrounded_mmhg"] > 15.0
        and flips["slab_15mmHg_gate_flip_at_kp_0p14"]["first_order_t50_60h_min_slab_unrounded_mmhg"] < 15.0,
        "First-order 15.0 mmHg slab gate must flip between 56 h and 60 h at k_p = 0.14",
    )
    ref_36_00944 = next(
        r
        for r in iip2["near_threshold_grid_refinement_audit_nz_51_101_201"]
        if r["release_kinetics_law"] == "first_order"
        and r["target_curve_midpoint_t50_h"] == 36.0
        and abs(r["coupled_internal_vascular_exchange_kp_s_inv"] - 0.0944) < 1e-6
    )
    _check(
        abs(ref_36_00944["mesh_refinement_nz_51_101_201"]["nz_51"]["min_slab_po2_unrounded_mmhg"] - 14.962351) < 1e-4
        and abs(ref_36_00944["mesh_refinement_nz_51_101_201"]["nz_101"]["min_slab_po2_unrounded_mmhg"] - 14.958606) < 1e-4
        and abs(ref_36_00944["mesh_refinement_nz_51_101_201"]["nz_201"]["min_slab_po2_unrounded_mmhg"] - 14.957812) < 1e-4
        and ref_36_00944["gate_stability_across_51_101_201"] is True,
        "36 h / k_p=0.0944 near-threshold refinement check mismatch",
    )
    zv_fo = iip2["prior_repository_corner_tests_and_zero_vessel_bound_context"]["exact_zero_vessel_bound_controls"]["first_order"][
        "exact_zero_internal_vascular_exchange_vperf_identically_zero"
    ]
    zv_sc = iip2["prior_repository_corner_tests_and_zero_vessel_bound_context"]["exact_zero_vessel_bound_controls"]["shrinking_core"][
        "exact_zero_internal_vascular_exchange_vperf_identically_zero"
    ]
    _check(
        zv_fo["zero_internal_vascular_state_enforced"] is True
        and abs(zv_fo["min_slab_po2_unrounded_mmhg"] - 0.111602) < 1e-4
        and abs(zv_sc["min_slab_po2_unrounded_mmhg"] - 0.071825) < 1e-4
        and zv_fo["gel_only_120um_core_po2_unrounded_mmhg"] < 1e-5
        and zv_sc["gel_only_120um_core_po2_unrounded_mmhg"] < 1e-5,
        "Exact zero-vessel bound check failed",
    )
    _check(
        iip2["open_question_title"] == "Does the oxygen supply hold up under pessimistic assumptions?",
        "IIP-0002 OQ-005 title mismatch",
    )
    iip2_json = json.dumps(iip2)
    _check(
        "PMID:14988090" in iip2_json and "PMID:31719926" in iip2_json and "31542725" not in iip2_json,
        "IIP-0002 provenance PMIDs mismatch",
    )
    req_txt = (REPO_ROOT / "requirements.txt").read_text(encoding="utf-8")
    _check("numpy>=2.0" in req_txt and "numpy>=1.24" not in req_txt, "requirements.txt must require numpy>=2.0 for np.trapezoid")

    # 2. Verify Layer 4 Part B 24h paracrine & metabolic ODE (AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)
    ode_part = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]
    _check("PMID:39169271" in ode_part["counter_evidence_and_scope_note_sol_v80_13"], "Missing PMID:39169271 counter-evidence")
    _check(len(ode_part["paracrine_submechanism_oat_ablation_table_sol_v80_r2_11"]) >= 5, "OAT table len < 5")
    _check(ode_part["low_glucose_challenge_below_4p2mM_oat1_verification_ultracode_v80_15"]["standard_protocol_oat1_exact_match_when_above_4p2mM"] is True, "OAT_1 standard match failed")
    _check(ode_part["low_glucose_challenge_below_4p2mM_oat1_verification_ultracode_v80_15"]["low_glucose_nadir_protection_delta_mM"] > 0.0, "OAT_1 low-G protection <= 0")
    _check(ode_part["zero_sst_and_common_constants_qualification_ultracode_v80_10"]["zero_sst_collapse_exact_match_to_beta_only"] is True, "Zero-SST collapse check failed")
    _check("non_overlapping_phase_volume_and_cad_ledger_7cm3_ultracode_v80_18" in alg, "Missing phase volume ledger")
    _check(len(pde_part["shrinking_core_coupled_bvp_safety_comparison_ultracode_v80_09"]["coupled_shrinking_core_regimes"]) == 4, "Shrinking-core regimes != 4")
    grafts = {g["graft_id"]: g for g in ode_part["graft_architectures_compared"]}
    tri = grafts["OPENISLET_60_30_10_TRIHORMONAL_G6PC2_WT"]
    mono = grafts["BETA_ONLY_MONOCULTURE_NO_ALPHA_DELTA"]
    g6ko = grafts["G6PC2_KO_DISRUPTED_GLYCOLYTIC_RHEOSTAT_GY002"]

    _check(tri["passes_clinical_glycemic_safety_targets"] is True, "Tri-hormonal graft should pass glycemic safety")
    _check(tri["time_in_tight_range_3p9_to_10p0_mM_pct"] >= 85.0, "Tri-hormonal TIR < 85%")
    _check(tri["time_severe_hypoglycemia_lt_3p0_mM_54_mg_dL_pct"] == 0.0, "Tri-hormonal severe hypo != 0%")
    _check(mono["passes_clinical_glycemic_safety_targets"] is False, "Beta-only graft should fail glycemic safety")
    _check(mono["min_glucose_mM"] < 3.9 and mono["time_below_range_lt_3p9_mM_70_mg_dL_pct"] > 4.0, "Beta-only hypo check failed")
    _check(g6ko["passes_clinical_glycemic_safety_targets"] is False, "G6PC2 KO graft should fail glycemic safety")
    _check(g6ko["min_glucose_mM"] < 3.0, "G6PC2 KO min glucose >= 3.0 mM")
    _check(
        tri["peripheral_insulin_action_clearance_p2_min_inv"]
        == mono["peripheral_insulin_action_clearance_p2_min_inv"]
        == g6ko["peripheral_insulin_action_clearance_p2_min_inv"],
        "Peripheral insulin clearance p2 differs across grafts",
    )
    _check(
        tri["insulin_secretion_shutoff_half_time_min"]
        < mono["insulin_secretion_shutoff_half_time_min"]
        < g6ko["insulin_secretion_shutoff_half_time_min"],
        "Shutoff half-time ordering failed",
    )

    for curve in ode_part["graft_architectures_compared"]:
        _check(curve["pde_coupling_applied"] is False, "Reference comparison falsely claims PDE coupling")
        _check(curve["source_pde_regime_id"] is None, "Reference comparison inherited an unrelated PDE")
        _check("viable_endocrine_fraction_coupled" not in curve, "Reference functional dose mislabeled coupled")
    for row in ode_part["matched_pde_ode_qualification_scenarios"]:
        _check(row["source_pde_regime_id"] == row["pde_regime_id"], "Coupled output lacks its exact PDE identifier")
        _check(row["pde_regime_id"] in regimes, "Coupled output refers to nonexistent PDE")
        _check(not row["coupled_complete_product_pass"] or regimes[row["pde_regime_id"]]["passes_all_oxygen_peroxide_and_ph_gates"], "Failed PDE produced full-product PASS")

    # 3. Verify Layer 5 Part A CRISPR gRNA compiler & live dbSNP/ClinVar records
    part_a = l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]
    guides = part_a["guides"]
    _check(len(guides) == 5, "Expected 5 CRISPR guides")
    _check(l5["summary"]["coding_off_targets_0_to_2_mismatches"] is None, "Unknown counts cannot sum to zero")
    for g in guides:
        _check(g["grch38_reference_exact_match_verified"] is True, f"Guide {g['guide_id']} GRCh38 match failed")
        _check(g["grch38_0_1_2_mismatch_coding_off_targets"] is None, f"Guide {g['guide_id']} unperformed genome-wide count must be unknown")
        _check(g["genomewide_off_target_enumeration_status"] == "NOT_ENUMERATED", "Unperformed scan has incorrect status")
        _check(g["cfd_mit_off_target_specificity_score"] >= 91.0, f"Guide {g['guide_id']} CFD < 91.0")
        _check(g["intact_target_present_in_any_hdr_donor_plasmid"] is False, f"Guide {g['guide_id']} present in HDR donor")
        _check(g["intact_target_present_in_post_hdr_edited_allele_either_strand"] is False, f"Guide {g['guide_id']} present in edited allele")
        _check(g["hdr_full_protospacer_disrupted_post_knockin"] is True, f"Guide {g['guide_id']} full protospacer not disrupted")
        _check(
            g["hdr_pam_and_seed_disrupted_post_knockin"]
            is (not g["raw_pam_and_seed_subsequence_retained_on_either_post_edit_junction"]),
            f"Guide {g['guide_id']} hdr_pam_and_seed_disrupted_post_knockin mismatch with junction audit",
        )
    g1 = next(g for g in guides if g["guide_id"].startswith("OI-GRNA-01"))
    g2 = next(g for g in guides if g["guide_id"].startswith("OI-GRNA-02"))
    g3 = next(g for g in guides if g["guide_id"].startswith("OI-GRNA-03"))
    g5 = next(g for g in guides if g["guide_id"].startswith("OI-GRNA-05"))
    _check(abs(g1["sugimoto_1995_rna_dna_hybrid_dg37_kcal_mol"] - (-36.09)) < 1e-2, "g1 Sugimoto dG37 mismatch")
    _check(abs(g3["sugimoto_1995_rna_dna_hybrid_dg37_kcal_mol"] - (-23.09)) < 1e-2, "g3 Sugimoto dG37 mismatch")
    _check(g2["target_locus_grch38"] == "chr16:10898963-10898989 (-)", "g2 CIITA Exon 5 locus mismatch")
    _check(len(g2["targeted_transcripts_pan_isoform"]) == 3, "g2 targeted transcripts != 3")
    _check("PEmax" in g5["nuclease"] and "ABE8e" not in g5["nuclease"], "g5 nuclease should be PEmax")
    peg = g5["epegRNA_extension_specification"]
    _check(peg["primer_binding_site_pbs_13nt_5to3"] == "CACGGTGAGAGCC", "g5 PBS mismatch")
    _check(peg["rt_template_rtt_12nt_5to3"] == "CCAGAGCGCCCC", "g5 RTT mismatch")
    _check(peg["rtt_written_plus_strand_10352438_10352449"] == "GGGGCGCTCTGG", "g5 plus strand mismatch")

    _check(len(part_a["live_ncbi_dbsnp_polymorphism_audit"]) == 4, "dbSNP audit count != 4")
    _check("rs34536443_TYK2_P1104A" in part_a["live_ncbi_dbsnp_polymorphism_audit"], "Missing rs34536443")

    # 4. Verify Layer 5 Part B Personalized Patient Compiler & 10,000-patient simulation
    part_b = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]
    _check(part_b["total_virtual_patients_simulated"] == 10000, "Total virtual patients != 10000")
    _check(len(part_b["ancestry_stratified_cohorts"]) == 5, "Ancestry cohorts != 5")
    _check(part_b["global_5yr_survival_unedited_allogeneic_pct"] == 0.0, "Unedited survival != 0.0%")
    _check(part_b["global_5yr_survival_hla_e_only_pct"] < 55.0, "HLA-E only survival >= 55%")
    _check(
        part_b["global_5yr_survival_ablated_compiler_control_pct"]
        < part_b["global_5yr_survival_standardized_openislet_plan_a_pct"]
        < part_b["global_5yr_survival_personalized_dna_compiled_openislet_pct"],
        "Survival ordering failed across Ablated < Plan A < Personalized",
    )
    _check(l5["summary"]["analytic_max_ancestry_disparity_gap_personalized_pct"] == part_b["analytic_max_ancestry_disparity_gap_personalized_pct_v80_065"], "Primary analytic ancestry summary drift")
    _check(l5["summary"]["sampled_max_ancestry_disparity_gap_personalized_pct"] == part_b["sampled_max_ancestry_disparity_gap_personalized_pct"], "Secondary sampled ancestry summary drift")
    _check(part_b["sampled_max_ancestry_disparity_gap_personalized_pct"] <= 2.5, "Sampled ancestry disparity > 2.5%")
    _check(len(part_b["representative_patient_compilations"]) == 4, "Representative patients != 4")
    pat_plan_b = next(p for p in part_b["representative_patient_compilations"] if p["patient_id"].startswith("PT-002"))
    _check("Plan B" in pat_plan_b["plan_architecture_prescribed"], "PT-002 should prescribe Plan B")
    _check(all(not g_id.startswith("OI-GRNA-01") for g_id in pat_plan_b["compiled_grna_panel"]), "Plan B must not use OI-GRNA-01")
    _check("PLAN_B_HLAA_EXON2_KO_GUIDE_INVESTIGATIONAL" in pat_plan_b["compiled_grna_panel"], "Plan B missing HLA-A guide")
    _check("PLAN_B_HLAB_EXON2_KO_GUIDE_INVESTIGATIONAL" in pat_plan_b["compiled_grna_panel"], "Plan B missing HLA-B guide")
    _check(not pat_plan_b["locus1_construct_prescribed"].startswith("OI-L2-DONOR-1"), "Plan B must not prescribe OI-L2-DONOR-1*")
    _check(len(part_b["correlated_frailty_copula_sensitivity_sol_v80_r2_07"]) == 4, "Copula sensitivity len != 4")

    # AUDIT-REF & AUDIT-REF Negative-control input validation probes on Layer 5 compiler
    rejected_invalid_sirpa = False
    try:
        l5_mod.compile_personalized_patient_blueprint({
            "patient_id": "TEST-BAD-SIRPA",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.35,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.2,
            "ifn1_ish_genetic_risk_percentile": 50.0,
            "sirpa_genotype": "NOT_A_REAL_ALLELE",
        })
    except ValueError:
        rejected_invalid_sirpa = True
    _check(rejected_invalid_sirpa, "AUDIT-REF FAILED: invalid SIRPA allele was not rejected!")

    rejected_nan_pct = False
    try:
        l5_mod.compile_personalized_patient_blueprint({
            "patient_id": "TEST-NAN-PCT",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": float("nan"),
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.2,
            "ifn1_ish_genetic_risk_percentile": 50.0,
            "sirpa_genotype": "v1/v1",
        })
    except ValueError:
        rejected_nan_pct = True
    _check(rejected_nan_pct, "AUDIT-REF FAILED: NaN percentage was not rejected!")

    rejected_bool_as_float = False
    try:
        l5_mod.compile_personalized_patient_blueprint({
            "patient_id": "TEST-BOOL-FLOAT",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": True,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.2,
            "ifn1_ish_genetic_risk_percentile": 50.0,
            "sirpa_genotype": "v1/v1",
        })
    except ValueError:
        rejected_bool_as_float = True
    _check(rejected_bool_as_float, "AUDIT-REF FAILED: boolean True passed as float fraction was not rejected!")

    rejected_bad_regimen = False
    try:
        l5_mod.compile_personalized_patient_blueprint({
            "patient_id": "TEST-BAD-REGIMEN",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.25,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.2,
            "ifn1_ish_genetic_risk_percentile": 50.0,
            "sirpa_genotype": "v1/v1",
            "carries_rare_snp_in_standard_grna": "false",
        })
    except ValueError:
        rejected_bad_regimen = True
    _check(rejected_bad_regimen, "AUDIT-REF FAILED: string 'false' for carries_rare_snp_in_standard_grna was not rejected!")

    # 5. Verify Closed-Loop Wet-Lab Conjugate Bayesian Updater & AUDIT-REF + AUDIT-REF + AUDIT-REF fail-closed probes
    _check(calib["calibration_mode"] == "SYNTHETIC_DRY_RUN_VERIFICATION_DO_NOT_CLAIM_AS_WETLAB_DATA", "Calibration mode mismatch")
    _check(len(calib["conjugate_posterior_updates"]) >= 4, "Posterior updates < 4")

    # AUDIT-REF: Independent asymmetric conjugate-normal golden oracle (prior=0.0, sd=1.0; obs=1.0, se=2.0 -> posterior_mean=0.2, NOT 0.8)
    asym_golden = wetlab.conjugate_normal_update(0.0, 1.0, 1.0, 2.0)
    _check(
        abs(asym_golden["posterior_mean_exact"] - 0.2) < 1e-9 and abs(asym_golden["posterior_mean"] - 0.2) < 1e-4,
        f"Asymmetric conjugate-normal golden check failed: expected 0.2, got {asym_golden['posterior_mean_exact']} (swapped prior/obs precision weights!)",
    )

    # AUDIT-REF: Two sequential near-floor observations (obs_se = 0.001 -> posterior_sd < 0.001) must succeed
    near_floor_1 = wetlab.run_calibration({
        "batch_id": "BATCH-FLOOR-01",
        "obs_stage7_yield_ratio": 0.96,
        "obs_stage7_yield_se": 0.001,
        "obs_cd47_fold": "UNAVAILABLE",
        "obs_serpinb9_log10_kass": "UNAVAILABLE",
        "obs_icasp9_log10_kill": "UNAVAILABLE",
    }, False)
    near_floor_2 = wetlab.run_calibration({
        "batch_id": "BATCH-FLOOR-02",
        "obs_stage7_yield_ratio": 0.95,
        "obs_stage7_yield_se": 0.001,
        "obs_cd47_fold": "UNAVAILABLE",
        "obs_serpinb9_log10_kass": "UNAVAILABLE",
        "obs_icasp9_log10_kill": "UNAVAILABLE",
    }, False, prior_posterior_report=near_floor_1)
    _check(near_floor_2["ingested_batch_ids"] == ["BATCH-FLOOR-01", "BATCH-FLOOR-02"], "Sequential near-floor update failed")

    # AUDIT-REF: Tampered posterior_mean_exact with unchanged conjugate_posterior_updates_sha256 must reject
    tampered_prior = json.loads(json.dumps(near_floor_1))
    next(iter(tampered_prior["conjugate_posterior_updates"].values()))["posterior_mean_exact"] = 0.123456
    rejected_tampered_prior = False
    try:
        wetlab.run_calibration({
            "batch_id": "BATCH-AFTER-TAMPER",
            "obs_stage7_yield_ratio": 0.95,
            "obs_stage7_yield_se": 0.03,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        }, False, prior_posterior_report=tampered_prior)
    except ValueError:
        rejected_tampered_prior = True
    _check(rejected_tampered_prior, "AUDIT-REF FAILED: tampered posterior_mean_exact was not rejected!")

    # AUDIT-REF / AUDIT-REF: Tampered endpoint status with recomputed SHA-256 must reject
    import hashlib
    tampered_status_prior = json.loads(json.dumps(near_floor_1))
    ep0 = next(iter(tampered_status_prior["conjugate_posterior_updates"].values()))
    ep0["status"] = "TAMPERED_BYPASS_STATUS"
    ep0["posterior_mean_exact"] = 0.123456
    tampered_status_prior["conjugate_posterior_updates_sha256"] = hashlib.sha256(
        json.dumps(tampered_status_prior["conjugate_posterior_updates"], sort_keys=True, allow_nan=False).encode("utf-8")
    ).hexdigest()
    rejected_tampered_status = False
    try:
        wetlab.run_calibration({
            "batch_id": "BATCH-AFTER-STATUS-TAMPER",
            "obs_stage7_yield_ratio": 0.95,
            "obs_stage7_yield_se": 0.03,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        }, False, prior_posterior_report=tampered_status_prior)
    except ValueError:
        rejected_tampered_status = True
    _check(rejected_tampered_status, "AUDIT-REF FAILED: unrecognized endpoint status with recomputed hash was not rejected!")

    rejected_empty_empirical = False
    try:
        wetlab.run_calibration({}, False)
    except ValueError:
        rejected_empty_empirical = True
    _check(rejected_empty_empirical, "AUDIT-REF FAILED: empty empirical wet-lab payload was not rejected!")

    rejected_out_of_range = False
    try:
        wetlab.run_calibration({
            "batch_id": "BATCH-TEST-01",
            "obs_stage7_yield_ratio": -0.25,
            "obs_stage7_yield_se": 0.04,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        }, False)
    except ValueError:
        rejected_out_of_range = True
    _check(rejected_out_of_range, "AUDIT-REF FAILED: negative yield ratio was not rejected!")

    seq_b1 = wetlab.run_calibration({
        "batch_id": "BATCH-SEQ-01",
        "obs_stage7_yield_ratio": 0.99,
        "obs_stage7_yield_se": 0.03,
        "obs_cd47_fold": "UNAVAILABLE",
        "obs_serpinb9_log10_kass": "UNAVAILABLE",
        "obs_icasp9_log10_kill": "UNAVAILABLE",
    }, False)
    seq_b2 = wetlab.run_calibration({
        "batch_id": "BATCH-SEQ-02",
        "obs_stage7_yield_ratio": 0.98,
        "obs_stage7_yield_se": 0.03,
        "obs_cd47_fold": "UNAVAILABLE",
        "obs_serpinb9_log10_kass": "UNAVAILABLE",
        "obs_icasp9_log10_kill": "UNAVAILABLE",
    }, False, prior_posterior_report=seq_b1)
    _check(seq_b2["ingested_batch_ids"] == ["BATCH-SEQ-01", "BATCH-SEQ-02"], "Sequential batch IDs mismatch")
    rejected_dup_batch = False
    try:
        wetlab.run_calibration({
            "batch_id": "BATCH-SEQ-01",
            "obs_stage7_yield_ratio": 0.98,
            "obs_stage7_yield_se": 0.03,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        }, False, prior_posterior_report=seq_b2)
    except ValueError:
        rejected_dup_batch = True
    _check(rejected_dup_batch, "AUDIT-REF FAILED: duplicate batch_id replay was not rejected!")

    # 6. Verify generate_release_manifest.py --export-public-bundle (AUDIT-REF) and isolated CLI --out flag (AUDIT-REF)
    orig_calib_bytes = WETLAB_CALIB.read_bytes()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        bundle_out = tmp_path / "public_bundle"
        manifest_mod.export_public_bundle(bundle_out)
        _check((bundle_out / "manifest.json").exists(), "Missing exported manifest.json")
        _check((bundle_out / "data" / "extracted" / "cross_study_discovery_graph.json").exists(), "Missing exported discovery graph")

        emp_in = tmp_path / "emp_batch.json"
        emp_out = tmp_path / "emp_posterior_out.json"
        emp_in.write_text(json.dumps({
            "batch_id": "BATCH-CLI-TMP-01",
            "assay_id": "ASSAY-DE1-TMP",
            "obs_stage7_yield_ratio": 0.97,
            "obs_stage7_yield_se": 0.03,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        }), encoding="utf-8")
        proc = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "ingest_wetlab_results.py"), "--input-json", str(emp_in), "--output", str(emp_out)],
            capture_output=True,
            text=True,
            check=False,
        )
        _check(proc.returncode == 0, f"ingest_wetlab_results.py CLI failed with returncode={proc.returncode}: {proc.stderr}")
        _check(emp_out.exists(), "CLI --output file was not created")
        _check(WETLAB_CALIB.read_bytes() == orig_calib_bytes, "AUDIT-REF FAILED: committed calibration artifact was modified!")

    print(
        "layer4_and_layer5_planes gate: PASS (336h PDE + isolated peroxide control, 24h paracrine ODE + OAT ablation, "
        "5/5 strand-exact CRISPR gRNAs + epegRNA round-trip, 4 live dbSNP/ClinVar records, 10,000-patient compiler simulation, "
        "sequential Bayesian updater, and atomic public bundle export verified)"
    )


run_layer4_and_layer5_suite = main


def test_layer4_and_layer5_planes() -> None:
    main()


if __name__ == "__main__":
    main()
