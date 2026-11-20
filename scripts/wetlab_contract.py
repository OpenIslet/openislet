"""One measurement contract for Python ingestion and generated MCP metadata."""
from __future__ import annotations

import ast
import inspect
from datetime import datetime, timezone
from decimal import Context, Decimal, localcontext
from functools import lru_cache
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import re

CONTRACT_RELATIVE_PATH = Path('data/protocols/wetlab_endpoint_contract_v0.8.0.json')
IIP0031_MIN_UNTREATED_CFU_PLATED_CELLS = 2
STAGE0_DECISION_TREE_RELATIVE_PATH = Path('data/protocols/stage0_decision_tree_v1.json')
DECIMAL = re.compile(r'^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$')
STAGE0_MAX_INDEPENDENT_LOT_REPLICATES = 10_000.0

STAGE0_DECISION_TREE_ENDPOINT_SPECS = [
    ('stage0_gate0b_t1d_pbmc_4digit_hlac_cd8_elispot_sfw_per_1e6', 'obs_gate0b_hlac_cd8_elispot_sfw', 'obs_gate0b_hlac_cd8_elispot_se', 0.0, 5000.0, 'IFNg_granzymeB_SFW_per_1e6_CD8_T_cells'),
    ('stage0_gate0b_hlac_sct_ciita_surface_molecules_per_cell', 'obs_gate0b_hlac_sct_molecules_per_cell', 'obs_gate0b_hlac_sct_molecules_se', 0.0, 1000000.0, 'QuantiBRITE_molecules_per_cell'),
    ('stage0_gate0c_hour0_resting_serpinb9_cytosolic_uM', 'obs_gate0c_hour0_serpinb9_uM', 'obs_gate0c_hour0_serpinb9_se', 0.0, 20.0, 'cytosolic_micromolar_uM'),
    ('stage0_gate0c_hour0_resting_pdl1_surface_molecules_per_cell', 'obs_gate0c_hour0_pdl1_molecules_per_cell', 'obs_gate0c_hour0_pdl1_molecules_se', 0.0, 1000000.0, 'QuantiBRITE_molecules_per_cell'),
    ('stage0_gate0c_ifng_tnfa_induction_t50_hours', 'obs_gate0c_induction_t50_h', 'obs_gate0c_induction_t50_se', 0.0, 72.0, 'hours_to_half_maximal_protein_plateau'),
    ('stage0_gate0c_2h_ctl_nk_lysis_pct_regime2', 'obs_gate0c_2h_lysis_pct', 'obs_gate0c_2h_lysis_se', 0.0, 100.0, 'specific_lysis_percentage_at_2h'),
    ('stage0_gate0e_hc10_surface_free_heavy_chain_mfi_over_isotype', 'obs_gate0e_hc10_surface_mfi_ratio', 'obs_gate0e_hc10_surface_mfi_se', 0.0, 100.0, 'MFI_fold_over_isotype_control'),
    ('stage0_gate0e_hc10_intracellular_retained_heavy_chain_fold_over_wt', 'obs_gate0e_hc10_intracellular_fold', 'obs_gate0e_hc10_intracellular_se', 0.0, 100.0, 'fold_over_unedited_wt_ifng'),
    ('stage0_gate0e_ifng_16p7mM_chop_ddit3_log2fc', 'obs_gate0e_chop_log2fc', 'obs_gate0e_chop_log2fc_se', -10.0, 15.0, 'log2_fold_change_vs_resting'),
    ('stage0_gate0e_ifng_16p7mM_gsis_si_retention_ratio', 'obs_gate0e_gsis_retention_ratio', 'obs_gate0e_gsis_retention_se', 0.0, 3.0, 'ratio_stressed_to_unstressed_GSIS_SI'),
    ('stage0_gate0m_perifusion_2p8mM_shutoff_t90_min', 'obs_gate0m_shutoff_t90_min', 'obs_gate0m_shutoff_t90_se', 0.0, 60.0, 'minutes_to_90pct_secretion_shutoff'),
    ('stage0_gate0f_d28_explant_collagen_capsule_thickness_um', 'obs_gate0f_capsule_thickness_um', 'obs_gate0f_capsule_thickness_se', 0.0, 500.0, 'micrometers_dense_collagen_capsule'),
    ('stage0_dxs_aavs1_locus4_joint_icasp9_cd47_retention_fraction', 'obs_dxs_locus4_joint_retention_fraction', 'obs_dxs_locus4_joint_retention_se', 0.0, 1.0, 'fraction_retaining_both_CD47_and_AP1903_kill'),
    # Canonical stage0_decision_tree_v1.json meanKey/seKey bindings + S9 hard safety stop endpoints (AUDIT-REF):
    ('locus3_serpinb9_hour0_cytosolic_uM', 'obs_serpinb9_hour0_uM', 'obs_serpinb9_hour0_se', 0.0, 20.0, 'hour0_unstimulated_cytosolic_serpinb9_uM_closed_interval_[0.0, 20.0]'),
    ('locus3_pdl1_hour0_surface_molecules_per_cell', 'obs_pdl1_hour0_mols', 'obs_pdl1_hour0_se', 0.0, 500000.0, 'hour0_unstimulated_surface_pdl1_molecules_per_cell_closed_interval_[0.0, 500000.0]'),
    ('regime_ii_uninduced_2h_specific_lysis_pct', 'obs_regimeII_lysis_2h_pct', 'obs_regimeII_lysis_2h_se', 0.0, 100.0, 'regime_ii_2h_specific_lysis_percentage_closed_interval_[0.0, 100.0]'),
    ('b2m_null_orphan_hc10_surface_fold_over_wt', 'obs_hc10_surface_fold_wt', 'obs_hc10_surface_fold_se', 0.0, 50.0, 'hc10_free_heavy_chain_surface_mfi_fold_over_wt_closed_interval_[0.0, 50.0]'),
    ('b2m_null_chop_ddit3_coincidence_fold_over_basal', 'obs_chop_coincidence_fold', 'obs_chop_coincidence_se', 0.0, 50.0, 'ddit3_chop_coincidence_arm_fold_over_basal_closed_interval_[0.0, 50.0]'),
    ('chronic_day28_90_capsule_thickness_um', 'obs_capsule_dz_um', 'obs_capsule_dz_se', 0.0, 500.0, 'physical_fibrotic_capsule_thickness_um_closed_interval_[0.0, 500.0]'),
    ('chronic_day28_90_capsule_o2_diffusivity_cm2_s', 'obs_capsule_d_o2_cm2_s', 'obs_capsule_d_o2_se', 1e-7, 5e-5, 'capsule_oxygen_diffusivity_cm2_s_closed_interval_[1e-7, 5e-5]'),
    ('stage7_basal_2p8mM_insulin_secretion_fraction_of_peak', 'obs_basal_insulin_frac_peak', 'obs_basal_insulin_frac_se', 0.0, 1.0, 'basal_2p8mM_secretion_fraction_of_16p7mM_peak_closed_interval_[0.0, 1.0]'),
    ('stage7_step_ramp_glucose_secretion_threshold_mM', 'obs_glucose_threshold_mM', 'obs_glucose_threshold_se', 1.0, 15.0, 'step_ramp_first_sustained_secretion_threshold_mM_closed_interval_[1.0, 15.0]'),
    ('b_prime_foreign_group_matched_cd8_specific_lysis_pct', 'obs_foreign_hlac_cd8_lysis_pct', 'obs_foreign_hlac_cd8_lysis_se', 0.0, 100.0, 'foreign_4digit_group_matched_cd8_specific_lysis_pct_closed_interval_[0.0, 100.0]'),
    # Gate_P05_kinetics (IIP-0010 / DEC-088): qualified per distinct CaO2:PDMS wt%/vol% recipe, PSD, catalase coating, and as-implanted state (never inherited across formulations):
    ('cao2_hydrolysis_kinetic_law_code', 'obs_cao2_kinetic_law', 'obs_cao2_kinetic_law_se', 0.0, 2.0, 'kinetic_law_code_0_first_order_1_shrinking_core_2_ginstling_brounshtein_closed_interval_[0.0, 2.0]'),
    ('stage7_lot_ocr_ratio_over_nominal_34_nmol_cm3_s', 'obs_lot_ocr_ratio', 'obs_lot_ocr_ratio_se', 0.2, 3.0, 'lot_ocr_vmax_ratio_over_nominal_34_nmol_cm3_s_closed_interval_[0.2, 3.0]'),
    ('safety_kill_s9_icasp9_log10_kill', 'obs_icasp9_log10_kill', 'obs_icasp9_log10_kill_se', 0.0, 10.0, 'log10_clonal_depletion_kill_depth'),
    ('safety_kill_s9_tp53_ppm1d_bcor_vaf', 'obs_tp53_ppm1d_bcor_vaf', 'obs_tp53_ppm1d_bcor_vaf_se', 0.0, 1.0, 'variant_allele_fraction_0_to_1'),
    ('safety_kill_s9_castseq_translocation_freq', 'obs_castseq_translocation_freq', 'obs_castseq_translocation_se', 0.0, 1.0, 'castseq_translocation_frequency_0_to_1'),
    # Additional unified JSON switch endpoints (AUDIT-REF / SAF-01):
    ('locus1_cd47_hour0_surface_molecules_per_cell', 'obs_cd47_hour0_molecules', 'obs_cd47_hour0_molecules_se', 0.0, 1000000.0, 'QuantiBRITE_CD47_molecules_per_cell'),
    ('locus1_macrophage_phagocytosis_index_pct', 'obs_phagocytosis_index_pct', 'obs_phagocytosis_index_se', 0.0, 100.0, 'primary_human_macrophage_phagocytosis_percentage'),
    ('locus1_hlae_hour0_surface_molecules_per_cell', 'obs_hlae_hour0_molecules', 'obs_hlae_hour0_molecules_se', 0.0, 1000000.0, 'QuantiBRITE_HLAE_molecules_per_cell'),
    ('locus1_nkg2a_positive_nk_lysis_pct', 'obs_nkg2a_lysis_pct', 'obs_nkg2a_lysis_se', 0.0, 100.0, 'nkg2a_positive_nk_specific_lysis_percentage'),
    ('locus3_serpinb9_stoichiometry_of_inhibition', 'obs_serpinb9_si', 'obs_serpinb9_si_se', 1.0, 10.0, 'moles_serpinb9_per_mole_gzmb_inactivated'),
    ('locus3_serpinb9_intracellular_half_life_h', 'obs_serpinb9_half_life_h', 'obs_serpinb9_half_life_se', 0.1, 72.0, 'cycloheximide_chase_half_life_hours'),
    ('locus3_serpinb9_gzmb_kass_M_inv_s_inv', 'obs_serpinb9_kass_M_s', 'obs_serpinb9_kass_M_s_se', 1e3, 1e9, 'second_order_association_rate_M_inv_s_inv'),
    ('plan_a_star_reference_cd8_specific_lysis_pct', 'obs_plan_a_star_cd8_lysis_pct', 'obs_plan_a_star_cd8_lysis_se', 0.0, 100.0, 'plan_a_star_comparator_cd8_lysis_percentage'),
    ('stage7_immaturity_corner_min_glucose_mM', 'obs_immaturity_corner_min_glucose_mM', 'obs_immaturity_corner_min_glucose_se', 0.5, 15.0, 'minimum_plasma_glucose_mM_under_exercise_challenge'),
    ('chronic_day90_aged_hydrogel_o2_diffusivity_cm2_s', 'obs_d_o2_gel_aged_cm2_s', 'obs_d_o2_gel_aged_se', 1e-7, 5e-5, 'aged_lamina_hydrogel_o2_diffusivity_cm2_s'),
    ('dxs_locus4_cd47low_icasp9low_double_loss_fraction', 'obs_dxs_cd47low_icasp9low_fraction', 'obs_dxs_cd47low_icasp9low_se', 0.0, 1.0, 'fraction_double_negative_cd47low_icasp9low'),
    ('dxs_locus4_ap1903_surviving_cd47high_escape_fraction', 'obs_dxs_ap1903_cd47high_escape_fraction', 'obs_dxs_ap1903_cd47high_escape_se', 0.0, 1.0, 'clonogenic_cd47high_survivor_fraction_post_ap1903'),
    # TYK2 x inducible-promoter epistasis & 3x2 satellite shield induction endpoints (AUDIT-REF; DEC-040 & AUDIT-REF; DEC-041):
    ('stage0_tyk2_x1_x4_hlae_epistasis_beta', 'obs_beta_x1_x4_epistasis', 'obs_beta_x1_x4_epistasis_se', -50.0, 50.0, 'standardized_or_pp_interaction_coefficient_X1_X4_half_effect'),
    ('stage0_tyk2_x3_x4_locus3_epistasis_beta', 'obs_beta_x3_x4_epistasis', 'obs_beta_x3_x4_epistasis_se', -50.0, 50.0, 'standardized_or_pp_interaction_coefficient_X3_X4_half_effect'),
    ('stage0_tyk2_p1104a_3x2_satellite_shield_induction_pct', 'obs_tyk2_p1104a_3x2_satellite_shield_induction_pct', 'obs_tyk2_p1104a_3x2_satellite_shield_induction_se', 0.0, 200.0, 'ifng_tnfa_shield_induction_percentage_of_tyk2_wt_plateau'),
    # Sorted NKG2C+ adaptive NK & KIR+NKG2A- licensed NK lysis + B'-SCT direct CD8 lysis + rare-event count depth endpoints (AUDIT-REF & S2-01; DEC-041):
    ('locus1_nkg2c_positive_adaptive_nk_lysis_pct', 'obs_nkg2c_adaptive_nk_lysis_pct', 'obs_nkg2c_adaptive_nk_lysis_se', 0.0, 100.0, 'cmv_imprinted_nkg2c_pos_nkg2a_neg_cd57_pos_nk_specific_lysis_percentage'),
    ('locus1_kir_licensed_nkg2a_negative_nk_lysis_pct', 'obs_kir_nkg2a_neg_nk_lysis_pct', 'obs_kir_nkg2a_neg_nk_lysis_se', 0.0, 100.0, 'kir2dl_kir3dl1_licensed_nkg2a_neg_nk_specific_lysis_percentage'),
    ('b_prime_sct_mismatched_recipient_cd8_specific_lysis_pct', 'obs_b_prime_sct_cd8_lysis_pct', 'obs_b_prime_sct_cd8_lysis_se', 0.0, 100.0, 'hlac_mismatched_recipient_cd8_specific_lysis_against_b_prime_sct_percentage'),
    ('safety_kill_s9_tp53_ppm1d_bcor_duplex_molecules_assayed', 'obs_tp53_ppm1d_bcor_duplex_depth', 'obs_tp53_ppm1d_bcor_duplex_depth_se', 100.0, 1e9, 'unique_duplex_consensus_molecules_interrogated'),
    ('safety_kill_s9_tp53_ppm1d_bcor_mutant_molecules_k', 'obs_tp53_ppm1d_bcor_mutant_molecules_k', 'obs_tp53_ppm1d_bcor_mutant_molecules_k_se', 0.0, 1e6, 'integer_mutant_duplex_molecules_observed'),
    ('safety_kill_s9_castseq_duplex_or_amplicon_depth', 'obs_castseq_duplex_or_amplicon_depth', 'obs_castseq_duplex_or_amplicon_depth_se', 100.0, 1e9, 'unique_castseq_amplicon_or_duplex_molecules_interrogated'),
    ('safety_kill_s9_castseq_mutant_junctions_k', 'obs_castseq_mutant_junctions_k', 'obs_castseq_mutant_junctions_k_se', 0.0, 1e6, 'integer_castseq_translocation_junctions_observed'),
    ('dxs_locus4_rare_event_screened_cell_count', 'obs_dxs_screened_cell_count', 'obs_dxs_screened_cell_count_se', 100.0, 1e11, 'total_cells_interrogated_for_dxs_escape_and_retention'),
    ('dxs_locus4_clonogenic_fraction_f', 'obs_dxs_clonogenic_fraction', 'obs_dxs_clonogenic_fraction_se', 1e-6, 1.0, 'fraction_of_screened_cells_capable_of_clonogenic_colony_formation'),
    ('dxs_locus4_spike_in_detection_efficiency_eta', 'obs_dxs_spike_in_detection_efficiency', 'obs_dxs_spike_in_detection_efficiency_se', 0.01, 1.0, 'empirical_or_clopper_pearson_lower95_spike_in_recovery_efficiency'),
    ('locus1_hlae_nkg2c_over_nkg2a_kd_selectivity_ratio', 'obs_hlae_nkg2c_over_nkg2a_kd_ratio', 'obs_hlae_nkg2c_over_nkg2a_kd_ratio_se', 0.1, 100.0, 'covalent_sct_spr_kd_nkg2c_over_kd_nkg2a_ratio'),
    # Pivotal combined-stress Stage-0 coupon gate endpoints (STAGE0_COMBINED_STRESS_KILL_CHAIN_COUPON_MUSE_R4; AUDIT-REF & AUDIT-REF; DEC-043/044):
    ('combined_stress_coupon_embodiment_a_star_lysis_pct', 'obs_combined_stress_embodiment_a_star_lysis_pct', 'obs_combined_stress_embodiment_a_star_lysis_se', 0.0, 100.0, 'combined_stress_embodiment_a_star_specific_lysis_percentage'),
    ('combined_stress_coupon_comp_r4_lysis_pct', 'obs_combined_stress_comp_r4_lysis_pct', 'obs_combined_stress_comp_r4_lysis_se', 0.0, 100.0, 'combined_stress_comp_r4_specific_lysis_percentage'),
    ('combined_stress_coupon_comp_r4_gsis_si', 'obs_combined_stress_comp_r4_gsis_si', 'obs_combined_stress_comp_r4_gsis_si_se', 0.0, 30.0, 'combined_stress_comp_r4_post_stress_gsis_stimulation_index'),
    ('combined_stress_coupon_comp_r4_viable_beta_mass_pct', 'obs_combined_stress_comp_r4_viable_beta_mass_pct', 'obs_combined_stress_comp_r4_viable_beta_mass_se', 0.0, 100.0, 'combined_stress_comp_r4_post_stress_viable_beta_mass_percentage'),
    ('combined_stress_coupon_comp_r4_shutoff_t90_min', 'obs_combined_stress_comp_r4_shutoff_t90_min', 'obs_combined_stress_comp_r4_shutoff_t90_se', 0.1, 30.0, 'combined_stress_comp_r4_post_stepdown_insulin_shutoff_t90_minutes'),
    ('combined_stress_coupon_comp_r4_max_h2o2_uM', 'obs_combined_stress_comp_r4_max_h2o2_uM', 'obs_combined_stress_comp_r4_max_h2o2_se', 0.0, 1000.0, 'combined_stress_comp_r4_peak_peroxide_uM'),
    ('combined_stress_coupon_comp_r4_max_ph', 'obs_combined_stress_comp_r4_max_ph', 'obs_combined_stress_comp_r4_max_ph_se', 6.5, 8.5, 'combined_stress_comp_r4_peak_alkaline_ph'),
    # MCB clonogenic escape, challenged alpha/delta counter-regulation, and indirect APC-CD4 autoimmunity endpoints (AUDIT-REF & AUDIT-REF; DEC-044/045):
    ('safety_kill_s9_mcb_effective_clonogenic_trials', 'obs_mcb_effective_clonogenic_trials', 'obs_mcb_effective_clonogenic_trials_se', 100.0, 1e11, 'effective_clonogenic_opportunities_interrogated_in_mcb_post_ap1903_outgrowth'),
    ('safety_kill_s9_mcb_clonogenic_survivor_count_k', 'obs_mcb_clonogenic_survivor_count_k', 'obs_mcb_clonogenic_survivor_count_k_se', 0.0, 1e6, 'integer_surviving_clonogenic_colonies_observed_post_ap1903_outgrowth'),
    ('stage7_alpha_glucagon_low_glucose_fold_change', 'obs_alpha_glucagon_low_glc_fold', 'obs_alpha_glucagon_low_glc_fold_se', 0.1, 50.0, 'alpha_cell_glucagon_secretion_fold_at_1mM_plus_arg_over_11mM_glc'),
    ('stage7_alpha_counterregulatory_gain_ratio', 'obs_alpha_counterregulatory_gain_ratio', 'obs_alpha_counterregulatory_gain_ratio_se', 0.0, 5.0, 'alpha_counterregulatory_glucagon_gain_fraction_of_mature_human_islet'),
    ('stage7_ucn3_delta_paracrine_maturity_ratio', 'obs_ucn3_maturity_ratio', 'obs_ucn3_maturity_ratio_se', 0.0, 5.0, 'ucn3_and_delta_somatostatin_paracrine_feedback_fraction_of_mature_islet'),
    ('gate0m_glucagon_2mM_counterregulatory_fold_increase', 'obs_gate0m_glucagon_2mM_fold_increase', 'obs_gate0m_glucagon_2mM_fold_increase_se', 0.1, 50.0, 'stepdown_16p7_to_2p0mM_glucose_glucagon_fold_increase'),
    ('gate0d_indirect_apc_dq2_dq8_cd4_stimulation_index', 'obs_indirect_apc_dq2_dq8_cd4_stimulation_index', 'obs_indirect_apc_dq2_dq8_cd4_stimulation_index_se', 0.1, 100.0, 'recipient_apc_pulsed_stressed_graft_lysate_hip_cd4_proliferation_stimulation_index'),
    ('gate0d_indirect_apc_cd4_ifng_stimulation_index', 'obs_indirect_apc_cd4_ifng_si', 'obs_indirect_apc_cd4_ifng_si_se', 0.1, 100.0, 'recipient_apc_pulsed_stressed_graft_lysate_hip_cd4_ifng_stimulation_index'),
    ('gate0d_sa_fasl_activated_cd4_apoptosis_pct', 'obs_gate0d_activated_cd4_fasl_apoptosis_pct', 'obs_gate0d_activated_cd4_fasl_apoptosis_se', 0.0, 100.0, 'sa_fasl_cuff_activated_autoreactive_cd4_apoptosis_percentage_at_24h'),
    ('gate0d_transwell_bystander_beta_annexin_v_pct', 'obs_gate0d_transwell_beta_annexin_v_pct', 'obs_gate0d_transwell_beta_annexin_v_se', 0.0, 100.0, 'indirect_apc_cd4_transwell_bystander_beta_cell_annexin_v_apoptosis_percentage'),
    ('gate0e_exogenous_b2m_w632_surface_fold_over_isotype', 'obs_exogenous_b2m_w632_surface_fold_isotype', 'obs_exogenous_b2m_w632_surface_fold_se', 0.0, 100.0, 'ifng_plus_exogenous_b2m_3ug_ml_w632_surface_mfi_fold_over_isotype'),
    ('gate0e_exogenous_b2m_t1d_ctl_lysis_pct', 'obs_exogenous_b2m_t1d_ctl_lysis_pct', 'obs_exogenous_b2m_t1d_ctl_lysis_se', 0.0, 100.0, 'ifng_plus_exogenous_b2m_3ug_ml_primary_t1d_autoreactive_cd8_ctl_lysis_percentage'),
    ('gate0f_outer_band_vascular_perfusion_rate_s_inv', 'obs_outer_band_kperf_s_inv', 'obs_outer_band_kperf_se', 0.0, 1.0, 'chronic_explant_outer_vascular_band_perfusion_rate_constant_s_inv'),
    # Pan-lineage CD47 QPCTL/complement/vascular shield, receptor-negative NKG2C NK, 12-pulse unprimed hypoxic serial lysis, and insulin clearance endpoints (AUDIT-REF; DEC-046):
    ('locus1_cd47_cc2c6_over_b6h12_mfi_ratio', 'obs_cd47_cc2c6_over_b6h12_ratio', 'obs_cd47_cc2c6_over_b6h12_se', 0.0, 2.0, 'qpctl_pyroglutamate_pglu19_cc2c6_over_pan_cd47_b6h12_mfi_ratio'),
    ('locus1_complement_c3b_opsonized_macrophage_phagocytosis_pct', 'obs_complement_c3b_opsonized_phagocytosis_pct', 'obs_complement_c3b_opsonized_phagocytosis_se', 0.0, 100.0, 'c3b_igg_opsonized_primary_human_macrophage_phagocytosis_percentage'),
    ('graft_iec_vascular_endothelial_support_viability_pct', 'obs_iec_vascular_support_viability_pct', 'obs_iec_vascular_support_viability_se', 0.0, 100.0, 'post_stress_iec_pericyte_vascular_support_lineage_viability_percentage'),
    ('locus1_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_se', 0.0, 100.0, 'receptor_negative_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_specific_lysis_percentage'),
    ('gate0c_12pulse_unprimed_hypoxic_serial_lysis_pct', 'obs_12pulse_unprimed_hypoxic_serial_lysis_pct', 'obs_12pulse_unprimed_hypoxic_serial_lysis_se', 0.0, 100.0, 'unprimed_hour0_hypoxic_12pulse_serial_nk_ctl_specific_lysis_percentage'),
    ('stage7_recipient_insulin_clearance_tau_min', 'obs_recipient_insulin_clearance_tau_min', 'obs_recipient_insulin_clearance_tau_se', 1.0, 60.0, 'systemic_or_peritoneal_insulin_clearance_time_constant_minutes'),
    # Residual pluripotency post-purge, 20q11.21 BCL2L1 CNV, hazard-state proliferative survivor, Gate 0S-Retrieval en-bloc recovery, AP1903 insulin dump, peri-implant fluid film, 72h hypoxic CD47/HLA-E retention, single-cell p10 CD47, IL-12/IL-18 SIRPa-low NK, 21d primed B'-SCT ELISpot, and SA-FasL neutrophil chemotaxis endpoints (AUDIT-REF; DEC-047):
    ('gate0s_residual_pluripotent_post_purge_fraction', 'obs_residual_psc_effective_fraction', 'obs_residual_psc_effective_se', 0.0, 1.0, 'post_dual_purge_effective_residual_tra160_oct4_fraction'),
    ('gate0s_chr20q11_bcl2l1_cnv_ratio', 'obs_20q11_bcl2l1_cnv_ratio', 'obs_20q11_bcl2l1_cnv_se', 0.0, 10.0, 'chr20q11_21_bcl2l1_ddpcr_copy_number_ratio_diploid_2p0'),
    ('gate0s_hazard_state_prolif_survivor_count', 'obs_hazard_state_prolif_survivor_count_k', 'obs_hazard_state_prolif_survivor_se', 0.0, 1000000.0, 'hypoxic_isr_sorted_ki67_proliferative_ap1903_survivor_count_k'),
    ('gate0s_retrieval_en_bloc_recovery_fraction', 'obs_explant_en_bloc_recovery_fraction', 'obs_explant_en_bloc_recovery_se', 0.0, 1.0, 'gate0s_retrieval_30d_90d_en_bloc_pouch_recovery_fraction'),
    ('gate0s_ap1903_supernatant_insulin_dump_fraction_24h', 'obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h', 'obs_gate0s_ap1903_supernatant_insulin_dump_se', 0.0, 1.0, 'ap1903_ablation_24h_supernatant_stored_insulin_release_fraction'),
    ('gate0f_periimplant_fluid_film_thickness_um', 'obs_periimplant_fluid_film_um', 'obs_periimplant_fluid_film_se', 0.0, 10000.0, 'high_frequency_ultrasound_periimplant_seroma_fluid_film_um'),
    ('gate0a_hypoxic_72h_cd47_retention_molecules', 'obs_hypoxic_72h_cd47_retention_molecules', 'obs_hypoxic_72h_cd47_retention_se', 0.0, 1000000.0, 'hypoxic_1pct_o2_72h_cd47_surface_molecules_per_cell'),
    ('gate0a_hypoxic_72h_hlae_retention_molecules', 'obs_hypoxic_72h_hlae_retention_molecules', 'obs_hypoxic_72h_hlae_retention_se', 0.0, 1000000.0, 'hypoxic_1pct_o2_72h_hlae_surface_molecules_per_cell'),
    ('gate0a_cd47_single_cell_p10_molecules', 'obs_cd47_p10_single_cell_molecules', 'obs_cd47_p10_single_cell_se', 0.0, 1000000.0, 'quantibrite_single_cell_10th_percentile_cd47_molecules_per_cell'),
    ('gate0a_il12_il18_sirpa_low_nk_lysis_pct', 'obs_il12_il18_sirpa_low_nk_lysis_pct', 'obs_il12_il18_sirpa_low_nk_lysis_se', 0.0, 100.0, 'il12_il18_primed_sirpa_low_cd94_blocked_nk_lysis_percentage'),
    ('gate0b_b_prime_sct_21d_primed_cd8_elispot_sfu_per_1e6', 'obs_b_prime_sct_21d_primed_cd8_elispot_sfu_per_1e6', 'obs_b_prime_sct_21d_primed_cd8_elispot_se', 0.0, 10000.0, 'autologous_dc_21d_primed_cd8_ifng_elispot_sfu_per_1e6'),
    ('gate0g_safasl_neutrophil_chemotaxis_index', 'obs_safasl_neutrophil_chemotaxis_index', 'obs_safasl_neutrophil_chemotaxis_se', 0.0, 100.0, 'primary_human_neutrophil_boyce_chamber_chemotaxis_fold_index'),
    # Audit R5 pulsed-PK AP1903 washout kill, Step-1 sub-lethal AP1903 survivor count, 15 mmHg hypoxic f_ISR translation ratio, iEC p10 CD47 & ABO antigen, post-stress NKG2D/DNAM-1 lysis, and p16/SASP senescence endpoints (MW-R5-S1-01/03/06/08/10/22 & S2-F06; DEC-048):
    ('gate0s_ap1903_pulsed_pk_washout_log10_kill', 'obs_ap1903_pulsed_pk_log10_kill', 'obs_ap1903_pulsed_pk_log10_kill_se', 0.0, 10.0, 'pulsed_pk_100nM_2h_washout_2x24h_icasp9_log10_kill'),
    ('gate0s_ap1903_step1_sublethal_survivor_count', 'obs_ap1903_step1_sublethal_survivor_count_k', 'obs_ap1903_step1_sublethal_survivor_se', 0.0, 1000000.0, 'step1_10nM_ap1903_resistant_survivor_colony_count_k'),
    ('gate0a_hypoxic_15mmhg_f_isr_translation_ratio', 'obs_hypoxic_15mmhg_f_isr_translation_ratio', 'obs_hypoxic_15mmhg_f_isr_translation_se', 0.0, 2.0, 'puromycin_opp_or_s35_met_translation_rate_at_15mmhg_16p7mM_glc_ifng_tnfa_over_normoxic'),
    ('gate0a_iec_cd47_single_cell_p10_molecules', 'obs_iec_cd47_p10_molecules', 'obs_iec_cd47_p10_se', 0.0, 1000000.0, 've_cadherin_pos_iec_single_cell_p10_cd47_molecules_per_cell'),
    ('gate0a_iec_abo_ab_antigen_mfi_ratio', 'obs_iec_abo_ab_antigen_mfi_ratio', 'obs_iec_abo_ab_antigen_mfi_se', 0.0, 100.0, 've_cadherin_pos_iec_blood_group_ab_antigen_mfi_over_o_isotype'),
    ('gate0a_post_stress_nkg2d_dnam1_attributed_lysis_pct', 'obs_post_stress_nkg2d_dnam1_attributed_lysis_pct', 'obs_post_stress_nkg2d_dnam1_attributed_lysis_se', 0.0, 100.0, 'post_stress_anti_nkg2d_anti_dnam1_reversible_nk_lysis_percentage'),
    ('gate0s_senescent_p16_sasp_positive_cell_pct', 'obs_senescent_p16_sasp_positive_cell_pct', 'obs_senescent_p16_sasp_positive_cell_se', 0.0, 100.0, 'day28_post_stress_p16ink4a_sabgal_sasp_positive_cell_percentage'),
    # AUDIT-REF & Terra-6.1 R6 rare-event effective trial denominators, AP1903 monotherapy pulsed-PK kill, 1-hour acute AP1903 insulin burst, and single-arm Embodiment A* combined-stress secondary endpoints (AUDIT-REF & AUDIT-REF; DEC-050):
    ('gate0s_residual_pluripotent_effective_trials', 'obs_residual_psc_effective_trials', 'obs_residual_psc_effective_trials_se', 100.0, 1.0e11, 'effective_clonogenic_or_ddpcr_spike_in_adjusted_psc_trials_interrogated'),
    ('gate0s_hazard_state_prolif_effective_trials', 'obs_hazard_state_prolif_effective_trials', 'obs_hazard_state_prolif_effective_trials_se', 100.0, 1.0e11, 'effective_clonogenic_hazard_state_proliferative_trials_interrogated'),
    ('gate0s_ap1903_step1_effective_trials', 'obs_ap1903_step1_effective_trials', 'obs_ap1903_step1_effective_trials_se', 100.0, 1.0e11, 'effective_clonogenic_step1_10nM_ap1903_trials_interrogated'),
    ('gate0s_ap1903_monotherapy_pulsed_pk_log10_kill', 'obs_ap1903_monotherapy_pulsed_pk_log10_kill', 'obs_ap1903_monotherapy_pulsed_pk_log10_kill_se', 0.0, 10.0, 'ap1903_monotherapy_pulsed_pk_100nM_2h_washout_2x24h_hazard_state_log10_kill'),
    ('gate0s_ap1903_1h_acute_insulin_burst_iu', 'obs_gate0s_ap1903_1h_insulin_burst_iu', 'obs_gate0s_ap1903_1h_insulin_burst_iu_se', 0.0, 200.0, 'ap1903_ablation_peak_1h_acute_supernatant_insulin_burst_iu_per_442k_ieq'),
    ('combined_stress_coupon_embodiment_a_star_viable_beta_mass_pct', 'obs_combined_stress_embodiment_a_star_viable_beta_mass_pct', 'obs_combined_stress_embodiment_a_star_viable_beta_mass_se', 0.0, 100.0, 'combined_stress_embodiment_a_star_post_stress_viable_beta_mass_percentage'),
    ('combined_stress_coupon_embodiment_a_star_shutoff_t90_min', 'obs_combined_stress_embodiment_a_star_shutoff_t90_min', 'obs_combined_stress_embodiment_a_star_shutoff_t90_se', 0.1, 30.0, 'combined_stress_embodiment_a_star_post_stepdown_insulin_shutoff_t90_minutes'),
    ('combined_stress_coupon_embodiment_a_star_max_h2o2_uM', 'obs_combined_stress_embodiment_a_star_max_h2o2_uM', 'obs_combined_stress_embodiment_a_star_max_h2o2_se', 0.0, 1000.0, 'combined_stress_embodiment_a_star_peak_peroxide_uM'),
    ('combined_stress_coupon_embodiment_a_star_max_ph', 'obs_combined_stress_embodiment_a_star_max_ph', 'obs_combined_stress_embodiment_a_star_max_ph_se', 6.5, 8.5, 'combined_stress_embodiment_a_star_peak_alkaline_ph'),
    # Audit R6 & Audit R6 post-72h hypoxic + ISR single-cell p10 CD47/HLA-E retention and single-cell 20q11.21 BCL2L1 mosaicism endpoints (AUDIT-REF & AUDIT-REF; DEC-051/052):
    ('gate0a_hypoxic_72h_cd47_single_cell_p10_molecules', 'obs_hypoxic_72h_cd47_p10_single_cell_molecules', 'obs_hypoxic_72h_cd47_p10_single_cell_se', 0.0, 2000000.0, 'post_72h_hypoxic_isr_single_cell_10th_percentile_cd47_molecules_per_cell'),
    ('gate0a_hypoxic_72h_hlae_single_cell_p10_molecules', 'obs_hypoxic_72h_hlae_p10_single_cell_molecules', 'obs_hypoxic_72h_hlae_p10_single_cell_se', 0.0, 1000000.0, 'post_72h_hypoxic_isr_single_cell_10th_percentile_hlae_molecules_per_cell'),
    ('gate0s_chr20q11_bcl2l1_single_cell_mosaic_fraction', 'obs_20q11_bcl2l1_single_cell_mosaic_fraction', 'obs_20q11_bcl2l1_single_cell_mosaic_se', 0.0, 1.0, 'single_cell_dna_fish_or_scdna_chr20q11_21_bcl2l1_trisomic_mosaic_cell_fraction'),
    # Audit Full R6 Gate 0V viral sanctuary, post-stress 20q11.21 outgrowth, assembled-cassette Chandler-loop IBMIR, 5-transgene joint coverage, efferocytosis retention, functional Nutlin-3a p53, pericyte coverage, and assembled-cluster lineage/drift endpoints (MW-R6-S1-01..16; DEC-053):
    ('gate0v_cvb_day14_log10_viral_reduction', 'obs_gate0v_cvb_day14_log10_viral_reduction', 'obs_gate0v_cvb_day14_log10_viral_reduction_se', 0.0, 10.0, 'cvb3_cvb4_hcmv_day14_log10_viral_genome_copy_reduction'),
    ('gate0v_day28_persistent_vp1_ie1_positive_cell_pct', 'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct', 'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_se', 0.0, 100.0, 'day28_persistent_cvb_vp1_or_hcmv_ie1_positive_beta_cell_percentage'),
    ('gate0s_post_stress_20q11_bcl2l1_outgrowth_fold', 'obs_post_stress_20q11_bcl2l1_outgrowth_fold', 'obs_post_stress_20q11_bcl2l1_outgrowth_se', 0.0, 100.0, 'post_hypoxia_cytokine_sublethal_ap1903_20q11_21_bcl2l1_outgrowth_enrichment_fold'),
    ('gate0f_assembled_cassette_chandler_loop_sc5b9_fold', 'obs_assembled_cassette_chandler_loop_sc5b9_fold', 'obs_assembled_cassette_chandler_loop_sc5b9_se', 0.0, 100.0, 'assembled_cassette_recalcified_whole_blood_chandler_loop_sc5b9_fold_over_baseline'),
    ('gate0f_assembled_cassette_chandler_loop_platelet_retention_pct', 'obs_assembled_cassette_chandler_loop_platelet_retention_pct', 'obs_assembled_cassette_chandler_loop_platelet_retention_se', 0.0, 100.0, 'assembled_cassette_recalcified_whole_blood_chandler_loop_platelet_retention_percentage'),
    ('gate0a_hypoxic_72h_five_transgene_joint_above_threshold_fraction', 'obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction', 'obs_hypoxic_72h_five_transgene_joint_above_threshold_se', 0.0, 1.0, 'post_72h_hypoxic_isr_single_cell_joint_five_transgene_above_threshold_fraction'),
    ('gate0a_apoptotic_effero_clearance_retention_pct_of_wt', 'obs_apoptotic_effero_clearance_retention_pct_of_wt', 'obs_apoptotic_effero_clearance_retention_se', 0.0, 200.0, 'annexin_v_positive_apoptotic_corpse_mertk_axl_tim4_efferocytosis_clearance_percentage_of_wt'),
    ('gate0s_nutlin3a_p21_mdm2_fold_induction', 'obs_nutlin3a_p21_mdm2_fold_induction', 'obs_nutlin3a_p21_mdm2_fold_induction_se', 0.0, 100.0, 'functional_p53_nutlin3a_10uM_8h_cdkn1a_p21_and_mdm2_mrna_fold_induction'),
    ('gate0f_vascular_pericyte_ng2_pdgfrb_coverage_pct', 'obs_vascular_pericyte_ng2_pdgfrb_coverage_pct', 'obs_vascular_pericyte_ng2_pdgfrb_coverage_se', 0.0, 100.0, 'abluminal_ng2_pdgfrb_pericyte_coverage_of_cd31_iec_microvessels_percentage'),
    ('gate0m_assembled_cluster_slc18a1_enterochromaffin_pct', 'obs_assembled_cluster_slc18a1_enterochromaffin_pct', 'obs_assembled_cluster_slc18a1_enterochromaffin_se', 0.0, 100.0, 'assembled_cluster_slc18a1_fev_tph1_enterochromaffin_like_off_target_cell_percentage'),
    ('gate0m_assembled_cluster_ins_gcg_polyhormonal_pct', 'obs_assembled_cluster_ins_gcg_polyhormonal_pct', 'obs_assembled_cluster_ins_gcg_polyhormonal_se', 0.0, 100.0, 'assembled_cluster_ins_pos_gcg_pos_polyhormonal_cell_percentage'),
    ('gate0s_assembled_cluster_ki67_progenitor_pct', 'obs_assembled_cluster_ki67_progenitor_pct', 'obs_assembled_cluster_ki67_progenitor_se', 0.0, 100.0, 'assembled_cluster_ki67_sox9_residual_ductal_acinar_progenitor_cell_percentage'),
    ('gate0m_post_stress_mafa_nkx61_beta_identity_retention_pct', 'obs_post_stress_mafa_nkx61_beta_identity_retention_pct', 'obs_post_stress_mafa_nkx61_beta_identity_retention_se', 0.0, 100.0, 'post_stress_72h_ins_pos_nkx61_pos_mafa_pos_mature_beta_identity_retention_percentage'),
    # Audit R7 & Audit Full R6 supplementary / co-registered endpoint specs for complete DEC-053/055 parity (AUDIT-REF; DEC-053/055):
    ('gate0v_cvb4_10moi_72h_viral_titer_log10_pfu', 'obs_cvb4_10moi_72h_viral_titer_log10_pfu', 'obs_cvb4_10moi_72h_viral_titer_log10_pfu_se', 0.0, 10.0, 'cvb4_10moi_72h_residual_viral_titer_log10_pfu_per_ml'),
    ('gate0s_nutlin3a_p21_cdkn1a_fold_induction', 'obs_nutlin3a_p21_cdkn1a_fold_induction', 'obs_nutlin3a_p21_cdkn1a_fold_induction_se', 0.0, 100.0, 'functional_p53_nutlin3a_10uM_8h_cdkn1a_p21_mrna_fold_induction'),
    ('gate0m_sc_atac_ec_chromatin_spillover_pct', 'obs_sc_atac_ec_chromatin_spillover_pct', 'obs_sc_atac_ec_chromatin_spillover_se', 0.0, 100.0, 'single_cell_atac_enterochromaffin_chromatin_spillover_percentage'),
    ('gate0m_endocrine_synaptophysin_pos_purity_pct', 'obs_endocrine_synaptophysin_pos_purity_pct', 'obs_endocrine_synaptophysin_pos_purity_se', 0.0, 100.0, 'assembled_cluster_synaptophysin_positive_endocrine_purity_percentage'),
    ('gate0s_ductal_sox9_krt19_double_pos_pct', 'obs_ductal_sox9_krt19_double_pos_pct', 'obs_ductal_sox9_krt19_double_pos_se', 0.0, 100.0, 'assembled_cluster_sox9_krt19_double_positive_ductal_progenitor_percentage'),
    ('gate0m_enterochromaffin_fev_tph1_pos_pct', 'obs_enterochromaffin_fev_tph1_pos_pct', 'obs_enterochromaffin_fev_tph1_pos_se', 0.0, 100.0, 'assembled_cluster_fev_tph1_positive_enterochromaffin_like_cell_percentage'),
    ('gate0f_c3b_cd55_cd59_complement_deposition_mfi_ratio', 'obs_c3b_cd55_cd59_complement_deposition_mfi_ratio', 'obs_c3b_cd55_cd59_complement_deposition_mfi_se', 0.0, 100.0, 'assembled_cassette_c3b_c3d_complement_deposition_mfi_ratio_over_isotype'),
    ('gate0a_nkp46_innate_nk_lysis_pct', 'obs_nkp46_innate_nk_lysis_pct', 'obs_nkp46_innate_nk_lysis_se', 0.0, 100.0, 'nkp46_positive_innate_nk_cell_4h_specific_lysis_percentage'),
    ('gate0g_sa_fasl_microgel_nk_degranulation_cd107a_pct', 'obs_sa_fasl_microgel_nk_degranulation_cd107a_pct', 'obs_sa_fasl_microgel_nk_degranulation_cd107a_se', -50.0, 100.0, 'sa_fasl_microgel_nkp46_nk_delta_cd107a_degranulation_percentage_points'),
    ('gate0g_safasl_delta_tunel_pp', 'obs_safasl_delta_tunel_pp', 'obs_safasl_delta_tunel_se', -50.0, 100.0, 'sa_fasl_microgel_cytokine_stimulated_beta_cell_delta_tunel_apoptosis_percentage_points'),
    ('gate0g_safasl_zb4_anti_fas_rescue_pp', 'obs_safasl_zb4_anti_fas_rescue_pp', 'obs_safasl_zb4_anti_fas_rescue_se', -50.0, 100.0, 'sa_fasl_microgel_zb4_neutralizing_anti_fas_mab_apoptosis_rescue_percentage_points'),
    ('gate0g_safasl_treg_viable_retention_pct_of_thbd', 'obs_safasl_treg_viable_retention_pct_of_thbd', 'obs_safasl_treg_viable_retention_se', 0.0, 200.0, 'sa_fasl_microgel_foxp3_treg_viable_retention_percentage_of_sa_thbd_only_control'),
    ('gate0a_effero_clearance_pct', 'obs_effero_clearance_pct', 'obs_effero_clearance_se', 0.0, 200.0, 'annexin_v_apoptotic_corpse_efferocytosis_clearance_percentage_alias'),
    ('gate0f_iec_anastomosis_t50_h', 'obs_iec_anastomosis_t50_h', 'obs_iec_anastomosis_t50_se', 0.0, 168.0, 'iec_pericyte_microvascular_anastomosis_half_time_hours_alias'),
    ('gate0f_day28_functional_beta_survival_fraction', 'obs_day28_functional_beta_survival_fraction', 'obs_day28_functional_beta_survival_se', 0.0, 1.0, 'day28_chronic_pouch_functional_beta_cell_survival_fraction'),
    ('gate0m_pouch_iapp_oligomer_a11_dot_blot_fold', 'obs_pouch_iapp_oligomer_a11_dot_blot_fold', 'obs_pouch_iapp_oligomer_a11_dot_blot_se', 0.0, 100.0, 'confined_pouch_72h_high_glucose_a11_iapp_toxic_oligomer_dot_blot_fold_over_baseline'),
    ('gate0m_8mM_to_2mM_insulin_shutoff_ratio', 'obs_gate0m_8mM_to_2mM_insulin_shutoff_ratio', 'obs_gate0m_8mM_to_2mM_insulin_shutoff_se', 0.0, 1.0, 'perifusion_8mM_to_2mM_stepdown_15min_residual_insulin_secretion_ratio'),
    # Scored-cell denominator endpoints for the 6 proportion/percentage rare-event gates so zero-event (k=0) observations enforce exact Clopper-Pearson/Poisson 95% upper confidence bounds (AUDIT-REF & F2-03; DEC-071):
    ('gate0s_chr20q11_bcl2l1_single_cell_mosaic_cells_scored', 'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored', 'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored_se', 10.0, 1e9, 'single_cell_dna_fish_or_scdna_cells_scored_for_20q11_21_bcl2l1_mosaicism'),
    ('gate0v_day28_persistent_vp1_ie1_cells_scored', 'obs_gate0v_day28_persistent_vp1_ie1_cells_scored', 'obs_gate0v_day28_persistent_vp1_ie1_cells_scored_se', 10.0, 1e9, 'day28_beta_cells_scored_for_persistent_vp1_or_ie1_positivity'),
    ('gate0s_assembled_cluster_ki67_progenitor_cells_scored', 'obs_assembled_cluster_ki67_progenitor_cells_scored', 'obs_assembled_cluster_ki67_progenitor_cells_scored_se', 10.0, 1e9, 'assembled_cluster_cells_scored_for_ki67_sox9_residual_progenitors'),
    ('gate0s_ductal_sox9_krt19_double_pos_cells_scored', 'obs_ductal_sox9_krt19_double_pos_cells_scored', 'obs_ductal_sox9_krt19_double_pos_cells_scored_se', 10.0, 1e9, 'assembled_cluster_cells_scored_for_sox9_krt19_double_positive_ductal_progenitors'),
    ('gate0m_assembled_cluster_slc18a1_enterochromaffin_cells_scored', 'obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored', 'obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored_se', 10.0, 1e9, 'assembled_cluster_cells_scored_for_slc18a1_fev_tph1_enterochromaffin_cells'),
    ('gate0s_senescent_p16_sasp_cells_scored', 'obs_senescent_p16_sasp_cells_scored', 'obs_senescent_p16_sasp_cells_scored_se', 10.0, 1e9, 'day28_post_stress_cells_scored_for_p16ink4a_sabgal_sasp_positivity'),
    # Post-adjunct functional-depletion rechallenge endpoints for Arm 2 (conditioned then functionally depleted) sibling-aliquot qualification (IIP-0012; DEC-091):
    ('gate0d_arm2_post_depletion_viable_beta_recovery_pct', 'obs_arm2_post_depletion_viable_beta_recovery_pct', 'obs_arm2_post_depletion_viable_beta_recovery_se', 0.0, 150.0, 'post_adjunct_depletion_rechallenge_viable_ins_nkx61_beta_recovery_percentage_of_sham'),
    ('gate0d_arm2_post_depletion_perifusion_gsis_si', 'obs_arm2_post_depletion_perifusion_gsis_si', 'obs_arm2_post_depletion_perifusion_gsis_si_se', 0.1, 30.0, 'post_adjunct_depletion_rechallenge_dynamic_perifusion_gsis_stimulation_index'),
    ('gate0f_arm2_post_depletion_nhs_c5b9_mac_pct', 'obs_arm2_post_depletion_nhs_c5b9_mac_pct', 'obs_arm2_post_depletion_nhs_c5b9_mac_se', 0.0, 100.0, 'post_adjunct_depletion_20pct_nhs_c5b9_mac_positive_percentage'),
    # Two-strand knock-in splice-challenger IVT-RNA recovery, UMI isoform census, TOST host-exon retention, and dsRNA/ISG qualification endpoints (IIP-0015; DEC-092):
    ('gate0s_challenger_ivt_relative_recovery_eta', 'obs_challenger_ivt_relative_recovery_eta', 'obs_challenger_ivt_relative_recovery_se', 0.01, 2.0, 'pre_extraction_chimeric_ivt_rna_relative_recovery_ratio_eta'),
    ('gate0s_challenger_umi_molecules_profiled', 'obs_challenger_umi_molecules_profiled', 'obs_challenger_umi_molecules_profiled_se', 0.0, 1e9, 'unique_consensus_umi_cdna_molecules_profiled'),
    ('gate0s_challenger_aberrant_splice_umi_count_k', 'obs_challenger_aberrant_splice_umi_count_k', 'obs_challenger_aberrant_splice_umi_count_se', 0.0, 1e6, 'observed_aberrant_chimeric_or_cryptic_splice_umi_molecule_count_k'),
    ('gate0s_challenger_host_exon_splice_retention_ratio', 'obs_challenger_host_exon_splice_retention_ratio', 'obs_challenger_host_exon_splice_retention_se', 0.0, 5.0, 'host_exon_canonical_splice_retention_ratio_vs_unedited_parental'),
    ('gate0s_challenger_dsrna_j2_or_isg_fold_wt', 'obs_challenger_dsrna_j2_or_isg_fold_wt', 'obs_challenger_dsrna_j2_or_isg_fold_se', 0.0, 50.0, 'anti_dsrna_j2_or_ifit1_oas1_isg15_fold_over_unedited_parental'),
    # Stage-0 Gate 0A (Rule S10) covalent-SCT SPR, paired resting + 24h IFN-gamma NK subset, UL40 CTL, CD8-coreceptor/Treg, and non-self seam endpoints (IIP-0016; DEC-094):
    ('stage0_s10_covalent_sct_nkg2a_kd_uM', 'obs_hlae_sct_nkg2a_kd_uM', 'obs_hlae_sct_nkg2a_kd_se_uM', 0.001, 200.0, 'covalent_sct_spr_cd94_nkg2a_kd_uM'),
    ('stage0_s10_covalent_sct_nkg2c_kd_uM', 'obs_hlae_sct_nkg2c_kd_uM', 'obs_hlae_sct_nkg2c_kd_se_uM', 0.001, 500.0, 'covalent_sct_spr_cd94_nkg2c_kd_uM'),
    ('stage0_s10_wt_sct_reference_nkg2a_kd_uM', 'obs_wt_sct_nkg2a_kd_uM', 'obs_wt_sct_nkg2a_kd_se_uM', 0.001, 200.0, 'parallel_wt_hlae_sct_reference_cd94_nkg2a_kd_uM'),
    ('stage0_s10_nkg2a_pos_resting_nk_lysis_pct', 'obs_nkg2a_pos_resting_nk_lysis_pct', 'obs_nkg2a_pos_resting_nk_lysis_se', 0.0, 100.0, 'regime_ii_resting_nkg2a_pos_nkg2c_neg_kir_neg_nk_lysis_percentage'),
    ('stage0_s10_nkg2a_pos_ifng24h_nk_lysis_pct', 'obs_nkg2a_pos_ifng24h_nk_lysis_pct', 'obs_nkg2a_pos_ifng24h_nk_lysis_se', 0.0, 100.0, 'regime_i_ifng24h_nkg2a_pos_nkg2c_neg_kir_neg_nk_lysis_percentage'),
    ('stage0_s10_nkg2a_pos_resting_cd107a_pct', 'obs_nkg2a_pos_resting_cd107a_pct', 'obs_nkg2a_pos_resting_cd107a_se', 0.0, 100.0, 'regime_ii_resting_nkg2a_pos_nkg2c_neg_nk_cd107a_degranulation_percentage'),
    ('stage0_s10_nkg2a_pos_ifng24h_cd107a_pct', 'obs_nkg2a_pos_ifng24h_cd107a_pct', 'obs_nkg2a_pos_ifng24h_cd107a_se', 0.0, 100.0, 'regime_i_ifng24h_nkg2a_pos_nkg2c_neg_nk_cd107a_degranulation_percentage'),
    ('stage0_s10_nkg2c_kir_pos_resting_lysis_pct', 'obs_nkg2c_kir_pos_resting_lysis_pct', 'obs_nkg2c_kir_pos_resting_lysis_se', 0.0, 100.0, 'regime_ii_resting_nkg2c_pos_nkg2a_neg_kir2dl_pos_nk_lysis_percentage'),
    ('stage0_s10_nkg2c_kir_pos_ifng24h_lysis_pct', 'obs_nkg2c_kir_pos_ifng24h_lysis_pct', 'obs_nkg2c_kir_pos_ifng24h_lysis_se', 0.0, 100.0, 'regime_i_ifng24h_nkg2c_pos_nkg2a_neg_kir2dl_pos_nk_lysis_percentage'),
    ('stage0_s10_double_neg_nkg2c_resting_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_se', 0.0, 100.0, 'regime_ii_resting_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_percentage'),
    ('stage0_s10_double_neg_nkg2c_ifng24h_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_se', 0.0, 100.0, 'regime_i_ifng24h_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_percentage'),
    ('stage0_s10_hcmv_pooled_3subset_resting_lysis_pct', 'obs_hcmv_pooled_3subset_resting_nk_lysis_pct', 'obs_hcmv_pooled_3subset_resting_nk_lysis_se', 0.0, 100.0, 'regime_ii_resting_total_3subset_sorted_hcmv_nk_lysis_percentage'),
    ('stage0_s10_hcmv_pooled_3subset_ifng24h_lysis_pct', 'obs_hcmv_pooled_3subset_ifng24h_nk_lysis_pct', 'obs_hcmv_pooled_3subset_ifng24h_nk_lysis_se', 0.0, 100.0, 'regime_i_ifng24h_total_3subset_sorted_hcmv_nk_lysis_percentage'),
    ('stage0_s10_nkg2c_cd107a_fold_vs_b2m_ko', 'obs_nkg2c_cd107a_degranulation_fold_vs_b2m_ko', 'obs_nkg2c_cd107a_degranulation_fold_se', 0.0, 20.0, 'nkg2a_neg_nkg2c_pos_cd107a_degranulation_fold_over_b2m_ko_control'),
    ('stage0_s10_cmv_ul40_hlae_cd8_ctl_lysis_pct', 'obs_cmv_ul40_hlae_cd8_lysis_pct', 'obs_cmv_ul40_hlae_cd8_lysis_se', 0.0, 100.0, 'cmv_ul40_hlae_restricted_cd8_ctl_4h_specific_lysis_percentage'),
    ('stage0_s10_hcmv_pos_donor_count', 'obs_hcmv_pos_donor_count', 'obs_hcmv_pos_donor_count_se', 0.0, 1000.0, 'independent_hcmv_seropositive_nk_donors_assayed_integer_count'),
    ('stage0_s10_nk_failing_donor_count', 'obs_s10_nk_failing_donor_count', 'obs_s10_nk_failing_donor_count_se', 0.0, 1000.0, 'hcmv_nk_donors_exceeding_individual_subset_lysis_threshold_integer_count'),
    ('stage0_s10_ul40_reactive_donor_count', 'obs_ul40_reactive_donor_count', 'obs_ul40_reactive_donor_count_se', 0.0, 1000.0, 'tetramer_or_peptide_confirmed_ul40_reactive_hcmv_cd8_ctl_donors_integer_count'),
    ('stage0_s10_ul40_ctl_failing_donor_count', 'obs_ul40_ctl_failing_donor_count', 'obs_ul40_ctl_failing_donor_count_se', 0.0, 1000.0, 'ul40_reactive_cd8_ctl_donors_exceeding_5pct_or_2pp_over_wt_integer_count'),
    ('stage0_s10_hlae_cd8aa_spr_kd_uM', 'obs_hlae_cd8aa_kd_uM', 'obs_hlae_cd8aa_kd_se_uM', 0.001, 50000.0, 'covalent_sct_spr_cd8alpha_alpha_homodimer_kd_uM'),
    ('stage0_s10_hlae_cd8ab_spr_kd_uM', 'obs_hlae_cd8ab_kd_uM', 'obs_hlae_cd8ab_kd_se_uM', 0.001, 50000.0, 'covalent_sct_spr_cd8alpha_beta_heterodimer_kd_uM'),
    ('stage0_s10_cd8_treg_suppression_retention_pct_of_wt', 'obs_cd8_treg_suppression_retention_pct_of_wt', 'obs_cd8_treg_suppression_retention_se', 0.0, 200.0, 'hlae_restricted_cd8_treg_suppression_retention_percentage_of_wt_sct'),
    ('stage0_s10_cd8_treg_suppression_delta_vs_wt_pp', 'obs_cd8_treg_suppression_delta_vs_wt_pp', 'obs_cd8_treg_suppression_delta_vs_wt_se', -100.0, 100.0, 'hlae_restricted_cd8_treg_suppression_delta_mutant_minus_wt_percentage_points'),
    ('stage0_s10_ul40_ctl_lysis_delta_vs_wt_pp', 'obs_cmv_ul40_hlae_cd8_lysis_delta_vs_wt_pp', 'obs_cmv_ul40_hlae_cd8_lysis_delta_vs_wt_se', -100.0, 100.0, 'cmv_ul40_cd8_ctl_4h_lysis_delta_d227k_t228a_minus_wt_sct_percentage_points'),
    ('stage0_s10_nonself_dq2_dq8_indirect_apc_si', 'obs_nonself_dq2_dq8_indirect_apc_si', 'obs_nonself_dq2_dq8_indirect_apc_si_se', 0.1, 100.0, 'recipient_dq2_dq8_apc_indirect_cd4_stimulation_index_on_nonself_seam_or_alpha3'),
    ('stage0_s10_nonself_host_class1_cross_pres_lysis_pct', 'obs_nonself_host_class1_cross_pres_lysis_pct', 'obs_nonself_host_class1_cross_pres_lysis_se', 0.0, 100.0, 'host_hla_abc_cross_presented_cd8_ctl_specific_lysis_percentage'),
    # Stage-0 Gate 0B Design B2 (ISLET1_B2_5LOCUS_TYK2_PE) retained-HLA-C TYK2 p.Pro1104Ala (SAPTKFLEL) stratified OR-hold endpoints (IIP-0017; DEC-095):
    ('stage0_gate0b_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp', 'obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp', 'obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_se', -100.0, 100.0, '4digit_hlac_matched_cd8_lysis_delta_design_b2_minus_isogenic_bwt_percentage_points'),
    ('stage0_gate0b_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6', 'obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6', 'obs_b2_tyk2_p1104a_primed_cd8_elispot_se', 0.0, 5000.0, 'day21_saptkflel_primed_autologous_cd8_ifng_elispot_sfu_per_1e6_cells'),
    ('stage0_gate0b_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp', 'obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp', 'obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_se', -100.0, 100.0, 'ifng_treated_iec_pericyte_cd8_lysis_delta_design_b2_minus_isogenic_bwt_percentage_points'),
    # Stage-0 Locus 3 Mode-B Multi-miR(IL1R1+FAS) splicing, Drosha/Dicer 5' seed fidelity, AGO2 miR-375 retention, shield protein retention, receptor knockdown, and matched absolute GSIS output endpoints (IIP-0018; DEC-096):
    ('stage0_l3_multimir_mature_spliced_fraction', 'obs_l3_multimir_spliced_fraction', 'obs_l3_multimir_spliced_fraction_se', 0.0, 1.0, 'mature_exon1_exon2_spliced_transcript_fraction_0_to_1'),
    ('stage0_l3_multimir_5prime_seed_fidelity_pct', 'obs_l3_multimir_5prime_seed_fidelity_pct', 'obs_l3_multimir_5prime_seed_fidelity_se', 0.0, 100.0, 'small_rna_seq_exact_5prime_drosha_dicer_seed_register_percentage'),
    ('stage0_l3_multimir_ago2_mir375_retention_pct', 'obs_l3_multimir_ago2_mir375_retention_pct', 'obs_l3_multimir_ago2_mir375_retention_se', 0.0, 200.0, 'ago2_rip_endogenous_mir375_occupancy_retention_percentage_of_arm1'),
    ('stage0_l3_multimir_shield_protein_retention_pct', 'obs_l3_multimir_shield_protein_retention_pct', 'obs_l3_multimir_shield_protein_retention_se', 0.0, 200.0, 'min_of_cytosolic_serpinb9_and_surface_pdl1_protein_retention_percentage_of_arm1'),
    ('stage0_l3_multimir_min_receptor_knockdown_pct', 'obs_l3_multimir_min_receptor_knockdown_pct', 'obs_l3_multimir_min_receptor_knockdown_se', 0.0, 100.0, 'surface_flow_knockdown_percentage_vs_arm2_scramble_for_targeted_receptors'),
    ('stage0_l3_multimir_matched_gsis_output_retention_pct', 'obs_l3_multimir_matched_gsis_output_retention_pct', 'obs_l3_multimir_matched_gsis_output_retention_se', 0.0, 200.0, 'matched_absolute_perifusion_insulin_output_retention_percentage_of_arm1'),
    # Stage-0 explicit-3D transport hazard qualification gates for Rule S6 lateral spacing / dose-matched 80 um sizing, Readout 5 (Gate P05) wet-sieved microbeads + co-cast catalase barrier films, and Design Z two-stage empty-lumen + loaded closed-state confirmation (IIP-0019; DEC-097):
    ('stage0_s6_p10_nearest_neighbor_edge_gap_um', 'obs_s6_p10_nearest_neighbor_edge_gap_um', 'obs_s6_p10_nearest_neighbor_edge_gap_se', 0.0, 500.0, 'tenth_percentile_nearest_neighbor_cluster_edge_to_edge_gap_micrometers'),
    ('stage0_s6_cluster_diam_p90_um', 'obs_s6_cluster_diam_p90_um', 'obs_s6_cluster_diam_p90_se', 20.0, 300.0, 'ninetieth_percentile_assembled_sc_islet_cluster_diameter_micrometers'),
    ('stage0_s6_divider_14d_hydrated_po2_retention_pct', 'obs_s6_divider_14d_hydrated_po2_retention_pct', 'obs_s6_divider_14d_hydrated_po2_retention_se', 0.0, 150.0, 'day14_37c_hydrated_breathable_divider_oxygen_permeability_retention_percentage'),
    ('stage0_s6_hydrogel_shear_modulus_gprime_pa', 'obs_s6_hydrogel_shear_modulus_gprime_pa', 'obs_s6_hydrogel_shear_modulus_gprime_se', 10.0, 20000.0, 'pericellular_hydrogel_storage_shear_modulus_g_prime_pascals'),
    ('stage0_p05_microbead_volume_d99_um', 'obs_p05_microbead_volume_d99_um', 'obs_p05_microbead_volume_d99_se', 1.0, 600.0, 'laser_diffraction_volume_weighted_d99_pdms_cao2_microbead_diameter_micrometers'),
    ('stage0_p05_catalase_barrier_film_um', 'obs_p05_catalase_barrier_film_um', 'obs_p05_catalase_barrier_film_se', 0.0, 100.0, 'cell_free_immobilized_catalase_barrier_film_thickness_per_side_micrometers'),
    ('stage0_p05_islet_plane_h2o2_max_uM', 'obs_p05_islet_plane_h2o2_max_uM', 'obs_p05_islet_plane_h2o2_max_se', 0.0, 200.0, 'interstitial_islet_plane_peak_hydrogen_peroxide_micromolar'),
    ('stage0_p05_day7_catalase_retention_pct', 'obs_p05_day7_catalase_retention_pct', 'obs_p05_day7_catalase_retention_se', 0.0, 150.0, 'day7_dynamic_washout_active_film_catalase_retention_percentage'),
    ('stage0_p05_336h_retained_catalase_iu_333cm2', 'obs_p05_336h_retained_catalase_iu_333cm2', 'obs_p05_336h_retained_catalase_iu_333cm2_se', 0.0, 300000.0, 'retained_active_catalase_inventory_across_333p33cm2_two_pouch_bridge_iu'),
    ('stage0_design_z_empty_lumen_5site_min_po2_mmhg', 'obs_design_z_empty_lumen_5site_min_po2_mmhg', 'obs_design_z_empty_lumen_5site_min_po2_se', 0.0, 120.0, 'five_position_minimum_empty_prevascularized_lumen_phosphorescence_po2_mmhg'),
    ('stage0_design_z_vessel_standoff_p90_um', 'obs_design_z_vessel_standoff_p90_um', 'obs_design_z_vessel_standoff_p90_se', 0.0, 250.0, 'ninetieth_percentile_nearest_perfused_capillary_to_slab_standoff_micrometers'),
    ('stage0_design_z_loaded_surface_po2_mmhg', 'obs_design_z_loaded_surface_po2_mmhg', 'obs_design_z_loaded_surface_po2_se', 0.0, 120.0, 'loaded_closed_state_pericellular_cluster_surface_po2_mmhg'),
    ('stage0_design_z_loaded_core_po2_mmhg', 'obs_design_z_loaded_core_po2_mmhg', 'obs_design_z_loaded_core_po2_se', 0.0, 120.0, 'loaded_closed_state_minimum_cluster_core_po2_mmhg'),
    ('stage0_design_z_loaded_perifusion_gsis_si', 'obs_design_z_loaded_perifusion_gsis_si', 'obs_design_z_loaded_perifusion_gsis_se', 0.1, 50.0, 'loaded_closed_state_dynamic_perifusion_gsis_stimulation_index_preserving_readout2_ge3'),
    # Stage-0 Gate P05 recombinant human CAT (P04040) species/endotoxin/shed-mass lock and anti-PEG / Day-28 anti-SA humoral, DPR-maleimide hydrolysis, and primary human hepatocyte safety gates (IIP-0020; DEC-099):
    ('stage0_p05_catalase_endotoxin_eu_per_mg', 'obs_p05_catalase_endotoxin_eu_per_mg', 'obs_p05_catalase_endotoxin_eu_per_mg_se', 0.0, 50.0, 'recombinant_human_catalase_endotoxin_burden_eu_per_mg_protein'),
    ('stage0_p05_total_implant_catalase_endotoxin_eu', 'obs_p05_total_implant_catalase_endotoxin_eu', 'obs_p05_total_implant_catalase_endotoxin_eu_se', 0.0, 50.0, 'total_implant_catalase_endotoxin_burden_eu_across_333p33cm2'),
    ('stage0_p05_day7_shed_catalase_mass_pct', 'obs_p05_day7_shed_catalase_mass_pct', 'obs_p05_day7_shed_catalase_mass_se', 0.0, 100.0, 'day7_human_cat_elisa_and_mass_balance_shed_protein_percentage'),
    ('stage0_cuff_surface_c3b_deposition_pct_pos_ctrl', 'obs_cuff_surface_c3b_deposition_pct_pos_ctrl', 'obs_cuff_surface_c3b_deposition_se', 0.0, 200.0, 'surface_c1q_c3b_flow_deposition_percentage_of_positive_control_in_donor_sera'),
    ('stage0_cuff_soluble_sc5b9_pct_pos_ctrl', 'obs_cuff_soluble_sc5b9_pct_pos_ctrl', 'obs_cuff_soluble_sc5b9_se', 0.0, 200.0, 'fluid_phase_soluble_sc5b9_percentage_of_positive_control_in_donor_sera'),
    ('stage0_cuff_macrophage_fcgammar_cytokine_fold', 'obs_cuff_macrophage_fcgammar_cytokine_fold', 'obs_cuff_macrophage_fcgammar_cytokine_se', 0.1, 50.0, 'primary_human_macrophage_fcgammar_tnfa_il1b_fold_over_blank_gel'),
    ('stage0_cuff_primed_serum_fasl_ec50_fold_shift', 'obs_cuff_primed_serum_fasl_ec50_fold_shift', 'obs_cuff_primed_serum_fasl_ec50_fold_se', 0.05, 50.0, 'day28_primed_serum_activated_cd8_fasl_ec50_fold_shift_two_sided_0p5_to_2p0'),
    ('stage0_cuff_primed_serum_fasl_emax_apoptosis_pct', 'obs_cuff_primed_serum_fasl_emax_apoptosis_pct', 'obs_cuff_primed_serum_fasl_emax_apoptosis_se', 0.0, 100.0, 'day28_primed_serum_activated_cd8_maximal_fasl_apoptosis_percentage'),
    ('stage0_cuff_primed_serum_thbd_apc_retention_pct', 'obs_cuff_primed_serum_thbd_apc_retention_pct', 'obs_cuff_primed_serum_thbd_apc_retention_se', 0.0, 150.0, 'day28_primed_serum_thrombomodulin_protein_c_activation_retention_percentage'),
    ('stage0_cuff_primed_serum_cfh_c3b_cleavage_retention_pct', 'obs_cuff_primed_serum_cfh_c3b_cleavage_retention_pct', 'obs_cuff_primed_serum_cfh_c3b_cleavage_retention_se', 0.0, 150.0, 'day28_primed_serum_factor_h_factor_i_c3b_ic3b_cleavage_retention_percentage'),
    ('stage0_cuff_dpr_ring_hydrolysis_pct', 'obs_cuff_dpr_ring_hydrolysis_pct', 'obs_cuff_dpr_ring_hydrolysis_se', 0.0, 100.0, 'lcms_thiosuccinimide_ring_hydrolysis_completion_percentage_pre_cell_contact'),
    ('stage0_cuff_14d_hsa_shed_bioactive_fasl_pct', 'obs_cuff_14d_hsa_shed_bioactive_fasl_pct', 'obs_cuff_14d_hsa_shed_bioactive_fasl_se', 0.0, 100.0, 'day14_37c_50mg_ml_hsa_total_free_and_polymer_fragment_bioactive_fasl_release_percentage'),
    ('stage0_cuff_h1a_sec_mals_hexamer_purity_pct', 'obs_cuff_h1a_sec_mals_hexamer_purity_pct', 'obs_cuff_h1a_sec_mals_hexamer_purity_se', 0.0, 100.0, 'sec_mals_hexameric_h1a_adipoq_fasl_oligomer_fraction_percentage'),
    ('stage0_cuff_unreacted_free_thiol_pct', 'obs_cuff_unreacted_free_thiol_pct', 'obs_cuff_unreacted_free_thiol_se', 0.0, 100.0, 'post_conjugation_unreacted_accessible_free_thiol_percentage'),
    ('stage0_cuff_hepatocyte_72h_caspase37_fold', 'obs_cuff_hepatocyte_72h_caspase37_fold', 'obs_cuff_hepatocyte_72h_caspase37_se', 0.1, 50.0, 'primary_human_hepatocyte_72h_caspase3_7_fold_over_vehicle_leachate_and_mmp_digest'),
    ('stage0_cuff_hepatocyte_72h_alt_ldh_fold', 'obs_cuff_hepatocyte_72h_alt_ldh_fold', 'obs_cuff_hepatocyte_72h_alt_ldh_se', 0.1, 50.0, 'primary_human_hepatocyte_72h_alt_ldh_release_fold_over_vehicle_leachate_and_mmp_digest'),
    ('stage0_cuff_serum_donor_count', 'obs_cuff_serum_donor_count', 'obs_cuff_serum_donor_count_se', 0.0, 500.0, 'unpooled_human_serum_biological_donor_count_across_humoral_panel'),
    ('stage0_cuff_serum_failing_donor_count', 'obs_cuff_serum_failing_donor_count', 'obs_cuff_serum_failing_donor_count_se', 0.0, 500.0, 'failing_individual_donor_count_in_humoral_or_hepatocyte_panel'),
    # Stage-0 120-h lineage-resolved metabolic-rest OCR and 4-arm Nec-1s x Ferrostatin-1 conditioning, pellet-LOQ washout & NQO1 gates (IIP-0021; DEC-100):
    ('stage0_rest_120h_graft_ocr_ratio', 'obs_rest_120h_graft_ocr_ratio', 'obs_rest_120h_graft_ocr_ratio_se', 0.05, 3.0, 'assembled_60_30_10_plus_iec_0_to_120h_simultaneous_graft_ocr_ratio_vs_vehicle'),
    ('stage0_rest_depot_release_half_life_h', 'obs_rest_depot_release_half_life_h', 'obs_rest_depot_release_half_life_se_h', 0.001, 500.0, 'encapsulated_metabolic_rest_microdepot_release_half_life_hours_in_assembled_gel'),
    ('stage0_rest_post_washout_gsis_si', 'obs_rest_post_washout_gsis_si', 'obs_rest_post_washout_gsis_se', 0.1, 50.0, 'post_washout_dynamic_perifusion_gsis_stimulation_index'),
    ('stage0_rest_post_washout_gsis_output_retention_pct', 'obs_rest_post_washout_gsis_output_retention_pct', 'obs_rest_post_washout_gsis_output_retention_se', 0.0, 200.0, 'post_washout_absolute_gsis_insulin_output_retention_pct_of_unconditioned_control'),
    ('stage0_rest_post_washout_alpha_glucagon_fold', 'obs_rest_post_washout_alpha_glucagon_fold', 'obs_rest_post_washout_alpha_glucagon_se', 0.1, 30.0, 'post_washout_low_glucose_alpha_cell_glucagon_fold_induction'),
    ('stage0_cond_sc_beta_viable_yield_pct', 'obs_cond_sc_beta_viable_yield_pct', 'obs_cond_sc_beta_viable_yield_se', 0.0, 150.0, 'post_72h_35mmhg_h2o2_cytokine_cpeptide_sc_beta_absolute_viable_yield_pct'),
    ('stage0_cond_iec_viable_yield_pct', 'obs_cond_iec_viable_yield_pct', 'obs_cond_iec_viable_yield_se', 0.0, 150.0, 'post_72h_35mmhg_h2o2_cytokine_cdh5_plvap_iec_absolute_viable_yield_pct'),
    ('stage0_cond_final_wash_supernatant_nM', 'obs_cond_final_wash_supernatant_nM', 'obs_cond_final_wash_supernatant_se_nM', 0.0, 10000.0, 'pre_implant_final_wash_supernatant_lcms_residual_conditioning_agent_nM'),
    ('stage0_cond_pellet_residual_over_loq_ratio', 'obs_cond_pellet_residual_over_loq_ratio', 'obs_cond_pellet_residual_over_loq_se', 0.0, 1000.0, 'washed_cell_pellet_extract_lcms_residual_concentration_divided_by_validated_loq'),
    ('stage0_cond_pellet_extraction_recovery_pct', 'obs_cond_pellet_extraction_recovery_pct', 'obs_cond_pellet_extraction_recovery_se', 0.0, 150.0, 'matrix_spiked_washed_cell_pellet_lcms_extraction_recovery_percentage'),
    ('stage0_cond_endocrine_nqo1_activity_pct', 'obs_cond_endocrine_nqo1_activity_pct', 'obs_cond_endocrine_nqo1_activity_se', 0.0, 200.0, 'undiluted_endocrine_compartment_cellular_nqo1_enzymatic_activity_pct_of_vehicle'),
    ('stage0_cond_vascular_nqo1_activity_pct', 'obs_cond_vascular_nqo1_activity_pct', 'obs_cond_vascular_nqo1_activity_se', 0.0, 200.0, 'undiluted_vascular_iec_compartment_cellular_nqo1_enzymatic_activity_pct_of_vehicle'),
    ('stage0_cond_post_washout_gsis_si', 'obs_cond_post_washout_gsis_si', 'obs_cond_post_washout_gsis_se', 0.1, 50.0, 'post_conditioning_and_washout_dynamic_perifusion_gsis_stimulation_index'),
    # Stage-0 Track-2 B'-SCT vs HLA-G1 NK rescue split, B2M-null PVR/CD226 NK sign-inversion check, and lineage-separated SC-beta vs iEC necroptosis attribution (IIP-0022; DEC-101):
    ('stage0_track2_bsct_kir_nkg2a_neg_nk_lysis_pct', 'obs_track2_bsct_kir_nkg2a_neg_nk_lysis_pct', 'obs_track2_bsct_kir_nkg2a_neg_nk_lysis_se', 0.0, 100.0, 'sorted_c1_c2_kir2dl_pos_nkg2a_neg_primary_nk_lysis_pct_on_track2_bsct'),
    ('stage0_track2_bsct_donor_count', 'obs_track2_bsct_donor_count', 'obs_track2_bsct_donor_count_se', 0.0, 500.0, 'independent_kir_typed_primary_nk_donor_count_for_track2_bsct'),
    ('stage0_track2_bsct_failing_donor_count', 'obs_track2_bsct_failing_donor_count', 'obs_track2_bsct_failing_donor_count_se', 0.0, 500.0, 'individual_donors_exceeding_8pct_lysis_on_track2_bsct'),
    ('stage0_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct', 'obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct', 'obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_se', 0.0, 100.0, 'sorted_lilrb1_neg_nkg2a_neg_kir_pos_primary_nk_lysis_pct_on_hlag1_challenger'),
    ('stage0_track2_hlag1_donor_count', 'obs_track2_hlag1_donor_count', 'obs_track2_hlag1_donor_count_se', 0.0, 500.0, 'independent_primary_nk_donor_count_for_hlag1_lilrb1_neg_subset'),
    ('stage0_track2_hlag1_failing_donor_count', 'obs_track2_hlag1_failing_donor_count', 'obs_track2_hlag1_failing_donor_count_se', 0.0, 500.0, 'individual_donors_exceeding_8pct_lysis_on_hlag1_lilrb1_neg_subset'),
    ('stage0_b2m_null_pvr_cd226_nk_lysis_pct', 'obs_b2m_null_pvr_cd226_nk_lysis_pct', 'obs_b2m_null_pvr_cd226_nk_lysis_se', 0.0, 100.0, 'sorted_tigit_low_cd226_high_nkg2a_neg_kir_pos_nk_lysis_pct_on_b2m_null_sc_islets'),
    # Stage-0 Design Y wildcard (Y-01..Y-04) fail-closed endpoints + rare-event/age/lineage denominators (IIP-0023; DEC-102):
    ('stage0_design_y01_max_lineage_escape_pct', 'obs_design_y01_max_lineage_escape_pct', 'obs_design_y01_max_lineage_escape_se', 0.0, 100.0, 'max_single_cell_classical_hla_abc_hla_ii_or_cxcl10_escape_pct_across_4_lineages_ge_60d'),
    ('stage0_design_y01_lineage_cells_scored', 'obs_design_y01_lineage_cells_scored', 'obs_design_y01_lineage_cells_scored_se', 0.0, 10000000.0, 'min_single_cells_scored_per_lineage_across_beta_alpha_delta_iec'),
    ('stage0_design_y01_cbe_offtarget_pct', 'obs_design_y01_cbe_offtarget_pct', 'obs_design_y01_cbe_offtarget_se', 0.0, 100.0, 'tadcbed_a_to_g_rescue_or_rhampseq_offtarget_editing_percentage'),
    ('stage0_design_y02_pre_enrichment_ki_pct', 'obs_design_y02_pre_enrichment_ki_pct', 'obs_design_y02_pre_enrichment_ki_se', 0.0, 100.0, 'pre_enrichment_full_length_attb_minicircle_knockin_percentage_in_hpscs'),
    ('stage0_design_y02_clonal_ddpcr_cn', 'obs_design_y02_clonal_ddpcr_cn', 'obs_design_y02_clonal_ddpcr_cn_se', 0.0, 10.0, 'selected_clone_biallelic_payload_ddpcr_copy_number'),
    ('stage0_design_y02_wt_b2m_residual_allele_frac', 'obs_design_y02_wt_b2m_residual_allele_frac', 'obs_design_y02_wt_b2m_residual_allele_se', 0.0, 1.0, 'out_out_junction_pcr_and_phased_ont_residual_wt_b2m_allele_fraction'),
    ('stage0_design_y02_backbone_kanr_vaf', 'obs_design_y02_backbone_kanr_vaf', 'obs_design_y02_backbone_kanr_vaf_se', 0.0, 1.0, 'bacterial_kanr_or_ori_plasmid_backbone_variant_allele_fraction'),
    ('stage0_design_y02_wgs_pseudosite_k', 'obs_design_y02_wgs_pseudosite_k', 'obs_design_y02_wgs_pseudosite_k_se', 0.0, 1000.0, 'observed_non_target_bxb1_pseudosite_insertions_across_ot1_ot2_ot4_and_longread_wgs'),
    ('stage0_design_y03_iec_beta_selectivity_ratio', 'obs_design_y03_iec_beta_selectivity_ratio', 'obs_design_y03_iec_beta_selectivity_se', 0.0, 10000.0, 'quantibrite_iec_to_sc_beta_cd200_or_cd39_expression_ratio'),
    ('stage0_design_y03_beta_mir375_retention_pct', 'obs_design_y03_beta_mir375_retention_pct', 'obs_design_y03_beta_mir375_retention_se', 0.0, 200.0, 'ago2_rip_endogenous_beta_mir375_occupancy_retention_percentage'),
    ('stage0_design_y03_coculture_gsis_si', 'obs_design_y03_coculture_gsis_si', 'obs_design_y03_coculture_gsis_se', 0.1, 50.0, 'dynamic_perifusion_gsis_stimulation_index_under_100uM_amp_atp_iec_coculture'),
    ('stage0_design_y03_gsis_retention_pct', 'obs_design_y03_gsis_retention_pct', 'obs_design_y03_gsis_retention_se', 0.0, 200.0, 'perifusion_gsis_si_retention_percentage_vs_unarmed_iec_control'),
    ('stage0_design_y03_supernatant_adenosine_nM', 'obs_design_y03_supernatant_adenosine_nM', 'obs_design_y03_supernatant_adenosine_se_nM', 0.0, 10000.0, 'perifusion_effluent_adenosine_concentration_nM_at_16p7mM_glucose'),
    ('stage0_design_y03_dpcpx_gsis_delta_pct', 'obs_design_y03_dpcpx_gsis_delta_pct', 'obs_design_y03_dpcpx_gsis_delta_se', 0.0, 200.0, 'dpcpx_adora1_antagonist_gsis_unmasking_percentage_shift'),
    ('stage0_design_y04_multicistron_retention_pct', 'obs_design_y04_multicistron_retention_pct', 'obs_design_y04_multicistron_retention_se', 0.0, 100.0, 'day60_single_cell_5cloak_beta_and_6cloak_iec_copositivity_percentage'),
    ('stage0_design_y04_nhs20_c5b9_mac_pct', 'obs_design_y04_nhs20_c5b9_mac_pct', 'obs_design_y04_nhs20_c5b9_mac_se', 0.0, 100.0, 'twenty_pct_normal_human_serum_c5b9_mac_deposition_percentage'),
    ('stage0_design_y04_indirect_apc_cd4_si', 'obs_design_y04_indirect_apc_cd4_si', 'obs_design_y04_indirect_apc_cd4_se', 0.1, 50.0, 'hla_dq2_dq8_host_apc_indirect_cd4_stimulation_index'),
    ('stage0_design_y04_ap1903_clonogenic_survival_frac', 'obs_design_y04_ap1903_clonogenic_survival_frac', 'obs_design_y04_ap1903_clonogenic_survival_se', 0.0, 1.0, 'tier6b_ap1903_hypoxic_cytokine_clonogenic_survival_fraction'),
    ('stage0_design_y_sampling_age_days', 'obs_design_y_sampling_age_days', 'obs_design_y_sampling_age_days_se', 0.0, 365.0, 'total_post_edit_differentiation_and_cytokine_rechallenge_age_days'),
    ('stage0_design_y02_ont_intact_homolog_pct', 'obs_design_y02_ont_intact_homolog_pct', 'obs_design_y02_ont_intact_homolog_se', 0.0, 100.0, 'phased_15_25kb_ont_full_length_intact_homolog_percentage'),
    ('stage0_design_y04_all4_lineages_scored_count', 'obs_design_y04_all4_lineages_scored_count', 'obs_design_y04_all4_lineages_scored_se', 0.0, 10.0, 'number_of_product_lineages_beta_alpha_delta_iec_passing_gy009_and_y01'),
    ('stage0_design_y04_ap1903_survivor_count_k', 'obs_design_y04_ap1903_survivor_count_k', 'obs_design_y04_ap1903_survivor_count_k_se', 0.0, 1000000.0, 'observed_ap1903_clonogenic_survivor_colonies_k'),
    ('stage0_design_y04_ap1903_effective_trials', 'obs_design_y04_ap1903_effective_trials', 'obs_design_y04_ap1903_effective_trials_se', 0.0, 1e12, 'effective_clonogenic_spike_adjusted_trials_n_times_f_clon_times_eta'),
    # Stage-0 atomic provenance, donor-bank kinship, bank/lot 1q32.1 MDM4, Readout 2 perifusion kinetics, CD16+ ADCC, catalase mass, Gate 0S dual shedding, Student-t replicates & 3-metric IAPP gates (IIP-0024; DEC-103):
    ('stage0_mcb_donor_pairwise_ibd_pi_hat', 'obs_mcb_donor_pairwise_ibd_pi_hat', 'obs_mcb_donor_pairwise_ibd_pi_hat_se', 0.0, 1.0, 'genome_wide_pairwise_ibd_pi_hat_between_independent_abo_oo_donor_seedlines'),
    ('stage0_mcb_1q32_mdm4_cnv_ratio', 'obs_mcb_1q32_mdm4_cnv_ratio', 'obs_mcb_1q32_mdm4_cnv_se', 0.0, 10.0, 'mcb_bulk_ddpcr_diploid_copy_number_at_1q32p1_mdm4_normal_2p0'),
    ('stage0_mcb_1q32_mdm4_single_cell_mosaic_fraction', 'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction', 'obs_mcb_1q32_mdm4_single_cell_mosaic_se', 0.0, 1.0, 'mcb_single_cell_dna_fish_1q32p1_mdm4_gain_mosaic_proportion'),
    ('stage0_mcb_1q32_mdm4_single_cell_mosaic_cells_scored', 'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored', 'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored_se', 1.0, 100000000.0, 'mcb_single_cells_scored_for_1q32p1_mdm4_mosaicism'),
    ('stage0_lot_1q32_mdm4_cnv_ratio', 'obs_1q32_mdm4_cnv_ratio', 'obs_1q32_mdm4_cnv_se', 0.0, 10.0, 'stage7_lot_bulk_ddpcr_diploid_copy_number_at_1q32p1_mdm4_normal_2p0'),
    ('stage0_lot_1q32_mdm4_single_cell_mosaic_fraction', 'obs_1q32_mdm4_single_cell_mosaic_fraction', 'obs_1q32_mdm4_single_cell_mosaic_se', 0.0, 1.0, 'stage7_lot_single_cell_dna_fish_1q32p1_mdm4_gain_mosaic_proportion'),
    ('stage0_lot_1q32_mdm4_single_cell_mosaic_cells_scored', 'obs_1q32_mdm4_single_cell_mosaic_cells_scored', 'obs_1q32_mdm4_single_cell_mosaic_cells_scored_se', 1.0, 100000000.0, 'stage7_lot_single_cells_scored_for_1q32p1_mdm4_mosaicism'),
    ('stage0_readout2_perifusion_gsis_si', 'obs_readout2_perifusion_gsis_si', 'obs_readout2_perifusion_gsis_si_se', 0.1, 100.0, 'readout2_dynamic_perifusion_gsis_stimulation_index'),
    ('stage0_readout2_first_phase_peak_to_basal_ratio', 'obs_readout2_first_phase_peak_to_basal_ratio', 'obs_readout2_first_phase_peak_to_basal_se', 0.1, 100.0, 'readout2_first_phase_peak_to_basal_secretion_ratio'),
    ('stage0_readout2_first_phase_time_to_peak_min', 'obs_readout2_first_phase_time_to_peak_min', 'obs_readout2_first_phase_time_to_peak_se', 0.0, 60.0, 'readout2_minutes_from_high_glucose_step_to_first_phase_peak'),
    ('stage0_readout2_stepdown_shutoff_t50_min', 'obs_readout2_stepdown_shutoff_t50_min', 'obs_readout2_stepdown_shutoff_t50_se', 0.0, 60.0, 'readout2_stepdown_secretion_arrest_half_time_t50_minutes_distinct_from_gate0m_t90'),
    ('stage0_exact_product_net_cd16_adcc_lysis_pct', 'obs_exact_product_net_cd16_adcc_lysis_pct', 'obs_exact_product_net_cd16_adcc_lysis_se', 0.0, 100.0, 'paired_igg_depleted_cd16_pos_nk_net_adcc_specific_lysis_percentage'),
    ('stage0_cd16_adcc_positive_control_lysis_pct', 'obs_cd16_adcc_positive_control_lysis_pct', 'obs_cd16_adcc_positive_control_lysis_se', 0.0, 100.0, 'opsonized_positive_control_cd16_pos_nk_adcc_specific_lysis_percentage'),
    ('stage0_p05_total_loaded_catalase_protein_mg', 'obs_p05_total_loaded_catalase_protein_mg', 'obs_p05_total_loaded_catalase_protein_mg_se', 0.0, 1000.0, 'total_loaded_recombinant_human_p04040_catalase_protein_mass_mg_across_333p33cm2'),
    ('stage0_gate0s_shedding_ddpcr_events_k', 'obs_gate0s_shedding_ddpcr_events_k', 'obs_gate0s_shedding_ddpcr_events_k_se', 0.0, 1000000.0, 'day90_sister_unit_ex_vivo_ddpcr_extra_device_shedding_positive_events_k'),
    ('stage0_gate0s_shedding_ddpcr_trials_n', 'obs_gate0s_shedding_ddpcr_trials_n', 'obs_gate0s_shedding_ddpcr_trials_n_se', 1.0, 10000000000.0, 'day90_sister_unit_ex_vivo_ddpcr_qualified_genome_equivalents_n'),
    ('stage0_gate0s_shedding_viable_events_k', 'obs_gate0s_shedding_viable_events_k', 'obs_gate0s_shedding_viable_events_k_se', 0.0, 1000000.0, 'day90_sister_unit_viable_outgrowth_extra_device_shedding_positive_events_k'),
    ('stage0_gate0s_shedding_viable_trials_n', 'obs_gate0s_shedding_viable_trials_n', 'obs_gate0s_shedding_viable_trials_n_se', 1.0, 10000000000.0, 'day90_sister_unit_viable_outgrowth_qualified_plated_equivalents_n'),
    ('stage0_independent_lot_replicates_n', 'obs_stage0_independent_lot_replicates_n', 'obs_stage0_independent_lot_replicates_n_se', 1.0, STAGE0_MAX_INDEPENDENT_LOT_REPLICATES, 'independent_differentiation_lot_replicates_n_for_exact_student_t_bounds'),
    ('stage0_iapp_3metric_thioflavin_s_area_pct', 'obs_iapp_3metric_thioflavin_s_area_pct', 'obs_iapp_3metric_thioflavin_s_area_se', 0.0, 100.0, 'colocalized_thioflavin_s_positive_amyloid_fibril_area_percentage'),
    ('stage0_iapp_3metric_soluble_oligomer_fold', 'obs_iapp_3metric_soluble_oligomer_fold', 'obs_iapp_3metric_soluble_oligomer_se', 0.0, 50.0, 'soluble_iapp_oligomer_elisa_fold_over_quantifiable_age_matched_reference'),
    ('stage0_iapp_3metric_insoluble_fibril_fold', 'obs_iapp_3metric_insoluble_fibril_fold', 'obs_iapp_3metric_insoluble_fibril_se', 0.0, 50.0, 'formic_acid_extracted_insoluble_iapp_fold_over_quantifiable_age_matched_reference'),
    # Stage-0 Hurdle-4 two-timepoint human P04040 catalase envelope, exact z=70.0 um gel-interface H2O2, separated Day-30 vs Day-90 Gate 0S en-bloc retrieval, and zero-primate Contrasts B1, B2, B3 endpoints (IIP-0025; DEC-104):
    ('stage0_p05_precure_active_catalase_input_u_333cm2', 'obs_p05_precure_active_catalase_input_u_333cm2', 'obs_p05_precure_active_catalase_input_se_u', 0.0, 1000000.0, 'precure_input_recombinant_human_p04040_catalase_units_across_333p33cm2'),
    ('stage0_p05_postcure_t0_active_catalase_iu_333cm2', 'obs_p05_postcure_t0_active_catalase_iu_333cm2', 'obs_p05_postcure_t0_active_catalase_se_iu', 0.0, 1000000.0, 'postcure_t0_immobilized_active_catalase_iu_across_333p33cm2'),
    ('stage0_p05_cure_immobilization_yield_pct', 'obs_p05_cure_immobilization_yield_pct', 'obs_p05_cure_immobilization_yield_se', 0.0, 100.0, 'postcure_t0_over_precure_input_catalase_immobilization_yield_percentage'),
    ('stage0_p05_336h_catalase_retention_pct', 'obs_p05_336h_catalase_retention_pct', 'obs_p05_336h_catalase_retention_se', 0.0, 150.0, 'hour336_day14_dynamic_washout_active_film_catalase_retention_percentage_of_postcure_t0'),
    ('stage0_p05_z70_gel_interface_h2o2_max_uM', 'obs_p05_z70_gel_interface_h2o2_max_uM', 'obs_p05_z70_gel_interface_h2o2_max_se', 0.0, 200.0, 'exact_z70um_catalase_film_to_gel_interface_peak_24h_repeat_burst_h2o2_micromolar'),
    ('stage0_gate0s_day30_en_bloc_recovery_fraction', 'obs_gate0s_day30_en_bloc_recovery_fraction', 'obs_gate0s_day30_en_bloc_recovery_se', 0.0, 1.0, 'gate0s_day30_en_bloc_vascularized_pouch_recovery_fraction'),
    ('stage0_gate0s_day90_en_bloc_recovery_fraction', 'obs_gate0s_day90_en_bloc_recovery_fraction', 'obs_gate0s_day90_en_bloc_recovery_se', 0.0, 1.0, 'gate0s_day90_en_bloc_vascularized_pouch_recovery_fraction'),
    ('stage0_contrast_b1_paired_min_face_po2_gain_mmhg', 'obs_contrast_b1_paired_min_face_po2_gain_mmhg', 'obs_contrast_b1_paired_min_face_po2_gain_se', -20.0, 40.0, 'contrast_b1_paired_80um_minus_100um_dose_matched_worst_face_po2_gain_mmhg'),
    ('stage0_contrast_b1_80um_worst_face_po2_mmhg', 'obs_contrast_b1_80um_worst_face_po2_mmhg', 'obs_contrast_b1_80um_worst_face_po2_se', 0.0, 120.0, 'contrast_b1_dose_matched_80um_worst_face_pericellular_po2_mmhg'),
    ('stage0_contrast_b1_80um_perifusion_gsis_si', 'obs_contrast_b1_80um_perifusion_gsis_si', 'obs_contrast_b1_80um_perifusion_gsis_se', 0.1, 50.0, 'contrast_b1_dose_matched_80um_post_hypoxia_perifusion_gsis_stimulation_index'),
    ('stage0_contrast_b2_human_exvivo_viable_window_h', 'obs_contrast_b2_human_exvivo_viable_window_h', 'obs_contrast_b2_human_exvivo_viable_window_se_h', 0.0, 336.0, 'contrast_b2_paired_declined_human_omental_suprafascial_viable_perfusion_window_hours'),
    ('stage0_contrast_b2_omental_vs_suprafascial_po2_gain_mmhg', 'obs_contrast_b2_omental_vs_suprafascial_po2_gain_mmhg', 'obs_contrast_b2_omental_vs_suprafascial_po2_gain_se', -30.0, 60.0, 'contrast_b2_paired_omental_minus_suprafascial_pericellular_po2_gain_mmhg'),
    ('stage0_contrast_b2_porcine_day90_en_bloc_recovery_fraction', 'obs_contrast_b2_porcine_day90_en_bloc_recovery_fraction', 'obs_contrast_b2_porcine_day90_en_bloc_recovery_se', 0.0, 1.0, 'contrast_b2_chronic_porcine_day90_en_bloc_surgical_retrieval_fraction'),
    ('stage0_contrast_b3_iec_omitted_surface_po2_mmhg', 'obs_contrast_b3_iec_omitted_surface_po2_mmhg', 'obs_contrast_b3_iec_omitted_surface_po2_se', 0.0, 120.0, 'contrast_b3_design_z_140um_donor_iec_omitted_loaded_worst_face_po2_mmhg'),
    ('stage0_contrast_b3_iec_omitted_outer_kperf_s_inv', 'obs_contrast_b3_iec_omitted_outer_kperf_s_inv', 'obs_contrast_b3_iec_omitted_outer_kperf_se', 0.0, 1.0, 'contrast_b3_design_z_140um_donor_iec_omitted_host_ingrowth_kperf_s_inv'),
    ('stage0_contrast_b3_iec_omitted_perifusion_gsis_si', 'obs_contrast_b3_iec_omitted_perifusion_gsis_si', 'obs_contrast_b3_iec_omitted_perifusion_gsis_se', 0.1, 50.0, 'contrast_b3_design_z_140um_donor_iec_omitted_loaded_perifusion_gsis_si'),
    # Stage-0 Design B endogenous HLA-C phased CFD + Corrah 2011 DT9/HLA-E specificity control, per-locus segmented BAF, Locus-5 RACE/UMI (>=3,689), and TYK2 PE4max MSI + piecewise repeat-indel endpoints (IIP-0026; DEC-105):
    ('stage0_design_b_hlac_phased_guide_cfd_max', 'obs_design_b_hlac_phased_guide_cfd_max', 'obs_design_b_hlac_phased_guide_cfd_max_se', 0.0, 1.0, 'phased_donor_hlac_exons2_3_max_wt_cfd_across_hla_ab_guides'),
    ('stage0_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio', 'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio', 'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_se', 0.0, 100.0, 'dt9_mfi_fold_over_matched_hlae_only_control_corrah_pmid_21248048'),
    ('stage0_post_edit_segmented_baf_chr19q', 'obs_post_edit_segmented_baf_chr19q', 'obs_post_edit_segmented_baf_chr19q_se', 0.0, 1.0, 'post_edit_phased_50snp_segmented_baf_chr19q_aavs1'),
    ('stage0_post_edit_segmented_baf_chr11p15', 'obs_post_edit_segmented_baf_chr11p15', 'obs_post_edit_segmented_baf_chr11p15_se', 0.0, 1.0, 'post_edit_phased_50snp_segmented_baf_chr11p15_imprinting_region'),
    ('stage0_post_edit_segmented_baf_chr6p21', 'obs_post_edit_segmented_baf_chr6p21', 'obs_post_edit_segmented_baf_chr6p21_se', 0.0, 1.0, 'post_edit_phased_50snp_segmented_baf_chr6p21_mhc_region'),
    ('stage0_design_b_locus5_5prime_race_coord_bp', 'obs_design_b_locus5_5prime_race_coord_bp', 'obs_design_b_locus5_5prime_race_coord_bp_se', 1.0, 10000.0, 'design_b_locus5_mature_5prime_race_donor_coordinate_bp_1based'),
    ('stage0_design_b_locus5_3prime_race_coord_bp', 'obs_design_b_locus5_3prime_race_coord_bp', 'obs_design_b_locus5_3prime_race_coord_bp_se', 1.0, 10000.0, 'design_b_locus5_mature_3prime_race_donor_coordinate_bp_1based'),
    ('stage0_design_b_locus5_race_umi_n_eff', 'obs_design_b_locus5_race_umi_n_eff', 'obs_design_b_locus5_race_umi_n_eff_se', 1.0, 1000000000.0, 'design_b_locus5_race_effective_independent_umi_trials_n_eff_min_3689'),
    ('stage0_design_b_locus5_aberrant_umi_k', 'obs_design_b_locus5_aberrant_umi_k', 'obs_design_b_locus5_aberrant_umi_k_se', 0.0, 1000000.0, 'design_b_locus5_aberrant_chimeric_or_upstream_kozak_umi_count_k'),
    ('stage0_tyk2_pe4max_pentaplex_msi_shifts_k', 'obs_tyk2_pe4max_pentaplex_msi_shifts_k', 'obs_tyk2_pe4max_pentaplex_msi_shifts_k_se', 0.0, 100.0, 'tyk2_pe4max_pentaplex_microsatellite_instability_shifted_markers_k'),
    ('stage0_tyk2_pe4max_wgs_repeat_indels_clone', 'obs_tyk2_pe4max_wgs_repeat_indels_clone', 'obs_tyk2_pe4max_wgs_repeat_indels_clone_se', 0.0, 1000000.0, 'tyk2_pe4max_clonal_wgs_callable_repeat_de_novo_indels_count'),
    ('stage0_tyk2_pe4max_wgs_repeat_indels_ctrl', 'obs_tyk2_pe4max_wgs_repeat_indels_ctrl', 'obs_tyk2_pe4max_wgs_repeat_indels_ctrl_se', 0.0, 1000000.0, 'tyk2_pe4max_matched_control_wgs_callable_repeat_de_novo_indels_count'),
    ('stage0_tyk2_pe4max_coding_or_splice_de_novo_indels', 'obs_tyk2_pe4max_coding_or_splice_de_novo_indels', 'obs_tyk2_pe4max_coding_or_splice_de_novo_indels_se', 0.0, 10000.0, 'tyk2_pe4max_coding_or_canonical_splice_de_novo_indels_count'),
    # Stage-0 arm-separated Gate 0d, paired intact-plus-increment CD16+ ADCC, measured-potency DEC-084 mixture, and 72h/d28/d84 dual-modality IAPP endpoints (IIP-0028; DEC-106):
    ('stage0_gate0d_arm1_unconditioned_cd4_si', 'obs_arm1_unconditioned_cd4_si', 'obs_arm1_unconditioned_cd4_si_se', 0.1, 100.0, 'arm1_never_conditioned_indirect_apc_dq2_dq8_cd4_stimulation_index'),
    ('stage0_gate0d_arm1_unconditioned_annexin_v_pct', 'obs_arm1_unconditioned_annexin_v_pct', 'obs_arm1_unconditioned_annexin_v_se', 0.0, 100.0, 'arm1_never_conditioned_transwell_beta_annexin_v_positive_percentage'),
    ('stage0_gate0d_parallel_arm3_control_cd4_apoptosis_pct', 'obs_parallel_arm3_control_cd4_apoptosis_pct', 'obs_parallel_arm3_control_cd4_apoptosis_se', 0.0, 100.0, 'contemporaneous_parallel_arm3_active_fasl_control_activated_cd4_apoptosis_percentage'),
    ('stage0_gate0d_post_depletion_input_normalized_cpep_ratio', 'obs_post_depletion_input_normalized_cpep_ratio', 'obs_post_depletion_input_normalized_cpep_ratio_se', 0.01, 5.0, 'arm2_post_depletion_input_normalized_stimulated_c_peptide_ratio_vs_sham'),
    ('stage0_fcr_cd16_intact_serum_total_lysis_pct', 'obs_fcr_cd16_intact_serum_total_lysis_pct', 'obs_fcr_cd16_intact_serum_total_lysis_se', 0.0, 100.0, 'intact_recipient_serum_fcgr3a_cd16_total_cellular_lysis_percentage'),
    ('stage0_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct', 'obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct', 'obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_se', -20.0, 100.0, 'paired_igg_attributable_fcgr3a_cd16_net_adcc_lysis_increment_vs_igg_depleted_serum_percentage'),
    ('stage0_adcc_isotype_background_lysis_pct', 'obs_adcc_isotype_background_lysis_pct', 'obs_adcc_isotype_background_lysis_se', 0.0, 100.0, 'matrix_matched_isotype_control_cd16_background_cellular_lysis_percentage'),
    ('stage0_igg_depletion_efficiency_pct', 'obs_igg_depletion_efficiency_pct', 'obs_igg_depletion_efficiency_se', 0.0, 100.0, 'paired_recipient_serum_protein_g_igg_depletion_efficiency_percentage'),
    ('stage0_dec084_starting_potency_ratio_out_b_over_out_a', 'obs_dec084_starting_potency_ratio_out_b_over_out_a', 'obs_dec084_starting_potency_ratio_out_b_over_out_a_se', 0.1, 5.0, 'measured_starting_per_beta_stimulated_c_peptide_output_ratio_candidate_over_reference'),
    ('stage0_dec084_ocr_ratio_b_over_a', 'obs_dec084_ocr_ratio_b_over_a', 'obs_dec084_ocr_ratio_b_over_a_se', 0.1, 5.0, 'equal_output_loaded_whole_graft_ocr_ratio_candidate_over_60_30_10_reference'),
    ('stage0_dec084_cpep_auc_ratio_b_over_a', 'obs_dec084_cpep_auc_ratio_b_over_a', 'obs_dec084_cpep_auc_ratio_b_over_a_se', 0.1, 5.0, 'equal_output_loaded_stimulated_c_peptide_auc_ratio_candidate_over_60_30_10_reference'),
    ('stage0_dec084_gcg_auc_ratio_b_over_a', 'obs_dec084_gcg_auc_ratio_b_over_a', 'obs_dec084_gcg_auc_ratio_b_over_a_se', 0.1, 5.0, 'equal_output_loaded_glucopenic_glucagon_auc_ratio_candidate_over_60_30_10_reference'),
    ('stage0_iapp_72h_immunocapture_oligomer_elisa_fold', 'obs_iapp_72h_immunocapture_oligomer_elisa_fold', 'obs_iapp_72h_immunocapture_oligomer_elisa_se', 0.0, 50.0, 'lot_scope_72h_iapp_immunocapture_oligomer_elisa_fold_over_age_matched_ref'),
    ('stage0_iapp_72h_colocalized_fibril_area_pct', 'obs_iapp_72h_colocalized_fibril_area_pct', 'obs_iapp_72h_colocalized_fibril_area_se', 0.0, 100.0, 'lot_scope_72h_anti_iapp_colocalized_fibril_area_percentage'),
    ('stage0_iapp_d28_immunocapture_oligomer_elisa_fold', 'obs_iapp_d28_immunocapture_oligomer_elisa_fold', 'obs_iapp_d28_immunocapture_oligomer_elisa_se', 0.0, 50.0, 'lot_scope_day28_iapp_immunocapture_oligomer_elisa_fold_over_age_matched_ref'),
    ('stage0_iapp_d28_colocalized_fibril_area_pct', 'obs_iapp_d28_colocalized_fibril_area_pct', 'obs_iapp_d28_colocalized_fibril_area_se', 0.0, 100.0, 'lot_scope_day28_anti_iapp_colocalized_fibril_area_percentage'),
    ('stage0_iapp_d84_immunocapture_oligomer_elisa_fold', 'obs_iapp_d84_immunocapture_oligomer_elisa_fold', 'obs_iapp_d84_immunocapture_oligomer_elisa_se', 0.0, 50.0, 'terminal_sister_device_day84_iapp_immunocapture_oligomer_elisa_fold_over_age_matched_ref'),
    ('stage0_iapp_d84_colocalized_fibril_area_pct', 'obs_iapp_d84_colocalized_fibril_area_pct', 'obs_iapp_d84_colocalized_fibril_area_se', 0.0, 100.0, 'terminal_sister_device_day84_anti_iapp_colocalized_fibril_area_percentage'),
    # Stage-0 resting (0 h) + interferon-stimulated (24 h) CXCL10 promoter-of-origin 5'-RACE TSS attribution & conditional BIO-06 AAV-ITR RNA-origin endpoints (IIP-0027; DEC-107):
    ('stage0_cxcl10_resting_race_umi_n_eff', 'obs_cxcl10_resting_race_umi_n_eff', 'obs_cxcl10_resting_race_umi_n_eff_se', 0.0, 1000000000.0, 'cxcl10_resting_0h_5prime_race_recovery_adjusted_effective_umi_n_eff'),
    ('stage0_cxcl10_resting_unexpected_tss_umi_k', 'obs_cxcl10_resting_unexpected_tss_umi_k', 'obs_cxcl10_resting_unexpected_tss_umi_k_se', 0.0, 1000000.0, 'cxcl10_resting_0h_unexpected_internal_promoter_tss_umi_count_k'),
    ('stage0_cxcl10_resting_unresolved_origin_umi_k', 'obs_cxcl10_resting_unresolved_origin_umi_k', 'obs_cxcl10_resting_unresolved_origin_umi_k_se', 0.0, 1000000.0, 'cxcl10_resting_0h_unresolved_tss_origin_umi_count_k'),
    ('stage0_cxcl10_ifng24h_race_umi_n_eff', 'obs_cxcl10_ifng24h_race_umi_n_eff', 'obs_cxcl10_ifng24h_race_umi_n_eff_se', 0.0, 1000000000.0, 'cxcl10_ifng24h_5prime_race_recovery_adjusted_effective_umi_n_eff'),
    ('stage0_cxcl10_ifng24h_unexpected_tss_umi_k', 'obs_cxcl10_ifng24h_unexpected_tss_umi_k', 'obs_cxcl10_ifng24h_unexpected_tss_umi_k_se', 0.0, 1000000.0, 'cxcl10_ifng24h_unexpected_internal_promoter_tss_umi_count_k'),
    ('stage0_cxcl10_ifng24h_unresolved_origin_umi_k', 'obs_cxcl10_ifng24h_unresolved_origin_umi_k', 'obs_cxcl10_ifng24h_unresolved_origin_umi_k_se', 0.0, 1000000.0, 'cxcl10_ifng24h_unresolved_tss_origin_umi_count_k'),
    ('stage0_bio06_retained_aav_itr_dna_pct', 'obs_bio06_retained_aav_itr_dna_pct', 'obs_bio06_retained_aav_itr_dna_pct_se', 0.0, 100.0, 'bio06_locus_junction_retained_aav_itr_dna_allele_percentage'),
    ('stage0_bio06_itr_origin_rna_umi_k', 'obs_bio06_itr_origin_rna_umi_k', 'obs_bio06_itr_origin_rna_umi_k_se', 0.0, 1000000.0, 'bio06_conditional_aav_itr_origin_transcript_5prime_race_umi_count_k'),
    # Stage-0 ALK5i-free two-clock maturation, two-pulse Readout 2 with adult-islet reference floor, inert-tracer tau_chamber, and >138.3 IU Gate 0S atomic escalation endpoints (IIP-0029; DEC-108):
    ('stage0_stage7_total_maturation_days', 'obs_stage7_total_maturation_days', 'obs_stage7_total_maturation_days_se', 0.0, 90.0, 'total_stage7_culture_duration_days_from_stage7_start_utc'),
    ('stage0_stage7_days_post_sort_reaggregation', 'obs_stage7_days_post_sort_reaggregation', 'obs_stage7_days_post_sort_reaggregation_se', 0.0, 90.0, 'days_elapsed_from_sort_reaggregation_utc_to_assay_harvest_utc'),
    ('stage0_stage7_eta_sort_viable_recovery', 'obs_stage7_eta_sort_viable_recovery', 'obs_stage7_eta_sort_viable_recovery_se', 0.0, 1.0, 'sorting_and_reaggregation_viable_recovery_fraction_on_common_denominator'),
    ('stage0_readout2_inert_tracer_tau_chamber_min', 'obs_readout2_inert_tracer_tau_chamber_min', 'obs_readout2_inert_tracer_tau_chamber_se', 0.0, 30.0, 'direct_fluorescent_inert_tracer_perifusion_chamber_mixing_time_constant_minutes'),
    ('stage0_readout2_pulse2_net_increment_z2', 'obs_readout2_pulse2_net_increment_z2', 'obs_readout2_pulse2_net_increment_z2_se', -50.0, 500.0, 'recovery_weighted_second_pulse_net_increment_z2_on_registered_unit_and_common_denominator'),
    ('stage0_readout2_adult_islet_ref_z2', 'obs_readout2_adult_islet_ref_z2', 'obs_readout2_adult_islet_ref_z2_se', -50.0, 500.0, 'contemporaneous_same_assay_primary_adult_human_islet_second_pulse_net_increment_ref_z2'),
    ('stage0_readout2_low2_frac_peak1', 'obs_readout2_low2_frac_peak1', 'obs_readout2_low2_frac_peak1_se', 0.0, 5.0, 'second_low_glucose_washout_secretion_as_fraction_of_first_phase_peak1'),
    ('stage0_readout2_low2_over_low1_ratio', 'obs_readout2_low2_over_low1_ratio', 'obs_readout2_low2_over_low1_ratio_se', 0.0, 20.0, 'second_low_glucose_washout_over_first_low_glucose_baseline_ratio'),
    ('stage0_stage7_s8_fieller_z2_diff', 'obs_stage7_s8_fieller_z2_diff', 'obs_stage7_s8_fieller_z2_diff_se', -500.0, 500.0, 'within_recipe_rule_s8_paired_fieller_contrast_z2_ext_minus_0p85_z2_early'),
    ('stage0_stage7_stored_mature_insulin_iu_implanted', 'obs_stage7_stored_mature_insulin_iu_implanted', 'obs_stage7_stored_mature_insulin_iu_implanted_se', 0.0, 1000.0, 'specific_mature_insulin_inventory_iu_at_actual_implanted_dose'),
    ('stage0_stage7_stored_mature_insulin_max_lot_iu', 'obs_stage7_stored_mature_insulin_max_lot_iu', 'obs_stage7_stored_mature_insulin_max_lot_iu_se', 0.0, 1000.0, 'maximum_individual_lot_specific_mature_insulin_inventory_iu_at_actual_implanted_dose'),
    ('stage0_stage7_proinsulin_to_cpeptide_molar_ratio', 'obs_stage7_proinsulin_to_cpeptide_molar_ratio', 'obs_stage7_proinsulin_to_cpeptide_molar_ratio_se', 0.0, 5.0, 'exploratory_stage0_in_vitro_intact_proinsulin_to_c_peptide_content_molar_ratio'),
    ('stage0_stage7_pre_iec_zm447439_residual_uM', 'obs_stage7_pre_iec_zm447439_residual_uM', 'obs_stage7_pre_iec_zm447439_residual_se', 0.0, 10.0, 'analytical_lc_ms_ms_residual_zm447439_concentration_micromolar_prior_to_iec_coassembly'),
    # Stage-0 Locus-2 CIITA Mode-B v3 complement/ADCC/efferocytosis/pericyte/AAV6 Rule S11 endpoints (IIP-0031; DEC-110):
    ('stage0_l2_v3_ap_c5b9_pct', 'obs_l2_v3_ap_c5b9_pct', 'obs_l2_v3_ap_c5b9_se', 0.0, 100.0, 'post_cuff_alternative_pathway_20pct_nhs_c5b9_deposition_pct'),
    ('stage0_l2_v3_classical_c4b2a_lysis_pct', 'obs_l2_v3_classical_c4b2a_lysis_pct', 'obs_l2_v3_classical_c4b2a_lysis_se', 0.0, 100.0, 'post_cuff_classical_c1q_sensitized_c4b2a_cdc_lysis_pct'),
    ('stage0_l2_v3_lectin_c4b2a_lysis_pct', 'obs_l2_v3_lectin_c4b2a_lysis_pct', 'obs_l2_v3_lectin_c4b2a_lysis_se', 0.0, 100.0, 'post_cuff_lectin_mbl_mannan_c4b2a_cdc_lysis_pct'),
    ('stage0_l2_v3_overall_20pct_nhs_lysis_pct', 'obs_l2_v3_overall_20pct_nhs_lysis_pct', 'obs_l2_v3_overall_20pct_nhs_lysis_se', 0.0, 100.0, 'post_cuff_playbook_overall_20pct_nhs_specific_lysis_pct'),
    ('stage0_l2_v3_intact_serum_adcc_lysis_pct', 'obs_l2_v3_intact_serum_adcc_lysis_pct', 'obs_l2_v3_intact_serum_adcc_lysis_se', 0.0, 100.0, 'post_cuff_gy023_intact_serum_cd16_adcc_total_lysis_pct'),
    ('stage0_l2_v3_igg_increment_adcc_lysis_pct', 'obs_l2_v3_igg_increment_adcc_lysis_pct', 'obs_l2_v3_igg_increment_adcc_lysis_se', -100.0, 100.0, 'post_cuff_gy023_paired_igg_attributable_cd16_adcc_increment_pct'),
    ('stage0_l2_v3_exact_stack_effero_abs_pct', 'obs_l2_v3_exact_stack_effero_abs_pct', 'obs_l2_v3_exact_stack_effero_abs_se', 0.0, 100.0, 'exact_stack_apoptotic_corpse_efferocytosis_absolute_clearance_pct'),
    ('stage0_l2_v3_exact_stack_effero_pct_of_wt', 'obs_l2_v3_exact_stack_effero_pct_of_wt', 'obs_l2_v3_exact_stack_effero_pct_of_wt_se', 0.0, 200.0, 'exact_stack_apoptotic_corpse_efferocytosis_pct_of_wt'),
    ('stage0_l2_v3_pericyte_coverage_pct', 'obs_l2_v3_pericyte_coverage_pct', 'obs_l2_v3_pericyte_coverage_se', 0.0, 100.0, 'gate0g_vascular_pericyte_ng2_pdgfrb_coverage_pct'),
    ('stage0_l2_v3_aav6_full_length_monomer_pct', 'obs_l2_v3_aav6_full_length_monomer_pct', 'obs_l2_v3_aav6_full_length_monomer_se', 0.0, 100.0, 'reconciled_aav6_alkaline_gel_smrt_ont_full_length_monomer_pct'),
    ('stage0_l2_v3_aav6_truncated_species_pct', 'obs_l2_v3_aav6_truncated_species_pct', 'obs_l2_v3_aav6_truncated_species_se', 0.0, 100.0, 'reconciled_aav6_truncated_genome_species_pct'),
    ('stage0_l2_v3_aav6_dimeric_or_intermediate_pct', 'obs_l2_v3_aav6_dimeric_or_intermediate_pct', 'obs_l2_v3_aav6_dimeric_or_intermediate_se', 0.0, 100.0, 'reconciled_aav6_dimeric_and_intermediate_species_pct'),
]

STAGE0_S10_V3_MEAN_KEYS: frozenset[str] = frozenset({
    'obs_hlae_sct_nkg2a_kd_uM',
    'obs_hlae_sct_nkg2c_kd_uM',
    'obs_wt_sct_nkg2a_kd_uM',
    'obs_nkg2a_pos_resting_nk_lysis_pct',
    'obs_nkg2a_pos_ifng24h_nk_lysis_pct',
    'obs_nkg2a_pos_resting_cd107a_pct',
    'obs_nkg2a_pos_ifng24h_cd107a_pct',
    'obs_nkg2c_kir_pos_resting_lysis_pct',
    'obs_nkg2c_kir_pos_ifng24h_lysis_pct',
    'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_pct',
    'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_pct',
    'obs_hcmv_pooled_3subset_resting_nk_lysis_pct',
    'obs_hcmv_pooled_3subset_ifng24h_nk_lysis_pct',
    'obs_nkg2c_cd107a_degranulation_fold_vs_b2m_ko',
    'obs_cmv_ul40_hlae_cd8_lysis_pct',
    'obs_hcmv_pos_donor_count',
    'obs_s10_nk_failing_donor_count',
    'obs_ul40_reactive_donor_count',
    'obs_ul40_ctl_failing_donor_count',
    'obs_hlae_cd8aa_kd_uM',
    'obs_hlae_cd8ab_kd_uM',
    'obs_cd8_treg_suppression_retention_pct_of_wt',
    'obs_cd8_treg_suppression_delta_vs_wt_pp',
    'obs_cmv_ul40_hlae_cd8_lysis_delta_vs_wt_pp',
    'obs_nonself_dq2_dq8_indirect_apc_si',
    'obs_nonself_host_class1_cross_pres_lysis_pct',
})

STAGE0_S10_V3_METADATA_KEYS: tuple[str, ...] = (
    's10_construct_id',
    's10_construct_donor_sha256',
    'current_lot_id',
    's10_qualified_lot_id',
    's10_hlae_floor_mode',
)

STAGE0_S10_V3_BOOLEAN_KEYS: tuple[str, ...] = (
    's10_is_replacement_challenger_lot',
    's10_race_5p_3p_junction_verified',
    'hla_a2_d227k_t228a_positive_control_spr_verified',
)

ALLOWED_S10_HLAE_FLOOR_MODES: frozenset[str] = frozenset({
    'GATE_0A_SCREENING_15K',
    'LOT_RELEASE_25K',
})

ALLOWED_B2_TYK2_RECIPIENT_GENOTYPE_STRATA: frozenset[str] = frozenset({
    'G_G',
    'G_C',
    'C_C',
})

ALLOWED_B2_TYK2_HLAC_ALLELES: frozenset[str] = frozenset({
    'HLA-C*01:02',
    'HLA-C*02:02',
    'HLA-C*03:04',
    'HLA-C*04:01',
    'HLA-C*05:01',
    'HLA-C*06:02',
    'HLA-C*07:01',
    'HLA-C*07:02',
    'HLA-C*08:01',
    'HLA-C*12:02',
    'HLA-C*12:03',
    'HLA-C*14:02',
    'HLA-C*15:02',
    'HLA-C*16:01',
})

ALLOWED_S6_CLUSTER_SIZE_ARMS: frozenset[str] = frozenset({
    'DXB_100UM_1500IEQ',
    'DOSE_MATCHED_80UM_1500IEQ',
    'DERATED_80UM_768IEQ',
})

ALLOWED_P05_BARRIER_SLAB_GEOMETRY_MODES: frozenset[str] = frozenset({
    'EXPANDED_390UM_2X125UM_BANDS',
    'FIXED_350UM_2X105UM_BANDS',
})

ALLOWED_DESIGN_Z_DIVIDER_ARMS: frozenset[str] = frozenset({
    'PLAIN_HYDROGEL_140UM',
    'EXPLORATORY_18PCT_DIVIDER_140UM',
})

# Stage-0 Gate P05 human catalase species lock and Mode-B streptavidin-free / DPR-maleimide humoral cuff qualification enums (IIP-0020; DEC-099; IIP-0025; DEC-104):
ALLOWED_P05_CATALASE_UNIPROT_ACCESSIONS: frozenset[str] = frozenset({
    'P04040',
    'P00432',
    'P00432_BOVINE_DISQUALIFIED',
    'P55303_FUNGAL_DISQUALIFIED',
    'UNSPECIFIED_XENO_DISQUALIFIED',
})
ALLOWED_CONTRAST_B1_GEOMETRY_ARMS: frozenset[str] = frozenset({
    'PRIMARY_EXPANDED_390UM_2X125UM',
    'SECONDARY_FIXED_350UM_2X105UM',
})
ALLOWED_CONTRAST_B2_ARMS: frozenset[str] = frozenset({
    'PAIRED_HUMAN_EXVIVO_AND_PORCINE_DAY30_90',
    'HUMAN_EXVIVO_ONLY_INCOMPLETE',
    'PORCINE_RETRIEVAL_ONLY_INCOMPLETE',
})
ALLOWED_CONTRAST_B3_DONOR_IEC_ARMS: frozenset[str] = frozenset({
    'WITH_DONOR_IEC_9E8',
    'OMIT_DONOR_IEC_8E8_CHALLENGER',
})
ALLOWED_TYK2_PE_MODALITIES: frozenset[str] = frozenset({
    'PE2MAX_LEAD',
    'PE4MAX_MSI_GATED',
})
ALLOWED_IIP0026_SCHEMA_ERAS: frozenset[str] = frozenset({
    'IIP0026_V1',
    'LEGACY_PRE_IIP0026',
})
ALLOWED_CUFF_BIOCONJUGATION_ARMS: frozenset[str] = frozenset({
    'CANONICAL_SA_BIOTIN_PEG4MAL_BASELINE',
    'DPR_HYDROLYZED_CFH_THBD_H1A_ADIPOQ_HEXAFASL',
    'DPR_HYDROLYZED_CFH_THBD_T0_TNC_FASL',
    'DPR_HYDROLYZED_CFH_THBD_FASL_FREE',
    'ZWITTERIONIC_ALGINATE_PEG_FREE_CUFF',
})
ALLOWED_CUFF_SERUM_CHALLENGE_STRATA: frozenset[str] = frozenset({
    'NORMAL_DONOR_SERA',
    'HIGH_ANTI_PEG_IGG_GT_500NG_ML',
    'ANTI_SA_POSITIVE_HUMAN_SERA',
    'DAY28_SA_FASL_PRIMED_SERA',
})

# Stage-0 120-h lineage-partitioned metabolic-rest OCR and 4-arm Nec-1s x Ferrostatin-1 conditioning / NQO1 / pellet-washout enums (IIP-0021; DEC-100):
ALLOWED_REST_OCR_COUPON_ARMS: frozenset[str] = frozenset({
    'VEHICLE_791P5MG_CAO2_BASELINE',
    'ENCAPSULATED_KATP_MICRODEPOT_CHALLENGER',
    'UNENCAPSULATED_FREE_DRUG_REJECTED',
})
ALLOWED_CONDITIONING_4ARM_ARMS: frozenset[str] = frozenset({
    'VEHICLE_CONTROL',
    'NEC1S_10UM_SELECTIVE',
    'FERROSTATIN1_2UM',
    'NEC1S_10UM_PLUS_FERROSTATIN1_2UM',
    'NEC1S_100UM_BENCHMARK_COMPARATOR',
    'NEC1_100UM_ANTIFERROPTOTIC_COMPARATOR',
})
ALLOWED_CONDITIONING_WASHOUT_ANALYTE_SCOPES: frozenset[str] = frozenset({
    'SUPERNATANT_AND_CELL_PELLET_PARENT_PLUS_METABOLITES',
    'AQUEOUS_SUPERNATANT_ONLY_INSUFFICIENT',
})

# Stage-0 Track-2 NK rescue split (HLA-C B'-SCT KIR vs held HLA-G1 LILRB1 challenger; IIP-0022; DEC-101):
ALLOWED_TRACK2_RESCUE_ROUTES: frozenset[str] = frozenset({
    'TRACK1_MHC_NULL_NO_TRACK2_SWITCH',
    'TRACK2_HLAC_B_SCT_KIR_RESCUE',
    'TRACK2_HLAG1_LILRB1_CHALLENGER_HELD',
})

# Stage-0 Design Y wildcard (<= 20% firewall) integrase, silencing, and payload arm enums (IIP-0023; DEC-102):
ALLOWED_DESIGN_Y_INTEGRASE_MODALITIES: frozenset[str] = frozenset({
    'EEPASSIGE_B2M_ATTP_GENOME_ATTB_MINICIRCLE',
    'LEGACY_ATTB_GENOME_ATTP_PLASMID_REJECTED',
})
ALLOWED_DESIGN_Y_SILENCING_MODALITIES: frozenset[str] = frozenset({
    'TADCBED_CIITA_W334X_CXCL10_M1I',
    'TADCBED_CIITA_Q53X_CXCL10_M1I_COMPARATOR',
    'CRISPROFF_EPIGENETIC_COMPARATOR',
})
ALLOWED_DESIGN_Y_PAYLOAD_ARMS: frozenset[str] = frozenset({
    'DESIGN_Y_CORE_CD200_IEC_11066BP',
    'DESIGN_Y_CD59_PLUS_CD200_IEC_11859BP',
    'DESIGN_Y_CD39_IEC_COMPARATOR_12555BP',
})

STAGE0_S10_V3_CONSTRUCT_REGISTRY: dict[str, dict[str, object]] = {
    'S10_CTRL_WT_VMAPRTLVL_ASTAR_4143': {
        'construct_id': 'S10_CTRL_WT_VMAPRTLVL_ASTAR_4143',
        'donor_sha256': 'aff2ec0d862b9a99f54518465221106444e3ce0f8b055b29dd0f41faac262225',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTLVL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'BASELINE_PARENTAL_CONTROL',
        'promotion_eligibility': 'BASELINE_REFERENCE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': False,
    },
    'S10_CTRL_WT_VMAPRTLVL_ARM1B_3974': {
        'construct_id': 'S10_CTRL_WT_VMAPRTLVL_ARM1B_3974',
        'donor_sha256': '3d516e2f0fa4408ffe60aa111183ad4e57f4a669c77cba2b2aad15e322d17964',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTLVL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'BASELINE_PARENTAL_CONTROL',
        'promotion_eligibility': 'BASELINE_REFERENCE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': False,
    },
    'S10_CTRL_POS_VMAPRTLFL_ASTAR_4143': {
        'construct_id': 'S10_CTRL_POS_VMAPRTLFL_ASTAR_4143',
        'donor_sha256': 'ba47b154207763f44c72baef3b27f35beecc84cfd37ae6be9eb009605202955b',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTLFL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'POSITIVE_NKG2C_AGONIST_RATIO_PARADOX_CONTROL',
        'promotion_eligibility': 'CALIBRATION_CONTROL_NEVER_PROMOTE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CTRL_POS_VMAPRTLFL_ARM1B_3974': {
        'construct_id': 'S10_CTRL_POS_VMAPRTLFL_ARM1B_3974',
        'donor_sha256': '6e6288704e9cdc243b1d97ad6c3d7404bd03579477cf6ef796137e853423cd2d',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTLFL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'POSITIVE_NKG2C_AGONIST_RATIO_PARADOX_CONTROL',
        'promotion_eligibility': 'CALIBRATION_CONTROL_NEVER_PROMOTE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CTRL_NEG_QMRPVSRVL_ASTAR_4143': {
        'construct_id': 'S10_CTRL_NEG_QMRPVSRVL_ASTAR_4143',
        'donor_sha256': 'a6fd61149f07a9c87e69526729448e44679e949521f8bc7dc4036d0b16fa6489',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'QMRPVSRVL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'GROOVE_STABILIZING_RECEPTOR_PROTECTION_NEGATIVE_CONTROL',
        'promotion_eligibility': 'CALIBRATION_CONTROL_NEVER_PROMOTE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CTRL_NEG_QMRPVSRVL_ARM1B_3974': {
        'construct_id': 'S10_CTRL_NEG_QMRPVSRVL_ARM1B_3974',
        'donor_sha256': '9e064dbfa7ca1ba55b7d7e55e3cb8b36869e7abae310e4d9ee7d51a95ee305d7',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'QMRPVSRVL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'GROOVE_STABILIZING_RECEPTOR_PROTECTION_NEGATIVE_CONTROL',
        'promotion_eligibility': 'CALIBRATION_CONTROL_NEVER_PROMOTE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143': {
        'construct_id': 'S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143',
        'donor_sha256': '797cf4472a096797cb46490ef133100866827de88420b02b5e6ac15195715765',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTVLL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'EXPLORATORY_WEAK_RECOGNITION_DUAL_ATTENUATED_SP6B_COMPARATOR',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL1_SP6B_VMAPRTVLL_ARM1B_3974': {
        'construct_id': 'S10_CHAL1_SP6B_VMAPRTVLL_ARM1B_3974',
        'donor_sha256': '663cb4bcdc5f0473f6b310ecd60d1a7d218cf4843162d5fbe9497b89e8845791',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTVLL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'EXPLORATORY_WEAK_RECOGNITION_DUAL_ATTENUATED_SP6B_COMPARATOR',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL2_SP1B_VTAPRTLLL_ASTAR_4143': {
        'construct_id': 'S10_CHAL2_SP1B_VTAPRTLLL_ASTAR_4143',
        'donor_sha256': '180ad50d74b035752ddd4fda4f370419db07c4bc470767d793aefbe68b1b992e',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VTAPRTLLL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'COVALENT_MINUS21T_BYPASS_WEAK_ENHANCER_CALIBRATION_COMPARATOR',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL2_SP1B_VTAPRTLLL_ARM1B_3974': {
        'construct_id': 'S10_CHAL2_SP1B_VTAPRTLLL_ARM1B_3974',
        'donor_sha256': 'b136a00e566e0dc7069536d734f76677904b224f518540ce616943d7d39de216',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VTAPRTLLL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'COVALENT_MINUS21T_BYPASS_WEAK_ENHANCER_CALIBRATION_COMPARATOR',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL3_SP5A_VMPPRTLLL_ASTAR_4143': {
        'construct_id': 'S10_CHAL3_SP5A_VMPPRTLLL_ASTAR_4143',
        'donor_sha256': '9e2de8f649d3d7b4b33bb37c6e10cb21b7ed587dba110243ba6e7a2dee53ac7b',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMPPRTLLL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'FUNCTIONAL_P3PRO_UL40_COHORT_NULL_WEAK_ENHANCER_COMPARATOR',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL3_SP5A_VMPPRTLLL_ARM1B_3974': {
        'construct_id': 'S10_CHAL3_SP5A_VMPPRTLLL_ARM1B_3974',
        'donor_sha256': 'f286ca0cf0cb9b399bd80b334022d0b57fe57594d750730577cfb6d378a8aae0',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMPPRTLLL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'FUNCTIONAL_P3PRO_UL40_COHORT_NULL_WEAK_ENHANCER_COMPARATOR',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL4_P8A_VMAPRTLAL_ASTAR_4143': {
        'construct_id': 'S10_CHAL4_P8A_VMAPRTLAL_ASTAR_4143',
        'donor_sha256': '334eddd5d77553825c39a17a451117b68433933807f283054bc15aa5c050a9e2',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTLAL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'SYNTHETIC_P8ALA_UL40_NULL_NKG2A_SELECTIVE_MECHANISTIC_PROBE',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL4_P8A_VMAPRTLAL_ARM1B_3974': {
        'construct_id': 'S10_CHAL4_P8A_VMAPRTLAL_ARM1B_3974',
        'donor_sha256': 'c76c54199c495281461cdc3237021313b9e15077639f5dbc6e4e385cf1568811',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTLAL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'SYNTHETIC_P8ALA_UL40_NULL_NKG2A_SELECTIVE_MECHANISTIC_PROBE',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL5A_WT_VL_D227K_T228A_ASTAR_4143': {
        'construct_id': 'S10_CHAL5A_WT_VL_D227K_T228A_ASTAR_4143',
        'donor_sha256': '8d9286c847228b19dde95e3113ea817b5622394646276c9c48b24edd3e9e0902',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTLVL',
        'alpha3_allele': 'MUT_D227K_T228A',
        'role': 'SINGLE_FACTOR_ALPHA3_CD8_CORECEPTOR_NULL_CHALLENGER',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL5A_WT_VL_D227K_T228A_ARM1B_3974': {
        'construct_id': 'S10_CHAL5A_WT_VL_D227K_T228A_ARM1B_3974',
        'donor_sha256': '80f60a9de0409438775df4a49bb6862d2d7cb0f87f34c6df5611210e92586ebe',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTLVL',
        'alpha3_allele': 'MUT_D227K_T228A',
        'role': 'SINGLE_FACTOR_ALPHA3_CD8_CORECEPTOR_NULL_CHALLENGER',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL5B_SP6B_VMAPRTVLL_D227K_T228A_ASTAR_4143': {
        'construct_id': 'S10_CHAL5B_SP6B_VMAPRTVLL_D227K_T228A_ASTAR_4143',
        'donor_sha256': 'dc87e10e98b4599486d34aa4e76a306072670ee2a67c7d8a9c5d8f5d0584003b',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTVLL',
        'alpha3_allele': 'MUT_D227K_T228A',
        'role': 'FACTORIAL_2X2_SP6B_PLUS_ALPHA3_CD8_NULL_CHALLENGER',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CHAL5B_SP6B_VMAPRTVLL_D227K_T228A_ARM1B_3974': {
        'construct_id': 'S10_CHAL5B_SP6B_VMAPRTVLL_D227K_T228A_ARM1B_3974',
        'donor_sha256': '2b2f2875f7a9bb8392d92e592667b4850edb89c77c3469b0ea1737317385a31b',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTVLL',
        'alpha3_allele': 'MUT_D227K_T228A',
        'role': 'FACTORIAL_2X2_SP6B_PLUS_ALPHA3_CD8_NULL_CHALLENGER',
        'promotion_eligibility': 'ELIGIBLE_MODE_B_CHALLENGER',
        'is_replacement_challenger': True,
        'requires_nonself_peptide_gate': True,
    },
    'S10_CTRL_HARDWT_VMAPRTLVL_ASTAR_4143': {
        'construct_id': 'S10_CTRL_HARDWT_VMAPRTLVL_ASTAR_4143',
        'donor_sha256': 'b5e0bee49a1c59d9368dbe083132485cd8da666fd6e311d22abcd103c589f1e9',
        'parent_construct_id': 'S3_L1_ASTAR_STOP_ONLY_4143',
        'nonamer_aa': 'VMAPRTLVL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'SYNONYMOUS_CTC_SPLICE_HARDENED_WT_SINGLE_FACTOR_CONTROL',
        'promotion_eligibility': 'BASELINE_REFERENCE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': False,
    },
    'S10_CTRL_HARDWT_VMAPRTLVL_ARM1B_3974': {
        'construct_id': 'S10_CTRL_HARDWT_VMAPRTLVL_ARM1B_3974',
        'donor_sha256': 'e3b3f04321dceccf0fb013d7572a7bbccbf9aa4ae10dbe47a14d659b90f0652f',
        'parent_construct_id': 'L1_ARM1B_STOP_ONLY_S3_PAM_SHIELDED_3974',
        'nonamer_aa': 'VMAPRTLVL',
        'alpha3_allele': 'WT_D227_T228',
        'role': 'SYNONYMOUS_CTC_SPLICE_HARDENED_WT_SINGLE_FACTOR_CONTROL',
        'promotion_eligibility': 'BASELINE_REFERENCE',
        'is_replacement_challenger': False,
        'requires_nonself_peptide_gate': False,
    },
}

STICKY_S9_PRIOR_VARIANT_SAFETY_STATUSES: frozenset[str] = frozenset({
    'HARD_SAFETY_STOP_REJECTED',
    'HARD_SAFETY_STOP_STICKY_PRIOR_S9_RETAINED',
    'REJECTED_BY_S9_HARD_SAFETY_STOP',
    'HARD_SAFETY_STOP_MCB_CLONOGENIC_SURVIVOR_DETECTED',
})

ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES: frozenset[str] = STICKY_S9_PRIOR_VARIANT_SAFETY_STATUSES | frozenset({
    'CLEARED',
    'FULL_5ENDPOINT_S9_PANEL_CLEARED',
    'CORE_3ENDPOINT_GENOMIC_SUICIDE_CLEARED_ER_GATE0E_PENDING',
    'HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL',
    'HOLD_INCOMPLETE_5ENDPOINT_S9_PANEL',
    'HOLD_INCOMPLETE_CORE_3ENDPOINT_S9_PANEL',
    'HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM',
    'HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM',
    'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING',
    'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM',
    'HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM',
    'HOLD_HAZARD_STATE_PROLIF_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM',
    'HOLD_AP1903_STEP1_EFFECTIVE_TRIALS_BELOW_29957323_MINIMUM',
    'HOLD_VASCULARIZED_EXPLANT_EN_BLOC_RECOVERY_UNVERIFIED',
})

ALLOWED_STICKY_S9_RESET_REASON_CODES: frozenset[str] = frozenset({
    'INSTRUMENT_CONTAMINATION_VERIFIED_NEGATIVE_RETEST',
    'SAMPLE_SWAP_STR_CONFIRMED_RETEST',
    'PIPELINE_BIOINFORMATIC_ARTIFACT_RAW_BAM_REVERIFIED',
    'ASSAY_ARTIFACT_CONFIRMED_BY_REPEAT_DUPLEX_SEQ',
})

ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES: frozenset[str] = frozenset({
    'FRACTION_HALF_EFFECT',
    'PP_HALF_EFFECT',
    'PP_FULL_CONTRAST',
    'DID_PP',
})

ALLOWED_ADJUNCT_CONDITIONING_HISTORY_ARMS: frozenset[str] = frozenset({
    'ARM1_NEVER_CONDITIONED',
    'ARM2_CONDITIONED_DEPLETED',
    'ARM3_CONTINUOUSLY_ACTIVE',
})

ALLOWED_ADJUNCT_DEPLETION_MODES: frozenset[str] = frozenset({
    'NATURAL_AGING_DECAY',
    'PHYSICAL_WASHOUT',
    'NOT_APPLICABLE',
})

ALLOWED_ADJUNCT_SA_FHL1_STATUSES: frozenset[str] = frozenset({
    'DEPLETED_VERIFIED',
    'RETAINED_ACTIVE',
    'ABSENT_FROM_LOT',
})

ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES: frozenset[str] = frozenset({
    'RETAINED_CONDITIONED_RESPONDERS',
    'FRESH_RECIPIENT_RESPONDERS',
    'PAIRED_BOTH',
})

UNSHIELDED_REFERENCE_SPLICE_CHALLENGER_ARMS: frozenset[str] = frozenset({
    'S1_L3_CAN_SYN2_4090',
    'S1_L3_CAN_SYN4_4090',
    'S1_L3_CAN_BGH_ONLY_4090',
    'S1_L3_CAN_WPRE_ONLY_4090',
    'S2_L3_SF_SYN4_ONLY_3799',
    'S2_L3_SF_SYN4_BGH_3799',
    'S2_L3_SF_SYN4_1xFII80_3879',
    'S2_L3_SF_SYN4_BGH_1xFII80_3879',
    'S4_L1_ASTAR_505AA_STOP_ONLY_4143',
    'S4_L1_ARM1B_500AA_STOP_ONLY_3974',
})

LOCUS3_MULTIMIR_MODEB_CONTROL_ARMS: frozenset[str] = frozenset({
    'ARM1_L3_WT_PARENT_4078',
    'ARM2_L3_V2_INTRON_NT_SCRAMBLE_4398',
    'ARM3_L3_V2_INTRON_SINGLE_IL1R1_4398',
    'ARM4_L3_V2_INTRON_SINGLE_FAS_4398',
})

LOCUS3_MULTIMIR_MODEB_DUPLEX_CHALLENGER_ARMS: frozenset[str] = frozenset({
    'ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398',
    'ARM6_L3_V2_INTRON_DUPLEX_IL1R1_FAS_C341S_TCA_4398',
})

LOCUS3_MULTIMIR_MODEB_ALL_ARMS: frozenset[str] = (
    LOCUS3_MULTIMIR_MODEB_CONTROL_ARMS | LOCUS3_MULTIMIR_MODEB_DUPLEX_CHALLENGER_ARMS
)

ALLOWED_KNOCKIN_SPLICE_CHALLENGER_ARMS: frozenset[str] = (
    UNSHIELDED_REFERENCE_SPLICE_CHALLENGER_ARMS
    | LOCUS3_MULTIMIR_MODEB_ALL_ARMS
    | frozenset({
        'PARENTAL_REFERENCE',
        'S1_L3_CAN_SYN2_PAM_SHIELDED_4090',
        'S1_L3_CAN_SYN4_PAM_SHIELDED_4090',
        'S1_L3_CAN_BGH_ONLY_PAM_SHIELDED_4090',
        'S1_L3_CAN_WPRE_ONLY_PAM_SHIELDED_4090',
        'S2_L3_SF_SYN4_ONLY_PAM_SHIELDED_3799',
        'S2_L3_SF_SYN4_BGH_PAM_SHIELDED_3799',
        'S2_L3_SF_SYN4_1xFII80_PAM_SHIELDED_3879',
        'S2_L3_SF_SYN4_BGH_1xFII80_PAM_SHIELDED_3879',
        'S2B_L3_BFB_SYN4_PAM_SHIELDED_4078',
        'S2B_L3_BFB_SYN4_BGH_PAM_SHIELDED_4078',
        'S3_L4_WINA_ONLY_3457',
        'S3_L4_EEF230C_ONLY_3457',
        'S3_L4_SPLICE_230C_3457',
        'S3_L4_CODIR_7NT_230C_3457',
        'S3_L4_CONVERGENT_WTCARGO_3457',
        'S3_L4_CONVERGENT_5NT_230C_3457',
        'S3_L4_CONVERGENT_2STRAND_230C_3457',
        'S3_L4_WINDOW_B_COMPARATOR',
        'S3_L4_NO_UCOE_COMPARATOR',
        'S4_L1_ASTAR_505AA_STOP_ONLY_PAM_SHIELDED_4143',
        'S4_L1_ARM1B_500AA_STOP_ONLY_PAM_SHIELDED_3974',
        'S4_L1_DXS_ARM1C_500AA_STOP_ONLY_3300',
        'S4_L1_DXS_ARM1C_500AA_STOP_ONLY_3459',
    })
)

ALLOWED_KNOCKIN_LOCUS4_ORIENTATIONS: frozenset[str] = frozenset({
    'CO_DIRECTIONAL_MINUS',
    'CONVERGENT_PLUS',
    'NOT_APPLICABLE',
})

STAGE0_CUSTOM_DEFAULT_SE: dict[str, float] = {
    'obs_tp53_ppm1d_bcor_vaf': 1e-5,
    'obs_castseq_translocation_freq': 1e-5,
    'obs_dxs_locus4_joint_retention_fraction': 1e-5,
    'obs_dxs_cd47low_icasp9low_fraction': 1e-7,
    'obs_dxs_ap1903_cd47high_escape_fraction': 1e-9,
    'obs_capsule_d_o2_cm2_s': 1e-7,
    'obs_d_o2_gel_aged_cm2_s': 1e-7,
    'obs_serpinb9_hour0_uM': 0.02,
    'obs_gate0c_hour0_serpinb9_uM': 0.02,
    'obs_pdl1_hour0_mols': 150.0,
    'obs_gate0c_hour0_pdl1_molecules_per_cell': 150.0,
    'obs_regimeII_lysis_2h_pct': 0.5,
    'obs_gate0c_2h_lysis_pct': 0.5,
    'obs_hc10_surface_fold_wt': 0.05,
    'obs_gate0e_hc10_surface_mfi_ratio': 0.05,
    'obs_chop_coincidence_fold': 0.05,
    'obs_capsule_dz_um': 0.5,
    'obs_gate0f_capsule_thickness_um': 0.5,
    'obs_basal_insulin_frac_peak': 0.01,
    'obs_glucose_threshold_mM': 0.1,
    'obs_foreign_hlac_cd8_lysis_pct': 0.5,
    'obs_cao2_kinetic_law': 0.05,
    'obs_lot_ocr_ratio': 0.015,
    'obs_residual_psc_effective_fraction': 1e-9,
    'obs_20q11_bcl2l1_cnv_ratio': 0.01,
    'obs_20q11_bcl2l1_single_cell_mosaic_fraction': 0.001,
    'obs_ap1903_pulsed_pk_log10_kill': 0.05,
    'obs_ap1903_monotherapy_pulsed_pk_log10_kill': 0.05,
    'obs_gate0v_cvb_day14_log10_viral_reduction': 0.05,
    'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct': 0.001,
    'obs_post_stress_20q11_bcl2l1_outgrowth_fold': 0.01,
    'obs_assembled_cassette_chandler_loop_sc5b9_fold': 0.05,
    'obs_assembled_cassette_chandler_loop_platelet_retention_pct': 0.5,
    'obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction': 0.01,
    'obs_apoptotic_effero_clearance_retention_pct_of_wt': 0.5,
    'obs_nutlin3a_p21_mdm2_fold_induction': 0.1,
    'obs_vascular_pericyte_ng2_pdgfrb_coverage_pct': 0.5,
    'obs_assembled_cluster_slc18a1_enterochromaffin_pct': 0.05,
    'obs_assembled_cluster_ins_gcg_polyhormonal_pct': 0.05,
    'obs_assembled_cluster_ki67_progenitor_pct': 0.02,
    'obs_post_stress_mafa_nkx61_beta_identity_retention_pct': 0.5,
    'obs_cvb4_10moi_72h_viral_titer_log10_pfu': 0.05,
    'obs_nutlin3a_p21_cdkn1a_fold_induction': 0.1,
    'obs_sc_atac_ec_chromatin_spillover_pct': 0.1,
    'obs_endocrine_synaptophysin_pos_purity_pct': 0.2,
    'obs_ductal_sox9_krt19_double_pos_pct': 0.01,
    'obs_enterochromaffin_fev_tph1_pos_pct': 0.05,
    'obs_c3b_cd55_cd59_complement_deposition_mfi_ratio': 0.05,
    'obs_nkp46_innate_nk_lysis_pct': 0.5,
    'obs_sa_fasl_microgel_nk_degranulation_cd107a_pct': 0.2,
    'obs_safasl_delta_tunel_pp': 0.5,
    'obs_safasl_zb4_anti_fas_rescue_pp': 0.5,
    'obs_safasl_treg_viable_retention_pct_of_thbd': 1.0,
    'obs_day28_functional_beta_survival_fraction': 0.01,
    'obs_pouch_iapp_oligomer_a11_dot_blot_fold': 0.05,
    'obs_gate0m_8mM_to_2mM_insulin_shutoff_ratio': 0.01,
    # Scale-aware conservative assay-calibrated default SEs for molecule-count, rate-constant, and large-span endpoints (R7-AUDIT-REF; DEC-062):
    'obs_cd47_hour0_molecules': 2500.0,
    'obs_hlae_hour0_molecules': 1000.0,
    'obs_cd47_p10_single_cell_molecules': 5000.0,
    'obs_iec_cd47_p10_molecules': 5000.0,
    'obs_hypoxic_72h_cd47_retention_molecules': 5000.0,
    'obs_hypoxic_72h_hlae_retention_molecules': 1000.0,
    'obs_hypoxic_72h_cd47_p10_single_cell_molecules': 5000.0,
    'obs_hypoxic_72h_hlae_p10_single_cell_molecules': 1000.0,
    'obs_gate0b_hlac_sct_molecules_per_cell': 1500.0,
    'obs_serpinb9_kass_M_s': 50000.0,
    'obs_gate0b_hlac_cd8_elispot_sfw': 5.0,
    'obs_b_prime_sct_21d_primed_cd8_elispot_sfu_per_1e6': 5.0,
    'obs_periimplant_fluid_film_um': 2.5,
    'obs_gate0s_ap1903_1h_insulin_burst_iu': 0.25,
    'obs_recipient_insulin_clearance_tau_min': 0.5,
    'obs_gate0c_induction_t50_h': 0.5,
    'obs_gate0e_hc10_intracellular_fold': 0.05,
    'obs_gate0e_chop_log2fc': 0.1,
    'obs_gate0m_shutoff_t90_min': 0.2,
    'obs_phagocytosis_index_pct': 0.85,
    'obs_nkg2a_lysis_pct': 1.25,
    'obs_serpinb9_si': 0.1,
    'obs_serpinb9_half_life_h': 0.5,
    'obs_plan_a_star_cd8_lysis_pct': 0.5,
    'obs_immaturity_corner_min_glucose_mM': 0.1,
    'obs_beta_x1_x4_epistasis': 0.05,
    'obs_beta_x3_x4_epistasis': 0.05,
    'obs_tyk2_p1104a_3x2_satellite_shield_induction_pct': 1.0,
    'obs_nkg2c_adaptive_nk_lysis_pct': 1.25,
    'obs_kir_nkg2a_neg_nk_lysis_pct': 1.25,
    'obs_b_prime_sct_cd8_lysis_pct': 0.5,
    'obs_hlae_nkg2c_over_nkg2a_kd_ratio': 0.25,
    'obs_combined_stress_embodiment_a_star_lysis_pct': 0.5,
    'obs_combined_stress_comp_r4_lysis_pct': 0.5,
    'obs_combined_stress_comp_r4_gsis_si': 0.15,
    'obs_combined_stress_comp_r4_viable_beta_mass_pct': 0.5,
    'obs_combined_stress_comp_r4_shutoff_t90_min': 0.2,
    'obs_combined_stress_comp_r4_max_h2o2_uM': 0.2,
    'obs_combined_stress_embodiment_a_star_viable_beta_mass_pct': 0.5,
    'obs_combined_stress_embodiment_a_star_shutoff_t90_min': 0.2,
    'obs_combined_stress_embodiment_a_star_max_h2o2_uM': 0.2,
    'obs_alpha_glucagon_low_glc_fold': 0.1,
    'obs_alpha_counterregulatory_gain_ratio': 0.05,
    'obs_ucn3_maturity_ratio': 0.05,
    'obs_gate0m_glucagon_2mM_fold_increase': 0.1,
    'obs_indirect_apc_dq2_dq8_cd4_stimulation_index': 0.1,
    'obs_indirect_apc_cd4_ifng_si': 0.1,
    'obs_gate0d_activated_cd4_fasl_apoptosis_pct': 0.5,
    'obs_gate0d_transwell_beta_annexin_v_pct': 0.5,
    'obs_exogenous_b2m_w632_surface_fold_isotype': 0.05,
    'obs_exogenous_b2m_t1d_ctl_lysis_pct': 0.5,
    'obs_complement_c3b_opsonized_phagocytosis_pct': 0.85,
    'obs_iec_vascular_support_viability_pct': 0.5,
    'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct': 1.25,
    'obs_12pulse_unprimed_hypoxic_serial_lysis_pct': 0.5,
    'obs_il12_il18_sirpa_low_nk_lysis_pct': 1.85,
    'obs_safasl_neutrophil_chemotaxis_index': 0.1,
    'obs_iec_abo_ab_antigen_mfi_ratio': 0.05,
    'obs_post_stress_nkg2d_dnam1_attributed_lysis_pct': 0.5,
    'obs_senescent_p16_sasp_positive_cell_pct': 0.1,
    'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored': 1.0,
    'obs_gate0v_day28_persistent_vp1_ie1_cells_scored': 1.0,
    'obs_assembled_cluster_ki67_progenitor_cells_scored': 1.0,
    'obs_ductal_sox9_krt19_double_pos_cells_scored': 1.0,
    'obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored': 1.0,
    'obs_senescent_p16_sasp_cells_scored': 1.0,
    'obs_arm2_post_depletion_viable_beta_recovery_pct': 1.0,
    'obs_arm2_post_depletion_perifusion_gsis_si': 0.15,
    'obs_arm2_post_depletion_nhs_c5b9_mac_pct': 0.5,
    'obs_challenger_ivt_relative_recovery_eta': 0.02,
    'obs_challenger_umi_molecules_profiled': 1.0,
    'obs_challenger_aberrant_splice_umi_count_k': 1.0,
    'obs_challenger_host_exon_splice_retention_ratio': 0.03,
    'obs_challenger_dsrna_j2_or_isg_fold_wt': 0.05,
    'obs_hlae_sct_nkg2a_kd_uM': 0.15,
    'obs_hlae_sct_nkg2c_kd_uM': 0.35,
    'obs_wt_sct_nkg2a_kd_uM': 0.15,
    'obs_nkg2a_pos_resting_nk_lysis_pct': 0.75,
    'obs_nkg2a_pos_ifng24h_nk_lysis_pct': 0.75,
    'obs_nkg2a_pos_resting_cd107a_pct': 0.85,
    'obs_nkg2a_pos_ifng24h_cd107a_pct': 0.85,
    'obs_nkg2c_kir_pos_resting_lysis_pct': 0.75,
    'obs_nkg2c_kir_pos_ifng24h_lysis_pct': 0.75,
    'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_pct': 0.75,
    'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_pct': 0.75,
    'obs_hcmv_pooled_3subset_resting_nk_lysis_pct': 0.85,
    'obs_hcmv_pooled_3subset_ifng24h_nk_lysis_pct': 0.85,
    'obs_nkg2c_cd107a_degranulation_fold_vs_b2m_ko': 0.08,
    'obs_cmv_ul40_hlae_cd8_lysis_pct': 0.65,
    'obs_hcmv_pos_donor_count': 1e-12,
    'obs_s10_nk_failing_donor_count': 1e-12,
    'obs_ul40_reactive_donor_count': 1e-12,
    'obs_ul40_ctl_failing_donor_count': 1e-12,
    'obs_hlae_cd8aa_kd_uM': 50.0,
    'obs_hlae_cd8ab_kd_uM': 50.0,
    'obs_cd8_treg_suppression_retention_pct_of_wt': 2.5,
    'obs_cd8_treg_suppression_delta_vs_wt_pp': 1.5,
    'obs_cmv_ul40_hlae_cd8_lysis_delta_vs_wt_pp': 0.5,
    'obs_nonself_dq2_dq8_indirect_apc_si': 0.10,
    'obs_nonself_host_class1_cross_pres_lysis_pct': 0.65,
    'obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp': 0.75,
    'obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6': 5.0,
    'obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp': 0.75,
    'obs_l3_multimir_spliced_fraction': 0.01,
    'obs_l3_multimir_5prime_seed_fidelity_pct': 1.0,
    'obs_l3_multimir_ago2_mir375_retention_pct': 1.5,
    'obs_l3_multimir_shield_protein_retention_pct': 1.5,
    'obs_l3_multimir_min_receptor_knockdown_pct': 1.5,
    'obs_l3_multimir_matched_gsis_output_retention_pct': 1.5,
    'obs_s6_p10_nearest_neighbor_edge_gap_um': 2.0,
    'obs_s6_cluster_diam_p90_um': 2.0,
    'obs_s6_divider_14d_hydrated_po2_retention_pct': 1.5,
    'obs_s6_hydrogel_shear_modulus_gprime_pa': 25.0,
    'obs_p05_microbead_volume_d99_um': 1.5,
    'obs_p05_catalase_barrier_film_um': 1.0,
    'obs_p05_islet_plane_h2o2_max_uM': 0.15,
    'obs_p05_day7_catalase_retention_pct': 1.5,
    'obs_p05_336h_retained_catalase_iu_333cm2': 500.0,
    'obs_design_z_empty_lumen_5site_min_po2_mmhg': 1.0,
    'obs_design_z_vessel_standoff_p90_um': 1.0,
    'obs_design_z_loaded_surface_po2_mmhg': 0.75,
    'obs_design_z_loaded_core_po2_mmhg': 0.50,
    'obs_design_z_loaded_perifusion_gsis_si': 0.15,
    'obs_p05_catalase_endotoxin_eu_per_mg': 0.01,
    'obs_p05_total_implant_catalase_endotoxin_eu': 0.03,
    'obs_p05_day7_shed_catalase_mass_pct': 1.5,
    'obs_cuff_surface_c3b_deposition_pct_pos_ctrl': 0.4,
    'obs_cuff_soluble_sc5b9_pct_pos_ctrl': 0.4,
    'obs_cuff_macrophage_fcgammar_cytokine_fold': 0.05,
    'obs_cuff_primed_serum_fasl_ec50_fold_shift': 0.10,
    'obs_cuff_primed_serum_fasl_emax_apoptosis_pct': 1.5,
    'obs_cuff_primed_serum_thbd_apc_retention_pct': 1.5,
    'obs_cuff_primed_serum_cfh_c3b_cleavage_retention_pct': 1.5,
    'obs_cuff_dpr_ring_hydrolysis_pct': 0.2,
    'obs_cuff_14d_hsa_shed_bioactive_fasl_pct': 0.3,
    'obs_cuff_h1a_sec_mals_hexamer_purity_pct': 0.8,
    'obs_cuff_unreacted_free_thiol_pct': 0.3,
    'obs_cuff_hepatocyte_72h_caspase37_fold': 0.03,
    'obs_cuff_hepatocyte_72h_alt_ldh_fold': 0.03,
    'obs_cuff_serum_donor_count': 1e-12,
    'obs_cuff_serum_failing_donor_count': 1e-12,
    'obs_rest_120h_graft_ocr_ratio': 0.02,
    'obs_rest_depot_release_half_life_h': 2.0,
    'obs_rest_post_washout_gsis_si': 0.15,
    'obs_rest_post_washout_gsis_output_retention_pct': 1.5,
    'obs_rest_post_washout_alpha_glucagon_fold': 0.10,
    'obs_cond_sc_beta_viable_yield_pct': 1.5,
    'obs_cond_iec_viable_yield_pct': 1.5,
    'obs_cond_final_wash_supernatant_nM': 0.5,
    'obs_cond_pellet_residual_over_loq_ratio': 0.05,
    'obs_cond_pellet_extraction_recovery_pct': 1.5,
    'obs_cond_endocrine_nqo1_activity_pct': 1.5,
    'obs_cond_vascular_nqo1_activity_pct': 1.5,
    'obs_cond_post_washout_gsis_si': 0.15,
    'obs_track2_bsct_kir_nkg2a_neg_nk_lysis_pct': 0.75,
    'obs_track2_bsct_donor_count': 1e-12,
    'obs_track2_bsct_failing_donor_count': 1e-12,
    'obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct': 0.75,
    'obs_track2_hlag1_donor_count': 1e-12,
    'obs_track2_hlag1_failing_donor_count': 1e-12,
    'obs_b2m_null_pvr_cd226_nk_lysis_pct': 0.75,
    'obs_design_y01_max_lineage_escape_pct': 0.10,
    'obs_design_y01_lineage_cells_scored': 1.0,
    'obs_design_y01_cbe_offtarget_pct': 0.01,
    'obs_design_y02_pre_enrichment_ki_pct': 0.25,
    'obs_design_y02_clonal_ddpcr_cn': 0.02,
    'obs_design_y02_wt_b2m_residual_allele_frac': 1e-12,
    'obs_design_y02_backbone_kanr_vaf': 1e-12,
    'obs_design_y02_wgs_pseudosite_k': 1e-12,
    'obs_design_y03_iec_beta_selectivity_ratio': 2.0,
    'obs_design_y03_beta_mir375_retention_pct': 1.5,
    'obs_design_y03_coculture_gsis_si': 0.15,
    'obs_design_y03_gsis_retention_pct': 1.5,
    'obs_design_y03_supernatant_adenosine_nM': 1.0,
    'obs_design_y03_dpcpx_gsis_delta_pct': 1.0,
    'obs_design_y04_multicistron_retention_pct': 0.5,
    'obs_design_y04_nhs20_c5b9_mac_pct': 0.4,
    'obs_design_y04_indirect_apc_cd4_si': 0.10,
    'obs_design_y04_ap1903_clonogenic_survival_frac': 1e-12,
    'obs_design_y_sampling_age_days': 1e-12,
    'obs_design_y02_ont_intact_homolog_pct': 0.2,
    'obs_design_y04_all4_lineages_scored_count': 1e-12,
    'obs_design_y04_ap1903_survivor_count_k': 1e-12,
    'obs_design_y04_ap1903_effective_trials': 1.0,
    'obs_mcb_donor_pairwise_ibd_pi_hat': 0.001,
    'obs_mcb_1q32_mdm4_cnv_ratio': 0.01,
    'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction': 0.001,
    'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored': 1.0,
    'obs_1q32_mdm4_cnv_ratio': 0.01,
    'obs_1q32_mdm4_single_cell_mosaic_fraction': 0.001,
    'obs_1q32_mdm4_single_cell_mosaic_cells_scored': 1.0,
    'obs_readout2_perifusion_gsis_si': 0.15,
    'obs_readout2_first_phase_peak_to_basal_ratio': 0.20,
    'obs_readout2_first_phase_time_to_peak_min': 0.25,
    'obs_readout2_stepdown_shutoff_t50_min': 0.10,
    'obs_exact_product_net_cd16_adcc_lysis_pct': 0.35,
    'obs_cd16_adcc_positive_control_lysis_pct': 2.0,
    'obs_p05_total_loaded_catalase_protein_mg': 0.05,
    'obs_gate0s_shedding_ddpcr_events_k': 1e-12,
    'obs_gate0s_shedding_ddpcr_trials_n': 1.0,
    'obs_gate0s_shedding_viable_events_k': 1e-12,
    'obs_gate0s_shedding_viable_trials_n': 1.0,
    'obs_stage0_independent_lot_replicates_n': 1e-12,
    'obs_iapp_3metric_thioflavin_s_area_pct': 0.04,
    'obs_iapp_3metric_soluble_oligomer_fold': 0.05,
    'obs_iapp_3metric_insoluble_fibril_fold': 0.05,
    'obs_p05_precure_active_catalase_input_u_333cm2': 500.0,
    'obs_p05_postcure_t0_active_catalase_iu_333cm2': 500.0,
    'obs_p05_cure_immobilization_yield_pct': 1.0,
    'obs_p05_336h_catalase_retention_pct': 1.5,
    'obs_p05_z70_gel_interface_h2o2_max_uM': 0.10,
    'obs_gate0s_day30_en_bloc_recovery_fraction': 0.002,
    'obs_gate0s_day90_en_bloc_recovery_fraction': 0.002,
    'obs_contrast_b1_paired_min_face_po2_gain_mmhg': 0.25,
    'obs_contrast_b1_80um_worst_face_po2_mmhg': 0.40,
    'obs_contrast_b1_80um_perifusion_gsis_si': 0.15,
    'obs_contrast_b2_human_exvivo_viable_window_h': 1.0,
    'obs_contrast_b2_omental_vs_suprafascial_po2_gain_mmhg': 0.50,
    'obs_contrast_b2_porcine_day90_en_bloc_recovery_fraction': 0.002,
    'obs_contrast_b3_iec_omitted_surface_po2_mmhg': 0.40,
    'obs_contrast_b3_iec_omitted_outer_kperf_s_inv': 0.005,
    'obs_contrast_b3_iec_omitted_perifusion_gsis_si': 0.15,
    'obs_design_b_hlac_phased_guide_cfd_max': 1e-6,
    'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio': 0.10,
    'obs_post_edit_segmented_baf_chr19q': 0.008,
    'obs_post_edit_segmented_baf_chr11p15': 0.008,
    'obs_post_edit_segmented_baf_chr6p21': 0.008,
    'obs_design_b_locus5_5prime_race_coord_bp': 1e-12,
    'obs_design_b_locus5_3prime_race_coord_bp': 1e-12,
    'obs_design_b_locus5_race_umi_n_eff': 1.0,
    'obs_design_b_locus5_aberrant_umi_k': 1e-12,
    'obs_tyk2_pe4max_pentaplex_msi_shifts_k': 1e-12,
    'obs_tyk2_pe4max_wgs_repeat_indels_clone': 1e-12,
    'obs_tyk2_pe4max_wgs_repeat_indels_ctrl': 1e-12,
    'obs_tyk2_pe4max_coding_or_splice_de_novo_indels': 1e-12,
    'obs_arm1_unconditioned_cd4_si': 0.08,
    'obs_arm1_unconditioned_annexin_v_pct': 1.0,
    'obs_parallel_arm3_control_cd4_apoptosis_pct': 1.5,
    'obs_post_depletion_input_normalized_cpep_ratio': 0.02,
    'obs_fcr_cd16_intact_serum_total_lysis_pct': 0.40,
    'obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct': 0.35,
    'obs_adcc_isotype_background_lysis_pct': 0.30,
    'obs_igg_depletion_efficiency_pct': 0.25,
    'obs_dec084_starting_potency_ratio_out_b_over_out_a': 0.02,
    'obs_dec084_ocr_ratio_b_over_a': 0.01,
    'obs_dec084_cpep_auc_ratio_b_over_a': 0.02,
    'obs_dec084_gcg_auc_ratio_b_over_a': 0.02,
    'obs_iapp_72h_immunocapture_oligomer_elisa_fold': 0.05,
    'obs_iapp_72h_colocalized_fibril_area_pct': 0.04,
    'obs_iapp_d28_immunocapture_oligomer_elisa_fold': 0.05,
    'obs_iapp_d28_colocalized_fibril_area_pct': 0.04,
    'obs_iapp_d84_immunocapture_oligomer_elisa_fold': 0.05,
    'obs_iapp_d84_colocalized_fibril_area_pct': 0.04,
    'obs_cxcl10_resting_race_umi_n_eff': 1e-12,
    'obs_cxcl10_resting_unexpected_tss_umi_k': 1e-12,
    'obs_cxcl10_resting_unresolved_origin_umi_k': 1e-12,
    'obs_cxcl10_ifng24h_race_umi_n_eff': 1e-12,
    'obs_cxcl10_ifng24h_unexpected_tss_umi_k': 1e-12,
    'obs_cxcl10_ifng24h_unresolved_origin_umi_k': 1e-12,
    'obs_bio06_retained_aav_itr_dna_pct': 0.02,
    'obs_bio06_itr_origin_rna_umi_k': 1e-12,
}

STAGE0_CUSTOM_MINIMUM_SE: dict[str, float] = {
    'obs_tp53_ppm1d_bcor_vaf': 1e-12,
    'obs_castseq_translocation_freq': 1e-12,
    'obs_dxs_locus4_joint_retention_fraction': 1e-12,
    'obs_dxs_cd47low_icasp9low_fraction': 1e-12,
    'obs_dxs_ap1903_cd47high_escape_fraction': 1e-12,
    'obs_capsule_d_o2_cm2_s': 1e-12,
    'obs_d_o2_gel_aged_cm2_s': 1e-12,
    'obs_residual_psc_effective_fraction': 1e-12,
    'obs_20q11_bcl2l1_single_cell_mosaic_fraction': 1e-6,
    'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction': 1e-6,
    'obs_1q32_mdm4_single_cell_mosaic_fraction': 1e-6,
    'obs_mcb_donor_pairwise_ibd_pi_hat': 1e-6,
    'obs_gate0s_shedding_ddpcr_events_k': 1e-12,
    'obs_gate0s_shedding_viable_events_k': 1e-12,
    'obs_stage0_independent_lot_replicates_n': 1e-12,
    'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct': 1e-12,
    'obs_cd47_hour0_molecules': 1.0,
    'obs_hlae_hour0_molecules': 1.0,
    'obs_cd47_p10_single_cell_molecules': 1.0,
    'obs_iec_cd47_p10_molecules': 1.0,
    'obs_hypoxic_72h_cd47_retention_molecules': 1.0,
    'obs_hypoxic_72h_hlae_retention_molecules': 1.0,
    'obs_hypoxic_72h_cd47_p10_single_cell_molecules': 1.0,
    'obs_hypoxic_72h_hlae_p10_single_cell_molecules': 1.0,
    'obs_gate0b_hlac_sct_molecules_per_cell': 1.0,
    'obs_pdl1_hour0_mols': 1.0,
    'obs_gate0c_hour0_pdl1_molecules_per_cell': 1.0,
    'obs_serpinb9_kass_M_s': 10.0,
    'obs_challenger_aberrant_splice_umi_count_k': 1e-12,
    'obs_hcmv_pos_donor_count': 1e-12,
    'obs_s10_nk_failing_donor_count': 1e-12,
    'obs_ul40_reactive_donor_count': 1e-12,
    'obs_ul40_ctl_failing_donor_count': 1e-12,
    'obs_cuff_serum_donor_count': 1e-12,
    'obs_cuff_serum_failing_donor_count': 1e-12,
    'obs_track2_bsct_donor_count': 1e-12,
    'obs_track2_bsct_failing_donor_count': 1e-12,
    'obs_track2_hlag1_donor_count': 1e-12,
    'obs_track2_hlag1_failing_donor_count': 1e-12,
    'obs_design_y01_cbe_offtarget_pct': 1e-6,
    'obs_design_y02_wt_b2m_residual_allele_frac': 1e-12,
    'obs_design_y02_backbone_kanr_vaf': 1e-12,
    'obs_design_y02_wgs_pseudosite_k': 1e-12,
    'obs_design_y04_ap1903_clonogenic_survival_frac': 1e-12,
    'obs_design_y_sampling_age_days': 1e-12,
    'obs_design_y04_all4_lineages_scored_count': 1e-12,
    'obs_design_y04_ap1903_survivor_count_k': 1e-12,
    'obs_design_b_hlac_phased_guide_cfd_max': 1e-12,
    'obs_design_b_locus5_5prime_race_coord_bp': 1e-12,
    'obs_design_b_locus5_3prime_race_coord_bp': 1e-12,
    'obs_design_b_locus5_race_umi_n_eff': 1e-12,
    'obs_design_b_locus5_aberrant_umi_k': 1e-12,
    'obs_tyk2_pe4max_pentaplex_msi_shifts_k': 1e-12,
    'obs_tyk2_pe4max_wgs_repeat_indels_clone': 1e-12,
    'obs_tyk2_pe4max_wgs_repeat_indels_ctrl': 1e-12,
    'obs_tyk2_pe4max_coding_or_splice_de_novo_indels': 1e-12,
    'obs_cxcl10_resting_race_umi_n_eff': 1e-12,
    'obs_cxcl10_resting_unexpected_tss_umi_k': 1e-12,
    'obs_cxcl10_resting_unresolved_origin_umi_k': 1e-12,
    'obs_cxcl10_ifng24h_race_umi_n_eff': 1e-12,
    'obs_cxcl10_ifng24h_unexpected_tss_umi_k': 1e-12,
    'obs_cxcl10_ifng24h_unresolved_origin_umi_k': 1e-12,
    'obs_bio06_retained_aav_itr_dna_pct': 1e-6,
    'obs_bio06_itr_origin_rna_umi_k': 1e-12,
}


def _scale_aware_stage0_default_se(mean_key: str, lo: float, hi: float) -> float:
    """Return an assay-calibrated or scale-proportional default SE without a unit-blind 0.05 cap on large-scale endpoints (R7-AUDIT-REF; DEC-062)."""
    if mean_key in STAGE0_CUSTOM_DEFAULT_SE:
        return float(STAGE0_CUSTOM_DEFAULT_SE[mean_key])
    span = float(hi) - float(lo)
    if span <= 5.0:
        return max(1e-5, min(0.05, span * 0.01))
    return round(span * 0.005, 6)

# Exact one-sided 95% Poisson upper bound table 0.5 * chi2_{0.95, 2*(k+1)} for k = 0..15 mutant duplex molecules (AUDIT-REF; DEC-041):
POISSON_UPPER95_LAMBDA_TABLE: tuple[float, ...] = (
    2.995732273553991,
    4.743864518390577,
    6.295793622168813,
    7.753656527939652,
    9.153519026707094,
    10.513034914451096,
    11.842395658958893,
    13.148113796999684,
    14.434649711103376,
    15.705216422284304,
    16.96221923565611,
    18.207514641439424,
    19.442569331048507,
    20.668569076007206,
    21.886485907105575,
    23.097129760760103,
)


def _regularized_gamma_q_iip0031(a: float, x: float) -> float:
    """Exact regularized upper incomplete gamma Q(a, x) = P(Poisson(x) <= a - 1) for a >= 1, x > 0."""
    if x <= 0.0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap = a
        term = 1.0 / a
        total = term
        for _ in range(50000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) <= abs(total) * 1e-15:
                break
        p = total * math.exp(-x + a * math.log(x) - gln)
        return max(0.0, min(1.0, 1.0 - p))
    fpmin = 1e-300
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, 50000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) <= 1e-15:
            break
    q = math.exp(-x + a * math.log(x) - gln) * h
    return max(0.0, min(1.0, q))


def _poisson_exact_upper95_lambda(k_int: int, alpha: float = 0.05) -> float:
    """Exact one-sided (1 - alpha) Poisson upper mean lambda solving P(Poisson(lambda) <= k_int) = alpha."""
    if k_int <= 0:
        return -math.log(alpha)
    a = float(k_int + 1)
    if k_int > 10_000_000:
        # High-order Cornish-Fisher expansion (relative error < 1e-13 for a > 1e7):
        z = 1.6448536269514722
        return a + z * math.sqrt(a) + (z * z - 1.0) / 3.0 + (z * z * z - 7.0 * z) / (36.0 * math.sqrt(a))
    lo = float(k_int)
    hi = max(a + 4.0 * math.sqrt(a) + 10.0, a * 1.5)
    while _regularized_gamma_q_iip0031(a, hi) > alpha:
        hi *= 2.0
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        if _regularized_gamma_q_iip0031(a, mid) > alpha:
            lo = mid
        else:
            hi = mid
    return hi


def _poisson_exact_upper95_rate(k_events: float, n_trials: float) -> float:
    """Return exact one-sided 95% Poisson upper rate bound lambda_upper95(k) / n_trials for all k >= 0."""
    if n_trials <= 0.0:
        return 1.0
    k_int = max(0, int(round(k_events)))
    lam_u95 = _poisson_exact_upper95_lambda(k_int, alpha=0.05)
    return lam_u95 / n_trials


# Co-supported equivalent Stage-0 endpoint pairs that must never disagree in the same payload (AUDIT-REF & AUDIT-REF; DEC-041):
# Note: obs_hc10_surface_fold_wt (fold over WT) and obs_gate0e_hc10_surface_mfi_ratio (MFI fold over isotype) are distinct physical denominators and are not forced equal.
CO_SUPPORTED_STAGE0_EQUIVALENT_KEY_PAIRS: tuple[tuple[str, str], ...] = (
    ('obs_serpinb9_hour0_uM', 'obs_gate0c_hour0_serpinb9_uM'),
    ('obs_serpinb9_hour0_se', 'obs_gate0c_hour0_serpinb9_se'),
    ('obs_pdl1_hour0_mols', 'obs_gate0c_hour0_pdl1_molecules_per_cell'),
    ('obs_pdl1_hour0_se', 'obs_gate0c_hour0_pdl1_molecules_se'),
    ('obs_regimeII_lysis_2h_pct', 'obs_gate0c_2h_lysis_pct'),
    ('obs_regimeII_lysis_2h_se', 'obs_gate0c_2h_lysis_se'),
    ('obs_capsule_dz_um', 'obs_gate0f_capsule_thickness_um'),
    ('obs_capsule_dz_se', 'obs_gate0f_capsule_thickness_se'),
    ('obs_d_o2_gel_aged_cm2_s', 'obs_capsule_d_o2_cm2_s'),
    ('obs_d_o2_gel_aged_se', 'obs_capsule_d_o2_se'),
)

# Bidirectional alias map between stage0_decision_tree_v1.json switch trigger strings and canonical contract meanKey/seKey (AUDIT-REF & AUDIT-REF):
STAGE0_JSON_TRIGGER_KEY_ALIASES: dict[str, str] = {
    'obs_regimeII_nkg2a_lysis_pct': 'obs_nkg2a_lysis_pct',
    'obs_regimeII_nkg2a_lysis_se': 'obs_nkg2a_lysis_se',
    'obs_astar_control_cd8_lysis_pct': 'obs_plan_a_star_cd8_lysis_pct',
    'obs_astar_control_cd8_lysis_se': 'obs_plan_a_star_cd8_lysis_se',
    'obs_immaturity_corner_min_g_mM': 'obs_immaturity_corner_min_glucose_mM',
    'obs_immaturity_corner_min_g_se': 'obs_immaturity_corner_min_glucose_se',
    'obs_hc10_surface_mfi_ratio': 'obs_gate0e_hc10_surface_mfi_ratio',
    'obs_hc10_surface_mfi_se': 'obs_gate0e_hc10_surface_mfi_se',
    'obs_granzyme_b_2h_lysis_pct': 'obs_regimeII_lysis_2h_pct',
    'obs_granzyme_b_2h_lysis_se': 'obs_regimeII_lysis_2h_se',
    'obs_chronic_capsule_thickness_um': 'obs_capsule_dz_um',
    'obs_chronic_capsule_thickness_se': 'obs_capsule_dz_se',
    'X1:X4': 'obs_beta_x1_x4_epistasis',
    'X1:X4_se': 'obs_beta_x1_x4_epistasis_se',
    'X3:X4': 'obs_beta_x3_x4_epistasis',
    'X3:X4_se': 'obs_beta_x3_x4_epistasis_se',
}


def _canonicalize_stage0_alias_keys(raw_dict: dict) -> dict:
    """Translate stage0_decision_tree_v1.json trigger alias keys to canonical contract meanKey/seKey, rejecting conflicting duplicates (AUDIT-REF, AUDIT-REF, & AUDIT-REF)."""
    out = dict(raw_dict)
    for alias_key, canon_key in STAGE0_JSON_TRIGGER_KEY_ALIASES.items():
        if alias_key in out:
            alias_val = out.pop(alias_key)
            if canon_key in out and out[canon_key] not in ('UNAVAILABLE', None) and alias_val not in ('UNAVAILABLE', None):
                if finite_decimal(out[canon_key], canon_key) != finite_decimal(alias_val, alias_key):
                    raise ValueError(
                        f'Conflicting Stage-0 observation values for alias {alias_key!r} ({alias_val!r}) '
                        f'and canonical key {canon_key!r} ({out[canon_key]!r}) (AUDIT-REF)'
                    )
            elif canon_key not in out or out[canon_key] in ('UNAVAILABLE', None):
                out[canon_key] = alias_val
    for k_primary, k_secondary in CO_SUPPORTED_STAGE0_EQUIVALENT_KEY_PAIRS:
        if (
            k_primary in out
            and k_secondary in out
            and out[k_primary] not in ('UNAVAILABLE', None)
            and out[k_secondary] not in ('UNAVAILABLE', None)
        ):
            if finite_decimal(out[k_primary], k_primary) != finite_decimal(out[k_secondary], k_secondary):
                raise ValueError(
                    f'Conflicting co-supported Stage-0 observation values for {k_primary!r} ({out[k_primary]!r}) '
                    f'and {k_secondary!r} ({out[k_secondary]!r}) (AUDIT-REF)'
                )
    # Cross-scale consistency checks for CHOP (fold vs 2**log2fc) and SERPINB9 k_ass (M^-1 s^-1 vs 10**log10_kass) (AUDIT-REF):
    if (
        'obs_chop_coincidence_fold' in out
        and 'obs_gate0e_chop_log2fc' in out
        and out['obs_chop_coincidence_fold'] not in ('UNAVAILABLE', None)
        and out['obs_gate0e_chop_log2fc'] not in ('UNAVAILABLE', None)
    ):
        chop_f = finite_decimal(out['obs_chop_coincidence_fold'], 'obs_chop_coincidence_fold')
        chop_l2 = finite_decimal(out['obs_gate0e_chop_log2fc'], 'obs_gate0e_chop_log2fc')
        if not math.isclose(chop_f, 2.0 ** chop_l2, rel_tol=1e-2, abs_tol=1e-2):
            raise ValueError(
                f"Conflicting CHOP fold ({chop_f}) and 2**log2FC ({2.0 ** chop_l2:.4f}) in Stage-0 observations (AUDIT-REF)"
            )
    if (
        'obs_serpinb9_kass_M_s' in out
        and 'obs_serpinb9_log10_kass' in out
        and out['obs_serpinb9_kass_M_s'] not in ('UNAVAILABLE', None)
        and out['obs_serpinb9_log10_kass'] not in ('UNAVAILABLE', None)
    ):
        kass_lin = finite_decimal(out['obs_serpinb9_kass_M_s'], 'obs_serpinb9_kass_M_s')
        kass_l10 = finite_decimal(out['obs_serpinb9_log10_kass'], 'obs_serpinb9_log10_kass')
        if not math.isclose(kass_lin, 10.0 ** kass_l10, rel_tol=1e-2, abs_tol=1.0):
            raise ValueError(
                f"Conflicting SERPINB9 k_ass ({kass_lin}) and 10**log10_kass ({10.0 ** kass_l10:.1f}) in Stage-0 observations (AUDIT-REF)"
            )
    return out

ARM1B_LOCUS1_COMPATIBLE_FAMILY_A_IDS = frozenset({
    'ISLET1_A4_ARM1B',
    'ISLET1_A4_ARM1B_TYK2_PE',
    'ISLET1_A4_X_MOL_RNLS_WT',
    'ISLET1_Z_X_MOL_140UM_ZERO_CAO2',
    'ISLET1_A4_L3_SCARFREE_3799',
    'ISLET1_A4_L3_BFB_BASALFLOOR',
    'ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500',
    'ISLET1_B_SCT_C1_CIITA',
    'ISLET1_B_SCT_C2_CIITA',
    'ISLET1_B_SCT_C1C2_CIITA',
    'ISLET1_B_SCT_C1_CIITA_TYK2_PE',
    'ISLET1_B_SCT_C2_CIITA_TYK2_PE',
    'ISLET1_B_SCT_C1C2_CIITA_TYK2_PE',
})

ARM1B_SCOPED_MEAN_KEYS = frozenset({
    'obs_stage7_yield_ratio',
    'obs_cd47_fold',
    'obs_ifng_hla_dr_dq_negative_pct',
})

GLOBAL_CANONICAL_ESTIMAND_SCALES = frozenset({
    'canonical_endpoint_contract_units_v0.8.0',
    'canonical_v0.8.0',
})

ALLOWED_ESTIMAND_SCALES_BY_MEAN_KEY = {
    'obs_stage7_yield_ratio': frozenset({
        'normalized_ratio_to_unedited_wt_control_closed_interval_[0.0, 2.0]',
        'normalized_ratio_to_unedited_wt_control',
        'ratio_to_wt_control',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_cd47_fold': frozenset({
        'fold_over_endogenous_baseline_closed_interval_[0.0, 100.0]',
        'fold_over_endogenous_baseline',
        'fold_over_wt_baseline',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_serpinb9_log10_kass': frozenset({
        'log10_second_order_rate_constant_M_inv_s_inv_closed_interval_[0.0, 10.0]',
        'log10_second_order_rate_constant_M_inv_s_inv',
        'log10_kass_M_inv_s_inv',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_icasp9_log10_kill': frozenset({
        'log10_clonal_depletion_kill_depth_closed_interval_[0.0, 10.0]',
        'log10_clonal_depletion_kill_depth',
        'log10_kill_depth',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_pouch_min_po2_mmhg': frozenset({
        'min_14d_slab_po2_mmhg_closed_interval_[0.0, 160.0]',
        'min_14d_slab_po2_mmhg',
        'mmHg',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_pouch_max_h2o2_uM': frozenset({
        'max_14d_slab_h2o2_uM_closed_interval_[0.0, 100.0]',
        'max_14d_slab_h2o2_uM',
        'uM',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_pouch_max_ph': frozenset({
        'max_14d_pouch_alkaline_ph_closed_interval_[6.0, 9.0]',
        'max_14d_pouch_alkaline_ph',
        'pH',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
    'obs_ifng_hla_dr_dq_negative_pct': frozenset({
        'ifng_48h_hla_dr_dq_negative_cell_percentage_closed_interval_[0.0, 100.0]',
        'ifng_48h_hla_dr_dq_negative_cell_percentage',
        'percentage_0_to_100',
    }) | GLOBAL_CANONICAL_ESTIMAND_SCALES,
}
ALL_ALLOWED_ESTIMAND_SCALES = frozenset().union(*ALLOWED_ESTIMAND_SCALES_BY_MEAN_KEY.values())


_PRODUCT_CONTRACT_MODULE_CACHE: dict[tuple[str, Any, Any], Any] = {}


def _load_product_contract_module(root: Path):
    root_res = root.resolve()
    py_path = root_res / 'scripts/product_contract.py'
    json_path = root_res / 'data/protocols/islet1_product_definitions_v0.8.0.json'
    try:
        st_py = py_path.stat()
        py_sig: Any = (st_py.st_mtime_ns, st_py.st_size)
    except OSError:
        py_sig = None
    try:
        st_js = json_path.stat()
        js_sig: Any = (st_js.st_mtime_ns, st_js.st_size)
    except OSError:
        js_sig = None
    cache_key = (str(root_res), py_sig, js_sig)
    cached = _PRODUCT_CONTRACT_MODULE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    spec = importlib.util.spec_from_file_location('openislet_product_contract', py_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if len(_PRODUCT_CONTRACT_MODULE_CACHE) > 32:
        _PRODUCT_CONTRACT_MODULE_CACHE.clear()
    _PRODUCT_CONTRACT_MODULE_CACHE[cache_key] = module
    return module


def validate_product_and_estimand_compatibility(
    product_variant_id: str | None,
    estimand_scale: str | None,
    observed_mean_keys: list[str] | set[str],
    *,
    root: Path | str | None = None,
) -> None:
    """Enforce registry membership & Arm-1B family compatibility for product_variant_id and canonical unit compatibility for estimand_scale (AUDIT-REF)."""
    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    obs_set = set(observed_mean_keys)
    if product_variant_id is not None:
        if not isinstance(product_variant_id, str) or not product_variant_id.strip() or product_variant_id != product_variant_id.strip():
            raise ValueError('product_variant_id must be a canonical non-whitespace-padded nonempty string (IIP-0024; DEC-103)')
        p_mod = _load_product_contract_module(root_path)
        product_def = p_mod.get_product_definition(product_variant_id, root=root_path)
        canon_id = product_def['product_variant_id']
        is_track2_b_sct = canon_id.startswith('ISLET1_B_SCT_')
        if (product_def.get('family') != 'A' and not is_track2_b_sct) or (
            bool(obs_set & ARM1B_SCOPED_MEAN_KEYS) and canon_id not in ARM1B_LOCUS1_COMPATIBLE_FAMILY_A_IDS
        ):
            raise ValueError(
                f'Incompatible product_variant_id {product_variant_id!r} (resolved to {canon_id!r}, family={product_def.get("family")!r}): '
                f'wet-lab Bayesian calibration priors are anchored to Family-A Arm-1B or Track-2 B\'-SCT Arm-1B products (AUDIT-REF; AUDIT-REF)'
            )
    if estimand_scale is not None:
        if not isinstance(estimand_scale, str) or not estimand_scale.strip():
            raise ValueError('estimand_scale must be a nonempty string')
        esc = estimand_scale.strip()
        if esc not in ALL_ALLOWED_ESTIMAND_SCALES:
            raise ValueError(
                f'Invalid estimand_scale {estimand_scale!r}: must match canonical endpoint unit specification (AUDIT-REF)'
            )
        for mk in sorted(obs_set):
            allowed_for_mk = ALLOWED_ESTIMAND_SCALES_BY_MEAN_KEY.get(mk, GLOBAL_CANONICAL_ESTIMAND_SCALES)
            if esc not in allowed_for_mk:
                raise ValueError(
                    f'Incompatible estimand_scale {estimand_scale!r} for observed endpoint {mk!r}; allowed scales: {sorted(allowed_for_mk)} (AUDIT-REF)'
                )


STAGE0_ORIGINAL_PRIMARY_MANDATORY_MEAN_KEYS: tuple[str, ...] = (
    'obs_stage7_yield_ratio',
    'obs_cd47_fold',
    'obs_serpinb9_log10_kass',
    'obs_icasp9_log10_kill',
    'obs_pouch_min_po2_mmhg',
    'obs_pouch_max_h2o2_uM',
    'obs_pouch_max_ph',
    'obs_ifng_hla_dr_dq_negative_pct',
)

STAGE0_RELEASE_CRITICAL_MOLECULE_AND_RATE_MEAN_KEYS: tuple[str, ...] = (
    'obs_cd47_hour0_molecules',
    'obs_hlae_hour0_molecules',
    'obs_cd47_p10_single_cell_molecules',
    'obs_hypoxic_72h_cd47_p10_single_cell_molecules',
    'obs_hypoxic_72h_hlae_p10_single_cell_molecules',
    'obs_hypoxic_72h_cd47_retention_molecules',
    'obs_hypoxic_72h_hlae_retention_molecules',
    'obs_iec_cd47_p10_molecules',
    'obs_gate0b_hlac_sct_molecules_per_cell',
    'obs_serpinb9_kass_M_s',
    'obs_nkg2a_lysis_pct',
    'obs_kir_nkg2a_neg_nk_lysis_pct',
    'obs_nkg2c_adaptive_nk_lysis_pct',
    'obs_il12_il18_sirpa_low_nk_lysis_pct',
    'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct',
    'obs_phagocytosis_index_pct',
    'obs_complement_c3b_opsonized_phagocytosis_pct',
)

STAGE0_DECISION_TREE_BOOLEAN_KEYS: tuple[str, ...] = (
    'require_full_5endpoint_s9_panel',
    'require_full_7metric_combined_stress_panel',
    's1_gzmb_kass_or_pdl1_scar_failed',
    's1_l3_bfb_replacement_qualified',
    's2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior',
    's3_arm1b_and_arm1x_below_15k_and_arm1c_splicing_ge_0p95',
    's4_rnls_wt_clears_gsis_shutoff_and_ocr_and_er_gate',
    's5_fasl_delta_tunel_gt_5pp',
    's5_sa_thbd_fhl1_only_replacement_qualified',
    's6_er_rescue_genotype_requalified',
    'b_prime_sct_recipient_cd8_mlc_cleared',
    'track2_b_sct_or_hlag1_lilrb1_rescue_qualified',
    'track2_b_sct_or_hlag1_product_switched',
    'cd47_playbook_200k_and_5pct_phago_waiver_approved',
    's8_fasting_leak_or_ucn3_failed',
    's8_staged_dose_and_omental_portal_drainage_qualified',
    's9_mineral_pore_occlusion_gt_15pct',
    's9_chronic_80um_lowload_or_omental_scaffold_qualified',
    's9_chronic_80um_lowload_qualified',
    's9_omental_revascularization_qualified',
    's9_staged_cassette_excision_and_d10w_clamp_qualified',
    'periimplant_closed_suction_drain_verified',
    'obs_vascularized_explant_en_bloc_recovery_verified',
    'obs_ap1903_pulsed_pk_washout_regimen_verified',
    'obs_masked_perfusion_confirmed',
    's10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct',
    'combined_stress_coupon_cleared_lysis_le_8pct_and_si_ge_2p5',
    'prior_batch_sticky_s9_safety_stop',
    'reset_prior_safety_stop_after_investigation',
    'oos_retest_investigation_approved',
    'allow_cross_clone_pooling_for_exploratory_research',
    'obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified',
    'obs_arm2_sa_thbd_bioactivity_depleted_verified',
    'obs_arm2_contemporaneous_sham_handling_matched',
    'obs_arm3_ongoing_adjunct_bioactivity_verified',
    'obs_challenger_ivt_rna_recovery_standard_verified',
    'obs_challenger_pam_seed_shield_coconversion_verified',
    'obs_challenger_host_isoforms_above_loq_verified',
    'obs_challenger_sv40pa_and_host_3prime_cleavage_checked',
    's10_is_replacement_challenger_lot',
    's10_race_5p_3p_junction_verified',
    'hla_a2_d227k_t228a_positive_control_spr_verified',
    'b2_tyk2_saptkflel_allele_assigned_ms_detected',
    'b2_tyk2_peptide_pulsed_positive_control_verified',
    'b2_tyk2_hlac_blockade_abrogation_verified',
    'l3_multimir_pre_aav6_2x2_transient_kd_falsifier_cleared',
    'l3_multimir_low_ifn_and_resting_shield_retention_verified',
    'l3_multimir_arm6_c341s_hlac_and_dq2_dq8_peptide_screen_cleared',
    's6_dose_matched_1500_ieq_verified',
    's6_divider_iso10993_leachables_cleared',
    'p05_full_336h_coupled_o2_h2o2_ca_ph_qualified',
    'design_z_empty_lumen_14d_persistence_verified',
    'design_z_loaded_closed_state_confirmed',
    'require_full_iip0020_biomaterial_and_humoral_panel',
    'p05_recombinant_human_catalase_verified',
    'cuff_dpr_ring_hydrolysis_lcms_verified',
    'cuff_conjugate_disulfide_and_site_occupancy_lcms_verified',
    'cuff_hexafasl_inventory_within_provisional_cap_verified',
    'require_full_iip0021_ocr_and_conditioning_panel',
    'rest_120h_longitudinal_simultaneous_coverage_verified',
    'rest_120h_viable_lineage_yield_and_during_exposure_function_verified',
    'conditioning_pellet_loq_and_extraction_recovery_verified',
    'conditioning_orthogonal_pathway_controls_verified',
    'track2_b_sct_kir_rescue_qualified',
    'track2_hlag1_lilrb1_rescue_qualified',
    'require_full_iip0022_track2_and_guardrail_panel',
    'lineage_necroptosis_time_resolved_viability_and_rescue_verified',
    'design_y_attlr_outside_orf_5p3p_race_verified',
    'design_y_zero_2a_scar_proteomics_verified',
    'design_y_dpcpx_adora1_antagonist_control_run',
    'design_y_anti_cd200_blockade_reversal_verified',
    'design_y_y04_full_coculture_and_clonogenic_confirmed',
    'require_full_iip0024_atomic_provenance_and_bank_panel',
    'cd16_adcc_igg_depletion_paired_control_verified',
    'iapp_age_matched_reference_above_loq_verified',
    'mcb_imprinting_ic1_ic2_igdmr_conjunctive_cleared',
    'stage7_lot_imprinting_ic1_ic2_igdmr_conjunctive_cleared',
    'mcb_and_lot_xci_46xx_3axis_cleared_or_46xy',
    'mcb_and_lot_mtdna_1000x_numt_filtered_cleared',
    'require_full_iip0025_hurdle4_transport_and_challenger_panel',
    'p05_two_timepoint_same_lot_envelope_verified',
    'contrast_b2_route_qualification_distinct_from_portfolio_credit_verified',
    'contrast_b3_zero_primate_host_ingrowth_verified',
    'require_full_iip0026_molecular_provenance_panel',
    'donor_seed_gnomad_target_site_wgs_verified',
    'post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified',
    'design_b_hla_ab_separated_pass_delivery_verified',
    'design_b_chr6p21_1p41mb_loh_sv_cleared',
    'design_b_hlac_dt9_positive',
    'design_b_hlae_3d12_positive',
    'design_b_hlac_allele_specific_surface_or_phased_mrna_verified',
    'b_sct_ciita_exon5_and_trimer_verified',
    'tyk2_phased_haplotype_sap_verified',
    'require_full_iip0028_postdepletion_adcc_mixture_and_iapp_panel',
    'washout_and_natural_decay_unpooled_verified',
    'carrier_and_soluble_adjunct_below_lod_verified',
    'fcr_effector_fcgr3a_158v_typed_verified',
    'dec084_measured_starting_composition_and_fresh_equal_output_verified',
    'iapp_sister_visit_lot_and_viable_beta_mass_concordance_verified',
    'require_full_iip0027_promoter_origin_and_itr_panel',
    'cxcl10_insert_matched_comparator_control_verified',
    'cxcl10_resting_expression_below_loq_confirmed',
    'bio06_itr_rna_origin_attribution_verified',
    'require_full_iip0029_twoclock_maturation_panel',
    'obs_alk5i_free_maturation_applied',
    'terminal_alk5i_tgfb_inhibitor_continued',
    'stage7_t3_10nM_and_nac_1mM_retained',
    'zm447439_washout_verified_pre_iec',
    'pre_coassembly_iec_vascular_function_verified',
    'readout2_both_40_and_100mmhg_tracer_aligned_verified',
    'require_full_iip0031_locus2_and_tier6b_panel',
    'ciita_ki_allele_informative_transcript_recovery_verified',
    'ciita_ko_allele_informative_transcript_recovery_verified',
    'ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified',
    'lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified',
    'l2_v3_same_background_intronless_control_linked_to_cargo',
    'l2_v3_same_background_scramble_control_linked_to_cargo',
    'phasing_verified',
    'cis_verified',
    'wt_inducibility_verified',
    'rnai_attribution_verified',
    'obs_l2_v3_chr16_ha_wgs_pass',
)


def derive_contract(root: Path | str | None = None) -> dict:
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    tree = ast.parse((root / 'scripts/ingest_wetlab_results.py').read_text())
    constants = {}
    for node in tree.body:
        name = node.targets[0].id if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name) else node.target.id if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) else None
        if name in {'ENDPOINT_SPECS', 'EXTENDED_STAGE0_STAGE1_SPECS', 'MIN_BIOLOGICAL_SE_FLOOR'}:
            constants[name] = ast.literal_eval(node.value)
    endpoints = []
    for row in constants['ENDPOINT_SPECS']:
        name, _pm, _ps, mean, se, lo, hi, unit = row
        endpoints.append(dict(endpoint_id=name, meanKey=mean, seKey=se, min=lo, max=hi, se_required=True, default_se=None, unit=unit))
    for row in constants['EXTENDED_STAGE0_STAGE1_SPECS']:
        name, _pm, _ps, mean, se, default, lo, hi, unit = row
        endpoints.append(dict(endpoint_id=name, meanKey=mean, seKey=se, min=lo, max=hi, se_required=False, default_se=default, unit=unit))
    primary_means = {e['meanKey'] for e in endpoints}
    stage0_tree_endpoints = []
    for eid, mean, se, lo, hi, unit in STAGE0_DECISION_TREE_ENDPOINT_SPECS:
        if mean in primary_means:
            continue
        default_se = _scale_aware_stage0_default_se(mean, lo, hi)
        min_se = STAGE0_CUSTOM_MINIMUM_SE.get(mean, constants['MIN_BIOLOGICAL_SE_FLOOR'])
        stage0_tree_endpoints.append(dict(endpoint_id=eid, meanKey=mean, seKey=se, min=lo, max=hi, se_required=False, default_se=default_se, minimum_se=min_se, unit=unit))
    return {'schema_version': '0.8.0', 'source': 'scripts/ingest_wetlab_results.py constant endpoint specifications',
            'minimum_se': constants['MIN_BIOLOGICAL_SE_FLOOR'], 'missing_mean': 'UNAVAILABLE',
            'unknown_keys': 'REJECT',
            'tier6b_sorted_lot': {
                'minimum_untreated_plated_cells_n': IIP0031_MIN_UNTREATED_CFU_PLATED_CELLS,
                'underpowered_untreated_cfu_reason': 'UNDERPOWERED_UNTREATED_CFU_CONTROL_DENOMINATOR',
            },
            'metadata_keys': [
                'batch_id',
                'assay_id',
                'product_variant_id',
                'estimand_scale',
                's8b_epistasis_estimand_scale',
                'obs_beta_x1_x4_epistasis_scale',
                'obs_beta_x3_x4_epistasis_scale',
                'specimen_or_clone_id',
                'specimen_or_clone_ids',
                'mcb_clone_id',
                'prior_variant_safety_status',
                'reset_prior_safety_stop_reason',
                'reset_prior_safety_stop_reason_code',
                'reset_prior_safety_stop_repeat_batch_ids',
                'reset_prior_safety_stop_independent_qp_approver_id',
                'oos_retest_assignable_cause_code',
                'adjunct_conditioning_history_arm',
                'adjunct_depletion_mode',
                'adjunct_sa_fhl1_status',
                'adjunct_immune_responder_history',
                'knockin_splice_challenger_arm',
                'knockin_locus4_orientation',
                's10_construct_id',
                's10_construct_donor_sha256',
                'current_lot_id',
                's10_qualified_lot_id',
                's10_hlae_floor_mode',
                'b2_tyk2_recipient_genotype_stratum',
                'b2_tyk2_hlac_allele',
                's6_cluster_size_arm',
                'p05_barrier_slab_geometry_mode',
                'design_z_divider_arm',
                'p05_catalase_uniprot_accession',
                'cuff_bioconjugation_arm',
                'cuff_serum_challenge_stratum',
                'rest_ocr_coupon_arm',
                'conditioning_4arm_arm',
                'conditioning_washout_analyte_scope',
                'track2_rescue_route',
                'design_y_integrase_modality',
                'design_y_silencing_modality',
                'design_y_payload_arm',
                'donor_seedline_id',
                'device_assembly_id',
                'physical_route_id',
                'assay_campaign_id',
                'reset_prior_safety_stop_target_category',
                'submitter_id',
                'contrast_b1_geometry_arm',
                'contrast_b2_perfusion_and_retrieval_arm',
                'contrast_b3_donor_iec_arm',
                'bound_c1_c5_pair_sha256',
                'tyk2_pe_modality',
                'iip0026_schema_era',
                'recipient_serum_pair_id',
                'intact_serum_sample_id',
                'igg_depleted_serum_sample_id',
                'positive_adcc_control_id',
                'isotype_adcc_control_id',
                'arm2_depleted_sample_id',
                'contemporaneous_arm3_control_sample_id',
                'sister_assembly_id',
                'cxcl10_tss_donor_id',
                'cxcl10_tss_transcription_unit_id',
                'construct_sha256',
                'post_sort_endocrine_zm447439_regimen',
                'stage7_start_utc',
                'sort_reaggregation_utc',
                'assay_harvest_utc',
                'registered_mixture_id',
                'locked_analysis_plan_id',
                'gate0s_elevated_store_scope_matched_device_lot_id',
                'background_id',
                'mcb_id',
                'l2_v3_donor_id',
                'tier6b_presort_lineage_panel',
                'tier6b_clonogenic_risk_unit',
            ],
            'boolean_policy_keys': list(STAGE0_DECISION_TREE_BOOLEAN_KEYS),
            'endpoints': endpoints,
            'stage0_decision_tree_endpoints': stage0_tree_endpoints}


def write_contract(root: Path | str | None = None) -> dict:
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    value = derive_contract(root)
    (root / CONTRACT_RELATIVE_PATH).write_text(json.dumps(value, indent=2) + '\n')
    return value


def finite_decimal(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f'{name} must be a finite decimal number')
    if isinstance(value, str) and not DECIMAL.fullmatch(value.strip()):
        raise ValueError(f'{name} must be a nonempty decimal/scientific number')
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f'{name} must be finite')
    # IEEE signed zero is one physical observation, not a distinct replay identity.
    return 0.0 if number == 0.0 else number


def normalize_measurements(payload: dict, *, require_observation: bool = True,
                           root: Path | str | None = None) -> tuple[dict, dict]:
    """Validate endpoints, retaining explicit missingness and declared SE defaults.

    Batch identity/replay is an additional stateful ingestor requirement; this
    shared measurement check is also used by the stateless draft validator.
    """
    if not isinstance(payload, dict):
        raise ValueError('Measurement payload must be an object')
    payload = _canonicalize_stage0_alias_keys(payload)
    contract = derive_contract(root)
    endpoints = contract['endpoints']
    stage0_endpoints = contract.get('stage0_decision_tree_endpoints', [])
    all_endpoints = list(endpoints) + list(stage0_endpoints)
    bool_policy_keys = set(contract.get('boolean_policy_keys', ['require_full_5endpoint_s9_panel']))
    allowed = set(contract['metadata_keys']) | bool_policy_keys | {e[k] for e in all_endpoints for k in ['meanKey', 'seKey']}
    unknown = payload.keys() - allowed
    if unknown:
        raise ValueError(f'Unknown measurement keys: {sorted(unknown)}')
    result = {}
    for bkey in bool_policy_keys:
        if bkey in payload and payload[bkey] is not None:
            if type(payload[bkey]) is not bool:
                raise ValueError(f'{bkey} must be a boolean when provided')
            result[bkey] = payload[bkey]
    mean_keys = {e['meanKey'] for e in all_endpoints}
    for key in contract['metadata_keys']:
        if key not in payload:
            continue
        if key == 'specimen_or_clone_ids':
            raw_map = payload[key]
            if not isinstance(raw_map, dict) or not raw_map:
                raise ValueError('specimen_or_clone_ids must be a non-empty dictionary mapping endpoint meanKey to specimen/clone ID')
            clean_map = {}
            for raw_mk, sv in raw_map.items():
                mk = STAGE0_JSON_TRIGGER_KEY_ALIASES.get(raw_mk, raw_mk)
                if mk not in mean_keys or payload.get(mk, 'UNAVAILABLE') == 'UNAVAILABLE':
                    raise ValueError(f'specimen_or_clone_ids key {raw_mk!r} must reference an observed endpoint meanKey')
                if not isinstance(sv, str) or not sv.strip() or len(sv) > 120:
                    raise ValueError(f'specimen_or_clone_ids[{raw_mk!r}] must be a nonempty string of at most 120 characters')
                if mk in clean_map and clean_map[mk] != sv.strip():
                    raise ValueError(f'Conflicting specimen_or_clone_ids for canonical endpoint {mk!r}')
                clean_map[mk] = sv.strip()
            result[key] = clean_map
        elif key == 'reset_prior_safety_stop_repeat_batch_ids':
            raw_list = payload[key]
            if (
                not isinstance(raw_list, list)
                or len(raw_list) < 2
                or not all(isinstance(x, str) and x.strip() and len(x) <= 120 for x in raw_list)
                or len({x.strip() for x in raw_list}) < 2
            ):
                raise ValueError('reset_prior_safety_stop_repeat_batch_ids must be a list of at least 2 distinct nonempty batch ID strings of at most 120 characters')
            result[key] = [x.strip() for x in raw_list]
        elif key in ('s8b_epistasis_estimand_scale', 'obs_beta_x1_x4_epistasis_scale', 'obs_beta_x3_x4_epistasis_scale'):
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES:
                raise ValueError(f'{key} must be one of {sorted(ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES)}')
            result[key] = payload[key].strip()
        elif key == 'prior_variant_safety_status':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES:
                raise ValueError(f'prior_variant_safety_status must be one of {sorted(ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES)} (AUDIT-REF)')
            result[key] = payload[key].strip()
        elif key == 'adjunct_conditioning_history_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_ADJUNCT_CONDITIONING_HISTORY_ARMS:
                raise ValueError(f'adjunct_conditioning_history_arm must be one of {sorted(ALLOWED_ADJUNCT_CONDITIONING_HISTORY_ARMS)} (IIP-0012; DEC-091)')
            result[key] = payload[key].strip()
        elif key == 'adjunct_depletion_mode':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_ADJUNCT_DEPLETION_MODES:
                raise ValueError(f'adjunct_depletion_mode must be one of {sorted(ALLOWED_ADJUNCT_DEPLETION_MODES)} (IIP-0012; DEC-091)')
            result[key] = payload[key].strip()
        elif key == 'adjunct_sa_fhl1_status':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_ADJUNCT_SA_FHL1_STATUSES:
                raise ValueError(f'adjunct_sa_fhl1_status must be one of {sorted(ALLOWED_ADJUNCT_SA_FHL1_STATUSES)} (IIP-0012; DEC-091)')
            result[key] = payload[key].strip()
        elif key == 'adjunct_immune_responder_history':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES:
                raise ValueError(f'adjunct_immune_responder_history must be one of {sorted(ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES)} (IIP-0012; DEC-091)')
            result[key] = payload[key].strip()
        elif key == 'knockin_splice_challenger_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_KNOCKIN_SPLICE_CHALLENGER_ARMS:
                raise ValueError(f'knockin_splice_challenger_arm must be one of {sorted(ALLOWED_KNOCKIN_SPLICE_CHALLENGER_ARMS)} (IIP-0015; DEC-092)')
            result[key] = payload[key].strip()
        elif key == 'knockin_locus4_orientation':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_KNOCKIN_LOCUS4_ORIENTATIONS:
                raise ValueError(f'knockin_locus4_orientation must be one of {sorted(ALLOWED_KNOCKIN_LOCUS4_ORIENTATIONS)} (IIP-0015; DEC-092)')
            result[key] = payload[key].strip()
        elif key == 's10_hlae_floor_mode':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_S10_HLAE_FLOOR_MODES:
                raise ValueError(f's10_hlae_floor_mode must be one of {sorted(ALLOWED_S10_HLAE_FLOOR_MODES)} (IIP-0016; DEC-094)')
            result[key] = payload[key].strip()
        elif key == 'b2_tyk2_recipient_genotype_stratum':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_B2_TYK2_RECIPIENT_GENOTYPE_STRATA:
                raise ValueError(f'b2_tyk2_recipient_genotype_stratum must be one of {sorted(ALLOWED_B2_TYK2_RECIPIENT_GENOTYPE_STRATA)} (IIP-0017; DEC-095)')
            result[key] = payload[key].strip()
        elif key == 'b2_tyk2_hlac_allele':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_B2_TYK2_HLAC_ALLELES:
                raise ValueError(f'b2_tyk2_hlac_allele must be one of {sorted(ALLOWED_B2_TYK2_HLAC_ALLELES)} (IIP-0017; DEC-095)')
            result[key] = payload[key].strip()
        elif key == 's6_cluster_size_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_S6_CLUSTER_SIZE_ARMS:
                raise ValueError(f's6_cluster_size_arm must be one of {sorted(ALLOWED_S6_CLUSTER_SIZE_ARMS)} (IIP-0019; DEC-097)')
            result[key] = payload[key].strip()
        elif key == 'p05_barrier_slab_geometry_mode':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_P05_BARRIER_SLAB_GEOMETRY_MODES:
                raise ValueError(f'p05_barrier_slab_geometry_mode must be one of {sorted(ALLOWED_P05_BARRIER_SLAB_GEOMETRY_MODES)} (IIP-0019; DEC-097)')
            result[key] = payload[key].strip()
        elif key == 'design_z_divider_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_DESIGN_Z_DIVIDER_ARMS:
                raise ValueError(f'design_z_divider_arm must be one of {sorted(ALLOWED_DESIGN_Z_DIVIDER_ARMS)} (IIP-0019; DEC-097)')
            result[key] = payload[key].strip()
        elif key == 'p05_catalase_uniprot_accession':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_P05_CATALASE_UNIPROT_ACCESSIONS:
                raise ValueError(f'p05_catalase_uniprot_accession must be one of {sorted(ALLOWED_P05_CATALASE_UNIPROT_ACCESSIONS)} (IIP-0020; DEC-099)')
            result[key] = payload[key].strip()
        elif key == 'cuff_bioconjugation_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CUFF_BIOCONJUGATION_ARMS:
                raise ValueError(f'cuff_bioconjugation_arm must be one of {sorted(ALLOWED_CUFF_BIOCONJUGATION_ARMS)} (IIP-0020; DEC-099)')
            result[key] = payload[key].strip()
        elif key == 'cuff_serum_challenge_stratum':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CUFF_SERUM_CHALLENGE_STRATA:
                raise ValueError(f'cuff_serum_challenge_stratum must be one of {sorted(ALLOWED_CUFF_SERUM_CHALLENGE_STRATA)} (IIP-0020; DEC-099)')
            result[key] = payload[key].strip()
        elif key == 'rest_ocr_coupon_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_REST_OCR_COUPON_ARMS:
                raise ValueError(f'rest_ocr_coupon_arm must be one of {sorted(ALLOWED_REST_OCR_COUPON_ARMS)} (IIP-0021; DEC-100)')
            result[key] = payload[key].strip()
        elif key == 'conditioning_4arm_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CONDITIONING_4ARM_ARMS:
                raise ValueError(f'conditioning_4arm_arm must be one of {sorted(ALLOWED_CONDITIONING_4ARM_ARMS)} (IIP-0021; DEC-100)')
            result[key] = payload[key].strip()
        elif key == 'conditioning_washout_analyte_scope':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CONDITIONING_WASHOUT_ANALYTE_SCOPES:
                raise ValueError(f'conditioning_washout_analyte_scope must be one of {sorted(ALLOWED_CONDITIONING_WASHOUT_ANALYTE_SCOPES)} (IIP-0021; DEC-100)')
            result[key] = payload[key].strip()
        elif key == 'track2_rescue_route':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_TRACK2_RESCUE_ROUTES:
                raise ValueError(f'track2_rescue_route must be one of {sorted(ALLOWED_TRACK2_RESCUE_ROUTES)} (IIP-0022; DEC-101)')
            result[key] = payload[key].strip()
        elif key == 'design_y_integrase_modality':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_DESIGN_Y_INTEGRASE_MODALITIES:
                raise ValueError(f'design_y_integrase_modality must be one of {sorted(ALLOWED_DESIGN_Y_INTEGRASE_MODALITIES)} (IIP-0023; DEC-102)')
            result[key] = payload[key].strip()
        elif key == 'design_y_silencing_modality':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_DESIGN_Y_SILENCING_MODALITIES:
                raise ValueError(f'design_y_silencing_modality must be one of {sorted(ALLOWED_DESIGN_Y_SILENCING_MODALITIES)} (IIP-0023; DEC-102)')
            result[key] = payload[key].strip()
        elif key == 'design_y_payload_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_DESIGN_Y_PAYLOAD_ARMS:
                raise ValueError(f'design_y_payload_arm must be one of {sorted(ALLOWED_DESIGN_Y_PAYLOAD_ARMS)} (IIP-0023; DEC-102)')
            result[key] = payload[key].strip()
        elif key == 'contrast_b1_geometry_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CONTRAST_B1_GEOMETRY_ARMS:
                raise ValueError(f'contrast_b1_geometry_arm must be one of {sorted(ALLOWED_CONTRAST_B1_GEOMETRY_ARMS)} (IIP-0025; DEC-104)')
            result[key] = payload[key].strip()
        elif key == 'contrast_b2_perfusion_and_retrieval_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CONTRAST_B2_ARMS:
                raise ValueError(f'contrast_b2_perfusion_and_retrieval_arm must be one of {sorted(ALLOWED_CONTRAST_B2_ARMS)} (IIP-0025; DEC-104)')
            result[key] = payload[key].strip()
        elif key == 'contrast_b3_donor_iec_arm':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_CONTRAST_B3_DONOR_IEC_ARMS:
                raise ValueError(f'contrast_b3_donor_iec_arm must be one of {sorted(ALLOWED_CONTRAST_B3_DONOR_IEC_ARMS)} (IIP-0025; DEC-104)')
            result[key] = payload[key].strip()
        elif key == 'tyk2_pe_modality':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_TYK2_PE_MODALITIES:
                raise ValueError(f'tyk2_pe_modality must be one of {sorted(ALLOWED_TYK2_PE_MODALITIES)} (IIP-0026; DEC-105)')
            result[key] = payload[key].strip()
        elif key == 'iip0026_schema_era':
            if not isinstance(payload[key], str) or payload[key].strip() not in ALLOWED_IIP0026_SCHEMA_ERAS:
                raise ValueError(f'iip0026_schema_era must be one of {sorted(ALLOWED_IIP0026_SCHEMA_ERAS)} (IIP-0026; DEC-105)')
            result[key] = payload[key].strip()
        elif key == 'bound_c1_c5_pair_sha256':
            if not isinstance(payload[key], str) or not re.fullmatch(r'[0-9a-f]{64}', payload[key].strip()):
                raise ValueError('bound_c1_c5_pair_sha256 must be a lowercase 64-character SHA-256 hex digest (IIP-0026; DEC-105)')
            result[key] = payload[key].strip()
        elif key == 'product_variant_id':
            if not isinstance(payload[key], str) or not payload[key].strip() or payload[key] != payload[key].strip() or len(payload[key]) > 120:
                raise ValueError('product_variant_id must be a canonical non-whitespace-padded nonempty string of at most 120 characters (IIP-0024; DEC-103)')
            result[key] = payload[key]
        else:
            if not isinstance(payload[key], str) or not payload[key].strip() or len(payload[key]) > 120:
                raise ValueError(f'{key} must be a nonempty string of at most 120 characters')
            result[key] = payload[key].strip()
    if 'specimen_or_clone_id' in result and 'specimen_or_clone_ids' in result:
        sid_scalar = result['specimen_or_clone_id']
        for mk, sv in result['specimen_or_clone_ids'].items():
            if sv != sid_scalar:
                raise ValueError(
                    f'Conflicting specimen_or_clone_id ({sid_scalar!r}) and specimen_or_clone_ids[{mk!r}] ({sv!r}) in measurement payload (AUDIT-REF)'
                )
    observed = []; defaulted = []
    for endpoint in endpoints:
        mk, sk = endpoint['meanKey'], endpoint['seKey']
        mean, se = payload.get(mk, 'UNAVAILABLE'), payload.get(sk, 'UNAVAILABLE')
        if mean == 'UNAVAILABLE':
            if se != 'UNAVAILABLE':
                raise ValueError(f'{sk} has no observed {mk}')
            result[mk] = 'UNAVAILABLE'
            continue
        number = finite_decimal(mean, mk)
        if not endpoint['min'] <= number <= endpoint['max']:
            raise ValueError(f'{mk} must be in [{endpoint["min"]}, {endpoint["max"]}]')
        if sk not in payload and not endpoint['se_required']:
            se = endpoint['default_se']; defaulted.append(sk)
        error = finite_decimal(se, sk)
        if error < contract['minimum_se']:
            raise ValueError(f'{sk} must be >= {contract["minimum_se"]}')
        result[mk] = number; result[sk] = error; observed.append(mk)
    for endpoint in stage0_endpoints:
        mk, sk = endpoint['meanKey'], endpoint['seKey']
        mean, se = payload.get(mk, 'UNAVAILABLE'), payload.get(sk, 'UNAVAILABLE')
        if mean == 'UNAVAILABLE':
            if se != 'UNAVAILABLE':
                raise ValueError(f'{sk} has no observed {mk}')
            continue
        number = finite_decimal(mean, mk)
        if not endpoint['min'] <= number <= endpoint['max']:
            raise ValueError(f'{mk} must be in [{endpoint["min"]}, {endpoint["max"]}]')
        if sk not in payload or se == 'UNAVAILABLE':
            se = endpoint['default_se']; defaulted.append(sk)
        error = finite_decimal(se, sk)
        min_se = endpoint.get('minimum_se', STAGE0_CUSTOM_MINIMUM_SE.get(mk, contract['minimum_se']))
        if error < min_se:
            raise ValueError(f'{sk} must be >= {min_se}')
        result[mk] = number; result[sk] = error; observed.append(mk)
    if require_observation and not observed:
        raise ValueError('At least one endpoint must be observed')
    validate_product_and_estimand_compatibility(
        result.get('product_variant_id'),
        result.get('estimand_scale'),
        observed,
        root=root,
    )
    return result, {'validated_endpoints': observed, 'defaulted_standard_error_keys': defaulted,
                    'default_se_status': 'DECLARED_ASSUMPTION_NOT_MEASURED' if defaulted else 'NOT_USED'}


def evaluate_s10_functional_and_replacement_gate_v3(
    obs: dict,
    construct_registry: dict[str, dict[str, object]] | None = None,
    *,
    defaulted_se_keys: list[str] | set[str] | None = None,
) -> dict:
    """Fail-closed Stage-0 Gate 0A (Rule S10) functional selectivity & comparator gate v3 (IIP-0016; DEC-094).

    Enforces:
      - Exact 20-construct registry lookup and hard block on calibration controls (VMAPRTLFL, QMRPVSRVL).
      - Unconditional cryptographic donor SHA-256 (s10_construct_donor_sha256) and lot-identity binding
        (current_lot_id == s10_qualified_lot_id), plus strict s10_is_replacement_challenger_lot class match
        and Stage-0 5'/3'-RACE junction verification (s10_race_5p_3p_junction_verified == True).
      - Explicit SE bounds on all continuous endpoints, positive covalent-SCT SPR K_D bounds, and paired
        resting (Regime II) + 24h IFN-gamma (Regime I) NK subset lysis/CD107a, pooled 3-subset lysis,
        NKG2C CD107a fold vs B2M-KO, and HCMV UL40 / HLA-E CD8+ CTL 4h lysis across n >= 6 donors (0 failing).
      - Mode-specific resting HLA-E floor ('GATE_0A_SCREENING_15K' -> >= 15,000 vs 'LOT_RELEASE_25K' ->
        >= 25,000 plus 72h hypoxic p10 >= 15,000).
      - Full Layer-3 D227K/T228A biophysics, parallel HLA-A*02:01 WT vs D227K/T228A SPR positive control,
        CD8+ Treg suppression retention (>= 80% of WT and delta >= -10.0 pp), and UL40 CTL delta vs WT (<= +2.0 pp).
      - Unfiltered non-self seam / alpha3 indirect DQ2/DQ8 APC CD4 SI <= 1.5 and host Class-I CTL lysis <= 5.0%.
    """
    registry = STAGE0_S10_V3_CONSTRUCT_REGISTRY if construct_registry is None else construct_registry
    defaulted_set = set(defaulted_se_keys or ())
    z95 = 1.6448536269514722

    cid = obs.get('s10_construct_id')
    if not isinstance(cid, str) or cid not in registry:
        return {
            'gate_pass': False,
            'status': 'HOLD_INVALID_OR_UNREGISTERED_S10_CONSTRUCT_ID',
            'failed_checks': ['s10_construct_id'],
        }
    meta = registry[cid]
    if meta['promotion_eligibility'] == 'CALIBRATION_CONTROL_NEVER_PROMOTE':
        return {
            'gate_pass': False,
            'status': 'HOLD_CALIBRATION_CONTROL_NEVER_ELIGIBLE_FOR_CLINICAL_PROMOTION',
            'failed_checks': ['calibration_control_hard_block'],
        }

    sha_in = obs.get('s10_construct_donor_sha256')
    lot_id = obs.get('current_lot_id')
    qual_lot_id = obs.get('s10_qualified_lot_id')
    is_repl = obs.get('s10_is_replacement_challenger_lot')
    race_ok = obs.get('s10_race_5p_3p_junction_verified')
    floor_mode = obs.get('s10_hlae_floor_mode', 'GATE_0A_SCREENING_15K')

    if (
        type(is_repl) is not bool
        or type(race_ok) is not bool
        or not isinstance(floor_mode, str)
        or floor_mode not in ALLOWED_S10_HLAE_FLOOR_MODES
    ):
        return {
            'gate_pass': False,
            'status': 'HOLD_INVALID_BOOLEAN_OR_FLOOR_MODE_METADATA',
            'failed_checks': ['s10_boolean_or_floor_mode_typecheck'],
        }

    spec_lookup: dict[str, tuple[float, float, float]] = {
        mk: (lo, hi, STAGE0_CUSTOM_MINIMUM_SE.get(mk, 0.001))
        for _eid, mk, _sk, lo, hi, _u in STAGE0_DECISION_TREE_ENDPOINT_SPECS
        if mk in STAGE0_S10_V3_MEAN_KEYS
    }
    spec_lookup['obs_hlae_hour0_molecules'] = (0.0, 1000000.0, 1.0)
    spec_lookup['obs_hypoxic_72h_hlae_p10_single_cell_molecules'] = (0.0, 1000000.0, 1.0)

    def _read_pair(mk: str, sk: str) -> tuple[float, float] | None:
        if mk in defaulted_set or sk in defaulted_set:
            return None
        v = obs.get(mk)
        s = obs.get(sk)
        if isinstance(v, bool) or isinstance(s, bool) or not isinstance(v, (int, float)) or not isinstance(s, (int, float)):
            return None
        vf, sf = float(v), float(s)
        if not math.isfinite(vf) or not math.isfinite(sf):
            return None
        lo, hi, min_se = spec_lookup[mk]
        if not (lo <= vf <= hi) or sf < min_se:
            return None
        return vf, sf

    def _read_exact_nonneg_int(mk: str, sk: str) -> int | None:
        v = obs.get(mk)
        s = obs.get(sk, 1e-12)
        if isinstance(v, bool) or isinstance(s, bool) or not isinstance(v, (int, float)) or not isinstance(s, (int, float)):
            return None
        vf, sf = float(v), float(s)
        if not math.isfinite(vf) or not math.isfinite(sf) or vf < 0.0 or vf > 1000.0 or abs(vf - round(vf)) > 1e-9 or sf < 1e-12:
            return None
        return int(round(vf))

    req_continuous_pairs = (
        ('obs_hlae_sct_nkg2a_kd_uM', 'obs_hlae_sct_nkg2a_kd_se_uM'),
        ('obs_hlae_sct_nkg2c_kd_uM', 'obs_hlae_sct_nkg2c_kd_se_uM'),
        ('obs_hlae_hour0_molecules', 'obs_hlae_hour0_molecules_se'),
        ('obs_nkg2a_pos_resting_nk_lysis_pct', 'obs_nkg2a_pos_resting_nk_lysis_se'),
        ('obs_nkg2a_pos_ifng24h_nk_lysis_pct', 'obs_nkg2a_pos_ifng24h_nk_lysis_se'),
        ('obs_nkg2a_pos_resting_cd107a_pct', 'obs_nkg2a_pos_resting_cd107a_se'),
        ('obs_nkg2a_pos_ifng24h_cd107a_pct', 'obs_nkg2a_pos_ifng24h_cd107a_se'),
        ('obs_nkg2c_kir_pos_resting_lysis_pct', 'obs_nkg2c_kir_pos_resting_lysis_se'),
        ('obs_nkg2c_kir_pos_ifng24h_lysis_pct', 'obs_nkg2c_kir_pos_ifng24h_lysis_se'),
        ('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_se'),
        ('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_pct', 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_se'),
        ('obs_hcmv_pooled_3subset_resting_nk_lysis_pct', 'obs_hcmv_pooled_3subset_resting_nk_lysis_se'),
        ('obs_hcmv_pooled_3subset_ifng24h_nk_lysis_pct', 'obs_hcmv_pooled_3subset_ifng24h_nk_lysis_se'),
        ('obs_nkg2c_cd107a_degranulation_fold_vs_b2m_ko', 'obs_nkg2c_cd107a_degranulation_fold_se'),
        ('obs_cmv_ul40_hlae_cd8_lysis_pct', 'obs_cmv_ul40_hlae_cd8_lysis_se'),
    )
    vals: dict[str, tuple[float, float]] = {}
    for mk, sk in req_continuous_pairs:
        pr = _read_pair(mk, sk)
        if pr is None:
            return {
                'gate_pass': False,
                'status': 'HOLD_MISSING_OR_OUT_OF_DOMAIN_S10_ENDPOINT_OR_SE',
                'failed_checks': [f'{mk}/{sk}'],
            }
        vals[mk] = pr

    donor_n = _read_exact_nonneg_int('obs_hcmv_pos_donor_count', 'obs_hcmv_pos_donor_count_se')
    nk_fail_n = _read_exact_nonneg_int('obs_s10_nk_failing_donor_count', 'obs_s10_nk_failing_donor_count_se')
    ul40_resp_n = _read_exact_nonneg_int('obs_ul40_reactive_donor_count', 'obs_ul40_reactive_donor_count_se')
    ul40_fail_n = _read_exact_nonneg_int('obs_ul40_ctl_failing_donor_count', 'obs_ul40_ctl_failing_donor_count_se')
    if donor_n is None or nk_fail_n is None or ul40_resp_n is None or ul40_fail_n is None:
        return {
            'gate_pass': False,
            'status': 'HOLD_MISSING_OR_INVALID_DONOR_COUNT_PANEL',
            'failed_checks': ['donor_panel_integer_counts'],
        }

    def _u95(mk: str) -> float:
        m, s = vals[mk]
        return m + z95 * s

    def _l95(mk: str) -> float:
        m, s = vals[mk]
        return m - z95 * s

    ka_u95 = _u95('obs_hlae_sct_nkg2a_kd_uM')
    kc_l95 = max(0.0, _l95('obs_hlae_sct_nkg2c_kd_uM'))
    sct_ratio_mean = vals['obs_hlae_sct_nkg2c_kd_uM'][0] / vals['obs_hlae_sct_nkg2a_kd_uM'][0]
    hlae_h0_l95 = _l95('obs_hlae_hour0_molecules')
    req_h0_floor = 25000.0 if floor_mode == 'LOT_RELEASE_25K' else 15000.0

    checks: dict[str, bool] = {
        'unconditional_construct_sha256_and_lot_identity_bound': (
            isinstance(sha_in, str)
            and sha_in == meta['donor_sha256']
            and isinstance(lot_id, str)
            and 0 < len(lot_id.strip()) <= 120
            and lot_id == qual_lot_id
        ),
        'replacement_challenger_flag_matches_construct_class': (
            is_repl is bool(meta['is_replacement_challenger'])
        ),
        'stage0_5p_3p_race_splice_junction_verified': (race_ok is True),
        'hcmv_pos_nk_donor_panel_ge_6_with_zero_failing_donors': (
            donor_n >= 6 and nk_fail_n == 0
        ),
        'ul40_reactive_ctl_donor_panel_ge_6_with_zero_failing_donors': (
            ul40_resp_n >= 6 and ul40_fail_n == 0
        ),
        'sct_nkg2a_kd_positive_and_upper95_le_10uM': (0.001 <= vals['obs_hlae_sct_nkg2a_kd_uM'][0] and ka_u95 <= 10.0),
        'sct_nkg2c_kd_positive': (vals['obs_hlae_sct_nkg2c_kd_uM'][0] >= 0.001),
        'resting_hlae_hour0_lower95_clears_mode_floor': (hlae_h0_l95 >= req_h0_floor),
        'resting_nkg2a_pos_nk_lysis_upper95_le_8pct': (_u95('obs_nkg2a_pos_resting_nk_lysis_pct') <= 8.0),
        'ifng24h_nkg2a_pos_nk_lysis_upper95_le_8pct': (_u95('obs_nkg2a_pos_ifng24h_nk_lysis_pct') <= 8.0),
        'resting_nkg2a_pos_nk_cd107a_upper95_le_10pct': (_u95('obs_nkg2a_pos_resting_cd107a_pct') <= 10.0),
        'ifng24h_nkg2a_pos_nk_cd107a_upper95_le_10pct': (_u95('obs_nkg2a_pos_ifng24h_cd107a_pct') <= 10.0),
        'resting_nkg2c_kir_pos_lysis_upper95_le_8pct': (_u95('obs_nkg2c_kir_pos_resting_lysis_pct') <= 8.0),
        'ifng24h_nkg2c_kir_pos_lysis_upper95_le_8pct': (_u95('obs_nkg2c_kir_pos_ifng24h_lysis_pct') <= 8.0),
        'resting_double_neg_nkg2c_pos_lysis_upper95_le_8pct': (_u95('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_pct') <= 8.0),
        'ifng24h_double_neg_nkg2c_pos_lysis_upper95_le_8pct': (_u95('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_pct') <= 8.0),
        'resting_pooled_3subset_hcmv_nk_lysis_upper95_le_10pct': (_u95('obs_hcmv_pooled_3subset_resting_nk_lysis_pct') <= 10.0),
        'ifng24h_pooled_3subset_hcmv_nk_lysis_upper95_le_10pct': (_u95('obs_hcmv_pooled_3subset_ifng24h_nk_lysis_pct') <= 10.0),
        'nkg2c_cd107a_degranulation_fold_upper95_le_1p25x': (_u95('obs_nkg2c_cd107a_degranulation_fold_vs_b2m_ko') <= 1.25),
        'cmv_ul40_hlae_cd8_ctl_lysis_upper95_le_5pct': (_u95('obs_cmv_ul40_hlae_cd8_lysis_pct') <= 5.0),
    }

    if floor_mode == 'LOT_RELEASE_25K':
        p10_pair = _read_pair('obs_hypoxic_72h_hlae_p10_single_cell_molecules', 'obs_hypoxic_72h_hlae_p10_single_cell_se')
        checks['release_mode_hypoxic_72h_hlae_p10_lower95_ge_15000'] = (
            p10_pair is not None and (p10_pair[0] - z95 * p10_pair[1]) >= 15000.0
        )

    if meta['alpha3_allele'] == 'MUT_D227K_T228A':
        wt_ka = _read_pair('obs_wt_sct_nkg2a_kd_uM', 'obs_wt_sct_nkg2a_kd_se_uM')
        cd8aa = _read_pair('obs_hlae_cd8aa_kd_uM', 'obs_hlae_cd8aa_kd_se_uM')
        cd8ab = _read_pair('obs_hlae_cd8ab_kd_uM', 'obs_hlae_cd8ab_kd_se_uM')
        treg_ret = _read_pair('obs_cd8_treg_suppression_retention_pct_of_wt', 'obs_cd8_treg_suppression_retention_se')
        treg_del = _read_pair('obs_cd8_treg_suppression_delta_vs_wt_pp', 'obs_cd8_treg_suppression_delta_vs_wt_se')
        ctl_del_wt = _read_pair('obs_cmv_ul40_hlae_cd8_lysis_delta_vs_wt_pp', 'obs_cmv_ul40_hlae_cd8_lysis_delta_vs_wt_se')
        hla_a2_ctrl = obs.get('hla_a2_d227k_t228a_positive_control_spr_verified')

        checks['d227k_t228a_layer3_nkg2a_kd_upper95_le_2p5uM_and_le_1p25x_wt'] = (
            wt_ka is not None and ka_u95 <= 2.5 and ka_u95 <= 1.25 * wt_ka[0]
        )
        checks['d227k_t228a_layer3_cd8aa_and_cd8ab_kd_lower95_gt_1000uM'] = (
            cd8aa is not None
            and cd8ab is not None
            and (cd8aa[0] - z95 * cd8aa[1]) > 1000.0
            and (cd8ab[0] - z95 * cd8ab[1]) > 1000.0
        )
        checks['d227k_t228a_layer3_hla_a2_positive_control_spr_verified'] = (hla_a2_ctrl is True)
        checks['d227k_t228a_layer3_cd8_treg_retention_ge_80pct_and_delta_ge_minus10pp'] = (
            treg_ret is not None
            and treg_del is not None
            and (treg_ret[0] - z95 * treg_ret[1]) >= 80.0
            and (treg_del[0] - z95 * treg_del[1]) >= -10.0
        )
        checks['d227k_t228a_layer3_ul40_ctl_delta_vs_wt_upper95_le_2pp'] = (
            ctl_del_wt is not None
            and (ctl_del_wt[0] + z95 * ctl_del_wt[1]) <= 2.0
        )

    if meta['requires_nonself_peptide_gate']:
        dq_si = _read_pair('obs_nonself_dq2_dq8_indirect_apc_si', 'obs_nonself_dq2_dq8_indirect_apc_si_se')
        host_ctl = _read_pair('obs_nonself_host_class1_cross_pres_lysis_pct', 'obs_nonself_host_class1_cross_pres_lysis_se')
        checks['nonself_seam_or_alpha3_dq2_dq8_and_class1_cleared'] = (
            dq_si is not None
            and host_ctl is not None
            and (dq_si[0] + z95 * dq_si[1]) <= 1.5
            and (host_ctl[0] + z95 * host_ctl[1]) <= 5.0
        )

    failed = [k for k, ok in checks.items() if not ok]
    gate_pass = len(failed) == 0
    if gate_pass:
        if sct_ratio_mean >= 4.0:
            mode = 'PASS_S10_STATIC_RATIO_AND_PAIRED_FUNCTIONAL_GATE_V3'
        else:
            mode = 'PASS_S10_FUNCTIONAL_SELECTIVITY_WITH_SUB4_SCT_RATIO_RECORDED_V3'
    else:
        mode = 'HOLD_S10_FUNCTIONAL_OR_BIOPHYSICAL_GATE_V3'

    return {
        'gate_pass': gate_pass,
        'status': mode,
        'construct_id': cid,
        'floor_mode': floor_mode,
        'recorded_sct_kd_nkg2a_upper95_uM': round(ka_u95, 4),
        'recorded_sct_kd_nkg2c_lower95_uM': round(kc_l95, 4),
        'recorded_sct_kc_over_ka_mean': round(sct_ratio_mean, 4),
        'checks': checks,
        'failed_checks': failed,
    }


def _betacf_iip0024(a, b, x):
    # The binomial inversion path uses extended precision and requires verified convergence.
    extended = isinstance(x, Decimal)
    max_it = 16384 if extended else 240
    eps = Decimal('1e-29') if extended else 3.0e-15
    fpmin = Decimal('1e-48') if extended else 1e-300
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1
    d = 1 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1 / d
    h = d
    for m in range(1, max_it + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1 / d
        del_v = d * c
        h *= del_v
        if abs(del_v - 1) < eps:
            return h
    if extended:
        raise ArithmeticError('incomplete beta continued fraction did not converge')
    return h


def _reg_inc_beta_iip0024(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf_iip0024(a, b, x) / a
    return 1.0 - bt * _betacf_iip0024(b, a, 1.0 - x) / b


def student_t_cdf(t_val: float, df: float) -> float:
    if df <= 0.0 or not math.isfinite(df) or not math.isfinite(t_val):
        raise ValueError('student_t_cdf requires finite positive df and finite t_val (IIP-0024; DEC-103)')
    x = df / (df + t_val * t_val)
    ib = _reg_inc_beta_iip0024(0.5 * df, 0.5, x)
    return 1.0 - 0.5 * ib if t_val >= 0.0 else 0.5 * ib


def student_t_ppf(p: float, df: float) -> float:
    if not (0.0 < p < 1.0) or not math.isfinite(df) or df <= 0.0:
        raise ValueError('student_t_ppf requires p in (0, 1) and finite df > 0 (IIP-0024; DEC-103)')
    if abs(p - 0.5) < 1e-15:
        return 0.0
    if p < 0.5:
        return -student_t_ppf(1.0 - p, df)
    lo, hi = 0.0, 200.0
    while student_t_cdf(hi, df) < p:
        hi *= 2.0
    for _ in range(96):
        mid = 0.5 * (lo + hi)
        if student_t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def student_t_ppf_one_sided_95(df: float) -> float:
    """Exact one-sided 95% Student-t critical value t_{0.95, df} (IIP-0024; DEC-103)."""
    return student_t_ppf(0.95, df)


def student_t_q95(n_independent_units: int) -> float:
    """Exact one-sided 95% Student-t critical value for n >= 2 independent units (df = n - 1; IIP-0024; DEC-103)."""
    if type(n_independent_units) is not int or n_independent_units < 2:
        raise ValueError('student_t_q95 requires strict integer n_independent_units >= 2 (IIP-0024; DEC-103)')
    return student_t_ppf(0.95, float(n_independent_units - 1))


@lru_cache(maxsize=512)
def _log_gamma_positive_integer(n: int) -> Decimal:
    # Called in the fixed-precision context below. Use exact factorials for small parameters
    # and the log-gamma Stirling expansion for larger ones. At n >= 64, the first
    # omitted term has magnitude below 2e-28.
    if n < 64:
        return Decimal(math.factorial(n - 1)).ln()
    z = Decimal(n)
    inv2 = 1 / (z * z)
    correction = Decimal(0)
    for numerator, denominator in (
        (-3617, 122400), (1, 156), (-691, 360360), (1, 1188),
        (-1, 1680), (1, 1260), (-1, 360), (1, 12),
    ):
        correction = correction * inv2 + Decimal(numerator) / denominator
    half_log_two_pi = Decimal('0.91893853320467274178032973640561763986139747363778')
    return (z - Decimal('0.5')) * z.ln() - z + half_log_two_pi + correction / z


@lru_cache(maxsize=1024)
def _cp_interior_limit(k: int, n: int, alpha: float, upper: bool) -> float:
    # Invert the requested beta tail. Evaluating its smaller tail directly avoids
    # losing a small alpha to subtraction from 1. Round the bracket outward cleanly.
    a, b = (k + 1, n - k) if upper else (k, n - k + 1)
    lo, hi = 0.0, 1.0
    with localcontext(Context(prec=96)):
        alpha_dec = Decimal.from_float(alpha)
        log_inv_beta = _log_gamma_positive_integer(a + b) - _log_gamma_positive_integer(a) - _log_gamma_positive_integer(b)
        split = Decimal(a + 1) / (a + b + 2)

        def _is_below(p_float: float) -> bool:
            x = Decimal.from_float(p_float)
            y = Decimal(1) - x
            factor = (log_inv_beta + a * x.ln() + b * y.ln()).exp()
            if x < split:
                tail = factor * _betacf_iip0024(Decimal(a), Decimal(b), x) / a if factor else Decimal(0)
                ok_range = Decimal(0) <= tail <= Decimal(1)
                if not ok_range:
                    raise ArithmeticError('regularized beta probability outside [0, 1]')
                return tail < (Decimal(1) - alpha_dec if upper else alpha_dec)
            tail = factor * _betacf_iip0024(Decimal(b), Decimal(a), y) / b if factor else Decimal(0)
            ok_range = Decimal(0) <= tail <= Decimal(1)
            if not ok_range:
                raise ArithmeticError('regularized beta probability outside [0, 1]')
            return tail > (alpha_dec if upper else Decimal(1) - alpha_dec)

        while True:
            mid = (lo + hi) * 0.5
            if mid == lo or mid == hi:
                cand = hi if upper else lo
                if 0.0 < cand < 1.0:
                    try:
                        # Ensure exact outward rounding at the final float64 endpoint without wasting an extra ULP:
                        if upper and _is_below(cand):
                            cand = min(1.0, math.nextafter(cand, 1.0))
                        elif (not upper) and not _is_below(cand):
                            cand = max(0.0, math.nextafter(cand, 0.0))
                    except ArithmeticError:
                        cand = min(1.0, math.nextafter(hi, 1.0)) if upper else max(0.0, math.nextafter(lo, 0.0))
                return cand
            try:
                below = _is_below(mid)
            except ArithmeticError:
                # Retain the already-established outer bracket on numerical failure.
                return hi if upper else lo
            if below:
                lo = mid
            else:
                hi = mid


def _clopper_pearson_limit(k: int, n: int, alpha: float, *, upper: bool) -> float:
    if type(k) is not int or type(n) is not int or n <= 0 or not 0 <= k <= n:
        raise ValueError('Clopper-Pearson requires integer 0 <= k <= n and n > 0')
    if type(alpha) not in (int, float) or not math.isfinite(float(alpha)) or not 0 < alpha < 1:
        raise ValueError('Clopper-Pearson requires a finite numeric alpha strictly between 0 and 1')
    if upper and k == n:
        return 1.0
    if not upper and k == 0:
        return 0.0
    if n.bit_length() > 53:
        return 1.0 if upper else 0.0
    if (upper and k in (0, n - 1)) or (not upper and k in (1, n)):
        from decimal import localcontext
        with localcontext() as ctx:
            ctx.prec = 96
            a_dec = Decimal.from_float(float(alpha))
            n_dec = Decimal(n)
            if upper and k == 0:
                exact = Decimal(1) - (a_dec.ln() / n_dec).exp()
            elif not upper and k == n:
                exact = (a_dec.ln() / n_dec).exp()
            elif upper and k == n - 1:
                exact = ((Decimal(1) - a_dec).ln() / n_dec).exp()
            else:
                exact = Decimal(1) - ((Decimal(1) - a_dec).ln() / n_dec).exp()
            if exact <= 0:
                return 0.0
            if exact >= 1:
                return 1.0
            res = float(exact)
            if upper:
                while Decimal.from_float(res) < exact and res < 1.0:
                    res = math.nextafter(res, 1.0)
            else:
                while Decimal.from_float(res) > exact and res > 0.0:
                    res = math.nextafter(res, 0.0)
            return res
    return _cp_interior_limit(k, n, float(alpha), upper)


def clopper_pearson_upper95_independent(k: int, n: int, alpha: float = 0.05) -> float:
    """One-sided binomial upper confidence limit, conservatively enclosed on numerical failure."""
    return _clopper_pearson_limit(k, n, alpha, upper=True)


def clopper_pearson_lower95_independent(k: int, n: int, alpha: float = 0.05) -> float:
    """One-sided binomial lower confidence limit, conservatively enclosed on numerical failure."""
    return _clopper_pearson_limit(k, n, alpha, upper=False)


S9_CORE_LOT_ENDPOINT_IDS = (
    'obs_icasp9_log10_kill',
    'obs_tp53_ppm1d_bcor_vaf',
    'obs_castseq_translocation_freq',
)
S9_FREQUENCY_ASSAY_KEYS = {
    'obs_tp53_ppm1d_bcor_vaf': ('obs_tp53_ppm1d_bcor_duplex_depth', 'obs_tp53_ppm1d_bcor_mutant_molecules_k'),
    'obs_castseq_translocation_freq': ('obs_castseq_duplex_or_amplicon_depth', 'obs_castseq_mutant_junctions_k'),
}
S9_REPLICATE_KEYS = {
    'default': ('obs_stage0_independent_lot_replicates_n', 'replicates_n', 'sample_n'),
    'obs_icasp9_log10_kill': ('obs_icasp9_log10_kill_n', 'obs_icasp9_log10_kill_sample_n', 'obs_icasp9_log10_kill_replicates_n'),
    'obs_tp53_ppm1d_bcor_vaf': ('obs_tp53_ppm1d_bcor_vaf_n', 'obs_tp53_ppm1d_bcor_vaf_sample_n', 'obs_tp53_ppm1d_bcor_vaf_replicates_n'),
    'obs_castseq_translocation_freq': ('obs_castseq_translocation_freq_n', 'obs_castseq_translocation_n', 'obs_castseq_translocation_sample_n', 'obs_castseq_translocation_replicates_n'),
}
_S9_ASSAY_RANGES = {
    key: (low, high) for _, key, _, low, high, _ in STAGE0_DECISION_TREE_ENDPOINT_SPECS
    if key in S9_CORE_LOT_ENDPOINT_IDS or any(key in pair for pair in S9_FREQUENCY_ASSAY_KEYS.values())
}


def _resolve_s9_replicate_counts_detailed_iip0031(
    containers: tuple[dict[str, Any], ...] | list[dict[str, Any]],
    *,
    flag_explicit_none: bool = False,
) -> tuple[dict[str, int | None], set[str], bool]:
    """
    Order-independent S9 replicate-count resolver across all containers and aliases
    in each S9_REPLICATE_KEYS group, unifying lot-default and endpoint-specific keys
    identically to Atomic Stage-0 Ledger (IIP-0031; DEC-110).
    Returns (resolved_by_group, conflicted_groups, any_conflict_or_invalid):
      - For `default`, inspects `S9_REPLICATE_KEYS['default']`; for each endpoint `eid`,
        inspects `(*S9_REPLICATE_KEYS['default'], *S9_REPLICATE_KEYS[eid])` so a conflict
        between `replicates_n` (e.g. 6) and `obs_icasp9_log10_kill_n` (e.g. 100) is handled
        identically across Direct, Atomic Ledger, Serialized Ingestion, and Chronological Replay.
      - Any group with multiple distinct supplied values (type-aware, unhashable-safe)
        or any invalid value (not an integral count in [6, 10000]) is flagged.
      - For S9 hard-stop evaluation, picks the most conservative valid integer count >= 2
        (`min(valid_ge2)`, capped at <= 6 when an invalid entry also appears in the group)
        so contradictory counts always evaluate symmetrically at the wider Student-t bound.
    """
    resolved: dict[str, int | None] = {}
    conflicted_groups: set[str] = set()
    group_key_map = {
        'default': S9_REPLICATE_KEYS['default'],
        **{
            eid: (*S9_REPLICATE_KEYS['default'], *S9_REPLICATE_KEYS[eid])
            for eid in S9_CORE_LOT_ENDPOINT_IDS
        },
    }
    for group_name, group_keys in group_key_map.items():
        supplied: list[Any] = []
        explicit_none_seen = False
        for container in containers:
            if not isinstance(container, dict):
                continue
            for key in group_keys:
                if key not in container:
                    continue
                val = container[key]
                if val is None:
                    if flag_explicit_none:
                        explicit_none_seen = True
                    continue
                supplied.append(val)
        if not supplied and not explicit_none_seen:
            continue
        first_val = supplied[0] if supplied else None
        group_conflict = False
        for v in supplied[1:]:
            if _is_integral_count_iip0026(v) and _is_integral_count_iip0026(first_val):
                if int(v) != int(first_val):
                    group_conflict = True
                    break
            elif type(v) is not type(first_val) or v != first_val:
                group_conflict = True
                break
        group_invalid = explicit_none_seen or any(
            not _is_integral_count_iip0026(v) or not (6 <= int(v) <= 10000)
            for v in supplied
        )
        if group_conflict or group_invalid:
            conflicted_groups.add(group_name)
        valid_ge2 = [
            (v if type(v) is int else int(round(v)))
            for v in supplied
            if (
                (type(v) is int and 2 <= v <= 10000)
                or (isinstance(v, float) and math.isfinite(v) and 2.0 <= v <= 10000.0 and v.is_integer())
            )
        ]
        if valid_ge2:
            resolved[group_name] = min(valid_ge2 + ([6] if group_invalid else []))
        else:
            resolved[group_name] = None
    return resolved, conflicted_groups, bool(conflicted_groups)


def _evaluate_s9_tp53_vaf_guard(
    vaf_val: Any,
    vaf_se: Any,
    duplex_depth: Any,
    *,
    sample_n: Any = None,
    mutant_molecules_k: Any = None,
    require_explicit_depth: bool = True,
    cap: float = 1.0e-4,
    min_depth: float = 62958.0,
    endpoint_label: str = 'obs_tp53_ppm1d_bcor_vaf',
) -> dict[str, Any]:
    """
    Single shared Rule S9 TP53/PPM1D/BCOR VAF & CAST-Seq translocation guard across ALL 4 entry paths
    (Direct product contract, Atomic Stage-0 Ledger, Serialized JSON Ingestion, and Chronological Replay;
    IIP-0031; DEC-110).
    - Never defaults omitted duplex depth to 62,958 on live paths (`require_explicit_depth=True` holds).
    - Computes the EXACT same one-sided 95% Student-t (df = sample_n - 1, defaulting to df = 5 / n = 6
      per DEC-104 when omitted), Wald normal, and exact Poisson upper95 bounds identically across all 4 paths.
    """
    if not _is_strict_finite_number_iip0024(vaf_val) or not _is_strict_finite_number_iip0024(vaf_se):
        return {
            'status': 'HOLD_S9_OUT_OF_RANGE_VAF_OR_SE',
            'hard_stop': False,
            'upper95': 1.0,
            'depth_explicit': False,
            'depth_underpowered': True,
        }
    val_f = float(vaf_val)
    se_f = float(vaf_se)
    if val_f < 0.0 or val_f > 1.0 or se_f < 0.0:
        return {
            'status': 'HOLD_S9_OUT_OF_RANGE_VAF_OR_SE',
            'hard_stop': False,
            'upper95': 1.0,
            'depth_explicit': False,
            'depth_underpowered': True,
        }
    if (
        sample_n is not None
        and type(sample_n) is not bool
        and isinstance(sample_n, (int, float))
        and math.isfinite(float(sample_n))
        and float(sample_n) >= 2.0
    ):
        t95_mult = student_t_ppf_one_sided_95(float(int(round(float(sample_n))) - 1))
    else:
        t95_mult = student_t_ppf_one_sided_95(5.0)
    t_u95 = val_f + t95_mult * se_f
    wald_u95 = val_f + 1.645 * se_f
    param_u95 = max(t_u95, wald_u95)

    depth_key, count_key = S9_FREQUENCY_ASSAY_KEYS[endpoint_label]
    bad_count = mutant_molecules_k is not None and (
        not _is_integral_count_iip0026(mutant_molecules_k)
        or not 0 <= mutant_molecules_k <= _S9_ASSAY_RANGES[count_key][1]
    )
    if bad_count:
        mutant_molecules_k = None
    depth_explicit = duplex_depth is not None
    if not depth_explicit:
        k_ref = val_f * min_depth
        if mutant_molecules_k is not None:
            k_ref = max(float(mutant_molecules_k), k_ref)
        pois_ref = _poisson_exact_upper95_rate(k_ref, min_depth)
        eff_u95 = max(param_u95, pois_ref)
        is_stop = bool(val_f >= cap or eff_u95 >= cap)
        if is_stop:
            return {
                'status': 'FAIL_HARD_STOP',
                'hard_stop': True,
                'upper95': eff_u95,
                'depth_explicit': False,
                'depth_underpowered': bool(require_explicit_depth),
            }
        return {
            'status': 'HOLD_S9_MISSING_OR_INVALID_DUPLEX_DEPTH',
            'hard_stop': False,
            'upper95': eff_u95,
            'depth_explicit': False,
            'depth_underpowered': bool(require_explicit_depth),
        }

    if (
        not _is_integral_count_iip0026(duplex_depth)
        or not _S9_ASSAY_RANGES[depth_key][0] <= duplex_depth <= _S9_ASSAY_RANGES[depth_key][1]
    ):
        is_stop = bool(val_f >= cap or param_u95 >= cap)
        return {
            'status': 'FAIL_HARD_STOP' if is_stop else 'HOLD_S9_MISSING_OR_INVALID_DUPLEX_DEPTH',
            'hard_stop': is_stop,
            'upper95': param_u95,
            'depth_explicit': True,
            'depth_underpowered': True,
        }

    depth_f = float(duplex_depth)
    k_from_rate = val_f * depth_f
    if mutant_molecules_k is not None:
        k_sup = float(mutant_molecules_k)
        if abs(k_sup - k_from_rate) > 0.51:
            bad_count = True
        k_est = max(k_sup, k_from_rate)
    else:
        k_est = k_from_rate
    pois_u95 = _poisson_exact_upper95_rate(k_est, depth_f)
    eff_u95 = max(param_u95, pois_u95)
    if val_f >= cap or eff_u95 >= cap:
        return {
            'status': 'FAIL_HARD_STOP',
            'hard_stop': True,
            'upper95': eff_u95,
            'depth_explicit': True,
            'depth_underpowered': bool(depth_f < min_depth),
        }
    if bad_count:
        return {
            'status': 'HOLD_S9_INVALID_OR_CONFLICTING_EVENT_COUNT',
            'hard_stop': False,
            'upper95': eff_u95,
            'depth_explicit': True,
            'depth_underpowered': depth_f < min_depth,
        }
    if depth_f < min_depth:
        return {
            'status': 'HOLD_S9_DUPLEX_DEPTH_UNDERPOWERED_BELOW_62958',
            'hard_stop': False,
            'upper95': eff_u95,
            'depth_explicit': True,
            'depth_underpowered': True,
        }
    return {
        'status': 'PASS',
        'hard_stop': False,
        'upper95': eff_u95,
        'depth_explicit': True,
        'depth_underpowered': False,
    }


STAGE0_STRICT_INTEGER_COUNT_AND_DENOMINATOR_KEYS: frozenset[str] = frozenset({
    'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored',
    'obs_17q_birc5_survivin_single_cell_mosaic_cells_scored',
    'obs_gate0v_day28_persistent_vp1_ie1_cells_scored',
    'obs_cd24_sox2_oct4_3d_escapee_cells_scored',
    'obs_off_target_gp2_amylase_acinar_cells_scored',
    'obs_off_target_krt19_cftr_ductal_cells_scored',
    'obs_assembled_cluster_ki67_progenitor_cells_scored',
    'obs_ductal_sox9_krt19_double_pos_cells_scored',
    'obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored',
    'obs_senescent_p16_sasp_cells_scored',
    'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored',
    'obs_1q32_mdm4_single_cell_mosaic_cells_scored',
    'obs_mcb_effective_clonogenic_trials',
    'obs_mcb_clonogenic_survivor_count_k',
    'obs_residual_psc_effective_trials',
    'obs_hazard_state_prolif_effective_trials',
    'obs_hazard_state_prolif_survivor_count_k',
    'obs_ap1903_step1_sublethal_survivor_count_k',
    'obs_dxs_screened_cell_count',
    'obs_gate0s_shedding_ddpcr_events_k',
    'obs_gate0s_shedding_ddpcr_trials_n',
    'obs_gate0s_shedding_viable_events_k',
    'obs_gate0s_shedding_viable_trials_n',
    'obs_stage0_independent_lot_replicates_n',
    'obs_cuff_serum_donor_count',
    'obs_cuff_serum_failing_donor_count',
    'obs_hcmv_pos_donor_count',
    'obs_s10_nk_failing_donor_count',
    'obs_ul40_reactive_donor_count',
    'obs_ul40_ctl_failing_donor_count',
    'obs_track2_bsct_donor_count',
    'obs_track2_bsct_failing_donor_count',
    'obs_track2_hlag1_donor_count',
    'obs_track2_hlag1_failing_donor_count',
    'obs_challenger_aberrant_splice_umi_count_k',
    'obs_design_y02_wgs_pseudosite_k',
    'obs_design_y04_all4_lineages_scored_count',
    'obs_design_y04_ap1903_survivor_count_k',
    'obs_design_b_locus5_5prime_race_coord_bp',
    'obs_design_b_locus5_3prime_race_coord_bp',
    'obs_design_b_locus5_race_umi_n_eff',
    'obs_design_b_locus5_aberrant_umi_k',
    'obs_tyk2_pe4max_pentaplex_msi_shifts_k',
    'obs_tyk2_pe4max_wgs_repeat_indels_clone',
    'obs_tyk2_pe4max_wgs_repeat_indels_ctrl',
    'obs_tyk2_pe4max_coding_or_splice_de_novo_indels',
    'obs_cxcl10_resting_race_umi_n_eff',
    'obs_cxcl10_resting_unexpected_tss_umi_k',
    'obs_cxcl10_resting_unresolved_origin_umi_k',
    'obs_cxcl10_ifng24h_race_umi_n_eff',
    'obs_cxcl10_ifng24h_unexpected_tss_umi_k',
    'obs_cxcl10_ifng24h_unresolved_origin_umi_k',
    'obs_bio06_itr_origin_rna_umi_k',
})


def evaluate_stage0_decision_tree(
    observations: dict,
    *,
    root: Path | str | None = None,
    require_full_5endpoint_s9_panel: bool = False,
    defaulted_se_keys: list[str] | set[str] | None = None,
    prior_batch_sticky_s9_safety_stop: bool = False,
    product_def: dict | None = None,
) -> dict:
    """Deterministically evaluate the machine-readable Stage-0 promotion/fallback decision tree (S1..S9; AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, & R7S3-F01/F02; DEC-058; IIP-0024; DEC-103)."""
    if not isinstance(observations, dict):
        raise ValueError('Stage-0 decision-tree observations must be a dictionary')
    if type(require_full_5endpoint_s9_panel) is not bool:
        raise ValueError('require_full_5endpoint_s9_panel must be a boolean')
    if type(prior_batch_sticky_s9_safety_stop) is not bool:
        raise ValueError('prior_batch_sticky_s9_safety_stop must be a boolean')
    observations = _canonicalize_stage0_alias_keys(observations)
    for bkey in STAGE0_DECISION_TREE_BOOLEAN_KEYS:
        if bkey in observations and observations[bkey] is not None:
            if type(observations[bkey]) is not bool:
                raise ValueError(f'{bkey} must be a boolean when provided')
    if observations.get('require_full_5endpoint_s9_panel') is True:
        require_full_5endpoint_s9_panel = True
    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    tree_doc = json.loads((root_path / STAGE0_DECISION_TREE_RELATIVE_PATH).read_text(encoding='utf-8'))
    p_mod = _load_product_contract_module(root_path)
    # Canonical product resolution (IIP-0024 & IIP-0028; DEC-103 & DEC-106): reject whitespace-padded IDs and caller-overridden lifecycle/execution_supported fields!
    obs_pvid_raw = observations.get('product_variant_id')
    if obs_pvid_raw is not None:
        if not isinstance(obs_pvid_raw, str) or not obs_pvid_raw.strip() or obs_pvid_raw != obs_pvid_raw.strip():
            raise ValueError(f'Invalid or whitespace-padded product_variant_id: {obs_pvid_raw!r} (IIP-0024; DEC-103)')
        canon_from_obs = p_mod.get_product_definition(obs_pvid_raw, root=root_path)
        if product_def is not None:
            if not isinstance(product_def, dict):
                raise ValueError('product_def must be a dictionary when provided (IIP-0024; DEC-103)')
            kw_pvid = product_def.get('product_variant_id')
            if not isinstance(kw_pvid, str) or kw_pvid != kw_pvid.strip() or p_mod.get_product_definition(kw_pvid, root=root_path)['product_variant_id'] != canon_from_obs['product_variant_id']:
                raise ValueError('Conflicting or whitespace-padded product_def and observations product_variant_id (IIP-0024 / IIP-0028; DEC-103 / DEC-106)')
            if 'execution_supported' in product_def and product_def['execution_supported'] is not canon_from_obs['execution_supported']:
                raise ValueError('Caller-overridden execution_supported in product_def rejected (IIP-0028; DEC-106)')
            if 'lifecycle_state' in product_def and product_def.get('lifecycle_state') != canon_from_obs.get('lifecycle_state'):
                raise ValueError('Caller-overridden lifecycle_state in product_def rejected (IIP-0028; DEC-106)')
        product_def = canon_from_obs
    elif product_def is not None:
        if (
            not isinstance(product_def, dict)
            or not isinstance(product_def.get('product_variant_id'), str)
            or not product_def['product_variant_id'].strip()
            or product_def['product_variant_id'] != product_def['product_variant_id'].strip()
        ):
            raise ValueError('product_def must contain a canonical unpadded string product_variant_id (IIP-0024 / IIP-0028; DEC-103 / DEC-106)')
        canon_kw = p_mod.get_product_definition(product_def['product_variant_id'], root=root_path)
        if 'execution_supported' in product_def and product_def['execution_supported'] is not canon_kw['execution_supported']:
            raise ValueError('Caller-overridden execution_supported in product_def rejected (IIP-0028; DEC-106)')
        if 'lifecycle_state' in product_def and product_def.get('lifecycle_state') != canon_kw.get('lifecycle_state'):
            raise ValueError('Caller-overridden lifecycle_state in product_def rejected (IIP-0028; DEC-106)')
        product_def = canon_kw

    contract = derive_contract(root_path)
    default_se_by_mean: dict[str, float] = {}
    for e in list(contract['endpoints']) + list(contract.get('stage0_decision_tree_endpoints', [])):
        if e.get('default_se') is not None:
            default_se_by_mean[e['meanKey']] = float(e['default_se'])
    primary_spec_by_mean = {
        e['meanKey']: (e['endpoint_id'], e['meanKey'], e['seKey'], e['min'], e['max'], e['unit'])
        for e in contract['endpoints']
    }
    spec_by_mean = dict(primary_spec_by_mean)
    for row in STAGE0_DECISION_TREE_ENDPOINT_SPECS:
        spec_by_mean[row[1]] = row
    spec_by_se = {row[2]: row for row in spec_by_mean.values()}
    allowed_boolean_and_meta_flags = {
        'batch_id',
        'assay_id',
        'product_variant_id',
        'estimand_scale',
        's8b_epistasis_estimand_scale',
        'obs_beta_x1_x4_epistasis_scale',
        'obs_beta_x3_x4_epistasis_scale',
        'specimen_or_clone_id',
        'specimen_or_clone_ids',
        'require_full_5endpoint_s9_panel',
        'require_full_7metric_combined_stress_panel',
        's1_gzmb_kass_or_pdl1_scar_failed',
        's1_l3_bfb_replacement_qualified',
        's2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior',
        's3_arm1b_and_arm1x_below_15k_and_arm1c_splicing_ge_0p95',
        's4_rnls_wt_clears_gsis_shutoff_and_ocr_and_er_gate',
        's5_fasl_delta_tunel_gt_5pp',
        's5_sa_thbd_fhl1_only_replacement_qualified',
        's6_er_rescue_genotype_requalified',
        'b_prime_sct_recipient_cd8_mlc_cleared',
        's8_fasting_leak_or_ucn3_failed',
        's8_staged_dose_and_omental_portal_drainage_qualified',
        's9_mineral_pore_occlusion_gt_15pct',
        's9_chronic_80um_lowload_or_omental_scaffold_qualified',
        's9_chronic_80um_lowload_qualified',
        's9_omental_revascularization_qualified',
        'obs_masked_perfusion_confirmed',
        's10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct',
        'combined_stress_coupon_cleared_lysis_le_8pct_and_si_ge_2p5',
        'prior_batch_sticky_s9_safety_stop',
        'prior_variant_safety_status',
        'reset_prior_safety_stop_after_investigation',
        'reset_prior_safety_stop_reason',
        'reset_prior_safety_stop_reason_code',
        'reset_prior_safety_stop_repeat_batch_ids',
        'reset_prior_safety_stop_independent_qp_approver_id',
        'verified_ledger_passing_s9_batch_ids',
        'track2_b_sct_or_hlag1_lilrb1_rescue_qualified',
        'track2_b_sct_or_hlag1_product_switched',
        'cd47_playbook_200k_and_5pct_phago_waiver_approved',
        's9_staged_cassette_excision_and_d10w_clamp_qualified',
        'periimplant_closed_suction_drain_verified',
        'obs_vascularized_explant_en_bloc_recovery_verified',
        'obs_ap1903_pulsed_pk_washout_regimen_verified',
        'mcb_clone_id',
        'oos_retest_investigation_approved',
        'oos_retest_assignable_cause_code',
        'allow_cross_clone_pooling_for_exploratory_research',
        'obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified',
        'obs_arm2_sa_thbd_bioactivity_depleted_verified',
        'obs_arm2_contemporaneous_sham_handling_matched',
        'obs_arm3_ongoing_adjunct_bioactivity_verified',
        'adjunct_conditioning_history_arm',
        'adjunct_depletion_mode',
        'adjunct_sa_fhl1_status',
        'adjunct_immune_responder_history',
        'obs_challenger_ivt_rna_recovery_standard_verified',
        'obs_challenger_pam_seed_shield_coconversion_verified',
        'obs_challenger_host_isoforms_above_loq_verified',
        'obs_challenger_sv40pa_and_host_3prime_cleavage_checked',
        'knockin_splice_challenger_arm',
        'knockin_locus4_orientation',
        's10_construct_id',
        's10_construct_donor_sha256',
        'current_lot_id',
        's10_qualified_lot_id',
        's10_hlae_floor_mode',
        's10_is_replacement_challenger_lot',
        's10_race_5p_3p_junction_verified',
        'hla_a2_d227k_t228a_positive_control_spr_verified',
        'b2_tyk2_recipient_genotype_stratum',
        'b2_tyk2_hlac_allele',
        'b2_tyk2_saptkflel_allele_assigned_ms_detected',
        'b2_tyk2_peptide_pulsed_positive_control_verified',
        'b2_tyk2_hlac_blockade_abrogation_verified',
        'l3_multimir_pre_aav6_2x2_transient_kd_falsifier_cleared',
        'l3_multimir_low_ifn_and_resting_shield_retention_verified',
        'l3_multimir_arm6_c341s_hlac_and_dq2_dq8_peptide_screen_cleared',
        's6_cluster_size_arm',
        'p05_barrier_slab_geometry_mode',
        'design_z_divider_arm',
        's6_dose_matched_1500_ieq_verified',
        's6_divider_iso10993_leachables_cleared',
        'p05_full_336h_coupled_o2_h2o2_ca_ph_qualified',
        'design_z_empty_lumen_14d_persistence_verified',
        'design_z_loaded_closed_state_confirmed',
        'p05_catalase_uniprot_accession',
        'cuff_bioconjugation_arm',
        'cuff_serum_challenge_stratum',
        'require_full_iip0020_biomaterial_and_humoral_panel',
        'p05_recombinant_human_catalase_verified',
        'cuff_dpr_ring_hydrolysis_lcms_verified',
        'cuff_conjugate_disulfide_and_site_occupancy_lcms_verified',
        'cuff_hexafasl_inventory_within_provisional_cap_verified',
        'rest_ocr_coupon_arm',
        'conditioning_4arm_arm',
        'conditioning_washout_analyte_scope',
        'require_full_iip0021_ocr_and_conditioning_panel',
        'rest_120h_longitudinal_simultaneous_coverage_verified',
        'rest_120h_viable_lineage_yield_and_during_exposure_function_verified',
        'conditioning_pellet_loq_and_extraction_recovery_verified',
        'conditioning_orthogonal_pathway_controls_verified',
        'track2_rescue_route',
        'track2_b_sct_kir_rescue_qualified',
        'track2_hlag1_lilrb1_rescue_qualified',
        'require_full_iip0022_track2_and_guardrail_panel',
        'lineage_necroptosis_time_resolved_viability_and_rescue_verified',
        'design_y_integrase_modality',
        'design_y_silencing_modality',
        'design_y_payload_arm',
        'design_y_attlr_outside_orf_5p3p_race_verified',
        'design_y_zero_2a_scar_proteomics_verified',
        'design_y_dpcpx_adora1_antagonist_control_run',
        'design_y_anti_cd200_blockade_reversal_verified',
        'design_y_y04_full_coculture_and_clonogenic_confirmed',
        'donor_seedline_id',
        'device_assembly_id',
        'physical_route_id',
        'assay_campaign_id',
        'reset_prior_safety_stop_target_category',
        'submitter_id',
        'require_full_iip0024_atomic_provenance_and_bank_panel',
        'cd16_adcc_igg_depletion_paired_control_verified',
        'iapp_age_matched_reference_above_loq_verified',
        'mcb_imprinting_ic1_ic2_igdmr_conjunctive_cleared',
        'stage7_lot_imprinting_ic1_ic2_igdmr_conjunctive_cleared',
        'mcb_and_lot_xci_46xx_3axis_cleared_or_46xy',
        'mcb_and_lot_mtdna_1000x_numt_filtered_cleared',
        'require_full_iip0025_hurdle4_transport_and_challenger_panel',
        'p05_two_timepoint_same_lot_envelope_verified',
        'contrast_b2_route_qualification_distinct_from_portfolio_credit_verified',
        'contrast_b3_zero_primate_host_ingrowth_verified',
        'contrast_b1_geometry_arm',
        'contrast_b2_perfusion_and_retrieval_arm',
        'contrast_b3_donor_iec_arm',
        'require_full_iip0026_molecular_provenance_panel',
        'donor_seed_gnomad_target_site_wgs_verified',
        'post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified',
        'design_b_hla_ab_separated_pass_delivery_verified',
        'design_b_chr6p21_1p41mb_loh_sv_cleared',
        'design_b_hlac_dt9_positive',
        'design_b_hlae_3d12_positive',
        'design_b_hlac_allele_specific_surface_or_phased_mrna_verified',
        'b_sct_ciita_exon5_and_trimer_verified',
        'tyk2_phased_haplotype_sap_verified',
        'bound_c1_c5_pair_sha256',
        'tyk2_pe_modality',
        'iip0026_schema_era',
        'require_full_iip0028_postdepletion_adcc_mixture_and_iapp_panel',
        'washout_and_natural_decay_unpooled_verified',
        'carrier_and_soluble_adjunct_below_lod_verified',
        'fcr_effector_fcgr3a_158v_typed_verified',
        'dec084_measured_starting_composition_and_fresh_equal_output_verified',
        'iapp_sister_visit_lot_and_viable_beta_mass_concordance_verified',
        'recipient_serum_pair_id',
        'intact_serum_sample_id',
        'igg_depleted_serum_sample_id',
        'positive_adcc_control_id',
        'isotype_adcc_control_id',
        'arm2_depleted_sample_id',
        'contemporaneous_arm3_control_sample_id',
        'sister_assembly_id',
        'require_full_iip0027_promoter_origin_and_itr_panel',
        'cxcl10_insert_matched_comparator_control_verified',
        'cxcl10_resting_expression_below_loq_confirmed',
        'bio06_itr_rna_origin_attribution_verified',
        'cxcl10_tss_donor_id',
        'cxcl10_tss_transcription_unit_id',
        'require_full_iip0029_twoclock_maturation_panel',
        'obs_alk5i_free_maturation_applied',
        'terminal_alk5i_tgfb_inhibitor_continued',
        'stage7_t3_10nM_and_nac_1mM_retained',
        'zm447439_washout_verified_pre_iec',
        'pre_coassembly_iec_vascular_function_verified',
        'readout2_both_40_and_100mmhg_tracer_aligned_verified',
        'construct_sha256',
        'post_sort_endocrine_zm447439_regimen',
        'stage7_start_utc',
        'sort_reaggregation_utc',
        'assay_harvest_utc',
        'registered_mixture_id',
        'locked_analysis_plan_id',
        'gate0s_elevated_store_scope_matched_device_lot_id',
        'iip0029_twoclock_record',
        'require_full_iip0031_locus2_and_tier6b_panel',
        'ciita_ki_allele_informative_transcript_recovery_verified',
        'ciita_ko_allele_informative_transcript_recovery_verified',
        'ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified',
        'lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified',
        'l2_v3_same_background_intronless_control_linked_to_cargo',
        'l2_v3_same_background_scramble_control_linked_to_cargo',
        'phasing_verified',
        'cis_verified',
        'wt_inducibility_verified',
        'rnai_attribution_verified',
        'obs_l2_v3_chr16_ha_wgs_pass',
        'obs_l2_v3_chr16_ha_wgs_variant_count',
        'background_id',
        'mcb_id',
        'l2_v3_donor_id',
        'donor_id',
        'sample_id',
        'sample_unit_id',
        'replicates_n',
        'sample_n',
        'ciita_locus2_genotype_copy_class',
        'tier6b_presort_lineage_panel',
        'tier6b_clonogenic_risk_unit',
        'tier6b_sorted_lot_observations',
        'iip0031_locus2_record',
        's11_replicate_n',
        'rnai_specific_benefit_claimed',
        'rnai_attribution_scramble_control_id',
        'rnai_attribution_scramble_background_id',
        'rnai_attribution_intronless_control_id',
        'rnai_attribution_intronless_background_id',
        'obs_l2_v3_plvap_iec_network_t50_h',
        'obs_l2_v3_plvap_iec_network_t50_se',
        'obs_l2_v3_plvap_iec_viability_pct',
        'obs_l2_v3_plvap_iec_viability_se',
        'obs_ciita_ki_allele_informative_reads_n',
        'obs_ciita_ko_allele_informative_reads_n',
        'obs_ciita_ki_ex4_ex6_aberrant_skip_pct',
        'obs_ciita_ki_ex4_ex6_aberrant_skip_se',
        'obs_ciita_ko_ex4_ex6_aberrant_skip_pct',
        'obs_ciita_ko_ex4_ex6_aberrant_skip_se',
        'obs_ciita_ki_allele1_informative_reads_n',
        'obs_ciita_ki_allele2_informative_reads_n',
        'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_pct',
        'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_se',
        'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_pct',
        'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_se',
        'obs_edited_72h_ifng_hla_dr_dp_dq_pos_pct',
        'obs_edited_72h_ifng_hla_dr_dp_dq_pos_se',
        'obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_pct',
        'obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_se',
        'obs_cis_socs1_tost90_ratio_lo',
        'obs_cis_socs1_tost90_ratio_hi',
        'obs_cis_dexi_tost90_ratio_lo',
        'obs_cis_dexi_tost90_ratio_hi',
        'obs_cis_clec16a_tost90_ratio_lo',
        'obs_cis_clec16a_tost90_ratio_hi',
        'obs_16p13_ddpcr_cn_ratio',
        'obs_16p13_phased_baf',
        'obs_l2_v3_unit1_eef1a1_5race_umi_n_eff',
        'obs_l2_v3_unit1_eef1a1_unexpected_tss_umi_k',
        'obs_l2_v3_unit1_eef1a1_resting_0h_5race_umi_n_eff',
        'obs_l2_v3_unit1_eef1a1_resting_0h_unexpected_tss_umi_k',
        'obs_l2_v3_unit1_eef1a1_stim_24h_5race_umi_n_eff',
        'obs_l2_v3_unit1_eef1a1_stim_24h_unexpected_tss_umi_k',
        'obs_l2_v3_unit2_resting_0h_5race_umi_n_eff',
        'obs_l2_v3_unit2_resting_0h_unexpected_tss_umi_k',
        'obs_l2_v3_unit2_stim_24h_5race_umi_n_eff',
        'obs_l2_v3_unit2_stim_24h_unexpected_tss_umi_k',
        'obs_l2_v3_unit2_tss_pos1_in_donor',
        'obs_l2_v3_stage0_eef1a1_cpg_methylation_pct',
        'obs_l2_v3_stage0_eef1a1_cpg_methylation_se',
        'obs_l2_v3_stage0_chs4_cpg_methylation_pct',
        'obs_l2_v3_stage0_chs4_cpg_methylation_se',
        'obs_l2_v3_stage0_dm1_cpg_methylation_pct',
        'obs_l2_v3_stage0_dm1_cpg_methylation_se',
        'obs_l2_v3_stage7_eef1a1_cpg_methylation_pct',
        'obs_l2_v3_stage7_eef1a1_cpg_methylation_se',
        'obs_l2_v3_stage7_chs4_cpg_methylation_pct',
        'obs_l2_v3_stage7_chs4_cpg_methylation_se',
        'obs_l2_v3_stage7_dm1_cpg_methylation_pct',
        'obs_l2_v3_stage7_dm1_cpg_methylation_se',
        'obs_l2_v3_chs4_ctcf_chip_fold_enrichment',
        'obs_l2_v3_chs4_ctcf_chip_fold_enrichment_se',
        'obs_l2_v3_dm1_ctcf_chip_fold_enrichment',
        'obs_l2_v3_dm1_ctcf_chip_fold_enrichment_se',
        'obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment',
        'obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment_se',
        'obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment',
        'obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment_se',
        'obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment',
        'obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment_se',
        'obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment',
        'obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment_se',
        'obs_l2_v3_igg_depleted_serum_adcc_lysis_pct',
        'obs_l2_v3_igg_depleted_serum_adcc_lysis_se',
        'obs_l2_v3_adcc_serum_donor_pair_id',
        'obs_l2_v3_adcc_intact_specimen_id',
        'obs_l2_v3_adcc_igg_depleted_specimen_id',
        'obs_l2_v3_barrier_matched_80bp_comparator_id',
        'obs_l2_v3_barrier_matched_80bp_comparator_sha256',
        'obs_l2_v3_barrier_matched_80bp_comparator_sample_id',
        'obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_pct',
        'obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_se',
        'obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment',
        'obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment_se',
        'obs_tp53_ppm1d_bcor_vaf_n',
        'obs_tp53_ppm1d_bcor_vaf_sample_n',
        'obs_tp53_ppm1d_bcor_vaf_replicates_n',
        'obs_castseq_translocation_freq_n',
        'obs_castseq_translocation_n',
        'obs_castseq_translocation_sample_n',
        'obs_castseq_translocation_replicates_n',
        'obs_icasp9_log10_kill_n',
        'obs_icasp9_log10_kill_sample_n',
        'obs_icasp9_log10_kill_replicates_n',
        'obs_stage0_independent_lot_replicates_n',
        'rnai_attribution_scramble_control_sha256',
        'rnai_attribution_scramble_sample_id',
        'rnai_attribution_intronless_control_sha256',
        'rnai_attribution_intronless_sample_id',
        'ciita_ki_allele1_informative_transcript_recovery_verified',
        'ciita_ki_allele2_informative_transcript_recovery_verified',
        'obs_l2_v3_tandem_polya_readthrough_fold',
        'obs_l2_v3_tandem_polya_readthrough_se',
        'obs_l2_v3_cfh_stp_gpi_quantibrite_molecules',
        'obs_l2_v3_cfh_stp_gpi_quantibrite_se',
        'obs_l2_v3_unit2_copy_expression_quantibrite_or_fold',
        'obs_l2_v3_unit2_copy_expression_se',
        'obs_l2_v3_small_rnaseq_mir375_retention_pct',
        'obs_l2_v3_small_rnaseq_mir375_retention_se',
        'obs_l2_v3_intron_splicing_efficiency_pct',
        'obs_l2_v3_intron_splicing_efficiency_se',
        'obs_l2_v3_cd276_kd_pct',
        'obs_l2_v3_cd276_kd_se',
        'obs_l2_v3_pvr_kd_pct',
        'obs_l2_v3_pvr_kd_se',
        'obs_l2_v3_joint_stratum_nk_lysis_pct',
        'obs_l2_v3_joint_stratum_nk_lysis_se',
        'obs_l2_v3_manf_locus3_nfkb_suppression_ratio',
        'obs_l2_v3_manf_locus3_nfkb_suppression_se',
        'obs_l2_v3_manf_ki67_prolif_pct',
        'obs_l2_v3_manf_ki67_prolif_se',
        'obs_l2_v3_manf_ocr_ratio_vs_ctrl',
        'obs_l2_v3_manf_ocr_ratio_vs_ctrl_se',
        'obs_l2_v3_ccl22_treg_to_effector_chemotaxis_ratio',
        'obs_l2_v3_ccl22_treg_to_effector_chemotaxis_se',
        'obs_l2_v3_ccl22_effector_chemotaxis_fold_ctrl',
        'obs_l2_v3_ccl22_effector_chemotaxis_fold_se',
        'syn_5xisre_cd200_in_design_y',
        'has_intronic_mir_cd276_pvr',
        'rnai_specific_benefit_claimed',
        's11_replicate_n',
        'per_locus_segmented_baf',
        'post_edit_per_locus_segmented_baf',
        'phasing_verified',
        'cis_verified',
        'wt_inducibility_verified',
        'rnai_attribution_verified',
        'ciita_locus2_genotype_copy_class',
        'ciita_ki_allele_informative_transcript_recovery_verified',
        'ciita_ko_allele_informative_transcript_recovery_verified',
        'ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified',
        'lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified',
        'l2_v3_same_background_intronless_control_linked_to_cargo',
        'l2_v3_same_background_scramble_control_linked_to_cargo',
    }
    unknown_keys = sorted(set(observations) - set(spec_by_mean) - set(spec_by_se) - allowed_boolean_and_meta_flags)
    if unknown_keys:
        raise ValueError(f'Unrecognized Stage-0 decision-tree observation keys: {unknown_keys}')
    for _str_meta_k in (
        'donor_seedline_id',
        'device_assembly_id',
        'physical_route_id',
        'assay_campaign_id',
        'reset_prior_safety_stop_target_category',
        'submitter_id',
        'recipient_serum_pair_id',
        'intact_serum_sample_id',
        'igg_depleted_serum_sample_id',
        'positive_adcc_control_id',
        'isotype_adcc_control_id',
        'arm2_depleted_sample_id',
        'contemporaneous_arm3_control_sample_id',
        'sister_assembly_id',
        'cxcl10_tss_donor_id',
        'cxcl10_tss_transcription_unit_id',
        'construct_sha256',
        'post_sort_endocrine_zm447439_regimen',
        'stage7_start_utc',
        'sort_reaggregation_utc',
        'assay_harvest_utc',
        'registered_mixture_id',
        'locked_analysis_plan_id',
        'gate0s_elevated_store_scope_matched_device_lot_id',
        'background_id',
        'mcb_id',
        'l2_v3_donor_id',
        'tier6b_presort_lineage_panel',
        'tier6b_clonogenic_risk_unit',
    ):
        if _str_meta_k in observations and observations[_str_meta_k] is not None:
            _sv = observations[_str_meta_k]
            if not isinstance(_sv, str) or not _sv.strip() or len(_sv) > 120:
                raise ValueError(f'{_str_meta_k} must be a nonempty string of at most 120 characters (IIP-0024 / IIP-0028; DEC-103 / DEC-106)')
    if 'bound_c1_c5_pair_sha256' in observations and observations['bound_c1_c5_pair_sha256'] is not None:
        _bsha = observations['bound_c1_c5_pair_sha256']
        if not isinstance(_bsha, str) or not re.fullmatch(r'[0-9a-f]{64}', _bsha.strip()):
            raise ValueError('bound_c1_c5_pair_sha256 must be a lowercase 64-character SHA-256 hex digest (IIP-0026; DEC-105)')

    # Pre-normalization strict integer check on count, scored-cell, donor, trial, and replicate denominators (IIP-0024; DEC-103):
    non_integer_scored_cell_or_trial_denominator = False
    for _int_k in STAGE0_STRICT_INTEGER_COUNT_AND_DENOMINATOR_KEYS:
        if _int_k in observations and observations[_int_k] not in ('UNAVAILABLE', None):
            _raw_iv = observations[_int_k]
            if isinstance(_raw_iv, bool):
                raise ValueError(f'{_int_k} must be an integer count, not boolean (IIP-0024; DEC-103)')
            _dec_iv = finite_decimal(_raw_iv, _int_k)
            if abs(_dec_iv - round(_dec_iv)) > 1e-9 or _dec_iv < 0.0:
                non_integer_scored_cell_or_trial_denominator = True
    if 'prior_variant_safety_status' in observations and observations['prior_variant_safety_status'] is not None:
        pvs_val = observations['prior_variant_safety_status']
        if not isinstance(pvs_val, str) or pvs_val.strip() not in ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES:
            raise ValueError(
                f"Invalid prior_variant_safety_status {pvs_val!r}; must be one of {sorted(ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES)} (AUDIT-REF; DEC-071)"
            )
    for _adj_k, _adj_allowed in (
        ('adjunct_conditioning_history_arm', ALLOWED_ADJUNCT_CONDITIONING_HISTORY_ARMS),
        ('adjunct_depletion_mode', ALLOWED_ADJUNCT_DEPLETION_MODES),
        ('adjunct_sa_fhl1_status', ALLOWED_ADJUNCT_SA_FHL1_STATUSES),
        ('adjunct_immune_responder_history', ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES),
        ('knockin_splice_challenger_arm', ALLOWED_KNOCKIN_SPLICE_CHALLENGER_ARMS),
        ('knockin_locus4_orientation', ALLOWED_KNOCKIN_LOCUS4_ORIENTATIONS),
        ('s10_hlae_floor_mode', ALLOWED_S10_HLAE_FLOOR_MODES),
        ('b2_tyk2_recipient_genotype_stratum', ALLOWED_B2_TYK2_RECIPIENT_GENOTYPE_STRATA),
        ('b2_tyk2_hlac_allele', ALLOWED_B2_TYK2_HLAC_ALLELES),
        ('s6_cluster_size_arm', ALLOWED_S6_CLUSTER_SIZE_ARMS),
        ('p05_barrier_slab_geometry_mode', ALLOWED_P05_BARRIER_SLAB_GEOMETRY_MODES),
        ('design_z_divider_arm', ALLOWED_DESIGN_Z_DIVIDER_ARMS),
        ('p05_catalase_uniprot_accession', ALLOWED_P05_CATALASE_UNIPROT_ACCESSIONS),
        ('cuff_bioconjugation_arm', ALLOWED_CUFF_BIOCONJUGATION_ARMS),
        ('cuff_serum_challenge_stratum', ALLOWED_CUFF_SERUM_CHALLENGE_STRATA),
        ('rest_ocr_coupon_arm', ALLOWED_REST_OCR_COUPON_ARMS),
        ('conditioning_4arm_arm', ALLOWED_CONDITIONING_4ARM_ARMS),
        ('conditioning_washout_analyte_scope', ALLOWED_CONDITIONING_WASHOUT_ANALYTE_SCOPES),
        ('track2_rescue_route', ALLOWED_TRACK2_RESCUE_ROUTES),
        ('design_y_integrase_modality', ALLOWED_DESIGN_Y_INTEGRASE_MODALITIES),
        ('design_y_silencing_modality', ALLOWED_DESIGN_Y_SILENCING_MODALITIES),
        ('design_y_payload_arm', ALLOWED_DESIGN_Y_PAYLOAD_ARMS),
        ('contrast_b1_geometry_arm', ALLOWED_CONTRAST_B1_GEOMETRY_ARMS),
        ('contrast_b2_perfusion_and_retrieval_arm', ALLOWED_CONTRAST_B2_ARMS),
        ('contrast_b3_donor_iec_arm', ALLOWED_CONTRAST_B3_DONOR_IEC_ARMS),
        ('tyk2_pe_modality', ALLOWED_TYK2_PE_MODALITIES),
        ('iip0026_schema_era', ALLOWED_IIP0026_SCHEMA_ERAS),
    ):
        if _adj_k in observations and observations[_adj_k] is not None:
            _adj_v = observations[_adj_k]
            if not isinstance(_adj_v, str) or _adj_v.strip() not in _adj_allowed:
                raise ValueError(
                    f"Invalid {_adj_k} {_adj_v!r}; must be one of {sorted(_adj_allowed)} (IIP-0012 / IIP-0015 / IIP-0016; DEC-091 / DEC-092 / DEC-094)"
                )

    defaulted_se_set = set(defaulted_se_keys or ())
    validated_stage0: dict[str, float] = {}
    validated_stage0_se: dict[str, float] = {}
    se_provenance_by_mean: dict[str, str] = {}
    for mk, (eid, _mk, sk, lo, hi, _unit) in spec_by_mean.items():
        if mk in observations and observations[mk] != 'UNAVAILABLE' and observations[mk] is not None:
            val = finite_decimal(observations[mk], mk)
            if not lo <= val <= hi:
                raise ValueError(f'{mk} ({eid}) must be in [{lo}, {hi}]')
            validated_stage0[mk] = val
            if sk in observations and observations[sk] not in ('UNAVAILABLE', None) and sk not in defaulted_se_set and mk not in defaulted_se_set:
                se_val = finite_decimal(observations[sk], sk)
                min_se = (
                    contract['minimum_se']
                    if mk in primary_spec_by_mean and primary_spec_by_mean[mk] and default_se_by_mean.get(mk) is None
                    else STAGE0_CUSTOM_MINIMUM_SE.get(mk, contract['minimum_se'])
                )
                if se_val < min_se:
                    raise ValueError(f'{sk} must be >= {min_se}')
                validated_stage0_se[sk] = se_val
                se_provenance_by_mean[mk] = 'EXPLICIT_SE_SUPPLIED'
            else:
                # Unified scale-aware SE fallback across direct evaluate_stage0_decision_tree and normalize_measurements (AUDIT-REF, AUDIT-REF, R7-AUDIT-REF & AUDIT-REF; DEC-041/050/062/071):
                imp_se = (
                    finite_decimal(observations[sk], sk)
                    if (sk in observations and observations[sk] not in ('UNAVAILABLE', None))
                    else STAGE0_CUSTOM_DEFAULT_SE.get(
                        mk,
                        default_se_by_mean.get(mk, _scale_aware_stage0_default_se(mk, lo, hi)),
                    )
                )
                min_se_imp = (
                    contract['minimum_se']
                    if mk in primary_spec_by_mean and primary_spec_by_mean[mk] and default_se_by_mean.get(mk) is None
                    else STAGE0_CUSTOM_MINIMUM_SE.get(mk, contract['minimum_se'])
                )
                if imp_se < min_se_imp:
                    raise ValueError(f'{sk} must be >= {min_se_imp}')
                validated_stage0_se[sk] = imp_se
                se_provenance_by_mean[mk] = 'DEFAULT_SE_IMPUTED_NOT_EXPLICIT'

    def _explicit_se(mean_key: str) -> float:
        if mean_key not in validated_stage0 or se_provenance_by_mean.get(mean_key) != 'EXPLICIT_SE_SUPPLIED':
            return 0.0
        row = spec_by_mean.get(mean_key)
        if not row:
            return 0.0
        return float(validated_stage0_se.get(row[2], 0.0))

    def _validated_se(mean_key: str) -> float:
        if mean_key not in validated_stage0:
            return 0.0
        row = spec_by_mean.get(mean_key)
        if not row:
            return 0.0
        return float(validated_stage0_se.get(row[2], 0.0))

    stage0_replicates_below_6_hold = False
    crit_95_mult = 1.644854
    if 'obs_stage0_independent_lot_replicates_n' in validated_stage0:
        _n_rep_raw = validated_stage0['obs_stage0_independent_lot_replicates_n']
        if abs(_n_rep_raw - round(_n_rep_raw)) > 1e-9 or int(round(_n_rep_raw)) < 6:
            stage0_replicates_below_6_hold = True
        else:
            crit_95_mult = student_t_ppf_one_sided_95(float(int(round(_n_rep_raw)) - 1))

    def _upper95(mean_key: str, default: float = 0.0) -> float:
        if mean_key not in validated_stage0:
            return default
        return float(validated_stage0[mean_key]) + crit_95_mult * _validated_se(mean_key)

    def _lower95(mean_key: str, default: float = 0.0) -> float:
        if mean_key not in validated_stage0:
            return default
        return float(validated_stage0[mean_key]) - crit_95_mult * _validated_se(mean_key)

    # Exact one-sided 95% Student-t multiplier defaulting to df = 5 (n = 6 minimum lot replicates) when replicate count is omitted (IIP-0025; DEC-104):
    _student_t_mult_df5_default = (
        crit_95_mult
        if ('obs_stage0_independent_lot_replicates_n' in validated_stage0 and not stage0_replicates_below_6_hold)
        else student_t_ppf_one_sided_95(5.0)
    )

    def _upper95_student_t(mean_key: str, default: float = 0.0) -> float:
        if mean_key not in validated_stage0:
            return default
        return float(validated_stage0[mean_key]) + _student_t_mult_df5_default * _validated_se(mean_key)

    def _lower95_student_t(mean_key: str, default: float = 0.0) -> float:
        if mean_key not in validated_stage0:
            return default
        return float(validated_stage0[mean_key]) - _student_t_mult_df5_default * _validated_se(mean_key)

    serpinb9_h0 = validated_stage0.get('obs_serpinb9_hour0_uM', validated_stage0.get('obs_gate0c_hour0_serpinb9_uM', 0.35))
    serpinb9_h0_lower95 = (
        _lower95('obs_serpinb9_hour0_uM', serpinb9_h0)
        if 'obs_serpinb9_hour0_uM' in validated_stage0
        else _lower95('obs_gate0c_hour0_serpinb9_uM', serpinb9_h0)
    )
    lysis_2h = validated_stage0.get('obs_regimeII_lysis_2h_pct', validated_stage0.get('obs_gate0c_2h_lysis_pct', 5.0))
    lysis_2h_upper95 = (
        _upper95('obs_regimeII_lysis_2h_pct', lysis_2h)
        if 'obs_regimeII_lysis_2h_pct' in validated_stage0
        else _upper95('obs_gate0c_2h_lysis_pct', lysis_2h)
    )
    chop_fold_candidates = []
    if 'obs_chop_coincidence_fold' in validated_stage0:
        chop_fold_candidates.append(validated_stage0['obs_chop_coincidence_fold'])
    if 'obs_gate0e_chop_log2fc' in validated_stage0:
        chop_fold_candidates.append(2.0 ** validated_stage0['obs_gate0e_chop_log2fc'])
    chop_fold = max(chop_fold_candidates) if chop_fold_candidates else (2.0 ** 0.4)

    hc10_fold_candidates = []
    if 'obs_hc10_surface_fold_wt' in validated_stage0:
        hc10_fold_candidates.append(validated_stage0['obs_hc10_surface_fold_wt'])
    if 'obs_gate0e_hc10_surface_mfi_ratio' in validated_stage0:
        hc10_fold_candidates.append(validated_stage0['obs_gate0e_hc10_surface_mfi_ratio'])
    hc10_fold = max(hc10_fold_candidates) if hc10_fold_candidates else 1.1
    # Default unmeasured capsule thickness to 10.0 um (<= 10.4 um) so S9_Chronic does not falsely fire on omission (AUDIT-REF; DEC-041):
    has_chronic_capsule_obs = ('obs_capsule_dz_um' in validated_stage0) or ('obs_gate0f_capsule_thickness_um' in validated_stage0)
    capsule_dz = validated_stage0.get('obs_capsule_dz_um', validated_stage0.get('obs_gate0f_capsule_thickness_um', 10.0))
    capsule_dz_upper95 = (
        _upper95('obs_capsule_dz_um', capsule_dz)
        if 'obs_capsule_dz_um' in validated_stage0
        else _upper95('obs_gate0f_capsule_thickness_um', capsule_dz)
    )

    has_icasp9_obs = 'obs_icasp9_log10_kill' in validated_stage0
    has_vaf_obs = 'obs_tp53_ppm1d_bcor_vaf' in validated_stage0
    has_cast_obs = 'obs_castseq_translocation_freq' in validated_stage0
    has_chop_obs = ('obs_chop_coincidence_fold' in validated_stage0) or ('obs_gate0e_chop_log2fc' in validated_stage0)
    has_hc10_obs = ('obs_hc10_surface_fold_wt' in validated_stage0) or ('obs_gate0e_hc10_surface_mfi_ratio' in validated_stage0)
    mandatory_s9_safety_panel_complete = bool(
        has_icasp9_obs and has_vaf_obs and has_cast_obs
    )
    full_5endpoint_s9_panel_complete = bool(
        mandatory_s9_safety_panel_complete and has_chop_obs and has_hc10_obs
    )

    # Single shared order-independent Rule S9 iCasp9 / TP53-PPM1D-BCOR VAF / CAST-Seq replicate-count lookup across ALL 4 entry paths (IIP-0031; DEC-110):
    _s9_n_resolved, _s9_n_conflict_groups, _s9_n_any_conflict = _resolve_s9_replicate_counts_detailed_iip0031(
        (validated_stage0, observations)
    )
    _s9_sample_n_default = _s9_n_resolved.get('default')
    _icasp9_sample_n = _s9_n_resolved.get('obs_icasp9_log10_kill')
    if _icasp9_sample_n is None:
        _icasp9_sample_n = _s9_sample_n_default
    _icasp9_t_mult = (
        student_t_ppf_one_sided_95(float(int(round(float(_icasp9_sample_n))) - 1))
        if (
            _icasp9_sample_n is not None
            and type(_icasp9_sample_n) is not bool
            and isinstance(_icasp9_sample_n, (int, float))
            and math.isfinite(float(_icasp9_sample_n))
            and float(_icasp9_sample_n) >= 2.0
        )
        else _student_t_mult_df5_default
    )
    icasp9_lower95 = (
        validated_stage0['obs_icasp9_log10_kill'] - _icasp9_t_mult * validated_stage0_se.get('obs_icasp9_log10_kill_se', 0.0)
        if has_icasp9_obs
        else 4.5
    )
    _vaf_sample_n = _s9_n_resolved.get('obs_tp53_ppm1d_bcor_vaf')
    if _vaf_sample_n is None:
        _vaf_sample_n = _s9_sample_n_default
    _cast_sample_n = _s9_n_resolved.get('obs_castseq_translocation_freq')
    if _cast_sample_n is None:
        _cast_sample_n = _s9_sample_n_default

    tp53_duplex_depth_explicitly_supplied = 'obs_tp53_ppm1d_bcor_duplex_depth' in validated_stage0
    vaf_duplex_depth = validated_stage0.get('obs_tp53_ppm1d_bcor_duplex_depth', 62958.0)
    if has_vaf_obs:
        _vaf_guard = _evaluate_s9_tp53_vaf_guard(
            validated_stage0['obs_tp53_ppm1d_bcor_vaf'],
            validated_stage0_se.get('obs_tp53_ppm1d_bcor_vaf_se', 0.0),
            validated_stage0.get('obs_tp53_ppm1d_bcor_duplex_depth'),
            sample_n=_vaf_sample_n,
            mutant_molecules_k=validated_stage0.get('obs_tp53_ppm1d_bcor_mutant_molecules_k'),
            require_explicit_depth=False,
            endpoint_label='obs_tp53_ppm1d_bcor_vaf',
        )
        vaf_upper95 = float(_vaf_guard['upper95'])
        tp53_duplex_depth_underpowered_below_62958 = bool(
            tp53_duplex_depth_explicitly_supplied and _vaf_guard['depth_underpowered']
        )
    else:
        vaf_upper95 = 0.0
        tp53_duplex_depth_underpowered_below_62958 = bool(
            tp53_duplex_depth_explicitly_supplied and vaf_duplex_depth < 62958.0
        )

    cast_duplex_depth_explicitly_supplied = 'obs_castseq_duplex_or_amplicon_depth' in validated_stage0
    cast_duplex_depth = validated_stage0.get('obs_castseq_duplex_or_amplicon_depth', 62958.0)
    if has_cast_obs:
        _cast_guard = _evaluate_s9_tp53_vaf_guard(
            validated_stage0['obs_castseq_translocation_freq'],
            validated_stage0_se.get('obs_castseq_translocation_se', 0.0),
            validated_stage0.get('obs_castseq_duplex_or_amplicon_depth'),
            sample_n=_cast_sample_n,
            mutant_molecules_k=validated_stage0.get('obs_castseq_mutant_junctions_k'),
            require_explicit_depth=False,
            endpoint_label='obs_castseq_translocation_freq',
        )
        cast_upper95 = float(_cast_guard['upper95'])
        cast_duplex_depth_underpowered_below_62958 = bool(
            cast_duplex_depth_explicitly_supplied and _cast_guard['depth_underpowered']
        )
    else:
        cast_upper95 = 0.0
        cast_duplex_depth_underpowered_below_62958 = bool(
            cast_duplex_depth_explicitly_supplied and cast_duplex_depth < 62958.0
        )
    mcb_clonogenic_trials_explicitly_supplied = 'obs_mcb_effective_clonogenic_trials' in validated_stage0
    mcb_clonogenic_trials = validated_stage0.get('obs_mcb_effective_clonogenic_trials', 0.0)
    mcb_clonogenic_trials_underpowered_below_299573228 = bool(
        mcb_clonogenic_trials_explicitly_supplied and mcb_clonogenic_trials < 299573228.0
    )
    mcb_survivors_k_explicitly_supplied = 'obs_mcb_clonogenic_survivor_count_k' in validated_stage0
    mcb_survivors_k = validated_stage0.get('obs_mcb_clonogenic_survivor_count_k', 0.0)
    mcb_clonogenic_survivor_detected = bool(
        mcb_survivors_k_explicitly_supplied and mcb_survivors_k > 0.0
    )
    mcb_clonogenic_outcome_missing_when_trials_supplied = bool(
        mcb_clonogenic_trials_explicitly_supplied and not mcb_survivors_k_explicitly_supplied
    )
    mcb_clonogenic_trials_missing_when_outcome_zero_supplied = bool(
        mcb_survivors_k_explicitly_supplied and mcb_survivors_k == 0.0 and not mcb_clonogenic_trials_explicitly_supplied
    )

    chop_upper95_candidates = []
    if 'obs_chop_coincidence_fold' in validated_stage0:
        chop_upper95_candidates.append(
            validated_stage0['obs_chop_coincidence_fold'] + 1.645 * validated_stage0_se.get('obs_chop_coincidence_se', 0.0)
        )
    if 'obs_gate0e_chop_log2fc' in validated_stage0:
        chop_upper95_candidates.append(
            2.0 ** (validated_stage0['obs_gate0e_chop_log2fc'] + 1.645 * validated_stage0_se.get('obs_gate0e_chop_log2fc_se', 0.0))
        )
    chop_upper95 = max(chop_upper95_candidates) if chop_upper95_candidates else chop_fold

    hc10_upper95_candidates = []
    if 'obs_hc10_surface_fold_wt' in validated_stage0:
        hc10_upper95_candidates.append(
            validated_stage0['obs_hc10_surface_fold_wt'] + 1.645 * validated_stage0_se.get('obs_hc10_surface_fold_se', 0.0)
        )
    if 'obs_gate0e_hc10_surface_mfi_ratio' in validated_stage0:
        hc10_upper95_candidates.append(
            validated_stage0['obs_gate0e_hc10_surface_mfi_ratio'] + 1.645 * validated_stage0_se.get('obs_gate0e_hc10_surface_mfi_se', 0.0)
        )
    hc10_upper95 = max(hc10_upper95_candidates) if hc10_upper95_candidates else hc10_fold

    # Structured QP governance + verified repeat-batch ledger evidence required to reset a sticky S9 hard safety stop (AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF; DEC-042/043/044/071):
    prior_s9_flag = bool(
        prior_batch_sticky_s9_safety_stop is True
        or observations.get('prior_batch_sticky_s9_safety_stop') is True
        or str(observations.get('prior_variant_safety_status', '')).strip() in STICKY_S9_PRIOR_VARIANT_SAFETY_STATUSES
    )
    repeat_batch_ids = observations.get('reset_prior_safety_stop_repeat_batch_ids')
    verified_ledger_passing = observations.get('verified_ledger_passing_s9_batch_ids')
    cur_batch_id_str = str(observations.get('batch_id', '')).strip()
    structured_s9_reset_valid = bool(
        observations.get('reset_prior_safety_stop_after_investigation') is True
        and (
            (
                isinstance(observations.get('reset_prior_safety_stop_reason_code'), str)
                and observations.get('reset_prior_safety_stop_reason_code') in ALLOWED_STICKY_S9_RESET_REASON_CODES
            )
            or (
                isinstance(observations.get('reset_prior_safety_stop_reason'), str)
                and observations.get('reset_prior_safety_stop_reason') in ALLOWED_STICKY_S9_RESET_REASON_CODES
            )
        )
        and isinstance(repeat_batch_ids, list)
        and len(repeat_batch_ids) >= 2
        and all(isinstance(x, str) and x.strip() for x in repeat_batch_ids)
        and len({x.strip() for x in repeat_batch_ids}) >= 2
        and (not cur_batch_id_str or cur_batch_id_str not in {x.strip() for x in repeat_batch_ids})
        and isinstance(observations.get('reset_prior_safety_stop_independent_qp_approver_id'), str)
        and len(observations['reset_prior_safety_stop_independent_qp_approver_id'].strip()) >= 3
        and (
            not observations.get('submitter_id')
            or str(observations.get('submitter_id')).strip() != observations['reset_prior_safety_stop_independent_qp_approver_id'].strip()
        )
        and isinstance(verified_ledger_passing, (list, tuple, set))
        and {x.strip() for x in repeat_batch_ids}.issubset({str(b).strip() for b in verified_ledger_passing if str(b).strip()})
    )
    sticky_prior_s9_stop = bool(prior_s9_flag and not structured_s9_reset_valid)
    prior_s9_reset_history_retained = bool(prior_s9_flag and structured_s9_reset_valid)

    # Rare-event trial denominators & exact Poisson/Clopper-Pearson upper95 bounds for residual PSC, hazard-state proliferative survivor, Step-1 10 nM AP1903 survivor, and proportion/percentage gates (AUDIT-REF, DEC-050/051/071; IIP-0024; DEC-103):
    def _prop_or_pct_u95(mean_key: str, n_key: str, is_pct: bool = True) -> float:
        nonlocal non_integer_scored_cell_or_trial_denominator
        w_u95 = _upper95(mean_key)
        if mean_key in validated_stage0 and n_key in validated_stage0 and validated_stage0[n_key] > 0.0:
            n_c = validated_stage0[n_key]
            if abs(n_c - round(n_c)) > 1e-9 or round(n_c) <= 0:
                non_integer_scored_cell_or_trial_denominator = True
                return 100.0 if is_pct else 1.0
            n_int = int(round(n_c))
            scale = 100.0 if is_pct else 1.0
            k_ev = (validated_stage0[mean_key] / scale) * float(n_int)
            if abs(k_ev - round(k_ev)) > 1e-6 or round(k_ev) < 0 or round(k_ev) > n_int:
                non_integer_scored_cell_or_trial_denominator = True
                return 100.0 if is_pct else 1.0
            k_int = int(round(k_ev))
            cp_rate = clopper_pearson_upper95_independent(k_int, n_int)
            return max(w_u95, cp_rate * scale)
        return w_u95

    # Clonal mosaicism and Day-28 Gate 0V zero-event counts NEVER accept an explicit SE without an integer cells_scored denominator (IIP-0024; DEC-103):
    strict_scored_cell_zero_keys = frozenset({
        'obs_20q11_bcl2l1_single_cell_mosaic_fraction',
        'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction',
        'obs_1q32_mdm4_single_cell_mosaic_fraction',
        'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct',
    })
    mosaic_explicit_se_missing_scored_cell_denominator = any(
        mk in validated_stage0
        and nk not in validated_stage0
        and se_provenance_by_mean.get(mk) == 'EXPLICIT_SE_SUPPLIED'
        for mk, nk in (
            ('obs_20q11_bcl2l1_single_cell_mosaic_fraction', 'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored'),
            ('obs_mcb_1q32_mdm4_single_cell_mosaic_fraction', 'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored'),
            ('obs_1q32_mdm4_single_cell_mosaic_fraction', 'obs_1q32_mdm4_single_cell_mosaic_cells_scored'),
        )
    )
    zero_event_proportion_missing_denominator_or_explicit_se = any(
        mk in validated_stage0
        and validated_stage0[mk] == 0.0
        and nk not in validated_stage0
        and (
            mk in strict_scored_cell_zero_keys
            or se_provenance_by_mean.get(mk) != 'EXPLICIT_SE_SUPPLIED'
            or _validated_se(mk) < 5e-4
        )
        for mk, nk in (
            ('obs_20q11_bcl2l1_single_cell_mosaic_fraction', 'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored'),
            ('obs_mcb_1q32_mdm4_single_cell_mosaic_fraction', 'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored'),
            ('obs_1q32_mdm4_single_cell_mosaic_fraction', 'obs_1q32_mdm4_single_cell_mosaic_cells_scored'),
            ('obs_17q_birc5_survivin_single_cell_mosaic_fraction', 'obs_17q_birc5_survivin_single_cell_mosaic_cells_scored'),
            ('obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct', 'obs_gate0v_day28_persistent_vp1_ie1_cells_scored'),
            ('obs_cd24_sox2_oct4_3d_escapee_pct', 'obs_cd24_sox2_oct4_3d_escapee_cells_scored'),
            ('obs_off_target_gp2_amylase_acinar_fraction', 'obs_off_target_gp2_amylase_acinar_cells_scored'),
            ('obs_off_target_krt19_cftr_ductal_fraction', 'obs_off_target_krt19_cftr_ductal_cells_scored'),
            ('obs_assembled_cluster_ki67_progenitor_pct', 'obs_assembled_cluster_ki67_progenitor_cells_scored'),
            ('obs_ductal_sox9_krt19_double_pos_pct', 'obs_ductal_sox9_krt19_double_pos_cells_scored'),
            ('obs_assembled_cluster_slc18a1_enterochromaffin_pct', 'obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored'),
            ('obs_senescent_p16_sasp_positive_cell_pct', 'obs_senescent_p16_sasp_cells_scored'),
        )
    )

    residual_psc_trials_supplied = 'obs_residual_psc_effective_trials' in validated_stage0
    residual_psc_trials = validated_stage0.get('obs_residual_psc_effective_trials', 0.0)
    residual_psc_validated_se = _validated_se('obs_residual_psc_effective_fraction')
    residual_psc_trials_underpowered = bool(
        (residual_psc_trials_supplied and residual_psc_trials < 299573228.0)
        or (
            'obs_residual_psc_effective_fraction' in validated_stage0
            and validated_stage0['obs_residual_psc_effective_fraction'] <= 1e-8
            and not residual_psc_trials_supplied
        )
    )
    residual_psc_wald_upper95 = validated_stage0.get('obs_residual_psc_effective_fraction', 0.0) + 1.644854 * residual_psc_validated_se
    residual_psc_poisson_upper95 = (
        _poisson_exact_upper95_rate(
            validated_stage0.get('obs_residual_psc_effective_fraction', 0.0) * residual_psc_trials,
            residual_psc_trials,
        )
        if residual_psc_trials_supplied and residual_psc_trials > 0.0
        else 0.0
    )
    residual_psc_upper95 = max(residual_psc_wald_upper95, residual_psc_poisson_upper95)
    residual_psc_failed = bool(
        'obs_residual_psc_effective_fraction' in validated_stage0
        and (validated_stage0['obs_residual_psc_effective_fraction'] > 1e-8 or residual_psc_upper95 > 1e-8)
    )
    chr20q11_se = (
        validated_stage0_se.get('obs_20q11_bcl2l1_cnv_se', 0.0)
        if 'obs_20q11_bcl2l1_cnv_ratio' in validated_stage0
        else 0.0
    )
    chr20q11_lower95 = validated_stage0.get('obs_20q11_bcl2l1_cnv_ratio', 2.0) - 1.644854 * chr20q11_se
    chr20q11_upper95 = validated_stage0.get('obs_20q11_bcl2l1_cnv_ratio', 2.0) + 1.644854 * chr20q11_se
    chr20q11_bulk_cnv_failed = bool(
        'obs_20q11_bcl2l1_cnv_ratio' in validated_stage0
        and (
            validated_stage0['obs_20q11_bcl2l1_cnv_ratio'] < 1.85
            or validated_stage0['obs_20q11_bcl2l1_cnv_ratio'] > 2.15
            or chr20q11_lower95 < 1.85
            or chr20q11_upper95 > 2.15
        )
    )
    chr20q11_mosaic_upper95 = _prop_or_pct_u95(
        'obs_20q11_bcl2l1_single_cell_mosaic_fraction',
        'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored',
        is_pct=False,
    )
    chr20q11_single_cell_mosaic_failed = bool(
        'obs_20q11_bcl2l1_single_cell_mosaic_fraction' in validated_stage0
        and (
            validated_stage0['obs_20q11_bcl2l1_single_cell_mosaic_fraction'] > 0.01
            or chr20q11_mosaic_upper95 > 0.01
        )
    )
    chr20q11_bcl2l1_cnv_failed = bool(chr20q11_bulk_cnv_failed or chr20q11_single_cell_mosaic_failed)
    hazard_prolif_trials_supplied = 'obs_hazard_state_prolif_effective_trials' in validated_stage0
    hazard_prolif_trials = validated_stage0.get('obs_hazard_state_prolif_effective_trials', 0.0)
    hazard_prolif_trials_underpowered = bool(
        (hazard_prolif_trials_supplied and hazard_prolif_trials < 299573228.0)
        or (
            'obs_hazard_state_prolif_survivor_count_k' in validated_stage0
            and validated_stage0['obs_hazard_state_prolif_survivor_count_k'] == 0.0
            and not hazard_prolif_trials_supplied
        )
    )
    hazard_state_prolif_survivor_detected = bool(
        ('obs_hazard_state_prolif_survivor_count_k' in validated_stage0 and validated_stage0['obs_hazard_state_prolif_survivor_count_k'] > 0.0)
        or (
            'obs_hazard_state_prolif_survivor_count_k' in validated_stage0
            and hazard_prolif_trials_supplied
            and hazard_prolif_trials > 0.0
            and _poisson_exact_upper95_rate(validated_stage0['obs_hazard_state_prolif_survivor_count_k'], hazard_prolif_trials) > 1e-8
            and not hazard_prolif_trials_underpowered
        )
    )
    ap1903_pulsed_pk_se = (
        validated_stage0_se.get('obs_ap1903_pulsed_pk_log10_kill_se', 0.0)
        if 'obs_ap1903_pulsed_pk_log10_kill' in validated_stage0
        else 0.0
    )
    ap1903_pulsed_pk_lower95 = validated_stage0.get('obs_ap1903_pulsed_pk_log10_kill', 4.5) - 1.644854 * ap1903_pulsed_pk_se
    ap1903_mono_pk_se = (
        validated_stage0_se.get('obs_ap1903_monotherapy_pulsed_pk_log10_kill_se', 0.0)
        if 'obs_ap1903_monotherapy_pulsed_pk_log10_kill' in validated_stage0
        else 0.0
    )
    ap1903_mono_pk_lower95 = validated_stage0.get('obs_ap1903_monotherapy_pulsed_pk_log10_kill', 4.5) - 1.644854 * ap1903_mono_pk_se
    ap1903_generic_pk_failed = bool(
        (
            'obs_ap1903_pulsed_pk_log10_kill' in validated_stage0
            and (validated_stage0['obs_ap1903_pulsed_pk_log10_kill'] < 4.0 or ap1903_pulsed_pk_lower95 < 4.0)
        )
        or observations.get('obs_ap1903_pulsed_pk_washout_regimen_verified') is False
    )
    ap1903_mono_pk_failed = bool(
        'obs_ap1903_monotherapy_pulsed_pk_log10_kill' in validated_stage0
        and (validated_stage0['obs_ap1903_monotherapy_pulsed_pk_log10_kill'] < 4.0 or ap1903_mono_pk_lower95 < 4.0)
    )
    ap1903_pulsed_pk_kill_failed = bool(ap1903_generic_pk_failed or ap1903_mono_pk_failed)
    ap1903_step1_trials_supplied = 'obs_ap1903_step1_effective_trials' in validated_stage0
    ap1903_step1_trials = validated_stage0.get('obs_ap1903_step1_effective_trials', 0.0)
    ap1903_step1_trials_underpowered = bool(
        (ap1903_step1_trials_supplied and ap1903_step1_trials < 29957322.8)
        or (
            'obs_ap1903_step1_sublethal_survivor_count_k' in validated_stage0
            and validated_stage0['obs_ap1903_step1_sublethal_survivor_count_k'] == 0.0
            and not ap1903_step1_trials_supplied
        )
    )
    ap1903_step1_sublethal_survivor_detected = bool(
        ('obs_ap1903_step1_sublethal_survivor_count_k' in validated_stage0 and validated_stage0['obs_ap1903_step1_sublethal_survivor_count_k'] > 0.0)
        or (
            'obs_ap1903_step1_sublethal_survivor_count_k' in validated_stage0
            and ap1903_step1_trials_supplied
            and ap1903_step1_trials > 0.0
            and _poisson_exact_upper95_rate(validated_stage0['obs_ap1903_step1_sublethal_survivor_count_k'], ap1903_step1_trials) > 1e-7
            and not ap1903_step1_trials_underpowered
        )
    )
    vascularized_explant_recovery_failed_s9 = bool(
        observations.get('obs_vascularized_explant_en_bloc_recovery_verified') is False
    )
    explant_retrieval_failed = bool(
        ('obs_explant_en_bloc_recovery_fraction' in validated_stage0 and _lower95_student_t('obs_explant_en_bloc_recovery_fraction') < 0.99)
        or ('obs_gate0s_day30_en_bloc_recovery_fraction' in validated_stage0 and _lower95_student_t('obs_gate0s_day30_en_bloc_recovery_fraction') < 0.99)
        or ('obs_gate0s_day90_en_bloc_recovery_fraction' in validated_stage0 and _lower95_student_t('obs_gate0s_day90_en_bloc_recovery_fraction') < 0.99)
        or vascularized_explant_recovery_failed_s9
    )
    post_stress_20q11_outgrowth_failed = bool(
        'obs_post_stress_20q11_bcl2l1_outgrowth_fold' in validated_stage0
        and (
            validated_stage0['obs_post_stress_20q11_bcl2l1_outgrowth_fold'] > 1.0
            or _upper95('obs_post_stress_20q11_bcl2l1_outgrowth_fold') > 1.20
        )
    )
    nutlin3a_p53_functional_failed = bool(
        (
            'obs_nutlin3a_p21_mdm2_fold_induction' in validated_stage0
            and (
                validated_stage0['obs_nutlin3a_p21_mdm2_fold_induction'] < 5.0
                or _lower95('obs_nutlin3a_p21_mdm2_fold_induction') < 5.0
            )
        )
        or (
            'obs_nutlin3a_p21_cdkn1a_fold_induction' in validated_stage0
            and (
                validated_stage0['obs_nutlin3a_p21_cdkn1a_fold_induction'] < 5.0
                or _lower95('obs_nutlin3a_p21_cdkn1a_fold_induction') < 5.0
            )
        )
    )
    s9_causal_failed_categories: list[str] = []
    if (
        (has_icasp9_obs and (validated_stage0['obs_icasp9_log10_kill'] < 4.0 or icasp9_lower95 < 4.0))
        or (has_vaf_obs and (validated_stage0['obs_tp53_ppm1d_bcor_vaf'] >= 1e-4 or vaf_upper95 >= 1e-4))
        or (has_cast_obs and (validated_stage0['obs_castseq_translocation_freq'] >= 1e-4 or cast_upper95 >= 1e-4))
    ):
        s9_causal_failed_categories.append('CORE_S9_VAF_OR_CASTSEQ')
    if (has_chop_obs and (chop_fold >= 2.0 or chop_upper95 >= 2.0)) or (has_hc10_obs and (hc10_fold > 2.0 or hc10_upper95 > 2.0)):
        s9_causal_failed_categories.append('ER_CHOP_HC10')
    if mcb_clonogenic_survivor_detected:
        s9_causal_failed_categories.append('MCB_CLONOGENIC')
    if residual_psc_failed:
        s9_causal_failed_categories.append('RESIDUAL_PSC')
    if chr20q11_bulk_cnv_failed:
        s9_causal_failed_categories.append('CNV_20Q11')
    if chr20q11_single_cell_mosaic_failed:
        s9_causal_failed_categories.append('SINGLE_CELL_MOSAIC_20Q11')
    if post_stress_20q11_outgrowth_failed:
        s9_causal_failed_categories.append('POST_STRESS_20Q11_OUTGROWTH')
    if nutlin3a_p53_functional_failed:
        s9_causal_failed_categories.append('FUNCTIONAL_P53_NUTLIN3A')
    if hazard_state_prolif_survivor_detected:
        s9_causal_failed_categories.append('HAZARD_PROLIF_SURVIVOR')
    if ap1903_generic_pk_failed:
        s9_causal_failed_categories.append('PULSED_PK_KILL')
    # Bank & Stage-7 lot 1q32.1 MDM4 diploid CN + >= 299-cell mosaicism and observed-zero vetoes on Day-28 Gate 0V & Day-90 Gate 0S dual shedding (IIP-0024; DEC-103):
    mcb_mdm4_cnv_failed = bool(
        'obs_mcb_1q32_mdm4_cnv_ratio' in validated_stage0
        and (
            validated_stage0['obs_mcb_1q32_mdm4_cnv_ratio'] < 1.85
            or validated_stage0['obs_mcb_1q32_mdm4_cnv_ratio'] > 2.15
            or _lower95('obs_mcb_1q32_mdm4_cnv_ratio') < 1.85
            or _upper95('obs_mcb_1q32_mdm4_cnv_ratio') > 2.15
        )
    )
    lot_mdm4_cnv_failed = bool(
        'obs_1q32_mdm4_cnv_ratio' in validated_stage0
        and (
            validated_stage0['obs_1q32_mdm4_cnv_ratio'] < 1.85
            or validated_stage0['obs_1q32_mdm4_cnv_ratio'] > 2.15
            or _lower95('obs_1q32_mdm4_cnv_ratio') < 1.85
            or _upper95('obs_1q32_mdm4_cnv_ratio') > 2.15
        )
    )
    mcb_mdm4_mosaic_u95 = _prop_or_pct_u95(
        'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction',
        'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored',
        is_pct=False,
    )
    lot_mdm4_mosaic_u95 = _prop_or_pct_u95(
        'obs_1q32_mdm4_single_cell_mosaic_fraction',
        'obs_1q32_mdm4_single_cell_mosaic_cells_scored',
        is_pct=False,
    )
    mdm4_mosaic_failed = bool(
        ('obs_mcb_1q32_mdm4_single_cell_mosaic_fraction' in validated_stage0 and (validated_stage0['obs_mcb_1q32_mdm4_single_cell_mosaic_fraction'] > 0.01 or mcb_mdm4_mosaic_u95 > 0.01))
        or ('obs_1q32_mdm4_single_cell_mosaic_fraction' in validated_stage0 and (validated_stage0['obs_1q32_mdm4_single_cell_mosaic_fraction'] > 0.01 or lot_mdm4_mosaic_u95 > 0.01))
    )
    gate0v_day28_observed_positive_veto = bool(
        'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct' in validated_stage0
        and validated_stage0['obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct'] > 0.0
    )
    gate0s_shedding_observed_positive_veto = bool(
        ('obs_gate0s_shedding_ddpcr_events_k' in validated_stage0 and validated_stage0['obs_gate0s_shedding_ddpcr_events_k'] > 0.0)
        or ('obs_gate0s_shedding_viable_events_k' in validated_stage0 and validated_stage0['obs_gate0s_shedding_viable_events_k'] > 0.0)
    )
    p05_z70_interface_h2o2_hard_stop_veto = bool(
        'obs_p05_z70_gel_interface_h2o2_max_uM' in validated_stage0
        and (
            validated_stage0['obs_p05_z70_gel_interface_h2o2_max_uM'] > 5.0
            or _upper95_student_t('obs_p05_z70_gel_interface_h2o2_max_uM') > 5.0
        )
    )
    if ap1903_mono_pk_failed:
        s9_causal_failed_categories.append('MONOTHERAPY_PULSED_PK_KILL')
    if ap1903_step1_sublethal_survivor_detected:
        s9_causal_failed_categories.append('STEP1_SUBLETHAL_SURVIVOR')
    if explant_retrieval_failed:
        s9_causal_failed_categories.append('EXPLANT_RECOVERY')
    if mcb_mdm4_cnv_failed or lot_mdm4_cnv_failed:
        s9_causal_failed_categories.append('CNV_1Q32_MDM4')
    if mdm4_mosaic_failed:
        s9_causal_failed_categories.append('SINGLE_CELL_MOSAIC_1Q32_MDM4')
    if gate0v_day28_observed_positive_veto:
        s9_causal_failed_categories.append('GATE0V_PERSISTENT_VIRAL_POSITIVE')
    if gate0s_shedding_observed_positive_veto:
        s9_causal_failed_categories.append('GATE0S_DUAL_SHEDDING_POSITIVE')
    if p05_z70_interface_h2o2_hard_stop_veto:
        s9_causal_failed_categories.append('P05_Z70_INTERFACE_H2O2_GT_5UM')

    current_batch_s9_endpoint_failed = bool(s9_causal_failed_categories)
    safety_kill_s9_triggered = bool(
        sticky_prior_s9_stop
        or current_batch_s9_endpoint_failed
    )

    nkg2a_explicit_se = _explicit_se('obs_nkg2a_lysis_pct')
    kir_explicit_se = _explicit_se('obs_kir_nkg2a_neg_nk_lysis_pct')
    nkg2c_explicit_se = _explicit_se('obs_nkg2c_adaptive_nk_lysis_pct')
    il12_il18_explicit_se = _explicit_se('obs_il12_il18_sirpa_low_nk_lysis_pct')
    phago_explicit_se = _explicit_se('obs_phagocytosis_index_pct')
    receptor_neg_nkg2c_explicit_se = _explicit_se('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct')
    c3b_phago_explicit_se = _explicit_se('obs_complement_c3b_opsonized_phagocytosis_pct')
    nkg2a_lysis_upper95 = _upper95('obs_nkg2a_lysis_pct', 2.0)
    kir_nkg2a_neg_upper95 = _upper95('obs_kir_nkg2a_neg_nk_lysis_pct', 2.0)
    nkg2c_adaptive_upper95 = _upper95('obs_nkg2c_adaptive_nk_lysis_pct', 2.0)
    il12_il18_sirpa_low_upper95 = _upper95('obs_il12_il18_sirpa_low_nk_lysis_pct', 2.0)
    phagocytosis_upper95 = _upper95('obs_phagocytosis_index_pct', 3.0)
    receptor_neg_nkg2c_upper95 = _upper95('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct', 2.0)
    c3b_phagocytosis_upper95 = _upper95('obs_complement_c3b_opsonized_phagocytosis_pct', 3.0)

    # Enforce explicit SE caps (<= 3.0% for NK lysis panels, <= 2.5% for macrophage phagocytosis) AND compare one-sided 95% UCBs directly against the biological release gates (8.0% lysis, 10.0% IL-12/IL-18 lysis, 5.0% Playbook phagocytosis; MW-R5-S2-F01 & AUDIT-REF; DEC-048/049):
    immune_se_cap_exceeded = bool(
        nkg2a_explicit_se > 3.0
        or kir_explicit_se > 3.0
        or nkg2c_explicit_se > 3.0
        or il12_il18_explicit_se > 3.0
        or receptor_neg_nkg2c_explicit_se > 3.0
        or phago_explicit_se > 2.5
        or c3b_phago_explicit_se > 2.5
    )
    excessive_immune_donor_variance_or_upper95_failed = bool(
        immune_se_cap_exceeded
        or ('obs_nkg2a_lysis_pct' in validated_stage0 and validated_stage0['obs_nkg2a_lysis_pct'] <= 8.0 and nkg2a_lysis_upper95 > 8.0)
        or ('obs_kir_nkg2a_neg_nk_lysis_pct' in validated_stage0 and validated_stage0['obs_kir_nkg2a_neg_nk_lysis_pct'] <= 8.0 and kir_nkg2a_neg_upper95 > 8.0)
        or ('obs_nkg2c_adaptive_nk_lysis_pct' in validated_stage0 and validated_stage0['obs_nkg2c_adaptive_nk_lysis_pct'] <= 8.0 and nkg2c_adaptive_upper95 > 8.0)
        or ('obs_il12_il18_sirpa_low_nk_lysis_pct' in validated_stage0 and validated_stage0['obs_il12_il18_sirpa_low_nk_lysis_pct'] <= 10.0 and il12_il18_sirpa_low_upper95 > 10.0)
        or ('obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct' in validated_stage0 and validated_stage0['obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct'] <= 8.0 and receptor_neg_nkg2c_upper95 > 8.0)
        or ('obs_phagocytosis_index_pct' in validated_stage0 and validated_stage0['obs_phagocytosis_index_pct'] <= 5.0 and phagocytosis_upper95 > 5.0)
        or ('obs_complement_c3b_opsonized_phagocytosis_pct' in validated_stage0 and validated_stage0['obs_complement_c3b_opsonized_phagocytosis_pct'] <= 5.0 and c3b_phagocytosis_upper95 > 5.0)
    )

    receptor_neg_nkg2c_nk_failed = bool(
        (
            'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct' in validated_stage0
            and receptor_neg_nkg2c_upper95 > 8.0
        )
        or (
            'obs_nkp46_innate_nk_lysis_pct' in validated_stage0
            and _upper95('obs_nkp46_innate_nk_lysis_pct') > 8.0
        )
    )
    post_stress_nkg2d_dnam1_lysis_failed = bool(
        'obs_post_stress_nkg2d_dnam1_attributed_lysis_pct' in validated_stage0
        and _upper95('obs_post_stress_nkg2d_dnam1_attributed_lysis_pct') > 5.0
    )
    s1_promoted = bool(
        (
            serpinb9_h0_lower95 < 0.30
            or lysis_2h_upper95 > 8.0
            or _upper95('obs_12pulse_unprimed_hypoxic_serial_lysis_pct', 5.0) > 8.0
            or receptor_neg_nkg2c_nk_failed
            or post_stress_nkg2d_dnam1_lysis_failed
            or _upper95('obs_serpinb9_si', 1.2) > 1.5
            or _lower95('obs_serpinb9_half_life_h', 6.0) < 4.0
            or _lower95('obs_serpinb9_kass_M_s', 5e6) < 1e6
            or _lower95('obs_serpinb9_log10_kass', 6.7) < 6.0
            or observations.get('s1_gzmb_kass_or_pdl1_scar_failed') is True
        )
        and lysis_2h_upper95 <= 25.0
    )
    catastrophic_nkg2a_lysis_exceeds_arm1x = bool(
        'obs_nkg2a_lysis_pct' in validated_stage0
        and (
            validated_stage0['obs_nkg2a_lysis_pct'] > 15.0
            or nkg2a_lysis_upper95 > 15.0
        )
    )
    # Uncertainty-aware upper95 KIR+NKG2A- (<= 8.0%) and IL-12/IL-18 SIRPa-low (<= 10.0%) escape routing & fail-closed Track-2 B'-SCT vs HLA-G1 split (IIP-0022; DEC-101):
    # Preserves provisional core_14_panel compatibility when subset endpoints are unmeasured while preventing raw-mean <= 8.0% with upper95 > 8.0% from misrouting to S2 Arm 1X.
    kir_nkg2a_neg_nk_exceeded_raw = bool(
        ('obs_kir_nkg2a_neg_nk_lysis_pct' in validated_stage0 and (validated_stage0['obs_kir_nkg2a_neg_nk_lysis_pct'] > 8.0 or kir_nkg2a_neg_upper95 > 8.0))
        or ('obs_il12_il18_sirpa_low_nk_lysis_pct' in validated_stage0 and (validated_stage0['obs_il12_il18_sirpa_low_nk_lysis_pct'] > 10.0 or il12_il18_sirpa_low_upper95 > 10.0))
    )
    _t2_route = observations.get('track2_rescue_route')
    if isinstance(_t2_route, str):
        _t2_route = _t2_route.strip()
    _legacy_conflated_t2 = observations.get('track2_b_sct_or_hlag1_lilrb1_rescue_qualified') is True
    _bsct_flag = observations.get('track2_b_sct_kir_rescue_qualified') is True
    _hlag1_flag = observations.get('track2_hlag1_lilrb1_rescue_qualified') is True
    _bsct_quant_present = any(
        _k in validated_stage0
        for _k in ('obs_track2_bsct_kir_nkg2a_neg_nk_lysis_pct', 'obs_track2_bsct_donor_count', 'obs_track2_bsct_failing_donor_count')
    )
    _hlag1_quant_present = any(
        _k in validated_stage0
        for _k in ('obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct', 'obs_track2_hlag1_donor_count', 'obs_track2_hlag1_failing_donor_count')
    )
    _bsct_quant_cleared = bool(
        'obs_track2_bsct_kir_nkg2a_neg_nk_lysis_pct' in validated_stage0
        and se_provenance_by_mean.get('obs_track2_bsct_kir_nkg2a_neg_nk_lysis_pct') == 'EXPLICIT_SE_SUPPLIED'
        and _upper95('obs_track2_bsct_kir_nkg2a_neg_nk_lysis_pct') <= 8.0
        and 'obs_track2_bsct_donor_count' in validated_stage0
        and abs(validated_stage0['obs_track2_bsct_donor_count'] - round(validated_stage0['obs_track2_bsct_donor_count'])) <= 1e-9
        and int(round(validated_stage0['obs_track2_bsct_donor_count'])) >= 6
        and 'obs_track2_bsct_failing_donor_count' in validated_stage0
        and abs(validated_stage0['obs_track2_bsct_failing_donor_count'] - round(validated_stage0['obs_track2_bsct_failing_donor_count'])) <= 1e-9
        and int(round(validated_stage0['obs_track2_bsct_failing_donor_count'])) == 0
    )
    _hlag1_quant_cleared_uniform_8pct = bool(
        'obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct' in validated_stage0
        and se_provenance_by_mean.get('obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct') == 'EXPLICIT_SE_SUPPLIED'
        and _upper95('obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct') <= 8.0
        and 'obs_track2_hlag1_donor_count' in validated_stage0
        and abs(validated_stage0['obs_track2_hlag1_donor_count'] - round(validated_stage0['obs_track2_hlag1_donor_count'])) <= 1e-9
        and int(round(validated_stage0['obs_track2_hlag1_donor_count'])) >= 6
        and 'obs_track2_hlag1_failing_donor_count' in validated_stage0
        and abs(validated_stage0['obs_track2_hlag1_failing_donor_count'] - round(validated_stage0['obs_track2_hlag1_failing_donor_count'])) <= 1e-9
        and int(round(validated_stage0['obs_track2_hlag1_failing_donor_count'])) == 0
    )
    # Deprecated conflated alias preserves legacy B'-SCT fallback candidate promotion ONLY when neither HLA-G1 route nor failing B'-SCT quant panel is present; never auto-clears HLA-G1!
    track2_bsct_route_promoted = bool(
        kir_nkg2a_neg_nk_exceeded_raw
        and not _hlag1_flag
        and _t2_route != 'TRACK2_HLAG1_LILRB1_CHALLENGER_HELD'
        and (
            (_bsct_flag and (not _bsct_quant_present or _bsct_quant_cleared))
            or (_legacy_conflated_t2 and not _bsct_flag and _t2_route is None and not _bsct_quant_present)
        )
    )
    track2_hlag1_route_attempted = bool(
        _hlag1_flag or _t2_route == 'TRACK2_HLAG1_LILRB1_CHALLENGER_HELD' or _hlag1_quant_present
    )
    track2_kir_or_hlag1_fallback_promoted = bool(track2_bsct_route_promoted)
    kir_nkg2a_neg_nk_failed = bool(
        kir_nkg2a_neg_nk_exceeded_raw
        and not track2_bsct_route_promoted
    )
    # S10 triggers on one-sided 95% bounds (CMV-imprinted NKG2C+NKG2A-CD57+ adaptive NK lysis upper95 > 8% OR receptor-negative NKG2A-KIR2DL-NKG2C+ NK lysis upper95 > 8% OR covalent-SCT NKG2C/NKG2A KD selectivity ratio lower95 < 4.0) or any failing S10 v3 functional/comparator panel (AUDIT-REF & IIP-0016; DEC-042/046/094):
    _s10_v3_active = bool(
        any(observations.get(_mk) is not None for _mk in STAGE0_S10_V3_METADATA_KEYS if _mk != 'current_lot_id')
        or any(observations.get(_bk) is not None for _bk in STAGE0_S10_V3_BOOLEAN_KEYS)
        or any(_ek in validated_stage0 for _ek in STAGE0_S10_V3_MEAN_KEYS)
    )
    _has_explicit_s10_obs = bool(
        _s10_v3_active
        or 'obs_nkg2c_adaptive_nk_lysis_pct' in validated_stage0
        or 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct' in validated_stage0
        or 'obs_hlae_nkg2c_over_nkg2a_kd_ratio' in validated_stage0
        or observations.get('s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct') is not None
    )
    s10_v3_eval: dict | None = None
    s10_v3_hold_status: str | None = None
    if _s10_v3_active:
        _v3_obs = dict(observations)
        for _mk_v, _val_v in validated_stage0.items():
            _v3_obs[_mk_v] = _val_v
        for _sk_v, _sval_v in validated_stage0_se.items():
            _v3_obs[_sk_v] = _sval_v
        s10_v3_eval = evaluate_s10_functional_and_replacement_gate_v3(
            _v3_obs,
            STAGE0_S10_V3_CONSTRUCT_REGISTRY,
            defaulted_se_keys=defaulted_se_set,
        )
        if not s10_v3_eval['gate_pass']:
            s10_v3_hold_status = str(s10_v3_eval['status'])
    s10_kd_ratio_lower95_sub4 = bool(_lower95('obs_hlae_nkg2c_over_nkg2a_kd_ratio', 4.5) < 4.0)
    s10_nkg2c_triggered = bool(
        ('obs_nkg2c_adaptive_nk_lysis_pct' in validated_stage0 and nkg2c_adaptive_upper95 > 8.0)
        or receptor_neg_nkg2c_nk_failed
        or s10_kd_ratio_lower95_sub4
        or (_s10_v3_active and s10_v3_eval is not None and not s10_v3_eval['gate_pass'])
    )
    s10_coupon_cleared_legacy = bool(
        observations.get('s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct') is True
        and not receptor_neg_nkg2c_nk_failed
        and ('obs_nkg2c_adaptive_nk_lysis_pct' not in validated_stage0 or nkg2c_adaptive_upper95 <= 8.0)
    )
    # Unified S10 clearance binding: when S10 v3 panel/metadata is active, clearance is strictly bound to s10_v3_eval['gate_pass'] (defeating spoofed legacy coupon=True); bare legacy coupon=True unmasks S2 screening only and cannot clear lot_release_readiness_tier without v3 LOT_RELEASE_25K qualification (IIP-0016; DEC-094):
    s10_coupon_cleared = bool(
        s10_v3_eval['gate_pass']
        if (_s10_v3_active and s10_v3_eval is not None)
        else s10_coupon_cleared_legacy
    )
    s10_v3_release_qualified = bool(
        _s10_v3_active
        and s10_v3_eval is not None
        and s10_v3_eval['gate_pass'] is True
        and str(observations.get('s10_hlae_floor_mode', 'GATE_0A_SCREENING_15K')).strip() == 'LOT_RELEASE_25K'
        and not receptor_neg_nkg2c_nk_failed
        and ('obs_nkg2c_adaptive_nk_lysis_pct' not in validated_stage0 or nkg2c_adaptive_upper95 <= 8.0)
    )
    # S2 triggers on low resting HLA-E (<15,000 molecules/cell) OR moderate CD94/NKG2A+ NK lysis in (8.0%, 15.0%], strictly when neither catastrophic (>15%) NKG2A+ lysis nor CD94-independent KIR+NKG2A- / IL-12/IL-18 SIRPa-low NK failure nor uncleared S10 NKG2C peptide-selectivity failure is active (MW-R5-S2-F01/F02 & AUDIT-REF; DEC-048/049):
    s2_arm1b_failed = bool(
        (
            _lower95('obs_hlae_hour0_molecules', 20000.0) < 15000.0
            or (8.0 < validated_stage0.get('obs_nkg2a_lysis_pct', 2.0) <= 15.0)
        )
        and not catastrophic_nkg2a_lysis_exceeds_arm1x
        and not kir_nkg2a_neg_nk_exceeded_raw
        and not (s10_nkg2c_triggered and not s10_coupon_cleared)
    )
    has_all_3_dxs_endpoints = bool(
        'obs_dxs_locus4_joint_retention_fraction' in validated_stage0
        and 'obs_dxs_cd47low_icasp9low_fraction' in validated_stage0
        and 'obs_dxs_ap1903_cd47high_escape_fraction' in validated_stage0
    )
    has_explicit_dxs_clonogenic_factors = bool(
        'obs_dxs_screened_cell_count' in validated_stage0
        and 'obs_dxs_clonogenic_fraction' in validated_stage0
        and 'obs_dxs_spike_in_detection_efficiency' in validated_stage0
    )
    dxs_screened_n = validated_stage0.get('obs_dxs_screened_cell_count', 0.0)
    dxs_f_clonogenic = validated_stage0.get('obs_dxs_clonogenic_fraction', 0.0)
    dxs_eta_detection = validated_stage0.get('obs_dxs_spike_in_detection_efficiency', 0.0)
    dxs_effective_clonogenic_n = (
        dxs_screened_n * dxs_f_clonogenic * dxs_eta_detection
        if has_explicit_dxs_clonogenic_factors
        else 0.0
    )
    # Rare-event DX-S verification ALWAYS requires explicit obs_dxs_screened_cell_count, explicit obs_dxs_clonogenic_fraction (> 0), explicit obs_dxs_spike_in_detection_efficiency (>= 0.30), and effective clonogenic trial count >= 299,573,228 (no silent f=1.0, eta=1.0 defaults; MW-R5-S2-F03; DEC-048):
    dxs_explicit_se_or_cell_count_supplied = bool(
        has_all_3_dxs_endpoints
        and has_explicit_dxs_clonogenic_factors
        and dxs_f_clonogenic > 0.0
        and dxs_eta_detection >= 0.30
        and dxs_effective_clonogenic_n >= 299573228.0
    )
    dxs_retention_lower95 = (
        min(
            validated_stage0['obs_dxs_locus4_joint_retention_fraction']
            - 1.645 * validated_stage0_se.get('obs_dxs_locus4_joint_retention_se', 0.0),
            1.0
            - _poisson_exact_upper95_rate(
                (1.0 - validated_stage0['obs_dxs_locus4_joint_retention_fraction']) * dxs_effective_clonogenic_n,
                dxs_effective_clonogenic_n,
            )
            if dxs_effective_clonogenic_n > 0.0
            else 1.0,
        )
        if 'obs_dxs_locus4_joint_retention_fraction' in validated_stage0
        else 0.0
    )
    dxs_double_loss_upper95 = (
        max(
            validated_stage0['obs_dxs_cd47low_icasp9low_fraction']
            + 1.645 * validated_stage0_se.get('obs_dxs_cd47low_icasp9low_se', 0.0),
            _poisson_exact_upper95_rate(
                validated_stage0['obs_dxs_cd47low_icasp9low_fraction'] * dxs_effective_clonogenic_n,
                dxs_effective_clonogenic_n,
            )
            if dxs_effective_clonogenic_n > 0.0
            else 0.0,
        )
        if 'obs_dxs_cd47low_icasp9low_fraction' in validated_stage0
        else 1.0
    )
    dxs_escape_upper95 = (
        max(
            validated_stage0['obs_dxs_ap1903_cd47high_escape_fraction']
            + 1.645 * validated_stage0_se.get('obs_dxs_ap1903_cd47high_escape_se', 0.0),
            _poisson_exact_upper95_rate(
                validated_stage0['obs_dxs_ap1903_cd47high_escape_fraction'] * dxs_effective_clonogenic_n,
                dxs_effective_clonogenic_n,
            )
            if dxs_effective_clonogenic_n > 0.0
            else 0.0,
        )
        if 'obs_dxs_ap1903_cd47high_escape_fraction' in validated_stage0
        else 1.0
    )
    is_family_a_or_unspecified = bool(
        product_def is None or str(product_def.get('family', 'A')).upper().startswith('A')
    )
    dxs_joint_retention_verified = bool(
        is_family_a_or_unspecified
        and has_all_3_dxs_endpoints
        and dxs_explicit_se_or_cell_count_supplied
        and validated_stage0['obs_dxs_locus4_joint_retention_fraction'] >= 0.9999
        and dxs_retention_lower95 >= 0.9999
        and validated_stage0['obs_dxs_cd47low_icasp9low_fraction'] <= 1e-6
        and dxs_double_loss_upper95 <= 1e-6
        and validated_stage0['obs_dxs_ap1903_cd47high_escape_fraction'] <= 1e-8
        and dxs_escape_upper95 <= 1e-8
    )
    dxs_hold_missing_explicit_se_or_counts = bool(
        has_all_3_dxs_endpoints and not dxs_explicit_se_or_cell_count_supplied
    )
    iec_vascular_support_failed = bool(
        ('obs_iec_vascular_support_viability_pct' in validated_stage0 and _lower95('obs_iec_vascular_support_viability_pct') < 85.0)
        or ('obs_iec_cd47_p10_molecules' in validated_stage0 and _lower95('obs_iec_cd47_p10_molecules') < 180000.0)
        or ('obs_iec_abo_ab_antigen_mfi_ratio' in validated_stage0 and _upper95('obs_iec_abo_ab_antigen_mfi_ratio') > 1.05)
        or ('obs_vascular_pericyte_ng2_pdgfrb_coverage_pct' in validated_stage0 and _lower95('obs_vascular_pericyte_ng2_pdgfrb_coverage_pct') < 70.0)
        or ('obs_iec_anastomosis_t50_h' in validated_stage0 and _upper95('obs_iec_anastomosis_t50_h') > 36.0)
    )
    effero_clearance_failed = bool(
        ('obs_apoptotic_effero_clearance_retention_pct_of_wt' in validated_stage0 and _lower95('obs_apoptotic_effero_clearance_retention_pct_of_wt') < 70.0)
        or ('obs_effero_clearance_pct' in validated_stage0 and _lower95('obs_effero_clearance_pct') < 80.0)
    )
    vp1_ie1_upper95 = _prop_or_pct_u95(
        'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct',
        'obs_gate0v_day28_persistent_vp1_ie1_cells_scored',
        is_pct=True,
    )
    gate0v_viral_sanctuary_failed = bool(
        ('obs_gate0v_cvb_day14_log10_viral_reduction' in validated_stage0 and _lower95('obs_gate0v_cvb_day14_log10_viral_reduction') < 3.0)
        or ('obs_cvb4_10moi_72h_viral_titer_log10_pfu' in validated_stage0 and _upper95('obs_cvb4_10moi_72h_viral_titer_log10_pfu') > 2.0)
        or (
            'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct' in validated_stage0
            and (
                validated_stage0['obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct'] > 0.0
                or vp1_ie1_upper95 > 0.05
            )
        )
    )
    assembled_cassette_chandler_ibmir_failed = bool(
        ('obs_assembled_cassette_chandler_loop_sc5b9_fold' in validated_stage0 and _upper95('obs_assembled_cassette_chandler_loop_sc5b9_fold') > 2.0)
        or ('obs_assembled_cassette_chandler_loop_platelet_retention_pct' in validated_stage0 and _lower95('obs_assembled_cassette_chandler_loop_platelet_retention_pct') < 70.0)
        or ('obs_c3b_cd55_cd59_complement_deposition_mfi_ratio' in validated_stage0 and _upper95('obs_c3b_cd55_cd59_complement_deposition_mfi_ratio') > 1.5)
    )
    slc18a1_ec_upper95 = _prop_or_pct_u95(
        'obs_assembled_cluster_slc18a1_enterochromaffin_pct',
        'obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored',
        is_pct=True,
    )
    sox9_krt19_upper95 = _prop_or_pct_u95(
        'obs_ductal_sox9_krt19_double_pos_pct',
        'obs_ductal_sox9_krt19_double_pos_cells_scored',
        is_pct=True,
    )
    ki67_prog_upper95 = _prop_or_pct_u95(
        'obs_assembled_cluster_ki67_progenitor_pct',
        'obs_assembled_cluster_ki67_progenitor_cells_scored',
        is_pct=True,
    )
    assembled_cluster_lineage_or_identity_failed = bool(
        ('obs_assembled_cluster_slc18a1_enterochromaffin_pct' in validated_stage0 and slc18a1_ec_upper95 > 1.0)
        or ('obs_enterochromaffin_fev_tph1_pos_pct' in validated_stage0 and _upper95('obs_enterochromaffin_fev_tph1_pos_pct') > 1.0)
        or ('obs_sc_atac_ec_chromatin_spillover_pct' in validated_stage0 and _upper95('obs_sc_atac_ec_chromatin_spillover_pct') > 5.0)
        or ('obs_endocrine_synaptophysin_pos_purity_pct' in validated_stage0 and _lower95('obs_endocrine_synaptophysin_pos_purity_pct') < 95.0)
        or ('obs_ductal_sox9_krt19_double_pos_pct' in validated_stage0 and sox9_krt19_upper95 > 0.1)
        or ('obs_assembled_cluster_ins_gcg_polyhormonal_pct' in validated_stage0 and _upper95('obs_assembled_cluster_ins_gcg_polyhormonal_pct') > 2.0)
        or ('obs_assembled_cluster_ki67_progenitor_pct' in validated_stage0 and ki67_prog_upper95 > 0.5)
        or ('obs_post_stress_mafa_nkx61_beta_identity_retention_pct' in validated_stage0 and _lower95('obs_post_stress_mafa_nkx61_beta_identity_retention_pct') < 85.0)
        or ('obs_pouch_iapp_oligomer_a11_dot_blot_fold' in validated_stage0 and _upper95('obs_pouch_iapp_oligomer_a11_dot_blot_fold') > 1.5)
    )
    cd47_or_phagocytosis_failed = bool(
        observations.get('s3_arm1b_and_arm1x_below_15k_and_arm1c_splicing_ge_0p95') is True
        or _lower95('obs_cd47_hour0_molecules', 220000.0) < 50000.0
        or validated_stage0.get('obs_phagocytosis_index_pct', 3.0) > 10.0
        or _lower95('obs_cd47_cc2c6_over_b6h12_ratio', 0.95) < 0.90
        or validated_stage0.get('obs_complement_c3b_opsonized_phagocytosis_pct', 3.0) > 5.0
    )
    # Separate soft resting CD47 200k / 5% phagocytosis band from non-waivable mandatory post-stress/hypoxic/single-cell p10 and retention floors (R7-AUDIT-REF; DEC-061):
    resting_cd47_200k_or_5pct_phago_soft_cliff_failed = bool(
        ('obs_cd47_hour0_molecules' in validated_stage0 and 50000.0 <= _lower95('obs_cd47_hour0_molecules') < 200000.0)
        or ('obs_phagocytosis_index_pct' in validated_stage0 and 5.0 < validated_stage0['obs_phagocytosis_index_pct'] <= 10.0)
    )
    mandatory_hypoxic_or_p10_retention_cliff_failed = bool(
        (
            'obs_cd47_p10_single_cell_molecules' in validated_stage0
            and (
                validated_stage0['obs_cd47_p10_single_cell_molecules'] < 180000.0
                or _lower95('obs_cd47_p10_single_cell_molecules') < 180000.0
            )
        )
        or (
            'obs_hypoxic_72h_cd47_p10_single_cell_molecules' in validated_stage0
            and (
                validated_stage0['obs_hypoxic_72h_cd47_p10_single_cell_molecules'] < 180000.0
                or _lower95('obs_hypoxic_72h_cd47_p10_single_cell_molecules') < 180000.0
            )
        )
        or (
            'obs_hypoxic_72h_hlae_p10_single_cell_molecules' in validated_stage0
            and (
                validated_stage0['obs_hypoxic_72h_hlae_p10_single_cell_molecules'] < 15000.0
                or _lower95('obs_hypoxic_72h_hlae_p10_single_cell_molecules') < 15000.0
            )
        )
        or (
            'obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction' in validated_stage0
            and (
                validated_stage0['obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction'] < 0.90
                or _lower95('obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction') < 0.90
            )
        )
        or (
            'obs_hypoxic_72h_cd47_retention_molecules' in validated_stage0
            and (
                validated_stage0['obs_hypoxic_72h_cd47_retention_molecules'] < 180000.0
                or _lower95('obs_hypoxic_72h_cd47_retention_molecules') < 180000.0
            )
        )
        or (
            'obs_hypoxic_72h_hlae_retention_molecules' in validated_stage0
            and (
                validated_stage0['obs_hypoxic_72h_hlae_retention_molecules'] < 15000.0
                or _lower95('obs_hypoxic_72h_hlae_retention_molecules') < 15000.0
            )
        )
        or (
            'obs_hypoxic_15mmhg_f_isr_translation_ratio' in validated_stage0
            and (
                validated_stage0['obs_hypoxic_15mmhg_f_isr_translation_ratio'] < 0.50
                or _lower95('obs_hypoxic_15mmhg_f_isr_translation_ratio') < 0.50
            )
        )
    )
    cd47_or_phago_playbook_cliff_failed = bool(
        resting_cd47_200k_or_5pct_phago_soft_cliff_failed
        or mandatory_hypoxic_or_p10_retention_cliff_failed
    )
    unmitigated_cd47_or_hypoxic_retention_cliff_failed = bool(
        not dxs_joint_retention_verified
        and (
            mandatory_hypoxic_or_p10_retention_cliff_failed
            or (
                resting_cd47_200k_or_5pct_phago_soft_cliff_failed
                and observations.get('cd47_playbook_200k_and_5pct_phago_waiver_approved') is not True
            )
        )
    )
    # S3 is evaluated independently of S2 (no `not s2_promoted` masking; AUDIT-REF; DEC-047):
    s3_promoted = bool(
        dxs_joint_retention_verified
        and (cd47_or_phagocytosis_failed or cd47_or_phago_playbook_cliff_failed)
    )
    s2_promoted = bool(
        (
            observations.get('s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior') is True
            or s2_arm1b_failed
        )
        and not catastrophic_nkg2a_lysis_exceeds_arm1x
        and not kir_nkg2a_neg_nk_exceeded_raw
        and not (s10_nkg2c_triggered and not s10_coupon_cleared)
        and not s3_promoted
    )
    lot_ocr_ratio = validated_stage0.get('obs_lot_ocr_ratio', 1.0)
    lot_ocr_upper95 = _upper95('obs_lot_ocr_ratio', lot_ocr_ratio)
    lot_ocr_exceeds_max_1p50x_ceiling = bool('obs_lot_ocr_ratio' in validated_stage0 and lot_ocr_upper95 > 1.50)
    s4_promoted = bool(
        observations.get('s4_rnls_wt_clears_gsis_shutoff_and_ocr_and_er_gate') is True
        or lot_ocr_upper95 > 1.097
    )
    s5_promoted = bool(
        observations.get('s5_fasl_delta_tunel_gt_5pp') is True
        or ('obs_safasl_delta_tunel_pp' in validated_stage0 and _upper95('obs_safasl_delta_tunel_pp') > 5.0)
        or ('obs_safasl_zb4_anti_fas_rescue_pp' in validated_stage0 and _upper95('obs_safasl_zb4_anti_fas_rescue_pp') > 5.0)
        or ('obs_safasl_treg_viable_retention_pct_of_thbd' in validated_stage0 and _lower95('obs_safasl_treg_viable_retention_pct_of_thbd') < 85.0)
        or ('obs_safasl_neutrophil_chemotaxis_index' in validated_stage0 and _upper95('obs_safasl_neutrophil_chemotaxis_index') > 1.5)
        or ('obs_sa_fasl_microgel_nk_degranulation_cd107a_pct' in validated_stage0 and _upper95('obs_sa_fasl_microgel_nk_degranulation_cd107a_pct') > 5.0)
    )
    exogenous_b2m_class1_rescue_failed = bool(
        ('obs_exogenous_b2m_w632_surface_fold_isotype' in validated_stage0 and _upper95('obs_exogenous_b2m_w632_surface_fold_isotype') > 1.2)
        or ('obs_exogenous_b2m_t1d_ctl_lysis_pct' in validated_stage0 and _upper95('obs_exogenous_b2m_t1d_ctl_lysis_pct') > 5.0)
    )
    s6_promoted = bool(
        ('obs_gate0e_chop_log2fc' in validated_stage0 and _upper95('obs_gate0e_chop_log2fc', 0.4) > math.log2(1.5))
        or chop_fold > 1.5
        or (has_chop_obs and chop_upper95 > 1.5)
        or _lower95('obs_gate0e_gsis_retention_ratio', 0.92) < 0.85
        or hc10_fold > 1.5
        or ('obs_hc10_surface_fold_wt' in validated_stage0 and _upper95('obs_hc10_surface_fold_wt') > 1.5)
        or ('obs_gate0e_hc10_surface_mfi_ratio' in validated_stage0 and _upper95('obs_gate0e_hc10_surface_mfi_ratio') > 1.5)
        or exogenous_b2m_class1_rescue_failed
    )
    b_prime_sct_cd8_upper95 = (
        validated_stage0['obs_b_prime_sct_cd8_lysis_pct'] + 1.645 * validated_stage0_se.get('obs_b_prime_sct_cd8_lysis_se', 0.0)
        if 'obs_b_prime_sct_cd8_lysis_pct' in validated_stage0
        else 0.0
    )
    b_prime_sct_direct_cd8_failed = bool(
        (
            'obs_b_prime_sct_cd8_lysis_pct' in validated_stage0
            and (
                validated_stage0['obs_b_prime_sct_cd8_lysis_pct'] > 5.0
                or b_prime_sct_cd8_upper95 > 10.0
                or (
                    'obs_plan_a_star_cd8_lysis_pct' in validated_stage0
                    and (validated_stage0['obs_b_prime_sct_cd8_lysis_pct'] - validated_stage0['obs_plan_a_star_cd8_lysis_pct']) > 5.0
                )
            )
        )
        or (
            'obs_b_prime_sct_21d_primed_cd8_elispot_sfu_per_1e6' in validated_stage0
            and _upper95('obs_b_prime_sct_21d_primed_cd8_elispot_sfu_per_1e6') > 50.0
        )
    )
    s7_allo_triggered = bool(
        _upper95('obs_gate0b_hlac_cd8_elispot_sfw', 10.0) > 20.0
        or _upper95('obs_foreign_hlac_cd8_lysis_pct', 0.0) > 10.0
        or (
            'obs_foreign_hlac_cd8_lysis_pct' in validated_stage0
            and 'obs_plan_a_star_cd8_lysis_pct' in validated_stage0
            and (validated_stage0['obs_foreign_hlac_cd8_lysis_pct'] - validated_stage0['obs_plan_a_star_cd8_lysis_pct']) > 5.0
        )
        or b_prime_sct_direct_cd8_failed
    )
    sct_density_explicitly_qualified = bool(
        'obs_gate0b_hlac_sct_molecules_per_cell' in validated_stage0
        and _lower95('obs_gate0b_hlac_sct_molecules_per_cell') >= 15000.0
    )
    sct_cd8_mlc_explicitly_cleared = bool(
        not b_prime_sct_direct_cd8_failed
        and (
            observations.get('b_prime_sct_recipient_cd8_mlc_cleared') is True
            or ('obs_b_prime_sct_cd8_lysis_pct' in validated_stage0 and validated_stage0['obs_b_prime_sct_cd8_lysis_pct'] <= 5.0 and b_prime_sct_cd8_upper95 <= 10.0)
        )
    )
    _strict_s7_mlc_and_se_required = bool(
        product_def is not None
        or observations.get('b_prime_sct_recipient_cd8_mlc_cleared') is False
        or observations.get('require_full_iip0026_molecular_provenance_panel') is True
    )
    _sct_density_se_explicit = bool(
        se_provenance_by_mean.get('obs_gate0b_hlac_sct_molecules_per_cell') == 'EXPLICIT_SE_SUPPLIED'
    )
    _has_explicit_s7_hlac_obs = bool(
        any(
            _sk in validated_stage0
            for _sk in (
                'obs_gate0b_hlac_cd8_elispot_sfw',
                'obs_foreign_hlac_cd8_lysis_pct',
                'obs_gate0b_hlac_sct_molecules_per_cell',
                'obs_b_prime_sct_cd8_lysis_pct',
            )
        )
        or observations.get('b_prime_sct_recipient_cd8_mlc_cleared') is not None
    )
    s7_promoted = bool(
        s7_allo_triggered
        and sct_density_explicitly_qualified
        and not b_prime_sct_direct_cd8_failed
        and (not _strict_s7_mlc_and_se_required or (sct_cd8_mlc_explicitly_cleared and _sct_density_se_explicit))
    )
    # Design B2 (ISLET1_B2_5LOCUS_TYK2_PE) retained-HLA-C TYK2 p.Pro1104Ala (SAPTKFLEL) empirical OR-hold gate (IIP-0017; DEC-095):
    _b2_tyk2_product_selected = bool(
        product_def is not None
        and str(product_def.get('product_variant_id', '')) == 'ISLET1_B2_5LOCUS_TYK2_PE'
    )
    _b2_tyk2_stratum = (
        observations['b2_tyk2_recipient_genotype_stratum'].strip()
        if isinstance(observations.get('b2_tyk2_recipient_genotype_stratum'), str)
        else None
    )
    _b2_tyk2_allele = (
        observations['b2_tyk2_hlac_allele'].strip()
        if isinstance(observations.get('b2_tyk2_hlac_allele'), str)
        else None
    )
    _b2_tyk2_gate_active = bool(
        _b2_tyk2_product_selected
        or _b2_tyk2_stratum is not None
        or _b2_tyk2_allele is not None
        or any(
            observations.get(_bk) is not None
            for _bk in (
                'b2_tyk2_saptkflel_allele_assigned_ms_detected',
                'b2_tyk2_peptide_pulsed_positive_control_verified',
                'b2_tyk2_hlac_blockade_abrogation_verified',
            )
        )
        or any(
            _bmk in validated_stage0
            for _bmk in (
                'obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp',
                'obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6',
                'obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp',
            )
        )
    )
    b2_tyk2_gate_status: str | None = None
    b2_tyk2_hold_status: str | None = None
    b2_tyk2_eval_summary: dict[str, object] | None = None
    if _b2_tyk2_gate_active:
        _b2_sc_hazard = bool(
            'obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp' in validated_stage0
            and _upper95('obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp') > 5.0
        )
        _b2_iec_hazard = bool(
            'obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp' in validated_stage0
            and _upper95('obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp') > 5.0
        )
        _b2_elispot_hazard = bool(
            'obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6' in validated_stage0
            and _upper95('obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6') > 50.0
        )
        _b2_ms_flag = observations.get('b2_tyk2_saptkflel_allele_assigned_ms_detected')
        _b2_pos_ctrl = observations.get('b2_tyk2_peptide_pulsed_positive_control_verified')
        _b2_block_ok = observations.get('b2_tyk2_hlac_blockade_abrogation_verified')
        _b2_quant_powered_complete = bool(
            _b2_tyk2_stratum in ALLOWED_B2_TYK2_RECIPIENT_GENOTYPE_STRATA
            and _b2_tyk2_allele in ALLOWED_B2_TYK2_HLAC_ALLELES
            and _b2_ms_flag is False
            and all(
                _bmk in validated_stage0 and se_provenance_by_mean.get(_bmk) == 'EXPLICIT_SE_SUPPLIED'
                for _bmk in (
                    'obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp',
                    'obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6',
                    'obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp',
                )
            )
            and validated_stage0_se.get('obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_se', 999.0) <= 2.5
            and validated_stage0_se.get('obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_se', 999.0) <= 2.5
            and validated_stage0_se.get('obs_b2_tyk2_p1104a_primed_cd8_elispot_se', 999.0) <= 20.0
        )
        if _b2_sc_hazard or _b2_iec_hazard or _b2_elispot_hazard:
            b2_tyk2_gate_status = 'HOLD_DESIGN_B2_TYK2_P1104A_HLAC_NEOANTIGEN_HAZARD_DETECTED'
            b2_tyk2_hold_status = b2_tyk2_gate_status
        elif _b2_ms_flag is True:
            b2_tyk2_gate_status = 'HOLD_DESIGN_B2_TYK2_P1104A_HLAC_MS_DETECTED_OR_DISCORDANT_REQUIRES_STRATIFIED_REVIEW'
            b2_tyk2_hold_status = b2_tyk2_gate_status
        elif _b2_pos_ctrl is False or _b2_block_ok is False:
            b2_tyk2_gate_status = 'HOLD_DESIGN_B2_TYK2_P1104A_POSITIVE_CONTROL_OR_CAUSAL_RESTRICTION_UNVERIFIED'
            b2_tyk2_hold_status = b2_tyk2_gate_status
        elif not _b2_quant_powered_complete:
            b2_tyk2_gate_status = 'HOLD_DESIGN_B2_TYK2_P1104A_INCOMPLETE_OR_UNDERPOWERED_STRATIFIED_PANEL'
            b2_tyk2_hold_status = b2_tyk2_gate_status
        elif _b2_pos_ctrl is not True or _b2_block_ok is not True:
            b2_tyk2_gate_status = 'HOLD_DESIGN_B2_TYK2_P1104A_POSITIVE_CONTROL_OR_CAUSAL_RESTRICTION_UNVERIFIED'
            b2_tyk2_hold_status = b2_tyk2_gate_status
        else:
            b2_tyk2_gate_status = 'CLEARED_DESIGN_B2_TYK2_P1104A_STRATIFIED_HLAC_GATE'
        b2_tyk2_eval_summary = {
            'gate_pass': b2_tyk2_hold_status is None,
            'status': b2_tyk2_gate_status,
            'recipient_genotype_stratum': _b2_tyk2_stratum,
            'hlac_allele': _b2_tyk2_allele,
            'saptkflel_allele_assigned_ms_detected': _b2_ms_flag,
            'peptide_pulsed_positive_control_verified': _b2_pos_ctrl,
            'hlac_blockade_abrogation_verified': _b2_block_ok,
        }
    s7_hold_unqualified_sct = bool(s7_allo_triggered and not sct_density_explicitly_qualified)
    s7_hold_unverified_sct_allo = bool(
        (
            s7_allo_triggered
            and sct_density_explicitly_qualified
            and (not sct_cd8_mlc_explicitly_cleared or (_strict_s7_mlc_and_se_required and not _sct_density_se_explicit))
        )
        or b_prime_sct_direct_cd8_failed
    )
    glucagon_2mM_lower95 = _lower95('obs_gate0m_glucagon_2mM_fold_increase', 2.5)
    prolonged_insulin_clearance_failed = bool(
        'obs_recipient_insulin_clearance_tau_min' in validated_stage0
        and _upper95('obs_recipient_insulin_clearance_tau_min') >= 8.0
        and observations.get('s8_staged_dose_and_omental_portal_drainage_qualified') is not True
    )
    s8_promoted = bool(
        _upper95('obs_gate0m_shutoff_t90_min', 2.4) > 3.0
        or ('obs_gate0m_8mM_to_2mM_insulin_shutoff_ratio' in validated_stage0 and _upper95('obs_gate0m_8mM_to_2mM_insulin_shutoff_ratio') > 0.20)
        or _upper95('obs_basal_insulin_frac_peak', 0.10) > 0.15
        or _lower95('obs_glucose_threshold_mM', 5.0) < 4.5
        or _lower95('obs_immaturity_corner_min_glucose_mM', 3.6) < 3.5
        or _lower95('obs_alpha_glucagon_low_glc_fold', 2.5) < 2.0
        or validated_stage0.get('obs_gate0m_glucagon_2mM_fold_increase', 2.5) < 2.0
        or glucagon_2mM_lower95 < 2.0
        or _lower95('obs_alpha_counterregulatory_gain_ratio', 1.0) < 0.80
        or _lower95('obs_ucn3_maturity_ratio', 1.0) < 0.70
        or prolonged_insulin_clearance_failed
        or observations.get('s8_fasting_leak_or_ucn3_failed') is True
    )
    _adj_arm_raw = observations.get('adjunct_conditioning_history_arm')
    adjunct_conditioning_history_arm = (
        _adj_arm_raw.strip() if isinstance(_adj_arm_raw, str) else 'ARM3_CONTINUOUSLY_ACTIVE'
    )
    adjunct_depletion_mode = (
        observations['adjunct_depletion_mode'].strip()
        if isinstance(observations.get('adjunct_depletion_mode'), str)
        else None
    )
    adjunct_sa_fhl1_status = (
        observations['adjunct_sa_fhl1_status'].strip()
        if isinstance(observations.get('adjunct_sa_fhl1_status'), str)
        else None
    )
    adjunct_immune_responder_history = (
        observations['adjunct_immune_responder_history'].strip()
        if isinstance(observations.get('adjunct_immune_responder_history'), str)
        else None
    )
    arm1_never_conditioned_reference_hold = bool(adjunct_conditioning_history_arm == 'ARM1_NEVER_CONDITIONED')
    arm3_ongoing_bioactivity_explicit_false = bool(
        adjunct_conditioning_history_arm == 'ARM3_CONTINUOUSLY_ACTIVE'
        and observations.get('obs_arm3_ongoing_adjunct_bioactivity_verified') is False
    )
    has_complete_arm2_6endpoint_panel = bool(
        'obs_indirect_apc_dq2_dq8_cd4_stimulation_index' in validated_stage0
        and 'obs_indirect_apc_cd4_ifng_si' in validated_stage0
        and 'obs_gate0d_transwell_beta_annexin_v_pct' in validated_stage0
        and 'obs_arm2_post_depletion_viable_beta_recovery_pct' in validated_stage0
        and 'obs_arm2_post_depletion_perifusion_gsis_si' in validated_stage0
        and 'obs_arm2_post_depletion_nhs_c5b9_mac_pct' in validated_stage0
    )
    if adjunct_conditioning_history_arm == 'ARM2_CONDITIONED_DEPLETED':
        # Exempt active-FasL >=80% CD4 apoptosis requirement on verified functionally depleted Arm 2 while enforcing shared post-depletion injury, activation, complement, recovery, and dynamic perifusion GSIS gates (IIP-0012; DEC-091):
        indirect_apc_cd4_failed = False
        arm2_depletion_verification_incomplete = not (
            observations.get('obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified') is True
            and observations.get('obs_arm2_sa_thbd_bioactivity_depleted_verified') is True
            and observations.get('obs_arm2_contemporaneous_sham_handling_matched') is True
            and observations.get('obs_arm3_ongoing_adjunct_bioactivity_verified') is True
            and adjunct_depletion_mode in ('NATURAL_AGING_DECAY', 'PHYSICAL_WASHOUT')
            and adjunct_sa_fhl1_status in ALLOWED_ADJUNCT_SA_FHL1_STATUSES
            and adjunct_immune_responder_history in ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES
        )
        arm2_rebound_injury_or_function_failed = bool(
            ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index' in validated_stage0 and _upper95('obs_indirect_apc_dq2_dq8_cd4_stimulation_index') > 1.5)
            or ('obs_indirect_apc_cd4_ifng_si' in validated_stage0 and _upper95('obs_indirect_apc_cd4_ifng_si') > 2.0)
            or ('obs_gate0d_transwell_beta_annexin_v_pct' in validated_stage0 and _upper95('obs_gate0d_transwell_beta_annexin_v_pct') > 15.0)
            or ('obs_arm2_post_depletion_viable_beta_recovery_pct' in validated_stage0 and _lower95('obs_arm2_post_depletion_viable_beta_recovery_pct') < 85.0)
            or ('obs_arm2_post_depletion_perifusion_gsis_si' in validated_stage0 and _lower95('obs_arm2_post_depletion_perifusion_gsis_si') < 2.5)
            or ('obs_arm2_post_depletion_nhs_c5b9_mac_pct' in validated_stage0 and _upper95('obs_arm2_post_depletion_nhs_c5b9_mac_pct') > 5.0)
        )
        arm2_post_depletion_panel_incomplete = bool(
            not arm2_depletion_verification_incomplete
            and not arm2_rebound_injury_or_function_failed
            and not has_complete_arm2_6endpoint_panel
        )
    else:
        arm2_depletion_verification_incomplete = False
        arm2_rebound_injury_or_function_failed = False
        arm2_post_depletion_panel_incomplete = False
        indirect_apc_cd4_failed = bool(
            ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index' in validated_stage0 and _upper95('obs_indirect_apc_dq2_dq8_cd4_stimulation_index') > 1.5)
            or ('obs_indirect_apc_cd4_ifng_si' in validated_stage0 and _upper95('obs_indirect_apc_cd4_ifng_si') > 2.0)
            or (
                adjunct_conditioning_history_arm == 'ARM3_CONTINUOUSLY_ACTIVE'
                and 'obs_gate0d_activated_cd4_fasl_apoptosis_pct' in validated_stage0
                and _lower95('obs_gate0d_activated_cd4_fasl_apoptosis_pct') < 80.0
            )
            or ('obs_gate0d_transwell_beta_annexin_v_pct' in validated_stage0 and _upper95('obs_gate0d_transwell_beta_annexin_v_pct') > 15.0)
        )

    # Two-strand knock-in splice-challenger Stage-0 gate evaluation (IIP-0015; DEC-092):
    _chall_arm_raw = observations.get('knockin_splice_challenger_arm')
    knockin_splice_challenger_arm = (
        _chall_arm_raw.strip() if isinstance(_chall_arm_raw, str) else 'PARENTAL_REFERENCE'
    )
    knockin_locus4_orientation = (
        observations['knockin_locus4_orientation'].strip()
        if isinstance(observations.get('knockin_locus4_orientation'), str)
        else None
    )
    _chall_any_input = bool(
        _chall_arm_raw is not None
        or knockin_locus4_orientation is not None
        or any(
            observations.get(_ck) is not None
            for _ck in (
                'obs_challenger_ivt_rna_recovery_standard_verified',
                'obs_challenger_pam_seed_shield_coconversion_verified',
                'obs_challenger_host_isoforms_above_loq_verified',
                'obs_challenger_sv40pa_and_host_3prime_cleavage_checked',
            )
        )
        or any(
            _ck in validated_stage0
            for _ck in (
                'obs_challenger_ivt_relative_recovery_eta',
                'obs_challenger_umi_molecules_profiled',
                'obs_challenger_aberrant_splice_umi_count_k',
                'obs_challenger_host_exon_splice_retention_ratio',
                'obs_challenger_dsrna_j2_or_isg_fold_wt',
            )
        )
    )
    challenger_umi_eff = (
        validated_stage0['obs_challenger_umi_molecules_profiled']
        * min(1.0, validated_stage0['obs_challenger_ivt_relative_recovery_eta'])
        if (
            'obs_challenger_umi_molecules_profiled' in validated_stage0
            and 'obs_challenger_ivt_relative_recovery_eta' in validated_stage0
        )
        else 0.0
    )
    challenger_hold_status: str | None = None
    if knockin_splice_challenger_arm != 'PARENTAL_REFERENCE' or _chall_any_input:
        _is_active_chall = bool(knockin_splice_challenger_arm != 'PARENTAL_REFERENCE')
        _tost_z95 = 1.6448536269514722
        _host_ret_mean = validated_stage0.get('obs_challenger_host_exon_splice_retention_ratio')
        _host_ret_se = _explicit_se('obs_challenger_host_exon_splice_retention_ratio')
        _host_ret_se_supplied = (
            se_provenance_by_mean.get('obs_challenger_host_exon_splice_retention_ratio') == 'EXPLICIT_SE_SUPPLIED'
        )
        _host_ret_tost_ok = bool(
            _host_ret_mean is not None
            and _host_ret_se_supplied
            and (_host_ret_mean - _tost_z95 * _host_ret_se) >= 0.80
            and (_host_ret_mean + _tost_z95 * _host_ret_se) <= 1.25
        )
        _orient_mismatch = bool(
            (
                knockin_splice_challenger_arm.startswith('S3_L4_CONVERGENT_')
                and knockin_locus4_orientation != 'CONVERGENT_PLUS'
            )
            or (
                knockin_splice_challenger_arm in (
                    'S3_L4_WINA_ONLY_3457',
                    'S3_L4_EEF230C_ONLY_3457',
                    'S3_L4_SPLICE_230C_3457',
                    'S3_L4_CODIR_7NT_230C_3457',
                )
                and knockin_locus4_orientation == 'CONVERGENT_PLUS'
            )
        )
        if (
            knockin_splice_challenger_arm in UNSHIELDED_REFERENCE_SPLICE_CHALLENGER_ARMS
            or observations.get('obs_challenger_pam_seed_shield_coconversion_verified') is False
            or (_is_active_chall and observations.get('obs_challenger_pam_seed_shield_coconversion_verified') is not True)
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_UNSHIELDED_OR_PAM_SEED_COCONVERSION_UNVERIFIED'
        elif _orient_mismatch:
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_LOCUS4_ORIENTATION_METADATA_MISMATCH'
        elif (
            observations.get('obs_challenger_ivt_rna_recovery_standard_verified') is False
            or (_is_active_chall and observations.get('obs_challenger_ivt_rna_recovery_standard_verified') is not True)
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_PREEXTRACTION_IVT_RNA_RECOVERY_STANDARD_UNVERIFIED'
        elif (
            observations.get('obs_challenger_host_isoforms_above_loq_verified') is False
            or (_is_active_chall and observations.get('obs_challenger_host_isoforms_above_loq_verified') is not True)
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_HOST_ISOFORM_BELOW_LOQ_OR_UNVERIFIED'
        elif (
            observations.get('obs_challenger_sv40pa_and_host_3prime_cleavage_checked') is False
            or (_is_active_chall and observations.get('obs_challenger_sv40pa_and_host_3prime_cleavage_checked') is not True)
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_SV40PA_AND_HOST_3PRIME_CLEAVAGE_UNCHECKED'
        elif (
            'obs_challenger_aberrant_splice_umi_count_k' in validated_stage0
            and validated_stage0['obs_challenger_aberrant_splice_umi_count_k'] > 0.0
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_ABERRANT_CHIMERIC_UMI_DETECTED'
        elif (
            ('obs_challenger_umi_molecules_profiled' in validated_stage0 or _is_active_chall)
            and (
                'obs_challenger_umi_molecules_profiled' not in validated_stage0
                or 'obs_challenger_ivt_relative_recovery_eta' not in validated_stage0
                or 'obs_challenger_aberrant_splice_umi_count_k' not in validated_stage0
                or challenger_umi_eff < 3689.0
            )
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_EFFECTIVE_UMI_BELOW_3689_RECOVERY_ADJUSTED_MINIMUM'
        elif ('obs_challenger_host_exon_splice_retention_ratio' in validated_stage0 or _is_active_chall) and not _host_ret_tost_ok:
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_HOST_EXON_RETENTION_OUTSIDE_90PCT_TOST_0P80_1P25'
        elif (
            ('obs_challenger_dsrna_j2_or_isg_fold_wt' in validated_stage0 or _is_active_chall)
            and (
                'obs_challenger_dsrna_j2_or_isg_fold_wt' not in validated_stage0
                or _upper95('obs_challenger_dsrna_j2_or_isg_fold_wt') > 1.5
            )
        ):
            challenger_hold_status = 'HOLD_CHALLENGER_SPLICE_CONVERGENT_DSRNA_J2_OR_ISG_EXCEEDS_1P5X_CEILING'

    # Locus 3 Mode-B Multi-miR(IL1R1+FAS) 6-arm challenger qualification gate (IIP-0018; DEC-096):
    _l3_multimir_mean_keys = (
        'obs_l3_multimir_spliced_fraction',
        'obs_l3_multimir_5prime_seed_fidelity_pct',
        'obs_l3_multimir_ago2_mir375_retention_pct',
        'obs_l3_multimir_shield_protein_retention_pct',
        'obs_l3_multimir_min_receptor_knockdown_pct',
        'obs_l3_multimir_matched_gsis_output_retention_pct',
    )
    _l3_multimir_gate_active = bool(
        knockin_splice_challenger_arm in LOCUS3_MULTIMIR_MODEB_ALL_ARMS
        or any(
            observations.get(_lmk) is not None
            for _lmk in (
                'l3_multimir_pre_aav6_2x2_transient_kd_falsifier_cleared',
                'l3_multimir_low_ifn_and_resting_shield_retention_verified',
                'l3_multimir_arm6_c341s_hlac_and_dq2_dq8_peptide_screen_cleared',
            )
        )
        or any(_lmk in validated_stage0 for _lmk in _l3_multimir_mean_keys)
    )
    l3_multimir_gate_status: str | None = None
    l3_multimir_hold_status: str | None = None
    l3_multimir_eval_summary: dict[str, object] | None = None
    if _l3_multimir_gate_active:
        _l3_is_arm6 = bool(knockin_splice_challenger_arm == 'ARM6_L3_V2_INTRON_DUPLEX_IL1R1_FAS_C341S_TCA_4398')
        _l3_is_design_b = bool(
            product_def is not None
            and str(product_def.get('family', 'A')).upper().startswith('B')
        )
        _l3_quant_complete = all(
            _lmk in validated_stage0 and se_provenance_by_mean.get(_lmk) == 'EXPLICIT_SE_SUPPLIED'
            for _lmk in _l3_multimir_mean_keys
        )
        if knockin_splice_challenger_arm in LOCUS3_MULTIMIR_MODEB_CONTROL_ARMS:
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_FACTORIAL_CONTROL_ARM_NOT_RELEASE_CHALLENGER'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif knockin_splice_challenger_arm not in LOCUS3_MULTIMIR_MODEB_DUPLEX_CHALLENGER_ARMS:
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_DUPLEX_CHALLENGER_ARM_ID_REQUIRED'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif observations.get('l3_multimir_pre_aav6_2x2_transient_kd_falsifier_cleared') is not True:
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_PRE_AAV6_2X2_TRANSIENT_FALSIFIER_UNVERIFIED'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif observations.get('l3_multimir_low_ifn_and_resting_shield_retention_verified') is not True:
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_RESTING_OR_LOW_IFN_SHIELD_RETENTION_UNVERIFIED'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif (
            (_l3_is_arm6 or observations.get('l3_multimir_arm6_c341s_hlac_and_dq2_dq8_peptide_screen_cleared') is False)
            and observations.get('l3_multimir_arm6_c341s_hlac_and_dq2_dq8_peptide_screen_cleared') is not True
        ):
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_ARM6_C341S_HLAC_OR_DQ2_DQ8_PEPTIDE_SCREEN_UNVERIFIED'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif not _l3_quant_complete:
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_INCOMPLETE_OR_DEFAULT_SE_QUALIFICATION_PANEL'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif (
            _lower95('obs_l3_multimir_spliced_fraction') < 0.85
            or _lower95('obs_l3_multimir_5prime_seed_fidelity_pct') < 85.0
            or _lower95('obs_l3_multimir_ago2_mir375_retention_pct') < 90.0
        ):
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_SPLICING_DROSHA_SEED_OR_AGO2_MIR375_FAILED'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif (
            _lower95('obs_l3_multimir_shield_protein_retention_pct') < 85.0
            or _lower95('obs_l3_multimir_min_receptor_knockdown_pct') < 70.0
            or _lower95('obs_l3_multimir_matched_gsis_output_retention_pct') < 85.0
        ):
            l3_multimir_gate_status = 'HOLD_L3_MULTIMIR_RECEPTOR_KD_SHIELD_RETENTION_OR_MATCHED_GSIS_FAILED'
            l3_multimir_hold_status = l3_multimir_gate_status
        elif challenger_hold_status is not None:
            l3_multimir_gate_status = challenger_hold_status
            l3_multimir_hold_status = challenger_hold_status
        else:
            l3_multimir_gate_status = 'CLEARED_L3_MULTIMIR_MODEB_STAGE0_CHALLENGER_GATE'
        l3_multimir_eval_summary = {
            'gate_pass': l3_multimir_hold_status is None,
            'status': l3_multimir_gate_status,
            'challenger_arm': knockin_splice_challenger_arm,
            'design_b_retained_hlac_background': _l3_is_design_b,
            'arm6_c341s_peptide_screen_required': _l3_is_arm6,
        }

    # Stage-0 explicit-3D transport qualification gates for Rule S6 lateral spacing / dose-matched 80 um sizing, Readout 5 (Gate P05) wet-sieved microbeads + co-cast catalase barrier films, and Design Z two-stage empty-lumen + loaded closed-state confirmation (IIP-0019; DEC-097; product-bound Gate P05 fix IIP-0025; DEC-104):
    _is_design_z_product = bool(
        product_def is not None
        and str(product_def.get('product_variant_id', '')) == 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2'
    )
    _req_iip0025_full = observations.get('require_full_iip0025_hurdle4_transport_and_challenger_panel') is True
    _s6_arm_val = observations.get('s6_cluster_size_arm')
    if isinstance(_s6_arm_val, str):
        _s6_arm_val = _s6_arm_val.strip()
    _p05_geom_val = observations.get('p05_barrier_slab_geometry_mode')
    if isinstance(_p05_geom_val, str):
        _p05_geom_val = _p05_geom_val.strip()
    _dz_arm_val = observations.get('design_z_divider_arm')
    if isinstance(_dz_arm_val, str):
        _dz_arm_val = _dz_arm_val.strip()

    _s6_mean_keys = (
        'obs_s6_p10_nearest_neighbor_edge_gap_um',
        'obs_s6_cluster_diam_p90_um',
    )
    _s6_divider_mean_keys = (
        'obs_s6_divider_14d_hydrated_po2_retention_pct',
        'obs_s6_hydrogel_shear_modulus_gprime_pa',
    )
    _p05_mean_keys = (
        'obs_p05_microbead_volume_d99_um',
        'obs_p05_catalase_barrier_film_um',
        'obs_p05_islet_plane_h2o2_max_uM',
        'obs_p05_day7_catalase_retention_pct',
        'obs_p05_336h_retained_catalase_iu_333cm2',
    )
    _dz_empty_mean_keys = (
        'obs_design_z_empty_lumen_5site_min_po2_mmhg',
        'obs_design_z_vessel_standoff_p90_um',
    )
    _dz_loaded_mean_keys = (
        'obs_design_z_loaded_surface_po2_mmhg',
        'obs_design_z_loaded_core_po2_mmhg',
        'obs_design_z_loaded_perifusion_gsis_si',
    )

    _divider_qc_active = bool(
        _dz_arm_val == 'EXPLORATORY_18PCT_DIVIDER_140UM'
        or observations.get('s6_divider_iso10993_leachables_cleared') is not None
        or any(_mk in validated_stage0 for _mk in _s6_divider_mean_keys)
    )
    _divider_qc_cleared = bool(
        observations.get('s6_divider_iso10993_leachables_cleared') is True
        and all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _s6_divider_mean_keys
        )
        and _lower95('obs_s6_divider_14d_hydrated_po2_retention_pct') >= 90.0
        and _upper95('obs_s6_hydrogel_shear_modulus_gprime_pa') <= 1000.0
    )

    _s6_active = bool(
        _s6_arm_val is not None
        or observations.get('s6_dose_matched_1500_ieq_verified') is not None
        or any(_mk in validated_stage0 for _mk in _s6_mean_keys)
        or (_divider_qc_active and _dz_arm_val is None)
    )
    _p05_active = bool(
        _p05_geom_val is not None
        or observations.get('p05_full_336h_coupled_o2_h2o2_ca_ph_qualified') is not None
        or any(_mk in validated_stage0 for _mk in _p05_mean_keys)
    )
    _dz_active = bool(
        (product_def is None or _is_design_z_product)
        and (
            _dz_arm_val is not None
            or observations.get('design_z_empty_lumen_14d_persistence_verified') is not None
            or observations.get('design_z_loaded_closed_state_confirmed') is not None
            or any(_mk in validated_stage0 for _mk in (*_dz_empty_mean_keys, *_dz_loaded_mean_keys))
        )
    )

    iip0019_s6_status: str | None = None
    iip0019_p05_status: str | None = None
    iip0019_dz_status: str | None = None
    iip0019_transport_hold_status: str | None = None

    if _s6_active:
        _s6_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _s6_mean_keys
        )
        if _s6_arm_val is None or not _s6_quant_ok:
            iip0019_s6_status = 'HOLD_S6_INCOMPLETE_OR_DEFAULT_SE_LATERAL_SPACING_PANEL'
        elif _s6_arm_val == 'DERATED_80UM_768IEQ' or observations.get('s6_dose_matched_1500_ieq_verified') is not True:
            iip0019_s6_status = 'HOLD_S6_DOSE_DERATED_80UM_768IEQ_NOT_SINGLE_POUCH_EQUIVALENT_REQUIRES_DOSE_MATCHED_1500IEQ_OR_AREA_EXPANSION'
        elif _divider_qc_active and not _divider_qc_cleared:
            iip0019_s6_status = 'HOLD_S6_EXPLORATORY_DIVIDER_HYDRATED_PO2_LEACHABLES_OR_MODULUS_QC_FAILED'
        else:
            _diam_cap = 80.0 if _s6_arm_val == 'DOSE_MATCHED_80UM_1500IEQ' else 100.0
            _gap_floor = 55.0 if _s6_arm_val == 'DOSE_MATCHED_80UM_1500IEQ' else 85.0
            if (
                _upper95('obs_s6_cluster_diam_p90_um') > _diam_cap
                or _lower95('obs_s6_p10_nearest_neighbor_edge_gap_um') < _gap_floor
            ):
                iip0019_s6_status = 'HOLD_S6_LATERAL_CLUMPING_EDGE_GAP_OR_CLUSTER_DIAMETER_CAP_FAILED'
            else:
                iip0019_s6_status = 'CLEARED_S6_LATERAL_SPACING_AND_DOSE_MATCHED_CLUSTER_CAP'
        if iip0019_s6_status != 'CLEARED_S6_LATERAL_SPACING_AND_DOSE_MATCHED_CLUSTER_CAP':
            iip0019_transport_hold_status = iip0019_s6_status

    if _p05_active:
        _p05_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _p05_mean_keys
        )
        if _p05_geom_val is None or not _p05_quant_ok:
            iip0019_p05_status = 'HOLD_P05_INCOMPLETE_OR_DEFAULT_SE_BARRIER_QUALIFICATION_PANEL'
        elif observations.get('p05_full_336h_coupled_o2_h2o2_ca_ph_qualified') is not True:
            iip0019_p05_status = 'HOLD_P05_FULL_336H_COUPLED_TRANSPORT_QUALIFICATION_UNVERIFIED'
        elif _upper95('obs_p05_islet_plane_h2o2_max_uM') > 5.0:
            iip0019_p05_status = 'HOLD_P05_ISLET_PLANE_H2O2_EXCEEDS_5P0UM_TIER2_ABORT_CEILING'
        elif _upper95('obs_p05_islet_plane_h2o2_max_uM') > 2.5:
            iip0019_p05_status = 'HOLD_P05_ISLET_PLANE_H2O2_EXCEEDS_2P5UM_TIER1_TARGET_REQUIRES_FULL_73800IU_RETENTION'
        elif (
            _upper95('obs_p05_microbead_volume_d99_um') > 50.0
            or _lower95('obs_p05_catalase_barrier_film_um') < 20.0
            or _lower95('obs_p05_day7_catalase_retention_pct') < 70.0
            or _lower95('obs_p05_336h_retained_catalase_iu_333cm2') < 73800.0
        ):
            iip0019_p05_status = 'HOLD_P05_MICROBEAD_D99_FILM_THICKNESS_OR_CATALASE_RETENTION_FAILED'
        else:
            iip0019_p05_status = 'CLEARED_P05_WET_SIEVED_MICROBEAD_AND_20UM_CATALASE_BARRIER_GATE'
        if iip0019_p05_status != 'CLEARED_P05_WET_SIEVED_MICROBEAD_AND_20UM_CATALASE_BARRIER_GATE' and iip0019_transport_hold_status is None:
            iip0019_transport_hold_status = iip0019_p05_status

    if _dz_active:
        _dz_empty_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _dz_empty_mean_keys
        )
        _dz_loaded_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _dz_loaded_mean_keys
        )
        if _dz_arm_val is None or not _dz_empty_ok:
            iip0019_dz_status = 'HOLD_DESIGN_Z_INCOMPLETE_OR_DEFAULT_SE_EMPTY_LUMEN_OPTODE_PANEL'
        elif _dz_arm_val == 'EXPLORATORY_18PCT_DIVIDER_140UM' and not _divider_qc_cleared:
            iip0019_dz_status = 'HOLD_S6_EXPLORATORY_DIVIDER_HYDRATED_PO2_LEACHABLES_OR_MODULUS_QC_FAILED'
        else:
            _lumen_floor = 35.0 if _dz_arm_val == 'EXPLORATORY_18PCT_DIVIDER_140UM' else 40.0
            _standoff_cap = 20.0 if _dz_arm_val == 'EXPLORATORY_18PCT_DIVIDER_140UM' else 15.0
            if (
                observations.get('design_z_empty_lumen_14d_persistence_verified') is not True
                or _lower95('obs_design_z_empty_lumen_5site_min_po2_mmhg') < _lumen_floor
                or _upper95('obs_design_z_vessel_standoff_p90_um') > _standoff_cap
            ):
                iip0019_dz_status = 'HOLD_DESIGN_Z_EMPTY_LUMEN_PO2_STANDOFF_OR_14D_PERSISTENCE_FAILED'
            elif observations.get('design_z_loaded_closed_state_confirmed') is not True or not _dz_loaded_ok:
                iip0019_dz_status = 'HOLD_DESIGN_Z_EMPTY_LUMEN_ONLY_REQUIRES_LOADED_CLOSED_STATE_CONFIRMATION'
            elif (
                _lower95('obs_design_z_loaded_surface_po2_mmhg') < 15.0
                or _lower95('obs_design_z_loaded_core_po2_mmhg') < 5.0
                or _lower95('obs_design_z_loaded_perifusion_gsis_si') < 3.0
            ):
                iip0019_dz_status = 'HOLD_DESIGN_Z_LOADED_CLOSED_STATE_PO2_OR_READOUT2_GSIS_SI3_FAILED'
            else:
                iip0019_dz_status = 'CLEARED_DESIGN_Z_TWO_STAGE_EMPTY_LUMEN_AND_LOADED_CLOSED_STATE_GATE'
        if iip0019_dz_status != 'CLEARED_DESIGN_Z_TWO_STAGE_EMPTY_LUMEN_AND_LOADED_CLOSED_STATE_GATE' and iip0019_transport_hold_status is None:
            iip0019_transport_hold_status = iip0019_dz_status

    iip0019_eval_summary: dict[str, object] | None = None
    if _s6_active or _p05_active or _dz_active:
        iip0019_eval_summary = {
            'gate_pass': iip0019_transport_hold_status is None,
            'hold_status': iip0019_transport_hold_status,
            's6_cluster_size_arm': _s6_arm_val,
            's6_lateral_spacing_status': iip0019_s6_status,
            'p05_barrier_slab_geometry_mode': _p05_geom_val,
            'p05_microbead_barrier_status': iip0019_p05_status,
            'design_z_divider_arm': _dz_arm_val,
            'design_z_two_stage_status': iip0019_dz_status,
        }

    # Stage-0 Gate P05 recombinant human CAT (P04040) species/shed-mass/endotoxin lock and anti-PEG / Day-28 anti-SA humoral + DPR-maleimide/hepatocyte cuff gate (IIP-0020; DEC-099):
    _req_iip0020_full = observations.get('require_full_iip0020_biomaterial_and_humoral_panel') is True
    _p05_cat_acc = observations.get('p05_catalase_uniprot_accession')
    if isinstance(_p05_cat_acc, str):
        _p05_cat_acc = _p05_cat_acc.strip()
    _cuff_arm = observations.get('cuff_bioconjugation_arm')
    if isinstance(_cuff_arm, str):
        _cuff_arm = _cuff_arm.strip()
    _cuff_stratum = observations.get('cuff_serum_challenge_stratum')
    if isinstance(_cuff_stratum, str):
        _cuff_stratum = _cuff_stratum.strip()

    _iip0020_p05_mean_keys = (
        'obs_p05_catalase_endotoxin_eu_per_mg',
        'obs_p05_total_implant_catalase_endotoxin_eu',
        'obs_p05_day7_shed_catalase_mass_pct',
    )
    _iip0020_cuff_core_mean_keys = (
        'obs_cuff_surface_c3b_deposition_pct_pos_ctrl',
        'obs_cuff_soluble_sc5b9_pct_pos_ctrl',
        'obs_cuff_macrophage_fcgammar_cytokine_fold',
        'obs_cuff_primed_serum_thbd_apc_retention_pct',
        'obs_cuff_primed_serum_cfh_c3b_cleavage_retention_pct',
        'obs_cuff_dpr_ring_hydrolysis_pct',
        'obs_cuff_hepatocyte_72h_caspase37_fold',
        'obs_cuff_hepatocyte_72h_alt_ldh_fold',
    )
    _iip0020_cuff_donor_count_keys = (
        'obs_cuff_serum_donor_count',
        'obs_cuff_serum_failing_donor_count',
    )
    _iip0020_fasl_mean_keys = (
        'obs_cuff_primed_serum_fasl_ec50_fold_shift',
        'obs_cuff_primed_serum_fasl_emax_apoptosis_pct',
        'obs_cuff_14d_hsa_shed_bioactive_fasl_pct',
    )
    _iip0020_p05_species_active = bool(
        _req_iip0020_full
        or (_req_iip0025_full and not _is_design_z_product)
        or (product_def is not None and not _is_design_z_product and _p05_active)
        or _p05_cat_acc is not None
        or observations.get('p05_recombinant_human_catalase_verified') is not None
        or any(_mk in validated_stage0 for _mk in _iip0020_p05_mean_keys)
    )
    _iip0020_cuff_active = bool(
        _req_iip0020_full
        or _cuff_arm is not None
        or _cuff_stratum is not None
        or observations.get('cuff_dpr_ring_hydrolysis_lcms_verified') is not None
        or observations.get('cuff_conjugate_disulfide_and_site_occupancy_lcms_verified') is not None
        or observations.get('cuff_hexafasl_inventory_within_provisional_cap_verified') is not None
        or any(_mk in validated_stage0 for _mk in (*_iip0020_cuff_core_mean_keys, *_iip0020_cuff_donor_count_keys, *_iip0020_fasl_mean_keys, 'obs_cuff_h1a_sec_mals_hexamer_purity_pct', 'obs_cuff_unreacted_free_thiol_pct'))
    )
    iip0020_p05_species_status: str | None = None
    iip0020_cuff_status: str | None = None
    iip0020_hold_status: str | None = None

    # Product-bind Gate P05 strictly to CaO2-containing routes (not _is_design_z_product) so extra Design Z fields never suppress species/barrier holds (IIP-0025; DEC-104):
    if _iip0020_p05_species_active and not _is_design_z_product:
        _iip20_p05_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _iip0020_p05_mean_keys
        )
        if _p05_cat_acc != 'P04040' or observations.get('p05_recombinant_human_catalase_verified') is not True:
            iip0020_p05_species_status = 'HOLD_P05_XENO_OR_UNVERIFIED_CATALASE_SPECIES_REQUIRES_HUMAN_P04040'
        elif not _iip20_p05_quant_ok:
            iip0020_p05_species_status = 'HOLD_P05_INCOMPLETE_OR_DEFAULT_SE_CATALASE_ENDOTOXIN_OR_SHED_MASS_PANEL'
        elif (
            _upper95('obs_p05_catalase_endotoxin_eu_per_mg') > 0.10
            or _upper95('obs_p05_total_implant_catalase_endotoxin_eu') > 0.50
            or _upper95('obs_p05_day7_shed_catalase_mass_pct') > 30.0
        ):
            iip0020_p05_species_status = 'HOLD_P05_CATALASE_ENDOTOXIN_OR_DAY7_SHED_MASS_EXCEEDED'
        else:
            iip0020_p05_species_status = 'CLEARED_P05_HUMAN_P04040_SPECIES_ENDOTOXIN_AND_SHED_MASS_GATE'
        if iip0020_p05_species_status != 'CLEARED_P05_HUMAN_P04040_SPECIES_ENDOTOXIN_AND_SHED_MASS_GATE':
            iip0020_hold_status = iip0020_p05_species_status

    if _iip0020_cuff_active:
        _is_fasl_free = _cuff_arm == 'DPR_HYDROLYZED_CFH_THBD_FASL_FREE'
        _is_h1a = _cuff_arm == 'DPR_HYDROLYZED_CFH_THBD_H1A_ADIPOQ_HEXAFASL'
        _is_dpr_covalent = _cuff_arm in (
            'DPR_HYDROLYZED_CFH_THBD_H1A_ADIPOQ_HEXAFASL',
            'DPR_HYDROLYZED_CFH_THBD_T0_TNC_FASL',
            'DPR_HYDROLYZED_CFH_THBD_FASL_FREE',
        )
        _req_cuff_keys = list(_iip0020_cuff_core_mean_keys)
        if not _is_fasl_free:
            _req_cuff_keys.extend(_iip0020_fasl_mean_keys)
        if _is_dpr_covalent:
            _req_cuff_keys.append('obs_cuff_unreacted_free_thiol_pct')
        if _is_h1a:
            _req_cuff_keys.append('obs_cuff_h1a_sec_mals_hexamer_purity_pct')
        _cuff_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _req_cuff_keys
        ) and all(
            _mk in validated_stage0 and abs(validated_stage0[_mk] - round(validated_stage0[_mk])) <= 1e-9
            for _mk in _iip0020_cuff_donor_count_keys
        )
        _cuff_donor_n = int(round(validated_stage0['obs_cuff_serum_donor_count'])) if 'obs_cuff_serum_donor_count' in validated_stage0 else -1
        _cuff_fail_n = int(round(validated_stage0['obs_cuff_serum_failing_donor_count'])) if 'obs_cuff_serum_failing_donor_count' in validated_stage0 else 999
        if _cuff_arm is None or _cuff_stratum is None or not _cuff_quant_ok:
            iip0020_cuff_status = 'HOLD_CUFF_INCOMPLETE_OR_DEFAULT_SE_HUMORAL_MALEIMIDE_HEPATOCYTE_PANEL'
        elif observations.get('cuff_dpr_ring_hydrolysis_lcms_verified') is not True or _lower95('obs_cuff_dpr_ring_hydrolysis_pct') < 98.0:
            iip0020_cuff_status = 'HOLD_CUFF_THIOSUCCINIMIDE_RING_HYDROLYSIS_LT_98PCT_RETRO_MICHAEL_HAZARD'
        elif _is_dpr_covalent and (
            observations.get('cuff_conjugate_disulfide_and_site_occupancy_lcms_verified') is not True
            or _upper95('obs_cuff_unreacted_free_thiol_pct') > 5.0
        ):
            iip0020_cuff_status = 'HOLD_CUFF_CONJUGATE_DISULFIDE_SITE_OCCUPANCY_OR_FREE_THIOL_FAILED'
        elif _is_h1a and _lower95('obs_cuff_h1a_sec_mals_hexamer_purity_pct') < 90.0:
            iip0020_cuff_status = 'HOLD_CUFF_H1A_HEXAMER_PURITY_LT_90PCT'
        elif not _is_fasl_free and (
            observations.get('cuff_hexafasl_inventory_within_provisional_cap_verified') is not True
            or _upper95('obs_cuff_14d_hsa_shed_bioactive_fasl_pct') > 5.0
            or _upper95('obs_cuff_primed_serum_fasl_ec50_fold_shift') > 2.0
            or _lower95('obs_cuff_primed_serum_fasl_ec50_fold_shift') < 0.50
            or _lower95('obs_cuff_primed_serum_fasl_emax_apoptosis_pct') < 80.0
        ):
            iip0020_cuff_status = 'HOLD_CUFF_FASL_CAP_SHEDDING_OR_PRIMED_SERUM_POTENCY_FAILED'
        elif (
            _cuff_donor_n < 12
            or _cuff_fail_n > 0
            or _upper95('obs_cuff_surface_c3b_deposition_pct_pos_ctrl') > 5.0
            or _upper95('obs_cuff_soluble_sc5b9_pct_pos_ctrl') > 5.0
            or _upper95('obs_cuff_macrophage_fcgammar_cytokine_fold') > 1.5
            or _lower95('obs_cuff_primed_serum_thbd_apc_retention_pct') < 80.0
            or _lower95('obs_cuff_primed_serum_cfh_c3b_cleavage_retention_pct') < 80.0
        ):
            iip0020_cuff_status = 'HOLD_CUFF_ANTI_PEG_OR_DAY28_ANTI_SA_HUMORAL_COMPLEMENT_OR_COFACTOR_FAILED'
        elif (
            _upper95('obs_cuff_hepatocyte_72h_caspase37_fold') > 1.20
            or _upper95('obs_cuff_hepatocyte_72h_alt_ldh_fold') > 1.20
        ):
            iip0020_cuff_status = 'HOLD_CUFF_PRIMARY_HUMAN_HEPATOCYTE_LEACHATE_OR_DIGEST_TOXICITY_EXCEEDED'
        else:
            iip0020_cuff_status = 'CLEARED_CUFF_HUMORAL_HYDROLYSIS_AND_HEPATOCYTE_SAFETY_GATE'
        if iip0020_cuff_status != 'CLEARED_CUFF_HUMORAL_HYDROLYSIS_AND_HEPATOCYTE_SAFETY_GATE' and iip0020_hold_status is None:
            iip0020_hold_status = iip0020_cuff_status

    iip0020_eval_summary: dict[str, object] | None = None
    if _iip0020_p05_species_active or _iip0020_cuff_active:
        iip0020_eval_summary = {
            'gate_pass': iip0020_hold_status is None,
            'hold_status': iip0020_hold_status,
            'p05_catalase_uniprot_accession': _p05_cat_acc,
            'p05_catalase_species_status': iip0020_p05_species_status,
            'cuff_bioconjugation_arm': _cuff_arm,
            'cuff_serum_challenge_stratum': _cuff_stratum,
            'cuff_humoral_and_safety_status': iip0020_cuff_status,
        }

    # Stage-0 120-h lineage-resolved metabolic-rest OCR coupon and 4-arm Nec-1s x Ferrostatin-1 pellet-LOQ/NQO1 conditioning gate (IIP-0021; DEC-100):
    _req_iip0021_full = observations.get('require_full_iip0021_ocr_and_conditioning_panel') is True
    _rest_arm = observations.get('rest_ocr_coupon_arm')
    if isinstance(_rest_arm, str):
        _rest_arm = _rest_arm.strip()
    _cond_arm = observations.get('conditioning_4arm_arm')
    if isinstance(_cond_arm, str):
        _cond_arm = _cond_arm.strip()
    _cond_scope = observations.get('conditioning_washout_analyte_scope')
    if isinstance(_cond_scope, str):
        _cond_scope = _cond_scope.strip()

    _iip0021_rest_keys = (
        'obs_rest_120h_graft_ocr_ratio',
        'obs_rest_depot_release_half_life_h',
        'obs_rest_post_washout_gsis_si',
        'obs_rest_post_washout_gsis_output_retention_pct',
        'obs_rest_post_washout_alpha_glucagon_fold',
    )
    _iip0021_cond_keys = (
        'obs_cond_sc_beta_viable_yield_pct',
        'obs_cond_iec_viable_yield_pct',
        'obs_cond_final_wash_supernatant_nM',
        'obs_cond_pellet_residual_over_loq_ratio',
        'obs_cond_pellet_extraction_recovery_pct',
        'obs_cond_endocrine_nqo1_activity_pct',
        'obs_cond_vascular_nqo1_activity_pct',
        'obs_cond_post_washout_gsis_si',
    )
    _iip0021_rest_active = bool(
        _req_iip0021_full
        or _rest_arm is not None
        or observations.get('rest_120h_longitudinal_simultaneous_coverage_verified') is not None
        or observations.get('rest_120h_viable_lineage_yield_and_during_exposure_function_verified') is not None
        or any(_mk in validated_stage0 for _mk in _iip0021_rest_keys)
    )
    _iip0021_cond_active = bool(
        _req_iip0021_full
        or _cond_arm is not None
        or _cond_scope is not None
        or observations.get('conditioning_pellet_loq_and_extraction_recovery_verified') is not None
        or observations.get('conditioning_orthogonal_pathway_controls_verified') is not None
        or any(_mk in validated_stage0 for _mk in _iip0021_cond_keys)
    )
    iip0021_rest_status: str | None = None
    iip0021_cond_status: str | None = None
    iip0021_hold_status: str | None = None

    if _iip0021_rest_active:
        if _rest_arm == 'VEHICLE_791P5MG_CAO2_BASELINE':
            iip0021_rest_status = 'BASELINE_791P5MG_CAO2_RETAINED_NO_METABOLIC_REST_DERATING'
        elif _rest_arm == 'UNENCAPSULATED_FREE_DRUG_REJECTED':
            iip0021_rest_status = 'HOLD_REST_OCR_UNENCAPSULATED_FREE_DRUG_RAPID_SLAB_WASHOUT_BARRED'
        else:
            _rest_quant_ok = all(
                _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
                for _mk in _iip0021_rest_keys
            )
            if _rest_arm != 'ENCAPSULATED_KATP_MICRODEPOT_CHALLENGER' or not _rest_quant_ok:
                iip0021_rest_status = 'HOLD_REST_OCR_INCOMPLETE_OR_DEFAULT_SE_120H_COUPON_PANEL'
            elif (
                observations.get('rest_120h_longitudinal_simultaneous_coverage_verified') is not True
                or observations.get('rest_120h_viable_lineage_yield_and_during_exposure_function_verified') is not True
            ):
                iip0021_rest_status = 'HOLD_REST_OCR_MISSING_LONGITUDINAL_COVERAGE_OR_VIABLE_LINEAGE_FUNCTION_CHECK'
            elif _lower95('obs_rest_depot_release_half_life_h') < 72.0 or _upper95('obs_rest_120h_graft_ocr_ratio') > 0.80:
                iip0021_rest_status = 'HOLD_REST_OCR_120H_RATIO_GT_0P80_OR_DEPOT_HALFLIFE_LT_72H_KEEP_791P5MG_CAO2'
            elif (
                _lower95('obs_rest_post_washout_gsis_si') < 3.0
                or _lower95('obs_rest_post_washout_gsis_output_retention_pct') < 85.0
                or _lower95('obs_rest_post_washout_alpha_glucagon_fold') < 2.0
            ):
                iip0021_rest_status = 'HOLD_REST_OCR_POST_WASHOUT_GSIS_ABSOLUTE_OUTPUT_OR_GLUCAGON_FAILED'
            else:
                iip0021_rest_status = 'CLEARED_REST_120H_OCR_COUPON_ELIGIBLE_FOR_STAGE2_COUPLED_CAO2_QUALIFICATION_ONLY'
        if iip0021_rest_status not in (
            'BASELINE_791P5MG_CAO2_RETAINED_NO_METABOLIC_REST_DERATING',
            'CLEARED_REST_120H_OCR_COUPON_ELIGIBLE_FOR_STAGE2_COUPLED_CAO2_QUALIFICATION_ONLY',
        ):
            iip0021_hold_status = iip0021_rest_status

    if _iip0021_cond_active:
        if _cond_arm in ('NEC1S_100UM_BENCHMARK_COMPARATOR', 'NEC1_100UM_ANTIFERROPTOTIC_COMPARATOR'):
            iip0021_cond_status = 'HOLD_CONDITIONING_BENCHMARK_OR_NEC1_ANTIFERROPTOTIC_COMPARATOR_NOT_RELEASE_ARM'
        elif _cond_scope == 'AQUEOUS_SUPERNATANT_ONLY_INSUFFICIENT':
            iip0021_cond_status = 'HOLD_CONDITIONING_AQUEOUS_SUPERNATANT_ONLY_REQUIRES_WASHED_CELL_PELLET_LCMS'
        else:
            _cond_quant_ok = all(
                _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
                for _mk in _iip0021_cond_keys
            )
            if _cond_arm is None or _cond_scope != 'SUPERNATANT_AND_CELL_PELLET_PARENT_PLUS_METABOLITES' or not _cond_quant_ok:
                iip0021_cond_status = 'HOLD_CONDITIONING_INCOMPLETE_OR_DEFAULT_SE_4ARM_PELLET_NQO1_PANEL'
            elif (
                observations.get('conditioning_pellet_loq_and_extraction_recovery_verified') is not True
                or observations.get('conditioning_orthogonal_pathway_controls_verified') is not True
                or _lower95('obs_cond_pellet_extraction_recovery_pct') < 80.0
            ):
                iip0021_cond_status = 'HOLD_CONDITIONING_PELLET_LOQ_EXTRACTION_RECOVERY_OR_PATHWAY_CONTROLS_UNVERIFIED'
            elif (
                _upper95('obs_cond_final_wash_supernatant_nM') >= 10.0
                or _upper95('obs_cond_pellet_residual_over_loq_ratio') > 1.0
                or _lower95('obs_cond_endocrine_nqo1_activity_pct') < 90.0
                or _lower95('obs_cond_vascular_nqo1_activity_pct') < 90.0
            ):
                iip0021_cond_status = 'HOLD_CONDITIONING_RESIDUAL_DRUG_OR_NQO1_INHIBITION_DETECTED'
            elif (
                _lower95('obs_cond_sc_beta_viable_yield_pct') < 85.0
                or _lower95('obs_cond_iec_viable_yield_pct') < 85.0
                or _lower95('obs_cond_post_washout_gsis_si') < 3.0
            ):
                iip0021_cond_status = 'HOLD_CONDITIONING_LINEAGE_VIABLE_YIELD_OR_POST_WASHOUT_GSIS_SI3_FAILED'
            else:
                iip0021_cond_status = 'CLEARED_CONDITIONING_4ARM_PELLET_LOQ_NQO1_AND_LINEAGE_GSIS_GATE'
        if iip0021_cond_status != 'CLEARED_CONDITIONING_4ARM_PELLET_LOQ_NQO1_AND_LINEAGE_GSIS_GATE' and iip0021_hold_status is None:
            iip0021_hold_status = iip0021_cond_status

    iip0021_eval_summary: dict[str, object] | None = None
    if _iip0021_rest_active or _iip0021_cond_active:
        iip0021_eval_summary = {
            'gate_pass': iip0021_hold_status is None,
            'hold_status': iip0021_hold_status,
            'rest_ocr_coupon_arm': _rest_arm,
            'rest_120h_ocr_status': iip0021_rest_status,
            'conditioning_4arm_arm': _cond_arm,
            'conditioning_washout_analyte_scope': _cond_scope,
            'conditioning_pellet_nqo1_status': iip0021_cond_status,
        }

    # Stage-0 Track-2 B'-SCT vs HLA-G1 split, B2M-null PVR/CD226 NK check, and lineage-separated necroptosis guardrail evaluation (IIP-0022; DEC-101):
    _req_iip0022_full = observations.get('require_full_iip0022_track2_and_guardrail_panel') is True
    _iip0022_active = bool(
        _req_iip0022_full
        or _t2_route is not None
        or _bsct_flag
        or _hlag1_flag
        or _bsct_quant_present
        or _hlag1_quant_present
        or 'obs_b2m_null_pvr_cd226_nk_lysis_pct' in validated_stage0
        or observations.get('lineage_necroptosis_time_resolved_viability_and_rescue_verified') is not None
    )
    iip0022_hold_status: str | None = None
    iip0022_track2_status: str | None = None
    if _iip0022_active:
        if _bsct_flag and _hlag1_flag:
            iip0022_track2_status = 'HOLD_TRACK2_MUTUALLY_EXCLUSIVE_BSCT_AND_HLAG1_FLAGS_BOTH_SET'
            iip0022_hold_status = iip0022_track2_status
        elif _legacy_conflated_t2 and (_t2_route is not None or _req_iip0022_full):
            iip0022_track2_status = 'HOLD_TRACK2_DEPRECATED_CONFLATED_RESCUE_FLAG_REQUIRES_EXPLICIT_SPLIT_ROUTE'
            iip0022_hold_status = iip0022_track2_status
        elif track2_hlag1_route_attempted:
            if not _hlag1_quant_cleared_uniform_8pct:
                iip0022_track2_status = 'HOLD_TRACK2_HLAG1_LILRB1_NEG_KIR_NK_EXCEEDS_8PCT_OR_UNDERPOWERED_DONOR_PANEL'
            else:
                # Even if sorted LILRB1- NKG2A- KIR+ primary NK clears <= 8.0% across n >= 6 donors, HLA-G1 remains HELD until registered in islet1_product_definitions_v0.8.0.json:
                iip0022_track2_status = 'HOLD_TRACK2_HLAG1_CLEARED_UNIFORM_8PCT_COUPON_PENDING_PRODUCT_DEFINITION_REGISTRATION'
            iip0022_hold_status = iip0022_track2_status
        elif _bsct_flag or _t2_route == 'TRACK2_HLAC_B_SCT_KIR_RESCUE' or _bsct_quant_present:
            if not _bsct_quant_cleared:
                iip0022_track2_status = 'HOLD_TRACK2_BSCT_KIR_NK_EXCEEDS_8PCT_OR_MISSING_EXPLICIT_6DONOR_PANEL'
                iip0022_hold_status = iip0022_track2_status
            else:
                iip0022_track2_status = 'CLEARED_TRACK2_BSCT_KIR_6DONOR_8PCT_COUPON'
        else:
            iip0022_track2_status = 'TRACK1_MHC_NULL_ACTIVE_NO_TRACK2_SWITCH'

        if iip0022_hold_status is None and ('obs_b2m_null_pvr_cd226_nk_lysis_pct' in validated_stage0 or _req_iip0022_full):
            if (
                'obs_b2m_null_pvr_cd226_nk_lysis_pct' not in validated_stage0
                or se_provenance_by_mean.get('obs_b2m_null_pvr_cd226_nk_lysis_pct') != 'EXPLICIT_SE_SUPPLIED'
                or _upper95('obs_b2m_null_pvr_cd226_nk_lysis_pct') > 8.0
            ):
                iip0022_hold_status = 'HOLD_B2M_NULL_TIGIT_LOW_CD226_HIGH_NK_LYSIS_EXCEEDS_8PCT_OR_MISSING_EXPLICIT_SE'
        if iip0022_hold_status is None and (observations.get('lineage_necroptosis_time_resolved_viability_and_rescue_verified') is False or (_req_iip0022_full and observations.get('lineage_necroptosis_time_resolved_viability_and_rescue_verified') is not True)):
            iip0022_hold_status = 'HOLD_LINEAGE_NECROPTOSIS_TIME_RESOLVED_VIABILITY_AND_RESCUE_UNVERIFIED'

    iip0022_eval_summary: dict[str, object] | None = None
    if _iip0022_active:
        iip0022_eval_summary = {
            'gate_pass': iip0022_hold_status is None,
            'hold_status': iip0022_hold_status,
            'track2_rescue_route': _t2_route,
            'track2_route_status': iip0022_track2_status,
            'deprecated_conflated_flag_supplied': _legacy_conflated_t2,
        }

    # Stage-0 Design Y wildcard fail-closed Y-01..Y-04 evaluation (IIP-0023; DEC-102):
    # Completely inactive when no design_y_* / obs_design_y* key is supplied so A/B/X behavior is 100% preserved.
    _dy_int_mod = observations.get('design_y_integrase_modality')
    if isinstance(_dy_int_mod, str):
        _dy_int_mod = _dy_int_mod.strip()
    _dy_sil_mod = observations.get('design_y_silencing_modality')
    if isinstance(_dy_sil_mod, str):
        _dy_sil_mod = _dy_sil_mod.strip()
    _dy_pay_arm = observations.get('design_y_payload_arm')
    if isinstance(_dy_pay_arm, str):
        _dy_pay_arm = _dy_pay_arm.strip()
    _dy_bool_keys = (
        'design_y_attlr_outside_orf_5p3p_race_verified',
        'design_y_zero_2a_scar_proteomics_verified',
        'design_y_dpcpx_adora1_antagonist_control_run',
        'design_y_anti_cd200_blockade_reversal_verified',
        'design_y_y04_full_coculture_and_clonogenic_confirmed',
    )
    _dy_quant_keys = (
        'obs_design_y01_max_lineage_escape_pct',
        'obs_design_y01_lineage_cells_scored',
        'obs_design_y01_cbe_offtarget_pct',
        'obs_design_y02_pre_enrichment_ki_pct',
        'obs_design_y02_clonal_ddpcr_cn',
        'obs_design_y02_wt_b2m_residual_allele_frac',
        'obs_design_y02_backbone_kanr_vaf',
        'obs_design_y02_wgs_pseudosite_k',
        'obs_design_y03_iec_beta_selectivity_ratio',
        'obs_design_y03_beta_mir375_retention_pct',
        'obs_design_y03_coculture_gsis_si',
        'obs_design_y03_gsis_retention_pct',
        'obs_design_y03_supernatant_adenosine_nM',
        'obs_design_y03_dpcpx_gsis_delta_pct',
        'obs_design_y04_multicistron_retention_pct',
        'obs_design_y04_nhs20_c5b9_mac_pct',
        'obs_design_y04_indirect_apc_cd4_si',
        'obs_design_y04_ap1903_clonogenic_survival_frac',
        'obs_design_y_sampling_age_days',
        'obs_design_y02_ont_intact_homolog_pct',
        'obs_design_y04_all4_lineages_scored_count',
        'obs_design_y04_ap1903_survivor_count_k',
        'obs_design_y04_ap1903_effective_trials',
    )
    _iip0023_dy_active = bool(
        _dy_int_mod is not None
        or _dy_sil_mod is not None
        or _dy_pay_arm is not None
        or any(observations.get(_bk) is not None for _bk in _dy_bool_keys)
        or any(_qk in validated_stage0 for _qk in _dy_quant_keys)
    )
    iip0023_dy_status = 'NOT_TESTED'
    iip0023_hold_status: str | None = None
    if _iip0023_dy_active:
        if _dy_int_mod == 'LEGACY_ATTB_GENOME_ATTP_PLASMID_REJECTED':
            iip0023_dy_status = 'HOLD_DESIGN_Y02_LEGACY_ATTB_IN_GENOME_POLARITY_OT4_PSEUDOSITE_RISK_REJECTED'
        elif _dy_pay_arm == 'DESIGN_Y_CD39_IEC_COMPARATOR_12555BP':
            iip0023_dy_status = 'HOLD_DESIGN_Y03_CD39_ECTONUCLEOTIDASE_COMPARATOR_NOT_PRIMARY_PAYLOAD'
        elif (
            _dy_int_mod != 'EEPASSIGE_B2M_ATTP_GENOME_ATTB_MINICIRCLE'
            or _dy_sil_mod is None
            or _dy_pay_arm is None
            or not all(observations.get(_bk) is True for _bk in _dy_bool_keys)
            or not all(_qk in validated_stage0 and se_provenance_by_mean.get(_qk) == 'EXPLICIT_SE_SUPPLIED' for _qk in _dy_quant_keys)
        ):
            iip0023_dy_status = 'HOLD_DESIGN_Y_INCOMPLETE_OR_DEFAULT_SE_Y01_TO_Y04_PANEL'
        elif (
            validated_stage0['obs_design_y_sampling_age_days'] < 60.0
            or abs(validated_stage0['obs_design_y04_all4_lineages_scored_count'] - 4.0) > 1e-9
            or _lower95('obs_design_y01_lineage_cells_scored') < 10000.0
            or _upper95('obs_design_y01_max_lineage_escape_pct') > 1.0
            or _upper95('obs_design_y01_cbe_offtarget_pct') > 0.10
        ):
            iip0023_dy_status = 'HOLD_DESIGN_Y01_60D_4LINEAGE_SILENCE_OR_TADCBED_HAPLOTYPE_FAILED'
        elif (
            _lower95('obs_design_y02_pre_enrichment_ki_pct') < 3.0
            or _lower95('obs_design_y02_clonal_ddpcr_cn') < 1.90
            or _upper95('obs_design_y02_clonal_ddpcr_cn') > 2.10
            or _lower95('obs_design_y02_ont_intact_homolog_pct') < 99.0
            or _upper95('obs_design_y02_wt_b2m_residual_allele_frac') > 0.001
            or _upper95('obs_design_y02_backbone_kanr_vaf') > 0.0001
            or abs(validated_stage0['obs_design_y02_wgs_pseudosite_k']) > 1e-9
        ):
            iip0023_dy_status = 'HOLD_DESIGN_Y02_EEPASSIGE_HOMOZYGOUS_CN2_KANR_OR_PSEUDOSITE_GATE_FAILED'
        elif (
            _lower95('obs_design_y03_iec_beta_selectivity_ratio') < 50.0
            or _lower95('obs_design_y03_beta_mir375_retention_pct') < 90.0
            or _lower95('obs_design_y03_coculture_gsis_si') < 2.5
            or _lower95('obs_design_y03_gsis_retention_pct') < 90.0
            or _upper95('obs_design_y03_supernatant_adenosine_nM') > 25.0
            or _upper95('obs_design_y03_dpcpx_gsis_delta_pct') > 10.0
        ):
            iip0023_dy_status = 'HOLD_DESIGN_Y03_ADORA1_PARACRINE_GSIS_OR_CD200_SELECTIVITY_FAILED'
        else:
            _dy_surv_k = validated_stage0['obs_design_y04_ap1903_survivor_count_k']
            _dy_neff = validated_stage0['obs_design_y04_ap1903_effective_trials']
            _dy_poisson_u95 = _poisson_exact_upper95_rate(_dy_surv_k, _dy_neff) if _dy_neff > 0.0 else 1.0
            if (
                _lower95('obs_design_y04_multicistron_retention_pct') < 95.0
                or _upper95('obs_design_y04_nhs20_c5b9_mac_pct') > 5.0
                or _upper95('obs_design_y04_indirect_apc_cd4_si') > 1.5
                or 'obs_kir_nkg2a_neg_nk_lysis_pct' not in validated_stage0
                or kir_nkg2a_neg_upper95 > 8.0
                or _dy_neff < 299573228.0
                or abs(_dy_surv_k) > 1e-9
                or _dy_poisson_u95 > 1e-8
                or _upper95('obs_design_y04_ap1903_clonogenic_survival_frac') > 1e-8
                or ('obs_gsis_si' in validated_stage0 and _lower95('obs_gsis_si') < 3.0)
            ):
                iip0023_dy_status = 'HOLD_DESIGN_Y04_MULTICISTRON_COMPLEMENT_KIR_OR_3E8_CLONOGENIC_GATE_FAILED'
            else:
                iip0023_dy_status = 'CLEARED_DESIGN_Y_WILDCARD_COUPONS_ELIGIBLE_FOR_SEPARATE_D2_CHALLENGE_PROPOSAL_ONLY'
        if iip0023_dy_status != 'CLEARED_DESIGN_Y_WILDCARD_COUPONS_ELIGIBLE_FOR_SEPARATE_D2_CHALLENGE_PROPOSAL_ONLY':
            iip0023_hold_status = iip0023_dy_status

    iip0023_eval_summary: dict[str, object] | None = None
    if _iip0023_dy_active:
        iip0023_eval_summary = {
            'gate_pass': iip0023_hold_status is None,
            'hold_status': iip0023_hold_status,
            'wildcard_status': iip0023_dy_status,
            'design_y_integrase_modality': _dy_int_mod,
            'design_y_silencing_modality': _dy_sil_mod,
            'design_y_payload_arm': _dy_pay_arm,
            'auto_promotes_into_core_blueprints': False,
        }

    # S8b TYK2 x inducible-promoter epistasis: support explicit enumerated scale {FRACTION_HALF_EFFECT, PP_HALF_EFFECT, PP_FULL_CONTRAST, DID_PP} with per-endpoint scale retention across batches, unified fractional half-effect threshold (point < -0.05 or lower95 < -0.10), and sign-independent underpowered SE hold (AUDIT-REF & AUDIT-REF; DEC-042/044):
    s8b_scale_enum = observations.get('s8b_epistasis_estimand_scale')
    if s8b_scale_enum is not None:
        if not isinstance(s8b_scale_enum, str) or s8b_scale_enum.strip() not in ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES:
            raise ValueError(
                f"Invalid s8b_epistasis_estimand_scale {s8b_scale_enum!r}; must be one of {sorted(ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES)} (AUDIT-REF)"
            )
        s8b_scale_enum = s8b_scale_enum.strip()

    def _eval_tyk2_beta(mean_key: str, se_key: str, scale_key: str) -> tuple[bool, bool]:
        if mean_key not in validated_stage0:
            return False, False
        raw_b = validated_stage0[mean_key]
        ep_scale_raw = observations.get(scale_key)
        if ep_scale_raw is not None:
            if not isinstance(ep_scale_raw, str) or ep_scale_raw.strip() not in ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES:
                raise ValueError(
                    f"Invalid {scale_key} {ep_scale_raw!r}; must be one of {sorted(ALLOWED_S8B_EPISTASIS_ESTIMAND_SCALES)} (AUDIT-REF)"
                )
            eff_scale = ep_scale_raw.strip()
        else:
            eff_scale = s8b_scale_enum
        if eff_scale == 'FRACTION_HALF_EFFECT':
            scale_div = 1.0
        elif eff_scale == 'PP_HALF_EFFECT':
            scale_div = 100.0
        elif eff_scale == 'PP_FULL_CONTRAST':
            scale_div = 200.0
        elif eff_scale == 'DID_PP':
            scale_div = 400.0
        else:
            is_pp_scale = abs(raw_b) > 1.0 or str(observations.get('estimand_scale', '')).endswith('_pp')
            scale_div = 100.0 if is_pp_scale else 1.0
        b_frac = raw_b / scale_div
        thresh_frac = -0.05
        if se_provenance_by_mean.get(mean_key) == 'EXPLICIT_SE_SUPPLIED':
            raw_s = validated_stage0_se.get(se_key, 0.0)
            s_frac = (raw_s / scale_div) if eff_scale is not None else (raw_s / scale_div if (scale_div > 1.0 or raw_s > 0.5) else raw_s)
            # If SE is so wide (> 0.15 fractional = 15 pp half-effect) that the test cannot rule out a >10 pp deficit, flag as underpowered/inconclusive hold regardless of point sign (AUDIT-REF & AUDIT-REF):
            if s_frac > 0.15:
                return False, True
        else:
            s_frac = float(STAGE0_CUSTOM_DEFAULT_SE.get(mean_key, 0.05))
        b_lower95 = b_frac - 1.645 * s_frac
        return bool(b_frac < thresh_frac or b_lower95 < -0.10), False

    x14_trig, x14_underpowered = _eval_tyk2_beta('obs_beta_x1_x4_epistasis', 'obs_beta_x1_x4_epistasis_se', 'obs_beta_x1_x4_epistasis_scale')
    x34_trig, x34_underpowered = _eval_tyk2_beta('obs_beta_x3_x4_epistasis', 'obs_beta_x3_x4_epistasis_se', 'obs_beta_x3_x4_epistasis_scale')
    sat_trig = bool(_lower95('obs_tyk2_p1104a_3x2_satellite_shield_induction_pct', 95.0) < 70.0)
    is_tyk2_product_or_unspecified = bool(
        product_def is None
        or 'TYK2' in str(product_def.get('product_variant_id', '')).upper()
        or any('TYK2' in str(e).upper() for e in product_def.get('edits', []))
    )
    s8_tyk2_epistasis_triggered = bool(
        is_tyk2_product_or_unspecified and (x14_trig or x34_trig or sat_trig)
    )
    s8_tyk2_epistasis_underpowered_hold = bool(
        is_tyk2_product_or_unspecified and (x14_underpowered or x34_underpowered)
    )
    subthreshold_outer_band_kperf = bool(
        'obs_outer_band_kperf_s_inv' in validated_stage0
        and 1e-5 < _lower95('obs_outer_band_kperf_s_inv') < 0.10
    )
    periimplant_fluid_film_gt_80_failed = bool(
        'obs_periimplant_fluid_film_um' in validated_stage0
        and _upper95('obs_periimplant_fluid_film_um') > 80.0
    )
    periimplant_fluid_film_50_to_80_unmitigated = bool(
        'obs_periimplant_fluid_film_um' in validated_stage0
        and 50.0 < _upper95('obs_periimplant_fluid_film_um') <= 80.0
        and observations.get('periimplant_closed_suction_drain_verified') is not True
    )
    periimplant_fluid_film_failed = bool(
        periimplant_fluid_film_gt_80_failed
        or periimplant_fluid_film_50_to_80_unmitigated
    )
    explant_retrieval_failed = bool(
        ('obs_explant_en_bloc_recovery_fraction' in validated_stage0 and _lower95_student_t('obs_explant_en_bloc_recovery_fraction') < 0.99)
        or ('obs_gate0s_day30_en_bloc_recovery_fraction' in validated_stage0 and _lower95_student_t('obs_gate0s_day30_en_bloc_recovery_fraction') < 0.99)
        or ('obs_gate0s_day90_en_bloc_recovery_fraction' in validated_stage0 and _lower95_student_t('obs_gate0s_day90_en_bloc_recovery_fraction') < 0.99)
        or observations.get('obs_vascularized_explant_en_bloc_recovery_verified') is False
    )
    ap1903_insulin_dump_failed = bool(
        ('obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h' in validated_stage0 and _upper95('obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h') > 0.15)
        or ('obs_gate0s_ap1903_1h_insulin_burst_iu' in validated_stage0 and _upper95('obs_gate0s_ap1903_1h_insulin_burst_iu') > 2.8)
    )
    ap1903_insulin_dump_unmitigated = bool(
        ap1903_insulin_dump_failed
        and observations.get('s9_staged_cassette_excision_and_d10w_clamp_qualified') is not True
    )
    senescent_p16_upper95 = _prop_or_pct_u95(
        'obs_senescent_p16_sasp_positive_cell_pct',
        'obs_senescent_p16_sasp_cells_scored',
        is_pct=True,
    )
    senescent_p16_sasp_failed = bool(
        'obs_senescent_p16_sasp_positive_cell_pct' in validated_stage0
        and senescent_p16_upper95 > 1.0
    )
    periimplant_fluid_film_gt_50_triggered = bool(
        'obs_periimplant_fluid_film_um' in validated_stage0
        and (
            validated_stage0['obs_periimplant_fluid_film_um'] > 50.0
            or _upper95('obs_periimplant_fluid_film_um') > 50.0
        )
    )
    s9_promoted = bool(
        (has_chronic_capsule_obs and capsule_dz_upper95 > 10.4)
        or _lower95('obs_d_o2_gel_aged_cm2_s', 1.5e-5) < 1.50e-5
        or _lower95('obs_capsule_d_o2_cm2_s', 1.5e-5) < 1.50e-5
        or _upper95('obs_cao2_kinetic_law', 0.0) >= 1.5
        or observations.get('obs_masked_perfusion_confirmed') is True
        or observations.get('s9_mineral_pore_occlusion_gt_15pct') is True
        or periimplant_fluid_film_gt_50_triggered
    )
    primary_biological_or_pouch_gate_failed = bool(
        ('obs_pouch_min_po2_mmhg' in validated_stage0 and _lower95('obs_pouch_min_po2_mmhg') < 15.0)
        or ('obs_pouch_max_h2o2_uM' in validated_stage0 and _upper95('obs_pouch_max_h2o2_uM') > 5.0)
        or ('obs_pouch_max_ph' in validated_stage0 and _upper95('obs_pouch_max_ph') > 7.45)
        or ('obs_ifng_hla_dr_dq_negative_pct' in validated_stage0 and _lower95('obs_ifng_hla_dr_dq_negative_pct') < 99.0)
        or ('obs_stage7_yield_ratio' in validated_stage0 and _lower95('obs_stage7_yield_ratio') < 0.70)
        or ('obs_cd47_fold' in validated_stage0 and _lower95('obs_cd47_fold') < 3.0)
        or ('obs_day28_functional_beta_survival_fraction' in validated_stage0 and _lower95('obs_day28_functional_beta_survival_fraction') < 0.80)
    )

    # Stage-0 atomic provenance, donor-bank kinship, Readout-2 kinetics, CD16+ ADCC, coupled catalase endotoxin, Gate 0S dual shedding, 3-metric IAPP, and Design Z route binding (IIP-0024; DEC-103):
    _req_iip0024_full = observations.get('require_full_iip0024_atomic_provenance_and_bank_panel') is True
    _iip0024_quant_or_flag_keys = (
        'obs_mcb_donor_pairwise_ibd_pi_hat',
        'obs_mcb_1q32_mdm4_cnv_ratio',
        'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction',
        'obs_1q32_mdm4_cnv_ratio',
        'obs_1q32_mdm4_single_cell_mosaic_fraction',
        'obs_readout2_perifusion_gsis_si',
        'obs_readout2_first_phase_peak_to_basal_ratio',
        'obs_readout2_first_phase_time_to_peak_min',
        'obs_readout2_stepdown_shutoff_t50_min',
        'obs_exact_product_net_cd16_adcc_lysis_pct',
        'obs_cd16_adcc_positive_control_lysis_pct',
        'obs_p05_total_loaded_catalase_protein_mg',
        'obs_gate0s_shedding_ddpcr_events_k',
        'obs_gate0s_shedding_ddpcr_trials_n',
        'obs_gate0s_shedding_viable_events_k',
        'obs_gate0s_shedding_viable_trials_n',
        'obs_stage0_independent_lot_replicates_n',
        'obs_iapp_3metric_thioflavin_s_area_pct',
        'obs_iapp_3metric_soluble_oligomer_fold',
        'obs_iapp_3metric_insoluble_fibril_fold',
    )
    _iip0024_bool_or_meta_keys = (
        'cd16_adcc_igg_depletion_paired_control_verified',
        'iapp_age_matched_reference_above_loq_verified',
        'mcb_imprinting_ic1_ic2_igdmr_conjunctive_cleared',
        'stage7_lot_imprinting_ic1_ic2_igdmr_conjunctive_cleared',
        'mcb_and_lot_xci_46xx_3axis_cleared_or_46xy',
        'mcb_and_lot_mtdna_1000x_numt_filtered_cleared',
        'physical_route_id',
        'donor_seedline_id',
        'device_assembly_id',
    )
    _is_design_z_product = bool(
        product_def is not None
        and str(product_def.get('product_variant_id', '')) == 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2'
    )
    _iip0024_active = bool(
        _req_iip0024_full
        or _is_design_z_product
        or any(_k in validated_stage0 for _k in _iip0024_quant_or_flag_keys)
        or any(observations.get(_bk) is not None for _bk in _iip0024_bool_or_meta_keys)
    )
    iip0024_hold_status: str | None = None
    _phys_route_in = observations.get('physical_route_id')
    if isinstance(_phys_route_in, str):
        _phys_route_in = _phys_route_in.strip()
    if _phys_route_in is not None:
        _is_z_route = _phys_route_in in ('Z_OMENTAL_140UM_ZERO_CAO2', 'DESIGN_Z_140UM_ZERO_CAO2_PLAIN_HYDROGEL', 'DESIGN_Z_140UM_ZERO_CAO2_EXPLORATORY_DIVIDER')
        _is_cao2_route = _phys_route_in == 'SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH'
        if not (_is_z_route or _is_cao2_route):
            iip0024_hold_status = 'HOLD_UNREGISTERED_PHYSICAL_ROUTE_ID'
        elif _is_design_z_product and not _is_z_route:
            iip0024_hold_status = 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH'
        elif product_def is not None and not _is_design_z_product and product_def.get('family') in ('A', 'B') and _is_z_route:
            iip0024_hold_status = 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH'
    if iip0024_hold_status is None and product_def is not None and not _is_design_z_product and (
        _dz_arm_val is not None or any(_k.startswith('obs_design_z_') for _k in validated_stage0)
    ):
        iip0024_hold_status = 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH'
    if iip0024_hold_status is None and _is_design_z_product and (
        observations.get('p05_catalase_uniprot_accession') is not None
        or any(_k.startswith('obs_p05_') for _k in validated_stage0)
    ):
        iip0024_hold_status = 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH'
    if iip0024_hold_status is None and _is_design_z_product:
        if iip0019_dz_status != 'CLEARED_DESIGN_Z_TWO_STAGE_EMPTY_LUMEN_AND_LOADED_CLOSED_STATE_GATE':
            iip0024_hold_status = iip0019_transport_hold_status or 'HOLD_DESIGN_Z_WILDCARD_TWO_STAGE_QUALIFICATION_REQUIRED'
    if iip0024_hold_status is None and stage0_replicates_below_6_hold:
        iip0024_hold_status = 'HOLD_STAGE0_INDEPENDENT_LOT_REPLICATES_BELOW_6_STUDENT_T_MINIMUM'
    if iip0024_hold_status is None and 'obs_mcb_donor_pairwise_ibd_pi_hat' in validated_stage0:
        if validated_stage0['obs_mcb_donor_pairwise_ibd_pi_hat'] >= 0.05 or _upper95('obs_mcb_donor_pairwise_ibd_pi_hat') >= 0.05:
            iip0024_hold_status = 'HOLD_MCB_DONOR_SEEDLINE_PAIRWISE_IBD_PI_HAT_GE_0P05'
    if iip0024_hold_status is None and any(
        _rk in validated_stage0
        for _rk in (
            'obs_readout2_perifusion_gsis_si',
            'obs_readout2_first_phase_peak_to_basal_ratio',
            'obs_readout2_first_phase_time_to_peak_min',
            'obs_readout2_stepdown_shutoff_t50_min',
        )
    ):
        if (
            ('obs_readout2_perifusion_gsis_si' in validated_stage0 and _lower95('obs_readout2_perifusion_gsis_si') < 3.0)
            or ('obs_readout2_first_phase_peak_to_basal_ratio' in validated_stage0 and _lower95('obs_readout2_first_phase_peak_to_basal_ratio') < 4.5)
            or ('obs_readout2_first_phase_time_to_peak_min' in validated_stage0 and _upper95('obs_readout2_first_phase_time_to_peak_min') > 8.0)
            or ('obs_readout2_stepdown_shutoff_t50_min' in validated_stage0 and _upper95('obs_readout2_stepdown_shutoff_t50_min') > 3.0)
            or (
                'obs_readout2_stepdown_shutoff_t50_min' in validated_stage0
                and 'obs_gate0m_shutoff_t90_min' in validated_stage0
                and validated_stage0['obs_readout2_stepdown_shutoff_t50_min'] > validated_stage0['obs_gate0m_shutoff_t90_min']
            )
        ):
            iip0024_hold_status = 'HOLD_READOUT2_PERIFUSION_KINETICS_OR_SHUTOFF_MONOTONICITY_FAILED'
    if iip0024_hold_status is None and (
        'obs_exact_product_net_cd16_adcc_lysis_pct' in validated_stage0
        or 'obs_cd16_adcc_positive_control_lysis_pct' in validated_stage0
        or observations.get('cd16_adcc_igg_depletion_paired_control_verified') is not None
    ):
        _has_any_net_adcc = (
            'obs_exact_product_net_cd16_adcc_lysis_pct' in validated_stage0
            or 'obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct' in validated_stage0
        )
        if (
            observations.get('cd16_adcc_igg_depletion_paired_control_verified') is not True
            or not _has_any_net_adcc
            or 'obs_cd16_adcc_positive_control_lysis_pct' not in validated_stage0
            or (
                'obs_exact_product_net_cd16_adcc_lysis_pct' in validated_stage0
                and _upper95('obs_exact_product_net_cd16_adcc_lysis_pct') > 5.0
            )
            or (
                'obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct' in validated_stage0
                and _upper95('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct') > 5.0
            )
            or _lower95('obs_cd16_adcc_positive_control_lysis_pct') < 25.0
        ):
            iip0024_hold_status = 'HOLD_PAIRED_IGG_DEPLETED_CD16_ADCC_GATE_FAILED_OR_UNVERIFIED'
    if iip0024_hold_status is None and 'obs_p05_total_loaded_catalase_protein_mg' in validated_stage0 and not _is_design_z_product:
        _m_cat = validated_stage0['obs_p05_total_loaded_catalase_protein_mg']
        if _m_cat <= 0.0:
            iip0024_hold_status = 'HOLD_P05_TOTAL_LOADED_CATALASE_MASS_NONPOSITIVE'
        elif 'obs_p05_catalase_endotoxin_eu_per_mg' in validated_stage0:
            _coupled_eu_cap = min(0.10, 0.50 / _m_cat)
            if _upper95_student_t('obs_p05_catalase_endotoxin_eu_per_mg') > _coupled_eu_cap:
                iip0024_hold_status = 'HOLD_P05_MASS_COUPLED_CATALASE_ENDOTOXIN_CAP_EXCEEDED'
            elif 'obs_p05_total_implant_catalase_endotoxin_eu' in validated_stage0:
                if _upper95_student_t('obs_p05_total_implant_catalase_endotoxin_eu') > 0.50:
                    iip0024_hold_status = 'HOLD_P05_TOTAL_IMPLANT_CATALASE_ENDOTOXIN_GT_0P50EU'
                else:
                    _exp_tot_eu = _m_cat * validated_stage0['obs_p05_catalase_endotoxin_eu_per_mg']
                    if abs(validated_stage0['obs_p05_total_implant_catalase_endotoxin_eu'] - _exp_tot_eu) > max(0.015, 0.05 * _exp_tot_eu):
                        iip0024_hold_status = 'HOLD_P05_COUPLED_TOTAL_ENDOTOXIN_MASS_PRODUCT_INCONSISTENT'
    if iip0024_hold_status is None and any(
        _sk in validated_stage0
        for _sk in (
            'obs_gate0s_shedding_ddpcr_events_k',
            'obs_gate0s_shedding_ddpcr_trials_n',
            'obs_gate0s_shedding_viable_events_k',
            'obs_gate0s_shedding_viable_trials_n',
        )
    ):
        if not all(
            _sk in validated_stage0
            for _sk in (
                'obs_gate0s_shedding_ddpcr_events_k',
                'obs_gate0s_shedding_ddpcr_trials_n',
                'obs_gate0s_shedding_viable_events_k',
                'obs_gate0s_shedding_viable_trials_n',
            )
        ):
            iip0024_hold_status = 'HOLD_GATE0S_DUAL_SHEDDING_PANEL_INCOMPLETE'
        else:
            _k_dd_f = validated_stage0['obs_gate0s_shedding_ddpcr_events_k']
            _n_dd_f = validated_stage0['obs_gate0s_shedding_ddpcr_trials_n']
            _k_vi_f = validated_stage0['obs_gate0s_shedding_viable_events_k']
            _n_vi_f = validated_stage0['obs_gate0s_shedding_viable_trials_n']
            if any(abs(_vf - round(_vf)) > 1e-9 for _vf in (_k_dd_f, _n_dd_f, _k_vi_f, _n_vi_f)):
                non_integer_scored_cell_or_trial_denominator = True
                iip0024_hold_status = 'HOLD_NON_INTEGER_SCORED_CELL_OR_TRIAL_DENOMINATOR'
            else:
                _spec_map_s0 = observations.get('specimen_or_clone_ids')
                if (
                    isinstance(_spec_map_s0, dict)
                    and 'obs_gate0s_shedding_ddpcr_events_k' in _spec_map_s0
                    and 'obs_gate0s_shedding_viable_events_k' in _spec_map_s0
                    and str(_spec_map_s0['obs_gate0s_shedding_ddpcr_events_k']).strip()
                    == str(_spec_map_s0['obs_gate0s_shedding_viable_events_k']).strip()
                ):
                    iip0024_hold_status = 'HOLD_GATE0S_DUAL_SHEDDING_REQUIRES_INDEPENDENT_SPECIMEN_UNITS'
                else:
                    _n_dd = int(round(_n_dd_f))
                    _n_vi = int(round(_n_vi_f))
                    if (
                        _n_dd < 2995731
                        or _n_vi < 2995731
                        or clopper_pearson_upper95_independent(0, _n_dd) >= 1e-6
                        or clopper_pearson_upper95_independent(0, _n_vi) >= 1e-6
                    ):
                        iip0024_hold_status = 'HOLD_GATE0S_DUAL_SHEDDING_TRIALS_BELOW_2995731_BINOMIAL_1E6_FLOOR'
    if iip0024_hold_status is None and (
        any(_ik in validated_stage0 for _ik in ('obs_iapp_3metric_thioflavin_s_area_pct', 'obs_iapp_3metric_soluble_oligomer_fold', 'obs_iapp_3metric_insoluble_fibril_fold'))
        or (
            observations.get('iapp_age_matched_reference_above_loq_verified') is not None
            and not any(
                _ik28 in validated_stage0
                for _ik28 in (
                    'obs_iapp_72h_immunocapture_oligomer_elisa_fold',
                    'obs_iapp_72h_colocalized_fibril_area_pct',
                    'obs_iapp_d28_immunocapture_oligomer_elisa_fold',
                    'obs_iapp_d28_colocalized_fibril_area_pct',
                    'obs_iapp_d84_immunocapture_oligomer_elisa_fold',
                    'obs_iapp_d84_colocalized_fibril_area_pct',
                )
            )
        )
    ):
        if (
            observations.get('iapp_age_matched_reference_above_loq_verified') is not True
            or 'obs_iapp_3metric_thioflavin_s_area_pct' not in validated_stage0
            or 'obs_iapp_3metric_soluble_oligomer_fold' not in validated_stage0
            or 'obs_iapp_3metric_insoluble_fibril_fold' not in validated_stage0
            or _upper95('obs_iapp_3metric_thioflavin_s_area_pct') > 1.0
            or _upper95('obs_iapp_3metric_soluble_oligomer_fold') > 1.5
            or _upper95('obs_iapp_3metric_insoluble_fibril_fold') > 1.5
        ):
            iip0024_hold_status = 'HOLD_IAPP_3METRIC_PANEL_OR_REFERENCE_LOQ_FAILED'
    if iip0024_hold_status is None:
        for _epi_flag_k, _epi_hold_code in (
            ('mcb_imprinting_ic1_ic2_igdmr_conjunctive_cleared', 'HOLD_MCB_CONJUNCTIVE_IMPRINTING_DMR_NOT_CLEARED'),
            ('stage7_lot_imprinting_ic1_ic2_igdmr_conjunctive_cleared', 'HOLD_STAGE7_LOT_CONJUNCTIVE_IMPRINTING_DMR_NOT_CLEARED'),
            ('mcb_and_lot_xci_46xx_3axis_cleared_or_46xy', 'HOLD_46XX_XCI_3AXIS_NOT_CLEARED'),
            ('mcb_and_lot_mtdna_1000x_numt_filtered_cleared', 'HOLD_SCOPED_1000X_NUMT_MTDNA_NOT_CLEARED'),
        ):
            if observations.get(_epi_flag_k) is False or (_req_iip0024_full and observations.get(_epi_flag_k) is not True):
                iip0024_hold_status = _epi_hold_code
                break
    if iip0024_hold_status is None and _req_iip0024_full:
        _req_24_keys = (
            'obs_mcb_donor_pairwise_ibd_pi_hat',
            'obs_mcb_1q32_mdm4_cnv_ratio',
            'obs_mcb_1q32_mdm4_single_cell_mosaic_fraction',
            'obs_mcb_1q32_mdm4_single_cell_mosaic_cells_scored',
            'obs_1q32_mdm4_cnv_ratio',
            'obs_1q32_mdm4_single_cell_mosaic_fraction',
            'obs_1q32_mdm4_single_cell_mosaic_cells_scored',
            'obs_readout2_perifusion_gsis_si',
            'obs_readout2_first_phase_peak_to_basal_ratio',
            'obs_readout2_first_phase_time_to_peak_min',
            'obs_readout2_stepdown_shutoff_t50_min',
            'obs_exact_product_net_cd16_adcc_lysis_pct',
            'obs_cd16_adcc_positive_control_lysis_pct',
            'obs_gate0s_shedding_ddpcr_events_k',
            'obs_gate0s_shedding_ddpcr_trials_n',
            'obs_gate0s_shedding_viable_events_k',
            'obs_gate0s_shedding_viable_trials_n',
            'obs_iapp_3metric_thioflavin_s_area_pct',
            'obs_iapp_3metric_soluble_oligomer_fold',
            'obs_iapp_3metric_insoluble_fibril_fold',
        )
        if not all(_k in validated_stage0 for _k in _req_24_keys) or (not _is_design_z_product and 'obs_p05_total_loaded_catalase_protein_mg' not in validated_stage0):
            iip0024_hold_status = 'HOLD_IIP0024_MANDATORY_ATOMIC_PROVENANCE_OR_BANK_PANEL_INCOMPLETE'
    iip0024_eval_summary: dict[str, object] | None = None
    if _iip0024_active:
        iip0024_eval_summary = {
            'gate_pass': iip0024_hold_status is None,
            'hold_status': iip0024_hold_status,
            'status': iip0024_hold_status if iip0024_hold_status is not None else 'CLEARED_IIP0024_STAGE0_ATOMIC_PROVENANCE_AND_BANK_GUARDRAILS',
            'student_t_critical_multiplier_used': round(crit_95_mult, 8),
        }

    # Stage-0 Hurdle-4 two-timepoint human P04040 catalase envelope hold, exact z=70.0 um gel-interface H2O2, separated Day-30/Day-90 Gate 0S retrieval, and zero-primate Contrasts B1, B2, B3 (IIP-0025; DEC-104):
    _b1_arm = observations.get('contrast_b1_geometry_arm')
    if isinstance(_b1_arm, str):
        _b1_arm = _b1_arm.strip()
    _b2_arm = observations.get('contrast_b2_perfusion_and_retrieval_arm')
    if isinstance(_b2_arm, str):
        _b2_arm = _b2_arm.strip()
    _b3_arm = observations.get('contrast_b3_donor_iec_arm')
    if isinstance(_b3_arm, str):
        _b3_arm = _b3_arm.strip()

    _iip0025_p05_new_keys = (
        'obs_p05_precure_active_catalase_input_u_333cm2',
        'obs_p05_postcure_t0_active_catalase_iu_333cm2',
        'obs_p05_cure_immobilization_yield_pct',
        'obs_p05_336h_catalase_retention_pct',
        'obs_p05_z70_gel_interface_h2o2_max_uM',
    )
    _iip0025_p05_all_req_keys = (
        *_iip0025_p05_new_keys,
        'obs_p05_day7_catalase_retention_pct',
        'obs_p05_day7_shed_catalase_mass_pct',
        'obs_p05_336h_retained_catalase_iu_333cm2',
        'obs_p05_total_loaded_catalase_protein_mg',
        'obs_p05_catalase_endotoxin_eu_per_mg',
        'obs_p05_total_implant_catalase_endotoxin_eu',
    )
    _iip0025_b1_keys = (
        'obs_contrast_b1_paired_min_face_po2_gain_mmhg',
        'obs_contrast_b1_80um_worst_face_po2_mmhg',
        'obs_contrast_b1_80um_perifusion_gsis_si',
    )
    _iip0025_b2_keys = (
        'obs_contrast_b2_human_exvivo_viable_window_h',
        'obs_contrast_b2_omental_vs_suprafascial_po2_gain_mmhg',
        'obs_contrast_b2_porcine_day90_en_bloc_recovery_fraction',
    )
    _iip0025_b3_keys = (
        'obs_contrast_b3_iec_omitted_surface_po2_mmhg',
        'obs_contrast_b3_iec_omitted_outer_kperf_s_inv',
        'obs_contrast_b3_iec_omitted_perifusion_gsis_si',
    )
    _iip0025_p05_env_active = bool(
        not _is_design_z_product
        and (
            _req_iip0025_full
            or observations.get('p05_two_timepoint_same_lot_envelope_verified') is not None
            or any(_k in validated_stage0 for _k in _iip0025_p05_new_keys)
        )
    )
    _iip0025_b1_active = bool(
        _b1_arm is not None
        or any(_k in validated_stage0 for _k in _iip0025_b1_keys)
    )
    _iip0025_b2_active = bool(
        _b2_arm is not None
        or observations.get('contrast_b2_route_qualification_distinct_from_portfolio_credit_verified') is not None
        or any(_k in validated_stage0 for _k in _iip0025_b2_keys)
    )
    _iip0025_b3_active = bool(
        _b3_arm is not None
        or observations.get('contrast_b3_zero_primate_host_ingrowth_verified') is not None
        or any(_k in validated_stage0 for _k in _iip0025_b3_keys)
    )
    _iip0025_day30_90_active = bool(
        'obs_gate0s_day30_en_bloc_recovery_fraction' in validated_stage0
        or 'obs_gate0s_day90_en_bloc_recovery_fraction' in validated_stage0
    )

    iip0025_p05_env_status: str | None = None
    iip0025_b1_status: str | None = None
    iip0025_b2_status: str | None = None
    iip0025_b3_status: str | None = None
    iip0025_hold_status: str | None = None

    if _is_design_z_product and (
        observations.get('p05_two_timepoint_same_lot_envelope_verified') is not None
        or any(_k in validated_stage0 for _k in _iip0025_p05_new_keys)
    ):
        iip0025_hold_status = 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH'

    if _iip0025_p05_env_active and iip0025_hold_status is None:
        _p05_25_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _iip0025_p05_all_req_keys
        )
        if _p05_cat_acc != 'P04040' or observations.get('p05_recombinant_human_catalase_verified') is not True:
            iip0025_p05_env_status = 'HOLD_P05_XENO_OR_UNVERIFIED_CATALASE_SPECIES_REQUIRES_HUMAN_P04040'
        elif stage0_replicates_below_6_hold:
            iip0025_p05_env_status = 'HOLD_STAGE0_INDEPENDENT_LOT_REPLICATES_BELOW_6_STUDENT_T_MINIMUM'
        elif p05_z70_interface_h2o2_hard_stop_veto:
            iip0025_p05_env_status = 'HOLD_P05_Z70_GEL_INTERFACE_H2O2_EXCEEDS_5P0UM_HARD_STOP'
        elif 'obs_p05_z70_gel_interface_h2o2_max_uM' in validated_stage0 and _upper95_student_t('obs_p05_z70_gel_interface_h2o2_max_uM') > 2.5:
            iip0025_p05_env_status = 'HOLD_P05_Z70_GEL_INTERFACE_H2O2_EXCEEDS_2P5UM_TIER1_TARGET'
        elif observations.get('p05_two_timepoint_same_lot_envelope_verified') is not True or not _p05_25_quant_ok:
            iip0025_p05_env_status = 'HOLD_P05_TWO_TIMEPOINT_ENVELOPE_INCOMPLETE_OR_DEFAULT_SE'
        else:
            _u_pre = validated_stage0['obs_p05_precure_active_catalase_input_u_333cm2']
            _u_post0 = validated_stage0['obs_p05_postcure_t0_active_catalase_iu_333cm2']
            _y_cure = validated_stage0['obs_p05_cure_immobilization_yield_pct']
            _r_d7 = validated_stage0['obs_p05_day7_catalase_retention_pct']
            _s_d7 = validated_stage0['obs_p05_day7_shed_catalase_mass_pct']
            _r_336 = validated_stage0['obs_p05_336h_catalase_retention_pct']
            _u_336 = validated_stage0['obs_p05_336h_retained_catalase_iu_333cm2']
            _m_cat_25 = validated_stage0['obs_p05_total_loaded_catalase_protein_mg']
            if _u_pre > 158300.0 or _upper95_student_t('obs_p05_precure_active_catalase_input_u_333cm2') > 158300.0:
                iip0025_p05_env_status = 'HOLD_P05_PRECURE_CATALASE_INPUT_EXCEEDS_158300U_CAP'
            elif (
                _u_pre <= 0.0
                or _u_post0 <= 0.0
                or _m_cat_25 <= 0.0
                or abs(_u_post0 - _u_pre * (_y_cure / 100.0)) > max(150.0, 0.02 * _u_post0)
                or abs(_u_336 - _u_post0 * (_r_336 / 100.0)) > max(150.0, 0.02 * max(1.0, _u_336))
                or _r_336 > _r_d7 + 1e-6
                or (_r_d7 + _s_d7) > 102.0
            ):
                iip0025_p05_env_status = 'HOLD_P05_TWO_TIMEPOINT_ALGEBRAIC_YIELD_OR_RETENTION_INCONSISTENT'
            elif (
                _lower95_student_t('obs_p05_day7_catalase_retention_pct') < 70.0
                or _upper95_student_t('obs_p05_day7_shed_catalase_mass_pct') > 30.0
                or _lower95_student_t('obs_p05_336h_retained_catalase_iu_333cm2') < 73800.0
            ):
                iip0025_p05_env_status = 'HOLD_P05_TWO_TIMEPOINT_DAY7_70PCT_OR_336H_73800IU_FAILED'
            else:
                _unrounded_eu_cap = min(0.10, 0.50 / _m_cat_25)
                _exp_tot_eu_25 = _m_cat_25 * validated_stage0['obs_p05_catalase_endotoxin_eu_per_mg']
                if _upper95_student_t('obs_p05_catalase_endotoxin_eu_per_mg') > _unrounded_eu_cap:
                    iip0025_p05_env_status = 'HOLD_P05_MASS_COUPLED_CATALASE_ENDOTOXIN_CAP_EXCEEDED'
                elif _upper95_student_t('obs_p05_total_implant_catalase_endotoxin_eu') > 0.50:
                    iip0025_p05_env_status = 'HOLD_P05_TOTAL_IMPLANT_CATALASE_ENDOTOXIN_GT_0P50EU'
                elif abs(validated_stage0['obs_p05_total_implant_catalase_endotoxin_eu'] - _exp_tot_eu_25) > max(0.015, 0.05 * _exp_tot_eu_25):
                    iip0025_p05_env_status = 'HOLD_P05_COUPLED_TOTAL_ENDOTOXIN_MASS_PRODUCT_INCONSISTENT'
                else:
                    iip0025_p05_env_status = 'CLEARED_P05_TWO_TIMEPOINT_P04040_ENVELOPE_AND_Z70_H2O2_GATE'
        if iip0025_p05_env_status != 'CLEARED_P05_TWO_TIMEPOINT_P04040_ENVELOPE_AND_Z70_H2O2_GATE':
            iip0025_hold_status = iip0025_p05_env_status

    if _iip0025_b1_active:
        _b1_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _iip0025_b1_keys
        )
        if _b1_arm not in ALLOWED_CONTRAST_B1_GEOMETRY_ARMS or not _b1_quant_ok:
            iip0025_b1_status = 'HOLD_CONTRAST_B1_INCOMPLETE_OR_DEFAULT_SE_PANEL'
        elif observations.get('s6_dose_matched_1500_ieq_verified') is not True:
            iip0025_b1_status = 'HOLD_CONTRAST_B1_REQUIRES_DOSE_MATCHED_1500IEQ_PITCH'
        elif (
            _lower95_student_t('obs_contrast_b1_80um_worst_face_po2_mmhg') < 15.0
            or _lower95_student_t('obs_contrast_b1_paired_min_face_po2_gain_mmhg') <= 0.0
            or _lower95_student_t('obs_contrast_b1_80um_perifusion_gsis_si') < 3.0
        ):
            iip0025_b1_status = 'HOLD_CONTRAST_B1_DOSE_MATCHED_80UM_MIN_FACE_OR_GSIS_FAILED'
        else:
            iip0025_b1_status = 'CLEARED_CONTRAST_B1_DOSE_MATCHED_80UM_VS_100UM_CHALLENGER'
        if iip0025_b1_status != 'CLEARED_CONTRAST_B1_DOSE_MATCHED_80UM_VS_100UM_CHALLENGER' and iip0025_hold_status is None:
            iip0025_hold_status = iip0025_b1_status

    if _iip0025_b2_active:
        _b2_quant_ok = all(
            _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
            for _mk in _iip0025_b2_keys
        )
        if _b2_arm != 'PAIRED_HUMAN_EXVIVO_AND_PORCINE_DAY30_90':
            iip0025_b2_status = 'HOLD_CONTRAST_B2_INCOMPLETE_EXVIVO_PERFUSION_OR_PORCINE_RETRIEVAL_ARM'
        elif observations.get('contrast_b2_route_qualification_distinct_from_portfolio_credit_verified') is not True:
            iip0025_b2_status = 'HOLD_CONTRAST_B2_ROUTE_QUALIFICATION_CONFLATED_WITH_PORTFOLIO_CREDIT'
        elif not _b2_quant_ok:
            iip0025_b2_status = 'HOLD_CONTRAST_B2_INCOMPLETE_OR_DEFAULT_SE_PANEL'
        elif (
            _lower95_student_t('obs_contrast_b2_human_exvivo_viable_window_h') < 24.0
            or _lower95_student_t('obs_contrast_b2_omental_vs_suprafascial_po2_gain_mmhg') <= 0.0
            or _lower95_student_t('obs_contrast_b2_porcine_day90_en_bloc_recovery_fraction') < 0.99
        ):
            iip0025_b2_status = 'HOLD_CONTRAST_B2_PERFUSION_WINDOW_PO2_GAIN_OR_PORCINE_RETRIEVAL_FAILED'
        else:
            iip0025_b2_status = 'CLEARED_CONTRAST_B2_OMENTAL_VS_SUPRAFASCIAL_ZERO_PRIMATE_ROUTE_GATE'
        if iip0025_b2_status != 'CLEARED_CONTRAST_B2_OMENTAL_VS_SUPRAFASCIAL_ZERO_PRIMATE_ROUTE_GATE' and iip0025_hold_status is None:
            iip0025_hold_status = iip0025_b2_status

    if _iip0025_b3_active:
        if _b3_arm not in ALLOWED_CONTRAST_B3_DONOR_IEC_ARMS:
            iip0025_b3_status = 'HOLD_CONTRAST_B3_MISSING_OR_INVALID_DONOR_IEC_ARM'
        elif observations.get('contrast_b3_zero_primate_host_ingrowth_verified') is not True:
            iip0025_b3_status = 'HOLD_CONTRAST_B3_HOST_VESSEL_INGROWTH_UNVERIFIED'
        elif _b3_arm == 'OMIT_DONOR_IEC_8E8_CHALLENGER':
            _b3_quant_ok = all(
                _mk in validated_stage0 and se_provenance_by_mean.get(_mk) == 'EXPLICIT_SE_SUPPLIED'
                for _mk in _iip0025_b3_keys
            )
            if not _b3_quant_ok:
                iip0025_b3_status = 'HOLD_CONTRAST_B3_INCOMPLETE_OR_DEFAULT_SE_PANEL'
            elif (
                _lower95_student_t('obs_contrast_b3_iec_omitted_surface_po2_mmhg') < 15.0
                or _lower95_student_t('obs_contrast_b3_iec_omitted_outer_kperf_s_inv') < 0.10
                or _lower95_student_t('obs_contrast_b3_iec_omitted_perifusion_gsis_si') < 3.0
            ):
                iip0025_b3_status = 'HOLD_CONTRAST_B3_DONOR_IEC_OMITTED_INGROWTH_PO2_OR_GSIS_FAILED'
            else:
                iip0025_b3_status = 'CLEARED_CONTRAST_B3_OMIT_DONOR_IEC_8E8_CHALLENGER'
        else:
            iip0025_b3_status = 'CLEARED_CONTRAST_B3_BASELINE_WITH_DONOR_IEC_9E8_REFERENCE'
        if not str(iip0025_b3_status).startswith('CLEARED_') and iip0025_hold_status is None:
            iip0025_hold_status = iip0025_b3_status

    iip0025_eval_summary: dict[str, object] | None = None
    if _iip0025_p05_env_active or _iip0025_b1_active or _iip0025_b2_active or _iip0025_b3_active or _iip0025_day30_90_active:
        iip0025_eval_summary = {
            'gate_pass': iip0025_hold_status is None and not explant_retrieval_failed,
            'hold_status': iip0025_hold_status,
            'status': (
                iip0025_hold_status
                if iip0025_hold_status is not None
                else (
                    'HOLD_S9_EXPLANT_EN_BLOC_RECOVERY_BELOW_99PCT'
                    if explant_retrieval_failed
                    else 'CLEARED_IIP0025_HURDLE4_SPATIAL_OXYGEN_CATALASE_AND_CHALLENGERS'
                )
            ),
            'p05_two_timepoint_envelope_status': iip0025_p05_env_status,
            'contrast_b1_status': iip0025_b1_status,
            'contrast_b2_status': iip0025_b2_status,
            'contrast_b3_status': iip0025_b3_status,
            'student_t_critical_multiplier_used': round(_student_t_mult_df5_default, 8),
        }

    # Stage-0 molecular provenance, 4-state Design B c1+c5 registry binding, Corrah 2011 DT9/HLA-E control, per-locus segmented BAF, Locus-5 RACE/UMI >=3,689, & TYK2 PE4max MSI/piecewise guardrails (IIP-0026; DEC-105):
    _iip0026_obs_keys = (
        'obs_design_b_hlac_phased_guide_cfd_max',
        'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio',
        'obs_post_edit_segmented_baf_chr19q',
        'obs_post_edit_segmented_baf_chr11p15',
        'obs_post_edit_segmented_baf_chr6p21',
        'obs_design_b_locus5_5prime_race_coord_bp',
        'obs_design_b_locus5_3prime_race_coord_bp',
        'obs_design_b_locus5_race_umi_n_eff',
        'obs_design_b_locus5_aberrant_umi_k',
        'obs_tyk2_pe4max_pentaplex_msi_shifts_k',
        'obs_tyk2_pe4max_wgs_repeat_indels_clone',
        'obs_tyk2_pe4max_wgs_repeat_indels_ctrl',
        'obs_tyk2_pe4max_coding_or_splice_de_novo_indels',
    )
    _iip0026_bool_or_meta_keys = (
        'require_full_iip0026_molecular_provenance_panel',
        'donor_seed_gnomad_target_site_wgs_verified',
        'post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified',
        'design_b_hla_ab_separated_pass_delivery_verified',
        'design_b_chr6p21_1p41mb_loh_sv_cleared',
        'design_b_hlac_dt9_positive',
        'design_b_hlae_3d12_positive',
        'design_b_hlac_allele_specific_surface_or_phased_mrna_verified',
        'b_sct_ciita_exon5_and_trimer_verified',
        'tyk2_phased_haplotype_sap_verified',
        'bound_c1_c5_pair_sha256',
        'tyk2_pe_modality',
        'iip0026_schema_era',
    )
    _req_iip0026_full = bool(observations.get('require_full_iip0026_molecular_provenance_panel') is True)
    _iip0026_active = bool(
        _req_iip0026_full
        or any(_k in validated_stage0 for _k in _iip0026_obs_keys)
        or any(observations.get(_k) is not None for _k in _iip0026_bool_or_meta_keys)
    )
    _pvid_s0_26 = str(product_def['product_variant_id']) if product_def is not None else ''
    _is_endogenous_hlac_b_26 = _pvid_s0_26 in ('ISLET1_B2_5LOCUS_BASE', 'ISLET1_B2_5LOCUS_TYK2_PE', 'ISLET1_B1_6LOCUS_TYK2_KO_INS')
    _is_bsct_26 = _pvid_s0_26.startswith('ISLET1_B_SCT_')
    _is_tyk2_pe_26 = (_pvid_s0_26 == 'ISLET1_B2_5LOCUS_TYK2_PE') or ('TYK2_PE' in _pvid_s0_26)
    iip0026_hold_status: str | None = None
    iip0026_hlac_status: str | None = None
    if _iip0026_active:
        _expected_b2_sha = (
            str(product_def.get('staged_design_b_donor_reference_binding', {}).get('bound_c1_c5_pair_sha256', '')).strip()
            if isinstance(product_def, dict) and isinstance(product_def.get('staged_design_b_donor_reference_binding'), dict)
            else ''
        )
        if observations.get('iip0026_schema_era') == 'LEGACY_PRE_IIP0026':
            iip0026_hold_status = 'HOLD_PROVISIONAL_LEGACY_SCHEMA_MISSING_IIP0026_EVIDENCE'
        elif _expected_b2_sha:
            _sup_sha = observations.get('bound_c1_c5_pair_sha256')
            if _sup_sha is None:
                if _req_iip0026_full:
                    iip0026_hold_status = 'HOLD_MISSING_BOUND_C1_C5_DONOR_PAIR_DIGEST'
            elif _sup_sha.strip() != _expected_b2_sha:
                iip0026_hold_status = 'HOLD_DONOR_REGISTRY_DIGEST_MISMATCH'
        if iip0026_hold_status is None and (
            _req_iip0026_full or observations.get('donor_seed_gnomad_target_site_wgs_verified') is False
        ):
            if observations.get('donor_seed_gnomad_target_site_wgs_verified') is not True:
                iip0026_hold_status = 'HOLD_MISSING_DONOR_TARGET_WGS'
        if iip0026_hold_status is None and (
            _req_iip0026_full
            or any(_bk in validated_stage0 for _bk in ('obs_post_edit_segmented_baf_chr19q', 'obs_post_edit_segmented_baf_chr11p15', 'obs_post_edit_segmented_baf_chr6p21'))
            or observations.get('post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified') is not None
        ):
            _req_baf_keys = (
                ('obs_post_edit_segmented_baf_chr19q', 'obs_post_edit_segmented_baf_chr11p15', 'obs_post_edit_segmented_baf_chr6p21')
                if _is_endogenous_hlac_b_26
                else ('obs_post_edit_segmented_baf_chr19q', 'obs_post_edit_segmented_baf_chr11p15')
            )
            if (
                any(_bk not in validated_stage0 for _bk in _req_baf_keys)
                or observations.get('post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified') is not True
            ):
                iip0026_hold_status = 'HOLD_MISSING_SEGMENTED_PHASED_BAF_SURVEILLANCE'
            elif max(abs(validated_stage0[_bk] - 0.50) for _bk in _req_baf_keys) > 0.05:
                iip0026_hold_status = 'HOLD_DISCORDANT_PER_LOCUS_SEGMENTED_BAF'
        if _is_endogenous_hlac_b_26:
            _cfd_26 = validated_stage0.get('obs_design_b_hlac_phased_guide_cfd_max')
            _dt9_26 = validated_stage0.get('obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio')
            _orth_26 = observations.get('design_b_hlac_allele_specific_surface_or_phased_mrna_verified')
            if (
                _cfd_26 is None
                or not (0.0 <= _cfd_26 < 0.05)
                or observations.get('design_b_hla_ab_separated_pass_delivery_verified') is not True
                or observations.get('design_b_chr6p21_1p41mb_loh_sv_cleared') is not True
                or observations.get('design_b_hlac_dt9_positive') is not True
                or observations.get('design_b_hlae_3d12_positive') is not True
            ):
                iip0026_hlac_status = 'HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE'
            elif _dt9_26 is None and _orth_26 is not True:
                iip0026_hlac_status = 'HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY'
            elif _dt9_26 is not None and _dt9_26 <= 0.0 and _orth_26 is not True:
                iip0026_hlac_status = 'HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY'
            elif _dt9_26 is not None and _dt9_26 < 2.0 and _orth_26 is not True:
                iip0026_hlac_status = 'HOLD_HLAC_SURFACE_LOSS_MASKED_BY_DT9_HLAE_CROSSREACTIVITY'
            else:
                iip0026_hlac_status = 'CLEARED_ENDOGENOUS_DESIGNB_HLAC_PHASED_AND_DT9_HLAE_CONTROL'
            if iip0026_hold_status is None and not str(iip0026_hlac_status).startswith('CLEARED_'):
                iip0026_hold_status = iip0026_hlac_status
            if iip0026_hold_status is None and (
                _req_iip0026_full
                or any(
                    _rk in validated_stage0
                    for _rk in (
                        'obs_design_b_locus5_5prime_race_coord_bp',
                        'obs_design_b_locus5_3prime_race_coord_bp',
                        'obs_design_b_locus5_race_umi_n_eff',
                        'obs_design_b_locus5_aberrant_umi_k',
                    )
                )
            ):
                _r5_26 = validated_stage0.get('obs_design_b_locus5_5prime_race_coord_bp')
                _r3_26 = validated_stage0.get('obs_design_b_locus5_3prime_race_coord_bp')
                _un_26 = validated_stage0.get('obs_design_b_locus5_race_umi_n_eff')
                _uk_26 = validated_stage0.get('obs_design_b_locus5_aberrant_umi_k')
                if _r5_26 is None or _r3_26 is None:
                    iip0026_hold_status = 'HOLD_MISSING_DESIGNB_LOCUS5_RACE_COORDINATES'
                elif int(round(_r5_26)) <= 0 or int(round(_r3_26)) <= 0 or int(round(_r5_26)) >= int(round(_r3_26)):
                    iip0026_hold_status = 'HOLD_INVALID_RACE_COORDINATES'
                elif _un_26 is None or int(round(_un_26)) < 3689 or _uk_26 is None:
                    iip0026_hold_status = 'HOLD_MISSING_DESIGNB_LOCUS5_TRANSCRIPT_UMI_3689'
                elif int(round(_uk_26)) > 0:
                    iip0026_hold_status = 'HOLD_DESIGNB_LOCUS5_CHIMERIC_TRANSCRIPT_DETECTED'
        elif _is_bsct_26 and (_req_iip0026_full or observations.get('b_sct_ciita_exon5_and_trimer_verified') is not None):
            if observations.get('b_sct_ciita_exon5_and_trimer_verified') is not True and iip0026_hold_status is None:
                iip0026_hold_status = 'HOLD_MISSING_BSCT_CIITA_EXON5_AND_TRIMER_EVIDENCE'
        if iip0026_hold_status is None and (
            _is_tyk2_pe_26
            and (
                _req_iip0026_full
                or observations.get('tyk2_pe_modality') is not None
                or observations.get('tyk2_phased_haplotype_sap_verified') is not None
                or any(
                    _tk in validated_stage0
                    for _tk in (
                        'obs_tyk2_pe4max_pentaplex_msi_shifts_k',
                        'obs_tyk2_pe4max_wgs_repeat_indels_clone',
                        'obs_tyk2_pe4max_wgs_repeat_indels_ctrl',
                        'obs_tyk2_pe4max_coding_or_splice_de_novo_indels',
                    )
                )
            )
        ):
            _pem_26 = observations.get('tyk2_pe_modality')
            _pem_str = _pem_26.strip() if isinstance(_pem_26, str) else None
            if _pem_str not in ALLOWED_TYK2_PE_MODALITIES or observations.get('tyk2_phased_haplotype_sap_verified') is not True:
                iip0026_hold_status = 'HOLD_MISSING_TYK2_PE_MODALITY_OR_PHASED_HAPLOTYPE'
            elif _pem_str == 'PE4MAX_MSI_GATED':
                _msi_26 = validated_stage0.get('obs_tyk2_pe4max_pentaplex_msi_shifts_k')
                _cc_26 = validated_stage0.get('obs_tyk2_pe4max_wgs_repeat_indels_clone')
                _ct_26 = validated_stage0.get('obs_tyk2_pe4max_wgs_repeat_indels_ctrl')
                _cd_26 = validated_stage0.get('obs_tyk2_pe4max_coding_or_splice_de_novo_indels')
                if _msi_26 is None:
                    iip0026_hold_status = 'HOLD_MISSING_MSI'
                elif int(round(_msi_26)) > 0:
                    iip0026_hold_status = 'HOLD_MSI_DETECTED_PE4MAX_REJECTED'
                elif _cc_26 is None or _ct_26 is None or _cd_26 is None:
                    iip0026_hold_status = 'HOLD_MISSING_PE4MAX_REPEAT_INDEL_EVIDENCE'
                elif int(round(_cd_26)) > 0:
                    iip0026_hold_status = 'HOLD_PE4MAX_CODING_OR_SPLICE_DE_NOVO_INDEL_DETECTED'
                else:
                    _ct_i = int(round(_ct_26))
                    _cap_26 = (_ct_i + 5) if _ct_i <= 25 else int(math.floor(1.2 * _ct_i))
                    if int(round(_cc_26)) > _cap_26:
                        iip0026_hold_status = 'HOLD_PE4MAX_REPEAT_INDEL_CAP_EXCEEDED'
    iip0026_eval_summary: dict[str, object] | None = None
    if _iip0026_active:
        iip0026_eval_summary = {
            'gate_pass': iip0026_hold_status is None,
            'hold_status': iip0026_hold_status,
            'endogenous_hlac_status': iip0026_hlac_status,
            'status': iip0026_hold_status if iip0026_hold_status is not None else 'CLEARED_IIP0026_MOLECULAR_PROVENANCE_AND_DESIGNB_GUARDRAILS',
        }

    # Stage-0 registry-bound arm-separated Gate 0d, paired intact-plus-increment CD16+ ADCC, measured-potency DEC-084 mixture, and visit/lot/sister-bound 72h/d28/d84 dual-modality IAPP gates (IIP-0028; DEC-106):
    _iip0028_gate0d_adcc_obs_keys = (
        'obs_arm1_unconditioned_cd4_si',
        'obs_arm1_unconditioned_annexin_v_pct',
        'obs_parallel_arm3_control_cd4_apoptosis_pct',
        'obs_post_depletion_input_normalized_cpep_ratio',
        'obs_fcr_cd16_intact_serum_total_lysis_pct',
        'obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct',
        'obs_adcc_isotype_background_lysis_pct',
        'obs_igg_depletion_efficiency_pct',
    )
    _iip0028_dec084_obs_keys = (
        'obs_dec084_starting_potency_ratio_out_b_over_out_a',
        'obs_dec084_ocr_ratio_b_over_a',
        'obs_dec084_cpep_auc_ratio_b_over_a',
        'obs_dec084_gcg_auc_ratio_b_over_a',
    )
    _iip0028_iapp_obs_keys = (
        'obs_iapp_72h_immunocapture_oligomer_elisa_fold',
        'obs_iapp_72h_colocalized_fibril_area_pct',
        'obs_iapp_d28_immunocapture_oligomer_elisa_fold',
        'obs_iapp_d28_colocalized_fibril_area_pct',
        'obs_iapp_d84_immunocapture_oligomer_elisa_fold',
        'obs_iapp_d84_colocalized_fibril_area_pct',
    )
    _iip0028_meta_or_bool_keys = (
        'require_full_iip0028_postdepletion_adcc_mixture_and_iapp_panel',
        'washout_and_natural_decay_unpooled_verified',
        'carrier_and_soluble_adjunct_below_lod_verified',
        'fcr_effector_fcgr3a_158v_typed_verified',
        'dec084_measured_starting_composition_and_fresh_equal_output_verified',
        'iapp_sister_visit_lot_and_viable_beta_mass_concordance_verified',
        'recipient_serum_pair_id',
        'intact_serum_sample_id',
        'igg_depleted_serum_sample_id',
        'positive_adcc_control_id',
        'isotype_adcc_control_id',
        'arm2_depleted_sample_id',
        'contemporaneous_arm3_control_sample_id',
        'sister_assembly_id',
    )
    _req_iip0028_full = bool(observations.get('require_full_iip0028_postdepletion_adcc_mixture_and_iapp_panel') is True)
    _iip0028_gate0d_adcc_active = bool(
        _req_iip0028_full
        or any(_k in validated_stage0 for _k in _iip0028_gate0d_adcc_obs_keys)
        or any(
            observations.get(_k) is not None
            for _k in (
                'washout_and_natural_decay_unpooled_verified',
                'carrier_and_soluble_adjunct_below_lod_verified',
                'fcr_effector_fcgr3a_158v_typed_verified',
                'recipient_serum_pair_id',
                'intact_serum_sample_id',
                'igg_depleted_serum_sample_id',
                'positive_adcc_control_id',
                'isotype_adcc_control_id',
                'arm2_depleted_sample_id',
                'contemporaneous_arm3_control_sample_id',
            )
        )
    )
    _iip0028_dec084_active = bool(
        any(_k in validated_stage0 for _k in _iip0028_dec084_obs_keys)
        or observations.get('dec084_measured_starting_composition_and_fresh_equal_output_verified') is not None
    )
    _iip0028_iapp_active = bool(
        _req_iip0028_full
        or any(_k in validated_stage0 for _k in _iip0028_iapp_obs_keys)
        or observations.get('iapp_sister_visit_lot_and_viable_beta_mass_concordance_verified') is not None
        or observations.get('sister_assembly_id') is not None
    )
    _iip0028_active = bool(
        _req_iip0028_full
        or _iip0028_gate0d_adcc_active
        or _iip0028_dec084_active
        or _iip0028_iapp_active
    )
    iip0028_hold_status: str | None = None
    iip0028_gate0d_adcc_status: str | None = None
    iip0028_dec084_status: str | None = None
    iip0028_iapp_status: str | None = None
    if _iip0028_active:
        def _explicit_se_ok_28(mk_28: str) -> bool:
            return mk_28 in validated_stage0 and se_provenance_by_mean.get(mk_28) == 'EXPLICIT_SE_SUPPLIED'

        def _unpadded_id_28(raw_id_28: object) -> bool:
            return isinstance(raw_id_28, str) and raw_id_28 == raw_id_28.strip() and len(raw_id_28) >= 4

        if product_def is None:
            iip0028_hold_status = 'HOLD_MISSING_PRODUCT_VARIANT_ID'
        elif str(product_def.get('product_variant_id', '')).startswith('ISLET1_B_SCT_'):
            iip0028_hold_status = 'HOLD_UNSUPPORTED_TRACK2_BSCT_REQUIRES_SEPARATE_GATE0B_SCT'

        if iip0028_hold_status is None and _iip0028_gate0d_adcc_active:
            if not isinstance(_adj_arm_raw, str) or _adj_arm_raw != _adj_arm_raw.strip() or _adj_arm_raw not in ALLOWED_ADJUNCT_CONDITIONING_HISTORY_ARMS:
                iip0028_gate0d_adcc_status = 'HOLD_GATE0D_MISSING_EXPLICIT_HISTORY_ARM'
            elif _adj_arm_raw == 'ARM1_NEVER_CONDITIONED':
                if not (_explicit_se_ok_28('obs_arm1_unconditioned_cd4_si') and _explicit_se_ok_28('obs_arm1_unconditioned_annexin_v_pct')):
                    iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM1_INVALID_OR_MISSING_CONTROL'
                else:
                    iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM1_NEVER_CONDITIONED_REFERENCE_CONTROL'
            else:
                _sp_id = observations.get('recipient_serum_pair_id')
                _in_sid = observations.get('intact_serum_sample_id')
                _dp_sid = observations.get('igg_depleted_serum_sample_id')
                _pos_cid = observations.get('positive_adcc_control_id')
                _iso_cid = observations.get('isotype_adcc_control_id')
                if not _unpadded_id_28(_sp_id):
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_MISSING_OR_INVALID_RECIPIENT_SERUM_PAIR_ID'
                elif not (_unpadded_id_28(_in_sid) and _unpadded_id_28(_dp_sid)) or _in_sid == _dp_sid:
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_MISSING_OR_DUPLICATE_PAIRED_SERUM_SAMPLE_IDS'
                elif not (_unpadded_id_28(_pos_cid) and _unpadded_id_28(_iso_cid)) or _pos_cid == _iso_cid:
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_MISSING_OR_DUPLICATE_POSITIVE_ISOTYPE_CONTROL_IDS'
                elif observations.get('cd16_adcc_igg_depletion_paired_control_verified') is not True:
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_MISSING_IGG_DEPLETED_MATCHED_SERUM_CONTROL'
                elif observations.get('fcr_effector_fcgr3a_158v_typed_verified') is not True:
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_UNVERIFIED_CD16_FCGR3A_EFFECTOR_COMPETENCE'
                elif not _explicit_se_ok_28('obs_igg_depletion_efficiency_pct') or validated_stage0['obs_igg_depletion_efficiency_pct'] < 95.0:
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_IGG_DEPLETION_EFFICIENCY_BELOW_95PCT'
                elif not (
                    _explicit_se_ok_28('obs_fcr_cd16_intact_serum_total_lysis_pct')
                    and _explicit_se_ok_28('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct')
                    and _explicit_se_ok_28('obs_cd16_adcc_positive_control_lysis_pct')
                    and _explicit_se_ok_28('obs_adcc_isotype_background_lysis_pct')
                ):
                    iip0028_gate0d_adcc_status = 'HOLD_ADCC_INCOMPLETE_OR_DEFAULT_SE_INTACT_AND_INCREMENT_PANEL'
                elif _upper95('obs_fcr_cd16_intact_serum_total_lysis_pct') > 5.0:
                    iip0028_gate0d_adcc_status = 'HOLD_EXCEEDED_UPPER95_obs_fcr_cd16_intact_serum_total_lysis_pct'
                elif _upper95('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct') > 5.0:
                    iip0028_gate0d_adcc_status = 'HOLD_EXCEEDED_UPPER95_obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct'
                elif _lower95('obs_cd16_adcc_positive_control_lysis_pct') < 25.0:
                    iip0028_gate0d_adcc_status = 'HOLD_BELOW_LOWER95_obs_cd16_adcc_positive_control_lysis_pct'
                elif _upper95('obs_adcc_isotype_background_lysis_pct') > 5.0:
                    iip0028_gate0d_adcc_status = 'HOLD_EXCEEDED_UPPER95_obs_adcc_isotype_background_lysis_pct'
                elif _adj_arm_raw == 'ARM3_CONTINUOUSLY_ACTIVE':
                    if observations.get('obs_arm3_ongoing_adjunct_bioactivity_verified') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM3_ONGOING_BIOACTIVITY_UNVERIFIED'
                    elif not all(
                        _explicit_se_ok_28(_mk3)
                        for _mk3 in (
                            'obs_gate0d_activated_cd4_fasl_apoptosis_pct',
                            'obs_indirect_apc_dq2_dq8_cd4_stimulation_index',
                            'obs_indirect_apc_cd4_ifng_si',
                            'obs_gate0d_transwell_beta_annexin_v_pct',
                        )
                    ):
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM3_INCOMPLETE_OR_DEFAULT_SE_PANEL'
                    elif (
                        _lower95('obs_gate0d_activated_cd4_fasl_apoptosis_pct') < 80.0
                        or _upper95('obs_indirect_apc_dq2_dq8_cd4_stimulation_index') > 1.5
                        or _upper95('obs_indirect_apc_cd4_ifng_si') > 2.0
                        or _upper95('obs_gate0d_transwell_beta_annexin_v_pct') > 15.0
                    ):
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM3_ACTIVE_ADJUNCT_BOUND_FAILED'
                    else:
                        iip0028_gate0d_adcc_status = 'CLEARED_ARM3_ACTIVE_ADJUNCT_AND_FCR_ADCC_PROFILE'
                else:  # ARM2_CONDITIONED_DEPLETED
                    _a2_sid = observations.get('arm2_depleted_sample_id')
                    _a3_cid = observations.get('contemporaneous_arm3_control_sample_id')
                    if adjunct_depletion_mode not in ('NATURAL_AGING_DECAY', 'PHYSICAL_WASHOUT'):
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_POOLED_OR_UNSPECIFIED_DEPLETION_MODE'
                    elif observations.get('washout_and_natural_decay_unpooled_verified') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_POOLED_WASHOUT_AND_AGING_ALIQUOTS'
                    elif observations.get('carrier_and_soluble_adjunct_below_lod_verified') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_RESIDUAL_ADJUNCT_ABOVE_LOD'
                    elif observations.get('obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_FASL_DEPLETION_UNVERIFIED'
                    elif observations.get('obs_arm2_sa_thbd_bioactivity_depleted_verified') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_SA_THBD_DEPLETION_UNVERIFIED'
                    elif observations.get('obs_arm2_contemporaneous_sham_handling_matched') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_CONTEMPORANEOUS_SHAM_UNMATCHED'
                    elif observations.get('obs_arm3_ongoing_adjunct_bioactivity_verified') is not True:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_PARALLEL_ARM3_BIOACTIVITY_CONTROL_UNVERIFIED'
                    elif not (_unpadded_id_28(_a2_sid) and _unpadded_id_28(_a3_cid)) or _a2_sid == _a3_cid:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_MISSING_OR_DUPLICATE_CONTEMPORANEOUS_ARM3_CONTROL_SAMPLE_ID'
                    elif adjunct_immune_responder_history not in ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES:
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_MISSING_IMMUNE_RESPONDER_HISTORY'
                    elif adjunct_sa_fhl1_status == 'RETAINED_ACTIVE':
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_PARTIALLY_DEPLETED_SA_FHL1_RETAINED_ACTIVE_NOT_COMPLETE_DEPLETION'
                    elif adjunct_sa_fhl1_status not in ('DEPLETED_VERIFIED', 'ABSENT_FROM_LOT'):
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_INVALID_OR_MISSING_SA_FHL1_STATUS'
                    elif not all(
                        _explicit_se_ok_28(_mk2)
                        for _mk2 in (
                            'obs_parallel_arm3_control_cd4_apoptosis_pct',
                            'obs_indirect_apc_dq2_dq8_cd4_stimulation_index',
                            'obs_indirect_apc_cd4_ifng_si',
                            'obs_gate0d_transwell_beta_annexin_v_pct',
                            'obs_arm2_post_depletion_nhs_c5b9_mac_pct',
                            'obs_arm2_post_depletion_viable_beta_recovery_pct',
                            'obs_arm2_post_depletion_perifusion_gsis_si',
                            'obs_post_depletion_input_normalized_cpep_ratio',
                        )
                    ):
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM2_INCOMPLETE_OR_DEFAULT_SE_POST_DEPLETION_PANEL'
                    elif (
                        _lower95('obs_parallel_arm3_control_cd4_apoptosis_pct') < 80.0
                        or _upper95('obs_indirect_apc_dq2_dq8_cd4_stimulation_index') > 1.5
                        or _upper95('obs_indirect_apc_cd4_ifng_si') > 2.0
                        or _upper95('obs_gate0d_transwell_beta_annexin_v_pct') > 15.0
                        or _upper95('obs_arm2_post_depletion_nhs_c5b9_mac_pct') > 5.0
                        or _lower95('obs_arm2_post_depletion_viable_beta_recovery_pct') < 85.0
                        or _lower95('obs_arm2_post_depletion_perifusion_gsis_si') < 2.5
                        or _lower95('obs_post_depletion_input_normalized_cpep_ratio') < 0.85
                    ):
                        iip0028_gate0d_adcc_status = 'HOLD_GATE0D_ARM2_POST_DEPLETION_STUDENT_T_BOUND_FAILED'
                    else:
                        iip0028_gate0d_adcc_status = 'CLEARED_ARM2_COMPLETE_DEPLETION_GATE0D_AND_FCR_ADCC_OVERLAY'
            if not str(iip0028_gate0d_adcc_status).startswith('CLEARED_'):
                iip0028_hold_status = iip0028_gate0d_adcc_status

        if iip0028_hold_status is None and _iip0028_dec084_active:
            if observations.get('dec084_measured_starting_composition_and_fresh_equal_output_verified') is not True:
                iip0028_dec084_status = 'HOLD_DEC084_UNVERIFIED_STARTING_COMPOSITION_OR_EQUAL_OUTPUT_DOSE_RESPONSE'
            elif not all(_explicit_se_ok_28(_dk) for _dk in _iip0028_dec084_obs_keys):
                iip0028_dec084_status = 'HOLD_DEC084_INCOMPLETE_OR_DEFAULT_SE_CONFIRMATION_PANEL'
            elif _upper95('obs_dec084_ocr_ratio_b_over_a') > 0.90:
                iip0028_dec084_status = 'HOLD_DEC084_OCR_RATIO_UPPER95_EXCEEDS_0P90'
            elif _lower95('obs_dec084_cpep_auc_ratio_b_over_a') < 0.85 or _lower95('obs_dec084_gcg_auc_ratio_b_over_a') < 0.85:
                iip0028_dec084_status = 'HOLD_DEC084_CPEP_OR_GLUCAGON_AUC_LOWER95_BELOW_0P85'
            else:
                iip0028_dec084_status = 'CLEARED_DEC084_MEASURED_POTENCY_OCR_AND_GLUCAGON_GATE'
            if not str(iip0028_dec084_status).startswith('CLEARED_'):
                iip0028_hold_status = iip0028_dec084_status

        if iip0028_hold_status is None and _iip0028_iapp_active:
            _sis_id_28 = observations.get('sister_assembly_id')
            _dev_id_28 = observations.get('device_assembly_id')
            _n_lots_28 = validated_stage0.get('obs_stage0_independent_lot_replicates_n')
            if observations.get('iapp_age_matched_reference_above_loq_verified') is not True or observations.get('iapp_sister_visit_lot_and_viable_beta_mass_concordance_verified') is not True:
                iip0028_iapp_status = 'HOLD_INCONCLUSIVE_IAPP_REFERENCE_LOQ_OR_VIABLE_BETA_CONCORDANCE_UNVERIFIED'
            elif _n_lots_28 is None or int(round(_n_lots_28)) < 6:
                iip0028_iapp_status = 'HOLD_IAPP_UNVERIFIED_INDEPENDENT_LOTS_LT_6'
            elif not _unpadded_id_28(_sis_id_28) or (_dev_id_28 is not None and _sis_id_28 == _dev_id_28):
                iip0028_iapp_status = 'HOLD_IAPP_D84_REQUIRES_AUTHORIZED_SISTER_ASSEMBLY'
            elif not all(_explicit_se_ok_28(_ik) for _ik in _iip0028_iapp_obs_keys):
                iip0028_iapp_status = 'HOLD_INCONCLUSIVE_IAPP_INCOMPLETE_OR_DEFAULT_SE_72H_D28_D84_PANEL'
            else:
                _tp_pairs_28 = [
                    (validated_stage0['obs_iapp_72h_immunocapture_oligomer_elisa_fold'], validated_stage0['obs_iapp_72h_colocalized_fibril_area_pct']),
                    (validated_stage0['obs_iapp_d28_immunocapture_oligomer_elisa_fold'], validated_stage0['obs_iapp_d28_colocalized_fibril_area_pct']),
                    (validated_stage0['obs_iapp_d84_immunocapture_oligomer_elisa_fold'], validated_stage0['obs_iapp_d84_colocalized_fibril_area_pct']),
                ]
                if len(set(_tp_pairs_28)) != 3:
                    iip0028_iapp_status = 'HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE'
                elif (
                    _lower95('obs_iapp_72h_immunocapture_oligomer_elisa_fold') > 1.5
                    or _lower95('obs_iapp_72h_colocalized_fibril_area_pct') > 1.0
                ):
                    iip0028_iapp_status = 'HOLD_FAIL_ACUTE_PROTEOSTASIS_72H'
                elif (
                    _upper95('obs_iapp_72h_immunocapture_oligomer_elisa_fold') > 1.5
                    or _upper95('obs_iapp_72h_colocalized_fibril_area_pct') > 1.0
                ):
                    iip0028_iapp_status = 'HOLD_INCONCLUSIVE_IAPP_EVIDENCE_72H_UPPER_BOUND'
                else:
                    _d28_e_l = _lower95('obs_iapp_d28_immunocapture_oligomer_elisa_fold') > 1.5
                    _d28_f_l = _lower95('obs_iapp_d28_colocalized_fibril_area_pct') > 1.0
                    _d84_e_l = _lower95('obs_iapp_d84_immunocapture_oligomer_elisa_fold') > 1.5
                    _d84_f_l = _lower95('obs_iapp_d84_colocalized_fibril_area_pct') > 1.0
                    if (_d28_e_l and _d28_f_l) or (_d84_e_l and _d84_f_l):
                        iip0028_iapp_status = 'HOLD_CONFIRMED_CHRONIC_IAPP_BREACH_P2_ITM2B_RESERVE_ELIGIBLE'
                    elif _d28_e_l or _d28_f_l or _d84_e_l or _d84_f_l:
                        iip0028_iapp_status = 'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION'
                    elif (
                        _upper95('obs_iapp_d28_immunocapture_oligomer_elisa_fold') > 1.5
                        or _upper95('obs_iapp_d28_colocalized_fibril_area_pct') > 1.0
                        or _upper95('obs_iapp_d84_immunocapture_oligomer_elisa_fold') > 1.5
                        or _upper95('obs_iapp_d84_colocalized_fibril_area_pct') > 1.0
                    ):
                        iip0028_iapp_status = 'HOLD_INCONCLUSIVE_IAPP_EVIDENCE'
                    else:
                        iip0028_iapp_status = 'CLEARED_DUAL_ORTHOGONAL_IAPP_72H_D28_D84'
            if not str(iip0028_iapp_status).startswith('CLEARED_'):
                iip0028_hold_status = iip0028_iapp_status

    iip0028_eval_summary: dict[str, object] | None = None
    if _iip0028_active:
        iip0028_eval_summary = {
            'gate_pass': iip0028_hold_status is None,
            'hold_status': iip0028_hold_status,
            'gate0d_and_adcc_status': iip0028_gate0d_adcc_status,
            'dec084_mixture_status': iip0028_dec084_status,
            'iapp_3way_status': iip0028_iapp_status,
            'p2_itm2b_active_weight_pct': 0.0,
            'status': (
                iip0028_hold_status
                if iip0028_hold_status is not None
                else 'CLEARED_IIP0028_POSTDEPLETION_ADCC_MIXTURE_AND_IAPP_GATES'
            ),
        }

    # Stage-0 resting (0 h) + interferon-stimulated (24 h) CXCL10 promoter-of-origin 5'-RACE TSS attribution & conditional BIO-06 AAV-ITR RNA-origin gate (IIP-0027; DEC-107):
    _iip0027_obs_keys = (
        'obs_cxcl10_resting_race_umi_n_eff',
        'obs_cxcl10_resting_unexpected_tss_umi_k',
        'obs_cxcl10_resting_unresolved_origin_umi_k',
        'obs_cxcl10_ifng24h_race_umi_n_eff',
        'obs_cxcl10_ifng24h_unexpected_tss_umi_k',
        'obs_cxcl10_ifng24h_unresolved_origin_umi_k',
        'obs_bio06_retained_aav_itr_dna_pct',
        'obs_bio06_itr_origin_rna_umi_k',
    )
    _iip0027_meta_or_bool_keys = (
        'require_full_iip0027_promoter_origin_and_itr_panel',
        'cxcl10_insert_matched_comparator_control_verified',
        'cxcl10_resting_expression_below_loq_confirmed',
        'bio06_itr_rna_origin_attribution_verified',
        'cxcl10_tss_donor_id',
        'cxcl10_tss_transcription_unit_id',
    )
    _iip0027_active = bool(
        any(_k in validated_stage0 for _k in _iip0027_obs_keys)
        or any(observations.get(_k) is not None for _k in _iip0027_meta_or_bool_keys)
    )
    iip0027_eval_summary: dict[str, object] | None = None
    iip0027_hold_status: str | None = None
    if _iip0027_active:
        iip0027_eval_summary = evaluate_cxcl10_promoter_origin_and_itr_gate_iip_0027(
            observations,
            validated_stage0=validated_stage0,
            non_integer_denominator_flagged=non_integer_scored_cell_or_trial_denominator,
        )
        if iip0027_eval_summary['gate_pass'] is not True:
            iip0027_hold_status = str(iip0027_eval_summary['status'])

    # Stage-0 ALK5i-free two-clock maturation, two-pulse Readout 2 with adult-islet reference floor,
    # inert-tracer tau_chamber, and >138.3 IU Gate 0S atomic escalation gate (IIP-0029; DEC-108):
    _iip0029_obs_keys = (
        'obs_stage7_total_maturation_days',
        'obs_stage7_days_post_sort_reaggregation',
        'obs_stage7_eta_sort_viable_recovery',
        'obs_readout2_inert_tracer_tau_chamber_min',
        'obs_readout2_pulse2_net_increment_z2',
        'obs_readout2_adult_islet_ref_z2',
        'obs_readout2_low2_frac_peak1',
        'obs_readout2_low2_over_low1_ratio',
        'obs_stage7_s8_fieller_z2_diff',
        'obs_stage7_stored_mature_insulin_iu_implanted',
        'obs_stage7_stored_mature_insulin_max_lot_iu',
        'obs_stage7_proinsulin_to_cpeptide_molar_ratio',
        'obs_stage7_pre_iec_zm447439_residual_uM',
    )
    _iip0029_meta_or_bool_keys = (
        'require_full_iip0029_twoclock_maturation_panel',
        'obs_alk5i_free_maturation_applied',
        'terminal_alk5i_tgfb_inhibitor_continued',
        'stage7_t3_10nM_and_nac_1mM_retained',
        'zm447439_washout_verified_pre_iec',
        'pre_coassembly_iec_vascular_function_verified',
        'readout2_both_40_and_100mmhg_tracer_aligned_verified',
        'post_sort_endocrine_zm447439_regimen',
        'stage7_start_utc',
        'sort_reaggregation_utc',
        'assay_harvest_utc',
        'registered_mixture_id',
        'locked_analysis_plan_id',
        'gate0s_elevated_store_scope_matched_device_lot_id',
        'iip0029_twoclock_record',
    )
    _iip0029_active = bool(
        any(_k in validated_stage0 for _k in _iip0029_obs_keys)
        or any(observations.get(_k) is not None for _k in _iip0029_meta_or_bool_keys)
    )
    iip0029_eval_summary: dict[str, object] | None = None
    iip0029_hold_status: str | None = None
    if _iip0029_active:
        iip0029_eval_summary = evaluate_stage0_twoclock_flat_or_structured_iip_0029(
            observations,
            validated_stage0=validated_stage0,
            se_provenance_by_mean=se_provenance_by_mean,
            upper95_fn=_upper95,
            lower95_fn=_lower95,
            root=root_path,
            product_def=product_def,
        )
        if iip0029_eval_summary['gate_pass'] is not True:
            iip0029_hold_status = str(iip0029_eval_summary['status'])

    _staged_l2_v3_ids = getattr(p_mod, 'STAGED_LOCUS2_V3_PRODUCT_IDS', frozenset())
    _has_explicit_pid_or_s11_dt = bool(
        product_def is not None
        or 'product_variant_id' in observations
        or 'iip0031_locus2_record' in observations
        or any(k in observations for k in _ALL_S11_RECOGNIZED_EVIDENCE_KEYS)
    )
    if _has_explicit_pid_or_s11_dt:
        iip0031_s11_eval = evaluate_locus2_ciita_modeb_challenger_v3_gate(
            observations,
            root=root_path,
            product_def=product_def,
            _trusted_p_mod=p_mod,
        )
    else:
        iip0031_s11_eval = {
            'gate_id': RULE_S11_GATE_ID_IIP0031,
            'gate_pass': True,
            'status': 'NOT_APPLICABLE_NON_LOCUS2_MODEB_VARIANT',
            'reasons': [],
            'product_variant_id': None,
            'is_staged_locus2_modeb_challenger': False,
        }
    iip0031_tier6b_eval: dict[str, Any] | None = None
    if 'tier6b_sorted_lot_observations' in observations and observations.get('tier6b_sorted_lot_observations') is not None:
        _l2b_dt = product_def.get('staged_locus2_modeb_donor_binding') if isinstance(product_def, dict) else None
        iip0031_tier6b_eval = evaluate_tier6b_per_fraction_thinning_iip0031(
            observations.get('tier6b_sorted_lot_observations'),
            target_scope={
                'product_variant_id': (product_def.get('product_variant_id') if isinstance(product_def, dict) else observations.get('product_variant_id')),
                'construct_sha256': observations.get('construct_sha256') or (_l2b_dt.get('locus2_donor_sha256') if isinstance(_l2b_dt, dict) else None),
                'background_id': observations.get('background_id') or (_l2b_dt.get('background_id') if isinstance(_l2b_dt, dict) else None),
                'donor_seedline_id': observations.get('donor_seedline_id', observations.get('seed_line_id')),
                'mcb_id': observations.get('mcb_id'),
                'lot_id': observations.get('lot_id', observations.get('current_lot_id')),
                'campaign_id': observations.get('campaign_id'),
            },
            root=root_path,
        )
    iip0031_hold_status: str | None = None
    if iip0031_tier6b_eval is not None and iip0031_tier6b_eval.get('hard_stop') is True:
        iip0031_hold_status = 'FAIL_TIER6B_OUTGROWTH_DETECTED'
    elif iip0031_s11_eval.get('gate_pass') is not True:
        iip0031_hold_status = str(iip0031_s11_eval.get('status'))
    elif iip0031_tier6b_eval is not None and iip0031_tier6b_eval.get('gate_pass') is not True:
        iip0031_hold_status = str(iip0031_tier6b_eval.get('status'))

    resolved_arch = dict(tree_doc.get('default_architecture_assignments', {}))
    if product_def is not None:
        _pvid_canon = str(product_def['product_variant_id'])
        resolved_arch['product_variant_id'] = _pvid_canon
        resolved_arch['product_execution_supported'] = bool(product_def.get('execution_supported', False))
        if _pvid_canon in _staged_l2_v3_ids:
            _l2b_ra = product_def.get('staged_locus2_modeb_donor_binding') or {}
            resolved_arch['locus2_construct'] = str(_l2b_ra.get('locus2_donor_id', _pvid_canon))
            if product_def.get('family') == 'B':
                resolved_arch['locus1_construct'] = f'{_pvid_canon}_HLAA_KO_HLAB_KO_ENDOGENOUS_HLAC'
        elif _pvid_canon.startswith('ISLET1_B_SCT_'):
            resolved_arch['locus1_construct'] = f'{_pvid_canon}_B2M_ARM1B_PLUS_CIITA_HLAC_SCT'
            resolved_arch['locus2_construct'] = f'{_pvid_canon}_CIITA_HLAC_SCT'
        elif _pvid_canon in ('ISLET1_B2_5LOCUS_BASE', 'ISLET1_B2_5LOCUS_TYK2_PE', 'ISLET1_B1_6LOCUS_TYK2_KO_INS'):
            resolved_arch['locus1_construct'] = f'{_pvid_canon}_HLAA_KO_HLAB_KO_ENDOGENOUS_HLAC'
        elif _pvid_canon == 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2':
            resolved_arch['locus1_construct'] = 'ISLET1_A4_X_MOL_RNLS_WT_ON_Z_140UM_ZERO_CAO2'
    if iip0024_eval_summary is not None:
        resolved_arch['stage0_atomic_provenance_and_bank_guardrails_iip_0024'] = iip0024_eval_summary['status']
        resolved_arch['stage0_atomic_provenance_and_bank_evaluation_iip_0024'] = iip0024_eval_summary
    else:
        resolved_arch['stage0_atomic_provenance_and_bank_guardrails_iip_0024'] = 'NOT_TESTED'
    if iip0025_eval_summary is not None:
        resolved_arch['hurdle4_spatial_oxygen_catalase_and_route_challengers_iip_0025'] = iip0025_eval_summary['status']
        resolved_arch['hurdle4_spatial_oxygen_catalase_and_route_evaluation_iip_0025'] = iip0025_eval_summary
    else:
        resolved_arch['hurdle4_spatial_oxygen_catalase_and_route_challengers_iip_0025'] = 'NOT_TESTED'
    if knockin_splice_challenger_arm != 'PARENTAL_REFERENCE' or _chall_any_input:
        resolved_arch['knockin_splice_challenger_arm'] = knockin_splice_challenger_arm
        if knockin_locus4_orientation is not None:
            resolved_arch['knockin_locus4_orientation'] = knockin_locus4_orientation
        if challenger_umi_eff > 0.0:
            resolved_arch['challenger_effective_umi_molecules'] = round(challenger_umi_eff, 3)
        resolved_arch['knockin_splice_challenger_qualification_status'] = (
            challenger_hold_status
            if challenger_hold_status is not None
            else 'CHALLENGER_SPLICE_AND_RECOVERY_ADJUSTED_UMI_PANEL_CLEARED'
        )
    if l3_multimir_gate_status is not None:
        resolved_arch['locus3_multimir_modeb_gate_status'] = l3_multimir_gate_status
        resolved_arch['locus3_multimir_modeb_evaluation'] = l3_multimir_eval_summary
    if iip0019_eval_summary is not None:
        resolved_arch['transport_3d_clumping_h2o2_and_design_z_evaluation'] = iip0019_eval_summary
        if iip0019_s6_status is not None:
            resolved_arch['s6_3d_lateral_spacing_gate_status'] = iip0019_s6_status
        if iip0019_p05_status is not None:
            resolved_arch['p05_microbead_catalase_barrier_gate_status'] = iip0019_p05_status
        if iip0019_dz_status is not None:
            resolved_arch['design_z_two_stage_optode_gate_status'] = iip0019_dz_status
    if iip0020_eval_summary is not None:
        resolved_arch['biomaterial_antigen_humoral_and_human_cuff_evaluation'] = iip0020_eval_summary
        if iip0020_p05_species_status is not None:
            resolved_arch['p05_catalase_species_gate_status'] = iip0020_p05_species_status
        if iip0020_cuff_status is not None:
            resolved_arch['cuff_humoral_and_hepatocyte_gate_status'] = iip0020_cuff_status
    if iip0021_eval_summary is not None:
        resolved_arch['metabolic_rest_ocr_and_conditioning_evaluation'] = iip0021_eval_summary
        if iip0021_rest_status is not None:
            resolved_arch['rest_120h_ocr_coupon_status'] = iip0021_rest_status
        if iip0021_cond_status is not None:
            resolved_arch['conditioning_4arm_nqo1_gate_status'] = iip0021_cond_status
    if iip0022_eval_summary is not None:
        resolved_arch['track2_nk_split_and_pvr_necroptosis_evaluation'] = iip0022_eval_summary
        if iip0022_track2_status is not None:
            resolved_arch['track2_split_route_status'] = iip0022_track2_status
    if iip0023_eval_summary is not None:
        resolved_arch['design_y_wildcard_evaluation'] = iip0023_eval_summary
        resolved_arch['design_y_wildcard_status'] = iip0023_dy_status
    if (
        _adj_arm_raw is not None
        or adjunct_depletion_mode is not None
        or adjunct_sa_fhl1_status is not None
        or adjunct_immune_responder_history is not None
        or any(
            observations.get(_k) is not None
            for _k in (
                'obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified',
                'obs_arm2_sa_thbd_bioactivity_depleted_verified',
                'obs_arm2_contemporaneous_sham_handling_matched',
                'obs_arm3_ongoing_adjunct_bioactivity_verified',
            )
        )
        or any(
            _k in validated_stage0
            for _k in (
                'obs_arm2_post_depletion_viable_beta_recovery_pct',
                'obs_arm2_post_depletion_perifusion_gsis_si',
                'obs_arm2_post_depletion_nhs_c5b9_mac_pct',
            )
        )
    ):
        resolved_arch['adjunct_conditioning_history_arm'] = adjunct_conditioning_history_arm
        if adjunct_depletion_mode is not None:
            resolved_arch['adjunct_depletion_mode'] = adjunct_depletion_mode
        if adjunct_immune_responder_history is not None:
            resolved_arch['adjunct_immune_responder_history'] = adjunct_immune_responder_history
        if adjunct_sa_fhl1_status is not None:
            resolved_arch['adjunct_sa_fhl1_status'] = adjunct_sa_fhl1_status
            resolved_arch['adjunct_depletion_truth_in_labeling_status'] = (
                'PARTIAL_FASL_THBD_DEPLETED_SA_FHL1_RETAINED_ACTIVE_NOT_COMPLETE_ADJUNCT_DEPLETION'
                if adjunct_sa_fhl1_status == 'RETAINED_ACTIVE'
                else (
                    'COMPLETE_DECLARED_ADJUNCT_BIOACTIVITY_DEPLETED'
                    if (
                        observations.get('obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified') is True
                        and observations.get('obs_arm2_sa_thbd_bioactivity_depleted_verified') is True
                    )
                    else 'UNVERIFIED_ADJUNCT_DEPLETION_STATE'
                )
            )
        if adjunct_conditioning_history_arm == 'ARM2_CONDITIONED_DEPLETED':
            resolved_arch['arm2_post_adjunct_depletion_qualification_status'] = (
                'HOLD_GATE0D_ARM2_ADJUNCT_DEPLETION_OR_SHAM_VERIFICATION_INCOMPLETE'
                if arm2_depletion_verification_incomplete
                else (
                    'HOLD_GATE0D_ARM2_POST_ADJUNCT_DEPLETION_REBOUND_INJURY_OR_FUNCTION_LOSS'
                    if arm2_rebound_injury_or_function_failed
                    else (
                        'HOLD_GATE0D_ARM2_POST_ADJUNCT_DEPLETION_PANEL_INCOMPLETE'
                        if arm2_post_depletion_panel_incomplete
                        else 'ARM2_POST_ADJUNCT_DEPLETION_RECHALLENGE_CLEARED'
                    )
                )
            )
    if prior_s9_reset_history_retained:
        resolved_arch['prior_s9_safety_stop_history_retained'] = True
    if 'obs_dxs_screened_cell_count' in validated_stage0:
        resolved_arch['dxs_effective_clonogenic_trials'] = dxs_effective_clonogenic_n
        resolved_arch['dxs_denominator_scope'] = (
            'EXPLICIT_CLONOGENIC_AND_SPIKEIN_ADJUSTED'
            if has_explicit_dxs_clonogenic_factors
            else 'MISSING_EXPLICIT_CLONOGENIC_FRACTION_OR_SPIKEIN_DETECTION_EFFICIENCY'
        )
    if s1_promoted:
        resolved_arch['locus3_construct'] = 'ISLET1_A4_L3_BFB_BASALFLOOR'
    if s3_promoted:
        resolved_arch['locus1_construct'] = 'ISLET1_A4_DXS_ARM1C_L4CD47'
    elif catastrophic_nkg2a_lysis_exceeds_arm1x:
        resolved_arch['locus1_construct'] = 'HOLD_CATASTROPHIC_NKG2A_LYSIS_EXCEEDS_ARM1X_RESCUE_ENVELOPE'
    elif kir_nkg2a_neg_nk_failed:
        resolved_arch['locus1_construct'] = 'HOLD_KIR_NKG2A_NEG_NK_OR_IL12_IL18_SIRPA_LOW_LYSIS_FAILED_REQUIRES_TRACK2_B_SCT_OR_HLAG1_CHALLENGER'
    elif track2_hlag1_route_attempted:
        resolved_arch['locus1_construct'] = 'HOLD_TRACK2_HLAG1_CHALLENGER_PENDING_PRODUCT_REGISTRATION'
        resolved_arch['locus2_construct'] = 'HOLD_TRACK2_HLAG1_CHALLENGER_PENDING_PRODUCT_REGISTRATION'
        resolved_arch['track2_kir_subbank_construct'] = 'HOLD_HLAG1_C42S_LILRB1_CHALLENGER_UNREGISTERED_PRODUCT'
        resolved_arch['track1_mhc_null_kir_escape_status'] = 'HOLD_HLAG1_LILRB1_ROUTE_PENDING_PRODUCT_REGISTRATION'
    elif track2_kir_or_hlag1_fallback_promoted:
        resolved_arch['locus1_construct'] = 'PROMOTE_TRACK2_B_SCT_KIR_RESCUE' if (_bsct_flag or _t2_route == 'TRACK2_HLAC_B_SCT_KIR_RESCUE') else 'PROMOTE_TRACK2_B_SCT_OR_HLAG1_LILRB1_RESCUE'
        resolved_arch['locus2_construct'] = 'ISLET1_B_SCT_C1_CIITA_OR_C2_CIITA_KIR_RESCUE'
        resolved_arch['track2_kir_subbank_construct'] = 'ISLET1_B_SCT_C1_OR_C2_CIITA_KIR_RESCUE' if (_bsct_flag or _t2_route == 'TRACK2_HLAC_B_SCT_KIR_RESCUE') else 'ISLET1_B_SCT_C1_OR_C2_CIITA_OR_HLAG1_C42S_CHALLENGER'
        resolved_arch['track1_mhc_null_kir_escape_status'] = 'PROMOTED_TO_TRACK2_HLAC_B_SCT_KIR_RESCUE' if (_bsct_flag or _t2_route == 'TRACK2_HLAC_B_SCT_KIR_RESCUE') else 'PROMOTED_TO_TRACK2_B_SCT_OR_HLAG1_LILRB1_CHALLENGER'
    elif s2_promoted:
        resolved_arch['locus1_construct'] = 'ISLET1_A4_X_MOL_RNLS_WT'
        if observations.get('s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior') is not True:
            resolved_arch['locus1_arm1x_qualification_status'] = 'CANDIDATE_PROMOTED_REQUIRES_ARM1X_GATE0A_QUALIFICATION'
    elif cd47_or_phagocytosis_failed and not dxs_joint_retention_verified:
        resolved_arch['locus1_construct'] = 'HOLD_LOCUS1_CD47_DEFICIENT_PENDING_DXS_QUALIFICATION'
    elif unmitigated_cd47_or_hypoxic_retention_cliff_failed:
        resolved_arch['locus1_construct'] = 'HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED'
    if s10_v3_eval is not None:
        resolved_arch['s10_v3_functional_selectivity_evaluation'] = s10_v3_eval
    if s10_nkg2c_triggered or (_s10_v3_active and s10_v3_eval is not None):
        if _s10_v3_active and s10_v3_eval is not None:
            if not s10_v3_eval['gate_pass']:
                s10_status_val = str(s10_v3_eval['status'])
            elif not s10_v3_release_qualified:
                s10_status_val = 'CLEARED_S10_FUNCTIONAL_SELECTIVITY_SCREENING_15K_V3'
            elif s10_kd_ratio_lower95_sub4 or float(s10_v3_eval['recorded_sct_kc_over_ka_mean']) < 4.0:
                s10_status_val = 'CLEARED_S10_FUNCTIONAL_SELECTIVITY_SUB4_SCT_V3'
            else:
                s10_status_val = 'CLEARED_S10_STATIC_RATIO_AND_PAIRED_FUNCTIONAL_GATE_V3'
        else:
            s10_status_val = (
                'CLEARED_SCT_SELECTIVITY_COUPON_KC_KA_GE_4'
                if s10_coupon_cleared
                else 'HOLD_S10_NKG2C_SELECTIVITY_FAILED_REQUIRES_SCT_PEPTIDE_SELECTIVITY_COUPON_OR_CMV_STRATIFIED_LABEL'
            )
        resolved_arch['nkg2c_selectivity_contingency_status'] = s10_status_val
        resolved_arch['nkg2c_adaptive_nk_selectivity_status'] = s10_status_val
    if s8_tyk2_epistasis_triggered:
        if not s3_promoted and resolved_arch.get('locus1_construct') != 'HOLD_LOCUS1_CD47_DEFICIENT_PENDING_DXS_QUALIFICATION':
            resolved_arch['locus1_construct'] = 'ISLET1_A4_X_MOL_RNLS_WT'
        resolved_arch['locus3_construct'] = 'ISLET1_A4_L3_BFB_BASALFLOOR'
        resolved_arch['tyk2_epistasis_contingency_status'] = 'COUPLE_TYK2_P1104A_STRICTLY_WITH_ARM1X_AND_L3_BFB_CONSTITUTIVE_BASAL_FLOOR'
    elif (x14_trig or x34_trig or sat_trig) and not is_tyk2_product_or_unspecified:
        resolved_arch['tyk2_epistasis_contingency_status'] = 'INACTIVE_NON_TYK2_PRODUCT_VARIANT'
    if s4_promoted and 'obs_lot_ocr_ratio' in validated_stage0 and lot_ocr_upper95 > 1.097:
        resolved_arch['lot_payload_derating_factor'] = round(1.0 / validated_stage0['obs_lot_ocr_ratio'], 4)
        resolved_arch['pouch_areal_loading_cap_ieq_cm2'] = 1500.0
    if s8_promoted or prolonged_insulin_clearance_failed or assembled_cluster_lineage_or_identity_failed:
        resolved_arch['prescribed_functional_dose_cap_ieq'] = 150000.0
        resolved_arch['sentinel_staged_cassette_dosing_status'] = 'STAGE1_SINGLE_100CM2_SENTINEL_CASSETTE_150K_IEQ_CAP_BEFORE_STEPUP'
    if s6_promoted:
        resolved_arch['er_rescue_genotype_pivot'] = (
            'PROMOTE_ISLET1_A4_K6_HLAA_HLAB_KO_EXOGENOUS_B2M_SHIELD'
            if exogenous_b2m_class1_rescue_failed
            else 'REQUALIFY_ISLET1_A4_ARM1B_TYK2_PE_OR_K6_HLAA_HLAB_KO'
        )
    if s7_promoted:
        resolved_arch['track2_kir_subbank_construct'] = 'ISLET1_B_SCT_C1_OR_C2_CIITA'
    elif b_prime_sct_direct_cd8_failed:
        resolved_arch['track2_kir_subbank_construct'] = 'HOLD_B_PRIME_SCT_DIRECT_CD8_ALLOREACTIVITY_FAILED'
    elif s7_hold_unqualified_sct:
        resolved_arch['track2_kir_subbank_construct'] = 'HOLD_PENDING_B_SCT_GATE0B_DENSITY_QUALIFICATION'
    if b2_tyk2_gate_status is not None:
        resolved_arch['design_b2_tyk2_p1104a_hlac_gate_status'] = b2_tyk2_gate_status
        resolved_arch['design_b2_tyk2_p1104a_hlac_evaluation'] = b2_tyk2_eval_summary
        if _b2_tyk2_stratum is not None:
            resolved_arch['b2_tyk2_recipient_genotype_stratum'] = _b2_tyk2_stratum
        if _b2_tyk2_allele is not None:
            resolved_arch['b2_tyk2_hlac_allele'] = _b2_tyk2_allele

    # Terminal NO-GO and unmitigated fallback HOLD states (AUDIT-REF, S2-01..03, AUDIT-REF, AUDIT-REF, AUDIT-REF, MW-R5-S2-F01..F06, AUDIT-REF & AUDIT-REF; DEC-041..050):
    if prior_s9_flag:
        resolved_arch['reset_current_batch_id_explicitly_supplied'] = bool(cur_batch_id_str)
    resolved_arch['tp53_denominator_provenance'] = (
        'EXPLICIT_OBSERVED_DEPTH'
        if tp53_duplex_depth_explicitly_supplied
        else 'ASSUMED_PROTOCOL_MINIMUM_62958_PENDING_EXPLICIT_TP53_DENOMINATOR'
    )
    resolved_arch['castseq_denominator_provenance'] = (
        'EXPLICIT_OBSERVED_DEPTH'
        if cast_duplex_depth_explicitly_supplied
        else 'ASSUMED_PROTOCOL_MINIMUM_62958_PENDING_EXPLICIT_CAST_DENOMINATOR'
    )
    resolved_arch['mcb_clonogenic_escape_qualification_status'] = (
        'HARD_SAFETY_STOP_MCB_CLONOGENIC_SURVIVOR_DETECTED'
        if mcb_clonogenic_survivor_detected
        else (
            'CLEARED_MCB_CLONOGENIC_EFFECTIVE_TRIALS_GE_299573228'
            if (
                mcb_clonogenic_trials_explicitly_supplied
                and mcb_clonogenic_trials >= 299573228.0
                and mcb_survivors_k_explicitly_supplied
                and mcb_survivors_k == 0.0
            )
            else (
                'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING'
                if mcb_clonogenic_outcome_missing_when_trials_supplied
                else (
                    'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM'
                    if (mcb_clonogenic_trials_underpowered_below_299573228 or mcb_clonogenic_trials_missing_when_outcome_zero_supplied)
                    else 'BULK_4LOG_KILL_ONLY_REQUIRES_SEPARATE_MCB_3E8_CLONOGENIC_QUALIFICATION'
                )
            )
        )
    )
    avascular_outer_band_dropout = bool(
        'obs_outer_band_kperf_s_inv' in validated_stage0
        and _lower95('obs_outer_band_kperf_s_inv') <= 1e-5
    )
    omental_revascularization_cleared = bool(
        observations.get('s9_omental_revascularization_qualified') is True
        or (
            observations.get('s9_chronic_80um_lowload_or_omental_scaffold_qualified') is True
            and observations.get('s9_chronic_80um_lowload_qualified') is not True
        )
    )
    chronic_30_to_50_cleared = bool(
        omental_revascularization_cleared
        or observations.get('s9_chronic_80um_lowload_qualified') is True
        or observations.get('s9_chronic_80um_lowload_or_omental_scaffold_qualified') is True
    )
    if (has_chronic_capsule_obs and capsule_dz_upper95 > 50.0) or avascular_outer_band_dropout:
        resolved_arch['chronic_capsule_50um_or_avascular_escalation_status'] = (
            'CLEARED_BY_OMENTAL_REVASCULARIZATION'
            if omental_revascularization_cleared
            else 'REQUIRES_OMENTAL_REVASCULARIZATION_80UM_LOWLOAD_INSUFFICIENT_ABOVE_50UM_OR_AVASCULAR_DROPOUT'
        )
    if lot_ocr_exceeds_max_1p50x_ceiling:
        resolved_arch['lot_maturation_status'] = 'NO_GO_LOT_OCR_EXCEEDS_1P50X_MAX_POUCH_ENVELOPE'
    elif has_chronic_capsule_obs and capsule_dz_upper95 > 100.0:
        resolved_arch['lot_maturation_status'] = 'NO_GO_CHRONIC_CAPSULE_EXCEEDS_100UM_VIABILITY_ENVELOPE'
    elif lysis_2h_upper95 > 25.0:
        resolved_arch['lot_maturation_status'] = 'NO_GO_ACUTE_REGIME_II_LYSIS_EXCEEDS_25PCT_CEILING'
    elif tp53_duplex_depth_underpowered_below_62958:
        resolved_arch['lot_maturation_status'] = 'HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
    elif cast_duplex_depth_underpowered_below_62958:
        resolved_arch['lot_maturation_status'] = 'HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
    elif mcb_clonogenic_outcome_missing_when_trials_supplied:
        resolved_arch['lot_maturation_status'] = 'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING'
    elif mcb_clonogenic_trials_underpowered_below_299573228 or mcb_clonogenic_trials_missing_when_outcome_zero_supplied:
        resolved_arch['lot_maturation_status'] = 'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM'
    elif residual_psc_trials_underpowered:
        resolved_arch['lot_maturation_status'] = 'HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
    elif hazard_prolif_trials_underpowered:
        resolved_arch['lot_maturation_status'] = 'HOLD_HAZARD_STATE_PROLIF_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
    elif ap1903_step1_trials_underpowered:
        resolved_arch['lot_maturation_status'] = 'HOLD_AP1903_STEP1_EFFECTIVE_TRIALS_BELOW_29957323_MINIMUM'
    elif non_integer_scored_cell_or_trial_denominator:
        resolved_arch['lot_maturation_status'] = 'HOLD_NON_INTEGER_SCORED_CELL_OR_TRIAL_DENOMINATOR'
    elif mosaic_explicit_se_missing_scored_cell_denominator:
        resolved_arch['lot_maturation_status'] = 'HOLD_SINGLE_CELL_MOSAIC_REQUIRES_INTEGER_CELLS_SCORED_DENOMINATOR'
    elif zero_event_proportion_missing_denominator_or_explicit_se:
        resolved_arch['lot_maturation_status'] = 'HOLD_ZERO_EVENT_PROPORTION_GATE_REQUIRES_EXPLICIT_SE_OR_SCORED_CELL_DENOMINATOR'
    elif ap1903_insulin_dump_failed and observations.get('s9_staged_cassette_excision_and_d10w_clamp_qualified') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0S_AP1903_INSULIN_DUMP_EXCEEDS_15PCT_REQUIRES_STAGED_CASSETTE_EXCISION_AND_CLAMP'
    elif periimplant_fluid_film_gt_80_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_PERIIMPLANT_FLUID_FILM_GT_100UM_SEROMA_HYPOXIA_BREACH'
    elif periimplant_fluid_film_50_to_80_unmitigated:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_PERIIMPLANT_FLUID_FILM_GT_50UM_REQUIRES_CLOSED_SUCTION_DRAIN_AND_LOWLOAD_1500'
    elif explant_retrieval_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_EXPLANT_EN_BLOC_RECOVERY_BELOW_99PCT'
    elif senescent_p16_sasp_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0S_SENESCENT_P16_SASP_FRACTION_EXCEEDS_1PCT_CEILING'
    elif gate0v_viral_sanctuary_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0V_VIRAL_SANCTUARY_PERSISTENCE_FAILED'
    elif assembled_cassette_chandler_ibmir_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_ASSEMBLED_CASSETTE_CHANDLER_LOOP_IBMIR_OR_COMPLEMENT_FAILED'
    elif effero_clearance_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_APOPTOTIC_EFFEROCYTOSIS_CLEARANCE_IMPAIRED_DAMPS_RISK'
    elif assembled_cluster_lineage_or_identity_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_ASSEMBLED_CLUSTER_OFF_TARGET_LINEAGE_OR_BETA_IDENTITY_DRIFT_FAILED'
    elif indirect_apc_cd4_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0D_INDIRECT_APC_DQ2_DQ8_CD4_AUTOIMMUNE_STIMULATION_EXCEEDED'
    elif arm3_ongoing_bioactivity_explicit_false:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0D_ARM3_ONGOING_ADJUNCT_BIOACTIVITY_UNVERIFIED'
    elif arm1_never_conditioned_reference_hold:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0D_ARM1_NEVER_CONDITIONED_REFERENCE_CONTROL_NOT_RELEASE_ARM'
    elif arm2_depletion_verification_incomplete:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0D_ARM2_ADJUNCT_DEPLETION_OR_SHAM_VERIFICATION_INCOMPLETE'
    elif arm2_rebound_injury_or_function_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0D_ARM2_POST_ADJUNCT_DEPLETION_REBOUND_INJURY_OR_FUNCTION_LOSS'
    elif arm2_post_depletion_panel_incomplete:
        resolved_arch['lot_maturation_status'] = 'HOLD_GATE0D_ARM2_POST_ADJUNCT_DEPLETION_PANEL_INCOMPLETE'
    elif l3_multimir_hold_status is not None:
        resolved_arch['lot_maturation_status'] = l3_multimir_hold_status
    elif iip0019_transport_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0019_transport_hold_status
    elif iip0020_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0020_hold_status
    elif iip0021_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0021_hold_status
    elif iip0022_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0022_hold_status
    elif iip0023_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0023_hold_status
    elif iip0024_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0024_hold_status
    elif iip0025_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0025_hold_status
    elif iip0026_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0026_hold_status
    elif iip0028_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0028_hold_status
    elif iip0027_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0027_hold_status
    elif iip0029_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0029_hold_status
    elif iip0031_hold_status is not None:
        resolved_arch['lot_maturation_status'] = iip0031_hold_status
    elif challenger_hold_status is not None:
        resolved_arch['lot_maturation_status'] = challenger_hold_status
    elif avascular_outer_band_dropout and not omental_revascularization_cleared:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_CHRONIC_CAPSULE_GT_50UM_OR_AVASCULAR_DROPOUT_REQUIRES_OMENTAL_REVASCULARIZATION'
    elif subthreshold_outer_band_kperf and not omental_revascularization_cleared:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_SUBTHRESHOLD_OUTER_BAND_PERFUSION_BELOW_0P10_S_INV_REQUIRES_OMENTAL_REVASCULARIZATION'
    elif has_chronic_capsule_obs and capsule_dz_upper95 > 50.0 and observations.get('s9_chronic_80um_lowload_qualified') is True and not omental_revascularization_cleared:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_CHRONIC_CAPSULE_GT_50UM_OR_AVASCULAR_DROPOUT_REQUIRES_OMENTAL_REVASCULARIZATION'
    elif has_chronic_capsule_obs and capsule_dz_upper95 > 30.0 and not chronic_30_to_50_cleared:
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_CHRONIC_CAPSULE_GT_30UM_REQUIRES_80UM_LOWLOAD_OR_OMENTAL_QUALIFICATION'
    elif prolonged_insulin_clearance_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_S8_PROLONGED_INSULIN_CLEARANCE_OR_SUBCUTANEOUS_DRAINAGE_REQUIRES_DOSE_STAGING_AND_OMENTAL_PORTAL_DRAINAGE'
    elif s8_promoted:
        resolved_arch['lot_maturation_status'] = 'HOLD_FOR_PLUS_3W_EXTENDED_STAGE7_MATURATION'
    elif primary_biological_or_pouch_gate_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_FAILED_PRIMARY_BIOLOGICAL_OR_POUCH_GATE'
    elif dxs_hold_missing_explicit_se_or_counts:
        resolved_arch['lot_maturation_status'] = 'HOLD_DXS_RARE_EVENT_REQUIRES_EXPLICIT_SE_OR_CELL_COUNT'
    elif iec_vascular_support_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_PAN_LINEAGE_CD47_QPCTL_OR_VASCULAR_SHIELD_FAILED'
    elif catastrophic_nkg2a_lysis_exceeds_arm1x:
        resolved_arch['lot_maturation_status'] = 'HOLD_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_EXCEEDS_ARM1X_RESCUE_ENVELOPE'
    elif kir_nkg2a_neg_nk_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_KIR_NKG2A_NEG_NK_OR_IL12_IL18_SIRPA_LOW_LYSIS_EXCEEDED_REQUIRES_TRACK2_B_SCT_OR_HLAG1_CHALLENGER'
    elif track2_kir_or_hlag1_fallback_promoted and product_def is not None and str(product_def.get('family', 'A')).upper().startswith('A'):
        resolved_arch['lot_maturation_status'] = 'HOLD_TRACK1_MHC_NULL_KIR_ESCAPE_REQUIRES_TRACK2_B_SCT_OR_HLAG1_PRODUCT_SWITCH'
    elif track2_kir_or_hlag1_fallback_promoted and observations.get('track2_b_sct_or_hlag1_product_switched') is not True:
        resolved_arch['lot_maturation_status'] = 'CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION'
    elif receptor_neg_nkg2c_nk_failed and observations.get('s1_l3_bfb_replacement_qualified') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_S10_RECEPTOR_NEGATIVE_NKG2C_NK_LYSIS_EXCEEDED_REQUIRES_L3_BFB_AND_FASL_CUFF_CLEARANCE'
    elif s10_v3_hold_status is not None:
        resolved_arch['lot_maturation_status'] = s10_v3_hold_status
    elif s10_nkg2c_triggered and not s10_coupon_cleared:
        resolved_arch['lot_maturation_status'] = 'HOLD_S10_NKG2C_SELECTIVITY_FAILED_REQUIRES_SCT_PEPTIDE_SELECTIVITY_COUPON_OR_CMV_STRATIFIED_LABEL'
    elif excessive_immune_donor_variance_or_upper95_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_EXCESSIVE_DONOR_VARIANCE_OR_95PCT_UPPER_CONFIDENCE_BOUND_EXCEEDED'
    elif cd47_or_phagocytosis_failed and not dxs_joint_retention_verified:
        resolved_arch['lot_maturation_status'] = 'HOLD_LOCUS1_CD47_DEFICIENT_PENDING_DXS_QUALIFICATION'
    elif unmitigated_cd47_or_hypoxic_retention_cliff_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED'
    elif post_stress_nkg2d_dnam1_lysis_failed and observations.get('s1_l3_bfb_replacement_qualified') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_POST_STRESS_NKG2D_DNAM1_ATTRIBUTED_LYSIS_EXCEEDS_5PCT_REQUIRES_L3_BFB_AND_FASL_CUFF'
    elif b_prime_sct_direct_cd8_failed:
        resolved_arch['lot_maturation_status'] = 'HOLD_B_PRIME_ALLO_FAILED_REQUIRES_4DIGIT_MATCH_OR_SCT_CD8_MLC_CLEARANCE'
    elif s7_hold_unqualified_sct:
        resolved_arch['lot_maturation_status'] = 'HOLD_PENDING_B_SCT_GATE0B_DENSITY_QUALIFICATION'
    elif s7_hold_unverified_sct_allo:
        resolved_arch['lot_maturation_status'] = 'HOLD_B_PRIME_ALLO_FAILED_REQUIRES_4DIGIT_MATCH_OR_SCT_CD8_MLC_CLEARANCE'
    elif b2_tyk2_hold_status is not None:
        resolved_arch['lot_maturation_status'] = b2_tyk2_hold_status
    elif s8_tyk2_epistasis_underpowered_hold:
        resolved_arch['lot_maturation_status'] = 'HOLD_INCONCLUSIVE_UNDERPOWERED_TYK2_EPISTASIS_SE'
    elif s8_tyk2_epistasis_triggered and observations.get('s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_PENDING_COUPLED_ARM1X_L3_BFB_TYK2_QUALIFICATION'
    elif s6_promoted and observations.get('s6_er_rescue_genotype_requalified') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_S6_ER_UPR_THRESHOLD_BREACHED_REQUIRES_GENOTYPE_REQUALIFICATION'
    elif s5_promoted and observations.get('s5_sa_thbd_fhl1_only_replacement_qualified') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_S5_FASL_BYSTANDER_TOXICITY_REQUIRES_THBD_FHL1_ONLY_CUFF_QUALIFICATION'
    elif s1_promoted and observations.get('s1_l3_bfb_replacement_qualified') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_PENDING_L3_BFB_LOCUS3_REPLACEMENT_QUALIFICATION'
    elif s2_arm1b_failed and observations.get('s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_PENDING_ARM1X_LOCUS1_REPLACEMENT_QUALIFICATION'
    elif s4_promoted and observations.get('s4_rnls_wt_clears_gsis_shutoff_and_ocr_and_er_gate') is not True:
        resolved_arch['lot_maturation_status'] = 'HOLD_S4_OCR_EXCEEDS_1P097X_REQUIRES_LOWLOAD_1500_AND_RNLS_WT_QUALIFICATION'

    # Product-to-evidence binding for Track-2 B'-SCT and composite COMP_R4 targets when product_def is supplied (AUDIT-REF; DEC-050):
    is_track2_sct_product = bool(
        product_def is not None
        and (
            (
                str(product_def.get('family', 'A')).upper().startswith('B')
                and str(product_def.get('product_variant_id', '')) not in ('ISLET1_B2_5LOCUS_TYK2_PE', 'ISLET1_B2_5LOCUS_BASE')
            )
            or 'SCT' in str(product_def.get('product_variant_id', '')).upper()
        )
    )
    is_comp_r4_product = bool(
        product_def is not None
        and 'COMP_R4' in str(product_def.get('product_variant_id', '')).upper()
    )
    if (
        is_track2_sct_product
        and (
            not sct_density_explicitly_qualified
            or se_provenance_by_mean.get('obs_gate0b_hlac_sct_molecules_per_cell') != 'EXPLICIT_SE_SUPPLIED'
            or not sct_cd8_mlc_explicitly_cleared
        )
        and resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY')
    ):
        resolved_arch['lot_maturation_status'] = 'HOLD_TRACK2_B_SCT_REQUIRES_EXPLICIT_GATE0B_SCT_DENSITY_AND_CD8_MLC_QUALIFICATION'

    if s9_promoted:
        resolved_arch['pouch_cluster_diameter_cap_um'] = (
            80.0
            if (capsule_dz_upper95 > 18.7 or observations.get('s9_chronic_80um_lowload_qualified') is True)
            else 100.0
        )
        resolved_arch['required_qualification_coupon_max_cluster_diameter_um'] = resolved_arch['pouch_cluster_diameter_cap_um']
        if has_chronic_capsule_obs and capsule_dz_upper95 > 30.0:
            resolved_arch['chronic_capsule_gsis_derating_status'] = (
                'HOLD_CHRONIC_CAPSULE_GT_30UM_BREACHES_15MMHG_SURFACE_GSIS_REQUIRES_80UM_VIABILITY_ONLY_OR_VASCULARIZED_OMENTAL_SCAFFOLD'
            )
        if (
            capsule_dz_upper95 >= 14.0
            or _lower95('obs_d_o2_gel_aged_cm2_s', 1.5e-5) < 1.5e-5
            or _lower95('obs_capsule_d_o2_cm2_s', 1.5e-5) < 1.5e-5
            or periimplant_fluid_film_gt_50_triggered
        ):
            resolved_arch['pouch_areal_loading_cap_ieq_cm2'] = 1500.0
        if periimplant_fluid_film_gt_50_triggered:
            resolved_arch['periimplant_fluid_film_cluster_cap_status'] = (
                'ENFORCE_DXB_100UM_OR_80UM_CLUSTER_CAP_AND_1500_LOWLOAD_WITH_VERIFIED_DRAIN'
            )
    if periimplant_fluid_film_gt_50_triggered:
        resolved_arch['pouch_areal_loading_cap_ieq_cm2'] = 1500.0

    required_panel_complete = full_5endpoint_s9_panel_complete if require_full_5endpoint_s9_panel else mandatory_s9_safety_panel_complete
    if safety_kill_s9_triggered:
        safety_status = (
            'HARD_SAFETY_STOP_STICKY_PRIOR_S9_RETAINED'
            if (sticky_prior_s9_stop and not current_batch_s9_endpoint_failed)
            else 'HARD_SAFETY_STOP_REJECTED'
        )
        s9_panel_completeness_tier = 'REJECTED_BY_S9_HARD_SAFETY_STOP'
        resolved_arch['variant_safety_status'] = 'HARD_SAFETY_STOP_REJECTED'
        resolved_arch['lot_maturation_status'] = 'HARD_SAFETY_STOP_REJECTED'
    elif tp53_duplex_depth_underpowered_below_62958:
        safety_status = 'HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
        s9_panel_completeness_tier = 'HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
        resolved_arch['variant_safety_status'] = 'HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
        resolved_arch['lot_maturation_status'] = 'HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
    elif cast_duplex_depth_underpowered_below_62958:
        safety_status = 'HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
        s9_panel_completeness_tier = 'HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
        resolved_arch['variant_safety_status'] = 'HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
        resolved_arch['lot_maturation_status'] = 'HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM'
    elif mcb_clonogenic_outcome_missing_when_trials_supplied:
        safety_status = 'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING'
        s9_panel_completeness_tier = 'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING'
        resolved_arch['variant_safety_status'] = 'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING'
        resolved_arch['lot_maturation_status'] = 'HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING'
    elif mcb_clonogenic_trials_underpowered_below_299573228 or mcb_clonogenic_trials_missing_when_outcome_zero_supplied:
        safety_status = 'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM'
        s9_panel_completeness_tier = 'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM'
        resolved_arch['variant_safety_status'] = 'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM'
        resolved_arch['lot_maturation_status'] = 'HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM'
    elif residual_psc_trials_underpowered:
        safety_status = 'HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
        s9_panel_completeness_tier = 'HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
        resolved_arch['variant_safety_status'] = 'HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
        resolved_arch['lot_maturation_status'] = 'HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
    elif hazard_prolif_trials_underpowered:
        safety_status = 'HOLD_HAZARD_STATE_PROLIF_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
        s9_panel_completeness_tier = 'HOLD_HAZARD_STATE_PROLIF_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
        resolved_arch['variant_safety_status'] = 'HOLD_HAZARD_STATE_PROLIF_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
        resolved_arch['lot_maturation_status'] = 'HOLD_HAZARD_STATE_PROLIF_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM'
    elif ap1903_step1_trials_underpowered:
        safety_status = 'HOLD_AP1903_STEP1_EFFECTIVE_TRIALS_BELOW_29957323_MINIMUM'
        s9_panel_completeness_tier = 'HOLD_AP1903_STEP1_EFFECTIVE_TRIALS_BELOW_29957323_MINIMUM'
        resolved_arch['variant_safety_status'] = 'HOLD_AP1903_STEP1_EFFECTIVE_TRIALS_BELOW_29957323_MINIMUM'
        resolved_arch['lot_maturation_status'] = 'HOLD_AP1903_STEP1_EFFECTIVE_TRIALS_BELOW_29957323_MINIMUM'
    elif observations.get('obs_vascularized_explant_en_bloc_recovery_verified') is False:
        safety_status = 'HOLD_VASCULARIZED_EXPLANT_EN_BLOC_RECOVERY_UNVERIFIED'
        s9_panel_completeness_tier = 'HOLD_VASCULARIZED_EXPLANT_EN_BLOC_RECOVERY_UNVERIFIED'
        resolved_arch['variant_safety_status'] = 'HOLD_VASCULARIZED_EXPLANT_EN_BLOC_RECOVERY_UNVERIFIED'
        resolved_arch['lot_maturation_status'] = 'HOLD_S9_EXPLANT_EN_BLOC_RECOVERY_BELOW_99PCT'
    elif not required_panel_complete:
        safety_status = 'HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL'
        s9_panel_completeness_tier = (
            'HOLD_INCOMPLETE_5ENDPOINT_S9_PANEL'
            if (require_full_5endpoint_s9_panel and mandatory_s9_safety_panel_complete)
            else 'HOLD_INCOMPLETE_CORE_3ENDPOINT_S9_PANEL'
        )
        resolved_arch['variant_safety_status'] = 'HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL'
        if not s8_promoted and not str(resolved_arch.get('lot_maturation_status', '')).startswith(('HOLD_', 'NO_GO_')):
            resolved_arch['lot_maturation_status'] = 'HOLD_INCOMPLETE_SAFETY_PANEL'
    else:
        safety_status = 'CLEARED'
        s9_panel_completeness_tier = (
            'FULL_5ENDPOINT_S9_PANEL_CLEARED'
            if full_5endpoint_s9_panel_complete
            else 'CORE_3ENDPOINT_GENOMIC_SUICIDE_CLEARED_ER_GATE0E_PENDING'
        )
        resolved_arch['variant_safety_status'] = 'CLEARED'
        if prior_s9_reset_history_retained and resolved_arch.get('lot_maturation_status') == 'READY_FOR_ENCAPSULATION':
            resolved_arch['lot_maturation_status'] = 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'

    # Pivotal combined-stress Stage-0 coupon gate evaluation (STAGE0_COMBINED_STRESS_KILL_CHAIN_COUPON_MUSE_R4; AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, MW-R5-S2-F04 & AUDIT-REF; DEC-043..050):
    # Prohibits positive defaults: PROMOTE of ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500 ALWAYS requires all 4 COMP_R4 secondary metrics, and single-arm RETAIN of Embodiment A* accepts either COMP_R4 or Embodiment A* 4-endpoint secondary panels.
    combined_stress_spec = dict(tree_doc.get('pivotal_combined_stress_stage0_coupon_gate_muse_r4', {}))
    has_a_star_cs = 'obs_combined_stress_embodiment_a_star_lysis_pct' in validated_stage0
    has_comp_cs_lysis = 'obs_combined_stress_comp_r4_lysis_pct' in validated_stage0
    has_comp_cs_si = 'obs_combined_stress_comp_r4_gsis_si' in validated_stage0
    has_all_4_comp_cs_secondary = bool(
        'obs_combined_stress_comp_r4_viable_beta_mass_pct' in validated_stage0
        and 'obs_combined_stress_comp_r4_shutoff_t90_min' in validated_stage0
        and 'obs_combined_stress_comp_r4_max_h2o2_uM' in validated_stage0
        and 'obs_combined_stress_comp_r4_max_ph' in validated_stage0
    )
    has_all_4_a_star_cs_secondary = bool(
        'obs_combined_stress_embodiment_a_star_viable_beta_mass_pct' in validated_stage0
        and 'obs_combined_stress_embodiment_a_star_shutoff_t90_min' in validated_stage0
        and 'obs_combined_stress_embodiment_a_star_max_h2o2_uM' in validated_stage0
        and 'obs_combined_stress_embodiment_a_star_max_ph' in validated_stage0
    )
    has_all_4_cs_secondary = bool(has_all_4_comp_cs_secondary or has_all_4_a_star_cs_secondary)
    cs_secondary_failed = bool(
        ('obs_combined_stress_comp_r4_viable_beta_mass_pct' in validated_stage0 and _lower95('obs_combined_stress_comp_r4_viable_beta_mass_pct') < 85.0)
        or ('obs_combined_stress_comp_r4_shutoff_t90_min' in validated_stage0 and _upper95('obs_combined_stress_comp_r4_shutoff_t90_min') > 3.0)
        or ('obs_combined_stress_comp_r4_max_h2o2_uM' in validated_stage0 and _upper95('obs_combined_stress_comp_r4_max_h2o2_uM') > 5.0)
        or ('obs_combined_stress_comp_r4_max_ph' in validated_stage0 and _upper95('obs_combined_stress_comp_r4_max_ph') > 7.50)
        or ('obs_combined_stress_embodiment_a_star_viable_beta_mass_pct' in validated_stage0 and _lower95('obs_combined_stress_embodiment_a_star_viable_beta_mass_pct') < 85.0)
        or ('obs_combined_stress_embodiment_a_star_shutoff_t90_min' in validated_stage0 and _upper95('obs_combined_stress_embodiment_a_star_shutoff_t90_min') > 3.0)
        or ('obs_combined_stress_embodiment_a_star_max_h2o2_uM' in validated_stage0 and _upper95('obs_combined_stress_embodiment_a_star_max_h2o2_uM') > 5.0)
        or ('obs_combined_stress_embodiment_a_star_max_ph' in validated_stage0 and _upper95('obs_combined_stress_embodiment_a_star_max_ph') > 7.50)
    )
    has_combined_stress_obs = bool(
        has_a_star_cs
        or has_comp_cs_lysis
        or has_comp_cs_si
        or 'obs_combined_stress_comp_r4_viable_beta_mass_pct' in validated_stage0
        or 'obs_combined_stress_comp_r4_shutoff_t90_min' in validated_stage0
        or 'obs_combined_stress_comp_r4_max_h2o2_uM' in validated_stage0
        or 'obs_combined_stress_comp_r4_max_ph' in validated_stage0
        or 'obs_combined_stress_embodiment_a_star_viable_beta_mass_pct' in validated_stage0
        or 'obs_combined_stress_embodiment_a_star_shutoff_t90_min' in validated_stage0
        or 'obs_combined_stress_embodiment_a_star_max_h2o2_uM' in validated_stage0
        or 'obs_combined_stress_embodiment_a_star_max_ph' in validated_stage0
        or observations.get('combined_stress_coupon_cleared_lysis_le_8pct_and_si_ge_2p5') is not None
    )
    combined_stress_spec['full_7metric_panel_complete'] = bool(
        (has_a_star_cs and has_comp_cs_lysis and has_comp_cs_si and has_all_4_comp_cs_secondary)
        or (has_a_star_cs and has_all_4_a_star_cs_secondary)
    )
    combined_stress_spec['secondary_panel_provenance'] = (
        'EXPLICIT_4SECONDARY_ENDPOINTS_OBSERVED'
        if has_all_4_cs_secondary
        else 'CORE_3COMPARATOR_ONLY_SECONDARY_VIABILITY_SHUTOFF_H2O2_PH_PENDING'
    )
    combined_stress_promoted = False
    if has_combined_stress_obs:
        if observations.get('combined_stress_coupon_cleared_lysis_le_8pct_and_si_ge_2p5') is False or cs_secondary_failed:
            combined_stress_coupon_status = 'HOLD_COMBINED_STRESS_COUPON_FAILED'
            if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                resolved_arch['lot_maturation_status'] = 'HOLD_COMBINED_STRESS_COUPON_FAILED'
        elif observations.get('require_full_7metric_combined_stress_panel') is True and not has_all_4_cs_secondary:
            combined_stress_coupon_status = 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
            if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                resolved_arch['lot_maturation_status'] = 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
        elif has_a_star_cs and has_comp_cs_lysis and has_comp_cs_si:
            a_star_cs_lysis_u95 = _upper95('obs_combined_stress_embodiment_a_star_lysis_pct')
            comp_cs_lysis_u95 = _upper95('obs_combined_stress_comp_r4_lysis_pct')
            comp_cs_si_l95 = _lower95('obs_combined_stress_comp_r4_gsis_si')
            if resolved_arch.get('lot_maturation_status') not in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                combined_stress_coupon_status = 'HOLD_COMPONENT_GATE_FAILURE_BLOCKS_COMPOSITE_PROMOTION'
            elif a_star_cs_lysis_u95 > 8.0 and comp_cs_lysis_u95 <= 8.0 and comp_cs_si_l95 >= 2.5:
                # Composite promotion ALWAYS requires all 4 COMP_R4 secondary endpoints (no waiver via require_full_7metric_combined_stress_panel=False; AUDIT-REF & AUDIT-REF; DEC-050/051):
                if not has_all_4_comp_cs_secondary:
                    hold_cs_tok = (
                        'HOLD_COMP_R4_PROMOTION_REQUIRES_COMPLETE_4METRIC_SECONDARY_PANEL'
                        if observations.get('require_full_7metric_combined_stress_panel') is False
                        else 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
                    )
                    combined_stress_coupon_status = hold_cs_tok
                    if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                        resolved_arch['lot_maturation_status'] = hold_cs_tok
                else:
                    combined_stress_coupon_status = 'PROMOTE_ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500'
                    combined_stress_promoted = True
                    resolved_arch['combined_stress_coupon_promoted_variant'] = 'ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500'
            elif a_star_cs_lysis_u95 <= 8.0 and comp_cs_lysis_u95 <= 8.0 and comp_cs_si_l95 >= 2.5:
                if observations.get('require_full_7metric_combined_stress_panel', True) is not False and not has_all_4_cs_secondary:
                    combined_stress_coupon_status = 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
                    if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                        resolved_arch['lot_maturation_status'] = 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
                else:
                    combined_stress_coupon_status = 'RETAIN_EMBODIMENT_A_STAR_4143BP'
            else:
                combined_stress_coupon_status = 'HOLD_COMBINED_STRESS_COUPON_FAILED'
                if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                    resolved_arch['lot_maturation_status'] = 'HOLD_COMBINED_STRESS_COUPON_FAILED'
        elif has_a_star_cs and _upper95('obs_combined_stress_embodiment_a_star_lysis_pct') <= 8.0 and not has_comp_cs_lysis and not has_comp_cs_si:
            if observations.get('require_full_7metric_combined_stress_panel', True) is not False and not has_all_4_cs_secondary:
                combined_stress_coupon_status = 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
                if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                    resolved_arch['lot_maturation_status'] = 'HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL'
            else:
                combined_stress_coupon_status = 'RETAIN_EMBODIMENT_A_STAR_4143BP'
        else:
            combined_stress_coupon_status = 'HOLD_COMBINED_STRESS_INCOMPLETE_COMPARATOR_EVIDENCE'
            if resolved_arch.get('lot_maturation_status') in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'):
                resolved_arch['lot_maturation_status'] = 'HOLD_COMBINED_STRESS_INCOMPLETE_COMPARATOR_EVIDENCE'
    else:
        combined_stress_coupon_status = 'PENDING_COMBINED_STRESS_COUPON_EXECUTION'
    combined_stress_spec['evaluation_status'] = combined_stress_coupon_status

    has_explicit_maturity_obs = bool(
        'obs_gate0m_shutoff_t90_min' in validated_stage0
        and 'obs_basal_insulin_frac_peak' in validated_stage0
        and 'obs_glucose_threshold_mM' in validated_stage0
    )
    has_explicit_immune_and_ocr_core_obs = bool(
        'obs_hlae_hour0_molecules' in validated_stage0
        and 'obs_cd47_hour0_molecules' in validated_stage0
        and 'obs_nkg2a_lysis_pct' in validated_stage0
        and ('obs_serpinb9_hour0_uM' in validated_stage0 or 'obs_gate0c_hour0_serpinb9_uM' in validated_stage0)
        and ('obs_regimeII_lysis_2h_pct' in validated_stage0 or 'obs_gate0c_2h_lysis_pct' in validated_stage0)
        and 'obs_lot_ocr_ratio' in validated_stage0
    )
    has_explicit_indirect_apc_gate0d_obs = bool(
        (
            has_complete_arm2_6endpoint_panel
            and not arm2_depletion_verification_incomplete
            and not arm2_rebound_injury_or_function_failed
        )
        if adjunct_conditioning_history_arm == 'ARM2_CONDITIONED_DEPLETED'
        else (
            adjunct_conditioning_history_arm == 'ARM3_CONTINUOUSLY_ACTIVE'
            and 'obs_indirect_apc_dq2_dq8_cd4_stimulation_index' in validated_stage0
            and 'obs_indirect_apc_cd4_ifng_si' in validated_stage0
            and 'obs_gate0d_activated_cd4_fasl_apoptosis_pct' in validated_stage0
            and 'obs_gate0d_transwell_beta_annexin_v_pct' in validated_stage0
            and not indirect_apc_cd4_failed
            and not arm3_ongoing_bioactivity_explicit_false
        )
    )
    has_explicit_exogenous_b2m_obs = bool(
        'obs_exogenous_b2m_w632_surface_fold_isotype' in validated_stage0
        and 'obs_exogenous_b2m_t1d_ctl_lysis_pct' in validated_stage0
        and not exogenous_b2m_class1_rescue_failed
    )
    has_explicit_glucagon_and_perfusion_obs = bool(
        ('obs_gate0m_glucagon_2mM_fold_increase' in validated_stage0 or 'obs_alpha_glucagon_low_glc_fold' in validated_stage0)
        and 'obs_outer_band_kperf_s_inv' in validated_stage0
        and _lower95('obs_outer_band_kperf_s_inv') >= 0.10
    )
    has_explicit_dec046_astra_r5_obs = bool(
        'obs_cd47_cc2c6_over_b6h12_ratio' in validated_stage0
        and 'obs_complement_c3b_opsonized_phagocytosis_pct' in validated_stage0
        and 'obs_iec_vascular_support_viability_pct' in validated_stage0
        and 'obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct' in validated_stage0
        and 'obs_12pulse_unprimed_hypoxic_serial_lysis_pct' in validated_stage0
        and 'obs_recipient_insulin_clearance_tau_min' in validated_stage0
        and 'obs_hlae_nkg2c_over_nkg2a_kd_ratio' in validated_stage0
        and 'obs_il12_il18_sirpa_low_nk_lysis_pct' in validated_stage0
        and 'obs_immaturity_corner_min_glucose_mM' in validated_stage0
        and 'obs_alpha_counterregulatory_gain_ratio' in validated_stage0
        and 'obs_ucn3_maturity_ratio' in validated_stage0
    )
    has_explicit_dec047_tumour_and_retention_obs = bool(
        'obs_residual_psc_effective_fraction' in validated_stage0
        and residual_psc_trials_supplied
        and residual_psc_trials >= 299573228.0
        and 'obs_20q11_bcl2l1_cnv_ratio' in validated_stage0
        and 'obs_20q11_bcl2l1_single_cell_mosaic_fraction' in validated_stage0
        and 'obs_hazard_state_prolif_survivor_count_k' in validated_stage0
        and hazard_prolif_trials_supplied
        and hazard_prolif_trials >= 299573228.0
        and 'obs_cd47_p10_single_cell_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_cd47_p10_single_cell_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and validated_stage0['obs_cd47_p10_single_cell_molecules'] >= 180000.0
        and _lower95('obs_cd47_p10_single_cell_molecules') >= 180000.0
        and 'obs_hypoxic_72h_cd47_p10_single_cell_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_hypoxic_72h_cd47_p10_single_cell_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and validated_stage0['obs_hypoxic_72h_cd47_p10_single_cell_molecules'] >= 180000.0
        and _lower95('obs_hypoxic_72h_cd47_p10_single_cell_molecules') >= 180000.0
        and 'obs_hypoxic_72h_hlae_p10_single_cell_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_hypoxic_72h_hlae_p10_single_cell_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and validated_stage0['obs_hypoxic_72h_hlae_p10_single_cell_molecules'] >= 15000.0
        and _lower95('obs_hypoxic_72h_hlae_p10_single_cell_molecules') >= 15000.0
        and 'obs_hypoxic_72h_cd47_retention_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_hypoxic_72h_cd47_retention_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and validated_stage0['obs_hypoxic_72h_cd47_retention_molecules'] >= 180000.0
        and _lower95('obs_hypoxic_72h_cd47_retention_molecules') >= 180000.0
        and 'obs_hypoxic_72h_hlae_retention_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_hypoxic_72h_hlae_retention_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and validated_stage0['obs_hypoxic_72h_hlae_retention_molecules'] >= 15000.0
        and _lower95('obs_hypoxic_72h_hlae_retention_molecules') >= 15000.0
        and 'obs_periimplant_fluid_film_um' in validated_stage0
        and 'obs_explant_en_bloc_recovery_fraction' in validated_stage0
        and 'obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h' in validated_stage0
        and 'obs_gate0s_ap1903_1h_insulin_burst_iu' in validated_stage0
        and _upper95('obs_gate0s_ap1903_1h_insulin_burst_iu') <= 2.8
        and not residual_psc_failed
        and not chr20q11_bcl2l1_cnv_failed
        and not hazard_state_prolif_survivor_detected
        and not periimplant_fluid_film_failed
        and not explant_retrieval_failed
        and not ap1903_insulin_dump_unmitigated
        and not mandatory_hypoxic_or_p10_retention_cliff_failed
        and not unmitigated_cd47_or_hypoxic_retention_cliff_failed
        and not zero_event_proportion_missing_denominator_or_explicit_se
    )
    has_explicit_dec048_muse_r5_obs = bool(
        'obs_ap1903_pulsed_pk_log10_kill' in validated_stage0
        and 'obs_ap1903_monotherapy_pulsed_pk_log10_kill' in validated_stage0
        and ap1903_mono_pk_lower95 >= 4.0
        and 'obs_ap1903_step1_sublethal_survivor_count_k' in validated_stage0
        and ap1903_step1_trials_supplied
        and ap1903_step1_trials >= 29957322.8
        and 'obs_hypoxic_15mmhg_f_isr_translation_ratio' in validated_stage0
        and _lower95('obs_hypoxic_15mmhg_f_isr_translation_ratio') >= 0.50
        and 'obs_iec_cd47_p10_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_iec_cd47_p10_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and _lower95('obs_iec_cd47_p10_molecules') >= 180000.0
        and 'obs_iec_abo_ab_antigen_mfi_ratio' in validated_stage0
        and _upper95('obs_iec_abo_ab_antigen_mfi_ratio') <= 1.05
        and 'obs_post_stress_nkg2d_dnam1_attributed_lysis_pct' in validated_stage0
        and _upper95('obs_post_stress_nkg2d_dnam1_attributed_lysis_pct') <= 5.0
        and 'obs_senescent_p16_sasp_positive_cell_pct' in validated_stage0
        and senescent_p16_upper95 <= 1.0
        and not ap1903_pulsed_pk_kill_failed
        and not ap1903_step1_sublethal_survivor_detected
        and not iec_vascular_support_failed
        and not post_stress_nkg2d_dnam1_lysis_failed
        and not senescent_p16_sasp_failed
        and not zero_event_proportion_missing_denominator_or_explicit_se
    )
    # Single authoritative policy table for all 13 canonical DEC-053/055/057 endpoints required directly (no supplementary alias substitution; R7S2-F01; DEC-057):
    has_explicit_dec053_muse_r6_obs = bool(
        'obs_gate0v_cvb_day14_log10_viral_reduction' in validated_stage0
        and validated_stage0['obs_gate0v_cvb_day14_log10_viral_reduction'] >= 3.0
        and _lower95('obs_gate0v_cvb_day14_log10_viral_reduction') >= 3.0
        and 'obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct' in validated_stage0
        and validated_stage0['obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct'] <= 0.0
        and vp1_ie1_upper95 <= 0.05
        and 'obs_post_stress_20q11_bcl2l1_outgrowth_fold' in validated_stage0
        and validated_stage0['obs_post_stress_20q11_bcl2l1_outgrowth_fold'] <= 1.0
        and _upper95('obs_post_stress_20q11_bcl2l1_outgrowth_fold') <= 1.20
        and 'obs_assembled_cassette_chandler_loop_sc5b9_fold' in validated_stage0
        and validated_stage0['obs_assembled_cassette_chandler_loop_sc5b9_fold'] <= 2.0
        and _upper95('obs_assembled_cassette_chandler_loop_sc5b9_fold') <= 2.0
        and 'obs_assembled_cassette_chandler_loop_platelet_retention_pct' in validated_stage0
        and validated_stage0['obs_assembled_cassette_chandler_loop_platelet_retention_pct'] >= 70.0
        and _lower95('obs_assembled_cassette_chandler_loop_platelet_retention_pct') >= 70.0
        and 'obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction' in validated_stage0
        and validated_stage0['obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction'] >= 0.90
        and _lower95('obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction') >= 0.90
        and 'obs_apoptotic_effero_clearance_retention_pct_of_wt' in validated_stage0
        and validated_stage0['obs_apoptotic_effero_clearance_retention_pct_of_wt'] >= 70.0
        and _lower95('obs_apoptotic_effero_clearance_retention_pct_of_wt') >= 70.0
        and 'obs_nutlin3a_p21_mdm2_fold_induction' in validated_stage0
        and validated_stage0['obs_nutlin3a_p21_mdm2_fold_induction'] >= 5.0
        and _lower95('obs_nutlin3a_p21_mdm2_fold_induction') >= 5.0
        and 'obs_vascular_pericyte_ng2_pdgfrb_coverage_pct' in validated_stage0
        and validated_stage0['obs_vascular_pericyte_ng2_pdgfrb_coverage_pct'] >= 70.0
        and _lower95('obs_vascular_pericyte_ng2_pdgfrb_coverage_pct') >= 70.0
        and 'obs_assembled_cluster_slc18a1_enterochromaffin_pct' in validated_stage0
        and validated_stage0['obs_assembled_cluster_slc18a1_enterochromaffin_pct'] <= 1.0
        and slc18a1_ec_upper95 <= 1.0
        and 'obs_assembled_cluster_ins_gcg_polyhormonal_pct' in validated_stage0
        and validated_stage0['obs_assembled_cluster_ins_gcg_polyhormonal_pct'] <= 2.0
        and _upper95('obs_assembled_cluster_ins_gcg_polyhormonal_pct') <= 2.0
        and 'obs_assembled_cluster_ki67_progenitor_pct' in validated_stage0
        and validated_stage0['obs_assembled_cluster_ki67_progenitor_pct'] <= 0.5
        and ki67_prog_upper95 <= 0.5
        and 'obs_post_stress_mafa_nkx61_beta_identity_retention_pct' in validated_stage0
        and validated_stage0['obs_post_stress_mafa_nkx61_beta_identity_retention_pct'] >= 85.0
        and _lower95('obs_post_stress_mafa_nkx61_beta_identity_retention_pct') >= 85.0
        and not post_stress_20q11_outgrowth_failed
        and not nutlin3a_p53_functional_failed
        and not gate0v_viral_sanctuary_failed
        and not assembled_cassette_chandler_ibmir_failed
        and not effero_clearance_failed
        and not assembled_cluster_lineage_or_identity_failed
        and not zero_event_proportion_missing_denominator_or_explicit_se
    )
    # Mandatory original exact-product primary biological & pouch endpoints required for READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE (R7S3-F01; DEC-058; Design Z Zero-CaO2 H2O2 exemption IIP-0024; DEC-103):
    # Prevents composite COMP-R4 combined-stress secondary coupons or omitted/UNAVAILABLE primary pouch oxygen/peroxide/pH/HLA-II/yield/CD47/SERPINB9/iCasp9 observations from awarding COMPLETE Stage-0 readiness.
    _orig_primary_keys_for_product = tuple(
        mk for mk in STAGE0_ORIGINAL_PRIMARY_MANDATORY_MEAN_KEYS
        if not (_is_design_z_product and mk in ('obs_pouch_min_po2_mmhg', 'obs_pouch_max_h2o2_uM', 'obs_pouch_max_ph'))
    )
    missing_original_primary_endpoints = [
        mk for mk in _orig_primary_keys_for_product if mk not in validated_stage0
    ]
    has_explicit_original_primary_obs = bool(
        not missing_original_primary_endpoints
        and all(
            se_provenance_by_mean.get(mk) == 'EXPLICIT_SE_SUPPLIED'
            for mk in _orig_primary_keys_for_product
        )
        and validated_stage0['obs_stage7_yield_ratio'] >= 0.70
        and _lower95('obs_stage7_yield_ratio') >= 0.70
        and validated_stage0['obs_cd47_fold'] >= 3.0
        and _lower95('obs_cd47_fold') >= 3.0
        and validated_stage0['obs_serpinb9_log10_kass'] >= 6.0
        and _lower95('obs_serpinb9_log10_kass') >= 6.0
        and validated_stage0['obs_icasp9_log10_kill'] >= 4.0
        and icasp9_lower95 >= 4.0
        and (_is_design_z_product or (validated_stage0['obs_pouch_min_po2_mmhg'] >= 15.0 and _lower95('obs_pouch_min_po2_mmhg') >= 15.0))
        and (_is_design_z_product or (validated_stage0['obs_pouch_max_h2o2_uM'] <= 5.0 and _upper95('obs_pouch_max_h2o2_uM') <= 5.0))
        and (_is_design_z_product or (validated_stage0['obs_pouch_max_ph'] <= 7.45 and _upper95('obs_pouch_max_ph') <= 7.45))
        and validated_stage0['obs_ifng_hla_dr_dq_negative_pct'] >= 99.0
        and _lower95('obs_ifng_hla_dr_dq_negative_pct') >= 99.0
        and not primary_biological_or_pouch_gate_failed
    )
    has_explicit_release_critical_molecule_and_rate_se = bool(
        'obs_cd47_hour0_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_cd47_hour0_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and 'obs_hlae_hour0_molecules' in validated_stage0
        and se_provenance_by_mean.get('obs_hlae_hour0_molecules') == 'EXPLICIT_SE_SUPPLIED'
        and all(
            se_provenance_by_mean.get(mk) == 'EXPLICIT_SE_SUPPLIED'
            for mk in STAGE0_RELEASE_CRITICAL_MOLECULE_AND_RATE_MEAN_KEYS
            if mk in validated_stage0
        )
    )
    has_explicit_extended_stage0_obs = bool(
        has_explicit_original_primary_obs
        and has_explicit_release_critical_molecule_and_rate_se
        and 'obs_nkg2c_adaptive_nk_lysis_pct' in validated_stage0
        and 'obs_kir_nkg2a_neg_nk_lysis_pct' in validated_stage0
        and 'obs_phagocytosis_index_pct' in validated_stage0
        and tp53_duplex_depth_explicitly_supplied
        and vaf_duplex_depth >= 62958.0
        and cast_duplex_depth_explicitly_supplied
        and cast_duplex_depth >= 62958.0
        and mcb_clonogenic_trials_explicitly_supplied
        and mcb_clonogenic_trials >= 299573228.0
        and mcb_survivors_k_explicitly_supplied
        and mcb_survivors_k == 0.0
        and not mcb_clonogenic_survivor_detected
        and has_chronic_capsule_obs
        and has_explicit_indirect_apc_gate0d_obs
        and has_explicit_exogenous_b2m_obs
        and has_explicit_glucagon_and_perfusion_obs
        and has_explicit_dec046_astra_r5_obs
        and has_explicit_dec047_tumour_and_retention_obs
        and has_explicit_dec048_muse_r5_obs
        and has_explicit_dec053_muse_r6_obs
        and combined_stress_spec.get('full_7metric_panel_complete') is True
        and combined_stress_coupon_status == 'RETAIN_EMBODIMENT_A_STAR_4143BP'
        and not is_comp_r4_product
        and (not _req_iip0020_full or (iip0020_hold_status is None and (_dz_arm_val is not None or iip0020_p05_species_status == 'CLEARED_P05_HUMAN_P04040_SPECIES_ENDOTOXIN_AND_SHED_MASS_GATE') and iip0020_cuff_status == 'CLEARED_CUFF_HUMORAL_HYDROLYSIS_AND_HEPATOCYTE_SAFETY_GATE'))
        and (not _req_iip0021_full or (iip0021_hold_status is None and iip0021_rest_status is not None and iip0021_cond_status == 'CLEARED_CONDITIONING_4ARM_PELLET_LOQ_NQO1_AND_LINEAGE_GSIS_GATE'))
        and (not _req_iip0022_full or iip0022_hold_status is None)
        and (not _req_iip0024_full or iip0024_hold_status is None)
        and not non_integer_scored_cell_or_trial_denominator
        and not mosaic_explicit_se_missing_scored_cell_denominator
    )
    any_fallback_switch_triggered = bool(
        s1_promoted
        or s2_promoted
        or s3_promoted
        or s4_promoted
        or s5_promoted
        or s6_promoted
        or s7_allo_triggered
        or s8_promoted
        or s8_tyk2_epistasis_triggered
        or s9_promoted
        or (s10_nkg2c_triggered and not s10_v3_release_qualified)
        or (_s10_v3_active and not s10_v3_release_qualified)
        or combined_stress_promoted
        or track2_kir_or_hlag1_fallback_promoted
    )
    any_resolved_component_on_hold = any(
        str(resolved_arch.get(k, '')).startswith('HOLD_')
        for k in ('locus1_construct', 'locus2_construct', 'locus3_construct', 'locus4_construct', 'track2_kir_subbank_construct')
    )
    _prod_non_executable = bool(product_def is not None and product_def.get('execution_supported') is False)
    _pvid_for_tier = str(product_def.get('product_variant_id', '')) if product_def is not None else ''
    if safety_kill_s9_triggered:
        lot_release_readiness_tier = 'HARD_SAFETY_STOP_REJECTED'
    elif iip0031_tier6b_eval is not None and iip0031_tier6b_eval.get('hard_stop') is True:
        lot_release_readiness_tier = 'HARD_SAFETY_STOP_REJECTED'
        resolved_arch['lot_maturation_status'] = 'FAIL_TIER6B_OUTGROWTH_DETECTED'
    elif str(resolved_arch.get('lot_maturation_status', '')).startswith('NO_GO_'):
        lot_release_readiness_tier = 'NO_GO_TERMINAL_BIOTRANSPORT_OR_IMMUNE_CEILING_EXCEEDED'
    elif resolved_arch.get('lot_maturation_status') not in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY') or any_resolved_component_on_hold:
        lot_release_readiness_tier = 'HOLD_FAILED_OR_INCOMPLETE_GATE'
    elif is_comp_r4_product and combined_stress_coupon_status == 'RETAIN_EMBODIMENT_A_STAR_4143BP':
        lot_release_readiness_tier = 'CANDIDATE_NOT_INDICATED_RETAIN_EMBODIMENT_A_STAR_4143BP'
    elif _prod_non_executable and (_pvid_for_tier == 'ISLET1_B2_5LOCUS_BASE' or _pvid_for_tier in _staged_l2_v3_ids):
        lot_release_readiness_tier = 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
        resolved_arch['lot_maturation_status'] = 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
    elif _prod_non_executable and _pvid_for_tier == 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2':
        if any_fallback_switch_triggered:
            lot_release_readiness_tier = 'CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION'
        elif full_5endpoint_s9_panel_complete and has_explicit_maturity_obs and has_explicit_immune_and_ocr_core_obs and has_explicit_extended_stage0_obs:
            lot_release_readiness_tier = 'WILDCARD_COUPON_QUALIFIED_NON_EXECUTABLE'
            resolved_arch['lot_maturation_status'] = 'WILDCARD_COUPON_QUALIFIED_NON_EXECUTABLE'
        else:
            lot_release_readiness_tier = 'HOLD_FAILED_OR_INCOMPLETE_GATE'
            resolved_arch['lot_maturation_status'] = 'HOLD_INCOMPLETE_WILDCARD_STAGE0_EVIDENCE'
    elif _prod_non_executable:
        lot_release_readiness_tier = 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD'
        resolved_arch['lot_maturation_status'] = 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD'
    elif any_fallback_switch_triggered:
        lot_release_readiness_tier = 'CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION'
    elif full_5endpoint_s9_panel_complete and has_explicit_maturity_obs and has_explicit_immune_and_ocr_core_obs and has_explicit_extended_stage0_obs:
        lot_release_readiness_tier = (
            'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'
            if prior_s9_reset_history_retained
            else 'READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE'
        )
    elif full_5endpoint_s9_panel_complete and has_explicit_maturity_obs and has_explicit_immune_and_ocr_core_obs:
        lot_release_readiness_tier = (
            'PROVISIONAL_CORE_14_PANEL_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY'
            if prior_s9_reset_history_retained
            else 'PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED'
        )
    elif full_5endpoint_s9_panel_complete and has_explicit_maturity_obs:
        lot_release_readiness_tier = (
            'PROVISIONAL_S9_AND_MATURITY_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY'
            if prior_s9_reset_history_retained
            else 'PROVISIONAL_S9_AND_MATURITY_CLEARED_IMMUNE_AND_OCR_GATES_UNMEASURED'
        )
    else:
        lot_release_readiness_tier = (
            'PROVISIONAL_CORE_S9_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY'
            if prior_s9_reset_history_retained
            else 'PROVISIONAL_CORE_S9_CLEARED_UNMEASURED_GATES_PENDING'
        )
    resolved_arch['lot_release_readiness_tier'] = lot_release_readiness_tier

    decisions = {
        'S1_Locus3_Kinetic_and_Topology_Gate_0C': (
            'NO_GO_ACUTE_REGIME_II_LYSIS_GT_25PCT_TERMINAL'
            if lysis_2h_upper95 > 25.0
            else ('PROMOTE_ISLET1_A4_L3_BFB_BASALFLOOR' if s1_promoted else 'RETAIN_CHAMPION_LOCUS3_OR_L3_SCARFREE_3799')
        ),
        'S2_Locus1_Arm1B_vs_Arm1X_Gate_0A': (
            'HOLD_S2_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_REQUIRES_MULTI_RECEPTOR_RESCUE'
            if catastrophic_nkg2a_lysis_exceeds_arm1x
            else (
                'PROMOTE_ISLET1_A4_X_MOL_RNLS_WT'
                if s2_promoted
                else 'RETAIN_EMBODIMENT_A_STAR_ARM1B'
            )
        ),
        'S3_Locus1_Constitutive_Split_Fallback_DX_S': (
            'PROMOTE_ISLET1_A4_DXS_ARM1C_L4CD47'
            if s3_promoted
            else (
                'HOLD_DXS_RARE_EVENT_REQUIRES_EXPLICIT_SE_OR_CELL_COUNT'
                if dxs_hold_missing_explicit_se_or_counts
                else 'INACTIVE'
            )
        ),
        'S4_Metabolic_Demand_and_Autoreactive_Shield_Gate_0D': 'LOCK_RNLS_WT_AND_TXNIP_KO' if s4_promoted else 'RETAIN_BASELINE_METABOLIC_CHASSIS',
        'S5_Microgel_Boundary_Cuff_Gate_0G': 'DROP_SA_FASL_RETAIN_SA_THBD_PLUS_SA_FHL1' if s5_promoted else 'RETAIN_PERIMETER_SA_FASL_CUFF',
        'S6_B2M_Orphan_Heavy_Chain_ERAD_Gate_0E': 'ACTIVATE_TYK2_P1104A_OR_RNLS_KO_OR_NLRC5_CRISPROFF' if s6_promoted else 'PASS_NO_ERAD_UPR_PENALTY',
        'S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B': (
            'PROMOTE_ISLET1_B_SCT_C1_OR_C2_CIITA'
            if s7_promoted
            else (
                'HOLD_B_PRIME_SCT_DIRECT_CD8_ALLOREACTIVITY_FAILED_RETAIN_PURE_MHC_NULL_TRACK1'
                if b_prime_sct_direct_cd8_failed
                else (
                    'HOLD_B_PRIME_ALLO_FAILED_PENDING_B_SCT_DENSITY_QUALIFICATION'
                    if s7_hold_unqualified_sct
                    else (
                        'HOLD_B_PRIME_ALLO_FAILED_REQUIRES_4DIGIT_MATCH_OR_SCT_CD8_MLC_CLEARANCE'
                        if s7_hold_unverified_sct_allo
                        else (
                            b2_tyk2_gate_status
                            if b2_tyk2_gate_status is not None
                            else (
                                (
                                    iip0026_hlac_status
                                    if iip0026_hlac_status is not None
                                    else 'HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE'
                                )
                                if _is_endogenous_hlac_b_26
                                else (
                                    'UNMEASURED_EXTENDED_S7_GATE0B_HLAC_PANEL'
                                    if not _has_explicit_s7_hlac_obs
                                    else 'RETAIN_4DIGIT_MATCHED_HLAC'
                                )
                            )
                        )
                    )
                )
            )
        ),
        'S8_SC_Islet_Paracrine_Maturity_Gate_0M': 'MANDATE_EXTENDED_S7W6_MATURATION_AND_80_100UM_CLUSTERS' if s8_promoted else 'PASS_MATURITY_GATE_0M',
        'S8b_TYK2_Inducible_Promoter_Epistasis_Gate_0E': (
            'COUPLE_TYK2_P1104A_STRICTLY_WITH_ARM1X_AND_L3_BFB_CONSTITUTIVE_BASAL_FLOOR'
            if s8_tyk2_epistasis_triggered
            else (
                'HOLD_INCONCLUSIVE_UNDERPOWERED_TYK2_EPISTASIS_SE'
                if s8_tyk2_epistasis_underpowered_hold
                else 'PASS_TYK2_INDUCIBLE_PROMOTER_COMPATIBILITY'
            )
        ),
        'S9_Chronic_Fibrotic_Capsule_and_Mineralization_Gate_0F': 'ENFORCE_80UM_CLUSTER_CAP_AND_CITRATE_CHELATION' if s9_promoted else 'PASS_CHRONIC_CAPSULE_GATE_0F',
        'S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A': (
            (
                str(s10_v3_eval['status'])
                if not s10_v3_eval['gate_pass']
                else (
                    'CLEARED_S10_FUNCTIONAL_SELECTIVITY_SCREENING_15K_V3'
                    if not s10_v3_release_qualified
                    else (
                        'CLEARED_S10_FUNCTIONAL_SELECTIVITY_SUB4_SCT_V3'
                        if (s10_kd_ratio_lower95_sub4 or float(s10_v3_eval['recorded_sct_kc_over_ka_mean']) < 4.0)
                        else 'CLEARED_S10_STATIC_RATIO_AND_PAIRED_FUNCTIONAL_GATE_V3'
                    )
                )
            )
            if (_s10_v3_active and s10_v3_eval is not None)
            else (
                'CLEARED_SCT_SELECTIVITY_COUPON_KC_KA_GE_4'
                if (s10_nkg2c_triggered and s10_coupon_cleared)
                else (
                    'STRATIFY_CMV_NKG2C_RECIPIENTS_AND_EVALUATE_SCT_SELECTIVITY_MUTANT_COUPON'
                    if s10_nkg2c_triggered
                    else (
                        'UNMEASURED_EXTENDED_S10_PANEL'
                        if not _has_explicit_s10_obs
                        else 'PASS_CMV_NKG2C_SELECTIVITY_GATE_0A'
                    )
                )
            )
        ),
        'S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate': str(iip0031_s11_eval['status']),
    }
    canonical_json_switch_id_map = {
        'JSON_S1_CD47_Split_Fallback_DX_S': decisions['S3_Locus1_Constitutive_Split_Fallback_DX_S'],
        'JSON_S2_HLAE_Sense_Promoter_Gate_0A': decisions['S2_Locus1_Arm1B_vs_Arm1X_Gate_0A'],
        'JSON_S3_Locus3_Hour0_Kinetic_Gate_0C': decisions['S1_Locus3_Kinetic_and_Topology_Gate_0C'],
        'JSON_S4_Locus3_P2A_Tail_Biophysics_Gate_0C': decisions['S1_Locus3_Kinetic_and_Topology_Gate_0C'],
        'JSON_S5_Plan_B_4Digit_Allo_vs_B_SCT_Gate_0B': decisions['S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B'],
        'JSON_S6_Chronic_Capsule_and_Mineralization_Gate_0F': decisions['S9_Chronic_Fibrotic_Capsule_and_Mineralization_Gate_0F'],
        'JSON_S7_Lot_OCR_Ratio_and_RNLS_Gate_0D': decisions['S4_Metabolic_Demand_and_Autoreactive_Shield_Gate_0D'],
        'JSON_S8_SC_Islet_Maturity_Gate_0M': decisions['S8_SC_Islet_Paracrine_Maturity_Gate_0M'],
        'JSON_S8b_TYK2_Inducible_Promoter_Epistasis': decisions['S8b_TYK2_Inducible_Promoter_Epistasis_Gate_0E'],
        'JSON_S9_Hard_Safety_Stop': safety_status,
        'JSON_S10_CMV_NKG2C_Adaptive_NK_Selectivity_Gate_0A': decisions['S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A'],
        'JSON_S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate': decisions['S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate'],
    }
    return {
        'schema_version': tree_doc['schema_version'],
        'strict_precedence_order': tree_doc['strict_precedence_order'],
        'mandatory_s9_safety_panel_complete': mandatory_s9_safety_panel_complete,
        'full_5endpoint_s9_panel_complete': full_5endpoint_s9_panel_complete,
        's9_5endpoint_er_safety_panel_complete': full_5endpoint_s9_panel_complete,
        's9_causal_failed_categories': s9_causal_failed_categories,
        's9_panel_completeness_tier': s9_panel_completeness_tier,
        's9_frequency_assay_holds': [
            guard['status'] for guard in (
                _vaf_guard if has_vaf_obs else {},
                _cast_guard if has_cast_obs else {},
            ) if guard.get('status', '').startswith('HOLD_')
        ],
        'lot_release_readiness_tier': lot_release_readiness_tier,
        'safety_kill_s9_triggered': safety_kill_s9_triggered,
        'sticky_cross_batch_s9_stop_active': sticky_prior_s9_stop,
        'primary_biological_or_pouch_gate_failed': primary_biological_or_pouch_gate_failed,
        'zero_event_proportion_missing_denominator_or_explicit_se': zero_event_proportion_missing_denominator_or_explicit_se,
        'has_explicit_original_primary_obs': has_explicit_original_primary_obs,
        'has_explicit_release_critical_molecule_and_rate_se': has_explicit_release_critical_molecule_and_rate_se,
        'missing_original_primary_endpoints': missing_original_primary_endpoints,
        'any_fallback_switch_triggered': any_fallback_switch_triggered,
        'safety_kill_s9_status': safety_status,
        'resolved_architecture_assignments': resolved_arch,
        'validated_stage0_endpoints': validated_stage0,
        'validated_stage0_se_endpoints': validated_stage0_se,
        'stage0_se_provenance': se_provenance_by_mean,
        's9_point_estimate_fallback_se_policy': 'SCALE_AWARE_DEFAULT_SE_IMPUTED_ON_SCALAR_OMISSION_TYPED_ATOM_REQUIRES_MEASURED_SE',
        's9_se_imputed_endpoints': sorted(
            mk for mk in (
                'obs_icasp9_log10_kill',
                'obs_tp53_ppm1d_bcor_vaf',
                'obs_castseq_translocation_freq',
                'obs_chop_coincidence_fold',
                'obs_gate0e_chop_log2fc',
                'obs_hc10_surface_fold_wt',
                'obs_gate0e_hc10_surface_mfi_ratio',
            )
            if se_provenance_by_mean.get(mk) == 'DEFAULT_SE_IMPUTED_NOT_EXPLICIT'
        ),
        'switch_decisions': decisions,
        'canonical_json_switch_id_map': canonical_json_switch_id_map,
        'pivotal_combined_stress_coupon_gate': combined_stress_spec,
        **(
            {'iip0026_molecular_provenance_evaluation': iip0026_eval_summary}
            if iip0026_eval_summary is not None
            else {}
        ),
        **(
            {'iip0028_postdepletion_adcc_mixture_and_iapp_evaluation': iip0028_eval_summary}
            if iip0028_eval_summary is not None
            else {}
        ),
        **(
            {'iip0027_promoter_origin_and_itr_evaluation': iip0027_eval_summary}
            if iip0027_eval_summary is not None
            else {}
        ),
        **(
            {'iip0029_twoclock_maturation_and_repeated_challenge_evaluation': iip0029_eval_summary}
            if iip0029_eval_summary is not None
            else {}
        ),
        **(
            {
                'iip0031_locus2_and_tier6b_evaluation': {
                    's11_evaluation': iip0031_s11_eval,
                    **({'tier6b_sorted_lot_evaluation': iip0031_tier6b_eval} if iip0031_tier6b_eval is not None else {}),
                }
            }
            if (iip0031_s11_eval.get('status') != 'NOT_APPLICABLE_NON_LOCUS2_MODEB_VARIANT' or iip0031_tier6b_eval is not None)
            else {}
        ),
    }


# ============================================================================
# Stage-0 Named Atomic Provenance Scopes, ENDPOINT_REGISTRY_SPEC, Exact
# Student-t / 6+6 Welch TOST, & Conjunctive DMR / XCI / mtDNA / MDM4 Validators
# (IIP-0024; DEC-103)
# ============================================================================
import dataclasses
import hashlib
from types import MappingProxyType
from typing import Any


def _is_strict_int_iip0024(val: Any) -> bool:
    return isinstance(val, int) and not isinstance(val, bool)


def _is_strict_finite_number_iip0024(val: Any) -> bool:
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        return False
    try:
        return math.isfinite(val)
    except OverflowError:
        return False


def _deep_freeze_iip0024(obj: Any) -> Any:
    if isinstance(obj, (dict, MappingProxyType)):
        for k in obj.keys():
            if not isinstance(k, str):
                raise TypeError(f'AtomicObservation mapping keys must be strings, got {type(k).__name__} (IIP-0024; DEC-103)')
        return MappingProxyType({k: _deep_freeze_iip0024(v) for k, v in sorted(obj.items(), key=lambda kv: kv[0])})
    if isinstance(obj, (list, tuple)):
        return tuple(_deep_freeze_iip0024(x) for x in obj)
    if isinstance(obj, (set, frozenset)):
        return frozenset(_deep_freeze_iip0024(x) for x in obj)
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    raise TypeError(f'Unsupported mutable type inside AtomicObservation: {type(obj).__name__} (IIP-0024; DEC-103)')


def _deep_to_jsonable_iip0024(obj: Any) -> Any:
    if isinstance(obj, MappingProxyType) or isinstance(obj, dict):
        return {str(k): _deep_to_jsonable_iip0024(obj[k]) for k in sorted(obj.keys(), key=str)}
    if isinstance(obj, (list, tuple)):
        return [_deep_to_jsonable_iip0024(x) for x in obj]
    if isinstance(obj, frozenset):
        return sorted(_deep_to_jsonable_iip0024(x) for x in obj)
    return obj


class MCBScopeKey:
    """4-key Master Cell Bank provenance scope tuple (IIP-0024; DEC-103)."""
    __slots__ = ('product_variant_id', 'construct_sha256', 'donor_seedline_id', 'mcb_id')

    def __init__(self, product_variant_id: str, construct_sha256: str, donor_seedline_id: str, mcb_id: str) -> None:
        for fname, val in (
            ('product_variant_id', product_variant_id),
            ('construct_sha256', construct_sha256),
            ('donor_seedline_id', donor_seedline_id),
            ('mcb_id', mcb_id),
        ):
            if not isinstance(val, str) or not val.strip() or val != val.strip():
                raise ValueError(f'MCBScopeKey.{fname} must be a canonical non-whitespace-padded string (IIP-0024; DEC-103)')
            object.__setattr__(self, fname, val)

    def __setattr__(self, name: str, value: Any) -> None:
        raise TypeError('MCBScopeKey is immutable (IIP-0024; DEC-103)')

    def __delattr__(self, name: str) -> None:
        raise TypeError('MCBScopeKey is immutable (IIP-0024; DEC-103)')

    def as_tuple(self) -> tuple[str, str, str, str]:
        return (self.product_variant_id, self.construct_sha256, self.donor_seedline_id, self.mcb_id)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, MCBScopeKey) and self.as_tuple() == other.as_tuple()

    def __hash__(self) -> int:
        return hash(('MCBScopeKey', self.as_tuple()))

    def __repr__(self) -> str:
        return f'MCBScopeKey{self.as_tuple()!r}'


class LotScopeKey:
    """5-key Stage-7 differentiation lot provenance scope tuple (IIP-0024; DEC-103)."""
    __slots__ = ('product_variant_id', 'construct_sha256', 'donor_seedline_id', 'mcb_id', 'lot_id')

    def __init__(
        self,
        product_variant_id: str,
        construct_sha256: str,
        donor_seedline_id: str,
        mcb_id: str,
        lot_id: str,
    ) -> None:
        for fname, val in (
            ('product_variant_id', product_variant_id),
            ('construct_sha256', construct_sha256),
            ('donor_seedline_id', donor_seedline_id),
            ('mcb_id', mcb_id),
            ('lot_id', lot_id),
        ):
            if not isinstance(val, str) or not val.strip() or val != val.strip():
                raise ValueError(f'LotScopeKey.{fname} must be a canonical non-whitespace-padded string (IIP-0024; DEC-103)')
            object.__setattr__(self, fname, val)

    def __setattr__(self, name: str, value: Any) -> None:
        raise TypeError('LotScopeKey is immutable (IIP-0024; DEC-103)')

    def __delattr__(self, name: str) -> None:
        raise TypeError('LotScopeKey is immutable (IIP-0024; DEC-103)')

    def as_tuple(self) -> tuple[str, str, str, str, str]:
        return (self.product_variant_id, self.construct_sha256, self.donor_seedline_id, self.mcb_id, self.lot_id)

    def to_mcb_key(self) -> MCBScopeKey:
        return MCBScopeKey(self.product_variant_id, self.construct_sha256, self.donor_seedline_id, self.mcb_id)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, LotScopeKey) and self.as_tuple() == other.as_tuple()

    def __hash__(self) -> int:
        return hash(('LotScopeKey', self.as_tuple()))

    def __repr__(self) -> str:
        return f'LotScopeKey{self.as_tuple()!r}'


class AssembledDeviceLotScopeKey:
    """7-key assembled device x lot x physical route provenance scope tuple (IIP-0024; DEC-103)."""
    __slots__ = (
        'product_variant_id',
        'construct_sha256',
        'donor_seedline_id',
        'mcb_id',
        'lot_id',
        'device_assembly_id',
        'physical_route_id',
    )

    def __init__(
        self,
        product_variant_id: str,
        construct_sha256: str,
        donor_seedline_id: str,
        mcb_id: str,
        lot_id: str,
        device_assembly_id: str,
        physical_route_id: str,
    ) -> None:
        for fname, val in (
            ('product_variant_id', product_variant_id),
            ('construct_sha256', construct_sha256),
            ('donor_seedline_id', donor_seedline_id),
            ('mcb_id', mcb_id),
            ('lot_id', lot_id),
            ('device_assembly_id', device_assembly_id),
            ('physical_route_id', physical_route_id),
        ):
            if not isinstance(val, str) or not val.strip() or val != val.strip():
                raise ValueError(f'AssembledDeviceLotScopeKey.{fname} must be a canonical non-whitespace-padded string (IIP-0024; DEC-103)')
            object.__setattr__(self, fname, val)

    def __setattr__(self, name: str, value: Any) -> None:
        raise TypeError('AssembledDeviceLotScopeKey is immutable (IIP-0024; DEC-103)')

    def __delattr__(self, name: str) -> None:
        raise TypeError('AssembledDeviceLotScopeKey is immutable (IIP-0024; DEC-103)')

    def as_tuple(self) -> tuple[str, str, str, str, str, str, str]:
        return (
            self.product_variant_id,
            self.construct_sha256,
            self.donor_seedline_id,
            self.mcb_id,
            self.lot_id,
            self.device_assembly_id,
            self.physical_route_id,
        )

    def to_lot_key(self) -> LotScopeKey:
        return LotScopeKey(
            self.product_variant_id,
            self.construct_sha256,
            self.donor_seedline_id,
            self.mcb_id,
            self.lot_id,
        )

    def to_mcb_key(self) -> MCBScopeKey:
        return MCBScopeKey(self.product_variant_id, self.construct_sha256, self.donor_seedline_id, self.mcb_id)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, AssembledDeviceLotScopeKey) and self.as_tuple() == other.as_tuple()

    def __hash__(self) -> int:
        return hash(('AssembledDeviceLotScopeKey', self.as_tuple()))

    def __repr__(self) -> str:
        return f'AssembledDeviceLotScopeKey{self.as_tuple()!r}'


ENDPOINT_REGISTRY_SPEC: MappingProxyType = _deep_freeze_iip0024({
    'registry_version': 'IIP-0024-v1.0',
    'decision_id': 'DEC-103',
    'mcb_mandatory_endpoint_ids': (
        'mcb_donor_pairwise_ibd_pi_hat',
        'mcb_wgs_30x_structural_haplotype_verified',
        'mcb_karyotype_gband_ogm_verified',
        'mcb_imprinting_ic1_record',
        'mcb_imprinting_ic2_record',
        'mcb_imprinting_ig_dmr_record',
        'mcb_xci_46xx_3axis_record',
        'mcb_mtdna_1000x_scoped_record',
        'mcb_1q32_mdm4_diploid_cn',
        'mcb_1q32_mdm4_single_cell_mosaic',
        'mcb_20q11_bcl2l1_diploid_cn',
        'mcb_20q11_bcl2l1_single_cell_mosaic',
        'mcb_tier6b_pure_hpsc_clonogenic_escape',
    ),
    'lot_mandatory_endpoint_ids': (
        'lot_tier6b_sorted_prolif_clonogenic_escape',
    ),
    'device_sister_terminal_endpoint_ids': (
        'obs_explant_en_bloc_recovery_fraction',
        'gate0s_shedding_ddpcr_assay',
        'gate0s_shedding_viable_outgrowth_assay',
        'obs_iapp_d84_immunocapture_oligomer_elisa_fold',
        'obs_iapp_d84_colocalized_fibril_area_pct',
    ),
    'count_endpoint_caps_and_zero_vetoes': {
        'mcb_1q32_mdm4_single_cell_mosaic': {'cap': 0.01, 'require_observed_zero': False},
        'mcb_20q11_bcl2l1_single_cell_mosaic': {'cap': 0.01, 'require_observed_zero': False},
        'mcb_tier6b_pure_hpsc_clonogenic_escape': {'cap': 1e-8, 'require_observed_zero': True, 'use_poisson_protocol_floor': True, 'forbid_implicit_1p0_defaults': True},
        'obs_tp53_ppm1d_bcor_vaf': {'cap': 1.0e-4, 'require_observed_zero': False, 'use_poisson_protocol_floor': True},
        'obs_castseq_translocation_freq': {'cap': 1.0e-4, 'require_observed_zero': False, 'use_poisson_protocol_floor': True},
        'obs_20q11_bcl2l1_single_cell_mosaic': {'cap': 0.01, 'require_observed_zero': False},
        'obs_1q32_mdm4_single_cell_mosaic': {'cap': 0.01, 'require_observed_zero': False},
        'obs_assembled_cluster_ki67_progenitor': {'cap': 0.005, 'require_observed_zero': False},
        'obs_ductal_sox9_krt19_double_pos': {'cap': 0.001, 'require_observed_zero': False},
        'obs_assembled_cluster_slc18a1_enterochromaffin': {'cap': 0.01, 'require_observed_zero': False},
        'obs_senescent_p16_sasp_positive': {'cap': 0.01, 'require_observed_zero': False},
        'obs_gate0v_day28_persistent_vp1_ie1': {'cap': 0.0005, 'require_observed_zero': True},
        'lot_tier6b_sorted_prolif_clonogenic_escape': {
            'cap': 3.0e-7,
            'require_observed_zero': True,
            'requires_per_fraction_thinning': True,
            'forbid_implicit_1p0_defaults': True,
            'detection_floor_d_l_min': 0.50,
            'scope_key_tuple': ('endpoint_id', 'scope_key', 'fraction_id'),
        },
        'obs_day28_sox2_oct4_nanog_residual_psc': {'cap': 0.0005, 'require_observed_zero': False},
        'obs_day28_mki67_top2a_prolif_escapee': {'cap': 0.0005, 'require_observed_zero': False},
        'gate0s_shedding_ddpcr_assay': {'cap': 1e-6, 'require_observed_zero': True},
        'gate0s_shedding_viable_outgrowth_assay': {'cap': 1e-6, 'require_observed_zero': True},
    },
    'physical_routes': {
        'SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH': {
            'requires_cao2_catalase_p05': True,
            'requires_design_z_two_stage': False,
            'requires_cuff_panel': True,
            'canonical_design_z_divider_arm': None,
        },
        'Z_OMENTAL_140UM_ZERO_CAO2': {
            'requires_cao2_catalase_p05': False,
            'requires_design_z_two_stage': True,
            'requires_cuff_panel': True,
            'canonical_design_z_divider_arm': 'PLAIN_HYDROGEL_140UM',
        },
        'DESIGN_Z_140UM_ZERO_CAO2_PLAIN_HYDROGEL': {
            'requires_cao2_catalase_p05': False,
            'requires_design_z_two_stage': True,
            'requires_cuff_panel': True,
            'canonical_design_z_divider_arm': 'PLAIN_HYDROGEL_140UM',
        },
        'DESIGN_Z_140UM_ZERO_CAO2_EXPLORATORY_DIVIDER': {
            'requires_cao2_catalase_p05': False,
            'requires_design_z_two_stage': True,
            'requires_cuff_panel': True,
            'canonical_design_z_divider_arm': 'EXPLORATORY_18PCT_DIVIDER_140UM',
        },
    },
})

REGISTERED_VALIDATED_CLUSTER_ESTIMATORS_IIP0024: frozenset[str] = frozenset({
    'BETA_BINOMIAL_EXACT_PROFILE_L95_V1',
    'DESIGN_BASED_HYPERGEOMETRIC_CLUSTER_BOUND_V1',
})

_ATOMIC_ENDPOINT_KINDS = {
    **dict.fromkeys(ENDPOINT_REGISTRY_SPEC['count_endpoint_caps_and_zero_vetoes'], 'EXACT_INTEGER_COUNT'),
    **dict.fromkeys(S9_CORE_LOT_ENDPOINT_IDS, 'CONTINUOUS_STUDENT_T'),
    **dict.fromkeys((
        'stage7_stored_mature_insulin_iu', 'gate0s_ap1903_1h_insulin_burst_iu',
        'gate0s_ap1903_24h_insulin_dump_fraction', 'obs_explant_en_bloc_recovery_fraction',
        'mcb_1q32_mdm4_diploid_cn', 'mcb_20q11_bcl2l1_diploid_cn',
        'mcb_donor_pairwise_ibd_pi_hat', 'obs_iapp_d84_immunocapture_oligomer_elisa_fold',
        'obs_iapp_d84_colocalized_fibril_area_pct',
    ), 'CONTINUOUS_STUDENT_T'),
    **dict.fromkeys((
        'mcb_wgs_30x_structural_haplotype_verified', 'mcb_karyotype_gband_ogm_verified',
    ), 'METADATA_CONTROL'),
    **dict.fromkeys((
        'mcb_imprinting_ic1_record', 'mcb_imprinting_ic2_record', 'mcb_imprinting_ig_dmr_record',
        'mcb_xci_46xx_3axis_record', 'mcb_mtdna_1000x_scoped_record',
        'stage0_twoclock_repeated_challenge_panel', 'stage0_locus2_modeb_v3_panel',
    ), 'STRUCTURED_PANEL'),
}

_ATOMIC_REQUIRED_SCOPES = {
    **dict.fromkeys(ENDPOINT_REGISTRY_SPEC['mcb_mandatory_endpoint_ids'], 'MCB_SCOPE'),
    **dict.fromkeys(ENDPOINT_REGISTRY_SPEC['device_sister_terminal_endpoint_ids'], 'ASSEMBLED_DEVICE_LOT_SCOPE'),
    **dict.fromkeys(('gate0s_ap1903_1h_insulin_burst_iu', 'gate0s_ap1903_24h_insulin_dump_fraction'), 'ASSEMBLED_DEVICE_LOT_SCOPE'),
    **dict.fromkeys((*S9_CORE_LOT_ENDPOINT_IDS, 'stage7_stored_mature_insulin_iu'), 'LOT_SCOPE'),
}


def _atomic_mcb_panel_status(validator, scope: str, payload: Any) -> str:
    if not isinstance(payload, dict):
        return 'HOLD_MALFORMED_STRUCTURED_PANEL'
    try:
        args = inspect.signature(validator).bind(endpoint_scope=scope, **payload)
    except TypeError:
        return 'HOLD_MALFORMED_STRUCTURED_PANEL'
    status = validator(*args.args, **args.kwargs)['status']
    return 'PASS' if status == 'NOT_APPLICABLE_46XY' else status


class AtomicObservation:
    """
    Deep-immutable Stage-0 evidence atom (MappingProxyType-frozen) derived strictly from
    ENDPOINT_REGISTRY_SPEC (forbidding caller-supplied canonical_projection or cap overrides;
    IIP-0024; DEC-103).
    """
    __slots__ = (
        'endpoint_id',
        'scope_type',
        'scope_key',
        'kind',
        'value_or_k',
        'se_or_n',
        'sample_n',
        'unit',
        'controls',
        'assay_version',
        'campaign_id',
        '_campaign_explicit',
        'timepoint',
        'sample_unit_id',
        'atom_sha256',
    )

    def __init__(
        self,
        *,
        endpoint_id: str,
        scope_type: str,
        scope_key: MCBScopeKey | LotScopeKey | AssembledDeviceLotScopeKey,
        kind: str,
        value_or_k: Any,
        se_or_n: Any = None,
        sample_n: int | None = None,
        unit: str = 'canonical_unit',
        controls: dict[str, Any] | MappingProxyType | None = None,
        assay_version: str = 'STAGE0_IIP0024_V1',
        campaign_id: str | None = None,
        timepoint: str = 'STAGE0',
        sample_unit_id: str = 'SAMPLE-01',
        **forbidden_kwargs: Any,
    ) -> None:
        if forbidden_kwargs:
            raise ValueError(
                f'Caller-supplied projection or threshold overrides are prohibited on AtomicObservation: '
                f'{sorted(forbidden_kwargs.keys())} (IIP-0024; DEC-103)'
            )
        if not isinstance(endpoint_id, str) or not endpoint_id.strip() or endpoint_id != endpoint_id.strip():
            raise ValueError('endpoint_id must be a canonical non-whitespace-padded nonempty string (IIP-0024; DEC-103)')
        if scope_type == 'MCB_SCOPE' and not isinstance(scope_key, MCBScopeKey):
            raise ValueError('MCB_SCOPE requires an MCBScopeKey instance (IIP-0024; DEC-103)')
        if scope_type == 'LOT_SCOPE' and not isinstance(scope_key, LotScopeKey):
            raise ValueError('LOT_SCOPE requires a LotScopeKey instance (IIP-0024; DEC-103)')
        if scope_type == 'ASSEMBLED_DEVICE_LOT_SCOPE' and not isinstance(scope_key, AssembledDeviceLotScopeKey):
            raise ValueError('ASSEMBLED_DEVICE_LOT_SCOPE requires an AssembledDeviceLotScopeKey instance (IIP-0024; DEC-103)')
        if scope_type not in ('MCB_SCOPE', 'LOT_SCOPE', 'ASSEMBLED_DEVICE_LOT_SCOPE'):
            raise ValueError(f'Unknown scope_type {scope_type!r} (IIP-0024; DEC-103)')
        if kind not in ('CONTINUOUS_STUDENT_T', 'SCALAR_EXACT', 'EXACT_INTEGER_COUNT', 'STRUCTURED_PANEL', 'METADATA_CONTROL'):
            raise ValueError(f'Unknown AtomicObservation kind {kind!r} (IIP-0024; DEC-103)')
        if kind == 'EXACT_INTEGER_COUNT':
            if not _is_strict_int_iip0024(value_or_k) or not _is_strict_int_iip0024(se_or_n):
                raise ValueError('EXACT_INTEGER_COUNT requires strict Python integers (k, n) pre-normalization (IIP-0024; DEC-103)')
            if se_or_n <= 0 or value_or_k < 0 or value_or_k > se_or_n:
                raise ValueError('EXACT_INTEGER_COUNT requires 0 <= k <= n and n > 0 (IIP-0024; DEC-103)')
        if kind == 'CONTINUOUS_STUDENT_T':
            if (
                not _is_strict_finite_number_iip0024(value_or_k)
                or not _is_strict_finite_number_iip0024(se_or_n)
                or float(se_or_n) <= 0.0
                or not _is_strict_int_iip0024(sample_n)
                or sample_n < 6
            ):
                raise ValueError('CONTINUOUS_STUDENT_T requires finite mean, SE > 0, and integer sample_n >= 6 (IIP-0024; DEC-103)')
        if controls is not None and not isinstance(controls, (dict, MappingProxyType)):
            raise ValueError('AtomicObservation controls must be a dict or MappingProxyType (IIP-0024; DEC-103)')
        ctrl_dict = dict(controls) if controls is not None else {}
        for bad_ck in ('canonical_projection', 'canonical_projection_keys', 'override_cap', 'override_requires_k_eq_0', 'cap'):
            if bad_ck in ctrl_dict:
                raise ValueError(f'Caller-supplied {bad_ck!r} inside controls is prohibited; caps derive from ENDPOINT_REGISTRY_SPEC (IIP-0024; DEC-103)')
        frozen_val = _deep_freeze_iip0024(value_or_k)
        frozen_ctrl = _deep_freeze_iip0024(ctrl_dict)
        campaign_explicit = campaign_id is not None
        eff_campaign_id = 'CAMP-DEFAULT-01' if campaign_id is None else campaign_id
        for s_name, s_val in (
            ('unit', unit),
            ('assay_version', assay_version),
            ('campaign_id', eff_campaign_id),
            ('timepoint', timepoint),
            ('sample_unit_id', sample_unit_id),
        ):
            if not isinstance(s_val, str) or not s_val.strip() or s_val != s_val.strip():
                raise ValueError(f'{s_name} must be a canonical non-whitespace-padded nonempty string (IIP-0024; DEC-103)')
        canonical_doc = {
            'registry_version': ENDPOINT_REGISTRY_SPEC['registry_version'],
            'endpoint_id': endpoint_id,
            'scope_type': scope_type,
            'scope_key': list(scope_key.as_tuple()),
            'kind': kind,
            'value_or_k': _deep_to_jsonable_iip0024(frozen_val),
            'se_or_n': se_or_n,
            'sample_n': sample_n,
            'unit': unit,
            'controls': _deep_to_jsonable_iip0024(frozen_ctrl),
            'assay_version': assay_version,
            'campaign_id': eff_campaign_id,
            'timepoint': timepoint,
            'sample_unit_id': sample_unit_id,
        }
        digest = hashlib.sha256(
            json.dumps(canonical_doc, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')
        ).hexdigest()
        object.__setattr__(self, 'endpoint_id', endpoint_id)
        object.__setattr__(self, 'scope_type', scope_type)
        object.__setattr__(self, 'scope_key', scope_key)
        object.__setattr__(self, 'kind', kind)
        object.__setattr__(self, 'value_or_k', frozen_val)
        object.__setattr__(self, 'se_or_n', se_or_n)
        object.__setattr__(self, 'sample_n', sample_n)
        object.__setattr__(self, 'unit', unit)
        object.__setattr__(self, 'controls', frozen_ctrl)
        object.__setattr__(self, 'assay_version', assay_version)
        object.__setattr__(self, 'campaign_id', eff_campaign_id)
        object.__setattr__(self, '_campaign_explicit', campaign_explicit)
        object.__setattr__(self, 'timepoint', timepoint)
        object.__setattr__(self, 'sample_unit_id', sample_unit_id)
        object.__setattr__(self, 'atom_sha256', digest)

    def __setattr__(self, name: str, value: Any) -> None:
        raise TypeError('AtomicObservation is strictly immutable (IIP-0024; DEC-103)')

    def __delattr__(self, name: str) -> None:
        raise TypeError('AtomicObservation is strictly immutable (IIP-0024; DEC-103)')

    def digest(self) -> str:
        return self.atom_sha256


def evaluate_rare_event_count_iip0024(
    *,
    gate_id: str,
    k: Any,
    n: Any,
    max_u95_fraction: float,
    require_observed_zero: bool = False,
    reported_rate: Any = None,
    rate_scale: float = 1.0,
    sampling_independence_model: str = 'INDEPENDENT_BIOLOGICAL_OPPORTUNITIES',
    clustered_sampling: bool = False,
    validated_clustered_estimator_id: str | None = None,
    eta_spike_recovery_l95: Any = None,
    f_clonogenic: Any = None,
    use_poisson_protocol_floor: bool = False,
    require_explicit_calibration: bool = False,
) -> dict[str, Any]:
    """
    Fail-closed rare-event integer count validator (IIP-0024; DEC-103; hardened IIP-0031; DEC-110):
      - Rejects bool, float (e.g. 298.9 or 300.0), NaN, Inf, < 0, k > n with structured HOLD.
      - Enforces k == 0 observed-zero hard veto BEFORE calibration/recovery holds so an attributable
        integer k >= 1 is never downgraded to a calibration hold.
      - Forbids implicit 1.0 defaults on missing recovery / f_clonogenic for Tier 6B / clonogenic
        or explicit-calibration gates.
      - Computes exact independent binomial-thinning sensitivity U95 = (1 - 0.05**(1/n)) / (f_clono * eta)
        (holding 0/920 at eta = 0.325 where U95 = 0.01000288 > 0.01, never floor(n * eta)).
      - Rejects unregistered cluster estimators and naive DEFF substitution.
    """
    is_tier6b_gate = ('TIER6B' in str(gate_id).upper()) or ('CLONOGENIC' in str(gate_id).upper())
    # Latch attributable integer k >= 1 BEFORE missing/invalid n, calibration holds, recovery bounds, or clustered checks (IIP-0031; DEC-110):
    if require_observed_zero and _is_strict_int_iip0024(k) and k > 0:
        u95_nom = (
            round(clopper_pearson_upper95_independent(k, n), 10)
            if (_is_strict_int_iip0024(n) and n >= k and n > 0)
            else None
        )
        return {
            'status': 'FAIL_HARD_STOP',
            'reason': 'FAIL_TIER6B_OUTGROWTH_DETECTED' if is_tier6b_gate else 'OBSERVED_POSITIVE_EVENT_VIOLATES_MANDATORY_ZERO_VETO',
            'gate_id': gate_id,
            'k': k,
            'n_nominal': n if _is_strict_int_iip0024(n) else None,
            'u95_nominal': u95_nom,
        }
    if k is None or n is None:
        return {
            'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_MISSING_ATOMIC_INTEGER_K_OR_N',
            'gate_id': gate_id,
        }
    if not _is_strict_int_iip0024(k) or not _is_strict_int_iip0024(n):
        return {
            'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_NON_INTEGER_COUNT_OR_DENOMINATOR_REJECTED',
            'gate_id': gate_id,
        }
    if n <= 0 or k < 0 or k > n:
        return {
            'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_INVALID_COUNT_BOUNDS',
            'gate_id': gate_id,
        }
    if reported_rate is not None:
        if not _is_strict_finite_number_iip0024(reported_rate):
            return {'status': 'HOLD_NONFINITE_REPORTED_RATE', 'gate_id': gate_id}
        expected_rate = rate_scale * (k / float(n))
        if abs(float(reported_rate) - expected_rate) > 1e-9:
            return {'status': 'HOLD_RATE_COUNT_INCONSISTENCY_NEVER_INFER_FROM_ROUNDED_RATE', 'gate_id': gate_id}
    if clustered_sampling or sampling_independence_model != 'INDEPENDENT_BIOLOGICAL_OPPORTUNITIES':
        if (
            not isinstance(validated_clustered_estimator_id, str)
            or validated_clustered_estimator_id.strip() not in REGISTERED_VALIDATED_CLUSTER_ESTIMATORS_IIP0024
        ):
            return {
                'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_CLUSTERED_OR_UNREGISTERED_ESTIMATOR_REJECTED',
                'gate_id': gate_id,
            }
    if is_tier6b_gate or require_explicit_calibration:
        if eta_spike_recovery_l95 is None or f_clonogenic is None:
            return {
                'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_MISSING_RECOVERY_OR_CLONOGENIC_FACTOR',
                'reason': 'MISSING_RECOVERY_OR_CLONOGENIC_FACTOR_NO_IMPLICIT_1P0_DEFAULT',
                'gate_id': gate_id,
            }
    else:
        if eta_spike_recovery_l95 is None:
            eta_spike_recovery_l95 = 1.0
        if f_clonogenic is None:
            f_clonogenic = 1.0
    if not _is_strict_finite_number_iip0024(eta_spike_recovery_l95) or not _is_strict_finite_number_iip0024(f_clonogenic):
        return {
            'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_NONFINITE_RECOVERY_OR_CLONOGENIC_FACTOR',
            'gate_id': gate_id,
        }
    eta = float(eta_spike_recovery_l95)
    f_clono = float(f_clonogenic)
    if not (0.0 < eta <= 1.0) or not (0.0 < f_clono <= 1.0):
        return {
            'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION' if is_tier6b_gate else 'HOLD_OUT_OF_RANGE_RECOVERY_OR_CLONOGENIC_FACTOR',
            'gate_id': gate_id,
        }
    if is_tier6b_gate and (eta < 0.50 or (f_clono * eta) < 0.50):
        return {'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION', 'reason': 'TIER6B_DETECTION_BELOW_0P50_FLOOR', 'gate_id': gate_id}
    if (eta < 1.0 or f_clono < 1.0) and eta < 0.30:
        return {'status': 'HOLD_SPIKE_RECOVERY_BELOW_0P30_QUALIFICATION_FLOOR', 'gate_id': gate_id}
    eff_prob = f_clono * eta
    if eff_prob < 1.0 and k > 0:
        return {
            'status': 'HOLD_NONZERO_RECOVERY_ADJUSTED_COUNT_REQUIRES_VALIDATED_JOINT_ESTIMATOR',
            'gate_id': gate_id,
        }
    if k == 0:
        # Exact independent per-trial detection probability P(obs=0) = (1 - p * f_clono * eta)**n <= 0.05:
        u95_bin = min(1.0, (-math.expm1(math.log(0.05) / float(n))) / eff_prob)
        u95_pois = min(1.0, (-math.log(0.05) / (float(n) * eff_prob)))
    else:
        u95_bin = clopper_pearson_upper95_independent(k, n)
        u95_pois = u95_bin
    u95_eval = u95_pois if use_poisson_protocol_floor else u95_bin
    exceeds = (u95_eval >= max_u95_fraction) if gate_id.startswith('GATE0S_SHEDDING_') else (u95_eval > max_u95_fraction)
    if exceeds:
        if k == 0:
            return {
                'status': 'HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR',
                'gate_id': gate_id,
                'k': 0,
                'n_nominal': n,
                'effective_detection_factor': eff_prob,
                'u95_binomial': u95_bin,
                'u95_poisson': u95_pois,
                'cap_fraction': max_u95_fraction,
            }
        return {
            'status': 'FAIL_HARD_STOP',
            'reason': 'CONFIRMED_COUNT_UPPER95_EXCEEDS_CAP',
            'gate_id': gate_id,
            'k': k,
            'n_nominal': n,
            'u95_binomial': u95_bin,
            'cap_fraction': max_u95_fraction,
        }
    return {
        'status': 'PASS',
        'gate_id': gate_id,
        'k': k,
        'n_nominal': n,
        'effective_detection_factor': eff_prob,
        'u95_binomial': u95_bin,
        'u95_poisson': u95_pois,
        'cap_fraction': max_u95_fraction,
    }


REGISTERED_DMR_PANEL_SPECS_IIP0024: MappingProxyType = _deep_freeze_iip0024({
    'IC1': {
        'panel_version_id': 'IC1_H19_CTCF_CORE_V1',
        'dmr_class': 'PRIMARY_GERMLINE_DMR',
        'expected_methylated_parent': 'PATERNAL',
        'grch38_window': 'chr11:1999500-2003500',
    },
    'IC2': {
        'panel_version_id': 'IC2_KVDMR1_PROMOTER_CORE_V1',
        'dmr_class': 'PRIMARY_GERMLINE_DMR',
        'expected_methylated_parent': 'MATERNAL',
        'grch38_window': 'chr11:2698800-2701199',
    },
    'IG_DMR': {
        'panel_version_id': 'IG_DMR_KAGAMI_CG4_6CPG_CG6_19CPG_V1',
        'dmr_class': 'PRIMARY_GERMLINE_DMR',
        'expected_methylated_parent': 'PATERNAL',
        'grch38_window': 'chr14:100809350-100811900',
        'amplicon_cpg_counts': {'CG4': 6, 'CG6': 19},
    },
    'MEG3_DMR': {
        'panel_version_id': 'MEG3_DMR_SOMATIC_STAGE7_V1',
        'dmr_class': 'SECONDARY_SOMATIC_DMR',
        'expected_methylated_parent': 'PATERNAL',
        'grch38_window': 'chr14:100824100-100828299',
    },
})


def evaluate_conjunctive_imprinting_dmr_iip0024(
    *,
    endpoint_scope: str,
    dmr_id: Any,
    panel_version_id: Any,
    per_cpg_callable_verified: Any,
    assay_controls_valid: Any,
    bulk_beta: Any,
    phased: Any,
    beta_paternal: Any = None,
    beta_maternal: Any = None,
    min_allele_molecule_depth: Any = None,
    mosaic_discordant_epiallele_frac: Any = None,
    stage7_meg3_dlk1_rtqpcr_and_perifusion_cleared: Any = None,
) -> dict[str, Any]:
    """NaN-safe conjunctive DMR validator for IC1, IC2/KvDMR1, IG-DMR (CG4/CG6), and Stage-7 MEG3-DMR diagnostic HOLD (IIP-0024; DEC-103)."""
    if endpoint_scope not in ('MCB_SCOPE', 'LOT_SCOPE'):
        return {'status': 'HOLD_INVALID_ENDPOINT_SCOPE'}
    if not isinstance(dmr_id, str) or dmr_id not in REGISTERED_DMR_PANEL_SPECS_IIP0024:
        return {'status': 'HOLD_UNKNOWN_DMR_ID', 'dmr_id': dmr_id}
    spec = REGISTERED_DMR_PANEL_SPECS_IIP0024[dmr_id]
    if panel_version_id != spec['panel_version_id'] or per_cpg_callable_verified is not True:
        return {'status': 'HOLD_UNQUALIFIED_OR_MISMATCHED_DMR_PANEL_VERSION', 'dmr_id': dmr_id}
    if assay_controls_valid is not True:
        return {'status': 'HOLD_INVALID_OR_MISSING_DMR_ASSAY_CONTROLS', 'dmr_id': dmr_id}
    if not _is_strict_finite_number_iip0024(bulk_beta) or not (0.0 <= float(bulk_beta) <= 1.0):
        return {'status': 'HOLD_NONFINITE_OR_OUT_OF_RANGE_BULK_BETA', 'dmr_id': dmr_id}
    if not isinstance(phased, bool):
        return {'status': 'HOLD_INVALID_PHASED_FLAG', 'dmr_id': dmr_id}
    b_bulk = float(bulk_beta)
    if dmr_id == 'MEG3_DMR':
        aberrant = not (0.35 <= b_bulk <= 0.65)
        if phased:
            if (
                not _is_strict_int_iip0024(min_allele_molecule_depth)
                or min_allele_molecule_depth < 30
                or not _is_strict_finite_number_iip0024(beta_paternal)
                or not _is_strict_finite_number_iip0024(beta_maternal)
                or not (0.0 <= float(beta_paternal) <= 1.0)
                or not (0.0 <= float(beta_maternal) <= 1.0)
            ):
                return {'status': 'HOLD_MEG3_DMR_INCOMPLETE_PHASED_METRICS', 'dmr_id': dmr_id}
            if (float(beta_paternal) - float(beta_maternal)) < 0.60:
                aberrant = True
        if aberrant:
            if endpoint_scope == 'LOT_SCOPE' and stage7_meg3_dlk1_rtqpcr_and_perifusion_cleared is False:
                return {'status': 'FAIL_HARD_STOP', 'dmr_id': dmr_id, 'reason': 'STAGE7_MEG3_SILENCING_AND_PERIFUSION_FAILED'}
            return {
                'status': 'HOLD_STAGE7_MEG3_DMR_REQUIRES_MEG3_DLK1_EXPRESSION_AND_PERIFUSION_CHECK',
                'dmr_id': dmr_id,
                'endpoint_scope': endpoint_scope,
            }
        return {'status': 'PASS', 'dmr_id': dmr_id, 'endpoint_scope': endpoint_scope}
    if b_bulk < 0.35 or b_bulk > 0.65:
        return {'status': 'FAIL_HARD_STOP', 'dmr_id': dmr_id, 'reason': 'BULK_METHYLATION_OUTSIDE_0P35_0P65'}
    if not phased or not _is_strict_int_iip0024(min_allele_molecule_depth) or min_allele_molecule_depth < 30:
        return {'status': 'HOLD_UNPHASED_OR_LOW_ALLELE_DEPTH', 'dmr_id': dmr_id}
    if (
        not _is_strict_finite_number_iip0024(beta_paternal)
        or not _is_strict_finite_number_iip0024(beta_maternal)
        or not (0.0 <= float(beta_paternal) <= 1.0)
        or not (0.0 <= float(beta_maternal) <= 1.0)
    ):
        return {'status': 'HOLD_NONFINITE_OR_MISSING_PARENTAL_ALLELE_BETA', 'dmr_id': dmr_id}
    if not _is_strict_finite_number_iip0024(mosaic_discordant_epiallele_frac) or not (0.0 <= float(mosaic_discordant_epiallele_frac) <= 1.0):
        return {'status': 'HOLD_MISSING_OR_INVALID_MOSAIC_EPIALLELE_FRACTION', 'dmr_id': dmr_id}
    b_pat, b_mat, mos = float(beta_paternal), float(beta_maternal), float(mosaic_discordant_epiallele_frac)
    b_meth = b_mat if spec['expected_methylated_parent'] == 'MATERNAL' else b_pat
    b_unmeth = b_pat if spec['expected_methylated_parent'] == 'MATERNAL' else b_mat
    signed_delta = b_meth - b_unmeth
    if signed_delta < 0.60 or b_meth < 0.75 or b_unmeth > 0.25 or mos > 0.10:
        return {
            'status': 'FAIL_HARD_STOP',
            'dmr_id': dmr_id,
            'signed_parental_delta': round(signed_delta, 4),
            'mosaic_discordant_epiallele_frac': mos,
        }
    return {'status': 'PASS', 'dmr_id': dmr_id, 'endpoint_scope': endpoint_scope, 'signed_parental_delta': round(signed_delta, 4)}


def evaluate_conjunctive_xci_46xx_iip0024(
    *,
    endpoint_scope: str,
    donor_sex_karyotype: Any,
    assay_controls_valid: Any,
    xci_promoter_panel_id: Any = None,
    xist_pos_k: Any = None,
    xist_nuclei_n: Any = None,
    xci_subject_promoter_cpg_retention: Any = None,
    allele_specific_x_rna_informative_genes: Any = None,
    allele_specific_x_rna_monoallelic_fraction: Any = None,
) -> dict[str, Any]:
    """Fail-closed 3-axis 46,XX XCI erosion validator (n >= 300 nuclei, k >= 279 -> L95 >= 0.90, promoter CpG >= 0.85, X-RNA >= 0.90 across >= 10 genes; IIP-0024; DEC-103)."""
    if endpoint_scope not in ('MCB_SCOPE', 'LOT_SCOPE'):
        return {'status': 'HOLD_INVALID_ENDPOINT_SCOPE'}
    if donor_sex_karyotype == '46,XY':
        return {'status': 'NOT_APPLICABLE_46XY', 'endpoint_scope': endpoint_scope}
    if donor_sex_karyotype != '46,XX':
        return {'status': 'HOLD_MISSING_OR_UNRECOGNIZED_KARYOTYPE'}
    if assay_controls_valid is not True:
        return {'status': 'HOLD_INVALID_OR_MISSING_XCI_ASSAY_CONTROLS'}
    if xci_promoter_panel_id != 'XCI_SUBJECT_PROMOTER_CPG_PANEL_V1':
        return {'status': 'HOLD_MISSING_OR_UNKNOWN_XCI_PROMOTER_PANEL_ID'}
    if not _is_strict_int_iip0024(xist_pos_k) or not _is_strict_int_iip0024(xist_nuclei_n) or xist_pos_k < 0 or xist_pos_k > xist_nuclei_n:
        return {'status': 'HOLD_NON_INTEGER_OR_OUT_OF_BOUNDS_XIST_COUNTS'}
    if xist_nuclei_n < 300:
        return {'status': 'HOLD_INSUFFICIENT_XIST_NUCLEI_REQUIRES_GE_300'}
    if not _is_strict_finite_number_iip0024(xci_subject_promoter_cpg_retention) or not (0.0 <= float(xci_subject_promoter_cpg_retention) <= 1.0):
        return {'status': 'HOLD_NONFINITE_OR_OUT_OF_RANGE_XCI_PROMOTER_RETENTION'}
    xist_l95 = clopper_pearson_lower95_independent(xist_pos_k, xist_nuclei_n)
    if xist_l95 < 0.90 or float(xci_subject_promoter_cpg_retention) < 0.85:
        return {'status': 'FAIL_HARD_STOP', 'xist_lower95': round(xist_l95, 6), 'promoter_cpg_retention': float(xci_subject_promoter_cpg_retention)}
    if not _is_strict_int_iip0024(allele_specific_x_rna_informative_genes) or allele_specific_x_rna_informative_genes < 10:
        return {'status': 'HOLD_INSUFFICIENT_INFORMATIVE_X_RNA_SNPS_REQUIRES_GE_10'}
    if not _is_strict_finite_number_iip0024(allele_specific_x_rna_monoallelic_fraction) or not (0.0 <= float(allele_specific_x_rna_monoallelic_fraction) <= 1.0):
        return {'status': 'HOLD_NONFINITE_OR_MISSING_X_RNA_MONOALLELIC_FRACTION'}
    if float(allele_specific_x_rna_monoallelic_fraction) < 0.90:
        return {'status': 'FAIL_HARD_STOP', 'reason': 'BIALLELIC_X_LINKED_RNA_REACTIVATION_DETECTED'}
    return {
        'status': 'PASS',
        'endpoint_scope': endpoint_scope,
        'xist_lower95': round(xist_l95, 6),
        'promoter_cpg_retention': float(xci_subject_promoter_cpg_retention),
        'monoallelic_rna_fraction': float(allele_specific_x_rna_monoallelic_fraction),
    }


def evaluate_scoped_mtdna_gate_iip0024(
    *,
    endpoint_scope: str,
    min_unique_chrm_molecule_depth: Any,
    numt_filter_and_snv_del_spikes_validated: Any,
    variants: Any,
    deletion_vaf_max: Any = None,
    seahorse_coupled_respiration_pass: Any = None,
) -> dict[str, Any]:
    """Fail-closed >= 1,000x NUMT-filtered chrM heteroplasmy, deletion, and respiration validator (IIP-0024; DEC-103)."""
    if endpoint_scope not in ('MCB_SCOPE', 'LOT_SCOPE'):
        return {'status': 'HOLD_INVALID_ENDPOINT_SCOPE'}
    if not _is_strict_int_iip0024(min_unique_chrm_molecule_depth) or min_unique_chrm_molecule_depth < 1000:
        return {'status': 'HOLD_CHRM_UNIQUE_DEPTH_BELOW_1000X'}
    if numt_filter_and_snv_del_spikes_validated is not True:
        return {'status': 'HOLD_NUMT_FILTER_OR_SNV_DEL_SPIKE_CONTROLS_UNVALIDATED'}
    if not _is_strict_finite_number_iip0024(deletion_vaf_max) or not (0.0 <= float(deletion_vaf_max) <= 1.0):
        return {'status': 'HOLD_MISSING_OR_INVALID_MTDNA_DELETION_VAF_CONTROL'}
    if not isinstance(seahorse_coupled_respiration_pass, bool):
        return {'status': 'HOLD_MISSING_MTDNA_SEAHORSE_RESPIRATION_CONTROL'}
    if not isinstance(variants, list):
        return {'status': 'HOLD_INVALID_MTDNA_VARIANTS_LIST'}
    if seahorse_coupled_respiration_pass is False or float(deletion_vaf_max) > 0.02:
        return {'status': 'FAIL_HARD_STOP', 'reason': 'MTDNA_DELETION_GT_2PCT_OR_RESPIRATION_FAILURE'}
    borderline_hold = False
    for v in variants:
        if not isinstance(v, dict):
            return {'status': 'HOLD_MALFORMED_MTDNA_VARIANT_RECORD'}
        vid = v.get('variant_id')
        vaf_raw = v.get('vaf')
        if not isinstance(vid, str) or not vid.strip() or not _is_strict_finite_number_iip0024(vaf_raw) or not (0.0 <= float(vaf_raw) <= 1.0):
            return {'status': 'HOLD_INVALID_MTDNA_VARIANT_ID_OR_VAF'}
        vaf = float(vaf_raw)
        b_hap = v.get('verified_benign_donor_germline_haplogroup_polymorphism')
        b_path = v.get('known_pathogenic_mitomap_clinvar')
        b_acq = v.get('acquired_denovo_or_private_coding_trna_rrna')
        if any(x is not None and not isinstance(x, bool) for x in (b_hap, b_path, b_acq)):
            return {'status': 'HOLD_NON_BOOLEAN_MTDNA_CLASSIFICATION_FLAG', 'variant_id': vid}
        is_hap = b_hap is True
        is_path = b_path is True
        is_acq = b_acq is True
        if is_hap and (is_path or is_acq):
            return {'status': 'HOLD_CONTRADICTORY_MTDNA_BENIGN_AND_ACQUIRED_OR_PATHOGENIC_FLAGS', 'variant_id': vid}
        if not (is_hap or is_path or is_acq):
            return {'status': 'HOLD_UNCLASSIFIED_MTDNA_VARIANT_REQUIRES_CURATION', 'variant_id': vid}
        if vid.strip() == 'm.3243A>G' and vaf > 0.01:
            return {'status': 'FAIL_HARD_STOP', 'reason': 'M3243AG_EXCEEDS_1PCT_CAP', 'vaf': vaf}
        if is_hap:
            continue
        if vaf > 0.05:
            return {'status': 'FAIL_HARD_STOP', 'reason': 'SCOPED_MTDNA_HETEROPLASMY_EXCEEDS_5PCT_CAP', 'variant_id': vid, 'vaf': vaf}
        if 0.03 <= vaf <= 0.05:
            borderline_hold = True
    if borderline_hold:
        return {'status': 'HOLD_BORDERLINE_3_TO_5PCT_ACQUIRED_MTDNA_REQUIRES_ORTHOGONAL_DDPCR'}
    return {'status': 'PASS', 'endpoint_scope': endpoint_scope}


def _norm_cdf_iip0024(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _chi2_pdf_iip0024(w: float, df: float) -> float:
    if w <= 0.0:
        return 0.0
    half_df = 0.5 * df
    return math.exp((half_df - 1.0) * math.log(w) - 0.5 * w - half_df * math.log(2.0) - math.lgamma(half_df))


def welch_tost_power_exact_2d_iip0024(
    m1: int,
    m2: int,
    sigma1: float,
    sigma2: float,
    delta_true: float,
    theta_margin: float,
    n_pts: int = 350,
) -> float:
    """Exact 2D chi-square quadrature over independent (S1^2, S2^2) for unpooled Welch 90% TOST power (IIP-0024; DEC-103)."""
    df1 = float(m1 - 1)
    df2 = float(m2 - 1)
    se_true = math.sqrt((sigma1 ** 2) / m1 + (sigma2 ** 2) / m2)
    u1_max = math.sqrt(df1 + 15.0 * math.sqrt(2.0 * df1))
    u2_max = math.sqrt(df2 + 15.0 * math.sqrt(2.0 * df2))
    h1 = u1_max / n_pts
    h2 = u2_max / n_pts
    u1_nodes = [
        ((i + 0.5) * h1, _chi2_pdf_iip0024(((i + 0.5) * h1) ** 2, df1) * 2.0 * ((i + 0.5) * h1) * h1)
        for i in range(n_pts)
    ]
    u2_nodes = [
        ((j + 0.5) * h2, _chi2_pdf_iip0024(((j + 0.5) * h2) ** 2, df2) * 2.0 * ((j + 0.5) * h2) * h2)
        for j in range(n_pts)
    ]
    c1 = (sigma1 ** 2) / (m1 * df1)
    c2 = (sigma2 ** 2) / (m2 * df2)
    df_lo = min(df1, df2)
    df_hi = df1 + df2
    n_tab = 4096
    df_step = (df_hi - df_lo) / n_tab
    t_tab = [student_t_ppf_one_sided_95(df_lo + idx * df_step) for idx in range(n_tab + 1)]
    inv_sqrt2 = 1.0 / math.sqrt(2.0)

    acc = 0.0
    for u1, wt1 in u1_nodes:
        a1 = c1 * (u1 * u1)
        a1_sq_df = (a1 * a1) / df1
        for u2, wt2 in u2_nodes:
            a2 = c2 * (u2 * u2)
            v_tot = a1 + a2
            se_w = math.sqrt(v_tot)
            df_w = (v_tot * v_tot) / (a1_sq_df + (a2 * a2) / df2)
            pos = (df_w - df_lo) / df_step
            idx = int(pos)
            if idx < 0:
                t_val = t_tab[0]
            elif idx >= n_tab:
                t_val = t_tab[n_tab]
            else:
                frac = pos - idx
                t_val = t_tab[idx] + frac * (t_tab[idx + 1] - t_tab[idx])
            half_w = t_val * se_w
            if half_w < theta_margin:
                hi_z = (theta_margin - half_w - delta_true) / se_true
                lo_z = (-theta_margin + half_w - delta_true) / se_true
                if hi_z > lo_z:
                    acc += 0.5 * (math.erf(hi_z * inv_sqrt2) - math.erf(lo_z * inv_sqrt2)) * wt1 * wt2
    return min(1.0, max(0.0, acc))


def evaluate_donor_blocked_contrast_and_welch_tost_iip0024(n_pts: int = 350) -> dict[str, Any]:
    """
    Deterministic audit of asymptotic z=1.644854 Type-I inflation vs Student-t,
    unblocked cross-donor false >25% shifts vs Donor-1 paired blocking, and
    6+3 vs 6+6 unpooled Welch 90% TOST on [0.75, 1.3333] (IIP-0024; DEC-103).
    """
    z95 = 1.6448536269514722
    type1_z_m3_pct = 100.0 * (1.0 - student_t_cdf(z95, 2.0))
    type1_z_m6_pct = 100.0 * (1.0 - student_t_cdf(z95, 5.0))
    type1_t_m6_pct = 100.0 * (1.0 - student_t_cdf(student_t_ppf_one_sided_95(5.0), 5.0))
    # False >25% (|Delta_log| > ln(1.25)) shift probability under unblocked cross-donor vs blocked Donor-1 contrast (m=6, sigma_batch=0.12):
    margin_25 = math.log(1.25)
    sigma_batch = 0.12
    sd_blocked = math.sqrt(2.0 * (sigma_batch ** 2) / 6.0)
    blocked_false_shift_pct = 100.0 * 2.0 * (1.0 - _norm_cdf_iip0024(margin_25 / sd_blocked))
    unblocked_false_shifts_pct: dict[str, float] = {}
    for s_don in (0.10, 0.20, 0.30):
        sd_unb = math.sqrt(2.0 * (s_don ** 2) + 2.0 * (sigma_batch ** 2) / 6.0)
        unblocked_false_shifts_pct[f'sigma_donor_{s_don:.2f}'] = round(100.0 * 2.0 * (1.0 - _norm_cdf_iip0024(margin_25 / sd_unb)), 2)

    theta_eq = math.log(4.0 / 3.0)
    w_63_0 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 3, 0.12, 0.12, 0.0, theta_eq, n_pts=n_pts)
    w_63_10 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 3, 0.12, 0.12, math.log(1.10), theta_eq, n_pts=n_pts)
    w_63_15 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 3, 0.12, 0.12, math.log(1.15), theta_eq, n_pts=n_pts)
    w_66_0 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 6, 0.12, 0.12, 0.0, theta_eq, n_pts=n_pts)
    w_66_10 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 6, 0.12, 0.12, math.log(1.10), theta_eq, n_pts=n_pts)
    w_66_15 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 6, 0.12, 0.12, math.log(1.15), theta_eq, n_pts=n_pts)
    w_66_type1_sd_006_015 = 100.0 * welch_tost_power_exact_2d_iip0024(6, 6, 0.06, 0.15, theta_eq, theta_eq, n_pts=n_pts)
    return {
        'proposal_id': 'IIP-0024',
        'decision_id': 'DEC-103',
        'asymptotic_z_type1_error_pct': {
            'm_eq_3_df2': round(type1_z_m3_pct, 3),
            'm_eq_6_df5': round(type1_z_m6_pct, 3),
            'student_t_m_eq_6_df5': round(type1_t_m6_pct, 3),
        },
        'construct_contrast_false_gt25pct_shift_pct': {
            'donor1_paired_blocked_m6': round(blocked_false_shift_pct, 2),
            'unblocked_cross_donor_m6': unblocked_false_shifts_pct,
        },
        'welch_90pct_tost_power_on_0p75_1p3333_pct': {
            'lots_6_plus_3': {
                'drift_0pct': round(w_63_0, 2),
                'drift_plus10pct': round(w_63_10, 2),
                'drift_plus15pct': round(w_63_15, 3),
            },
            'lots_6_plus_6': {
                'drift_0pct': round(w_66_0, 2),
                'drift_plus10pct': round(w_66_10, 2),
                'drift_plus15pct': round(w_66_15, 2),
                'boundary_type1_sd_0p06_vs_0p15_pct': round(w_66_type1_sd_006_015, 3),
            },
        },
    }


def evaluate_iip0024_atomic_stage0_ledger(
    batches: list[dict[str, Any]],
    target_provenance: dict[str, Any],
    *,
    target_campaign_id: str,
    root: Path | str | None = None,
    specimen_lifecycle_stage: str = 'FULL_LOT_RELEASE_WITH_SISTER_UNIT_TERMINAL',
    authorized_sister_assembly_ids: frozenset[str] = frozenset(),
) -> dict[str, Any]:
    """
    Fail-closed atomic Stage-0 ledger evaluator enforcing 4/5/7-key provenance scopes,
    adverse observation latching across batch merges until an authorized same-scope QP reset,
    distinct specimen units on Gate 0S dual shedding, and non-executable product holds (IIP-0024; DEC-103).
    """
    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    p_mod = _load_product_contract_module(root_path)
    try:
        mcb_k = MCBScopeKey(
            target_provenance['product_variant_id'],
            target_provenance['construct_sha256'],
            target_provenance['donor_seedline_id'],
            target_provenance['mcb_id'],
        )
        lot_k = LotScopeKey(
            mcb_k.product_variant_id,
            mcb_k.construct_sha256,
            mcb_k.donor_seedline_id,
            mcb_k.mcb_id,
            target_provenance['lot_id'],
        )
    except (KeyError, TypeError, ValueError) as exc:
        return {'readiness_tier': 'HOLD_MISSING_OR_INVALID_IMMUTABLE_PROVENANCE_KEY', 'overall_status': 'HOLD_MISSING_OR_INVALID_IMMUTABLE_PROVENANCE_KEY', 'hard_stop_precedence_winner': None, 'execution_supported': False, 'holds': [str(exc)]}

    try:
        pdef = p_mod.get_product_definition(mcb_k.product_variant_id, root=root_path)
    except ValueError as exc:
        return {'readiness_tier': 'HOLD_UNKNOWN_PRODUCT_VARIANT', 'overall_status': 'HOLD_UNKNOWN_PRODUCT_VARIANT', 'hard_stop_precedence_winner': None, 'execution_supported': False, 'holds': [str(exc)]}

    pvid_canon = str(pdef['product_variant_id'])
    staged_l2_v3_ids = getattr(p_mod, 'STAGED_LOCUS2_V3_PRODUCT_IDS', frozenset())
    is_non_exec_product = (pdef.get('execution_supported') is False)
    is_z_prod = (pvid_canon == 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2')
    is_staged_research_comparator = (pvid_canon == 'ISLET1_B2_5LOCUS_BASE' or pvid_canon in staged_l2_v3_ids)

    mcb_k_c = MCBScopeKey(pvid_canon, mcb_k.construct_sha256, mcb_k.donor_seedline_id, mcb_k.mcb_id)
    lot_k_c = LotScopeKey(pvid_canon, mcb_k.construct_sha256, mcb_k.donor_seedline_id, mcb_k.mcb_id, lot_k.lot_id)

    early_device_or_route_hold: tuple[str, str] | None = None
    dev_k_c: AssembledDeviceLotScopeKey | None = None
    try:
        dev_k = AssembledDeviceLotScopeKey(
            mcb_k.product_variant_id,
            mcb_k.construct_sha256,
            mcb_k.donor_seedline_id,
            mcb_k.mcb_id,
            lot_k.lot_id,
            target_provenance['device_assembly_id'],
            target_provenance['physical_route_id'],
        )
        dev_k_c = AssembledDeviceLotScopeKey(
            pvid_canon,
            mcb_k.construct_sha256,
            mcb_k.donor_seedline_id,
            mcb_k.mcb_id,
            lot_k.lot_id,
            dev_k.device_assembly_id,
            dev_k.physical_route_id,
        )
        routes = ENDPOINT_REGISTRY_SPEC['physical_routes']
        if dev_k.physical_route_id not in routes:
            early_device_or_route_hold = ('HOLD_UNREGISTERED_PHYSICAL_ROUTE_ID', dev_k.physical_route_id)
        else:
            route_spec = routes[dev_k.physical_route_id]
            if is_z_prod and not route_spec['requires_design_z_two_stage']:
                early_device_or_route_hold = ('HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH', dev_k.physical_route_id)
            elif not is_z_prod and route_spec['requires_design_z_two_stage']:
                early_device_or_route_hold = ('HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH', dev_k.physical_route_id)
    except (KeyError, TypeError, ValueError) as exc:
        early_device_or_route_hold = ('HOLD_MISSING_OR_INVALID_IMMUTABLE_PROVENANCE_KEY', str(exc))

    def _canon_scope_tuple_iip0031(sk: Any) -> tuple[str, ...]:
        raw_t = sk.as_tuple()
        c_pid, _ = _canonicalize_product_id_iip0031(raw_t[0], root_path)
        return (c_pid or raw_t[0], *raw_t[1:])

    mcb_store: dict[str, AtomicObservation] = {}
    lot_store: dict[str, AtomicObservation] = {}
    dev_store: dict[str, AtomicObservation] = {}
    atom_status_by_store_key: dict[tuple[str, str], str] = {}
    invalidated_stale_keys: set[str] = set()
    latched_adverse: dict[tuple[str, tuple[str, ...]], dict[str, Any]] = {}
    retained_reset_audit_log: list[dict[str, Any]] = []
    structural_holds: list[str] = []
    malformed_batch_holds: list[str] = []
    tier6b_atoms_collected: list[Any] = []
    s11_panel_eval_collected: dict[str, Any] | None = None
    latched_s11_fail_eval: dict[str, Any] | None = None
    latched_s11_fail_b_idx: int | None = None
    latched_s11_all_fail_reasons: list[str] = []
    elevated_store_gate0s_escalated_on_lot = False
    iip0029_lot_source_fingerprint: tuple[Any, ...] | None = None
    iip0029_lot_store_sample_unit_id: str | None = None

    count_caps = ENDPOINT_REGISTRY_SPEC['count_endpoint_caps_and_zero_vetoes']
    sister_terminal_ids = set(ENDPOINT_REGISTRY_SPEC['device_sister_terminal_endpoint_ids'])
    target_scope_tuples = {mcb_k_c.as_tuple(), lot_k_c.as_tuple()}
    if dev_k_c is not None and early_device_or_route_hold is None:
        target_scope_tuples.add(dev_k_c.as_tuple())

    for b_idx, batch in enumerate(batches):
        if not isinstance(batch, dict):
            malformed_batch_holds.append('HOLD_MALFORMED_BATCH_ATOM_LIST')
            continue
        has_t6b_container = False
        for _nested_name in ('scope', 'measurements', 'iip0026_molecular_observations'):
            if _nested_name in batch and batch[_nested_name] is not None and not isinstance(batch[_nested_name], dict):
                malformed_batch_holds.append('HOLD_MALFORMED_BATCH_ATOM_LIST')
                tier6b_atoms_collected.append(None)
        _b_t6_ctrs = (
            batch,
            batch.get('scope') if isinstance(batch.get('scope'), dict) else {},
            batch.get('iip0026_molecular_observations') if isinstance(batch.get('iip0026_molecular_observations'), dict) else {},
            batch.get('measurements') if isinstance(batch.get('measurements'), dict) else {},
        )
        _b_t6_null_seen = False
        for _b_src in _b_t6_ctrs:
            if 'tier6b_sorted_lot_observations' in _b_src:
                if _b_src['tier6b_sorted_lot_observations'] is None:
                    _b_t6_null_seen = True
                else:
                    has_t6b_container = True
                    _b_t6 = _b_src['tier6b_sorted_lot_observations']
                    if isinstance(_b_t6, (list, tuple)):
                        tier6b_atoms_collected.extend(_b_t6)
                    else:
                        tier6b_atoms_collected.append(None)
        if _b_t6_null_seen and has_t6b_container:
            tier6b_atoms_collected.append(None)
        atoms = batch.get('atoms', [] if (has_t6b_container or 'tier6b_sorted_lot_observations' in batch) else None)
        if not isinstance(atoms, list):
            malformed_batch_holds.append('HOLD_MALFORMED_BATCH_ATOM_LIST')
            continue
        for atom in atoms:
            if not isinstance(atom, AtomicObservation):
                malformed_batch_holds.append('HOLD_NON_ATOMIC_OBSERVATION_OBJECT')
                if isinstance(atom, dict) and (
                    atom.get('endpoint_id') == 'lot_tier6b_sorted_prolif_clonogenic_escape'
                    or (isinstance(atom.get('controls'), dict) and atom['controls'].get('endpoint_id') == 'lot_tier6b_sorted_prolif_clonogenic_escape')
                ):
                    tier6b_atoms_collected.append(atom)
                else:
                    tier6b_atoms_collected.append(None)
                continue
            sk_tup = _canon_scope_tuple_iip0031(atom.scope_key)
            # Route sorted-lot Tier 6B observations exclusively through the IIP-0031 per-fraction thinning evaluator:
            nested_tier6b = atom.controls.get('endpoint_id') == 'lot_tier6b_sorted_prolif_clonogenic_escape'
            if nested_tier6b or atom.endpoint_id == 'lot_tier6b_sorted_prolif_clonogenic_escape':
                tier6b_atoms_collected.append(atom)
            if atom.endpoint_id == 'lot_tier6b_sorted_prolif_clonogenic_escape':
                if atom.scope_type == 'LOT_SCOPE' and sk_tup == lot_k_c.as_tuple():
                    if atom.campaign_id != target_campaign_id:
                        structural_holds.append(f'HOLD_UNLINKED_CAMPAIGN:{atom.endpoint_id}')
                    else:
                        lot_store[atom.endpoint_id] = atom
                elif atom.scope_type == 'LOT_SCOPE' and sk_tup[:4] == mcb_k_c.as_tuple():
                    invalidated_stale_keys.add(atom.endpoint_id)
                continue
            # Evaluate atom validity and latch any adverse failure on its exact scope:
            expected_kind = _ATOMIC_ENDPOINT_KINDS.get(atom.endpoint_id)
            if expected_kind is None:
                atom_status = 'HOLD_UNKNOWN_ATOMIC_ENDPOINT'
            elif atom.kind != expected_kind:
                atom_status = 'HOLD_INVALID_ATOM_KIND_OR_NON_NUMERIC_VALUE'
            elif atom.scope_type != _ATOMIC_REQUIRED_SCOPES.get(atom.endpoint_id, atom.scope_type):
                atom_status = 'HOLD_INVALID_ATOMIC_ENDPOINT_SCOPE'
            elif atom.endpoint_id in S9_CORE_LOT_ENDPOINT_IDS and not (
                _S9_ASSAY_RANGES[atom.endpoint_id][0] <= atom.value_or_k <= _S9_ASSAY_RANGES[atom.endpoint_id][1]
                and STAGE0_CUSTOM_MINIMUM_SE.get(atom.endpoint_id, 0.001) <= atom.se_or_n <= 1e12
            ):
                atom_status = 'HOLD_S9_OUT_OF_RANGE_MEASUREMENT_OR_REPLICATE_COUNT'
            elif atom.kind == 'CONTINUOUS_STUDENT_T' and atom.endpoint_id not in S9_CORE_LOT_ENDPOINT_IDS and atom.sample_n > 10000:
                atom_status = 'HOLD_S9_OUT_OF_RANGE_MEASUREMENT_OR_REPLICATE_COUNT'
            elif atom.kind == 'METADATA_CONTROL':
                atom_status = 'PASS' if atom.value_or_k is True else 'HOLD_UNVERIFIED_METADATA_CONTROL'
            elif atom.kind == 'EXACT_INTEGER_COUNT':
                spec_c = count_caps[atom.endpoint_id]
                ctrl = dict(atom.controls)
                c_res = evaluate_rare_event_count_iip0024(
                    gate_id=('GATE0S_SHEDDING_' + atom.endpoint_id.upper() if atom.endpoint_id.startswith('gate0s_shedding_') else atom.endpoint_id.upper()),
                    k=atom.value_or_k,
                    n=atom.se_or_n,
                    max_u95_fraction=float(spec_c['cap']),
                    require_observed_zero=bool(spec_c['require_observed_zero']),
                    sampling_independence_model=str(ctrl.get('sampling_independence_model', 'INDEPENDENT_BIOLOGICAL_OPPORTUNITIES')),
                    clustered_sampling=bool(ctrl.get('clustered_sampling', False)),
                    validated_clustered_estimator_id=ctrl.get('validated_clustered_estimator_id'),
                    eta_spike_recovery_l95=ctrl.get('eta_spike_recovery_l95'),
                    f_clonogenic=ctrl.get('f_clonogenic'),
                    use_poisson_protocol_floor=bool(spec_c.get('use_poisson_protocol_floor', False)),
                )
                atom_status = c_res['status']
            elif atom.kind == 'CONTINUOUS_STUDENT_T':
                # The immutable constructor validated the numeric mean, SE and sample_n >= 6.
                atom_status = 'HOLD_UNEVALUATED_CONTINUOUS_ENDPOINT'
                _eff_cont_n = min(atom.sample_n, 10000)
                t_95_a = student_t_ppf_one_sided_95(float(_eff_cont_n - 1))
                u95_a = float(atom.value_or_k) + t_95_a * float(atom.se_or_n)
                l95_a = float(atom.value_or_k) - t_95_a * float(atom.se_or_n)
                if atom.endpoint_id in ('mcb_1q32_mdm4_diploid_cn', 'mcb_20q11_bcl2l1_diploid_cn'):
                    atom_status = 'PASS' if 1.85 <= l95_a <= u95_a <= 2.15 else 'FAIL_HARD_STOP'
                elif atom.endpoint_id == 'mcb_donor_pairwise_ibd_pi_hat':
                    atom_status = 'PASS' if 0 <= atom.value_or_k and u95_a < 0.05 else 'HOLD_DONOR_PAIRWISE_IBD'
                elif atom.endpoint_id == 'obs_iapp_d84_immunocapture_oligomer_elisa_fold':
                    atom_status = 'PASS' if 0 <= atom.value_or_k and u95_a <= 1.5 else 'HOLD_INCONCLUSIVE_IAPP_EVIDENCE'
                elif atom.endpoint_id == 'obs_iapp_d84_colocalized_fibril_area_pct':
                    atom_status = 'PASS' if 0 <= atom.value_or_k and u95_a <= 1.0 else 'HOLD_INCONCLUSIVE_IAPP_EVIDENCE'
                elif atom.endpoint_id == 'stage7_stored_mature_insulin_iu':
                    ctrl_a = dict(atom.controls)
                    max_lot_iu_a = ctrl_a.get('max_lot_iu', float(atom.value_or_k))
                    if not _is_strict_finite_number_iip0024(max_lot_iu_a) or float(max_lot_iu_a) < 0.0:
                        atom_status = 'HOLD_INVALID_STORED_MATURE_INSULIN_MAX_LOT_IU'
                    else:
                        crossed_138p3 = (
                            u95_a > float(STAGE0_IIP0029_CONSTANTS['GATE0S_ABSOLUTE_STORED_INSULIN_CEILING_IU'])
                            or float(max_lot_iu_a) > float(STAGE0_IIP0029_CONSTANTS['GATE0S_ABSOLUTE_STORED_INSULIN_CEILING_IU'])
                        )
                        if atom.scope_type == 'LOT_SCOPE' and sk_tup == lot_k_c.as_tuple():
                            if elevated_store_gate0s_escalated_on_lot and not crossed_138p3:
                                atom_status = 'FAIL_HARD_STOP'
                                structural_holds.append('FAIL_IIP0029_IMMUTABLE_SOURCE_OR_THRESHOLD_REPLACEMENT_DETECTED')
                            elif iip0029_lot_store_sample_unit_id is not None and atom.sample_unit_id != iip0029_lot_store_sample_unit_id and not crossed_138p3 and elevated_store_gate0s_escalated_on_lot:
                                atom_status = 'FAIL_HARD_STOP'
                                structural_holds.append('FAIL_IIP0029_IMMUTABLE_SOURCE_OR_THRESHOLD_REPLACEMENT_DETECTED')
                            else:
                                iip0029_lot_store_sample_unit_id = atom.sample_unit_id
                                if crossed_138p3:
                                    elevated_store_gate0s_escalated_on_lot = True
                                atom_status = 'PASS'
                elif atom.endpoint_id == 'gate0s_ap1903_1h_insulin_burst_iu':
                    atom_status = 'PASS' if u95_a <= 2.8 else 'FAIL_HARD_STOP'
                elif atom.endpoint_id == 'gate0s_ap1903_24h_insulin_dump_fraction':
                    atom_status = 'PASS' if u95_a <= 0.15 else 'FAIL_HARD_STOP'
                elif atom.endpoint_id == 'obs_explant_en_bloc_recovery_fraction':
                    atom_status = 'PASS' if l95_a >= 0.99 else 'FAIL_HARD_STOP'
                elif atom.endpoint_id == 'obs_icasp9_log10_kill':
                    ctrl_ic9 = dict(atom.controls)
                    _b_n_ctrs = (
                        ctrl_ic9,
                        batch,
                        batch.get('scope') if isinstance(batch.get('scope'), dict) else {},
                        batch.get('measurements') if isinstance(batch.get('measurements'), dict) else {},
                        batch.get('iip0026_molecular_observations') if isinstance(batch.get('iip0026_molecular_observations'), dict) else {},
                        target_provenance,
                    )
                    _ic9_n_vals = [atom.sample_n] + [
                        _c_n[k]
                        for _c_n in _b_n_ctrs
                        for k in (*S9_REPLICATE_KEYS['default'], *S9_REPLICATE_KEYS['obs_icasp9_log10_kill'])
                        if k in _c_n
                    ]
                    _ic9_n_bad = any(
                        not _is_integral_count_iip0026(v) or not (6 <= int(v) <= 10000) or int(v) != atom.sample_n
                        for v in _ic9_n_vals
                    )
                    _ic9_valid_ge2 = [
                        (v if type(v) is int else int(round(v)))
                        for v in _ic9_n_vals
                        if (
                            (type(v) is int and 2 <= v <= 10000)
                            or (isinstance(v, float) and math.isfinite(v) and 2.0 <= v <= 10000.0 and v.is_integer())
                        )
                    ]
                    _ic9_eff_n = min(_ic9_valid_ge2 + ([6] if _ic9_n_bad else [])) if _ic9_valid_ge2 else 6
                    _ic9_t95 = student_t_ppf_one_sided_95(float(_ic9_eff_n - 1))
                    _ic9_l95 = float(atom.value_or_k) - _ic9_t95 * float(atom.se_or_n)
                    if _ic9_l95 < 4.0:
                        atom_status = 'FAIL_HARD_STOP'
                    elif _ic9_n_bad:
                        atom_status = 'HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT'
                    else:
                        atom_status = 'PASS'
                elif atom.endpoint_id in ('obs_tp53_ppm1d_bcor_vaf', 'obs_castseq_translocation_freq'):
                    ctrl_s9 = dict(atom.controls)
                    _b_n_ctrs = (
                        ctrl_s9,
                        batch,
                        batch.get('scope') if isinstance(batch.get('scope'), dict) else {},
                        batch.get('measurements') if isinstance(batch.get('measurements'), dict) else {},
                        batch.get('iip0026_molecular_observations') if isinstance(batch.get('iip0026_molecular_observations'), dict) else {},
                        target_provenance,
                    )
                    _s9_n_vals = [atom.sample_n] + [
                        _c_n[k]
                        for _c_n in _b_n_ctrs
                        for k in (*S9_REPLICATE_KEYS['default'], *S9_REPLICATE_KEYS[atom.endpoint_id])
                        if k in _c_n
                    ]
                    _s9_n_bad = any(
                        not _is_integral_count_iip0026(v) or not (6 <= int(v) <= 10000) or int(v) != atom.sample_n
                        for v in _s9_n_vals
                    )
                    _s9_valid_ge2 = [
                        (v if type(v) is int else int(round(v)))
                        for v in _s9_n_vals
                        if (
                            (type(v) is int and 2 <= v <= 10000)
                            or (isinstance(v, float) and math.isfinite(v) and 2.0 <= v <= 10000.0 and v.is_integer())
                        )
                    ]
                    _s9_eff_n = min(_s9_valid_ge2 + ([6] if _s9_n_bad else [])) if _s9_valid_ge2 else 6
                    depth_key, count_key = S9_FREQUENCY_ASSAY_KEYS[atom.endpoint_id]
                    depth_raw = ctrl_s9.get(depth_key, ctrl_s9.get('duplex_depth'))
                    s9_res = _evaluate_s9_tp53_vaf_guard(
                        atom.value_or_k,
                        atom.se_or_n,
                        depth_raw,
                        sample_n=_s9_eff_n,
                        mutant_molecules_k=ctrl_s9.get(count_key),
                        require_explicit_depth=True,
                        endpoint_label=atom.endpoint_id,
                    )
                    atom_status = str(s9_res['status'])
                    if not s9_res['hard_stop'] and _s9_n_bad:
                        atom_status = 'HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT'
                    elif not s9_res['hard_stop'] and 'duplex_depth' in ctrl_s9 and depth_key in ctrl_s9 and (
                        not _is_integral_count_iip0026(ctrl_s9['duplex_depth'])
                        or ctrl_s9['duplex_depth'] != ctrl_s9[depth_key]
                    ):
                        atom_status = 'HOLD_S9_CONFLICTING_DUPLEX_DEPTH'
            elif atom.kind == 'STRUCTURED_PANEL':
                p_dict = _deep_to_jsonable_iip0024(atom.value_or_k)
                if 'imprinting_' in atom.endpoint_id or 'meg3_dmr' in atom.endpoint_id:
                    atom_status = _atomic_mcb_panel_status(evaluate_conjunctive_imprinting_dmr_iip0024, atom.scope_type, p_dict)
                elif 'xci_46xx' in atom.endpoint_id:
                    atom_status = _atomic_mcb_panel_status(evaluate_conjunctive_xci_46xx_iip0024, atom.scope_type, p_dict)
                elif 'mtdna_1000x' in atom.endpoint_id:
                    atom_status = _atomic_mcb_panel_status(evaluate_scoped_mtdna_gate_iip0024, atom.scope_type, p_dict)
                elif atom.endpoint_id == 'stage0_locus2_modeb_v3_panel':
                    if atom.scope_type != 'LOT_SCOPE':
                        atom_status = 'HOLD_S11_PANEL_REQUIRES_LOT_SCOPE'
                        if sk_tup in target_scope_tuples:
                            structural_holds.append('HOLD_S11_PANEL_REQUIRES_LOT_SCOPE')
                    elif not isinstance(p_dict, dict):
                        atom_status = 'HOLD_S11_MALFORMED_RECORD'
                        if sk_tup == lot_k_c.as_tuple():
                            structural_holds.append('HOLD_S11_MALFORMED_RECORD')
                    else:
                        outer_ctx_31 = {
                            'product_variant_id': pvid_canon,
                            'construct_sha256': atom.scope_key.construct_sha256,
                            'donor_seedline_id': atom.scope_key.donor_seedline_id,
                            'mcb_id': atom.scope_key.mcb_id,
                            'lot_id': atom.scope_key.lot_id,
                            'campaign_id': atom.campaign_id,
                            'sample_id': atom.sample_unit_id,
                            'scope_type': atom.scope_type,
                            **({'background_id': target_provenance['background_id']} if 'background_id' in target_provenance else {}),
                            **({'l2_v3_donor_id': target_provenance['l2_v3_donor_id']} if 'l2_v3_donor_id' in target_provenance else {}),
                            **({'donor_id': target_provenance['donor_id']} if 'donor_id' in target_provenance else {}),
                        }
                        res_31 = evaluate_locus2_ciita_modeb_challenger_v3_gate(
                            {**outer_ctx_31, 'iip0031_locus2_record': dict(p_dict)},
                            root=root_path,
                            product_def=pdef,
                            _trusted_p_mod=p_mod,
                        )
                        st_31 = str(res_31.get('status'))
                        if sk_tup == lot_k_c.as_tuple() and atom.campaign_id == target_campaign_id:
                            if st_31.startswith('FAIL_'):
                                for _rf31 in res_31.get('reasons', [st_31]):
                                    if str(_rf31).startswith('FAIL_') and str(_rf31) not in latched_s11_all_fail_reasons:
                                        latched_s11_all_fail_reasons.append(str(_rf31))
                                if latched_s11_fail_eval is None or (
                                    b_idx == latched_s11_fail_b_idx
                                    and _pick_headline_s11_fail_iip0031([str(res_31['status']), str(latched_s11_fail_eval['status'])]) == str(res_31['status'])
                                ):
                                    latched_s11_fail_eval = dict(res_31)
                                    latched_s11_fail_b_idx = b_idx
                            s11_panel_eval_collected = res_31
                        atom_status = 'PASS' if res_31.get('gate_pass') is True else st_31
                elif atom.endpoint_id == 'stage0_twoclock_repeated_challenge_panel':
                    if atom.scope_type != 'LOT_SCOPE':
                        atom_status = 'HOLD_IIP0029_PANEL_REQUIRES_LOT_SCOPE'
                    else:
                        rec_29 = dict(p_dict) if isinstance(p_dict, dict) else {}
                        rec_29.setdefault(
                            'lot_scope_key',
                            {
                                'product_variant_id': atom.scope_key.product_variant_id,
                                'construct_sha256': atom.scope_key.construct_sha256,
                                'donor_seedline_id': atom.scope_key.donor_seedline_id,
                                'mcb_id': atom.scope_key.mcb_id,
                                'lot_id': atom.scope_key.lot_id,
                            },
                        )
                        res_29 = evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(rec_29, root=root_path)
                        fp_29 = _extract_iip0029_source_fingerprint(rec_29)
                        if sk_tup == lot_k_c.as_tuple():
                            if iip0029_lot_source_fingerprint is not None and fp_29 != iip0029_lot_source_fingerprint:
                                atom_status = 'FAIL_HARD_STOP'
                                structural_holds.append('FAIL_IIP0029_IMMUTABLE_SOURCE_OR_THRESHOLD_REPLACEMENT_DETECTED')
                            else:
                                iip0029_lot_source_fingerprint = fp_29
                                if res_29.get('gate0s_escalation_required') is True:
                                    elevated_store_gate0s_escalated_on_lot = True
                                if res_29.get('subgate_pass') is True:
                                    atom_status = 'PASS'
                                elif res_29.get('gate0s_escalation_required') is True and res_29.get('all_non_store_ceiling_checks_passed') is True:
                                    atom_status = 'PASS_PENDING_GATE0S_ESCALATION'
                                else:
                                    atom_status = str(res_29.get('status', 'HOLD_IIP0029_SUBGATE_FAILED'))
                        else:
                            atom_status = 'PASS' if res_29.get('subgate_pass') is True else str(res_29.get('status', 'HOLD_IIP0029_SUBGATE_FAILED'))
                else:
                    atom_status = 'HOLD_UNKNOWN_STRUCTURED_PANEL'

            latch_key = (atom.endpoint_id, sk_tup)
            _campaign_linked_for_latch = (
                atom.scope_type == 'MCB_SCOPE'
                or atom.campaign_id == target_campaign_id
            )
            if atom_status == 'FAIL_HARD_STOP':
                if _campaign_linked_for_latch:
                    latched_adverse[latch_key] = {'stop_idx': b_idx, 'repeats': {}}
                elif sk_tup in target_scope_tuples:
                    structural_holds.append(f'HOLD_UNLINKED_CAMPAIGN_ADVERSE_OBSERVATION:{atom.endpoint_id}:{atom.campaign_id}')
            elif (
                latch_key in latched_adverse
                and atom_status == 'PASS'
                and b_idx > latched_adverse[latch_key]['stop_idx']
                and batch.get('reset_prior_safety_stop_after_investigation') is not True
            ):
                raw_bid = batch.get('batch_id')
                if (
                    dict(atom.controls).get('washout_and_retest_verified') is True
                    and atom.campaign_id == target_campaign_id
                    and isinstance(raw_bid, str)
                    and raw_bid.strip()
                    and raw_bid == raw_bid.strip()
                ):
                    latched_adverse[latch_key]['repeats'][raw_bid] = atom.sample_unit_id

            if atom.scope_type == 'MCB_SCOPE':
                if sk_tup == mcb_k_c.as_tuple():
                    mcb_store[atom.endpoint_id] = atom
                    atom_status_by_store_key[('MCB_SCOPE', atom.endpoint_id)] = atom_status
            elif atom.scope_type == 'LOT_SCOPE':
                if sk_tup == lot_k_c.as_tuple():
                    if atom.campaign_id != target_campaign_id:
                        structural_holds.append(f'HOLD_UNLINKED_CAMPAIGN:{atom.endpoint_id}')
                    else:
                        lot_store[atom.endpoint_id] = atom
                        atom_status_by_store_key[('LOT_SCOPE', atom.endpoint_id)] = atom_status
                elif sk_tup[:4] == mcb_k_c.as_tuple():
                    invalidated_stale_keys.add(atom.endpoint_id)
            elif atom.scope_type == 'ASSEMBLED_DEVICE_LOT_SCOPE':
                same_or_sister = (
                    dev_k_c is not None
                    and early_device_or_route_hold is None
                    and sk_tup[:5] == lot_k_c.as_tuple()
                    and atom.scope_key.physical_route_id == dev_k_c.physical_route_id
                    and (
                        atom.scope_key.device_assembly_id == dev_k_c.device_assembly_id
                        or atom.scope_key.device_assembly_id in authorized_sister_assembly_ids
                    )
                )
                if same_or_sister:
                    if atom.scope_key.device_assembly_id != dev_k_c.device_assembly_id and atom.endpoint_id not in sister_terminal_ids:
                        structural_holds.append(f'HOLD_SISTER_UNIT_NON_TERMINAL_TRANSFER_FORBIDDEN:{atom.endpoint_id}')
                    elif atom.campaign_id != target_campaign_id:
                        structural_holds.append(f'HOLD_UNLINKED_CAMPAIGN:{atom.endpoint_id}')
                    else:
                        dev_store[atom.endpoint_id] = atom
                        atom_status_by_store_key[('ASSEMBLED_DEVICE_LOT_SCOPE', atom.endpoint_id)] = atom_status
                elif sk_tup[:4] == mcb_k_c.as_tuple():
                    invalidated_stale_keys.add(atom.endpoint_id)

        if batch.get('reset_prior_safety_stop_after_investigation') is True:
            req_eid = batch.get('reset_prior_safety_stop_target_category')
            rep_ids = batch.get('reset_prior_safety_stop_repeat_batch_ids')
            qp_id = batch.get('reset_prior_safety_stop_independent_qp_approver_id')
            sub_id = batch.get('submitter_id')
            r_code = batch.get('reset_prior_safety_stop_reason_code')
            cur_bid = batch.get('batch_id')
            if 'reset_prior_safety_stop_target_scope_tuple' in batch and batch['reset_prior_safety_stop_target_scope_tuple'] is not None:
                req_scope = batch['reset_prior_safety_stop_target_scope_tuple']
                if isinstance(req_scope, (list, tuple)) and all(isinstance(x, str) and x.strip() and x == x.strip() for x in req_scope):
                    allowed_scopes_for_reset = {tuple(req_scope)}
                else:
                    allowed_scopes_for_reset = set()
            else:
                allowed_scopes_for_reset = target_scope_tuples
            if (
                isinstance(req_eid, str)
                and req_eid.strip()
                and req_eid == req_eid.strip()
                and isinstance(qp_id, str)
                and qp_id.strip()
                and qp_id == qp_id.strip()
                and isinstance(sub_id, str)
                and sub_id.strip()
                and sub_id == sub_id.strip()
                and qp_id != sub_id
                and isinstance(r_code, str)
                and r_code in ALLOWED_STICKY_S9_RESET_REASON_CODES
                and isinstance(rep_ids, list)
                and all(isinstance(x, str) and x.strip() and x == x.strip() for x in rep_ids)
                and len(set(rep_ids)) >= 2
                and (not isinstance(cur_bid, str) or cur_bid not in rep_ids)
                and bool(allowed_scopes_for_reset)
            ):
                for lk_cand in list(latched_adverse.keys()):
                    if lk_cand[0] == req_eid and lk_cand[1] in allowed_scopes_for_reset and lk_cand[1] in target_scope_tuples:
                        reps_map = latched_adverse[lk_cand]['repeats']
                        if set(rep_ids).issubset(set(reps_map.keys())) and len({reps_map[r] for r in rep_ids if str(reps_map[r]).strip()}) >= 2:
                            retained_reset_audit_log.append({
                                'cleared_endpoint_id': req_eid,
                                'scope_tuple': list(lk_cand[1]),
                                'qp_approver_id': qp_id,
                                'submitter_id': sub_id,
                                'repeat_batch_ids': list(rep_ids),
                                'reason_code': r_code,
                            })
                            del latched_adverse[lk_cand]

    tier6b_eval = evaluate_tier6b_per_fraction_thinning_iip0031(
        tier6b_atoms_collected,
        target_scope={**target_provenance, 'product_variant_id': pvid_canon, 'campaign_id': target_campaign_id},
        root=root_path,
    )
    if latched_s11_fail_eval is not None:
        s11_eval = dict(latched_s11_fail_eval)
        merged_s11_reasons = list(latched_s11_all_fail_reasons) + [
            r for r in s11_eval.get('reasons', []) if r not in latched_s11_all_fail_reasons
        ]
        s11_eval['reasons'] = merged_s11_reasons
    elif s11_panel_eval_collected is not None:
        s11_eval = s11_panel_eval_collected
    else:
        s11_eval = evaluate_locus2_ciita_modeb_challenger_v3_gate(
            {
                'product_variant_id': pvid_canon,
                'construct_sha256': mcb_k_c.construct_sha256,
                'donor_seedline_id': lot_k_c.donor_seedline_id,
                'mcb_id': lot_k_c.mcb_id,
                'lot_id': lot_k_c.lot_id,
                **({'background_id': target_provenance['background_id']} if 'background_id' in target_provenance else {}),
                **({'l2_v3_donor_id': target_provenance['l2_v3_donor_id']} if 'l2_v3_donor_id' in target_provenance else {}),
                **({'donor_id': target_provenance['donor_id']} if 'donor_id' in target_provenance else {}),
            },
            root=root_path,
            product_def=pdef,
            _trusted_p_mod=p_mod,
        )
    _t6b_meta = {
        'tier6b_sorted_lot_evaluation': tier6b_eval,
        'tier6b_disposition': tier6b_eval['disposition'],
        'tier6b_reason': tier6b_eval['reason_code'],
        's11_evaluation': s11_eval,
        's9_point_estimate_fallback_se_policy': 'SCALE_AWARE_DEFAULT_SE_IMPUTED_ON_SCALAR_OMISSION_TYPED_ATOM_REQUIRES_MEASURED_SE',
        's9_se_imputed_endpoints': [],
    }

    holds: list[str] = list(structural_holds)
    if early_device_or_route_hold is not None:
        holds.append(f'{early_device_or_route_hold[0]}:{early_device_or_route_hold[1]}')
    if elevated_store_gate0s_escalated_on_lot:
        _g0s_req_ids = (
            'gate0s_ap1903_1h_insulin_burst_iu',
            'gate0s_ap1903_24h_insulin_dump_fraction',
            'obs_explant_en_bloc_recovery_fraction',
        )
        _g0s_all_ok = all(
            (eid in dev_store) and (atom_status_by_store_key.get(('ASSEMBLED_DEVICE_LOT_SCOPE', eid)) == 'PASS')
            for eid in _g0s_req_ids
        )
        if _g0s_all_ok:
            if atom_status_by_store_key.get(('LOT_SCOPE', 'stage0_twoclock_repeated_challenge_panel')) == 'PASS_PENDING_GATE0S_ESCALATION':
                atom_status_by_store_key[('LOT_SCOPE', 'stage0_twoclock_repeated_challenge_panel')] = 'PASS'
        else:
            holds.append('HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW')
    for (_sc_t, _eid_t), _st_t in sorted(atom_status_by_store_key.items()):
        if _st_t == 'PASS_PENDING_GATE0S_ESCALATION':
            if 'HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW' not in holds:
                holds.append('HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW')
        elif _eid_t == 'stage0_locus2_modeb_v3_panel' and _st_t.startswith(('HOLD_S11_', 'FAIL_S11_')):
            continue
        elif _st_t not in ('PASS', 'FAIL_HARD_STOP'):
            holds.append(f'{_st_t}:{_eid_t}')

    # PRECEDENCE 1: Independent S9 Hard Safety Stop ALWAYS keeps first precedence:
    active_target_latches = [eid for (eid, skt) in latched_adverse if skt in target_scope_tuples]
    if active_target_latches:
        return {
            'readiness_tier': 'HARD_SAFETY_STOP_REJECTED',
            'overall_status': 'HARD_SAFETY_STOP_REJECTED',
            'hard_stop_precedence_winner': 'S9_HARD_SAFETY_STOP',
            'execution_supported': False,
            'holds': [f'latched_adverse:{e}' for e in sorted(active_target_latches)] + holds,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }

    # PRECEDENCE 2: Attributable sorted-lot Tier 6B integer k_i >= 1 latches FAIL_HARD_STOP BEFORE
    # calibration holds, malformed records, completeness holds, or non-executable research holds:
    if tier6b_eval.get('hard_stop') is True:
        return {
            'readiness_tier': 'HARD_SAFETY_STOP_REJECTED',
            'overall_status': 'FAIL_HARD_STOP',
            'hard_stop_precedence_winner': 'FAIL_TIER6B_OUTGROWTH_DETECTED',
            'execution_supported': False,
            'holds': ['FAIL_TIER6B_OUTGROWTH_DETECTED'] + holds,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }

    # PRECEDENCE 3: Rule S11 hard architectural/safety stops evaluate BEFORE non-executable research holds:
    if str(s11_eval.get('status', '')).startswith('FAIL_'):
        return {
            'readiness_tier': str(s11_eval['status']),
            'overall_status': str(s11_eval['status']),
            'hard_stop_precedence_winner': str(s11_eval['status']),
            'execution_supported': False,
            'holds': list(s11_eval.get('reasons', [s11_eval['status']])) + holds,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }

    if early_device_or_route_hold is not None:
        return {
            'readiness_tier': early_device_or_route_hold[0],
            'overall_status': early_device_or_route_hold[0],
            'hard_stop_precedence_winner': None,
            'execution_supported': False,
            'holds': [early_device_or_route_hold[1]] + [h for h in holds if h != f'{early_device_or_route_hold[0]}:{early_device_or_route_hold[1]}'],
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }

    # PRECEDENCE 4: Malformed batch or non-AtomicObservation holds (only when no S9/Tier 6B/S11 hard stop latched):
    if malformed_batch_holds:
        return {
            'readiness_tier': malformed_batch_holds[0],
            'overall_status': malformed_batch_holds[0],
            'hard_stop_precedence_winner': None,
            'execution_supported': False,
            'holds': malformed_batch_holds + holds,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }

    if 'gate0s_shedding_ddpcr_assay' in dev_store and 'gate0s_shedding_viable_outgrowth_assay' in dev_store:
        if dev_store['gate0s_shedding_ddpcr_assay'].sample_unit_id == dev_store['gate0s_shedding_viable_outgrowth_assay'].sample_unit_id:
            holds.append('HOLD_GATE0S_DUAL_SHEDDING_REQUIRES_INDEPENDENT_SPECIMEN_UNITS')

    # PRECEDENCE 5: Mandatory core S9 safety observations, Rule S11 holds, and mandatory sorted-lot Tier 6B per-fraction thinning
    # (matching evaluate_stage0_product_contract hold ordering identically across all 4 entry paths):
    missing_s9_lot = [
        eid for eid in S9_CORE_LOT_ENDPOINT_IDS
        if eid not in lot_store or lot_store[eid].kind != 'CONTINUOUS_STUDENT_T'
        or atom_status_by_store_key.get(('LOT_SCOPE', eid)) != 'PASS'
    ]
    if missing_s9_lot and 'HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS' not in holds:
        holds.append('HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS')
    if s11_eval.get('gate_pass') is not True:
        _s11_h = str(s11_eval['status'])
        if _s11_h not in holds:
            holds.append(_s11_h)
    if tier6b_eval.get('gate_pass') is not True:
        _t6b_h = str(tier6b_eval['status'])
        if _t6b_h not in holds:
            holds.append(_t6b_h)

    missing_mcb = [mk for mk in ENDPOINT_REGISTRY_SPEC['mcb_mandatory_endpoint_ids'] if mk not in mcb_store]
    if holds:
        return {
            'readiness_tier': 'HOLD_FAILED_STAGE0_GATE',
            'overall_status': holds[0],
            'hard_stop_precedence_winner': None,
            'execution_supported': False,
            'holds': holds + missing_mcb + missing_s9_lot,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }
    if missing_mcb or missing_s9_lot or (not lot_store and tier6b_eval.get('gate_pass') is not True) or (specimen_lifecycle_stage != 'AUTHORIZED_TEST_ARTICLE_PRE_IMPLANT' and any(tk not in dev_store for tk in sister_terminal_ids)):
        unresolved_stale_keys = {k for k in invalidated_stale_keys if k not in lot_store and k not in dev_store}
        tier = (
            'HOLD_LOT_OR_DEVICE_SCOPED_EVIDENCE_MISSING_AFTER_PROVENANCE_CHANGE'
            if unresolved_stale_keys
            else 'PROVISIONAL_PARTIAL_LOT_SCREENING_ONLY'
        )
        return {
            'readiness_tier': tier,
            'overall_status': tier,
            'hard_stop_precedence_winner': None,
            'execution_supported': False,
            'holds': missing_mcb or ['incomplete_lot_or_device_panel'],
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }
    if is_staged_research_comparator:
        return {
            'readiness_tier': 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC',
            'overall_status': 'RESEARCH_COMPARATOR_STAGE0_PANEL_CLEARED_PENDING_STATE4_PROMOTION',
            'hard_stop_precedence_winner': None,
            'holds': [f'HOLD_NON_EXECUTABLE_STAGE0_RESEARCH_COMPARATOR:{pvid_canon}'],
            'execution_supported': False,
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }
    if is_non_exec_product and not is_z_prod:
        return {
            'readiness_tier': 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD',
            'overall_status': 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD',
            'hard_stop_precedence_winner': None,
            'holds': [f'non_executable_product:{pvid_canon}'],
            'execution_supported': False,
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }
    if specimen_lifecycle_stage == 'AUTHORIZED_TEST_ARTICLE_PRE_IMPLANT':
        return {
            'readiness_tier': 'AUTHORIZED_TEST_ARTICLE_ASSEMBLED_PENDING_LONGITUDINAL_READOUT',
            'overall_status': 'AUTHORIZED_TEST_ARTICLE_ASSEMBLED_PENDING_LONGITUDINAL_READOUT',
            'hard_stop_precedence_winner': None,
            'holds': [],
            'execution_supported': not is_non_exec_product,
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }
    if is_z_prod:
        return {
            'readiness_tier': 'WILDCARD_COUPON_QUALIFIED_NON_EXECUTABLE',
            'overall_status': 'WILDCARD_COUPON_QUALIFIED_NON_EXECUTABLE',
            'hard_stop_precedence_winner': None,
            'holds': [],
            'execution_supported': False,
            'wildcard_coupon_qualified': True,
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
            **_t6b_meta,
        }
    final_tier = (
        'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'
        if retained_reset_audit_log
        else 'READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE'
    )
    return {
        'readiness_tier': final_tier,
        'overall_status': final_tier,
        'hard_stop_precedence_winner': None,
        'holds': [],
        'execution_supported': True,
        'transferred_mcb_count': len(mcb_store),
        'retained_reset_audit_log': retained_reset_audit_log,
        **_t6b_meta,
    }


def evaluate_coupled_catalase_endotoxin_gate_iip0025(
    *,
    m_total_cat_mg: Any,
    eu_per_mg_mean: Any,
    eu_per_mg_se: Any,
    total_implant_eu_mean: Any,
    total_implant_eu_se: Any,
    replicates_n: Any = 6,
) -> dict[str, Any]:
    """
    Fail-closed unrounded mass-coupled human P04040 catalase endotoxin validator enforcing BOTH
    U95_t(EU/mg) <= min(0.10, 0.50 / m_total_CAT_mg) (rejecting 0.05505 at 9.0830769 mg where
    exact cap is 0.0550474255) AND canonical total implant endotoxin U95_t <= 0.50 EU (IIP-0025; DEC-104).
    """
    if not _is_strict_int_iip0024(replicates_n) or replicates_n < 6:
        return {'gate_pass': False, 'status': 'HOLD_STAGE0_INDEPENDENT_LOT_REPLICATES_BELOW_6_STUDENT_T_MINIMUM'}
    for v_num in (m_total_cat_mg, eu_per_mg_mean, eu_per_mg_se, total_implant_eu_mean, total_implant_eu_se):
        if not _is_strict_finite_number_iip0024(v_num):
            return {'gate_pass': False, 'status': 'HOLD_NONFINITE_CATALASE_ENDOTOXIN_INPUT'}
    m_mg = float(m_total_cat_mg)
    c_mean, c_se = float(eu_per_mg_mean), float(eu_per_mg_se)
    tot_mean, tot_se = float(total_implant_eu_mean), float(total_implant_eu_se)
    if m_mg <= 0.0 or c_mean < 0.0 or c_se <= 0.0 or tot_mean < 0.0 or tot_se <= 0.0:
        return {'gate_pass': False, 'status': 'HOLD_INVALID_CATALASE_MASS_OR_SE_BOUNDS'}
    t95 = student_t_ppf_one_sided_95(float(replicates_n - 1))
    u95_eu_per_mg = c_mean + t95 * c_se
    u95_tot_eu = tot_mean + t95 * tot_se
    exact_cap = min(0.10, 0.50 / m_mg)
    expected_tot = m_mg * c_mean
    if u95_eu_per_mg > exact_cap:
        return {
            'gate_pass': False,
            'status': 'HOLD_P05_MASS_COUPLED_CATALASE_ENDOTOXIN_CAP_EXCEEDED',
            'exact_unrounded_cap_eu_per_mg': exact_cap,
            'u95_eu_per_mg': u95_eu_per_mg,
        }
    if u95_tot_eu > 0.50:
        return {
            'gate_pass': False,
            'status': 'HOLD_P05_TOTAL_IMPLANT_CATALASE_ENDOTOXIN_GT_0P50EU',
            'u95_total_implant_eu': u95_tot_eu,
        }
    if abs(tot_mean - expected_tot) > max(0.015, 0.05 * expected_tot):
        return {
            'gate_pass': False,
            'status': 'HOLD_P05_COUPLED_TOTAL_ENDOTOXIN_MASS_PRODUCT_INCONSISTENT',
            'expected_total_eu': expected_tot,
            'reported_total_eu': tot_mean,
        }
    return {
        'gate_pass': True,
        'status': 'CLEARED_COUPLED_UNROUNDED_CATALASE_ENDOTOXIN_GATE',
        'exact_unrounded_cap_eu_per_mg': exact_cap,
        'u95_eu_per_mg': u95_eu_per_mg,
        'u95_total_implant_eu': u95_tot_eu,
    }


def evaluate_dual_assay_gate0s_shedding_iip0025(
    *,
    ddpcr_k: Any,
    ddpcr_n: Any,
    ddpcr_eta_l95: Any = 1.0,
    ddpcr_specimen_unit_id: str,
    viable_k: Any,
    viable_n: Any,
    viable_eta_l95: Any = 1.0,
    viable_f_clonogenic: Any = 1.0,
    viable_specimen_unit_id: str,
    day90_recovery_mean: Any = None,
    day90_recovery_se: Any = None,
    day90_replicates_n: Any = 6,
    **forbidden_kwargs: Any,
) -> dict[str, Any]:
    """
    Fail-closed Day-90 Gate 0S dual-assay shedding & separated en-bloc recovery validator:
    forbids caller cap overrides (fixed 1e-6), requires distinct specimen units, strict integer
    observed k == 0 with exact independent sensitivity U95 < 1e-6 (N_eff >= 2,995,731), and keeps
    Student-t L95 >= 0.99 physical retrieval distinct from cellular shedding (IIP-0025; DEC-104).
    """
    if forbidden_kwargs:
        return {'gate_pass': False, 'status': 'HOLD_CALLER_CAP_OR_PARAMETER_OVERRIDE_FORBIDDEN'}
    if (
        not isinstance(ddpcr_specimen_unit_id, str)
        or not isinstance(viable_specimen_unit_id, str)
        or not ddpcr_specimen_unit_id.strip()
        or not viable_specimen_unit_id.strip()
        or ddpcr_specimen_unit_id.strip() == viable_specimen_unit_id.strip()
    ):
        return {'gate_pass': False, 'status': 'HOLD_GATE0S_DUAL_SHEDDING_REQUIRES_INDEPENDENT_SPECIMEN_UNITS'}
    dd_res = evaluate_rare_event_count_iip0024(
        gate_id='GATE0S_SHEDDING_DDPCR',
        k=ddpcr_k,
        n=ddpcr_n,
        max_u95_fraction=1e-6,
        require_observed_zero=True,
        eta_spike_recovery_l95=ddpcr_eta_l95,
    )
    if dd_res['status'] != 'PASS':
        return {'gate_pass': False, 'status': dd_res['status'], 'failed_assay': 'ddPCR', 'detail': dd_res}
    vi_res = evaluate_rare_event_count_iip0024(
        gate_id='GATE0S_SHEDDING_VIABLE_OUTGROWTH',
        k=viable_k,
        n=viable_n,
        max_u95_fraction=1e-6,
        require_observed_zero=True,
        eta_spike_recovery_l95=viable_eta_l95,
        f_clonogenic=viable_f_clonogenic,
    )
    if vi_res['status'] != 'PASS':
        return {'gate_pass': False, 'status': vi_res['status'], 'failed_assay': 'viable_outgrowth', 'detail': vi_res}
    rec_l95: float | None = None
    if day90_recovery_mean is not None or day90_recovery_se is not None:
        if (
            not _is_strict_int_iip0024(day90_replicates_n)
            or day90_replicates_n < 6
            or not _is_strict_finite_number_iip0024(day90_recovery_mean)
            or not _is_strict_finite_number_iip0024(day90_recovery_se)
            or float(day90_recovery_se) <= 0.0
        ):
            return {'gate_pass': False, 'status': 'HOLD_INVALID_DAY90_EN_BLOC_RECOVERY_INPUT'}
        t95_r = student_t_ppf_one_sided_95(float(day90_replicates_n - 1))
        rec_l95 = float(day90_recovery_mean) - t95_r * float(day90_recovery_se)
        if rec_l95 < 0.99:
            return {
                'gate_pass': False,
                'status': 'HOLD_S9_EXPLANT_EN_BLOC_RECOVERY_BELOW_99PCT',
                'day90_recovery_lower95': rec_l95,
            }
    return {
        'gate_pass': True,
        'status': 'CLEARED_GATE0S_DUAL_SHEDDING_AND_SEPARATED_RECOVERY',
        'ddpcr_u95': dd_res['u95_binomial'],
        'viable_u95': vi_res['u95_binomial'],
        'day90_recovery_lower95': rec_l95,
    }


# ============================================================================
# Shared Registry-Bound Stage-0 Research & Molecular Guardrail Evaluator
# Decoupled from Bayesian Calibration Admission (IIP-0026; DEC-105)
# ============================================================================
DESIGN_B_FASTA_SHA256 = '4bb20389773c76b2705031ba675487825dfbd94d920a692cde75b7d56af22aa4'
DESIGN_B_C1_C5_PAIR_SHA256 = '82c351af082e0f97361e715c36e2ab5872b8d097efa1c9f59dfda53ea60496ea'

IIP0026_REQUIRED_CORE_S9_OBS_KEYS: tuple[str, ...] = (
    'obs_icasp9_log10_kill',
    'obs_icasp9_log10_kill_se',
    'obs_tp53_ppm1d_bcor_vaf',
    'obs_tp53_ppm1d_bcor_vaf_se',
    'obs_tp53_ppm1d_bcor_duplex_depth',
    'obs_castseq_translocation_freq',
    'obs_castseq_translocation_se',
    'obs_castseq_duplex_or_amplicon_depth',
    'obs_mcb_clonogenic_survivor_count_k',
    'obs_mcb_effective_clonogenic_trials',
)

IIP0026_FORBIDDEN_CALLER_OVERRIDE_KEYS: frozenset[str] = frozenset({
    'lifecycle_state',
    'iip0026_lifecycle_state',
    'execution_supported',
    'override_execution_supported',
    'force_state4_promotion',
    'bypass_hlac_phased_gate',
    'bypass_pe4max_msi_gate',
})

IIP0026_ALLOWED_MOLECULAR_OBS_KEYS: frozenset[str] = frozenset({
    'bound_c1_c5_pair_sha256',
    'donor_seed_gnomad_target_site_wgs_verified',
    'post_edit_per_locus_segmented_baf',
    'per_locus_segmented_baf',
    'obs_post_edit_segmented_baf_chr19q',
    'obs_post_edit_segmented_baf_chr11p15',
    'obs_post_edit_segmented_baf_chr6p21',
    'post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified',
    'design_b_hlac_phased_guide_cfd_max',
    'obs_design_b_hlac_phased_guide_cfd_max',
    'design_b_hla_ab_separated_pass_delivery_verified',
    'design_b_chr6p21_1p41mb_loh_sv_cleared',
    'design_b_hlac_dt9_positive',
    'design_b_hlae_3d12_positive',
    'design_b_hlac_dt9_matched_hlae_only_control_delta_ratio',
    'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio',
    'design_b_hlac_allele_specific_surface_or_phased_mrna_verified',
    'design_b_locus5_5prime_race_coord_bp',
    'obs_design_b_locus5_5prime_race_coord_bp',
    'design_b_locus5_3prime_race_coord_bp',
    'obs_design_b_locus5_3prime_race_coord_bp',
    'design_b_locus5_pdl1_5prime_race_coord_bp',
    'design_b_locus5_race_umi_n_eff',
    'obs_design_b_locus5_race_umi_n_eff',
    'design_b_locus5_aberrant_umi_k',
    'obs_design_b_locus5_aberrant_umi_k',
    'b_sct_ciita_exon5_and_trimer_verified',
    'tyk2_pe_modality',
    'tyk2_phased_haplotype_sap_verified',
    'tyk2_pe4max_pentaplex_msi_shifts_k',
    'obs_tyk2_pe4max_pentaplex_msi_shifts_k',
    'tyk2_pe4max_wgs_repeat_indels_clone',
    'obs_tyk2_pe4max_wgs_repeat_indels_clone',
    'tyk2_pe4max_wgs_repeat_indels_ctrl',
    'obs_tyk2_pe4max_wgs_repeat_indels_ctrl',
    'tyk2_pe4max_coding_or_splice_de_novo_indels',
    'obs_tyk2_pe4max_coding_or_splice_de_novo_indels',
    'schema_era',
    'iip0026_schema_era',
    'adjunct_conditioning_history_arm',
    'adjunct_depletion_mode',
    'adjunct_sa_fhl1_status',
    'adjunct_immune_responder_history',
    'device_assembly_id',
    'recipient_serum_pair_id',
    'intact_serum_sample_id',
    'igg_depleted_serum_sample_id',
    'positive_adcc_control_id',
    'isotype_adcc_control_id',
    'arm2_depleted_sample_id',
    'contemporaneous_arm3_control_sample_id',
    'sister_assembly_id',
    'construct_sha256',
    'post_sort_endocrine_zm447439_regimen',
    'stage7_start_utc',
    'sort_reaggregation_utc',
    'assay_harvest_utc',
    'registered_mixture_id',
    'locked_analysis_plan_id',
    'gate0s_elevated_store_scope_matched_device_lot_id',
    'iip0029_twoclock_record',
    'background_id',
    'mcb_id',
    'donor_seedline_id',
    'seed_line_id',
    'l2_v3_donor_id',
    'donor_id',
    'sample_id',
    'sample_unit_id',
    'replicates_n',
    'sample_n',
    's11_replicate_n',
    'tier6b_presort_lineage_panel',
    'tier6b_clonogenic_risk_unit',
    'tier6b_sorted_lot_observations',
    'iip0031_locus2_record',
    'rnai_specific_benefit_claimed',
    'phasing_verified',
    'cis_verified',
    'wt_inducibility_verified',
    'rnai_attribution_verified',
    'obs_l2_v3_chr16_ha_wgs_pass',
    'obs_l2_v3_chr16_ha_wgs_variant_count',
    'ciita_locus2_genotype_copy_class',
    'rnai_attribution_scramble_control_id',
    'rnai_attribution_scramble_background_id',
    'rnai_attribution_intronless_control_id',
    'rnai_attribution_intronless_background_id',
    'has_intronic_mir_cd276_pvr',
    'syn_5xisre_cd200_in_design_y',
})


def evaluate_pe4max_piecewise_repeat_indel_cap_iip0026(
    c_clone: Any,
    c_ctrl: Any,
    coding_or_splice_de_novo_indels: Any,
) -> str:
    """
    Evaluate held PE4max clonal WGS repeat-indel cap:
      - 0 coding/splice de novo indels,
      - C_clone <= C_ctrl + 5 for C_ctrl <= 25,
      - C_clone <= floor(1.2 * C_ctrl) for C_ctrl > 25 (IIP-0026; DEC-105).
    """
    for v in (c_clone, c_ctrl, coding_or_splice_de_novo_indels):
        if not _is_strict_int_iip0024(v) or v < 0:
            return 'HOLD_MISSING_PE4MAX_REPEAT_INDEL_EVIDENCE'
    if coding_or_splice_de_novo_indels != 0:
        return 'FAIL_PE4MAX_CODING_OR_SPLICE_DE_NOVO_INDEL_DETECTED'
    cap = (c_ctrl + 5) if c_ctrl <= 25 else int(math.floor(1.2 * c_ctrl))
    return 'PASS' if c_clone <= cap else 'FAIL_PE4MAX_REPEAT_INDEL_CAP_EXCEEDED'


def _is_integral_count_iip0026(val: Any) -> bool:
    return _is_strict_int_iip0024(val) or (
        isinstance(val, float) and math.isfinite(val) and val.is_integer()
    )


def _derive_versioned_product_registry_entry_iip0026(
    pdef: dict[str, Any],
) -> dict[str, Any]:
    pid = str(pdef['product_variant_id'])
    exec_sup = bool(pdef.get('execution_supported', False))
    lifecycle_state = str(
        pdef.get(
            'product_lifecycle_state',
            pdef.get(
                'lifecycle_state',
                (
                    'EXECUTION_SUPPORTED_QUALIFIED_PRODUCT'
                    if exec_sup
                    else (
                        'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
                        if pid in ('ISLET1_B2_5LOCUS_BASE', 'ISLET1_B2_5LOCUS_TYK2_PE') or isinstance(pdef.get('staged_locus2_modeb_donor_binding'), dict)
                        else 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD'
                    )
                ),
            ),
        )
    )
    bound_sha: str | None = None
    binding_block = pdef.get('staged_design_b_donor_reference_binding')
    if isinstance(binding_block, dict) and isinstance(binding_block.get('bound_c1_c5_pair_sha256'), str):
        bound_sha = binding_block['bound_c1_c5_pair_sha256'].strip()
    elif pid in ('ISLET1_B2_5LOCUS_BASE', 'ISLET1_B2_5LOCUS_TYK2_PE'):
        bound_sha = DESIGN_B_C1_C5_PAIR_SHA256
    is_b_sct = pid.startswith('ISLET1_B_SCT_')
    is_endogenous_hlac_b = (pdef.get('family') == 'B') and (not is_b_sct)
    is_tyk2_pe = ('TYK2_PE' in pid) or ('Pro1104Ala' in json.dumps(pdef.get('alleles', {})))
    is_design_y = pid.startswith('ISLET1_Y_') or (pdef.get('family') == 'Y')
    return {
        'product_variant_id': pid,
        'family': pdef.get('family'),
        'execution_supported': exec_sup,
        'lifecycle_state': lifecycle_state,
        'bound_c1_c5_pair_sha256': bound_sha,
        'requires_endogenous_hlac_and_locus5_gates': is_endogenous_hlac_b,
        'requires_b_sct_gates': is_b_sct,
        'requires_tyk2_pe_gates': is_tyk2_pe,
        'is_held_design_y': is_design_y,
        'required_baf_loci': ('chr19q', 'chr11p15', 'chr6p21') if is_endogenous_hlac_b else ('chr19q', 'chr11p15'),
    }


_IIP0026_DESIGNB_TYK2_ALIAS_PAIRS: tuple[tuple[str, str, bool], ...] = (
    ('design_b_hlac_phased_guide_cfd_max', 'obs_design_b_hlac_phased_guide_cfd_max', False),
    ('design_b_hlac_dt9_matched_hlae_only_control_delta_ratio', 'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio', False),
    ('design_b_locus5_5prime_race_coord_bp', 'obs_design_b_locus5_5prime_race_coord_bp', True),
    ('design_b_locus5_3prime_race_coord_bp', 'obs_design_b_locus5_3prime_race_coord_bp', True),
    ('design_b_locus5_race_umi_n_eff', 'obs_design_b_locus5_race_umi_n_eff', True),
    ('design_b_locus5_aberrant_umi_k', 'obs_design_b_locus5_aberrant_umi_k', True),
    ('tyk2_pe4max_pentaplex_msi_shifts_k', 'obs_tyk2_pe4max_pentaplex_msi_shifts_k', True),
    ('tyk2_pe4max_wgs_repeat_indels_clone', 'obs_tyk2_pe4max_wgs_repeat_indels_clone', True),
    ('tyk2_pe4max_wgs_repeat_indels_ctrl', 'obs_tyk2_pe4max_wgs_repeat_indels_ctrl', True),
    ('tyk2_pe4max_coding_or_splice_de_novo_indels', 'obs_tyk2_pe4max_coding_or_splice_de_novo_indels', True),
)


def _baf_locus_discordant_iip0031(val: Any) -> bool:
    """Exact Decimal check for |BAF - 0.50| > 0.05 (inclusive at 0.45 and 0.55; IIP-0031; DEC-110)."""
    if not _is_strict_finite_number_iip0024(val):
        return False
    d_val = Decimal(str(val))
    if d_val < Decimal('0') or d_val > Decimal('1'):
        return False
    return abs(d_val - Decimal('0.50')) > Decimal('0.05')


def _extract_per_locus_baf_dict_iip0026(payload: dict[str, Any]) -> dict[str, Any] | None:
    raw_dict = payload.get('post_edit_per_locus_segmented_baf', payload.get('per_locus_segmented_baf'))
    if (
        'post_edit_per_locus_segmented_baf' in payload
        and 'per_locus_segmented_baf' in payload
        and payload['post_edit_per_locus_segmented_baf'] != payload['per_locus_segmented_baf']
    ):
        return None
    loci_pairs = (
        ('chr19q', 'obs_post_edit_segmented_baf_chr19q'),
        ('chr11p15', 'obs_post_edit_segmented_baf_chr11p15'),
        ('chr6p21', 'obs_post_edit_segmented_baf_chr6p21'),
    )
    if isinstance(raw_dict, dict):
        for loc, ok in loci_pairs:
            if ok in payload:
                ov = payload[ok]
                if ov is None or ov == 'UNAVAILABLE' or (loc in raw_dict and raw_dict[loc] != ov):
                    return None
        return dict(raw_dict)
    if raw_dict is not None:
        return None
    assembled: dict[str, Any] = {}
    for loc, ok in loci_pairs:
        if ok in payload and payload[ok] is not None and payload[ok] != 'UNAVAILABLE':
            assembled[loc] = payload[ok]
    return assembled if assembled else None


def _reconcile_iip0026_aliases_and_baf_iip0031(
    containers: tuple[dict[str, Any], ...] | list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, list[Any]], dict[str, Any] | None, dict[str, list[Any]], bool]:
    """
    Canonicalize all 10 Design-B/TYK2 alias pairs and all BAF representations across containers.
    Never lets a preferred 0-valued spelling shadow a positive adverse spelling, and marks any
    alias/null/type disagreement as `has_conflict=True` while retaining all valid candidates.
    """
    clean_aliases: dict[str, Any] = {}
    alias_cands: dict[str, list[Any]] = {}
    has_conflict = False

    def _val_sig(x: Any) -> str | None:
        if type(x) is bool:
            return f'bool:{x}'
        if _is_strict_int_iip0024(x):
            if abs(x).bit_length() > 53:
                return None
            return f'num:{Decimal(x).normalize()}'
        if _is_strict_finite_number_iip0024(x):
            return f'num:{Decimal(str(x)).normalize()}'
        try:
            return 'json:' + json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False)
        except (TypeError, ValueError, OverflowError, RecursionError):
            return None

    for short_k, obs_k, req_int in _IIP0026_DESIGNB_TYK2_ALIAS_PAIRS:
        raw_vals: list[Any] = []
        for c in containers:
            if not isinstance(c, dict):
                continue
            for ak in (obs_k, short_k):
                if ak in c:
                    raw_vals.append(c[ak])
        if not raw_vals:
            continue
        has_null = any(v is None or v == 'UNAVAILABLE' for v in raw_vals)
        non_null = [v for v in raw_vals if v is not None and v != 'UNAVAILABLE']
        valid_cands: list[Any] = []
        seen_sigs: set[str] = set()
        bad_val = False
        for v in non_null:
            if req_int:
                if not _is_strict_int_iip0024(v) or abs(v).bit_length() > 53:
                    bad_val = True
                    continue
            elif not _is_strict_finite_number_iip0024(v):
                bad_val = True
                continue
            sig = _val_sig(v)
            if sig is None:
                bad_val = True
                continue
            if sig not in seen_sigs:
                seen_sigs.add(sig)
                valid_cands.append(int(v) if req_int else float(v))
        if valid_cands:
            alias_cands[obs_k] = valid_cands
        if bad_val or (has_null and bool(non_null)) or len(valid_cands) > 1:
            has_conflict = True
        elif len(valid_cands) == 1:
            clean_aliases[obs_k] = valid_cands[0]
            clean_aliases[short_k] = valid_cands[0]

    loci_pairs = (
        ('chr19q', 'obs_post_edit_segmented_baf_chr19q'),
        ('chr11p15', 'obs_post_edit_segmented_baf_chr11p15'),
        ('chr6p21', 'obs_post_edit_segmented_baf_chr6p21'),
    )
    baf_locus_cands: dict[str, list[Any]] = {loc: [] for loc, _ in loci_pairs}
    baf_bad = False
    baf_any_present = False
    locus_raw: dict[str, list[Any]] = {loc: [] for loc, _ in loci_pairs}
    for c in containers:
        if not isinstance(c, dict):
            continue
        for map_k in ('post_edit_per_locus_segmented_baf', 'per_locus_segmented_baf'):
            if map_k in c:
                baf_any_present = True
                mv = c[map_k]
                if not isinstance(mv, dict):
                    baf_bad = True
                else:
                    for mk_sub, mv_sub in mv.items():
                        if not isinstance(mk_sub, str):
                            baf_bad = True
                        elif mk_sub in locus_raw:
                            locus_raw[mk_sub].append(mv_sub)
        for loc, ok in loci_pairs:
            if ok in c:
                baf_any_present = True
                locus_raw[loc].append(c[ok])

    clean_baf: dict[str, Any] = {}
    for loc, _ in loci_pairs:
        lvals = locus_raw[loc]
        if not lvals:
            continue
        l_null = any(v is None or v == 'UNAVAILABLE' for v in lvals)
        l_non = [v for v in lvals if v is not None and v != 'UNAVAILABLE']
        l_valid: list[Any] = []
        l_sigs: set[str] = set()
        for v in l_non:
            if not _is_strict_finite_number_iip0024(v) or not (0.0 <= float(v) <= 1.0):
                baf_bad = True
                continue
            sig = _val_sig(v)
            if sig is None:
                baf_bad = True
                continue
            if sig not in l_sigs:
                l_sigs.add(sig)
                l_valid.append(v)
        baf_locus_cands[loc] = l_valid
        if (l_null and bool(l_non)) or len(l_valid) > 1:
            baf_bad = True
        elif len(l_valid) == 1:
            clean_baf[loc] = l_valid[0]

    if baf_bad:
        has_conflict = True
        clean_baf_out: dict[str, Any] | None = None
    else:
        clean_baf_out = clean_baf if (baf_any_present and clean_baf) else None

    return clean_aliases, alias_cands, clean_baf_out, baf_locus_cands, has_conflict


def _evaluate_iip0026_candidate_fails_iip0031(
    reg: dict[str, Any],
    containers: tuple[dict[str, Any], ...] | list[dict[str, Any]],
    alias_cands: dict[str, list[Any]],
    baf_locus_cands: dict[str, list[Any]],
) -> list[str]:
    """Evaluate all non-null candidate values across containers so adverse Design-B/TYK2/BAF stops are never shadowed."""
    cand_fails: list[str] = []

    def _add(code: str) -> None:
        if code not in cand_fails:
            cand_fails.append(code)

    for loc in reg['required_baf_loci']:
        for bv in baf_locus_cands.get(loc, ()):
            if _baf_locus_discordant_iip0031(bv):
                _add('FAIL_DISCORDANT_PER_LOCUS_SEGMENTED_BAF')

    def _any_bool_true(k: str) -> bool:
        return any(isinstance(c, dict) and c.get(k) is True for c in containers)

    def _all_bool_not_true(k: str) -> bool:
        return not _any_bool_true(k)

    if reg['requires_endogenous_hlac_and_locus5_gates']:
        if _any_bool_true('design_b_hlac_dt9_positive') and _any_bool_true('design_b_hlae_3d12_positive') and _all_bool_not_true('design_b_hlac_allele_specific_surface_or_phased_mrna_verified'):
            for dv in alias_cands.get('obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio', ()):
                if 0.0 < float(dv) < 2.0:
                    _add('FAIL_HLAC_SURFACE_LOSS_MASKED_BY_DT9_HLAE_CROSSREACTIVITY')
        max_donor_bp = 4078 if reg['bound_c1_c5_pair_sha256'] == DESIGN_B_C1_C5_PAIR_SHA256 else 10000
        for _rc_c in containers:
            if not isinstance(_rc_c, dict):
                continue
            for _rk_chk in (
                'design_b_locus5_5prime_race_coord_bp', 'obs_design_b_locus5_5prime_race_coord_bp',
                'design_b_locus5_3prime_race_coord_bp', 'obs_design_b_locus5_3prime_race_coord_bp',
            ):
                if _rk_chk in _rc_c and _rc_c[_rk_chk] is not None and _rc_c[_rk_chk] != 'UNAVAILABLE':
                    _rv_chk = _rc_c[_rk_chk]
                    if not _is_strict_int_iip0024(_rv_chk) or _rv_chk <= 0 or _rv_chk > max_donor_bp:
                        _add('FAIL_INVALID_RACE_COORDINATES')
        r5_list = alias_cands.get('obs_design_b_locus5_5prime_race_coord_bp', [])
        r3_list = alias_cands.get('obs_design_b_locus5_3prime_race_coord_bp', [])
        for r5 in r5_list:
            for r3 in r3_list:
                if not (0 < int(r5) < int(r3) <= max_donor_bp):
                    _add('FAIL_INVALID_RACE_COORDINATES')
        for c in containers:
            if isinstance(c, dict) and 'design_b_locus5_pdl1_5prime_race_coord_bp' in c:
                pdl1_v = c['design_b_locus5_pdl1_5prime_race_coord_bp']
                if pdl1_v is not None and pdl1_v != 'UNRESOLVED':
                    if not _is_strict_int_iip0024(pdl1_v) or pdl1_v <= 0 or pdl1_v > 2381:
                        _add('FAIL_INVALID_RACE_COORDINATES')
                    elif pdl1_v <= 2308:
                        _add('FAIL_MATURE_TRANSCRIPT_INCLUDES_UPSTREAM_OUT_OF_FRAME_KOZAK_ORF')
        for uk in alias_cands.get('obs_design_b_locus5_aberrant_umi_k', ()):
            if _is_strict_int_iip0024(uk) and uk > 0:
                _add('FAIL_DESIGNB_LOCUS5_CHIMERIC_TRANSCRIPT_DETECTED')

    if reg['requires_tyk2_pe_gates']:
        has_pe4 = any(isinstance(c, dict) and c.get('tyk2_pe_modality') == 'PE4MAX_MSI_GATED' for c in containers)
        if has_pe4 and _any_bool_true('tyk2_phased_haplotype_sap_verified'):
            for mk in alias_cands.get('obs_tyk2_pe4max_pentaplex_msi_shifts_k', ()):
                if _is_strict_int_iip0024(mk) and mk > 0:
                    _add('FAIL_MSI_DETECTED')
            for cd in alias_cands.get('obs_tyk2_pe4max_coding_or_splice_de_novo_indels', ()):
                if _is_strict_int_iip0024(cd) and cd > 0:
                    _add('FAIL_PE4MAX_CODING_OR_SPLICE_DE_NOVO_INDEL_DETECTED')
            for cc in alias_cands.get('obs_tyk2_pe4max_wgs_repeat_indels_clone', ()):
                for ct in alias_cands.get('obs_tyk2_pe4max_wgs_repeat_indels_ctrl', ()):
                    if _is_strict_int_iip0024(cc) and cc >= 0 and _is_strict_int_iip0024(ct) and ct >= 0:
                        cap = (ct + 5) if ct <= 25 else int(math.floor(1.2 * ct))
                        if cc > cap:
                            _add('FAIL_PE4MAX_REPEAT_INDEL_CAP_EXCEEDED')
    return cand_fails


_S11_MEAN_SE_PAIRS_IIP0031: tuple[tuple[str, str], ...] = (
    ('obs_l2_v3_ap_c5b9_pct', 'obs_l2_v3_ap_c5b9_se'),
    ('obs_l2_v3_classical_c4b2a_lysis_pct', 'obs_l2_v3_classical_c4b2a_lysis_se'),
    ('obs_l2_v3_lectin_c4b2a_lysis_pct', 'obs_l2_v3_lectin_c4b2a_lysis_se'),
    ('obs_l2_v3_overall_20pct_nhs_lysis_pct', 'obs_l2_v3_overall_20pct_nhs_lysis_se'),
    ('obs_l2_v3_intact_serum_adcc_lysis_pct', 'obs_l2_v3_intact_serum_adcc_lysis_se'),
    ('obs_l2_v3_igg_depleted_serum_adcc_lysis_pct', 'obs_l2_v3_igg_depleted_serum_adcc_lysis_se'),
    ('obs_l2_v3_igg_increment_adcc_lysis_pct', 'obs_l2_v3_igg_increment_adcc_lysis_se'),
    ('obs_l2_v3_exact_stack_effero_abs_pct', 'obs_l2_v3_exact_stack_effero_abs_se'),
    ('obs_l2_v3_exact_stack_effero_pct_of_wt', 'obs_l2_v3_exact_stack_effero_pct_of_wt_se'),
    ('obs_l2_v3_pericyte_coverage_pct', 'obs_l2_v3_pericyte_coverage_se'),
    ('obs_l2_v3_plvap_iec_network_t50_h', 'obs_l2_v3_plvap_iec_network_t50_se'),
    ('obs_l2_v3_plvap_iec_viability_pct', 'obs_l2_v3_plvap_iec_viability_se'),
    ('obs_l2_v3_tandem_polya_readthrough_fold', 'obs_l2_v3_tandem_polya_readthrough_se'),
    ('obs_l2_v3_cfh_stp_gpi_quantibrite_molecules', 'obs_l2_v3_cfh_stp_gpi_quantibrite_se'),
    ('obs_l2_v3_unit2_copy_expression_quantibrite_or_fold', 'obs_l2_v3_unit2_copy_expression_se'),
    ('obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_se'),
    ('obs_ciita_allele1_ex4_ex6_aberrant_skip_pct', 'obs_ciita_allele1_ex4_ex6_aberrant_skip_se'),
    ('obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_se'),
    ('obs_ciita_allele2_ex4_ex6_aberrant_skip_pct', 'obs_ciita_allele2_ex4_ex6_aberrant_skip_se'),
    ('obs_ciita_ki_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ki_ex4_ex6_aberrant_skip_se'),
    ('obs_ciita_ko_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ko_ex4_ex6_aberrant_skip_se'),
    ('obs_edited_72h_ifng_hla_dr_dp_dq_pos_pct', 'obs_edited_72h_ifng_hla_dr_dp_dq_pos_se'),
    ('obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_pct', 'obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_se'),
    ('obs_l2_v3_stage0_eef1a1_cpg_methylation_pct', 'obs_l2_v3_stage0_eef1a1_cpg_methylation_se'),
    ('obs_l2_v3_stage0_chs4_cpg_methylation_pct', 'obs_l2_v3_stage0_chs4_cpg_methylation_se'),
    ('obs_l2_v3_stage0_dm1_cpg_methylation_pct', 'obs_l2_v3_stage0_dm1_cpg_methylation_se'),
    ('obs_l2_v3_stage7_eef1a1_cpg_methylation_pct', 'obs_l2_v3_stage7_eef1a1_cpg_methylation_se'),
    ('obs_l2_v3_stage7_chs4_cpg_methylation_pct', 'obs_l2_v3_stage7_chs4_cpg_methylation_se'),
    ('obs_l2_v3_stage7_dm1_cpg_methylation_pct', 'obs_l2_v3_stage7_dm1_cpg_methylation_se'),
    ('obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment_se'),
    ('obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment_se'),
    ('obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment_se'),
    ('obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment_se'),
    ('obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_pct', 'obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_se'),
    ('obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment', 'obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment_se'),
    ('obs_l2_v3_manf_locus3_nfkb_suppression_ratio', 'obs_l2_v3_manf_locus3_nfkb_suppression_se'),
    ('obs_l2_v3_manf_ki67_prolif_pct', 'obs_l2_v3_manf_ki67_prolif_se'),
    ('obs_l2_v3_manf_ocr_ratio_vs_ctrl', 'obs_l2_v3_manf_ocr_ratio_vs_ctrl_se'),
    ('obs_l2_v3_ccl22_treg_to_effector_chemotaxis_ratio', 'obs_l2_v3_ccl22_treg_to_effector_chemotaxis_se'),
    ('obs_l2_v3_ccl22_effector_chemotaxis_fold_ctrl', 'obs_l2_v3_ccl22_effector_chemotaxis_fold_se'),
    ('obs_l2_v3_small_rnaseq_mir375_retention_pct', 'obs_l2_v3_small_rnaseq_mir375_retention_se'),
    ('obs_l2_v3_intron_splicing_efficiency_pct', 'obs_l2_v3_intron_splicing_efficiency_se'),
    ('obs_l2_v3_cd276_kd_pct', 'obs_l2_v3_cd276_kd_se'),
    ('obs_l2_v3_pvr_kd_pct', 'obs_l2_v3_pvr_kd_se'),
    ('obs_l2_v3_joint_stratum_nk_lysis_pct', 'obs_l2_v3_joint_stratum_nk_lysis_se'),
)


def _merge_stage0_container_fields_iip0031(
    containers: tuple[dict[str, Any], ...] | list[dict[str, Any]],
    *,
    excluded_keys: frozenset[str] = frozenset(),
    key_predicate: Any = None,
) -> tuple[dict[str, Any], list[dict[str, Any]], bool]:
    """
    Conflict-aware, null-safe Stage-0 field merger shared across Direct, Serialized Ingestion,
    and Chronological Replay (IIP-0031; DEC-110).
    - Never lets an explicit `None` / `'UNAVAILABLE'` overwrite a non-null measurement.
    - Flags ANY key where `None`/`'UNAVAILABLE'` co-occurs with a non-null value OR where
      multiple containers supply distinct non-null values (strict JSON wire-type signature),
      OR where any container contains a non-string dictionary key.
    - Builds strictly per-container S9 probe candidates (NEVER cross-synthesizing `mean` from
      container A with `se` from container B when their `mean`s disagree, while pairing a bare
      repeated `mean` with an explicit valid `se` from another container that reported the exact
      same `mean`).
    """
    all_keys: list[str] = []
    seen_k: set[str] = set()
    has_conflict = False
    for c in containers:
        if not isinstance(c, dict):
            continue
        for raw_k in c:
            if not isinstance(raw_k, str):
                has_conflict = True
                continue
            ks = raw_k
            if ks in excluded_keys or ks in seen_k:
                continue
            if key_predicate is not None and not key_predicate(ks):
                continue
            seen_k.add(ks)
            all_keys.append(ks)
    clean_merged: dict[str, Any] = {}
    conflict_cands_by_key: dict[str, list[Any]] = {}
    for k in all_keys:
        vals = [c[k] for c in containers if isinstance(c, dict) and k in c]
        has_null = any(v is None or v == 'UNAVAILABLE' for v in vals)
        non_null = [v for v in vals if v is not None and v != 'UNAVAILABLE']
        distinct_non_null: list[Any] = []
        seen_sigs: set[str] = set()
        malformed_val = False
        for v in non_null:
            if _is_strict_int_iip0024(v) and abs(v).bit_length() > 53:
                malformed_val = True
                continue
            try:
                sv = json.dumps(v, sort_keys=True, separators=(',', ':'), allow_nan=False)
            except (TypeError, ValueError, OverflowError, RecursionError):
                malformed_val = True
                distinct_non_null.append(v)
                continue
            if sv not in seen_sigs:
                seen_sigs.add(sv)
                distinct_non_null.append(v)
        if malformed_val or (has_null and bool(non_null)) or len(distinct_non_null) > 1:
            has_conflict = True
            if distinct_non_null:
                conflict_cands_by_key[k] = distinct_non_null
        elif distinct_non_null:
            clean_merged[k] = distinct_non_null[0]

    _all_mean_se_pairs = (
        ('obs_icasp9_log10_kill', 'obs_icasp9_log10_kill_se'),
        ('obs_tp53_ppm1d_bcor_vaf', 'obs_tp53_ppm1d_bcor_vaf_se'),
        ('obs_castseq_translocation_freq', 'obs_castseq_translocation_se'),
        ('obs_chop_coincidence_fold', 'obs_chop_coincidence_se'),
        ('obs_gate0e_chop_log2fc', 'obs_gate0e_chop_log2fc_se'),
        ('obs_hc10_surface_fold_wt', 'obs_hc10_surface_fold_se'),
        ('obs_gate0e_hc10_surface_mfi_ratio', 'obs_gate0e_hc10_surface_mfi_se'),
        *_S11_MEAN_SE_PAIRS_IIP0031,
    )
    probe_candidates: list[dict[str, Any]] = []
    if conflict_cands_by_key:
        for c in containers:
            if not isinstance(c, dict):
                continue
            cand_c = dict(clean_merged)
            c_contributed = False
            for k in conflict_cands_by_key:
                if k in c and c[k] is not None and c[k] != 'UNAVAILABLE':
                    if _is_strict_int_iip0024(c[k]) and abs(c[k]).bit_length() > 53:
                        continue
                    cand_c[k] = c[k]
                    c_contributed = True
            if not c_contributed:
                continue
            # Never cross-breed mean from container A with SE from container B when means disagree,
            # and do not impute a default SE onto a bare repeated mean when another container already
            # supplied an explicit valid SE for that exact same mean (or a mean-free shared SE) (Finding #9):
            for mk, sk in _all_mean_se_pairs:
                if mk not in conflict_cands_by_key and sk not in conflict_cands_by_key:
                    continue
                if mk in cand_c and sk not in c:
                    if _is_strict_finite_number_iip0024(cand_c[mk]):
                        m_f = float(cand_c[mk])
                        matching_ses = [
                            c_other[sk]
                            for c_other in containers
                            if isinstance(c_other, dict)
                            and mk in c_other
                            and _is_strict_finite_number_iip0024(c_other[mk])
                            and float(c_other[mk]) == m_f
                            and sk in c_other
                            and _is_strict_finite_number_iip0024(c_other[sk])
                            and float(c_other[sk]) > 0.0
                        ]
                        mean_free_ses = [
                            c_other[sk]
                            for c_other in containers
                            if isinstance(c_other, dict)
                            and mk not in c_other
                            and sk in c_other
                            and _is_strict_finite_number_iip0024(c_other[sk])
                            and float(c_other[sk]) > 0.0
                        ]
                        if matching_ses:
                            cand_c[sk] = max(matching_ses, key=float)
                        elif mean_free_ses:
                            cand_c[sk] = max(mean_free_ses, key=float)
                        else:
                            cand_c.pop(sk, None)
                    else:
                        cand_c.pop(sk, None)
                elif sk in cand_c and mk not in cand_c and mk in conflict_cands_by_key:
                    # Standalone SE without its mean in a container when mean is conflicted across containers:
                    # drop orphan SE so it cannot attach to a non-existent mean.
                    cand_c.pop(sk, None)
            probe_candidates.append(cand_c)
    return clean_merged, probe_candidates, has_conflict


def _eval_stage0_tree_preserving_s9_iip0031(
    safe_obs: dict[str, Any],
    *,
    root_path: Path,
) -> tuple[dict[str, Any], bool]:
    """
    Evaluate `evaluate_stage0_decision_tree` while isolating unrelated co-supported Stage-0 alias
    collisions (e.g. `obs_regimeII_lysis_2h_pct` vs `obs_gate0c_2h_lysis_pct`, CHOP fold vs `2**log2FC`,
    or SERPINB9 `k_ass` vs `10**log10_kass`) so valid S9 hard stops ALWAYS survive without `ValueError`
    while any alias collision forces a non-clearing `HOLD_CONFLICTING_STAGE0_MEASUREMENT` (Finding #18).
    """
    try:
        return evaluate_stage0_decision_tree(safe_obs, root=root_path), False
    except ValueError:
        pass

    conflict_pairs: list[tuple[str, str]] = []
    for ak, ck in STAGE0_JSON_TRIGGER_KEY_ALIASES.items():
        if ak in safe_obs and ck in safe_obs:
            conflict_pairs.append((ak, ck))
    for k1, k2 in CO_SUPPORTED_STAGE0_EQUIVALENT_KEY_PAIRS:
        if k1 in safe_obs and k2 in safe_obs:
            conflict_pairs.append((k1, k2))
    for k1, k2 in (
        ('obs_chop_coincidence_fold', 'obs_gate0e_chop_log2fc'),
        ('obs_serpinb9_kass_M_s', 'obs_serpinb9_log10_kass'),
    ):
        if k1 in safe_obs and k2 in safe_obs:
            conflict_pairs.append((k1, k2))

    # Strip non-S9 conflicted pairs and probe both sides of any S9-relevant conflicted pair (CHOP / HC10):
    s9_pair_keys = {
        ('obs_chop_coincidence_fold', 'obs_gate0e_chop_log2fc'),
        ('obs_hc10_surface_fold_wt', 'obs_gate0e_hc10_surface_mfi_ratio'),
    }
    stripped = dict(safe_obs)
    active_s9_pairs: list[tuple[str, str]] = []
    for k1, k2 in conflict_pairs:
        if (k1, k2) in s9_pair_keys or (k2, k1) in s9_pair_keys:
            active_s9_pairs.append((k1, k2))
        else:
            stripped.pop(k1, None)
            stripped.pop(k2, None)

    variants: list[dict[str, Any]] = [stripped]
    for k1, k2 in active_s9_pairs:
        next_v: list[dict[str, Any]] = []
        for d in variants:
            d1 = dict(d)
            d1.pop(k2, None)
            d2 = dict(d)
            d2.pop(k1, None)
            next_v.extend((d1, d2))
        variants = next_v

    fallback_tree: dict[str, Any] | None = None
    for vd in variants:
        try:
            te = evaluate_stage0_decision_tree(vd, root=root_path)
        except ValueError:
            continue
        if te['safety_kill_s9_triggered']:
            return te, True
        if fallback_tree is None:
            fallback_tree = te
    if fallback_tree is not None:
        return fallback_tree, True
    # Final minimal S9-only projection if another non-S9 metadata field caused ValueError:
    s9_only_keys = set(IIP0026_REQUIRED_CORE_S9_OBS_KEYS) | set(
        k for ks in S9_REPLICATE_KEYS.values() for k in ks
    ) | {
        'product_variant_id', 'current_lot_id', 'assay_campaign_id',
        'obs_tp53_ppm1d_bcor_mutant_molecules_k', 'obs_castseq_mutant_junctions_k',
        'obs_hazard_state_prolif_survivor_count_k', 'obs_hazard_state_prolif_effective_trials',
        'obs_ap1903_step1_survivor_count_k', 'obs_ap1903_step1_effective_trials',
        'obs_residual_psc_effective_fraction', 'obs_residual_psc_effective_trials',
    }
    min_obs = {k: v for k, v in safe_obs.items() if k in s9_only_keys}
    return evaluate_stage0_decision_tree(min_obs, root=root_path), True


def _pick_headline_s11_fail_iip0031(fail_codes: list[str] | set[str]) -> str:
    """Deterministic S11 headline failure selector: substantive stops outrank schema-only unrecognized-key errors."""
    substantive = sorted(r for r in fail_codes if r != 'FAIL_S11_UNRECOGNIZED_RECORD_KEYS')
    if substantive:
        return substantive[0]
    return sorted(fail_codes)[0]


def _evaluate_all_s11_containers_iip0031(
    base_env: dict[str, Any],
    containers: tuple[dict[str, Any], ...] | list[dict[str, Any]],
    *,
    root_path: Path,
    product_def: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Evaluate both nested `iip0031_locus2_record` panels AND flat Rule S11 observations across
    all supported Stage-0 containers (`top`, `scope`, `iip0026_molecular_observations`,
    `measurements`) without letting a `None`, malformed, or conflicting sibling erase an
    independently attributable `FAIL_S11_*` stop (Findings #3, #16, #46; Patch v11).
    """
    env_clean = {
        k: v for k, v in base_env.items()
        if isinstance(k, str) and k not in ('scope', 'measurements', 'iip0026_molecular_observations', 'iip0031_locus2_record')
    }
    present_vals: list[Any] = [
        c['iip0031_locus2_record']
        for c in containers
        if isinstance(c, dict) and 'iip0031_locus2_record' in c
    ]

    _s11_rep_alias_keys = ('obs_stage0_independent_lot_replicates_n', 'replicates_n', 's11_replicate_n')

    def _is_s11_domain_key(k: str) -> bool:
        return bool(
            k.startswith(('obs_l2_v3_', 'obs_ciita_', 'obs_edited_72h_', 'obs_wt_lineage_', 'obs_cis_', 'obs_16p13_'))
            or k in _ALL_S11_RECOGNIZED_EVIDENCE_KEYS
            or k in ('l2_v3_donor_id', 'donor_id')
        )

    def _is_s11_flat_key(k: str) -> bool:
        return bool(
            _is_s11_domain_key(k)
            or k in _s11_rep_alias_keys
            or k == 'background_id'
        )

    flat_s11_keys: list[str] = []
    seen_flat_k: set[str] = set()
    any_s11_domain_key_present = False
    for c in containers:
        if not isinstance(c, dict):
            continue
        for k in c:
            if not isinstance(k, str) or k in seen_flat_k:
                continue
            if _is_s11_flat_key(k):
                seen_flat_k.add(k)
                flat_s11_keys.append(k)
                if _is_s11_domain_key(k):
                    any_s11_domain_key_present = True

    flat_s11_conflicted_keys: set[str] = set()
    for k in flat_s11_keys:
        vals = [c[k] for c in containers if isinstance(c, dict) and k in c]
        has_null = any(v is None or v == 'UNAVAILABLE' for v in vals)
        non_null = [v for v in vals if v is not None and v != 'UNAVAILABLE']
        seen_sigs: set[str] = set()
        malformed_val = False
        for v in non_null:
            if _is_strict_int_iip0024(v) and abs(v).bit_length() > 53:
                malformed_val = True
                continue
            try:
                sv = json.dumps(v, sort_keys=True, separators=(',', ':'), allow_nan=False)
            except (TypeError, ValueError, OverflowError, RecursionError):
                malformed_val = True
                continue
            seen_sigs.add(sv)
        if malformed_val or (has_null and bool(non_null)) or len(seen_sigs) > 1:
            flat_s11_conflicted_keys.add(k)

    rep_vals_all = [
        c[rk]
        for c in containers
        if isinstance(c, dict)
        for rk in _s11_rep_alias_keys
        if rk in c
    ]
    rep_has_null = any(v is None or v == 'UNAVAILABLE' for v in rep_vals_all)
    rep_non_null = [v for v in rep_vals_all if v is not None and v != 'UNAVAILABLE']
    rep_sigs: set[str] = set()
    rep_malformed = False
    for rv in rep_non_null:
        if not _is_strict_int_iip0024(rv) or abs(rv).bit_length() > 53:
            rep_malformed = True
            continue
        rep_sigs.add(str(rv))
    s11_rep_cross_conflict = bool(
        rep_vals_all and (rep_malformed or (rep_has_null and bool(rep_non_null)) or len(rep_sigs) > 1)
    )

    is_l2_product_def = bool(
        isinstance(product_def, dict) and isinstance(product_def.get('staged_locus2_modeb_donor_binding'), dict)
    )
    has_flat_s11_domain_conflict = any(_is_s11_domain_key(k) for k in flat_s11_conflicted_keys)
    has_flat_s11_conflict = bool(
        has_flat_s11_domain_conflict
        or (s11_rep_cross_conflict and (is_l2_product_def or any_s11_domain_key_present or bool(present_vals)))
        or ('background_id' in flat_s11_conflicted_keys and (is_l2_product_def or any_s11_domain_key_present or bool(present_vals)))
    )

    flat_probe_envs: list[dict[str, Any]] = []
    if has_flat_s11_conflict:
        for c in containers:
            if not isinstance(c, dict):
                continue
            cand = dict(env_clean)
            c_contributed = False
            for k in flat_s11_conflicted_keys:
                if k in c and c[k] is not None and c[k] != 'UNAVAILABLE':
                    if _is_strict_int_iip0024(c[k]) and abs(c[k]).bit_length() > 53:
                        continue
                    cand[k] = c[k]
                    c_contributed = True
            c_rep_keys = [rk for rk in _s11_rep_alias_keys if rk in c and c[rk] is not None and c[rk] != 'UNAVAILABLE']
            if s11_rep_cross_conflict and c_rep_keys:
                for rk in _s11_rep_alias_keys:
                    cand.pop(rk, None)
                for rk in c_rep_keys:
                    if not (_is_strict_int_iip0024(c[rk]) and abs(c[rk]).bit_length() > 53):
                        cand[rk] = c[rk]
                        c_contributed = True
            if not c_contributed:
                continue
            # Never cross-breed mean from container A with SE from container B when means disagree:
            for mk, sk in _S11_MEAN_SE_PAIRS_IIP0031:
                if mk not in flat_s11_conflicted_keys and sk not in flat_s11_conflicted_keys:
                    continue
                if mk in cand and sk not in c:
                    if _is_strict_finite_number_iip0024(cand[mk]):
                        m_f = float(cand[mk])
                        matching_ses = [
                            c_o[sk]
                            for c_o in containers
                            if isinstance(c_o, dict)
                            and mk in c_o
                            and _is_strict_finite_number_iip0024(c_o[mk])
                            and float(c_o[mk]) == m_f
                            and sk in c_o
                            and _is_strict_finite_number_iip0024(c_o[sk])
                            and float(c_o[sk]) > 0.0
                        ]
                        mean_free_ses = [
                            c_o[sk]
                            for c_o in containers
                            if isinstance(c_o, dict)
                            and mk not in c_o
                            and sk in c_o
                            and _is_strict_finite_number_iip0024(c_o[sk])
                            and float(c_o[sk]) > 0.0
                        ]
                        if matching_ses:
                            cand[sk] = max(matching_ses, key=float)
                        elif mean_free_ses:
                            cand[sk] = max(mean_free_ses, key=float)
                        else:
                            cand.pop(sk, None)
                    else:
                        cand.pop(sk, None)
                elif sk in cand and mk not in cand and mk in flat_s11_conflicted_keys:
                    cand.pop(sk, None)
            flat_probe_envs.append(cand)

    if not present_vals and not flat_probe_envs:
        return evaluate_locus2_ciita_modeb_challenger_v3_gate(env_clean, root=root_path, product_def=product_def)

    has_null_or_non_dict = any(v is None or not isinstance(v, dict) for v in present_vals)
    dict_panels: list[dict[str, Any]] = []
    seen_panel_sigs: set[str] = set()
    for v in present_vals:
        if isinstance(v, dict):
            try:
                sig = json.dumps(v, sort_keys=True, separators=(',', ':'), allow_nan=False)
            except (TypeError, ValueError, OverflowError, RecursionError):
                has_null_or_non_dict = True
                dict_panels.append(v)
                continue
            if sig not in seen_panel_sigs:
                seen_panel_sigs.add(sig)
                dict_panels.append(v)
    has_conflicting_dicts = len(dict_panels) > 1

    evals: list[dict[str, Any]] = []
    if dict_panels:
        for dp in dict_panels:
            evals.append(
                evaluate_locus2_ciita_modeb_challenger_v3_gate(
                    {**env_clean, 'iip0031_locus2_record': dp},
                    root=root_path,
                    product_def=product_def,
                )
            )
            for f_env in flat_probe_envs:
                evals.append(
                    evaluate_locus2_ciita_modeb_challenger_v3_gate(
                        {**f_env, 'iip0031_locus2_record': dp},
                        root=root_path,
                        product_def=product_def,
                    )
                )
    elif present_vals:
        # Only null / non-dict values supplied: still probe flat S11 fields for architectural/quantitative stops,
        # then enforce HOLD_S11_MALFORMED_RECORD if no FAIL_S11_* fired!
        flat_ev = evaluate_locus2_ciita_modeb_challenger_v3_gate(
            {**env_clean, 'iip0031_locus2_record': '__MALFORMED_NON_DICT__'},
            root=root_path,
            product_def=product_def,
        )
        flat_only_ev = evaluate_locus2_ciita_modeb_challenger_v3_gate(
            env_clean,
            root=root_path,
            product_def=product_def,
        )
        evals.extend((flat_only_ev, flat_ev))
    else:
        evals.append(
            evaluate_locus2_ciita_modeb_challenger_v3_gate(
                env_clean,
                root=root_path,
                product_def=product_def,
            )
        )

    for f_env in flat_probe_envs:
        evals.append(
            evaluate_locus2_ciita_modeb_challenger_v3_gate(
                f_env,
                root=root_path,
                product_def=product_def,
            )
        )

    all_fails: list[str] = []
    all_holds: list[str] = []
    for ev in evals:
        st = str(ev.get('status', ''))
        if st.startswith('FAIL_') and st not in all_fails:
            all_fails.append(st)
        for r in ev.get('reasons', []):
            rs = str(r)
            if rs.startswith('FAIL_'):
                if rs not in all_fails:
                    all_fails.append(rs)
            elif rs not in all_holds:
                all_holds.append(rs)

    has_any_s11_conflict = bool(has_conflicting_dicts or has_flat_s11_conflict)
    if has_any_s11_conflict and 'HOLD_S11_CONFLICTING_NESTED_EVIDENCE' not in all_holds:
        all_holds.insert(0, 'HOLD_S11_CONFLICTING_NESTED_EVIDENCE')
    if has_null_or_non_dict and 'HOLD_S11_MALFORMED_RECORD' not in all_holds:
        all_holds.insert(0, 'HOLD_S11_MALFORMED_RECORD')

    base_out = dict(evals[0])
    if all_fails:
        win_fail = _pick_headline_s11_fail_iip0031(all_fails)
        ordered_fails = [win_fail] + sorted(f for f in all_fails if f != win_fail)
        base_out['status'] = win_fail
        base_out['gate_pass'] = False
        base_out['reasons'] = ordered_fails + all_holds
        return base_out
    if has_any_s11_conflict:
        base_out['status'] = 'HOLD_S11_CONFLICTING_NESTED_EVIDENCE'
        base_out['gate_pass'] = False
        base_out['reasons'] = all_holds
        return base_out
    if has_null_or_non_dict:
        base_out['status'] = 'HOLD_S11_MALFORMED_RECORD'
        base_out['gate_pass'] = False
        base_out['reasons'] = all_holds
        return base_out
    return base_out


def _build_safe_stage0_tree_obs_iip0026(
    payload: dict[str, Any],
    *,
    root_path: Path,
) -> dict[str, Any]:
    """Project valid Stage-0 scalar, boolean, and metadata keys onto evaluate_stage0_decision_tree for branch parity."""
    contract = derive_contract(root_path)
    valid_mean_se: dict[str, tuple[float, float]] = {}
    for e in list(contract['endpoints']) + list(contract.get('stage0_decision_tree_endpoints', [])):
        valid_mean_se[e['meanKey']] = (float(e['min']), float(e['max']))
        valid_mean_se[e['seKey']] = (float(e.get('minimum_se') or contract['minimum_se']), 1e12)
    out: dict[str, Any] = {}
    pid = payload.get('product_variant_id')
    if isinstance(pid, str) and pid.strip() and pid == pid.strip():
        out['product_variant_id'] = pid
    lot_v = payload.get('current_lot_id', payload.get('lot_id'))
    if isinstance(lot_v, str) and lot_v.strip() and len(lot_v.strip()) <= 120:
        out['current_lot_id'] = lot_v.strip()
    camp_v = payload.get('assay_campaign_id', payload.get('campaign_id'))
    if isinstance(camp_v, str) and camp_v.strip() and len(camp_v.strip()) <= 120:
        out['assay_campaign_id'] = camp_v.strip()
    bsha = payload.get('bound_c1_c5_pair_sha256')
    if isinstance(bsha, str) and re.fullmatch(r'[0-9a-f]{64}', bsha.strip()):
        out['bound_c1_c5_pair_sha256'] = bsha.strip()
    pem = payload.get('tyk2_pe_modality')
    if isinstance(pem, str) and pem.strip() in ALLOWED_TYK2_PE_MODALITIES:
        out['tyk2_pe_modality'] = pem.strip()
    sera = payload.get('iip0026_schema_era', payload.get('schema_era'))
    if sera == 'LEGACY_PRE_IIP0026':
        out['iip0026_schema_era'] = 'LEGACY_PRE_IIP0026'
    elif sera == 'IIP0026_V1':
        out['iip0026_schema_era'] = 'IIP0026_V1'
    for enum_k, enum_allowed in (
        ('adjunct_conditioning_history_arm', ALLOWED_ADJUNCT_CONDITIONING_HISTORY_ARMS),
        ('adjunct_depletion_mode', ALLOWED_ADJUNCT_DEPLETION_MODES),
        ('adjunct_sa_fhl1_status', ALLOWED_ADJUNCT_SA_FHL1_STATUSES),
        ('adjunct_immune_responder_history', ALLOWED_ADJUNCT_IMMUNE_RESPONDER_HISTORIES),
    ):
        ev = payload.get(enum_k)
        if isinstance(ev, str) and ev == ev.strip() and ev in enum_allowed:
            out[enum_k] = ev
    for str_id_k in (
        'device_assembly_id',
        'recipient_serum_pair_id',
        'intact_serum_sample_id',
        'igg_depleted_serum_sample_id',
        'positive_adcc_control_id',
        'isotype_adcc_control_id',
        'arm2_depleted_sample_id',
        'contemporaneous_arm3_control_sample_id',
        'sister_assembly_id',
        'construct_sha256',
        'post_sort_endocrine_zm447439_regimen',
        'stage7_start_utc',
        'sort_reaggregation_utc',
        'assay_harvest_utc',
        'registered_mixture_id',
        'locked_analysis_plan_id',
        'gate0s_elevated_store_scope_matched_device_lot_id',
        'background_id',
        'mcb_id',
        'donor_seedline_id',
        'seed_line_id',
        'l2_v3_donor_id',
        'donor_id',
        'sample_id',
        'sample_unit_id',
        'tier6b_presort_lineage_panel',
        'tier6b_clonogenic_risk_unit',
    ):
        sv = payload.get(str_id_k)
        if isinstance(sv, str) and sv == sv.strip() and 2 <= len(sv) <= 120:
            out[str_id_k] = sv
    _safe_s9_n_resolved, _safe_s9_n_conflicts, _ = _resolve_s9_replicate_counts_detailed_iip0031(
        (payload,),
        flag_explicit_none=True,
    )
    for _grp_name, _grp_keys in S9_REPLICATE_KEYS.items():
        if _grp_name in _safe_s9_n_conflicts:
            if _safe_s9_n_resolved.get(_grp_name) is not None:
                out[_grp_keys[0]] = _safe_s9_n_resolved[_grp_name]
            continue
        for _n_k in _grp_keys:
            _nv = payload.get(_n_k)
            if _is_strict_finite_number_iip0024(_nv) and 1.0 <= float(_nv) <= 10000.0:
                out[_n_k] = int(_nv) if type(_nv) is int else float(_nv)
    if 'iip0029_twoclock_record' in payload and payload['iip0029_twoclock_record'] is not None:
        out['iip0029_twoclock_record'] = payload['iip0029_twoclock_record']
    if 'iip0031_locus2_record' in payload and payload['iip0031_locus2_record'] is not None:
        out['iip0031_locus2_record'] = payload['iip0031_locus2_record']
    if 'tier6b_sorted_lot_observations' in payload and payload['tier6b_sorted_lot_observations'] is not None:
        out['tier6b_sorted_lot_observations'] = payload['tier6b_sorted_lot_observations']
    for extra_bool_k in ('rnai_specific_benefit_claimed', 'has_intronic_mir_cd276_pvr', 'syn_5xisre_cd200_in_design_y'):
        if type(payload.get(extra_bool_k)) is bool:
            out[extra_bool_k] = payload[extra_bool_k]
    for bkey in STAGE0_DECISION_TREE_BOOLEAN_KEYS:
        if type(payload.get(bkey)) is bool:
            out[bkey] = payload[bkey]
    baf_map = _extract_per_locus_baf_dict_iip0026(payload)
    if isinstance(baf_map, dict):
        for loc, ok in (
            ('chr19q', 'obs_post_edit_segmented_baf_chr19q'),
            ('chr11p15', 'obs_post_edit_segmented_baf_chr11p15'),
            ('chr6p21', 'obs_post_edit_segmented_baf_chr6p21'),
        ):
            v = baf_map.get(loc)
            if _is_strict_finite_number_iip0024(v) and 0.0 <= float(v) <= 1.0:
                out[ok] = float(v)
    for short_k, obs_k, req_int in _IIP0026_DESIGNB_TYK2_ALIAS_PAIRS:
        raw_v = payload.get(obs_k, payload.get(short_k))
        if raw_v is None:
            continue
        if req_int:
            if not _is_strict_int_iip0024(raw_v) or abs(raw_v).bit_length() > 53:
                continue
        elif not _is_strict_finite_number_iip0024(raw_v):
            continue
        lo, hi = valid_mean_se[obs_k]
        fv = float(raw_v)
        if lo <= fv <= hi:
            out[obs_k] = int(raw_v) if req_int else fv
    for k, v in payload.items():
        if not isinstance(k, str):
            continue
        if k in valid_mean_se and k not in out and _is_strict_finite_number_iip0024(v):
            if k in STAGE0_STRICT_INTEGER_COUNT_AND_DENOMINATOR_KEYS and not _is_integral_count_iip0026(v):
                continue
            lo, hi = valid_mean_se[k]
            fv = float(v)
            if lo <= fv <= hi:
                out[k] = fv
    return out


def evaluate_stage0_research_and_molecular_guardrails(
    payload: dict[str, Any],
    *,
    root: Path | str | None = None,
) -> dict[str, Any]:
    """
    Single authoritative registry-bound Stage-0 research & molecular guardrail evaluator (IIP-0026; DEC-105).
    Decoupled from Bayesian calibration admission (`validate_product_and_estimand_compatibility`),
    derives `lifecycle_state` and `execution_supported` strictly from the versioned product registry
    (rejecting caller override keys fail-closed), and produces identical branch (`switch_decisions`)
    and lot decisions across Direct evaluation, Serialized JSON Ingestion, and Chronological Replay.
    """
    if not isinstance(payload, dict):
        raise ValueError('Stage-0 research guardrail payload must be a dictionary (IIP-0026; DEC-105)')
    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    p_mod = _load_product_contract_module(root_path)
    prod_doc = p_mod.load_product_contracts(root_path)
    prod_by_id = {str(p['product_variant_id']): p for p in prod_doc['products']}
    reg_by_id = {
        pid_k: _derive_versioned_product_registry_entry_iip0026(p)
        for pid_k, p in prod_by_id.items()
    }

    _nested_meas_d = payload.get('measurements') if isinstance(payload.get('measurements'), dict) else {}
    _nested_scope_d = payload.get('scope') if isinstance(payload.get('scope'), dict) else {}
    _nested_mol_d = payload.get('iip0026_molecular_observations') if isinstance(payload.get('iip0026_molecular_observations'), dict) else {}
    _direct_containers = (payload, _nested_meas_d, _nested_scope_d, _nested_mol_d)
    _direct_nonstring_key = any(
        isinstance(_c, dict) and any(not isinstance(_rk, str) for _rk in _c)
        for _c in _direct_containers
    )
    _direct_malformed_container = any(
        name in payload and not isinstance(payload[name], dict)
        for name in ('scope', 'measurements', 'iip0026_molecular_observations')
    )
    _direct_scope_conflict = False
    for _sc_names in (('product_variant_id',), ('lot_id', 'current_lot_id'), ('campaign_id', 'assay_campaign_id')):
        _sc_vals = [
            _ctr[_k]
            for _ctr in (payload, _nested_scope_d)
            for _k in _sc_names
            if _k in _ctr
        ]
        if len(_sc_vals) >= 2 and not all(_v == _sc_vals[0] for _v in _sc_vals):
            if _sc_names == ('product_variant_id',) and all(isinstance(_v, str) for _v in _sc_vals):
                _canon_set = {_canonicalize_product_id_iip0031(_v, root_path)[0] or _v for _v in _sc_vals}
                if len(_canon_set) == 1:
                    continue
            _direct_scope_conflict = True

    pid_raw = payload.get('product_variant_id') if 'product_variant_id' in payload else _nested_scope_d.get('product_variant_id')
    if pid_raw is None and 'product_variant_id' not in payload:
        pid_raw = _nested_scope_d.get('product_variant_id')
    pid, pdef_resolved = _canonicalize_product_id_iip0031(pid_raw, root_path)
    if pid is None or pdef_resolved is None or pid not in reg_by_id:
        empty_s0 = evaluate_stage0_decision_tree({}, root=root_path)
        return {
            'overall_status': 'HOLD_UNREGISTERED_OR_MISSING_PRODUCT_VARIANT_ID',
            'lot_release_readiness_tier': 'HOLD_UNREGISTERED_OR_MISSING_PRODUCT_VARIANT_ID',
            'iip0026_reasons': ['HOLD_UNREGISTERED_OR_MISSING_PRODUCT_VARIANT_ID'],
            'hard_stop_precedence_winner': None,
            'lifecycle_state': 'UNREGISTERED',
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': empty_s0['switch_decisions'],
        }

    reg = reg_by_id[pid]
    pdef_cur = prod_by_id[pid]
    lifecycle_state = reg['lifecycle_state']
    registry_exec_supported = bool(reg['execution_supported'])

    _all_s9_rep_d = frozenset(k for ks in S9_REPLICATE_KEYS.values() for k in ks)
    _all_alias_d = frozenset(k for sh, ob, _ in _IIP0026_DESIGNB_TYK2_ALIAS_PAIRS for k in (sh, ob)) | frozenset({
        'post_edit_per_locus_segmented_baf', 'per_locus_segmented_baf',
        'obs_post_edit_segmented_baf_chr19q', 'obs_post_edit_segmented_baf_chr11p15', 'obs_post_edit_segmented_baf_chr6p21',
    })
    _excluded_merge_d = _all_s9_rep_d | _all_alias_d | frozenset({
        'scope', 'measurements', 'iip0026_molecular_observations',
        'tier6b_sorted_lot_observations', 'iip0031_locus2_record', 'iip0029_twoclock_record',
        '_s9_replicate_count_conflict', '_stage0_conflicting_measurement',
    })
    _d_s9_resolved, _d_s9_conflicts, _d_s9_any_conflict = _resolve_s9_replicate_counts_detailed_iip0031(
        _direct_containers,
        flag_explicit_none=True,
    )
    _clean_aliases_d, _alias_cands_d, _clean_baf_d, _baf_cands_d, _alias_or_baf_conflict_d = _reconcile_iip0026_aliases_and_baf_iip0031(
        _direct_containers
    )
    _direct_clean_obs, _direct_probes, _direct_meas_conflict = _merge_stage0_container_fields_iip0031(
        _direct_containers,
        excluded_keys=_excluded_merge_d,
        key_predicate=lambda k: k.startswith('obs_') or k in IIP0026_ALLOWED_MOLECULAR_OBS_KEYS or k in _ALL_S11_RECOGNIZED_EVIDENCE_KEYS or k in STAGE0_DECISION_TREE_BOOLEAN_KEYS,
    )
    canonical_payload: dict[str, Any] = {
        k: v for _src_c in (_nested_scope_d, payload) for k, v in _src_c.items() if isinstance(k, str)
    }
    canonical_payload['product_variant_id'] = pid
    for _ck_rm in list(canonical_payload.keys()):
        if (
            _ck_rm not in _excluded_merge_d
            and (_ck_rm.startswith('obs_') or _ck_rm in IIP0026_ALLOWED_MOLECULAR_OBS_KEYS or _ck_rm in _ALL_S11_RECOGNIZED_EVIDENCE_KEYS or _ck_rm in STAGE0_DECISION_TREE_BOOLEAN_KEYS)
        ):
            canonical_payload.pop(_ck_rm, None)
        elif _ck_rm in _all_alias_d or _ck_rm in _all_s9_rep_d:
            canonical_payload.pop(_ck_rm, None)
    canonical_payload.update(_direct_clean_obs)
    canonical_payload.update(_clean_aliases_d)
    if _clean_baf_d is not None:
        canonical_payload['post_edit_per_locus_segmented_baf'] = dict(_clean_baf_d)
    for _grp_d, _gkeys_d in S9_REPLICATE_KEYS.items():
        if _d_s9_resolved.get(_grp_d) is not None:
            canonical_payload[_gkeys_d[0]] = _d_s9_resolved[_grp_d]
    if _d_s9_any_conflict:
        canonical_payload['_s9_replicate_count_conflict'] = True

    safe_s0_obs = _build_safe_stage0_tree_obs_iip0026(canonical_payload, root_path=root_path)
    tree_eval, _tree_alias_conflict = _eval_stage0_tree_preserving_s9_iip0031(safe_s0_obs, root_path=root_path)
    if not tree_eval['safety_kill_s9_triggered'] and _direct_probes:
        for _p_cand in _direct_probes:
            _p_full = {**canonical_payload, **_p_cand}
            _p_safe = _build_safe_stage0_tree_obs_iip0026(_p_full, root_path=root_path)
            _p_tree, _p_ac = _eval_stage0_tree_preserving_s9_iip0031(_p_safe, root_path=root_path)
            if _p_ac:
                _tree_alias_conflict = True
            if _p_tree['safety_kill_s9_triggered']:
                tree_eval = _p_tree
                break
    switch_decisions = dict(tree_eval['switch_decisions'])

    _l2_bind_cur = pdef_cur.get('staged_locus2_modeb_donor_binding') if isinstance(pdef_cur, dict) else None
    _target_scope_6b = {
        'product_variant_id': pid,
        'construct_sha256': canonical_payload.get('construct_sha256'),
        'background_id': canonical_payload.get('background_id') or (_l2_bind_cur.get('background_id') if isinstance(_l2_bind_cur, dict) else None),
        'donor_seedline_id': canonical_payload.get('donor_seedline_id', canonical_payload.get('seed_line_id')),
        'mcb_id': canonical_payload.get('mcb_id'),
        'lot_id': canonical_payload.get('lot_id', canonical_payload.get('current_lot_id')),
        'campaign_id': canonical_payload.get('campaign_id', canonical_payload.get('assay_campaign_id')),
    }
    _direct_t6b_merged: list[Any] = []
    _direct_t6b_present = False
    for _d_src in (payload, _nested_scope_d, _nested_mol_d, _nested_meas_d):
        if 'tier6b_sorted_lot_observations' in _d_src:
            _direct_t6b_present = True
            _dv = _d_src['tier6b_sorted_lot_observations']
            if isinstance(_dv, (list, tuple)):
                _direct_t6b_merged.extend(_dv)
            else:
                _direct_t6b_merged.append('__EXPLICIT_NULL_OR_MALFORMED_TIER6B_CONTAINER__')
    tier6b_eval = evaluate_tier6b_per_fraction_thinning_iip0031(
        _direct_t6b_merged if _direct_t6b_present else [],
        target_scope=_target_scope_6b,
        root=root_path,
    )
    s11_eval = _evaluate_all_s11_containers_iip0031(
        canonical_payload,
        (payload, _nested_scope_d, _nested_mol_d, _nested_meas_d),
        root_path=root_path,
        product_def=pdef_cur,
    )
    switch_decisions[RULE_S11_GATE_ID_IIP0031] = str(s11_eval['status'])
    switch_decisions[f'JSON_{RULE_S11_GATE_ID_IIP0031}'] = str(s11_eval['status'])
    _t6b_common = {
        'tier6b_sorted_lot_evaluation': tier6b_eval,
        'tier6b_disposition': tier6b_eval['disposition'],
        'tier6b_reason': tier6b_eval['reason_code'],
        's11_evaluation': s11_eval,
        's9_point_estimate_fallback_se_policy': tree_eval.get(
            's9_point_estimate_fallback_se_policy',
            'SCALE_AWARE_DEFAULT_SE_IMPUTED_ON_SCALAR_OMISSION_TYPED_ATOM_REQUIRES_MEASURED_SE',
        ),
        's9_se_imputed_endpoints': list(tree_eval.get('s9_se_imputed_endpoints', [])),
    }

    _caller_override_rejected = any(
        isinstance(k, str) and (
            (k in IIP0026_FORBIDDEN_CALLER_OVERRIDE_KEYS)
            or k.endswith('_lifecycle_state')
            or k.endswith('_execution_supported')
        )
        for _ovr_ctr in (payload, _nested_scope_d, _nested_mol_d)
        for k in _ovr_ctr
    )

    # PRECEDENCE 1: Independent S9 Hard Safety Stop ALWAYS keeps first precedence:
    if tree_eval['safety_kill_s9_triggered']:
        _s9_reasons = ['HARD_SAFETY_STOP_S9']
        if _caller_override_rejected:
            _s9_reasons.append('FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED')
        if _direct_scope_conflict:
            _s9_reasons.append('FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING')
        return {
            'overall_status': 'HARD_SAFETY_STOP_REJECTED',
            'lot_release_readiness_tier': 'HARD_SAFETY_STOP_REJECTED',
            'iip0026_reasons': _s9_reasons,
            'hard_stop_precedence_winner': 'S9_HARD_SAFETY_STOP',
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    # PRECEDENCE 2: Attributable sorted-lot Tier 6B integer k_i >= 1 latches FAIL_HARD_STOP BEFORE
    # calibration holds, completeness holds, malformed records, caller overrides, or non-executable holds:
    if tier6b_eval.get('hard_stop') is True:
        _t6_reasons = ['FAIL_TIER6B_OUTGROWTH_DETECTED']
        if _caller_override_rejected:
            _t6_reasons.append('FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED')
        if _direct_scope_conflict:
            _t6_reasons.append('FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING')
        return {
            'overall_status': 'FAIL_HARD_STOP',
            'lot_release_readiness_tier': 'HARD_SAFETY_STOP_REJECTED',
            'iip0026_reasons': _t6_reasons,
            'hard_stop_precedence_winner': 'FAIL_TIER6B_OUTGROWTH_DETECTED',
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    if _direct_scope_conflict:
        return {
            'overall_status': 'FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING',
            'lot_release_readiness_tier': 'FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING',
            'iip0026_reasons': ['FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING'],
            'hard_stop_precedence_winner': None,
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    lot_id = canonical_payload.get('lot_id', canonical_payload.get('current_lot_id'))
    campaign_id = canonical_payload.get('campaign_id', canonical_payload.get('assay_campaign_id'))
    _missing_lot_or_campaign = bool(
        not isinstance(lot_id, str)
        or not lot_id.strip()
        or not isinstance(campaign_id, str)
        or not campaign_id.strip()
    )

    # PRECEDENCE 3: Rule S11 hard architectural/safety stops evaluate BEFORE caller overrides, missing-obs, or research holds:
    if not _missing_lot_or_campaign and str(s11_eval.get('status', '')).startswith('FAIL_'):
        _s11_fail = str(s11_eval['status'])
        _s11_reasons_out = list(s11_eval.get('reasons', [_s11_fail]))
        if _caller_override_rejected and 'FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED' not in _s11_reasons_out:
            _s11_reasons_out.append('FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED')
        if (
            payload.get('_stage0_conflicting_measurement') is True
            or _direct_meas_conflict
            or _alias_or_baf_conflict_d
            or _tree_alias_conflict
        ) and 'HOLD_CONFLICTING_STAGE0_MEASUREMENT' not in _s11_reasons_out:
            _s11_reasons_out.append('HOLD_CONFLICTING_STAGE0_MEASUREMENT')
        return {
            'overall_status': _s11_fail,
            'lot_release_readiness_tier': _s11_fail,
            'iip0026_reasons': _s11_reasons_out,
            'hard_stop_precedence_winner': _s11_fail,
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    # Reject any caller attempt to supply or override lifecycle_state or execution_supported:
    if _caller_override_rejected:
        return {
            'overall_status': 'FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED',
            'lot_release_readiness_tier': 'FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED',
            'iip0026_reasons': ['FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED'],
            'hard_stop_precedence_winner': None,
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    if _missing_lot_or_campaign:
        return {
            'overall_status': 'HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE',
            'lot_release_readiness_tier': 'HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE',
            'iip0026_reasons': ['HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE'],
            'hard_stop_precedence_winner': None,
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    if reg['bound_c1_c5_pair_sha256'] is not None:
        supplied_sha = canonical_payload.get('bound_c1_c5_pair_sha256')
        if supplied_sha != reg['bound_c1_c5_pair_sha256']:
            reason_sha = (
                'HOLD_MISSING_BOUND_C1_C5_DONOR_PAIR_DIGEST'
                if supplied_sha is None
                else 'FAIL_DONOR_REGISTRY_DIGEST_MISMATCH'
            )
            return {
                'overall_status': reason_sha,
                'lot_release_readiness_tier': reason_sha,
                'iip0026_reasons': [reason_sha],
                'hard_stop_precedence_winner': None,
                'lifecycle_state': lifecycle_state,
                'execution_supported': False,
                'effective_execution_permitted': False,
                'switch_decisions': switch_decisions,
                **_t6b_common,
            }

    # Evaluate candidate Design-B / TYK2 / BAF adverse failures across all containers BEFORE returning
    # conflict holds so a contradictory 0-alias can NEVER suppress an attributable positive failure (#13, #15):
    _cand_iip0026_fails = _evaluate_iip0026_candidate_fails_iip0031(
        reg,
        _direct_containers,
        _alias_cands_d,
        _baf_cands_d,
    )

    _, _s9_rep_conflict_groups, _s9_rep_has_conflict = _resolve_s9_replicate_counts_detailed_iip0031(
        _direct_containers,
        flag_explicit_none=True,
    )
    s9_replicate_count_conflict = bool(payload.get('_s9_replicate_count_conflict') is True or _s9_rep_has_conflict)
    stage0_conflicting_measurement = bool(
        payload.get('_stage0_conflicting_measurement') is True
        or _direct_meas_conflict
        or _direct_malformed_container
        or _direct_nonstring_key
        or _alias_or_baf_conflict_d
        or _tree_alias_conflict
    )
    if stage0_conflicting_measurement:
        _confl_reasons: list[str] = list(_cand_iip0026_fails)
        if 'HOLD_CONFLICTING_STAGE0_MEASUREMENT' not in _confl_reasons:
            _confl_reasons.append('HOLD_CONFLICTING_STAGE0_MEASUREMENT')
        if s9_replicate_count_conflict and 'HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT' not in _confl_reasons:
            _confl_reasons.append('HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT')
        _confl_status = _cand_iip0026_fails[0] if _cand_iip0026_fails else 'HOLD_CONFLICTING_STAGE0_MEASUREMENT'
        return {
            'overall_status': _confl_status,
            'lot_release_readiness_tier': _confl_status,
            'iip0026_reasons': _confl_reasons,
            'hard_stop_precedence_winner': None,
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }
    missing_obs = [k for k in IIP0026_REQUIRED_CORE_S9_OBS_KEYS if k not in safe_s0_obs]
    if (
        missing_obs
        or float(canonical_payload['obs_tp53_ppm1d_bcor_duplex_depth']) <= 0.0
        or float(canonical_payload['obs_castseq_duplex_or_amplicon_depth']) <= 0.0
        or not _is_integral_count_iip0026(canonical_payload['obs_mcb_clonogenic_survivor_count_k'])
        or not _is_integral_count_iip0026(canonical_payload['obs_mcb_effective_clonogenic_trials'])
        or int(canonical_payload['obs_mcb_effective_clonogenic_trials']) <= 0
    ):
        _missing_s9_status = (
            'HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT'
            if s9_replicate_count_conflict
            else 'HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS'
        )
        _missing_s9_reasons = (
            ['HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT', 'HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS']
            if s9_replicate_count_conflict
            else ['HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS']
        )
        return {
            'overall_status': _missing_s9_status,
            'lot_release_readiness_tier': _missing_s9_status,
            'iip0026_reasons': _missing_s9_reasons,
            'hard_stop_precedence_winner': None,
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
            **_t6b_common,
        }

    reasons: list[str] = list(tree_eval['s9_frequency_assay_holds'])
    for _, count_key in S9_FREQUENCY_ASSAY_KEYS.values():
        if any(isinstance(_c, dict) and count_key in _c for _c in _direct_containers) and count_key not in safe_s0_obs:
            reasons.append('HOLD_S9_INVALID_OR_CONFLICTING_EVENT_COUNT')
            break
    if s9_replicate_count_conflict:
        reasons.append('HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT')
    sera_val = canonical_payload.get('iip0026_schema_era', canonical_payload.get('schema_era'))
    if isinstance(sera_val, str) and sera_val.startswith('LEGACY_PRE_'):
        reasons.append('HOLD_PROVISIONAL_LEGACY_SCHEMA_MISSING_IIP0026_EVIDENCE')

    if canonical_payload.get('donor_seed_gnomad_target_site_wgs_verified') is not True:
        reasons.append('HOLD_MISSING_DONOR_TARGET_WGS')

    per_locus_baf = _clean_baf_d if _clean_baf_d is not None else _extract_per_locus_baf_dict_iip0026(canonical_payload)
    sv_recip = canonical_payload.get('post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified')
    req_loci = reg['required_baf_loci']
    if (
        not isinstance(per_locus_baf, dict)
        or any(loc not in per_locus_baf or not _is_strict_finite_number_iip0024(per_locus_baf[loc]) for loc in req_loci)
        or sv_recip is not True
    ):
        reasons.append('HOLD_MISSING_SEGMENTED_PHASED_BAF_SURVEILLANCE')
    else:
        if any(_baf_locus_discordant_iip0031(per_locus_baf[loc]) for loc in req_loci):
            reasons.append('FAIL_DISCORDANT_PER_LOCUS_SEGMENTED_BAF')

    if reg['requires_endogenous_hlac_and_locus5_gates']:
        cfd_max = canonical_payload.get('obs_design_b_hlac_phased_guide_cfd_max', canonical_payload.get('design_b_hlac_phased_guide_cfd_max'))
        sep_pass = canonical_payload.get('design_b_hla_ab_separated_pass_delivery_verified')
        chr6_loh = canonical_payload.get('design_b_chr6p21_1p41mb_loh_sv_cleared')
        if (
            not _is_strict_finite_number_iip0024(cfd_max)
            or not (Decimal('0') <= Decimal(str(cfd_max)) < Decimal('0.05'))
            or sep_pass is not True
            or chr6_loh is not True
        ):
            reasons.append('HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE')

        dt9_pos = canonical_payload.get('design_b_hlac_dt9_positive')
        hlae_3d12_pos = canonical_payload.get('design_b_hlae_3d12_positive')
        dt9_delta = canonical_payload.get(
            'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio',
            canonical_payload.get('design_b_hlac_dt9_matched_hlae_only_control_delta_ratio'),
        )
        hlac_orth = canonical_payload.get('design_b_hlac_allele_specific_surface_or_phased_mrna_verified')
        if dt9_pos is not True or hlae_3d12_pos is not True:
            if 'HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE' not in reasons:
                reasons.append('HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE')
        elif dt9_delta is None and hlac_orth is not True:
            reasons.append('HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY')
        elif dt9_delta is not None:
            if not _is_strict_finite_number_iip0024(dt9_delta) or Decimal(str(dt9_delta)) <= Decimal('0'):
                reasons.append('HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY')
            elif Decimal(str(dt9_delta)) < Decimal('2') and hlac_orth is not True:
                reasons.append('FAIL_HLAC_SURFACE_LOSS_MASKED_BY_DT9_HLAE_CROSSREACTIVITY')

        r5_bp = canonical_payload.get('obs_design_b_locus5_5prime_race_coord_bp', canonical_payload.get('design_b_locus5_5prime_race_coord_bp'))
        r3_bp = canonical_payload.get('obs_design_b_locus5_3prime_race_coord_bp', canonical_payload.get('design_b_locus5_3prime_race_coord_bp'))
        max_donor_bp = 4078 if reg['bound_c1_c5_pair_sha256'] == DESIGN_B_C1_C5_PAIR_SHA256 else 10000
        if r5_bp is None or r3_bp is None:
            reasons.append('HOLD_MISSING_DESIGNB_LOCUS5_RACE_COORDINATES')
        elif not (
            _is_strict_int_iip0024(r5_bp)
            and _is_strict_int_iip0024(r3_bp)
            and 0 < r5_bp < r3_bp <= max_donor_bp
        ):
            reasons.append('FAIL_INVALID_RACE_COORDINATES')

        pdl1_r5 = canonical_payload.get('design_b_locus5_pdl1_5prime_race_coord_bp')
        if reg['bound_c1_c5_pair_sha256'] == DESIGN_B_C1_C5_PAIR_SHA256 and (pdl1_r5 is None or pdl1_r5 == 'UNRESOLVED'):
            reasons.append('HOLD_PENDING_5PRIME_RACE_TSS')
        elif pdl1_r5 == 'UNRESOLVED':
            reasons.append('HOLD_PENDING_5PRIME_RACE_TSS')
        elif pdl1_r5 is not None:
            if not _is_strict_int_iip0024(pdl1_r5) or pdl1_r5 <= 0 or pdl1_r5 > 2381:
                reasons.append('FAIL_INVALID_RACE_COORDINATES')
            elif pdl1_r5 <= 2308:
                reasons.append('FAIL_MATURE_TRANSCRIPT_INCLUDES_UPSTREAM_OUT_OF_FRAME_KOZAK_ORF')

        b_umi_n = canonical_payload.get('obs_design_b_locus5_race_umi_n_eff', canonical_payload.get('design_b_locus5_race_umi_n_eff'))
        b_umi_k = canonical_payload.get('obs_design_b_locus5_aberrant_umi_k', canonical_payload.get('design_b_locus5_aberrant_umi_k'))
        if not _is_strict_int_iip0024(b_umi_n) or b_umi_n < 3689 or not _is_strict_int_iip0024(b_umi_k) or b_umi_k < 0:
            reasons.append('HOLD_MISSING_DESIGNB_LOCUS5_TRANSCRIPT_UMI_3689')
        elif b_umi_k > 0:
            reasons.append('FAIL_DESIGNB_LOCUS5_CHIMERIC_TRANSCRIPT_DETECTED')
    elif reg['requires_b_sct_gates']:
        if canonical_payload.get('b_sct_ciita_exon5_and_trimer_verified') is not True:
            reasons.append('HOLD_MISSING_BSCT_CIITA_EXON5_AND_TRIMER_EVIDENCE')

    if reg['requires_tyk2_pe_gates']:
        pe_mod = canonical_payload.get('tyk2_pe_modality')
        if not isinstance(pe_mod, str) or pe_mod not in ALLOWED_TYK2_PE_MODALITIES or canonical_payload.get('tyk2_phased_haplotype_sap_verified') is not True:
            reasons.append('HOLD_MISSING_TYK2_PE_MODALITY_OR_PHASED_HAPLOTYPE')
        elif pe_mod == 'PE4MAX_MSI_GATED':
            msi_k = canonical_payload.get('obs_tyk2_pe4max_pentaplex_msi_shifts_k', canonical_payload.get('tyk2_pe4max_pentaplex_msi_shifts_k'))
            if msi_k is None or not _is_strict_int_iip0024(msi_k) or msi_k < 0:
                reasons.append('HOLD_MISSING_MSI')
            elif msi_k > 0:
                reasons.append('FAIL_MSI_DETECTED')
            pe_st = evaluate_pe4max_piecewise_repeat_indel_cap_iip0026(
                canonical_payload.get('obs_tyk2_pe4max_wgs_repeat_indels_clone', canonical_payload.get('tyk2_pe4max_wgs_repeat_indels_clone')),
                canonical_payload.get('obs_tyk2_pe4max_wgs_repeat_indels_ctrl', canonical_payload.get('tyk2_pe4max_wgs_repeat_indels_ctrl')),
                canonical_payload.get('obs_tyk2_pe4max_coding_or_splice_de_novo_indels', canonical_payload.get('tyk2_pe4max_coding_or_splice_de_novo_indels')),
            )
            if pe_st != 'PASS':
                if not (msi_k is None and pe_st == 'HOLD_MISSING_PE4MAX_REPEAT_INDEL_EVIDENCE'):
                    reasons.append(pe_st)

    for _cf in _cand_iip0026_fails:
        if _cf not in reasons:
            reasons.append(_cf)

    if reg['is_held_design_y']:
        reasons.append('HOLD_DESIGNY_EXPLORATORY_UNBIASED_CAPTURE_AND_63K_OR_90K_LOD_REQUIRED')

    # Propagate IIP-0029 two-clock/two-pulse & >138.3 IU Gate 0S escalation holds identically on all variants:
    _iip0029_tree_eval = tree_eval.get('iip0029_twoclock_maturation_and_repeated_challenge_evaluation')
    if isinstance(_iip0029_tree_eval, dict) and _iip0029_tree_eval.get('gate_pass') is not True:
        _r29 = str(_iip0029_tree_eval.get('status', 'HOLD_IIP0029_SUBGATE_FAILED'))
        if _r29 not in reasons:
            reasons.append(_r29)

    # Propagate IIP-0031 Rule S11 holds and mandatory sorted-lot Tier 6B per-fraction thinning qualification:
    if s11_eval.get('gate_pass') is not True:
        _r31_s11 = str(s11_eval.get('status', 'HOLD_S11_LOCUS2_GATE_FAILED'))
        if _r31_s11 not in reasons:
            reasons.append(_r31_s11)
    if tier6b_eval.get('gate_pass') is not True:
        _r31_6b = str(tier6b_eval.get('status', 'HOLD_TIER6B_UNQUALIFIED_FRACTION'))
        if _r31_6b not in reasons:
            reasons.append(_r31_6b)

    # Preserve stricter applicable Stage-0 tree holds on non-research-comparator executable variants:
    tree_tier = str(tree_eval['lot_release_readiness_tier'])
    if (
        not reasons
        and lifecycle_state != 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
        and registry_exec_supported
        and tree_tier.startswith(('HOLD_', 'NO_GO_', 'HARD_SAFETY_STOP'))
    ):
        reasons.append(str(tree_eval['resolved_architecture_assignments'].get('lot_maturation_status', tree_tier)))

    fails = [r for r in reasons if r.startswith('FAIL_')]
    if fails:
        status = fails[0]
    elif reasons:
        status = reasons[0]
    else:
        status = (
            'RESEARCH_COMPARATOR_STAGE0_PANEL_CLEARED_PENDING_STATE4_PROMOTION'
            if lifecycle_state == 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
            else 'STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED'
        )
    effective_exec = bool(
        registry_exec_supported
        and len(reasons) == 0
        and lifecycle_state != 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
    )
    return {
        'overall_status': status,
        'lot_release_readiness_tier': tree_tier if status.startswith('STAGE0_') else status,
        'iip0026_reasons': reasons,
        'hard_stop_precedence_winner': None,
        'lifecycle_state': lifecycle_state,
        'execution_supported': registry_exec_supported if len(reasons) == 0 else False,
        'effective_execution_permitted': effective_exec,
        'endogenous_hlac_gate_applied': reg['requires_endogenous_hlac_and_locus5_gates'],
        'b_sct_route_separated': reg['requires_b_sct_gates'],
        'switch_decisions': switch_decisions,
        **_t6b_common,
        **(
            {'iip0029_twoclock_maturation_and_repeated_challenge_evaluation': _iip0029_tree_eval}
            if isinstance(_iip0029_tree_eval, dict)
            else {}
        ),
    }


# ============================================================================
# Post-Depletion Gate 0d, Paired Intact-Plus-Increment CD16+ ADCC, Upstream B2
# Research-Comparator Authorization, Measured-Potency DEC-084 Mixture, and
# Visit/Lot/Sister-Bound 72h/d28/d84 Dual-Modality IAPP Evaluators
# (IIP-0028; DEC-106)
# ============================================================================
ITM2B_NM021999_5_CDS_SHA256 = '6fc8d3d7a9b4851a9b50c026d44f02e7d2539e122a0b727d2a5cc4ee93df8903'
ITM2B_NM021999_5_FASTA_SHA256 = '4a4cb15d0544884a62d86ebf4137c5284edad8210b6e32e74c2eb0aff005d13e'


def _valid_domain_scalar_iip0028(val: Any, lo: float, hi: float) -> bool:
    return _is_strict_finite_number_iip0024(val) and (lo <= float(val) <= hi)


def _valid_positive_se_iip0028(se: Any) -> bool:
    return _is_strict_finite_number_iip0024(se) and (float(se) > 0.0)


def _valid_strict_unpadded_id_iip0028(val: Any, min_len: int = 4) -> bool:
    return isinstance(val, str) and (val == val.strip()) and (len(val) >= min_len)


def build_valid_upstream_b2_record_iip_0028(
    mcb_id: str = 'MCB-B2-001',
    lot_id: str = 'LOT-B2-001',
) -> dict[str, Any]:
    """Build canonical upstream ISLET1_B2_5LOCUS_BASE research-comparator authorization record (IIP-0028; DEC-106)."""
    return {
        'product_variant_id': 'ISLET1_B2_5LOCUS_BASE',
        'mcb_id': mcb_id,
        'lot_id': lot_id,
        'lifecycle_state': 'STAGE0_RESEARCH_COMPARATOR_AUTHORIZED_FOR_PAIRED_STUDY',
        'execution_supported': False,
        'design_b_fasta_sha256': DESIGN_B_FASTA_SHA256,
        'bound_c1_c5_pair_sha256': DESIGN_B_C1_C5_PAIR_SHA256,
        'verified_qc_atoms': {
            'mcb_hlac_phased_allele_and_ab_ko_qc': {
                'atom_id': f'QCATOM-MCB-HLAC-{mcb_id}',
                'mcb_id': mcb_id,
                'status': 'PASS',
                'biallelic_hla_ab_ciita_ko_pct': 99.4,
                'phased_hlac_allele': 'HLA-C*03:04:01',
            },
            'mcb_c1_c5_junction_and_5p3p_race_qc': {
                'atom_id': f'QCATOM-MCB-C1C5-{mcb_id}',
                'mcb_id': mcb_id,
                'status': 'PASS',
                'c1_sha256_prefix': '3392af83',
                'c5_sha256_prefix': 'ba264b14',
                'five_prime_three_prime_race_exact': True,
                'biallelic_c1_c5_hdr_pct': 99.1,
            },
            'lot_hlac_dt9_and_hlae_3d12_specificity_qc': {
                'atom_id': f'QCATOM-LOT-DT9-3D12-{lot_id}',
                'lot_id': lot_id,
                'status': 'PASS',
                'dt9_hlac_fold_over_hlae_only_ctrl': 8.5,
                'hlae_3d12_molecules': 28500.0,
            },
        },
    }


def verify_upstream_b2_paired_study_authorization_iip_0028(
    auth_rec: Any,
    *,
    expected_mcb_id: str,
    expected_lot_id: str,
) -> tuple[bool, str]:
    """
    Strict upstream registry and QC-atom verifier for ISLET1_B2_5LOCUS_BASE (IIP-0028; DEC-106):
      - Never trusts a caller boolean;
      - Binds DESIGN_B_FASTA_SHA256 (4bb20389...) and DESIGN_B_C1_C5_PAIR_SHA256 (82c351af...);
      - Enforces execution_supported is False while authorizing paired Stage-0 ex-vivo responder
        comparisons only when all three MCB/lot allele, c1+c5 RACE/junction, and DT9+3D12
        HLA-E-controlled surface specificity QC atoms pass.
    """
    if not isinstance(auth_rec, dict):
        return (False, 'HOLD_B2_MISSING_UPSTREAM_AUTHORIZATION_RECORD')
    if auth_rec.get('product_variant_id') != 'ISLET1_B2_5LOCUS_BASE':
        return (False, 'HOLD_B2_UPSTREAM_PRODUCT_MISMATCH')
    if auth_rec.get('mcb_id') != expected_mcb_id or auth_rec.get('lot_id') != expected_lot_id:
        return (False, 'HOLD_B2_UPSTREAM_MCB_OR_LOT_MISMATCH')
    if auth_rec.get('lifecycle_state') != 'STAGE0_RESEARCH_COMPARATOR_AUTHORIZED_FOR_PAIRED_STUDY':
        return (False, 'HOLD_B2_LIFECYCLE_NOT_AUTHORIZED_FOR_PAIRED_STUDY')
    if auth_rec.get('execution_supported') is not False:
        return (False, 'HOLD_B2_EXECUTION_SUPPORTED_MUST_REMAIN_FALSE')
    if auth_rec.get('design_b_fasta_sha256') != DESIGN_B_FASTA_SHA256:
        return (False, 'HOLD_B2_FASTA_DIGEST_MISMATCH')
    if auth_rec.get('bound_c1_c5_pair_sha256') != DESIGN_B_C1_C5_PAIR_SHA256:
        return (False, 'HOLD_B2_C1_C5_PAIR_DIGEST_MISMATCH')

    qc_atoms = auth_rec.get('verified_qc_atoms')
    if not isinstance(qc_atoms, dict):
        return (False, 'HOLD_B2_MISSING_UPSTREAM_QC_ATOMS')

    a1 = qc_atoms.get('mcb_hlac_phased_allele_and_ab_ko_qc')
    if (
        not isinstance(a1, dict)
        or not _valid_strict_unpadded_id_iip0028(a1.get('atom_id'))
        or a1.get('mcb_id') != expected_mcb_id
        or a1.get('status') != 'PASS'
        or not _valid_domain_scalar_iip0028(a1.get('biallelic_hla_ab_ciita_ko_pct'), 98.0, 100.0)
        or not _valid_strict_unpadded_id_iip0028(a1.get('phased_hlac_allele'), min_len=6)
    ):
        return (False, 'HOLD_B2_UNVERIFIED_MCB_HLAC_ALLELE_AND_AB_KO_QC_ATOM')

    a2 = qc_atoms.get('mcb_c1_c5_junction_and_5p3p_race_qc')
    if (
        not isinstance(a2, dict)
        or not _valid_strict_unpadded_id_iip0028(a2.get('atom_id'))
        or a2.get('mcb_id') != expected_mcb_id
        or a2.get('status') != 'PASS'
        or a2.get('c1_sha256_prefix') != '3392af83'
        or a2.get('c5_sha256_prefix') != 'ba264b14'
        or a2.get('five_prime_three_prime_race_exact') is not True
        or not _valid_domain_scalar_iip0028(a2.get('biallelic_c1_c5_hdr_pct'), 98.0, 100.0)
    ):
        return (False, 'HOLD_B2_UNVERIFIED_MCB_C1_C5_JUNCTION_AND_RACE_QC_ATOM')

    a3 = qc_atoms.get('lot_hlac_dt9_and_hlae_3d12_specificity_qc')
    if (
        not isinstance(a3, dict)
        or not _valid_strict_unpadded_id_iip0028(a3.get('atom_id'))
        or a3.get('lot_id') != expected_lot_id
        or a3.get('status') != 'PASS'
        or not _valid_domain_scalar_iip0028(a3.get('dt9_hlac_fold_over_hlae_only_ctrl'), 2.0, 1000.0)
        or not _valid_domain_scalar_iip0028(a3.get('hlae_3d12_molecules'), 15000.0, 1000000.0)
    ):
        return (False, 'HOLD_B2_UNVERIFIED_LOT_HLAC_DT9_HLAE_3D12_SPECIFICITY_QC_ATOM')

    return (True, 'AUTHORIZED_B2_STAGE0_PAIRED_STUDY_EXECUTION_SUPPORTED_FALSE')


def _read_se_alias_iip0028(obs: dict[str, Any], primary_se_key: str, *alias_se_keys: str) -> Any:
    if primary_se_key in obs and obs[primary_se_key] is not None:
        return obs[primary_se_key]
    for ak in alias_se_keys:
        if ak in obs and obs[ak] is not None:
            return obs[ak]
    return None


def evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
    obs: dict[str, Any],
    *,
    root: Path | str | None = None,
    upstream_b2_authorization_record: dict[str, Any] | None = None,
) -> tuple[bool, str]:
    """
    Fail-closed registry-bound Gate 0d + intact-and-increment CD16+ ADCC evaluator (IIP-0028; DEC-106):
      - Rejects omitted, None, empty, whitespace-padded, unregistered, caller-overridden, or Track-2 B'-SCT IDs;
      - Gates ISLET1_B2_5LOCUS_BASE on verified upstream 82c351af... c1+c5 authorization + QC atoms (execution_supported=False);
      - Separates Arm 1 reference, Arm 2 post-depletion (preserving all 6 Arm-2 explicit-SE endpoints + contemporaneous
        Arm-3 active-FasL control + unpooled washout/aging + segregated RETAINED_ACTIVE SA-FHL-1), and Arm 3 profiles;
      - Caps BOTH intact-serum total cellular killing AND paired IgG-attributable CD16+ ADCC increment at U95_t <= 5.0%
        (preventing high Fc-independent background injury such as 30% intact vs 29% depleted from clearing on a 1% delta),
        alongside positive control L95_t >= 25.0% and isotype background U95_t <= 5.0%.
    """
    if not isinstance(obs, dict):
        return (False, 'HOLD_INVALID_OBSERVATION_PAYLOAD')
    if any(k in obs for k in ('product_defs', 'override_product_defs', 'execution_supported', 'override_execution_supported')):
        return (False, 'HOLD_CALLER_PRODUCT_DEF_OR_EXECUTION_OVERRIDE_REJECTED')
    if 'product_variant_id' not in obs or obs['product_variant_id'] is None:
        return (False, 'HOLD_MISSING_PRODUCT_VARIANT_ID')
    prod_id = obs['product_variant_id']
    if not isinstance(prod_id, str) or prod_id == '':
        return (False, 'HOLD_MISSING_PRODUCT_VARIANT_ID')
    if prod_id != prod_id.strip():
        return (False, 'HOLD_WHITESPACE_PADDED_PRODUCT_VARIANT_ID')
    if prod_id.startswith('ISLET1_B_SCT_'):
        return (False, 'HOLD_UNSUPPORTED_TRACK2_BSCT_REQUIRES_SEPARATE_GATE0B_SCT')

    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    p_mod = _load_product_contract_module(root_path)
    prod_doc = p_mod.load_product_contracts(root_path)
    reg_by_id = {str(p['product_variant_id']): p for p in prod_doc['products']}
    reg_entry = reg_by_id.get(prod_id)
    if reg_entry is None:
        return (False, 'HOLD_UNREGISTERED_PRODUCT_VARIANT_ID')

    mcb_id = obs.get('mcb_id', obs.get('mcb_clone_id'))
    lot_id = obs.get('lot_id', obs.get('current_lot_id'))
    if not _valid_strict_unpadded_id_iip0028(mcb_id) or not _valid_strict_unpadded_id_iip0028(lot_id):
        return (False, 'HOLD_GATE0D_MISSING_OR_PADDED_MCB_OR_LOT_ID')

    if prod_id == 'ISLET1_B2_5LOCUS_BASE':
        if upstream_b2_authorization_record is None:
            if (
                obs.get('registered_and_qc_cleared') is True
                or obs.get('comparator_registered_and_qc_cleared') is True
                or any(str(k).endswith('_registered_and_qc_cleared') and obs.get(k) is True for k in obs)
            ):
                return (False, 'HOLD_B2_CALLER_BOOLEAN_REJECTED_MISSING_UPSTREAM_REGISTRY_RECORD')
            return (False, 'HOLD_B2_PENDING_UPSTREAM_COMPARATOR_REGISTRATION')
        b2_ok, b2_code = verify_upstream_b2_paired_study_authorization_iip_0028(
            upstream_b2_authorization_record,
            expected_mcb_id=str(mcb_id),
            expected_lot_id=str(lot_id),
        )
        if not b2_ok:
            return (False, b2_code)
    elif reg_entry.get('execution_supported') is False:
        return (False, 'HOLD_NON_EXECUTABLE_RESEARCH_COMPARATOR_VARIANT')

    arm = obs.get('adjunct_conditioning_history_arm')
    if arm not in ('ARM1_NEVER_CONDITIONED', 'ARM2_CONDITIONED_DEPLETED', 'ARM3_CONTINUOUSLY_ACTIVE'):
        return (False, 'HOLD_GATE0D_MISSING_EXPLICIT_HISTORY_ARM')

    n_rep = obs.get('n_independent_blocks', obs.get('obs_stage0_independent_lot_replicates_n'))
    if not _is_strict_int_iip0024(n_rep) or n_rep < 6:
        return (False, 'HOLD_GATE0D_INSUFFICIENT_INDEPENDENT_BLOCKS_LT_6')
    t_crit = student_t_ppf(0.95, n_rep - 1)

    if arm == 'ARM1_NEVER_CONDITIONED':
        for mk, sk1, sk2, lo, hi in (
            ('obs_arm1_unconditioned_cd4_si', 'obs_arm1_unconditioned_cd4_si_se', 'se_arm1_cd4_si', 0.1, 100.0),
            ('obs_arm1_unconditioned_annexin_v_pct', 'obs_arm1_unconditioned_annexin_v_se', 'se_arm1_annexin_v', 0.0, 100.0),
        ):
            s_v = _read_se_alias_iip0028(obs, sk1, sk2)
            if not _valid_domain_scalar_iip0028(obs.get(mk), lo, hi) or not _valid_positive_se_iip0028(s_v):
                return (False, f'HOLD_GATE0D_ARM1_INVALID_OR_MISSING_CONTROL_{mk}')
        return (False, 'HOLD_GATE0D_ARM1_NEVER_CONDITIONED_REFERENCE_CONTROL')

    if not _valid_strict_unpadded_id_iip0028(obs.get('recipient_serum_pair_id')):
        return (False, 'HOLD_ADCC_MISSING_OR_INVALID_RECIPIENT_SERUM_PAIR_ID')
    intact_sid = obs.get('intact_serum_sample_id')
    dep_sid = obs.get('igg_depleted_serum_sample_id')
    if not _valid_strict_unpadded_id_iip0028(intact_sid) or not _valid_strict_unpadded_id_iip0028(dep_sid) or intact_sid == dep_sid:
        return (False, 'HOLD_ADCC_MISSING_OR_DUPLICATE_PAIRED_SERUM_SAMPLE_IDS')
    pos_cid = obs.get('positive_adcc_control_id')
    iso_cid = obs.get('isotype_adcc_control_id')
    if not _valid_strict_unpadded_id_iip0028(pos_cid) or not _valid_strict_unpadded_id_iip0028(iso_cid) or pos_cid == iso_cid:
        return (False, 'HOLD_ADCC_MISSING_OR_DUPLICATE_POSITIVE_ISOTYPE_CONTROL_IDS')
    if obs.get('igg_depleted_matched_serum_control_verified', obs.get('cd16_adcc_igg_depletion_paired_control_verified')) is not True:
        return (False, 'HOLD_ADCC_MISSING_IGG_DEPLETED_MATCHED_SERUM_CONTROL')
    if not _valid_domain_scalar_iip0028(obs.get('obs_igg_depletion_efficiency_pct'), 95.0, 100.0):
        return (False, 'HOLD_ADCC_IGG_DEPLETION_EFFICIENCY_BELOW_95PCT')
    if obs.get('fcr_effector_fcgr3a_158v_typed_verified') is not True:
        return (False, 'HOLD_ADCC_UNVERIFIED_CD16_FCGR3A_EFFECTOR_COMPETENCE')

    # Enforce BOTH intact-serum total cellular killing U95_t <= 5.0% AND paired IgG-attributable increment U95_t <= 5.0%:
    intact_m = obs.get('obs_fcr_cd16_intact_serum_total_lysis_pct')
    intact_s = _read_se_alias_iip0028(obs, 'obs_fcr_cd16_intact_serum_total_lysis_se', 'se_intact_adcc')
    dep_m = obs.get('obs_fcr_cd16_igg_depleted_serum_lysis_pct')
    if not _valid_domain_scalar_iip0028(intact_m, 0.0, 100.0) or not _valid_positive_se_iip0028(intact_s):
        return (False, 'HOLD_MISSING_OR_INVALID_obs_fcr_cd16_intact_serum_total_lysis_pct')
    if dep_m is not None:
        if not _valid_domain_scalar_iip0028(dep_m, 0.0, 100.0):
            return (False, 'HOLD_MISSING_OR_INVALID_obs_fcr_cd16_igg_depleted_serum_lysis_pct')
        net_check = obs.get('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct')
        if _is_strict_finite_number_iip0024(net_check) and abs((float(intact_m) - float(dep_m)) - float(net_check)) > 1e-4:
            return (False, 'HOLD_ADCC_INTACT_MINUS_DEPLETED_ARITHMETIC_MISMATCH')
    if float(intact_m) + t_crit * float(intact_s) > 5.0:
        return (False, 'HOLD_EXCEEDED_UPPER95_obs_fcr_cd16_intact_serum_total_lysis_pct')

    pos_adcc_m = obs.get('obs_adcc_positive_control_lysis_pct', obs.get('obs_cd16_adcc_positive_control_lysis_pct'))
    pos_adcc_s = _read_se_alias_iip0028(obs, 'obs_cd16_adcc_positive_control_lysis_se', 'se_pos_adcc')
    if not _valid_domain_scalar_iip0028(pos_adcc_m, 0.0, 100.0) or not _valid_positive_se_iip0028(pos_adcc_s):
        return (False, 'HOLD_MISSING_OR_INVALID_obs_adcc_positive_control_lysis_pct')
    if float(pos_adcc_m) - t_crit * float(pos_adcc_s) < 25.0:
        return (False, 'HOLD_BELOW_LOWER95_obs_adcc_positive_control_lysis_pct')

    if arm == 'ARM3_CONTINUOUSLY_ACTIVE':
        if obs.get('obs_arm3_ongoing_adjunct_bioactivity_verified') is not True:
            return (False, 'HOLD_GATE0D_ARM3_ONGOING_BIOACTIVITY_UNVERIFIED')
        arm3_bounds = [
            ('obs_gate0d_activated_cd4_fasl_apoptosis_pct', ('obs_gate0d_activated_cd4_fasl_apoptosis_se', 'se_arm3_cd4_apoptosis'), 'lower', 80.0, 0.0, 100.0),
            ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index', ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index_se', 'obs_indirect_apc_dq2_dq8_cd4_stimulation_se', 'se_cd4_si'), 'upper', 1.5, 0.1, 100.0),
            ('obs_indirect_apc_cd4_ifng_si', ('obs_indirect_apc_cd4_ifng_si_se', 'obs_indirect_apc_cd4_ifng_se', 'se_ifng_si'), 'upper', 2.0, 0.1, 100.0),
            ('obs_gate0d_transwell_beta_annexin_v_pct', ('obs_gate0d_transwell_beta_annexin_v_se', 'se_annexin_v'), 'upper', 15.0, 0.0, 100.0),
            ('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct', ('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_se', 'se_net_adcc'), 'upper', 5.0, -20.0, 100.0),
            ('obs_adcc_isotype_background_lysis_pct', ('obs_adcc_isotype_background_lysis_se', 'se_isotype_adcc'), 'upper', 5.0, 0.0, 100.0),
        ]
        for mk, sks, direction, thresh, lo, hi in arm3_bounds:
            m = obs.get(mk)
            s = _read_se_alias_iip0028(obs, sks[0], *sks[1:])
            if not _valid_domain_scalar_iip0028(m, lo, hi) or not _valid_positive_se_iip0028(s):
                return (False, f'HOLD_ARM3_INVALID_OR_MISSING_{mk}')
            mf, sf = float(m), float(s)
            if direction == 'upper' and (mf + t_crit * sf) > thresh:
                return (False, f'HOLD_ARM3_EXCEEDED_UPPER95_{mk}')
            if direction == 'lower' and (mf - t_crit * sf) < thresh:
                return (False, f'HOLD_ARM3_BELOW_LOWER95_{mk}')
        if prod_id == 'ISLET1_B2_5LOCUS_BASE':
            return (True, 'RESEARCH_ONLY_PAIRED_STUDY_ARM3_CLEARED_B2_EXECUTION_SUPPORTED_FALSE')
        return (True, 'CLEARED_ARM3_ACTIVE_ADJUNCT_AND_FCR_ADCC_PROFILE')

    dep_mode = obs.get('adjunct_depletion_mode')
    if dep_mode not in ('PHYSICAL_WASHOUT', 'NATURAL_AGING_DECAY'):
        return (False, 'HOLD_GATE0D_POOLED_OR_UNSPECIFIED_DEPLETION_MODE')
    if obs.get('washout_and_natural_decay_unpooled_verified') is not True:
        return (False, 'HOLD_GATE0D_POOLED_WASHOUT_AND_AGING_ALIQUOTS')
    if obs.get('carrier_and_soluble_adjunct_below_lod_verified') is not True:
        return (False, 'HOLD_GATE0D_RESIDUAL_ADJUNCT_ABOVE_LOD')
    if obs.get('obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified') is not True:
        return (False, 'HOLD_GATE0D_FASL_DEPLETION_UNVERIFIED')
    if obs.get('obs_arm2_sa_thbd_bioactivity_depleted_verified') is not True:
        return (False, 'HOLD_GATE0D_SA_THBD_DEPLETION_UNVERIFIED')
    if obs.get('obs_arm2_contemporaneous_sham_handling_matched') is not True:
        return (False, 'HOLD_GATE0D_CONTEMPORANEOUS_SHAM_UNMATCHED')
    if obs.get('obs_arm3_ongoing_adjunct_bioactivity_verified') is not True:
        return (False, 'HOLD_GATE0D_PARALLEL_ARM3_BIOACTIVITY_CONTROL_UNVERIFIED')
    arm2_sid = obs.get('arm2_depleted_sample_id')
    arm3_ctrl_sid = obs.get('contemporaneous_arm3_control_sample_id')
    if not _valid_strict_unpadded_id_iip0028(arm2_sid) or not _valid_strict_unpadded_id_iip0028(arm3_ctrl_sid) or arm2_sid == arm3_ctrl_sid:
        return (False, 'HOLD_GATE0D_MISSING_OR_DUPLICATE_CONTEMPORANEOUS_ARM3_CONTROL_SAMPLE_ID')
    resp_hist = obs.get('adjunct_immune_responder_history')
    if resp_hist not in ('RETAINED_CONDITIONED_RESPONDERS', 'FRESH_RECIPIENT_RESPONDERS', 'PAIRED_BOTH'):
        return (False, 'HOLD_GATE0D_MISSING_IMMUNE_RESPONDER_HISTORY')
    fhl1_status = obs.get('adjunct_sa_fhl1_status')
    if fhl1_status == 'RETAINED_ACTIVE':
        return (False, 'HOLD_GATE0D_PARTIALLY_DEPLETED_SA_FHL1_RETAINED_ACTIVE_NOT_COMPLETE_DEPLETION')
    if fhl1_status not in ('DEPLETED_VERIFIED', 'ABSENT_FROM_LOT'):
        return (False, 'HOLD_GATE0D_INVALID_OR_MISSING_SA_FHL1_STATUS')

    arm2_bounds = [
        ('obs_parallel_arm3_control_cd4_apoptosis_pct', ('obs_parallel_arm3_control_cd4_apoptosis_se', 'se_parallel_arm3_cd4_apoptosis'), 'lower', 80.0, 0.0, 100.0),
        ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index', ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index_se', 'obs_indirect_apc_dq2_dq8_cd4_stimulation_se', 'se_cd4_si'), 'upper', 1.5, 0.1, 100.0),
        ('obs_indirect_apc_cd4_ifng_si', ('obs_indirect_apc_cd4_ifng_si_se', 'obs_indirect_apc_cd4_ifng_se', 'se_ifng_si'), 'upper', 2.0, 0.1, 100.0),
        ('obs_gate0d_transwell_beta_annexin_v_pct', ('obs_gate0d_transwell_beta_annexin_v_se', 'se_annexin_v'), 'upper', 15.0, 0.0, 100.0),
        ('obs_arm2_post_depletion_nhs_c5b9_mac_pct', ('obs_arm2_post_depletion_nhs_c5b9_mac_se', 'se_c5b9'), 'upper', 5.0, 0.0, 100.0),
        ('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct', ('obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_se', 'se_net_adcc'), 'upper', 5.0, -20.0, 100.0),
        ('obs_adcc_isotype_background_lysis_pct', ('obs_adcc_isotype_background_lysis_se', 'se_isotype_adcc'), 'upper', 5.0, 0.0, 100.0),
        ('obs_arm2_post_depletion_viable_beta_recovery_pct', ('obs_arm2_post_depletion_viable_beta_recovery_se', 'se_beta_rec'), 'lower', 85.0, 0.0, 150.0),
        ('obs_arm2_post_depletion_perifusion_gsis_si', ('obs_arm2_post_depletion_perifusion_gsis_si_se', 'se_gsis_si'), 'lower', 2.5, 0.1, 100.0),
        ('obs_post_depletion_input_normalized_cpep_ratio', ('obs_post_depletion_input_normalized_cpep_ratio_se', 'se_cpep_ratio'), 'lower', 0.85, 0.01, 5.0),
    ]
    for mk, sks, direction, thresh, lo, hi in arm2_bounds:
        m = obs.get(mk)
        s = _read_se_alias_iip0028(obs, sks[0], *sks[1:])
        if not _valid_domain_scalar_iip0028(m, lo, hi) or not _valid_positive_se_iip0028(s):
            return (False, f'HOLD_MISSING_OR_INVALID_{mk}')
        mf, sf = float(m), float(s)
        if direction == 'upper' and (mf + t_crit * sf) > thresh:
            return (False, f'HOLD_EXCEEDED_UPPER95_{mk}')
        if direction == 'lower' and (mf - t_crit * sf) < thresh:
            return (False, f'HOLD_BELOW_LOWER95_{mk}')
    if prod_id == 'ISLET1_B2_5LOCUS_BASE':
        return (True, 'RESEARCH_ONLY_PAIRED_STUDY_ARM2_CLEARED_B2_EXECUTION_SUPPORTED_FALSE')
    return (True, 'CLEARED_ARM2_COMPLETE_DEPLETION_GATE0D_AND_FCR_ADCC_OVERLAY')


def evaluate_paired_adcc_and_postdepletion_report_iip_0028(
    obs: dict[str, Any],
    *,
    root: Path | str | None = None,
    upstream_b2_authorization_record: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Structured report wrapper around `evaluate_paired_adcc_and_postdepletion_gate_iip_0028` preserving `execution_supported` (IIP-0028; DEC-106)."""
    ok, code = evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        obs,
        root=root,
        upstream_b2_authorization_record=upstream_b2_authorization_record,
    )
    prod_id = obs.get('product_variant_id') if isinstance(obs, dict) else None
    is_b2 = prod_id == 'ISLET1_B2_5LOCUS_BASE'
    return {
        'gate_pass': ok,
        'status': code,
        'product_variant_id': prod_id if isinstance(prod_id, str) else None,
        'execution_supported': False if (not ok or is_b2) else True,
        'paired_study_research_only': bool(ok and is_b2),
    }


def evaluate_dec084_measured_potency_mixture_gate_iip_0028(
    obs: dict[str, Any],
) -> dict[str, Any]:
    """
    Measured-potency equal-output DEC-084 mixture (70:20:10 / 70:25:5 vs 60:30:10) evaluator (IIP-0028; DEC-106):
      - Requires n >= 8 independent confirmation blocks disjoint from n >= 3 fixed-input selection blocks;
      - Requires measured starting lineage fractions, viable counts, and fresh-lot equal-output dose-response;
      - Enforces implanted loading == (6/7) * (out_A / out_B) * 1500.0 IEQ/cm^2, U95_t(OCR_B/A) <= 0.90,
        dual-mesh (40x40x70 & 60x60x105) worst-face pO2 >= 15.0 mmHg, L95_t(Cpep_AUC_B/A) >= 0.85,
        L95_t(GCG_AUC_B/A) >= 0.85, Day-28 L95_t(beta_recovery) >= 0.80, and Readout 2 insulin/glucagon bounds.
    """
    if not isinstance(obs, dict):
        return {'gate_pass': False, 'status': 'HOLD_DEC084_INVALID_PAYLOAD'}
    if obs.get('selection_blocks_excluded_from_confirmation') is not True:
        return {'gate_pass': False, 'status': 'HOLD_DEC084_SELECTION_BLOCKS_NOT_EXCLUDED'}
    if obs.get('measured_starting_composition_and_viable_counts_verified') is not True:
        return {'gate_pass': False, 'status': 'HOLD_DEC084_UNVERIFIED_STARTING_COMPOSITION_OR_VIABLE_COUNTS'}
    n_conf = obs.get('n_confirmation_blocks')
    if not _is_strict_int_iip0024(n_conf) or n_conf < 8:
        return {'gate_pass': False, 'status': 'HOLD_DEC084_CONFIRMATION_BLOCKS_LT_8'}
    t_crit = student_t_ppf(0.95, n_conf - 1)

    pot = obs.get('obs_dec084_starting_potency_ratio_out_b_over_out_a')
    pot_se = _read_se_alias_iip0028(obs, 'obs_dec084_starting_potency_ratio_out_b_over_out_a_se', 'se_starting_potency_ratio')
    load = obs.get('implanted_areal_loading_ieq_cm2')
    if not _valid_domain_scalar_iip0028(pot, 0.2, 5.0) or not _valid_positive_se_iip0028(pot_se) or not _valid_domain_scalar_iip0028(load, 100.0, 5000.0):
        return {'gate_pass': False, 'status': 'HOLD_DEC084_MISSING_OR_INVALID_POTENCY_OR_LOADING'}
    expected_load = (6.0 / 7.0) * (1.0 / float(pot)) * 1500.0
    if abs(float(load) - expected_load) > 0.25:
        return {'gate_pass': False, 'status': 'HOLD_DEC084_LOADING_DOES_NOT_MATCH_MEASURED_POTENCY_FORMULA'}

    po2_nx40 = obs.get('dual_band_worst_face_po2_nx40_mmhg')
    po2_nx60 = obs.get('dual_band_worst_face_po2_nx60_mmhg')
    if not _valid_domain_scalar_iip0028(po2_nx40, 0.0, 120.0) or not _valid_domain_scalar_iip0028(po2_nx60, 0.0, 120.0):
        return {'gate_pass': False, 'status': 'HOLD_DEC084_MISSING_MULTIMESH_WORST_FACE_OXYGEN'}
    if min(float(po2_nx40), float(po2_nx60)) < 15.0:
        if max(float(po2_nx40), float(po2_nx60)) >= 15.0:
            return {'gate_pass': False, 'status': 'HOLD_DEC084_NUMERICALLY_MESH_SENSITIVE_WORST_FACE_OXYGEN_BELOW_15MMHG'}
        return {'gate_pass': False, 'status': 'HOLD_DEC084_WORST_FACE_OXYGEN_BELOW_15MMHG'}

    bounds = [
        ('obs_dec084_ocr_ratio_b_over_a', ('obs_dec084_ocr_ratio_b_over_a_se', 'se_ocr_ratio'), 'upper', 0.90, 0.1, 5.0),
        ('obs_dec084_cpep_auc_ratio_b_over_a', ('obs_dec084_cpep_auc_ratio_b_over_a_se', 'se_cpep_auc_ratio'), 'lower', 0.85, 0.1, 5.0),
        ('obs_dec084_gcg_auc_ratio_b_over_a', ('obs_dec084_gcg_auc_ratio_b_over_a_se', 'se_gcg_auc_ratio'), 'lower', 0.85, 0.1, 5.0),
        ('obs_day28_beta_recovery_ratio', ('obs_day28_beta_recovery_se', 'se_day28_beta_recovery'), 'lower', 0.80, 0.0, 2.0),
        ('obs_readout2_perifusion_gsis_si', ('obs_readout2_perifusion_gsis_se', 'se_readout2_si'), 'lower', 3.0, 0.1, 50.0),
        ('obs_readout2_first_phase_peak_to_basal_ratio', ('obs_readout2_first_phase_peak_to_basal_se', 'se_readout2_peak_basal'), 'lower', 4.5, 0.1, 100.0),
        ('obs_readout2_first_phase_time_to_peak_min', ('obs_readout2_first_phase_time_to_peak_se', 'se_readout2_tpeak'), 'upper', 8.0, 0.1, 60.0),
        ('obs_readout2_stepdown_shutoff_t50_min', ('obs_readout2_stepdown_shutoff_t50_se', 'se_readout2_t50'), 'upper', 3.0, 0.05, 60.0),
        ('obs_gate0m_glucagon_2mM_fold_increase', ('obs_gate0m_glucagon_2mM_fold_se', 'se_readout2_gcg_fold'), 'lower', 2.0, 0.1, 50.0),
    ]
    for mk, sks, direction, thresh, lo, hi in bounds:
        m = obs.get(mk)
        s = _read_se_alias_iip0028(obs, sks[0], *sks[1:])
        if not _valid_domain_scalar_iip0028(m, lo, hi) or not _valid_positive_se_iip0028(s):
            return {'gate_pass': False, 'status': f'HOLD_DEC084_MISSING_OR_INVALID_{mk}'}
        mf, sf = float(m), float(s)
        if direction == 'upper' and (mf + t_crit * sf) > thresh:
            return {'gate_pass': False, 'status': f'HOLD_DEC084_EXCEEDED_UPPER95_{mk}'}
        if direction == 'lower' and (mf - t_crit * sf) < thresh:
            return {'gate_pass': False, 'status': f'HOLD_DEC084_BELOW_LOWER95_{mk}'}

    return {
        'gate_pass': True,
        'status': 'CLEARED_DEC084_MEASURED_POTENCY_LOW_OCR_AND_GLUCAGON_CONFIRMATION',
        'computed_equal_output_loading_ieq_cm2': round(expected_load, 4),
    }


def _validate_iapp_control_block_iip0028(rec: dict[str, Any], sample_id: str, tp_label: str) -> tuple[bool, str]:
    loq = rec.get('loq_signal')
    if not _is_strict_finite_number_iip0024(loq) or float(loq) <= 0.0:
        return (False, f'HOLD_INCONCLUSIVE_IAPP_INVALID_LOQ_{tp_label}')
    loq_f = float(loq)

    ko_id = rec.get('ko_blank_id')
    ko_sig = rec.get('ko_blank_signal')
    ins_sig = rec.get('insulin_fibril_blank_signal')
    if not _valid_strict_unpadded_id_iip0028(ko_id):
        return (False, f'HOLD_INCONCLUSIVE_IAPP_MISSING_KO_BLANK_ID_{tp_label}')
    if not _is_strict_finite_number_iip0024(ko_sig) or not (0.0 <= float(ko_sig) < loq_f):
        return (False, f'HOLD_INCONCLUSIVE_IAPP_BLANK_ABOVE_LOQ_{tp_label}')
    if not _is_strict_finite_number_iip0024(ins_sig) or not (0.0 <= float(ins_sig) < loq_f):
        return (False, f'HOLD_INCONCLUSIVE_IAPP_INSULIN_FIBRIL_BLANK_ABOVE_LOQ_{tp_label}')

    ref_id = rec.get('age_matched_ref_control_id')
    ref_sig = rec.get('ref_signal')
    if not _valid_strict_unpadded_id_iip0028(ref_id) or ref_id in (ko_id, sample_id):
        return (False, f'HOLD_INCONCLUSIVE_IAPP_INVALID_OR_KO_REUSED_REF_ID_{tp_label}')
    if not _is_strict_finite_number_iip0024(ref_sig) or not (float(ref_sig) > loq_f > 0.0):
        return (False, f'HOLD_INCONCLUSIVE_IAPP_REFERENCE_NOT_ABOVE_LOQ_{tp_label}')

    spike_pct = rec.get('spike_recovery_iapp_1_37_pct')
    if not _valid_domain_scalar_iip0028(spike_pct, 80.0, 120.0):
        return (False, f'HOLD_INCONCLUSIVE_IAPP_SPIKE_RECOVERY_OUT_OF_RANGE_{tp_label}')

    return (True, 'OK')


def evaluate_orthogonal_iapp_3way_gate_iip_0028(
    obs: dict[str, Any],
    *,
    root: Path | str | None = None,
) -> dict[str, Any]:
    """
    Visit/lot/sister-bound & concordance/replication-gated 72-h / Day-28 / Day-84 dual-modality IAPP evaluator (IIP-0028; DEC-106):
      - Validates product_variant_id against the versioned product registry;
      - Derives n >= 6 strictly from distinct verified contributing_lot_ids (rejecting free integer n);
      - Binds 72h (timepoint_h == 72) and d28 (timepoint_d == 28) to LOT_SCOPE and destructive d84 (timepoint_d == 84)
        to TERMINAL_SISTER_DEVICE_SCOPE on an authorized_sister_assembly_id != implant_device_assembly_id, rejecting
        duplicate sample_ids OR copied raw measurement vectors across visits as HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE;
      - Evaluates viable beta mass (>= 80.0% L95_t) and stimulated GSIS SI (>= 2.5 L95_t) alongside IAPP so a falling
        amyloid signal caused by beta-cell attrition cannot clear as clean proteostasis;
      - Requires within-lot dual-modality L95_t co-exceedance OR disjoint >= 6-lot repeat confirmation before emitting
        CONFIRMED_CHRONIC_P2_ELIGIBLE, holding archival NM_021999.5 ITM2B at 0% active weight.
    """
    def _fail(disp: str, **extra: Any) -> dict[str, Any]:
        return {
            'core_cleared': False,
            'disposition': disp,
            'p2_itm2b_reserve_eligible': False,
            'p2_active_weight_pct': 0.0,
            **extra,
        }

    if not isinstance(obs, dict):
        return _fail('HOLD_INCONCLUSIVE_IAPP_INVALID_PAYLOAD')

    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    p_mod = _load_product_contract_module(root_path)
    prod_doc = p_mod.load_product_contracts(root_path)
    reg_ids = {str(p['product_variant_id']) for p in prod_doc['products']}

    prod_id = obs.get('product_variant_id')
    if not _valid_strict_unpadded_id_iip0028(prod_id) or prod_id not in reg_ids:
        return _fail('HOLD_INCONCLUSIVE_IAPP_UNREGISTERED_OR_MISSING_PRODUCT_ID')

    lot_list = obs.get('contributing_lot_ids')
    if not isinstance(lot_list, (list, tuple)) or len(lot_list) < 6:
        return _fail('HOLD_IAPP_UNVERIFIED_INDEPENDENT_LOTS_LT_6')
    if not all(_valid_strict_unpadded_id_iip0028(x) for x in lot_list) or len(set(lot_list)) != len(lot_list):
        return _fail('HOLD_IAPP_DUPLICATE_OR_INVALID_LOT_IDS')
    n_rep = len(lot_list)
    if 'n_independent_blocks' in obs:
        nb = obs['n_independent_blocks']
        if not _is_strict_int_iip0024(nb) or nb != n_rep:
            return _fail('HOLD_IAPP_DECLARED_N_MISMATCH_WITH_UNIQUE_LOTS')

    primary_lot_id = obs.get('primary_lot_id')
    if not _valid_strict_unpadded_id_iip0028(primary_lot_id) or primary_lot_id not in lot_list:
        return _fail('HOLD_INCONCLUSIVE_IAPP_MISSING_OR_UNLINKED_PRIMARY_LOT_ID')

    implant_dev_id = obs.get('implant_device_assembly_id')
    if not _valid_strict_unpadded_id_iip0028(implant_dev_id):
        return _fail('HOLD_INCONCLUSIVE_IAPP_MISSING_IMPLANT_ASSEMBLY_ID')

    sister_ids = obs.get('authorized_sister_assembly_ids')
    if (
        not isinstance(sister_ids, (list, tuple))
        or len(sister_ids) == 0
        or not all(_valid_strict_unpadded_id_iip0028(s) for s in sister_ids)
        or implant_dev_id in sister_ids
    ):
        return _fail('HOLD_IAPP_MISSING_OR_SELF_REFERENTIAL_SISTER_ASSEMBLIES')

    tp_records = obs.get('timepoint_records')
    if not isinstance(tp_records, dict):
        return _fail('HOLD_INCONCLUSIVE_IAPP_MISSING_TIMEPOINT_RECORDS')

    seen_sample_ids: set[str] = set()
    seen_raw_vectors: set[tuple[float, ...]] = set()
    t_crit = student_t_ppf(0.95, n_rep - 1)
    inconclusive_flags: list[str] = []
    beta_attrition_confounded = False
    chronic_tp_modality_breaches: dict[str, dict[str, list[str]]] = {
        'd28': {'mod1_elisa': [], 'mod2_fibril': []},
        'd84': {'mod1_elisa': [], 'mod2_fibril': []},
    }

    for tp in ('72h', 'd28', 'd84'):
        rec = tp_records.get(tp)
        if not isinstance(rec, dict):
            return _fail(f'HOLD_INCONCLUSIVE_IAPP_MISSING_TIMEPOINT_{tp}')
        if rec.get('product_variant_id') != prod_id or rec.get('lot_id') != primary_lot_id:
            return _fail(f'HOLD_INCONCLUSIVE_IAPP_PROVENANCE_MISMATCH_{tp}')

        sid = rec.get('sample_id')
        if not _valid_strict_unpadded_id_iip0028(sid) or sid in seen_sample_ids:
            return _fail('HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE')
        seen_sample_ids.add(str(sid))

        if tp == '72h':
            tph = rec.get('timepoint_h')
            if rec.get('scope') != 'LOT_SCOPE' or not _is_strict_int_iip0024(tph) or tph != 72 or ('timepoint_d' in rec) or ('sister_assembly_id' in rec):
                return _fail('HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE')
        elif tp == 'd28':
            tpd = rec.get('timepoint_d')
            tph = rec.get('timepoint_h')
            if (
                rec.get('scope') != 'LOT_SCOPE'
                or not _is_strict_int_iip0024(tpd)
                or tpd != 28
                or (tph is not None and (not _is_strict_int_iip0024(tph) or tph != 672))
                or ('sister_assembly_id' in rec)
            ):
                return _fail('HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE')
        else:
            tpd = rec.get('timepoint_d')
            tph = rec.get('timepoint_h')
            sis_id = rec.get('sister_assembly_id')
            if rec.get('scope') != 'TERMINAL_SISTER_DEVICE_SCOPE':
                return _fail('HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE')
            if not _is_strict_int_iip0024(tpd) or tpd != 84 or (tph is not None and (not _is_strict_int_iip0024(tph) or tph != 2016)):
                return _fail('HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE')
            if not _valid_strict_unpadded_id_iip0028(sis_id) or sis_id == implant_dev_id or sis_id not in sister_ids:
                return _fail('HOLD_IAPP_D84_REQUIRES_AUTHORIZED_SISTER_ASSEMBLY')

        ctrl_ok, ctrl_code = _validate_iapp_control_block_iip0028(rec, str(sid), tp)
        if not ctrl_ok:
            return _fail(ctrl_code)

        # Require viable beta mass (>= 80% L95_t) & stimulated secretion (>= 2.5 SI L95_t) alongside IAPP:
        vbm = rec.get('viable_beta_mass_pct')
        vbms = rec.get('se_viable_beta_mass')
        gsi = rec.get('stimulated_gsis_si')
        gsis = rec.get('se_gsis_si')
        if not _valid_domain_scalar_iip0028(vbm, 0.0, 150.0) or not _valid_positive_se_iip0028(vbms):
            return _fail(f'HOLD_INCONCLUSIVE_IAPP_MISSING_VIABLE_BETA_MASS_{tp}')
        if not _valid_domain_scalar_iip0028(gsi, 0.1, 100.0) or not _valid_positive_se_iip0028(gsis):
            return _fail(f'HOLD_INCONCLUSIVE_IAPP_MISSING_STIMULATED_SECRETION_{tp}')
        if (float(vbm) - t_crit * float(vbms)) < 80.0 or (float(gsi) - t_crit * float(gsis)) < 2.5:
            beta_attrition_confounded = True

        metrics = [
            ('elisa_oligomer_fold', 'se_elisa_fold', 1.5, 0.0, 100.0, 'mod1_elisa'),
            ('colocalized_fibril_fold', 'se_fibril_fold', 1.5, 0.0, 100.0, 'mod2_fibril'),
            ('colocalized_amyloid_area_pct', 'se_area_pct', 1.0, 0.0, 100.0, 'mod2_fibril'),
        ]
        vec_vals: list[float] = []
        for mk, sk, cap, lo, hi, mod_group in metrics:
            m, s = rec.get(mk), rec.get(sk)
            if not _valid_domain_scalar_iip0028(m, lo, hi) or not _valid_positive_se_iip0028(s):
                return _fail(f'HOLD_INCONCLUSIVE_IAPP_INVALID_METRIC_{tp}_{mk}')
            mf, sf = float(m), float(s)
            vec_vals.extend([round(mf, 9), round(sf, 9)])
            l95 = mf - t_crit * sf
            u95 = mf + t_crit * sf
            if l95 > cap:
                if tp == '72h':
                    return _fail(f'FAIL_ACUTE_PROTEOSTASIS_72H_{mk}')
                chronic_tp_modality_breaches[tp][mod_group].append(f'{tp}_{mk}')
            elif u95 > cap:
                inconclusive_flags.append(f'{tp}_{mk}')

        raw_vec = tuple(vec_vals)
        if raw_vec in seen_raw_vectors:
            return _fail('HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE')
        seen_raw_vectors.add(raw_vec)

    if any(f.startswith('72h_') for f in inconclusive_flags):
        return _fail('HOLD_INCONCLUSIVE_IAPP_EVIDENCE_72H_UPPER_BOUND', inconclusive_upper_bounds=inconclusive_flags)

    concordant_chronic_tps: list[str] = []
    discordant_chronic_tps: list[tuple[str, str, list[str]]] = []
    all_chronic_breaches: list[str] = []
    for tp in ('d28', 'd84'):
        m1 = chronic_tp_modality_breaches[tp]['mod1_elisa']
        m2 = chronic_tp_modality_breaches[tp]['mod2_fibril']
        all_chronic_breaches.extend(m1 + m2)
        if m1 and m2:
            if tp_records[tp].get('within_lot_dual_modality_coexceedance_verified', True) is True:
                concordant_chronic_tps.append(tp)
            else:
                discordant_chronic_tps.append((tp, 'mod1_elisa', m1 + m2))
        elif m1 or m2:
            which_mod = 'mod1_elisa' if m1 else 'mod2_fibril'
            discordant_chronic_tps.append((tp, which_mod, m1 if m1 else m2))

    if concordant_chronic_tps:
        return {
            'core_cleared': False,
            'disposition': 'CONFIRMED_CHRONIC_P2_ELIGIBLE',
            'confirmation_route': 'CONCORDANT_WITHIN_LOT_DUAL_MODALITY_L95_BREACH',
            'concordant_timepoints': concordant_chronic_tps,
            'confirmed_chronic_breaches': all_chronic_breaches,
            'p2_itm2b_reserve_eligible': True,
            'p2_active_weight_pct': 0.0,
            'p2_b2_branch_hold_pending_b2_registry_qc': True,
        }

    if discordant_chronic_tps:
        rep = obs.get('independent_repeat_confirmation')
        if not isinstance(rep, dict):
            return _fail(
                'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION',
                discordant_chronic_breaches=[b for _, _, bs in discordant_chronic_tps for b in bs],
            )
        rep_proto = rep.get('preregistered_replication_protocol_id')
        rep_lots = rep.get('repeat_lot_ids')
        rep_sid = rep.get('repeat_sample_id')
        rep_tp = rep.get('timepoint')
        if (
            not _valid_strict_unpadded_id_iip0028(rep_proto)
            or not isinstance(rep_lots, (list, tuple))
            or len(rep_lots) < 6
            or not all(_valid_strict_unpadded_id_iip0028(x) for x in rep_lots)
            or len(set(rep_lots)) != len(rep_lots)
            or not set(rep_lots).isdisjoint(set(lot_list))
            or not _valid_strict_unpadded_id_iip0028(rep_sid)
            or rep_sid in seen_sample_ids
        ):
            return _fail(
                'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION',
                discordant_chronic_breaches=[b for _, _, bs in discordant_chronic_tps for b in bs],
                replication_rejection_reason='MISSING_OR_NON_DISJOINT_REPEAT_LOTS_OR_SAMPLE',
            )
        target_tp, target_mod, target_breaches = discordant_chronic_tps[-1]
        expected_tp_d = 28 if target_tp == 'd28' else 84
        if rep_tp != target_tp or rep.get('timepoint_d') != expected_tp_d:
            return _fail(
                'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION',
                discordant_chronic_breaches=target_breaches,
                replication_rejection_reason='REPEAT_TIMEPOINT_MISMATCH',
            )
        if target_tp == 'd84':
            rep_sis = rep.get('repeat_sister_assembly_id')
            if not _valid_strict_unpadded_id_iip0028(rep_sis) or rep_sis == implant_dev_id or rep_sis in sister_ids:
                return _fail(
                    'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION',
                    discordant_chronic_breaches=target_breaches,
                    replication_rejection_reason='REPEAT_D84_REQUIRES_INDEPENDENT_SISTER_ASSEMBLY',
                )
        rep_ctrl_ok, _ = _validate_iapp_control_block_iip0028(rep, str(rep_sid), f'repeat_{target_tp}')
        if not rep_ctrl_ok:
            return _fail(
                'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION',
                discordant_chronic_breaches=target_breaches,
                replication_rejection_reason='REPEAT_CONTROL_BLOCK_INVALID',
            )
        t_rep = student_t_ppf(0.95, len(rep_lots) - 1)
        rep_confirmed = False
        if target_mod == 'mod1_elisa':
            rm, rs = rep.get('elisa_oligomer_fold'), rep.get('se_elisa_fold')
            if _valid_domain_scalar_iip0028(rm, 0.0, 100.0) and _valid_positive_se_iip0028(rs) and (float(rm) - t_rep * float(rs)) > 1.5:
                rep_confirmed = True
        else:
            rfm, rfs = rep.get('colocalized_fibril_fold'), rep.get('se_fibril_fold')
            ram, ras = rep.get('colocalized_amyloid_area_pct'), rep.get('se_area_pct')
            if _valid_domain_scalar_iip0028(rfm, 0.0, 100.0) and _valid_positive_se_iip0028(rfs) and (float(rfm) - t_rep * float(rfs)) > 1.5:
                rep_confirmed = True
            if _valid_domain_scalar_iip0028(ram, 0.0, 100.0) and _valid_positive_se_iip0028(ras) and (float(ram) - t_rep * float(ras)) > 1.0:
                rep_confirmed = True

        if not rep_confirmed:
            return _fail(
                'HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION',
                discordant_chronic_breaches=target_breaches,
                replication_rejection_reason='REPEAT_L95_DID_NOT_CONFIRM_BREACH',
            )

        return {
            'core_cleared': False,
            'disposition': 'CONFIRMED_CHRONIC_P2_ELIGIBLE',
            'confirmation_route': 'INDEPENDENT_DISJOINT_REPEAT_LOT_BATCH',
            'confirmed_chronic_breaches': target_breaches,
            'repeat_lot_count': len(rep_lots),
            'p2_itm2b_reserve_eligible': True,
            'p2_active_weight_pct': 0.0,
            'p2_b2_branch_hold_pending_b2_registry_qc': True,
        }

    if beta_attrition_confounded:
        return _fail('HOLD_IAPP_CONFOUNDED_BY_VIABLE_BETA_MASS_OR_SECRETION_LOSS')

    if inconclusive_flags:
        return _fail('HOLD_INCONCLUSIVE_IAPP_EVIDENCE', inconclusive_upper_bounds=inconclusive_flags)

    return {
        'core_cleared': True,
        'disposition': 'CLEARED_DUAL_ORTHOGONAL_IAPP_72H_D28_D84',
        'verified_unique_lots': n_rep,
        'p2_itm2b_reserve_eligible': False,
        'p2_active_weight_pct': 0.0,
    }


# ============================================================================
# Separate Resting (0 h) & Interferon-Stimulated (24 h) CXCL10 Promoter-of-
# Origin 5'-RACE TSS Attribution & Conditional BIO-06 AAV-ITR Origin Gate
# (IIP-0027; DEC-107)
# ============================================================================
def evaluate_cxcl10_promoter_origin_and_itr_gate_iip_0027(
    observations: dict[str, Any],
    *,
    validated_stage0: dict[str, float] | None = None,
    non_integer_denominator_flagged: bool = False,
) -> dict[str, Any]:
    """
    Fail-closed Stage-0 resting (0 h) + interferon-stimulated (24 h) CXCL10 promoter-of-origin
    5'-RACE TSS attribution and conditional BIO-06 AAV-ITR RNA-origin evaluator (IIP-0027; DEC-107).

    Enforces:
      - Named donor (cxcl10_tss_donor_id) and transcription unit (cxcl10_tss_transcription_unit_id);
      - Insert-matched comparator control verification (Giannakis et al. 2003 PMID:12951775);
      - Separate resting (0 h) AND interferon-stimulated (24 h) 5'-RACE TSS attribution classifying
        expected host / designed internal vs unexpected internal (k == 0) and unresolved (k == 0)
        origins, with recovery-adjusted N_UMI_eff >= 3,689 wherever expressed (3,688 exact two-sided
        95% Clopper-Pearson minimum for < 0.10% at k = 0; Kivioja et al. 2011 PMID:22101854);
      - Conditional Readout 1 BIO-06 AAV-ITR RNA-origin attribution (Earley et al. 2020 PMID:31914802;
        Haberman et al. 2000 PMID:10954575; Flotte et al. 1993 PMID:7679117) whenever locus-specific
        junction sequencing documents retained ITR DNA (obs_bio06_retained_aav_itr_dna_pct > 0%,
        alongside the < 0.5% DNA capture ceiling from Hanlon et al. 2019 PMID:31570731).
    """
    if not isinstance(observations, dict):
        raise ValueError('IIP-0027 promoter-origin observations must be a dictionary (IIP-0027; DEC-107)')

    def _get_num(key: str) -> float | None:
        if validated_stage0 is not None and key in validated_stage0:
            return float(validated_stage0[key])
        raw = observations.get(key)
        if raw is None or raw == 'UNAVAILABLE':
            return None
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            return None
        val = float(raw)
        return val if math.isfinite(val) else None

    def _get_nonneg_int(key: str) -> tuple[bool, int | None]:
        raw = observations.get(key)
        if raw is None or raw == 'UNAVAILABLE':
            if validated_stage0 is not None and key in validated_stage0:
                vf = float(validated_stage0[key])
                if vf >= 0.0 and abs(vf - round(vf)) <= 1e-9:
                    return True, int(round(vf))
                return False, None
            return True, None
        if isinstance(raw, bool):
            return False, None
        if isinstance(raw, int):
            return (raw >= 0), (raw if raw >= 0 else None)
        if isinstance(raw, float) and math.isfinite(raw) and raw >= 0.0 and abs(raw - round(raw)) <= 1e-9:
            return True, int(round(raw))
        return False, None

    def _unpadded_id(val: object) -> bool:
        return isinstance(val, str) and val == val.strip() and 3 <= len(val) <= 120

    def _result(status: str, gate_pass: bool, **extra: Any) -> dict[str, Any]:
        return {
            'gate_pass': gate_pass,
            'hold_status': None if gate_pass else status,
            'status': status,
            'cxcl10_tss_donor_id': observations.get('cxcl10_tss_donor_id'),
            'cxcl10_tss_transcription_unit_id': observations.get('cxcl10_tss_transcription_unit_id'),
            'oq011_closed_by_computational_motif_count': False,
            'requires_empirical_race_per_condition_and_unit': True,
            **extra,
        }

    # Conditional BIO-06 ITR checks can also fire independently if retained ITR DNA or ITR RNA is observed:
    itr_dna_pct = _get_num('obs_bio06_retained_aav_itr_dna_pct')
    itr_rna_ok, itr_rna_k = _get_nonneg_int('obs_bio06_itr_origin_rna_umi_k')
    if not itr_rna_ok:
        return _result('HOLD_IIP0027_INVALID_NON_INTEGER_UMI_COUNT', False)
    if itr_dna_pct is not None and itr_dna_pct >= 0.5:
        return _result('HOLD_FAIL_BIO06_AAV_ITR_DNA_CAPTURE_GE_0P5PCT', False, retained_itr_dna_pct=itr_dna_pct)
    if itr_rna_k is not None and itr_rna_k > 0:
        return _result('HOLD_FAIL_IIP0027_AAV_ITR_ORIGIN_TRANSCRIPT_DETECTED', False, itr_origin_rna_umi_k=itr_rna_k)
    if itr_dna_pct is not None and itr_dna_pct > 0.0:
        if observations.get('bio06_itr_rna_origin_attribution_verified') is not True or itr_rna_k is None:
            return _result(
                'HOLD_IIP0027_RETAINED_AAV_ITR_DNA_REQUIRES_RNA_ORIGIN_ATTRIBUTION',
                False,
                retained_itr_dna_pct=itr_dna_pct,
            )

    # Check if only clean BIO-06 ITR check (0% DNA or verified 0 ITR RNA) was supplied without CXCL10 promoter-origin fields:
    cxcl10_keys_present = any(
        observations.get(k) is not None or (validated_stage0 is not None and k in validated_stage0)
        for k in (
            'require_full_iip0027_promoter_origin_and_itr_panel',
            'cxcl10_insert_matched_comparator_control_verified',
            'cxcl10_resting_expression_below_loq_confirmed',
            'cxcl10_tss_donor_id',
            'cxcl10_tss_transcription_unit_id',
            'obs_cxcl10_resting_race_umi_n_eff',
            'obs_cxcl10_resting_unexpected_tss_umi_k',
            'obs_cxcl10_resting_unresolved_origin_umi_k',
            'obs_cxcl10_ifng24h_race_umi_n_eff',
            'obs_cxcl10_ifng24h_unexpected_tss_umi_k',
            'obs_cxcl10_ifng24h_unresolved_origin_umi_k',
        )
    )
    if not cxcl10_keys_present:
        return _result('CLEARED_IIP0027_BIO06_CONDITIONAL_AAV_ITR_ORIGIN_CHECK', True, retained_itr_dna_pct=itr_dna_pct, itr_origin_rna_umi_k=itr_rna_k)

    if non_integer_denominator_flagged:
        return _result('HOLD_IIP0027_INVALID_NON_INTEGER_UMI_COUNT', False)
    if not (_unpadded_id(observations.get('cxcl10_tss_donor_id')) and _unpadded_id(observations.get('cxcl10_tss_transcription_unit_id'))):
        return _result('HOLD_IIP0027_MISSING_NAMED_DONOR_OR_TRANSCRIPTION_UNIT', False)
    if observations.get('cxcl10_insert_matched_comparator_control_verified') is not True:
        return _result('HOLD_IIP0027_MISSING_INSERT_MATCHED_COMPARATOR_CONTROL', False)

    r_n_ok, r_n = _get_nonneg_int('obs_cxcl10_resting_race_umi_n_eff')
    r_u_ok, r_u = _get_nonneg_int('obs_cxcl10_resting_unexpected_tss_umi_k')
    r_q_ok, r_q = _get_nonneg_int('obs_cxcl10_resting_unresolved_origin_umi_k')
    s_n_ok, s_n = _get_nonneg_int('obs_cxcl10_ifng24h_race_umi_n_eff')
    s_u_ok, s_u = _get_nonneg_int('obs_cxcl10_ifng24h_unexpected_tss_umi_k')
    s_q_ok, s_q = _get_nonneg_int('obs_cxcl10_ifng24h_unresolved_origin_umi_k')
    if not (r_n_ok and r_u_ok and r_q_ok and s_n_ok and s_u_ok and s_q_ok):
        return _result('HOLD_IIP0027_INVALID_NON_INTEGER_UMI_COUNT', False)

    if r_n is None or r_u is None or r_q is None:
        return _result('HOLD_IIP0027_MISSING_RESTING_0H_PROMOTER_ORIGIN_5PRIME_RACE', False)
    if s_n is None or s_u is None or s_q is None:
        return _result('HOLD_IIP0027_MISSING_IFNG24H_STIMULATED_PROMOTER_ORIGIN_5PRIME_RACE', False)

    if r_u > 0:
        return _result('HOLD_FAIL_IIP0027_RESTING_UNEXPECTED_INTERNAL_TSS_DETECTED', False, resting_unexpected_tss_umi_k=r_u)
    if s_u > 0:
        return _result('HOLD_FAIL_IIP0027_IFNG24H_UNEXPECTED_INTERNAL_TSS_DETECTED', False, ifng24h_unexpected_tss_umi_k=s_u)
    if r_q > 0:
        return _result('HOLD_IIP0027_RESTING_UNRESOLVED_TSS_ORIGIN', False, resting_unresolved_origin_umi_k=r_q)
    if s_q > 0:
        return _result('HOLD_IIP0027_IFNG24H_UNRESOLVED_TSS_ORIGIN', False, ifng24h_unresolved_origin_umi_k=s_q)

    resting_below_loq = observations.get('cxcl10_resting_expression_below_loq_confirmed') is True
    if not resting_below_loq and r_n < 3689:
        return _result('HOLD_IIP0027_RESTING_UMI_BELOW_3689', False, resting_race_umi_n_eff=r_n)
    if s_n < 3689:
        return _result('HOLD_IIP0027_IFNG24H_UMI_BELOW_3689', False, ifng24h_race_umi_n_eff=s_n)

    return _result(
        'CLEARED_IIP0027_RESTING_AND_STIMULATED_PROMOTER_ORIGIN_AND_ITR_GATES',
        True,
        resting_race_umi_n_eff=r_n,
        resting_below_loq_confirmed=resting_below_loq,
        ifng24h_race_umi_n_eff=s_n,
        retained_itr_dna_pct=itr_dna_pct,
        itr_origin_rna_umi_k=itr_rna_k,
    )


# ============================================================================
# Stage-0 ALK5i-Free Two-Clock Maturation, Two-Pulse Readout 2 with Adult-Islet
# Reference Floor, Inert-Tracer tau_chamber, and Non-Bypassable >138.3 IU Gate 0S
# Atomic-Ledger Escalation (IIP-0029; DEC-108)
# ============================================================================
_HEX64_RE_IIP0029 = re.compile(r'^[0-9a-fA-F]{64}$')

STAGE0_IIP0029_CONSTANTS: MappingProxyType = _deep_freeze_iip0024({
    'MIN_INDEPENDENT_LOTS_N': 6,
    'MIN_DEC084_MIXTURE_LOTS_N': 8,
    'MIN_PILOT_SELECTION_LOTS_N': 3,
    'STAGE7_MIN_TOTAL_DAYS': 7.0,
    'STAGE7_MAX_TOTAL_DAYS': 49.0,
    'SORTED_MIN_SORT_DAY': 7.0,
    'SORTED_MIN_POST_SORT_DAYS': 5.0,
    'BASELINE_WINDOW_MAX_STAGE7_DAYS': 28.0,
    'BASELINE_WINDOW_MAX_POST_SORT_DAYS': 14.0,
    'TIMESTAMP_DAY_TOLERANCE': 0.25,
    'SORTED_ETA_MIN': 0.05,
    'SORTED_ETA_MAX_EXCLUSIVE': 0.95,
    'CONTINUOUS_ZM447439_MIN_UM': 0.45,
    'CONTINUOUS_ZM447439_MAX_UM': 0.55,
    'PRE_IEC_ZM447439_MAX_LLOQ_UM': 0.005,
    'PRE_IEC_MIN_WASH_CYCLES': 3,
    'PRE_IEC_SPIKE_RECOVERY_MIN_PCT': 80.0,
    'PRE_IEC_SPIKE_RECOVERY_MAX_PCT': 120.0,
    'ALLOWED_ASSAY_PO2_MMHG': (40.0, 100.0),
    'INTEGRATED_SI1_MIN': 3.0,
    'DYNAMIC_PEAK1_OVER_BASAL1_MIN': 4.5,
    'DYNAMIC_TIME_TO_PEAK1_MAX_MIN': 8.0,
    'GATE0M_BASAL1_FRAC_PEAK1_MAX': 0.15,
    'DYNAMIC_LOW2_FRAC_PEAK1_MAX': 0.15,
    'LOW2_OVER_LOW1_WASHOUT_MAX': 1.50,
    'STATIC_LOW_OVER_KCL_MAX': 0.25,
    'GATE0M_GLUCAGON_FOLD_MIN': 2.0,
    'GATE0M_GLUCOSE_THRESHOLD_MIN_MM': 4.5,
    'GATE0M_UCN3_RATIO_MIN': 0.70,
    'GATE0M_COUNTERREG_GAIN_MIN': 0.80,
    'TRACER_TAU_CHAMBER_MAX_MIN': 0.50,
    'MIN_UNCENSORED_LOW2_DURATION_MIN': 10.0,
    'READOUT2_T50_MAX_MIN': 3.0,
    'GATE0M_T90_MAX_MIN': 3.0,
    'ADULT_ISLET_REF_PULSE2_MIN_FRACTION': 0.50,
    'REGISTERED_HORMONE_UNIT_SPECS': {
        'pmol_cpep_per_1e6_presort_cells_per_min': {
            'z_min_floor': 0.50,
            'required_denominator': '1e6_presort_viable_endocrine_cells',
        },
        'uIU_insulin_per_ug_dna_per_min': {
            'z_min_floor': 0.50,
            'required_denominator': 'ug_presort_viable_endocrine_dna',
        },
        'ng_insulin_per_1e5_presort_cells_per_30min': {
            'z_min_floor': 0.50,
            'required_denominator': '1e5_presort_viable_endocrine_cells',
        },
    },
    'REGISTERED_MIXTURES': {
        'REF_60_30_10': {'beta': 0.60, 'alpha': 0.30, 'delta': 0.10, 'is_reference': True},
        'DEC084_M1_70_20_10': {'beta': 0.70, 'alpha': 0.20, 'delta': 0.10, 'is_reference': False},
        'DEC084_M2_70_25_5': {'beta': 0.70, 'alpha': 0.25, 'delta': 0.05, 'is_reference': False},
        'DEC084_M3_80_15_5': {'beta': 0.80, 'alpha': 0.15, 'delta': 0.05, 'is_reference': False},
        'DEC084_M4_80_10_10': {'beta': 0.80, 'alpha': 0.10, 'delta': 0.10, 'is_reference': False},
    },
    'STORED_INSULIN_FLOOR_AT_300K_BETA_IU': 41.5,
    'GATE0S_ABSOLUTE_STORED_INSULIN_CEILING_IU': 138.3,
    'PROINSULIN_TO_CPEPTIDE_CONTENT_MOLAR_MAX': 0.10,
    'S8_WITHIN_RECIPE_FIELLER_THETA0': 0.85,
})

_FORBIDDEN_IIP0029_OVERRIDE_KEYS: frozenset[str] = frozenset({
    'execution_supported',
    'clinical_release_eligible',
    'gate0s_elevated_stored_insulin_clamp_reverified',
    'override_gate0s_ceiling',
    'bypass_atomic_ledger_gate0s',
})


def _strict_finite_float_iip0029(val: Any) -> float | None:
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        return None
    fv = float(val)
    return fv if math.isfinite(fv) else None


def _strict_pos_int_iip0029(val: Any) -> int | None:
    if isinstance(val, bool) or not isinstance(val, int) or val <= 0:
        return None
    return int(val)


def _strict_unpadded_str_iip0029(val: Any, min_len: int = 1) -> str | None:
    if not isinstance(val, str):
        return None
    if len(val) < min_len or val != val.strip() or not val.strip():
        return None
    return val


def _valid_construct_sha256_iip0029(val: Any) -> bool:
    if not isinstance(val, str) or not _HEX64_RE_IIP0029.fullmatch(val):
        return False
    return val != ('0' * 64)


def _parse_utc_iso_iip0029(val: Any) -> datetime | None:
    s = _strict_unpadded_str_iip0029(val, 16)
    if s is None:
        return None
    try:
        iso = s[:-1] + '+00:00' if s.endswith('Z') else s
        dt = datetime.fromisoformat(iso)
    except ValueError:
        return None
    if dt.tzinfo is None:
        return None
    return dt.astimezone(timezone.utc)


def _coerce_lot_scope_key_iip0029(raw: Any) -> LotScopeKey | None:
    if isinstance(raw, LotScopeKey):
        return raw
    if isinstance(raw, (dict, MappingProxyType)):
        req = ('product_variant_id', 'construct_sha256', 'donor_seedline_id', 'mcb_id', 'lot_id')
        if set(raw.keys()) != set(req):
            return None
        try:
            return LotScopeKey(
                raw['product_variant_id'],
                raw['construct_sha256'],
                raw['donor_seedline_id'],
                raw['mcb_id'],
                raw['lot_id'],
            )
        except (TypeError, ValueError):
            return None
    return None


def _one_sided_t_bounds_iip0029(
    vals: list[float],
    t_crit: float,
    *,
    extra_within_var: float = 0.0,
) -> tuple[float, float, float]:
    n = len(vals)
    if n < 2:
        return (vals[0] if vals else float('nan'), float('-inf'), float('inf'))
    mu = sum(vals) / float(n)
    s2 = sum((x - mu) ** 2 for x in vals) / float(n - 1)
    if s2 <= 0.0:
        # Degenerate zero-between-lot-variance arrays cannot pass uncertainty gates:
        return mu, float('-inf'), float('inf')
    se = math.sqrt(max(0.0, s2 / float(n) + max(0.0, float(extra_within_var))))
    return mu, mu - t_crit * se, mu + t_crit * se


def _extract_iip0029_source_fingerprint(record: Any) -> tuple[Any, ...]:
    if not isinstance(record, (dict, MappingProxyType)):
        return ()
    sk = _coerce_lot_scope_key_iip0029(record.get('lot_scope_key'))
    csha = sk.construct_sha256 if sk is not None else None
    lids = tuple(record.get('confirmation_lot_ids') or ()) if isinstance(record.get('confirmation_lot_ids'), (list, tuple)) else ()
    lot_obs = record.get('lot_observations') or ()
    tids: list[str] = []
    if isinstance(lot_obs, (list, tuple)):
        for it in lot_obs:
            if isinstance(it, (dict, MappingProxyType)):
                tids.append(str(it.get('trace_id', '')))
                if 'matched_100mmhg_trace_id' in it:
                    tids.append(str(it.get('matched_100mmhg_trace_id', '')))
    return (
        csha,
        record.get('stage7_start_utc'),
        record.get('sort_reaggregation_utc'),
        record.get('assay_harvest_utc'),
        record.get('locked_analysis_plan_id'),
        lids,
        tuple(tids),
    )


def evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(
    record: dict[str, Any] | MappingProxyType,
    *,
    root: Path | str | None = None,
) -> dict[str, Any]:
    """
    Fail-closed Stage-0 two-clock & repeated-challenge empirical qualification subgate (IIP-0029; DEC-108).
    NEVER grants standalone clinical release authority (`release_authority` remains strictly
    'DELEGATED_TO_LIVE_ATOMIC_STAGE0_LEDGER', preserving registry `execution_supported` holds on
    ISLET1_B2_5LOCUS_BASE and ISLET1_Z_X_MOL_140UM_ZERO_CAO2).

    Enforces:
      1. Unit-bound absolute increment floors (Z1, Z2) PLUS contemporaneous same-assay primary adult
         human islet second-pulse floor (L95_t(RefAdult_Z2) >= z_min and L95_t(Z2 - 0.50 * RefAdult_Z2) >= 0)
         both inside and outside Rule S8, plus immutable within-recipe Rule S8 Fieller (theta_0 = 0.85)
         with exact Bonferroni across K windows or disjoint pilot with locked_analysis_plan_id.
      2. Preserved canonical Readout 2 (SI_1 >= 3.0, Peak1/Basal1 >= 4.5 <= 8.0 min, t50 <= 3.0 min)
         and full Gate 0M (t90 <= 3.0 min, t50_b < t90_b, Basal1/Peak1 <= 0.15 explicitly separated from
         nominal Readout 2 4.5x with window definitions, glucagon >= 2.0x, threshold >= 4.5 mM,
         UCN3 >= 0.70, counterregulation >= 0.80).
      3. Non-overridable escalation to full scope-matched Gate 0S atomic-ledger review whenever
         U95_t(stored_insulin_iu) > 138.3 IU OR max_lot > 138.3 IU.
      4. Locked UTC timestamps for dual clocks (phenotype failure remains S8 hold trigger; clocks classify
         extended exposure), independent differentiation_block_id coverage at BOTH 40 mmHg AND 100 mmHg,
         per-trace tracer alignment (tau_chamber_b <= 0.50 min AND U95_t <= 0.50 min) & censored short
         low-interval handling, common-denominator sorting+reaggregation recovery, registered-mixture +
         measured (beta, alpha, delta) DEC-084 linkage, pre-iEC LC-MS/MS washout + vascular function
         verification, mature-insulin-specific inventory assay, and explicit exploratory labeling of
         Keymeulen proinsulin:C-peptide <= 0.10.
    """
    C = STAGE0_IIP0029_CONSTANTS

    def _hold(status: str, **extra: Any) -> dict[str, Any]:
        out: dict[str, Any] = {
            'subgate_pass': False,
            'release_authority': 'DELEGATED_TO_LIVE_ATOMIC_STAGE0_LEDGER',
            'gate0s_escalation_required': False,
            'all_non_store_ceiling_checks_passed': False,
            'status': status,
        }
        out.update(extra)
        return out

    if not isinstance(record, (dict, MappingProxyType)):
        return _hold('HOLD_INVALID_RECORD_TYPE')

    # Reject any caller attempt to override registry execution_supported, release authority, or Gate 0S ceiling:
    if any(k in record for k in _FORBIDDEN_IIP0029_OVERRIDE_KEYS):
        return _hold('HOLD_FORBIDDEN_CALLER_RELEASE_OR_EXECUTION_OVERRIDE')

    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    p_mod = _load_product_contract_module(root_path)
    prod_doc = p_mod.load_product_contracts(root_path)
    product_registry = {str(p['product_variant_id']): p for p in prod_doc['products']}

    # 1. Registry-bound LotScopeKey & strict 64-hex non-zero construct_sha256:
    scope_key = _coerce_lot_scope_key_iip0029(record.get('lot_scope_key'))
    if scope_key is None:
        return _hold('HOLD_MISSING_OR_INVALID_IIP0024_LOT_SCOPE_KEY')
    pvid = _strict_unpadded_str_iip0029(scope_key.product_variant_id, 4)
    if pvid is None or pvid not in product_registry:
        return _hold('HOLD_UNREGISTERED_OR_UNSUPPORTED_PRODUCT_VARIANT_ID')
    reg_entry = product_registry[pvid]
    registry_exec_supported = bool(reg_entry.get('execution_supported', False))

    if not _valid_construct_sha256_iip0029(scope_key.construct_sha256):
        return _hold('HOLD_INVALID_OR_ZERO_CONSTRUCT_SHA256')
    if (
        _strict_unpadded_str_iip0029(scope_key.donor_seedline_id, 2) is None
        or _strict_unpadded_str_iip0029(scope_key.mcb_id, 2) is None
        or _strict_unpadded_str_iip0029(scope_key.lot_id, 2) is None
    ):
        return _hold('HOLD_INVALID_SCOPE_PROVENANCE_IDENTIFIERS')

    # 2. Registered mixture + measured (beta, alpha, delta) composition linkage to DEC-084 (Sol BOX2 Fix #4):
    reg_mix_id = _strict_unpadded_str_iip0029(record.get('registered_mixture_id'), 4)
    if reg_mix_id is None or reg_mix_id not in C['REGISTERED_MIXTURES']:
        return _hold('HOLD_MISSING_OR_UNREGISTERED_MIXTURE_ID')
    mix_spec = C['REGISTERED_MIXTURES'][reg_mix_id]
    beta_frac = _strict_finite_float_iip0029(record.get('beta_cell_fraction'))
    alpha_frac = _strict_finite_float_iip0029(record.get('alpha_cell_fraction'))
    delta_frac = _strict_finite_float_iip0029(record.get('delta_cell_fraction'))
    imp_ieq = _strict_finite_float_iip0029(record.get('implanted_total_ieq'))
    if (
        beta_frac is None
        or alpha_frac is None
        or delta_frac is None
        or imp_ieq is None
        or imp_ieq <= 0.0
        or not (0.0 < beta_frac <= 1.0)
        or not (0.0 <= alpha_frac <= 1.0)
        or not (0.0 <= delta_frac <= 1.0)
    ):
        return _hold('HOLD_MISSING_OR_INVALID_MEASURED_COMPOSITION_OR_IEQ')
    comp_sum = beta_frac + alpha_frac + delta_frac
    if not (0.95 <= comp_sum <= 1.0 + 1e-6):
        return _hold('HOLD_INVALID_MEASURED_ENDOCRINE_COMPOSITION_SUM')

    comp_within_registered_band = (
        abs(beta_frac - float(mix_spec['beta'])) <= 0.05 + 1e-9
        and abs(alpha_frac - float(mix_spec['alpha'])) <= 0.05 + 1e-9
        and abs(delta_frac - float(mix_spec['delta'])) <= 0.03 + 1e-9
    )
    if mix_spec['is_reference']:
        alpha_delta_ratio = alpha_frac / max(1e-9, delta_frac)
        if not comp_within_registered_band or not (2.0 <= alpha_delta_ratio <= 4.5):
            return _hold(
                'HOLD_MEASURED_COMPOSITION_DRIFT_FROM_REGISTERED_MIXTURE_REQUIRES_DEC084_QUALIFICATION',
                measured_composition={'beta': beta_frac, 'alpha': alpha_frac, 'delta': delta_frac},
            )
    elif not comp_within_registered_band:
        return _hold(
            'HOLD_MEASURED_COMPOSITION_OUTSIDE_REGISTERED_DEC084_TOLERANCE',
            measured_composition={'beta': beta_frac, 'alpha': alpha_frac, 'delta': delta_frac},
        )

    is_dec084_non_ref_mixture = not bool(mix_spec['is_reference'])
    n_lots = _strict_pos_int_iip0029(record.get('n_independent_lots'))
    min_n_req = int(C['MIN_DEC084_MIXTURE_LOTS_N'] if is_dec084_non_ref_mixture else C['MIN_INDEPENDENT_LOTS_N'])
    if n_lots is None or n_lots < min_n_req:
        return _hold(
            'HOLD_DEC084_MIXTURE_REQUIRES_N_GE_8_CONFIRMATION_BLOCKS'
            if is_dec084_non_ref_mixture
            else 'HOLD_INSUFFICIENT_INDEPENDENT_LOTS_REQUIRES_N_GE_6'
        )
    if is_dec084_non_ref_mixture:
        dec084_obs = record.get('dec084_mixture_observations')
        if isinstance(dec084_obs, (dict, MappingProxyType)):
            dec084_res = evaluate_dec084_measured_potency_mixture_gate_iip_0028(dict(dec084_obs))
            if dec084_res.get('gate_pass') is not True:
                return _hold(
                    'HOLD_DEC084_REQUIRES_LINKED_MATCHED_AGE_MIXTURE_GATE',
                    dec084_gate_status=dec084_res.get('status'),
                )
        elif record.get('dec084_matched_age_mixture_gate_passed') is not True:
            return _hold('HOLD_DEC084_REQUIRES_LINKED_MATCHED_AGE_MIXTURE_GATE')

    lot_ids_raw = record.get('confirmation_lot_ids')
    if not isinstance(lot_ids_raw, (list, tuple)) or len(lot_ids_raw) != n_lots:
        return _hold('HOLD_DUPLICATE_OR_INVALID_CONFIRMATION_LOT_IDS')
    if any(_strict_unpadded_str_iip0029(x, 2) is None for x in lot_ids_raw):
        return _hold('HOLD_DUPLICATE_OR_INVALID_CONFIRMATION_LOT_IDS')
    lot_ids = [str(x) for x in lot_ids_raw]
    if len(set(lot_ids)) != n_lots:
        return _hold('HOLD_DUPLICATE_OR_INVALID_CONFIRMATION_LOT_IDS')

    # 3. Locked analysis plan, explicit K, pilot disjointness, and exact Bonferroni multiplicity:
    locked_plan_id = _strict_unpadded_str_iip0029(record.get('locked_analysis_plan_id'), 4)
    if locked_plan_id is None:
        return _hold('HOLD_MISSING_LOCKED_ANALYSIS_PLAN_ID')
    if 'candidate_windows_evaluated_k' not in record:
        return _hold('HOLD_MISSING_EXPLICIT_CANDIDATE_WINDOWS_K')
    k_cand = _strict_pos_int_iip0029(record.get('candidate_windows_evaluated_k'))
    if k_cand is None or k_cand > 5:
        return _hold('HOLD_INVALID_CANDIDATE_WINDOW_COUNT')

    win_design = record.get('window_selection_design')
    pilot_block_ids_set: set[str] = set()
    if win_design == 'SINGLE_PRESPECIFIED_PRIMARY':
        if k_cand != 1:
            return _hold('HOLD_WINNERS_CURSE_UNADJUSTED_MULTIWINDOW_OR_LOT_OVERLAP')
        t_crit = student_t_ppf(0.95, float(n_lots - 1))
    elif win_design == 'PILOT_SCREEN_THEN_DISJOINT_CONFIRMATION':
        pilot_raw = record.get('pilot_selection_lot_ids')
        pilot_blk_raw = record.get('pilot_differentiation_block_ids')
        if (
            not isinstance(pilot_raw, (list, tuple))
            or len(pilot_raw) < int(C['MIN_PILOT_SELECTION_LOTS_N'])
            or any(_strict_unpadded_str_iip0029(x, 2) is None for x in pilot_raw)
        ):
            return _hold('HOLD_WINNERS_CURSE_UNADJUSTED_MULTIWINDOW_OR_LOT_OVERLAP')
        if len(set(pilot_raw)) != len(pilot_raw) or (set(pilot_raw) & set(lot_ids)):
            return _hold('HOLD_WINNERS_CURSE_UNADJUSTED_MULTIWINDOW_OR_LOT_OVERLAP')
        if (
            not isinstance(pilot_blk_raw, (list, tuple))
            or len(pilot_blk_raw) < int(C['MIN_PILOT_SELECTION_LOTS_N'])
            or any(_strict_unpadded_str_iip0029(x, 2) is None for x in pilot_blk_raw)
            or len(set(pilot_blk_raw)) != len(pilot_blk_raw)
        ):
            return _hold('HOLD_WINNERS_CURSE_UNADJUSTED_MULTIWINDOW_OR_LOT_OVERLAP')
        pilot_block_ids_set = {str(x) for x in pilot_blk_raw}
        t_crit = student_t_ppf(0.95, float(n_lots - 1))
    elif win_design == 'MULTI_WINDOW_CONFIRMATION_BONFERRONI':
        if k_cand < 2:
            return _hold('HOLD_INVALID_BONFERRONI_SINGLE_WINDOW_CONFIG')
        t_crit = student_t_ppf(1.0 - 0.05 / float(k_cand), float(n_lots - 1))
    else:
        return _hold('HOLD_WINNERS_CURSE_UNADJUSTED_MULTIWINDOW_OR_LOT_OVERLAP')

    # 4. Registered hormone unit, common-denominator recovery basis, and fluidics mode (Sol Fix #4):
    hormone_unit = _strict_unpadded_str_iip0029(record.get('hormone_unit'), 3)
    norm_denom = _strict_unpadded_str_iip0029(record.get('normalization_denominator'), 3)
    recov_basis = _strict_unpadded_str_iip0029(record.get('recovery_denominator_basis'), 3)
    if hormone_unit is None or hormone_unit not in C['REGISTERED_HORMONE_UNIT_SPECS'] or norm_denom is None or recov_basis is None:
        return _hold('HOLD_MISSING_OR_UNREGISTERED_HORMONE_UNIT_OR_DENOMINATOR')
    unit_spec = C['REGISTERED_HORMONE_UNIT_SPECS'][hormone_unit]
    if norm_denom != str(unit_spec['required_denominator']) or recov_basis != norm_denom:
        return _hold('HOLD_RECOVERY_DENOMINATOR_BASIS_MISMATCH')
    if record.get('eta_sort_includes_sorting_and_reaggregation_losses') is not True:
        return _hold('HOLD_ETA_SORT_MUST_INCLUDE_SORTING_AND_REAGGREGATION_LOSSES')
    z_min_floor = float(unit_spec['z_min_floor'])

    fluidics = record.get('assay_fluidics_mode')
    if fluidics not in ('DYNAMIC_PERIFUSION', 'STATIC_BATCH_BENCHMARK'):
        return _hold('HOLD_INVALID_ASSAY_FLUIDICS_MODE')

    # 5. Source-bound per-lot records with independent differentiation_block_ids and matched 40 + 100 mmHg coverage:
    lot_obs_list = record.get('lot_observations')
    if not isinstance(lot_obs_list, (list, tuple)) or len(lot_obs_list) != n_lots:
        return _hold('HOLD_MISSING_SOURCE_BOUND_PER_LOT_OBSERVATIONS')

    seen_obs_lot_ids: list[str] = []
    seen_block_ids: list[str] = []
    seen_trace_ids: list[str] = []
    block_by_lot: dict[str, str] = {}

    eta_v: list[float] = []
    basal1_v: list[float] = []
    low1_v: list[float] = []
    high1_v: list[float] = []
    low2_v: list[float] = []
    high2_v: list[float] = []
    kcl_v: list[float] = []
    adult_ref_z2_v: list[float] = []
    stored_iu_v: list[float] = []
    proins_ratio_v: list[float] = []

    peak1_v: list[float] = []
    tpeak1_v: list[float] = []
    t50_v: list[float] = []
    t90_v: list[float] = []
    tau_ch_v: list[float] = []
    gcg_fold_v: list[float] = []
    glc_thresh_v: list[float] = []
    ucn3_v: list[float] = []
    counterreg_v: list[float] = []

    # Matched 100 mmHg sibling vectors for DYNAMIC_PERIFUSION:
    h1_100_v: list[float] = []
    l1_100_v: list[float] = []
    h2_100_v: list[float] = []
    l2_100_v: list[float] = []
    tau_100_v: list[float] = []

    for item in lot_obs_list:
        if not isinstance(item, (dict, MappingProxyType)):
            return _hold('HOLD_INVALID_PER_LOT_OBSERVATION_RECORD')
        lsk = _coerce_lot_scope_key_iip0029(item.get('lot_scope_key'))
        if lsk is None:
            return _hold('HOLD_PER_LOT_RECORD_MISSING_LOT_SCOPE_KEY')
        if (
            lsk.product_variant_id != scope_key.product_variant_id
            or lsk.construct_sha256 != scope_key.construct_sha256
            or lsk.donor_seedline_id != scope_key.donor_seedline_id
            or lsk.mcb_id != scope_key.mcb_id
        ):
            return _hold('HOLD_PER_LOT_SCOPE_KEY_PROVENANCE_MISMATCH')
        lid = _strict_unpadded_str_iip0029(lsk.lot_id, 2)
        bid = _strict_unpadded_str_iip0029(item.get('differentiation_block_id'), 2)
        tid = _strict_unpadded_str_iip0029(item.get('trace_id'), 3)
        po2 = _strict_finite_float_iip0029(item.get('assay_po2_mmhg'))
        if lid is None or bid is None or tid is None or po2 is None or po2 not in C['ALLOWED_ASSAY_PO2_MMHG']:
            return _hold('HOLD_PER_LOT_MISSING_BLOCK_TRACE_ID_OR_INVALID_ASSAY_PO2')
        seen_obs_lot_ids.append(lid)
        seen_block_ids.append(bid)
        seen_trace_ids.append(tid)
        block_by_lot[lid] = bid

        vals: dict[str, float] = {}
        for k in (
            'eta_sort_viable_recovery',
            'basal1_pre_pulse',
            'low1_mean',
            'high1_mean',
            'low2_mean',
            'high2_mean',
            'kcl_mean',
            'adult_islet_ref_z2',
            'stored_mature_insulin_iu_implanted',
            'proinsulin_to_cpeptide_content_molar_ratio',
        ):
            fv = _strict_finite_float_iip0029(item.get(k))
            if fv is None or fv < 0.0:
                return _hold(f'HOLD_MISSING_OR_INVALID_PER_LOT_{k.upper()}')
            vals[k] = fv
        eta_v.append(vals['eta_sort_viable_recovery'])
        basal1_v.append(vals['basal1_pre_pulse'])
        low1_v.append(vals['low1_mean'])
        high1_v.append(vals['high1_mean'])
        low2_v.append(vals['low2_mean'])
        high2_v.append(vals['high2_mean'])
        kcl_v.append(vals['kcl_mean'])
        adult_ref_z2_v.append(vals['adult_islet_ref_z2'])
        stored_iu_v.append(vals['stored_mature_insulin_iu_implanted'])
        proins_ratio_v.append(vals['proinsulin_to_cpeptide_content_molar_ratio'])

        if fluidics == 'DYNAMIC_PERIFUSION':
            # Primary trace must be boundary-matched 40 mmHg paired with a matched 100 mmHg sibling trace:
            if abs(po2 - 40.0) > 1e-6:
                return _hold('HOLD_DYNAMIC_PRIMARY_TRACE_MUST_BE_40MMHG_PAIRED_WITH_100MMHG_SIBLING')
            # Censored short low-interval handling (Sol Fix #4):
            low2_dur = _strict_finite_float_iip0029(item.get('low2_washout_duration_min'))
            if low2_dur is None or low2_dur < float(C['MIN_UNCENSORED_LOW2_DURATION_MIN']) or item.get('effluent_t90_shutoff_min') is None:
                return _hold('HOLD_CENSORED_SHORT_LOW_INTERVAL_SHUTOFF_UNRESOLVED')
            # Per-trace tracer validity, alignment, and individual chamber ceiling <= 0.50 min:
            if item.get('tracer_response_valid') is not True or item.get('tracer_time_aligned') is not True:
                return _hold('HOLD_PER_TRACE_TRACER_RESPONSE_OR_TIME_ALIGNMENT_UNVERIFIED')
            for k, dst in (
                ('peak1_instantaneous', peak1_v),
                ('time_to_peak1_min', tpeak1_v),
                ('effluent_t50_shutoff_min', t50_v),
                ('effluent_t90_shutoff_min', t90_v),
                ('inert_tracer_tau_chamber_min', tau_ch_v),
                ('glucagon_2mM_fold_increase', gcg_fold_v),
                ('glucose_threshold_mM', glc_thresh_v),
                ('ucn3_maturity_ratio', ucn3_v),
                ('alpha_counterregulatory_gain_ratio', counterreg_v),
            ):
                fv = _strict_finite_float_iip0029(item.get(k))
                if fv is None or fv <= 0.0:
                    return _hold('HOLD_MISSING_DYNAMIC_PERIFUSION_PEAK_SHUTOFF_TRACER_OR_GATE0M_METRIC')
                dst.append(fv)
            if tau_ch_v[-1] > float(C['TRACER_TAU_CHAMBER_MAX_MIN']):
                return _hold('HOLD_INDIVIDUAL_CHAMBER_DISPERSION_TAU_EXCEEDED_0P50_MIN')

            # Matched 100 mmHg sibling trace coverage on every independent block:
            tid_100 = _strict_unpadded_str_iip0029(item.get('matched_100mmhg_trace_id'), 3)
            po2_100 = _strict_finite_float_iip0029(item.get('matched_100mmhg_assay_po2_mmhg'))
            tau_100 = _strict_finite_float_iip0029(item.get('matched_100mmhg_tracer_tau_chamber_min'))
            h1_100 = _strict_finite_float_iip0029(item.get('matched_100mmhg_high1_mean'))
            l1_100 = _strict_finite_float_iip0029(item.get('matched_100mmhg_low1_mean'))
            h2_100 = _strict_finite_float_iip0029(item.get('matched_100mmhg_high2_mean'))
            l2_100 = _strict_finite_float_iip0029(item.get('matched_100mmhg_low2_mean'))
            if (
                tid_100 is None
                or po2_100 is None
                or abs(po2_100 - 100.0) > 1e-6
                or tau_100 is None
                or tau_100 <= 0.0
                or h1_100 is None or h1_100 < 0.0
                or l1_100 is None or l1_100 < 0.0
                or h2_100 is None or h2_100 < 0.0
                or l2_100 is None or l2_100 < 0.0
                or item.get('matched_100mmhg_tracer_response_valid') is not True
                or item.get('matched_100mmhg_tracer_time_aligned') is not True
            ):
                return _hold('HOLD_MISSING_MATCHED_100MMHG_OXYGEN_COVERAGE_OR_TRACER_ALIGNMENT')
            if tau_100 > float(C['TRACER_TAU_CHAMBER_MAX_MIN']):
                return _hold('HOLD_INDIVIDUAL_CHAMBER_DISPERSION_TAU_EXCEEDED_0P50_MIN')
            seen_trace_ids.append(tid_100)
            h1_100_v.append(h1_100)
            l1_100_v.append(l1_100)
            h2_100_v.append(h2_100)
            l2_100_v.append(l2_100)
            tau_100_v.append(tau_100)

    expected_traces = (2 * n_lots) if fluidics == 'DYNAMIC_PERIFUSION' else n_lots
    if seen_obs_lot_ids != lot_ids or len(set(seen_trace_ids)) != expected_traces:
        return _hold('HOLD_PER_LOT_IDS_OR_TRACE_IDS_MISMATCH_CONFIRMATION_SET')
    if len(set(seen_block_ids)) != n_lots or (pilot_block_ids_set & set(seen_block_ids)):
        return _hold('HOLD_NON_INDEPENDENT_DIFFERENTIATION_BLOCKS_OR_PILOT_BLOCK_OVERLAP')

    # 6. Locked UTC Timestamps, Dual Clocks & Process Mode (Sol Fix #4):
    proc_mode = record.get('preparation_process_mode')
    sorting_performed = record.get('sorting_performed')
    if proc_mode not in ('SORTED_REAGGREGATED', 'UNSORTED_CONTROL') or type(sorting_performed) is not bool:
        return _hold('HOLD_MISSING_OR_INVALID_PROCESS_SORTING_MODE')
    if (proc_mode == 'SORTED_REAGGREGATED' and sorting_performed is not True) or (
        proc_mode == 'UNSORTED_CONTROL' and sorting_performed is not False
    ):
        return _hold('HOLD_SORTING_MODE_BOOLEAN_CONTRADICTION')

    dt_s7_start = _parse_utc_iso_iip0029(record.get('stage7_start_utc'))
    dt_sort = _parse_utc_iso_iip0029(record.get('sort_reaggregation_utc'))
    dt_harvest = _parse_utc_iso_iip0029(record.get('assay_harvest_utc'))
    if dt_s7_start is None or dt_sort is None or dt_harvest is None:
        return _hold('HOLD_MISSING_OR_INVALID_LOCKED_CLOCK_TIMESTAMPS')
    if not (dt_s7_start <= dt_sort <= dt_harvest):
        return _hold('HOLD_NON_MONOTONIC_CLOCK_TIMESTAMPS')

    s7_total = _strict_finite_float_iip0029(record.get('stage7_total_days'))
    s7_sort = _strict_finite_float_iip0029(record.get('sort_reaggregation_stage7_day'))
    post_sort = _strict_finite_float_iip0029(record.get('days_post_sort_reaggregation'))
    if s7_total is None or s7_sort is None or post_sort is None:
        return _hold('HOLD_MISSING_OR_INVALID_TWOCLOCK_METADATA')
    if not (float(C['STAGE7_MIN_TOTAL_DAYS']) <= s7_total <= float(C['STAGE7_MAX_TOTAL_DAYS'])):
        return _hold('HOLD_STAGE7_TOTAL_DAYS_OUT_OF_QUALIFIED_RANGE')
    if abs((s7_total - s7_sort) - post_sort) > 1e-6:
        return _hold('HOLD_TWOCLOCK_ARITHMETIC_MISMATCH')

    ts_total_d = (dt_harvest - dt_s7_start).total_seconds() / 86400.0
    ts_sort_d = (dt_sort - dt_s7_start).total_seconds() / 86400.0
    ts_post_d = (dt_harvest - dt_sort).total_seconds() / 86400.0
    tol_d = float(C['TIMESTAMP_DAY_TOLERANCE'])
    if abs(ts_total_d - s7_total) > tol_d or abs(ts_sort_d - s7_sort) > tol_d or abs(ts_post_d - post_sort) > tol_d:
        return _hold('HOLD_TWOCLOCK_TIMESTAMP_DAY_MISMATCH')

    if proc_mode == 'SORTED_REAGGREGATED':
        if s7_sort < float(C['SORTED_MIN_SORT_DAY']) or post_sort < float(C['SORTED_MIN_POST_SORT_DAYS']):
            return _hold('HOLD_SORTED_CLOCK_BELOW_MINIMUM_MATURATION_DAYS')
        if any(e < float(C['SORTED_ETA_MIN']) or e >= float(C['SORTED_ETA_MAX_EXCLUSIVE']) for e in eta_v):
            return _hold('HOLD_IMPOSSIBLE_OR_UNQUALIFIED_SORTED_ETA_RECOVERY')
    else:
        if abs(post_sort) > 1e-9 or abs(s7_sort - s7_total) > 1e-9 or any(abs(e - 1.0) > 1e-9 for e in eta_v):
            return _hold('HOLD_UNSORTED_CONTROL_CLOCK_OR_ETA_MISMATCH')

    # 7. Stage-7 Maturation Pharmacology (ALK5i-free + T3/NAC + ZM447439 tolerance + LLOQ washout + vascular function):
    if record.get('obs_alk5i_free_maturation_applied') is not True or record.get('terminal_alk5i_tgfb_inhibitor_continued') is not False:
        return _hold('HOLD_PROHIBITED_OR_UNVERIFIED_ALK5I_FREE_MATURATION')
    if record.get('stage7_t3_10nM_and_nac_1mM_retained') is not True:
        return _hold('HOLD_STAGE7_MISSING_REQUIRED_10NM_T3_OR_1MM_NAC')

    zm_regimen = record.get('post_sort_endocrine_zm447439_regimen')
    zm_uM = _strict_finite_float_iip0029(record.get('post_sort_zm447439_uM'))
    if zm_regimen not in ('S7W0_ONLY_0UM_POST_SORT', 'CONTINUOUS_0P5UM_BALBOA_S7') or zm_uM is None:
        return _hold('HOLD_MISSING_OR_UNREGISTERED_ZM447439_COVARIATE')
    if zm_regimen == 'S7W0_ONLY_0UM_POST_SORT' and abs(zm_uM) > 1e-9:
        return _hold('HOLD_ZM447439_CONCENTRATION_REGIMEN_MISMATCH')
    if zm_regimen == 'CONTINUOUS_0P5UM_BALBOA_S7' and not (
        float(C['CONTINUOUS_ZM447439_MIN_UM']) <= zm_uM <= float(C['CONTINUOUS_ZM447439_MAX_UM'])
    ):
        return _hold('HOLD_ZM447439_CONCENTRATION_REGIMEN_MISMATCH')

    iec_coassembled = record.get('vascular_iec_coassembled')
    if type(iec_coassembled) is not bool:
        return _hold('HOLD_MISSING_VASCULAR_IEC_COASSEMBLY_DECLARATION')
    if iec_coassembled is True:
        coassembly_zm = _strict_finite_float_iip0029(record.get('coassembly_zm447439_uM'))
        coassembly_lloq = _strict_finite_float_iip0029(record.get('coassembly_zm447439_lloq_uM'))
        wash_method = _strict_unpadded_str_iip0029(record.get('zm447439_washout_assay_method'), 3)
        wash_cycles = _strict_pos_int_iip0029(record.get('zm447439_wash_cycles'))
        spike_rec = _strict_finite_float_iip0029(record.get('zm447439_matrix_spike_recovery_pct'))
        if (
            coassembly_zm is None
            or coassembly_lloq is None
            or wash_method is None
            or wash_cycles is None
            or spike_rec is None
            or wash_cycles < int(C['PRE_IEC_MIN_WASH_CYCLES'])
            or not (0.0 < coassembly_lloq <= float(C['PRE_IEC_ZM447439_MAX_LLOQ_UM']))
            or not (0.0 <= coassembly_zm <= coassembly_lloq)
            or not (float(C['PRE_IEC_SPIKE_RECOVERY_MIN_PCT']) <= spike_rec <= float(C['PRE_IEC_SPIKE_RECOVERY_MAX_PCT']))
            or record.get('zm447439_washout_verified_pre_iec') is not True
            or record.get('pre_coassembly_iec_vascular_function_verified') is not True
        ):
            return _hold('HOLD_UNWASHED_OR_UNVERIFIED_LLOQ_ZM447439_OR_VASCULAR_FUNCTION_PRE_IEC')

    # 8. Integrated SI_1, KCl Reserve, Separate Readout 2 vs Gate 0M Basal/Peak/Shutoff & Full Gate 0M Paracrine (Sol Fix #2):
    _, low1_l95, _ = _one_sided_t_bounds_iip0029(low1_v, t_crit)
    _, si1_l95, _ = _one_sided_t_bounds_iip0029([h - float(C['INTEGRATED_SI1_MIN']) * lo for h, lo in zip(high1_v, low1_v)], t_crit)
    _, kcl_vs_h2_l95, _ = _one_sided_t_bounds_iip0029([k - h2 for k, h2 in zip(kcl_v, high2_v)], t_crit)
    _, kcl_vs_l2_l95, _ = _one_sided_t_bounds_iip0029([k - l2 for k, l2 in zip(kcl_v, low2_v)], t_crit)
    if low1_l95 <= 0.0 or si1_l95 < 0.0 or kcl_vs_h2_l95 <= 0.0 or kcl_vs_l2_l95 <= 0.0:
        return _hold('HOLD_INTEGRATED_SI1_OR_KCL_RESERVE_FAILED')

    basal_peak_stats: dict[str, Any] | None = None
    if fluidics == 'DYNAMIC_PERIFUSION':
        _, _, tau_ch_u95 = _one_sided_t_bounds_iip0029(tau_ch_v, t_crit)
        _, _, tau_100_u95 = _one_sided_t_bounds_iip0029(tau_100_v, t_crit)
        if tau_ch_u95 > float(C['TRACER_TAU_CHAMBER_MAX_MIN']) or tau_100_u95 > float(C['TRACER_TAU_CHAMBER_MAX_MIN']):
            return _hold('HOLD_CHAMBER_DISPERSION_TAU_EXCEEDED_REQUIRES_FAST_CHIP')
        if any(a >= b for a, b in zip(t50_v, t90_v)):
            return _hold('HOLD_SHUTOFF_T50_GE_T90_MONOTONICITY_VIOLATION')

        _, _, t50_u95 = _one_sided_t_bounds_iip0029(t50_v, t_crit)
        if t50_u95 > float(C['READOUT2_T50_MAX_MIN']):
            return _hold('HOLD_READOUT2_STEPDOWN_T50_U95_EXCEEDED')
        _, _, t90_u95 = _one_sided_t_bounds_iip0029(t90_v, t_crit)
        if t90_u95 > float(C['GATE0M_T90_MAX_MIN']):
            return _hold('HOLD_GATE0M_S8_STEPDOWN_T90_U95_EXCEEDED')

        # Canonical Readout 2 (Peak1/Basal1 >= 4.5 within <= 8.0 min) explicitly distinguished from
        # stricter Gate 0M / extended-maturation low-basal check (Basal1/Peak1 <= 0.15 => Peak1/Basal1 >= 6.67x):
        _, p1_4p5_l95, _ = _one_sided_t_bounds_iip0029([p - float(C['DYNAMIC_PEAK1_OVER_BASAL1_MIN']) * b for p, b in zip(peak1_v, basal1_v)], t_crit)
        _, _, tpeak1_u95 = _one_sided_t_bounds_iip0029(tpeak1_v, t_crit)
        if p1_4p5_l95 < 0.0 or tpeak1_u95 > float(C['DYNAMIC_TIME_TO_PEAK1_MAX_MIN']):
            return _hold('HOLD_READOUT2_FIRST_PHASE_PEAK_OR_LATENCY_FAILED')

        _, g0m_b1_l95, _ = _one_sided_t_bounds_iip0029([float(C['GATE0M_BASAL1_FRAC_PEAK1_MAX']) * p - b for p, b in zip(peak1_v, basal1_v)], t_crit)
        _, g0m_l2_l95, _ = _one_sided_t_bounds_iip0029([float(C['DYNAMIC_LOW2_FRAC_PEAK1_MAX']) * p - l2 for p, l2 in zip(peak1_v, low2_v)], t_crit)
        _, wash_l95, _ = _one_sided_t_bounds_iip0029([float(C['LOW2_OVER_LOW1_WASHOUT_MAX']) * l1 - l2 for l1, l2 in zip(low1_v, low2_v)], t_crit)
        if g0m_b1_l95 < 0.0 or g0m_l2_l95 < 0.0 or wash_l95 < 0.0:
            return _hold('HOLD_GATE0M_BASAL_FRAC_PEAK_OR_LOW2_WASHOUT_FAILED')

        # Full canonical Gate 0M paracrine checks (glucagon, threshold, UCN3, counterregulation):
        _, gcg_l95, _ = _one_sided_t_bounds_iip0029(gcg_fold_v, t_crit)
        _, thresh_l95, _ = _one_sided_t_bounds_iip0029(glc_thresh_v, t_crit)
        _, ucn3_l95, _ = _one_sided_t_bounds_iip0029(ucn3_v, t_crit)
        _, creg_l95, _ = _one_sided_t_bounds_iip0029(counterreg_v, t_crit)
        if (
            gcg_l95 < float(C['GATE0M_GLUCAGON_FOLD_MIN'])
            or thresh_l95 < float(C['GATE0M_GLUCOSE_THRESHOLD_MIN_MM'])
            or ucn3_l95 < float(C['GATE0M_UCN3_RATIO_MIN'])
            or creg_l95 < float(C['GATE0M_COUNTERREG_GAIN_MIN'])
        ):
            return _hold('HOLD_GATE0M_PARACRINE_GLUCAGON_THRESHOLD_UCN3_OR_COUNTERREG_FAILED')

        # Matched 100 mmHg sibling SI_1 and Z2 positivity verification:
        _, si1_100_l95, _ = _one_sided_t_bounds_iip0029([h - float(C['INTEGRATED_SI1_MIN']) * lo for h, lo in zip(h1_100_v, l1_100_v)], t_crit)
        z2_100_vec = [e * (h2 - l2) for e, h2, l2 in zip(eta_v, h2_100_v, l2_100_v)]
        _, z2_100_l95, _ = _one_sided_t_bounds_iip0029(z2_100_vec, t_crit)
        if si1_100_l95 < 0.0 or z2_100_l95 < z_min_floor:
            return _hold('HOLD_MATCHED_100MMHG_SIBLING_SI1_OR_PULSE2_INCREMENT_FAILED')

        b1_p1_ratios = [b / max(1e-12, p) for b, p in zip(basal1_v, peak1_v)]
        p1_b1_ratios = [p / max(1e-12, b) for b, p in zip(basal1_v, peak1_v)]
        basal_peak_stats = {
            'basal1_pre_pulse_window_min': [-10.0, 0.0],
            'peak1_window_min': [0.0, 8.0],
            'readout2_nominal_peak1_over_basal1_min': float(C['DYNAMIC_PEAK1_OVER_BASAL1_MIN']),
            'readout2_peak1_over_basal1_mean': round(sum(p1_b1_ratios) / float(n_lots), 4),
            'readout2_peak1_minus_4p5_basal1_l95': round(p1_4p5_l95, 4),
            'gate0m_strict_basal1_over_peak1_max': float(C['GATE0M_BASAL1_FRAC_PEAK1_MAX']),
            'gate0m_basal1_over_peak1_mean': round(sum(b1_p1_ratios) / float(n_lots), 4),
            'gate0m_0p15_peak1_minus_basal1_l95': round(g0m_b1_l95, 4),
            'relationship_note': 'Gate 0M Basal1/Peak1 <= 0.15 (Peak1/Basal1 >= 6.67x) is stricter than nominal Readout 2 Peak1/Basal1 >= 4.5x on the same pre-step [-10, 0] min Basal1 and [0, 8] min Peak1 windows.',
        }
    else:
        _, l1_kcl_l95, _ = _one_sided_t_bounds_iip0029([float(C['STATIC_LOW_OVER_KCL_MAX']) * k - l1 for k, l1 in zip(kcl_v, low1_v)], t_crit)
        _, l2_kcl_l95, _ = _one_sided_t_bounds_iip0029([float(C['STATIC_LOW_OVER_KCL_MAX']) * k - l2 for k, l2 in zip(kcl_v, low2_v)], t_crit)
        if l1_kcl_l95 < 0.0 or l2_kcl_l95 < 0.0:
            return _hold('HOLD_STATIC_BASAL_LEAK_OVER_KCL_EXCEEDED')

    # 9. Lot-Paired Pre-Sort-Normalized Net Increments Z1, Z2 + Contemporaneous Adult-Islet Reference Floor (Sol Fix #1):
    if 'within_lot_technical_var_z2' not in record:
        return _hold('HOLD_MISSING_EXPLICIT_WITHIN_LOT_TECHNICAL_VARIANCE')
    tech_var_z2 = _strict_finite_float_iip0029(record.get('within_lot_technical_var_z2'))
    if tech_var_z2 is None or tech_var_z2 < 0.0:
        return _hold('HOLD_INVALID_WITHIN_LOT_TECHNICAL_VARIANCE')

    z1_vec = [e * (h1 - l1) for e, h1, l1 in zip(eta_v, high1_v, low1_v)]
    z2_vec = [e * (h2 - l2) for e, h2, l2 in zip(eta_v, high2_v, low2_v)]
    z1_mean, z1_l95, _ = _one_sided_t_bounds_iip0029(z1_vec, t_crit)
    z2_mean, z2_l95, _ = _one_sided_t_bounds_iip0029(z2_vec, t_crit, extra_within_var=tech_var_z2 / float(n_lots))
    if z1_l95 < z_min_floor or z2_l95 < z_min_floor:
        return _hold('HOLD_NET_INCREMENT_PER_PRESORT_CELL_L95_INSUFFICIENT')

    # Contemporaneous same-assay adult-islet reference floor (blocks tiny 0.001 increments clearing on ratio alone!):
    ref_mean, ref_z2_l95, _ = _one_sided_t_bounds_iip0029(adult_ref_z2_v, t_crit)
    if ref_z2_l95 < z_min_floor:
        return _hold('HOLD_CONTEMPORANEOUS_ADULT_ISLET_REFERENCE_PULSE2_INSUFFICIENT')
    ref_frac = float(C['ADULT_ISLET_REF_PULSE2_MIN_FRACTION'])
    _, z2_vs_ref_l95, _ = _one_sided_t_bounds_iip0029([z2 - ref_frac * rv for z2, rv in zip(z2_vec, adult_ref_z2_v)], t_crit)
    if z2_vs_ref_l95 < 0.0:
        return _hold(
            'HOLD_SECOND_PULSE_BELOW_ADULT_ISLET_REFERENCE_FLOOR',
            z2_l95=round(z2_l95, 4),
            adult_ref_z2_mean=round(ref_mean, 4),
            z2_minus_half_adult_ref_l95=round(z2_vs_ref_l95, 4),
        )

    # 10. Clock-Classified Extended-Maturation Exposure -> Paired Within-Recipe Rule S8 Fieller (theta_0 = 0.85):
    is_extended_s8 = (s7_total > float(C['BASELINE_WINDOW_MAX_STAGE7_DAYS'])) or (
        post_sort > float(C['BASELINE_WINDOW_MAX_POST_SORT_DAYS'])
    )
    s8_fieller_l95: float | None = None
    if is_extended_s8:
        early_map = record.get('paired_earlier_window_z2_by_lot')
        if not isinstance(early_map, (dict, MappingProxyType)) or set(early_map.keys()) != set(lot_ids):
            return _hold('HOLD_S8_EXTENDED_CLOCK_REQUIRES_SOURCE_LINKED_EARLIER_WINDOW')
        early_vals: list[float] = []
        early_100_vals: list[float] = []
        early_traces_seen: set[str] = set()
        for lid in lot_ids:
            raw_e = early_map[lid]
            if isinstance(raw_e, (dict, MappingProxyType)):
                ebid = _strict_unpadded_str_iip0029(raw_e.get('differentiation_block_id'), 2)
                et40 = _strict_unpadded_str_iip0029(raw_e.get('trace_id_40mmhg'), 3)
                ez40 = _strict_finite_float_iip0029(raw_e.get('z2_40mmhg'))
                dt_eh = _parse_utc_iso_iip0029(raw_e.get('early_harvest_utc'))
                if (
                    ebid != block_by_lot[lid]
                    or et40 is None
                    or ez40 is None
                    or dt_eh is None
                    or dt_eh >= dt_harvest
                    or et40 in seen_trace_ids
                    or et40 in early_traces_seen
                ):
                    return _hold('HOLD_S8_EXTENDED_CLOCK_EARLIER_WINDOW_BLOCK_OR_TIMESTAMP_MISMATCH')
                early_traces_seen.add(et40)
                early_vals.append(ez40)
                if fluidics == 'DYNAMIC_PERIFUSION':
                    et100 = _strict_unpadded_str_iip0029(raw_e.get('trace_id_100mmhg'), 3)
                    ez100 = _strict_finite_float_iip0029(raw_e.get('z2_100mmhg'))
                    if et100 is None or ez100 is None or et100 in seen_trace_ids or et100 in early_traces_seen:
                        return _hold('HOLD_S8_EXTENDED_CLOCK_MISSING_MATCHED_100MMHG_EARLIER_ARM')
                    early_traces_seen.add(et100)
                    early_100_vals.append(ez100)
            else:
                return _hold('HOLD_S8_EXTENDED_CLOCK_REQUIRES_BLOCK_BOUND_EARLIER_RECORD')
        _, early_z2_l95, _ = _one_sided_t_bounds_iip0029(early_vals, t_crit)
        if early_z2_l95 < z_min_floor:
            return _hold('HOLD_S8_EARLIER_REFERENCE_WINDOW_FAILED_POSITIVITY_FLOOR')
        theta0 = float(C['S8_WITHIN_RECIPE_FIELLER_THETA0'])
        _, s8_fieller_l95, _ = _one_sided_t_bounds_iip0029([z_e - theta0 * z_0 for z_e, z_0 in zip(z2_vec, early_vals)], t_crit)
        if s8_fieller_l95 < 0.0:
            return _hold(
                'HOLD_S8_EXTENDED_MATURATION_FIELLER_DEGRADED_KEEP_EARLIER_WINDOW',
                s8_fieller_diff_l95=round(s8_fieller_l95, 4),
                theta0_applied=theta0,
            )
        if fluidics == 'DYNAMIC_PERIFUSION':
            _, early_100_l95, _ = _one_sided_t_bounds_iip0029(early_100_vals, t_crit)
            z2_100_v = [e * (h2 - l2) for e, h2, l2 in zip(eta_v, h2_100_v, l2_100_v)]
            _, s8_100_fieller_l95, _ = _one_sided_t_bounds_iip0029([z_e - theta0 * z_0 for z_e, z_0 in zip(z2_100_v, early_100_vals)], t_crit)
            if early_100_l95 < z_min_floor or s8_100_fieller_l95 < 0.0:
                return _hold(
                    'HOLD_S8_EXTENDED_MATURATION_100MMHG_FIELLER_DEGRADED',
                    s8_100mmhg_fieller_diff_l95=round(s8_100_fieller_l95, 4),
                )

    # 11. Specific Mature-Insulin Inventory Conservation, Exploratory Proinsulin:C-Peptide, & >138.3 IU Gate 0S Trigger:
    if record.get('insulin_inventory_assay_specificity') != 'MATURE_INSULIN_SPECIFIC_LOW_PROINSULIN_CROSSREACTIVITY':
        return _hold('HOLD_CROSS_REACTIVE_OR_UNVERIFIED_MATURE_INSULIN_INVENTORY_ASSAY')
    if 'assay_cumulative_secreted_iu_equiv' not in record:
        return _hold('HOLD_MISSING_EXPLICIT_ASSAY_CUMULATIVE_SECRETED_IU')
    cum_assay_iu = _strict_finite_float_iip0029(record.get('assay_cumulative_secreted_iu_equiv'))
    if cum_assay_iu is None or cum_assay_iu < 0.0:
        return _hold('HOLD_INVALID_ASSAY_CUMULATIVE_SECRETED_IU')

    beta_scale = (beta_frac / 0.60) * (imp_ieq / 500_000.0)
    min_phys_store_iu = max(float(C['STORED_INSULIN_FLOOR_AT_300K_BETA_IU']) * beta_scale, cum_assay_iu)
    _, stored_l95, stored_u95 = _one_sided_t_bounds_iip0029(stored_iu_v, t_crit)
    if stored_l95 < min_phys_store_iu:
        return _hold('HOLD_STORED_INSULIN_L95_BELOW_PHYSIOLOGICAL_OR_MASS_FLOOR')

    _, _, proins_u95 = _one_sided_t_bounds_iip0029(proins_ratio_v, t_crit)
    if proins_u95 > float(C['PROINSULIN_TO_CPEPTIDE_CONTENT_MOLAR_MAX']):
        return _hold(
            'HOLD_ELEVATED_PROINSULIN_TO_CPEPTIDE_CONTENT_RATIO',
            proinsulin_to_cpeptide_u95=round(proins_u95, 4),
            proinsulin_to_cpeptide_marker_classification='EXPLORATORY_STAGE0_IN_VITRO_PROCESSING_MARKER_PENDING_CALIBRATION',
        )

    max_lot_store_iu = max(stored_iu_v)
    ceiling_iu = float(C['GATE0S_ABSOLUTE_STORED_INSULIN_CEILING_IU'])
    if stored_u95 > ceiling_iu or max_lot_store_iu > ceiling_iu:
        return _hold(
            'HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW',
            gate0s_escalation_required=True,
            all_non_store_ceiling_checks_passed=True,
            stored_insulin_l95_iu=round(stored_l95, 4),
            stored_insulin_u95_iu=round(stored_u95, 4),
            stored_insulin_max_lot_iu=round(max_lot_store_iu, 4),
            registry_execution_supported=registry_exec_supported,
            proinsulin_to_cpeptide_u95=round(proins_u95, 4),
            proinsulin_to_cpeptide_marker_classification='EXPLORATORY_STAGE0_IN_VITRO_PROCESSING_MARKER_PENDING_CALIBRATION',
        )

    common_payload: dict[str, Any] = {
        'subgate_pass': True,
        'release_authority': 'DELEGATED_TO_LIVE_ATOMIC_STAGE0_LEDGER',
        'gate0s_escalation_required': False,
        'all_non_store_ceiling_checks_passed': True,
        'registry_execution_supported': registry_exec_supported,
        'is_extended_s8_exposure': is_extended_s8,
        'registered_mixture_id': reg_mix_id,
        't_crit_applied': round(t_crit, 6),
        'z1_mean': round(z1_mean, 4),
        'z1_l95': round(z1_l95, 4),
        'z2_mean': round(z2_mean, 4),
        'z2_l95': round(z2_l95, 4),
        'adult_islet_ref_z2_l95': round(ref_z2_l95, 4),
        'z2_minus_half_adult_ref_l95': round(z2_vs_ref_l95, 4),
        'stored_insulin_l95_iu': round(stored_l95, 4),
        'stored_insulin_u95_iu': round(stored_u95, 4),
        'stored_insulin_max_lot_iu': round(max_lot_store_iu, 4),
        'proinsulin_to_cpeptide_u95': round(proins_u95, 4),
        'proinsulin_to_cpeptide_marker_classification': 'EXPLORATORY_STAGE0_IN_VITRO_PROCESSING_MARKER_PENDING_CALIBRATION',
        'proinsulin_provenance_note': (
            'Keymeulen et al. 2024 (PMID:38012450) measured post-hoc plasma proinsulin:C-peptide in a selected responder '
            'and does not validate 0.10 as a universal in-vitro cutoff; retained strictly as an exploratory Stage-0 processing marker.'
        ),
        's8_fieller_diff_l95': None if s8_fieller_l95 is None else round(s8_fieller_l95, 4),
        'post_sort_endocrine_zm447439_regimen': zm_regimen,
        'basal_peak_statistics_and_windows': basal_peak_stats,
        'diagnostic_high2_over_high1_mean_nonblocking': round(sum(h2 / max(1e-12, h1) for h2, h1 in zip(high2_v, high1_v)) / float(n_lots), 4),
    }
    if fluidics == 'STATIC_BATCH_BENCHMARK' or proc_mode == 'UNSORTED_CONTROL':
        common_payload['status'] = 'COMPARATIVE_BENCHMARK_QUALIFIED_NON_RELEASE'
        return common_payload
    if not registry_exec_supported:
        common_payload['status'] = 'RESEARCH_COMPARATOR_QUALIFIED_NON_EXECUTABLE'
        return common_payload
    common_payload['status'] = 'QUALIFIED_TWOCLOCK_REPEATED_CHALLENGE_SUBGATE_EVIDENCE_ONLY'
    return common_payload


def evaluate_stage0_twoclock_flat_or_structured_iip_0029(
    observations: dict[str, Any],
    *,
    validated_stage0: dict[str, float],
    se_provenance_by_mean: dict[str, str],
    upper95_fn: Any,
    lower95_fn: Any,
    root: Path | str | None = None,
    product_def: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Unified Stage-0 IIP-0029 evaluator wired identically into evaluate_stage0_decision_tree,
    evaluate_stage0_research_and_molecular_guardrails, ingest_stage0_research_or_wetlab_batch,
    and replay_stage0_chronological_batches (DEC-108).
    Enforces identical fail-closed handling for missing records, ALK5i/washout violations,
    wrong-lot scope mismatches, and non-overridable >138.3 IU Gate 0S escalations.
    """
    C = STAGE0_IIP0029_CONSTANTS

    def _res(status: str, gate_pass: bool, **extra: Any) -> dict[str, Any]:
        out: dict[str, Any] = {
            'gate_pass': gate_pass,
            'status': status,
            'release_authority': 'DELEGATED_TO_LIVE_ATOMIC_STAGE0_LEDGER',
            'decision_id': 'DEC-108',
            'contribution_id': 'IIP-0029',
        }
        out.update(extra)
        return out

    def _check_scope_matched_gate0s_escalation() -> tuple[bool, str]:
        req_g0s = (
            'obs_gate0s_ap1903_1h_insulin_burst_iu',
            'obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h',
            'obs_explant_en_bloc_recovery_fraction',
        )
        if not all(k in validated_stage0 and se_provenance_by_mean.get(k) == 'EXPLICIT_SE_SUPPLIED' for k in req_g0s):
            return False, 'HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW'
        outer_lot = observations.get('current_lot_id', observations.get('lot_id'))
        outer_dev = observations.get('device_assembly_id')
        matched_tag = observations.get('gate0s_elevated_store_scope_matched_device_lot_id')
        if outer_lot is not None or outer_dev is not None:
            if not isinstance(matched_tag, str) or not matched_tag.strip():
                return False, 'HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW'
            if outer_lot is not None and str(outer_lot) not in matched_tag:
                return False, 'FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH'
            if outer_dev is not None and str(outer_dev) not in matched_tag:
                return False, 'FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH'
        if upper95_fn('obs_gate0s_ap1903_1h_insulin_burst_iu') > 2.8:
            return False, 'FAIL_GATE0S_ELEVATED_STORE_BURST_OR_DUMP_BREACHED'
        if upper95_fn('obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h') > 0.15:
            return False, 'FAIL_GATE0S_ELEVATED_STORE_BURST_OR_DUMP_BREACHED'
        if lower95_fn('obs_explant_en_bloc_recovery_fraction') < 0.99:
            return False, 'FAIL_GATE0S_ELEVATED_STORE_BURST_OR_DUMP_BREACHED'
        return True, 'CLEARED_ESCALATED_SCOPE_MATCHED_GATE0S_REVIEW'

    structured = observations.get('iip0029_twoclock_record')
    if structured is not None:
        if not isinstance(structured, (dict, MappingProxyType)):
            return _res('HOLD_INVALID_RECORD_TYPE', False)
        # Verify lot & product provenance match outer payload scope (rejects wrong-lot borrowing!):
        rec_sk = _coerce_lot_scope_key_iip0029(structured.get('lot_scope_key'))
        if rec_sk is not None:
            outer_pvid = product_def.get('product_variant_id') if isinstance(product_def, dict) else observations.get('product_variant_id')
            outer_lot = observations.get('current_lot_id', observations.get('lot_id'))
            outer_csha = observations.get('construct_sha256')
            if outer_pvid is not None and rec_sk.product_variant_id != str(outer_pvid):
                return _res('FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH', False)
            if outer_lot is not None and rec_sk.lot_id != str(outer_lot):
                return _res('FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH', False)
            if outer_csha is not None and rec_sk.construct_sha256 != str(outer_csha):
                return _res('FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH', False)

        sub = evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(structured, root=root)
        if sub.get('gate0s_escalation_required') is True and sub.get('all_non_store_ceiling_checks_passed') is True:
            g0s_ok, g0s_status = _check_scope_matched_gate0s_escalation()
            if not g0s_ok:
                return _res(g0s_status, False, subgate=sub, gate0s_escalation_required=True)
            return _res(
                'CLEARED_IIP0029_WITH_ESCALATED_SCOPE_MATCHED_GATE0S_REVIEW',
                True,
                subgate=sub,
                gate0s_escalation_required=True,
                gate0s_escalation_status=g0s_status,
            )
        return _res(
            str(sub['status']),
            bool(sub['subgate_pass']),
            subgate=sub,
            gate0s_escalation_required=bool(sub.get('gate0s_escalation_required', False)),
        )

    # Flat endpoint evaluation path (identical fail-closed rules when structured record is omitted):
    if observations.get('obs_alk5i_free_maturation_applied') is not True or observations.get('terminal_alk5i_tgfb_inhibitor_continued') is not False:
        return _res('HOLD_PROHIBITED_OR_UNVERIFIED_ALK5I_FREE_MATURATION', False)
    if observations.get('stage7_t3_10nM_and_nac_1mM_retained') is not True:
        return _res('HOLD_STAGE7_MISSING_REQUIRED_10NM_T3_OR_1MM_NAC', False)

    zm_reg = observations.get('post_sort_endocrine_zm447439_regimen')
    if zm_reg not in ('S7W0_ONLY_0UM_POST_SORT', 'CONTINUOUS_0P5UM_BALBOA_S7'):
        return _res('HOLD_MISSING_OR_UNREGISTERED_ZM447439_COVARIATE', False)
    if (
        observations.get('zm447439_washout_verified_pre_iec') is not True
        or observations.get('pre_coassembly_iec_vascular_function_verified') is not True
        or 'obs_stage7_pre_iec_zm447439_residual_uM' not in validated_stage0
        or upper95_fn('obs_stage7_pre_iec_zm447439_residual_uM') > float(C['PRE_IEC_ZM447439_MAX_LLOQ_UM'])
    ):
        return _res('HOLD_UNWASHED_OR_UNVERIFIED_LLOQ_ZM447439_OR_VASCULAR_FUNCTION_PRE_IEC', False)

    if observations.get('readout2_both_40_and_100mmhg_tracer_aligned_verified') is not True:
        return _res('HOLD_MISSING_MATCHED_100MMHG_OXYGEN_COVERAGE_OR_TRACER_ALIGNMENT', False)
    if _strict_unpadded_str_iip0029(observations.get('locked_analysis_plan_id'), 4) is None:
        return _res('HOLD_MISSING_LOCKED_ANALYSIS_PLAN_ID', False)
    reg_mix = observations.get('registered_mixture_id')
    if not isinstance(reg_mix, str) or reg_mix not in C['REGISTERED_MIXTURES']:
        return _res('HOLD_MISSING_OR_UNREGISTERED_MIXTURE_ID', False)

    dt_s = _parse_utc_iso_iip0029(observations.get('stage7_start_utc'))
    dt_r = _parse_utc_iso_iip0029(observations.get('sort_reaggregation_utc'))
    dt_h = _parse_utc_iso_iip0029(observations.get('assay_harvest_utc'))
    if dt_s is None or dt_r is None or dt_h is None or not (dt_s <= dt_r <= dt_h):
        return _res('HOLD_MISSING_OR_INVALID_LOCKED_CLOCK_TIMESTAMPS', False)

    req_flat_means = (
        'obs_stage7_total_maturation_days',
        'obs_stage7_days_post_sort_reaggregation',
        'obs_stage7_eta_sort_viable_recovery',
        'obs_readout2_inert_tracer_tau_chamber_min',
        'obs_readout2_pulse2_net_increment_z2',
        'obs_readout2_adult_islet_ref_z2',
        'obs_readout2_low2_frac_peak1',
        'obs_readout2_low2_over_low1_ratio',
        'obs_stage7_stored_mature_insulin_iu_implanted',
        'obs_stage7_stored_mature_insulin_max_lot_iu',
        'obs_stage7_proinsulin_to_cpeptide_molar_ratio',
    )
    if any(k not in validated_stage0 or se_provenance_by_mean.get(k) != 'EXPLICIT_SE_SUPPLIED' for k in req_flat_means):
        return _res('HOLD_IIP0029_INCOMPLETE_OR_DEFAULT_SE_FLAT_PANEL', False)

    if upper95_fn('obs_readout2_inert_tracer_tau_chamber_min') > float(C['TRACER_TAU_CHAMBER_MAX_MIN']):
        return _res('HOLD_CHAMBER_DISPERSION_TAU_EXCEEDED_REQUIRES_FAST_CHIP', False)
    if lower95_fn('obs_readout2_pulse2_net_increment_z2') < 0.50 or lower95_fn('obs_readout2_adult_islet_ref_z2') < 0.50:
        return _res('HOLD_NET_INCREMENT_PER_PRESORT_CELL_L95_INSUFFICIENT', False)
    if lower95_fn('obs_readout2_pulse2_net_increment_z2') < float(C['ADULT_ISLET_REF_PULSE2_MIN_FRACTION']) * validated_stage0['obs_readout2_adult_islet_ref_z2']:
        return _res('HOLD_SECOND_PULSE_BELOW_ADULT_ISLET_REFERENCE_FLOOR', False)
    if upper95_fn('obs_readout2_low2_frac_peak1') > float(C['DYNAMIC_LOW2_FRAC_PEAK1_MAX']) or upper95_fn('obs_readout2_low2_over_low1_ratio') > float(C['LOW2_OVER_LOW1_WASHOUT_MAX']):
        return _res('HOLD_GATE0M_BASAL_FRAC_PEAK_OR_LOW2_WASHOUT_FAILED', False)

    is_ext_flat = (
        validated_stage0['obs_stage7_total_maturation_days'] > float(C['BASELINE_WINDOW_MAX_STAGE7_DAYS'])
        or validated_stage0['obs_stage7_days_post_sort_reaggregation'] > float(C['BASELINE_WINDOW_MAX_POST_SORT_DAYS'])
    )
    if is_ext_flat:
        if 'obs_stage7_s8_fieller_z2_diff' not in validated_stage0 or se_provenance_by_mean.get('obs_stage7_s8_fieller_z2_diff') != 'EXPLICIT_SE_SUPPLIED':
            return _res('HOLD_S8_EXTENDED_CLOCK_REQUIRES_SOURCE_LINKED_EARLIER_WINDOW', False)
        if lower95_fn('obs_stage7_s8_fieller_z2_diff') < 0.0:
            return _res('HOLD_S8_EXTENDED_MATURATION_FIELLER_DEGRADED_KEEP_EARLIER_WINDOW', False)

    if lower95_fn('obs_stage7_stored_mature_insulin_iu_implanted') < float(C['STORED_INSULIN_FLOOR_AT_300K_BETA_IU']):
        return _res('HOLD_STORED_INSULIN_L95_BELOW_PHYSIOLOGICAL_OR_MASS_FLOOR', False)
    if upper95_fn('obs_stage7_proinsulin_to_cpeptide_molar_ratio') > float(C['PROINSULIN_TO_CPEPTIDE_CONTENT_MOLAR_MAX']):
        return _res('HOLD_ELEVATED_PROINSULIN_TO_CPEPTIDE_CONTENT_RATIO', False)

    u95_store_flat = upper95_fn('obs_stage7_stored_mature_insulin_iu_implanted')
    max_lot_flat = validated_stage0['obs_stage7_stored_mature_insulin_max_lot_iu']
    ceiling_iu = float(C['GATE0S_ABSOLUTE_STORED_INSULIN_CEILING_IU'])
    if u95_store_flat > ceiling_iu or max_lot_flat > ceiling_iu:
        g0s_ok, g0s_status = _check_scope_matched_gate0s_escalation()
        if not g0s_ok:
            return _res(
                g0s_status,
                False,
                gate0s_escalation_required=True,
                stored_insulin_u95_iu=round(u95_store_flat, 4),
                stored_insulin_max_lot_iu=round(max_lot_flat, 4),
            )
        return _res(
            'CLEARED_IIP0029_WITH_ESCALATED_SCOPE_MATCHED_GATE0S_REVIEW',
            True,
            gate0s_escalation_required=True,
            gate0s_escalation_status=g0s_status,
            stored_insulin_u95_iu=round(u95_store_flat, 4),
            stored_insulin_max_lot_iu=round(max_lot_flat, 4),
        )

    return _res(
        'CLEARED_IIP0029_TWOCLOCK_REPEATED_CHALLENGE_FLAT_PANEL',
        True,
        gate0s_escalation_required=False,
        is_extended_s8_exposure=is_ext_flat,
        stored_insulin_u95_iu=round(u95_store_flat, 4),
        stored_insulin_max_lot_iu=round(max_lot_flat, 4),
        proinsulin_to_cpeptide_marker_classification='EXPLORATORY_STAGE0_IN_VITRO_PROCESSING_MARKER_PENDING_CALIBRATION',
    )


# ============================================================================
# Stage-0 Locus-2 CIITA Mode-B v3 Comparators, Rule S11, & Auditable GY-018
# Per-Fraction Binomial-Thinning Tier 6B Evaluator (IIP-0031; DEC-110)
# ============================================================================
TIER6B_PER_FRACTION_LOT_CAP_IIP0031: float = 3.0e-7
TIER6B_DETECTION_FLOOR_IIP0031: float = 0.50
TIER6B_ALLOWED_PRESORT_LINEAGE_PANEL_IIP0031: str = 'CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY'
TIER6B_ALLOWED_ENDOGENOUS_BASES_IIP0031: frozenset[str] = frozenset({
    'MATCHED_UNTREATED_CLONOGENIC_CONTROL',
    'VALIDATED_CONSERVATIVE_BOUND',
})
TIER6B_ALLOWED_RISK_UNITS_IIP0031: frozenset[str] = frozenset({
    'SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY',
    'CD31_NEG_CD144_NEG_CLONOGENIC_OPPORTUNITY',
})
TIER6B_ALLOWED_TIMEPOINTS_IIP0031: frozenset[str] = frozenset({
    'STAGE7_PRE_ASSEMBLY',
    'PRE_ASSEMBLY_STAGE7',
    'STAGE0_PRE_ASSEMBLY_STAGE7',
    'STAGE0_STAGE7_PRE_ASSEMBLY',
})
TIER6B_ALLOWED_STAGE_LABELS_IIP0031: frozenset[str] = frozenset({
    'STAGE7_PRE_ASSEMBLY',
    'PRE_ASSEMBLY_STAGE7',
    'STAGE0_PRE_ASSEMBLY_STAGE7',
    'STAGE0_STAGE7_PRE_ASSEMBLY',
})
TIER6B_ALLOWED_ASSAY_VERSIONS_IIP0031: frozenset[str] = frozenset({
    'STAGE0_IIP0031_V1',
    'STAGE0_IIP0031_TIER6B_V1',
    'TIER6B_SORTED_LOT_V1',
})
TIER6B_ALLOWED_UNITS_IIP0031: frozenset[str] = frozenset({
    'clonogenic_opportunities',
    'sorted_clonogenic_opportunities',
})
TIER6B_FORBIDDEN_RISK_UNITS_IIP0031: frozenset[str] = frozenset({
    'TOTAL_CELLS',
    'TOTAL_CELL_COUNT',
    'SPIKE_IN_ONLY',
    'RAW_CELLS',
    'UNNORMALIZED_CELLS',
    'WHOLE_LOT_CELLS',
})
RULE_S11_GATE_ID_IIP0031: str = 'S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate'
IIP0031_BARRIER_CHS4_250_CLASSIFICATION: str = 'COMBINED_LENGTH_AND_SEQUENCE_CHANGE_250BP_PLUS_U78775_2_C79T_C80A'


def _canonicalize_product_id_detailed_iip0031(
    pvid: Any,
    root_path: Path,
    *,
    _trusted_p_mod: Any = None,
) -> tuple[str | None, dict[str, Any] | None, bool]:
    """
    Resolve `pvid` against the verified product registry, distinguishing:
      - `registry_ok=False`: product registry or pinned construct FASTA failed integrity validation;
      - `registry_ok=True, pdef=None`: registry is valid, `pvid` is malformed/missing/unregistered;
      - `registry_ok=True, pdef=dict`: verified canonical product definition.
    """
    p_mod = _trusted_p_mod if _trusted_p_mod is not None else _load_product_contract_module(root_path)
    try:
        p_mod.load_product_contracts(root_path)
    except (ValueError, OSError):
        return None, None, False
    if not isinstance(pvid, str) or not pvid.strip() or pvid != pvid.strip():
        return None, None, True
    try:
        pdef = p_mod.get_product_definition(pvid, root=root_path)
        return str(pdef['product_variant_id']), pdef, True
    except ValueError:
        return pvid, None, True


def _canonicalize_product_id_iip0031(pvid: Any, root_path: Path) -> tuple[str | None, dict[str, Any] | None]:
    canon_id, pdef, _reg_ok = _canonicalize_product_id_detailed_iip0031(pvid, root_path)
    return canon_id, pdef


def _exact_one_sided_cp_lower_iip0031(r_k: int, r_n: int, alpha_each: float) -> float:
    return clopper_pearson_lower95_independent(r_k, r_n, alpha=alpha_each)


def _conservative_int_product_floor_iip0031(int_count: int, cp_lower_float: float) -> int:
    """Exact integer floor of `int_count * cp_lower_float` without binary float multiplication rounding (#24)."""
    if int_count <= 0 or cp_lower_float <= 0.0:
        return 0
    num, den = cp_lower_float.as_integer_ratio()
    return (int_count * num) // den


def _worst_case_weighted_u95_over_bounds_iip0031(
    bounds: dict[str, tuple[float, float]],
    p_u_by_frac: dict[str, float],
) -> tuple[float, dict[str, float]]:
    """Maximize sum_i w_i * p_i_U subject to w_i_lo <= w_i <= w_i_hi and sum_i w_i == 1."""
    w_alloc = {fid: lo for fid, (lo, _hi) in bounds.items()}
    rem = max(0.0, 1.0 - sum(w_alloc.values()))
    for fid in sorted(bounds.keys(), key=lambda f: p_u_by_frac[f], reverse=True):
        if rem <= 0.0:
            break
        lo, hi = bounds[fid]
        add = min(rem, max(0.0, hi - lo))
        w_alloc[fid] += add
        rem -= add
    p_lot = sum(w_alloc[fid] * p_u_by_frac[fid] for fid in bounds)
    return float(p_lot), w_alloc


def compute_tier6b_per_fraction_thinning_bound_math_only(
    fractions: list[dict[str, Any]],
    *,
    lot_cap: float = TIER6B_PER_FRACTION_LOT_CAP_IIP0031,
    detection_floor: float = TIER6B_DETECTION_FLOOR_IIP0031,
) -> dict[str, Any]:
    """
    Standalone pure-arithmetic per-fraction binomial-thinning boundary helper (IIP-0031; DEC-110).
    Strictly non-authorizing (`release_clearance_permitted=False`): NEVER callable by live release,
    product-contract, atomic-ledger, JSON-ingestion, or replay evaluators, which ALWAYS require raw
    integer dual-calibrant counts `(r_k_hpsc, r_n_hpsc)` and `(r_k_diff, r_n_diff)` per observation.
    """
    if not isinstance(fractions, (list, tuple)) or len(fractions) < 2:
        raise ValueError('fractions must be a sequence of at least 2 fraction specifications')
    if not _is_strict_finite_number_iip0024(lot_cap) or not (0.0 < float(lot_cap) <= 1.0):
        raise ValueError('lot_cap must be a finite numeric bound in (0, 1]')
    if not _is_strict_finite_number_iip0024(detection_floor) or not (0.0 < float(detection_floor) <= 1.0):
        raise ValueError('detection_floor must be a finite numeric bound in (0, 1]')
    m = len(fractions)
    alpha_obs = 0.05 / (2.0 * float(m))
    p_i_u: dict[str, float] = {}
    d_by_f: dict[str, float] = {}
    w_by_f: dict[str, float] = {}
    n_by_f: dict[str, int] = {}
    for entry in fractions:
        if not isinstance(entry, dict):
            raise ValueError('each fraction specification must be a dictionary')
        fid = str(entry.get('fraction_id', '')).strip()
        n_i = entry.get('N')
        w_i = entry.get('w_i')
        d_l = entry.get('d_i_L')
        k_i = entry.get('k', 0)
        if (
            not fid
            or fid in p_i_u
            or type(k_i) is not int
            or k_i != 0
            or type(n_i) is not int
            or n_i < 1
            or n_i.bit_length() > 53
            or not _is_strict_finite_number_iip0024(w_i)
            or not (0.0 < float(w_i) <= 1.0)
            or not _is_strict_finite_number_iip0024(d_l)
            or not (0.0 < float(d_l) <= 1.0)
        ):
            raise ValueError(f'invalid math-only fraction specification: {entry!r}')
        w_f = float(w_i)
        d_f = float(d_l)
        w_by_f[fid] = w_f
        d_by_f[fid] = d_f
        n_by_f[fid] = n_i
        p_i_u[fid] = min(1.0, (-math.expm1(math.log(alpha_obs) / float(n_i))) / d_f)
    if abs(sum(w_by_f.values()) - 1.0) > 1e-9:
        raise ValueError('fraction weights w_i must sum to 1.0')
    p_lot_u95 = float(sum(w_by_f[fid] * p_i_u[fid] for fid in w_by_f))
    floor_ok = all(d_by_f[fid] >= detection_floor for fid in d_by_f)
    return {
        'math_only_non_authorizing': True,
        'release_clearance_permitted': False,
        'm_fractions': m,
        'alpha_obs': alpha_obs,
        'N_by_fraction': n_by_f,
        'w_by_fraction': w_by_f,
        'd_L_by_fraction': d_by_f,
        'p_i_U_by_fraction': p_i_u,
        'p_lot_U95': p_lot_u95,
        'threshold': lot_cap,
        'detection_floor_cleared': floor_ok,
        'within_numerical_cap': bool(floor_ok and p_lot_u95 <= lot_cap),
    }


_TIER6B_SAMPLE_ALIASES = (
    'sample_id', 'sample_unit_id', 'source_sample_id', 'aliquot_id',
    'specimen_id', 'original_sample_id', 'parent_sample_id',
)
_TIER6B_RECORD_ID_ALIASES = (
    'record_id', 'source_record_id', 'observation_id', 'measurement_id', 'entry_id', 'row_id',
    'atom_id', 'event_id', 'uid', 'id',
)
_TIER6B_FIELD_ALIASES = {
    'value_or_k': 'k', 'outgrowth_count': 'k', 'survivors_k': 'k', 'colonies_k': 'k',
    'outgrowth_k': 'k', 'outgrowth_count_k': 'k', 'survivor_count_k': 'k',
    'se_or_n': 'N', 'n': 'N', 'seed_line_id': 'donor_seedline_id',
    'treated_plated_cells_n': 'plated_treated_cells_n',
}
_TIER6B_PROVENANCE_KEYS = (
    'product_variant_id', 'construct_sha256', 'background_id',
    'donor_seedline_id', 'seed_line_id', 'mcb_id',
    'lot_id', 'current_lot_id', 'campaign_id', 'assay_campaign_id', 'scope_type',
)


def _tier6b_field_value(key: str, value: Any) -> Any:
    """Normalize equivalent wire encodings without coercing strings, booleans or counts."""
    if key == 'declared_fraction_roster' and isinstance(value, (list, tuple)):
        if all(isinstance(label, str) for label in value):
            return sorted(value)
    if key in ('hpsc_calibrant', 'edited_differentiated_calibrant'):
        if isinstance(value, dict) and set(value) == {'r_k', 'r_n'}:
            return [value['r_k'], value['r_n']]
    return value


def _tier6b_signature(value: Any) -> str:
    # Python container equality equates True, 1 and 1.0. Counts must retain their wire types.
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


_Tier6BFields = dict[str, Any]
_Tier6BRawAliquot = tuple[_Tier6BFields, Any, list[Any]]


def _extract_tier6b_outer_and_full_entry_iip0031(
    item: Any,
) -> tuple[_Tier6BFields, _Tier6BFields, Any, list[Any], int, str | None] | None:
    """Preserve the outer provenance, all raw counts and identity aliases before merging any fields."""
    if isinstance(item, AtomicObservation):
        sk = item.scope_key
        outer: dict[str, Any] = {
            'endpoint_id': item.endpoint_id, 'scope_type': item.scope_type, 'kind': item.kind,
            'k': item.value_or_k, 'N': item.se_or_n, 'sample_id': item.sample_unit_id, 'sample_unit_id': item.sample_unit_id,
            'product_variant_id': sk.product_variant_id, 'construct_sha256': sk.construct_sha256,
            'donor_seedline_id': sk.donor_seedline_id, 'mcb_id': sk.mcb_id,
            'sample_n': item.sample_n, 'unit': item.unit, 'assay_version': item.assay_version,
            'timepoint': item.timepoint,
        }
        if getattr(item, '_campaign_explicit', True):
            outer['campaign_id'] = item.campaign_id
        if isinstance(sk, (LotScopeKey, AssembledDeviceLotScopeKey)):
            outer['lot_id'] = sk.lot_id
        controls = _deep_to_jsonable_iip0024(item.controls)
    elif isinstance(item, dict):
        if any(not isinstance(rk, str) for rk in item):
            outer = {key: value for key, value in item.items() if isinstance(key, str) and key != 'controls'}
            controls = item.get('controls')
            _nonstring_outer = True
        else:
            outer = {key: value for key, value in item.items() if key != 'controls'}
            controls = item.get('controls')
            _nonstring_outer = False
    else:
        return None

    if 'donor_seedline_id' not in outer and 'seed_line_id' in outer:
        outer['donor_seedline_id'] = outer['seed_line_id']
    if 'lot_id' not in outer and 'current_lot_id' in outer:
        outer['lot_id'] = outer['current_lot_id']
    if 'campaign_id' not in outer and 'assay_campaign_id' in outer:
        outer['campaign_id'] = outer['assay_campaign_id']

    conflict: str | None = 'MALFORMED_TIER6B_RECORD_REJECTED' if (isinstance(item, dict) and _nonstring_outer) else None
    if controls is None:
        controls = {}
    elif not isinstance(controls, dict):
        controls = {}
        conflict = 'MALFORMED_TIER6B_RECORD_REJECTED'
    elif any(not isinstance(rk, str) for rk in controls):
        controls = {k: v for k, v in controls.items() if isinstance(k, str)}
        conflict = 'MALFORMED_TIER6B_RECORD_REJECTED'

    # Resolve nested endpoint_id and provenance conflicts BEFORE collecting positive counts (#6, #8, T9):
    nested_prov_conflict = False
    if 'endpoint_id' in controls and controls['endpoint_id'] is not None:
        if outer.get('endpoint_id') is None:
            outer['endpoint_id'] = controls['endpoint_id']
        elif outer['endpoint_id'] != controls['endpoint_id']:
            nested_prov_conflict = True
            conflict = 'CONFLICTING_OUTER_AND_NESTED_PROVENANCE'

    _prov_equiv_groups = (
        ('product_variant_id',),
        ('construct_sha256',),
        ('background_id',),
        ('donor_seedline_id', 'seed_line_id'),
        ('mcb_id',),
        ('lot_id', 'current_lot_id'),
        ('campaign_id', 'assay_campaign_id'),
        ('scope_type',),
    )
    for grp in _prov_equiv_groups:
        o_vals = [outer[k] for k in grp if k in outer and outer[k] is not None]
        c_vals = [controls[k] for k in grp if k in controls and controls[k] is not None]
        if o_vals and c_vals and any(cv != o_vals[0] for cv in c_vals):
            nested_prov_conflict = True
            conflict = 'CONFLICTING_OUTER_AND_NESTED_PROVENANCE'
        elif not o_vals and c_vals:
            if len(set(repr(cv) for cv in c_vals)) > 1:
                nested_prov_conflict = True
                conflict = 'CONFLICTING_OUTER_AND_NESTED_PROVENANCE'
            else:
                outer[grp[0]] = c_vals[0]

    for id_k in _TIER6B_SAMPLE_ALIASES + _TIER6B_RECORD_ID_ALIASES:
        if id_k in outer and id_k in controls and outer[id_k] != controls[id_k]:
            conflict = f'CONFLICTING_OUTER_AND_NESTED_IDENTITY:{id_k}'

    fields: dict[str, Any] = {}
    aliases: list[Any] = []
    outer_positive = 0
    controls_positive = 0
    for src_idx, source in enumerate((outer, controls)):
        for raw_key, value in source.items():
            if raw_key in _TIER6B_SAMPLE_ALIASES:
                if value is not None:
                    aliases.append(value)
                continue
            key = _TIER6B_FIELD_ALIASES.get(raw_key, raw_key)
            if key == 'k' and type(value) is int and value > 0:
                if src_idx == 0:
                    outer_positive = max(outer_positive, value)
                else:
                    controls_positive = max(controls_positive, value)
            value = _tier6b_field_value(key, value)
            if key in fields:
                try:
                    if _tier6b_signature(fields[key]) != _tier6b_signature(value):
                        if conflict is None:
                            conflict = 'CONFLICTING_TIER6B_OBSERVATION_FIELDS'
                except (TypeError, ValueError, OverflowError, RecursionError):
                    conflict = 'MALFORMED_TIER6B_RECORD_REJECTED'
            else:
                fields[key] = value

    # Never let a nested positive with conflicting/foreign nested provenance inherit an outer zero's target lot (#6):
    positive_count = outer_positive if nested_prov_conflict else max(outer_positive, controls_positive)

    # Background / fraction / record_id carried inside controls for typed atoms use the same binding:
    for key in ('fraction_id', 'record_id', 'background_id'):
        if key not in outer and key in fields:
            outer[key] = fields[key]
    primary_id = next((source[key] for source in (outer, controls)
                       for key in ('sample_id', 'sample_unit_id') if source.get(key) is not None), None)
    outer['sample_id'] = primary_id
    return outer, fields, primary_id, aliases, positive_count, conflict


def _tier6b_distinct_aliquots(
    records: list[_Tier6BRawAliquot],
) -> tuple[list[tuple[_Tier6BFields, list[str]]], list[str], str | None]:
    """Check the complete payload of every member of each transitive physical-identity group."""
    parents: dict[str, str] = {}

    def representative(token: str) -> str:
        parents.setdefault(token, token)
        while token != parents[token]:
            parents[token] = parents[parents[token]]
            token = parents[token]
        return token

    for fields, primary_id, aliases in records:
        if not isinstance(primary_id, str) or not primary_id or primary_id.strip() != primary_id:
            return [], [], 'MISSING_OR_INVALID_PROVENANCE_FIELD:sample_id'
        for token in aliases:
            if not isinstance(token, str) or not token or token.strip() != token:
                return [], [], 'INVALID_SAMPLE_OR_SOURCE_SAMPLE_ID'
            left, right = representative(primary_id), representative(token)
            parents[max(left, right)] = min(left, right)

    by_sample = {}
    by_record_id: dict[tuple[str, str], tuple[str, str]] = {}
    duplicate_ids = set()
    try:
        for fields, primary_id, aliases in records:
            canonical_id = representative(primary_id)
            rec_id_vals: list[tuple[str, str]] = []
            for rk_alias in _TIER6B_RECORD_ID_ALIASES:
                if rk_alias in fields and fields[rk_alias] is not None:
                    rv = fields[rk_alias]
                    if not isinstance(rv, str) or not rv.strip() or rv != rv.strip():
                        return [], [], f'INVALID_TIER6B_RECORD_ID:{rk_alias}'
                    rec_id_vals.append((rk_alias, rv))
                elif rk_alias in fields and fields[rk_alias] is None:
                    return [], [], f'INVALID_TIER6B_RECORD_ID:{rk_alias}'
            if len({rv for _, rv in rec_id_vals}) > 1:
                return [], [], 'CONFLICTING_TIER6B_RECORD_ID_ALIASES'
            signature = _tier6b_signature({k: v for k, v in fields.items() if k not in _TIER6B_RECORD_ID_ALIASES})
            for rk_alias, rv in rec_id_vals:
                for ns_key in (('any_record_id', rv), (rk_alias, rv)):
                    if ns_key in by_record_id:
                        prev_canon_sid, prev_rec_sig = by_record_id[ns_key]
                        if canonical_id != prev_canon_sid or signature != prev_rec_sig:
                            return [], [], f'CONFLICTING_OR_RELABELED_DUPLICATE_RECORD_ID:{rv}'
                    else:
                        by_record_id[ns_key] = (canonical_id, signature)
            if canonical_id in by_sample:
                previous, previous_signature, component_aliases = by_sample[canonical_id]
                if fields.get('fraction_id') != previous.get('fraction_id'):
                    return [], [], f'DUPLICATE_SAMPLE_ID_ACROSS_FRACTIONS:{canonical_id}'
                if signature != previous_signature:
                    return [], [], f'CONFLICTING_REPLAYED_SAMPLE_ID:{canonical_id}'
                component_aliases.update(aliases)
                duplicate_ids.add(canonical_id)
            else:
                by_sample[canonical_id] = (fields, signature, set(aliases))
    except (TypeError, ValueError, OverflowError, RecursionError):
        return [], [], 'MALFORMED_TIER6B_RECORD_REJECTED'
    distinct = [
        (dict(fields, sample_id=sid), sorted(aliases))
        for sid, (fields, _, aliases) in sorted(by_sample.items())
    ]
    return distinct, sorted(duplicate_ids), None


def _is_valid_sha256_hex_iip0031(val: Any) -> bool:
    return isinstance(val, str) and len(val) == 64 and all(c in '0123456789abcdef' for c in val)


def _is_attributable_to_target_lot_iip0031(
    rec: dict[str, Any],
    *,
    target_canon_pvid: str | None,
    target_construct_sha256: str | None,
    target_background_id: str | None,
    target_seedline_id: str | None,
    target_mcb_id: str | None,
    target_lot_id: str | None,
    target_campaign_id: str | None,
    target_scope_complete: bool,
    root_path: Path,
) -> tuple[bool, str]:
    if rec.get('endpoint_id') != 'lot_tier6b_sorted_prolif_clonogenic_escape':
        return False, 'non_sorted_lot_tier6b_endpoint'
    st = rec.get('scope_type', 'LOT_SCOPE')
    if st != 'LOT_SCOPE':
        return False, f'wrong_scope_type:{st}'
    if not target_scope_complete:
        return False, 'incomplete_or_invalid_target_scope'
    r_pvid_raw = rec.get('product_variant_id')
    if r_pvid_raw is None:
        return False, 'missing_product_variant_id'
    r_pvid, r_pdef = _canonicalize_product_id_iip0031(r_pvid_raw, root_path)
    if r_pdef is None or r_pvid != target_canon_pvid:
        return False, 'product_variant_mismatch'
    r_lot = rec.get('lot_id', rec.get('current_lot_id'))
    if not isinstance(r_lot, str) or not r_lot.strip() or r_lot != r_lot.strip():
        return False, 'missing_lot_id'
    if r_lot != target_lot_id:
        return False, 'wrong_lot_id'
    r_sha = rec.get('construct_sha256')
    if not _is_valid_sha256_hex_iip0031(r_sha) or r_sha != target_construct_sha256:
        return False, 'wrong_construct_sha256'
    r_seed = rec.get('donor_seedline_id', rec.get('seed_line_id'))
    if not isinstance(r_seed, str) or not r_seed.strip() or r_seed != r_seed.strip() or r_seed != target_seedline_id:
        return False, 'wrong_seed_line_id'
    r_mcb = rec.get('mcb_id')
    if not isinstance(r_mcb, str) or not r_mcb.strip() or r_mcb != r_mcb.strip() or r_mcb != target_mcb_id:
        return False, 'wrong_mcb_id'
    r_bg = rec.get('background_id')
    if r_bg is not None:
        if not isinstance(r_bg, str) or not r_bg.strip() or r_bg != r_bg.strip():
            return False, 'wrong_background_id'
        if target_background_id is not None and r_bg != target_background_id:
            return False, 'wrong_background_id'
    if 'campaign_id' in rec or 'assay_campaign_id' in rec:
        r_camp = rec.get('campaign_id', rec.get('assay_campaign_id'))
        if r_camp is not None:
            if not isinstance(r_camp, str) or not r_camp.strip() or r_camp != r_camp.strip():
                return False, 'invalid_campaign_id'
            if target_campaign_id is not None and r_camp != target_campaign_id:
                return False, 'wrong_campaign_id'
        else:
            return False, 'invalid_campaign_id'
    return True, 'matched'


def evaluate_tier6b_per_fraction_thinning_iip0031(
    records: Any,
    *,
    target_scope: dict[str, Any],
    root: Path | str | None = None,
) -> dict[str, Any]:
    """
    Auditable GY-018-compliant pre-assembly Stage-7 CD31-/CD144- sorted-lot Tier 6B
    per-fraction binomial-thinning evaluator (IIP-0031; DEC-110).
    """
    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    target_scope_malformed = not isinstance(target_scope, dict)
    if target_scope_malformed:
        target_scope = {}
    t_pvid, t_pdef = _canonicalize_product_id_iip0031(target_scope.get('product_variant_id'), root_path)
    t_sha_raw = target_scope.get('construct_sha256')
    t_sha = t_sha_raw if _is_valid_sha256_hex_iip0031(t_sha_raw) else None
    t_bg_raw = target_scope.get('background_id')
    if (not isinstance(t_bg_raw, str) or not t_bg_raw.strip()) and isinstance(t_pdef, dict):
        l2b = t_pdef.get('staged_locus2_modeb_donor_binding')
        if isinstance(l2b, dict) and isinstance(l2b.get('background_id'), str):
            t_bg_raw = l2b['background_id']
    t_bg = t_bg_raw if isinstance(t_bg_raw, str) and t_bg_raw.strip() and t_bg_raw == t_bg_raw.strip() else None
    t_seed_raw = target_scope.get('donor_seedline_id', target_scope.get('seed_line_id'))
    t_seed = t_seed_raw if isinstance(t_seed_raw, str) and t_seed_raw.strip() and t_seed_raw == t_seed_raw.strip() else None
    t_mcb_raw = target_scope.get('mcb_id')
    t_mcb = t_mcb_raw if isinstance(t_mcb_raw, str) and t_mcb_raw.strip() and t_mcb_raw == t_mcb_raw.strip() else None
    t_lot_raw = target_scope.get('lot_id', target_scope.get('current_lot_id'))
    t_lot = t_lot_raw if isinstance(t_lot_raw, str) and t_lot_raw.strip() and t_lot_raw == t_lot_raw.strip() else None
    t_camp_raw = target_scope.get('campaign_id', target_scope.get('assay_campaign_id'))
    t_camp_invalid = (
        ('campaign_id' in target_scope or 'assay_campaign_id' in target_scope)
        and t_camp_raw is not None
        and (not isinstance(t_camp_raw, str) or not t_camp_raw.strip() or t_camp_raw != t_camp_raw.strip())
    )
    t_camp = t_camp_raw if isinstance(t_camp_raw, str) and t_camp_raw.strip() and t_camp_raw == t_camp_raw.strip() else None

    t_l2_sha_mismatch = False
    if isinstance(t_pdef, dict) and isinstance(t_pdef.get('staged_locus2_modeb_donor_binding'), dict):
        exp_l2_sha = t_pdef['staged_locus2_modeb_donor_binding'].get('locus2_donor_sha256')
        if isinstance(exp_l2_sha, str) and t_sha != exp_l2_sha:
            t_l2_sha_mismatch = True

    target_scope_complete = bool(
        not target_scope_malformed
        and not t_camp_invalid
        and not t_l2_sha_mismatch
        and t_pvid is not None
        and t_pdef is not None
        and t_sha is not None
        and t_seed is not None
        and t_mcb is not None
        and t_lot is not None
    )

    def _hold_unqual(reason: str, **extra: Any) -> dict[str, Any]:
        return {
            'gate_pass': False,
            'status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION',
            'disposition': 'HOLD_TIER6B_UNQUALIFIED_FRACTION',
            'reason_code': reason,
            'tier6b_status': 'HOLD_TIER6B_UNQUALIFIED_FRACTION',
            'hard_stop': False,
            'p_lot_U95': extra.pop('p_lot_U95', None),
            'threshold': TIER6B_PER_FRACTION_LOT_CAP_IIP0031,
            'd_L_by_fraction': extra.pop('d_L_by_fraction', {}),
            'd_L_details': extra.pop('d_L_details', {}),
            'detection_floor_cleared': extra.pop('detection_floor_cleared', False),
            'deduplicated_sample_ids': extra.pop('deduplicated_sample_ids', []),
            **extra,
        }

    if not isinstance(records, (list, tuple)):
        return _hold_unqual('TIER6B_RECORDS_MUST_BE_SEQUENCE')

    malformed_seen = False
    same_lot_provenance_conflict = False
    unlinked_campaign_positive = False
    latched_outgrowths: list[dict[str, Any]] = []
    unlinked_outgrowths: list[dict[str, Any]] = []
    record_conflicts: set[str] = set()
    target_recs: list[_Tier6BRawAliquot] = []
    ignored_out_of_scope: list[dict[str, Any]] = []

    # PASS 1: Scan EVERY entry's outer lot identity + strict integer k_i >= 1 BEFORE parsing ancillary controls.
    # Foreign-lot records are excluded BEFORE updating target-lot `malformed_seen` or `record_conflicts` (#23, T5)!
    for idx, raw_item in enumerate(records):
        observation = _extract_tier6b_outer_and_full_entry_iip0031(raw_item)
        if observation is None:
            malformed_seen = True
            continue
        outer_rec, full_rec, primary_id, aliases, positive_count, conflict = observation
        is_attr, attr_reason = _is_attributable_to_target_lot_iip0031(
            outer_rec,
            target_canon_pvid=t_pvid,
            target_construct_sha256=t_sha,
            target_background_id=t_bg,
            target_seedline_id=t_seed,
            target_mcb_id=t_mcb,
            target_lot_id=t_lot,
            target_campaign_id=t_camp,
            target_scope_complete=target_scope_complete,
            root_path=root_path,
        )
        if not is_attr:
            if positive_count > 0:
                unlinked_outgrowths.append({
                    'index': idx,
                    'reason': attr_reason,
                    'lot_id': outer_rec.get('lot_id'),
                    'campaign_id': outer_rec.get('campaign_id'),
                    'k': positive_count,
                })
            if (
                outer_rec.get('endpoint_id') == 'lot_tier6b_sorted_prolif_clonogenic_escape'
                and outer_rec.get('scope_type', 'LOT_SCOPE') == 'LOT_SCOPE'
                and t_lot is not None
                and outer_rec.get('lot_id') == t_lot
            ):
                if attr_reason == 'wrong_campaign_id' and positive_count > 0:
                    unlinked_campaign_positive = True
                elif attr_reason in ('wrong_construct_sha256', 'wrong_seed_line_id', 'wrong_mcb_id', 'wrong_background_id', 'invalid_campaign_id'):
                    same_lot_provenance_conflict = True
            elif attr_reason.startswith('wrong_scope_type:') and outer_rec.get('lot_id') == t_lot:
                same_lot_provenance_conflict = True
            elif outer_rec.get('endpoint_id') != 'lot_tier6b_sorted_prolif_clonogenic_escape' and outer_rec.get('lot_id') in (None, t_lot) and not outer_rec.get('endpoint_id'):
                malformed_seen = True
            ignored_out_of_scope.append({'index': idx, 'reason': attr_reason, 'endpoint_id': outer_rec.get('endpoint_id')})
            continue

        if conflict == 'MALFORMED_TIER6B_RECORD_REJECTED':
            malformed_seen = True
        elif conflict:
            record_conflicts.add(conflict)

        if positive_count:
            latched_outgrowths.append({
                'index': idx,
                'fraction_id': outer_rec.get('fraction_id'),
                'sample_id': outer_rec.get('sample_id'),
                'record_id': outer_rec.get('record_id'),
                'k': positive_count,
                'N': full_rec.get('N'),
                'lot_id': outer_rec.get('lot_id'),
                'campaign_id': outer_rec.get('campaign_id', t_camp),
            })
        full_rec.setdefault('sample_n', None)
        target_recs.append((full_rec, primary_id, aliases))

    if latched_outgrowths:
        return {
            'gate_pass': False,
            'status': 'FAIL_HARD_STOP',
            'disposition': 'FAIL_HARD_STOP',
            'reason_code': 'FAIL_TIER6B_OUTGROWTH_DETECTED',
            'tier6b_status': 'FAIL_TIER6B_OUTGROWTH_DETECTED',
            'hard_stop': True,
            'latched_outgrowth_events': latched_outgrowths,
            'unlinked_or_foreign_outgrowth_events': unlinked_outgrowths,
            'malformed_records_observed': malformed_seen,
            'p_lot_U95': None,
            'threshold': TIER6B_PER_FRACTION_LOT_CAP_IIP0031,
            'd_L_by_fraction': {},
            'd_L_details': {},
            'detection_floor_cleared': False,
            'deduplicated_sample_ids': [],
        }

    # PASS 2: Strict provenance, completeness, per-observation calibration/CFU/mass-balance, and exact zero-outgrowth math:
    if t_l2_sha_mismatch:
        return _hold_unqual('OUTER_CONSTRUCT_SHA256_MISMATCH_WITH_LOCUS2_BINDING', unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)
    if not target_scope_complete:
        return _hold_unqual('MISSING_OR_INVALID_OUTER_TARGET_IMMUTABLE_PROVENANCE_SCOPE', unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)
    if unlinked_campaign_positive:
        return _hold_unqual('UNLINKED_CAMPAIGN_TIER6B_POSITIVE_EVENT', unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)
    if malformed_seen:
        return _hold_unqual('MALFORMED_TIER6B_RECORD_REJECTED', unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)
    if record_conflicts:
        return _hold_unqual(min(record_conflicts), unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)
    if same_lot_provenance_conflict:
        return _hold_unqual('INCONSISTENT_PROVENANCE_ACROSS_SAME_LOT_TIER6B_RECORDS', unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)
    if not target_recs:
        return _hold_unqual('MISSING_SORTED_LOT_TIER6B_EVIDENCE', ignored_out_of_scope=ignored_out_of_scope, unlinked_or_foreign_outgrowth_events=unlinked_outgrowths)

    # Examine every lot declaration, including duplicate aliquots, before deduplication (including lot-consistent timepoint/assay_version/unit, T11):
    lot_fields = {
        'presort_total_input_cells': 'INCONSISTENT_PRESORT_TOTAL_INPUT_CELLS_ACROSS_RECORDS',
        'presort_reconciled_loss_cells': 'INCONSISTENT_PRESORT_RECONCILED_LOSS_CELLS_ACROSS_RECORDS',
        'unreconciled_loss_fraction': 'INCONSISTENT_UNRECONCILED_LOSS_FRACTION_ACROSS_RECORDS',
        'declared_fraction_roster': 'INCONSISTENT_FRACTION_ROSTER_ACROSS_RECORDS',
        'timepoint': 'INCONSISTENT_TIER6B_TIMEPOINT_ACROSS_LOT_RECORDS',
        'assay_version': 'INCONSISTENT_TIER6B_ASSAY_VERSION_ACROSS_LOT_RECORDS',
        'unit': 'INCONSISTENT_TIER6B_UNIT_ACROSS_LOT_RECORDS',
    }
    try:
        lot_conflicts = [
            reason for key, reason in lot_fields.items()
            if len({_tier6b_signature(rec.get(key)) for rec, _, _ in target_recs}) > 1
        ]
    except (TypeError, ValueError, OverflowError, RecursionError):
        return _hold_unqual('MALFORMED_TIER6B_RECORD_REJECTED')
    distinct_records, deduplicated_sample_ids, duplicate_error = _tier6b_distinct_aliquots(target_recs)
    if duplicate_error:
        return _hold_unqual(duplicate_error)
    # Validate per-row metadata BEFORE lot-level consistency checks so invalid stage/assay/unit labels surface their specific reason:
    lot_roster: list[str] | None = None
    lot_risk_unit: str | None = None
    lot_bg_seen: str | None = t_bg
    lot_presort_total_cells: int | None = None
    lot_reconciled_loss_cells: int | None = None
    seen_untr_sid_to_frac: dict[str, str] = {}
    seen_cal_batch_to_frac: dict[str, str] = {}
    seen_cal_specimen_to_frac: dict[str, str] = {}
    all_treated_sample_ids: set[str] = set()
    all_untreated_sample_ids: set[str] = set()
    all_calibrant_specimen_ids: set[str] = set()
    agg_by_frac: dict[str, dict[str, Any]] = {}

    _allowed_untr_keys = frozenset({
        'untreated_sample_id', 'untreated_cfu_colonies_k', 'untreated_plated_cells_n',
        'treated', 'ap1903_treated', 'is_untreated', 'control_role',
        *_TIER6B_SAMPLE_ALIASES, *_TIER6B_RECORD_ID_ALIASES,
    })
    _allowed_cal_dict_keys = frozenset({
        'r_k', 'r_n', 'cal_batch_id', 'calibrant_role', 'treated', 'ap1903_treated',
        *_TIER6B_SAMPLE_ALIASES, *_TIER6B_RECORD_ID_ALIASES,
    })

    for rec, aliquot_tokens in distinct_records:
        if t_camp is not None and rec.get('campaign_id') != t_camp:
            return _hold_unqual('INCONSISTENT_RECORD_CAMPAIGN_WITH_TARGET_LOT')
        if rec.get('scope_type') != 'LOT_SCOPE':
            return _hold_unqual('MISSING_OR_INVALID_TIER6B_LOT_SCOPE')
        if rec.get('kind') != 'EXACT_INTEGER_COUNT':
            return _hold_unqual('NON_COUNT_OR_METADATA_TIER6B_SUBSTITUTE_REJECTED')
        tp_val = rec.get('timepoint')
        if not isinstance(tp_val, str) or not tp_val.strip() or tp_val != tp_val.strip() or tp_val not in TIER6B_ALLOWED_TIMEPOINTS_IIP0031:
            return _hold_unqual('INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT')
        for _stage_k in ('differentiation_stage', 'stage', 'cell_stage', 'developmental_stage'):
            if _stage_k in rec and rec[_stage_k] is not None:
                _sv = rec[_stage_k]
                if not isinstance(_sv, str) or not _sv.strip() or _sv != _sv.strip() or _sv not in TIER6B_ALLOWED_STAGE_LABELS_IIP0031:
                    return _hold_unqual('INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT')
        if 'assembly_state' in rec and rec['assembly_state'] is not None and rec['assembly_state'] != 'PRE_ASSEMBLY':
            return _hold_unqual('INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT')
        if 'pre_assembly' in rec and rec['pre_assembly'] is not None and rec['pre_assembly'] is not True:
            return _hold_unqual('INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT')
        av_val = rec.get('assay_version')
        if not isinstance(av_val, str) or not av_val.strip() or av_val != av_val.strip() or av_val not in TIER6B_ALLOWED_ASSAY_VERSIONS_IIP0031:
            return _hold_unqual('UNVALIDATED_TIER6B_ASSAY_VERSION')
        unit_val = rec.get('unit')
        if not isinstance(unit_val, str) or not unit_val.strip() or unit_val != unit_val.strip() or unit_val not in TIER6B_ALLOWED_UNITS_IIP0031:
            return _hold_unqual('INCOMPATIBLE_TIER6B_MEASUREMENT_UNIT')
        if rec.get('sample_n') is not None:
            return _hold_unqual('INVALID_TIER6B_SAMPLE_N_FOR_EXACT_COUNT')
        for req_id_k in ('product_variant_id', 'construct_sha256', 'donor_seedline_id', 'mcb_id', 'lot_id', 'sample_id', 'fraction_id'):
            v_id = rec.get(req_id_k)
            if not isinstance(v_id, str) or not v_id.strip() or v_id != v_id.strip():
                return _hold_unqual(f'MISSING_OR_INVALID_PROVENANCE_FIELD:{req_id_k}')
        rec_pvid_c, _ = _canonicalize_product_id_iip0031(rec['product_variant_id'], root_path)
        if (
            rec_pvid_c != t_pvid
            or rec['construct_sha256'] != t_sha
            or rec['donor_seedline_id'] != t_seed
            or rec['mcb_id'] != t_mcb
            or rec['lot_id'] != t_lot
        ):
            return _hold_unqual('INCONSISTENT_RECORD_PROVENANCE_WITH_TARGET_LOT')

        r_bg_val = rec.get('background_id')
        if t_bg is not None and r_bg_val != t_bg:
            return _hold_unqual('MISSING_OR_MISMATCHED_BACKGROUND_ID')
        if r_bg_val is not None:
            if not isinstance(r_bg_val, str) or not r_bg_val.strip() or r_bg_val != r_bg_val.strip():
                return _hold_unqual('INVALID_BACKGROUND_ID')
            if lot_bg_seen is None:
                lot_bg_seen = r_bg_val
            elif lot_bg_seen != r_bg_val:
                return _hold_unqual('INCONSISTENT_BACKGROUND_ID_ACROSS_FRACTIONS')

        # Strictly prohibit calculation-only supplied d_L, single-fraction bypass, and caller weight overrides/unverified bounds:
        if rec.get('calculation_only_supplied_d_L') not in (None, False) or 'supplied_d_i_L' in rec:
            return _hold_unqual('CALCULATION_ONLY_SUPPLIED_D_L_PROHIBITED_ON_LIVE_PATHS')
        if rec.get('single_fraction_unsorted_lot') not in (None, False):
            return _hold_unqual('SINGLE_FRACTION_UNSORTED_BYPASS_PROHIBITED_IN_SORTED_LOT_TIER6B')
        for _forb_wk in ('w_i_lo', 'w_i_hi', 'w_i_upper', 'weight_override', 'w_override'):
            if _forb_wk in rec and rec[_forb_wk] is not None:
                return _hold_unqual('CALLER_WEIGHT_OVERRIDE_OR_UNVERIFIED_BOUNDS_PROHIBITED')
        # Validate per-row `w_i` type and bounds BEFORE aggregation so lexical row order never masks Boolean/invalid `w_i` (#25):
        if 'w_i' in rec:
            w_row_raw = rec.get('w_i')
            if not _is_strict_finite_number_iip0024(w_row_raw):
                return _hold_unqual(f'INVALID_OR_MISSING_FRACTION_WEIGHT:{rec["fraction_id"]}')
            if not (0.0 < float(w_row_raw) <= 1.0):
                return _hold_unqual(f'OUT_OF_RANGE_FRACTION_WEIGHT:{rec["fraction_id"]}')

        k_raw = rec.get('k')
        n_raw = rec.get('N')
        if type(k_raw) is not int or type(n_raw) is not int:
            return _hold_unqual('NON_INTEGER_OR_BOOLEAN_TIER6B_COUNT_OR_DENOMINATOR')
        if n_raw < 1 or k_raw != 0 or k_raw > n_raw:
            return _hold_unqual('INVALID_TIER6B_DENOMINATOR_OR_COUNT_BOUNDS')

        if rec.get('clustered_sampling') is True or (
            rec.get('clustered_sampling') is not None and type(rec.get('clustered_sampling')) is not bool
        ):
            return _hold_unqual('UNSUPPORTED_CLUSTERED_SAMPLING_REJECTED')
        samp_model = rec.get('sampling_independence_model', 'INDEPENDENT_BIOLOGICAL_OPPORTUNITIES')
        if samp_model != 'INDEPENDENT_BIOLOGICAL_OPPORTUNITIES':
            return _hold_unqual('UNSUPPORTED_SAMPLING_INDEPENDENCE_MODEL')

        panel = rec.get('presort_lineage_panel')
        if panel != TIER6B_ALLOWED_PRESORT_LINEAGE_PANEL_IIP0031:
            return _hold_unqual('INVALID_OR_BLINDED_PRESORT_LINEAGE_PANEL')

        risk_unit = rec.get('clonogenic_risk_unit')
        if (
            not isinstance(risk_unit, str)
            or not risk_unit.strip()
            or risk_unit != risk_unit.strip()
            or risk_unit not in TIER6B_ALLOWED_RISK_UNITS_IIP0031
        ):
            return _hold_unqual('INVALID_OR_UNSUBSTANTIATED_CLONOGENIC_RISK_UNIT')
        if lot_risk_unit is None:
            lot_risk_unit = risk_unit
        elif lot_risk_unit != risk_unit:
            return _hold_unqual('INCONSISTENT_CLONOGENIC_RISK_UNIT_ACROSS_FRACTIONS')

        if rec.get('endogenous_denominator_verified') is not True:
            return _hold_unqual('OMITTED_OR_UNVERIFIED_ENDOGENOUS_DENOMINATOR')
        endo_basis = rec.get('endogenous_denominator_basis')
        if not isinstance(endo_basis, str) or endo_basis not in TIER6B_ALLOWED_ENDOGENOUS_BASES_IIP0031:
            return _hold_unqual('INVALID_ENDOGENOUS_DENOMINATOR_BASIS')

        # Require closed-schema source-linked disjoint untreated (no-AP1903) CFU control tuple on EVERY observation (#27, T7):
        untr_ctrl = rec.get('untreated_cfu_control')
        if not isinstance(untr_ctrl, dict):
            return _hold_unqual('MISSING_SOURCE_LINKED_UNTREATED_CFU_CONTROL')
        if any(not isinstance(uk_k, str) or uk_k not in _allowed_untr_keys for uk_k in untr_ctrl):
            return _hold_unqual('INVALID_UNTREATED_CFU_CONTROL_SCHEMA')
        for _tr_flag_k in ('treated', 'ap1903_treated'):
            if _tr_flag_k in untr_ctrl and untr_ctrl[_tr_flag_k] is not False:
                return _hold_unqual('CONTRADICTORY_TREATED_FLAG_IN_UNTREATED_CFU_CONTROL')
        if 'is_untreated' in untr_ctrl and untr_ctrl['is_untreated'] is not True:
            return _hold_unqual('CONTRADICTORY_TREATED_FLAG_IN_UNTREATED_CFU_CONTROL')
        if 'control_role' in untr_ctrl and untr_ctrl['control_role'] != 'UNTREATED_CFU_CONTROL':
            return _hold_unqual('CONTRADICTORY_CONTROL_ROLE_IN_UNTREATED_CFU_CONTROL')
        untr_sid = untr_ctrl.get('untreated_sample_id')
        untr_k = untr_ctrl.get('untreated_cfu_colonies_k')
        untr_n = untr_ctrl.get('untreated_plated_cells_n')
        if (
            not isinstance(untr_sid, str)
            or not untr_sid.strip()
            or untr_sid != untr_sid.strip()
            or type(untr_k) is not int
            or type(untr_n) is not int
            or untr_n < 1
            or not (1 <= untr_k <= untr_n)
        ):
            return _hold_unqual('INVALID_OR_ZERO_GROWTH_UNTREATED_CFU_CONTROL')
        untr_all_sids: list[str] = [untr_sid]
        for _ua_k in _TIER6B_SAMPLE_ALIASES + _TIER6B_RECORD_ID_ALIASES:
            if _ua_k in untr_ctrl:
                _ua_v = untr_ctrl[_ua_k]
                if not isinstance(_ua_v, str) or not _ua_v.strip() or _ua_v != _ua_v.strip():
                    return _hold_unqual('INVALID_OR_ZERO_GROWTH_UNTREATED_CFU_CONTROL')
                if _ua_v not in untr_all_sids:
                    untr_all_sids.append(_ua_v)
        if untr_n < IIP0031_MIN_UNTREATED_CFU_PLATED_CELLS and untr_n.bit_length() <= 53:
            return _hold_unqual('UNDERPOWERED_UNTREATED_CFU_CONTROL_DENOMINATOR')

        # Live untreated CFU assays are ALWAYS finite samples; forbid opting out of 3-way Clopper-Pearson coverage:
        if 'statistically_estimated_denominator_or_share' in rec and rec['statistically_estimated_denominator_or_share'] is not None:
            stat_est_flag = rec['statistically_estimated_denominator_or_share']
            if type(stat_est_flag) is not bool:
                return _hold_unqual('INVALID_STATISTICAL_ESTIMATION_FLAG')
            if stat_est_flag is False:
                return _hold_unqual('EXACT_FINITE_CFU_OPT_OUT_PROHIBITED')

        # Require explicit integer pre-sort cell balance, disposal attestation for sort losses, & physical plated-cell bounds:
        if rec.get('presort_mass_balance_reconciled') is not True:
            return _hold_unqual('UNRECONCILED_PRESORT_MASS_BALANCE')
        unrec_loss = rec.get('unreconciled_loss_fraction')
        if not _is_strict_finite_number_iip0024(unrec_loss) or float(unrec_loss) < 0.0 or float(unrec_loss) > 1e-12:
            return _hold_unqual('UNRECONCILED_SORTING_OR_PROCESSING_LOSSES')
        ps_tot = rec.get('presort_total_input_cells')
        ps_loss = rec.get('presort_reconciled_loss_cells')
        ps_frac = rec.get('presort_fraction_cells')
        plated_n = rec.get('plated_treated_cells_n')
        if (
            type(ps_tot) is not int
            or ps_tot < 1
            or type(ps_loss) is not int
            or ps_loss < 0
            or type(ps_frac) is not int
            or ps_frac < 1
            or type(plated_n) is not int
            or plated_n < n_raw
        ):
            return _hold_unqual('MISSING_OR_INVALID_INTEGER_PRESORT_CELL_MASS_BALANCE')
        loss_disp_flag = rec.get('excluded_or_lost_material_discarded_not_in_graft')
        if loss_disp_flag is not None and type(loss_disp_flag) is not bool:
            return _hold_unqual('UNVERIFIED_DISPOSAL_OF_RECONCILED_SORT_LOSSES')
        if loss_disp_flag is False or (ps_loss > 0 and loss_disp_flag is not True):
            return _hold_unqual('UNVERIFIED_DISPOSAL_OF_RECONCILED_SORT_LOSSES')
        lot_presort_total_cells = ps_tot
        lot_reconciled_loss_cells = ps_loss

        # Validate EVERY observation's separate dual calibrants BEFORE fraction aggregation (closed schema + disjoint IDs, #27, T8):
        if rec.get('pooled_calibrants') is True or (
            rec.get('pooled_calibrants') is not None and type(rec.get('pooled_calibrants')) is not bool
        ):
            return _hold_unqual(f'POOLED_CALIBRANTS_PROHIBITED:{rec["fraction_id"]}')
        c_hpsc = rec.get('hpsc_calibrant')
        c_diff = rec.get('edited_differentiated_calibrant')
        if not isinstance(c_hpsc, (list, tuple, dict)) or not isinstance(c_diff, (list, tuple, dict)):
            return _hold_unqual(f'MISSING_SEPARATE_DUAL_CALIBRANTS:{rec["fraction_id"]}')
        hpsc_cal_ids: list[str] = []
        diff_cal_ids: list[str] = []
        if isinstance(c_hpsc, dict):
            if any(not isinstance(ck_k, str) or ck_k not in _allowed_cal_dict_keys for ck_k in c_hpsc):
                return _hold_unqual(f'INVALID_HPSC_CALIBRANT_PAIR:{rec["fraction_id"]}')
            if 'calibrant_role' in c_hpsc and c_hpsc['calibrant_role'] != 'HPSC_CALIBRANT':
                return _hold_unqual(f'INVALID_HPSC_CALIBRANT_PAIR:{rec["fraction_id"]}')
            for _tr_ck in ('treated', 'ap1903_treated'):
                if _tr_ck in c_hpsc and c_hpsc[_tr_ck] is not True:
                    return _hold_unqual(f'INVALID_HPSC_CALIBRANT_PAIR:{rec["fraction_id"]}')
            for _cid_k in ('cal_batch_id',) + _TIER6B_SAMPLE_ALIASES + _TIER6B_RECORD_ID_ALIASES:
                if _cid_k in c_hpsc:
                    _cid_v = c_hpsc[_cid_k]
                    if not isinstance(_cid_v, str) or not _cid_v.strip() or _cid_v != _cid_v.strip():
                        return _hold_unqual(f'INVALID_HPSC_CALIBRANT_PAIR:{rec["fraction_id"]}')
                    if _cid_v not in hpsc_cal_ids:
                        hpsc_cal_ids.append(_cid_v)
            rk_h, rn_h = c_hpsc.get('r_k'), c_hpsc.get('r_n')
        elif len(c_hpsc) == 2:
            rk_h, rn_h = c_hpsc[0], c_hpsc[1]
        else:
            return _hold_unqual(f'INVALID_HPSC_CALIBRANT_PAIR:{rec["fraction_id"]}')
        if isinstance(c_diff, dict):
            if any(not isinstance(ck_k, str) or ck_k not in _allowed_cal_dict_keys for ck_k in c_diff):
                return _hold_unqual(f'INVALID_DIFF_CALIBRANT_PAIR:{rec["fraction_id"]}')
            if 'calibrant_role' in c_diff and c_diff['calibrant_role'] != 'EDITED_DIFFERENTIATED_CALIBRANT':
                return _hold_unqual(f'INVALID_DIFF_CALIBRANT_PAIR:{rec["fraction_id"]}')
            for _tr_ck in ('treated', 'ap1903_treated'):
                if _tr_ck in c_diff and c_diff[_tr_ck] is not True:
                    return _hold_unqual(f'INVALID_DIFF_CALIBRANT_PAIR:{rec["fraction_id"]}')
            for _cid_k in ('cal_batch_id',) + _TIER6B_SAMPLE_ALIASES + _TIER6B_RECORD_ID_ALIASES:
                if _cid_k in c_diff:
                    _cid_v = c_diff[_cid_k]
                    if not isinstance(_cid_v, str) or not _cid_v.strip() or _cid_v != _cid_v.strip():
                        return _hold_unqual(f'INVALID_DIFF_CALIBRANT_PAIR:{rec["fraction_id"]}')
                    if _cid_v not in diff_cal_ids:
                        diff_cal_ids.append(_cid_v)
            rk_d, rn_d = c_diff.get('r_k'), c_diff.get('r_n')
        elif len(c_diff) == 2:
            rk_d, rn_d = c_diff[0], c_diff[1]
        else:
            return _hold_unqual(f'INVALID_DIFF_CALIBRANT_PAIR:{rec["fraction_id"]}')
        if set(hpsc_cal_ids).intersection(diff_cal_ids):
            return _hold_unqual(f'NON_DISJOINT_HPSC_AND_DIFF_CALIBRANT_SPECIMEN_ID:{rec["fraction_id"]}')
        if (
            type(rk_h) is not int
            or type(rn_h) is not int
            or type(rk_d) is not int
            or type(rn_d) is not int
            or rn_h < 1
            or rn_d < 1
            or not (0 <= rk_h <= rn_h)
            or not (0 <= rk_d <= rn_d)
        ):
            return _hold_unqual(f'NON_INTEGER_OR_OUT_OF_BOUNDS_CALIBRANT_COUNTS:{rec["fraction_id"]}')
        if any(count.bit_length() > 53 for count in (n_raw, ps_tot, ps_loss, ps_frac, plated_n, untr_n, rn_h, rn_d)):
            return _hold_unqual('TIER6B_COUNTS_EXCEED_SUPPORTED_NUMERIC_PRECISION')
        if n_raw > (plated_n * untr_k) // untr_n:
            return _hold_unqual(
                'UNTREATED_CFU_DOES_NOT_SUPPORT_CLAIMED_OPPORTUNITIES',
                untreated_cfu_max_supported_N=(plated_n * untr_k) // untr_n,
                claimed_N=n_raw,
            )
        cal_tuple_obs = (rk_h, rn_h, rk_d, rn_d)

        cal_batch_raw = rec.get('cal_batch_id')
        if cal_batch_raw is not None:
            if not isinstance(cal_batch_raw, str) or not cal_batch_raw.strip() or cal_batch_raw != cal_batch_raw.strip():
                return _hold_unqual(f'INVALID_CAL_BATCH_ID:{rec["fraction_id"]}')

        roster_raw = rec.get('declared_fraction_roster')
        if (
            not isinstance(roster_raw, (list, tuple))
            or len(roster_raw) < 2
            or any(not isinstance(f, str) or not f.strip() or f != f.strip() for f in roster_raw)
            or len(set(roster_raw)) != len(roster_raw)
            or len({re.sub(r'[^a-z0-9]+', '_', f.lower()).strip('_') for f in roster_raw}) != len(roster_raw)
        ):
            return _hold_unqual('INVALID_OR_MISSING_EXHAUSTIVE_FRACTION_ROSTER')
        lot_roster = list(roster_raw)

        fid = rec['fraction_id']
        sid = rec['sample_id']
        if fid not in lot_roster:
            return _hold_unqual(f'UNDECLARED_FRACTION_ID:{fid}')

        for _usid_tok in untr_all_sids:
            if _usid_tok in aliquot_tokens:
                return _hold_unqual(f'NON_DISJOINT_UNTREATED_CFU_SAMPLE_ID:{sid}')
            if _usid_tok in seen_untr_sid_to_frac:
                if seen_untr_sid_to_frac[_usid_tok] != fid:
                    return _hold_unqual(f'NON_DISJOINT_UNTREATED_CFU_SAMPLE_ID_ACROSS_FRACTIONS:{_usid_tok}')
                return _hold_unqual(f'BORROWED_OR_SHARED_UNTREATED_CFU_CONTROL_ACROSS_OBSERVATIONS:{_usid_tok}')
        for _usid_tok in untr_all_sids:
            seen_untr_sid_to_frac[_usid_tok] = fid

        if cal_batch_raw is not None:
            if cal_batch_raw in seen_cal_batch_to_frac:
                if seen_cal_batch_to_frac[cal_batch_raw] != fid:
                    return _hold_unqual(f'NON_DISJOINT_CALIBRATION_BATCH_ID_ACROSS_FRACTIONS:{cal_batch_raw}')
                return _hold_unqual(f'BORROWED_OR_SHARED_CALIBRATION_BATCH_ACROSS_OBSERVATIONS:{cal_batch_raw}')
            seen_cal_batch_to_frac[cal_batch_raw] = fid

        for _cspec_id in hpsc_cal_ids + diff_cal_ids:
            if _cspec_id in aliquot_tokens or _cspec_id in untr_all_sids:
                return _hold_unqual(f'NON_DISJOINT_CALIBRANT_SPECIMEN_ID:{_cspec_id}')
            if _cspec_id in seen_cal_specimen_to_frac:
                if seen_cal_specimen_to_frac[_cspec_id] != fid:
                    return _hold_unqual(f'NON_DISJOINT_CALIBRATION_BATCH_ID_ACROSS_FRACTIONS:{_cspec_id}')
                return _hold_unqual(f'BORROWED_OR_SHARED_CALIBRATION_BATCH_ACROSS_OBSERVATIONS:{_cspec_id}')
        for _cspec_id in hpsc_cal_ids + diff_cal_ids:
            seen_cal_specimen_to_frac[_cspec_id] = fid

        all_treated_sample_ids.update(aliquot_tokens)
        all_untreated_sample_ids.update(untr_all_sids)
        all_calibrant_specimen_ids.update(hpsc_cal_ids)
        all_calibrant_specimen_ids.update(diff_cal_ids)
        canon_tok = sid
        cal_key_obs = ('cal_batch_id', cal_batch_raw) if cal_batch_raw is not None else ('obs_sample_id', canon_tok)

        _allowed_tier6b_roles = frozenset({'ENRICHED', 'FLOW_THROUGH_OR_EXCLUDED'})
        _explicit_frole: str | None = None
        if 'fraction_role' in rec:
            _raw_frole = rec.get('fraction_role')
            if (
                not isinstance(_raw_frole, str)
                or not _raw_frole.strip()
                or _raw_frole != _raw_frole.strip()
                or _raw_frole not in _allowed_tier6b_roles
            ):
                return _hold_unqual(f'INVALID_TIER6B_FRACTION_ROLE:{fid}')
            _explicit_frole = _raw_frole

        if fid not in agg_by_frac:
            agg_by_frac[fid] = {
                'N': n_raw,
                'k': k_raw,
                'plated_treated_cells_sum': plated_n,
                'presort_fraction_cells': ps_frac,
                'declared_clonogenic_opportunities': [],
                'sample_ids': [canon_tok],
                'records': [rec],
                'obs_identities': [(rec, list(aliquot_tokens))],
                'cal_by_key': {cal_key_obs: cal_tuple_obs},
                'untr_by_id': {untr_sid: (untr_k, untr_n)},
                'obs_list': [(n_raw, plated_n, untr_k, untr_n)],
                'fraction_roles': ({_explicit_frole} if _explicit_frole is not None else set()),
                'rec': rec,
            }
        else:
            prev_r = agg_by_frac[fid]['rec']
            if _explicit_frole is not None:
                agg_by_frac[fid]['fraction_roles'].add(_explicit_frole)
                if len(agg_by_frac[fid]['fraction_roles']) > 1:
                    return _hold_unqual(f'CONFLICTING_TIER6B_FRACTION_ROLE:{fid}')
            if (
                _tier6b_signature(prev_r.get('w_i')) != _tier6b_signature(rec.get('w_i'))
                or agg_by_frac[fid]['presort_fraction_cells'] != ps_frac
            ):
                return _hold_unqual(f'CONFLICTING_WEIGHTS_OR_PRESORT_CELLS_WITHIN_FRACTION:{fid}')
            if cal_key_obs in agg_by_frac[fid]['cal_by_key']:
                if agg_by_frac[fid]['cal_by_key'][cal_key_obs] != cal_tuple_obs:
                    return _hold_unqual(f'CONFLICTING_CALIBRATION_COUNTS_FOR_SAME_CAL_BATCH_ID:{fid}')
            else:
                agg_by_frac[fid]['cal_by_key'][cal_key_obs] = cal_tuple_obs
            if untr_sid in agg_by_frac[fid]['untr_by_id']:
                if agg_by_frac[fid]['untr_by_id'][untr_sid] != (untr_k, untr_n):
                    return _hold_unqual(f'CONFLICTING_UNTREATED_CFU_CONTROL_COUNTS:{untr_sid}')
            else:
                agg_by_frac[fid]['untr_by_id'][untr_sid] = (untr_k, untr_n)
            agg_by_frac[fid]['N'] += n_raw
            agg_by_frac[fid]['plated_treated_cells_sum'] += plated_n
            agg_by_frac[fid]['sample_ids'].append(canon_tok)
            agg_by_frac[fid]['records'].append(rec)
            agg_by_frac[fid]['obs_identities'].append((rec, list(aliquot_tokens)))
            agg_by_frac[fid]['obs_list'].append((n_raw, plated_n, untr_k, untr_n))
        if rec.get('presort_clonogenic_opportunities') is not None:
            agg_by_frac[fid]['declared_clonogenic_opportunities'].append(rec['presort_clonogenic_opportunities'])

    if lot_conflicts:
        return _hold_unqual(lot_conflicts[0])

    if not all_untreated_sample_ids.isdisjoint(all_treated_sample_ids):
        return _hold_unqual('NON_DISJOINT_UNTREATED_CFU_SAMPLE_ID_ACROSS_LOT')
    if not all_calibrant_specimen_ids.isdisjoint(all_treated_sample_ids) or not all_calibrant_specimen_ids.isdisjoint(all_untreated_sample_ids):
        return _hold_unqual('NON_DISJOINT_CALIBRANT_SPECIMEN_ID_ACROSS_LOT')

    if lot_roster is None or lot_presort_total_cells is None or lot_reconciled_loss_cells is None:
        return _hold_unqual('MISSING_SORTED_LOT_TIER6B_EVIDENCE')
    m_fracs = len(lot_roster)

    def _norm_frac_label(lbl: str) -> tuple[str, list[str], set[str]]:
        words = re.findall(r'[a-z0-9]+', lbl.lower())
        norm = '_'.join(words)
        return norm, words, set(words)

    _negation_words = frozenset({
        'not', 'non', 'no', 'un', 'without', 'anti', 'never', 'neither', 'except', 'minus', 'negative',
        'mixed', 'hybrid', 'pooled', 'pool', 'unsorted', 'unfractionated', 'bulk', 'crude', 'raw',
        'whole', 'total', 'impure', 'contaminated', 'contaminant', 'unknown', 'unclassified',
        'ambiguous', 'uncertain', 'partial', 'other', 'bypass', 'skipped', 'mock', 'blank',
        'reverse', 'inverted', 'failed', 'unpurified', 'unselected', 'pseudo', 'quasi',
    })
    _negated_substrings = (
        'unenriched', 'nonenriched', 'notenriched',
        'unflowthrough', 'nonflowthrough', 'notflowthrough',
        'neg_enriched', 'enriched_neg',
        'neg_flow_through', 'flow_through_neg',
        'neg_flowthrough', 'flowthrough_neg',
        'neg_ft', 'ft_neg',
        'unsortedenriched', 'mixedenriched', 'bulkenriched',
    )
    _ft_or_excl_words = frozenset({
        'flowthrough', 'ft', 'excluded', 'exclude', 'depleted', 'deplete', 'depletion',
        'discard', 'discarded', 'waste', 'rejected', 'reject', 'residual', 'remainder',
        'remaining', 'leftover', 'removed', 'stripped', 'unbound', 'unretained',
        'wash', 'eluate', 'eluted', 'sidestream', 'nontarget', 'offtarget',
    })
    _hybrid_concat_substrings = (
        'depletedenriched', 'enricheddepleted',
        'excludedenriched', 'enrichedexcluded',
        'discardenriched', 'enricheddiscard',
        'wasteenriched', 'enrichedwaste',
        'flowthroughenriched', 'enrichedflowthrough',
        'rejectedenriched', 'enrichedrejected',
    )

    def _is_negated_lbl(lbl: str) -> bool:
        norm, words, word_set = _norm_frac_label(lbl)
        if word_set.intersection(_negation_words):
            return True
        if any(sub in norm for sub in _negated_substrings):
            return True
        for i in range(len(words) - 1):
            if words[i] == 'neg' and words[i + 1] in ('enriched', 'flow', 'flowthrough', 'ft'):
                return True
            if words[i] in ('enriched', 'through', 'flowthrough', 'ft') and words[i + 1] == 'neg':
                return True
        return any(
            w.startswith(('non', 'not', 'un')) and any(w[len(p):].startswith(('enriched', 'flow', 'ft')) for p in ('non', 'not', 'un') if w.startswith(p))
            for w in words
        )

    if any(_is_negated_lbl(f) for f in lot_roster):
        return _hold_unqual('AMBIGUOUS_OR_NEGATED_FRACTION_ID')

    def _has_enr_token(lbl: str) -> bool:
        norm, _words, word_set = _norm_frac_label(lbl)
        return ('enriched' in word_set) or ('enriched' in norm)

    def _is_ft_lbl(lbl: str) -> bool:
        _norm, words, word_set = _norm_frac_label(lbl)
        return (
            ('flowthrough' in word_set)
            or ('ft' in word_set)
            or any(words[i] == 'flow' and words[i + 1] == 'through' for i in range(len(words) - 1))
        )

    def _has_ft_or_excl_token(lbl: str) -> bool:
        norm, words, word_set = _norm_frac_label(lbl)
        return (
            _is_ft_lbl(lbl)
            or bool(word_set.intersection(_ft_or_excl_words))
            or any(sub in norm for sub in _hybrid_concat_substrings)
            or any(
                (words[i], words[i + 1]) in (('neg', 'sort'), ('side', 'stream'), ('non', 'target'), ('off', 'target'))
                for i in range(len(words) - 1)
            )
        )

    # Strict closed positive token grammar for unambiguous sorted-lot ENRICHED vs FLOW_THROUGH_OR_EXCLUDED roles:
    _enr_neg_markers = frozenset({'cd31', 'cd144', 'cd34', 'pecam1', 'cdh5', 'kdr', 'vegfr2', 'cd276', 'b7h3', 'cd45'})
    _enr_allowed_qualifiers = frozenset({
        'fraction', 'subfraction', 'arm', 'part', 'sort', 'sorted', 'stage7',
        'endocrine', 'islet', 'beta', 'pure', 'target',
    })
    _ft_excl_allowed_tokens = frozenset({
        'flow', 'through', 'flowthrough', 'ft', 'and', 'or',
        'excluded', 'exclude', 'discard', 'discarded', 'waste', 'depleted', 'deplete',
        'fraction', 'subfraction', 'arm', 'part', 'sort', 'sorted',
        'neg', 'pos', 'cd31', 'cd144', 'cd34', 'pecam1', 'cdh5', 'kdr', 'vegfr2', 'cd45',
        'vascular', 'endothelial', 'stage7', 'pre', 'assembly',
    })

    def _is_strict_enr_lbl(lbl: str) -> bool:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', lbl):
            return False
        _norm, words, _word_set = _norm_frac_label(lbl)
        if not words or words[0] != 'enriched' or _has_ft_or_excl_token(lbl):
            return False
        idx_w = 1
        while idx_w < len(words):
            w = words[idx_w]
            if w in _enr_neg_markers and idx_w + 1 < len(words) and words[idx_w + 1] == 'neg':
                idx_w += 2
                continue
            if w == 'pre' and idx_w + 1 < len(words) and words[idx_w + 1] == 'assembly':
                idx_w += 2
                continue
            if w in _enr_allowed_qualifiers or re.fullmatch(r'[0-9]+', w):
                idx_w += 1
                continue
            return False
        return True

    def _is_strict_ft_or_excl_lbl(lbl: str) -> bool:
        if not re.fullmatch(r'[A-Za-z0-9_-]+', lbl):
            return False
        _norm, words, _word_set = _norm_frac_label(lbl)
        if not words or _has_enr_token(lbl):
            return False
        starts_ft_excl = (
            words[0] in ('flowthrough', 'ft', 'excluded', 'exclude', 'discard', 'discarded', 'waste', 'depleted', 'deplete')
            or (len(words) >= 2 and (words[0], words[1]) in (('flow', 'through'), ('neg', 'sort')))
        )
        if not starts_ft_excl:
            return False
        return all((w in _ft_excl_allowed_tokens or re.fullmatch(r'[0-9]+', w)) for w in words)

    hybrid_fracs = [f for f in lot_roster if _has_enr_token(f) and _has_ft_or_excl_token(f)]
    if hybrid_fracs:
        return _hold_unqual('AMBIGUOUS_HYBRID_ENRICHED_FLOW_THROUGH_FRACTION_ID')
    # Every declared fraction in the roster must classify into an unambiguous strict sorted-lot role (#26, T6; Patch v11):
    unclassified_fracs = [f for f in lot_roster if not _is_strict_enr_lbl(f) and not _is_strict_ft_or_excl_lbl(f)]
    if unclassified_fracs:
        return _hold_unqual(f'UNCLASSIFIED_OR_AMBIGUOUS_FRACTION_ROLE:{",".join(unclassified_fracs)}')
    pure_enriched = [f for f in lot_roster if _is_strict_enr_lbl(f)]
    pure_ft = [f for f in lot_roster if _is_strict_ft_or_excl_lbl(f) and _is_ft_lbl(f)]
    if m_fracs < 2 or not pure_enriched or not pure_ft or set(pure_enriched).intersection(pure_ft):
        return _hold_unqual('OMITTED_ENRICHED_OR_FLOW_THROUGH_FRACTION_IN_SORTED_LOT')
    missing_fracs = [f for f in lot_roster if f not in agg_by_frac]
    if missing_fracs:
        return _hold_unqual(f'MISSING_DECLARED_FRACTION_EVIDENCE:{",".join(missing_fracs)}')
    # Cross-check optional explicit typed `fraction_role` against the unambiguous label classification:
    for f in lot_roster:
        declared_roles = agg_by_frac[f]['fraction_roles']
        if declared_roles:
            expected_role = 'ENRICHED' if _is_strict_enr_lbl(f) else 'FLOW_THROUGH_OR_EXCLUDED'
            if declared_roles != {expected_role}:
                return _hold_unqual(f'CONFLICTING_TIER6B_FRACTION_ROLE_AND_LABEL:{f}')

    # Enforce exact integer pre-sort mass balance & disjoint treated + untreated cell budget <= fraction physical cells:
    total_sorted_cells = sum(int(agg_by_frac[f]['presort_fraction_cells']) for f in lot_roster)
    if total_sorted_cells + lot_reconciled_loss_cells != lot_presort_total_cells:
        return _hold_unqual('PRESORT_INTEGER_MASS_BALANCE_MISMATCH')
    for f in lot_roster:
        untr_n_f_sum = sum(n_u for _k_u, n_u in agg_by_frac[f]['untr_by_id'].values())
        if int(agg_by_frac[f]['plated_treated_cells_sum']) > int(agg_by_frac[f]['presort_fraction_cells']):
            return _hold_unqual(f'PLATED_CELLS_EXCEED_FRACTION_PHYSICAL_CELLS:{f}')
        if int(agg_by_frac[f]['plated_treated_cells_sum']) + untr_n_f_sum > int(agg_by_frac[f]['presort_fraction_cells']):
            return _hold_unqual(f'PLATED_AND_UNTREATED_CELLS_EXCEED_FRACTION_PHYSICAL_CELLS:{f}')

    # Verify point-estimate pre-sort clonogenic-risk-unit shares w_i = E_clono_i / sum(E_clono_j) if w_i is declared:
    e_clono_by_frac: dict[str, Fraction] = {}
    for f in lot_roster:
        uk_sum = sum(k_u for k_u, _n_u in agg_by_frac[f]['untr_by_id'].values())
        un_sum = sum(n_u for _k_u, n_u in agg_by_frac[f]['untr_by_id'].values())
        e_clono_f = Fraction(int(agg_by_frac[f]['presort_fraction_cells']), 1) * Fraction(uk_sum, un_sum)
        for sup_e in agg_by_frac[f]['declared_clonogenic_opportunities']:
            if not _is_strict_finite_number_iip0024(sup_e) or float(sup_e) <= 0.0:
                return _hold_unqual(f'INVALID_PRESORT_CLONOGENIC_OPPORTUNITIES:{f}')
            if abs(float(sup_e) - float(e_clono_f)) > max(1.0, 1e-6 * float(e_clono_f)):
                return _hold_unqual(f'PRESORT_CLONOGENIC_OPPORTUNITIES_MISMATCH_WITH_UNTREATED_CFU:{f}')
        e_clono_by_frac[f] = e_clono_f
    total_e_clono = sum(e_clono_by_frac.values(), Fraction(0, 1))
    if total_e_clono <= 0:
        return _hold_unqual('ZERO_TOTAL_PRE_SORT_CLONOGENIC_RISK_UNITS')

    w_point: dict[str, float] = {}
    for f in lot_roster:
        w_verified_f = float(e_clono_by_frac[f] / total_e_clono)
        r_f = agg_by_frac[f]['rec']
        if 'w_i' in r_f:
            w_raw = r_f.get('w_i')
            if not _is_strict_finite_number_iip0024(w_raw):
                return _hold_unqual(f'INVALID_OR_MISSING_FRACTION_WEIGHT:{f}')
            w_v = float(w_raw)
            if not (0.0 < w_v <= 1.0):
                return _hold_unqual(f'OUT_OF_RANGE_FRACTION_WEIGHT:{f}')
            if abs(w_v - w_verified_f) > 1e-9:
                return _hold_unqual(f'FRACTION_WEIGHT_MISMATCH_WITH_VERIFIED_CLONOGENIC_RISK_UNIT_SHARE:{f}')
        w_point[f] = w_verified_f
    if abs(sum(w_point.values()) - 1.0) > 1e-9:
        return _hold_unqual('INVALID_FRACTION_WEIGHT_SUM_NOT_EQUAL_TO_1')

    # UNCONDITIONAL 3-way simultaneous Clopper-Pearson finite-sample 95% coverage across ALL m fractions:
    alpha_obs = 0.05 / (3.0 * float(m_fracs))
    alpha_cal = 0.05 / (3.0 * float(m_fracs))
    alpha_cfu = 0.05 / (3.0 * float(m_fracs))
    alpha_each = alpha_cal / 2.0

    d_L_by_frac: dict[str, float] = {}
    d_L_details: dict[str, dict[str, Any]] = {}
    p_i_U_by_frac: dict[str, float] = {}
    N_by_frac: dict[str, int] = {}
    e_lo_by_frac: dict[str, float] = {}
    e_hi_by_frac: dict[str, float] = {}

    for f in sorted(lot_roster):
        N_i = int(agg_by_frac[f]['N'])
        plated_sum_f = int(agg_by_frac[f]['plated_treated_cells_sum'])
        untr_items = sorted(agg_by_frac[f]['untr_by_id'].items())
        u_f_cnt = len(untr_items)
        alpha_cfu_tail = alpha_cfu / (2.0 * float(u_f_cnt))
        f_lo_vals = [_exact_one_sided_cp_lower_iip0031(uk_v, un_v, alpha_cfu_tail) for _, (uk_v, un_v) in untr_items]
        f_hi_vals = [clopper_pearson_upper95_independent(uk_v, un_v, alpha=alpha_cfu_tail) for _, (uk_v, un_v) in untr_items]
        f_lo_f = min(f_lo_vals)
        f_hi_f = max(f_hi_vals)
        for (n_sub, plated_sub, uk_sub, un_sub) in agg_by_frac[f]['obs_list']:
            f_sub_l = _exact_one_sided_cp_lower_iip0031(uk_sub, un_sub, alpha_cfu_tail)
            max_sup_sub = _conservative_int_product_floor_iip0031(plated_sub, f_sub_l)
            if n_sub > max_sup_sub:
                return _hold_unqual(
                    'UNTREATED_CFU_DOES_NOT_SUPPORT_CLAIMED_OPPORTUNITIES',
                    untreated_cfu_max_supported_N=max_sup_sub,
                    claimed_N=n_sub,
                    fraction_id=f,
                )
        max_sup_frac = _conservative_int_product_floor_iip0031(plated_sum_f, f_lo_f)
        if N_i > max_sup_frac:
            return _hold_unqual(
                'UNTREATED_CFU_DOES_NOT_SUPPORT_CLAIMED_OPPORTUNITIES',
                untreated_cfu_max_supported_N=max_sup_frac,
                claimed_N=N_i,
                fraction_id=f,
            )
        c_f = float(agg_by_frac[f]['presort_fraction_cells'])
        e_lo_by_frac[f] = c_f * f_lo_f
        e_hi_by_frac[f] = c_f * f_hi_f
        N_by_frac[f] = N_i

        cal_items = sorted(agg_by_frac[f]['cal_by_key'].items())
        b_i_count = len(cal_items)
        alpha_each_f = alpha_cal / (2.0 * float(b_i_count))
        per_cal_bounds = []
        for (_ck, (rk_h, rn_h, rk_d, rn_d)) in cal_items:
            d_l_h = _exact_one_sided_cp_lower_iip0031(rk_h, rn_h, alpha_each_f)
            d_l_d = _exact_one_sided_cp_lower_iip0031(rk_d, rn_d, alpha_each_f)
            per_cal_bounds.append((min(d_l_h, d_l_d), rk_h, rn_h, d_l_h, rk_d, rn_d, d_l_d))
        worst_cal = min(per_cal_bounds, key=lambda x: x[0])
        d_i_L = worst_cal[0]
        d_L_details[f] = {
            'mode': 'EXACT_DUAL_SEPARATE_CALIBRANTS',
            'distinct_calibration_batches_in_fraction': b_i_count,
            'alpha_per_calibrant_in_fraction': alpha_each_f,
            'r_k_hpsc': worst_cal[1],
            'r_n_hpsc': worst_cal[2],
            'd_L_hpsc': worst_cal[3],
            'r_k_diff': worst_cal[4],
            'r_n_diff': worst_cal[5],
            'd_L_diff': worst_cal[6],
            'd_i_L': d_i_L,
            'f_cfu_L': f_lo_f,
            'f_cfu_U': f_hi_f,
        }
        d_L_by_frac[f] = d_i_L
        if d_i_L > 0.0 and f_lo_f > 0.0:
            p_exact_plated = min(1.0, (-math.expm1(math.log(alpha_obs) / float(plated_sum_f))) / (f_lo_f * d_i_L))
            p_decl_n = min(1.0, (-math.expm1(math.log(alpha_obs) / float(N_i))) / d_i_L)
            p_i_U_by_frac[f] = max(p_exact_plated, p_decl_n)
        else:
            p_i_U_by_frac[f] = 1.0

    # Exact global maximizer of sum_i (E_i / sum_j E_j) * p_{i,U} over the simultaneous box E_i in [E_lo_i, E_hi_i]:
    sorted_by_p = sorted(lot_roster, key=lambda fk: (-p_i_U_by_frac[fk], fk))
    best_p_lot = -1.0
    w_used: dict[str, float] = dict(w_point)
    for r_cut in range(len(sorted_by_p) + 1):
        e_vec = {
            fk: (e_hi_by_frac[fk] if idx_r < r_cut else e_lo_by_frac[fk])
            for idx_r, fk in enumerate(sorted_by_p)
        }
        tot_e_vec = sum(e_vec.values())
        if tot_e_vec <= 0.0:
            continue
        cand_w = {fk: e_vec[fk] / tot_e_vec for fk in lot_roster}
        cand_p = float(sum(cand_w[fk] * p_i_U_by_frac[fk] for fk in lot_roster))
        if cand_p > best_p_lot:
            best_p_lot = cand_p
            w_used = cand_w
    p_lot_U95 = best_p_lot if best_p_lot >= 0.0 else 1.0

    detection_floor_cleared = all(d_L_by_frac[f] >= TIER6B_DETECTION_FLOOR_IIP0031 for f in lot_roster)
    if not detection_floor_cleared:
        return _hold_unqual(
            'DETECTION_LOWER_BOUND_BELOW_0P50_FLOOR',
            p_lot_U95=p_lot_U95,
            d_L_by_fraction=d_L_by_frac,
            d_L_details=d_L_details,
            p_i_U_by_fraction=p_i_U_by_frac,
            N_by_fraction=N_by_frac,
            w_by_fraction=w_used,
            alpha_obs=alpha_obs,
            alpha_cal=alpha_cal,
            alpha_each=alpha_each,
            m_fractions=m_fracs,
            detection_floor_cleared=False,
            deduplicated_sample_ids=sorted(deduplicated_sample_ids),
        )

    # Exact unrounded comparison against 3.0e-7 (exact equality passes; never round before comparing):
    if p_lot_U95 > TIER6B_PER_FRACTION_LOT_CAP_IIP0031:
        return {
            'gate_pass': False,
            'status': 'HOLD_TIER6B_INSUFFICIENT_DENOMINATOR',
            'disposition': 'HOLD_TIER6B_INSUFFICIENT_DENOMINATOR',
            'reason_code': 'HOLD_TIER6B_LOT_U95_EXCEEDS_3E_7',
            'tier6b_status': 'HOLD_TIER6B_INSUFFICIENT_DENOMINATOR',
            'hard_stop': False,
            'p_lot_U95': p_lot_U95,
            'threshold': TIER6B_PER_FRACTION_LOT_CAP_IIP0031,
            'd_L_by_fraction': d_L_by_frac,
            'd_L_details': d_L_details,
            'p_i_U_by_fraction': p_i_U_by_frac,
            'N_by_fraction': N_by_frac,
            'w_by_fraction': w_used,
            'alpha_obs': alpha_obs,
            'alpha_cal': alpha_cal,
            'alpha_each': alpha_each,
            'm_fractions': m_fracs,
            'detection_floor_cleared': True,
            'deduplicated_sample_ids': sorted(deduplicated_sample_ids),
        }

    # Summing opportunities across >= 2 observations in the same fraction requires explicit, distinct
    # source-owned trial & calibration identifiers on every summed row so relabelled duplicate payloads
    # without stable source identity can NEVER pass Tier 6B (#5):
    for f in sorted(lot_roster):
        f_obs_ids = agg_by_frac[f]['obs_identities']
        if len(f_obs_ids) >= 2:
            seen_rids: set[str] = set()
            seen_srcs: set[str] = set()
            seen_cals: set[str] = set()
            for r_obs, r_aliquot_tokens in f_obs_ids:
                rid_val = next((r_obs[rk] for rk in _TIER6B_RECORD_ID_ALIASES if rk in r_obs and isinstance(r_obs[rk], str) and r_obs[rk].strip()), None)
                src_tokens = [tok for tok in r_aliquot_tokens[1:] if isinstance(tok, str) and tok.strip()]
                cal_val = r_obs.get('cal_batch_id')
                if (
                    rid_val is None or rid_val in seen_rids
                    or not src_tokens or any(tok in seen_srcs for tok in src_tokens)
                    or not isinstance(cal_val, str) or not cal_val.strip() or cal_val in seen_cals
                ):
                    return _hold_unqual(
                        'MISSING_SOURCE_OWNED_TRIAL_OR_CALIBRANT_IDENTITY_FOR_MULTI_OBSERVATION_SUM',
                        p_lot_U95=p_lot_U95,
                        d_L_by_fraction=d_L_by_frac,
                        d_L_details=d_L_details,
                        p_i_U_by_fraction=p_i_U_by_frac,
                        N_by_fraction=N_by_frac,
                        w_by_fraction=w_used,
                        alpha_obs=alpha_obs,
                        alpha_cal=alpha_cal,
                        alpha_each=alpha_each,
                        m_fractions=m_fracs,
                        detection_floor_cleared=True,
                        deduplicated_sample_ids=sorted(deduplicated_sample_ids),
                    )
                seen_rids.add(rid_val)
                seen_srcs.update(src_tokens)
                seen_cals.add(cal_val)

    return {
        'gate_pass': True,
        'status': 'PASS',
        'disposition': 'PASS',
        'reason_code': 'PASS_TIER6B_PER_FRACTION_THINNING',
        'tier6b_status': 'PASS_TIER6B_PER_FRACTION_THINNING',
        'hard_stop': False,
        'p_lot_U95': p_lot_U95,
        'threshold': TIER6B_PER_FRACTION_LOT_CAP_IIP0031,
        'd_L_by_fraction': d_L_by_frac,
        'd_L_details': d_L_details,
        'p_i_U_by_fraction': p_i_U_by_frac,
        'N_by_fraction': N_by_frac,
        'w_by_fraction': w_used,
        'alpha_obs': alpha_obs,
        'alpha_cal': alpha_cal,
        'alpha_each': alpha_each,
        'm_fractions': m_fracs,
        'detection_floor_cleared': True,
        'deduplicated_sample_ids': sorted(deduplicated_sample_ids),
    }


def _paired_difference_se_consistent(arm_a: Any, arm_b: Any, difference: Any) -> bool:
    """Check the sharp paired-SE bounds using the decimal values carried by JSON."""
    values = (arm_a, arm_b, difference)
    if not all(_is_strict_finite_number_iip0024(v) and v >= 0 for v in values):
        return False
    a, b, d = (Fraction(str(v)) for v in values)
    return abs(a - b) <= d <= a + b


_ALL_S11_RECOGNIZED_EVIDENCE_KEYS: frozenset[str] = frozenset({
    'ciita_locus2_genotype_copy_class',
    'ciita_ki_allele_informative_transcript_recovery_verified',
    'ciita_ko_allele_informative_transcript_recovery_verified',
    'ciita_ki_allele1_informative_transcript_recovery_verified',
    'ciita_ki_allele2_informative_transcript_recovery_verified',
    'ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified',
    'lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified',
    'phasing_verified',
    'cis_verified',
    'wt_inducibility_verified',
    'rnai_attribution_verified',
    'rnai_specific_benefit_claimed',
    'rnai_attribution_scramble_control_id',
    'rnai_attribution_scramble_background_id',
    'rnai_attribution_scramble_control_sha256',
    'rnai_attribution_scramble_sample_id',
    'rnai_attribution_intronless_control_id',
    'rnai_attribution_intronless_background_id',
    'rnai_attribution_intronless_control_sha256',
    'rnai_attribution_intronless_sample_id',
    'l2_v3_same_background_intronless_control_linked_to_cargo',
    'l2_v3_same_background_scramble_control_linked_to_cargo',
    'has_intronic_mir_cd276_pvr',
    'syn_5xisre_cd200_in_design_y',
    's11_replicate_n',
    'obs_edited_72h_ifng_hla_dr_dp_dq_pos_pct',
    'obs_edited_72h_ifng_hla_dr_dp_dq_pos_se',
    'obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_pct',
    'obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_se',
    'obs_cis_socs1_tost90_ratio_lo',
    'obs_cis_socs1_tost90_ratio_hi',
    'obs_cis_dexi_tost90_ratio_lo',
    'obs_cis_dexi_tost90_ratio_hi',
    'obs_cis_clec16a_tost90_ratio_lo',
    'obs_cis_clec16a_tost90_ratio_hi',
    'obs_16p13_ddpcr_cn_ratio',
    'obs_16p13_phased_baf',
})

_S11_KIKO_ONLY_FIELDS: frozenset[str] = frozenset({
    'obs_ciita_ki_allele_informative_reads_n',
    'obs_ciita_ki_ex4_ex6_aberrant_skip_reads_k',
    'obs_ciita_ki_ex4_ex6_aberrant_skip_pct',
    'obs_ciita_ki_ex4_ex6_aberrant_skip_se',
    'obs_ciita_ko_allele_informative_reads_n',
    'obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k',
    'obs_ciita_ko_ex4_ex6_aberrant_skip_pct',
    'obs_ciita_ko_ex4_ex6_aberrant_skip_se',
})

_S11_KIKI_ONLY_FIELDS: frozenset[str] = frozenset({
    'obs_ciita_ki_allele1_informative_reads_n',
    'obs_ciita_ki_allele2_informative_reads_n',
    'obs_ciita_allele1_informative_reads_n',
    'obs_ciita_allele2_informative_reads_n',
    'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_reads_k',
    'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_reads_k',
    'obs_ciita_allele1_ex4_ex6_aberrant_skip_reads_k',
    'obs_ciita_allele2_ex4_ex6_aberrant_skip_reads_k',
    'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_pct',
    'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_pct',
    'obs_ciita_allele1_ex4_ex6_aberrant_skip_pct',
    'obs_ciita_allele2_ex4_ex6_aberrant_skip_pct',
    'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_se',
    'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_se',
    'obs_ciita_allele1_ex4_ex6_aberrant_skip_se',
    'obs_ciita_allele2_ex4_ex6_aberrant_skip_se',
})


def evaluate_locus2_ciita_modeb_challenger_v3_gate(
    payload: dict[str, Any],
    *,
    root: Path | str | None = None,
    product_def: dict[str, Any] | None = None,
    _trusted_p_mod: Any = None,
) -> dict[str, Any]:
    """
    Fail-closed Locus-2 CIITA Mode-B v3 & Rule S11 evaluator (IIP-0031; DEC-110).
    """
    root_path = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    p_mod = _trusted_p_mod if _trusted_p_mod is not None else _load_product_contract_module(root_path)

    payload_malformed = not isinstance(payload, dict)
    if payload_malformed:
        payload = {}
    raw_rec_obj = payload.get('iip0031_locus2_record')
    rec_malformed = payload_malformed or (raw_rec_obj is not None and not isinstance(raw_rec_obj, dict))
    rec_in: dict[str, Any] = dict(raw_rec_obj) if isinstance(raw_rec_obj, dict) else {}

    outer_only = {k: v for k, v in payload.items() if isinstance(k, str) and k != 'iip0031_locus2_record' and v is not None}
    pvid_supplied = ('product_variant_id' in payload) or ('product_variant_id' in rec_in)
    pvid_outer_raw = payload.get('product_variant_id') if 'product_variant_id' in payload else rec_in.get('product_variant_id')
    if not pvid_supplied and isinstance(product_def, dict) and _trusted_p_mod is not None:
        pvid_outer_raw = product_def.get('product_variant_id')
        pvid_supplied = True
    pvid_canon, pdef_loaded, registry_ok = _canonicalize_product_id_detailed_iip0031(
        pvid_outer_raw,
        root_path,
        _trusted_p_mod=p_mod,
    )

    # Never let an unverified caller `product_def` override the authoritative registry binding (#14):
    if pdef_loaded is not None:
        if isinstance(product_def, dict) and hasattr(p_mod, 'validate_product_definition_dict'):
            try:
                pdef = p_mod.validate_product_definition_dict(product_def, root=root_path)
                if pdef.get('product_variant_id') != pdef_loaded.get('product_variant_id'):
                    pdef = pdef_loaded
            except ValueError:
                pdef = pdef_loaded
        else:
            pdef = pdef_loaded
    else:
        pdef = None

    l2_binding = pdef.get('staged_locus2_modeb_donor_binding') if isinstance(pdef, dict) else None
    expected_donor_id = str(l2_binding.get('locus2_donor_id')) if isinstance(l2_binding, dict) and l2_binding.get('locus2_donor_id') else None
    expected_donor_sha = str(l2_binding.get('locus2_donor_sha256')) if isinstance(l2_binding, dict) and l2_binding.get('locus2_donor_sha256') else None
    expected_bg_id = str(l2_binding.get('background_id')) if isinstance(l2_binding, dict) and l2_binding.get('background_id') else None
    supplied_donor_ids = [
        v for v in (
            outer_only.get('l2_v3_donor_id'),
            outer_only.get('donor_id'),
            rec_in.get('l2_v3_donor_id'),
            rec_in.get('donor_id'),
        )
        if v is not None
    ]
    supplied_bg_ids = [
        v for v in (outer_only.get('background_id'), rec_in.get('background_id'))
        if v is not None
    ]
    donor_id = (
        expected_donor_id
        if expected_donor_id is not None
        else (supplied_donor_ids[0] if supplied_donor_ids and isinstance(supplied_donor_ids[0], str) else None)
    )
    _all_candidate_donors = {
        d for d in ([expected_donor_id] + supplied_donor_ids) if isinstance(d, str) and d.strip()
    }
    _all_candidate_bgs = [
        b.strip() for b in ([expected_bg_id] + supplied_bg_ids) if isinstance(b, str) and b.strip()
    ]
    _bg_outer = outer_only.get('background_id')
    _bg_nested = rec_in.get('background_id')
    bg_id_raw = (
        _bg_outer.strip() if isinstance(_bg_outer, str) and _bg_outer.strip()
        else (expected_bg_id or (_bg_nested.strip() if isinstance(_bg_nested, str) else ''))
    )
    bg_id = expected_bg_id if expected_bg_id is not None else bg_id_raw
    is_l2_variant = bool(pvid_canon in getattr(p_mod, 'STAGED_LOCUS2_V3_PRODUCT_IDS', frozenset())) or bool(donor_id) or bool(supplied_donor_ids)
    merged = {**rec_in, **outer_only}

    def _out(status: str, gate_pass: bool, reasons: list[str], **extra: Any) -> dict[str, Any]:
        return {
            'rule_id': RULE_S11_GATE_ID_IIP0031,
            'decision_id': 'DEC-110',
            'contribution_id': 'IIP-0031',
            'product_variant_id': pvid_canon,
            'l2_v3_donor_id': donor_id,
            'status': status,
            'gate_pass': gate_pass,
            'reasons': reasons,
            'execution_supported': False,
            **(
                {'barrier_intervention_classification': IIP0031_BARRIER_CHS4_250_CLASSIFICATION}
                if donor_id == 'L2_V3_BARRIER_CHS4_250' or pvid_canon == 'ISLET1_A4_L2_V3_BARRIER_CHS4_250'
                else {}
            ),
            **extra,
        }

    if not registry_ok:
        return _out(
            'HOLD_S11_UNVERIFIED_OR_CORRUPTED_PRODUCT_REGISTRY',
            False,
            ['HOLD_S11_UNVERIFIED_OR_CORRUPTED_PRODUCT_REGISTRY'],
        )

    has_any_s11_fields = bool(rec_in) or any(
        isinstance(k, str) and (
            k.startswith(('obs_l2_v3_', 'obs_ciita_'))
            or k in _ALL_S11_RECOGNIZED_EVIDENCE_KEYS
        )
        for k in outer_only
    )

    # 1. Hard architectural compatibility blocks (evaluated uniformly across flat & nested inputs, #12, #14):
    if isinstance(pvid_canon, str) and pvid_canon.startswith('ISLET1_B_SCT_') and (
        donor_id is not None
        or bool(supplied_donor_ids)
        or has_any_s11_fields
        or rec_malformed
        or merged.get('require_full_iip0031_locus2_and_tier6b_panel') is True
    ):
        return _out(
            'FAIL_S11_TRACK2_BSCT_INCOMPATIBLE_AT_LOCUS_2',
            False,
            ['FAIL_S11_TRACK2_BSCT_INCOMPATIBLE_AT_LOCUS_2', 'INCOMPATIBLE_AT_LOCUS_2'],
        )

    is_design_y_any = bool(
        (isinstance(pvid_canon, str) and pvid_canon.startswith('ISLET1_Y_'))
        or (isinstance(pvid_outer_raw, str) and pvid_outer_raw == pvid_outer_raw.strip() and pvid_outer_raw.startswith('ISLET1_Y_'))
        or any(b.upper().startswith('DESIGN_Y') for b in _all_candidate_bgs)
    )
    has_syn5_cd200_any = bool(
        _all_candidate_donors.intersection({
            'L2_V3_CO_PRIMARY_CD200',
            'L2_V3_DESIGN_B_CD200',
            'L2_V3_CTRL_NT_SCRAMBLE',
            'L2_V3_BARRIER_CHS4_250',
        })
        or (isinstance(l2_binding, dict) and l2_binding.get('unit2_promoter') == 'Syn_5xISRE_minP' and l2_binding.get('unit2_cargo') == 'CD200')
        or any(c.get('syn_5xisre_cd200_in_design_y') is True for c in (outer_only, rec_in))
    )
    if is_design_y_any and has_syn5_cd200_any:
        return _out('FAIL_S11_SYN5XISRE_CD200_BARRED_IN_DESIGN_Y_Y03', False, ['FAIL_S11_SYN5XISRE_CD200_BARRED_IN_DESIGN_Y_Y03'])

    has_mir_active = bool(
        (isinstance(l2_binding, dict) and l2_binding.get('has_intronic_mir') is True and l2_binding.get('mir_class') == 'ACTIVE_CD276_PVR')
        or any(c.get('has_intronic_mir_cd276_pvr') is True for c in (outer_only, rec_in))
        or bool(_all_candidate_donors.intersection({
            'L2_V3_CO_PRIMARY_CD200',
            'L2_V3_CO_PRIMARY_MANF',
            'L2_V3_ALT2_CCL22',
            'L2_V3_BARRIER_CHS4_250',
            'L2_V3_CD274P_ZERO_TATATA',
        }))
    )
    is_design_b_bg = bool(
        (isinstance(pdef, dict) and pdef.get('family') == 'B' and not str(pvid_canon or '').startswith('ISLET1_B_SCT_'))
        or any(b.upper().startswith('DESIGN_B') for b in _all_candidate_bgs)
    )
    if is_design_b_bg and has_mir_active:
        return _out('FAIL_S11_MIR_CD276_PVR_BARRED_ON_DESIGN_B_HLAC_INTACT', False, ['FAIL_S11_MIR_CD276_PVR_BARRED_ON_DESIGN_B_HLAC_INTACT'])

    if pvid_supplied and pdef_loaded is None:
        return _out(
            'HOLD_S11_UNREGISTERED_OR_INVALID_PRODUCT_VARIANT_ID',
            False,
            ['HOLD_S11_UNREGISTERED_OR_INVALID_PRODUCT_VARIANT_ID'],
        )

    if not pvid_supplied and not is_l2_variant and not has_any_s11_fields and not rec_malformed:
        return _out(
            'HOLD_S11_UNREGISTERED_OR_INVALID_PRODUCT_VARIANT_ID',
            False,
            ['HOLD_S11_UNREGISTERED_OR_INVALID_PRODUCT_VARIANT_ID'],
        )

    if (
        pdef_loaded is not None
        and not is_l2_variant
        and not has_any_s11_fields
        and not rec_malformed
        and merged.get('require_full_iip0031_locus2_and_tier6b_panel') is not True
    ):
        return _out('NOT_APPLICABLE_NON_LOCUS2_MODEB_VARIANT', True, [])

    if rec_malformed:
        return _out('HOLD_S11_MALFORMED_RECORD', False, ['HOLD_S11_MALFORMED_RECORD'])

    # 2. Enforce LOT_SCOPE, strict Donor-ID + SHA-256 + Background Registry Binding, & forbid nested provenance overrides:
    for sc_check_v in (outer_only.get('scope_type'), rec_in.get('scope_type')):
        if sc_check_v is not None and sc_check_v != 'LOT_SCOPE':
            return _out('HOLD_S11_PANEL_REQUIRES_LOT_SCOPE', False, ['HOLD_S11_PANEL_REQUIRES_LOT_SCOPE'])

    if isinstance(l2_binding, dict):
        if expected_donor_id is not None:
            for sup_did in supplied_donor_ids:
                if not isinstance(sup_did, str) or sup_did != expected_donor_id:
                    return _out('HOLD_S11_DONOR_MISMATCH', False, ['HOLD_S11_DONOR_MISMATCH'])
        if expected_donor_sha is not None:
            for sup_sha in (outer_only.get('construct_sha256'), rec_in.get('construct_sha256')):
                if sup_sha is not None and sup_sha != expected_donor_sha:
                    return _out('HOLD_S11_DONOR_MISMATCH', False, ['HOLD_S11_DONOR_MISMATCH'])
        if expected_bg_id is not None:
            for sup_bg in supplied_bg_ids:
                if not isinstance(sup_bg, str) or sup_bg != expected_bg_id:
                    return _out('HOLD_S11_BACKGROUND_MISMATCH', False, ['HOLD_S11_BACKGROUND_MISMATCH'])

    if 'product_variant_id' in outer_only and 'product_variant_id' in rec_in and rec_in['product_variant_id'] is not None:
        nested_pvid_c, _, _ = _canonicalize_product_id_detailed_iip0031(rec_in['product_variant_id'], root_path, _trusted_p_mod=p_mod)
        if nested_pvid_c != pvid_canon:
            return _out('FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH', False, ['FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH'])
    for id_k in ('construct_sha256', 'background_id', 'donor_seedline_id', 'seed_line_id', 'mcb_id', 'lot_id', 'current_lot_id', 'sample_id', 'sample_unit_id', 'l2_v3_donor_id', 'donor_id'):
        if id_k in outer_only and id_k in rec_in and rec_in[id_k] is not None and rec_in[id_k] != outer_only[id_k]:
            return _out('FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH', False, ['FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH'])
    _all_s11_sids = [
        v for v in (
            outer_only.get('sample_id'),
            outer_only.get('sample_unit_id'),
            rec_in.get('sample_id'),
            rec_in.get('sample_unit_id'),
        )
        if v is not None
    ]
    if any(not isinstance(v, str) or not v.strip() or v != v.strip() for v in _all_s11_sids):
        return _out('HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING', False, ['HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING'])
    if len(set(_all_s11_sids)) > 1:
        return _out('FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH', False, ['FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH'])
    _all_s11_lids = [
        v for v in (
            outer_only.get('lot_id'),
            outer_only.get('current_lot_id'),
            rec_in.get('lot_id'),
            rec_in.get('current_lot_id'),
        )
        if v is not None
    ]
    if any(not isinstance(v, str) or not v.strip() or v != v.strip() for v in _all_s11_lids):
        return _out('HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING', False, ['HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING'])
    if len(set(_all_s11_lids)) > 1:
        return _out('FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH', False, ['FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH'])
    _all_s11_cids = [
        v for v in (
            outer_only.get('campaign_id'),
            outer_only.get('assay_campaign_id'),
            rec_in.get('campaign_id'),
            rec_in.get('assay_campaign_id'),
        )
        if v is not None
    ]
    if any(not isinstance(v, str) or not v.strip() or v != v.strip() for v in _all_s11_cids) or len(set(_all_s11_cids)) > 1:
        return _out('HOLD_S11_CONFLICTING_NESTED_EVIDENCE', False, ['HOLD_S11_CONFLICTING_NESTED_EVIDENCE'])
    _all_s11_seeds = [
        v for v in (
            outer_only.get('donor_seedline_id'),
            outer_only.get('seed_line_id'),
            rec_in.get('donor_seedline_id'),
            rec_in.get('seed_line_id'),
        )
        if v is not None
    ]
    if any(not isinstance(v, str) or not v.strip() or v != v.strip() for v in _all_s11_seeds):
        return _out('HOLD_S11_MALFORMED_RECORD', False, ['HOLD_S11_MALFORMED_RECORD'])
    if len(set(_all_s11_seeds)) > 1:
        return _out('FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH', False, ['FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH'])
    _all_s11_dids = [
        v for v in (
            outer_only.get('l2_v3_donor_id'),
            outer_only.get('donor_id'),
            rec_in.get('l2_v3_donor_id'),
            rec_in.get('donor_id'),
        )
        if v is not None
    ]
    if any(not isinstance(v, str) or not v.strip() or v != v.strip() for v in _all_s11_dids):
        return _out('HOLD_S11_MALFORMED_RECORD', False, ['HOLD_S11_MALFORMED_RECORD'])
    if len(set(_all_s11_dids)) > 1:
        return _out('FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH', False, ['FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH'])
    for _prov_k in ('construct_sha256', 'background_id', 'mcb_id', 'product_variant_id'):
        for _prov_v in (outer_only.get(_prov_k), rec_in.get(_prov_k)):
            if _prov_v is not None and (not isinstance(_prov_v, str) or not _prov_v.strip() or _prov_v != _prov_v.strip()):
                return _out('HOLD_S11_MALFORMED_RECORD', False, ['HOLD_S11_MALFORMED_RECORD'])

    # Omitted S11 evidence on an L2_V3_* comparator MUST NEVER clear:
    if not has_any_s11_fields:
        return _out('HOLD_S11_INCOMPLETE_EVIDENCE', False, ['HOLD_S11_INCOMPLETE_EVIDENCE'])

    if 'LOT_SCOPE' not in (payload.get('scope_type'), rec_in.get('scope_type')):
        return _out('HOLD_S11_PANEL_REQUIRES_LOT_SCOPE', False, ['HOLD_S11_PANEL_REQUIRES_LOT_SCOPE'])

    # Endpoint-granular outer/nested conflict detection (#4, s11_extra.py):
    # Preserve independently attributable quantitative FAIL_S11_* stops while flagging conflicting fields as HOLD_S11_CONFLICTING_NESTED_EVIDENCE.
    conflicting_overlap_keys: set[str] = set()
    for key in sorted((payload.keys() & rec_in.keys()) - {'product_variant_id'}):
        try:
            outer_value = json.dumps(payload[key], sort_keys=True, allow_nan=False)
            nested_value = json.dumps(rec_in[key], sort_keys=True, allow_nan=False)
        except (TypeError, ValueError, OverflowError):
            return _out('HOLD_S11_MALFORMED_RECORD', False, ['HOLD_S11_MALFORMED_RECORD'])
        if outer_value != nested_value:
            conflicting_overlap_keys.add(key)

    eff_lot_id = _all_s11_lids[0] if _all_s11_lids else None
    eff_sample_id = _all_s11_sids[0] if _all_s11_sids else None

    def _eval_s11_quantitative_view(view: dict[str, Any]) -> tuple[list[str], list[str]]:
        v_fails: list[str] = []
        v_holds: list[str] = []

        def _add_f(code: str) -> None:
            if code not in v_fails:
                v_fails.append(code)

        def _add_h(code: str) -> None:
            if code not in v_holds:
                v_holds.append(code)

        if isinstance(l2_binding, dict) and outer_only.get('construct_sha256') is None and rec_in.get('construct_sha256') is None:
            _add_h('HOLD_S11_MISSING_OUTER_CONSTRUCT_SHA256')

        if (
            not isinstance(eff_lot_id, str)
            or not eff_lot_id.strip()
            or eff_lot_id != eff_lot_id.strip()
            or not isinstance(eff_sample_id, str)
            or not eff_sample_id.strip()
            or eff_sample_id != eff_sample_id.strip()
        ):
            _add_h('HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING')

        geno_class_pre = view.get('ciita_locus2_genotype_copy_class')
        _is_kiki_pre = geno_class_pre in ('KI/KI', '2x', '2X')
        _is_kiko_pre = geno_class_pre in ('KI/KO', '1x', '1X')

        phasing_any_bad = any(
            _pflag in view and view[_pflag] is not True
            for _pflag in (
                'phasing_verified',
                'ciita_ki_allele_informative_transcript_recovery_verified',
                'ciita_ko_allele_informative_transcript_recovery_verified',
                'ciita_ki_allele1_informative_transcript_recovery_verified',
                'ciita_ki_allele2_informative_transcript_recovery_verified',
            )
        )
        has_kiko_recov_keys = any(
            k in view
            for k in (
                'ciita_ki_allele_informative_transcript_recovery_verified',
                'ciita_ko_allele_informative_transcript_recovery_verified',
            )
        )
        has_kiki_recov_keys = any(
            k in view
            for k in (
                'ciita_ki_allele1_informative_transcript_recovery_verified',
                'ciita_ki_allele2_informative_transcript_recovery_verified',
            )
        )
        if _is_kiki_pre:
            if has_kiko_recov_keys:
                _add_h('HOLD_S11_GENOTYPE_CLASS_ALLELE_FIELD_MISMATCH')
            allele_recovery_attested = (
                view.get('ciita_ki_allele1_informative_transcript_recovery_verified') is True
                and view.get('ciita_ki_allele2_informative_transcript_recovery_verified') is True
            )
        elif _is_kiko_pre:
            if has_kiki_recov_keys:
                _add_h('HOLD_S11_GENOTYPE_CLASS_ALLELE_FIELD_MISMATCH')
            allele_recovery_attested = (
                view.get('ciita_ki_allele_informative_transcript_recovery_verified') is True
                and view.get('ciita_ko_allele_informative_transcript_recovery_verified') is True
            )
        else:
            allele_recovery_attested = False
        phasing_attested = (view.get('phasing_verified') is True) and allele_recovery_attested
        if phasing_any_bad or not phasing_attested:
            _add_h('HOLD_S11_CIITA_BIALLELIC_PHASING_UNVERIFIED')

        cis_any_bad = any(
            _cflag in view and view[_cflag] is not True
            for _cflag in ('cis_verified', 'ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified')
        )
        cis_attested = (
            view.get('cis_verified') is True
            or view.get('ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified') is True
        )
        if cis_any_bad or not cis_attested:
            _add_h('HOLD_S11_16P13_CIS_SENTINELS_UNVERIFIED')

        wt_any_bad = any(
            _wflag in view and view[_wflag] is not True
            for _wflag in ('wt_inducibility_verified', 'lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified')
        )
        wt_attested = (
            view.get('wt_inducibility_verified') is True
            or view.get('lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified') is True
        )
        if wt_any_bad or not wt_attested:
            _add_h('HOLD_S11_WT_HLAII_INDUCIBILITY_UNVERIFIED')

        if (
            ('rnai_attribution_verified' in view and view['rnai_attribution_verified'] is not True)
            or (has_mir_active and view.get('rnai_attribution_verified') is not True)
        ):
            _add_h('HOLD_S11_RNAI_ATTRIBUTION_UNVERIFIED')

        # 3. Strict replicate count validation across ALL 3 supported aliases (including `s11_replicate_n`):
        _rep_alias_keys = ('obs_stage0_independent_lot_replicates_n', 'replicates_n', 's11_replicate_n')
        replicate_aliases = [view[k] for k in _rep_alias_keys if k in view]
        t95: float | None = None
        rep_n_raw = replicate_aliases[0] if replicate_aliases else None
        if not replicate_aliases or any(
            not _is_strict_int_iip0024(v) or not 3 <= v <= STAGE0_MAX_INDEPENDENT_LOT_REPLICATES or v != rep_n_raw
            for v in replicate_aliases
        ):
            _add_h('HOLD_S11_INVALID_OR_MISSING_REPLICATES_N')
        else:
            t95 = student_t_ppf_one_sided_95(float(rep_n_raw - 1))
            if rep_n_raw < 6:
                _add_h('HOLD_S11_REPLICATES_BELOW_6_LOT_MINIMUM')

        def _interval_t(
            mean_k: str,
            se_k: str,
            *,
            lo_dom: float = 0.0,
            hi_dom: float = 100.0,
        ) -> tuple[bool, float | None, float | None]:
            mv, sv = view.get(mean_k), view.get(se_k)
            if mv is None and sv is None:
                return False, None, None
            if (
                not _is_strict_finite_number_iip0024(mv)
                or not _is_strict_finite_number_iip0024(sv)
                or float(sv) <= 0.0
                or float(sv) > (hi_dom - lo_dom)
                or not (lo_dom <= float(mv) <= hi_dom)
            ):
                _add_h('HOLD_S11_INVALID_SE_OR_OUT_OF_RANGE_PHYSIOLOGICAL_READOUT')
                return False, None, None
            if t95 is None:
                return False, None, None
            mf, sf = float(mv), float(sv)
            return True, mf - t95 * sf, mf + t95 * sf

        def _pair_num_ok(mean_k: str, se_k: str, lo_dom: float, hi_dom: float) -> bool:
            mv, sv = view.get(mean_k), view.get(se_k)
            return bool(
                _is_strict_finite_number_iip0024(mv)
                and _is_strict_finite_number_iip0024(sv)
                and float(sv) > 0.0
                and float(sv) <= (hi_dom - lo_dom)
                and (lo_dom <= float(mv) <= hi_dom)
            )

        # 4. Post-cuff complement, conjunctive GY-023 ADCC, PLVAP+ iEC, efferocytosis & Gate 0G pericytes:
        ok_ap, ap_l95, ap_u95 = _interval_t('obs_l2_v3_ap_c5b9_pct', 'obs_l2_v3_ap_c5b9_se')
        ok_cls, cls_l95, cls_u95 = _interval_t('obs_l2_v3_classical_c4b2a_lysis_pct', 'obs_l2_v3_classical_c4b2a_lysis_se')
        ok_lec, lec_l95, lec_u95 = _interval_t('obs_l2_v3_lectin_c4b2a_lysis_pct', 'obs_l2_v3_lectin_c4b2a_lysis_se')
        ok_nhs, _nhs_l95, nhs_u95 = _interval_t('obs_l2_v3_overall_20pct_nhs_lysis_pct', 'obs_l2_v3_overall_20pct_nhs_lysis_se')
        ok_atot, _at_l95, adcc_tot_u95 = _interval_t('obs_l2_v3_intact_serum_adcc_lysis_pct', 'obs_l2_v3_intact_serum_adcc_lysis_se')
        ok_adep, _ad_l95, _ad_u95 = _interval_t('obs_l2_v3_igg_depleted_serum_adcc_lysis_pct', 'obs_l2_v3_igg_depleted_serum_adcc_lysis_se')
        ok_ainc, _ai_l95, adcc_inc_u95 = _interval_t('obs_l2_v3_igg_increment_adcc_lysis_pct', 'obs_l2_v3_igg_increment_adcc_lysis_se', lo_dom=-100.0, hi_dom=100.0)

        # Validate paired GY-023 ADCC mean equality (exact Decimal) and Cauchy-Schwarz SE bounds EVEN WHEN t95 is None (#20, #30):
        num_atot = _pair_num_ok('obs_l2_v3_intact_serum_adcc_lysis_pct', 'obs_l2_v3_intact_serum_adcc_lysis_se', 0.0, 100.0)
        num_adep = _pair_num_ok('obs_l2_v3_igg_depleted_serum_adcc_lysis_pct', 'obs_l2_v3_igg_depleted_serum_adcc_lysis_se', 0.0, 100.0)
        num_ainc = _pair_num_ok('obs_l2_v3_igg_increment_adcc_lysis_pct', 'obs_l2_v3_igg_increment_adcc_lysis_se', -100.0, 100.0)
        adcc_pair_invalid = False
        adcc_se_invalid = False
        if num_atot and num_ainc:
            _intact_d = Decimal(str(view['obs_l2_v3_intact_serum_adcc_lysis_pct']))
            _inc_d = Decimal(str(view['obs_l2_v3_igg_increment_adcc_lysis_pct']))
            _implied_dep_d = _intact_d - _inc_d
            if _implied_dep_d < Decimal('0') or _implied_dep_d > Decimal('100'):
                adcc_pair_invalid = True
                _add_h('HOLD_S11_INVALID_SE_OR_OUT_OF_RANGE_PHYSIOLOGICAL_READOUT')
            _dep_sup = view.get('obs_l2_v3_igg_depleted_serum_adcc_lysis_pct')
            _adcc_pair_id = view.get('obs_l2_v3_adcc_serum_donor_pair_id')
            _adcc_intact_sid = view.get('obs_l2_v3_adcc_intact_specimen_id')
            _adcc_dep_sid = view.get('obs_l2_v3_adcc_igg_depleted_specimen_id')
            if (
                not num_adep
                or _dep_sup is None
                or not _is_strict_finite_number_iip0024(_dep_sup)
                or Decimal(str(_dep_sup)) < Decimal('0')
                or Decimal(str(_dep_sup)) > Decimal('100')
                or Decimal(str(_dep_sup)) != _implied_dep_d
                or not isinstance(_adcc_pair_id, str)
                or not _adcc_pair_id.strip()
                or _adcc_pair_id != _adcc_pair_id.strip()
                or not isinstance(_adcc_intact_sid, str)
                or not _adcc_intact_sid.strip()
                or _adcc_intact_sid != _adcc_intact_sid.strip()
                or not isinstance(_adcc_dep_sid, str)
                or not _adcc_dep_sid.strip()
                or _adcc_dep_sid != _adcc_dep_sid.strip()
                or _adcc_intact_sid == _adcc_dep_sid
            ):
                adcc_pair_invalid = True
                _add_h('HOLD_S11_UNVERIFIED_PAIRED_GY023_ADCC_SPECIMENS')
        if num_atot and num_adep and num_ainc and not _paired_difference_se_consistent(
            view['obs_l2_v3_intact_serum_adcc_lysis_se'],
            view['obs_l2_v3_igg_depleted_serum_adcc_lysis_se'],
            view['obs_l2_v3_igg_increment_adcc_lysis_se'],
        ):
            adcc_se_invalid = True
            _add_h('HOLD_S11_IMPOSSIBLE_PAIRED_ADCC_UNCERTAINTY')
        adcc_pair_valid = bool(ok_atot and ok_adep and ok_ainc and not adcc_pair_invalid and not adcc_se_invalid)

        ok_eabs, eff_abs_l95, _ea_u95 = _interval_t('obs_l2_v3_exact_stack_effero_abs_pct', 'obs_l2_v3_exact_stack_effero_abs_se')
        ok_ewt, eff_wt_l95, _ew_u95 = _interval_t('obs_l2_v3_exact_stack_effero_pct_of_wt', 'obs_l2_v3_exact_stack_effero_pct_of_wt_se', lo_dom=0.0, hi_dom=200.0)
        ok_peri, peri_l95, _pe_u95 = _interval_t('obs_l2_v3_pericyte_coverage_pct', 'obs_l2_v3_pericyte_coverage_se')
        ok_t50, t50_l95, t50_u95 = _interval_t('obs_l2_v3_plvap_iec_network_t50_h', 'obs_l2_v3_plvap_iec_network_t50_se', lo_dom=0.1, hi_dom=240.0)
        ok_ivia, ivia_l95, ivia_u95 = _interval_t('obs_l2_v3_plvap_iec_viability_pct', 'obs_l2_v3_plvap_iec_viability_se')

        if ok_ap and ap_u95 is not None and ap_u95 > 5.0:
            _add_f('FAIL_S11_POST_CUFF_AP_C5B9_EXCEEDS_5PCT')
        if ok_cls and cls_l95 is not None and cls_u95 is not None:
            if cls_l95 > 10.0:
                _add_f('FAIL_S11_POST_CUFF_CLASSICAL_C4B2A_LYSIS_EXCEEDS_10PCT')
            elif cls_u95 > 10.0:
                _add_h('HOLD_S11_INCONCLUSIVE_CLASSICAL_C4B2A_STRADDLES_10PCT')
        if ok_lec and lec_l95 is not None and lec_u95 is not None:
            if lec_l95 > 10.0:
                _add_f('FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT')
            elif lec_u95 > 10.0:
                _add_h('HOLD_S11_INCONCLUSIVE_LECTIN_C4B2A_STRADDLES_10PCT')
        if ok_nhs and nhs_u95 is not None and nhs_u95 > 5.0:
            _add_f('FAIL_S11_PLAYBOOK_OVERALL_20PCT_NHS_LYSIS_EXCEEDS_5PCT')
        if (ok_atot and adcc_tot_u95 is not None and adcc_tot_u95 > 5.0) or (
            adcc_pair_valid and adcc_inc_u95 is not None and adcc_inc_u95 > 5.0
        ):
            _add_f('FAIL_S11_CONJUNCTIVE_GY023_ADCC_EXCEEDS_5PCT')
        if (ok_eabs and eff_abs_l95 is not None and eff_abs_l95 < 80.0) or (ok_ewt and eff_wt_l95 is not None and eff_wt_l95 < 70.0):
            _add_f('FAIL_S11_EXACT_STACK_EFFEROCYTOSIS_BELOW_80ABS_OR_70WT')
        if ok_peri and peri_l95 is not None and peri_l95 < 70.0:
            _add_f('FAIL_S11_GATE0G_PERICYTE_COVERAGE_BELOW_70PCT')
        if ok_t50 and t50_l95 is not None and t50_u95 is not None:
            if t50_l95 > 36.0:
                _add_f('FAIL_S11_PLVAP_IEC_NETWORK_MORPHOGENESIS_EXCEEDS_36H')
            elif t50_u95 > 36.0:
                _add_h('HOLD_S11_INCONCLUSIVE_PLVAP_IEC_NETWORK_T50_STRADDLES_36H')
        if ok_ivia and ivia_l95 is not None and ivia_u95 is not None:
            if ivia_u95 < 85.0:
                _add_f('FAIL_S11_PLVAP_IEC_VIABILITY_BELOW_85PCT')
            elif ivia_l95 < 85.0:
                _add_h('HOLD_S11_INCONCLUSIVE_PLVAP_IEC_VIABILITY_STRADDLES_85PCT')

        if not all((ok_ap, ok_cls, ok_lec, ok_nhs, ok_atot, ok_ainc, ok_eabs, ok_ewt, ok_peri, ok_t50, ok_ivia)):
            _add_h('HOLD_S11_INCOMPLETE_COMPLEMENT_ADCC_VASCULAR_EFFEROCYTOSIS_OR_PERICYTE_PANEL')

        # 5. Reconciled AAV6 species fractions & quantitative chr16 HA WGS + tandem poly(A) readthrough + copy-specific expression:
        aav_m = view.get('obs_l2_v3_aav6_full_length_monomer_pct')
        aav_t = view.get('obs_l2_v3_aav6_truncated_species_pct')
        aav_d = view.get('obs_l2_v3_aav6_dimeric_or_intermediate_pct')
        if any(x is None for x in (aav_m, aav_t, aav_d)):
            if any(x is not None for x in (aav_m, aav_t, aav_d)):
                _add_h('HOLD_S11_AAV6_SPECIES_DENOMINATOR_UNRECONCILED')
            else:
                _add_h('HOLD_S11_MISSING_RECONCILED_AAV6_SPECIES_FRACTIONS')
        elif not all(_is_strict_finite_number_iip0024(x) for x in (aav_m, aav_t, aav_d)):
            _add_h('HOLD_S11_AAV6_SPECIES_DENOMINATOR_UNRECONCILED')
        else:
            fm, ft, fd = float(aav_m), float(aav_t), float(aav_d)
            if fm < 0.0 or ft < 0.0 or fd < 0.0 or fm > 100.0 or ft > 100.0 or fd > 100.0 or abs((fm + ft + fd) - 100.0) > 1e-6:
                _add_h('HOLD_S11_AAV6_SPECIES_DENOMINATOR_UNRECONCILED')
            elif fm < 80.0 or ft >= 10.0 or fd > 15.0:
                _add_f('FAIL_S11_AAV6_HOMOGENEITY_SPEC_BREACHED')

        wgs_pass = view.get('obs_l2_v3_chr16_ha_wgs_pass')
        wgs_vcount = view.get('obs_l2_v3_chr16_ha_wgs_variant_count')
        if wgs_pass is False or (_is_strict_int_iip0024(wgs_vcount) and wgs_vcount > 0):
            _add_f('FAIL_S11_CHR16_HA_WGS_VARIANT_DETECTED')
        elif wgs_pass is not True or not _is_strict_int_iip0024(wgs_vcount) or wgs_vcount != 0:
            _add_h('HOLD_S11_MISSING_CHR16_HA_WGS')

        ok_rt, rt_l95, rt_u95 = _interval_t('obs_l2_v3_tandem_polya_readthrough_fold', 'obs_l2_v3_tandem_polya_readthrough_se', lo_dom=0.0, hi_dom=100.0)
        if not ok_rt:
            _add_h('HOLD_S11_MISSING_TANDEM_POLYA_READTHROUGH_EVIDENCE')
        elif rt_l95 is not None and rt_l95 > 2.0:
            _add_f('FAIL_S11_TANDEM_POLYA_READTHROUGH_EXCEEDS_2X')
        elif rt_u95 is not None and rt_u95 > 2.0:
            _add_h('HOLD_S11_INCONCLUSIVE_TANDEM_POLYA_READTHROUGH_STRADDLES_2X')

        geno_class = view.get('ciita_locus2_genotype_copy_class')
        if geno_class not in ('KI/KO', '1x', '1X', 'KI/KI', '2x', '2X'):
            _add_h('HOLD_S11_MISSING_OR_INVALID_CIITA_GENOTYPE_COPY_CLASS')
            is_2x_kiki = False
        else:
            is_2x_kiki = geno_class in ('KI/KI', '2x', '2X')

        cfh_qb_min = 15000.0 if is_2x_kiki else 10000.0
        ok_cfh_qb, cfh_qb_l95, cfh_qb_u95 = _interval_t('obs_l2_v3_cfh_stp_gpi_quantibrite_molecules', 'obs_l2_v3_cfh_stp_gpi_quantibrite_se', lo_dom=0.0, hi_dom=2000000.0)
        ok_u2_qb, u2_qb_l95, _u2_qb_u95 = _interval_t('obs_l2_v3_unit2_copy_expression_quantibrite_or_fold', 'obs_l2_v3_unit2_copy_expression_se', lo_dom=0.0, hi_dom=2000000.0)
        if not ok_cfh_qb or not ok_u2_qb:
            _add_h('HOLD_S11_MISSING_COPY_SPECIFIC_QUANTIBRITE_EXPRESSION')
        elif (cfh_qb_l95 is not None and cfh_qb_l95 < cfh_qb_min) or (cfh_qb_u95 is not None and cfh_qb_u95 > 500000.0) or (u2_qb_l95 is not None and u2_qb_l95 < 1.0):
            _add_f('FAIL_S11_COPY_SPECIFIC_QUANTIBRITE_EXPRESSION_OUT_OF_BOUNDS')

        # 6. Quantitative phased biallelic CIITA Ex4->Ex6: independent per-allele exact CP validation,
        #    full KI/KI vs KI/KO field binding, `bit_length() <= 53` overflow guard, and modern/legacy alias conflict preservation:
        has_kiko_fields = any(k in view for k in _S11_KIKO_ONLY_FIELDS)
        has_kiki_fields = any(k in view for k in _S11_KIKI_ONLY_FIELDS)

        def _collect_alias_vals(keys: tuple[str, ...]) -> tuple[list[Any], bool]:
            vals = [view[k] for k in keys if k in view]
            if not vals:
                return [], False
            sigs: list[str] = []
            distinct: list[Any] = []
            bad_alias = False
            for v in vals:
                if _is_strict_int_iip0024(v) and abs(v).bit_length() > 53:
                    bad_alias = True
                    continue
                try:
                    sv = json.dumps(v, sort_keys=True, allow_nan=False)
                except (TypeError, ValueError, OverflowError):
                    bad_alias = True
                    continue
                if sv not in sigs:
                    sigs.append(sv)
                    distinct.append(v)
            return distinct, bool(bad_alias or len(distinct) > 1)

        allele_specs: list[tuple[list[Any], list[Any], list[str]]] = []
        if is_2x_kiki:
            if has_kiko_fields or not has_kiki_fields:
                _add_h('HOLD_S11_GENOTYPE_CLASS_ALLELE_FIELD_MISMATCH')
            for idx_a in (1, 2):
                n_vals, n_conf = _collect_alias_vals((f'obs_ciita_ki_allele{idx_a}_informative_reads_n', f'obs_ciita_allele{idx_a}_informative_reads_n'))
                k_vals, k_conf = _collect_alias_vals((f'obs_ciita_ki_allele{idx_a}_ex4_ex6_aberrant_skip_reads_k', f'obs_ciita_allele{idx_a}_ex4_ex6_aberrant_skip_reads_k'))
                if n_conf or k_conf:
                    _add_h('HOLD_S11_CONFLICTING_CIITA_ALLELE_ALIASES')
                pct_ks = [pk for pk in (f'obs_ciita_ki_allele{idx_a}_ex4_ex6_aberrant_skip_pct', f'obs_ciita_allele{idx_a}_ex4_ex6_aberrant_skip_pct') if pk in view]
                allele_specs.append((k_vals, n_vals, pct_ks))
            ok_a1_m, a1_l_m, a1_u_m = _interval_t('obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ki_allele1_ex4_ex6_aberrant_skip_se')
            ok_a1_l, a1_l_l, a1_u_l = _interval_t('obs_ciita_allele1_ex4_ex6_aberrant_skip_pct', 'obs_ciita_allele1_ex4_ex6_aberrant_skip_se')
            ok_a2_m, a2_l_m, a2_u_m = _interval_t('obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ki_allele2_ex4_ex6_aberrant_skip_se')
            ok_a2_l, a2_l_l, a2_u_l = _interval_t('obs_ciita_allele2_ex4_ex6_aberrant_skip_pct', 'obs_ciita_allele2_ex4_ex6_aberrant_skip_se')
            t_sk_intervals = [(ok_a1_m, a1_l_m, a1_u_m), (ok_a1_l, a1_l_l, a1_u_l), (ok_a2_m, a2_l_m, a2_u_m), (ok_a2_l, a2_l_l, a2_u_l)]
            ok_a1_sk = ok_a1_m or ok_a1_l
            ok_a2_sk = ok_a2_m or ok_a2_l
        else:
            if has_kiki_fields:
                _add_h('HOLD_S11_GENOTYPE_CLASS_ALLELE_FIELD_MISMATCH')
            a1_n_v = [view['obs_ciita_ki_allele_informative_reads_n']] if 'obs_ciita_ki_allele_informative_reads_n' in view else []
            a1_k_v = [view['obs_ciita_ki_ex4_ex6_aberrant_skip_reads_k']] if 'obs_ciita_ki_ex4_ex6_aberrant_skip_reads_k' in view else []
            a2_n_v = [view['obs_ciita_ko_allele_informative_reads_n']] if 'obs_ciita_ko_allele_informative_reads_n' in view else []
            a2_k_v = [view['obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k']] if 'obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k' in view else []
            allele_specs = [
                (a1_k_v, a1_n_v, ['obs_ciita_ki_ex4_ex6_aberrant_skip_pct'] if 'obs_ciita_ki_ex4_ex6_aberrant_skip_pct' in view else []),
                (a2_k_v, a2_n_v, ['obs_ciita_ko_ex4_ex6_aberrant_skip_pct'] if 'obs_ciita_ko_ex4_ex6_aberrant_skip_pct' in view else []),
            ]
            ok_a1_sk, a1_sk_l95, a1_sk_u95 = _interval_t('obs_ciita_ki_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ki_ex4_ex6_aberrant_skip_se')
            ok_a2_sk, a2_sk_l95, a2_sk_u95 = _interval_t('obs_ciita_ko_ex4_ex6_aberrant_skip_pct', 'obs_ciita_ko_ex4_ex6_aberrant_skip_se')
            t_sk_intervals = [(ok_a1_sk, a1_sk_l95, a1_sk_u95), (ok_a2_sk, a2_sk_l95, a2_sk_u95)]

        ok_ed_hla, ed_hla_l95, ed_hla_u95 = _interval_t('obs_edited_72h_ifng_hla_dr_dp_dq_pos_pct', 'obs_edited_72h_ifng_hla_dr_dp_dq_pos_se')
        ok_wt_hla, wt_hla_l95, _wt_hla_u95 = _interval_t('obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_pct', 'obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_se')

        if any(ok_i and l_i is not None and l_i > 1.0 for ok_i, l_i, _u_i in t_sk_intervals):
            _add_f('FAIL_S11_CIITA_ABERRANT_EX4_EX6_SKIP_EXCEEDS_1PCT')
        elif any(ok_i and u_i is not None and u_i > 1.0 for ok_i, _l_i, u_i in t_sk_intervals):
            _add_h('HOLD_S11_INCONCLUSIVE_CIITA_EX4_EX6_SKIP_STRADDLES_1PCT')

        if ok_ed_hla and ed_hla_l95 is not None and ed_hla_l95 > 1.0:
            _add_f('FAIL_S11_EDITED_CELL_HLA_II_ESCAPE_DETECTED')
        elif ok_ed_hla and ed_hla_u95 is not None and ed_hla_u95 > 1.0:
            _add_h('HOLD_S11_INCONCLUSIVE_EDITED_CELL_HLA_II_STRADDLES_1PCT')

        # Evaluate EACH allele's exact binomial counts independently BEFORE checking full-panel completeness:
        all_alleles_counts_complete = True
        for k_cands, n_cands, pct_keys in allele_specs:
            valid_pairs: list[tuple[int, int]] = []
            allele_any_invalid = False
            for kv in k_cands:
                for nv in n_cands:
                    if (
                        _is_strict_int_iip0024(nv)
                        and nv >= 300
                        and nv.bit_length() <= 53
                        and _is_strict_int_iip0024(kv)
                        and 0 <= kv <= nv
                    ):
                        valid_pairs.append((kv, nv))
                    else:
                        allele_any_invalid = True
            if not valid_pairs or allele_any_invalid or len(k_cands) != 1 or len(n_cands) != 1:
                all_alleles_counts_complete = False
            for _sk_k, _rn_k in valid_pairs:
                _emp_pct_d = (Decimal(_sk_k) * Decimal('100')) / Decimal(_rn_k)
                for _pk_name in pct_keys:
                    _raw_p = view.get(_pk_name)
                    if _is_strict_finite_number_iip0024(_raw_p):
                        if abs(Decimal(str(_raw_p)) - _emp_pct_d) > Decimal('0.05'):
                            _add_h('HOLD_S11_CIITA_SKIP_COUNT_PERCENTAGE_MISMATCH')
                _cp_l95_pct = 100.0 * _exact_one_sided_cp_lower_iip0031(_sk_k, _rn_k, 0.05)
                _cp_u95_pct = 100.0 * clopper_pearson_upper95_independent(_sk_k, _rn_k, alpha=0.05)
                if _cp_l95_pct > 1.0:
                    _add_f('FAIL_S11_CIITA_ABERRANT_EX4_EX6_SKIP_EXCEEDS_1PCT')
                elif _cp_u95_pct > 1.0:
                    _add_h('HOLD_S11_INCONCLUSIVE_CIITA_EX4_EX6_SKIP_STRADDLES_1PCT')

        if not all_alleles_counts_complete or not all((ok_a1_sk, ok_a2_sk, ok_ed_hla, ok_wt_hla)):
            _add_h('HOLD_S11_MISSING_QUANTITATIVE_CIITA_BIALLELIC_SPLICE_OR_HLAII_EVIDENCE')
        elif wt_hla_l95 is None or wt_hla_l95 < 50.0:
            _add_h('HOLD_S11_UNINDUCIBLE_WT_HLAII_POSITIVE_CONTROL')

        for cis_gene in ('socs1', 'dexi', 'clec16a'):
            lo_g = view.get(f'obs_cis_{cis_gene}_tost90_ratio_lo')
            hi_g = view.get(f'obs_cis_{cis_gene}_tost90_ratio_hi')
            if (
                not _is_strict_finite_number_iip0024(lo_g)
                or not _is_strict_finite_number_iip0024(hi_g)
                or float(lo_g) <= 0.0
                or float(lo_g) >= float(hi_g)
            ):
                _add_h('HOLD_S11_MISSING_OR_INVALID_16P13_CIS_SENTINEL_TOST')
            else:
                flo, fhi = float(lo_g), float(hi_g)
                if flo > 1.25 or fhi < 0.80:
                    _add_f('FAIL_S11_16P13_CIS_SENTINEL_TOST_BREACHED')
                elif flo < 0.80 or fhi > 1.25:
                    _add_h('HOLD_S11_INCONCLUSIVE_16P13_CIS_SENTINEL_TOST')

        cn16 = view.get('obs_16p13_ddpcr_cn_ratio')
        baf16 = view.get('obs_16p13_phased_baf')
        if not _is_strict_finite_number_iip0024(cn16) or not _is_strict_finite_number_iip0024(baf16):
            _add_h('HOLD_S11_MISSING_16P13_CN_OR_PHASED_BAF')
        elif not (
            Decimal('1.85') <= Decimal(str(cn16)) <= Decimal('2.15')
            and Decimal('0.45') <= Decimal(str(baf16)) <= Decimal('0.55')
        ):
            _add_f('FAIL_S11_16P13_LOH_OR_CN_BREACHED')

        # 7. Quantitative Locus-2 DEC-107 5'-RACE across ALL 4/4 promoter-time datasets (reconciling legacy aliases & 0 <= k <= n):
        eef0_n_vals, eef0_n_conf = _collect_alias_vals(('obs_l2_v3_unit1_eef1a1_resting_0h_5race_umi_n_eff', 'obs_l2_v3_unit1_eef1a1_5race_umi_n_eff'))
        eef0_k_vals, eef0_k_conf = _collect_alias_vals(('obs_l2_v3_unit1_eef1a1_resting_0h_unexpected_tss_umi_k', 'obs_l2_v3_unit1_eef1a1_unexpected_tss_umi_k'))
        if eef0_n_conf or eef0_k_conf:
            _add_h('HOLD_S11_CONFLICTING_RACE_ALIASES')
        race_4_specs = (
            (eef0_n_vals, eef0_k_vals),
            (
                [view['obs_l2_v3_unit1_eef1a1_stim_24h_5race_umi_n_eff']] if 'obs_l2_v3_unit1_eef1a1_stim_24h_5race_umi_n_eff' in view else [],
                [view['obs_l2_v3_unit1_eef1a1_stim_24h_unexpected_tss_umi_k']] if 'obs_l2_v3_unit1_eef1a1_stim_24h_unexpected_tss_umi_k' in view else [],
            ),
            (
                [view['obs_l2_v3_unit2_resting_0h_5race_umi_n_eff']] if 'obs_l2_v3_unit2_resting_0h_5race_umi_n_eff' in view else [],
                [view['obs_l2_v3_unit2_resting_0h_unexpected_tss_umi_k']] if 'obs_l2_v3_unit2_resting_0h_unexpected_tss_umi_k' in view else [],
            ),
            (
                [view['obs_l2_v3_unit2_stim_24h_5race_umi_n_eff']] if 'obs_l2_v3_unit2_stim_24h_5race_umi_n_eff' in view else [],
                [view['obs_l2_v3_unit2_stim_24h_unexpected_tss_umi_k']] if 'obs_l2_v3_unit2_stim_24h_unexpected_tss_umi_k' in view else [],
            ),
        )
        for rn_cands, rk_cands in race_4_specs:
            if not rn_cands or not rk_cands:
                _add_h('HOLD_S11_MISSING_LOCUS2_DUAL_PROMOTER_5PRIME_RACE_3689')
                continue
            ds_valid_any = False
            ds_invalid_any = False
            for rn_v in rn_cands:
                for rk_v in rk_cands:
                    if (
                        not _is_strict_int_iip0024(rn_v)
                        or rn_v < 3689
                        or rn_v.bit_length() > 53
                        or not _is_strict_int_iip0024(rk_v)
                        or rk_v < 0
                        or rk_v > rn_v
                    ):
                        ds_invalid_any = True
                        if _is_strict_int_iip0024(rn_v) and _is_strict_int_iip0024(rk_v) and rk_v > rn_v:
                            _add_h('HOLD_S11_INVALID_RACE_COUNT_EXCEEDS_DENOMINATOR')
                    else:
                        ds_valid_any = True
                        if rk_v > 0:
                            _add_f('FAIL_S11_LOCUS2_UNEXPECTED_5PRIME_RACE_TSS')
            if not ds_valid_any or ds_invalid_any:
                _add_h('HOLD_S11_MISSING_LOCUS2_DUAL_PROMOTER_5PRIME_RACE_3689')

        is_syn5_donor = (
            l2_binding.get('unit2_promoter') == 'Syn_5xISRE_minP'
            if isinstance(l2_binding, dict)
            else (donor_id != 'L2_V3_CD274P_ZERO_TATATA')
        )
        if is_syn5_donor:
            u2_tss_pos = view.get('obs_l2_v3_unit2_tss_pos1_in_donor')
            if not _is_strict_int_iip0024(u2_tss_pos):
                _add_h('HOLD_S11_PENDING_SYN5XISRE_5PRIME_RACE_TSS_COORDINATE')
            elif u2_tss_pos <= 626 or u2_tss_pos > 698:
                _add_f('FAIL_S11_SYN5XISRE_UPSTREAM_ATG_INCLUDED_OR_INVALID_TSS')

        # 8. Quantitative LONGITUDINAL Stage-0 AND Stage-7 bisulfite CpG methylation & CTCF ChIP:
        ok_me0_eef, me0_eef_l95, me0_eef_u95 = _interval_t('obs_l2_v3_stage0_eef1a1_cpg_methylation_pct', 'obs_l2_v3_stage0_eef1a1_cpg_methylation_se')
        ok_me0_chs4, me0_chs4_l95, me0_chs4_u95 = _interval_t('obs_l2_v3_stage0_chs4_cpg_methylation_pct', 'obs_l2_v3_stage0_chs4_cpg_methylation_se')
        ok_me0_dm1, me0_dm1_l95, me0_dm1_u95 = _interval_t('obs_l2_v3_stage0_dm1_cpg_methylation_pct', 'obs_l2_v3_stage0_dm1_cpg_methylation_se')
        ok_me7_eef, me7_eef_l95, me7_eef_u95 = _interval_t('obs_l2_v3_stage7_eef1a1_cpg_methylation_pct', 'obs_l2_v3_stage7_eef1a1_cpg_methylation_se')
        ok_me7_chs4, me7_chs4_l95, me7_chs4_u95 = _interval_t('obs_l2_v3_stage7_chs4_cpg_methylation_pct', 'obs_l2_v3_stage7_chs4_cpg_methylation_se')
        ok_me7_dm1, me7_dm1_l95, me7_dm1_u95 = _interval_t('obs_l2_v3_stage7_dm1_cpg_methylation_pct', 'obs_l2_v3_stage7_dm1_cpg_methylation_se')
        ok_ch0_chs4, ch0_chs4_l95, _ch0_chs4_u95 = _interval_t('obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment_se', lo_dom=0.0, hi_dom=1000.0)
        ok_ch0_dm1, ch0_dm1_l95, _ch0_dm1_u95 = _interval_t('obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment_se', lo_dom=0.0, hi_dom=1000.0)
        ok_ch7_chs4, ch7_chs4_l95, _ch7_chs4_u95 = _interval_t('obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment_se', lo_dom=0.0, hi_dom=1000.0)
        ok_ch7_dm1, ch7_dm1_l95, _ch7_dm1_u95 = _interval_t('obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment', 'obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment_se', lo_dom=0.0, hi_dom=1000.0)

        meth_tuples = (
            (ok_me0_eef, me0_eef_l95, me0_eef_u95),
            (ok_me0_chs4, me0_chs4_l95, me0_chs4_u95),
            (ok_me0_dm1, me0_dm1_l95, me0_dm1_u95),
            (ok_me7_eef, me7_eef_l95, me7_eef_u95),
            (ok_me7_chs4, me7_chs4_l95, me7_chs4_u95),
            (ok_me7_dm1, me7_dm1_l95, me7_dm1_u95),
        )
        if any(ok_m and ml is not None and ml > 20.0 for (ok_m, ml, _mu) in meth_tuples):
            _add_f('FAIL_S11_LOCUS2_CPG_HYPERMETHYLATION_EXCEEDS_20PCT')
        elif any(ok_m and mu is not None and mu > 20.0 for (ok_m, _ml, mu) in meth_tuples):
            _add_h('HOLD_S11_INCONCLUSIVE_LOCUS2_CPG_METHYLATION_STRADDLES_20PCT')

        if not all((ok_me0_eef, ok_me0_chs4, ok_me0_dm1, ok_me7_eef, ok_me7_chs4, ok_me7_dm1, ok_ch0_chs4, ok_ch0_dm1, ok_ch7_chs4, ok_ch7_dm1)):
            _add_h('HOLD_S11_MISSING_LONGITUDINAL_STAGE0_AND_STAGE7_METHYLATION_OR_CTCF_PANEL')
        elif any(cl is not None and cl < 5.0 for cl in (ch0_chs4_l95, ch0_dm1_l95, ch7_chs4_l95, ch7_dm1_l95)):
            _add_h('HOLD_S11_INSUFFICIENT_CTCF_CHIP_OCCUPANCY_BELOW_5X')

        if donor_id == 'L2_V3_BARRIER_CHS4_250' or pvid_canon == 'ISLET1_A4_L2_V3_BARRIER_CHS4_250':
            _bar_cmp_id = view.get('obs_l2_v3_barrier_matched_80bp_comparator_id')
            _bar_cmp_sha = view.get('obs_l2_v3_barrier_matched_80bp_comparator_sha256')
            _bar_cmp_sid = view.get('obs_l2_v3_barrier_matched_80bp_comparator_sample_id')
            ok_bar_me7, _bm7_l, _bm7_u = _interval_t('obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_pct', 'obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_se')
            ok_bar_ch7, _bc7_l, _bc7_u = _interval_t('obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment', 'obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment_se', lo_dom=0.0, hi_dom=1000.0)
            _exp_80bp_id = getattr(p_mod, 'LOCUS2_V3_BARRIER_MATCHED_80BP_COMPARATOR_ID', 'L2_V3_CO_PRIMARY_CD200')
            _exp_80bp_sha = '723e052172adfe5ae6c74e35e44cf727202d2716010772041d8c12abd76529c8'
            if (
                _bar_cmp_id != _exp_80bp_id
                or _bar_cmp_sha != _exp_80bp_sha
                or not isinstance(_bar_cmp_sid, str)
                or not _bar_cmp_sid.strip()
                or _bar_cmp_sid != _bar_cmp_sid.strip()
                or (isinstance(eff_sample_id, str) and _bar_cmp_sid == eff_sample_id)
                or (isinstance(view.get('rnai_attribution_scramble_sample_id'), str) and _bar_cmp_sid == view.get('rnai_attribution_scramble_sample_id'))
                or (isinstance(view.get('rnai_attribution_intronless_sample_id'), str) and _bar_cmp_sid == view.get('rnai_attribution_intronless_sample_id'))
                or not ok_bar_me7
                or not ok_bar_ch7
            ):
                _add_h('HOLD_S11_MISSING_LINKED_BARRIER_80BP_COMPARATOR_MEASUREMENT')

        # 9. Conditional MANF / CCL22 cargo-specific quantitative gates:
        unit2_cargo = l2_binding.get('unit2_cargo') if isinstance(l2_binding, dict) else None
        if unit2_cargo == 'MANF':
            ok_mnf_l3, mnf_l3_l95, _mnf_l3_u95 = _interval_t('obs_l2_v3_manf_locus3_nfkb_suppression_ratio', 'obs_l2_v3_manf_locus3_nfkb_suppression_se', lo_dom=0.0, hi_dom=10.0)
            ok_mnf_ki, mnf_ki_l95, mnf_ki_u95 = _interval_t('obs_l2_v3_manf_ki67_prolif_pct', 'obs_l2_v3_manf_ki67_prolif_se')
            ok_mnf_oc, mnf_oc_l95, _mnf_oc_u95 = _interval_t('obs_l2_v3_manf_ocr_ratio_vs_ctrl', 'obs_l2_v3_manf_ocr_ratio_vs_ctrl_se', lo_dom=0.0, hi_dom=10.0)
            if not all((ok_mnf_l3, ok_mnf_ki, ok_mnf_oc)):
                _add_h('HOLD_S11_MISSING_CONDITIONAL_MANF_GY001_PROLIF_OR_OCR_PANEL')
            elif (mnf_l3_l95 is not None and mnf_l3_l95 < 0.70) or (mnf_ki_u95 is not None and mnf_ki_u95 > 2.0) or (mnf_oc_l95 is not None and mnf_oc_l95 < 0.85):
                _add_f('FAIL_S11_CONDITIONAL_MANF_GY001_PROLIF_OR_OCR_BREACHED')
        elif unit2_cargo == 'CCL22':
            ok_ccl_r, ccl_r_l95, _ccl_r_u95 = _interval_t('obs_l2_v3_ccl22_treg_to_effector_chemotaxis_ratio', 'obs_l2_v3_ccl22_treg_to_effector_chemotaxis_se', lo_dom=0.0, hi_dom=100.0)
            ok_ccl_e, _ccl_e_l95, ccl_e_u95 = _interval_t('obs_l2_v3_ccl22_effector_chemotaxis_fold_ctrl', 'obs_l2_v3_ccl22_effector_chemotaxis_fold_se', lo_dom=0.0, hi_dom=100.0)
            if not all((ok_ccl_r, ok_ccl_e)):
                _add_h('HOLD_S11_MISSING_CONDITIONAL_CCL22_CHEMOTAXIS_PANEL')
            elif (ccl_r_l95 is not None and ccl_r_l95 < 2.0) or (ccl_e_u95 is not None and ccl_e_u95 > 1.25):
                _add_f('FAIL_S11_CONDITIONAL_CCL22_TREG_EFFECTOR_CHEMOTAXIS_BREACHED')

        # 10. Quantitative active-RNAi (miR) processing/KD/joint-stratum NK lysis & matched same-background B2M-/- attribution controls:
        if has_mir_active:
            ok_m375, m375_l95, m375_u95 = _interval_t('obs_l2_v3_small_rnaseq_mir375_retention_pct', 'obs_l2_v3_small_rnaseq_mir375_retention_se', lo_dom=0.0, hi_dom=200.0)
            ok_spl, spl_l95, spl_u95 = _interval_t('obs_l2_v3_intron_splicing_efficiency_pct', 'obs_l2_v3_intron_splicing_efficiency_se')
            ok_c276, c276_l95, c276_u95 = _interval_t('obs_l2_v3_cd276_kd_pct', 'obs_l2_v3_cd276_kd_se')
            ok_pvr, pvr_l95, pvr_u95 = _interval_t('obs_l2_v3_pvr_kd_pct', 'obs_l2_v3_pvr_kd_se')
            ok_nk, nk_l95, nk_u95 = _interval_t('obs_l2_v3_joint_stratum_nk_lysis_pct', 'obs_l2_v3_joint_stratum_nk_lysis_se')

            if ok_m375 and m375_u95 is not None and m375_u95 < 90.0:
                _add_f('FAIL_S11_ACTIVE_MIR_MIR375_DEPLETION_BELOW_90PCT')
            elif ok_m375 and m375_l95 is not None and m375_l95 < 90.0:
                _add_h('HOLD_S11_INCONCLUSIVE_MIR375_RETENTION_STRADDLES_90PCT')
            if ok_spl and spl_u95 is not None and spl_u95 < 85.0:
                _add_f('FAIL_S11_ACTIVE_MIR_INTRON_SPLICING_BELOW_85PCT')
            elif ok_spl and spl_l95 is not None and spl_l95 < 85.0:
                _add_h('HOLD_S11_INCONCLUSIVE_INTRON_SPLICING_STRADDLES_85PCT')
            if (ok_c276 and c276_u95 is not None and c276_u95 < 70.0) or (ok_pvr and pvr_u95 is not None and pvr_u95 < 70.0):
                _add_f('FAIL_S11_ACTIVE_MIR_CD276_OR_PVR_KD_BELOW_70PCT')
            elif (ok_c276 and c276_l95 is not None and c276_l95 < 70.0) or (ok_pvr and pvr_l95 is not None and pvr_l95 < 70.0):
                _add_h('HOLD_S11_INCONCLUSIVE_CD276_OR_PVR_KD_STRADDLES_70PCT')
            if ok_nk and nk_l95 is not None and nk_l95 > 8.0:
                _add_f('FAIL_S11_ACTIVE_MIR_JOINT_STRATUM_NK_LYSIS_EXCEEDS_8PCT')
            elif ok_nk and nk_u95 is not None and nk_u95 > 8.0:
                _add_h('HOLD_S11_INCONCLUSIVE_JOINT_STRATUM_NK_LYSIS_STRADDLES_8PCT')
            if not all((ok_m375, ok_spl, ok_c276, ok_pvr, ok_nk)):
                _add_h('HOLD_S11_MISSING_ACTIVE_MIR_PROCESSING_KD_OR_JOINT_STRATUM_NK_PANEL')

        scr_id = view.get('rnai_attribution_scramble_control_id')
        scr_bg = view.get('rnai_attribution_scramble_background_id')
        scr_sha = view.get('rnai_attribution_scramble_control_sha256')
        scr_sid = view.get('rnai_attribution_scramble_sample_id')
        intr_id = view.get('rnai_attribution_intronless_control_id')
        intr_bg = view.get('rnai_attribution_intronless_background_id')
        intr_sha = view.get('rnai_attribution_intronless_control_sha256')
        intr_sid = view.get('rnai_attribution_intronless_sample_id')
        rnai_claimed = view.get('rnai_specific_benefit_claimed') is True
        rnai_ctrl_flags_bad = (
            view.get('l2_v3_same_background_intronless_control_linked_to_cargo') is False
            or view.get('l2_v3_same_background_scramble_control_linked_to_cargo') is False
        )
        if rnai_claimed or has_mir_active or rnai_ctrl_flags_bad:
            _exp_scr_id = p_mod.LOCUS2_V3_SCRAMBLE_CONTROL_ID
            _pinned_intr_map = p_mod.LOCUS2_V3_PINNED_B2MKO_INTRONLESS_CONTROL_IDS
            _exp_intr_id = _pinned_intr_map.get(str(unit2_cargo)) if unit2_cargo is not None else None
            registered = {binding['locus2_donor_id']: binding for binding in p_mod.LOCUS2_V3_DONOR_BINDINGS.values()}
            scramble_reference = registered.get(_exp_scr_id, {})
            _exp_scr_sha = scramble_reference.get('locus2_donor_sha256')
            intronless_reference = registered.get(_exp_intr_id, {})
            _hex64_ok = lambda s: isinstance(s, str) and bool(re.fullmatch(r'[0-9a-f]{64}', s))
            _sid_ok = lambda s: isinstance(s, str) and bool(s.strip()) and s == s.strip()
            matched_ctrls_valid = (
                isinstance(scr_id, str)
                and scr_id == _exp_scr_id
                and scramble_reference.get('unit2_cargo') == unit2_cargo
                and _hex64_ok(_exp_scr_sha)
                and scr_sha == _exp_scr_sha
                and isinstance(scr_bg, str)
                and bool(bg_id)
                and scr_bg == bg_id
                and _sid_ok(scr_sid)
                and isinstance(intr_id, str)
                and _exp_intr_id is not None
                and intr_id == _exp_intr_id
                and _hex64_ok(intr_sha)
                and intr_sha != _exp_scr_sha
                and (expected_donor_sha is None or intr_sha != expected_donor_sha)
                and isinstance(intr_bg, str)
                and intr_bg == bg_id
                and not intr_id.startswith('L2_V3_DESIGN_B_')
                and not intr_bg.upper().startswith('DESIGN_B')
                and _sid_ok(intr_sid)
                and intr_sid != scr_sid
                and (not isinstance(eff_sample_id, str) or (scr_sid != eff_sample_id and intr_sid != eff_sample_id))
            )
            if rnai_claimed or has_mir_active:
                if (
                    view.get('l2_v3_same_background_intronless_control_linked_to_cargo') is not True
                    or view.get('l2_v3_same_background_scramble_control_linked_to_cargo') is not True
                    or not matched_ctrls_valid
                ):
                    _add_h('HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS')
            elif rnai_ctrl_flags_bad:
                _add_h('HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS')
            if (rnai_claimed or has_mir_active) and not (
                _hex64_ok(intr_sha)
                and len(set(intr_sha)) > 1
                and intr_sha == intronless_reference.get('locus2_donor_sha256')
                and intronless_reference.get('unit2_cargo') == unit2_cargo
                and intronless_reference.get('background_id') == bg_id
            ):
                _add_h('HOLD_S11_UNVERIFIED_INTRONLESS_CONTROL_DIGEST')

        _add_h('HOLD_S11_UNVERIFIED_ATTRIBUTION_CONTROL_MEASUREMENTS')
        return v_fails, v_holds

    if conflicting_overlap_keys:
        base_unconflicted = {
            k: v for k, v in merged.items() if k not in conflicting_overlap_keys
        }
        view_outer = dict(base_unconflicted)
        view_nested = dict(base_unconflicted)
        for ck in conflicting_overlap_keys:
            ov = payload.get(ck)
            nv = rec_in.get(ck)
            if ov is not None and ov != 'UNAVAILABLE':
                view_outer[ck] = ov
            if nv is not None and nv != 'UNAVAILABLE':
                view_nested[ck] = nv
        f_out, h_out = _eval_s11_quantitative_view(view_outer)
        f_nes, h_nes = _eval_s11_quantitative_view(view_nested)
        fails = list(dict.fromkeys(f_out + f_nes))
        holds = ['HOLD_S11_CONFLICTING_NESTED_EVIDENCE'] + [
            h for h in dict.fromkeys(h_out + h_nes) if h != 'HOLD_S11_CONFLICTING_NESTED_EVIDENCE'
        ]
    else:
        fails, holds = _eval_s11_quantitative_view(merged)

    if fails:
        if 'HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING' in holds:
            return _out(
                'HOLD_S11_UNLINKED_ADVERSE_SIGNAL_MISSING_LOT_OR_SAMPLE_BINDING',
                False,
                ['HOLD_S11_UNLINKED_ADVERSE_SIGNAL_MISSING_LOT_OR_SAMPLE_BINDING', *fails, *holds],
                unlinked_adverse_signals=list(fails),
                provenance_attributed=False,
            )
        return _out(fails[0], False, fails + holds, provenance_attributed=True)
    return _out(holds[0], False, holds, provenance_attributed=('HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING' not in holds))


def evaluate_stage0_product_contract(
    payload: dict[str, Any],
    *,
    root: Path | str | None = None,
) -> dict[str, Any]:
    """
    Direct Stage-0 product contract evaluator (Entry Point 1; unified with
    evaluate_stage0_research_and_molecular_guardrails, Atomic Stage-0 Ledger,
    Serialized JSON Ingestion, and Chronological Replay; IIP-0031; DEC-110).
    """
    return evaluate_stage0_research_and_molecular_guardrails(payload, root=root)


