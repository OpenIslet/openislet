"""One measurement contract for Python ingestion and generated MCP metadata."""
from __future__ import annotations

import ast
import importlib.util
import json
import math
from pathlib import Path
import re

CONTRACT_RELATIVE_PATH = Path('data/protocols/wetlab_endpoint_contract_v0.8.0.json')
STAGE0_DECISION_TREE_RELATIVE_PATH = Path('data/protocols/stage0_decision_tree_v1.json')
DECIMAL = re.compile(r'^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$')

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
    ('stage0_independent_lot_replicates_n', 'obs_stage0_independent_lot_replicates_n', 'obs_stage0_independent_lot_replicates_n_se', 1.0, 10000.0, 'independent_differentiation_lot_replicates_n_for_exact_student_t_bounds'),
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


def _poisson_exact_upper95_rate(k_events: float, n_trials: float) -> float:
    """Return exact one-sided 95% Poisson upper rate bound lambda_upper95(k) / n_trials (AUDIT-REF)."""
    if n_trials <= 0.0:
        return 1.0
    k_int = max(0, int(round(k_events)))
    if k_int < len(POISSON_UPPER95_LAMBDA_TABLE):
        lam_u95 = POISSON_UPPER95_LAMBDA_TABLE[k_int]
    else:
        # Wilson-Hilferty chi-square approximation for df = 2*(k+1) at z_0.95 = 1.6448536269514722:
        df = 2.0 * (k_int + 1)
        w = 1.0 - 2.0 / (9.0 * df) + 1.6448536269514722 * math.sqrt(2.0 / (9.0 * df))
        lam_u95 = 0.5 * df * (w ** 3)
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


def _load_product_contract_module(root: Path):
    spec = importlib.util.spec_from_file_location('openislet_product_contract', root / 'scripts/product_contract.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
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


def _betacf_iip0024(a: float, b: float, x: float) -> float:
    max_it = 240
    eps = 3.0e-15
    fpmin = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_it + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        del_v = d * c
        h *= del_v
        if abs(del_v - 1.0) < eps:
            break
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


def clopper_pearson_upper95_independent(k: int, n: int, alpha: float = 0.05) -> float:
    """Exact Clopper-Pearson one-sided (1 - alpha) upper bound for independent (k, n) counts (IIP-0024; DEC-103)."""
    if type(k) is not int or type(n) is not int or n <= 0 or k < 0 or k > n:
        raise ValueError('clopper_pearson_upper95_independent requires integer 0 <= k <= n, n > 0 (IIP-0024; DEC-103)')
    if k == n:
        return 1.0
    if k == 0:
        return -math.expm1(math.log(alpha) / float(n))
    lo, hi = 0.0, 1.0
    for _ in range(95):
        mid = 0.5 * (lo + hi)
        if _reg_inc_beta_iip0024(float(k + 1), float(n - k), mid) < (1.0 - alpha):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def clopper_pearson_lower95_independent(k: int, n: int, alpha: float = 0.05) -> float:
    if type(k) is not int or type(n) is not int or n <= 0 or k < 0 or k > n:
        raise ValueError('clopper_pearson_lower95_independent requires integer 0 <= k <= n, n > 0 (IIP-0024; DEC-103)')
    if k == 0:
        return 0.0
    return 1.0 - clopper_pearson_upper95_independent(n - k, n, alpha=alpha)


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
    # Canonical product resolution (IIP-0024; DEC-103): reject whitespace-padded IDs and always re-resolve from registry so caller-supplied product_def dicts cannot override execution_supported!
    obs_pvid_raw = observations.get('product_variant_id')
    if obs_pvid_raw is not None:
        if not isinstance(obs_pvid_raw, str) or not obs_pvid_raw.strip() or obs_pvid_raw != obs_pvid_raw.strip():
            raise ValueError(f'Invalid or whitespace-padded product_variant_id: {obs_pvid_raw!r} (IIP-0024; DEC-103)')
        canon_from_obs = p_mod.get_product_definition(obs_pvid_raw, root=root_path)
        if product_def is not None:
            if not isinstance(product_def, dict):
                raise ValueError('product_def must be a dictionary when provided (IIP-0024; DEC-103)')
            kw_pvid = product_def.get('product_variant_id')
            if not isinstance(kw_pvid, str) or p_mod.get_product_definition(kw_pvid, root=root_path)['product_variant_id'] != canon_from_obs['product_variant_id']:
                raise ValueError('Conflicting product_def and observations product_variant_id (IIP-0024; DEC-103)')
        product_def = canon_from_obs
    elif product_def is not None:
        if not isinstance(product_def, dict) or not isinstance(product_def.get('product_variant_id'), str):
            raise ValueError('product_def must contain a canonical string product_variant_id (IIP-0024; DEC-103)')
        product_def = p_mod.get_product_definition(product_def['product_variant_id'], root=root_path)

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
    }
    unknown_keys = sorted(set(observations) - set(spec_by_mean) - set(spec_by_se) - allowed_boolean_and_meta_flags)
    if unknown_keys:
        raise ValueError(f'Unrecognized Stage-0 decision-tree observation keys: {unknown_keys}')
    for _str_meta_k in ('donor_seedline_id', 'device_assembly_id', 'physical_route_id', 'assay_campaign_id', 'reset_prior_safety_stop_target_category', 'submitter_id'):
        if _str_meta_k in observations and observations[_str_meta_k] is not None:
            _sv = observations[_str_meta_k]
            if not isinstance(_sv, str) or not _sv.strip() or len(_sv) > 120:
                raise ValueError(f'{_str_meta_k} must be a nonempty string of at most 120 characters (IIP-0024; DEC-103)')
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

    icasp9_lower95 = (
        validated_stage0['obs_icasp9_log10_kill'] - 1.645 * validated_stage0_se.get('obs_icasp9_log10_kill_se', 0.0)
        if has_icasp9_obs
        else 4.5
    )
    # Exact one-sided 95% Poisson bound enforced for ALL k >= 0 (including k = 0 and regardless of whether an explicit Wald SE was supplied; AUDIT-REF; DEC-042):
    tp53_duplex_depth_explicitly_supplied = 'obs_tp53_ppm1d_bcor_duplex_depth' in validated_stage0
    vaf_duplex_depth = validated_stage0.get('obs_tp53_ppm1d_bcor_duplex_depth', 62958.0)
    tp53_duplex_depth_underpowered_below_62958 = bool(
        tp53_duplex_depth_explicitly_supplied and vaf_duplex_depth < 62958.0
    )
    vaf_wald_upper95 = (
        validated_stage0['obs_tp53_ppm1d_bcor_vaf'] + 1.645 * validated_stage0_se.get('obs_tp53_ppm1d_bcor_vaf_se', 0.0)
        if has_vaf_obs
        else 0.0
    )
    if has_vaf_obs:
        vaf_k_from_rate = validated_stage0['obs_tp53_ppm1d_bcor_vaf'] * vaf_duplex_depth
        if 'obs_tp53_ppm1d_bcor_mutant_molecules_k' in validated_stage0:
            k_supplied = validated_stage0['obs_tp53_ppm1d_bcor_mutant_molecules_k']
            if tp53_duplex_depth_explicitly_supplied and abs(k_supplied - vaf_k_from_rate) > 0.51:
                raise ValueError(
                    f"Conflicting obs_tp53_ppm1d_bcor_mutant_molecules_k ({k_supplied}) and "
                    f"obs_tp53_ppm1d_bcor_vaf * obs_tp53_ppm1d_bcor_duplex_depth ({vaf_k_from_rate:.4f}) (AUDIT-REF)"
                )
            vaf_k_est = max(k_supplied, vaf_k_from_rate)
        else:
            vaf_k_est = vaf_k_from_rate
        vaf_poisson_upper95 = _poisson_exact_upper95_rate(vaf_k_est, vaf_duplex_depth)
        vaf_upper95 = max(vaf_wald_upper95, vaf_poisson_upper95)
    else:
        vaf_upper95 = 0.0

    cast_duplex_depth_explicitly_supplied = 'obs_castseq_duplex_or_amplicon_depth' in validated_stage0
    cast_duplex_depth = validated_stage0.get('obs_castseq_duplex_or_amplicon_depth', 62958.0)
    cast_duplex_depth_underpowered_below_62958 = bool(
        cast_duplex_depth_explicitly_supplied and cast_duplex_depth < 62958.0
    )
    cast_wald_upper95 = (
        validated_stage0['obs_castseq_translocation_freq'] + 1.645 * validated_stage0_se.get('obs_castseq_translocation_se', 0.0)
        if has_cast_obs
        else 0.0
    )
    if has_cast_obs:
        cast_k_from_rate = validated_stage0['obs_castseq_translocation_freq'] * cast_duplex_depth
        if 'obs_castseq_mutant_junctions_k' in validated_stage0:
            cast_k_est = max(validated_stage0['obs_castseq_mutant_junctions_k'], cast_k_from_rate)
        else:
            cast_k_est = cast_k_from_rate
        cast_poisson_upper95 = _poisson_exact_upper95_rate(cast_k_est, cast_duplex_depth)
        cast_upper95 = max(cast_wald_upper95, cast_poisson_upper95)
    else:
        cast_upper95 = 0.0
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
            observations.get('reset_prior_safety_stop_reason_code') in ALLOWED_STICKY_S9_RESET_REASON_CODES
            or observations.get('reset_prior_safety_stop_reason') in ALLOWED_STICKY_S9_RESET_REASON_CODES
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
        any(observations.get(_mk) is not None for _mk in STAGE0_S10_V3_METADATA_KEYS)
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
        if (
            observations.get('cd16_adcc_igg_depletion_paired_control_verified') is not True
            or 'obs_exact_product_net_cd16_adcc_lysis_pct' not in validated_stage0
            or 'obs_cd16_adcc_positive_control_lysis_pct' not in validated_stage0
            or _upper95('obs_exact_product_net_cd16_adcc_lysis_pct') > 5.0
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
        or observations.get('iapp_age_matched_reference_above_loq_verified') is not None
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

    resolved_arch = dict(tree_doc.get('default_architecture_assignments', {}))
    if product_def is not None:
        _pvid_canon = str(product_def['product_variant_id'])
        resolved_arch['product_variant_id'] = _pvid_canon
        resolved_arch['product_execution_supported'] = bool(product_def.get('execution_supported', False))
        if _pvid_canon.startswith('ISLET1_B_SCT_'):
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
        and _lower95('obs_icasp9_log10_kill') >= 4.0
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
    elif str(resolved_arch.get('lot_maturation_status', '')).startswith('NO_GO_'):
        lot_release_readiness_tier = 'NO_GO_TERMINAL_BIOTRANSPORT_OR_IMMUNE_CEILING_EXCEEDED'
    elif resolved_arch.get('lot_maturation_status') not in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY') or any_resolved_component_on_hold:
        lot_release_readiness_tier = 'HOLD_FAILED_OR_INCOMPLETE_GATE'
    elif is_comp_r4_product and combined_stress_coupon_status == 'RETAIN_EMBODIMENT_A_STAR_4143BP':
        lot_release_readiness_tier = 'CANDIDATE_NOT_INDICATED_RETAIN_EMBODIMENT_A_STAR_4143BP'
    elif _prod_non_executable and _pvid_for_tier == 'ISLET1_B2_5LOCUS_BASE':
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
    }
    return {
        'schema_version': tree_doc['schema_version'],
        'strict_precedence_order': tree_doc['strict_precedence_order'],
        'mandatory_s9_safety_panel_complete': mandatory_s9_safety_panel_complete,
        'full_5endpoint_s9_panel_complete': full_5endpoint_s9_panel_complete,
        's9_5endpoint_er_safety_panel_complete': full_5endpoint_s9_panel_complete,
        's9_causal_failed_categories': s9_causal_failed_categories,
        's9_panel_completeness_tier': s9_panel_completeness_tier,
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
        'switch_decisions': decisions,
        'canonical_json_switch_id_map': canonical_json_switch_id_map,
        'pivotal_combined_stress_coupon_gate': combined_stress_spec,
        **(
            {'iip0026_molecular_provenance_evaluation': iip0026_eval_summary}
            if iip0026_eval_summary is not None
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
    return math.isfinite(float(val))


def _deep_freeze_iip0024(obj: Any) -> Any:
    if isinstance(obj, MappingProxyType):
        return obj
    if isinstance(obj, dict):
        return MappingProxyType({str(k): _deep_freeze_iip0024(v) for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))})
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
        'mcb_tier6b_pure_hpsc_clonogenic_escape': {'cap': 1e-8, 'require_observed_zero': True, 'use_poisson_protocol_floor': True},
        'obs_20q11_bcl2l1_single_cell_mosaic': {'cap': 0.01, 'require_observed_zero': False},
        'obs_1q32_mdm4_single_cell_mosaic': {'cap': 0.01, 'require_observed_zero': False},
        'obs_assembled_cluster_ki67_progenitor': {'cap': 0.005, 'require_observed_zero': False},
        'obs_ductal_sox9_krt19_double_pos': {'cap': 0.001, 'require_observed_zero': False},
        'obs_assembled_cluster_slc18a1_enterochromaffin': {'cap': 0.01, 'require_observed_zero': False},
        'obs_senescent_p16_sasp_positive': {'cap': 0.01, 'require_observed_zero': False},
        'obs_gate0v_day28_persistent_vp1_ie1': {'cap': 0.0005, 'require_observed_zero': True},
        'lot_tier6b_sorted_prolif_clonogenic_escape': {'cap': 3.0e-7, 'require_observed_zero': True},
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
        campaign_id: str = 'CAMP-DEFAULT-01',
        timepoint: str = 'STAGE0',
        sample_unit_id: str = 'SAMPLE-01',
        **forbidden_kwargs: Any,
    ) -> None:
        if forbidden_kwargs:
            raise ValueError(
                f'Caller-supplied projection or threshold overrides are prohibited on AtomicObservation: '
                f'{sorted(forbidden_kwargs.keys())} (IIP-0024; DEC-103)'
            )
        if not isinstance(endpoint_id, str) or not endpoint_id.strip():
            raise ValueError('endpoint_id must be a nonempty string (IIP-0024; DEC-103)')
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
        ctrl_dict = dict(controls or {})
        for bad_ck in ('canonical_projection', 'canonical_projection_keys', 'override_cap', 'override_requires_k_eq_0', 'cap'):
            if bad_ck in ctrl_dict:
                raise ValueError(f'Caller-supplied {bad_ck!r} inside controls is prohibited; caps derive from ENDPOINT_REGISTRY_SPEC (IIP-0024; DEC-103)')
        frozen_val = _deep_freeze_iip0024(value_or_k)
        frozen_ctrl = _deep_freeze_iip0024(ctrl_dict)
        for s_name, s_val in (
            ('unit', unit),
            ('assay_version', assay_version),
            ('campaign_id', campaign_id),
            ('timepoint', timepoint),
            ('sample_unit_id', sample_unit_id),
        ):
            if not isinstance(s_val, str) or not s_val.strip():
                raise ValueError(f'{s_name} must be a nonempty string (IIP-0024; DEC-103)')
        canonical_doc = {
            'registry_version': ENDPOINT_REGISTRY_SPEC['registry_version'],
            'endpoint_id': endpoint_id.strip(),
            'scope_type': scope_type,
            'scope_key': list(scope_key.as_tuple()),
            'kind': kind,
            'value_or_k': _deep_to_jsonable_iip0024(frozen_val),
            'se_or_n': se_or_n,
            'sample_n': sample_n,
            'unit': unit.strip(),
            'controls': _deep_to_jsonable_iip0024(frozen_ctrl),
            'assay_version': assay_version.strip(),
            'campaign_id': campaign_id.strip(),
            'timepoint': timepoint.strip(),
            'sample_unit_id': sample_unit_id.strip(),
        }
        digest = hashlib.sha256(
            json.dumps(canonical_doc, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')
        ).hexdigest()
        object.__setattr__(self, 'endpoint_id', endpoint_id.strip())
        object.__setattr__(self, 'scope_type', scope_type)
        object.__setattr__(self, 'scope_key', scope_key)
        object.__setattr__(self, 'kind', kind)
        object.__setattr__(self, 'value_or_k', frozen_val)
        object.__setattr__(self, 'se_or_n', se_or_n)
        object.__setattr__(self, 'sample_n', sample_n)
        object.__setattr__(self, 'unit', unit.strip())
        object.__setattr__(self, 'controls', frozen_ctrl)
        object.__setattr__(self, 'assay_version', assay_version.strip())
        object.__setattr__(self, 'campaign_id', campaign_id.strip())
        object.__setattr__(self, 'timepoint', timepoint.strip())
        object.__setattr__(self, 'sample_unit_id', sample_unit_id.strip())
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
    eta_spike_recovery_l95: Any = 1.0,
    f_clonogenic: Any = 1.0,
    use_poisson_protocol_floor: bool = False,
) -> dict[str, Any]:
    """
    Fail-closed rare-event integer count validator (IIP-0024; DEC-103):
      - Rejects bool, float (e.g. 298.9 or 300.0), NaN, Inf, < 0, k > n with structured HOLD.
      - Enforces k == 0 observed-zero hard veto where specified even when U95(k, n) <= cap.
      - Computes exact independent binomial-thinning sensitivity U95 = (1 - 0.05**(1/n)) / (f_clono * eta)
        (holding 0/920 at eta = 0.325 where U95 = 0.01000288 > 0.01, never floor(n * eta)).
      - Rejects unregistered cluster estimators and naive DEFF substitution.
    """
    if k is None or n is None:
        return {'status': 'HOLD_MISSING_ATOMIC_INTEGER_K_OR_N', 'gate_id': gate_id}
    if not _is_strict_int_iip0024(k) or not _is_strict_int_iip0024(n):
        return {'status': 'HOLD_NON_INTEGER_COUNT_OR_DENOMINATOR_REJECTED', 'gate_id': gate_id}
    if n <= 0 or k < 0 or k > n:
        return {'status': 'HOLD_INVALID_COUNT_BOUNDS', 'gate_id': gate_id}
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
                'status': 'HOLD_CLUSTERED_OR_UNREGISTERED_ESTIMATOR_REJECTED',
                'gate_id': gate_id,
            }
    if not _is_strict_finite_number_iip0024(eta_spike_recovery_l95) or not _is_strict_finite_number_iip0024(f_clonogenic):
        return {'status': 'HOLD_NONFINITE_RECOVERY_OR_CLONOGENIC_FACTOR', 'gate_id': gate_id}
    eta = float(eta_spike_recovery_l95)
    f_clono = float(f_clonogenic)
    if not (0.0 < eta <= 1.0) or not (0.0 < f_clono <= 1.0):
        return {'status': 'HOLD_OUT_OF_RANGE_RECOVERY_OR_CLONOGENIC_FACTOR', 'gate_id': gate_id}
    if (eta < 1.0 or f_clono < 1.0) and eta < 0.30:
        return {'status': 'HOLD_SPIKE_RECOVERY_BELOW_0P30_QUALIFICATION_FLOOR', 'gate_id': gate_id}
    if require_observed_zero and k > 0:
        return {
            'status': 'FAIL_HARD_STOP',
            'reason': 'OBSERVED_POSITIVE_EVENT_VIOLATES_MANDATORY_ZERO_VETO',
            'gate_id': gate_id,
            'k': k,
            'n_nominal': n,
            'u95_nominal': round(clopper_pearson_upper95_independent(k, n), 10),
        }
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
        dev_k = AssembledDeviceLotScopeKey(
            mcb_k.product_variant_id,
            mcb_k.construct_sha256,
            mcb_k.donor_seedline_id,
            mcb_k.mcb_id,
            lot_k.lot_id,
            target_provenance['device_assembly_id'],
            target_provenance['physical_route_id'],
        )
    except (KeyError, TypeError, ValueError) as exc:
        return {'readiness_tier': 'HOLD_MISSING_OR_INVALID_IMMUTABLE_PROVENANCE_KEY', 'holds': [str(exc)]}

    try:
        pdef = p_mod.get_product_definition(mcb_k.product_variant_id, root=root_path)
    except ValueError as exc:
        return {'readiness_tier': 'HOLD_UNKNOWN_PRODUCT_VARIANT', 'holds': [str(exc)]}

    pvid_canon = pdef['product_variant_id']
    if pvid_canon == 'ISLET1_B2_5LOCUS_BASE':
        return {
            'readiness_tier': 'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC',
            'holds': ['HOLD_STAGED_B2_5LOCUS_BASE_EXECUTION_SUPPORTED_FALSE'],
            'execution_supported': False,
        }
    if pdef.get('execution_supported') is False and pvid_canon != 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2':
        return {
            'readiness_tier': 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD',
            'holds': [f'non_executable_product:{pvid_canon}'],
            'execution_supported': False,
        }

    routes = ENDPOINT_REGISTRY_SPEC['physical_routes']
    if dev_k.physical_route_id not in routes:
        return {'readiness_tier': 'HOLD_UNREGISTERED_PHYSICAL_ROUTE_ID', 'holds': [dev_k.physical_route_id]}
    route_spec = routes[dev_k.physical_route_id]
    is_z_prod = (pvid_canon == 'ISLET1_Z_X_MOL_140UM_ZERO_CAO2')
    if is_z_prod and not route_spec['requires_design_z_two_stage']:
        return {'readiness_tier': 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH', 'holds': [dev_k.physical_route_id]}
    if not is_z_prod and route_spec['requires_design_z_two_stage']:
        return {'readiness_tier': 'HOLD_BIDIRECTIONAL_PRODUCT_PHYSICAL_ROUTE_MISMATCH', 'holds': [dev_k.physical_route_id]}

    mcb_store: dict[str, AtomicObservation] = {}
    lot_store: dict[str, AtomicObservation] = {}
    dev_store: dict[str, AtomicObservation] = {}
    atom_status_by_store_key: dict[tuple[str, str], str] = {}
    invalidated_stale_keys: set[str] = set()
    latched_adverse: dict[tuple[str, tuple[str, ...]], dict[str, Any]] = {}
    retained_reset_audit_log: list[dict[str, Any]] = []
    structural_holds: list[str] = []

    count_caps = ENDPOINT_REGISTRY_SPEC['count_endpoint_caps_and_zero_vetoes']
    sister_terminal_ids = set(ENDPOINT_REGISTRY_SPEC['device_sister_terminal_endpoint_ids'])
    target_scope_tuples = {mcb_k.as_tuple(), lot_k.as_tuple(), dev_k.as_tuple()}

    for b_idx, batch in enumerate(batches):
        atoms = batch.get('atoms')
        if not isinstance(atoms, list):
            return {'readiness_tier': 'HOLD_MALFORMED_BATCH_ATOM_LIST', 'holds': ['atoms_must_be_list']}
        for atom in atoms:
            if not isinstance(atom, AtomicObservation):
                return {'readiness_tier': 'HOLD_NON_ATOMIC_OBSERVATION_OBJECT', 'holds': ['raw_dict_rejected']}
            sk_tup = atom.scope_key.as_tuple()
            # Evaluate atom validity and latch any adverse failure on its exact scope:
            atom_status = 'PASS'
            if atom.kind == 'EXACT_INTEGER_COUNT':
                spec_c = count_caps.get(atom.endpoint_id, {'cap': 0.01, 'require_observed_zero': False})
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
                    eta_spike_recovery_l95=ctrl.get('eta_spike_recovery_l95', 1.0),
                    f_clonogenic=ctrl.get('f_clonogenic', 1.0),
                    use_poisson_protocol_floor=bool(spec_c.get('use_poisson_protocol_floor', False)),
                )
                atom_status = c_res['status']
            elif atom.kind == 'STRUCTURED_PANEL':
                p_dict = dict(atom.value_or_k)
                if 'imprinting_' in atom.endpoint_id or 'meg3_dmr' in atom.endpoint_id:
                    atom_status = evaluate_conjunctive_imprinting_dmr_iip0024(endpoint_scope=atom.scope_type, **p_dict)['status']
                elif 'xci_46xx' in atom.endpoint_id:
                    st_x = evaluate_conjunctive_xci_46xx_iip0024(endpoint_scope=atom.scope_type, **p_dict)['status']
                    atom_status = 'PASS' if st_x in ('PASS', 'NOT_APPLICABLE_46XY') else st_x
                elif 'mtdna_1000x' in atom.endpoint_id:
                    atom_status = evaluate_scoped_mtdna_gate_iip0024(endpoint_scope=atom.scope_type, **p_dict)['status']
                else:
                    atom_status = 'HOLD_UNKNOWN_STRUCTURED_PANEL'

            latch_key = (atom.endpoint_id, sk_tup)
            if atom_status == 'FAIL_HARD_STOP':
                latched_adverse[latch_key] = {'stop_idx': b_idx, 'repeats': {}}
            elif latch_key in latched_adverse and atom_status == 'PASS' and b_idx > latched_adverse[latch_key]['stop_idx']:
                if dict(atom.controls).get('washout_and_retest_verified') is True and atom.campaign_id == target_campaign_id:
                    latched_adverse[latch_key]['repeats'][str(batch.get('batch_id', b_idx))] = atom.sample_unit_id

            if atom.scope_type == 'MCB_SCOPE':
                if atom.scope_key == mcb_k:
                    mcb_store[atom.endpoint_id] = atom
                    atom_status_by_store_key[('MCB_SCOPE', atom.endpoint_id)] = atom_status
            elif atom.scope_type == 'LOT_SCOPE':
                if atom.scope_key == lot_k:
                    if atom.campaign_id != target_campaign_id:
                        structural_holds.append(f'HOLD_UNLINKED_CAMPAIGN:{atom.endpoint_id}')
                    else:
                        lot_store[atom.endpoint_id] = atom
                        atom_status_by_store_key[('LOT_SCOPE', atom.endpoint_id)] = atom_status
                elif atom.scope_key.to_mcb_key() == mcb_k:
                    invalidated_stale_keys.add(atom.endpoint_id)
            elif atom.scope_type == 'ASSEMBLED_DEVICE_LOT_SCOPE':
                same_or_sister = (
                    atom.scope_key.to_lot_key() == lot_k
                    and atom.scope_key.physical_route_id == dev_k.physical_route_id
                    and (
                        atom.scope_key.device_assembly_id == dev_k.device_assembly_id
                        or atom.scope_key.device_assembly_id in authorized_sister_assembly_ids
                    )
                )
                if same_or_sister:
                    if atom.scope_key.device_assembly_id != dev_k.device_assembly_id and atom.endpoint_id not in sister_terminal_ids:
                        structural_holds.append(f'HOLD_SISTER_UNIT_NON_TERMINAL_TRANSFER_FORBIDDEN:{atom.endpoint_id}')
                    elif atom.campaign_id != target_campaign_id:
                        structural_holds.append(f'HOLD_UNLINKED_CAMPAIGN:{atom.endpoint_id}')
                    else:
                        dev_store[atom.endpoint_id] = atom
                        atom_status_by_store_key[('ASSEMBLED_DEVICE_LOT_SCOPE', atom.endpoint_id)] = atom_status
                elif atom.scope_key.to_mcb_key() == mcb_k:
                    invalidated_stale_keys.add(atom.endpoint_id)

        if batch.get('reset_prior_safety_stop_after_investigation') is True:
            req_eid = batch.get('reset_prior_safety_stop_target_category')
            rep_ids = batch.get('reset_prior_safety_stop_repeat_batch_ids') or []
            qp_id = str(batch.get('reset_prior_safety_stop_independent_qp_approver_id') or '').strip()
            sub_id = str(batch.get('submitter_id') or '').strip()
            r_code = batch.get('reset_prior_safety_stop_reason_code')
            req_scope = batch.get('reset_prior_safety_stop_target_scope_tuple')
            allowed_scopes_for_reset = {tuple(req_scope)} if isinstance(req_scope, (list, tuple)) else target_scope_tuples
            if (
                isinstance(req_eid, str)
                and qp_id
                and sub_id
                and qp_id != sub_id
                and r_code in ALLOWED_STICKY_S9_RESET_REASON_CODES
                and isinstance(rep_ids, list)
                and len(set(rep_ids)) >= 2
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

    holds: list[str] = list(structural_holds)
    for (_sc_t, _eid_t), _st_t in sorted(atom_status_by_store_key.items()):
        if _st_t not in ('PASS', 'FAIL_HARD_STOP'):
            holds.append(f'{_st_t}:{_eid_t}')

    active_target_latches = [eid for (eid, skt) in latched_adverse if skt in target_scope_tuples]
    if active_target_latches:
        return {
            'readiness_tier': 'HARD_SAFETY_STOP_REJECTED',
            'holds': [f'latched_adverse:{e}' for e in sorted(active_target_latches)] + holds,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
        }

    if 'gate0s_shedding_ddpcr_assay' in dev_store and 'gate0s_shedding_viable_outgrowth_assay' in dev_store:
        if dev_store['gate0s_shedding_ddpcr_assay'].sample_unit_id == dev_store['gate0s_shedding_viable_outgrowth_assay'].sample_unit_id:
            holds.append('HOLD_GATE0S_DUAL_SHEDDING_REQUIRES_INDEPENDENT_SPECIMEN_UNITS')

    missing_mcb = [mk for mk in ENDPOINT_REGISTRY_SPEC['mcb_mandatory_endpoint_ids'] if mk not in mcb_store]
    if holds:
        return {
            'readiness_tier': 'HOLD_FAILED_STAGE0_GATE',
            'holds': holds + missing_mcb,
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'retained_reset_audit_log': retained_reset_audit_log,
        }
    if missing_mcb or not lot_store or (specimen_lifecycle_stage != 'AUTHORIZED_TEST_ARTICLE_PRE_IMPLANT' and any(tk not in dev_store for tk in sister_terminal_ids)):
        tier = (
            'HOLD_LOT_OR_DEVICE_SCOPED_EVIDENCE_MISSING_AFTER_PROVENANCE_CHANGE'
            if invalidated_stale_keys
            else 'PROVISIONAL_PARTIAL_LOT_SCREENING_ONLY'
        )
        return {
            'readiness_tier': tier,
            'holds': missing_mcb or ['incomplete_lot_or_device_panel'],
            'invalidated_stale_keys': sorted(invalidated_stale_keys),
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
        }
    if specimen_lifecycle_stage == 'AUTHORIZED_TEST_ARTICLE_PRE_IMPLANT':
        return {
            'readiness_tier': 'AUTHORIZED_TEST_ARTICLE_ASSEMBLED_PENDING_LONGITUDINAL_READOUT',
            'holds': [],
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
        }
    if is_z_prod:
        return {
            'readiness_tier': 'WILDCARD_COUPON_QUALIFIED_NON_EXECUTABLE',
            'holds': [],
            'execution_supported': False,
            'wildcard_coupon_qualified': True,
            'transferred_mcb_count': len(mcb_store),
            'retained_reset_audit_log': retained_reset_audit_log,
        }
    return {
        'readiness_tier': (
            'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY'
            if retained_reset_audit_log
            else 'READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE'
        ),
        'holds': [],
        'execution_supported': True,
        'transferred_mcb_count': len(mcb_store),
        'retained_reset_audit_log': retained_reset_audit_log,
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
            (
                'EXECUTION_SUPPORTED_QUALIFIED_PRODUCT'
                if exec_sup
                else (
                    'STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC'
                    if pid in ('ISLET1_B2_5LOCUS_BASE', 'ISLET1_B2_5LOCUS_TYK2_PE')
                    else 'RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD'
                )
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


def _extract_per_locus_baf_dict_iip0026(payload: dict[str, Any]) -> dict[str, Any] | None:
    raw_dict = payload.get('post_edit_per_locus_segmented_baf')
    if isinstance(raw_dict, dict):
        return dict(raw_dict)
    if raw_dict is not None:
        return None
    assembled: dict[str, Any] = {}
    for loc, ok in (
        ('chr19q', 'obs_post_edit_segmented_baf_chr19q'),
        ('chr11p15', 'obs_post_edit_segmented_baf_chr11p15'),
        ('chr6p21', 'obs_post_edit_segmented_baf_chr6p21'),
    ):
        if ok in payload and payload[ok] is not None and payload[ok] != 'UNAVAILABLE':
            assembled[loc] = payload[ok]
    return assembled if assembled else None


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
    alias_pairs = (
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
    for short_k, obs_k, req_int in alias_pairs:
        raw_v = payload.get(obs_k, payload.get(short_k))
        if raw_v is None:
            continue
        if req_int:
            if not _is_strict_int_iip0024(raw_v):
                continue
        elif not _is_strict_finite_number_iip0024(raw_v):
            continue
        lo, hi = valid_mean_se[obs_k]
        fv = float(raw_v)
        if lo <= fv <= hi:
            out[obs_k] = int(raw_v) if req_int else fv
    for k, v in payload.items():
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
    reg_by_id = {
        str(p['product_variant_id']): _derive_versioned_product_registry_entry_iip0026(p)
        for p in prod_doc['products']
    }

    pid = payload.get('product_variant_id')
    if not isinstance(pid, str) or not pid.strip() or pid != pid.strip() or pid not in reg_by_id:
        empty_s0 = evaluate_stage0_decision_tree({}, root=root_path)
        return {
            'overall_status': 'HOLD_UNREGISTERED_OR_MISSING_PRODUCT_VARIANT_ID',
            'lot_release_readiness_tier': 'HOLD_UNREGISTERED_OR_MISSING_PRODUCT_VARIANT_ID',
            'iip0026_reasons': ['HOLD_UNREGISTERED_OR_MISSING_PRODUCT_VARIANT_ID'],
            'lifecycle_state': 'UNREGISTERED',
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': empty_s0['switch_decisions'],
        }

    reg = reg_by_id[pid]
    lifecycle_state = reg['lifecycle_state']
    registry_exec_supported = bool(reg['execution_supported'])
    safe_s0_obs = _build_safe_stage0_tree_obs_iip0026(payload, root_path=root_path)
    tree_eval = evaluate_stage0_decision_tree(safe_s0_obs, root=root_path)
    switch_decisions = dict(tree_eval['switch_decisions'])

    # Reject any caller attempt to supply or override lifecycle_state or execution_supported:
    if any(
        (k in IIP0026_FORBIDDEN_CALLER_OVERRIDE_KEYS) or str(k).endswith('_lifecycle_state') or str(k).endswith('_execution_supported')
        for k in payload
    ):
        return {
            'overall_status': 'FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED',
            'lot_release_readiness_tier': 'FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED',
            'iip0026_reasons': ['FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED'],
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
        }

    lot_id = payload.get('lot_id', payload.get('current_lot_id'))
    campaign_id = payload.get('campaign_id', payload.get('assay_campaign_id'))
    if (
        not isinstance(lot_id, str)
        or not lot_id.strip()
        or not isinstance(campaign_id, str)
        or not campaign_id.strip()
    ):
        return {
            'overall_status': 'HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE',
            'lot_release_readiness_tier': 'HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE',
            'iip0026_reasons': ['HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE'],
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
        }

    if reg['bound_c1_c5_pair_sha256'] is not None:
        supplied_sha = payload.get('bound_c1_c5_pair_sha256')
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
                'lifecycle_state': lifecycle_state,
                'execution_supported': False,
                'effective_execution_permitted': False,
                'switch_decisions': switch_decisions,
            }

    missing_obs = [k for k in IIP0026_REQUIRED_CORE_S9_OBS_KEYS if not _is_strict_finite_number_iip0024(payload.get(k))]
    if (
        missing_obs
        or float(payload['obs_tp53_ppm1d_bcor_duplex_depth']) <= 0.0
        or float(payload['obs_castseq_duplex_or_amplicon_depth']) <= 0.0
        or not _is_integral_count_iip0026(payload['obs_mcb_clonogenic_survivor_count_k'])
        or not _is_integral_count_iip0026(payload['obs_mcb_effective_clonogenic_trials'])
        or int(payload['obs_mcb_effective_clonogenic_trials']) <= 0
    ):
        return {
            'overall_status': 'HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS',
            'lot_release_readiness_tier': 'HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS',
            'iip0026_reasons': ['HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS'],
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
        }

    if tree_eval['safety_kill_s9_triggered']:
        return {
            'overall_status': 'HARD_SAFETY_STOP_REJECTED',
            'lot_release_readiness_tier': 'HARD_SAFETY_STOP_REJECTED',
            'iip0026_reasons': ['HARD_SAFETY_STOP_S9'],
            'lifecycle_state': lifecycle_state,
            'execution_supported': False,
            'effective_execution_permitted': False,
            'switch_decisions': switch_decisions,
        }

    reasons: list[str] = []
    sera_val = payload.get('iip0026_schema_era', payload.get('schema_era'))
    if isinstance(sera_val, str) and sera_val.startswith('LEGACY_PRE_'):
        reasons.append('HOLD_PROVISIONAL_LEGACY_SCHEMA_MISSING_IIP0026_EVIDENCE')

    if payload.get('donor_seed_gnomad_target_site_wgs_verified') is not True:
        reasons.append('HOLD_MISSING_DONOR_TARGET_WGS')

    per_locus_baf = _extract_per_locus_baf_dict_iip0026(payload)
    sv_recip = payload.get('post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified')
    req_loci = reg['required_baf_loci']
    if (
        not isinstance(per_locus_baf, dict)
        or any(loc not in per_locus_baf or not _is_strict_finite_number_iip0024(per_locus_baf[loc]) for loc in req_loci)
        or sv_recip is not True
    ):
        reasons.append('HOLD_MISSING_SEGMENTED_PHASED_BAF_SURVEILLANCE')
    else:
        max_abs_locus_baf_shift = max(abs(float(per_locus_baf[loc]) - 0.50) for loc in req_loci)
        if max_abs_locus_baf_shift > 0.05:
            reasons.append('FAIL_DISCORDANT_PER_LOCUS_SEGMENTED_BAF')

    if reg['requires_endogenous_hlac_and_locus5_gates']:
        cfd_max = payload.get('obs_design_b_hlac_phased_guide_cfd_max', payload.get('design_b_hlac_phased_guide_cfd_max'))
        sep_pass = payload.get('design_b_hla_ab_separated_pass_delivery_verified')
        chr6_loh = payload.get('design_b_chr6p21_1p41mb_loh_sv_cleared')
        if (
            not _is_strict_finite_number_iip0024(cfd_max)
            or not (0.0 <= float(cfd_max) < 0.05)
            or sep_pass is not True
            or chr6_loh is not True
        ):
            reasons.append('HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE')

        dt9_pos = payload.get('design_b_hlac_dt9_positive')
        hlae_3d12_pos = payload.get('design_b_hlae_3d12_positive')
        dt9_delta = payload.get(
            'obs_design_b_hlac_dt9_matched_hlae_only_control_delta_ratio',
            payload.get('design_b_hlac_dt9_matched_hlae_only_control_delta_ratio'),
        )
        hlac_orth = payload.get('design_b_hlac_allele_specific_surface_or_phased_mrna_verified')
        if dt9_pos is not True or hlae_3d12_pos is not True:
            if 'HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE' not in reasons:
                reasons.append('HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE')
        elif dt9_delta is None and hlac_orth is not True:
            reasons.append('HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY')
        elif dt9_delta is not None:
            if not _is_strict_finite_number_iip0024(dt9_delta) or float(dt9_delta) <= 0.0:
                reasons.append('HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY')
            elif float(dt9_delta) < 2.0 and hlac_orth is not True:
                reasons.append('FAIL_HLAC_SURFACE_LOSS_MASKED_BY_DT9_HLAE_CROSSREACTIVITY')

        r5_bp = payload.get('obs_design_b_locus5_5prime_race_coord_bp', payload.get('design_b_locus5_5prime_race_coord_bp'))
        r3_bp = payload.get('obs_design_b_locus5_3prime_race_coord_bp', payload.get('design_b_locus5_3prime_race_coord_bp'))
        max_donor_bp = 4078 if reg['bound_c1_c5_pair_sha256'] == DESIGN_B_C1_C5_PAIR_SHA256 else 10000
        if r5_bp is None or r3_bp is None:
            reasons.append('HOLD_MISSING_DESIGNB_LOCUS5_RACE_COORDINATES')
        elif not (
            _is_strict_int_iip0024(r5_bp)
            and _is_strict_int_iip0024(r3_bp)
            and 0 < r5_bp < r3_bp <= max_donor_bp
        ):
            reasons.append('FAIL_INVALID_RACE_COORDINATES')

        pdl1_r5 = payload.get('design_b_locus5_pdl1_5prime_race_coord_bp')
        if pdl1_r5 == 'UNRESOLVED':
            reasons.append('HOLD_PENDING_5PRIME_RACE_TSS')
        elif pdl1_r5 is not None:
            if not _is_strict_int_iip0024(pdl1_r5) or pdl1_r5 <= 0 or pdl1_r5 > 2381:
                reasons.append('FAIL_INVALID_RACE_COORDINATES')
            elif pdl1_r5 <= 2308:
                reasons.append('FAIL_MATURE_TRANSCRIPT_INCLUDES_UPSTREAM_OUT_OF_FRAME_KOZAK_ORF')

        b_umi_n = payload.get('obs_design_b_locus5_race_umi_n_eff', payload.get('design_b_locus5_race_umi_n_eff'))
        b_umi_k = payload.get('obs_design_b_locus5_aberrant_umi_k', payload.get('design_b_locus5_aberrant_umi_k'))
        if not _is_strict_int_iip0024(b_umi_n) or b_umi_n < 3689 or not _is_strict_int_iip0024(b_umi_k) or b_umi_k < 0:
            reasons.append('HOLD_MISSING_DESIGNB_LOCUS5_TRANSCRIPT_UMI_3689')
        elif b_umi_k > 0:
            reasons.append('FAIL_DESIGNB_LOCUS5_CHIMERIC_TRANSCRIPT_DETECTED')
    elif reg['requires_b_sct_gates']:
        if payload.get('b_sct_ciita_exon5_and_trimer_verified') is not True:
            reasons.append('HOLD_MISSING_BSCT_CIITA_EXON5_AND_TRIMER_EVIDENCE')

    if reg['requires_tyk2_pe_gates']:
        pe_mod = payload.get('tyk2_pe_modality')
        if pe_mod not in ALLOWED_TYK2_PE_MODALITIES or payload.get('tyk2_phased_haplotype_sap_verified') is not True:
            reasons.append('HOLD_MISSING_TYK2_PE_MODALITY_OR_PHASED_HAPLOTYPE')
        elif pe_mod == 'PE4MAX_MSI_GATED':
            msi_k = payload.get('obs_tyk2_pe4max_pentaplex_msi_shifts_k', payload.get('tyk2_pe4max_pentaplex_msi_shifts_k'))
            if msi_k is None or not _is_strict_int_iip0024(msi_k) or msi_k < 0:
                reasons.append('HOLD_MISSING_MSI')
            elif msi_k > 0:
                reasons.append('FAIL_MSI_DETECTED')
            pe_st = evaluate_pe4max_piecewise_repeat_indel_cap_iip0026(
                payload.get('obs_tyk2_pe4max_wgs_repeat_indels_clone', payload.get('tyk2_pe4max_wgs_repeat_indels_clone')),
                payload.get('obs_tyk2_pe4max_wgs_repeat_indels_ctrl', payload.get('tyk2_pe4max_wgs_repeat_indels_ctrl')),
                payload.get('obs_tyk2_pe4max_coding_or_splice_de_novo_indels', payload.get('tyk2_pe4max_coding_or_splice_de_novo_indels')),
            )
            if pe_st != 'PASS':
                if not (msi_k is None and pe_st == 'HOLD_MISSING_PE4MAX_REPEAT_INDEL_EVIDENCE'):
                    reasons.append(pe_st)

    if reg['is_held_design_y']:
        reasons.append('HOLD_DESIGNY_EXPLORATORY_UNBIASED_CAPTURE_AND_63K_OR_90K_LOD_REQUIRED')

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
        'lifecycle_state': lifecycle_state,
        'execution_supported': registry_exec_supported if len(reasons) == 0 else False,
        'effective_execution_permitted': effective_exec,
        'endogenous_hlac_gate_applied': reg['requires_endogenous_hlac_and_locus5_gates'],
        'b_sct_route_separated': reg['requires_b_sct_gates'],
        'switch_decisions': switch_decisions,
    }
