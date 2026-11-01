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
]

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
        if not isinstance(product_variant_id, str) or not product_variant_id.strip():
            raise ValueError('product_variant_id must be a nonempty string')
        p_mod = _load_product_contract_module(root_path)
        product_def = p_mod.get_product_definition(product_variant_id.strip(), root=root_path)
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


def evaluate_stage0_decision_tree(
    observations: dict,
    *,
    root: Path | str | None = None,
    require_full_5endpoint_s9_panel: bool = False,
    defaulted_se_keys: list[str] | set[str] | None = None,
    prior_batch_sticky_s9_safety_stop: bool = False,
    product_def: dict | None = None,
) -> dict:
    """Deterministically evaluate the machine-readable Stage-0 promotion/fallback decision tree (S1..S9; AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, & R7S3-F01/F02; DEC-058)."""
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
    if product_def is None and 'product_variant_id' in observations and observations['product_variant_id'] is not None:
        product_def = p_mod.get_product_definition(str(observations['product_variant_id']).strip(), root=root_path)

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
    }
    unknown_keys = sorted(set(observations) - set(spec_by_mean) - set(spec_by_se) - allowed_boolean_and_meta_flags)
    if unknown_keys:
        raise ValueError(f'Unrecognized Stage-0 decision-tree observation keys: {unknown_keys}')
    if 'prior_variant_safety_status' in observations and observations['prior_variant_safety_status'] is not None:
        pvs_val = observations['prior_variant_safety_status']
        if not isinstance(pvs_val, str) or pvs_val.strip() not in ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES:
            raise ValueError(
                f"Invalid prior_variant_safety_status {pvs_val!r}; must be one of {sorted(ALLOWED_PRIOR_VARIANT_SAFETY_STATUSES)} (AUDIT-REF; DEC-071)"
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

    def _upper95(mean_key: str, default: float = 0.0) -> float:
        if mean_key not in validated_stage0:
            return default
        return float(validated_stage0[mean_key]) + 1.644854 * _validated_se(mean_key)

    def _lower95(mean_key: str, default: float = 0.0) -> float:
        if mean_key not in validated_stage0:
            return default
        return float(validated_stage0[mean_key]) - 1.644854 * _validated_se(mean_key)

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
        and isinstance(verified_ledger_passing, (list, tuple, set))
        and {x.strip() for x in repeat_batch_ids}.issubset({str(b).strip() for b in verified_ledger_passing if str(b).strip()})
    )
    sticky_prior_s9_stop = bool(prior_s9_flag and not structured_s9_reset_valid)
    prior_s9_reset_history_retained = bool(prior_s9_flag and structured_s9_reset_valid)

    # Rare-event trial denominators & exact Poisson/Clopper-Pearson upper95 bounds for residual PSC, hazard-state proliferative survivor, Step-1 10 nM AP1903 survivor, and 6 proportion/percentage gates (AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF; DEC-050/051/071):
    def _prop_or_pct_u95(mean_key: str, n_key: str, is_pct: bool = True) -> float:
        w_u95 = _upper95(mean_key)
        if mean_key in validated_stage0 and n_key in validated_stage0 and validated_stage0[n_key] > 0.0:
            n_c = validated_stage0[n_key]
            scale = 100.0 if is_pct else 1.0
            k_ev = (validated_stage0[mean_key] / scale) * n_c
            cp_rate = (
                (1.0 - (0.05 ** (1.0 / n_c)))
                if round(k_ev) <= 0
                else min(1.0, _poisson_exact_upper95_rate(k_ev, n_c))
            )
            return max(w_u95, cp_rate * scale)
        return w_u95

    zero_event_proportion_missing_denominator_or_explicit_se = any(
        mk in validated_stage0
        and validated_stage0[mk] == 0.0
        and nk not in validated_stage0
        and (
            se_provenance_by_mean.get(mk) != 'EXPLICIT_SE_SUPPLIED'
            or _validated_se(mk) < 5e-4
        )
        for mk, nk in (
            ('obs_20q11_bcl2l1_single_cell_mosaic_fraction', 'obs_20q11_bcl2l1_single_cell_mosaic_cells_scored'),
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
        ('obs_explant_en_bloc_recovery_fraction' in validated_stage0 and _lower95('obs_explant_en_bloc_recovery_fraction') < 0.99)
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
    if ap1903_mono_pk_failed:
        s9_causal_failed_categories.append('MONOTHERAPY_PULSED_PK_KILL')
    if ap1903_step1_sublethal_survivor_detected:
        s9_causal_failed_categories.append('STEP1_SUBLETHAL_SURVIVOR')
    if explant_retrieval_failed:
        s9_causal_failed_categories.append('EXPLANT_RECOVERY')

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
            and validated_stage0['obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct'] > 8.0
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
    kir_nkg2a_neg_nk_exceeded_raw = bool(
        validated_stage0.get('obs_kir_nkg2a_neg_nk_lysis_pct', 2.0) > 8.0
        or ('obs_il12_il18_sirpa_low_nk_lysis_pct' in validated_stage0 and validated_stage0['obs_il12_il18_sirpa_low_nk_lysis_pct'] > 10.0)
    )
    track2_kir_or_hlag1_fallback_promoted = bool(
        kir_nkg2a_neg_nk_exceeded_raw
        and observations.get('track2_b_sct_or_hlag1_lilrb1_rescue_qualified') is True
    )
    kir_nkg2a_neg_nk_failed = bool(
        kir_nkg2a_neg_nk_exceeded_raw
        and not track2_kir_or_hlag1_fallback_promoted
    )
    # S10 triggers on CMV-imprinted NKG2C+NKG2A-CD57+ adaptive NK lysis > 8% OR receptor-negative NKG2A-KIR2DL-NKG2C+ NK lysis > 8% OR covalent-SCT NKG2C/NKG2A KD selectivity ratio < 4.0 (AUDIT-REF & AUDIT-REF; DEC-042/046):
    s10_nkg2c_triggered = bool(
        validated_stage0.get('obs_nkg2c_adaptive_nk_lysis_pct', 2.0) > 8.0
        or receptor_neg_nkg2c_nk_failed
        or _lower95('obs_hlae_nkg2c_over_nkg2a_kd_ratio', 4.5) < 4.0
    )
    s10_coupon_cleared = bool(
        observations.get('s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct') is True
        and not receptor_neg_nkg2c_nk_failed
        and validated_stage0.get('obs_nkg2c_adaptive_nk_lysis_pct', 2.0) <= 8.0
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
    s7_promoted = bool(s7_allo_triggered and sct_density_explicitly_qualified and not b_prime_sct_direct_cd8_failed)
    s7_hold_unqualified_sct = bool(s7_allo_triggered and not sct_density_explicitly_qualified)
    s7_hold_unverified_sct_allo = bool((s7_allo_triggered and sct_density_explicitly_qualified and not sct_cd8_mlc_explicitly_cleared) or b_prime_sct_direct_cd8_failed)
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
    indirect_apc_cd4_failed = bool(
        ('obs_indirect_apc_dq2_dq8_cd4_stimulation_index' in validated_stage0 and _upper95('obs_indirect_apc_dq2_dq8_cd4_stimulation_index') > 1.5)
        or ('obs_indirect_apc_cd4_ifng_si' in validated_stage0 and _upper95('obs_indirect_apc_cd4_ifng_si') > 2.0)
        or ('obs_gate0d_activated_cd4_fasl_apoptosis_pct' in validated_stage0 and _lower95('obs_gate0d_activated_cd4_fasl_apoptosis_pct') < 80.0)
        or ('obs_gate0d_transwell_beta_annexin_v_pct' in validated_stage0 and _upper95('obs_gate0d_transwell_beta_annexin_v_pct') > 15.0)
    )

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
        ('obs_explant_en_bloc_recovery_fraction' in validated_stage0 and _lower95('obs_explant_en_bloc_recovery_fraction') < 0.99)
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

    resolved_arch = dict(tree_doc.get('default_architecture_assignments', {}))
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
    elif track2_kir_or_hlag1_fallback_promoted:
        resolved_arch['locus1_construct'] = 'PROMOTE_TRACK2_B_SCT_OR_HLAG1_LILRB1_RESCUE'
        resolved_arch['locus2_construct'] = 'ISLET1_B_SCT_C1_CIITA_OR_C2_CIITA_KIR_RESCUE'
        resolved_arch['track2_kir_subbank_construct'] = 'ISLET1_B_SCT_C1_OR_C2_CIITA_OR_HLAG1_C42S_CHALLENGER'
        resolved_arch['track1_mhc_null_kir_escape_status'] = 'PROMOTED_TO_TRACK2_B_SCT_OR_HLAG1_LILRB1_CHALLENGER'
    elif s2_promoted:
        resolved_arch['locus1_construct'] = 'ISLET1_A4_X_MOL_RNLS_WT'
        if observations.get('s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior') is not True:
            resolved_arch['locus1_arm1x_qualification_status'] = 'CANDIDATE_PROMOTED_REQUIRES_ARM1X_GATE0A_QUALIFICATION'
    elif cd47_or_phagocytosis_failed and not dxs_joint_retention_verified:
        resolved_arch['locus1_construct'] = 'HOLD_LOCUS1_CD47_DEFICIENT_PENDING_DXS_QUALIFICATION'
    elif unmitigated_cd47_or_hypoxic_retention_cliff_failed:
        resolved_arch['locus1_construct'] = 'HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED'
    if s10_nkg2c_triggered:
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
            str(product_def.get('family', 'A')).upper().startswith('B')
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
        'obs_indirect_apc_dq2_dq8_cd4_stimulation_index' in validated_stage0
        and 'obs_indirect_apc_cd4_ifng_si' in validated_stage0
        and 'obs_gate0d_activated_cd4_fasl_apoptosis_pct' in validated_stage0
        and 'obs_gate0d_transwell_beta_annexin_v_pct' in validated_stage0
        and not indirect_apc_cd4_failed
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
    # Mandatory original exact-product primary biological & pouch endpoints required for READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE (R7S3-F01; DEC-058):
    # Prevents composite COMP-R4 combined-stress secondary coupons or omitted/UNAVAILABLE primary pouch oxygen/peroxide/pH/HLA-II/yield/CD47/SERPINB9/iCasp9 observations from awarding COMPLETE Stage-0 readiness.
    missing_original_primary_endpoints = [
        mk for mk in STAGE0_ORIGINAL_PRIMARY_MANDATORY_MEAN_KEYS if mk not in validated_stage0
    ]
    has_explicit_original_primary_obs = bool(
        not missing_original_primary_endpoints
        and all(
            se_provenance_by_mean.get(mk) == 'EXPLICIT_SE_SUPPLIED'
            for mk in STAGE0_ORIGINAL_PRIMARY_MANDATORY_MEAN_KEYS
        )
        and validated_stage0['obs_stage7_yield_ratio'] >= 0.70
        and _lower95('obs_stage7_yield_ratio') >= 0.70
        and validated_stage0['obs_cd47_fold'] >= 3.0
        and _lower95('obs_cd47_fold') >= 3.0
        and validated_stage0['obs_serpinb9_log10_kass'] >= 6.0
        and _lower95('obs_serpinb9_log10_kass') >= 6.0
        and validated_stage0['obs_icasp9_log10_kill'] >= 4.0
        and _lower95('obs_icasp9_log10_kill') >= 4.0
        and validated_stage0['obs_pouch_min_po2_mmhg'] >= 15.0
        and _lower95('obs_pouch_min_po2_mmhg') >= 15.0
        and validated_stage0['obs_pouch_max_h2o2_uM'] <= 5.0
        and _upper95('obs_pouch_max_h2o2_uM') <= 5.0
        and validated_stage0['obs_pouch_max_ph'] <= 7.45
        and _upper95('obs_pouch_max_ph') <= 7.45
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
        or s10_nkg2c_triggered
        or combined_stress_promoted
        or track2_kir_or_hlag1_fallback_promoted
    )
    any_resolved_component_on_hold = any(
        str(resolved_arch.get(k, '')).startswith('HOLD_')
        for k in ('locus1_construct', 'locus2_construct', 'locus3_construct', 'locus4_construct', 'track2_kir_subbank_construct')
    )
    if safety_kill_s9_triggered:
        lot_release_readiness_tier = 'HARD_SAFETY_STOP_REJECTED'
    elif str(resolved_arch.get('lot_maturation_status', '')).startswith('NO_GO_'):
        lot_release_readiness_tier = 'NO_GO_TERMINAL_BIOTRANSPORT_OR_IMMUNE_CEILING_EXCEEDED'
    elif resolved_arch.get('lot_maturation_status') not in ('READY_FOR_ENCAPSULATION', 'READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY') or any_resolved_component_on_hold:
        lot_release_readiness_tier = 'HOLD_FAILED_OR_INCOMPLETE_GATE'
    elif is_comp_r4_product and combined_stress_coupon_status == 'RETAIN_EMBODIMENT_A_STAR_4143BP':
        lot_release_readiness_tier = 'CANDIDATE_NOT_INDICATED_RETAIN_EMBODIMENT_A_STAR_4143BP'
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
                    else 'RETAIN_4DIGIT_MATCHED_HLAC'
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
            'CLEARED_SCT_SELECTIVITY_COUPON_KC_KA_GE_4'
            if (s10_nkg2c_triggered and s10_coupon_cleared)
            else (
                'STRATIFY_CMV_NKG2C_RECIPIENTS_AND_EVALUATE_SCT_SELECTIVITY_MUTANT_COUPON'
                if s10_nkg2c_triggered
                else 'PASS_CMV_NKG2C_SELECTIVITY_GATE_0A'
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
    }
