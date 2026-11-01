"""Adversarial two-pass controls; these test code integrity, not physical qualification."""
from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import ingest_wetlab_results as w
from scripts import run_layer4_tissue_and_metabolic_twin as l4
from scripts import run_layer5_personalized_patient_compiler as l5


def require(test, message):
    if not test:
        raise AssertionError(message)


def observed(batch="first"):
    return dict(batch_id=batch, obs_stage7_yield_ratio=.91, obs_stage7_yield_se=.04,
                obs_pouch_min_po2_mmhg=23., obs_pouch_min_po2_se=.8)


@pytest.mark.parametrize("remove", ["ingested_batch_payload_sha256s", "batch_payload_sha256",
                                    "optional_pouch_and_36run_screen_updates_v80_070",
                                    "calibration_state_sha256"])
def test_prior_state_cannot_drop_replay_history_or_extended_endpoints(remove):
    prior = w.run_calibration(observed(), False)
    prior.pop(remove)
    with pytest.raises(ValueError):
        w.run_calibration(observed("second"), False, prior)


def test_all_eight_endpoint_states_preserved_across_partial_batches():
    prior = w.run_calibration(observed(), False)
    nxt = w.run_calibration(dict(batch_id="second", obs_cd47_fold=6., obs_cd47_fold_se=.4), False, prior)
    require(len(nxt["optional_pouch_and_36run_screen_updates_v80_070"]) == 4, "extended state absent")
    for name, ep in prior["optional_pouch_and_36run_screen_updates_v80_070"].items():
        require(ep["posterior_mean_exact"] == nxt["optional_pouch_and_36run_screen_updates_v80_070"][name]["posterior_mean_exact"], "missing endpoint reset")
    require(nxt["observed_endpoints_updated_count"] == 1, "missing endpoints treated as observations")


def test_histories_remain_required_even_if_state_digest_recomputed():
    prior = w.run_calibration(observed(), False)
    prior["ingested_batch_payload_sha256s"] = []
    prior["calibration_state_sha256"] = w._calibration_state_sha256(prior)
    with pytest.raises(ValueError):
        w.run_calibration(observed("second"), False, prior)


def test_extended_deletion_rejected_even_with_recomputed_digest():
    prior = w.run_calibration(observed(), False)
    prior["optional_pouch_and_36run_screen_updates_v80_070"] = {}
    prior["calibration_state_sha256"] = w._calibration_state_sha256(prior)
    with pytest.raises(ValueError):
        w.run_calibration(dict(batch_id="second", obs_cd47_fold=6., obs_cd47_fold_se=.4), False, prior)


@pytest.mark.parametrize("dt", [0., -1., float("nan"), float("inf"), True])
def test_rk4_invalid_steps_fail_before_loop(dt):
    with pytest.raises(ValueError):
        l4.rk4_integrate(lambda t,y:-y, (0.,1.), np.array([1.]), [1.], dt_max=dt)


@pytest.mark.parametrize("times", [[.5,.1], [-1.], [2.], [float("nan")]])
def test_rk4_invalid_times_fail_closed(times):
    with pytest.raises(ValueError):
        l4.rk4_integrate(lambda t,y:-y, (0.,1.), np.array([1.]), times)


def test_rk4_nonfinite_derivative_fails_instead_of_emitting_nan():
    with pytest.raises(ValueError):
        l4.rk4_integrate(lambda t,y:np.array([float("nan")]), (0.,1.), np.array([1.]), [1.])


@pytest.mark.parametrize("update", [{"release_kinetics_law":"typo"}, {"cao2_half_life_h":0.},
    {"f_chem_initial":-1.}, {"catalase_present":"false"}, {"iec_max_perfusion_fraction":1.01},
    {"demand_multiplier_extra":float("nan")}, {"boundary_po2_mmhg":float("inf")}])
def test_pde_physical_domains_fail_before_solver(update):
    args=dict(regime_id="negative", f_chem_initial=.45, spatial_mode="central", cao2_half_life_h=68.,
              catalase_present=True, iec_max_perfusion_fraction=.97)
    args.update(update)
    with pytest.raises(ValueError):
        l4._run_single_pouch_pde_regime(**args)


def test_layout_b_mass_porosity_and_cluster_fit_are_not_fabricated():
    led=l4.verify_exact_slab_source_distribution_algebra()["non_overlapping_phase_volume_and_cad_ledger_7cm3_ultracode_v80_18"]
    b=led["layout_option_b_31wt_pct_cao2_confined_central_lamina_microbeads"]
    require(abs(b["central_lamina_void_fraction_from_mass_and_density"]-.1754)<.0001,"mass-derived porosity mismatch")
    require(.70 < b["eps_0p42_inventory_retention_at_same_formulation_and_lamina"] < .71,"changed porosity silently keeps same inventory")
    require(abs(b["outer_band_thickness_um_each"]-115.5)<1e-10,"outer band geometry wrong")
    require(b["diameter_120um_full_outer_band_exclusion_feasible"] is False,"oversize sphere declared feasible")


def test_exact_intact_junction_cannot_be_reported_absent_or_immune():
    result=l5.compute_split_junction_residual_cleavage_metrics("A"*20,"AGG","HiFi-SpCas9", "A"*20+"AGG", "A"*20+"AGG")
    require(result["full_20_to_23nt_target_absent_on_both_junctions"] is False,"intact target mislabeled absent")
    require(result["cleavage_immune_post_hdr_hifi_r691a"] is None,"heuristic mislabeled empirical immunity")
    require(result["max_junction_cfd_cleavage_probability"] is None,"uncalibrated score mislabeled probability")


@pytest.mark.parametrize("fn", [l5.revcomp,l5.compute_sugimoto_1995_rna_dna_dg37_kcal_mol,l5.compute_duplex_dg37_kcal_mol,l5.gc_percent])
@pytest.mark.parametrize("bad", ["", "ACNT", "AG T", 42])
def test_ambiguous_sequences_cannot_receive_thermodynamic_scores(fn,bad):
    with pytest.raises(ValueError):fn(bad)


def test_fabricated_offtarget_receipt_is_not_authenticated_by_hash_shape():
    with pytest.raises(ValueError):
        l5.compute_panel_and_reference_off_target_specificity("test", "HiFi-SpCas9", 0, 0, [], None,
          genomewide_enumeration_receipt=dict(tool="x",version="1",reference_build="GRCh38",hit_loci=[],sha256="a"*64))


def test_autologous_guide_selection_ignores_unrelated_donor_mcb_flag():
    patient=dict(patient_id="test",ancestry="EUR",nkg2a_neg_kir_educated_nk_fraction=.4,
       pre_existing_gad65_ia2_znt8_t_cell_zscore=1.,ifn1_ish_genetic_risk_percentile=40.,sirpa_genotype="v1/v1",
       requested_track="TRACK_3_AUTOLOGOUS_IPSC",carries_rare_snp_in_standard_grna=False,
       donor_mcb_carries_rare_snp_in_standard_grna=True)
    result=l5.compile_personalized_patient_blueprint(patient)
    require(result["edited_cell_line_requires_custom_allele_grna"] is False,"unrelated donor SNP contaminated autologous design")


def test_clipped_ci_and_heuristic_are_not_labeled_exact_alternative_posteriors():
    result=w.conjugate_normal_update(.1,.2,0.,.1,min_valid=0.,max_valid=1.)
    q=result["interval_and_robustness_qualification"]
    require("NOT_A_95PCT_TRUNCATED" in q["physical_clipped_fields"],"clipped interval mislabeled")
    require("NOT_A_STUDENT_T_POSTERIOR" in q["student_t_df5_legacy_field"],"heuristic mislabeled")


def test_hwe_inverse_cdf_uses_cumulative_heterozygote_mass():
    draws=[l5._hwe_alternate_alleles((i+.5)/10000., .5) for i in range(10000)]
    require([draws.count(i) for i in (0,1,2)] == [2500,5000,2500], "HWE genotype masses are not (1-p)^2,2p(1-p),p^2")


def test_coupled_verdict_cannot_bypass_failed_conservative_core():
    row=dict(regime_id="strict_failure_with_favorable_upper_bound", min_120um_cluster_core_po2_mmhg=9.,
        min_slab_po2_across_14_days_mmhg=20., max_interstitial_h2o2_uM=1.,max_interstitial_ph=7.4,
        coupled_dynamic_14d_viable_endocrine_fraction_v80_066=1., passes_all_oxygen_peroxide_and_ph_gates=True,
        passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate=False,
        wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001=dict(two_zone_fv_gel_only_supply_120um_core_po2_unrounded_mmhg=1.99))
    result=l4.solve_24h_paracrine_and_whole_body_metabolic_ode([row])["matched_pde_ode_qualification_scenarios"][0]
    require(result["coupled_complete_product_pass"] is False,"uniform-supply PASS bypassed conservative qualification")
    require(result["core_coupling_source"] == "gel_only_unrounded","coupling used favorable upper-bound core")
