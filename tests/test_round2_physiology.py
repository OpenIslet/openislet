"""Round-2 numerical-boundary controls; all fixtures are synthetic, not experiments."""
from __future__ import annotations
import copy
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import run_layer4_tissue_and_metabolic_twin as l4
from scripts import run_layer5_personalized_patient_compiler as l5
from scripts import ingest_wetlab_results as w


def require(condition, message):
    if not condition:raise AssertionError(message)


def patient(**updates):
    row=dict(patient_id='SYNTHETIC-ROUND2',ancestry='EUR',nkg2a_neg_kir_educated_nk_fraction=.2,
        pre_existing_gad65_ia2_znt8_t_cell_zscore=1.,ifn1_ish_genetic_risk_percentile=40.,sirpa_genotype='v1/v1')
    row.update(updates);return row


@pytest.mark.parametrize('rhs',[lambda t,y:1.,lambda t,y:np.array([1.]),lambda t,y:np.array([[1.,1.]]),
    lambda t,y:np.array([True,False]),lambda t,y:np.array([1.+2j,1.]),lambda t,y:np.array(['1','2'])])
def test_rk4_derivative_cannot_broadcast_or_coerce(rhs):
    with pytest.raises(ValueError,match='RHS'):
        l4.rk4_integrate(rhs,(0.,1.),np.array([1.,2.]),[1.])


@pytest.mark.parametrize('times',[[True],['1'],[np.bool_(True)],[1.+0j]])
def test_rk4_output_times_are_real_numeric_not_coerced(times):
    with pytest.raises(ValueError):l4.rk4_integrate(lambda t,y:-y,(0.,1.),np.array([1.]),times)


@pytest.mark.parametrize('state',[[True],['1'],np.array([1.+2j])])
def test_rk4_initial_state_rejects_nonreal_or_boolean_values(state):
    with pytest.raises(ValueError):l4.rk4_integrate(lambda t,y:-y,(0.,1.),state,[1.])


def test_rk4_vector_positive_control_matches_independent_exponential():
    actual=l4.rk4_integrate(lambda t,y:np.array([-y[0],-2.*y[1]]),(0.,1.),np.array([1.,2.]),[0.,.5,1.],dt_max=.01)
    expected=np.array([[1.,math.exp(-.5),math.exp(-1.)],[2.,2.*math.exp(-1.),2.*math.exp(-2.)]])
    require(np.max(np.abs(actual-expected))<1e-8,'vector RK4 changed the declared ODE')


@pytest.mark.parametrize('d',[True,'1',float('nan'),float('inf'),0.,-1.])
def test_diffusivity_override_has_strict_physical_type_and_domain(d):
    with pytest.raises(ValueError):l4._solve_two_zone_fv_wigner_seitz_bvp(20.,.006,n_per_zone=5,d_hydrogel_override=d)


@pytest.mark.parametrize('dose',[True,False,float('nan'),float('inf'),-1.,0.,49999.,2000001.,'500000'])
def test_standalone_evaluator_rejects_invalid_compiled_doses(dose):
    p=patient();reg=l5.compile_personalized_patient_blueprint(p)
    reg['host_modifier_genotypes_and_dosing']['prescribed_curative_dose_ieq']=dose
    with pytest.raises(ValueError):l5.evaluate_regimen_survival_probability(p,reg)


@pytest.mark.parametrize('dose',[True,float('nan'),float('inf'),-1.])
def test_standalone_evaluator_also_validates_patient_dose_even_when_compiled_value_exists(dose):
    p=patient();reg=l5.compile_personalized_patient_blueprint(p);p['prescribed_dose_ieq']=dose
    with pytest.raises(ValueError):l5.evaluate_regimen_survival_probability(p,reg)


@pytest.mark.parametrize('key',['donor_tyk2_rs34536443_genotype','recipient_tyk2_rs34536443_genotype'])
def test_standalone_evaluator_rejects_invalid_tyk2_genotypes(key):
    p=patient();reg=l5.compile_personalized_patient_blueprint(p);p[key]='INVALID'
    with pytest.raises(ValueError):l5.evaluate_regimen_survival_probability(p,reg)


@pytest.mark.parametrize('location,key',[(None,'final_graft_tyk2_p1104a_protective_state'),
                                      ('host_modifier_genotypes_and_dosing','already_homozygous_tyk2_p1104a')])
def test_derived_flags_cannot_forge_natural_tyk2_protection(location,key):
    p=patient();reg=l5.compile_personalized_patient_blueprint(p)
    block=reg if location is None else reg[location];block[key]=True
    p.update(ifn1_ish_genetic_risk_percentile=90.,pre_existing_gad65_ia2_znt8_t_cell_zscore=3.)
    with pytest.raises(ValueError,match='Inconsistent derived TYK2'):
        l5.evaluate_regimen_survival_probability(p,reg)


@pytest.mark.parametrize('graft', ['G/G','G/C','C/C'])
def test_valid_compiler_output_still_evaluates_for_all_starting_tyk2_states(graft):
    p=patient(donor_tyk2_rs34536443_genotype=graft,ifn1_ish_genetic_risk_percentile=90.)
    score=l5.evaluate_regimen_survival_probability(p,l5.compile_personalized_patient_blueprint(p))
    require(math.isfinite(score) and .05<=score<=.988,'valid compiler output failed')


def test_huge_finite_measurement_error_has_finite_no_information_limit():
    result=w.conjugate_normal_update(.9,.1,.95,1e200,min_valid=0.,max_valid=1.)
    require(result['posterior_mean_exact']==.9 and result['posterior_sd_exact']==.1,'uninformative finite observation moved prior')
    require(result['precision_gain_ratio']==1. and result['kl_divergence_post_vs_prior_nats']==0.,'huge SE created information')
    json.dumps(result,allow_nan=False)


@pytest.mark.parametrize('sd',[1e200,1e308])
def test_unrepresentable_prior_precision_fails_explicitly(sd):
    with pytest.raises(ValueError,match='representable numeric range'):
        w.conjugate_normal_update(.9,sd,.95,.1,min_valid=0.,max_valid=1.)


def test_regular_conjugate_update_matches_independent_closed_form():
    pm,psd,om,ose=.9,.1,.7,.03
    v=1./(1./psd**2+1./ose**2);m=v*(pm/psd**2+om/ose**2)
    result=w.conjugate_normal_update(pm,psd,om,ose,min_valid=0.,max_valid=1.)
    require(abs(result['posterior_mean_exact']-m)<1e-14 and abs(result['posterior_sd_exact']-math.sqrt(v))<1e-14,'conjugate arithmetic drift')


def test_huge_se_partial_batch_can_chain_with_full_state_integrity():
    first=w.run_calibration(dict(batch_id='SYNTHETIC-LARGE-SE',obs_stage7_yield_ratio=.9,obs_stage7_yield_se=1e200),False)
    second=w.run_calibration(dict(batch_id='SYNTHETIC-NEXT',obs_stage7_yield_ratio=.92,obs_stage7_yield_se=.05),False,first)
    require(second['observed_endpoints_updated_count']==1,'stable update could not chain')
    json.dumps(second,allow_nan=False)


def test_far_future_inosculation_is_an_explicit_unvascularized_control():
    result=l4._run_single_pouch_pde_regime('SYNTHETIC_DELAYED_CONTROL',.45,'central',68.,True,.97,
        inosculation_t50_h=10000.,_assess_numerical_gate=False)
    require(result['passes_all_oxygen_peroxide_and_ph_gates'] is False,'delayed vascularization unexpectedly passes')
    require(result['passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate'] is False,'delayed control passed strict qualification')
    json.dumps(result,allow_nan=False)


def test_nonconverged_favorable_oxygen_cannot_pass_product_qualification(monkeypatch):
    original=l4.solve_spherical_cluster_po2_bvp
    def failed_solver(*args,**kwargs):
        result=original(*args,**kwargs)
        result['spherical_bvp_solver_diagnostics']['converged']=False
        return result
    monkeypatch.setattr(l4,'solve_spherical_cluster_po2_bvp',failed_solver)
    result=l4._run_single_pouch_pde_regime('SYNTHETIC_SOLVER_FAILURE',.45,'central',68.,True,.97,_assess_numerical_gate=False)
    require(result['min_slab_po2_across_14_days_mmhg']>=15.,'fixture does not isolate convergence from oxygen level')
    require(result['required_120um_oxygen_solvers_converged'] is False,'solver failure not propagated')
    require(result['passes_all_oxygen_peroxide_and_ph_gates'] is False,'unconverged upper-bound oxygen passed')
    require(result['passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate'] is False,'unconverged strict oxygen passed')
    require(result['strict_gel_only_numerical_gate_assessment']['status']=='NUMERICAL_SOLVER_UNRESOLVED','failure falsely labeled physical certainty')


@pytest.mark.parametrize('pm,om',[(-1.,1.),(1.,3.)])
def test_resealed_prior_cannot_hide_out_of_domain_history_with_in_range_posterior(pm,om):
    report=w.run_calibration(dict(batch_id='SYNTHETIC-DOMAIN-FIRST',obs_stage7_yield_ratio=.9,obs_stage7_yield_se=.1),False)
    ep=report['conjugate_posterior_updates']['stage7_sc_beta_yield_fraction_arm1b']
    # A complete algebraically correct update under the wrong, wider support.
    ep.update(w.conjugate_normal_update(pm,.1,om,.1,min_valid=-2.,max_valid=4.))
    report['conjugate_posterior_updates_sha256']=hashlib.sha256(json.dumps(report['conjugate_posterior_updates'],sort_keys=True,allow_nan=False).encode()).hexdigest()
    report['calibration_state_sha256']=w._calibration_state_sha256(report)
    with pytest.raises(ValueError,match='physical domain'):
        w.run_calibration(dict(batch_id='SYNTHETIC-DOMAIN-NEXT',obs_stage7_yield_ratio=.95,obs_stage7_yield_se=.05),False,report)


def test_resealed_retained_state_cannot_underflow_prior_precision():
    report=w.run_calibration(dict(batch_id='SYNTHETIC-PRECISION-FIRST',obs_stage7_yield_ratio=.9,obs_stage7_yield_se=.1),False)
    ep=report['conjugate_posterior_updates']['cd47_igv_hour0_fold_over_baseline_arm1b']
    for key in ['prior_sd','prior_sd_exact','posterior_sd','posterior_sd_exact']:ep[key]=1e200
    ep['posterior_precision_tau_exact']=0.
    report['conjugate_posterior_updates_sha256']=hashlib.sha256(json.dumps(report['conjugate_posterior_updates'],sort_keys=True,allow_nan=False).encode()).hexdigest()
    report['calibration_state_sha256']=w._calibration_state_sha256(report)
    with pytest.raises(ValueError,match='representable numeric range'):
        w.run_calibration(dict(batch_id='SYNTHETIC-PRECISION-NEXT',obs_stage7_yield_ratio=.95,obs_stage7_yield_se=.05),False,report)


@pytest.mark.parametrize('values',[[5.,5.,5.],[0.,0.,0.],[5.,6.,7.]])
def test_no_secretion_decline_does_not_become_instant_or_arbitrary_shutoff(values):
    require(l4._trajectory_half_drop_time(np.array([448.,449.,450.]),np.array(values),448.,450.) is None,
            'constant/rising secretion was described as a shutoff')


def test_half_drop_uses_real_window_and_linear_time_interpolation():
    elapsed=l4._trajectory_half_drop_time(np.array([448.,449.,450.]),np.array([10.,8.,6.]),448.,450.)
    require(abs(elapsed-1.)<1e-12,'wrong intervention reference clock')
    irregular=l4._trajectory_half_drop_time(np.array([0.,.1,1.]),np.array([10.,8.,0.]),0.,1.)
    require(abs(irregular-.4375)<1e-12,'metric assumed uniform time spacing')


def test_low_glucose_challenge_metrics_use_the_actual_exercise_schedule(monkeypatch):
    original=l4._trajectory_half_drop_time;windows=[]
    def observed_window(times,secretion,start,end):
        windows.append((start,end));return original(times,secretion,start,end)
    monkeypatch.setattr(l4,'_trajectory_half_drop_time',observed_window)
    result=l4.solve_24h_paracrine_and_whole_body_metabolic_ode()
    require(windows.count((448.,536.))==3,'the three low-glucose interventions used the standard protocol metric window')
    for row in result['graft_architectures_compared']:
        require(row['insulin_secretion_shutoff_measurement_window_min']==[600.,660.],'reference protocol event moved')
        require(row['glucagon_peak_measurement_window_including_30min_recovery_min']==[600.,690.],'reference glucagon window moved')



def model_atlas_fixture(tmp_path,monkeypatch,**values):
    atlas=json.loads((ROOT/'data/extracted/layer2_regulatory_and_perturbation_atlas.json').read_text())
    profiles={row['genotype_id']:row for row in atlas['stage0_to_stage7_perturbation_plane']['perturbation_profiles']}
    for name,updates in values.items():profiles[name].update(updates)
    path=tmp_path/'synthetic-model-atlas.json';path.write_text(json.dumps(atlas))
    monkeypatch.setattr(w,'LAYER2_ATLAS_PATH',path)
    return atlas,path


def test_zero_reference_model_ratio_does_not_block_independent_assay_ingestion(tmp_path,monkeypatch):
    model_atlas_fixture(tmp_path,monkeypatch,WT_UNEDITED_CONTROL={'hour0_basal_cd47_molecules_per_cell':0},
                        SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS={'hour0_basal_cd47_molecules_per_cell':0})
    metadata=w.current_model_reference_metadata()['cd47_igv_hour0_fold_over_baseline_arm1b']
    require(metadata['layer2_or_layer3_model_reference_value'] is None,'undefined model fold became a number')
    require(metadata['model_reference_status']=='UNDEFINED_ZERO_REFERENCE_CD47_DENSITY','undefined ratio lacks explicit status')
    result=w.run_calibration(dict(batch_id='SYNTHETIC-ZERO-MODEL',obs_stage7_yield_ratio=.9,obs_stage7_yield_se=.05),False)
    prior=result['conjugate_posterior_updates']['cd47_igv_hour0_fold_over_baseline_arm1b']
    require(prior['prior_mean_exact']==7.81 and prior['status']=='UNAVAILABLE_PRIOR_RETAINED','independent assay prior was retuned')
    require(result['observed_endpoints_updated_count']==1,'valid unrelated assay was blocked')
    json.dumps(result,allow_nan=False)


def test_zero_yield_reference_also_stays_undefined_without_changing_assay_prior(tmp_path,monkeypatch):
    model_atlas_fixture(tmp_path,monkeypatch,WT_UNEDITED_CONTROL={'stage7_sc_beta_yield_fraction':0})
    result=w.run_calibration(dict(batch_id='SYNTHETIC-ZERO-YIELD-MODEL',obs_cd47_fold=6.,obs_cd47_fold_se=.2),False)
    metadata=result['layer2_and_ranker_reconciliation_note_claude_v80_013_017']['current_model_references']['stage7_sc_beta_yield_fraction_arm1b']
    require(metadata['layer2_or_layer3_model_reference_value'] is None,'undefined yield model ratio fabricated')
    require(metadata['model_reference_status']=='UNDEFINED_ZERO_REFERENCE_STAGE7_YIELD','undefined yield status absent')
    require(result['conjugate_posterior_updates']['stage7_sc_beta_yield_fraction_arm1b']['prior_mean_exact']==.985,'yield target prior changed')


@pytest.mark.parametrize('profile,key,value',[
    ('WT_UNEDITED_CONTROL','hour0_basal_cd47_molecules_per_cell',True),
    ('WT_UNEDITED_CONTROL','hour0_basal_cd47_molecules_per_cell',-1.),
    ('WT_UNEDITED_CONTROL','hour0_basal_cd47_molecules_per_cell',float('nan')),
    ('SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS','hour0_basal_cd47_molecules_per_cell',float('inf')),
    ('SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS','stage7_sc_beta_yield_fraction',True),
    ('WT_UNEDITED_CONTROL','stage7_sc_beta_yield_fraction',1.01),
])
def test_model_reference_inputs_cannot_be_coerced_or_nonfinite(tmp_path,monkeypatch,profile,key,value):
    model_atlas_fixture(tmp_path,monkeypatch,**{profile:{key:value}})
    with pytest.raises(ValueError):w.current_model_reference_metadata()


def test_model_reference_profile_ids_are_unambiguous(tmp_path,monkeypatch):
    atlas,path=model_atlas_fixture(tmp_path,monkeypatch)
    atlas['stage0_to_stage7_perturbation_plane']['perturbation_profiles'].append(copy.deepcopy(atlas['stage0_to_stage7_perturbation_plane']['perturbation_profiles'][0]))
    path.write_text(json.dumps(atlas))
    with pytest.raises(ValueError,match='unique genotype'):w.current_model_reference_metadata()
