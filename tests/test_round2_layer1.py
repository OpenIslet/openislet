"""Round-2 boundary regressions; all non-corpus measurements are synthetic fixtures."""
from __future__ import annotations
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
from unittest import mock
import pytest
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import extract_ring1_corpus as ext
from scripts import rank_islet1_stacks as rank


def semantic_fixture():
    source={'pmid':'12345678','title':'Synthetic B2M knockout unit fixture',
            'abstract':'B2M knockout improved cell viability to 80% compared with 20% in unedited control cells (n = 6).',
            'fulltext_data_sections':{}}
    return ext.extract_one_dossier(source),source


def hazard_fixture():
    return {'study_metadata':{'pmid':'12345678','source_text_type':'Abstract_Only','max_biological_fidelity_tier':'Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse'},
            'arms':[{'gate2_3_4_genetic_and_epitope_edits':[{'gene_symbol_hgnc':'B2M','edit_modality':'CRISPR_Cas9_KO'}]}],
            'atomic_effects':[{'effect_id':'unit1','focal_factor_X':'B2M knockout','verbatim_evidence_quote':'B2M knockout reduced specific lysis to 10% versus 40% in controls while viability was maintained.',
             'outcome_variable_Y':'specific lysis','effect_direction':'Protective_or_Improved','secretome_or_differentiation_penalty_observed':False,
             'log_response_ratio_or_smd':math.log(.25),'endpoint_polarity_registry':{'polarity_class':'Lower_Is_Better_Hazard_Endpoint','orientation_sign':-1}}]}


def public_fixture(source):
    public=copy.deepcopy(source);effect=public['atomic_effects'][0]
    digest=hashlib.sha256(effect['verbatim_evidence_quote'].encode()).hexdigest()
    effect['verbatim_evidence_quote']=f'[REDACTED_NON_OA_EXCERPT: see PMID:12345678 field=verbatim_evidence_quote genes=B2M pert=KO_LOF quote_sha256={digest}]'
    return public


def utility_args(profiles=None):
    return ['ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO','AUTO_NONE','CYTO_CXCL10_SWAP_PDL1_SERPINB9','ER_NONE',
            next(iter(rank.NICHE_VASCULAR_MODULES)),'SAF_AAVS1_CAG_ICASP9_BIALLELIC',profiles if profiles is not None else {'INS':{'n_empirical_studies':1}}]


@pytest.mark.parametrize('field',['log_response_ratio_or_smd','treatment_mean','control_mean','treatment_n','p_value_bound'])
@pytest.mark.parametrize('value',[float('nan'),float('inf'),-float('inf'),True])
def test_semantic_numerical_evidence_rejects_nonfinite_or_boolean(field,value):
    study,source=semantic_fixture();study['atomic_effects'][0][field]=value
    with pytest.raises((ValueError,TypeError)):ext.validate_study_semantics(study,source)


@pytest.mark.parametrize('orientation',[True,False,1.0,-1.0,'1'])
def test_semantic_endpoint_orientation_requires_exact_integer(orientation):
    study,source=semantic_fixture();study['atomic_effects'][0]['endpoint_polarity_registry']['orientation_sign']=orientation
    with pytest.raises(ValueError):ext.validate_study_semantics(study,source)


@pytest.mark.parametrize('collection',['arms','atomic_effects'])
def test_semantic_duplicate_identifiers_rejected(collection):
    study,source=semantic_fixture();study[collection].append(copy.deepcopy(study[collection][0]))
    with pytest.raises(ValueError,match='unique'):ext.validate_study_semantics(study,source)


def test_semantic_span_bounds_cannot_be_satisfied_by_python_slice_clamping():
    study,source=semantic_fixture();ext.validate_study_semantics(study,source)
    study['atomic_effects'][0]['source_span']='abstract[0:999999]|norm=exact'
    with pytest.raises(ValueError,match='offsets'):ext.validate_study_semantics(study,source)
    study['atomic_effects'][0]['source_span']='fulltext_data_sections.results_paragraphs[999][0:5]|norm=exact'
    with pytest.raises(ValueError,match='section index'):ext.validate_study_semantics(study,source)


def test_private_and_public_numeric_hazard_use_same_typed_endpoint():
    private=hazard_fixture();public=public_fixture(private)
    # Mentioning maintained viability elsewhere in the sentence does not change
    # the explicitly measured lysis endpoint's welfare direction.
    a=rank.build_target_evidence_profiles([private],{})['B2M']
    b=rank.build_target_evidence_profiles([public],{})['B2M']
    assert a==b
    assert a['posterior_mean_log_rr']==pytest.approx(math.log(4),abs=5e-5)
    assert a['n_deleterious_or_penalty_effects']==0
    inverted=copy.deepcopy(private);inverted['atomic_effects'][0]['log_response_ratio_or_smd']*=-1
    c=rank.build_target_evidence_profiles([inverted],{})['B2M']
    assert c['posterior_mean_log_rr']==-a['posterior_mean_log_rr']


@pytest.mark.parametrize('change',['missing_digest','wrong_pmid','wrong_field','unknown_gene','trailing_prose'])
def test_ranker_authenticates_opaque_metadata_before_attribution(change):
    study=public_fixture(hazard_fixture());effect=study['atomic_effects'][0];q=effect['verbatim_evidence_quote']
    if change=='missing_digest':q=q[:q.index(' quote_sha256=')]+']'
    elif change=='wrong_pmid':q=q.replace('12345678','87654321')
    elif change=='wrong_field':q=q.replace('field=verbatim_evidence_quote','field=note')
    elif change=='unknown_gene':q=q.replace('genes=B2M','genes=SOURCE_PROSE')
    else:q+=' Additional source prose.'
    effect['verbatim_evidence_quote']=q
    with pytest.raises(ValueError,match='opaque'):rank.build_target_evidence_profiles([study],{})


def test_duplicate_study_and_unknown_source_cannot_gain_empirical_weight():
    study=hazard_fixture()
    with pytest.raises(ValueError,match='Duplicate study'):rank.build_target_evidence_profiles([study,copy.deepcopy(study)],{})
    study['study_metadata']['source_text_type']='UNREVIEWED_UNKNOWN'
    with pytest.raises(ValueError,match='source_text_type'):rank.build_target_evidence_profiles([study],{})


@pytest.mark.parametrize('value',[float('nan'),float('inf'),-float('inf'),True])
def test_profile_numeric_lrr_and_omics_reject_nonfinite_boolean(value):
    study=hazard_fixture();study['atomic_effects'][0]['log_response_ratio_or_smd']=value
    with pytest.raises(ValueError):rank.build_target_evidence_profiles([study],{})
    with pytest.raises(ValueError):rank.build_target_evidence_profiles([],{"B2M":{"gnomad_pli":value}})


@pytest.mark.parametrize('value',[float('nan'),float('inf'),-float('inf'),True,'1'])
def test_nonfinite_or_coerced_signal_cannot_saturate_to_positive_score(value):
    with pytest.raises(ValueError):rank.sigmoid(value)
    with pytest.raises(ValueError):rank.signed_gene_signal({'B2M':{'net_evidence_and_constraint_score':value}},['B2M'])


@pytest.mark.parametrize('value',[True,False,1.5,-1,'1'])
def test_empirical_count_requires_nonnegative_integer(value):
    with pytest.raises(ValueError):rank.compute_series_system_utility(*utility_args({'INS':{'n_empirical_studies':value}}))


@pytest.mark.parametrize('value',[float('nan'),float('inf'),True,-.1])
def test_heuristic_dispersion_must_be_finite_and_nonnegative(value):
    with pytest.raises(ValueError):rank.compute_series_system_utility(*utility_args({'B2M':{'n_empirical_studies':1,'posterior_se':value}}))


def test_unknown_recipient_context_is_not_silently_adult():
    with pytest.raises(ValueError,match='recipient context'):rank.compute_series_system_utility(*utility_args(),recipient_context='Pediatric_Early_Onset_T1D_TYPO')
    assert rank.compute_series_system_utility(*utility_args(),recipient_context='Adult_Onset_LADA')['series_system_joint_survival_probability']>0


def test_bottleneck_is_selected_before_display_rounding():
    with mock.patch.object(rank,'sigmoid',side_effect=[.60004,.60001,.8,.9,.85,.75]):
        result=rank.compute_series_system_utility(*utility_args())
    assert result['limiting_bottleneck_gate']=='Gate2A_Allo_Adaptive_T_and_DSA_Stealth'
    assert result['limiting_bottleneck_score']==.6


def test_finite_sigmoid_and_missing_evidence_positive_controls():
    assert rank.sigmoid(0)==.5
    assert rank.sigmoid(1000)==pytest.approx(1/(1+math.exp(-15)))
    assert rank.run_combinatorial_multi_gate_search(rank.build_target_evidence_profiles([],{}))[0]['composite_series_system_utility']==0


def test_warning_metadata_describes_executable_discount_not_phantom_ceiling(monkeypatch):
    semantics=json.loads(rank.GRAVEYARD_PATH.read_text())['severity_level_semantics']['WARN']
    assert 'Gate 2B indirect-APC ceiling (0.84)' not in semantics
    assert '0.84 ceiling description is not an implemented hard probability ceiling' in semantics
    assert '-0.25' in semantics
    # A source-registry mutation demonstrates the executable model has no0.84
    # clamp; this is a synthetic sensitivity control, not a benefit prediction.
    entry=rank.COMPLETE_GATE_COEFFICIENT_REGISTRY['Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth']['base_intercept_pediatric']
    monkeypatch.setitem(entry,'logit_value',100.)
    result=rank.compute_series_system_utility(*utility_args())
    assert result['gate_scores']['Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth']>.84
