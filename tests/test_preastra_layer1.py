"""Category-A Layer 1 negative controls. All assay fixtures here are synthetic tests."""
from __future__ import annotations
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from unittest import mock
import pytest
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from scripts import extract_ring1_corpus as ext
from scripts import rank_islet1_stacks as rank
from scripts import build_study_dossiers as dossiers
from scripts import ingest_wetlab_results as ingestor


def load_acceptance():
    spec = importlib.util.spec_from_file_location('preastra_acceptance', ROOT/'tests/test_adversarial_acceptance.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def empirical_fixture(**observations):
    return ingestor.run_calibration({'batch_id':'SYNTHETIC-UNIT-TEST-NOT-WETLAB', **{spec[3]:'UNAVAILABLE' for spec in ingestor.ENDPOINT_SPECS}, **observations}, is_dry_run=False)


def test_a01_current_enum_consumes_validated_snapshot_and_unmapped_stays_unapplied():
    # Deliberately exercise empirical-mode plumbing; no physical observation is asserted.
    report = empirical_fixture(obs_stage7_yield_ratio=.8, obs_stage7_yield_se=.1)
    status = rank.prepare_wetlab_evaluation_context(report)
    assert status['applied_to_production_priors']
    assert set(status['mapped_nonzero_parameter_deltas']) == {'stage7_trihormonal_120um_mvf_cao2_niche_prior'}
    delta = rank._pdelta('stage7_trihormonal_120um_mvf_cao2_niche_prior', .45)
    with mock.patch.object(Path,'exists',side_effect=AssertionError('hot-loop I/O')), mock.patch.object(Path,'stat',side_effect=AssertionError('hot-loop I/O')), mock.patch.object(Path,'read_text',side_effect=AssertionError('hot-loop I/O')):
        assert rank._pdelta('stage7_trihormonal_120um_mvf_cao2_niche_prior', .45) == delta
        rank._rdelta('Gate1_Maturation_GSIS_and_Paracrine','base_intercept', .4)
    for mutation in ({'epistemic_lineage_status':'MIXED_LINEAGE'}, {'calibration_mode':'SYNTHETIC_DRY_RUN_VERIFICATION_DO_NOT_CLAIM_AS_WETLAB_DATA'}, {'conjugate_posterior_updates_sha256':'0'*64}):
        bad = {**report, **mutation}
        assert not rank.prepare_wetlab_evaluation_context(bad)['applied_to_production_priors']
        assert not rank._load_active_wetlab_posterior_deltas()
    h2o2 = empirical_fixture(obs_pouch_max_h2o2_uM=.5, obs_pouch_max_h2o2_se=.1)
    status = rank.prepare_wetlab_evaluation_context(h2o2)
    assert not status['applied_to_production_priors']
    assert status['mapped_nonzero_parameter_deltas'] == {}
    assert status['unmapped_endpoint_deltas']['pouch_max_peroxide_14d_uM'] != 0
    rank.prepare_wetlab_evaluation_context({'calibration_mode':'NONE'})


def test_a01_one_calibration_read_per_evaluation():
    report = empirical_fixture(obs_stage7_yield_ratio=.8, obs_stage7_yield_se=.1)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td)/'posterior.json'; path.write_text(json.dumps(report))
        original = Path.read_text; reads=[]
        def reader(p,*args,**kwargs):
            if p == path: reads.append(p)
            return original(p,*args,**kwargs)
        with mock.patch.object(Path,'read_text',reader):
            rank.prepare_wetlab_evaluation_context(path=path)
            for _ in range(1000): rank._pdelta('stage7_trihormonal_120um_mvf_cao2_niche_prior', .45)
        assert len(reads) == 1
    rank.prepare_wetlab_evaluation_context({'calibration_mode':'NONE'})


@pytest.mark.parametrize('sentence',[
    'There was not a significant increase in apoptosis.',
    'Knockout failed to improve graft survival.',
    'Deletion did not significantly improve graft survival.',
    'Knockout did not prolong graft survival.',
    'No evidence that B2M knockout improved survival was observed.',
    'B2M deletion was not associated with increased apoptosis.',
    'Neither B2M deletion nor CIITA deletion significantly changed survival.',
])
def test_a02_scoped_null_variants(sentence):
    assert ext.classify_sentence_direction_and_penalty(sentence) == ('Neutral_or_NonInferior',False)


def test_a02_production_multi_target_clause_binding():
    q='While downregulation of SLC30A8 had no effect on beta cell survival, loss of UTP23 markedly reduced cell viability.'
    edits=[{'gene_symbol_hgnc':g,'edit_modality':'CRISPR_Cas9_KO','functional_gate':'Gate4'} for g in ['SLC30A8','UTP23']]
    effects,*_ = ext.extract_atomic_effects_and_modifiers('UNIT',{'pmid':'UNIT','title':'Gene perturbation in human islets','abstract':q,'fulltext_data_sections':{}},edits,'hESC','Stage7_Mature_SC_Islet','None_InVitro_Only','Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse')
    bound = {e['focal_factor_X'].split()[0]:e for e in effects}
    assert bound['SLC30A8']['effect_direction'] == 'Neutral_or_NonInferior'
    assert bound['UTP23']['effect_direction'] == 'Deleterious_or_Impaired'
    assert not bound['SLC30A8']['secretome_or_differentiation_penalty_observed']
    assert bound['UTP23']['secretome_or_differentiation_penalty_observed']
    assert all(e['verbatim_evidence_quote'] in q for e in effects)


def test_a03_decimal_exact_roles_and_typed_inequalities():
    assert ext._num_in_quote(.05,'Cell viability increased (P = .05).','p_value')
    assert not ext._num_in_quote(.0504,'Cell viability increased (P = .05).','p_value')
    assert not ext._num_in_quote(9e-7,'Cell viability increased (P = 1e-7).','p_value')
    assert ext.parse_pvalue_and_n('Viability increased (P < .05, n = 6).') == (None,6)
    assert not ext._num_in_quote(99,'Medium catalog 99-603-CV was used.','sample_size_n')
    assert not ext._num_in_quote(17,'Insulin remained stable until week 17.','treatment_mean')
    assert ext._num_in_quote(6,'Measured in n = 6 independent donors.','treatment_n')
    quote='Insulin secretion was 80 ng/mL/min versus 20 ng/mL/min in controls.'
    assert ext._num_in_quote(80,quote,'treatment_mean')
    assert not ext._num_in_quote(20,quote,'treatment_mean')
    assert ext._num_in_quote(20,quote,'control_mean')


def test_a04_source_identifier_is_not_ignored():
    # Use regenerated production extraction, then alter only the source identifier.
    if not (ROOT/'data/dossiers/38661000.json').exists(): pytest.skip('Private source dossiers intentionally excluded from public export')
    dossier = json.loads((ROOT/'data/dossiers/38661000.json').read_text())
    doc = ext.extract_one_dossier(dossier)
    ext.validate_study_semantics(doc,dossier)
    bad = copy.deepcopy(doc); bad['atomic_effects'][0]['quote_source_dossier_pmid'] = '99999999'
    with pytest.raises(ValueError,match='quote_source_dossier_pmid'): ext.validate_study_semantics(bad,dossier)
    bad['atomic_effects'][0]['quote_source_dossier_pmid']='27163171'
    with pytest.raises(ValueError,match='quote_source_dossier_pmid'): ext.validate_study_semantics(bad,dossier)


def test_a05_recipient_sirpa_matching_does_not_become_cd47_donor_gof():
    doc=json.loads((ROOT/'data/extracted/40668892.json').read_text())
    profile=rank.build_target_evidence_profiles([doc],{'CD47':{}})['CD47']
    assert profile['n_gof_oe_effects'] == 0
    assert profile['gof_weighted_protective_score'] == 0
    assert profile['ko_weighted_protective_score'] == 0
    assert profile.get('contextual_untyped_effects')


def test_a06_execution_rows_and_unevaluated_shortcut():
    protocol=json.loads((ROOT/'data/protocols/hip_core_fractional_factorial_v0.6.2.json').read_text())
    contract=protocol['execution_row_contract']
    assert contract['n_run_rows']==36 and len(contract['run_ids'])==36
    assert contract['n_unique_genotypes']==25 and contract['execute_repeat_rows'] is True
    assert 'all 36 coded run rows' in protocol['biological_replication_and_blocking_plan']['randomization_and_blocking']
    alternative=protocol['unevaluated_61_clone_alternative']
    assert alternative['standard_errors'] is None and alternative['minimum_detectable_effects'] is None
    assert list(protocol['column_sums_across_36_runs_by_factor'].values())==[-10,-10,-6,-10,-6,-10]


def test_a07_plausible_cache_mutants_are_rejected():
    acceptance=load_acceptance()
    original={'mc':{'n_monte_carlo_draws':100,'scenario_a_rank_max':20},'variance':{'share':.991}}
    for keys,val in [(('mc','n_monte_carlo_draws'),101),(('mc','scenario_a_rank_max'),100),(('variance','share'),.9911)]:
        mutant=copy.deepcopy(original); mutant[keys[0]][keys[1]]=val
        assert not acceptance.sensitivity_artifact_matches(mutant,original)
    assert acceptance.sensitivity_artifact_matches(original,copy.deepcopy(original))
    summary=rank._rank_precision_summary([1,2,3,4],[17,18,19,20],'unit fixture')
    lo,hi=summary['scenario_a_rank1_frequency_wilson95']; assert lo<.25<hi
    assert summary['scenario_a_rank_quantiles_linear']['0.5']==2.5


def test_a08_failed_forced_fetch_preserves_rich_cache():
    with tempfile.TemporaryDirectory() as td, mock.patch.object(dossiers,'DOSSIER_DIR',Path(td).resolve()), mock.patch.object(dossiers,'fetch_epmc_fulltext_xml',return_value=None):
        path=Path(td)/'12345678.json'
        old={'pmid':'12345678','abstract':'Existing full scientific abstract','has_parsed_fulltext_xml':True,'fulltext_data_sections':{'results_paragraphs':['Original experimental source sentence']}}
        path.write_text(json.dumps(old))
        _,has_ft,_=dossiers.build_one_dossier('12345678',{'pmid':'12345678','pmcid':'PMC12345678','title':'Metadata only'},[],force=True)
        after=json.loads(path.read_text())
        assert after['abstract']==old['abstract'] and after['fulltext_data_sections']==old['fulltext_data_sections'] and has_ft
    from scripts import identify_and_stratify_corpus as stratifier
    assert stratifier.main.__defaults__ == (False,)


def test_a09_artifact_ledger_facts_are_computed():
    if not (ROOT/'data/extracted/islet1_candidate_stack_rankings.json').exists(): pytest.skip('Precomputed ranking cache intentionally excluded from public export')
    artifact=json.loads((ROOT/'data/extracted/islet1_candidate_stack_rankings.json').read_text())
    ledger=artifact['layer1_ultracode_audit_reconciliation_ledger']
    studies=[json.loads(p.read_text()) for p in (ROOT/'data/extracted').glob('*.json') if p.stem.isdigit()]
    assert ledger['source_counts']['studies']==len(studies)
    assert ledger['source_counts']['atomic_effects']==sum(len(study['atomic_effects']) for study in studies)
    assert ledger['search']['all_stacks']==len(artifact['all_evaluated_combinatorial_stacks'])
    assert ledger['profile_counts']['independently_adjudicated_direction_effects']==0
    assert ledger['registry']['auxiliary_scalar_leaves']==len(ledger['registry']['auxiliary_scalar_paths'])


def test_a10_sintov_direct_c10g_clause_has_no_survival_mean():
    if not (ROOT/'data/curated/36055241.json').exists(): pytest.skip('Private curated source intentionally excluded from public export')
    doc=json.loads((ROOT/'data/curated/36055241.json').read_text()); effect=doc['atomic_effects'][0]
    assert 'C10G SC-islet graft insulin levels remained stable' in effect['verbatim_evidence_quote']
    assert effect['treatment_mean'] is None and effect['log_response_ratio_or_smd'] is None
    dossier=json.loads((ROOT/'data/dossiers/36055241.json').read_text())
    ext._populate_effect_provenance(effect,dossier,ext._compute_dossier_sha256('36055241',dossier),'InVitro_Only')
    assert effect['contrast_metadata']['contrast_type']=='right_censored_followup_observation'
    assert '17-week' in effect['contrast_metadata']['timepoint']
    assert effect['verbatim_evidence_quote'] in dossier['fulltext_data_sections']['results_paragraphs'][59]


@pytest.mark.parametrize('alias',['CD274','PDL1','PD-L1'])
def test_a11_alias_readouts_excluded_actual_oe_retained(alias):
    readout=f'IRF1 silencing led to a significant decrease in IFNalpha-induced {alias} mRNA expression.'
    assert 'CD274' not in {e['gene_symbol_hgnc'] for e in ext.detect_genetic_edits(readout)}
    positive=f'{alias} overexpression protected edited human beta cells from T cell killing.'
    assert 'CD274' in {e['gene_symbol_hgnc'] for e in ext.detect_genetic_edits(positive)}


def test_a05_genuine_donor_gof_and_ko_channels_survive():
    for modality,channel in [("Lentiviral_or_Transposon_OE","gof"),("CRISPR_Cas9_KO","ko")]:
        doc={"study_metadata":{"pmid":"UNIT-POSITIVE-"+channel,"source_text_type":"FullText_XML_or_PDF","max_biological_fidelity_tier":"Tier_2B_Human_Cells_InVitro_HumanImmune_CoCulture"},
             "arms":[{"gate2_3_4_genetic_and_epitope_edits":[{"gene_symbol_hgnc":"CD47","edit_modality":modality}]}],
             "atomic_effects":[{"effect_id":"UNIT-POSITIVE","focal_factor_X":"Donor CD47 "+("overexpression" if channel=="gof" else "knockout"),"verbatim_evidence_quote":"Edited donor cells showed improved cell survival.","effect_direction":"Protective_or_Improved","secretome_or_differentiation_penalty_observed":False,"log_response_ratio_or_smd":None}]}
        profile=rank.build_target_evidence_profiles([doc],{})['CD47']
        assert profile[channel+'_weighted_protective_score']>0
        assert profile[('ko' if channel=='gof' else 'gof')+'_weighted_protective_score']==0


@pytest.mark.parametrize('alias',['CD274','PDL1','PD-L1'])
def test_a11_b2m_readout_subject_and_true_alias_ko(alias):
    assert [(e['gene_symbol_hgnc'],e['edit_modality']) for e in ext.detect_genetic_edits(f'B2M knockout reduced {alias} expression.')] == [('B2M','CRISPR_Cas9_KO')]
    assert ('CD274','CRISPR_Cas9_KO') in [(e['gene_symbol_hgnc'],e['edit_modality']) for e in ext.detect_genetic_edits(f'{alias} knockout impaired cell viability.')]
    assert not {'CD274','HLA-E'} & {e['gene_symbol_hgnc'] for e in ext.detect_genetic_edits('IRF1 silencing led to a significant decrease in IFNalpha-induced PDL1 and HLA-E mRNA expression.')}


def test_a03_real_catalog_number_cannot_be_sample_size():
    if not (ROOT/'data/dossiers/27163171.json').exists(): pytest.skip('Private source dossier intentionally excluded')
    dossier=json.loads((ROOT/'data/dossiers/27163171.json').read_text())
    doc=ext.extract_one_dossier(dossier)
    doc['arms'][0]['sample_size_n']=99
    with pytest.raises(ValueError,match='sample_size_n'): ext.validate_study_semantics(doc,dossier)


def test_a07_throwing_recompute_is_exercised_not_replaced_by_cached_json():
    acceptance=load_acceptance()
    with mock.patch.object(acceptance.rank,'compute_direct_output_logit_sensitivity_analysis',side_effect=RuntimeError('RECOMPUTE-MUST-EXECUTE')) as probe:
        with pytest.raises(RuntimeError,match='RECOMPUTE-MUST-EXECUTE'):
            acceptance.replay_sensitivity_artifact({},[],[],{}, {'plausible_cached_result':True})
        assert probe.call_count == 1


def test_a08_successful_reacquisition_updates_source_hash():
    import hashlib
    source_xml = b'<article><front><article-meta><article-id pub-id-type="pmid">12345678</article-id><article-id pub-id-type="pmc">PMC12345678</article-id></article-meta></front></article>'
    with tempfile.TemporaryDirectory() as td, mock.patch.object(dossiers,'DOSSIER_DIR',Path(td).resolve()), mock.patch.object(dossiers,'fetch_epmc_fulltext_xml',return_value=source_xml), mock.patch.object(dossiers,'parse_jats_fulltext',return_value={'results_paragraphs':['New measured source text, fixture only']}):
        path=Path(td)/'12345678.json'; path.write_text(json.dumps({'pmid':'12345678','abstract':'Old abstract','fulltext_data_sections':{'results_paragraphs':['Old source text']}}))
        before=hashlib.sha256(path.read_bytes()).hexdigest()
        dossiers.build_one_dossier('12345678',{'pmid':'12345678','pmcid':'PMC12345678','abstract':'New source abstract'},[],force=True)
        after=json.loads(path.read_text())
        assert after['abstract']=='New source abstract' and after['fulltext_data_sections']['results_paragraphs']==['New measured source text, fixture only']
        assert before != hashlib.sha256(path.read_bytes()).hexdigest()


def test_a01_only_consumed_endpoint_changes_search_and_snapshot_is_stable():
    _, studies, omics=rank.load_inputs()
    profiles=rank.build_target_evidence_profiles(studies,omics)
    rank.prepare_wetlab_evaluation_context({'calibration_mode':'NONE'})
    baseline=rank.run_combinatorial_multi_gate_search(profiles)
    report=empirical_fixture(obs_stage7_yield_ratio=.8,obs_stage7_yield_se=.1)
    rank.prepare_wetlab_evaluation_context(report)
    shifted=rank.run_combinatorial_multi_gate_search(profiles)
    assert shifted != baseline
    assert shifted == rank.run_combinatorial_multi_gate_search(profiles)
    rank.prepare_wetlab_evaluation_context(empirical_fixture(obs_pouch_max_h2o2_uM=.5,obs_pouch_max_h2o2_se=.1))
    assert rank.run_combinatorial_multi_gate_search(profiles) == baseline
    import hashlib
    tampered=copy.deepcopy(report)
    tampered['conjugate_posterior_updates']['stage7_sc_beta_yield_fraction_arm1b']['posterior_mean_exact'] += .1
    tampered['conjugate_posterior_updates_sha256']=hashlib.sha256(json.dumps(tampered['conjugate_posterior_updates'],sort_keys=True,allow_nan=False).encode()).hexdigest()
    assert not rank.prepare_wetlab_evaluation_context(tampered)['applied_to_production_priors']


def test_a03_inequality_metadata_contains_bound_not_observed_p():
    quote='Edited cells maintained viability (P < .05, n = 6 donors).'
    dossier={'pmid':'UNIT','abstract':quote}
    effect={'verbatim_evidence_quote':quote,'effect_id':'UNIT','focal_factor_X':'Edited cells','outcome_domain':'Cell_Viability_or_Apoptosis','conditioning_context_Z':{}}
    ext._populate_effect_provenance(effect,dossier,ext._compute_dossier_sha256('UNIT',dossier),'InVitro_Only')
    assert effect['p_value'] is None and effect['p_value_relation']=='<' and effect['p_value_bound']==.05


def test_a02_graph_penalty_metadata_does_not_claim_independent_adjudication():
    graph=json.loads((ROOT/'data/extracted/cross_study_discovery_graph.json').read_text())
    diagnostics=graph['corpus_epistemic_and_fidelity_diagnostics']
    assert diagnostics['heuristic_penalty_effects_count']==diagnostics['adjudicated_penalty_effects_count']
    assert diagnostics['deprecated_field_aliases']['adjudicated_penalty_effects_count']['alias_of']=='heuristic_penalty_effects_count'
    assert 'no independent manual accuracy dataset' in diagnostics['adjudication_semantics']
    assert sum(diagnostics['direction_and_penalty_method_counts'].values())==diagnostics['active_deduplicated_effects_count']


def test_producer_jsonl_exports_are_current_and_match_manifest_format():
    from scripts.study_artifact_exports import write_study_exports
    study={'study_metadata':{'pmid':'12345679','citation_short':'Test fixture'},'atomic_effects':[{'effect_id':'one','focal_factor_X':'B2M'}]}
    with tempfile.TemporaryDirectory() as td:
        root=Path(td).resolve();studies_path,effects_path,count=write_study_exports(root,[study])
        assert count==1
        assert studies_path.read_text()==json.dumps(study,ensure_ascii=False)+'\n'
        expected={'pmid':'12345679','citation_short':'Test fixture','source_text_type':None,'superseded_by_pmid':None,**study['atomic_effects'][0]}
        assert effects_path.read_text()==json.dumps(expected,ensure_ascii=False)+'\n'
        study['atomic_effects'][0]['effect_id']='new-current-source'
        _,_,count=write_study_exports(root,[study])
        assert json.loads(effects_path.read_text())['effect_id']=='new-current-source'
        assert '"effect_id": "one"' not in effects_path.read_text()


def test_hu_written_recipient_count_is_source_backed_not_catalog_number():
    source='FLuc+ wtallo and HIPallo islets were transplanted into one allogeneic rhesus macaque each.'
    assert ext._num_in_quote(1,source,'sample_size_n')
    assert not ext._num_in_quote(2,source,'sample_size_n')
    assert not ext._num_in_quote(1,'One well received 1 microgram of Cas9 reagent.','sample_size_n')


def test_lrr_shrinkage_source_fixture_changes_branch_without_modifying_persisted_source():
    acceptance=load_acceptance()
    doc=json.loads((ROOT/'data/extracted/35918471.json').read_text())
    before=copy.deepcopy(doc)
    probe=acceptance.source_bound_lrr_shrinkage_probe([doc],{})
    assert probe['verified'] and probe['source_effect_id']=='pmid35918471_eff1'
    assert probe['persisted_study_modified'] is False and doc==before
    assert probe['before_shrunk_lrr'] != probe['after_shrunk_lrr']

# Comprehensive two-pass audit: boundary and challenger negative controls.
def test_master_graph_schema_rejects_nonfinite_bounds_and_ambiguous_oneof():
    from scripts.build_cross_study_graph import validate_against_schema
    cases = [(float('nan'), {'type':'number'}), (float('inf'), {'type':'number'}),
             (-.1, {'type':'number','minimum':0}), (1.1, {'type':'number','maximum':1}),
             (2, {'oneOf':[{'type':'number'},{'type':'integer'}]})]
    for value, schema in cases:
        with pytest.raises(ValueError): validate_against_schema(value, schema)
    validate_against_schema(.5, {'type':'number','minimum':0,'maximum':1})


def test_master_actual_schema_enforces_pvalue_bound():
    from scripts.build_cross_study_graph import validate_against_schema
    schema = json.loads((ROOT/'schemas/extraction_schema_v0.2.json').read_text())
    node = schema['properties']['atomic_effects']['items']['properties']['p_value_bound']
    for value in [-.1, 1.1, float('nan'), float('inf')]:
        with pytest.raises(ValueError): validate_against_schema(value,node)


def test_master_hip_gate_survives_optimization_and_rejects_coercion():
    protocol=json.loads((ROOT/'data/protocols/hip_core_fractional_factorial_v0.6.2.json').read_text())
    rows=protocol['runs']
    with pytest.raises(ValueError):rank._hip_design_summary([rows[0]]*36)
    with pytest.raises(ValueError):rank._hip_design_summary([])
    for value in [True,1.,'1',0,2]:
        mutant=copy.deepcopy(rows);mutant[0]['coded_levels']['X1']=value
        with pytest.raises(ValueError):rank._hip_model_rows(mutant)
    mutant=copy.deepcopy(rows);mutant[0]['stage_block']=True
    with pytest.raises(ValueError):rank._hip_model_rows(mutant)


def test_master_graveyard_missing_malformed_or_unsupported_rules_fail_closed(tmp_path,monkeypatch):
    original=json.loads(rank.GRAVEYARD_PATH.read_text())
    path=tmp_path/'graveyard.json';monkeypatch.setattr(rank,'GRAVEYARD_PATH',path)
    with pytest.raises(ValueError):rank._load_graveyard_entries_by_id()
    mutations=[[],{}, {'graveyard_entries':[]}, {'graveyard_entries':original['graveyard_entries'][:-1]},
               {'graveyard_entries':original['graveyard_entries']+[original['graveyard_entries'][0]]}]
    for mutate in mutations:
        path.write_text(json.dumps(mutate));rank._GRAVEYARD_CACHE['entries']=None
        with pytest.raises(ValueError):rank._load_graveyard_entries_by_id()
    for key,val in [('severity_level','APPROVED'),('trigger_pattern',{'requires_any_module':['UNKNOWN']})]:
        bad=copy.deepcopy(original);bad['graveyard_entries'][0][key]=val;path.write_text(json.dumps(bad));rank._GRAVEYARD_CACHE['entries']=None
        with pytest.raises(ValueError):rank._load_graveyard_entries_by_id()
    path.write_text(json.dumps(original));rank._GRAVEYARD_CACHE['entries']=None
    assert len(rank._load_graveyard_entries_by_id())==10


def test_master_ins_edit_does_not_clear_indirect_apc_warning_or_understate_edits():
    selected=['ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO','AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT','CYTO_CXCL10_SWAP_PDL1_SERPINB9','ER_NONE','NICHE_STAGE7_TRIHORMONAL_120UM_MVF_CAO2_MICROGELS','SAF_AAVS1_CAG_ICASP9_BIALLELIC']
    # Use registry keys rather than duplicated labels for the niche fixture.
    selected[4]=next(iter(rank.NICHE_VASCULAR_MODULES))
    flags=rank.evaluate_graveyard_rules_for_stack(*selected,1)
    assert any(r['graveyard_id'].startswith('GRAVEYARD_010') for r in flags)
    for value in [0,True,1.,-1]:
        with pytest.raises(ValueError):rank.evaluate_graveyard_rules_for_stack(*selected,value)


def test_master_export_rejects_stale_pin_before_any_output(tmp_path,monkeypatch):
    import hashlib
    from scripts import generate_release_manifest as export
    root=tmp_path/'root';root.mkdir();(root/'README.md').write_text('original')
    manifest={'artifacts':{'README.md':{'sha256':hashlib.sha256(b'original').hexdigest(),'bytes':len(b'original')}},'public_redistribution_allowlist_spec':{'public_package_included_paths':['README.md','manifest.json'],'local_research_cache_excluded_from_public_bundle':[]}}
    (root/'manifest.json').write_text(json.dumps(manifest));(root/'README.md').write_text('modified bytes')
    monkeypatch.setattr(export,'ROOT',root)
    dest=tmp_path/'public';unrelated=tmp_path/'.public.staging';unrelated.mkdir();(unrelated/'preserve').write_text('owned by somebody else')
    with pytest.raises(ValueError,match='Pinned source hash mismatch'):export.export_public_bundle(dest)
    assert not dest.exists() and (unrelated/'preserve').read_text()=='owned by somebody else'
    (root/'README.md').write_text('original');export._verify_pinned_source(root,'README.md',manifest)
    with pytest.raises(ValueError):export._verify_pinned_source(root,'../README.md',manifest)


def test_master_portfolio_separates_grid_and_unmaterialized_molecular_hypothesis():
    if not (ROOT/'data/extracted/islet1_candidate_stack_rankings.json').exists(): pytest.skip('Precomputed ranking cache intentionally excluded from public export')
    document=json.loads((ROOT/'data/extracted/islet1_candidate_stack_rankings.json').read_text())
    for calibration in [False,True]:
        portfolio=rank.build_champion_challenger_portfolio(document['ranked_candidate_architectures'],document['all_evaluated_combinatorial_stacks'],{'applied_to_production_priors':calibration})
        grid=portfolio['slots']['X'];molecular=portfolio['blueprint_x_stage0_molecular_challenger_spec'];promotion=portfolio['two_step_promotion_gate']
        assert grid['B']==1 and grid['molecular_proposal_included_in_score'] is False
        assert molecular['B']==0 and molecular['execution_supported'] is False
        for key in ['physical_manufacturing_yield','exact_donor_length_bp','packaged_length_bp','resting_hlae_surface_density_molecules_per_cell','resting_cd47_surface_density_molecules_per_cell','composite_series_system_utility']:
            assert molecular[key] is None
        assert promotion['promotion_authorized'] is False and promotion['exact_product_empirical_gate_results'] is None


def test_master_omics_annotations_do_not_assert_exact_product_immunity_or_amyloid_cure():
    from scripts.summarize_omics_baseline import ENGINEERING_ANNOTATIONS
    assert '100%' not in ENGINEERING_ANNOTATIONS['IAPP']['rationale']
    assert 'unmeasured exact-product' in ENGINEERING_ANNOTATIONS['IAPP']['rationale']
    assert 'indirect presentation' in ENGINEERING_ANNOTATIONS['CIITA']['rationale']
    assert 'some CMV UL40 isolates' in ENGINEERING_ANNOTATIONS['HLA-E']['rationale']
    assert 'SERPINB9_P2A_PDL1' in ENGINEERING_ANNOTATIONS['CXCL10']['recommended_edit']


def test_master_export_failure_does_not_delete_unowned_staging_sibling(tmp_path,monkeypatch):
    import hashlib
    from scripts import generate_release_manifest as export
    root=tmp_path/'root';root.mkdir();(root/'README.md').write_text('original')
    manifest={'artifacts':{'README.md':{'sha256':hashlib.sha256(b'original').hexdigest(),'bytes':len(b'original')}},'public_redistribution_allowlist_spec':{'public_package_included_paths':['README.md','manifest.json'],'local_research_cache_excluded_from_public_bundle':[]}}
    (root/'manifest.json').write_text(json.dumps(manifest));monkeypatch.setattr(export,'ROOT',root)
    sibling=tmp_path/'.public.staging';sibling.mkdir();(sibling/'owned.txt').write_text('preserve')
    # This intentionally incomplete fixture fails after staging/copying, not
    # before. Only this export's unique temporary directory may be removed.
    with pytest.raises(FileNotFoundError):export.export_public_bundle(tmp_path/'public')
    assert (sibling/'owned.txt').read_text()=='preserve'
    assert not list(tmp_path.glob('.public.staging-*'))


def test_master_portfolio_watchlist_never_substitutes_nonqualifying_candidate():
    if not (ROOT/'data/extracted/islet1_candidate_stack_rankings.json').exists(): pytest.skip('Precomputed ranking cache intentionally excluded from public export')
    doc=json.loads((ROOT/'data/extracted/islet1_candidate_stack_rankings.json').read_text())
    grid=doc['all_evaluated_combinatorial_stacks'];single=[copy.deepcopy(grid[1])]
    single[0]['triggered_graveyard_rules']=[{'severity_level':'WARN','graveyard_id':'fixture'}]
    p=rank.build_champion_challenger_portfolio(doc['ranked_candidate_architectures'],single,{})
    assert p['pareto_challenger_watchlist']['zero_graveyard_warn_challenger']['combination_id'] is None
    p=rank.build_champion_challenger_portfolio(doc['ranked_candidate_architectures'],list(reversed(grid)),{})
    assert p['slots']['X']['global_rank']==2


def test_master_review_metadata_cannot_become_primary_experimental_arm():
    if not (ROOT/'data/dossiers/35482056.json').exists(): pytest.skip('Private dossier omitted in public bundle')
    from scripts.build_cross_study_graph import validate_against_schema
    dossier=json.loads((ROOT/'data/dossiers/35482056.json').read_text())
    doc=ext.extract_one_dossier(dossier)
    assert doc['study_metadata']['source_text_type']=='Review_or_Commentary'
    for field in ['arms','atomic_effects','causal_mechanism_edges','effect_modifier_interactions']:
        assert doc[field]==[]
    assert doc['screening_decision']['is_crispr_or_high_throughput_screen'] is False
    assert doc['screening_decision']['exclusion_reason']=='Review_Editorial_or_Commentary'
    assert doc['screening_decision']['included'] is False
    validate_against_schema(doc,json.loads((ROOT/'schemas/extraction_schema_v0.2.json').read_text()))


def test_master_source_eligibility_crosses_effect_export_boundary(tmp_path):
    from scripts.study_artifact_exports import write_study_exports
    studies=[{'study_metadata':{'pmid':'12345679','source_text_type':'Preprint_Superseded_By_Peer_Reviewed','superseded_by_pmid':'12345678'},'atomic_effects':[{'effect_id':'one'}]}]
    _,effects,_=write_study_exports(tmp_path,studies);row=json.loads(effects.read_text())
    assert row['source_text_type']=='Preprint_Superseded_By_Peer_Reviewed' and row['superseded_by_pmid']=='12345678'



def test_master_short_verbatim_fragments_cannot_evade_aggregate_quote_budget():
    from scripts import generate_release_manifest as export
    source = "first exact copied source phrase followed by second exact copied source phrase"
    fragments = {"arm_label": "first exact copied source phrase", "focal_factor_X": "second exact copied source phrase"}
    expected=sum(map(len,fragments.values()))
    assert export._max_verbatim_shingle_chars(fragments,source)==expected
    assert export._max_verbatim_shingle_chars("unrelated short descriptor",source)==0


@pytest.mark.parametrize("policy",[
    {"non_oa_per_article_aggregate_verbatim_ceiling_characters":241},
    {"non_oa_per_article_aggregate_verbatim_ceiling_characters":True},
    {"non_oa_per_article_max_abstract_share":.26},
    {"non_oa_per_article_max_abstract_share":float('nan')},
])
def test_master_manifest_cannot_weaken_public_quote_policy(tmp_path,monkeypatch,policy):
    from scripts import generate_release_manifest as export
    root=tmp_path/'root';root.mkdir();(root/'manifest.json').write_text(json.dumps({'artifacts':{},'public_redistribution_allowlist_spec':policy}))
    monkeypatch.setattr(export,'ROOT',root)
    with pytest.raises(ValueError):export.export_public_bundle(tmp_path/'public')



def test_master_discovery_reharvest_cannot_publish_full_abstract():
    import hashlib
    from scripts.generate_release_manifest import _public_discovery_record
    private={'pmid':'12345678','title':'Citation title','abstract':'Synthetic source text. '*100}
    public=_public_discovery_record(private)
    assert 'abstract' not in public and public['has_abstract'] is True
    assert public['abstract_sha256']==hashlib.sha256(private['abstract'].encode()).hexdigest()
    assert 'abstract' in private and _public_discovery_record(public)==public


def test_master_stratifier_preserves_two_primary_reports_as_empirical():
    from scripts import identify_and_stratify_corpus as strat
    assert strat.PROJECT_SUMMARY_PMIDS=={'42430558'}
    assert not {'41617975','41993274'} & strat.PROJECT_SUMMARY_PMIDS


def test_master_renderer_refreshes_all_current_count_formats_and_is_idempotent():
    from scripts.render_documentation import derive_facts, render_text
    facts=derive_facts(ROOT, len(json.loads((ROOT/'manifest.json').read_text())['artifacts']))
    text='# Test\n\n- **`3,101` atomic quantitative & qualitative effects** with **`3,101 / 3,101` (`100.0%`) exact (`norm=exact`) single-field substring match**.\n- Corpus (`563` studies, `594` arms, `530` strict peer-reviewed empirical studies).\n'
    rendered=render_text(text,facts,'README.md')
    count=facts['corpus']['atomic_effects'];arms=facts['corpus']['arms'];strict=facts['corpus']['ring1_corpus_count_reconciliation']['strict_active_peer_reviewed_empirical_studies_total']
    assert f'`{count:,} / {count:,}`'in rendered
    assert f'`{arms}` arms'in rendered and f'`{strict}` strict peer-reviewed empirical'in rendered
    assert render_text(rendered,facts,'README.md')==rendered
    historical=render_text('# Prior design\n\nHistorical interpretation.',facts,'ISLET-1-SYNTHESIS-REPORT-v0.4.md')
    assert historical.count('**Earlier-revision rationale:**')==1
    assert render_text(historical,facts,'ISLET-1-SYNTHESIS-REPORT-v0.4.md')==historical


def test_master_quote_counter_and_exporter_share_key_bound_metadata_scope():
    from scripts import generate_release_manifest as export
    source = 'The kidney capsule received cells. None were lost. The graft remained functional.'
    doc = {'arms': [{'transplant_site': 'Kidney_Capsule'}] * 30,
           'atomic_effects': [{'conditioning_context_Z': {'Z6_microenvironment_oxygen_and_vascular_niche': 'Kidney_Capsule'}}] * 30}
    assert export._max_verbatim_shingle_chars(doc, source) == 0
    # A source sentence/code lookalike in free text cannot use the enum exemption.
    for field in ('verbatim_evidence_quote', 'arm_label', 'focal_factor_X', 'outcome_variable_Y'):
        assert not export._is_fixed_public_schema_metadata('None', field)
        assert export._max_verbatim_shingle_chars({field: 'None'}, source) == 4
        assert export._max_verbatim_shingle_chars({field: 'Kidney_Capsule'}, source) == len('Kidney_Capsule')
    # Case changes and prose appended to a typed metadata key are not exact codes.
    for value in ('kidney capsule', 'The graft remained functional'):
        assert not export._is_fixed_public_schema_metadata(value, 'transplant_site')
        assert export._max_verbatim_shingle_chars({'transplant_site': value}, source) == len(value)
    fragments = [{'verbatim_evidence_quote': 'The graft remained functional'}] * 10
    assert export._max_verbatim_shingle_chars(fragments, source) > 240
