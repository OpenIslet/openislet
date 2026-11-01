"""Independent source identity and numeric-domain controls at supporting producers."""
import copy
import json
import os
from pathlib import Path
from unittest.mock import Mock

import pytest

from scripts import build_study_dossiers as dossiers
from scripts import summarize_omics_baseline as omics


@pytest.fixture
def omics_source(tmp_path, monkeypatch):
    root=tmp_path/'source';folder=root/'data/omics';folder.mkdir(parents=True);(root/'data/extracted').mkdir()
    raw={'data':{'gene':{'gene_id':'ENSG00000166710','symbol':'B2M','gnomad_constraint':{'pli':.7,'oe_lof':.2,'oe_lof_lower':.1,'oe_lof_upper':1.5}}}}
    (folder/'gnomad_B2M.json').write_text(json.dumps(raw))
    monkeypatch.setattr(omics,'BASE_DIR',root);monkeypatch.setattr(omics,'OMICS_DIR',folder)
    for p in [folder/'target_prior_table.json',root/'data/extracted/omics_baseline_summary.json']:p.write_text('{"unchanged":true}\n')
    return root,folder,raw


def assert_omics_not_published(root,folder):
    assert (folder/'target_prior_table.json').read_text()=='{"unchanged":true}\n'
    assert (root/'data/extracted/omics_baseline_summary.json').read_text()=='{"unchanged":true}\n'


@pytest.mark.parametrize('field,value', [('pli',True),('pli','0.7'),('pli',float('nan')),('pli',float('inf')),('pli',-.1),('pli',1.1),('pli',10**400),('oe_lof',True),('oe_lof',float('nan')),('oe_lof_upper',float('inf')),('oe_lof_upper',-1.)])
def test_omics_invalid_numeric_input_cannot_publish_prior(omics_source,field,value):
    root,folder,raw=omics_source;raw['data']['gene']['gnomad_constraint'][field]=value;(folder/'gnomad_B2M.json').write_text(json.dumps(raw))
    with pytest.raises(ValueError):omics.main()
    assert_omics_not_published(root,folder)


def test_omics_allows_missing_metrics_and_loeuf_above_one(omics_source):
    _,folder,raw=omics_source;raw['data']['gene']['gnomad_constraint']['pli']=None;(folder/'gnomad_B2M.json').write_text(json.dumps(raw));omics.main()
    record=json.loads((folder/'target_prior_table.json').read_text())['targets'][0]
    assert record['gnomad_pli'] is None and record['gnomad_loeuf_upper']==1.5


@pytest.mark.parametrize('field,value',[('data',False),('gene',False),('gnomad_constraint',False),('gnomad_constraint',[]),('gnomad_constraint',0)])
def test_omics_false_container_is_not_silently_missing_data(omics_source,field,value):
    root,folder,raw=omics_source
    if field=='data':raw['data']=value
    elif field=='gene':raw['data']['gene']=value
    else:raw['data']['gene']['gnomad_constraint']=value
    (folder/'gnomad_B2M.json').write_text(json.dumps(raw))
    with pytest.raises(ValueError):omics.main()
    assert_omics_not_published(root,folder)


@pytest.mark.parametrize('mutation',['gene-identity','interval'])
def test_omics_rejects_inconsistent_source_identity_or_interval(omics_source,mutation):
    root,folder,raw=omics_source
    if mutation=='gene-identity':raw['data']['gene']['symbol']='CD47'
    else:raw['data']['gene']['gnomad_constraint']['oe_lof_lower']=2.
    (folder/'gnomad_B2M.json').write_text(json.dumps(raw))
    with pytest.raises(ValueError):omics.main()
    assert_omics_not_published(root,folder)


def test_omics_preserves_pLI_precedence_over_rounded_legacy_alias(omics_source):
    _,folder,raw=omics_source
    raw['data']['gene']['gnomad_constraint'].update({'pLI':.9894400238990784,'pli':.9894})
    (folder/'gnomad_B2M.json').write_text(json.dumps(raw));omics.main()
    record=json.loads((folder/'target_prior_table.json').read_text())['targets'][0]
    assert record['gnomad_pli']==.9894400238990784


@pytest.mark.parametrize('gene_id',[None,True,123,'','ENST00000166710','ENSGbad'])
def test_omics_gene_identifier_requires_human_gene_syntax(omics_source,gene_id):
    root,folder,raw=omics_source;raw['data']['gene']['gene_id']=gene_id
    (folder/'gnomad_B2M.json').write_text(json.dumps(raw))
    with pytest.raises(ValueError):omics.main()
    assert_omics_not_published(root,folder)


@pytest.mark.parametrize('existing',[False,True])
def test_omics_second_publication_failure_restores_both_outputs(omics_source,monkeypatch,existing):
    root,folder,_=omics_source;prior=folder/'target_prior_table.json';summary=root/'data/extracted/omics_baseline_summary.json'
    if existing:
        prior.chmod(0o640);summary.chmod(0o600)
    else:
        prior.unlink();summary.unlink()
    io_os=omics.publish_byte_artifacts.__globals__['os'];replace=io_os.replace;injected=False
    def fail_summary_once(source,destination,*args,**kwargs):
        nonlocal injected
        if Path(destination)==summary and not injected:
            injected=True;raise OSError('INJECTED_SECOND_OMICS_COMMIT_FAILURE')
        return replace(source,destination,*args,**kwargs)
    monkeypatch.setattr(io_os,'replace',fail_summary_once)
    with pytest.raises(OSError,match='INJECTED_SECOND_OMICS_COMMIT_FAILURE'):omics.main()
    assert injected and not (root/'data/.omics-publication.lock').exists()
    if existing:
        assert_omics_not_published(root,folder)
        assert prior.stat().st_mode&0o777==0o640 and summary.stat().st_mode&0o777==0o600
    else:
        assert not prior.exists() and not summary.exists()


@pytest.mark.parametrize('mutation',['disease','duplicate','boolean','missing-score','negative','above-one'])
def test_opentargets_identity_and_score_are_checked_before_publication(omics_source,mutation):
    root,folder,_=omics_source;row={'target':{'approvedSymbol':'B2M'},'score':.3};data={'disease':{'id':'MONDO_0005147','associatedTargets':{'rows':[row]}}}
    if mutation=='disease':data['disease']['id']='MONDO_0005148'
    elif mutation=='duplicate':data['disease']['associatedTargets']['rows'].append(copy.deepcopy(row))
    elif mutation=='boolean':row['score']=True
    elif mutation=='missing-score':row.pop('score')
    elif mutation=='negative':row['score']=-.1
    else:row['score']=1.1
    (folder/'t1d_associated_targets_top100.json').write_text(json.dumps(data))
    with pytest.raises(ValueError):omics.main()
    assert_omics_not_published(root,folder)


def test_opentargets_valid_probability_is_preserved(omics_source):
    _,folder,_=omics_source;(folder/'t1d_associated_targets_top100.json').write_text(json.dumps({'disease':{'id':'MONDO_0005147','associatedTargets':{'rows':[{'target':{'approvedSymbol':'B2M'},'score':.3}]}}}));omics.main()
    record=json.loads((folder/'target_prior_table.json').read_text())['targets'][0]
    assert record['opentargets_t1d_score']==.3 and record['opentargets_score_source']=='OpenTargets_MONDO_0005147_Direct_GraphQL_Top25'


@pytest.mark.parametrize('pmid',['../escape','/absolute','01','0','1234567890','１２３','123\n',True])
def test_dossier_identifier_cannot_escape_cache_or_alias_another_identity(tmp_path,monkeypatch,pmid):
    folder=tmp_path/'dossiers';monkeypatch.setattr(dossiers,'DOSSIER_DIR',folder);fetch=Mock(side_effect=AssertionError('No network allowed'));monkeypatch.setattr(dossiers,'fetch_epmc_fulltext_xml',fetch)
    with pytest.raises(ValueError):dossiers.build_one_dossier(pmid,{'pmid':pmid},[])
    assert not folder.exists() and not (tmp_path/'escape.json').exists();fetch.assert_not_called()


@pytest.mark.parametrize('force',[False,True])
def test_wrong_cached_pmid_is_never_reused_or_merged(tmp_path,monkeypatch,force):
    monkeypatch.setattr(dossiers,'DOSSIER_DIR',tmp_path);p=tmp_path/'40757665.json';original=json.dumps({'pmid':'38508139','abstract':'Unrelated paper'});p.write_text(original)
    with pytest.raises(ValueError):dossiers.build_one_dossier('40757665',{'pmid':'40757665'},[],force)
    assert p.read_text()==original


@pytest.mark.parametrize('record',[{'pmid':'38508139'},{'pmid':'40757665','pmcid':'PMC-UNIT'}, {'pmid':'40757665','pmcid':'../other'}])
def test_mismatched_or_noncanonical_source_metadata_fails_before_fetch(tmp_path,monkeypatch,record):
    monkeypatch.setattr(dossiers,'DOSSIER_DIR',tmp_path);fetch=Mock();monkeypatch.setattr(dossiers,'fetch_epmc_fulltext_xml',fetch)
    with pytest.raises(ValueError):dossiers.build_one_dossier('40757665',record,[])
    fetch.assert_not_called();assert not list(tmp_path.iterdir())


def test_dossier_ordinary_cached_identity_is_preserved(tmp_path,monkeypatch):
    monkeypatch.setattr(dossiers,'DOSSIER_DIR',tmp_path);p=tmp_path/'40757665.json';original=json.dumps({'pmid':'40757665','abstract':'A cached scientific source','has_parsed_fulltext_xml':False});p.write_text(original)
    assert dossiers.build_one_dossier('40757665',{'pmid':'40757665'},[])==('40757665',False,0)
    assert p.read_text()==original


def test_dossier_truthy_string_cannot_claim_fulltext_availability(tmp_path,monkeypatch):
    monkeypatch.setattr(dossiers,'DOSSIER_DIR',tmp_path);p=tmp_path/'40757665.json';original=json.dumps({'pmid':'40757665','has_parsed_fulltext_xml':'false'});p.write_text(original)
    with pytest.raises(ValueError):dossiers.build_one_dossier('40757665',{'pmid':'40757665'},[])
    assert p.read_text()==original


def test_dossier_rejects_a_symlink_above_the_cache_directory(tmp_path,monkeypatch):
    outside=tmp_path/'outside';outside.mkdir();link=tmp_path/'alias';link.symlink_to(outside,target_is_directory=True)
    monkeypatch.setattr(dossiers,'DOSSIER_DIR',link/'dossiers')
    with pytest.raises(ValueError):dossiers.build_one_dossier('40757665',{'pmid':'40757665'},[])
    assert not (outside/'dossiers').exists()


@pytest.mark.parametrize('kind',['symlink','hardlink','fifo','directory'])
def test_dossier_rejects_linked_or_special_cached_files_before_fetch(tmp_path,monkeypatch,kind):
    monkeypatch.setattr(dossiers,'DOSSIER_DIR',tmp_path);p=tmp_path/'40757665.json';original=tmp_path/'original.json';original.write_text(json.dumps({'pmid':'40757665'}))
    if kind=='symlink':p.symlink_to(original)
    elif kind=='hardlink':os.link(original,p)
    elif kind=='fifo':os.mkfifo(p)
    else:p.mkdir()
    fetch=Mock();monkeypatch.setattr(dossiers,'fetch_epmc_fulltext_xml',fetch)
    with pytest.raises(ValueError):dossiers.build_one_dossier('40757665',{'pmid':'40757665'},[])
    fetch.assert_not_called()


def test_step1_convergence_swarm_invariants(tmp_path):
    """Regression coverage for Step-1 (Audit) 7-round gauntlet adjudication fixes."""
    from scripts import rank_islet1_stacks as ranker
    from scripts import ingest_wetlab_results as validator

    # 1. GRAVEYARD_003 triggers on SAF_G6PC2_ENDOGENOUS_ICASP9 as well as AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR (Opus F-19 / Audit M-01)
    gy_map = ranker._load_graveyard_entries_by_id()
    gy003 = gy_map["GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK"]
    req_mods = gy003["trigger_pattern"]["requires_any_module"]
    assert "SAF_G6PC2_ENDOGENOUS_ICASP9" in req_mods
    assert "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR" in req_mods
    for r in gy_map.values():
        assert "evaluator_logit_delta" in r and "evaluator_target_gate" in r
    trig_g6pc2_saf = ranker.evaluate_graveyard_rules_for_stack(
        "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO",
        "AUTO_NONE",
        "CYTO_CXCL10_SWAP_PDL1_SERPINB9",
        "ER_NONE",
        "NICHE_TRIHORMONAL_120UM_MVF_CAO2",
        "SAF_G6PC2_ENDOGENOUS_ICASP9",
        0,
    )
    trig_ids = {t["graveyard_id"] for t in trig_g6pc2_saf}
    assert "GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH" in trig_ids
    assert "GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK" in trig_ids

    # 2. Multi-batch wet-lab chaining compares final posterior against genesis priors & preserves Batch-1 shifts (Audit W-C1)
    ROOT = Path(__file__).resolve().parent.parent
    rep1 = validator.run_calibration(
        {
            "batch_id": "BATCH_01",
            "obs_stage7_yield_ratio": 1.045,
            "obs_stage7_yield_se": 0.020,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        },
        is_dry_run=False,
    )
    rep2 = validator.run_calibration(
        {
            "batch_id": "BATCH_02",
            "obs_stage7_yield_ratio": "UNAVAILABLE",
            "obs_cd47_fold": 8.40,
            "obs_cd47_fold_se": 0.30,
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        },
        is_dry_run=False,
        prior_posterior_report=rep1,
    )
    assert rep2["ingested_batch_ids"] == ["BATCH_01", "BATCH_02"]
    try:
        ctx1 = ranker.prepare_wetlab_evaluation_context(report=rep1)
        ctx2 = ranker.prepare_wetlab_evaluation_context(report=rep2)
        assert ctx1["raw_endpoint_posterior_deltas"]["stage7_sc_beta_yield_fraction_arm1b"] > 0.02
        # Batch 2 did not re-measure stage7 yield, yet its cumulative shift vs genesis prior (0.985) is preserved!
        assert abs(
            ctx2["raw_endpoint_posterior_deltas"]["stage7_sc_beta_yield_fraction_arm1b"]
            - ctx1["raw_endpoint_posterior_deltas"]["stage7_sc_beta_yield_fraction_arm1b"]
        ) < 1e-12
        assert ctx2["mapped_nonzero_parameter_deltas"]["stage7_trihormonal_120um_mvf_cao2_niche_prior"] > 0.02
        assert ctx2["mapped_nonzero_parameter_deltas"]["resting_b2m_promoter_hour0_density_discount"] > 0.01
    finally:
        ranker.prepare_wetlab_evaluation_context()

    # 3. Layer 4 Design X 100 um thinned central lamina (a = 100/350) & Catalase multi-denominator ledger (DEC-014 AUDIT-REF / Opus L4-13)
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pde = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    dx = pde["design_x_100um_lamina_transport_qualification_v80_018"]
    assert abs(dx["lamina_half_width_fraction_a"] - (100.0 / 350.0)) < 1e-6
    assert abs(dx["first_order_nz201"]["min_slab_po2_mmhg"] - 23.75) < 0.05
    assert abs(dx["first_order_nz201"]["gel_only_120um_core_po2_mmhg"] - 7.23) < 0.05
    cat_ledger = pde["catalase_multi_denominator_iu_ledger_v80_093"]
    assert abs(cat_ledger["denominators"]["lamina_aqueous_pore_volume_layout_b_0p41745_cm3"]["total_pouch_catalase_iu"] - 11271.2) < 1.0

    # 4. Layer 3 D227/T228 distance audit & SERPINB9 full mobile loop / hinge distances (Opus L3-01, L3-09)
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    d227_audit = c1["hlae_alpha3_d227_t228_distance_audit"]
    assert abs(d227_audit["min_heavy_atom_distance_range_to_cd94_nkg2a_heterodimer_angstrom"][0] - 47.00) < 0.05
    assert abs(d227_audit["min_heavy_atom_distance_range_across_individual_cd94_and_nkg2a_subunits_angstrom"][1] - 52.06) < 0.05
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    sb9 = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]
    assert abs(sb9["min_heavy_atom_distance_cterm_372_376_to_rcl_335_350_angstrom"] - 9.77) < 0.05
    assert abs(sb9["min_heavy_atom_distance_cterm_372_376_to_proximal_hinge_318_334_angstrom"] - 7.96) < 0.05
    assert abs(sb9["pro376_only_min_heavy_atom_distance_to_rcl_335_345_angstrom"] - 26.13) < 0.05


def test_sol61_convergence_round1_32_findings_invariants():
    """Verify all 32 Audit convergence round 1 findings (AUDIT-REF)."""
    import math
    from scripts import rank_islet1_stacks as ranker
    from scripts import ingest_wetlab_results as validator

    ROOT = Path(__file__).resolve().parent.parent

    # 1. AUDIT-REF: Spoofed genesis_prior_mean_exact or identity-posterior inflation is rejected
    rep1 = validator.run_calibration(
        {
            "batch_id": "BATCH_01",
            "obs_stage7_yield_ratio": 1.045,
            "obs_stage7_yield_se": 0.020,
            "obs_cd47_fold": "UNAVAILABLE",
            "obs_serpinb9_log10_kass": "UNAVAILABLE",
            "obs_icasp9_log10_kill": "UNAVAILABLE",
        },
        is_dry_run=False,
    )
    import hashlib
    spoofed = json.loads(json.dumps(rep1))
    spoofed["conjugate_posterior_updates"]["stage7_sc_beta_yield_fraction_arm1b"]["genesis_prior_mean_exact"] = 0.0
    spoofed["conjugate_posterior_updates_sha256"] = hashlib.sha256(
        json.dumps(spoofed["conjugate_posterior_updates"], sort_keys=True, allow_nan=False).encode("utf-8")
    ).hexdigest()
    spoofed["calibration_state_sha256"] = validator._calibration_state_sha256(spoofed)
    try:
        ctx_spoofed = ranker.prepare_wetlab_evaluation_context(report=spoofed)
        assert ctx_spoofed["applied_to_production_priors"] is False
        assert "genesis_prior_mean_exact" in ctx_spoofed["status_reason"]
    finally:
        ranker.prepare_wetlab_evaluation_context()

    try:
        validator._validate_prior_posterior_report(spoofed, False, "BATCH_02", "b" * 64)
        raise AssertionError("Expected ValueError for spoofed genesis_prior_mean_exact in ingest_wetlab_results")
    except ValueError as exc:
        assert "genesis_prior_mean_exact" in str(exc)

    # 2. AUDIT-REF: sigmoid_unclipped preserves IEEE-754 subnormal/normal floats in [-745, -700)
    assert ranker.sigmoid_unclipped(-705.0) > 0.0
    assert abs(ranker.sigmoid_unclipped(-705.0) - math.exp(-705.0)) < 1e-305

    # 3. AUDIT-REF: Layer 2 SV40pA 122/122, HNRNPA2B1 167 bp overlap, Arm 1B vs Arm 1X codon table, Chandra 2022 separation
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    cv = l2["construct_sequence_and_grammar_verification"]
    sv40_ver = cv["regulatory_element_audits"]["SV40_late_polyA_122bp_J02400"]["j02400_alignment_verification"]
    assert sv40_ver["exact_identity_bp"] == 122
    assert sv40_ver["query_length_bp"] == 122
    assert sv40_ver["range_1based_closed"] == "2580..2701"
    ucoe_audit = cv["regulatory_element_audits"]["CBX3_A2UCOE_core_455bp_R17"]["gencode_v46_enst00000354667_8_exon1_overlap_audit"]
    assert ucoe_audit["exact_overlap_length_bp"] == 167
    arm1b_entry = next(c for c in cv["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    codon_tbl = arm1b_entry["arm1b_vs_arm1x_codon_coordinate_table"]
    assert codon_tbl["cd47_cys15_mature_cys33_precursor_1based_span_in_arm1b_and_arm1x"]["donor_1based_coordinates"] == "1152..1154"
    assert codon_tbl["hlae_d227_t228_mature_aa389_390_precursor_in_arm1b"]["donor_1based_coordinates_in_3974bp_arm1b"] == "2814..2819"
    assert "2896..2901" in codon_tbl["hlae_d227_t228_mature_aa389_390_precursor_in_arm1b"]["donor_1based_coordinates_in_4056bp_arm1x_with_230bp_sense_promoter"]
    chandra_anchor = l2["stage0_to_stage7_perturbation_plane"]["external_benchmark_reference_anchors_not_fitted"][0]
    assert "published_source_measurements" in chandra_anchor and "rescaled_model_prior_baselines" in chandra_anchor

    # 4. AUDIT-REF: Layer 3 3CDG substitution table, 4ZQK Ala18..Ala132 audit, SERPINB9 PAE & Pro376.C distances
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    subs = c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis"]["nkg2c_3cdg_homolog_substitution_audit"]["ectodomain_substitutions_resolved_in_3cdg_chain_k"]
    assert any(s["nkg2a_res"] == "Ser167" and s["nkg2c_res"] == "Ala165" for s in subs)
    c4 = next(c for c in l3["complexes"] if "PD1" in c["complex_id"])
    assert "Ala18..Ala132" in c4["pdl1_to_pd1_interface_4ZQK"]["chain_a_resolved_residue_range_audit"]["pdb_4zqk_chain_a_resolved_span"]
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    sb9 = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]
    assert abs(sb9["pro376_c_peptide_bond_anchor_to_rcl_335_345_angstrom"] - 28.57) < 0.05
    assert abs(sb9["pro376_c_peptide_bond_anchor_to_rcl_335_350_angstrom"] - 15.07) < 0.05

    # 5. AUDIT-REF, 017, 028, 030: Layer 4 spatial localization audit, coupled 15 mmHg GSIS gate, lamina fraction guard
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pde = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    sp_audit = pde["catalase_multi_denominator_iu_ledger_v80_093"]["spatial_enzyme_localization_and_saturation_equivalence_audit"]
    assert abs(sp_audit["hour0_newton_bvp_spatial_and_dose_controls_uM"]["uniform_slab_k1_18p0_s_inv_189000_iu_bulk"] - 0.7249) < 0.01
    assert abs(sp_audit["hour0_newton_bvp_spatial_and_dose_controls_uM"]["uniform_slab_k1_1p0734_s_inv_11271_iu_spread_over_7cm3"] - 9.8158) < 0.01
    for row in l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["zero_sst_and_common_constants_qualification_ultracode_v80_10"]["coupled_pde_viability_to_metabolic_ode_scenarios"]:
        assert "coupled_survival_and_chemistry_screen_pass" in row
        assert "coupled_full_gsis_cluster_surface_15mmHg_qualified_pass" in row

    # 6. AUDIT-REF: Layer 5 unrounded analytic gap, dbSNP study allele count separation, Track 3 vs Track 1/2 SNP gate
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    p_sim = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]
    unrounded_rates = [c["analytic_expected_5yr_survival_personalized_compiled_unrounded_pct"] for c in p_sim["ancestry_stratified_cohorts"]]
    assert p_sim["analytic_max_ancestry_disparity_gap_personalized_pct_v80_065"] == round(max(unrounded_rates) - min(unrounded_rates), 2)
    dbsnp_sum = l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]["live_ncbi_dbsnp_polymorphism_audit"]
    for rec in dbsnp_sum.values():
        assert "untyped_study_allele_count_observations" in rec
        assert "derived_folded_minor_allele_frequency_summary" in rec


def test_kimi31_convergence_round1_findings_invariants():
    """Verify all Audit convergence round 1 findings (AUDIT-REF) and DEC-017."""
    from scripts import run_layer5_personalized_patient_compiler as l5_compiler

    ROOT = Path(__file__).resolve().parent.parent

    # 1. AUDIT-REF, 008, 012: Layer 2 D390/T391 SCT precursor key, dynamic polyA identity, single-variable challengers
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    cv = l2["construct_sequence_and_grammar_verification"]
    arm1b_entry = next(c for c in cv["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    codon_tbl = arm1b_entry["arm1b_vs_arm1x_codon_coordinate_table"]
    assert "hlae_d227_t228_mature_aa390_391_precursor_in_arm1b" in codon_tbl
    assert codon_tbl["hlae_d227_t228_mature_aa390_391_precursor_in_arm1b"]["donor_1based_coordinates_in_3974bp_arm1b"] == "2814..2819"
    assert "single_variable_challenger_arms_kimi31_q1" in codon_tbl
    challengers = codon_tbl["single_variable_challenger_arms_kimi31_q1"]
    assert "Arm_1B_CD47_C15G_single_nucleotide_challenger" in challengers
    assert "Arm_1B_L1L2_Synonymous_Linker_Diversified_challenger" in challengers
    assert (ROOT / "docs/STAGE-0-EXECUTION-PLAYBOOK.md").exists()

    # 2. AUDIT-REF, 009: Layer 3 reconciled SPR dictionary, Kaiser 0.1 SD, dynamic NKG2C distances
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    sel_analysis = c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis"]
    assert sel_analysis["literature_spr_reference_affinities"]["kaiser_2008_table1_kd_mutant_nkg2c_165_168siis_vmaprtlfl_sd_uM"] == 0.1
    sel_proof = c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_proof"]
    assert sel_proof["kd_cd94_nkg2a_vmaprtlvl_uM"] is None
    assert sel_proof["kd_cd94_nkg2a_vmaprtlfl_uM"] == 0.8
    assert sel_proof["kd_cd94_nkg2c_vmaprtlfl_uM"] == 5.2
    c2 = next(c for c in l3["complexes"] if "CD47" in c["complex_id"])
    assert "GRAVEYARD_009" in c2["gpi_anchor_and_sirpa_polymorphism_verification"]["graveyard_rule_verified"]

    # 3. AUDIT-REF, 004, 010: Layer 4 unrounded surface pO2, Design X coupled ODE scenarios & localized Catalase regime
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pde = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    dx = pde["design_x_100um_lamina_transport_qualification_v80_018"]
    assert "first_order_nz201_localized_catalase_22680_iu_pore_water_dose" in dx
    assert dx["first_order_nz201_localized_catalase_22680_iu_pore_water_dose"]["max_interstitial_h2o2_uM"] <= 5.0
    assert abs(dx["required_dry_cao2_weight_percent_in_pdms_beads"] - 46.39) < 0.05
    sp_audit = pde["catalase_multi_denominator_iu_ledger_v80_093"]["spatial_enzyme_localization_and_saturation_equivalence_audit"]
    assert "central_lamina_masked_k1_3p1572_s_inv_11271_iu_pore_water_dose" in sp_audit["hour0_newton_bvp_spatial_and_dose_controls_uM"]

    ode_coupled = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["zero_sst_and_common_constants_qualification_ultracode_v80_10"]["coupled_pde_viability_to_metabolic_ode_scenarios"]
    lead_238 = next(r for r in ode_coupled if r["pde_regime"] == "LEAD_CENTRAL_LAMINA_PDMS_CAO2_F045_PLUS_CATALASE")
    assert lead_238["coupled_survival_and_chemistry_screen_pass"] is True
    assert lead_238["coupled_complete_product_pass"] is False
    dx_201 = next(r for r in ode_coupled if r["pde_regime"] == "DESIGN_X_100UM_LAMINA_EPS042_FO_201")
    assert dx_201["coupled_complete_product_pass"] is True
    dx_loc = next(r for r in ode_coupled if r["pde_regime"] == "DESIGN_X_100UM_LAMINA_EPS042_FO_201_LOCALIZED_CATALASE")
    assert dx_loc["coupled_complete_product_pass"] is True

    # 4. AUDIT-REF: Layer 5 sirpa_genotype causal wiring and fail-closed validation
    base_pt = {
        "patient_id": "PT_TEST",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.20,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 50.0,
        "sirpa_genotype": "v1/v1",
        "cd47_sirpa_binding_compromised": True,
    }
    reg = l5_compiler.compile_personalized_patient_blueprint(base_pt)
    p_v1 = l5_compiler.evaluate_regimen_survival_probability(base_pt, reg)
    pt_v2 = dict(base_pt, sirpa_genotype="v2/v2")
    p_v2 = l5_compiler.evaluate_regimen_survival_probability(pt_v2, reg)
    assert p_v2 < p_v1
    try:
        l5_compiler.evaluate_regimen_survival_probability(dict(base_pt, sirpa_genotype="INVALID_SIRPA"), reg)
        raise AssertionError("Expected ValueError for invalid sirpa_genotype")
    except ValueError as exc:
        assert "sirpa_genotype" in str(exc)


def test_astra_convergence_round1_remediations_and_blueprint_upgrades() -> None:
    """Verify Step-4 Audit remediations and concrete blueprint upgrades."""
    import copy
    import json
    from pathlib import Path
    from scripts import ingest_wetlab_results as wetlab_ingest
    from scripts import run_layer4_tissue_and_metabolic_twin as l4_mod
    from scripts import run_layer5_personalized_patient_compiler as l5_compiler

    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Q2 2^3 Arm 1X factorial, Q3 40-run bridge option, D2-F04 noncentral-t power, D2-F05 stage_block_B SE/MDE
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx_spec = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        assert "stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1" in bx_spec
        sub_23 = bx_spec["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]
        assert sub_23["runs"] == 8 and sub_23["total_clonal_evaluations"] == 24
    hip_proto = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    assert "hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1" in hip_proto
    assert hip_proto["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]["full_29_column_model_rank"] == "29/29"
    ooc_pw = hip_proto["statistical_analysis_and_power_plan"]["organ_on_chip_stage2_6channel_power_qualification"]["two_sided_t_test_df10_ncp_2_887"]
    assert abs(ooc_pw["power_at_two_sided_alpha_0_05"] - 0.7413) < 1e-3
    assert abs(ooc_pw["power_at_n8_df14_two_sided_alpha_0_05"] - 0.8736) < 1e-3
    sb_col = (
        hip_proto["statistical_analysis_and_power_plan"]["minimum_detectable_effects_at_80pct_power"][
            "mode1_108_independent_run_nested_clones_exact_iid_run_covariance"
        ]["two_stage_23col_model_with_additive_stage_block_B"]["stage_block_B"]
    )
    assert abs(sb_col["full_stage_difference_beta_B_se_pp"] - 4.5055) < 1e-3
    assert abs(sb_col["full_stage_difference_beta_B_mde_80pct_alpha05_two_sided_df13_pp"] - 13.66) < 0.05

    # 2. Layer 2: L2-001..004 & BIO-01..02 (unthresholded longest_exact_repeat_bp, Sense Terminator Shield, Scar-Free Locus 3)
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    cv = l2["construct_sequence_and_grammar_verification"]
    loc3 = next(c for c in cv["constructs"] if c["construct_id"] == "LOCUS3_CXCL10_SERPINB9_TPT2A_PDL1_R01")
    assert loc3["donor_longest_exact_repeat_bp"] == 14
    assert loc3["cds_only_longest_exact_repeat_bp"] == 14
    arm1b = next(c for c in cv["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    challengers = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    assert "Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1" in challengers
    assert challengers["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]["donor_length_bp"] == 4128
    assert challengers["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]["legacy_unspaced_122bp_sv40_donor_length_bp"] == 4096
    assert "Arm_1A_lite_AAV6_compatible_challenger_sonnet55_r1" in challengers
    assert challengers["Arm_1A_lite_AAV6_compatible_challenger_sonnet55_r1"]["donor_length_bp"] == 3891
    assert "DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1" in challengers
    assert "B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1" in challengers
    assert "Locus3_ScarFree_Dual_Promoter_SERPINB9_PDL1_challenger_astra_r1" in challengers
    assert challengers["Locus3_ScarFree_Dual_Promoter_SERPINB9_PDL1_challenger_astra_r1"]["donor_length_bp"] == 3799
    c15g_desc = challengers["Arm_1B_CD47_C15G_single_nucleotide_challenger"]["plus_strand_mutation"]
    assert "A->C transversion at 1-based donor nucleotide 1154" in c15g_desc

    # 3. Layer 3: S2-02 Kaiser 2008 WT NKG2C SD = 0.2, D4-F01 Cheng 2013 soluble PD-L1:CD80, BIO-03 stalk reach
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    assert c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis"]["literature_spr_reference_affinities"]["kaiser_2008_table1_kd_wt_cd94_nkg2c_vmaprtlfl_sd_uM"] == 0.2
    assert c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_proof"]["kd_cd94_nkg2c_vmaprtlfl_sd_uM"] == 0.2
    c2 = next(c for c in l3["complexes"] if "CD47" in c["complex_id"])
    assert "Morrissey et al. 2020 PMID:32768386" in json.dumps(c2["gpi_anchor_and_sirpa_polymorphism_verification"])
    c4 = next(c for c in l3["complexes"] if "PDL1" in c["complex_id"])
    assert "Cheng et al. 2013 (PMID:23417675 Table 2)" in json.dumps(c4)

    # 4. Layer 4: D5-01 mesh refinement spatial catalase propagation, D5-02 fail-closed surface pO2 gate, BIO-04/05
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pde = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    dx = pde["design_x_100um_lamina_transport_qualification_v80_018"]
    assert "fixed_inventory_two_scale_ocr_headroom_sensitivity_astra_r1" in dx
    assert "cao2_72h_hydrolysis_calcium_alkali_stoichiometry_and_calcification_gate_astra_r1" in dx
    assert "dynamic_viability_trapezoidal_vs_right_endpoint_sensitivity_astra_r1" in dx
    # Verify D5-02: deleting passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only on a regime with surface pO2 < 15.0 mmHg still fails
    lead_reg_copy = copy.deepcopy(next(r for r in pde["regimes_simulated"] if r["regime_id"] == "LEAD_CENTRAL_LAMINA_PDMS_CAO2_F045_PLUS_CATALASE"))
    lead_reg_copy.pop("passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only", None)
    ode_forged = l4_mod.solve_24h_paracrine_and_whole_body_metabolic_ode(pde_regimes=[lead_reg_copy])
    forged_row = ode_forged["zero_sst_and_common_constants_qualification_ultracode_v80_10"]["coupled_pde_viability_to_metabolic_ode_scenarios"][0]
    assert forged_row["coupled_full_gsis_cluster_surface_15mmHg_qualified_pass"] is False
    assert forged_row["coupled_complete_product_pass"] is False

    # 5. Layer 5 & Wet-Lab Ingest: S2-03 non-string sirpa_genotype, D6-F01 Plan B product_definition, D6-F02 replay check
    base_pt = {
        "patient_id": "PT_ASTRA",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.52,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 50.0,
        "sirpa_genotype": "v1/v1",
    }
    for bad_sirpa in (123, ["V1/V1"], {"g": "V1/V1"}):
        try:
            l5_compiler.compile_personalized_patient_blueprint(dict(base_pt, sirpa_genotype=bad_sirpa))
            raise AssertionError(f"Expected ValueError for non-string sirpa_genotype={bad_sirpa!r}")
        except ValueError as exc:
            assert "sirpa_genotype" in str(exc)
    reg_b = l5_compiler.compile_personalized_patient_blueprint(base_pt)
    assert reg_b["selected_fda_deployment_track"] == "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE"
    assert "Plan B" in reg_b["plan_architecture_prescribed"]
    p_b_ok = l5_compiler.evaluate_regimen_survival_probability(base_pt, reg_b)
    assert p_b_ok > 0.90
    try:
        l5_compiler.evaluate_regimen_survival_probability(base_pt, dict(reg_b, product_definition="B2M_KO_WITHOUT_CD47"))
        raise AssertionError("Expected ValueError for non-dict product_definition")
    except ValueError as exc:
        assert "product_definition" in str(exc)


def test_sonnet55_convergence_round1_upgrades_and_dec020_tournament() -> None:
    """Verify Step-5 Audit upgrades and DEC-020 A/B/X/Y/Z tournament registry."""
    import json
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: DEC-020 A/B/X/Y/Z tournament registry & Q5 dual-clock governance
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        abxyz = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]
        assert "80%" in abxyz["review_focus_allocation_protocol"] and "20%" in abxyz["review_focus_allocation_protocol"]
        assert "Y_Radical_ZeroDSB_and_Active_Indirect_Cloak_Wildcard" in abxyz["tier2_radical_outside_the_box_wildcards_20pct_focus"]
        assert "Z_Radical_ZeroCaO2_Microvascular_Ribbon_Wildcard" in abxyz["tier2_radical_outside_the_box_wildcards_20pct_focus"]
        assert "q5_dual_clock_governance_and_screen_confirm_firewall_sonnet55_r1" in abxyz

    # 2. Layer 2: bidirectional SV40 polyA orientation note & non-palindromic sense terminator
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    sv40_j = l2["construct_sequence_and_grammar_verification"]["regulatory_element_audits"]["SV40_late_polyA_122bp_J02400"]["j02400_alignment_verification"]
    assert "bidirectional_sv40_polya_orientation_note_sonnet55_r1" in sv40_j
    assert "1 canonical AATAAA hexamer at local 1-based positions 78..83" in sv40_j["bidirectional_sv40_polya_orientation_note_sonnet55_r1"]
    assert sv40_j["donor_plus_strand_canonical_AATAAA_1based_spans"] == ["65..70", "94..99"]
    assert sv40_j["transcribed_antisense_cd47_minus_strand_canonical_AATAAA_1based_spans"] == ["78..83"]

    # 3. Layer 3: excluded-volume surface geodesic L1 linker audit & 20-aa (G4S)4 challenger
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert abs(geo["excluded_volume_surface_geodesic_3cdg_copy1_angstrom"] - 35.6) < 0.1
    assert geo["l1_20aa_g4s4_challenger_contour_length_at_3p5A_angstrom"] == 70.0

    # 4. Layer 4: host-gel avascular stand-off + masked-lamina perfusion (DX-b <=100 um cap), D_O2 time-lag gate, Erlang lag
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    standoff = dx["host_gel_avascular_standoff_and_masked_lamina_perfusion_sensitivity_sonnet55_r1"]
    assert standoff["dx_b_100um_cluster_diameter_cap_masked_perfusion_surface_po2_mmhg"] == 15.6
    assert "d_o2_hydrogel_provenance_and_stage0_time_lag_coupon_gate_sonnet55_r1" in dx
    assert "design_z_radical_zero_cao2_microvascular_pfc_ribbon_wildcard_spec_dec020" in dx
    ode_qual = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["parameter_literature_provenance_table_v80_091"]
    assert "extravascular_erlang_transport_lag_and_hepatic_partition_sensitivity_sonnet55_r1" in ode_qual

    # 5. Graveyard 005 & Layer 5: Lin et al. 2023 HLA-C*07:02 SP-3C VMAPRALLL non-functional leader-peptide note
    gy = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    g5 = next(r for r in gy["graveyard_entries"] if r["graveyard_id"].startswith("GRAVEYARD_005"))
    assert "VMAPRALLL" in g5["why_it_fails"]
    assert "PMC10690437" in g5["why_it_fails"]
    assert "Nature Immunology" in g5["why_it_fails"]


def test_muse_wonder_convergence_round1_remediations_and_upgrades() -> None:
    """Verify Step-6 Audit Full Stage-1 & Stage-2 findings and DEC-022 upgrades across Layers 1-5."""
    import json
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: 40-run explicit coded_levels, 6 bridge pairs, Var(beta_B)/sigma^2 = 0.302425, 29/29 rank, Arm 1X 4056 vs 4225 bp, off-grid Design X projection
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        assert bx["candidate_arithmetic_donor_length_derivation_bp"] == 4056
        assert bx["fully_loaded_arm1x_donor_length_with_spa_maz4_and_20aa_l1_bp"] == 4225
        assert bx["off_grid_molecular_design_x_locus1_and_lamina_counterfactual_projection_muse_r1"]["off_grid_design_x_constitutive_hlae_composite_utility_projection"] == 0.1162
        assert bx["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]["two_batch_hip_core_matched_sigma_run_pp"] == 5.4314

    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    opt40 = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]
    assert opt40["matched_stage0_stage1_bridge_genotype_pairs"] == 6
    assert abs(opt40["var_beta_B_over_sigma2_23col_in_40run_bridge"] - 0.302425) < 1e-5
    assert opt40["full_29_column_model_rank"] == "29/29"
    assert len(opt40["added_stage1_bridge_runs_37_to_40"]) == 4
    assert all(len(r["coded_levels"]) == 6 for r in opt40["added_stage1_bridge_runs_37_to_40"])

    # 2. Layer 2: DX-S iCasp9-GSG-P2A-CD47-IgV-GPI cistron order, 45-bp L1-in-L2 repeat, assembled canonical-part SHA-256 digests
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    arm1b = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    assert ch["Arm_1B_L1L2_Synonymous_Linker_Diversified_challenger"]["canonical_l1_in_l2_exact_identity_bp"] == 45
    assert ch["Arm_1B_L1L2_Synonymous_Linker_Diversified_challenger"]["diversified_l1_l2_longest_exact_repeat_bp"] == 9
    dxs = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]
    assert "iCasp9(nostop, 1,188 bp) - GSG-P2A(66 bp) - CD47-IgV-GPI(stop, 528 bp)" in dxs["cistron_topology_order_requirement_muse_r1"]
    assert len(dxs["locus4_icasp9_p2a_cd47_donor_sha256"]) == 64
    l3_sf = ch["Locus3_ScarFree_Dual_Promoter_SERPINB9_PDL1_challenger_astra_r1"]
    assert l3_sf["canonical_unspaced_part_sum_donor_length_bp"] == 3717
    assert l3_sf["donor_length_bp"] == 3799

    # 3. Layer 3: 8-copy chord range [17.96, 19.55], max-copy 43.7 A WLC tension 101.441 pN vs 22.955 pN
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert geo["straight_chord_range_across_8_crystal_copies_angstrom"] == [17.96, 19.55]
    assert abs(geo["l1_15aa_geodesic_wlc_metrics_at_max_copy_43p7A"]["marko_siggia_entropic_tension_pN_at_310K_lp4A"] - 101.441) < 0.1
    assert abs(geo["l1_20aa_g4s4_challenger_geodesic_wlc_metrics_at_max_copy_43p7A"]["marko_siggia_entropic_tension_pN_at_310K_lp4A"] - 22.955) < 0.1

    # 4. Layer 4: Two-tier D_O2 gate, joint adverse-corner stress audit, Design Z effective-medium PFOB transport, Control 5 2.0458 vs 2.0459 uM
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pde = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    dx = pde["design_x_100um_lamina_transport_qualification_v80_018"]
    assert "joint_adverse_corner_stress_audit_muse_r1" in dx
    dz = dx["design_z_radical_zero_cao2_microvascular_pfc_ribbon_wildcard_spec_dec020"]
    assert dz["effective_medium_transport_parameterization_muse_r1"]["effective_oxygen_permeability_multiplier_p_eff_over_p_aq"] == 2.50
    cat_ctrl = pde["catalase_multi_denominator_iu_ledger_v80_093"]["spatial_enzyme_localization_and_saturation_equivalence_audit"]["hour0_newton_bvp_spatial_and_dose_controls_uM"]
    assert abs(cat_ctrl["design_x_100um_lamina_masked_k1_7p5600_s_inv_22680_iu_pore_water_dose_hour0_nz51"] - 2.0458) < 5e-3
    assert abs(cat_ctrl["design_x_100um_lamina_masked_k1_7p5600_s_inv_22680_iu_336h_pde_nz201_max"] - 2.046) < 1e-2

    # 5. Layer 5: Gate 0 30x donor WGS & pre-cloning junction re-cleavage gate
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    assert "gate_zero_donor_wgs_and_precloning_junction_recleavage_gate_muse_r1" in l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]["post_hdr_junction_cfd_and_split_protospacer_audit_claude_v80_004"]


def test_argon_r2_convergence_round2_step1_remediations() -> None:
    """Verify Audit Convergence Round 2 Step 1 (DEC-023) remediations across Layers 1-5."""
    import json
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: D1-R2-01..03 (Runs 37-40 schema keys, Mode A vs Mode B kappa_2, 29-col residual_df=11 MDEs)
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    opt40 = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]
    for run in opt40["added_stage1_bridge_runs_37_to_40"]:
        for key in ("X1", "X2", "X3", "X4", "X5", "X6", "stage_block", "decoded_genotype"):
            assert key in run
        assert len(run["coded_levels"]) == 6
        assert "X4" in run["coded_levels"]
        assert "X4_TYK2_KO" in run["named_coded_levels"]
        assert "X5_TXNIP_KO" in run["named_coded_levels"]
    modes = opt40["operating_modes_and_condition_number_reconciliation_argon_r2"]
    assert modes["mode_a_primary_binary_bridge_Z3_plus1"]["rank_23col_with_stage_block_B"] == "23/23"
    assert modes["mode_b_quadratic_sentinel_dose_curvature_Z3_0"]["condition_number_24col_column_scaled"] == 7.42
    pwr29 = opt40["twenty_nine_column_interaction_power_diagnostics_argon_r2"]
    assert pwr29["residual_df_29col"] == 11
    assert pwr29["pure_error_df"] == 9
    assert pwr29["lack_of_fit_df_29col"] == 2

    # 2. Layer 2: D2-R2-01..04 (B2M_AsCas12a_S3 4-nt wobble coords, diversified L1/L2 9-bp repeat, Arm 1A-lite 225-bp bGH_polyA, DX-S 1.50x post-skip)
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    arm1b = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    s3_shield = ch["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]
    assert any("chr15:44711591T>C" in s for s in s3_shield["mutated_codons_in_b2m_exon1_leader"])
    assert any("seed position +4" in s for s in s3_shield["mutated_codons_in_b2m_exon1_leader"])
    assert "198..203" in s3_shield["eef1a1_5prime_top_tract_optimization"]
    a1a_lite = ch["Arm_1A_lite_AAV6_compatible_challenger_sonnet55_r1"]
    assert a1a_lite["cross_donor_locus3_longest_shared_block_canonical_parts_bp"] == 225
    assert a1a_lite["cross_donor_locus3_longest_shared_block_with_orthogonal_polya_and_diversified_2a_bp"] == 11
    dxs = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]
    assert "1.50x post-2A-skip" in dxs["scientific_rationale"]

    # 3. Layer 3: D3-R2-01..03, 06 (8 crystal copies explicitly listed, summary geodesic screen booleans, P02 PCA1 -17.03 Da check)
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    assert l3["summary"]["canonical_l1_15aa_surface_geodesic_screen_d_lt_0p65_Lc_passed"] is False
    assert l3["summary"]["challenger_l1_20aa_g4s4_surface_geodesic_screen_d_lt_0p65_Lc_passed"] is True
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert len(geo["per_copy_p9_c_to_b2m_ile1_n_angstrom"]) == 8
    assert geo["mean_straight_chord_across_8_crystal_copies_angstrom"] == 19.00
    assert "PCA1 pyroglutamate cyclization" in geo["stage0_intact_sct_and_signal_peptidase_gate_p02"]

    # 4. Layer 4: D4-R2-01..05 (Both nz51 and nz201 Hour-0 controls, pde_verified_corner_matrix_nz201_argon_r2, Pe_axial=12500 & Gz=51.3)
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pde = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    cat_ctrl = pde["catalase_multi_denominator_iu_ledger_v80_093"]["spatial_enzyme_localization_and_saturation_equivalence_audit"]["hour0_newton_bvp_spatial_and_dose_controls_uM"]
    assert abs(cat_ctrl["uniform_slab_k1_1p0734_s_inv_11271_iu_spread_over_7cm3"] - 9.8210) < 1e-3
    assert abs(cat_ctrl["nz201_fine_grid_controls_uM"]["uniform_slab_k1_1p0734_s_inv_11271_iu_spread_over_7cm3_nz201"] - 9.8158) < 1e-3
    assert abs(cat_ctrl["design_x_100um_lamina_masked_k1_7p5600_s_inv_22680_iu_pore_water_dose_hour0_nz201"] - 2.0459) < 1e-3
    dx = pde["design_x_100um_lamina_transport_qualification_v80_018"]
    corner_mat = dx["joint_adverse_corner_stress_audit_muse_r1"]["pde_verified_corner_matrix_nz201_argon_r2"]
    assert abs(corner_mat["baseline_pb40_cat1x_fo_d_floor_9p4e6"]["min_slab_mmhg"] - 23.70) < 0.05
    assert abs(corner_mat["masked_perf_standoff30um_ocr1p10x_pb40_d_target_1p5e5"]["min_slab_mmhg"] - 18.01) < 0.05
    dz = dx["design_z_radical_zero_cao2_microvascular_pfc_ribbon_wildcard_spec_dec020"]
    graetz = dz["effective_medium_transport_parameterization_muse_r1"]["axial_microvascular_advection_graetz_verification_argon_r2"]
    assert graetz["axial_peclet_number_pe_axial"] == 12500.0
    assert graetz["graetz_number_gz"] == 51.3

    # 5. Layer 5: D5-R2-01..05 (OI-GRNA-03 right-junction 14-bp spacer + AGG intact, CFD=0.247822, Plan B vs Plan B-prime structured keys)
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    j_audit = l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]["post_hdr_junction_cfd_and_split_protospacer_audit_claude_v80_004"]
    r_gate = j_audit["gate_zero_donor_wgs_and_precloning_junction_recleavage_gate_muse_r1"]
    assert "0.247822" in r_gate["pre_cloning_cxcl10_right_junction_recleavage_gate"]
    assert "0.143472" in r_gate["pre_cloning_cxcl10_right_junction_recleavage_gate"]
    assert "stage0_challenger_junction_equivalence_audit_argon_r2" in r_gate
    mfg_audit = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    assert mfg_audit["unmaterialized_plan_b_without_cd47_hlae_mean_survival_pct"] == 91.3735
    assert mfg_audit["conditional_on_release_qc_mean_personalized_survival_pct"] == 97.4963


def test_kimi31_r2_convergence_round2_step3_remediations() -> None:
    """Verify Audit Convergence Round 2 Step 3 (DEC-025) remediations across Layers 1-5."""
    import json
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: exact df=16 MDEs (5.7312 / 11.4624 pp), r=9 (0.603) vs r=15 (0.8228) power, Var(beta_B)/sigma^2 = 3/4 in 46-run, pan-lineage 93,350 sync
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        sub = bx["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]
        assert sub["main_effect_mde_exact_df16_noncentral_t_at_5p4314_pp"] == 5.7312
        assert sub["full_biological_contrast_2beta_mde_exact_df16_noncentral_t_at_5p4314_pp"] == 11.4624
        assert sub["clones_per_cell_exact_power_at_r9_under_two_batch_variance_model"] == 0.6030
        assert sub["clones_per_cell_exact_power_at_r15_df112_under_two_batch_variance_model"] == 0.8228
        assert sub["clones_per_cell_required_for_exact_80pct_power_under_two_batch_variance_model_r15"] == 15

        t_gate = l1["whole_product_clonogenic_safety_model"]
        assert t_gate["whole_product_proliferating_cells_total_N_prolif"] == 93350
        assert t_gate["whole_product_expected_surviving_clonogenic_cells_lambda"] == 0.0009335

    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    opt40 = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]
    ext46 = opt40["twenty_nine_column_interaction_power_diagnostics_argon_r2"]["unit_leverage_diagnostic_and_optional_46run_extension_sol61_r2"]["optional_46run_duplicated_bridge_extension"]
    assert ext46["stage_block_B_var_over_sigma2_29col_in_46run_exact"] == "3/4"
    modes = opt40["operating_modes_and_condition_number_reconciliation_argon_r2"]
    assert modes["mode_a_primary_binary_bridge_Z3_plus1"]["condition_number_22col_unscaled_with_intercept"] == 11.2549
    assert modes["mode_b_quadratic_sentinel_dose_curvature_Z3_0"]["condition_number_24col_B_pm1_raw"] == 15.9055

    # 2. Layer 2: All 7 challengers carry full sequences & live stem scan; pure S3 vs pure 5'TOP vs combined hashes; 66-bp post-insertion uORF audit
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    arm1b = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    s3 = ch["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]
    assert s3["pure_b2m_s3_4nt_only_donor_3974bp_sha256"] != s3["assembled_arm1b_s3_3974bp_donor_sha256"]
    assert s3["pure_eef1a1_5top_cagaga_only_donor_3974bp_sha256"] != s3["assembled_arm1b_s3_3974bp_donor_sha256"]
    term = ch["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]
    assert term["post_insertion_b2m_sense_uorf_audit_kimi31_r2"]["post_insertion_uorf_length_bp"] == 87
    assert term["post_insertion_b2m_sense_uorf_audit_kimi31_r2"]["post_insertion_uorf_peptide_aa"] == 28

    # 3. Layer 3: P9:C->B2M(Gln2:N) 22.56 A, V1 4CMM Asp101 insertion within 4.5A contact shell, exact_product_sirpa_v8_kd_uM is None
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert geo["per_copy_chord_and_4p0A_collision_audit_sol61_r2"]["3CDG_copy1_ABP"]["p9_c_to_b2m_gln2_n_angstrom"] == 22.5632
    c2 = next(c for c in l3["complexes"] if "CD47" in c["complex_id"])
    poly = c2["gpi_anchor_and_sirpa_polymorphism_verification"]
    assert poly["v1_4cmm_asp101_insertion_within_4p5A_cd47_contact_shell"] is True
    assert poly["exact_product_sirpa_v8_kd_uM"] is None

    # 4. Layer 4: Aligned [0.30 L, 1.20 L] Ca2+ retention arrays ([19.0331, 4.7583] mM & [96.85, 87.39]%)
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    eff = dx["cao2_72h_hydrolysis_calcium_alkali_stoichiometry_and_calcification_gate_astra_r1"]["effluent_mass_balance_and_pdms_retention_requirement_sol61_r2"]
    assert eff["sink_exchange_volume_over_72h_at_0p5_to_2p0_mL_cm2_day_liters"] == [0.30, 1.20]
    assert eff["unretained_100pct_export_added_ca2_plus_mM"] == [19.0331, 4.7583]
    assert eff["minimum_required_pdms_or_precipitate_ca2_plus_retention_pct"] == [96.85, 87.39]

    # 5. Layer 5: Discounted Plan B vs discounted Plan B-prime delta = -5.1901 pp
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    assert mfg["unmaterialized_plan_b_minus_materialized_plan_b_prime_complexity_discounted_pp"] == -5.1901
    assert mfg["unpenalized_plan_b_minus_discounted_plan_b_prime_cross_scale_delta_pp"] == -1.4628


def test_astra_r2_convergence_round2_step4_remediations() -> None:
    """Verify DEC-026 / AUDIT-REF Convergence Round 2 Step 4 remediations across Layers 1-5."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Primary r=15 clone requirement, 15-aa vs 20-aa L1 satellite, continuous Z3 matrix, and 93,350 common-mode bound
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        sub = bx["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]
        assert sub["clones_per_cell_required_for_standalone_80pct_power_on_5pp_full_contrast"] == 15
        assert sub["standalone_two_sided_power_for_5pp_full_contrast_at_n8_two_batch"] == 0.2558
        assert "matched_l1_linker_15aa_vs_20aa_satellite_contrast_astra_r2" in sub
        assert l1["whole_product_clonogenic_safety_model"]["quantitative_common_mode_failure_floor"]["maximum_common_mode_q_to_keep_pan_lineage_93350_risk_below_1e_3"] == 0.0000670

    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    opt40 = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]
    modes = opt40["operating_modes_and_condition_number_reconciliation_argon_r2"]
    assert modes["mode_a_primary_binary_bridge_Z3_plus1"]["condition_number_23col_B_pm1_raw"] == 14.7461
    assert modes["mode_a_primary_binary_bridge_Z3_plus1"]["max_hat_leverage_23col_exact_fraction"] == "8111/9566"
    assert modes["mode_b_quadratic_sentinel_dose_curvature_Z3_0"]["continuous_dose_polynomial_Z3_matrix_astra_r2"]["condition_number_24col_continuous_Z3_raw_B01"] == 16.2370
    assert opt40["twenty_nine_column_interaction_power_diagnostics_argon_r2"]["clone_level_df_120_clones_mode1_iid_120_minus_29"] == 91

    # 2. Layer 2: 87-bp (28-aa) uORF in Arm 1B Sense Terminator Shield & double-uAUG-ablated Arm 1C in DX-S
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    arm1b = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    uorf = ch["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]["post_insertion_b2m_sense_uorf_audit_kimi31_r2"]
    assert uorf["post_insertion_uorf_length_bp"] == 87
    assert uorf["post_insertion_uorf_peptide_sequence"] == "MSRSVALAVLALLSLSGLEAIIKDLYFH"
    dxs = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]
    assert dxs["promoter_trapped_post_splice_uorf_audit_and_rescue_astra_r2"]["post_splice_upstream_atg_count_between_tss509_and_hlae_atg610"] == 0

    # 3. Layer 3: Vales-Gomez 1999 10.77 uM NKG2C + VMAPRTLFL & Hatherley 2014 soluble 37C SPR affinities
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    spr = c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis"]["literature_spr_reference_affinities"]
    assert spr["vales_gomez_1999_table1_kd_cd94_nkg2c_vmaprtlfl_uM"] == 10.77
    c2 = next(c for c in l3["complexes"] if "CD47" in c["complex_id"])
    poly = c2["gpi_anchor_and_sirpa_polymorphism_verification"]
    assert poly["hatherley_2014_table2_soluble_spr_37C_kd_sirpa_v1_uM"] == 0.74
    assert poly["hatherley_2014_table2_soluble_spr_37C_kd_sirpa_v2_uM"] == 0.64

    # 4. Layer 4: Propagated Hour-16 adverse corner extrema and joint v=0.60 + Erlang lag qualification
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    t16 = dx["joint_adverse_corner_stress_audit_muse_r1"]["two_scale_propagated_and_full_trajectory_hour16_extrema_astra_r2"]
    assert t16["masked_perf_standoff30um_ocr1p10x_pb40_d_target_1p5e5_two_scale_propagated_t16h"]["r40_80um_surf_core_mmhg_at_t16h"] == [16.63, 11.07]
    bridge = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]
    assert bridge["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]["low_glucose_exercise_challenge_all_flux_scaled_v0p60_T9p5332_min_glucose_mM"] == 3.2734

    # 5. Layer 5: Dynamic patient-level complexity discount (87.6462% / -5.1901 pp)
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    assert mfg["unmaterialized_plan_b_complexity_discounted_mean_utility_pct"] == 87.6462
    assert mfg["unmaterialized_plan_b_minus_materialized_plan_b_prime_complexity_discounted_pp"] == -5.1901


def test_sonnet55_r2_convergence_round2_step5_remediations() -> None:
    """Verify DEC-027 / AUDIT-REF Convergence Round 2 Step 5 & live AlphaGenome/AlphaFold/UniProt remediations."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Embodiment A*, B' 2x2 satellite, off-grid conditionality (0.1085 vs 0.1162), exact 24col fractions (16180/19049 & 90371537/106195448), and 0.0000670 common-mode q
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        proj = bx["off_grid_molecular_design_x_locus1_and_lamina_counterfactual_projection_muse_r1"]
        assert proj["unconditional_composite_utility_if_graveyard_009_retained"] == 0.1085
        sub = bx["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]
        assert "arm1x_tier_endpoint_role_distinction_sonnet55_r2" in sub
        assert sub["mandatory_satellite_status"] == "MANDATORY_STAGE0_LINKER_DECONFOUNDING_CONTROL"
        triad = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]
        assert "embodiment_a_star_platform_upgrade_and_conditionality_sonnet55_r2" in triad["A_Champion_Islet1"]
        assert triad["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]["b_prime_stage0_isogenic_2x2_satellite_arm_sonnet55_r2"]["total_clones"] == 12
        assert l1["whole_product_clonogenic_safety_model"]["pan_lineage_proliferating_cell_budget_sol61_r2"]["maximum_common_mode_q_to_keep_pan_lineage_total_risk_below_1e_3"] == 0.0000670

    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    opt40 = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]
    modes = opt40["operating_modes_and_condition_number_reconciliation_argon_r2"]
    assert modes["mode_a_primary_binary_bridge_Z3_plus1"]["condition_number_23col_B_pm1_unit_column_norm_exact"] == 14.7461
    assert modes["mode_b_quadratic_sentinel_dose_curvature_Z3_0"]["max_hat_leverage_24col_mode_b_indicator_exact_fraction"] == "16180/19049"
    assert modes["mode_b_quadratic_sentinel_dose_curvature_Z3_0"]["continuous_dose_polynomial_Z3_matrix_astra_r2"]["max_hat_leverage_24col_continuous_Z3_exact_fraction"] == "90371537/106195448"

    # 2. Layer 2: U1 snRNA nt 600 T->G splice-donor rescue (TAG|GTAAGT 8/9 match, 0 uORFs) & unconditional 3-locus PAM shield + AlphaGenome live verification
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    arm1b = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    u1_rescue = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]["promoter_trapped_post_splice_uorf_audit_and_rescue_astra_r2"]["u1_snrna_consensus_600t_to_g_splice_donor_rescue_sonnet55_r2"]
    assert "TAG|GTAAGT" in u1_rescue["preferred_single_nt_600t_to_g_rescue_mutations"]
    assert u1_rescue["u1_consensus_600t_to_g_post_splice_upstream_atg_count_between_tss509_and_hlae_atg610"] == 0
    shield3 = ch["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]["unconditional_3locus_recut_immunity_pam_shield_bundle_sonnet55_r2"]
    assert shield3["deployment_status"] == "UNCONDITIONAL_PREVENTIVE_DONOR_DEFAULT_ACROSS_ALL_3_HDR_LOCI"
    assert "alphagenome_avi_live_verification_2026_10_27" in shield3

    # 3. Layer 3: Exact 3CDG copy-1 nearest non-peptide neighbors (A:Ser143:OG 3.27 A & A:Asp119:O 3.67 A), geodesic radius sensitivity, and RCL PAE qualification
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    nn = geo["covalent_anchor_nearest_non_peptide_neighbors_3cdg_copy1_sonnet55_r2"]
    assert nn["P9_Leu9_C_nearest_heavy_chain_A_atoms_angstrom"]["A:Ser143:OG"] == 3.27
    assert nn["B2M_Ile1_N_nearest_non_ile1_atoms_angstrom"]["A:Asp119:O_heavy_chain"] == 3.67
    assert "excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2" in geo
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    assert "rcl_low_plddt_and_high_pae_coordinate_uncertainty_qualification_sonnet55_r2" in c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]

    # 4. Layer 4: Ginstling-Brounshtein nz=201 localized-catalase PDE run (20.13 mmHg / 14.88 mmHg / 3.303 uM), per-lot OCR gate, and paracrine-preserved v=0.60 sweep
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    gb = dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]
    assert gb["min_slab_po2_mmhg"] == 20.13
    assert gb["gel_only_cluster_surface_po2_mmhg"] == 14.88
    assert gb["max_interstitial_h2o2_uM"] == 3.303
    assert dx["fixed_inventory_two_scale_ocr_headroom_sensitivity_astra_r1"]["per_lot_ocr_release_gate_sonnet55_r2"]["max_allowable_ocr_lot_over_ocr_model_for_120um_clusters"] == 1.097
    bridge = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]
    para_v060 = bridge["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]["paracrine_preserved_systemic_scaled_v060_erlang_sweep_sonnet55_r2"]
    assert para_v060[1]["low_glucose_challenge_min_glucose_0p25min_unrounded_mM"] == 3.5706
    assert para_v060[2]["low_glucose_challenge_min_glucose_0p25min_unrounded_mM"] == 3.3802

    # 5. Layer 5 & Wet-Lab Ingestor: Hsu 2013 MIT score (0.0037 / 0.0026), 4-corner prior box vs Plan A, and genesis-omission rejection
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    corners = mfg["four_corner_prior_box_vs_plan_a_comparison_sonnet55_r2"]
    assert len(corners) == 4
    assert corners[2]["personalized_exceeds_standardized_plan_a"] is False
    assert corners[3]["personalized_exceeds_standardized_plan_a"] is False

    from scripts.ingest_wetlab_results import run_calibration
    b1 = run_calibration({"batch_id": "BATCH-GEN-01", "obs_cd47_fold": 4.90, "obs_cd47_fold_se": 0.12}, is_dry_run=False)
    with pytest.raises(ValueError, match="AUDIT-REF"):
        run_calibration(
            {"batch_id": "BATCH-GEN-02", "obs_cd47_fold": 4.95, "obs_cd47_fold_se": 0.12, "product_variant_id": "ARM_1B"},
            is_dry_run=False,
            prior_posterior_report=b1,
        )
    with pytest.raises(ValueError, match="AUDIT-REF"):
        run_calibration(
            {"batch_id": "BATCH-GEN-02", "obs_cd47_fold": 4.95, "obs_cd47_fold_se": 0.12, "specimen_or_clone_ids": {"obs_cd47_fold": "CLONE-01"}},
            is_dry_run=False,
            prior_posterior_report=b1,
        )


def test_muse_wonder_r2_convergence_round2_step6_remediations() -> None:
    """Verify DEC-028 / AUDIT-REF Convergence Round 2 Step 6 fixed-point remediations."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Embodiment A* defended Champion status (4,143 / 4,433 bp) and standalone power approximation note (0.2318 / 0.2558 / 0.2562)
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        sub = bx["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]
        p_note = sub["standalone_two_sided_power_approximation_note_muse_r2"]
        assert p_note["exact_df16_noncentral_t_power"] == 0.2318
        assert p_note["normal_approx_power_unrounded_sigma_run_sqrt_177_over_6_5p431390"] == 0.2558
        assert p_note["normal_approx_power_rounded_sigma_run_5p43"] == 0.2559
        triad = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]
        a_star = triad["A_Champion_Islet1"]["embodiment_a_star_platform_upgrade_and_conditionality_sonnet55_r2"]
        assert a_star["defended_champion_status"] == "EMBODIMENT_A_STAR_PROMOTED_AS_DEFENDED_CHAMPION_CONFIGURATION"
        assert a_star["embodiment_a_star_locus1_donor_length_bp"] == 4143
        assert a_star["embodiment_a_star_locus1_packaged_aav6_length_bp"] == 4433

    # 2. Layer 2: 4-arm (12-clone) Arm 1C splice panel including constitutive 3,300-bp reference & AlphaGenome single-source provenance pin
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    arm1b = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = arm1b["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    u1_rescue = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]["promoter_trapped_post_splice_uorf_audit_and_rescue_astra_r2"]["u1_snrna_consensus_600t_to_g_splice_donor_rescue_sonnet55_r2"]
    assert u1_rescue["stage0_3arm_9clone_splicing_and_mfi_gate_sonnet55_r2"]["n_arms_including_constitutive_3300bp_reference_muse_r2"] == 4
    assert u1_rescue["stage0_3arm_9clone_splicing_and_mfi_gate_sonnet55_r2"]["n_clones_total_including_constitutive_3300bp_reference_muse_r2"] == 12
    avi = ch["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]["unconditional_3locus_recut_immunity_pam_shield_bundle_sonnet55_r2"]["alphagenome_avi_live_verification_2026_10_27"]
    assert avi["provenance_and_evidence_class_muse_r2"]["evidence_class"] == "LIVE_VERIFIED_BY_PRIVATE_API_RECEIPT_9_EDITS_PLUS_50_ROW_SAMPLE"
    assert avi["chr15:44711591:T>C_B2M_S3_Leu15_PAM"]["verification_status"] == "verified against the private receipt"
    receipt_path = ROOT / "data/raw/alphagenome_avi_live_receipt_2026_10_28.json"
    if receipt_path.is_file():
        rec = json.loads(receipt_path.read_text(encoding="utf-8"))
        assert rec["tool_provenance"]["archive_6decimal_match_summary"]["scoped_9_edits_matched"] == "9/9"
        assert rec["tool_provenance"]["archive_6decimal_match_summary"]["random_50_archive_rows_matched"] == "50/50"

    # 3. Layer 3: Structured PMID-anchored SIRPa SPR ledger & coordinate-derived RCL pLDDTs
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c2 = next(c for c in l3["complexes"] if "CD47" in c["complex_id"])
    sirpa_rec = c2["gpi_anchor_and_sirpa_polymorphism_verification"]["sirpa_temperature_and_immobilization_affinity_reconciliation_sonnet55_r2"]
    assert len(sirpa_rec["pmid_anchored_spr_measurement_ledger"]) == 4
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    assert "53.30" in c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]["rcl_low_plddt_and_high_pae_coordinate_uncertainty_qualification_sonnet55_r2"]

    # 4. Layer 4: Ginstling-Brounshtein 120-um surface failure flag, joint Ginstling corner, tightened t50<=4.35 min gate, and RK4 step-halving audit
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    gb = dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]
    assert gb["passes_15mmHg_surface_gsis_gate_at_120um"] is False
    assert gb["dxb_100um_cluster_surface_and_core_po2_gel_only_mmhg"] == [16.48, 8.64]
    assert "joint_ginstling_ocr1p097x_pb32p5_corner_muse_r2" in gb
    bridge = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]
    assert bridge["tightened_stage2_washout_target_for_v060_exercise_challenge_t50_le_4p35_min"]["max_washout_t50_min"] == 4.35
    assert bridge["rk4_step_halving_refinement_audit_muse_r2"]["verified_nadir_discretization_drift_within_0p0002_mM"] is True

    # 5. Layer 5: Off-target heuristic freeze (hsu_2013_mit_off_target_score is None), HLA-Cw*0402 >=75.30 uM lower bound, SP-2C inference tag, and SP-3C sensitivity grid
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    j_audit = l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]["post_hdr_junction_cfd_and_split_protospacer_audit_claude_v80_004"]
    cxcl10_r = j_audit["OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9"]["ISLET1_LOCUS3_CXCL10_SERPINB9_GSGP2A_PDL1_NO_FURIN_4090bp"]["three_prime_insert_ha_right_junction"]
    assert cxcl10_r["hsu_2013_mit_off_target_score"] is None
    assert cxcl10_r["published_cfd_or_hsu_mit_matrix_status"] == "UNAVAILABLE_UNCALIBRATED_POSITIONAL_HEURISTIC_ONLY"
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    hlac_tbl = mfg["hla_c_allele_leader_peptide_and_nkg2a_kd_table_sonnet55_r2"]
    assert hlac_tbl[1]["vales_gomez_1999_cd94_nkg2a_kd_lower_bound_uM"] == ">= 75.30"
    assert hlac_tbl[2]["endogenous_hla_e_nkg2a_engageable_evidence_status"] == "DIRECT_SOLUBLE_HLAE_SPR_LIN_2023_SUPP_TABLE1_P6_THR_P8_LEU_INTACT_SCT_UNMEASURED"
    assert len(mfg["sp3c_uncloaked_survival_penalty_sensitivity_grid_muse_r2"]) == 4


def test_argon_r3_convergence_round3_step1_remediations() -> None:
    """Verify DEC-029 (Audit 5-Domain Swarm Review, Convergence Round 3 Step 1) across Layers 1-5 & ingest_wetlab_results."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Exact composite and pediatric counterfactual trajectories on Blueprint B
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        reg = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]
        assert "DEC-029" in reg["decision_ids"]
        triad = reg["tier1_stage0_execution_triad_80pct_focus"]
        cf = triad["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]["graveyard_005_and_010_counterfactual_sensitivity_sonnet55_r2"]
        assert cf["exact_composite_trajectory_50_50_ped_lada"] == [0.0673, 0.0623, 0.0589]
        assert cf["exact_pediatric_trajectory"] == [0.0664, 0.0612, 0.0579]

    # 2. Layer 2: Unconditional 3-locus PAM+seed shield bundle HA_R coordinate precision (nt 17 CXCL10, nt 14 AAVS1)
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    loci_s = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]["unconditional_3locus_recut_immunity_pam_shield_bundle_sonnet55_r2"]["loci_shielded"]
    assert "CXCL10_HA_R_600bp_minus nt 17" in loci_s[1]
    assert "AAVS1_HA_R_600bp_minus nt 14" in loci_s[2]

    # 3. Layer 3: Exact 3CDG copy-1 B:Ile1:N neighbors (3.58 A B:Gln2:N), deduplicated capsule collisions (40; [36, 42]), and RCL 335-345 min/max
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    nbrs = geo["covalent_anchor_nearest_non_peptide_neighbors_3cdg_copy1_sonnet55_r2"]
    assert nbrs["B2M_Ile1_N_nearest_non_ile1_atoms_angstrom"]["B:Gln2:N_adjacent_backbone"] == 3.58
    assert nbrs["ile1_resolved_6_copies_3cdg_3bze_chord_range_angstrom"] == [19.03, 19.55]
    assert geo["per_copy_chord_and_4p0A_collision_audit_sol61_r2"]["3CDG_copy1_ABP"]["capsule_0_le_t_le_1_4p0A_heavy_atoms_excluding_p9_and_b2m_0_1"] == 40
    assert "40 deduplicated heavy atoms in 3CDG copy 1 (range 36-42 across the 8 copies" in geo["excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2"]["open_chord_cylinder_vs_capsule_collision_convention_note"]
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    rcl_note = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]["rcl_low_plddt_and_high_pae_coordinate_uncertainty_qualification_sonnet55_r2"]
    assert "max 59.00 at Ser345" in rcl_note
    assert "26.13 A [28.57 A from Pro376.C]" in rcl_note

    # 4. Layer 4: Exact BVP vs linear DX-b pO2, Ginstling coupled PDE regime, Blueprint A 100-um BVP, and v=0.60 standard-protocol table sync
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    gb = dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]
    assert gb["dxb_100um_cluster_surface_and_core_po2_gel_only_exact_bvp_mmhg"] == [16.44, 8.63]
    assert gb["joint_ginstling_ocr1p097x_pb32p5_corner_muse_r2"]["scaled_1p097x_bvp_surface_and_core_po2_120_100_80um_mmhg"]["100um"] == [6.51, 0.21]
    bpa_rec = dx["blueprint_a_outer_band_115p5um_vs_100um_cluster_cap_reconciliation_sonnet55_r2"]
    assert bpa_rec["blueprint_a_100um_cluster_surface_and_core_po2_outer_annulus_layout_b_mmhg"] == [20.62, 12.72]
    part_b = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]
    coupled_ids = [r["pde_regime_id"] for r in part_b["zero_sst_and_common_constants_qualification_ultracode_v80_10"]["coupled_pde_viability_to_metabolic_ode_scenarios"]]
    assert "DESIGN_X_100UM_LAMINA_EPS042_GB_201_LOCALIZED_CATALASE" in coupled_ids
    bridge = part_b["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]
    assert bridge["standard_protocol_v0p60_T9p5332_min_glucose_mM"] == 4.59
    assert bridge["standard_protocol_v0p60_T9p5332_max_postprandial_glucose_mM"] == 8.25
    assert bridge["coupled_ginstling_surface_gsis_times_erlang_t50_8min_exercise_corner_muse_r2"]["linear_threshold_15mmHg_gsis_fraction_at_14p88_mmhg"] == 0.992

    # 5. Layer 5: Exact 10k-patient simulated KIR routing counts (4982, 3583, 2367, 1148) & population SP-3C sensitivity grid
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    kir_sens = mfg["plan_b_kir_nk_threshold_sensitivity_sonnet55_r2"]
    assert [row["plan_b_routed_patients_of_10000"] for row in kir_sens] == [4982, 3583, 2367, 1148]
    sp3c_grid = mfg["sp3c_uncloaked_survival_penalty_sensitivity_grid_muse_r2"]
    assert [row["population_10k_uncloaked_plan_b_mean_survival_pct"] for row in sp3c_grid] == [90.9001, 90.4267, 89.4799, 87.8230]


def test_sol61_convergence_round3_step2_dec030_remediations() -> None:
    """Verify DEC-030 (Audit Convergence Round 3 Step 2) across Layers 1-5 & wet-lab contract."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Split-batch S1-A/S1-B GLS covariance, 0.1161 vs 0.1162 utility, 0.2559 rounded power, 898.2M physical cells, and 7.4225 vs 16.7269 condition numbers
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        bx = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        assert bx["off_grid_molecular_design_x_locus1_and_lamina_counterfactual_projection_muse_r1"]["off_grid_design_x_constitutive_hlae_composite_utility_projection_exact_unrounded_inputs"] == 0.1161
        sub8 = bx["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]["standalone_two_sided_power_approximation_note_muse_r2"]
        assert sub8["normal_approx_power_rounded_sigma_run_5p43_exact"] == 0.2559
        assert sub8["clone_only_sigma_eps_zero_sigma_run_sqrt_64_over_3_pp"] == 4.6188
        assert sub8["single_batch_with_assay_residual_sigma_run_sqrt_113_over_3_pp"] == 6.1373
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    mcb = hip_p["whole_product_clonogenic_safety_model"]
    assert mcb["required_physical_cells_plated_for_literal_300M_effective_cells_at_33p4pct"] == 898203593
    assert mcb["quantitative_common_mode_failure_floor"]["nonparametric_zero_failure_batches_required_for_95pct_upper_bound_6p70e_5"] == 44713
    mode_b = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]["operating_modes_and_condition_number_reconciliation_argon_r2"]["mode_b_quadratic_sentinel_dose_curvature_Z3_0"]
    assert mode_b["condition_number_15col_2fi_subblock_exact"] == 7.4225
    assert mode_b["condition_number_24col_euclidean_unit_column_norm_scaled"] == 16.7269
    b_sweep = hip_p["statistical_analysis_and_power_plan"]["minimum_detectable_effects_at_80pct_power"]["stage_by_main_interaction_and_batch_covariance_qualification_ultracode_v80_12"]["differentiation_batch_random_intercept_sensitivity_sweep"]["batch_sd_2p5_pp"]
    assert "stage1_split_batch_s1a_vs_s1b_incidence_2batch_residual_gls_se_sol61_r3" in b_sweep

    # 2. Layer 2: Materialized 4,143-bp Embodiment A* donor & SHA-256, splicing efficiency SE = S/(S+U) >= 0.85, and local-archive AlphaGenome tag
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    a_star = ch["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]
    assert len(a_star["embodiment_a_star_combined_4143bp_donor_sequence"]) == 4143
    assert len(a_star["embodiment_a_star_combined_4143bp_donor_sha256"]) == 64
    u1_gate = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]["promoter_trapped_post_splice_uorf_audit_and_rescue_astra_r2"]["u1_snrna_consensus_600t_to_g_splice_donor_rescue_sonnet55_r2"]["stage0_3arm_9clone_splicing_and_mfi_gate_sonnet55_r2"]
    assert u1_gate["splicing_efficiency_se_s_over_s_plus_u_minimum"] == 0.85
    assert u1_gate["splicing_efficiency_se_high_benchmark_for_lt_5pct_intron_retention"] == 0.95
    ag_prov = ch["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]["unconditional_3locus_recut_immunity_pam_shield_bundle_sonnet55_r2"]["alphagenome_avi_live_verification_2026_10_27"]["provenance_and_evidence_class_muse_r2"]
    assert ag_prov["evidence_class"] == "LIVE_VERIFIED_BY_PRIVATE_API_RECEIPT_9_EDITS_PLUS_50_ROW_SAMPLE"

    # 3. Layer 3: 1MHE copy-2 18.00 A chord, 37+5+2=44 vs 40 capsule collision breakdown, directional PAE, and PCA1/PGP note
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert "17.96 / 18.00 A [17.9601 A in 1MHE copy 1 and 18.0036 A in 1MHE copy 2]" in geo["stage0_intact_sct_and_signal_peptidase_gate_p02"]
    assert "37 + 5 + 2 = 44 [range 39-46 across the 8 copies]" in geo["excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2"]["open_chord_cylinder_vs_capsule_collision_convention_note"]
    c2 = next(c for c in l3["complexes"] if "CD47" in c["complex_id"])
    assert "crystal_2JJS_chainC_pca1_edman_and_lcms_note_sol61_r3" in c2["designed_vs_crystallized_cd47_sequence_identity_audit_ultracode_v80_11"]
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    assert "15.20-27.29 A AlphaFold directional expected positional error (PAE" in c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]["rcl_low_plddt_and_high_pae_coordinate_uncertainty_qualification_sonnet55_r2"]

    # 4. Layer 4: Forwarded d_hydrogel_override BVP (0.36 mmHg at D_floor, 11.35 mmHg at D_target), GB stock/source conservation reconciliation, and strict 4.34 min washout
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    gb = dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]
    assert abs(gb["stock_vs_regularized_source_conservation_reconciliation_sol61_r3"]["early_hour_clamp_stock_minus_source_delta_fraction"] - 0.070295) < 1e-6
    ext16 = dx["joint_adverse_corner_stress_audit_muse_r1"]["two_scale_propagated_and_full_trajectory_hour16_extrema_astra_r2"]
    assert ext16["masked_perf_standoff30um_ocr1p10x_pb40_d_floor_9p4e6_two_scale_propagated_t16h"]["r40_80um_surf_core_mmhg_at_t16h"] == [6.57, 1.79]
    assert ext16["masked_perf_standoff30um_ocr1p10x_pb40_d_target_1p5e5_two_scale_propagated_t16h"]["r40_80um_surf_core_mmhg_at_t16h"] == [16.63, 11.07]
    bridge = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]
    assert bridge["low_glucose_exercise_challenge_all_flux_scaled_v0p60_T0_min_glucose_mM"] == 3.8104
    assert bridge["tightened_stage2_washout_target_for_v060_exercise_challenge_t50_le_4p35_min"]["strict_unrounded_all_flux_max_washout_t50_min"] == 4.34

    # 5. Layer 5 & Wet-Lab Contract: HLA-C*06:02 in SP-1C, U1 conflicting specimen guard, and U2 historical measurement metadata guard
    from scripts import ingest_wetlab_results as wetlab

    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    hlac_tbl = mfg["hla_c_allele_leader_peptide_and_nkg2a_kd_table_sonnet55_r2"]
    assert "HLA-C*06:02" in hlac_tbl[0]["hla_c_allele_or_group"]
    assert "HLA-C*06:02" not in hlac_tbl[2]["hla_c_allele_or_group"]
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-U1-CONFLICT",
                "specimen_or_clone_id": "CLONE-A",
                "specimen_or_clone_ids": {"obs_stage7_yield_ratio": "CLONE-B"},
                "obs_stage7_yield_ratio": 0.90,
                "obs_stage7_yield_se": 0.05,
            },
            is_dry_run=True,
        )
        raise AssertionError("Expected ValueError for conflicting specimen_or_clone_id and specimen_or_clone_ids")
    except ValueError as exc:
        assert "AUDIT-REF" in str(exc)


def test_kimi31_convergence_round3_step3_dec031_remediations() -> None:
    """Verify DEC-031 (Audit Convergence Round 3 Step 3) remediations across Layers 1-5."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Mode 1 (108 independent clones) vs Mode 2 (75 shared-genotype clones) split-batch GLS SEs, live counterfactuals, and 44713 vs 44711 batches
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    b_sweep = hip_p["statistical_analysis_and_power_plan"]["minimum_detectable_effects_at_80pct_power"]["stage_by_main_interaction_and_batch_covariance_qualification_ultracode_v80_12"]["differentiation_batch_random_intercept_sensitivity_sweep"]["batch_sd_2p5_pp"]
    m2_2b = b_sweep["stage1_split_batch_s1a_vs_s1b_incidence_2batch_residual_gls_se_sol61_r3"]
    m1_2b = b_sweep["stage1_split_batch_s1a_vs_s1b_mode1_independent_run_nested_clones_2batch_residual_gls_se_kimi31_r3"]
    assert m1_2b["main_effect_half_effect_beta_gls_se_range_pp"] == [1.2359, 1.3086]
    assert m1_2b["two_factor_interaction_half_effect_beta_gls_se_range_pp"] == [1.4499, 4.6648]
    assert m1_2b["stage_block_B_gls_se_pp"] == 5.154
    assert m2_2b["main_effect_half_effect_beta_gls_se_range_pp"] == [1.2655, 1.3461]
    assert m2_2b["two_factor_interaction_half_effect_beta_gls_se_range_pp"] == [1.8054, 4.0872]
    assert m2_2b["stage_block_B_gls_se_pp"] == 3.7139
    cm = hip_p["whole_product_clonogenic_safety_model"]["quantitative_common_mode_failure_floor"]
    assert cm["nonparametric_zero_failure_batches_required_for_95pct_upper_bound_exact_qmax"] == 44713
    assert cm["nonparametric_zero_failure_batches_required_for_95pct_upper_bound_literal_6p70e_5"] == 44711

    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        pf_x = l1["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]
        assert pf_x["off_grid_molecular_design_x_locus1_and_lamina_counterfactual_projection_muse_r1"]["computation_mode"] == "LIVE_EVALUATOR_DERIVED_KIMI31_R3"
        triad = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]
        assert triad["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]["graveyard_005_and_010_counterfactual_sensitivity_sonnet55_r2"]["computation_mode"] == "LIVE_EVALUATOR_DERIVED_KIMI31_R3"
        assert "isogenic_leader_swap_c1_c2_kir_deconfounding_control_kimi31_r3" in triad["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]["b_prime_stage0_isogenic_2x2_satellite_arm_sonnet55_r2"]

    # 2. Layer 2: Materialized post-splice 101-nt 5'-UTR and 2,218-nt mRNA for Arm 1C-599_600_CC and Arm 1C-U1-600G
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    spliced_node = ch["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]["promoter_trapped_post_splice_uorf_audit_and_rescue_astra_r2"]["u1_snrna_consensus_600t_to_g_splice_donor_rescue_sonnet55_r2"]["materialized_post_splice_transcripts_and_retired_642nt_note_kimi31_r3"]
    assert len(spliced_node["arm1c_u1_600g_spliced_5utr_101nt_sequence"]) == 101
    assert len(spliced_node["arm1c_599_600_cc_spliced_5utr_101nt_sequence"]) == 101
    assert len(spliced_node["arm1c_u1_600g_spliced_mrna_2218nt_sequence"]) == 2218
    assert len(spliced_node["arm1c_599_600_cc_spliced_mrna_2218nt_sequence"]) == 2218

    # 3. Layer 3: Executable 26-connected 3D voxel Dijkstra geodesic solver on 3CDG copy 1
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    vox = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]["excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2"]
    assert vox["executable_3d_voxel_dijkstra_solver_verified_kimi31_r3"] is True
    assert abs(vox["live_3d_voxel_dijkstra_3cdg_copy1_r_obs_2p5A_angstrom"] - 33.8) < 0.05

    # 4. Layer 4: Exact piecewise mass-conserving Ginstling-Brounshtein PDE run, z=50 um interface H2O2, and live 9-state ODE washout runs
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    gb = dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]
    cons = gb["stock_vs_regularized_source_conservation_reconciliation_sol61_r3"]
    assert cons["early_window_regularization_deficit_fraction_of_336h_classical_release"] == 0.078471
    assert cons["exact_piecewise_conserving_pde_run_kimi31_r3"]["stock_minus_source_conservation_error"] == 0.0
    assert abs(gb["max_interstitial_h2o2_at_lamina_interface_z50um_uM"] - 2.96) < 0.05
    bridge = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]
    t_gate = bridge["tightened_stage2_washout_target_for_v060_exercise_challenge_t50_le_4p35_min"]
    assert t_gate["strict_raw_glucose_washout_gate_passed_at_T5p17_t50_4p34"] is True
    assert t_gate["strict_raw_glucose_washout_gate_passed_at_T5p179_t50_4p35"] is False

    # 5. Layer 5 & Graveyard: Same-scale K=4,B=0 discounted Plan A comparison, product_contract_scope_kimi31_r3, and C*15:05 in GRAVEYARD_005
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    box = mfg["four_corner_prior_box_vs_plan_a_comparison_sonnet55_r2"]
    assert box[0]["standardized_plan_a_k4_b0_discount_factor"] == 1.0
    assert box[0]["standardized_plan_a_discounted_utility_pct"] == 92.5019
    assert "product_contract_scope_kimi31_r3" in mfg["sp3c_uncloaked_survival_penalty_sensitivity_grid_muse_r2"][0]
    gy = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    gy5 = next(e for e in gy["graveyard_entries"] if e["graveyard_id"].startswith("GRAVEYARD_005"))
    assert "C*15:05" in json.dumps(gy5)


def test_astra_ultracode_convergence_round3_step4_dec033_remediations() -> None:
    """Verify DEC-033 (Audit Convergence Round 3 Step 4) remediations across Layers 1-5 & wet-lab contract."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: Runs 37-40 decoded_genotype & locus1_physical_cassette_id, two-sample external control SE, and Clopper-Pearson eta_L95 MCB cells
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    added40 = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]["added_stage1_bridge_runs_37_to_40"]
    assert [r["locus1_physical_cassette_id"] for r in added40] == [
        "LOC1_X1P_X2P_HEF1A_HLAE_TPT2A_CD47GPI_BICISTRONIC",
        "LOC1_X1P_X2M_B2M_HLAE_TPT2A_CD47WT_BICISTRONIC",
        "LOC1_X1M_X2P_HEF1A_CD47GPI_MONOCISTRONIC",
        "LOC1_X1P_X2M_B2M_HLAE_TPT2A_CD47WT_BICISTRONIC",
    ]
    assert all("AAVS1_CD47" not in r["decoded_genotype"] and "X2_CD47" in r["decoded_genotype"] for r in added40)
    mcb_d = hip_p["whole_product_clonogenic_safety_model"]["clonogenic_validation_tier_distinction"]
    assert mcb_d["one_time_mcb_qualification"]["required_physical_cells_plated_at_clopper_pearson_lower_95pct_efficiency_26p961pct"] == 1111135446
    assert mcb_d["routine_per_lot_sub_lot_screen"]["required_physical_cells_plated_at_33p4pct_detection_efficiency"] == 29940120

    # 2. Layer 2: Materialized 4,128-bp SPA+S3+15aaL1 satellite control donor and corrected sense HLA-E coordinates (1801..3318)
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    a_star = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]
    assert len(a_star["matched_4128bp_spa_s3_15aa_l1_satellite_control_donor_sequence"]) == 4128
    assert a_star["sense_hla_e_sct_coordinate_table_in_4143bp_and_4128bp_donors_astra_r3"]["HLAE_SCT_ORF_4143bp_20aa_L1_1based"] == "1801..3318"

    # 3. Layer 3: Official AlphaFold PAE row=aligned / col=scored convention, Lin 2023 Supp Table 1 soluble HLA-E SPR, and Hatherley 2014 [0.64-0.78 uM]
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    spr_lin = c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_proof"]["lin_2023_supp_table1_soluble_peptide_hlae0103_spr_reference_astra_r3"]
    assert spr_lin["soluble_hlae0103_vmaprtlvl_cd94_nkg2a_kd_uM"] == 2.10
    assert spr_lin["soluble_hlae0103_vmaprtlll_sp2c_cd94_nkg2a_kd_uM"] == 1.70
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    pae_a3 = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]["alphafold_pae_uncertainty_summary_claude_v80_065"]["pae_row_aligned_col_scored_convention_astra_r3"]
    assert pae_a3["pae_rows_rcl_335_345_aligned_cols_core_1_322_scored_mean_angstrom"] == 27.29
    assert pae_a3["pae_rows_core_1_322_aligned_cols_rcl_335_345_scored_mean_angstrom"] == 15.20

    # 4. Layer 4: Exact 0.825513 clamped GB source integral (0.070295 delta), coupled_complete_product_pass gating on physical v=0.60, and RK4 non-smooth RHS audit
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    cons = dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]["stock_vs_regularized_source_conservation_reconciliation_sol61_r3"]
    assert cons["early_hour_pdms_shell_clamp_336h_cumulative_source_integral_fraction"] == 0.825513
    assert cons["early_hour_clamp_stock_minus_source_delta_fraction"] == 0.070295
    coupled_rows = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["zero_sst_and_common_constants_qualification_ultracode_v80_10"]["coupled_pde_viability_to_metabolic_ode_scenarios"]
    dx_gb_row = next(r for r in coupled_rows if r["pde_regime_id"] == "DESIGN_X_100UM_LAMINA_EPS042_GB_201_LOCALIZED_CATALASE")
    assert dx_gb_row["coupled_complete_product_pass"] is False
    assert dx_gb_row["surface_gsis_secretion_coupled_trajectory_sonnet55_r2"]["physical_300k_beta_eq_paracrine_preserved_passes_clinical_glycemic_safety_targets"] is True

    # 5. Layer 5 & Wet-Lab Contract: Reject non-existent / Family-B product_variant_id and wrong estimand_scale, C1_C1 underscore normalization
    from scripts import ingest_wetlab_results as wetlab
    from scripts import run_layer5_personalized_patient_compiler as l5_mod

    with pytest.raises(ValueError, match="Unknown product variant"):
        wetlab.run_calibration({"batch_id": "B-BAD-PROD", "product_variant_id": "ISLET1_DOES_NOT_EXIST", "obs_stage7_yield_ratio": 0.95, "obs_stage7_yield_se": 0.02}, is_dry_run=True)
    with pytest.raises(ValueError, match="AUDIT-REF"):
        wetlab.run_calibration({"batch_id": "B-FAM-B", "product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE", "obs_stage7_yield_ratio": 0.95, "obs_stage7_yield_se": 0.02}, is_dry_run=True)
    with pytest.raises(ValueError, match="AUDIT-REF"):
        wetlab.run_calibration({"batch_id": "B-BAD-SCALE", "estimand_scale": "percentage_of_unedited_WT_control", "obs_stage7_yield_ratio": 0.95, "obs_stage7_yield_se": 0.02}, is_dry_run=True)

    pt_c2 = {
        "patient_id": "PT-TEST-C2",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "ifih1_rs1990760_genotype": "C/T",
        "ptpn22_rs2476601_genotype": "G/G",
        "recipient_hla_e_genotype": "*01:03/*01:03",
        "tyk2_rs34536443_genotype": "G/G",
        "recipient_kir_ligand_group": "C2/C2",
        "recipient_bmi_kg_m2": 24.0,
        "carries_rare_snp_in_standard_grna": False,
    }
    reg_b = l5_mod.compile_personalized_patient_blueprint(pt_c2)
    reg_b_mismatched = dict(reg_b)
    reg_b_mismatched["locus1_construct_prescribed"] = "PLAN_B_SELECTIVE_HLAA_HLAB_KO_RETAINED_HLAC_SUBLINE [C1_C1]"
    reg_b_mismatched["product_definition"] = dict(reg_b["product_definition"])
    reg_b_mismatched["product_definition"]["alleles"] = dict(reg_b["product_definition"]["alleles"])
    reg_b_mismatched["product_definition"]["alleles"]["HLA-C"] = "RETENTION_REQUIRED_KIR_MATCHED_C1_C1"
    assert l5_mod.evaluate_regimen_survival_probability(pt_c2, reg_b_mismatched) < 0.50


def test_sonnet55_and_astra_specialists_dec034_remediations() -> None:
    """Verify DEC-034 (Audit Step 5 + 3 Parallel Audit Specialists) remediations across Layers 1-5."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: B-prime K=5, B=1 counterfactual, Mode B 29 genotypes / 87 founders, and deep-copy graveyard cache
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        b_proj = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]["graveyard_005_and_010_counterfactual_sensitivity_sonnet55_r2"]
        assert b_proj["design_b_prime_k5_b1_counterfactual_projection_sonnet55_r3"]["projected_composite_utility_at_k5_b1"] == 0.0782
        assert b_proj["design_b_prime_k5_b1_counterfactual_projection_sonnet55_r3"]["projected_pediatric_utility_at_k5_b1"] == 0.0771
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    assert "mode_b_genotype_and_clone_sharing_note_astra_s1" in hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]
    from scripts import rank_islet1_stacks as ranker
    gy_copy = ranker._load_graveyard_entries_by_id()
    gy_copy["GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER"]["severity"] = "WARN"
    assert ranker._load_graveyard_entries_by_id()["GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER"]["severity"] == "Fatal_Gate2_and_Gate4_Failure"

    # 2. Layer 2: L1=1888..1947 vs mature B2M 1948..2007, matched 4128bp aliases, and AVI legacy alias
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    a_star = ch["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]
    coord_tbl = a_star["sense_hla_e_sct_coordinate_table_in_4143bp_and_4128bp_donors_astra_r3"]
    assert coord_tbl["L1_linker_20aa_60bp_in_4143bp_donor_1based"] == "1888..1947"
    assert coord_tbl["mature_B2M_first_20aa_in_4143bp_donor_1based"] == "1948..2007"
    assert len(a_star["arm1b_spa_s3_15aal1_matched_linker_control_4128bp_donor_sequence"]) == 4128
    avi_prov = ch["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]["unconditional_3locus_recut_immunity_pam_shield_bundle_sonnet55_r2"]["alphagenome_avi_live_verification_2026_10_27"]["provenance_and_evidence_class_muse_r2"]
    assert avi_prov["legacy_evidence_class_alias"] == "SUPERSEDED_BY_NINE_SCOPED_ABX_SNVS_CROSS_CHECKED_AGAINST_LOCAL_3118_ARCHIVE"

    # 3. Layer 3: Like-for-like 8-copy Ile1 chord mean 19.41 A, 0.25A r_obs=2.55A jump 34.0 A, and PAE terminal prose 21.27 / 25.45 A
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert geo["covalent_anchor_nearest_non_peptide_neighbors_3cdg_copy1_sonnet55_r2"]["like_for_like_8_copies_with_1mhe_ile1a_n_chord_mean_angstrom"] == 19.41
    assert geo["excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2"]["grid_spacing_discretization_sensitivity_3cdg_copy1_r_obs_2p5A_astra_r3"]["spacing_0p25A_r_obs_2p55A_percolation_jump_angstrom"] == 34.0
    c3 = next(c for c in l3["complexes"] if "SERPINB9" in c["complex_id"])
    pae_note = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]["alphafold_pae_uncertainty_summary_claude_v80_065"]["pae_row_aligned_col_scored_convention_astra_r3"]["convention_note"]
    assert "21.27 A" in pae_note and "25.45 A" in pae_note and "27.95" not in pae_note

    # 4. Layer 4: Strict sub-hourly 15.0341 mmHg at 1.096x OCR, two-scale OCR-propagated BVP [15.0006, 3.0808] / [14.431, 1.969], r_clamp=0.0161938861, and DOP853 3.500446
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    pa = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    dx = pa["design_x_100um_lamina_transport_qualification_v80_018"]
    ocr_sens = dx["fixed_inventory_two_scale_ocr_headroom_sensitivity_astra_r1"]
    assert ocr_sens["ocr_multiplier_1p096x_strict_sub_hourly_cluster_surface_po2_mmhg"] == 15.0341
    assert ocr_sens["ocr_multiplier_1p097x_design_x_two_scale_propagated_120um_surface_and_core_mmhg"] == [15.0006, 3.0808]
    assert ocr_sens["ocr_multiplier_1p097x_design_a_19p92mmhg_slab_with_radial_1p097x_120um_surface_and_core_mmhg"] == [14.431, 1.969]
    assert "0.0161938861" in dx["ginstling_brounshtein_nz201_localized_catalase_sonnet55_r2"]["early_hour_burst_and_pdms_barrier_note"]
    assert "matched_a_vs_x_geometry_and_kinetics_comparison_sonnet55_r3" in dx["blueprint_a_outer_band_115p5um_vs_100um_cluster_cap_reconciliation_sonnet55_r2"]
    pb = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]
    bridge_d34 = pb["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]
    assert bridge_d34["rk4_step_halving_refinement_audit_muse_r2"]["dop853_continuous_reference_nadir_at_T5p17_mM"] == 3.500446
    from scripts import run_layer4_tissue_and_metabolic_twin as l4_mod
    with pytest.raises(ValueError):
        l4_mod._validate_lamina_half_width_fraction(complex(0.5, 0.0))

    # 5. Layer 5 & Wet-Lab: Threshold table [4982, 3583, 2367, 1148], C2/C1 normalization, C0702 C1 matching, NONE_* cloak sentinel, WETLAB-01/02
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    t_rows = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["plan_b_kir_nk_threshold_sensitivity_sonnet55_r2"]
    assert [r["plan_b_routed_patients_of_10000"] for r in t_rows] == [4982, 3583, 2367, 1148]
    from scripts import ingest_wetlab_results as wetlab
    from scripts import run_layer5_personalized_patient_compiler as l5_mod
    pt_het = {
        "patient_id": "PT-HET",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.54,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 45.0,
        "sirpa_genotype": "v1/v1",
        "ifih1_rs1990760_genotype": None,
        "ptpn22_rs2476601_genotype": None,
        "recipient_hla_e_genotype": None,
        "tyk2_rs34536443_genotype": "G/G",
        "recipient_kir_ligand_group": "C2/C1",
        "recipient_bmi_kg_m2": 24.0,
        "carries_rare_snp_in_standard_grna": False,
    }
    reg_het = l5_mod.compile_personalized_patient_blueprint(pt_het)
    assert l5_mod.evaluate_regimen_survival_probability(pt_het, reg_het) > 0.95
    reg_c1_mis = dict(reg_het, product_definition=dict(reg_het["product_definition"], alleles=dict(reg_het["product_definition"]["alleles"], **{"HLA-C": "RETAINED_HLA_C0702_HOMOZYGOUS_SP3C"})))
    assert l5_mod.evaluate_regimen_survival_probability(pt_het, reg_c1_mis) < 0.50
    reg_no_cloak = dict(reg_het, product_definition=dict(reg_het["product_definition"], cassette_topology=dict(reg_het["product_definition"]["cassette_topology"], CD47_HLAE="NONE_B_PRIME_CLOAK_NOT_MATERIALIZED_IN_PRODUCT_DEFINITIONS")))
    assert l5_mod.evaluate_regimen_survival_probability(pt_het, reg_no_cloak) < 0.80

    emp1 = wetlab.run_calibration({"batch_id": "B-EMP-1", "obs_stage7_yield_ratio": 0.95, "obs_stage7_yield_se": 0.02}, is_dry_run=False)
    emp1_synth = dict(emp1, contains_synthetic_records=True)
    emp1_synth["calibration_state_sha256"] = wetlab._calibration_state_sha256(emp1_synth)
    with pytest.raises(ValueError, match="WETLAB-01"):
        wetlab.run_calibration({"batch_id": "B-EMP-2", "obs_stage7_yield_ratio": 0.96, "obs_stage7_yield_se": 0.02}, is_dry_run=False, prior_posterior_report=emp1_synth)


def test_sol61_specialists_dec035_remediations(monkeypatch) -> None:
    """Verify DEC-035 (3 Parallel GPT-6.1 Sol Max Specialists) remediations across Layers 1-5."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Product contract full SHA-256 verification and negative control rejection
    from scripts import product_contract as pc
    a_prod = pc.get_product_definition("ISLET1_A4_ARM1B")
    assert a_prod["locus1_donor_variants_metadata"]["canonical_arm1b_3974bp"]["donor_sha256"] == "f74728767987c18bdcb8c9d49773ee8dbe10e6006e57aa622638e4fb07062b0b"
    assert a_prod["locus1_donor_variants_metadata"]["embodiment_a_star_4143bp"]["donor_sha256"] == "0995f8e35791d219dcc905df4bd50145c9b2a43ddb84efdb502b41dc1f2ae6a8"
    bad_prod = json.loads(json.dumps(a_prod))
    bad_prod["locus1_donor_variants_metadata"]["embodiment_a_star_4143bp"]["donor_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="Invalid or mismatched donor_sha256"):
        pc.validate_product_definition(bad_prod)

    # 2. Layer 1: cells_plated_n effective-opportunity disambiguation and condition_number_22col_raw == 11.2549
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        c_tiers = l1["whole_product_clonogenic_safety_model"]["clonogenic_validation_tier_distinction"]
        assert c_tiers["master_cell_bank_pooled_validation"]["cells_plated_n_is_effective_clonogenic_opportunities_legacy_alias"] is True
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    mode_a = hip_p["hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1"]["operating_modes_and_condition_number_reconciliation_argon_r2"]["mode_a_primary_binary_bridge_Z3_plus1"]
    assert mode_a["condition_number_22col_raw"] == 11.2549

    # 3. Layer 2 & 3: 10 shared codons note, 19.4068875 4-dec mean, common Ile1/Ile1A voxel ranges, citation cache fixes
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    assert "10 shared codons" in ch["Arm_1B_Sense_Readthrough_Terminator_Shield_challenger_astra_r1"]["sense_hla_e_sct_coordinate_table_in_4143bp_and_4128bp_donors_astra_r3"]["l1_codon_diversification_note"]
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c1 = next(c for c in l3["complexes"] if "HLAE" in c["complex_id"])
    geo = c1["single_chain_trimer_linker_physics"]["excluded_volume_surface_geodesic_l1_audit_sonnet55_r1"]
    assert geo["covalent_anchor_nearest_non_peptide_neighbors_3cdg_copy1_sonnet55_r2"]["like_for_like_8_copies_with_1mhe_ile1a_n_chord_mean_from_4dec_rounded_inputs_angstrom"] == 19.4068875
    vox8 = geo["excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2"]["live_3d_voxel_dijkstra_8_copy_summary_at_0p5A_grid_astra_r3"]
    assert vox8["r_obs_2p5A"]["range_open_copies_common_ile1_ile1a_anchor_angstrom"] == [22.4, 39.8]
    cit = json.loads((ROOT / "data/raw/doc_citation_metadata_cache.json").read_text(encoding="utf-8"))["pmids"]
    assert cit["37703830"]["journal"] == "Immunity"
    assert cit["18585785"]["pubdate_year"] == "2008" and cit["18585785"]["journal"] == "Mol Immunol"

    # 4. Layer 4: Primary bvp_r60 propagates demand_multiplier_extra into surface_gsis_15mmhg_gate_passed
    from scripts import run_layer4_tissue_and_metabolic_twin as l4_mod
    r_110 = l4_mod._run_single_pouch_pde_regime(
        "TEST_DX_OCR_1P10X", 0.45, "central", 68.0, True, 0.97,
        lamina_half_width_fraction=100.0 / 350.0, demand_spatial_mode="outer_annulus",
        lamina_d_o2_scale=0.42 ** 1.5, catalase_spatial_mode="central_lamina_masked",
        catalase_k1_override_s_inv=7.56, nz_override=201, demand_multiplier_extra=1.10, _assess_numerical_gate=False
    )
    assert r_110["passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only"] is False
    assert r_110["wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001"]["two_zone_fv_gel_only_supply_surface_po2_mmhg"] == 14.90

    # 5. Layer 5: Dynamic prior coupling in plan_b_kir_nk_threshold_sensitivity_sonnet55_r2 and secondary_safe_harbor_cloak="KO"
    from scripts import rank_islet1_stacks as ranker
    from scripts import run_layer5_personalized_patient_compiler as l5_mod
    monkeypatch.setitem(ranker.ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["locus_dsb_penalty_rate_lambda_K"], "point_estimate", 0.25)
    monkeypatch.setitem(ranker.ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["base_edit_nickase_penalty_factor_gamma_B"], "point_estimate", 0.92)
    sim_alt = l5_mod.simulate_10000_patient_global_population_trial()
    t_disc_alt = [r["complexity_discounted_personalized_utility_pct"] for r in sim_alt["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["plan_b_kir_nk_threshold_sensitivity_sonnet55_r2"]]
    assert t_disc_alt == [83.7095, 86.6191, 89.1551, 91.5311]


def test_muse_wonder_r3_dec036_high_impact_biological_upgrades() -> None:
    """Verify DEC-036 (Audit Full Round 3 Step 6 High-Impact Biological & Clinical Failure-Mode Red-Team) across Layers 1-5 & Stage-0 Decision Tree."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1: W_BF Basal-Floor satellite, B'-SCT promotion, pass-aware tuple, DX-S correlation pricing, and TYK2xISRE interaction estimands
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        pf = l1["champion_challenger_portfolio"]
        sub8 = pf["blueprint_x_stage0_molecular_challenger_spec"]["stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1"]
        assert sub8["locus3_w_bf_basal_floor_2arm_satellite_muse_r3"]["n_arms"] == 2
        assert sub8["locus3_w_bf_basal_floor_2arm_satellite_muse_r3"]["total_clones"] == 6
        triad = pf["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]
        assert "b2m_null_orphan_hc_erad_status_muse_r3" in triad["A_Champion_Islet1"]
        b_node = triad["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]
        assert b_node["b_prime_4digit_hlac_rearming_and_b_prime_sct_promotion_muse_r3"]["primary_product_variant_ids"] == [
            "ISLET1_B_SCT_C1_CIITA",
            "ISLET1_B_SCT_C2_CIITA",
            "ISLET1_B_SCT_C1C2_CIITA",
        ]
        assert b_node["b_prime_4digit_hlac_rearming_and_b_prime_sct_promotion_muse_r3"]["pass_aware_editing_complexity_tuple_muse_r3"]["B_prime_SCT_primary"]["P_passes"] == 2
        assert "dxs_single_locus_shield_plus_switch_correlation_pricing_muse_r3" in triad["X_Locked_Number1_Evolutionary_Contender"]
    hip_p = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    tyk2_int = hip_p["alias_structure"]["co_primary_tyk2_inducible_promoter_interaction_estimands_muse_r3"]
    assert tyk2_int["mode1_23col_half_effect_se_pp_at_sigma_run_5p4314"] == {
        "X1:X4": 1.7797,
        "X3:X4": 2.6407,
        "mean_across_15_2fi_pairs": 2.7669,
        "arithmetic_mean_across_15_2fi_pairs": 2.5924,
        "arithmetic_mean_across_15_2fi_pairs_exact": 2.5915,
        "mean_across_6_main_effects": 1.2912,
        "arithmetic_mean_across_6_main_effects": 1.2910,
    }
    assert tyk2_int["co_primary_interactions"] == ["X1:X4", "X3:X4"]

    # 2. Layer 2: L3-BF-B (4,078 bp) & L3-BF-C (4,160 bp) Basal-Floor donors, Hour-0 kinetic race fields, F2 ER heavy-chain flux, and F7 TYK2 epistasis
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    bf_node = ch["Locus3_BasalFloor_ScarFree_SERPINB9_PDL1_challenger_muse_r3"]
    assert bf_node["l3_bfb_donor_length_bp"] == 4078
    assert bf_node["l3_bfc_dual_basal_floor_donor_length_bp"] == 4160
    p_plane = l2["stage0_to_stage7_perturbation_plane"]
    p0 = p_plane["perturbation_profiles"][0]
    assert "hour0_kinetic_race_status" in p0
    assert "er_heavy_chain_flux_and_upr_hazard_muse_r3" in p_plane
    assert "tyk2_inducible_promoter_epistasis_prior_band_muse_r3" in p_plane

    # 3. Layer 4: Chronic zero-source & mineralization sweep (F4), SC-islet immaturity paracrine sweep (F5), and ER load coupling note (F2)
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    f4_node = dx["joint_adverse_corner_stress_audit_muse_r1"]["chronic_zero_source_and_mineralization_extension_muse_r3"]
    assert "cao2_exhaustion_day14_to_day90_regime" in f4_node
    part_b = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]
    assert "er_load_viability_discount_coupling_note_muse_r3" in part_b
    f5_node = part_b["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]
    assert len(f5_node["corners"]) == 6

    # 4. Layer 5 & Stage-0 Decision Tree: B'-SCT product definitions, 4-digit HLA-C alloreactivity penalty, and deterministic S1..S9 precedence
    from scripts import ingest_wetlab_results as wetlab
    from scripts import product_contract as pc
    from scripts import run_layer5_personalized_patient_compiler as l5_mod

    for vid in (
        "ISLET1_A4_X_MOL_RNLS_WT",
        "ISLET1_A4_DXS_ARM1C_L4CD47",
        "ISLET1_A4_L3_SCARFREE_3799",
        "ISLET1_A4_L3_BFB_BASALFLOOR",
        "ISLET1_B_SCT_C1_CIITA",
        "ISLET1_B_SCT_C2_CIITA",
    ):
        pdef = pc.get_product_definition(vid)
        pc.validate_product_definition(pdef)

    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    mfg = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]
    assert mfg["b_prime_4digit_hlac_alloreactivity_and_b_prime_sct_pivot_audit_muse_r3"]["endogenous_hlac_retention_product_demotion"]["status"] == "EXPLORATORY_REQUIRES_4DIGIT_BANK"
    assert mfg["hour0_kinetic_race_and_basal_floor_audit_muse_r3"]["basal_floor_challenger_product_variant_id"] == "ISLET1_A4_L3_BFB_BASALFLOOR"

    pt_allo = {
        "patient_id": "PT-4DIGIT-ALLO",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.2,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C1",
        "hlac_4digit_relation": "KIR_MATCH_FOREIGN_4DIGIT_ALLO",
    }
    reg_b = l5_mod.compile_personalized_patient_blueprint(pt_allo)
    assert reg_b["b_prime_4digit_hlac_and_b_prime_sct_routing_audit_muse_r3"]["canonical_track2_kir_subbank_product_variant_id"] == "ISLET1_B_SCT_C1_CIITA"
    assert l5_mod.evaluate_regimen_survival_probability(pt_allo, reg_b) < 0.75

    # Verify deterministic Stage-0 decision tree precedence S1..S9
    dt_res = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_gate0c_hour0_serpinb9_uM": 0.12,
            "obs_gate0c_2h_lysis_pct": 14.5,
            "obs_gate0b_hlac_cd8_elispot_sfw": 150.0,
            "obs_gate0b_hlac_sct_molecules_per_cell": 25000.0,
        }
    )
    assert dt_res["strict_precedence_order"][0] == "1_SAFETY_KILL_S9"
    assert len(dt_res["switch_decisions"]) == 11
    assert dt_res["switch_decisions"]["S1_Locus3_Kinetic_and_Topology_Gate_0C"] == "PROMOTE_ISLET1_A4_L3_BFB_BASALFLOOR"
    assert dt_res["switch_decisions"]["S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B"] == "PROMOTE_ISLET1_B_SCT_C1_OR_C2_CIITA"


def test_argon_r4_dec037_biological_redteam_upgrades() -> None:
    """Verify DEC-037 (Audit / Gemini 3.1 Pro Round 4 Step 1 5-Agent Biological & Biotransport Red-Team Swarm) across Layers 1-5 & Stage-0 Contract."""
    ROOT = Path(__file__).resolve().parent.parent

    # 1. Layer 1 & Layer 3: GAVDPLLAL (PDB:1EFX) & QYDDAVYKL (PDB:1IM9) human self-peptide HLA-C SCTs + Tyr84Ala F-pocket lateral channel + B'-SCT counterfactual projection
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        b_node = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]["B_Stratified_KIR_SubBank_PlanB_and_Bprime"]["b_prime_4digit_hlac_rearming_and_b_prime_sct_promotion_muse_r3"]
        assert "GAVDPLLAL" in b_node["tethered_human_self_peptides_and_f_pocket_clearance_argon_r4"]["C1_HLA_C0304_self_peptide"]
        assert "QYDDAVYKL" in b_node["tethered_human_self_peptides_and_f_pocket_clearance_argon_r4"]["C2_HLA_C0401_self_peptide"]
        assert b_node["design_b_prime_sct_k4_b0_counterfactual_projection_argon_r4"]["engineering_complexity_discount"] == 1.0

    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    hlac_audit = l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]
    assert hlac_audit["c1_subbank_hlac0304_kir2dl2_1EFX"]["self_peptide_9mer"] == "GAVDPLLAL"
    assert hlac_audit["c2_subbank_hlac0401_kir2dl1_1IM9"]["self_peptide_9mer"] == "QYDDAVYKL"
    assert hlac_audit["f_pocket_lateral_channel_mutation_required_for_hlac_sct"]["engineered_groove_mutation"] == "Tyr84Ala (Y84A)"

    # 2. Layer 2: L3-BF-B hybrid basal+inducible promoter (Syn_3xISRE_2xNFkB_hPGK1_Sp1_basalfloor_230) & L3-BF-C repeat-free dual-promoter
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    bf_node = ch["Locus3_BasalFloor_ScarFree_SERPINB9_PDL1_challenger_muse_r3"]
    assert bf_node["l3_bfb_donor_length_bp"] == 4078
    assert bf_node["l3_bfb_5top_oligopyrimidine_tract_present"] is False
    assert bf_node["l3_bfc_dual_basal_floor_donor_length_bp"] == 4160
    assert bf_node["l3_bfc_single_promoter_bicistronic_3915bp_donor_length_bp"] == 3915
    er_node = l2["stage0_to_stage7_perturbation_plane"]["er_heavy_chain_flux_and_upr_hazard_muse_r3"]
    assert "zero free soluble B2M in trans" in er_node["cis_tethered_hlae_sct_cannot_chaperone_orphan_hc_in_trans_note_argon_r4"]

    # 3. Layer 4: Physical f_source=0.34 chronic zero-source audit (D_O2_gel >= 1.50e-5 cm^2/s) & 12-state bidirectional glucagon/ISF-glucose lag ODE
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    f4_phys = dx["joint_adverse_corner_stress_audit_muse_r1"]["chronic_zero_source_and_mineralization_extension_muse_r3"]["physical_outer_annulus_masked_chronic_zero_source_audit_argon_r4"]
    assert f4_phys["outer_annulus_masked_s_chem_0_pb40_d_o2_high_1p50e_5"]["delta_z_56p0um"]["100um_dxb_surf_core_mmhg"][0] >= 15.0
    part_b = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]
    bidir = part_b["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]["bidirectional_glucagon_isf_glucose_lag_and_portal_vs_subcutaneous_partition_argon_r4"]
    assert bidir["omental_portal_draining_site_partition_f_sys_ins_1p00_f_hep_gcg_1p00"]["T_ins_4p00min_t50_3p36min_tightened_bidirectional_gate"]["clears_3p50mM_at_v1p00_local_I1"] is True

    # 4. Layer 5 & Wet-Lab Contract: B'-SCT product variants pass evaluate_regimen_survival_probability and S9 safety kill triggers hard stop
    from scripts import ingest_wetlab_results as wetlab
    from scripts import product_contract as pc
    from scripts import run_layer5_personalized_patient_compiler as l5_mod

    for vid in (
        "ISLET1_B_SCT_C1_CIITA_TYK2_PE",
        "ISLET1_B_SCT_C2_CIITA_TYK2_PE",
        "ISLET1_B_SCT_C1C2_CIITA",
        "ISLET1_B_SCT_C1C2_CIITA_TYK2_PE",
    ):
        pc.validate_product_definition(pc.get_product_definition(vid))

    pt_c1 = {
        "patient_id": "PT-B-SCT-C1",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C1",
        "hlac_4digit_relation": "KIR_MATCH_FOREIGN_4DIGIT_ALLO",
    }
    reg_sct = {
        "selected_fda_deployment_track": "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE",
        "product_variant_id": "ISLET1_B_SCT_C1_CIITA",
        "product_definition": pc.get_product_definition("ISLET1_B_SCT_C1_CIITA"),
        "locus1_construct_prescribed": "OI-L2-DONOR-1B_B2M_EEF1A1_HLAE_CD47IGV + CIITA_EXON5_B_SCT_C0304_GAVDPLLAL",
        "plan_architecture_prescribed": "Plan B-SCT (B2M::Arm1B + CIITA::HLA-C*03:04[GAVDPLLAL] SCT)",
        "include_tyk2_p1104a_rs34536443_shield": False,
        "compiled_grna_panel": [
            "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "OI-GRNA-02_CIITA_EXON5_ASCAS12A",
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        ],
    }
    assert abs(l5_mod.evaluate_regimen_survival_probability(pt_c1, reg_sct) - 0.9017) < 1e-6
    assert abs(l5_mod.evaluate_regimen_survival_probability({**pt_c1, "b_prime_sct_recipient_crossmatch_and_mlc_cleared": True}, reg_sct) - 0.9817) < 1e-6

    dt_kill = wetlab.evaluate_stage0_decision_tree({"obs_icasp9_log10_kill": 3.5})
    assert dt_kill["safety_kill_s9_triggered"] is True
    assert dt_kill["resolved_architecture_assignments"]["variant_safety_status"] == "HARD_SAFETY_STOP_REJECTED"


def test_sol61_ultracode_r4_dec038_biological_redteam_upgrades() -> None:
    """Verify DEC-038 (GPT AUDIT-REF Audit Round 4 Step 2 Biological & Clinical Failure-Mode Red-Team) across Layers 1-5 & Wet-Lab Contract."""
    ROOT = Path(__file__).resolve().parent.parent
    from scripts import ingest_wetlab_results as wetlab
    from scripts import product_contract as pc
    from scripts import run_layer5_personalized_patient_compiler as l5_mod
    from scripts import wetlab_contract as wc

    # 1. R4-01 (SAF-01): Fail-closed Stage-0 safety panel & uncertainty-aware 95% bounds + ingestion normalization
    s0_mean_keys = {t[1] for t in wc.STAGE0_DECISION_TREE_ENDPOINT_SPECS}
    assert "obs_tp53_ppm1d_bcor_vaf" in s0_mean_keys
    assert "obs_castseq_translocation_freq" in s0_mean_keys
    norm_s9, _ = wc.normalize_measurements(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_icasp9_log10_kill_se": 0.1,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_serpinb9_hour0_uM": 0.35,
        }
    )
    assert norm_s9["obs_tp53_ppm1d_bcor_vaf"] == 0.0

    # Empty or incomplete safety observations MUST fail closed on HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL
    dt_empty = wetlab.evaluate_stage0_decision_tree({})
    assert dt_empty["mandatory_s9_safety_panel_complete"] is False
    assert dt_empty["safety_kill_s9_status"] == "HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL"
    assert dt_empty["resolved_architecture_assignments"]["variant_safety_status"] == "HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL"
    assert dt_empty["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_INCOMPLETE_SAFETY_PANEL"

    # Complete clean safety panel clears S9
    dt_clean = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_icasp9_log10_kill_se": 0.1,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
            "obs_castseq_translocation_freq": 0.0,
            "obs_castseq_translocation_se": 1e-6,
        }
    )
    assert dt_clean["mandatory_s9_safety_panel_complete"] is True
    assert dt_clean["safety_kill_s9_triggered"] is False
    assert dt_clean["resolved_architecture_assignments"]["variant_safety_status"] == "CLEARED"
    assert dt_clean["resolved_architecture_assignments"]["lot_maturation_status"] == "READY_FOR_ENCAPSULATION"

    # High SE pushing one-sided 95% bound across safety threshold triggers S9 hard stop
    dt_se_fail = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.10,
            "obs_icasp9_log10_kill_se": 0.10,  # lower95 = 4.10 - 0.1645 = 3.9355 < 4.0
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
        }
    )
    assert dt_se_fail["safety_kill_s9_triggered"] is True
    assert dt_se_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HARD_SAFETY_STOP_REJECTED"

    # S3 DX-S promotion requires observed obs_dxs_locus4_joint_retention_fraction >= 0.9999
    dt_s3_unqualified = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_cd47_hour0_molecules": 40000.0,
        }
    )
    assert dt_s3_unqualified["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] == "INACTIVE"
    dt_s3_missing_counts = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_cd47_hour0_molecules": 40000.0,
            "obs_dxs_locus4_joint_retention_fraction": 0.99995,
            "obs_dxs_cd47low_icasp9low_fraction": 0.0,
            "obs_dxs_ap1903_cd47high_escape_fraction": 0.0,
        }
    )
    assert dt_s3_missing_counts["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] == "HOLD_DXS_RARE_EVENT_REQUIRES_EXPLICIT_SE_OR_CELL_COUNT"
    dt_s3_qualified = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_cd47_hour0_molecules": 40000.0,
            "obs_dxs_locus4_joint_retention_fraction": 0.99995,
            "obs_dxs_cd47low_icasp9low_fraction": 0.0,
            "obs_dxs_ap1903_cd47high_escape_fraction": 0.0,
            "obs_dxs_screened_cell_count": 300000000,
            "obs_dxs_clonogenic_fraction": 1.0,
            "obs_dxs_spike_in_detection_efficiency": 1.0,
        }
    )
    assert dt_s3_qualified["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] == "PROMOTE_ISLET1_A4_DXS_ARM1C_L4CD47"

    # 2. R4-02, R4-06, R4-10, R4-18: B'-SCT anti-HLA-C DSA, KIR2DS1 non-activation gate, and 1EFX/1IM9 Tyr84-Asp183 contact audit
    pt_dsa = {
        "patient_id": "PT-DSA-POS",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C2/C2",
        "hlac_4digit_relation": "KIR_MATCH_FOREIGN_4DIGIT_ALLO",
        "recipient_anti_hla_c_dsa": True,
    }
    reg_c2_sct = {
        "selected_fda_deployment_track": "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE",
        "product_variant_id": "ISLET1_B_SCT_C2_CIITA",
        "product_definition": pc.get_product_definition("ISLET1_B_SCT_C2_CIITA"),
        "locus1_construct_prescribed": "OI-L2-DONOR-1B_B2M_EEF1A1_HLAE_CD47IGV + CIITA_EXON5_B_SCT_C0401_QYDDAVYKL",
        "plan_architecture_prescribed": "Plan B-SCT (B2M::Arm1B + CIITA::HLA-C*04:01[QYDDAVYKL] SCT)",
        "include_tyk2_p1104a_rs34536443_shield": False,
        "compiled_grna_panel": [
            "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "OI-GRNA-02_CIITA_EXON5_ASCAS12A",
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        ],
    }
    assert l5_mod.evaluate_regimen_survival_probability(pt_dsa, reg_c2_sct) < 0.60
    bp_dsa = l5_mod.compile_personalized_patient_blueprint(pt_dsa, track2_subbank_mode="B_PRIME_SCT")
    assert bp_dsa["selected_product_id"] == "ISLET1_B_SCT_C2_CIITA"
    assert bp_dsa["b_prime_4digit_hlac_and_b_prime_sct_routing_audit_muse_r3"]["track2_hlac_4digit_match_status"] == "HOLD_INCOMPATIBLE_ANTI_HLAC_DSA_CONTRAINDICATED_REQUIRES_CROSSMATCH_OR_TRACK3"

    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    y84_audit = l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]
    assert abs(y84_audit["pdb_1efx_hlac0304_gavdpllal_kir2dl2_native_geometry"]["tyr84_oh_to_kir2dl2_asp183_od2_distance_angstrom"] - 3.5584) < 1e-4
    assert abs(y84_audit["pdb_1im9_hlac0401_qyddavykl_kir2dl1_native_geometry"]["tyr84_ce1_to_kir2dl1_asp183_od2_distance_angstrom"] - 3.2425) < 1e-4

    # 3. R4-03, R4-04, R4-16, R4-17: Layer 4 1,500 IEQ/cm2 low-load envelope, CaO2->CaCO3 expansion (+0.1336 cm3), 360k IEQ immaturity counterexamples & 13-state ODE
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    dx = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]
    f4_phys = dx["joint_adverse_corner_stress_audit_muse_r1"]["chronic_zero_source_and_mineralization_extension_muse_r3"]["physical_outer_annulus_masked_chronic_zero_source_audit_argon_r4"]
    lowload = f4_phys["joint_chronic_source_free_areal_loading_capsule_ocr_envelope_sol61_r4"]["scenarios_at_d_o2_1p50e_5_outer_annulus_s_chem_0"]
    assert lowload["load_2500_pb40_cap14um_phys_42um_eq_ocr1p096"]["clears_15mmhg_surface_and_5mmhg_core_for_100um_dxb"] is False
    assert lowload["load_1500_dxb_lowload_pb35_cap15um_phys_45um_eq_ocr1p097"]["clears_15mmhg_surface_and_5mmhg_core_for_100um_dxb"] is True
    assert abs(f4_phys["phase_resolved_cao2_to_caco3_net_solid_expansion_audit_sol61_r4"]["net_solid_volume_expansion_cao2_to_caco3_cm3"] - 0.1336) < 1e-3

    part_b = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]
    bidir = part_b["parameter_literature_provenance_table_v80_091"]["total_endocrine_ieq_to_beta_equivalent_ieq_dose_bridge_sol61_r2"]["joint_actual_dose_v060_and_erlang_lag_qualification_astra_r2"]["bidirectional_glucagon_isf_glucose_lag_and_portal_vs_subcutaneous_partition_argon_r4"]
    sol_ode = bidir["joint_360k_immaturity_counterexample_and_13state_hepatic_peripheral_split_sol61_r4"]
    assert sol_ode["immaturity_360k_ieq_v0p432_counterexamples"]["single_defect_alpha_gain_0p50_v0p432_T_ins_5p17min_min_glucose_mM"] < 3.50
    assert sol_ode["immaturity_360k_ieq_v0p432_counterexamples"]["full_immaturity_corner_v0p432_T_ins_5p17min_min_glucose_mM"] < 3.00
    assert sol_ode["immaturity_360k_ieq_v0p432_counterexamples"]["full_immaturity_corner_v0p432_T_ins_2p00min_portal_ref_min_glucose_mM"] < 3.50

    # 4. R4-07, R4-08, R4-09, R4-15: Layer 1 & Layer 2 ER taxonomy, SERPINB9 serial-pulse ODE & IAPP challenger
    hip_proto = json.loads((ROOT / "data/protocols/hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    co_prim = hip_proto["alias_structure"]["co_primary_tyk2_inducible_promoter_interaction_estimands_muse_r3"]
    assert co_prim["isogenic_tyk2_p1104a_by_l3_bfb_hybrid_promoter_satellite_sol61_r4"]["design_id"] == "STAGE0_TYK2_GENOTYPE_BY_LOCUS3_PROMOTER_3X2_SATELLITE"

    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    er_tax = l2["stage0_to_stage7_perturbation_plane"]["er_heavy_chain_flux_and_upr_hazard_muse_r3"]["four_mechanism_er_taxonomy_and_tyk2_p1104a_context_sensitivity_sol61_r4"]
    assert er_tax["tyk2_p1104a_ifna_orphan_hc_flux_sensitivity_branches_vs_wt"]["human_leukocyte_wt_like_0_to_20pct_ifna_attenuation_branch_pmid_30578352_42268925"] == [0.8, 1.0]

    # 5. AUDIT-REF (DEC-039):
    # - Fail-closed string boolean rejection & separated anti-HLA-C DSA vs MLC clearance
    for bad_bool_key in ("recipient_anti_hla_c_dsa", "b_prime_sct_recipient_crossmatch_and_mlc_cleared"):
        try:
            l5_mod.evaluate_regimen_survival_probability(dict(pt_dsa, **{bad_bool_key: "false"}), reg_c2_sct)
            raise AssertionError(f"Expected ValueError for string boolean {bad_bool_key}")
        except ValueError:
            pass
    # MLC-only flag cannot clear anti-HLA-C DSA without exact serum crossmatch negative:
    pt_dsa_mlc_only = dict(pt_dsa, b_prime_sct_recipient_crossmatch_and_mlc_cleared=True)
    assert l5_mod.evaluate_regimen_survival_probability(pt_dsa_mlc_only, reg_c2_sct) < 0.65
    pt_dsa_serum_and_mlc = dict(
        pt_dsa,
        b_prime_sct_recipient_crossmatch_and_mlc_cleared=True,
        b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative=True,
    )
    assert l5_mod.evaluate_regimen_survival_probability(pt_dsa_serum_and_mlc, reg_c2_sct) > 0.90

    # - Stage-0 decision-tree JSON trigger aliases & CHOP/HC10/DX-S 95% confidence bounds
    dt_alias = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_hlae_hour0_molecules": 12000.0,
            "obs_regimeII_nkg2a_lysis_pct": 12.0,
            "obs_chop_coincidence_fold": 1.85,
            "obs_chop_coincidence_se": 0.15,
        }
    )
    assert dt_alias["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "PROMOTE_ISLET1_A4_X_MOL_RNLS_WT"
    assert dt_alias["safety_kill_s9_triggered"] is True  # chop_upper95 = 1.85 + 1.645 * 0.15 = 2.09675 >= 2.0

    # - Layer 3 31/22 open-chord obstacles, DT-1 vs DT-2, and CD47 TSP-1/shedding audit
    assert y84_audit["pdb_1efx_hlac0304_gavdpllal_kir2dl2_native_geometry"]["l1_open_chord_heavy_atom_obstacles_within_4p0A_excl_endpoints"]["y84a_all_heavy_atoms"] == 31
    assert y84_audit["pdb_1im9_hlac0401_qyddavykl_kir2dl1_native_geometry"]["l1_open_chord_heavy_atom_obstacles_within_4p0A_excl_endpoints"]["y84a_sidechain_plus_cbeta_atoms"] == 22
    assert "endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4" in l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]

    # - Layer 4 nz=201 monotone chronic PDE & ACC expansion
    nz201_audit = f4_phys["joint_chronic_source_free_areal_loading_capsule_ocr_envelope_sol61_r4"]["live_nz201_monotone_chronic_pde_and_postural_folded_occlusion_kimi31_r4"]
    assert nz201_audit["strict_ocr_monotonicity_at_2500_pb35_cap15um_nz201"]["monotonically_decreasing_in_ocr"] is True
    assert f4_phys["phase_resolved_cao2_to_caco3_net_solid_expansion_audit_sol61_r4"]["amorphous_acc_and_porous_mineral_phase_branches_kimi31_r4"]["hydrated_amorphous_calcium_carbonate_acc_1p75_g_cm3_pmid_38343991"]["solid_volume_cm3"] > 0.60

    # 6. AUDIT-REF (DEC-040):
    # - WET-R4-01..04: Scale-matched STAGE0_CUSTOM_DEFAULT_SE, conflicting alias rejection, S8b TYK2 epistasis, and lot_release_readiness_tier
    assert wc.STAGE0_CUSTOM_DEFAULT_SE["obs_dxs_cd47low_icasp9low_fraction"] == 1e-7
    assert wc.STAGE0_CUSTOM_DEFAULT_SE["obs_dxs_ap1903_cd47high_escape_fraction"] == 1e-9
    try:
        wc.normalize_measurements(
            {
                "obs_nkg2a_lysis_pct": 6.0,
                "obs_regimeII_nkg2a_lysis_pct": 12.0,
            }
        )
        raise AssertionError("Expected ValueError on conflicting co-supplied Stage-0 alias keys")
    except ValueError as exc:
        assert "Conflicting" in str(exc)
    dt_epistasis = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_beta_x1_x4_epistasis": -0.14,
            "product_variant_id": "ISLET1_A4_ARM1B_TYK2_PE",
        }
    )
    assert dt_epistasis["switch_decisions"]["S8b_TYK2_Inducible_Promoter_Epistasis_Gate_0E"] == "COUPLE_TYK2_P1104A_STRICTLY_WITH_ARM1X_AND_L3_BFB_CONSTITUTIVE_BASAL_FLOOR"
    assert dt_epistasis["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # - L5-S2-01..03: C2 B'-SCT un-cleared CD8+ MLC hold, strict hlac_4digit_relation validation, and 333.33 cm2 chronic low-load area
    try:
        l5_mod.compile_personalized_patient_blueprint(dict(pt_dsa, hlac_4digit_relation="UNKNOWN_RELATION"))
        raise AssertionError("Expected ValueError for invalid hlac_4digit_relation")
    except ValueError:
        pass
    pt_c2_no_mlc = {
        "patient_id": "PT-C2-NOMLC",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C2/C2",
        "hlac_4digit_relation": "SELF_4DIGIT_MATCH",
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared": False,
    }
    bp_c2_no_mlc = l5_mod.compile_personalized_patient_blueprint(pt_c2_no_mlc, track2_subbank_mode="B_PRIME_SCT")
    assert bp_c2_no_mlc["product_qualification_status"] == "HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_CLEARANCE"
    assert bp_c2_no_mlc["selected_product_qualification_status"] == "HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_CLEARANCE"
    assert bp_c2_no_mlc["lot_release_cleared"] is False
    assert abs(bp_c2_no_mlc["host_modifier_genotypes_and_dosing"]["required_pouch_exchange_area_cm2_at_chronic_lowload_1500_ieq_cm2"] - 333.33) < 0.05

    # - COV-01 & L2-D: Exact (X^T X)^-1 diagonals 347/3232 & 191/808 and canonical 3,996-bp Syn_5xISRE_minP comparator
    assert co_prim["mode1_23col_var_beta_over_sigma2"]["X1:X4"] == "347/3232 (0.107364)"
    assert co_prim["mode1_23col_var_beta_over_sigma2"]["X3:X4"] == "191/808 (0.236386)"
    syn_comp = co_prim["isogenic_tyk2_p1104a_by_l3_bfb_hybrid_promoter_satellite_sol61_r4"]["canonical_syn5xisre_comparator_construct_bindings_astra6_r4"]
    c_b2m = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    ch = c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    bf_node = ch["Locus3_BasalFloor_ScarFree_SERPINB9_PDL1_challenger_muse_r3"]
    assert syn_comp["matched_2afree_comparator"]["donor_length_bp"] == 3996
    assert syn_comp["matched_2afree_comparator"]["donor_sequence_sha256"] == bf_node["l3_syn5xisre_matched_2afree_comparator_3996bp_donor_sha256"] == "afae235efdf6bd013f6bb690ac53f7d0ebbf5a6590e441d225debb1bb42f85ba"
    assert syn_comp["canonical_scarfree_3799_comparator"]["donor_sequence_sha256"] == ch["Locus3_ScarFree_Dual_Promoter_SERPINB9_PDL1_challenger_astra_r1"]["assembled_3799bp_spaced_donor_sha256"] == "3d6bf68e2f0b70911cfea69255e4733db2862ce3baa2b92ae4206561e284bcf7"

    # - L4 live chronic BVP audit & splanchnic first-pass extraction sensitivity
    assert nz201_audit["computation_mode"] == "LIVE_NZ201_NEWTON_BVP_DERIVED_ASTRA6_R4"
    meal_audit = part_b["meal_mass_balance_and_distribution_volume_provenance_v80_008"]
    assert "splanchnic_first_pass_extraction_sensitivity_audit_astra6_r4" in meal_audit

    # 7. AUDIT-REF (DEC-041):
    # - S2-01: Exact Poisson upper95 bound at 62,958 duplex depth stops k >= 3 (VAF = 3/62958 = 4.765e-5 -> upper95 = 1.231e-4 > 1e-4)
    dt_vaf_k3 = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 3.0 / 62958.0,
            "obs_castseq_translocation_freq": 0.0,
        }
    )
    assert dt_vaf_k3["safety_kill_s9_triggered"] is True
    # - S2-01: DX-S zero-event row without explicit SE or >= 3e8 cell count does NOT promote DX-S
    dt_dxs_no_count = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_cd47_hour0_molecules": 40000.0,
            "obs_dxs_locus4_joint_retention_fraction": 1.0,
            "obs_dxs_cd47low_icasp9low_fraction": 0.0,
            "obs_dxs_ap1903_cd47high_escape_fraction": 0.0,
        }
    )
    assert dt_dxs_no_count["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] != "PROMOTE_ISLET1_A4_DXS_ARM1C_L4CD47"
    # - S2-05: C2/C1 normalized to C1/C2 and triggers C2 holds
    bp_c2c1 = l5_mod.compile_personalized_patient_blueprint(
        dict(pt_c2_no_mlc, recipient_kir_ligand_group="C2/C1"),
        track2_subbank_mode="B_PRIME_SCT",
    )
    assert bp_c2c1["product_qualification_status"] == "HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_CLEARANCE"
    assert "HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_CLEARANCE" in bp_c2c1["all_active_qualification_holds"]
    # - S2-04 & D-a4: P_b = 25 mmHg areal density sweep & 1EFX/1IM9 SHA-pinned live geometry
    assert "pb25_standoff100_areal_density_sweep_sonnet55_r4" in nz201_audit["postural_and_folded_pouch_contact_occlusion_sensitivity"]
    assert "dxb_lowload_1500_333cm2_full_336h_transient_pde_chain_sonnet55_r4" in nz201_audit["postural_and_folded_pouch_contact_occlusion_sensitivity"]
    assert l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]["pdb_1efx_hlac0304_gavdpllal_kir2dl2_native_geometry"]["lys146_nz_to_kir2dl2_asp183_od2_salt_bridge_angstrom"] == 2.4601
    assert "nkg2a_vs_nkg2c_2d_synaptic_selectivity_bound_sonnet55_r4" in l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]

    # 8. AUDIT-REF (DEC-042):
    # - R-01: Shallow depth (5,000x) at vaf=0.0 triggers S9 (2.99573/5000 = 5.99e-4 >= 1e-4);
    #         depth 40,000x (<62,958x) at vaf=0.0 holds on HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM;
    #         explicit SE=1e-6 with k=3 at 62,958x still triggers S9 via max(Wald, exact Poisson).
    dt_vaf_shallow = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_tp53_ppm1d_bcor_duplex_depth": 5000.0,
            "obs_castseq_translocation_freq": 0.0,
        }
    )
    assert dt_vaf_shallow["safety_kill_s9_triggered"] is True
    dt_vaf_middepth = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_tp53_ppm1d_bcor_duplex_depth": 40000.0,
            "obs_castseq_translocation_freq": 0.0,
        }
    )
    assert dt_vaf_middepth["safety_kill_s9_triggered"] is False
    assert dt_vaf_middepth["safety_kill_s9_status"] == "HOLD_TP53_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM"
    dt_vaf_k3_with_se = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 3.0 / 62958.0,
            "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
            "obs_castseq_translocation_freq": 0.0,
        }
    )
    assert dt_vaf_k3_with_se["safety_kill_s9_triggered"] is True

    # - R-02: S10 separates NKG2C+ adaptive NK peptide-selectivity failure from S2 Arm 1X promotion
    dt_s10_nkg2c = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_hlae_hour0_molecules": 18000.0,
            "obs_nkg2a_lysis_pct": 6.0,
            "obs_nkg2c_adaptive_nk_lysis_pct": 14.0,
        }
    )
    assert dt_s10_nkg2c["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "RETAIN_EMBODIMENT_A_STAR_ARM1B"
    assert dt_s10_nkg2c["switch_decisions"]["S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A"] == "STRATIFY_CMV_NKG2C_RECIPIENTS_AND_EVALUATE_SCT_SELECTIVITY_MUTANT_COUPON"
    assert dt_s10_nkg2c["resolved_architecture_assignments"]["nkg2c_adaptive_nk_selectivity_status"] == "HOLD_S10_NKG2C_SELECTIVITY_FAILED_REQUIRES_SCT_PEPTIDE_SELECTIVITY_COUPON_OR_CMV_STRATIFIED_LABEL"
    assert dt_s10_nkg2c["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # - R-03: 14-endpoint core panel returns PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED;
    #         full extended panel (with duplex depth >= 62958, NKG2C, KIR, phagocytosis, and chronic capsule) returns READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE.
    core_14_panel = {
        "obs_icasp9_log10_kill": 4.6,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_castseq_translocation_freq": 0.0,
        "obs_chop_coincidence_fold": 1.15,
        "obs_hc10_surface_fold_wt": 1.05,
        "obs_gate0m_shutoff_t90_min": 2.2,
        "obs_basal_insulin_frac_peak": 0.09,
        "obs_glucose_threshold_mM": 5.1,
        "obs_hlae_hour0_molecules": 18000.0,
        "obs_cd47_hour0_molecules": 220000.0,
        "obs_nkg2a_lysis_pct": 5.5,
        "obs_serpinb9_hour0_uM": 0.35,
        "obs_regimeII_lysis_2h_pct": 5.0,
        "obs_lot_ocr_ratio": 1.02,
    }
    dt_core14 = wetlab.evaluate_stage0_decision_tree(core_14_panel)
    assert dt_core14["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
    full_ext_input = {
        **core_14_panel,
        "obs_hlae_hour0_molecules_se": 200.0,
        "obs_cd47_hour0_molecules_se": 1000.0,
        "obs_stage7_yield_ratio": 0.985,
        "obs_stage7_yield_se": 0.01,
        "obs_cd47_fold": 8.0,
        "obs_cd47_fold_se": 0.01,
        "obs_serpinb9_log10_kass": 6.4,
        "obs_serpinb9_log10_kass_se": 0.01,
        "obs_icasp9_log10_kill_se": 0.05,
        "obs_pouch_min_po2_mmhg": 23.0,
        "obs_pouch_min_po2_se": 0.01,
        "obs_pouch_max_h2o2_uM": 2.0,
        "obs_pouch_max_h2o2_se": 0.01,
        "obs_pouch_max_ph": 7.40,
        "obs_pouch_max_ph_se": 0.001,
        "obs_ifng_hla_dr_dq_negative_pct": 99.99,
        "obs_ifng_hla_dr_dq_negative_se": 0.001,
        "obs_tp53_ppm1d_bcor_duplex_depth": 62958.0,
        "obs_castseq_duplex_or_amplicon_depth": 62958.0,
        "obs_mcb_effective_clonogenic_trials": 299573228.0,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_nkg2a_lysis_se": 0.5,
        "obs_nkg2c_adaptive_nk_lysis_pct": 5.0,
        "obs_nkg2c_adaptive_nk_lysis_se": 0.5,
        "obs_kir_nkg2a_neg_nk_lysis_pct": 5.5,
        "obs_kir_nkg2a_neg_nk_lysis_se": 0.5,
        "obs_phagocytosis_index_pct": 3.5,
        "obs_phagocytosis_index_se": 0.5,
        "obs_chronic_capsule_thickness_um": 9.5,
        "obs_outer_band_kperf_s_inv": 0.1358,
        "obs_indirect_apc_cd4_ifng_si": 1.4,
        "obs_indirect_apc_dq2_dq8_cd4_stimulation_index": 1.2,
        "obs_gate0d_activated_cd4_fasl_apoptosis_pct": 85.0,
        "obs_gate0d_transwell_beta_annexin_v_pct": 6.0,
        "obs_exogenous_b2m_w632_surface_fold_isotype": 1.05,
        "obs_exogenous_b2m_t1d_ctl_lysis_pct": 2.0,
        "obs_gate0m_glucagon_2mM_fold_increase": 2.4,
        "obs_residual_psc_effective_fraction": 1.0e-10,
        "obs_residual_psc_effective_trials": 299573228.0,
        "obs_20q11_bcl2l1_cnv_ratio": 2.00,
        "obs_20q11_bcl2l1_single_cell_mosaic_fraction": 0.002,
        "obs_hazard_state_prolif_survivor_count_k": 0.0,
        "obs_hazard_state_prolif_effective_trials": 299573228.0,
        "obs_cd47_p10_single_cell_molecules": 210000.0,
        "obs_cd47_p10_single_cell_se": 2000.0,
        "obs_hypoxic_72h_cd47_p10_single_cell_molecules": 182000.0,
        "obs_hypoxic_72h_cd47_p10_single_cell_se": 500.0,
        "obs_hypoxic_72h_hlae_p10_single_cell_molecules": 15200.0,
        "obs_hypoxic_72h_hlae_p10_single_cell_se": 50.0,
        "obs_hypoxic_72h_cd47_retention_molecules": 183354.0,
        "obs_hypoxic_72h_cd47_retention_se": 500.0,
        "obs_hypoxic_72h_hlae_retention_molecules": 15250.0,
        "obs_hypoxic_72h_hlae_retention_se": 50.0,
        "obs_periimplant_fluid_film_um": 35.0,
        "obs_explant_en_bloc_recovery_fraction": 0.995,
        "obs_explant_en_bloc_recovery_se": 0.001,
        "obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h": 0.08,
        "obs_gate0s_ap1903_1h_insulin_burst_iu": 1.8,
        "obs_ap1903_pulsed_pk_log10_kill": 4.6,
        "obs_ap1903_monotherapy_pulsed_pk_log10_kill": 4.5,
        "obs_ap1903_step1_sublethal_survivor_count_k": 0.0,
        "obs_ap1903_step1_effective_trials": 29957323.0,
        "obs_hypoxic_15mmhg_f_isr_translation_ratio": 0.78,
        "obs_iec_cd47_p10_molecules": 195000.0,
        "obs_iec_cd47_p10_se": 2000.0,
        "obs_iec_abo_ab_antigen_mfi_ratio": 1.01,
        "obs_iec_abo_ab_antigen_mfi_se": 0.005,
        "obs_post_stress_nkg2d_dnam1_attributed_lysis_pct": 2.5,
        "obs_senescent_p16_sasp_positive_cell_pct": 0.3,
        "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct": 5.2,
        "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_se": 0.5,
        "obs_il12_il18_sirpa_low_nk_lysis_pct": 6.4,
        "obs_il12_il18_sirpa_low_nk_lysis_se": 0.6,
        "obs_cd47_cc2c6_over_b6h12_ratio": 0.94,
        "obs_complement_c3b_opsonized_phagocytosis_pct": 3.2,
        "obs_complement_c3b_opsonized_phagocytosis_se": 0.4,
        "obs_iec_vascular_support_viability_pct": 91.5,
        "obs_recipient_insulin_clearance_tau_min": 6.5,
        "obs_12pulse_unprimed_hypoxic_serial_lysis_pct": 5.8,
        "obs_hlae_nkg2c_over_nkg2a_kd_ratio": 4.8,
        "obs_immaturity_corner_min_glucose_mM": 3.85,
        "obs_alpha_counterregulatory_gain_ratio": 0.92,
        "obs_ucn3_maturity_ratio": 0.88,
        "obs_combined_stress_embodiment_a_star_lysis_pct": 5.8,
        "obs_combined_stress_comp_r4_lysis_pct": 5.2,
        "obs_combined_stress_comp_r4_gsis_si": 2.9,
        "obs_combined_stress_comp_r4_viable_beta_mass_pct": 91.0,
        "obs_combined_stress_comp_r4_shutoff_t90_min": 2.1,
        "obs_combined_stress_comp_r4_max_h2o2_uM": 2.8,
        "obs_combined_stress_comp_r4_max_ph": 7.38,
        "obs_gate0v_cvb_day14_log10_viral_reduction": 4.2,
        "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct": 0.0,
        "obs_gate0v_day28_persistent_vp1_ie1_cells_scored": 10000.0,
        "obs_post_stress_20q11_bcl2l1_outgrowth_fold": 0.92,
        "obs_assembled_cassette_chandler_loop_sc5b9_fold": 1.25,
        "obs_assembled_cassette_chandler_loop_platelet_retention_pct": 86.0,
        "obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction": 0.95,
        "obs_apoptotic_effero_clearance_retention_pct_of_wt": 92.0,
        "obs_nutlin3a_p21_mdm2_fold_induction": 8.4,
        "obs_vascular_pericyte_ng2_pdgfrb_coverage_pct": 78.0,
        "obs_assembled_cluster_slc18a1_enterochromaffin_pct": 0.3,
        "obs_assembled_cluster_ins_gcg_polyhormonal_pct": 0.8,
        "obs_assembled_cluster_ki67_progenitor_pct": 0.1,
        "obs_post_stress_mafa_nkx61_beta_identity_retention_pct": 91.0,
    }
    dt_full_ext = wetlab.evaluate_stage0_decision_tree(full_ext_input)
    assert dt_full_ext["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

    # - R-04: C1/C1 B'-SCT without explicit hlac_4digit_relation or MLC clearance holds on HOLD_B_PRIME_SCT_C1_UNSPECIFIED_4DIGIT_RELATION_REQUIRES_ALLELE_TYPING_AND_CD8_MLC_CLEARANCE
    pt_c1_untyped = {
        "patient_id": "PT-C1-UNTYPED",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C1",
    }
    bp_c1_untyped = l5_mod.compile_personalized_patient_blueprint(pt_c1_untyped, track2_subbank_mode="B_PRIME_SCT")
    assert bp_c1_untyped["product_qualification_status"] == "HOLD_B_PRIME_SCT_C1_UNSPECIFIED_4DIGIT_RELATION_REQUIRES_ALLELE_TYPING_AND_CD8_MLC_CLEARANCE"
    assert bp_c1_untyped["routing_holds_cleared"] is False
    assert bp_c1_untyped["lot_release_cleared"] is False

    # - R-05 & R-06: Literature-provenanced HLA-E KD table (Lin 2023 PMID:37264229 1.70/4.97 uM = 2.9235x for VMAPRTLLL) & 137/2424 main-effect diagonal
    nkg2c_audit = l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]["nkg2a_vs_nkg2c_2d_synaptic_selectivity_bound_sonnet55_r4"]
    assert nkg2c_audit["soluble_3d_spr_kd_reference"]["vmaprtlll_hla_a1_c_leader_challenger_kd_uM"]["nkg2a_cd94_kd_uM"] == 1.70
    assert nkg2c_audit["soluble_3d_spr_kd_reference"]["vmaprtlll_hla_a1_c_leader_challenger_kd_uM"]["nkg2c_cd94_kd_uM"] == 4.97
    assert nkg2c_audit["soluble_3d_spr_kd_reference"]["covalent_sct_linker_selectivity_qualification_status"] == "UNPROVEN_IN_SCT_FORMAT_REQUIRES_GATE_0A_SPR"
    assert co_prim["mode1_23col_mean_main_effect_var_beta_over_sigma2_across_6_factors"].startswith("137/2424")

    # - R-07: Sticky S9 reset requires 3-field QP governance and retains prior_s9_safety_stop_history_retained
    dt_sticky_bad_reset = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "prior_variant_safety_status": "HARD_SAFETY_STOP_REJECTED",
            "reset_prior_safety_stop_after_investigation": True,
        }
    )
    assert dt_sticky_bad_reset["safety_kill_s9_triggered"] is True
    assert dt_sticky_bad_reset["safety_kill_s9_status"] == "HARD_SAFETY_STOP_STICKY_PRIOR_S9_RETAINED"
    dt_sticky_qp_cleared = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "prior_variant_safety_status": "HARD_SAFETY_STOP_REJECTED",
            "reset_prior_safety_stop_after_investigation": True,
            "reset_prior_safety_stop_reason_code": "ASSAY_ARTIFACT_CONFIRMED_BY_REPEAT_DUPLEX_SEQ",
            "reset_prior_safety_stop_repeat_batch_ids": ["BATCH-R1", "BATCH-R2"],
            "verified_ledger_passing_s9_batch_ids": ["BATCH-R1", "BATCH-R2"],
            "reset_prior_safety_stop_independent_qp_approver_id": "QP-001",
        }
    )
    assert dt_sticky_qp_cleared["safety_kill_s9_triggered"] is False
    assert dt_sticky_qp_cleared["resolved_architecture_assignments"]["prior_s9_safety_stop_history_retained"] is True
    assert dt_sticky_qp_cleared["resolved_architecture_assignments"]["lot_maturation_status"] == "READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY"

    # 9. AUDIT-REF (DEC-043):
    # - AUDIT-REF: Composite 4-way lead uses HLA-C*03:04 (not GRAVEYARD_005 HLA-C*07:02) and AsCas12a Ultra, and combined-stress coupon gate is wired into evaluate_stage0_decision_tree
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        l1 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        comp_r4_spec = l1["champion_challenger_portfolio"]["abxyz_champion_contender_wildcard_tournament_registry_dec020"]["tier1_stage0_execution_triad_80pct_focus"]["X_Locked_Number1_Evolutionary_Contender"]["composite_4way_research_lead_specification_muse_r4"]
        assert "HLA-C*03:04" in comp_r4_spec["locus_2_ciita"]
        assert "AsCas12a Ultra" in comp_r4_spec["locus_2_ciita"]
    dt_cs_promote = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_combined_stress_embodiment_a_star_lysis_pct": 11.5,
            "obs_combined_stress_comp_r4_lysis_pct": 6.2,
            "obs_combined_stress_comp_r4_gsis_si": 2.9,
            "obs_combined_stress_comp_r4_viable_beta_mass_pct": 91.0,
            "obs_combined_stress_comp_r4_shutoff_t90_min": 2.1,
            "obs_combined_stress_comp_r4_max_h2o2_uM": 2.8,
            "obs_combined_stress_comp_r4_max_ph": 7.38,
        }
    )
    assert dt_cs_promote["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "PROMOTE_ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500"
    assert dt_cs_promote["resolved_architecture_assignments"]["combined_stress_coupon_promoted_variant"] == "ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500"
    # - D5-01 & D5-02: Duplicate repeat_batch_ids ["B1", "B1"] fails to clear sticky S9, and FRACTION_HALF_EFFECT uses -0.05 threshold
    dt_dup_batches = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "prior_variant_safety_status": "HARD_SAFETY_STOP_REJECTED",
            "reset_prior_safety_stop_reason_code": "INSTRUMENT_CONTAMINATION_VERIFIED_NEGATIVE_RETEST",
            "reset_prior_safety_stop_repeat_batch_ids": ["BATCH-DUP", "BATCH-DUP"],
            "reset_prior_safety_stop_independent_qp_approver_id": "QP-002",
        }
    )
    assert dt_dup_batches["safety_kill_s9_triggered"] is True
    for scale_name, beta_val in (
        ("FRACTION_HALF_EFFECT", -0.06),
        ("PP_HALF_EFFECT", -6.0),
        ("PP_FULL_CONTRAST", -12.0),
        ("DID_PP", -24.0),
    ):
        dt_scale = wetlab.evaluate_stage0_decision_tree(
            {
                "obs_icasp9_log10_kill": 4.5,
                "obs_tp53_ppm1d_bcor_vaf": 0.0,
                "obs_castseq_translocation_freq": 0.0,
                "s8b_epistasis_estimand_scale": scale_name,
                "obs_beta_x1_x4_epistasis": beta_val,
                "obs_beta_x1_x4_epistasis_se": abs(beta_val) / 6.0,
            }
        )
        assert dt_scale["switch_decisions"]["S8b_TYK2_Inducible_Promoter_Epistasis_Gate_0E"] == "COUPLE_TYK2_P1104A_STRICTLY_WITH_ARM1X_AND_L3_BFB_CONSTITUTIVE_BASAL_FLOOR"

    # 10. AUDIT-REF [DEC-044]:
    # - R5-01: Combined-stress coupon fails closed when comparator lysis/SI are missing, when COMP-R4 secondary viability/shutoff/chemistry fail, or when a component fallback switch triggers
    dt_cs_missing_comparators = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_combined_stress_comp_r4_lysis_pct": 6.0,
            "obs_combined_stress_comp_r4_gsis_si": 3.0,
        }
    )
    assert dt_cs_missing_comparators["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMBINED_STRESS_INCOMPLETE_COMPARATOR_EVIDENCE"
    assert dt_cs_missing_comparators["resolved_architecture_assignments"].get("combined_stress_coupon_promoted_variant") is None

    dt_cs_comp_fail = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_hlae_hour0_molecules": 8000.0,  # triggers S2 hold without Arm 1X qualification
            "obs_combined_stress_embodiment_a_star_lysis_pct": 12.0,
            "obs_combined_stress_comp_r4_lysis_pct": 6.0,
            "obs_combined_stress_comp_r4_gsis_si": 3.0,
        }
    )
    assert dt_cs_comp_fail["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMPONENT_GATE_FAILURE_BLOCKS_COMPOSITE_PROMOTION"
    assert dt_cs_comp_fail["resolved_architecture_assignments"].get("combined_stress_coupon_promoted_variant") is None

    dt_cs_sec_fail = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_combined_stress_embodiment_a_star_lysis_pct": 12.0,
            "obs_combined_stress_comp_r4_lysis_pct": 6.0,
            "obs_combined_stress_comp_r4_gsis_si": 3.0,
            "obs_combined_stress_comp_r4_viable_beta_mass_pct": 80.0,  # < 85%
        }
    )
    assert dt_cs_sec_fail["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMBINED_STRESS_COUPON_FAILED"

    # - R5-02: Sticky S9 reset requires reset_prior_safety_stop_after_investigation=True, non-self repeat_batch_ids, and verified ledger passing batches
    dt_sticky_flag_false = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "prior_variant_safety_status": "HARD_SAFETY_STOP_REJECTED",
            "reset_prior_safety_stop_after_investigation": False,
            "reset_prior_safety_stop_reason_code": "ASSAY_ARTIFACT_CONFIRMED_BY_REPEAT_DUPLEX_SEQ",
            "reset_prior_safety_stop_repeat_batch_ids": ["BATCH-R1", "BATCH-R2"],
            "verified_ledger_passing_s9_batch_ids": ["BATCH-R1", "BATCH-R2"],
            "reset_prior_safety_stop_independent_qp_approver_id": "QP-001",
        }
    )
    assert dt_sticky_flag_false["safety_kill_s9_triggered"] is True

    dt_sticky_self_batch = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "batch_id": "BATCH-R1",
            "prior_variant_safety_status": "HARD_SAFETY_STOP_REJECTED",
            "reset_prior_safety_stop_after_investigation": True,
            "reset_prior_safety_stop_reason_code": "ASSAY_ARTIFACT_CONFIRMED_BY_REPEAT_DUPLEX_SEQ",
            "reset_prior_safety_stop_repeat_batch_ids": ["BATCH-R1", "BATCH-R2"],
            "verified_ledger_passing_s9_batch_ids": ["BATCH-R1", "BATCH-R2"],
            "reset_prior_safety_stop_independent_qp_approver_id": "QP-001",
        }
    )
    assert dt_sticky_self_batch["safety_kill_s9_triggered"] is True

    dt_sticky_unverified_ledger = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "prior_variant_safety_status": "HARD_SAFETY_STOP_REJECTED",
            "reset_prior_safety_stop_after_investigation": True,
            "reset_prior_safety_stop_reason_code": "ASSAY_ARTIFACT_CONFIRMED_BY_REPEAT_DUPLEX_SEQ",
            "reset_prior_safety_stop_repeat_batch_ids": ["FABRICATED-B1", "FABRICATED-B2"],
            "verified_ledger_passing_s9_batch_ids": ["BATCH-R1"],
            "reset_prior_safety_stop_independent_qp_approver_id": "QP-001",
        }
    )
    assert dt_sticky_unverified_ledger["safety_kill_s9_triggered"] is True

    # - R5-03: Per-endpoint epistasis scale tracking and unified -0.05 fractional threshold
    dt_epistasis_mixed_scale = wetlab.evaluate_stage0_decision_tree(
        {
            "obs_icasp9_log10_kill": 4.5,
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_castseq_translocation_freq": 0.0,
            "obs_beta_x1_x4_epistasis": -0.08,
            "obs_beta_x1_x4_epistasis_se": 0.02,
            "obs_beta_x1_x4_epistasis_scale": "FRACTION_HALF_EFFECT",
            "obs_beta_x3_x4_epistasis": -1.0,
            "obs_beta_x3_x4_epistasis_se": 1.0,
            "obs_beta_x3_x4_epistasis_scale": "PP_HALF_EFFECT",
            "s8b_epistasis_estimand_scale": "PP_HALF_EFFECT",
        }
    )
    assert dt_epistasis_mixed_scale["switch_decisions"]["S8b_TYK2_Inducible_Promoter_Epistasis_Gate_0E"] == "COUPLE_TYK2_P1104A_STRICTLY_WITH_ARM1X_AND_L3_BFB_CONSTITUTIVE_BASAL_FLOOR"

    # - R5-04: Layer-5 compiler separates in_silico_candidate_routing_cleared from lot_release_cleared
    l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    for p_entry in l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["representative_patient_compilations"]:
        assert p_entry["in_silico_candidate_routing_cleared"] == (
            p_entry["selected_product_qualification_status"] == "SEQUENCE_COMPLETE_IN_SILICO_CANDIDATE_REQUIRING_WETLAB_RELEASE_QC"
        )
        assert p_entry["lot_release_cleared"] is False

    # - R5-05 & R5-S2-01: CAST-seq depth < 62,958 hold, MCB clonogenic trials < 299,573,228 hold, and DX-S requires >= 299,573,227 effective clonogenic trials
    dt_cast_low_depth = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_castseq_mutant_junctions_k": 0,
            "obs_castseq_duplex_or_amplicon_depth": 35000,
        }
    )
    assert dt_cast_low_depth["safety_kill_s9_triggered"] is False
    assert dt_cast_low_depth["safety_kill_s9_status"] == "HOLD_CASTSEQ_DUPLEX_DEPTH_BELOW_62958_PROTOCOL_MINIMUM"

    dt_mcb_low_trials = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_mcb_effective_clonogenic_trials": 1.0e7,
            "obs_mcb_clonogenic_survivor_count_k": 0,
        }
    )
    assert dt_mcb_low_trials["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_MCB_CLONOGENIC_TRIALS_BELOW_299573228_MINIMUM"

    dt_dxs_missing_n = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_dxs_locus4_joint_retention_fraction": 0.99995,
            "obs_dxs_cd47low_icasp9low_fraction": 5.0e-7,
            "obs_dxs_ap1903_cd47high_escape_fraction": 5.0e-9,
        }
    )
    assert dt_dxs_missing_n["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] == "HOLD_DXS_RARE_EVENT_REQUIRES_EXPLICIT_SE_OR_CELL_COUNT"

    # - R5-06..R5-14: Verify Layer 2, Layer 3, and Layer 4 Round-5 qualification blocks
    l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))

    # R5-06: Chronic capsule thickness sweep 15..100 um + k_perf = 0
    cap_sweep = (
        l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
        ["design_x_100um_lamina_transport_qualification_v80_018"]
        ["joint_adverse_corner_stress_audit_muse_r1"]
        ["chronic_zero_source_and_mineralization_extension_muse_r3"]
        ["physical_outer_annulus_masked_chronic_zero_source_audit_argon_r4"]
        ["joint_chronic_source_free_areal_loading_capsule_ocr_envelope_sol61_r4"]
        ["chronic_fibrotic_capsule_thickness_sweep_50_to_100um_sol61_r5"]
    )
    assert cap_sweep["computation_mode"] == "LIVE_NZ201_NEWTON_BVP_AND_SPHERICAL_CLUSTER_DERIVED_SOL61_R5"
    assert cap_sweep["core_viability_5mmhg_max_capsule_um_100um_cluster"] == 30.0
    assert cap_sweep["core_viability_5mmhg_max_capsule_um_80um_cluster"] == 100.0
    dt_cap_40 = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_capsule_dz_um": 40.0})
    assert dt_cap_40["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S9_CHRONIC_CAPSULE_GT_30UM_REQUIRES_80UM_LOWLOAD_OR_OMENTAL_QUALIFICATION"
    dt_cap_90 = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_capsule_dz_um": 90.0, "s9_chronic_80um_lowload_qualified": True})
    assert dt_cap_90["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S9_CHRONIC_CAPSULE_GT_50UM_OR_AVASCULAR_DROPOUT_REQUIRES_OMENTAL_REVASCULARIZATION"

    # R5-07: 12-pulse hypoxic coupon reconciliation and synaptic suppression gate
    arm1b_construct = next(c for c in l2["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    r5_07_block = arm1b_construct["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]["Locus3_BasalFloor_ScarFree_SERPINB9_PDL1_challenger_muse_r3"]["hour0_kinetic_race_and_stoichiometric_bolus_model_muse_r3"]["serpinb9_suicide_substrate_serial_pulse_ode_and_perk_eif2a_isr_audit_sol61_r4"]["twelve_pulse_hypoxic_coupon_reconciliation_and_synaptic_suppression_gate_sol61_r5"]
    assert r5_07_block["finding_id"] == "AUDIT-REF"
    assert r5_07_block["required_upstream_synaptic_degranulation_probability_cap_p_degran_max"] == 0.15

    # R5-08 & R5-09 & R5-13: Layer 3 pure NKG2A- NKG2C+ audit, 4-construct SCT comparator dict, and QPCTL / complement phagocytosis gate
    hlac_r4 = l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]
    assert "pure_nkg2a_neg_nkg2c_pos_subset_and_hlae_d227k_t228a_ctl_audit_sol61_r5" in hlac_r4["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]["nkg2a_vs_nkg2c_2d_synaptic_selectivity_bound_sonnet55_r4"]
    sct_4way = hlac_r4["f_pocket_lateral_channel_mutation_required_for_hlac_sct"]["disulfide_trap_comparators_disentangled_kimi31_r4"]
    assert len(sct_4way) == 4
    assert "L3_B_SCT_Y84A_L1_20AA" in sct_4way
    assert "3.5169 A" in sct_4way["L3_B_SCT_Y84A_L1_20AA"]
    assert "qpctl_pyroglutamate_and_opsonized_complement_phagocytosis_gate_sol61_r5" in hlac_r4["endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4"]

    # R5-10, R5-11, R5-12, R5-14: Indirect APC, common-mode iCasp9, exercise counter-regulatory margin, and B2M free heavy-chain rescue
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        assert "indirect_apc_hip_cd4_gate_0d_sol61_r5" in comp_r4_spec
    assert "b2m_ko_free_heavy_chain_surface_leakage_and_exogenous_b2m_rescue_gate_sol61_r5" in l2["stage0_to_stage7_perturbation_plane"]["er_heavy_chain_flux_and_upr_hazard_muse_r3"]
    assert "challenged_v0p60_300k_beta_eq_exercise_margin_audit_sol61_r5" in l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]

    # 11. AUDIT-REF [DEC-045]:
    # - AUDIT-REF: MCB clonogenic survivor count > 0 triggers HARD_SAFETY_STOP_MCB_CLONOGENIC_SURVIVOR_DETECTED, and missing CAST-seq/MCB trials keeps readiness tier provisional
    dt_mcb_survivor = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_mcb_clonogenic_survivor_count_k": 1})
    assert dt_mcb_survivor["safety_kill_s9_triggered"] is True
    assert dt_mcb_survivor["resolved_architecture_assignments"]["mcb_clonogenic_escape_qualification_status"] == "HARD_SAFETY_STOP_MCB_CLONOGENIC_SURVIVOR_DETECTED"
    dt_missing_mcb = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_tp53_ppm1d_bcor_duplex_depth": 62958.0,
            "obs_castseq_duplex_or_amplicon_depth": 62958.0,
            "obs_nkg2c_adaptive_nk_lysis_pct": 5.0,
            "obs_kir_nkg2a_neg_nk_lysis_pct": 6.0,
            "obs_kir_nkg2a_neg_nk_lysis_se": 0.5,
            "obs_phagocytosis_index_pct": 3.5,
            "obs_chronic_capsule_thickness_um": 9.5,
            "obs_capsule_dz_se": 0.2,
        }
    )
    assert dt_missing_mcb["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"

    # - AUDIT-REF: Numeric obs_b_prime_sct_cd8_lysis_pct > 10% triggers S7 directly and overrides b_prime_sct_recipient_cd8_mlc_cleared=True
    dt_sct_cd8 = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_b_prime_sct_cd8_lysis_pct": 14.0,
            "b_prime_sct_recipient_cd8_mlc_cleared": True,
        }
    )
    assert dt_sct_cd8["switch_decisions"]["S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B"] == "HOLD_B_PRIME_SCT_DIRECT_CD8_ALLOREACTIVITY_FAILED_RETAIN_PURE_MHC_NULL_TRACK1"
    assert dt_sct_cd8["resolved_architecture_assignments"]["track2_kir_subbank_construct"] == "HOLD_B_PRIME_SCT_DIRECT_CD8_ALLOREACTIVITY_FAILED"
    assert dt_sct_cd8["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_B_PRIME_ALLO_FAILED_REQUIRES_4DIGIT_MATCH_OR_SCT_CD8_MLC_CLEARANCE"

    # - AUDIT-REF: require_full_7metric_combined_stress_panel fails closed when secondary panel metrics are omitted, and S5 bystander failure blocks composite promotion
    dt_cs_7metric_hold = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "require_full_7metric_combined_stress_panel": True,
            "obs_combined_stress_embodiment_a_star_lysis_pct": 11.5,
            "obs_combined_stress_comp_r4_lysis_pct": 6.2,
            "obs_combined_stress_comp_r4_gsis_si": 2.9,
        }
    )
    assert dt_cs_7metric_hold["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL"
    dt_s5_hold = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "s5_fasl_delta_tunel_gt_5pp": True})
    assert dt_s5_hold["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S5_FASL_BYSTANDER_TOXICITY_REQUIRES_THBD_FHL1_ONLY_CUFF_QUALIFICATION"

    # - AUDIT-REF: Forced Track 1 on high-KIR recipient holds without track1_sorted_kir_nkg2a_neg_nk_lysis_cleared=True, and Bw4/C1C2 holds are enforced
    pt_high_kir_forced_t1 = {
        "patient_id": "PT-HIGH-KIR-FORCED-T1",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C1",
        "requested_track": "TRACK_1_OFF_THE_SHELF",
    }
    bp_forced_t1 = l5_mod.compile_personalized_patient_blueprint(pt_high_kir_forced_t1)
    assert bp_forced_t1["product_qualification_status"] == "HOLD_FORCED_TRACK1_OVERRIDES_HIGH_KIR_NK_REQUIRES_SORTED_KIR_NK_CLEARANCE"
    assert bp_forced_t1["routing_holds_cleared"] is False

    # - AUDIT-REF: Exogenous B2M rescue, indirect APC Gate 0D, and glucagon 2mM fold increase holds
    dt_exo_b2m = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_exogenous_b2m_w632_surface_fold_isotype": 2.4})
    assert dt_exo_b2m["resolved_architecture_assignments"]["er_rescue_genotype_pivot"] == "PROMOTE_ISLET1_A4_K6_HLAA_HLAB_KO_EXOGENOUS_B2M_SHIELD"
    assert dt_exo_b2m["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S6_ER_UPR_THRESHOLD_BREACHED_REQUIRES_GENOTYPE_REQUALIFICATION"
    dt_gate0d = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_indirect_apc_cd4_ifng_si": 2.8})
    assert dt_gate0d["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_GATE0D_INDIRECT_APC_DQ2_DQ8_CD4_AUTOIMMUNE_STIMULATION_EXCEEDED"
    dt_gcg = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_gate0m_glucagon_2mM_fold_increase": 1.4})
    assert dt_gcg["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_FOR_PLUS_3W_EXTENDED_STAGE7_MATURATION"

    # 12. AUDIT-REF [DEC-046]:
    # - R5-S1-01: Supplying obs_mcb_effective_clonogenic_trials >= 299,573,228 without obs_mcb_clonogenic_survivor_count_k fails closed
    dt_mcb_missing_k = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_mcb_effective_clonogenic_trials": 300000000.0,
        }
    )
    assert dt_mcb_missing_k["resolved_architecture_assignments"]["mcb_clonogenic_escape_qualification_status"] == "HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING"
    assert dt_mcb_missing_k["safety_kill_s9_status"] == "HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING"
    assert dt_mcb_missing_k["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_MCB_CLONOGENIC_SURVIVOR_OUTCOME_MISSING"

    # - R5-S1-02: 3-metric combined-stress input without secondary metrics holds by default
    dt_cs_3metric_default_hold = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_combined_stress_embodiment_a_star_lysis_pct": 11.5,
            "obs_combined_stress_comp_r4_lysis_pct": 6.2,
            "obs_combined_stress_comp_r4_gsis_si": 2.9,
        }
    )
    assert dt_cs_3metric_default_hold["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL"
    assert dt_cs_3metric_default_hold["resolved_architecture_assignments"].get("combined_stress_coupon_promoted_variant") is None

    # - R5-S1-07: C1/C2 B'-SCT without explicit c1c2_heteroallelic_phased_longread_ddpcr_verified fails closed
    pt_c1c2_unphased = {
        "patient_id": "PT-C1C2-UNPHASED",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C2",
        "hlac_4digit_relation": "SELF_4DIGIT_MATCH",
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared": True,
        "b_prime_sct_kir2ds1_non_activation_cleared": True,
    }
    bp_c1c2_unphased = l5_mod.compile_personalized_patient_blueprint(pt_c1c2_unphased, track2_subbank_mode="B_PRIME_SCT")
    assert "HOLD_C1C2_HETEROALLELIC_CIITA_SCT_REQUIRES_PHASED_LONGREAD_AND_DDPCR_VERIFICATION" in bp_c1c2_unphased["all_active_qualification_holds"]
    assert bp_c1c2_unphased["routing_holds_cleared"] is False
    assert bp_c1c2_unphased["b_prime_4digit_hlac_and_b_prime_sct_routing_audit_muse_r3"]["c1c2_heteroallelic_phased_longread_ddpcr_verified"] is None
    assert bp_c1c2_unphased["b_prime_4digit_hlac_and_b_prime_sct_routing_audit_muse_r3"]["c1c2_heteroallelic_phased_verification_status"] == "HOLD_UNMEASURED_MISSING_PHASING_EVIDENCE"

    # 13. AUDIT-REF [DEC-047]:
    # - P-01 / P-13 / L16: Residual PSC > 1e-8, 20q11.21 BCL2L1 CNV outside [1.85, 2.15], and hazard-state proliferative survivor > 0 trigger hard S9 safety stops
    dt_psc_fail = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_residual_psc_effective_fraction": 5e-7})
    assert dt_psc_fail["safety_kill_s9_triggered"] is True
    assert dt_psc_fail["safety_kill_s9_status"] == "HARD_SAFETY_STOP_REJECTED"
    assert dt_psc_fail["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"
    dt_cnv_fail = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_20q11_bcl2l1_cnv_ratio": 2.85})
    assert dt_cnv_fail["safety_kill_s9_triggered"] is True
    assert dt_cnv_fail["safety_kill_s9_status"] == "HARD_SAFETY_STOP_REJECTED"
    assert dt_cnv_fail["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"
    dt_haz_surv = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_hazard_state_prolif_survivor_count_k": 1.0})
    assert dt_haz_surv["safety_kill_s9_triggered"] is True
    assert dt_haz_surv["safety_kill_s9_status"] == "HARD_SAFETY_STOP_REJECTED"
    assert dt_haz_surv["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"

    # - P-02 / P-04 / P-09 / P-12 / L07 / L18 / L20:
    #   (a) One-sided 95% upper confidence bound on NKG2A+ lysis (7.5% + 1.644854 * 3.2% = 12.76% > 12.0%) holds lot
    dt_upper95_nk = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_nkg2a_lysis_pct": 7.5, "obs_nkg2a_lysis_se": 3.2}
    )
    assert dt_upper95_nk["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_EXCESSIVE_DONOR_VARIANCE_OR_95PCT_UPPER_CONFIDENCE_BOUND_EXCEEDED"
    #   (b) NKG2A- KIR+ NK lysis > 8% does NOT mis-route to S2 Arm 1X HLA-E; instead triggers Track 2 B'-SCT / HLA-G1 challenger hold
    dt_kir_neg_nk = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_kir_nkg2a_neg_nk_lysis_pct": 11.0})
    assert dt_kir_neg_nk["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "RETAIN_EMBODIMENT_A_STAR_ARM1B"
    assert dt_kir_neg_nk["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_KIR_NKG2A_NEG_NK_OR_IL12_IL18_SIRPA_LOW_LYSIS_EXCEEDED_REQUIRES_TRACK2_B_SCT_OR_HLAG1_CHALLENGER"
    #   (c) Simultaneous HLA-E < 15k and CD47 < 50k unmasks S3 (DX-S split cassette) when DX-S rare-event panel is qualified
    dt_dxs_unmasked = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_hlae_hour0_molecules": 11000.0,
            "obs_cd47_hour0_molecules": 35000.0,
            "obs_dxs_locus4_joint_retention_fraction": 0.99999,
            "obs_dxs_locus4_joint_retention_se": 1e-6,
            "obs_dxs_cd47low_icasp9low_fraction": 0.0,
            "obs_dxs_cd47low_icasp9low_se": 1e-7,
            "obs_dxs_ap1903_cd47high_escape_fraction": 0.0,
            "obs_dxs_ap1903_cd47high_escape_se": 1e-9,
            "obs_dxs_screened_cell_count": 900000000.0,
            "obs_dxs_clonogenic_fraction": 1.0,
            "obs_dxs_spike_in_detection_efficiency": 0.35,
        }
    )
    assert dt_dxs_unmasked["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] == "PROMOTE_ISLET1_A4_DXS_ARM1C_L4CD47"
    assert dt_dxs_unmasked["resolved_architecture_assignments"]["locus1_construct"] == "ISLET1_A4_DXS_ARM1C_L4CD47"
    #   (d) Playbook CD47 cliff (55,000 < 200,000) or phagocytosis cliff (7.5% > 5.0%) holds lot
    dt_cd47_cliff = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_cd47_hour0_molecules": 55000.0})
    assert dt_cd47_cliff["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED"

    # - P-03 / P-07 / P-08 / L12: Subthreshold k_perf (0.04 < 0.10 s^-1), peri-implant fluid film (> 100 um), explant recovery (< 0.99), and AP1903 24h insulin dump (> 0.15) hold lot
    dt_kperf_sub = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_outer_band_kperf_s_inv": 0.04})
    assert dt_kperf_sub["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S9_SUBTHRESHOLD_OUTER_BAND_PERFUSION_BELOW_0P10_S_INV_REQUIRES_OMENTAL_REVASCULARIZATION"
    dt_film_fail = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_periimplant_fluid_film_um": 135.0})
    assert dt_film_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S9_PERIIMPLANT_FLUID_FILM_GT_100UM_SEROMA_HYPOXIA_BREACH"
    dt_dump_fail = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h": 0.28})
    assert dt_dump_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_GATE0S_AP1903_INSULIN_DUMP_EXCEEDS_15PCT_REQUIRES_STAGED_CASSETTE_EXCISION_AND_CLAMP"

    # - Layer 1..5 AUDIT-REF audit blocks present and Layer-5 modular cassette & hazard flags verified
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        assert "unconditional_patient_level_tumor_harm_and_two_stage_psc_purge_model_sonnet55_r5" in l1
    assert "coupled_five_transgene_hypoxic_isr_surface_decay_and_leader_rescue_audit_sonnet55_r5" in l2["stage0_to_stage7_perturbation_plane"]
    assert "hlag1_sct_lilrb1_challenger_and_il12_il18_sirpa_low_nk_audit_sonnet55_r5" in l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4"]
    assert "bed_level_convective_o2_budget_and_periimplant_fluid_film_audit_sonnet55_r5" in l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]
    assert "icasp9_ap1903_ablation_insulin_dump_and_staged_cassette_excision_ode_audit_sonnet55_r5" in l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]
    assert "population_simulation_with_cmv_nkg2c_and_sensitization_sonnet55_r5" in l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["b_prime_4digit_hlac_alloreactivity_and_b_prime_sct_pivot_audit_muse_r3"]["hazard_active_cohort_sensitivity_audit_sonnet55_r4"]

    # 14. AUDIT-REF [DEC-048]:
    # - R1 (MW-R5-S2-F01): Direct biological gate UCBs (8% lysis, 10% IL12/18, 5% phago) & catastrophic NKG2A > 15% precedence over Arm 1X / UCB
    dt_ucb_nkg2a_6p5_se1p5 = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_nkg2a_lysis_pct": 6.5, "obs_nkg2a_lysis_se": 1.5}
    )  # 6.5 + 1.644854 * 1.5 = 8.967% > 8.0%
    assert dt_ucb_nkg2a_6p5_se1p5["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_EXCESSIVE_DONOR_VARIANCE_OR_95PCT_UPPER_CONFIDENCE_BOUND_EXCEEDED"
    dt_cat_nkg2a_18_hlae_10k = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_nkg2a_lysis_pct": 18.0, "obs_hlae_hour0_molecules": 10000.0}
    )
    assert dt_cat_nkg2a_18_hlae_10k["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "HOLD_S2_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_REQUIRES_MULTI_RECEPTOR_RESCUE"
    assert dt_cat_nkg2a_18_hlae_10k["resolved_architecture_assignments"]["locus1_construct"] == "HOLD_CATASTROPHIC_NKG2A_LYSIS_EXCEEDS_ARM1X_RESCUE_ENVELOPE"
    assert dt_cat_nkg2a_18_hlae_10k["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_EXCEEDS_ARM1X_RESCUE_ENVELOPE"

    # - R2 (MW-R5-S2-F02): L3-BF-B or DX-S alone cannot clear sorted KIR+ NKG2A- NK failure; requires track2_b_sct_or_hlag1_lilrb1_rescue_qualified=True
    dt_kir_l3bfb_only = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_kir_nkg2a_neg_nk_lysis_pct": 12.0, "s1_l3_bfb_replacement_qualified": True}
    )
    assert dt_kir_l3bfb_only["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_KIR_NKG2A_NEG_NK_OR_IL12_IL18_SIRPA_LOW_LYSIS_EXCEEDED_REQUIRES_TRACK2_B_SCT_OR_HLAG1_CHALLENGER"
    dt_kir_track2_qualified = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_kir_nkg2a_neg_nk_lysis_pct": 12.0, "track2_b_sct_or_hlag1_lilrb1_rescue_qualified": True}
    )
    assert dt_kir_track2_qualified["resolved_architecture_assignments"]["locus1_construct"] == "PROMOTE_TRACK2_B_SCT_OR_HLAG1_LILRB1_RESCUE"
    assert dt_kir_track2_qualified["resolved_architecture_assignments"]["lot_maturation_status"] == "CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION"
    dt_kir_track2_switched = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_kir_nkg2a_neg_nk_lysis_pct": 12.0,
            "track2_b_sct_or_hlag1_lilrb1_rescue_qualified": True,
            "track2_b_sct_or_hlag1_product_switched": True,
        }
    )
    assert dt_kir_track2_switched["resolved_architecture_assignments"]["lot_maturation_status"] == "READY_FOR_ENCAPSULATION"

    # - R3 (MW-R5-S2-F03): DX-S requires explicit obs_dxs_clonogenic_fraction and obs_dxs_spike_in_detection_efficiency >= 0.30
    dt_dxs_missing_eta = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_cd47_hour0_molecules": 35000.0,
            "obs_dxs_locus4_joint_retention_fraction": 0.99999,
            "obs_dxs_cd47low_icasp9low_fraction": 0.0,
            "obs_dxs_ap1903_cd47high_escape_fraction": 0.0,
            "obs_dxs_screened_cell_count": 300000000.0,
        }
    )
    assert dt_dxs_missing_eta["switch_decisions"]["S3_Locus1_Constitutive_Split_Fallback_DX_S"] == "HOLD_DXS_RARE_EVENT_REQUIRES_EXPLICIT_SE_OR_CELL_COUNT"

    # - R4 (MW-R5-S2-F04): Combined-stress RETAIN branch enforces require_full_7metric_combined_stress_panel, and omitting DEC-047 endpoints keeps readiness tier provisional
    dt_cs_retain_missing_sec = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_combined_stress_embodiment_a_star_lysis_pct": 5.5,
            "obs_combined_stress_comp_r4_lysis_pct": 5.0,
            "obs_combined_stress_comp_r4_gsis_si": 3.0,
        }
    )
    assert dt_cs_retain_missing_sec["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMBINED_STRESS_INCOMPLETE_SECONDARY_PANEL"

    # - R5 & R6 (MW-R5-S2-F05/F06 & MW-R5-S1-18): Normalize all new boolean/meta flags, verify Ser147/Ser353 and 0.4295/0.3395 lognormal fractions, and verify 50 um fluid-film gate
    norm_flags, _ = wc.normalize_measurements(
        {
            "obs_icasp9_log10_kill": 4.6,
            "obs_icasp9_log10_kill_se": 0.1,
            "track2_b_sct_or_hlag1_lilrb1_rescue_qualified": True,
            "cd47_playbook_200k_and_5pct_phago_waiver_approved": False,
            "s9_staged_cassette_excision_and_d10w_clamp_qualified": True,
            "periimplant_closed_suction_drain_verified": True,
            "obs_vascularized_explant_en_bloc_recovery_verified": True,
            "obs_ap1903_pulsed_pk_washout_regimen_verified": True,
            "mcb_clone_id": "MCB-CLONE-01",
            "oos_retest_investigation_approved": True,
            "oos_retest_assignable_cause_code": "PIPETTE_CALIBRATION_FAULT_CONFIRMED",
        }
    )
    assert norm_flags["mcb_clone_id"] == "MCB-CLONE-01"
    dt_film_65_no_drain = wetlab.evaluate_stage0_decision_tree({**core_14_panel, "obs_periimplant_fluid_film_um": 65.0})
    assert dt_film_65_no_drain["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S9_PERIIMPLANT_FLUID_FILM_GT_50UM_REQUIRES_CLOSED_SUCTION_DRAIN_AND_LOWLOAD_1500"
    dt_film_65_with_drain = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_periimplant_fluid_film_um": 65.0, "periimplant_closed_suction_drain_verified": True}
    )
    assert dt_film_65_with_drain["resolved_architecture_assignments"]["lot_maturation_status"] == "READY_FOR_ENCAPSULATION"
    assert dt_film_65_with_drain["resolved_architecture_assignments"]["pouch_areal_loading_cap_ieq_cm2"] == 1500.0
    isr_r5 = l2["stage0_to_stage7_perturbation_plane"]["coupled_five_transgene_hypoxic_isr_surface_decay_and_leader_rescue_audit_sonnet55_r5"]
    omega_r5 = l2["construct_sequence_and_grammar_verification"]["hlae_0101_vs_0103_gly107arg_and_cd47_daf_omega_site_sequence_audit_sonnet55_r5"]
    assert "Ser147" in omega_r5["cd47_igv_daf_175aa_omega_site_numbering_reconciliation"]["omega_cleavage_site_residue"]
    assert isr_r5["single_cell_lognormal_variegation_and_p10_floor_gate"]["fraction_of_cells_below_200k_at_bulk_mean_240k_sigma_0p45"] == 0.4285
    assert isr_r5["single_cell_lognormal_variegation_and_p10_floor_gate"]["fraction_of_cells_below_180k_at_bulk_mean_240k_sigma_0p45"] == 0.3393
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        assert l1["non_hla_neoepitope_and_complement_flow_crossmatch_audit_sonnet55_r5"]["chimeric_surface_junction_neoepitope_inventory_count"] == 306
        assert l1["non_hla_neoepitope_and_complement_flow_crossmatch_audit_sonnet55_r5"]["exact_boundary_constrained_8_to_11mer_window_count"] == 326

    # 15. AUDIT-REF [DEC-049]:
    # - AUDIT-REF: Uncleared S10 blocks S2 Arm 1X override, whereas cleared S10 coupon unmasks S2 Arm 1B failure
    dt_s10_blocks_s2 = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_nkg2c_adaptive_nk_lysis_pct": 12.0,
            "s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior": True,
        }
    )
    assert dt_s10_blocks_s2["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "RETAIN_EMBODIMENT_A_STAR_ARM1B"
    assert dt_s10_blocks_s2["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S10_NKG2C_SELECTIVITY_FAILED_REQUIRES_SCT_PEPTIDE_SELECTIVITY_COUPON_OR_CMV_STRATIFIED_LABEL"
    dt_s10_cleared_unmasks_s2 = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_hlae_hour0_molecules": 10800.0,
            "obs_hlae_nkg2c_over_nkg2a_kd_ratio": 3.5,
            "s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct": True,
        }
    )
    assert dt_s10_cleared_unmasks_s2["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "PROMOTE_ISLET1_A4_X_MOL_RNLS_WT"
    assert dt_s10_cleared_unmasks_s2["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_PENDING_ARM1X_LOCUS1_REPLACEMENT_QUALIFICATION"
    # - AUDIT-REF: 4-hazard competing-risks decomposition (85.7575% additive, 86.5048% multiplicative)
    pop_r5 = l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["b_prime_4digit_hlac_alloreactivity_and_b_prime_sct_pivot_audit_muse_r3"]["hazard_active_cohort_sensitivity_audit_sonnet55_r4"]["population_simulation_with_cmv_nkg2c_and_sensitization_sonnet55_r5"]
    assert pop_r5["four_hazard_competing_risks_decomposition_from_100pct"]["pre_clearance_unscreened_4hazard_additive_survival_pct"] == 85.7575
    assert pop_r5["pre_clearance_independent_competing_risks_multiplicative_survival_pct"] == 86.5048
    # - AUDIT-REF: HLA-E*01:01 (Arg107 in 3CDG) vs HLA-E*01:03 (Gly107 in pinned donors) and SERPINB9 t_half = 6.0h
    assert "Arg107Gly" in omega_r5["hlae_0101_gly107_vs_0103_arg107_audit"]["hlae_0103_canonical_substitution"]
    assert isr_r5["five_transgene_72h_hypoxic_isr_decay_trajectories_at_f_isr_0p20"]["SERPINB9_Locus3_BFB"]["t_half_hours"] == 6.0

    # 16. AUDIT-REF & AUDIT-REF [DEC-050]:
    # - R6-01 & S2-01: Catastrophic NKG2A > 15% and Regime II 2h lysis > 25% cannot be waived by S2/S1 flags, and Track-2 KIR escape cannot clear a Family-A product
    dt_cat_nkg2a_s2_override_blocked = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_nkg2a_lysis_pct": 19.0,
            "obs_hlae_hour0_molecules": 10000.0,
            "s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior": True,
        }
    )
    assert dt_cat_nkg2a_s2_override_blocked["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "HOLD_S2_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_REQUIRES_MULTI_RECEPTOR_RESCUE"
    assert dt_cat_nkg2a_s2_override_blocked["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_EXCEEDS_ARM1X_RESCUE_ENVELOPE"

    dt_cat_reg2_s1_override_blocked = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_regimeII_lysis_2h_pct": 28.0,
            "s1_l3_bfb_replacement_qualified": True,
        }
    )
    assert dt_cat_reg2_s1_override_blocked["switch_decisions"]["S1_Locus3_Kinetic_and_Topology_Gate_0C"] == "NO_GO_ACUTE_REGIME_II_LYSIS_GT_25PCT_TERMINAL"
    assert dt_cat_reg2_s1_override_blocked["resolved_architecture_assignments"]["lot_maturation_status"] == "NO_GO_ACUTE_REGIME_II_LYSIS_EXCEEDS_25PCT_CEILING"

    from scripts.product_contract import get_product_definition
    prod_a_star = get_product_definition("ISLET1_A4_ARM1B")
    dt_kir_family_a_mismatch = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_kir_nkg2a_neg_nk_lysis_pct": 14.0,
            "track2_b_sct_or_hlag1_lilrb1_rescue_qualified": True,
        },
        product_def=prod_a_star,
    )
    assert dt_kir_family_a_mismatch["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_TRACK1_MHC_NULL_KIR_ESCAPE_REQUIRES_TRACK2_B_SCT_OR_HLAG1_PRODUCT_SWITCH"

    # - R6-02/R6-06 & S1-02/S1-09: Rare-event denominators required for residual PSC, hazard-state survivor, and Step-1 survivor; vascularized explant flag wired
    dt_psc_no_denom = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_residual_psc_effective_fraction": 0.0}
    )
    assert dt_psc_no_denom["safety_kill_s9_status"] == "HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM"
    dt_vasc_explant_false = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_explant_en_bloc_recovery_fraction": 0.995,
            "obs_explant_en_bloc_recovery_se": 0.001,
            "obs_vascularized_explant_en_bloc_recovery_verified": False,
        }
    )
    assert dt_vasc_explant_false["safety_kill_s9_triggered"] is True
    assert "EXPLANT_RECOVERY" in dt_vasc_explant_false["s9_causal_failed_categories"]
    assert dt_vasc_explant_false["safety_kill_s9_status"] == "HARD_SAFETY_STOP_REJECTED"

    # - R6-03/S1-03/S2-04: One-sided 95% confidence bounds enforced on extended endpoints and COMP_R4 promotion requires all 4 COMP_R4 secondary metrics
    dt_20q_ci_breach = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_20q11_bcl2l1_cnv_ratio": 2.12, "obs_20q11_bcl2l1_cnv_se": 0.05}
    )
    assert dt_20q_ci_breach["safety_kill_s9_triggered"] is True
    assert "CNV_20Q11" in dt_20q_ci_breach["s9_causal_failed_categories"]
    dt_comp_r4_waived_missing_sec = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_combined_stress_embodiment_a_star_lysis_pct": 14.0,
            "obs_combined_stress_comp_r4_lysis_pct": 6.0,
            "obs_combined_stress_comp_r4_gsis_si": 3.0,
            "require_full_7metric_combined_stress_panel": False,
        }
    )
    assert dt_comp_r4_waived_missing_sec["pivotal_combined_stress_coupon_gate"]["evaluation_status"] == "HOLD_COMP_R4_PROMOTION_REQUIRES_COMPLETE_4METRIC_SECONDARY_PANEL"

    # - R6-12/S2-03/S2-05: Layer-5 autoantigen z >= 2.0 holds without gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared=True, and explicit crossmatch=False holds at any cPRA
    pt3_ref = next(
        p for p in l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["representative_patient_compilations"]
        if p["patient_id"].startswith("PT-003")
    )
    assert "HOLD_TYK2_INDEPENDENT_IFNG_TNFA_AND_INDIRECT_DR3_DQ8_CLEARANCE_REQUIRED" in pt3_ref["all_active_qualification_holds"]
    assert pt3_ref["in_silico_candidate_routing_cleared"] is False
    pt_low_cpra_failed_xmatch = {
        "patient_id": "PT-LOW-CPRA-XMATCH-FAIL",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.32,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_cpra_pct": 10.0,
        "exact_product_non_hla_and_sct_flow_crossmatch_negative": False,
    }
    bp_failed_xmatch = l5_mod.compile_personalized_patient_blueprint(pt_low_cpra_failed_xmatch)
    assert "HOLD_EXACT_PRODUCT_NON_HLA_AND_SCT_FLOW_CROSSMATCH_FAILED" in bp_failed_xmatch["all_active_qualification_holds"]
    assert l5_mod.evaluate_regimen_survival_probability(pt_low_cpra_failed_xmatch, bp_failed_xmatch) < 0.75

    # 17. AUDIT-REF [DEC-051]:
    # - S2-D1: Catastrophic NKG2A+ lysis upper95 > 15.0% blocks Arm 1X promotion even when nkg2a_explicit_se > 3.0%
    dt_cat_nkg2a_high_se = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_nkg2a_lysis_pct": 14.0,
            "obs_nkg2a_lysis_se": 4.0,  # upper95 = 20.579% > 15.0%
            "s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior": True,
        }
    )
    assert dt_cat_nkg2a_high_se["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "HOLD_S2_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_REQUIRES_MULTI_RECEPTOR_RESCUE"
    assert dt_cat_nkg2a_high_se["resolved_architecture_assignments"]["lot_maturation_status"] != "READY_FOR_ENCAPSULATION"

    # - S2-D2: Residual PSC fraction <= 1e-8 with explicit SE but missing obs_residual_psc_effective_trials still holds
    dt_psc_explicit_se_no_trials = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_residual_psc_effective_fraction": 0.0,
            "obs_residual_psc_effective_se": 1e-9,
        }
    )
    assert dt_psc_explicit_se_no_trials["safety_kill_s9_status"] == "HOLD_RESIDUAL_PSC_EFFECTIVE_TRIALS_BELOW_299573228_MINIMUM"

    # - S2-L2: Post-72h hypoxic + ISR single-cell p10 CD47 endpoint wired into Gate 0A Playbook cliff & exact lognormal p10 math in Layer 2
    dt_72h_p10_fail = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_hypoxic_72h_cd47_p10_single_cell_molecules": 150000.0,
            "obs_hypoxic_72h_cd47_p10_single_cell_se": 5000.0,
        }
    )
    assert dt_72h_p10_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED"
    p10_72h_audit = l2["stage0_to_stage7_perturbation_plane"]["coupled_five_transgene_hypoxic_isr_surface_decay_and_leader_rescue_audit_sonnet55_r5"]["single_cell_lognormal_variegation_and_p10_floor_gate"]
    assert p10_72h_audit["joint_72h_hypoxic_isr_single_cell_p10_cd47_audit_sol61_r6"]["stressed_72h_bulk_mean_183354_at_sigma_0p45_p10_molecules"] == 93081
    assert p10_72h_audit["joint_72h_hypoxic_isr_single_cell_p10_cd47_audit_sol61_r6"]["required_72h_stressed_bulk_mean_for_p10_ge_180k_at_sigma_0p45_molecules"] == 354570
    assert p10_72h_audit["mandatory_stage0_gate_0a_pass_criteria"]["obs_hypoxic_72h_cd47_p10_single_cell_molecules_min"] == 180000

    # - S2-L3 & S1-L5: Layer-3 Wooldridge 2005/2010 D227K/T228A + SA-THBD common dose d* & Layer-5 explicit KIR ligand group requirement for B_PRIME_SCT
    assert "PMID:15837791" in l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]["nkg2a_vs_nkg2c_2d_synaptic_selectivity_bound_sonnet55_r4"]["pure_nkg2a_neg_nkg2c_pos_subset_and_hlae_d227k_t228a_ctl_audit_sol61_r5"]["hlae_d227k_t228a_cd8_null_and_cmv_ul40_ctl_gate_sol61_r6"]["primary_citations"]
    assert l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4"]["non_hla_chimeric_junction_neoepitope_and_flow_crossmatch_audit_sonnet55_r5"]["sa_thbd_tafi_fibrin_scaffold_vs_ibmir_patency_gate_sol61_r6"]["common_selected_dose_d_star_required"] is True
    bp_b_sct_no_kir = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-BSCT-NOKIR",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.55,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 40.0,
            "sirpa_genotype": "v1/v1",
        },
        track2_subbank_mode="B_PRIME_SCT",
    )
    assert "HOLD_B_PRIME_SCT_REQUIRES_EXPLICIT_RECIPIENT_KIR_LIGAND_GROUP_TYPING" in bp_b_sct_no_kir["all_active_qualification_holds"]
    assert bp_b_sct_no_kir["lot_release_cleared"] is False

    # 18. AUDIT-REF [DEC-052]:
    # - R6-B01 & B02: Dose-level 8-PSC escape with 1% 20q11.21 mosaic is 8.0792e-6 (not double-counted 8.079e-14), and single-cell 20q11.21 mosaic > 1% triggers S9 SINGLE_CELL_MOSAIC_20Q11
    import math
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        tumor_r5 = l1["unconditional_patient_level_tumor_harm_and_two_stage_psc_purge_model_sonnet55_r5"]
        assert math.isclose(tumor_r5["computed_dose_level_8psc_escape_if_1pct_mosaic_20q11_pre_retrieval"], 8.0792e-6, rel_tol=1e-4)
        assert math.isclose(tumor_r5["computed_dose_level_8psc_harm_with_1pct_mosaic_20q11_subclone_bound"], 8.0892e-6, rel_tol=1e-4)
    dt_mosaic_fail = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_20q11_bcl2l1_cnv_ratio": 2.02,
            "obs_20q11_bcl2l1_single_cell_mosaic_fraction": 0.02,
        }
    )
    assert dt_mosaic_fail["safety_kill_s9_triggered"] is True
    assert "SINGLE_CELL_MOSAIC_20Q11" in dt_mosaic_fail["s9_causal_failed_categories"]

    # - R6-B03 & B04: Stressed 72h HLA-E p10 < 15,000 triggers Playbook cliff hold, and Layer-2 exposes 12,622 / 26,144 stressed HLA-E p10 math
    dt_72h_hlae_p10_fail = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_hypoxic_72h_hlae_p10_single_cell_molecules": 12622.0,
        }
    )
    assert dt_72h_hlae_p10_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED"
    assert p10_72h_audit["joint_72h_hypoxic_isr_single_cell_p10_cd47_audit_sol61_r6"]["hlae_arm1x_hour0_22000_at_f_isr_0p85_and_sigma_0p28_post_72h_p10_molecules"] == 12622
    assert p10_72h_audit["joint_72h_hypoxic_isr_single_cell_p10_cd47_audit_sol61_r6"]["required_hlae_hour0_bulk_mean_for_post_72h_f_isr_0p85_p10_ge_15k_at_sigma_0p28_molecules"] == 26144
    assert p10_72h_audit["mandatory_stage0_gate_0a_pass_criteria"]["obs_hypoxic_72h_hlae_p10_single_cell_molecules_min"] == 15000

    # - R6-D01: Default/imputed SE is retained in _upper95 / _lower95 (e.g., NKG2A 14% with defaulted SE=4.0% triggers catastrophic hold, and retrieval 1.0 with default SE=0.01 triggers EXPLANT_RECOVERY)
    dt_defaulted_nkg2a_cat = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_nkg2a_lysis_pct": 14.0,
            "obs_nkg2a_lysis_se": 4.0,
            "s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior": True,
        },
        defaulted_se_keys=["obs_nkg2a_lysis_se"],
    )
    assert dt_defaulted_nkg2a_cat["switch_decisions"]["S2_Locus1_Arm1B_vs_Arm1X_Gate_0A"] == "HOLD_S2_CATASTROPHIC_NKG2A_LYSIS_GT_15PCT_REQUIRES_MULTI_RECEPTOR_RESCUE"
    dt_default_se_retrieval = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_explant_en_bloc_recovery_fraction": 1.0,  # default SE = 0.01 -> lower95 = 0.98355 < 0.99
        }
    )
    assert dt_default_se_retrieval["safety_kill_s9_triggered"] is True
    assert "EXPLANT_RECOVERY" in dt_default_se_retrieval["s9_causal_failed_categories"]

    # - R6-D04: Explicit exact_product_non_hla_and_sct_flow_crossmatch_negative=False holds even on Track 3 autologous iPSC
    bp_track3_failed_xmatch = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-TRACK3-XMATCH-FAIL",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.32,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 40.0,
            "sirpa_genotype": "v1/v1",
            "requested_track": "TRACK_3_AUTOLOGOUS_IPSC",
            "exact_product_non_hla_and_sct_flow_crossmatch_negative": False,
        }
    )
    assert "HOLD_EXACT_PRODUCT_NON_HLA_AND_SCT_FLOW_CROSSMATCH_FAILED" in bp_track3_failed_xmatch["all_active_qualification_holds"]
    assert bp_track3_failed_xmatch["lot_release_cleared"] is False

    # - R6-B13 & B16: Layer-4 coupled chronic BVP + spherical cluster + seroma film table and portal lag x clearance exercise envelope
    film_limit = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["bed_level_convective_o2_budget_and_periimplant_fluid_film_audit_sonnet55_r5"]["periimplant_stagnant_fluid_film_transport_limit"]
    assert film_limit["coupled_bvp_and_series_fluid_film_core_po2_table_r6"]["live_chronic_bvp_and_spherical_cluster_series_film_table_astra6_r6"]["cluster_100um_core_po2_mmhg_by_film_um"]["80"] == 4.71
    assert film_limit["coupled_bvp_and_series_fluid_film_core_po2_table_r6"]["live_chronic_bvp_and_spherical_cluster_series_film_table_astra6_r6"]["cluster_80um_dxb_core_po2_mmhg_by_film_um"]["80"] == 8.95
    ex_env = l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]["challenged_v0p60_300k_beta_eq_exercise_margin_audit_sol61_r5"]["portal_drainage_lag_by_clearance_exercise_envelope_astra6_r6"]
    assert ex_env["lag_5p17min_tau_12min_min_glucose_mM"] == 2.7160

    # 19. AUDIT-REF [DEC-053] & AUDIT-REF [DEC-054] (Meta Audit Full AUDIT-REF):
    # - S2-R1: Acute Regime II 2h lysis > 25% sets S1_Locus3_Kinetic_and_Topology_Gate_0C to NO_GO_ACUTE_REGIME_II_LYSIS_GT_25PCT_TERMINAL
    assert dt_cat_reg2_s1_override_blocked["switch_decisions"]["S1_Locus3_Kinetic_and_Topology_Gate_0C"] == "NO_GO_ACUTE_REGIME_II_LYSIS_GT_25PCT_TERMINAL"
    assert dt_cat_reg2_s1_override_blocked["resolved_architecture_assignments"]["locus3_construct"] == "ISLET1_A4_ARM1B"

    # - S2-R2: Sorted KIR+ NKG2A- NK or IL-12/IL-18 SIRPa-low lysis failure uses Track-2 B'-SCT/HLA-G1 hold string (no residual DX-S/L3-BFB reference)
    assert dt_kir_neg_nk["resolved_architecture_assignments"]["locus1_construct"] == "HOLD_KIR_NKG2A_NEG_NK_OR_IL12_IL18_SIRPA_LOW_LYSIS_FAILED_REQUIRES_TRACK2_B_SCT_OR_HLAG1_CHALLENGER"
    assert dt_kir_neg_nk["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_KIR_NKG2A_NEG_NK_OR_IL12_IL18_SIRPA_LOW_LYSIS_EXCEEDED_REQUIRES_TRACK2_B_SCT_OR_HLAG1_CHALLENGER"

    # - S2-R3: Layer-3 HLA-A2 D227K/T228A positive control SPR required, donor_level_pass_criterion, and treg_primary_estimand_hierarchy exposed
    d227k_gate = l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4"]["nkg2a_vs_nkg2c_2d_synaptic_selectivity_bound_sonnet55_r4"]["pure_nkg2a_neg_nkg2c_pos_subset_and_hlae_d227k_t228a_ctl_audit_sol61_r5"]["hlae_d227k_t228a_cd8_null_and_cmv_ul40_ctl_gate_sol61_r6"]
    assert d227k_gate["hla_a2_d227k_t228a_positive_control_spr_required"] is True
    assert "donor_level_pass_criterion" in d227k_gate
    fasl_gate = l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4"]["non_hla_chimeric_junction_neoepitope_and_flow_crossmatch_audit_sonnet55_r5"]["sa_fasl_microgel_nkp46_innate_nk_neutrality_and_teff_selectivity_gate_sol61_r6"]
    assert "treg_primary_estimand_hierarchy" in fasl_gate
    assert "membrane_complement_regulator_cd55_scr14_cd59_and_assembled_cassette_ibmir_gate_muse_r6" in l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4"]["non_hla_chimeric_junction_neoepitope_and_flow_crossmatch_audit_sonnet55_r5"]

    # - S2-R4: Live chronic BVP + spherical cluster series film table computed dynamically and matches all 16 live solver values
    live_film_tbl = film_limit["coupled_bvp_and_series_fluid_film_core_po2_table_r6"]["live_chronic_bvp_and_spherical_cluster_series_film_table_astra6_r6"]
    assert live_film_tbl["computation_mode"] == "LIVE_NZ201_CHRONIC_BVP_AND_SPHERICAL_CLUSTER_SOLVER_MW_R6_S2_R4"
    assert live_film_tbl["cluster_100um_surface_po2_mmhg_by_film_um"] == {"0": 16.52, "35": 14.70, "50": 14.10, "80": 13.09}
    assert live_film_tbl["cluster_100um_core_po2_mmhg_by_film_um"] == {"0": 7.96, "35": 6.23, "50": 5.65, "80": 4.71}
    assert live_film_tbl["cluster_80um_dxb_surface_po2_mmhg_by_film_um"] == {"0": 17.89, "35": 16.06, "50": 15.45, "80": 14.43}
    assert live_film_tbl["cluster_80um_dxb_core_po2_mmhg_by_film_um"] == {"0": 12.35, "35": 10.55, "50": 9.95, "80": 8.95}

    # - S2-R5: Explicit Track 2 B_ENDORETAINED_HLAC without kir_ligand_group triggers HOLD_B_PRIME_SCT_REQUIRES_EXPLICIT_RECIPIENT_KIR_LIGAND_GROUP_TYPING
    bp_b2_no_kir = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-B2-NOKIR",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.55,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 40.0,
            "sirpa_genotype": "v1/v1",
        },
        track2_subbank_mode="B_ENDORETAINED_HLAC",
    )
    assert "HOLD_B_PRIME_SCT_REQUIRES_EXPLICIT_RECIPIENT_KIR_LIGAND_GROUP_TYPING" in bp_b2_no_kir["all_active_qualification_holds"]
    assert bp_b2_no_kir["lot_release_cleared"] is False

    # - Stage-1 Audit Full R6 gates: Gate 0V viral persistence, post-stress 20q11 outgrowth S9, Nutlin-3a p53 S9, Chandler-loop IBMIR, efferocytosis clearance, and lineage drift
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        assert math.isclose(tumor_r5["computed_dose_level_8psc_harm_with_1pct_mosaic_20q11_pre_retrieval_bound"], 9.0792e-6, rel_tol=1e-4)
        assert math.isclose(tumor_r5["computed_post_retrieval_1pct_residual_8psc_harm_with_1pct_mosaic_20q11_including_1e8_metastasis_floor"], 9.0792e-8, rel_tol=1e-4)
        assert math.isclose(tumor_r5["computed_post_retrieval_1pct_residual_8psc_harm_with_1pct_mosaic_20q11_with_retrieval_attenuated_metastasis"], 8.0892e-8, rel_tol=1e-4)
    dt_gate0v_fail = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_gate0v_cvb_day14_log10_viral_reduction": 1.2}
    )
    assert dt_gate0v_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_GATE0V_VIRAL_SANCTUARY_PERSISTENCE_FAILED"
    dt_nutlin_fail = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_nutlin3a_p21_mdm2_fold_induction": 2.1}
    )
    assert dt_nutlin_fail["safety_kill_s9_triggered"] is True
    assert "FUNCTIONAL_P53_NUTLIN3A" in dt_nutlin_fail["s9_causal_failed_categories"]
    dt_outgrowth_null_pass = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_post_stress_20q11_bcl2l1_outgrowth_fold": 1.0, "obs_post_stress_20q11_bcl2l1_outgrowth_se": 0.01}
    )
    assert "POST_STRESS_20Q11_OUTGROWTH" not in dt_outgrowth_null_pass["s9_causal_failed_categories"]
    dt_outgrowth_fail = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_post_stress_20q11_bcl2l1_outgrowth_fold": 1.8}
    )
    assert dt_outgrowth_fail["safety_kill_s9_triggered"] is True
    assert "POST_STRESS_20Q11_OUTGROWTH" in dt_outgrowth_fail["s9_causal_failed_categories"]
    dt_chandler_fail = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_assembled_cassette_chandler_loop_sc5b9_fold": 2.4}
    )
    assert dt_chandler_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_ASSEMBLED_CASSETTE_CHANDLER_LOOP_IBMIR_OR_COMPLEMENT_FAILED"
    dt_effero_fail = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_apoptotic_effero_clearance_retention_pct_of_wt": 55.0}
    )
    assert dt_effero_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_APOPTOTIC_EFFEROCYTOSIS_CLEARANCE_IMPAIRED_DAMPS_RISK"
    dt_lineage_fail = wetlab.evaluate_stage0_decision_tree(
        {**core_14_panel, "obs_assembled_cluster_slc18a1_enterochromaffin_pct": 4.5}
    )
    assert dt_lineage_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_ASSEMBLED_CLUSTER_OFF_TARGET_LINEAGE_OR_BETA_IDENTITY_DRIFT_FAILED"
    assert dt_lineage_fail["resolved_architecture_assignments"]["prescribed_functional_dose_cap_ieq"] == 150000.0

    # 20. AUDIT-REF & AUDIT-REF [DEC-055 & DEC-056]:
    # - R7-F01: Leave-one-out / None / UNAVAILABLE on any of the 13 DEC-053 endpoints blocks READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE
    full_ext_payload = dict(full_ext_input)
    dec053_required_13 = [
        "obs_gate0v_cvb_day14_log10_viral_reduction",
        "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct",
        "obs_post_stress_20q11_bcl2l1_outgrowth_fold",
        "obs_assembled_cassette_chandler_loop_sc5b9_fold",
        "obs_assembled_cassette_chandler_loop_platelet_retention_pct",
        "obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction",
        "obs_apoptotic_effero_clearance_retention_pct_of_wt",
        "obs_nutlin3a_p21_mdm2_fold_induction",
        "obs_vascular_pericyte_ng2_pdgfrb_coverage_pct",
        "obs_assembled_cluster_slc18a1_enterochromaffin_pct",
        "obs_assembled_cluster_ins_gcg_polyhormonal_pct",
        "obs_assembled_cluster_ki67_progenitor_pct",
        "obs_post_stress_mafa_nkx61_beta_identity_retention_pct",
    ]
    for k_req in dec053_required_13:
        loo_payload = {k: v for k, v in full_ext_payload.items() if k != k_req}
        loo_res = wetlab.evaluate_stage0_decision_tree(loo_payload)
        assert loo_res["lot_release_readiness_tier"] == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
        unavail_res = wetlab.evaluate_stage0_decision_tree({**full_ext_payload, k_req: "UNAVAILABLE"})
        assert unavail_res["lot_release_readiness_tier"] == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"

    # - R7-F02: Film 50..80 um with verified drain triggers S9_Chronic, enforces <= 100 um (or 80 um) cluster diameter cap and 1500 IEQ/cm2 loading cap, and returns CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION
    dt_film_79p5_drain = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_payload,
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "obs_periimplant_fluid_film_um": 79.5,
            "obs_periimplant_fluid_film_se": 0.01,
            "periimplant_closed_suction_drain_verified": True,
        }
    )
    assert dt_film_79p5_drain["resolved_architecture_assignments"]["pouch_cluster_diameter_cap_um"] == 100.0
    assert dt_film_79p5_drain["resolved_architecture_assignments"]["pouch_areal_loading_cap_ieq_cm2"] == 1500.0
    assert dt_film_79p5_drain["switch_decisions"]["S9_Chronic_Fibrotic_Capsule_and_Mineralization_Gate_0F"] == "ENFORCE_80UM_CLUSTER_CAP_AND_CITRATE_CHELATION"
    assert dt_film_79p5_drain["lot_release_readiness_tier"] == "CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION"

    # 21. AUDIT-REF [DEC-057]:
    # - R7S2-F01: Supplementary aliases (obs_nutlin3a_p21_cdkn1a_fold_induction, obs_enterochromaffin_fev_tph1_pos_pct) cannot substitute for omitted/None/UNAVAILABLE canonical DEC-053 keys
    alias_sub_nutlin = {
        k: v for k, v in full_ext_payload.items() if k != "obs_nutlin3a_p21_mdm2_fold_induction"
    }
    alias_sub_nutlin["obs_nutlin3a_p21_cdkn1a_fold_induction"] = 8.4
    assert (
        wetlab.evaluate_stage0_decision_tree(alias_sub_nutlin)["lot_release_readiness_tier"]
        == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
    )
    assert (
        wetlab.evaluate_stage0_decision_tree(
            {**alias_sub_nutlin, "obs_nutlin3a_p21_mdm2_fold_induction": "UNAVAILABLE"}
        )["lot_release_readiness_tier"]
        == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
    )
    alias_sub_ec = {
        k: v for k, v in full_ext_payload.items() if k != "obs_assembled_cluster_slc18a1_enterochromaffin_pct"
    }
    alias_sub_ec["obs_enterochromaffin_fev_tph1_pos_pct"] = 0.3
    assert (
        wetlab.evaluate_stage0_decision_tree(alias_sub_ec)["lot_release_readiness_tier"]
        == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
    )
    assert (
        wetlab.evaluate_stage0_decision_tree(
            {**alias_sub_ec, "obs_assembled_cluster_slc18a1_enterochromaffin_pct": "UNAVAILABLE"}
        )["lot_release_readiness_tier"]
        == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
    )

    # - R7S2-F01: obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction enforces >= 0.90 (and lower95 >= 0.90)
    dt_joint_805 = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_payload,
            "obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction": 0.805,
            "obs_hypoxic_72h_five_transgene_joint_above_threshold_se": 0.001,
        }
    )
    assert (
        dt_joint_805["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED"
    )
    assert dt_joint_805["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    dt_joint_905_wide_se = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_payload,
            "obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction": 0.905,
            "obs_hypoxic_72h_five_transgene_joint_above_threshold_se": 0.01,  # lower95 = 0.88855 < 0.90
        }
    )
    assert dt_joint_905_wide_se["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    assert (
        l2["stage0_to_stage7_perturbation_plane"][
            "five_transgene_joint_coverage_serpinb9_stoichiometry_and_lineage_drift_audit_muse_r6"
        ]["joint_five_transgene_72h_hypoxic_isr_single_cell_coverage_gate"][
            "min_joint_five_transgene_coverage_fraction"
        ]
        == 0.90
    )

    # - R7S2-F02 & Q1/Q2: Historical ledger strict boolean validation (27 persisted boolean flags, 31 ledger metadata booleans) & Layer-4 5.71 mmHg prose alignment
    assert len(wetlab.STAGE0_PERSISTED_BOOLEAN_FLAGS) == 27
    assert "require_full_5endpoint_s9_panel" in wetlab.STAGE0_PERSISTED_BOOLEAN_FLAGS
    assert "require_full_7metric_combined_stress_panel" in wetlab.STAGE0_PERSISTED_BOOLEAN_FLAGS
    assert len(wetlab.HISTORICAL_LEDGER_BOOLEAN_METADATA_KEYS) == 31
    assert "5.71 mmHg at 100 um clusters with matched D_caps = D_gel" in film_limit["surgical_mitigation"]

    batch1_r7s2 = {
        **full_ext_payload,
        "batch_id": "BATCH-R7S2-01",
        "mcb_clone_id": "CLONE-R7S2",
        "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
        "obs_icasp9_log10_kill_se": 0.05,
        "obs_vascularized_explant_en_bloc_recovery_verified": False,
    }
    rep1_r7s2 = wetlab.run_calibration(batch1_r7s2, is_dry_run=False)
    assert rep1_r7s2["stage0_decision_tree_evaluation"]["safety_kill_s9_triggered"] is True
    for bad_bool_val in ("false", "False", 0):
        tampered_rep1 = json.loads(json.dumps(rep1_r7s2))
        tampered_rep1["ordered_batch_observation_ledger"][0][
            "obs_vascularized_explant_en_bloc_recovery_verified"
        ] = bad_bool_val
        tampered_rep1["stage0_decision_tree_evaluation"]["safety_kill_s9_triggered"] = False
        tampered_rep1["stage0_decision_tree_evaluation"][
            "resolved_architecture_assignments"
        ]["variant_safety_status"] = "CLEARED"
        tampered_rep1["calibration_state_sha256"] = wetlab._calibration_state_sha256(
            tampered_rep1
        )
        batch2_r7s2 = {
            "batch_id": "BATCH-R7S2-02",
            "mcb_clone_id": "CLONE-R7S2",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "obs_stage7_yield_ratio": 0.92,
            "obs_stage7_yield_se": 0.04,
        }
        try:
            wetlab.run_calibration(batch2_r7s2, False, tampered_rep1)
            raise AssertionError(f"Expected ValueError for bad historical boolean {bad_bool_val!r}")
        except ValueError as exc:
            assert "obs_vascularized_explant_en_bloc_recovery_verified" in str(exc)

    # =========================================================================
    # 22. AUDIT-REF [DEC-058] STAGE-0 COMPLETE-EVIDENCE ORIGINAL PRIMARY
    #     ENDPOINTS (R7S3-F01), DIRECT EVALUATOR BOOLEAN VALIDATION (R7S3-F02),
    #     AND BOUNDED QUALIFICATIONS (Q1..Q5)
    # =========================================================================
    # - R7S3-F01: Leave-one-out / None / UNAVAILABLE on any of the 8 original primary biological & pouch endpoints blocks READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE
    primary_8_keys = [
        ("obs_stage7_yield_ratio", "obs_stage7_yield_se"),
        ("obs_cd47_fold", "obs_cd47_fold_se"),
        ("obs_serpinb9_log10_kass", "obs_serpinb9_log10_kass_se"),
        ("obs_icasp9_log10_kill", "obs_icasp9_log10_kill_se"),
        ("obs_pouch_min_po2_mmhg", "obs_pouch_min_po2_se"),
        ("obs_pouch_max_h2o2_uM", "obs_pouch_max_h2o2_se"),
        ("obs_pouch_max_ph", "obs_pouch_max_ph_se"),
        ("obs_ifng_hla_dr_dq_negative_pct", "obs_ifng_hla_dr_dq_negative_se"),
    ]
    for mk, sk in primary_8_keys:
        omitted_p = {k: v for k, v in full_ext_payload.items() if k not in (mk, sk)}
        dt_omitted = wetlab.evaluate_stage0_decision_tree(omitted_p)
        assert dt_omitted["has_explicit_original_primary_obs"] is False, f"Expected False for omitted {mk}"
        assert mk in dt_omitted["missing_original_primary_endpoints"]
        assert dt_omitted["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

        none_p = {**omitted_p, mk: None}
        dt_none = wetlab.evaluate_stage0_decision_tree(none_p)
        assert dt_none["has_explicit_original_primary_obs"] is False, f"Expected False for None {mk}"
        assert dt_none["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

        unavail_p = {**omitted_p, mk: "UNAVAILABLE"}
        dt_unavail = wetlab.evaluate_stage0_decision_tree(unavail_p)
        assert dt_unavail["has_explicit_original_primary_obs"] is False, f"Expected False for UNAVAILABLE {mk}"
        assert dt_unavail["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

        if mk != "obs_icasp9_log10_kill":
            cal_omitted = wetlab.run_calibration(
                {**omitted_p, "batch_id": f"BATCH-R7S3-OMIT-{mk}"}, is_dry_run=False
            )
            assert (
                cal_omitted["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
                != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
            )

    # - R7S3-F02: Both wetlab_contract.evaluate_stage0_decision_tree and ingest_wetlab_results.evaluate_stage0_decision_tree reject non-boolean policy values
    from scripts import wetlab_contract as wc_mod

    assert len(wc_mod.STAGE0_DECISION_TREE_BOOLEAN_KEYS) == 31
    for bkey in wc_mod.STAGE0_DECISION_TREE_BOOLEAN_KEYS:
        for bad_val in ("false", "False", 0, "true", "True", 1):
            for eval_fn in (wc_mod.evaluate_stage0_decision_tree, wetlab.evaluate_stage0_decision_tree):
                try:
                    eval_fn({**full_ext_payload, bkey: bad_val})
                    raise AssertionError(f"Expected ValueError for {bkey}={bad_val!r} in {eval_fn}")
                except ValueError as exc:
                    assert bkey in str(exc)

    # - Qualifications Q2..Q4: Gate 0V PMIDs, Layer-5 UNSPECIFIED 4-digit relation & Track 3 crossmatch measurement, and keyword require_full_5endpoint_s9_panel persistence
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        assert tumor_r5["antiviral_sanctuary_and_cell_intrinsic_viral_fuse_gate_0v_muse_r6"]["primary_citations"] == [
            "PMID:19266182",
            "PMID:23974921",
        ]
    pt_c1_short_unspec = {
        **pt_c1_untyped,
        "patient_id": "PT-C1-SHORT-UNSPEC",
        "hlac_4digit_relation": "UNSPECIFIED",
    }
    comp_c1_short = l5_mod.compile_personalized_patient_blueprint(pt_c1_short_unspec, track2_subbank_mode="B_PRIME_SCT")
    assert (
        "HOLD_B_PRIME_SCT_C1_UNSPECIFIED_4DIGIT_RELATION_REQUIRES_ALLELE_TYPING_AND_CD8_MLC_CLEARANCE"
        in comp_c1_short["all_active_qualification_holds"]
    )
    pt_track3_no_xm = {
        "patient_id": "PT-TRACK3-NO-XM",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.25,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 35.0,
        "sirpa_genotype": "v1/v1",
        "cmv_nkg2c_adaptive_nk_fraction": 0.04,
        "recipient_bw4_kir3dl1_educated": False,
        "recipient_anti_hla_e_dsa": False,
        "requested_track": "TRACK_3_AUTOLOGOUS_IPSC",
        "wetlab_stage0_complete_lot_release_verified": True,
    }
    comp_t3_no_xm = l5_mod.compile_personalized_patient_blueprint(pt_track3_no_xm)
    assert comp_t3_no_xm["selected_fda_deployment_track"] == "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC"
    assert comp_t3_no_xm["recipient_immune_hazard_panel_explicitly_measured"] is False
    comp_t3_with_xm = l5_mod.compile_personalized_patient_blueprint(
        {
            **pt_track3_no_xm,
            "exact_product_non_hla_and_sct_flow_crossmatch_negative": True,
        }
    )
    assert comp_t3_with_xm["recipient_immune_hazard_panel_explicitly_measured"] is True

    rep_kw_5ep = wetlab.run_calibration(
        {
            "batch_id": "BATCH-R7S3-KW5EP",
            "obs_stage7_yield_ratio": 0.95,
            "obs_stage7_yield_se": 0.02,
        },
        is_dry_run=False,
        require_full_5endpoint_s9_panel=True,
    )
    assert (
        rep_kw_5ep["ordered_batch_observation_ledger"][0].get("require_full_5endpoint_s9_panel")
        is True
    )

    # 23. AUDIT-REF [DEC-059]:
    # - L2-Q1: exact 6-decimal Bernoulli breakthrough literal 0.026922
    assert (
        '"basal_0p30uM_start_immediate_licensed_hypoxic_replenishment_0p4005uM_prob_breakthrough": 0.026922'
        in (ROOT / "scripts/run_layer2_perturbation_and_regulatory_plane.py").read_text(encoding="utf-8")
    )
    # - ING-Q1: reject unknown top-level keys in prior stage0_decision_tree_evaluation even after recomputing digest
    import copy as _copy_r7s4

    tampered_eval_rep = _copy_r7s4.deepcopy(rep_kw_5ep)
    tampered_eval_rep["stage0_decision_tree_evaluation"]["unrecognized_policy"] = True
    tampered_eval_rep["calibration_state_sha256"] = wetlab._calibration_state_sha256(tampered_eval_rep)
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-R7S4-UNKNOWN-EVAL",
                "obs_stage7_yield_ratio": 0.96,
                "obs_stage7_yield_se": 0.02,
            },
            is_dry_run=False,
            prior_posterior_report=tampered_eval_rep,
        )
        raise AssertionError("Expected ValueError for unrecognized key in prior stage0_decision_tree_evaluation")
    except ValueError as exc:
        assert "unrecognized_policy" in str(exc)

    # - L4-Q1: direct 2-decimal and 4-decimal 100-um core pO2 exposed alongside two-step rounded table
    assert live_film_tbl["cluster_100um_core_po2_4dec_mmhg_by_film_um"]["0"] == 7.965
    assert live_film_tbl["cluster_100um_core_po2_direct_2dec_mmhg_by_film_um"]["0"] == 7.97

    # 24. AUDIT-REF [DEC-060]:
    l3_txt_r7s5 = (ROOT / "scripts/run_layer3_structural_binder_plane.py").read_text(encoding="utf-8")
    assert '"apoptotic_beta_cell_efferocytosis_clearance_gate_muse_r6"' in l3_txt_r7s5
    assert "MerTK/AXL/TIM-4" in l3_txt_r7s5
    assert "annexin_v_positive_apoptotic_corpse_mertk_axl_tim4_efferocytosis_clearance_percentage_of_wt" in (
        ROOT / "scripts/wetlab_contract.py"
    ).read_text(encoding="utf-8")

    # 25. AUDIT-REF [DEC-061] — R7-AUDIT-REF:
    # cd47_playbook_200k_and_5pct_phago_waiver_approved=True must NEVER waive mandatory post-stress/hypoxic/single-cell p10 or retention gates
    _se_map_r7s6 = {row[1]: row[2] for row in wc.STAGE0_DECISION_TREE_ENDPOINT_SPECS}
    for _mk_r7s6, _val_r7s6 in (
        ("obs_cd47_p10_single_cell_molecules", 180050.0),
        ("obs_hypoxic_72h_cd47_p10_single_cell_molecules", 180050.0),
        ("obs_hypoxic_72h_hlae_p10_single_cell_molecules", 15050.0),
        ("obs_hypoxic_72h_cd47_retention_molecules", 180050.0),
        ("obs_hypoxic_72h_hlae_retention_molecules", 15050.0),
    ):
        for _waiver_r7s6 in (False, True):
            _payload_r7s6 = {
                **full_ext_input,
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                "mcb_clone_id": "CLONE-AUDIT-REF",
                "batch_id": f"BATCH-{_mk_r7s6}-{_waiver_r7s6}",
                "cd47_playbook_200k_and_5pct_phago_waiver_approved": _waiver_r7s6,
                "obs_hypoxic_72h_five_transgene_joint_above_threshold_fraction": 0.92,
                "obs_hypoxic_72h_five_transgene_joint_above_threshold_se": 0.005,
                _mk_r7s6: _val_r7s6,
                _se_map_r7s6[_mk_r7s6]: 100.0,
            }
            _dt_r7s6 = wetlab.evaluate_stage0_decision_tree(_payload_r7s6)
            assert _dt_r7s6["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
            assert (
                _dt_r7s6["resolved_architecture_assignments"]["lot_maturation_status"]
                == "HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED"
            )
            assert (
                _dt_r7s6["resolved_architecture_assignments"]["locus1_construct"]
                == "HOLD_PLAYBOOK_CD47_200K_OR_PHAGO_5PCT_OR_HYPOXIC_RETENTION_CLIFF_FAILED"
            )
            _rep1_r7s6 = wetlab.run_calibration(_payload_r7s6, is_dry_run=False)
            assert (
                _rep1_r7s6["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
                == "HOLD_FAILED_OR_INCOMPLETE_GATE"
            )
            _rep2_r7s6 = wetlab.run_calibration(
                {
                    "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                    "mcb_clone_id": "CLONE-AUDIT-REF",
                    "batch_id": f"BATCH-NEXT-{_mk_r7s6}-{_waiver_r7s6}",
                    "obs_stage7_yield_ratio": 0.96,
                    "obs_stage7_yield_se": 0.01,
                },
                is_dry_run=False,
                prior_posterior_report=_rep1_r7s6,
            )
            assert (
                _rep2_r7s6["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
                == "HOLD_FAILED_OR_INCOMPLETE_GATE"
            )
    # Soft resting CD47 185k (>= 50k and < 200k) with all hypoxic/p10 floors passing IS waived when waiver=True:
    _soft_cd47_no_waiver = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_input,
            "obs_cd47_hour0_molecules": 185000.0,
            _se_map_r7s6["obs_cd47_hour0_molecules"]: 100.0,
            "cd47_playbook_200k_and_5pct_phago_waiver_approved": False,
        }
    )
    assert _soft_cd47_no_waiver["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    _soft_cd47_with_waiver = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_input,
            "obs_cd47_hour0_molecules": 185000.0,
            _se_map_r7s6["obs_cd47_hour0_molecules"]: 100.0,
            "cd47_playbook_200k_and_5pct_phago_waiver_approved": True,
        }
    )
    assert _soft_cd47_with_waiver["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    assert _soft_cd47_with_waiver["resolved_architecture_assignments"]["locus1_construct"] == "ISLET1_A4_ARM1B"

    # 26. AUDIT-REF [DEC-062] — R7-AUDIT-REF:
    # Scale-aware conservative default SEs + explicit-SE requirement on p10/retention molecule-count endpoints prevent omitted-SE bypass
    _contract_r7s7 = wc.derive_contract(ROOT)
    for _ep_r7s7 in _contract_r7s7["stage0_decision_tree_endpoints"]:
        if _ep_r7s7["max"] - _ep_r7s7["min"] >= 100000.0 and "molecules_per_cell" in _ep_r7s7["unit"].lower():
            assert _ep_r7s7["default_se"] >= 150.0, f"Expected scale-aware default_se >= 150.0 for {_ep_r7s7['meanKey']}"
            assert _ep_r7s7["minimum_se"] >= 1.0, f"Expected minimum_se >= 1.0 for {_ep_r7s7['meanKey']}"
    assert wc.STAGE0_CUSTOM_DEFAULT_SE["obs_serpinb9_kass_M_s"] >= 50000.0

    for _mk_r7s7, _sk_r7s7, _borderline_val_r7s7, _high_val_r7s7 in (
        ("obs_cd47_p10_single_cell_molecules", "obs_cd47_p10_single_cell_se", 180500.0, 210000.0),
        ("obs_hypoxic_72h_cd47_p10_single_cell_molecules", "obs_hypoxic_72h_cd47_p10_single_cell_se", 180500.0, 210000.0),
        ("obs_hypoxic_72h_hlae_p10_single_cell_molecules", "obs_hypoxic_72h_hlae_p10_single_cell_se", 15050.0, 25000.0),
        ("obs_hypoxic_72h_cd47_retention_molecules", "obs_hypoxic_72h_cd47_retention_se", 180500.0, 210000.0),
        ("obs_hypoxic_72h_hlae_retention_molecules", "obs_hypoxic_72h_hlae_retention_se", 15050.0, 25000.0),
        ("obs_iec_cd47_p10_molecules", "obs_iec_cd47_p10_se", 180500.0, 210000.0),
    ):
        # Leave-one-SE-out at borderline mean must fail _lower95 via conservative default SE and return HOLD_FAILED_OR_INCOMPLETE_GATE:
        _borderline_omit_se = {k: v for k, v in full_ext_input.items() if k != _sk_r7s7}
        _borderline_omit_se[_mk_r7s7] = _borderline_val_r7s7
        _dt_borderline_omit = wetlab.evaluate_stage0_decision_tree(_borderline_omit_se)
        assert _dt_borderline_omit["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE", (
            f"Expected HOLD_FAILED_OR_INCOMPLETE_GATE for omitted {_sk_r7s7} at borderline {_mk_r7s7}={_borderline_val_r7s7}"
        )
        _cal_borderline_omit = wetlab.run_calibration(
            {
                **_borderline_omit_se,
                "batch_id": f"BATCH-R7S7-BORDER-{_mk_r7s7}",
                "mcb_clone_id": "CLONE-R7S7",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            },
            is_dry_run=False,
        )
        assert (
            _cal_borderline_omit["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
            == "HOLD_FAILED_OR_INCOMPLETE_GATE"
        )
        # Leave-one-SE-out even at high passing mean must downgrade from COMPLETE to PROVISIONAL (because se_provenance != EXPLICIT_SE_SUPPLIED):
        _high_omit_se = {k: v for k, v in full_ext_input.items() if k != _sk_r7s7}
        _high_omit_se[_mk_r7s7] = _high_val_r7s7
        _dt_high_omit = wetlab.evaluate_stage0_decision_tree(_high_omit_se)
        assert (
            _dt_high_omit["lot_release_readiness_tier"]
            == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
        )

    # Also verify resting hour0 CD47 (200,500) and HLA-E (15,050) with omitted SE hold at borderline means:
    assert (
        wetlab.evaluate_stage0_decision_tree(
            {
                **{k: v for k, v in full_ext_input.items() if k != "obs_cd47_hour0_molecules_se"},
                "obs_cd47_hour0_molecules": 200500.0,
            }
        )["lot_release_readiness_tier"]
        == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    )
    assert (
        wetlab.evaluate_stage0_decision_tree(
            {
                **{k: v for k, v in full_ext_input.items() if k != "obs_hlae_hour0_molecules_se"},
                "obs_hlae_hour0_molecules": 15050.0,
            }
        )["lot_release_readiness_tier"]
        == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    )
    # Also verify that omitting any of the 4 EXTENDED_STAGE0_STAGE1_SPECS SE keys in run_calibration() propagates DEFAULT_SE_IMPUTED_NOT_EXPLICIT and blocks COMPLETE readiness (DEC-062 / DEC-063):
    for _ext_se_key in ("obs_pouch_min_po2_se", "obs_pouch_max_h2o2_se", "obs_pouch_max_ph_se", "obs_ifng_hla_dr_dq_negative_se"):
        _cal_ext_se_omit = wetlab.run_calibration(
            {
                **{k: v for k, v in full_ext_input.items() if k != _ext_se_key},
                "batch_id": f"BATCH-R7S8-EXT-SE-OMIT-{_ext_se_key}",
                "mcb_clone_id": "MCB-CLONE-R7S8-EXT-SE",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            },
            is_dry_run=False,
        )
        assert (
            _cal_ext_se_omit["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
            == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
        )
        assert (
            _cal_ext_se_omit["stage0_decision_tree_evaluation"]["has_explicit_original_primary_obs"]
            is False
        )

    # 27. AUDIT-REF [DEC-064] — R8-F01 (Layer-4 Transient Low-Load Density vs Intrinsic OCR Separation) & R8-F02 (Resting Molecule & Rate Explicit-SE Enforcement):
    from scripts import run_layer4_tissue_and_metabolic_twin as l4_mod

    _r8_seroma50_pb30 = l4_mod._run_single_pouch_pde_regime(
        "R8_SEROMA50_PB30_LOWLOAD_1500",
        0.45 * 0.60,
        "central",
        68.0,
        True,
        0.97,
        boundary_po2_mmhg=30.0,
        demand_multiplier_extra=0.60,
        areal_loading_ieq_cm2=1500.0,
        cell_intrinsic_ocr_multiplier=1.0,
        prescribed_cluster_diameter_um=100.0,
        nz_override=201,
        demand_spatial_mode="outer_annulus",
        lamina_d_o2_scale=0.42 ** 1.5,
        lamina_half_width_fraction=100.0 / 350.0,
        catalase_spatial_mode="central_lamina_masked",
        catalase_k1_override_s_inv=18.0 * 0.42,
        d_o2_hydrogel_override_cm2_s=1.5e-5,
        perfusion_spatial_mode="outer_annulus_masked_all_species",
        external_standoff_um=50.0 * 1.5e-5 / 1.584e-5,
        _assess_numerical_gate=False,
    )
    assert math.isclose(_r8_seroma50_pb30["min_slab_po2_unrounded_mmhg"], 17.114450, abs_tol=1e-4)
    assert math.isclose(
        _r8_seroma50_pb30["spatial_control_volume_ledger"]["local_islet_volume_fraction_for_unit_cell"],
        0.10602875,
        abs_tol=1e-6,
    )
    assert math.isclose(
        _r8_seroma50_pb30["wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001"][
            "two_zone_fv_gel_only_supply_120um_core_po2_unrounded_mmhg"
        ],
        1.714933,
        abs_tol=1e-3,
    )
    assert _r8_seroma50_pb30["passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate"] is False
    assert math.isclose(
        _r8_seroma50_pb30["wigner_seitz_100um_dxb_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_100um_core_po2_unrounded_mmhg"
        ],
        5.933834,
        abs_tol=1e-3,
    )
    assert (
        _r8_seroma50_pb30["wigner_seitz_prescribed_cluster_pericellular_shell_audit_dec064"][
            "passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate"
        ]
        is True
    )
    assert (
        _r8_seroma50_pb30["wigner_seitz_prescribed_cluster_pericellular_shell_audit_dec064"][
            "passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only"
        ]
        is False
    )

    for _mk_r8, _sk_r8, _val_r8, _extra_r8 in (
        ("obs_cd47_hour0_molecules", "obs_cd47_hour0_molecules_se", 220000.0, {}),
        ("obs_hlae_hour0_molecules", "obs_hlae_hour0_molecules_se", 18000.0, {}),
        ("obs_serpinb9_kass_M_s", "obs_serpinb9_kass_M_s_se", 5000000.0, {"obs_serpinb9_log10_kass": math.log10(5000000.0)}),
        ("obs_gate0b_hlac_sct_molecules_per_cell", "obs_gate0b_hlac_sct_molecules_se", 30000.0, {}),
    ):
        _omit_r8 = {k: v for k, v in full_ext_input.items() if k != _sk_r8}
        _omit_r8[_mk_r8] = _val_r8
        _omit_r8.update(_extra_r8)
        _dt_r8 = wetlab.evaluate_stage0_decision_tree(_omit_r8)
        assert _dt_r8["has_explicit_release_critical_molecule_and_rate_se"] is False
        assert _dt_r8["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
        _cal_r8 = wetlab.run_calibration(
            {
                **_omit_r8,
                "batch_id": f"BATCH-R8-OMIT-{_mk_r8}",
                "mcb_clone_id": "CLONE-R8-OMIT",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            },
            is_dry_run=False,
        )
        assert (
            _cal_r8["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
            != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
        )

    # 28. AUDIT-REF [DEC-065] & AUDIT-REF [DEC-066] — Layer-4 80-um Seroma Note Provenance & _extract_canonical_measurements Scale-Aware Default SE:
    _r8_seroma80_pb35 = l4_mod._run_single_pouch_pde_regime(
        "R8_SEROMA80_PB35_LOWLOAD_1500",
        0.45 * 0.60,
        "central",
        68.0,
        True,
        0.97,
        boundary_po2_mmhg=35.0,
        demand_multiplier_extra=0.60,
        areal_loading_ieq_cm2=1500.0,
        cell_intrinsic_ocr_multiplier=1.0,
        prescribed_cluster_diameter_um=100.0,
        nz_override=201,
        demand_spatial_mode="outer_annulus",
        lamina_d_o2_scale=0.42 ** 1.5,
        lamina_half_width_fraction=100.0 / 350.0,
        catalase_spatial_mode="central_lamina_masked",
        catalase_k1_override_s_inv=18.0 * 0.42,
        d_o2_hydrogel_override_cm2_s=1.5e-5,
        perfusion_spatial_mode="outer_annulus_masked_all_species",
        external_standoff_um=80.0 * 1.5e-5 / 1.584e-5,
        _assess_numerical_gate=False,
    )
    assert math.isclose(_r8_seroma80_pb35["min_slab_po2_unrounded_mmhg"], 18.470753, abs_tol=1e-4)
    assert (
        _r8_seroma80_pb35["wigner_seitz_prescribed_cluster_pericellular_shell_audit_dec064"][
            "passes_strict_gel_only_two_zone_fv_2p0mmHg_core_gate"
        ]
        is True
    )
    _canon_meas_r8 = wetlab._extract_canonical_measurements(
        {
            "obs_stage7_yield_ratio": 0.90,
            "obs_stage7_yield_se": 0.03,
            "obs_serpinb9_kass_M_s": 5000000.0,
            "obs_gate0b_hlac_sct_molecules_per_cell": 30000.0,
        }
    )
    assert math.isclose(_canon_meas_r8["obs_serpinb9_kass_M_s_se"], 50000.0, rel_tol=1e-9)
    assert math.isclose(_canon_meas_r8["obs_gate0b_hlac_sct_molecules_se"], 1500.0, rel_tol=1e-9)

    # 29. AUDIT-REF [DEC-067] — Exact 80-um Seroma (75.7576 um -> 18.470753 mmHg -> [14.9767, 7.2140] & [16.2198, 11.1907]) Provenance:
    assert math.isclose(
        _r8_seroma80_pb35["wigner_seitz_100um_dxb_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg"
        ],
        14.976720,
        abs_tol=1e-4,
    )
    assert math.isclose(
        _r8_seroma80_pb35["wigner_seitz_100um_dxb_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_100um_core_po2_unrounded_mmhg"
        ],
        7.213984,
        abs_tol=1e-4,
    )
    assert math.isclose(
        _r8_seroma80_pb35["wigner_seitz_80um_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg"
        ],
        16.219850,
        abs_tol=1e-4,
    )
    assert math.isclose(
        _r8_seroma80_pb35["wigner_seitz_80um_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_80um_core_po2_unrounded_mmhg"
        ],
        11.190726,
        abs_tol=1e-4,
    )
    assert "stage0_custom_defaults_hist.get" in (ROOT / "scripts/ingest_wetlab_results.py").read_text(encoding="utf-8")
    # S2-B1: Pin direct 80.0 um external_standoff_um case (c) -> 17.922458 mmHg, [14.434758, 6.693928], [15.674277, 10.652800]:
    _r8_seroma80_direct = l4_mod._run_single_pouch_pde_regime(
        "R8_SEROMA80_DIRECT_PB35_LOWLOAD_1500",
        0.45 * 0.60,
        "central",
        68.0,
        True,
        0.97,
        boundary_po2_mmhg=35.0,
        demand_multiplier_extra=0.60,
        areal_loading_ieq_cm2=1500.0,
        cell_intrinsic_ocr_multiplier=1.0,
        prescribed_cluster_diameter_um=80.0,
        nz_override=201,
        demand_spatial_mode="outer_annulus",
        lamina_d_o2_scale=0.42 ** 1.5,
        lamina_half_width_fraction=100.0 / 350.0,
        catalase_spatial_mode="central_lamina_masked",
        catalase_k1_override_s_inv=18.0 * 0.42,
        d_o2_hydrogel_override_cm2_s=1.5e-5,
        perfusion_spatial_mode="outer_annulus_masked_all_species",
        external_standoff_um=80.0,
        _assess_numerical_gate=True,
    )
    assert math.isclose(_r8_seroma80_direct["min_slab_po2_unrounded_mmhg"], 17.922458, abs_tol=1e-4)
    assert math.isclose(
        _r8_seroma80_direct["wigner_seitz_100um_dxb_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg"
        ],
        14.434758,
        abs_tol=1e-4,
    )
    assert math.isclose(
        _r8_seroma80_direct["wigner_seitz_100um_dxb_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_100um_core_po2_unrounded_mmhg"
        ],
        6.693928,
        abs_tol=1e-4,
    )
    assert math.isclose(
        _r8_seroma80_direct["wigner_seitz_80um_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_surface_po2_unrounded_mmhg"
        ],
        15.674277,
        abs_tol=1e-4,
    )
    assert math.isclose(
        _r8_seroma80_direct["wigner_seitz_80um_cluster_pericellular_shell_audit_dec064"][
            "two_zone_fv_gel_only_supply_80um_core_po2_unrounded_mmhg"
        ],
        10.652800,
        abs_tol=1e-4,
    )
    assert (
        _r8_seroma80_direct["wigner_seitz_prescribed_cluster_pericellular_shell_audit_dec064"][
            "prescribed_cluster_diameter_um"
        ]
        == 80.0
    )
    assert (
        _r8_seroma80_direct["wigner_seitz_prescribed_cluster_pericellular_shell_audit_dec064"][
            "passes_cluster_surface_15mmHg_full_first_phase_gsis_gate_gel_only"
        ]
        is True
    )

    # 30. AUDIT-REF [DEC-071] — F2-01 (Rare-Event Assay Groups Atomic Pooling), F1-02 (7 Capped Immune Endpoints Explicit SE & Conservative Defaults), F1-01 (Sticky S9 Vocabulary & Single-Cell Mosaic Floor), F1-03..F1-07, F2-02..F2-05, L-01, D-01, A2-01, B2-01..02, D2-01..02, E-01..04, S2-B1:
    # - F2-01: Multi-batch rare-event denominator pooling must NOT pair stale Batch-1 depth (e.g. 500,000x) with a new Batch-2 VAF (e.g. 5e-5) when Batch-2 omits depth
    _b1_f201 = wetlab.run_calibration(
        {
            **full_ext_input,
            "batch_id": "BATCH-F201-1",
            "mcb_clone_id": "CLONE-F201",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "obs_tp53_ppm1d_bcor_vaf": 0.0,
            "obs_tp53_ppm1d_bcor_duplex_depth": 500000.0,
        },
        is_dry_run=False,
    )
    assert (
        _b1_f201["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
        == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    )
    _b2_f201 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-F201-2",
            "mcb_clone_id": "CLONE-F201",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "obs_tp53_ppm1d_bcor_vaf": 5e-5,  # at 62,958x (~3.15 mutants) -> Poisson upper95 > 1e-4 -> S9 stop!
        },
        is_dry_run=False,
        prior_posterior_report=_b1_f201,
    )
    assert _b2_f201["stage0_decision_tree_evaluation"]["safety_kill_s9_triggered"] is True
    assert (
        _b2_f201["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
        == "HARD_SAFETY_STOP_REJECTED"
    )
    # Also verify CAST-seq stale depth is dropped when obs_castseq_translocation_freq is re-observed without depth:
    _b1_cast_f201 = wetlab.run_calibration(
        {
            **full_ext_input,
            "batch_id": "BATCH-F201-CAST-1",
            "mcb_clone_id": "CLONE-F201-CAST",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "obs_castseq_translocation_freq": 0.0,
            "obs_castseq_duplex_or_amplicon_depth": 500000.0,
        },
        is_dry_run=False,
    )
    _b2_cast_f201 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-F201-CAST-2",
            "mcb_clone_id": "CLONE-F201-CAST",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "obs_castseq_translocation_freq": 0.0,  # without depth -> HOLD_CASTSEQ_DUPLEX_DEPTH_UNMEASURED
        },
        is_dry_run=False,
        prior_posterior_report=_b1_cast_f201,
    )
    assert (
        _b2_cast_f201["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"]
        != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    )

    # - F1-02: Omitting explicit SE on any of the 7 capped immune lysis/phagocytosis endpoints blocks COMPLETE readiness, and borderline means fail under conservative default SEs:
    for _imm_mk, _imm_sk, _border_val in (
        ("obs_phagocytosis_index_pct", "obs_phagocytosis_index_se", 4.2),
        ("obs_nkg2a_lysis_pct", "obs_nkg2a_lysis_se", 7.0),
        ("obs_nkg2c_adaptive_nk_lysis_pct", "obs_nkg2c_adaptive_nk_lysis_se", 7.0),
        ("obs_kir_nkg2a_neg_nk_lysis_pct", "obs_kir_nkg2a_neg_nk_lysis_se", 7.0),
        ("obs_complement_c3b_opsonized_phagocytosis_pct", "obs_complement_c3b_opsonized_phagocytosis_se", 4.2),
        ("obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct", "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_se", 7.0),
        ("obs_il12_il18_sirpa_low_nk_lysis_pct", "obs_il12_il18_sirpa_low_nk_lysis_se", 8.5),
    ):
        _omit_imm = {k: v for k, v in full_ext_input.items() if k != _imm_sk}
        _dt_omit_imm = wetlab.evaluate_stage0_decision_tree(_omit_imm)
        assert _dt_omit_imm["has_explicit_release_critical_molecule_and_rate_se"] is False
        assert (
            _dt_omit_imm["lot_release_readiness_tier"]
            == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED"
        )
        _border_imm = {**_omit_imm, _imm_mk: _border_val}
        _dt_border_imm = wetlab.evaluate_stage0_decision_tree(_border_imm)
        assert _dt_border_imm["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # - F1-01 & F2-02: Sticky S9 vocabulary validation rejects unknown prior_variant_safety_status strings and latches all 4 S9 stop tokens:
    for _sticky_tok in wc.STICKY_S9_PRIOR_VARIANT_SAFETY_STATUSES:
        _dt_sticky_tok = wetlab.evaluate_stage0_decision_tree(
            {**full_ext_input, "prior_variant_safety_status": _sticky_tok}
        )
        assert _dt_sticky_tok["safety_kill_s9_triggered"] is True
        assert _dt_sticky_tok["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"
    with pytest.raises(ValueError, match="prior_variant_safety_status"):
        wetlab.evaluate_stage0_decision_tree(
            {**full_ext_input, "prior_variant_safety_status": "STOPPED_IN_BATCH_1"}
        )
    # Single-cell 20q11 mosaic fraction floor is 1e-6 (not 1e-3):
    _norm_mosaic, _ = wc.normalize_measurements(
        {
            "obs_20q11_bcl2l1_single_cell_mosaic_fraction": 0.002,
            "obs_20q11_bcl2l1_single_cell_mosaic_se": 0.0002,
        }
    )
    assert math.isclose(_norm_mosaic["obs_20q11_bcl2l1_single_cell_mosaic_se"], 0.0002, rel_tol=1e-9)

    # - F1-04 & F2-03: Zero-event proportion/percentage gates require explicit SE or >= 1,000 cells scored:
    _dt_vp1_zero_no_denom = wetlab.evaluate_stage0_decision_tree(
        {
            k: v
            for k, v in full_ext_input.items()
            if k not in ("obs_gate0v_day28_persistent_vp1_ie1_cells_scored", "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_se")
        }
    )
    assert _dt_vp1_zero_no_denom["zero_event_proportion_missing_denominator_or_explicit_se"] is True
    assert (
        _dt_vp1_zero_no_denom["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_ZERO_EVENT_PROPORTION_GATE_REQUIRES_EXPLICIT_SE_OR_SCORED_CELL_DENOMINATOR"
    )
    assert _dt_vp1_zero_no_denom["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    # With 100 cells scored, 0.0% has Clopper-Pearson upper95 = 2.95% > 0.1% -> Gate 0V fails closed:
    _dt_vp1_100cells = wetlab.evaluate_stage0_decision_tree(
        {**full_ext_input, "obs_gate0v_day28_persistent_vp1_ie1_cells_scored": 100.0}
    )
    assert _dt_vp1_100cells["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # - F1-05 & F1-06 & F1-07 & L-01 & F2-05:
    #   Track-2 B'-SCT variants pass validate_product_and_estimand_compatibility;
    #   residual PSC uses validated SE floor;
    #   obs_gate0e_chop_log2fc S6 threshold is log2(1.5) (~0.585);
    #   obs_pouch_max_h2o2_uM > 5.0 uM fails closed;
    #   stale s8b_epistasis_estimand_scale is not persisted across batches.
    wc.validate_product_and_estimand_compatibility(
        "ISLET1_B_SCT_C1_CIITA", None, ["obs_stage7_yield_ratio", "obs_cd47_fold"]
    )
    wc.validate_product_and_estimand_compatibility(
        "ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500", None, ["obs_stage7_yield_ratio", "obs_cd47_fold"]
    )
    _chop_input = {k: v for k, v in full_ext_input.items() if k != "obs_chop_coincidence_fold"}
    _dt_chop_065 = wetlab.evaluate_stage0_decision_tree(
        {**_chop_input, "obs_gate0e_chop_log2fc": 0.65, "obs_gate0e_chop_log2fc_se": 0.01}
    )
    assert (
        _dt_chop_065["switch_decisions"]["S6_B2M_Orphan_Heavy_Chain_ERAD_Gate_0E"]
        != "RETAIN_STANDARD_STOICHIOMETRY"
    )
    _dt_h2o2_7 = wetlab.evaluate_stage0_decision_tree(
        {**full_ext_input, "obs_pouch_max_h2o2_uM": 7.0, "obs_pouch_max_h2o2_se": 0.1}
    )
    assert _dt_h2o2_7["primary_biological_or_pouch_gate_failed"] is True
    assert _dt_h2o2_7["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # - E-01..E-04: Layer-5 UNRESOLVED cloak token, schema typo rejection, CD47 SIRPa-compromised hold, extreme KIR > 0.70 hold, and rare_snp_pam_disrupted hold:
    with pytest.raises(ValueError, match="Unrecognized patient input keys"):
        l5_mod.compile_personalized_patient_blueprint(
            {
                "patient_id": "PT-TYPO",
                "ancestry": "EUR",
                "nkg2a_neg_kir_educated_nk_fraction": 0.25,
                "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
                "ifn1_ish_genetic_risk_percentile": 35.0,
                "sirpa_genotype": "v1/v1",
                "recipient_anti_hlae_dsa": True,  # typo!
            }
        )
    _bp_pam_disrupted = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-PAM-DISRUPT",
            "ancestry": "AFR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.25,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 35.0,
            "sirpa_genotype": "v1/v1",
            "carries_rare_snp_in_standard_grna": True,
            "rare_snp_pam_disrupted": True,
        }
    )
    assert (
        "HOLD_RARE_SNP_PAM_DISRUPTED_REQUIRES_SHIFTED_PROTOSPACER_OR_ORTHOGONAL_NUCLEASE_QUALIFICATION"
        in _bp_pam_disrupted["all_active_qualification_holds"]
    )
    assert _bp_pam_disrupted["lot_release_cleared"] is False

    # 31. Gemini 3.1 Audit 10-Worker Swarm (Convergence Round 8, Step 9; DEC-072):
    # - AUDIT-REF: rare_snp_pam_disrupted=True on Track 1/2 with donor_mcb_has_protospacer_snp=True
    #   still triggers HOLD_RARE_SNP_PAM_DISRUPTED, omits PATIENT_ALLELE_MATCHED_CUSTOM_SPACER_SYNTHESIS,
    #   and applies the -0.35 survival penalty in evaluate_patient_regimen_survival.
    _pt_track1_pam_mcb = {
        "patient_id": "PT-TRACK1-PAM-MCB",
        "ancestry": "AFR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.25,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 35.0,
        "sirpa_genotype": "v1/v1",
        "donor_mcb_carries_rare_snp_in_standard_grna": True,
        "rare_snp_pam_disrupted": True,
    }
    _bp_t1_pam = l5_mod.compile_personalized_patient_blueprint(_pt_track1_pam_mcb)
    assert (
        "HOLD_RARE_SNP_PAM_DISRUPTED_REQUIRES_SHIFTED_PROTOSPACER_OR_ORTHOGONAL_NUCLEASE_QUALIFICATION"
        in _bp_t1_pam["all_active_qualification_holds"]
    )
    assert "PATIENT_ALLELE_MATCHED_CUSTOM_SPACER_SYNTHESIS" not in _bp_t1_pam["compiled_grna_panel"]
    _surv_t1_pam = l5_mod.evaluate_regimen_survival_probability(_pt_track1_pam_mcb, _bp_t1_pam)
    assert _surv_t1_pam < 0.75
    # - AUDIT-REF: SIRPA v8 allele triggers HOLD_CD47_SIRPA_BINDING_COMPROMISED_REQUIRES_ALTERNATIVE_STEALTH_QUALIFICATION
    _bp_v8 = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-SIRPA-V8",
            "ancestry": "EAS",
            "nkg2a_neg_kir_educated_nk_fraction": 0.25,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 35.0,
            "sirpa_genotype": "v1/v8",
        }
    )
    assert (
        "HOLD_CD47_SIRPA_BINDING_COMPROMISED_REQUIRES_ALTERNATIVE_STEALTH_QUALIFICATION"
        in _bp_v8["all_active_qualification_holds"]
    )
    # - AUDIT-REF & B01: zero_event_proportion_missing_denominator_or_explicit_se exported in return dict
    #   and sub-5e-4 floor SE without *_cells_scored on k=0 does not bypass denominator requirement:
    assert "zero_event_proportion_missing_denominator_or_explicit_se" in _dt_chop_065
    _zero_ev_subfloor_input = {
        k: v
        for k, v in full_ext_input.items()
        if k not in ("obs_gate0v_day28_persistent_vp1_ie1_cells_scored",)
    }
    _dt_zero_subfloor = wetlab.evaluate_stage0_decision_tree(
        {
            **_zero_ev_subfloor_input,
            "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct": 0.0,
            "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_se": 1e-6,
        }
    )
    assert _dt_zero_subfloor["zero_event_proportion_missing_denominator_or_explicit_se"] is True
    assert _dt_zero_subfloor["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    # - AUDIT-REF: Both obs_hc10_surface_fold_wt (1.1, passing) and obs_gate0e_hc10_surface_mfi_ratio (2.1, failing > 1.8)
    #   supplied together must trigger S9/S6 via max() rather than masking the failing alias via elif:
    _dt_both_hc10 = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_input,
            "obs_hc10_surface_fold_wt": 1.1,
            "obs_hc10_surface_fold_se": 0.05,
            "obs_gate0e_hc10_surface_mfi_ratio": 2.1,
            "obs_gate0e_hc10_surface_mfi_se": 0.05,
        }
    )
    assert _dt_both_hc10["safety_kill_s9_triggered"] is True
    assert "ER_CHOP_HC10" in _dt_both_hc10["s9_causal_failed_categories"]

    # 32. DEC-073 (GPT AUDIT-REF Audit Round 9 Step 1 Layers 1-5 Audit) regression checks:
    # - R9-S2-L1-01: X-grid role and blueprint prose accurately identify #2 (ER_IAPP_PRAMLINTIDE_RNLS_WT) vs #3 (AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT) K=4, B=1 lexical tie
    _r9_rank_atlas = (
        json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists()
        else None
    )
    if _r9_rank_atlas is not None:
        _r9_x_role = _r9_rank_atlas["champion_challenger_portfolio"]["slots"]["X"]["role"]
        assert "ER_IAPP_PRAMLINTIDE_RNLS_WT" in _r9_x_role and "AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT" in _r9_x_role
    # - R9-S2-L2-01..03: Arm 1C 3,300 bp constitutive + 3,459 bp promoter-trapped in product registry,
    #   3,915-bp bicistronic control uses authentic 135-bp wpre_mut6_589[:135], and CD47 note maps 112/139/148 links
    _r9_arm1c_prod = pc.get_product_definition("ISLET1_A4_DXS_ARM1C_L4CD47")
    assert "3,300 bp constitutive" in _r9_arm1c_prod["cassette_topology"]["B2M"]
    assert "3,459 bp promoter-trapped" in _r9_arm1c_prod["cassette_topology"]["B2M"]
    _r9_l2_atlas = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    _r9_c_b2m = next(c for c in _r9_l2_atlas["construct_sequence_and_grammar_verification"]["constructs"] if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01")
    _r9_bf_node = _r9_c_b2m["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]["Locus3_BasalFloor_ScarFree_SERPINB9_PDL1_challenger_muse_r3"]
    assert "WPREmut6_5prime_core_135bp" in _r9_bf_node["ssa_direct_repeat_elimination_audit_argon_r4"]["l3_bfc_3915bp_single_promoter_bicistron"]
    assert _r9_bf_node["ssa_direct_repeat_elimination_audit_argon_r4"]["exact_50mer_window_repeat_audit_kimi31_r4"]["l3_bfc_3915bp_bicistronic_duplicated_50mer_windows"] == 0
    _r9_cd47_note = json.loads((ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text(encoding="utf-8"))["gtex_v8_v10_human_tissue_expression"]["CD47"]["coordinate_convention_note"]
    assert "0.756757) yields Arm 1B = 238,004" in _r9_cd47_note and "4.966803x" in _r9_cd47_note
    assert "283,828" in _r9_cd47_note and "299,103" in _r9_cd47_note
    # - R9-S2-L5-01 & R9-S1-L5-01: NO_SCT and UNMATERIALIZED trigger the absent-cloak penalty on both field paths,
    #   and rare_snp_pam_disrupted=True unconditionally triggers HOLD_RARE_SNP_PAM_DISRUPTED even when carrier flags are omitted/False
    _r9_pt_plan_b = {
        "patient_id": "PT-R9-PLANB-CLOAK",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 0.0,
        "ifn1_ish_genetic_risk_percentile": 50.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C1",
        "hlac_4digit_relation": "SELF_4DIGIT_MATCH",
        "track2_subbank_mode": "ENDOGENOUS_HLAC_EXPLORATORY",
    }
    _r9_bp_plan_b = l5_mod.compile_personalized_patient_blueprint(_r9_pt_plan_b)
    _r9_score_absent = l5_mod.evaluate_regimen_survival_probability(
        _r9_pt_plan_b, {**_r9_bp_plan_b, "secondary_safe_harbor_cloak": "ABSENT"}
    )
    for _r9_tok in ("NO_SCT", "UNMATERIALIZED"):
        assert (
            abs(
                l5_mod.evaluate_regimen_survival_probability(
                    _r9_pt_plan_b, {**_r9_bp_plan_b, "secondary_safe_harbor_cloak": _r9_tok}
                )
                - _r9_score_absent
            )
            < 1e-9
        )
    _r9_bp_pam_orphan = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-R9-ORPHAN-PAM",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.20,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 0.5,
            "ifn1_ish_genetic_risk_percentile": 30.0,
            "sirpa_genotype": "v1/v1",
            "rare_snp_pam_disrupted": True,
        }
    )
    assert (
        "HOLD_RARE_SNP_PAM_DISRUPTED_REQUIRES_SHIFTED_PROTOSPACER_OR_ORTHOGONAL_NUCLEASE_QUALIFICATION"
        in _r9_bp_pam_orphan["all_active_qualification_holds"]
    )
    assert _r9_bp_pam_orphan["routing_holds_cleared"] is False

    # 33. DEC-074 (GPT AUDIT-REF Audit Round 9 Step 2 Layers 1-5 Audit) regression checks:
    # - R9-S1-L1B-01: Whole-product 9.335-PSC vs 8.0-PSC vs 93.35-PSC tumor harm bounds explicitly distinguished
    if _r9_rank_atlas is not None:
        _r9_s2_harm = _r9_rank_atlas["whole_product_clonogenic_safety_model"]["unconditional_patient_level_tumor_harm_and_two_stage_psc_purge_model_sonnet55_r5"]
        assert abs(_r9_s2_harm["computed_whole_product_9p335psc_harm_with_1pct_mosaic_20q11_true_pre_retrieval"] - 1.04274e-5) < 1e-9
        assert abs(_r9_s2_harm["computed_whole_product_9p335psc_harm_with_1pct_mosaic_20q11_and_99pct_en_bloc_retrieval"] - 1.04274e-7) < 1e-11
        assert abs(_r9_s2_harm["computed_intermediate_93p35psc_10x_mosaic_harm_pre_retrieval"] - 2.0175e-6) < 1e-11
        assert abs(_r9_s2_harm["computed_intermediate_93p35psc_10x_mosaic_harm_with_99pct_retrieval"] - 2.0175e-8) < 1e-12
    # - R9-S1-L3B-01: ISLET-1-BLUEPRINT.md uses the matched-control incremental SA-FasL fratricide rule (Delta_TUNEL+ <= +5.0 percentage points)
    _r9_s2_bp_md = (ROOT / "ISLET-1-BLUEPRINT.md").read_text(encoding="utf-8")
    assert "Delta\\text{TUNEL}^+ > +5.0\\text{ percentage points}" in _r9_s2_bp_md
    assert "1.0427 \\times 10^{-5}" in _r9_s2_bp_md and "1.0427 \\times 10^{-7}" in _r9_s2_bp_md
    # - R9-S1-L4B-01: Layer 4 exercise margin summary_note and bidirectional_13state dict match the live 9-state (2.7554 / 2.1294 mM at tau=18, 3.1257 mM at tau=12) and 13-state (2.7160 mM at tau=12, 3.1775 mM at baseline tau=7.14) ODE outputs
    _r9_s2_l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    _r9_s2_ex = _r9_s2_l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]["challenged_v0p60_300k_beta_eq_exercise_margin_audit_sol61_r5"]
    assert _r9_s2_ex["bidirectional_13state_and_site_drainage_partition_audit_astra6_r5"]["symmetric_13state_bidirectional_v0p60_alpha1p0_tau12_min_glucose_mM"] == 2.716
    assert _r9_s2_ex["bidirectional_13state_and_site_drainage_partition_audit_astra6_r5"]["symmetric_13state_bidirectional_v0p60_alpha1p0_baseline_tau7p14_min_glucose_mM"] == 3.1775
    assert _r9_s2_ex["bidirectional_13state_and_site_drainage_partition_audit_astra6_r5"]["single_compartment_9state_v0p60_alpha1p0_tau18_min_glucose_mM"] == 2.7554
    assert _r9_s2_ex["bidirectional_13state_and_site_drainage_partition_audit_astra6_r5"]["single_compartment_9state_v0p60_alpha0p50_tau18_min_glucose_mM"] == 2.1294
    assert "2.7554 mM at tau_clear = 18.0 min" in _r9_s2_ex["summary_note"]
    assert "2.1294 mM at tau_clear = 18.0 min" in _r9_s2_ex["summary_note"]
    # - Layer 2 & Layer 5 <5.0/10 notes: SERPINB9 pulse-4 hypoxia-equilibrated unneutralized GzmB = 0.038049 uM, and TYK2 PE Track-2 status label reflects K4_B1_3PASS
    _r9_s2_sb9 = _r9_bf_node["hour0_kinetic_race_and_stoichiometric_bolus_model_muse_r3"]["serpinb9_suicide_substrate_serial_pulse_ode_and_perk_eif2a_isr_audit_sol61_r4"]
    assert _r9_s2_sb9["sustained_12_pulse_22h_depletion_and_stationary_reserve_audit_astra6_r4"]["first_pulse_with_unneutralized_gzmb_bolus_from_hypoxia_equilibrated_start_fISR_times_S0_sonnet55_r4"]["licensed_S_start_0p4005_fISR_0p45_SI_1p0"]["unneutralized_gzmb_uM"] == 0.038049
    assert _r9_s2_sb9["twelve_pulse_hypoxic_coupon_reconciliation_and_synaptic_suppression_gate_sol61_r5"]["unattenuated_0p15uM_gzmb_per_pulse_depletion_envelope"]["basal_0p30uM_under_f_isr_0p45_first_exhaustion_pulse"] == 2
    _r9_s2_bp_tyk2 = l5_mod.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-R9-TRACK2-TYK2",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.55,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 95.0,
            "sirpa_genotype": "v1/v1",
            "recipient_kir_ligand_group": "C1/C1",
        }
    )
    assert (
        _r9_s2_bp_tyk2["b_prime_4digit_hlac_and_b_prime_sct_routing_audit_muse_r3"]["b_prime_sct_status"]
        == "CANONICAL_TRACK2_KIR_ENGAGING_SUBBANK_K4_B1_3PASS_TYK2_PE_SATELLITE"
    )

    # 34. DEC-075 (Audit Convergence Round 9 Step 3 Layers 1-5 Audit) regression checks:
    # - AUDIT-REF & AUDIT-REF: 2A-free CD55(SCR1-4)-GPI + CD59-GPI complement challenger (no upstream GPI-P2A defect) and 10 experimental PDB structures (8 mmCIF + 2 legacy PDB)
    _r9_s3_l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    assert _r9_s3_l3["summary"]["total_experimental_pdb_structures_verified"] == 10
    assert _r9_s3_l3["summary"]["experimental_rcsb_mmcif_structures_verified"] == 8
    assert _r9_s3_l3["summary"]["experimental_rcsb_legacy_pdb_structures_verified"] == 2
    _r9_s3_chandler = _r9_s3_l3["hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4"]["endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4"]["non_hla_chimeric_junction_neoepitope_and_flow_crossmatch_audit_sonnet55_r5"]["membrane_complement_regulator_cd55_scr14_cd59_and_assembled_cassette_ibmir_gate_muse_r6"]["assembled_cassette_chandler_loop_ibmir_and_complement_gate"]["cell_intrinsic_cd55_scr14_cd59_challenger"]
    assert "SV40pA-CD55(SCR1-4)-GPI(antisense) <- EEF1A1-core | hPGK1-min(sense) -> CD59-GPI-bGHpA" in _r9_s3_chandler
    assert "explicitly rejecting 2A-linked CD55-P2A-CD59" in _r9_s3_chandler
    # - AUDIT-REF / R9-S3-L1-01..02 & AUDIT-REF: exact 10x mosaic q_esc = 1.09e-8 (P_pre = 2.0175e-6, P_post = 2.0175e-8), whole-product 9.335e8 denominator for 93.35 PSCs, and CLAUDE.md schema_version 0.7.0
    assert "2.0175 \\times 10^{-6}" in _r9_s2_bp_md and "2.0175 \\times 10^{-8}" in _r9_s2_bp_md
    assert "1.09 \\times 10^{-8}" in _r9_s2_bp_md
    if _r9_rank_atlas is not None:
        assert "2.0175e-6" in _r9_s2_harm["pre_vs_post_retrieval_metastasis_term_reconciliation_muse_r6"]
    if (ROOT / "CLAUDE.md").exists():
        assert "schema_version 0.7.0" in (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    # - AUDIT-REF / s2_03: Kazmi 2026 direct islet surface biotinylation (not microgels) and OR ZB4 rescue > 5.0 pp rule
    assert "Kazmi et al. 2026 JCI Insight (PMID:41842939)" in _r9_s2_bp_md
    assert "STAGE0_DECISION_TREE_ENDPOINT_SPECS" in _r9_s2_bp_md
    assert "ZB4` rescue $\\le 5.0\\text{ pp}$" in _r9_s2_bp_md
    # - R9-S1-L4B-01-RESIDUAL: 0.25-min RK4 adult reference 3.5710 mM (1-min sample 3.5715 mM; literal tau=7.14 3.5713 mM) synchronized
    assert _r9_s2_ex["challenged_v0p60_300k_beta_eq_adult_ref_min_glucose_mM"] == 3.571
    assert _r9_s2_ex["challenged_v0p60_300k_beta_eq_alpha_gain_0p50_min_glucose_mM"] == 3.1277
    assert _r9_s2_ex["challenged_v0p60_300k_beta_eq_combined_immaturity_min_glucose_mM"] == 2.9523

    # 35. DEC-076 (Audit Convergence Round 9 Step 4 Layers 1-5 Audit) regression checks:
    # - AUDIT-REF: experimental_pdb_accession_inventory matches the exact 10 authenticated experimental inputs (4CMM, 4ZQK, 1IAU instead of 2VSC, 1G1Q, 1I8L)
    assert _r9_s3_l3["summary"]["experimental_pdb_accession_inventory"] == [
        "1MHE", "3CDG", "3BZE", "2JJS", "4CMM", "4KJY", "4ZQK", "1IAU", "1EFX", "1IM9"
    ]
    # - AUDIT-REF: Decision Experiment #2 separates Headen/Lei SA-FasL rapamycin from Turan 2024 syngeneic SA-THBD+SA-CD47 without rapamycin
    if _r9_rank_atlas is not None:
        _r9_s4_exp2_why = _r9_rank_atlas["step7_high_value_of_information_gaps"][1]["why_needed"]
        assert "Turan et al. 2024 PMID:38709159 evaluated direct islet surface biotinylation co-displaying SA-THBD [SA-TM] and SA-CD47 in a murine syngeneic intraportal IBMIR/engraftment model without rapamycin" in _r9_s4_exp2_why
    # - AUDIT-REF & L1B-N01: ISLET-1-BLUEPRINT.md binds 0.7249 uM (0.725 uM) to 189,000 IU uniform bulk-slab, 4.0607 uM to localized 11,271.2 IU pore-water dose, and 9.8158 uM to dispersed 11,271.2 IU, plus 2.00 x 10^6 vs 1.92 x 10^6 seeding convention
    assert "4.0607 µM" in _r9_s2_bp_md and "9.8158 µM" in _r9_s2_bp_md and "189,000 IU" in _r9_s2_bp_md
    assert "2.00 × 10^6` total input endocrine cells" in _r9_s2_bp_md
    # - AUDIT-REF & L4B-N1: insulin_content_per_beta_eeq_mU_range == [0.138, 0.461] and 1-min subsample minima == 3.5715 / 3.1277 / 2.9523 mM
    _r9_s4_pool = _r9_s2_l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]["icasp9_ap1903_ablation_insulin_dump_and_staged_cassette_excision_ode_audit_sonnet55_r5"]["intracellular_granule_insulin_pool_300k_beta_eeq_to_500k_total_ieq"]
    assert _r9_s4_pool["insulin_content_per_total_ieq_mU_range"] == [0.083, 0.277]
    assert _r9_s4_pool["insulin_content_per_beta_eeq_mU_range"] == [0.138, 0.461]
    assert _r9_s2_ex["challenged_v0p60_300k_beta_eq_adult_ref_1min_sample_min_glucose_mM"] == 3.5715
    assert _r9_s2_ex["challenged_v0p60_300k_beta_eq_alpha_gain_0p50_1min_sample_min_glucose_mM"] == 3.1277
    assert _r9_s2_ex["challenged_v0p60_300k_beta_eq_combined_immaturity_1min_sample_min_glucose_mM"] == 2.9523
    # - AUDIT-REF: Layer 5 baseline_97p50_conditioning_scope dynamically formats 96.00% K=4 [96.0046%]
    _r9_s4_l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    _r9_s4_scope = _r9_s4_l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["b_prime_4digit_hlac_alloreactivity_and_b_prime_sct_pivot_audit_muse_r3"]["hazard_active_cohort_sensitivity_audit_sonnet55_r4"]["baseline_97p50_conditioning_scope"]
    assert "96.00% K=4 [96.0046%] complexity-discounted" in _r9_s4_scope
    assert "96.65%" not in _r9_s4_scope

    # 36. DEC-077 (Audit Convergence Round 9 Step 5 Layers 1-5 Audit) regression checks:
    # - L3B-01/02/03: Headen et al. 2018 Nat Mater PMID:29867165 replaces bladder cancer ID 2935_5840 and splits SA-FasL vs SA-THBD citations
    assert ("2935" + "5840") not in _r9_s2_bp_md
    assert "Headen et al. 2018 (PMID:29867165)" in _r9_s2_bp_md or "PMID:29867165" in _r9_s2_bp_md
    if _r9_rank_atlas is not None:
        assert ("2935" + "5840") not in _r9_rank_atlas["ranked_candidate_architectures"][0]["locus1_and_locus4_engineering_audit_notes"]["sa_fasl_sa_thbd_immunosuppression_and_fratricide_note"]
        assert "PMID:29867165" in _r9_rank_atlas["ranked_candidate_architectures"][0]["locus1_and_locus4_engineering_audit_notes"]["sa_fasl_sa_thbd_immunosuppression_and_fratricide_note"]
        assert "PMID:29867165" in _r9_s4_exp2_why and ("2935" + "5840") not in _r9_s4_exp2_why
    # - L1A-01 & GF-01: ISLET-1-BLUEPRINT.md documents 19.99-23.07 pp t_13 MDE for 38-run 4-bridge high-VIF 2FI, sigma_run = 5.43 pp (5.4314 pp) vs 4.90 pp, and 0.1041 (0.104090 unrounded) for Scenario A
    assert "19.99–23.07 pp" in _r9_s2_bp_md and "5.43\\text{ pp}" in _r9_s2_bp_md
    assert "**`0.1041`** [`0.104090` unrounded" in _r9_s2_bp_md
    # - L1B-F1, L1B-F2, R9-S1-L3A-01: Unified Slot -> Grid ID -> Stage-0 Product Variant ID roster, IEEE-754 binary64 cohort-rounded average + reverse=True tie-break, and WT alpha3 Asp227/Thr228 conditional default vs D227K/T228A satellite coupon
    if _r9_rank_atlas is not None:
        _r9_s5_slots = _r9_rank_atlas["champion_challenger_portfolio"]["slots"]
        assert "ISLET1_A_STAR_4LOCUS" in _r9_s5_slots["A"]["stage0_product_variant_ids"]
        assert "ISLET1_B_SCT_C1_CIITA" in _r9_s5_slots["B"]["stage0_product_variant_ids"]
        assert "ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500" in _r9_s5_slots["X"]["stage0_molecular_product_variant_ids"]
        assert "reverse=True" in _r9_s5_slots["X"]["role"]
        assert "WT alpha3 Asp227/Thr228 conditional default" in _r9_rank_atlas["champion_challenger_portfolio"]["blueprint_x_stage0_molecular_challenger_spec"]["locus1_arm_1x_challenger"]
    assert "R9-S1-L3A-01 [DEC-077]" in _r9_s3_l3["complexes"][0]["hlae_alpha3_d227_t228_distance_audit"]["structural_interpretation_note"]
    # - L2A-01 & L2A-02: Layer 2 materializes both Arm 1C 3300 bp donors (S3+uORF+5'TOP b1e17584... and constitutive native-parts 259d1654...) and Arm 1A-S3 4-nt PAM-shielded 4880 bp donor
    _r9_s5_l2 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    _r9_s5_q1 = _r9_s5_l2["construct_sequence_and_grammar_verification"]["constructs"][1]["arm1b_vs_arm1x_codon_coordinate_table"]["single_variable_challenger_arms_kimi31_q1"]
    _r9_s5_dxs = _r9_s5_q1["DX_S_Split_Cassette_Arm1C_HLAE_plus_Locus4_CD47_challenger_sonnet55_r1"]
    _r9_s5_sh = _r9_s5_q1["B2M_AsCas12a_S3_4nt_synonymous_PAM_seed_shield_sonnet55_r1"]["unconditional_3locus_recut_immunity_pam_shield_bundle_sonnet55_r2"]["materialized_shielded_homology_arms_and_donors_kimi31_r3"]
    assert _r9_s5_dxs["locus1_arm1c_constitutive_eef1a1_intron_donor_sha256"].startswith("b1e17584")
    assert _r9_s5_dxs["locus1_arm1c_constitutive_native_parts_3300bp_donor_sha256"].startswith("259d1654")
    assert len(_r9_s5_dxs["locus1_arm1c_constitutive_native_parts_3300bp_donor_sequence"]) == 3300
    assert len(_r9_s5_sh["arm1a_4880bp_s3_4nt_pam_shielded_donor_sequence"]) == 4880
    assert "canonical_fasta_vs_shielded_synthesis_reconciliation_note_muse_r9" in _r9_s5_sh
    # - L4A-F01 & L4A-F02: Layer 4 reconciles 0.42 cm3 local sentinel patch OH- generation flux (8.2561e-9 mol/(cm3*s), 2.91x clearance ratio) and two-sided vs per-face 1500 IEQ/cm2 OCR fluxes (9.0124e-11 vs 4.5062e-11 mol/(cm2*s))
    _r9_s5_muse_r6 = _r9_s2_l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["bed_level_convective_o2_budget_and_periimplant_fluid_film_audit_sonnet55_r5"]["periimplant_stagnant_fluid_film_transport_limit"]["coupled_bvp_and_series_fluid_film_core_po2_table_r6"]["cao2_hydroxide_peroxide_inventory_and_peak_ocr_seroma_fick_audit_muse_r6"]
    _r9_s5_cao2 = _r9_s5_muse_r6["cao2_total_alkaline_and_peroxide_stoichiometric_inventory_7cm3_pouch"]
    assert _r9_s5_cao2["sentinel_patch_volume_cm3_at_12cm2_350um"] == 0.42
    assert _r9_s5_cao2["sentinel_local_bicarbonate_washout_clearance_ratio_at_100p78mg_in_0p42cm3"] == 2.91
    _r9_s5_fick = _r9_s5_muse_r6["peak_stimulated_vs_reference_ocr_seroma_fick_reconciliation"]
    assert abs(_r9_s5_fick["pde_two_sided_total_pure_islet_vmax_flux_1500_ieq_cm2_mol_cm2_s"] - 9.0124e-11) < 1e-14
    assert abs(_r9_s5_fick["pde_per_face_pure_islet_vmax_flux_1500_ieq_cm2_mol_cm2_s"] - 4.5062e-11) < 1e-14
    # - L5A-01, L5A-02, M1, M3, M4: Layer 5 evaluate_regimen_survival_probability enforces Bw4/KIR3DL1 -0.18 penalty & SIRPA v8 CD47 loss, fails closed when Plan B omits product_definition, and normalizes uppercase SIRPA genotype
    _r9_s5_pt = {
        "patient_id": "PT-R9-UPPER-V8",
        "ancestry": "AFR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.45,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 50.0,
        "sirpa_genotype": "V1/V8",
        "recipient_kir_ligand_group": "C1/C1",
        "recipient_bw4_kir3dl1_educated": True,
        "bw4_kir3dl1_sorted_nk_lysis_cleared": False,
    }
    _r9_s5_upper_v8 = l5_mod.compile_personalized_patient_blueprint(_r9_s5_pt)
    assert (
        "HOLD_CD47_SIRPA_BINDING_COMPROMISED_REQUIRES_ALTERNATIVE_STEALTH_QUALIFICATION"
        in _r9_s5_upper_v8["all_active_qualification_holds"]
    )
    _r9_s5_reg_no_pd = dict(_r9_s5_upper_v8)
    _r9_s5_reg_no_pd.pop("product_definition", None)
    _r9_s5_bw4_surv = l5_mod.evaluate_regimen_survival_probability(_r9_s5_pt, _r9_s5_reg_no_pd)
    _r9_s5_bw4_cleared_surv = l5_mod.evaluate_regimen_survival_probability(
        {**_r9_s5_pt, "bw4_kir3dl1_sorted_nk_lysis_cleared": True}, _r9_s5_reg_no_pd
    )
    assert _r9_s5_bw4_surv < 0.50
    assert abs((_r9_s5_bw4_cleared_surv - _r9_s5_bw4_surv) - 0.18) < 1e-9

    # 37. DEC-078 (Audit Convergence Round 9 Step 6 Layers 1-5 Audit + Public Bundle Privacy Guard) regression checks:
    # - Privacy Guard: verify_public_bundle_privacy_guard passes on live repo and fails closed on synthetic forbidden leaks
    from scripts import generate_release_manifest as grm_mod
    import tempfile

    _priv_res = grm_mod.verify_public_bundle_privacy_guard(ROOT, grm_mod.PUBLIC_PACKAGE_INCLUDED_PATHS)
    assert _priv_res["passed"] is True
    assert _priv_res["violation_count"] == 0
    with tempfile.TemporaryDirectory() as _tmp_bundle:
        _tb = Path(_tmp_bundle)
        (_tb / "README.md").write_text("Clean file with Kev only.\n", encoding="utf-8")
        assert grm_mod.verify_public_bundle_privacy_guard(_tb, ["README.md"])["passed"] is True
        for _leak_snippet in (
            "/" + "Users/alice/secret",
            "contact: test" + "@gmail.com",
            "token: ghp_" + "A" * 36,
        ):
            (_tb / "README.md").write_text(f"Leak: {_leak_snippet}\n", encoding="utf-8")
            try:
                grm_mod.verify_public_bundle_privacy_guard(_tb, ["README.md"])
                raise AssertionError(f"Expected RuntimeError for leak snippet {_leak_snippet!r}")
            except RuntimeError as _exc:
                assert "PRIVACY GUARD FAILED" in str(_exc)
        # Verify planted CLAUDE.md is blocked as a forbidden public path (A):
        (_tb / "README.md").write_text("Clean file with Kev only.\n", encoding="utf-8")
        (_tb / "CLAUDE.md").write_text("Internal instructions.\n", encoding="utf-8")
        try:
            grm_mod.verify_public_bundle_privacy_guard(_tb, ["README.md", "CLAUDE.md"])
            raise AssertionError("Expected RuntimeError for planted CLAUDE.md")
        except RuntimeError as _exc_cl:
            assert "PRIVACY GUARD FAILED" in str(_exc_cl) and "forbidden path" in str(_exc_cl)
        (_tb / "CLAUDE.md").unlink()
        # Verify fail-closed when require_private_pattern_file=True and private pattern file is missing (B1 / R18):
        try:
            grm_mod.verify_public_bundle_privacy_guard(
                _tb, ["README.md"], source_root=_tb, require_private_pattern_file=True
            )
            raise AssertionError("Expected RuntimeError when private_export_guard_patterns.json is missing")
        except RuntimeError as _exc_missing_priv:
            assert "Missing required private policy file" in str(_exc_missing_priv)
        # Verify machine patterns from private_export_guard_patterns.json are enforced when present (B1):
        _priv_pat_file = ROOT / "data/raw/private_export_guard_patterns.json"
        if _priv_pat_file.exists():
            _priv_cfg = json.loads(_priv_pat_file.read_text(encoding="utf-8"))
            _mach_pat = _priv_cfg.get("regex_patterns", {}).get("machines_and_network", "")
            assert bool(_mach_pat)
        _syn_priv_dir = _tb / "data/raw"
        _syn_priv_dir.mkdir(parents=True, exist_ok=True)
        (_syn_priv_dir / "private_export_guard_patterns.json").write_text(
            '{"regex_patterns":{"machines_and_network":"SYNTHETIC_PRIVATE_MACHINE_TOKEN"}}',
            encoding="utf-8",
        )
        (_tb / "README.md").write_text("Host: SYNTHETIC_PRIVATE_MACHINE_TOKEN\n", encoding="utf-8")
        try:
            grm_mod.verify_public_bundle_privacy_guard(
                _tb, ["README.md"], source_root=_tb, require_private_pattern_file=True
            )
            raise AssertionError("Expected RuntimeError for private machines_and_network pattern")
        except RuntimeError as _exc_mach:
            assert "PRIVACY GUARD FAILED" in str(_exc_mach) and "machines_and_network" in str(_exc_mach)
        (_syn_priv_dir / "private_export_guard_patterns.json").unlink()

    # - G-CITES-01 & L3A-NUM-02: doc_citation_metadata_cache.json has all 6 added classic PMIDs and off-topic 15972683/29403493 removed from Layer 3; KPNA2 is P52292 residues 204-212
    _cite_cache = json.loads((ROOT / "data/raw/doc_citation_metadata_cache.json").read_text(encoding="utf-8"))["pmids"]
    for _pmid_chk in ("2124002", "1986378", "1512546", "8120026", "27807284", "26239055"):
        assert _pmid_chk in _cite_cache and bool(_cite_cache[_pmid_chk].get("first_author"))
    assert "Vitiello" in _cite_cache["2124002"]["first_author"]
    assert "Dendrou" in _cite_cache["27807284"]["first_author"]
    assert "Marroqui" in _cite_cache["26239055"]["first_author"]
    _l3_raw_txt = (ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8")
    assert ("1597" + "2683") not in _l3_raw_txt and ("2940" + "3493") not in _l3_raw_txt
    _prod_defs = json.loads((ROOT / "data/protocols/islet1_product_definitions_v0.8.0.json").read_text(encoding="utf-8"))
    assert "UniProt P52292 residues 204-212" in json.dumps(_prod_defs)

    # - L2A-PAR-02 / S2-B1 & L2B-COD-01: unshielded 4048-bp DX-S donor (79433ac7...) materialized alongside PAM-shielded (fbbd5b8f...) and 5'TOP (a6eb1de8...); Layer-2 free heavy-chain citations verified
    assert _r9_s5_dxs["locus4_icasp9_p2a_cd47_native_parts_unshielded_donor_sha256"].startswith("79433ac7")
    assert len(_r9_s5_dxs["locus4_icasp9_p2a_cd47_native_parts_unshielded_donor_sequence"]) == 4048
    assert _r9_s5_sh["dxs_locus4_4048bp_unshielded_5top_cagaga_donor_sha256"].startswith("a6eb1de8")
    _l2_fhc_cites = _r9_s5_l2["stage0_to_stage7_perturbation_plane"]["er_heavy_chain_flux_and_upr_hazard_muse_r3"]["b2m_ko_free_heavy_chain_surface_leakage_and_exogenous_b2m_rescue_gate_sol61_r5"]["primary_citations"]
    assert _l2_fhc_cites == ["PMID:2124002", "PMID:1986378", "PMID:1512546", "PMID:8120026"]

    # - L3B-COD-01 & L3B-COD-02: SA-FasL Decision Experiment #2 endpoints wired in STAGE0_DECISION_TREE_ENDPOINT_SPECS and efferocytosis/iEC vascular holds decoupled from CD47/phagocytosis
    for _safasl_k in ("obs_safasl_delta_tunel_pp", "obs_safasl_zb4_anti_fas_rescue_pp", "obs_safasl_treg_viable_retention_pct_of_thbd"):
        assert _safasl_k in s0_mean_keys
    _dt_safasl_fail = wetlab.evaluate_stage0_decision_tree({**full_ext_input, "obs_safasl_delta_tunel_pp": 7.5, "obs_safasl_delta_tunel_se": 0.2})
    assert _dt_safasl_fail["switch_decisions"]["S5_Microgel_Boundary_Cuff_Gate_0G"] == "DROP_SA_FASL_RETAIN_SA_THBD_PLUS_SA_FHL1"
    _dt_effero_fail = wetlab.evaluate_stage0_decision_tree({**full_ext_input, "obs_effero_clearance_pct": 65.0, "obs_effero_clearance_se": 1.0})
    assert _dt_effero_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_APOPTOTIC_EFFEROCYTOSIS_CLEARANCE_IMPAIRED_DAMPS_RISK"
    _dt_iec_fail = wetlab.evaluate_stage0_decision_tree({**full_ext_input, "obs_iec_anastomosis_t50_h": 42.0, "obs_iec_anastomosis_t50_se": 1.0})
    assert _dt_iec_fail["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_PAN_LINEAGE_CD47_QPCTL_OR_VASCULAR_SHIELD_FAILED"

    # - L4A-F2, L4B-F02, L5B-NUM-01..02, L5A-G3-01: Layer-4 core_gate_basis_note and coupled 13-state immaturity/lag audit; Layer-5 live challenger junction metrics, 37-nt tevopreQ1 + 5-nt AUUUU segmentation, and Track-independent moderate KIR hold on Arm 1B
    _r9_s6_l4 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    _r9_s6_sens = _r9_s6_l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["lead_regime_4_scenario_sensitivity_table_sol_v80_12"]
    assert all("core_gate_basis_note" in s for s in _r9_s6_sens)
    _r9_s6_env = _r9_s6_l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]["challenged_v0p60_300k_beta_eq_exercise_margin_audit_sol61_r5"]["portal_drainage_lag_by_clearance_exercise_envelope_astra6_r6"]
    assert "coupled_13state_portal_vs_subcutaneous_immaturity_and_lag_audit_sonnet55_r9" in _r9_s6_env
    _r9_s6_l5 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    _r9_s6_epeg = _r9_s6_l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]["guides"][4]["epegRNA_extension_specification"]
    assert _r9_s6_epeg["separate_intervening_linker_length_nt"] == 0
    assert len(_r9_s6_epeg["tevopreq1_pseudoknot_37nt_core_rna_5to3"]) == 37
    assert len(_r9_s6_epeg["u6_termination_tail_5nt_rna_5to3"]) == 5
    _r9_s6_junc = _r9_s6_l5["part_a_crispr_grna_compiler_and_polymorphism_audit"]["post_hdr_junction_cfd_and_split_protospacer_audit_claude_v80_004"]["gate_zero_donor_wgs_and_precloning_junction_recleavage_gate_muse_r1"]["stage0_challenger_junction_equivalence_audit_argon_r2"]
    assert "0.173249" in _r9_s6_junc["arm1c_constitutive_and_arm1a_lite_locus1"]
    assert "0.157781" in _r9_s6_junc["arm1b_cd47_c15g_and_l1l2_synonymous_locus1"]
    assert "0.114255" in _r9_s6_junc["arm1c_promoter_trapped_and_arm1b_sense_terminator_shield_locus1"]
    assert "0.068317" in _r9_s6_junc["arm1c_promoter_trapped_and_arm1b_sense_terminator_shield_shield_locus1"] if "arm1c_promoter_trapped_and_arm1b_sense_terminator_shield_shield_locus1" in _r9_s6_junc else "0.068317" in _r9_s6_junc["arm1c_promoter_trapped_and_arm1b_sense_terminator_shield_locus1"]

    # 38. DEC-080 (GPT AUDIT-REF Audit Convergence Round 9 Step 8 Final Closure):
    # - R9S8-L3B-01: obs_safasl_treg_viable_retention_pct_of_thbd 95% lower-CI < 85.0% triggers S5_Microgel_Boundary_Cuff_Gate_0G = DROP_SA_FASL_RETAIN_SA_THBD_PLUS_SA_FHL1
    _dt_safasl_treg80 = wetlab.evaluate_stage0_decision_tree(
        {**full_ext_input, "obs_safasl_treg_viable_retention_pct_of_thbd": 80.0, "obs_safasl_treg_viable_retention_se": 1.0}
    )
    assert _dt_safasl_treg80["switch_decisions"]["S5_Microgel_Boundary_Cuff_Gate_0G"] == "DROP_SA_FASL_RETAIN_SA_THBD_PLUS_SA_FHL1"
    assert _dt_safasl_treg80["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S5_FASL_BYSTANDER_TOXICITY_REQUIRES_THBD_FHL1_ONLY_CUFF_QUALIFICATION"

    # - R9S8-L5B-01: OpenAI project key (sk-proj-...) blocked by verify_public_bundle_privacy_guard
    with tempfile.TemporaryDirectory() as _tb_sk:
        _tb_sk_p = Path(_tb_sk)
        (_tb_sk_p / "README.md").write_text("Leak: " + "sk-proj-" + "A" * 156 + "\n", encoding="utf-8")
        try:
            grm_mod.verify_public_bundle_privacy_guard(_tb_sk_p, ["README.md"])
            raise AssertionError("Expected RuntimeError for synthetic sk-proj- token")
        except RuntimeError as _exc_sk:
            assert "PRIVACY GUARD FAILED" in str(_exc_sk) and "private_keys_and_tokens" in str(_exc_sk)

    # - R9S8-L4B-01 & L4-A: coupled_13state_portal_vs_subcutaneous_immaturity_and_lag_audit_sonnet55_r9, staged_sentinel_functional_dosing_and_paediatric_exercise_derating_muse_r6, and hour0_newton_bvp_spatial_and_dose_controls_uM dynamically computed from live solvers
    _r9_s8_13st = _r9_s6_env["coupled_13state_portal_vs_subcutaneous_immaturity_and_lag_audit_sonnet55_r9"]
    assert abs(_r9_s8_13st["standard_protocol_12state_v0p72_360k_beta_eq_lag10_alpha0p50_min_glucose_mM"] - 3.5059) < 0.001
    assert abs(_r9_s8_13st["standard_protocol_portal_13state_v0p72_360k_beta_eq_lag10_alpha0p50_canonical_kI0p14_min_glucose_mM"] - 3.5059) < 0.001
    assert abs(_r9_s8_13st["standard_protocol_portal_13state_v0p72_360k_beta_eq_lag10_alpha0p50_tau8p0_min_glucose_mM"] - 3.3399) < 0.001
    assert abs(_r9_s8_13st["low_glucose_challenge_13state_portal_v0p72_360k_beta_eq_lag10_alpha0p50_min_glucose_mM"] - 1.8651) < 0.001
    _r9_s8_sent = _r9_s6_l4["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]["icasp9_ap1903_ablation_insulin_dump_and_staged_cassette_excision_ode_audit_sonnet55_r5"]["staged_sentinel_functional_dosing_and_paediatric_exercise_derating_muse_r6"]
    assert abs(_r9_s8_sent["simulated_v0p18_portal_tau8_default_g0_min_glucose_mM"] - 4.6372) < 0.001
    assert abs(_r9_s8_sent["simulated_v0p36_portal_tau8_default_g0_min_glucose_mM"] - 3.3217) < 0.001
    assert abs(_r9_s8_sent["simulated_v0p60_portal_tau8_min_glucose_mM"] - 3.0722) < 0.001
    _r9_s8_cat = _r9_s6_l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["catalase_multi_denominator_iu_ledger_v80_093"]["spatial_enzyme_localization_and_saturation_equivalence_audit"]["hour0_newton_bvp_spatial_and_dose_controls_uM"]["nz201_fine_grid_controls_uM"]
    assert abs(_r9_s8_cat["design_x_100um_source_119um_mask_k1_6p3529_s_inv_22680_iu_diluted_hour0_nz201"] - 2.3806) < 0.001
    assert abs(_r9_s8_cat["layout_b_119um_source_119um_mask_k1_6p3529_s_inv_22680_iu_hour0_nz201"] - 2.0476) < 0.001
    assert abs(_r9_s8_cat["design_x_100um_source_100um_mask_k1_6p3529_s_inv_19059_iu_hour0_nz201"] - 2.4288) < 0.001


























