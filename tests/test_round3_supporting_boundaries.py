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
        assert bx["off_grid_molecular_design_x_locus1_and_lamina_counterfactual_projection_muse_r1"]["off_grid_design_x_constitutive_hlae_composite_utility_projection"] == 0.1160
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
        assert proj["unconditional_composite_utility_if_graveyard_009_retained"] == 0.1084
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
    assert len(dt_res["switch_decisions"]) == 12
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

    # - R7S2-F02 & Q1/Q2: Historical ledger strict boolean validation (118 persisted boolean flags, 122 ledger metadata booleans) & Layer-4 5.71 mmHg prose alignment
    assert len(wetlab.STAGE0_PERSISTED_BOOLEAN_FLAGS) == 118
    assert "require_full_5endpoint_s9_panel" in wetlab.STAGE0_PERSISTED_BOOLEAN_FLAGS
    assert "require_full_7metric_combined_stress_panel" in wetlab.STAGE0_PERSISTED_BOOLEAN_FLAGS
    assert len(wetlab.HISTORICAL_LEDGER_BOOLEAN_METADATA_KEYS) == 122
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

    assert len(wc_mod.STAGE0_DECISION_TREE_BOOLEAN_KEYS) == 122
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

    # 34. DEC-075 & DEC-090 (IIP-0011) Layer 3 structural reference inventory checks:
    # - AUDIT-REF & IIP-0011: 2A-free CD55(SCR1-4)-GPI + CD59-GPI complement challenger (no upstream GPI-P2A defect) and 13 experimental PDB structures (11 mmCIF + 2 legacy PDB, adding 1BL4, 1JXQ, 1NW9 for Locus 4 iCasp9)
    _r9_s3_l3 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    assert _r9_s3_l3["summary"]["total_experimental_pdb_structures_verified"] == 15
    assert _r9_s3_l3["summary"]["experimental_rcsb_mmcif_structures_verified"] == 11
    assert _r9_s3_l3["summary"]["experimental_rcsb_legacy_pdb_structures_verified"] == 4
    assert _r9_s3_l3["summary"]["canonical_experimental_pdb_structures_verified"] == 13
    assert _r9_s3_l3["summary"]["challenger_iip0031_experimental_pdb_structures_verified"] == 2
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
        "1MHE", "3CDG", "3BZE", "2JJS", "4CMM", "4KJY", "4ZQK", "1IAU", "1EFX", "1IM9", "1BL4", "1JXQ", "1NW9", "2WII", "5O32"
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

    # 39. DEC-081 / IIP-0003 (Re-judged four-hurdle conditional scenario specification, score/rank measures & OQ-016; reassessed 2026-11-15 after IIP-0016..IIP-0019):
    from scripts import rank_islet1_stacks as ranker
    _iip3 = ranker._build_odds_composition_and_profile_scope_audit_iip_0003()
    assert _iip3["proposal_id"] == "IIP-0003" and _iip3["open_question_id"] == "OQ-016"
    assert _iip3["initial_reassessment_date"] == "2026-11-05" and _iip3["reassessment_date"] == "2026-11-16"
    assert "Linus Schumacher (ORCID 0000-0003-0797-3406, GitHub @ljschumacher, University of Edinburgh)" in _iip3["credit"]
    assert len(_iip3["reassessed_hurdles"]) == 4
    assert [h["one_design_design_a"] for h in _iip3["reassessed_hurdles"]] == [0.88, 0.78, 0.57, 0.69]
    assert [h["one_design_design_b"] for h in _iip3["reassessed_hurdles"]] == [0.88, 0.79, 0.64, 0.70]
    assert [h["one_design_design_x"] for h in _iip3["reassessed_hurdles"]] == [0.88, 0.82, 0.73, 0.69]
    assert [h["three_designs_parallel"] for h in _iip3["reassessed_hurdles"]] == [0.88, 0.87, 0.86, 0.80]
    assert abs(_iip3["multiplied_conditional_totals"]["design_a_exact"] - 0.26996112) < 1e-8
    assert abs(_iip3["multiplied_conditional_totals"]["design_b_exact"] - 0.3114496) < 1e-8
    assert abs(_iip3["multiplied_conditional_totals"]["design_x_exact"] - 0.36346992) < 1e-8
    assert abs(_iip3["multiplied_conditional_totals"]["three_designs_exact"] - 0.5267328) < 1e-8
    assert _iip3["multiplied_conditional_totals"]["prior_2026_11_15_dec098_baseline"]["three_designs_rounded_pct"] == 51
    assert _iip3["multiplied_conditional_totals"]["prior_2026_11_05_baseline"]["one_design_rounded_pct"] == 21
    assert _iip3["multiplied_conditional_totals"]["prior_2026_11_05_baseline"]["three_designs_rounded_pct"] == 37
    assert _iip3["section3_score_and_rank_measure_reconciliation"]["rank_column_measure"] == "overall_50_50_composite_score"
    assert _iip3["new_open_question"]["id"] == "OQ-016"
    if _r9_rank_atlas is not None:
        assert _r9_rank_atlas["odds_composition_and_profile_scope_audit_iip_0003"]["proposal_id"] == "IIP-0003"
        assert any(v.get("open_question_id") == "OQ-016" for v in _r9_rank_atlas["step7_high_value_of_information_gaps"])
    assert "IIP-0003" in _r9_s2_bp_md and "OQ-016" in _r9_s2_bp_md and "Linus Schumacher" in _r9_s2_bp_md
    assert "(Gate 5 oxygen bridge)" in _r9_s2_bp_md

    # 40. DEC-085..DEC-088 / IIP-0007..IIP-0010 (OQ-017, Xu 2024 Tier_2B + Hu primary-islet CD47 arm, Design X scope mapping, and per-formulation CaO2 t_1/2 sensitivity):
    _xu_rec = json.loads((ROOT / "data/extracted/38508139.json").read_text(encoding="utf-8"))
    assert _xu_rec["study_metadata"]["max_biological_fidelity_tier"] == "Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture"
    assert len(_xu_rec["arms"]) == 4
    assert all(a["gate1_substrate_and_differentiation"]["cell_source_class"] == "NonIslet_Immune_Endothelial_or_Hematopoietic_Model" for a in _xu_rec["arms"])
    _hu_rec = json.loads((ROOT / "data/extracted/37156915.json").read_text(encoding="utf-8"))
    assert len(_hu_rec["arms"]) == 5
    assert any(
        a["arm_id"] in ("hu2024_human_primary_islets_cd47_threshold", "pmid37156915_arm5")
        and a["gate1_substrate_and_differentiation"]["cell_source_class"] == "Primary_Human_Islet"
        for a in _hu_rec["arms"]
    )
    if _r9_rank_atlas is not None:
        assert any(v.get("open_question_id") == "OQ-017" and "IIP-0007" in v.get("why_needed", "") for v in _r9_rank_atlas["step7_high_value_of_information_gaps"])
        _gy9 = next(r for r in _r9_rank_atlas["top_15_combinatorial_search_results"][0]["triggered_graveyard_rules"] if r["graveyard_id"].startswith("GRAVEYARD_009"))
        assert "38508139" not in _gy9["reason"]
    _iip9_map = _r9_s6_l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["design_x_100um_lamina_transport_qualification_v80_018"]["blueprint_a_outer_band_115p5um_vs_100um_cluster_cap_reconciliation_sonnet55_r2"]["representation_and_score_scope_mapping_iip_0009"]
    assert _iip9_map["iip_id"] == "IIP-0009" and set(_iip9_map["open_question_ids"]) == {"OQ-014", "OQ-015"}
    _iip10_sens = _r9_s6_l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["formulation_specific_cao2_half_life_sensitivity_iip_0010"]
    assert _iip10_sens["iip_id"] == "IIP-0010" and _iip10_sens["open_question_id"] == "OQ-005"
    assert len(_iip10_sens["grid_rows"]) == 84
    _cad_p5 = _r9_s6_l4["part_a_exact_oxygen_boundary_value_and_volume_algebra"]["non_overlapping_phase_volume_and_cad_ledger_7cm3_ultracode_v80_18"]
    assert _cad_p5["cao2_crystal_particle_size_distribution_um"] == "unknown_pending_stage0_measurement"
    assert _cad_p5["pdms_cao2_composite_microbead_diameter_distribution_um"] == "unknown_pending_stage0_measurement"
    assert "IIP-0008" in _r9_s2_bp_md and "IIP-0009" in _r9_s2_bp_md and "IIP-0010" in _r9_s2_bp_md

    # 41. DEC-090 / IIP-0011 (Locus 3/4 cis/UCOE guardrails, SLC30A8 hold, iCasp9 1BL4/1JXQ/1NW9 structural audit, donor TYK2 sensitivity & donor IFIH1 provenance invariance):
    _r9_l2_iip11 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    assert _r9_l2_iip11["stage0_core_construct_and_secretion_guardrails_iip_0011"]["proposal_id"] == "IIP-0011"
    _r9_l3_iip11 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    _icasp9_aud = _r9_l3_iip11["icasp9_locus4_structural_reference_audit"]
    assert _icasp9_aud["proposal_id"] == "IIP-0011" and _icasp9_aud["qualification_status"] == "FULL_FUSION_UNRESOLVED"
    assert _icasp9_aud["fkbp12_f36v_domain_reference_1BL4"]["identical_residues"] == 107
    assert _icasp9_aud["caspase9_catalytic_homodimer_reference_1JXQ"]["all_model1_polymer_heavy_atoms_including_zero_occupancy"]["heavy_atom_contact_pairs_le_4p5A"] == 674
    assert _icasp9_aud["caspase9_catalytic_homodimer_reference_1JXQ"]["all_model1_polymer_heavy_atoms_including_zero_occupancy"]["chain_A_interface_residues_insertion_code_aware"] == 44
    assert _icasp9_aud["caspase9_catalytic_homodimer_reference_1JXQ"]["all_model1_polymer_heavy_atoms_including_zero_occupancy"]["chain_B_interface_residues_insertion_code_aware"] == 42
    assert _icasp9_aud["caspase9_catalytic_homodimer_reference_1JXQ"]["positive_occupancy_only_polymer_heavy_atoms"]["heavy_atom_contact_pairs_le_4p5A"] == 667
    assert _icasp9_aud["caspase9_xiap_bir3_inhibitory_complex_reference_1NW9"]["heavy_atom_contact_pairs_le_4p5A"] == 363
    assert _icasp9_aud["caspase9_xiap_bir3_inhibitory_complex_reference_1NW9"]["cleaved_atpf_316_319_heavy_atom_contact_pairs_le_4p5A"] == 145
    assert _icasp9_aud["caspase9_xiap_bir3_inhibitory_complex_reference_1NW9"]["shared_canonical_casp9_interface_residues_with_1JXQ_homodimer"]["vs_1JXQ_chain_A_count"] == 17
    assert _icasp9_aud["caspase9_xiap_bir3_inhibitory_complex_reference_1NW9"]["shared_canonical_casp9_interface_residues_with_1JXQ_homodimer"]["vs_1JXQ_chain_B_count"] == 19
    _tyk2_sens = _r9_s6_l5["part_b_personalized_patient_compiler_and_10k_population_trial"]["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["donor_seedline_tyk2_genotype_sensitivity_summary"]
    assert _tyk2_sens["proposal_id"] == "IIP-0011" and _tyk2_sens["assigned_survival_invariance_verified"] is True
    _tyk2_rows = {r["donor_tyk2_rs34536443_genotype"]: r for r in _tyk2_sens["donor_genotype_comparison_table"]}
    assert _tyk2_rows["G/G"]["kb_complexity_counts"] == {"K4_B0": 4286, "K4_B1": 3347, "K5_B0": 1340, "K5_B1": 1027}
    assert _tyk2_rows["G/C"]["kb_complexity_counts"] == {"K4_B0": 4286, "K4_B1": 3347, "K5_B0": 1340, "K5_B1": 1027}
    assert _tyk2_rows["C/C"]["kb_complexity_counts"] == {"K4_B0": 7633, "K4_B1": 0, "K5_B0": 2367, "K5_B1": 0}
    assert _tyk2_rows["C/C"]["recipients_prescribed_tyk2_prime_editing"] == 0
    from scripts import run_layer5_personalized_patient_compiler as _l5_mod
    _pat_base = {
        "patient_id": "TEST-IFIH1-PROV",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.20,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 85.0,
        "sirpa_genotype": "v1/v2",
        "tyk2_rs34536443_genotype": "G/G",
        "ifih1_rs1990760_genotype": "T/T",
        "ptpn22_rs2476601_genotype": "G/G",
        "recipient_bmi_kg_m2": 22.0,
    }
    _comp_none = _l5_mod.compile_personalized_patient_blueprint(dict(_pat_base))
    assert _comp_none["host_modifier_genotypes_and_dosing"]["optional_field_input_provenance"]["donor_ifih1_rs1990760_genotype"] == "DEFAULT_IMPUTED_NOT_MEASURED"
    assert _comp_none["host_modifier_genotypes_and_dosing"]["donor_ifih1_rs1990760_provenance"]["input_status"] == "UNMEASURED_PROVENANCE_ONLY"
    _comp_unm = _l5_mod.compile_personalized_patient_blueprint({**_pat_base, "donor_ifih1_rs1990760_genotype": "UNMEASURED"})
    assert _comp_unm["host_modifier_genotypes_and_dosing"]["donor_ifih1_rs1990760_provenance"]["input_status"] == "UNMEASURED_PROVENANCE_ONLY"
    for _d_ifih1 in ("C/C", "C/T", "T/C", "T/T"):
        _comp_var = _l5_mod.compile_personalized_patient_blueprint({**_pat_base, "donor_ifih1_rs1990760_genotype": _d_ifih1})
        assert _comp_var["host_modifier_genotypes_and_dosing"]["optional_field_input_provenance"]["donor_ifih1_rs1990760_genotype"] == "MEASURED_PATIENT_INPUT"
        assert _comp_var["host_modifier_genotypes_and_dosing"]["donor_ifih1_rs1990760_provenance"]["input_status"] == "MEASURED_DONOR_INPUT"
        assert _comp_var["editing_complexity_and_manufacturing_accounting_v80_009"]["k_dsb_or_indel_chromosomal_loci"] == _comp_unm["editing_complexity_and_manufacturing_accounting_v80_009"]["k_dsb_or_indel_chromosomal_loci"]
        assert _comp_var["editing_complexity_and_manufacturing_accounting_v80_009"]["b_prime_or_base_edit_loci"] == _comp_unm["editing_complexity_and_manufacturing_accounting_v80_009"]["b_prime_or_base_edit_loci"]
        assert _comp_var["product_variant_id"] == _comp_unm["product_variant_id"]
        assert _l5_mod.evaluate_regimen_survival_probability({**_pat_base, "donor_ifih1_rs1990760_genotype": _d_ifih1}, _comp_var) == _l5_mod.evaluate_regimen_survival_probability(_pat_base, _comp_unm)
    try:
        _l5_mod.compile_personalized_patient_blueprint({**_pat_base, "donor_ifih1_rs1990760_genotype": "G/G"})
        raise AssertionError("Expected ValueError for invalid donor_ifih1_rs1990760_genotype")
    except ValueError:
        pass
    assert "IIP-0011" in _r9_s2_bp_md and "DEC-090" in _r9_s2_bp_md

    # =========================================================================
    # 42. DEC-091 / IIP-0012 POST-ADJUNCT THREE-HISTORY RECHALLENGE CONTRACT,
    #     SA-FHL-1 TRUTH-IN-LABELING, UNCONFOUNDED CD55 SATELLITE & HELD
    #     STAGE-0 CHALLENGER FALSIFIER GATES
    # =========================================================================
    _dt_v1_iip12 = json.loads((ROOT / "data/protocols/stage0_decision_tree_v1.json").read_text(encoding="utf-8"))
    assert "post_adjunct_depletion_three_history_rechallenge_and_cd55_gate_iip_0012" in _dt_v1_iip12
    assert "stage0_held_challenger_falsifier_gates_iip_0012" in _dt_v1_iip12
    _held_12 = _dt_v1_iip12["stage0_held_challenger_falsifier_gates_iip_0012"]
    assert _held_12["initial_state"] == "HOLD_FOR_PROSPECTIVE_FALSIFIER"
    assert len(_held_12["held_challenger_families"]) == 4

    assert "stage0_post_adjunct_and_held_challenger_falsifier_gates_iip_0012" in _r9_l2_iip11
    assert (
        _r9_l2_iip11["stage0_post_adjunct_and_held_challenger_falsifier_gates_iip_0012"]["core_scoring_impact"]
        == "ZERO_LAYER1_SCORE_CHANGE"
    )
    def _find_nested_dict_key(node, target_key):
        if isinstance(node, dict):
            if target_key in node:
                return node[target_key]
            for val in node.values():
                hit = _find_nested_dict_key(val, target_key)
                if hit is not None:
                    return hit
        elif isinstance(node, list):
            for item in node:
                hit = _find_nested_dict_key(item, target_key)
                if hit is not None:
                    return hit
        return None

    _cd55_l3 = _find_nested_dict_key(
        _r9_l3_iip11,
        "ectopic_cd55_scr14_unconfounded_causal_attribution_satellite_iip_0012",
    )
    assert _cd55_l3 is not None and _cd55_l3["decision_id"] == "DEC-091"

    # Active Arm 3 default still fails when activated CD4 FasL apoptosis < 80%:
    _arm3_low_fasl = {
        **full_ext_payload,
        "obs_indirect_apc_dq2_dq8_cd4_stimulation_index": 1.1,
        "obs_indirect_apc_cd4_ifng_si": 1.2,
        "obs_gate0d_activated_cd4_fasl_apoptosis_pct": 12.0,
        "obs_gate0d_transwell_beta_annexin_v_pct": 6.0,
    }
    _eval_arm3_low = wc_mod.evaluate_stage0_decision_tree(_arm3_low_fasl)
    assert (
        _eval_arm3_low["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_GATE0D_INDIRECT_APC_DQ2_DQ8_CD4_AUTOIMMUNE_STIMULATION_EXCEEDED"
    )

    # Verified Arm 2 exempts <80% FasL CD4 apoptosis, enforces SA-FHL-1 truth-in-labeling, and clears shared 6-endpoint post-depletion bounds:
    _arm2_pass = {
        **_arm3_low_fasl,
        "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED",
        "adjunct_depletion_mode": "NATURAL_AGING_DECAY",
        "adjunct_sa_fhl1_status": "RETAINED_ACTIVE",
        "adjunct_immune_responder_history": "PAIRED_BOTH",
        "obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified": True,
        "obs_arm2_sa_thbd_bioactivity_depleted_verified": True,
        "obs_arm2_contemporaneous_sham_handling_matched": True,
        "obs_arm3_ongoing_adjunct_bioactivity_verified": True,
        "obs_arm2_post_depletion_viable_beta_recovery_pct": 92.0,
        "obs_arm2_post_depletion_perifusion_gsis_si": 3.4,
        "obs_arm2_post_depletion_nhs_c5b9_mac_pct": 2.5,
    }
    _eval_arm2_pass = wc_mod.evaluate_stage0_decision_tree(_arm2_pass)
    assert (
        _eval_arm2_pass["resolved_architecture_assignments"]["arm2_post_adjunct_depletion_qualification_status"]
        == "ARM2_POST_ADJUNCT_DEPLETION_RECHALLENGE_CLEARED"
    )
    assert (
        _eval_arm2_pass["resolved_architecture_assignments"]["adjunct_depletion_truth_in_labeling_status"]
        == "PARTIAL_FASL_THBD_DEPLETED_SA_FHL1_RETAINED_ACTIVE_NOT_COMPLETE_ADJUNCT_DEPLETION"
    )

    # Incomplete Arm 2 verification & rebound injury holds + non-pooling of NATURAL_AGING_DECAY vs PHYSICAL_WASHOUT:
    _eval_arm2_unver = wc_mod.evaluate_stage0_decision_tree(
        {**_arm2_pass, "obs_arm2_sa_thbd_bioactivity_depleted_verified": False}
    )
    assert (
        _eval_arm2_unver["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_GATE0D_ARM2_ADJUNCT_DEPLETION_OR_SHAM_VERIFICATION_INCOMPLETE"
    )
    _eval_arm2_rebound = wc_mod.evaluate_stage0_decision_tree(
        {**_arm2_pass, "obs_arm2_post_depletion_nhs_c5b9_mac_pct": 9.5}
    )
    assert (
        _eval_arm2_rebound["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_GATE0D_ARM2_POST_ADJUNCT_DEPLETION_REBOUND_INJURY_OR_FUNCTION_LOSS"
    )
    _cal_dm1 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP12-DM1",
            "mcb_clone_id": "CLONE-IIP12",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED",
            "adjunct_depletion_mode": "NATURAL_AGING_DECAY",
            "obs_stage7_yield_ratio": 0.91,
            "obs_stage7_yield_se": 0.03,
        },
        is_dry_run=False,
    )
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-IIP12-DM2",
                "mcb_clone_id": "CLONE-IIP12",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED",
                "adjunct_depletion_mode": "PHYSICAL_WASHOUT",
                "obs_cd47_fold": 4.2,
                "obs_cd47_fold_se": 0.2,
            },
            False,
            _cal_dm1,
        )
        raise AssertionError("Expected ValueError when pooling NATURAL_AGING_DECAY and PHYSICAL_WASHOUT batches")
    except ValueError as _exc_dm:
        assert "adjunct_depletion_mode" in str(_exc_dm)
    assert "IIP-0012" in _r9_s2_bp_md and "DEC-091" in _r9_s2_bp_md

    # =========================================================================
    # 43. DEC-092 / IIP-0015: TWO-STRAND KNOCK-IN SPLICE CENSUS, PAM-SHIELDED
    #     MODE-B CHALLENGER MANIFEST, SEPARATED 505/500-AA LOCUS-1 SEAMS, AND
    #     RECOVERY-CALIBRATED UMI / 90% TOST / LOQ / NON-POOLING GATES
    # =========================================================================
    import hashlib as _hl_iip15

    assert (
        _hl_iip15.sha256(
            (ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta").read_bytes()
        ).hexdigest()
        == "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"
    )
    _l2_iip15 = _r9_l2_iip11["stage0_two_strand_splice_census_and_challenger_manifest_iip_0015"]
    assert _l2_iip15["proposal_id"] == "IIP-0015" and _l2_iip15["decision_id"] == "DEC-092"
    assert len(_l2_iip15["validator_fail_closed_adversarial_selftests"]) == 11
    assert all(
        row["rejected_as_expected"] is True
        for row in _l2_iip15["validator_fail_closed_adversarial_selftests"]
    )
    _ref_hash_by_id = {
        row["construct_id"]: row["donor_sha256"][:8]
        for row in _l2_iip15["qualified_mode_b_reference_derivations"]
    }
    for _cid, _pfx in (
        ("S3_L3_CAN_PDL1_SYN2_4090", "3cb5ce72"),
        ("S3_L3_CAN_PDL1_SYN4_4090", "ccd86d6a"),
        ("S3_L3_CAN_BGH_ONLY_4090", "304a08d9"),
        ("S3_L3_SF_PDL1_SYN4_3799", "a678fade"),
        ("S3_L3_SF_SYN4_1xFII80_3879", "6d594329"),
        ("S3_L4_SPLICE_230C_NOSPA_3457", "3bc15171"),
        ("S3_L4_CONVERGENT_2STRAND_230C_3457", "70424667"),
        ("S3_L1_ASTAR_STOP_ONLY_4143", "aff2ec0d"),
    ):
        assert _ref_hash_by_id[_cid] == _pfx
    _l4_mat_arms = _l2_iip15["pam_seed_shielded_synthesis_ready_manifest"][
        "locus4_aavs1_pam_shielded_and_unconfounded_2x2_orientation_matrix"
    ]["arms"]
    assert len(_l4_mat_arms) >= 5
    _conv_row = next(
        r for r in _l4_mat_arms if r.get("alias_construct_id") == "S3_L4_CONVERGENT_2STRAND_230C_3457"
    )
    assert "70424667" in json.dumps(_conv_row) and "977495f1" in json.dumps(_conv_row)
    _ab_l1 = _l2_iip15["aberrant_locus1_sct_to_b2m_exon2_cryptic_splice_audit"]
    assert _ab_l1["embodiment_a_star_505aa_sct1514_audit"]["sct_precursor_aa"] == 505
    assert _ab_l1["embodiment_a_star_505aa_sct1514_audit"]["total_chimeric_polypeptide_aa"] == 537
    assert _ab_l1["arm1b_and_arm1c_500aa_sct1499_audit"]["sct_precursor_aa"] == 500
    assert _ab_l1["arm1b_and_arm1c_500aa_sct1499_audit"]["total_chimeric_polypeptide_aa"] == 532
    assert len(
        _ab_l1["embodiment_a_star_505aa_sct1514_audit"]["appended_b2m_exon2_plus1_frame_tail_32aa"]
    ) == 32
    _rec_l1 = _l2_iip15["portfolio_chimeric_boundary_reconciliation"]
    assert _rec_l1["nine_boundaries_8_11mer"] == 326
    assert _rec_l1["nine_boundaries_15mer_full_polypeptide"] == 112
    assert _rec_l1["ten_boundaries_with_hlac_l2_alpha1_15mer_full_polypeptide"] == 126
    assert "two_strand_knockin_splice_challenger_and_umi_isoform_gate_iip_0015" in json.loads(
        (ROOT / "data/protocols/stage0_decision_tree_v1.json").read_text(encoding="utf-8")
    )

    _chall_base_pass = {
        **full_ext_payload,
        "knockin_splice_challenger_arm": "S1_L3_CAN_SYN4_PAM_SHIELDED_4090",
        "knockin_locus4_orientation": "NOT_APPLICABLE",
        "obs_challenger_pam_seed_shield_coconversion_verified": True,
        "obs_challenger_ivt_rna_recovery_standard_verified": True,
        "obs_challenger_host_isoforms_above_loq_verified": True,
        "obs_challenger_sv40pa_and_host_3prime_cleavage_checked": True,
        "obs_challenger_ivt_relative_recovery_eta": 0.90,
        "obs_challenger_ivt_relative_recovery_se": 0.02,
        "obs_challenger_umi_molecules_profiled": 4200.0,
        "obs_challenger_umi_molecules_profiled_se": 1.0,
        "obs_challenger_aberrant_splice_umi_count_k": 0.0,
        "obs_challenger_aberrant_splice_umi_count_se": 1e-6,
        "obs_challenger_host_exon_splice_retention_ratio": 1.01,
        "obs_challenger_host_exon_splice_retention_se": 0.03,
        "obs_challenger_dsrna_j2_or_isg_fold_wt": 1.08,
        "obs_challenger_dsrna_j2_or_isg_fold_se": 0.05,
    }
    _eval_chall_pass = wc_mod.evaluate_stage0_decision_tree(_chall_base_pass)
    assert (
        _eval_chall_pass["resolved_architecture_assignments"][
            "knockin_splice_challenger_qualification_status"
        ]
        == "CHALLENGER_SPLICE_AND_RECOVERY_ADJUSTED_UMI_PANEL_CLEARED"
    )
    assert (
        _eval_chall_pass["resolved_architecture_assignments"][
            "challenger_effective_umi_molecules"
        ]
        == 3780.0
    )

    # Unshielded reference derivation holds from synthesis/release
    _eval_unsh = wc_mod.evaluate_stage0_decision_tree(
        {**_chall_base_pass, "knockin_splice_challenger_arm": "S1_L3_CAN_SYN4_4090"}
    )
    assert (
        _eval_unsh["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_CHALLENGER_SPLICE_UNSHIELDED_OR_PAM_SEED_COCONVERSION_UNVERIFIED"
    )

    # 10% under-recovery (eta=0.90) on nominal N=3700 -> N_UMI_eff=3330 < 3689 holds!
    _eval_low_umi = wc_mod.evaluate_stage0_decision_tree(
        {**_chall_base_pass, "obs_challenger_umi_molecules_profiled": 3700.0}
    )
    assert (
        _eval_low_umi["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_CHALLENGER_SPLICE_EFFECTIVE_UMI_BELOW_3689_RECOVERY_ADJUSTED_MINIMUM"
    )

    # TOST outside [0.80, 1.25], LOQ=False, SV40pA unchecked, and dsRNA > 1.5x holds
    assert (
        wc_mod.evaluate_stage0_decision_tree(
            {**_chall_base_pass, "obs_challenger_host_exon_splice_retention_ratio": 0.81, "obs_challenger_host_exon_splice_retention_se": 0.04}
        )["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_CHALLENGER_SPLICE_HOST_EXON_RETENTION_OUTSIDE_90PCT_TOST_0P80_1P25"
    )
    assert (
        wc_mod.evaluate_stage0_decision_tree(
            {**_chall_base_pass, "obs_challenger_host_isoforms_above_loq_verified": False}
        )["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_CHALLENGER_SPLICE_HOST_ISOFORM_BELOW_LOQ_OR_UNVERIFIED"
    )
    assert (
        wc_mod.evaluate_stage0_decision_tree(
            {**_chall_base_pass, "obs_challenger_sv40pa_and_host_3prime_cleavage_checked": False}
        )["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_CHALLENGER_SPLICE_SV40PA_AND_HOST_3PRIME_CLEAVAGE_UNCHECKED"
    )
    assert (
        wc_mod.evaluate_stage0_decision_tree(
            {**_chall_base_pass, "obs_challenger_dsrna_j2_or_isg_fold_wt": 1.48, "obs_challenger_dsrna_j2_or_isg_fold_se": 0.05}
        )["resolved_architecture_assignments"]["lot_maturation_status"]
        == "HOLD_CHALLENGER_SPLICE_CONVERGENT_DSRNA_J2_OR_ISG_EXCEEDS_1P5X_CEILING"
    )

    # Cross-batch non-pooling guard across distinct challenger arms or Locus-4 orientations
    _cal_ch1 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP15-C1",
            "mcb_clone_id": "CLONE-IIP15",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "knockin_splice_challenger_arm": "S3_L4_SPLICE_230C_3457",
            "knockin_locus4_orientation": "CO_DIRECTIONAL_MINUS",
            "obs_stage7_yield_ratio": 0.92,
            "obs_stage7_yield_se": 0.03,
        },
        is_dry_run=False,
    )
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-IIP15-C2",
                "mcb_clone_id": "CLONE-IIP15",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                "knockin_splice_challenger_arm": "S3_L4_CONVERGENT_2STRAND_230C_3457",
                "knockin_locus4_orientation": "CONVERGENT_PLUS",
                "obs_cd47_fold": 4.2,
                "obs_cd47_fold_se": 0.2,
            },
            False,
            _cal_ch1,
        )
        raise AssertionError("Expected ValueError when pooling distinct knockin_splice_challenger_arm / orientation batches")
    except ValueError as _exc_ch:
        assert "knockin_splice_challenger_arm" in str(_exc_ch) or "knockin_locus4_orientation" in str(_exc_ch)
    assert "IIP-0015" in _r9_s2_bp_md and "DEC-092" in _r9_s2_bp_md

    # 31. IIP-0016 / DEC-094: Reconciled Lin 2023 HLA-E*01:03 6-leader SPR table, 20-construct Mode-B S10 comparator panel, and C1-C4 v3 contract repairs
    import hashlib as _hl_iip16
    assert (
        _hl_iip16.sha256((ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta").read_bytes()).hexdigest()
        == "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"
    )
    _s10_fasta_bytes = (ROOT / "data/constructs/islet1_locus1_s10_modeb_comparators_v0.8.0.fasta").read_bytes()
    assert _hl_iip16.sha256(_s10_fasta_bytes).hexdigest() == "a589ea72144cf1ad2cfa82e9b08c22e075746c63a9264e276a3bd4ec3f896854"

    _l2_s10 = _r9_l2_iip11["stage0_s10_hlae_comparator_panel_manifest_iip_0016"]
    assert _l2_s10["proposal_id"] == "IIP-0016" and _l2_s10["decision_id"] == "DEC-094"
    assert len(_l2_s10["verified_20_constructs"]) == 20
    _s10_by_id = {c["construct_id"]: c for c in _l2_s10["verified_20_constructs"]}
    assert _s10_by_id["S10_CTRL_HARDWT_VMAPRTLVL_ASTAR_4143"]["donor_sha256"].startswith("b5e0bee49a1c")
    assert _s10_by_id["S10_CTRL_HARDWT_VMAPRTLVL_ARM1B_3974"]["donor_sha256"].startswith("e3b3f04321dc")
    assert _s10_by_id["S10_CTRL_HARDWT_VMAPRTLVL_ASTAR_4143"]["sct_precursor_aa_diffs_vs_wt_parent"] == []
    assert _s10_by_id["S10_CTRL_HARDWT_VMAPRTLVL_ASTAR_4143"]["cargo_acceptors_delta_vs_parent"] == -2
    assert _s10_by_id["S10_CHAL5B_SP6B_VMAPRTVLL_D227K_T228A_ASTAR_4143"]["sct_precursor_aa_diffs_vs_wt_parent"] == ["L27V", "V28L", "D395K", "T396A"]
    assert _s10_by_id["S10_CHAL5B_SP6B_VMAPRTVLL_D227K_T228A_ARM1B_3974"]["sct_precursor_aa_diffs_vs_wt_parent"] == ["L27V", "V28L", "D390K", "T391A"]
    for _c16 in _l2_s10["verified_20_constructs"]:
        assert _c16["created_donors_cv70_or_u1_7"] == 0
        assert _c16["created_new_acceptor_coordinates"] == 0
        assert _c16["junction_8_to_11mer_windows"] == 168
        assert _c16["junction_15mer_windows"] == 60

    _l3_s10 = _r9_l3_iip11["complexes"][0]["s10_hlae_nkg2a_nkg2c_selectivity_and_14pdb_audit_iip_0016"]
    assert len(_l3_s10["lin_2023_pmid_37264229_hlae_0103_complete_6_leader_spr_table"]) == 6
    assert all(r["kc_over_ka_point_estimate"] < 4.0 for r in _l3_s10["lin_2023_pmid_37264229_hlae_0103_complete_6_leader_spr_table"])
    assert len(_l3_s10["pdb_14_structure_reconciled_audit"]["peptide_bound_13_structures"]) == 13
    assert "3BDW" in _l3_s10["pdb_14_structure_reconciled_audit"]["unliganded_receptor_structure_excluded_from_chords"]

    # C1: Unmeasured S10 on core-14 emits UNMEASURED_EXTENDED_S10_PANEL
    _dt_core14_s10 = wetlab.evaluate_stage0_decision_tree(core_14_panel)
    assert _dt_core14_s10["switch_decisions"]["S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A"] == "UNMEASURED_EXTENDED_S10_PANEL"
    assert _dt_core14_s10["canonical_json_switch_id_map"]["JSON_S10_CMV_NKG2C_Adaptive_NK_Selectivity_Gate_0A"] == "UNMEASURED_EXTENDED_S10_PANEL"

    # C2: Raw mean 7.5% with SE=1.0 -> upper95=9.1449% > 8.0% triggers local S10
    _dt_c2 = wetlab.evaluate_stage0_decision_tree(
        {
            **core_14_panel,
            "obs_nkg2c_adaptive_nk_lysis_pct": 7.5,
            "obs_nkg2c_adaptive_nk_lysis_se": 1.0,
        }
    )
    assert _dt_c2["switch_decisions"]["S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A"] == "STRATIFY_CMV_NKG2C_RECIPIENTS_AND_EVALUATE_SCT_SELECTIVITY_MUTANT_COUPON"

    # C3 & C4: Calibration controls blocked; spoofed legacy coupon=True rejected for complete release; valid v3 challenger advances to READY
    _sp6b_sha = wc.STAGE0_S10_V3_CONSTRUCT_REGISTRY["S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143"]["donor_sha256"]
    _fl_sha = wc.STAGE0_S10_V3_CONSTRUCT_REGISTRY["S10_CTRL_POS_VMAPRTLFL_ASTAR_4143"]["donor_sha256"]
    _dt_cal_blocked = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_payload,
            "s10_construct_id": "S10_CTRL_POS_VMAPRTLFL_ASTAR_4143",
            "s10_construct_donor_sha256": _fl_sha,
            "current_lot_id": "LOT-FL-1",
            "s10_qualified_lot_id": "LOT-FL-1",
            "s10_is_replacement_challenger_lot": True,
            "s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct": True,
        }
    )
    assert _dt_cal_blocked["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CALIBRATION_CONTROL_NEVER_ELIGIBLE_FOR_CLINICAL_PROMOTION"
    assert _dt_cal_blocked["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    _dt_spoof_legacy = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_payload,
            "obs_hlae_nkg2c_over_nkg2a_kd_ratio": 2.68,
            "obs_hlae_nkg2c_over_nkg2a_kd_ratio_se": 0.10,
            "s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct": True,
        }
    )
    assert _dt_spoof_legacy["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "CANDIDATE_FALLBACK_PROMOTED_PENDING_REPLACEMENT_QUALIFICATION"

    _v3_valid_payload = {
        **full_ext_payload,
        "obs_hlae_nkg2c_over_nkg2a_kd_ratio": 2.68,
        "obs_hlae_nkg2c_over_nkg2a_kd_ratio_se": 0.10,
        "s10_construct_id": "S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143",
        "s10_construct_donor_sha256": _sp6b_sha,
        "current_lot_id": "LOT-SP6B-01",
        "s10_qualified_lot_id": "LOT-SP6B-01",
        "s10_hlae_floor_mode": "LOT_RELEASE_25K",
        "s10_is_replacement_challenger_lot": True,
        "s10_race_5p_3p_junction_verified": True,
        "obs_hlae_hour0_molecules": 28000.0,
        "obs_hlae_hour0_molecules_se": 500.0,
        "obs_hypoxic_72h_hlae_p10_single_cell_molecules": 18000.0,
        "obs_hypoxic_72h_hlae_p10_single_cell_se": 400.0,
        "obs_hlae_sct_nkg2a_kd_uM": 2.98,
        "obs_hlae_sct_nkg2a_kd_se_uM": 0.12,
        "obs_hlae_sct_nkg2c_kd_uM": 8.00,
        "obs_hlae_sct_nkg2c_kd_se_uM": 0.20,
        "obs_nkg2a_pos_resting_nk_lysis_pct": 3.1,
        "obs_nkg2a_pos_resting_nk_lysis_se": 0.6,
        "obs_nkg2a_pos_ifng24h_nk_lysis_pct": 3.5,
        "obs_nkg2a_pos_ifng24h_nk_lysis_se": 0.6,
        "obs_nkg2a_pos_resting_cd107a_pct": 4.2,
        "obs_nkg2a_pos_resting_cd107a_se": 0.6,
        "obs_nkg2a_pos_ifng24h_cd107a_pct": 4.5,
        "obs_nkg2a_pos_ifng24h_cd107a_se": 0.6,
        "obs_nkg2c_kir_pos_resting_lysis_pct": 4.1,
        "obs_nkg2c_kir_pos_resting_lysis_se": 0.6,
        "obs_nkg2c_kir_pos_ifng24h_lysis_pct": 4.6,
        "obs_nkg2c_kir_pos_ifng24h_lysis_se": 0.6,
        "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_pct": 4.4,
        "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_resting_lysis_se": 0.6,
        "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_pct": 4.9,
        "obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_ifng24h_lysis_se": 0.6,
        "obs_hcmv_pooled_3subset_resting_nk_lysis_pct": 4.5,
        "obs_hcmv_pooled_3subset_resting_nk_lysis_se": 0.6,
        "obs_hcmv_pooled_3subset_ifng24h_nk_lysis_pct": 5.0,
        "obs_hcmv_pooled_3subset_ifng24h_nk_lysis_se": 0.6,
        "obs_nkg2c_cd107a_degranulation_fold_vs_b2m_ko": 1.02,
        "obs_nkg2c_cd107a_degranulation_fold_se": 0.05,
        "obs_hcmv_pos_donor_count": 6.0,
        "obs_hcmv_pos_donor_count_se": 0.001,
        "obs_s10_nk_failing_donor_count": 0.0,
        "obs_s10_nk_failing_donor_count_se": 0.001,
        "obs_cmv_ul40_hlae_cd8_lysis_pct": 1.8,
        "obs_cmv_ul40_hlae_cd8_lysis_se": 0.4,
        "obs_ul40_reactive_donor_count": 6.0,
        "obs_ul40_reactive_donor_count_se": 0.001,
        "obs_ul40_ctl_failing_donor_count": 0.0,
        "obs_ul40_ctl_failing_donor_count_se": 0.001,
        "obs_nonself_dq2_dq8_indirect_apc_si": 1.12,
        "obs_nonself_dq2_dq8_indirect_apc_si_se": 0.08,
        "obs_nonself_host_class1_cross_pres_lysis_pct": 1.5,
        "obs_nonself_host_class1_cross_pres_lysis_se": 0.4,
    }
    _dt_v3_ready = wetlab.evaluate_stage0_decision_tree(_v3_valid_payload)
    assert _dt_v3_ready["switch_decisions"]["S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A"] == "CLEARED_S10_FUNCTIONAL_SELECTIVITY_SUB4_SCT_V3"
    assert _dt_v3_ready["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

    # Screening 15k mode holds complete release even when v3 screening passes; unrelated failure also holds
    _dt_v3_screen_only = wetlab.evaluate_stage0_decision_tree({**_v3_valid_payload, "s10_hlae_floor_mode": "GATE_0A_SCREENING_15K"})
    assert _dt_v3_screen_only["switch_decisions"]["S10_CMV_NKG2C_Adaptive_NK_Peptide_Selectivity_Gate_0A"] == "CLEARED_S10_FUNCTIONAL_SELECTIVITY_SCREENING_15K_V3"
    assert _dt_v3_screen_only["resolved_architecture_assignments"]["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    _dt_v3_unrelated_hold = wetlab.evaluate_stage0_decision_tree({**_v3_valid_payload, "obs_phagocytosis_index_pct": 9.0, "obs_phagocytosis_index_se": 0.5})
    assert _dt_v3_unrelated_hold["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # Cross-batch non-pooling guard on s10_construct_id
    _cal_s10_a = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP16-1",
            "mcb_clone_id": "CLONE-IIP16",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "s10_construct_id": "S10_CTRL_WT_VMAPRTLVL_ASTAR_4143",
            "obs_stage7_yield_ratio": 0.94,
            "obs_stage7_yield_se": 0.02,
        },
        is_dry_run=False,
    )
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-IIP16-2",
                "mcb_clone_id": "CLONE-IIP16",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                "s10_construct_id": "S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143",
                "obs_cd47_fold": 4.4,
                "obs_cd47_fold_se": 0.2,
            },
            False,
            _cal_s10_a,
        )
        raise AssertionError("Expected ValueError when pooling distinct s10_construct_id batches")
    except ValueError as _exc_s10:
        assert "s10_construct_id" in str(_exc_s10)
    assert "IIP-0016" in _r9_s2_bp_md and "DEC-094" in _r9_s2_bp_md

    # --- IIP-0017 (DEC-095) Species-Separated Proteome Census, TYK2 P1104A HLA-C OR-Hold Gate & AlphaMissense v3 Provenance ---
    _l2_iip17 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    _l3_iip17 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    _l5_iip17 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    _c17 = _l2_iip17["species_separated_proteome_subtracted_peptide_census_iip_0017"]
    assert _c17["decision_id"] == "DEC-095" and _c17["contribution_id"] == "IIP-0017"
    assert _c17["canonical_fasta_and_ranking_invariance"]["canonical_scenario_a_fasta_sha256"] == "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"
    assert _c17["rule_c3_and_reference_proteome_provenance"]["uniprot_human_swissprot_snapshot"]["fasta_gz_sha256"] == "40a8856b4e3f3620b11635143b501ef701558ed97309b0cdb4423ef09e4618d0"
    assert len(_c17["rule_c3_and_reference_proteome_provenance"]["mhcflurry_2_3_28_and_v2_3_0_weights_provenance"]["pinned_model_file_sha256s"]) == 10
    _iedb17 = _c17["rule_c3_and_reference_proteome_provenance"]["iedb_and_arber_2013_t2a_vs_p2a_reconciliation"]
    assert _iedb17["total_study_iedb_records"] == 44 and len(_iedb17["t2a_records"]) == 12
    assert sum(1 for r in _iedb17["t2a_records"] if r["source_organism_name"] == "Thosea asigna virus (Setothosea asigna virus)") == 4
    assert sum(1 for r in _iedb17["t2a_records"] if r["source_organism_name"] is None) == 8
    assert len(_c17["species_catalog_22_molecular_species"]["per_species_ledger"]) == 22
    _arch17 = _c17["causally_separated_architecture_census"]["architectures"]
    assert _arch17["DesignA_EmbodimentAstar_1xP2A_mature_only_WT_TYK2"]["unique_unmatched_8_11mer_count"] == 319
    assert _arch17["DesignA_EmbodimentAstar_1xP2A_mature_only_plus_optional_TYK2_P1104A"]["unique_unmatched_8_11mer_count"] == 357
    assert _arch17["DesignA_EmbodimentAstar_plus_existing_L3BFB_4078bp_0xP2A_mature_only_WT_TYK2"]["unique_unmatched_8_11mer_count"] == 235
    assert _arch17["DesignA_EmbodimentAstar_plus_existing_L3BFB_4078bp_0xP2A_mature_only_plus_optional_TYK2_P1104A"]["unique_unmatched_8_11mer_count"] == 273
    assert _arch17["DesignB2_TYK2_P1104A_plus_Canonical_P2A_Locus3"]["unique_unmatched_8_11mer_count"] == 243
    assert _arch17["DesignB2_TYK2_P1104A_plus_Existing_ScarFree_L3BFB_4078bp_or_DXL3_3799bp"]["unique_unmatched_8_11mer_count"] == 121
    assert len(_c17["tyk2_p1104a_32row_hlac_and_hlae_presentation_table"]["rows"]) == 32
    _am17 = _l3_iip17["alphamissense_v3_ccby4_provenance_iip_0017"]
    assert _am17["gating_authority"] == "NON_GATING_DESCRIPTIVE_ANNOTATION_ONLY"
    assert _am17["dataset_and_license"]["zenodo_v3_record_doi"] == "10.5281/zenodo.10813168"
    assert _am17["tyk2_exon23_p1104a_window_audit"]["exonic_missense_snv_count"] == 96
    assert _am17["tyk2_exon23_p1104a_window_audit"]["likely_pathogenic_gt_0564_count"] == 38
    assert _l5_iip17["alphamissense_v3_ccby4_provenance_iip_0017"]["dataset_and_license"]["zenodo_v3_record_doi"] == "10.5281/zenodo.10813168"

    _b2_valid_base = {
        **full_ext_payload,
        "product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE",
        "b2_tyk2_recipient_genotype_stratum": "G_G",
        "b2_tyk2_hlac_allele": "HLA-C*12:03",
        "b2_tyk2_saptkflel_allele_assigned_ms_detected": False,
        "b2_tyk2_peptide_pulsed_positive_control_verified": True,
        "b2_tyk2_hlac_blockade_abrogation_verified": True,
        "obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp": 0.8,
        "obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_se": 0.6,
        "obs_b2_tyk2_p1104a_primed_cd8_elispot_sfu_per_1e6": 12.0,
        "obs_b2_tyk2_p1104a_primed_cd8_elispot_se": 4.0,
        "obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp": 0.9,
        "obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_se": 0.6,
    }
    _dt_b2_clear = wetlab.evaluate_stage0_decision_tree(_b2_valid_base)
    assert _dt_b2_clear["resolved_architecture_assignments"]["design_b2_tyk2_p1104a_hlac_gate_status"] == "CLEARED_DESIGN_B2_TYK2_P1104A_STRATIFIED_HLAC_GATE"

    # MS negative + functional killing positive in C_C stratum -> still holds on NEOANTIGEN_HAZARD_DETECTED (OR-hold across all strata)
    _dt_b2_cc_hazard = wetlab.evaluate_stage0_decision_tree({
        **_b2_valid_base,
        "b2_tyk2_recipient_genotype_stratum": "C_C",
        "obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_pp": 4.8,
        "obs_b2_tyk2_p1104a_iec_cd8_lysis_delta_vs_bwt_se": 0.5,
    })
    assert _dt_b2_cc_hazard["resolved_architecture_assignments"]["design_b2_tyk2_p1104a_hlac_gate_status"] == "HOLD_DESIGN_B2_TYK2_P1104A_HLAC_NEOANTIGEN_HAZARD_DETECTED"
    assert _dt_b2_cc_hazard["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_DESIGN_B2_TYK2_P1104A_HLAC_NEOANTIGEN_HAZARD_DETECTED"

    # MS positive + functional killing negative -> holds on MS_DETECTED_OR_DISCORDANT
    _dt_b2_ms_pos = wetlab.evaluate_stage0_decision_tree({
        **_b2_valid_base,
        "b2_tyk2_saptkflel_allele_assigned_ms_detected": True,
    })
    assert _dt_b2_ms_pos["resolved_architecture_assignments"]["design_b2_tyk2_p1104a_hlac_gate_status"] == "HOLD_DESIGN_B2_TYK2_P1104A_HLAC_MS_DETECTED_OR_DISCORDANT_REQUIRES_STRATIFIED_REVIEW"

    # Missing positive control / blockade -> POSITIVE_CONTROL_OR_CAUSAL_RESTRICTION_UNVERIFIED
    _dt_b2_unver = wetlab.evaluate_stage0_decision_tree({
        **_b2_valid_base,
        "b2_tyk2_hlac_blockade_abrogation_verified": None,
    })
    assert _dt_b2_unver["resolved_architecture_assignments"]["design_b2_tyk2_p1104a_hlac_gate_status"] == "HOLD_DESIGN_B2_TYK2_P1104A_POSITIVE_CONTROL_OR_CAUSAL_RESTRICTION_UNVERIFIED"

    # Underpowered SE -> INCOMPLETE_OR_UNDERPOWERED_STRATIFIED_PANEL
    _dt_b2_underpow = wetlab.evaluate_stage0_decision_tree({
        **_b2_valid_base,
        "obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_pp": 0.1,
        "obs_b2_tyk2_p1104a_isogenic_cd8_lysis_delta_vs_bwt_se": 2.9,
    })
    assert _dt_b2_underpow["resolved_architecture_assignments"]["design_b2_tyk2_p1104a_hlac_gate_status"] == "HOLD_DESIGN_B2_TYK2_P1104A_INCOMPLETE_OR_UNDERPOWERED_STRATIFIED_PANEL"

    # Cross-batch non-pooling guard on b2_tyk2_recipient_genotype_stratum and b2_tyk2_hlac_allele
    _cal_b2_1 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP17-1",
            "mcb_clone_id": "CLONE-IIP17",
            "b2_tyk2_recipient_genotype_stratum": "G_G",
            "b2_tyk2_hlac_allele": "HLA-C*12:03",
            "obs_stage7_yield_ratio": 0.94,
            "obs_stage7_yield_se": 0.02,
        },
        is_dry_run=False,
    )
    for _bad_b2_meta, _expected_err_sub in (
        ({"b2_tyk2_recipient_genotype_stratum": "C_C", "b2_tyk2_hlac_allele": "HLA-C*12:03"}, "b2_tyk2_recipient_genotype_stratum"),
        ({"b2_tyk2_recipient_genotype_stratum": "G_G", "b2_tyk2_hlac_allele": "HLA-C*04:01"}, "b2_tyk2_hlac_allele"),
    ):
        try:
            wetlab.run_calibration(
                {
                    "batch_id": "BATCH-IIP17-2",
                    "mcb_clone_id": "CLONE-IIP17",
                    **_bad_b2_meta,
                    "obs_cd47_fold": 4.4,
                    "obs_cd47_fold_se": 0.2,
                },
                False,
                _cal_b2_1,
            )
            raise AssertionError(f"Expected ValueError when pooling {_expected_err_sub}")
        except ValueError as _exc_b2:
            assert _expected_err_sub in str(_exc_b2)
    _r9_s2_bp_md_17 = (ROOT / "ISLET-1-BLUEPRINT.md").read_text(encoding="utf-8")
    assert "IIP-0017" in _r9_s2_bp_md_17 and "DEC-095" in _r9_s2_bp_md_17

    # --- IIP-0018 (DEC-096) 6-Arm Locus 3 Multi-miR(IL1R1+FAS) Mode-B Suite, Editor-Scoped CBE Needles & Stage-0 Gate ---
    import hashlib

    _f18_path = ROOT / "data/constructs/islet1_locus3_multimir_modeb_comparators_v0.8.0.fasta"
    assert _f18_path.exists()
    assert hashlib.sha256(_f18_path.read_bytes()).hexdigest() == "22df54a27ef3f68b646b63dc036b11f5a66a9fcadb50378dfb8057e108c79ffa"
    assert hashlib.sha256((ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta").read_bytes()).hexdigest() == "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"
    _l2_iip18 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))["stage0_locus3_multimir_and_cbe_needle_manifest_iip_0018"]
    assert _l2_iip18["proposal_id"] == "IIP-0018" and _l2_iip18["decision_id"] == "DEC-096"
    _arms18 = {a["construct_id"]: a for a in _l2_iip18["verified_6_arm_constructs"]}
    assert len(_arms18) == 6
    assert _arms18["ARM1_L3_WT_PARENT_4078"]["donor_length_bp"] == 4078
    assert _arms18["ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398"]["donor_length_bp"] == 4398
    assert _arms18["ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398"]["packaged_aav6_length_bp"] == 4688
    assert _arms18["ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398"]["novel_8_to_11mers_vs_arm1"] == 0
    assert _arms18["ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398"]["novel_15mers_vs_arm1"] == 0
    assert _arms18["ARM6_L3_V2_INTRON_DUPLEX_IL1R1_FAS_C341S_TCA_4398"]["novel_8_to_11mers_vs_arm1"] == 38
    assert _arms18["ARM6_L3_V2_INTRON_DUPLEX_IL1R1_FAS_C341S_TCA_4398"]["novel_15mers_vs_arm1"] == 15
    assert len(_l2_iip18["standalone_cbe_comparator_registry_ensembl_116"]["candidate_guides_9_across_7_variants"]) == 9
    _l3_iip18 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))["hurdle3_receptor_and_compartment_reconciliation_iip_0018"]
    assert _l3_iip18["proposal_id"] == "IIP-0018" and _l3_iip18["decision_id"] == "DEC-096"
    _c4dep = _l3_iip18["il1_signaling_vs_decoy_crystal_contacts_4p5A"]["pdb_4dep_il1beta_il1r1_il1racp_pmid_22426547"]
    assert _c4dep["copy_1_chains_A_B_C"]["IL1beta_to_IL1R1"] == 504 and _c4dep["copy_1_chains_A_B_C"]["binary_complex_to_IL1RAcP"] == 250
    assert _c4dep["copy_2_chains_D_E_F"]["IL1beta_to_IL1R1"] == 618 and _c4dep["copy_2_chains_D_E_F"]["binary_complex_to_IL1RAcP"] == 352
    _c4msv = _l3_iip18["fasl_dcr3_c3_assembly_geometry_pdb_4msv"]["c3_symmetry_ca_triangle_edges_angstrom"]
    assert abs(_c4msv["dcr3_chainB_cys193_cterm_edge_A"] - 61.5796) < 1e-3 and abs(_c4msv["fasl_chainA_leu143_resolved_nterm_stalk_edge_A"] - 19.6245) < 1e-3

    _l3_mm_valid = {
        **full_ext_payload,
        "knockin_splice_challenger_arm": "ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398",
        "obs_challenger_pam_seed_shield_coconversion_verified": True,
        "obs_challenger_ivt_rna_recovery_standard_verified": True,
        "obs_challenger_host_isoforms_above_loq_verified": True,
        "obs_challenger_sv40pa_and_host_3prime_cleavage_checked": True,
        "obs_challenger_ivt_relative_recovery_eta": 0.96,
        "obs_challenger_ivt_relative_recovery_se": 0.01,
        "obs_challenger_umi_molecules_profiled": 5000.0,
        "obs_challenger_umi_molecules_profiled_se": 10.0,
        "obs_challenger_aberrant_splice_umi_count_k": 0.0,
        "obs_challenger_aberrant_splice_umi_count_se": 0.001,
        "obs_challenger_host_exon_splice_retention_ratio": 1.01,
        "obs_challenger_host_exon_splice_retention_se": 0.03,
        "obs_challenger_dsrna_j2_or_isg_fold_wt": 1.05,
        "obs_challenger_dsrna_j2_or_isg_fold_se": 0.05,
        "l3_multimir_pre_aav6_2x2_transient_kd_falsifier_cleared": True,
        "l3_multimir_low_ifn_and_resting_shield_retention_verified": True,
        "obs_l3_multimir_spliced_fraction": 0.93,
        "obs_l3_multimir_spliced_fraction_se": 0.01,
        "obs_l3_multimir_5prime_seed_fidelity_pct": 92.0,
        "obs_l3_multimir_5prime_seed_fidelity_se": 1.0,
        "obs_l3_multimir_ago2_mir375_retention_pct": 95.0,
        "obs_l3_multimir_ago2_mir375_retention_se": 1.0,
        "obs_l3_multimir_shield_protein_retention_pct": 91.0,
        "obs_l3_multimir_shield_protein_retention_se": 1.5,
        "obs_l3_multimir_min_receptor_knockdown_pct": 78.0,
        "obs_l3_multimir_min_receptor_knockdown_se": 1.5,
        "obs_l3_multimir_matched_gsis_output_retention_pct": 93.0,
        "obs_l3_multimir_matched_gsis_output_retention_se": 1.5,
    }
    _dt_mm_ok = wetlab.evaluate_stage0_decision_tree(_l3_mm_valid)
    assert _dt_mm_ok["resolved_architecture_assignments"]["locus3_multimir_modeb_gate_status"] == "CLEARED_L3_MULTIMIR_MODEB_STAGE0_CHALLENGER_GATE"
    assert _dt_mm_ok["resolved_architecture_assignments"]["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

    assert wetlab.evaluate_stage0_decision_tree({**_l3_mm_valid, "knockin_splice_challenger_arm": "ARM2_L3_V2_INTRON_NT_SCRAMBLE_4398"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_L3_MULTIMIR_FACTORIAL_CONTROL_ARM_NOT_RELEASE_CHALLENGER"
    assert wetlab.evaluate_stage0_decision_tree({**_l3_mm_valid, "l3_multimir_pre_aav6_2x2_transient_kd_falsifier_cleared": None})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_L3_MULTIMIR_PRE_AAV6_2X2_TRANSIENT_FALSIFIER_UNVERIFIED"
    assert wetlab.evaluate_stage0_decision_tree({**_l3_mm_valid, "knockin_splice_challenger_arm": "ARM6_L3_V2_INTRON_DUPLEX_IL1R1_FAS_C341S_TCA_4398"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_L3_MULTIMIR_ARM6_C341S_HLAC_OR_DQ2_DQ8_PEPTIDE_SCREEN_UNVERIFIED"
    assert wetlab.evaluate_stage0_decision_tree({**_l3_mm_valid, "obs_l3_multimir_ago2_mir375_retention_pct": 85.0})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_L3_MULTIMIR_SPLICING_DROSHA_SEED_OR_AGO2_MIR375_FAILED"
    assert wetlab.evaluate_stage0_decision_tree({**_l3_mm_valid, "obs_l3_multimir_matched_gsis_output_retention_pct": 80.0})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_L3_MULTIMIR_RECEPTOR_KD_SHIELD_RETENTION_OR_MATCHED_GSIS_FAILED"

    _cal_mm_1 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP18-1",
            "mcb_clone_id": "CLONE-IIP18",
            "knockin_splice_challenger_arm": "ARM5_L3_V2_INTRON_DUPLEX_IL1R1_FAS_WT_SB9_4398",
            "obs_stage7_yield_ratio": 0.94,
            "obs_stage7_yield_se": 0.02,
        },
        is_dry_run=False,
    )
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-IIP18-2",
                "mcb_clone_id": "CLONE-IIP18",
                "knockin_splice_challenger_arm": "ARM6_L3_V2_INTRON_DUPLEX_IL1R1_FAS_C341S_TCA_4398",
                "obs_cd47_fold": 4.4,
                "obs_cd47_fold_se": 0.2,
            },
            False,
            _cal_mm_1,
        )
        raise AssertionError("Expected ValueError when pooling ARM5 and ARM6")
    except ValueError as _exc_mm:
        assert "knockin_splice_challenger_arm" in str(_exc_mm)
    _r9_s2_bp_md_18 = (ROOT / "ISLET-1-BLUEPRINT.md").read_text(encoding="utf-8")
    assert "IIP-0018" in _r9_s2_bp_md_18 and "DEC-096" in _r9_s2_bp_md_18

    # IIP-0019 / DEC-097: Layer-4 explicit-3D transport audit + Stage-0 S6 / P05 / Design-Z fail-closed qualification gates
    _l4_iip19 = l4["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["explicit_3d_clumping_and_microbead_h2o2_audit_iip_0019"]
    assert _l4_iip19["iip_id"] == "IIP-0019" and _l4_iip19["decision_id"] == "DEC-097"
    _g3dv = _l4_iip19["3d_solver_verification_gate"]
    assert _g3dv["gate_id"] == "Gate_3D_V" and _g3dv["distinct_from_seed_line_gate_0"] is True
    assert _g3dv["planar_1d_o2_rel_diff_pct"] < 0.5 and _g3dv["planar_1d_h2o2_rel_diff_pct"] < 1.0
    assert _g3dv["spherical_120um_closure_reconciliation"]["dirichlet_3d_vs_1d_rel_diff_pct"] < 0.5
    assert abs(_g3dv["spherical_120um_closure_reconciliation"]["layer4_coupled_wigner_seitz_uniform_supply_core_mmhg"] - 4.0762) < 1e-3
    assert abs(_g3dv["spherical_120um_closure_reconciliation"]["layer4_coupled_wigner_seitz_gel_only_lower_bound_core_mmhg"] - 2.1468) < 1e-3
    assert _g3dv["negative_controls"]["all_3_negative_controls_passed"] is True
    _s6_3d = _l4_iip19["rule_s6_3d_lateral_clumping_capsule_and_dose_matching_audit"]
    assert _s6_3d["exact_operator_equivalence_15um_phys_collagen_vs_45um_geleq_max_diff_mmhg"] == 0.0
    assert _s6_3d["physical_collagen_vs_hydrogel_equivalent_standoff_comparison"]["100um_1500ieq_15um_physical_collagen_45um_hydrogel_equiv"]["fv_harmonic_face_mean_mmhg"] > 17.0
    assert _s6_3d["physical_collagen_vs_hydrogel_equivalent_standoff_comparison"]["100um_1500ieq_45um_physical_collagen_135um_hydrogel_equiv"]["fv_harmonic_face_mean_mmhg"] < 13.0
    assert _s6_3d["cluster_80um_dose_matched_1500ieq_vs_pitch_held_768ieq_comparison"]["dose_matched_80um_1500ieq_142p22um_square_152p83um_hex_pitch"]["square_pitch_um"] == 142.22
    assert _s6_3d["joint_plus10pct_100um_knife_edge_discretization_audit"]["classification"] == "NUMERICALLY_KNIFE_EDGE_UNRESOLVED"
    _rsa_ens = _s6_3d["stochastic_3d_rsa_10_tile_160_cluster_ensemble"]
    assert len(_rsa_ens["tiles"]) == 10 and _rsa_ens["fail_shell_0p95_1p06_lt15_count"] == 67 and _rsa_ens["clusters_with_any_voxel_lt0p5mmhg_severe_anoxia_count"] == 0
    _h2o2_blk = _l4_iip19["discrete_microbead_h2o2_and_barrier_film_cross_section_audit"]
    assert _h2o2_blk["tier1_2p5uM_target_vs_tier2_5p0uM_ceiling_and_catalase_budget"]["zero_added_enzyme_split_tier1_2p5uM_status_at_24h_or_50pct_washout"] == "EXCEEDS_2P5UM_TIER1_TARGET_CLEARS_ONLY_5P0UM_TIER2_CEILING"
    assert abs(_h2o2_blk["tier1_2p5uM_target_vs_tier2_5p0uM_ceiling_and_catalase_budget"]["required_post_cure_t0_active_catalase_at_70pct_day7_retention_iu_333cm2"] - 105428.6) < 0.2

    _t3d_valid = {
        **full_ext_payload,
        "s6_cluster_size_arm": "DOSE_MATCHED_80UM_1500IEQ",
        "s6_dose_matched_1500_ieq_verified": True,
        "obs_s6_p10_nearest_neighbor_edge_gap_um": 62.0,
        "obs_s6_p10_nearest_neighbor_edge_gap_se": 1.5,
        "obs_s6_cluster_diam_p90_um": 76.0,
        "obs_s6_cluster_diam_p90_se": 1.0,
        "p05_barrier_slab_geometry_mode": "EXPANDED_390UM_2X125UM_BANDS",
        "p05_full_336h_coupled_o2_h2o2_ca_ph_qualified": True,
        "obs_p05_microbead_volume_d99_um": 45.0,
        "obs_p05_microbead_volume_d99_se": 1.0,
        "obs_p05_catalase_barrier_film_um": 22.0,
        "obs_p05_catalase_barrier_film_se": 0.5,
        "obs_p05_islet_plane_h2o2_max_uM": 1.4,
        "obs_p05_islet_plane_h2o2_max_se": 0.15,
        "obs_p05_day7_catalase_retention_pct": 78.0,
        "obs_p05_day7_catalase_retention_se": 1.5,
        "obs_p05_336h_retained_catalase_iu_333cm2": 78000.0,
        "obs_p05_336h_retained_catalase_iu_333cm2_se": 500.0,
        "design_z_divider_arm": "PLAIN_HYDROGEL_140UM",
        "design_z_empty_lumen_14d_persistence_verified": True,
        "obs_design_z_empty_lumen_5site_min_po2_mmhg": 43.0,
        "obs_design_z_empty_lumen_5site_min_po2_se": 1.0,
        "obs_design_z_vessel_standoff_p90_um": 12.0,
        "obs_design_z_vessel_standoff_p90_se": 0.5,
        "design_z_loaded_closed_state_confirmed": True,
        "obs_design_z_loaded_surface_po2_mmhg": 19.5,
        "obs_design_z_loaded_surface_po2_se": 0.5,
        "obs_design_z_loaded_core_po2_mmhg": 11.2,
        "obs_design_z_loaded_core_po2_se": 0.4,
        "obs_design_z_loaded_perifusion_gsis_si": 3.4,
        "obs_design_z_loaded_perifusion_gsis_se": 0.1,
    }
    _dt_t3d_ok = wetlab.evaluate_stage0_decision_tree(_t3d_valid)
    assert _dt_t3d_ok["resolved_architecture_assignments"]["transport_3d_clumping_h2o2_and_design_z_evaluation"]["gate_pass"] is True
    assert _dt_t3d_ok["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

    assert wetlab.evaluate_stage0_decision_tree({**_t3d_valid, "s6_cluster_size_arm": "DERATED_80UM_768IEQ"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_S6_DOSE_DERATED_80UM_768IEQ_NOT_SINGLE_POUCH_EQUIVALENT_REQUIRES_DOSE_MATCHED_1500IEQ_OR_AREA_EXPANSION"
    assert wetlab.evaluate_stage0_decision_tree({**_t3d_valid, "obs_p05_islet_plane_h2o2_max_uM": 2.815})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_P05_ISLET_PLANE_H2O2_EXCEEDS_2P5UM_TIER1_TARGET_REQUIRES_FULL_73800IU_RETENTION"
    assert wetlab.evaluate_stage0_decision_tree({**_t3d_valid, "obs_p05_islet_plane_h2o2_max_uM": 6.88})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_P05_ISLET_PLANE_H2O2_EXCEEDS_5P0UM_TIER2_ABORT_CEILING"
    assert wetlab.evaluate_stage0_decision_tree({**_t3d_valid, "p05_full_336h_coupled_o2_h2o2_ca_ph_qualified": False})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_P05_FULL_336H_COUPLED_TRANSPORT_QUALIFICATION_UNVERIFIED"
    assert wetlab.evaluate_stage0_decision_tree({**_t3d_valid, "design_z_loaded_closed_state_confirmed": False})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_DESIGN_Z_EMPTY_LUMEN_ONLY_REQUIRES_LOADED_CLOSED_STATE_CONFIRMATION"
    assert wetlab.evaluate_stage0_decision_tree({**_t3d_valid, "obs_design_z_loaded_perifusion_gsis_si": 2.6})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_DESIGN_Z_LOADED_CLOSED_STATE_PO2_OR_READOUT2_GSIS_SI3_FAILED"

    _cal_t3d_1 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP19-1",
            "mcb_clone_id": "CLONE-IIP19",
            "s6_cluster_size_arm": "DXB_100UM_1500IEQ",
            "obs_stage7_yield_ratio": 0.94,
            "obs_stage7_yield_se": 0.02,
        },
        is_dry_run=False,
    )
    try:
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-IIP19-2",
                "mcb_clone_id": "CLONE-IIP19",
                "s6_cluster_size_arm": "DOSE_MATCHED_80UM_1500IEQ",
                "obs_cd47_fold": 4.4,
                "obs_cd47_fold_se": 0.2,
            },
            False,
            _cal_t3d_1,
        )
        raise AssertionError("Expected ValueError when pooling distinct s6_cluster_size_arm conditions")
    except ValueError as _exc_t3d:
        assert "s6_cluster_size_arm" in str(_exc_t3d)
    assert "IIP-0019" in _r9_s2_bp_md_18 and "DEC-097" in _r9_s2_bp_md_18

    # 45. IIP-0020 / DEC-099: Extracellular biomaterial census, Gate P05 P04040 lock, humoral/DPR/hepatocyte cuff gates & OQ-016 reassessment:
    _l2_iip20 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))["extracellular_biomaterial_foreign_peptide_census_iip_0020"]
    assert _l2_iip20["decision_id"] == "DEC-099" and _l2_iip20["pinned_human_swissprot_protein_count"] == 20429
    assert _l2_iip20["audited_extracellular_chains"]["Human_CAT_P04040"]["positional_foreign_8_11mers"] == "0/2074"
    assert _l2_iip20["audited_extracellular_chains"]["Human_CAT_P04040_MetCleaved_2_527"]["positional_foreign_8_11mers"] == "0/2070"
    assert _l2_iip20["audited_extracellular_chains"]["Bovine_CAT_P00432"]["positional_foreign_8_11mers"] == "1049/2074"
    assert _l2_iip20["audited_extracellular_chains"]["CFH_19_448"]["positional_foreign_8_11mers"] == "0/1686"
    _l3_iip20 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))["biomaterial_cuff_and_catalase_structural_audit_iip_0020"]
    assert _l3_iip20["registered_pdb_entries"]["2FHA"]["status"] == "REJECTED_GEOMETRIC_AND_CYSTEINE_MISMATCH"

    _iip20_valid = {
        **_t3d_valid,
        "design_z_divider_arm": None,
        "require_full_iip0020_biomaterial_and_humoral_panel": True,
        "p05_catalase_uniprot_accession": "P04040",
        "p05_recombinant_human_catalase_verified": True,
        "obs_p05_catalase_endotoxin_eu_per_mg": 0.04,
        "obs_p05_catalase_endotoxin_eu_per_mg_se": 0.01,
        "obs_p05_total_implant_catalase_endotoxin_eu": 0.18,
        "obs_p05_total_implant_catalase_endotoxin_eu_se": 0.03,
        "obs_p05_day7_shed_catalase_mass_pct": 12.0,
        "obs_p05_day7_shed_catalase_mass_se": 1.5,
        "cuff_bioconjugation_arm": "DPR_HYDROLYZED_CFH_THBD_H1A_ADIPOQ_HEXAFASL",
        "cuff_serum_challenge_stratum": "HIGH_ANTI_PEG_IGG_GT_500NG_ML",
        "cuff_dpr_ring_hydrolysis_lcms_verified": True,
        "cuff_conjugate_disulfide_and_site_occupancy_lcms_verified": True,
        "cuff_hexafasl_inventory_within_provisional_cap_verified": True,
        "obs_cuff_surface_c3b_deposition_pct_pos_ctrl": 2.1,
        "obs_cuff_surface_c3b_deposition_se": 0.4,
        "obs_cuff_soluble_sc5b9_pct_pos_ctrl": 2.4,
        "obs_cuff_soluble_sc5b9_se": 0.4,
        "obs_cuff_macrophage_fcgammar_cytokine_fold": 1.15,
        "obs_cuff_macrophage_fcgammar_cytokine_se": 0.05,
        "obs_cuff_primed_serum_fasl_ec50_fold_shift": 1.25,
        "obs_cuff_primed_serum_fasl_ec50_fold_se": 0.10,
        "obs_cuff_primed_serum_fasl_emax_apoptosis_pct": 88.0,
        "obs_cuff_primed_serum_fasl_emax_apoptosis_se": 1.5,
        "obs_cuff_primed_serum_thbd_apc_retention_pct": 89.0,
        "obs_cuff_primed_serum_thbd_apc_retention_se": 1.5,
        "obs_cuff_primed_serum_cfh_c3b_cleavage_retention_pct": 90.0,
        "obs_cuff_primed_serum_cfh_c3b_cleavage_retention_se": 1.5,
        "obs_cuff_dpr_ring_hydrolysis_pct": 99.1,
        "obs_cuff_dpr_ring_hydrolysis_se": 0.2,
        "obs_cuff_14d_hsa_shed_bioactive_fasl_pct": 1.8,
        "obs_cuff_14d_hsa_shed_bioactive_fasl_se": 0.3,
        "obs_cuff_h1a_sec_mals_hexamer_purity_pct": 94.5,
        "obs_cuff_h1a_sec_mals_hexamer_purity_se": 0.8,
        "obs_cuff_unreacted_free_thiol_pct": 1.9,
        "obs_cuff_unreacted_free_thiol_se": 0.3,
        "obs_cuff_hepatocyte_72h_caspase37_fold": 1.04,
        "obs_cuff_hepatocyte_72h_caspase37_se": 0.03,
        "obs_cuff_hepatocyte_72h_alt_ldh_fold": 1.05,
        "obs_cuff_hepatocyte_72h_alt_ldh_se": 0.03,
        "obs_cuff_serum_donor_count": 14.0,
        "obs_cuff_serum_donor_count_se": 0.001,
        "obs_cuff_serum_failing_donor_count": 0.0,
        "obs_cuff_serum_failing_donor_count_se": 0.001,
    }
    _clean_iip20_valid = {k: v for k, v in _iip20_valid.items() if not k.startswith("design_z_") and not k.startswith("obs_design_z_")}
    _dt_iip20_ok = wetlab.evaluate_stage0_decision_tree(_clean_iip20_valid)
    assert _dt_iip20_ok["resolved_architecture_assignments"]["biomaterial_antigen_humoral_and_human_cuff_evaluation"]["gate_pass"] is True
    assert _dt_iip20_ok["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    assert wetlab.evaluate_stage0_decision_tree({**_clean_iip20_valid, "p05_catalase_uniprot_accession": "P00432_BOVINE_DISQUALIFIED"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_P05_XENO_OR_UNVERIFIED_CATALASE_SPECIES_REQUIRES_HUMAN_P04040"
    assert wetlab.evaluate_stage0_decision_tree({**full_ext_payload, "require_full_iip0020_biomaterial_and_humoral_panel": True})["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    assert wetlab.evaluate_stage0_decision_tree({**_clean_iip20_valid, "obs_cuff_dpr_ring_hydrolysis_pct": 95.0})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CUFF_THIOSUCCINIMIDE_RING_HYDROLYSIS_LT_98PCT_RETRO_MICHAEL_HAZARD"
    assert wetlab.evaluate_stage0_decision_tree({**_clean_iip20_valid, "obs_cuff_hepatocyte_72h_caspase37_fold": 1.35})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CUFF_PRIMARY_HUMAN_HEPATOCYTE_LEACHATE_OR_DIGEST_TOXICITY_EXCEEDED"

    # 46. IIP-0021 / DEC-100: Two-scale lineage OCR helper, parameterized Design Z / washout / M101 audit, and Stage-0 120-h OCR + 4-arm conditioning gates:
    from scripts import run_layer4_tissue_and_metabolic_twin as l4_iip21
    _l4_iip21 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["lineage_partitioned_metabolic_rest_and_m101_exclusion_audit"]
    assert _l4_iip21["decision_id"] == "DEC-100" and _l4_iip21["baseline_cao2_mass_mg_preserved"] == 791.5
    _ts_res = l4_iip21.run_lineage_partitioned_intrinsic_ocr_regime(
        regime_id="TEST_IIP21_TWO_SCALE",
        areal_loading_ieq_cm2=2500.0,
        endocrine_intrinsic_ocr_multiplier=0.55,
        vascular_intrinsic_ocr_multiplier=1.0,
    )
    assert abs(_ts_res["slab_normalized_ocr_multiplier"] - (0.55 + 0.15) / 1.15) < 1e-6
    assert abs(_ts_res["cluster_endocrine_ocr_multiplier"] - 0.55) < 1e-6
    assert abs(_l4_iip21["m101_closed_slab_standalone_reservoir_bounds"]["bounds"]["HostBed_40_to_15mmHg_n2p5"]["releasable_umol_above_floor"] - 0.25707) < 1e-4

    _iip21_valid = {
        **_clean_iip20_valid,
        "require_full_iip0021_ocr_and_conditioning_panel": True,
        "rest_ocr_coupon_arm": "ENCAPSULATED_KATP_MICRODEPOT_CHALLENGER",
        "rest_120h_longitudinal_simultaneous_coverage_verified": True,
        "rest_120h_viable_lineage_yield_and_during_exposure_function_verified": True,
        "obs_rest_120h_graft_ocr_ratio": 0.72,
        "obs_rest_120h_graft_ocr_ratio_se": 0.02,
        "obs_rest_depot_release_half_life_h": 84.0,
        "obs_rest_depot_release_half_life_se_h": 2.0,
        "obs_rest_post_washout_gsis_si": 3.5,
        "obs_rest_post_washout_gsis_se": 0.12,
        "obs_rest_post_washout_gsis_output_retention_pct": 92.0,
        "obs_rest_post_washout_gsis_output_retention_se": 1.5,
        "obs_rest_post_washout_alpha_glucagon_fold": 2.5,
        "obs_rest_post_washout_alpha_glucagon_se": 0.10,
        "conditioning_4arm_arm": "FERROSTATIN1_2UM",
        "conditioning_washout_analyte_scope": "SUPERNATANT_AND_CELL_PELLET_PARENT_PLUS_METABOLITES",
        "conditioning_pellet_loq_and_extraction_recovery_verified": True,
        "conditioning_orthogonal_pathway_controls_verified": True,
        "obs_cond_sc_beta_viable_yield_pct": 91.0,
        "obs_cond_sc_beta_viable_yield_se": 1.2,
        "obs_cond_iec_viable_yield_pct": 90.0,
        "obs_cond_iec_viable_yield_se": 1.2,
        "obs_cond_final_wash_supernatant_nM": 2.0,
        "obs_cond_final_wash_supernatant_se_nM": 0.4,
        "obs_cond_pellet_residual_over_loq_ratio": 0.45,
        "obs_cond_pellet_residual_over_loq_se": 0.05,
        "obs_cond_pellet_extraction_recovery_pct": 88.0,
        "obs_cond_pellet_extraction_recovery_se": 1.5,
        "obs_cond_endocrine_nqo1_activity_pct": 96.0,
        "obs_cond_endocrine_nqo1_activity_se": 1.2,
        "obs_cond_vascular_nqo1_activity_pct": 95.0,
        "obs_cond_vascular_nqo1_activity_se": 1.2,
        "obs_cond_post_washout_gsis_si": 3.4,
        "obs_cond_post_washout_gsis_se": 0.12,
    }
    _dt_iip21_ok = wetlab.evaluate_stage0_decision_tree(_iip21_valid)
    assert _dt_iip21_ok["resolved_architecture_assignments"]["metabolic_rest_ocr_and_conditioning_evaluation"]["gate_pass"] is True
    assert _dt_iip21_ok["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    assert wetlab.evaluate_stage0_decision_tree({**_iip21_valid, "rest_ocr_coupon_arm": "UNENCAPSULATED_FREE_DRUG_REJECTED"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_REST_OCR_UNENCAPSULATED_FREE_DRUG_RAPID_SLAB_WASHOUT_BARRED"
    assert wetlab.evaluate_stage0_decision_tree({**_iip21_valid, "conditioning_4arm_arm": "NEC1_100UM_ANTIFERROPTOTIC_COMPARATOR"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CONDITIONING_BENCHMARK_OR_NEC1_ANTIFERROPTOTIC_COMPARATOR_NOT_RELEASE_ARM"
    assert wetlab.evaluate_stage0_decision_tree({**_iip21_valid, "obs_cond_endocrine_nqo1_activity_pct": 82.0})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_CONDITIONING_RESIDUAL_DRUG_OR_NQO1_INHIBITION_DETECTED"

    # 47. IIP-0022 / DEC-101: Fail-closed Track-2 B'-SCT vs HLA-G1 split, upper95 KIR routing, off-grid Graveyard guardrails, Layer-3 PVR/HLA-G1 audit, and 2^12 common-shock generator:
    from scripts import rank_islet1_stacks as rk_iip22
    _tree_dec098 = rk_iip22.compute_common_shock_2pow12_portfolio_tree(
        [0.88, 0.78, 0.56, 0.68],
        [0.88, 0.79, 0.64, 0.70],
        [0.88, 0.81, 0.70, 0.68],
        [0.88, 0.87, 0.88, 0.85],
    )
    assert _tree_dec098["total_enumerated_states"] == 4096 and _tree_dec098["positive_mass_states"] == 1024
    assert [round(x, 6) for x in _tree_dec098["surviving_portfolio_conditionals_exact"]] == [0.88, 0.869429, 0.847262, 0.794158]
    assert round(_tree_dec098["locked_3slot_2pow12_exact"], 6) == 0.514804
    assert round(_tree_dec098["legacy_modular_3slot_cell_x_pouch_exact"], 6) == 0.547113
    assert round(_tree_dec098["modular_2route_shared_ax_pouch_exact"], 6) == 0.531555
    assert round(_tree_dec098["locked_2route_shared_ax_pouch_exact"], 6) == 0.497144

    _gy_iip22 = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    assert len(_gy_iip22["graveyard_entries"]) == 10
    assert _gy_iip22["off_grid_molecular_dead_end_guardrails"]["decision_id"] == "DEC-101"

    _l3_iip22 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))["pvr_tigit_cd226_cd96_and_hlag1_structural_audit_iip_0022"]
    assert _l3_iip22["authenticated_rcsb_mmcif_records"]["3UDW"]["pvr_ala67_contacts_le_4p5A"] == 0

    # Passing mean (7.0%) with failing upper95 (8.645%) on KIR+NKG2A- NK lysis must route to KIR hold rather than misrouting to S2 Arm 1X:
    _kir_u95_fail = wetlab.evaluate_stage0_decision_tree({
        "obs_hlae_hour0_molecules": 12000.0,
        "obs_hlae_hour0_molecules_se": 500.0,
        "obs_kir_nkg2a_neg_nk_lysis_pct": 7.0,
        "obs_kir_nkg2a_neg_nk_lysis_se": 1.0,
    })
    assert _kir_u95_fail["resolved_architecture_assignments"]["locus1_construct"].startswith("HOLD_KIR_NKG2A_NEG_NK")

    # HLA-G1 route remains HELD pending product registration and never receives B'-SCT locus2_construct:
    _hlag1_attempt = wetlab.evaluate_stage0_decision_tree({
        **_iip21_valid,
        "track2_rescue_route": "TRACK2_HLAG1_LILRB1_CHALLENGER_HELD",
        "track2_hlag1_lilrb1_rescue_qualified": True,
        "obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_pct": 4.0,
        "obs_track2_hlag1_lilrb1_neg_kir_nk_lysis_se": 0.8,
        "obs_track2_hlag1_donor_count": 6.0,
        "obs_track2_hlag1_donor_count_se": 1e-12,
        "obs_track2_hlag1_failing_donor_count": 0.0,
        "obs_track2_hlag1_failing_donor_count_se": 1e-12,
    })
    assert _hlag1_attempt["resolved_architecture_assignments"]["locus2_construct"] == "HOLD_TRACK2_HLAG1_CHALLENGER_PENDING_PRODUCT_REGISTRATION"
    assert _hlag1_attempt["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # 48. IIP-0023 / DEC-102: Design Y wildcard comparator FASTA, Layer-2 CGI/TadCBEd/adenosine manifest, and fail-closed Y-01..Y-04 contract (with zero-event Poisson/SE floor check & unchanged A/B/X when absent):
    assert "design_y_wildcard_status" not in _dt_iip21_ok["resolved_architecture_assignments"]
    _l2_iip23 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))["stage0_design_y_wildcard_manifest_iip_0023"]
    assert _l2_iip23["fasta_artifact"]["file_sha256"] == "c7f54b824a73cf9d40247d4fd251cb5c28a8eca4a1530b9661284200bb4cbabc"
    assert [c["length_bp"] for c in _l2_iip23["fasta_artifact"]["constructs"]] == [11066, 11859, 12555]

    _dy_valid = {
        **_iip21_valid,
        "obs_kir_nkg2a_neg_nk_lysis_pct": 3.5,
        "obs_kir_nkg2a_neg_nk_lysis_se": 0.6,
        "design_y_integrase_modality": "EEPASSIGE_B2M_ATTP_GENOME_ATTB_MINICIRCLE",
        "design_y_silencing_modality": "TADCBED_CIITA_W334X_CXCL10_M1I",
        "design_y_payload_arm": "DESIGN_Y_CORE_CD200_IEC_11066BP",
        "design_y_attlr_outside_orf_5p3p_race_verified": True,
        "design_y_zero_2a_scar_proteomics_verified": True,
        "design_y_dpcpx_adora1_antagonist_control_run": True,
        "design_y_anti_cd200_blockade_reversal_verified": True,
        "design_y_y04_full_coculture_and_clonogenic_confirmed": True,
        "obs_design_y01_max_lineage_escape_pct": 0.25,
        "obs_design_y01_max_lineage_escape_se": 0.05,
        "obs_design_y01_lineage_cells_scored": 15000.0,
        "obs_design_y01_lineage_cells_scored_se": 10.0,
        "obs_design_y01_cbe_offtarget_pct": 0.02,
        "obs_design_y01_cbe_offtarget_se": 0.005,
        "obs_design_y02_pre_enrichment_ki_pct": 5.2,
        "obs_design_y02_pre_enrichment_ki_se": 0.3,
        "obs_design_y02_clonal_ddpcr_cn": 2.00,
        "obs_design_y02_clonal_ddpcr_cn_se": 0.02,
        "obs_design_y02_wt_b2m_residual_allele_frac": 0.0,
        "obs_design_y02_wt_b2m_residual_allele_se": 1e-5,
        "obs_design_y02_backbone_kanr_vaf": 0.0,
        "obs_design_y02_backbone_kanr_vaf_se": 1e-6,
        "obs_design_y02_wgs_pseudosite_k": 0.0,
        "obs_design_y02_wgs_pseudosite_k_se": 1e-12,
        "obs_design_y03_iec_beta_selectivity_ratio": 85.0,
        "obs_design_y03_iec_beta_selectivity_se": 4.0,
        "obs_design_y03_beta_mir375_retention_pct": 97.0,
        "obs_design_y03_beta_mir375_retention_se": 1.0,
        "obs_design_y03_coculture_gsis_si": 3.4,
        "obs_design_y03_coculture_gsis_se": 0.12,
        "obs_design_y03_gsis_retention_pct": 96.0,
        "obs_design_y03_gsis_retention_se": 1.2,
        "obs_design_y03_supernatant_adenosine_nM": 8.0,
        "obs_design_y03_supernatant_adenosine_se_nM": 1.0,
        "obs_design_y03_dpcpx_gsis_delta_pct": 3.0,
        "obs_design_y03_dpcpx_gsis_delta_se": 0.8,
        "obs_design_y04_multicistron_retention_pct": 97.5,
        "obs_design_y04_multicistron_retention_se": 0.5,
        "obs_design_y04_nhs20_c5b9_mac_pct": 2.1,
        "obs_design_y04_nhs20_c5b9_mac_se": 0.3,
        "obs_design_y04_indirect_apc_cd4_si": 1.15,
        "obs_design_y04_indirect_apc_cd4_se": 0.06,
        "obs_design_y04_ap1903_clonogenic_survival_frac": 0.0,
        "obs_design_y04_ap1903_clonogenic_survival_se": 1e-12,
        "obs_design_y_sampling_age_days": 62.0,
        "obs_design_y_sampling_age_days_se": 1e-12,
        "obs_design_y02_ont_intact_homolog_pct": 99.6,
        "obs_design_y02_ont_intact_homolog_se": 0.1,
        "obs_design_y04_all4_lineages_scored_count": 4.0,
        "obs_design_y04_all4_lineages_scored_se": 1e-12,
        "obs_design_y04_ap1903_survivor_count_k": 0.0,
        "obs_design_y04_ap1903_survivor_count_k_se": 1e-12,
        "obs_design_y04_ap1903_effective_trials": 3.05e8,
        "obs_design_y04_ap1903_effective_trials_se": 1.0,
    }
    _dy_res_ok = wetlab.evaluate_stage0_decision_tree(_dy_valid)
    assert _dy_res_ok["resolved_architecture_assignments"]["design_y_wildcard_evaluation"]["gate_pass"] is True
    assert wetlab.evaluate_stage0_decision_tree({**_dy_valid, "obs_design_y_sampling_age_days": 49.0})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_DESIGN_Y01_60D_4LINEAGE_SILENCE_OR_TADCBED_HAPLOTYPE_FAILED"
    assert wetlab.evaluate_stage0_decision_tree({**_dy_valid, "design_y_integrase_modality": "LEGACY_ATTB_GENOME_ATTP_PLASMID_REJECTED"})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_DESIGN_Y02_LEGACY_ATTB_IN_GENOME_POLARITY_OT4_PSEUDOSITE_RISK_REJECTED"
    assert wetlab.evaluate_stage0_decision_tree({**_dy_valid, "obs_design_y04_ap1903_effective_trials": 1.0e7})["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_DESIGN_Y04_MULTICISTRON_COMPLEMENT_KIR_OR_3E8_CLONOGENIC_GATE_FAILED"

    # 49. IIP-0024 / DEC-103: Atomic provenance scopes, cross-lot & mean-only SE isolation, integer (k, n) & k==0 vetoes,
    # non-executable comparator holds (ISLET1_Z_X_MOL_140UM_ZERO_CAO2 / B2_BASE / B'-SCT), Student-t & 6+6 Welch TOST,
    # conjunctive DMR/XCI/MDM4/mtDNA guardrails, and Layers 1, 3, 4 & 5 hardening:
    from scripts import product_contract as pc_iip24
    from scripts import wetlab_contract as wc_iip24
    from scripts import run_layer3_structural_binder_plane as l3_iip24
    from scripts import run_layer5_personalized_patient_compiler as l5_iip24

    # (a) Product catalog registration & whitespace rejection:
    _prod_z = pc_iip24.get_product_definition("ISLET1_Z_X_MOL_140UM_ZERO_CAO2", root=ROOT)
    assert _prod_z["execution_supported"] is False and _prod_z["physical_route_id"] == "Z_OMENTAL_140UM_ZERO_CAO2"
    assert _prod_z["cao2_and_p04040_catalase_exempt"] is True
    with pytest.raises(ValueError):
        pc_iip24.get_product_definition(" ISLET1_A4_ARM1B_EMB_A_STAR_4143 ", root=ROOT)

    # (b) Deep-immutable AtomicObservation & typed 4/5/7-key scope envelopes:
    _c_sha = "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"
    _mcb_sc = wc_iip24.MCBScopeKey("ISLET1_A4_ARM1B_EMB_A_STAR_4143", _c_sha, "SEED-OO-01", "MCB-01")
    _lot_sc = wc_iip24.LotScopeKey("ISLET1_A4_ARM1B_EMB_A_STAR_4143", _c_sha, "SEED-OO-01", "MCB-01", "LOT-01")
    _dev_sc = wc_iip24.AssembledDeviceLotScopeKey(
        "ISLET1_A4_ARM1B_EMB_A_STAR_4143", _c_sha, "SEED-OO-01", "MCB-01", "LOT-01", "DEV-01", "A_STAR_350UM_CAO2_P04040"
    )
    assert len(_mcb_sc.as_tuple()) == 4 and len(_lot_sc.as_tuple()) == 5 and len(_dev_sc.as_tuple()) == 7
    _atom = wc_iip24.AtomicObservation(
        endpoint_id="stage0_sc_beta_yield_ratio",
        scope_type="LOT_SCOPE",
        scope_key=_lot_sc,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=0.97,
        se_or_n=0.01,
        sample_n=6,
        controls={"batch_id": "B-ATOM-1"},
    )
    with pytest.raises(TypeError):
        _atom.controls["batch_id"] = "TAMPERED"  # type: ignore[index]
    with pytest.raises(ValueError):
        wc_iip24.AtomicObservation(
            endpoint_id="stage0_sc_beta_yield_ratio",
            scope_type="LOT_SCOPE",
            scope_key=_lot_sc,
            kind="CONTINUOUS_STUDENT_T",
            value_or_k=0.97,
            se_or_n=0.01,
            sample_n=6,
            canonical_projection={"spoof": True},
        )

    # (c) Cross-lot evidence borrowing block & mean-only SE inheritance block in ingest_wetlab_results:
    _b_lot1 = wetlab.run_calibration(
        {
            **_v3_valid_payload,
            "batch_id": "BATCH-IIP24-LOT1",
            "mcb_clone_id": "MCB-IIP24-01",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "current_lot_id": "LOT-SP6B-01",
        },
        is_dry_run=False,
    )
    assert _b_lot1["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    _b_lot_untested = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP24-LOT-UNTESTED",
            "mcb_clone_id": "MCB-IIP24-01",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "s10_construct_id": "S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143",
            "current_lot_id": "LOT-UNTESTED",
            "obs_stage7_yield_ratio": 0.98,
            "obs_stage7_yield_se": 0.01,
        },
        False,
        _b_lot1,
    )
    assert _b_lot_untested["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    assert _b_lot_untested["stage0_decision_tree_evaluation"]["mandatory_s9_safety_panel_complete"] is False

    _b_mean_only_se = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP24-MEAN-ONLY",
            "mcb_clone_id": "MCB-IIP24-01",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "s10_construct_id": "S10_CHAL1_SP6B_VMAPRTVLL_ASTAR_4143",
            "current_lot_id": "LOT-SP6B-01",
            "s10_qualified_lot_id": "LOT-SP6B-01",
            "obs_hlae_hour0_molecules": 27500.0,
        },
        False,
        _b_lot1,
    )
    assert _b_mean_only_se["stage0_decision_tree_evaluation"]["stage0_se_provenance"]["obs_hlae_hour0_molecules"] == "DEFAULT_SE_IMPUTED_NOT_EXPLICIT"
    assert _b_mean_only_se["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"

    # (d) Execution-disabled product variants & bidirectional route-mismatch holds:
    _dz_fields = {
        "design_z_divider_arm": "PLAIN_HYDROGEL_140UM",
        "design_z_empty_lumen_14d_persistence_verified": True,
        "design_z_loaded_closed_state_confirmed": True,
        "obs_design_z_empty_lumen_5site_min_po2_mmhg": 44.0,
        "obs_design_z_empty_lumen_5site_min_po2_se": 0.5,
        "obs_design_z_vessel_standoff_p90_um": 11.0,
        "obs_design_z_vessel_standoff_p90_se": 0.5,
        "obs_design_z_loaded_surface_po2_mmhg": 18.5,
        "obs_design_z_loaded_surface_po2_se": 0.5,
        "obs_design_z_loaded_core_po2_mmhg": 8.2,
        "obs_design_z_loaded_core_po2_se": 0.3,
        "obs_design_z_loaded_perifusion_gsis_si": 3.5,
        "obs_design_z_loaded_perifusion_gsis_se": 0.1,
    }
    _dz_only_base = {
        k: v
        for k, v in _clean_iip20_valid.items()
        if not k.startswith("p05_")
        and not k.startswith("obs_p05_")
        and k != "require_full_iip0020_biomaterial_and_humoral_panel"
    }
    _dt_z_ok = wetlab.evaluate_stage0_decision_tree({**_dz_only_base, **_dz_fields, "product_variant_id": "ISLET1_Z_X_MOL_140UM_ZERO_CAO2"})
    assert _dt_z_ok["lot_release_readiness_tier"] == "WILDCARD_COUPON_QUALIFIED_NON_EXECUTABLE"
    assert wetlab.evaluate_stage0_decision_tree({**_dz_only_base, "product_variant_id": "ISLET1_Z_X_MOL_140UM_ZERO_CAO2"})["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    assert wetlab.evaluate_stage0_decision_tree({**_clean_iip20_valid, **_dz_fields, "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143"})["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    assert wetlab.evaluate_stage0_decision_tree({**full_ext_payload, "product_variant_id": "ISLET1_B2_5LOCUS_BASE"})["lot_release_readiness_tier"] == "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"
    _dt_bsct = wetlab.evaluate_stage0_decision_tree(
        {
            **full_ext_payload,
            "product_variant_id": "ISLET1_B_SCT_C1_CIITA",
            "obs_gate0b_hlac_sct_molecules_per_cell": 25000.0,
            "obs_gate0b_hlac_sct_molecules_se": 1000.0,
            "b_prime_sct_recipient_cd8_mlc_cleared": True,
        }
    )
    assert _dt_bsct["lot_release_readiness_tier"] == "RESEARCH_COMPARATOR_UNSUPPORTED_VARIANT_HOLD"
    assert _dt_bsct["resolved_architecture_assignments"]["locus1_construct"] != "ISLET1_A4_ARM1B"

    # (e) Strict integer (k, n) counts, sensitivity U95 = (1 - 0.05**(1/n))/eta, k==0 Gate 0V veto, & monotonic S9 latching:
    assert wetlab.evaluate_stage0_decision_tree({**full_ext_payload, "obs_20q11_bcl2l1_single_cell_mosaic_cells_scored": 298.9})["lot_release_readiness_tier"] in ("HOLD_FAILED_OR_INCOMPLETE_GATE", "HARD_SAFETY_STOP_REJECTED")
    _no_20q_denom = {k: v for k, v in full_ext_payload.items() if k != "obs_20q11_bcl2l1_single_cell_mosaic_cells_scored"}
    assert wetlab.evaluate_stage0_decision_tree({**_no_20q_denom, "obs_20q11_bcl2l1_single_cell_mosaic_fraction": 0.005, "obs_20q11_bcl2l1_single_cell_mosaic_se": 0.001})["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    assert wc_iip24.evaluate_rare_event_count_iip0024(gate_id="G", k=0, n=298, max_u95_fraction=0.01)["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
    assert wc_iip24.evaluate_rare_event_count_iip0024(gate_id="G", k=0, n=299, max_u95_fraction=0.01)["status"] == "PASS"
    _re_920 = wc_iip24.evaluate_rare_event_count_iip0024(gate_id="G", k=0, n=920, max_u95_fraction=0.01, eta_spike_recovery_l95=0.325)
    assert _re_920["status"].startswith("HOLD") and abs(_re_920["u95_binomial"] - 0.01000288) < 1e-6
    assert wc_iip24.evaluate_rare_event_count_iip0024(
        gate_id="G", k=0, n=1200, max_u95_fraction=0.01, clustered_sampling=True, validated_clustered_estimator_id="UNREGISTERED_DEFF_SHORTCUT"
    )["status"].startswith("HOLD")
    _dt_0v_veto = wetlab.evaluate_stage0_decision_tree(
        {
            **_clean_iip20_valid,
            "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct": 0.01,
            "obs_gate0v_day28_persistent_vp1_ie1_cells_scored": 10000,
        }
    )
    assert _dt_0v_veto["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"

    _s9_stop_batch = wetlab.run_calibration(
        {
            **_clean_iip20_valid,
            "batch_id": "BATCH-IIP24-0V-STOP",
            "mcb_clone_id": "MCB-IIP24-LATCH",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "current_lot_id": "LOT-LATCH-01",
            "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct": 0.01,
            "obs_gate0v_day28_persistent_vp1_ie1_cells_scored": 10000,
        },
        is_dry_run=False,
    )
    assert _s9_stop_batch["stage0_decision_tree_evaluation"]["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"
    _s9_pass1 = wetlab.run_calibration(
        {
            **_clean_iip20_valid,
            "batch_id": "BATCH-IIP24-0V-PASS1",
            "mcb_clone_id": "MCB-IIP24-LATCH",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "current_lot_id": "LOT-LATCH-01",
            "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct": 0.0,
            "obs_gate0v_day28_persistent_vp1_ie1_cells_scored": 10000,
        },
        False,
        _s9_stop_batch,
    )
    assert _s9_pass1["stage0_decision_tree_evaluation"]["safety_kill_s9_triggered"] is True
    with pytest.raises(ValueError):
        wetlab.run_calibration(
            {
                **_clean_iip20_valid,
                "batch_id": "BATCH-IIP24-0V-BADCAT",
                "mcb_clone_id": "MCB-IIP24-LATCH",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                "current_lot_id": "LOT-LATCH-01",
                "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct": 0.0,
                "obs_gate0v_day28_persistent_vp1_ie1_cells_scored": 10000,
                "reset_prior_safety_stop_after_investigation": True,
                "reset_prior_safety_stop_justification": "ROOT_CAUSE_WASHOUT_VERIFIED_WRONG_CATEGORY",
                "reset_prior_safety_stop_approver_id": "QP-INDEPENDENT-01",
                "reset_prior_safety_stop_target_category": "DUPLEX_DRIVER_VAF",
            },
            False,
            _s9_pass1,
        )

    # (f) Exact one-sided Student-t (df=5, n=6) vs normal z, and 6+6 vs 6+3 Welch 90% TOST:
    assert abs(wc_iip24.student_t_ppf_one_sided_95(5) - 2.01504837) < 1e-6
    _dt_t_borderline = wetlab.evaluate_stage0_decision_tree(
        {
            **_clean_iip20_valid,
            "obs_hlae_hour0_molecules": 15850.0,
            "obs_hlae_hour0_molecules_se": 450.0,
            "obs_stage0_independent_lot_replicates_n": 6.0,
            "obs_stage0_independent_lot_replicates_n_se": 1e-12,
        }
    )
    assert _dt_t_borderline["lot_release_readiness_tier"] != "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    _tost_res = wc_iip24.evaluate_donor_blocked_contrast_and_welch_tost_iip0024(n_pts=80)
    assert _tost_res["welch_90pct_tost_power_on_0p75_1p3333_pct"]["lots_6_plus_3"]["drift_plus15pct"] < 50.0
    assert _tost_res["welch_90pct_tost_power_on_0p75_1p3333_pct"]["lots_6_plus_6"]["drift_0pct"] >= 96.0
    assert _tost_res["welch_90pct_tost_power_on_0p75_1p3333_pct"]["lots_6_plus_6"]["boundary_type1_sd_0p06_vs_0p15_pct"] <= 5.2

    # (g) Conjunctive DMR, 46,XX XCI (278/300 vs 279/300), MDM4 stop, and >=1,000x mtDNA validators:
    assert wc_iip24.evaluate_conjunctive_imprinting_dmr_iip0024(
        endpoint_scope="LOT_SCOPE",
        dmr_id="MEG3_DMR",
        panel_version_id="MEG3_DMR_SOMATIC_STAGE7_V1",
        per_cpg_callable_verified=True,
        assay_controls_valid=True,
        bulk_beta=0.18,
        phased=False,
        mosaic_discordant_epiallele_frac=0.02,
    )["status"].startswith("HOLD")
    assert wc_iip24.evaluate_conjunctive_xci_46xx_iip0024(
        endpoint_scope="MCB_SCOPE",
        donor_sex_karyotype="46,XX",
        assay_controls_valid=True,
        xci_promoter_panel_id="XCI_SUBJECT_PROMOTER_CPG_PANEL_V1",
        xist_pos_k=278,
        xist_nuclei_n=300,
        xci_subject_promoter_cpg_retention=0.95,
        allele_specific_x_rna_informative_genes=12,
        allele_specific_x_rna_monoallelic_fraction=0.96,
    )["status"].startswith("FAIL")
    assert wc_iip24.evaluate_conjunctive_xci_46xx_iip0024(
        endpoint_scope="MCB_SCOPE",
        donor_sex_karyotype="46,XX",
        assay_controls_valid=True,
        xci_promoter_panel_id="XCI_SUBJECT_PROMOTER_CPG_PANEL_V1",
        xist_pos_k=279,
        xist_nuclei_n=300,
        xci_subject_promoter_cpg_retention=0.95,
        allele_specific_x_rna_informative_genes=12,
        allele_specific_x_rna_monoallelic_fraction=0.96,
    )["status"] == "PASS"
    assert wc_iip24.evaluate_scoped_mtdna_gate_iip0024(
        endpoint_scope="MCB_SCOPE",
        min_unique_chrm_molecule_depth=1500,
        numt_filter_and_snv_del_spikes_validated=True,
        variants=[],
        deletion_vaf_max=None,
        seahorse_coupled_respiration_pass=True,
    )["status"].startswith("HOLD")
    assert wetlab.evaluate_stage0_decision_tree(
        {**_clean_iip20_valid, "obs_1q32_mdm4_cnv_ratio": 2.35, "obs_1q32_mdm4_cnv_se": 0.02}
    )["lot_release_readiness_tier"] == "HARD_SAFETY_STOP_REJECTED"

    # (h) Shared Layers 1, 3, 4 & 5 hardening checks:
    with pytest.raises(ValueError):
        rk_iip22.compute_common_shock_2pow12_portfolio_tree(
            [0.88, 0.87 + 0.5e-12, 0.56, 0.68], [0.88, 0.79, 0.64, 0.70], [0.88, 0.81, 0.70, 0.68], [0.88, 0.87, 0.88, 0.85]
        )
    with pytest.raises(ValueError):
        rk_iip22.compute_common_shock_2pow12_portfolio_tree(
            [True, 0.78, 0.56, 0.68], [0.88, 0.79, 0.64, 0.70], [0.88, 0.81, 0.70, 0.68], [0.88, 0.87, 0.88, 0.85]  # type: ignore[list-item]
        )
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        _rk_json_24 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        assert _rk_json_24["odds_composition_and_profile_scope_audit_iip_0003"]["multiplied_conditional_totals"]["iip_0024_donor_bank_h1_sensitivity_audit"]["decision_id"] == "DEC-103"

    _l3_json_24 = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    _l3_mf_24 = _l3_json_24["product_function_manifest_resolution_iip_0024"]
    assert _l3_mf_24["staged_design_b2_c1_c5_pair_sha256_reference"] == "82c351af082e0f97361e715c36e2ab5872b8d097efa1c9f59dfda53ea60496ea"
    assert _l3_mf_24["function_counts_by_route"]["ISLET1_A4_ARM1B"] == 5
    assert _l3_mf_24["function_counts_by_route"]["ISLET1_B_SCT_C1_CIITA"] == 6
    _l3_ic_24 = _l3_json_24["icasp9_1nw9_1jxq_insertion_code_audit_iip_0024"]["policies"]
    assert _l3_ic_24["default_exclude_insertion_codes"]["xiap_bir3_heavy_atoms_lt_2p2A"] == 565
    assert _l3_ic_24["default_exclude_insertion_codes"]["xiap_bir3_heavy_atoms_total"] == 746
    assert _l3_ic_24["retain_insertion_codes"]["xiap_bir3_heavy_atoms_lt_2p2A"] == 562
    assert _l3_ic_24["retain_insertion_codes"]["xiap_bir3_heavy_atoms_total"] == 746
    assert len(l3_iip24.derive_required_product_functions_iip_0024(_prod_z)) == 5

    # Exact linear threshold-crossing elapsed TBR catching <4.0% sample vs >4.0% elapsed failure:
    _ts_grid = [i * (1440.0 / 600.0) for i in range(601)]
    _gl_grid = [3.80 if 10 <= i < 34 else (3.91 if i in (9, 34) else 5.5) for i in range(601)]
    _el_m = l4_iip21.compute_exact_elapsed_glycemic_metrics_iip_0024(_ts_grid, _gl_grid)
    assert _el_m["sample_fraction_lt_3p9_pct"] < 4.0 and _el_m["tbr_lt_3p9_pct"] > 4.0 and _el_m["passes_layer4_model_gate"] is False
    with pytest.raises(ValueError):
        l4_iip21.compute_exact_elapsed_glycemic_metrics_iip_0024([0.0, 20.0] + [20.0 + i * 5.0 for i in range(1, 290)], [5.5] * 291)

    # Layer 5 explicit-false B'-SCT component overriding summary True while Design A ignores B'-SCT flags:
    _pt_bsct = {
        "patient_id": "PT-IIP24-BSCT",
        "ancestry": "EUR",
        "nkg2a_neg_kir_educated_nk_fraction": 0.55,
        "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
        "ifn1_ish_genetic_risk_percentile": 40.0,
        "sirpa_genotype": "v1/v1",
        "recipient_kir_ligand_group": "C1/C1",
        "hlac_4digit_relation": "KIR_MATCH_FOREIGN_4DIGIT_ALLO",
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared": True,
        "b_prime_sct_recipient_cd8_mlc_negative": False,
        "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative": True,
    }
    _bp_bsct = l5_iip24.compile_personalized_patient_blueprint(_pt_bsct, track2_subbank_mode="B_PRIME_SCT")
    assert _bp_bsct["routing_holds_cleared"] is False
    _eval_a = l5_iip24.evaluate_product_scoped_routing_and_stealth_axes_iip_0024(
        {**_pt_bsct, "exact_hlae_sct_serum_crossmatch_negative": True, "exact_product_non_hla_and_sct_flow_crossmatch_negative": True},
        target_product_def=pc_iip24.get_product_definition("ISLET1_A4_ARM1B", root=ROOT),
    )
    _eval_bsct = l5_iip24.evaluate_product_scoped_routing_and_stealth_axes_iip_0024(
        {**_pt_bsct, "exact_hlae_sct_serum_crossmatch_negative": True, "exact_product_non_hla_and_sct_flow_crossmatch_negative": True},
        target_product_def=pc_iip24.get_product_definition("ISLET1_B_SCT_C1_CIITA", root=ROOT),
    )
    assert _eval_a["scheduler_routing_eligible"] is True and _eval_bsct["scheduler_routing_eligible"] is False
    _eval_unk = l5_iip24.evaluate_product_scoped_routing_and_stealth_axes_iip_0024(
        {"sirpa_genotype": None},
        target_product_def=pc_iip24.get_product_definition("ISLET1_A4_ARM1B", root=ROOT),
    )
    assert _eval_unk["independent_stealth_axes"]["cd47_sirpa_axis"] == "UNKNOWN"

    # =========================================================================
    # 50. Hurdle-4 Spatial-Minimum Oxygen Gates, Product-Bound Two-Timepoint
    #     Human P04040 Catalase Envelope Hold, Unified Student-t Retrieval,
    #     Dual-Assay Shedding & Contrasts B1-B3 (IIP-0025; DEC-104)
    # =========================================================================
    from scripts import wetlab_contract as wc_mod
    _l4_full_25 = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    _h4_25 = _l4_full_25["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["hurdle4_spatial_oxygen_catalase_and_route_challengers_iip_0025"]
    assert _h4_25["contribution_id"] == "IIP-0025" and _h4_25["decision_id"] == "DEC-104"
    assert _h4_25["graveyard_card_id"] == "GY-020"
    assert "PMID:21788575" in _h4_25["primary_citations"] and "PMID:21788581" not in _h4_25["primary_citations"]

    _uc_25 = _h4_25["ordered_unit_cell_multimesh_dual_band_gate"]
    assert _uc_25["nominal_100um_350um_no_film"]["passes_unrounded_multimesh_dual_band_iip_0025"] is False
    assert _uc_25["nominal_100um_350um_no_film"]["fv_face_min_by_nx_40_60_80_mmhg"] == [12.7988, 12.6776, 13.0542]
    assert _uc_25["layout_a_350um_barrier_100um"]["passes_unrounded_multimesh_dual_band_iip_0025"] is False
    assert _uc_25["layout_a_350um_barrier_100um"]["fv_face_min_by_nx_40_60_80_mmhg"] == [12.6447, 12.8978, 12.8766]
    assert _uc_25["nominal_80um_350um_no_film"]["passes_unrounded_multimesh_dual_band_iip_0025"] is True
    assert _uc_25["nominal_80um_350um_no_film"]["fv_face_min_by_nx_40_60_80_mmhg"] == [16.1611, 16.1164, 16.3162]
    assert _uc_25["layout_a_350um_barrier_80um"]["passes_unrounded_multimesh_dual_band_iip_0025"] is True
    assert _uc_25["layout_a_350um_barrier_80um"]["fv_face_min_by_nx_40_60_80_mmhg"] == [15.7101, 15.9496, 15.8598]
    assert _uc_25["plus10pct_joint_stress_80um_layout_a_350um_barrier"]["passes_unrounded_multimesh_dual_band_iip_0025"] is False
    assert _uc_25["nx60_nz105_boundary_plane_symmetry_diagnostics"]["lamina_boundary_asymmetry_voxels"] == 7200
    assert _uc_25["wigner_seitz_unrounded_0p497_boundary_check"]["passes_unrounded_0p5mmhg_gate"] is False
    assert _uc_25["contrast_b1_paired_effect_390um_primary_nx40_mmhg"] == 2.6781
    assert _uc_25["contrast_b1_paired_effect_350um_sensitivity_nx40_mmhg"] == 3.0654

    _dz_25 = _h4_25["design_z_and_contrast_b3_unrounded_min_face_audit"]
    assert _dz_25["pb30_100um_d_chronic"]["passes_unrounded_min_face_iip_0025"] is False
    assert _dz_25["pb30_100um_d_chronic"]["fv_harmonic_face_min_mmhg"] == 14.2854
    assert 30.70 <= _dz_25["unrounded_bed_po2_cutoffs_mmhg"]["mode_a_mixed_100um_d_chronic_72x72x42"] <= 30.75
    assert 35.55 <= _dz_25["unrounded_bed_po2_cutoffs_mmhg"]["mode_a_mixed_100um_d_chronic_div3_collagen_48x48x28"] <= 35.65
    assert _dz_25["pfc_5volpct_conditional_20x_solubility_reservoir"]["standalone_buffer_seconds"] == 4.525
    assert _dz_25["contrast_b3_donor_iec_omitted_pb30_intact_kperf_0p1358"]["100um"]["passes_unrounded_min_face_iip_0025"] is True
    assert _dz_25["contrast_b3_donor_iec_omitted_pb30_halved_kperf_0p0679"]["100um"]["passes_unrounded_min_face_iip_0025"] is False

    _p05_25 = _h4_25["gate_p05_two_timepoint_p04040_and_dual_basis_h2o2_audit"]
    assert [
        _p05_25["stationary_aged_repeat_burst_dual_basis_h2o2"]["basis_1500_ieq_cm2"][k]["gel_interface_z70p0um_max_h2o2_uM"]
        for k in ("retained_73800_iu", "undersized_51660_iu", "half_36900_iu", "undersized_26357_iu")
    ] == [0.7019, 1.297, 2.2076, 3.5984]
    assert [
        _p05_25["stationary_aged_repeat_burst_dual_basis_h2o2"]["basis_2500_ieq_cm2"][k]["gel_interface_z70p0um_max_h2o2_uM"]
        for k in ("retained_73800_iu", "undersized_51660_iu", "half_36900_iu", "undersized_26357_iu")
    ] == [1.17, 2.162, 3.68, 5.9989]

    # Gate P05 subgate bypass fix: adding design_z_divider_arm to a Design A/X payload NEVER suppresses bovine P00432 species hold:
    _bovine_bypass_attempt = wc_mod.evaluate_stage0_decision_tree({
        "product_variant_id": "ISLET1_A4_ARM1B",
        "p05_catalase_uniprot_accession": "P00432",
        "p05_recombinant_human_catalase_verified": False,
        "design_z_divider_arm": "PLAIN_HYDROGEL_140UM",
    })
    assert _bovine_bypass_attempt["resolved_architecture_assignments"]["p05_catalase_species_gate_status"] == "HOLD_P05_XENO_OR_UNVERIFIED_CATALASE_SPECIES_REQUIRES_HUMAN_P04040"
    assert _bovine_bypass_attempt["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_P05_XENO_OR_UNVERIFIED_CATALASE_SPECIES_REQUIRES_HUMAN_P04040"

    # Registered CaO2 product running P05 barrier panel without P04040 species panel cannot clear P05 subgate:
    _omit_species_attempt = wc_mod.evaluate_stage0_decision_tree({
        "product_variant_id": "ISLET1_A4_X_MOL_RNLS_WT",
        "p05_barrier_slab_geometry_mode": "FIXED_350UM_2X105UM_BANDS",
        "p05_full_336h_coupled_o2_h2o2_ca_ph_qualified": True,
        "obs_p05_microbead_volume_d99_um": 42.0,
        "obs_p05_microbead_volume_d99_se": 1.0,
        "obs_p05_catalase_barrier_film_um": 24.0,
        "obs_p05_catalase_barrier_film_se": 0.5,
        "obs_p05_islet_plane_h2o2_max_uM": 1.1,
        "obs_p05_islet_plane_h2o2_max_se": 0.05,
        "obs_p05_day7_catalase_retention_pct": 80.0,
        "obs_p05_day7_catalase_retention_se": 1.0,
        "obs_p05_336h_retained_catalase_iu_333cm2": 77171.25,
        "obs_p05_336h_retained_catalase_iu_333cm2_se": 300.0,
    })
    assert _omit_species_attempt["resolved_architecture_assignments"]["p05_catalase_species_gate_status"] == "HOLD_P05_XENO_OR_UNVERIFIED_CATALASE_SPECIES_REQUIRES_HUMAN_P04040"

    # Two-timepoint human P04040 catalase envelope + unrounded endotoxin (0.0550474255 vs 0.05505) + z=70.0 um H2O2:
    _m_9083 = 158300.0 / 17428.015  # 9.0830769... mg
    _t95_df5 = wc_mod.student_t_ppf_one_sided_95(5.0)
    _eu_fail_se = 0.001
    _eu_fail_mean = 0.05505 - _t95_df5 * _eu_fail_se  # U95_t == 0.05505 > 0.0550474255!
    assert wc_mod.evaluate_coupled_catalase_endotoxin_gate_iip0025(
        m_total_cat_mg=_m_9083,
        eu_per_mg_mean=_eu_fail_mean,
        eu_per_mg_se=_eu_fail_se,
        total_implant_eu_mean=_m_9083 * _eu_fail_mean,
        total_implant_eu_se=0.005,
        replicates_n=6,
    )["status"] == "HOLD_P05_MASS_COUPLED_CATALASE_ENDOTOXIN_CAP_EXCEEDED"
    assert wc_mod.evaluate_coupled_catalase_endotoxin_gate_iip0025(
        m_total_cat_mg=_m_9083,
        eu_per_mg_mean=0.040,
        eu_per_mg_se=0.002,
        total_implant_eu_mean=_m_9083 * 0.040,
        total_implant_eu_se=0.010,
        replicates_n=6,
    )["gate_pass"] is True

    _valid_p05_25 = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "obs_stage0_independent_lot_replicates_n": 6.0,
        "obs_stage0_independent_lot_replicates_n_se": 1e-12,
        "p05_catalase_uniprot_accession": "P04040",
        "p05_recombinant_human_catalase_verified": True,
        "p05_two_timepoint_same_lot_envelope_verified": True,
        "obs_p05_precure_active_catalase_input_u_333cm2": 156000.0,
        "obs_p05_precure_active_catalase_input_se_u": 200.0,
        "obs_p05_cure_immobilization_yield_pct": 65.0,
        "obs_p05_cure_immobilization_yield_se": 0.5,
        "obs_p05_postcure_t0_active_catalase_iu_333cm2": 101400.0,
        "obs_p05_postcure_t0_active_catalase_se_iu": 200.0,
        "obs_p05_day7_catalase_retention_pct": 80.0,
        "obs_p05_day7_catalase_retention_se": 1.0,
        "obs_p05_day7_shed_catalase_mass_pct": 18.0,
        "obs_p05_day7_shed_catalase_mass_se": 1.0,
        "obs_p05_336h_catalase_retention_pct": 76.0,
        "obs_p05_336h_catalase_retention_se": 1.0,
        "obs_p05_336h_retained_catalase_iu_333cm2": 77064.0,
        "obs_p05_336h_retained_catalase_iu_333cm2_se": 300.0,
        "obs_p05_z70_gel_interface_h2o2_max_uM": 1.17,
        "obs_p05_z70_gel_interface_h2o2_max_se": 0.05,
        "obs_p05_total_loaded_catalase_protein_mg": 8.95,
        "obs_p05_total_loaded_catalase_protein_mg_se": 0.05,
        "obs_p05_catalase_endotoxin_eu_per_mg": 0.040,
        "obs_p05_catalase_endotoxin_eu_per_mg_se": 0.002,
        "obs_p05_total_implant_catalase_endotoxin_eu": 0.358,
        "obs_p05_total_implant_catalase_endotoxin_eu_se": 0.010,
    }
    _ev_p05_ok = wc_mod.evaluate_stage0_decision_tree(_valid_p05_25)
    assert _ev_p05_ok["resolved_architecture_assignments"]["hurdle4_spatial_oxygen_catalase_and_route_evaluation_iip_0025"]["gate_pass"] is True
    assert wc_mod.evaluate_stage0_decision_tree({
        **_valid_p05_25,
        "obs_p05_precure_active_catalase_input_u_333cm2": 295200.0,
        "obs_p05_postcure_t0_active_catalase_iu_333cm2": 191880.0,
        "obs_p05_336h_retained_catalase_iu_333cm2": 145828.8,
    })["resolved_architecture_assignments"]["hurdle4_spatial_oxygen_catalase_and_route_challengers_iip_0025"] == "HOLD_P05_PRECURE_CATALASE_INPUT_EXCEEDS_158300U_CAP"
    assert wc_mod.evaluate_stage0_decision_tree({
        **_valid_p05_25,
        "obs_p05_z70_gel_interface_h2o2_max_uM": 5.9989,
    })["safety_kill_s9_triggered"] is True

    # Unified Student-t Day-90 Gate 0S retrieval across wetlab_contract.py AND ingest_wetlab_results.py + dual-assay shedding:
    assert wetlab._explant_recovery_clears_student_t_0p99({
        "obs_explant_en_bloc_recovery_fraction": 0.995,
        "obs_explant_en_bloc_recovery_se": 0.0028,
    }) is False
    assert wetlab._explant_recovery_clears_student_t_0p99({
        "obs_explant_en_bloc_recovery_fraction": 0.995,
        "obs_explant_en_bloc_recovery_se": 0.001,
    }) is True
    assert wc_mod.evaluate_stage0_decision_tree({
        "obs_explant_en_bloc_recovery_fraction": 0.995,
        "obs_explant_en_bloc_recovery_se": 0.0028,
    })["safety_kill_s9_triggered"] is True
    assert wc_mod.evaluate_dual_assay_gate0s_shedding_iip0025(
        ddpcr_k=0,
        ddpcr_n=2995731,
        ddpcr_specimen_unit_id="SISTER-DDPCR-01",
        viable_k=0,
        viable_n=2995731,
        viable_specimen_unit_id="SISTER-VIABLE-02",
        day90_recovery_mean=0.995,
        day90_recovery_se=0.001,
        day90_replicates_n=6,
    )["gate_pass"] is True
    assert wc_mod.evaluate_dual_assay_gate0s_shedding_iip0025(
        ddpcr_k=1,
        ddpcr_n=5000000,
        ddpcr_specimen_unit_id="SISTER-DDPCR-01",
        viable_k=0,
        viable_n=5000000,
        viable_specimen_unit_id="SISTER-VIABLE-02",
    )["status"] == "FAIL_HARD_STOP"
    assert wc_mod.evaluate_dual_assay_gate0s_shedding_iip0025(
        ddpcr_k=0,
        ddpcr_n=2995731,
        ddpcr_specimen_unit_id="SAME-UNIT",
        viable_k=0,
        viable_n=2995731,
        viable_specimen_unit_id="SAME-UNIT",
    )["status"] == "HOLD_GATE0S_DUAL_SHEDDING_REQUIRES_INDEPENDENT_SPECIMEN_UNITS"
    assert wc_mod.evaluate_dual_assay_gate0s_shedding_iip0025(
        ddpcr_k=0,
        ddpcr_n=2995731,
        ddpcr_specimen_unit_id="SISTER-DDPCR-01",
        viable_k=0,
        viable_n=2995731,
        viable_specimen_unit_id="SISTER-VIABLE-02",
        max_u95_fraction=1e-4,
    )["status"] == "HOLD_CALLER_CAP_OR_PARAMETER_OVERRIDE_FORBIDDEN"

    # Contrasts B1, B2, B3 & GY-020:
    assert wc_mod.evaluate_stage0_decision_tree({
        "contrast_b2_perfusion_and_retrieval_arm": "PAIRED_HUMAN_EXVIVO_AND_PORCINE_DAY30_90",
        "contrast_b2_route_qualification_distinct_from_portfolio_credit_verified": False,
        "obs_contrast_b2_human_exvivo_viable_window_h": 36.0,
        "obs_contrast_b2_human_exvivo_viable_window_se_h": 1.0,
        "obs_contrast_b2_omental_vs_suprafascial_po2_gain_mmhg": 8.0,
        "obs_contrast_b2_omental_vs_suprafascial_po2_gain_se": 0.5,
        "obs_contrast_b2_porcine_day90_en_bloc_recovery_fraction": 0.996,
        "obs_contrast_b2_porcine_day90_en_bloc_recovery_se": 0.001,
    })["resolved_architecture_assignments"]["hurdle4_spatial_oxygen_catalase_and_route_challengers_iip_0025"] == "HOLD_CONTRAST_B2_ROUTE_QUALIFICATION_CONFLATED_WITH_PORTFOLIO_CREDIT"

    _gy_25 = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    assert len(_gy_25["graveyard_entries"]) == 10
    _g20 = [g for g in _gy_25["off_grid_molecular_dead_end_guardrails"]["guardrails"] if g["graveyard_card_id"] == "GY-020"]
    assert len(_g20) == 1 and _g20[0]["severity_level"] == "BLOCK" and _g20[0]["contribution_id"] == "IIP-0025"

    # =========================================================================
    # 51. 69-Key gnomAD v4.1 Join, 17-Record HLA-C Scoping (114/115 vs 106/108),
    #     4-State 2A-Free Design B Comparator FASTA, Transcribed-Interval Layer-2
    #     ORF Scanner, 3-Path Stage-0 Parity & INS 3'-UTR Retirement (IIP-0026; DEC-105)
    # =========================================================================
    from scripts import run_layer2_perturbation_and_regulatory_plane as l2_mod_26
    from scripts import run_layer5_personalized_patient_compiler as l5_mod_26

    # (a) Staged Design B FASTA bytes, 4-state registry binding, & canonical Scenario A FASTA invariance:
    _db_fasta_path = ROOT / "data/constructs/islet1_design_b_complete_constructs_v0.8.0.fasta"
    _db_fasta_bytes = _db_fasta_path.read_bytes()
    assert len(_db_fasta_bytes) == 25802
    assert hashlib.sha256(_db_fasta_bytes).hexdigest() == "4bb20389773c76b2705031ba675487825dfbd94d920a692cde75b7d56af22aa4"
    assert hashlib.sha256((ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta").read_bytes()).hexdigest() == "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"
    _pdefs_26 = pc.load_product_contracts(ROOT)
    _fam_b_26 = [p for p in _pdefs_26["products"] if p.get("family") == "B"]
    assert len(_fam_b_26) == 11 and all(p["execution_supported"] is False for p in _fam_b_26)
    for _b2_id in ("ISLET1_B2_5LOCUS_BASE", "ISLET1_B2_5LOCUS_TYK2_PE"):
        _b2_def = pc.get_product_definition(_b2_id, root=ROOT)
        assert _b2_def["lifecycle_state"] == "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"
        assert _b2_def["staged_design_b_donor_reference_binding"]["bound_c1_c5_pair_sha256"] == "82c351af082e0f97361e715c36e2ab5872b8d097efa1c9f59dfda53ea60496ea"

    # (b) Layer-2 24-segment provenance (c6 C72T/C73A = U78775.2 C79T/C80A), 975-nt c1/c2 counterexample, & OQ-011 retirement:
    _l2_full_26 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    _m26_l2 = _l2_full_26["stage0_design_b_comparator_provenance_and_ins_retirement_manifest_iip_0026"]
    assert _m26_l2["contribution_id"] == "IIP-0026" and _m26_l2["decision_id"] == "DEC-105"
    assert len(_m26_l2["per_segment_provenance_ledger_24"]) == 24
    _c6_chs4_seg = next(s for s in _m26_l2["per_segment_provenance_ledger_24"] if s["segment_id"] == "chs4_c72t_c73a_250")
    assert "c6" in _c6_chs4_seg["used_in"] and "C79T/C80A" in _c6_chs4_seg["transformations"]
    assert _m26_l2["c6_two_factor_comparator_audit"]["exact_two_factor_length_identity"] is True
    _sp_rec_26 = _m26_l2["swissprot_unique_string_vs_positional_reconciliation"]
    assert _sp_rec_26["arm6_serpinb9_c341s_missense_benchmark"] == {
        "legacy_integer_count_delta_8_11mer": 0,
        "unique_8_11mers_added_absent_from_swissprot": 38,
        "unique_8_11mers_removed_vs_wt": 38,
        "unique_15mers_added_absent_from_swissprot": 15,
    }
    assert _sp_rec_26["complete_islet1_b2_5locus_base_mature_8_11_and_15"] == [108, 39]
    assert _sp_rec_26["complete_islet1_b2_5locus_base_precursor_8_11_and_15"] == [109, 44]
    assert _sp_rec_26["complete_islet1_b2_5locus_tyk2_pe_mature_8_11_and_15"] == [146, 54]
    assert _sp_rec_26["complete_islet1_b2_5locus_tyk2_pe_precursor_8_11_and_15"] == [147, 59]
    assert len(_m26_l2["iip0015_pdl1_syn2_syn4_versioned_disclosure_note"]["iip0018_multimir_6_descendants_minus_orfs_ge_300nt"]) == 6
    _oq11_26 = _m26_l2["oq011_ins_3utr_eif4g2_knockin_retirement_audit"]
    assert _oq11_26["graveyard_card_id"] == "GY-022" and _oq11_26["expansion_fold_vs_73nt"] == 37.14
    assert abs(_oq11_26["curative_500k_ieq_non_beta_unshielded_cell_fraction"] - 0.485806) < 1e-5

    _c1_rec = next(r for r in _m26_l2["design_b_construct_audits"] if r["short_id"] == "c1")
    assert _c1_rec["p1_to_p4"]["gate_P2"]["upstream_regulatory_kozak_gate_status"] == "HOLD_PENDING_5PRIME_RACE_TSS"
    assert _c1_rec["p1_to_p4"]["gate_P2"]["upstream_regulatory_kozak_runs_gated_on_5prime_race"][0]["coords_1based"] == [2308, 3282]
    assert _c1_rec["p1_to_p4"]["gate_P2"]["upstream_regulatory_kozak_runs_gated_on_5prime_race"][0]["span_nt"] == 975

    # Direct Layer-2 validator probes on c1 DNA (wrong SHA, upstream TSS <= 2308 FAIL vs > 2308 PASS, internal 357-nt +1 ORF FAIL):
    _c1_dna = next(
        seq for hdr, seq in l2_mod_26._parse_fasta_records_iip0026(_db_fasta_path)
        if hdr.startswith("DESIGNB_L5A_CXCL10_BFB_SYN6_ACC_4078bp")
    )
    _c1_sha = hashlib.sha256(_c1_dna.encode("ascii")).hexdigest()
    _l2_ref_raw = json.loads(l2_mod_26.REF_PATH.read_text(encoding="utf-8"))
    _sb9_aa = _l2_ref_raw["uniprot_protein_sequences"]["P50453"]["sequence"]
    _pdl1_aa = _l2_ref_raw["uniprot_protein_sequences"]["Q9NZQ7"]["sequence"]
    _c1_orfs = [
        {"name": "SERPINB9_P50453", "start_0based": 888, "end_0based": 2019, "strand": "+", "has_stop": True, "expected_precursor_aa": _sb9_aa, "mature_slice_0based": (0, 376), "reference_key": "P50453_SERPINB9"},
        {"name": "CD274_PDL1_Q9NZQ7", "start_0based": 2380, "end_0based": 3253, "strand": "+", "has_stop": True, "expected_precursor_aa": _pdl1_aa, "mature_slice_0based": (18, 290), "reference_key": "Q9NZQ7_PDL1"},
    ]
    _c1_ti = _c1_rec["p1_to_p4"]["gate_P1"]["transcribed_intervals"]
    _c1_pin = {"P50453_SERPINB9": _sb9_aa, "Q9NZQ7_PDL1": _pdl1_aa}
    with pytest.raises(ValueError, match="FAIL_CONSTRUCT_SHA256_MISMATCH"):
        l2_mod_26.validate_design_b_transcribed_construct_iip_0026(
            design_id="C1_BAD_SHA",
            donor_seq=_c1_dna,
            cargo_bounds_0based=(600, len(_c1_dna) - 600),
            orf_specs=_c1_orfs,
            parent_mature_proteins=[_sb9_aa, _pdl1_aa[18:]],
            swissprot_k8_11=set(),
            swissprot_k15=set(),
            pinned_layer3_sequences=_c1_pin,
            expected_donor_sha256="0" * 64,
            transcribed_intervals=_c1_ti,
        )
    with pytest.raises(ValueError, match="FAIL_TRANSCRIBED_UPSTREAM_KOZAK_ORF_BY_5PRIME_RACE"):
        l2_mod_26.validate_design_b_transcribed_construct_iip_0026(
            design_id="C1_RACE_UPSTREAM_2300",
            donor_seq=_c1_dna,
            cargo_bounds_0based=(600, len(_c1_dna) - 600),
            orf_specs=_c1_orfs,
            parent_mature_proteins=[_sb9_aa, _pdl1_aa[18:]],
            swissprot_k8_11=set(),
            swissprot_k15=set(),
            pinned_layer3_sequences=_c1_pin,
            expected_donor_sha256=_c1_sha,
            transcribed_intervals=_c1_ti,
            observed_5prime_race_end_bp=2300,
            observed_3prime_race_end_bp=3478,
            readthrough_umi_n_eff=3689,
            readthrough_umi_k=0,
        )
    _c1_race_ok = l2_mod_26.validate_design_b_transcribed_construct_iip_0026(
        design_id="C1_RACE_DOWNSTREAM_2350",
        donor_seq=_c1_dna,
        cargo_bounds_0based=(600, len(_c1_dna) - 600),
        orf_specs=_c1_orfs,
        parent_mature_proteins=[_sb9_aa, _pdl1_aa[18:]],
        swissprot_k8_11=set(),
        swissprot_k15=set(),
        pinned_layer3_sequences=_c1_pin,
        expected_donor_sha256=_c1_sha,
        transcribed_intervals=_c1_ti,
        observed_5prime_race_end_bp=2350,
        observed_3prime_race_end_bp=3478,
        readthrough_umi_n_eff=3689,
        readthrough_umi_k=0,
    )
    assert _c1_race_ok["gate_P2"]["upstream_regulatory_kozak_gate_status"] == "CLEARED_DOWNSTREAM_5PRIME_RACE_TSS_VERIFIED"
    assert _c1_race_ok["gate_P1"]["transcript_boundary_gate_status"] == "VERIFIED_OBSERVED_ENDS"
    # Synthetic internal out-of-frame 357-nt +1 [AG]NNATGG ORF wholly inside a 600-nt CDS must fail closed:
    _syn_cds_600 = "ATG" + ("GCT" * 3) + "GAGAATGGG" + ("CTG" * 116) + "GTAACT" + ("GCT" * 74) + "TAA"
    assert len(_syn_cds_600) == 600
    _syn_donor = ("A" * 600) + _syn_cds_600 + ("A" * 600)
    _syn_aa = l2_mod_26._translate_strict_cds_iip0026(_syn_cds_600, True)
    with pytest.raises(ValueError, match="FAIL_UNDECLARED_KOZAK_ORF"):
        l2_mod_26.validate_design_b_transcribed_construct_iip_0026(
            design_id="SYN_INTERNAL_PLUS1_357NT",
            donor_seq=_syn_donor,
            cargo_bounds_0based=(600, 1200),
            orf_specs=[{"name": "SYN_CDS", "start_0based": 600, "end_0based": 1200, "strand": "+", "has_stop": True, "expected_precursor_aa": _syn_aa, "mature_slice_0based": (0, len(_syn_aa)), "reference_key": "SYN_REF"}],
            parent_mature_proteins=[_syn_aa],
            swissprot_k8_11=set(),
            swissprot_k15=set(),
            pinned_layer3_sequences={"SYN_REF": _syn_aa},
            expected_donor_sha256=hashlib.sha256(_syn_donor.encode("ascii")).hexdigest(),
        )

    # (c) Layer-5 69-key gnomAD v4.1 join, 17-record HLA-C (114/115 vs 106/108), typed CIITA Exon-5 routing, & outward odds invariance:
    _l5_full_26 = json.loads((ROOT / "data/extracted/layer5_personalized_patient_and_crispr_atlas.json").read_text(encoding="utf-8"))
    _a26_l5 = _l5_full_26["variant_aware_gnomad_hlac_pe4max_and_designy_audit_iip_0026"]
    assert _a26_l5["contribution_id"] == "IIP-0026" and _a26_l5["decision_id"] == "DEC-105"
    assert _a26_l5["quarantined_doench2016_weight_table_entries_bundled"] == 0
    _g69 = _a26_l5["gnomad_v4_1_69key_normalized_join"]
    assert _g69["total_dbsnp155_loci_audited"] == 56 and _g69["total_normalized_alleles_audited"] == 69
    assert _g69["either_modality_denominators"] == {"resolved_in_either_wgs_or_exome_alleles": 49, "unresolved_absent_in_both_modalities_alleles": 20}
    assert _g69["wgs_only_denominators"] == {"pass_wgs_alleles": 26, "flagged_wgs_alleles": 1, "unresolved_wgs_alleles": 42}
    assert _g69["exome_only_denominators"] == {"pass_exome_alleles": 43, "flagged_exome_alleles": 2, "unresolved_exome_alleles": 24}

    assert l5_mod_26.left_normalize_vcf_allele_iip_0026("chr1", 104, 106, "AA", "A", 101, "CAAAAACT")[2] == 101
    assert l5_mod_26.left_normalize_vcf_allele_iip_0026("chr1", 104, 106, "AA", "A", 101, "AAAAAACT")[0] == "UNRESOLVED_BOUNDARY"
    _m_lt_1e4 = l5_mod_26.summarize_gnomad_modality_iip_0026({"ac": 1, "an": 10001, "filters": []})
    _m_gt_1e4 = l5_mod_26.summarize_gnomad_modality_iip_0026({"ac": 1, "an": 9999, "filters": []})
    assert _m_lt_1e4["global_af"] == 0.0001 and (_m_lt_1e4["ac"] * 10000 >= _m_lt_1e4["an"]) is False
    assert _m_gt_1e4["global_af"] == 0.0001 and (_m_gt_1e4["ac"] * 10000 >= _m_gt_1e4["an"]) is True
    assert l5_mod_26.summarize_gnomad_modality_iip_0026({"ac": 5, "an": 10000})["resolution_status"] == "UNRESOLVED"

    _hlac_17 = _a26_l5["hlac_17_record_and_hla_ab_exon23_candidate_audit"]
    assert _hlac_17["total_sequence_records_evaluated"] == 17
    assert _hlac_17["exon_boundary_counting_conventions"]["HLA_A"] == {"nominal_spcas9_cut_in_exon23_count": 114, "full_23mer_wholly_inside_exon23_count": 106}
    assert _hlac_17["exon_boundary_counting_conventions"]["HLA_B"] == {"nominal_spcas9_cut_in_exon23_count": 115, "full_23mer_wholly_inside_exon23_count": 108}
    assert _hlac_17["hla_a_exon23_all_nnn_audit"]["exon23_window_lt005_across_17_sampled_records_count"] == 11
    assert _hlac_17["hla_b_exon23_all_nnn_audit"]["exon23_window_lt005_across_17_sampled_records_count"] == 6
    assert len(_hlac_17["full_record_10_to_11_mismatch_non_exon_hits"]) == 3
    assert "21248048" in _hlac_17["dt9_hlae_crossreactivity_specificity_control_rule"]
    assert _a26_l5["mode_a_vs_held_mode_b_pam_shields_and_12_alignments"]["nominated_grch38_alignments_12_across_11_windows"]["total_alignments"] == 12
    assert _a26_l5["mode_a_vs_held_mode_b_pam_shields_and_12_alignments"]["nominated_grch38_alignments_12_across_11_windows"]["unique_genomic_windows"] == 11
    assert round(_a26_l5["oq016_2pow12_portfolio_invariance"]["design_a_standalone_pct"]) == 27
    assert round(_a26_l5["oq016_2pow12_portfolio_invariance"]["design_b_standalone_pct"]) == 31
    assert round(_a26_l5["oq016_2pow12_portfolio_invariance"]["design_x_standalone_pct"]) == 36
    assert round(_a26_l5["oq016_2pow12_portfolio_invariance"]["three_slot_locked_pct"]) == 53

    _comp_b_26 = l5_mod_26.compile_personalized_patient_blueprint(
        {
            "patient_id": "PT-IIP0026-PLANB",
            "ancestry": "EUR",
            "nkg2a_neg_kir_educated_nk_fraction": 0.55,
            "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.0,
            "ifn1_ish_genetic_risk_percentile": 40.0,
            "sirpa_genotype": "v1/v1",
            "recipient_kir_ligand_group": "C1/C1",
            "hlac_4digit_relation": "SELF_4DIGIT_MATCH",
        },
        track2_subbank_mode="ENDOGENOUS_HLAC_EXPLORATORY",
    )
    assert "OI-GRNA-02_CIITA_EXON5_ASCAS12A" in _comp_b_26["compiled_grna_panel"]
    assert "OI-GRNA-02_CIITA_EXON1_ASCAS12A" not in _comp_b_26["compiled_grna_panel"]
    assert _comp_b_26["mhc_locus_routing_contract_iip_0026"]["ciita_target_locus_id"] == "OI-GRNA-02_CIITA_EXON5_ASCAS12A"
    assert _l5_full_26["part_b_personalized_patient_compiler_and_10k_population_trial"]["global_5yr_survival_personalized_dna_compiled_openislet_pct"] == 97.4

    # (d) S7 missing-evidence & SCT-without-MLC fixes + 19-case 3-Path Direct/Ingest/Replay branch & lot parity:
    assert wc_mod.evaluate_stage0_decision_tree({})["switch_decisions"]["S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B"] == "UNMEASURED_EXTENDED_S7_GATE0B_HLAC_PANEL"
    assert wc_mod.evaluate_stage0_decision_tree({"product_variant_id": "ISLET1_B2_5LOCUS_BASE"})["switch_decisions"]["S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B"] == "HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE"
    _sct_no_mlc = wc_mod.evaluate_stage0_decision_tree({
        "product_variant_id": "ISLET1_B_SCT_C1_CIITA",
        "obs_gate0b_hlac_cd8_elispot_sfw": 150.0,
        "obs_gate0b_hlac_cd8_elispot_se": 5.0,
        "obs_gate0b_hlac_sct_molecules_per_cell": 25000.0,
        "obs_gate0b_hlac_sct_molecules_se": 500.0,
        "b_prime_sct_recipient_cd8_mlc_cleared": False,
    })
    assert _sct_no_mlc["switch_decisions"]["S7_Design_B_Prime_4Digit_Allo_vs_B_SCT_Gate_0B"] == "HOLD_B_PRIME_ALLO_FAILED_REQUIRES_4DIGIT_MATCH_OR_SCT_CD8_MLC_CLEARANCE"

    _base_s9_26 = {
        "obs_icasp9_log10_kill": 4.6,
        "obs_icasp9_log10_kill_se": 0.08,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000.0,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000.0,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300000000,
    }
    def _valid_t6b_26(pvid: str, lot: str, csha: str = wc_mod.DESIGN_B_C1_C5_PAIR_SHA256, camp: str = "CAMP_S0_IIP0026") -> list:
        return [
            {
                "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
                "scope_type": "LOT_SCOPE",
                "kind": "EXACT_INTEGER_COUNT",
                "timepoint": "STAGE7_PRE_ASSEMBLY",
                "assay_version": "TIER6B_SORTED_LOT_V1",
                "unit": "clonogenic_opportunities",
                "product_variant_id": pvid,
                "construct_sha256": csha,
                "donor_seedline_id": "SEED_IIP0026_01",
                "mcb_id": "MCB_IIP0026_01",
                "lot_id": lot,
                "campaign_id": camp,
                "sample_id": f"SMP_{lot}_ENR",
                "fraction_id": "enriched",
                "declared_fraction_roster": ["enriched", "flow_through"],
                "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
                "clonogenic_risk_unit": "CD31_NEG_CD144_NEG_CLONOGENIC_OPPORTUNITY",
                "endogenous_denominator_verified": True,
                "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
                "untreated_cfu_control": {"untreated_sample_id": f"UNTR_{lot}_ENR", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
                "presort_mass_balance_reconciled": True,
                "unreconciled_loss_fraction": 0.0,
                "presort_total_input_cells": 3_000_000_000,
                "presort_reconciled_loss_cells": 0,
                "presort_fraction_cells": 2_999_700_000,
                "plated_treated_cells_n": 80_000_000,
                "w_i": 0.9999,
                "N": 30_000_000,
                "k": 0,
                "hpsc_calibrant": {"r_k": 20, "r_n": 20},
                "edited_differentiated_calibrant": {"r_k": 20, "r_n": 20},
            },
            {
                "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
                "scope_type": "LOT_SCOPE",
                "kind": "EXACT_INTEGER_COUNT",
                "timepoint": "STAGE7_PRE_ASSEMBLY",
                "assay_version": "TIER6B_SORTED_LOT_V1",
                "unit": "clonogenic_opportunities",
                "product_variant_id": pvid,
                "construct_sha256": csha,
                "donor_seedline_id": "SEED_IIP0026_01",
                "mcb_id": "MCB_IIP0026_01",
                "lot_id": lot,
                "campaign_id": camp,
                "sample_id": f"SMP_{lot}_FT",
                "fraction_id": "flow_through",
                "declared_fraction_roster": ["enriched", "flow_through"],
                "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
                "clonogenic_risk_unit": "CD31_NEG_CD144_NEG_CLONOGENIC_OPPORTUNITY",
                "endogenous_denominator_verified": True,
                "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
                "untreated_cfu_control": {"untreated_sample_id": f"UNTR_{lot}_FT", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
                "presort_mass_balance_reconciled": True,
                "unreconciled_loss_fraction": 0.0,
                "presort_total_input_cells": 3_000_000_000,
                "presort_reconciled_loss_cells": 0,
                "presort_fraction_cells": 300_000,
                "plated_treated_cells_n": 260_000,
                "w_i": 0.0001,
                "N": 100_000,
                "k": 0,
                "hpsc_calibrant": {"r_k": 10, "r_n": 10},
                "edited_differentiated_calibrant": {"r_k": 10, "r_n": 10},
            },
        ]

    _complete_b2_26 = {
        **_base_s9_26,
        "product_variant_id": "ISLET1_B2_5LOCUS_BASE",
        "lot_id": "LOT_B2_IIP0026",
        "campaign_id": "CAMP_S0_IIP0026",
        "bound_c1_c5_pair_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
        "construct_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
        "donor_seedline_id": "SEED_IIP0026_01",
        "mcb_id": "MCB_IIP0026_01",
        "tier6b_sorted_lot_observations": _valid_t6b_26("ISLET1_B2_5LOCUS_BASE", "LOT_B2_IIP0026"),
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_per_locus_segmented_baf": {"chr19q": 0.502, "chr11p15": 0.499, "chr6p21": 0.504},
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "design_b_hlac_phased_guide_cfd_max": 0.0188,
        "design_b_hla_ab_separated_pass_delivery_verified": True,
        "design_b_chr6p21_1p41mb_loh_sv_cleared": True,
        "design_b_hlac_dt9_positive": True,
        "design_b_hlae_3d12_positive": True,
        "design_b_hlac_dt9_matched_hlae_only_control_delta_ratio": 4.25,
        "design_b_hlac_allele_specific_surface_or_phased_mrna_verified": True,
        "design_b_locus5_5prime_race_coord_bp": 880,
        "design_b_locus5_3prime_race_coord_bp": 3478,
        "design_b_locus5_pdl1_5prime_race_coord_bp": 2350,
        "design_b_locus5_race_umi_n_eff": 3689,
        "design_b_locus5_aberrant_umi_k": 0,
    }

    def _assert_3path_parity_26(case_payload: dict) -> dict:
        d_out = wc_mod.evaluate_stage0_research_and_molecular_guardrails(dict(case_payload), root=ROOT)
        scope_26 = {
            "product_variant_id": case_payload.get("product_variant_id"),
            "lot_id": case_payload.get("lot_id"),
            "campaign_id": case_payload.get("campaign_id"),
            **({"iip0026_schema_era": case_payload["iip0026_schema_era"]} if "iip0026_schema_era" in case_payload else {}),
        }
        meas_26 = {k: v for k, v in case_payload.items() if str(k).startswith("obs_") and v is not None}
        mol_26 = {
            k: v for k, v in case_payload.items()
            if k not in ("product_variant_id", "lot_id", "campaign_id", "iip0026_schema_era") and not str(k).startswith("obs_")
        }
        i_out = wetlab.ingest_stage0_research_or_wetlab_batch(
            json.dumps({"scope": scope_26, "measurements": meas_26, "iip0026_molecular_observations": mol_26}),
            root=ROOT,
        )
        mol_items = list(mol_26.items())
        mid_i = len(mol_items) // 2
        r_out = wetlab.replay_stage0_chronological_batches(
            [
                {"seq_index": 3, "utc_timestamp": "2026-11-18T10:00:00Z", "scope": dict(scope_26), "measurements": {}, "iip0026_molecular_observations": dict(mol_items[mid_i:])},
                {"seq_index": 1, "utc_timestamp": "2026-11-18T08:00:00Z", "scope": dict(scope_26), "measurements": dict(meas_26), "iip0026_molecular_observations": {}},
                {"seq_index": 2, "utc_timestamp": "2026-11-18T09:00:00Z", "scope": dict(scope_26), "measurements": {}, "iip0026_molecular_observations": dict(mol_items[:mid_i])},
            ],
            root=ROOT,
        )
        assert d_out["overall_status"] == i_out["overall_status"] == r_out["overall_status"]
        assert d_out["lot_release_readiness_tier"] == i_out["lot_release_readiness_tier"] == r_out["lot_release_readiness_tier"]
        assert d_out["iip0026_reasons"] == i_out["iip0026_reasons"] == r_out["iip0026_reasons"]
        assert d_out["lifecycle_state"] == i_out["lifecycle_state"] == r_out["lifecycle_state"]
        assert d_out["execution_supported"] == i_out["execution_supported"] == r_out["execution_supported"]
        assert d_out["effective_execution_permitted"] == i_out["effective_execution_permitted"] == r_out["effective_execution_permitted"]
        assert d_out["switch_decisions"] == i_out["switch_decisions"] == r_out["switch_decisions"]
        return d_out

    _res_b2_ok = _assert_3path_parity_26(_complete_b2_26)
    assert _res_b2_ok["overall_status"] == "RESEARCH_COMPARATOR_STAGE0_PANEL_CLEARED_PENDING_STATE4_PROMOTION"
    assert _res_b2_ok["execution_supported"] is False and _res_b2_ok["effective_execution_permitted"] is False
    assert _assert_3path_parity_26({**_complete_b2_26, "design_b_locus5_pdl1_5prime_race_coord_bp": 2300})["overall_status"] == "FAIL_MATURE_TRANSCRIPT_INCLUDES_UPSTREAM_OUT_OF_FRAME_KOZAK_ORF"
    assert _assert_3path_parity_26({k: v for k, v in _complete_b2_26.items() if k != "donor_seed_gnomad_target_site_wgs_verified"})["overall_status"] == "HOLD_MISSING_DONOR_TARGET_WGS"
    assert _assert_3path_parity_26({k: v for k, v in _complete_b2_26.items() if k != "design_b_hlac_phased_guide_cfd_max"})["overall_status"] == "HOLD_MISSING_ENDOGENOUS_DESIGNB_HLAC_PHASED_CLEARANCE"
    assert _assert_3path_parity_26({**_complete_b2_26, "post_edit_per_locus_segmented_baf": {"chr19q": 0.30, "chr11p15": 0.70, "chr6p21": 0.50}})["overall_status"] == "FAIL_DISCORDANT_PER_LOCUS_SEGMENTED_BAF"
    assert _assert_3path_parity_26({**_complete_b2_26, "design_b_locus5_race_umi_n_eff": 3687})["overall_status"] == "HOLD_MISSING_DESIGNB_LOCUS5_TRANSCRIPT_UMI_3689"
    assert _assert_3path_parity_26({**_complete_b2_26, "design_b_locus5_5prime_race_coord_bp": True, "design_b_locus5_3prime_race_coord_bp": -2})["overall_status"] == "FAIL_INVALID_RACE_COORDINATES"
    assert _assert_3path_parity_26({k: v for k, v in _complete_b2_26.items() if k not in ("design_b_hlac_dt9_matched_hlae_only_control_delta_ratio", "design_b_hlac_allele_specific_surface_or_phased_mrna_verified")})["overall_status"] == "HOLD_MISSING_HLAC_SPECIFICITY_CONTROL_FOR_DT9_HLAE_CROSSREACTIVITY"
    assert _assert_3path_parity_26({**_complete_b2_26, "design_b_hlac_dt9_matched_hlae_only_control_delta_ratio": 1.04, "design_b_hlac_allele_specific_surface_or_phased_mrna_verified": False})["overall_status"] == "FAIL_HLAC_SURFACE_LOSS_MASKED_BY_DT9_HLAE_CROSSREACTIVITY"
    assert _assert_3path_parity_26({**_complete_b2_26, "bound_c1_c5_pair_sha256": "0" * 64})["overall_status"] == "FAIL_DONOR_REGISTRY_DIGEST_MISMATCH"
    assert _assert_3path_parity_26({**_complete_b2_26, "iip0026_lifecycle_state": "State 4"})["overall_status"] == "FAIL_CALLER_LIFECYCLE_OR_EXECUTION_OVERRIDE_REJECTED"
    assert _assert_3path_parity_26({**_complete_b2_26, "lot_id": "  "})["overall_status"] == "HOLD_MISSING_LOT_OR_CAMPAIGN_SCOPE"
    assert _assert_3path_parity_26({**_complete_b2_26, "product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE", "tyk2_pe_modality": "PE4MAX_MSI_GATED", "tyk2_phased_haplotype_sap_verified": True})["overall_status"] == "HOLD_MISSING_MSI"
    _pe4_msi_fail = _assert_3path_parity_26({
        **_complete_b2_26,
        "product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE",
        "tyk2_pe_modality": "PE4MAX_MSI_GATED",
        "tyk2_phased_haplotype_sap_verified": True,
        "tyk2_pe4max_pentaplex_msi_shifts_k": 2,
        "tyk2_pe4max_wgs_repeat_indels_clone": 121,
        "tyk2_pe4max_wgs_repeat_indels_ctrl": 100,
        "tyk2_pe4max_coding_or_splice_de_novo_indels": 0,
    })
    assert _pe4_msi_fail["overall_status"] == "FAIL_MSI_DETECTED" and "FAIL_PE4MAX_REPEAT_INDEL_CAP_EXCEEDED" in _pe4_msi_fail["iip0026_reasons"]
    assert _assert_3path_parity_26({
        **_complete_b2_26,
        "product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE",
        "tyk2_pe_modality": "PE4MAX_MSI_GATED",
        "tyk2_phased_haplotype_sap_verified": True,
        "tyk2_pe4max_pentaplex_msi_shifts_k": 0,
        "tyk2_pe4max_wgs_repeat_indels_clone": 121,
        "tyk2_pe4max_wgs_repeat_indels_ctrl": 100,
        "tyk2_pe4max_coding_or_splice_de_novo_indels": 0,
    })["overall_status"] == "FAIL_PE4MAX_REPEAT_INDEL_CAP_EXCEEDED"
    assert _assert_3path_parity_26({
        **_base_s9_26,
        "product_variant_id": "ISLET1_B_SCT_C1_CIITA",
        "lot_id": "LOT_BSCT_26",
        "campaign_id": "CAMP_S0_26",
        "construct_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
        "donor_seedline_id": "SEED_IIP0026_01",
        "mcb_id": "MCB_IIP0026_01",
        "tier6b_sorted_lot_observations": _valid_t6b_26("ISLET1_B_SCT_C1_CIITA", "LOT_BSCT_26", camp="CAMP_S0_26"),
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_per_locus_segmented_baf": {"chr19q": 0.498, "chr11p15": 0.503},
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "b_sct_ciita_exon5_and_trimer_verified": True,
    })["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED"

    # Cross-lot mixing rejection & monotonic adverse MSI latching in replay_stage0_chronological_batches:
    assert wetlab.replay_stage0_chronological_batches(
        [
            {"seq_index": 1, "utc_timestamp": "2026-11-18T08:00:00Z", "scope": {"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "lot_id": "LOT_1", "campaign_id": "C1"}, "measurements": dict(_base_s9_26), "iip0026_molecular_observations": {}},
            {"seq_index": 2, "utc_timestamp": "2026-11-18T09:00:00Z", "scope": {"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "lot_id": "LOT_2", "campaign_id": "C1"}, "measurements": {}, "iip0026_molecular_observations": {}},
        ],
        root=ROOT,
    )["overall_status"] == "FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING"
    _latched_msi_26 = wetlab.replay_stage0_chronological_batches(
        [
            {
                "seq_index": 1,
                "utc_timestamp": "2026-11-18T08:00:00Z",
                "scope": {"product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE", "lot_id": "LOT_PE_LATCH", "campaign_id": "CAMP_PE"},
                "measurements": dict(_base_s9_26),
                "iip0026_molecular_observations": {
                    **{k: v for k, v in _complete_b2_26.items() if not str(k).startswith("obs_") and k not in ("product_variant_id", "lot_id", "campaign_id")},
                    "tyk2_pe_modality": "PE4MAX_MSI_GATED",
                    "tyk2_phased_haplotype_sap_verified": True,
                    "tyk2_pe4max_pentaplex_msi_shifts_k": 2,
                    "tyk2_pe4max_wgs_repeat_indels_clone": 12,
                    "tyk2_pe4max_wgs_repeat_indels_ctrl": 10,
                    "tyk2_pe4max_coding_or_splice_de_novo_indels": 0,
                },
            },
            {
                "seq_index": 2,
                "utc_timestamp": "2026-11-18T11:00:00Z",
                "scope": {"product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE", "lot_id": "LOT_PE_LATCH", "campaign_id": "CAMP_PE"},
                "measurements": {},
                "iip0026_molecular_observations": {"tyk2_pe4max_pentaplex_msi_shifts_k": 0},
            },
        ],
        root=ROOT,
    )
    assert _latched_msi_26["overall_status"] == "FAIL_MSI_DETECTED" and _latched_msi_26["effective_execution_permitted"] is False

    # (e) Graveyard GY-021 & GY-022 (10 combinatorial grid + >=12 off-grid):
    _gy_26 = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    _og_26 = _gy_26["off_grid_molecular_dead_end_guardrails"]["guardrails"]
    assert len(_gy_26["graveyard_entries"]) == 10 and len(_og_26) >= 12
    _gy21 = [g for g in _og_26 if g["graveyard_card_id"] == "GY-021"]
    _gy22 = [g for g in _og_26 if g["graveyard_card_id"] == "GY-022"]
    assert len(_gy21) == 1 and _gy21[0]["guardrail_id"].startswith("OFFGRID_GUARD_011_") and _gy21[0]["contribution_id"] == "IIP-0026"
    assert len(_gy22) == 1 and _gy22[0]["guardrail_id"].startswith("OFFGRID_GUARD_012_") and _gy22[0]["contribution_id"] == "IIP-0026"

    # =========================================================================
    # 52. DEC-106 / IIP-0028: POST-DEPLETION GATE 0D, INTACT-PLUS-INCREMENT
    #     CD16+ ADCC, UPSTREAM B2 PAIRED MCNEMAR / PHI_AB CONCORDANCE,
    #     MEASURED-POTENCY DEC-084 MIXTURE & 72H/D28/D84 DUAL-MODALITY IAPP
    # =========================================================================
    _nm_fasta_28 = ROOT / "data/constructs/nm_021999_5_archival_mrna.fasta"
    assert _nm_fasta_28.is_file()
    assert hashlib.sha256(_nm_fasta_28.read_bytes()).hexdigest() == wc_mod.ITM2B_NM021999_5_FASTA_SHA256
    _nm_mrna_28 = "".join(ln.strip().upper() for ln in _nm_fasta_28.read_text(encoding="utf-8").splitlines()[1:] if ln.strip())
    assert len(_nm_mrna_28) == 10089
    assert hashlib.sha256(_nm_mrna_28[155:956].encode("ascii")).hexdigest() == wc_mod.ITM2B_NM021999_5_CDS_SHA256
    assert hashlib.sha256((ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta").read_bytes()).hexdigest() == "b5ba4e2b4bbbcc7d6cf1138975c5fae5c6ab074dda43973f071c671cd8cd1d88"

    _l2_28_doc = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    _itm2b_28 = _l2_28_doc["stage0_p2_itm2b_nm021999_reserve_and_brichos_audit_iip_0028"]
    assert _itm2b_28["contribution_id"] == "IIP-0028" and _itm2b_28["decision_id"] == "DEC-106"
    assert _itm2b_28["active_weight_pct"] == 0.0
    assert _itm2b_28["archival_mrna_provenance"]["cds_ascii_sha256"] == wc_mod.ITM2B_NM021999_5_CDS_SHA256

    # (a) Registry-bound Gate 0d + intact-and-increment CD16+ ADCC + upstream B2 responder contract:
    _base_arm2_28 = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "mcb_id": "MCB-A4-001",
        "lot_id": "LOT-A4-001",
        "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED",
        "adjunct_depletion_mode": "NATURAL_AGING_DECAY",
        "washout_and_natural_decay_unpooled_verified": True,
        "carrier_and_soluble_adjunct_below_lod_verified": True,
        "obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified": True,
        "obs_arm2_sa_thbd_bioactivity_depleted_verified": True,
        "obs_arm2_contemporaneous_sham_handling_matched": True,
        "obs_arm3_ongoing_adjunct_bioactivity_verified": True,
        "arm2_depleted_sample_id": "SMP-ARM2-DEP-01",
        "contemporaneous_arm3_control_sample_id": "SMP-ARM3-POSCTRL-01",
        "obs_parallel_arm3_control_cd4_apoptosis_pct": 88.5,
        "se_parallel_arm3_cd4_apoptosis": 1.8,
        "adjunct_immune_responder_history": "PAIRED_BOTH",
        "adjunct_sa_fhl1_status": "DEPLETED_VERIFIED",
        "n_independent_blocks": 6,
        "recipient_serum_pair_id": "SERUM-PAIR-T1D-001",
        "intact_serum_sample_id": "SERUM-INTACT-001",
        "igg_depleted_serum_sample_id": "SERUM-IGGDEP-001",
        "positive_adcc_control_id": "ADCC-POSCTRL-MAB-01",
        "isotype_adcc_control_id": "ADCC-ISOCTRL-IGG1-01",
        "obs_igg_depletion_efficiency_pct": 98.4,
        "igg_depleted_matched_serum_control_verified": True,
        "fcr_effector_fcgr3a_158v_typed_verified": True,
        "obs_indirect_apc_dq2_dq8_cd4_stimulation_index": 1.18, "se_cd4_si": 0.08,
        "obs_indirect_apc_cd4_ifng_si": 1.42, "se_ifng_si": 0.12,
        "obs_gate0d_transwell_beta_annexin_v_pct": 7.5, "se_annexin_v": 1.2,
        "obs_arm2_post_depletion_nhs_c5b9_mac_pct": 2.4, "se_c5b9": 0.5,
        "obs_fcr_cd16_intact_serum_total_lysis_pct": 3.1, "se_intact_adcc": 0.6,
        "obs_fcr_cd16_igg_depleted_serum_lysis_pct": 0.5,
        "obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct": 2.6, "se_net_adcc": 0.7,
        "obs_adcc_positive_control_lysis_pct": 42.0, "se_pos_adcc": 3.0,
        "obs_adcc_isotype_background_lysis_pct": 1.4, "se_isotype_adcc": 0.4,
        "obs_arm2_post_depletion_viable_beta_recovery_pct": 91.5, "se_beta_rec": 1.8,
        "obs_arm2_post_depletion_perifusion_gsis_si": 3.2, "se_gsis_si": 0.2,
        "obs_post_depletion_input_normalized_cpep_ratio": 0.94, "se_cpep_ratio": 0.03,
    }
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(_base_arm2_28, root=ROOT) == (
        True,
        "CLEARED_ARM2_COMPLETE_DEPLETION_GATE0D_AND_FCR_ADCC_OVERLAY",
    )
    _b2_auth_28 = wc_mod.build_valid_upstream_b2_record_iip_0028("MCB-B2-001", "LOT-B2-001")
    _b2_arm2_28 = {**_base_arm2_28, "product_variant_id": "ISLET1_B2_5LOCUS_BASE", "mcb_id": "MCB-B2-001", "lot_id": "LOT-B2-001"}
    _b2_rep_28 = wc_mod.evaluate_paired_adcc_and_postdepletion_report_iip_0028(
        _b2_arm2_28, root=ROOT, upstream_b2_authorization_record=_b2_auth_28
    )
    assert _b2_rep_28["gate_pass"] is True
    assert _b2_rep_28["status"] == "RESEARCH_ONLY_PAIRED_STUDY_ARM2_CLEARED_B2_EXECUTION_SUPPORTED_FALSE"
    assert _b2_rep_28["execution_supported"] is False and _b2_rep_28["paired_study_research_only"] is True

    # Negative counterexamples including 30% intact vs 29% depleted (1% delta) high background injury:
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {k: v for k, v in _base_arm2_28.items() if k != "product_variant_id"}, root=ROOT
    ) == (False, "HOLD_MISSING_PRODUCT_VARIANT_ID")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {**_base_arm2_28, "product_variant_id": " ISLET1_B2_5LOCUS_BASE "}, root=ROOT
    ) == (False, "HOLD_WHITESPACE_PADDED_PRODUCT_VARIANT_ID")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {**_base_arm2_28, "product_variant_id": "BOGUS"}, root=ROOT
    ) == (False, "HOLD_UNREGISTERED_PRODUCT_VARIANT_ID")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {**_base_arm2_28, "product_variant_id": "ISLET1_B_SCT_C1_CIITA"}, root=ROOT
    ) == (False, "HOLD_UNSUPPORTED_TRACK2_BSCT_REQUIRES_SEPARATE_GATE0B_SCT")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {**_b2_arm2_28, "registered_and_qc_cleared": True}, root=ROOT
    ) == (False, "HOLD_B2_CALLER_BOOLEAN_REJECTED_MISSING_UPSTREAM_REGISTRY_RECORD")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        _b2_arm2_28, root=ROOT, upstream_b2_authorization_record={**_b2_auth_28, "execution_supported": True}
    ) == (False, "HOLD_B2_EXECUTION_SUPPORTED_MUST_REMAIN_FALSE")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        _b2_arm2_28, root=ROOT, upstream_b2_authorization_record={**_b2_auth_28, "bound_c1_c5_pair_sha256": "f" * 64}
    ) == (False, "HOLD_B2_C1_C5_PAIR_DIGEST_MISMATCH")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {
            **_base_arm2_28,
            "obs_fcr_cd16_intact_serum_total_lysis_pct": 30.0,
            "obs_fcr_cd16_igg_depleted_serum_lysis_pct": 29.0,
            "obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct": 1.0,
        },
        root=ROOT,
    ) == (False, "HOLD_EXCEEDED_UPPER95_obs_fcr_cd16_intact_serum_total_lysis_pct")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {**_base_arm2_28, "adjunct_sa_fhl1_status": "RETAINED_ACTIVE"}, root=ROOT
    ) == (False, "HOLD_GATE0D_PARTIALLY_DEPLETED_SA_FHL1_RETAINED_ACTIVE_NOT_COMPLETE_DEPLETION")
    assert wc_mod.evaluate_paired_adcc_and_postdepletion_gate_iip_0028(
        {**_base_arm2_28, "obs_arm3_ongoing_adjunct_bioactivity_verified": False}, root=ROOT
    ) == (False, "HOLD_GATE0D_PARALLEL_ARM3_BIOACTIVITY_CONTROL_UNVERIFIED")
    with pytest.raises(ValueError, match="execution_supported"):
        wc_mod.evaluate_stage0_decision_tree(
            {"product_variant_id": "ISLET1_B2_5LOCUS_BASE"},
            root=ROOT,
            product_def={"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "execution_supported": True},
        )

    # Cross-arm batch pooling block in run_calibration and replay_stage0_chronological_batches:
    _cal_arm3_28 = wetlab.run_calibration(
        {
            "batch_id": "BATCH-IIP28-ARM3",
            "mcb_clone_id": "CLONE-IIP28",
            "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
            "adjunct_conditioning_history_arm": "ARM3_CONTINUOUSLY_ACTIVE",
            "obs_stage7_yield_ratio": 0.92,
            "obs_stage7_yield_se": 0.02,
        },
        is_dry_run=False,
    )
    with pytest.raises(ValueError, match="adjunct_conditioning_history_arm"):
        wetlab.run_calibration(
            {
                "batch_id": "BATCH-IIP28-ARM2",
                "mcb_clone_id": "CLONE-IIP28",
                "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
                "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED",
                "obs_cd47_fold": 4.5,
                "obs_cd47_fold_se": 0.2,
            },
            False,
            _cal_arm3_28,
        )
    assert wetlab.replay_stage0_chronological_batches(
        [
            {"seq_index": 1, "utc_timestamp": "2026-11-18T08:00:00Z", "scope": {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": "L28", "campaign_id": "C28", "adjunct_conditioning_history_arm": "ARM3_CONTINUOUSLY_ACTIVE"}, "measurements": dict(_base_s9_26), "iip0026_molecular_observations": {}},
            {"seq_index": 2, "utc_timestamp": "2026-11-18T09:00:00Z", "scope": {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": "L28", "campaign_id": "C28", "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED"}, "measurements": {}, "iip0026_molecular_observations": {}},
        ],
        root=ROOT,
    )["overall_status"] == "FAIL_CROSS_CONDITIONING_ARM_REPLAY_MIXING"

    # (b) Stage-0 tree integration when require_full_iip0028_postdepletion_adcc_mixture_and_iapp_panel=True:
    _s0_28_valid = {
        **_clean_iip20_valid,
        "product_variant_id": "ISLET1_A4_ARM1B_EMB_A_STAR_4143",
        "device_assembly_id": "DEV-IMPLANT-28",
        "sister_assembly_id": "DEV-SISTER-28",
        "obs_stage0_independent_lot_replicates_n": 12,
        "require_full_iip0028_postdepletion_adcc_mixture_and_iapp_panel": True,
        "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED",
        "adjunct_depletion_mode": "NATURAL_AGING_DECAY",
        "adjunct_sa_fhl1_status": "DEPLETED_VERIFIED",
        "adjunct_immune_responder_history": "PAIRED_BOTH",
        "washout_and_natural_decay_unpooled_verified": True,
        "carrier_and_soluble_adjunct_below_lod_verified": True,
        "obs_arm2_bound_and_supernatant_fasl_bioactivity_depleted_verified": True,
        "obs_arm2_sa_thbd_bioactivity_depleted_verified": True,
        "obs_arm2_contemporaneous_sham_handling_matched": True,
        "obs_arm3_ongoing_adjunct_bioactivity_verified": True,
        "arm2_depleted_sample_id": "SMP-ARM2-DEP-28",
        "contemporaneous_arm3_control_sample_id": "SMP-ARM3-CTRL-28",
        "obs_parallel_arm3_control_cd4_apoptosis_pct": 89.0,
        "obs_parallel_arm3_control_cd4_apoptosis_se": 1.5,
        "obs_indirect_apc_dq2_dq8_cd4_stimulation_index": 1.18,
        "obs_indirect_apc_dq2_dq8_cd4_stimulation_index_se": 0.06,
        "obs_indirect_apc_cd4_ifng_si": 1.42,
        "obs_indirect_apc_cd4_ifng_si_se": 0.08,
        "obs_gate0d_transwell_beta_annexin_v_pct": 7.5,
        "obs_gate0d_transwell_beta_annexin_v_se": 0.9,
        "obs_arm2_post_depletion_nhs_c5b9_mac_pct": 2.2,
        "obs_arm2_post_depletion_nhs_c5b9_mac_se": 0.4,
        "obs_arm2_post_depletion_viable_beta_recovery_pct": 92.0,
        "obs_arm2_post_depletion_viable_beta_recovery_se": 1.2,
        "obs_arm2_post_depletion_perifusion_gsis_si": 3.3,
        "obs_arm2_post_depletion_perifusion_gsis_si_se": 0.15,
        "obs_post_depletion_input_normalized_cpep_ratio": 0.94,
        "obs_post_depletion_input_normalized_cpep_ratio_se": 0.02,
        "dec084_measured_starting_composition_and_fresh_equal_output_verified": True,
        "obs_dec084_starting_potency_ratio_out_b_over_out_a": 1.0,
        "obs_dec084_starting_potency_ratio_out_b_over_out_a_se": 0.02,
        "obs_dec084_ocr_ratio_b_over_a": 0.87,
        "obs_dec084_ocr_ratio_b_over_a_se": 0.01,
        "obs_dec084_cpep_auc_ratio_b_over_a": 0.96,
        "obs_dec084_cpep_auc_ratio_b_over_a_se": 0.02,
        "obs_dec084_gcg_auc_ratio_b_over_a": 0.92,
        "obs_dec084_gcg_auc_ratio_b_over_a_se": 0.02,
        "recipient_serum_pair_id": "SERUM-PAIR-28",
        "intact_serum_sample_id": "SERUM-INTACT-28",
        "igg_depleted_serum_sample_id": "SERUM-DEP-28",
        "positive_adcc_control_id": "ADCC-POS-28",
        "isotype_adcc_control_id": "ADCC-ISO-28",
        "cd16_adcc_igg_depletion_paired_control_verified": True,
        "fcr_effector_fcgr3a_158v_typed_verified": True,
        "obs_igg_depletion_efficiency_pct": 98.5,
        "obs_igg_depletion_efficiency_se": 0.3,
        "obs_fcr_cd16_intact_serum_total_lysis_pct": 2.8,
        "obs_fcr_cd16_intact_serum_total_lysis_se": 0.5,
        "obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_pct": 2.3,
        "obs_fcr_cd16_net_adcc_lysis_vs_igg_depleted_se": 0.5,
        "obs_cd16_adcc_positive_control_lysis_pct": 41.0,
        "obs_cd16_adcc_positive_control_lysis_se": 2.5,
        "obs_adcc_isotype_background_lysis_pct": 1.3,
        "obs_adcc_isotype_background_lysis_se": 0.3,
        "iapp_age_matched_reference_above_loq_verified": True,
        "iapp_sister_visit_lot_and_viable_beta_mass_concordance_verified": True,
        "obs_iapp_72h_immunocapture_oligomer_elisa_fold": 1.04,
        "obs_iapp_72h_immunocapture_oligomer_elisa_se": 0.05,
        "obs_iapp_72h_colocalized_fibril_area_pct": 0.14,
        "obs_iapp_72h_colocalized_fibril_area_se": 0.03,
        "obs_iapp_d28_immunocapture_oligomer_elisa_fold": 1.13,
        "obs_iapp_d28_immunocapture_oligomer_elisa_se": 0.06,
        "obs_iapp_d28_colocalized_fibril_area_pct": 0.27,
        "obs_iapp_d28_colocalized_fibril_area_se": 0.05,
        "obs_iapp_d84_immunocapture_oligomer_elisa_fold": 1.21,
        "obs_iapp_d84_immunocapture_oligomer_elisa_se": 0.07,
        "obs_iapp_d84_colocalized_fibril_area_pct": 0.41,
        "obs_iapp_d84_colocalized_fibril_area_se": 0.06,
    }
    _dt_28_ok = wc_mod.evaluate_stage0_decision_tree(_s0_28_valid, root=ROOT)
    assert _dt_28_ok["lot_release_readiness_tier"] == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE"
    assert _dt_28_ok["iip0028_postdepletion_adcc_mixture_and_iapp_evaluation"]["gate_pass"] is True
    assert wc_mod.evaluate_stage0_decision_tree(
        {k: v for k, v in _s0_28_valid.items() if k != "obs_arm2_post_depletion_viable_beta_recovery_se"}, root=ROOT
    )["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"
    assert wc_mod.evaluate_stage0_decision_tree(
        {**_s0_28_valid, "obs_fcr_cd16_intact_serum_total_lysis_pct": 30.0}, root=ROOT
    )["lot_release_readiness_tier"] == "HOLD_FAILED_OR_INCOMPLETE_GATE"

    # (c) Layer-1 paired McNemar + phi_AB audit & Layer-4 measured-potency DEC-084 mixture audit:
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        _rk_28_doc = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        _mc_28 = _rk_28_doc["odds_composition_and_profile_scope_audit_iip_0003"]["multiplied_conditional_totals"]["iip_0028_paired_mcnemar_and_phi_ab_concordance_audit"]
        assert _mc_28["proposal_id"] == "IIP-0028"
        _nuis_28 = _mc_28["exact_conditional_mcnemar_summary"]["n24_global_nuisance_size_profile"]
        assert abs(_nuis_28["global_worst_case_size_pct"] - 3.6699) < 1e-3
        assert abs(_nuis_28["argmax_total_discordance_q"] - 0.96335) < 1e-3
        assert abs(_nuis_28["size_at_total_discordance_q_1p0_pct"] - 3.1957) < 1e-3

    _l4_28_doc = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    _mix_28 = _l4_28_doc["part_a_336h_spatial_reaction_diffusion_pouch_pde"]["measured_potency_dec084_mixture_and_glucagon_audit_iip_0028"]
    assert _mix_28["contribution_id"] == "IIP-0028" and _mix_28["decision_id"] == "DEC-106"
    assert _mix_28["measured_potency_loading_and_ocr"]["ocr_ratio_0p95_rounded_5dec"] == [0.915, 0.90226]
    assert _mix_28["dual_band_worst_face_3d_fvm_oxygen_sweep"]["boundary_po2_34mmhg_narrow_margin_80um_plus10_stress_mmhg"]["passes_15mmhg_floor_at_34mmhg_fixed_vascular"] is False
    assert _mix_28["rk4_24h_exercise_ode_elapsed_tbr_and_glucagon_sensitivities"]["reference_60_30_10_preserved_paracrine"]["tbr_below_3p9_elapsed_pct"] == 1.1652

    _dec84_in_28 = {
        "selection_blocks_excluded_from_confirmation": True,
        "measured_starting_composition_and_viable_counts_verified": True,
        "n_confirmation_blocks": 8,
        "obs_dec084_starting_potency_ratio_out_b_over_out_a": 1.0,
        "se_starting_potency_ratio": 0.02,
        "implanted_areal_loading_ieq_cm2": 1285.7143,
        "dual_band_worst_face_po2_nx40_mmhg": 15.3699,
        "dual_band_worst_face_po2_nx60_mmhg": 15.6504,
        "obs_dec084_ocr_ratio_b_over_a": 0.8758, "se_ocr_ratio": 0.008,
        "obs_dec084_cpep_auc_ratio_b_over_a": 0.96, "se_cpep_auc_ratio": 0.02,
        "obs_dec084_gcg_auc_ratio_b_over_a": 0.92, "se_gcg_auc_ratio": 0.02,
        "obs_day28_beta_recovery_ratio": 0.88, "se_day28_beta_recovery": 0.02,
        "obs_readout2_perifusion_gsis_si": 3.4, "se_readout2_si": 0.12,
        "obs_readout2_first_phase_peak_to_basal_ratio": 5.2, "se_readout2_peak_basal": 0.2,
        "obs_readout2_first_phase_time_to_peak_min": 5.5, "se_readout2_tpeak": 0.3,
        "obs_readout2_stepdown_shutoff_t50_min": 2.1, "se_readout2_t50": 0.15,
        "obs_gate0m_glucagon_2mM_fold_increase": 2.4, "se_readout2_gcg_fold": 0.12,
    }
    assert wc_mod.evaluate_dec084_measured_potency_mixture_gate_iip_0028(_dec84_in_28)["gate_pass"] is True
    assert wc_mod.evaluate_dec084_measured_potency_mixture_gate_iip_0028(
        {**_dec84_in_28, "obs_dec084_starting_potency_ratio_out_b_over_out_a": 0.95, "implanted_areal_loading_ieq_cm2": 1353.3835, "obs_dec084_ocr_ratio_b_over_a": 0.915}
    )["status"] == "HOLD_DEC084_EXCEEDED_UPPER95_obs_dec084_ocr_ratio_b_over_a"
    assert wc_mod.evaluate_dec084_measured_potency_mixture_gate_iip_0028(
        {**_dec84_in_28, "dual_band_worst_face_po2_nx40_mmhg": 14.7514, "dual_band_worst_face_po2_nx60_mmhg": 15.0465}
    )["status"] == "HOLD_DEC084_NUMERICALLY_MESH_SENSITIVE_WORST_FACE_OXYGEN_BELOW_15MMHG"

    # (d) Visit/lot/sister-bound 72h/d28/d84 dual-modality IAPP + viable beta mass + GY-023:
    _lots_28 = [f"LOT-A4-28-{i:03d}" for i in range(1, 9)]
    _rep_lots_28 = [f"LOT-A4-REP28-{i:03d}" for i in range(1, 7)]

    def _mk_iapp_tp_28(tp: str, ef: float, es: float, ff: float, fs: float, ap: float, aps: float, **kw: object) -> dict:
        rec = {
            "product_variant_id": "ISLET1_A4_ARM1B",
            "lot_id": "LOT-A4-28-001",
            "sample_id": kw.get("sample_id", f"SMP-28-{tp.upper()}-01"),
            "loq_signal": kw.get("loq_signal", 1.0),
            "ko_blank_id": f"KO-28-{tp.upper()}",
            "ko_blank_signal": 0.15,
            "insulin_fibril_blank_signal": 0.20,
            "age_matched_ref_control_id": f"REF-28-{tp.upper()}",
            "ref_signal": kw.get("ref_signal", 12.0),
            "spike_recovery_iapp_1_37_pct": 93.0,
            "viable_beta_mass_pct": kw.get("viable_beta_mass_pct", 90.0),
            "se_viable_beta_mass": 1.5,
            "stimulated_gsis_si": kw.get("stimulated_gsis_si", 3.2),
            "se_gsis_si": 0.15,
            "elisa_oligomer_fold": ef, "se_elisa_fold": es,
            "colocalized_fibril_fold": ff, "se_fibril_fold": fs,
            "colocalized_amyloid_area_pct": ap, "se_area_pct": aps,
        }
        if "within_lot_dual_modality_coexceedance_verified" in kw:
            rec["within_lot_dual_modality_coexceedance_verified"] = kw["within_lot_dual_modality_coexceedance_verified"]
        if tp == "72h":
            rec["scope"] = "LOT_SCOPE"
            rec["timepoint_h"] = 72
        elif tp == "d28":
            rec["scope"] = "LOT_SCOPE"
            rec["timepoint_d"] = 28
        else:
            rec["scope"] = "TERMINAL_SISTER_DEVICE_SCOPE"
            rec["timepoint_d"] = 84
            rec["sister_assembly_id"] = kw.get("sister_id", "DEV-SISTER-28-01")
        return rec

    _base_iapp_28 = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "primary_lot_id": "LOT-A4-28-001",
        "contributing_lot_ids": _lots_28,
        "n_independent_blocks": 8,
        "implant_device_assembly_id": "DEV-IMPLANT-28-01",
        "authorized_sister_assembly_ids": ["DEV-SISTER-28-01", "DEV-SISTER-28-02"],
        "timepoint_records": {
            "72h": _mk_iapp_tp_28("72h", 1.05, 0.06, 1.02, 0.05, 0.15, 0.04),
            "d28": _mk_iapp_tp_28("d28", 1.14, 0.08, 1.12, 0.07, 0.28, 0.06),
            "d84": _mk_iapp_tp_28("d84", 1.22, 0.09, 1.19, 0.08, 0.42, 0.08),
        },
    }
    _ia_ok_28 = wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(_base_iapp_28, root=ROOT)
    assert _ia_ok_28["core_cleared"] is True and _ia_ok_28["p2_itm2b_reserve_eligible"] is False and _ia_ok_28["p2_active_weight_pct"] == 0.0
    # Reused raw measurement vector under a fresh sample_id must fail closed:
    assert wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(
        {
            **_base_iapp_28,
            "timepoint_records": {
                **_base_iapp_28["timepoint_records"],
                "d28": _mk_iapp_tp_28("d28", 1.05, 0.06, 1.02, 0.05, 0.15, 0.04, sample_id="SMP-FRESH-REUSED-RAW"),
            },
        },
        root=ROOT,
    )["disposition"] == "HOLD_DUPLICATE_OR_WRONG_TIMEPOINT_SAMPLE"
    assert wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(
        {
            **_base_iapp_28,
            "timepoint_records": {
                **_base_iapp_28["timepoint_records"],
                "d84": _mk_iapp_tp_28("d84", 1.22, 0.09, 1.19, 0.08, 0.42, 0.08, sister_id="DEV-IMPLANT-28-01"),
            },
        },
        root=ROOT,
    )["disposition"] == "HOLD_IAPP_D84_REQUIRES_AUTHORIZED_SISTER_ASSEMBLY"
    # Beta-cell attrition confounding low amyloid signal must fail closed:
    assert wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(
        {
            **_base_iapp_28,
            "timepoint_records": {
                **_base_iapp_28["timepoint_records"],
                "d84": _mk_iapp_tp_28("d84", 1.08, 0.06, 1.05, 0.05, 0.20, 0.04, viable_beta_mass_pct=62.0),
            },
        },
        root=ROOT,
    )["disposition"] == "HOLD_IAPP_CONFOUNDED_BY_VIABLE_BETA_MASS_OR_SECRETION_LOSS"
    # Disjoint-lot-only dual-modality breaches without within-lot co-exceedance require replication:
    assert wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(
        {
            **_base_iapp_28,
            "timepoint_records": {
                **_base_iapp_28["timepoint_records"],
                "d84": _mk_iapp_tp_28("d84", 1.85, 0.08, 1.72, 0.07, 1.35, 0.08, within_lot_dual_modality_coexceedance_verified=False),
            },
        },
        root=ROOT,
    )["disposition"] == "HOLD_DISCORDANT_CHRONIC_NEEDS_REPLICATION"
    _ia_conc_28 = wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(
        {
            **_base_iapp_28,
            "timepoint_records": {
                **_base_iapp_28["timepoint_records"],
                "d84": _mk_iapp_tp_28("d84", 1.85, 0.08, 1.72, 0.07, 1.35, 0.08),
            },
        },
        root=ROOT,
    )
    assert _ia_conc_28["disposition"] == "CONFIRMED_CHRONIC_P2_ELIGIBLE" and _ia_conc_28["p2_active_weight_pct"] == 0.0
    _ia_rep_28 = wc_mod.evaluate_orthogonal_iapp_3way_gate_iip_0028(
        {
            **_base_iapp_28,
            "timepoint_records": {
                **_base_iapp_28["timepoint_records"],
                "d84": _mk_iapp_tp_28("d84", 1.85, 0.08, 1.19, 0.08, 0.42, 0.08),
            },
            "independent_repeat_confirmation": {
                "preregistered_replication_protocol_id": "PROTO-REP-28",
                "repeat_lot_ids": _rep_lots_28,
                "repeat_sample_id": "SMP-REP-D84-28",
                "repeat_sister_assembly_id": "DEV-SISTER-REP-28",
                "timepoint": "d84",
                "timepoint_d": 84,
                "loq_signal": 1.0,
                "ko_blank_id": "KO-REP-28",
                "ko_blank_signal": 0.14,
                "insulin_fibril_blank_signal": 0.19,
                "age_matched_ref_control_id": "REF-REP-28",
                "ref_signal": 12.2,
                "spike_recovery_iapp_1_37_pct": 91.0,
                "elisa_oligomer_fold": 1.84, "se_elisa_fold": 0.08,
            },
        },
        root=ROOT,
    )
    assert _ia_rep_28["disposition"] == "CONFIRMED_CHRONIC_P2_ELIGIBLE" and _ia_rep_28["p2_active_weight_pct"] == 0.0

    # (e) Graveyard GY-023 (10 combinatorial grid + >=13 off-grid; 15 after IIP-0031 GY-025):
    _gy_28 = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    _og_28 = _gy_28["off_grid_molecular_dead_end_guardrails"]["guardrails"]
    assert len(_gy_28["graveyard_entries"]) == 10 and len(_og_28) == 15
    _gy23 = [g for g in _og_28 if g["graveyard_card_id"] == "GY-023"]
    assert len(_gy23) == 1 and _gy23[0]["guardrail_id"].startswith("OFFGRID_GUARD_013_") and _gy23[0]["contribution_id"] == "IIP-0028"


def test_53_iip0027_promoter_origin_tss_6pwm_ledger_and_c1_c2_race_hold() -> None:
    """
    Verify IIP-0027 / DEC-107 across:
      (a) Layer-2 part-resolved 6-PWM JASPAR annotation ledger across canonical constructs,
          DX-L3, L3-BF-B (including Unit-1 53-bp repeat attribution), DEC-105 c1..c6,
          IIP-0018 ARM1, and 145-bp AAV2 ITR (0/0/0), rejecting motif-count closure of OQ-011;
      (b) Tightened omitted/null/UNRESOLVED c1/c2 PD-L1 5'-RACE HOLD_PENDING_5PRIME_RACE_TSS
          parity across Direct, Serialized JSON Ingestion, and Chronological Replay paths;
      (c) Fail-closed Stage-0 resting (0 h) + interferon-stimulated (24 h) CXCL10 promoter-of-origin
          5'-RACE TSS attribution and conditional BIO-06 AAV-ITR RNA-origin attribution;
      (d) Contract & decision-tree sync and Layer-1 2,160-stack invariance.
    """
    ROOT = Path(__file__).resolve().parent.parent
    from scripts import ingest_wetlab_results as wetlab
    from scripts import wetlab_contract as wc_mod

    # (a) Layer-2 6-PWM ledger + 53-bp repeat resolution + 145-bp AAV2 ITR negative sensitivity:
    _l2_27 = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    _aud_27 = _l2_27["stage0_locus3_initiation_6pwm_and_aav_itr_audit_iip_0027"]
    assert _aud_27["decision_id"] == "DEC-107" and _aud_27["contribution_id"] == "IIP-0027"
    assert _aud_27["motif_count_stopping_rule_status"] == "REJECTED_DESCRIPTIVE_ANNOTATION_ONLY"
    assert _aud_27["oq011_status"] == "OPEN_REQUIRES_RESTING_AND_STIMULATED_PROMOTER_OF_ORIGIN_5PRIME_RACE"
    _scans_27 = _aud_27["scanned_constructs_and_itr"]
    _expected_triplets = {
        "L1_ARM1A_4880": [3, 3, 3],
        "L1_ARM1B_3974": [12, 4, 6],
        "L3_CANONICAL_4090": [4, 4, 2],
        "L4_3457": [1, 1, 0],
        "L3_SCARFREE_3799": [10, 3, 5],
        "L3_BFB_BASALFLOOR_4078": [12, 5, 5],
        "DEC105_c1": [12, 5, 5],
        "DEC105_c2": [10, 3, 5],
        "DEC105_c3": [3, 2, 2],
        "DEC105_c4": [1, 1, 0],
        "DEC105_c5": [0, 0, 0],
        "DEC105_c6": [1, 1, 0],
        "IIP0018_ARM1_L3_WT_PARENT_4078": [12, 5, 5],
        "AAV2_ITR_145BP_J01901_1_145": [0, 0, 0],
    }
    for _cid_27, _exp_t in _expected_triplets.items():
        for _pc_lbl in ("pseudocount_0p8_production", "pseudocount_0p5_comparator"):
            assert _scans_27[_cid_27]["jaspar_6pwm_by_pseudocount"][_pc_lbl]["counts_window_merged_region6_total"] == _exp_t
    for _k2a in (
        "L3_SCARFREE_3799", "L3_BFB_BASALFLOOR_4078",
        "DEC105_c1", "DEC105_c2", "DEC105_c3", "DEC105_c4", "DEC105_c5", "DEC105_c6",
        "IIP0018_ARM1_L3_WT_PARENT_4078",
    ):
        assert _scans_27[_k2a]["jaspar_6pwm_by_pseudocount"]["pseudocount_0p8_production"]["non_promoter_cargo_significant_window_hits"] == 0
    _rep53 = _aud_27["summary_findings"]["l3_bfb_and_c1_53bp_direct_repeat_resolution"]
    assert (_rep53["repeat_length_bp"], _rep53["pos1_a_1based"], _rep53["pos1_b_1based"]) == (53, 650, 680)
    assert _rep53["wholly_inside_declared_unit1_promoter"] is True and _rep53["identical_in_dec105_c1_650_879"] is True

    # (b) Tightened omitted/null/UNRESOLVED c1/c2 PD-L1 5'-RACE hold across all 3 paths:
    _base_b2_27 = {
        "obs_icasp9_log10_kill": 4.6,
        "obs_icasp9_log10_kill_se": 0.08,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000.0,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000.0,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300000000,
        "product_variant_id": "ISLET1_B2_5LOCUS_BASE",
        "lot_id": "LOT_B2_IIP0027",
        "campaign_id": "CAMP_S0_IIP0027",
        "bound_c1_c5_pair_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
        "construct_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
        "donor_seedline_id": "SEED_IIP0027_01",
        "mcb_id": "MCB_IIP0027_01",
        "tier6b_sorted_lot_observations": [
            {
                "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
                "scope_type": "LOT_SCOPE",
                "kind": "EXACT_INTEGER_COUNT",
                "timepoint": "STAGE7_PRE_ASSEMBLY",
                "assay_version": "TIER6B_SORTED_LOT_V1",
                "unit": "clonogenic_opportunities",
                "product_variant_id": "ISLET1_B2_5LOCUS_BASE",
                "construct_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
                "donor_seedline_id": "SEED_IIP0027_01",
                "mcb_id": "MCB_IIP0027_01",
                "lot_id": "LOT_B2_IIP0027",
                "campaign_id": "CAMP_S0_IIP0027",
                "sample_id": "SMP_27_ENR",
                "fraction_id": "enriched",
                "declared_fraction_roster": ["enriched", "flow_through"],
                "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
                "clonogenic_risk_unit": "CD31_NEG_CD144_NEG_CLONOGENIC_OPPORTUNITY",
                "endogenous_denominator_verified": True,
                "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_27_ENR", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
                "presort_mass_balance_reconciled": True,
                "unreconciled_loss_fraction": 0.0,
                "presort_total_input_cells": 3_000_000_000,
                "presort_reconciled_loss_cells": 0,
                "presort_fraction_cells": 2_999_700_000,
                "plated_treated_cells_n": 80_000_000,
                "w_i": 0.9999,
                "N": 30_000_000,
                "k": 0,
                "hpsc_calibrant": {"r_k": 20, "r_n": 20},
                "edited_differentiated_calibrant": {"r_k": 20, "r_n": 20},
            },
            {
                "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
                "scope_type": "LOT_SCOPE",
                "kind": "EXACT_INTEGER_COUNT",
                "timepoint": "STAGE7_PRE_ASSEMBLY",
                "assay_version": "TIER6B_SORTED_LOT_V1",
                "unit": "clonogenic_opportunities",
                "product_variant_id": "ISLET1_B2_5LOCUS_BASE",
                "construct_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
                "donor_seedline_id": "SEED_IIP0027_01",
                "mcb_id": "MCB_IIP0027_01",
                "lot_id": "LOT_B2_IIP0027",
                "campaign_id": "CAMP_S0_IIP0027",
                "sample_id": "SMP_27_FT",
                "fraction_id": "flow_through",
                "declared_fraction_roster": ["enriched", "flow_through"],
                "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
                "clonogenic_risk_unit": "CD31_NEG_CD144_NEG_CLONOGENIC_OPPORTUNITY",
                "endogenous_denominator_verified": True,
                "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_27_FT", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
                "presort_mass_balance_reconciled": True,
                "unreconciled_loss_fraction": 0.0,
                "presort_total_input_cells": 3_000_000_000,
                "presort_reconciled_loss_cells": 0,
                "presort_fraction_cells": 300_000,
                "plated_treated_cells_n": 260_000,
                "w_i": 0.0001,
                "N": 100_000,
                "k": 0,
                "hpsc_calibrant": {"r_k": 10, "r_n": 10},
                "edited_differentiated_calibrant": {"r_k": 10, "r_n": 10},
            },
        ],
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_per_locus_segmented_baf": {"chr19q": 0.502, "chr11p15": 0.499, "chr6p21": 0.504},
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "design_b_hlac_phased_guide_cfd_max": 0.0188,
        "design_b_hla_ab_separated_pass_delivery_verified": True,
        "design_b_chr6p21_1p41mb_loh_sv_cleared": True,
        "design_b_hlac_dt9_positive": True,
        "design_b_hlae_3d12_positive": True,
        "design_b_hlac_dt9_matched_hlae_only_control_delta_ratio": 4.25,
        "design_b_hlac_allele_specific_surface_or_phased_mrna_verified": True,
        "design_b_locus5_5prime_race_coord_bp": 880,
        "design_b_locus5_3prime_race_coord_bp": 3478,
        "design_b_locus5_pdl1_5prime_race_coord_bp": 2350,
        "design_b_locus5_race_umi_n_eff": 3689,
        "design_b_locus5_aberrant_umi_k": 0,
    }

    def _parity_3path_27(case_p: dict) -> dict:
        d_out = wc_mod.evaluate_stage0_research_and_molecular_guardrails(dict(case_p), root=ROOT)
        sc_27 = {k: case_p.get(k) for k in ("product_variant_id", "lot_id", "campaign_id")}
        m_27 = {k: v for k, v in case_p.items() if str(k).startswith("obs_") and v is not None}
        mol_27 = {k: v for k, v in case_p.items() if k not in sc_27 and not str(k).startswith("obs_")}
        i_out = wetlab.ingest_stage0_research_or_wetlab_batch(
            json.dumps({"scope": sc_27, "measurements": m_27, "iip0026_molecular_observations": mol_27}),
            root=ROOT,
        )
        r_out = wetlab.replay_stage0_chronological_batches(
            [
                {"seq_index": 1, "utc_timestamp": "2026-11-18T12:00:00Z", "scope": dict(sc_27), "measurements": dict(m_27), "iip0026_molecular_observations": {}},
                {"seq_index": 2, "utc_timestamp": "2026-11-18T13:00:00Z", "scope": dict(sc_27), "measurements": {}, "iip0026_molecular_observations": dict(mol_27)},
            ],
            root=ROOT,
        )
        assert d_out["overall_status"] == i_out["overall_status"] == r_out["overall_status"]
        assert d_out["iip0026_reasons"] == i_out["iip0026_reasons"] == r_out["iip0026_reasons"]
        return d_out

    assert _parity_3path_27(_base_b2_27)["overall_status"] == "RESEARCH_COMPARATOR_STAGE0_PANEL_CLEARED_PENDING_STATE4_PROMOTION"
    _omitted_27 = {k: v for k, v in _base_b2_27.items() if k != "design_b_locus5_pdl1_5prime_race_coord_bp"}
    assert _parity_3path_27(_omitted_27)["overall_status"] == "HOLD_PENDING_5PRIME_RACE_TSS"
    assert _parity_3path_27({**_base_b2_27, "design_b_locus5_pdl1_5prime_race_coord_bp": None})["overall_status"] == "HOLD_PENDING_5PRIME_RACE_TSS"
    assert _parity_3path_27({**_base_b2_27, "design_b_locus5_pdl1_5prime_race_coord_bp": "UNRESOLVED"})["overall_status"] == "HOLD_PENDING_5PRIME_RACE_TSS"
    assert _parity_3path_27({**_base_b2_27, "design_b_locus5_pdl1_5prime_race_coord_bp": 2300})["overall_status"] == "FAIL_MATURE_TRANSCRIPT_INCLUDES_UPSTREAM_OUT_OF_FRAME_KOZAK_ORF"

    # (c) Stage-0 resting + stimulated CXCL10 promoter-of-origin 5'-RACE & conditional BIO-06 ITR gate:
    _clean_27 = {
        "obs_icasp9_log10_kill": 4.6,
        "obs_icasp9_log10_kill_se": 0.08,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "require_full_iip0027_promoter_origin_and_itr_panel": True,
        "cxcl10_tss_donor_id": "L3_BFB_BASALFLOOR_4078",
        "cxcl10_tss_transcription_unit_id": "UNIT2_PDL1",
        "cxcl10_insert_matched_comparator_control_verified": True,
        "obs_cxcl10_resting_race_umi_n_eff": 3689,
        "obs_cxcl10_resting_unexpected_tss_umi_k": 0,
        "obs_cxcl10_resting_unresolved_origin_umi_k": 0,
        "obs_cxcl10_ifng24h_race_umi_n_eff": 3689,
        "obs_cxcl10_ifng24h_unexpected_tss_umi_k": 0,
        "obs_cxcl10_ifng24h_unresolved_origin_umi_k": 0,
        "obs_bio06_retained_aav_itr_dna_pct": 0.2,
        "bio06_itr_rna_origin_attribution_verified": True,
        "obs_bio06_itr_origin_rna_umi_k": 0,
    }
    _res_clean_27 = wc_mod.evaluate_stage0_decision_tree(dict(_clean_27), root=ROOT)
    assert _res_clean_27["iip0027_promoter_origin_and_itr_evaluation"]["status"] == "CLEARED_IIP0027_RESTING_AND_STIMULATED_PROMOTER_ORIGIN_AND_ITR_GATES"
    assert _res_clean_27["iip0027_promoter_origin_and_itr_evaluation"]["oq011_closed_by_computational_motif_count"] is False

    # Resting below LoQ (N=0) allowed ONLY when confirmed AND stimulated 24 h still clears N >= 3,689:
    _res_loq_27 = wc_mod.evaluate_stage0_decision_tree(
        {**_clean_27, "obs_cxcl10_resting_race_umi_n_eff": 0, "cxcl10_resting_expression_below_loq_confirmed": True},
        root=ROOT,
    )
    assert _res_loq_27["iip0027_promoter_origin_and_itr_evaluation"]["gate_pass"] is True

    # Resting expressed with N = 3688 < 3689 holds; resting-only holds; stimulated unexpected TSS holds:
    assert wc_mod.evaluate_stage0_decision_tree({**_clean_27, "obs_cxcl10_resting_race_umi_n_eff": 3688}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_IIP0027_RESTING_UMI_BELOW_3689"
    assert wc_mod.evaluate_stage0_decision_tree({k: v for k, v in _clean_27.items() if "ifng24h" not in k}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_IIP0027_MISSING_IFNG24H_STIMULATED_PROMOTER_ORIGIN_5PRIME_RACE"
    assert wc_mod.evaluate_stage0_decision_tree({k: v for k, v in _clean_27.items() if "resting_race" not in k}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_IIP0027_MISSING_RESTING_0H_PROMOTER_ORIGIN_5PRIME_RACE"
    assert wc_mod.evaluate_stage0_decision_tree({**_clean_27, "cxcl10_insert_matched_comparator_control_verified": False}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_IIP0027_MISSING_INSERT_MATCHED_COMPARATOR_CONTROL"
    assert wc_mod.evaluate_stage0_decision_tree({**_clean_27, "obs_cxcl10_ifng24h_unexpected_tss_umi_k": 1}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_FAIL_IIP0027_IFNG24H_UNEXPECTED_INTERNAL_TSS_DETECTED"
    assert wc_mod.evaluate_stage0_decision_tree({**_clean_27, "obs_cxcl10_ifng24h_unresolved_origin_umi_k": 2}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_IIP0027_IFNG24H_UNRESOLVED_TSS_ORIGIN"

    # Conditional BIO-06 ITR RNA-origin checks:
    _no_itr_rna = {k: v for k, v in _clean_27.items() if k not in ("bio06_itr_rna_origin_attribution_verified", "obs_bio06_itr_origin_rna_umi_k")}
    assert wc_mod.evaluate_stage0_decision_tree(_no_itr_rna, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_IIP0027_RETAINED_AAV_ITR_DNA_REQUIRES_RNA_ORIGIN_ATTRIBUTION"
    assert wc_mod.evaluate_stage0_decision_tree({**_clean_27, "obs_bio06_itr_origin_rna_umi_k": 1}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_FAIL_IIP0027_AAV_ITR_ORIGIN_TRANSCRIPT_DETECTED"
    assert wc_mod.evaluate_stage0_decision_tree({**_clean_27, "obs_bio06_retained_aav_itr_dna_pct": 0.5}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_FAIL_BIO06_AAV_ITR_DNA_CAPTURE_GE_0P5PCT"

    # (d) Contract, decision-tree & Layer-1 2,160-stack invariance:
    _ct_27 = json.loads((ROOT / "data/protocols/wetlab_endpoint_contract_v0.8.0.json").read_text(encoding="utf-8"))
    assert len(_ct_27["stage0_decision_tree_endpoints"]) == 366
    assert len(_ct_27["metadata_keys"]) == 77
    assert len(_ct_27["boolean_policy_keys"]) == 122
    _dt_27 = json.loads((ROOT / "data/protocols/stage0_decision_tree_v1.json").read_text(encoding="utf-8"))
    assert _dt_27["promoter_origin_tss_and_conditional_aav_itr_iip_0027"]["decision_id"] == "DEC-107"
    if (ROOT / "data/extracted/islet1_candidate_stack_rankings.json").exists():
        _rk_27 = json.loads((ROOT / "data/extracted/islet1_candidate_stack_rankings.json").read_text(encoding="utf-8"))
        assert _rk_27["total_valid_combinatorial_stacks_evaluated"] == 2160


def test_54_iip0029_alk5i_free_maturation_twoclock_gates_and_mode_b_v2() -> None:
    """
    Verify IIP-0029 / DEC-108 across:
      (a) Cross-file Rule S8 ALK5i-free reconciliation (Layer 4 pivot_rule_S8, decision tree,
          Playbook Row 7 S8 / Row 8 S8b, Master Spec, HIP-Core protocol, and Blueprint);
      (b) Source-bound two-clock/two-pulse empirical qualification subgate with ALL independent
          review hardenings (absolute + contemporaneous adult-islet reference second-pulse floors,
          separate Readout 2 vs Gate 0M Basal1/Peak1 <= 0.15 + full Gate 0M paracrine checks,
          locked UTC timestamps, independent differentiation blocks with matched 40 + 100 mmHg
          coverage, per-trace tracer alignment & censored short low-interval handling, common-
          denominator sorting+reaggregation recovery, registered-mixture + measured-composition
          DEC-084 linkage, pre-iEC LC-MS/MS washout + vascular function verification, exploratory
          Keymeulen proinsulin:C-peptide labeling, and zero standalone clinical release bypass);
      (c) Identical non-overridable >138.3 IU Gate 0S escalation (U95_t or max_lot), wrong-lot
          rejection, and immutable source-replacement latching across Flat, Direct, Serialized
          JSON Ingestion, Chronological Replay, and Atomic Stage-0 Ledger paths;
      (d) Opt-in exploratory Mode B v2 granule-store sensitivity comparator (domain guards,
          exact Mode A default parity, multi-state SciPy DOP853 verification, positivity/mass
          conservation, censored short-interval shutoff, hepatic sweep, and ceteris-paribus
          isolated tau_refill meal sensitivity);
      (e) Idea Graveyard GY-024 (OFFGRID_GUARD_014_*, 10 grid + 14 off-grid = 24 total) and
          updated contract counts (354 endpoints / 72 metadata keys / 110 boolean policy keys).
    """
    import importlib.util
    import json
    from pathlib import Path

    ROOT = Path(__file__).resolve().parents[1]
    from scripts import ingest_wetlab_results as ing_mod
    from scripts import wetlab_contract as wc_mod

    def _load_mod(name: str, path: Path):
        spec = importlib.util.spec_from_file_location(name, path)
        assert spec and spec.loader
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    l4_mod = _load_mod("l4_iip0029_test", ROOT / "scripts" / "run_layer4_tissue_and_metabolic_twin.py")

    # (a) Cross-file Rule S8 ALK5i-free reconciliation:
    l4_atlas = json.loads((ROOT / "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    corner = l4_atlas["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["sc_islet_immaturity_paracrine_corner_sweep_muse_r3"]["stage0_gate_0m_maturity_pass_criteria"]
    assert "pivot_rule_S6" not in corner and "pivot_rule_S8" in corner
    assert "ALK5i-free" in corner["pivot_rule_S8"] and "zero terminal ALK5i II" in corner["pivot_rule_S8"]
    assert "+7-14 days thyroid/TGFb-inhibitor" not in corner["pivot_rule_S8"]

    dt_doc = json.loads((ROOT / "data/protocols/stage0_decision_tree_v1.json").read_text(encoding="utf-8"))
    s8_row = next(s for s in dt_doc["switches"] if s["switch_id"] == "S8")
    assert "ALK5i-free" in s8_row["action"] and "PMID:30661993" in s8_row["action"]
    assert dt_doc["twoclock_alk5i_free_maturation_and_repeated_challenge_iip_0029"]["decision_id"] == "DEC-108"

    pb_txt = (ROOT / "docs/STAGE-0-EXECUTION-PLAYBOOK.md").read_text(encoding="utf-8")
    assert "| **7** | **`S8`** |" in pb_txt and "| **8** | **`S8b`** |" in pb_txt
    assert "#### 4.1.16" in pb_txt and "DEC-108" in pb_txt
    bp_txt = (ROOT / "ISLET-1-BLUEPRINT.md").read_text(encoding="utf-8")
    assert "IIP-0029" in bp_txt and "DEC-108" in bp_txt
    ms_txt = (ROOT / "docs/ISLET-1-MASTER-SPEC-v0.1.md").read_text(encoding="utf-8")
    assert "`ALK5i_II` withdrawal" in ms_txt

    # (b) Source-bound two-clock / two-pulse empirical qualification subgate + ALL 4 review hardenings:
    valid_sha = "f74728767987c18bdcb8c9d49773ee8dbe10e6006e57aa622638e4fb07062b0b"

    def _make_record(
        pvid: str = "ISLET1_A4_ARM1B",
        csha: str = valid_sha,
        lot_id: str = "LOT_S7_2026_01",
        n_lots: int = 6,
        reg_mix: str = "REF_60_30_10",
        beta_f: float = 0.60,
        alpha_f: float = 0.30,
        delta_f: float = 0.10,
        imp_ieq: float = 500_000.0,
        proc_mode: str = "SORTED_REAGGREGATED",
        fluidics: str = "DYNAMIC_PERIFUSION",
        **overrides,
    ):
        lot_ids = [f"CONF_{i+1:02d}" for i in range(n_lots)]
        block_ids = [f"DIFF_BLK_{i+1:02d}" for i in range(n_lots)]
        base_eta = [1.0] * n_lots if proc_mode == "UNSORTED_CONTROL" else [0.64, 0.68, 0.61, 0.66, 0.63, 0.67, 0.65, 0.62][:n_lots]
        base_basal1 = [0.40, 0.43, 0.39, 0.42, 0.41, 0.39, 0.42, 0.40][:n_lots]
        base_low1 = [0.42, 0.45, 0.40, 0.44, 0.43, 0.41, 0.44, 0.42][:n_lots]
        base_high1 = [2.25, 2.35, 2.15, 2.30, 2.28, 2.20, 2.32, 2.24][:n_lots]
        base_peak1 = [3.80, 4.05, 3.65, 3.90, 3.85, 3.75, 3.95, 3.82][:n_lots]
        base_tpeak1 = [3.2, 3.5, 3.0, 3.4, 3.3, 3.1, 3.4, 3.2][:n_lots]
        base_low2 = [0.44, 0.47, 0.42, 0.46, 0.45, 0.43, 0.46, 0.44][:n_lots]
        base_high2 = [1.85, 1.94, 1.78, 1.90, 1.84, 1.81, 1.92, 1.83][:n_lots]
        base_kcl = [4.90, 5.10, 4.75, 5.05, 4.95, 4.85, 5.00, 4.92][:n_lots]
        base_adult_ref = [1.20, 1.26, 1.15, 1.24, 1.19, 1.18, 1.25, 1.21][:n_lots]
        base_t50 = [1.05, 1.15, 0.98, 1.12, 1.08, 1.02, 1.10, 1.04][:n_lots]
        base_t90 = [2.35, 2.48, 2.28, 2.42, 2.38, 2.31, 2.44, 2.34][:n_lots]
        base_tau = [0.32, 0.36, 0.30, 0.35, 0.34, 0.33, 0.35, 0.31][:n_lots]
        base_store = [86.0, 92.0, 84.0, 90.0, 88.0, 91.0, 89.0, 87.0][:n_lots]
        base_proins = [0.042, 0.048, 0.039, 0.045, 0.044, 0.041, 0.046, 0.043][:n_lots]
        base_gcg = [3.1, 3.3, 2.9, 3.2, 3.0, 3.15, 3.25, 3.05][:n_lots]
        base_thr = [5.1, 5.3, 4.9, 5.2, 5.0, 5.15, 5.25, 5.05][:n_lots]
        base_ucn3 = [0.86, 0.89, 0.84, 0.88, 0.85, 0.87, 0.88, 0.86][:n_lots]
        base_creg = [0.92, 0.95, 0.90, 0.94, 0.91, 0.93, 0.94, 0.92][:n_lots]

        per_lot_ov = overrides.pop("_per_lot_overrides", {})
        lot_obs = []
        for i, lid in enumerate(lot_ids):
            it = {
                "lot_scope_key": {
                    "product_variant_id": pvid,
                    "construct_sha256": csha,
                    "donor_seedline_id": "DONOR1_H1_GMP",
                    "mcb_id": "MCB_D1_001",
                    "lot_id": lid,
                },
                "differentiation_block_id": block_ids[i],
                "trace_id": f"TR_40_{lid}",
                "assay_po2_mmhg": 40.0,
                "eta_sort_viable_recovery": base_eta[i],
                "basal1_pre_pulse": base_basal1[i],
                "low1_mean": base_low1[i],
                "high1_mean": base_high1[i],
                "peak1_instantaneous": base_peak1[i],
                "time_to_peak1_min": base_tpeak1[i],
                "low2_mean": base_low2[i],
                "high2_mean": base_high2[i],
                "low2_washout_duration_min": 15.0,
                "kcl_mean": base_kcl[i],
                "adult_islet_ref_z2": base_adult_ref[i],
                "effluent_t50_shutoff_min": base_t50[i],
                "effluent_t90_shutoff_min": base_t90[i],
                "inert_tracer_tau_chamber_min": base_tau[i],
                "tracer_response_valid": True,
                "tracer_time_aligned": True,
                "glucagon_2mM_fold_increase": base_gcg[i],
                "glucose_threshold_mM": base_thr[i],
                "ucn3_maturity_ratio": base_ucn3[i],
                "alpha_counterregulatory_gain_ratio": base_creg[i],
                "matched_100mmhg_trace_id": f"TR_100_{lid}",
                "matched_100mmhg_assay_po2_mmhg": 100.0,
                "matched_100mmhg_low1_mean": base_low1[i],
                "matched_100mmhg_high1_mean": base_high1[i] * 1.05,
                "matched_100mmhg_low2_mean": base_low2[i],
                "matched_100mmhg_high2_mean": base_high2[i] * 1.05,
                "matched_100mmhg_tracer_tau_chamber_min": base_tau[i],
                "matched_100mmhg_tracer_response_valid": True,
                "matched_100mmhg_tracer_time_aligned": True,
                "stored_mature_insulin_iu_implanted": base_store[i],
                "proinsulin_to_cpeptide_content_molar_ratio": base_proins[i],
            }
            for ko, vo in per_lot_ov.items():
                it[ko] = vo[i] if isinstance(vo, list) else vo
            lot_obs.append(it)

        rec = {
            "lot_scope_key": {
                "product_variant_id": pvid,
                "construct_sha256": csha,
                "donor_seedline_id": "DONOR1_H1_GMP",
                "mcb_id": "MCB_D1_001",
                "lot_id": lot_id,
            },
            "n_independent_lots": n_lots,
            "locked_analysis_plan_id": "SAP_IIP0029_V1",
            "pilot_selection_lot_ids": ["PILOT_01", "PILOT_02", "PILOT_03"],
            "pilot_differentiation_block_ids": ["PILOT_BLK_01", "PILOT_BLK_02", "PILOT_BLK_03"],
            "confirmation_lot_ids": lot_ids,
            "window_selection_design": "PILOT_SCREEN_THEN_DISJOINT_CONFIRMATION",
            "candidate_windows_evaluated_k": 1,
            "preparation_process_mode": proc_mode,
            "sorting_performed": (proc_mode == "SORTED_REAGGREGATED"),
            "stage7_start_utc": "2026-10-01T00:00:00Z",
            "sort_reaggregation_utc": "2026-10-15T00:00:00Z",
            "assay_harvest_utc": "2026-10-22T00:00:00Z" if proc_mode == "SORTED_REAGGREGATED" else "2026-10-15T00:00:00Z",
            "stage7_total_days": 21.0 if proc_mode == "SORTED_REAGGREGATED" else 14.0,
            "sort_reaggregation_stage7_day": 14.0,
            "days_post_sort_reaggregation": 7.0 if proc_mode == "SORTED_REAGGREGATED" else 0.0,
            "obs_alk5i_free_maturation_applied": True,
            "terminal_alk5i_tgfb_inhibitor_continued": False,
            "stage7_t3_10nM_and_nac_1mM_retained": True,
            "post_sort_endocrine_zm447439_regimen": "CONTINUOUS_0P5UM_BALBOA_S7",
            "post_sort_zm447439_uM": 0.50,
            "vascular_iec_coassembled": True,
            "coassembly_zm447439_uM": 0.001,
            "coassembly_zm447439_lloq_uM": 0.005,
            "zm447439_washout_assay_method": "LC-MS/MS",
            "zm447439_wash_cycles": 4,
            "zm447439_matrix_spike_recovery_pct": 98.5,
            "zm447439_washout_verified_pre_iec": True,
            "pre_coassembly_iec_vascular_function_verified": True,
            "assay_fluidics_mode": fluidics,
            "hormone_unit": "pmol_cpep_per_1e6_presort_cells_per_min",
            "normalization_denominator": "1e6_presort_viable_endocrine_cells",
            "recovery_denominator_basis": "1e6_presort_viable_endocrine_cells",
            "eta_sort_includes_sorting_and_reaggregation_losses": True,
            "within_lot_technical_var_z2": 0.0004,
            "registered_mixture_id": reg_mix,
            "beta_cell_fraction": beta_f,
            "alpha_cell_fraction": alpha_f,
            "delta_cell_fraction": delta_f,
            "implanted_total_ieq": imp_ieq,
            "insulin_inventory_assay_specificity": "MATURE_INSULIN_SPECIFIC_LOW_PROINSULIN_CROSSREACTIVITY",
            "assay_cumulative_secreted_iu_equiv": 2.4,
            "lot_observations": lot_obs,
        }
        rec.update(overrides)
        return rec

    r_ok = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(_make_record(), root=ROOT)
    assert r_ok["subgate_pass"] is True
    assert "clinical_release_eligible" not in r_ok
    assert r_ok["release_authority"] == "DELEGATED_TO_LIVE_ATOMIC_STAGE0_LEDGER"
    assert r_ok["status"] == "QUALIFIED_TWOCLOCK_REPEATED_CHALLENGE_SUBGATE_EVIDENCE_ONLY"
    assert r_ok["proinsulin_to_cpeptide_marker_classification"] == "EXPLORATORY_STAGE0_IN_VITRO_PROCESSING_MARKER_PENDING_CALIBRATION"
    assert r_ok["basal_peak_statistics_and_windows"]["gate0m_strict_basal1_over_peak1_max"] == 0.15
    assert r_ok["basal_peak_statistics_and_windows"]["readout2_nominal_peak1_over_basal1_min"] == 4.5

    # Realistic analytical measurement noise around 60:30:10 stays on REF_60_30_10 (not misclassified as new DEC-084 recipe):
    r_noise = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(_make_record(beta_f=0.61, alpha_f=0.29, delta_f=0.10), root=ROOT)
    assert r_noise["subgate_pass"] is True

    # alpha:delta shift at constant beta=0.60 IS caught!
    r_ad_shift = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(_make_record(beta_f=0.60, alpha_f=0.15, delta_f=0.25), root=ROOT)
    assert r_ad_shift["subgate_pass"] is False and r_ad_shift["status"] == "HOLD_MEASURED_COMPOSITION_DRIFT_FROM_REGISTERED_MIXTURE_REQUIRES_DEC084_QUALIFICATION"

    # Registry execution_supported holds preserved on B2 and Design Z:
    r_b2 = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(_make_record(pvid="ISLET1_B2_5LOCUS_BASE"), root=ROOT)
    r_z = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(_make_record(pvid="ISLET1_Z_X_MOL_140UM_ZERO_CAO2"), root=ROOT)
    assert r_b2["registry_execution_supported"] is False and r_b2["status"] == "RESEARCH_COMPARATOR_QUALIFIED_NON_EXECUTABLE"
    assert r_z["registry_execution_supported"] is False and r_z["status"] == "RESEARCH_COMPARATOR_QUALIFIED_NON_EXECUTABLE"

    # Extended Rule S8 passes when paired by block + both 40/100 mmHg + early timestamp:
    early_paired = {
        f"CONF_{i+1:02d}": {
            "differentiation_block_id": f"DIFF_BLK_{i+1:02d}",
            "trace_id_40mmhg": f"TR_EARLY_40_CONF_{i+1:02d}",
            "z2_40mmhg": v,
            "trace_id_100mmhg": f"TR_EARLY_100_CONF_{i+1:02d}",
            "z2_100mmhg": v * 1.04,
            "early_harvest_utc": "2026-10-22T00:00:00Z",
        }
        for i, v in enumerate([0.82, 0.86, 0.78, 0.85, 0.80, 0.83])
    }
    r_s8 = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(
        _make_record(
            assay_harvest_utc="2026-11-05T00:00:00Z",
            stage7_total_days=35.0,
            days_post_sort_reaggregation=21.0,
            paired_earlier_window_z2_by_lot=early_paired,
        ),
        root=ROOT,
    )
    assert r_s8["subgate_pass"] is True and r_s8["is_extended_s8_exposure"] is True

    # Adversarial negative controls covering all 4 independent review hardenings:
    nc_records = [
        ("TINY_0P001_INCREMENTS_BLOCKED", _make_record(_per_lot_overrides={"low2_mean": [0.400, 0.401, 0.399, 0.400, 0.402, 0.398], "high2_mean": [0.401, 0.402, 0.400, 0.401, 0.403, 0.399]})),
        ("BELOW_HALF_ADULT_REF_Z2_BLOCKED", _make_record(_per_lot_overrides={"adult_islet_ref_z2": [2.8, 3.0, 2.7, 2.9, 2.85, 2.95]})),
        ("TINY_ADULT_REF_BLOCKED", _make_record(_per_lot_overrides={"adult_islet_ref_z2": [0.001, 0.002, 0.001, 0.0015, 0.0012, 0.0018]})),
        ("GATE0M_GLUCAGON_FAILED", _make_record(_per_lot_overrides={"glucagon_2mM_fold_increase": [1.4, 1.5, 1.3, 1.45, 1.55, 1.35]})),
        ("GATE0M_UCN3_FAILED", _make_record(_per_lot_overrides={"ucn3_maturity_ratio": [0.55, 0.58, 0.52, 0.57, 0.56, 0.54]})),
        ("TERMINAL_ALK5I_CONTINUED", _make_record(terminal_alk5i_tgfb_inhibitor_continued=True)),
        ("MISSING_PRE_IEC_VASCULAR_FUNCTION", _make_record(pre_coassembly_iec_vascular_function_verified=False)),
        ("BAD_ZM_SPIKE_RECOVERY", _make_record(zm447439_matrix_spike_recovery_pct=62.0)),
        ("TIMESTAMP_DAY_MISMATCH", _make_record(assay_harvest_utc="2026-10-30T00:00:00Z")),
        ("DUPLICATE_DIFF_BLOCK", _make_record(_per_lot_overrides={"differentiation_block_id": ["DIFF_BLK_01"] * 6})),
        ("MISSING_100MMHG_SIBLING", _make_record(_per_lot_overrides={"matched_100mmhg_assay_po2_mmhg": [40.0] * 6})),
        ("SINGLE_SLOW_CHAMBER_0P55", _make_record(_per_lot_overrides={"inert_tracer_tau_chamber_min": [0.22, 0.24, 0.23, 0.25, 0.21, 0.55]})),
        ("UNALIGNED_TRACER", _make_record(_per_lot_overrides={"tracer_time_aligned": [True, True, True, False, True, True]})),
        ("SHORT_2MIN_CENSORED_LOW2", _make_record(_per_lot_overrides={"low2_washout_duration_min": [2.0] * 6})),
        ("PER_CELL_ETA_ON_PER_DNA", _make_record(hormone_unit="uIU_insulin_per_ug_dna_per_min", normalization_denominator="ug_presort_viable_endocrine_dna", recovery_denominator_basis="1e6_presort_viable_endocrine_cells")),
        ("CROSS_REACTIVE_INSULIN_ELISA", _make_record(insulin_inventory_assay_specificity="GENERIC_CROSS_REACTIVE_ELISA")),
        ("CALLER_GATE0S_BOOLEAN_OVERRIDE", _make_record(gate0s_elevated_stored_insulin_clamp_reverified=True)),
        ("DEC084_70_20_10_WITH_N6", _make_record(n_lots=6, reg_mix="DEC084_M1_70_20_10", beta_f=0.70, alpha_f=0.20, delta_f=0.10, dec084_matched_age_mixture_gate_passed=True)),
        ("STORED_MAX_LOT_POOLING_BYPASS_150IU", _make_record(_per_lot_overrides={"stored_mature_insulin_iu_implanted": [84.0, 88.0, 90.0, 92.0, 94.0, 150.0]})),
    ]
    for nc_name, nc_rec in nc_records:
        nc_out = wc_mod.evaluate_stage0_twoclock_and_repeated_challenge_iip_0029(nc_rec, root=ROOT)
        assert nc_out["subgate_pass"] is False, f"Expected {nc_name} to fail closed, got {nc_out}"

    # (c) Identical >138.3 IU Gate 0S escalation, wrong-lot rejection & source-replacement latching across ALL 4 paths:
    rec_hi_store = _make_record(_per_lot_overrides={"stored_mature_insulin_iu_implanted": [144.0, 147.0, 145.0, 149.0, 146.0, 148.0]})
    base_core_s9 = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "current_lot_id": "LOT_S7_2026_01",
        "assay_campaign_id": "CAMP_S0_29",
        "construct_sha256": valid_sha,
        "device_assembly_id": "DEV_29_01",
        "obs_icasp9_log10_kill": 4.8,
        "obs_icasp9_log10_kill_se": 0.12,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300000000,
        "obs_chop_coincidence_fold": 1.1,
        "obs_chop_coincidence_se": 0.05,
        "obs_hc10_surface_fold_wt": 1.05,
        "obs_hc10_surface_fold_se": 0.05,
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "obs_post_edit_segmented_baf_chr19q": 0.50,
        "obs_post_edit_segmented_baf_chr11p15": 0.50,
    }
    # 1. Missing Gate 0S when >138.3 IU -> identical HOLD across Flat, Direct Guardrails, Ingestor, and Replay:
    flat_hi = wc_mod.evaluate_stage0_decision_tree({**base_core_s9, "iip0029_twoclock_record": rec_hi_store}, root=ROOT)
    assert flat_hi["resolved_architecture_assignments"]["lot_maturation_status"] == "HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW"
    dir_hi = wc_mod.evaluate_stage0_research_and_molecular_guardrails({**base_core_s9, "lot_id": "LOT_S7_2026_01", "campaign_id": "CAMP_S0_29", "iip0029_twoclock_record": rec_hi_store}, root=ROOT)
    assert dir_hi["overall_status"] == "HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW"
    ing_doc_hi = {
        "scope": {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": "LOT_S7_2026_01", "campaign_id": "CAMP_S0_29", "construct_sha256": valid_sha, "device_assembly_id": "DEV_29_01"},
        "measurements": {k: v for k, v in base_core_s9.items() if k.startswith("obs_")},
        "iip0026_molecular_observations": {"donor_seed_gnomad_target_site_wgs_verified": True, "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True},
        "iip0029_twoclock_record": rec_hi_store,
    }
    ing_hi = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(ing_doc_hi), root=ROOT)
    assert ing_hi["overall_status"] == "HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW"

    # 2. Wrong lot inside iip0029_twoclock_record -> FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH identically:
    rec_wrong_lot = _make_record(lot_id="OTHER_LOT_99")
    assert wc_mod.evaluate_stage0_decision_tree({**base_core_s9, "iip0029_twoclock_record": rec_wrong_lot}, root=ROOT)["resolved_architecture_assignments"]["lot_maturation_status"] == "FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH"
    assert ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps({**ing_doc_hi, "iip0029_twoclock_record": rec_wrong_lot}), root=ROOT)["overall_status"] == "FAIL_IIP0029_WRONG_LOT_OR_PRODUCT_SCOPE_MISMATCH"

    # 3. Scope-matched passing Gate 0S clears the >138.3 IU escalation; breached >2.8 IU burst fails closed:
    g0s_pass = {
        "obs_gate0s_ap1903_1h_insulin_burst_iu": 1.6,
        "obs_gate0s_ap1903_1h_insulin_burst_iu_se": 0.2,
        "obs_gate0s_ap1903_supernatant_insulin_dump_fraction_24h": 0.08,
        "obs_gate0s_ap1903_supernatant_insulin_dump_se": 0.01,
        "obs_explant_en_bloc_recovery_fraction": 0.996,
        "obs_explant_en_bloc_recovery_se": 0.001,
        "gate0s_elevated_store_scope_matched_device_lot_id": "LOT_S7_2026_01__DEV_29_01",
    }
    _t6b_pass_54 = [
        {
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "kind": "EXACT_INTEGER_COUNT",
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "product_variant_id": "ISLET1_A4_ARM1B",
            "construct_sha256": valid_sha,
            "donor_seedline_id": "DONOR1_H1_GMP",
            "mcb_id": "MCB_D1_001",
            "lot_id": "LOT_S7_2026_01",
            "campaign_id": "CAMP_S0_29",
            "fraction_id": "enriched_cd31_neg_cd144_neg",
            "declared_fraction_roster": ["enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"],
            "sample_id": "T54_S1",
            "k": 0,
            "N": 30_000_000,
            "w_i": 0.9999,
            "hpsc_calibrant": (20, 20),
            "edited_differentiated_calibrant": (20, 20),
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "untreated_cfu_control": {"untreated_sample_id": "T54_UNTR_1", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": 3_000_000_000,
            "presort_reconciled_loss_cells": 0,
            "presort_fraction_cells": 2_999_700_000,
            "plated_treated_cells_n": 80_000_000,
        },
        {
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "kind": "EXACT_INTEGER_COUNT",
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "product_variant_id": "ISLET1_A4_ARM1B",
            "construct_sha256": valid_sha,
            "donor_seedline_id": "DONOR1_H1_GMP",
            "mcb_id": "MCB_D1_001",
            "lot_id": "LOT_S7_2026_01",
            "campaign_id": "CAMP_S0_29",
            "fraction_id": "flow_through_and_excluded",
            "declared_fraction_roster": ["enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"],
            "sample_id": "T54_S2",
            "k": 0,
            "N": 100_000,
            "w_i": 0.0001,
            "hpsc_calibrant": (10, 10),
            "edited_differentiated_calibrant": (10, 10),
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "untreated_cfu_control": {"untreated_sample_id": "T54_UNTR_2", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": 3_000_000_000,
            "presort_reconciled_loss_cells": 0,
            "presort_fraction_cells": 300_000,
            "plated_treated_cells_n": 260_000,
        },
    ]
    ing_doc_cleared = {
        **ing_doc_hi,
        "scope": {**ing_doc_hi["scope"], "donor_seedline_id": "DONOR1_H1_GMP", "mcb_id": "MCB_D1_001"},
        "measurements": {**ing_doc_hi["measurements"], **{k: v for k, v in g0s_pass.items() if k.startswith("obs_")}},
        "iip0026_molecular_observations": {**ing_doc_hi["iip0026_molecular_observations"], "gate0s_elevated_store_scope_matched_device_lot_id": "LOT_S7_2026_01__DEV_29_01"},
        "tier6b_sorted_lot_observations": _t6b_pass_54,
    }
    assert ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(ing_doc_cleared), root=ROOT)["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED"
    ing_doc_breach = {
        **ing_doc_cleared,
        "measurements": {**ing_doc_cleared["measurements"], "obs_gate0s_ap1903_1h_insulin_burst_iu": 3.2, "obs_gate0s_ap1903_1h_insulin_burst_iu_se": 0.2},
    }
    assert ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(ing_doc_breach), root=ROOT)["overall_status"] == "FAIL_GATE0S_ELEVATED_STORE_BURST_OR_DUMP_BREACHED"

    # 4. Chronological Replay latches immutable source-replacement & threshold-lowering mutation:
    rep_mut = ing_mod.replay_stage0_chronological_batches(
        [
            {"utc_timestamp": "2026-10-22T01:00:00Z", "seq_index": 1, **ing_doc_hi},
            {"utc_timestamp": "2026-10-22T02:00:00Z", "seq_index": 2, **{**ing_doc_hi, "iip0029_twoclock_record": _make_record()}},
        ],
        root=ROOT,
    )
    assert rep_mut["overall_status"] == "FAIL_IIP0029_IMMUTABLE_SOURCE_OR_THRESHOLD_REPLACEMENT_DETECTED"

    # 5. Live Atomic Stage-0 Ledger enforces identical >138.3 IU escalation, scope-matched Gate 0S resolution, and source-replacement hold:
    mcb_k = wc_mod.MCBScopeKey("ISLET1_A4_ARM1B", valid_sha, "DONOR1_H1_GMP", "MCB_D1_001")
    lot_k = wc_mod.LotScopeKey("ISLET1_A4_ARM1B", valid_sha, "DONOR1_H1_GMP", "MCB_D1_001", "LOT_S7_2026_01")
    dev_k = wc_mod.AssembledDeviceLotScopeKey("ISLET1_A4_ARM1B", valid_sha, "DONOR1_H1_GMP", "MCB_D1_001", "LOT_S7_2026_01", "DEV_29_01", "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH")
    mcb_atoms = [
        wc_mod.AtomicObservation(endpoint_id=eid, scope_type="MCB_SCOPE", scope_key=mcb_k, kind="METADATA_CONTROL", value_or_k=True, campaign_id="CAMP_S0_29", sample_unit_id=f"MCB_U_{idx}")
        for idx, eid in enumerate(wc_mod.ENDPOINT_REGISTRY_SPEC["mcb_mandatory_endpoint_ids"])
    ]
    t6b_atoms_54 = [
        wc_mod.AtomicObservation(
            endpoint_id="lot_tier6b_sorted_prolif_clonogenic_escape",
            scope_type="LOT_SCOPE",
            scope_key=lot_k,
            kind="EXACT_INTEGER_COUNT",
            value_or_k=0,
            se_or_n=r["N"],
            timepoint=r["timepoint"],
            assay_version=r["assay_version"],
            unit=r["unit"],
            campaign_id="CAMP_S0_29",
            sample_unit_id=r["sample_id"],
            controls={k: v for k, v in r.items() if k not in ("endpoint_id", "scope_type", "kind", "k", "N", "timepoint", "assay_version", "unit")},
        )
        for r in _t6b_pass_54
    ]
    sister_term_atoms = [
        wc_mod.AtomicObservation(endpoint_id=eid, scope_type="ASSEMBLED_DEVICE_LOT_SCOPE", scope_key=dev_k, kind="METADATA_CONTROL", value_or_k=True, campaign_id="CAMP_S0_29", sample_unit_id=f"DEV_U_{idx}")
        for idx, eid in enumerate(wc_mod.ENDPOINT_REGISTRY_SPEC["device_sister_terminal_endpoint_ids"])
        if eid != "obs_explant_en_bloc_recovery_fraction"
    ] + [
        wc_mod.AtomicObservation(endpoint_id="obs_explant_en_bloc_recovery_fraction", scope_type="ASSEMBLED_DEVICE_LOT_SCOPE", scope_key=dev_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.996, se_or_n=0.001, sample_n=6, campaign_id="CAMP_S0_29", sample_unit_id="DEV_RETR_01")
    ]
    panel_hi_atom = wc_mod.AtomicObservation(endpoint_id="stage0_twoclock_repeated_challenge_panel", scope_type="LOT_SCOPE", scope_key=lot_k, kind="STRUCTURED_PANEL", value_or_k=rec_hi_store, campaign_id="CAMP_S0_29", sample_unit_id="LOT_PANEL_01")
    s9_lot_atoms_54 = [
        wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=4.8, se_or_n=0.12, sample_n=6, campaign_id="CAMP_S0_29", sample_unit_id="LOT29_IC9"),
        wc_mod.AtomicObservation(endpoint_id="obs_tp53_ppm1d_bcor_vaf", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6, campaign_id="CAMP_S0_29", sample_unit_id="LOT29_TP53", controls={"obs_tp53_ppm1d_bcor_duplex_depth": 70000}),
        wc_mod.AtomicObservation(endpoint_id="obs_castseq_translocation_freq", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6, campaign_id="CAMP_S0_29", sample_unit_id="LOT29_CAST", controls={"obs_castseq_duplex_or_amplicon_depth": 70000}),
    ]
    target_prov = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "construct_sha256": valid_sha,
        "donor_seedline_id": "DONOR1_H1_GMP",
        "mcb_id": "MCB_D1_001",
        "lot_id": "LOT_S7_2026_01",
        "device_assembly_id": "DEV_29_01",
        "physical_route_id": "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH",
    }
    led_missing_g0s = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"batch_id": "B1", "atoms": mcb_atoms + s9_lot_atoms_54 + t6b_atoms_54 + [panel_hi_atom] + sister_term_atoms}], target_prov, target_campaign_id="CAMP_S0_29", root=ROOT)
    assert led_missing_g0s["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE"
    assert "HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW" in led_missing_g0s["holds"]

    g0s_atoms_pass = [
        wc_mod.AtomicObservation(endpoint_id="gate0s_ap1903_1h_insulin_burst_iu", scope_type="ASSEMBLED_DEVICE_LOT_SCOPE", scope_key=dev_k, kind="CONTINUOUS_STUDENT_T", value_or_k=1.6, se_or_n=0.2, sample_n=6, campaign_id="CAMP_S0_29", sample_unit_id="DEV_BURST_01"),
        wc_mod.AtomicObservation(endpoint_id="gate0s_ap1903_24h_insulin_dump_fraction", scope_type="ASSEMBLED_DEVICE_LOT_SCOPE", scope_key=dev_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.08, se_or_n=0.01, sample_n=6, campaign_id="CAMP_S0_29", sample_unit_id="DEV_DUMP_01"),
    ]
    led_cleared = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"batch_id": "B1", "atoms": mcb_atoms + s9_lot_atoms_54 + t6b_atoms_54 + [panel_hi_atom] + sister_term_atoms + g0s_atoms_pass}], target_prov, target_campaign_id="CAMP_S0_29", root=ROOT)
    assert led_cleared["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE"
    assert "HOLD_STORED_INSULIN_EXCEEDS_138P3_IU_REQUIRES_FULL_LINKED_GATE0S_LEDGER_REVIEW" not in led_cleared["holds"]
    assert "HOLD_INVALID_ATOM_KIND_OR_NON_NUMERIC_VALUE:mcb_1q32_mdm4_single_cell_mosaic" in led_cleared["holds"]
    assert "HOLD_INVALID_ATOM_KIND_OR_NON_NUMERIC_VALUE:gate0s_shedding_ddpcr_assay" in led_cleared["holds"]

    # (d) Exploratory Mode B v2 sensitivity comparator & Mode A default parity:
    mode_a_def = l4_mod.simulate_metabolic_twin_granule_sensitivity_iip_0029(granule_conservation_mode="mode_a_default")
    assert mode_a_def["tir_3p9_to_10_elapsed_pct"] == 100.0 and mode_a_def["min_glucose_mM"] == 4.59
    for bad_fn in (
        lambda: l4_mod.ModeBV2GranuleStoreParams(f_rrp=1.01),
        lambda: l4_mod.ModeBV2GranuleStoreParams(tau_chamber_min=0.002),
        lambda: l4_mod.build_mode_b_v2_granule_store_system_iip_0029(l4_mod.ModeBV2GranuleStoreParams(), extravascular_insulin_lag_total_min=0.01),
        lambda: l4_mod.simulate_metabolic_twin_granule_sensitivity_iip_0029(granule_conservation_mode="bogus_mode"),
    ):
        try:
            bad_fn()
            raise AssertionError("Expected ValueError on invalid Mode B v2 parameter/mode")
        except ValueError:
            pass

    mb_audit = l4_atlas["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["exploratory_two_clock_granule_store_sensitivity_audit_iip_0029"]
    assert mb_audit["evidence_status"] == "EXPLORATORY_OPT_IN_SENSITIVITY_COMPARATOR_NOT_A_RELEASE_GATE"
    assert mb_audit["default_mode_a_unchanged"] is True
    assert mb_audit["real_scipy_dop853_multi_state_verification"]["verified"] is True
    assert mb_audit["real_scipy_dop853_multi_state_verification"]["dop853_0_to_1440min_max_abs_g_diff_mM"] < 1e-5
    assert mb_audit["real_scipy_dop853_multi_state_verification"]["dop853_0_to_1440min_max_abs_qp_diff_iu"] < 1e-6
    assert mb_audit["real_scipy_dop853_multi_state_verification"]["dop853_0_to_1440min_max_abs_qr_diff_iu"] < 1e-5
    assert mb_audit["real_scipy_dop853_multi_state_verification"]["dop853_0_to_1440min_max_abs_cum_prehep_diff_iu"] < 1e-5
    assert mb_audit["effluent_shutoff_and_short_interval_censoring_audit"]["tau_chamber_0p60_min_2min_low_t90_censored"] is None
    assert mb_audit["effluent_shutoff_and_short_interval_censoring_audit"]["tau_chamber_0p60_min_60min_low_t90_min"] > 3.0
    assert mb_audit["positivity_mass_conservation_and_biosyn_audit"]["c0_0p25_iu_min_qp_iu"] > 0.0
    assert mb_audit["positivity_mass_conservation_and_biosyn_audit"]["c0_0p25_iu_relative_mass_error"] < 5e-12
    assert mb_audit["ceteris_paribus_isolated_tau_refill_meal_sensitivity"]["canonical_5h_7h_250g_tir_tau16p5_pct"] >= 85.0
    assert mb_audit["ceteris_paribus_isolated_tau_refill_meal_sensitivity"]["canonical_5h_7h_250g_tir_tau45p0_pct"] < 80.0

    # (e) Graveyard GY-024 (10 grid + 15 off-grid = 25 total after IIP-0031 GY-025):
    gy_doc = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    og_list = gy_doc["off_grid_molecular_dead_end_guardrails"]["guardrails"]
    assert len(gy_doc["graveyard_entries"]) == 10 and len(og_list) == 15
    gy24 = [g for g in og_list if g["graveyard_card_id"] == "GY-024"]
    assert len(gy24) == 1 and gy24[0]["guardrail_id"].startswith("OFFGRID_GUARD_014_") and gy24[0]["contribution_id"] == "IIP-0029"


@pytest.mark.parametrize("malformed", [None, [], False, 0, "", "null"])
def test_s11_nonobject_panel_payload_holds(malformed):
    from scripts.wetlab_contract import evaluate_locus2_ciita_modeb_challenger_v3_gate

    result = evaluate_locus2_ciita_modeb_challenger_v3_gate(malformed)
    assert result["gate_pass"] is False
    assert result["status"] == "HOLD_S11_MALFORMED_RECORD"


@pytest.mark.parametrize("invalid", [None, [], {}, 11, True, "", " lot "])
@pytest.mark.parametrize("key", ["lot_id", "current_lot_id", "sample_id", "sample_unit_id"])
def test_s11_invalid_identity_alias_holds(key, invalid):
    from scripts.wetlab_contract import evaluate_locus2_ciita_modeb_challenger_v3_gate

    record = {"scope_type": "LOT_SCOPE", key: invalid}
    result = evaluate_locus2_ciita_modeb_challenger_v3_gate({"iip0031_locus2_record": record})
    assert result["gate_pass"] is False
    assert result["status"] == "HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING"


@pytest.mark.parametrize("outer,nested", [(3, 4), (4, 3), (True, 1), (1, True), (None, 1), (1, None)])
def test_s11_conflicting_nested_evidence_holds(outer, nested):
    from scripts.wetlab_contract import evaluate_locus2_ciita_modeb_challenger_v3_gate

    key = "obs_stage0_independent_lot_replicates_n"
    payload = {key: outer, "iip0031_locus2_record": {"scope_type": "LOT_SCOPE", key: nested}}
    result = evaluate_locus2_ciita_modeb_challenger_v3_gate(payload)
    assert result["status"] == "HOLD_S11_CONFLICTING_NESTED_EVIDENCE"
    assert result["gate_pass"] is False


@pytest.mark.parametrize("reported,alias", [
    (6, 6), (10_000, 10_000), (10_001, 10_001), (10**308, 10**308),
    (10**350, 10**350), (6, 100), (100, 6), (6, True), (6, 6.0),
])
def test_s11_replicates_respect_count_type_aliases_and_contract_range(reported, alias):
    from scripts.wetlab_contract import evaluate_locus2_ciita_modeb_challenger_v3_gate

    rec = {"scope_type": "LOT_SCOPE", "obs_stage0_independent_lot_replicates_n": reported, "replicates_n": alias}
    result = evaluate_locus2_ciita_modeb_challenger_v3_gate({"iip0031_locus2_record": rec})
    assert result["gate_pass"] is False
    valid = type(alias) is int and reported == alias and 6 <= reported <= 10_000
    assert ("HOLD_S11_INVALID_OR_MISSING_REPLICATES_N" in result["reasons"]) == (not valid)


@pytest.mark.parametrize("payload", [
    {"iip0031_locus2_record": {"lot_id": "lot"}},
    {"obs_l2_v3_adcc_intact_specimen_id": "sample"},
])
def test_s11_requires_explicit_scope(payload):
    from scripts.wetlab_contract import evaluate_locus2_ciita_modeb_challenger_v3_gate

    result = evaluate_locus2_ciita_modeb_challenger_v3_gate(payload)
    assert result["status"] == "HOLD_S11_PANEL_REQUIRES_LOT_SCOPE"
    assert result["gate_pass"] is False


def test_55_iip0031_locus2_modeb_comparators_rule_s11_and_tier6b_per_fraction_thinning():
    """
    IIP-0031 (DEC-110; GY-025 / OFFGRID_GUARD_015_*; Rule S11):
    Verifies all 7 exact Tier 6B per-fraction thinning acceptance fixtures across ALL FOUR entry points
    (direct product contract, atomic Stage-0 ledger, serialized JSON ingestion, chronological replay),
    all 8 non-executable L2_V3_* products, Rule S11 hard/numerical/attribution gates, Layer-2 v3 manifest
    & <=1 mm >=20 bp disjoint IR scanner, Layer-3 2WII/5O32/Val62 provenance, GY-025, and bit-for-bit
    Layer-1 2,160/2,160 & OQ-016 invariance.
    """
    import hashlib
    import importlib.util
    import json
    from pathlib import Path
    from typing import Any

    ROOT = Path(__file__).resolve().parents[1]

    def _load_mod(name: str, path: Path):
        spec = importlib.util.spec_from_file_location(name, path)
        assert spec and spec.loader
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    wc_mod = _load_mod("wc_iip0031_55", ROOT / "scripts/wetlab_contract.py")
    ing_mod = _load_mod("ing_iip0031_55", ROOT / "scripts/ingest_wetlab_results.py")
    pc_mod = _load_mod("pc_iip0031_55", ROOT / "scripts/product_contract.py")
    l2_mod = _load_mod("l2_iip0031_55", ROOT / "scripts/run_layer2_perturbation_and_regulatory_plane.py")

    p_defs = pc_mod.load_product_contracts(ROOT)
    by_id = {p["product_variant_id"]: p for p in p_defs["products"]}
    assert len(p_defs["products"]) == 27
    l2_v3_ids = sorted(pc_mod.STAGED_LOCUS2_V3_PRODUCT_IDS)
    assert len(l2_v3_ids) == 8

    # Verify locked v3 FASTA SHA-256 and every per-record digest:
    v3_fasta_path = ROOT / "data/constructs/islet1_locus2_ciita_modeb_comparators_v0.8.0.fasta"
    assert v3_fasta_path.exists()
    assert hashlib.sha256(v3_fasta_path.read_bytes()).hexdigest() == "25ffcb10e2bad77a4b34725c7eff0532990ebf7b80864cf7217d33f5ccbe5b6d"
    fasta_seqs = {}
    cur_h = None
    cur_buf = []
    for raw_ln in v3_fasta_path.read_text(encoding="utf-8").splitlines():
        s_ln = raw_ln.strip()
        if not s_ln:
            continue
        if s_ln.startswith(">"):
            if cur_h is not None:
                fasta_seqs[cur_h] = "".join(cur_buf).upper()
            cur_h = s_ln[1:].split()[0]
            cur_buf = []
        else:
            cur_buf.append(s_ln)
    if cur_h is not None:
        fasta_seqs[cur_h] = "".join(cur_buf).upper()
    assert len(fasta_seqs) == 8

    for pvid in l2_v3_ids:
        pdef = by_id[pvid]
        assert pdef["execution_supported"] is False
        assert pdef["lifecycle_state"] == "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"
        l2b = pdef["staged_locus2_modeb_donor_binding"]
        did = l2b["locus2_donor_id"]
        assert did in fasta_seqs
        assert l2b["locus2_donor_sha256"] == hashlib.sha256(fasta_seqs[did].encode("ascii")).hexdigest()
        assert l2b["rule_s11_gate_id"] == "S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate"

    # Helper to construct clean baseline executable A4 lot across ALL FOUR entry points:
    a4_sha = "f74728767987c18bdcb8c9d49773ee8dbe10e6006e57aa622638e4fb07062b0b"
    LOT_ID = "LOT_S7_IIP0031_01"
    CAMP_ID = "CAMP_S0_31"
    SEED_ID = "DONOR1_H1_GMP"
    MCB_ID = "MCB_D1_001"
    DEV_ID = "DEV_31_01"
    ROUTE_ID = "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH"

    mcb_k = wc_mod.MCBScopeKey("ISLET1_A4_ARM1B", a4_sha, SEED_ID, MCB_ID)
    lot_k = wc_mod.LotScopeKey("ISLET1_A4_ARM1B", a4_sha, SEED_ID, MCB_ID, LOT_ID)
    dev_k = wc_mod.AssembledDeviceLotScopeKey("ISLET1_A4_ARM1B", a4_sha, SEED_ID, MCB_ID, LOT_ID, DEV_ID, ROUTE_ID)
    target_prov = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "construct_sha256": a4_sha,
        "donor_seedline_id": SEED_ID,
        "mcb_id": MCB_ID,
        "lot_id": LOT_ID,
        "device_assembly_id": DEV_ID,
        "physical_route_id": ROUTE_ID,
    }
    mcb_atoms_base = [
        wc_mod.AtomicObservation(endpoint_id=eid, scope_type="MCB_SCOPE", scope_key=mcb_k, kind="METADATA_CONTROL", value_or_k=True, campaign_id=CAMP_ID, sample_unit_id=f"MCB31_{i}")
        for i, eid in enumerate(wc_mod.ENDPOINT_REGISTRY_SPEC["mcb_mandatory_endpoint_ids"])
    ]
    dev_atoms_base = [
        wc_mod.AtomicObservation(endpoint_id=eid, scope_type="ASSEMBLED_DEVICE_LOT_SCOPE", scope_key=dev_k, kind="METADATA_CONTROL", value_or_k=True, campaign_id=CAMP_ID, sample_unit_id=f"DEV31_{i}")
        for i, eid in enumerate(wc_mod.ENDPOINT_REGISTRY_SPEC["device_sister_terminal_endpoint_ids"])
        if eid != "obs_explant_en_bloc_recovery_fraction"
    ] + [
        wc_mod.AtomicObservation(endpoint_id="obs_explant_en_bloc_recovery_fraction", scope_type="ASSEMBLED_DEVICE_LOT_SCOPE", scope_key=dev_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.996, se_or_n=0.001, sample_n=6, campaign_id=CAMP_ID, sample_unit_id="DEV31_RETR")
    ]
    base_meas_31 = {
        "obs_icasp9_log10_kill": 4.8,
        "obs_icasp9_log10_kill_se": 0.12,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300000000,
        "obs_chop_coincidence_fold": 1.1,
        "obs_chop_coincidence_se": 0.05,
        "obs_hc10_surface_fold_wt": 1.05,
        "obs_hc10_surface_fold_se": 0.05,
        "obs_post_edit_segmented_baf_chr19q": 0.50,
        "obs_post_edit_segmented_baf_chr11p15": 0.50,
    }
    base_mol_31 = {
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
    }
    base_scope_31 = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "lot_id": LOT_ID,
        "campaign_id": CAMP_ID,
        "construct_sha256": a4_sha,
        "donor_seedline_id": SEED_ID,
        "mcb_id": MCB_ID,
        "device_assembly_id": DEV_ID,
    }

    def _mk_frac(
        fid: str,
        sid: str,
        *,
        k: Any = 0,
        N: Any = 30_000_000,
        w_i: Any = 0.9999,
        hpsc: Any = (20, 20),
        diff: Any = (20, 20),
        roster: Any = ("enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"),
        lot_id: str = LOT_ID,
        construct_sha256: str = a4_sha,
        donor_seedline_id: str = SEED_ID,
        mcb_id: str = MCB_ID,
        product_variant_id: str = "ISLET1_A4_ARM1B",
        endpoint_id: str = "lot_tier6b_sorted_prolif_clonogenic_escape",
        scope_type: str = "LOT_SCOPE",
        kind: str = "EXACT_INTEGER_COUNT",
        extra: Any = None,
    ) -> dict:
        tot_cells_default = 3_000_000_000
        if isinstance(w_i, (int, float)) and not isinstance(w_i, bool) and 0.0 < float(w_i) <= 1.0:
            frac_cells_default = max(1, int(round(float(w_i) * float(tot_cells_default))))
        else:
            frac_cells_default = 2_999_700_000 if "enriched" in str(fid) else 300_000
        plated_default = int((int(N) * 8 + 2) // 3) if (type(N) is int and N >= 1) else 80_000_000
        rec = {
            "endpoint_id": endpoint_id,
            "scope_type": scope_type,
            "kind": kind,
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "product_variant_id": product_variant_id,
            "construct_sha256": construct_sha256,
            "donor_seedline_id": donor_seedline_id,
            "mcb_id": mcb_id,
            "lot_id": lot_id,
            "campaign_id": CAMP_ID,
            "fraction_id": fid,
            "declared_fraction_roster": list(roster) if isinstance(roster, (list, tuple)) else roster,
            "sample_id": sid,
            "k": k,
            "N": N,
            "w_i": w_i,
            "hpsc_calibrant": list(hpsc) if isinstance(hpsc, tuple) else hpsc,
            "edited_differentiated_calibrant": list(diff) if isinstance(diff, tuple) else diff,
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "untreated_cfu_control": {
                "untreated_sample_id": f"{sid}_UNTR_CFU",
                "untreated_cfu_colonies_k": 4000,
                "untreated_plated_cells_n": 10000,
            },
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": tot_cells_default,
            "presort_reconciled_loss_cells": 0,
            "presort_fraction_cells": frac_cells_default,
            "plated_treated_cells_n": plated_default,
        }
        if extra:
            rec.update(extra)
        return rec

    def _to_atom(r: Any) -> Any:
        if not isinstance(r, dict):
            return r
        if "controls" in r and r["controls"] is not None and not isinstance(r["controls"], dict):
            raise ValueError("malformed non-dict controls routed to ledger raw fallback lane")
        st = r.get("scope_type", "LOT_SCOPE")
        pvid_r = r.get("product_variant_id", "ISLET1_A4_ARM1B")
        sha_r = r.get("construct_sha256", a4_sha)
        seed_r = r.get("donor_seedline_id", SEED_ID)
        mcb_r = r.get("mcb_id", MCB_ID)
        lot_r = r.get("lot_id", LOT_ID)
        if st == "MCB_SCOPE":
            sk = wc_mod.MCBScopeKey(pvid_r, sha_r, seed_r, mcb_r)
        elif st == "ASSEMBLED_DEVICE_LOT_SCOPE":
            sk = wc_mod.AssembledDeviceLotScopeKey(pvid_r, sha_r, seed_r, mcb_r, lot_r, DEV_ID, ROUTE_ID)
        else:
            sk = wc_mod.LotScopeKey(pvid_r, sha_r, seed_r, mcb_r, lot_r)
        atom_meta_kw = {}
        for _mk in ("timepoint", "assay_version", "unit"):
            if _mk in r and isinstance(r[_mk], str) and r[_mk].strip():
                atom_meta_kw[_mk] = r[_mk]
        ctrl = {
            k: v for k, v in r.items()
            if k not in (
                "endpoint_id", "scope_type", "kind", "k", "N",
                "product_variant_id", "construct_sha256", "donor_seedline_id",
                "mcb_id", "lot_id", "controls", "timepoint", "assay_version",
                "unit", "campaign_id", "sample_n",
            )
        }
        if isinstance(r.get("controls"), dict):
            ctrl.update(r["controls"])
        return wc_mod.AtomicObservation(
            endpoint_id=r.get("endpoint_id", "lot_tier6b_sorted_prolif_clonogenic_escape"),
            scope_type=st,
            scope_key=sk,
            kind=r.get("kind", "EXACT_INTEGER_COUNT"),
            value_or_k=r.get("k"),
            se_or_n=r.get("N"),
            campaign_id=CAMP_ID,
            sample_unit_id=str(r.get("sample_id", "S_UNK")),
            controls=ctrl,
            **atom_meta_kw,
        )

    def _eval_4paths(
        t6b_records: list,
        *,
        scope_override: Any = None,
        target_prov_override: Any = None,
        meas_override: Any = None,
        s11_record: Any = None,
        include_s11_key: bool = False,
        s11_ledger_scope: str = "LOT_SCOPE",
        extra_ledger_atoms: Any = None,
        extra_ledger_batches: Any = None,
        extra_replay_records: Any = None,
        split_replay_batches: bool = True,
    ) -> dict:
        scope_use = {k: v for k, v in {**base_scope_31, **(scope_override or {})}.items() if v is not None}
        tprov_use = {k: v for k, v in {**target_prov, **(target_prov_override or {})}.items() if v is not None}
        meas_use = {**base_meas_31, **(meas_override or {})}
        s11_dict_part = {"iip0031_locus2_record": s11_record} if (include_s11_key or s11_record is not None) else {}

        # Path 1: Direct product contract
        p1 = wc_mod.evaluate_stage0_product_contract(
            {**scope_use, **meas_use, **base_mol_31, **s11_dict_part, "tier6b_sorted_lot_observations": list(t6b_records)},
            root=ROOT,
        )

        # Path 2: Atomic Stage-0 Ledger (strict AtomicObservation where constructible; raw tier6b_sorted_lot_observations
        # when testing invalid non-int/bool/zero-N/malformed count payloads rejected by AtomicObservation.__init__):
        pvid_l2 = tprov_use.get("product_variant_id", "ISLET1_A4_ARM1B")
        sha_l2 = tprov_use.get("construct_sha256", a4_sha)
        seed_l2 = tprov_use.get("donor_seedline_id", SEED_ID)
        mcb_l2 = tprov_use.get("mcb_id", MCB_ID)
        lot_l2 = tprov_use.get("lot_id", LOT_ID)
        mcb_k_use = wc_mod.MCBScopeKey(pvid_l2, sha_l2, seed_l2, mcb_l2)
        lot_k_use = wc_mod.LotScopeKey(pvid_l2, sha_l2, seed_l2, mcb_l2, lot_l2)
        dev_k_use = wc_mod.AssembledDeviceLotScopeKey(pvid_l2, sha_l2, seed_l2, mcb_l2, lot_l2, DEV_ID, ROUTE_ID)
        extra_eids = {getattr(a, "endpoint_id", None) for a in (extra_ledger_atoms or [])}
        _vaf_n_l2 = int(meas_use.get("obs_tp53_ppm1d_bcor_vaf_n", meas_use.get("obs_tp53_ppm1d_bcor_vaf_sample_n", meas_use.get("sample_n", 6))))
        _cast_n_l2 = int(meas_use.get("obs_castseq_translocation_freq_n", meas_use.get("obs_castseq_translocation_n", meas_use.get("sample_n", 6))))
        s9_lot_atoms_use = []
        if "obs_icasp9_log10_kill" not in extra_eids and meas_use.get("obs_icasp9_log10_kill") is not None:
            s9_lot_atoms_use.append(
                wc_mod.AtomicObservation(
                    endpoint_id="obs_icasp9_log10_kill",
                    scope_type="LOT_SCOPE",
                    scope_key=lot_k_use,
                    kind="CONTINUOUS_STUDENT_T",
                    value_or_k=meas_use["obs_icasp9_log10_kill"],
                    se_or_n=meas_use["obs_icasp9_log10_kill_se"],
                    sample_n=6,
                    campaign_id=CAMP_ID,
                    sample_unit_id="LOT31_IC9",
                )
            )
        if "obs_tp53_ppm1d_bcor_vaf" not in extra_eids and meas_use.get("obs_tp53_ppm1d_bcor_vaf") is not None:
            _vaf_ctrl_l2 = {}
            if meas_use.get("obs_tp53_ppm1d_bcor_duplex_depth") is not None:
                _vaf_ctrl_l2["obs_tp53_ppm1d_bcor_duplex_depth"] = meas_use["obs_tp53_ppm1d_bcor_duplex_depth"]
            s9_lot_atoms_use.append(
                wc_mod.AtomicObservation(
                    endpoint_id="obs_tp53_ppm1d_bcor_vaf",
                    scope_type="LOT_SCOPE",
                    scope_key=lot_k_use,
                    kind="CONTINUOUS_STUDENT_T",
                    value_or_k=meas_use["obs_tp53_ppm1d_bcor_vaf"],
                    se_or_n=meas_use["obs_tp53_ppm1d_bcor_vaf_se"],
                    sample_n=_vaf_n_l2,
                    campaign_id=CAMP_ID,
                    sample_unit_id="LOT31_TP53",
                    controls=_vaf_ctrl_l2,
                )
            )
        if "obs_castseq_translocation_freq" not in extra_eids and meas_use.get("obs_castseq_translocation_freq") is not None:
            _cast_ctrl_l2 = {}
            if meas_use.get("obs_castseq_duplex_or_amplicon_depth") is not None:
                _cast_ctrl_l2["obs_castseq_duplex_or_amplicon_depth"] = meas_use["obs_castseq_duplex_or_amplicon_depth"]
            s9_lot_atoms_use.append(
                wc_mod.AtomicObservation(
                    endpoint_id="obs_castseq_translocation_freq",
                    scope_type="LOT_SCOPE",
                    scope_key=lot_k_use,
                    kind="CONTINUOUS_STUDENT_T",
                    value_or_k=meas_use["obs_castseq_translocation_freq"],
                    se_or_n=meas_use["obs_castseq_translocation_se"],
                    sample_n=_cast_n_l2,
                    campaign_id=CAMP_ID,
                    sample_unit_id="LOT31_CAST",
                    controls=_cast_ctrl_l2,
                )
            )
        s11_atoms_use = []
        if include_s11_key or s11_record is not None:
            sk_s11 = mcb_k_use if s11_ledger_scope == "MCB_SCOPE" else (dev_k_use if s11_ledger_scope == "ASSEMBLED_DEVICE_LOT_SCOPE" else lot_k_use)
            s11_atoms_use.append(
                wc_mod.AtomicObservation(
                    endpoint_id="stage0_locus2_modeb_v3_panel",
                    scope_type=s11_ledger_scope,
                    scope_key=sk_s11,
                    kind="STRUCTURED_PANEL",
                    value_or_k=s11_record,
                    campaign_id=CAMP_ID,
                    sample_unit_id=str(scope_use.get("sample_id", scope_use.get("sample_unit_id", "S11_PANEL_01"))),
                )
            )
        # These inputs contain lot evidence only. Unrelated bank and device
        # endpoints stay missing on each path.
        batches_l2 = [{"batch_id": "B_BASE", "atoms": s9_lot_atoms_use + s11_atoms_use + (extra_ledger_atoms or [])}]
        for idx_a, raw_r in enumerate(t6b_records):
            if not isinstance(raw_r, dict):
                batches_l2.append({"batch_id": f"B_T6_{idx_a}", "atoms": [raw_r]})
                continue
            try:
                a_obj = _to_atom(raw_r)
                batches_l2.append({"batch_id": f"B_T6_{idx_a}", "atoms": [a_obj]})
            except (ValueError, TypeError):
                batches_l2.append({"batch_id": f"B_T6_RAW_{idx_a}", "atoms": [], "tier6b_sorted_lot_observations": [raw_r]})
        if extra_ledger_batches:
            batches_l2.extend(extra_ledger_batches)
        p2 = wc_mod.evaluate_iip0024_atomic_stage0_ledger(batches_l2, tprov_use, target_campaign_id=CAMP_ID, root=ROOT)

        # Path 3: Serialized JSON Ingestion
        doc3 = {
            "scope": dict(scope_use),
            "measurements": dict(meas_use),
            "iip0026_molecular_observations": dict(base_mol_31),
            "tier6b_sorted_lot_observations": list(t6b_records),
            **s11_dict_part,
        }
        p3 = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(doc3), root=ROOT)

        # Path 4: Chronological Replay
        replay_recs = []
        if split_replay_batches and len(t6b_records) >= 2:
            replay_recs.append({
                "utc_timestamp": "2026-11-19T01:00:00Z",
                "seq_index": 1,
                "scope": dict(scope_use),
                "measurements": dict(meas_use),
                "iip0026_molecular_observations": dict(base_mol_31),
                "tier6b_sorted_lot_observations": [t6b_records[0]],
                **s11_dict_part,
            })
            for r_i, item in enumerate(t6b_records[1:], start=2):
                if isinstance(item, dict):
                    replay_recs.append({
                        "utc_timestamp": f"2026-11-19T0{min(r_i, 9)}:00:00Z",
                        "seq_index": r_i,
                        "scope": dict(scope_use),
                        "tier6b_sorted_lot_observations": [item],
                    })
                else:
                    replay_recs.append({"utc_timestamp": f"2026-11-19T0{min(r_i, 9)}:00:00Z", "seq_index": r_i, "malformed_record": True})
        else:
            replay_recs.append({
                "utc_timestamp": "2026-11-19T01:00:00Z",
                "seq_index": 1,
                "scope": dict(scope_use),
                "measurements": dict(meas_use),
                "iip0026_molecular_observations": dict(base_mol_31),
                "tier6b_sorted_lot_observations": list(t6b_records),
                **s11_dict_part,
            })
        if extra_replay_records:
            replay_recs.extend(extra_replay_records)
        p4 = ing_mod.replay_stage0_chronological_batches(replay_recs, root=ROOT)
        return {"direct": p1, "ledger": p2, "ingest": p3, "replay": p4}

    # --- FIXTURE 1: Raw-count pass across all 4 entry points (unconditional 3-way simultaneous Clopper-Pearson) + exact 2-way math_only benchmark ---
    m_fix1_2way = wc_mod.compute_tier6b_per_fraction_thinning_bound_math_only([
        {"fraction_id": "enriched", "N": 30_000_000, "w_i": 0.9999, "d_i_L": wc_mod._exact_one_sided_cp_lower_iip0031(20, 20, 0.05 / 8.0)},
        {"fraction_id": "flow_through", "N": 100_000, "w_i": 0.0001, "d_i_L": wc_mod._exact_one_sided_cp_lower_iip0031(10, 10, 0.05 / 8.0)},
    ])
    assert m_fix1_2way["release_clearance_permitted"] is False and abs(m_fix1_2way["p_lot_U95"] - 1.9552099522e-7) < 1e-15
    assert abs(wc_mod._exact_one_sided_cp_lower_iip0031(20, 20, 0.05 / 8.0) - 0.7758790062) < 1e-9
    assert abs(wc_mod._exact_one_sided_cp_lower_iip0031(10, 10, 0.05 / 8.0) - 0.6019882323) < 1e-9
    m_fix2_2way = wc_mod.compute_tier6b_per_fraction_thinning_bound_math_only([
        {"fraction_id": "enriched", "N": 9_999_500, "w_i": 0.9999, "d_i_L": wc_mod._exact_one_sided_cp_lower_iip0031(20, 20, 0.05 / 8.0)},
        {"fraction_id": "flow_through", "N": 1_000, "w_i": 0.0001, "d_i_L": wc_mod._exact_one_sided_cp_lower_iip0031(10, 10, 0.05 / 8.0)},
    ])
    assert m_fix2_2way["within_numerical_cap"] is False and abs(m_fix2_2way["p_lot_U95"] - 1.2910869131e-6) < 1e-15
    assert abs(wc_mod._exact_one_sided_cp_lower_iip0031(7, 7, 0.05 / 8.0) - 0.4843125430) < 1e-9
    assert abs(wc_mod._exact_one_sided_cp_lower_iip0031(8, 8, 0.05 / 8.0) - 0.5302552806) < 1e-9

    fix1_recs = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F1_ENRICH", k=0, N=30_000_000, w_i=0.9999, hpsc=(20, 20), diff=(20, 20)),
        _mk_frac("flow_through_and_excluded", "F1_FT", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    res1 = _eval_4paths(fix1_recs)
    for nm, out in res1.items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["disposition"] == "PASS", (nm, ev)
        assert abs(ev["d_L_by_fraction"]["enriched_cd31_neg_cd144_neg"] - 0.7603077860) < 1e-9, (nm, ev)
        assert abs(ev["d_L_by_fraction"]["flow_through_and_excluded"] - 0.5780679294) < 1e-9, (nm, ev)
        assert ev["p_lot_U95"] <= 3.0e-7 and abs(ev["p_lot_U95"] - 2.1870701159e-7) < 1e-11, (nm, ev)
    assert res1["direct"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED"
    assert res1["ingest"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED"
    assert res1["replay"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED"
    assert res1["ledger"]["readiness_tier"] == "PROVISIONAL_PARTIAL_LOT_SCREENING_ONLY"
    assert set(res1["ledger"]["holds"]) == set(wc_mod.ENDPOINT_REGISTRY_SPEC["mcb_mandatory_endpoint_ids"])

    # --- FIXTURE 2: Zero counts, insufficient evidence across all 4 entry points (including single-dose ~80k scale) ---
    fix2_recs = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F2_ENRICH", k=0, N=9_999_500, w_i=0.9999, hpsc=(20, 20), diff=(20, 20)),
        _mk_frac("flow_through_and_excluded", "F2_FT", k=0, N=1_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    res2 = _eval_4paths(fix2_recs)
    for nm, out in res2.items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, ev)
        assert ev["p_lot_U95"] > 1.4e-6 and abs(ev["p_lot_U95"] - 1.5112233424e-6) < 1e-10, (nm, ev)
        assert out["overall_status"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, out)
    # Single-graft-dose (~80,000 total clonogenic opportunities across 2 fractions) fails closed on all 4 paths (> 5.98419e-5 single-fraction alpha=0.05/6 floor):
    fix2_single_dose_80k = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F2_80K_1", k=0, N=79_000, w_i=79_000 / 80_000, hpsc=(20, 20), diff=(20, 20)),
        _mk_frac("flow_through_and_excluded", "F2_80K_2", k=0, N=1_000, w_i=1_000 / 80_000, hpsc=(20, 20), diff=(20, 20)),
    ]
    for nm, out in _eval_4paths(fix2_single_dose_80k).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, out)
        assert out["tier6b_sorted_lot_evaluation"]["p_lot_U95"] > 5.98419e-5, (nm, out)
        assert out["tier6b_sorted_lot_evaluation"]["p_lot_U95"] > 1.19e-4, (nm, out)

    # --- FIXTURE 3: Detection floor & separate calibrants across all 4 entry points ---
    fix3_7 = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F3_7_ENRICH", k=0, N=9_999_500, w_i=0.9999, hpsc=(20, 20), diff=(7, 7)),
        _mk_frac("flow_through_and_excluded", "F3_7_FT", k=0, N=1_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    res3_7 = _eval_4paths(fix3_7)
    for nm, out in res3_7.items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, ev)
        assert ev["detection_floor_cleared"] is False
        assert abs(ev["d_L_by_fraction"]["enriched_cd31_neg_cd144_neg"] - 0.4570564295) < 1e-9, (nm, ev)
        assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
    # Pooling 27/27 must NOT clear:
    fix3_pool = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F3_P_ENRICH", k=0, N=30_000_000, w_i=0.9999, hpsc=(27, 27), diff=(27, 27), extra={"pooled_calibrants": True}),
        _mk_frac("flow_through_and_excluded", "F3_P_FT", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for nm, out in _eval_4paths(fix3_pool).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
    # Changing 7/7 -> 8/8 clears detection floor ONLY (0.5040499718 >= 0.50), still holds on N:
    fix3_8 = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F3_8_ENRICH", k=0, N=9_999_500, w_i=0.9999, hpsc=(20, 20), diff=(8, 8)),
        _mk_frac("flow_through_and_excluded", "F3_8_FT", k=0, N=1_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    res3_8 = _eval_4paths(fix3_8)
    for nm, out in res3_8.items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["detection_floor_cleared"] is True, (nm, ev)
        assert abs(ev["d_L_by_fraction"]["enriched_cd31_neg_cd144_neg"] - 0.5040499718) < 1e-9, (nm, ev)
        assert ev["disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, ev)
        assert out["overall_status"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, out)

    # --- FIXTURE 4: Pure-math unrounded 29,213,508 vs 29,213,509 boundary vs live rejection of supplied d_L ---
    m_hold = wc_mod.compute_tier6b_per_fraction_thinning_bound_math_only([
        {"fraction_id": "enriched", "N": 29_213_508, "w_i": 0.5, "d_i_L": 0.50},
        {"fraction_id": "flow_through", "N": 29_213_508, "w_i": 0.5, "d_i_L": 0.50},
    ])
    assert m_hold["math_only_non_authorizing"] is True and m_hold["release_clearance_permitted"] is False
    assert m_hold["within_numerical_cap"] is False and abs(m_hold["p_lot_U95"] - 3.0000000725841394e-7) < 1e-18
    m_pass = wc_mod.compute_tier6b_per_fraction_thinning_bound_math_only([
        {"fraction_id": "enriched", "N": 29_213_509, "w_i": 0.5, "d_i_L": 0.50},
        {"fraction_id": "flow_through", "N": 29_213_509, "w_i": 0.5, "d_i_L": 0.50},
    ])
    assert m_pass["math_only_non_authorizing"] is True and m_pass["release_clearance_permitted"] is False
    assert m_pass["within_numerical_cap"] is True and abs(m_pass["p_lot_U95"] - 2.9999999698919284e-7) < 1e-18

    # COUNTEREXAMPLE 1 (Sol R6): Live paths MUST REJECT calculation_only_supplied_d_L / supplied_d_i_L (absent OR present calibrants)!
    for sup_lbl, sup_recs in (
        (
            "calibrants_absent_supplied_d_L",
            [
                _mk_frac("enriched_cd31_neg_cd144_neg", "CE1_A1", k=0, N=29_213_509, w_i=0.5, hpsc=None, diff=None, extra={"calculation_only_supplied_d_L": True, "supplied_d_i_L": 0.50}),
                _mk_frac("flow_through_and_excluded", "CE1_A2", k=0, N=29_213_509, w_i=0.5, hpsc=None, diff=None, extra={"calculation_only_supplied_d_L": True, "supplied_d_i_L": 0.50}),
            ],
        ),
        (
            "calibrants_present_but_supplied_d_L_flag",
            [
                _mk_frac("enriched_cd31_neg_cd144_neg", "CE1_B1", k=0, N=29_213_509, w_i=0.5, extra={"calculation_only_supplied_d_L": True, "supplied_d_i_L": 0.50}),
                _mk_frac("flow_through_and_excluded", "CE1_B2", k=0, N=29_213_509, w_i=0.5, extra={"calculation_only_supplied_d_L": True, "supplied_d_i_L": 0.50}),
            ],
        ),
    ):
        for nm, out in _eval_4paths(sup_recs).items():
            assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (sup_lbl, nm, out)
            assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (sup_lbl, nm, out)

    # --- FIXTURE 5 & COUNTEREXAMPLES 2, 3, 6 (Missing/invalid evidence, single_fraction_unsorted_lot, order-invariance, CFU & mass-balance) ---
    obs_15m_cal = _mk_frac("enriched_cd31_neg_cd144_neg", "CE3_CAL_15M", k=0, N=15_000_000, w_i=0.9999, hpsc=(20, 20), diff=(20, 20))
    obs_15m_uncal = _mk_frac("enriched_cd31_neg_cd144_neg", "CE3_UNCAL_15M", k=0, N=15_000_000, w_i=0.9999, hpsc=None, diff=None)
    obs_15m_low7 = _mk_frac("enriched_cd31_neg_cd144_neg", "CE3_LOW7_15M", k=0, N=15_000_000, w_i=0.9999, hpsc=(20, 20), diff=(7, 7))
    obs_ft_ok = _mk_frac("flow_through_and_excluded", "CE3_FT_OK", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10))

    fix5_cases = {
        "omitted_flow_through": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_NOFT", k=0, N=30_000_000, w_i=1.0, roster=("enriched_cd31_neg_cd144_neg",)),
        ],
        "ce2_single_fraction_unsorted_lot_bypass_attempt": [
            _mk_frac(
                "enriched_cd31_neg_cd144_neg",
                "CE2_UNSORT_BYPASS",
                k=0,
                N=30_000_000,
                w_i=1.0,
                roster=("enriched_cd31_neg_cd144_neg",),
                extra={"single_fraction_unsorted_lot": True},
            ),
        ],
        "ce3_partial_cal_calibrated_first": [obs_15m_cal, obs_15m_uncal, obs_ft_ok],
        "ce3_partial_cal_uncalibrated_first": [obs_15m_uncal, obs_15m_cal, obs_ft_ok],
        "ce3_mixed_7of7_cal_order_ab": [obs_15m_cal, obs_15m_low7, obs_ft_ok],
        "ce3_mixed_7of7_cal_order_ba": [obs_15m_low7, obs_15m_cal, obs_ft_ok],
        "unreconciled_losses": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_LOSS1", k=0, N=30_000_000, w_i=0.9999, extra={"unreconciled_loss_fraction": 0.02}),
            _mk_frac("flow_through_and_excluded", "F5_LOSS2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "unreconciled_integer_presort_cell_balance": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_IBAL1", k=0, N=30_000_000, w_i=0.9999, extra={"presort_fraction_cells": 2_900_000_000}),
            _mk_frac("flow_through_and_excluded", "F5_IBAL2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "plated_exceeds_fraction_physical_cells": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_PL1", k=0, N=30_000_000, w_i=0.9999),
            _mk_frac("flow_through_and_excluded", "F5_PL2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), extra={"plated_treated_cells_n": 400_000}),
        ],
        "missing_untreated_cfu_control": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_UCFU1", k=0, N=30_000_000, w_i=0.9999, extra={"untreated_cfu_control": None}),
            _mk_frac("flow_through_and_excluded", "F5_UCFU2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "non_disjoint_untreated_cfu_sample_id": [
            _mk_frac(
                "enriched_cd31_neg_cd144_neg",
                "F5_DISJ1",
                k=0,
                N=30_000_000,
                w_i=0.9999,
                extra={"untreated_cfu_control": {"untreated_sample_id": "F5_DISJ2", "untreated_cfu_colonies_k": 100, "untreated_plated_cells_n": 500}},
            ),
            _mk_frac("flow_through_and_excluded", "F5_DISJ2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "omitted_endogenous_denominator": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_ENDO1", k=0, N=30_000_000, w_i=0.9999, extra={"endogenous_denominator_verified": False}),
            _mk_frac("flow_through_and_excluded", "F5_ENDO2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "omitted_calibrant": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_CAL1", k=0, N=30_000_000, w_i=0.9999, diff=None),
            _mk_frac("flow_through_and_excluded", "F5_CAL2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "invalid_weights": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_W1", k=0, N=30_000_000, w_i=0.70),
            _mk_frac("flow_through_and_excluded", "F5_W2", k=0, N=100_000, w_i=0.20, hpsc=(10, 10), diff=(10, 10)),
        ],
        "zero_N": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_N0_1", k=0, N=0, w_i=0.9999),
            _mk_frac("flow_through_and_excluded", "F5_N0_2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "boolean_count": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_BOOL1", k=False, N=30_000_000, w_i=0.9999),
            _mk_frac("flow_through_and_excluded", "F5_BOOL2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "float_count": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_FLT1", k=0.0, N=30_000_000, w_i=0.9999),
            _mk_frac("flow_through_and_excluded", "F5_FLT2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "clustered_sampling": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F5_CL1", k=0, N=30_000_000, w_i=0.9999, extra={"clustered_sampling": True}),
            _mk_frac("flow_through_and_excluded", "F5_CL2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
    }
    res_order_cache: dict[str, dict] = {}
    for sub_name, sub_recs in fix5_cases.items():
        res_sub = _eval_4paths(sub_recs)
        res_order_cache[sub_name] = res_sub
        for nm, out in res_sub.items():
            assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (sub_name, nm, out)
            assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (sub_name, nm, out)
    # Verify exact order invariance on reversed sub-observation sequences:
    for nm in ("direct", "ledger", "ingest", "replay"):
        assert res_order_cache["ce3_partial_cal_calibrated_first"][nm]["tier6b_reason"] == res_order_cache["ce3_partial_cal_uncalibrated_first"][nm]["tier6b_reason"]
        assert res_order_cache["ce3_mixed_7of7_cal_order_ab"][nm]["tier6b_sorted_lot_evaluation"]["d_L_by_fraction"] == res_order_cache["ce3_mixed_7of7_cal_order_ba"][nm]["tier6b_sorted_lot_evaluation"]["d_L_by_fraction"]

    # Also verify evaluate_rare_event_count_iip0024 NEVER defaults missing eta_spike_recovery_l95 / f_clonogenic to 1.0:
    re_miss = wc_mod.evaluate_rare_event_count_iip0024(
        gate_id="TIER6B_SORTED_CLONOGENIC_CHECK",
        k=0,
        n=30_000_000,
        max_u95_fraction=3.0e-7,
        require_observed_zero=True,
        eta_spike_recovery_l95=None,
        f_clonogenic=None,
    )
    assert re_miss["status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION"
    re_pos_uncal = wc_mod.evaluate_rare_event_count_iip0024(
        gate_id="TIER6B_SORTED_CLONOGENIC_CHECK",
        k=1,
        n=0,
        max_u95_fraction=3.0e-7,
        require_observed_zero=True,
        eta_spike_recovery_l95=None,
        f_clonogenic=None,
    )
    assert re_pos_uncal["status"] == "FAIL_HARD_STOP"

    # --- FIXTURE 6 & COUNTEREXAMPLES 4, 5 (Positive outgrowth with malformed controls + S9 VAF 0.0002 parity across all 4 paths) ---
    fix6_valid_cal = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F6_POS_V", k=1, N=30_000_000, w_i=0.9999, hpsc=(20, 20), diff=(20, 20)),
        "MALFORMED_NON_DICT_RECORD",
        _mk_frac("enriched_cd31_neg_cd144_neg", "F6_Z1", k=0, N=30_000_000, w_i=0.9999, hpsc=(20, 20), diff=(20, 20)),
        _mk_frac("flow_through_and_excluded", "F6_Z2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    fix6_invalid_cal = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "F6_POS_BAD", k=1, N=0, w_i=None, hpsc=None, diff=None, extra={"endogenous_denominator_verified": False}),
        "MALFORMED_NON_DICT_RECORD",
        _mk_frac("enriched_cd31_neg_cd144_neg", "F6_Z3", k=0, N=30_000_000, w_i=0.9999, hpsc=(20, 20), diff=(20, 20)),
        _mk_frac("flow_through_and_excluded", "F6_Z4", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    # COUNTEREXAMPLE 4 (Sol R6): Same-lot k=1 with valid outer identity and malformed ancillary `controls`, followed by valid zeros!
    ce4_malformed_controls_pos = [
        {
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "product_variant_id": "ISLET1_A4_ARM1B",
            "construct_sha256": a4_sha,
            "donor_seedline_id": SEED_ID,
            "mcb_id": MCB_ID,
            "lot_id": LOT_ID,
            "fraction_id": "enriched_cd31_neg_cd144_neg",
            "sample_id": "CE4_POS_MALCTRL",
            "k": 1,
            "N": 30_000_000,
            "controls": "MALFORMED_ANCILLARY_CONTROLS_STRING",
        },
        *fix1_recs,
    ]
    for label6, recs6 in (
        ("valid_cal", fix6_valid_cal),
        ("missing_invalid_cal", fix6_invalid_cal),
        ("ce4_malformed_controls_on_positive", ce4_malformed_controls_pos),
    ):
        res6 = _eval_4paths(recs6)
        for nm, out in res6.items():
            assert out["tier6b_disposition"] == "FAIL_HARD_STOP", (label6, nm, out)
            assert out["tier6b_reason"] == "FAIL_TIER6B_OUTGROWTH_DETECTED", (label6, nm, out)
            assert out["overall_status"] == "FAIL_HARD_STOP", (label6, nm, out)
            assert out["hard_stop_precedence_winner"] == "FAIL_TIER6B_OUTGROWTH_DETECTED", (label6, nm, out)
    # Also verify raw dict inside ledger batch["atoms"] (fallback non-AtomicObservation lane) latches FAIL_HARD_STOP on k=1:
    led_raw_atoms_pos = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
        [
            {"batch_id": "B_BASE", "atoms": mcb_atoms_base + dev_atoms_base},
            {"batch_id": "B_RAW_ATOM_LANE", "atoms": [ce4_malformed_controls_pos[0]]},
            {"batch_id": "B_CLEAN_ZEROS", "atoms": [_to_atom(r) for r in fix1_recs]},
        ],
        target_prov,
        target_campaign_id=CAMP_ID,
        root=ROOT,
    )
    assert led_raw_atoms_pos["overall_status"] == "FAIL_HARD_STOP"
    assert led_raw_atoms_pos["hard_stop_precedence_winner"] == "FAIL_TIER6B_OUTGROWTH_DETECTED"

    # Independent S9 hard stop keeps FIRST precedence over Tier 6B k=1 across all 4 entry points:
    s9_atom_fail = wc_mod.AtomicObservation(
        endpoint_id="obs_icasp9_log10_kill",
        scope_type="LOT_SCOPE",
        scope_key=lot_k,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=2.0,
        se_or_n=0.1,
        sample_n=6,
        campaign_id=CAMP_ID,
        sample_unit_id="S9_KILL_FAIL",
    )
    res6_s9 = _eval_4paths(
        fix6_valid_cal,
        meas_override={"obs_icasp9_log10_kill": 2.0, "obs_icasp9_log10_kill_se": 0.1},
        extra_ledger_atoms=[s9_atom_fail],
    )
    for nm, out in res6_s9.items():
        assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (nm, out)
        assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (nm, out)

    # COUNTEREXAMPLE 5 (Sol R6): Qualifying Tier 6B with TP53/PPM1D/BCOR VAF mean 0.0002, SE 0.000001, depth 70,000
    # MUST hard-stop under S9 (>= 1e-4) across ALL FOUR routes INCLUDING the atomic ledger!
    s9_vaf_0p0002_atom = wc_mod.AtomicObservation(
        endpoint_id="obs_tp53_ppm1d_bcor_vaf",
        scope_type="LOT_SCOPE",
        scope_key=lot_k,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=0.0002,
        se_or_n=0.000001,
        sample_n=6,
        campaign_id=CAMP_ID,
        sample_unit_id="S9_VAF_2E4",
        controls={"obs_tp53_ppm1d_bcor_duplex_depth": 70000},
    )
    res_ce5_vaf = _eval_4paths(
        fix1_recs,
        meas_override={
            "obs_tp53_ppm1d_bcor_vaf": 0.0002,
            "obs_tp53_ppm1d_bcor_vaf_se": 0.000001,
            "obs_tp53_ppm1d_bcor_duplex_depth": 70000,
        },
        extra_ledger_atoms=[s9_vaf_0p0002_atom],
    )
    for nm, out in res_ce5_vaf.items():
        assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (nm, out)
        assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (nm, out)
    # Underpowered VAF depth (< 62,958) at VAF=0.0 holds on the atomic ledger:
    s9_vaf_underdepth_atom = wc_mod.AtomicObservation(
        endpoint_id="obs_tp53_ppm1d_bcor_vaf",
        scope_type="LOT_SCOPE",
        scope_key=lot_k,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=0.0,
        se_or_n=1e-6,
        sample_n=6,
        campaign_id=CAMP_ID,
        sample_unit_id="S9_VAF_LOWDEPTH",
        controls={"obs_tp53_ppm1d_bcor_duplex_depth": 30000},
    )
    led_lowdepth = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
        [{"batch_id": "B1", "atoms": mcb_atoms_base + dev_atoms_base + [s9_vaf_underdepth_atom] + [_to_atom(r) for r in fix1_recs]}],
        target_prov,
        target_campaign_id=CAMP_ID,
        root=ROOT,
    )
    assert led_lowdepth["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE"

    # --- FIXTURE 7 & COUNTEREXAMPLE 6 (Scope, omitted outer provenance, mixed fraction identities, sample_id dedup, cross-lot non-latch) ---
    fix7_rejects = {
        "wrong_lot": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F7_WL1", lot_id="OTHER_LOT_99"),
            _mk_frac("flow_through_and_excluded", "F7_WL2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), lot_id="OTHER_LOT_99"),
        ],
        "wrong_construct": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F7_WC1", construct_sha256="0" * 64),
            _mk_frac("flow_through_and_excluded", "F7_WC2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), construct_sha256="0" * 64),
        ],
        "wrong_bank": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F7_WB1", mcb_id="MCB_OTHER_999"),
            _mk_frac("flow_through_and_excluded", "F7_WB2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), mcb_id="MCB_OTHER_999"),
        ],
        "ce6_inconsistent_identities_across_two_fractions": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "CE6_MIX1", construct_sha256=a4_sha, mcb_id=MCB_ID),
            _mk_frac("flow_through_and_excluded", "CE6_MIX2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), construct_sha256="b" * 64, mcb_id="MCB_OTHER_MIX"),
        ],
        "metadata_instead_of_counts": [
            _mk_frac("enriched_cd31_neg_cd144_neg", "F7_MD1", kind="METADATA_CONTROL"),
            _mk_frac("flow_through_and_excluded", "F7_MD2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ],
        "pure_hpsc_mcb_instead_of_sorted_lot": [
            _mk_frac(
                "enriched_cd31_neg_cd144_neg",
                "F7_MCB1",
                N=900_000_000,
                scope_type="MCB_SCOPE",
                endpoint_id="mcb_tier6b_pure_hpsc_clonogenic_escape",
                extra={"eta_spike_recovery_l95": 0.80, "f_clonogenic": 0.80},
            ),
        ],
        "assembled_product_instead_of_sorted_lot": [
            _mk_frac(
                "enriched_cd31_neg_cd144_neg",
                "F7_DEV1",
                N=300_000_000,
                scope_type="ASSEMBLED_DEVICE_LOT_SCOPE",
                endpoint_id="gate0s_shedding_viable_outgrowth_assay",
                extra={"eta_spike_recovery_l95": 0.80, "f_clonogenic": 0.80},
            ),
        ],
    }
    for sub7, recs7 in fix7_rejects.items():
        for nm, out in _eval_4paths(recs7).items():
            assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (sub7, nm, out)
            assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (sub7, nm, out)

    # COUNTEREXAMPLE 6b (Sol R6): Omitting outer construct_sha256, donor_seedline_id, and mcb_id MUST hold on ALL 4 paths!
    mixed_two_fracs = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "CE6B_1", construct_sha256="a" * 64, donor_seedline_id="SEED_A", mcb_id="MCB_A"),
        _mk_frac("flow_through_and_excluded", "CE6B_2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), construct_sha256="b" * 64, donor_seedline_id="SEED_B", mcb_id="MCB_B"),
    ]
    omitted_outer_scope = {"construct_sha256": None, "donor_seedline_id": None, "mcb_id": None}
    for recs_om in (fix1_recs, mixed_two_fracs):
        p1_om = wc_mod.evaluate_stage0_product_contract(
            {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": LOT_ID, "campaign_id": CAMP_ID, **base_meas_31, **base_mol_31, "tier6b_sorted_lot_observations": recs_om},
            root=ROOT,
        )
        p2_om = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
            [{"batch_id": "B1", "atoms": [_to_atom(r) for r in recs_om]}],
            {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": LOT_ID, "device_assembly_id": DEV_ID, "physical_route_id": ROUTE_ID},
            target_campaign_id=CAMP_ID,
            root=ROOT,
        )
        p3_om = ing_mod.ingest_stage0_research_or_wetlab_batch(
            json.dumps({
                "scope": {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": LOT_ID, "campaign_id": CAMP_ID},
                "measurements": dict(base_meas_31),
                "iip0026_molecular_observations": dict(base_mol_31),
                "tier6b_sorted_lot_observations": recs_om,
            }),
            root=ROOT,
        )
        p4_om = ing_mod.replay_stage0_chronological_batches(
            [{
                "utc_timestamp": "2026-11-19T01:00:00Z",
                "seq_index": 1,
                "scope": {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": LOT_ID, "campaign_id": CAMP_ID},
                "measurements": dict(base_meas_31),
                "iip0026_molecular_observations": dict(base_mol_31),
                "tier6b_sorted_lot_observations": recs_om,
            }],
            root=ROOT,
        )
        assert p1_om["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION" and p1_om["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", p1_om
        assert p3_om["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION" and p3_om["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", p3_om
        assert p4_om["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION" and p4_om["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", p4_om
        assert p2_om["readiness_tier"] == "HOLD_MISSING_OR_INVALID_IMMUTABLE_PROVENANCE_KEY", p2_om
        assert omitted_outer_scope["construct_sha256"] is None

    # 7b. Replayed sample_id under a new record_id MUST NOT inflate N_i:
    fix7_dup = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "DUP_SID_01", k=0, N=15_000_000, w_i=0.9999, extra={"record_id": "REC_A"}),
        _mk_frac("enriched_cd31_neg_cd144_neg", "DUP_SID_01", k=0, N=15_000_000, w_i=0.9999, extra={"record_id": "REC_B"}),
        _mk_frac("flow_through_and_excluded", "F7_FT_DUP", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    res7_dup = _eval_4paths(fix7_dup)
    for nm, out in res7_dup.items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["N_by_fraction"]["enriched_cd31_neg_cd144_neg"] == 15_000_000, (nm, ev)
        assert "DUP_SID_01" in ev["deduplicated_sample_ids"], (nm, ev)
        assert ev["disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, ev)
    # Distinct calibrated sample_ids with independent source-owned trial & calibrant IDs summing to 30,000,000 DO pass:
    fix7_distinct = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "DIST_SID_01", k=0, N=15_000_000, w_i=0.9999, extra={"record_id": "DIST_REC_01", "source_sample_id": "DIST_SRC_01", "cal_batch_id": "DIST_CAL_01"}),
        _mk_frac("enriched_cd31_neg_cd144_neg", "DIST_SID_02", k=0, N=15_000_000, w_i=0.9999, extra={"record_id": "DIST_REC_02", "source_sample_id": "DIST_SRC_02", "cal_batch_id": "DIST_CAL_02"}),
        _mk_frac("flow_through_and_excluded", "F7_FT_DIST", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for nm, out in _eval_4paths(fix7_distinct).items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["N_by_fraction"]["enriched_cd31_neg_cd144_neg"] == 30_000_000, (nm, ev)
        assert ev["disposition"] == "PASS", (nm, ev)
    # Without source-owned record_id / source_sample_id / cal_batch_id, multi-observation summing MUST hold (#5):
    fix7_missing_src_id = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "DIST_NOID_01", k=0, N=15_000_000, w_i=0.9999),
        _mk_frac("enriched_cd31_neg_cd144_neg", "DIST_NOID_02", k=0, N=15_000_000, w_i=0.9999),
        _mk_frac("flow_through_and_excluded", "F7_FT_NOID", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for nm, out in _eval_4paths(fix7_missing_src_id).items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, ev)
        assert ev["reason_code"] == "MISSING_SOURCE_OWNED_TRIAL_OR_CALIBRANT_IDENTITY_FOR_MULTI_OBSERVATION_SUM", (nm, ev)

    # 7c. Adverse observation (k=1) from a DIFFERENT lot does NOT latch target lot across all 4 entry points:
    fix7_other_lot_pos = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "OTHER_POS_1", k=1, N=30_000_000, w_i=0.9999, lot_id="OTHER_LOT_99"),
        *fix1_recs,
    ]
    res7_ol = _eval_4paths(fix7_other_lot_pos)
    for nm, out in res7_ol.items():
        assert out["tier6b_disposition"] == "PASS", (nm, out)
    assert res7_ol["direct"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED"
    assert res7_ol["ledger"]["readiness_tier"] == "PROVISIONAL_PARTIAL_LOT_SCREENING_ONLY"

    # --- RULE S11 QUANTITATIVE PANEL, NON-EXECUTABLE L2_V3_* PRECEDENCE, & SOL R6 S11 COUNTEREXAMPLES ---
    co_cd200_def = by_id["ISLET1_A4_L2_V3_CO_PRIMARY_CD200"]
    co_cd200_sha = co_cd200_def["staged_locus2_modeb_donor_binding"]["locus2_donor_sha256"]
    s11_passing_panel = {
        "scope_type": "LOT_SCOPE",
        "sample_id": "S11_PANEL_01",
        "ciita_locus2_genotype_copy_class": "KI/KO",
        "obs_stage0_independent_lot_replicates_n": 6,
        "phasing_verified": True,
        "cis_verified": True,
        "wt_inducibility_verified": True,
        "rnai_attribution_verified": True,
        "ciita_ki_allele_informative_transcript_recovery_verified": True,
        "ciita_ko_allele_informative_transcript_recovery_verified": True,
        "ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified": True,
        "lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified": True,
        "l2_v3_same_background_intronless_control_linked_to_cargo": True,
        "l2_v3_same_background_scramble_control_linked_to_cargo": True,
        "obs_l2_v3_chr16_ha_wgs_pass": True,
        "obs_l2_v3_chr16_ha_wgs_variant_count": 0,
        "obs_l2_v3_tandem_polya_readthrough_fold": 1.15, "obs_l2_v3_tandem_polya_readthrough_se": 0.10,
        "obs_l2_v3_cfh_stp_gpi_quantibrite_molecules": 85000.0, "obs_l2_v3_cfh_stp_gpi_quantibrite_se": 4000.0,
        "obs_l2_v3_unit2_copy_expression_quantibrite_or_fold": 45000.0, "obs_l2_v3_unit2_copy_expression_se": 2500.0,
        "obs_ciita_ki_allele_informative_reads_n": 500,
        "obs_ciita_ko_allele_informative_reads_n": 500,
        "obs_ciita_ki_ex4_ex6_aberrant_skip_reads_k": 0,
        "obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k": 0,
        "obs_ciita_ki_ex4_ex6_aberrant_skip_pct": 0.0, "obs_ciita_ki_ex4_ex6_aberrant_skip_se": 0.05,
        "obs_ciita_ko_ex4_ex6_aberrant_skip_pct": 0.0, "obs_ciita_ko_ex4_ex6_aberrant_skip_se": 0.05,
        "obs_edited_72h_ifng_hla_dr_dp_dq_pos_pct": 0.2, "obs_edited_72h_ifng_hla_dr_dp_dq_pos_se": 0.05,
        "obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_pct": 75.0, "obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_se": 2.0,
        "obs_cis_socs1_tost90_ratio_lo": 0.92, "obs_cis_socs1_tost90_ratio_hi": 1.08,
        "obs_cis_dexi_tost90_ratio_lo": 0.90, "obs_cis_dexi_tost90_ratio_hi": 1.10,
        "obs_cis_clec16a_tost90_ratio_lo": 0.91, "obs_cis_clec16a_tost90_ratio_hi": 1.09,
        "obs_16p13_ddpcr_cn_ratio": 2.00,
        "obs_16p13_phased_baf": 0.50,
        "obs_l2_v3_unit1_eef1a1_resting_0h_5race_umi_n_eff": 4200, "obs_l2_v3_unit1_eef1a1_resting_0h_unexpected_tss_umi_k": 0,
        "obs_l2_v3_unit1_eef1a1_stim_24h_5race_umi_n_eff": 4300, "obs_l2_v3_unit1_eef1a1_stim_24h_unexpected_tss_umi_k": 0,
        "obs_l2_v3_unit2_resting_0h_5race_umi_n_eff": 4100, "obs_l2_v3_unit2_resting_0h_unexpected_tss_umi_k": 0,
        "obs_l2_v3_unit2_stim_24h_5race_umi_n_eff": 4500, "obs_l2_v3_unit2_stim_24h_unexpected_tss_umi_k": 0,
        "obs_l2_v3_unit2_tss_pos1_in_donor": 665,
        "obs_l2_v3_stage0_eef1a1_cpg_methylation_pct": 3.5, "obs_l2_v3_stage0_eef1a1_cpg_methylation_se": 0.4,
        "obs_l2_v3_stage0_chs4_cpg_methylation_pct": 4.2, "obs_l2_v3_stage0_chs4_cpg_methylation_se": 0.5,
        "obs_l2_v3_stage0_dm1_cpg_methylation_pct": 4.8, "obs_l2_v3_stage0_dm1_cpg_methylation_se": 0.5,
        "obs_l2_v3_stage7_eef1a1_cpg_methylation_pct": 4.0, "obs_l2_v3_stage7_eef1a1_cpg_methylation_se": 0.5,
        "obs_l2_v3_stage7_chs4_cpg_methylation_pct": 5.0, "obs_l2_v3_stage7_chs4_cpg_methylation_se": 0.6,
        "obs_l2_v3_stage7_dm1_cpg_methylation_pct": 6.0, "obs_l2_v3_stage7_dm1_cpg_methylation_se": 0.7,
        "obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment": 12.4, "obs_l2_v3_stage0_chs4_ctcf_chip_fold_enrichment_se": 1.0,
        "obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment": 11.8, "obs_l2_v3_stage0_dm1_ctcf_chip_fold_enrichment_se": 1.0,
        "obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment": 12.0, "obs_l2_v3_stage7_chs4_ctcf_chip_fold_enrichment_se": 1.0,
        "obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment": 11.5, "obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment_se": 1.0,
        "obs_l2_v3_ap_c5b9_pct": 1.8, "obs_l2_v3_ap_c5b9_se": 0.2,
        "obs_l2_v3_classical_c4b2a_lysis_pct": 4.0, "obs_l2_v3_classical_c4b2a_lysis_se": 0.5,
        "obs_l2_v3_lectin_c4b2a_lysis_pct": 4.2, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.5,
        "obs_l2_v3_overall_20pct_nhs_lysis_pct": 2.0, "obs_l2_v3_overall_20pct_nhs_lysis_se": 0.2,
        "obs_l2_v3_intact_serum_adcc_lysis_pct": 1.5, "obs_l2_v3_intact_serum_adcc_lysis_se": 0.2,
        "obs_l2_v3_igg_depleted_serum_adcc_lysis_pct": 0.3, "obs_l2_v3_igg_depleted_serum_adcc_lysis_se": 0.1,
        "obs_l2_v3_igg_increment_adcc_lysis_pct": 1.2, "obs_l2_v3_igg_increment_adcc_lysis_se": 0.2,
        "obs_l2_v3_adcc_serum_donor_pair_id": "SERUM_PAIR_D1_01",
        "obs_l2_v3_adcc_intact_specimen_id": "ADCC_INTACT_SPEC_01",
        "obs_l2_v3_adcc_igg_depleted_specimen_id": "ADCC_IGG_DEP_SPEC_01",
        "obs_l2_v3_plvap_iec_network_t50_h": 24.0, "obs_l2_v3_plvap_iec_network_t50_se": 1.5,
        "obs_l2_v3_plvap_iec_viability_pct": 92.0, "obs_l2_v3_plvap_iec_viability_se": 1.0,
        "obs_l2_v3_exact_stack_effero_abs_pct": 88.0, "obs_l2_v3_exact_stack_effero_abs_se": 1.0,
        "obs_l2_v3_exact_stack_effero_pct_of_wt": 82.0, "obs_l2_v3_exact_stack_effero_pct_of_wt_se": 1.5,
        "obs_l2_v3_pericyte_coverage_pct": 82.0, "obs_l2_v3_pericyte_coverage_se": 1.5,
        "obs_l2_v3_aav6_full_length_monomer_pct": 86.0,
        "obs_l2_v3_aav6_truncated_species_pct": 5.0,
        "obs_l2_v3_aav6_dimeric_or_intermediate_pct": 9.0,
        "obs_l2_v3_small_rnaseq_mir375_retention_pct": 96.5, "obs_l2_v3_small_rnaseq_mir375_retention_se": 1.0,
        "obs_l2_v3_intron_splicing_efficiency_pct": 94.0, "obs_l2_v3_intron_splicing_efficiency_se": 1.2,
        "obs_l2_v3_cd276_kd_pct": 84.0, "obs_l2_v3_cd276_kd_se": 2.0,
        "obs_l2_v3_pvr_kd_pct": 82.0, "obs_l2_v3_pvr_kd_se": 2.0,
        "obs_l2_v3_joint_stratum_nk_lysis_pct": 5.2, "obs_l2_v3_joint_stratum_nk_lysis_se": 0.8,
        "rnai_attribution_scramble_control_id": "L2_V3_CTRL_NT_SCRAMBLE",
        "rnai_attribution_scramble_control_sha256": "d1c8d1bbe01e89dac2256d48b51d2dac03b7cba41d076a0e9bbed6cc249c43ff",
        "rnai_attribution_scramble_background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "rnai_attribution_scramble_sample_id": "S11_SCRAMBLE_CTRL_SPEC_01",
        "rnai_attribution_intronless_control_id": "L2_V3_B2MKO_INTRONLESS_CD200_CTRL",
        "rnai_attribution_intronless_control_sha256": "1" * 64,
        "rnai_attribution_intronless_background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "rnai_attribution_intronless_sample_id": "S11_INTRONLESS_CTRL_SPEC_01",
    }
    l2_scope_cd200 = {
        "product_variant_id": "ISLET1_A4_L2_V3_CO_PRIMARY_CD200",
        "construct_sha256": co_cd200_sha,
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
    }
    l2_t6b_pass = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "L2_P1", product_variant_id="ISLET1_A4_L2_V3_CO_PRIMARY_CD200", construct_sha256=co_cd200_sha, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
        _mk_frac("flow_through_and_excluded", "L2_P2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), product_variant_id="ISLET1_A4_L2_V3_CO_PRIMARY_CD200", construct_sha256=co_cd200_sha, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
    ]

    def assert_only_unverified_controls(nm, out):
        missing_digest = "HOLD_S11_UNVERIFIED_INTRONLESS_CONTROL_DIGEST"
        assert out["s11_evaluation"]["reasons"] == [
            missing_digest,
            "HOLD_S11_UNVERIFIED_ATTRIBUTION_CONTROL_MEASUREMENTS",
        ], (nm, out)
        assert out["s11_evaluation"]["gate_pass"] is False, (nm, out)
        assert out["s11_evaluation"]["status"] == missing_digest, (nm, out)
        assert out["overall_status"] == missing_digest, (nm, out)

    # A specimen label and a submitted digest cannot establish the external controls.
    res_l2_clean = _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
        s11_record=s11_passing_panel,
    )
    for nm, out in res_l2_clean.items():
        assert out["execution_supported"] is False, (nm, out)
        assert_only_unverified_controls(nm, out)

    impossible_se = "HOLD_S11_IMPOSSIBLE_PAIRED_ADCC_UNCERTAINTY"
    for difference_se, is_possible in (
        (0.01, False), (0.09999999999999999, False), (0.10, True), (0.30, True),
        (0.30000000000000004, False), (0.35, False),
    ):
        for nm, out in _eval_4paths(
            l2_t6b_pass, scope_override=l2_scope_cd200, target_prov_override=l2_scope_cd200,
            s11_record={**s11_passing_panel, "obs_l2_v3_igg_increment_adcc_lysis_se": difference_se},
        ).items():
            if is_possible:
                assert_only_unverified_controls(nm, out)
            else:
                assert out["s11_evaluation"]["status"] == impossible_se, (nm, out)
                assert out["overall_status"] == impossible_se, (nm, out)

    # COUNTEREXAMPLE 7 (Sol R6): CD200 comparator with qualifying Tier 6B and OMITTED S11 MUST HOLD on ALL 4 paths!
    res_l2_omitted_s11 = _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
    )
    for nm, out in res_l2_omitted_s11.items():
        assert out["s11_evaluation"]["gate_pass"] is False, (nm, out)
        assert out["s11_evaluation"]["status"] == "HOLD_S11_INCOMPLETE_EVIDENCE", (nm, out)
        assert out["overall_status"] == "HOLD_S11_INCOMPLETE_EVIDENCE", (nm, out)
    assert res_l2_omitted_s11["ledger"]["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE"

    # COUNTEREXAMPLE 8 (Sol R6): S11 panel with SEs removed, invalid replicate count, -100% AP deposition, nested override, or MCB_SCOPE MUST NOT clear on any path!
    panel_no_se = {k: v for k, v in s11_passing_panel.items() if not k.endswith("_se")}
    panel_neg_ap = {**s11_passing_panel, "obs_l2_v3_ap_c5b9_pct": -100.0}
    panel_bad_rep = {**s11_passing_panel, "obs_stage0_independent_lot_replicates_n": 2}
    panel_nested_ovr = {**s11_passing_panel, "lot_id": "NESTED_FOREIGN_LOT_99"}
    for bad_s11_lbl, bad_s11_rec, exp_prefix in (
        ("missing_all_ses", panel_no_se, "HOLD_S11_"),
        ("negative_100pct_ap_deposition", panel_neg_ap, "HOLD_S11_"),
        ("invalid_replicates_n_2", panel_bad_rep, "HOLD_S11_"),
        ("nested_lot_id_override", panel_nested_ovr, "FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH"),
    ):
        res_bad_s11 = _eval_4paths(
            l2_t6b_pass,
            scope_override=l2_scope_cd200,
            target_prov_override=l2_scope_cd200,
            s11_record=bad_s11_rec,
        )
        for nm, out in res_bad_s11.items():
            assert str(out["overall_status"]).startswith(exp_prefix), (bad_s11_lbl, nm, out)
            assert out["overall_status"] != "RESEARCH_COMPARATOR_STAGE0_PANEL_CLEARED_PENDING_STATE4_PROMOTION", (bad_s11_lbl, nm, out)

    # MCB-scoped S11 panel MUST be rejected as lot evidence across all 4 paths:
    res_mcb_s11 = _eval_4paths(
        l2_t6b_pass,
        scope_override={**l2_scope_cd200, "scope_type": "MCB_SCOPE"},
        target_prov_override=l2_scope_cd200,
        s11_record={**s11_passing_panel, "scope_type": "MCB_SCOPE"},
        s11_ledger_scope="MCB_SCOPE",
    )
    for nm, out in res_mcb_s11.items():
        assert out["overall_status"] != "RESEARCH_COMPARATOR_STAGE0_PANEL_CLEARED_PENDING_STATE4_PROMOTION", (nm, out)
        assert "HOLD_S11_" in str(out["overall_status"]) or "HOLD_TIER6B_" in str(out["overall_status"]), (nm, out)
    # Specifically test ledger with valid LOT_SCOPE Tier 6B + MCB_SCOPE S11 panel:
    res_mcb_s11_ledger_only = _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
        s11_record=s11_passing_panel,
        s11_ledger_scope="MCB_SCOPE",
    )["ledger"]
    assert res_mcb_s11_ledger_only["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE"
    assert "HOLD_S11_PANEL_REQUIRES_LOT_SCOPE" in res_mcb_s11_ledger_only["holds"]

    # COUNTEREXAMPLE 9 (Sol R6): 18% lectin failure followed by a later clean S11 batch on the SAME lot MUST stay failed on ALL 4 paths!
    panel_18pct_lectin = {**s11_passing_panel, "obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.5}
    lot_k_cd200 = wc_mod.LotScopeKey("ISLET1_A4_L2_V3_CO_PRIMARY_CD200", co_cd200_sha, SEED_ID, MCB_ID, LOT_ID)
    clean_s11_atom_later = wc_mod.AtomicObservation(
        endpoint_id="stage0_locus2_modeb_v3_panel",
        scope_type="LOT_SCOPE",
        scope_key=lot_k_cd200,
        kind="STRUCTURED_PANEL",
        value_or_k=s11_passing_panel,
        campaign_id=CAMP_ID,
        sample_unit_id="S11_CLEAN_RETRY_02",
    )
    res_18_lectin_latch = _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
        s11_record=panel_18pct_lectin,
        extra_ledger_batches=[{"batch_id": "B_LATER_CLEAN_S11", "atoms": [clean_s11_atom_later]}],
        extra_replay_records=[{
            "utc_timestamp": "2026-11-19T08:00:00Z",
            "seq_index": 99,
            "scope": {**base_scope_31, **l2_scope_cd200},
            "iip0031_locus2_record": s11_passing_panel,
        }],
    )
    for nm, out in res_18_lectin_latch.items():
        assert out["overall_status"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", (nm, out)
        assert out["hard_stop_precedence_winner"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", (nm, out)
        assert out["s11_evaluation"]["status"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", (nm, out)
        if "switch_decisions" in out:
            assert out["switch_decisions"][wc_mod.RULE_S11_GATE_ID_IIP0031] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", (nm, out["switch_decisions"])
            assert out["switch_decisions"][f"JSON_{wc_mod.RULE_S11_GATE_ID_IIP0031}"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", (nm, out["switch_decisions"])

    # Adverse Tier 6B k=1 on non-executable L2_V3_* evaluates BEFORE research hold across all 4 paths:
    res_l2_pos = _eval_4paths(
        [
            _mk_frac("enriched_cd31_neg_cd144_neg", "L2_POS1", k=1, product_variant_id="ISLET1_A4_L2_V3_CO_PRIMARY_CD200", construct_sha256=co_cd200_sha, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
        ],
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
    )
    for nm, out in res_l2_pos.items():
        assert out["overall_status"] == "FAIL_HARD_STOP", (nm, out)
        assert out["hard_stop_precedence_winner"] == "FAIL_TIER6B_OUTGROWTH_DETECTED", (nm, out)

    # Rule S11 hard architectural blocks evaluate BEFORE research hold:
    assert wc_mod.evaluate_stage0_product_contract(
        {**base_scope_31, "product_variant_id": "ISLET1_B_SCT_C1_CIITA", "l2_v3_donor_id": "L2_V3_CO_PRIMARY_CD200", **base_meas_31},
        root=ROOT,
    )["overall_status"] == "FAIL_S11_TRACK2_BSCT_INCOMPATIBLE_AT_LOCUS_2"
    assert wc_mod.evaluate_stage0_product_contract(
        {**base_scope_31, "product_variant_id": "ISLET1_B2_5LOCUS_BASE", "l2_v3_donor_id": "L2_V3_CO_PRIMARY_CD200", **base_meas_31},
        root=ROOT,
    )["overall_status"] == "FAIL_S11_MIR_CD276_PVR_BARRED_ON_DESIGN_B_HLAC_INTACT"
    assert wc_mod.evaluate_stage0_product_contract(
        {**base_scope_31, "product_variant_id": "ISLET1_A4_L2_V3_CO_PRIMARY_CD200", "background_id": "DESIGN_Y_Y03", **base_meas_31},
        root=ROOT,
    )["overall_status"] == "FAIL_S11_SYN5XISRE_CD200_BARRED_IN_DESIGN_Y_Y03"
    assert wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate(
        {"product_variant_id": "ISLET1_Y_RESEARCH", "l2_v3_donor_id": "L2_V3_CO_PRIMARY_CD200"},
        root=ROOT,
    )["status"] == "FAIL_S11_SYN5XISRE_CD200_BARRED_IN_DESIGN_Y_Y03"

    # COUNTEREXAMPLE 10 (Sol R6) + Quantitative Rule S11 stops & straddlers across ALL 4 paths:
    for s11_override, exp_st in (
        # Classical & lectin interval straddlers (11.0%, SE 0.5, n=6 -> L95=9.9925% <= 10.0% < U95=12.0075% -> HOLD):
        ({"obs_l2_v3_classical_c4b2a_lysis_pct": 11.0, "obs_l2_v3_classical_c4b2a_lysis_se": 0.5}, "HOLD_S11_INCONCLUSIVE_CLASSICAL_C4B2A_STRADDLES_10PCT"),
        ({"obs_l2_v3_lectin_c4b2a_lysis_pct": 11.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.5}, "HOLD_S11_INCONCLUSIVE_LECTIN_C4B2A_STRADDLES_10PCT"),
        # Classical & lectin true stops (12.0%, SE 0.5, n=6 -> L95=10.9925% > 10.0% -> FAIL):
        ({"obs_l2_v3_classical_c4b2a_lysis_pct": 12.0, "obs_l2_v3_classical_c4b2a_lysis_se": 0.5}, "FAIL_S11_POST_CUFF_CLASSICAL_C4B2A_LYSIS_EXCEEDS_10PCT"),
        ({"obs_l2_v3_lectin_c4b2a_lysis_pct": 12.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.5}, "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT"),
        ({"obs_l2_v3_ap_c5b9_pct": 5.5, "obs_l2_v3_ap_c5b9_se": 0.2}, "FAIL_S11_POST_CUFF_AP_C5B9_EXCEEDS_5PCT"),
        ({"obs_l2_v3_intact_serum_adcc_lysis_pct": 6.2, "obs_l2_v3_igg_depleted_serum_adcc_lysis_pct": 1.0, "obs_l2_v3_igg_increment_adcc_lysis_pct": 5.2, "obs_l2_v3_igg_increment_adcc_lysis_se": 0.2}, "FAIL_S11_CONJUNCTIVE_GY023_ADCC_EXCEEDS_5PCT"),
        ({"obs_l2_v3_plvap_iec_network_t50_h": 42.0, "obs_l2_v3_plvap_iec_network_t50_se": 1.0}, "FAIL_S11_PLVAP_IEC_NETWORK_MORPHOGENESIS_EXCEEDS_36H"),
        ({"obs_l2_v3_plvap_iec_viability_pct": 75.0, "obs_l2_v3_plvap_iec_viability_se": 1.0}, "FAIL_S11_PLVAP_IEC_VIABILITY_BELOW_85PCT"),
        ({"obs_l2_v3_exact_stack_effero_abs_pct": 75.0, "obs_l2_v3_exact_stack_effero_abs_se": 1.0}, "FAIL_S11_EXACT_STACK_EFFEROCYTOSIS_BELOW_80ABS_OR_70WT"),
        ({"obs_l2_v3_pericyte_coverage_pct": 65.0, "obs_l2_v3_pericyte_coverage_se": 1.0}, "FAIL_S11_GATE0G_PERICYTE_COVERAGE_BELOW_70PCT"),
        ({"obs_l2_v3_aav6_full_length_monomer_pct": 75.0, "obs_l2_v3_aav6_truncated_species_pct": 12.0, "obs_l2_v3_aav6_dimeric_or_intermediate_pct": 13.0}, "FAIL_S11_AAV6_HOMOGENEITY_SPEC_BREACHED"),
        ({"obs_l2_v3_aav6_full_length_monomer_pct": 85.0, "obs_l2_v3_aav6_truncated_species_pct": 5.0, "obs_l2_v3_aav6_dimeric_or_intermediate_pct": 5.0}, "HOLD_S11_AAV6_SPECIES_DENOMINATOR_UNRECONCILED"),
        ({"obs_ciita_ko_allele_informative_reads_n": 25}, "HOLD_S11_MISSING_QUANTITATIVE_CIITA_BIALLELIC_SPLICE_OR_HLAII_EVIDENCE"),
        ({"obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k": 20, "obs_ciita_ko_ex4_ex6_aberrant_skip_pct": 4.0, "obs_ciita_ko_ex4_ex6_aberrant_skip_se": 0.3}, "FAIL_S11_CIITA_ABERRANT_EX4_EX6_SKIP_EXCEEDS_1PCT"),
        ({"obs_edited_72h_ifng_hla_dr_dp_dq_pos_pct": 4.0, "obs_edited_72h_ifng_hla_dr_dp_dq_pos_se": 0.3}, "FAIL_S11_EDITED_CELL_HLA_II_ESCAPE_DETECTED"),
        ({"obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_pct": 20.0, "obs_wt_lineage_matched_72h_ifng_hla_dr_dp_dq_pos_se": 1.5}, "HOLD_S11_UNINDUCIBLE_WT_HLAII_POSITIVE_CONTROL"),
        ({"obs_cis_dexi_tost90_ratio_lo": 1.30, "obs_cis_dexi_tost90_ratio_hi": 1.50}, "FAIL_S11_16P13_CIS_SENTINEL_TOST_BREACHED"),
        ({"obs_cis_dexi_tost90_ratio_lo": 0.75, "obs_cis_dexi_tost90_ratio_hi": 1.05}, "HOLD_S11_INCONCLUSIVE_16P13_CIS_SENTINEL_TOST"),
        ({"obs_l2_v3_unit2_stim_24h_5race_umi_n_eff": 2000}, "HOLD_S11_MISSING_LOCUS2_DUAL_PROMOTER_5PRIME_RACE_3689"),
        ({"obs_l2_v3_unit2_tss_pos1_in_donor": 620}, "FAIL_S11_SYN5XISRE_UPSTREAM_ATG_INCLUDED_OR_INVALID_TSS"),
        ({"obs_l2_v3_stage7_dm1_cpg_methylation_pct": 35.0, "obs_l2_v3_stage7_dm1_cpg_methylation_se": 1.5}, "FAIL_S11_LOCUS2_CPG_HYPERMETHYLATION_EXCEEDS_20PCT"),
        ({"obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment": 2.0, "obs_l2_v3_stage7_dm1_ctcf_chip_fold_enrichment_se": 0.3}, "HOLD_S11_INSUFFICIENT_CTCF_CHIP_OCCUPANCY_BELOW_5X"),
        ({"rnai_specific_benefit_claimed": True, "l2_v3_same_background_intronless_control_linked_to_cargo": False}, "HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS"),
    ):
        res_s11_4p = _eval_4paths(
            l2_t6b_pass,
            scope_override=l2_scope_cd200,
            target_prov_override=l2_scope_cd200,
            s11_record={**s11_passing_panel, **s11_override},
        )
        for nm, out in res_s11_4p.items():
            assert out["s11_evaluation"]["status"] == exp_st, (s11_override, nm, out)
            assert out["overall_status"] == exp_st, (s11_override, nm, out)

    # --- ROUND-7 DIFFERENTIAL CONFORMANCE CASES (ALL 4 ENTRY PATHS) ---
    # R7-1: Unconditional untreated CFU quantitative bound (1/30M & 1/100k claiming 30M/100k MUST HOLD on all 4 paths):
    r7_1_untr_low = [
        _mk_frac(
            "enriched_cd31_neg_cd144_neg",
            "R7_UCFU_1",
            k=0,
            N=30_000_000,
            w_i=0.9999,
            extra={
                "plated_treated_cells_n": 30_000_000,
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_R7_1A", "untreated_cfu_colonies_k": 1, "untreated_plated_cells_n": 30_000_000},
            },
        ),
        _mk_frac(
            "flow_through_and_excluded",
            "R7_UCFU_2",
            k=0,
            N=100_000,
            w_i=0.0001,
            hpsc=(10, 10),
            diff=(10, 10),
            extra={
                "plated_treated_cells_n": 100_000,
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_R7_1B", "untreated_cfu_colonies_k": 1, "untreated_plated_cells_n": 100_000},
            },
        ),
    ]
    for nm, out in _eval_4paths(r7_1_untr_low).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert "UNTREATED_CFU_DOES_NOT_SUPPORT_CLAIMED_OPPORTUNITIES" in str(out["tier6b_reason"]), (nm, out)
        assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)

    # R7-2: Forbid caller-supplied w=(0.99999, 0.00001) overriding 0.5/0.5 shares AND forbid weight_override/w_i_lo/w_i_hi shortcuts:
    r7_2_w_spoof = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "R7_W_1", k=0, N=30_000_000, w_i=0.99999, extra={"presort_total_input_cells": 600_000_000, "presort_fraction_cells": 300_000_000}),
        _mk_frac("flow_through_and_excluded", "R7_W_2", k=0, N=100_000, w_i=0.00001, hpsc=(10, 10), diff=(10, 10), extra={"presort_total_input_cells": 600_000_000, "presort_fraction_cells": 300_000_000}),
    ]
    for nm, out in _eval_4paths(r7_2_w_spoof).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert "FRACTION_WEIGHT_MISMATCH_WITH_VERIFIED_CLONOGENIC_RISK_UNIT_SHARE" in str(out["tier6b_reason"]), (nm, out)
        assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
    for ovr_k, ovr_extra in (
        ("weight_override", {"weight_override": 0.99999}),
        ("w_i_bounds_shortcut", {"w_i_lo": 0.99998, "w_i_hi": 0.99999}),
        ("w_i_upper_shortcut", {"w_i_upper": 0.99999}),
    ):
        r7_2_ovr = [
            _mk_frac("enriched_cd31_neg_cd144_neg", f"R7_OVR1_{ovr_k}", k=0, N=30_000_000, w_i=0.9999, extra=ovr_extra),
            _mk_frac("flow_through_and_excluded", f"R7_OVR2_{ovr_k}", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ]
        for nm, out in _eval_4paths(r7_2_ovr).items():
            assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (ovr_k, nm, out)
            assert "CALLER_WEIGHT_OVERRIDE_OR_UNVERIFIED_BOUNDS_PROHIBITED" in str(out["tier6b_reason"]), (ovr_k, nm, out)
    # Clonogenic-risk-unit share w_i = E_clono_i / sum(E_clono_j) consistency:
    # Enriched: 100M physical * 0.04 CFU = 4M clono (w=4/5=0.8); FT: 2.5M physical * 0.40 CFU = 1M clono (w=1/5=0.2)
    r7_2_clono_mismatch = [
        _mk_frac(
            "enriched_cd31_neg_cd144_neg",
            "R7_CLM1",
            k=0,
            N=3_000_000,
            w_i=100.0 / 102.5,  # raw physical share (wrong risk unit!)
            extra={
                "presort_total_input_cells": 102_500_000,
                "presort_fraction_cells": 100_000_000,
                "plated_treated_cells_n": 75_000_000,
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_CLM1", "untreated_cfu_colonies_k": 20, "untreated_plated_cells_n": 500},
            },
        ),
        _mk_frac(
            "flow_through_and_excluded",
            "R7_CLM2",
            k=0,
            N=500_000,
            w_i=2.5 / 102.5,
            hpsc=(10, 10),
            diff=(10, 10),
            extra={
                "presort_total_input_cells": 102_500_000,
                "presort_fraction_cells": 2_500_000,
                "plated_treated_cells_n": 1_250_000,
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_CLM2", "untreated_cfu_colonies_k": 200, "untreated_plated_cells_n": 500},
            },
        ),
    ]
    for nm, out in _eval_4paths(r7_2_clono_mismatch).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert "FRACTION_WEIGHT_MISMATCH_WITH_VERIFIED_CLONOGENIC_RISK_UNIT_SHARE" in str(out["tier6b_reason"]), (nm, out)

    # R7-3: Calibration batch multiplicity across two distinct 8/8 batches in one fraction -> d_L = 0.4622158621 < 0.50 -> HOLD on all 4 paths:
    for cal_mode, cal_a_extra, cal_b_extra in (
        ("distinct_cal_batch_ids", {"cal_batch_id": "CAL_BATCH_A"}, {"cal_batch_id": "CAL_BATCH_B"}),
        ("omitted_cal_batch_ids_distinct_samples", {}, {}),
    ):
        r7_3_cal_mult = [
            _mk_frac("enriched_cd31_neg_cd144_neg", f"R7_CM1_{cal_mode}", k=0, N=15_000_000, w_i=0.9999, hpsc=(20, 20), diff=(8, 8), extra=cal_a_extra),
            _mk_frac("enriched_cd31_neg_cd144_neg", f"R7_CM2_{cal_mode}", k=0, N=15_000_000, w_i=0.9999, hpsc=(20, 20), diff=(8, 8), extra=cal_b_extra),
            _mk_frac("flow_through_and_excluded", f"R7_CM_FT_{cal_mode}", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ]
        for nm, out in _eval_4paths(r7_3_cal_mult).items():
            ev = out["tier6b_sorted_lot_evaluation"]
            assert ev["disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (cal_mode, nm, ev)
            assert ev["detection_floor_cleared"] is False, (cal_mode, nm, ev)
            assert abs(ev["d_L_by_fraction"]["enriched_cd31_neg_cd144_neg"] - 0.4622158621) < 1e-9, (cal_mode, nm, ev)
            assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (cal_mode, nm, out)

    # R7-4: Conflicting outer vs nested controls['sample_id'] MUST HOLD unqualified (T9), while nested aliquot_id or shared source_sample_id deduplicates to 15M on ALL 4 paths (never clears READY!):
    r7_4_dup_nested_conflict = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "R7_DUP_15M", k=0, N=15_000_000, w_i=0.9999, extra={"controls": {"sample_id": "WRAP_NESTED_1"}}),
        _mk_frac("enriched_cd31_neg_cd144_neg", "R7_DUP_15M", k=0, N=15_000_000, w_i=0.9999, extra={"controls": {"sample_id": "WRAP_NESTED_2"}}),
        _mk_frac("flow_through_and_excluded", "R7_DUP_FT1", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    res_dup_conflict = _eval_4paths(r7_4_dup_nested_conflict)
    for nm, out in res_dup_conflict.items():
        ev = out["tier6b_sorted_lot_evaluation"]
        assert ev["disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, ev)
        assert "CONFLICTING_OUTER_AND_NESTED_IDENTITY" in str(ev["reason_code"]), (nm, ev)
    r7_4_dup_nested_alias = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "R7_DUP_15M", k=0, N=15_000_000, w_i=0.9999, extra={"controls": {"aliquot_id": "WRAP_NESTED_1"}}),
        _mk_frac("enriched_cd31_neg_cd144_neg", "R7_DUP_15M", k=0, N=15_000_000, w_i=0.9999, extra={"controls": {"aliquot_id": "WRAP_NESTED_2"}}),
        _mk_frac("flow_through_and_excluded", "R7_DUP_FT1", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    r7_4_dup_source = [
        _mk_frac(
            "enriched_cd31_neg_cd144_neg",
            "R7_WRAP_OUTER_1",
            k=0,
            N=15_000_000,
            w_i=0.9999,
            extra={
                "source_sample_id": "PHYS_ALIQUOT_15M_A",
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_PHYS_15M_A", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            },
        ),
        _mk_frac(
            "enriched_cd31_neg_cd144_neg",
            "R7_WRAP_OUTER_2",
            k=0,
            N=15_000_000,
            w_i=0.9999,
            extra={
                "source_sample_id": "PHYS_ALIQUOT_15M_A",
                "untreated_cfu_control": {"untreated_sample_id": "UNTR_PHYS_15M_A", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            },
        ),
        _mk_frac("flow_through_and_excluded", "R7_DUP_FT2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for dup_lbl, dup_recs in (("nested_aliquot_ids", r7_4_dup_nested_alias), ("shared_source_sample_id", r7_4_dup_source)):
        res_dup = _eval_4paths(dup_recs)
        for nm, out in res_dup.items():
            ev = out["tier6b_sorted_lot_evaluation"]
            assert ev["N_by_fraction"]["enriched_cd31_neg_cd144_neg"] == 15_000_000, (dup_lbl, nm, ev)
            assert ev["disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (dup_lbl, nm, ev)
            assert out["overall_status"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (dup_lbl, nm, out)
        assert res_dup["ledger"]["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE", (dup_lbl, res_dup["ledger"])

    # R7-5: Strict Donor-ID + SHA-256 Registry Binding (CD200 + correct SHA + spoofed CD274P donor label + missing Unit-2 TSS -> HOLD_S11_DONOR_MISMATCH on all 4 paths):
    panel_no_tss_spoof_cd274p = {k: v for k, v in s11_passing_panel.items() if k != "obs_l2_v3_unit2_tss_pos1_in_donor"}
    l2_t6b_bad_sha = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "L2_BS1", product_variant_id="ISLET1_A4_L2_V3_CO_PRIMARY_CD200", construct_sha256="f" * 64, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
        _mk_frac("flow_through_and_excluded", "L2_BS2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), product_variant_id="ISLET1_A4_L2_V3_CO_PRIMARY_CD200", construct_sha256="f" * 64, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
    ]
    for spoof_lbl, t6_use, sc_spoof, rec_spoof in (
        ("outer_spoofed_cd274p_donor_id", l2_t6b_pass, {**l2_scope_cd200, "l2_v3_donor_id": "L2_V3_CD274P_ZERO_TATATA"}, panel_no_tss_spoof_cd274p),
        ("nested_spoofed_cd274p_donor_id", l2_t6b_pass, l2_scope_cd200, {**panel_no_tss_spoof_cd274p, "donor_id": "L2_V3_CD274P_ZERO_TATATA"}),
        ("mismatched_donor_sha256", l2_t6b_bad_sha, {**l2_scope_cd200, "construct_sha256": "f" * 64}, s11_passing_panel),
    ):
        res_spoof = _eval_4paths(
            t6_use,
            scope_override=sc_spoof,
            target_prov_override=sc_spoof,
            s11_record=rec_spoof,
        )
        for nm, out in res_spoof.items():
            assert out["s11_evaluation"]["status"] == "HOLD_S11_DONOR_MISMATCH", (spoof_lbl, nm, out)
            assert out["overall_status"] == "HOLD_S11_DONOR_MISMATCH", (spoof_lbl, nm, out)

    # R7-6: Reject False verification flags, missing 4th RACE, missing Stage-0 methylation, n_reads=100 (< 300), missing WGS/readthrough/QuantiBRITE/active-RNAi/matched B2M-/- controls, or missing direct lot/sample ID:
    r7_6_s11_holds = [
        ("phasing_verified_false", {**s11_passing_panel, "phasing_verified": False}, "HOLD_S11_CIITA_BIALLELIC_PHASING_UNVERIFIED"),
        ("cis_verified_false", {**s11_passing_panel, "cis_verified": False}, "HOLD_S11_16P13_CIS_SENTINELS_UNVERIFIED"),
        ("wt_inducibility_verified_false", {**s11_passing_panel, "wt_inducibility_verified": False}, "HOLD_S11_WT_HLAII_INDUCIBILITY_UNVERIFIED"),
        ("rnai_attribution_verified_false", {**s11_passing_panel, "rnai_attribution_verified": False}, "HOLD_S11_RNAI_ATTRIBUTION_UNVERIFIED"),
        ("missing_4th_race_eef1a1_24h", {k: v for k, v in s11_passing_panel.items() if k != "obs_l2_v3_unit1_eef1a1_stim_24h_5race_umi_n_eff"}, "HOLD_S11_MISSING_LOCUS2_DUAL_PROMOTER_5PRIME_RACE_3689"),
        ("missing_stage0_methylation", {k: v for k, v in s11_passing_panel.items() if not k.startswith("obs_l2_v3_stage0_eef1a1_cpg_methylation_")}, "HOLD_S11_MISSING_LONGITUDINAL_STAGE0_AND_STAGE7_METHYLATION_OR_CTCF_PANEL"),
        ("ciita_reads_100_below_300", {**s11_passing_panel, "obs_ciita_ko_allele_informative_reads_n": 100}, "HOLD_S11_MISSING_QUANTITATIVE_CIITA_BIALLELIC_SPLICE_OR_HLAII_EVIDENCE"),
        ("missing_chr16_ha_wgs", {k: v for k, v in s11_passing_panel.items() if k != "obs_l2_v3_chr16_ha_wgs_pass"}, "HOLD_S11_MISSING_CHR16_HA_WGS"),
        ("missing_tandem_polya_readthrough", {k: v for k, v in s11_passing_panel.items() if not k.startswith("obs_l2_v3_tandem_polya_readthrough_")}, "HOLD_S11_MISSING_TANDEM_POLYA_READTHROUGH_EVIDENCE"),
        ("missing_quantibrite_cfh", {k: v for k, v in s11_passing_panel.items() if not k.startswith("obs_l2_v3_cfh_stp_gpi_quantibrite_")}, "HOLD_S11_MISSING_COPY_SPECIFIC_QUANTIBRITE_EXPRESSION"),
        ("missing_small_rnaseq_mir375", {k: v for k, v in s11_passing_panel.items() if not k.startswith("obs_l2_v3_small_rnaseq_mir375_retention_")}, "HOLD_S11_MISSING_ACTIVE_MIR_PROCESSING_KD_OR_JOINT_STRATUM_NK_PANEL"),
        ("design_b_hlac_intact_wrong_rnai_ctrl", {**s11_passing_panel, "rnai_attribution_intronless_control_id": "L2_V3_DESIGN_B_CD200"}, "HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS"),
    ]
    for r7_lbl, r7_rec, r7_exp in r7_6_s11_holds:
        for nm, out in _eval_4paths(
            l2_t6b_pass,
            scope_override=l2_scope_cd200,
            target_prov_override=l2_scope_cd200,
            s11_record=r7_rec,
        ).items():
            assert out["s11_evaluation"]["status"] == r7_exp, (r7_lbl, nm, out)
            assert out["overall_status"] == r7_exp, (r7_lbl, nm, out)
    # Direct S11 without lot_id or sample_id MUST HOLD (`HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING`):
    s11_no_lot_direct = wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate(
        {"product_variant_id": "ISLET1_A4_L2_V3_CO_PRIMARY_CD200", "construct_sha256": co_cd200_sha, "iip0031_locus2_record": {k: v for k, v in s11_passing_panel.items() if k != "sample_id"}},
        root=ROOT,
    )
    assert s11_no_lot_direct["gate_pass"] is False and s11_no_lot_direct["status"] == "HOLD_S11_MISSING_LOT_OR_SAMPLE_BINDING", s11_no_lot_direct
    # Genotype-matched allele evidence leaves exactly the unresolved control reasons.
    s11_kiki_panel = {
        **{
            k: v
            for k, v in s11_passing_panel.items()
            if not k.startswith("obs_ciita_k")
            and k not in ("ciita_ki_allele_informative_transcript_recovery_verified", "ciita_ko_allele_informative_transcript_recovery_verified")
        },
        "ciita_locus2_genotype_copy_class": "KI/KI",
        "ciita_ki_allele1_informative_transcript_recovery_verified": True,
        "ciita_ki_allele2_informative_transcript_recovery_verified": True,
        "obs_ciita_allele1_informative_reads_n": 400,
        "obs_ciita_allele2_informative_reads_n": 400,
        "obs_ciita_allele1_ex4_ex6_aberrant_skip_reads_k": 0,
        "obs_ciita_allele2_ex4_ex6_aberrant_skip_reads_k": 0,
        "obs_ciita_allele1_ex4_ex6_aberrant_skip_pct": 0.0,
        "obs_ciita_allele1_ex4_ex6_aberrant_skip_se": 0.05,
        "obs_ciita_allele2_ex4_ex6_aberrant_skip_pct": 0.0,
        "obs_ciita_allele2_ex4_ex6_aberrant_skip_se": 0.05,
    }
    for nm, out in _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
        s11_record=s11_kiki_panel,
    ).items():
        assert_only_unverified_controls(nm, out)

    # R7-7: Zero VAF with missing duplex depth MUST HOLD on ALL 4 paths (including atomic ledger — no silent 62,958 default!):
    s9_vaf_missing_depth_atom = wc_mod.AtomicObservation(
        endpoint_id="obs_tp53_ppm1d_bcor_vaf",
        scope_type="LOT_SCOPE",
        scope_key=lot_k,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=0.0,
        se_or_n=1e-6,
        sample_n=6,
        campaign_id=CAMP_ID,
        sample_unit_id="S9_VAF_NODEPTH",
        controls={},
    )
    meas_no_depth = {k: v for k, v in base_meas_31.items() if k != "obs_tp53_ppm1d_bcor_duplex_depth"}
    meas_no_depth["obs_tp53_ppm1d_bcor_vaf"] = 0.0
    meas_no_depth["obs_tp53_ppm1d_bcor_vaf_se"] = 1e-6
    p1_nd = wc_mod.evaluate_stage0_product_contract(
        {**base_scope_31, **meas_no_depth, **base_mol_31, "tier6b_sorted_lot_observations": list(fix1_recs)},
        root=ROOT,
    )
    p2_nd = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
        [{"batch_id": "B_ND", "atoms": mcb_atoms_base + dev_atoms_base + [s9_vaf_missing_depth_atom] + [_to_atom(r) for r in fix1_recs]}],
        target_prov,
        target_campaign_id=CAMP_ID,
        root=ROOT,
    )
    p3_nd = ing_mod.ingest_stage0_research_or_wetlab_batch(
        json.dumps({"scope": dict(base_scope_31), "measurements": meas_no_depth, "iip0026_molecular_observations": dict(base_mol_31), "tier6b_sorted_lot_observations": list(fix1_recs)}),
        root=ROOT,
    )
    p4_nd = ing_mod.replay_stage0_chronological_batches(
        [{"utc_timestamp": "2026-11-19T01:00:00Z", "seq_index": 1, "scope": dict(base_scope_31), "measurements": meas_no_depth, "iip0026_molecular_observations": dict(base_mol_31), "tier6b_sorted_lot_observations": list(fix1_recs)}],
        root=ROOT,
    )
    assert p1_nd["overall_status"].startswith("HOLD_"), p1_nd
    assert p3_nd["overall_status"].startswith("HOLD_"), p3_nd
    assert p4_nd["overall_status"].startswith("HOLD_"), p4_nd
    assert p2_nd["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE" and any("HOLD_S9_MISSING_OR_INVALID_DUPLEX_DEPTH" in h for h in p2_nd["holds"]), p2_nd

    # R7-8: Single shared S9 VAF guard: VAF=4.0e-5, SE=3.1e-5, depth=500,000, n=6 -> IDENTICAL HARD_SAFETY_STOP_REJECTED on ALL 4 paths!
    s9_vaf_border_atom = wc_mod.AtomicObservation(
        endpoint_id="obs_tp53_ppm1d_bcor_vaf",
        scope_type="LOT_SCOPE",
        scope_key=lot_k,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=4.0e-5,
        se_or_n=3.1e-5,
        sample_n=6,
        campaign_id=CAMP_ID,
        sample_unit_id="S9_VAF_BORDER_T",
        controls={"obs_tp53_ppm1d_bcor_duplex_depth": 500000},
    )
    res_r7_8_s9 = _eval_4paths(
        fix1_recs,
        meas_override={
            "obs_tp53_ppm1d_bcor_vaf": 4.0e-5,
            "obs_tp53_ppm1d_bcor_vaf_se": 3.1e-5,
            "obs_tp53_ppm1d_bcor_duplex_depth": 500000,
            "obs_tp53_ppm1d_bcor_vaf_sample_n": 6,
        },
        extra_ledger_atoms=[s9_vaf_border_atom],
    )
    for nm, out in res_r7_8_s9.items():
        assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (nm, out)
        assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (nm, out)

    # A valid signed paired difference leaves exactly the unresolved control reasons.
    res_r7_10_signed_adcc = _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
        s11_record={**s11_passing_panel, "obs_l2_v3_igg_depleted_serum_adcc_lysis_pct": 1.8, "obs_l2_v3_igg_increment_adcc_lysis_pct": -0.3, "obs_l2_v3_igg_increment_adcc_lysis_se": 0.2},
    )
    for nm, out in res_r7_10_signed_adcc.items():
        assert_only_unverified_controls(nm, out)
    res_r7_10_neg_intact = _eval_4paths(
        l2_t6b_pass,
        scope_override=l2_scope_cd200,
        target_prov_override=l2_scope_cd200,
        s11_record={**s11_passing_panel, "obs_l2_v3_intact_serum_adcc_lysis_pct": -0.3, "obs_l2_v3_intact_serum_adcc_lysis_se": 0.2},
    )
    for nm, out in res_r7_10_neg_intact.items():
        assert out["s11_evaluation"]["status"] == "HOLD_S11_INVALID_SE_OR_OUT_OF_RANGE_PHYSIOLOGICAL_READOUT", (nm, out)

    # --- ROUND-8 DIFFERENTIAL CONFORMANCE CASES (ALL 4 ENTRY PATHS + PERMUTATION & ADVERSARIAL FUZZ) ---
    # R8-1: Finite untreated CFU (1/1 -> 17k/500M AND even 70/2M) MUST HOLD unconditionally on all 4 paths;
    #       explicit statistically_estimated_denominator_or_share=False opt-out MUST HOLD;
    #       90% reconciled sort losses without explicit disposal attestation MUST HOLD:
    r8_1_finite_1of1 = [
        _mk_frac(
            "enriched_cd31_neg_cd144_neg",
            "R8_1OF1_E",
            k=0,
            N=17_000,
            w_i=None,
            extra={
                "presort_total_input_cells": 500_017_002,
                "presort_fraction_cells": 17_001,
                "plated_treated_cells_n": 17_000,
                "untreated_cfu_control": {"untreated_sample_id": "R8_U1E", "untreated_cfu_colonies_k": 1, "untreated_plated_cells_n": 1},
            },
        ),
        _mk_frac(
            "flow_through_and_excluded",
            "R8_1OF1_F",
            k=0,
            N=500_000_000,
            w_i=None,
            extra={
                "presort_total_input_cells": 500_017_002,
                "presort_fraction_cells": 500_000_001,
                "plated_treated_cells_n": 500_000_000,
                "untreated_cfu_control": {"untreated_sample_id": "R8_U1F", "untreated_cfu_colonies_k": 1, "untreated_plated_cells_n": 1},
            },
        ),
    ]
    for nm, out in _eval_4paths(r8_1_finite_1of1).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert out["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
    # Even if caller lowers N to 70 & 2,083,333 on 1/1 CFU, worst-case E-box bound still >> 3e-7 -> HOLD on all 4 paths:
    r8_1_finite_1of1_lowN = [
        {**r8_1_finite_1of1[0], "sample_id": "R8_1OF1_E_LOW", "N": 70},
        {**r8_1_finite_1of1[1], "sample_id": "R8_1OF1_F_LOW", "N": 2_083_333},
    ]
    for nm, out in _eval_4paths(r8_1_finite_1of1_lowN).items():
        assert out["tier6b_disposition"] in ("HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", "HOLD_TIER6B_UNQUALIFIED_FRACTION"), (nm, out)
        assert out["overall_status"] != "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED", (nm, out)
    # Opting out via statistically_estimated_denominator_or_share=False MUST HOLD on all 4 paths:
    r8_1_opt_out = [
        {**fix1_recs[0], "sample_id": "R8_OPT1", "statistically_estimated_denominator_or_share": False},
        {**fix1_recs[1], "sample_id": "R8_OPT2", "statistically_estimated_denominator_or_share": False},
    ]
    for nm, out in _eval_4paths(r8_1_opt_out).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert out["tier6b_reason"] == "EXACT_FINITE_CFU_OPT_OUT_PROHIBITED", (nm, out)
    # 90% reconciled sort losses without explicit excluded_or_lost_material_discarded_not_in_graft=True MUST HOLD:
    r8_1_loss90 = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "R8_L90_E", k=0, N=30_000_000, w_i=0.9999, extra={"presort_total_input_cells": 30_000_000_000, "presort_reconciled_loss_cells": 27_000_000_000, "presort_fraction_cells": 2_999_700_000}),
        _mk_frac("flow_through_and_excluded", "R8_L90_F", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), extra={"presort_total_input_cells": 30_000_000_000, "presort_reconciled_loss_cells": 27_000_000_000, "presort_fraction_cells": 300_000}),
    ]
    for nm, out in _eval_4paths(r8_1_loss90).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert out["tier6b_reason"] == "UNVERIFIED_DISPOSAL_OF_RECONCILED_SORT_LOSSES", (nm, out)

    # R8-2a: Order-independent transitive alias dedup across ALL 6 permutations of [A, B, Link(A,B)] -> 15M HOLD on ALL 4 paths!
    import itertools as _it_r8
    rec_alias_a = _mk_frac("enriched_cd31_neg_cd144_neg", "R8_ALI_A", k=0, N=15_000_000, w_i=0.9999, extra={"untreated_cfu_control": {"untreated_sample_id": "R8_ALI_U", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000}})
    rec_alias_b = _mk_frac("enriched_cd31_neg_cd144_neg", "R8_ALI_B", k=0, N=15_000_000, w_i=0.9999, extra={"untreated_cfu_control": {"untreated_sample_id": "R8_ALI_U", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000}})
    rec_alias_link = _mk_frac("enriched_cd31_neg_cd144_neg", "R8_ALI_A", k=0, N=15_000_000, w_i=0.9999, extra={"source_sample_id": "R8_ALI_B", "untreated_cfu_control": {"untreated_sample_id": "R8_ALI_U", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000}})
    rec_alias_ft = _mk_frac("flow_through_and_excluded", "R8_ALI_FT", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10))
    for perm_idx, perm_t in enumerate(_it_r8.permutations((rec_alias_a, rec_alias_b, rec_alias_link))):
        res_perm = _eval_4paths([*perm_t, rec_alias_ft])
        for nm, out in res_perm.items():
            ev = out["tier6b_sorted_lot_evaluation"]
            assert ev["N_by_fraction"]["enriched_cd31_neg_cd144_neg"] == 15_000_000, (perm_idx, nm, ev)
            assert ev["disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (perm_idx, nm, ev)
            assert out["overall_status"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (perm_idx, nm, out)
        assert res_perm["ledger"]["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE", (perm_idx, res_perm["ledger"])

    # R8-2b: Dual-location Tier 6B (clean top-level + nested k=1 in mol_obs) MUST FAIL_HARD_STOP on Direct, Ingest, AND Replay identically!
    k1_hidden = _mk_frac("enriched_cd31_neg_cd144_neg", "R8_HIDDEN_K1", k=1, N=100_000, w_i=0.9999)
    dual_doc = {
        "scope": dict(base_scope_31),
        "measurements": dict(base_meas_31),
        "iip0026_molecular_observations": {**base_mol_31, "tier6b_sorted_lot_observations": [k1_hidden]},
        "tier6b_sorted_lot_observations": list(fix1_recs),
    }
    dual_direct = wc_mod.evaluate_stage0_product_contract(
        {**base_scope_31, **base_meas_31, **base_mol_31, "iip0026_molecular_observations": {"tier6b_sorted_lot_observations": [k1_hidden]}, "tier6b_sorted_lot_observations": list(fix1_recs)},
        root=ROOT,
    )
    dual_ingest = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(dual_doc), root=ROOT)
    dual_replay = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 1, **dual_doc}], root=ROOT)
    for nm_d, out_d in (("direct", dual_direct), ("ingest", dual_ingest), ("replay", dual_replay)):
        assert out_d["overall_status"] == "FAIL_HARD_STOP", (nm_d, out_d)
        assert out_d["hard_stop_precedence_winner"] == "FAIL_TIER6B_OUTGROWTH_DETECTED", (nm_d, out_d)

    # R8-3: S11 Round-8 catches across ALL 4 paths (1/300 CP U95=1.57% > 1%, omitted attestations, KI/KI+KI/KO mismatch, reused RNAi control ID, arbitrary background_id, 103% implied ADCC):
    r8_3_s11_cases = [
        ("ciita_1_of_300_cp_u95_gt_1pct", {**s11_passing_panel, "obs_ciita_ko_allele_informative_reads_n": 300, "obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k": 1, "obs_ciita_ko_ex4_ex6_aberrant_skip_pct": 0.333, "obs_ciita_ko_ex4_ex6_aberrant_skip_se": 0.01}, "HOLD_S11_INCONCLUSIVE_CIITA_EX4_EX6_SKIP_STRADDLES_1PCT"),
        ("omitted_phasing_attestations", {k: v for k, v in s11_passing_panel.items() if k not in ("phasing_verified", "ciita_ki_allele_informative_transcript_recovery_verified", "ciita_ko_allele_informative_transcript_recovery_verified")}, "HOLD_S11_CIITA_BIALLELIC_PHASING_UNVERIFIED"),
        ("omitted_cis_attestations", {k: v for k, v in s11_passing_panel.items() if k not in ("cis_verified", "ciita_cis_sentinels_socs1_dexi_clec16a_tost_verified")}, "HOLD_S11_16P13_CIS_SENTINELS_UNVERIFIED"),
        ("omitted_wt_attestations", {k: v for k, v in s11_passing_panel.items() if k not in ("wt_inducibility_verified", "lineage_matched_wt_stage7_hla_dr_dp_dq_inducibility_verified")}, "HOLD_S11_WT_HLAII_INDUCIBILITY_UNVERIFIED"),
        ("omitted_rnai_attestation", {k: v for k, v in s11_passing_panel.items() if k != "rnai_attribution_verified"}, "HOLD_S11_RNAI_ATTRIBUTION_UNVERIFIED"),
        ("kiki_label_with_kiko_fields", {**s11_passing_panel, "ciita_locus2_genotype_copy_class": "KI/KI"}, "HOLD_S11_GENOTYPE_CLASS_ALLELE_FIELD_MISMATCH"),
        ("reused_fictitious_rnai_ctrl_id", {**s11_passing_panel, "rnai_attribution_scramble_control_id": "FAKE_CTRL_SAME", "rnai_attribution_intronless_control_id": "FAKE_CTRL_SAME"}, "HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS"),
        ("paired_adcc_103pct_depleted_arm", {**s11_passing_panel, "obs_l2_v3_intact_serum_adcc_lysis_pct": 3.0, "obs_l2_v3_intact_serum_adcc_lysis_se": 0.2, "obs_l2_v3_igg_increment_adcc_lysis_pct": -100.0, "obs_l2_v3_igg_increment_adcc_lysis_se": 0.2}, "HOLD_S11_INVALID_SE_OR_OUT_OF_RANGE_PHYSIOLOGICAL_READOUT"),
    ]
    for r8_lbl, r8_rec, r8_exp in r8_3_s11_cases:
        for nm, out in _eval_4paths(
            l2_t6b_pass,
            scope_override=l2_scope_cd200,
            target_prov_override=l2_scope_cd200,
            s11_record=r8_rec,
        ).items():
            assert out["s11_evaluation"]["status"] == r8_exp, (r8_lbl, nm, out)
            assert out["overall_status"] == r8_exp, (r8_lbl, nm, out)
    # Arbitrary background_id ("FAKE_BG") overriding registry background MUST HOLD with HOLD_S11_BACKGROUND_MISMATCH on all 4 paths:
    fake_bg_scope = {**l2_scope_cd200, "background_id": "FAKE_BG"}
    fake_bg_rec = {**s11_passing_panel, "rnai_attribution_scramble_background_id": "FAKE_BG", "rnai_attribution_intronless_background_id": "FAKE_BG"}
    fake_bg_t6b = [{**r, "background_id": "FAKE_BG"} for r in l2_t6b_pass]
    for nm, out in _eval_4paths(fake_bg_t6b, scope_override=fake_bg_scope, target_prov_override=fake_bg_scope, s11_record=fake_bg_rec).items():
        assert out["s11_evaluation"]["status"] == "HOLD_S11_BACKGROUND_MISMATCH", (nm, out)
        assert out["overall_status"] == "HOLD_S11_BACKGROUND_MISMATCH", (nm, out)

    # R8-4a: TP53 VAF=5e-5, SE=1e-6 with duplex depth OMITTED -> IDENTICAL HARD_SAFETY_STOP_REJECTED on ALL 4 paths!
    s9_5e5_no_depth_atom = wc_mod.AtomicObservation(
        endpoint_id="obs_tp53_ppm1d_bcor_vaf",
        scope_type="LOT_SCOPE",
        scope_key=lot_k,
        kind="CONTINUOUS_STUDENT_T",
        value_or_k=5.0e-5,
        se_or_n=1e-6,
        sample_n=6,
        campaign_id=CAMP_ID,
        sample_unit_id="S9_5E5_NODEPTH",
        controls={},
    )
    meas_5e5_nd = {k: v for k, v in base_meas_31.items() if k != "obs_tp53_ppm1d_bcor_duplex_depth"}
    meas_5e5_nd["obs_tp53_ppm1d_bcor_vaf"] = 5.0e-5
    meas_5e5_nd["obs_tp53_ppm1d_bcor_vaf_se"] = 1e-6
    p1_5e5 = wc_mod.evaluate_stage0_product_contract({**base_scope_31, **meas_5e5_nd, **base_mol_31, "tier6b_sorted_lot_observations": list(fix1_recs)}, root=ROOT)
    p2_5e5 = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"batch_id": "B_5E5", "atoms": mcb_atoms_base + dev_atoms_base + [s9_5e5_no_depth_atom] + [_to_atom(r) for r in fix1_recs]}], target_prov, target_campaign_id=CAMP_ID, root=ROOT)
    p3_5e5 = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps({"scope": dict(base_scope_31), "measurements": meas_5e5_nd, "iip0026_molecular_observations": dict(base_mol_31), "tier6b_sorted_lot_observations": list(fix1_recs)}), root=ROOT)
    p4_5e5 = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 1, "scope": dict(base_scope_31), "measurements": meas_5e5_nd, "iip0026_molecular_observations": dict(base_mol_31), "tier6b_sorted_lot_observations": list(fix1_recs)}], root=ROOT)
    assert p1_5e5["overall_status"] == p3_5e5["overall_status"] == p4_5e5["overall_status"] == "HARD_SAFETY_STOP_REJECTED"
    assert p2_5e5["readiness_tier"] == "HARD_SAFETY_STOP_REJECTED", p2_5e5

    # R8-4b: Earlier S9 hard stop + later S11 18% lectin failure -> overall retains S9 precedence WHILE switch_decisions preserves latched S11 FAIL!
    replay_s9_then_s11 = ing_mod.replay_stage0_chronological_batches(
        [
            {
                "utc_timestamp": "2026-11-20T01:00:00Z",
                "seq_index": 1,
                "scope": {**base_scope_31, **l2_scope_cd200},
                "measurements": {**base_meas_31, "obs_tp53_ppm1d_bcor_vaf": 2.5e-4, "obs_tp53_ppm1d_bcor_vaf_se": 1e-6},
                "iip0026_molecular_observations": dict(base_mol_31),
                "tier6b_sorted_lot_observations": list(l2_t6b_pass),
                "iip0031_locus2_record": s11_passing_panel,
            },
            {
                "utc_timestamp": "2026-11-20T02:00:00Z",
                "seq_index": 2,
                "scope": {**base_scope_31, **l2_scope_cd200},
                "measurements": dict(base_meas_31),
                "iip0031_locus2_record": panel_18pct_lectin,
            },
            {
                "utc_timestamp": "2026-11-20T03:00:00Z",
                "seq_index": 3,
                "scope": {**base_scope_31, **l2_scope_cd200},
                "measurements": dict(base_meas_31),
                "iip0031_locus2_record": s11_passing_panel,
            },
        ],
        root=ROOT,
    )
    assert replay_s9_then_s11["overall_status"] == "HARD_SAFETY_STOP_REJECTED"
    assert replay_s9_then_s11["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP"
    assert replay_s9_then_s11["s11_evaluation"]["status"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT"
    assert replay_s9_then_s11["switch_decisions"][wc_mod.RULE_S11_GATE_ID_IIP0031] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT"
    assert replay_s9_then_s11["switch_decisions"][f"JSON_{wc_mod.RULE_S11_GATE_ID_IIP0031}"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT"

    # R8-5: Removing jaspar_pwms from Layer-2 reference MUST raise ValueError (never silently report 0 hits / empty CTCF!):
    raw_ref_no_pwm = json.loads((ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text(encoding="utf-8"))
    raw_ref_no_pwm.pop("jaspar_pwms", None)
    try:
        l2_mod.build_locus2_ciita_modeb_comparators_v3_manifest(ref=raw_ref_no_pwm)
        raise AssertionError("Expected ValueError when jaspar_pwms is missing")
    except ValueError as exc_pwm:
        assert "JASPAR" in str(exc_pwm) or "PWM" in str(exc_pwm), str(exc_pwm)

    # --- ROUND-9 DIFFERENTIAL 4-PATH REPLAY REGRESSIONS (Tier 6B nested aliases & hybrid labels, S9 completeness & n=100 parity, S11 strict attribution & longitudinal/barrier/ADCC provenance) ---
    # R9-1a: Conflicting outer vs nested controls['sample_unit_id'] MUST HOLD unqualified (T9), while shared nested aliquot_id deduplicates identically to 15M HOLD across ALL 4 paths!
    r9_conflict_enr_1 = _mk_frac(
        "enriched_cd31_neg_cd144_neg",
        "R9_OUTER_A",
        k=0,
        N=15_000_000,
        w_i=0.9999,
        hpsc=(20, 20),
        diff=(20, 20),
        extra={
            "sample_unit_id": "R9_OUTER_A",
            "untreated_cfu_control": {"untreated_sample_id": "R9_SHARED_UNTR_CFU", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            "controls": {"sample_id": "R9_OUTER_A", "sample_unit_id": "R9_SHARED_ALIQUOT"},
        },
    )
    r9_alias_enr_1 = _mk_frac(
        "enriched_cd31_neg_cd144_neg",
        "R9_OUTER_A",
        k=0,
        N=15_000_000,
        w_i=0.9999,
        hpsc=(20, 20),
        diff=(20, 20),
        extra={
            "sample_unit_id": "R9_OUTER_A",
            "untreated_cfu_control": {"untreated_sample_id": "R9_SHARED_UNTR_CFU", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            "controls": {"sample_id": "R9_OUTER_A", "aliquot_id": "R9_SHARED_ALIQUOT"},
        },
    )
    r9_alias_enr_2 = _mk_frac(
        "enriched_cd31_neg_cd144_neg",
        "R9_SHARED_ALIQUOT",
        k=0,
        N=15_000_000,
        w_i=0.9999,
        hpsc=(20, 20),
        diff=(20, 20),
        extra={
            "untreated_cfu_control": {"untreated_sample_id": "R9_SHARED_UNTR_CFU", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
        },
    )
    r9_ft = _mk_frac("flow_through_and_excluded", "R9_FT", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10))
    for nm_c, out_c in _eval_4paths([r9_conflict_enr_1, r9_alias_enr_2, r9_ft]).items():
        ev_c = out_c["tier6b_sorted_lot_evaluation"]
        assert ev_c["disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm_c, ev_c)
        assert "CONFLICTING_OUTER_AND_NESTED_IDENTITY" in str(ev_c["reason_code"]), (nm_c, ev_c)
    for perm_r9 in ([r9_alias_enr_1, r9_alias_enr_2, r9_ft], [r9_alias_enr_2, r9_alias_enr_1, r9_ft]):
        res_r9_al = _eval_4paths(perm_r9)
        for nm, out in res_r9_al.items():
            ev_al = out["tier6b_sorted_lot_evaluation"]
            assert ev_al["disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, ev_al)
            assert ev_al["N_by_fraction"]["enriched_cd31_neg_cd144_neg"] == 15_000_000, (nm, ev_al)
            assert abs(ev_al["p_lot_U95"] - 4.2857727075e-7) < 1e-11, (nm, ev_al)

    # R9-1b: Hybrid fraction label ("enriched_flow_through" + "secondary_enriched") MUST HOLD across ALL 4 paths!
    r9_hybrid_roster = ("enriched_flow_through", "secondary_enriched")
    r9_hybrid_recs = [
        _mk_frac("enriched_flow_through", "R9_H1", k=0, N=30_000_000, w_i=0.9999, roster=r9_hybrid_roster),
        _mk_frac("secondary_enriched", "R9_H2", k=0, N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), roster=r9_hybrid_roster),
    ]
    for nm, out in _eval_4paths(r9_hybrid_recs).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert out["tier6b_reason"] == "AMBIGUOUS_HYBRID_ENRICHED_FLOW_THROUGH_FRACTION_ID", (nm, out)

    # R9-2a: Omitting core S9 lot observations MUST HOLD across ALL 4 paths (including atomic ledger — NEVER READY_FOR_ENCAPSULATION)!
    p2_no_s9 = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
        [{"batch_id": "B_NOS9", "atoms": [_to_atom(r) for r in fix1_recs]}],
        target_prov,
        target_campaign_id=CAMP_ID,
        root=ROOT,
    )
    assert p2_no_s9["readiness_tier"] == "HOLD_FAILED_STAGE0_GATE", p2_no_s9
    assert p2_no_s9["overall_status"] == "HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS", p2_no_s9
    meas_no_tp53 = {k: v for k, v in base_meas_31.items() if k != "obs_tp53_ppm1d_bcor_vaf"}
    p1_no_s9 = wc_mod.evaluate_stage0_product_contract(
        {**base_scope_31, **meas_no_tp53, **base_mol_31, "tier6b_sorted_lot_observations": list(fix1_recs)},
        root=ROOT,
    )
    assert p1_no_s9["overall_status"] == "HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS", p1_no_s9

    # R9-2b: Endpoint replicate count n=100 vs n=6 on synthetic VAF (vaf=2e-5, se=4.5e-5, depth=70000) MUST match across ALL 4 paths!
    # At n=100: t95(99)=1.66039 -> 2e-5 + 1.66039*4.5e-5 = 9.47e-5 < 1e-4 (PASS); at n=6: t95(5)=2.01505 -> 1.107e-4 >= 1e-4 (HARD STOP)!
    res_r9_n100 = _eval_4paths(
        fix1_recs,
        meas_override={"obs_tp53_ppm1d_bcor_vaf": 2.0e-5, "obs_tp53_ppm1d_bcor_vaf_se": 4.5e-5, "obs_tp53_ppm1d_bcor_vaf_n": 100},
    )
    assert res_r9_n100["direct"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED", res_r9_n100["direct"]
    assert res_r9_n100["ingest"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED", res_r9_n100["ingest"]
    assert res_r9_n100["replay"]["overall_status"] == "STAGE0_IIP0026_MOLECULAR_GUARDRAILS_CLEARED", res_r9_n100["replay"]
    assert res_r9_n100["ledger"]["readiness_tier"] == "PROVISIONAL_PARTIAL_LOT_SCREENING_ONLY", res_r9_n100["ledger"]
    res_r9_n6 = _eval_4paths(
        fix1_recs,
        meas_override={"obs_tp53_ppm1d_bcor_vaf": 2.0e-5, "obs_tp53_ppm1d_bcor_vaf_se": 4.5e-5, "obs_tp53_ppm1d_bcor_vaf_n": 6},
    )
    for nm, out in res_r9_n6.items():
        ov_n6 = out["readiness_tier"] if nm == "ledger" else out["overall_status"]
        assert ov_n6 == "HARD_SAFETY_STOP_REJECTED", (nm, out)

    # R9-3a/b/c/d: Unregistered intronless control ID, KI/KO with KI/KI recovery attestations, conflicting outer vs nested sample_id,
    # untimed CTCF ChIP, missing/identical GY-023 IgG-depleted specimens, and L2_V3_BARRIER_CHS4_250 linked 80-bp comparator!
    r9_s11_cases = [
        ("unregistered_intronless_ctrl_id", {**s11_passing_panel, "rnai_attribution_intronless_control_id": "UNREGISTERED_FAKE_CTRL_99"}, "HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS"),
        ("wrong_cargo_intronless_ctrl_id", {**s11_passing_panel, "rnai_attribution_intronless_control_id": "L2_V3_B2MKO_INTRONLESS_MANF_CTRL"}, "HOLD_S11_RNAI_ATTRIBUTION_REQUIRES_SAME_BACKGROUND_CARGO_LINKED_INTRONLESS_AND_SCRAMBLE_CONTROLS"),
        ("kiko_with_only_kiki_recovery_attestations", {**{k: v for k, v in s11_passing_panel.items() if k not in ("ciita_ki_allele_informative_transcript_recovery_verified", "ciita_ko_allele_informative_transcript_recovery_verified")}, "ciita_ki_allele1_informative_transcript_recovery_verified": True, "ciita_ki_allele2_informative_transcript_recovery_verified": True}, "HOLD_S11_GENOTYPE_CLASS_ALLELE_FIELD_MISMATCH"),
        ("missing_stage0_ctcf_chip_untimed_only", {**{k: v for k, v in s11_passing_panel.items() if not k.startswith("obs_l2_v3_stage0_chs4_ctcf_")}, "obs_l2_v3_chs4_ctcf_chip_fold_enrichment": 12.0, "obs_l2_v3_chs4_ctcf_chip_fold_enrichment_se": 1.0}, "HOLD_S11_MISSING_LONGITUDINAL_STAGE0_AND_STAGE7_METHYLATION_OR_CTCF_PANEL"),
        ("missing_gy023_igg_depleted_arm", {k: v for k, v in s11_passing_panel.items() if not k.startswith("obs_l2_v3_igg_depleted_")}, "HOLD_S11_UNVERIFIED_PAIRED_GY023_ADCC_SPECIMENS"),
        ("identical_gy023_intact_and_depleted_specimen_ids", {**s11_passing_panel, "obs_l2_v3_adcc_igg_depleted_specimen_id": "ADCC_INTACT_SPEC_01"}, "HOLD_S11_UNVERIFIED_PAIRED_GY023_ADCC_SPECIMENS"),
    ]
    for r9_lbl, r9_rec, r9_exp in r9_s11_cases:
        for nm, out in _eval_4paths(
            l2_t6b_pass,
            scope_override=l2_scope_cd200,
            target_prov_override=l2_scope_cd200,
            s11_record=r9_rec,
        ).items():
            assert out["s11_evaluation"]["status"] == r9_exp, (r9_lbl, nm, out)
            assert out["overall_status"] == r9_exp, (r9_lbl, nm, out)

    # R9-3c: Conflicting outer sample_id / sample_unit_id vs nested sample_id MUST FAIL identically across ALL 4 paths!
    for nm, out in _eval_4paths(
        l2_t6b_pass,
        scope_override={**l2_scope_cd200, "sample_id": "OUTER_S11_CONFLICT_ID"},
        target_prov_override=l2_scope_cd200,
        s11_record=s11_passing_panel,
    ).items():
        assert out["s11_evaluation"]["status"] == "FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH", (nm, out)
        assert out["overall_status"] == "FAIL_S11_NESTED_PROVENANCE_SCOPE_MISMATCH", (nm, out)

    # Linking the barrier measurement clears only its own missing-evidence reason.
    bar250_def = by_id["ISLET1_A4_L2_V3_BARRIER_CHS4_250"]
    bar250_sha = bar250_def["staged_locus2_modeb_donor_binding"]["locus2_donor_sha256"]
    bar250_scope = {"product_variant_id": "ISLET1_A4_L2_V3_BARRIER_CHS4_250", "construct_sha256": bar250_sha, "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}
    bar250_t6b = [
        _mk_frac("enriched_cd31_neg_cd144_neg", "B250_P1", product_variant_id="ISLET1_A4_L2_V3_BARRIER_CHS4_250", construct_sha256=bar250_sha, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
        _mk_frac("flow_through_and_excluded", "B250_P2", N=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), product_variant_id="ISLET1_A4_L2_V3_BARRIER_CHS4_250", construct_sha256=bar250_sha, extra={"background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}),
    ]
    for nm, out in _eval_4paths(bar250_t6b, scope_override=bar250_scope, target_prov_override=bar250_scope, s11_record=s11_passing_panel).items():
        assert out["s11_evaluation"]["status"] == "HOLD_S11_MISSING_LINKED_BARRIER_80BP_COMPARATOR_MEASUREMENT", (nm, out)
    bar250_panel_linked = {
        **s11_passing_panel,
        "obs_l2_v3_barrier_matched_80bp_comparator_id": "L2_V3_CO_PRIMARY_CD200",
        "obs_l2_v3_barrier_matched_80bp_comparator_sha256": co_cd200_sha,
        "obs_l2_v3_barrier_matched_80bp_comparator_sample_id": "S11_LINKED_80BP_COMP_SPEC_01",
        "obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_pct": 4.5,
        "obs_l2_v3_barrier_matched_80bp_stage7_eef1a1_cpg_methylation_se": 0.5,
        "obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment": 11.2,
        "obs_l2_v3_barrier_matched_80bp_stage7_chs4_ctcf_chip_fold_enrichment_se": 0.9,
    }
    for nm, out in _eval_4paths(bar250_t6b, scope_override=bar250_scope, target_prov_override=bar250_scope, s11_record=bar250_panel_linked).items():
        assert out["execution_supported"] is False, (nm, out)
        assert_only_unverified_controls(nm, out)

    # Composite cHS4-250 disclosure classification:
    chs4_ev = wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate(
        {"product_variant_id": "L2_V3_BARRIER_CHS4_250"},
        root=ROOT,
    )
    assert chs4_ev["barrier_intervention_classification"] == "COMBINED_LENGTH_AND_SEQUENCE_CHANGE_250BP_PLUS_U78775_2_C79T_C80A"
    assert chs4_ev["gate_pass"] is False and chs4_ev["status"] == "HOLD_S11_INCOMPLETE_EVIDENCE"

    # Verify S8b preserved for TYK2 and S11 present in stage0_decision_tree_v1.json:
    dt_doc = json.loads((ROOT / "data/protocols/stage0_decision_tree_v1.json").read_text(encoding="utf-8"))
    smap = dt_doc["switch_id_to_contract_decision_key_map"]
    assert smap["S8b"] == "S8b_TYK2_Inducible_Promoter_Epistasis_Gate_0E"
    assert smap["S11"] == "S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate"

    # --- LAYER 2 MANIFEST & <= 1 mm >= 20 bp DISJOINT IR SCANNER ---
    stem_a = "ACGTTGCATGCAACGTTGCA"
    rc_map = str.maketrans("ACGT", "TGCA")
    stem_b_exact = stem_a.translate(rc_map)[::-1]
    stem_b_1mm = stem_b_exact[:9] + ("A" if stem_b_exact[9] != "A" else "C") + stem_b_exact[10:]
    synth_seq = "GGGGGG" + stem_a + ("A" * 25) + stem_b_1mm + "CCCCCC"
    synth_hits = l2_mod.scan_disjoint_inverted_repeats_with_mismatches(synth_seq, min_k=20, max_mismatches=1)
    assert len(synth_hits) >= 1 and synth_hits[0]["mismatches"] == 1 and synth_hits[0]["length_bp"] >= 20

    l2_atlas = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    l2_v3_man = l2_atlas["stage0_locus2_ciita_modeb_comparators_v3_manifest_iip_0031"]
    assert l2_v3_man["proposal_id"] == "IIP-0031" and l2_v3_man["decision_id"] == "DEC-110"
    assert l2_v3_man["provenance_computation_mode"] == "DYNAMICALLY_RECOMPUTED_END_TO_END_FROM_BUNDLED_FASTA_AND_REFERENCE_JSON"
    assert l2_v3_man["mode_b_comparator_fasta_sha256"] == "25ffcb10e2bad77a4b34725c7eff0532990ebf7b80864cf7217d33f5ccbe5b6d"
    assert l2_v3_man["repeat_summary"]["max_disjoint_dr_across_8_bp"] == 29
    assert l2_v3_man["repeat_summary"]["overlapping_dr_syn5xisre_bp"] == 52
    assert l2_v3_man["repeat_summary"]["overlapping_dr_cd274p_maz4_bp"] == 35
    assert l2_v3_man["repeat_summary"]["max_dr_outside_canonical_tandems_bp"] == 16
    assert l2_v3_man["repeat_summary"]["disjoint_ir_ge_16bp_total_across_8"] == 0
    assert l2_v3_man["repeat_summary"]["disjoint_ir_ge_20bp_le_1mm_total_across_8"] == 0
    assert l2_v3_man["intronic_mir_and_rnai_attribution_scope"]["candidate_seed_pairs_screened"] == 296
    assert l2_v3_man["intronic_mir_and_rnai_attribution_scope"]["replayed_8mer_hits_total"] == 0
    assert l2_v3_man["intronic_mir_and_rnai_attribution_scope"]["replayed_7mer_m8_hits_total"] == 0
    assert len(l2_v3_man["barrier_chs4_250_disclosure"]["shared_80bp_diffs"]) == 2
    assert [d["u78775_2_pos_1based"] for d in l2_v3_man["barrier_chs4_250_disclosure"]["shared_80bp_diffs"]] == [79, 80]
    # Round-7 Layer-2 provenance completeness: exact 8 CIITA isoforms, spatial 6-PWM JASPAR audit, 10/10 backup guides vs 15 OI proxies + 22 UTRs, zero bare assert:
    assert l2_v3_man["ciita_all_8_refseq_isoforms_audit"]["isoforms_audited_count"] == 8
    jaspar_sp = l2_v3_man["descriptive_6pwm_jaspar_spatial_audit"]
    assert jaspar_sp["wholly_non_promoter_6pwm_hits"] == 0
    assert jaspar_sp["promoter_boundary_irf1_nt538_558_6pwm_hits_count"] == 3
    assert jaspar_sp["promoter_boundary_irf1_nt538_558_1pwm_hits_count"] == 7
    ctcf_hits = jaspar_sp["co_primary_cd200_ctcf_insulator_hits"]
    assert abs(ctcf_hits["cHS4"]["p_adj_1pwm"] - 0.045041) < 1e-4
    assert abs(ctcf_hits["cHS4"]["p_adj_6pwm"] - 0.270243) < 1e-4
    assert abs(ctcf_hits["DM1"]["p_adj_1pwm"] - 0.011812) < 1e-4
    assert abs(ctcf_hits["DM1"]["p_adj_6pwm"] - 0.070875) < 1e-4
    bk_ledger = l2_v3_man["intronic_mir_and_rnai_attribution_scope"]["backup_guide_screen_ledger_15oi_22utr"]
    assert len(bk_ledger) == 10
    assert l2_v3_man["intronic_mir_and_rnai_attribution_scope"]["backup_guides_screened_against_all_15_openislet_proxies_and_22_utrs"] is True
    assert all(g["openislet_15_proxies_8mer_or_7mer_m8_hits"] == 0 and g["human_22_refseq_3utrs_8mer_or_7mer_m8_hits"] == 0 for g in bk_ledger)
    # Missing any 1 of the 8 CIITA isoforms MUST raise ValueError (never silently accept 7!):
    raw_ref_copy = json.loads((ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text(encoding="utf-8"))
    raw_ref_copy["iip0031_locus2_modeb_v3_reference_bundle"]["ciita_refseq_isoforms_8"].pop("NM_001379334.1")
    try:
        l2_mod.build_locus2_ciita_modeb_comparators_v3_manifest(ref=raw_ref_copy)
        raise AssertionError("Expected ValueError when only 7 CIITA isoforms are present")
    except ValueError as exc_ciita:
        assert "8" in str(exc_ciita), str(exc_ciita)

    # Zero bare `assert` statements in Layer-2 and Layer-3 production scripts (survives python3 -O / test_master_molecular_validation.py):
    import ast as _ast_mod
    for prod_rel in ("scripts/run_layer2_perturbation_and_regulatory_plane.py", "scripts/run_layer3_structural_binder_plane.py"):
        prod_tree = _ast_mod.parse((ROOT / prod_rel).read_text(encoding="utf-8"))
        assert not [n for n in _ast_mod.walk(prod_tree) if isinstance(n, _ast_mod.Assert)], prod_rel

    # --- LAYER 3 STRUCTURAL PROVENANCE (PDB 2WII, 5O32, Val62) ---
    l3_atlas = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    l3_prov = l3_atlas["locus2_cfh_stp_gpi_structural_provenance_iip_0031"]
    assert l3_prov["proposal_id"] == "IIP-0031" and l3_prov["decision_id"] == "DEC-110"
    assert l3_prov["coordinate_computation_mode"] == "DYNAMICALLY_RECOMPUTED_FROM_BUNDLED_PDB_SNAPSHOTS_ON_DISK"
    assert l3_prov["pdb_2wii"]["cfh_c3b_heavy_atom_pairs_le_4p5A"] == 475
    assert l3_prov["pdb_2wii"]["post_glu264_deposited_coordinate_audit"]["native_residue_265_geometry_claimed"] is False
    assert l3_prov["pdb_2wii"]["post_glu264_deposited_coordinate_audit"]["exit_vector_from_residue_265_claimed"] is False
    assert l3_prov["pdb_5o32"]["asu_copy_1"]["last_cfh_residue_contacting_factor_i_le_4p5A"] == 197
    assert l3_prov["pdb_5o32"]["asu_copy_2"]["last_cfh_residue_contacting_factor_i_le_4p5A"] == 197
    assert l3_prov["pdb_2wii"]["engineered_v3_allotype"] == "Val62 (canonical UniProt P08603)"
    assert l3_prov["val62_vs_ile62_swissprot_census"]["record_provenance_label"] == "archived_uniprot_2026_04_descriptive_record"
    assert l3_prov["val62_vs_ile62_swissprot_census"]["recomputed_from_bundled_disk_inputs"] is False

    # --- GY-025 (Provisional Stage-0 Guardrail) & DEC-110 ---
    gy_doc = json.loads((ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    og_list = gy_doc["off_grid_molecular_dead_end_guardrails"]["guardrails"]
    assert len(gy_doc["graveyard_entries"]) == 10 and len(og_list) == 15
    gy25 = [g for g in og_list if g["graveyard_card_id"] == "GY-025"]
    assert len(gy25) == 1 and gy25[0]["guardrail_id"].startswith("OFFGRID_GUARD_015_") and gy25[0]["contribution_id"] == "IIP-0031" and gy25[0]["severity_level"] == "BLOCK"
    assert gy25[0]["curation_status"] == "provisional_stage0_computational_guardrail"
    assert gy25[0]["status"].startswith("PROVISIONAL_STAGE0_GUARDRAIL_")
    assert "permanently barred" not in gy25[0]["summary"].lower()
    assert "DEC-110" in (ROOT / "ISLET-1-BLUEPRINT.md").read_text(encoding="utf-8")
    pdl_path = ROOT / "docs/PROJECT-DECISIONS-LOG.md"
    if pdl_path.exists():
        assert "DEC-110" in pdl_path.read_text(encoding="utf-8")

    # --- LAYER 1 2,160/2,160 & OQ-016 INVARIANCE ---
    pareto_path = ROOT / "data/extracted/islet1_combinatorial_pareto_and_ladder_v0.8.0.json"
    assert hashlib.sha256(pareto_path.read_bytes()).hexdigest() == "591ddeb9a706c0bc95fc509765b512aadad8a403287a796e63dd52e7db548575"
    l1_pareto = json.loads(pareto_path.read_text(encoding="utf-8"))
    assert len(l1_pareto["evaluated_stacks"]) == 2160
    rk_mod = _load_mod("rk_iip0031_55", ROOT / "scripts/rank_islet1_stacks.py")
    oq16 = rk_mod._build_odds_composition_and_profile_scope_audit_iip_0003()
    assert oq16["open_question_id"] == "OQ-016"
    assert hashlib.sha256(json.dumps(oq16, sort_keys=True).encode("utf-8")).hexdigest() == "9ea7e8be4ee7e6a368d4fa3c3cc6385c2bbeadb88c4e9f37e4407a4a77c906a7"


S9_TEST_ENDPOINTS = {
    "obs_icasp9_log10_kill": ("obs_icasp9_log10_kill_se", None),
    "obs_tp53_ppm1d_bcor_vaf": ("obs_tp53_ppm1d_bcor_vaf_se", "obs_tp53_ppm1d_bcor_duplex_depth"),
    "obs_castseq_translocation_freq": ("obs_castseq_translocation_se", "obs_castseq_duplex_or_amplicon_depth"),
}


def s9_incomplete_lot_results(product_id, measurements, atom_overrides=None):
    """Evaluate safety data in a lot that lacks the independent evidence needed for clearance."""
    from scripts import ingest_wetlab_results as ingest
    from scripts import wetlab_contract as contract

    scope = {
        "product_variant_id": product_id,
        "construct_sha256": "7" * 64,
        "donor_seedline_id": "SEED-S9",
        "mcb_id": "MCB-S9",
        "lot_id": "LOT-S9",
        "campaign_id": "CAMPAIGN-S9",
        "device_assembly_id": "DEVICE-S9",
        "physical_route_id": "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH",
    }
    lot_key = contract.LotScopeKey(
        product_id, scope["construct_sha256"], scope["donor_seedline_id"], scope["mcb_id"], scope["lot_id"]
    )
    atoms = []
    for endpoint, (se_key, depth_key) in S9_TEST_ENDPOINTS.items():
        fields = {
            "kind": "CONTINUOUS_STUDENT_T",
            "value_or_k": measurements[endpoint],
            "se_or_n": measurements[se_key],
            "sample_n": measurements.get(endpoint + "_n", 6),
            "controls": {depth_key: measurements[depth_key]} if depth_key in measurements else {},
        }
        fields.update((atom_overrides or {}).get(endpoint, {}))
        atoms.append(contract.AtomicObservation(
            endpoint_id=endpoint, scope_type="LOT_SCOPE", scope_key=lot_key,
            campaign_id=scope["campaign_id"], sample_unit_id=endpoint, **fields,
        ))
    molecular = {"bound_c1_c5_pair_sha256": contract.DESIGN_B_C1_C5_PAIR_SHA256} if product_id == "ISLET1_B2_5LOCUS_BASE" else {}
    batch = {"scope": scope, "measurements": measurements, "iip0026_molecular_observations": molecular}
    return {
        "direct": contract.evaluate_stage0_product_contract({**scope, **measurements, **molecular}),
        "ledger": contract.evaluate_iip0024_atomic_stage0_ledger(
            [{"atoms": atoms}], scope, target_campaign_id=scope["campaign_id"]
        ),
        "json": ingest.ingest_stage0_research_or_wetlab_batch(json.dumps(batch)),
        "replay": ingest.replay_stage0_chronological_batches([{"seq_index": 1, **batch}]),
    }


@pytest.fixture
def s9_measurements():
    return {
        "obs_icasp9_log10_kill": 4.8,
        "obs_icasp9_log10_kill_se": 0.12,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70_000,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70_000,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300_000_000,
    }


@pytest.mark.parametrize("product_id", ["ISLET1_A4_ARM1B", "ISLET1_B2_5LOCUS_BASE"])
@pytest.mark.parametrize("depth_key", ["obs_tp53_ppm1d_bcor_duplex_depth", "obs_castseq_duplex_or_amplicon_depth"])
@pytest.mark.parametrize("depth", [None, 60_000, 62_957, 62_958])
def test_s9_input_integrity_across_entry_paths(s9_measurements, product_id, depth_key, depth):
    if depth is None:
        del s9_measurements[depth_key]
    else:
        s9_measurements[depth_key] = depth
    for path, result in s9_incomplete_lot_results(product_id, s9_measurements).items():
        reasons = result["holds"] if path == "ledger" else result["iip0026_reasons"]
        assert result["overall_status"].startswith("HOLD_"), (path, reasons)
        depth_hold = any("DUPLEX" in reason and "62958" in reason for reason in reasons)
        if depth is None:
            assert any(
                reason.startswith(("HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS", "HOLD_S9_MISSING_OR_INVALID_DUPLEX_DEPTH"))
                for reason in reasons
            ), (path, reasons)
        else:
            assert depth_hold == (depth < 62_958), (path, reasons)


@pytest.mark.parametrize("endpoint", S9_TEST_ENDPOINTS)
@pytest.mark.parametrize(("kind", "value", "count"), [
    ("METADATA_CONTROL", True, None),
    ("SCALAR_EXACT", 1, None),
    ("EXACT_INTEGER_COUNT", 0, 70_000),
    ("STRUCTURED_PANEL", {}, None),
])
def test_s9_input_integrity_wrong_kind(s9_measurements, endpoint, kind, value, count):
    # Scalar transports carry a Boolean in the numeric field. Ledger represents
    # the equivalent failure as a legal atom kind for the wrong evidence type.
    s9_measurements[endpoint] = True
    overrides = {endpoint: {"kind": kind, "value_or_k": value, "se_or_n": count, "sample_n": None}}
    for path, result in s9_incomplete_lot_results("ISLET1_A4_ARM1B", s9_measurements, overrides).items():
        assert result["overall_status"].startswith("HOLD_"), path
        if path == "ledger":
            assert any("KIND" in reason and endpoint in reason for reason in result["holds"]), result["holds"]
        else:
            assert "HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS" in result["iip0026_reasons"], path


@pytest.mark.parametrize("endpoint", S9_TEST_ENDPOINTS)
@pytest.mark.parametrize("n", [6, 100])
def test_s9_input_integrity_replicate_bounds(s9_measurements, endpoint, n):
    se_key, _depth_key = S9_TEST_ENDPOINTS[endpoint]
    s9_measurements[endpoint], s9_measurements[se_key] = (
        (4.18, 0.10) if endpoint == "obs_icasp9_log10_kill" else (2e-5, 4.5e-5)
    )
    s9_measurements[endpoint + "_n"] = n
    for path, result in s9_incomplete_lot_results("ISLET1_A4_ARM1B", s9_measurements).items():
        if n == 6:
            assert result["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", path
            assert result["overall_status"] == "HARD_SAFETY_STOP_REJECTED", path
        else:
            assert result["hard_stop_precedence_winner"] is None, path
            assert result["overall_status"].startswith("HOLD_"), path


@pytest.mark.parametrize("endpoint", S9_TEST_ENDPOINTS)
def test_s9_input_integrity_replay_replacement_cannot_borrow_assay_fields(s9_measurements, endpoint):
    from scripts import ingest_wetlab_results as ingest

    se_key, depth_key = S9_TEST_ENDPOINTS[endpoint]
    scope = {"product_variant_id": "ISLET1_A4_ARM1B", "lot_id": "LOT-S9", "campaign_id": "CAMPAIGN-S9"}
    replacement = {endpoint: s9_measurements[endpoint]}
    if depth_key is not None:
        replacement[se_key] = s9_measurements[se_key]
    result = ingest.replay_stage0_chronological_batches([
        {"scope": scope, "seq_index": 1, "measurements": s9_measurements},
        {"scope": scope, "seq_index": 2, "measurements": replacement},
    ])
    assert "HOLD_MISSING_CORE_S9_MOLECULAR_OBSERVATIONS" in result["iip0026_reasons"]














def test_iip0031_binomial_confidence_limits_preserve_probability_tails(monkeypatch):
    from math import comb, isclose
    from scripts import wetlab_contract as contract

    lower = contract.clopper_pearson_lower95_independent
    upper = contract.clopper_pearson_upper95_independent
    for n in (5, 13, 29):
        for k in range(n + 1):
            for alpha in (0.005, 0.05, 0.75):
                p_lo, p_hi = lower(k, n, alpha), upper(k, n, alpha)
                if k:
                    right = sum(comb(n, i) * p_lo ** i * (1-p_lo) ** (n-i) for i in range(k, n+1))
                    assert isclose(right, alpha, rel_tol=3e-12)
                else:
                    assert p_lo == 0
                if k != n:
                    left = sum(comb(n, i) * p_hi ** i * (1-p_hi) ** (n-i) for i in range(k+1))
                    assert isclose(left, alpha, rel_tol=3e-12)
                else:
                    assert p_hi == 1

    # Independently evaluated with SciPy's inverse regularized beta functions.
    assert abs(lower(500_000_000, 1_000_000_000, 0.05/6) - 0.49996214735588584) < 1e-13
    assert abs(upper(500_000_000, 1_000_000_000, 0.05/6) - 0.5000378526441142) < 1e-13
    assert isclose(lower(1, 1_000_000_000, 1e-15), 1e-24, rel_tol=1e-14)
    assert 0 < lower(2, 10, 1e-150) < 1e-74
    for limit in (lower, upper):
        for k, n in ((True, 2), (1, False), (0, 0), (1.0, 4), (-1, 5), (8, 7)):
            with pytest.raises(ValueError):
                limit(k, n)
        for k in (0, 1, 3):
            for alpha in (True, False, '0.05', None, float('nan'), float('inf'), 0, 1):
                with pytest.raises(ValueError):
                    limit(k, 3, alpha)
    assert lower(2, 2**60) == 0 and upper(2, 2**60) == 1

    # Force a fresh inversion to exhaust its numeric solver before either bracket moves.
    def failed_cf(*args):
        raise ArithmeticError('unresolved tail')

    contract._cp_interior_limit.cache_clear()
    monkeypatch.setattr(contract, '_betacf_iip0024', failed_cf)
    try:
        assert lower(3, 11) == 0 and upper(3, 11) == 1
    finally:
        contract._cp_interior_limit.cache_clear()


def test_56_iip0031_four_path_nested_flat_poisson_and_identity_parity() -> None:
    """Four-path conformance regressions for nested Tier 6B containers, flat replay replacements, exact Poisson boundaries, negated fraction labels, and S11 campaign conflicts."""
    import json
    from scripts import ingest_wetlab_results as ing_mod
    from scripts import product_contract as pc_mod
    from scripts import wetlab_contract as wc_mod

    ROOT = Path(__file__).resolve().parent.parent
    prods = {p["product_variant_id"]: p for p in pc_mod.load_product_contracts(ROOT)["products"]}
    pvid_a = "ISLET1_A4_ARM1B"
    pvid_l2 = "ISLET1_A4_L2_V3_CO_PRIMARY_CD200"
    sha_a = "a" * 64
    sha_l2 = prods[pvid_l2]["staged_locus2_modeb_donor_binding"]["locus2_donor_sha256"]
    seed_id, mcb_id, lot_id, camp_id = "SEED_R11", "MCB_R11", "LOT_R11", "CAMP_R11"
    dev_id, route_id = "DEV_R11", "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH"

    def _mk_t6(fid: str, sid: str, *, k: int = 0, n: int = 30_000_000, w_i: float = 0.9999,
               roster: list[str] | None = None, tot_cells: int = 3_000_000_000, frac_cells: int | None = None,
               hpsc: list[int] | None = None, diff: list[int] | None = None) -> dict:
        r_use = roster if roster is not None else ["enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"]
        fc = frac_cells if frac_cells is not None else (2_999_700_000 if fid == r_use[0] else 300_000)
        return {
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "kind": "EXACT_INTEGER_COUNT",
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "product_variant_id": pvid_a,
            "construct_sha256": sha_a,
            "donor_seedline_id": seed_id,
            "mcb_id": mcb_id,
            "lot_id": lot_id,
            "campaign_id": camp_id,
            "fraction_id": fid,
            "declared_fraction_roster": list(r_use),
            "sample_id": sid,
            "k": k,
            "N": n,
            "w_i": w_i,
            "hpsc_calibrant": list(hpsc or [20, 20]),
            "edited_differentiated_calibrant": list(diff or [10, 10]),
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "untreated_cfu_control": {
                "untreated_sample_id": f"{sid}_UNTR",
                "untreated_cfu_colonies_k": 4000,
                "untreated_plated_cells_n": 10000,
            },
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": tot_cells,
            "presort_reconciled_loss_cells": 0,
            "presort_fraction_cells": fc,
            "plated_treated_cells_n": int((n * 8 + 2) // 3),
        }

    clean_t6 = [
        _mk_t6("enriched_cd31_neg_cd144_neg", "ENR_1"),
        _mk_t6("flow_through_and_excluded", "FT_1", n=100_000, w_i=0.0001, hpsc=[10, 10], diff=[10, 10]),
    ]
    base_scope = {
        "product_variant_id": pvid_a, "construct_sha256": sha_a,
        "donor_seedline_id": seed_id, "mcb_id": mcb_id, "lot_id": lot_id, "campaign_id": camp_id,
    }
    base_meas = {
        "obs_stage7_yield_ratio": 0.86, "obs_stage7_yield_se": 0.02,
        "obs_cd47_fold": 4.2, "obs_cd47_se": 0.15,
        "obs_serpinb9_log10_kass": 6.5, "obs_serpinb9_log10_kass_se": 0.05,
        "obs_icasp9_log10_kill": 4.8, "obs_icasp9_log10_kill_se": 0.10,
        "obs_pouch_min_po2_mmhg": 24.0, "obs_pouch_min_po2_se": 1.0,
        "obs_pouch_max_h2o2_uM": 1.8, "obs_pouch_max_h2o2_se": 0.2,
        "obs_pouch_max_ph": 7.30, "obs_pouch_max_ph_se": 0.02,
        "obs_ifng_hla_dr_dq_negative_pct": 99.8, "obs_ifng_hla_dr_dq_negative_se": 0.05,
        "obs_cd47_hour0_molecules": 180000.0, "obs_cd47_hour0_molecules_se": 5000.0,
        "obs_hlae_hour0_molecules": 45000.0, "obs_hlae_hour0_molecules_se": 1500.0,
        "obs_tp53_ppm1d_bcor_vaf": 0.0, "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0,
        "obs_castseq_translocation_freq": 0.0, "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000, "obs_castseq_mutant_junctions_k": 0,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 350_000_000,
        "obs_chop_coincidence_fold": 1.1, "obs_chop_coincidence_se": 0.05,
        "obs_hc10_surface_fold_wt": 1.05, "obs_hc10_surface_fold_se": 0.05,
        "obs_post_edit_segmented_baf_chr19q": 0.50,
        "obs_post_edit_segmented_baf_chr11p15": 0.50,
    }
    base_mol = {
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
    }
    lot_k = wc_mod.LotScopeKey(pvid_a, sha_a, seed_id, mcb_id, lot_id)
    target_prov = {**base_scope, "device_assembly_id": dev_id, "physical_route_id": route_id}
    s9_atoms = [
        wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=4.8, se_or_n=0.10, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="S9_IC9"),
        wc_mod.AtomicObservation(endpoint_id="obs_tp53_ppm1d_bcor_vaf", scope_type="LOT_SCOPE", scope_key=lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="S9_VAF",
                                 controls={"duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0}),
        wc_mod.AtomicObservation(endpoint_id="obs_castseq_translocation_freq", scope_type="LOT_SCOPE", scope_key=lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="S9_CAST",
                                 controls={"duplex_depth": 70000, "obs_castseq_mutant_junctions_k": 0}),
    ]

    # 1. Atomic Ledger MUST inspect nested measurements["tier6b_sorted_lot_observations"] identically to Direct/Ingest/Replay:
    pos_nested = _mk_t6("enriched_cd31_neg_cd144_neg", "ENR_1", k=1)
    conf_tot_nested = _mk_t6("enriched_cd31_neg_cd144_neg", "ENR_1", tot_cells=3_500_000_000)
    for container_key in ("measurements", "iip0026_molecular_observations"):
        l_pos = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
            [{"batch_id": "B1", "atoms": list(s9_atoms), "tier6b_sorted_lot_observations": list(clean_t6),
              container_key: {"tier6b_sorted_lot_observations": [pos_nested]}}],
            target_prov, target_campaign_id=camp_id, root=ROOT,
        )
        assert l_pos["overall_status"] == "FAIL_HARD_STOP", (container_key, l_pos)
        assert l_pos["tier6b_reason"] == "FAIL_TIER6B_OUTGROWTH_DETECTED", (container_key, l_pos)

        l_conf = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
            [{"batch_id": "B1", "atoms": list(s9_atoms), "tier6b_sorted_lot_observations": list(clean_t6),
              container_key: {"tier6b_sorted_lot_observations": [conf_tot_nested]}}],
            target_prov, target_campaign_id=camp_id, root=ROOT,
        )
        assert l_conf["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (container_key, l_conf)
        assert l_conf["overall_status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (container_key, l_conf)

    # 2. Chronological Replay MUST process flat top-level obs_* S9 replacements identically to nested measurements:
    rp_flat_adverse = ing_mod.replay_stage0_chronological_batches([
        {"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": dict(base_scope),
         "measurements": dict(base_meas), "iip0026_molecular_observations": dict(base_mol),
         "tier6b_sorted_lot_observations": list(clean_t6)},
        {"utc_timestamp": "2026-11-20T02:00:00Z", "seq_index": 1, "scope": dict(base_scope),
         "obs_tp53_ppm1d_bcor_vaf": 2e-4, "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
         "obs_tp53_ppm1d_bcor_duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 14},
    ], root=ROOT)
    assert rp_flat_adverse["overall_status"] == "HARD_SAFETY_STOP_REJECTED", rp_flat_adverse
    assert rp_flat_adverse["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", rp_flat_adverse

    # Flat adverse replacement followed by clean structured batch MUST stay latched:
    rp_flat_then_clean = ing_mod.replay_stage0_chronological_batches([
        {"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": dict(base_scope),
         "measurements": dict(base_meas), "iip0026_molecular_observations": dict(base_mol),
         "tier6b_sorted_lot_observations": list(clean_t6)},
        {"utc_timestamp": "2026-11-20T02:00:00Z", "seq_index": 1, "scope": dict(base_scope),
         "obs_tp53_ppm1d_bcor_vaf": 2e-4, "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
         "obs_tp53_ppm1d_bcor_duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 14},
        {"utc_timestamp": "2026-11-20T03:00:00Z", "seq_index": 2, "scope": dict(base_scope),
         "measurements": dict(base_meas)},
    ], root=ROOT)
    assert rp_flat_then_clean["overall_status"] == "HARD_SAFETY_STOP_REJECTED", rp_flat_then_clean

    # 3. Exact Poisson upper bound boundary at k=16, n=243,000 MUST exceed 1e-4 and hard-stop on ALL 4 paths:
    u95_16 = wc_mod._poisson_exact_upper95_rate(16, 243_000)
    assert abs(u95_16 - 1.000048711261e-4) < 1e-10 and u95_16 > 1.0e-4, u95_16
    meas_16 = {
        **base_meas,
        "obs_tp53_ppm1d_bcor_vaf": 16.0 / 243_000.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 243_000,
        "obs_tp53_ppm1d_bcor_mutant_molecules_k": 16,
    }
    d_16 = wc_mod.evaluate_stage0_product_contract({**base_scope, **meas_16, **base_mol, "tier6b_sorted_lot_observations": list(clean_t6)}, root=ROOT)
    i_16 = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps({"scope": base_scope, "measurements": meas_16, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": clean_t6}), root=ROOT)
    r_16 = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": base_scope, "measurements": meas_16, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": clean_t6}], root=ROOT)
    s9_atoms_16 = [
        s9_atoms[0],
        wc_mod.AtomicObservation(endpoint_id="obs_tp53_ppm1d_bcor_vaf", scope_type="LOT_SCOPE", scope_key=lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=16.0 / 243_000.0, se_or_n=1e-6, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="S9_VAF16",
                                 controls={"duplex_depth": 243_000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 16}),
        s9_atoms[2],
    ]
    l_16 = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"batch_id": "B1", "atoms": s9_atoms_16, "tier6b_sorted_lot_observations": list(clean_t6)}], target_prov, target_campaign_id=camp_id, root=ROOT)
    for nm, out in (("direct", d_16), ("ingest", i_16), ("replay", r_16), ("ledger", l_16)):
        assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (nm, out)

    # 4a. Negated Tier 6B fraction labels MUST be rejected on ALL 4 paths:
    for bad_enr, bad_ft in (
        ("not_enriched", "not_flow_through"),
        ("non_enriched", "flow_through"),
        ("unenriched", "flow_through"),
        ("enriched_neg", "flow_through"),
        ("enriched", "no_flow_through"),
    ):
        neg_t6 = [
            _mk_t6(bad_enr, "NEG_E", roster=[bad_enr, bad_ft]),
            _mk_t6(bad_ft, "NEG_F", n=100_000, w_i=0.0001, roster=[bad_enr, bad_ft], hpsc=[10, 10], diff=[10, 10]),
        ]
        ev_neg = wc_mod.evaluate_tier6b_per_fraction_thinning_iip0031(neg_t6, target_scope=base_scope, root=ROOT)
        assert ev_neg["status"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (bad_enr, bad_ft, ev_neg)
        assert ev_neg["reason_code"] == "AMBIGUOUS_OR_NEGATED_FRACTION_ID", (bad_enr, bad_ft, ev_neg)

    # 4b. Adverse S11 panel with conflicting campaign_id OR assay_campaign_id MUST hold identically across ALL 4 paths:
    l2_scope = {
        "product_variant_id": pvid_l2, "construct_sha256": sha_l2,
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "donor_seedline_id": seed_id, "mcb_id": mcb_id, "lot_id": lot_id, "campaign_id": camp_id,
        "sample_id": "S11_CAMP_SID", "scope_type": "LOT_SCOPE",
    }
    l2_lot_k = wc_mod.LotScopeKey(pvid_l2, sha_l2, seed_id, mcb_id, lot_id)
    l2_tprov = {**l2_scope, "device_assembly_id": dev_id, "physical_route_id": route_id}
    l2_s9_atoms = [
        wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=l2_lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=4.8, se_or_n=0.10, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="L2_IC9"),
        wc_mod.AtomicObservation(endpoint_id="obs_tp53_ppm1d_bcor_vaf", scope_type="LOT_SCOPE", scope_key=l2_lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="L2_VAF",
                                 controls={"duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0}),
        wc_mod.AtomicObservation(endpoint_id="obs_castseq_translocation_freq", scope_type="LOT_SCOPE", scope_key=l2_lot_k,
                                 kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                                 campaign_id=camp_id, sample_unit_id="L2_CAST",
                                 controls={"duplex_depth": 70000, "obs_castseq_mutant_junctions_k": 0}),
    ]
    l2_clean_t6 = [
        {**r, "product_variant_id": pvid_l2, "construct_sha256": sha_l2, "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL"}
        for r in clean_t6
    ]
    for camp_key in ("campaign_id", "assay_campaign_id"):
        adverse_s11 = {
            "scope_type": "LOT_SCOPE",
            "sample_id": "S11_CAMP_SID",
            camp_key: "OTHER_CAMP_CONFLICT",
            "obs_l2_v3_ap_c5b9_pct": 25.0,
            "obs_l2_v3_ap_c5b9_se": 0.5,
            "replicates_n": 6,
        }
        d_c = wc_mod.evaluate_stage0_product_contract({**l2_scope, **base_meas, **base_mol, "tier6b_sorted_lot_observations": l2_clean_t6, "iip0031_locus2_record": adverse_s11}, root=ROOT)
        i_c = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps({"scope": l2_scope, "measurements": base_meas, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": l2_clean_t6, "iip0031_locus2_record": adverse_s11}), root=ROOT)
        r_c = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": l2_scope, "measurements": base_meas, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": l2_clean_t6, "iip0031_locus2_record": adverse_s11}], root=ROOT)
        s11_atom = wc_mod.AtomicObservation(endpoint_id="stage0_locus2_modeb_v3_panel", scope_type="LOT_SCOPE", scope_key=l2_lot_k,
                                            kind="STRUCTURED_PANEL", value_or_k=adverse_s11, campaign_id=camp_id, sample_unit_id="S11_CAMP_SID")
        l_c = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"batch_id": "B1", "atoms": [*l2_s9_atoms, s11_atom], "tier6b_sorted_lot_observations": l2_clean_t6}], l2_tprov, target_campaign_id=camp_id, root=ROOT)
        for nm, out in (("direct", d_c), ("ingest", i_c), ("replay", r_c), ("ledger", l_c)):
            assert out["s11_evaluation"]["status"] == "HOLD_S11_CONFLICTING_NESTED_EVIDENCE", (camp_key, nm, out)
            assert out["overall_status"] == "HOLD_S11_CONFLICTING_NESTED_EVIDENCE", (camp_key, nm, out)

    # 4c. Decision-tree JSON and Graveyard JSON contain ZERO stale 3.001e7 / 30.01M / 0.05/(2m)-live / >=100-reads wording:
    dt_txt = (ROOT / "data/protocols/stage0_decision_tree_v1.json").read_text(encoding="utf-8")
    gy_txt = (ROOT / "data/graveyard/idea_graveyard_v0.5.json").read_text(encoding="utf-8")
    assert "3.001e7" not in dt_txt and "30.01" not in dt_txt and "informative reads (>= 100)" not in dt_txt
    assert "30.01" not in gy_txt and "3.001e7" not in gy_txt


def test_57_iip0031_round12_metadata_borrowing_unhashable_and_replicate_conflict_parity():
    """
    IIP-0031 Round-12 4-path conformance regressions:
      1. Tier 6B strict Stage-7 pre-assembly timepoint, assay_version, unit, and canonical
         clonogenic_risk_unit allowlists + prohibition on borrowing cal_batch_id or
         untreated_sample_id across distinct treated aliquots (all 4 paths).
      2. Unhashable/malformed JSON values (lists/dicts/oversized ints) in S11 provenance
         fields never raise TypeError/OverflowError and never mask an independent S9 stop.
      3. Contradictory S9 replicate counts across containers (measurements/scope/doc/mol_obs)
         and aliases resolve symmetrically at the conservative minimum valid n for hard stops
         and flag HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT when no hard stop fires.
      4. Playbook single-dose 80k bound reflects live two-fraction alpha_obs = 0.05/6 (5.98419e-5).
    """
    import hashlib
    import json
    import math
    from pathlib import Path
    import scripts.ingest_wetlab_results as ing_mod
    import scripts.product_contract as pc_mod
    import scripts.wetlab_contract as wc_mod

    ROOT = Path(__file__).resolve().parents[1]
    p_doc = pc_mod.load_product_contracts(ROOT)
    p_by_id = {p["product_variant_id"]: p for p in p_doc["products"]}
    a4_sha = hashlib.sha256(b"ISLET1_A4_ARM1B_R12_TEST").hexdigest()
    l2_pvid = "ISLET1_A4_L2_V3_CO_PRIMARY_CD200"
    l2_sha = p_by_id[l2_pvid]["staged_locus2_modeb_donor_binding"]["locus2_donor_sha256"]

    seed_id, mcb_id, lot_id, camp_id = "SEED-R12-01", "MCB-R12-01", "LOT-R12-01", "CAMP-R12-01"
    dev_id, route_id = "DEV-R12-01", "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH"

    base_scope = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "construct_sha256": a4_sha,
        "donor_seedline_id": seed_id,
        "mcb_id": mcb_id,
        "lot_id": lot_id,
        "campaign_id": camp_id,
    }
    base_meas = {
        "obs_icasp9_log10_kill": 4.80,
        "obs_icasp9_log10_kill_se": 0.10,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000,
        "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000,
        "obs_castseq_mutant_junctions_k": 0,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300_000_000,
    }
    base_mol = {
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "per_locus_segmented_baf": {"chr19q": 0.50, "chr11p15": 0.50},
    }
    lot_k = wc_mod.LotScopeKey("ISLET1_A4_ARM1B", a4_sha, seed_id, mcb_id, lot_id)
    target_prov = {**base_scope, "device_assembly_id": dev_id, "physical_route_id": route_id}

    def _s9_atoms(lk, ic9_mean=4.80, ic9_se=0.10, ic9_n=6, ic9_controls=None):
        return [
            wc_mod.AtomicObservation(
                endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lk,
                kind="CONTINUOUS_STUDENT_T", value_or_k=ic9_mean, se_or_n=ic9_se, sample_n=ic9_n,
                campaign_id=camp_id, sample_unit_id="R12_IC9", controls=ic9_controls or {},
            ),
            wc_mod.AtomicObservation(
                endpoint_id="obs_tp53_ppm1d_bcor_vaf", scope_type="LOT_SCOPE", scope_key=lk,
                kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                campaign_id=camp_id, sample_unit_id="R12_VAF",
                controls={"duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0},
            ),
            wc_mod.AtomicObservation(
                endpoint_id="obs_castseq_translocation_freq", scope_type="LOT_SCOPE", scope_key=lk,
                kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                campaign_id=camp_id, sample_unit_id="R12_CAST",
                controls={"duplex_depth": 70000, "obs_castseq_mutant_junctions_k": 0},
            ),
        ]

    def _mk_t6(fid, sid, *, n, w_i, hpsc=(20, 20), diff=(20, 20), extra=None):
        tot = 3_000_000_000
        frac_cells = max(1, int(round(w_i * tot)))
        plated = int((n * 8 + 2) // 3)
        rec = {
            **base_scope,
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "kind": "EXACT_INTEGER_COUNT",
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "fraction_id": fid,
            "declared_fraction_roster": ["enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"],
            "sample_id": sid,
            "k": 0,
            "N": n,
            "w_i": w_i,
            "hpsc_calibrant": list(hpsc),
            "edited_differentiated_calibrant": list(diff),
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "untreated_cfu_control": {
                "untreated_sample_id": f"{sid}_UNTR",
                "untreated_cfu_colonies_k": 4000,
                "untreated_plated_cells_n": 10000,
            },
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": tot,
            "presort_reconciled_loss_cells": 0,
            "presort_fraction_cells": frac_cells,
            "plated_treated_cells_n": plated,
        }
        if extra:
            rec.update(extra)
        return rec

    def _to_t6_atom(rec):
        atom_kwargs = {}
        for ak in ("timepoint", "assay_version", "unit"):
            if ak in rec and isinstance(rec[ak], str) and rec[ak].strip():
                atom_kwargs[ak] = rec[ak]
        ctrl = {
            k: v for k, v in rec.items()
            if k not in (
                "endpoint_id", "scope_type", "kind", "k", "N", "product_variant_id",
                "construct_sha256", "donor_seedline_id", "mcb_id", "lot_id", "campaign_id",
                "timepoint", "assay_version", "unit", "sample_n",
            )
        }
        return wc_mod.AtomicObservation(
            endpoint_id=rec["endpoint_id"],
            scope_type=rec["scope_type"],
            scope_key=lot_k,
            kind=rec["kind"],
            value_or_k=rec["k"],
            se_or_n=rec["N"],
            sample_n=rec.get("sample_n"),
            campaign_id=rec.get("campaign_id", camp_id),
            sample_unit_id=rec["sample_id"],
            controls=ctrl,
            **atom_kwargs,
        )

    def _run_t6_4paths(t6_recs):
        d = wc_mod.evaluate_stage0_product_contract(
            {**base_scope, **base_meas, **base_mol, "tier6b_sorted_lot_observations": list(t6_recs)},
            root=ROOT,
        )
        i = ing_mod.ingest_stage0_research_or_wetlab_batch(
            json.dumps({"scope": base_scope, "measurements": base_meas, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": list(t6_recs)}),
            root=ROOT,
        )
        r = ing_mod.replay_stage0_chronological_batches(
            [{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": base_scope, "measurements": base_meas, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": list(t6_recs)}],
            root=ROOT,
        )
        l = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
            [{"batch_id": "B1", "atoms": [*_s9_atoms(lot_k), *[_to_t6_atom(x) for x in t6_recs]]}],
            target_prov,
            target_campaign_id=camp_id,
            root=ROOT,
        )
        return {"direct": d, "ingest": i, "replay": r, "ledger": l}

    # 1a. Explicit Stage-7 pre-assembly valid metadata passes identically across all 4 paths:
    valid_s7_meta = {"timepoint": "STAGE7_PRE_ASSEMBLY", "assay_version": "TIER6B_SORTED_LOT_V1", "unit": "clonogenic_opportunities"}
    t6_valid_s7 = [
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_OK_E", n=30_000_000, w_i=0.9999, extra=valid_s7_meta),
        _mk_t6("flow_through_and_excluded", "R12_OK_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), extra=valid_s7_meta),
    ]
    for nm, out in _run_t6_4paths(t6_valid_s7).items():
        assert out["tier6b_disposition"] == "PASS", (nm, out)
        assert abs(out["tier6b_sorted_lot_evaluation"]["p_lot_U95"] - 2.1870701159e-7) < 1e-12, (nm, out)

    # 1b. Wrong stage, unvalidated assay version, incompatible unit, and unsupported risk unit MUST HOLD on all 4 paths:
    bad_t6_cases = [
        ("stage2_timepoint", {"timepoint": "STAGE2"}, "INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT"),
        ("post_assembly_timepoint", {"timepoint": "POST_ASSEMBLY"}, "INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT"),
        ("stage2_diff_stage", {"differentiation_stage": "STAGE2"}, "INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT"),
        ("post_assembly_state", {"assembly_state": "POST_ASSEMBLY"}, "INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT"),
        ("unvalidated_assay_version", {"assay_version": "UNVALIDATED_V99"}, "UNVALIDATED_TIER6B_ASSAY_VERSION"),
        ("incompatible_unit", {"unit": "percent"}, "INCOMPATIBLE_TIER6B_MEASUREMENT_UNIT"),
        ("unsupported_whole_lot_risk_unit", {"clonogenic_risk_unit": "WHOLE_LOT_CELLS"}, "INVALID_OR_UNSUBSTANTIATED_CLONOGENIC_RISK_UNIT"),
        ("arbitrary_risk_unit", {"clonogenic_risk_unit": "CUSTOM_UNVERIFIED_UNIT"}, "INVALID_OR_UNSUBSTANTIATED_CLONOGENIC_RISK_UNIT"),
    ]
    for lbl, bad_extra, exp_reason in bad_t6_cases:
        recs_bad = [
            _mk_t6("enriched_cd31_neg_cd144_neg", f"R12_{lbl}_E", n=30_000_000, w_i=0.9999, extra=bad_extra),
            _mk_t6("flow_through_and_excluded", f"R12_{lbl}_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
        ]
        for nm, out in _run_t6_4paths(recs_bad).items():
            assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (lbl, nm, out)
            assert out["tier6b_reason"] == exp_reason, (lbl, nm, out)

    # 1c. Two distinct treated aliquots in the same fraction borrowing cal_batch_id OR untreated_sample_id MUST HOLD on all 4 paths:
    t6_borrow_cal = [
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_BC_E1", n=15_000_000, w_i=0.9999, extra={"cal_batch_id": "SHARED_CAL_01"}),
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_BC_E2", n=15_000_000, w_i=0.9999, extra={"cal_batch_id": "SHARED_CAL_01"}),
        _mk_t6("flow_through_and_excluded", "R12_BC_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for nm, out in _run_t6_4paths(t6_borrow_cal).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert "BORROWED_OR_SHARED_CALIBRATION_BATCH_ACROSS_OBSERVATIONS" in str(out["tier6b_reason"]), (nm, out)

    shared_untr = {"untreated_sample_id": "SHARED_UNTR_01", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000}
    t6_borrow_untr = [
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_BU_E1", n=15_000_000, w_i=0.9999, extra={"untreated_cfu_control": dict(shared_untr)}),
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_BU_E2", n=15_000_000, w_i=0.9999, extra={"untreated_cfu_control": dict(shared_untr)}),
        _mk_t6("flow_through_and_excluded", "R12_BU_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for nm, out in _run_t6_4paths(t6_borrow_untr).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert "BORROWED_OR_SHARED_UNTREATED_CFU_CONTROL_ACROSS_OBSERVATIONS" in str(out["tier6b_reason"]), (nm, out)

    t6_indep_two_aliquots = [
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_IND_E1", n=15_000_000, w_i=0.9999, extra={"record_id": "R12_REC_A", "source_sample_id": "R12_SRC_A", "cal_batch_id": "IND_CAL_A"}),
        _mk_t6("enriched_cd31_neg_cd144_neg", "R12_IND_E2", n=15_000_000, w_i=0.9999, extra={"record_id": "R12_REC_B", "source_sample_id": "R12_SRC_B", "cal_batch_id": "IND_CAL_B"}),
        _mk_t6("flow_through_and_excluded", "R12_IND_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    for nm, out in _run_t6_4paths(t6_indep_two_aliquots).items():
        assert out["tier6b_disposition"] == "PASS", (nm, out)

    # 2. Malformed/unhashable JSON types in S11 provenance fields NEVER raise TypeError and NEVER mask an S9 hard stop:
    l2_scope = {
        "product_variant_id": l2_pvid,
        "construct_sha256": l2_sha,
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "donor_seedline_id": seed_id,
        "mcb_id": mcb_id,
        "lot_id": lot_id,
        "campaign_id": camp_id,
        "sample_id": "R12_S11_SID",
        "scope_type": "LOT_SCOPE",
    }
    l2_lot_k = wc_mod.LotScopeKey(l2_pvid, l2_sha, seed_id, mcb_id, lot_id)
    l2_tprov = {**l2_scope, "device_assembly_id": dev_id, "physical_route_id": route_id}
    s9_fail_meas = {**base_meas, "obs_icasp9_log10_kill": 3.50, "obs_icasp9_log10_kill_se": 0.10}
    for bad_key in ("donor_seedline_id", "seed_line_id", "l2_v3_donor_id", "donor_id", "background_id", "sample_id", "lot_id", "campaign_id", "mcb_id", "construct_sha256"):
        for bad_val in (["x"], {"nested": 1}):
            bad_s11_rec = {"scope_type": "LOT_SCOPE", bad_key: bad_val}
            d_u = wc_mod.evaluate_stage0_product_contract({**l2_scope, **s9_fail_meas, **base_mol, "iip0031_locus2_record": bad_s11_rec}, root=ROOT)
            i_u = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps({"scope": l2_scope, "measurements": s9_fail_meas, "iip0026_molecular_observations": base_mol, "iip0031_locus2_record": bad_s11_rec}), root=ROOT)
            r_u = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": l2_scope, "measurements": s9_fail_meas, "iip0026_molecular_observations": base_mol, "iip0031_locus2_record": bad_s11_rec}], root=ROOT)
            s11_a = wc_mod.AtomicObservation(
                endpoint_id="stage0_locus2_modeb_v3_panel", scope_type="LOT_SCOPE", scope_key=l2_lot_k,
                kind="STRUCTURED_PANEL", value_or_k=bad_s11_rec, campaign_id=camp_id, sample_unit_id="R12_S11_SID",
            )
            l_u = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
                [{"batch_id": "B_U", "atoms": [*_s9_atoms(l2_lot_k, ic9_mean=3.50, ic9_se=0.10), s11_a]}],
                l2_tprov, target_campaign_id=camp_id, root=ROOT,
            )
            for nm, out in (("direct", d_u), ("ingest", i_u), ("replay", r_u), ("ledger", l_u)):
                assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (bad_key, bad_val, nm, out)
                assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (bad_key, bad_val, nm, out)
            # Without S9 hard stop, malformed nested provenance also fails closed cleanly without raising:
            clean_direct = wc_mod.evaluate_stage0_product_contract({**l2_scope, **base_meas, **base_mol, "iip0031_locus2_record": bad_s11_rec}, root=ROOT)
            assert clean_direct["s11_evaluation"]["gate_pass"] is False, (bad_key, bad_val, clean_direct)
    # Outer scope list in donor_seedline_id (direct/ingest/replay) alongside S9 hard stop also never raises:
    outer_list_scope = {**l2_scope, "donor_seedline_id": ["x"]}
    for nm, out in (
        ("direct", wc_mod.evaluate_stage0_product_contract({**outer_list_scope, **s9_fail_meas, **base_mol}, root=ROOT)),
        ("ingest", ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps({"scope": outer_list_scope, "measurements": s9_fail_meas, "iip0026_molecular_observations": base_mol}), root=ROOT)),
        ("replay", ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": outer_list_scope, "measurements": s9_fail_meas, "iip0026_molecular_observations": base_mol}], root=ROOT)),
    ):
        assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (nm, out)
        assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (nm, out)

    # 3. Contradictory S9 replicate counts (100 vs 6) across containers & aliases are order-independent on ALL 4 paths:
    #    - Adverse iCasp9=4.18, se=0.10 (passes at n=100 [4.014>=4.0], hard-stops at n=6 [3.9785<4.0]) -> HARD_SAFETY_STOP_REJECTED in BOTH orders!
    #    - Clean iCasp9=4.80, se=0.10 -> HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT in BOTH orders (even when tp53/castseq are omitted)!
    for key_m, key_s in (
        ("obs_icasp9_log10_kill_n", "obs_icasp9_log10_kill_n"),
        ("replicates_n", "replicates_n"),
        ("obs_icasp9_log10_kill_n", "obs_icasp9_log10_kill_sample_n"),
        ("replicates_n", "sample_n"),
    ):
        for n_meas, n_scope in ((100, 6), (6, 100)):
            # 3a. Adverse partial submission (only iCasp9 supplied, missing tp53/castseq):
            partial_doc_adv = {
                "scope": {**base_scope, key_s: n_scope},
                "measurements": {"obs_icasp9_log10_kill": 4.18, "obs_icasp9_log10_kill_se": 0.10, key_m: n_meas},
            }
            ing_adv = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(partial_doc_adv), root=ROOT)
            rep_adv = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, **partial_doc_adv}], root=ROOT)
            assert ing_adv["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (key_m, key_s, n_meas, n_scope, ing_adv)
            assert rep_adv["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (key_m, key_s, n_meas, n_scope, rep_adv)

            # 3b. Clean partial submission (iCasp9=4.80, missing tp53/castseq) -> MUST flag conflict explicitly, never hide behind missing-core only:
            partial_doc_ok = {
                "scope": {**base_scope, key_s: n_scope},
                "measurements": {"obs_icasp9_log10_kill": 4.80, "obs_icasp9_log10_kill_se": 0.10, key_m: n_meas},
            }
            ing_ok = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(partial_doc_ok), root=ROOT)
            rep_ok = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, **partial_doc_ok}], root=ROOT)
            assert ing_ok["overall_status"] == "HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT", (key_m, key_s, n_meas, n_scope, ing_ok)
            assert rep_ok["overall_status"] == "HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT", (key_m, key_s, n_meas, n_scope, rep_ok)

            # 3c. Clean full submission with contradictory replicate counts -> HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT:
            full_doc_ok = {
                "scope": {**base_scope, key_s: n_scope},
                "measurements": {**base_meas, key_m: n_meas},
                "iip0026_molecular_observations": base_mol,
                "tier6b_sorted_lot_observations": t6_valid_s7,
            }
            ing_f = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(full_doc_ok), root=ROOT)
            rep_f = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, **full_doc_ok}], root=ROOT)
            assert ing_f["overall_status"] == "HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT", (key_m, key_s, n_meas, n_scope, ing_f)
            assert rep_f["overall_status"] == "HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT", (key_m, key_s, n_meas, n_scope, rep_f)

    # 3d. Direct conflicting aliases & Atomic Ledger conflicting controls match identically:
    for d_n1, d_n2 in ((100, 6), (6, 100)):
        d_adv = wc_mod.evaluate_stage0_product_contract(
            {**base_scope, "obs_icasp9_log10_kill": 4.18, "obs_icasp9_log10_kill_se": 0.10, "obs_icasp9_log10_kill_n": d_n1, "obs_icasp9_log10_kill_sample_n": d_n2},
            root=ROOT,
        )
        assert d_adv["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (d_n1, d_n2, d_adv)
        d_ok = wc_mod.evaluate_stage0_product_contract(
            {**base_scope, "obs_icasp9_log10_kill": 4.80, "obs_icasp9_log10_kill_se": 0.10, "obs_icasp9_log10_kill_n": d_n1, "obs_icasp9_log10_kill_sample_n": d_n2},
            root=ROOT,
        )
        assert d_ok["overall_status"] == "HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT", (d_n1, d_n2, d_ok)
    # Oversized int (10**400) in replicate count never raises OverflowError:
    d_huge = wc_mod.evaluate_stage0_product_contract(
        {**base_scope, "obs_icasp9_log10_kill": 4.18, "obs_icasp9_log10_kill_se": 0.10, "obs_icasp9_log10_kill_n": 100, "obs_icasp9_log10_kill_sample_n": 10**400},
        root=ROOT,
    )
    assert d_huge["overall_status"] == "HARD_SAFETY_STOP_REJECTED", d_huge
    # Atomic Ledger conflicting controls on iCasp9 atom:
    l_adv = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
        [{"batch_id": "B_CONFL_ADV", "atoms": [*_s9_atoms(lot_k, ic9_mean=4.18, ic9_se=0.10, ic9_n=100, ic9_controls={"obs_icasp9_log10_kill_n": 6})]}],
        target_prov, target_campaign_id=camp_id, root=ROOT,
    )
    assert l_adv["overall_status"] == "HARD_SAFETY_STOP_REJECTED", l_adv
    l_ok = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
        [{"batch_id": "B_CONFL_OK", "atoms": [*_s9_atoms(lot_k, ic9_mean=4.80, ic9_se=0.10, ic9_n=100, ic9_controls={"obs_icasp9_log10_kill_n": 6}), *[_to_t6_atom(x) for x in t6_valid_s7]]}],
        target_prov, target_campaign_id=camp_id, root=ROOT,
    )
    assert any("HOLD_S9_INVALID_OR_CONFLICTING_REPLICATE_COUNT" in h for h in l_ok["holds"]), l_ok

    # 4. Playbook single-dose bound uses live alpha_obs = 0.05/6 (5.98419e-5), never stale 5.48e-5:
    pb_txt = (ROOT / "docs/STAGE-0-EXECUTION-PLAYBOOK.md").read_text(encoding="utf-8")
    assert "5.98419" in pb_txt and "5.48 \\times 10^{-5}" not in pb_txt
    assert abs((-math.expm1(math.log(0.05 / 6.0) / 80_000.0)) - 5.9841856189e-5) < 1e-12


def test_58_iip0031_round13_duplicate_record_id_explicit_metadata_null_conflict_and_22gene_pins():
    """
    IIP-0031 Round-13 4-path conformance regressions:
      1. Copying a Tier 6B observation with the same `record_id` (or alias) and relabeled
         treated/untreated/calibrant identifiers MUST NOT earn extra trial credit on ANY path
         (`CONFLICTING_OR_RELABELED_DUPLICATE_RECORD_ID` -> `HOLD_TIER6B_UNQUALIFIED_FRACTION`).
      2. Omitting explicit Stage-7 pre-assembly `timepoint`, `assay_version`, or `unit` MUST NOT
         receive generic-default clearance credit on ANY path (`HOLD_TIER6B_UNQUALIFIED_FRACTION`).
      3. An adverse S9 record paired with a nested or top-level `null` for the same mean/SE OR
         with cross-group replicate-count conflict (`replicates_n=100` vs `obs_icasp9_log10_kill_n=6`
         in either order) MUST retain `HARD_SAFETY_STOP_REJECTED` (`S9_HARD_SAFETY_STOP`) across
         ALL 4 paths (and non-adverse conflicts MUST hold rather than clearing).
      4. The 22-gene RefSeq off-target panel pins accessions, SHA-256 digests, ATG..STOP CDS
         bounds, and nonempty 3'-UTRs, rejecting out-of-range slices, swapped accessions, or
         mutated sequences while verifying true CDS 7mer-m8 hits separately from 5'-UTR/full-mRNA.
    """
    import copy
    import hashlib
    import json
    from pathlib import Path
    import pytest
    import scripts.ingest_wetlab_results as ing_mod
    import scripts.product_contract as pc_mod
    import scripts.run_layer2_perturbation_and_regulatory_plane as l2_mod
    import scripts.wetlab_contract as wc_mod

    ROOT = Path(__file__).resolve().parents[1]
    p_doc = pc_mod.load_product_contracts(ROOT)
    p_by_id = {p["product_variant_id"]: p for p in p_doc["products"]}
    a4_sha = hashlib.sha256(b"ISLET1_A4_ARM1B_R13_TEST").hexdigest()
    l2_pvid = "ISLET1_A4_L2_V3_CO_PRIMARY_CD200"
    l2_sha = p_by_id[l2_pvid]["staged_locus2_modeb_donor_binding"]["locus2_donor_sha256"]

    seed_id, mcb_id, lot_id, camp_id = "SEED-R13-01", "MCB-R13-01", "LOT-R13-01", "CAMP-R13-01"
    dev_id, route_id = "DEV-R13-01", "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH"
    base_scope = {
        "product_variant_id": l2_pvid,
        "construct_sha256": l2_sha,
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "donor_seedline_id": seed_id,
        "mcb_id": mcb_id,
        "lot_id": lot_id,
        "campaign_id": camp_id,
    }
    base_meas = {
        "obs_icasp9_log10_kill": 4.80,
        "obs_icasp9_log10_kill_se": 0.10,
        "obs_tp53_ppm1d_bcor_vaf": 0.0,
        "obs_tp53_ppm1d_bcor_vaf_se": 1e-6,
        "obs_tp53_ppm1d_bcor_duplex_depth": 70000,
        "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0,
        "obs_castseq_translocation_freq": 0.0,
        "obs_castseq_translocation_se": 1e-6,
        "obs_castseq_duplex_or_amplicon_depth": 70000,
        "obs_castseq_mutant_junctions_k": 0,
        "obs_mcb_clonogenic_survivor_count_k": 0,
        "obs_mcb_effective_clonogenic_trials": 300_000_000,
    }
    base_mol = {
        "donor_seed_gnomad_target_site_wgs_verified": True,
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "per_locus_segmented_baf": {"chr19q": 0.50, "chr11p15": 0.50},
    }
    lot_k = wc_mod.LotScopeKey(l2_pvid, l2_sha, seed_id, mcb_id, lot_id)
    target_prov = {**base_scope, "device_assembly_id": dev_id, "physical_route_id": route_id}

    def _s9_atoms(lk, ic9_mean=4.80, ic9_se=0.10, ic9_n=6, ic9_controls=None):
        return [
            wc_mod.AtomicObservation(
                endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lk,
                kind="CONTINUOUS_STUDENT_T", value_or_k=ic9_mean, se_or_n=ic9_se, sample_n=ic9_n,
                campaign_id=camp_id, sample_unit_id="R13_IC9", controls=ic9_controls or {},
            ),
            wc_mod.AtomicObservation(
                endpoint_id="obs_tp53_ppm1d_bcor_vaf", scope_type="LOT_SCOPE", scope_key=lk,
                kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                campaign_id=camp_id, sample_unit_id="R13_VAF",
                controls={"duplex_depth": 70000, "obs_tp53_ppm1d_bcor_mutant_molecules_k": 0},
            ),
            wc_mod.AtomicObservation(
                endpoint_id="obs_castseq_translocation_freq", scope_type="LOT_SCOPE", scope_key=lk,
                kind="CONTINUOUS_STUDENT_T", value_or_k=0.0, se_or_n=1e-6, sample_n=6,
                campaign_id=camp_id, sample_unit_id="R13_CAST",
                controls={"duplex_depth": 70000, "obs_castseq_mutant_junctions_k": 0},
            ),
        ]

    def _mk_t6_r13(fid, sid, *, n, w_i, hpsc=(20, 20), diff=(20, 20), extra=None, omit=()):
        tot = 3_000_000_000
        frac_cells = 2_999_700_000 if "enriched" in fid else 300_000
        plated = int((n * 8 + 2) // 3)
        rec = {
            **base_scope,
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "kind": "EXACT_INTEGER_COUNT",
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "fraction_id": fid,
            "declared_fraction_roster": ["enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"],
            "sample_id": sid,
            "k": 0,
            "N": n,
            "w_i": w_i,
            "hpsc_calibrant": list(hpsc),
            "edited_differentiated_calibrant": list(diff),
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "untreated_cfu_control": {
                "untreated_sample_id": f"{sid}_UNTR",
                "untreated_cfu_colonies_k": 4000,
                "untreated_plated_cells_n": 10000,
            },
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": tot,
            "presort_reconciled_loss_cells": 0,
            "presort_fraction_cells": frac_cells,
            "plated_treated_cells_n": plated,
        }
        if extra:
            rec.update(extra)
        for ok in omit:
            rec.pop(ok, None)
        return rec

    def _to_atom_r13(rec):
        akw = {}
        for k in ("timepoint", "assay_version", "unit"):
            if k in rec and isinstance(rec[k], str) and rec[k].strip():
                akw[k] = rec[k]
        ctrl = {
            k: v for k, v in rec.items()
            if k not in (
                "endpoint_id", "scope_type", "kind", "k", "N",
                "product_variant_id", "construct_sha256", "donor_seedline_id",
                "mcb_id", "lot_id", "campaign_id", "timepoint", "assay_version",
                "unit", "sample_n",
            )
        }
        return wc_mod.AtomicObservation(
            endpoint_id=rec["endpoint_id"],
            scope_type=rec.get("scope_type", "LOT_SCOPE"),
            scope_key=lot_k,
            kind=rec.get("kind", "EXACT_INTEGER_COUNT"),
            value_or_k=rec["k"],
            se_or_n=rec["N"],
            campaign_id=camp_id,
            sample_unit_id=rec["sample_id"],
            controls=ctrl,
            **akw,
        )

    def _eval_4_r13(t6_recs):
        d = wc_mod.evaluate_stage0_product_contract(
            {**base_scope, **base_meas, **base_mol, "tier6b_sorted_lot_observations": list(t6_recs)},
            root=ROOT,
        )
        i = ing_mod.ingest_stage0_research_or_wetlab_batch(
            json.dumps({"scope": base_scope, "measurements": base_meas, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": list(t6_recs)}),
            root=ROOT,
        )
        r = ing_mod.replay_stage0_chronological_batches(
            [{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, "scope": base_scope, "measurements": base_meas, "iip0026_molecular_observations": base_mol, "tier6b_sorted_lot_observations": list(t6_recs)}],
            root=ROOT,
        )
        l = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
            [{"batch_id": "B_R13", "atoms": [*_s9_atoms(lot_k), *[_to_atom_r13(x) for x in t6_recs]]}],
            target_prov,
            target_campaign_id=camp_id,
            root=ROOT,
        )
        return {"direct": d, "ingest": i, "replay": r, "ledger": l}

    # 1. Same-source `record_id` copy with relabeled sample/untreated/calibrant IDs MUST NOT double N on ANY path:
    underpowered_base = [
        _mk_t6_r13("enriched_cd31_neg_cd144_neg", "R13_ENR_1", n=18_000_000, w_i=0.9999, extra={"record_id": "SRC_REC_001", "cal_batch_id": "CAL_E1"}),
        _mk_t6_r13("flow_through_and_excluded", "R13_FT_1", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), extra={"record_id": "SRC_REC_002", "cal_batch_id": "CAL_F1"}),
    ]
    for nm, out in _eval_4_r13(underpowered_base).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_INSUFFICIENT_DENOMINATOR", (nm, out)
        assert abs(out["tier6b_sorted_lot_evaluation"]["p_lot_U95"] - 3.586205214e-7) < 1e-12, (nm, out)
        assert out["s11_evaluation"]["status"] == "HOLD_S11_INCOMPLETE_EVIDENCE", (nm, out)

    copied_same_rec_id = [
        underpowered_base[0],
        _mk_t6_r13("enriched_cd31_neg_cd144_neg", "R13_ENR_2_RELABELED", n=18_000_000, w_i=0.9999, extra={"record_id": "SRC_REC_001", "cal_batch_id": "CAL_E2_RELABELED"}),
        underpowered_base[1],
    ]
    for nm, out in _eval_4_r13(copied_same_rec_id).items():
        assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (nm, out)
        assert "CONFLICTING_OR_RELABELED_DUPLICATE_RECORD_ID" in str(out["tier6b_reason"]), (nm, out)
        assert out["s11_evaluation"]["status"] == "HOLD_S11_INCOMPLETE_EVIDENCE", (nm, out)

    # 2. Missing explicit `timepoint`, `assay_version`, or `unit` MUST HOLD on all 4 paths (no generic defaults):
    for omitted_fields, expected_reason in (
        (("timepoint",), "INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT"),
        (("assay_version",), "UNVALIDATED_TIER6B_ASSAY_VERSION"),
        (("unit",), "INCOMPATIBLE_TIER6B_MEASUREMENT_UNIT"),
        (("timepoint", "assay_version", "unit"), "INVALID_OR_WRONG_STAGE_TIER6B_TIMEPOINT"),
    ):
        missing_meta_recs = [
            _mk_t6_r13("enriched_cd31_neg_cd144_neg", "R13_MM_E", n=30_000_000, w_i=0.9999, omit=omitted_fields),
            _mk_t6_r13("flow_through_and_excluded", "R13_MM_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10), omit=omitted_fields),
        ]
        for nm, out in _eval_4_r13(missing_meta_recs).items():
            assert out["tier6b_disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (omitted_fields, nm, out)
            assert out["tier6b_reason"] == expected_reason, (omitted_fields, nm, out)

    # 3a. Adverse S9 paired with nested/top-level `null` for the same mean/SE MUST retain S9_HARD_SAFETY_STOP on Direct, Ingest, and Replay:
    valid_t6_r13 = [
        _mk_t6_r13("enriched_cd31_neg_cd144_neg", "R13_OK_E", n=30_000_000, w_i=0.9999),
        _mk_t6_r13("flow_through_and_excluded", "R13_OK_F", n=100_000, w_i=0.0001, hpsc=(10, 10), diff=(10, 10)),
    ]
    null_adverse_packets = [
        {
            "scope": base_scope,
            "obs_icasp9_log10_kill": 3.50,
            "obs_icasp9_log10_kill_se": 0.10,
            "measurements": {"obs_icasp9_log10_kill": None},
        },
        {
            "scope": base_scope,
            "obs_icasp9_log10_kill": None,
            "measurements": {"obs_icasp9_log10_kill": 3.50, "obs_icasp9_log10_kill_se": 0.10},
        },
        {
            "scope": {**base_scope, "obs_icasp9_log10_kill": None},
            "measurements": {"obs_icasp9_log10_kill": 3.50, "obs_icasp9_log10_kill_se": 0.10},
        },
    ]
    for pkt in null_adverse_packets:
        d_na = wc_mod.evaluate_stage0_product_contract(pkt, root=ROOT)
        i_na = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(pkt), root=ROOT)
        r_na = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, **pkt}], root=ROOT)
        for nm, out in (("direct", d_na), ("ingest", i_na), ("replay", r_na)):
            assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (nm, pkt, out)
            assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (nm, pkt, out)
            assert out["s11_evaluation"]["status"].startswith("HOLD_S11_"), (nm, pkt, out)

    # Non-adverse clean S9 paired with nested `null` MUST hold (`HOLD_CONFLICTING_STAGE0_MEASUREMENT`) and never clear:
    a4_scope = {
        "product_variant_id": "ISLET1_A4_ARM1B",
        "construct_sha256": a4_sha,
        "donor_seedline_id": seed_id,
        "mcb_id": mcb_id,
        "lot_id": lot_id,
        "campaign_id": camp_id,
    }
    a4_t6 = [{**r, "product_variant_id": "ISLET1_A4_ARM1B", "construct_sha256": a4_sha} for r in valid_t6_r13]
    for r_x in a4_t6:
        r_x.pop("background_id", None)
    clean_plus_null_pkt = {
        "scope": a4_scope,
        "obs_icasp9_log10_kill": 4.80,
        "measurements": {**base_meas, "obs_icasp9_log10_kill": None},
        "iip0026_molecular_observations": base_mol,
        "tier6b_sorted_lot_observations": a4_t6,
    }
    for nm, out in (
        ("direct", wc_mod.evaluate_stage0_product_contract(clean_plus_null_pkt, root=ROOT)),
        ("ingest", ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(clean_plus_null_pkt), root=ROOT)),
        ("replay", ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, **clean_plus_null_pkt}], root=ROOT)),
    ):
        assert out["overall_status"] == "HOLD_CONFLICTING_STAGE0_MEASUREMENT", (nm, out)
        assert out["effective_execution_permitted"] is False, (nm, out)

    # 3b. Cross-group replicate-count conflict (`replicates_n` vs `obs_icasp9_log10_kill_n`: 100 vs 6) MUST produce S9_HARD_SAFETY_STOP on ALL 4 paths:
    for n_def, n_ep in ((100, 6), (6, 100)):
        pkt_xgrp = {
            "scope": {**base_scope, "replicates_n": n_def},
            "measurements": {"obs_icasp9_log10_kill": 4.18, "obs_icasp9_log10_kill_se": 0.10, "obs_icasp9_log10_kill_n": n_ep},
        }
        d_xg = wc_mod.evaluate_stage0_product_contract(pkt_xgrp, root=ROOT)
        i_xg = ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(pkt_xgrp), root=ROOT)
        r_xg = ing_mod.replay_stage0_chronological_batches([{"utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 0, **pkt_xgrp}], root=ROOT)
        l_xg = wc_mod.evaluate_iip0024_atomic_stage0_ledger(
            [{"batch_id": "B_XG", "replicates_n": n_def, "atoms": [*_s9_atoms(lot_k, ic9_mean=4.18, ic9_se=0.10, ic9_n=n_ep, ic9_controls={"replicates_n": n_def})]}],
            target_prov, target_campaign_id=camp_id, root=ROOT,
        )
        for nm, out in (("direct", d_xg), ("ingest", i_xg), ("replay", r_xg), ("ledger", l_xg)):
            assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED", (n_def, n_ep, nm, out)
            assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (n_def, n_ep, nm, out)
            assert out["s11_evaluation"]["status"].startswith("HOLD_S11_"), (n_def, n_ep, nm, out)

    # 4. 22-gene RefSeq off-target panel pins accessions, CDS bounds, SHA-256s, and nonempty 3'-UTRs:
    raw_ref = json.loads((ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text(encoding="utf-8"))
    ref_b = raw_ref["iip0031_locus2_modeb_v3_reference_bundle"]
    u3, cds_m, u5, full_m = l2_mod._verified_22gene_refseq_panel(ref_b)
    assert len(u3) == 22 and all(len(v) >= 1 for v in u3.values())
    assert len(cds_m) == 22 and all(v.startswith("ATG") and v[-3:] in ("TAA", "TAG", "TGA") for v in cds_m.values())

    bad_utr_b = copy.deepcopy(ref_b)
    bad_utr_b["human_22_gene_refseq_panel"]["INS"]["cds_end_1based"] = len(bad_utr_b["human_22_gene_refseq_panel"]["INS"]["mrna_sequence"])
    with pytest.raises(ValueError):
        l2_mod._verified_22gene_refseq_panel(bad_utr_b)

    bad_acc_b = copy.deepcopy(ref_b)
    bad_acc_b["human_22_gene_refseq_panel"]["INS"]["accession"] = "NM_000000.0"
    with pytest.raises(ValueError):
        l2_mod._verified_22gene_refseq_panel(bad_acc_b)

    bad_sha_b = copy.deepcopy(ref_b)
    bad_sha_b["human_22_gene_refseq_panel"]["INS"]["mrna_sequence"] = "ATG" + bad_sha_b["human_22_gene_refseq_panel"]["INS"]["mrna_sequence"][3:]
    with pytest.raises(ValueError):
        l2_mod._verified_22gene_refseq_panel(bad_sha_b)


def test_59_iip0031_hardened_4path_tier6b_s11_layer3_conformance():
    """
    Comprehensive conformance suite for IIP-0031 hardened 4-path evidence validation,
    Tier 6B source-owned multi-observation identity & exact integer/CP bounds, Rule S11
    provenance/ADCC/CIITA/RACE isolation, and Layer-3 2WII DBREF Gln1013 (33.04 A) / 15-PDB inventory.
    """
    import copy
    import json
    from types import MappingProxyType
    from pathlib import Path
    from scripts import ingest_wetlab_results as ing_mod
    from scripts import product_contract as pc_mod
    from scripts import wetlab_contract as wc_mod

    ROOT = Path(__file__).resolve().parents[1]
    l2_pdef = pc_mod.get_product_definition("L2_V3_CO_PRIMARY_CD200", root=ROOT)
    l2_pvid = l2_pdef["product_variant_id"]
    l2_sha = l2_pdef["staged_locus2_modeb_donor_binding"]["locus2_donor_sha256"]
    bsct_scope = {
        "product_variant_id": "ISLET1_B_SCT_C1_CIITA",
        "construct_sha256": "7" * 64,
        "donor_seedline_id": "S_QA",
        "mcb_id": "M_QA",
        "lot_id": "L_QA",
        "campaign_id": "C_QA",
        "device_assembly_id": "D_QA",
        "physical_route_id": "SUPRAFASCIAL_350UM_CAO2_CATALASE_POUCH",
    }
    l2_scope = {**bsct_scope, "product_variant_id": l2_pvid, "construct_sha256": l2_sha, "scope_type": "LOT_SCOPE", "sample_id": "S_L2"}
    s9_adv = {"obs_icasp9_log10_kill": 0.1, "obs_icasp9_log10_kill_se": 0.1}
    t6_pos = {**bsct_scope, "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape", "scope_type": "LOT_SCOPE", "sample_id": "T6_POS", "k": 1}
    lectin_fail = {"obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.1, "replicates_n": 6}
    classical_fail = {"obs_l2_v3_classical_c4b2a_lysis_pct": 18.0, "obs_l2_v3_classical_c4b2a_lysis_se": 0.1, "replicates_n": 6}

    # 1. Malformed container + 10**400 scalar never suppress independently attributable S9 / Tier 6B / S11 hard stops across all 4 dict/wire/replay paths:
    for exp_w, doc in (
        ("S9_HARD_SAFETY_STOP", {**bsct_scope, **s9_adv, "measurements": []}),
        ("S9_HARD_SAFETY_STOP", {"scope": bsct_scope, "measurements": {**s9_adv, "obs_stage7_yield_ratio": 10**400}}),
        ("FAIL_TIER6B_OUTGROWTH_DETECTED", {**bsct_scope, "tier6b_sorted_lot_observations": [t6_pos], "scope": []}),
        ("FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", {"scope": l2_scope, "iip0031_locus2_record": None, "iip0026_molecular_observations": {"iip0031_locus2_record": lectin_fail}}),
        ("FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", {"scope": l2_scope, "iip0031_locus2_record": {}, "iip0026_molecular_observations": {"iip0031_locus2_record": lectin_fail}}),
        ("FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", {"scope": l2_scope, "iip0031_locus2_record": 10**400, "iip0026_molecular_observations": {"iip0031_locus2_record": lectin_fail}}),
    ):
        for fn in (
            lambda d: wc_mod.evaluate_stage0_product_contract(d, root=ROOT),
            lambda d: ing_mod.ingest_stage0_research_or_wetlab_batch(d, root=ROOT),
            lambda d: ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(d), root=ROOT),
            lambda d: ing_mod.replay_stage0_chronological_batches([d], root=ROOT),
        ):
            out = fn(doc)
            assert out["hard_stop_precedence_winner"] == exp_w, (exp_w, out)
            assert out["effective_execution_permitted"] is False

    # 2. Mixed-scope / cross-arm / ambiguous-envelope / incomplete-scope replay preserves attributable hard stops in both orders:
    adv_s9_pkt = {"scope": bsct_scope, "measurements": s9_adv}
    for peer in (
        {"scope": {**bsct_scope, "lot_id": "OTHER"}},
        {"scope": bsct_scope, "lot_id": None},
        {"scope": bsct_scope, "adjunct_conditioning_history_arm": "ARM2_CONDITIONED_DEPLETED"},
    ):
        a_use = {**adv_s9_pkt, "adjunct_conditioning_history_arm": "ARM1_NEVER_CONDITIONED"} if "adjunct_conditioning_history_arm" in peer else adv_s9_pkt
        for order in ([a_use, peer], [peer, a_use]):
            r_mix = ing_mod.replay_stage0_chronological_batches(order, root=ROOT)
            assert r_mix["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", (peer, r_mix)
            assert r_mix["effective_execution_permitted"] is False
    inc_s = {k: v for k, v in bsct_scope.items() if k != "campaign_id"}
    for order in ([{"scope": inc_s, "measurements": s9_adv}, {"scope": inc_s, "measurements": {"obs_icasp9_log10_kill": 10**400}}],
                  [{"scope": inc_s, "measurements": {"obs_icasp9_log10_kill": 5.0, "obs_icasp9_log10_kill_se": 0.1}}, {"scope": inc_s, "measurements": s9_adv}]):
        r_inc = ing_mod.replay_stage0_chronological_batches(order, root=ROOT)
        assert r_inc["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", r_inc

    # 3. ISO-8601 UTC timestamp ordering & tie determinism + S11 precedence over generic failures:
    t_before = {"scope": l2_scope, "iip0031_locus2_record": lectin_fail, "utc_timestamp": "2026-11-20T00:00:00Z", "seq_index": 1}
    t_after = {"scope": l2_scope, "iip0031_locus2_record": classical_fail, "utc_timestamp": "2026-11-20T00:00:00.001Z", "seq_index": 2}
    for order in ([t_before, t_after], [t_after, t_before]):
        r_ts = ing_mod.replay_stage0_chronological_batches(order, root=ROOT)
        assert r_ts["overall_status"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", r_ts
        assert "FAIL_S11_POST_CUFF_CLASSICAL_C4B2A_LYSIS_EXCEEDS_10PCT" in r_ts["s11_evaluation"]["reasons"]
    tie_a = {"scope": l2_scope, "iip0031_locus2_record": lectin_fail, "seq_index": 1}
    tie_b = {"scope": l2_scope, "iip0031_locus2_record": classical_fail, "seq_index": 1}
    r_tie1 = ing_mod.replay_stage0_chronological_batches([tie_a, tie_b], root=ROOT)
    r_tie2 = ing_mod.replay_stage0_chronological_batches([tie_b, tie_a], root=ROOT)
    assert r_tie1["overall_status"] == r_tie2["overall_status"] == "FAIL_S11_POST_CUFF_CLASSICAL_C4B2A_LYSIS_EXCEEDS_10PCT"
    s11_bsct = {"scope": bsct_scope, "iip0031_locus2_record": {"test": True}}
    ovr_pkt = {"scope": bsct_scope, "execution_supported": False}
    for order in ([ovr_pkt, s11_bsct], [s11_bsct, ovr_pkt]):
        r_p7 = ing_mod.replay_stage0_chronological_batches(order, root=ROOT)
        assert r_p7["overall_status"] == "FAIL_S11_TRACK2_BSCT_INCOMPATIBLE_AT_LOCUS_2", r_p7

    # 4. Shared BAF / Design-B / TYK2 alias conflict preservation across all 4 containers & paths:
    b2_tyk2_base = {
        "scope": {**bsct_scope, "product_variant_id": "ISLET1_B2_5LOCUS_TYK2_PE"},
        "obs_icasp9_log10_kill": 5.0, "obs_icasp9_log10_kill_se": 0.1,
        "obs_tp53_ppm1d_bcor_vaf": 0.0, "obs_tp53_ppm1d_bcor_vaf_se": 1e-6, "obs_tp53_ppm1d_bcor_duplex_depth": 70000,
        "obs_castseq_translocation_freq": 0.0, "obs_castseq_translocation_se": 1e-6, "obs_castseq_duplex_or_amplicon_depth": 70000,
        "obs_mcb_clonogenic_survivor_count_k": 0, "obs_mcb_effective_clonogenic_trials": 300000000,
        "bound_c1_c5_pair_sha256": wc_mod.DESIGN_B_C1_C5_PAIR_SHA256,
        "post_edit_reciprocal_mosaic_sv_or_clonal_zygosity_verified": True,
        "design_b_locus5_race_umi_n_eff": 4000, "tyk2_pe_modality": "PE4MAX_MSI_GATED", "tyk2_phased_haplotype_sap_verified": True,
    }
    for doc_a, exp_reason in (
        ({**b2_tyk2_base, "measurements": {"per_locus_segmented_baf": {"chr19q": 0.25, "chr11p15": 0.25, "chr6p21": 0.25}}}, "FAIL_DISCORDANT_PER_LOCUS_SEGMENTED_BAF"),
        ({**b2_tyk2_base, "design_b_locus5_aberrant_umi_k": 1, "obs_design_b_locus5_aberrant_umi_k": 0}, "FAIL_DESIGNB_LOCUS5_CHIMERIC_TRANSCRIPT_DETECTED"),
        ({**b2_tyk2_base, "tyk2_pe4max_pentaplex_msi_shifts_k": 1, "obs_tyk2_pe4max_pentaplex_msi_shifts_k": 0}, "FAIL_MSI_DETECTED"),
    ):
        for fn in (
            lambda d: wc_mod.evaluate_stage0_product_contract(d, root=ROOT),
            lambda d: ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(d), root=ROOT),
            lambda d: ing_mod.replay_stage0_chronological_batches([d], root=ROOT),
        ):
            assert exp_reason in fn(doc_a)["iip0026_reasons"], (exp_reason, doc_a)

    # 5. Atomic ledger wrong-campaign S9 isolation, 10**400 atom isolation, QP reset hardening, and MappingProxyType immutability:
    lot_k = wc_mod.LotScopeKey(bsct_scope["product_variant_id"], bsct_scope["construct_sha256"], bsct_scope["donor_seedline_id"], bsct_scope["mcb_id"], bsct_scope["lot_id"])
    mcb_k = wc_mod.MCBScopeKey(bsct_scope["product_variant_id"], bsct_scope["construct_sha256"], bsct_scope["donor_seedline_id"], bsct_scope["mcb_id"])
    a_bad_camp = wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.1, se_or_n=0.1, sample_n=6, campaign_id="OTHER", sample_unit_id="S1")
    a_ok = wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=5.0, se_or_n=0.1, sample_n=6, campaign_id="C_QA", sample_unit_id="S2")
    l_wc = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"atoms": [a_bad_camp]}, {"atoms": [a_ok]}], bsct_scope, target_campaign_id="C_QA", root=ROOT)
    assert l_wc["hard_stop_precedence_winner"] is None and any("HOLD_UNLINKED_CAMPAIGN" in h for h in l_wc["holds"]), l_wc
    a_s9 = wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=0.1, se_or_n=0.1, sample_n=6, campaign_id="C_QA", sample_unit_id="S3")
    a_huge = wc_mod.AtomicObservation(endpoint_id="mcb_1q32_mdm4_diploid_cn", scope_type="MCB_SCOPE", scope_key=mcb_k, kind="CONTINUOUS_STUDENT_T", value_or_k=2.0, se_or_n=0.01, sample_n=10**400, campaign_id="C_QA", sample_unit_id="S4")
    l_huge = wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"atoms": [a_s9]}, {"atoms": [a_huge]}], bsct_scope, target_campaign_id="C_QA", root=ROOT)
    assert l_huge["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP", l_huge
    # MappingProxyType deep-freeze cannot be mutated post-construction:
    mut_d = {"k_val": 4000}
    a_mp = wc_mod.AtomicObservation(endpoint_id="obs_icasp9_log10_kill", scope_type="LOT_SCOPE", scope_key=lot_k, kind="CONTINUOUS_STUDENT_T", value_or_k=5.0, se_or_n=0.1, sample_n=6, campaign_id="C_QA", sample_unit_id="S5", controls={"nested": MappingProxyType(mut_d)})
    d_pre = a_mp.digest()
    mut_d["k_val"] = 1
    assert a_mp.digest() == d_pre and a_mp.controls["nested"]["k_val"] == 4000

    # 6. Rule S11 fail-closed unknown ID, public product_def override resistance, CIITA 10**400 guard, paired ADCC exact decimal check, & 5'-RACE k>n check:
    assert wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate({"product_variant_id": "ISLET1_UNKNOWN_999"}, root=ROOT)["gate_pass"] is False
    fake_pdef = copy.deepcopy(pc_mod.get_product_definition("ISLET1_A4_ARM1B", root=ROOT))
    s11_ovr_test = wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate({"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "has_intronic_mir_cd276_pvr": True}, root=ROOT, product_def=fake_pdef)
    assert s11_ovr_test["status"] == "FAIL_S11_MIR_CD276_PVR_BARRED_ON_DESIGN_B_HLAC_INTACT"
    s11_adcc_bad = wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate(
        {**l2_scope, "iip0031_locus2_record": {"obs_l2_v3_intact_serum_adcc_lysis_pct": 2.0, "obs_l2_v3_igg_depleted_serum_adcc_lysis_pct": 1.0, "obs_l2_v3_igg_increment_adcc_lysis_pct": 9.0, "obs_l2_v3_igg_increment_adcc_lysis_se": 0.1, "replicates_n": 6}},
        root=ROOT,
    )
    assert s11_adcc_bad["status"].startswith("HOLD_S11_") and "FAIL_S11_CONJUNCTIVE_GY023_ADCC_EXCEEDS_5PCT" not in s11_adcc_bad["reasons"], s11_adcc_bad
    s11_ciita_huge = wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate(
        {**l2_scope, "iip0031_locus2_record": {"replicates_n": 6, "obs_ciita_ki_allele_informative_reads_n": 10**400, "obs_ciita_ki_ex4_ex6_aberrant_skip_reads_k": 0, "obs_ciita_ko_allele_informative_reads_n": 500, "obs_ciita_ko_ex4_ex6_aberrant_skip_reads_k": 20}},
        root=ROOT,
    )
    assert s11_ciita_huge["status"] == "FAIL_S11_CIITA_ABERRANT_EX4_EX6_SKIP_EXCEEDS_1PCT", s11_ciita_huge
    s11_race_oob = wc_mod.evaluate_locus2_ciita_modeb_challenger_v3_gate(
        {**l2_scope, "iip0031_locus2_record": {"replicates_n": 6, "obs_l2_v3_unit1_eef1a1_resting_0h_5race_umi_n_eff": 4000, "obs_l2_v3_unit1_eef1a1_resting_0h_unexpected_tss_umi_k": 4001}},
        root=ROOT,
    )
    assert s11_race_oob["status"].startswith("HOLD_S11_") and "FAIL_S11_LOCUS2_UNEXPECTED_5PRIME_RACE_START_DETECTED" not in s11_race_oob["reasons"], s11_race_oob

    # 7. Math-only parameter validation & exact Layer-3 2WII DBREF Gln1013 (33.04 A) / 15-PDB inventory:
    for bad_cap in (-1, float("inf"), float("nan"), False, "0.01"):
        with pytest.raises(ValueError):
            wc_mod.compute_tier6b_per_fraction_thinning_bound_math_only(
                [{"fraction_id": "a", "N": 10_000_000, "w_i": 0.5, "d_i_L": 0.5}, {"fraction_id": "b", "N": 10_000_000, "w_i": 0.5, "d_i_L": 0.5}],
                lot_cap=bad_cap,
            )
    l3_atlas = json.loads((ROOT / "data/extracted/layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    assert l3_atlas["summary"]["total_experimental_pdb_structures_verified"] == 15
    assert l3_atlas["locus2_cfh_stp_gpi_structural_provenance_iip_0031"]["pdb_2wii"]["glu264_ca_to_c3b_ted_q1013_ca_angstrom"] == 33.04

    # 8. Patch v11 table-driven 4-path mutation conformance:
    # (a) Conflicting flat Rule S11 measurements across supported containers MUST NEVER erase an attributable FAIL_S11_* hard stop
    # across Direct, Dict Ingestion, Wire Ingestion, Chronological Replay (single- and multi-record), and Atomic Ledger!
    l2_lot_k = wc_mod.LotScopeKey(
        l2_scope["product_variant_id"],
        l2_scope["construct_sha256"],
        l2_scope["donor_seedline_id"],
        l2_scope["mcb_id"],
        l2_scope["lot_id"],
    )
    flat_s11_conflict_cases = [
        (
            "lectin_18_vs_4",
            {"obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.1, "replicates_n": 6},
            {"obs_l2_v3_lectin_c4b2a_lysis_pct": 4.0},
            "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT",
        ),
        (
            "lectin_18_vs_null",
            {"obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.1, "replicates_n": 6},
            {"obs_l2_v3_lectin_c4b2a_lysis_pct": None},
            "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT",
        ),
        (
            "classical_18_vs_4",
            {"obs_l2_v3_classical_c4b2a_lysis_pct": 18.0, "obs_l2_v3_classical_c4b2a_lysis_se": 0.1, "replicates_n": 6},
            {"obs_l2_v3_classical_c4b2a_lysis_pct": 4.0, "obs_l2_v3_classical_c4b2a_lysis_se": 0.2},
            "FAIL_S11_POST_CUFF_CLASSICAL_C4B2A_LYSIS_EXCEEDS_10PCT",
        ),
        (
            "ap_c5b9_8_vs_1",
            {"obs_l2_v3_ap_c5b9_pct": 8.0, "obs_l2_v3_ap_c5b9_se": 0.1, "replicates_n": 6},
            {"obs_l2_v3_ap_c5b9_pct": 1.0},
            "FAIL_S11_POST_CUFF_AP_C5B9_EXCEEDS_5PCT",
        ),
        (
            "intact_adcc_12_vs_2",
            {"obs_l2_v3_intact_serum_adcc_lysis_pct": 12.0, "obs_l2_v3_intact_serum_adcc_lysis_se": 0.1, "replicates_n": 6},
            {"obs_l2_v3_intact_serum_adcc_lysis_pct": 2.0},
            "FAIL_S11_CONJUNCTIVE_GY023_ADCC_EXCEEDS_5PCT",
        ),
        (
            "lectin_with_conflicting_bad_replicates_n",
            {"obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.1, "replicates_n": 6},
            {"obs_l2_v3_lectin_c4b2a_lysis_pct": 4.0, "replicates_n": 2},
            "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT",
        ),
        (
            "chr16_ha_wgs_false_vs_true",
            {"obs_l2_v3_chr16_ha_wgs_pass": False, "obs_l2_v3_chr16_ha_wgs_variant_count": 2},
            {"obs_l2_v3_chr16_ha_wgs_pass": True, "obs_l2_v3_chr16_ha_wgs_variant_count": 0},
            "FAIL_S11_CHR16_HA_WGS_VARIANT_DETECTED",
        ),
    ]
    for case_lbl, adv_block, safe_block, exp_s11 in flat_s11_conflict_cases:
        # Test BOTH container orientations (adverse in measurements vs safe in mol, AND adverse in mol vs safe in measurements):
        for meas_blk, mol_blk in ((adv_block, safe_block), (safe_block, adv_block)):
            pkt = {
                "scope": dict(l2_scope),
                "measurements": dict(meas_blk),
                "iip0026_molecular_observations": dict(mol_blk),
            }
            for p_nm, res_p in (
                ("direct", wc_mod.evaluate_stage0_product_contract(pkt, root=ROOT)),
                ("ingest_dict", ing_mod.ingest_stage0_research_or_wetlab_batch(pkt, root=ROOT)),
                ("ingest_wire", ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(pkt), root=ROOT)),
                ("replay_1rec", ing_mod.replay_stage0_chronological_batches([pkt], root=ROOT)),
            ):
                assert res_p["overall_status"] == exp_s11, (case_lbl, p_nm, res_p)
                assert res_p["hard_stop_precedence_winner"] == exp_s11, (case_lbl, p_nm, res_p)
                assert res_p["s11_evaluation"]["status"] == exp_s11, (case_lbl, p_nm, res_p)
                assert "HOLD_CONFLICTING_STAGE0_MEASUREMENT" in res_p["iip0026_reasons"], (case_lbl, p_nm, res_p)
        # Atomic ledger parity for contradictory structured panels:
        a_adv = wc_mod.AtomicObservation(
            endpoint_id="stage0_locus2_modeb_v3_panel",
            scope_type="LOT_SCOPE",
            scope_key=l2_lot_k,
            kind="STRUCTURED_PANEL",
            value_or_k={**adv_block, "scope_type": "LOT_SCOPE"},
            campaign_id="C_QA",
            sample_unit_id="S1",
        )
        a_safe = wc_mod.AtomicObservation(
            endpoint_id="stage0_locus2_modeb_v3_panel",
            scope_type="LOT_SCOPE",
            scope_key=l2_lot_k,
            kind="STRUCTURED_PANEL",
            value_or_k={**safe_block, "scope_type": "LOT_SCOPE"},
            campaign_id="C_QA",
            sample_unit_id="S2",
        )
        for batch_seq in ([{"atoms": [a_adv]}, {"atoms": [a_safe]}], [{"atoms": [a_safe]}, {"atoms": [a_adv]}]):
            res_l = wc_mod.evaluate_iip0024_atomic_stage0_ledger(batch_seq, l2_scope, target_campaign_id="C_QA", root=ROOT)
            assert res_l["overall_status"] == exp_s11, (case_lbl, res_l)
            assert res_l["hard_stop_precedence_winner"] == exp_s11, (case_lbl, res_l)
            assert res_l["s11_evaluation"]["status"] == exp_s11, (case_lbl, res_l)

    # Multi-record chronological replay where SE & replicates_n arrive in Record 1 and conflicting flat lectin means arrive in Record 2:
    r_multi_1 = {"scope": dict(l2_scope), "measurements": {"obs_l2_v3_lectin_c4b2a_lysis_se": 0.1, "replicates_n": 6}, "utc_timestamp": "2026-11-20T01:00:00Z", "seq_index": 1}
    r_multi_2 = {"scope": dict(l2_scope), "measurements": {"obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0}, "iip0026_molecular_observations": {"obs_l2_v3_lectin_c4b2a_lysis_pct": 4.0}, "utc_timestamp": "2026-11-20T01:01:00Z", "seq_index": 2}
    res_rm = ing_mod.replay_stage0_chronological_batches([r_multi_1, r_multi_2], root=ROOT)
    assert res_rm["overall_status"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", res_rm
    assert res_rm["hard_stop_precedence_winner"] == "FAIL_S11_POST_CUFF_LECTIN_C4B2A_LYSIS_EXCEEDS_10PCT", res_rm

    # Anti-cross-breeding control: safe mean+SE in measurements (4.0, SE=0.1) + bare adverse mean (18.0, NO SE) in mol MUST NOT fabricate SE=0.1 onto 18.0:
    pkt_no_xb = {
        "scope": dict(l2_scope),
        "measurements": {"obs_l2_v3_lectin_c4b2a_lysis_pct": 4.0, "obs_l2_v3_lectin_c4b2a_lysis_se": 0.1, "replicates_n": 6},
        "iip0026_molecular_observations": {"obs_l2_v3_lectin_c4b2a_lysis_pct": 18.0},
    }
    res_no_xb = wc_mod.evaluate_stage0_product_contract(pkt_no_xb, root=ROOT)
    assert res_no_xb["overall_status"] == "HOLD_CONFLICTING_STAGE0_MEASUREMENT", res_no_xb
    assert res_no_xb["s11_evaluation"]["status"] == "HOLD_S11_CONFLICTING_NESTED_EVIDENCE", res_no_xb

    # (b) Tier 6B strict fraction-role classification & optional typed `fraction_role` across ALL 4 paths:
    def _t6b_pair_v11(enr_lbl: str, ft_lbl: str = "flow_through_and_excluded", *, enr_extra: dict | None = None, ft_extra: dict | None = None) -> list[dict]:
        roster = [enr_lbl, ft_lbl]
        base_meta = {
            "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape",
            "scope_type": "LOT_SCOPE",
            "kind": "EXACT_INTEGER_COUNT",
            "product_variant_id": l2_scope["product_variant_id"],
            "construct_sha256": l2_scope["construct_sha256"],
            "background_id": l2_pdef["staged_locus2_modeb_donor_binding"]["background_id"],
            "donor_seedline_id": l2_scope["donor_seedline_id"],
            "mcb_id": l2_scope["mcb_id"],
            "lot_id": l2_scope["lot_id"],
            "campaign_id": l2_scope["campaign_id"],
            "timepoint": "STAGE7_PRE_ASSEMBLY",
            "assay_version": "TIER6B_SORTED_LOT_V1",
            "unit": "clonogenic_opportunities",
            "presort_lineage_panel": "CD31_NEG_CD144_NEG_ENDOCRINE_PRE_ASSEMBLY",
            "clonogenic_risk_unit": "SORTED_SINGLE_CELL_EQUIVALENT_CLONOGENIC_OPPORTUNITY",
            "endogenous_denominator_verified": True,
            "endogenous_denominator_basis": "MATCHED_UNTREATED_CLONOGENIC_CONTROL",
            "presort_mass_balance_reconciled": True,
            "unreconciled_loss_fraction": 0.0,
            "presort_total_input_cells": 3_000_000_000,
            "presort_reconciled_loss_cells": 0,
            "excluded_or_lost_material_discarded_not_in_graft": True,
            "declared_fraction_roster": roster,
            "k": 0,
        }
        r_enr = {
            **base_meta,
            "fraction_id": enr_lbl,
            "sample_id": "V11_ENR_S1",
            "N": 30_000_000,
            "plated_treated_cells_n": 80_000_000,
            "presort_fraction_cells": 2_999_700_000,
            "w_i": 0.9999,
            "hpsc_calibrant": {"r_k": 20, "r_n": 20, "cal_batch_id": "V11_CAL_E_H"},
            "edited_differentiated_calibrant": {"r_k": 20, "r_n": 20, "cal_batch_id": "V11_CAL_E_D"},
            "untreated_cfu_control": {"untreated_sample_id": "V11_UNTR_E", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            **(enr_extra or {}),
        }
        r_ft = {
            **base_meta,
            "fraction_id": ft_lbl,
            "sample_id": "V11_FT_S1",
            "N": 100_000,
            "plated_treated_cells_n": 266_667,
            "presort_fraction_cells": 300_000,
            "w_i": 0.0001,
            "hpsc_calibrant": {"r_k": 10, "r_n": 10, "cal_batch_id": "V11_CAL_F_H"},
            "edited_differentiated_calibrant": {"r_k": 10, "r_n": 10, "cal_batch_id": "V11_CAL_F_D"},
            "untreated_cfu_control": {"untreated_sample_id": "V11_UNTR_F", "untreated_cfu_colonies_k": 4000, "untreated_plated_cells_n": 10000},
            **(ft_extra or {}),
        }
        return [r_enr, r_ft]

    def _eval_t6b_4paths_v11(recs: list[dict]) -> dict[str, dict]:
        pkt = {"scope": dict(l2_scope), "tier6b_sorted_lot_observations": recs}
        atoms = [
            wc_mod.AtomicObservation(
                endpoint_id=r["endpoint_id"],
                scope_type="LOT_SCOPE",
                scope_key=l2_lot_k,
                kind="EXACT_INTEGER_COUNT",
                value_or_k=r["k"],
                se_or_n=r["N"],
                campaign_id=r["campaign_id"],
                sample_unit_id=r["sample_id"],
                unit=r["unit"],
                assay_version=r["assay_version"],
                timepoint=r["timepoint"],
                controls={k: v for k, v in r.items() if k not in ("endpoint_id", "scope_type", "kind", "k", "N", "product_variant_id", "construct_sha256", "donor_seedline_id", "mcb_id", "lot_id", "campaign_id", "sample_id", "unit", "assay_version", "timepoint")},
            )
            for r in recs
        ]
        return {
            "direct": wc_mod.evaluate_stage0_product_contract(pkt, root=ROOT)["tier6b_sorted_lot_evaluation"],
            "ingest_wire": ing_mod.ingest_stage0_research_or_wetlab_batch(json.dumps(pkt), root=ROOT)["tier6b_sorted_lot_evaluation"],
            "replay": ing_mod.replay_stage0_chronological_batches([pkt], root=ROOT)["tier6b_sorted_lot_evaluation"],
            "ledger": wc_mod.evaluate_iip0024_atomic_stage0_ledger([{"atoms": atoms}], l2_scope, target_campaign_id="C_QA", root=ROOT)["tier6b_sorted_lot_evaluation"],
        }

    # Positive controls (with and without explicit typed fraction_role) MUST PASS with exact 2.1870701159e-7 on ALL 4 paths:
    for pos_recs in (
        _t6b_pair_v11("enriched_cd31_neg_cd144_neg", "flow_through_and_excluded"),
        _t6b_pair_v11("enriched_cd31_neg_cd144_neg", "flow_through_and_excluded", enr_extra={"fraction_role": "ENRICHED"}, ft_extra={"fraction_role": "FLOW_THROUGH_OR_EXCLUDED"}),
    ):
        for p_nm, ev_pos in _eval_t6b_4paths_v11(pos_recs).items():
            assert ev_pos["disposition"] == "PASS", (p_nm, ev_pos)
            assert abs(ev_pos["p_lot_U95"] - 2.187070115946855e-7) < 1e-12, (p_nm, ev_pos)

    # Ambiguous/contradictory enriched labels (including Box 2's exact `depleted_enriched_cd31_neg_cd144_neg`, even with explicit `fraction_role="ENRICHED"`!) MUST HOLD on ALL 4 paths:
    for bad_enr_lbl, bad_enr_extra in (
        ("depleted_enriched_cd31_neg_cd144_neg", None),
        ("depleted_enriched_cd31_neg_cd144_neg", {"fraction_role": "ENRICHED"}),
        ("enriched_depleted_cd31_neg_cd144_neg", None),
        ("excluded_enriched_cd31_neg_cd144_neg", None),
        ("discard_enriched", None),
        ("waste_enriched", None),
        ("rejected_enriched", None),
        ("residual_enriched", None),
        ("unselected_enriched", None),
        ("mixed_enriched", None),
        ("unsorted_enriched", None),
        ("bulk_enriched", None),
        ("pseudo_enriched_cd31_neg_cd144_neg", None),
        ("enriched_cd31_pos_cd144_neg", None),
        ("enriched_cd31_neg_cd144_neg", {"fraction_role": "FLOW_THROUGH_OR_EXCLUDED"}),
        ("enriched_cd31_neg_cd144_neg", {"fraction_role": "UNKNOWN_ROLE"}),
    ):
        bad_recs = _t6b_pair_v11(bad_enr_lbl, "flow_through_and_excluded", enr_extra=bad_enr_extra)
        for p_nm, ev_bad in _eval_t6b_4paths_v11(bad_recs).items():
            assert ev_bad["disposition"] == "HOLD_TIER6B_UNQUALIFIED_FRACTION", (bad_enr_lbl, bad_enr_extra, p_nm, ev_bad)
            assert ev_bad["gate_pass"] is False, (bad_enr_lbl, bad_enr_extra, p_nm, ev_bad)



