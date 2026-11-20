"""Root policy negative controls. All assay inputs are synthetic test fixtures.

Only pure normalization/calibration APIs run here; canonical assay artifacts are
never written. Explicit require checks remain active with Python optimization.
"""
from __future__ import annotations
import copy
import hashlib
import json
import math
import sys
import subprocess
import time
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import engineering_priors as engineering
from scripts import generate_release_manifest as exporter
from scripts import ingest_wetlab_results as ingest
from scripts import product_contract as products
from scripts import public_excerpt_policy as excerpts
from scripts import wetlab_contract as measurements


def require(condition,message):
    if not condition:raise AssertionError(message)


def token(quote="Synthetic source fixture for hashing only.",*,pmid="12345678",metadata="field=verbatim_evidence_quote genes=B2M,CD47 pert=KO_LOF"):
    return f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} {metadata} quote_sha256={hashlib.sha256(quote.encode()).hexdigest()}]"


def test_opaque_token_requires_actual_source_digest_and_enclosing_pmid():
    quote="Synthetic source fixture for hashing only.";valid=token(quote)
    require(excerpts.validate_token(valid,expected_pmid="12345678",source_quote=quote),"genuine generated token failed")
    require(not excerpts.validate_token(valid,expected_pmid="87654321",source_quote=quote),"wrong enclosing PMID accepted")
    require(not excerpts.validate_token(valid,source_quote=quote+" changed"),"digest not source-bound")
    require(excerpts.validate_token(valid,expected_pmid="12345678"),"bounded opaque reference failed when source unavailable")
    require(exporter._validate_public_redaction_token(valid,expected_pmid="12345678",source_quote=quote),"exporter bypassed shared parser")


@pytest.mark.parametrize("metadata",[
    "field=verbatim_evidence_quote field=note","field=verbatim_evidence_quote unknown=INS",
    "field=note span=source,prose","field=source,prose,encoded",
    "field=note genes=INS,this,is,source,prose","field=note genes=B2M,B2M","field=note genes=NOT_A_GENE",
    "field=note mod=source,prose","field=note pert=source,prose","field=note pert=GOF_OE pert=KO_LOF",
    "field=note genes=","field=note=more","genes=INS pert=NONE",
])
def test_opaque_tokens_reject_duplicate_unknown_and_source_encoded_fields(metadata):
    require(not excerpts.validate_token(token(metadata=metadata),expected_pmid="12345678"),"untrusted metadata exemption")


@pytest.mark.parametrize("digest",["0"*64,"f"*63,"g"*64,"F"*64,"not-a-digest"])
def test_redaction_rejects_invalid_or_source_unmatched_hash(digest):
    valid=token();altered=valid[:valid.index("quote_sha256=")]+"quote_sha256="+digest+"]"
    require(not excerpts.validate_token(altered,source_quote="Synthetic source fixture for hashing only."),"fake hash accepted")


@pytest.mark.parametrize("source_path",["abstract","fulltext_data_sections.results_paragraphs[0]"])
def test_import_token_resolves_exact_source_path_and_rejects_mutation(source_path):
    source="prefix synthetic quote suffix";start,end=7,22
    dossier={"abstract":source,"fulltext_data_sections":{"results_paragraphs":[source]}}
    record={"source_span":f"{source_path}[{start}:{end}]|norm=exact","verbatim_evidence_quote":token(source[start:end])}
    exporter._validate_imported_redaction_tokens(record,"12345678",dossier)
    wrong=copy.deepcopy(record);wrong["verbatim_evidence_quote"]=token("another source")
    with pytest.raises(ValueError):exporter._validate_imported_redaction_tokens(wrong,"12345678",dossier)
    wrong=copy.deepcopy(record);wrong["source_span"]="unrelated.path[0:1]|norm=exact"
    with pytest.raises(ValueError):exporter._validate_imported_redaction_tokens(wrong,"12345678",dossier)


def test_unknown_fields_and_trailing_prose_are_not_exemptions():
    require(excerpts.safe_field("verbatim_quote")=="verbatim_quote","valid bounded field missing")
    require(excerpts.safe_field("this is source-bearing prose")=="omitted_text","free-text exemption")
    # Bound the negative control externally: a parser regression must fail CI,
    # never hang the entire verification run through catastrophic backtracking.
    code="import json,sys;from scripts.public_excerpt_policy import validate_token;sys.exit(1 if validate_token(json.load(sys.stdin)) else 0)"
    started=time.monotonic()
    proc=subprocess.run([sys.executable,"-B","-c",code],input=json.dumps(token()+" appended source prose"),text=True,capture_output=True,cwd=ROOT,timeout=2)
    require(time.monotonic()-started<2.,"bounded token rejection exceeded the 2-second whole-process deadline")
    require(proc.returncode==0,"trailing prose accepted by the opaque-token validator")
    require(not excerpts.validate_token(token(metadata="field=note genes="+",".join(["INS"]*33))),"unbounded repeated genes accepted")


@pytest.fixture(scope="module")
def product_rows():return products.load_product_contracts()["products"]


def test_exact_products_separate_six_locus_comparator_and_five_locus_bank(product_rows):
    by_id={row["product_variant_id"]:row for row in product_rows}
    for identifier,k,b in [("ISLET1_B1_6LOCUS_TYK2_KO_INS",6,1),("ISLET1_B2_5LOCUS_BASE",5,0),("ISLET1_B2_5LOCUS_TYK2_PE",5,1)]:
        require((by_id[identifier]["K"],by_id[identifier]["B"])==(k,b),"distinct B variant silently substituted")
    for row in product_rows:
        require(row["K"]==len(set(row["dsb_or_indel_loci"])) and row["B"]==len(set(row["base_or_prime_edit_loci"])),"K/B and explicit loci diverge")
        require(row["physical_qualification_status"]=="PREREGISTERED_PENDING" and row["physical_manufacturing_yield"] is None,"unmeasured product/yield qualified")
        require(row["donor_lineage"] and row["hla_retention"] and row["pending_assay_gates"],"donor/HLA/assay scope omitted")
        if row["family"]=="B":
            require(row["execution_supported"] is False,"unmaterialized B design executable")
            require("retained" in row["alleles"]["HLA-C"],"B bank HLA-C retention omitted")


@pytest.mark.parametrize("key,value",[
    ("K",4.),("K",True),("K",-1),("K",5),("B",0.),("B",True),
    ("execution_supported","true"),("family","UNKNOWN"),("product_variant_id","source-bearing prose"),
    ("dsb_or_indel_loci",["B2M","B2M","CXCL10","AAVS1"]),("dsb_or_indel_loci",["B2M","CIITA","CXCL10","not a locus"]),
    ("alleles",{}),("cassette_topology",{}),("donor_lineage",{}),("hla_retention",{}),("pending_assay_gates",[]),
    ("aliases",["source prose"]),("physical_qualification_status","VALIDATED"),("physical_manufacturing_yield",.99),
    ("sequence_materialization_status","NOT_MATERIALIZED_SEPARATE_DONOR_LEDGER_REQUIRED"),
])
def test_product_rejects_incomplete_or_mistyped_variants(product_rows,key,value):
    row=copy.deepcopy(product_rows[0]);row[key]=value
    with pytest.raises(ValueError):products.validate_product_definition(row)


def test_product_requires_edited_alleles_and_pending_scope(product_rows):
    for key in ("alleles","donor_lineage","hla_retention","pending_assay_gates"):
        row=copy.deepcopy(product_rows[0]);del row[key]
        with pytest.raises(ValueError):products.validate_product_definition(row)
    row=copy.deepcopy(product_rows[0]);del row["alleles"]["B2M"]
    with pytest.raises(ValueError):products.validate_product_definition(row)
    row=copy.deepcopy(next(r for r in product_rows if r["family"]=="B"));row["execution_supported"]=True
    with pytest.raises(ValueError):products.validate_product_definition(row)


@pytest.mark.parametrize("collision",["id","alias_alias","alias_id"])
def test_registry_rejects_global_identifier_alias_collisions(tmp_path,product_rows,collision):
    rows=copy.deepcopy(product_rows)
    if collision=="id":rows[1]["product_variant_id"]=rows[0]["product_variant_id"]
    elif collision=="alias_alias":rows[0]["aliases"]=["SHARED_ALIAS"];rows[1]["aliases"]=["SHARED_ALIAS"]
    else:rows[0]["aliases"]=[rows[1]["product_variant_id"]]
    path=tmp_path/products.CONTRACT_PATH;path.parent.mkdir(parents=True);path.write_text(json.dumps({"schema_version":"0.8.0","products":rows}))
    with pytest.raises(ValueError):products.load_product_contracts(tmp_path)


def test_product_lookup_is_copy_isolated_and_unknown_ids_fail(product_rows):
    identifier=product_rows[0]["product_variant_id"];row=products.get_product_definition(identifier);row["alleles"]["B2M"]="MUTATED"
    require(products.get_product_definition(identifier)["alleles"]["B2M"]!="MUTATED","lookup poisoned contract")
    with pytest.raises(ValueError):products.get_product_definition("ISLET1_UNKNOWN_VARIANT")


@pytest.mark.parametrize("K,B",[(k,b)for k in range(4,7)for b in range(2)])
def test_engineering_utility_independent_formula(K,B):
    for lam,gamma in ((.15,.965),(.08,.92),(.25,.985),(0.,1.)):
        got=engineering.engineering_complexity_discount(K,B,lam,gamma)
        require(abs(got-math.exp(-lam*max(0,K-4))*gamma**B)<1e-14,"shared utility formula drift")
        require(0<got<=1,"utility outside dimensionless domain")


@pytest.mark.parametrize("args",[(-1,0),(4,-1),(4.,0),(4,0.),(True,0),(4,False),(4,0,-.1,.965),(4,0,float('nan'),.965),(4,0,float('inf'),.965),(4,0,.15,0),(4,0,.15,1.1),(4,0,.15,float('nan')),(4,0,.15,True)])
def test_engineering_utility_rejects_bad_domains(args):
    with pytest.raises(ValueError):engineering.engineering_complexity_discount(*args)


CONTRACT=measurements.derive_contract();ENDPOINTS=CONTRACT['endpoints'];VALUES=[.95,6.2,6.,4.4,25.,.02,7.42,99.5]


def payload_for(index,**overrides):
    e=ENDPOINTS[index];p={'batch_id':f'QA-ROOT-ENDPOINT-{index}',e['meanKey']:VALUES[index],e['seKey']:max(.001,e['default_se']or.03)};p.update(overrides);return p


def test_generated_endpoint_contract_matches_current_constants():
    require(len(ENDPOINTS)==8 and len({e['meanKey']for e in ENDPOINTS})==8,'8-endpoint parity lost')
    require([e['se_required']for e in ENDPOINTS]==[True]*4+[False]*4,'SE policy drift')
    require([e['default_se']for e in ENDPOINTS[4:]]==[.85,.005,.020,.25],'extended default drift')
    require(CONTRACT['minimum_se']==.001,'SE floor drift')
    require(json.loads((ROOT/measurements.CONTRACT_RELATIVE_PATH).read_text())==CONTRACT,'generated endpoint metadata stale')


@pytest.mark.parametrize('index',range(8))
def test_all_eight_endpoints_normalize_and_reach_pure_ingestor(index,monkeypatch):
    e=ENDPOINTS[index];payload=payload_for(index);normalized,v=measurements.normalize_measurements(payload)
    require(v['validated_endpoints']==[e['meanKey']]and normalized[e['meanKey']]==VALUES[index],'measurement lost or changed')
    def no_write(*_a,**_k):raise AssertionError('pure calibration attempted disk write')
    monkeypatch.setattr(ingest,'_atomic_write_json',no_write);r=ingest.run_calibration(payload,False)
    require(r['observed_endpoints_updated_count']==1,'single observed endpoint lost')
    table=r['conjugate_posterior_updates']if index<4 else r['optional_pouch_and_36run_screen_updates_v80_070']
    require(table[e['endpoint_id']]['status']=='UPDATED_FROM_OBSERVED_READOUT','observation retained prior')


@pytest.mark.parametrize('index',range(8))
def test_unavailable_retains_missingness_and_orphan_se_rejects(index):
    e=ENDPOINTS[index];n,v=measurements.normalize_measurements({e['meanKey']:'UNAVAILABLE'},require_observation=False)
    require(n[e['meanKey']]=='UNAVAILABLE'and e['seKey']not in n and v['validated_endpoints']==[],'missingness fabricated a measurement')
    for payload in ({e['meanKey']:'UNAVAILABLE',e['seKey']:.03},{e['seKey']:.03}):
        with pytest.raises(ValueError):measurements.normalize_measurements(payload,require_observation=False)


@pytest.mark.parametrize('index',range(8))
def test_required_and_default_se_rules_match_ingestor(index):
    e=ENDPOINTS[index];payload=payload_for(index);del payload[e['seKey']]
    if index<4:
        with pytest.raises(ValueError):measurements.normalize_measurements(payload)
        with pytest.raises(ValueError):ingest.run_calibration(payload,False)
    else:
        n,v=measurements.normalize_measurements(payload)
        require(n[e['seKey']]==e['default_se']and v['defaulted_standard_error_keys']==[e['seKey']],'default SE absent or unqualified')
        r=ingest.run_calibration(payload,False)
        require(r['measurement_validation']['default_se_status']=='DECLARED_ASSUMPTION_NOT_MEASURED','default falsely measured')


@pytest.mark.parametrize('bad',['',' ','NaN','Infinity','0x10','1_0','1,0','١','１２',True,None,{},[],float('nan'),float('inf')])
def test_blank_non_decimal_and_nonfinite_reject_before_conversion(bad):
    with pytest.raises(ValueError):measurements.finite_decimal(bad,'probe')
    with pytest.raises(ValueError):measurements.normalize_measurements(payload_for(0,obs_stage7_yield_ratio=bad))
    with pytest.raises(ValueError):ingest.run_calibration(payload_for(0,obs_stage7_yield_ratio=bad),False)


@pytest.mark.parametrize('index',range(8))
def test_bounds_se_floor_and_unknown_fields_reject_in_both_paths(index):
    e=ENDPOINTS[index]
    for patch in ({e['meanKey']:e['min']-.01},{e['meanKey']:e['max']+.01},{e['seKey']:.000999},{e['seKey']:''},{'obs_unknown_endpoint':1}):
        payload=payload_for(index,**patch)
        with pytest.raises(ValueError):measurements.normalize_measurements(payload)
        with pytest.raises(ValueError):ingest.run_calibration(payload,False)
    n,_=measurements.normalize_measurements(payload_for(index,**{e['seKey']:'1e-3'}))
    require(n[e['seKey']]==.001,'valid scientific SE rejected')


@pytest.mark.parametrize('value,expected',[('+1.25e-2',.0125),('  1E+2  ',100.),('.5',.5),('1.',1.),('-0',0.)])
def test_decimal_scientific_positive_controls(value,expected):require(measurements.finite_decimal(value,'probe')==expected,'valid numeric string drift')


def test_metadata_unknown_keys_and_duplicate_json_keys_reject():
    for patch in ({'batch_id':''},{'batch_id':'  '},{'batch_id':False},{'assay_id':' '},{'unrecognized_metadata':'x'}):
        with pytest.raises(ValueError):measurements.normalize_measurements(payload_for(0,**patch))
    with pytest.raises(ValueError):ingest.loads_strict_json('{"obs_stage7_yield_ratio":0.9,"obs_stage7_yield_ratio":1}')
    with pytest.raises(ValueError):measurements.normalize_measurements({},require_observation=True)


@pytest.mark.parametrize("value_pair", ["0,1", "1,0"])
@pytest.mark.parametrize("container", [
    '{"k":%s,"k":%s}', '{"measurements":{"k":%s,"k":%s}}',
    '{"tier6b_sorted_lot_observations":[{"controls":{"k":%s,"k":%s}}]}',
])
def test_stage0_serialized_ingestion_rejects_duplicate_keys(container, value_pair):
    raw = container % tuple(value_pair.split(","))
    with pytest.raises(ValueError, match="Duplicate JSON key"):
        ingest.ingest_stage0_research_or_wetlab_batch(raw)


@pytest.mark.parametrize("token", ["NaN", "Infinity", "-Infinity", "1e9999", "-1e9999"])
def test_stage0_serialized_ingestion_rejects_nonfinite_numbers(token):
    raw = '{"measurements":{"value":' + token + '}}'
    with pytest.raises(ValueError, match="finite"):
        ingest.ingest_stage0_research_or_wetlab_batch(raw)


@pytest.mark.parametrize("value", [None, [], False, "", ["key", "value"]])
@pytest.mark.parametrize("name", ["scope", "measurements", "iip0026_molecular_observations"])
def test_stage0_ingestion_rejects_malformed_structured_blocks(value, name):
    packet = {name: value}
    for document in (packet, json.dumps(packet)):
        with pytest.raises(ValueError, match="must be an object"):
            ingest.ingest_stage0_research_or_wetlab_batch(document)


@pytest.mark.parametrize("value", [[], [["lot_id", "one"]], None, 1, True])
def test_stage0_ingestion_requires_object_root(value):
    for document in (value, json.dumps(value)):
        with pytest.raises(ValueError, match="must decode to a JSON object"):
            ingest.ingest_stage0_research_or_wetlab_batch(document)


@pytest.mark.parametrize("malformed", [
    {}, {"scope": []}, {"scope": {"lot_id": "orphan"}}, None,
    {"measurements": ["malformed"]}, {"iip0026_molecular_observations": False},
    {"measurements": {"obs_icasp9_log10_kill": True}},
    {"measurements": {"obs_icasp9_log10_kill": 10**600}},
])
@pytest.mark.parametrize("adverse_first", [True, False])
def test_stage0_replay_keeps_s9_stop_across_malformed_records(malformed, adverse_first):
    scope = {"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "lot_id": "LOT-TEST", "campaign_id": "CAMP-TEST"}
    adverse = {"scope": scope, "measurements": {"obs_icasp9_log10_kill": 0.1, "obs_icasp9_log10_kill_se": 0.1}}
    if malformed and "scope" not in malformed:
        malformed = {"scope": scope, **malformed}
    records = [adverse, malformed] if adverse_first else [malformed, adverse]
    out = ingest.replay_stage0_chronological_batches(records)
    assert out["overall_status"] == "HARD_SAFETY_STOP_REJECTED"
    assert out["hard_stop_precedence_winner"] == "S9_HARD_SAFETY_STOP"
    assert out["effective_execution_permitted"] is False


def test_stage0_replay_rejects_two_attributable_lots():
    scope = {"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "lot_id": "A", "campaign_id": "C"}
    records = [{"scope": scope}, {"scope": {**scope, "lot_id": "B"}}]
    for ordered in (records, records[::-1]):
        out = ingest.replay_stage0_chronological_batches(ordered)
        assert out["overall_status"] == "FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING"
        assert out["effective_execution_permitted"] is False


@pytest.mark.parametrize("key", ["lot_id", "current_lot_id", "campaign_id", "assay_campaign_id", "product_variant_id"])
def test_stage0_rejects_ambiguous_envelope_provenance(key):
    packet = {"scope": {"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "lot_id": "A", "campaign_id": "C"},
              key: "OTHER"}
    with pytest.raises(ValueError, match="Conflicting Stage-0 batch scope"):
        ingest.ingest_stage0_research_or_wetlab_batch(json.dumps(packet))
    out = ingest.replay_stage0_chronological_batches([packet])
    assert out["overall_status"] == "FAIL_CROSS_LOT_OR_CAMPAIGN_REPLAY_MIXING"
    assert out["effective_execution_permitted"] is False


@pytest.mark.parametrize("explicit_malformed", [False, True])
@pytest.mark.parametrize("adverse_first", [True, False])
def test_stage0_replay_keeps_attributable_outgrowth_with_malformed_peers(explicit_malformed, adverse_first):
    scope = {"product_variant_id": "ISLET1_B2_5LOCUS_BASE", "lot_id": "LOT-TEST", "campaign_id": "CAMP-TEST",
             "construct_sha256": hashlib.sha256(b"transport-test").hexdigest(),
             "donor_seedline_id": "TEST-SEED", "mcb_id": "TEST-MCB"}
    count = {**scope, "endpoint_id": "lot_tier6b_sorted_prolif_clonogenic_escape", "scope_type": "LOT_SCOPE",
             "sample_id": "transport-test", "k": 1}
    adverse = {"scope": scope, "tier6b_sorted_lot_observations": [count], "malformed_record": explicit_malformed}
    zero = {"scope": scope, "tier6b_sorted_lot_observations": [{**count, "k": 0}]}
    records = [adverse, {}, zero] if adverse_first else [zero, {}, adverse]
    out = ingest.replay_stage0_chronological_batches(records)
    assert out["overall_status"] == "FAIL_HARD_STOP"
    assert out["hard_stop_precedence_winner"] == "FAIL_TIER6B_OUTGROWTH_DETECTED"
    assert out["effective_execution_permitted"] is False


def test_report_integrity_lineage_and_replay_remain_fail_closed():
    report=ingest.run_calibration(payload_for(0),False);next_payload=payload_for(0,batch_id='QA-NEXT',obs_stage7_yield_ratio=.96)
    for lineage in ingest.SUPPORTED_EMPIRICAL_LINEAGES:
        legacy=copy.deepcopy(report);legacy['epistemic_lineage_status']=lineage
        # A changed label must not bypass the complete-state integrity gate.
        if lineage != report['epistemic_lineage_status']:
            with pytest.raises(ValueError,match='complete calibration-state digest mismatch'):
                ingest.run_calibration(next_payload,False,prior_posterior_report=legacy)
        # Explicitly sealed in-memory positive fixture; this authenticates no assay.
        legacy['calibration_state_sha256']=ingest._calibration_state_sha256(legacy)
        require(ingest.run_calibration(next_payload,False,prior_posterior_report=legacy)['observed_endpoints_updated_count']==1,'valid legacy lineage rejected')
    for lineage in (None,'SYNTHETIC_DRY_RUN_LINEAGE','MIXED_LINEAGE'):
        invalid=copy.deepcopy(report);invalid['epistemic_lineage_status']=lineage
        invalid['calibration_state_sha256']=ingest._calibration_state_sha256(invalid)
        with pytest.raises(ValueError):ingest.run_calibration(next_payload,False,prior_posterior_report=invalid)
    missing=copy.deepcopy(report);del missing['conjugate_posterior_updates_sha256']
    missing['calibration_state_sha256']=ingest._calibration_state_sha256(missing)
    with pytest.raises(ValueError):ingest.run_calibration(next_payload,False,prior_posterior_report=missing)
    bad=copy.deepcopy(report);next(iter(bad['conjugate_posterior_updates'].values()))['posterior_mean_exact']=.123
    bad['calibration_state_sha256']=ingest._calibration_state_sha256(bad)
    with pytest.raises(ValueError):ingest.run_calibration(next_payload,False,prior_posterior_report=bad)
    with pytest.raises(ValueError):ingest.run_calibration(payload_for(0,batch_id='QA-RENAMED-SAME-CONTENT'),False,prior_posterior_report=report)
    dry=ingest.run_calibration(payload_for(0,batch_id='DRY-RUN-QA'),True)
    with pytest.raises(ValueError):ingest.run_calibration(next_payload,False,prior_posterior_report=dry)


def test_generated_complement_gene_tokens_are_bounded_positive_controls():
    quote = "Complement reference fixture"
    digest = hashlib.sha256(quote.encode()).hexdigest()
    marker = f"[REDACTED_NON_OA_EXCERPT: see PMID:41385368 field=verbatim_evidence_quote genes=CD55,CD59,THBD pert=NONE quote_sha256={digest}]"
    require(excerpts.validate_token(marker, expected_pmid="41385368", source_quote=quote), "known complement symbols must roundtrip")


@pytest.mark.parametrize("value,field", [
    ("Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture", "max_biological_fidelity_tier"),
    ("Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture", "Z8_species_age_and_sex_context"),
    ("fulltext_data_sections.figure_legends[59].caption[411:588]|norm=exact", "source_span"),
])
def test_public_reexport_preserves_closed_schema_codes_and_typed_offsets(value,field):
    require(exporter._is_fixed_public_schema_metadata(value,field), "fixed schema annotation must survive source-free re-export")


@pytest.mark.parametrize("value,field", [
    ("Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture copied source words", "max_biological_fidelity_tier"),
    ("fulltext_data_sections.encoded_source_words_here[59][411:588]|norm=exact", "source_span"),
    ("any six word source sentence here", "source_span"),
])
def test_public_metadata_exemption_cannot_accept_prose_or_unbounded_source_paths(value,field):
    require(not exporter._is_fixed_public_schema_metadata(value,field), "source prose is not fixed schema metadata")


def test_opaque_arm_reference_rewrite_is_bijective_and_preserves_foreign_keys():
    old = "pmid41385368_arm1_some_descriptive_source_words_here"
    new = "pmid41385368_arm1"
    value = {"arms":[{"arm_id":old}], "effects":[{"focal_arm_id":old,"comparator_arm_id":None}], "unrelated":"retained"}
    rewritten = exporter._rewrite_exact_arm_references(value,{old:new})
    require(rewritten["arms"][0]["arm_id"]==rewritten["effects"][0]["focal_arm_id"]==new,"arm foreign-key identity lost")
    require(rewritten["effects"][0]["comparator_arm_id"] is None and rewritten["unrelated"]=="retained","unrelated metadata changed")
    require(exporter._rewrite_exact_arm_references(rewritten,{new:new})==rewritten,"opaque IDs are not idempotent")


def test_public_composite_gene_code_exemption_is_bounded():
    require(exporter._is_fixed_public_schema_metadata("HLA-B_HLA-C_HLA-A24","gene_symbol_hgnc"),"known legacy gene tuple lost")
    require(not exporter._is_fixed_public_schema_metadata("some_arbitrary_source_words_encoded_here","gene_symbol_hgnc"),"arbitrary source prose treated as a gene tuple")
