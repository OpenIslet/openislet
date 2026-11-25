from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_module(rel_path: str, module_name: str):
    path = REPO_ROOT / rel_path
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_context_evidence_graph_integrity_and_replay() -> None:
    ceg_mod = _load_module("scripts/build_context_evidence_graph.py", "openislet_ceg_builder")
    manifest_mod = _load_module("scripts/generate_release_manifest.py", "openislet_release_manifest")

    ceg_path = REPO_ROOT / ceg_mod.CEG_REL_PATH
    atlas_path = REPO_ROOT / ceg_mod.ATLAS_REL_PATH
    assert ceg_path.is_file(), f"Missing {ceg_path}"
    assert atlas_path.is_file(), f"Missing {atlas_path}"

    # Static asset ceiling (< 25 MiB) for direct web hosting and public bundle downloads
    assert ceg_path.stat().st_size < 25 * 1024 * 1024
    assert atlas_path.stat().st_size < 2 * 1024 * 1024

    ceg = json.loads(ceg_path.read_text(encoding="utf-8"))
    atlas_bytes = atlas_path.read_bytes()
    atlas = json.loads(atlas_bytes.decode("utf-8"))

    counts = ceg_mod.verify_ceg_payload(ceg, repo_root=REPO_ROOT)
    assert counts == {
        "studies": 1176,
        "arms": 7740,
        "outcomes": 14167,
        "edges": 2483,
        "effectors": 20,
        "graveyard_index": 10,
    }

    funnel = ceg["screening_funnel"]
    assert funnel["screened_records_total"] == 3938
    assert funnel["ring1_existing_records_screened"] + funnel["new_candidate_records_screened"] == 3938
    assert (
        funnel["candidates_advanced_total"]
        + funnel["contextual_reviews_parked_total"]
        + funnel["off_topic_records_excluded_total"]
        == 3938
    )
    assert funnel["candidates_advanced_total"] == 3263
    assert funnel["wave1_ceg_records_extracted_total"] == 1176
    assert funnel["wave1_primary_empirical_extracted_total"] == 1150
    assert funnel["wave1_trial_registry_protocols_extracted"] == 11
    assert funnel["wave1_regulatory_or_guidance_records_extracted"] == 7
    assert funnel["wave1_computational_models_extracted"] == 4
    assert funnel["wave1_correspondence_or_secondary_extracted"] == 1
    assert funnel["wave1_superseded_preprint_duplicates_extracted"] == 3
    assert (
        funnel["wave1_primary_empirical_extracted_total"]
        + funnel["wave1_trial_registry_protocols_extracted"]
        + funnel["wave1_regulatory_or_guidance_records_extracted"]
        + funnel["wave1_computational_models_extracted"]
        + funnel["wave1_correspondence_or_secondary_extracted"]
        + funnel["wave1_superseded_preprint_duplicates_extracted"]
        == 1176
    )
    assert (
        funnel["wave1_ceg_records_extracted_total"]
        + funnel["remaining_supporting_catalog_queued_for_later_waves"]
        == 3263
    )
    assert (
        funnel["wave1_ring1_empirical_core_extracted"] + funnel["wave1_new_empirical_studies_extracted"]
        == 1150
    )
    assert (
        funnel["wave1_ring1_core_records_extracted_total"] + funnel["wave1_non_r1_records_extracted_total"]
        == 1176
    )
    # Negative controls: computational/modeling-primary, correspondence, all 3 superseded preprints, VC-01 & Sweet M3020:O1
    for _cm_qid in ("M0223", "M0263", "M0372", "R1_0367"):
        assert ceg["studies"][_cm_qid]["record_class"] == "computational_model_or_simulation"
    assert ceg["studies"]["M0263"].get("has_experimental_subcomponent") is True
    assert ceg["studies"]["R1_0367"].get("has_experimental_subcomponent") is True
    assert ceg["studies"]["R1_0381"]["record_class"] == "correspondence_or_secondary_summary"
    assert ceg["studies"]["R1_0381"].get("same_patient_followup_of_queue_id") == "R1_0334"
    assert ceg["outcomes"]["R1_0381:O1"]["source_pin"]["source_type"] == "title_or_citation_metadata_only"
    assert ceg["outcomes"]["R1_0381:O2"]["source_pin"]["source_type"] == "title_or_citation_metadata_only"
    for _dup_qid in ("R1_0258", "R1_0298", "M0198"):
        assert ceg["studies"][_dup_qid]["record_class"] == "superseded_preprint_duplicate"
        assert all(e["primary_certification_eligible"] is False for e in ceg["edges"] if e["queue_id"] == _dup_qid)
    assert ceg["outcomes"]["M3020:O1"]["direction"] == "descriptive_or_neutral"
    assert ceg["outcomes"]["M0119:O2"]["stat_type"] == "descriptive_scalar"
    assert ceg["outcomes"]["M0119:O2"]["source_pin"]["source_type"] == "trial_registry_results"
    assert "114" in ceg["outcomes"]["M0119:O2"]["quantitative_summary"] and "0 participants analyzed" in ceg["outcomes"]["M0119:O2"]["quantitative_summary"]
    assert ceg["outcomes"]["M0119:O3"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["M0119:O3"]["source_pin"]["source_type"] == "trial_registry_results"
    assert ceg["outcomes"]["M0119:O4"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["M0119:O4"]["source_pin"]["source_type"] == "trial_registry_status"
    assert "whyStopped" in str(ceg["outcomes"]["M0119:O4"]["source_pin"].get("fulltext_section") or "")
    # Live recomputation of all 3 public-domain NCT02239354 VC-01 hashes + source-backed BANDIT M0155:O4/O6 hashes
    import hashlib as _hl_v15
    _vc01_canonical_sources = {
        "M0119:O2": '{"classes":[{"categories":[{"measurements":[{"groupId":"OG000","value":"114"}]}]}],"denoms":[{"counts":[{"groupId":"OG000","value":"19"}],"units":"Participants"}],"groups":[{"description":"Up to 2 VC-01-250\u2122 Combination Product implants\\n\\nUp to 6 VC-01-20 sentinel units\\n\\nVC-01\u2122 Combination Product: Biologic and Device","id":"OG000","title":"Cohort 1"}],"paramType":"NUMBER","populationDescription":"The Safety Analysis Set (SAS) will include all T1DM subjects who were enrolled into the study and in whom an implant surgery was attempted, regardless if any VC-01-250 units or sentinel-sized units were actually implanted. This includes all subjects in both Cohort 1 and Cohort 2.","reportingStatus":"POSTED","timeFrame":"Thru the Month 24 Visit","title":"Number of Adverse Events Reported During the Study.","type":"PRIMARY","unitOfMeasure":"events"}',
        "M0119:O3": '{"denoms":[{"counts":[{"groupId":"OG000","value":"0"}],"units":"Participants"}],"description":"The Full Analysis Set (FAS) is the intent-to-treat (ITT) set of subjects. This set is defined as all T1DM subjects who were enrolled into the study and received implantation of at least one VC-01-250 or sentinel unit on Study Day 1. The FAS (Cohort 2 T1DM subjects who meet the FAS criteria) will be used to analyze the primary efficacy endpoint.","groups":[{"description":"Up to 4 or 6 VC-01\u2122 Combination Product implants\\n\\nUp to 3 VC-01-20 sentinel units\\n\\nVC-01\u2122 Combination Product: Biologic and Device","id":"OG000","title":"Cohort 2"}],"populationDescription":"No subjects were enrolled in Cohort 2. Therefore, no data was collected for this endpoint.","reportingStatus":"POSTED","timeFrame":"Baseline to the Month 6 Visit","title":"Change in C-peptide","type":"PRIMARY"}',
        "M0119:O4": '{"overallStatus":"TERMINATED","whyStopped":"Insufficient functional product engraftment"}',
    }
    for _vc_oid, _vc_bytes_str in _vc01_canonical_sources.items():
        assert ceg["outcomes"][_vc_oid]["quote_sha256"] == _hl_v15.sha256(_vc_bytes_str.encode("utf-8")).hexdigest()
    assert ceg["outcomes"]["M0155:O4"]["source_pin"]["source_type"] == "abstract"
    assert ceg["outcomes"]["M0155:O6"]["source_pin"]["source_type"] == "abstract"
    assert ceg["outcomes"]["M0155:O6"]["stat_type"] == "qualitative_no_pvalue"
    _d42_path = REPO_ROOT / "data" / "dossiers" / "42627334.json"
    if _d42_path.is_file():
        _d42_ab = json.loads(_d42_path.read_text(encoding="utf-8"))["abstract"]
        for _b_oid in ("M0155:O4", "M0155:O6"):
            _m_sp = re.fullmatch(r"abstract\[(\d+):(\d+)\]", str(ceg["outcomes"][_b_oid]["source_pin"].get("dossier_span") or ""))
            assert _m_sp is not None, f"Missing dossier_span on {_b_oid}"
            _b_slice = _d42_ab[int(_m_sp.group(1)):int(_m_sp.group(2))]
            assert ceg["outcomes"][_b_oid]["quote_sha256"] == _hl_v15.sha256(_b_slice.encode("utf-8")).hexdigest()
    _vc01_edge = next(e for e in ceg["edges"] if e["edge_id"] == "M0119:L1")
    assert _vc01_edge["primary_certification_eligible"] is False
    assert "M0119:O2" in _vc01_edge["cited_outcome_tuples"]
    assert "March 24, 2022" in _vc01_edge["claim_summary"]
    # Shi M0144:O2 & M0349:O1 title-only stat_type & Pham M0060:O13 (mean_unspecified_spread) / M0060:L1 precision
    assert ceg["outcomes"]["M0144:O2"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["M0144:O2"]["source_pin"]["source_type"] == "title_or_citation_metadata_only"
    assert ceg["outcomes"]["M0349:O1"]["source_pin"]["source_type"] == "title_or_citation_metadata_only"
    assert ceg["outcomes"]["M0060:O13"]["stat_type"] == "mean_unspecified_spread"
    assert "243.53" in ceg["outcomes"]["M0060:O13"]["quantitative_summary"]
    assert ceg["outcomes"]["M0237:O4"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["M0237:O4"]["direction"] == "unfavorable_or_toxic"
    assert ceg["outcomes"]["M0237:O19"]["direction"] == "favorable_to_intervention"
    assert ceg["outcomes"]["M0237:O5"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["PMID:40976248:O1"]["stat_type"] == "adjusted_difference_95ci"
    assert ceg["outcomes"]["PMID:10531307:O2"]["source_pin"]["table_id"] == "Table V"
    assert "Table VI" not in ceg["outcomes"]["PMID:10531307:O2"]["quantitative_summary"]
    assert "2006-2009" in ceg["outcomes"]["PMID:29862647:O1"]["quantitative_summary"]
    assert "2011-2012" not in ceg["outcomes"]["PMID:29862647:O1"]["quantitative_summary"]
    assert "M2895:O2" in ceg["outcomes"] and "M2895:O3" in ceg["outcomes"]
    assert "TNF-alpha" in ceg["outcomes"]["M2895:O2"]["endpoint_name"]
    assert "LHR-A" in ceg["outcomes"]["PMID:10531307:O1"]["quantitative_summary"] and "sCR1" in ceg["outcomes"]["PMID:10531307:O1"]["quantitative_summary"]
    assert len(ceg["studies"]["PMID:16177096"]["variable_arm_ids"]) == 4
    assert "PMID:24222351:O2" in ceg["outcomes"] and "PMID:29862647:O2" in ceg["outcomes"]
    _edge_map = {e["edge_id"]: e for e in ceg["edges"]}
    assert "comparable to normoxia" not in _edge_map["M0060:L1"]["claim_summary"]
    assert _edge_map["M0349:L1"]["primary_certification_eligible"] is False
    for _sim_eid in (
        "R1_0365:L1", "R1_0365:L2", "M0035:L1", "M0035:L2", "M0105:L1",
        "M0397:L2", "M0099:L3", "M0139:L2", "M0216:L2", "M0226:L2", "M0245:L1",
        "M0277:L3", "M0454:L2", "M0693:L2", "M1490:L2", "M1746:L1", "M1756:L2", "M2700:L2",
    ):
        assert _edge_map[_sim_eid]["primary_certification_eligible"] is False
    for _pos_emp_eid in ("M0035:L3", "M0397:L1", "M0099:L2", "M0139:L1", "M1756:L1"):
        assert _edge_map[_pos_emp_eid]["primary_certification_eligible"] is True
    assert ceg["outcomes"]["R1_0365:O1"]["evidence_provenance"] == "mixed_empirical_and_computational"
    assert ceg["outcomes"]["M0397:O11"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0099:O24"]["evidence_provenance"] == "computational_model_or_simulation"
    assert "M0237:O19" in _edge_map["M0237:L1"]["cited_outcome_tuples"] and "M0237:O4" not in _edge_map["M0237:L1"]["cited_outcome_tuples"]
    assert "M0237:O4" in _edge_map["M0237:L2"]["cited_outcome_tuples"]
    for _dm_eid in ("M0067:L1", "M0067:L2", "M0067:L3", "M0086:L1", "M0086:L2", "M0086:L3", "M0086:L4", "PMID:41543848:L1"):
        assert _edge_map[_dm_eid]["applicability"] == "domain_mismatch"
        assert _edge_map[_dm_eid]["primary_certification_eligible"] is False
    assert "13 oncology" in _edge_map["M0067:L1"]["claim_summary"] and "1 autoimmune" in _edge_map["M0067:L1"]["claim_summary"]
    for _np_oid in ("M0086:O1", "M0086:O2", "M0086:O3", "M0067:O1", "M0067:O2", "M0067:O3"):
        assert ceg["outcomes"][_np_oid]["stat_type"] == "qualitative_no_pvalue"
    # M0172 modeled outcomes, all domain_mismatch & all noninformative edges must be ineligible
    assert ceg["outcomes"]["M0172:O9"]["evidence_provenance"] == "computational_model_or_simulation"
    assert _edge_map["M0172:L1"]["primary_certification_eligible"] is False
    assert _edge_map["M0172:L2"]["primary_certification_eligible"] is False
    assert all(e["primary_certification_eligible"] is False for e in ceg["edges"] if e.get("applicability") == "domain_mismatch")
    assert all(e["primary_certification_eligible"] is False for e in ceg["edges"] if e.get("relation") == "noninformative")
    # Gate 5 safety-switch edges M0069:L1 and M0088:L1 must map strictly to Gate 5 (not Gate 4)
    assert ceg_mod._classify_gates_and_designs(_edge_map["M0069:L1"]["target_gate_or_module"])[0] == ["5"]
    assert ceg_mod._classify_gates_and_designs(_edge_map["M0088:L1"]["target_gate_or_module"])[0] == ["5"]
    # Simulated-transfer edges on mixed studies (M0200:L1/L2, M1716:L2) must be ineligible while empirical edges stay eligible
    assert _edge_map["M0200:L1"]["primary_certification_eligible"] is False
    assert _edge_map["M0200:L2"]["primary_certification_eligible"] is False
    assert _edge_map["M0200:L3"]["primary_certification_eligible"] is True
    assert _edge_map["M1716:L1"]["primary_certification_eligible"] is True
    assert _edge_map["M1716:L2"]["primary_certification_eligible"] is False
    # Sweet M3020:O2/O3 dispersion & high-glucose absolute OCR precision
    assert ceg["outcomes"]["M3020:O2"]["stat_type"] == "mean_unspecified_spread"
    assert "at high glucose" in ceg["outcomes"]["M3020:O3"]["quantitative_summary"]
    # Verifier negative controls (IIP-0032 eligibility, provenance, and supersession guards)
    import copy as _copy
    import pytest as _pytest
    import tempfile as _tempfile
    _bad_narr_pin = _copy.deepcopy(ceg)
    _bad_narr_pin["outcomes"]["M0349:O1"]["source_pin"] = {"source_type": "abstract"}
    with _pytest.raises(ValueError, match="title-only"):
        ceg_mod.verify_ceg_payload(_bad_narr_pin, repo_root=REPO_ROOT)
    _bad_untag_sim = _copy.deepcopy(ceg)
    _bad_untag_sim["outcomes"]["M0105:O3"].pop("evidence_provenance", None)
    with _pytest.raises(ValueError, match="computational/simulated"):
        ceg_mod.verify_ceg_payload(_bad_untag_sim, repo_root=REPO_ROOT)
    for _v11_sim_oid in ("R1_0365:O1", "M0397:O11", "M0099:O24", "M0139:O4", "M1756:O8"):
        _bad_v11_untag = _copy.deepcopy(ceg)
        _bad_v11_untag["outcomes"][_v11_sim_oid].pop("evidence_provenance", None)
        with _pytest.raises(ValueError, match="computational/simulated"):
            ceg_mod.verify_ceg_payload(_bad_v11_untag, repo_root=REPO_ROOT)
    for _v11_sim_eid in ("R1_0365:L1", "M0397:L2", "M0099:L3"):
        _bad_v11_elig = _copy.deepcopy(ceg)
        _e_v11 = next(e for e in _bad_v11_elig["edges"] if e["edge_id"] == _v11_sim_eid)
        _e_v11["primary_certification_eligible"] = True
        _e_v11.pop("evidence_provenance", None)
        with _pytest.raises(ValueError, match="computational/simulated|mixed provenance"):
            ceg_mod.verify_ceg_payload(_bad_v11_elig, repo_root=REPO_ROOT)
    _bad_struct_stat = _copy.deepcopy(ceg)
    _bad_struct_stat["outcomes"]["M0086:O2"]["quantitative_summary"] = _bad_struct_stat["outcomes"]["M0086:O2"]["quantitative_summary"].replace("stat=qualitative_no_pvalue", "stat=qualitative_with_pvalue")
    with _pytest.raises(ValueError, match="disagrees with"):
        ceg_mod.verify_ceg_payload(_bad_struct_stat, repo_root=REPO_ROOT)
    _bad_sup = _copy.deepcopy(ceg)
    _bad_sup["studies"]["R1_0258"]["superseded_by_queue_id"] = "NO_SUCH_STUDY"
    with _pytest.raises(ValueError, match="superseded_by_queue_id"):
        ceg_mod.verify_ceg_payload(_bad_sup, repo_root=REPO_ROOT)
    _bad_pmid = _copy.deepcopy(ceg)
    _bad_pmid["studies"]["R1_0258"]["superseded_by_pmid"] = "00000000"
    with _pytest.raises(ValueError, match="superseded_by_pmid"):
        ceg_mod.verify_ceg_payload(_bad_pmid, repo_root=REPO_ROOT)
    _bad_cycle = _copy.deepcopy(ceg)
    _bad_cycle["studies"]["R1_0282"]["superseded_by_queue_id"] = "R1_0258"
    _bad_cycle["studies"]["R1_0282"]["superseded_by_pmid"] = "37905040"
    _bad_cycle["studies"]["R1_0282"]["record_class"] = "superseded_preprint_duplicate"
    with _pytest.raises(ValueError, match="supersession cycle"):
        ceg_mod.verify_ceg_payload(_bad_cycle, repo_root=REPO_ROOT)
    _bad_cite = _copy.deepcopy(ceg)
    _emp_edge = next(e for e in _bad_cite["edges"] if e["primary_certification_eligible"] is True)
    _emp_edge["cited_outcome_tuples"].append("R1_0381:O1")
    _emp_edge["cited_outcome_ids"].append("R1_0381:O1")
    with _pytest.raises(ValueError, match="non-empirical/duplicate"):
        ceg_mod.verify_ceg_payload(_bad_cite, repo_root=REPO_ROOT)
    _bad_sim = _copy.deepcopy(ceg)
    _emp_edge2 = next(e for e in _bad_sim["edges"] if e["primary_certification_eligible"] is True)
    _emp_edge2["cited_outcome_tuples"].append("M0200:O12")
    _emp_edge2["cited_outcome_ids"].append("M0200:O12")
    with _pytest.raises(ValueError, match="computational/simulated"):
        ceg_mod.verify_ceg_payload(_bad_sim, repo_root=REPO_ROOT)
    _bad_mixed = _copy.deepcopy(ceg)
    _emp_edge_m = next(e for e in _bad_mixed["edges"] if e["primary_certification_eligible"] is True and e["cited_outcome_tuples"])
    _bad_mixed["outcomes"][_emp_edge_m["cited_outcome_tuples"][0]]["evidence_provenance"] = "mixed_measured_and_simulated_model"
    with _pytest.raises(ValueError, match="mixed provenance"):
        ceg_mod.verify_ceg_payload(_bad_mixed, repo_root=REPO_ROOT)
    _bad_noninf = _copy.deepcopy(ceg)
    next(e for e in _bad_noninf["edges"] if e["edge_id"] == "M0144:L1")["primary_certification_eligible"] = True
    with _pytest.raises(ValueError, match="Noninformative edge"):
        ceg_mod.verify_ceg_payload(_bad_noninf, repo_root=REPO_ROOT)
    _bad_title = _copy.deepcopy(ceg)
    _m144_e = next(e for e in _bad_title["edges"] if e["edge_id"] == "M0144:L1")
    _m144_e["relation"] = "contextual_constraint"
    _m144_e["primary_certification_eligible"] = True
    with _pytest.raises(ValueError, match="needs_primary_fulltext|title-only"):
        ceg_mod.verify_ceg_payload(_bad_title, repo_root=REPO_ROOT)
    _bad_m349_elig = _copy.deepcopy(ceg)
    _m349_e = next(e for e in _bad_m349_elig["edges"] if e["edge_id"] == "M0349:L1")
    _m349_e["primary_certification_eligible"] = True
    with _pytest.raises(ValueError, match="Noninformative edge"):
        ceg_mod.verify_ceg_payload(_bad_m349_elig, repo_root=REPO_ROOT)
    _m349_e["relation"] = "contextual_constraint"
    with _pytest.raises(ValueError, match="title-only"):
        ceg_mod.verify_ceg_payload(_bad_m349_elig, repo_root=REPO_ROOT)
    _bad_dm = _copy.deepcopy(ceg)
    next(e for e in _bad_dm["edges"] if e["edge_id"] == "M0067:L1")["primary_certification_eligible"] = True
    with _pytest.raises(ValueError, match="domain_mismatch"):
        ceg_mod.verify_ceg_payload(_bad_dm, repo_root=REPO_ROOT)
    _bad_m172 = _copy.deepcopy(ceg)
    _m172_e = next(e for e in _bad_m172["edges"] if e["edge_id"] == "M0172:L1")
    _m172_e["primary_certification_eligible"] = True
    _m172_e.pop("evidence_provenance", None)
    with _pytest.raises(ValueError, match="computational/simulated"):
        ceg_mod.verify_ceg_payload(_bad_m172, repo_root=REPO_ROOT)
    _bad_bool = _copy.deepcopy(ceg)
    next(e for e in _bad_bool["edges"] if e["edge_id"] == "M0200:L1")["primary_certification_eligible"] = "false"
    with _pytest.raises(ValueError, match="strict bool"):
        ceg_mod.verify_ceg_payload(_bad_bool, repo_root=REPO_ROOT)
    _bad_loc = _copy.deepcopy(ceg)
    _bad_loc["edges"][0]["cited_outcome_ids"].append("O9999")
    with _pytest.raises(ValueError, match="unresolved outcome"):
        ceg_mod.verify_ceg_payload(_bad_loc, repo_root=REPO_ROOT)
    _bad_sp = _copy.deepcopy(ceg)
    _bad_sp["outcomes"][next(iter(_bad_sp["outcomes"]))]["source_pin"] = "abstract"
    with _pytest.raises(ValueError, match="source_pin must be a dict"):
        ceg_mod.verify_ceg_payload(_bad_sp, repo_root=REPO_ROOT)
    _bad_idx = _copy.deepcopy(ceg)
    _bad_idx["edges"][0]["link_index"] = 99
    with _pytest.raises(ValueError, match="link_index"):
        ceg_mod.verify_ceg_payload(_bad_idx, repo_root=REPO_ROOT)
    _bad_sf = _copy.deepcopy(ceg)
    _bad_sf["screening_funnel"]["screened_records_total"] = 3900
    with _pytest.raises(ValueError, match="screening_funnel mismatch"):
        ceg_mod.verify_ceg_payload(_bad_sf, repo_root=REPO_ROOT)
    # Importer negative controls (--import-raw-graph): non-bool "false" and broken/empty screening ledger
    with _tempfile.TemporaryDirectory() as _td:
        _tdp = Path(_td)
        _raw_g = {
            "study_contexts": ceg["studies"],
            "variable_nodes": ceg["arms"],
            "outcome_nodes": ceg["outcomes"],
            "edges": _copy.deepcopy(ceg["edges"]),
            "effector_nodes": ceg["effectors"],
            "graveyard_index": ceg["graveyard_index"],
        }
        _good_ledger = {
            "summary": {
                "total_studies": 3938,
                "ring1_existing_records_screened": 594,
                "new_candidate_records_screened": 3344,
                "overall_counts": {"ADVANCE": 3263, "PARK_CONTEXT": 522, "EXCLUDE_OFF_TOPIC": 153},
                "candidate_extraction_waves": {"Wave_2B_Supporting_Catalog_P2_P3": 2087},
            }
        }
        _pkts = _tdp / "pkts"
        _pkts.mkdir()
        (_tdp / "good_ledger.json").write_text(json.dumps(_good_ledger), encoding="utf-8")
        (_tdp / "empty_ledger.json").write_text("{}", encoding="utf-8")
        _bad_ledger = _copy.deepcopy(_good_ledger)
        _bad_ledger["summary"]["total_studies"] = 3900
        (_tdp / "bad_ledger.json").write_text(json.dumps(_bad_ledger), encoding="utf-8")
        (_tdp / "raw_good.json").write_text(json.dumps(_raw_g), encoding="utf-8")
        with _pytest.raises(ValueError, match="screening_ledger"):
            ceg_mod.import_and_materialize(_tdp / "raw_good.json", _tdp / "bad_ledger.json", _pkts, repo_root=REPO_ROOT)
        with _pytest.raises(ValueError, match="screening_ledger"):
            ceg_mod.import_and_materialize(_tdp / "raw_good.json", _tdp / "empty_ledger.json", _pkts, repo_root=REPO_ROOT)
        _raw_str_bool = _copy.deepcopy(_raw_g)
        next(e for e in _raw_str_bool["edges"] if e["edge_id"] == "M0200:L3")["primary_certification_eligible"] = "false"
        (_tdp / "raw_str_bool.json").write_text(json.dumps(_raw_str_bool), encoding="utf-8")
        with _pytest.raises(ValueError, match="strict bool"):
            ceg_mod.import_and_materialize(_tdp / "raw_str_bool.json", _tdp / "good_ledger.json", _pkts, repo_root=REPO_ROOT)
        _raw_st_pin = _copy.deepcopy(_raw_g)
        _raw_st_pin["study_contexts"] = _copy.deepcopy(ceg["studies"])
        _raw_st_pin["study_contexts"]["M0200"]["source_pin"] = {"source_type": "title_or_citation_metadata_only"}
        (_tdp / "raw_st_pin.json").write_text(json.dumps(_raw_st_pin), encoding="utf-8")
        _imp_ceg, _ = ceg_mod.import_and_materialize(_tdp / "raw_st_pin.json", _tdp / "good_ledger.json", _pkts, repo_root=REPO_ROOT)
        assert _imp_ceg["studies"]["M0200"]["source_pin"] == {"source_type": "title_or_citation_metadata_only"}
        assert next(e for e in _imp_ceg["edges"] if e["edge_id"] == "M0200:L3")["primary_certification_eligible"] is False
        # v12 Bug #2: preserve explicit computational record_class on M0223 across --import-raw-graph even under prose redaction
        _raw_m223 = _copy.deepcopy(_raw_g)
        _raw_m223["study_contexts"] = _copy.deepcopy(ceg["studies"])
        _raw_m223["study_contexts"]["M0223"]["biological_model"] = "[Structured biological_model (M0223) | text_sha256=" + "a" * 64 + "]"
        _raw_m223["study_contexts"]["M0223"]["immune_context"] = "[Structured immune_context (M0223) | text_sha256=" + "b" * 64 + "]"
        _raw_m223["study_contexts"]["M0223"]["follow_up"] = "[Structured follow_up (M0223) | text_sha256=" + "c" * 64 + "]"
        (_tdp / "raw_m223.json").write_text(json.dumps(_raw_m223), encoding="utf-8")
        _imp_m223, _ = ceg_mod.import_and_materialize(_tdp / "raw_m223.json", _tdp / "good_ledger.json", _pkts, repo_root=REPO_ROOT)
        assert _imp_m223["studies"]["M0223"]["record_class"] == "computational_model_or_simulation"
        assert all(e["primary_certification_eligible"] is False for e in _imp_m223["edges"] if e["queue_id"] == "M0223")
        ceg_mod.verify_ceg_payload(_imp_m223, repo_root=REPO_ROOT)
    # v12 Bug #2 shipped fixes & verifier negative controls (M0174:L2/O7-O8, M0080:O4/O11, M0155:L1, empty citations, backlink symmetry, title-only parent cross-link)
    assert ceg["outcomes"]["M0174:O7"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0174:O8"]["evidence_provenance"] == "computational_model_or_simulation"
    assert _edge_map["M0174:L2"]["primary_certification_eligible"] is False
    assert _edge_map["M0174:L1"]["primary_certification_eligible"] is False
    assert _edge_map["M0930:L1"]["primary_certification_eligible"] is False
    assert _edge_map["M0065:L1"]["primary_certification_eligible"] is False
    assert _edge_map["M1511:L2"]["primary_certification_eligible"] is False
    assert ceg["outcomes"]["M0174:O3"]["evidence_provenance"] == "mixed_empirical_and_computational"
    assert ceg["outcomes"]["M0930:O7"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0930:O12"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0065:O13"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M1511:O3"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0080:O4"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["M0080:O11"]["stat_type"] == "qualitative_no_pvalue"
    assert "week-24" in _edge_map["M0155:L1"]["claim_summary"] and "week-48" in _edge_map["M0155:L1"]["claim_summary"]
    assert "transient immunomodulation of T1D autoimmunity does not give durable protection" not in _edge_map["M0155:L1"]["claim_summary"]
    assert len(ceg["studies"]["PMID:24222351"]["variable_arm_ids"]) == 5 and "PMID:24222351:O3" in ceg["outcomes"]
    assert len(ceg["studies"]["PMID:10531307"]["variable_arm_ids"]) == 4
    _bad_m174 = _copy.deepcopy(ceg)
    _m174_e = next(e for e in _bad_m174["edges"] if e["edge_id"] == "M0174:L2")
    _m174_e["primary_certification_eligible"] = True
    _m174_e.pop("evidence_provenance", None)
    with _pytest.raises(ValueError, match="computational/simulated"):
        ceg_mod.verify_ceg_payload(_bad_m174, repo_root=REPO_ROOT)
    _bad_empty_elig = _copy.deepcopy(ceg)
    _m200_l3 = next(e for e in _bad_empty_elig["edges"] if e["edge_id"] == "M0200:L3")
    _m200_l3["cited_outcome_ids"] = []
    _m200_l3["cited_outcome_tuples"] = []
    with _pytest.raises(ValueError, match="must cite at least one outcome"):
        ceg_mod.verify_ceg_payload(_bad_empty_elig, repo_root=REPO_ROOT)
    _bad_stale_bwd = _copy.deepcopy(ceg)
    _m200_l3_b = next(e for e in _bad_stale_bwd["edges"] if e["edge_id"] == "M0200:L3")
    _m200_l3_b["primary_certification_eligible"] = False
    _m200_l3_b["cited_outcome_ids"] = []
    _m200_l3_b["cited_outcome_tuples"] = []
    with _pytest.raises(ValueError, match="asymmetric"):
        ceg_mod.verify_ceg_payload(_bad_stale_bwd, repo_root=REPO_ROOT)
    _bad_xlink_title = _copy.deepcopy(ceg)
    _bad_xlink_title["studies"]["M0155"]["source_pin"] = {"source_type": "title_or_citation_metadata_only"}
    for _e_m155 in _bad_xlink_title["edges"]:
        if _e_m155["queue_id"] == "M0155":
            _e_m155["primary_certification_eligible"] = False
    _m200_xl = next(e for e in _bad_xlink_title["edges"] if e["edge_id"] == "M0200:L3")
    _m200_xl["cited_outcome_ids"].append("M0155:O1")
    _m200_xl["cited_outcome_tuples"].append("M0155:O1")
    _bad_xlink_title["outcomes"]["M0155:O1"].setdefault("cited_by_edge_ids", []).append("M0200:L3")
    with _pytest.raises(ValueError, match="title-only study"):
        ceg_mod.verify_ceg_payload(_bad_xlink_title, repo_root=REPO_ROOT)
    # v13 Bug #2, #4, #5: Krych-Goldberg 1999/2005, M0155:O6, M1439:O3, duplicate edge_id, duplicate backlink, DOI-backed regulatory import after redaction
    assert "one intervening CCP" in ceg["arms"]["VAR:PMID:10531307:A4_ShortenedSpacer_OneCCP_Ctrl"]["label"]
    assert "unspaced" not in ceg["arms"]["VAR:PMID:10531307:A4_ShortenedSpacer_OneCCP_Ctrl"]["label"].lower()
    assert "58.0 +/- 0.2 SD" in ceg["outcomes"]["PMID:16177096:O1"]["quantitative_summary"] and "highest observed arm mean" in ceg["outcomes"]["PMID:16177096:O1"]["quantitative_summary"]
    assert "Fig. 2A" in ceg["outcomes"]["PMID:16177096:O2"]["quantitative_summary"]
    assert "week-24" in ceg["outcomes"]["M0155:O6"]["comparison_arms"] and "week-48" in ceg["outcomes"]["M0155:O6"]["comparison_arms"]
    assert "in vitro" in ceg["outcomes"]["M1439:O3"]["quantitative_summary"].lower() and ceg["outcomes"]["M1439:O3"]["source_pin"].get("figure_id") == "SI Appendix Fig. S7"
    # v14 Bug #2, #5, #7: Nair M0237:O19 separation, M0049:L1 & R1_0142:L2 computational exclusion, Yang M1439 4 arms & untreated diabetic Lewis rats, PMID:16177096:O1 p_value_or_ci null
    assert "matching Ring-1" not in ceg["outcomes"]["M0237:O19"]["quantitative_summary"]
    assert ceg["outcomes"]["M0049:O1"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0049:O2"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["M0049:O15"]["evidence_provenance"] == "computational_model_or_simulation"
    assert ceg["outcomes"]["R1_0142:O10"]["evidence_provenance"] == "computational_model_or_simulation"
    assert _edge_map["M0049:L1"]["primary_certification_eligible"] is False
    assert _edge_map["M0049:L2"]["primary_certification_eligible"] is True
    assert _edge_map["R1_0142:L2"]["primary_certification_eligible"] is False
    assert _edge_map["R1_0142:L1"]["primary_certification_eligible"] is True
    assert len(ceg["studies"]["M1439"]["variable_arm_ids"]) == 5 and "VAR:M1439:A4_single_EqC_layout" in ceg["arms"] and "VAR:M1439:A5_untreated_diabetic_lewis_nodevice" in ceg["arms"]
    assert "M1439:O4" in ceg["outcomes"] and "A5_untreated_diabetic_lewis_nodevice" in ceg["outcomes"]["M1439:O4"]["comparison_arms"]
    assert "A2_human_islet_lewis_noflow" in ceg["outcomes"]["M1439:O2"]["comparison_arms"] and "alpha < 0.001" not in str(ceg["outcomes"]["M1439:O2"].get("p_value_or_ci") or "")
    assert "alpha < 0.001" in str(ceg["outcomes"]["M1439:O4"].get("p_value_or_ci") or "")
    assert "A4_single_EqC_layout" in ceg["outcomes"]["M1439:O3"]["comparison_arms"] and "dextran" in ceg["outcomes"]["M1439:O3"]["quantitative_summary"].lower()
    assert "solute release from the second eqc" not in ceg["outcomes"]["M1439:O3"]["quantitative_summary"].lower()
    assert len(ceg["studies"]["M2895"]["variable_arm_ids"]) == 6 and "VAR:M2895:A6_BMDM_hIAPP_Fig4D" in ceg["arms"]
    assert "A6_BMDM_hIAPP_Fig4D" in ceg["outcomes"]["M2895:O3"]["comparison_arms"] and "A4_BMDM_hIAPP_Ctrl" not in ceg["outcomes"]["M2895:O3"]["comparison_arms"]
    assert "il-1r antagonist (il-1ra) anakinra" not in json.dumps(ceg["arms"]["VAR:M2895:A3_BMDM_hIAPP_Anakinra"]).lower()
    assert "cell-free" not in json.dumps(ceg["studies"]["PMID:10531307"]).lower() and "sheep erythrocytes" in json.dumps(ceg["studies"]["PMID:10531307"]).lower()
    assert "cell-free" not in json.dumps(ceg["studies"]["PMID:16177096"]).lower() and "sheep erythrocytes" in json.dumps(ceg["studies"]["PMID:16177096"]).lower()
    assert ceg["outcomes"]["M0237:O6"]["comparison_arms"].startswith("A1_eBC vs A3_NEC")
    assert "60 baricitinib, 31 placebo" in str(ceg["outcomes"]["M0155:O2"].get("n_summary") or "") and "0.35" in str(ceg["outcomes"]["M0155:O2"].get("comparator_value") or "")
    assert ceg["outcomes"]["M0002:O1"]["manuscript_tier"] == "T5_Rodent_or_CellLine"
    assert ceg["outcomes"]["M0002:O6"]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["PMID:16177096:O1"].get("p_value_or_ci") is None
    # Negative controls: flipping ONLY M0155:L3 to eligible MUST fail; dropping provenance on curated model outcomes MUST fail
    _bad_m155_l3 = _copy.deepcopy(ceg)
    next(e for e in _bad_m155_l3["edges"] if e["edge_id"] == "M0155:L3")["primary_certification_eligible"] = True
    with _pytest.raises(ValueError, match="Curated-ineligible edge M0155:L3"):
        ceg_mod.verify_ceg_payload(_bad_m155_l3, repo_root=REPO_ROOT)
    for _mod_oid in ("M0049:O1", "R1_0142:O10", "M0174:O3", "M0930:O7", "M0065:O13", "M1511:O3"):
        _bad_mod_oc = _copy.deepcopy(ceg)
        _bad_mod_oc["outcomes"][_mod_oid].pop("evidence_provenance", None)
        with _pytest.raises(ValueError, match="computational/simulated"):
            ceg_mod.verify_ceg_payload(_bad_mod_oc, repo_root=REPO_ROOT)
    # Negative control: zero _SourceGuard >=6-word overlaps on M2895 & M1439 against dossiers when present
    for _sg_qid, _sg_pmid in (("M2895", "21813778"), ("M1439", "34504013")):
        _sg_dp = REPO_ROOT / "data" / "dossiers" / f"{_sg_pmid}.json"
        if _sg_dp.is_file():
            _sg_dos = json.loads(_sg_dp.read_text(encoding="utf-8"))
            _sg_ft = _sg_dos.get("fulltext_data_sections") or {}
            _sg_txt = " ".join([
                _sg_dos.get("title") or "", _sg_dos.get("abstract") or "",
                *[p for p in (_sg_ft.get("results_paragraphs") or []) if isinstance(p, str)],
                *[p for p in (_sg_ft.get("methods_paragraphs") or []) if isinstance(p, str)],
                *[f.get("caption") or "" for f in (_sg_ft.get("figure_legends") or []) if isinstance(f, dict)],
            ])
            _sg_ws = ceg_mod._words(_sg_txt)
            _sg_sh = {tuple(_sg_ws[i:i+6]) for i in range(max(0, len(_sg_ws) - 5))}
            _sg_guard = ceg_mod._SourceGuard({_sg_qid: _sg_sh}, [_sg_txt])
            for _ak, _av in ceg["arms"].items():
                if _av.get("queue_id") == _sg_qid:
                    for _iv in (_av.get("edits_or_interventions") or []):
                        assert not _sg_guard.overlaps(_sg_qid, str(_iv)), f"Overlap in {_ak}: {_iv}"
                    assert not _sg_guard.overlaps(_sg_qid, str(_av.get("label") or "")), f"Overlap in {_ak} label"
            for _ok, _ov in ceg["outcomes"].items():
                if _ov.get("queue_id") == _sg_qid:
                    assert not _sg_guard.overlaps(_sg_qid, str(_ov.get("quantitative_summary") or "")), f"Overlap in {_ok}"
    _bad_dup_edge = _copy.deepcopy(ceg)
    _m69_e = next(e for e in _bad_dup_edge["edges"] if e["queue_id"] == "M0069")
    _m69_clone = _copy.deepcopy(_m69_e)
    _m69_clone["link_index"] = _bad_dup_edge["studies"]["M0069"]["links_count"]
    _bad_dup_edge["edges"].append(_m69_clone)
    _bad_dup_edge["studies"]["M0069"]["links_count"] += 1
    _bad_dup_edge["counts"]["edges"] += 1
    with _pytest.raises(ValueError, match="Duplicate edge_id|Expected 2,483"):
        ceg_mod.verify_ceg_payload(_bad_dup_edge, repo_root=REPO_ROOT)
    # Also test duplicate edge_id preserving exact total edge count (2,483)
    _bad_dup_edge2 = _copy.deepcopy(ceg)
    _e_idx_drop = next(i for i, e in enumerate(_bad_dup_edge2["edges"]) if e["edge_id"] == "M0088:L1")
    _dropped = _bad_dup_edge2["edges"].pop(_e_idx_drop)
    _bad_dup_edge2["studies"]["M0088"]["links_count"] -= 1
    for _t_d in _dropped.get("cited_outcome_tuples") or []:
        _bad_dup_edge2["outcomes"][_t_d]["cited_by_edge_ids"] = [x for x in _bad_dup_edge2["outcomes"][_t_d].get("cited_by_edge_ids", []) if x != _dropped["edge_id"]]
    _m69_e2 = next(e for e in _bad_dup_edge2["edges"] if e["queue_id"] == "M0069")
    _m69_c2 = _copy.deepcopy(_m69_e2)
    _m69_c2["link_index"] = _bad_dup_edge2["studies"]["M0069"]["links_count"]
    _bad_dup_edge2["edges"].append(_m69_c2)
    _bad_dup_edge2["studies"]["M0069"]["links_count"] += 1
    with _pytest.raises(ValueError, match="Duplicate edge_id"):
        ceg_mod.verify_ceg_payload(_bad_dup_edge2, repo_root=REPO_ROOT)
    _bad_dup_bwd = _copy.deepcopy(ceg)
    _oid_bwd = next(oid for oid, o in _bad_dup_bwd["outcomes"].items() if o.get("cited_by_edge_ids"))
    _bad_dup_bwd["outcomes"][_oid_bwd]["cited_by_edge_ids"].append(_bad_dup_bwd["outcomes"][_oid_bwd]["cited_by_edge_ids"][0])
    with _pytest.raises(ValueError, match="duplicate edge_id"):
        ceg_mod.verify_ceg_payload(_bad_dup_bwd, repo_root=REPO_ROOT)
    # v13 Bug #4: DOI-backed regulatory_or_guidance_document survives --import-raw-graph prose redaction and passes verify_ceg_payload
    import json as _json_reg, tempfile as _tf_reg
    with _tf_reg.TemporaryDirectory() as _td_reg:
        _tdp_reg = Path(_td_reg)
        _pkts_reg = _tdp_reg / "pkts"
        _pkts_reg.mkdir()
        (_pkts_reg / "test_packets.json").write_text(_json_reg.dumps([{"queue_id": "M0223", "abstract": "Regulatory guidance on synthetic islet device validation for clinical trial submission."}]), encoding="utf-8")
        _raw_reg = _copy.deepcopy(ceg)
        _raw_reg["study_contexts"] = _raw_reg.pop("studies")
        _raw_reg["variable_nodes"] = _raw_reg.pop("arms")
        _raw_reg["outcome_nodes"] = _raw_reg.pop("outcomes")
        _raw_reg["effector_nodes"] = _raw_reg.pop("effectors")
        _reg_qid = "M0223"
        _raw_reg["study_contexts"][_reg_qid]["id_type"] = "DOI"
        _raw_reg["study_contexts"][_reg_qid]["id_value"] = "10.1000/synthetic.regulatory.guidance"
        _raw_reg["study_contexts"][_reg_qid]["record_class"] = "regulatory_or_guidance_document"
        _raw_reg["study_contexts"][_reg_qid]["biological_model"] = "Regulatory guidance on synthetic islet device validation for clinical trial submission."
        _rp_reg = _tdp_reg / "raw_reg.json"
        _lp_reg = _tdp_reg / "ledger_reg.json"
        _rp_reg.write_text(_json_reg.dumps(_raw_reg), encoding="utf-8")
        _lp_reg.write_text(_json_reg.dumps({"summary": {"total_studies": 3938, "ring1_existing_records_screened": 594, "new_candidate_records_screened": 3344, "overall_counts": {"ADVANCE": 3263, "PARK_CONTEXT": 522, "EXCLUDE_OFF_TOPIC": 153}, "candidate_extraction_waves": {"Wave_2B_Supporting_Catalog_P2_P3": 2087}}}), encoding="utf-8")
        _imp_reg, _ = ceg_mod.import_and_materialize(_rp_reg, _lp_reg, _pkts_reg, repo_root=REPO_ROOT)
        assert _imp_reg["studies"][_reg_qid]["biological_model"].startswith("[Structured ")
        assert _imp_reg["studies"][_reg_qid]["record_class"] == "regulatory_or_guidance_document"
        assert all(e["primary_certification_eligible"] is False for e in _imp_reg["edges"] if e["queue_id"] == _reg_qid)
        ceg_mod.verify_ceg_payload(_imp_reg, repo_root=REPO_ROOT)

    # Rights & token hygiene: zero inline verbatim quotes, valid SHA-256 digests, zero invalid redaction markers
    for oid, out_node in ceg["outcomes"].items():
        assert "verbatim_quote" not in out_node or out_node["verbatim_quote"] is None, oid
        if "quote_sha256" in out_node:
            assert re.fullmatch(r"[0-9a-f]{64}", out_node["quote_sha256"]), oid

    assert manifest_mod.public_token_errors(ceg) == []
    assert manifest_mod.public_token_errors(atlas) == []

    # Deterministic source-free replay of the Layer-1 CEG atlas from context_evidence_graph_v1.json
    replayed_atlas = ceg_mod.build_layer1_atlas_from_ceg(ceg)
    replayed_bytes = (json.dumps(replayed_atlas, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    assert atlas_bytes == replayed_bytes
    assert atlas["summary"]["ceg_records_total"] == 1176
    assert atlas["summary"]["primary_empirical_studies_extracted"] == 1150
    assert atlas["summary"]["trial_registry_protocols_in_ceg"] == 11
    assert atlas["summary"]["regulatory_or_guidance_records_in_ceg"] == 7
    assert atlas["summary"]["computational_models_in_ceg"] == 4
    assert atlas["summary"]["correspondence_or_secondary_records_in_ceg"] == 1
    assert atlas["summary"]["superseded_preprint_duplicates_in_ceg"] == 3
    assert atlas["summary"]["ring1_core_records_in_ceg_total"] == 404
    assert atlas["summary"]["ring1_empirical_core_studies_in_ceg"] == 400
    assert atlas["summary"]["new_expansion_empirical_studies_in_ceg"] == 750
    assert atlas["summary"]["non_r1_ceg_records_total"] == 772
    assert atlas["record_class_distribution"] == {
        "computational_model_or_simulation": 4,
        "correspondence_or_secondary_summary": 1,
        "empirical_primary_study": 1150,
        "regulatory_or_guidance_document": 7,
        "superseded_preprint_duplicate": 3,
        "trial_registry_protocol": 11,
    }
    assert atlas["summary"]["experimental_and_control_arms_total"] == 7740
    assert atlas["summary"]["structured_outcomes_total"] == 14167
    assert atlas["summary"]["blueprint_links_total"] == 2483
    assert atlas["summary"]["ring1_scoring_baseline_studies_unchanged"] == 594
    assert atlas["summary"]["inline_verbatim_quotes_shipped"] == 0
    assert atlas["summary"]["unresolved_outcome_citations"] == 0

    # Manifest alignment when manifest.json is present
    manifest_path = REPO_ROOT / "manifest.json"
    if manifest_path.is_file():
        mf = json.loads(manifest_path.read_text(encoding="utf-8"))
        mf_counts = mf["counts"]
        assert mf_counts["extracted_studies_total"] == 594
        assert mf_counts["extracted_atomic_effects_total"] == 3158
        assert mf_counts["ceg_studies_total"] == 1176
        assert mf_counts["ceg_empirical_primary_studies_total"] == 1150
        assert mf_counts["ceg_trial_registry_protocols_total"] == 11
        assert mf_counts["ceg_regulatory_or_guidance_records_total"] == 7
        assert mf_counts["ceg_computational_models_total"] == 4
        assert mf_counts["ceg_correspondence_secondary_total"] == 1
        assert mf_counts["ceg_superseded_preprint_duplicates_total"] == 3
        assert mf_counts["ceg_arms_total"] == 7740
        assert mf_counts["ceg_outcomes_total"] == 14167
        assert mf_counts["ceg_blueprint_links_total"] == 2483
