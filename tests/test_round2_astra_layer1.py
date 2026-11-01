"""Independent third-pass Layer-1 controls; synthetic cases are not measurements."""
from __future__ import annotations

import copy
import json
import math
import os
from pathlib import Path
import random
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import build_cross_study_graph as graph
from scripts import extract_ring1_corpus as ext
from scripts import rank_islet1_stacks as rank


def study_fixture():
    return {
        "study_metadata": {"pmid": "12345678", "source_text_type": "Abstract_Only", "max_biological_fidelity_tier": "Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse"},
        "arms": [{"gate2_3_4_genetic_and_epitope_edits": [{"gene_symbol_hgnc": "B2M", "edit_modality": "CRISPR_Cas9_KO"}]}],
        "atomic_effects": [{"effect_id": "fixture1", "focal_factor_X": "B2M knockout", "verbatim_evidence_quote": "B2M knockout improved cell viability.", "outcome_variable_Y": "cell viability", "effect_direction": "Protective_or_Improved", "secretome_or_differentiation_penalty_observed": False, "log_response_ratio_or_smd": None}],
    }


def utility_args(profiles):
    return ["ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", "AUTO_NONE", "CYTO_CXCL10_SWAP_PDL1_SERPINB9", "ER_NONE", next(iter(rank.NICHE_VASCULAR_MODULES)), "SAF_AAVS1_CAG_ICASP9_BIALLELIC", profiles]


@pytest.mark.parametrize("qualifier", [">", "<", "≥", "≤", "≈", "~", "approximately ", "roughly ", "at least ", "at most ", "more than ", "less than ", "about ", "around ", "up to ", "nearly ", "almost ", "over ", "under "])
@pytest.mark.parametrize("magnitude", ["2-fold", "2.5 ± 0.3-fold"])
def test_qualified_fold_is_not_promoted_to_exact_lrr(qualifier, magnitude):
    result = ext.parse_quantitative_numbers(f"B2M knockout improved cell viability by {qualifier}{magnitude} compared with controls.", "Protective_or_Improved")
    assert result[:4] == (None, None, None, None)
    assert result[-1] is None
    assert "quarantined" in result[4]


def test_exact_fold_positive_control_and_current_approximate_corpus_effect():
    exact = ext.parse_quantitative_numbers("B2M knockout increased cell viability by 2-fold compared with controls (p < 0.05).", "Protective_or_Improved")
    assert exact[0] == 2
    assert exact[-1] == pytest.approx(math.log(2), abs=5e-5)
    document = json.loads((ROOT / "data/extracted/31883920.json").read_text())
    approximate = next(e for e in document["atomic_effects"] if "roughly 2-fold" in e["verbatim_evidence_quote"])
    assert ext.parse_quantitative_numbers(approximate["verbatim_evidence_quote"], approximate["effect_direction"])[-1] is None


@pytest.mark.parametrize("source", ["p=1e-1000", "p<1e-1000", "p=1 × 10^-1000", "p=1 × 10^99999999"])
def test_unrepresentable_p_value_stays_unknown(source):
    assert ext._pvalue_tokens(source) == []
    assert ext.parse_pvalue_and_n(source)[0] is None


@pytest.mark.parametrize("source,expected", [("p=0.03", .03), ("p=1e-200", 1e-200), ("p=0", 0.)])
def test_representable_source_p_values_preserved(source, expected):
    assert ext.parse_pvalue_and_n(source)[0] == expected


@pytest.mark.parametrize("field,value", [
    ("effect_id", ""), ("effect_id", None), ("effect_direction", "Protective_TYPO"),
    ("secretome_or_differentiation_penalty_observed", "false"),
    ("secretome_or_differentiation_penalty_observed", 0),
    ("focal_factor_X", []), ("verbatim_evidence_quote", {}), ("outcome_variable_Y", 17),
])
def test_ranker_rejects_malformed_evidence_before_gene_attribution(field, value):
    fixture = study_fixture()
    fixture["atomic_effects"][0][field] = value
    with pytest.raises(ValueError):
        rank.build_target_evidence_profiles([fixture], {})


def test_duplicate_effect_identity_cannot_multiply_counts():
    fixture = study_fixture()
    fixture["atomic_effects"].append(copy.deepcopy(fixture["atomic_effects"][0]))
    with pytest.raises(ValueError, match="unique"):
        rank.build_target_evidence_profiles([fixture], {})


def test_all_profile_counts_are_validated_independent_of_mapping_order():
    for rows in [[("INS", {"n_empirical_studies": 1}), ("UNUSED", {"n_empirical_studies": "false"})], [("UNUSED", {"n_empirical_studies": "false"}), ("INS", {"n_empirical_studies": 1})]]:
        with pytest.raises(ValueError, match="Empirical study counts"):
            rank.compute_series_system_utility(*utility_args(dict(rows)))


def test_very_large_finite_uncertainty_yields_zero_floor_not_overflow():
    result = rank.compute_series_system_utility(*utility_args({"B2M": {"n_empirical_studies": 1, "posterior_se": 1e200}}))
    assert result["lower_10pct_cred_bound"] == 0
    assert math.isfinite(result["series_system_joint_survival_probability_unrounded"])


def test_overflowing_evidence_cannot_turn_cancellation_into_positive_saturation():
    profiles = {"X": {"net_evidence_and_constraint_score": 1e308}, "Y": {"net_evidence_and_constraint_score": -1e308}}
    for genes in [["X", "Y"], ["Y", "X"]]:
        with pytest.raises(ValueError, match="finite"):
            rank.signed_gene_signal(profiles, genes, scale=100)


@pytest.mark.parametrize("value", [True, "0.15", float("nan"), float("inf"), 10**500])
def test_physical_prior_domain_rejected_before_float_coercion(monkeypatch, value):
    monkeypatch.setitem(rank.ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["locus_dsb_penalty_rate_lambda_K"], "point_estimate", value)
    with pytest.raises(ValueError):
        rank.compute_series_system_utility(*utility_args({"INS": {"n_empirical_studies": 1}}))


@pytest.mark.parametrize("value", [True, "0.35", float("nan"), -1.])
def test_floor_registry_rejects_invalid_dispersion_and_discount(monkeypatch, value):
    params = rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"]["lower_10pct_cred_bound_hyperparameters"]["parameter_value"]
    monkeypatch.setitem(params, "default_gene_se", value)
    with pytest.raises(ValueError):
        rank.compute_series_system_utility(*utility_args({"INS": {"n_empirical_studies": 1}}))


@pytest.mark.parametrize("value", [True, "0.2", float("nan"), -1.])
def test_cd47_penalty_is_not_truthy_coerced_or_negative(value):
    with pytest.raises(ValueError):
        rank.compute_series_system_utility(*utility_args({"CD47": {"n_empirical_studies": 1, "weighted_deleterious_penalty_score": value}}))


def test_exact_direct_score_replays_all_feasible_production_stacks(monkeypatch):
    # Deliberately asymmetric cohort weights catch a hard-coded 0.5 / 0.5.
    secondary = rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["secondary_ranking_tie_break_and_display_constants"]
    monkeypatch.setitem(secondary["composite_utility_pediatric_weight"], "parameter_value", .2)
    monkeypatch.setitem(secondary["composite_utility_adult_weight"], "parameter_value", .8)
    for candidate in rank.run_combinatorial_multi_gate_search({"INS": {"n_empirical_studies": 1}}):
        expected = candidate["composite_series_system_utility_unrounded"]
        # The display field is intentionally corrupted. It must not contribute.
        candidate["utility_pediatric_early_onset_t1d"]["manufacturing_locus_penalty_factor"] = .1234
        assert rank._direct_output_score(candidate, [0.] * 6) == pytest.approx(expected, rel=2e-15, abs=1e-17)


def test_direct_score_zero_evidence_and_misaligned_gates_fail_closed():
    candidate = rank.run_combinatorial_multi_gate_search({})[0]
    assert rank._direct_output_score(candidate, [0.] * 6) == 0
    with pytest.raises(ValueError, match="aligned"):
        rank._direct_output_score(candidate, [0.] * 5)


def test_opentargets_provenance_survives_profile_construction():
    profile = rank.build_target_evidence_profiles([], {"B2M": {"opentargets_t1d_l2g_score": .02, "opentargets_score_source": "SYNTHETIC_SOURCE_IDENTIFIER"}})["B2M"]
    assert profile["opentargets_prior_source"] == "SYNTHETIC_SOURCE_IDENTIFIER"


def test_graveyard_same_size_mutation_cannot_hide_behind_restored_mtime(monkeypatch, tmp_path):
    target = tmp_path / "graveyard.json"
    data = rank.GRAVEYARD_PATH.read_text()
    target.write_text(data)
    monkeypatch.setattr(rank, "GRAVEYARD_PATH", target)
    rank._load_graveyard_entries_by_id()
    stat = target.stat()
    changed = data.replace('"severity_level": "BLOCK"', '"severity_level": "INFO "', 1)
    assert len(changed) == len(data) and changed != data
    target.write_text(changed)
    os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    with pytest.raises(ValueError):
        rank._load_graveyard_entries_by_id()


@pytest.mark.parametrize("source,stage", [(None, None), ("Unspecified_In_Source_Abstract", "Primary_Mature_Islet"), ("Mixed_Arm_Contexts", "Stage7_MetabolicallyMature_SC_Islet"), ("hESC", None)])
def test_equal_unknown_contexts_do_not_establish_biological_match(source, stage):
    context = {"cell_source_class": source, "canonical_differentiation_stage_evaluated": stage}
    assert not graph._has_matched_cell_source_and_stage(context, context)


def test_known_context_match_positive_control():
    context = {"cell_source_class": "hESC", "canonical_differentiation_stage_evaluated": "Stage7_MetabolicallyMature_SC_Islet"}
    assert graph._has_matched_cell_source_and_stage(context, dict(context))




def test_graph_rejects_filename_identity_mismatch(monkeypatch, tmp_path):
    record = json.loads((ROOT / "data/extracted/38508139.json").read_text())
    (tmp_path / "99999999.json").write_text(json.dumps(record))
    monkeypatch.setattr(graph, "EXTRACTED_DIR", tmp_path)
    with pytest.raises(ValueError, match="filename and study PMID"):
        graph.main()
    assert not (tmp_path / "cross_study_discovery_graph.json").exists()


@pytest.mark.parametrize("checkpoint", ['{"pmid":"1234"}\nnot-json\n', '{"pmid":"1234"}\n{"pmid":"1234"}\n', '{"pmid":null}\n'])
def test_harvester_rejects_broken_resume_before_network_or_writes(monkeypatch, tmp_path, checkpoint):
    from scripts import harvest_corpus as harvest
    target = tmp_path / "checkpoint.jsonl"
    target.write_text(checkpoint)
    def forbidden_network(*args, **kwargs):
        pytest.fail("Network must not run against an invalid resume checkpoint")
    monkeypatch.setattr(harvest.urllib.request, "urlopen", forbidden_network)
    with pytest.raises(ValueError, match="Cannot resume malformed corpus"):
        harvest.efetch_metadata(["1234", "5678"], {}, target)
    assert target.read_text() == checkpoint


def test_harvester_resume_repairs_only_missing_jsonl_separator(monkeypatch, tmp_path):
    from scripts import harvest_corpus as harvest
    import io
    target = tmp_path / "checkpoint.jsonl"
    target.write_text('{"pmid":"1234"}')
    xml = b'<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>5678</PMID><Article><ArticleTitle>Fixture</ArticleTitle></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>'
    monkeypatch.setattr(harvest.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(xml))
    monkeypatch.setattr(harvest.time, "sleep", lambda *a: None)
    assert harvest.efetch_metadata(["1234", "5678"], {}, target) == 2
    assert [json.loads(line)["pmid"] for line in target.read_text().splitlines()] == ["1234", "5678"]


@pytest.mark.parametrize("second_page", [{"count":"3", "idlist": []}, {"count":"3", "idlist": ["1"]}])
def test_harvester_incomplete_or_repeated_pagination_cannot_claim_complete(monkeypatch, second_page):
    from scripts import harvest_corpus as harvest
    import io
    pages = iter([{"count": "3", "idlist": ["1", "2"]}, second_page])
    monkeypatch.setattr(harvest.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(json.dumps({"esearchresult": next(pages)}).encode()))
    monkeypatch.setattr(harvest.time, "sleep", lambda *a: None)
    with pytest.raises(ValueError, match="pagination"):
        harvest.esearch_pmids("synthetic query", retmax=2)


def test_harvester_complete_pagination_positive_control(monkeypatch):
    from scripts import harvest_corpus as harvest
    import io
    pages = iter([{"count": "3", "idlist": ["1", "2"]}, {"count": "3", "idlist": ["3"]}])
    monkeypatch.setattr(harvest.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(json.dumps({"esearchresult": next(pages)}).encode()))
    monkeypatch.setattr(harvest.time, "sleep", lambda *a: None)
    assert harvest.esearch_pmids("synthetic query", retmax=2) == ["1", "2", "3"]
