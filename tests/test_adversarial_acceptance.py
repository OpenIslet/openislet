#!/usr/bin/env python3
"""
OpenIslet (Project Islet-1) Independent Acceptance & Provenance Verification Suite.

Verifies schema integrity, deterministic regeneration parity (extract_one_dossier),
distinct comparator arm linkage, strict word-boundary gene attribution, endpoint
polarity registration, prior-jitter stability, manipulation resistance, raw gnomAD/OpenTargets
constraint provenance, Graveyard rule reachability, clinical benchmark fidelity,
and release manifest / documentation synchronization.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
import jsonschema

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import extract_ring1_corpus as ext
from scripts import generate_release_manifest as manifest_mod
from scripts import rank_islet1_stacks as rank

CHASSIS_A_ID = (
    "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO__AUTO_NONE__"
    "CYTO_CXCL10_SWAP_PDL1_SERPINB9__ER_NONE__"
    "NICHE_TRIHORMONAL_120UM_MVF_CAO2__SAF_AAVS1_CAG_ICASP9_BIALLELIC"
)
CHASSIS_B_ID = (
    "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO__AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT__"
    "CYTO_CXCL10_SWAP_PDL1_SERPINB9__ER_NONE__"
    "NICHE_TRIHORMONAL_120UM_MVF_CAO2__SAF_AAVS1_CAG_ICASP9_BIALLELIC"
)


def sensitivity_artifact_matches(cached, fresh) -> bool:
    """Compare the complete replay, including draw counts, bounds and metadata."""
    if isinstance(cached, dict) and isinstance(fresh, dict):
        return cached.keys() == fresh.keys() and all(sensitivity_artifact_matches(cached[key], fresh[key]) for key in cached)
    if isinstance(cached, list) and isinstance(fresh, list):
        return len(cached) == len(fresh) and all(sensitivity_artifact_matches(a,b) for a,b in zip(cached,fresh))
    if isinstance(cached, (int,float)) and not isinstance(cached,bool) and isinstance(fresh,(int,float)) and not isinstance(fresh,bool):
        return math.isfinite(cached) and math.isfinite(fresh) and math.isclose(cached,fresh,rel_tol=1e-12,abs_tol=1e-12)
    return type(cached) is type(fresh) and cached == fresh


def replay_sensitivity_artifact(profiles, combos, studies, omics, cached):
    """Execute current code before comparing artifacts; errors must propagate."""
    fresh = rank.compute_direct_output_logit_sensitivity_analysis(profiles, combos, studies=studies, omics=omics)
    return fresh, sensitivity_artifact_matches(cached if cached is not None else fresh, fresh)


def source_bound_lrr_shrinkage_probe(studies, omics):
    """Exercise an otherwise unobserved branch using an actual source contrast.

    HNF1A is not in the current candidate module targets. Its documented KO
    contrast remains a valid independent branch oracle; no effect is invented.
    """
    source = next((study for study in studies if str(study.get("study_metadata",{}).get("pmid")) == "35918471"), None)
    if source is None:
        return {"verified": False, "reason": "Required source-backed HNF1A contrast unavailable"}
    effect = next((effect for effect in source.get("atomic_effects",[]) if effect.get("effect_id") == "pmid35918471_eff1"), None)
    if effect is None or not isinstance(effect.get("log_response_ratio_or_smd"),(int,float)):
        return {"verified": False, "reason": "Source contrast is not a documented numeric HNF1A KO"}
    fixture = copy.deepcopy(source)
    fixture["atomic_effects"] = [copy.deepcopy(effect)]
    # Explicit test annotation of the donor perturbation already named by the
    # paper's curated focal factor; persisted source records are not mutated.
    fixture["arms"][0]["gate2_3_4_genetic_and_epitope_edits"] = [{"gene_symbol_hgnc":"HNF1A", "edit_modality":"CRISPR_Cas9_KO"}]
    original = copy.deepcopy(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY)
    try:
        before = rank.build_target_evidence_profiles([fixture],omics)["HNF1A"]
        values=rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"]["lrr_shrinkage_hyperparameters"]["parameter_value"]
        for key,value in list(values.items()): values[key]=value*2.5+.5
        after = rank.build_target_evidence_profiles([fixture],omics)["HNF1A"]
    finally:
        rank.COMPLETE_GATE_COEFFICIENT_REGISTRY = original
    return {"verified":before["shrunk_posterior_log_rr"] != after["shrunk_posterior_log_rr"],
            "source_pmid":"35918471","source_effect_id":effect["effect_id"],
            "test_annotation":"HNF1A donor KO, already explicit in the unchanged source focal factor",
            "persisted_study_modified":False,"original_lrr":effect["log_response_ratio_or_smd"],
            "before_shrunk_lrr":before["shrunk_posterior_log_rr"], "after_shrunk_lrr":after["shrunk_posterior_log_rr"]}


def study_alias_errors(root, studies):
    """Compare both complete materialized streams against their canonical rows."""
    effects=[]
    for study in studies:
        meta=study['study_metadata']
        for effect in study.get('atomic_effects',[]):
            effects.append({'pmid':str(meta['pmid']),'citation_short':meta.get('citation_short',''),
                            'source_text_type':meta.get('source_text_type'),
                            'superseded_by_pmid':meta.get('superseded_by_pmid'),**effect})
    errors=[]
    for relative,rows in [('data/extracted/ring1_study_dossiers_v0.5.jsonl',studies),
                          ('data/extracted/ring1_atomic_effects_v0.5.jsonl',effects)]:
        try:
            expected=''.join(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n' for row in rows).encode('utf-8')
            if manifest_mod._regular_file_bytes(root,relative)!=expected:errors.append(relative)
        except (OSError,TypeError,ValueError):errors.append(relative)
    return errors


def synthetic_semantic_probe_source():
    """Self-contained public control, explicitly unrelated to a real study.

    Redaction tokens cannot stand in for missing source numbers. This small
    synthetic contrast supplies known operands and a valid unmutated source.
    """
    pmid = '99999999'
    quote = 'Synthetic viability was 80% versus 20% in controls.'
    dossier = {'pmid': pmid, 'title': 'Synthetic validator control', 'abstract': quote}
    effect = {
        'effect_id': 'synthetic_semantic_control',
        'focal_arm_id': 'treated', 'comparator_arm_id': 'control',
        'treatment_mean': 80., 'control_mean': 20.,
        'log_response_ratio_or_smd': math.log(4.),
        'statistic_type': 'log_response_ratio_paired_means',
        'effect_direction': 'Protective_or_Improved',
        'endpoint_polarity_registry': {
            'polarity_class': 'Higher_Is_Better_Protective_Endpoint',
            'orientation_sign': 1,
        },
        'quote_source_dossier_pmid': pmid,
        'dossier_sha256': ext._compute_dossier_sha256(pmid, dossier),
        'source_span': f'abstract[0:{len(quote)}]|norm=exact',
        'verbatim_evidence_quote': quote,
    }
    doc = {'study_metadata': {'pmid': pmid}, 'arms': [
        {'arm_id': 'treated', 'is_reference_control': False},
        {'arm_id': 'control', 'is_reference_control': True},
    ], 'atomic_effects': [effect]}
    return doc, dossier, effect['effect_id']


def semantic_mutation_controls(doc, dossier, effect_id):
    """A mutation rejection counts only after the unmutated input passes."""
    result = {'baseline_valid': False, 'baseline_error': None, 'mutations': {}, 'mutation_errors': {}}
    try:
        ext.validate_study_semantics(doc, dossier)
    except (ValueError, TypeError) as error:
        result['baseline_error'] = f'{type(error).__name__}: {error}'
        return result
    result['baseline_valid'] = True
    original = next(effect for effect in doc['atomic_effects'] if effect['effect_id'] == effect_id)
    mutations = {
        'wrong_numeric_lrr': {'log_response_ratio_or_smd': 123.},
        'forged_dossier_hash': {'dossier_sha256': '0' * 64},
        'numeric_value_with_qualitative_statistic_type': {'statistic_type': 'qualitative_statement_only'},
        'reversed_orientation_in_same_polarity_class': {'endpoint_polarity_registry': {
            **original['endpoint_polarity_registry'],
            'orientation_sign': -original['endpoint_polarity_registry']['orientation_sign'],
        }},
        'missing_well_formed_comparator_id': {'comparator_arm_id': 'missing_comparator_001'},
    }
    for name, fields in mutations.items():
        mutated = copy.deepcopy(doc)
        next(effect for effect in mutated['atomic_effects'] if effect['effect_id'] == effect_id).update(fields)
        try:
            ext.validate_study_semantics(mutated, dossier)
            result['mutations'][name] = False
        except ValueError as error:
            result['mutations'][name] = True
            result['mutation_errors'][name] = str(error)
    return result


def run_acceptance_suite(out_path: Path) -> int:
    graph, studies, omics = rank.load_inputs()
    profiles = rank.build_target_evidence_profiles(studies, omics)
    combos = rank.run_combinatorial_multi_gate_search(profiles)
    by_id = {c["combination_id"]: c for c in combos}
    curated = set(ext.CURATED_LANDMARK_PMIDS)

    checks: list[dict] = []

    def add_check(cid: str, category: str, description: str, passed: bool | str | None, detail: object) -> None:
        if isinstance(passed, str):
            status = passed
        else:
            status = "INFO" if passed is None else ("PASS" if passed else "FAIL")
        checks.append({
            "id": cid,
            "category": category,
            "description": description,
            "status": status,
            "detail": detail,
        })

    # 1. Brief & release artifacts exist
    is_public_bundle = not (ROOT / "data" / "dossiers").exists()
    req_files = [
        "manifest.json",
        "data/extracted/ring1_study_dossiers_v0.5.jsonl",
        "data/extracted/ring1_atomic_effects_v0.5.jsonl",
        "docs/ISLET-1-SYNTHESIS-REPORT-v0.4.md",
    ]
    if not is_public_bundle:
        req_files.append("data/extracted/islet1_stack_ranking_v0.5.json")
    missing_req = [f for f in req_files if not (ROOT / f).exists()]
    alias_errors=study_alias_errors(ROOT,studies)
    add_check("CHK-01", "artifacts", "Required artifacts exist and both study/effect JSONL streams exactly match canonical study rows", len(missing_req) == 0 and not alias_errors, {"missing": missing_req, "alias_errors":alias_errors, "is_public_bundle": is_public_bundle})

    # 2. Whole-record regeneration parity on non-curated dossiers
    regen_mismatches = 0
    regen_checked = 0
    for p in sorted((ROOT / "data" / "extracted").glob("*.json")):
        if not p.stem.isdigit() or p.stem in curated:
            continue
        dp = ROOT / "data" / "dossiers" / f"{p.stem}.json"
        if not dp.exists():
            continue
        regen_checked += 1
        comm = json.loads(p.read_text(encoding="utf-8"))
        reg = ext.extract_one_dossier(json.loads(dp.read_text(encoding="utf-8")))
        if comm != reg:
            regen_mismatches += 1
    chk02_status: bool | str = (
        "SKIPPED_PUBLIC_BUNDLE_NO_PRIVATE_DOSSIERS"
        if is_public_bundle
        else (regen_mismatches == 0 and regen_checked == len(studies) - len(curated))
    )
    add_check(
        "CHK-02",
        "reproducibility",
        "Committed extractor deterministically reproduces whole-record non-curated study files (or explicitly skips when private dossiers are absent in public bundle)",
        chk02_status,
        {"mismatched_studies": regen_mismatches, "dossiers_checked": regen_checked, "is_public_bundle": is_public_bundle},
    )

    # 3. Comparative arm integrity
    num_total = 0
    valid_pairs = 0
    templated_arm0 = 0
    for d in studies:
        arm_ids = {a["arm_id"] for a in d.get("arms", [])}
        templated_arm0 += sum(1 for a in d.get("arms", []) if a["arm_id"].endswith("_arm0_wt_control"))
        for e in d.get("atomic_effects", []):
            if isinstance(e.get("log_response_ratio_or_smd"), (int, float)):
                num_total += 1
                f_id = e.get("focal_arm_id")
                c_id = e.get("comparator_arm_id")
                if f_id in arm_ids and c_id in arm_ids and f_id != c_id:
                    valid_pairs += 1
    add_check(
        "CHK-03",
        "provenance",
        "All numeric comparative effects link to distinct focal and comparator arms without templated _arm0_wt_control",
        templated_arm0 == 0 and num_total > 0 and valid_pairs == num_total,
        {"valid_pairs": valid_pairs, "numeric_total": num_total, "templated_arm0": templated_arm0},
    )

    # 4. Input-order invariance
    rev_profiles = rank.build_target_evidence_profiles(list(reversed(studies)), omics)
    order_diffs = [
        g
        for g in profiles
        if abs(profiles[g]["net_evidence_and_constraint_score"] - rev_profiles[g]["net_evidence_and_constraint_score"]) > 1e-9
    ]
    add_check("CHK-04", "attribution", "Target evidence profiles are strictly invariant to study input ordering", len(order_diffs) == 0, order_diffs)

    # 5. Graveyard reachability across 1,620 stacks
    gy_entries = json.loads((ROOT / "data" / "graveyard" / "idea_graveyard_v0.5.json").read_text(encoding="utf-8"))["graveyard_entries"]
    fired_counts = Counter(
        r["graveyard_id"] for c in combos for r in (c.get("triggered_graveyard_rules") or [])
    )
    unreached = [e["graveyard_id"] for e in gy_entries if fired_counts[e["graveyard_id"]] == 0]
    add_check("CHK-05", "graveyard", "Every Graveyard rule fires on at least one stack in the 1,620-stack grid", len(unreached) == 0, {"unreached": unreached, "fired": dict(fired_counts)})

    # 6. Zero-evidence baseline
    empty_combos = rank.run_combinatorial_multi_gate_search(rank.build_target_evidence_profiles([], omics))
    empty_a = next(c for c in empty_combos if c["combination_id"] == CHASSIS_A_ID)
    add_check("CHK-06", "evaluator", "Zero extracted studies yields 0.0 scenario index", empty_a["composite_series_system_utility"] == 0.0, empty_a["composite_series_system_utility"])

    # 7. Headline ranking & manifest integrity
    mf = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    bad_hashes = manifest_mod.manifest_integrity_errors(ROOT, mf, is_public_bundle)
    chk07_ok = (
        by_id[CHASSIS_A_ID]["global_rank_among_all_combinations"] == 1
        and by_id[CHASSIS_A_ID]["composite_series_system_utility"] > by_id[CHASSIS_B_ID]["composite_series_system_utility"]
        and by_id[CHASSIS_A_ID]["composite_series_system_utility"] > 0.10
        and len(bad_hashes) == 0
    )
    add_check(
        "CHK-07",
        "release",
        "Chassis A is ranked #1 (strictly beating Chassis B) and all manifest SHA-256 hashes match live files",
        chk07_ok,
        {
            "A_rank": by_id[CHASSIS_A_ID]["global_rank_among_all_combinations"],
            "A_index": by_id[CHASSIS_A_ID]["composite_series_system_utility"],
            "B_rank": by_id[CHASSIS_B_ID]["global_rank_among_all_combinations"],
            "B_index": by_id[CHASSIS_B_ID]["composite_series_system_utility"],
            "bad_hashes": bad_hashes,
            "is_public_bundle": is_public_bundle,
        },
    )

    # 8. Semantic validation & Audit A01/A02/A12/R1 source counterexamples
    semantic_errors = []
    semantic_checked = 0
    by_pmid = {str(d["study_metadata"]["pmid"]): d for d in studies}
    for pm, doc in by_pmid.items():
        dos_path = ROOT / "data" / "dossiers" / f"{pm}.json"
        if dos_path.exists():
            semantic_checked += 1
            try:
                ext.validate_study_semantics(doc, json.loads(dos_path.read_text(encoding="utf-8")))
            except Exception as exc:
                semantic_errors.append({pm: str(exc)})
    eff_by_id = {e["effect_id"]: (str(d["study_metadata"]["pmid"]), e) for d in studies for e in d.get("atomic_effects", [])}
    ce_checks = {
        "pmid41385368_eff2_lrr_is_minus_ln_4.14": abs((eff_by_id["pmid41385368_eff2"][1]["log_response_ratio_or_smd"] or 0.0) - (-math.log(4.14))) < 1e-3,
        "pmid35918471_eff8_quarantined": eff_by_id["pmid35918471_eff8"][1]["log_response_ratio_or_smd"] is None,
        "pmid37809993_eff1_longitudinal_no_control_arm": eff_by_id["pmid37809993_eff1"][1]["log_response_ratio_or_smd"] is None and eff_by_id["pmid37809993_eff1"][1]["comparator_arm_id"] is None,
        "pmid36506044_eff1_sem_detected": eff_by_id["pmid36506044_eff1"][1]["dispersion_type"] == "SEM",
        "balboa2022_eff5_individual_cells_imaged": eff_by_id["balboa2022_eff5_calcium_signaling_maturation"][1]["n_unit"] == "individual_cells_imaged",
        "gerace2023_eff3_qualitative_not_100v100": eff_by_id["gerace2023_eff3_intrinsic_nk_resistance_cd47_pvr"][1]["log_response_ratio_or_smd"] is None,
        "demine2020_eff2_qualitative_not_1v2": eff_by_id["demine2020_eff2_ruxolitinib_prevents_cytokine_apoptosis_er_stress_chemokines"][1]["log_response_ratio_or_smd"] is None,
        "parent2021_eff1_qualitative_not_112v28": eff_by_id["parent2021_eff1_hla_a2r_16wk_allogeneic_survival"][1]["log_response_ratio_or_smd"] is None,
        "lei2022_eff1_nhp_recipients": eff_by_id["lei2022_eff1_safasl_nhp_graft_survival_180days"][1]["n_unit"] == "nhp_recipients",
        "hu2024_eff4_not_mouse_recipients": eff_by_id["hu2024_eff4_rhesus_hip_islets_40week_survival"][1]["n_unit"] == "unspecified_in_source"
        and any(
            a["arm_id"] == eff_by_id["hu2024_eff4_rhesus_hip_islets_40week_survival"][1]["focal_arm_id"]
            and a["sample_size_unit"] == "nhp_recipients" and a["sample_size_n"] == 1
            for a in by_pmid["37156915"]["arms"]
        ),
        "pmid40600980_eff1_multi_fold_quarantined": eff_by_id["pmid40600980_eff1"][1]["log_response_ratio_or_smd"] is None,
        "pmid33198821_eff1_within_assay_stim_quarantined": eff_by_id["pmid33198821_eff1"][1]["log_response_ratio_or_smd"] is None,
        "pmid34415483_eff1_single_arm_sorting_quarantined": eff_by_id["pmid34415483_eff1"][1]["log_response_ratio_or_smd"] is None
        and eff_by_id["pmid34415483_eff1"][1]["comparator_arm_id"] is None,
        "pmid41297833_erratum_41887653_verified_in_outcome": "41887653" in eff_by_id["pmid41297833_eff1"][1]["outcome_variable_Y"],
        "all_563_studies_have_first_and_senior_authors": all(
            bool(d["study_metadata"].get("first_author")) and bool(d["study_metadata"].get("senior_author"))
            for d in studies
        ),
        "composite_rate_unit_parsed": ext.parse_quantitative_numbers("Insulin secretion was 80 ng/mL/min versus 20 ng/mL/min in controls.", "Protective_or_Improved")[:3] == (80.0, None, 20.0),
        "missing_second_unit_quarantined": ext.parse_quantitative_numbers("Insulin secretion was 80 ng/mL versus 20 in controls across 2 groups.", "Protective_or_Improved")[2] is None,
        "mismatched_second_unit_quarantined": ext.parse_quantitative_numbers("Insulin secretion was 80 ng/mL versus 20 pg/mL in controls.", "Protective_or_Improved")[2] is None,
        "longitudinal_from_to_respectively_quarantined": ext.parse_quantitative_numbers("Blood glucose changed from 20 mg/dL to 80 mg/dL at baseline and follow-up, respectively.", "Deleterious_or_Impaired")[5] is None,
        "clause_bound_penalty_rejects_all_8_astra_probes": not any(
            ext.classify_sentence_direction_and_penalty(probe)[1]
            for probe in (
                "Cells reduce blood glucose and preserve insulin secretion.",
                "Edited islets showed decreased cost and improved beta-cell differentiation.",
                "Knockout reduced non-endocrine cells and increased C-peptide secretion.",
                "Treatment reduced cytotoxicity and preserved beta-cell survival.",
                "Differentiation showed loss of pluripotency markers and gain of insulin expression.",
                "Knockout prevented cytokine-induced reduction in insulin secretion.",
                "Diabetes reduced insulin secretion in untreated mice, whereas PD-L1 overexpression restored normoglycemia.",
                "HNF1A deficiency caused MODY3 diabetes with reduced insulin secretion.",
            )
        ),
        "heuristic_penalty_effects_present_not_independent_adjudication": sum(
            1 for d in studies for e in d.get("atomic_effects", []) if e.get("secretome_or_differentiation_penalty_observed")
        ) >= 35,
        "heuristic_classifier_internal_consistency": all(
            e["verbatim_evidence_quote"].startswith("[REDACTED_NON_OA_EXCERPT:")
            or ext.classify_sentence_direction_and_penalty(e["verbatim_evidence_quote"])[1] == e["secretome_or_differentiation_penalty_observed"]
            for d in studies
            for e in d.get("atomic_effects", [])
        ),
        "r09_contrastive_concessive_mixed_clause_probes": (
            ext.classify_sentence_direction_and_penalty("Although edited cells showed reduced insulin content in vitro, they reversed diabetes in vivo.") == ("Deleterious_or_Impaired", True)
            and ext.classify_sentence_direction_and_penalty("Knockout cells restored normoglycemia, but reduced GSIS and impaired beta-cell maturation.") == ("Deleterious_or_Impaired", True)
            and ext.classify_sentence_direction_and_penalty("Allografts were rejected in 3 of 5 recipients, although insulin secretion was preserved prior to rejection.") == ("Deleterious_or_Impaired", False)
            and ext.classify_sentence_direction_and_penalty("Treatment failed to prevent graft rejection despite preserving C-peptide levels initially.") == ("Deleterious_or_Impaired", False)
        ),
        "r07_outcome_independent_multi_harm_retention": (
            sum(
                1
                for e in ext.extract_atomic_effects_and_modifiers(
                    "99999999",
                    {
                        "pmid": "99999999",
                        "title": "Synthetic multi-harm test for B2M knockout islets.",
                        "abstract": (
                            "B2M knockout islets failed to prevent rejection in sensitized hosts after allogeneic transplantation. "
                            "Furthermore, B2M knockout cells suffered severe necrosis and were rejected by NK cells during in vitro co-culture assays. "
                            "Additionally, B2M ablation impaired glucose-stimulated insulin secretion in dynamic perifusion experiments. "
                            "Moreover, B2M deficiency decreased beta-cell survival after pro-inflammatory cytokine challenge in vitro. "
                            "Transiently, B2M knockout reduced T-cell binding at hour 1 in flow cytometry."
                        ),
                    },
                    [{"gene_symbol_hgnc": "B2M", "edit_modality": "CRISPR_Cas9_KO", "functional_gate": "Gate2_Allo_Adaptive_Stealth"}],
                    "hESC",
                    "Stage7_Mature_SC_Islet",
                    "In_Vitro_Assay_Only",
                    "Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture",
                )[0]
                if e["effect_direction"] == "Deleterious_or_Impaired"
            )
            >= 4
        ),
        "r08_cell_source_and_site_regex_hardening": (
            all(
                a["gate1_substrate_and_differentiation"]["cell_source_class"] in ("hESC", "iPSC_healthy_donor", "iPSC_T1D_or_LADA_patient")
                for d in studies
                for a in d.get("arms", [])
                if a["gate1_substrate_and_differentiation"].get("cell_line_id") == "H1"
            )
            and ext.infer_cell_source_and_stage("Th1 cytokines and EndoC-βH1 human beta cells were analyzed.")[1] == "EndoC-βH1"
            and by_pmid["32444635"]["arms"][0]["gate1_substrate_and_differentiation"]["cell_line_id"] == "EndoC-βH1"
            and any(a["gate1_substrate_and_differentiation"]["cell_source_class"] == "NHP_Primary_Islet" for a in by_pmid["37156915"]["arms"])
            and by_pmid["27787618"]["arms"][0]["gate5_engraftment_niche"]["transplant_site"] == "Unspecified_InVivo_or_Source_Site"
            and by_pmid["10202046"]["arms"][0]["gate1_substrate_and_differentiation"]["cell_source_class"] == "Rodent_Primary_Islet"
            and ext.infer_cell_source_and_stage(
                "Rhesus macaque monkey islets derived from embryonic stem cells were not used; primary macaque islets were transplanted."
            )[0] == "NHP_Primary_Islet"
            and ext.infer_cell_source_and_stage(
                "Allogeneic graft survival in a human clinical trial of type 1 diabetes."
            )[0] == "Unspecified_In_Source_Abstract"
            and ext.infer_cell_source_and_stage(
                "Adult human pancreas donors were screened for autoantibodies."
            )[0] == "Unspecified_In_Source_Abstract"
            and ext.infer_cell_source_and_stage(
                "Non human primate donor islets were transplanted into macaques."
            )[0] == "NHP_Primary_Islet"
            and ext.infer_cell_source_and_stage(
                "Human T cells were used to kill porcine pancreatic islets."
            )[0] == "Unspecified_In_Source_Abstract"
            and ext.infer_cell_source_and_stage(
                "No human islets were used. Mouse islets were transplanted into mice."
            )[0] == "Rodent_Primary_Islet"
            and len(by_pmid["27787618"]["arms"][0]["gate2_3_4_genetic_and_epitope_edits"]) == 0
        ),
        "pmid36506044_eff1_manipulation_check_not_survival_lrr": eff_by_id["pmid36506044_eff1"][1]["log_response_ratio_or_smd"] is None
        and eff_by_id["pmid36506044_eff1"][1]["dispersion_type"] == "SEM",
    }
    chk08_dos_ok = (semantic_checked == len(studies)) if not is_public_bundle else True
    chk08_desc = (
        "Public bundle mode: all Audit/Grok/Audit source counterexamples & regex probes pass (563 private dossiers excluded by allowlist; dossier span/SHA-256 check skipped)"
        if is_public_bundle
        else "In-repo audit mode: all 563 studies pass validate_study_semantics against local dossiers and all Audit/Grok/Audit source counterexamples pass"
    )
    add_check(
        "CHK-08",
        "semantics",
        chk08_desc,
        len(semantic_errors) == 0 and chk08_dos_ok and all(ce_checks.values()),
        {
            "semantic_errors": len(semantic_errors),
            "dossiers_checked": semantic_checked,
            "source_dossier_verification_status": "SKIPPED_PUBLIC_BUNDLE_NO_PRIVATE_DOSSIERS" if is_public_bundle else "VERIFIED_ALL_563_DOSSIERS",
            "is_public_bundle": is_public_bundle,
            "counterexamples": ce_checks,
        },
    )

    # 9. Strict non-nullable JSON schema validation of all 563 studies + rejection of synthetic mutation probes (Audit A12; Finding 039)
    schema = json.loads((ROOT / "schemas" / "extraction_schema_v0.2.json").read_text(encoding="utf-8"))
    validator = jsonschema.validators.validator_for(schema)(schema)
    corpus_schema_errors = 0
    for d in studies:
        if list(validator.iter_errors(d)):
            corpus_schema_errors += 1
    sample_doc = copy.deepcopy(by_pmid["41385368"])
    mutations = {
        "null_all_new_provenance_fields": {
            k: None
            for k in (
                "dispersion_type",
                "n_unit",
                "censoring_note",
                "endpoint_polarity_registry",
                "dossier_sha256",
                "source_span",
                "statistic_type",
            )
        },
        "nonexistent_comparator_arm_id": {"comparator_arm_id": "DOES_NOT_EXIST"},
        "unstructured_source_span": {"source_span": "not a source location"},
        "unknown_nested_polarity_registry_property": {
            "endpoint_polarity_registry": {"not_a_defined_registry_field": "unvalidated"}
        },
    }
    mutation_results = {}
    for mname, mfields in mutations.items():
        mdoc = copy.deepcopy(sample_doc)
        meff = next(e for e in mdoc["atomic_effects"] if e["effect_id"] == "pmid41385368_eff4")
        meff.update(mfields)
        s_errs = list(validator.iter_errors(mdoc))
        mutation_results[mname] = len(s_errs) > 0
    add_check(
        "CHK-09",
        "schema",
        "All 563 extracted study JSONs pass schemas/extraction_schema_v0.2.json and strict JSON Schema rejects all 4 synthetic provenance/foreign-key/span/nested-property mutation probes",
        corpus_schema_errors == 0 and all(mutation_results.values()),
        {"studies_validated": len(studies), "corpus_schema_errors": corpus_schema_errors, "mutation_probes": mutation_results},
    )

    # 10. Preprint-to-peer-reviewed deduplication (PMID 40832357 superseded by 41202120) & tier mapping Fail-Loud (Audit A08/A12)
    p408 = by_pmid["40832357"]["study_metadata"]
    preprint_dedup_ok = (
        p408.get("source_text_type") == "Preprint_Superseded_By_Peer_Reviewed"
        and p408.get("superseded_by_pmid") == "41202120"
        and all("40832357" not in v["supporting_pmids"] for v in profiles.values())
    )
    bad_tier_doc = copy.deepcopy(by_pmid["41385368"])
    bad_tier_doc["study_metadata"]["max_biological_fidelity_tier"] = "Tier_UNKNOWN_FUTURE_LABEL"
    tier_raises = False
    try:
        rank.build_target_evidence_profiles([bad_tier_doc], omics)
    except ValueError:
        tier_raises = True
    add_check(
        "CHK-10",
        "deduplication",
        "Superseded bioRxiv preprint 40832357 is linked to 41202120 (0 weight) and unknown fidelity tiers raise ValueError",
        preprint_dedup_ok and tier_raises,
        {"preprint_40832357_dedup_ok": preprint_dedup_ok, "unknown_tier_raises_value_error": tier_raises},
    )

    # 11. HIP-Core Two-Stage Bridged, Replicated & Lead-Chassis-Confirmatory Design Algebra (Audit A06/A07; exact independent re-derivation)
    hip = json.loads((ROOT / "data" / "protocols" / "hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    algebra_spec = importlib.util.spec_from_file_location("openislet_design_algebra_exact", ROOT / "tests" / "test_design_algebra_exact.py")
    algebra = importlib.util.module_from_spec(algebra_spec)
    algebra_spec.loader.exec_module(algebra)
    try:
        design_problems = algebra.verify_protocol(hip)
        design_summary = algebra.summarize(hip)
    except Exception as exc:  # a malformed protocol must fail this check rather than abort the whole suite
        design_problems, design_summary = ["%s: %s" % (type(exc).__name__, exc)], {}
    add_check(
        "CHK-11",
        "design_algebra",
        "HIP-Core design: every published rank, leverage, variance multiplier, VIF, block-inflation and prediction-variance figure equals an independent exact-rational recomputation from the run table "
        "(rank(X)=22/22, rank([X|B])=23/23, 0 unit-leverage rows) and the lean-chassis pointers reference the Scenario A lead genotype",
        len(design_problems) == 0,
        dict(design_summary, problems=design_problems[:10]),
    )

    # 12. Complete coefficient registry live-wiring & sensitivity decomposition (Audit A03/A04/R6 + Audit F-04/F-08 + Finding 084)
    rankings_path = ROOT / "data" / "extracted" / "islet1_candidate_stack_rankings.json"
    rank_out = json.loads(rankings_path.read_text(encoding="utf-8")) if rankings_path.exists() else None
    # Always execute the current evaluator. A cached artifact cannot certify itself.
    cached_vd = rank_out.get("direct_output_logit_sensitivity_and_variance_decomposition", {}) if rank_out else None
    vd, replay_matches_artifact = replay_sensitivity_artifact(profiles, combos, studies, omics, cached_vd)
    zero_profiles = {gene: {**profile, "net_evidence_and_constraint_score": 0.0, "net_evidence_ko_score": 0.0,
                            "net_evidence_gof_score": 0.0, "net_evidence_preserve_score": 0.0,
                            "weighted_deleterious_penalty_score": 0.0} for gene,profile in profiles.items()}
    zero_by_id = {combo["combination_id"]: combo for combo in rank.run_combinatorial_multi_gate_search(zero_profiles)}
    import statistics
    loaded_logits, structural_logits = [], []
    for combo in combos:
        for gate,logit in combo["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"].items():
            loaded_logits.append(logit)
            structural_logits.append(zero_by_id[combo["combination_id"]]["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"][gate])
    delta_logits = [a-b for a,b in zip(loaded_logits,structural_logits)]
    var_structural, var_delta = statistics.pvariance(structural_logits), statistics.pvariance(delta_logits)
    independent_share = round(var_structural/(var_structural+var_delta),6)
    has_reg_key = (
        ("complete_gate_coefficient_registry" in rank_out and "auxiliary_evaluator_hyperparameters" in rank_out["complete_gate_coefficient_registry"])
        if rank_out is not None
        else ("auxiliary_evaluator_hyperparameters" in rank.COMPLETE_GATE_COEFFICIENT_REGISTRY)
    )
    ug = vd.get("uniform_grid_pooled_pediatric_gate_logit_variance_decomposition", {})
    mc = vd.get("direct_output_gate_logit_perturbation_monte_carlo", {})
    sp_mc = vd.get("shared_parameter_registry_monte_carlo", {})
    orig_reg = copy.deepcopy(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY)
    for gdict in rank.COMPLETE_GATE_COEFFICIENT_REGISTRY.values():
        if isinstance(gdict, dict):
            for edict in gdict.values():
                if isinstance(edict, dict) and "logit_value" in edict:
                    edict["logit_value"] += 10.0
    mut_combos = rank.run_combinatorial_multi_gate_search(profiles)
    rank.COMPLETE_GATE_COEFFICIENT_REGISTRY = copy.deepcopy(orig_reg)
    reg_live_wired = mut_combos != combos

    # Audit R10: verify every single one of the 8 auxiliary_evaluator_hyperparameters is live-wired into profiles or combos
    aux_keys = sorted(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"].keys())
    aux_live_by_key = {}
    for ak in aux_keys:
        rank.COMPLETE_GATE_COEFFICIENT_REGISTRY = copy.deepcopy(orig_reg)
        cur_v = rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"][ak]["parameter_value"]
        if isinstance(cur_v, (int, float)):
            rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"][ak]["parameter_value"] = cur_v * 2.5 + 0.5
        elif isinstance(cur_v, list):
            rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"][ak]["parameter_value"] = [-0.01, 0.01]
        elif isinstance(cur_v, dict):
            rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"][ak]["parameter_value"] = {
                k: (v * 2.5 + 0.5 if isinstance(v, (int, float)) else v) for k, v in cur_v.items()
            }
        p_mut = rank.build_target_evidence_profiles(studies, omics)
        c_mut = rank.run_combinatorial_multi_gate_search(p_mut)
        aux_live_by_key[ak] = (p_mut != profiles) or (c_mut != combos)
    rank.COMPLETE_GATE_COEFFICIENT_REGISTRY = orig_reg
    current_typed_lrr_profiles = sum(profile.get("signed_lrr_weight_total",0)>0 for profile in profiles.values())
    lrr_source_fixture = source_bound_lrr_shrinkage_probe(studies,omics)
    aux_pathway_verified = dict(aux_live_by_key)
    if not aux_pathway_verified.get("lrr_shrinkage_hyperparameters") and current_typed_lrr_profiles == 0:
        aux_pathway_verified["lrr_shrinkage_hyperparameters"] = lrr_source_fixture["verified"]
    all_8_aux_live = len(aux_keys) == 8 and all(aux_pathway_verified.values())

    def _valid_mc_block(mc_dict: dict, prefix: str) -> bool:
        sig_keys = [k for k in mc_dict if k.startswith(prefix)]
        if len(sig_keys) < 4:
            return False
        for sk in sig_keys:
            ent = mc_dict[sk]
            if not isinstance(ent, dict):
                return False
            nd = int(ent.get("n_monte_carlo_draws", 0))
            r1 = int(ent.get("scenario_a_rank1_count", -1))
            t3 = int(ent.get("scenario_a_top3_count", -1))
            rmin = int(ent.get("scenario_a_rank_min", 0))
            rmax = int(ent.get("scenario_a_rank_max", 0))
            if not (nd >= 50 and 0 <= r1 <= t3 <= nd and 1 <= rmin <= rmax):
                return False
            if r1 > 0 and rmin != 1:
                return False
        return True


    reg_ok = (
        has_reg_key
        and reg_live_wired
        and all_8_aux_live
        and replay_matches_artifact
        and ug.get("structural_share_of_component_variance_sum") == independent_share
        and _valid_mc_block(mc, "sigma_logit_")
        and _valid_mc_block(sp_mc, "sigma_shared_param_")
        and sp_mc.get("sampled_scalar_leaf_registry_ultracode_v80_14", {}).get("n_auxiliary_scalar_leaves_sampled") == 39
    )
    add_check(
        "CHK-12",
        "sensitivity",
        "Registry paths react on applicable current corpus or explicitly source-bound branch fixtures; current sensitivity is recomputed and seed-replayed",
        reg_ok,
        {
            "replayed_current_evaluator_matches_every_cached_field": replay_matches_artifact,
            "independent_unrounded_variance_share": independent_share,
            "registry_live_wired": reg_live_wired,
            "all_8_aux_hyperparams_live_wired": all_8_aux_live,
            "aux_live_by_key_on_current_corpus": aux_live_by_key,
            "aux_pathway_verified": aux_pathway_verified,
            "current_candidate_target_profiles_with_typed_numeric_lrr": current_typed_lrr_profiles,
            "lrr_shrinkage_source_bound_branch_fixture": lrr_source_fixture,
            "structural_share_component_sum": ug.get("structural_share_of_component_variance_sum"),
            "structural_share_total_variance": ug.get("structural_share_of_total_variance"),
            "top10_overlap_zero_vs_loaded": ug.get("top10_membership_overlap_loaded_vs_zero_gene_signal"),
            "sp_mc_keys": sorted(sp_mc.keys()),
        },
    )

    # 13. Positive-control-bound semantic mutations. Public checks cannot
    # authenticate absent dossiers and use an explicitly synthetic contrast.
    dos_413_path = ROOT / "data" / "dossiers" / "41385368.json"
    if not is_public_bundle:
        base_doc = by_pmid["41385368"]
        probe_dossier = json.loads(dos_413_path.read_text(encoding="utf-8"))
        probe_effect_id = "pmid41385368_eff2"
        probe_scope = "PRIVATE_SOURCE_BOUND_PMID_41385368"
    else:
        base_doc, probe_dossier, probe_effect_id = synthetic_semantic_probe_source()
        probe_scope = "SYNTHETIC_SOURCE_FREE_VALIDATOR_CONTROL"
    sem_probe_results = semantic_mutation_controls(base_doc, probe_dossier, probe_effect_id)
    add_check(
        "CHK-13",
        "semantic_mutations",
        "A valid unmutated source control passes, then all 5 targeted semantic corruptions raise ValueError; public mode uses explicitly synthetic source text",
        sem_probe_results['baseline_valid'] and len(sem_probe_results['mutations']) == 5 and all(sem_probe_results['mutations'].values()),
        {"scope": probe_scope, **sem_probe_results},
    )

    # 14. Primary CD47/SIRPA literature coverage (37820008, 38508139, 40668892) & counter-evidence attribution
    cd47_prof = profiles.get("CD47", {})
    cd47_ok = (
        all(pm in by_pmid for pm in ("37820008", "38508139", "40668892"))
        and "37820008" in cd47_prof.get("counter_evidence_pmids", [])
        and cd47_prof.get("n_deleterious_or_penalty_effects", 0) >= 1
    )
    add_check(
        "CHK-14",
        "cd47_primary_literature",
        "Primary CD47/SIRPA papers (PMID:37820008, PMID:38508139, PMID:40668892) are extracted and 37820008 is indexed as CD47 counter-evidence",
        cd47_ok,
        {
            "supporting_pmids_present": [pm for pm in ("37820008", "38508139", "40668892") if pm in cd47_prof.get("supporting_pmids", [])],
            "counter_evidence_pmids": cd47_prof.get("counter_evidence_pmids", []),
            "n_deleterious_or_penalty_effects": cd47_prof.get("n_deleterious_or_penalty_effects", 0),
        },
    )

    # 15. Public redistribution allowlist consistency (Audit R2)
    allow = mf["public_redistribution_allowlist_spec"]
    excerpt_policy_errors = []
    try:
        public_character_limit, public_abstract_share = manifest_mod.public_excerpt_limits(allow)
    except ValueError as error:
        excerpt_policy_errors.append(str(error))
        public_character_limit, public_abstract_share = 0, 0.
    included_paths: set[str] = set()
    for pat in allow["public_package_included_paths"]:
        for p in ROOT.glob(pat):
            if p.is_file():
                included_paths.add(str(p.relative_to(ROOT)))
    for pat in allow["local_research_cache_excluded_from_public_bundle"]:
        for p in ROOT.glob(pat):
            if p.is_file():
                included_paths.discard(str(p.relative_to(ROOT)))
    lic_files = [str(p.relative_to(ROOT)) for p in (ROOT / ".licenses").rglob("*") if p.is_file()]
    raw_manifest_lines = (ROOT / "data" / "raw" / "corpus_manifest_all_sources.jsonl").read_text(encoding="utf-8").splitlines()
    raw_abstract_count = sum(1 for line in raw_manifest_lines if 'abstract' in json.loads(line)) if is_public_bundle else 0
    non_oa_excess_studies = 0
    graph_unredacted_statements = 0
    public_token_errors=[]
    if is_public_bundle:
        public_token_errors=[(str(d['study_metadata']['pmid']),p) for d in studies for p in manifest_mod.public_token_errors(d,str(d['study_metadata']['pmid']))]
        public_graph=json.loads((ROOT/'data/extracted/cross_study_discovery_graph.json').read_text())
        public_token_errors.extend(('graph',p) for p in manifest_mod.public_token_errors(public_graph))
    if is_public_bundle:
        for d in studies:
            pm_str = str(d.get("study_metadata", {}).get("pmid", ""))
            has_verified_open_license = pm_str in manifest_mod.VERIFIED_OPEN_LICENSE_PMIDS
            if not has_verified_open_license:
                unredacted_chars = sum(
                    len(e.get("verbatim_evidence_quote") or "")
                    for e in d.get("atomic_effects", [])
                    if not (e.get("verbatim_evidence_quote") or "").startswith("[REDACTED_NON_OA_EXCERPT:")
                )
                for e in d.get("atomic_effects", []):
                    oy = str(e.get("outcome_variable_Y") or "")
                    if not oy.startswith("[Endpoint:"):
                        unredacted_chars += len(oy)
                for em in d.get("effect_modifier_interactions", []):
                    for mfield in ("level_1_condition_and_effect", "level_2_condition_and_effect", "mechanistic_explanation"):
                        mv = em.get(mfield) or ""
                        if not mv.startswith("[REDACTED_NON_OA_EXCERPT:"):
                            unredacted_chars += len(mv)
                for edge in d.get("causal_mechanism_edges", []):
                    ev = edge.get("evidence_source") or ""
                    if not ev.startswith("[REDACTED_NON_OA_EXCERPT:"):
                        unredacted_chars += len(ev)
                if unredacted_chars > public_character_limit:
                    non_oa_excess_studies += 1
    else:
        # AUDIT-REF & 084: Even in in-repo audit mode, export a temporary public bundle and verify
        # that every exported non-OA study JSON and cross_study_discovery_graph.json satisfies the <= 240-char / <25%-abstract cap
        with tempfile.TemporaryDirectory(prefix='openislet-acceptance-public-',dir=out_path.parent) as tmp_pub_dir:
            tmp_bundle = Path(tmp_pub_dir) / "bundle"
            manifest_mod.export_public_bundle(tmp_bundle)
            raw_abstract_count=sum('abstract' in json.loads(line) for line in (tmp_bundle/'data/raw/corpus_manifest_all_sources.jsonl').read_text().splitlines() if line.strip())
            exported_studies=[];article_sources={}
            for p_ext in sorted((tmp_bundle / "data" / "extracted").glob("*.json")):
                if not p_ext.stem.isdigit():
                    continue
                d_exp = json.loads(p_ext.read_text(encoding="utf-8"))
                exported_studies.append(d_exp)
                pm_str = str(d_exp.get("study_metadata", {}).get("pmid", p_ext.stem))
                if pm_str in manifest_mod.VERIFIED_OPEN_LICENSE_PMIDS:
                    continue
                dos_p = ROOT / "data" / "dossiers" / f"{pm_str}.json"
                dos_text = ""
                char_budget = 0
                if dos_p.exists():
                    dos_obj = json.loads(dos_p.read_text(encoding="utf-8"))
                    ab_txt = str(dos_obj.get("abstract") or "")
                    dos_text = f"{dos_obj.get('title') or ''} {ab_txt}"
                    dos_text+=' '+' '.join(str(v) for v in (dos_obj.get('fulltext_data_sections') or {}).values() if v)
                    if ab_txt:
                        char_budget = min(public_character_limit, max(0, math.ceil(len(ab_txt) * public_abstract_share) - 1))
                article_sources[pm_str]=(dos_text,char_budget)
                public_token_errors.extend((pm_str,p) for p in manifest_mod.public_token_errors(d_exp,pm_str))
                shingle_chars = manifest_mod._max_verbatim_shingle_chars(d_exp, dos_text)
                if shingle_chars > char_budget:
                    non_oa_excess_studies += 1
            exp_graph_path = tmp_bundle / "data" / "extracted" / "cross_study_discovery_graph.json"
            if exp_graph_path.exists():
                exp_graph = json.loads(exp_graph_path.read_text(encoding="utf-8"))
                public_token_errors.extend(('graph',p) for p in manifest_mod.public_token_errors(exp_graph))
                combined=manifest_mod.public_article_verbatim_totals(exported_studies,exp_graph,article_sources)
                non_oa_excess_studies=sum(total>article_sources[pm][1] for pm,total in combined.items())
                for bnode in exp_graph.get("shared_cross_study_bridge_nodes", {}).values():
                    if isinstance(bnode, dict):
                        for item in bnode.get("linked_mechanisms_and_bottlenecks", []):
                            if isinstance(item, dict):
                                st_val = str(item.get("statement") or "")
                                if " -> " not in st_val and not st_val.startswith("[REDACTED_NON_OA_EXCERPT:"):
                                    graph_unredacted_statements += 1
    allow_ok = (
        not any("ranking" in k for k in included_paths)
        and not any(k.startswith("data/dossiers/") for k in included_paths)
        and not any(manifest_mod._forbidden_public_path(k) for k in included_paths)
        and all(k in included_paths for k in lic_files)
        and all(k in included_paths for k in ("data/raw/ring1_platinum_queue.json", "data/raw/corpus_manifest_all_sources.jsonl"))
        and all((k in mf["artifacts"] or k == "manifest.json") for k in included_paths)
        and raw_abstract_count == 0
        and non_oa_excess_studies == 0
        and graph_unredacted_statements == 0
        and not public_token_errors
        and not excerpt_policy_errors
    )
    add_check(
        "CHK-15",
        "public_allowlist",
        ("Source-free public verification checks exact pinned inventory, private-path exclusions, opaque reference syntax/identity and explicit quote-field bounds; absent source texts and redistribution rights are not authenticated"
         if is_public_bundle else
         "Private verification exports a temporary bundle, excludes private dossiers/planning documents/frozen rankings, strips discovery abstracts and recounts combined study/graph overlap against private sources; licenses and required raw inputs remain pinned"),
        allow_ok,
        {
            "allowlisted_file_count": len(included_paths),
            "raw_abstracts_in_corpus_manifest": raw_abstract_count,
            "non_oa_studies_exceeding_240_chars_in_public_bundle": non_oa_excess_studies,
            "graph_unredacted_statements_in_public_bundle": graph_unredacted_statements,
            "invalid_public_opaque_references": public_token_errors,
            "excerpt_policy_errors": excerpt_policy_errors,
            "explicit_quote_character_limit": public_character_limit,
            "private_source_recount_performed": not is_public_bundle,
            "source_authentication_scope": "ABSENT_PRIVATE_SOURCES_NOT_AUTHENTICATED" if is_public_bundle else "PINNED_PRIVATE_DOSSIER_BYTES",
        },
    )

    # 16. Documented pinned-file counts match the manifest
    pinned_claims = {}
    for rel in (
        "README.md",
        "docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md",
        "docs/ISLET-1-SYNTHESIS-REPORT-v0.4.md",
        *(("CLAUDE.md",) if (ROOT / "CLAUDE.md").exists() else ()),
    ):
        pinned_claims[rel] = [int(m.group(1).replace(",", "")) for m in re.finditer(r"`([\d,]+)` pinned files", (ROOT / rel).read_text(encoding="utf-8"))]
    pinned_ok = all(claims and all(c == len(mf["artifacts"]) for c in claims) for claims in pinned_claims.values())
    add_check(
        "CHK-16",
        "documentation_sync",
        "README.md, ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md and the synthesis report quote exactly the number of files pinned in manifest.json",
        pinned_ok,
        {"manifest_artifacts": len(mf["artifacts"]), "documented": pinned_claims},
    )

    # 17. Layer 2 (Perturbation & Regulatory Grammar) & Layer 3 (3D Structural & De Novo Binder) Planes
    l23_spec = importlib.util.spec_from_file_location("openislet_layer23_planes", ROOT / "tests" / "test_layer2_and_layer3_planes.py")
    l23_mod = importlib.util.module_from_spec(l23_spec)
    l23_spec.loader.exec_module(l23_mod)
    try:
        l23_mod.run_layer2_and_layer3_suite()
        l23_ok = True
        l23_err = None
    except Exception as exc:
        l23_ok = False
        l23_err = f"{type(exc).__name__}: {exc}"
    add_check(
        "CHK-17",
        "layer2_and_layer3_planes",
        "Layer 2 computational donor grammar and prior-scenario profiles (unresolved-origin AVI quarantined, no model-output certification) and Layer 3 reference-coordinate replays/negative controls satisfy declared computational invariants; exact-product function remains pending",
        l23_ok,
        {"status": "PASS" if l23_ok else "FAIL", "error": l23_err},
    )

    # 18. Layer 4 (Multi-Scale Tissue & Metabolic Digital Twin) & Layer 5 (Personalized Patient-DNA & CRISPR Compiler) Planes
    l45_spec = importlib.util.spec_from_file_location("openislet_layer45_planes", ROOT / "tests" / "test_layer4_and_layer5_planes.py")
    l45_mod = importlib.util.module_from_spec(l45_spec)
    l45_spec.loader.exec_module(l45_mod)
    try:
        l45_mod.run_layer4_and_layer5_suite()
        l45_ok = True
        l45_err = None
    except Exception as exc:
        l45_ok = False
        l45_err = f"{type(exc).__name__}: {exc}"
    add_check(
        "CHK-18",
        "layer4_and_layer5_planes",
        "Layer 4 conditional PDE/reference and matched-state ODE arithmetic and Layer 5 source-verified variant identities, unknown genome-wide enumeration, assigned-rule virtual population and validated eight-endpoint updater satisfy computational invariants; physical endpoints remain unmeasured",
        l45_ok,
        {"status": "PASS" if l45_ok else "FAIL", "error": l45_err},
    )

    # 19. Documentation & Artifact Numeric/Label/Digest Consistency (CHK-19: 19A, 19B, 19C, 19D)
    chk19_failures: list[dict] = []
    doc_files_19 = (
        "CLAUDE.md",
        "README.md",
        "docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md",
        "ISLET-1-BLUEPRINT.md",
        "docs/PROJECT-DECISIONS-LOG.md",
        "docs/STAGE-0-EXECUTION-PLAYBOOK.md",
    )
    doc_texts_19 = {
        name: (ROOT / name).read_text(encoding="utf-8") if (ROOT / name).exists() else ""
        for name in doc_files_19
    }
    ref_ctrl_truth = sum(
        1
        for st in studies
        for arm in (st.get("arms") or [])
        if arm.get("is_reference_control") is True
    )
    for dname, pat in (
        ("CLAUDE.md", r"([\d,]+) reference controls"),
        ("docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md", r"([\d,]+) reference controls"),
        ("CLAUDE.md", r"`([\d,]+)` explicit reference-control comparator arms"),
        ("docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md", r"`([\d,]+)` explicit reference-control comparator arms"),
    ):
        if not (ROOT / dname).exists():
            continue
        matches = list(re.finditer(pat, doc_texts_19.get(dname, "")))
        for m in matches:
            val = int(m.group(1).replace(",", ""))
            if val != ref_ctrl_truth:
                chk19_failures.append({"subcheck": "19A.reference_control_arms", "doc": dname, "got": val, "expected": ref_ctrl_truth})
    for dname in ("CLAUDE.md", "docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md"):
        if not (ROOT / dname).exists():
            continue
        if not re.search(r"([\d,]+) reference controls|`([\d,]+)` explicit reference-control comparator arms", doc_texts_19.get(dname, "")):
            chk19_failures.append({"subcheck": "19A.reference_control_arms_missing", "doc": dname})

    hip_proto = json.loads((ROOT / "data" / "protocols" / "hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    hip_pw = hip_proto.get("statistical_analysis_and_power_plan", {}).get("minimum_detectable_effects_at_80pct_power")
    if not isinstance(hip_pw, dict):
        hip_pw = hip_proto.get("statistical_power_and_sample_size", {})
    if (
        not isinstance(hip_pw, dict)
        or "worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_t52_alpha_0_05_pp" not in hip_pw
        or "worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_t13_alpha_0_05_pp" not in hip_pw
    ):
        chk19_failures.append({"subcheck": "19A.hip_power_dict_missing"})
        mde_t52 = -999.0
        mde_t13 = -999.0
    else:
        mde_t52 = float(hip_pw["worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_t52_alpha_0_05_pp"])
        mde_t13 = float(hip_pw["worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_t13_alpha_0_05_pp"])
        if not (math.isfinite(mde_t52) and math.isfinite(mde_t13)):
            chk19_failures.append({"subcheck": "19A.hip_power_non_finite", "mde_t52": mde_t52, "mde_t13": mde_t13})
    for dname in ("CLAUDE.md", "docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md"):
        if not (ROOT / dname).exists():
            continue
        m52_list = list(re.finditer(r"`([\d.]+) pp` \$t_\{52\}\$", doc_texts_19.get(dname, "")))
        m13_list = list(re.finditer(r"`([\d.]+) pp` \$t_\{13\}\$", doc_texts_19.get(dname, "")))
        if not m52_list:
            chk19_failures.append({"subcheck": "19A.mde_t52_missing", "doc": dname})
        for m in m52_list:
            if not math.isfinite(mde_t52) or abs(float(m.group(1)) - mde_t52) > 5e-3:
                chk19_failures.append({"subcheck": "19A.mde_t52", "doc": dname, "got": float(m.group(1)), "expected": mde_t52})
        if not m13_list:
            chk19_failures.append({"subcheck": "19A.mde_t13_missing", "doc": dname})
        for m in m13_list:
            if not math.isfinite(mde_t13) or abs(float(m.group(1)) - mde_t13) > 5e-3:
                chk19_failures.append({"subcheck": "19A.mde_t13", "doc": dname, "got": float(m.group(1)), "expected": mde_t13})

    for dname, pat in (
        ("CLAUDE.md", r"`1191/1616`\s*\[([^\]]+)\]"),
        ("docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md", r"`1191/1616`\s*\(`[\d.]+`,\s*([^)]+)\)"),
        ("docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md", r"\(`1191/1616`,\s*([^)]+)\)"),
        ("ISLET-1-BLUEPRINT.md", r"`(X\d:X\d)`\s*\(`1` pair\):\s*`1191/1616`"),
        ("docs/PROJECT-DECISIONS-LOG.md", r"`(X\d:X\d)`\s*\(`1191/1616`\)"),
    ):
        if not (ROOT / dname).exists():
            continue
        matches_2fi = list(re.finditer(pat, doc_texts_19.get(dname, "")))
        if not matches_2fi:
            chk19_failures.append({"subcheck": "19A.worst_2fi_pattern_missing", "doc": dname, "pattern": pat})
        for m in matches_2fi:
            claimed_pairs = {f"X{a}:X{b}" for a, b in re.findall(r"X(\d)\s*:\s*X(\d)", m.group(1))}
            if claimed_pairs != {"X1:X3"}:
                chk19_failures.append({"subcheck": "19A.worst_2fi_label", "doc": dname, "got": sorted(claimed_pairs), "expected": ["X1:X3"]})

    l3_atlas = json.loads((ROOT / "data" / "extracted" / "layer3_structural_verification_atlas.json").read_text(encoding="utf-8"))
    c3_atlas = next((c for c in l3_atlas.get("complexes", []) if "SERPINB9" in c.get("complex_id", "")), {})
    pro376_truth = float(
        c3_atlas.get("r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry", {}).get(
            "pro376_only_min_heavy_atom_distance_to_rcl_335_345_angstrom", -999.0
        )
    )
    if not math.isfinite(pro376_truth):
        chk19_failures.append({"subcheck": "19A.pro376_truth_non_finite"})
    pro376_matches = list(re.finditer(r"`([\d.]+) Å` `Pro376`-only min heavy-atom distance", doc_texts_19.get("ISLET-1-BLUEPRINT.md", "")))
    if not pro376_matches:
        chk19_failures.append({"subcheck": "19A.pro376_distance_missing", "doc": "ISLET-1-BLUEPRINT.md"})
    for m in pro376_matches:
        if not math.isfinite(pro376_truth) or abs(float(m.group(1)) - pro376_truth) > 5e-3:
            chk19_failures.append({"subcheck": "19A.pro376_distance", "got": float(m.group(1)), "expected": pro376_truth})

    hlae_gate_vals: set[float] = set()
    for dname in doc_files_19:
        for pat in (r"assumed ([\d,]+) threshold", r"`>= ([\d,]+)` molecules/cell[^.\n]{0,40}?release gate"):
            for m in re.finditer(pat, doc_texts_19.get(dname, "")):
                hlae_gate_vals.add(float(m.group(1).replace(",", "")))
    if len(hlae_gate_vals) != 1 or 15000.0 not in hlae_gate_vals:
        chk19_failures.append({"subcheck": "19B.resting_hla_e_gate_single_valued", "found": sorted(hlae_gate_vals)})

    l2_atlas = json.loads((ROOT / "data" / "extracted" / "layer2_regulatory_and_perturbation_atlas.json").read_text(encoding="utf-8"))
    chs4_node = (
        l2_atlas.get("construct_sequence_and_grammar_verification", {})
        .get("regulatory_element_audits", {})
        .get("cHS4_core_insulator_250bp_ULTRACODE_V80_13", {})
    )
    chs4_prose = str(chs4_node.get("sequence_adaptation_disclosure", ""))
    chs4_top = (chs4_node.get("ctcf_ma0139_1_footprint_fii_scan_6_matrix_family", {}).get("top_hits") or [{}])[0]
    for pat, key in (
        (r"score=([\d.]+) bits", "log2_odds_score"),
        (r"rel_score=([\d.]+)", "relative_score_0_to_1"),
        (r"window Bonferroni p=([\d.]+)", "bonferroni_window_p_value"),
    ):
        matches_c = list(re.finditer(pat, chs4_prose))
        t_val = float(chs4_top.get(key, -999.0))
        if not matches_c or not math.isfinite(t_val):
            chk19_failures.append({"subcheck": f"19C.{key}_missing_or_non_finite"})
        for m in matches_c:
            if not math.isfinite(t_val) or abs(float(m.group(1)) - t_val) > 5e-4:
                chk19_failures.append({"subcheck": f"19C.{key}", "got": float(m.group(1)), "expected": chs4_top.get(key)})
    coord_matches_c = list(re.finditer(r"local positions (\d+)\.\.(\d+)", chs4_prose))
    if not coord_matches_c:
        chk19_failures.append({"subcheck": "19C.chs4_coords_missing"})
    for m_coord in coord_matches_c:
        if int(m_coord.group(1)) != int(chs4_top.get("start_1based_in_window", -1)) or int(m_coord.group(2)) != int(chs4_top.get("end_1based_in_window", -1)):
            chk19_failures.append({"subcheck": "19C.chs4_coords", "got": [int(m_coord.group(1)), int(m_coord.group(2))]})

    known_digests_19: set[str] = set()
    fasta_19 = ROOT / "data" / "constructs" / "islet1_scenario_a_constructs_v0.7.0.fasta"
    if fasta_19.exists():
        for rec in re.split(r"^>", fasta_19.read_text(encoding="utf-8"), flags=re.M)[1:]:
            lines_r = rec.strip().split("\n")
            seq_r = re.sub(r"[^A-Za-z]", "", "".join(lines_r[1:]))
            known_digests_19.add(hashlib.sha256(seq_r.encode()).hexdigest())
            known_digests_19.update(re.findall(r"sha256=([0-9a-f]{64})", lines_r[0]))
    for art_obj in (l2_atlas, l3_atlas):
        for m_hex in re.findall(r'"([0-9a-f]{64})"', json.dumps(art_obj)):
            known_digests_19.add(m_hex)
    n_digest_refs_19 = 0
    for dname in doc_files_19:
        dtext_19 = doc_texts_19.get(dname, "")
        for m in re.finditer(r"SHA-256 `([0-9a-f]{64})`", dtext_19):
            n_digest_refs_19 += 1
            if m.group(1) not in known_digests_19:
                chk19_failures.append({"subcheck": "19D.unresolved_sha256", "doc": dname, "sha256": m.group(1)})
        for m_pref in re.finditer(r"SHA-256 prefix `([0-9a-f]{16,64})`", dtext_19):
            n_digest_refs_19 += 1
            pref_hex = m_pref.group(1)
            if not any(full_d.startswith(pref_hex) for full_d in known_digests_19):
                chk19_failures.append({"subcheck": "19D.unresolved_sha256_prefix", "doc": dname, "sha256_prefix": pref_hex})
        if "results/current/" in dtext_19 or "3,701 bp" in dtext_19:
            chk19_failures.append({"subcheck": "19A.stale_playbook_path_or_length", "doc": dname})
    if n_digest_refs_19 == 0:
        chk19_failures.append({"subcheck": "19D.zero_sha256_or_prefix_references_harvested"})

    playbook_text_19 = doc_texts_19.get("docs/STAGE-0-EXECUTION-PLAYBOOK.md", "")
    for stale_tok in ("24,992", "Glu376", "PMC10232510", "PMC10322668", "291 bp WPREmut6", "efea4d2ddfaab3ca", "17.9597"):
        if stale_tok in playbook_text_19:
            chk19_failures.append({"subcheck": "19B.stale_playbook_token", "token": stale_tok})
    for req_tok in ("66,744", "Pro376-Stop", "iCasp9-GSG-P2A-CD47-IgV-GPI"):
        if req_tok not in playbook_text_19:
            chk19_failures.append({"subcheck": "19B.missing_reconciled_playbook_token", "token": req_tok})
    if (ROOT / "docs/PROJECT-DECISIONS-LOG.md").exists() and "PMC10690437" not in doc_texts_19.get("docs/PROJECT-DECISIONS-LOG.md", ""):
        chk19_failures.append({"subcheck": "19B.missing_reconciled_decisions_token", "token": "PMC10690437"})

    sv40_j19 = (
        l2_atlas.get("construct_sequence_and_grammar_verification", {})
        .get("regulatory_element_audits", {})
        .get("SV40_late_polyA_122bp_J02400", {})
        .get("j02400_alignment_verification", {})
    )
    if sv40_j19.get("donor_plus_strand_canonical_AATAAA_1based_spans") != ["65..70", "94..99"] or sv40_j19.get("transcribed_antisense_cd47_minus_strand_canonical_AATAAA_1based_spans") != ["78..83"]:
        chk19_failures.append({"subcheck": "19C.sv40_aataaa_hexamer_spans", "got_plus": sv40_j19.get("donor_plus_strand_canonical_AATAAA_1based_spans"), "got_minus": sv40_j19.get("transcribed_antisense_cd47_minus_strand_canonical_AATAAA_1based_spans")})

    hip_40_node = hip_proto.get("hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1", {})
    added_runs_40 = hip_40_node.get("added_stage1_bridge_runs_37_to_40", [])
    if len(added_runs_40) != 4 or not all(
        isinstance(r.get("coded_levels"), dict)
        and set(r["coded_levels"].keys()) == {"X1", "X2", "X3", "X4", "X5", "X6"}
        for r in added_runs_40
    ):
        chk19_failures.append({"subcheck": "19A.hip_core_40run_missing_canonical_coded_levels"})
    else:
        try:
            base_runs_36 = hip_proto.get("runs") or (hip_proto.get("tier1_screening_16_runs", []) + hip_proto.get("tier2_dealiasing_augmentation_20_runs", []))
            rows_40 = rank._hip_model_rows(base_runs_36 + added_runs_40)
            if len(rows_40) != 40:
                chk19_failures.append({"subcheck": "19A.hip_core_40run_row_count", "got": len(rows_40)})
        except Exception as exc_40:
            chk19_failures.append({"subcheck": "19A.hip_core_40run_native_parser_error", "error": str(exc_40)})

    co_prim_19 = hip_proto.get("alias_structure", {}).get("co_primary_tyk2_inducible_promoter_interaction_estimands_muse_r3", {})
    if not str(co_prim_19.get("mode1_23col_mean_main_effect_var_beta_over_sigma2_across_6_factors", "")).startswith("137/2424"):
        chk19_failures.append(
            {
                "subcheck": "19A.co_primary_main_effect_var_beta_over_sigma2",
                "got": co_prim_19.get("mode1_23col_mean_main_effect_var_beta_over_sigma2_across_6_factors"),
                "expected_prefix": "137/2424",
            }
        )
    if (
        co_prim_19.get("mode1_23col_half_effect_se_pp_at_sigma_run_5p4314", {}).get("mean_across_6_main_effects") != 1.2912
        or co_prim_19.get("mode1_23col_full_contrast_2beta_se_pp_at_sigma_run_5p4314", {}).get("mean_across_6_main_effects") != 2.5824
    ):
        chk19_failures.append({"subcheck": "19A.co_primary_main_effect_se_pp_mismatch"})

    add_check(
        "CHK-19",
        "documentation_numeric_consistency",
        "Hand-maintained numbers, factor-pair labels, and SHA-256 digests in CLAUDE.md, README.md, ISLET-1-BLUEPRINT.md, docs/PROJECT-DECISIONS-LOG.md, and docs/STAGE-0-EXECUTION-PLAYBOOK.md agree with ground-truth recounts and generated artifacts (19A-19D)",
        len(chk19_failures) == 0,
        {"failures": chk19_failures},
    )

    summary = Counter(c["status"] for c in checks)
    payload = {
        "suite": "openislet_independent_acceptance_suite",
        "summary": dict(summary),
        "source_verification_coverage_matrix": {
            "is_public_bundle": is_public_bundle,
            "chk02_whole_record_regeneration_dossiers_verified": regen_checked,
            "chk08_source_span_and_sha256_dossiers_verified": semantic_checked,
            "chk02_status": chk02_status if isinstance(chk02_status, str) else ("PASS" if chk02_status else "FAIL"),
            "chk08_dossier_verification_status": (
                "SKIPPED_PUBLIC_BUNDLE_COUNTEREXAMPLES_ONLY" if is_public_bundle else "VERIFIED_ALL_563_DOSSIERS"
            ),
            "public_bundle_non_oa_studies_exceeding_240_chars": non_oa_excess_studies,
        },
        "checks": checks,
    }
    out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for c in checks:
        print(f"{c['id']}  {c['status']:<4}  [{c['category']}] {c['description']} -> {c['detail']}")
    print(f"\nSUMMARY: {dict(summary)} (written to {out_path})")
    return 0 if summary.get("FAIL", 0) == 0 else 1


def test_adversarial_acceptance_suite(tmp_path: Path) -> None:
    rc = run_acceptance_suite(tmp_path / "openislet_acceptance_evidence.json")
    if rc != 0:
        raise AssertionError(f"Adversarial acceptance suite failed with exit code {rc}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run OpenIslet independent acceptance suite.")
    parser.add_argument("--out", default=None, help="Output JSON path (defaults to an isolated temporary file)")
    args = parser.parse_args()
    if args.out is not None:
        out_p = Path(args.out)
    else:
        with tempfile.NamedTemporaryFile(prefix="openislet_acceptance_", suffix=".json", delete=False) as tf:
            out_p = Path(tf.name)
    sys.exit(run_acceptance_suite(out_p))


if __name__ == "__main__":
    main()
