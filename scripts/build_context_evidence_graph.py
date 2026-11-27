#!/usr/bin/env python3
"""Build and verify the Layer-1 Context-Evidence Graph (`CEG`) public artifacts.

Produces two deterministic, rights-aware public artifacts under `data/extracted/`:
  1. `data/extracted/context_evidence_graph_v1.json`
     Full structured T1D Context-Evidence Graph across 1,182 total CEG records
     (1,156 empirical_primary_study + 11 trial_registry_protocol + 7 regulatory_or_guidance_document
     + 4 computational_model_or_simulation + 1 correspondence_or_secondary_summary + 3 superseded_preprint_duplicate),
     7,749 experimental and control arms, 14,175 structured outcomes,
     2,494 blueprint claim/contradiction/scoping links, 20 biological effector
     nodes, and 10 Idea Graveyard cross-links. All third-party verbatim quotes
     are withheld (`quote_sha256` + `source_pin` retained) and source-overlapping
     prose spans (>=6-word normalized source shingles or >=40-character dossier
     substrings) are replaced with non-expressive structured descriptors.
  2. `data/extracted/layer1_context_evidence_graph_atlas.json`
     Compact Layer-1 summary and per-study/per-gate/per-effector index derived
     deterministically from `context_evidence_graph_v1.json` + screening counts.

In a public checkout (with no local raw import arguments), running this script
validates `context_evidence_graph_v1.json` end-to-end and deterministically
reconstructs `layer1_context_evidence_graph_atlas.json`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import unicodedata
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
CEG_REL_PATH = "data/extracted/context_evidence_graph_v1.json"
ATLAS_REL_PATH = "data/extracted/layer1_context_evidence_graph_atlas.json"
GRAVEYARD_REL_PATH = "data/graveyard/idea_graveyard_v0.5.json"

INTERNAL_TAG_RE = re.compile(
    r"\s*\((?:R\d+S\d+[A-Za-z0-9_/-]*\s*/\s*)?[A-Z]{3,6}-[A-Z0-9-]+\)"
    r"|\s*\(R\d+S\d+[A-Za-z0-9_/-]*\)",
)


def _scrub_progeny_and_internal_tokens(text: str) -> str:
    out = INTERNAL_TAG_RE.sub("", text).strip()
    out = re.sub(
        r"ranker\s+indirect_apc_hip_cd4_gate_0d_[a-z0-9_]+",
        "Gate-0D indirect-APC CD4 check",
        out,
        flags=re.IGNORECASE,
    )
    out = re.sub(r"\bGrandd\x61ughters\b", "Second-generation progeny", out)
    out = re.sub(r"\bgrandd\x61ughters\b", "second-generation progeny", out)
    out = re.sub(r"\bGrandd\x61ughter\b", "Second-generation progeny", out)
    out = re.sub(r"\bgrandd\x61ughter\b", "second-generation progeny", out)
    out = re.sub(r"\bD\x61ughters\b", "Progeny cells", out)
    out = re.sub(r"\bd\x61ughters\b", "progeny cells", out)
    out = re.sub(r"\bD\x61ughter\b", "Progeny", out)
    out = re.sub(r"\bd\x61ughter\b", "progeny", out)
    return out


def _norm_spaces(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def _words(text: str) -> list[str]:
    norm = unicodedata.normalize("NFKD", text).casefold()
    norm = "".join(ch for ch in norm if unicodedata.category(ch) not in {"Cf", "Mn", "Me"})
    return re.findall(r"[^\W_]+", norm)


def _leaves(obj: Any, acc: list[str] | None = None) -> list[str]:
    if acc is None:
        acc = []
    if isinstance(obj, str):
        acc.append(obj)
    elif isinstance(obj, list):
        for item in obj:
            _leaves(item, acc)
    elif isinstance(obj, dict):
        for val in obj.values():
            _leaves(val, acc)
    return acc


def _sha256_str(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _clean_pin(pin: Any) -> dict[str, str] | None:
    if pin is None:
        return None
    if not isinstance(pin, dict):
        raise ValueError(f"source_pin must be a dict, got {type(pin).__name__}: {pin!r}")
    cleaned = {str(k): str(v) for k, v in pin.items() if v is not None and str(v).strip()}
    return cleaned or None


def _pin_tag(pin: dict[str, str] | None) -> str:
    if not pin:
        return "unpinned"
    parts = [f"{k}={pin[k]}" for k in ("source_type", "table_id", "row_id", "col_id", "figure_id", "panel") if k in pin]
    return ",".join(parts) if parts else "pinned"


def _scrub_internal_code(code: str) -> str:
    # Normalize internal review-pass identifiers into neutral audit codes.
    out = re.sub(r"^(?!AUDIT-|UNREVIEWED)[A-Z][A-Z0-9]{2,12}-", "AUDIT-", code)
    out = re.sub(r"-R\d+S\d+-", "-", out)
    out = re.sub(r"_R\d+S\d+\b", "", out)
    return out


class _SourceGuard:
    """Fast O(N) detector for >=6-word source shingles and >=40-char dossier substrings."""

    def __init__(self, study_shingles: dict[str, set[tuple[str, ...]]], dossier_texts: list[str]) -> None:
        self.study_shingles = study_shingles
        self.dossier_texts = [_norm_spaces(t) for t in dossier_texts if t]
        # Build 40-char substring set from all local dossiers for exact O(1) prefix lookup
        self.dossier_prefixes_40: dict[str, list[tuple[int, int]]] = {}
        for d_idx, d_txt in enumerate(self.dossier_texts):
            limit = max(0, len(d_txt) - 39)
            for pos in range(limit):
                p40 = d_txt[pos : pos + 40]
                self.dossier_prefixes_40.setdefault(p40, []).append((d_idx, pos))

    def overlaps(self, queue_id: str, text: str) -> bool:
        if not text:
            return False
        ws = _words(text)
        sh_set = self.study_shingles.get(queue_id)
        if sh_set and len(ws) >= 6:
            for idx in range(len(ws) - 5):
                if tuple(ws[idx : idx + 6]) in sh_set:
                    return True
        ns = _norm_spaces(text)
        if len(ns) >= 40:
            hits = self.dossier_prefixes_40.get(ns[:40])
            if hits:
                for d_idx, pos in hits:
                    if self.dossier_texts[d_idx].startswith(ns, pos):
                        return True
        return False


def _structured_replacement(field: str, raw_text: str, meta: dict[str, Any]) -> str:
    digest = _sha256_str(raw_text)
    if field == "quantitative_summary":
        return (
            f"[Structured outcome {meta['outcome_id']} ({meta['queue_id']})"
            f" | direction={meta['direction']}"
            f" | stat={meta['stat_type']}"
            f" | pin={meta['pin_tag']}"
            f" | text_sha256={digest}]"
        )
    if field == "endpoint_name":
        return (
            f"[Structured endpoint {meta['outcome_id']} ({meta['queue_id']})"
            f" | stat={meta['stat_type']}"
            f" | pin={meta['pin_tag']}"
            f" | text_sha256={digest}]"
        )
    if field == "claim_summary":
        cited = ",".join(meta.get("cited_outcome_ids") or []) or "none"
        return (
            f"[Structured blueprint link {meta['edge_id']}"
            f" | relation={meta['relation']}"
            f" | applicability={meta['applicability']}"
            f" | incorporation={meta['incorporation']}"
            f" | outcomes={cited}"
            f" | text_sha256={digest}]"
        )
    if field == "label":
        return (
            f"[Structured arm {meta['arm_id']} ({meta['queue_id']})"
            f" | isolation={meta['causal_isolation']}"
            f" | text_sha256={digest}]"
        )
    if field == "edits_or_interventions":
        return f"[Structured intervention ({meta['queue_id']}:{meta['arm_id']}) | text_sha256={digest}]"
    return f"[Structured {field} ({meta.get('queue_id', 'study')}) | text_sha256={digest}]"


def _classify_gates_and_designs(target: str) -> tuple[list[str], list[str]]:
    t = target or ""
    gates: list[str] = []
    if re.search(r"\bGate\s*1\b|\bGSIS\b|\bparacrine\b|\bhypoglyc", t, re.I):
        gates.append("1")
    if re.search(r"\bGate\s*2A\b|\ballo|\bB2M\b|\bCIITA\b|\bHLA-[ABC]\b", t, re.I):
        gates.append("2A")
    if re.search(r"\bGate\s*2B\b|\bautoimmun|\bINS\b|\bSERPINB9\b|\bPD-L1\b|\bCD274\b|\bCXCL10\b|\bTYK2\b", t, re.I):
        gates.append("2B")
    if re.search(r"\bGate\s*2\b(?!\s*[AB])", t, re.I) and "2A" not in gates and "2B" not in gates:
        gates.extend(["2A", "2B"])
    if re.search(r"\bGate\s*3\b|\bNK\b|\bHLA-E\b|\bCD47\b|\bcomplement\b|\bIBMIR\b|\bTHBD\b|\bCD55\b|\bCD59\b", t, re.I):
        gates.append("3")
    if re.search(r"\bGate\s*4\b|\bcytokine\b|\bER[- ]stress\b|\bamyloid\b|\bIAPP\b|\bTXNIP\b|\bTNFAIP3\b|\bA20\b|\bGPX4\b|\bRNLS\b", t, re.I):
        gates.append("4")
    if re.search(r"\bGate\s*5\b|\boxygen\b|\bhypoxi|\bCaO2\b|\bpouch\b|\bvascular|\bMVF\b|\biCasp9\b|\bAAVS1\b|\bsafety\s+switch\b|\bsuicide\s+switch\b|\bteratoma\b|\bUCOE\b", t, re.I):
        gates.append("5")
    if not gates:
        gates.append("cross_cutting")

    designs: list[str] = []
    if re.search(r"\bDesign\s*A\b|\bScenario\s*A\b|\bChassis\s*A\b|\bArm\s*1[AB]\b|\bLocus\s*[1234]\b", t, re.I):
        designs.append("A")
    if re.search(r"\bDesign\s*B\b|\bScenario\s*B\b|\bChassis\s*B\b|\bHLA-C\b|\bTYK2\b", t, re.I):
        designs.append("B")
    if re.search(r"\bDesign\s*X\b|\bIslet-X\b|\bChallenger\s*X\b|\bINS\b.*\bsignal|\bFasL\b", t, re.I):
        designs.append("X")
    if re.search(r"\bDesign\s*Y\b|\bDesign\s*Y-FullStack\b", t, re.I):
        designs.append("Y")
    if re.search(r"\bD1b(?:-Auto-C)?\b|\bDesign\s*D1b\b", t, re.I):
        designs.append("D1b")
    if not designs:
        designs = ["A", "B", "X"]
    return sorted(set(gates)), sorted(set(designs))


SUPERSEDED_PREPRINT_PMIDS_IN_CEG = frozenset({"37745505", "37905040", "39605713"})
CORRESPONDENCE_SECONDARY_PMIDS_IN_CEG = frozenset({"42430558"})
VALID_CEG_RECORD_CLASSES = frozenset({
    "empirical_primary_study",
    "trial_registry_protocol",
    "regulatory_or_guidance_document",
    "computational_model_or_simulation",
    "correspondence_or_secondary_summary",
    "superseded_preprint_duplicate",
})
CURATED_INELIGIBLE_EDGE_IDS = frozenset({
    "DOI:10.1016/j.celbio.2026.100455:L1",
    "M0060:L4",
    "M0079:L2",
    "M0097:L2",
    "M0105:L3",
    "M0155:L3",
    "M0178:L1",
    "M0178:L2",
    "M0382:L1",
    "M0528:L1",
    "M0545:L1",
    "M1430:L3",
    "M1439:L1",
    "M1566:L3",
    "M1867:L1",
    "M2895:L1",
    "M3020:L1",
    "PMID:10531307:L1",
    "PMID:16177096:L1",
    "PMID:24222351:L1",
    "PMID:25491100:L1",
    "PMID:38844131:L1",
    "R1_0134:L1",
    "R1_0134:L2",
    "R1_0137:L1",
    "R1_0137:L2",
    "R1_0189:L1",
    "R1_0224:L4",
    "R1_0238:L3",
    "R1_0317:L4",
    "R1_0378:L3",
    "R1_0390:L1",
    "R1_0393:L3",
    "M1421:L1",
    "M1422:L1",
    "M1422:L2",
    "M1453:L1",
    "M2036:L1",
    "M2449:L1",
    "PMID:19555051:L1",
    "M1522:L4",
    "PMID:17375075:L1",
    "PMID:42020736:L1",
    "PMID:42056515:L1",
    "PMID:42749803:L1",
    "M0189:L1",
    "M1407:L1",
    "PMID:36944335:L1",
    "PMID:39753706:L1",
    "PMID:41436438:L1",
    "PMID:41999746:L1",
    "M0051:L4",
    "M0154:L3",
    "M0159:L3",
    "M0206:L3",
    "M0241:L3",
    # Legacy edges downgraded under IIP-0035 after all cited outcomes reconciled to qualitative_no_pvalue (plus M0051:L3 & M0154:L2 unvalidated_transfer and M0206:L2 untested screening-policy novel_upgrade_hypothesis)
    "M0051:L1",
    "M0051:L2",
    "M0051:L3",
    "M0154:L1",
    "M0154:L2",
    "M0159:L1",
    "M0159:L2",
    "M0206:L2",
    "M0241:L1",
    "M0241:L2",
})
CURATED_NON_EMPIRICAL_OUTCOME_IDS = frozenset({
    "M0049:O1", "M0049:O2", "M0049:O3", "M0049:O4", "M0049:O5",
    "M0049:O6", "M0049:O7", "M0049:O8", "M0049:O9", "M0049:O10", "M0049:O15",
    "R1_0142:O10",
    "M0174:O3", "M0174:O7", "M0174:O8",
    "M0930:O7", "M0930:O12",
    "M0065:O13",
    "M1511:O3",
})


def classify_ceg_study_record(study: dict[str, Any]) -> str:
    """Deterministically classify a CEG study record into a typed record class."""
    id_type = str(study.get("id_type") or "").strip().lower()
    id_val = str(study.get("id_value") or "").strip().lower()
    pmid = str(study.get("pmid") or (study.get("id_value") if id_type == "pmid" else "") or "").strip()
    bio_model = str(study.get("biological_model") or "").strip().lower()
    imm_ctx = str(study.get("immune_context") or "").strip().lower()
    follow_up = str(study.get("follow_up") or "").strip().lower()
    ft = str(study.get("fidelity_tier") or "").strip()
    explicit_rc = str(study.get("record_class") or "").strip()
    if id_type == "nct" or id_val.startswith("nct"):
        return "trial_registry_protocol"
    if (
        id_type in ("url", "spl set id")
        or any(dom in id_val for dom in ("fda.gov", "federalregister.gov", "ema.europa.eu"))
        or bio_model.startswith(("regulatory guidance", "regulatory label", "regulatory text"))
        or explicit_rc == "regulatory_or_guidance_document"
    ):
        return "regulatory_or_guidance_document"
    if (
        study.get("superseded_by_queue_id")
        or study.get("superseded_by_pmid")
        or pmid in SUPERSEDED_PREPRINT_PMIDS_IN_CEG
    ):
        return "superseded_preprint_duplicate"
    if (
        study.get("publication_type") == "correspondence_or_secondary_summary"
        or pmid in CORRESPONDENCE_SECONDARY_PMIDS_IN_CEG
    ):
        return "correspondence_or_secondary_summary"
    if ft == "T6_Biophysical_Model_or_Regulatory" and (
        explicit_rc == "computational_model_or_simulation"
        or any(tok in imm_ctx for tok in ("computational model", "none modeled", "not applicable (computational"))
        or any(tok in follow_up for tok in ("not applicable (computational model)", "pseudo-steady-state model"))
        or (
            bio_model.startswith((
                "mathematical ",
                "two-dimensional axisymmetric finite-element",
                "three-dimensional finite-element",
            ))
            and "validation" not in bio_model
            and "fibroblast" not in bio_model
        )
    ):
        return "computational_model_or_simulation"
    if explicit_rc in VALID_CEG_RECORD_CLASSES and explicit_rc != "empirical_primary_study":
        return explicit_rc
    return "empirical_primary_study"


def build_layer1_atlas_from_ceg(ceg: dict[str, Any]) -> dict[str, Any]:
    """Deterministically derive the compact Layer-1 CEG atlas from `context_evidence_graph_v1.json`."""
    studies = ceg["studies"]
    arms = ceg["arms"]
    outcomes = ceg["outcomes"]
    edges = ceg["edges"]
    effectors = ceg["effectors"]
    graveyard = ceg["graveyard_index"]

    tier_studies = Counter(s.get("fidelity_tier") or "unspecified" for s in studies.values())
    tier_outcomes = Counter(o.get("manuscript_tier") or "unspecified" for o in outcomes.values())
    id_types = Counter(s.get("id_type") or "unspecified" for s in studies.values())
    record_class_counts = Counter(classify_ceg_study_record(s) for s in studies.values())
    r1_total_count = sum(1 for q in studies if q.startswith("R1_"))
    r1_empirical_count = sum(
        1 for q, s in studies.items() if q.startswith("R1_") and classify_ceg_study_record(s) == "empirical_primary_study"
    )
    non_r1_empirical_count = sum(
        1 for q, s in studies.items() if not q.startswith("R1_") and classify_ceg_study_record(s) == "empirical_primary_study"
    )
    non_r1_total_count = sum(1 for q in studies if not q.startswith("R1_"))

    direction_counts = Counter(o.get("direction") or "unspecified" for o in outcomes.values())
    stat_counts = Counter(o.get("stat_type") or "unspecified" for o in outcomes.values())
    cited_outcomes = sum(1 for o in outcomes.values() if o.get("cited_by_edge_ids"))

    relation_counts = Counter(e.get("relation") or "unspecified" for e in edges)
    incorporation_counts = Counter(e.get("incorporation") or "unspecified" for e in edges)
    applicability_counts = Counter(e.get("applicability") or "unspecified" for e in edges)

    by_gate: dict[str, dict[str, int]] = {}
    by_design: dict[str, dict[str, int]] = {}
    study_relations: dict[str, Counter[str]] = {q: Counter() for q in studies}
    study_gates: dict[str, set[str]] = {q: set() for q in studies}

    for e in edges:
        rel = str(e.get("relation") or "unspecified")
        qid = e["queue_id"]
        study_relations[qid][rel] += 1
        gates, designs = _classify_gates_and_designs(str(e.get("target_gate_or_module") or ""))
        for g in gates:
            study_gates[qid].add(g)
            bucket = by_gate.setdefault(g, {"total_links": 0})
            bucket["total_links"] += 1
            bucket[rel] = bucket.get(rel, 0) + 1
        for d_id in designs:
            db = by_design.setdefault(d_id, {"total_links": 0})
            db["total_links"] += 1
            db[rel] = db.get(rel, 0) + 1

    effector_summaries = []
    for eff_id in sorted(effectors):
        eff = effectors[eff_id]
        effector_summaries.append({
            "effector_id": eff["effector_id"],
            "topic_id": eff["topic_id"],
            "label": eff["label"],
            "category": eff["category"],
            "causal_isolation": eff["causal_isolation"],
            "mediation_status": eff["mediation_status"],
            "linked_studies_count": len(eff.get("linked_study_ids") or []),
            "linked_arms_count": len(eff.get("linked_variable_ids") or []),
            "linked_outcomes_count": len(eff.get("linked_outcome_tuples") or []),
        })

    graveyard_summaries = [
        {
            "graveyard_id": g["graveyard_id"],
            "category": g.get("category"),
        }
        for g in graveyard
    ]

    study_index = []
    for qid in sorted(studies):
        s = studies[qid]
        rec: dict[str, Any] = {
            "queue_id": qid,
            "record_class": classify_ceg_study_record(s),
            "id_type": s["id_type"],
            "id_value": s["id_value"],
            "fidelity_tier": s.get("fidelity_tier"),
            "arms_count": len(s.get("variable_arm_ids") or []),
            "outcomes_count": s.get("outcomes_count", 0),
            "links_count": s.get("links_count", 0),
            "gates_touched": sorted(study_gates.get(qid) or []),
            "link_relations": dict(sorted(study_relations[qid].items())),
        }
        if s.get("pub_year") is not None:
            rec["pub_year"] = s["pub_year"]
        for opt in ("pmid", "pmcid", "doi"):
            if s.get(opt):
                rec[opt] = s[opt]
        if s.get("candidate_effector_topics"):
            rec["effector_topics"] = s["candidate_effector_topics"]
        study_index.append(rec)

    return {
        "schema_version": "openislet_layer1_context_evidence_graph_atlas_v1.0",
        "status": "COMPLETE_LAYER1_CONTEXT_EVIDENCE_GRAPH_ATLAS_V1",
        "evidence_class": "STRUCTURED_LITERATURE_CONTEXT_EVIDENCE_GRAPH",
        "scope_and_scoring_separation_note": (
            f"Maps {len(studies):,} total CEG records ({record_class_counts['empirical_primary_study']:,} empirical primary studies, "
            f"{record_class_counts['trial_registry_protocol']:,} trial registry protocols, "
            f"{record_class_counts['regulatory_or_guidance_document']:,} regulatory/guidance documents, "
            f"{record_class_counts['computational_model_or_simulation']:,} computational/modeling-primary records, "
            f"{record_class_counts['correspondence_or_secondary_summary']:,} correspondence/secondary summary, and "
            f"{record_class_counts['superseded_preprint_duplicate']:,} superseded preprint duplicates), "
            f"{len(arms):,} experimental/control arms, {len(outcomes):,} structured outcomes, and {len(edges):,} blueprint links "
            f"across all five Islet-1 biological gates and 20 effector pathways, screened from "
            f"{ceg['screening_funnel']['screened_records_total']:,} deduplicated candidate records. "
            "The Layer-1 2,160-stack prior-anchored scenario evaluator retains its 604-study Ring-1 scoring baseline "
            "unchanged (+0.00 pp under IIP-0035) unless an explicit Islet Improvement Proposal "
            "updates a module parameter or Idea Graveyard rule."
        ),
        "summary": {
            "screened_candidate_records_total": ceg["screening_funnel"]["screened_records_total"],
            "candidates_advanced_total": ceg["screening_funnel"]["candidates_advanced_total"],
            "ceg_records_total": len(studies),
            "primary_empirical_studies_extracted": record_class_counts["empirical_primary_study"],
            "trial_registry_protocols_in_ceg": record_class_counts["trial_registry_protocol"],
            "regulatory_or_guidance_records_in_ceg": record_class_counts["regulatory_or_guidance_document"],
            "computational_models_in_ceg": record_class_counts["computational_model_or_simulation"],
            "correspondence_or_secondary_records_in_ceg": record_class_counts["correspondence_or_secondary_summary"],
            "superseded_preprint_duplicates_in_ceg": record_class_counts["superseded_preprint_duplicate"],
            "ring1_core_records_in_ceg_total": r1_total_count,
            "ring1_empirical_core_studies_in_ceg": r1_empirical_count,
            "new_expansion_empirical_studies_in_ceg": non_r1_empirical_count,
            "non_r1_ceg_records_total": non_r1_total_count,
            "ring1_scoring_baseline_studies_unchanged": 604,
            "experimental_and_control_arms_total": len(arms),
            "structured_outcomes_total": len(outcomes),
            "outcomes_cited_by_blueprint_links": cited_outcomes,
            "background_arm_outcomes_retained": len(outcomes) - cited_outcomes,
            "blueprint_links_total": len(edges),
            "biological_effector_nodes_total": len(effectors),
            "graveyard_rules_indexed_total": len(graveyard),
            "unresolved_outcome_citations": 0,
            "inline_verbatim_quotes_shipped": 0,
        },
        "screening_funnel": ceg["screening_funnel"],
        "record_class_distribution": dict(sorted(record_class_counts.items())),
        "study_identifier_types": dict(sorted(id_types.items())),
        "fidelity_tier_distribution": {
            "studies_by_manuscript_tier": dict(sorted(tier_studies.items())),
            "outcomes_by_manuscript_tier": dict(sorted(tier_outcomes.items())),
        },
        "outcome_distributions": {
            "by_direction": dict(sorted(direction_counts.items())),
            "by_stat_type": dict(sorted(stat_counts.items())),
        },
        "blueprint_link_distributions": {
            "by_relation": dict(sorted(relation_counts.items())),
            "by_incorporation": dict(sorted(incorporation_counts.items())),
            "by_applicability": dict(sorted(applicability_counts.items())),
        },
        "gate_and_design_claim_balance": {
            "by_gate": {k: dict(sorted(v.items())) for k, v in sorted(by_gate.items())},
            "by_design": {k: dict(sorted(v.items())) for k, v in sorted(by_design.items())},
        },
        "effector_atlas": effector_summaries,
        "graveyard_evidence_index": graveyard_summaries,
        "study_index": study_index,
    }


def verify_ceg_payload(ceg: dict[str, Any], repo_root: Path = REPO_ROOT, *, require_sources: bool = False) -> dict[str, int]:
    """Fail-closed structural, referential, and rights verification of `context_evidence_graph_v1.json`."""
    studies = ceg.get("studies")
    arms = ceg.get("arms")
    outcomes = ceg.get("outcomes")
    edges = ceg.get("edges")
    effectors = ceg.get("effectors")
    graveyard = ceg.get("graveyard_index")

    if not isinstance(studies, dict) or len(studies) != 1182:
        raise ValueError(f"Expected 1,182 studies in CEG, got {len(studies) if isinstance(studies, dict) else type(studies)}")
    if not isinstance(arms, dict) or len(arms) != 7749:
        raise ValueError(f"Expected 7,749 arms in CEG, got {len(arms) if isinstance(arms, dict) else type(arms)}")
    if not isinstance(outcomes, dict) or len(outcomes) != 14175:
        raise ValueError(f"Expected 14,175 outcomes in CEG, got {len(outcomes) if isinstance(outcomes, dict) else type(outcomes)}")
    if not isinstance(edges, list) or len(edges) != 2494:
        raise ValueError(f"Expected 2,494 edges in CEG, got {len(edges) if isinstance(edges, list) else type(edges)}")
    if not isinstance(effectors, dict) or len(effectors) != 20:
        raise ValueError(f"Expected 20 effectors in CEG, got {len(effectors) if isinstance(effectors, dict) else type(effectors)}")
    if not isinstance(graveyard, list) or len(graveyard) != 10:
        raise ValueError(f"Expected 10 graveyard entries in CEG, got {len(graveyard) if isinstance(graveyard, list) else type(graveyard)}")

    # Typed record_class integrity, supersession/same-patient linkage, and fail-closed certification gate
    def _study_pmid(st_obj: dict[str, Any]) -> str:
        return str(
            st_obj.get("pmid")
            or (st_obj.get("id_value") if str(st_obj.get("id_type") or "").strip().lower() == "pmid" else "")
            or ""
        ).strip()

    expected_rc_counts: Counter[str] = Counter()
    for qid, s in studies.items():
        expected_rc = classify_ceg_study_record(s)
        if s.get("record_class") != expected_rc:
            raise ValueError(f"Study {qid} has record_class={s.get('record_class')!r}, expected {expected_rc!r}")
        expected_rc_counts[expected_rc] += 1
        sup_qid = s.get("superseded_by_queue_id")
        sup_pmid = s.get("superseded_by_pmid")
        if expected_rc == "superseded_preprint_duplicate" and (not sup_qid or not sup_pmid):
            raise ValueError(f"Superseded duplicate study {qid} must declare superseded_by_queue_id and superseded_by_pmid")
        if sup_qid is not None:
            if not isinstance(sup_qid, str) or sup_qid == qid or sup_qid not in studies:
                raise ValueError(f"Study {qid} superseded_by_queue_id={sup_qid!r} does not resolve to another existing study")
            target_pmid = _study_pmid(studies[sup_qid])
            if not isinstance(sup_pmid, str) or not sup_pmid.strip() or sup_pmid.strip() != target_pmid:
                raise ValueError(
                    f"Study {qid} superseded_by_pmid={sup_pmid!r} mismatches target study {sup_qid} PMID {target_pmid!r}"
                )
            _seen_sup = {qid}
            _cur_sup = sup_qid
            while _cur_sup is not None:
                if _cur_sup in _seen_sup:
                    raise ValueError(f"Study {qid} superseded_by_queue_id chain forms a supersession cycle at {_cur_sup!r}")
                _seen_sup.add(_cur_sup)
                _nxt = studies[_cur_sup].get("superseded_by_queue_id") if _cur_sup in studies else None
                _cur_sup = _nxt if isinstance(_nxt, str) and _nxt in studies else None
        elif sup_pmid is not None:
            raise ValueError(f"Study {qid} declares superseded_by_pmid={sup_pmid!r} without superseded_by_queue_id")
        _st_sp = s.get("source_pin")
        if _st_sp is not None and not isinstance(_st_sp, dict):
            raise ValueError(f"Study {qid} source_pin must be a dict, got {type(_st_sp).__name__}")
        for _link_field in ("same_patient_followup_of_queue_id", "same_patient_linked_queue_id"):
            _lqid = s.get(_link_field)
            if _lqid is not None:
                if not isinstance(_lqid, str) or _lqid == qid or _lqid not in studies:
                    raise ValueError(f"Study {qid} {_link_field}={_lqid!r} does not resolve to another existing study")
                if _link_field == "same_patient_followup_of_queue_id" and s.get("same_patient_followup_of_pmid") is not None:
                    _lpmid = str(s.get("same_patient_followup_of_pmid") or "").strip()
                    _tpmid = _study_pmid(studies[_lqid])
                    if _lpmid != _tpmid:
                        raise ValueError(f"Study {qid} same_patient_followup_of_pmid={_lpmid!r} mismatches target {_lqid} PMID {_tpmid!r}")
    if (ceg.get("counts") or {}).get("record_class_counts") != dict(sorted(expected_rc_counts.items())):
        raise ValueError("CEG counts.record_class_counts disagrees with deterministic study record_class classification")
    _sf = ceg.get("screening_funnel") or {}
    if not isinstance(_sf, dict):
        raise ValueError("CEG screening_funnel must be a dict")
    _sf_total = _sf.get("screened_records_total")
    _sf_r1 = _sf.get("ring1_existing_records_screened")
    _sf_new = _sf.get("new_candidate_records_screened")
    _sf_adv = _sf.get("candidates_advanced_total")
    _sf_park = _sf.get("contextual_reviews_parked_total")
    _sf_excl = _sf.get("off_topic_records_excluded_total")
    _sf_w1 = _sf.get("wave1_ceg_records_extracted_total")
    _sf_rem = _sf.get("remaining_supporting_catalog_queued_for_later_waves")
    if not all(type(v) is int for v in (_sf_total, _sf_r1, _sf_new, _sf_adv, _sf_park, _sf_excl, _sf_w1, _sf_rem)):
        raise ValueError("CEG screening_funnel missing required integer totals")
    if _sf_r1 + _sf_new != _sf_total:
        raise ValueError(f"CEG screening_funnel mismatch: ring1 ({_sf_r1}) + new ({_sf_new}) != screened_records_total ({_sf_total})")
    if _sf_adv + _sf_park + _sf_excl != _sf_total:
        raise ValueError(f"CEG screening_funnel mismatch: advanced ({_sf_adv}) + parked ({_sf_park}) + excluded ({_sf_excl}) != screened_records_total ({_sf_total})")
    if _sf_w1 != len(studies) or _sf_w1 + _sf_rem != _sf_adv:
        raise ValueError(f"CEG screening_funnel mismatch: wave1 ({_sf_w1}) + remaining ({_sf_rem}) != candidates_advanced_total ({_sf_adv})")

    # Foreign-key integrity
    for qid, s in studies.items():
        for vid in s.get("variable_arm_ids") or []:
            if vid not in arms:
                raise ValueError(f"Study {qid} references missing arm {vid}")
    for vid, v in arms.items():
        if v.get("queue_id") not in studies:
            raise ValueError(f"Arm {vid} references missing study {v.get('queue_id')}")
    for oid, o in outcomes.items():
        if o.get("queue_id") not in studies:
            raise ValueError(f"Outcome {oid} references missing study {o.get('queue_id')}")
        if "verbatim_quote" in o and o["verbatim_quote"] is not None:
            raise ValueError(f"Outcome {oid} retains a raw verbatim_quote in public CEG")
        if o.get("quote_sha256") is not None and not re.fullmatch(r"[0-9a-f]{64}", str(o["quote_sha256"])):
            raise ValueError(f"Outcome {oid} has malformed quote_sha256")
        _o_sp = o.get("source_pin")
        if not isinstance(_o_sp, dict) or not isinstance(_o_sp.get("source_type"), str) or not _o_sp["source_type"].strip():
            raise ValueError(f"Outcome {oid} source_pin must be a dict with non-empty source_type, got {_o_sp!r}")
        _m_ep_stat = re.search(r"\[Structured endpoint [^\]]*\bstat=([^,|\s\]]+)", str(o.get("endpoint_name") or ""))
        if _m_ep_stat and _m_ep_stat.group(1) != o.get("stat_type"):
            raise ValueError(f"Outcome {oid} endpoint_name structured stat={_m_ep_stat.group(1)!r} disagrees with stat_type={o.get('stat_type')!r}")
        _m_qs_ds = re.search(r"\[Structured outcome [^\]]*\bdirection=([^,|\s\]]+)\s*(?:\||,)\s*stat=([^,|\s\]]+)", str(o.get("quantitative_summary") or ""))
        if _m_qs_ds and (_m_qs_ds.group(1) != o.get("direction") or _m_qs_ds.group(2) != o.get("stat_type")):
            raise ValueError(
                f"Outcome {oid} quantitative_summary structured (direction, stat)=({_m_qs_ds.group(1)!r}, {_m_qs_ds.group(2)!r}) "
                f"disagrees with ({o.get('direction')!r}, {o.get('stat_type')!r})"
            )
        if studies[o["queue_id"]].get("record_class") == "empirical_primary_study":
            _o_narr = f"{o.get('endpoint_name', '')} {o.get('quantitative_summary', '')}"
            if re.search(
                r"(?i)(\btitle[- ]only\b|\btitle\s+claim\b|\bqualitative\s+title\s+statement\s+only\b"
                r"|\bfrom\s+title\b|\bstated\s+only\s+in\s+the\s+title\b|\bread\s+from\s+the\s+title\s+alone\b"
                r"|\bthe\s+title\s+reports\b.*\bno\s+.*data\s+are\s+stated\b|\babstract\s+text\s+itself\s+does\s+not\s+report\b)",
                _o_narr,
            ) and _o_sp["source_type"] not in ("title_or_citation_metadata_only", "title"):
                raise ValueError(
                    f"Outcome {oid} narrative states title-only / unmeasured content while source_pin.source_type={_o_sp['source_type']!r}"
                )
            _o_sim_blob = f"{_o_narr} {o.get('comparison_arms', '')} {o.get('intervention_value', '')}"
            if (
                re.search(
                    r"(?i)(\bfinite[- ]element\s+(?:prediction|model)\b|\bmodel[- ]predicted\b"
                    r"|\bmodel\s+prediction,\s*not\s+a\s+measurement\b|\bsimulation\s+only\b|\bsimulation\s+parameter\b"
                    r"|\bsimulated\s+(?:core|critical|central|viability|oxygen|average|fraction|hypoxic|human|control|sonic|islet)\b"
                    r"|\bcfd\s+model\b|\bmathematical\s+safety\s+model\b|\bdescriptive\s+modeled\s+bank\b"
                    r"|\bmodelled\s+estimate\b|\blogistic[- ]regression\s+prediction\b"
                    r"|\bpo2_survival\b|\bmodels\s+barcode[- ]sharing\s+probabilities\b|\bmodeled\s+probability\b"
                    r"|\bpredicted\s+stability\b|\bextrapolation,\s*not\s+long[- ]term\s+observation\b"
                    r"|\bmonte\s+carlo\s+simulation\b|\bthiele[- ]modulus\s+calculation\b|\bcalculated\s+volumetric\s+ocr\b"
                    r"|\buse\s+a\s+simple\s+model\s+to\s+illustrate\b|\bmodeled\s+sequencing[- ]depth\s+saturation\b"
                    r"|\bchao2[- ]estimated\b|\btheoretical\s+limit\s+of\s+detection\b|\ban\s+estimate,\s*not\s+a\s+spike[- ]in\s+measurement\b"
                    r"|\bsuggested,\s*not\s+measured\b|\btheoretical\s+calculation\b"
                    r"|\b11[- ]gene\s+multiple\s+regression\s+model\b"
                    r"|\bpredicted\s+protein[- ]coding\s+consequence\b|\bestimated\s+test\s+sample\s+size\b)",
                    _o_sim_blob,
                )
                or re.search(r"\bA\d+_(?:model|sim|FEM|CFD|Model)_", str(o.get("comparison_arms") or ""))
                or oid in CURATED_NON_EMPIRICAL_OUTCOME_IDS
            ) and o.get("evidence_provenance") not in ("computational_model_or_simulation", "mixed_measured_and_simulated_model", "mixed_empirical_and_computational"):
                raise ValueError(
                    f"Outcome {oid} describes computational/simulated prediction without non-empirical evidence_provenance ({o.get('evidence_provenance')!r})"
                )
    _non_emp_prov = frozenset({"computational_model_or_simulation", "mixed_measured_and_simulated_model", "mixed_empirical_and_computational"})
    _title_only_pins = frozenset({"title_or_citation_metadata_only", "title"})
    _link_indices_by_study: dict[str, list[int]] = {}
    _fwd_cited_by: dict[str, set[str]] = {}
    _seen_edge_ids: set[str] = set()
    for e in edges:
        _eid = e.get("edge_id")
        if not isinstance(_eid, str) or not _eid.strip():
            raise ValueError(f"Edge has missing or non-string edge_id: {_eid!r}")
        if _eid in _seen_edge_ids:
            raise ValueError(f"Duplicate edge_id {_eid!r} in CEG edges")
        _seen_edge_ids.add(_eid)
        qid = e.get("queue_id")
        if qid not in studies:
            raise ValueError(f"Edge {e.get('edge_id')} references missing study {qid}")
        _lidx = e.get("link_index")
        if type(_lidx) is not int or _lidx < 0:
            raise ValueError(f"Edge {e.get('edge_id')} has invalid link_index={_lidx!r}")
        _link_indices_by_study.setdefault(qid, []).append(_lidx)
        _pce = e.get("primary_certification_eligible")
        if type(_pce) is not bool:
            raise ValueError(
                f"Edge {e.get('edge_id')} primary_certification_eligible must be strict bool, got {_pce!r} ({type(_pce).__name__})"
            )
        if studies[qid].get("record_class") != "empirical_primary_study" and _pce is not False:
            raise ValueError(f"Non-empirical study edge {e.get('edge_id')} ({studies[qid].get('record_class')}) cannot be primary_certification_eligible=True")
        if _pce is True and e.get("applicability") == "domain_mismatch":
            raise ValueError(f"Domain-mismatched edge {e.get('edge_id')} (applicability='domain_mismatch') cannot be primary_certification_eligible=True")
        if _pce is True and e.get("relation") == "noninformative":
            raise ValueError(f"Noninformative edge {e.get('edge_id')} (relation='noninformative') cannot be primary_certification_eligible=True")
        if _pce is True and e.get("readiness") == "needs_primary_fulltext":
            raise ValueError(f"Unauthenticated/title-only edge {e.get('edge_id')} (readiness='needs_primary_fulltext') cannot be primary_certification_eligible=True")
        _st_pin = studies[qid].get("source_pin")
        if _pce is True and isinstance(_st_pin, dict) and _st_pin.get("source_type") in _title_only_pins:
            raise ValueError(f"Title-only study edge {e.get('edge_id')} cannot be primary_certification_eligible=True")
        if _pce is True and e.get("evidence_provenance") in _non_emp_prov:
            raise ValueError(f"Simulated/modeled edge {e.get('edge_id')} ({e.get('evidence_provenance')}) cannot be primary_certification_eligible=True")
        if _pce is True and _eid in CURATED_INELIGIBLE_EDGE_IDS:
            raise ValueError(f"Curated-ineligible edge {_eid} cannot be primary_certification_eligible=True")
        if _pce is True and re.search(r"(?i)(\btitle[- ]only\b|\bno\s+measured\s+result\b|\bread\s+from\s+the\s+title\s+alone\b)", str(e.get("claim_summary") or "")):
            raise ValueError(f"Edge {e.get('edge_id')} narrative states title-only / no measured result while primary_certification_eligible=True")
        _cited_tuples = e.get("cited_outcome_tuples") or []
        if len(_cited_tuples) != len(set(_cited_tuples)):
            raise ValueError(f"Edge {e.get('edge_id')} cited_outcome_tuples contains duplicate outcome(s): {_cited_tuples!r}")
        if _pce is True and not _cited_tuples:
            raise ValueError(f"Primary-certification-eligible edge {e.get('edge_id')} must cite at least one outcome (cited_outcome_tuples is empty)")
        for tok in _cited_tuples:
            if not isinstance(tok, str) or tok not in outcomes:
                raise ValueError(f"Edge {e.get('edge_id')} cites unresolved outcome {tok}")
            _fwd_cited_by.setdefault(tok, set()).add(str(e.get("edge_id")))
            if _pce is True:
                _oc = outcomes[tok]
                _oc_qid = _oc.get("queue_id")
                if _oc_qid not in studies or studies[_oc_qid].get("record_class") != "empirical_primary_study":
                    raise ValueError(
                        f"Primary-certification-eligible edge {e.get('edge_id')} cites outcome {tok} owned by non-empirical/duplicate study {_oc_qid} ({studies.get(_oc_qid, {}).get('record_class')})"
                    )
                _oc_st_pin = studies[_oc_qid].get("source_pin")
                if isinstance(_oc_st_pin, dict) and _oc_st_pin.get("source_type") in _title_only_pins:
                    raise ValueError(
                        f"Primary-certification-eligible edge {e.get('edge_id')} cites outcome {tok} owned by title-only study {_oc_qid} (source_type={_oc_st_pin.get('source_type')!r})"
                    )
                if _oc.get("evidence_provenance") in _non_emp_prov:
                    raise ValueError(
                        f"Primary-certification-eligible edge {e.get('edge_id')} cites computational/simulated or mixed provenance outcome {tok} ({_oc.get('evidence_provenance')})"
                    )
                _oc_pin_t = (_oc.get("source_pin") or {}).get("source_type") if isinstance(_oc.get("source_pin"), dict) else None
                if _oc_pin_t in _title_only_pins:
                    raise ValueError(
                        f"Primary-certification-eligible edge {e.get('edge_id')} cites title-only outcome {tok} (source_type={_oc_pin_t!r})"
                    )
        _cited_local = e.get("cited_outcome_ids")
        if not isinstance(_cited_local, list):
            raise ValueError(f"Edge {e.get('edge_id')} cited_outcome_ids must be a list")
        if len(_cited_local) != len(set(_cited_local)):
            raise ValueError(f"Edge {e.get('edge_id')} cited_outcome_ids contains duplicate outcome(s): {_cited_local!r}")
        if _pce is True and not _cited_local:
            raise ValueError(f"Primary-certification-eligible edge {e.get('edge_id')} must cite at least one outcome (cited_outcome_ids is empty)")
        _resolved_local: set[str] = set()
        for _loc in _cited_local:
            if not isinstance(_loc, str) or not _loc.strip():
                raise ValueError(f"Edge {e.get('edge_id')} has invalid cited_outcome_id {_loc!r}")
            _full_loc = _loc if ":" in _loc else f"{qid}:{_loc}"
            if _full_loc not in outcomes:
                raise ValueError(f"Edge {e.get('edge_id')} cited_outcome_id {_loc!r} resolves to unresolved outcome {_full_loc!r}")
            _resolved_local.add(_full_loc)
        if _resolved_local != set(_cited_tuples):
            raise ValueError(
                f"Edge {e.get('edge_id')} cited_outcome_ids {sorted(_resolved_local)} mismatches cited_outcome_tuples {sorted(set(_cited_tuples))}"
            )
    for oid, _oc_sym in outcomes.items():
        _bwd_edges = _oc_sym.get("cited_by_edge_ids")
        if _bwd_edges is not None and not isinstance(_bwd_edges, list):
            raise ValueError(f"Outcome {oid} cited_by_edge_ids must be a list or omitted, got {type(_bwd_edges).__name__}")
        if _bwd_edges is not None and len(_bwd_edges) != len(set(_bwd_edges)):
            raise ValueError(f"Outcome {oid} cited_by_edge_ids contains duplicate edge_id(s): {_bwd_edges!r}")
        _bwd_set = set(_bwd_edges or [])
        _fwd_set = _fwd_cited_by.get(oid, set())
        if _bwd_set != _fwd_set:
            raise ValueError(
                f"Outcome {oid} cited_by_edge_ids {sorted(_bwd_set)} is asymmetric with citing edges {sorted(_fwd_set)}"
            )
    for qid, _idxs in _link_indices_by_study.items():
        if sorted(_idxs) != list(range(len(_idxs))):
            raise ValueError(f"Study {qid} has non-contiguous or out-of-range link_index sequence {sorted(_idxs)} (expected 0..{len(_idxs)-1})")
        if studies[qid].get("links_count") != len(_idxs):
            raise ValueError(f"Study {qid} links_count={studies[qid].get('links_count')} disagrees with {len(_idxs)} edges")
    for eff_id, eff in effectors.items():
        for sq in eff.get("linked_study_ids") or []:
            if sq not in studies:
                raise ValueError(f"Effector {eff_id} references missing study {sq}")
        for sv in eff.get("linked_variable_ids") or []:
            if sv not in arms:
                raise ValueError(f"Effector {eff_id} references missing arm {sv}")
        for tok in eff.get("linked_outcome_tuples") or []:
            if not isinstance(tok, str) or tok not in outcomes:
                raise ValueError(f"Effector {eff_id} references missing outcome {tok}")

    gy_path = repo_root / GRAVEYARD_REL_PATH
    if gy_path.is_file():
        gy_ids = {row["graveyard_id"] for row in json.loads(gy_path.read_text(encoding="utf-8")).get("graveyard_entries", [])}
        for g in graveyard:
            if g.get("graveyard_id") not in gy_ids:
                raise ValueError(f"CEG graveyard entry {g.get('graveyard_id')} missing from {GRAVEYARD_REL_PATH}")

    # Global study<->arm ownership symmetry, outcomes_count, and same-study explicit outcome arm FKs
    _actual_out_counts: dict[str, int] = {}
    for _oid_g, _o_g in outcomes.items():
        _oqid_g = _o_g.get("queue_id")
        _actual_out_counts[_oqid_g] = _actual_out_counts.get(_oqid_g, 0) + 1
        for _af in ("comparator_arm_id", "focal_arm_id", "intervention_arm_id"):
            _av = _o_g.get(_af)
            if _av is not None:
                if not isinstance(_av, str) or not _av.strip():
                    raise ValueError(f"Outcome {_oid_g} has invalid {_af}={_av!r}")
                _full_av = _av if _av.startswith("VAR:") else f"VAR:{_oqid_g}:{_av}"
                if _full_av not in arms:
                    raise ValueError(f"Outcome {_oid_g} {_af}={_av!r} does not resolve in CEG arms")
                if arms[_full_av].get("queue_id") != _oqid_g:
                    raise ValueError(f"Outcome {_oid_g} {_af}={_av!r} belongs to study {arms[_full_av].get('queue_id')!r}, not {_oqid_g!r}")
    for _qid_g, _st_g in studies.items():
        _st_vids = _st_g.get("variable_arm_ids") or []
        if len(_st_vids) != len(set(_st_vids)):
            raise ValueError(f"Study {_qid_g} has duplicate entries in variable_arm_ids")
        for _vid_g in _st_vids:
            if arms[_vid_g].get("queue_id") != _qid_g:
                raise ValueError(f"Study {_qid_g} variable_arm_ids references arm {_vid_g} owned by {arms[_vid_g].get('queue_id')!r}")
        if _st_g.get("outcomes_count") != _actual_out_counts.get(_qid_g, 0):
            raise ValueError(f"Study {_qid_g} outcomes_count={_st_g.get('outcomes_count')} disagrees with actual {_actual_out_counts.get(_qid_g, 0)}")
    for _vid_g, _a_g in arms.items():
        _aq_g = _a_g.get("queue_id")
        if _vid_g not in set(studies[_aq_g].get("variable_arm_ids") or []):
            raise ValueError(f"Arm {_vid_g} owned by {_aq_g} is missing from study variable_arm_ids")

    # Strict IIP-0034 & IIP-0035 new-record typed arm FKs, same-study joins, effector refs, no-partial-dossier span audit & full-leaf rolling-40 _SourceGuard
    _iip0034_qids = frozenset({"PMID:17375075", "PMID:42020736", "PMID:42056515", "PMID:42749803"})
    _iip0034_eids = frozenset({"PMID:17375075:L1", "PMID:42020736:L1", "PMID:42056515:L1", "PMID:42749803:L1", "M1522:L4"})
    _iip0034_req_pmids = ("17375075", "39696686", "42020736", "42056515", "42749803")
    _iip0035_qids = frozenset({"M0189", "M1407", "PMID:36944335", "PMID:39753706", "PMID:41436438", "PMID:41999746"})
    _iip0035_expected_edge_bindings = {
        "M0189:L1": ("M0189", ("M0189:O1", "M0189:O2"), "9489a2d260775e0591af72799e3a2297260f4b2df14dc29cdc59acbad01e7c7b"),
        "M1407:L1": ("M1407", ("M1407:O1",), "97edc6fff8dca53a997b4d6c544eebad6253eba29f0cd720cf335994947f451a"),
        "PMID:36944335:L1": ("PMID:36944335", ("PMID:36944335:O1",), "8177ff0a7f4837d6609155d02f26064aa674e6d6433d5519344421888b80d8cb"),
        "PMID:39753706:L1": ("PMID:39753706", ("PMID:39753706:O1", "PMID:39753706:O2"), "7f8df9e1e6c03a9d386e2e57c0592bc28fe88d8694487daf36d324ef1efe40cf"),
        "PMID:41436438:L1": ("PMID:41436438", ("PMID:41436438:O1",), "05ce92c5032f16e4b53d438cdcf9b63b4251d57a8769661be819ce1a104379d6"),
        "PMID:41999746:L1": ("PMID:41999746", ("PMID:41999746:O1",), "9a30fdc5656e490bf78c93d018a740d8b67eb27c2059a73a0059279db2db9e33"),
        "M0051:L4": ("M0051", ("M0051:O1", "M0051:O5"), "fa7fb021669e6a71f22cc4078dd477bda9904c47e819b2c250a9e2081f9962bf"),
        "M0154:L3": ("M0154", ("M0154:O1", "M0154:O3", "M0154:O4"), "dfc1fb8e65f7db6084552b462c21532d9cab4d326ca13eef0f280e50e4cd5f2c"),
        "M0159:L3": ("M0159", ("M0159:O2",), "1c18c658d8351a5128abb0b5a21846fa2a15fe4001b0b733cee5d650d51de8cc"),
        "M0206:L3": ("M0206", ("M0206:O1", "M0206:O2", "M0206:O9", "M0206:O11", "M0206:O12", "M0206:O17"), "9c3f1e834d57b4683fcb135274c0cc00e482f0752ba61636b1dd595b75d796b5"),
        "M0241:L3": ("M0241", ("M0241:O1", "M0241:O2", "M0241:O4"), "60a3e30568b5bd329c5595fbfb5fa779c344f597f0ba8f0bc84c049f9f24c00f"),
    }
    _iip0035_eids = frozenset(_iip0035_expected_edge_bindings)
    _iip0035_adjudicated_inelig_legacy_outcomes = {
        "M0051:L1": (("M0051:O1", "M0051:O5"), "fda2fee5e96b1ae2f699dea5a5919eef81acb185440cfdb446157be83fff582b"),
        "M0051:L2": (("M0051:O1", "M0051:O5"), "b0b8d71b045b709de53cc8d8975a44f88bf343886f9944a8775f4f107b0acb29"),
        "M0051:L3": (("M0051:O4",), "a1b650ca58d7ee3f65d7ce05ad380e3bcbe2de5515f58b2ad8ce2f5b9ce6e147"),
        "M0154:L1": (("M0154:O1", "M0154:O2"), "b0b04b722d7ae4282cf253721cc6c9a536b39ba352dea2d69f2711eb947c5966"),
        "M0154:L2": (("M0154:O1", "M0154:O3", "M0154:O4"), "54844b79b54e02a1ba799f3bc4cae9366b949d7ad737ab0ff317aec213f36667"),
        "M0159:L1": (("M0159:O2", "M0159:O3"), "eec8033b2915fc19c425683bff865193b98e0534e520b52a50e07bb2cd980081"),
        "M0159:L2": (("M0159:O1", "M0159:O2", "M0159:O3"), "c4ab223e2234ae63414688b3c3c0110cfba5f2f15c189331945cc3afccee5258"),
        "M0206:L2": (("M0206:O9", "M0206:O12", "M0206:O10", "M0206:O14", "M0206:O17"), "bf2d25dc5e33e302309c364e95c85c87934028065912bf4673e67b27c91de241"),
        "M0241:L1": (("M0241:O2", "M0241:O3"), "0415b7d3450a99d88ec2b370284472391a5b9c4399f1f72e05497c6eaac7f255"),
        "M0241:L2": (("M0241:O1", "M0241:O4"), "fe8c8f2f7c3d40ba478a86ce88e9a9c49d7ca24f6df3a2f1b9a6f3fa807554bb"),
    }
    _iip0035_m206_l1_expected_sha = "c8f8cb8f03df508b94717777577b0e2b76754a49c0b55140cf7f02f2d57224c1"
    _iip0035_expected_arm_digests = {
        "VAR:M0189:A1": "9020bdec3108bc526a05cf1e13398c37e0b7d1c01d2f7a0b8a5129e8ff7dbd19",
        "VAR:M0189:A2_allogeneic_ctrl": "d9a3ab49db3a1550ec72677e231e8f78c0669498b1ce87f9b3cc15809d9af47e",
        "VAR:M1407:A1": "80003667e1218301f14566c81b6b226f12f08b72e8ebff76ea177c81e2e11ee8",
        "VAR:PMID:36944335:A1": "223640329c6632b6a2bd0adb299ac60ed8629b709f0bd573a93abfd9335663d1",
        "VAR:PMID:39753706:A1": "12df2f6220286f5d8e98d4f9f950e2dfc95f30e207cb367fb44e8a483ab195c5",
        "VAR:PMID:39753706:A2_original_protocol_ctrl": "97ce0849c5fe173b9707d41139d258cc2c2038dff728c2301838e51fb7a776ad",
        "VAR:PMID:41436438:A1": "47e244e5b004fbfad064117b9e2becb09f127252c705e773f4b59a30927e2776",
        "VAR:PMID:41999746:A1": "f050b9f19a1ca640dc14637faf81855973cb5dc3c6376c0429e9af0371c9db03",
        "VAR:PMID:41999746:A2_p53_suppressed": "aff8562531c232b9f4f0ec8901053ba94fe42ce918a62eadf4ced84ddfff0c50",
    }
    _iip0035_expected_study_digests = {
        "M0051": "6531ee32a16d308915a5d42f3a5cf97e74260b6ed4406c14042c218aad55ff94",
        "M0154": "493ea738a959291d4fe95f298556f4e22a447a91889d38c39b97a9322d1d0f22",
        "M0159": "02ab4902a6c9f52b7033855e4c6e112fb774c899d936ddfe54f899663a8299da",
        "M0189": "9d8e7e951be7eecbe9761771166fbf91610828a6f468b91c7851a8ed963c652d",
        "M0206": "5b6984efae346da92c8d6c8f322c88eb0d3f171dcc23d841b7d7f2357449ae9f",
        "M0241": "0d02ded3504df0ca0df1743c25cec1529a160873e4b7b486b8c37b9b91216bba",
        "M1407": "2b1fe5f1b6c0b09e6d81c7e0519e2862d9015865c8eeeed23f90e6bc6b567274",
        "PMID:36944335": "27afa1ef3895cb36621504ab3dda8ac703b2416d5b5675b5c60dc6a093db2105",
        "PMID:39753706": "619cf75493ec231fb729348c03785c095019eaa3e5211376abde09505eaa594c",
        "PMID:41436438": "e7a680c29f0df638f85075d1475170c1a923241cd309910a34f5da982a338935",
        "PMID:41999746": "b370082130c9c8fedd438db6fa718e0c5978832116262baf92c07cbf7bb67e94",
    }
    _iip0035_expected_dossier_digests = {
        "2642862": "a463f7f1770ed78ca7c892db54c59b2970cc83226cd48e71302218b9b48d89fa",
        "31427818": "f6593b87fda7d0bc1211d17d30288dbb86f6b408af77482adf4022e7531c1b98",
        "35953586": "6ff5efa024e5614f8a98125f51fbb3958b90c36dc1426293328b6537fa3d678b",
        "36944335": "ef775e8802dc187a5e04949fb88c2dda258044d9a205849a08699ac01a3eb03a",
        "38843830": "cc6ba867300b9c420070dd7bb25a7a59e7c5395269388652b50ca567f24611d1",
        "39753706": "948c18dddd97935f2c0cb6fa73fd905c8dc2aff31e7b64dfbe39d574059c552e",
        "41201433": "184652ef2fdbc682abf441f25930a58bb990f238e4f3e8846bf4c229ae8dd544",
        "41436438": "9db3a80db0a50ded4b7d8054c451b9ed017ffe3e17201808eef81b7d84b2e787",
        "41690310": "b55f46afb7307ff6793ad519b06a5735aee77a0b610cdeda0c8ed0d03fec3d5f",
        "41764741": "e7095034b47de98edcf0c2c8a46a1659600e2299a8928880b56b81ec6a58dbe6",
        "41999746": "dd2f6aac6a4812460da51ebae389980cf025d626ea41f810bd985d5193c67aa1",
    }
    _iip0035_adjudicated_inelig_legacy = frozenset(_iip0035_adjudicated_inelig_legacy_outcomes)
    _iip0035_req_pmids = ("2642862", "31427818", "35953586", "36944335", "38843830", "39753706", "41201433", "41436438", "41690310", "41764741", "41999746")
    _iip0035_crosslink_eids = frozenset({"M0051:L4", "M0154:L3", "M0159:L3", "M0206:L3", "M0241:L3"})
    _strict_study_expected_pmid = {
        "PMID:17375075": "17375075", "PMID:42020736": "42020736", "PMID:42056515": "42056515", "PMID:42749803": "42749803",
        "M0189": "38843830", "M1407": "2642862", "PMID:36944335": "36944335", "PMID:39753706": "39753706",
        "PMID:41436438": "41436438", "PMID:41999746": "41999746",
    }
    _strict_new_qids = _iip0034_qids | _iip0035_qids
    _strict_new_eids = _iip0034_eids | _iip0035_eids
    _strict_req_pmids = _iip0034_req_pmids + _iip0035_req_pmids
    _edges_by_eid_35 = {e.get("edge_id"): e for e in edges if isinstance(e, dict) and isinstance(e.get("edge_id"), str)}
    _missing_35_eids = sorted(eid for eid in _iip0035_expected_edge_bindings if eid not in _edges_by_eid_35)
    if _missing_35_eids:
        raise ValueError(f"Missing required IIP-0035 CEG edges: {_missing_35_eids}")
    for _eid35, (_exp_owner35, _exp_outs35, _exp_esha35) in sorted(_iip0035_expected_edge_bindings.items()):
        _eobj35 = _edges_by_eid_35[_eid35]
        if _eobj35.get("queue_id") != _exp_owner35:
            raise ValueError(f"IIP-0035 edge {_eid35} owner queue_id={_eobj35.get('queue_id')!r} mismatches expected {_exp_owner35!r}")
        if tuple(_eobj35.get("cited_outcome_tuples") or ()) != _exp_outs35:
            raise ValueError(f"IIP-0035 edge {_eid35} cited_outcome_tuples={_eobj35.get('cited_outcome_tuples')!r} mismatches expected {list(_exp_outs35)!r}")
        if _eobj35.get("primary_certification_eligible") is not False:
            raise ValueError(f"IIP-0035 edge {_eid35} must keep primary_certification_eligible=False")
        if _eid35 in _iip0035_crosslink_eids:
            _pn35 = _eobj35.get("provenance_note")
            if (
                not isinstance(_pn35, str)
                or len(_pn35.strip()) < 40
                or not re.search(r"(?i)\b(?:reused|legacy|prior)\b", _pn35)
                or not re.search(r"(?i)\b(?:outside|non-certifying|primary_certification_eligible\s*=\s*false|scoping)\b", _pn35)
            ):
                raise ValueError(f"IIP-0035 cross-link edge {_eid35} missing grounded provenance_note disclosing reused legacy outcomes and non-certifying scoping status")
    for _leg_eid, (_leg_outs, _leg_sha) in sorted(_iip0035_adjudicated_inelig_legacy_outcomes.items()):
        if _leg_eid not in _edges_by_eid_35:
            raise ValueError(f"Missing adjudicated legacy edge {_leg_eid}")
        _le_o = _edges_by_eid_35[_leg_eid]
        if _le_o.get("primary_certification_eligible") is not False:
            raise ValueError(f"Adjudicated legacy ineligible edge {_leg_eid} must have primary_certification_eligible=False")
        if tuple(_le_o.get("cited_outcome_tuples") or ()) != _leg_outs:
            raise ValueError(f"Adjudicated legacy edge {_leg_eid} cited_outcome_tuples={_le_o.get('cited_outcome_tuples')!r} mismatches pinned {list(_leg_outs)!r}")
        if not isinstance(_le_o.get("legacy_certification_adjudication_note"), str) or len(_le_o["legacy_certification_adjudication_note"].strip()) < 25:
            raise ValueError(f"Adjudicated legacy edge {_leg_eid} missing required legacy_certification_adjudication_note")
    _iip0035_reclass_qual_oids = frozenset({
        "M0051:O1", "M0051:O2", "M0051:O3", "M0051:O4", "M0051:O5", "M0051:O6",
        "M0154:O1", "M0154:O2", "M0154:O3", "M0154:O4",
        "M0159:O1", "M0159:O2", "M0159:O3", "M0206:O16", "M0206:O17",
        "M0241:O1", "M0241:O2", "M0241:O3", "M0241:O4",
    })
    _QUANTITATIVE_CEG_STAT_TYPES = frozenset({
        "adjusted_difference_95ci", "descriptive_scalar", "fold_change", "hazard_or_odds_ratio",
        "mean_sd", "mean_sem", "mean_unspecified_spread", "median_iqr", "median_no_iqr",
        "proportion_or_rate", "qualitative_with_pvalue", "threshold_or_limit",
    })
    def _is_strictly_quantitative_outcome(_o_obj: object) -> bool:
        if not isinstance(_o_obj, dict):
            return False
        _st_t = _o_obj.get("stat_type")
        if not isinstance(_st_t, str) or _st_t not in _QUANTITATIVE_CEG_STAT_TYPES:
            return False
        _qs_s = _o_obj.get("quantitative_summary")
        if not isinstance(_qs_s, str) or not _qs_s.strip():
            return False
        if _qs_s.startswith("[Structured outcome "):
            _m_st = re.search(r"\| stat=([^|\]]+)", _qs_s)
            _m_dr = re.search(r"\| direction=([^|\]]+)", _qs_s)
            _m_sh = re.search(r"\| text_sha256=([0-9a-f]{64})\]", _qs_s)
            return bool(
                _m_st and _m_st.group(1).strip() == _st_t
                and _m_dr and _m_dr.group(1).strip() == str(_o_obj.get("direction") or "")
                and _m_sh and _m_sh.group(1) != ("0" * 64)
            )
        _qs_clean = re.sub(r"(?i)\b(?:19|20)\d{2}\b|\bfig(?:ure|\.)?\s*\d+[a-z]?\b|\bno\s+(?:statistically\s+)?significant\b", "", _qs_s)
        return bool(re.search(r"(?i)\d|\b(?:seven|eight|nine|ten)\b|\b(?:p|q)\s*[<>=]", _qs_clean))
    _m206_l1 = _edges_by_eid_35.get("M0206:L1")
    if _m206_l1 is None:
        raise ValueError("Missing adjudicated legacy edge M0206:L1")
    if _m206_l1.get("primary_certification_eligible") is not True:
        raise ValueError("Adjudicated legacy edge M0206:L1 must retain primary_certification_eligible=True")
    if not isinstance(_m206_l1.get("legacy_certification_adjudication_note"), str) or len(_m206_l1["legacy_certification_adjudication_note"].strip()) < 25:
        raise ValueError("Adjudicated legacy edge M0206:L1 missing required legacy_certification_adjudication_note")
    _m206_l1_cited = set(_m206_l1.get("cited_outcome_tuples") or ())
    for _mq1_oid in ("M0206:O1", "M0206:O4", "M0206:O5", "M0206:O9", "M0206:O12", "M0206:O14"):
        if _mq1_oid not in _m206_l1_cited or not _is_strictly_quantitative_outcome(outcomes.get(_mq1_oid)):
            raise ValueError("Adjudicated eligible edge M0206:L1 lost required quantitative supporting outcomes (M0206:O1/O4/O5/O9/O12/O14)")
    _m206_quant_pins = {
        "M0206:O9": "1fb9a20a547190d51e307c1166d1fa6740f56170a1780c32dfe3afe8fd6b67a0",
        "M0206:O10": "03ebbd6bb92b56428a1794c5559276cc3309a89d7cecaab7689deb00c40b4a02",
        "M0206:O12": "4086954e70c1a349884f93c82fce65a0386a125c54f6fbffcb785588250a6634",
        "M0206:O13": "f268b00678dfa3016807ec43a9f303505c248d23d77add088063323d87eecab5",
        "M0206:O14": "15d5092f8f123c4987a331cd186296b398e552ead87417cd8aac9ab8a511ac9b",
    }
    _iip0035_reclass_qual_full_pins = {
        "M0051:O1": "3c280612321fcd0a2c1702a2e245f02715a8b4af2a22245afcda8690ca45b391",
        "M0051:O2": "e947f8ed906f1f1900fe00e44be9adaa3bd51ded06a7c6d3864d80dc27f30ede",
        "M0051:O3": "e91c31a485df12a9fdb5f4e63623b7bcac9226b13612c1db3a4f448f3472bbc0",
        "M0051:O4": "4e485c64057d1c8c6b721e55afe6f83a7e3983c424b93250ec5da6904964578c",
        "M0051:O5": "e6be27c3d21e2c23e83c56be1f8daafa2d62b8a2e481a18c32eadc751ab59e6a",
        "M0051:O6": "381e4d62f9acbd2e727d6e48e2f17847e2b4beffc71852bf91460ecbe15093cc",
        "M0154:O1": "ccc97f2925d55996613d686c019577a2516dcaa8117e09e738acae00b8127fa8",
        "M0154:O2": "d0e8b9d6db11a449654b43e432df2ebab81910a16f9b03c33697abfe0f0fd296",
        "M0154:O3": "390f6b65652ddbe8ada191f2255171f0c553a74e40c87319824be5b11d20273a",
        "M0154:O4": "ef6231cf1cb2443b08777ad8966a11390cae6c5f1aa42fd51d42879e73ea4fec",
        "M0159:O1": "71d9e24fc182d73ad25926ec59cebeaa919350e0892ae86b456a5d71c6d3e66d",
        "M0159:O2": "553d2c8471937ee648dca3d48ccb9b739078d38e913a7697900703945ad81108",
        "M0159:O3": "5ff62feeb8b62d420a14e38214b9dc49ce0217872cf7920808e3813247258ab3",
        "M0206:O16": "48a1c081cf79d9594f51247a51c4d49886e843bab4ff46a6f9d447e8e665727c",
        "M0206:O17": "05c8575a22e162ef2ce26b238485246762d3a1287e3486485e66e5d0441d0630",
        "M0241:O1": "31436d5fdd95b9180016ddfd7aa5a048034af6b17963b3967b41d25db758392c",
        "M0241:O2": "d7c3513790f8c23e57d0df00846f61f48b54222b16454d599a9b683dd2ce32f9",
        "M0241:O3": "eecf009c9fb794a349020f96af4ef19910e3bdd9db69d32540d8be881f9a5cdf",
        "M0241:O4": "b4f3b1f2529d4ec7b4170177d0ca3bd0b2cf52cb15c23638648161fcd94183b9",
    }
    _iip0035_extra_touched_legacy_pins = {
        "M0206:O1": "cf8d9f6e1a0b1d2775456eb21f408f6f4d8c00fde2a02e5a13ea65ebb55d2a8e",
        "M0206:O2": "0ce6f72a53bfb54a8250c42a68c012a2f5a4dd7fe373d913721227390a1b62e1",
        "M0206:O11": "1c45715bbc92a8526f97d331eb79c3e3916eddb383a9e164a7141e012a6d1498",
    }
    for _e_rc in edges:
        if _e_rc.get("primary_certification_eligible") is True:
            _rc_cites = list(_e_rc.get("cited_outcome_tuples") or ())
            if any(tok in _iip0035_reclass_qual_oids for tok in _rc_cites):
                if not any(_is_strictly_quantitative_outcome(outcomes.get(tok)) for tok in _rc_cites):
                    raise ValueError(f"Primary-certification-eligible edge {_e_rc.get('edge_id')} citing reclassified qualitative outcomes retains zero quantitative outcomes")
    _iip0035_expected_outcome_semantics = {
        "M0189:O1": (("VAR:M0189:A1",), ("VAR:M0189:A2_allogeneic_ctrl",), "favorable_to_intervention", "2ec0a977b228e707bc2625dbecf2e40f2cfe382374a8e468e4b8820d637f6db4", ("8", "14", "autologous")),
        "M0189:O2": (("VAR:M0189:A2_allogeneic_ctrl",), (), "unfavorable_or_toxic", "42085854453c51bb47668076b3a22186a9f3c504532f88a6399cdecddee78fc5", ("8", "weeks", "allogeneic")),
        "M1407:O1": (("VAR:M1407:A1",), (), "unfavorable_or_toxic", "43fb4a54bd27bc337ee9e68e7f1e186093d7c11f24ad9a6d03481ef0ff5c36bc", ("200", "8", "3", "5")),
        "PMID:36944335:O1": (("VAR:PMID:36944335:A1",), (), "favorable_to_intervention", "01bb7a582cf7547daeb954f55a95499f7281452043f89c1bfa9a197a65932e85", ("16", "17")),
        "PMID:39753706:O1": (("VAR:PMID:39753706:A1",), ("VAR:PMID:39753706:A2_original_protocol_ctrl",), "favorable_to_intervention", "58bb03f7d4906ac6fc056f34de4ead6a76048797c435967179b231ee97df7a29", ("20-fold", "16", "resistant")),
        "PMID:39753706:O2": (("VAR:PMID:39753706:A1",), (), "favorable_to_intervention", "a7af9f243994a87a0abef115cbc4d55bd86076cc0988b22ffe82025757c46b65", ("10", "day", "15", "100%")),
        "PMID:41436438:O1": (("VAR:PMID:41436438:A1",), (), "favorable_to_intervention", "7f41cbfdb15b08013b43f55920367b85e54073eb9bbd56069406c76bca7b7191", ("14", "100")),
        "PMID:41999746:O1": (("VAR:PMID:41999746:A1",), ("VAR:PMID:41999746:A2_p53_suppressed",), "favorable_to_intervention", "2a8b7f3d956ade3d8d9076539ed4e3f26c98cc17d7634b0e7411aa9783495856", ("p53", "impairs")),
    }
    for _sem_oid in sorted(_iip0035_expected_outcome_semantics):
        if not isinstance(outcomes.get(_sem_oid), dict):
            raise ValueError(f"Missing required IIP-0035 outcome {_sem_oid}")
    for _chk_qid, _chk_pm in sorted(_strict_study_expected_pmid.items()):
        _chk_st = studies.get(_chk_qid)
        if isinstance(_chk_st, dict):
            if _study_pmid(_chk_st) != _chk_pm or str(_chk_st.get("id_value") or "") != _chk_pm:
                raise ValueError(f"Study {_chk_qid} owner PMID/id_value mismatch: got pmid={_study_pmid(_chk_st)!r}, id_value={_chk_st.get('id_value')!r}, expected {_chk_pm!r}")
            if _chk_qid.startswith("PMID:") and _chk_qid.split(":", 1)[1] != _chk_pm:
                raise ValueError(f"Study {_chk_qid} queue_id PMID prefix disagrees with expected {_chk_pm!r}")
    _sf_led = (ceg.get("screening_funnel") or {}).get("iip0035_record_level_screening_ledger")
    if not isinstance(_sf_led, list) or len(_sf_led) != 11:
        raise ValueError("CEG screening_funnel missing 11-record iip0035_record_level_screening_ledger")
    _led_expected_rows = {
        "M0051": ("41201433", "existing_wave1_ceg_candidate_crosslinked_only_deferred_from_ring1_to_preserve_layer1_baseline", (0, 0, 0, 0, 0, 0, 0, 0)),
        "M0154": ("41690310", "existing_wave1_ceg_candidate_promoted_to_ring1_and_crosslinked", (0, 1, -1, 0, 0, 0, 0, 0)),
        "M0159": ("41764741", "existing_wave1_ceg_candidate_promoted_to_ring1_and_crosslinked", (0, 1, -1, 0, 0, 0, 0, 0)),
        "M0206": ("35953586", "existing_wave1_ceg_candidate_promoted_to_ring1_and_crosslinked", (0, 1, -1, 0, 0, 0, 0, 0)),
        "M0241": ("31427818", "existing_wave1_ceg_candidate_promoted_to_ring1_and_crosslinked", (0, 1, -1, 0, 0, 0, 0, 0)),
        "M0189": ("38843830", "existing_wave2b_queued_candidate_extracted_to_ring1_and_wave1_ceg", (0, 1, -1, 0, 0, 0, 1, -1)),
        "M1407": ("2642862", "existing_parked_context_record_promoted_to_ring1_and_wave1_ceg", (0, 1, -1, 1, -1, 0, 1, 0)),
        "PMID:36944335": ("36944335", "new_literature_screened_advanced_to_ring1_and_wave1_ceg", (1, 1, 0, 1, 0, 0, 1, 0)),
        "PMID:39753706": ("39753706", "new_literature_screened_advanced_to_ring1_and_wave1_ceg", (1, 1, 0, 1, 0, 0, 1, 0)),
        "PMID:41436438": ("41436438", "new_literature_screened_advanced_to_ring1_and_wave1_ceg", (1, 1, 0, 1, 0, 0, 1, 0)),
        "PMID:41999746": ("41999746", "new_literature_screened_advanced_to_ring1_and_wave1_ceg", (1, 1, 0, 1, 0, 0, 1, 0)),
    }
    _led_delta_keys = ("d_screened_total", "d_ring1_screened", "d_new_screened", "d_advanced", "d_parked", "d_excluded", "d_wave1_extracted", "d_remaining_catalog")
    _seen_led_qids: set[str] = set()
    for _lr in _sf_led:
        if not isinstance(_lr, dict):
            raise ValueError("CEG iip0035_record_level_screening_ledger contains non-dict row")
        _lq = _lr.get("queue_id")
        if not isinstance(_lq, str) or _lq not in _led_expected_rows or _lq in _seen_led_qids:
            raise ValueError(f"CEG screening ledger has unexpected or duplicate queue_id={_lq!r}")
        _seen_led_qids.add(_lq)
        _exp_lpm, _exp_ltr, _exp_dtup = _led_expected_rows[_lq]
        if _lq not in studies or _study_pmid(studies[_lq]) != _exp_lpm or type(_lr.get("pmid")) is not str or _lr.get("pmid") != _exp_lpm:
            raise ValueError(f"CEG screening ledger row {_lq} PMID mismatch (row_pmid={_lr.get('pmid')!r}, study_pmid={_study_pmid(studies[_lq]) if _lq in studies else None!r}, expected={_exp_lpm!r})")
        if _lq.startswith("PMID:") and _lr.get("discovery_provenance") != "targeted_iip0035_proposal_screening_not_in_frozen_wave0_discovered_pmids_all_sources":
            raise ValueError(f"CEG screening ledger row {_lq} missing required targeted discovery_provenance tag (got {_lr.get('discovery_provenance')!r})")
        if _lr.get("transition") != _exp_ltr:
            raise ValueError(f"CEG screening ledger row {_lq} transition={_lr.get('transition')!r} mismatches expected {_exp_ltr!r}")
        _got_dtup = []
        for _dk in _led_delta_keys:
            _dv = _lr.get(_dk)
            if type(_dv) is not int:
                raise ValueError(f"CEG screening ledger row {_lq} field {_dk}={_dv!r} must be strict int")
            _got_dtup.append(_dv)
        if tuple(_got_dtup) != _exp_dtup:
            raise ValueError(f"CEG screening ledger row {_lq} deltas={tuple(_got_dtup)!r} mismatch expected {_exp_dtup!r}")
        if _got_dtup[1] + _got_dtup[2] != _got_dtup[0] or _got_dtup[3] + _got_dtup[4] + _got_dtup[5] != _got_dtup[0] or _got_dtup[6] + _got_dtup[7] != _got_dtup[3]:
            raise ValueError(f"CEG screening ledger row {_lq} violates per-row conservation: {tuple(_got_dtup)!r}")
    if _seen_led_qids != set(_led_expected_rows):
        raise ValueError(f"CEG screening ledger missing expected records: {sorted(set(_led_expected_rows) - _seen_led_qids)}")
    _led_base = {"d_screened_total": 3938, "d_ring1_screened": 594, "d_new_screened": 3344, "d_advanced": 3263, "d_parked": 522, "d_excluded": 153, "d_wave1_extracted": 1176, "d_remaining_catalog": 2087}
    _led_map = {"d_screened_total": "screened_records_total", "d_ring1_screened": "ring1_existing_records_screened", "d_new_screened": "new_candidate_records_screened", "d_advanced": "candidates_advanced_total", "d_parked": "contextual_reviews_parked_total", "d_excluded": "off_topic_records_excluded_total", "d_wave1_extracted": "wave1_ceg_records_extracted_total", "d_remaining_catalog": "remaining_supporting_catalog_queued_for_later_waves"}
    _ceg_counts_35 = ceg.get("counts") if isinstance(ceg.get("counts"), dict) else {}
    for _dk, _sfk in _led_map.items():
        _calc = _led_base[_dk] + sum(r[_dk] for r in _sf_led)
        if ceg["screening_funnel"].get(_sfk) != _calc:
            raise ValueError(f"CEG screening_funnel {_sfk}={ceg['screening_funnel'].get(_sfk)} disagrees with 11-record ledger sum {_calc}")
        if _sfk in _ceg_counts_35 and _ceg_counts_35.get(_sfk) != _calc:
            raise ValueError(f"CEG counts.{_sfk}={_ceg_counts_35.get(_sfk)} disagrees with screening_funnel {_calc}")
    _arm_tok_re = re.compile(r"\b(A\d+(?:_[A-Za-z0-9_]+)?)\b")
    _span_re = re.compile(r"^([A-Za-z0-9_\.\[\]]+)\[(\d+):(\d+)\]\|norm=exact$")

    for _qid34 in sorted(_strict_new_qids):
        _st34 = studies.get(_qid34)
        if isinstance(_st34, dict) and _qid34 in _strict_study_expected_pmid:
            _exp_pm = _strict_study_expected_pmid[_qid34]
            if _study_pmid(_st34) != _exp_pm or str(_st34.get("id_value") or "") != _exp_pm:
                raise ValueError(f"Study {_qid34} owner PMID/id_value mismatch: got pmid={_study_pmid(_st34)!r}, id_value={_st34.get('id_value')!r}, expected {_exp_pm!r}")
            if _qid34.startswith("PMID:") and _qid34.split(":", 1)[1] != _exp_pm:
                raise ValueError(f"Study {_qid34} queue_id PMID prefix disagrees with expected {_exp_pm!r}")
        if not isinstance(_st34, dict):
            raise ValueError(f"Missing required IIP-0034 study {_qid34}")
        for _ef in (_st34.get("candidate_effector_topics") or []):
            if _ef not in effectors:
                raise ValueError(f"IIP-0034 study {_qid34} references unresolved effector {_ef!r}")
    for _vid34, _a34 in arms.items():
        _aq34 = _a34.get("queue_id")
        if _aq34 not in _strict_new_qids:
            continue
        _exp_vid = f"VAR:{_aq34}:{_a34.get('arm_id')}"
        if _vid34 != _exp_vid or _a34.get("variable_id") != _exp_vid:
            raise ValueError(f"IIP-0034 arm {_vid34} mismatches embedded variable_id/arm_id ({_exp_vid})")
        for _ef in (_a34.get("candidate_effector_topics") or []):
            if _ef not in effectors:
                raise ValueError(f"IIP-0034 arm {_vid34} references unresolved effector {_ef!r}")
    _legacy_m1522_eids = frozenset({"M1522:L1", "M1522:L2"})
    for _e34 in edges:
        _eq34 = _e34.get("queue_id")
        _eid34 = str(_e34.get("edge_id") or "")
        if _eid34 in _iip0035_crosslink_eids:
            if _e34.get("primary_certification_eligible") is not False:
                raise ValueError(f"IIP-0035 cross-link edge {_eid34} must keep primary_certification_eligible=False")
            _pnote = _e34.get("provenance_note")
            if (
                not isinstance(_pnote, str)
                or len(_pnote.strip()) < 40
                or not re.search(r"(?i)\b(?:reused|legacy|prior)\b", _pnote)
                or not re.search(r"(?i)\b(?:outside|non-certifying|primary_certification_eligible\s*=\s*false|scoping)\b", _pnote)
            ):
                raise ValueError(
                    f"IIP-0035 cross-link edge {_eid34} missing grounded provenance_note disclosing reused legacy outcomes and non-certifying scoping status"
                )
        if _e34.get("primary_certification_eligible") is True and _e34.get("applicability") == "unvalidated_transfer":
            if _eq34 in _strict_new_qids or _eid34 in _strict_new_eids or (_eq34 == "M1522" and _eid34 not in _legacy_m1522_eids):
                raise ValueError(
                    f"Edge {_eid34} ({_eq34}) with applicability='unvalidated_transfer' cannot be primary_certification_eligible=True"
                )
        if _eid34 in _strict_new_eids or _eq34 in _strict_new_qids:
            for _tok34 in (_e34.get("cited_outcome_tuples") or []):
                if outcomes[_tok34].get("queue_id") != _eq34:
                    raise ValueError(f"IIP-0034 edge {_eid34} cites cross-study outcome {_tok34} (owned by {outcomes[_tok34].get('queue_id')!r}, not {_eq34!r})")
            for _ef in (_e34.get("candidate_effector_topics") or []):
                if _ef not in effectors:
                    raise ValueError(f"IIP-0034 edge {_eid34} references unresolved effector {_ef!r}")

    _dossiers_dir = repo_root / "data" / "dossiers"
    _has_dossiers_dir = _dossiers_dir.is_dir()
    _loaded_dos_by_pmid: dict[str, dict[str, Any]] = {}
    if _has_dossiers_dir or require_sources:
        _present_34 = [pm for pm in _iip0034_req_pmids if (_dossiers_dir / f"{pm}.json").is_file()]
        _present_35 = [pm for pm in _iip0035_req_pmids if (_dossiers_dir / f"{pm}.json").is_file()]
        _req_cohorts: list[str] = []
        if require_sources:
            _req_cohorts = list(_strict_req_pmids)
        else:
            if _present_34:
                _req_cohorts.extend(_iip0034_req_pmids)
            if _present_35:
                _req_cohorts.extend(_iip0035_req_pmids)
            if not _req_cohorts:
                _req_cohorts = list(_strict_req_pmids)
        _missing_dos = [pm for pm in _req_cohorts if not (_dossiers_dir / f"{pm}.json").is_file()]
        if _missing_dos:
            raise ValueError(
                f"Incomplete local dossier set for source-level CEG verification: missing {_missing_dos} under {_dossiers_dir}"
            )
        for _pm34 in _req_cohorts:
            _dobj_34 = json.loads((_dossiers_dir / f"{_pm34}.json").read_text(encoding="utf-8"))
            _int_pm34 = str(_dobj_34.get("pmid") or "").strip()
            if _int_pm34 != _pm34:
                raise ValueError(f"Dossier {_pm34}.json internal pmid={_int_pm34!r} mismatches filename stem {_pm34!r}")
            if _pm34 in _iip0035_expected_dossier_digests:
                _got_dsha35 = _sha256_str(json.dumps(_dobj_34, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
                if _got_dsha35 != _iip0035_expected_dossier_digests[_pm34]:
                    raise ValueError(f"IIP-0035 dossier {_pm34}.json full-object digest {_got_dsha35} mismatches pinned {_iip0035_expected_dossier_digests[_pm34]}")
            _loaded_dos_by_pmid[_pm34] = _dobj_34

    def _resolve_dossier_span(_dos_obj: dict[str, Any], _path_expr: str) -> str:
        _cur: Any = _dos_obj
        for _part in _path_expr.split("."):
            _m_idx = re.fullmatch(r"([A-Za-z0-9_]+)\[(\d+)\]", _part)
            if _m_idx:
                _cur = _cur[_m_idx.group(1)][int(_m_idx.group(2))]
            else:
                _cur = _cur[_part]
        if not isinstance(_cur, str):
            raise ValueError(f"Dossier span path {_path_expr!r} did not resolve to string")
        return _cur

    def _verify_pin_obj(_oid_p: str, _pm_p: str, _pin_obj: dict[str, Any], _expected_sha: Any) -> None:
        _sp_str = str((_pin_obj or {}).get("span") or "")
        _m_sp = _span_re.fullmatch(_sp_str)
        if not _m_sp:
            raise ValueError(f"Outcome {_oid_p} has malformed span={_sp_str!r}")
        _p_expr, _s_i, _e_i = _m_sp.group(1), int(_m_sp.group(2)), int(_m_sp.group(3))
        if _s_i < 0 or _e_i <= _s_i:
            raise ValueError(f"Outcome {_oid_p} has invalid span coordinates [{_s_i}:{_e_i}]")
        if not isinstance(_expected_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", _expected_sha):
            raise ValueError(f"Outcome {_oid_p} missing valid 64-hex quote_sha256")
        if _pm_p in _loaded_dos_by_pmid:
            _src_txt = _resolve_dossier_span(_loaded_dos_by_pmid[_pm_p], _p_expr)
            if _e_i > len(_src_txt):
                raise ValueError(f"Outcome {_oid_p} span [{_s_i}:{_e_i}] exceeds dossier length {len(_src_txt)}")
            _slice_sha = _sha256_str(_src_txt[_s_i:_e_i])
            if _slice_sha != _expected_sha:
                raise ValueError(f"Outcome {_oid_p} span digest {_slice_sha} mismatches quote_sha256 {_expected_sha!r}")

    _verified_iip0034_span_oids: list[str] = []
    _verified_iip0035_span_oids: list[str] = []
    for oid, o in outcomes.items():
        _oqid = o.get("queue_id")
        if _oqid not in _strict_new_qids:
            continue
        _exp_oid = f"{_oqid}:{o.get('outcome_id')}"
        if oid != _exp_oid or o.get("outcome_node_id") != _exp_oid:
            raise ValueError(f"IIP-0034 outcome {oid} mismatches embedded outcome_node_id/outcome_id ({_exp_oid})")
        for _ef in (o.get("candidate_effector_topics") or []):
            if _ef not in effectors:
                raise ValueError(f"IIP-0034 outcome {oid} references unresolved effector {_ef!r}")
        _iv_arms = o.get("intervention_arm_ids")
        _cmp_arms = o.get("comparator_arm_ids")
        if not isinstance(_iv_arms, list) or not _iv_arms or not isinstance(_cmp_arms, list):
            raise ValueError(f"IIP-0034 outcome {oid} missing required typed intervention_arm_ids/comparator_arm_ids lists")
        if len(_iv_arms) != len(set(_iv_arms)) or len(_cmp_arms) != len(set(_cmp_arms)):
            raise ValueError(f"IIP-0034 outcome {oid} contains duplicate entries in intervention_arm_ids or comparator_arm_ids")
        if not set(_iv_arms).isdisjoint(set(_cmp_arms)):
            raise ValueError(f"IIP-0034 outcome {oid} intervention_arm_ids and comparator_arm_ids overlap")
        for _t_vid in (*_iv_arms, *_cmp_arms):
            if not isinstance(_t_vid, str) or _t_vid not in arms:
                raise ValueError(f"IIP-0034 outcome {oid} typed arm ref {_t_vid!r} does not resolve in CEG arms")
            if arms[_t_vid].get("queue_id") != _oqid:
                raise ValueError(f"IIP-0034 outcome {oid} typed arm ref {_t_vid!r} belongs to {arms[_t_vid].get('queue_id')!r}, not {_oqid!r}")
        _opt_cmp = o.get("comparator_arm_id")
        if _opt_cmp is not None:
            _full_opt_cmp = _opt_cmp if _opt_cmp.startswith("VAR:") else f"VAR:{_oqid}:{_opt_cmp}"
            if _full_opt_cmp not in _cmp_arms:
                raise ValueError(f"IIP-0034 outcome {oid} comparator_arm_id {_opt_cmp!r} disagrees with comparator_arm_ids={_cmp_arms!r}")
        _ca_raw = o.get("comparison_arms")
        if not isinstance(_ca_raw, str) or not _ca_raw.strip():
            raise ValueError(f"IIP-0034 outcome {oid} comparison_arms must be a non-empty string")
        _ca_tokens = _arm_tok_re.findall(_ca_raw)
        _ca_vids: set[str] = set()
        for _atok in _ca_tokens:
            _cand_vid = f"VAR:{_oqid}:{_atok}"
            if _cand_vid not in arms or arms[_cand_vid].get("queue_id") != _oqid:
                raise ValueError(f"Outcome {oid} comparison_arms references unresolved or cross-study arm {_cand_vid}")
            _ca_vids.add(_cand_vid)
        if _ca_vids != (set(_iv_arms) | set(_cmp_arms)):
            raise ValueError(
                f"IIP-0034 outcome {oid} typed arm lists ({sorted(set(_iv_arms) | set(_cmp_arms))}) disagree with comparison_arms ({sorted(_ca_vids)})"
            )
        _vs_split = re.split(r"\bvs\b", _ca_raw, maxsplit=1, flags=re.IGNORECASE)
        if len(_vs_split) == 2:
            _lhs_vids = {f"VAR:{_oqid}:{t}" for t in _arm_tok_re.findall(_vs_split[0])}
            _rhs_vids = {f"VAR:{_oqid}:{t}" for t in _arm_tok_re.findall(_vs_split[1])}
            if _lhs_vids != set(_iv_arms) or _rhs_vids != set(_cmp_arms):
                raise ValueError(
                    f"IIP-0034 outcome {oid} typed intervention/comparator roles disagree with comparison_arms LHS/RHS"
                )
        elif _cmp_arms:
            raise ValueError(f"IIP-0034 outcome {oid} lists comparator_arm_ids={_cmp_arms!r} without a 'vs' contrast in comparison_arms")
        _opt_iv = o.get("intervention_arm_id")
        if _opt_iv is not None:
            _full_opt_iv = _opt_iv if str(_opt_iv).startswith("VAR:") else f"VAR:{_oqid}:{_opt_iv}"
            if _full_opt_iv != _iv_arms[0] or _full_opt_iv in _cmp_arms:
                raise ValueError(f"Outcome {oid} singular intervention_arm_id={_opt_iv!r} disagrees with primary intervention_arm_ids={_iv_arms!r} / comparator_arm_ids={_cmp_arms!r}")
        if _oqid in _iip0035_qids:
            _fa = o.get("focal_arm_id")
            if not isinstance(_fa, str) or _fa != _iv_arms[0] or _fa in _cmp_arms:
                raise ValueError(f"IIP-0035 outcome {oid} focal_arm_id={_fa!r} must match primary intervention arm {_iv_arms[0]!r} and never appear in comparator_arm_ids={_cmp_arms!r}")
            if not isinstance(_opt_iv, str) or _opt_iv != _iv_arms[0] or _opt_iv in _cmp_arms:
                raise ValueError(f"IIP-0035 outcome {oid} required singular intervention_arm_id={_opt_iv!r} must match primary intervention arm {_iv_arms[0]!r} and never appear in comparator_arm_ids={_cmp_arms!r}")
            if oid in _iip0035_expected_outcome_semantics:
                _sem_iv, _sem_cmp, _sem_dir, _sem_csha, _sem_toks = _iip0035_expected_outcome_semantics[oid]
                if tuple(_iv_arms) != _sem_iv or tuple(_cmp_arms) != _sem_cmp:
                    raise ValueError(f"IIP-0035 outcome {oid} arm-role semantics ({_iv_arms!r}, {_cmp_arms!r}) mismatch expected ({list(_sem_iv)!r}, {list(_sem_cmp)!r})")
                if o.get("direction") != _sem_dir:
                    raise ValueError(f"IIP-0035 outcome {oid} direction={o.get('direction')!r} mismatches expected {_sem_dir!r}")
                _joined_o = " ".join(str(o.get(k) or "") for k in ("endpoint_name", "intervention_value", "quantitative_summary")).lower()
                for _tk in _sem_toks:
                    if _tk.lower() not in _joined_o:
                        raise ValueError(f"IIP-0035 outcome {oid} missing required structured endpoint/magnitude token {_tk!r}")
        _pm_o = _study_pmid(studies[_oqid])
        _exp_pm_o = _strict_study_expected_pmid.get(_oqid)
        if _exp_pm_o is not None and _pm_o != _exp_pm_o:
            raise ValueError(f"Outcome {oid} resolved PMID {_pm_o!r} mismatches expected owner {_exp_pm_o!r}")
        if _loaded_dos_by_pmid and _oqid in _iip0035_qids and all(pm in _loaded_dos_by_pmid for pm in _iip0035_req_pmids):
            if _pm_o not in _loaded_dos_by_pmid:
                raise ValueError(f"Outcome {oid} owner PMID {_pm_o!r} missing from loaded IIP-0035 dossiers (refusing silent span skip)")
        if _loaded_dos_by_pmid and _oqid in _iip0034_qids and all(pm in _loaded_dos_by_pmid for pm in _iip0034_req_pmids):
            if _pm_o not in _loaded_dos_by_pmid:
                raise ValueError(f"Outcome {oid} owner PMID {_pm_o!r} missing from loaded IIP-0034 dossiers (refusing silent span skip)")
        _sp_primary = o.get("source_pin") or {}
        if isinstance(_sp_primary, dict) and "quote_sha256" in _sp_primary and _sp_primary["quote_sha256"] != o.get("quote_sha256"):
            raise ValueError(
                f"Outcome {oid} nested source_pin.quote_sha256 {_sp_primary['quote_sha256']!r} conflicts with top-level quote_sha256 {o.get('quote_sha256')!r}"
            )
        _verify_pin_obj(oid, _pm_o, _sp_primary, o.get("quote_sha256"))
        if _pm_o in _loaded_dos_by_pmid:
            if _oqid in _iip0035_qids:
                if oid in _iip0035_expected_outcome_semantics:
                    _m_t = _span_re.fullmatch(str((_sp_primary or {}).get("span") or ""))
                    if _m_t:
                        _sl_txt = _resolve_dossier_span(_loaded_dos_by_pmid[_pm_o], _m_t.group(1))[int(_m_t.group(2)):int(_m_t.group(3))].lower()
                        for _tk in _iip0035_expected_outcome_semantics[oid][4]:
                            if _tk.lower() not in _sl_txt:
                                raise ValueError(f"IIP-0035 outcome {oid} pinned source span lacks required magnitude/endpoint token {_tk!r}")
                _verified_iip0035_span_oids.append(oid)
            elif _oqid in _iip0034_qids:
                _verified_iip0034_span_oids.append(oid)
        for _sec_pin in (o.get("secondary_source_pins") or []):
            if not isinstance(_sec_pin, dict):
                raise ValueError(f"Outcome {oid} malformed secondary_source_pins entry")
            _verify_pin_obj(oid, _pm_o, _sec_pin, _sec_pin.get("quote_sha256"))

    if all(pm in _loaded_dos_by_pmid for pm in _iip0035_req_pmids) and len(_verified_iip0035_span_oids) != 8:
        raise ValueError(f"Expected exactly 8 exact-span-verified IIP-0035 outcomes, got {len(_verified_iip0035_span_oids)}: {sorted(_verified_iip0035_span_oids)}")
    if all(pm in _loaded_dos_by_pmid for pm in _iip0034_req_pmids) and len(_verified_iip0034_span_oids) != 8:
        raise ValueError(f"Expected exactly 8 exact-span-verified IIP-0034 outcomes, got {len(_verified_iip0034_span_oids)}: {sorted(_verified_iip0034_span_oids)}")

    if _loaded_dos_by_pmid:
        _union_shingles: set[tuple[str, ...]] = set()
        _rolling40_set: set[str] = set()
        for _d34 in _loaded_dos_by_pmid.values():
            _dt34 = " ".join(_leaves(_d34))
            _ws34 = _words(_dt34)
            _union_shingles.update(tuple(_ws34[i : i + 6]) for i in range(max(0, len(_ws34) - 5)))
            _nd34 = _norm_spaces(_dt34).lower()
            _rolling40_set.update(_nd34[i : i + 40] for i in range(max(0, len(_nd34) - 39)))

        _SG_SKIP_KEYS = frozenset({
            "doi", "pmid", "pmcid", "queue_id", "context_id", "variable_id", "outcome_node_id", "edge_id",
            "arm_id", "outcome_id", "quote_sha256", "span", "source_type", "intervention_arm_ids",
            "comparator_arm_ids", "cited_outcome_ids", "cited_outcome_tuples", "cited_by_edge_ids",
            "variable_arm_ids", "candidate_effector_topics", "record_class", "fidelity_tier", "manuscript_tier",
            "direction", "stat_type", "relation", "applicability", "incorporation", "readiness",
            "causal_isolation", "id_type", "id_value", "role", "figure_id", "panel", "table_id", "row_id", "col_id",
        })

        def _scan_sg(_obj: Any, _path_lbl: str) -> None:
            if isinstance(_obj, dict):
                for _k, _v in _obj.items():
                    if _k in _SG_SKIP_KEYS and not (isinstance(_v, str) and len(_v) > 32 and _k in {"figure_id", "panel", "table_id", "row_id", "col_id"}):
                        continue
                    _scan_sg(_v, f"{_path_lbl}.{_k}")
            elif isinstance(_obj, list):
                for _i, _v in enumerate(_obj):
                    _scan_sg(_v, f"{_path_lbl}[{_i}]")
            elif isinstance(_obj, str) and _obj:
                _ws_t = _words(_obj)
                if len(_ws_t) >= 6:
                    for _i in range(len(_ws_t) - 5):
                        if tuple(_ws_t[_i : _i + 6]) in _union_shingles:
                            raise ValueError(f"CEG {_path_lbl} violates _SourceGuard (6-word union shingle)")
                _ns_t = _norm_spaces(_obj).lower()
                if len(_ns_t) >= 40:
                    for _i in range(len(_ns_t) - 39):
                        if _ns_t[_i : _i + 40] in _rolling40_set:
                            raise ValueError(f"CEG {_path_lbl} violates _SourceGuard (rolling 40-char dossier substring)")

        _sg_qids: set[str] = set()
        _sg_eids: set[str] = set()
        if all(pm in _loaded_dos_by_pmid for pm in _iip0034_req_pmids):
            _sg_qids.update(_iip0034_qids)
            _sg_eids.update(_iip0034_eids)
        if all(pm in _loaded_dos_by_pmid for pm in _iip0035_req_pmids):
            _sg_qids.update(_iip0035_qids)
            _sg_eids.update(_iip0035_eids)
        for _qid34 in sorted(_sg_qids):
            _scan_sg(studies[_qid34], f"study.{_qid34}")
        for _vid, _a in arms.items():
            if _a.get("queue_id") in _sg_qids:
                _scan_sg(_a, f"arm.{_vid}")
        for _oid, _o in outcomes.items():
            if (
                _o.get("queue_id") in _sg_qids
                or (all(pm in _loaded_dos_by_pmid for pm in _iip0034_req_pmids) and _oid in {"M1522:O6", "M1522:O7", "M1522:O15", "M1522:O16"})
                or (all(pm in _loaded_dos_by_pmid for pm in _iip0035_req_pmids) and (_oid in _iip0035_reclass_qual_full_pins or _oid in _m206_quant_pins or _oid in _iip0035_extra_touched_legacy_pins))
            ):
                _scan_sg(_o, f"outcome.{_oid}")
        for _e in edges:
            if _e.get("edge_id") in _sg_eids:
                _scan_sg(_e, f"edge.{_e.get('edge_id')}")
            elif all(pm in _loaded_dos_by_pmid for pm in _iip0035_req_pmids) and (_e.get("edge_id") in _iip0035_adjudicated_inelig_legacy or _e.get("edge_id") in {"M0206:L1", "M0206:L2"}):
                _scan_sg(_e, f"edge.{_e.get('edge_id')}")
    # Run FULL-OBJECT IIP-0035 digest checks AFTER _SourceGuard so verbatim source leaks always trip _SourceGuard first!
    for _st_qid35, _exp_stsha35 in sorted(_iip0035_expected_study_digests.items()):
        _st_o35 = studies.get(_st_qid35)
        if not isinstance(_st_o35, dict):
            raise ValueError(f"Missing required IIP-0035 study {_st_qid35}")
        _got_stsha35 = _sha256_str(json.dumps(_st_o35, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _got_stsha35 != _exp_stsha35:
            raise ValueError(f"IIP-0035 study {_st_qid35} full-object digest {_got_stsha35} mismatches pinned {_exp_stsha35}")
    for _arm_vid35, _exp_asha35 in sorted(_iip0035_expected_arm_digests.items()):
        _arm_o35 = arms.get(_arm_vid35)
        if not isinstance(_arm_o35, dict):
            raise ValueError(f"Missing required IIP-0035 arm {_arm_vid35}")
        _got_asha35 = _sha256_str(json.dumps(_arm_o35, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _got_asha35 != _exp_asha35:
            raise ValueError(f"IIP-0035 arm {_arm_vid35} full-object digest {_got_asha35} mismatches pinned {_exp_asha35}")
    for _oid35_c, _sem_t35 in sorted(_iip0035_expected_outcome_semantics.items()):
        _oo35 = outcomes[_oid35_c]
        _got_csha35 = _sha256_str(json.dumps(_oo35, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _got_csha35 != _sem_t35[3]:
            raise ValueError(f"IIP-0035 outcome {_oid35_c} full-object structured endpoint/magnitude digest {_got_csha35} mismatches pinned {_sem_t35[3]}")
    for _mq_oid, _mq_sha in sorted(_m206_quant_pins.items()):
        _mq_o = outcomes.get(_mq_oid)
        _mq_got = _sha256_str(json.dumps(_mq_o, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _mq_got != _mq_sha:
            raise ValueError(f"Touched M0206 quantitative outcome {_mq_oid} full-object digest {_mq_got} mismatches pinned {_mq_sha}")
    for _rq_oid, _rq_sha in sorted({**_iip0035_reclass_qual_full_pins, **_iip0035_extra_touched_legacy_pins}.items()):
        _rq_o = outcomes.get(_rq_oid)
        if not isinstance(_rq_o, dict):
            raise ValueError(f"Missing proposal-modified legacy outcome {_rq_oid}")
        _rq_got = _sha256_str(json.dumps(_rq_o, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _rq_got != _rq_sha:
            raise ValueError(f"Proposal-modified legacy outcome {_rq_oid} full-object digest {_rq_got} mismatches pinned {_rq_sha}")
    for _eid35, (_exp_owner35, _exp_outs35, _exp_esha35) in sorted(_iip0035_expected_edge_bindings.items()):
        _eobj35 = _edges_by_eid_35[_eid35]
        _got_esha35 = _sha256_str(json.dumps(_eobj35, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _got_esha35 != _exp_esha35:
            raise ValueError(f"IIP-0035 edge {_eid35} full-object digest {_got_esha35} mismatches pinned {_exp_esha35}")
    for _leg_eid, (_leg_outs, _leg_sha) in sorted(_iip0035_adjudicated_inelig_legacy_outcomes.items()):
        _le_o = _edges_by_eid_35[_leg_eid]
        _got_lsha = _sha256_str(json.dumps(_le_o, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
        if _got_lsha != _leg_sha:
            raise ValueError(f"Adjudicated legacy edge {_leg_eid} full-object digest {_got_lsha} mismatches pinned {_leg_sha}")
    _got_m206_l1_sha = _sha256_str(json.dumps(_m206_l1, sort_keys=True, ensure_ascii=False, separators=(",", ":")))
    if _got_m206_l1_sha != _iip0035_m206_l1_expected_sha:
        raise ValueError(f"Adjudicated legacy edge M0206:L1 full-object digest {_got_m206_l1_sha} mismatches pinned {_iip0035_m206_l1_expected_sha}")

    return {
        "studies": len(studies),
        "arms": len(arms),
        "outcomes": len(outcomes),
        "edges": len(edges),
        "effectors": len(effectors),
        "graveyard_index": len(graveyard),
    }


def verify_ceg_payload_detailed(
    ceg: dict[str, Any],
    repo_root: Path = REPO_ROOT,
    *,
    require_sources: bool = False,
) -> dict[str, Any]:
    """Run verify_ceg_payload and return explicit verification_mode alongside counts."""
    counts = verify_ceg_payload(ceg, repo_root=repo_root, require_sources=require_sources)
    _dos_dir = repo_root / "data" / "dossiers"
    _p34 = ["17375075", "39696686", "42020736", "42056515", "42749803"]
    _p35 = ["2642862", "31427818", "35953586", "36944335", "38843830", "39753706", "41201433", "41436438", "41690310", "41764741", "41999746"]
    _q35 = frozenset({"M0189", "M1407", "PMID:36944335", "PMID:39753706", "PMID:41436438", "PMID:41999746"})
    has_src_34 = _dos_dir.is_dir() and all((_dos_dir / f"{pm}.json").is_file() for pm in _p34)
    has_src_35 = _dos_dir.is_dir() and all((_dos_dir / f"{pm}.json").is_file() for pm in _p35)
    _verified_35_oids = (
        sorted(oid for oid, o in (ceg.get("outcomes") or {}).items() if isinstance(o, dict) and o.get("queue_id") in _q35)
        if has_src_35
        else []
    )
    _mode = (
        "iip0035_11dossiers_8new_exact_spans_verified"
        if (has_src_35 and len(_verified_35_oids) == 8)
        else ("iip0034_5dossiers_8new_exact_spans_verified" if has_src_34 else "structural_public_only")
    )
    _e35_ids = sorted([
        "M0051:L4", "M0154:L3", "M0159:L3", "M0189:L1", "M0206:L3", "M0241:L3",
        "M1407:L1", "PMID:36944335:L1", "PMID:39753706:L1", "PMID:41436438:L1", "PMID:41999746:L1",
    ])
    return {
        "counts": counts,
        "verification_mode": _mode,
        "iip0034_dossiers_verified": _p34 if has_src_34 else [],
        "iip0035_dossiers_verified": _p35 if has_src_35 else [],
        "iip0035_verified_exact_span_outcome_ids": _verified_35_oids,
        "iip0035_verified_exact_span_outcome_count": len(_verified_35_oids),
        "iip0035_sourceguard_executed": bool(has_src_35),
        "iip0035_sourceguard_scanned_study_ids": sorted(_q35) if has_src_35 else [],
        "iip0035_sourceguard_scanned_edge_ids": _e35_ids if has_src_35 else [],
        "iip0035_sourceguard_scanned_touched_legacy_outcome_ids": (
            sorted([
                "M0051:O1", "M0051:O2", "M0051:O3", "M0051:O4", "M0051:O5", "M0051:O6",
                "M0154:O1", "M0154:O2", "M0154:O3", "M0154:O4",
                "M0159:O1", "M0159:O2", "M0159:O3", "M0206:O1", "M0206:O2", "M0206:O9", "M0206:O10",
                "M0206:O11", "M0206:O12", "M0206:O13", "M0206:O14", "M0206:O16", "M0206:O17",
                "M0241:O1", "M0241:O2", "M0241:O3", "M0241:O4",
            ])
            if has_src_35
            else []
        ),
        "iip0035_sourceguard_scanned_adjudicated_legacy_edge_ids": (
            sorted(["M0051:L1", "M0051:L2", "M0051:L3", "M0154:L1", "M0154:L2", "M0159:L1", "M0159:L2", "M0206:L1", "M0206:L2", "M0241:L1", "M0241:L2"])
            if has_src_35
            else []
        ),
        "iip0035_sourceguard_scan_scope": (
            "executed_on_6_new_study_families_11_new_edges_27_touched_legacy_outcomes_and_11_adjudicated_legacy_edges"
            if has_src_35
            else "skipped_public_bundle_no_local_dossiers"
        ),
        "reused_m1522_outcomes_exact_span_verified": False,
        "reused_legacy_crosslink_outcomes_exact_span_verified": False,
        "izadifar_dossier_contains_figure_legends": False,
        "m1522_reused_outcomes_note": "M1522:O6/O7/O15/O16 use figure/panel + digest pins outside the exact-span pass (Izadifar dossier carries no figure legends) and remain ineligible scoping references on M1522:L4.",
        "iip0035_reused_crosslink_outcomes_note": "Reused legacy cross-linked outcomes on M0051:L4, M0154:L3, M0159:L3, M0206:L3, and M0241:L3 carry prior CEG source pins outside the IIP-0035 8-new-outcome exact-span pass and remain primary_certification_eligible=False scoping qualifications only.",
    }


def import_and_materialize(
    raw_graph_path: Path,
    screening_ledger_path: Path,
    source_packets_dir: Path,
    repo_root: Path = REPO_ROOT,
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw = json.loads(raw_graph_path.read_text(encoding="utf-8"))
    ledger = json.loads(screening_ledger_path.read_text(encoding="utf-8"))

    study_shingles: dict[str, set[tuple[str, ...]]] = {}
    study_years: dict[str, int] = {}
    for pkt_file in sorted(source_packets_dir.glob("*packets*.json")):
        data = json.loads(pkt_file.read_text(encoding="utf-8"))
        if not isinstance(data, list):
            continue
        for pkt in data:
            if not isinstance(pkt, dict) or "queue_id" not in pkt:
                continue
            qid = str(pkt["queue_id"])
            yr = pkt.get("year")
            if isinstance(yr, int) and not isinstance(yr, bool) and 1900 <= yr <= 2100:
                study_years[qid] = yr
            ws = _words(" ".join(_leaves(pkt)))
            study_shingles.setdefault(qid, set()).update(
                tuple(ws[i : i + 6]) for i in range(max(0, len(ws) - 5))
            )

    pmid_to_qid = {
        str(s.get("pmid")): qid
        for qid, s in raw["study_contexts"].items()
        if s.get("pmid")
    }
    dossier_texts: list[str] = []
    dossier_dir = repo_root / "data" / "dossiers"
    if dossier_dir.is_dir():
        for fp in sorted(dossier_dir.glob("*.json")):
            d_obj = json.loads(fp.read_text(encoding="utf-8"))
            d_txt = " ".join(_leaves(d_obj))
            dossier_texts.append(d_txt)
            qid = pmid_to_qid.get(fp.stem)
            if qid:
                ws = _words(d_txt)
                study_shingles.setdefault(qid, set()).update(
                    tuple(ws[i : i + 6]) for i in range(max(0, len(ws) - 5))
                )

    guard = _SourceGuard(study_shingles, dossier_texts)

    def clean_text(qid: str, field: str, value: Any, meta: dict[str, Any]) -> Any:
        if not isinstance(value, str):
            return value
        text = _scrub_progeny_and_internal_tokens(value)
        if not text:
            return text
        if guard.overlaps(qid, value) or guard.overlaps(qid, text):
            return _structured_replacement(field, text, meta)
        return text

    studies: dict[str, Any] = {}
    for qid in sorted(raw["study_contexts"]):
        s = raw["study_contexts"][qid]
        meta = {"queue_id": qid}
        st: dict[str, Any] = {
            "context_id": s["context_id"],
            "queue_id": qid,
            "id_type": s["id_type"],
            "id_value": s["id_value"],
            "fidelity_tier": s.get("fidelity_tier"),
            "biological_model": clean_text(qid, "biological_model", s.get("biological_model"), meta),
            "immune_context": clean_text(qid, "immune_context", s.get("immune_context"), meta),
            "donor_or_cell_source": clean_text(qid, "donor_or_cell_source", s.get("donor_or_cell_source"), meta),
            "follow_up": clean_text(qid, "follow_up", s.get("follow_up"), meta),
            "variable_arm_ids": list(s.get("variable_arm_ids") or []),
            "blueprint_targets": [clean_text(qid, "blueprint_target", bt, meta) for bt in (s.get("blueprint_targets") or [])],
            "outcomes_count": int(s.get("outcomes_count") or 0),
            "links_count": int(s.get("links_count") or 0),
        }
        if qid in study_years:
            st["pub_year"] = study_years[qid]
        elif isinstance(s.get("pub_year"), int) and not isinstance(s.get("pub_year"), bool):
            st["pub_year"] = s["pub_year"]
        for opt in (
            "pmid",
            "pmcid",
            "doi",
            "superseded_by_queue_id",
            "superseded_by_pmid",
            "publication_type",
            "same_patient_followup_of_queue_id",
            "same_patient_followup_of_pmid",
            "same_patient_linked_queue_id",
            "correspondence_pin_status",
        ):
            if s.get(opt):
                st[opt] = s[opt]
        if "has_experimental_subcomponent" in s and s["has_experimental_subcomponent"] is not None:
            st["has_experimental_subcomponent"] = bool(s["has_experimental_subcomponent"])
        if s.get("experimental_subcomponent_note"):
            st["experimental_subcomponent_note"] = clean_text(
                qid, "experimental_subcomponent_note", s.get("experimental_subcomponent_note"), meta
            )
        if s.get("tier_qualification_note"):
            st["tier_qualification_note"] = clean_text(
                qid, "tier_qualification_note", s.get("tier_qualification_note"), meta
            )
        if s.get("candidate_effector_topics"):
            st["candidate_effector_topics"] = list(s["candidate_effector_topics"])
        st_pin = _clean_pin(s.get("source_pin"))
        if st_pin:
            st["source_pin"] = st_pin
        _raw_rc = s.get("record_class")
        if _raw_rc is not None and _raw_rc not in VALID_CEG_RECORD_CLASSES:
            raise ValueError(f"Raw study {qid} has invalid record_class={_raw_rc!r}")
        _pre_rc = classify_ceg_study_record(s)
        _post_rc = classify_ceg_study_record(st)
        st["record_class"] = (
            _pre_rc
            if _pre_rc != "empirical_primary_study"
            else (_raw_rc if _raw_rc in VALID_CEG_RECORD_CLASSES and _raw_rc != "empirical_primary_study" else _post_rc)
        )
        studies[qid] = st

    arms: dict[str, Any] = {}
    for vid in sorted(raw["variable_nodes"]):
        v = raw["variable_nodes"][vid]
        qid = v["queue_id"]
        meta = {
            "queue_id": qid,
            "arm_id": v["arm_id"],
            "causal_isolation": v.get("causal_isolation") or "NOT_ESTABLISHED",
        }
        ar: dict[str, Any] = {
            "variable_id": vid,
            "queue_id": qid,
            "context_id": v["context_id"],
            "arm_id": v["arm_id"],
            "label": clean_text(qid, "label", v.get("label"), meta),
            "edits_or_interventions": [
                clean_text(qid, "edits_or_interventions", ed, meta)
                for ed in (v.get("edits_or_interventions") or [])
            ],
            "n_donors_or_animals": clean_text(qid, "n_donors_or_animals", v.get("n_donors_or_animals"), meta),
            "causal_isolation": v.get("causal_isolation"),
        }
        if v.get("candidate_effector_topics"):
            ar["candidate_effector_topics"] = list(v["candidate_effector_topics"])
        arms[vid] = ar

    outcomes: dict[str, Any] = {}
    for oid in sorted(raw["outcome_nodes"]):
        o = raw["outcome_nodes"][oid]
        qid = o["queue_id"]
        pin = _clean_pin(o.get("source_pin"))
        vq = str(o.get("verbatim_quote") or "")
        meta = {
            "queue_id": qid,
            "outcome_id": o["outcome_id"],
            "direction": o.get("direction") or "unspecified",
            "stat_type": o.get("stat_type") or "unspecified",
            "pin_tag": _pin_tag(pin),
        }
        reach = o.get("outcome_reachability") or {}
        cited_edges = list(reach.get("cited_by_blueprint_link_edge_ids") or o.get("cited_by_edge_ids") or [])
        oc: dict[str, Any] = {
            "outcome_node_id": oid,
            "queue_id": qid,
            "outcome_id": o["outcome_id"],
            "endpoint_name": clean_text(qid, "endpoint_name", o.get("endpoint_name"), meta),
            "comparison_arms": clean_text(qid, "comparison_arms", o.get("comparison_arms"), meta),
            "direction": o.get("direction"),
            "stat_type": o.get("stat_type"),
            "quantitative_summary": clean_text(qid, "quantitative_summary", o.get("quantitative_summary"), meta),
            "manuscript_tier": o.get("manuscript_tier"),
        }
        for val_key in ("intervention_value", "comparator_value", "p_value_or_ci", "n_summary"):
            raw_v = o.get(val_key)
            if raw_v is not None and str(raw_v).strip() and str(raw_v).strip().lower() != "not reported":
                oc[val_key] = clean_text(qid, val_key, str(raw_v), meta)
        if o.get("requires_figure_digitization"):
            oc["requires_figure_digitization"] = True
        if pin:
            oc["source_pin"] = pin
        if vq:
            oc["quote_sha256"] = _sha256_str(vq)
        elif o.get("quote_sha256"):
            oc["quote_sha256"] = str(o["quote_sha256"])
        if o.get("evidence_provenance"):
            oc["evidence_provenance"] = str(o["evidence_provenance"])
        if cited_edges:
            oc["cited_by_edge_ids"] = cited_edges
        if o.get("candidate_effector_topics"):
            oc["candidate_effector_topics"] = list(o["candidate_effector_topics"])
        if isinstance(o.get("outcome_context_override"), dict):
            oco_in = o["outcome_context_override"]
            oc["outcome_context_override"] = {
                "study_manuscript_fidelity_tier": oco_in.get("study_manuscript_fidelity_tier"),
                "outcome_effective_assay_tier": oco_in.get("outcome_effective_assay_tier"),
                "override_reason": clean_text(qid, "override_reason", oco_in.get("override_reason"), meta),
                "override_provenance": _scrub_internal_code(str(oco_in.get("override_provenance") or "")),
            }
        if isinstance(o.get("comparator_audit_override"), dict):
            cao_in = o["comparator_audit_override"]
            keep_cao_keys = (
                "corrected_comparison_arms",
                "focal_intervention_arm",
                "misbound_master_arm_id",
                "baseline_comparator_arm",
                "assay_context",
                "audit_note",
            )
            oc["comparator_audit_override"] = {
                k: clean_text(qid, k, cao_in[k], meta)
                for k in keep_cao_keys
                if k in cao_in and cao_in[k] is not None
            }
        if isinstance(o.get("focal_arm_id"), str) and o["focal_arm_id"].strip():
            oc["focal_arm_id"] = str(o["focal_arm_id"]).strip()
        if isinstance(o.get("intervention_arm_id"), str) and o["intervention_arm_id"].strip():
            oc["intervention_arm_id"] = str(o["intervention_arm_id"]).strip()
        if isinstance(o.get("intervention_arm_ids"), list):
            oc["intervention_arm_ids"] = list(o["intervention_arm_ids"])
        if isinstance(o.get("comparator_arm_ids"), list):
            oc["comparator_arm_ids"] = list(o["comparator_arm_ids"])
        if isinstance(o.get("secondary_source_pins"), list):
            oc["secondary_source_pins"] = [dict(sp) for sp in o["secondary_source_pins"] if isinstance(sp, dict)]
        outcomes[oid] = oc

    edges: list[dict[str, Any]] = []
    for e in sorted(raw["edges"], key=lambda item: (item["queue_id"], item["link_index"])):
        qid = e["queue_id"]
        meta = {
            "queue_id": qid,
            "edge_id": e["edge_id"],
            "relation": e.get("relation") or "unspecified",
            "applicability": e.get("applicability") or "unspecified",
            "incorporation": e.get("incorporation") or "unspecified",
            "cited_outcome_ids": list(e.get("cited_outcome_ids") or []),
        }
        cited_tuples = [
            f"{t[0]}:{t[1]}" if isinstance(t, (list, tuple)) and len(t) == 2 else str(t)
            for t in (e.get("cited_outcome_tuples") or [])
        ]
        _raw_pce = e.get("primary_certification_eligible")
        if type(_raw_pce) is not bool:
            raise ValueError(
                f"Raw edge {e.get('edge_id')} primary_certification_eligible must be strict bool, got {_raw_pce!r} ({type(_raw_pce).__name__})"
            )
        _edge_prov = e.get("evidence_provenance")
        _imp_non_emp_prov = frozenset({"computational_model_or_simulation", "mixed_measured_and_simulated_model", "mixed_empirical_and_computational"})
        _imp_title_pins = frozenset({"title_or_citation_metadata_only", "title"})
        _st_pin_imp_t = (studies[qid].get("source_pin") or {}).get("source_type") if isinstance(studies[qid].get("source_pin"), dict) else None
        _cites_nonemp_or_sim = any(
            studies.get(outcomes[tok]["queue_id"], {}).get("record_class") != "empirical_primary_study"
            or outcomes[tok].get("evidence_provenance") in _imp_non_emp_prov
            or (isinstance(outcomes[tok].get("source_pin"), dict) and outcomes[tok]["source_pin"].get("source_type") in _imp_title_pins)
            or (
                isinstance(studies.get(outcomes[tok]["queue_id"], {}).get("source_pin"), dict)
                and studies[outcomes[tok]["queue_id"]]["source_pin"].get("source_type") in _imp_title_pins
            )
            for tok in cited_tuples
            if tok in outcomes
        )
        ed: dict[str, Any] = {
            "edge_id": e["edge_id"],
            "queue_id": qid,
            "link_index": int(e["link_index"]),
            "target_gate_or_module": clean_text(qid, "target_gate_or_module", e.get("target_gate_or_module"), meta),
            "relation": e.get("relation"),
            "applicability": e.get("applicability"),
            "incorporation": e.get("incorporation"),
            "readiness": e.get("readiness"),
            "primary_certification_eligible": (
                False
                if (
                    studies[qid].get("record_class") != "empirical_primary_study"
                    or e.get("applicability") == "domain_mismatch"
                    or e.get("relation") == "noninformative"
                    or e.get("readiness") == "needs_primary_fulltext"
                    or _st_pin_imp_t in _imp_title_pins
                    or _edge_prov in _imp_non_emp_prov
                    or e["edge_id"] in CURATED_INELIGIBLE_EDGE_IDS
                    or (
                        e.get("applicability") == "unvalidated_transfer"
                        and (
                            qid in {"PMID:17375075", "PMID:42020736", "PMID:42056515", "PMID:42749803"}
                            or (qid == "M1522" and e["edge_id"] not in {"M1522:L1", "M1522:L2"})
                        )
                    )
                    or not cited_tuples
                    or _cites_nonemp_or_sim
                )
                else _raw_pce
            ),
            "cited_outcome_ids": list(e.get("cited_outcome_ids") or []),
            "cited_outcome_tuples": cited_tuples,
            "claim_summary": clean_text(qid, "claim_summary", e.get("claim_summary"), meta),
        }
        if _edge_prov:
            ed["evidence_provenance"] = str(_edge_prov)
        if e.get("provenance_note"):
            ed["provenance_note"] = clean_text(qid, "provenance_note", e.get("provenance_note"), meta)
        if e.get("legacy_certification_adjudication_note"):
            ed["legacy_certification_adjudication_note"] = clean_text(qid, "legacy_certification_adjudication_note", e.get("legacy_certification_adjudication_note"), meta)
        if e.get("candidate_effector_topics"):
            ed["candidate_effector_topics"] = list(e["candidate_effector_topics"])
        edges.append(ed)

    effectors: dict[str, Any] = {}
    for eff_id in sorted(raw["effector_nodes"]):
        eff = raw["effector_nodes"][eff_id]
        effectors[eff_id] = {
            "effector_id": eff["effector_id"],
            "topic_id": eff["topic_id"],
            "label": eff["label"],
            "category": eff["category"],
            "causal_isolation": eff["causal_isolation"],
            "mediation_status": eff["mediation_status"],
            "semantic_status": eff.get("semantic_status"),
            "linked_study_ids": sorted(eff.get("linked_study_ids") or []),
            "linked_variable_ids": sorted(eff.get("linked_variable_ids") or []),
            "linked_outcome_tuples": [
                f"{t[0]}:{t[1]}" if isinstance(t, (list, tuple)) and len(t) == 2 else str(t)
                for t in (eff.get("linked_outcome_tuples") or [])
            ],
        }

    # Build graveyard_index directly from the canonical public graveyard file for zero source-prose drift
    gy_src = json.loads((repo_root / GRAVEYARD_REL_PATH).read_text(encoding="utf-8"))["graveyard_entries"]
    gy_by_id = {row["graveyard_id"]: row for row in gy_src}
    graveyard_index: list[dict[str, Any]] = []
    for g in raw["graveyard_index"]:
        gid = g["graveyard_id"]
        canon = gy_by_id.get(gid, {})
        graveyard_index.append({
            "graveyard_id": gid,
            "category": g.get("Category") or canon.get("severity") or "BLOCK",
            "severity": canon.get("severity"),
        })

    l_sum = ledger.get("summary")
    if not isinstance(l_sum, dict) or "total_studies" not in l_sum:
        raise ValueError("screening_ledger.summary must be a dict containing explicit total_studies and overall_counts")
    _oc_counts = l_sum.get("overall_counts")
    _wave_counts = l_sum.get("candidate_extraction_waves")
    if not isinstance(_oc_counts, dict) or not all(k in _oc_counts for k in ("ADVANCE", "PARK_CONTEXT", "EXCLUDE_OFF_TOPIC")):
        raise ValueError("screening_ledger.summary.overall_counts must specify ADVANCE, PARK_CONTEXT, and EXCLUDE_OFF_TOPIC")
    if not isinstance(_wave_counts, dict) or "Wave_2B_Supporting_Catalog_P2_P3" not in _wave_counts:
        raise ValueError("screening_ledger.summary.candidate_extraction_waves must specify Wave_2B_Supporting_Catalog_P2_P3")
    _tot_screened = l_sum["total_studies"]
    _r1_screened = l_sum.get("ring1_existing_records_screened", 604)
    _new_screened = l_sum.get("new_candidate_records_screened", 3338)
    _adv_total = _oc_counts["ADVANCE"]
    _park_total = _oc_counts["PARK_CONTEXT"]
    _excl_total = _oc_counts["EXCLUDE_OFF_TOPIC"]
    _rem_total = _wave_counts["Wave_2B_Supporting_Catalog_P2_P3"]
    if not all(type(v) is int and v >= 0 for v in (_tot_screened, _r1_screened, _new_screened, _adv_total, _park_total, _excl_total, _rem_total)):
        raise ValueError("screening_ledger counts must be non-negative integers")
    if _r1_screened + _new_screened != _tot_screened:
        raise ValueError(f"Inconsistent screening_ledger totals: ring1 ({_r1_screened}) + new ({_new_screened}) != total_studies ({_tot_screened})")
    if _adv_total + _park_total + _excl_total != _tot_screened:
        raise ValueError(f"Inconsistent screening_ledger totals: ADVANCE ({_adv_total}) + PARK_CONTEXT ({_park_total}) + EXCLUDE_OFF_TOPIC ({_excl_total}) != total_studies ({_tot_screened})")
    if len(studies) + _rem_total != _adv_total:
        raise ValueError(f"Inconsistent screening_ledger totals: extracted ({len(studies)}) + queued ({_rem_total}) != ADVANCE ({_adv_total})")
    rc_import = Counter(classify_ceg_study_record(s) for s in studies.values())
    screening_funnel = {
        "screened_records_total": _tot_screened,
        "ring1_existing_records_screened": _r1_screened,
        "new_candidate_records_screened": _new_screened,
        "candidates_advanced_total": _adv_total,
        "contextual_reviews_parked_total": _park_total,
        "off_topic_records_excluded_total": _excl_total,
        "wave1_ceg_records_extracted_total": len(studies),
        "wave1_primary_empirical_extracted_total": rc_import["empirical_primary_study"],
        "wave1_trial_registry_protocols_extracted": rc_import["trial_registry_protocol"],
        "wave1_regulatory_or_guidance_records_extracted": rc_import["regulatory_or_guidance_document"],
        "wave1_computational_models_extracted": rc_import["computational_model_or_simulation"],
        "wave1_correspondence_or_secondary_extracted": rc_import["correspondence_or_secondary_summary"],
        "wave1_superseded_preprint_duplicates_extracted": rc_import["superseded_preprint_duplicate"],
        "wave1_ring1_core_records_extracted_total": sum(1 for q in studies if q.startswith("R1_")),
        "wave1_ring1_empirical_core_extracted": sum(1 for q, s in studies.items() if q.startswith("R1_") and classify_ceg_study_record(s) == "empirical_primary_study"),
        "wave1_new_empirical_studies_extracted": sum(1 for q, s in studies.items() if not q.startswith("R1_") and classify_ceg_study_record(s) == "empirical_primary_study"),
        "wave1_non_r1_records_extracted_total": sum(1 for q in studies if not q.startswith("R1_")),
        "remaining_supporting_catalog_queued_for_later_waves": _rem_total,
    }
    _raw_sf_led = (raw.get("screening_funnel") or {}).get("iip0035_record_level_screening_ledger")
    _led_sf_led = l_sum.get("iip0035_record_level_screening_ledger")
    if isinstance(_led_sf_led, list) and len(_led_sf_led) == 11:
        screening_funnel["iip0035_record_level_screening_ledger"] = list(_led_sf_led)
    elif isinstance(_raw_sf_led, list) and len(_raw_sf_led) == 11:
        screening_funnel["iip0035_record_level_screening_ledger"] = list(_raw_sf_led)

    ceg_payload: dict[str, Any] = {
        "schema_version": "openislet_context_evidence_graph_v1.0",
        "graph_name": "OpenIslet T1D Context-Evidence Graph",
        "graph_acronym": "CEG",
        "license": "CC-BY-4.0",
        "quote_and_rights_policy": {
            "verbatim_quote_status": "WITHHELD_IN_PUBLIC_CEG_SEE_QUOTE_SHA256_AND_SOURCE_PIN",
            "non_oa_prose_overlap_policy": (
                "All inline verbatim_quote strings are withheld (preserving 64-hex quote_sha256 and source_pin), "
                "and any string field sharing a >=6-word normalized source shingle or >=40-character local dossier "
                "substring is replaced with a deterministic non-expressive structured descriptor."
            ),
            "omitted_outcome_field_defaults": {
                "verbatim_quote": None,
                "intervention_value": "not reported",
                "comparator_value": "not reported",
                "p_value_or_ci": "not reported",
                "n_summary": "not reported",
                "requires_figure_digitization": False,
                "cited_by_edge_ids": [],
            },
        },
        "evidence_hierarchy_and_scoring_note": (
            "Context-indexed evidence graph across 1,182 total CEG records (1,156 empirical_primary_study + 11 trial_registry_protocol + 7 regulatory_or_guidance_document + 4 computational_model_or_simulation + 1 correspondence_or_secondary_summary + 3 superseded_preprint_duplicate), 7,749 arms, 14,175 outcomes, "
            "and 2,494 blueprint links. Layer-1 2,160-stack prior-anchored scenario scores remain computed on the "
            "604-study Ring-1 baseline (+0.00 pp under IIP-0035) unless explicitly updated through an accepted Islet Improvement Proposal."
        ),
        "counts": {
            "screened_records_total": screening_funnel["screened_records_total"],
            "candidates_advanced_total": screening_funnel["candidates_advanced_total"],
            "studies": len(studies),
            "record_class_counts": dict(sorted(rc_import.items())),
            "arms": len(arms),
            "outcomes": len(outcomes),
            "edges": len(edges),
            "effectors": len(effectors),
            "graveyard_entries_indexed": len(graveyard_index),
            "unresolved_exact_outcome_mentions": 0,
        },
        "screening_funnel": screening_funnel,
        "studies": studies,
        "arms": arms,
        "outcomes": outcomes,
        "edges": edges,
        "effectors": effectors,
        "graveyard_index": graveyard_index,
    }

    verify_ceg_payload(ceg_payload, repo_root=repo_root)
    atlas_payload = build_layer1_atlas_from_ceg(ceg_payload)
    return ceg_payload, atlas_payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Build or verify OpenIslet Layer-1 Context-Evidence Graph artifacts.")
    parser.add_argument("--root", type=Path, default=REPO_ROOT, help="OpenIslet repository root")
    parser.add_argument("--verify", action="store_true", help="Verify existing artifacts without modifying files")
    parser.add_argument("--import-raw-graph", type=Path, default=None, help="Optional raw CEG JSON to import")
    parser.add_argument("--import-screening-ledger", type=Path, default=None, help="Optional screening ledger JSON")
    parser.add_argument("--source-packets-dir", type=Path, default=None, help="Optional directory of source packets for shingle scrubbing")
    args = parser.parse_args()

    root = args.root.resolve()
    ceg_path = root / CEG_REL_PATH
    atlas_path = root / ATLAS_REL_PATH

    if args.import_raw_graph is not None:
        if args.import_screening_ledger is None or args.source_packets_dir is None:
            raise ValueError("--import-raw-graph requires --import-screening-ledger and --source-packets-dir")
        ceg_payload, atlas_payload = import_and_materialize(
            args.import_raw_graph.resolve(),
            args.import_screening_ledger.resolve(),
            args.source_packets_dir.resolve(),
            repo_root=root,
        )
        ceg_bytes = (json.dumps(ceg_payload, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n").encode("utf-8")
        atlas_bytes = (json.dumps(atlas_payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
        ceg_path.parent.mkdir(parents=True, exist_ok=True)
        ceg_path.write_bytes(ceg_bytes)
        atlas_path.write_bytes(atlas_bytes)
        print(
            f"Materialized {CEG_REL_PATH} ({len(ceg_bytes) / 1024 / 1024:.2f} MB) "
            f"and {ATLAS_REL_PATH} ({len(atlas_bytes) / 1024 / 1024:.2f} MB)."
        )
        return

    if not ceg_path.is_file():
        raise FileNotFoundError(f"Missing {ceg_path}")
    import math as _math
    def _bceg_no_dup(pairs):
        out = {}
        for k, v in pairs:
            if k in out:
                raise ValueError(f"Duplicate JSON key in CEG: {k!r}")
            out[k] = v
        return out
    def _bceg_no_nan(tok):
        raise ValueError(f"Non-standard JSON constant in CEG: {tok!r}")
    def _bceg_finite_float(tok):
        v = float(tok)
        if not _math.isfinite(v):
            raise ValueError(f"Non-finite JSON float overflow in CEG: {tok!r}")
        return v
    ceg_payload = json.loads(ceg_path.read_text(encoding="utf-8"), object_pairs_hook=_bceg_no_dup, parse_constant=_bceg_no_nan, parse_float=_bceg_finite_float)
    detail = verify_ceg_payload_detailed(ceg_payload, repo_root=root)
    counts = detail["counts"]
    expected_atlas = (json.dumps(build_layer1_atlas_from_ceg(ceg_payload), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if args.verify:
        if not atlas_path.is_file() or atlas_path.read_bytes() != expected_atlas:
            raise ValueError(f"{ATLAS_REL_PATH} is missing or out of sync with {CEG_REL_PATH}")
        print(f"Verified CEG ({counts['studies']} studies, {counts['outcomes']} outcomes, {counts['edges']} links; mode={detail['verification_mode']}, exact_spans_verified={detail['iip0035_verified_exact_span_outcome_count']}, sourceguard_executed={detail['iip0035_sourceguard_executed']}) and Layer-1 atlas.")
    else:
        atlas_path.write_bytes(expected_atlas)
        print(f"Verified CEG and wrote {ATLAS_REL_PATH} ({len(expected_atlas) / 1024 / 1024:.2f} MB).")


if __name__ == "__main__":
    main()
