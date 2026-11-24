#!/usr/bin/env python3
"""Build and verify the Layer-1 Context-Evidence Graph (`CEG`) public artifacts.

Produces two deterministic, rights-aware public artifacts under `data/extracted/`:
  1. `data/extracted/context_evidence_graph_v1.json`
     Full structured T1D Context-Evidence Graph across 1,172 total CEG records
     (1,146 empirical_primary_study + 11 trial_registry_protocol + 7 regulatory_or_guidance_document
     + 4 computational_model_or_simulation + 1 correspondence_or_secondary_summary + 3 superseded_preprint_duplicate),
     7,725 experimental and control arms, 14,159 structured outcomes,
     2,478 blueprint claim/contradiction/scoping links, 20 biological effector
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
            "The Layer-1 2,160-stack prior-anchored scenario evaluator retains its 589-study Ring-1 scoring baseline "
            "unchanged (+0.00 pp under IIP-0033) unless an explicit Islet Improvement Proposal "
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
            "ring1_scoring_baseline_studies_unchanged": 589,
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


def verify_ceg_payload(ceg: dict[str, Any], repo_root: Path = REPO_ROOT) -> dict[str, int]:
    """Fail-closed structural, referential, and rights verification of `context_evidence_graph_v1.json`."""
    studies = ceg.get("studies")
    arms = ceg.get("arms")
    outcomes = ceg.get("outcomes")
    edges = ceg.get("edges")
    effectors = ceg.get("effectors")
    graveyard = ceg.get("graveyard_index")

    if not isinstance(studies, dict) or len(studies) != 1172:
        raise ValueError(f"Expected 1,172 studies in CEG, got {len(studies) if isinstance(studies, dict) else type(studies)}")
    if not isinstance(arms, dict) or len(arms) != 7725:
        raise ValueError(f"Expected 7,725 arms in CEG, got {len(arms) if isinstance(arms, dict) else type(arms)}")
    if not isinstance(outcomes, dict) or len(outcomes) != 14159:
        raise ValueError(f"Expected 14,159 outcomes in CEG, got {len(outcomes) if isinstance(outcomes, dict) else type(outcomes)}")
    if not isinstance(edges, list) or len(edges) != 2478:
        raise ValueError(f"Expected 2,478 edges in CEG, got {len(edges) if isinstance(edges, list) else type(edges)}")
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

    return {
        "studies": len(studies),
        "arms": len(arms),
        "outcomes": len(outcomes),
        "edges": len(edges),
        "effectors": len(effectors),
        "graveyard_index": len(graveyard),
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
    _r1_screened = l_sum.get("ring1_existing_records_screened", 589)
    _new_screened = l_sum.get("new_candidate_records_screened", 3345)
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
            "Context-indexed evidence graph across 1,172 total CEG records (1,146 empirical_primary_study + 11 trial_registry_protocol + 7 regulatory_or_guidance_document + 4 computational_model_or_simulation + 1 correspondence_or_secondary_summary + 3 superseded_preprint_duplicate), 7,725 arms, 14,159 outcomes, "
            "and 2,478 blueprint links. Layer-1 2,160-stack prior-anchored scenario scores remain computed on the "
            "589-study Ring-1 baseline (+0.00 pp under IIP-0033) unless explicitly updated through an accepted Islet Improvement Proposal."
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
    ceg_payload = json.loads(ceg_path.read_text(encoding="utf-8"))
    counts = verify_ceg_payload(ceg_payload, repo_root=root)
    expected_atlas = (json.dumps(build_layer1_atlas_from_ceg(ceg_payload), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if args.verify:
        if not atlas_path.is_file() or atlas_path.read_bytes() != expected_atlas:
            raise ValueError(f"{ATLAS_REL_PATH} is missing or out of sync with {CEG_REL_PATH}")
        print(f"Verified CEG ({counts['studies']} studies, {counts['outcomes']} outcomes, {counts['edges']} links) and Layer-1 atlas.")
    else:
        atlas_path.write_bytes(expected_atlas)
        print(f"Verified CEG and wrote {ATLAS_REL_PATH} ({len(expected_atlas) / 1024 / 1024:.2f} MB).")


if __name__ == "__main__":
    main()
