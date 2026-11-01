#!/usr/bin/env python3
"""
Step 4: Full Recursive Schema Validator & Cross-Study Causal / Effect-Modifier Graph Builder
1. Recursively validates all extracted JSON files in data/extracted/*.json against schemas/extraction_schema_v0.2.json
   (checking required fields, types, and enum values across all nested objects and arrays).
2. Builds the unified Cross-Study Causal & Effect-Modifier Graph (data/extracted/cross_study_discovery_graph.json):
   - Aggregates all atomic_effects (X -> Y | Z) across studies
   - Aggregates all effect_modifier_interactions (dEffect/dZ != 0) by modifier class Z1..Z8
   - Extracts all canonical biological concept nodes per edge/finding to link directed causal_mechanism_edges
     and cross_study_discovery_bridges across disconnected papers
   - Synthesizes explicit multi-study discovery chains (Study A -> Shared Biological Node -> Study B)
"""

from __future__ import annotations

import json
import re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parent.parent
SCHEMA_PATH = BASE_DIR / "schemas" / "extraction_schema_v0.2.json"
EXTRACTED_DIR = BASE_DIR / "data" / "extracted"

CANONICAL_CONCEPT_PATTERNS: dict[str, str] = {
    "HLA_DR3_DR4_T1DE1_T1DE2_Endotype_Architecture": r"\b(HLA-DR3|HLA-DR4|DR3-DQ2|DR4-DQ8|DR3-T1D|DR4-T1D|T1DE1|T1DE2|IAA-first|GADA-first|endotype)\b",
    "MastCell_Histamine_CPA3_TPSAB1_Axis": r"\b(mast cell|histamine|tryptase|CPA3|TPSAB1|FcεRI|cromolyn|ketotifen|cetirizine)\b",
    "PMN_Neutrophil_NETosis_PAD4_CD99_CD200_Axis": r"\b(polymorphonuclear|\bPMN\b|\bPMNs\b|NETosis|neutrophil extracellular trap|PAD4|PADI4|CD99|CD200|myeloperoxidase)\b",
    "RIPK1_RIPK3_MLKL_Necroptosis_and_Ferroptosis_Axis": r"\b(RIPK1|RIPK3|MLKL|pMLKL|necroptosis|necrostatin|Nec-1|ferroptosis|GPX4|CXCL16|OxLDL)\b",
    "B2M_HLA_Class_I_Complex": r"\b(B2M|beta-2-microglobulin|HLA-A|HLA-B|HLA-C|HLA-ABC|MHC class I\b|MHC-I\b|HLA class I)\b",
    "CIITA_HLA_Class_II_Complex": r"\b(CIITA|HLA-DR|HLA-DQ|HLA-DP|MHC class II\b|MHC-II\b|HLA class II)\b",
    "CD47_SIRPa_DonotEatMe_Axis": r"\b(CD47|SIRP\s*alpha|SIRPα)\b",
    "HLA_E_NKG2A_Inhibitory_Axis": r"\b(HLA-E|HLA-G|Qa-1b|NKG2A|NKG2C|VMAPRTLVL|B2M-HLA-E)\b",
    "PVR_TIGIT_MICA_MICB_NKG2D_Axis": r"\b(PVR|CD155|TIGIT|MICA|MICB|ULBP1|ULBP2|NKG2D|NCR3LG1)\b",
    "PD_L1_CD274_Checkpoint": r"\b(PD-L1|PDL1|CD274|PD-1)\b",
    "HVEM_CD40_BTN3A1_Costimulation": r"\b(HVEM|TNFRSF14|CD40|CD80|CD86|BTN3A1|CTLA4|CTLA-4)\b",
    "FasL_FASLG_Treg_Tolerance_Signal": r"\b(FasL|SA-FasL|FASLG|CD95L)\b",
    "IL2_TGFb_IL10_FOXP3_Treg_Induction": r"\b(FOXP3|FoxP3|Treg|regulatory T cell|IL-10|TGF-beta|TGFB1|IL-2 mutein|2B10)\b",
    "CXCL10_CXCR3_Chemokine_Recruitment": r"\b(CXCL10|IP-10|CXCL9|CXCL11|CCL5|CXCR3)\b",
    "IFNg_IFNa_IL1b_TNFa_Cytokine_Challenge": r"\b(IFN-gamma|IFNγ|IFN-alpha|IFNα|IL-1beta|IL-1β|TNF-alpha|TNFα|pro-inflammatory cytokine)\b",
    "TYK2_JAK1_STAT1_IRF1_Signaling": r"\b(TYK2|TYK2i|BMS-986165|JAK1|JAK2|STAT1|STAT2|IRF1|IRF-1|ruxolitinib|SOCS1|PSMB8|TAP1)\b",
    "TNFAIP3_A20_NFkB_JNK_Brake": r"\b(TNFAIP3|\bA20\b|NF-kB|NF-κB)\b",
    "PTPN2_Phosphatase_STAT1_Brake": r"\b(PTPN2|TC-PTP)\b",
    "RNLS_Glycolytic_NADPH_GSH_Redox_Shield": r"\b(RNLS|Rnls|renalase|pargyline|NADPH|GSH:GSSG|glutathione)\b",
    "ER_Stress_UPR_IRE1_PERK_CHOP_TXNIP": r"\b(ER stress|unfolded protein response|UPR|IRE1|XBP1|sXBP1|PERK|ATF4|ATF6|CHOP|DDIT3|TXNIP|BiP|HSPA5|MANF|thapsigargin|tunicamycin|DTT|Akita)\b",
    "SERPINB9_GranzymeB_Perforin_Apoptosis_Axis": r"\b(SERPINB9|PI-9|granzyme|perforin|caspase-3|cleaved caspase)\b",
    "SLC30A8_ZnT8_Granule_Zinc_Autoantigen": r"\b(SLC30A8|ZnT8|ZnT8A)\b",
    "G6PC2_IGRP_Glycolytic_Autoantigen": r"\b(G6PC2|IGRP)\b",
    "GAD2_GAD65_GABA_Paracrine_Autoantigen": r"\b(GAD2|GAD65|GADA|\bGABA\b)\b",
    "INS_Preproinsulin_Signal_Peptide_DRiP": r"\b(preproinsulin|insulin signal peptide|PPI15-24|INS-DRiP|insulin B9-23)\b",
    "CHGA_WE14_Hybrid_Insulin_Peptide_HIP": r"\b(CHGA|chromogranin A|WE14|hybrid insulin peptide|\bHIP\b|TGM2)\b",
    "IAPP_Amyloid_Oligomer_Inflammasome": r"\b(IAPP|hIAPP|islet amyloid|amylin|pramlintide|NLRP3)\b",
    "CXADR_IFIH1_Enterovirus_CVB_Infection": r"\b(CXADR|CAR-SIV|Coxsackie|CVB3|CVB4|enterovirus|IFIH1|MDA5)\b",
    "CD55_CD59_THBD_Complement_IBMIR_Axis": r"\b(CD55|CD59|CD46|CD39|ENTPD1|THBD|thrombomodulin|IBMIR|complement|tissue factor)\b",
    "Stage7_3D_Suspension_Cytoarchitecture": r"\b(Stage 7|3D suspension|3D aggregat|reaggregat|cluster diameter|microwell|AggreWell|organoid)\b",
    "MAFA_NKX6_1_PDX1_BetaCell_Identity": r"\b(MAFA|NKX6\.1|NKX6-1|PDX1|NEUROG3|NKX2\.2|ISL1|UCN3)\b",
    "WNT4_ESRRG_Mitochondrial_OXPHOS_Switch": r"\b(WNT4|ESRRG|ERR\s*gamma|ERRγ|OXPHOS|mitochondrial respiration|OCR|ECAR|LDHA|HK2|SLC16A1|MCT1|pyruvate)\b",
    "Biphasic_GSIS_and_Calcium_Coupling": r"\b(GSIS|glucose-stimulated insulin secretion|biphasic|first-phase|perifusion|C-peptide|calcium influx|KATP)\b",
    "Alpha_GCG_Delta_SST_Paracrine_CounterRegulation": r"\b(glucagon|\bGCG\b|somatostatin|\bSST\b|alpha cell|delta cell|paracrine)\b",
    "OffTarget_SLC18A1_Enterochromaffin_Depletion": r"\b(SLC18A1|enterochromaffin|SC-EC|TPH1|LMX1A)\b",
    "ZM447439_Aurora_G1_CellCycle_Exit": r"\b(ZM447439|aurora kinase|Ki-67|G0/G1|cell cycle exit|teratoma|iCasp9|HSV-TK)\b",
    "CD31_CD146_Microvessel_Inosculation": r"\b(microvessel|microvascular fragment|CD31|CD146|inosculation|anastomos|revasculariz|endothelial|HUVEC|iEC|VEGFA)\b",
    "Hypoxia_72h_Ischemia_and_O2_Scaffolds": r"\b(hypoxia|ischemi|oxygen tension|pO2|CaO2|calcium peroxide|perfluorocarbon|subcutaneous|omental|intramuscular)\b",
    "Alginate_or_Macrodevice_Encapsulation_FBR": r"\b(alginate|macroencapsul|microencapsul|TheraCyte|PEC-Direct|PEC-Encap|hydrogel|fibrotic|foreign body)\b",
    "Alloreactive_CD8_CD4_TCell_or_DSA_Attack": r"\b(alloreactive|allogeneic T|allograft rejection|PBMC|donor-specific antibod|DSA|CD8\+ T|CD4\+ T|cytotoxic T)\b",
    "Primary_NK_or_Macrophage_Innate_Cytotoxicity": r"\b(NK cell|natural killer|missing-self|missing self|macrophage|phagocytos)\b",
    "BetaCell_Survival_and_Glycemic_Cure": r"\b(normoglycemi|diabetes reversal|graft survival|beta-cell survival|cell viability|apoptosis|cell death|specific lysis)\b",
}


def extract_canonical_concepts(text: str) -> list[str]:
    found = []
    for concept, pat in CANONICAL_CONCEPT_PATTERNS.items():
        if re.search(pat, text, flags=re.IGNORECASE):
            found.append(concept)
    return found


def validate_against_schema(instance: Any, schema_node: dict[str, Any], path: str = "$") -> None:
    """Apply the complete declared JSON Schema and reject non-JSON numbers.

    Python's JSON decoder accepts NaN and Infinity by default. They are not
    scientific observations and must fail before any graph aggregation.
    """
    import math
    from jsonschema import Draft202012Validator
    def finite_json(value, location):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f"{location}: non-finite number is not valid JSON evidence")
        if isinstance(value, dict):
            for key, child in value.items():
                finite_json(child, f"{location}.{key}")
        elif isinstance(value, list):
            for index, child in enumerate(value):
                finite_json(child, f"{location}[{index}]")
    finite_json(instance, path)
    errors = list(Draft202012Validator(schema_node).iter_errors(instance))
    if errors:
        error = errors[0]
        suffix = "".join(f"[{part}]" if isinstance(part, int) else f".{part}" for part in error.absolute_path)
        raise ValueError(f"{path}{suffix}: {error.message}")


def _has_matched_cell_source_and_stage(left: dict, right: dict) -> bool:
    """Unknown labels are missing information, not evidence of compatibility."""
    for field in ("cell_source_class", "canonical_differentiation_stage_evaluated"):
        value = left.get(field)
        if not isinstance(value, str) or not value.strip() or value != right.get(field):
            return False
        if any(token in value.lower() for token in ("unknown", "unspecified", "mixed", "not_reported")):
            return False
    return True


def _validate_graph_references(doc: dict, path: str) -> None:
    for collection, identifier in (("arms", "arm_id"), ("atomic_effects", "effect_id")):
        identifiers = [row[identifier] for row in doc[collection]]
        if any(not isinstance(value, str) or not value for value in identifiers) or len(identifiers) != len(set(identifiers)):
            raise ValueError(f"{path}: {identifier} values must be nonempty and unique")
    arms = {arm["arm_id"]: arm for arm in doc["arms"]}
    for effect in doc["atomic_effects"]:
        focal = effect["focal_arm_id"]
        comparator = effect["comparator_arm_id"]
        if focal not in arms:
            raise ValueError(f"{path}: effect {effect['effect_id']} references an unknown focal arm")
        if comparator is not None and (comparator not in arms or comparator == focal):
            raise ValueError(f"{path}: effect {effect['effect_id']} references an invalid comparator arm")
        if effect.get("log_response_ratio_or_smd") is not None:
            if comparator is None or arms[focal].get("is_reference_control", False) or not arms[comparator].get("is_reference_control", False):
                raise ValueError(f"{path}: numeric effect {effect['effect_id']} requires treatment and reference-control arms")
        elif comparator is not None:
            raise ValueError(f"{path}: non-comparative effect {effect['effect_id']} cannot declare a comparator arm")


def _validate_graph_corpus_identity(studies: list[dict]) -> None:
    for collection, identifier in (("arms", "arm_id"), ("atomic_effects", "effect_id")):
        seen: dict[str, str] = {}
        for study in studies:
            pmid = str(study["study_metadata"]["pmid"])
            for row in study[collection]:
                value = row[identifier]
                namespace = re.match(r"^pmid([0-9]+)_", value, flags=re.IGNORECASE)
                if namespace is not None and namespace[1] != pmid:
                    raise ValueError(f"PMID {pmid}: {identifier} {value!r} declares a different study namespace")
                if value in seen:
                    raise ValueError(f"{identifier} {value!r} is reused across studies {seen[value]} and {pmid}")
                seen[value] = pmid


def main() -> None:
    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    with open(SCHEMA_PATH, "r", encoding="utf-8") as f:
        schema = json.load(f)

    studies = []
    for p in sorted(EXTRACTED_DIR.glob("*.json")):
        if not p.stem.isdigit():
            continue
        with open(p, "r", encoding="utf-8") as f:
            doc = json.load(f)
        validate_against_schema(doc, schema, path=p.name)
        if str(doc["study_metadata"]["pmid"]) != p.stem:
            raise ValueError(f"{p.name}: source filename and study PMID differ")
        _validate_graph_references(doc, p.name)
        studies.append(doc)
    _validate_graph_corpus_identity(studies)

    total_arms = 0
    total_atomic_effects = 0
    total_modifier_interactions = 0
    total_causal_edges = 0

    concept_to_studies: dict[str, set[str]] = defaultdict(set)
    concept_to_evidence: dict[str, list[dict[str, str]]] = defaultdict(list)
    gene_to_studies: dict[str, set[str]] = defaultdict(set)
    gene_to_gates: dict[str, set[str]] = defaultdict(set)
    edges_all = []
    modifiers_by_class: dict[str, list[dict]] = defaultdict(list)
    residual_escapes = []
    pleiotropic_bonuses = []
    regime_thresholds = []
    superseded_preprint_records_excluded = []
    project_summary_records_excluded = []
    review_records_excluded = []

    for doc in studies:
        meta = doc["study_metadata"]
        pmid = str(meta.get("pmid"))
        label = f"PMID:{pmid} ({meta.get('first_author', '')} {meta.get('pub_year', '')})"

        total_arms += len(doc.get("arms", []))
        total_atomic_effects += len(doc.get("atomic_effects", []))
        total_modifier_interactions += len(doc.get("effect_modifier_interactions", []))
        total_causal_edges += len(doc.get("causal_mechanism_edges", []))

        src_type = meta.get("source_text_type", "")
        sup_by = meta.get("superseded_by_pmid")
        if src_type == "Preprint_Superseded_By_Peer_Reviewed" or bool(sup_by):
            superseded_preprint_records_excluded.append(
                {
                    "pmid": pmid,
                    "superseded_by_pmid": sup_by,
                    "exclusion_reason": (
                        f"Superseded_By_Peer_Reviewed_PMID:{sup_by}"
                        if sup_by
                        else "Preprint_Superseded_By_Peer_Reviewed_Enum"
                    ),
                }
            )
            continue
        if src_type == "Review_or_Commentary":
            review_records_excluded.append({"pmid": pmid, "exclusion_reason": "Bibliographic review/context, not an independent experimental study"})
            continue
        if src_type == "Project_Summary":
            project_summary_records_excluded.append(pmid)
            continue

        for arm in doc.get("arms", []):
            for edit in arm.get("gate2_3_4_genetic_and_epitope_edits", []):
                gsym = edit.get("gene_symbol_hgnc")
                fgate = edit.get("functional_gate")
                if gsym:
                    gene_to_studies[gsym].add(label)
                    if fgate:
                        gene_to_gates[gsym].add(fgate)

        for eff in doc.get("atomic_effects", []):
            eff_text = f"{eff.get('focal_factor_X', '')} {eff.get('outcome_variable_Y', '')} {eff.get('verbatim_evidence_quote', '')}"
            for c in extract_canonical_concepts(eff_text):
                concept_to_studies[c].add(label)

        for em in doc.get("effect_modifier_interactions", []):
            mclass = em.get("modifier_class", "Unknown")
            modifiers_by_class[mclass].append({
                "study": label,
                "focal_factor_X": em.get("focal_factor_X"),
                "outcome_Y": em.get("outcome_Y"),
                "modifier_variable_Z": em.get("modifier_variable_Z"),
                "modification_type": em.get("modification_type"),
                "level_1": em.get("level_1_condition_and_effect"),
                "level_2": em.get("level_2_condition_and_effect"),
                "mechanistic_explanation": em.get("mechanistic_explanation"),
            })
            for c in extract_canonical_concepts(
                f"{em.get('focal_factor_X', '')} {em.get('modifier_variable_Z', '')} {em.get('mechanistic_explanation', '')}"
            ):
                concept_to_studies[c].add(label)

        arms_list = doc.get("arms", []) or []
        arm_sources = sorted({
            (a.get("gate1_substrate_and_differentiation") or {}).get("cell_source_class")
            for a in arms_list
            if (a.get("gate1_substrate_and_differentiation") or {}).get("cell_source_class")
        })
        arm_stages = sorted({
            (a.get("gate1_substrate_and_differentiation") or {}).get("canonical_differentiation_stage_evaluated")
            for a in arms_list
            if (a.get("gate1_substrate_and_differentiation") or {}).get("canonical_differentiation_stage_evaluated")
        })
        for edge in doc.get("causal_mechanism_edges", []):
            u = edge["upstream_node"]
            v = edge["downstream_node"]
            ev_text = f"{u} {v} {edge.get('evidence_source', '')}".lower().replace("β", "beta").replace("α", "alpha").replace("δ", "delta").replace("γ", "gamma")
            matched_source = arm_sources[0] if len(arm_sources) == 1 else None
            matched_stage = arm_stages[0] if len(arm_stages) == 1 else None
            has_rodent_token = any(
                k in ev_text
                for k in ("c57bl/6", "c57bl6", "mouse", "mice", "murine", "rodent", "rat ", "rats", "nod/scid", "nsg", "balb/c", "nod ")
            )
            has_explicit_rodent_target = any(
                k in ev_text
                for k in (
                    "murine beta", "mouse beta", "murine islet", "mouse islet", "rat islet", "rat beta",
                    "rodent islet", "rodent beta", "min6", "ins-1", "nit-1", "beta-tc",
                )
            )
            has_explicit_human_target = any(
                k in ev_text
                for k in (
                    "sc-islet", "sc-beta", "psc-islet", "ipsc", "hipsc", "hesc", "hues8", "cyt49", "mel1",
                    "human beta", "human islet", "human sc-", "patient-derived", "cadaveric", "t1d sc-", "znt8 lof sc-",
                )
            )
            has_human_effector = any(
                k in ev_text
                for k in ("human nk", "human pbmc", "human t cell", "human cd8", "human cd4", "human macrophage", "humanized")
            )
            has_rodent_effector = any(
                k in ev_text
                for k in ("mouse nk", "murine nk", "mouse t cell", "murine t cell", "mouse splenocyte", "murine splenocyte", "mouse macrophage")
            )
            has_human_reagent = any(
                k in ev_text
                for k in ("human interferon", "human ifn", "human il-", "human tnf", "human cytokine")
            )
            has_rodent_host = any(
                k in ev_text
                for k in (
                    "in mice", "into mice", "in nsg", "nsg mice", "nod/scid", "diabetic mice",
                    "recipient mice", "grafted mice", "transplanted in mice", "transplanted into mice",
                    "prevent diabetes in mice", "mice grafted", "animals transplanted",
                )
            )
            has_generic_human = "human" in ev_text and not (has_explicit_rodent_target and (has_human_effector or has_human_reagent))
            has_human_target = has_explicit_human_target or (has_generic_human and not has_explicit_rodent_target)

            if has_explicit_rodent_target and has_explicit_human_target:
                target_cell_species = "Human_and_Rodent_Mixed_Context"
                matched_source = matched_source or ("Multi_Arm_Mixed_Context:" + "|".join(arm_sources) if arm_sources else "Human_and_Rodent_Mixed_Context")
                matched_stage = matched_stage or "Mixed_Human_and_Rodent_Assay"
            elif (has_explicit_rodent_target or has_rodent_token) and not has_human_target:
                target_cell_species = "Rodent_or_Murine"
                matched_source = "Rodent_or_Murine_Sentence_Specific_Context"
                matched_stage = "InVivo_or_ExVivo_Rodent_Assay"
            else:
                target_cell_species = (
                    "Human"
                    if (has_human_target or any(s in ("Primary_Human_Islet", "iPSC_healthy_donor", "iPSC_T1D_patient", "hESC") for s in arm_sources))
                    else ("NHP" if "NHP_Primary_Islet" in arm_sources else "Unspecified")
                )
                if len(arm_sources) > 1:
                    if any(k in ev_text for k in ("rhesus", "macaque", "cynomolgus", "nhp", "monkey")) and "NHP_Primary_Islet" in arm_sources:
                        matched_source = "NHP_Primary_Islet"
                    elif any(k in ev_text for k in ("sc-islet", "sc-beta", "ipsc", "hipsc")) and "iPSC_healthy_donor" in arm_sources:
                        matched_source = "iPSC_healthy_donor"
                    elif any(k in ev_text for k in ("hesc", "hues8", "cyt49")) and "hESC" in arm_sources:
                        matched_source = "hESC"
                    else:
                        matched_source = "Multi_Arm_Mixed_Context:" + "|".join(arm_sources)

            host_recipient_species = (
                "Rodent_or_Murine"
                if (has_rodent_host or (has_rodent_token and target_cell_species == "Human" and not has_rodent_effector))
                else ("NHP" if any(k in ev_text for k in ("rhesus", "macaque", "cynomolgus", "baboon")) else "InVitro_or_Same_As_Target")
            )
            immune_effector_species = (
                "Human"
                if has_human_effector
                else ("Rodent_or_Murine" if has_rodent_effector else "Unspecified_or_Abiotic")
            )
            reagent_species = "Human" if has_human_reagent else "Unspecified"

            if len(arm_stages) > 1 and matched_stage is None:
                matched_stage = "Multi_Arm_Mixed_Stage:" + "|".join(arm_stages)
            edge_ctx = {
                "cell_source_class": matched_source,
                "target_cell_species": target_cell_species,
                "host_recipient_species": host_recipient_species,
                "immune_effector_species": immune_effector_species,
                "reagent_species": reagent_species,
                "canonical_differentiation_stage_evaluated": matched_stage,
                "study_arm_cell_source_classes": arm_sources,
                "study_arm_differentiation_stages": arm_stages,
                "max_biological_fidelity_tier": meta.get("max_biological_fidelity_tier"),
                "epistemic_class": "Sentence_Heuristic_Candidate_Association",
            }
            u_concepts = extract_canonical_concepts(u)
            v_concepts = extract_canonical_concepts(v)
            for c in sorted(set(u_concepts + v_concepts)):
                concept_to_studies[c].add(label)
                concept_to_evidence[c].append({
                    "study": label,
                    "role": "causal_edge",
                    "statement": f"{u} --[{edge['causal_relation']}]--> {v}",
                })
            edges_all.append({
                "study": label,
                "study_context": edge_ctx,
                "upstream_node": u,
                "upstream_concepts": u_concepts,
                "causal_relation": edge["causal_relation"],
                "downstream_node": v,
                "downstream_concepts": v_concepts,
                "evidence_source": edge["evidence_source"],
            })

        bridges = doc.get("cross_study_discovery_bridges", {})
        for pfind in bridges.get("collateral_and_pleiotropic_findings", []):
            pleiotropic_bonuses.append({"study": label, **pfind})
            for c in extract_canonical_concepts(pfind.get("observation", "")):
                concept_to_studies[c].add(label)
                concept_to_evidence[c].append({
                    "study": label,
                    "role": f"pleiotropy ({pfind.get('is_beneficial_synergy_or_deleterious_tradeoff')})",
                    "statement": pfind.get("observation", ""),
                })
        for rvec in bridges.get("residual_failure_and_escape_vectors", []):
            residual_escapes.append({"study": label, **rvec})
            for c in extract_canonical_concepts(
                f"{rvec.get('residual_loss_or_limitation', '')} {rvec.get('escape_mechanism_identified', '')}"
            ):
                concept_to_studies[c].add(label)
                concept_to_evidence[c].append({
                    "study": label,
                    "role": "residual_bottleneck",
                    "statement": rvec.get("residual_loss_or_limitation", ""),
                })
        for rthresh in bridges.get("regime_thresholds", []):
            regime_thresholds.append({"study": label, **rthresh})

    shared_bridge_nodes = {
        k: {
            "studies_count": len(v),
            "studies": sorted(list(v)),
            "linked_mechanisms_and_bottlenecks": concept_to_evidence[k],
        }
        for k, v in sorted(concept_to_studies.items(), key=lambda x: (-len(x[1]), x[0]))
        if len(v) >= 2
    }

    target_gene_atlas = {
        g: {
            "studies_count": len(s_set),
            "functional_gates": sorted(list(gene_to_gates[g])),
            "studies": sorted(list(s_set)),
        }
        for g, s_set in sorted(gene_to_studies.items(), key=lambda x: (-len(x[1]), x[0]))
    }

    # Build Multi-Hop Cross-Study Causal Chains: Study A (Node_U -> Node_M) + Study B (Node_M -> Node_V)
    # where Study A != Study B and Node_U != Node_M != Node_V
    incoming_by_concept: dict[str, list[dict[str, Any]]] = defaultdict(list)
    outgoing_by_concept: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for edge in edges_all:
        u_nodes = sorted(set(edge["upstream_concepts"]))
        v_nodes = sorted(set(edge["downstream_concepts"]))
        for un in u_nodes:
            for vn in v_nodes:
                if un != vn:
                    outgoing_by_concept[un].append({"mid": vn, **edge})
                    incoming_by_concept[vn].append({"start": un, **edge})

    TERMINAL_OUTCOME_ONLY_NODES = {
        "BetaCell_Survival_and_Glycemic_Cure",
        "Biphasic_GSIS_and_Calcium_Coupling",
        "Alloreactive_CD8_CD4_TCell_or_DSA_Attack",
        "Primary_NK_or_Macrophage_Innate_Cytotoxicity",
    }

    candidates_by_bridge: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen_chains = set()
    for mid_node, in_edges in sorted(incoming_by_concept.items()):
        if mid_node in TERMINAL_OUTCOME_ONLY_NODES:
            continue
        out_edges = outgoing_by_concept.get(mid_node, [])
        for e1 in in_edges[:25]:
            for e2 in out_edges[:25]:
                if e1["study"] == e2["study"]:
                    continue
                v_end = e2["mid"]
                if len({e1["start"], mid_node, v_end}) < 3:
                    continue
                chain_key = (e1["start"], mid_node, v_end, e1["study"], e2["study"])
                if chain_key in seen_chains:
                    continue
                seen_chains.add(chain_key)
                c1 = e1.get("study_context") or {}
                c2 = e2.get("study_context") or {}
                same_ctx = _has_matched_cell_source_and_stage(c1, c2)
                candidates_by_bridge[mid_node].append({
                    "origin_node": e1["start"],
                    "hop1_study": e1["study"],
                    "hop1_relation": e1["causal_relation"],
                    "hop1_context": c1,
                    "bridge_node": mid_node,
                    "hop2_study": e2["study"],
                    "hop2_relation": e2["causal_relation"],
                    "hop2_context": c2,
                    "terminal_node": v_end,
                    "cross_context_compatibility_note": (
                        "Matched cell source and differentiation stage across Hop 1 and Hop 2"
                        if same_ctx
                        else "Cross-context hypothesis path (requires experimental validation across differing cell sources/stages/species)"
                    ),
                    "hop1_evidence": e1["evidence_source"],
                    "hop2_evidence": e2["evidence_source"],
                })

    for mid_node, ch_list in candidates_by_bridge.items():
        ch_list.sort(
            key=lambda ch: (
                0 if str(ch.get("cross_context_compatibility_note", "")).startswith("Matched cell source") else 1,
                -(
                    len(concept_to_studies.get(ch["origin_node"], ()))
                    + len(concept_to_studies.get(ch["bridge_node"], ()))
                    + len(concept_to_studies.get(ch["terminal_node"], ()))
                ),
                ch["origin_node"],
                ch["terminal_node"],
                ch["hop1_study"],
                ch["hop2_study"],
            )
        )

    untruncated_bridge_counts = {k: len(v) for k, v in sorted(candidates_by_bridge.items())}
    chains_total_before_cap = sum(untruncated_bridge_counts.values())

    multi_hop_chains = []
    bridge_keys = sorted(candidates_by_bridge.keys())
    rr_idx = 0
    while len(multi_hop_chains) < 500:
        added_in_round = False
        for bk in bridge_keys:
            if rr_idx < len(candidates_by_bridge[bk]):
                multi_hop_chains.append(candidates_by_bridge[bk][rr_idx])
                added_in_round = True
                if len(multi_hop_chains) >= 500:
                    break
        if not added_in_round:
            break
        rr_idx += 1

    exported_bridge_counts: dict[str, int] = defaultdict(int)
    for ch in multi_hop_chains:
        exported_bridge_counts[ch["bridge_node"]] += 1

    tier_counts: dict[str, int] = defaultdict(int)
    source_type_counts: dict[str, int] = defaultdict(int)
    direction_counts: dict[str, int] = defaultdict(int)
    orientation_sign_counts: dict[str, int] = defaultdict(int)
    penalty_effects_count = 0
    exact_span_verified_count = 0
    quantitative_lrr_count = 0
    active_arms_count = 0
    active_effects_count = 0
    active_modifiers_count = 0
    active_causal_edges_count = 0
    non_superseded_arms_count = 0
    non_superseded_effects_count = 0
    non_superseded_modifiers_count = 0
    non_superseded_causal_edges_count = 0
    for doc in studies:
        m = doc.get("study_metadata", {})
        src_type = m.get("source_text_type", "")
        is_sup = src_type == "Preprint_Superseded_By_Peer_Reviewed" or bool(m.get("superseded_by_pmid"))
        if not is_sup:
            non_superseded_arms_count += len(doc.get("arms", []))
            non_superseded_effects_count += len(doc.get("atomic_effects", []))
            non_superseded_modifiers_count += len(doc.get("effect_modifier_interactions", []))
            non_superseded_causal_edges_count += len(doc.get("causal_mechanism_edges", []))
        if is_sup or src_type == "Project_Summary":
            continue
        active_arms_count += len(doc.get("arms", []))
        active_modifiers_count += len(doc.get("effect_modifier_interactions", []))
        active_causal_edges_count += len(doc.get("causal_mechanism_edges", []))
        tier_counts[str(m.get("max_biological_fidelity_tier", "Unknown"))] += 1
        source_type_counts[str(src_type or "Unknown")] += 1
        for eff in doc.get("atomic_effects", []):
            active_effects_count += 1
            direction_counts[str(eff.get("effect_direction", "Unknown"))] += 1
            if (
                eff.get("is_functional_or_developmental_penalty") is True
                or eff.get("secretome_or_differentiation_penalty_observed") in (True, "True", "Borderline")
            ):
                penalty_effects_count += 1
            sp_val = eff.get("source_span")
            if (isinstance(sp_val, dict) and sp_val.get("norm") == "exact") or ("norm=exact" in str(sp_val or "")):
                exact_span_verified_count += 1
            if eff.get("log_response_ratio_or_smd") is not None:
                quantitative_lrr_count += 1
            osign = (eff.get("endpoint_polarity_registry") or {}).get("orientation_sign")
            if osign is not None:
                orientation_sign_counts[str(osign)] += 1

    matched_ctx_chains = sum(
        1 for ch in multi_hop_chains
        if str(ch.get("cross_context_compatibility_note", "")).startswith("Matched cell source")
    )

    active_studies_list = [
        doc for doc in studies
        if doc.get("study_metadata", {}).get("source_text_type") not in ("Preprint_Superseded_By_Peer_Reviewed", "Project_Summary")
        and not bool(doc.get("study_metadata", {}).get("superseded_by_pmid"))
    ]
    review_pmids_dynamic = sorted(
        str(doc["study_metadata"]["pmid"])
        for doc in active_studies_list
        if doc["study_metadata"].get("source_text_type") == "Review_or_Commentary"
        or "review" in str((doc.get("screening_decision") or {}).get("exclusion_reason") or "").lower()
    )
    active_preprint_pmids_dynamic = sorted(
        str(doc["study_metadata"]["pmid"])
        for doc in active_studies_list
        if any(
            tok in str(doc["study_metadata"].get("journal") or "").lower()
            for tok in ("biorxiv", "medrxiv", "preprint")
        )
    )
    n_active_dedup = len(active_studies_list)
    n_review_dyn = len(review_pmids_dynamic)
    n_preprint_dyn = len(active_preprint_pmids_dynamic)
    n_active_empirical = n_active_dedup - n_review_dyn
    n_strict_peer_reviewed_empirical = n_active_dedup - n_review_dyn - n_preprint_dyn

    graph_out = {
        "schema_version": "0.2.0",
        "validated_studies_count": len(studies),
        "superseded_preprint_records_excluded_from_graph": superseded_preprint_records_excluded,
        "project_summary_records_excluded_from_graph": project_summary_records_excluded,
        "review_records_excluded_from_graph": review_records_excluded,
        "total_arms_extracted": total_arms,
        "total_atomic_effects_extracted": total_atomic_effects,
        "total_effect_modifier_interactions_extracted": total_modifier_interactions,
        "total_causal_mechanism_edges_extracted": total_causal_edges,
        "multi_hop_cross_study_causal_chains_count": len(multi_hop_chains),
        "corpus_epistemic_and_fidelity_diagnostics": {
            "active_deduplicated_studies_count": n_active_dedup,
            "active_empirical_studies_count": n_active_empirical,
            "non_superseded_studies_count": len(studies) - len(superseded_preprint_records_excluded),
            "superseded_preprint_studies_excluded_count": len(superseded_preprint_records_excluded),
            "project_summary_studies_excluded_count": len(project_summary_records_excluded),
            "review_studies_in_corpus_count": n_review_dyn,
            "review_studies_pmids": review_pmids_dynamic,
            "active_preprints_without_in_corpus_peer_reviewed_twin_pmids": active_preprint_pmids_dynamic,
            "ring1_corpus_count_reconciliation": {
                "validated_studies_total": len(studies),
                "superseded_preprints_excluded": len(superseded_preprint_records_excluded),
                "project_summaries_excluded": len(project_summary_records_excluded),
                "non_superseded_studies_total": len(studies) - len(superseded_preprint_records_excluded),
                "active_deduplicated_studies_total": n_active_dedup,
                "review_records_flagged": n_review_dyn,
                "active_preprints_without_in_corpus_twin_flagged": n_preprint_dyn,
                "strict_active_peer_reviewed_empirical_studies_total": n_strict_peer_reviewed_empirical,
            },
            "active_deduplicated_arms_count": active_arms_count,
            "active_deduplicated_effects_count": active_effects_count,
            "active_deduplicated_modifier_interactions_count": active_modifiers_count,
            "active_deduplicated_causal_edges_count": active_causal_edges_count,
            "non_superseded_arms_count": non_superseded_arms_count,
            "non_superseded_effects_count": non_superseded_effects_count,
            "non_superseded_modifier_interactions_count": non_superseded_modifiers_count,
            "non_superseded_causal_edges_count": non_superseded_causal_edges_count,
            "studies_by_biological_fidelity_tier": dict(sorted(tier_counts.items())),
            "studies_by_source_text_type": dict(sorted(source_type_counts.items())),
            "exact_source_span_verified_effects_count": exact_span_verified_count,
            "quantitative_lrr_effects_count": quantitative_lrr_count,
            "heuristic_penalty_effects_count": penalty_effects_count,
            "direction_and_penalty_method_counts": dict(Counter(
                effect.get("direction_adjudication_method", "rule_based_heuristic")
                for study in active_studies_list for effect in study.get("atomic_effects", [])
            )),
            "adjudication_semantics": "Rule-based heuristic extraction; no independent manual accuracy dataset is attached. Classifier agreement is not measured sensitivity or specificity.",
            "adjudicated_penalty_effects_count": penalty_effects_count,
            "deprecated_field_aliases": {
                "adjudicated_penalty_effects_count": {
                    "alias_of": "heuristic_penalty_effects_count",
                    "reason": "Legacy label; these flags were not independently adjudicated accuracy observations."
                }
            },
            "effects_by_direction": dict(sorted(direction_counts.items())),
            "effects_by_polarity_orientation_sign": dict(sorted(orientation_sign_counts.items())),
            "matched_context_multi_hop_chains_count": matched_ctx_chains,
            "cross_context_hypothesis_chains_count": len(multi_hop_chains) - matched_ctx_chains,
        },
        "multi_hop_hypothesis_paths_enumeration_policy": {
            "estimand_role": "Heuristic cross-study hypothesis-generation paths (NOT used as quantitative input weights in stack ranking)",
            "per_bridge_node_incoming_cap": 25,
            "per_bridge_node_outgoing_cap": 25,
            "global_exported_chain_cap": 500,
            "chains_total_before_cap": chains_total_before_cap,
            "active_non_terminal_bridge_concepts_represented_in_export": len(exported_bridge_counts),
            "ordering_basis": "Balanced round-robin selection across all non-terminal bridge concepts (prioritizing matched cell-source/stage context and total supporting study count within each bridge concept)",
            "untruncated_per_bridge_node_chain_counts": untruncated_bridge_counts,
            "exported_per_bridge_node_chain_counts": dict(sorted(exported_bridge_counts.items())),
        },
        "target_gene_atlas": target_gene_atlas,
        "shared_cross_study_bridge_nodes": shared_bridge_nodes,
        "multi_hop_cross_study_causal_chains": multi_hop_chains,
        "effect_modifiers_by_canonical_class": dict(modifiers_by_class),
        "cross_gate_pleiotropic_discoveries": pleiotropic_bonuses,
        "residual_failure_and_complementary_lever_pairs": residual_escapes,
        "regime_thresholds": regime_thresholds,
        "all_causal_edges": edges_all,
    }

    out_file = EXTRACTED_DIR / "cross_study_discovery_graph.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(graph_out, f, indent=2, ensure_ascii=False)
        f.write("\n")

    print(
        f"[GRAPH VALIDATED] {len(studies)} studies | {total_arms} arms | "
        f"{total_atomic_effects} atomic effects | {total_modifier_interactions} effect modifiers | "
        f"{total_causal_edges} causal edges | {len(shared_bridge_nodes)} shared cross-study bridge nodes"
    )


if __name__ == "__main__":
    main()
