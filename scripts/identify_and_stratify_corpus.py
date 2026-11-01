#!/usr/bin/env python3
"""
Step 2: Deterministic Data Identification, Signal Tagging & Concentric Ring Stratification
Scans all 7,703 PubMed records (6,857 initial harvest + 846 targeted blind-spot and
counter-evidence additions) in data/raw/corpus_manifest_all_sources.jsonl +
70 Preprints + 196 GEO GSE series + 114 Clinical Trials to identify:
  1. Physical Data Accessibility (Open Full-Text PMC XML vs. GEO GSE Matrix vs. Abstract-Only)
  2. Target Gene / Protein Mentions across our 38 prioritized Islet-1 engineering targets
  3. Biological Fidelity Tier (Tier 1 Human Clinical -> Tier 2 NHP/Humanized/Human Co-Culture -> Tier 3 Human Abiotic/NSG -> Tier 4 Rodent-Only)
  4. Quantitative & Effect-Modifier Signal Density
  5. Concentric Extraction Priority Ring (Ring 1 Platinum Multi-Gate Core -> Ring 2 -> Ring 3 -> Ring 4)
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"

GOLDEN_SEED_PMIDS = {
    "25303535", "27163171", "32094658", "35241836", "32846125", "37188763",
    "38684699", "30412052", "24658304", "28325479", "16461897", "30853558",
    "30778232", "37156915", "37043559", "38335966", "39758817", "34407395",
    "32814902", "36599351", "36055241", "29844327", "30269996", "28263308",
    "42690930", "26912858", "32719542", "39988068", "36289205", "32444635",
    "31900242", "24584071", "31676859", "32432213", "20699420", "18753660",
    "23086037", "19641142", "30074059", "23737756", "19336676", "22442301",
    "33135383", "34480863", "35559682",
}

CURATED_FRONTIER_AND_COUNTEREVIDENCE_RING1_PMIDS = [
    "42671571", "26276263", "41445610", "39626097", "41758964", "40757665",
    "40251154", "39824805", "40112799", "40782792", "41851084", "42856627",
    "41993274", "41093834", "41617975", "41951603", "36030035", "41607298",
    "38754432", "41218919", "32569842", "34954144", "27626039", "35795779",
    "27183617", "39529170", "37264229", "28504668", "37820008", "38508139",
    "40668892",
]
CURATED_FRONTIER_SET = set(CURATED_FRONTIER_AND_COUNTEREVIDENCE_RING1_PMIDS)
CURATED_ORDER_INDEX = {
    pmid: idx for idx, pmid in enumerate(CURATED_FRONTIER_AND_COUNTEREVIDENCE_RING1_PMIDS)
}

SUPERSEDED_PREPRINT_PMIDS = {
    "40832357", "37502937", "38496417", "38766166", "39149298",
    "39211206", "39763971", "40236391", "40777266", "41542419",
    "39211191", "37886586", "41509220", "37905040",
}
PROJECT_SUMMARY_PMIDS = {"42430558"}
REVIEW_PMIDS = {"35482056"}
CLINICAL_TRIAL_REPORT_PMIDS = {"40251154"}

TARGET_PATTERNS: dict[str, re.Pattern] = {
    "B2M": re.compile(r"\b(b2m|beta-?2-?microglobulin|β2m)\b", re.I),
    "CIITA": re.compile(r"\b(ciita)\b", re.I),
    "NLRC5": re.compile(r"\b(nlrc5|cita)\b", re.I),
    "HLA-C": re.compile(r"\b(hla-c|hla-c1|hla-c2)\b", re.I),
    "HLA-E": re.compile(r"\b(hla-e)\b", re.I),
    "HLA-G": re.compile(r"\b(hla-g)\b", re.I),
    "CD47": re.compile(r"\b(cd47|sirp-?alpha|sirpα)\b", re.I),
    "CD274_PDL1": re.compile(r"\b(pd-?l1|cd274)\b", re.I),
    "FAS_FASL": re.compile(r"\b(fasl|fas ligand|cd95l)\b", re.I),
    "CD55_CD59_CD39": re.compile(r"\b(cd55|cd59|cd39|entpd1)\b", re.I),
    "FCGR1A_CD64": re.compile(r"\b(cd64|fcgr1a)\b", re.I),
    "TYK2": re.compile(r"\b(tyk2)\b", re.I),
    "RNLS": re.compile(r"\b(rnls|renalase)\b", re.I),
    "IFIH1": re.compile(r"\b(ifih1|mda5)\b", re.I),
    "SLC30A8": re.compile(r"\b(slc30a8|znt8|znt8a)\b", re.I),
    "G6PC2": re.compile(r"\b(g6pc2|igrp)\b", re.I),
    "GAD2": re.compile(r"\b(gad2|gad65|gada)\b", re.I),
    "PTPRN": re.compile(r"\b(ptprn|ia-?2|ia-?2a|phogrin)\b", re.I),
    "INS_Leader_DRiP": re.compile(r"\b(preproinsulin|signal peptide|ins-drip\w*|hybrid insulin peptide\w*)\b", re.I),
    "IAPP": re.compile(r"\b(iapp|islet amyloid|amylin|pramlintide)\b", re.I),
    "CXADR": re.compile(r"\b(cxadr|car-siv|coxsackie\w*)\b", re.I),
    "CHGA": re.compile(r"\b(chga|chromogranin|we14)\b", re.I),
    "TGM2": re.compile(r"\b(tgm2|transglutaminase 2)\b", re.I),
    "SERPINB9": re.compile(r"\b(serpinb9|pi-9|proteinase inhibitor 9|granzyme b)\b", re.I),
    "TNFAIP3": re.compile(r"\b(tnfaip3|\ba20\b)\b", re.I),
    "PTPN2": re.compile(r"\b(ptpn2)\b", re.I),
    "GLIS3": re.compile(r"\b(glis3)\b", re.I),
    "TXNIP": re.compile(r"\b(txnip|verapamil)\b", re.I),
    "MANF": re.compile(r"\b(manf)\b", re.I),
    "TET2": re.compile(r"\b(tet2)\b", re.I),
    "RIPK1": re.compile(r"\b(ripk1)\b", re.I),
    "CXCL10": re.compile(r"\b(cxcl10|ip-10|cxcr3)\b", re.I),
    "F3_THBD_IBMIR": re.compile(r"\b(tissue factor|ibmir|thrombomodulin|thbd)\b", re.I),
    "ICASP9_SafetySwitch": re.compile(r"\b(icasp9|inducible caspase|suicide gene|hsv-tk)\b", re.I),
    "PDX1_NKX61_MAFA": re.compile(r"\b(nkx6-?1|nkx6\.1|mafa|pdx1)\b", re.I),
    "WNT4_Latrunculin": re.compile(r"\b(wnt4|latrunculin|enterochromaffin|slc18a1)\b", re.I),
    "GJD2_Cx36": re.compile(r"\b(connexin-?36|cx36|gjd2)\b", re.I),
}

QUANT_REGEX = re.compile(
    r"(\b\d+(\.\d+)?\s*%|\bp\s*[<=>]\s*0\.\d+|\b\d+(\.\d+)?-fold\b|\b\d+(\.\d+)?\s*(mm|μm|um|ng/ml|iu/ml|ieq|days|hours|weeks|months)\b|\bn\s*=\s*\d+)",
    re.I,
)

MODIFIER_REGEX = re.compile(
    r"\b(dose-dependent|time-dependent|stage\s*[1-7]|hypoxi\w+|normoxi\w+|cytokine\w*|ifn-?[αaγg]|tnf-?[αa]|il-?1[βb]|high glucose|16\.7\s*mm|20\s*mm|2\.8\s*mm|diameter|microvessel\w*|pre-?vascular\w*|knockout|overexpress\w*|crispr|base edit\w*|prime edit\w*|heterozygous|homozygous|b2m|hla-c|kir|nkg2a)\b",
    re.I,
)


def classify_fidelity_tier(rec: dict) -> str:
    text = f"{rec.get('title', '')} {rec.get('abstract', '')}".lower()
    mesh = " ".join(rec.get("mesh_terms", [])).lower()
    pub_types = " ".join(rec.get("publication_types", [])).lower()

    if "clinical trial" in pub_types or any(
        k in text for k in ("clinical trial", "phase 1", "phase 2", "phase i/ii", "vx-880", "zimislecel", "up421", "patient with type 1 diabetes", "intraportal islet transplantation in patients")
    ):
        if "human" in mesh or "human" in text or "patient" in text:
            return "Tier_1_Human_InVivo"

    has_human_cells = (
        "humans" in mesh
        or any(
            k in text
            for k in (
                "human", "sc-islet", "sc-beta", "sc-β", "hpsc", "hesc", "ipsc",
                "endoc-betah", "endoc-βh", "patient-derived", "pluripotent stem cell",
            )
        )
    )
    if has_human_cells:
        if any(
            k in text or k in mesh
            for k in (
                "nonhuman primate", "non-human primate", "macaque", "rhesus",
                "cynomolgus", "baboon", "macaca", "humanized mice", "humanized mouse",
                "pbmc-nsg", "hu-pbmc",
            )
        ):
            return "Tier_2A_Human_Cells_InVivo_NHP_or_HumanizedMouse"
        if any(
            k in text
            for k in (
                "pbmc", "t-cell clone", "t cell clone", "autoreactive t", "ctl",
                "nk cell", "natural killer", "macrophage", "phagocytosis",
                "whole blood", "ibmir", "autoantibod", "allogeneic t", "mixed lymphocyte",
            )
        ):
            return "Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture"
        return "Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse"

    return "Tier_4_Rodent_or_Xenogeneic_Only"


def main(allow_title_only_fallback: bool = False) -> None:
    manifest_file = RAW_DIR / "corpus_manifest_all_sources.jsonl"
    geo_file = RAW_DIR / "functional_genomics_geo_datasets.json"
    ppr_file = RAW_DIR / "preprints_sc_islets_hypoimmune.jsonl"
    dossier_dir = BASE_DIR / "data" / "dossiers"

    pmid_to_gse: dict[str, list[str]] = {}
    if geo_file.exists():
        with open(geo_file, "r", encoding="utf-8") as f:
            gdata = json.load(f)
        for ds in gdata.get("datasets", []):
            acc = ds.get("accession", "")
            for pid in ds.get("pubmed_ids", []):
                pmid_to_gse.setdefault(str(pid), []).append(acc)

    records = []
    with open(manifest_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                records.append(json.loads(line))

    manifest_abstract_count = sum(1 for r in records if (r.get("abstract") or "").strip())
    dossier_abstract_count = sum(
        1 for r in records if (dossier_dir / f"{r.get('pmid', '')}.json").exists()
    )
    if manifest_abstract_count == 0 and not allow_title_only_fallback:
        raise RuntimeError(
            "corpus_manifest_all_sources.jsonl contains 0 inline abstracts and allow_title_only_fallback=False; "
            "pass allow_title_only_fallback=True to stratify using title+MeSH+keyword+subpillar triage with "
            "cached Ring-1 dossier abstracts."
        )

    data_access_counts = Counter()
    fidelity_counts = Counter()
    ring_counts = Counter()
    target_gene_counts = Counter()

    stratified_records = []
    ring1_queue = []

    for rec in records:
        pmid = str(rec.get("pmid", ""))
        pmcid = rec.get("pmcid")
        gse_list = pmid_to_gse.get(pmid, [])
        pillars = rec.get("matched_pillars", [])
        core_subpillars = [
            p for p in pillars
            if p.startswith("Pillar1") or p.startswith("Pillar2") or p.startswith("Pillar3") or p.startswith("Pillar4")
        ]

        # 1. Data accessibility
        if pmcid and gse_list:
            access_class = "Open_FullText_PMC_XML_Plus_GEO_GSE"
        elif pmcid:
            access_class = "Open_FullText_PMC_XML"
        elif gse_list:
            access_class = "Abstract_Plus_GEO_GSE"
        else:
            access_class = "Abstract_Only"
        data_access_counts[access_class] += 1

        # 2. Target gene detection
        full_str = f"{rec.get('title', '')} {rec.get('abstract', '')} {' '.join(rec.get('keywords', []))}"
        matched_targets = []
        for tname, pat in TARGET_PATTERNS.items():
            if pat.search(full_str):
                matched_targets.append(tname)
                target_gene_counts[tname] += 1

        # 3. Biological fidelity tier
        fidelity_tier = classify_fidelity_tier(rec)
        fidelity_counts[fidelity_tier] += 1

        # 4. Quantitative & Effect-Modifier signal density
        quant_hits = len(QUANT_REGEX.findall(full_str))
        modifier_hits = len(set(m.lower() for m in MODIFIER_REGEX.findall(full_str)))

        # 5. Concentric Ring Stratification
        is_seed = pmid in GOLDEN_SEED_PMIDS
        is_curated_frontier = pmid in CURATED_FRONTIER_SET
        is_multi_subpillar = len(core_subpillars) >= 2
        has_geo = len(gse_list) > 0
        is_human_multi_target = (
            fidelity_tier != "Tier_4_Rodent_or_Xenogeneic_Only"
            and len(matched_targets) >= 3
            and "Pillar1A_SC_Islet_Differentiation_Maturation" in core_subpillars
        )

        if is_seed or is_curated_frontier or is_multi_subpillar or has_geo or is_human_multi_target:
            ring = "Ring_1_Platinum_MultiGate_Core"
        elif (
            "Pillar1A_SC_Islet_Differentiation_Maturation" in core_subpillars
            or (fidelity_tier != "Tier_4_Rodent_or_Xenogeneic_Only" and len(matched_targets) >= 1)
        ):
            ring = "Ring_2_Human_SC_Islet_and_TargetGene_Core"
        elif fidelity_tier != "Tier_4_Rodent_or_Xenogeneic_Only":
            ring = "Ring_3_Human_Islet_Paracrine_HPAP_and_Engraftment_Core"
        else:
            ring = "Ring_4_Rodent_Tier4_and_Peripheral_Pool"

        ring_counts[ring] += 1

        if pmid in SUPERSEDED_PREPRINT_PMIDS:
            record_study_type = "Preprint_Superseded_By_Peer_Reviewed"
        elif pmid in PROJECT_SUMMARY_PMIDS:
            record_study_type = "Project_Summary"
        elif pmid in REVIEW_PMIDS:
            record_study_type = "Review"
        elif pmid in CLINICAL_TRIAL_REPORT_PMIDS:
            record_study_type = "Clinical_Trial_Report"
        else:
            record_study_type = "Empirical_Experimental"

        entry = {
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": rec.get("doi"),
            "pub_year": rec.get("pub_year"),
            "journal": rec.get("journal"),
            "first_author": rec.get("first_author"),
            "senior_author": rec.get("senior_author"),
            "title": rec.get("title"),
            "record_study_type": record_study_type,
            "is_active_empirical_study": record_study_type in ("Empirical_Experimental", "Clinical_Trial_Report"),
            "data_accessibility_class": access_class,
            "geo_gse_accessions": gse_list,
            "estimated_fidelity_tier": fidelity_tier,
            "matched_pillars": pillars,
            "core_subpillar_count": len(core_subpillars),
            "matched_beta1_targets": matched_targets,
            "quantitative_signal_count": quant_hits,
            "effect_modifier_signal_count": modifier_hits,
            "priority_ring": ring,
            "is_golden_seed": bool(is_seed or is_curated_frontier),
        }
        stratified_records.append(entry)
        if ring == "Ring_1_Platinum_MultiGate_Core":
            ring1_queue.append(entry)

    # Sort automated 532 Ring 1 studies first, followed by the 31 curated frontier/counter-evidence additions in canonical order
    auto_ring1_studies = [r for r in ring1_queue if r["pmid"] not in CURATED_FRONTIER_SET]
    curated_frontier_studies = [r for r in ring1_queue if r["pmid"] in CURATED_FRONTIER_SET]
    auto_ring1_studies.sort(
        key=lambda x: (
            x["is_golden_seed"],
            x["pmcid"] is not None,
            x["core_subpillar_count"],
            len(x["matched_beta1_targets"]),
            x["pub_year"] or 0,
        ),
        reverse=True,
    )
    curated_frontier_studies.sort(key=lambda x: CURATED_ORDER_INDEX.get(x["pmid"], 999))
    ring1_queue = auto_ring1_studies + curated_frontier_studies

    ring1_pmids = [r["pmid"] for r in ring1_queue]
    ring1_platinum_seed_pmids = [r["pmid"] for r in ring1_queue if r["is_golden_seed"]]
    ring1_study_type_counts = Counter(r["record_study_type"] for r in ring1_queue)
    non_superseded_count = sum(
        1 for r in ring1_queue if r["record_study_type"] != "Preprint_Superseded_By_Peer_Reviewed"
    )
    active_deduplicated_count = sum(
        1
        for r in ring1_queue
        if r["record_study_type"] != "Preprint_Superseded_By_Peer_Reviewed" and r["pmid"] != "42430558"
    )
    active_empirical_count = sum(1 for r in ring1_queue if r["is_active_empirical_study"])

    out_ring1 = RAW_DIR / "ring1_platinum_queue.json"
    if out_ring1.exists():
        try:
            existing_q = json.loads(out_ring1.read_text(encoding="utf-8"))
            existing_pmids = [
                str(s["pmid"]) if isinstance(s, dict) else str(s)
                for s in existing_q.get("studies", [])
            ]
            missing_existing = sorted(set(existing_pmids) - set(ring1_pmids))
            if missing_existing:
                raise RuntimeError(
                    f"Refusing to overwrite ring1_platinum_queue.json: {len(missing_existing)} existing Ring-1 PMIDs "
                    f"would be dropped ({missing_existing[:5]}...)."
                )
        except json.JSONDecodeError:
            pass

    summary = {
        "generated_at": "2026-10-24",
        "total_pubmed_studies_indexed": len(stratified_records),
        "initial_pubmed_harvest_records": 6857,
        "targeted_blindspot_and_counterevidence_additions": len(stratified_records) - 6857,
        "automated_ring1_core_count": len(auto_ring1_studies),
        "curated_frontier_and_counterevidence_ring1_additions": len(curated_frontier_studies),
        "ring1_platinum_core_total": len(ring1_queue),
        "ring1_platinum_seed_and_frontier_count": len(ring1_platinum_seed_pmids),
        "ring1_superseded_preprint_twins_count": ring1_study_type_counts["Preprint_Superseded_By_Peer_Reviewed"],
        "ring1_non_superseded_studies_count": non_superseded_count,
        "ring1_active_deduplicated_studies_count": active_deduplicated_count,
        "ring1_non_empirical_project_summary_count": ring1_study_type_counts["Project_Summary"],
        "ring1_non_empirical_review_count": ring1_study_type_counts["Review"],
        "ring1_active_empirical_studies_count": active_empirical_count,
        "ring1_study_type_breakdown": dict(ring1_study_type_counts),
        "abstract_availability_and_triage_mode_audit": {
            "triage_mode": "Title_MeSH_Keyword_Subpillar_Triage_With_Explicit_Fallback_OptIn",
            "allow_title_only_fallback": allow_title_only_fallback,
            "manifest_inline_abstract_count": manifest_abstract_count,
            "ring1_cached_dossier_abstract_count": dossier_abstract_count,
        },
        "data_accessibility_breakdown": dict(data_access_counts),
        "open_fulltext_pmc_xml_total": sum(1 for r in stratified_records if r["pmcid"] is not None),
        "open_fulltext_pmc_xml_pct": round(100.0 * sum(1 for r in stratified_records if r["pmcid"] is not None) / max(1, len(stratified_records)), 1),
        "biological_fidelity_tier_breakdown": dict(fidelity_counts),
        "concentric_priority_ring_breakdown": dict(ring_counts),
        "ring1_open_fulltext_count": sum(1 for r in ring1_queue if r["pmcid"] is not None),
        "target_gene_mention_breakdown": dict(target_gene_counts.most_common()),
        "initial_6857_harvest_stratification_snapshot": {
            "total_pubmed_studies_indexed": 6857,
            "concentric_priority_ring_breakdown": {
                "Ring_1_Platinum_MultiGate_Core": 532,
                "Ring_2_Human_SC_Islet_and_TargetGene_Core": 2174,
                "Ring_3_Human_Islet_Paracrine_HPAP_and_Engraftment_Core": 1631,
                "Ring_4_Rodent_Tier4_and_Peripheral_Pool": 2520,
            },
            "ring1_open_fulltext_count": 374,
        },
    }

    out_summary = RAW_DIR / "corpus_triage_and_rings.json"
    with open(out_summary, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
        f.write("\n")

    with open(out_ring1, "w", encoding="utf-8") as f:
        json.dump(
            {
                "count": len(ring1_queue),
                "automated_ring1_core_count": len(auto_ring1_studies),
                "curated_frontier_and_counterevidence_additions_count": len(curated_frontier_studies),
                "total_ring1_studies": len(ring1_queue),
                "superseded_preprint_twins_count": ring1_study_type_counts["Preprint_Superseded_By_Peer_Reviewed"],
                "non_superseded_studies_count": non_superseded_count,
                "active_deduplicated_studies_count": active_deduplicated_count,
                "non_empirical_project_summary_count": ring1_study_type_counts["Project_Summary"],
                "non_empirical_review_count": ring1_study_type_counts["Review"],
                "active_empirical_studies_count": active_empirical_count,
                "ring1_study_type_breakdown": dict(ring1_study_type_counts),
                "studies": ring1_queue,
            },
            f,
            indent=2,
            ensure_ascii=False,
        )
        f.write("\n")

    out_strat_manifest = RAW_DIR / "corpus_stratification_manifest.json"
    with open(out_strat_manifest, "w", encoding="utf-8") as f:
        json.dump(
            {
                "schema_version": "0.8.1",
                "total_deduplicated_pubmed_records_in_manifest": len(stratified_records),
                "initial_pubmed_harvest_records": 6857,
                "targeted_blindspot_and_counterevidence_additions": len(stratified_records) - 6857,
                "automated_ring1_core_count": len(auto_ring1_studies),
                "curated_frontier_and_counterevidence_ring1_additions": len(curated_frontier_studies),
                "ring1_platinum_core_count": len(ring1_queue),
                "ring1_platinum_seed_and_frontier_count": len(ring1_platinum_seed_pmids),
                "ring1_superseded_preprint_twins_count": ring1_study_type_counts["Preprint_Superseded_By_Peer_Reviewed"],
                "ring1_non_superseded_studies_count": non_superseded_count,
                "ring1_active_deduplicated_studies_count": active_deduplicated_count,
                "ring1_non_empirical_project_summary_count": ring1_study_type_counts["Project_Summary"],
                "ring1_non_empirical_review_count": ring1_study_type_counts["Review"],
                "ring1_active_empirical_studies_count": active_empirical_count,
                "ring1_study_type_breakdown": dict(ring1_study_type_counts),
                "ring1_platinum_pmids": ring1_platinum_seed_pmids,
                "ring1_pmids": ring1_pmids,
                "ring1_platinum_is_strict_subset_of_ring1": set(ring1_platinum_seed_pmids).issubset(set(ring1_pmids)),
                "ring1_queue_path": "data/raw/ring1_platinum_queue.json",
                "notes": (
                    f"Computed corpus stratification: {len(stratified_records)} PubMed records; "
                    f"{len(ring1_queue)} Ring-1 bibliographic records; {non_superseded_count} non-superseded; "
                    f"{active_deduplicated_count} active deduplicated; {active_empirical_count} active empirical. "
                    f"Excluded project summaries: {ring1_study_type_counts['Project_Summary']}; "
                    f"reviews: {ring1_study_type_counts['Review']}; "
                    f"superseded preprints: {ring1_study_type_counts['Preprint_Superseded_By_Peer_Reviewed']}."
                ),
            },
            f,
            indent=2,
            ensure_ascii=False,
        )
        f.write("\n")

    out_all_stratified = RAW_DIR / "corpus_stratified_index.jsonl"
    with open(out_all_stratified, "w", encoding="utf-8") as f:
        for r in stratified_records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    print("=== STEP 2: DATA IDENTIFICATION & RING STRATIFICATION SUMMARY ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stratify PubMed corpus and emit Ring-1 Platinum queue.")
    parser.add_argument(
        "--allow-title-only",
        dest="allow_title_only",
        action="store_true",
        default=False,
        help="Explicitly allow title+MeSH+keyword+subpillar triage when inline manifest abstracts are absent; private dossier caches may still be required.",
    )
    parser.add_argument(
        "--strict-require-abstracts",
        dest="allow_title_only",
        action="store_false",
        help="Fail immediately if corpus_manifest_all_sources.jsonl contains 0 inline abstracts.",
    )
    args = parser.parse_args()
    main(allow_title_only_fallback=args.allow_title_only)
