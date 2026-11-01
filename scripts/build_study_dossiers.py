#!/usr/bin/env python3
"""
Step 3 & Step 5A: Deterministic Full-Text JATS XML "Study Data Dossier" Builder
For any list of PMIDs (--pilot-5, --ring1, or --pmids), pulls open-access Full-Text JATS XML
from Europe PMC concurrently, strips narrative noise (Introduction, Discussion, References), and extracts:
  1. Parsed XML <table-wrap> data tables (converted to Markdown + footers)
  2. All <fig> captions and panel legends (where sample sizes n, mean±SD, p-values, and modifiers Z live)
  3. Methods paragraphs (cell lines, differentiation recipes, genetic edits, transplant niche)
  4. Quantitative Results paragraphs (containing numbers, comparisons, and effect modifiers)
  5. Linked NCBI GEO (GSE) functional genomics metadata
Writes compact, verifiable dossiers to data/dossiers/{pmid}.json.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
DOSSIER_DIR = BASE_DIR / "data" / "dossiers"

PILOT_5_PMIDS = [
    "35241836",  # Balboa/Otonkoski 2022 Nat Biotechnol
    "37156915",  # Hu/Schrepfer 2024 Nat Biotechnol
    "42690930",  # Kobaisi/Scharfmann 2026 Cell Reports
    "32719542",  # Cai/Kissler 2020 Nat Metab
    "34480863",  # Aghazadeh/Nostro 2021 Cell Stem Cell
]


def validate_pmid(value: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[1-9][0-9]{0,8}", value) is None:
        raise ValueError("PMID must be a canonical positive identifier with at most nine digits")
    return value


def validate_dossier_identity(pmid: str, record: dict, label: str) -> None:
    if not isinstance(record, dict) or str(record.get("pmid", "")) != pmid:
        raise ValueError(f"{label} PMID does not match requested dossier {pmid}")
    pmcid = record.get("pmcid")
    if pmcid is not None and (not isinstance(pmcid, str) or re.fullmatch(r"PMC[1-9][0-9]*", pmcid) is None):
        raise ValueError(f"{label} PMCID must be a canonical PMC identifier")
    if "has_parsed_fulltext_xml" in record and type(record["has_parsed_fulltext_xml"]) is not bool:
        raise ValueError(f"{label} full-text availability must be boolean")


def clean_text(el: ET.Element | None) -> str:
    if el is None:
        return ""
    return " ".join("".join(el.itertext()).split())


def _normalize_doi(doi) -> str:
    return re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", str(doi).strip(), flags=re.I).casefold()


def validate_jats_source_identity(xml_bytes: bytes, pmid: str, record: dict) -> None:
    """Bind fetched article identifiers to the requested source before caching.

    Identifiers in references or embedded sub-articles are not identities of
    the fetched article. An agreeing PMID or PMCID is required; every supplied
    identifier of either type must agree, as must a DOI present on both sides.
    """
    try:
        article = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise ValueError("Fetched source is not valid JATS XML") from exc
    local = lambda element: element.tag.rsplit("}", 1)[-1]
    if local(article) != "article":
        raise ValueError("Fetched JATS source must be one article")
    fronts = [child for child in article if local(child) == "front"]
    metadata = [child for front in fronts for child in front if local(child) == "article-meta"]
    if len(fronts) != 1 or len(metadata) != 1:
        raise ValueError("Fetched JATS source must have one article identity")
    matched = False
    for element in metadata[0]:
        if local(element) != "article-id":
            continue
        kind = element.get("pub-id-type", "").lower()
        value = clean_text(element)
        if kind == "pmid":
            if value != pmid:
                raise ValueError("Fetched JATS PMID conflicts with requested source")
            matched = True
        elif kind in {"pmc", "pmcid", "pmc-uid"}:
            normalized = value if value.startswith("PMC") else "PMC" + value
            if normalized != record.get("pmcid"):
                raise ValueError("Fetched JATS PMCID conflicts with requested source")
            matched = True
        elif kind == "doi" and record.get("doi"):
            if _normalize_doi(value) != _normalize_doi(record["doi"]):
                raise ValueError("Fetched JATS DOI conflicts with requested source")
    if not matched:
        raise ValueError("Fetched JATS source lacks a matching PMID or PMCID")


def parse_table_wrap(tw: ET.Element) -> dict:
    label = clean_text(tw.find("label"))
    caption = clean_text(tw.find("caption"))
    foot = clean_text(tw.find("table-wrap-foot"))
    rows_md = []
    table_el = tw.find(".//table")
    if table_el is not None:
        for tr in table_el.findall(".//tr"):
            cells = [clean_text(c) for c in list(tr) if c.tag in ("th", "td")]
            if any(cells):
                rows_md.append("| " + " | ".join(cells) + " |")
                if len(rows_md) == 1:
                    rows_md.append("| " + " | ".join(["---"] * len(cells)) + " |")
    return {
        "label": label,
        "caption": caption,
        "markdown_table": "\n".join(rows_md[:60]),
        "footnotes": foot,
    }


def parse_jats_fulltext(xml_bytes: bytes) -> dict | None:
    try:
        root = ET.fromstring(xml_bytes)
    except Exception:
        return None
    tables = [parse_table_wrap(tw) for tw in root.findall(".//table-wrap")]
    figures = []
    for fig in root.findall(".//fig"):
        flabel = clean_text(fig.find("label"))
        fcap = clean_text(fig.find("caption"))
        if flabel or fcap:
            figures.append({"label": flabel, "caption": fcap})

    methods_paras = []
    results_paras = []

    body = root.find(".//body")
    if body is not None:
        for sec in body.findall("sec"):
            sec_title = clean_text(sec.find("title")).lower()
            sec_type = (sec.attrib.get("sec-type", "") or "").lower()
            if any(
                k in sec_title or k in sec_type
                for k in ("intro", "background", "discussion", "acknowledg", "competing", "declaration", "author")
            ):
                continue
            is_method = any(
                k in sec_title or k in sec_type
                for k in ("method", "material", "experimental", "procedure", "protocol")
            )
            is_result = any(k in sec_title or k in sec_type for k in ("result", "finding")) or not is_method

            for p in sec.findall(".//p"):
                ptxt = clean_text(p)
                if len(ptxt) < 40:
                    continue
                if is_method:
                    methods_paras.append(ptxt)
                elif is_result:
                    results_paras.append(ptxt)

    if not tables and not figures and not methods_paras and not results_paras:
        return None

    return {
        "raw_xml_char_length": len(xml_bytes),
        "tables_count": len(tables),
        "figures_count": len(figures),
        "methods_paragraphs_count": len(methods_paras),
        "results_paragraphs_count": len(results_paras),
        "xml_tables": tables,
        "figure_legends": figures,
        "methods_paragraphs": methods_paras[:60],
        "results_paragraphs": results_paras[:80],
    }


def fetch_epmc_fulltext_xml(pmcid: str) -> bytes | None:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    for attempt in range(2):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ProjectBeta1-Atlas/0.2"})
            with urllib.request.urlopen(req, timeout=25) as resp:
                data = resp.read()
                if b"<body" in data or b"<fig" in data or b"<sec" in data:
                    return data
        except Exception:
            pass
    return None


def build_one_dossier(pmid: str, rec: dict, geo_list: list[dict], force: bool = False) -> tuple[str, bool, int]:
    validate_pmid(pmid)
    validate_dossier_identity(pmid, rec, "Metadata")
    if type(force) is not bool:
        raise ValueError("force must be boolean")
    if not isinstance(geo_list, list) or any(not isinstance(row, dict) for row in geo_list):
        raise ValueError("GEO metadata must be a list of records")
    cache_path = DOSSIER_DIR.absolute()
    if any(part.is_symlink() for part in (cache_path, *cache_path.parents)):
        raise ValueError("Dossier directory and its ancestors may not be symbolic links")
    out_path = DOSSIER_DIR / f"{pmid}.json"
    if out_path.is_symlink() or (out_path.exists() and (not out_path.is_file() or out_path.stat().st_nlink != 1)):
        raise ValueError("Existing dossier must be a single-link regular file")
    existing = {}
    if out_path.exists():
        try:
            with open(out_path, "r", encoding="utf-8") as f:
                existing = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Existing dossier {pmid} cannot be read as JSON") from exc
        validate_dossier_identity(pmid, existing, "Cached")
        # A failed refresh may retain cached source text. It must never pair
        # those old bytes with a different incoming article identity.
        for key in ("pmcid", "doi"):
            old, new = existing.get(key), rec.get(key)
            if old and new:
                equal = _normalize_doi(old) == _normalize_doi(new) if key == "doi" else old == new
                if not equal:
                    raise ValueError(f"Incoming {key.upper()} conflicts with cached dossier source")
        if not force:
            return pmid, bool(existing.get("has_parsed_fulltext_xml")), 0

    pmcid = rec.get("pmcid")
    fulltext_block = None
    if pmcid:
        xml_bytes = fetch_epmc_fulltext_xml(pmcid)
        if xml_bytes:
            validate_jats_source_identity(xml_bytes, pmid, rec)
            fulltext_block = parse_jats_fulltext(xml_bytes)

    dossier = {
        "pmid": pmid,
        "pmcid": pmcid,
        "doi": rec.get("doi"),
        "pub_year": rec.get("pub_year"),
        "journal": rec.get("journal"),
        "first_author": rec.get("first_author"),
        "senior_author": rec.get("senior_author"),
        "title": rec.get("title"),
        "abstract": rec.get("abstract"),
        "publication_types": rec.get("publication_types", []),
        "mesh_terms": rec.get("mesh_terms", []),
        "keywords": rec.get("keywords", []),
        "matched_pillars": rec.get("matched_pillars", []),
        "linked_geo_datasets": geo_list,
        "has_parsed_fulltext_xml": fulltext_block is not None,
        "fulltext_data_sections": fulltext_block,
    }
    # --force requests reacquisition, not destruction of scientific source data
    # when a stripped manifest or a transient network failure supplies no text.
    for key, value in existing.items():
        if value not in (None, "", [], {}) and dossier.get(key) in (None, "", [], {}):
            dossier[key] = value
    if fulltext_block is None and existing.get("fulltext_data_sections"):
        dossier["fulltext_data_sections"] = existing["fulltext_data_sections"]
        dossier["has_parsed_fulltext_xml"] = bool(existing.get("has_parsed_fulltext_xml"))
    DOSSIER_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=DOSSIER_DIR, encoding="utf-8", delete=False) as fout:
        tmp_name = fout.name
        json.dump(dossier, fout, indent=2, ensure_ascii=False)
        fout.write("\n")
    os.replace(tmp_name, out_path)
    return pmid, bool(dossier.get("has_parsed_fulltext_xml")), 1


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Study Data Dossiers for Project Islet-1 (OpenIslet.org).")
    parser.add_argument("--pilot-5", action="store_true", help="Build dossiers for the 5 Calibration Pilot studies")
    parser.add_argument("--ring1", action="store_true", help="Build dossiers for all 532 Ring 1 Platinum studies")
    parser.add_argument("--pmids", nargs="*", help="Specific PMIDs to build dossiers for")
    parser.add_argument("--workers", type=int, default=12, help="Concurrent Europe PMC download threads")
    parser.add_argument("--force", action="store_true", help="Rebuild existing dossiers")
    args = parser.parse_args()

    if args.ring1:
        with open(RAW_DIR / "ring1_platinum_queue.json", "r", encoding="utf-8") as f:
            r1 = json.load(f)
        target_pmids = [str(s["pmid"]) for s in r1["studies"]]
    elif args.pilot_5:
        target_pmids = PILOT_5_PMIDS
    else:
        target_pmids = args.pmids or PILOT_5_PMIDS

    if args.workers < 1:
        parser.error("--workers must be positive")
    for pmid in target_pmids:
        validate_pmid(pmid)
    if len(target_pmids) != len(set(target_pmids)):
        raise ValueError("Requested dossier PMIDs must be unique")
    target_set = set(target_pmids)

    # Load metadata index
    meta_by_pmid: dict[str, dict] = {}
    manifest_file = RAW_DIR / "corpus_manifest_all_sources.jsonl"
    with open(manifest_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                rec = json.loads(line)
                if str(rec.get("pmid")) in target_set:
                    if str(rec["pmid"]) in meta_by_pmid:
                        raise ValueError("Duplicate requested PMID in metadata manifest")
                    meta_by_pmid[str(rec["pmid"])] = rec

    # Load GEO datasets
    pmid_to_geo: dict[str, list[dict]] = {}
    geo_file = RAW_DIR / "functional_genomics_geo_datasets.json"
    if geo_file.exists():
        with open(geo_file, "r", encoding="utf-8") as f:
            gdata = json.load(f)
        for ds in gdata.get("datasets", []):
            for pid in ds.get("pubmed_ids", []):
                pmid_to_geo.setdefault(str(pid), []).append(ds)

    print(f"[DOSSIER BUILDER] Building {len(target_pmids)} Study Data Dossiers with {args.workers} workers...")
    fulltext_count = 0
    abstract_only_count = 0
    newly_built = 0

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {
            pool.submit(
                build_one_dossier,
                pmid,
                meta_by_pmid.get(pmid, {"pmid": pmid}),
                pmid_to_geo.get(pmid, []),
                args.force,
            ): pmid
            for pmid in target_pmids
        }
        for idx, fut in enumerate(as_completed(futs), 1):
            pmid, has_ft, is_new = fut.result()
            if has_ft:
                fulltext_count += 1
            else:
                abstract_only_count += 1
            newly_built += is_new
            if idx % 50 == 0 or idx == len(target_pmids):
                print(
                    f"  Progress: {idx}/{len(target_pmids)} dossiers ready "
                    f"(FullText XML parsed: {fulltext_count}, Abstract+Metadata: {abstract_only_count}, Newly built: {newly_built})"
                )

    print(
        f"[COMPLETE] {len(target_pmids)} dossiers in {DOSSIER_DIR} | "
        f"FullText XML parsed: {fulltext_count} | Abstract+GEO/MeSH: {abstract_only_count}"
    )


if __name__ == "__main__":
    main()
