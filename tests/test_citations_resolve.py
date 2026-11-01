#!/usr/bin/env python3
"""
OpenIslet (Project Islet-1) Citation-Resolution Gate.

Every `PMID:<n>` written in the prose / code / protocol layer (README, CLAUDE.md, docs, protocol, graveyard and omics-prior JSON, scripts, tests)
must (a) resolve to a real PubMed record and (b) when it is written next to an "Author et al. YYYY" label, the label's surname and
year must match the record's first author and publication year.  The extraction pipeline already verifies the 563-study corpus; this
gate covers the hand-written layer, where digit-level PMID errors (e.g. a thermoelectrics paper cited as a PD-L1 palmitoylation
reference) are otherwise invisible to every other check.

Offline by default: PubMed fields are checked against `data/raw/doc_citation_metadata_cache.json` (a snapshot of NCBI E-utilities
`esummary`).  To refresh it (network required) after editing any citation:

    python3 tests/test_citations_resolve.py --refresh-cache

`--prune-cache` (offline) drops cache entries that no scanned file references any more.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE_PATH = ROOT / "data" / "raw" / "doc_citation_metadata_cache.json"
SCAN_GLOBS = (
    "ISLET-1-BLUEPRINT.md",
    "README.md",
    "CLAUDE.md",
    "docs/*.md",
    "data/protocols/*.json",
    "data/graveyard/*.json",
    "data/omics/target_prior_table.json",
    "scripts/*.py",
    "tests/*.py",
)
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="

PMID_RE = re.compile(r"PMID:? ?(\d{6,8})\b")
# "Lei et al. 2022 PMID:...", "Saunders et al. 2019 Cell Metab PMID:...", "Wu 2024 PMID:...", "Lin et al. 2023 (PMID:...)"
LABEL_RE = re.compile(
    r"([A-Z][A-Za-zÀ-ſ\-' ]{1,40}?)(?: et al\.?| & [A-Z][A-Za-zÀ-ſ\-]+)?,? ?\(?((?:19|20)\d{2})\)?"
    r"[\s,;:\[\(`]{0,6}(?:[A-Za-z0-9 ./\-]{0,24}?)PMID:? ?(\d{6,8})\b"
)
MALFORMED_YEAR_RE = re.compile(
    r"([A-Z][A-Za-zÀ-ſ\-' ]{1,40}?)(?: et al\.?| & [A-Z][A-Za-zÀ-ſ\-]+),? ?\(?((?:19|20)[0-9A-Za-z]{2}|NaN|XXXX|TBD)\)?[\s,;:\[\(`]{0,6}(?:[A-Za-z0-9 ./\-]{0,24}?)PMID:? ?(\d{6,8})\b"
)


def _fold(text: str) -> str:
    """Accent-fold, lower-case and keep letters only (Lizée == Lizee, Müller-Kuller == Muller Kuller)."""
    decomposed = unicodedata.normalize("NFKD", text)
    return re.sub(r"[^a-z]", "", "".join(c for c in decomposed if not unicodedata.combining(c)).lower())


def _surname_of(pubmed_author: str) -> str:
    """'Lissanu Deribe Y' -> 'Lissanu Deribe'; 'Zanta-Boussif MA' -> 'Zanta-Boussif' (PubMed lists 'Surname INITIALS')."""
    tokens = pubmed_author.split()
    while len(tokens) > 1 and tokens[-1].isupper() and len(tokens[-1]) <= 3:
        tokens.pop()
    return " ".join(tokens)


def scan_files() -> list[Path]:
    seen: dict[Path, None] = {}
    for pattern in SCAN_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            if path.is_file():
                seen[path] = None
    return list(seen)


def harvest() -> tuple[dict[str, list[str]], list[tuple[str, str, str, str]]]:
    """Return ({pmid: [file:line, ...]}, [(file:line, label_surname, label_year, pmid), ...])."""
    where: dict[str, list[str]] = defaultdict(list)
    labels: list[tuple[str, str, str, str]] = []
    for path in scan_files():
        rel = path.relative_to(ROOT).as_posix()
        for lineno, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            for m in PMID_RE.finditer(line):
                where[m.group(1)].append(f"{rel}:{lineno}")
            for m in LABEL_RE.finditer(line):
                name = m.group(1).strip().split()[-1]
                labels.append((f"{rel}:{lineno}", name, m.group(2), m.group(3)))
            if path.name != "test_citations_resolve.py":
                for m_bad in MALFORMED_YEAR_RE.finditer(line):
                    if not m_bad.group(2).isdigit():
                        name = m_bad.group(1).strip().split()[-1]
                        labels.append((f"{rel}:{lineno}", name, m_bad.group(2), m_bad.group(3)))
    return where, labels


def load_cache() -> dict:
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    return {"_meta": {}, "pmids": {}}


def fetch_esummary(pmids: list[str]) -> dict[str, dict]:
    """Fetch first author / year / journal / title for each PMID (batches of 80, with retries)."""
    out: dict[str, dict] = {}
    for start in range(0, len(pmids), 80):
        chunk = pmids[start:start + 80]
        result = None
        for attempt in range(4):
            try:
                req = urllib.request.Request(ESUMMARY_URL + ",".join(chunk), headers={"User-Agent": "openislet-citation-gate/1.0"})
                with urllib.request.urlopen(req, timeout=60) as resp:
                    result = json.load(resp).get("result", {})
                break
            except Exception:  # network hiccup: back off and retry
                time.sleep(1.5 * (attempt + 1))
        if result is None:
            raise RuntimeError("NCBI E-utilities unreachable; cache not refreshed")
        for pmid in chunk:
            rec = result.get(pmid)
            if not rec or "error" in rec or not rec.get("title"):
                out[pmid] = {"error": "not_found_in_pubmed"}
                continue
            authors = rec.get("authors") or []
            out[pmid] = {
                "first_author": authors[0]["name"] if authors else "",
                "pubdate_year": str(rec.get("pubdate", ""))[:4],
                "journal": rec.get("source", ""),
                "title": rec.get("title", ""),
            }
        time.sleep(0.4)
    return out


def check(where: dict[str, list[str]], labels: list[tuple[str, str, str, str]], pmids_meta: dict[str, dict]) -> dict[str, list]:
    problems: dict[str, list] = {"not_in_cache": [], "not_in_pubmed": [], "label_mismatch": []}
    for pmid, locs in sorted(where.items()):
        meta = pmids_meta.get(pmid)
        if meta is None:
            problems["not_in_cache"].append((pmid, locs[:3]))
        elif not isinstance(meta, dict) or meta.get("error"):
            problems["not_in_pubmed"].append((pmid, locs[:3]))
        elif (
            not isinstance(meta.get("first_author"), str)
            or not meta["first_author"].strip()
            or not isinstance(meta.get("title"), str)
            or not meta["title"].strip()
            or not isinstance(meta.get("pubdate_year"), str)
            or not meta["pubdate_year"].isdigit()
            or len(meta["pubdate_year"]) != 4
        ):
            problems["not_in_pubmed"].append((pmid, locs[:3]))
    for loc, name, year, pmid in labels:
        meta = pmids_meta.get(pmid)
        if not isinstance(meta, dict) or meta.get("error"):
            continue
        surname = _fold(_surname_of(str(meta.get("first_author", ""))))
        label = _fold(name)
        name_ok = bool(label) and bool(surname) and (label == surname or surname.endswith(label) or label.endswith(surname))
        try:
            year_ok = abs(int(year) - int(str(meta.get("pubdate_year", "")))) <= 1
        except (ValueError, TypeError):
            year_ok = False
        if not (name_ok and year_ok):
            problems["label_mismatch"].append((loc, f"'{name} {year}'", f"PMID:{pmid} = {meta.get('first_author')} {meta.get('pubdate_year')} ({str(meta.get('title', ''))[:60]})"))
    return problems


def run(refresh: bool = False, prune: bool = False) -> int:
    where, labels = harvest()
    cache = load_cache()
    pmids_meta: dict[str, dict] = cache.setdefault("pmids", {})
    if refresh:
        missing = sorted(set(where) - set(pmids_meta) | {p for p, m in pmids_meta.items() if m.get("error")})
        if missing:
            pmids_meta.update(fetch_esummary(missing))
        cache["_meta"] = {
            "source": "NCBI E-utilities esummary (db=pubmed)",
            "retrieved": date.today().isoformat(),
            "purpose": "offline snapshot for tests/test_citations_resolve.py; bibliographic facts only (first author, year, journal, title)",
        }
        prune = True
    if prune:
        for stale in sorted(set(pmids_meta) - set(where)):
            del pmids_meta[stale]
        CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        CACHE_PATH.write_text(json.dumps(cache, indent=1, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    problems = check(where, labels, pmids_meta)
    n_labelled = len({(loc, pmid) for loc, _, _, pmid in labels})
    print(f"citation gate: {len(where)} unique PMIDs in {len(scan_files())} files; {n_labelled} author/year-labelled occurrences checked")
    for kind, items in problems.items():
        if items:
            print(f"  FAIL {kind}: {len(items)}")
            for item in items[:40]:
                print("     ", *item)
    if not any(problems.values()):
        print("  PASS: every cited PMID resolves, and every adjacent author/year label matches PubMed")
        return 0
    if problems["not_in_cache"]:
        print("  (new PMID? run: python3 tests/test_citations_resolve.py --refresh-cache)")
    return 1


def test_every_cited_pmid_resolves_and_matches_its_label() -> None:  # pytest entry point
    if run() != 0:
        raise AssertionError("Citation resolution gate failed")


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify that every PMID cited in prose/code resolves and matches its author/year label.")
    parser.add_argument("--refresh-cache", action="store_true", help="fetch missing PMIDs from NCBI E-utilities and rewrite the cache")
    parser.add_argument("--prune-cache", action="store_true", help="offline: drop cache entries no scanned file references")
    args = parser.parse_args()
    sys.exit(run(refresh=args.refresh_cache, prune=args.prune_cache))


if __name__ == "__main__":
    main()
