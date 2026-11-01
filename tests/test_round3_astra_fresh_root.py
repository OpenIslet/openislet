"""Source-identity mutation probes at the acquisition-to-dossier boundary."""
import json

import pytest

from scripts import build_study_dossiers as dossiers
from scripts import harvest_corpus as harvest


RECORD = {"pmid": "38508139", "pmcid": "PMC10983038", "doi": "10.1016/j.xcrm.2024.101450"}


def article(ids, extra=""):
    return ("<article><front><article-meta>" + ids + "</article-meta></front>"
            "<body><sec><title>Results</title><p>A sufficiently long fixture result describes the source being acquired.</p></sec></body>"
            + extra + "</article>").encode()


@pytest.mark.parametrize("ids,extra", [
    ('<article-id pub-id-type="pmid">99999999</article-id>', ""),
    ('<article-id pub-id-type="pmc">PMC9999999</article-id>', ""),
    ('<article-id pub-id-type="pmid">38508139</article-id><article-id pub-id-type="pmc">9999999</article-id>', ""),
    ('<article-id pub-id-type="pmid">38508139</article-id><article-id pub-id-type="pmid">99999999</article-id>', ""),
    ('<article-id pub-id-type="pmid">38508139</article-id><article-id pub-id-type="doi">10.1234/unrelated</article-id>', ""),
    ("", '<back><ref-list><ref><article-id pub-id-type="pmid">38508139</article-id></ref></ref-list></back>'),
    ("", '<sub-article><front><article-meta><article-id pub-id-type="pmid">38508139</article-id></article-meta></front></sub-article>'),
])
@pytest.mark.parametrize("existing", [False, True])
def test_fetched_wrong_article_cannot_be_published_or_replace_cache(tmp_path, monkeypatch, ids, extra, existing):
    monkeypatch.setattr(dossiers, "DOSSIER_DIR", tmp_path)
    monkeypatch.setattr(dossiers, "fetch_epmc_fulltext_xml", lambda _: article(ids, extra))
    target = tmp_path / "38508139.json"
    before = json.dumps({**RECORD, "abstract": "Original authenticated abstract", "has_parsed_fulltext_xml": False})
    if existing:
        target.write_text(before)
    with pytest.raises(ValueError, match="JATS"):
        dossiers.build_one_dossier("38508139", RECORD, [], force=True)
    assert target.read_text() == before if existing else not target.exists()


@pytest.mark.parametrize("kind,value", [("pmid", "38508139"), ("pmc", "10983038"), ("pmcid", "PMC10983038"), ("pmc-uid", "10983038")])
def test_agreeing_fetched_article_is_published(tmp_path, monkeypatch, kind, value):
    monkeypatch.setattr(dossiers, "DOSSIER_DIR", tmp_path)
    source = article(f'<article-id pub-id-type="{kind}">{value}</article-id><article-id pub-id-type="doi">https://doi.org/10.1016/J.XCRM.2024.101450</article-id>')
    monkeypatch.setattr(dossiers, "fetch_epmc_fulltext_xml", lambda _: source)
    assert dossiers.build_one_dossier("38508139", RECORD, []) == ("38508139", True, 1)
    data = json.loads((tmp_path / "38508139.json").read_text())
    assert data["fulltext_data_sections"]["results_paragraphs"] == ["A sufficiently long fixture result describes the source being acquired."]


def test_namespaced_article_identity_is_supported():
    source = b'<article xmlns="urn:jats"><front><article-meta><article-id pub-id-type="pmid">38508139</article-id></article-meta></front></article>'
    dossiers.validate_jats_source_identity(source, "38508139", RECORD)


@pytest.mark.parametrize("source", [b"invalid XML", b"<article/>", b'<collection><article><front><article-meta><article-id pub-id-type="pmid">38508139</article-id></article-meta></front></article></collection>'])
def test_malformed_or_unbound_fetched_source_is_rejected(source):
    with pytest.raises(ValueError, match="JATS"):
        dossiers.validate_jats_source_identity(source, "38508139", RECORD)


@pytest.mark.parametrize("force", [False, True])
@pytest.mark.parametrize("key,value", [("pmcid", "PMC9999999"), ("doi", "10.1234/unrelated")])
def test_cached_source_cannot_be_relabelled_by_conflicting_metadata(tmp_path, monkeypatch, force, key, value):
    monkeypatch.setattr(dossiers, "DOSSIER_DIR", tmp_path)
    target = tmp_path / "38508139.json"
    before = json.dumps({**RECORD, "has_parsed_fulltext_xml": True,
                         "fulltext_data_sections": {"results_paragraphs": ["Original source text"]}})
    target.write_text(before)
    def must_not_fetch(_):
        pytest.fail("Conflicting cached identity must be rejected before acquisition")
    monkeypatch.setattr(dossiers, "fetch_epmc_fulltext_xml", must_not_fetch)
    with pytest.raises(ValueError, match="conflicts with cached dossier source"):
        dossiers.build_one_dossier("38508139", {**RECORD, key: value}, [], force=force)
    assert target.read_text() == before


@pytest.mark.parametrize("doi", ["https://doi.org/10.1016/J.XCRM.2024.101450", "doi: 10.1016/j.xcrm.2024.101450"])
def test_equivalent_doi_allows_failed_refresh_to_preserve_source(tmp_path, monkeypatch, doi):
    monkeypatch.setattr(dossiers, "DOSSIER_DIR", tmp_path)
    monkeypatch.setattr(dossiers, "fetch_epmc_fulltext_xml", lambda _: None)
    target = tmp_path / "38508139.json"
    source = {"results_paragraphs": ["Original source text"]}
    target.write_text(json.dumps({**RECORD, "has_parsed_fulltext_xml": True, "fulltext_data_sections": source}))
    dossiers.build_one_dossier("38508139", {**RECORD, "doi": doi}, [], force=True)
    after = json.loads(target.read_text())
    assert after["fulltext_data_sections"] == source
    assert after["has_parsed_fulltext_xml"] is True


def test_missing_refresh_identifiers_preserve_cached_identity_and_source(tmp_path, monkeypatch):
    monkeypatch.setattr(dossiers, "DOSSIER_DIR", tmp_path)
    target = tmp_path / "38508139.json"
    source = {"results_paragraphs": ["Original source text"]}
    target.write_text(json.dumps({**RECORD, "has_parsed_fulltext_xml": True, "fulltext_data_sections": source}))
    dossiers.build_one_dossier("38508139", {"pmid": "38508139"}, [], force=True)
    after = json.loads(target.read_text())
    assert {key: after[key] for key in RECORD} == RECORD
    assert after["fulltext_data_sections"] == source


class JsonResponse:
    def __init__(self, payload):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self):
        return json.dumps(self.payload).encode()


def install_pages(monkeypatch, pages):
    iterator = iter(pages)
    monkeypatch.setattr(harvest.urllib.request, "urlopen", lambda *args, **kwargs: JsonResponse(next(iterator)))
    monkeypatch.setattr(harvest.time, "sleep", lambda _: None)


@pytest.mark.parametrize("name", ["europepmc_search_ids", "europepmc_search_patents"])
@pytest.mark.parametrize("next_cursor", ["next", "*", None])
def test_discovery_cannot_return_incomplete_europepmc_pages(monkeypatch, name, next_cursor):
    install_pages(monkeypatch, [{"hitCount": 2, "resultList": {"result": [{"pmid": "123", "id": "123", "source": "PAT"}]}, "nextCursorMark": next_cursor}])
    with pytest.raises(ValueError, match="Europe PMC"):
        getattr(harvest, name)("synthetic", max_pages=1)


@pytest.mark.parametrize("name", ["europepmc_search_ids", "europepmc_search_patents"])
def test_europepmc_exact_budget_completion_keeps_all_results(monkeypatch, name):
    rows = [{"pmid": str(i), "id": str(i), "source": "PAT"} for i in (123, 124)]
    install_pages(monkeypatch, [
        {"hitCount": 2, "resultList": {"result": rows[:1]}, "nextCursorMark": "next"},
        {"hitCount": 2, "resultList": {"result": rows[1:]}, "nextCursorMark": "last"},
    ])
    result = getattr(harvest, name)("synthetic", max_pages=2)
    assert (result[0] if name.endswith("ids") else {item["id"] for item in result}) == {"123", "124"}


@pytest.mark.parametrize("name", ["europepmc_search_ids", "europepmc_search_patents"])
@pytest.mark.parametrize("budget", [0, -1, True, 1.5])
def test_discovery_rejects_invalid_page_budgets(name, budget):
    with pytest.raises(ValueError, match="page budget"):
        getattr(harvest, name)("synthetic", max_pages=budget)


def test_europepmc_count_change_and_empty_tail_are_not_completion(monkeypatch):
    first = {"hitCount": 2, "resultList": {"result": [{"pmid": "123"}]}, "nextCursorMark": "next"}
    for tail in [{"hitCount": 3, "resultList": {"result": []}}, {"hitCount": 2, "resultList": {"result": []}}]:
        install_pages(monkeypatch, [first, tail])
        with pytest.raises(ValueError, match="Europe PMC"):
            harvest.europepmc_search_ids("synthetic", max_pages=2)


def test_europepmc_complete_empty_search_remains_valid(monkeypatch):
    install_pages(monkeypatch, [{"hitCount": 0, "resultList": {"result": []}}])
    assert harvest.europepmc_search_ids("synthetic") == (set(), [])


def test_geo_truncated_search_fails_before_metadata_requests(monkeypatch):
    install_pages(monkeypatch, [{"esearchresult": {"count": "2001", "idlist": [str(i) for i in range(1, 2001)]}}])
    with pytest.raises(ValueError, match="GEO search retrieval is incomplete"):
        harvest.geo_gds_search_datasets({"synthetic": "fixture"})


def test_geo_complete_search_preserves_matched_dataset(monkeypatch):
    install_pages(monkeypatch, [
        {"esearchresult": {"count": "1", "idlist": ["123"]}},
        {"result": {"123": {"title": "Human islet fixture", "accession": "GSE123", "pubmedids": [456]}}},
    ])
    datasets, pmids = harvest.geo_gds_search_datasets({"synthetic": "fixture"})
    assert len(datasets) == 1 and datasets[0]["accession"] == "GSE123" and pmids == {"456"}


@pytest.mark.parametrize("name", ["europepmc_search_ids", "europepmc_search_patents"])
def test_europepmc_new_cursors_cannot_hide_repeated_record_identity(monkeypatch, name):
    row = {"pmid": "123", "id": "123", "source": "PAT"}
    install_pages(monkeypatch, [
        {"hitCount": 2, "resultList": {"result": [row]}, "nextCursorMark": "first"},
        {"hitCount": 2, "resultList": {"result": [{**row, "title": "Changed metadata"}]}, "nextCursorMark": "second"},
    ])
    with pytest.raises(ValueError, match="repeated a record identity"):
        getattr(harvest, name)("synthetic", max_pages=2)


@pytest.mark.parametrize("summary", [{}, {"123": {"error": "Document summary unavailable"}}])
def test_geo_missing_metadata_is_not_a_relevance_exclusion(monkeypatch, summary):
    install_pages(monkeypatch, [
        {"esearchresult": {"count": "1", "idlist": ["123"]}},
        {"result": summary},
    ])
    with pytest.raises(ValueError, match="GEO metadata response is missing or rejected"):
        harvest.geo_gds_search_datasets({"synthetic": "fixture"})


def test_complete_europepmc_discovery_can_filter_reviews(monkeypatch):
    install_pages(monkeypatch, [{"hitCount": 2, "resultList": {"result": [
        {"source": "MED", "id": "123", "pmid": "123", "pubType": "Review"},
        {"source": "MED", "id": "124", "pmid": "124", "pubType": "Journal Article"},
    ]}}])
    assert harvest.europepmc_search_ids("synthetic", max_pages=1) == ({"124"}, [])


@pytest.mark.parametrize("pmid", [True, False, 123, "0123", "１２３", "not-a-pmid"])
def test_europepmc_does_not_coerce_invalid_pmid_types(monkeypatch, pmid):
    install_pages(monkeypatch, [{"hitCount": 1, "resultList": {"result": [{"source": "MED", "id": "123", "pmid": pmid}]}}])
    with pytest.raises(ValueError, match="invalid PMID identity"):
        harvest.europepmc_search_ids("synthetic")


@pytest.mark.parametrize("payload", [[], {"hitCount": 0, "resultList": []}])
def test_europepmc_malformed_containers_fail_at_provider_boundary(monkeypatch, payload):
    install_pages(monkeypatch, [payload])
    with pytest.raises(ValueError, match="invalid result container"):
        harvest.europepmc_search_ids("synthetic")
