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
        "studies": 1151,
        "arms": 7659,
        "outcomes": 14105,
        "edges": 2447,
        "effectors": 20,
        "graveyard_index": 10,
    }

    funnel = ceg["screening_funnel"]
    assert funnel["screened_records_total"] == 3922
    assert funnel["empirical_candidates_advanced_total"] == 3246
    assert funnel["wave1_primary_empirical_extracted_total"] == 1151
    assert (
        funnel["wave1_ring1_empirical_core_extracted"] + funnel["wave1_new_empirical_studies_extracted"]
        == 1151
    )

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
    assert atlas["summary"]["primary_empirical_studies_extracted"] == 1151
    assert atlas["summary"]["experimental_and_control_arms_total"] == 7659
    assert atlas["summary"]["structured_outcomes_total"] == 14105
    assert atlas["summary"]["blueprint_links_total"] == 2447
    assert atlas["summary"]["ring1_scoring_baseline_studies_unchanged"] == 563
    assert atlas["summary"]["inline_verbatim_quotes_shipped"] == 0
    assert atlas["summary"]["unresolved_outcome_citations"] == 0

    # Manifest alignment when manifest.json is present
    manifest_path = REPO_ROOT / "manifest.json"
    if manifest_path.is_file():
        mf = json.loads(manifest_path.read_text(encoding="utf-8"))
        mf_counts = mf["counts"]
        assert mf_counts["extracted_studies_total"] == 563
        assert mf_counts["extracted_atomic_effects_total"] == 3094
        assert mf_counts["ceg_studies_total"] == 1151
        assert mf_counts["ceg_arms_total"] == 7659
        assert mf_counts["ceg_outcomes_total"] == 14105
        assert mf_counts["ceg_blueprint_links_total"] == 2447
