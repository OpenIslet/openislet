#!/usr/bin/env python3
"""
Generate the Wave-0 / P0-8 Release & Provenance Manifest (`manifest.json`) and
synchronize the three brief-referenced JSONL/JSON artifact aliases (`AT-00`):
  - data/extracted/ring1_study_dossiers_v0.5.jsonl
  - data/extracted/ring1_atomic_effects_v0.5.jsonl
  - data/extracted/islet1_stack_ranking_v0.5.json
"""

from __future__ import annotations

import hashlib
import json
import os
import stat
import fnmatch
import threading
from functools import lru_cache
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXTRACTED_DIR = ROOT / "data" / "extracted"
DOSSIERS_DIR = ROOT / "data" / "dossiers"
OMICS_DIR = ROOT / "data" / "omics"
RAW_DIR = ROOT / "data" / "raw"
GRAVEYARD_PATH = ROOT / "data" / "graveyard" / "idea_graveyard_v0.5.json"
_PUBLIC_POLICY_SNAPSHOT = ContextVar('openislet_public_policy_snapshot', default=None)
_GENERATION_LEASE = ContextVar('openislet_manifest_generation_lease', default=None)
_GENERATION_JOURNAL = ContextVar('openislet_manifest_generation_journal', default=None)
_EXPORT_PROCESS_LOCK = threading.RLock()
_NEVER_PUBLIC = (
    'data/dossiers/*', 'data/curated/*',
    # DEC-011 keeps detailed downstream plans private. The decisions log also
    # embeds those plans; neither file is a public blueprint artifact.
    'docs/OPENISLET-2026-2030-MASTER-TIMELINE.md',
    'docs/PROJECT-DECISIONS-LOG.md',
    'data/raw/alphagenome_avi_3118_snvs_full_archive.json',
    'data/raw/alphagenome_avi_26_seed_anchors.json',
    'data/raw/alphagenome_avi_live_receipt_2026_10_28.json',
    'data/raw/private_export_guard_patterns.json',
    'data/extracted/islet1_candidate_stack_rankings.json',
    'data/extracted/islet1_stack_ranking_v0.5.json',
    'CLAUDE.md',
    '.git', '.git/*', '.claude', '.claude/*', '.claude/**', '.env', '.env.*',
)


PUBLIC_PACKAGE_INCLUDED_PATHS: tuple[str, ...] = (
    "LICENSE",
    "LICENSE-DATA",
    "CODE_OF_CONDUCT.md",
    "CONTRIBUTING.md",
    "SECURITY.md",
    ".gitignore",
    ".github/openislet-mark.svg",
    ".github/openislet-mark-dark.svg",
    ".github/ISSUE_TEMPLATE/config.yml",
    ".github/ISSUE_TEMPLATE/files-problem.yml",
    "requirements.txt",
    "pyproject.toml",
    ".licenses/*.txt",
    "ISLET-1-BLUEPRINT.md",
    "README.md",
    "manifest.json",
    "docs/*.md",
    "schemas/extraction_schema_v0.2.json",
    "scripts/*.py",
    "tests/*.py",
    "data/raw/ring1_platinum_queue.json",
    "data/raw/corpus_manifest_all_sources.jsonl",
    "data/raw/corpus_stratification_manifest.json",
    "data/raw/doc_citation_metadata_cache.json",
    "data/raw/layer2_genomic_and_jaspar_reference.json",
    "data/raw/layer5_dbsnp_clinvar_reference.json",
    "data/constructs/*.fasta",
    "data/extracted/*.json",
    "data/extracted/*.jsonl",
    "data/graveyard/idea_graveyard_v0.5.json",
    "data/omics/*.json",
    "data/protocols/*.json",
    "data/structures/*.cif",
    "data/structures/*.pdb",
    "data/structures/*.json",
    "data/structures/*.md",
    "data/structures/pymol/*",
)


def _forbidden_public_path(relative):
    # Match portable identities too: a case alias can address the same private
    # file on APFS even when the manifest spells its name differently.
    key = relative.casefold()
    return any(part in {'__pycache__','.pytest_cache','.claude','.git'} or part == '.env' or part.startswith('.env.') for part in key.split('/')) or any(fnmatch.fnmatchcase(key, pattern.casefold()) for pattern in _NEVER_PUBLIC)


def verify_public_bundle_privacy_guard(
    bundle_root: Path,
    included_paths: list[str] | tuple[str, ...] | None = None,
    *,
    source_root: Path | None = None,
    require_private_pattern_file: bool | None = None,
) -> dict[str, int | bool]:
    """Fail-closed privacy and credential scan over a staged public bundle or allowlisted file set."""
    import re as _re
    base = source_root if source_root is not None else ROOT
    if require_private_pattern_file is None:
        require_private_pattern_file = not (base / ".openislet_public_bundle").exists()
    guard_cfg_path = base / "data" / "raw" / "private_export_guard_patterns.json"
    if not guard_cfg_path.is_file() and require_private_pattern_file:
        raise RuntimeError(
            "PRIVACY GUARD FAILED: Missing required private policy file data/raw/private_export_guard_patterns.json"
        )
    raw_patterns: dict[str, str] = {
        "home_and_temp_paths": (
            r"(?<![A-Za-z0-9._-])(?:/" + r"Users/|/" + r"home/[a-z0-9_-]+/|/"
            + r"private/tmp|/" + r"var/folders|~/\.?[A-Za-z0-9_-]{2,}/)"
        ),
        "private_keys_and_tokens": (
            r"(?:-----BEGIN [A-Z0-9 ]*PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{16,}|github_pat_[A-Za-z0-9_]{16,}|sk-ant-[A-Za-z0-9_-]{8,}|sk-proj-[A-Za-z0-9_-]{16,}|sk-svcacct-[A-Za-z0-9_-]{16,}|(?<![A-Za-z0-9_-])sk-[A-Za-z0-9]{20,}|AIza[0-9A-Za-z_-]{20,}|hf_[A-Za-z0-9]{20,}|xox[baprs]-[A-Za-z0-9-]{10,}|AKIA[0-9A-Z]{16})"
        ),
    }
    allowed_email_domains: set[str] = {"openislet.org", "users.noreply.github.com"}
    allowed_file_emails: dict[str, set[str]] = {
        "data/structures/pdb_00003cdg.cif": {"sw-help" + "@rcsb.rutgers.edu", "PDAdams" + "@lbl.gov"},
    }
    cfg = json.loads(guard_cfg_path.read_text(encoding="utf-8")) if guard_cfg_path.is_file() else {}
    for k, v in (cfg.get("regex_patterns") or {}).items():
        if isinstance(k, str) and isinstance(v, str):
            raw_patterns[k] = v
    for dom in cfg.get("allowed_email_domains") or ():
        if isinstance(dom, str):
            allowed_email_domains.add(dom.lower())
    for rel_k, em_list in (cfg.get("allowed_file_specific_emails") or {}).items():
        if isinstance(rel_k, str) and isinstance(em_list, list):
            allowed_file_emails.setdefault(rel_k, set()).update(str(e) for e in em_list)
    compiled_patterns = {k: _re.compile(v) for k, v in raw_patterns.items()}
    home_name = Path.home().name
    if len(home_name) >= 4 and home_name.lower() not in ("root", "admin", "user", "runner", "ubuntu"):
        compiled_patterns["local_account_name"] = _re.compile(_re.escape(home_name), _re.I)
    email_re = _re.compile(r"\b([A-Za-z0-9._%+-]+@([A-Za-z0-9.-]+\.[A-Za-z]{2,}))\b")
    candidate_files: list[tuple[Path, bool]] = []
    if included_paths is None:
        candidate_files = [(fp, False) for fp in sorted(bundle_root.rglob("*")) if fp.is_file()]
    else:
        seen_paths: set[Path] = set()
        for pat in included_paths:
            is_glob = any(ch in pat for ch in ("*", "?", "["))
            for fp in sorted(bundle_root.glob(pat)):
                if fp.is_file() and fp not in seen_paths:
                    seen_paths.add(fp)
                    candidate_files.append((fp, is_glob))
    files_scanned = 0
    for fp, from_glob in candidate_files:
        rel = fp.relative_to(bundle_root).as_posix()
        if rel != ".openislet_public_bundle" and _forbidden_public_path(rel):
            if from_glob:
                continue
            raise RuntimeError(f"PRIVACY GUARD FAILED: Public export privacy guard blocked forbidden path: {rel}")
        is_binary = rel.lower().endswith((".png", ".pse", ".jpg", ".jpeg", ".gif", ".webp", ".ico"))
        if is_binary:
            text = b"\n".join(_re.findall(rb"[ -~]{6,}", fp.read_bytes())).decode("ascii", errors="ignore")
        else:
            text = fp.read_text(encoding="utf-8", errors="ignore")
        for label, rx in compiled_patterns.items():
            if rx.search(text):
                raise RuntimeError(f"PRIVACY GUARD FAILED: Public export privacy guard blocked {label} in {rel}")
        if not is_binary:
            file_allowed = {e.lower() for e in allowed_file_emails.get(rel, set())}
            for full_email, domain in email_re.findall(text):
                if domain.lower() not in allowed_email_domains and full_email.lower() not in file_allowed:
                    raise RuntimeError(f"PRIVACY GUARD FAILED: Public export privacy guard blocked unauthorized email in {rel}")
        files_scanned += 1
    return {"passed": True, "violation_count": 0, "files_scanned": files_scanned}


def _canonical_relative_path(value: str) -> bool:
    import re
    return isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9._/-]+", value) is not None and all(part not in ("", ".", "..") for part in value.split("/"))


def _artifact_inventory_errors(artifacts) -> list[str]:
    """Validate even retained private pins before any path is dereferenced."""
    if not isinstance(artifacts, dict) or not artifacts:
        return ['invalid-or-empty-artifact-inventory']
    errors = []
    identities = set()
    for relative, pin in artifacts.items():
        if not _canonical_relative_path(relative) or relative in ('manifest.json', '.openislet_public_bundle'):
            errors.append('invalid-artifact-path:' + str(relative))
            continue
        identity = relative.casefold()
        if identity in identities:
            errors.append('duplicate-portable-artifact:' + relative)
        identities.add(identity)
        if (not isinstance(pin, dict) or not isinstance(pin.get('sha256'), str)
                or _re.fullmatch('[0-9a-f]{64}', pin['sha256']) is None
                or type(pin.get('bytes')) is not int or pin['bytes'] < 0):
            errors.append('invalid-artifact-pin:' + relative)
    return errors


def public_excerpt_limits(allow) -> tuple[int, float]:
    """Resolve one validated policy for exporter and acceptance verification."""
    import math
    characters = allow.get('non_oa_per_article_aggregate_verbatim_ceiling_characters', 240)
    share = allow.get('non_oa_per_article_max_abstract_share', .25)
    if type(characters) is not int or not 0 <= characters <= 240:
        raise ValueError('Public non-OA character policy may not exceed 240 characters')
    if isinstance(share, bool) or not isinstance(share, (int, float)) or not math.isfinite(share) or not 0 < share <= .25:
        raise ValueError('Public non-OA abstract share policy must be finite in (0, .25]')
    return characters, share


def _regular_file_bytes(root: Path, relative: str) -> bytes:
    """Read one single-link regular file without following release-tree links."""
    if not _canonical_relative_path(relative):
        raise ValueError("Release file path must be canonical and relative")
    current = root
    parts = relative.split("/")
    for i, part in enumerate(parts):
        current = current / part
        info = current.lstat()
        if stat.S_ISLNK(info.st_mode):
            raise ValueError("Release files and their parents may not be symbolic links")
        if i < len(parts)-1 and not stat.S_ISDIR(info.st_mode):
            raise ValueError("Release parent must be a directory")
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError("Release input must be a regular file with one hard link")
    fd = os.open(current, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        opened = os.fstat(fd)
        if not stat.S_ISREG(opened.st_mode) or opened.st_nlink != 1 or (opened.st_dev, opened.st_ino) != (info.st_dev, info.st_ino):
            raise ValueError("Release input identity changed while opening")
        with os.fdopen(fd, "rb", closefd=False) as stream:
            data = stream.read()
        after = os.fstat(fd)
        if (after.st_size, after.st_mtime_ns, after.st_ctime_ns) != (len(data), opened.st_mtime_ns, opened.st_ctime_ns):
            raise ValueError("Release input changed while reading")
        return data
    finally:
        os.close(fd)


def sha256_file(path: Path) -> str:
    base = ROOT if path.is_relative_to(ROOT) else path.parent
    return hashlib.sha256(_regular_file_bytes(base, path.relative_to(base).as_posix())).hexdigest()


@lru_cache(maxsize=1)
def _publication_support():
    """Load our sibling from this source location, independent of caller cwd.

    The connector loads this exporter by file spec, without a package or the
    science scripts on sys.path. Neither that caller path nor an unrelated
    module with our sibling's basename is authority for publication code.
    """
    import types
    path=Path(__file__).resolve().with_name('study_artifact_exports.py')
    source=_regular_file_bytes(path.parent,path.name)
    helper=types.ModuleType('_openislet_release_publication_support')
    helper.__file__=str(path)
    exec(compile(source,str(path),'exec'),helper.__dict__)
    return helper


@contextmanager
def _publication_lock(path: Path):
    with _publication_support().PublicationLease(path) as lease:yield lease


def _require_generation_lease():
    lease=_GENERATION_LEASE.get()
    if lease is not None:lease.require_owned()


def main() -> None:
    """Serialize generators and restore their previous outputs on caught failure.

    The manifest is still a sequence of per-file replacements. This rollback
    does not assert crash atomicity or snapshot isolation for concurrent readers.
    """
    with _publication_lock(ROOT/'.openislet-manifest.lock') as lease:
        journal=_publication_support().PublicationJournal()
        token=_GENERATION_LEASE.set(lease)
        journal_token=_GENERATION_JOURNAL.set(journal)
        try:
            _generate_release_manifest();lease.require_owned()
        except BaseException:
            if not lease.owns_lock():raise
            journal.rollback(lease.require_owned)
            raise
        finally:
            _GENERATION_JOURNAL.reset(journal_token)
            _GENERATION_LEASE.reset(token)


def _publish_generation_artifacts(payloads,lock_path):
    _require_generation_lease()
    _publication_support().publish_byte_artifacts(payloads,lock_path=lock_path,rollback_journal=_GENERATION_JOURNAL.get())


def _generate_release_manifest() -> None:
    from wetlab_contract import derive_contract, CONTRACT_RELATIVE_PATH
    from render_documentation import write_documentation
    _require_generation_lease()
    contract_bytes=(json.dumps(derive_contract(ROOT),indent=2,allow_nan=False)+'\n').encode('utf-8')
    _publish_generation_artifacts({ROOT/CONTRACT_RELATIVE_PATH:contract_bytes},ROOT/'.contract-publication.lock')
    _require_generation_lease()
    write_documentation(ROOT, manifest_pin_count=0, render=False,rollback_journal=_GENERATION_JOURNAL.get())
    study_files = sorted(
        [p for p in EXTRACTED_DIR.glob("*.json") if p.stem.isdigit()],
        key=lambda p: p.stem,
    )
    studies = [json.loads(p.read_text(encoding="utf-8")) for p in study_files]

    # Materialize producer exports through the same function used by extraction.
    from study_artifact_exports import write_study_exports
    _require_generation_lease()
    dossiers_jsonl, effects_jsonl, n_effects = write_study_exports(ROOT, studies,rollback_journal=_GENERATION_JOURNAL.get())

    # 3. Synchronize data/extracted/islet1_stack_ranking_v0.5.json
    canonical_ranking = EXTRACTED_DIR / "islet1_candidate_stack_rankings.json"
    brief_ranking = EXTRACTED_DIR / "islet1_stack_ranking_v0.5.json"
    ranking_payload = json.loads(canonical_ranking.read_text(encoding="utf-8"))
    ranking_payload["publication_freeze_status"] = {
        "outward_rankings_frozen": False,
        "status": "UNFROZEN_PRIOR_ANCHORED_BASELINE_v0.8.0",
        "classification": "PRIOR_ANCHORED_SCENARIO_OUTPUTS",
        "retired_chassis_pending_definition": ["Chassis C", "Chassis D"],
        "unfreeze_gate": "Explicit human approval from Kev after P0-1..P0-8 pass (UNFROZEN in v0.8.0 per DEC-003; precomputed ranking JSON caches remain excluded from public zip bundle allowlist while live rank_islet1_stacks.py is shipped)",
    }
    _require_generation_lease()
    ranking_bytes=(json.dumps(ranking_payload,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode('utf-8')
    _publish_generation_artifacts({canonical_ranking:ranking_bytes,brief_ranking:ranking_bytes},EXTRACTED_DIR/'.ranking-publication.lock')

    # 4. Build root manifest.json (Layers 1-5 and Provenance Baseline)
    tracked_artifacts = [
        "LICENSE",
        "LICENSE-DATA",
        "CODE_OF_CONDUCT.md",
        "CONTRIBUTING.md",
        "SECURITY.md",
        ".gitignore",
        ".github/openislet-mark.svg",
        ".github/openislet-mark-dark.svg",
        ".github/ISSUE_TEMPLATE/config.yml",
        ".github/ISSUE_TEMPLATE/files-problem.yml",
        "requirements.txt",
        "pyproject.toml",
        ".licenses/alphafold_database_fetch_and_analyze_LICENSE.txt",
        ".licenses/alphagenome_single_variant_analysis_LICENSE.txt",
        ".licenses/clinical_trials_database_LICENSE.txt",
        ".licenses/clinvar_database_LICENSE.txt",
        ".licenses/dbsnp_database_LICENSE.txt",
        ".licenses/encode_ccres_database_LICENSE.txt",
        ".licenses/ensembl_database_LICENSE.txt",
        ".licenses/foldseek_structural_search_LICENSE.txt",
        ".licenses/gnomad_database_LICENSE.txt",
        ".licenses/gtex_database_LICENSE.txt",
        ".licenses/human_protein_atlas_database_LICENSE.txt",
        ".licenses/jaspar_database_LICENSE.txt",
        ".licenses/literature_search_europepmc_LICENSE.txt",
        ".licenses/opentargets_database_LICENSE.txt",
        ".licenses/pdb_database_LICENSE.txt",
        ".licenses/pubmed_database_LICENSE.txt",
        ".licenses/pymol_LICENSE.txt",
        ".licenses/ucsc_conservation_and_tfbs_LICENSE.txt",
        ".licenses/uniprot_database_LICENSE.txt",
        "ISLET-1-BLUEPRINT.md",
        "README.md",
        "docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md",
        "docs/ISLET-1-MASTER-SPEC-v0.1.md",
        "docs/ISLET-1-SYNTHESIS-REPORT-v0.4.md",
        "docs/ISLET-1-CANDIDATE-BLUEPRINTS-AND-OPTIMIZER-v0.2.md",
        "docs/PROJECT-ISLET-1-RING1-DISCOVERY-ATLAS-AND-STRATEGY.md",
        "docs/PROJECT-DECISIONS-LOG.md",
        "docs/OPENISLET-2026-2030-MASTER-TIMELINE.md",
        "docs/STAGE-0-EXECUTION-PLAYBOOK.md",
        "data/raw/corpus_manifest_all_sources.jsonl",
        "data/raw/ring1_platinum_queue.json",
        "data/raw/corpus_stratification_manifest.json",
        "data/raw/doc_citation_metadata_cache.json",
        "data/raw/layer2_genomic_and_jaspar_reference.json",
        "data/raw/alphagenome_avi_3118_snvs_full_archive.json",
        "data/raw/alphagenome_avi_26_seed_anchors.json",
        "data/raw/alphagenome_avi_live_receipt_2026_10_28.json",
        "data/raw/layer5_dbsnp_clinvar_reference.json",
        "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta",
        "data/graveyard/idea_graveyard_v0.5.json",
        "data/omics/target_prior_table.json",
        "data/omics/t1d_associated_targets_top100.json",
        "data/omics/t1d_disease_search.json",
        "data/extracted/cross_study_discovery_graph.json",
        "data/extracted/islet1_candidate_stack_rankings.json",
        "data/extracted/islet1_stack_ranking_v0.5.json",
        "data/extracted/ring1_study_dossiers_v0.5.jsonl",
        "data/extracted/ring1_atomic_effects_v0.5.jsonl",
        "data/extracted/layer2_regulatory_and_perturbation_atlas.json",
        "data/extracted/layer3_structural_verification_atlas.json",
        "data/extracted/layer4_tissue_and_metabolic_twin_atlas.json",
        "data/extracted/layer5_personalized_patient_and_crispr_atlas.json",
        "data/extracted/omics_baseline_summary.json",
        "data/extracted/wetlab_bayesian_posterior_calibration_v0.8.0.json",
        "data/extracted/wetlab_bayesian_posterior_calibration.json",
        "data/protocols/hip_core_fractional_factorial_v0.6.2.json",
        "data/structures/AF-P10144-F1-metadata.json",
        "data/structures/AF-P10144-F1-model_v6.cif",
        "data/structures/AF-P10144-F1-predicted_aligned_error_v6.json",
        "data/structures/AF-P50453-F1-metadata.json",
        "data/structures/AF-P50453-F1-model_v6.cif",
        "data/structures/AF-P50453-F1-predicted_aligned_error_v6.json",
        "data/structures/1EFX.pdb",
        "data/structures/1IM9.pdb",
        "data/structures/pdb_00001iau.cif",
        "data/structures/pdb_00001mhe.cif",
        "data/structures/pdb_00002jjs.cif",
        "data/structures/pdb_00003bze.cif",
        "data/structures/pdb_00003cdg.cif",
        "data/structures/pdb_00004cmm.cif",
        "data/structures/pdb_00004kjy.cif",
        "data/structures/pdb_00004zqk.cif",
        "data/structures/pdb_metadata.json",
        "data/structures/pymol/layer3_hlae_cd94_nkg2a_3cdg.png",
        "data/structures/pymol/openislet_layer3_all_complexes_v0.7.0.pse",
        "data/structures/serpinb9_foldseek_results.md",
        "schemas/extraction_schema_v0.2.json",
        "scripts/study_artifact_exports.py",
        "scripts/engineering_priors.py",
        "scripts/product_contract.py",
        "scripts/public_excerpt_policy.py",
        "scripts/wetlab_contract.py",
        "scripts/render_documentation.py",
        "data/protocols/islet1_product_definitions_v0.8.0.json",
        "data/protocols/wetlab_endpoint_contract_v0.8.0.json",
        "data/protocols/stage0_decision_tree_v1.json",
        "data/extracted/documentation_facts_v0.8.0.json",
        "tests/test_preastra_layer1.py",
        "tests/test_preastra_molecular.py",
        "tests/test_preastra_physiology_compiler.py",
        "tests/test_preastra_root.py",
        "tests/test_master_molecular_validation.py",
        "tests/test_master_physiology_audit.py",
        "tests/test_round2_layer1.py",
        "tests/test_round2_molecular.py",
        "tests/test_round2_physiology.py",
        "tests/test_round2_root.py",
        "tests/test_round2_astra_root.py",
        "tests/test_round2_astra_layer1.py",
        "tests/test_round2_astra_molecular.py",
        "tests/test_round2_astra_physiology.py",
        "scripts/build_cross_study_graph.py",
        "scripts/build_study_dossiers.py",
        "scripts/extract_ring1_corpus.py",
        "scripts/generate_release_manifest.py",
        "scripts/harvest_corpus.py",
        "scripts/identify_and_stratify_corpus.py",
        "scripts/ingest_wetlab_results.py",
        "scripts/rank_islet1_stacks.py",
        "scripts/run_layer2_perturbation_and_regulatory_plane.py",
        "scripts/run_layer3_structural_binder_plane.py",
        "scripts/run_layer4_tissue_and_metabolic_twin.py",
        "scripts/run_layer5_personalized_patient_compiler.py",
        "scripts/summarize_omics_baseline.py",
        "tests/test_adversarial_acceptance.py",
        "tests/test_citations_resolve.py",
        "tests/test_design_algebra_exact.py",
        "tests/test_layer2_and_layer3_planes.py",
        "tests/test_layer4_and_layer5_planes.py",
    ]
    artifact_entries = {}
    for rel in tracked_artifacts:
        p = ROOT / rel
        if not p.exists():
            raise ValueError(f"Required release artifact is missing: {rel}")
        artifact_entries[rel] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
        }
    for p in study_files:
        rel = str(p.relative_to(ROOT))
        artifact_entries[rel] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
        }

    gnomad_files = sorted(OMICS_DIR.glob("gnomad_*.json"))
    for p in gnomad_files:
        rel = str(p.relative_to(ROOT))
        artifact_entries[rel] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
        }

    dossier_files = sorted(DOSSIERS_DIR.glob("*.json"))
    for p in dossier_files:
        rel = str(p.relative_to(ROOT))
        artifact_entries[rel] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
        }

    curated_files = sorted((ROOT / "data" / "curated").glob("*.json"))
    for p in curated_files:
        rel = str(p.relative_to(ROOT))
        artifact_entries[rel] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
        }
    for p in sorted((ROOT/'tests').glob('test_round3_*.py')):
        rel=p.relative_to(ROOT).as_posix()
        artifact_entries[rel]={"sha256":sha256_file(p),"bytes":p.stat().st_size}
    graveyard_data = json.loads(GRAVEYARD_PATH.read_text(encoding="utf-8"))
    graph_data = json.loads(
        (EXTRACTED_DIR / "cross_study_discovery_graph.json").read_text(encoding="utf-8")
    )

    n_ft_dossiers = 0
    for p in dossier_files:
        d_obj = json.loads(p.read_text(encoding="utf-8"))
        if d_obj.get("has_parsed_fulltext_xml") or d_obj.get("has_fulltext_xml"):
            n_ft_dossiers += 1

    layer4_atlas = json.loads((EXTRACTED_DIR / "layer4_tissue_and_metabolic_twin_atlas.json").read_text(encoding="utf-8"))
    n_layer4_pde_regimes = len(
        (layer4_atlas.get("part_a_336h_spatial_reaction_diffusion_pouch_pde") or {}).get("regimes_simulated") or {}
    )
    if n_layer4_pde_regimes == 0:
        raise ValueError("Layer 4 must contain generated PDE regimes")
    layer2_atlas = json.loads((EXTRACTED_DIR / "layer2_regulatory_and_perturbation_atlas.json").read_text())
    layer3_atlas = json.loads((EXTRACTED_DIR / "layer3_structural_verification_atlas.json").read_text())
    layer5_atlas = json.loads((EXTRACTED_DIR / "layer5_personalized_patient_and_crispr_atlas.json").read_text())

    manifest = {
        "schema_version": "openislet_release_manifest_v0.8.0",
        "release_date": "2026-11-01",
        "release_stage": "v0.8.0 Complete 5-Layer Stack & Provenance Baseline (2026-11-01 Public Release)",
        "publication_freeze": {
            "status": "UNFROZEN_PRIOR_ANCHORED_BASELINE_v0.8.0",
            "outward_rankings_frozen": False,
            "retired_chassis": ["Chassis C", "Chassis D"],
            "pilot_experiment_gate": "HIP-core 36-run two-stage bridged, replicated & lead-chassis-confirmatory fractional-factorial protocol specified (requires biological replicate blocking and human sign-off before execution)",
            "mcp_mode": "MODE_A_READ_AND_IIP_PROPOSAL_QUEUE_OPEN_V0_8_0",
        },
        "preprint_publication_deduplication_registry": {
            "40832357": {
                "preprint_doi": "10.1101/2025.08.12.669867",
                "superseded_by_pmid": "41202120",
                "published_doi": "10.1126/sciadv.adx9755",
                "source_text_type": "Preprint_Superseded_By_Peer_Reviewed",
                "quantitative_evidence_weight": 0.0,
            }
        },
        "auxiliary_non_ring1_dossiers_registry": {
            "note": "Two Ring-2 review dossiers cached in data/dossiers/ are outside the 563-study Ring-1 primary experimental extracted corpus (563 + 2 = 565 total dossiers):",
            "26649319": "Akter et al. 2016 (Journal of Diabetes Research, doi:10.1155/2016/2798269): 'Islet Amyloid Polypeptide: Structure, Function, and Pathophysiology.' (Ring-2 review dossier retained in local dossier cache; not part of 563-study Ring-1 extracted queue)",
            "35650334": "Nekoua, Alidjinou & Hober 2022 (Nature Reviews Endocrinology, doi:10.1038/s41574-022-00688-1): 'Persistent coxsackievirus B infection and pathogenesis of type 1 diabetes mellitus.' (Ring-2 review dossier retained in local dossier cache; not part of 563-study Ring-1 extracted queue)",
        },
        "public_redistribution_allowlist_spec": {
            "policy_note": (
                "Separates the local research audit repository tree (which caches data/dossiers/*.json and frozen scenario-ranking JSONs locally for SHA-256 and character-offset verification) "
                "from the rights-aware public distribution package allowlist. Availability in PMC does not by itself grant bulk full-text redistribution rights "
                "for non-CC-BY articles; public export bundles exclude data/dossiers/*.json and frozen outward ranking files while retaining their SHA-256 verification digests in manifest.json, "
                "enforce a strict per-article aggregate verbatim quote ceiling (<=240 chars and <25% of abstract length) across all string fields for all studies lacking an explicit machine-verified CC-BY/CC0 license tag, "
                "and include all third-party database license notices (.licenses/*.txt), Layers 2-5 constructs/structures/atlases, and required raw queue/manifest inputs."
            ),
            "cc_by_4_0_license_scope": (
                "Applies strictly to OpenIslet-authored structured annotations, coded arms, polarity registries, effect sizes, character-span offsets, "
                "SHA-256 digests, and cross-study graphs; third-party verbatim quotes are carved out in LICENSE and redacted above the 240-char / <25% per-article cap "
                "for all studies without an explicit CC-BY/CC0 license tag in --export-public-bundle."
            ),
            "excerpt_ceiling_characters": 350,
            "non_oa_per_article_aggregate_verbatim_ceiling_characters": 240,
            "non_oa_per_article_max_abstract_share": 0.25,
            "excerpt_ceiling_legal_note": "350-character per-effect quote ceiling (and <=240-character / <25%-abstract per-article aggregate ceiling in public bundles) is an internal repository design limit for factual attribution, not a claim of a universal statutory copyright safe harbor.",
            "source_rights_and_redistribution_classification_summary": {
                "total_local_dossiers_indexed": len(dossier_files),
                "pmc_xml_parsed_dossiers_local_cache_only": n_ft_dossiers,
                "pubmed_abstract_and_metadata_only_dossiers": len(dossier_files) - n_ft_dossiers,
                "full_dossier_cache_redistribution_status": "LOCAL_AUDIT_CACHE_ONLY_EXCLUDED_FROM_PUBLIC_BUNDLE",
                "extracted_study_json_excerpt_policy": "Local audit tree: factual metadata + short verbatim evidence quotes (<=350 chars) + SHA-256 pointer to source dossier. Exported public bundle: all studies lacking a machine-verified CC-BY/CC0 license field are capped at <=240 aggregate verbatim chars (<25% of abstract) across all string fields with excess quotes replaced by [REDACTED_NON_OA_EXCERPT: see PMID + source_span + quote_sha256].",
            },
            "public_package_included_paths": [
                "LICENSE",
                "LICENSE-DATA",
                "CODE_OF_CONDUCT.md",
                "CONTRIBUTING.md",
                "SECURITY.md",
                ".gitignore",
                ".github/openislet-mark.svg",
                ".github/openislet-mark-dark.svg",
                ".github/ISSUE_TEMPLATE/config.yml",
                ".github/ISSUE_TEMPLATE/files-problem.yml",
                "requirements.txt",
                "pyproject.toml",
                ".licenses/*.txt",
                "ISLET-1-BLUEPRINT.md",
                "README.md",
                "manifest.json",
                "docs/*.md",
                "schemas/extraction_schema_v0.2.json",
                "scripts/*.py",
                "tests/*.py",
                "data/raw/ring1_platinum_queue.json",
                "data/raw/corpus_manifest_all_sources.jsonl",
                "data/raw/corpus_stratification_manifest.json",
                "data/raw/doc_citation_metadata_cache.json",
                "data/raw/layer2_genomic_and_jaspar_reference.json",
                "data/raw/layer5_dbsnp_clinvar_reference.json",
                "data/constructs/*.fasta",
                "data/extracted/*.json",
                "data/extracted/*.jsonl",
                "data/graveyard/idea_graveyard_v0.5.json",
                "data/omics/*.json",
                "data/protocols/*.json",
                "data/structures/*.cif",
                "data/structures/*.pdb",
                "data/structures/*.json",
                "data/structures/*.md",
                "data/structures/pymol/*",
            ],
            "local_research_cache_excluded_from_public_bundle": [
                "CLAUDE.md",
                "data/dossiers/*.json",
                "data/curated/*.json",
                "docs/OPENISLET-2026-2030-MASTER-TIMELINE.md",
                "docs/PROJECT-DECISIONS-LOG.md",
                "data/raw/alphagenome_avi_3118_snvs_full_archive.json",
                "data/raw/alphagenome_avi_26_seed_anchors.json",
                "data/raw/alphagenome_avi_live_receipt_2026_10_28.json",
                "data/raw/private_export_guard_patterns.json",
                "data/extracted/islet1_candidate_stack_rankings.json",
                "data/extracted/islet1_stack_ranking_v0.5.json",
            ],
        },
        "counts": {
            "corpus_manifest_rows": sum(
                1
                for _ in (RAW_DIR / "corpus_manifest_all_sources.jsonl").open(
                    encoding="utf-8"
                )
            ),
            "ring1_dossiers_total": len(studies),
            "dossiers_cache_total_including_auxiliary": len(dossier_files),
            "auxiliary_non_ring1_dossiers_total": len(dossier_files) - len(studies),
            "extracted_studies_total": len(studies),
            "extracted_arms_total": sum(len(s.get("arms", [])) for s in studies),
            "extracted_atomic_effects_total": n_effects,
            "extracted_thresholds_total": len(graph_data.get("regime_thresholds", [])),
            "extracted_failure_modes_total": len(
                graph_data.get("residual_failure_and_complementary_lever_pairs", [])
            ),
            "graveyard_entries_total": len(graveyard_data.get("graveyard_entries", [])),
            "omics_gnomad_raw_files_total": len(gnomad_files),
            "candidate_stacks_evaluated": ranking_payload.get(
                "total_valid_combinatorial_stacks_evaluated", 0
            ),
            "graph_target_nodes": len(graph_data.get("target_gene_atlas", [])),
            "layer2_construct_cassettes_total": layer2_atlas["summary"]["construct_cassettes_synthesized_and_audited"],
            "layer2_perturbation_arms_simulated": layer2_atlas["summary"]["perturbation_profiles_total"],
            "layer2_alphagenome_avi_snvs_archived": None,
            "layer2_avi_predictions_available": False,
            "layer2_unverified_private_archive_rows": layer2_atlas["summary"]["archival_table_supplied_rows_private_quarantine"],
            "layer3_pdb_experimental_structures_verified": layer3_atlas["summary"]["total_experimental_pdb_structures_verified"],
            "layer3_alphafold_structures_verified": layer3_atlas["summary"]["total_alphafold_v6_models_verified"],
            "layer3_immune_synapse_complexes_verified": len(layer3_atlas["complexes"]),
            "layer4_pde_regimes_simulated": n_layer4_pde_regimes,
            "layer4_paracrine_ode_grafts_simulated": len(layer4_atlas["part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode"]["graft_architectures_compared"]),
            "layer5_crispr_guides_compiled": layer5_atlas["summary"]["total_guides_compiled"],
            "layer5_virtual_patients_simulated": layer5_atlas["summary"]["total_virtual_patients_simulated"],
        },
        "artifacts": artifact_entries,
    }

    _require_generation_lease()
    write_documentation(ROOT, manifest_pin_count=len(artifact_entries), render=True,rollback_journal=_GENERATION_JOURNAL.get())
    for rel in artifact_entries:
        p = ROOT / rel
        artifact_entries[rel] = {"sha256": sha256_file(p), "bytes": p.stat().st_size}

    manifest_path = ROOT / "manifest.json"
    _require_generation_lease()
    _publish_generation_artifacts({manifest_path:(json.dumps(manifest,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode('utf-8')},ROOT/'.manifest-file-publication.lock')
    print(
        f"Wrote manifest.json and brief exports: {len(studies)} studies, {n_effects} atomic effects."
    )


# Explicit machine-verified open-license PMID allowlist; empty by default so all 563 studies enforce the <=240-char cap.
VERIFIED_OPEN_LICENSE_PMIDS: set[str] = set()

import re as _re

VALID_REDACTION_MARKER_RE = _re.compile(
    r"^\[REDACTED_NON_OA_EXCERPT: see PMID:\d{7,8}"
    r"(?: span=[A-Za-z0-9_.\[\]:|-]+)?"
    r"(?: field=[A-Za-z0-9_]+)?"
    r"(?: mod=[A-Za-z0-9_,-]+)?"
    r"(?: genes=[A-Za-z0-9_,-]+)?"
    r"(?: pert=(?:GOF_OE|KO_LOF|BOTH|NONE))?"
    r" quote_sha256=[0-9a-f]{64}\]$"
)


def _validate_public_redaction_token(text: str, *, expected_pmid: str | None = None, source_quote: str | None = None, expected_field: str | None = None) -> bool:
    import importlib.util
    helper_path = Path(__file__).parent / "public_excerpt_policy.py"
    # Loaded once per exporter module, including isolated fixture checkouts.
    global _REDACTION_HELPER
    if _REDACTION_HELPER is None:
        spec = importlib.util.spec_from_file_location("openislet_public_excerpt_policy", helper_path)
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        _REDACTION_HELPER = helper
    snapshot = _PUBLIC_POLICY_SNAPSHOT.get()
    helper = snapshot['helper'] if snapshot is not None else _REDACTION_HELPER
    return helper.validate_token(text, expected_pmid=expected_pmid, source_quote=source_quote, expected_field=expected_field)


_REDACTION_HELPER = None


def _public_redaction_field(name: str) -> str:
    _validate_public_redaction_token("")
    snapshot = _PUBLIC_POLICY_SNAPSHOT.get()
    return (snapshot['helper'] if snapshot is not None else _REDACTION_HELPER).safe_field(name)


def _is_authentic_redaction_marker(text: str) -> bool:
    """Compatibility name for bounded opaque syntax, not source acquisition proof."""
    return _validate_public_redaction_token(text)


_PUBLIC_SCHEMA_METADATA_CACHE = {}


def _is_bounded_public_metadata(value,field,*,bibliographic=False):
    _validate_public_redaction_token("")
    snapshot=_PUBLIC_POLICY_SNAPSHOT.get()
    helper=snapshot['helper'] if snapshot is not None else _REDACTION_HELPER
    return helper.is_bounded_metadata(value,field,bibliographic=bibliographic)


def _is_fixed_public_schema_metadata(value: str, field: str = "") -> bool:
    """Recognize finite schema codes only in their documented metadata fields.

    Free-text fields never gain an exemption merely by equalling a code from
    another property. Explicit aliases below are typed copies in conditioning
    contexts and graph projections, not arbitrary prose exemptions.
    """
    if not isinstance(value, str):
        return False
    # These fields carry bounded machine identity, not source prose. Arm IDs
    # are the deterministic identifiers minted by this exporter itself.
    identifiers = {
        'arm_id': r'pmid[1-9][0-9]{0,8}_arm[1-9][0-9]{0,5}',
        'focal_arm_id': r'pmid[1-9][0-9]{0,8}_arm[1-9][0-9]{0,5}',
        'comparator_arm_id': r'pmid[1-9][0-9]{0,8}_arm[1-9][0-9]{0,5}',
        'dossier_sha256': r'[0-9a-f]{64}',
        'quote_source_dossier_pmid': r'[1-9][0-9]{0,8}',
        'superseded_by_pmid': r'[1-9][0-9]{0,8}',
        'reconciled_from_superseded_preprint_pmid': r'[1-9][0-9]{0,8}',
    }
    if field in identifiers and _re.fullmatch(identifiers[field], value):
        return True
    snapshot = _PUBLIC_POLICY_SNAPSHOT.get()
    schema_bytes = snapshot['schema'] if snapshot is not None else _regular_file_bytes(ROOT, 'schemas/extraction_schema_v0.2.json')
    cache_key = hashlib.sha256(schema_bytes).digest()
    if cache_key not in _PUBLIC_SCHEMA_METADATA_CACHE:
        schema = json.loads(schema_bytes)
        by_field: dict[str, set[str]] = {}
        def collect(node, current_field=""):
            if isinstance(node, dict):
                by_field.setdefault(current_field, set()).update(v for v in node.get("enum", []) if isinstance(v, str))
                for key, child in node.items():
                    if key == "properties":
                        for name, prop in child.items(): collect(prop, name)
                    else: collect(child, current_field)
            elif isinstance(node, list):
                for child in node: collect(child, current_field)
        collect(schema)
        aliases = {
            "Z1_developmental_stage_and_timing": ("canonical_differentiation_stage_evaluated",),
            "Z2_genetic_background_and_allelic_dose": ("cell_source_class",),
            "Z3_immune_effector_and_stressor_pathway": ("outcome_domain", "challenge_type"),
            "Z6_microenvironment_oxygen_and_vascular_niche": ("transplant_site",),
            "Z7_promoter_and_epigenetic_context": ("promoter_or_locus_context",),
            "Z8_species_age_and_sex_context": ("max_biological_fidelity_tier",),
            "endpoint_name": ("outcome_domain",), "outcome_Y": ("outcome_domain",),
            "primary_gate_perturbed": ("functional_gate",),
            "secondary_gate_affected": ("functional_gate",),
            "functional_gates": ("functional_gate",),
            "hop1_relation": ("causal_relation",), "hop2_relation": ("causal_relation",),
            "study_arm_cell_source_classes": ("cell_source_class",),
            "study_arm_differentiation_stages": ("canonical_differentiation_stage_evaluated",),
        }
        for alias, origins in aliases.items():
            for origin in origins: by_field.setdefault(alias, set()).update(by_field.get(origin, set()))
        for name in ("systemic_immunosuppression_regimen", "encapsulation_or_biomaterial", "co_transplanted_cells_or_oxygen_source"):
            by_field.setdefault(name, set()).add("None")
        for name in ("measurement_unit", "units"):
            by_field.setdefault(name, set()).add("quarantined_ambiguous_numeric")
        span_pattern = schema["properties"]["atomic_effects"]["items"]["properties"]["source_span"]["pattern"]
        if len(_PUBLIC_SCHEMA_METADATA_CACHE) >= 16:
            _PUBLIC_SCHEMA_METADATA_CACHE.pop(next(iter(_PUBLIC_SCHEMA_METADATA_CACHE)))
        _PUBLIC_SCHEMA_METADATA_CACHE[cache_key] = ({k: frozenset(v) for k, v in by_field.items()}, _re.compile(span_pattern))
    by_field, span_pattern = _PUBLIC_SCHEMA_METADATA_CACHE[cache_key]
    if field == "gene_symbol_hgnc":
        _validate_public_redaction_token("")
        helper=snapshot['helper'] if snapshot is not None else _REDACTION_HELPER
        if value in helper.GENES or value in {"HLA-B_HLA-C_HLA-A24", "IL10_TGFB1_IL2_N88D", "JAK1_JAK2"}:
            return True
    return value in by_field.get(field, ()) or (field == "source_span" and span_pattern.fullmatch(value) is not None)


def _source_free_effect_identifier(value, field):
    # Source-free replay retains this bounded legacy identifier; it cannot
    # authenticate its original prose. With a private source available it is
    # deliberately counted for overlap, just like any source-derived slug.
    return field == 'effect_id' and isinstance(value, str) and _re.fullmatch(r'[a-z][a-z0-9_]{5,80}', value) is not None


def _rewrite_exact_arm_references(obj: object, identifier_map: dict[str, str]) -> object:
    """Publish opaque arm IDs, preserving every exact foreign-key reference."""
    if isinstance(obj, dict):
        return {k: _rewrite_exact_arm_references(v, identifier_map) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_rewrite_exact_arm_references(v, identifier_map) for v in obj]
    if isinstance(obj, str):
        return identifier_map.get(obj, obj)
    return obj


def _validate_imported_redaction_tokens(obj: object, pmid: str, dossier: dict | None, parent: dict | None = None, parent_key: str = "") -> None:
    """Input-origin references must not acquire trust merely by matching a regex."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            _validate_imported_redaction_tokens(value, pmid, dossier, obj, key)
    elif isinstance(obj, list):
        for value in obj:
            _validate_imported_redaction_tokens(value, pmid, dossier, parent, parent_key)
    elif isinstance(obj, str) and "[REDACTED_NON_OA_EXCERPT:" in obj:
        source_quote = None
        if dossier is not None:
            span = str((parent or {}).get("source_span") or "")
            m = _re.fullmatch(r"(.+)\[(\d+):(\d+)\]\|norm=exact", span)
            if m:
                value: object = dossier
                try:
                    tokens = _re.findall(r"([A-Za-z_][A-Za-z0-9_]*)|\[([0-9]+)\]", m[1])
                    rebuilt = ""
                    for key, index in tokens:
                        rebuilt += ("." if rebuilt and key else "") + key if key else f"[{index}]"
                        value = value[key] if key else value[int(index)]
                    if rebuilt != m[1]:
                        raise KeyError("Unsupported source path syntax")
                    if isinstance(value, str) and 0 <= int(m[2]) < int(m[3]) <= len(value):
                        source_quote = value[int(m[2]):int(m[3])]
                except (KeyError, TypeError, IndexError):
                    pass
            if source_quote is None:
                raise ValueError(f"Opaque input token for PMID {pmid} cannot be bound to the available source")
        if not _validate_public_redaction_token(obj, expected_pmid=pmid, source_quote=source_quote, expected_field=parent_key):
            raise ValueError(f"Invalid or unauthenticated public redaction token for PMID {pmid}")


def _has_symlink_component(path: Path, stop_at: Path | None = None) -> bool:
    """Return True if `path` or any parent up to `stop_at` is a symbolic link."""
    cur = path
    stop_res = stop_at.resolve() if stop_at is not None and not stop_at.is_symlink() else stop_at
    while True:
        if cur.is_symlink():
            return True
        if stop_res is not None and cur == stop_res:
            break
        if cur.parent == cur:
            break
        cur = cur.parent
    return False


def _source_words(text: str) -> list[str]:
    import unicodedata
    normalized = unicodedata.normalize('NFKD', text).casefold()
    normalized = ''.join(c for c in normalized if unicodedata.category(c) not in {'Cf','Mn','Me'})
    return _re.findall(r'[^\W_]+', normalized)


@lru_cache(maxsize=16)
def _normalized_source(text: str) -> str:
    # Complete immutable source content is the cache key, never a path/mtime.
    return " " + " ".join(_source_words(text)) + " "


def _max_verbatim_shingle_chars(text: object, source_text: str, shingle_words: int = 6, parent_key: str = "", title_str: str = "", bibliographic: bool = False) -> int:
    """Count source-overlapping prose, excluding only bounded metadata in its typed field.

    Exact short fragments count too. This is the same metadata scope used by
    enforcement and its post-condition; free-text enum lookalikes still count.
    """
    if not source_text or text is None:
        return 0
    if isinstance(text, dict):
        t_str = title_str or str((text.get("study_metadata") or {}).get("title") or "")
        return sum(_max_verbatim_shingle_chars(v, source_text, shingle_words, k, t_str, parent_key=='study_metadata') for k, v in text.items())
    if isinstance(text, list):
        return sum(_max_verbatim_shingle_chars(v, source_text, shingle_words, parent_key, title_str,bibliographic) for v in text)
    if not isinstance(text, str):
        return 0
    if _is_bounded_public_metadata(text,parent_key,bibliographic=bibliographic):
        return 0
    if not text or _is_authentic_redaction_marker(text) or _is_fixed_public_schema_metadata(text, parent_key):
        return 0
    w_txt = _source_words(text)
    src_norm = _normalized_source(source_text)
    if len(w_txt) < shingle_words:
        # A whole short source fragment is still verbatim. Otherwise source
        # sentences split across <=5-word leaves evade the aggregate budget.
        short_fragment = " " + " ".join(w_txt) + " "
        return len(text) if w_txt and short_fragment in src_norm else 0
    for i in range(len(w_txt) - shingle_words + 1):
        shingle = " " + " ".join(w_txt[i : i + shingle_words]) + " "
        if shingle in src_norm:
            return len(text)
    return 0


def _public_discovery_record(record: dict) -> dict:
    """Discovery metadata ships identifiers and abstract digests, never abstracts.

    A future live harvest may restore inline source text in the private raw
    cache. Export must independently enforce this boundary on every run.
    """
    if not isinstance(record, dict):
        raise ValueError("Discovery records must be JSON objects")
    # A shallow allowlist lets an adapter hide full text in authors/keywords.
    # Preserve the harvest's flat bibliographic schema; unexpected shapes fail
    # before publication rather than being silently converted or truncated.
    scalar_limits={'title':2000,'journal':1000,'first_author':512,'senior_author':512,'source_db':128}
    list_limits={'authors':(10000,512),'mesh_terms':(256,256),'keywords':(256,256),
                 'publication_types':(256,128),'matched_pillars':(256,128)}
    identifiers={'pmid':r'[1-9][0-9]{0,8}','pmcid':r'PMC[1-9][0-9]{0,9}',
                 'doi':r'10\.[0-9]{4,9}/[^\s\x00-\x1f\x7f]+','abstract_sha256':r'[0-9a-f]{64}'}
    flags={'has_abstract','has_open_pmc_fulltext'}
    metadata_fields=set(scalar_limits)|set(list_limits)|set(identifiers)|flags|{'pub_year'}
    def bounded_text(value, limit):
        return isinstance(value,str) and len(value)<=limit and _re.search(r'[\x00-\x1f\x7f]',value) is None
    public={}
    for key,value in record.items():
        if key not in metadata_fields:
            continue
        if value is None:
            valid=True
        elif key in scalar_limits:
            valid=bounded_text(value,scalar_limits[key])
        elif key in list_limits:
            count,limit=list_limits[key]
            valid=isinstance(value,list) and len(value)<=count and all(bounded_text(item,limit) for item in value)
        elif key in identifiers:
            valid=isinstance(value,str) and len(value)<=256 and _re.fullmatch(identifiers[key],value) is not None
        elif key in flags:
            valid=type(value) is bool
        else:
            valid=type(value) is int and 1<=value<=9999
        if not valid:
            raise ValueError(f'Discovery metadata field {key} has invalid type, shape or bounds')
        public[key]=value
    if "abstract" in record:
        abstract = record['abstract']
        if abstract is not None and not isinstance(abstract, str):
            raise ValueError("Discovery abstract must be text or null")
        public["has_abstract"] = bool(abstract)
        if abstract:
            public["abstract_sha256"] = hashlib.sha256(abstract.encode("utf-8")).hexdigest()
    return public


def public_token_errors(obj, expected_pmid=None, parent_key=''):
    """Check every opaque reference, including graph hops, without private text.

    This is a syntax/identity/field check only. Source-free public verification
    cannot authenticate the absent private source or its scientific meaning.
    """
    errors=[]
    def visit(value,pmid,key,path):
        if isinstance(value,dict):
            label=value.get('pmid') or value.get('study')
            if label:
                match=_re.search(r'(?<![0-9])([0-9]{7,8})(?![0-9])',str(label))
                if match:pmid=match[1]
            for name,child in value.items():
                child_pmid=pmid
                hop='hop1_study' if name in {'hop1_from','hop1_to','hop1_evidence','origin_node','bridge_node'} else 'hop2_study' if name in {'hop2_to','hop2_evidence','terminal_node'} else None
                if hop and value.get(hop):
                    match=_re.search(r'(?<![0-9])([0-9]{7,8})(?![0-9])',str(value[hop]))
                    if match:child_pmid=match[1]
                visit(child,child_pmid,name,path+'.'+name)
        elif isinstance(value,list):
            for index,child in enumerate(value):visit(child,pmid,key,path+f'[{index}]')
        elif isinstance(value,str) and '[REDACTED_NON_OA_EXCERPT:' in value:
            if not _validate_public_redaction_token(value,expected_pmid=pmid,expected_field=key):errors.append(path)
    visit(obj,expected_pmid,parent_key,'$')
    return errors


def public_article_verbatim_totals(studies, graph, sources):
    """Recount canonical study/graph fields together against pinned sources.

    JSONL aliases repeat these same canonical rows and are checked for exact
    materialization separately; their serialization is not a new annotation.
    """
    totals={pmid:0 for pmid in sources}
    normalized_sources=None
    unbound_matches={}
    for study in studies:
        pmid=str(study.get('study_metadata',{}).get('pmid',''))
        if pmid in sources:
            source,_=sources[pmid]
            totals[pmid]+=_max_verbatim_shingle_chars(study,source)
    def visit(value,pmid='',key=''):
        nonlocal normalized_sources
        if isinstance(value,dict):
            label=value.get('pmid') or value.get('study')
            if label:
                match=_re.search(r'(?<![0-9])([0-9]{7,8})(?![0-9])',str(label))
                if match:pmid=match[1]
            for name,child in value.items():
                child_pmid=pmid
                hop='hop1_study' if name in {'hop1_from','hop1_to','hop1_evidence','origin_node','bridge_node'} else 'hop2_study' if name in {'hop2_to','hop2_evidence','terminal_node'} else None
                if hop and value.get(hop):
                    match=_re.search(r'(?<![0-9])([0-9]{7,8})(?![0-9])',str(value[hop]))
                    if match:child_pmid=match[1]
                visit(child,child_pmid,name)
        elif isinstance(value,list):
            for child in value:visit(child,pmid,key)
        elif isinstance(value,str):
            if pmid in sources:
                source,_=sources[pmid]
                totals[pmid]+=_max_verbatim_shingle_chars(value,source,parent_key=key)
            else:
                # An invalid/missing graph identifier cannot hide prose from
                # the independent recount of source-overlapping strings.
                # Normalize each immutable source once for this recount; the
                # small process-wide LRU would thrash across a full corpus.
                if not value or _is_bounded_public_metadata(value,key) or _is_authentic_redaction_marker(value) or _is_fixed_public_schema_metadata(value,key):return
                if (key,value) not in unbound_matches:
                    words=_source_words(value);width=min(6,len(words))
                    if normalized_sources is None:normalized_sources={p:' '+' '.join(_source_words(source))+' ' for p,(source,_) in sources.items()}
                    shingles=[' '+' '.join(words[i:i+width])+' ' for i in range(len(words)-width+1)] if width else []
                    unbound_matches[key,value]=[p for p,source in normalized_sources.items() if any(shingle in source for shingle in shingles)]
                for source_pmid in unbound_matches[key,value]:totals[source_pmid]+=len(value)
    visit(graph)
    return totals


def _verify_pinned_source(root: Path, relative_path: str, manifest: dict) -> None:
    """Never relabel changed source bytes as the pinned release during export."""
    if relative_path == "manifest.json":
        return  # A manifest cannot contain its own digest.
    if not _canonical_relative_path(relative_path):
        raise ValueError(f"Unsafe manifest source path: {relative_path}")
    pin = manifest.get("artifacts", {}).get(relative_path, {})
    expected = pin.get("sha256")
    data = _regular_file_bytes(root, relative_path)
    actual = hashlib.sha256(data).hexdigest()
    if expected != actual or type(pin.get("bytes")) is not int or pin["bytes"] != len(data):
        raise ValueError(f"Pinned source hash mismatch for {relative_path}; regenerate the canonical manifest before export")


def manifest_integrity_errors(root: Path, manifest: dict, public: bool = False) -> list[str]:
    """Validate live identities, byte counts and exact declared public membership."""
    errors=[]
    artifacts=manifest.get('artifacts')
    artifact_errors = _artifact_inventory_errors(artifacts)
    if artifact_errors:
        return sorted(artifact_errors)
    inventory=manifest.get('public_export_file_inventory')
    if public:
        if (not isinstance(inventory,list) or not inventory or
                any(not _canonical_relative_path(rel) for rel in inventory) or
                len(set(inventory))!=len(inventory) or
                any(rel not in artifacts or rel in ('manifest.json','.openislet_public_bundle') for rel in inventory)):
            return ['invalid-public-inventory']
        required=set(inventory)
        errors.extend('forbidden:'+rel for rel in required if _forbidden_public_path(rel))
    else:
        required=set(artifacts)
    for rel in sorted(required):
        try:
            if rel=='manifest.json': raise ValueError('self pin')
            _verify_pinned_source(root,rel,manifest)
        except (OSError,ValueError,TypeError,KeyError): errors.append(rel)
    # Every file in a public tree must be declared. A cache directory name is
    # never an exemption for an undeclared payload; verification can use -B and
    # disable pytest's cache provider to keep runtime caches outside the bundle.
    if public:
        if (root / '.git').is_symlink():
            errors.append('unsafe:.git')
        for path in root.rglob('*'):
            rel=path.relative_to(root).as_posix()
            if rel == '.git' or rel.startswith('.git/'):
                continue
            if path.is_symlink() or (not path.is_dir() and (not stat.S_ISREG(path.lstat().st_mode) or path.lstat().st_nlink != 1)):
                errors.append('unsafe:'+rel); continue
            if not path.is_dir() and rel not in required|{'manifest.json','.openislet_public_bundle'}:
                errors.append('unexpected:'+rel)
        try: _regular_file_bytes(root,'.openislet_public_bundle')
        except (OSError,ValueError): errors.append('invalid-public-marker')
    return sorted(set(errors))


def export_public_bundle(dest_dir: Path) -> list[str]:
    # Snapshot imports temporarily affect Python's process-global import path
    # and bytecode setting. Serialize same-process exporters as well as taking
    # the destination lease used by independent processes.
    with _EXPORT_PROCESS_LOCK:return _export_public_bundle_serialized(dest_dir)


def _export_public_bundle_serialized(dest_dir: Path) -> list[str]:
    # Serialize publishers of the same destination without modifying source data.
    destination = Path(dest_dir).absolute()
    allowed_system_links = {Path('/tmp'),Path('/var'),Path('/etc')}
    if any(part.is_symlink() and part not in allowed_system_links for part in [destination,*destination.parents]):
        raise ValueError('Destination and ancestors may not be symbolic links')
    resolved=destination.resolve(); root=ROOT.resolve()
    if resolved == root or resolved in root.parents or root in resolved.parents:
        raise ValueError('Public export must be outside the repository tree')
    # Approved OS aliases (/var -> /private/var on macOS) were checked before
    # resolution. The lease itself always receives the canonical real path.
    destination=resolved
    destination.parent.mkdir(parents=True, exist_ok=True)
    with _publication_lock(destination.parent / ("." + destination.name + ".publish-lock")) as lease:
        return _export_public_bundle_locked(destination,publication_lease=lease)


def _export_public_bundle_locked(dest_dir: Path, publication_lease=None) -> list[str]:
    """Materializes a clean, rights-aware public distribution directory containing strictly the allowlisted files
    and enforcing the per-article non-OA aggregate verbatim quote ceiling (<=240 chars and <25% of abstract length)."""
    import re
    import shutil
    import tempfile
    def require_lease():
        if publication_lease is not None:publication_lease.require_owned()

    # Check symlink status on dest_dir and ALL ancestor components BEFORE .resolve() so symlinked destinations are rejected
    raw_dest = Path(dest_dir)
    allowed_sys_symlinks = {
        Path("/tmp"), Path("/var"), Path("/etc"),
    }
    for anc in [raw_dest, *raw_dest.parents]:
        if anc.is_symlink() and anc not in allowed_sys_symlinks:
            raise ValueError(f"Destination or ancestor {anc} must not be a symlink.")

    root_resolved = ROOT.resolve()
    final_dest_dir = raw_dest.resolve()
    if (
        final_dest_dir == root_resolved
        or final_dest_dir in root_resolved.parents
        or root_resolved in final_dest_dir.parents
    ):
        raise ValueError(
            f"Refusing to export public bundle into or over repository tree: {final_dest_dir} (ROOT={root_resolved})"
        )
    if final_dest_dir.exists():
        if not final_dest_dir.is_dir() or final_dest_dir.is_symlink():
            raise ValueError(f"Destination {final_dest_dir} must be a regular directory, not a file or symlink.")
        entries = list(final_dest_dir.iterdir())
        if entries and not (final_dest_dir / ".openislet_public_bundle").exists():
            raise ValueError(
                f"Refusing to overwrite non-empty directory {final_dest_dir} lacking '.openislet_public_bundle' marker."
            )
        if entries:
            _regular_file_bytes(final_dest_dir,'.openislet_public_bundle')

    manifest_bytes = _regular_file_bytes(ROOT, "manifest.json")
    manifest = json.loads(manifest_bytes)
    artifact_errors = _artifact_inventory_errors(manifest.get('artifacts'))
    if artifact_errors:
        raise ValueError('Invalid release artifact inventory: ' + ', '.join(artifact_errors))
    pinned_artifacts = set(manifest.get("artifacts", {}).keys()) | {"manifest.json"}
    allow = manifest["public_redistribution_allowlist_spec"]
    max_non_oa_chars, max_non_oa_share = public_excerpt_limits(allow)
    included: set[str] = set()
    # A glob is not proof that every file in the release still exists. Derive
    # the required set from the pinned inventory before enumerating disk.
    patterns = allow["public_package_included_paths"]
    excluded = allow["local_research_cache_excluded_from_public_bundle"]
    required = {rel for rel in pinned_artifacts if any(fnmatch.fnmatchcase(rel, pat) for pat in patterns)
                and not any(fnmatch.fnmatchcase(rel, pat) for pat in excluded)}
    if any(_forbidden_public_path(rel) for rel in required):
        raise ValueError('Public policy may not include private caches, planning documents or frozen ranking artifacts')
    for rel in required:
        _regular_file_bytes(ROOT, rel)
    for pat in allow["public_package_included_paths"]:
        if Path(pat).is_absolute() or ".." in Path(pat).parts:
            raise ValueError(f"Unsafe public include pattern: {pat}")
        # Also check non-glob prefix directory components for symlinks before globbing
        pat_parts = Path(pat).parts
        cur_check = ROOT
        for part in pat_parts:
            if "*" in part or "?" in part or "[" in part:
                break
            cur_check = cur_check / part
            if cur_check.is_symlink():
                raise ValueError(f"Symlink directory not permitted in public bundle allowlist: {cur_check}")
        for p in ROOT.glob(pat):
            if _has_symlink_component(p, stop_at=ROOT):
                raise ValueError(f"Symlink or symlinked ancestor not permitted in public bundle allowlist: {p}")
            info = p.lstat()
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise ValueError(f"Public include matched a non-regular or multiply-linked file: {p}")
            if p.is_file():
                rel_str = str(p.relative_to(ROOT))
                if rel_str not in pinned_artifacts:
                    raise ValueError(
                        f"Unpinned file {rel_str!r} matched public bundle glob {pat!r} but is absent from manifest['artifacts']"
                    )
                included.add(rel_str)
    for pat in allow["local_research_cache_excluded_from_public_bundle"]:
        if not isinstance(pat,str) or Path(pat).is_absolute() or ".." in Path(pat).parts:
            raise ValueError("Unsafe public exclusion pattern")
        for p in ROOT.glob(pat):
            if p.is_file() or p.is_symlink():
                included.discard(str(p.relative_to(ROOT)))

    if included != required:
        raise ValueError("Public allowlist does not exactly match the pinned release inventory")
    # Private source text controls the quote policy, so its pins matter even
    # though those dossiers never ship. Recheck on every actual read below.
    if 'public_export_file_inventory' not in manifest:
        for rel in manifest.get('artifacts',{}):
            if rel.startswith('data/dossiers/'):
                _verify_pinned_source(ROOT,rel,manifest)
    dossier_inputs = sorted((ROOT / "data" / "dossiers").glob("*.json"))
    for source in dossier_inputs:
        _verify_pinned_source(ROOT, source.relative_to(ROOT).as_posix(), manifest)

    def read_dossier(source):
        rel = source.relative_to(ROOT).as_posix()
        data = _regular_file_bytes(ROOT, rel)
        if hashlib.sha256(data).hexdigest() != manifest["artifacts"][rel]["sha256"]:
            raise ValueError("Private source changed during public export")
        return json.loads(data)

    # Validate source pins before creating any staging output. Unique siblings
    # prevent deleting unrelated directories or colliding with another export.
    for rel in sorted(included):
        _verify_pinned_source(ROOT, rel, manifest)
    final_dest_dir.parent.mkdir(parents=True, exist_ok=True)
    require_lease()
    staging_dir = Path(tempfile.mkdtemp(prefix=f".{final_dest_dir.name}.staging-", dir=final_dest_dir.parent))

    policy_token = None
    import sys
    previous_bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    previous_sys_path = list(sys.path)
    try:
        (staging_dir / ".openislet_public_bundle").write_text(
            "openislet_public_bundle_v0.8.0\n", encoding="utf-8"
        )
        work_dir = staging_dir
        for rel in sorted(included):
            require_lease()
            src = ROOT / rel
            if _has_symlink_component(src, stop_at=ROOT):
                raise ValueError(f"Symlink not permitted in public bundle export: {src}")
            dst = work_dir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            _verify_pinned_source(work_dir, rel, manifest)

        import types
        public_helper = types.ModuleType('openislet_export_snapshot_policy')
        public_helper.__file__ = str(work_dir / 'scripts/public_excerpt_policy.py')
        exec(compile(_regular_file_bytes(work_dir, 'scripts/public_excerpt_policy.py'), public_helper.__file__, 'exec'), public_helper.__dict__)
        schema_bytes = _regular_file_bytes(work_dir, 'schemas/extraction_schema_v0.2.json')
        policy_token = _PUBLIC_POLICY_SNAPSHOT.set({'helper':public_helper, 'schema':schema_bytes})

        discovery_path = work_dir / "data/raw/corpus_manifest_all_sources.jsonl"
        if discovery_path.exists():
            discovery_rows = [json.loads(line) for line in discovery_path.read_text(encoding="utf-8").splitlines() if line.strip()]
            public_rows = [_public_discovery_record(row) for row in discovery_rows]
            if public_rows != discovery_rows:
                discovery_path.write_text("".join(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n" for row in public_rows), encoding="utf-8")

        # Collect known gene symbols and aliases so redacted non-OA quotes preserve exact gene attribution in rank_islet1_stacks.py
        import importlib.util as _ilu
        _rank_spec = _ilu.spec_from_file_location("openislet_rank_helper", work_dir / "scripts" / "rank_islet1_stacks.py")
        _rank_mod = _ilu.module_from_spec(_rank_spec)
        _rank_spec.loader.exec_module(_rank_mod)
        _strict_gene_match = _rank_mod.strict_gene_match

        def _expand_gene_symbol(raw_g: str) -> list[str]:
            if not raw_g:
                return []
            if raw_g == "HLA-B_HLA-C_HLA-A24":
                return ["HLA-B", "HLA-C", "HLA-A"]
            if raw_g == "IL10_TGFB1_IL2_N88D":
                return ["IL10", "TGFB1", "IL2"]
            if raw_g == "JAK1_JAK2":
                return ["JAK1", "JAK2"]
            if "_" in raw_g:
                return [tok for tok in raw_g.split("_") if tok and not re.match(r"^[A-Z]\d+[A-Z]$", tok)]
            return [raw_g]

        omics_raw = json.loads(_regular_file_bytes(work_dir, 'data/omics/target_prior_table.json'))
        omics_targets = omics_raw.get("targets", []) if isinstance(omics_raw, dict) else []
        omics_genes = {t["gene_symbol"] for t in omics_targets if isinstance(t, dict) and t.get("gene_symbol")} | set((omics_raw.get("target_priors") or {}).keys())
        alias_map = {
            "ZN-T8": "SLC30A8", "ZNT8": "SLC30A8", "IGRP": "G6PC2", "GAD65": "GAD2",
            "IA-2": "PTPRN", "PD-L1": "CD274", "PDL1": "CD274", "A20": "TNFAIP3",
            "PI-9": "SERPINB9", "MDA5": "IFIH1", "CAR": "CXADR", "ICASP9": "CASP9", "IND-CASP9": "CASP9",
            "IL-10": "IL10",
        }
        known_genes = omics_genes | set(alias_map.values())
        study_paths = sorted(
            [p for p in (work_dir / "data" / "extracted").glob("*.json") if p.stem.isdigit()],
            key=lambda p: p.stem,
        )
        for sp in study_paths:
            doc = json.loads(sp.read_text(encoding="utf-8"))
            for arm in doc.get("arms", []):
                for ed in arm.get("gate2_3_4_genetic_and_epitope_edits", []):
                    for g_atom in _expand_gene_symbol(str(ed.get("gene_symbol_hgnc") or "")):
                        known_genes.add(g_atom)

        _validate_public_redaction_token("")
        known_genes &= public_helper.GENES

        # Identifier slugs must not become an encoded source-prose channel.
        # The private IDs remain untouched; publish a deterministic bijection.
        arm_identifier_map = {}
        for study_path in study_paths:
            source_study = json.loads(study_path.read_text(encoding="utf-8"))
            source_pmid = str(source_study["study_metadata"]["pmid"])
            for index, arm in enumerate(source_study.get("arms", []), 1):
                old_id, new_id = arm["arm_id"], f"pmid{source_pmid}_arm{index}"
                if old_id in arm_identifier_map and arm_identifier_map[old_id] != new_id:
                    raise ValueError("Non-unique study arm identifier cannot be safely exported")
                arm_identifier_map[old_id] = new_id

        marker_readout_genes = {"INS", "GCG", "SST", "CHGA", "NKX6-1", "PDX1", "MAFA"}
        gof_tokens = (
            "OVEREXPRESSION", "OVEREXPRESS", "KNOCK-IN", "KNOCKIN", "TRANSGENE",
            "TRANSGENIC", "ECTOPIC", "LENTIVIRAL", "GAIN-OF-FUNCTION", "GAIN OF FUNCTION",
            "FORCED EXPRESSION", "CONSTITUTIVE EXPRESSION", "INDUCIBLE EXPRESSION",
        )
        ko_tokens = (
            "KNOCKOUT", " KO ", " KO,", " KO.", "-/-", "−/−", "DEFICIENT", "DEFICIENCY",
            "DELETION", "NULL", "ABLATION", "SILENCING", "SHRNA", "SIRNA", "LOSS-OF-FUNCTION",
            "LOSS OF FUNCTION", "INHIBITION", "INHIBITOR", "BLOCKADE",
        )

        redacted_studies = []
        non_oa_pmids: set[str] = set()
        seen_pmids: set[str] = set()
        dossier_source_cache: dict[str, tuple[str, int]] = {}
        # Shared per-PMID verbatim shingle counter across BOTH study JSONs and cross_study_discovery_graph.json
        pmid_shingle_used: dict[str, int] = {}
        for sp in study_paths:
            doc = _rewrite_exact_arm_references(json.loads(sp.read_text(encoding="utf-8")), arm_identifier_map)
            meta = doc.get("study_metadata", {})
            pmid = str(meta.get("pmid", sp.stem))
            if pmid != sp.stem or pmid in seen_pmids or re.fullmatch(r'[1-9][0-9]{0,8}',pmid) is None:
                raise ValueError('Study filename, PMID and unique article identity must agree')
            seen_pmids.add(pmid)
            source_path = (ROOT / "data" / "dossiers") / f"{pmid}.json"
            imported_source = read_dossier(source_path) if source_path.exists() else None
            _validate_imported_redaction_tokens(doc, pmid, imported_source)
            has_verified_open_license = (
                pmid in VERIFIED_OPEN_LICENSE_PMIDS
                and meta.get("license") in ("CC-BY", "CC-BY-4.0", "CC0", "Public_Domain")
            )
            if not has_verified_open_license:
                non_oa_pmids.add(pmid)
                dos_path = (ROOT / "data" / "dossiers") / f"{pmid}.json"
                full_src_text = ""
                if dos_path.exists():
                    dos_obj = read_dossier(dos_path)
                    ab_str = dos_obj.get("abstract") or ""
                    ab_len = max(1, len(ab_str))
                    ft_secs = dos_obj.get("fulltext_data_sections") or {}
                    full_src_text = " ".join(
                        [ str(dos_obj.get("title") or ""), ab_str ]
                        + [ str(v) for v in ft_secs.values() if v ]
                    )
                    char_budget = min(max_non_oa_chars, max(0, int((max_non_oa_share - 0.01) * ab_len)))
                else:
                    # Fail closed (budget = 0) when a non-OA study has no local dossier
                    char_budget = 0
                dossier_source_cache[pmid] = (full_src_text, char_budget)
                pmid_shingle_used[pmid] = 0
                used_chars = 0
                study_genes = {
                    g_atom
                    for arm in doc.get("arms", [])
                    for ed in arm.get("gate2_3_4_genetic_and_epitope_edits", [])
                    for g_atom in _expand_gene_symbol(str(ed.get("gene_symbol_hgnc") or ""))
                }
                eff_redaction_meta: dict[str, tuple[str, str]] = {}
                eff_focal_redaction_meta: dict[str, tuple[str, str]] = {}
                for eff in doc.get("atomic_effects", []):
                    focal_only = eff.get("focal_factor_X") or ""
                    focal_only_u = focal_only.upper()
                    q_orig = eff.get("verbatim_evidence_quote") or ""
                    fx_orig = focal_only + " " + q_orig
                    fx_orig_u = fx_orig.upper()
                    matched_g = set()
                    for g in known_genes:
                        txt_g = focal_only if (g in marker_readout_genes and g not in study_genes) else fx_orig
                        txt_g_u = focal_only_u if (g in marker_readout_genes and g not in study_genes) else fx_orig_u
                        if g.upper() in txt_g_u and _strict_gene_match(g, txt_g):
                            matched_g.add(g)
                    for al, can in alias_map.items():
                        txt_a = focal_only if (can in marker_readout_genes and can not in study_genes) else fx_orig
                        txt_a_u = focal_only_u if (can in marker_readout_genes and can not in study_genes) else fx_orig_u
                        if al in txt_a_u and _strict_gene_match(al, txt_a):
                            matched_g.add(can)
                    if not matched_g and len(study_genes) == 1:
                        matched_g = set(study_genes)
                    has_gof = any(tok in focal_only_u or tok in fx_orig_u for tok in gof_tokens)
                    has_ko = any(tok in focal_only_u or tok in fx_orig_u for tok in ko_tokens)
                    pert_tag = "BOTH" if (has_gof and has_ko) else ("GOF_OE" if has_gof else ("KO_LOF" if has_ko else "NONE"))
                    genes_val = ",".join(sorted(matched_g)) if matched_g else "NONE"
                    effect_key = str(eff.get("effect_id") or "")
                    eff_redaction_meta[effect_key] = (genes_val, pert_tag)
                    # Focal metadata must not acquire additional readout genes
                    # from the evidence quote when public wording is redacted.
                    focal_genes = {g for g in known_genes if _strict_gene_match(g, focal_only)}
                    focal_genes.update(can for alias, can in alias_map.items()
                                       if can in known_genes and _strict_gene_match(alias, focal_only))
                    focal_genes_val = ",".join(sorted(focal_genes)) if focal_genes else "NONE"
                    eff_focal_redaction_meta[effect_key] = (focal_genes_val, pert_tag)

                    # Replace verbatim sentence slice in outcome_variable_Y with structured endpoint descriptor
                    out_y = eff.get("outcome_variable_Y") or ""
                    pol_reg = eff.get("endpoint_polarity_registry") or {}
                    pol_cls = pol_reg.get("polarity_class", "higher_is_better")
                    orient_sign = pol_reg.get("orientation_sign", 1)
                    if out_y.startswith("[Endpoint:"):
                        # Curated structured endpoint summary: verify if it shares any 6-word verbatim shingle with source dossier
                        v_len = _max_verbatim_shingle_chars(out_y, full_src_text) if full_src_text else len(out_y)
                        if v_len > 0:
                            if used_chars + len(out_y) <= char_budget:
                                used_chars += len(out_y)
                            else:
                                m_ep = re.match(r"^(\[Endpoint:[^\]]+\])", out_y)
                                ep_prefix = m_ep.group(1) if m_ep else f"[Endpoint: {eff.get('outcome_domain', 'General')}]"
                                err_note = " (verified against 2026 author/affiliation erratum PMID:41887653)" if "41887653" in out_y else ""
                                eff["outcome_variable_Y"] = f"{ep_prefix}{err_note} [polarity={pol_cls} orientation_sign={orient_sign}]"
                    else:
                        eff["outcome_variable_Y"] = (
                            f"[Endpoint: {eff.get('outcome_domain', 'General')} | polarity={pol_cls} | orientation_sign={orient_sign} | unit={eff.get('measurement_unit') or 'qualitative'}]"
                        )

                    comp_x0 = eff.get("comparator_X0") or ""
                    if (not full_src_text and len(comp_x0) > 0) or _max_verbatim_shingle_chars(comp_x0, full_src_text) > 0:
                        if used_chars + len(comp_x0) <= char_budget:
                            used_chars += len(comp_x0)
                        else:
                            eff["comparator_X0"] = "Control / comparator condition (see structured arm definitions)"

                    q = q_orig
                    if _is_authentic_redaction_marker(q):
                        # A validated opaque pointer already contains no quote;
                        # retain its source digest rather than hashing the token.
                        pass
                    elif used_chars + len(q) <= char_budget and len(q) > 0:
                        used_chars += len(q)
                    else:
                        q_sha = hashlib.sha256(q.encode("utf-8")).hexdigest()
                        span_tok = re.sub(r"[^A-Za-z0-9_.\[\]:|-]", "_", str(eff.get("source_span") or "none"))
                        eff["verbatim_evidence_quote"] = (
                            f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field=verbatim_evidence_quote genes={genes_val} pert={pert_tag} quote_sha256={q_sha}]"
                        )
                # Enforce article-wide non-OA verbatim cap across all other source-derived text channels
                for arm in doc.get("arms", []):
                    for ed in arm.get("gate2_3_4_genetic_and_epitope_edits", []):
                        smc = ed.get("sub_molecular_coordinate") or {}
                        desc = smc.get("exact_variant_or_swap_description") or ""
                        if desc and not desc.startswith("Experimental perturbation of ") and not _is_authentic_redaction_marker(desc):
                            if used_chars + len(desc) <= char_budget:
                                used_chars += len(desc)
                            else:
                                d_sha = hashlib.sha256(desc.encode("utf-8")).hexdigest()
                                raw_mod = str(ed.get("edit_modality") or "none")
                                mod_tok = raw_mod if raw_mod in public_helper.MODALITIES else "none"
                                g_tok = ",".join(sorted(set(_expand_gene_symbol(str(ed.get("gene_symbol_hgnc") or ""))) & known_genes)) or "NONE"
                                smc["exact_variant_or_swap_description"] = (
                                    f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field=exact_variant_or_swap_description mod={mod_tok} genes={g_tok} quote_sha256={d_sha}]"
                                )
                for em in doc.get("effect_modifier_interactions", []):
                    for mfield in ("level_1_condition_and_effect", "level_2_condition_and_effect", "mechanistic_explanation"):
                        mval = em.get(mfield) or ""
                        if mval and not _is_authentic_redaction_marker(mval):
                            if used_chars + len(mval) <= char_budget:
                                used_chars += len(mval)
                            else:
                                m_sha = hashlib.sha256(mval.encode("utf-8")).hexdigest()
                                fig_tok = re.sub(r"[^A-Za-z0-9_.\[\]:|-]", "_", str(em.get("source_figure_or_table") or "none"))
                                em[mfield] = (
                                    f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field={mfield} quote_sha256={m_sha}]"
                                )
                for edge in doc.get("causal_mechanism_edges", []):
                    for nfield in ("upstream_node", "downstream_node"):
                        nval = edge.get(nfield) or ""
                        if "[" in nval and nval.endswith("]"):
                            nval = nval.split("[", 1)[0].strip()
                        if "(" in nval and nval.endswith(")") and (not full_src_text or _max_verbatim_shingle_chars(nval, full_src_text) > 0):
                            nval = nval.split("(", 1)[0].strip()
                        edge[nfield] = nval
                    esrc = edge.get("evidence_source") or ""
                    if esrc and not _is_authentic_redaction_marker(esrc):
                        if used_chars + len(esrc) <= char_budget:
                            used_chars += len(esrc)
                        else:
                            e_sha = hashlib.sha256(esrc.encode("utf-8")).hexdigest()
                            edge["evidence_source"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field=evidence_source quote_sha256={e_sha}]"
                bridges = doc.get("cross_study_discovery_bridges") or {}
                for rv in bridges.get("residual_failure_and_escape_vectors", []):
                    for rfield in ("residual_loss_or_limitation", "escape_mechanism_identified"):
                        rval = rv.get(rfield) or ""
                        if rval and not _is_authentic_redaction_marker(rval):
                            if used_chars + len(rval) <= char_budget:
                                used_chars += len(rval)
                            else:
                                r_sha = hashlib.sha256(rval.encode("utf-8")).hexdigest()
                                rv[rfield] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field={rfield} quote_sha256={r_sha}]"
                for pf in bridges.get("collateral_and_pleiotropic_findings", []):
                    pval = pf.get("observation") or ""
                    if pval and not _is_authentic_redaction_marker(pval):
                        if used_chars + len(pval) <= char_budget:
                            used_chars += len(pval)
                        else:
                            p_sha = hashlib.sha256(pval.encode("utf-8")).hexdigest()
                            pf["observation"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field=observation quote_sha256={p_sha}]"
                for rt in bridges.get("regime_thresholds", []):
                    tval = rt.get("verbatim_quote") or ""
                    if tval and not _is_authentic_redaction_marker(tval):
                        if used_chars + len(tval) <= char_budget:
                            used_chars += len(tval)
                        else:
                            t_sha = hashlib.sha256(tval.encode("utf-8")).hexdigest()
                            rt["verbatim_quote"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field=verbatim_quote quote_sha256={t_sha}]"
                title_str = str((doc.get("study_metadata") or {}).get("title") or "")

                def _enforce_doc_shingle_budget(obj: object, parent_key: str = "", cur_eff_id: str = "", bibliographic: bool = False) -> object:
                    if isinstance(obj, dict):
                        eff_id_here = str(obj.get("effect_id") or cur_eff_id)
                        for k in list(obj.keys()):
                            obj[k] = _enforce_doc_shingle_budget(obj[k], k, eff_id_here,parent_key=='study_metadata')
                        return obj
                    if isinstance(obj, list):
                        for idx in range(len(obj)):
                            obj[idx] = _enforce_doc_shingle_budget(obj[idx], parent_key, cur_eff_id,bibliographic)
                        return obj
                    if isinstance(obj, str):
                        if _is_bounded_public_metadata(obj,parent_key,bibliographic=bibliographic):
                            return obj
                        if not full_src_text and _source_free_effect_identifier(obj, parent_key):
                            return obj
                        if _is_authentic_redaction_marker(obj) or _is_fixed_public_schema_metadata(obj, parent_key):
                            return obj
                        f_tok = _public_redaction_field(parent_key or "field")
                        gp_suffix = ""
                        if parent_key in ("focal_factor_X", "verbatim_evidence_quote") and cur_eff_id in eff_redaction_meta:
                            metadata = eff_focal_redaction_meta if parent_key == "focal_factor_X" else eff_redaction_meta
                            g_v, p_v = metadata[cur_eff_id]
                            gp_suffix = f" genes={g_v} pert={p_v}"
                        # Spoofed redaction prefix or missing dossier fails closed immediately
                        if obj.startswith("[REDACTED_NON_OA_EXCERPT:"):
                            s_sha = hashlib.sha256(obj.encode("utf-8")).hexdigest()
                            return f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field={f_tok}{gp_suffix} quote_sha256={s_sha}]"
                        v_len = _max_verbatim_shingle_chars(obj, full_src_text) if full_src_text else len(obj)
                        if v_len > 0:
                            if pmid_shingle_used[pmid] + v_len <= char_budget:
                                pmid_shingle_used[pmid] += v_len
                            else:
                                s_sha = hashlib.sha256(obj.encode("utf-8")).hexdigest()
                                return f"[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field={f_tok}{gp_suffix} quote_sha256={s_sha}]"
                        return obj
                    return obj

                _enforce_doc_shingle_budget(doc)
                sp.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            redacted_studies.append(doc)

        # Synchronize exported JSONL aliases and cross_study_discovery_graph.json in work_dir
        dossiers_jsonl = work_dir / "data" / "extracted" / "ring1_study_dossiers_v0.5.jsonl"
        with dossiers_jsonl.open("w", encoding="utf-8") as fh:
            for st in redacted_studies:
                fh.write(json.dumps(st, ensure_ascii=False) + "\n")

        effects_jsonl = work_dir / "data" / "extracted" / "ring1_atomic_effects_v0.5.jsonl"
        with effects_jsonl.open("w", encoding="utf-8") as fh:
            for st in redacted_studies:
                meta = st.get("study_metadata", {})
                pmid = str(meta.get("pmid", ""))
                cit = meta.get("citation_short", "")
                for eff in st.get("atomic_effects", []):
                    fh.write(json.dumps({"pmid": pmid, "citation_short": cit, "source_text_type": meta.get("source_text_type"), "superseded_by_pmid": meta.get("superseded_by_pmid"), **eff}, ensure_ascii=False) + "\n")

        graph_path = work_dir / "data" / "extracted" / "cross_study_discovery_graph.json"
        if graph_path.exists():
            g_obj = _rewrite_exact_arm_references(json.loads(graph_path.read_text(encoding="utf-8")), arm_identifier_map)
            def _extract_pmid(label_or_id: str) -> str:
                m_p = re.search(r"(?:PMID:)?(\d{7,8})", str(label_or_id or ""))
                return m_p.group(1) if m_p else ""

            def _strip_node_quote_suffix(node_str: str, pm: str) -> str:
                s = node_str
                if "[" in s and s.endswith("]"):
                    s = s.split("[", 1)[0].strip()
                if "(" in s and s.endswith(")"):
                    full_src_t, _ = dossier_source_cache.get(pm, ("", 0))
                    if not full_src_t or _max_verbatim_shingle_chars(s, full_src_t) > 0:
                        s = s.split("(", 1)[0].strip()
                return s

            for item in g_obj.get("regime_thresholds", []):
                pm = _extract_pmid(item.get("pmid") or item.get("study", ""))
                if pm in non_oa_pmids and item.get("verbatim_quote") and not _is_authentic_redaction_marker(item["verbatim_quote"]):
                    vq = item["verbatim_quote"]
                    item["verbatim_quote"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field=verbatim_quote quote_sha256={hashlib.sha256(vq.encode('utf-8')).hexdigest()}]"
            for item in g_obj.get("residual_failure_and_complementary_lever_pairs", []):
                pm = _extract_pmid(item.get("pmid") or item.get("study", ""))
                if pm in non_oa_pmids:
                    for rf_name in ("verbatim_quote", "residual_loss_or_limitation", "escape_mechanism_identified"):
                        vq = item.get(rf_name) or ""
                        if vq and not _is_authentic_redaction_marker(vq):
                            item[rf_name] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field={rf_name} quote_sha256={hashlib.sha256(vq.encode('utf-8')).hexdigest()}]"
            for mclass, mod_list in (g_obj.get("effect_modifiers_by_canonical_class") or g_obj.get("effect_modifiers_by_class_Z1_to_Z8") or {}).items():
                for item in mod_list:
                    pm = _extract_pmid(item.get("study", ""))
                    if pm in non_oa_pmids:
                        for mf_name in ("level_1", "level_2", "mechanistic_explanation"):
                            mv = item.get(mf_name) or ""
                            if mv and not _is_authentic_redaction_marker(mv):
                                item[mf_name] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field={mf_name} quote_sha256={hashlib.sha256(mv.encode('utf-8')).hexdigest()}]"
            for item in g_obj.get("cross_gate_pleiotropic_discoveries", []):
                pm = _extract_pmid(item.get("study", ""))
                if pm in non_oa_pmids:
                    obs = item.get("observation") or ""
                    if obs and not _is_authentic_redaction_marker(obs):
                        item["observation"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field=observation quote_sha256={hashlib.sha256(obs.encode('utf-8')).hexdigest()}]"
            bridge_nodes_obj = g_obj.get("shared_cross_study_bridge_nodes", {})
            for node in (bridge_nodes_obj.values() if isinstance(bridge_nodes_obj, dict) else bridge_nodes_obj):
                for lm in node.get("linked_mechanisms_and_bottlenecks", []):
                    pm = _extract_pmid(lm.get("study", ""))
                    if pm in non_oa_pmids:
                        stmt = lm.get("statement") or ""
                        if " -> " in stmt:
                            parts = stmt.split(" -> ", 1)
                            lhs = _strip_node_quote_suffix(parts[0].strip(), pm)
                            rhs = _strip_node_quote_suffix(parts[1].strip(), pm)
                            lm["statement"] = f"{lhs} -> {rhs}"
                        elif stmt and not _is_authentic_redaction_marker(stmt):
                            lm["statement"] = (
                                f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field=statement quote_sha256={hashlib.sha256(stmt.encode('utf-8')).hexdigest()}]"
                            )
            for ch in g_obj.get("multi_hop_cross_study_causal_chains", []):
                for h_idx in (1, 2):
                    pm = _extract_pmid(ch.get(f"hop{h_idx}_study", ""))
                    if pm in non_oa_pmids:
                        for nf in ("hop1_from", "hop1_to", "hop2_to"):
                            if nf in ch and isinstance(ch[nf], str):
                                ch[nf] = _strip_node_quote_suffix(ch[nf], pm)
                        ev = ch.get(f"hop{h_idx}_evidence") or ""
                        if ev and not _is_authentic_redaction_marker(ev):
                            ch[f"hop{h_idx}_evidence"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field=hop{h_idx}_evidence quote_sha256={hashlib.sha256(ev.encode('utf-8')).hexdigest()}]"
            for edge in g_obj.get("all_causal_edges") or g_obj.get("causal_edges_all") or []:
                pm = _extract_pmid(edge.get("study", ""))
                if pm in non_oa_pmids:
                    for nfield in ("upstream_node", "downstream_node"):
                        nval = edge.get(nfield) or ""
                        edge[nfield] = _strip_node_quote_suffix(nval, pm)
                    esrc = edge.get("evidence_source") or ""
                    if esrc and not _is_authentic_redaction_marker(esrc):
                        edge["evidence_source"] = f"[REDACTED_NON_OA_EXCERPT: see PMID:{pm} field=evidence_source quote_sha256={hashlib.sha256(esrc.encode('utf-8')).hexdigest()}]"

            # Recursive graph-wide shingle budget enforcement sharing `pmid_shingle_used[pm]` with study JSONs:
            normalized_sources={pm:_normalized_source(source) for pm,(source,_) in dossier_source_cache.items() if source}
            unbound_match_cache={}
            def _unbound_source_match(value,key):
                if _is_bounded_public_metadata(value,key):
                    return ''
                if _is_authentic_redaction_marker(value) or _is_fixed_public_schema_metadata(value,key):return ''
                words=_source_words(value)
                if not words:return ''
                cache_key=tuple(words)
                if cache_key not in unbound_match_cache:
                    width=min(6,len(words))
                    shingles=[' '+' '.join(words[i:i+width])+' ' for i in range(len(words)-width+1)]
                    unbound_match_cache[cache_key]=next((pm for pm,source in normalized_sources.items() if any(shingle in source for shingle in shingles)),'')
                return unbound_match_cache[cache_key]
            def _enforce_graph_shingle_budget(obj: object, current_pm: str = "", parent_key: str = "") -> object:
                if isinstance(obj, dict):
                    for key in obj:
                        if len(_source_words(key))>=6 and _unbound_source_match(key,''):
                            raise ValueError('Graph field name contains source-derived prose')
                    pm_here = current_pm
                    lab = obj.get("pmid") or obj.get("study")
                    if lab:
                        cand_pm = _extract_pmid(str(lab))
                        if cand_pm in non_oa_pmids:
                            pm_here = cand_pm
                    for k in list(obj.keys()):
                        k_pm = pm_here
                        if k in ("hop1_from", "hop1_to", "hop1_evidence", "origin_node", "bridge_node") and obj.get("hop1_study"):
                            c1 = _extract_pmid(str(obj.get("hop1_study")))
                            if c1 in non_oa_pmids:
                                k_pm = c1
                        elif k in ("hop2_to", "hop2_evidence", "terminal_node") and obj.get("hop2_study"):
                            c2 = _extract_pmid(str(obj.get("hop2_study")))
                            if c2 in non_oa_pmids:
                                k_pm = c2
                        obj[k] = _enforce_graph_shingle_budget(obj[k], k_pm, k)
                    return obj
                if isinstance(obj, list):
                    for idx in range(len(obj)):
                        obj[idx] = _enforce_graph_shingle_budget(obj[idx], current_pm, parent_key)
                    return obj
                if isinstance(obj,str) and current_pm not in non_oa_pmids:
                    source_pm=_unbound_source_match(obj,parent_key)
                    if source_pm:
                        digest=hashlib.sha256(obj.encode('utf-8')).hexdigest()
                        field=_public_redaction_field(parent_key or 'field')
                        return f'[REDACTED_NON_OA_EXCERPT: see PMID:{source_pm} field={field} quote_sha256={digest}]'
                if isinstance(obj, str) and current_pm in non_oa_pmids:
                    if _is_bounded_public_metadata(obj,parent_key):
                        return obj
                    if not dossier_source_cache.get(current_pm, ('', 0))[0] and _source_free_effect_identifier(obj, parent_key):
                        return obj
                    if _is_authentic_redaction_marker(obj) or _is_fixed_public_schema_metadata(obj, parent_key):
                        return obj
                    f_tok = _public_redaction_field(parent_key or "field")
                    if obj.startswith("[REDACTED_NON_OA_EXCERPT:"):
                        s_sha = hashlib.sha256(obj.encode("utf-8")).hexdigest()
                        return f"[REDACTED_NON_OA_EXCERPT: see PMID:{current_pm} field={f_tok} quote_sha256={s_sha}]"
                    full_src_t, c_budget = dossier_source_cache.get(current_pm, ("", 0))
                    v_len = _max_verbatim_shingle_chars(obj, full_src_t) if full_src_t else len(obj)
                    if v_len > 0:
                        if pmid_shingle_used.get(current_pm, 0) + v_len <= c_budget:
                            pmid_shingle_used[current_pm] = pmid_shingle_used.get(current_pm, 0) + v_len
                        else:
                            s_sha = hashlib.sha256(obj.encode("utf-8")).hexdigest()
                            return f"[REDACTED_NON_OA_EXCERPT: see PMID:{current_pm} field={f_tok} quote_sha256={s_sha}]"
                return obj

            _enforce_graph_shingle_budget(g_obj)
            graph_path.write_text(json.dumps(g_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

        public_graph=json.loads(graph_path.read_text(encoding='utf-8')) if graph_path.exists() else {}
        recomputed_totals=public_article_verbatim_totals(redacted_studies,public_graph,dossier_source_cache)
        if any(total>dossier_source_cache[pmid][1] for pmid,total in recomputed_totals.items()):
            raise ValueError('Combined study/graph excerpt budget post-condition failed')

        from jsonschema import Draft202012Validator
        export_schema = Draft202012Validator(json.loads(schema_bytes))
        for study_path in study_paths:
            exported_study = json.loads(study_path.read_text(encoding="utf-8"))
            export_schema.validate(exported_study)
            arm_ids = {arm["arm_id"] for arm in exported_study.get("arms", [])}
            for effect in exported_study.get("atomic_effects", []):
                if effect.get("focal_arm_id") not in arm_ids or (effect.get("comparator_arm_id") is not None and effect["comparator_arm_id"] not in arm_ids):
                    raise ValueError("Public export broke a study arm foreign key")

        # Post-condition shingle verification across all exported non-OA study JSONs
        for sp in study_paths:
            doc = json.loads(sp.read_text(encoding="utf-8"))
            pmid = str((doc.get("study_metadata") or {}).get("pmid", sp.stem))
            if pmid not in non_oa_pmids:
                continue
            full_src_text, char_budget = dossier_source_cache.get(pmid, ("", 0))
            title_str = str((doc.get("study_metadata") or {}).get("title") or "")

            def _walk_verbatim_chars(obj: object, parent_key: str = "",bibliographic: bool = False) -> int:
                if isinstance(obj, dict):
                    return sum(_walk_verbatim_chars(v, k,parent_key=='study_metadata') for k, v in obj.items())
                if isinstance(obj, list):
                    return sum(_walk_verbatim_chars(v, parent_key,bibliographic) for v in obj)
                if isinstance(obj, str):
                    if _is_bounded_public_metadata(obj,parent_key,bibliographic=bibliographic):
                        return 0
                    if not full_src_text and _source_free_effect_identifier(obj, parent_key):
                        return 0
                    if _is_authentic_redaction_marker(obj) or _is_fixed_public_schema_metadata(obj, parent_key):
                        return 0
                    if not full_src_text:
                        return len(obj)
                    return _max_verbatim_shingle_chars(obj, full_src_text)
                return 0

            _validate_imported_redaction_tokens(doc, pmid, None)
            total_v_chars = _walk_verbatim_chars(doc)
            if total_v_chars > char_budget:
                raise RuntimeError(
                    f"Public bundle post-condition violation for PMID:{pmid}: {total_v_chars} verbatim chars > budget {char_budget}"
                )

        # Rebind the prior plane to the actual public effects bytes, requiring all
        # numerical/semantic outputs to survive rights redaction unchanged.
        public_l2_path = work_dir / "data/extracted/layer2_regulatory_and_perturbation_atlas.json"
        public_effects_path = work_dir / "data/extracted/ring1_atomic_effects_v0.5.jsonl"
        if public_l2_path.exists() and public_effects_path.exists():
            import importlib.util
            plane_spec = importlib.util.spec_from_file_location("openislet_public_layer2_replay", work_dir / "scripts/run_layer2_perturbation_and_regulatory_plane.py")
            plane_module = importlib.util.module_from_spec(plane_spec)
            plane_spec.loader.exec_module(plane_module)
            plane_module.REPO_ROOT = work_dir
            public_l2 = json.loads(public_l2_path.read_text(encoding="utf-8"))
            original_plane = public_l2["stage0_to_stage7_perturbation_plane"]
            replayed_plane = plane_module.simulate_stage0_to_stage7_perturbation_plane(
                json.loads((work_dir / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text(encoding="utf-8")),
                public_l2["construct_sequence_and_grammar_verification"],
                public_l2["jaspar_pwm_regulatory_motif_scans"],
                effects_path=public_effects_path,
            )
            semantic = lambda plane: {k: v for k, v in plane.items() if k != "input_provenance"}
            if semantic(original_plane) != semantic(replayed_plane):
                raise ValueError("Rights transformation changed Layer 2 prior-plane semantics; refusing publication")
            public_l2["stage0_to_stage7_perturbation_plane"] = replayed_plane
            public_l2.setdefault("public_effects_provenance_rebinding", {
                "canonical_original_effects_sha256": original_plane["input_provenance"]["ring1_atomic_effects"]["sha256"],
                "status": "REPLAYED_PUBLIC_BYTES_NUMERICAL_AND_SEMANTIC_OUTPUTS_IDENTICAL",
                "qualification": "Source redaction changes byte identity; only input fingerprints are rebound. No new prediction or assay is claimed."
            })
            public_l2_path.write_text(json.dumps(public_l2, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

        # Update artifact digests in work_dir/manifest.json for the exported public bundle files
        pub_manifest = json.loads((work_dir / "manifest.json").read_text(encoding="utf-8"))
        pub_manifest["public_export_file_inventory"] = sorted(included - {"manifest.json"})
        pub_manifest["public_export_inventory_scope"] = "Exact shipped relative file paths; manifest.json and .openislet_public_bundle marker are the only omitted self/marker files"
        pub_manifest["public_arm_identifier_policy"] = {
            "format": "pmid<PMID>_arm<one-based arm position>",
            "arm_count": len(arm_identifier_map),
            "scope": "Deterministic opaque IDs; all exact foreign-key references transformed bijectively. Private source identifiers remain in the local audit tree."
        }
        for rel, entry in pub_manifest.get("artifacts", {}).items():
            fp = work_dir / rel
            if fp.exists():
                new_sha = sha256_file(fp)
                if new_sha != entry["sha256"]:
                    entry["local_audit_tree_sha256"] = entry["sha256"]
                    entry["sha256"] = new_sha
                    entry["bytes"] = fp.stat().st_size
        (work_dir / "manifest.json").write_text(
            json.dumps(pub_manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        if manifest_integrity_errors(work_dir,pub_manifest,public=True):
            raise ValueError('Transformed public bundle failed its exact inventory and digest post-condition')
        for study in redacted_studies:
            if public_token_errors(study,str(study['study_metadata']['pmid'])):
                raise ValueError('Public study contains an invalid opaque reference')
        if graph_path.exists() and public_token_errors(json.loads(graph_path.read_text(encoding='utf-8'))):
            raise ValueError('Public graph contains an invalid opaque reference')
        verify_public_bundle_privacy_guard(
            work_dir,
            source_root=ROOT,
            require_private_pattern_file=("public_export_file_inventory" not in manifest),
        )
    except Exception:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)
        raise
    finally:
        if policy_token is not None:_PUBLIC_POLICY_SNAPSHOT.reset(policy_token)
        sys.dont_write_bytecode = previous_bytecode
        sys.path[:] = previous_sys_path

    # Refuse publication if source inputs changed after their first validation.
    try:
        if _regular_file_bytes(ROOT, "manifest.json") != manifest_bytes:
            raise ValueError("Canonical manifest changed during export")
        for rel in sorted(included - {"manifest.json"}):
            _verify_pinned_source(ROOT, rel, manifest)
        for source in dossier_inputs:
            _verify_pinned_source(ROOT, source.relative_to(ROOT).as_posix(), manifest)
    except Exception:
        shutil.rmtree(staging_dir)
        raise

    backup_dir = None
    try:
        require_lease()
        backup_dir = Path(tempfile.mkdtemp(prefix=f".{final_dest_dir.name}.rollback-", dir=final_dest_dir.parent))
        backup_dir.rmdir()
        had_existing_bundle = final_dest_dir.exists()
        if had_existing_bundle:
            final_dest_dir.rename(backup_dir)
            require_lease()
        staging_dir.rename(final_dest_dir)
        require_lease()
    except Exception:
        owned=publication_lease is None or publication_lease.owns_lock()
        if owned and backup_dir is not None and backup_dir.exists() and not final_dest_dir.exists():
            backup_dir.rename(final_dest_dir)
        if staging_dir.exists():
            shutil.rmtree(staging_dir)
        raise
    if backup_dir is not None and backup_dir.exists():
        require_lease()
        shutil.rmtree(backup_dir)
    return sorted(included)


def verify_public_bundle_snapshot(bundle_root: Path | None = None) -> dict[str, int | bool]:
    """Verify shipped public bundle files against manifest.json without requiring private dossiers or decision logs."""
    root = Path(bundle_root) if bundle_root is not None else ROOT
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"manifest.json not found at {manifest_path}")
    mf = json.loads(manifest_path.read_text(encoding="utf-8"))
    is_pub = (root / ".openislet_public_bundle").exists() or "public_export_file_inventory" in mf
    errs = manifest_integrity_errors(root, mf, public=is_pub)
    if errs:
        raise ValueError(f"Manifest snapshot integrity failed ({len(errs)} error(s)): {errs[:5]}")
    n_verified = len(mf.get("public_export_file_inventory") or mf.get("artifacts") or {})
    return {"passed": True, "verified_files": n_verified, "public_bundle_mode": is_pub}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate release manifest and optional public bundle.")
    parser.add_argument("--export-public-bundle", type=str, default=None, help="Directory to export rights-aware public bundle")
    parser.add_argument(
        "--regenerate-repo-manifest",
        action="store_true",
        help="Also regenerate in-repo manifest.json and JSONL aliases when --export-public-bundle is specified",
    )
    parser.add_argument(
        "--verify-public-bundle",
        action="store_true",
        help="Verify pinned public bundle file hashes and inventory against manifest.json without requiring private caches",
    )
    args = parser.parse_args()
    if args.verify_public_bundle or (
        not args.export_public_bundle
        and not args.regenerate_repo_manifest
        and ((ROOT / ".openislet_public_bundle").exists() or not (ROOT / "docs" / "PROJECT-DECISIONS-LOG.md").exists())
    ):
        res = verify_public_bundle_snapshot(ROOT)
        print(
            f"Verified {res['verified_files']} pinned artifacts against manifest.json "
            f"(public_bundle_mode={res['public_bundle_mode']}). "
            "Note: Full manifest regeneration from raw dossiers requires the private audit cache "
            "(data/dossiers/*.json, docs/PROJECT-DECISIONS-LOG.md)."
        )
    else:
        if not args.export_public_bundle or args.regenerate_repo_manifest or not (ROOT / "manifest.json").exists():
            main()
        if args.export_public_bundle:
            exported_files = export_public_bundle(Path(args.export_public_bundle))
            print(f"Exported {len(exported_files)} allowlisted files to {args.export_public_bundle}")
