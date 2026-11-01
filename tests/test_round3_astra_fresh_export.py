"""Fresh release-boundary, acceptance positive-control and Markdown regressions."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import extract_ring1_corpus as extraction
from scripts import generate_release_manifest as export
from scripts import render_documentation as renderer
from tests import test_adversarial_acceptance as acceptance


def require(condition, message):
    if not condition:
        raise AssertionError(message)


PRIVATE_PLANS = (
    'docs/OPENISLET-2026-2030-MASTER-TIMELINE.md',
    'docs/PROJECT-DECISIONS-LOG.md',
)


def facts():
    return json.loads((ROOT / 'data/extracted/documentation_facts_v0.8.0.json').read_text())


@pytest.fixture
def source(tmp_path, monkeypatch):
    """Self-contained source fixture works in both private and public checkouts."""
    root = tmp_path.resolve() / 'source'
    files = [
        'scripts/rank_islet1_stacks.py', 'scripts/engineering_priors.py',
        'scripts/public_excerpt_policy.py', 'schemas/extraction_schema_v0.2.json',
        'data/omics/target_prior_table.json',
    ]
    for relative in files:
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, target)
    dossier = {
        'pmid': '99999999', 'title': 'Synthetic B2M knockout islet assay',
        'abstract': 'B2M knockout islets showed viability of 80 percent versus 20 percent in controls.',
        'pub_year': 2026, 'authors': ['Synthetic Tester'],
        'first_author': 'Synthetic Tester', 'senior_author': 'Synthetic Tester',
        'journal': 'Synthetic fixture',
    }
    study = extraction.extract_one_dossier(dossier)
    payloads = {
        'LICENSE': 'Synthetic license',
        'docs/PUBLIC-SCIENCE.md': '# Synthetic public science\n',
        PRIVATE_PLANS[0]: '# Private timeline sentinel\n',
        PRIVATE_PLANS[1]: '# Private decisions sentinel\n',
        'data/dossiers/99999999.json': json.dumps(dossier),
        'data/extracted/99999999.json': json.dumps(study),
        'data/extracted/cross_study_discovery_graph.json': '{}',
        'data/extracted/ring1_study_dossiers_v0.5.jsonl': '',
        'data/extracted/ring1_atomic_effects_v0.5.jsonl': '',
    }
    for relative, content in payloads.items():
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content)
        files.append(relative)
    priv_pat = root / 'data/raw/private_export_guard_patterns.json'
    priv_pat.parent.mkdir(parents=True, exist_ok=True)
    priv_pat.write_text('{"regex_patterns":{"synthetic_guard":"SYNTHETIC_GUARD_TOKEN"}}')
    manifest = {'artifacts': {}, 'public_redistribution_allowlist_spec': {
        'public_package_included_paths': [p for p in files if not p.startswith(('docs/', 'data/dossiers/'))] + ['docs/*.md', 'manifest.json'],
        'local_research_cache_excluded_from_public_bundle': ['data/dossiers/*.json', *PRIVATE_PLANS],
    }}
    def seal():
        for relative in files:
            payload = (root / relative).read_bytes()
            manifest['artifacts'][relative] = {'sha256': hashlib.sha256(payload).hexdigest(), 'bytes': len(payload)}
        (root / 'manifest.json').write_text(json.dumps(manifest))
    seal()
    monkeypatch.setattr(export, 'ROOT', root)
    return root, manifest, seal


def test_explicit_private_plans_are_excluded_by_real_export_and_source_free_reexport(source, tmp_path, monkeypatch):
    root, _, _ = source
    originals = {p: (root / p).read_bytes() for p in PRIVATE_PLANS}
    destination = tmp_path.resolve() / 'public'
    export.export_public_bundle(destination)
    require(all(not (destination / p).exists() for p in PRIVATE_PLANS), 'private planning document shipped')
    require((destination / 'docs/PUBLIC-SCIENCE.md').is_file(), 'public science was incorrectly excluded')
    require(all((root / p).read_bytes() == value for p, value in originals.items()), 'source document changed')
    monkeypatch.setattr(export, 'ROOT', destination)
    again = tmp_path.resolve() / 'again'
    export.export_public_bundle(again)
    require(all(not (again / p).exists() for p in PRIVATE_PLANS), 'source-free re-export re-enabled private planning')
    manifest = json.loads((again / 'manifest.json').read_text())
    require(export.manifest_integrity_errors(again, manifest, True) == [], 're-export is not a complete pinned bundle')


@pytest.mark.parametrize('relative', PRIVATE_PLANS)
@pytest.mark.parametrize('public_source', [False, True])
def test_policy_mutation_cannot_reinclude_private_plans(source, tmp_path, monkeypatch, relative, public_source):
    root, manifest, _ = source
    if public_source:
        public = tmp_path.resolve() / 'public'
        export.export_public_bundle(public)
        (public / relative).write_bytes((root / relative).read_bytes())
        root = public
        manifest = json.loads((root / 'manifest.json').read_text())
        monkeypatch.setattr(export, 'ROOT', root)
    manifest['public_redistribution_allowlist_spec']['local_research_cache_excluded_from_public_bundle'] = ['data/dossiers/*.json']
    (root / 'manifest.json').write_text(json.dumps(manifest))
    destination = tmp_path.resolve() / 'rejected'
    with pytest.raises(ValueError, match='planning documents'):
        export.export_public_bundle(destination)
    require(not destination.exists(), 'forbidden file reached publication')


@pytest.mark.parametrize('relative', [*PRIVATE_PLANS, 'DOCS/project-decisions-log.MD', 'DATA/DOSSIERS/private.JSON'])
def test_forbidden_public_paths_include_case_aliases(relative):
    require(export._forbidden_public_path(relative), 'case alias escaped the fixed private boundary')


@pytest.mark.parametrize('relative', ['../secret', '/absolute', 'docs/../secret', 'docs//secret', 'docs/./secret', 'docs\\secret', 'docs/é.txt'])
def test_unshipped_manifest_paths_cannot_escape_the_release_tree(source, tmp_path, relative):
    root, manifest, _ = source
    manifest['artifacts'][relative] = {'sha256': '0' * 64, 'bytes': 0}
    (root / 'manifest.json').write_text(json.dumps(manifest))
    require(bool(export.manifest_integrity_errors(root, manifest, True)), 'unshipped bad path passed inventory verification')
    with pytest.raises(ValueError, match='artifact inventory'):
        export.export_public_bundle(tmp_path.resolve() / 'rejected')


def test_duplicate_case_manifest_paths_have_no_distinct_portable_identity(source, tmp_path):
    root, manifest, _ = source
    manifest['artifacts']['license'] = copy.deepcopy(manifest['artifacts']['LICENSE'])
    (root / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='duplicate-portable'):
        export.export_public_bundle(tmp_path.resolve() / 'rejected')


def test_source_free_export_cannot_publish_short_unverifiable_prose(source, tmp_path):
    root, manifest, seal = source
    graph = {'study': 'PMID:' + '9' * 8, 'note': '秘' * 300}
    (root / 'data/extracted/cross_study_discovery_graph.json').write_text(json.dumps(graph))
    seal()
    shutil.rmtree(root / 'data/dossiers')
    manifest['artifacts'].pop('data/dossiers/99999999.json')
    manifest['public_export_file_inventory'] = [p for p in manifest['artifacts'] if p not in PRIVATE_PLANS]
    (root / 'manifest.json').write_text(json.dumps(manifest))
    destination = tmp_path.resolve() / 'public'
    export.export_public_bundle(destination)
    output = json.loads((destination / 'data/extracted/cross_study_discovery_graph.json').read_text())
    require(output['note'] != graph['note'], 'under-six-word prose escaped the zero-source budget')
    require(not export.public_token_errors(output), 'redaction token lost its field or PMID binding')


@pytest.mark.parametrize('text', ['café βeta γamma delta epsilon zeta', 'cafe\u0301 βeta γamma delta epsilon zeta', 'ＣＡＦＥ βeta γamma delta epsilon zeta', 'ca\u200dfe βeta γamma delta epsilon zeta'])
def test_unicode_shingles_count_original_codepoints(text):
    source = 'café βeta γamma delta epsilon zeta'
    require(export._max_verbatim_shingle_chars({'note': text}, source) == len(text), 'normalization changed the counted character budget')


def test_source_backed_effect_identifier_slug_still_counts_source_overlap():
    value = 'alpha_beta_gamma_delta_epsilon_zeta'
    require(export._max_verbatim_shingle_chars({'effect_id': value}, value.replace('_', ' ')) == len(value),
            'source-free identifier syntax became an exemption against available source text')


@pytest.mark.parametrize('field,value', [
    ('non_oa_per_article_aggregate_verbatim_ceiling_characters', True),
    ('non_oa_per_article_aggregate_verbatim_ceiling_characters', 241),
    ('non_oa_per_article_max_abstract_share', False),
    ('non_oa_per_article_max_abstract_share', float('nan')),
    ('non_oa_per_article_max_abstract_share', float('inf')),
    ('non_oa_per_article_max_abstract_share', .26),
])
def test_exporter_and_acceptance_share_finite_tight_excerpt_limits(field, value):
    with pytest.raises(ValueError):
        export.public_excerpt_limits({field: value})
    require(export.public_excerpt_limits({'non_oa_per_article_aggregate_verbatim_ceiling_characters': 12,
                                        'non_oa_per_article_max_abstract_share': .125}) == (12, .125),
            'stricter release policy was silently relaxed')


def test_semantic_mutations_require_a_valid_baseline(monkeypatch):
    doc, dossier, effect_id = acceptance.synthetic_semantic_probe_source()
    def always_invalid(*args):
        raise ValueError('independent baseline defect')
    monkeypatch.setattr(acceptance.ext, 'validate_study_semantics', always_invalid)
    result = acceptance.semantic_mutation_controls(doc, dossier, effect_id)
    require(not result['baseline_valid'] and not result['mutations'], 'unrelated baseline failure counted as successful mutation rejection')


def test_public_semantic_positive_control_and_all_five_distinct_mutations():
    doc, dossier, effect_id = acceptance.synthetic_semantic_probe_source()
    result = acceptance.semantic_mutation_controls(doc, dossier, effect_id)
    require(result['baseline_valid'], 'synthetic public control is not valid')
    require(len(result['mutations']) == 5 and all(result['mutations'].values()), 'targeted corruption passed')
    errors = result['mutation_errors']
    for key, token in [('wrong_numeric_lrr', 'ln('), ('forged_dossier_hash', 'dossier_sha256'),
                       ('numeric_value_with_qualitative_statistic_type', 'statistic_type'),
                       ('reversed_orientation_in_same_polarity_class', 'orientation_sign'),
                       ('missing_well_formed_comparator_id', 'comparator_arm_id')]:
        require(token in errors[key], 'mutation rejected for an unrelated reason: ' + key)


@pytest.mark.parametrize('literal', [
    '- ```text\n  GSG-tPT2A MVDGTLLLL\n  ```',
    '1. ~~~text\n   GSG-tPT2A MVDGTLLLL\n   ~~~~',
    '> - ```text\n>   GSG-tPT2A MVDGTLLLL\n>   ```',
    '- > ```text\n  > CpG-free\n  > ```',
    '10. > ```text\n    > CpG-free\n    > ```',
    '- > - ```text\n  >   CpG-free\n  >   ```',
    '> - > ```text\n>   > CpG-free\n>   > ```',
    '<script>\nconst alias = "GSG-tPT2A";\n</script>',
    '<pre>\nGSG-tPT2A MVDGTLLLL\n</pre>',
])
def test_container_code_and_raw_html_literals_survive_repeat_render(literal):
    text = '# Test\n\n' + literal + '\n\nAfter.\n'
    rendered = renderer.render_text(text, facts(), 'fixture.md')
    require(literal in rendered and 'After.' in rendered, 'literal bytes or following authored content changed')
    require(renderer.render_text(rendered, facts(), 'fixture.md') == rendered, 'literal block changed on second render')


@pytest.mark.parametrize('heading', ['Title\n=====\n', 'Multi-line\nheading\n-------\n'])
def test_setext_heading_is_not_split_by_generated_blocks(heading):
    rendered = renderer.render_text(heading + '\nBody.\n', facts(), 'fixture.md')
    require(heading in rendered, 'generated facts split a setext heading from its underline')
    require(renderer.render_text(rendered, facts(), 'fixture.md') == rendered, 'setext heading render is not idempotent')


@pytest.mark.parametrize('link', [
    '[linked](https://example.org/CpG-free)',
    '[linked](https://example.org/CpG-free_(nested) "GSG-tPT2A")',
    '[linked](<https://example.org/CpG-free>)',
    '<https://example.org/CpG-free>',
    '[reference]: https://example.org/CpG-free "MVDGTLLLL"',
    '[linked](https://example.org/233/233)',
    '<https://example.org/233/233>',
    '[reference]:\n  https://example.org/CpG-free',
    '[reference]: https://example.org/path\n  "CpG-free title"',
    '[reference]:\n  <https://example.org/CpG-free_(nested)>\n  "CpG-free title"',
    '[reference]: https://example.org/CpG-free_(nested)\n"CpG-free title"',
    '[linked](https://example.org/233/233 "<a href=\'CpG-free\'>")',
    '<a href="https://example.org/CpG-free/233/233">source</a>',
    '<a\n href="https://example.org/CpG-free/233/233" title="GSG-tPT2A">source</a>',
    '<!--\nCpG-free 233/233\n-->',
])
def test_link_destinations_are_not_scientific_replacement_inputs(link):
    rendered = renderer.render_text('# Test\n\n' + link + '\n', facts(), 'fixture.md')
    require(link in rendered, 'rendering rewrote a URL, title or reference definition')
    require(renderer.render_text(rendered, facts(), 'fixture.md') == rendered, 'literal syntax changed on second render')


def test_reference_definition_protection_does_not_hide_following_prose():
    text = '# Test\n\n[ref]: https://example.org/CpG-free\n\nCpG-free\n'
    rendered = renderer.render_text(text, facts(), 'fixture.md')
    require('[ref]: https://example.org/CpG-free' in rendered, 'reference destination changed')
    require('\nCpG-containing candidate; anti-silencing pending\n' in rendered, 'ordinary prose was incorrectly protected')


@pytest.mark.parametrize('field,value', [
    ('authors', [{'name': 'Synthetic Author', 'full_text': 'private source' * 50}]),
    ('keywords', {'attachment': 'private source' * 50}),
    ('authors', [['private source']]),
    ('mesh_terms', [False]),
    ('title', {'full_text': 'private source'}),
    ('has_abstract', {'full_text': 'private source'}),
    ('pub_year', True),
    ('abstract_sha256', 'private source'),
    ('keywords', ['private source' * 50]),
])
def test_discovery_metadata_requires_bounded_flat_bibliographic_schema(field, value):
    with pytest.raises(ValueError, match='Discovery metadata field ' + field):
        export._public_discovery_record({'pmid': '99999999', field: value})


def test_discovery_metadata_retains_valid_typed_values_and_is_idempotent():
    source = {'pmid': '99999999', 'pmcid': None, 'doi': '10.1234/example', 'pub_year': 2026,
              'title': 'Synthetic scientific title', 'journal': 'Synthetic journal',
              'first_author': 'Synthetic Author', 'senior_author': 'Synthetic Consortium',
              'authors': ['Synthetic Author', 'Synthetic Consortium'],
              'mesh_terms': ['Islets of Langerhans'], 'keywords': ['beta cells'],
              'publication_types': ['Journal Article'], 'matched_pillars': ['Synthetic pillar'],
              'has_open_pmc_fulltext': False, 'source_db': 'Synthetic database',
              'abstract': 'Private synthetic abstract.'}
    original = copy.deepcopy(source)
    public = export._public_discovery_record(source)
    require(source == original, 'metadata projection mutated its input')
    require(all(public[key] == value for key, value in source.items() if key != 'abstract'), 'valid bibliographic metadata was dropped or changed')
    require('abstract' not in public and public['has_abstract'] is True, 'source abstract survived')
    require(public['abstract_sha256'] == hashlib.sha256(source['abstract'].encode()).hexdigest(), 'abstract digest was not retained')
    require(export._public_discovery_record(public) == public, 'source-free metadata projection is not idempotent')


@pytest.mark.parametrize('field,value', [
    ('authors', [{'name': 'Synthetic Author', 'full_text': 'Private source. ' * 50}]),
    ('keywords', {'attachment': 'Private source. ' * 50}),
])
def test_actual_export_rejects_nested_discovery_text_without_replacing_previous_bundle(source, tmp_path, field, value):
    root, manifest, _ = source
    destination = tmp_path.resolve() / 'public'
    export.export_public_bundle(destination)
    previous = {str(p.relative_to(destination)): p.read_bytes() for p in destination.rglob('*') if p.is_file()}
    relative = 'data/raw/corpus_manifest_all_sources.jsonl'
    target = root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps({'pmid': '99999999', field: value}) + '\n')
    manifest['artifacts'][relative] = {'sha256': hashlib.sha256(target.read_bytes()).hexdigest(), 'bytes': target.stat().st_size}
    manifest['public_redistribution_allowlist_spec']['public_package_included_paths'].append(relative)
    (root / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='Discovery metadata field ' + field):
        export.export_public_bundle(destination)
    after = {str(p.relative_to(destination)): p.read_bytes() for p in destination.rglob('*') if p.is_file()}
    require(after == previous, 'failed discovery projection changed the published bundle')


def test_inline_alias_replacement_preserves_code_span_delimiters():
    rendered = renderer.render_text('# Test\n\n`minimal 455 bp CBX3-A2UCOE`\n', facts(), 'fixture.md')
    require('`' + renderer.CANDIDATE + ' (legacy CBX3-A2UCOE alias; anti-silencing function pending)`' in rendered,
            'inline alias replacement consumed a code delimiter')


@pytest.mark.parametrize('field', ['superseded_by_pmid', 'reconciled_from_superseded_preprint_pmid'])
@pytest.mark.parametrize('value,valid', [('38661000', True), ('1', True), ('9' * 9, True),
                                       ('0', False), ('038661000', False), ('３８６６１０００', False),
                                       ('38661000 private prose', False), ('9' * 10, False),
                                       (38661000, False), (True, False)])
def test_supersession_metadata_exempts_only_canonical_pmid_values(field, value, valid):
    require(export._is_fixed_public_schema_metadata(value, field) is valid,
            'supersession metadata exemption accepted malformed identity or lost valid identity')


def test_source_free_reexport_preserves_canonical_supersession_metadata(source, tmp_path, monkeypatch):
    root, _, seal = source
    relative = 'data/extracted/99999999.json'
    study = json.loads((root / relative).read_text())
    study['study_metadata']['superseded_by_pmid'] = '8' * 8
    study['study_metadata']['reconciled_from_superseded_preprint_pmid'] = '7' * 8
    (root / relative).write_text(json.dumps(study))
    seal()
    public = tmp_path.resolve() / 'public'
    export.export_public_bundle(public)
    monkeypatch.setattr(export, 'ROOT', public)
    again = tmp_path.resolve() / 'again'
    export.export_public_bundle(again)
    output = json.loads((again / relative).read_text())['study_metadata']
    require(output['superseded_by_pmid'] == '8' * 8, 'canonical supersession PMID was redacted on source-free re-export')
    require(output['reconciled_from_superseded_preprint_pmid'] == '7' * 8, 'reconciliation PMID was redacted on source-free re-export')
    for line in (again / 'data/extracted/ring1_atomic_effects_v0.5.jsonl').read_text().splitlines():
        require(json.loads(line)['superseded_by_pmid'] == '8' * 8, 'public JSONL lost typed supersession identity')


@pytest.mark.parametrize('unknown', [None, 'UNAVAILABLE'])
def test_public_renderer_preserves_unavailable_rank_facts_without_a_zero_default(unknown):
    current = facts()
    current['layer1']['audit_ledger'] = {'search': {'all_stacks': unknown}}
    current['layer1']['ranked_architectures'] = []
    current['layer1']['sensitivity'] = {'uniform_grid_pooled_pediatric_gate_logit_variance_decomposition':
                                     {'structural_share_of_component_variance_sum': unknown}}
    text = renderer.render_text('# Public facts\n\nBody.\n', current, 'fixture.md')
    label = 'unavailable' if unknown is None else unknown
    require('Layer 1 evaluates ' + label + ' prior-anchored stacks' in text, 'unknown stack count was coerced')
    require('Structural share of component variance is ' + label + '%' in text, 'unknown variance share became numeric zero')
    require(renderer.render_text(text, current, 'fixture.md') == text, 'public unavailable render is not idempotent')


def test_derive_and_render_without_private_ranking_cache(tmp_path):
    root = tmp_path / 'source-free'
    path = root / 'data/raw/corpus_manifest_all_sources.jsonl'
    path.parent.mkdir(parents=True)
    path.write_text('')
    current = renderer.derive_facts(root, 0)
    text = renderer.render_text('# Public facts\n\nBody.\n', current, 'fixture.md')
    require('Layer 1 evaluates unavailable prior-anchored stacks' in text, 'missing private ranking cache crashed or became a zero count')
    require('Structural share of component variance is unavailable%' in text, 'missing private sensitivity was represented as measured zero')


def test_synthetic_graph_fixture_does_not_add_a_fabricated_bibliographic_citation():
    from tests.test_citations_resolve import PMID_RE
    require(not any(match[1] == '9' * 8 for match in PMID_RE.finditer(Path(__file__).read_text())),
            'synthetic graph fixture was written as a real bibliographic citation')
