"""Independent numeric and provenance controls; synthetic inputs are not study labels."""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import io
import json
import math
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import build_cross_study_graph as graph
from scripts import extract_ring1_corpus as ext
from scripts import harvest_corpus as harvest
from scripts import rank_islet1_stacks as rank


LARGE = '1' + '0' * 308
SMALL = '0.' + '0' * 307 + '1'


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('dispersions', [False, True])
def test_log_ratio_extremes_match_decimal_oracle_without_dividing(reverse, dispersions):
    treatment, control = (SMALL, LARGE) if reverse else (LARGE, SMALL)
    endpoint = 'cell death' if reverse else 'cell viability'
    suffix = ' ± 0' if dispersions else ''
    quote = f'B2M knockout improved {endpoint}: {treatment}{suffix} ng/mL versus {control}{suffix} ng/mL.'
    result = ext.parse_quantitative_numbers(quote, 'Protective_or_Improved')
    with localcontext() as ctx:
        ctx.prec = 80
        expected = float((Decimal(treatment) / Decimal(control)).ln())
    assert result[0] == float(treatment)
    assert result[2] == float(control)
    assert math.isfinite(result[5])
    assert result[5] == pytest.approx(expected, abs=5e-5)


@pytest.mark.parametrize('token', ['1' + '0' * 309, '0.' + '0' * 400 + '1'])
@pytest.mark.parametrize('template', [
    'B2M knockout remained detectable at {token} ng/mL.',
    'B2M knockout improved cell viability: {token} ± 1 ng/mL versus 1 ± 1 ng/mL.',
])
def test_unrepresentable_source_operands_are_quarantined(token, template):
    result = ext.parse_quantitative_numbers(template.format(token=token), 'Protective_or_Improved')
    assert result == (None, None, None, None, 'quarantined_nonfinite_numeric', None)


def test_extreme_source_bound_lrr_passes_semantic_validation():
    source = {'pmid': '12345678', 'title': 'Synthetic B2M knockout arithmetic fixture',
              'abstract': 'B2M knockout improved cell viability: 80% versus 20%.', 'fulltext_data_sections': {}}
    doc = ext.extract_one_dossier(source)
    quote = f'B2M knockout improved cell viability: {LARGE} ng/mL versus {SMALL} ng/mL.'
    source['abstract'] = quote
    effect = doc['atomic_effects'][0]
    control = copy.deepcopy(doc['arms'][0])
    control.update(arm_id='synthetic_explicit_control', is_reference_control=True)
    doc['arms'].append(control)
    effect.update(verbatim_evidence_quote=quote, treatment_mean=1e308, control_mean=1e-308,
                  log_response_ratio_or_smd=1418.3924, measurement_unit='ng/mL',
                  comparator_arm_id='synthetic_explicit_control', statistic_type='log_response_ratio_paired_means',
                  source_span=f'abstract[0:{len(quote)}]|norm=exact',
                  dossier_sha256=ext._compute_dossier_sha256('12345678', source))
    ext.validate_study_semantics(doc, source)
    effect['treatment_mean'] = 10 ** 500
    with pytest.raises(ValueError, match='finite numeric domain'):
        ext.validate_study_semantics(doc, source)


@pytest.mark.parametrize('value', [0, -1, True, float('nan'), float('inf'), 10 ** 500])
def test_log_ratio_domain_has_explicit_value_errors(value):
    with pytest.raises(ValueError):
        ext._log_response_ratio(value, 1)


def _xml(pmids):
    return ('<PubmedArticleSet>' + ''.join(
        f'<PubmedArticle><MedlineCitation><PMID>{pid}</PMID></MedlineCitation></PubmedArticle>'
        for pid in pmids) + '</PubmedArticleSet>').encode()


@pytest.mark.parametrize('returned', [[], ['1234'], ['1234', '1234'], ['1234', '9999'], ['1234', '5678', '9999']])
def test_efetch_bad_batch_never_advances_or_repairs_checkpoint(monkeypatch, tmp_path, returned):
    target = tmp_path / 'checkpoint.jsonl'
    original = b'{"pmid":"1111"}'
    target.write_bytes(original)
    monkeypatch.setattr(harvest.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(_xml(returned)))
    monkeypatch.setattr(harvest.time, 'sleep', lambda *a: None)
    with pytest.raises(ValueError, match='exactly the requested unique PMIDs'):
        harvest.efetch_metadata(['1111', '1234', '5678'], {}, target)
    assert target.read_bytes() == original


def test_efetch_valid_reordered_batch_preserves_source_membership(monkeypatch, tmp_path):
    target = tmp_path / 'checkpoint.jsonl'
    target.write_text('{"pmid":"1111"}')
    monkeypatch.setattr(harvest.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(_xml(['5678', '1234'])))
    monkeypatch.setattr(harvest.time, 'sleep', lambda *a: None)
    assert harvest.efetch_metadata(['1111', '1234', '5678'], {'5678': ['unit-pillar']}, target) == 3
    rows = [json.loads(line) for line in target.read_text().splitlines()]
    assert [row['pmid'] for row in rows] == ['1111', '5678', '1234']
    assert rows[1]['matched_pillars'] == ['unit-pillar']


def test_efetch_unrequested_checkpoint_rejected_before_network(monkeypatch, tmp_path):
    target = tmp_path / 'checkpoint.jsonl'
    target.write_text('{"pmid":"9999"}\n')
    def no_network(*args, **kwargs):
        raise AssertionError('network reached before checkpoint validation')
    monkeypatch.setattr(harvest.urllib.request, 'urlopen', no_network)
    with pytest.raises(ValueError, match='outside the requested query snapshot'):
        harvest.efetch_metadata(['1234'], {}, target)
    assert target.read_text() == '{"pmid":"9999"}\n'


def test_efetch_commits_only_complete_prior_batches(monkeypatch, tmp_path):
    target = tmp_path / 'checkpoint.jsonl'
    responses = iter([_xml([str(i) for i in range(1, 201)]), _xml([])])
    monkeypatch.setattr(harvest.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(next(responses)))
    monkeypatch.setattr(harvest.time, 'sleep', lambda *a: None)
    with pytest.raises(ValueError, match='checkpoint not advanced'):
        harvest.efetch_metadata([str(i) for i in range(1, 202)], {}, target)
    assert [json.loads(line)['pmid'] for line in target.read_text().splitlines()] == [str(i) for i in range(1, 201)]


@pytest.mark.parametrize('second', [{'count': '2', 'idlist': []}, {'count': '4', 'idlist': ['3']}, {'count': '-1', 'idlist': []}, {'count': True, 'idlist': []}])
def test_esearch_changed_or_invalid_count_cannot_signal_completion(monkeypatch, second):
    responses = iter([{'count': '3', 'idlist': ['1', '2']}, second])
    monkeypatch.setattr(harvest.urllib.request, 'urlopen', lambda *a, **k: io.BytesIO(json.dumps({'esearchresult': next(responses)}).encode()))
    monkeypatch.setattr(harvest.time, 'sleep', lambda *a: None)
    with pytest.raises(ValueError):
        harvest.esearch_pmids('synthetic query', retmax=2)


@pytest.mark.parametrize('change', ['duplicate_arm', 'duplicate_effect', 'unknown_focal', 'self_comparator'])
def test_schema_valid_but_ambiguous_graph_records_fail_before_output(monkeypatch, tmp_path, change):
    doc = json.loads((ROOT/'data/extracted/10728379.json').read_text())
    if change == 'duplicate_arm':
        doc['arms'].append(copy.deepcopy(doc['arms'][0]))
    elif change == 'duplicate_effect':
        doc['atomic_effects'].append(copy.deepcopy(doc['atomic_effects'][0]))
    elif change == 'unknown_focal':
        doc['atomic_effects'][0]['focal_arm_id'] = 'missing_arm'
    else:
        doc['atomic_effects'][0]['comparator_arm_id'] = doc['atomic_effects'][0]['focal_arm_id']
    schema = json.loads((ROOT/'schemas/extraction_schema_v0.2.json').read_text())
    graph.validate_against_schema(doc, schema)
    (tmp_path/'10728379.json').write_text(json.dumps(doc))
    target = tmp_path/'cross_study_discovery_graph.json'
    target.write_text('PRIOR VERIFIED GRAPH')
    monkeypatch.setattr(graph, 'EXTRACTED_DIR', tmp_path)
    with pytest.raises(ValueError):
        graph.main()
    assert target.read_text() == 'PRIOR VERIFIED GRAPH'


def test_current_corpus_graph_identifiers_and_relations_are_valid():
    studies = []
    for path in sorted((ROOT/'data/extracted').glob('*.json')):
        if path.stem.isdigit():
            study = json.loads(path.read_text())
            graph._validate_graph_references(study, path.name)
            studies.append(study)
    graph._validate_graph_corpus_identity(studies)


def _graph_study_pair():
    first = json.loads((ROOT/'data/extracted/10728379.json').read_text())
    second = json.loads(json.dumps(first).replace('pmid10728379_', 'pmid10728380_'))
    second['study_metadata']['pmid'] = '10728380'
    return first, second


def _set_graph_identifier(study, collection, key, value):
    previous = study[collection][0][key]
    study[collection][0][key] = value
    if collection == 'arms':
        for effect in study['atomic_effects']:
            for reference in ('focal_arm_id', 'comparator_arm_id'):
                if effect[reference] == previous:
                    effect[reference] = value


@pytest.mark.parametrize('collection,key', [('arms', 'arm_id'), ('atomic_effects', 'effect_id')])
@pytest.mark.parametrize('defect', ['foreign_namespace', 'custom_id_reuse'])
def test_cross_study_graph_identifiers_cannot_reuse_or_claim_another_study(monkeypatch, tmp_path, collection, key, defect):
    first, second = _graph_study_pair()
    if defect == 'foreign_namespace':
        _set_graph_identifier(second, collection, key, first[collection][0][key])
    else:
        for study in (first, second):
            _set_graph_identifier(study, collection, key, f'curated_shared_{key}')
    schema = json.loads((ROOT/'schemas/extraction_schema_v0.2.json').read_text())
    for study in (first, second):
        graph.validate_against_schema(study, schema)
        graph._validate_graph_references(study, 'synthetic graph fixture')
        (tmp_path/f"{study['study_metadata']['pmid']}.json").write_text(json.dumps(study))
    output = tmp_path/'cross_study_discovery_graph.json'
    output.write_text('PRIOR VERIFIED GRAPH')
    monkeypatch.setattr(graph, 'EXTRACTED_DIR', tmp_path)
    with pytest.raises(ValueError, match='namespace|reused across studies'):
        graph.main()
    assert output.read_text() == 'PRIOR VERIFIED GRAPH'


def test_distinct_curated_and_namespaced_graph_identifiers_preserve_valid_studies(monkeypatch, tmp_path):
    first, second = _graph_study_pair()
    _set_graph_identifier(first, 'atomic_effects', 'effect_id', 'curated_first_study_effect')
    _set_graph_identifier(second, 'atomic_effects', 'effect_id', 'curated_second_study_effect')
    for study in (first, second):
        (tmp_path/f"{study['study_metadata']['pmid']}.json").write_text(json.dumps(study))
    monkeypatch.setattr(graph, 'EXTRACTED_DIR', tmp_path)
    graph.main()
    output = json.loads((tmp_path/'cross_study_discovery_graph.json').read_text())
    assert output['validated_studies_count'] == 2
    assert output['total_atomic_effects_extracted'] == sum(len(study['atomic_effects']) for study in (first, second))


@pytest.mark.parametrize('profiles', [{}, {'INS': {'n_empirical_studies': 1}}])
def test_exact_tie_order_is_independent_of_module_insertion_order(monkeypatch, profiles):
    before = rank.run_combinatorial_multi_gate_search(profiles)
    for name in ('ALLO_CLOAK_MODULES', 'AUTOANTIGEN_MODULES', 'CYTOKINE_GRANZYME_MODULES', 'ER_AMYLOID_MODULES', 'NICHE_VASCULAR_MODULES', 'SAFETY_SWITCH_MODULES'):
        monkeypatch.setattr(rank, name, dict(reversed(list(getattr(rank, name).items()))))
    after = rank.run_combinatorial_multi_gate_search(profiles)
    assert after == before


@pytest.mark.parametrize('weights', [(-.1, 1.1), (2., 2.), (0., 0.), (.2, .7)])
def test_invalid_cohort_simplex_cannot_produce_a_utility(monkeypatch, weights):
    secondary = rank.COMPLETE_GATE_COEFFICIENT_REGISTRY['secondary_ranking_tie_break_and_display_constants']
    for cohort, value in zip(('pediatric', 'adult'), weights):
        monkeypatch.setitem(secondary[f'composite_utility_{cohort}_weight'], 'parameter_value', value)
    with pytest.raises(ValueError, match='sum to one'):
        rank.run_combinatorial_multi_gate_search({})


def test_negative_quadratic_discount_cannot_increase_the_floor(monkeypatch):
    secondary = rank.COMPLETE_GATE_COEFFICIENT_REGISTRY['secondary_ranking_tie_break_and_display_constants']
    monkeypatch.setitem(secondary['lower_10pct_cred_bound_b_edit_quadratic_discount'], 'parameter_value', -.5)
    args = ['ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO', 'AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT',
            'CYTO_CXCL10_SWAP_PDL1_SERPINB9', 'ER_NONE', 'NICHE_TRIHORMONAL_120UM_MVF_CAO2',
            'SAF_AAVS1_CAG_ICASP9_BIALLELIC', {'INS': {'n_empirical_studies': 1}}]
    with pytest.raises(ValueError, match='discount must be nonnegative'):
        rank.compute_series_system_utility(*args)


@pytest.mark.parametrize('a,b', [([], []), ([1], []), ([True], [1]), ([0], [1]), ([1.], [1])])
def test_rank_precision_requires_valid_paired_draws(a, b):
    with pytest.raises(ValueError, match='aligned nonempty'):
        rank._rank_precision_summary(a, b, 'synthetic')


def test_rank_precision_positive_control_has_expected_finite_precision():
    result = rank._rank_precision_summary([2] * 100, [4] * 100, 'synthetic')
    interval = result['scenario_a_rank1_frequency_wilson95']
    assert interval == pytest.approx([0, .03699349820698568], abs=2e-15)
    assert result['scenario_a_rank_quantiles_linear']['0.5'] == 2


@pytest.mark.parametrize('value', [2161, 10 ** 309, 10 ** 500])
def test_rank_precision_rejects_outside_grid_before_float_conversion(value):
    with pytest.raises(ValueError, match='declared population'):
        rank._rank_precision_summary([value], [1], 'synthetic')


@pytest.mark.parametrize('max_rank', [0, -1, True, 1.5, '360', 2161, 10 ** 500])
def test_rank_precision_rejects_invalid_population_bound(max_rank):
    with pytest.raises(ValueError, match='Rank population'):
        rank._rank_precision_summary([1], [1], 'synthetic', max_rank=max_rank)


def test_rank_precision_respects_explicit_subset_and_preserves_half_integer_medians():
    values = [1, 2, 3, 4]
    expected = rank._rank_precision_summary(values, list(reversed(values)), 'synthetic')
    actual = rank._rank_precision_summary(values, list(reversed(values)), 'synthetic', max_rank=4)
    assert actual == expected
    assert actual['scenario_a_rank_quantiles_linear']['0.5'] == 2.5
    with pytest.raises(ValueError, match='declared population'):
        rank._rank_precision_summary(values, list(reversed(values)), 'synthetic', max_rank=3)


def test_rank_precision_default_bound_tracks_current_grid_dimensions(monkeypatch):
    for name in ('ALLO_CLOAK_MODULES', 'AUTOANTIGEN_MODULES', 'CYTOKINE_GRANZYME_MODULES', 'ER_AMYLOID_MODULES', 'NICHE_VASCULAR_MODULES', 'SAFETY_SWITCH_MODULES'):
        modules = getattr(rank, name)
        first_key = next(iter(modules))
        monkeypatch.setattr(rank, name, {first_key: modules[first_key]})
    assert rank._rank_precision_summary([1], [1], 'synthetic')['scenario_a_rank_quantiles_linear']['0.5'] == 1
    with pytest.raises(ValueError, match='declared population'):
        rank._rank_precision_summary([2], [1], 'synthetic')


def test_sensitivity_rejects_duplicate_grid_before_estimating_variances():
    profiles = {'INS': {'n_empirical_studies': 1}}
    grid = rank.run_combinatorial_multi_gate_search(profiles)
    grid[-1] = copy.deepcopy(grid[0])
    registry = copy.deepcopy(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY)
    with pytest.raises(ValueError, match='each evaluated grid combination exactly once'):
        rank.compute_direct_output_logit_sensitivity_analysis(profiles, grid)
    assert rank.COMPLETE_GATE_COEFFICIENT_REGISTRY == registry
