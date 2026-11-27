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


def test_iip_0030_curated_corpus_carrier_o2_hpsc_sources_and_atomic_split():
    """Regression & adversarial controls for IIP-0030 (DEC-109) Ring-1 corpus repairs."""
    # 1. Schema enum count (11 -> 12) with hPSC_Derived_Unspecified_Line and additionalProperties: false
    schema = json.loads((ROOT / 'schemas/extraction_schema_v0.2.json').read_text(encoding='utf-8'))
    cs_enum = schema['properties']['arms']['items']['properties']['gate1_substrate_and_differentiation']['properties']['cell_source_class']['enum']
    assert len(cs_enum) == 15 and 'hPSC_Derived_Unspecified_Line' in cs_enum and 'None_Acellular_Device_or_Scaffold' in cs_enum
    assert schema['properties']['arms']['items']['properties']['gate1_substrate_and_differentiation']['additionalProperties'] is False
    assert schema['properties']['arms']['items']['properties']['gate5_engraftment_niche']['additionalProperties'] is False
    assert schema['properties']['atomic_effects']['items']['additionalProperties'] is False

    # 2. Full 48-case synthetic adversarial carrier precision battery (30 TP, 18 TN)
    tp_cases = [
        ('hb_hemoxcell_commercial', 'Islets were cultured with 1 g/L HEMOXCell to maintain oxygenation.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_m101_standalone', 'Supplementing islet preservation medium with M101 (1 g/L) improved islet viability.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_hemoxcell_m101', 'Islet preservation with HEMOXCell (M101) at 5 g/L improved graft survival.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_arenicola_m101', 'Arenicola marina M101 oxygen delivery preserved human islets for 24 h.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_vesicles_abbrev', 'Islets co-encapsulated with Hb-vesicles (HbV) survived 72h hypoxia.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_hboc_201', 'Perfusion with HBOC-201 maintained pO2 above 30 mmHg.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_hboc_phrase', 'Islets were supplemented with hemoglobin-based oxygen carrier (HBOC) in culture.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_oxyglobin', 'Oxyglobin supplementation improved marginal islet survival.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_alginate_microcapsule', 'Encapsulation in hemoglobin-loaded alginate prevented central necrosis.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_myoglobin_peg_gel', 'Myoglobin-conjugated PEG-4MAL gel buffered pericellular oxygen.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_crosslinked_hb', 'Infusion of cross-linked hemoglobin improved graft oxygenation.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_loaded_hydrogel', 'Islets were encapsulated in a hemoglobin-loaded hydrogel oxygen carrier for 72h.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('hb_m101_arenicola_hb', 'Supplementation with marine worm Arenicola marina hemoglobin (M101) preserved islets.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('pfc_oxygen_preceding', 'Hydrogel loaded with 5 vol% oxygen-carrying PFC prevented hypoxia.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_alginate_preceding', 'Alginate encapsulation with 10% emulsified PFC improved islet pO2.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_reversed_order', 'An oxygen-carrying emulsion containing 10% PFC improved islet oxygenation.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_perfluorodecalin', 'Islets cultured over a perfluorodecalin liquid layer showed higher OCR.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_pfob', 'PFOB nanoemulsions increased oxygen solubility in the scaffold.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_oxyfluor', 'Oxyfluor emulsion was added to the culture medium during islet shipment.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_emulsion', 'A perfluorocarbon emulsion oxygen carrier was added to the hydrogel.', 'Perfluorocarbon_O2_Carrier'),
        ('pfc_droplets', 'Loaded 5 vol% PFC droplets enhance oxygen transport in the slab.', 'Perfluorocarbon_O2_Carrier'),
        ('exog_o2_macrodevice', 'Wireless battery-free oxygen generating O2-Macrodevice enables immunosuppression-free islet transplantation.', 'Exogenous_O2_Chamber_or_Bioelectronic'),
        ('exog_beta_air_unicode', 'Transplantation of macroencapsulated human islets within the bioartificial pancreas \u03b2Air to patients with type 1 diabetes.', 'Exogenous_O2_Chamber_or_Bioelectronic'),
        ('tp_contrast_pfc_with_hemoxcell', 'Unlike prior perfluorocarbon emulsions, HEMOXCell was administered to preserve islet grafts.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('tp_pfd_without_hb', 'Islets were encapsulated in oxygen-loaded perfluorodecalin droplets, without hemoglobin.', 'Perfluorocarbon_O2_Carrier'),
        ('tp_cell_free_hb_perfusate', 'Cold pancreas perfusion with cell-free hemoglobin perfusate improved islet oxygenation.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('tp_peg_hemoglobin', 'Co-encapsulation in PEG-hemoglobin hydrogels sustained islet viability.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('tp_pegylated_hb', 'Islets supplemented with pegylated hemoglobin (PEG-Hb) survived hypoxia.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('tp_british_polymerised_haemoglobin', 'Cross-linked polymerised haemoglobin improved oxygenation during islet isolation.', 'Hemoglobin_or_Myoglobin_Hydrogel'),
        ('tp_lowercase_pfc_emulsion', 'Adding 10% pfc emulsion to the encapsulation hydrogel reduced hypoxic death.', 'Perfluorocarbon_O2_Carrier'),
    ]
    tn_cases = [
        ('tn_prefrontal_cortex_oxygen', 'Metabolic imaging of the PFC revealed altered oxygen consumption in T1D patients.'),
        ('tn_plaque_forming_cell_hydrogel', 'Splenic PFC response was measured after hydrogel implantation.'),
        ('tn_negated_pfc_emulsion', 'Unlike prior perfluorocarbon emulsions, our hydrogel-free graft uses no oxygen carrier.'),
        ('tn_negated_pfc_scaffold', 'Unlike prior perfluorocarbon scaffolds, our hydrogel-free approach used no oxygen carrier.'),
        ('tn_hba1c', 'Accompanied by a decrease in glycated hemoglobin (HbA1c) to 5%.'),
        ('tn_stat_pfc', 'Data were analyzed by two-sided t-test (pRA) and Kruskal-Wallis test (pFC).'),
        ('tn_hbo2_oximetry_vesicle', 'Measured HbO2 saturation near extracellular vesicle membranes.'),
        ('tn_clone_m101_hydrogel', 'Clone M101 hydrogel control without hemoglobin.'),
        ('tn_perfluoro_elastomer', 'Compressible perfluoro elastomer ferrules were used to seal PTFE tubing.'),
        ('tn_negated_hboc', 'Unlike prior HBOC-201 perfusion scaffolds, our device operates without hemoglobin.'),
        ('tn_no_hemoxcell_administered', 'No HEMOXCell was administered to any graft in this trial.'),
        ('tn_negated_cao2_loading', 'Control scaffolds were not loaded with calcium peroxide (CaO2-free PDMS).'),
        ('tn_no_cao2_added', 'No CaO2 was added to the baseline alginate device.'),
        ('tn_mri_tracer_pfc', 'Noninvasive tracking used a perfluorocarbon 19F MRI tracer to image graft location.'),
        ('tn_bare_pfc_imaging_contrast', 'Perfluorocarbon nanoparticles were evaluated strictly as a contrast agent for MRI imaging.'),
        ('tn_ultrasound_pfc_microbubble', 'Perfluorocarbon microbubble ultrasound contrast enabled vascular flow quantification.'),
        ('tn_hb_omitted', 'HEMOXCell was omitted from all control islets.'),
        ('tn_hba1c_with_islet_viability', 'Patients showed improved HbA1c (glycated hemoglobin) alongside graft viability.'),
    ]
    assert len(tp_cases) == 30 and len(tn_cases) == 18
    for cid, text, expected in tp_cases:
        got, _ = ext.classify_oxygen_bridge_carrier_hardened(text)
        assert got == expected, f'{cid}: got {got}, expected {expected}'
    for cid, text in tn_cases:
        got, _ = ext.classify_oxygen_bridge_carrier_hardened(text)
        assert got == 'None', f'{cid}: got {got}, expected None'

    # 3. Comparator-arm exogenous O2 inheritance guard
    assert ext._infer_comparator_oxygen_carrier('Exogenous_O2_Chamber_or_Bioelectronic', 'Unoxygenated control macrodevice') == 'None'
    assert ext._infer_comparator_oxygen_carrier('Exogenous_O2_Chamber_or_Bioelectronic', 'Reference control arm (device details unspecified)') is None
    assert ext._infer_comparator_oxygen_carrier('None', 'Vehicle control arm') == 'None'

    # 4. Live extracted corpus assertions across repaired records
    ext_dir = ROOT / 'data/extracted'
    by_pmid = {p.stem: json.loads(p.read_text(encoding='utf-8')) for p in ext_dir.glob('*.json') if p.stem.isdigit()}
    assert len(by_pmid) == 604
    carrier_arms = [
        (pmid, a['gate5_engraftment_niche']['oxygen_bridge_72h_carrier'])
        for pmid, doc in by_pmid.items()
        for a in doc['arms']
        if a['gate5_engraftment_niche']['oxygen_bridge_72h_carrier'] not in (None, 'None')
    ]
    assert sorted(carrier_arms) == [
        ('29288549', 'Exogenous_O2_Chamber_or_Bioelectronic'),
        ('40789835', 'Exogenous_O2_Chamber_or_Bioelectronic'),
        ('41595986', 'Perfluorocarbon_O2_Carrier'),
        ('42182975', 'Exogenous_O2_Chamber_or_Bioelectronic'),
    ]

    # Exact-line & unspecified-line stem-cell repairs (15 total: 7 exact + 8 unspecified)
    expected_exact = {
        '34345868': ('hESC', 'HUES8', 'InVivo_Engrafted_Maturation'),
        '36497105': ('hESC', 'MEL1', 'Stage4_PancreaticProgenitor_PDX1_NKX6_1'),
        '36624157': ('CiPSC_chemically_induced', 'hCiPSC-8# (hPSC-8)', 'InVivo_Engrafted_Maturation'),
        '37040771': ('hESC', 'H1/WA01', None),
        '37809993': ('iPSC_healthy_donor', 'TkDN4-M', 'Stage6_Immature_SC_Islet'),
        '40544428': ('hESC', 'Undisclosed proprietary hESC line (zimislecel / VX-880; USAN-documented hESC origin)', 'InVivo_Engrafted_Maturation'),
        '40700651': ('hESC', 'WA01/H1', 'Stage6_Immature_SC_Islet'),
    }
    for pmid, (exp_src, exp_line, exp_stg) in expected_exact.items():
        g1 = by_pmid[pmid]['arms'][0]['gate1_substrate_and_differentiation']
        assert g1['cell_source_class'] == exp_src, pmid
        assert g1['cell_line_id'] == exp_line, pmid
        assert g1['canonical_differentiation_stage_evaluated'] == exp_stg, pmid

    expected_unspec = ['31756031', '34861146', '36111962', '39211191', '40173261', '40987274', '41932324', '42244553']
    for pmid in expected_unspec:
        assert by_pmid[pmid]['arms'][0]['gate1_substrate_and_differentiation']['cell_source_class'] == 'hPSC_Derived_Unspecified_Line', pmid

    # Liang 36624157 atomic split (eff2 C-peptide 2.0 ng/mL + eff3 HbA1c -44%) and 3,094 total effects
    liang = by_pmid['36624157']
    assert liang['arms'][0]['gate5_engraftment_niche']['transplant_site'] == 'Intramuscular_Forearm_or_Rectus'
    liang_effs = {e['effect_id']: e for e in liang['atomic_effects']}
    assert set(liang_effs) == {'pmid36624157_eff1', 'pmid36624157_eff2', 'pmid36624157_eff3'}
    assert liang_effs['pmid36624157_eff2']['treatment_mean'] == 2.0
    assert liang_effs['pmid36624157_eff2']['measurement_unit'] == 'ng/mL'
    assert liang_effs['pmid36624157_eff2']['log_response_ratio_or_smd'] is None
    assert liang_effs['pmid36624157_eff2']['endpoint_polarity_registry']['orientation_sign'] == 1
    assert liang_effs['pmid36624157_eff3']['treatment_mean'] == 44.0
    assert liang_effs['pmid36624157_eff3']['measurement_unit'] == '%_relative_decrease'
    assert liang_effs['pmid36624157_eff3']['log_response_ratio_or_smd'] is None
    assert liang_effs['pmid36624157_eff3']['endpoint_polarity_registry']['orientation_sign'] == -1

    # Wang 39326417 & Reichman 40544428 curated fields (including mortality domain repair on 40544428 eff2)
    wang = by_pmid['39326417']
    assert wang['study_metadata']['nct_id'] == 'ChiCTR2300072200'
    assert wang['causal_mechanism_edges'] == []
    assert wang['arms'][0]['gate1_substrate_and_differentiation']['canonical_differentiation_stage_evaluated'] == 'InVivo_Engrafted_Maturation'
    reichman = by_pmid['40544428']
    assert reichman['study_metadata']['nct_id'] == 'NCT04786262'
    r_eff2 = next(e for e in reichman['atomic_effects'] if e['effect_id'] == 'pmid40544428_eff2')
    assert r_eff2['outcome_domain'] != 'Differentiation_Yield_and_Lineage_Purity'
    assert r_eff2['contrast_metadata']['endpoint_name'] == 'Clinical_Safety_and_Mortality_Adverse_Events'
    assert sum(len(d['atomic_effects']) for d in by_pmid.values()) == 3174
    # Box-2 Round-1 & v12 Ring-1 comparator/timepoint/stage/cohort repairs & clinical validator negative controls
    assert by_pmid['38993891']['atomic_effects'][0]['contrast_metadata']['timepoint'] is not None
    assert '<=' not in by_pmid['38993891']['atomic_effects'][0]['outcome_variable_Y']
    assert len(by_pmid['30710150']['arms']) == 4 and len(by_pmid['30710150']['atomic_effects']) == 1
    assert by_pmid['30710150']['arms'][3]['gate1_substrate_and_differentiation']['canonical_differentiation_stage_evaluated'] == 'Primary_Mature_Islet'
    assert by_pmid['30710150']['atomic_effects'][0]['comparator_arm_id'] in ('pmid30710150_arm4_adult_human_islets', 'pmid30710150_arm4')
    assert by_pmid['30710150']['atomic_effects'][0]['effect_direction'] == 'Descriptive_or_Unoriented'
    assert len(by_pmid['24222351']['arms']) == 5 and len(by_pmid['24222351']['atomic_effects']) == 3
    assert by_pmid['24222351']['atomic_effects'][0]['comparator_arm_id'] is None
    assert by_pmid['24222351']['atomic_effects'][1]['comparator_arm_id'] is not None
    assert by_pmid['24222351']['atomic_effects'][2]['comparator_arm_id'] is not None
    assert by_pmid['24222351']['atomic_effects'][1]['comparator_arm_id'] != by_pmid['24222351']['atomic_effects'][2]['comparator_arm_id']
    assert len(by_pmid['29862647']['arms']) == 3 and len(by_pmid['29862647']['atomic_effects']) == 2
    assert by_pmid['29862647']['atomic_effects'][0]['comparator_arm_id'] == by_pmid['29862647']['arms'][1]['arm_id']
    assert by_pmid['29862647']['atomic_effects'][1]['comparator_arm_id'] == by_pmid['29862647']['arms'][1]['arm_id']
    assert by_pmid['34504013']['atomic_effects'][1]['effect_id'] == 'pmid34504013_eff2_dual_eqc_dextran_tracer_release'
    assert 'cell-free' not in json.dumps(by_pmid['10531307']).lower() and 'cell-free' not in json.dumps(by_pmid['16177096']).lower()
    assert len(by_pmid['10531307']['arms']) == 4 and len(by_pmid['10531307']['atomic_effects']) == 2
    assert by_pmid['10531307']['atomic_effects'][0]['comparator_arm_id'] != by_pmid['10531307']['atomic_effects'][1]['comparator_arm_id']
    assert all(a['gate5_engraftment_niche'].get('clinical_conditioning_regimen_class') == 'Lymphodepleting_or_Myeloablative_Conditioning' for a in by_pmid['40812299']['arms'])
    assert len(by_pmid['41595986']['arms']) == 3
    assert by_pmid['41595986']['atomic_effects'][0]['source_figure_or_table'] == 'Abstract and Figs. 3, 5, and 6'
    assert by_pmid['41595986']['atomic_effects'][1]['comparator_arm_id'] is None
    assert by_pmid['41595986']['atomic_effects'][2]['comparator_arm_id'] is not None
    assert by_pmid['41595986']['atomic_effects'][0]['comparator_arm_id'] != by_pmid['41595986']['atomic_effects'][2]['comparator_arm_id']
    assert sum(len(e['verbatim_evidence_quote']) for e in by_pmid['41595986']['atomic_effects']) <= 220
    assert '[REDACTED_' not in json.dumps(by_pmid['41595986'])
    assert by_pmid['41595986']['atomic_effects'][2]['effect_direction'] == 'Protective_or_Improved'
    assert by_pmid['41595986']['atomic_effects'][2]['source_figure_or_table'] == 'Fig. 8'
    assert by_pmid['41595986']['atomic_effects'][2]['source_span'].startswith('fulltext_data_sections.results_paragraphs[27]')
    assert '24 h' in (by_pmid['41595986']['atomic_effects'][2]['contrast_metadata']['timepoint'] or '')
    assert by_pmid['42627334']['atomic_effects'][0]['statistic_type'] == 'log_response_ratio_paired_means'
    assert by_pmid['42627334']['atomic_effects'][0]['treatment_mean'] == 0.43 and by_pmid['42627334']['atomic_effects'][0]['control_mean'] == 0.35
    assert by_pmid['41543848']['atomic_effects'][1]['treatment_mean'] == 67.27 and by_pmid['41543848']['atomic_effects'][1]['control_mean'] == 44.95
    assert by_pmid['41543848']['atomic_effects'][4]['treatment_mean'] == 80.73 and by_pmid['41543848']['atomic_effects'][4]['control_mean'] == 66.62
    assert by_pmid['40320524']['arms'][0]['gate2_3_4_genetic_and_epitope_edits'][0]['edit_modality'] == 'CRISPR_Cas12a_Multiplex_KO'
    assert by_pmid['40320524']['arms'][0]['gate5_engraftment_niche']['transplant_site'] == 'Kidney_Capsule'
    assert len(by_pmid['21813778']['arms']) == 6 and len(by_pmid['21813778']['atomic_effects']) == 3
    assert by_pmid['18021279']['arms'][0]['sample_size_unit'] == 'human_islet_preparations'
    assert by_pmid['18021279']['atomic_effects'][1]['n_unit'] == 'human_islet_preparations'
    assert by_pmid['16177096']['atomic_effects'][0]['comparator_arm_id'] != by_pmid['16177096']['atomic_effects'][1]['comparator_arm_id']
    assert by_pmid['41870990']['arms'][0]['gate1_substrate_and_differentiation']['canonical_differentiation_stage_evaluated'] is None
    # v13 Bug #3: public-safe synthetic-fixture validator negative controls (runs in BOTH OI and public_export pytest!)
    import copy as _c_r3_s, importlib.util as _ilu_r3_s, pytest as _pt_r3_s
    _sp_ext_s = _ilu_r3_s.spec_from_file_location('openislet_ext_r3_pub_check', ROOT / 'scripts' / 'extract_ring1_corpus.py')
    _mod_ext_s = _ilu_r3_s.module_from_spec(_sp_ext_s)
    _sp_ext_s.loader.exec_module(_mod_ext_s)
    _syn_q_mathieu = 'At 12 months, the mean ln(AUC C-peptide + 1) in the 0\u00b75 mg/kg ATG group was 0\u00b7513 nmol/L per min (SD 0\u00b7032), with a mean baseline-adjusted difference from placebo of 0\u00b7102 nmol/L per min (95% CI 0\u00b7021-0\u00b7183; p=0\u00b7014).'
    _syn_dos_m = {'pmid': '40976248', 'title': 'MELD-ATG', 'abstract': _syn_q_mathieu, 'fulltext_data_sections': {}}
    _syn_doc_m = _c_r3_s.deepcopy(by_pmid['40976248'])
    for _a_sm in _syn_doc_m['arms']:
        _a_sm['sample_size_n'] = None
    _syn_doc_m['atomic_effects'] = [_syn_doc_m['atomic_effects'][0]]
    _syn_doc_m['atomic_effects'][0]['verbatim_evidence_quote'] = _syn_q_mathieu
    _syn_doc_m['atomic_effects'][0]['source_span'] = f'abstract[0:{len(_syn_q_mathieu)}]|norm=exact'
    _syn_doc_m['atomic_effects'][0]['dossier_sha256'] = _mod_ext_s._compute_dossier_sha256('40976248', _syn_dos_m)
    _mod_ext_s.validate_study_semantics(_syn_doc_m, dossier=_syn_dos_m)
    _syn_m_p_ci = _c_r3_s.deepcopy(_syn_doc_m)
    _syn_m_p_ci['atomic_effects'][0]['ci_95_lower'] = 0.014
    with _pt_r3_s.raises(ValueError, match='p-value clause|parsed CI span'):
        _mod_ext_s.validate_study_semantics(_syn_m_p_ci, dossier=_syn_dos_m)
    _syn_q_moran = 'between the anakinra and the placebo groups at 9 months was 0\u00b702 nmol/L (-0\u00b709 to 0\u00b715; p=0\u00b771).'
    _syn_dos_mo = {'pmid': '23562090', 'title': 'AIDA', 'abstract': _syn_q_moran, 'fulltext_data_sections': {}}
    _syn_doc_mo = _c_r3_s.deepcopy(by_pmid['23562090'])
    for _a_smo in _syn_doc_mo['arms']:
        _a_smo['sample_size_n'] = None
    _syn_doc_mo['atomic_effects'] = [_syn_doc_mo['atomic_effects'][0]]
    _syn_doc_mo['atomic_effects'][0]['verbatim_evidence_quote'] = _syn_q_moran
    _syn_doc_mo['atomic_effects'][0]['source_span'] = f'abstract[0:{len(_syn_q_moran)}]|norm=exact'
    _syn_doc_mo['atomic_effects'][0]['dossier_sha256'] = _mod_ext_s._compute_dossier_sha256('23562090', _syn_dos_mo)
    _mod_ext_s.validate_study_semantics(_syn_doc_mo, dossier=_syn_dos_mo)
    _syn_mo_p_ci = _c_r3_s.deepcopy(_syn_doc_mo)
    _syn_mo_p_ci['atomic_effects'][0]['ci_95_upper'] = 0.71
    with _pt_r3_s.raises(ValueError, match='p-value clause|parsed CI span'):
        _mod_ext_s.validate_study_semantics(_syn_mo_p_ci, dossier=_syn_dos_mo)
    _syn_q_wr = 'Clodronate liposome-mediated depletion of islet macrophages improved glucose tolerance and blocked proinflammatory gene expression'
    _syn_dos_wr = {'pmid': '24222351', 'title': 'WR2014', 'abstract': _syn_q_wr, 'fulltext_data_sections': {}}
    _syn_doc_wr = _c_r3_s.deepcopy(by_pmid['24222351'])
    for _a_swr in _syn_doc_wr['arms']:
        _a_swr['sample_size_n'] = None
    _syn_doc_wr['arms'][1]['arm_label'] = 'Clodronate-liposome macrophage-depleted high-fat-fed FVB hIAPP transgenic mice (Fig. 6)'
    _syn_doc_wr['arms'][2]['arm_label'] = 'Control-liposome macrophage-intact high-fat-fed FVB hIAPP transgenic mice comparator (Fig. 6)'
    _syn_doc_wr['arms'][3]['arm_label'] = 'Clodronate-liposome macrophage-depleted obese crossbred hIAPP mice (Fig. 7)'
    _syn_doc_wr['arms'][4]['arm_label'] = 'Control-liposome macrophage-intact obese crossbred hIAPP mice comparator (Fig. 7)'
    _syn_eff_wr = _syn_doc_wr['atomic_effects'][1]
    _syn_eff_wr['outcome_variable_Y'] = '[Endpoint: Fig. 6 high-fat-fed FVB hIAPP mouse glucose tolerance improvement and blockade of islet proinflammatory transcripts after clodronate macrophage depletion]'
    _syn_eff_wr['verbatim_evidence_quote'] = _syn_q_wr
    _syn_eff_wr['source_span'] = f'abstract[0:{len(_syn_q_wr)}]|norm=exact'
    _syn_eff_wr['dossier_sha256'] = _mod_ext_s._compute_dossier_sha256('24222351', _syn_dos_wr)
    _syn_doc_wr['atomic_effects'] = [_syn_eff_wr]
    _mod_ext_s.validate_study_semantics(_syn_doc_wr, dossier=_syn_dos_wr)
    _syn_wr_noc = _c_r3_s.deepcopy(_syn_doc_wr)
    _syn_wr_noc['atomic_effects'][0]['comparator_arm_id'] = None
    with _pt_r3_s.raises(ValueError, match='requires non-null comparator_arm_id'):
        _mod_ext_s.validate_study_semantics(_syn_wr_noc, dossier=_syn_dos_wr)
    _syn_wr_both = _c_r3_s.deepcopy(_syn_doc_wr)
    _syn_wr_both['atomic_effects'][0]['focal_arm_id'] = _syn_wr_both['arms'][3]['arm_id']
    _syn_wr_both['atomic_effects'][0]['comparator_arm_id'] = _syn_wr_both['arms'][4]['arm_id']
    _syn_wr_both['atomic_effects'][0]['contrast_metadata']['intervention_label'] = _syn_wr_both['arms'][3]['arm_label']
    _syn_wr_both['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _syn_wr_both['arms'][4]['arm_label']
    with _pt_r3_s.raises(ValueError, match='cohort_id|figure/table ref|cohort tokens conflict'):
        _mod_ext_s.validate_study_semantics(_syn_wr_both, dossier=_syn_dos_wr)
    # v14 Bug #6: strict interior CI (reject degenerate lower==point==upper) & build_cross_study_graph.py parity on BOTH OI and public_export!
    _sp_gr_s = _ilu_r3_s.spec_from_file_location('openislet_gr_r3_pub_check', ROOT / 'scripts' / 'build_cross_study_graph.py')
    _mod_gr_s = _ilu_r3_s.module_from_spec(_sp_gr_s)
    _sp_gr_s.loader.exec_module(_mod_gr_s)
    _syn_q_deg = 'between active and placebo was 0.10 nmol/L (95% CI 0.10 to 0.10; p=0.014).'
    _syn_dos_deg = {'pmid': '40976248', 'title': 'DEG', 'abstract': _syn_q_deg, 'fulltext_data_sections': {}}
    _syn_doc_deg = _c_r3_s.deepcopy(_syn_doc_m)
    _syn_doc_deg['atomic_effects'][0]['verbatim_evidence_quote'] = _syn_q_deg
    _syn_doc_deg['atomic_effects'][0]['source_span'] = f'abstract[0:{len(_syn_q_deg)}]|norm=exact'
    _syn_doc_deg['atomic_effects'][0]['dossier_sha256'] = _mod_ext_s._compute_dossier_sha256('40976248', _syn_dos_deg)
    _syn_doc_deg['atomic_effects'][0]['point_estimate'] = 0.10
    _syn_doc_deg['atomic_effects'][0]['ci_95_lower'] = 0.10
    _syn_doc_deg['atomic_effects'][0]['ci_95_upper'] = 0.10
    with _pt_r3_s.raises(ValueError, match='degenerate 95% CI|strict interior'):
        _mod_ext_s.validate_study_semantics(_syn_doc_deg, dossier=_syn_dos_deg)
    with _pt_r3_s.raises(ValueError, match='degenerate 95% CI|strict interior'):
        _mod_gr_s._validate_graph_references(_syn_doc_deg, '40976248.json')
    # Both validators reject Mathieu null exact p_value
    _syn_m_nullp = _c_r3_s.deepcopy(_syn_doc_m)
    _syn_m_nullp['atomic_effects'][0]['p_value'] = None
    with _pt_r3_s.raises(ValueError, match='disagrees with p_value_bound|requires a grounded exact p-value'):
        _mod_ext_s.validate_study_semantics(_syn_m_nullp, dossier=_syn_dos_m)
    with _pt_r3_s.raises(ValueError, match='disagrees with p_value_bound|requires a grounded exact p-value'):
        _mod_gr_s._validate_graph_references(_syn_m_nullp, '40976248.json')
    # Both validators reject Moran arm-pair swap
    _syn_mo_swap = _c_r3_s.deepcopy(_syn_doc_mo)
    _syn_mo_swap['atomic_effects'][0]['focal_arm_id'] = _syn_mo_swap['arms'][2]['arm_id']
    _syn_mo_swap['atomic_effects'][0]['comparator_arm_id'] = _syn_mo_swap['arms'][3]['arm_id']
    _syn_mo_swap['atomic_effects'][0]['contrast_metadata']['intervention_label'] = _syn_mo_swap['arms'][2]['arm_label']
    _syn_mo_swap['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _syn_mo_swap['arms'][3]['arm_label']
    with _pt_r3_s.raises(ValueError, match='cohort_id|conflicts with focal_arm_id|cohort tokens conflict'):
        _mod_ext_s.validate_study_semantics(_syn_mo_swap, dossier=_syn_dos_mo)
    with _pt_r3_s.raises(ValueError, match='cohort_id|conflicts with focal_arm_id|cohort tokens conflict'):
        _mod_gr_s._validate_graph_references(_syn_mo_swap, '23562090.json')
    # Both validators reject Westwell-Roper Fig. 6 -> Fig. 7 arm-pair swap
    with _pt_r3_s.raises(ValueError, match='cohort_id|figure/table ref|cohort tokens conflict'):
        _mod_gr_s._validate_graph_references(_syn_wr_both, '24222351.json')
    # Public/source-free semantic regression suite against LIVE shipped records (executes cleanly even when dossiers are absent!)
    # 1. Westwell-Roper 21813778:eff2: switch Fig. 5E comparator arm4 -> rodent-IAPP arm6 (even with comparator_label updated)
    _pub_wr5e = _c_r3_s.deepcopy(by_pmid['21813778'])
    _pub_wr5e['atomic_effects'][1]['comparator_arm_id'] = _pub_wr5e['arms'][5]['arm_id']
    _pub_wr5e['atomic_effects'][1]['contrast_metadata']['comparator_label'] = _pub_wr5e['arms'][5]['arm_label']
    with _pt_r3_s.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
        _mod_gr_s._validate_graph_references(_pub_wr5e, '21813778.json')
    # 2. Mathieu 40976248:eff1: replace CI lower 0.021 with p-value 0.014; remove exact p_value; collapse CI to point_estimate
    _pub_math_ci = _c_r3_s.deepcopy(by_pmid['40976248'])
    _pub_math_ci['atomic_effects'][0]['ci_95_lower'] = 0.014
    with _pt_r3_s.raises(ValueError, match='p-value clause|parsed CI span'):
        _mod_gr_s._validate_graph_references(_pub_math_ci, '40976248.json')
    _pub_math_nullp = _c_r3_s.deepcopy(by_pmid['40976248'])
    _pub_math_nullp['atomic_effects'][0]['p_value'] = None
    with _pt_r3_s.raises(ValueError, match='disagrees with p_value_bound|requires a grounded exact p-value'):
        _mod_gr_s._validate_graph_references(_pub_math_nullp, '40976248.json')
    _pub_math_col = _c_r3_s.deepcopy(by_pmid['40976248'])
    _pub_math_col['atomic_effects'][0]['ci_95_lower'] = 0.102
    _pub_math_col['atomic_effects'][0]['ci_95_upper'] = 0.102
    with _pt_r3_s.raises(ValueError, match='degenerate 95% CI|strict interior'):
        _mod_gr_s._validate_graph_references(_pub_math_col, '40976248.json')
    # 3. BANDIT 42627334:eff1: change exact p_value 0.336 -> 0.995 leaving bound & quote at 0.336
    _pub_bandit_p = _c_r3_s.deepcopy(by_pmid['42627334'])
    _pub_bandit_p['atomic_effects'][0]['p_value'] = 0.995
    with _pt_r3_s.raises(ValueError, match='disagrees with p_value_bound|does not match quote token'):
        _mod_gr_s._validate_graph_references(_pub_bandit_p, '42627334.json')
    # v12 Bug #1: rank_islet1_stacks must yield 0 target-profile diffs even if 40812299 free-text regimen is redacted (source-free pub2pub)
    import copy as _c_r3_p, importlib.util as _ilu_r3_p
    _sp_rnk = _ilu_r3_p.spec_from_file_location('openislet_rnk_r3_check', ROOT / 'scripts' / 'rank_islet1_stacks.py')
    _mod_rnk = _ilu_r3_p.module_from_spec(_sp_rnk)
    _sp_rnk.loader.exec_module(_mod_rnk)
    _, _rnk_studies, _rnk_omics = _mod_rnk.load_inputs()
    _prof_base = _mod_rnk.build_target_evidence_profiles(_rnk_studies, _rnk_omics)
    _rnk_redacted = _c_r3_p.deepcopy(_rnk_studies)
    for _sd_r in _rnk_redacted:
        if str((_sd_r.get('study_metadata') or {}).get('pmid')) == '40812299':
            for _a_r in _sd_r.get('arms', []):
                _a_r['gate5_engraftment_niche']['systemic_immunosuppression_regimen'] = '[REDACTED_NON_OA_EXCERPT: see PMID:40812299 field=omitted_text quote_sha256=' + 'a' * 64 + ']'
    _prof_redacted = _mod_rnk.build_target_evidence_profiles(_rnk_redacted, _rnk_omics)
    assert _prof_base == _prof_redacted, f'Redacting 40812299 regimen mutated profiles: {[k for k in _prof_base if _prof_base[k] != _prof_redacted.get(k)]}'
    if not (ROOT / '.openislet_public_bundle').is_file():
        import copy as _c_r3, importlib.util as _ilu_r3, json as _js_r3, pytest as _pt_r3, re as _re_r3
        assert 'IL1RN antagonist' not in by_pmid['21813778']['arms'][0]['arm_label']
        assert 'TNF-alpha' in by_pmid['21813778']['atomic_effects'][1]['outcome_variable_Y']
        assert 'LHR-A' in by_pmid['10531307']['arms'][0]['arm_label']
        assert 'Table V' in by_pmid['10531307']['arms'][2]['arm_label'] and 'Table VI' not in by_pmid['10531307']['arms'][2]['arm_label']
        assert '2006-2009' in by_pmid['29862647']['arms'][1]['arm_label'] and '2011-2012' not in _js_r3.dumps(by_pmid['29862647'])
        assert 'not formal Stage-7' not in _js_r3.dumps(by_pmid['41870990']) and '7-stage protocol' in _js_r3.dumps(by_pmid['41870990'])
        assert '13 oncology' in by_pmid['40812299']['atomic_effects'][0]['contrast_metadata']['population_or_model'] and '1 autoimmune' in by_pmid['40812299']['atomic_effects'][0]['contrast_metadata']['population_or_model']
        assert 'non-significant superiority' in (by_pmid['42627334']['atomic_effects'][0].get('effect_notes') or '')
        _sp_ext = _ilu_r3.spec_from_file_location('openislet_ext_r3_check', ROOT / 'scripts' / 'extract_ring1_corpus.py')
        _mod_ext = _ilu_r3.module_from_spec(_sp_ext)
        _sp_ext.loader.exec_module(_mod_ext)
        _moran = by_pmid['23562090']
        _mod_ext.validate_study_semantics(_moran)
        _m_fab = _c_r3.deepcopy(_moran)
        _m_fab['atomic_effects'][0]['point_estimate'] = 9.8
        with _pt_r3.raises(ValueError, match='outside|not grounded'):
            _mod_ext.validate_study_semantics(_m_fab)
        _m_rev = _c_r3.deepcopy(_moran)
        _m_rev['atomic_effects'][0]['ci_95_lower'] = 0.15
        _m_rev['atomic_effects'][0]['ci_95_upper'] = -0.09
        with _pt_r3.raises(ValueError, match='reversed 95% CI'):
            _mod_ext.validate_study_semantics(_m_rev)
        _m_noc = _c_r3.deepcopy(_moran)
        _m_noc['atomic_effects'][0]['comparator_arm_id'] = None
        with _pt_r3.raises(ValueError, match='non-null comparator_arm_id'):
            _mod_ext.validate_study_semantics(_m_noc)
        _m_swp_arm = _c_r3.deepcopy(_moran)
        _m_swp_arm['atomic_effects'][0]['comparator_arm_id'] = 'pmid23562090_arm4_cana_placebo'
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts with contrast_metadata|cross-trial|cross-wired'):
            _mod_ext.validate_study_semantics(_m_swp_arm)
        _m_swp_edit_lbl = _c_r3.deepcopy(_moran)
        _m_swp_edit_lbl['atomic_effects'][0]['comparator_arm_id'] = 'pmid23562090_arm4_cana_placebo'
        _m_swp_edit_lbl['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _m_swp_edit_lbl['arms'][3]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|cross-trial|cross-wired|conflicts'):
            _mod_ext.validate_study_semantics(_m_swp_edit_lbl)
        _m_swp_no_nct = _c_r3.deepcopy(_m_swp_edit_lbl)
        for _a_m in _m_swp_no_nct['arms']:
            _a_m['arm_label'] = _re_r3.sub(r'\s*\(NCT\d+\)', '', _a_m['arm_label'])
        _m_swp_no_nct['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _m_swp_no_nct['arms'][3]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|cross-wired|conflicts'):
            _mod_ext.validate_study_semantics(_m_swp_no_nct)
        # v12 Bug #4: Moran whole-pair reassignment (anakinra contrast -> canakinumab arm3 + arm4)
        _m_whole = _c_r3.deepcopy(_moran)
        _m_whole['atomic_effects'][0]['focal_arm_id'] = 'pmid23562090_arm3_canakinumab'
        _m_whole['atomic_effects'][0]['comparator_arm_id'] = 'pmid23562090_arm4_cana_placebo'
        _m_whole['atomic_effects'][0]['contrast_metadata']['intervention_label'] = _m_whole['arms'][2]['arm_label']
        _m_whole['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _m_whole['arms'][3]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts with focal_arm_id'):
            _mod_ext.validate_study_semantics(_m_whole)
        _m_swp_p = _c_r3.deepcopy(_moran)
        _m_swp_p['atomic_effects'][0]['p_value'] = 0.86
        with _pt_r3.raises(ValueError, match='not grounded|disagrees with p_value_bound'):
            _mod_ext.validate_study_semantics(_m_swp_p)
        _mathieu = by_pmid['40976248']
        _mathieu_drop_p = _c_r3.deepcopy(_mathieu)
        _mathieu_drop_p['atomic_effects'][0]['p_value'] = None
        with _pt_r3.raises(ValueError, match='disagrees with p_value_bound|requires a grounded exact p-value'):
            _mod_ext.validate_study_semantics(_mathieu_drop_p)
        # v12 Bug #4: Mathieu CI lower bound substituted for point_estimate
        _mathieu_ci_pe = _c_r3.deepcopy(_mathieu)
        _mathieu_ci_pe['atomic_effects'][0]['point_estimate'] = 0.021
        with _pt_r3.raises(ValueError, match='equals a 95% CI endpoint'):
            _mod_ext.validate_study_semantics(_mathieu_ci_pe)
        # v12 Bug #4: Synthetic inequality-only p-value with null exact p_value on adjusted_mean_difference_with_95ci_and_pvalue
        _mathieu_ineq = _c_r3.deepcopy(_mathieu)
        _mathieu_ineq_dos = _js_r3.loads((ROOT / 'data' / 'dossiers' / '40976248.json').read_text(encoding='utf-8'))
        _mathieu_ineq_dos['abstract'] = _mathieu_ineq_dos['abstract'].replace('p=0\u00b7014', 'p<0\u00b705')
        _mathieu_ineq['atomic_effects'][0]['verbatim_evidence_quote'] = _mathieu_ineq['atomic_effects'][0]['verbatim_evidence_quote'].replace('p=0\u00b7014', 'p<0\u00b705')
        _mathieu_ineq['atomic_effects'][0]['source_span'] = f"abstract[{_mathieu_ineq_dos['abstract'].index(_mathieu_ineq['atomic_effects'][0]['verbatim_evidence_quote'])}:{_mathieu_ineq_dos['abstract'].index(_mathieu_ineq['atomic_effects'][0]['verbatim_evidence_quote']) + len(_mathieu_ineq['atomic_effects'][0]['verbatim_evidence_quote'])}]|norm=exact"
        _mathieu_ineq['atomic_effects'][0]['dossier_sha256'] = _mod_ext._compute_dossier_sha256('40976248', _mathieu_ineq_dos)
        _mathieu_ineq['atomic_effects'][0]['p_value'] = None
        _mathieu_ineq['atomic_effects'][0]['p_value_relation'] = '<'
        _mathieu_ineq['atomic_effects'][0]['p_value_bound'] = 0.05
        with _pt_r3.raises(ValueError, match='requires a grounded exact p-value'):
            _mod_ext.validate_study_semantics(_mathieu_ineq, dossier=_mathieu_ineq_dos)
        # v12 Bug #4: Gou WT/GTKO comparator swaps & Westwell-Roper Fig. 5E -> rodent-IAPP & Fig. 6 -> Fig. 7 & Naziruddin CTL -> ETA
        _gou_wt = _c_r3.deepcopy(by_pmid['41543848'])
        _gou_wt['atomic_effects'][0]['comparator_arm_id'] = 'pmid41543848_arm4_gtko_paec_nc'
        _gou_wt['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _gou_wt['arms'][3]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
            _mod_ext.validate_study_semantics(_gou_wt)
        _gou_gtko = _c_r3.deepcopy(by_pmid['41543848'])
        _gou_gtko['atomic_effects'][3]['comparator_arm_id'] = 'pmid41543848_arm2_paec_nc'
        _gou_gtko['atomic_effects'][3]['contrast_metadata']['comparator_label'] = _gou_gtko['arms'][1]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
            _mod_ext.validate_study_semantics(_gou_gtko)
        _wr5e = _c_r3.deepcopy(by_pmid['21813778'])
        _wr5e['atomic_effects'][1]['comparator_arm_id'] = 'pmid21813778_arm6_bmdm_rodent_iapp_polypeptide_ctrl'
        _wr5e['atomic_effects'][1]['contrast_metadata']['comparator_label'] = _wr5e['arms'][5]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
            _mod_ext.validate_study_semantics(_wr5e)
        _wr67 = _c_r3.deepcopy(by_pmid['24222351'])
        _wr67['atomic_effects'][1]['comparator_arm_id'] = 'pmid24222351_arm5_control_liposome_obese_crossbred'
        _wr67['atomic_effects'][1]['contrast_metadata']['comparator_label'] = _wr67['arms'][4]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
            _mod_ext.validate_study_semantics(_wr67)
        _naz = _c_r3.deepcopy(by_pmid['29862647'])
        _naz['atomic_effects'][0]['comparator_arm_id'] = 'pmid29862647_arm3_tpiat_eta_only'
        _naz['atomic_effects'][0]['contrast_metadata']['comparator_label'] = _naz['arms'][2]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
            _mod_ext.validate_study_semantics(_naz)
        # v13 Bug #2, #3, #5 full-dossier negative controls & precision assertions
        assert 'unspaced' not in _js_r3.dumps(by_pmid['10531307']).lower() and 'one intervening CCP' in by_pmid['10531307']['arms'][3]['arm_label']
        assert '58.0 +/- 0.2 SD' in by_pmid['16177096']['atomic_effects'][0]['outcome_variable_Y'] and 'Fig. 2A' in by_pmid['16177096']['atomic_effects'][1]['source_figure_or_table']
        _mathieu_p_ci = _c_r3.deepcopy(_mathieu)
        _mathieu_p_ci['atomic_effects'][0]['ci_95_lower'] = 0.014
        with _pt_r3.raises(ValueError, match='p-value clause|parsed CI span'):
            _mod_ext.validate_study_semantics(_mathieu_p_ci)
        _moran_p_ci = _c_r3.deepcopy(_moran)
        _moran_p_ci['atomic_effects'][0]['ci_95_upper'] = 0.71
        with _pt_r3.raises(ValueError, match='p-value clause|parsed CI span'):
            _mod_ext.validate_study_semantics(_moran_p_ci)
        _wr67_noc = _c_r3.deepcopy(by_pmid['24222351'])
        _wr67_noc['atomic_effects'][1]['comparator_arm_id'] = None
        with _pt_r3.raises(ValueError, match='requires non-null comparator_arm_id'):
            _mod_ext.validate_study_semantics(_wr67_noc)
        _wr67_both = _c_r3.deepcopy(by_pmid['24222351'])
        _wr67_both['atomic_effects'][1]['focal_arm_id'] = 'pmid24222351_arm4_clodronate_obese_crossbred'
        _wr67_both['atomic_effects'][1]['comparator_arm_id'] = 'pmid24222351_arm5_control_liposome_obese_crossbred'
        _wr67_both['atomic_effects'][1]['contrast_metadata']['intervention_label'] = _wr67_both['arms'][3]['arm_label']
        _wr67_both['atomic_effects'][1]['contrast_metadata']['comparator_label'] = _wr67_both['arms'][4]['arm_label']
        with _pt_r3.raises(ValueError, match='cohort_id|figure/table ref|cohort tokens conflict'):
            _mod_ext.validate_study_semantics(_wr67_both)
        # v14 Bug #6: full-dossier parity check on build_cross_study_graph._validate_graph_references
        _sp_gr_f = _ilu_r3.spec_from_file_location('openislet_gr_r3_full_check', ROOT / 'scripts' / 'build_cross_study_graph.py')
        _mod_gr_f = _ilu_r3.module_from_spec(_sp_gr_f)
        _sp_gr_f.loader.exec_module(_mod_gr_f)
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts with focal_arm_id|cohort tokens conflict'):
            _mod_gr_f._validate_graph_references(_m_whole, '23562090.json')
        with _pt_r3.raises(ValueError, match='cohort_id|conflicts|cross-wired'):
            _mod_gr_f._validate_graph_references(_wr5e, '21813778.json')
        with _pt_r3.raises(ValueError, match='cohort_id|figure/table ref|cohort tokens conflict'):
            _mod_gr_f._validate_graph_references(_wr67_both, '24222351.json')
        with _pt_r3.raises(ValueError, match='disagrees with p_value_bound|requires a grounded exact p-value'):
            _mod_gr_f._validate_graph_references(_mathieu_drop_p, '40976248.json')
        # Case 3 (BANDIT) & case 4 (Mathieu in-memory dossier tamper with original on-disk digest) in validate_study_semantics
        _bandit_bad_p = _c_r3.deepcopy(by_pmid['42627334'])
        _bandit_bad_p['atomic_effects'][0]['p_value'] = 0.995
        with _pt_r3.raises(ValueError, match='not grounded|disagrees with p_value_bound'):
            _mod_ext.validate_study_semantics(_bandit_bad_p)
        _mathieu_mem_tamper = _c_r3.deepcopy(_mathieu)
        _mathieu_mem_dos = _js_r3.loads((ROOT / 'data' / 'dossiers' / '40976248.json').read_text(encoding='utf-8'))
        _mathieu_mem_dos['abstract'] = _mathieu_mem_dos['abstract'].replace('p=0\u00b7014', 'p=0\u00b7019')
        _mathieu_mem_tamper['atomic_effects'][0]['verbatim_evidence_quote'] = _mathieu_mem_tamper['atomic_effects'][0]['verbatim_evidence_quote'].replace('p=0\u00b7014', 'p=0\u00b7019')
        _mathieu_mem_tamper['atomic_effects'][0]['p_value'] = 0.019
        _mathieu_mem_tamper['atomic_effects'][0]['p_value_bound'] = 0.019
        # Deliberately retain original on-disk dossier_sha256: MUST FAIL because digest must bind exact same bytes as offsets!
        with _pt_r3.raises(ValueError, match='dossier_sha256 does not match declared source'):
            _mod_ext.validate_study_semantics(_mathieu_mem_tamper, dossier=_mathieu_mem_dos)

    # Public export 240-char aggregate verbatim quote ceiling check when running inside an exported bundle
    if (ROOT / '.openislet_public_bundle').is_file():
        verbatim_len = sum(
            len(e['verbatim_evidence_quote'])
            for e in liang['atomic_effects']
            if not e['verbatim_evidence_quote'].startswith(('[REDACTED_NON_OA_EXCERPT:', '[REDACTED_CONSERVATIVE_EXCERPT:'))
        )
        assert verbatim_len <= 240
        assert any(
            e['verbatim_evidence_quote'].startswith('[REDACTED_NON_OA_EXCERPT: see PMID:36624157')
            for e in liang['atomic_effects']
        )

