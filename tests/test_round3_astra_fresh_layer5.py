"""Fresh Stage-0 Layer5 arithmetic and publication counterexamples.

All assay inputs and publication faults below are synthetic. No live artifacts
are regenerated. Thermodynamic oracles use SantaLucia 1998 Eq. 1/Table 1.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import ingest_wetlab_results as w
from scripts import run_layer5_personalized_patient_compiler as l5


def require(condition, message):
    if not condition:
        raise AssertionError(message)


@pytest.mark.parametrize('sequence,expected', [
    ('CGTTGA', -5.35),  # Paper Figure 1: unlike terminal pairs; non-self-complementary.
    ('GCGC', -4.26),   # Two G:C ends and self-complementarity correction.
    ('ATAT', .15),     # Two A:T ends plus self-complementarity.
    ('AAAA', -.94),    # Two A:T ends without self-complementarity.
])
def test_santalucia_terminal_and_symmetry_terms_match_source_oracles(sequence, expected):
    require(l5.compute_duplex_dg37_kcal_mol(sequence) == expected,
            'DNA duplex omitted the terminal-pair or symmetry contribution')
    require(l5.compute_duplex_dg37_kcal_mol(l5.revcomp(sequence)) == expected,
            'DNA duplex energy changed under reverse complement')


@pytest.mark.parametrize('boundary', [2., 6., 100.])
def test_equal_boundary_means_remain_exactly_inside_the_physical_domain(boundary):
    result = w.conjugate_normal_update(boundary, .001, boundary, .003,
                                       min_valid=0., max_valid=boundary)
    require(result['posterior_mean_exact'] == boundary,
            'a convex combination of equal boundary means escaped its domain')
    with localcontext() as context:
        context.prec = 70
        expected_sd = float((1 / (1 / Decimal('.001')**2 + 1 / Decimal('.003')**2)).sqrt())
    require(math.isclose(result['posterior_sd_exact'], expected_sd, rel_tol=2e-15),
            'convex mean repair changed the conjugate precision update')


@pytest.mark.parametrize('left,right', [(1e308, -1e308), (-1e308, 1e308)])
def test_convex_mean_handles_opposite_extreme_finite_endpoints(left, right):
    require(w._normal_normal_moments(left, .001, right, .001)[0] == 0.,
            'stable convex update must not overflow the endpoint difference')


def test_existing_in_hull_display_rounding_tie_can_continue_as_version_two_state():
    # This legal genesis update rounded down under the existing v2 arithmetic.
    # Reordering every weighted sum can move it one ULP and invalidate its
    # display field during continuation, despite unchanged Gaussian moments.
    first = w.run_calibration(dict(batch_id='SYNTHETIC-TIE-A',
                                   obs_stage7_yield_ratio=.4529111328125,
                                   obs_stage7_yield_se=.01), True)
    endpoint = first['conjugate_posterior_updates']['stage7_sc_beta_yield_fraction_arm1b']
    require(endpoint['posterior_mean_exact'] == .50025 and endpoint['posterior_mean'] == .5002,
            'in-hull v2 arithmetic changed at a four-decimal display tie')
    second = w.run_calibration(dict(batch_id='SYNTHETIC-TIE-B',
                                    obs_stage7_yield_ratio=.8,
                                    obs_stage7_yield_se=.02), True, first)
    require(len(second['ingested_batch_ids']) == 2, 'valid v2 continuation failed')


@pytest.fixture
def publication(tmp_path, monkeypatch):
    monkeypatch.setattr(w.tempfile, 'gettempdir', lambda: str(tmp_path))
    canonical, alias = tmp_path/'canonical.json', tmp_path/'alias.json'
    monkeypatch.setattr(w, 'DEFAULT_OUT', canonical)
    monkeypatch.setattr(w, 'LEGACY_ALIAS_PATH', alias)
    report = {'calibration_mode': 'SYNTHETIC_DRY_RUN_VERIFICATION_DO_NOT_CLAIM_AS_WETLAB_DATA'}
    return canonical, alias, report


@pytest.mark.parametrize('existing', [False, True])
@pytest.mark.parametrize('intervening', [b'{"newer_owner":true}\n', b'partial external write'])
def test_alias_failure_preserves_intervening_canonical_bytes(publication, monkeypatch, existing, intervening):
    canonical, alias, report = publication
    old = b'{"previous":true}\n' if existing else None
    if old is not None:
        canonical.write_bytes(old)
        alias.write_bytes(old)
    original = w._atomic_write_json
    def fail_after_external_change(path, data):
        if path == alias:
            canonical.write_bytes(intervening)
            raise OSError('synthetic alias failure')
        return original(path, data)
    monkeypatch.setattr(w, '_atomic_write_json', fail_after_external_change)
    with pytest.raises(RuntimeError, match='manual recovery'):
        w._publish_calibration_report(canonical, report, {canonical: old, alias: old})
    require(canonical.read_bytes() == intervening, 'rollback destroyed intervening canonical content')
    require((alias.read_bytes() if alias.exists() else None) == old, 'unwritten alias was modified')


def test_intervening_alias_change_is_detected_before_second_write(publication, monkeypatch):
    canonical, alias, report = publication
    newer = b'{"newer_alias":true}\n'
    original = w._atomic_write_json
    def alter_alias_after_first_write(path, data):
        original(path, data)
        if path == canonical:
            alias.write_bytes(newer)
    monkeypatch.setattr(w, '_atomic_write_json', alter_alias_after_first_write)
    with pytest.raises(ValueError, match='changed during publication'):
        w._publish_calibration_report(canonical, report, {canonical: None, alias: None})
    require(alias.read_bytes() == newer, 'a later target was overwritten after its entry snapshot changed')
    require(not canonical.exists(), 'own first write was not rolled back')


def test_lock_replacement_revokes_existing_context(publication):
    canonical, _, _ = publication
    lock = canonical.with_name('.' + canonical.name + '.lock')
    with pytest.raises(ValueError, match='revoked'):
        with w._calibration_output_lock(canonical) as verify:
            lock.unlink()
            lock.write_bytes(b'new owner')
            verify()
    require(lock.read_bytes() == b'new owner', 'old context removed replacement lock')


def test_revoked_lock_stops_pair_publication_and_preserves_new_owner(publication, monkeypatch):
    canonical, alias, report = publication
    lock = canonical.with_name('.' + canonical.name + '.lock')
    newer = b'{"new_owner_state":true}\n'
    original = w._atomic_write_json
    def revoke_after_first_write(path, data):
        original(path, data)
        if path == canonical:
            lock.unlink()
            lock.write_bytes(b'new owner')
            canonical.write_bytes(newer)
    monkeypatch.setattr(w, '_atomic_write_json', revoke_after_first_write)
    with pytest.raises(RuntimeError, match='manual recovery'):
        w._publish_calibration_report(canonical, report, {canonical: None, alias: None})
    require(canonical.read_bytes() == newer and not alias.exists(), 'revoked owner continued publication or rollback')
    require(lock.read_bytes() == b'new owner', 'revoked owner changed the new lock')


def test_partial_multibatch_posterior_matches_original_data_decimal_oracle():
    specifications = [*w.ENDPOINT_SPECS, *w.EXTENDED_STAGE0_STAGE1_SPECS]
    history = {row[0]: [] for row in specifications}
    state = None
    for batch in range(12):
        payload = {'batch_id': f'SYNTHETIC-FRESH-L5-{batch}'}
        for index, row in enumerate(specifications):
            if (batch + index) % 3 == 0:
                low, high = row[-3:-1]
                mean = low + (high - low) * ((batch + 1) / 14.)
                se = [.001, .019, 1e100, .3][(batch + index) % 4]
                payload[row[3]] = mean
                payload[row[4]] = se
                history[row[0]].append((mean, se))
        state = w.run_calibration(payload, True, state)
    all_states = {**state['conjugate_posterior_updates'],
                  **state['optional_pouch_and_36run_screen_updates_v80_070']}
    for row in specifications:
        with localcontext() as context:
            context.prec = 100
            d = lambda x: Decimal(str(x))
            precision = 1 / d(row[2])**2
            moment = precision * d(row[1])
            for mean, se in history[row[0]]:
                weight = 1 / d(se)**2
                precision += weight
                moment += weight * d(mean)
            mean, sd = float(moment / precision), float((1 / precision).sqrt())
        actual = all_states[row[0]]
        require(math.isclose(actual['posterior_mean_exact'], mean, rel_tol=2e-14, abs_tol=2e-14),
                f'partial batch sequence changed the sufficient statistic for {row[0]}')
        require(math.isclose(actual['posterior_sd_exact'], sd, rel_tol=2e-14),
                f'partial batch sequence changed precision for {row[0]}')
