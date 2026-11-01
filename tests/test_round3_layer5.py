"""Independent Round-3 Layer-5 numerical and state-boundary oracles.

All observations, patients and publication failures are synthetic controls.
"""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import stat
from statistics import NormalDist
import sys
import unicodedata

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import ingest_wetlab_results as w
from scripts import product_contract as products
from scripts import run_layer5_personalized_patient_compiler as l5
from scripts import wetlab_contract as contract


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def observation(batch="SYNTHETIC-R3-A", value=.91, se=.001):
    return dict(batch_id=batch, obs_stage7_yield_ratio=value, obs_stage7_yield_se=se)


def patient(**changes):
    record = dict(patient_id="SYNTHETIC-R3", ancestry="EUR", nkg2a_neg_kir_educated_nk_fraction=.36,
                  pre_existing_gad65_ia2_znt8_t_cell_zscore=1., ifn1_ish_genetic_risk_percentile=45.,
                  sirpa_genotype="v1/v1")
    record.update(changes)
    return record


def reseal(report):
    report["conjugate_posterior_updates_sha256"] = hashlib.sha256(json.dumps(
        report["conjugate_posterior_updates"], sort_keys=True, allow_nan=False).encode()).hexdigest()
    report["calibration_state_sha256"] = w._calibration_state_sha256(report)


def decimal_oracle(pm, psd, observations):
    # Accumulate sufficient statistics from original data, never production posteriors.
    with localcontext() as context:
        context.prec = 100
        d = lambda value: Decimal(str(value))
        tau = 1 / d(psd) ** 2
        weighted_sum = tau * d(pm)
        for mean, se in observations:
            weight = 1 / d(se) ** 2
            tau += weight
            weighted_sum += weight * d(mean)
        return float(weighted_sum / tau), float((1 / tau).sqrt())


@pytest.mark.parametrize("endpoint_group,endpoint,field,bad", [
    ("conjugate_posterior_updates", "stage7_sc_beta_yield_fraction_arm1b", "posterior_precision_tau_exact", 1.),
    ("optional_pouch_and_36run_screen_updates_v80_070", "pouch_max_alkaline_ph_14d", "posterior_precision_tau_exact", 1.),
    ("conjugate_posterior_updates", "cd47_igv_hour0_fold_over_baseline_arm1b", "observed_mean_exact", 99.),
    ("optional_pouch_and_36run_screen_updates_v80_070", "pouch_max_alkaline_ph_14d", "observed_se_exact", .01),
    ("conjugate_posterior_updates", "stage7_sc_beta_yield_fraction_arm1b", "posterior_mean", .1),
    ("conjugate_posterior_updates", "stage7_sc_beta_yield_fraction_arm1b", "prior_data_conflict_detected", True),
])
def test_resealed_endpoint_cannot_contradict_exact_or_retained_state(endpoint_group, endpoint, field, bad):
    report = w.run_calibration(observation(value=.98), False)
    report[endpoint_group][endpoint][field] = bad
    reseal(report)
    with pytest.raises(ValueError, match="integrity|inconsistent|mismatch|retained|observed"):
        w.run_calibration(observation("SYNTHETIC-R3-B", .93), False, report)


@pytest.mark.parametrize("offset", [9e-8, -9e-8])
def test_resealed_posterior_sd_cannot_use_absolute_tolerance_as_precision_budget(offset):
    report = w.run_calibration(observation(), False)
    report["conjugate_posterior_updates"]["stage7_sc_beta_yield_fraction_arm1b"]["posterior_sd_exact"] += offset
    reseal(report)
    with pytest.raises(ValueError, match="integrity|mismatch"):
        w.run_calibration(observation("SYNTHETIC-R3-B", .93), False, report)


@pytest.mark.parametrize("parent", ["unrelated", "0"*64])
def test_genesis_state_cannot_claim_a_parent_after_resealing(parent):
    report = w.run_calibration(observation(), False)
    report["parent_report_sha256"] = parent
    reseal(report)
    with pytest.raises(ValueError, match="parent|lineage"):
        w.run_calibration(observation("SYNTHETIC-R3-B", .93), False, report)


def test_state_digest_binds_parent_identity():
    report = w.run_calibration(observation(), False)
    changed = copy.deepcopy(report)
    changed["parent_report_sha256"] = "0"*64
    require(w._calibration_state_sha256(report) != w._calibration_state_sha256(changed), "parent identity lies outside state digest")


def test_parent_bound_state_version_is_explicit_and_legacy_empirical_state_is_rejected():
    report=w.run_calibration(observation(),False)
    require(report["calibration_state_version"]==2,"parent-bound digest silently reused legacy state semantics")
    report["calibration_state_version"]=1
    reseal(report)
    with pytest.raises(ValueError,match="version=2.*original ordered observations"):
        w.run_calibration(observation("SYNTHETIC-R3-B",.93),False,report)


def test_resealed_algebraic_observation_must_match_latest_measurement_digest():
    report=w.run_calibration(observation(),False)
    endpoint=report["conjugate_posterior_updates"]["stage7_sc_beta_yield_fraction_arm1b"]
    # Recompute a completely valid conjugate update, but retain the hash of .91.
    endpoint.update(w.conjugate_normal_update(.985,.032,.95,.001,
        unit_and_support=endpoint["unit_and_support"],min_valid=0.,max_valid=2.))
    reseal(report)
    with pytest.raises(ValueError,match="measurement|payload"):
        w.run_calibration(observation("SYNTHETIC-R3-B",.95),False,report)


@pytest.mark.parametrize("group,endpoint,prior_sd", [
    ("conjugate_posterior_updates","cd47_igv_hour0_fold_over_baseline_arm1b",1e-8),
    ("optional_pouch_and_36run_screen_updates_v80_070","pouch_max_alkaline_ph_14d",1e-8),
])
def test_genesis_retained_prior_cannot_be_replaced_by_fabricated_precision(group,endpoint,prior_sd):
    report=w.run_calibration(observation(),False)
    state=report[group][endpoint]
    state.update(prior_sd_exact=prior_sd,prior_sd=prior_sd,posterior_sd_exact=prior_sd,posterior_sd=prior_sd,
                 posterior_precision_tau_exact=1/prior_sd**2,
                 posterior_95ci_low=round(state["prior_mean_exact"]-1.96*prior_sd,4),
                 posterior_95ci_high=round(state["prior_mean_exact"]+1.96*prior_sd,4))
    reseal(report)
    with pytest.raises(ValueError,match="genesis|precision|prior"):
        w.run_calibration(observation("SYNTHETIC-R3-B",.93),False,report)


@pytest.mark.parametrize("increment,forged_mean", [(0.,99.),(1e-8,7.8101),(1e-8,7.8099)])
def test_retained_multibatch_mean_must_be_reachable_from_genesis_precision(increment,forged_mean):
    first=w.run_calibration(observation(),False)
    report=w.run_calibration(observation("SYNTHETIC-R3-B",.93),False,first)
    endpoint=report["conjugate_posterior_updates"]["cd47_igv_hour0_fold_over_baseline_arm1b"]
    tau0=(1/.58)**2;sd=1/math.sqrt(tau0*(1+increment))
    endpoint.update(prior_mean_exact=forged_mean,posterior_mean_exact=forged_mean,
        prior_mean=round(forged_mean,4),posterior_mean=round(forged_mean,4),
        prior_sd_exact=sd,posterior_sd_exact=sd,prior_sd=round(sd,4),posterior_sd=round(sd,4),
        posterior_precision_tau_exact=(1/sd)**2,
        posterior_95ci_low=round(forged_mean-1.96*sd,4),posterior_95ci_high=round(forged_mean+1.96*sd,4))
    reseal(report)
    with pytest.raises(ValueError,match="reachable|history"):
        w.run_calibration(dict(batch_id="SYNTHETIC-R3-CD47-NEXT",obs_cd47_fold=7.8,obs_cd47_fold_se=.58),False,report)


@pytest.mark.parametrize("value,se", [(100.,5800.),(0.,5800.),(100.,1e8),(0.,1e8),(100.,1e300)])
def test_weak_bounded_observations_preserve_reachable_mean_roundoff(value,se):
    first=w.run_calibration(dict(batch_id="SYNTHETIC-R3-WEAK-A",obs_cd47_fold=value,obs_cd47_fold_se=se),False)
    second=w.run_calibration(observation("SYNTHETIC-R3-WEAK-B",.93),False,first)
    third=w.run_calibration(dict(batch_id="SYNTHETIC-R3-WEAK-C",obs_cd47_fold=7.8,obs_cd47_fold_se=.58),False,second)
    expected_mean,expected_sd=decimal_oracle(7.81,.58,[(value,se),(7.8,.58)])
    state=third["conjugate_posterior_updates"]["cd47_igv_hour0_fold_over_baseline_arm1b"]
    require(math.isclose(state["posterior_mean_exact"],expected_mean,rel_tol=2e-15),"legal weak observation failed or lost its mean")
    require(math.isclose(state["posterior_sd_exact"],expected_sd,rel_tol=2e-15),"legal weak observation failed precision oracle")


@pytest.mark.parametrize("value", [None, True, .5, "not-a-digest"])
def test_multibatch_state_requires_a_valid_parent_digest(value):
    first = w.run_calibration(observation(), False)
    report = w.run_calibration(observation("SYNTHETIC-R3-B", .93), False, first)
    report["parent_report_sha256"] = value
    reseal(report)
    with pytest.raises(ValueError, match="parent|lineage"):
        w.run_calibration(observation("SYNTHETIC-R3-C", .94), False, report)


def test_sequential_updates_match_decimal_original_data_sufficient_statistics():
    report = None
    samples = [(.90, .1), (.94, .005), (.92, .03), (1.01, .001), (.99, 1e200)]
    for index, (value, se) in enumerate(samples):
        report = w.run_calibration(observation(f"SYNTHETIC-R3-SEQ-{index}", value, se), False, report)
        actual = report["conjugate_posterior_updates"]["stage7_sc_beta_yield_fraction_arm1b"]
        mean, sd = decimal_oracle(.985, .032, samples[:index+1])
        require(math.isclose(actual["posterior_mean_exact"], mean, rel_tol=2e-15), "sequential posterior mean differs from original-data oracle")
        require(math.isclose(actual["posterior_sd_exact"], sd, rel_tol=2e-15), "sequential posterior SD differs from original-data oracle")


def test_all_eight_endpoint_chains_preserve_partial_missingness_against_original_data():
    specifications=[*w.ENDPOINT_SPECS,*w.EXTENDED_STAGE0_STAGE1_SPECS]
    observed={row[0]:[] for row in specifications}
    report=None
    for batch in range(3):
        payload={"batch_id":f"SYNTHETIC-R3-EIGHT-{batch}"}
        for index,row in enumerate(specifications):
            endpoint,mean,sd,mk,sk=row[:5]
            lower,upper=row[-3:-1]
            if (index+batch)%3==0:
                payload[mk]="UNAVAILABLE"
                continue
            value=min(upper,max(lower,mean+(batch-1)*sd/10))
            se=max(.001,sd/2)
            payload[mk]=value;payload[sk]=se;observed[endpoint].append((value,se))
        report=w.run_calibration(payload,False,report)
        states={**report["conjugate_posterior_updates"],**report["optional_pouch_and_36run_screen_updates_v80_070"]}
        for row in specifications:
            endpoint,mean,sd=row[:3]
            expected_mean,expected_sd=decimal_oracle(mean,sd,observed[endpoint])
            actual=states[endpoint]
            require(math.isclose(actual["posterior_mean_exact"],expected_mean,rel_tol=2e-15),f"{endpoint} mean lost its observed history")
            require(math.isclose(actual["posterior_sd_exact"],expected_sd,rel_tol=2e-15),f"{endpoint} missingness reset accumulated precision")


def test_many_legal_near_floor_batches_keep_adding_precision_without_rounding_state():
    report=None;samples=[]
    for index in range(128):
        value=.9+index*1e-6;samples.append((value,.001))
        report=w.run_calibration(observation(f"SYNTHETIC-R3-FLOOR-{index}",value),False,report)
    mean,sd=decimal_oracle(.985,.032,samples)
    actual=report["conjugate_posterior_updates"]["stage7_sc_beta_yield_fraction_arm1b"]
    require(math.isclose(actual["posterior_mean_exact"],mean,rel_tol=3e-15),"long chain accumulated mean rounding")
    require(math.isclose(actual["posterior_sd_exact"],sd,rel_tol=3e-15) and sd<.0001,"long chain lost legal near-floor precision")


def test_conflicting_batch_order_changes_only_final_step_heuristic_not_gaussian_state():
    endpoints=[]
    for sequence in [[.2,1.8],[1.8,.2]]:
        report=None
        for index,value in enumerate(sequence):
            report=w.run_calibration(observation(f"SYNTHETIC-R3-ORDER-{index}",value,.01),False,report)
        endpoints.append(report["conjugate_posterior_updates"]["stage7_sc_beta_yield_fraction_arm1b"])
    mean,sd=decimal_oracle(.985,.032,[(.2,.01),(1.8,.01)])
    for endpoint in endpoints:
        require(math.isclose(endpoint["posterior_mean_exact"],mean,rel_tol=2e-15),"conflict attenuation leaked into exact sequential Gaussian mean")
        require(math.isclose(endpoint["posterior_sd_exact"],sd,rel_tol=2e-15),"conflict attenuation changed declared observation precision")
        require("NOT_A_STUDENT_T_POSTERIOR" in endpoint["interval_and_robustness_qualification"]["student_t_df5_legacy_field"],"one-step comparator was mislabeled exact")
    require(endpoints[0]["conflict_attenuated_heuristic_mean"] != endpoints[1]["conflict_attenuated_heuristic_mean"],"order probe failed to distinguish final-step comparator")


@pytest.mark.parametrize("prior_sd,obs_se", [(1e-10, .001), (.5, .001), (1e140, .1), (.5, 1e300)])
def test_gaussian_and_conflict_heuristic_match_independent_scaled_oracles(prior_sd, obs_se):
    actual = w.conjugate_normal_update(.2, prior_sd, .8, obs_se, min_valid=0., max_valid=1.)
    mean, sd = decimal_oracle(.2, prior_sd, [(.8, obs_se)])
    require(math.isclose(actual["posterior_mean_exact"], mean, rel_tol=2e-15), "Gaussian mean differs from Decimal oracle")
    require(math.isclose(actual["posterior_sd_exact"], sd, rel_tol=2e-15), "Gaussian SD differs from Decimal oracle")
    with localcontext() as context:
        context.prec=100
        p, s, o, e = map(lambda v: Decimal(str(v)), (.2, prior_sd, .8, obs_se))
        z2 = (o-p)**2/(s*s+e*e)
        weight_prior = 1/(s*s)
        weight_obs = 1/(e*e*(1+z2/5))
        expected = float((weight_prior*p+weight_obs*o)/(weight_prior+weight_obs))
    heuristic=actual["conflict_attenuated_heuristic_mean"]
    require(heuristic == round(expected, 4), "conflict attenuation differs from independent formula")
    require(.2 <= heuristic <= actual["posterior_mean_exact"]+5.1e-5, "attenuation exceeds ordinary Gaussian displacement")
    json.dumps(actual, allow_nan=False)


def test_equal_huge_finite_means_do_not_overflow_precision_weighted_heuristic():
    result = w.conjugate_normal_update(1e308, 1., 1e308, 1., min_valid=0., max_valid=1.1e308)
    require(result["posterior_mean_exact"] == 1e308, "equal means moved")
    require(result["conflict_attenuated_heuristic_mean"] == 1e308, "finite convex mean overflowed")
    require(math.isclose(result["posterior_sd_exact"], math.sqrt(.5), rel_tol=1e-15), "precision addition changed")
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("family,spacer,pam", [("SpCas9", "A"*20, "CCC"), ("AsCas12a", "A"*23, "TTTT"),
                                               ("AsCas12aSpCas9", "A"*23, "TTTA")])
def test_declared_reference_pam_must_match_one_supported_nuclease_family(family, spacer, pam):
    target = pam+spacer if "Cas12a" in family else spacer+pam
    with pytest.raises(ValueError, match="PAM|nuclease"):
        l5.compute_split_junction_residual_cleavage_metrics(spacer, pam, family, target, target)


@pytest.mark.parametrize("family,spacer,pam", [("SpCas9", "ACGT"*5, "TGG"), ("AsCas12a", "ACGT"*5+"ACG", "TTTA")])
def test_junction_pam_position_and_seed_boundary_match_direct_mismatch_count(family, spacer, pam):
    target = pam+spacer if "Cas12a" in family else spacer+pam
    mutated = list(spacer)
    position = 18 if "Cas12a" in family else 9
    mutated[position-1] = "C" if mutated[position-1] != "C" else "A"
    mutant = "".join(mutated)
    changed = pam+mutant if "Cas12a" in family else mutant+pam
    result=l5.compute_split_junction_residual_cleavage_metrics(spacer,pam,family,target,changed)
    left,right=result["five_prime_ha_left_insert_junction"],result["three_prime_insert_ha_right_junction"]
    require(left["total_mismatches_vs_wt_spacer"] == 0 and right["mismatch_positions_1based"] == [position], "PAM orientation changed spacer coordinates")
    require(right["seed_region_mismatches"] == 1, "seed boundary off by one")
    require(result["full_20_to_23nt_target_absent_on_both_junctions"] is False, "exact window incorrectly called absent")


@pytest.mark.parametrize("score", [True, "65", -1., 101., float("inf"), float("nan")])
def test_reference_specificity_cannot_escape_finite_percentage_domain(score):
    with pytest.raises(ValueError, match="specificity|finite|numeric"):
        l5.compute_panel_and_reference_off_target_specificity("SYNTHETIC", "SpCas9", None, 0, [], score)


def test_assigned_panel_score_matches_independent_denominator_and_unknown_scan():
    result=l5.compute_panel_and_reference_off_target_specificity("SYNTHETIC", "SpCas9", None, 3, [.1,.2,.3], 65.)
    require(result["window_and_exome_panel_cfd_specificity_score"] == 62.5, "panel denominator changed")
    require(result["coding_exome_0_1_2_mismatch_hits"] is None, "unknown genome count became zero")


@pytest.mark.parametrize("nuclease,alternative_pam", [("Cas12a", "TTTA"), ("SpCas9", "TGG")])
@pytest.mark.parametrize("reverse", [False, True])
def test_donor_scan_rejects_exact_spacer_with_alternative_canonical_pam(monkeypatch,nuclease,alternative_pam,reverse):
    atlas=l5.build_part_a_crispr_grna_and_polymorphism_atlas()
    guide=next(g for g in atlas["guides"] if nuclease in g["nuclease"] and "PEmax" not in g["nuclease"])
    spacer=guide["protospacer_dna_5to3"]
    target=alternative_pam+spacer if nuclease=="Cas12a" else spacer+alternative_pam
    require(target != guide["full_target_site_5to3"],"probe does not vary the original PAM")
    if reverse:target=l5.revcomp(target)
    path=ROOT/"data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta"
    from scripts.run_layer2_perturbation_and_regulatory_plane import parse_construct_fasta_text
    records=parse_construct_fasta_text(path.read_text())
    donor=list(records)[-1];sequence=records[donor]
    records[donor]=sequence[:1500]+target+sequence[1500+len(target):]
    changed="".join(f">{key}\n{value}\n" for key,value in records.items())
    original=Path.read_text
    def altered(path_obj,*args,**kwargs):
        return changed if path_obj==path else original(path_obj,*args,**kwargs)
    monkeypatch.setattr(Path,"read_text",altered)
    with pytest.raises((ValueError,RuntimeError),match="target|Guide|donor"):
        l5.build_part_a_crispr_grna_and_polymorphism_atlas()


@pytest.mark.parametrize("aliases", [["DUPLICATE","DUPLICATE"], ["ISLET1_A4_ARM1A"]])
def test_product_standalone_validator_rejects_duplicate_identity_aliases(aliases):
    product=products.get_product_definition("ISLET1_A4_ARM1A")
    product["aliases"]=aliases
    with pytest.raises(ValueError, match="unique|alias"):
        products.validate_product_definition(product)


@pytest.mark.parametrize("panel", [[None], [True], [""], ["SPACER", "SPACER"], [{"guide":"SPACER"}]])
def test_regimen_panel_is_a_unique_nonempty_string_vector(panel):
    record=patient();regimen=l5.compile_personalized_patient_blueprint(record)
    regimen["compiled_grna_panel"]=panel
    with pytest.raises(ValueError, match="guide|panel"):
        l5.evaluate_regimen_survival_probability(record, regimen)


def test_hwe_inverse_cdf_matches_exact_grid_masses_and_boundary_sides():
    p=.2
    u=(np.arange(10000)+.5)/10000
    counts=np.bincount([l5._hwe_alternate_alleles(float(q),p) for q in u],minlength=3)
    require(counts.tolist() == [6400,3200,400], "HWE masses do not equal (1-p)^2, 2p(1-p), p^2")
    for boundary in (p*p,p*p+2*p*(1-p)):
        below=l5._hwe_alternate_alleles(math.nextafter(boundary,-math.inf),p)
        at=l5._hwe_alternate_alleles(boundary,p)
        require(below == at+1, "closed cumulative threshold lost one genotype boundary")


@pytest.fixture(scope="module")
def population_receipts():
    counts={};edit_counts={};copula=[[] for _ in range(4)];copula_calls=0
    def capture(record):
        nonlocal copula_calls
        regimen=l5.compile_personalized_patient_blueprint(record)
        if record["patient_id"].startswith("SIM-"):
            ancestry=record["ancestry"];counts[ancestry]=counts.get(ancestry,0)+1
            product=regimen["product_definition"];key=f'K{product["K"]}_B{product["B"]}'
            edit_counts[key]=edit_counts.get(key,0)+1
        elif record["patient_id"].startswith("COP-"):
            # Recover the IFN latent normal from its percentile; other traits are
            # measured affine normals with the compiler's declared end clipping.
            copula[copula_calls//10000].append((
                record["nkg2a_neg_kir_educated_nk_fraction"],
                record["pre_existing_gad65_ia2_znt8_t_cell_zscore"],
                NormalDist().inv_cdf(record["ifn1_ish_genetic_risk_percentile"]/100)))
            copula_calls+=1
        return regimen
    population=l5.simulate_10000_patient_global_population_trial(capture)
    return population,counts,edit_counts,copula


def test_population_allocations_sum_to_actual_10000_compiler_calls(population_receipts):
    population,counts,edit_counts,_=population_receipts
    require(counts == dict.fromkeys(["EUR","AFR","EAS","SAS","AMR"],2000),"population omitted or duplicated actual patient calls")
    actual=population["editing_complexity_and_manufacturing_sensitivity_audit_v80_009"]["kb_edit_count_distribution_across_10000_patients"]
    require(actual==edit_counts and sum(actual.values())==10000,"aggregate K/B counts do not match actual compiled product identities")
    for cohort in population["ancestry_stratified_cohorts"]:
        allocation=cohort["compiler_regimen_allocation_counts"]
        require(allocation["track1_off_the_shelf"]+allocation["track2_stratified_subline"]==2000,"track allocations fail conservation")


def test_copula_has_unit_normal_marginals_and_declared_pairwise_correlation(population_receipts):
    population,_,_,copula=population_receipts
    for index,loading in enumerate([0.,.25,.5,.75]):
        sample=np.array(copula[index])
        require(sample.shape==(10000,3),"copula sensitivity omitted compiler calls")
        latent_ifn=sample[:,2]
        require(abs(float(np.mean(latent_ifn)))<.04 and abs(float(np.var(latent_ifn))-1)<.07,"shared factor changed marginal normal distribution")
        correlations=np.corrcoef(sample,rowvar=False)
        # Clipping plus ancestry mean mixing attenuates observed correlations;
        # this tolerance includes sampling error, without requiring exact draws.
        pairwise=correlations[np.triu_indices(3,k=1)]
        require(float(np.max(abs(pairwise-loading**2)))<.065,"covariate correlation is inconsistent with shared-factor loading squared")
        row=population["correlated_frailty_copula_sensitivity_sol_v80_r2_07"][index]
        require(row["analytic_pairwise_trait_correlation_rho_r_squared"]==loading**2,"reported analytic covariance differs from Gaussian factor model")


def test_measurement_hash_identifies_decimal_equivalence_and_metadata_independence():
    first=w.run_calibration(observation(value=.9),False)
    with pytest.raises(ValueError, match="replay"):
        w.run_calibration(dict(batch_id="SYNTHETIC-R3-REPLAY",assay_id="renamed",obs_stage7_yield_ratio="9e-1",obs_stage7_yield_se="1e-3"),False,first)


@pytest.fixture
def filesystem_root(tmp_path, monkeypatch):
    # Keeps audit-run test files inside its explicit --basetemp evidence directory.
    monkeypatch.setattr(w.tempfile,"gettempdir",lambda:str(tmp_path))
    return tmp_path


def test_legacy_alias_cannot_be_published_as_an_independent_posterior(filesystem_root, monkeypatch):
    canonical=filesystem_root/"canonical.json";alias=filesystem_root/"legacy.json"
    monkeypatch.setattr(w,"DEFAULT_OUT",canonical)
    monkeypatch.setattr(w,"LEGACY_ALIAS_PATH",alias)
    first=w.run_calibration(observation(),False)
    with pytest.raises(ValueError,match="alias|canonical"):
        w._publish_calibration_report(alias,first,{alias:None})
    require(not alias.exists() and not canonical.exists(),"independent alias write changed publication state")


@pytest.mark.parametrize("existing",[False,True])
@pytest.mark.parametrize("spelling",["LEGACY.JSON",unicodedata.normalize("NFD","légacy.json")])
def test_reserved_alias_identity_covers_case_and_unicode_before_creation(filesystem_root,monkeypatch,existing,spelling):
    canonical=filesystem_root/"canonical.json"
    alias=filesystem_root/("légacy.json" if "é" in unicodedata.normalize("NFC",spelling) else "legacy.json")
    variant=filesystem_root/spelling
    monkeypatch.setattr(w,"DEFAULT_OUT",canonical);monkeypatch.setattr(w,"LEGACY_ALIAS_PATH",alias)
    first=w.run_calibration(observation(),False)
    before=None
    if existing:
        before=(json.dumps(first,sort_keys=True,allow_nan=False)+"\n").encode()
        canonical.write_bytes(before);alias.write_bytes(before)
    child=w.run_calibration(observation("SYNTHETIC-R3-CASE-B",.93),False,first)
    with pytest.raises(ValueError,match="alias|canonical|identity"):
        w._publish_calibration_report(variant,child,{variant:before})
    require((canonical.read_bytes() if canonical.exists() else None)==before,"alias variant altered canonical")
    require((alias.read_bytes() if alias.exists() else None)==before,"alias variant altered legacy output")


@pytest.mark.parametrize("existing",[False,True])
def test_reserved_canonical_case_variant_cannot_choose_independent_writer(filesystem_root,monkeypatch,existing):
    canonical=filesystem_root/"canonical.json";alias=filesystem_root/"legacy.json";variant=filesystem_root/"CANONICAL.JSON"
    monkeypatch.setattr(w,"DEFAULT_OUT",canonical);monkeypatch.setattr(w,"LEGACY_ALIAS_PATH",alias)
    first=w.run_calibration(observation(),False);before=None
    if existing:
        before=(json.dumps(first,sort_keys=True,allow_nan=False)+"\n").encode()
        canonical.write_bytes(before);alias.write_bytes(before)
    child=w.run_calibration(observation("SYNTHETIC-R3-CASE-B",.93),False,first)
    with pytest.raises(ValueError,match="canonical|identity"):
        w._publish_calibration_report(variant,child,{variant:before})
    require((canonical.read_bytes() if canonical.exists() else None)==before,"canonical variant selected an independent writer")
    require((alias.read_bytes() if alias.exists() else None)==before,"canonical variant split legacy state")


def test_rollback_of_new_output_fsyncs_directory_after_unlink(filesystem_root,monkeypatch):
    target=filesystem_root/"canonical.json";alias=filesystem_root/"legacy.json"
    monkeypatch.setattr(w,"DEFAULT_OUT",target);monkeypatch.setattr(w,"LEGACY_ALIAS_PATH",alias)
    first=w.run_calibration(observation(),False)
    events=[]
    original_fsync,original_unlink,original_write=w.os.fsync,w.os.unlink,w._atomic_write_json
    def trace_fsync(fd):
        if stat.S_ISDIR(os.fstat(fd).st_mode):events.append("dir_fsync")
        return original_fsync(fd)
    def trace_unlink(name,*args,**kwargs):
        result=original_unlink(name,*args,**kwargs)
        if name==target.name:events.append("rollback_unlink")
        return result
    def fail_alias(path,data):
        if path==alias:raise OSError("SYNTHETIC alias failure")
        return original_write(path,data)
    monkeypatch.setattr(w.os,"fsync",trace_fsync);monkeypatch.setattr(w.os,"unlink",trace_unlink)
    monkeypatch.setattr(w,"_atomic_write_json",fail_alias)
    with pytest.raises(OSError,match="SYNTHETIC"):
        w._publish_calibration_report(target,first,{target:None,alias:None})
    require(not target.exists() and not alias.exists(),"failed initial publication left an output")
    index=events.index("rollback_unlink")
    require("dir_fsync" in events[index+1:],"rollback deletion was not durable at directory boundary")


def test_dec008_bridge_is_absent_from_current_compiler_output():
    result=l5.compile_personalized_patient_blueprint(patient())
    text=json.dumps(result).lower()
    require("anakinra" not in text and "etanercept" not in text,"removed systemic combination reappeared in compiler output")
