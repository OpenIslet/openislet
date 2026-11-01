"""Cold-start arithmetic and provenance regression controls using synthetic evidence."""
from __future__ import annotations

import copy
import hashlib
import math
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import extract_ring1_corpus as ext
from scripts import rank_islet1_stacks as rank


@pytest.mark.parametrize("measurement", [
    ".5% versus 10%", "-20% versus 10%", "−20% versus 10%",
    "- 20% versus 10%", "1,00% versus 10%", "1e2-fold", "1e-2-fold",
    "1E+2 fold", "80 ± -2% versus 20 ± 1%", ".5 ng/mL",
    "1e2 ng/mL", "2.5e1 ng/mL", "1,00 ng/mL", "+20% versus 10%",
])
def test_unsupported_measurement_lexemes_cannot_become_different_numbers(measurement):
    parsed = ext.parse_quantitative_numbers(
        f"B2M knockout improved cell viability: {measurement} compared with controls.",
        "Protective_or_Improved",
    )
    assert parsed[:4] == (None, None, None, None)
    assert parsed[5] is None
    assert parsed[4] == "quarantined_unsupported_numeric_notation"


@pytest.mark.parametrize("quote,expected", [
    ("B2M knockout impaired cell viability: 0.5% versus 10%.", (0.5, None, 10.0, None, "%", -2.9957)),
    ("B2M knockout improved cell viability: 80% versus 20%.", (80.0, None, 20.0, None, "%", 1.3863)),
    ("B2M knockout improved cell viability: 1,000 ± 2 ng/mL versus 100 ± 1 ng/mL.", (1000.0, 2.0, 100.0, 1.0, "ng/mL", 2.3026)),
    ("B2M knockout improved cell viability by 2-fold compared with controls (p=1e-2).", (2.0, None, None, None, "fold_change_ratio", 0.6931)),
])
def test_canonical_measurements_and_scientific_p_values_remain_unchanged(quote, expected):
    direction = ext.classify_sentence_direction_and_penalty(quote)[0]
    assert ext.parse_quantitative_numbers(quote, direction) == expected


def _harmful_study():
    return {
        "study_metadata": {"pmid": "12345678", "source_text_type": "Abstract_Only",
                           "max_biological_fidelity_tier": "Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse"},
        "arms": [{"gate2_3_4_genetic_and_epitope_edits": [{"gene_symbol_hgnc": "B2M", "edit_modality": "CRISPR_Cas9_KO"}]}],
        "atomic_effects": [{"effect_id": "synthetic_harm", "focal_factor_X": "B2M knockout",
                            "verbatim_evidence_quote": "B2M knockout impaired cell viability.",
                            "outcome_variable_Y": "cell viability", "effect_direction": "Deleterious_or_Impaired",
                            "secretome_or_differentiation_penalty_observed": False,
                            "log_response_ratio_or_smd": None}],
    }


@pytest.mark.parametrize("group,key", [
    ("effect_direction_weights", "deleterious"),
    ("effect_direction_weights", "protective"),
    ("effect_direction_weights", "secretome_or_differentiation_penalty_addon"),
    ("effect_direction_weights", "shrunk_lrr_multiplier"),
    ("lrr_shrinkage_hyperparameters", "prior_weight"),
    ("lrr_shrinkage_hyperparameters", "per_empirical_study_slope"),
    ("lower_10pct_cred_bound_hyperparameters", "default_gene_se"),
])
@pytest.mark.parametrize("value", [-1.0, True, "0.75", float("nan"), float("inf")])
def test_profile_hyperparameters_fail_before_silent_coercion_or_sign_reversal(monkeypatch, group, key, value):
    original = rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"]
    mutant = copy.deepcopy(original)
    mutant[group]["parameter_value"][key] = value
    monkeypatch.setitem(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY, "auxiliary_evaluator_hyperparameters", mutant)
    with pytest.raises(ValueError):
        rank.build_target_evidence_profiles([_harmful_study()], {})


def test_harmful_evidence_remains_negative_with_valid_nonnegative_weights():
    profile = rank.build_target_evidence_profiles([_harmful_study()], {})["B2M"]
    assert profile["weighted_deleterious_penalty_score"] > 0
    assert profile["net_evidence_ko_score"] < 0
    assert math.isfinite(profile["posterior_se"]) and profile["posterior_se"] > 0


@pytest.mark.parametrize("canonical,alias", [("CD274", "PD-L1"), ("CASP9", "ICASP9")])
@pytest.mark.parametrize("duplicate_alias", [False, True])
def test_redaction_preserves_alias_target_identity_and_contrast_multiplicity(canonical, alias, duplicate_alias):
    study = _harmful_study()
    study["arms"][0]["gate2_3_4_genetic_and_epitope_edits"] = [
        {"gene_symbol_hgnc": canonical, "edit_modality": "Lentiviral_or_Transposon_OE"}
    ]
    effect = study["atomic_effects"][0]
    effect["focal_factor_X"] = f"{alias} overexpression"
    quote = f"{alias} overexpression impaired cell viability."
    effect["verbatim_evidence_quote"] = quote
    private = rank.build_target_evidence_profiles([study], {})
    token_genes = f"{canonical},{alias}" if duplicate_alias else alias
    effect["verbatim_evidence_quote"] = (
        f"[REDACTED_NON_OA_EXCERPT: see PMID:12345678 field=verbatim_evidence_quote "
        f"genes={token_genes} pert=GOF_OE quote_sha256={hashlib.sha256(quote.encode()).hexdigest()}]"
    )
    assert rank.validate_public_excerpt_token(effect["verbatim_evidence_quote"],
                                             expected_pmid="12345678", source_quote=quote)
    public = rank.build_target_evidence_profiles([study], {})
    assert public == private


@pytest.mark.parametrize("separator", ["–", "—", " – ", " — ", " ", "\u202f", "\u00a0", "\u2009", "_"])
def test_numeric_ranges_and_unsupported_groupings_cannot_be_scalar_suffixes(separator):
    parsed = ext.parse_quantitative_numbers(
        f"B2M knockout improved cell viability: 1{separator}500 ng/mL versus 100 ng/mL.",
        "Protective_or_Improved",
    )
    assert parsed == (None, None, None, None, "quarantined_unsupported_numeric_notation", None)


@pytest.mark.parametrize("qualifier", [">", ">=", "≥", "<", "<=", "≤", "~", "≈", "at least ", "more than ", "approximately "])
@pytest.mark.parametrize("position", ["treatment", "control", "dispersion"])
def test_qualified_paired_operands_never_become_exact_means_or_ratios(qualifier, position):
    if position == "treatment":
        measurement = f"{qualifier}80% versus 20%"
    elif position == "control":
        measurement = f"80% versus {qualifier}20%"
    else:
        measurement = f"80 ± {qualifier}2% versus 20 ± 1%"
    parsed = ext.parse_quantitative_numbers(
        f"B2M knockout improved cell viability: {measurement}.", "Protective_or_Improved",
    )
    assert parsed == (None, None, None, None, "quarantined_approximate_or_bounded_measurement", None)


@pytest.mark.parametrize("suffix", ["(p<0.05)", "(p < .05)", "(p≥0.05)", "(p=1e-2)", "after 2 days", "in 3 independent assays"])
def test_unrelated_numeric_context_does_not_quarantine_exact_operands(suffix):
    parsed = ext.parse_quantitative_numbers(
        f"B2M knockout improved cell viability: 80% versus 20% {suffix}.", "Protective_or_Improved",
    )
    assert parsed == (80.0, None, 20.0, None, "%", 1.3863)


def test_assay_ratio_before_decimal_fold_is_not_a_grouped_integer():
    assert ext.parse_quantitative_numbers(
        "B2M knockout impaired cell viability: 5:1 0.49-fold compared with controls.",
        "Deleterious_or_Impaired",
    ) == (0.49, None, None, None, "fold_change_ratio", -0.7133)


@pytest.mark.parametrize("prefix", ["for", "over", "up to", "more than", ">"])
def test_descriptive_censored_followup_remains_available_without_ratio(prefix):
    assert ext.parse_quantitative_numbers(
        f"B2M knockout prolonged graft survival {prefix} 100 days.", "Protective_or_Improved",
    ) == (100.0, None, None, None, "days", None)


@pytest.mark.parametrize("weight", ["deleterious", "protective", "secretome_or_differentiation_penalty_addon"])
def test_finite_weights_cannot_overflow_profile_accumulation(monkeypatch, weight):
    registry = copy.deepcopy(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"])
    registry["effect_direction_weights"]["parameter_value"][weight] = 1e308
    monkeypatch.setitem(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY, "auxiliary_evaluator_hyperparameters", registry)
    studies = []
    for index in range(100):
        study = _harmful_study()
        study["study_metadata"]["pmid"] = str(12345000 + index)
        effect = study["atomic_effects"][0]
        if weight == "protective":
            effect["effect_direction"] = "Protective_or_Improved"
            effect["verbatim_evidence_quote"] = "B2M knockout improved cell viability."
        effect["secretome_or_differentiation_penalty_observed"] = weight == "secretome_or_differentiation_penalty_addon"
        studies.append(study)
    with pytest.raises(ValueError, match="finite"):
        rank.build_target_evidence_profiles(studies, {})


@pytest.mark.parametrize("case", ["shrinkage_denominator", "omics_product", "omics_sum", "dispersion"])
def test_finite_inputs_cannot_hide_nonfinite_intermediate_profile_values(monkeypatch, case):
    registry = copy.deepcopy(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"])
    omics = {}
    if case == "shrinkage_denominator":
        registry["lrr_shrinkage_hyperparameters"]["parameter_value"].update(prior_weight=1e308, per_empirical_study_slope=1e308)
    elif case == "omics_product":
        registry["omics_constraint_weights"]["parameter_value"]["gnomad_loeuf_ko_slope"] = 1e308
        omics = {"B2M": {"gnomad_loeuf_upper": 2.0}}
    elif case == "omics_sum":
        registry["omics_constraint_weights"]["parameter_value"].update(gnomad_pli_ko_slope=1e308, opentargets_t1d_slope=1e308)
        omics = {"B2M": {"gnomad_pli": 1.0, "opentargets_t1d_score": 1.0}}
    else:
        registry["lower_10pct_cred_bound_hyperparameters"]["parameter_value"]["default_gene_se"] = 1e308
    monkeypatch.setitem(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY, "auxiliary_evaluator_hyperparameters", registry)
    study = _harmful_study()
    study["atomic_effects"][0]["log_response_ratio_or_smd"] = -1.0
    study["atomic_effects"][0]["endpoint_polarity_registry"] = {
        "polarity_class": "Higher_Is_Better_Protective_Endpoint", "orientation_sign": 1,
    }
    with pytest.raises(ValueError, match="finite"):
        rank.build_target_evidence_profiles([study], omics)


def test_large_finite_weights_remain_valid_when_derived_profile_is_finite(monkeypatch):
    registry = copy.deepcopy(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"])
    registry["effect_direction_weights"]["parameter_value"]["deleterious"] = 1e100
    monkeypatch.setitem(rank.COMPLETE_GATE_COEFFICIENT_REGISTRY, "auxiliary_evaluator_hyperparameters", registry)
    profile = rank.build_target_evidence_profiles([_harmful_study()], {})["B2M"]
    assert profile["weighted_deleterious_penalty_score"] == pytest.approx(7.2e98)
    assert all(math.isfinite(value) for value in profile.values() if isinstance(value, float))
