#!/usr/bin/env python3
"""
Closed-Loop Wet-Lab Results Ingestor & Bayesian Posterior Updater (`scripts/ingest_wetlab_results.py`, v0.8.0)
==============================================================================================================
Ingests empirical wet-lab bench readouts (JSON) from:
  1. The 36-Run `HIP-Core` Fractional Factorial Perturbation Screen (Stage-7 SC-beta yield, GSIS index,
     PBMC/NK allogeneic & autoreactive T-cell cytotoxicity, and 24h AP1903 iCasp9 clearance logs).
  2. SERPINB9-2A C-terminal tail Granzyme B second-order inhibition rate constant `k_ass` (`M^-1 s^-1`).
  3. 350 um PDMS-CaO2 + Catalase hydrogel pouch oxygen/peroxide/pH qualification readouts.

Fails closed on empty, missing, non-finite (`NaN`/`Inf`), non-positive or sub-floor SE (`se < 1e-3`),
duplicate JSON keys, duplicate `batch_id`s, or content-identical measurement payload replays
(`batch_payload_sha256`), and reports 3-sigma prior-data conflict diagnostics (`z_conflict > 3.0`)
alongside reconciled Layer-2/Layer-3 reference model values (closing AUDIT-REF, AUDIT-REF,
AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, and
AUDIT-REF).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import tempfile
import stat
import fcntl
import sys
import uuid
import unicodedata
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
DEFAULT_OUT = REPO_ROOT / "data" / "extracted" / "wetlab_bayesian_posterior_calibration_v0.8.0.json"
LEGACY_ALIAS_PATH = REPO_ROOT / "data" / "extracted" / "wetlab_bayesian_posterior_calibration.json"
LAYER2_ATLAS_PATH = REPO_ROOT / "data" / "extracted" / "layer2_regulatory_and_perturbation_atlas.json"

# Minimum biological standard-error floor (AUDIT-REF): prevents precision/IG_lab gaming via tiny SEs on observations.
MIN_BIOLOGICAL_SE_FLOOR = 1e-3
# Positivity floor on prior_sd so valid sequential precision accumulation (posterior_sd < 1e-3 after multiple legal obs_se >= 1e-3 updates) never fails closed (AUDIT-REF):
MIN_PRIOR_SD_POSITIVITY_FLOOR = 1e-12
PRIOR_DATA_CONFLICT_Z_THRESHOLD = 3.0
# Version 2 additionally binds parent identity and validates measurement/state algebra.
# Legacy empirical states require reconstruction from their ordered original observations.
CALIBRATION_STATE_VERSION = 2

_MEASUREMENT_CONTRACT_MODULE = None

SUPPORTED_EMPIRICAL_LINEAGES = frozenset({"PURE_EMPIRICAL_WETLAB_LINEAGE", "PURE_EMPIRICAL_LINEAGE"})


ENDPOINT_SPECS = [
    (
        "stage7_sc_beta_yield_fraction_arm1b",
        0.985,
        0.032,
        "obs_stage7_yield_ratio",
        "obs_stage7_yield_se",
        0.0,
        2.0,
        "normalized_ratio_to_unedited_wt_control_closed_interval_[0.0, 2.0]",
    ),
    (
        "cd47_igv_hour0_fold_over_baseline_arm1b",
        7.81,
        0.58,
        "obs_cd47_fold",
        "obs_cd47_fold_se",
        0.0,
        100.0,
        "fold_over_endogenous_baseline_closed_interval_[0.0, 100.0]",
    ),
    (
        "serpinb9_p2a_gzmb_log10_kass_M_inv_s_inv",
        6.18,
        0.25,
        "obs_serpinb9_log10_kass",
        "obs_serpinb9_log10_kass_se",
        0.0,
        10.0,
        "log10_second_order_rate_constant_M_inv_s_inv_closed_interval_[0.0, 10.0]",
    ),
    (
        "aa_vs1_ucoe_icasp9_ap1903_24h_log10_kill",
        4.42,
        0.18,
        "obs_icasp9_log10_kill",
        "obs_icasp9_log10_kill_se",
        0.0,
        10.0,
        "log10_clonal_depletion_kill_depth_closed_interval_[0.0, 10.0]",
    ),
]

EXTENDED_STAGE0_STAGE1_SPECS = [
    (
        "pouch_min_slab_po2_14d_mmhg",
        24.22,
        1.80,
        "obs_pouch_min_po2_mmhg",
        "obs_pouch_min_po2_se",
        0.85,
        0.0,
        160.0,
        "min_14d_slab_po2_mmhg_closed_interval_[0.0, 160.0]",
    ),
    (
        "pouch_max_peroxide_14d_uM",
        0.0175,
        0.010,
        "obs_pouch_max_h2o2_uM",
        "obs_pouch_max_h2o2_se",
        0.005,
        0.0,
        100.0,
        "max_14d_slab_h2o2_uM_closed_interval_[0.0, 100.0]",
    ),
    (
        "pouch_max_alkaline_ph_14d",
        7.400,
        0.040,
        "obs_pouch_max_ph",
        "obs_pouch_max_ph_se",
        0.020,
        6.0,
        9.0,
        "max_14d_pouch_alkaline_ph_closed_interval_[6.0, 9.0]",
    ),
    (
        "ifng_48h_hla_dr_dq_negative_pct_arm1b",
        99.60,
        0.35,
        "obs_ifng_hla_dr_dq_negative_pct",
        "obs_ifng_hla_dr_dq_negative_se",
        0.25,
        0.0,
        100.0,
        "ifng_48h_hla_dr_dq_negative_cell_percentage_closed_interval_[0.0, 100.0]",
    ),
]

ENDPOINT_PROVENANCE_AND_MODEL_RECONCILIATION: Dict[str, Dict[str, Any]] = {'stage7_sc_beta_yield_fraction_arm1b': {'prior_source': 'Independent preregistered engineering target '
                                                         'prior; not an observed result for this exact '
                                                         'product.',
                                         'layer2_or_layer3_model_reference_value': None,
                                         'mapped_ranker_prior_key': 'Gate_1_Endocrine_Function_and_Metabolic_Fitness.stage7_sc_beta_relative_yield_arm1b',
                                         'model_reference_status': 'UNAVAILABLE_OR_PREREGISTERED_PENDING'},
 'cd47_igv_hour0_fold_over_baseline_arm1b': {'prior_source': 'Independent preregistered engineering target '
                                                             'prior; not an observed result for this exact '
                                                             'product.',
                                             'layer2_or_layer3_model_reference_value': None,
                                             'layer2_arm1b_molecules_per_cell': None,
                                             'layer2_wt_molecules_per_cell': None,
                                             'mapped_ranker_prior_key': 'Gate_3_Innate_NK_Macrophage_and_Complement_Stealth.cd47_igv_hour0_fold_over_wt_baseline',
                                             'model_reference_status': 'UNAVAILABLE_OR_PREREGISTERED_PENDING'},
 'serpinb9_p2a_gzmb_log10_kass_M_inv_s_inv': {'prior_source': 'Independent preregistered engineering target '
                                                              'prior; not an observed result for this exact '
                                                              'product.',
                                              'layer2_or_layer3_model_reference_value': None,
                                              'mapped_ranker_prior_key': 'p2a_c_terminal_tail_serpinb9_activity_retention_factor_r01',
                                              'model_reference_status': 'UNAVAILABLE_OR_PREREGISTERED_PENDING'},
 'aa_vs1_ucoe_icasp9_ap1903_24h_log10_kill': {'prior_source': 'Independent preregistered engineering target '
                                                              'prior; not an observed result for this exact '
                                                              'product.',
                                              'layer2_or_layer3_model_reference_value': None,
                                              'mapped_ranker_prior_key': 'Gate_5_Clonal_Safety_Genomic_Stability_and_Failsafe_Lineage_Control.icasp9_24h_ap1903_log10_kill',
                                              'model_reference_status': 'UNAVAILABLE_OR_PREREGISTERED_PENDING'}}

VALID_CALIBRATION_MODES = {
    "EMPIRICAL_WETLAB_INGESTION",
    "SYNTHETIC_DRY_RUN_VERIFICATION_DO_NOT_CLAIM_AS_WETLAB_DATA",
}


def reject_duplicate_json_keys(pairs: List[Tuple[str, Any]]) -> Dict[str, Any]:
    """
    Strict JSON `object_pairs_hook` that rejects duplicate keys within any JSON object
    (preventing last-wins key shadowing in wet-lab payloads; closes AUDIT-REF).
    """
    out: Dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise ValueError(f"Duplicate JSON key {k!r} detected in payload (AUDIT-REF)")
        out[k] = v
    return out


def loads_strict_json(raw_text: str) -> Any:
    """Parse JSON rejecting duplicate keys and NaN/Inf constants."""
    def _reject_bad_const(c: str) -> Any:
        raise ValueError(f"Non-finite JSON constant {c!r} is prohibited")

    return json.loads(raw_text, object_pairs_hook=reject_duplicate_json_keys, parse_constant=_reject_bad_const)


def _require_finite_float(val: Any, name: str) -> float:
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        raise ValueError(f"{name} must be a finite numeric value, got {val!r}")
    fval = float(val)
    if not math.isfinite(fval):
        raise ValueError(f"{name} must be finite (not NaN/Inf), got {val!r}")
    return fval


def _safe_round_positive_sd(sd_val: float) -> float:
    r4 = round(sd_val, 4)
    return r4 if r4 > 0.0 else float(sd_val)


def _get_stage0_decision_tree_specs() -> list:
    global _MEASUREMENT_CONTRACT_MODULE
    if _MEASUREMENT_CONTRACT_MODULE is None:
        spec = importlib.util.spec_from_file_location("openislet_measurement_contract", REPO_ROOT / "scripts/wetlab_contract.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        _MEASUREMENT_CONTRACT_MODULE = helper
    return list(getattr(_MEASUREMENT_CONTRACT_MODULE, "STAGE0_DECISION_TREE_ENDPOINT_SPECS", []))


CANONICAL_MEASUREMENT_KEYS: frozenset[str] = frozenset(
    {r[3] for r in ENDPOINT_SPECS}
    | {r[4] for r in ENDPOINT_SPECS}
    | {r[3] for r in EXTENDED_STAGE0_STAGE1_SPECS}
    | {r[4] for r in EXTENDED_STAGE0_STAGE1_SPECS}
    | {r[1] for r in _get_stage0_decision_tree_specs()}
    | {r[2] for r in _get_stage0_decision_tree_specs()}
)

# Rare-event count/rate/denominator assay groups that must be updated or cleared atomically across multi-batch pooling (AUDIT-REF & F2-03; DEC-071):
RARE_EVENT_ASSAY_GROUPS: tuple[tuple[str, ...], ...] = (
    (
        "obs_tp53_ppm1d_bcor_vaf",
        "obs_tp53_ppm1d_bcor_duplex_depth",
        "obs_tp53_ppm1d_bcor_mutant_molecules_k",
        "obs_tp53_vaf",
        "obs_tp53_duplex_depth",
        "obs_tp53_mutant_molecules_k",
    ),
    (
        "obs_castseq_translocation_freq",
        "obs_castseq_duplex_or_amplicon_depth",
        "obs_castseq_mutant_junctions_k",
    ),
    (
        "obs_mcb_effective_clonogenic_trials",
        "obs_mcb_clonogenic_survivor_count_k",
    ),
    (
        "obs_residual_psc_effective_fraction",
        "obs_residual_psc_effective_trials",
    ),
    (
        "obs_hazard_state_prolif_survivor_count_k",
        "obs_hazard_state_prolif_effective_trials",
    ),
    (
        "obs_ap1903_step1_sublethal_survivor_count_k",
        "obs_ap1903_step1_effective_trials",
    ),
    (
        "obs_dxs_locus4_joint_retention_fraction",
        "obs_dxs_cd47low_icasp9low_fraction",
        "obs_dxs_ap1903_cd47high_escape_fraction",
        "obs_dxs_screened_cell_count",
        "obs_dxs_clonogenic_fraction",
        "obs_dxs_spike_in_detection_efficiency",
    ),
    (
        "obs_20q11_bcl2l1_single_cell_mosaic_fraction",
        "obs_20q11_bcl2l1_single_cell_mosaic_cells_scored",
    ),
    (
        "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct",
        "obs_gate0v_day28_persistent_vp1_ie1_cells_scored",
    ),
    (
        "obs_assembled_cluster_ki67_progenitor_pct",
        "obs_assembled_cluster_ki67_progenitor_cells_scored",
    ),
    (
        "obs_ductal_sox9_krt19_double_pos_pct",
        "obs_ductal_sox9_krt19_double_pos_cells_scored",
    ),
    (
        "obs_assembled_cluster_slc18a1_enterochromaffin_pct",
        "obs_assembled_cluster_slc18a1_enterochromaffin_cells_scored",
    ),
    (
        "obs_senescent_p16_sasp_positive_cell_pct",
        "obs_senescent_p16_sasp_cells_scored",
    ),
)


def _drop_unpaired_rare_event_members(pooled: Dict[str, Any], batch_meas: Dict[str, Any]) -> None:
    """When any member of a rare-event count/rate/denominator assay group is updated in batch_meas, drop any stale sibling members from prior batches that were not re-supplied in the same batch."""
    mean_to_se = {r[1]: r[2] for r in _get_stage0_decision_tree_specs()}
    for grp in RARE_EVENT_ASSAY_GROUPS:
        supplied_in_batch = {
            k for k in grp if k in batch_meas and batch_meas.get(k) not in ("UNAVAILABLE", None)
        }
        if supplied_in_batch:
            # Keep TP53 aliases synchronized if either canonical or alias key was supplied:
            if grp[0] == "obs_tp53_ppm1d_bcor_vaf":
                vaf_supplied = bool({"obs_tp53_ppm1d_bcor_vaf", "obs_tp53_vaf"} & supplied_in_batch)
                depth_supplied = bool({"obs_tp53_ppm1d_bcor_duplex_depth", "obs_tp53_duplex_depth"} & supplied_in_batch)
                k_supplied = bool({"obs_tp53_ppm1d_bcor_mutant_molecules_k", "obs_tp53_mutant_molecules_k"} & supplied_in_batch)
                if not vaf_supplied:
                    for k_drop in ("obs_tp53_ppm1d_bcor_vaf", "obs_tp53_vaf"):
                        pooled.pop(k_drop, None)
                        pooled.pop(mean_to_se.get(k_drop, ""), None)
                if not depth_supplied:
                    for k_drop in ("obs_tp53_ppm1d_bcor_duplex_depth", "obs_tp53_duplex_depth"):
                        pooled.pop(k_drop, None)
                        pooled.pop(mean_to_se.get(k_drop, ""), None)
                if not k_supplied:
                    for k_drop in ("obs_tp53_ppm1d_bcor_mutant_molecules_k", "obs_tp53_mutant_molecules_k"):
                        pooled.pop(k_drop, None)
                        pooled.pop(mean_to_se.get(k_drop, ""), None)
            else:
                for k in grp:
                    if k not in supplied_in_batch:
                        pooled.pop(k, None)
                        if k in mean_to_se:
                            pooled.pop(mean_to_se[k], None)


def _extract_canonical_measurements(input_payload: Dict[str, Any]) -> Dict[str, Any]:
    measurement_canonical: Dict[str, Any] = {}
    primary_means = {r[3] for r in ENDPOINT_SPECS} | {r[3] for r in EXTENDED_STAGE0_STAGE1_SPECS}
    for _ep_name, _def_m, _def_sd, mean_key, se_key, _min_v, _max_v, _unit in ENDPOINT_SPECS:
        m_val = input_payload.get(mean_key)
        s_val = input_payload.get(se_key)
        if m_val in ("UNAVAILABLE", None):
            measurement_canonical[mean_key] = "UNAVAILABLE"
            measurement_canonical[se_key] = "UNAVAILABLE"
        else:
            measurement_canonical[mean_key] = float(m_val) if isinstance(m_val, (int, float)) and not isinstance(m_val, bool) else m_val
            measurement_canonical[se_key] = float(s_val) if isinstance(s_val, (int, float)) and not isinstance(s_val, bool) else s_val
    for _ep_name, _def_m, _def_sd, mean_key, se_key, default_se, _min_v, _max_v, _unit in EXTENDED_STAGE0_STAGE1_SPECS:
        if mean_key in input_payload and input_payload.get(mean_key) not in ("UNAVAILABLE", None):
            m_val = input_payload.get(mean_key)
            s_val = input_payload.get(se_key, default_se)
            measurement_canonical[mean_key] = float(m_val) if isinstance(m_val, (int, float)) and not isinstance(m_val, bool) else m_val
            measurement_canonical[se_key] = float(s_val) if isinstance(s_val, (int, float)) and not isinstance(s_val, bool) else s_val
    stage0_custom_defaults = getattr(_MEASUREMENT_CONTRACT_MODULE, "STAGE0_CUSTOM_DEFAULT_SE", {})
    _scale_aware_s0_se = getattr(
        _MEASUREMENT_CONTRACT_MODULE,
        "_scale_aware_stage0_default_se",
        lambda mk, l_b, h_b: stage0_custom_defaults.get(mk, max(1e-3, (h_b - l_b) * 0.01)),
    )
    for _eid, mean_key, se_key, lo, hi, _unit in _get_stage0_decision_tree_specs():
        if mean_key in primary_means:
            continue
        if mean_key in input_payload and input_payload.get(mean_key) not in ("UNAVAILABLE", None):
            m_val = input_payload.get(mean_key)
            default_se = _scale_aware_s0_se(mean_key, lo, hi)
            s_val = input_payload.get(se_key, default_se)
            measurement_canonical[mean_key] = float(m_val) if isinstance(m_val, (int, float)) and not isinstance(m_val, bool) else m_val
            measurement_canonical[se_key] = float(s_val) if isinstance(s_val, (int, float)) and not isinstance(s_val, bool) else s_val
    return {
        key: (0.0 if isinstance(value, (int, float)) and not isinstance(value, bool) and value == 0.0 else value)
        for key, value in measurement_canonical.items()
    }


def compute_measurement_payload_sha256(input_payload: Dict[str, Any]) -> str:
    """
    Compute content-based SHA-256 over the canonical measurement fields only (excluding
    `batch_id`/`assay_id` label strings) so content-identical measurement replays under a
    new `batch_id` are deterministically detected and rejected (closing AUDIT-REF
    and AUDIT-REF across all 8 canonical Stage-0/Stage-1 endpoints).
    """
    measurement_canonical = _extract_canonical_measurements(input_payload)
    raw_measurement_bytes = json.dumps(measurement_canonical, sort_keys=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw_measurement_bytes).hexdigest()


def _convex_mean(left: float, right: float, left_weight: float, right_weight: float) -> float:
    """Preserve the exact convex hull, including equal physical-boundary means."""
    # Retain the previous arithmetic for in-hull values so existing version-2
    # reports at display-rounding ties remain valid continuation inputs.
    mean = left_weight * left + right_weight * right
    if not math.isfinite(mean):
        mean = left + right_weight * (right - left)
    if not math.isfinite(mean):
        raise ValueError("Normal update mean is outside representable numeric range")
    return min(max(left, right), max(min(left, right), mean))


def _normal_normal_moments(pm: float, psd: float, om: float, ose: float) -> tuple[float, float, float, float, float]:
    """Finite precision arithmetic without squaring a potentially huge finite SE."""
    tau_prior = (1.0 / psd) ** 2
    tau_obs = (1.0 / ose) ** 2
    tau_post = tau_prior + tau_obs
    if tau_prior <= 0.0 or tau_post <= 0.0 or not all(math.isfinite(v) for v in (tau_prior, tau_obs, tau_post)):
        raise ValueError("Normal update precision is outside representable numeric range")
    post_mean = _convex_mean(pm, om, tau_prior / tau_post, tau_obs / tau_post)
    post_sd = 1.0 / math.sqrt(tau_post)
    if not math.isfinite(post_mean) or not math.isfinite(post_sd) or post_sd <= 0.0:
        raise ValueError("Normal update moments are outside representable numeric range")
    return post_mean, post_sd, tau_prior, tau_obs, tau_post


def conjugate_normal_update(
    prior_mean: float,
    prior_sd: float,
    obs_mean: float,
    obs_se: float,
    unit_and_support: str = "",
    min_valid: float = 0.0,
    max_valid: float = 100.0,
) -> Dict[str, Any]:
    """
    Exact Normal-Normal conjugate Bayesian update with strict finite & biological SE-floor
    validation (`obs_se >= 1e-3`, `prior_sd > 1e-12`), 3-sigma prior-data conflict detection (`z_conflict > 3.0`),
    explicitly labeled clipping display (`[min_valid, max_valid]`, not a truncated-Normal posterior),
    heuristic conflict attenuation comparator (not a fitted Student-t posterior), and full-precision machine state
    preservation (closing AUDIT-REF, AUDIT-REF, 087, 103, and AUDIT-REF).
    """
    pm = _require_finite_float(prior_mean, "prior_mean")
    psd = _require_finite_float(prior_sd, "prior_sd")
    om = _require_finite_float(obs_mean, "obs_mean")
    ose = _require_finite_float(obs_se, "obs_se")
    lower = _require_finite_float(min_valid, "min_valid")
    upper = _require_finite_float(max_valid, "max_valid")
    if lower >= upper or not lower <= om <= upper or not lower <= pm <= upper:
        raise ValueError("prior and observed means must lie in the declared ordered physical domain")
    if psd <= MIN_PRIOR_SD_POSITIVITY_FLOOR:
        raise ValueError(
            f"prior_sd must be > positivity floor {MIN_PRIOR_SD_POSITIVITY_FLOOR}, got {psd} (AUDIT-REF)"
        )
    if ose < MIN_BIOLOGICAL_SE_FLOOR:
        raise ValueError(
            f"obs_se must be >= biological floor {MIN_BIOLOGICAL_SE_FLOOR}, got {ose} (AUDIT-REF)"
        )

    post_mean, post_sd, tau_prior, tau_obs, tau_post = _normal_normal_moments(pm, psd, om, ose)
    kl_nats = math.log(psd / post_sd) + 0.5 * ((post_sd / psd) ** 2 + ((post_mean - pm) / psd) ** 2 - 1.0)
    kl_bits = kl_nats / math.log(2.0)
    precision_gain_ratio = tau_post / tau_prior
    shannon_entropy_reduction_bits = 0.5 * math.log2(precision_gain_ratio)
    pred_sd = math.hypot(psd, ose)
    conflict_z = abs(om - pm) / pred_sd
    conflict_p = math.erfc(conflict_z / math.sqrt(2.0))
    # Student-t (df=5) likelihood weight attenuation under prior-data conflict (AUDIT-REF):
    t_weight_obs = tau_obs / (1.0 + (conflict_z ** 2) / 5.0)
    robust_tau = tau_prior + t_weight_obs
    robust_t_post_mean = _convex_mean(pm, om, tau_prior / robust_tau, t_weight_obs / robust_tau)
    diagnostics = (kl_nats, kl_bits, precision_gain_ratio, shannon_entropy_reduction_bits,
                   pred_sd, conflict_z, conflict_p, robust_t_post_mean)
    if not all(math.isfinite(value) for value in diagnostics):
        raise ValueError("Normal update diagnostics are outside representable numeric range")
    ci_low_raw = post_mean - 1.96 * post_sd
    ci_high_raw = post_mean + 1.96 * post_sd
    return {
        "status": "UPDATED_FROM_OBSERVED_READOUT",
        "unit_and_support": unit_and_support,
        "minimum_biological_se_floor_enforced": MIN_BIOLOGICAL_SE_FLOOR,
        "prior_mean": round(pm, 4),
        "prior_sd": _safe_round_positive_sd(psd),
        "prior_mean_exact": float(pm),
        "prior_sd_exact": float(psd),
        "observed_mean": round(om, 4),
        "observed_se": _safe_round_positive_sd(ose),
        "observed_mean_exact": float(om),
        "observed_se_exact": float(ose),
        "posterior_mean": round(post_mean, 4),
        "posterior_sd": _safe_round_positive_sd(post_sd),
        "posterior_mean_exact": float(post_mean),
        "posterior_sd_exact": float(post_sd),
        "posterior_precision_tau_exact": float(tau_post),
        "posterior_95ci_low": round(ci_low_raw, 4),
        "posterior_95ci_high": round(ci_high_raw, 4),
        "posterior_95ci_low_physical_clipped_v80_071": round(max(min_valid, ci_low_raw), 4),
        "posterior_95ci_high_physical_clipped_v80_071": round(min(max_valid, ci_high_raw), 4),
        "student_t_df5_conflict_robust_posterior_mean_v80_087": round(robust_t_post_mean, 4),
        "conflict_attenuated_heuristic_mean": round(robust_t_post_mean, 4),
        "interval_and_robustness_qualification": {
            "posterior": "UNBOUNDED_NORMAL_NORMAL_CONJUGATE_MODEL",
            "physical_clipped_fields": "DISPLAY_INTERSECTION_ONLY_NOT_A_95PCT_TRUNCATED_NORMAL_CREDIBLE_INTERVAL",
            "student_t_df5_legacy_field": "DEPRECATED_ALIAS_FOR_ONE_STEP_PRECISION_ATTENUATION_HEURISTIC_NOT_A_STUDENT_T_POSTERIOR",
        },
        "kl_divergence_post_vs_prior_nats": round(kl_nats, 4),
        "kl_divergence_post_vs_prior_bits": round(kl_bits, 4),
        "precision_gain_ratio": round(precision_gain_ratio, 4),
        "shannon_entropy_reduction_bits": round(shannon_entropy_reduction_bits, 4),
        "z_conflict": round(conflict_z, 4),
        "prior_data_conflict_z": round(conflict_z, 4),
        "prior_data_conflict_p_twosided": round(conflict_p, 4),
        "prior_data_conflict_flag_p_lt_0p05": bool(conflict_z >= 1.959964),
        "prior_data_conflict_detected": bool(conflict_z > PRIOR_DATA_CONFLICT_Z_THRESHOLD),
    }


def _calibration_state_sha256(report: Dict[str, Any]) -> str:
    """Bind all eight endpoint states and replay histories; integrity is not source authentication."""
    keys = ("calibration_state_version", "calibration_mode", "epistemic_lineage_status", "batch_id",
            "batch_payload_sha256", "ingested_batch_ids", "ingested_batch_payload_sha256s",
            "parent_report_sha256", "conjugate_posterior_updates", "optional_pouch_and_36run_screen_updates_v80_070")
    if any(k not in report for k in keys):
        raise ValueError("Incomplete calibration state; regenerate legacy reports from original observations")
    state_obj = {k: report[k] for k in keys}
    if "ordered_batch_observation_ledger" in report:
        state_obj["ordered_batch_observation_ledger"] = report["ordered_batch_observation_ledger"]
    if "stage0_decision_tree_evaluation" in report:
        state_obj["stage0_decision_tree_evaluation"] = report["stage0_decision_tree_evaluation"]
    for opt_lineage_key in ("contains_synthetic_records", "calibration_lineage_class"):
        if opt_lineage_key in report:
            state_obj[opt_lineage_key] = report[opt_lineage_key]
    return hashlib.sha256(json.dumps(state_obj, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _verify_endpoint_algebraic_integrity(ep_name: str, ep_obj: Dict[str, Any], pm: float, psd: float, *, min_valid: float, max_valid: float) -> None:
    """
    Verify that a prior report's `posterior_mean_exact` and `posterior_sd_exact` match the exact
    algebraic Normal-Normal conjugate update of its recorded `prior_*` and `observed_*` fields
    so tampered `posterior_mean_exact` values or forged status strings are rejected even if hashes
    are recomputed (AUDIT-REF, AUDIT-REF, & AUDIT-REF).
    """
    required_exact = {"prior_mean_exact", "prior_sd_exact", "observed_mean_exact", "observed_se_exact",
                      "posterior_mean_exact", "posterior_sd_exact", "posterior_precision_tau_exact"}
    if not required_exact <= ep_obj.keys():
        raise ValueError(f"Prior endpoint {ep_name!r} integrity requires complete exact machine state")
    status = ep_obj.get("status")
    if status == "UPDATED_FROM_OBSERVED_READOUT":
        pm0 = _require_finite_float(ep_obj["prior_mean_exact"], f"prior.{ep_name}.prior_mean_exact")
        psd0 = _require_finite_float(ep_obj["prior_sd_exact"], f"prior.{ep_name}.prior_sd_exact")
        om0 = _require_finite_float(ep_obj["observed_mean_exact"], f"prior.{ep_name}.observed_mean_exact")
        ose0 = _require_finite_float(ep_obj["observed_se_exact"], f"prior.{ep_name}.observed_se_exact")
        if not (min_valid <= pm0 <= max_valid and min_valid <= om0 <= max_valid):
            raise ValueError(f"Prior endpoint {ep_name!r} contains a prior/observed mean outside its declared physical domain")
        if psd0 <= MIN_PRIOR_SD_POSITIVITY_FLOOR or ose0 < MIN_BIOLOGICAL_SE_FLOOR:
            raise ValueError(f"Prior endpoint {ep_name!r} has invalid prior_sd={psd0} or observed_se={ose0}")
        expected = conjugate_normal_update(pm0, psd0, om0, ose0, min_valid=min_valid, max_valid=max_valid)
        expected.pop("unit_and_support")
    elif status == "UNAVAILABLE_PRIOR_RETAINED":
        pm0 = _require_finite_float(ep_obj["prior_mean_exact"], f"prior.{ep_name}.prior_mean_exact")
        psd0 = _require_finite_float(ep_obj["prior_sd_exact"], f"prior.{ep_name}.prior_sd_exact")
        if not min_valid <= pm0 <= max_valid:
            raise ValueError(f"Prior endpoint {ep_name!r} contains a prior mean outside its declared physical domain")
        if psd0 <= MIN_PRIOR_SD_POSITIVITY_FLOOR:
            raise ValueError(f"Prior endpoint {ep_name!r} has invalid retained prior_sd")
        expected = {
            "prior_mean": round(pm0, 4), "prior_sd": _safe_round_positive_sd(psd0),
            "posterior_mean_exact": pm0, "posterior_sd_exact": psd0,
            "posterior_mean": round(pm0, 4), "posterior_sd": _safe_round_positive_sd(psd0),
            "observed_mean": None, "observed_se": None, "observed_mean_exact": None, "observed_se_exact": None,
            "posterior_precision_tau_exact": (1.0 / psd0) ** 2,
            "posterior_95ci_low": round(pm0 - 1.96 * psd0, 4),
            "posterior_95ci_high": round(pm0 + 1.96 * psd0, 4),
            "kl_divergence_post_vs_prior_nats": 0.0, "kl_divergence_post_vs_prior_bits": 0.0,
            "precision_gain_ratio": 1.0, "shannon_entropy_reduction_bits": 0.0,
            "z_conflict": 0.0, "prior_data_conflict_z": 0.0, "prior_data_conflict_p_twosided": 1.0,
            "prior_data_conflict_flag_p_lt_0p05": False, "prior_data_conflict_detected": False,
        }
    else:
        raise ValueError(
            f"Prior report integrity violation on {ep_name!r}: invalid or unrecognized endpoint status {status!r}; "
            "must be 'UPDATED_FROM_OBSERVED_READOUT' or 'UNAVAILABLE_PRIOR_RETAINED' (AUDIT-REF)"
        )
    # Exact fields use relative arithmetic tolerance without an absolute SD budget:
    # 1e-7 absolute error can dominate small positive SDs and change sequential precision.
    for key, expected_value in expected.items():
        value = ep_obj.get(key)
        if isinstance(expected_value, bool):
            matches = type(value) is bool and value == expected_value
        elif isinstance(expected_value, (int, float)):
            numeric = _require_finite_float(value, f"prior.{ep_name}.{key}")
            matches = (math.isclose(numeric, expected_value, rel_tol=1e-12, abs_tol=0.0)
                       if key.endswith("_exact") else numeric == expected_value)
        else:
            matches = value == expected_value
        if not matches:
            raise ValueError(f"Prior endpoint {ep_name!r} integrity mismatch in {key!r}: {value!r} != {expected_value!r}")


# Single authoritative 27-flag persisted Stage-0 boolean registry shared across writer, helper replay, and production replay (R7S2-F02; DEC-057):
STAGE0_PERSISTED_BOOLEAN_FLAGS: tuple[str, ...] = (
    "require_full_5endpoint_s9_panel",
    "require_full_7metric_combined_stress_panel",
    "s1_gzmb_kass_or_pdl1_scar_failed",
    "s1_l3_bfb_replacement_qualified",
    "s2_arm1x_clears_15k_resting_hlae_and_cd47_c15g_noninferior",
    "s3_arm1b_and_arm1x_below_15k_and_arm1c_splicing_ge_0p95",
    "s4_rnls_wt_clears_gsis_shutoff_and_ocr_and_er_gate",
    "s5_fasl_delta_tunel_gt_5pp",
    "s5_sa_thbd_fhl1_only_replacement_qualified",
    "s6_er_rescue_genotype_requalified",
    "b_prime_sct_recipient_cd8_mlc_cleared",
    "s8_fasting_leak_or_ucn3_failed",
    "s8_staged_dose_and_omental_portal_drainage_qualified",
    "s9_mineral_pore_occlusion_gt_15pct",
    "s9_chronic_80um_lowload_or_omental_scaffold_qualified",
    "s9_chronic_80um_lowload_qualified",
    "s9_omental_revascularization_qualified",
    "obs_masked_perfusion_confirmed",
    "s10_sct_selectivity_coupon_cleared_kc_ka_ge_4_and_lysis_le_8pct",
    "combined_stress_coupon_cleared_lysis_le_8pct_and_si_ge_2p5",
    "track2_b_sct_or_hlag1_lilrb1_rescue_qualified",
    "track2_b_sct_or_hlag1_product_switched",
    "cd47_playbook_200k_and_5pct_phago_waiver_approved",
    "s9_staged_cassette_excision_and_d10w_clamp_qualified",
    "periimplant_closed_suction_drain_verified",
    "obs_vascularized_explant_en_bloc_recovery_verified",
    "obs_ap1903_pulsed_pk_washout_regimen_verified",
)

HISTORICAL_LEDGER_BOOLEAN_METADATA_KEYS: tuple[str, ...] = (
    *STAGE0_PERSISTED_BOOLEAN_FLAGS,
    "prior_batch_sticky_s9_safety_stop",
    "reset_prior_safety_stop_after_investigation",
    "allow_cross_clone_pooling_for_exploratory_research",
    "oos_retest_investigation_approved",
)

ALLOWED_ORDERED_BATCH_LEDGER_KEYS: frozenset[str] = frozenset({
    "batch_id",
    "batch_payload_sha256",
    "measurements",
    "defaulted_se_keys",
    "product_variant_id",
    "estimand_scale",
    "s8b_epistasis_estimand_scale",
    "obs_beta_x1_x4_epistasis_scale",
    "obs_beta_x3_x4_epistasis_scale",
    "specimen_or_clone_ids",
    "prior_variant_safety_status",
    "reset_prior_safety_stop_reason",
    "reset_prior_safety_stop_reason_code",
    "reset_prior_safety_stop_repeat_batch_ids",
    "reset_prior_safety_stop_independent_qp_approver_id",
    "mcb_clone_id",
    "oos_retest_assignable_cause_code",
    *HISTORICAL_LEDGER_BOOLEAN_METADATA_KEYS,
})


def _validate_prior_posterior_report(
    prior_report: Dict[str, Any],
    is_dry_run: bool,
    new_batch_id: str,
    new_batch_payload_sha256: str,
) -> tuple[list[str], list[str], Dict[str, tuple[float, float]]]:
    """
    Strictly validate `prior_posterior_report` schema, epistemic lineage, batch ID history,
    content-hash (`batch_payload_sha256`) history, cryptographic `conjugate_posterior_updates_sha256`,
    and algebraic integrity of all canonical endpoint posteriors before any branch executes
    (closing AUDIT-REF, AUDIT-REF, and AUDIT-REF).
    """
    if not isinstance(prior_report, dict):
        raise ValueError("prior_posterior_report must be a dictionary")
    if (type(prior_report.get("calibration_state_version")) is not int
            or prior_report["calibration_state_version"] != CALIBRATION_STATE_VERSION):
        raise ValueError("Prior report requires calibration_state_version=2; reconstruct legacy empirical history from original ordered observations")
    if prior_report.get("calibration_state_sha256") != _calibration_state_sha256(prior_report):
        raise ValueError("Prior report complete calibration-state digest mismatch")

    for req_key in ("calibration_mode", "ingested_batch_ids", "conjugate_posterior_updates"):
        if req_key not in prior_report:
            raise ValueError(f"prior_posterior_report missing required schema key {req_key!r}")

    prior_mode = prior_report["calibration_mode"]
    if not isinstance(prior_mode, str) or prior_mode not in VALID_CALIBRATION_MODES:
        raise ValueError(f"Invalid prior calibration_mode {prior_mode!r}")

    prior_lineage = prior_report.get("epistemic_lineage_status")
    if ((prior_mode == "EMPIRICAL_WETLAB_INGESTION" and prior_lineage not in SUPPORTED_EMPIRICAL_LINEAGES)
            or (prior_mode != "EMPIRICAL_WETLAB_INGESTION" and prior_lineage != "SYNTHETIC_DRY_RUN_LINEAGE")):
        raise ValueError("Prior calibration mode has inconsistent epistemic lineage")
    if "contains_synthetic_records" in prior_report:
        if type(prior_report["contains_synthetic_records"]) is not bool:
            raise ValueError("prior_posterior_report['contains_synthetic_records'] must be a boolean when provided (WETLAB-01)")
        if (prior_mode == "EMPIRICAL_WETLAB_INGESTION" and prior_report["contains_synthetic_records"] is True) or (
            prior_mode != "EMPIRICAL_WETLAB_INGESTION" and prior_report["contains_synthetic_records"] is False
        ):
            raise ValueError("Prior calibration mode has inconsistent contains_synthetic_records declaration (WETLAB-01)")
        if not is_dry_run and prior_report["contains_synthetic_records"] is not False:
            raise ValueError("Epistemic lineage violation: cannot chain prior report with contains_synthetic_records=True into empirical ingestion (WETLAB-01)")
    if "calibration_lineage_class" in prior_report:
        if not isinstance(prior_report["calibration_lineage_class"], str):
            raise ValueError("prior_posterior_report['calibration_lineage_class'] must be a string when provided (WETLAB-01)")
        if (prior_mode == "EMPIRICAL_WETLAB_INGESTION" and prior_report["calibration_lineage_class"] not in SUPPORTED_EMPIRICAL_LINEAGES) or (
            prior_mode != "EMPIRICAL_WETLAB_INGESTION" and prior_report["calibration_lineage_class"] != "SYNTHETIC_DRY_RUN_LINEAGE"
        ):
            raise ValueError("Prior calibration mode has inconsistent calibration_lineage_class declaration (WETLAB-01)")
    if not is_dry_run and (
        prior_mode != "EMPIRICAL_WETLAB_INGESTION"
        or prior_lineage not in SUPPORTED_EMPIRICAL_LINEAGES
    ):
        raise ValueError(
            f"Epistemic lineage violation: cannot chain prior report with calibration_mode={prior_mode!r}, "
            f"epistemic_lineage_status={prior_lineage!r} into an empirical wet-lab ingestion (AUDIT-REF)."
        )

    raw_batches = prior_report["ingested_batch_ids"]
    if not isinstance(raw_batches, list) or len(raw_batches) == 0:
        raise ValueError("prior_posterior_report['ingested_batch_ids'] must be a non-empty list of batch ID strings")
    norm_batches: list[str] = []
    for b in raw_batches:
        if not isinstance(b, str) or not b.strip() or b != b.strip() or len(b) > 120:
            raise ValueError(f"Invalid batch ID {b!r} in prior_posterior_report['ingested_batch_ids']")
        b_clean = b.strip()
        if not is_dry_run and b_clean.upper().startswith("DRY-RUN"):
            raise ValueError(
                f"Epistemic lineage violation: prior report contains synthetic batch_id {b_clean!r} in empirical mode"
            )
        if b_clean in norm_batches:
            raise ValueError(f"Duplicate batch ID {b_clean!r} inside prior_posterior_report['ingested_batch_ids']")
        norm_batches.append(b_clean)

    parent_digest = prior_report.get("parent_report_sha256")
    if ((len(norm_batches) == 1 and parent_digest is not None)
            or (len(norm_batches) > 1 and (not isinstance(parent_digest, str)
                or re.fullmatch(r"[0-9a-f]{64}", parent_digest) is None))):
        raise ValueError("Prior report parent lineage must be null at genesis and a SHA-256 for sequential batches")

    if new_batch_id in norm_batches:
        raise ValueError(f"Duplicate batch_id {new_batch_id!r} already ingested in prior posterior report!")

    # Validate content-level measurement payload hashes (closing AUDIT-REF):
    raw_prev_hashes = prior_report.get("ingested_batch_payload_sha256s")
    if (not isinstance(raw_prev_hashes, list) or len(raw_prev_hashes) != len(norm_batches)
            or any(not isinstance(h, str) or re.fullmatch(r"[0-9a-f]{64}", h) is None for h in raw_prev_hashes)
            or len(set(raw_prev_hashes)) != len(raw_prev_hashes)):
        raise ValueError("Prior payload digest history must contain one unique lowercase SHA-256 per batch")
    prev_payload_hashes = list(raw_prev_hashes)
    single_prev_hash = prior_report.get("batch_payload_sha256")
    if single_prev_hash != prev_payload_hashes[-1] or prior_report.get("batch_id") != norm_batches[-1]:
        raise ValueError("Latest calibration batch must match the final replay-history entry")

    if new_batch_payload_sha256 in prev_payload_hashes:
        raise ValueError(
            f"Duplicate measurement payload content replay detected: batch_payload_sha256={new_batch_payload_sha256} "
            f"(submitted under batch_id={new_batch_id!r}) matches an already-ingested batch (AUDIT-REF)!"
        )

    raw_updates = prior_report["conjugate_posterior_updates"]
    if not isinstance(raw_updates, dict) or set(raw_updates) != {row[0] for row in ENDPOINT_SPECS}:
        raise ValueError("prior report must contain exactly the canonical core endpoint states")

    # Verify cryptographic digest of conjugate_posterior_updates if present (or require it when exact fields exist; AUDIT-REF):
    expected_updates_sha256 = hashlib.sha256(
        json.dumps(raw_updates, sort_keys=True, allow_nan=False).encode("utf-8")
    ).hexdigest()
    recorded_updates_sha256 = prior_report.get("conjugate_posterior_updates_sha256")
    if recorded_updates_sha256 is None:
        raise ValueError("Prior report missing required conjugate_posterior_updates_sha256")
    if recorded_updates_sha256 is not None:
        if not isinstance(recorded_updates_sha256, str) or recorded_updates_sha256 != expected_updates_sha256:
            raise ValueError(
                f"Prior report integrity violation: conjugate_posterior_updates_sha256={recorded_updates_sha256!r} "
                f"does not match recomputed SHA-256 {expected_updates_sha256!r} (tampered prior_posterior_report rejected; AUDIT-REF)"
            )

    validated_priors: Dict[str, tuple[float, float]] = {}
    for ep_name, _def_m, _def_sd, _mk, _sk, min_valid, max_valid, _unit in ENDPOINT_SPECS:
        if ep_name not in raw_updates:
            raise ValueError(
                f"prior_posterior_report missing required canonical endpoint {ep_name!r} "
                "(silent reset to default prior is prohibited; AUDIT-REF)"
            )
        ep_obj = raw_updates[ep_name]
        if not isinstance(ep_obj, dict):
            raise ValueError(f"Prior endpoint {ep_name!r} must be a dictionary")
        raw_pm = ep_obj.get("posterior_mean_exact", ep_obj.get("posterior_mean"))
        raw_psd = ep_obj.get("posterior_sd_exact", ep_obj.get("posterior_sd"))
        pm = _require_finite_float(raw_pm, f"prior.{ep_name}.posterior_mean")
        psd = _require_finite_float(raw_psd, f"prior.{ep_name}.posterior_sd")
        if psd > 0.0 and (1.0 / psd) ** 2 == 0.0:
            raise ValueError(f"prior.{ep_name}.precision is outside representable numeric range")
        if psd <= MIN_PRIOR_SD_POSITIVITY_FLOOR:
            raise ValueError(
                f"prior.{ep_name}.posterior_sd must be > positivity floor {MIN_PRIOR_SD_POSITIVITY_FLOOR}, got {psd}"
            )
        if pm < min_valid or pm > max_valid:
            raise ValueError(
                f"prior.{ep_name}.posterior_mean={pm} outside valid physical domain [{min_valid}, {max_valid}]"
            )
        _verify_endpoint_algebraic_integrity(ep_name, ep_obj, pm, psd, min_valid=min_valid, max_valid=max_valid)
        validated_priors[ep_name] = (pm, psd)

    raw_ext_updates = prior_report.get("optional_pouch_and_36run_screen_updates_v80_070")
    expected_ext = {row[0] for row in EXTENDED_STAGE0_STAGE1_SPECS}
    if not isinstance(raw_ext_updates, dict) or set(raw_ext_updates) != expected_ext:
        raise ValueError("Prior report must preserve every extended endpoint, including unavailable priors")
    if any(not isinstance(v, dict) for v in raw_ext_updates.values()):
        raise ValueError("Every extended endpoint state must be an object")
    if isinstance(raw_ext_updates, dict):
        for ep_name, _def_m, _def_sd, _mk, _sk, _def_se, min_valid, max_valid, _unit in EXTENDED_STAGE0_STAGE1_SPECS:
            if ep_name in raw_ext_updates and isinstance(raw_ext_updates[ep_name], dict):
                ep_obj = raw_ext_updates[ep_name]
                raw_pm = ep_obj.get("posterior_mean_exact", ep_obj.get("posterior_mean"))
                raw_psd = ep_obj.get("posterior_sd_exact", ep_obj.get("posterior_sd"))
                pm = _require_finite_float(raw_pm, f"prior.{ep_name}.posterior_mean")
                psd = _require_finite_float(raw_psd, f"prior.{ep_name}.posterior_sd")
                if psd > 0.0 and (1.0 / psd) ** 2 == 0.0:
                    raise ValueError(f"prior.{ep_name}.precision is outside representable numeric range")
                if psd <= MIN_PRIOR_SD_POSITIVITY_FLOOR:
                    raise ValueError(f"prior.{ep_name}.posterior_sd must be > {MIN_PRIOR_SD_POSITIVITY_FLOOR}, got {psd}")
                if pm < min_valid or pm > max_valid:
                    raise ValueError(f"prior.{ep_name}.posterior_mean={pm} outside [{min_valid}, {max_valid}]")
                _verify_endpoint_algebraic_integrity(ep_name, ep_obj, pm, psd, min_valid=min_valid, max_valid=max_valid)
                validated_priors[ep_name] = (pm, psd)

    all_updates = [*raw_updates.values(), *raw_ext_updates.values()]
    actual_count = sum(ep["status"] == "UPDATED_FROM_OBSERVED_READOUT" for ep in all_updates)
    if (type(prior_report.get("observed_endpoints_updated_count")) is not int
            or prior_report["observed_endpoints_updated_count"] != actual_count):
        raise ValueError("Prior report observed endpoint count is inconsistent with endpoint state")
    actual_conflict = any(ep["prior_data_conflict_detected"] for ep in all_updates)
    if (type(prior_report.get("any_prior_data_conflict_detected_gt_3sigma")) is not bool
            or prior_report["any_prior_data_conflict_detected_gt_3sigma"] != actual_conflict):
        raise ValueError("Prior report conflict summary is inconsistent with endpoint state")
    latest_measurements = {}
    for specifications, endpoint_states in ((ENDPOINT_SPECS, raw_updates),
                                             (EXTENDED_STAGE0_STAGE1_SPECS, raw_ext_updates)):
        for row in specifications:
            ep_name, default_mean, default_sd, mean_key, se_key = row[:5]
            endpoint = endpoint_states[ep_name]
            if endpoint["status"] == "UPDATED_FROM_OBSERVED_READOUT":
                latest_measurements[mean_key] = endpoint["observed_mean_exact"]
                latest_measurements[se_key] = endpoint["observed_se_exact"]
            else:
                latest_measurements[mean_key] = "UNAVAILABLE"
            initial_precision = (1.0 / default_sd) ** 2
            prior_precision = (1.0 / endpoint["prior_sd_exact"]) ** 2
            # At most one legal independent observation per endpoint per batch;
            # no history can create more precision than n / minimum_SE**2.
            upper_prior_precision = initial_precision + (len(norm_batches) - 1) / MIN_BIOLOGICAL_SE_FLOOR**2
            if (prior_precision < initial_precision * (1.0 - 1e-12)
                    or prior_precision > upper_prior_precision * (1.0 + 1e-12)):
                raise ValueError(f"Prior endpoint {ep_name!r} precision is inconsistent with its batch history")
            # Added sufficient-statistic mass can only be a precision-weighted
            # mean of bounded observations. Zero added precision cannot move
            # the genesis mean. Centering the convex bounds avoids subtracting
            # nearly equal weighted means for genuinely weak observations.
            observation_weight = max(0.0, (prior_precision - initial_precision) / prior_precision)
            lower, upper = row[-3:-1]
            reachable_low = default_mean + observation_weight * (lower - default_mean)
            reachable_high = default_mean + observation_weight * (upper - default_mean)
            # A conservative per-batch floating-point budget for reciprocal,
            # square/root and convex-sum roundoff, scaled in ULPs rather than a
            # fixed physical mean/SD tolerance that could dominate weak data.
            scale = max(abs(lower), abs(upper), abs(default_mean))
            mean_roundoff = 32 * len(norm_batches) * math.ulp(scale)
            if not reachable_low - mean_roundoff <= endpoint["prior_mean_exact"] <= reachable_high + mean_roundoff:
                raise ValueError(f"Prior endpoint {ep_name!r} mean is outside its reachable bounded-observation history")
            if len(norm_batches) == 1 and (
                    not math.isclose(endpoint["prior_mean_exact"], default_mean, rel_tol=1e-12, abs_tol=0.0)
                    or not math.isclose(endpoint["prior_sd_exact"], default_sd, rel_tol=1e-12, abs_tol=0.0)):
                raise ValueError(f"Prior endpoint {ep_name!r} genesis prior differs from the declared initial prior")
            if "genesis_prior_mean_exact" in endpoint:
                gm = _require_finite_float(endpoint["genesis_prior_mean_exact"], f"prior.{ep_name}.genesis_prior_mean_exact")
                if not math.isclose(gm, float(default_mean), rel_tol=1e-12, abs_tol=0.0):
                    raise ValueError(
                        f"Prior endpoint {ep_name!r} genesis_prior_mean_exact={gm} does not match canonical initial prior mean {default_mean} (AUDIT-REF)"
                    )
            if "genesis_prior_sd_exact" in endpoint:
                gsd = _require_finite_float(endpoint["genesis_prior_sd_exact"], f"prior.{ep_name}.genesis_prior_sd_exact")
                if not math.isclose(gsd, float(default_sd), rel_tol=1e-12, abs_tol=0.0):
                    raise ValueError(
                        f"Prior endpoint {ep_name!r} genesis_prior_sd_exact={gsd} does not match canonical initial prior sd {default_sd} (AUDIT-REF)"
                    )
    # Reconstruct Stage-0 decision-tree observations from the prior report's stage0_decision_tree_evaluation / latest ledger entry (AUDIT-REF; DEC-040):
    global _MEASUREMENT_CONTRACT_MODULE
    if _MEASUREMENT_CONTRACT_MODULE is None:
        spec = importlib.util.spec_from_file_location("openislet_measurement_contract", REPO_ROOT / "scripts/wetlab_contract.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        _MEASUREMENT_CONTRACT_MODULE = helper
    primary_means_set = {r[3] for r in ENDPOINT_SPECS} | {r[3] for r in EXTENDED_STAGE0_STAGE1_SPECS}
    prior_stage0_eval = prior_report.get("stage0_decision_tree_evaluation")
    prior_validated_stage0 = (
        prior_stage0_eval.get("validated_stage0_endpoints", {})
        if isinstance(prior_stage0_eval, dict) and isinstance(prior_stage0_eval.get("validated_stage0_endpoints"), dict)
        else {}
    )
    prior_validated_stage0_se = (
        prior_stage0_eval.get("validated_stage0_se_endpoints", {})
        if isinstance(prior_stage0_eval, dict) and isinstance(prior_stage0_eval.get("validated_stage0_se_endpoints"), dict)
        else {}
    )
    prior_ledger_list = prior_report.get("ordered_batch_observation_ledger")
    latest_ledger_meas = (
        prior_ledger_list[-1].get("measurements", {})
        if isinstance(prior_ledger_list, list) and len(prior_ledger_list) > 0 and isinstance(prior_ledger_list[-1], dict) and isinstance(prior_ledger_list[-1].get("measurements"), dict)
        else {}
    )
    stage0_custom_defaults_hist = getattr(_MEASUREMENT_CONTRACT_MODULE, "STAGE0_CUSTOM_DEFAULT_SE", {})
    scale_aware_default_se_fn = getattr(
        _MEASUREMENT_CONTRACT_MODULE,
        "_scale_aware_stage0_default_se",
        lambda mk, lo_v, hi_v: stage0_custom_defaults_hist.get(mk, max(1e-3, (hi_v - lo_v) * 0.01)),
    )
    for _eid, mean_key, se_key, lo, hi, _unit in _get_stage0_decision_tree_specs():
        if mean_key in primary_means_set:
            continue
        if isinstance(prior_ledger_list, list) and len(prior_ledger_list) > 1:
            if mean_key in latest_ledger_meas and latest_ledger_meas[mean_key] not in ("UNAVAILABLE", None):
                m_val = _require_finite_float(latest_ledger_meas[mean_key], f"prior.stage0.{mean_key}")
                if mean_key in prior_validated_stage0 and not math.isclose(
                    _require_finite_float(prior_validated_stage0[mean_key], f"prior.stage0.{mean_key}"),
                    m_val,
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                ):
                    raise ValueError("Prior endpoint observations do not match the latest measurement payload digest")
                default_se = scale_aware_default_se_fn(mean_key, lo, hi)
                raw_se = latest_ledger_meas.get(se_key, prior_validated_stage0_se.get(se_key, default_se))
                s_val = _require_finite_float(raw_se, f"prior.stage0.{se_key}")
                latest_measurements[mean_key] = m_val
                latest_measurements[se_key] = s_val
        elif mean_key in prior_validated_stage0 and prior_validated_stage0[mean_key] not in ("UNAVAILABLE", None):
            m_val = _require_finite_float(prior_validated_stage0[mean_key], f"prior.stage0.{mean_key}")
            default_se = scale_aware_default_se_fn(mean_key, lo, hi)
            raw_se = prior_validated_stage0_se.get(se_key, latest_ledger_meas.get(se_key, default_se))
            s_val = _require_finite_float(raw_se, f"prior.stage0.{se_key}")
            latest_measurements[mean_key] = m_val
            latest_measurements[se_key] = s_val
    if compute_measurement_payload_sha256(latest_measurements) != single_prev_hash:
        raise ValueError("Prior endpoint observations do not match the latest measurement payload digest")
    if len(norm_batches) >= 1 or "ordered_batch_observation_ledger" in prior_report:
        obs_ledger = prior_report.get("ordered_batch_observation_ledger")
        if not isinstance(obs_ledger, list) or len(obs_ledger) != len(norm_batches):
            raise ValueError("Prior report requires a complete ordered_batch_observation_ledger matching its batch history")
        normalize_hist_measurements = _MEASUREMENT_CONTRACT_MODULE.normalize_measurements
        seen_specimen_keys: set[tuple[str, str]] = set()
        seen_specimen_endpoints: set[str] = set()
        seen_untagged_specimen_endpoints: set[str] = set()
        seen_product_variant_ids: set[str] = set()
        seen_estimand_scales: set[str] = set()
        for idx_b, entry in enumerate(obs_ledger):
            if (not isinstance(entry, dict)
                    or entry.get("batch_id") != norm_batches[idx_b]
                    or entry.get("batch_payload_sha256") != prev_payload_hashes[idx_b]
                    or not isinstance(entry.get("measurements"), dict)
                    or compute_measurement_payload_sha256(entry["measurements"]) != prev_payload_hashes[idx_b]):
                raise ValueError("ordered_batch_observation_ledger entry does not match ingested batch history or payload SHA-256")
            # Validate every historical measurement dict AND persisted boolean policy flag against the strict contract domain (AUDIT-REF, AUDIT-REF & R7S2-F02; DEC-057):
            non_canonical_hist_keys = set(entry["measurements"].keys()) - CANONICAL_MEASUREMENT_KEYS
            if non_canonical_hist_keys:
                raise ValueError(
                    f"ordered_batch_observation_ledger measurements must contain only canonical numerical mean/SE keys, not metadata keys {sorted(non_canonical_hist_keys)} (AUDIT-REF)"
                )
            if "specimen_or_clone_id" in entry:
                raise ValueError(
                    "ordered_batch_observation_ledger entries must normalize singular specimen_or_clone_id into specimen_or_clone_ids (AUDIT-REF)"
                )
            unknown_entry_keys = set(entry.keys()) - ALLOWED_ORDERED_BATCH_LEDGER_KEYS
            if unknown_entry_keys:
                raise ValueError(
                    f"ordered_batch_observation_ledger entry contains unsupported metadata keys {sorted(unknown_entry_keys)} (R7S2-F02; DEC-057)"
                )
            for bf_k in HISTORICAL_LEDGER_BOOLEAN_METADATA_KEYS:
                if bf_k in entry and entry[bf_k] is not None and type(entry[bf_k]) is not bool:
                    raise ValueError(
                        f"{bf_k} must be a boolean when provided in ordered_batch_observation_ledger (R7S2-F02; DEC-057)"
                    )
            merged_hist_payload_for_norm = dict(entry["measurements"])
            for bf_k in STAGE0_PERSISTED_BOOLEAN_FLAGS:
                if bf_k in entry and entry[bf_k] is not None:
                    merged_hist_payload_for_norm[bf_k] = entry[bf_k]
            normalize_hist_measurements(merged_hist_payload_for_norm, require_observation=False)
            obs_active_keys = {
                mk for mk, mv in entry["measurements"].items()
                if not mk.endswith("_se") and mv not in ("UNAVAILABLE", None)
            }
            if "product_variant_id" in entry and entry["product_variant_id"] is not None:
                if idx_b > 0 and not seen_product_variant_ids:
                    raise ValueError(
                        "Cannot introduce product_variant_id in ordered_batch_observation_ledger after genesis batch omitted it (AUDIT-REF)"
                    )
                pvid = entry["product_variant_id"]
                if not isinstance(pvid, str) or not pvid.strip() or pvid != pvid.strip() or len(pvid) > 120:
                    raise ValueError("ordered_batch_observation_ledger product_variant_id must be a canonical stripped nonempty string")
                from scripts.product_contract import get_product_definition
                canon_pvid = get_product_definition(pvid.strip())["product_variant_id"]
                seen_product_variant_ids.add(canon_pvid)
                if len(seen_product_variant_ids) > 1:
                    raise ValueError("Cross-batch product_variant_id mismatch inside ordered_batch_observation_ledger (AUDIT-REF)")
            elif seen_product_variant_ids:
                raise ValueError("Cannot omit product_variant_id in ordered_batch_observation_ledger after prior batch declared one (AUDIT-REF)")
            if "estimand_scale" in entry and entry["estimand_scale"] is not None:
                if idx_b > 0 and not seen_estimand_scales:
                    raise ValueError(
                        "Cannot introduce estimand_scale in ordered_batch_observation_ledger after genesis batch omitted it (AUDIT-REF)"
                    )
                esc = entry["estimand_scale"]
                if not isinstance(esc, str) or not esc.strip() or esc != esc.strip() or len(esc) > 120:
                    raise ValueError("ordered_batch_observation_ledger estimand_scale must be a canonical stripped nonempty string")
                seen_estimand_scales.add(esc.strip())
                if len(seen_estimand_scales) > 1:
                    raise ValueError("Cross-batch estimand_scale mismatch inside ordered_batch_observation_ledger (AUDIT-REF)")
            elif seen_estimand_scales:
                raise ValueError("Cannot omit estimand_scale in ordered_batch_observation_ledger after prior batch declared one (AUDIT-REF)")
            _MEASUREMENT_CONTRACT_MODULE.validate_product_and_estimand_compatibility(
                entry.get("product_variant_id"),
                entry.get("estimand_scale"),
                obs_active_keys,
                root=REPO_ROOT,
            )
            if "specimen_or_clone_id" in entry:
                raise ValueError(
                    "ordered_batch_observation_ledger entries must normalize singular specimen_or_clone_id into specimen_or_clone_ids (AUDIT-REF)"
                )
            if "specimen_or_clone_ids" in entry and entry["specimen_or_clone_ids"] is not None:
                smap = entry["specimen_or_clone_ids"]
                if not isinstance(smap, dict):
                    raise ValueError("ordered_batch_observation_ledger specimen_or_clone_ids must be a dictionary")
                for mk, sid in smap.items():
                    if not isinstance(mk, str) or mk not in obs_active_keys or not isinstance(sid, str) or not sid.strip() or sid != sid.strip() or len(sid) > 120:
                        raise ValueError("Invalid specimen_or_clone_ids entry in ordered_batch_observation_ledger")
                    if mk in seen_untagged_specimen_endpoints:
                        raise ValueError(
                            f"Cannot introduce specimen_or_clone_id for observed endpoint {mk!r} in ordered_batch_observation_ledger after prior batch observed {mk!r} without one (AUDIT-REF)"
                        )
                    pair = (mk, sid.strip())
                    if pair in seen_specimen_keys:
                        raise ValueError(
                            f"Duplicate specimen_or_clone_id {sid!r} for endpoint {mk!r} inside prior ordered_batch_observation_ledger (AUDIT-REF)"
                        )
                    seen_specimen_keys.add(pair)
                for req_mk in obs_active_keys & seen_specimen_endpoints:
                    if req_mk not in smap:
                        raise ValueError(
                            f"Cannot omit specimen_or_clone_id for observed endpoint {req_mk!r} in ordered_batch_observation_ledger after prior batch declared one (AUDIT-REF)"
                        )
                seen_specimen_endpoints.update(smap.keys())
                seen_untagged_specimen_endpoints.update(obs_active_keys - set(smap.keys()))
            elif obs_active_keys & seen_specimen_endpoints:
                raise ValueError(
                    "Cannot omit specimen_or_clone_ids in ordered_batch_observation_ledger after prior batch declared one for an observed endpoint (AUDIT-REF)"
                )
            else:
                seen_untagged_specimen_endpoints.update(obs_active_keys)
        for specifications, endpoint_states in ((ENDPOINT_SPECS, raw_updates),
                                                 (EXTENDED_STAGE0_STAGE1_SPECS, raw_ext_updates)):
            for row in specifications:
                ep_name, default_mean, default_sd, mean_key, se_key = row[:5]
                endpoint = endpoint_states[ep_name]
                rep_m, rep_sd = float(default_mean), float(default_sd)
                for entry in obs_ledger[:-1]:
                    m_dict = entry["measurements"]
                    if mean_key in m_dict and m_dict.get(mean_key) not in ("UNAVAILABLE", None):
                        om_i = _require_finite_float(m_dict[mean_key], mean_key)
                        os_i = _require_finite_float(m_dict[se_key], se_key)
                        rep_m, rep_sd, _, _, _ = _normal_normal_moments(rep_m, rep_sd, om_i, os_i)
                if (not math.isclose(endpoint["prior_mean_exact"], rep_m, rel_tol=1e-12, abs_tol=1e-12)
                        or not math.isclose(endpoint["prior_sd_exact"], rep_sd, rel_tol=1e-12, abs_tol=1e-12)):
                    raise ValueError(
                        f"Prior endpoint {ep_name!r} prior state is inconsistent with its ordered batch observation history (AUDIT-REF)"
                    )
                last_m_dict = obs_ledger[-1]["measurements"]
                post_m, post_sd = rep_m, rep_sd
                if mean_key in last_m_dict and last_m_dict.get(mean_key) not in ("UNAVAILABLE", None):
                    om_last = _require_finite_float(last_m_dict[mean_key], mean_key)
                    os_last = _require_finite_float(last_m_dict[se_key], se_key)
                    post_m, post_sd, _, _, _ = _normal_normal_moments(rep_m, rep_sd, om_last, os_last)
                if (not math.isclose(endpoint["posterior_mean_exact"], post_m, rel_tol=1e-12, abs_tol=1e-12)
                        or not math.isclose(endpoint["posterior_sd_exact"], post_sd, rel_tol=1e-12, abs_tol=1e-12)):
                    raise ValueError(
                        f"Prior endpoint {ep_name!r} posterior state is inconsistent with its ordered batch observation history (AUDIT-REF)"
                    )
        if isinstance(prior_stage0_eval, dict):
            allowed_stage0_eval_keys = {
                "schema_version",
                "strict_precedence_order",
                "mandatory_s9_safety_panel_complete",
                "full_5endpoint_s9_panel_complete",
                "s9_5endpoint_er_safety_panel_complete",
                "s9_causal_failed_categories",
                "s9_panel_completeness_tier",
                "lot_release_readiness_tier",
                "safety_kill_s9_triggered",
                "sticky_cross_batch_s9_stop_active",
                "primary_biological_or_pouch_gate_failed",
                "has_explicit_original_primary_obs",
                "has_explicit_release_critical_molecule_and_rate_se",
                "missing_original_primary_endpoints",
                "any_fallback_switch_triggered",
                "safety_kill_s9_status",
                "resolved_architecture_assignments",
                "validated_stage0_endpoints",
                "validated_stage0_se_endpoints",
                "stage0_se_provenance",
                "switch_decisions",
                "canonical_json_switch_id_map",
                "pivotal_combined_stress_coupon_gate",
                "prior_s9_safety_stop_cleared_after_investigation",
                "zero_event_proportion_missing_denominator_or_explicit_se",
            }
            unknown_stage0_eval_keys = sorted(set(prior_stage0_eval.keys()) - allowed_stage0_eval_keys)
            if unknown_stage0_eval_keys:
                raise ValueError(
                    f"Prior report stage0_decision_tree_evaluation contains unsupported keys {unknown_stage0_eval_keys} (ING-Q1; DEC-059)"
                )
            for s0_bool_k in (
                "safety_kill_s9_triggered",
                "mandatory_s9_safety_panel_complete",
                "full_5endpoint_s9_panel_complete",
                "s9_5endpoint_er_safety_panel_complete",
                "has_explicit_original_primary_obs",
                "has_explicit_release_critical_molecule_and_rate_se",
                "prior_s9_safety_stop_cleared_after_investigation",
                "zero_event_proportion_missing_denominator_or_explicit_se",
            ):
                if s0_bool_k in prior_stage0_eval and prior_stage0_eval[s0_bool_k] is not None and type(prior_stage0_eval[s0_bool_k]) is not bool:
                    raise ValueError(
                        f"Prior report stage0_decision_tree_evaluation.{s0_bool_k} must be a boolean when provided (R7S2-F02; DEC-057)"
                    )
            running_sticky_s9 = False
            stop_required_mcb_clonogenic_clearance = False
            stop_required_5endpoint_er_clearance = False
            stop_required_s9_categories: set[str] = set()
            post_stop_passing_s9_batch_ids: list[str] = []
            pooled_hist_for_s0: Dict[str, Any] = {}
            pooled_hist_def_se: set[str] = set()
            se_to_mean_hist = {
                **{
                    se_k: mean_k
                    for _eid, _pm0, _psd0, mean_k, se_k, _lo, _hi, _u in ENDPOINT_SPECS
                },
                **{
                    se_k: mean_k
                    for _eid, _pm0, _psd0, mean_k, se_k, _dse, _lo, _hi, _u in EXTENDED_STAGE0_STAGE1_SPECS
                },
                **{
                    se_k: mean_k
                    for _eid, mean_k, se_k, _lo, _hi, _u in _get_stage0_decision_tree_specs()
                },
            }
            helper_persisted_boolean_flags = STAGE0_PERSISTED_BOOLEAN_FLAGS
            for entry in obs_ledger:
                hist_meas_for_s0 = dict(entry["measurements"])
                for bf_k in helper_persisted_boolean_flags:
                    if bf_k in entry and entry[bf_k] is not None:
                        if type(entry[bf_k]) is not bool:
                            raise ValueError(
                                f"{bf_k} must be a boolean when provided in ordered_batch_observation_ledger (R7S2-F02; DEC-057)"
                            )
                        hist_meas_for_s0[bf_k] = entry[bf_k]
                hist_def_se = set(entry.get("defaulted_se_keys") or [])
                _drop_unpaired_rare_event_members(pooled_hist_for_s0, hist_meas_for_s0)
                for k, v in hist_meas_for_s0.items():
                    if v not in ("UNAVAILABLE", None):
                        pooled_hist_for_s0[k] = v
                        if k in se_to_mean_hist:
                            if k in hist_def_se or se_to_mean_hist[k] in hist_def_se:
                                pooled_hist_def_se.add(k)
                            else:
                                pooled_hist_def_se.discard(k)
                                pooled_hist_def_se.discard(se_to_mean_hist[k])
                for meta_k in ("product_variant_id", "estimand_scale"):
                    if entry.get(meta_k) is not None:
                        pooled_hist_for_s0[meta_k] = entry[meta_k]
                if entry.get("s8b_epistasis_estimand_scale") is not None:
                    pooled_hist_for_s0["s8b_epistasis_estimand_scale"] = entry["s8b_epistasis_estimand_scale"]
                else:
                    pooled_hist_for_s0.pop("s8b_epistasis_estimand_scale", None)
                if "obs_beta_x1_x4_epistasis" in hist_meas_for_s0 and hist_meas_for_s0.get("obs_beta_x1_x4_epistasis") not in ("UNAVAILABLE", None):
                    ep_sc_14 = entry.get("obs_beta_x1_x4_epistasis_scale") or entry.get("s8b_epistasis_estimand_scale")
                    if ep_sc_14 is not None:
                        pooled_hist_for_s0["obs_beta_x1_x4_epistasis_scale"] = ep_sc_14
                    else:
                        pooled_hist_for_s0.pop("obs_beta_x1_x4_epistasis_scale", None)
                if "obs_beta_x3_x4_epistasis" in hist_meas_for_s0 and hist_meas_for_s0.get("obs_beta_x3_x4_epistasis") not in ("UNAVAILABLE", None):
                    ep_sc_34 = entry.get("obs_beta_x3_x4_epistasis_scale") or entry.get("s8b_epistasis_estimand_scale")
                    if ep_sc_34 is not None:
                        pooled_hist_for_s0["obs_beta_x3_x4_epistasis_scale"] = ep_sc_34
                    else:
                        pooled_hist_for_s0.pop("obs_beta_x3_x4_epistasis_scale", None)
                single_batch_s0 = evaluate_stage0_decision_tree(
                    dict(hist_meas_for_s0),
                    defaulted_se_keys=sorted(hist_def_se),
                    prior_batch_sticky_s9_safety_stop=False,
                )
                if single_batch_s0.get("safety_kill_s9_triggered") is True:
                    post_stop_passing_s9_batch_ids.clear()
                    for cat in (single_batch_s0.get("s9_causal_failed_categories") or []):
                        stop_required_s9_categories.add(str(cat))
                    if (
                        single_batch_s0.get("resolved_architecture_assignments", {}).get("mcb_clonogenic_escape_qualification_status")
                        == "HARD_SAFETY_STOP_MCB_CLONOGENIC_SURVIVOR_DETECTED"
                        or "MCB_CLONOGENIC" in stop_required_s9_categories
                    ):
                        stop_required_mcb_clonogenic_clearance = True
                    if (
                        "ER_CHOP_HC10" in stop_required_s9_categories
                        or any(
                            "obs_chop_coincidence_fold" in r
                            or "obs_gate0e_chop_log2fc" in r
                            or "obs_hc10_surface_fold_wt" in r
                            or "obs_gate0e_hc10_surface_mfi_ratio" in r
                            for r in (single_batch_s0.get("triggered_rule_ids") or [])
                        )
                    ):
                        stop_required_5endpoint_er_clearance = True
                elif (
                    running_sticky_s9
                    and single_batch_s0.get("mandatory_s9_safety_panel_complete") is True
                    and single_batch_s0.get("safety_kill_s9_triggered") is False
                    and single_batch_s0.get("safety_kill_s9_status") == "CLEARED"
                    and (
                        "CORE_S9_VAF_OR_CASTSEQ" not in stop_required_s9_categories
                        or (
                            float(hist_meas_for_s0.get("obs_tp53_ppm1d_bcor_duplex_depth") or hist_meas_for_s0.get("obs_tp53_duplex_depth") or 0.0) >= 62958.0
                            and float(hist_meas_for_s0.get("obs_castseq_duplex_or_amplicon_depth") or 0.0) >= 62958.0
                        )
                    )
                    and (
                        not stop_required_mcb_clonogenic_clearance
                        or single_batch_s0.get("resolved_architecture_assignments", {}).get("mcb_clonogenic_escape_qualification_status")
                        == "CLEARED_MCB_CLONOGENIC_EFFECTIVE_TRIALS_GE_299573228"
                    )
                    and (
                        not stop_required_5endpoint_er_clearance
                        or single_batch_s0.get("s9_5endpoint_er_safety_panel_complete") is True
                        or single_batch_s0.get("full_5endpoint_s9_panel_complete") is True
                    )
                    and (
                        "RESIDUAL_PSC" not in stop_required_s9_categories
                        or (
                            hist_meas_for_s0.get("obs_residual_psc_effective_fraction") not in ("UNAVAILABLE", None)
                            and float(hist_meas_for_s0.get("obs_residual_psc_effective_trials") or 0.0) >= 299573228.0
                        )
                    )
                    and (
                        "CNV_20Q11" not in stop_required_s9_categories
                        or hist_meas_for_s0.get("obs_20q11_bcl2l1_cnv_ratio") not in ("UNAVAILABLE", None)
                    )
                    and (
                        "SINGLE_CELL_MOSAIC_20Q11" not in stop_required_s9_categories
                        or hist_meas_for_s0.get("obs_20q11_bcl2l1_single_cell_mosaic_fraction") not in ("UNAVAILABLE", None)
                    )
                    and (
                        "POST_STRESS_20Q11_OUTGROWTH" not in stop_required_s9_categories
                        or hist_meas_for_s0.get("obs_post_stress_20q11_bcl2l1_outgrowth_fold") not in ("UNAVAILABLE", None)
                    )
                    and (
                        "FUNCTIONAL_P53_NUTLIN3A" not in stop_required_s9_categories
                        or hist_meas_for_s0.get("obs_nutlin3a_p21_mdm2_fold_induction") not in ("UNAVAILABLE", None)
                    )
                    and (
                        "HAZARD_PROLIF_SURVIVOR" not in stop_required_s9_categories
                        or (
                            hist_meas_for_s0.get("obs_hazard_state_prolif_survivor_count_k") not in ("UNAVAILABLE", None)
                            and float(hist_meas_for_s0.get("obs_hazard_state_prolif_effective_trials") or 0.0) >= 299573228.0
                        )
                    )
                    and (
                        "PULSED_PK_KILL" not in stop_required_s9_categories
                        or (
                            hist_meas_for_s0.get("obs_ap1903_pulsed_pk_log10_kill") not in ("UNAVAILABLE", None)
                            and entry.get("obs_ap1903_pulsed_pk_washout_regimen_verified") is not False
                        )
                    )
                    and (
                        "MONOTHERAPY_PULSED_PK_KILL" not in stop_required_s9_categories
                        or (
                            hist_meas_for_s0.get("obs_ap1903_monotherapy_pulsed_pk_log10_kill") not in ("UNAVAILABLE", None)
                            and entry.get("obs_ap1903_pulsed_pk_washout_regimen_verified") is not False
                        )
                    )
                    and (
                        "STEP1_SUBLETHAL_SURVIVOR" not in stop_required_s9_categories
                        or (
                            hist_meas_for_s0.get("obs_ap1903_step1_sublethal_survivor_count_k") not in ("UNAVAILABLE", None)
                            and float(hist_meas_for_s0.get("obs_ap1903_step1_effective_trials") or 0.0) >= 29957322.8
                        )
                    )
                    and (
                        "EXPLANT_RECOVERY" not in stop_required_s9_categories
                        or (
                            entry.get("obs_vascularized_explant_en_bloc_recovery_verified") is True
                            and float(hist_meas_for_s0.get("obs_explant_en_bloc_recovery_fraction") or 0.0) >= 0.99
                        )
                    )
                    and isinstance(entry.get("batch_id"), str)
                    and entry["batch_id"].strip()
                ):
                    post_stop_passing_s9_batch_ids.append(entry["batch_id"].strip())
                step_s0_input = dict(pooled_hist_for_s0)
                if isinstance(entry.get("batch_id"), str) and entry["batch_id"].strip():
                    step_s0_input["batch_id"] = entry["batch_id"].strip()
                step_s0_input["verified_ledger_passing_s9_batch_ids"] = list(post_stop_passing_s9_batch_ids)
                for qp_k in (
                    "prior_batch_sticky_s9_safety_stop",
                    "prior_variant_safety_status",
                    "reset_prior_safety_stop_after_investigation",
                    "reset_prior_safety_stop_reason",
                    "reset_prior_safety_stop_reason_code",
                    "reset_prior_safety_stop_repeat_batch_ids",
                    "reset_prior_safety_stop_independent_qp_approver_id",
                ):
                    if entry.get(qp_k) is not None:
                        step_s0_input[qp_k] = entry[qp_k]
                hist_s0 = evaluate_stage0_decision_tree(
                    step_s0_input,
                    defaulted_se_keys=sorted(pooled_hist_def_se),
                    prior_batch_sticky_s9_safety_stop=running_sticky_s9,
                )
                running_sticky_s9 = bool(hist_s0.get("safety_kill_s9_triggered") is True) or bool(single_batch_s0.get("safety_kill_s9_triggered") is True)
                if not running_sticky_s9:
                    stop_required_mcb_clonogenic_clearance = False
                    stop_required_5endpoint_er_clearance = False
                    stop_required_s9_categories.clear()
            if running_sticky_s9 and prior_stage0_eval.get("safety_kill_s9_triggered") is not True:
                raise ValueError(
                    "Prior report stage0_decision_tree_evaluation.safety_kill_s9_triggered=False contradicts S9 hard safety stop in ordered_batch_observation_ledger (AUDIT-REF)"
                )
    return norm_batches, prev_payload_hashes, validated_priors


def current_model_reference_metadata() -> Dict[str, Dict[str, Any]]:
    """Derived comparisons are separate from independent target priors and assays."""
    metadata = {key: dict(value) for key, value in ENDPOINT_PROVENANCE_AND_MODEL_RECONCILIATION.items()}
    if not LAYER2_ATLAS_PATH.exists():
        return metadata
    atlas = loads_strict_json(LAYER2_ATLAS_PATH.read_text(encoding="utf-8"))
    try:
        rows = atlas["stage0_to_stage7_perturbation_plane"]["perturbation_profiles"]
        if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
            raise ValueError("Layer2 model profiles must be a list of objects")
        profiles = {}
        for row in rows:
            identifier = row.get("genotype_id")
            if not isinstance(identifier, str) or not identifier or identifier in profiles:
                raise ValueError("Layer2 model profiles require unique genotype identifiers")
            profiles[identifier] = row
        wt = profiles["WT_UNEDITED_CONTROL"]
        arm = profiles["SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS"]
    except (KeyError, TypeError) as exc:
        raise ValueError("Missing or malformed Layer2 model comparison profiles") from exc

    def model_number(profile: Dict[str, Any], key: str, label: str, maximum: float | None = None) -> float:
        try:
            value = _require_finite_float(profile[key], label)
        except (KeyError, OverflowError) as exc:
            raise ValueError(f"Invalid Layer2 model comparison value: {label}") from exc
        if value < 0.0 or (maximum is not None and value > maximum):
            raise ValueError(f"Layer2 model comparison value outside physical domain: {label}")
        return value

    wt_density = model_number(wt, "hour0_basal_cd47_molecules_per_cell", "WT CD47 density")
    arm_density = model_number(arm, "hour0_basal_cd47_molecules_per_cell", "Arm1B CD47 density")
    fold = arm_density / wt_density if wt_density > 0.0 else None
    cd = metadata["cd47_igv_hour0_fold_over_baseline_arm1b"]
    cd.update(layer2_or_layer3_model_reference_value=fold,
              layer2_arm1b_molecules_per_cell=arm_density,
              layer2_wt_molecules_per_cell=wt_density,
              model_reference_status=("CONDITIONAL_LAYER2_PRIOR_COMPARISON_NOT_MEASURED" if fold is not None
                                      else "UNDEFINED_ZERO_REFERENCE_CD47_DENSITY"))
    wt_yield = model_number(wt, "stage7_sc_beta_yield_fraction", "WT Stage7 yield", 1.0)
    arm_yield = model_number(arm, "stage7_sc_beta_yield_fraction", "Arm1B Stage7 yield", 1.0)
    yield_ratio = arm_yield / wt_yield if wt_yield > 0.0 else None
    metadata["stage7_sc_beta_yield_fraction_arm1b"].update(
        layer2_or_layer3_model_reference_value=yield_ratio,
        model_reference_status=("CONDITIONAL_LAYER2_PRIOR_COMPARISON_NOT_MEASURED" if yield_ratio is not None
                                else "UNDEFINED_ZERO_REFERENCE_STAGE7_YIELD"))
    return metadata


def evaluate_stage0_decision_tree(
    observations: Dict[str, Any],
    *,
    require_full_5endpoint_s9_panel: bool = False,
    defaulted_se_keys: Any = None,
    prior_batch_sticky_s9_safety_stop: bool = False,
    product_def: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Evaluate the machine-readable Stage-0 decision tree (`data/protocols/stage0_decision_tree_v1.json`, `S1-S9`) in strict precedence order (`DEC-036`, `DEC-039`, `DEC-040`, & `DEC-041`)."""
    global _MEASUREMENT_CONTRACT_MODULE
    if _MEASUREMENT_CONTRACT_MODULE is None:
        spec = importlib.util.spec_from_file_location("openislet_measurement_contract", REPO_ROOT / "scripts/wetlab_contract.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        _MEASUREMENT_CONTRACT_MODULE = helper
    return _MEASUREMENT_CONTRACT_MODULE.evaluate_stage0_decision_tree(
        observations,
        root=REPO_ROOT,
        require_full_5endpoint_s9_panel=require_full_5endpoint_s9_panel,
        defaulted_se_keys=defaulted_se_keys,
        prior_batch_sticky_s9_safety_stop=prior_batch_sticky_s9_safety_stop,
        product_def=product_def,
    )


def run_calibration(
    input_payload: Dict[str, Any],
    is_dry_run: bool,
    prior_posterior_report: Dict[str, Any] | None = None,
    fail_on_conflict_gt_3sigma: bool = False,
    *,
    require_full_5endpoint_s9_panel: bool = False,
) -> Dict[str, Any]:
    """
    Execute conjugate Bayesian calibration (closes AUDIT-REF, AUDIT-REF, AUDIT-REF,
    AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF,
    AUDIT-REF & 07, and AUDIT-REF / EXT-A).
    """
    if not isinstance(input_payload, dict) or not input_payload:
        raise ValueError("input_payload must be a non-empty dictionary")
    if not isinstance(is_dry_run, bool) or not isinstance(fail_on_conflict_gt_3sigma, bool) or not isinstance(require_full_5endpoint_s9_panel, bool):
        raise ValueError("Calibration mode, conflict policy, and S9 panel policy must be booleans")

    global _MEASUREMENT_CONTRACT_MODULE
    if _MEASUREMENT_CONTRACT_MODULE is None:
        spec = importlib.util.spec_from_file_location("openislet_measurement_contract", REPO_ROOT / "scripts/wetlab_contract.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        _MEASUREMENT_CONTRACT_MODULE = helper
    normalize_measurements = _MEASUREMENT_CONTRACT_MODULE.normalize_measurements
    received_payload_sha256 = hashlib.sha256(json.dumps(input_payload, sort_keys=True, allow_nan=False).encode("utf-8")).hexdigest()
    input_payload, measurement_validation = normalize_measurements(input_payload, require_observation=not is_dry_run)
    if input_payload.get("require_full_5endpoint_s9_panel") is True:
        require_full_5endpoint_s9_panel = True
    batch_id = input_payload.get("batch_id", input_payload.get("assay_id"))
    if not is_dry_run:
        if not isinstance(batch_id, str) or not batch_id.strip():
            raise ValueError("Empirical wet-lab payload must include a non-empty string 'batch_id' or 'assay_id'.")
        batch_id = batch_id.strip()
        if batch_id.upper().startswith("DRY-RUN"):
            raise ValueError(f"Empirical wet-lab payload cannot use a DRY-RUN batch_id ({batch_id!r}).")
    else:
        if batch_id is not None and (not isinstance(batch_id, str) or not batch_id.strip()):
            raise ValueError("batch_id must be a non-empty string when provided")
        batch_id = str(batch_id or "DRY-RUN-CALIBRATION-001").strip()

    batch_payload_sha256 = compute_measurement_payload_sha256(input_payload)

    prev_batches: list[str] = []
    prev_payload_hashes: list[str] = []
    validated_priors: Dict[str, tuple[float, float]] = {}
    parent_report_sha256: str | None = None
    if prior_posterior_report is not None:
        prev_batches, prev_payload_hashes, validated_priors = _validate_prior_posterior_report(
            prior_posterior_report,
            is_dry_run=is_dry_run,
            new_batch_id=batch_id,
            new_batch_payload_sha256=batch_payload_sha256,
        )
        parent_report_sha256 = hashlib.sha256(
            json.dumps(prior_posterior_report, sort_keys=True, allow_nan=False).encode("utf-8")
        ).hexdigest()

    has_any_extended_stage0_stage1_obs = any(
        (mean_key in input_payload and input_payload.get(mean_key) not in ("UNAVAILABLE", None))
        for _ep, _dm, _dsd, mean_key, _sk, _dse, _minv, _maxv, _u in EXTENDED_STAGE0_STAGE1_SPECS
    ) or any(
        (mean_key in input_payload and input_payload.get(mean_key) not in ("UNAVAILABLE", None))
        for _eid, mean_key, _sk, _lo, _hi, _u in _get_stage0_decision_tree_specs()
    )

    model_reference_metadata = current_model_reference_metadata()
    updates: Dict[str, Any] = {}
    observed_count = 0

    for ep_name, default_prior_m, default_prior_sd, mean_key, se_key, min_valid, max_valid, unit_desc in ENDPOINT_SPECS:
        if ep_name in validated_priors:
            prior_m, prior_sd = validated_priors[ep_name]
        else:
            prior_m, prior_sd = float(default_prior_m), float(default_prior_sd)

        prov_meta = model_reference_metadata.get(ep_name, {})

        # Allow pouch-only or Stage-0/Stage-1 partial batches to omit core keys when at least one extended endpoint is present
        if input_payload.get(mean_key) == "UNAVAILABLE" or (mean_key not in input_payload and has_any_extended_stage0_stage1_obs):
            updates[ep_name] = {
                "status": "UNAVAILABLE_PRIOR_RETAINED",
                "unit_and_support": unit_desc,
                "minimum_biological_se_floor_enforced": MIN_BIOLOGICAL_SE_FLOOR,
                "prior_provenance": prov_meta.get("prior_source", ""),
                "layer2_or_layer3_model_reference_value": prov_meta.get("layer2_or_layer3_model_reference_value"),
                "mapped_ranker_prior_key": prov_meta.get("mapped_ranker_prior_key"),
                "prior_mean": round(prior_m, 4),
                "prior_sd": _safe_round_positive_sd(prior_sd),
                "prior_mean_exact": float(prior_m),
                "prior_sd_exact": float(prior_sd),
                "observed_mean": None,
                "observed_se": None,
                "observed_mean_exact": None,
                "observed_se_exact": None,
                "posterior_mean": round(prior_m, 4),
                "posterior_sd": _safe_round_positive_sd(prior_sd),
                "posterior_mean_exact": float(prior_m),
                "posterior_sd_exact": float(prior_sd),
                "posterior_precision_tau_exact": float((1.0 / prior_sd) ** 2),
                "posterior_95ci_low": round(prior_m - 1.96 * prior_sd, 4),
                "posterior_95ci_high": round(prior_m + 1.96 * prior_sd, 4),
                "kl_divergence_post_vs_prior_nats": 0.0,
                "kl_divergence_post_vs_prior_bits": 0.0,
                "precision_gain_ratio": 1.0,
                "shannon_entropy_reduction_bits": 0.0,
                "z_conflict": 0.0,
                "prior_data_conflict_z": 0.0,
                "prior_data_conflict_p_twosided": 1.0,
                "prior_data_conflict_flag_p_lt_0p05": False,
                "prior_data_conflict_detected": False,
            }
            continue

        if mean_key not in input_payload or se_key not in input_payload:
            raise ValueError(
                f"Missing required empirical readout keys ({mean_key!r}, {se_key!r}) for endpoint {ep_name!r}; "
                "pass explicit 'UNAVAILABLE' if an assay endpoint was not yet acquired."
            )

        obs_m = _require_finite_float(input_payload[mean_key], mean_key)
        obs_s = _require_finite_float(input_payload[se_key], se_key)
        # Allow valid 0.0 zero-effect / assay-failure readouts on closed interval [min_valid, max_valid] (AUDIT-REF):
        if obs_m < min_valid or obs_m > max_valid:
            raise ValueError(
                f"Out-of-range physical observation for {mean_key}: {obs_m} (must be in [{min_valid}, {max_valid}])"
            )
        ep_update = conjugate_normal_update(
            prior_mean=prior_m,
            prior_sd=prior_sd,
            obs_mean=obs_m,
            obs_se=obs_s,
            unit_and_support=unit_desc,
            min_valid=min_valid,
            max_valid=max_valid,
        )
        ep_update["prior_provenance"] = prov_meta.get("prior_source", "")
        ep_update["layer2_or_layer3_model_reference_value"] = prov_meta.get("layer2_or_layer3_model_reference_value")
        ep_update["mapped_ranker_prior_key"] = prov_meta.get("mapped_ranker_prior_key")
        if fail_on_conflict_gt_3sigma and ep_update["prior_data_conflict_detected"]:
            raise ValueError(
                f"Prior-data conflict z={ep_update['z_conflict']} exceeds 3.0 sigma threshold on {ep_name!r}"
            )
        updates[ep_name] = ep_update
        observed_count += 1

    # Stage-0/Stage-1 Layer-4 pouch oxygen/peroxide/pH & 48h IFN-g HLA-DR/DQ negative % ingestion (closing AUDIT-REF & AUDIT-REF):
    optional_pouch_and_screen_updates: Dict[str, Any] = {}
    for ep_name, default_prior_m, default_prior_sd, mean_key, se_key, default_se, min_valid, max_valid, unit_desc in EXTENDED_STAGE0_STAGE1_SPECS:
        prior_m, prior_sd = validated_priors.get(ep_name, (float(default_prior_m), float(default_prior_sd)))
        if mean_key in input_payload and input_payload[mean_key] not in ("UNAVAILABLE", None):
            obs_val = _require_finite_float(input_payload[mean_key], mean_key)
            obs_se_val = _require_finite_float(input_payload.get(se_key, default_se), se_key)
            if obs_val < min_valid or obs_val > max_valid:
                raise ValueError(f"{mean_key}={obs_val} outside [{min_valid}, {max_valid}]")
            ext_up = conjugate_normal_update(
                prior_mean=prior_m,
                prior_sd=prior_sd,
                obs_mean=obs_val,
                obs_se=obs_se_val,
                unit_and_support=unit_desc,
                min_valid=min_valid,
                max_valid=max_valid,
            )
            if fail_on_conflict_gt_3sigma and ext_up["prior_data_conflict_detected"]:
                raise ValueError(
                    f"Prior-data conflict z={ext_up['z_conflict']} exceeds 3.0 sigma threshold on {ep_name!r}"
                )
            optional_pouch_and_screen_updates[ep_name] = ext_up
            observed_count += 1
        else:
            optional_pouch_and_screen_updates[ep_name] = {
                "status": "UNAVAILABLE_PRIOR_RETAINED",
                "unit_and_support": unit_desc,
                "minimum_biological_se_floor_enforced": MIN_BIOLOGICAL_SE_FLOOR,
                "prior_mean": round(prior_m, 4),
                "prior_sd": _safe_round_positive_sd(prior_sd),
                "prior_mean_exact": float(prior_m),
                "prior_sd_exact": float(prior_sd),
                "observed_mean": None,
                "observed_se": None,
                "observed_mean_exact": None,
                "observed_se_exact": None,
                "posterior_mean": round(prior_m, 4),
                "posterior_sd": _safe_round_positive_sd(prior_sd),
                "posterior_mean_exact": float(prior_m),
                "posterior_sd_exact": float(prior_sd),
                "posterior_precision_tau_exact": float((1.0 / prior_sd) ** 2),
                "posterior_95ci_low": round(prior_m - 1.96 * prior_sd, 4),
                "posterior_95ci_high": round(prior_m + 1.96 * prior_sd, 4),
                "kl_divergence_post_vs_prior_nats": 0.0,
                "kl_divergence_post_vs_prior_bits": 0.0,
                "precision_gain_ratio": 1.0,
                "shannon_entropy_reduction_bits": 0.0,
                "z_conflict": 0.0,
                "prior_data_conflict_z": 0.0,
                "prior_data_conflict_p_twosided": 1.0,
                "prior_data_conflict_flag_p_lt_0p05": False,
                "prior_data_conflict_detected": False,
            }

    all_candidate_mean_keys = list(dict.fromkeys(
        [r[3] for r in ENDPOINT_SPECS]
        + [r[3] for r in EXTENDED_STAGE0_STAGE1_SPECS]
        + [r[1] for r in _get_stage0_decision_tree_specs()]
    ))
    all_observed_mean_keys = [
        mk for mk in all_candidate_mean_keys
        if mk in input_payload and input_payload.get(mk) not in ("UNAVAILABLE", None)
    ]
    if not is_dry_run and observed_count == 0 and not all_observed_mean_keys:
        raise ValueError("Empirical wet-lab ingestion requires at least one observed numeric endpoint.")

    raw_bytes = json.dumps(input_payload, sort_keys=True, allow_nan=False).encode("utf-8")
    updates_bytes = json.dumps(updates, sort_keys=True, allow_nan=False).encode("utf-8")
    prev_obs_ledger = (
        list(prior_posterior_report.get("ordered_batch_observation_ledger", []))
        if prior_posterior_report is not None
        else []
    )
    batch_specimen_map: Dict[str, str] = {}
    if "specimen_or_clone_id" in input_payload and input_payload["specimen_or_clone_id"] is not None:
        raw_sid = input_payload["specimen_or_clone_id"]
        if not isinstance(raw_sid, str) or not raw_sid.strip():
            raise ValueError("specimen_or_clone_id must be a nonempty string when provided")
        for mk in all_observed_mean_keys:
            batch_specimen_map[mk] = raw_sid.strip()
    if "specimen_or_clone_ids" in input_payload and input_payload["specimen_or_clone_ids"] is not None:
        raw_smap = input_payload["specimen_or_clone_ids"]
        if not isinstance(raw_smap, dict):
            raise ValueError("specimen_or_clone_ids must be a dictionary when provided")
        for mk, sid in raw_smap.items():
            if not isinstance(mk, str) or mk not in all_observed_mean_keys or not isinstance(sid, str) or not sid.strip():
                raise ValueError("specimen_or_clone_ids keys must be observed endpoints and values must be nonempty strings")
            if mk in batch_specimen_map and batch_specimen_map[mk] != sid.strip():
                raise ValueError(
                    f"Conflicting specimen_or_clone_id ({batch_specimen_map[mk]!r}) and specimen_or_clone_ids[{mk!r}] ({sid.strip()!r}) in measurement payload (AUDIT-REF)"
                )
            batch_specimen_map[mk] = sid.strip()
    for prev_entry in prev_obs_ledger:
        if prev_entry.get("product_variant_id") is not None:
            if "product_variant_id" not in input_payload or input_payload["product_variant_id"] is None:
                raise ValueError(
                    f"Cross-batch product_variant_id omission: prior batch declared {prev_entry['product_variant_id']!r} "
                    f"but current batch omitted product_variant_id (AUDIT-REF)"
                )
            from scripts.product_contract import get_product_definition
            prev_canon_pvid = get_product_definition(str(prev_entry["product_variant_id"]).strip())["product_variant_id"]
            curr_canon_pvid = get_product_definition(str(input_payload["product_variant_id"]).strip())["product_variant_id"]
            if prev_canon_pvid != curr_canon_pvid:
                raise ValueError(
                    f"Cross-batch product_variant_id mismatch: prior batch used {prev_entry['product_variant_id']!r} "
                    f"vs current batch {input_payload['product_variant_id']!r} (AUDIT-REF)"
                )
        elif "product_variant_id" in input_payload and input_payload["product_variant_id"] is not None:
            raise ValueError(
                f"Cross-batch product_variant_id genesis omission: prior batch {prev_entry.get('batch_id')!r} omitted "
                f"product_variant_id but current batch declared {input_payload['product_variant_id']!r} (AUDIT-REF)"
            )
        if prev_entry.get("estimand_scale") is not None:
            if "estimand_scale" not in input_payload or input_payload["estimand_scale"] is None:
                raise ValueError(
                    f"Cross-batch estimand_scale omission: prior batch declared {prev_entry['estimand_scale']!r} "
                    f"but current batch omitted estimand_scale (AUDIT-REF)"
                )
            if prev_entry["estimand_scale"] != input_payload["estimand_scale"]:
                raise ValueError(
                    f"Cross-batch estimand_scale mismatch: prior batch used {prev_entry['estimand_scale']!r} "
                    f"vs current batch {input_payload['estimand_scale']!r} (AUDIT-REF)"
                )
        elif "estimand_scale" in input_payload and input_payload["estimand_scale"] is not None:
            raise ValueError(
                f"Cross-batch estimand_scale genesis omission: prior batch {prev_entry.get('batch_id')!r} omitted "
                f"estimand_scale but current batch declared {input_payload['estimand_scale']!r} (AUDIT-REF)"
            )
        prev_smap = prev_entry.get("specimen_or_clone_ids") or {}
        prev_meas = prev_entry.get("measurements") or {}
        for mk in all_observed_mean_keys:
            if mk in prev_smap and mk not in batch_specimen_map:
                raise ValueError(
                    f"Cross-batch specimen_or_clone_id omission for endpoint {mk!r}: prior batch declared {prev_smap[mk]!r} "
                    f"but current batch omitted specimen_or_clone_id (AUDIT-REF)"
                )
        for mk, sid in batch_specimen_map.items():
            if prev_meas.get(mk) not in ("UNAVAILABLE", None) and mk not in prev_smap:
                raise ValueError(
                    f"Cross-batch specimen_or_clone_id genesis omission for endpoint {mk!r}: prior batch {prev_entry.get('batch_id')!r} "
                    f"observed {mk!r} without specimen_or_clone_id but current batch declared {sid!r} (AUDIT-REF)"
                )
            if prev_smap.get(mk) == sid:
                raise ValueError(
                    f"Duplicate specimen_or_clone_id {sid!r} for endpoint {mk!r} across sequential batches "
                    f"(pseudoreplication rejected; AUDIT-REF)"
                )
    current_obs_entry: Dict[str, Any] = {
        "batch_id": batch_id,
        "batch_payload_sha256": batch_payload_sha256,
        "measurements": _extract_canonical_measurements(input_payload),
    }
    if measurement_validation.get("defaulted_standard_error_keys"):
        current_obs_entry["defaulted_se_keys"] = list(measurement_validation["defaulted_standard_error_keys"])
    if "product_variant_id" in input_payload:
        current_obs_entry["product_variant_id"] = input_payload["product_variant_id"]
    if "estimand_scale" in input_payload:
        current_obs_entry["estimand_scale"] = input_payload["estimand_scale"]
    if "s8b_epistasis_estimand_scale" in input_payload:
        current_obs_entry["s8b_epistasis_estimand_scale"] = input_payload["s8b_epistasis_estimand_scale"]
    if "obs_beta_x1_x4_epistasis" in input_payload and input_payload.get("obs_beta_x1_x4_epistasis") not in ("UNAVAILABLE", None):
        sc_14_curr = input_payload.get("obs_beta_x1_x4_epistasis_scale") or input_payload.get("s8b_epistasis_estimand_scale")
        if sc_14_curr is not None:
            current_obs_entry["obs_beta_x1_x4_epistasis_scale"] = sc_14_curr
    if "obs_beta_x3_x4_epistasis" in input_payload and input_payload.get("obs_beta_x3_x4_epistasis") not in ("UNAVAILABLE", None):
        sc_34_curr = input_payload.get("obs_beta_x3_x4_epistasis_scale") or input_payload.get("s8b_epistasis_estimand_scale")
        if sc_34_curr is not None:
            current_obs_entry["obs_beta_x3_x4_epistasis_scale"] = sc_34_curr
    stage0_persisted_boolean_flags = STAGE0_PERSISTED_BOOLEAN_FLAGS
    for qp_meta_k in (
        "prior_batch_sticky_s9_safety_stop",
        "prior_variant_safety_status",
        "reset_prior_safety_stop_after_investigation",
        "reset_prior_safety_stop_reason",
        "reset_prior_safety_stop_reason_code",
        "reset_prior_safety_stop_repeat_batch_ids",
        "reset_prior_safety_stop_independent_qp_approver_id",
        "mcb_clone_id",
        "allow_cross_clone_pooling_for_exploratory_research",
        "oos_retest_investigation_approved",
        "oos_retest_assignable_cause_code",
        *stage0_persisted_boolean_flags,
    ):
        if qp_meta_k in input_payload and input_payload[qp_meta_k] is not None:
            current_obs_entry[qp_meta_k] = input_payload[qp_meta_k]
    if require_full_5endpoint_s9_panel is True:
        current_obs_entry["require_full_5endpoint_s9_panel"] = True
    if batch_specimen_map:
        current_obs_entry["specimen_or_clone_ids"] = batch_specimen_map

    # Pool historical Stage-0 endpoint observations across chained batches and latch any prior S9 hard safety stop, cross-clone pooling mismatch, or uninvestigated non-S9 OOS overwrite (AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF; DEC-041/043/044/047/050):
    pooled_stage0_payload: Dict[str, Any] = {}
    pooled_defaulted_se_set: set[str] = set()
    prior_sticky_s9 = False
    stop_required_mcb_clonogenic_clearance = False
    stop_required_5endpoint_er_clearance = False
    stop_required_s9_categories: set[str] = set()
    prior_s9_history_retained_across_ledger = False
    post_stop_passing_s9_batch_ids: list[str] = []
    seen_mcb_clone_ids: set[str] = set()
    any_batch_omitted_mcb_clone_id = False
    prior_failed_stage0_gate_keys: set[str] = set()
    uninvestigated_non_s9_oos_overwrite_hold = False
    s9_3_endpoint_keys = {
        "obs_icasp9_log10_kill",
        "obs_icasp9_log10_kill_se",
        "obs_tp53_ppm1d_bcor_vaf",
        "obs_tp53_ppm1d_bcor_vaf_se",
        "obs_castseq_translocation_freq",
        "obs_castseq_translocation_se",
    }
    se_to_mean_lookup = {
        **{
            se_k: mean_k
            for _eid, _pm0, _psd0, mean_k, se_k, _lo, _hi, _u in ENDPOINT_SPECS
        },
        **{
            se_k: mean_k
            for _eid, _pm0, _psd0, mean_k, se_k, _dse, _lo, _hi, _u in EXTENDED_STAGE0_STAGE1_SPECS
        },
        **{
            se_k: mean_k
            for _eid, mean_k, se_k, _lo, _hi, _u in _get_stage0_decision_tree_specs()
        },
    }
    for hist_entry in prev_obs_ledger:
        if isinstance(hist_entry.get("mcb_clone_id"), str) and hist_entry["mcb_clone_id"].strip():
            seen_mcb_clone_ids.add(hist_entry["mcb_clone_id"].strip())
        else:
            any_batch_omitted_mcb_clone_id = True
        hist_meas = dict(hist_entry.get("measurements") or {})
        for bf_k in stage0_persisted_boolean_flags:
            if bf_k in hist_entry and hist_entry[bf_k] is not None:
                if type(hist_entry[bf_k]) is not bool:
                    raise ValueError(
                        f"{bf_k} must be a boolean when provided in ordered_batch_observation_ledger (R7S2-F02; DEC-057)"
                    )
                hist_meas[bf_k] = hist_entry[bf_k]
        hist_def_se = set(hist_entry.get("defaulted_se_keys") or [])
        hist_stage0_noncore_keys = {
            k for k, v in hist_meas.items()
            if k not in s9_3_endpoint_keys and not k.endswith("_se") and v not in ("UNAVAILABLE", None)
        }
        if hist_stage0_noncore_keys & prior_failed_stage0_gate_keys:
            if not (
                hist_entry.get("oos_retest_investigation_approved") is True
                and isinstance(hist_entry.get("oos_retest_assignable_cause_code"), str)
                and hist_entry["oos_retest_assignable_cause_code"].strip()
            ):
                uninvestigated_non_s9_oos_overwrite_hold = True
            else:
                prior_failed_stage0_gate_keys -= hist_stage0_noncore_keys
        _drop_unpaired_rare_event_members(pooled_stage0_payload, hist_meas)
        for k, v in hist_meas.items():
            if v not in ("UNAVAILABLE", None):
                pooled_stage0_payload[k] = v
                if k in se_to_mean_lookup:
                    if k in hist_def_se or se_to_mean_lookup[k] in hist_def_se:
                        pooled_defaulted_se_set.add(k)
                    else:
                        pooled_defaulted_se_set.discard(k)
                        pooled_defaulted_se_set.discard(se_to_mean_lookup[k])
        for meta_k in ("product_variant_id", "estimand_scale"):
            if hist_entry.get(meta_k) is not None:
                pooled_stage0_payload[meta_k] = hist_entry[meta_k]
        if hist_entry.get("s8b_epistasis_estimand_scale") is not None:
            pooled_stage0_payload["s8b_epistasis_estimand_scale"] = hist_entry["s8b_epistasis_estimand_scale"]
        else:
            pooled_stage0_payload.pop("s8b_epistasis_estimand_scale", None)
        if "obs_beta_x1_x4_epistasis" in hist_meas and hist_meas.get("obs_beta_x1_x4_epistasis") not in ("UNAVAILABLE", None):
            ep_sc_14 = hist_entry.get("obs_beta_x1_x4_epistasis_scale") or hist_entry.get("s8b_epistasis_estimand_scale")
            if ep_sc_14 is not None:
                pooled_stage0_payload["obs_beta_x1_x4_epistasis_scale"] = ep_sc_14
            else:
                pooled_stage0_payload.pop("obs_beta_x1_x4_epistasis_scale", None)
        if "obs_beta_x3_x4_epistasis" in hist_meas and hist_meas.get("obs_beta_x3_x4_epistasis") not in ("UNAVAILABLE", None):
            ep_sc_34 = hist_entry.get("obs_beta_x3_x4_epistasis_scale") or hist_entry.get("s8b_epistasis_estimand_scale")
            if ep_sc_34 is not None:
                pooled_stage0_payload["obs_beta_x3_x4_epistasis_scale"] = ep_sc_34
            else:
                pooled_stage0_payload.pop("obs_beta_x3_x4_epistasis_scale", None)
        single_hist_s0 = evaluate_stage0_decision_tree(
            dict(hist_meas),
            defaulted_se_keys=sorted(hist_def_se),
            prior_batch_sticky_s9_safety_stop=False,
        )
        single_hist_mat_status = str(single_hist_s0.get("resolved_architecture_assignments", {}).get("lot_maturation_status", ""))
        if (
            single_hist_s0.get("any_fallback_switch_triggered") is True
            or single_hist_s0.get("primary_biological_or_pouch_gate_failed") is True
            or single_hist_mat_status.startswith("NO_GO_")
            or (
                single_hist_mat_status.startswith("HOLD_")
                and single_hist_mat_status not in ("HOLD_INCOMPLETE_SAFETY_PANEL", "HOLD_UNMEASURED_MANDATORY_SAFETY_PANEL")
            )
            or hist_entry.get("combined_stress_coupon_cleared_lysis_le_8pct_and_si_ge_2p5") is False
            or hist_entry.get("obs_vascularized_explant_en_bloc_recovery_verified") is False
            or hist_entry.get("obs_ap1903_pulsed_pk_washout_regimen_verified") is False
        ):
            prior_failed_stage0_gate_keys.update(hist_stage0_noncore_keys)
        if single_hist_s0.get("safety_kill_s9_triggered") is True:
            post_stop_passing_s9_batch_ids.clear()
            for cat in (single_hist_s0.get("s9_causal_failed_categories") or []):
                stop_required_s9_categories.add(str(cat))
            if (
                single_hist_s0.get("resolved_architecture_assignments", {}).get("mcb_clonogenic_escape_qualification_status")
                == "HARD_SAFETY_STOP_MCB_CLONOGENIC_SURVIVOR_DETECTED"
                or "MCB_CLONOGENIC" in stop_required_s9_categories
            ):
                stop_required_mcb_clonogenic_clearance = True
            if (
                "ER_CHOP_HC10" in stop_required_s9_categories
                or any(
                    "obs_chop_coincidence_fold" in r
                    or "obs_gate0e_chop_log2fc" in r
                    or "obs_hc10_surface_fold_wt" in r
                    or "obs_gate0e_hc10_surface_mfi_ratio" in r
                    for r in (single_hist_s0.get("triggered_rule_ids") or [])
                )
            ):
                stop_required_5endpoint_er_clearance = True
        elif (
            prior_sticky_s9
            and single_hist_s0.get("mandatory_s9_safety_panel_complete") is True
            and single_hist_s0.get("safety_kill_s9_triggered") is False
            and single_hist_s0.get("safety_kill_s9_status") == "CLEARED"
            and (
                "CORE_S9_VAF_OR_CASTSEQ" not in stop_required_s9_categories
                or (
                    float(hist_meas.get("obs_tp53_ppm1d_bcor_duplex_depth") or hist_meas.get("obs_tp53_duplex_depth") or 0.0) >= 62958.0
                    and float(hist_meas.get("obs_castseq_duplex_or_amplicon_depth") or 0.0) >= 62958.0
                )
            )
            and (
                not stop_required_mcb_clonogenic_clearance
                or single_hist_s0.get("resolved_architecture_assignments", {}).get("mcb_clonogenic_escape_qualification_status")
                == "CLEARED_MCB_CLONOGENIC_EFFECTIVE_TRIALS_GE_299573228"
            )
            and (
                not stop_required_5endpoint_er_clearance
                or single_hist_s0.get("s9_5endpoint_er_safety_panel_complete") is True
                or single_hist_s0.get("full_5endpoint_s9_panel_complete") is True
            )
            and (
                "RESIDUAL_PSC" not in stop_required_s9_categories
                or (
                    hist_meas.get("obs_residual_psc_effective_fraction") not in ("UNAVAILABLE", None)
                    and float(hist_meas.get("obs_residual_psc_effective_trials") or 0.0) >= 299573228.0
                )
            )
            and (
                "CNV_20Q11" not in stop_required_s9_categories
                or hist_meas.get("obs_20q11_bcl2l1_cnv_ratio") not in ("UNAVAILABLE", None)
            )
            and (
                "SINGLE_CELL_MOSAIC_20Q11" not in stop_required_s9_categories
                or hist_meas.get("obs_20q11_bcl2l1_single_cell_mosaic_fraction") not in ("UNAVAILABLE", None)
            )
            and (
                "POST_STRESS_20Q11_OUTGROWTH" not in stop_required_s9_categories
                or hist_meas.get("obs_post_stress_20q11_bcl2l1_outgrowth_fold") not in ("UNAVAILABLE", None)
            )
            and (
                "FUNCTIONAL_P53_NUTLIN3A" not in stop_required_s9_categories
                or hist_meas.get("obs_nutlin3a_p21_mdm2_fold_induction") not in ("UNAVAILABLE", None)
            )
            and (
                "HAZARD_PROLIF_SURVIVOR" not in stop_required_s9_categories
                or (
                    hist_meas.get("obs_hazard_state_prolif_survivor_count_k") not in ("UNAVAILABLE", None)
                    and float(hist_meas.get("obs_hazard_state_prolif_effective_trials") or 0.0) >= 299573228.0
                )
            )
            and (
                "PULSED_PK_KILL" not in stop_required_s9_categories
                or (
                    hist_meas.get("obs_ap1903_pulsed_pk_log10_kill") not in ("UNAVAILABLE", None)
                    and hist_entry.get("obs_ap1903_pulsed_pk_washout_regimen_verified") is not False
                )
            )
            and (
                "MONOTHERAPY_PULSED_PK_KILL" not in stop_required_s9_categories
                or (
                    hist_meas.get("obs_ap1903_monotherapy_pulsed_pk_log10_kill") not in ("UNAVAILABLE", None)
                    and hist_entry.get("obs_ap1903_pulsed_pk_washout_regimen_verified") is not False
                )
            )
            and (
                "STEP1_SUBLETHAL_SURVIVOR" not in stop_required_s9_categories
                or (
                    hist_meas.get("obs_ap1903_step1_sublethal_survivor_count_k") not in ("UNAVAILABLE", None)
                    and float(hist_meas.get("obs_ap1903_step1_effective_trials") or 0.0) >= 29957322.8
                )
            )
            and (
                "EXPLANT_RECOVERY" not in stop_required_s9_categories
                or (
                    hist_entry.get("obs_vascularized_explant_en_bloc_recovery_verified") is True
                    and float(hist_meas.get("obs_explant_en_bloc_recovery_fraction") or 0.0) >= 0.99
                )
            )
            and isinstance(hist_entry.get("batch_id"), str)
            and hist_entry["batch_id"].strip()
        ):
            post_stop_passing_s9_batch_ids.append(hist_entry["batch_id"].strip())
        step_hist_payload = dict(pooled_stage0_payload)
        if isinstance(hist_entry.get("batch_id"), str) and hist_entry["batch_id"].strip():
            step_hist_payload["batch_id"] = hist_entry["batch_id"].strip()
        step_hist_payload["verified_ledger_passing_s9_batch_ids"] = list(post_stop_passing_s9_batch_ids)
        for qp_k in (
            "prior_batch_sticky_s9_safety_stop",
            "prior_variant_safety_status",
            "reset_prior_safety_stop_after_investigation",
            "reset_prior_safety_stop_reason",
            "reset_prior_safety_stop_reason_code",
            "reset_prior_safety_stop_repeat_batch_ids",
            "reset_prior_safety_stop_independent_qp_approver_id",
        ):
            if hist_entry.get(qp_k) is not None:
                step_hist_payload[qp_k] = hist_entry[qp_k]
        hist_s0_eval = evaluate_stage0_decision_tree(
            step_hist_payload,
            defaulted_se_keys=sorted(pooled_defaulted_se_set),
            prior_batch_sticky_s9_safety_stop=prior_sticky_s9,
        )
        prior_sticky_s9 = bool(hist_s0_eval.get("safety_kill_s9_triggered") is True) or bool(single_hist_s0.get("safety_kill_s9_triggered") is True)
        if not prior_sticky_s9:
            stop_required_mcb_clonogenic_clearance = False
            stop_required_5endpoint_er_clearance = False
            stop_required_s9_categories.clear()
        if hist_s0_eval.get("resolved_architecture_assignments", {}).get("prior_s9_safety_stop_history_retained") is True:
            prior_s9_history_retained_across_ledger = True
    if (
        isinstance((prior_posterior_report or {}).get("stage0_decision_tree_evaluation"), dict)
        and prior_posterior_report["stage0_decision_tree_evaluation"].get("safety_kill_s9_triggered") is True
    ):
        prior_sticky_s9 = True
    if isinstance(input_payload.get("mcb_clone_id"), str) and input_payload["mcb_clone_id"].strip():
        seen_mcb_clone_ids.add(input_payload["mcb_clone_id"].strip())
    elif prev_obs_ledger:
        any_batch_omitted_mcb_clone_id = True
    curr_meas_for_oos = dict(current_obs_entry.get("measurements") or {})
    for bf_k in stage0_persisted_boolean_flags:
        if bf_k in input_payload and input_payload[bf_k] is not None:
            curr_meas_for_oos[bf_k] = input_payload[bf_k]
    curr_stage0_noncore_keys = {
        k for k, v in curr_meas_for_oos.items()
        if k not in s9_3_endpoint_keys and not k.endswith("_se") and v not in ("UNAVAILABLE", None)
    }
    if curr_stage0_noncore_keys & prior_failed_stage0_gate_keys:
        if not (
            input_payload.get("oos_retest_investigation_approved") is True
            and isinstance(input_payload.get("oos_retest_assignable_cause_code"), str)
            and input_payload["oos_retest_assignable_cause_code"].strip()
        ):
            uninvestigated_non_s9_oos_overwrite_hold = True
    multi_batch_clone_unverified_or_mixed = bool(
        len(seen_mcb_clone_ids) > 1
        or (len(prev_obs_ledger) >= 1 and (any_batch_omitted_mcb_clone_id or len(seen_mcb_clone_ids) == 0))
    )
    cross_clone_pooling_hold = bool(
        multi_batch_clone_unverified_or_mixed
        and input_payload.get("allow_cross_clone_pooling_for_exploratory_research") is not True
    )
    curr_def_se = set(measurement_validation.get("defaulted_standard_error_keys") or [])
    _drop_unpaired_rare_event_members(pooled_stage0_payload, input_payload)
    saved_x14_scale = pooled_stage0_payload.get("obs_beta_x1_x4_epistasis_scale")
    saved_x34_scale = pooled_stage0_payload.get("obs_beta_x3_x4_epistasis_scale")
    pooled_stage0_payload.pop("s8b_epistasis_estimand_scale", None)
    for k, v in input_payload.items():
        if k in ("obs_beta_x1_x4_epistasis_scale", "obs_beta_x3_x4_epistasis_scale"):
            continue
        if v not in ("UNAVAILABLE", None):
            pooled_stage0_payload[k] = v
            if k in se_to_mean_lookup:
                if k in curr_def_se or se_to_mean_lookup[k] in curr_def_se:
                    pooled_defaulted_se_set.add(k)
                else:
                    pooled_defaulted_se_set.discard(k)
                    pooled_defaulted_se_set.discard(se_to_mean_lookup[k])
    if "obs_beta_x1_x4_epistasis" in input_payload and input_payload.get("obs_beta_x1_x4_epistasis") not in ("UNAVAILABLE", None):
        curr_sc_14 = input_payload.get("obs_beta_x1_x4_epistasis_scale") or input_payload.get("s8b_epistasis_estimand_scale")
        if curr_sc_14 is not None:
            pooled_stage0_payload["obs_beta_x1_x4_epistasis_scale"] = curr_sc_14
        else:
            pooled_stage0_payload.pop("obs_beta_x1_x4_epistasis_scale", None)
    elif saved_x14_scale is not None:
        pooled_stage0_payload["obs_beta_x1_x4_epistasis_scale"] = saved_x14_scale
    if "obs_beta_x3_x4_epistasis" in input_payload and input_payload.get("obs_beta_x3_x4_epistasis") not in ("UNAVAILABLE", None):
        curr_sc_34 = input_payload.get("obs_beta_x3_x4_epistasis_scale") or input_payload.get("s8b_epistasis_estimand_scale")
        if curr_sc_34 is not None:
            pooled_stage0_payload["obs_beta_x3_x4_epistasis_scale"] = curr_sc_34
        else:
            pooled_stage0_payload.pop("obs_beta_x3_x4_epistasis_scale", None)
    elif saved_x34_scale is not None:
        pooled_stage0_payload["obs_beta_x3_x4_epistasis_scale"] = saved_x34_scale
    pooled_stage0_payload["batch_id"] = batch_id
    pooled_stage0_payload["verified_ledger_passing_s9_batch_ids"] = list(post_stop_passing_s9_batch_ids)
    stage0_eval_report = evaluate_stage0_decision_tree(
        pooled_stage0_payload,
        require_full_5endpoint_s9_panel=require_full_5endpoint_s9_panel,
        defaulted_se_keys=sorted(pooled_defaulted_se_set),
        prior_batch_sticky_s9_safety_stop=prior_sticky_s9,
    )
    if (
        prior_s9_history_retained_across_ledger
        and stage0_eval_report.get("safety_kill_s9_triggered") is False
        and isinstance(stage0_eval_report.get("resolved_architecture_assignments"), dict)
    ):
        stage0_eval_report["resolved_architecture_assignments"]["prior_s9_safety_stop_history_retained"] = True
        if stage0_eval_report["resolved_architecture_assignments"].get("lot_maturation_status") == "READY_FOR_ENCAPSULATION":
            stage0_eval_report["resolved_architecture_assignments"]["lot_maturation_status"] = (
                "READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY"
            )
        if stage0_eval_report.get("lot_release_readiness_tier") == "READY_FOR_ENCAPSULATION_COMPLETE_STAGE0_EVIDENCE":
            stage0_eval_report["lot_release_readiness_tier"] = "READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = (
                "READY_WITH_RETAINED_PRIOR_S9_INVESTIGATION_HISTORY"
            )
        elif stage0_eval_report.get("lot_release_readiness_tier") == "PROVISIONAL_CORE_14_PANEL_CLEARED_EXTENDED_NK_PHAGO_DEPTH_CHRONIC_UNMEASURED":
            stage0_eval_report["lot_release_readiness_tier"] = "PROVISIONAL_CORE_14_PANEL_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = (
                "PROVISIONAL_CORE_14_PANEL_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY"
            )
        elif stage0_eval_report.get("lot_release_readiness_tier") == "PROVISIONAL_S9_AND_MATURITY_CLEARED_IMMUNE_AND_OCR_GATES_UNMEASURED":
            stage0_eval_report["lot_release_readiness_tier"] = "PROVISIONAL_S9_AND_MATURITY_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = (
                "PROVISIONAL_S9_AND_MATURITY_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY"
            )
        elif stage0_eval_report.get("lot_release_readiness_tier") == "PROVISIONAL_CORE_S9_CLEARED_UNMEASURED_GATES_PENDING":
            stage0_eval_report["lot_release_readiness_tier"] = "PROVISIONAL_CORE_S9_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = (
                "PROVISIONAL_CORE_S9_CLEARED_WITH_RETAINED_PRIOR_S9_HISTORY"
            )
    if cross_clone_pooling_hold and isinstance(stage0_eval_report.get("resolved_architecture_assignments"), dict):
        stage0_eval_report["resolved_architecture_assignments"]["cross_clone_pooled_evidence_status"] = (
            "HOLD_CROSS_CLONE_POOLED_EVIDENCE_REQUIRES_SINGLE_CLONE_LOT_QUALIFICATION"
        )
        if not stage0_eval_report.get("safety_kill_s9_triggered"):
            stage0_eval_report["resolved_architecture_assignments"]["lot_maturation_status"] = (
                "HOLD_CROSS_CLONE_POOLED_EVIDENCE_REQUIRES_SINGLE_CLONE_LOT_QUALIFICATION"
            )
            stage0_eval_report["lot_release_readiness_tier"] = "HOLD_FAILED_OR_INCOMPLETE_GATE"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = "HOLD_FAILED_OR_INCOMPLETE_GATE"
    elif (
        multi_batch_clone_unverified_or_mixed
        and input_payload.get("allow_cross_clone_pooling_for_exploratory_research") is True
        and isinstance(stage0_eval_report.get("resolved_architecture_assignments"), dict)
    ):
        stage0_eval_report["resolved_architecture_assignments"]["cross_clone_pooled_evidence_status"] = (
            "EXPLORATORY_CROSS_CLONE_POOLED_EVIDENCE_NON_RELEASE_ONLY"
        )
        if not stage0_eval_report.get("safety_kill_s9_triggered"):
            stage0_eval_report["resolved_architecture_assignments"]["lot_maturation_status"] = (
                "EXPLORATORY_CROSS_CLONE_POOLED_NON_RELEASE"
            )
            stage0_eval_report["lot_release_readiness_tier"] = "EXPLORATORY_CROSS_CLONE_POOLED_NON_RELEASE"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = (
                "EXPLORATORY_CROSS_CLONE_POOLED_NON_RELEASE"
            )
    if uninvestigated_non_s9_oos_overwrite_hold and isinstance(stage0_eval_report.get("resolved_architecture_assignments"), dict):
        stage0_eval_report["resolved_architecture_assignments"]["non_s9_oos_retest_governance_status"] = (
            "HOLD_UNINVESTIGATED_PRIOR_BATCH_OOS_OVERWRITE"
        )
        if not stage0_eval_report.get("safety_kill_s9_triggered"):
            stage0_eval_report["resolved_architecture_assignments"]["lot_maturation_status"] = (
                "HOLD_UNINVESTIGATED_PRIOR_BATCH_OOS_OVERWRITE"
            )
            stage0_eval_report["lot_release_readiness_tier"] = "HOLD_FAILED_OR_INCOMPLETE_GATE"
            stage0_eval_report["resolved_architecture_assignments"]["lot_release_readiness_tier"] = "HOLD_FAILED_OR_INCOMPLETE_GATE"
    result_report = {
        "schema_version": "0.8.0",
        "calibration_state_version": CALIBRATION_STATE_VERSION,
        "ingestor": "OpenIslet Closed-Loop Wet-Lab Conjugate Bayesian Updater",
        "calibration_mode": (
            "SYNTHETIC_DRY_RUN_VERIFICATION_DO_NOT_CLAIM_AS_WETLAB_DATA"
            if is_dry_run
            else "EMPIRICAL_WETLAB_INGESTION"
        ),
        "epistemic_lineage_status": (
            "SYNTHETIC_DRY_RUN_LINEAGE"
            if is_dry_run
            else "PURE_EMPIRICAL_WETLAB_LINEAGE"
        ),
        "batch_id": batch_id,
        "ingested_batch_ids": prev_batches + [batch_id],
        "batch_payload_sha256": batch_payload_sha256,
        "ingested_batch_payload_sha256s": prev_payload_hashes + [batch_payload_sha256],
        "ordered_batch_observation_ledger": prev_obs_ledger + [current_obs_entry],
        "parent_report_sha256": parent_report_sha256,
        "observed_endpoints_updated_count": observed_count,
        "minimum_biological_se_floor_enforced": MIN_BIOLOGICAL_SE_FLOOR,
        "minimum_prior_sd_positivity_floor_enforced": MIN_PRIOR_SD_POSITIVITY_FLOOR,
        "prior_data_conflict_z_threshold": PRIOR_DATA_CONFLICT_Z_THRESHOLD,
        "any_prior_data_conflict_detected_gt_3sigma": any(
            bool(u.get("prior_data_conflict_detected", False))
            for u in list(updates.values()) + list(optional_pouch_and_screen_updates.values())
        ),
        "supported_8_canonical_endpoints_v80": [
            "obs_stage7_yield_ratio",
            "obs_cd47_fold",
            "obs_serpinb9_log10_kass",
            "obs_icasp9_log10_kill",
            "obs_pouch_min_po2_mmhg",
            "obs_pouch_max_h2o2_uM",
            "obs_pouch_max_ph",
            "obs_ifng_hla_dr_dq_negative_pct",
        ],
        "input_payload_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "input_payload_received_sha256": received_payload_sha256,
        "measurement_validation": measurement_validation,
        "conjugate_posterior_updates_sha256": hashlib.sha256(updates_bytes).hexdigest(),
        "layer2_and_ranker_reconciliation_note_claude_v80_013_017": {
            "status": "TARGET_PRIORS_SEPARATE_FROM_MODEL_REFERENCES_AND_OBSERVATIONS",
            "current_model_references": model_reference_metadata,
            "exact_product_assays": "PREREGISTERED_PENDING",
            "synthetic_calibration_is_empirical": False,
        },
        "conjugate_posterior_updates": updates,
        "optional_pouch_and_36run_screen_updates_v80_070": optional_pouch_and_screen_updates,
        "stage0_decision_tree_evaluation": stage0_eval_report,
    }
    result_report["calibration_state_sha256"] = _calibration_state_sha256(result_report)
    result_report["state_integrity_scope"] = "Checks complete-state consistency and replay history; hashes and PURE_EMPIRICAL labels do not authenticate bench records or donor/clone independence."
    # Strictly verify output has zero NaN/Inf values before returning (AUDIT-REF):
    json.dumps(result_report, allow_nan=False)
    return result_report


def _validate_safe_cli_path(raw_path: str, label: str) -> Path:
    """Validate CLI file path against symlink traversal outside REPO_ROOT or system temp directory (AUDIT-REF)."""
    p = Path(raw_path).expanduser().absolute()
    if p.is_symlink():
        raise ValueError(f"{label} path {p} must not be a symbolic link (AUDIT-REF)")
    macos_system_root_symlinks = {Path("/var"), Path("/tmp"), Path("/etc")}
    for ancestor in p.parents:
        if ancestor in macos_system_root_symlinks:
            continue
        if ancestor.is_symlink():
            raise ValueError(f"{label} ancestor {ancestor} must not be a symbolic link (AUDIT-REF)")
    resolved = p.resolve()
    allowed_roots = (
        REPO_ROOT.resolve(),
        Path(tempfile.gettempdir()).resolve(),
        Path("/tmp").resolve(),
    )
    if not any(resolved == root or root in resolved.parents for root in allowed_roots):
        raise ValueError(f"{label} path {resolved} must reside inside repository or temporary directory (AUDIT-REF)")
    if p.exists():
        info = p.lstat()
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError(f"{label} must be a regular file with one hardlink")
    return resolved


@contextmanager
def _safe_parent_directory(path: Path, *, create: bool = False):
    """Pin every directory component without following attacker-swapped links."""
    resolved = _validate_safe_cli_path(str(path), "JSON path")
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    current = os.open("/", flags)
    try:
        for part in resolved.parent.parts[1:]:
            try:
                nxt = os.open(part, flags, dir_fd=current)
            except FileNotFoundError:
                if not create:
                    raise
                try:
                    os.mkdir(part, dir_fd=current)
                except FileExistsError:
                    pass
                nxt = os.open(part, flags, dir_fd=current)
            os.close(current)
            current = nxt
        yield current, resolved.name
    finally:
        os.close(current)


def _regular_single_link_at(parent: int, name: str, *, missing_ok: bool = False) -> None:
    try:
        info = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        if missing_ok:
            return
        raise
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError("JSON path must be a regular file with one hardlink")


def _read_safe_json_file(path: Path) -> tuple[Any, bytes]:
    """Open without following a swapped link or blocking on a special file."""
    with _safe_parent_directory(path) as (parent, name):
        flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        fd = os.open(name, flags, dir_fd=parent)
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise ValueError("JSON input must be a regular file with one hardlink")
            with os.fdopen(fd, "rb") as stream:
                fd = -1
                raw = stream.read()
            return loads_strict_json(raw.decode("utf-8")), raw
        finally:
            if fd >= 0:
                os.close(fd)


def _existing_json_bytes(path: Path) -> bytes | None:
    _validate_safe_cli_path(str(path), "JSON snapshot")
    try:
        return _read_safe_json_file(path)[1]
    except FileNotFoundError:
        return None


def _atomic_write_bytes(target_path: Path, content: bytes) -> None:
    with _safe_parent_directory(target_path, create=True) as (parent, name):
        _regular_single_link_at(parent, name, missing_ok=True)
        temporary = f".{name}.{uuid.uuid4().hex}.tmp"
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
        fd = os.open(temporary, flags, 0o600, dir_fd=parent)
        try:
            with os.fdopen(fd, "wb") as stream:
                fd = -1
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            _regular_single_link_at(parent, name, missing_ok=True)
            os.replace(temporary, name, src_dir_fd=parent, dst_dir_fd=parent)
            os.fsync(parent)
        finally:
            if fd >= 0:
                os.close(fd)
            try:
                os.unlink(temporary, dir_fd=parent)
            except FileNotFoundError:
                pass


def _atomic_write_json(target_path: Path, data: Dict[str, Any]) -> None:
    _atomic_write_bytes(target_path, (json.dumps(data, indent=2, allow_nan=False) + "\n").encode("utf-8"))


@contextmanager
def _calibration_output_lock(path: Path):
    # Keep the small lock inode: unlinking it after unlocking permits two later
    # writers to lock different inodes for the same output. Dry-run/no-write
    # execution never enters this context and never creates it.
    lock_path = path.with_name(f".{path.name}.lock")
    with _safe_parent_directory(lock_path, create=True) as (parent, name):
        flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        fd = os.open(name, flags, 0o600, dir_fd=parent)
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise ValueError("Calibration lock must be a regular file with one hardlink")
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError("Another calibration writer holds this output lock") from exc

            def verify_ownership():
                # An operator/revoker can unlink or replace a flock inode. The
                # old descriptor remains locked but no longer excludes writers
                # opening the replacement pathname.
                try:
                    with _safe_parent_directory(lock_path) as (current_parent, current_name):
                        named = os.stat(current_name, dir_fd=current_parent, follow_symlinks=False)
                    held = os.fstat(fd)
                except (OSError, ValueError) as exc:
                    raise ValueError("Calibration writer lock ownership was revoked") from exc
                if (not stat.S_ISREG(named.st_mode) or named.st_nlink != 1 or held.st_nlink != 1
                        or (named.st_dev, named.st_ino) != (held.st_dev, held.st_ino)):
                    raise ValueError("Calibration writer lock ownership was revoked")

            verify_ownership()
            yield verify_ownership
            verify_ownership()
        finally:
            os.close(fd)


def _publish_calibration_report(
    output_path: Path, report: Dict[str, Any], expected_snapshots: Dict[Path, bytes | None],
    *, force_overwrite_empirical: bool = False,
) -> None:
    """Reject stale writers, protect empirical lineage, and roll back alias failures.

    Each rename is atomic; the canonical/legacy pair is rollback protected under
    one writer lock. This is not a multi-file transaction for uncooperative readers.
    """
    output_path = _validate_safe_cli_path(str(output_path), "output")
    canonical_path, alias_path = DEFAULT_OUT.resolve(), LEGACY_ALIAS_PATH.resolve()
    def portable_identity(path: Path) -> Path:
        # Same NFC/casefold contract as study_artifact_exports._portable_path_key;
        # it also protects reserved targets before either file exists.
        return Path(unicodedata.normalize("NFC", str(path)).casefold())
    def same_destination(left: Path, right: Path) -> bool:
        if portable_identity(left) == portable_identity(right):
            return True
        try:
            return left.samefile(right)
        except FileNotFoundError:
            return False
    if same_destination(canonical_path, alias_path):
        raise ValueError("Canonical calibration and legacy alias require distinct portable identities")
    if same_destination(output_path, alias_path):
        raise ValueError("The legacy alias cannot be an independent output; publish through the canonical calibration path")
    if same_destination(output_path, canonical_path) and output_path != canonical_path:
        raise ValueError("Use the canonical calibration path spelling so every writer shares its publication lock")
    targets = [output_path]
    if output_path == canonical_path:
        targets.append(_validate_safe_cli_path(str(LEGACY_ALIAS_PATH), "legacy alias"))
    if set(expected_snapshots) != set(targets):
        raise ValueError("Output publication requires an entry snapshot for every target")
    published_bytes = (json.dumps(report, indent=2, allow_nan=False) + "\n").encode("utf-8")
    with _calibration_output_lock(output_path) as verify_ownership:
        for target in targets:
            current = _existing_json_bytes(target)
            if current != expected_snapshots[target]:
                raise ValueError("Calibration output changed since this update began; reload the latest posterior")
            if current is None:
                continue
            existing = loads_strict_json(current.decode("utf-8"))
            if isinstance(existing, dict) and existing.get("calibration_mode") == "EMPIRICAL_WETLAB_INGESTION":
                if report["calibration_mode"] != "EMPIRICAL_WETLAB_INGESTION":
                    if not force_overwrite_empirical:
                        raise ValueError("Refusing to overwrite empirical calibration with a synthetic dry-run")
                elif report.get("parent_report_sha256") != hashlib.sha256(
                        json.dumps(existing, sort_keys=True, allow_nan=False).encode("utf-8")).hexdigest():
                    raise ValueError("Empirical overwrite must extend the currently published posterior")
        published: list[Path] = []
        try:
            for target in targets:
                verify_ownership()
                if _existing_json_bytes(target) != expected_snapshots[target]:
                    raise ValueError("Calibration output changed during publication; preserving intervening bytes")
                # Track attempts before writing: rename may have succeeded even
                # when a subsequent directory fsync raises an exception.
                published.append(target)
                _atomic_write_json(target, report)
                verify_ownership()
        except BaseException as original:
            rollback_errors = []
            for target in reversed(published):
                try:
                    verify_ownership()
                    old = expected_snapshots[target]
                    current = _existing_json_bytes(target)
                    if current == old:
                        continue  # The attempted write failed before replacement.
                    if current != published_bytes:
                        raise ValueError(f"Preserving intervening calibration bytes at {target}")
                    if old is None:
                        with _safe_parent_directory(target) as (parent, name):
                            try:
                                os.unlink(name, dir_fd=parent)
                            except FileNotFoundError:
                                pass
                            os.fsync(parent)
                    else:
                        _atomic_write_bytes(target, old)
                except BaseException as failure:
                    rollback_errors.append(str(failure))
            if rollback_errors:
                raise RuntimeError("Calibration publication failed and rollback requires manual recovery: "
                                   + "; ".join(rollback_errors)) from original
            raise


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest wet-lab JSON readouts to update OpenIslet Bayesian priors.")
    parser.add_argument("--input-json", type=str, default="", help="Path to empirical wet-lab readout JSON")
    parser.add_argument("--prior-posterior-json", type=str, default="", help="Optional prior posterior JSON for sequential Bayesian updating")
    parser.add_argument(
        "--dry-run-calibration",
        "--synthetic-dry-run",
        dest="dry_run_calibration",
        action="store_true",
        help="Run synthetic dry-run pipeline verification",
    )
    parser.add_argument("--output", type=str, default=str(DEFAULT_OUT), help="Output JSON path")
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Validate and compute calibration report in-memory without writing output files to disk",
    )
    parser.add_argument(
        "--force-overwrite-empirical",
        action="store_true",
        help="Allow a synthetic dry-run to overwrite an existing EMPIRICAL_WETLAB_INGESTION artifact",
    )
    parser.add_argument(
        "--fail-on-conflict-gt-3sigma",
        action="store_true",
        help="Fail closed if any observed endpoint has z_conflict > 3.0",
    )
    args = parser.parse_args()

    if args.input_json and args.dry_run_calibration:
        raise ValueError("Cannot specify both --input-json and --dry-run-calibration/--synthetic-dry-run simultaneously.")
    if not args.input_json and not args.dry_run_calibration:
        raise ValueError("Must specify either --input-json <path> or --dry-run-calibration (--synthetic-dry-run).")

    prior_report = None
    prior_path = None
    prior_raw = None
    if args.prior_posterior_json:
        prior_path = _validate_safe_cli_path(args.prior_posterior_json, "prior_posterior_json")
        prior_report, prior_raw = _read_safe_json_file(prior_path)

    if args.input_json:
        in_path = _validate_safe_cli_path(args.input_json, "input_json")
        payload, _input_raw = _read_safe_json_file(in_path)
        is_dry = False
    else:
        payload = {
            "batch_id": "DRY-RUN-CALIBRATION-001",
            "obs_stage7_yield_ratio": 0.982,
            "obs_stage7_yield_se": 0.018,
            "obs_cd47_fold": 7.95,
            "obs_cd47_fold_se": 0.35,
            "obs_serpinb9_log10_kass": 6.21,
            "obs_serpinb9_log10_kass_se": 0.09,
            "obs_icasp9_log10_kill": 4.48,
            "obs_icasp9_log10_kill_se": 0.11,
        }
        is_dry = True

    expected_snapshots = {}
    if not args.no_write:
        out_path = _validate_safe_cli_path(args.output, "output")
        targets = [out_path] + ([LEGACY_ALIAS_PATH.resolve()] if out_path == DEFAULT_OUT.resolve() else [])
        for target in targets:
            # Preserve the exact state actually used for an in-place update,
            # including a race occurring between the input and output reads.
            expected_snapshots[target] = prior_raw if target == prior_path else _existing_json_bytes(target)

    report = run_calibration(
        payload,
        is_dry_run=is_dry,
        prior_posterior_report=prior_report,
        fail_on_conflict_gt_3sigma=args.fail_on_conflict_gt_3sigma,
    )
    if args.no_write:
        print(f"[Wet-Lab Ingestor v0.8.0] Verified in-memory (--no-write, mode={report['calibration_mode']})")
        return

    _publish_calibration_report(out_path, report, expected_snapshots,
                                force_overwrite_empirical=args.force_overwrite_empirical)
    print(f"[Wet-Lab Ingestor v0.8.0] Wrote {out_path} (mode={report['calibration_mode']})")


if __name__ == "__main__":
    main()
