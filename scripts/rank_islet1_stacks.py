#!/usr/bin/env python3
"""
Step 5 & 6: Project Islet-1 Sign-Sensitive, Endpoint-Oriented & Locus-Penalized Evaluator (v0.6.9).

Upgrades in v0.6.0 (incorporating Red-Team Reviews):
1. Zero-Evidence Baseline = 0.0:
   - If target_profiles is empty (`{}`), `run_combinatorial_multi_gate_search({})` returns 0.0
     across all gates and composite scores.
2. Broadcast Attribution Eliminated (Audit §5.2.1):
   - Atomic effects are credited ONLY to genes explicitly named in `focal_factor_X` or
     `verbatim_evidence_quote` (or to the single target gene in single-gene studies).
     Multi-gene studies never broadcast every effect to all co-mentioned genes.
3. Endpoint-Orientation-Aware Sign & Counter-Evidence Sensitivity (Audit §2.2 + Audit §4.2/§5.2):
   - Preserves algebraic signs of `log_response_ratio_or_smd` (`sign_flip` changes profiles & rankings).
   - Orients lower-is-better hazard endpoints (`specific lysis`, `apoptosis`, `cell death`) vs
     higher-is-better functional endpoints (`GSIS`, `viability`, `survival`, `C-peptide`).
   - Separates `effect_direction == "Deleterious_or_Impaired"` from independent
     `secretome_or_differentiation_penalty_observed`.
4. Exact Specification Parity & Base-Edit Penalty (Audit §1.3):
   - `search_cuts == named_cuts` and `search_base_edits_B == named_base_edits_B` across all architectures.
   - Implements both chromosomal DSB translocation penalty `exp(-0.15 * max(0, K - 4))` and
     nickase base/prime-edit penalty `0.965 ** B` (`PMID:39529170`).
   - Unifies `Chassis B` (`HLA-A/B/CIITA-/- [retaining -21M HLA-C] + TYK2-/- + CXCL10::PD-L1-P2A-SERPINB9 + AAVS1-CAG-iCasp9 + INS base edit`, `K=6, B=1`).
   - Locks `HLA-E*01:03` single-chain trimer nonamer to `VMAPRTLVL` (HLA-A*02:01 leader, minimizing NKG2C+ adaptive NK activation) and locks suicide switch to biallelic `AAVS1-CAG-iCasp9` (`<10^-8` escape bound).
5. Explicit Engineering Priors, Posterior Uncertainty & Executable Graveyard Rules:
   - Exposes `engineering_priors_and_sensitivity_table`, `pareto_frontier_by_locus_budget_K`,
     `global_rank_among_all_combinations`, `lower_10pct_cred_bound`, and `triggered_graveyard_rules`.
"""

from __future__ import annotations

import copy
import json
import hashlib
import importlib.util
from types import MappingProxyType
import math
import re
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXTRACTED_DIR = PROJECT_ROOT / "data" / "extracted"
OMICS_DIR = PROJECT_ROOT / "data" / "omics"
GRAVEYARD_PATH = PROJECT_ROOT / "data" / "graveyard" / "idea_graveyard_v0.5.json"
GRAPH_PATH = EXTRACTED_DIR / "cross_study_discovery_graph.json"
OUT_PATH = EXTRACTED_DIR / "islet1_candidate_stack_rankings.json"

def _load_local_helper(name):
    spec = importlib.util.spec_from_file_location("openislet_" + name, PROJECT_ROOT / "scripts" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


engineering_complexity_discount = _load_local_helper("engineering_priors").engineering_complexity_discount
validate_public_excerpt_token = _load_local_helper("public_excerpt_policy").validate_token


def _finite_real(value, label: str) -> float:
    """Reject unknown/non-real data before arithmetic, conversion or clipping."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite real number")
    try:
        number = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{label} is outside the finite numeric domain") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} must be a finite real number")
    return number


def _profile_has_empirical_studies(profile: dict) -> bool:
    if not isinstance(profile, dict):
        raise ValueError("Target profiles must contain explicit gene mappings")
    count = profile.get("n_empirical_studies", profile.get("n_studies_total", 0))
    if type(count) is not int or count < 0:
        raise ValueError("Empirical study counts must be nonnegative integers")
    return count > 0


def _profile_standard_error(profile: dict, default: float) -> float:
    if not isinstance(profile, dict):
        raise ValueError("Target profiles must contain explicit gene mappings")
    value = _finite_real(profile.get("posterior_se", default), "Gene evidence dispersion")
    if value < 0:
        raise ValueError("Gene evidence dispersion must be nonnegative")
    return value


def _effect_orientation(effect: dict) -> int:
    """One typed endpoint orientation for private and redacted representations."""
    polarity = effect.get("endpoint_polarity_registry")
    if polarity is None:
        # Compatibility for deliberately minimal qualitative fixtures. Numeric
        # contrasts require the explicit measured-endpoint polarity contract.
        if effect.get("log_response_ratio_or_smd") is not None:
            raise ValueError("Numeric evidence requires an endpoint_polarity_registry")
        return int(get_endpoint_orientation_sign((effect.get("outcome_variable_Y") or "") + " " + (effect.get("verbatim_evidence_quote") or "")))
    expected = {"Higher_Is_Better_Protective_Endpoint": 1, "Lower_Is_Better_Hazard_Endpoint": -1, "Qualitative_or_Unoriented_Endpoint": 0}
    if not isinstance(polarity, dict):
        raise ValueError("endpoint_polarity_registry must be a dictionary")
    kind = polarity.get("polarity_class"); sign = polarity.get("orientation_sign")
    if not isinstance(kind, str) or kind not in expected or type(sign) is not int or sign != expected[kind]:
        raise ValueError("Endpoint polarity class and exact integer orientation must agree")
    if effect.get("log_response_ratio_or_smd") is not None and sign == 0:
        raise ValueError("Numeric evidence requires an oriented endpoint")
    return sign


FIDELITY_WEIGHTS = {
    "Tier_1_Human_InVivo": 1.00,
    "Tier_2A_Human_Cells_InVivo_NHP_or_HumanizedMouse": 0.85,
    "Tier_2B_Human_Cells_InVitro_HumanImmune_CoCulture": 0.75,
    "Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture": 0.75,
    "Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse": 0.60,
    "Tier_4_Rodent_or_Xenogeneic_Only": 0.25,
}

# Explicit, auditable table of structural & microenvironmental priors with sensitivity ranges (Audit §5.4/§5.6 + Audit A03/A07)
ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE = {
    "locus_dsb_penalty_rate_lambda_K": {
        "point_estimate": 0.15,
        "plausible_range": [0.08, 0.25],
        "formula": "exp(-lambda_K * max(0, K - 4))",
        "source_rationale": "Clonal attrition, chromosomal translocation, and p53 selection penalty for multiplex Cas9/Cas12a DSBs above K=4.",
    },
    "base_edit_nickase_penalty_factor_gamma_B": {
        "point_estimate": 0.965,
        "plausible_range": [0.920, 0.985],
        "formula": "gamma_B ** B",
        "source_rationale": "Conservative engineering prior for multiplex nickase single-strand break conversion to DSBs/deletions and deaminase bystander/structural-variant risk per base-edited locus B (extrapolated from Wu et al. 2024 PMID:39529170, which attributed elevated structural variants under ABE to deaminase activity vs Cas9/PE in single-locus cell-line edits).",
    },
    "indirect_apc_presentation_residual_hazard_cap": {
        "point_estimate": -0.25,
        "plausible_range": [-0.60, -0.10],
        "formula": "Applied to Gate 2B logit for MHC-null designs lacking indirect APC stealth validation",
        "source_rationale": "Donor B2M-/- + CIITA-/- eliminates direct presentation but not host DC phagocytosis of shed secretome/debris presented on host HLA-DQ2/DQ8 (Audit Review §3.4).",
    },
    "resting_b2m_promoter_hour0_density_discount": {
        "point_estimate": -0.20,
        "plausible_range": [-0.55, 0.00],
        "formula": "Applied to Gate 3 logit for endogenous B2M-promoter bicistronic HLA-E-P2A-CD47 pending QuantiBRITE Hour-0 verification",
        "source_rationale": "Resting Stage 7 beta-cells have low basal B2M transcription prior to IFN-gamma exposure (Audit Review §3.1).",
    },
    "stage7_trihormonal_120um_mvf_cao2_niche_prior": {
        "point_estimate": 0.45,
        "plausible_range": [0.20, 0.70],
        "formula": "Applied to Gate 1 and Gate 5 when cluster diameter <=120 um + Stage 7 maturation + microvessels are used",
        "source_rationale": "Balboa 2022 (PMID:35241836), Aghazadeh 2021 (PMID:34480863), and Hill/Krogh R_crit = 74.6 um radial diffusion limit.",
    },
    "tyk2_stage0_differentiation_yield_penalty": {
        "point_estimate": -0.32,
        "plausible_range": [-0.55, -0.15],
        "formula": "Applied to Gate 1 logit when clonal TYK2-/- is carried through Stage 0->7 hPSC endocrine differentiation",
        "source_rationale": "Chandra et al. 2022 (PMID:36289205) demonstrated ~40-60% compromised endocrine progenitor emergence in TYK2-/- hPSC differentiation.",
    },
}

# Complete registry of all 46 prior gate-logit entries (corresponding to 48 float literals in compute_series_system_utility because 2 composite entries decompose into +positive - discount pairs) across Gates 1..5, plus all 8 live-wired auxiliary_evaluator_hyperparameters (Audit A03, Grok P2-3/P2-8, Audit R10)
COMPLETE_GATE_COEFFICIENT_REGISTRY = {
    "pw_base_scale": {
        "logit_value": 0.08,
        "plausible_range": [0.05, 0.12],
        "expression": "Global multiplicative base scale factor pw applied to all 46 structural gate-logit multipliers",
        "rationale": "Global structural logit scaling constant (default 0.08) mapping dimensionless module multipliers to log-odds increments (Finding 083)",
    },
    "cross_gate_niche_coupling_note": (
        "NICHE_TRIHORMONAL_120UM_MVF_CAO2 contributes to Gate 1 (+0.45), Gate 4 (+0.25 when er_k == ER_NONE), "
        "and Gate 5 (+0.45) because <=120 um Stage 7 tri-hormonal vascularized clusters physically couple "
        "paracrine GSIS maturation, lower basal ER/hypoxic stress, and 72h oxygen engraftment across the series system."
    ),
    "Gate1_Maturation_GSIS_and_Paracrine": {
        "base_intercept": {"logit_value": 0.40, "expression": "0.08 * 5.0", "rationale": "Baseline hPSC-islet Stage 6/7 endocrine maturation prior"},
        "NICHE_TRIHORMONAL_120UM_MVF_CAO2": {"logit_value": 0.45, "expression": "0.08 * 5.625 + _pdelta('stage7_trihormonal_120um_mvf_cao2_niche_prior', 0.45)", "rationale": "Stage 7 ZM+NAC+T3 tri-hormonal maturation + <=120 um diameter (PMID:35241836)"},
        "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION": {"logit_value": -0.45, "expression": "-0.08 * 5.625", "rationale": "Diffusional lag and central hypoxia in cell-impermeable macrodevice"},
        "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO": {"logit_value": -0.32, "expression": "-0.08 * 4.0 + _pdelta('tyk2_stage0_differentiation_yield_penalty', -0.32)", "rationale": "Clonal TYK2-/- compromises Stage 0->7 endocrine progenitor emergence by ~40-60% (PMID:36289205)"},
        "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR": {"logit_value": -0.55, "expression": "-0.08 * 6.875", "rationale": "G6PC2-/- fasting glucose set-point downshift + high-pLI PTPN2/GAD2 disruption (PMID:32569842, PMID:32432213)"},
        "ER_RNLS_KO_or_PRAMLINTIDE": {"logit_value": -0.38, "expression": "-0.08 * 4.75", "rationale": "Human iPSC-beta RNLS-/- elevated basal OCR / metabolic precaution prior (poor dynamic GSIS persisted in both WT and RNLS KO in PMID:39988068 Fig. 6)"},
    },
    "Gate2A_Allo_Adaptive_T_and_DSA_Stealth": {
        "base_intercept": {"logit_value": 0.10, "expression": "0.08 * 1.25", "rationale": "Baseline allogeneic T-cell/DSA survival prior"},
        "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO": {"logit_value": 0.55, "expression": "0.08 * 6.875", "rationale": "Complete classical HLA-A/B/C and HLA-II deletion + single-chain HLA-E/CD47 (PMID:37156915, PMID:28504668)"},
        "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO": {"logit_value": 0.40, "expression": "0.08 * 5.0", "rationale": "HLA-A/B/CIITA KO retaining 1 recipient-matched HLA-C allele + TYK2-/- (PMID:30853558)"},
        "ALLO_B2M_CIITA_KO_CD47_LENTI": {"logit_value": 0.38, "expression": "0.08 * 4.75", "rationale": "B2M-/- + CIITA-/- + lentiviral CD47 OE (UP421 / HIP benchmark, PMID:37156915, PMID:40757665)"},
        "ALLO_B2M_KO_PDL1_A20_KI": {"logit_value": -0.10, "expression": "0.08 * -1.25", "rationale": "B2M-/- + PD-L1/A20 KI leaving CIITA (HLA-II) intact (VCTX210 benchmark)"},
        "ALLO_NONE_MACRODEVICE": {"logit_value": -0.60, "expression": "0.08 * -7.5", "rationale": "Unedited allogeneic HLA inside cell-impermeable macrodevice (shed alloantigen DSA priming)"},
        "ALLO_NONE_NAKED": {"logit_value": -1.80, "expression": "0.08 * -22.5", "rationale": "Unedited allogeneic HLA without encapsulation or immunosuppression"},
    },
    "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth": {
        "base_intercept_pediatric": {"logit_value": 0.20, "expression": "0.08 * 2.5", "rationale": "Pediatric early-onset T1D high autoreactive memory baseline"},
        "base_intercept_lada": {"logit_value": 0.30, "expression": "0.08 * 3.75", "rationale": "Adult-onset LADA lower autoreactive memory baseline"},
        "B2M_NULL_DIRECT_CLOAK_PLUS_INDIRECT_APC_DISCOUNT": {"logit_value": 0.30, "expression": "0.08 * 6.875 - 0.08 * 3.125 + _pdelta('indirect_apc_presentation_residual_hazard_cap', -0.25)", "rationale": "Eliminates direct CD8+ autoantigen presentation (+0.55) minus host HLA-DQ2/DQ8 DC indirect presentation hazard (-0.25)"},
        "HLA_C_RETAINED_AUTO_NONE": {"logit_value": -0.55, "expression": "-0.08 * 6.875", "rationale": "Retained HLA-C presents unedited insulin/IGRP/ZnT8 epitopes to T1D memory CD8+ T-cells"},
        "HLA_C_RETAINED_AUTO_INS_EDIT": {"logit_value": 0.25, "expression": "0.08 * 3.125", "rationale": "INS leader/DRiP base edit reduces preproinsulin/DRiP presentation on retained HLA-C (PMID:42690930)"},
        "HLA_C_RETAINED_AUTO_FULL_KO": {"logit_value": 0.20, "expression": "0.08 * 2.5", "rationale": "Multi-autoantigen KO/base-edit on retained HLA-C"},
        "ALLO_NONE": {"logit_value": -1.10, "expression": "-0.08 * 13.75", "rationale": "Unedited HLA-I/II fully exposed to T1D autoreactive memory"},
        "CYTO_CXCL10_SWAP_PDL1_SERPINB9": {"logit_value": 0.35, "expression": "0.08 * 4.375", "rationale": "CXCL10 chemokine deletion + inducible PD-L1/SERPINB9 local checkpoint/granzyme shield (PMID:36055241, PMID:19641142, PMID:30269996)"},
        "CYTO_CONSTITUTIVE_IL10_VIRAL_OE": {"logit_value": -0.75, "expression": "-0.08 * 9.375", "rationale": "Constitutive viral IL-10 pleiotropic mast-cell/granulocyte hazard (Graveyard 008)"},
        "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1": {"logit_value": -0.85, "expression": "-0.08 * 10.625", "rationale": "Epistatic deadlock: STAT1-/- silences GAS/ISRE-inducible PD-L1 (Graveyard 001)"},
    },
    "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield": {
        "base_intercept": {"logit_value": 0.20, "expression": "0.08 * 2.5", "rationale": "Baseline innate NK/macrophage survival prior"},
        "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO": {"logit_value": 0.35, "expression": "0.08 * 6.875 - 0.08 * 2.5 + _pdelta('resting_b2m_promoter_hour0_density_discount', -0.20)", "rationale": "Dual CD94/NKG2A (HLA-E*01:03 VMAPRTLVL; non-covalent soluble 3D SPR CD94/NKG2A K_D = 2.10 uM vs CD94/NKG2C K_D = 6.52 uM -> K_C/K_A = 3.1048x < 4.0 [theta_C = 0.2436 at C_eff = 2.10 uM] per Lin et al. 2023 PMID:37264229, vs HLA-G leader VMAPRTLFL K_D(NKG2A) = 0.80 +/- 0.055 uM, K_D(NKG2C) = 5.20 +/- 0.20 uM -> K_C/K_A = 6.50x per Kaiser et al. 2008 PMID:18448674 Table 1 [Vales-Gomez et al. 1999 PMID:10428963 Table I reports 2.33 uM vs 10.77 uM = 4.62x for VMAPRTLFL, 7.35 uM for VMAPRTLIL, 11.23 uM for VMAPRTVLL]; covalent SCT K_C/K_A >= 4.0 selectivity and non-agonism on NKG2C+ adaptive NK subsets require Stage-0 Gate 0A-NK SPR verification and CMV/NKG2C+ stratification per Switch S10; DEC-078) + SIRPa (CD47) shield (+0.55) minus resting B2M-promoter Hour-0 density discount (-0.20)"},
        "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO": {"logit_value": 0.35, "expression": "0.08 * 4.375", "rationale": "Retained single HLA-C allele engages KIR2DL and supplies -21M leader peptide (subject to Graveyard 005 ~10-fold lower HLA-C basal expression check after HLA-A/B deletion; PMID:29302013, PMID:37264229)"},
        "ALLO_B2M_CIITA_KO_CD47_LENTI": {"logit_value": 0.30, "expression": "0.08 * 3.75", "rationale": "Constitutive lentiviral CD47 OE (>180k molecules/cell; PMID:37156915)"},
        "ALLO_B2M_KO_PDL1_A20_KI": {"logit_value": -0.75, "expression": "-0.08 * 9.375", "rationale": "B2M-/- without CD47 or HLA-E triggers NK missing-self lysis"},
        "ALLO_NONE": {"logit_value": -0.25, "expression": "-0.08 * 3.125", "rationale": "Unedited allogeneic graft innate/complement vulnerability"},
    },
    "Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience": {
        "base_intercept": {"logit_value": 0.10, "expression": "0.08 * 1.25", "rationale": "Baseline intrinsic stress resilience prior"},
        "CYTO_CXCL10_SWAP_PDL1_SERPINB9": {"logit_value": 0.45, "expression": "0.08 * 5.625", "rationale": "CXCL10 KO + SERPINB9 Granzyme-B inhibition + PD-L1 (PMID:19641142, PMID:32814902)"},
        "CYTO_CONSTITUTIVE_IL10_VIRAL_OE": {"logit_value": -1.00, "expression": "-0.08 * 12.5", "rationale": "Viral IL-10 overexpression failure mode (Graveyard 008)"},
        "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1": {"logit_value": -1.35, "expression": "-0.08 * 16.875", "rationale": "Epistatic promoter-pathway failure (Graveyard 001)"},
        "CYTO_NONE_WITH_VCTX210_A20": {"logit_value": 0.15, "expression": "0.08 * 1.875", "rationale": "TNFAIP3 (A20) NF-kB brake in VCTX210"},
        "CYTO_NONE_DEFAULT": {"logit_value": -0.65, "expression": "-0.08 * 8.125", "rationale": "Unprotected against IFN-g/IL-1b/TNF-a + Granzyme B"},
        "ER_RNLS_KO": {"logit_value": 0.18, "expression": "0.08 * 2.25", "rationale": "RNLS-/- ER stress / thapsigargin protection (PMID:32719542)"},
        "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE": {"logit_value": 0.25, "expression": "0.08 * 3.125", "rationale": "RNLS-/- + non-amyloidogenic IAPP proline edit (PMID:25140605, PMID:18753660)"},
        "ER_NONE_WITH_STAGE7_120UM_NICHE": {"logit_value": 0.25, "expression": "0.08 * 3.125", "rationale": "Stage 7 OXPHOS maturation + <=120 um oxygenation lowers basal ER/UPR stress"},
    },
    "Gate5_72h_Vascular_Oxygen_and_80yr_Safety": {
        "base_intercept": {"logit_value": 0.15, "expression": "0.08 * 1.875", "rationale": "Baseline engraftment & safety prior"},
        "NICHE_TRIHORMONAL_120UM_MVF_CAO2": {"logit_value": 0.45, "expression": "0.08 * 5.625 + _pdelta('stage7_trihormonal_120um_mvf_cao2_niche_prior', 0.45)", "rationale": "<=120 um clusters (<149 um Krogh anoxic limit) + iEC microvessels + 72h CaO2 bridge (PMID:34480863)"},
        "NICHE_NAKED_AVASCULAR_TRANSPLANT": {"logit_value": -0.55, "expression": "-0.08 * 6.875", "rationale": "50-70% 72h avascular ischemic core necrosis"},
        "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION": {"logit_value": -0.95, "expression": "-0.08 * 11.875", "rationale": "Severe avascular diffusion barrier + fibrotic foreign-body capsule"},
        "SAF_AAVS1_CAG_ICASP9_BIALLELIC": {"logit_value": 0.38, "expression": "0.08 * 4.75", "rationale": "Biallelic safe-harbor pan-lineage iCasp9 switch (PMID:35795779)"},
        "SAF_G6PC2_ENDOGENOUS_ICASP9": {"logit_value": -0.45, "expression": "0.08 * -5.625", "rationale": "Lineage-restricted G6PC2::iCasp9 fails to ablate teratoma/SLC18A1+ cells (Graveyard 002)"},
        "SAF_NONE_PEDIATRIC": {"logit_value": -0.45, "expression": "-0.08 * 5.625", "rationale": "No prodrug erasability over 70-80 year pediatric horizon"},
        "SAF_NONE_LADA": {"logit_value": -0.25, "expression": "-0.08 * 3.125", "rationale": "No prodrug erasability over adult LADA horizon"},
    },
    "auxiliary_evaluator_hyperparameters": {
        "signed_gene_signal_slope": {"parameter_value": 0.42, "rationale": "Linear slope multiplying per-gene net_evidence_and_constraint_score inside signed_gene_signal()"},
        "signed_gene_signal_clip_bounds": {"parameter_value": [-4.5, 4.5], "rationale": "Symmetric saturation clamp on cumulative per-branch gene signal"},
        "gene_signal_scales_by_gate_and_branch": {
            "parameter_value": {
                "Gate1_NICHE_INS_GAD2": 0.20,
                "Gate1_ALLO_HLA_C_TYK2_KO": 0.20,
                "Gate1_AUTO_FULL_KO": 0.35,
                "Gate1_ER_RNLS_KO": 0.25,
                "Gate2A_ALLO_B2M_HLAE_CD47": 0.30,
                "Gate2A_ALLO_HLA_C_TYK2": 0.28,
                "Gate2A_ALLO_UP421": 0.28,
                "Gate2A_ALLO_VCTX210": 0.25,
                "Gate2B_ALLO_B2M_NULL": 0.22,
                "Gate2B_HLA_C_INS_EDIT": 0.28,
                "Gate2B_HLA_C_FULL_KO": 0.25,
                "Gate2B_CYTO_CXCL10_SWAP": 0.25,
                "Gate2B_CYTO_IL10_OE": 0.25,
                "Gate3_ALLO_B2M_HLAE_CD47": 0.35,
                "Gate3_ALLO_HLA_C_TYK2": 0.30,
                "Gate3_ALLO_UP421": 0.28,
                "Gate4_CYTO_CXCL10_SWAP": 0.32,
                "Gate4_CYTO_IL10_OE": 0.30,
                "Gate4_CYTO_STAT1_KO": 0.30,
                "Gate4_CYTO_VCTX210_A20": 0.25,
                "Gate4_ER_RNLS_KO": 0.25,
                "Gate5_SAF_ICASP9": 0.25,
            },
            "rationale": "Branch-specific scaling weights coupling empirical target_profiles into each gate logit (F-04/R10 live-wired parameter inventory)",
        },
        "cd47_gsis_penalty_scale": {"parameter_value": -0.25, "rationale": "Multiplier routing CD47 weighted_deleterious_penalty_score into Gate 1 (x0.5 for signaling-attenuated CD47-IgV-GPI vs x1.0 for lentiviral CD47-WT; PMID:37820008, PMID:38508139)"},
        "lrr_shrinkage_hyperparameters": {
            "parameter_value": {"prior_weight": 4.50, "per_empirical_study_slope": 0.50},
            "formula": "signed_lrr_weighted_sum / (w_tot + prior_weight + per_empirical_study_slope * n_emp)",
            "rationale": "Empirical-Bayes-style shrinkage damping single-study outlier LRRs toward zero",
        },
        "effect_direction_weights": {"parameter_value": {"protective": 0.50, "deleterious": 0.75, "secretome_or_differentiation_penalty_addon": 0.35, "shrunk_lrr_multiplier": 1.50}, "rationale": "Asymmetric counter-evidence weighting in build_target_evidence_profiles()"},
        "omics_constraint_weights": {"parameter_value": {"gnomad_pli_ko_slope": -0.45, "gnomad_loeuf_ko_slope": 0.12, "opentargets_t1d_slope": 0.15}, "rationale": "Perturbation-direction-aware human genetic constraint weights"},
        "lower_10pct_cred_bound_hyperparameters": {"parameter_value": {"z_multiplier": 1.2816, "locus_cut_linear_discount": 0.008, "default_gene_se": 0.35, "no_gene_se": 0.45}, "rationale": "Heuristic conservative floor index parameters"},
    },
    "secondary_ranking_tie_break_and_display_constants": {
        "lower_10pct_cred_bound_b_edit_quadratic_discount": {
            "parameter_value": 0.005,
            "formula": "max(0.0, penalized_index * exp(-z_multiplier * se_combined) - locus_cut_linear_discount * K - lower_10pct_cred_bound_b_edit_quadratic_discount * (B ** 2))",
            "rationale": "Secondary display/audit quadratic base-edit floor discount (0.005 * B^2) reported alongside the linear K-only search-grid floor index (Finding 013)",
        },
        "composite_utility_pediatric_weight": {
            "parameter_value": 0.50,
            "rationale": "Equal 50% weight on Pediatric_Early_Onset_T1D in composite_series_system_utility (Finding 013)",
        },
        "composite_utility_adult_weight": {
            "parameter_value": 0.50,
            "rationale": "Equal 50% weight on Adult_Onset_LADA in composite_series_system_utility (Finding 013)",
        },
    },
}


def _sec_param(key: str, default: float) -> float:
    """Returns parameter_value for `key` from COMPLETE_GATE_COEFFICIENT_REGISTRY['secondary_ranking_tie_break_and_display_constants']."""
    sec = COMPLETE_GATE_COEFFICIENT_REGISTRY.get("secondary_ranking_tie_break_and_display_constants", {})
    if isinstance(sec, dict):
        entry = sec.get(key, {})
        if isinstance(entry, dict) and "parameter_value" in entry:
            return _finite_real(entry["parameter_value"], f"Secondary parameter {key}")
    return default


def _composite_utility_weights() -> tuple[float, float]:
    pediatric = _sec_param("composite_utility_pediatric_weight", 0.50)
    adult = _sec_param("composite_utility_adult_weight", 0.50)
    if not 0 <= pediatric <= 1 or not 0 <= adult <= 1 or not math.isclose(pediatric + adult, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("Composite cohort weights must be nonnegative and sum to one")
    return pediatric, adult


def _aux_param(key: str, default):
    """Returns parameter_value for `key` from COMPLETE_GATE_COEFFICIENT_REGISTRY['auxiliary_evaluator_hyperparameters'] so mutations propagate live (Audit R10)."""
    aux = COMPLETE_GATE_COEFFICIENT_REGISTRY.get("auxiliary_evaluator_hyperparameters", {})
    if isinstance(aux, dict):
        entry = aux.get(key, {})
        if isinstance(entry, dict) and "parameter_value" in entry:
            return entry["parameter_value"]
    return default


def _branch_scale(branch_key: str, default_scale: float) -> float:
    scales = _aux_param("gene_signal_scales_by_gate_and_branch", {})
    if isinstance(scales, dict) and branch_key in scales:
        return _finite_real(scales[branch_key], f"Gene-signal branch scale {branch_key}")
    return default_scale


def _pdelta(
    key: str,
    default: float,
    gate_key: str | None = None,
    entry_key: str | None = None,
    reg_default: float | None = None,
) -> float:
    """
    Returns (point_estimate - default) from ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE plus any active
    PURE_EMPIRICAL_LINEAGE wet-lab posterior shift so prior-table and wet-lab calibration mutations
    dynamically propagate into gate logits (Finding 013/081/082).
    """
    entry = ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.get(key, {})
    w_deltas = _load_active_wetlab_posterior_deltas()
    return (_finite_real(entry.get("point_estimate", default), f"Engineering prior {key}") - default) + _finite_real(w_deltas.get(key, 0.0), f"Wet-lab shift {key}")


def _rdelta(gate_key: str, entry_key: str, default_val: float) -> float:
    """Returns (logit_value - default_val) from COMPLETE_GATE_COEFFICIENT_REGISTRY plus any active PURE_EMPIRICAL_LINEAGE wet-lab posterior shift."""
    w_deltas = _load_active_wetlab_posterior_deltas()
    w_shift = _finite_real(w_deltas.get(f"{gate_key}::{entry_key}", 0.0), "Wet-lab gate shift")
    g_dict = COMPLETE_GATE_COEFFICIENT_REGISTRY.get(gate_key, {})
    if isinstance(g_dict, dict):
        e_dict = g_dict.get(entry_key, {})
        if isinstance(e_dict, dict) and "logit_value" in e_dict:
            return (_finite_real(e_dict["logit_value"], f"Gate coefficient {gate_key}.{entry_key}") - default_val) + w_shift
    return w_shift

ALLO_CLOAK_MODULES = {
    "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO": {
        "label": "B2M-ISRE::[EEF1A1-core+intron]-HLA-E*01:03(VMAPRTLVL)-GSG-P2A-CD47-IgV-GPI-WPREmut6 (or 2A-Free Bidirectional Arm 1B) + CIITA-/-",
        "cuts": 2,
        "base_edits": 0,
        "genes": ["B2M", "CIITA", "HLA-E", "CD47"],
        "description": "2 chromosomal cuts (B2M exon 1 knock-in via single-AAV6 [Arm 1B 2A-free bidirectional: 3,974 bp donor / 4,264 bp with ITRs <= 4,700 bp limit] or non-viral linear dsDNA/ceDNA HDR [Arm 1A bicistronic: 4,880 bp donor / 5,170 bp with ITRs; Shy et al. 2023 PMID:36008610 supports a different nonviral editing context; exact ceDNA/Arm 1A delivery in hPSCs remains unvalidated; note that bacterial IS110 Bridge-RNA recombinases in Durrant et al. 2024 PMID:38926615 and Pelea et al. 2026 PMID:41642947 are exploratory and not validated in human iPSCs] + CIITA biallelic KO). Eliminates classical HLA-A/B/C and HLA-II direct presentation while co-expressing single-chain HLA-E*01:03 (HLA-A2 VMAPRTLVL leader; ~12% overlap with clinical HCMV UL40 isolates per Lin et al. 2023 PMID:37264229, vs optional cytoplasmic KW-motif retention-signal removal per He et al. 2025 PMC12582296) and CD47-IgV-GPI placed at the 3' cistron so the 37-aa DAF GPI signal resides at the true extreme C-terminus without a 2A peptide scar (note membrane topology per Audit R01: because HLA-E*01:03's C-terminal cytoplasmic tail resides in the cytosol whereas Furin's catalytic domain is lumenal inside the trans-Golgi network, Furin is omitted and HLA-E*01:03 retains the 21-aa cytosolic -GSG-P2A peptide tail in Arm 1A, compared head-to-head in Stage-0 Decision Experiment #1 against scarless 2A-free bidirectional Arm 1B [SV40pA-CD47-IgV-GPI(antisense) <- EEF1A1-core | Syn_5xISRE_minP -> HLA-E*01:03(VMAPRTLVL)-bGHpA]; PMID:38508139, PMID:37156915, PMID:37820008). Subject to Hour-0 resting B2M-promoter surface density verification (Graveyard 009; Han et al. 2019 PMID:31040209).",
    },
    "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO": {
        "label": "HLA-A-/- + HLA-B-/- + CIITA-/- (Retaining 1 matched HLA-C allele) + TYK2-/-",
        "cuts": 4,
        "base_edits": 0,
        "genes": ["HLA-A", "HLA-B", "CIITA", "HLA-C", "TYK2"],
        "description": "4 chromosomal cuts (HLA-A KO + HLA-B KO + CIITA KO + TYK2 KO). Retains 1 recipient-matched HLA-C allele (Xu et al. 2019 PMID:30853558) for KIR2DL C1/C2 NK self-recognition and viral surveillance while blocking IFN-alpha-induced MHC-I hyper-expression via TYK2-/- (PMID:36289205); note that whereas the -21M/-21T leader dimorphism is primarily an HLA-B polymorphism and HLA-C typically encodes -21M (Ramsuran et al. 2018 PMID:29302013), deleting HLA-A and HLA-B removes the dominant HLA-A/-21M HLA-B leader supply so endogenous HLA-E stabilization must be verified empirically (Lin et al. 2023 PMID:37264229), and co-targeting linked HLA-A (~chr6:29.94 Mb) and HLA-B (~chr6:31.35 Mb) in cis carries a ~1.41 Mb interstitial deletion/LOH hazard across the intervening HLA-C locus (~31.27 Mb) unless executed sequentially or verified by long-read/ddPCR copy-number gating.",
    },
    "ALLO_B2M_CIITA_KO_CD47_LENTI": {
        "label": "B2M-/- + CIITA-/- + Constitutive Lentiviral CD47 OE (Sana UP421 / HIP benchmark)",
        "cuts": 3,
        "base_edits": 0,
        "genes": ["B2M", "CIITA", "CD47"],
        "description": "2 CRISPR nuclease DSB knockouts (B2M-/-, CIITA-/-) + 1 semi-random lentiviral CD47 integration (counted as K_total_genomic_modifications = 3 in the uniform locus budget, with K_crispr_dsb = 2 and K_lentiviral_integration = 1; Finding 025). Validated in NHP and n=1 human forearm intramuscular sentinel implant (PMID:37156915, PMID:40757665; 14-month follow-up correspondence in Carlsson et al. 2026 PMID:42430558, which carries zero quantitative weight in target_profiles because PubMed indexes NEJM correspondence without an abstract).",
    },
    "ALLO_B2M_KO_PDL1_A20_KI": {
        "label": "B2M-/- + PD-L1 KI + TNFAIP3 (A20) KI (CRISPR/ViaCyte VCTX210 benchmark, CIITA WT)",
        "cuts": 3,
        "base_edits": 0,
        "genes": ["B2M", "CD274", "TNFAIP3"],
        "description": "3 genomic edits (B2M KO + PD-L1 KI + TNFAIP3/A20 KI) without CIITA KO or CD47/HLA-E NK shield (VCTX210 clinical trial NCT05210530).",
    },
    "ALLO_NONE_UNEDITED_OR_SYSTEMIC_IMMUNOSUPPRESSION": {
        "label": "Unedited Wild-Type Allogeneic HLA (Vertex VX-264 macrodevice or zimislecel / VX-880 with chronic systemic immunosuppression)",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "0 genetic edits (K=0). Represents the unedited allogeneic SC-islet baseline: inside the immunosuppression-free evaluator, unedited allogeneic cells require a physical macro-encapsulation barrier (VX-264, NCT05791201; terminated March 2025 due to insufficient C-peptide across diffusion barrier), whereas clinical zimislecel (VX-880, NCT04786262 / NCT06832410) achieves insulin independence only via chronic systemic calcineurin/mTOR/anti-thymocyte immunosuppression outside the immunosuppression-free target product profile (Finding 025).",
    },
}

AUTOANTIGEN_MODULES = {
    "AUTO_NONE": {
        "label": "No Intracellular Autoantigen Edits (Lean Default for Classical-MHC-Null Chassis A)",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "Zero autoantigen edits. Eliminates direct classical HLA-A/B/C CD8+ autoantigen presentation when paired with B2M-/-, avoiding G6PC2-/- fasting glucose set-point downshift and SLC30A8-/- dense-core defects; residual indirect APC presentation is modeled explicitly.",
    },
    "AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT": {
        "label": "Base-Edit of INS Signal Peptide (PPI) + Silent INS-DRiP AUG Removal",
        "cuts": 0,
        "base_edits": 1,
        "genes": ["INS"],
        "description": "1 multiplex base-edit locus (B=1) at INS removing out-of-frame INS-DRiP initiation and/or replacing PPI signal peptide (PMID:28263308, PMID:42690930).",
    },
    "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR": {
        "label": "Gross Pan-Autoantigen/GWAS KO (G6PC2-/- + SLC30A8-/- + PTPN2-/-, 3 cuts) + 4-Locus Base/Prime Edits (INS, GAD2, IAPP, CXADR)",
        "cuts": 3,
        "base_edits": 4,
        "genes": ["G6PC2", "SLC30A8", "INS", "GAD2", "IAPP", "CXADR", "PTPN2"],
        "description": "High-burden 7-target autoantigen/GWAS modification: K=3 chromosomal DSB knockouts (G6PC2-/-, SLC30A8-/-, and high-pLI negative-feedback brake PTPN2-/- [gnomAD v4 pLI = 0.990]) + B=4 nickase base/prime-edit loci (INS, high-pLI paracrine GABA enzyme GAD2 [gnomAD v4 pLI = 1.000], IAPP, CXADR). Incurs G6PC2-/- fasting glucose set-point downshift hazard (PMID:32569842, PMID:34954144), Graveyard 007 high-pLI KO/coding-disruption block on PTPN2-/- and multiplex GAD2 disruption (PMID:19336676, PMID:32432213), and multiplex nickase DSB/translocation penalty (PMID:39529170).",
    },
}

CYTOKINE_GRANZYME_MODULES = {
    "CYTO_CXCL10_SWAP_PDL1_SERPINB9": {
        "label": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1 (or CXCL10::PD-L1-P2A-SERPINB9 Endogenous CXCL10 KO + Inflammation-Inducible Swap-In)",
        "cuts": 1,
        "base_edits": 0,
        "genes": ["CXCL10", "CD274", "SERPINB9"],
        "description": "1 chromosomal cut at CXCL10 exon 1: deletes graft CXCL10 chemokine and expresses Granzyme-B inhibitor SERPINB9 + PD-L1 under native IFN/NF-kB inflammation control (comparing inverted SERPINB9-GSG-tPT2A-PD-L1 against PD-L1-GSG-tPT2A-SERPINB9 in Decision Experiment #3; Yoshihara et al. 2020 PMID:32814902; Sintov et al. 2022 PMID:36055241; Parent et al. 2021 PMID:34407395; Liu et al. 2017 PMID:28526819).",
    },
    "CYTO_CONSTITUTIVE_IL10_VIRAL_OE": {
        "label": "Constitutive Viral IL-10 Overexpression (1 locus — Graveyard 008 Mast-Cell/Granulocyte Trap)",
        "cuts": 1,
        "base_edits": 0,
        "genes": ["IL10"],
        "description": "1 viral/transgenic locus constitutively secreting IL-10. Fails to provide durable alloprotection while risking local mast-cell/granulocyte expansion (mechanistic mast-cell expansion context: Polukort et al. 2016 PMID:27183617).",
    },
    "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1": {
        "label": "Clonal STAT1-/- + GAS/ISRE-Inducible PD-L1 (Epistatic Trap — Graveyard 001)",
        "cuts": 1,
        "base_edits": 0,
        "genes": ["STAT1", "CD274"],
        "description": "1 cut (STAT1 KO) paired with a STAT1-dependent GAS/ISRE promoter (silences inducible PD-L1; Sintov et al. 2022 PMID:36055241).",
    },
    "CYTO_NONE": {
        "label": "No Intrinsic Cytokine or Granzyme Shield",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "Zero intrinsic protection against IFN-gamma/IL-1beta/TNF-alpha or Granzyme B bystander killing.",
    },
}

ER_AMYLOID_MODULES = {
    "ER_NONE": {
        "label": "No Intrinsic ER/Amyloid Gene Edit (Relies on Stage 7 Metabolic Maturation)",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "0 extra cuts. Avoids the human iPSC-beta RNLS-/- elevated basal OCR / metabolic precaution prior (noting that poor dynamic GSIS persisted in both WT and RNLS KO in MacDonald et al. 2025 Fig. 6, PMID:39988068).",
    },
    "ER_RNLS_KO": {
        "label": "Biallelic RNLS-/- Knockout (1 cut)",
        "cuts": 1,
        "base_edits": 0,
        "genes": ["RNLS"],
        "description": "1 chromosomal cut (RNLS KO). Protects beta-cells against thapsigargin ER stress and CTL killing in mouse NIT-1 / SC-beta screens (Cai et al. 2020 PMID:32719542), while carrying a conservative human iPSC-beta basal OCR / metabolic precaution prior (poor dynamic GSIS persisted in both WT and RNLS KO in MacDonald et al. 2025 Fig. 6, PMID:39988068).",
    },
    "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE": {
        "label": "RNLS-/- (1 cut) + IAPP (A25P/S28P/S29P) Prime-Edit (1 base/prime edit)",
        "cuts": 1,
        "base_edits": 1,
        "genes": ["RNLS", "IAPP"],
        "description": "1 cut (RNLS KO) + 1 prime-edit locus (IAPP proline conversion, PMID:18753660).",
    },
    "ER_IAPP_PRAMLINTIDE_RNLS_WT": {
        "label": "RNLS+/+ (WT) + Non-Amyloidogenic IAPP (A25P/S28P/S29P) Prime-Edit (0 cuts, 1 prime edit)",
        "cuts": 0,
        "base_edits": 1,
        "genes": ["IAPP"],
        "description": "0 chromosomal DSB cuts + 1 scarless prime-edit converting amyloidogenic human IAPP20-29 to non-amyloidogenic proline-substituted pramlintide-like IAPP(A25P/S28P/S29P) while preserving wild-type RNLS metabolic fitness (eliminating the human iPSC-beta RNLS-/- OCR/GSIS precaution prior in MacDonald et al. 2025 PMID:39988068 while directly blocking human IAPP oligomer toxicity in primary/SC-islets per Potter et al. 2010 PMID:20160085 and Thorngren et al. 2024 PMID:38603470, and noting that A25P/S28P/S29P converts soluble monomers while seeding/conformational effects require empirical Thioflavin-T/Congo-red qualification per Wang et al. 2014 PMID:25140605; AUDIT-REF & AUDIT-REF).",
        "molecular_and_immunological_qualification_gate_kimi31_r4": {
            "mature_iapp_residues_modified": "Mature human IAPP residues A25P, S28P, S29P (precursor preproIAPP positions 58, 61, 62 on chr12p12.1)",
            "required_assays": [
                "pegRNA spacer/PBS/RTT sequence registration and off-target rhAmpSeq/amplicon-seq verification",
                "Intact-mass LC-MS verification of C-terminally amidated mature 37-aa IAPP(A25P/S28P/S29P) processing by PC1/3, PC2, CPE, and PAM",
                "Thioflavin-T and Congo-red amyloid fibril inhibition assay under 7-day 16.7 mM high-glucose challenge (accounting for residual seeding kinetics in Wang et al. 2014 PMID:25140605)",
                "Host HLA-DQ2/DQ8/DR4 dendritic-cell indirect presentation screen confirming no de novo neo-epitope CD4+ T-cell priming from the 3-proline substitution",
            ],
            "gate4_scoring_note": "In the 46-entry structural prior registry, ER_IAPP_PRAMLINTIDE_RNLS_WT shares the +0.25 NICHE_TRIHORMONAL_120UM_MVF_CAO2 Gate-4 baseline when paired with the Stage-7 120-um niche while avoiding ER_RNLS_KO's -0.38 Gate-1 precaution prior (incurring only the B=1 prime-edit factor 0.9650); empirical promotion of the IAPP base edit requires positive ΔGate4 in Stage-0 high-glucose amyloid challenge.",
        },
    },
}

NICHE_VASCULAR_MODULES = {
    "NICHE_TRIHORMONAL_120UM_MVF_CAO2": {
        "label": "Stage 7 Tri-Hormonal (60/30/10, <=120 um) + iEC Microvessels + 14d Pre-Vasc Pouch + 72h PDMS-CaO2/Catalase",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "Vascularized (non-retrievable without surgical excision / iCasp9 ablation) micro-pocket with <=120 um clusters (R <= 60 um, strictly below both the conservative rat islet 74.6 um R_crit and the Stage 7 human SC-islet R_crit = sqrt(6*D*C_0/V_max) = 92.6 um [185.2 um diameter] at V_max = 0.034 mol/m^3/s; PMID:35241836, PMID:34480863) and hydrophobic PDMS-CaO2 (OxySite) + catalase oxygen bridge (noting that 0.7915 g [791.5 mg / 10.981 mmol] is the total loaded pure CaO2 inventory across 200 cm^2 [of which 0.4116 g / 52.0% hydrolyzes over the 72-h window at t_1/2 = 68 h, embedded in 3.166 g / 2.533 g / 1.706 g total PDMS-CaO2 composite mass at 25.0 / 31.25 / 46.39 wt% for Layouts A / B / Design X; delivering 0.4116 g reacted CaO2 at an empirical 72-h utilization fraction f_72h = 0.45 requires 0.915 g pure CaO2 [2.033 g composite at 45 wt%], whereas Pedraza et al. 2012 PMID:22371586 validated multi-week bulk PDMS-CaO2 disk release [f_72h ~= 0.0888 -> 4.635 g pure CaO2 for 0.4116 g reacted, or 8.91 g if scaling the full 0.7915 g inventory by 1/0.0888] and Montazeri et al. 2016 PMID:26970510 used PVP/H2O2-PLGA microparticles, so 72-h burst-matched <=150 um PDMS-CaO2 kinetics require Stage-1 empirical calibration; DEC-078]).",
    },
    "NICHE_NAKED_AVASCULAR_TRANSPLANT": {
        "label": "Standard Avascular Islet Cluster Transplant (~150-250 um, No Pre-Vasc or O2 Bridge)",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "Acute avascular injection suffering 50-70% 72h ischemic core loss.",
    },
    "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION": {
        "label": "Avascular Cell-Impermeable Macro-Encapsulation Device (VX-264 / PEC-Encap style)",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "Retrievable physical barrier (>150-200 um diffusion distance) causing central hypoxia and fibrotic foreign-body reaction.",
    },
}

SAFETY_SWITCH_MODULES = {
    "SAF_AAVS1_CAG_ICASP9_BIALLELIC": {
        "label": "Biallelic AAVS1-[455bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE)]-EEF1A1-core-iCasp9 (or AAVS1-CAG-iCasp9 / GAPDH-P2A-iCasp9) Pan-Lineage Safe-Harbor Switch (1 locus)",
        "cuts": 1,
        "base_edits": 0,
        "genes": ["CASP9"],
        "description": "1 safe-harbor locus knocked in biallelically (3,457 bp donor / 3,747 bp with ITRs <= 4,700 bp single-AAV6 limit; insulated with 455 bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE alias; insulation function pending) + cHS4 boundaries driving EEF1A1-core-iCasp9; Muller-Kuller et al. 2015 PMID:25605798; Zhang et al. 2017 PMID:28860464; Pfaff et al. 2013 PMID:23307570; Ordovas et al. 2015 PMID:26455413) so undifferentiated iPSCs, endocrine cells, and off-target SLC18A1+ cells are erasable by AP1903/rimiducid (<10^-8 target escape bound per colony-forming proliferating cell; Wunderlich et al. 2022 PMID:35795779).",
    },
    "SAF_G6PC2_ENDOGENOUS_ICASP9": {
        "label": "G6PC2::iCasp9 Mature-Beta-Restricted Switch (1 locus — Graveyard 002 Escape Risk)",
        "cuts": 1,
        "base_edits": 0,
        "genes": ["CASP9", "G6PC2"],
        "description": "1 locus (G6PC2 KO + iCasp9 KI). Fails to ablate dedifferentiated, teratoma, or SLC18A1+ off-target cells (Yagyu 2016 PMID:27626039).",
    },
    "SAF_NONE": {
        "label": "No Suicide Switch",
        "cuts": 0,
        "base_edits": 0,
        "genes": [],
        "description": "Zero pharmacological erasability in the event of viral infection or oncogenic transformation.",
    },
}


def load_inputs():
    prepare_wetlab_evaluation_context()
    with open(GRAPH_PATH, "r", encoding="utf-8") as f:
        graph = json.load(f)
    studies = []
    for p in sorted(EXTRACTED_DIR.glob("*.json")):
        if not p.stem.isdigit():
            continue
        with open(p, "r", encoding="utf-8") as f:
            studies.append(json.load(f))
    with open(OMICS_DIR / "target_prior_table.json", "r", encoding="utf-8") as f:
        omics_raw = json.load(f)
    omics = {t["gene_symbol"]: t for t in omics_raw["targets"]}
    return graph, studies, omics


def get_endpoint_orientation_sign(quote_and_var: str) -> float:
    """
    Returns -1.0 if the measured variable in the text is a lower-is-better hazard metric
    (e.g., specific lysis, cell death, apoptosis, rejection, hyperglycemia) without a
    positive survival/viability/GSIS term, so that a negative LRR (reduction in lysis)
    maps to a positive welfare gain, while multiplying lrr by -1 (Sol's sign_flip probe)
    always negates the resulting welfare contribution.
    """
    ql = quote_and_var.lower()
    is_hazard = any(
        k in ql
        for k in (
            "specific lysis",
            "cytotoxicity",
            "apoptosis",
            "cell death",
            "apoptotic",
            "necrotic core",
            "blood glucose",
            "hyperglycemia",
            "amyloid area",
            "thioflavin",
        )
    )
    has_pos_metric = any(
        k in ql
        for k in (
            "survival",
            "viability",
            "gsis",
            "stimulation index",
            "c-peptide",
            "insulin secretion",
            "normoglycemia",
            "vascular density",
        )
    )
    if is_hazard and not has_pos_metric:
        return -1.0
    return 1.0


def strict_gene_match(sym: str, text: str) -> bool:
    """Word-boundary gene-symbol match (prevents INS/CAR/HLA substring false positives)."""
    s = sym.upper()
    t = text.upper()
    if s == "CAR":
        return bool(re.search(r"(?<![A-Z0-9])CAR(?![A-Z0-9])", t)) and "COXSACKIE" in t
    if s == "INS":
        return bool(
            re.search(
                r"(?<!0\.005\s)\[INS\]|\bINS\s*\((?:SIGNAL|BASE|CRISPR|KO|KNOCK)|\bINS\s*(?:KO\b|KNOCKOUT\b|KNOCK-OUT\b|-/-|−/−|BASE[ -]EDIT|LEADER\b|SIGNAL\s+PEPTIDE|DRIP\b|3'UTR|FUNCTIONAL\s+EVALUATION)|^INS$",
                t,
            )
        )
    if s.startswith("HLA-"):
        return bool(re.search(r"(?<![A-Z0-9])" + re.escape(s) + r"(?:W(?=[\d*])|(?![A-Z]))", t))
    return bool(re.search(r"(?<![A-Z0-9])" + re.escape(s) + r"(?![A-Z0-9])", t))


def MATCH(g: str, fx: str, fx_upper: str) -> bool:
    return strict_gene_match(g, fx)


def AMATCH(alias: str, fx: str, fx_upper: str) -> bool:
    return strict_gene_match(alias, fx)


def expand_gene_symbol_tokens(raw_g: str) -> list[str]:
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


def build_target_evidence_profiles(studies: list[dict], omics: dict) -> dict:
    """
    Computes sign-sensitive, counter-evidence-penalized, omics-constrained, and
    non-broadcasting target profiles.
    - Attribution Rule (Audit §5.2.1): Effects are attributed ONLY to genes
      explicitly named in `focal_factor_X` or `verbatim_evidence_quote`, or to the
      sole edited gene if `len(study_genes) == 1`. Multi-gene studies never broadcast
      unmatched effects across all study genes.
    """
    if not isinstance(studies, list) or any(not isinstance(study, dict) for study in studies) or not isinstance(omics, dict):
        raise TypeError("Evidence profiles require a study list and omics mapping")
    # These declared masses, weights and dispersions cannot be negative or
    # truthy-coerced: a negative harm weight reverses adverse evidence, while
    # NaN can otherwise leave a seemingly complete profile with unknown scores.
    def profile_parameters(group: str, defaults: dict, *, nonnegative: bool = True) -> dict:
        raw = _aux_param(group, {})
        if not isinstance(raw, dict):
            raise ValueError(f"Profile parameter group {group} must be a mapping")
        values = {key: _finite_real(raw.get(key, default), f"{group}.{key}")
                  for key, default in defaults.items()}
        if nonnegative and any(value < 0 for value in values.values()):
            raise ValueError(f"Profile parameter group {group} requires nonnegative values")
        return values

    direction_weights = profile_parameters("effect_direction_weights", {
        "deleterious": 0.75, "secretome_or_differentiation_penalty_addon": 0.35,
        "protective": 0.50, "shrunk_lrr_multiplier": 1.50,
    })
    shrink_parameters = profile_parameters("lrr_shrinkage_hyperparameters", {
        "prior_weight": 4.50, "per_empirical_study_slope": 0.50,
    })
    omics_weights = profile_parameters("omics_constraint_weights", {
        "gnomad_pli_ko_slope": -0.45, "gnomad_loeuf_ko_slope": 0.12,
        "opentargets_t1d_slope": 0.15,
    }, nonnegative=False)
    dispersion_parameters = profile_parameters("lower_10pct_cred_bound_hyperparameters", {
        "default_gene_se": 0.35,
    })

    def accumulate(st: dict, key: str, amount: float) -> None:
        # Finite input coefficients can overflow after repeated additions.
        # Reject the arithmetic failure instead of emitting an infinite score.
        st[key] = _finite_real(st[key] + amount, f"{st['gene_symbol']}.{key} accumulation")

    source_types = {"FullText_XML_or_PDF", "FullText_Plus_Supplement_or_GEO", "Abstract_Only", "Project_Summary", "Preprint_Superseded_By_Peer_Reviewed", "Review_or_Commentary"}
    pmids = []
    for study in studies:
        metadata = study.get("study_metadata")
        if not isinstance(metadata, dict) or not isinstance(metadata.get("pmid"), str) or not metadata["pmid"]:
            raise ValueError("Every evidence study requires an explicit PMID identifier")
        if metadata.get("source_text_type") not in source_types:
            raise ValueError(f"Unmapped source_text_type in PMID {metadata['pmid']}")
        effects = study.get("atomic_effects", [])
        if not isinstance(effects, list) or any(not isinstance(effect, dict) for effect in effects):
            raise ValueError("Atomic evidence must be a list of explicit effect records")
        effect_ids = [effect.get("effect_id") for effect in effects]
        if any(not isinstance(eid, str) or not eid for eid in effect_ids) or len(set(effect_ids)) != len(effect_ids):
            raise ValueError("Effect identifiers must be nonempty and unique within each study")
        for effect in effects:
            if effect.get("effect_direction") not in {"Protective_or_Improved", "Neutral_or_NonInferior", "Deleterious_or_Impaired"}:
                raise ValueError("Unknown effect_direction cannot be interpreted as neutral evidence")
            if type(effect.get("secretome_or_differentiation_penalty_observed", False)) is not bool:
                raise ValueError("Evidence penalty flags must be boolean")
            for key in ("focal_factor_X", "verbatim_evidence_quote", "outcome_variable_Y"):
                if key in effect and not isinstance(effect[key], str):
                    raise ValueError(f"Evidence {key} must be a string")
        pmids.append(metadata["pmid"])
    if len(pmids) != len(set(pmids)):
        raise ValueError("Duplicate study PMID would double-count empirical evidence")
    for gene, annotation in omics.items():
        if not isinstance(annotation, dict):
            raise ValueError(f"Omics annotation for {gene} must be a dictionary")
        for key in ("gnomad_pli", "gnomad_loeuf_upper", "opentargets_t1d_score", "curated_literature_prior_t1d_score"):
            value = annotation.get(key)
            if value is not None:
                number = _finite_real(value, f"{gene}.{key}")
                if number < 0 or (key != "gnomad_loeuf_upper" and number > 1):
                    raise ValueError(f"{gene}.{key} is outside its declared domain")
    target_stats: dict[str, dict] = {}

    def ensure_gene(g: str):
        if g not in target_stats:
            om = omics.get(g, {})
            rec_edit = str(om.get("recommended_edit") or "")
            is_gof = (
                g in {"CD47", "HLA-E", "CD274", "SERPINB9", "CASP9", "TNFAIP3", "IL10"}
                or ("knock-in" in rec_edit.lower() and "ko" not in rec_edit.lower())
            )
            target_stats[g] = {
                "gene_symbol": g,
                "gnomad_pli": om.get("gnomad_pli"),
                "gnomad_loeuf_upper": om.get("gnomad_loeuf_upper"),
                "opentargets_t1d_score": om.get("opentargets_t1d_score"),
                "curated_literature_prior_t1d_score": om.get("curated_literature_prior_t1d_score"),
                "opentargets_prior_source": om.get("opentargets_score_source", om.get("opentargets_score_provenance", om.get("opentargets_prior_source", "Unmapped_or_Zero_T1D_Association"))),
                "is_engineered_genetic_target": bool(om.get("is_engineered_genetic_target", g not in {"GCG", "SST", "PECAM1", "GCK", "PCSK1"})),
                "functional_gate": om.get("functional_gate"),
                "recommended_edit": om.get("recommended_edit"),
                "default_perturbation_direction": "Gain_of_Function_or_KnockIn" if is_gof else "Knockout_or_BaseEdit",
                "arm_edit_modalities_observed": [],
                "n_studies_total": 0,
                "n_empirical_studies": 0,
                "n_project_summary_records_excluded_from_weights": 0,
                "n_review_records_excluded_from_weights": 0,
                "n_superseded_preprint_records_excluded_from_weights": 0,
                "n_human_tier1_to_3_studies": 0,
                "n_rodent_tier4_studies": 0,
                "n_protective_effects": 0,
                "n_deleterious_or_penalty_effects": 0,
                "n_neutral_effects": 0,
                "n_ko_lof_effects": 0,
                "n_gof_oe_effects": 0,
                "n_associational_or_untyped_effects": 0,
                "weighted_protective_score": 0.0,
                "weighted_deleterious_penalty_score": 0.0,
                "ko_weighted_protective_score": 0.0,
                "ko_weighted_deleterious_score": 0.0,
                "gof_weighted_protective_score": 0.0,
                "gof_weighted_deleterious_score": 0.0,
                "signed_lrr_weighted_sum": 0.0,
                "signed_lrr_weight_total": 0.0,
                "ko_signed_lrr_weighted_sum": 0.0,
                "ko_signed_lrr_weight_total": 0.0,
                "gof_signed_lrr_weighted_sum": 0.0,
                "gof_signed_lrr_weight_total": 0.0,
                "supporting_pmids": [],
                "counter_evidence_pmids": [],
            }

    for g in omics:
        ensure_gene(g)

    alias_map = {
        "ZN-T8": "SLC30A8",
        "ZNT8": "SLC30A8",
        "IGRP": "G6PC2",
        "GAD65": "GAD2",
        "IA-2": "PTPRN",
        "PD-L1": "CD274",
        "PDL1": "CD274",
        "A20": "TNFAIP3",
        "PI-9": "SERPINB9",
        "MDA5": "IFIH1",
        "CAR": "CXADR",
        "ICASP9": "CASP9",
        "IND-CASP9": "CASP9",
        "IL-10": "IL10",
    }
    for canonical in sorted(set(alias_map.values())):
        ensure_gene(canonical)

    _expand_gene_symbol = expand_gene_symbol_tokens

    # Pre-register all study_genes across studies in sorted order so target_stats is 100% invariant to input study order (AT-03)
    all_study_genes = sorted({
        g_atom
        for doc in studies
        for arm in doc.get("arms", [])
        for ed in arm.get("gate2_3_4_genetic_and_epitope_edits", [])
        for g_atom in _expand_gene_symbol(str(ed.get("gene_symbol_hgnc") or ""))
    })
    for g in all_study_genes:
        ensure_gene(g)

    ordered_studies = sorted(
        studies, key=lambda d: str(d.get("study_metadata", {}).get("pmid", ""))
    )
    for doc in ordered_studies:
        meta = doc.get("study_metadata", {})
        pmid = str(meta.get("pmid", ""))
        src_type = meta.get("source_text_type", "")
        is_project_summary = src_type == "Project_Summary"
        is_review = src_type == "Review_or_Commentary"
        is_superseded_preprint = (
            src_type == "Preprint_Superseded_By_Peer_Reviewed"
            or bool(meta.get("superseded_by_pmid"))
        )
        is_non_weighted_record = is_project_summary or is_superseded_preprint or is_review
        tier = meta.get("max_biological_fidelity_tier", "Tier_4_Rodent_or_Xenogeneic_Only")
        if tier not in FIDELITY_WEIGHTS:
            raise ValueError(f"Unmapped max_biological_fidelity_tier {tier!r} in PMID {pmid}")
        w = FIDELITY_WEIGHTS[tier]
        is_human = tier != "Tier_4_Rodent_or_Xenogeneic_Only"

        study_genes = set()
        study_gene_edit_types: dict[str, set[str]] = {}
        for arm in doc.get("arms", []):
            for ed in arm.get("gate2_3_4_genetic_and_epitope_edits", []):
                raw_g = str(ed.get("gene_symbol_hgnc") or "")
                for g in _expand_gene_symbol(raw_g):
                    study_genes.add(g)
                    ensure_gene(g)
                    et = str(ed.get("edit_modality") or ed.get("edit_type") or "").strip()
                    if et:
                        study_gene_edit_types.setdefault(g, set()).add(et)
                        if et not in target_stats[g]["arm_edit_modalities_observed"]:
                            target_stats[g]["arm_edit_modalities_observed"].append(et)

        for g in study_genes:
            st = target_stats[g]
            st["n_studies_total"] += 1
            if is_project_summary:
                st["n_project_summary_records_excluded_from_weights"] += 1
            elif is_superseded_preprint:
                st["n_superseded_preprint_records_excluded_from_weights"] += 1
            elif is_review:
                st["n_review_records_excluded_from_weights"] += 1
            else:
                st["n_empirical_studies"] += 1
            if is_human:
                st["n_human_tier1_to_3_studies"] += 1
            else:
                st["n_rodent_tier4_studies"] += 1
            if pmid and not is_non_weighted_record and pmid not in st["supporting_pmids"]:
                st["supporting_pmids"].append(pmid)

        study_effs = doc.get("atomic_effects", [])
        n_effs_in_study = max(1, len(study_effs))
        # Project_Summary (42430558) and Preprint_Superseded_By_Peer_Reviewed (40832357 -> 41202120) are excluded from quantitative evidence weights (A02/A08/A10 + AUDIT-REF)
        eff_weight = 0.0 if is_non_weighted_record else (0.12 * w) / float(n_effs_in_study)

        marker_readout_genes = {"INS", "GCG", "SST", "CHGA", "NKX6-1", "PDX1", "MAFA"}
        for eff in study_effs:
            focal_only = eff.get("focal_factor_X") or ""
            focal_only_upper = focal_only.upper()
            q_str = eff.get("verbatim_evidence_quote") or ""
            for field_name, field_text in (("focal_factor_X", focal_only), ("verbatim_evidence_quote", q_str)):
                if "[REDACTED_NON_OA_EXCERPT:" in field_text:
                    declared_field = re.search(r"\bfield=([A-Za-z0-9_]+)", field_text)
                    if not validate_public_excerpt_token(field_text, expected_pmid=pmid) or declared_field is None or declared_field[1] != field_name:
                        raise ValueError(f"{pmid}:{eff.get('effect_id')}: invalid source-bound opaque {field_name}")
            fx = focal_only + " " + q_str
            fx_upper = fx.upper()
            redacted_genes_matches = re.findall(r"\bgenes=([A-Za-z0-9_,-]+)", fx) if "[REDACTED_NON_OA_EXCERPT:" in fx else []
            if redacted_genes_matches:
                g_tag_val = redacted_genes_matches[-1]
                matched_genes = set()
                if g_tag_val != "NONE":
                    for g_tok in g_tag_val.split(","):
                        g_tok = g_tok.strip()
                        if g_tok:
                            # Opaque metadata permits the same public aliases
                            # as source text. Preserve the target identity and
                            # deduplicate aliases before accumulating evidence.
                            g_tok = alias_map.get(g_tok, g_tok)
                            ensure_gene(g_tok)
                            matched_genes.add(g_tok)
            else:
                matched_genes = set()
                for g in list(target_stats.keys()):
                    text_for_g = focal_only if (g in marker_readout_genes and g not in study_genes) else fx
                    text_for_g_upper = focal_only_upper if (g in marker_readout_genes and g not in study_genes) else fx_upper
                    if g.upper() in text_for_g_upper:
                        if MATCH(g, text_for_g, text_for_g_upper):
                            matched_genes.add(g)
                for alias, canonical in alias_map.items():
                    text_for_a = focal_only if (canonical in marker_readout_genes and canonical not in study_genes) else fx
                    text_for_a_upper = focal_only_upper if (canonical in marker_readout_genes and canonical not in study_genes) else fx_upper
                    if alias in text_for_a_upper:
                        if AMATCH(alias, text_for_a, text_for_a_upper):
                            ensure_gene(canonical)
                            matched_genes.add(canonical)
            # Only attribute to study_genes if the study perturbs a single target gene (no multi-gene broadcast!)
            if not matched_genes and len(study_genes) == 1:
                matched_genes = set(study_genes)

            direction = eff.get("effect_direction", "Neutral_or_NonInferior")
            penalty = bool(eff.get("secretome_or_differentiation_penalty_observed", False))
            lrr = eff.get("log_response_ratio_or_smd")
            if lrr is not None:
                lrr = _finite_real(lrr, f"{pmid}:{eff.get('effect_id')} LRR")
            orient_sign = _effect_orientation(eff)

            redacted_pert_matches = re.findall(r"\bpert=(GOF_OE|KO_LOF|BOTH|NONE)\b", fx) if "[REDACTED_NON_OA_EXCERPT:" in fx else []
            if redacted_pert_matches:
                p_tag_val = redacted_pert_matches[-1]
                has_explicit_gof_tokens = p_tag_val in ("GOF_OE", "BOTH")
                has_explicit_ko_tokens = p_tag_val in ("KO_LOF", "BOTH")
            else:
                has_explicit_gof_tokens = any(
                    tok in focal_only_upper or tok in fx_upper
                    for tok in (
                        "OVEREXPRESSION", "OVEREXPRESS", "KNOCK-IN", "KNOCKIN", "TRANSGENE",
                        "TRANSGENIC", "ECTOPIC", "LENTIVIRAL", "GAIN-OF-FUNCTION", "GAIN OF FUNCTION",
                        "FORCED EXPRESSION", "CONSTITUTIVE EXPRESSION", "INDUCIBLE EXPRESSION",
                    )
                )
                has_explicit_ko_tokens = any(
                    tok in focal_only_upper or tok in fx_upper
                    for tok in (
                        "KNOCKOUT", " KO ", " KO,", " KO.", "-/-", "−/−", "DEFICIENT", "DEFICIENCY",
                        "DELETION", "NULL", "ABLATION", "SILENCING", "SHRNA", "SIRNA", "LOSS-OF-FUNCTION",
                        "LOSS OF FUNCTION", "INHIBITION", "INHIBITOR", "BLOCKADE",
                    )
                )

            for g in matched_genes:
                st = target_stats[g]
                if pmid and not is_non_weighted_record and pmid not in st["supporting_pmids"]:
                    st["supporting_pmids"].append(pmid)

                arm_ets = " ".join(sorted(study_gene_edit_types.get(g, set()))).upper()
                is_gof_default = st.get("default_perturbation_direction") == "Gain_of_Function_or_KnockIn"
                if has_explicit_gof_tokens and not has_explicit_ko_tokens:
                    pert_dir = "GOF_OE"
                elif has_explicit_ko_tokens and not has_explicit_gof_tokens:
                    pert_dir = "KO_LOF"
                elif any(
                    t in arm_ets
                    for t in (
                        "OVEREXPRESSION", "KNOCKIN", "KNOCK-IN", "TRANSGENE", "GAIN",
                        "LENTIVIRAL_OR_TRANSPOSON_OE", "SAFEHARBOR_KI_AAVS1_GAPDH_CLYBL",
                        "CRISPRA_ACTIVATION", "TRANSIENT_MRNA_OR_SURFACE_TETHERING",
                    )
                ) and not any(t in arm_ets for t in ("KNOCKOUT", "CRISPR_CAS9_KO", "CRISPR_CAS12A_MULTIPLEX_KO", "DELETION", "LOSS")):
                    pert_dir = "GOF_OE"
                elif any(t in arm_ets for t in ("KNOCKOUT", "KO", "DELETION", "BASEEDIT", "PRIMEEDIT", "SHRNA", "LOSS")):
                    pert_dir = "KO_LOF"
                else:
                    pert_dir = "ASSOCIATIONAL_OR_UNTYPED"

                if pert_dir == "GOF_OE":
                    st["n_gof_oe_effects"] += 1
                elif pert_dir == "KO_LOF":
                    st["n_ko_lof_effects"] += 1
                else:
                    st["n_associational_or_untyped_effects"] += 1

                explicit_dir = (
                    "GOF_OE"
                    if (has_explicit_gof_tokens and not has_explicit_ko_tokens)
                    else ("KO_LOF" if (has_explicit_ko_tokens and not has_explicit_gof_tokens) else None)
                )
                eff_bucket = pert_dir
                if eff_bucket == "ASSOCIATIONAL_OR_UNTYPED":
                    st.setdefault("contextual_untyped_effects", []).append({
                        "pmid": pmid, "effect_id": eff.get("effect_id"), "direction": direction,
                        "reason": "No donor intervention type; not consumed as KO/GOF efficacy"})
                    continue

                # Compute oriented signed LRR (preserving exact sign-flip sensitivity)
                oriented_lrr = None
                if isinstance(lrr, (int, float)):
                    clipped_raw = max(-3.0, min(3.0, float(lrr)))
                    oriented_lrr = orient_sign * clipped_raw
                    accumulate(st, "signed_lrr_weighted_sum", eff_weight * oriented_lrr)
                    accumulate(st, "signed_lrr_weight_total", eff_weight)
                    if eff_bucket == "GOF_OE":
                        accumulate(st, "gof_signed_lrr_weighted_sum", eff_weight * oriented_lrr)
                        accumulate(st, "gof_signed_lrr_weight_total", eff_weight)
                    else:
                        accumulate(st, "ko_signed_lrr_weighted_sum", eff_weight * oriented_lrr)
                        accumulate(st, "ko_signed_lrr_weight_total", eff_weight)

                is_harmful = (
                    (direction == "Deleterious_or_Impaired")
                    or penalty
                    or (oriented_lrr is not None and oriented_lrr < 0)
                )

                if is_harmful:
                    del_w = direction_weights["deleterious"]
                    pen_w = direction_weights["secretome_or_differentiation_penalty_addon"]
                    st["n_deleterious_or_penalty_effects"] += 1
                    harm_weight = del_w * eff_weight + (pen_w * eff_weight if penalty else 0.0)
                    accumulate(st, "weighted_deleterious_penalty_score", harm_weight)
                    if eff_bucket == "GOF_OE":
                        accumulate(st, "gof_weighted_deleterious_score", harm_weight)
                    else:
                        accumulate(st, "ko_weighted_deleterious_score", harm_weight)
                    if pmid and not is_non_weighted_record and pmid not in st["counter_evidence_pmids"]:
                        st["counter_evidence_pmids"].append(pmid)
                elif direction == "Protective_or_Improved":
                    prot_w = direction_weights["protective"]
                    st["n_protective_effects"] += 1
                    prot_val = prot_w * eff_weight
                    accumulate(st, "weighted_protective_score", prot_val)
                    if eff_bucket == "GOF_OE":
                        accumulate(st, "gof_weighted_protective_score", prot_val)
                    else:
                        accumulate(st, "ko_weighted_protective_score", prot_val)
                else:
                    st["n_neutral_effects"] += 1

    # Finalize posterior means and integrate perturbation-direction-aware gnomAD + OpenTargets constraints
    prior_w = shrink_parameters["prior_weight"]
    emp_slope = shrink_parameters["per_empirical_study_slope"]
    shrunk_mult = direction_weights["shrunk_lrr_multiplier"]
    pli_slope = omics_weights["gnomad_pli_ko_slope"]
    loeuf_slope = omics_weights["gnomad_loeuf_ko_slope"]
    ot_slope = omics_weights["opentargets_t1d_slope"]

    for g, st in target_stats.items():
        w_tot = st["signed_lrr_weight_total"]
        n_emp = st.get("n_empirical_studies", st["n_studies_total"])
        raw_mean_lrr = (st["signed_lrr_weighted_sum"] / w_tot) if w_tot > 0 else 0.0
        def shrunk_mean(weight: float, signed_sum: float) -> float:
            if weight <= 0:
                return 0.0
            denominator = _finite_real(weight + prior_w + emp_slope * n_emp,
                                       f"{g} shrinkage denominator")
            return _finite_real(signed_sum / denominator, f"{g} shrunk mean")

        shrunk_lrr = shrunk_mean(w_tot, st["signed_lrr_weighted_sum"])
        ko_w_tot = st["ko_signed_lrr_weight_total"]
        gof_w_tot = st["gof_signed_lrr_weight_total"]
        ko_shrunk_lrr = shrunk_mean(ko_w_tot, st["ko_signed_lrr_weighted_sum"])
        gof_shrunk_lrr = shrunk_mean(gof_w_tot, st["gof_signed_lrr_weighted_sum"])

        pli = st["gnomad_pli"]
        loeuf = st["gnomad_loeuf_upper"]
        ot = (
            st["opentargets_t1d_score"]
            if st.get("opentargets_t1d_score") is not None
            else st.get("curated_literature_prior_t1d_score")
        )

        pli_penalty = pli_slope * float(pli) if isinstance(pli, (int, float)) else 0.0
        loeuf_bonus = loeuf_slope * min(2.0, float(loeuf)) if isinstance(loeuf, (int, float)) else 0.0
        ot_bonus = ot_slope * float(ot) if isinstance(ot, (int, float)) else 0.0
        omics_ko_constraint_score = round(pli_penalty + loeuf_bonus + ot_bonus, 4)
        omics_gof_constraint_score = round(ot_bonus, 4)

        ko_base_ev = (
            st["ko_weighted_protective_score"]
            - st["ko_weighted_deleterious_score"]
            + shrunk_mult * ko_shrunk_lrr
        )
        gof_base_ev = (
            st["gof_weighted_protective_score"]
            - st["gof_weighted_deleterious_score"]
            + shrunk_mult * gof_shrunk_lrr
        )
        is_gof_default = st.get("default_perturbation_direction") == "Gain_of_Function_or_KnockIn"
        net_ko_score = ko_base_ev + omics_ko_constraint_score
        net_gof_score = gof_base_ev + omics_gof_constraint_score
        net_preserve_score = (
            st["weighted_protective_score"]
            - st["weighted_deleterious_penalty_score"]
            + shrunk_mult * shrunk_lrr
            + omics_gof_constraint_score
        )
        net_default_score = net_gof_score if is_gof_default else net_ko_score

        st["posterior_mean_log_rr"] = round(raw_mean_lrr, 4)
        st["shrunk_posterior_log_rr"] = round(shrunk_lrr, 4)
        st["omics_ko_constraint_adjustment"] = omics_ko_constraint_score
        st["omics_gof_constraint_adjustment"] = omics_gof_constraint_score
        st["net_evidence_ko_score"] = round(net_ko_score, 4)
        st["net_evidence_gof_score"] = round(net_gof_score, 4)
        st["net_evidence_preserve_score"] = round(net_preserve_score, 4)
        st["net_evidence_and_constraint_score"] = round(net_default_score, 4)
        st["total_evidence_contrasts_in_ring1"] = (
            st["n_protective_effects"] + st["n_deleterious_or_penalty_effects"] + st["n_neutral_effects"]
        )
        def_g_se_p = dispersion_parameters["default_gene_se"]
        heur_se = round((0.40 * (def_g_se_p / 0.35)) / math.sqrt(max(0.25, w_tot + 0.25 * n_emp)), 4)
        st["posterior_se"] = heur_se
        st["heuristic_evidence_sparsity_se"] = heur_se
        st["uncertainty_field_note"] = (
            "Heuristic inverse-root evidence-weight dispersion index (NOT a calibrated Bayesian posterior standard error)."
        )
        st["supporting_pmids"] = sorted(st["supporting_pmids"])[:25]
        st["counter_evidence_pmids"] = sorted(st["counter_evidence_pmids"])[:25]
        for key, value in st.items():
            if isinstance(value, float):
                _finite_real(value, f"{g}.{key} derived profile")

    return target_stats


_GRAVEYARD_CACHE: dict = {"mtime": None, "entries": None}


EXPECTED_GRAVEYARD_IDS = frozenset(['GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER', 'GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH', 'GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK', 'GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF', 'GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT', 'GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS', 'GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2', 'GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION', 'GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT', 'GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT'])


def _load_graveyard_entries_by_id() -> dict[str, dict]:
    """Load the complete implemented rule set; malformed policy fails closed."""
    try:
        stat = GRAVEYARD_PATH.stat()
        identity = (str(GRAVEYARD_PATH.resolve()), stat.st_dev, stat.st_ino, stat.st_ctime_ns, stat.st_mtime_ns, stat.st_size)
        if _GRAVEYARD_CACHE.get("identity") != identity or _GRAVEYARD_CACHE["entries"] is None:
            raw = json.loads(GRAVEYARD_PATH.read_text(encoding="utf-8"))
            rows = raw.get("graveyard_entries") if isinstance(raw, dict) else None
            if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
                raise ValueError("Graveyard must contain a rule-object list")
            identifiers = [row.get("graveyard_id") for row in rows]
            if any(not isinstance(identifier, str) for identifier in identifiers) or len(set(identifiers)) != len(identifiers) or set(identifiers) != EXPECTED_GRAVEYARD_IDS:
                raise ValueError("Graveyard rules must exactly cover all implemented unique identifiers")
            known_modules = set().union(ALLO_CLOAK_MODULES, AUTOANTIGEN_MODULES, CYTOKINE_GRANZYME_MODULES, ER_AMYLOID_MODULES, NICHE_VASCULAR_MODULES, SAFETY_SWITCH_MODULES)
            expected_logit_deltas = {
                "GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER": None,
                "GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH": None,
                "GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK": 0.0,
                "GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF": 0.0,
                "GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT": 0.0,
                "GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS": 0.0,
                "GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2": None,
                "GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION": None,
                "GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT": -0.20,
                "GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT": -0.25,
            }
            expected_severities = {
                "GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER": "BLOCK",
                "GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH": "BLOCK",
                "GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK": "WARN",
                "GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF": "WARN",
                "GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT": "WARN",
                "GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS": "WARN",
                "GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2": "BLOCK",
                "GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION": "BLOCK",
                "GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT": "WARN",
                "GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT": "WARN",
            }
            expected_trigger_modules = {
                "GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER": {"CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1"},
                "GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH": {"SAF_G6PC2_ENDOGENOUS_ICASP9"},
                "GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK": {"AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR", "SAF_G6PC2_ENDOGENOUS_ICASP9"},
                "GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF": {"ER_RNLS_KO", "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE"},
                "GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT": {"ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO"},
                "GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS": {
                    "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR",
                    "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE",
                    "AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT",
                    "ER_IAPP_PRAMLINTIDE_RNLS_WT",
                },
                "GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2": {"AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR"},
                "GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION": {"CYTO_CONSTITUTIVE_IL10_VIRAL_OE"},
                "GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT": {"ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO"},
                "GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT": {"ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", "ALLO_B2M_CIITA_KO_CD47_LENTI"},
            }
            expected_target_gates = {
                "GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER": "Gate2B_and_Gate4_hard_BLOCK",
                "GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH": "Gate5_hard_BLOCK",
                "GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK": "Gate1_precautionary_WARN",
                "GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF": "Gate1_and_Gate5_precautionary_WARN",
                "GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT": "Gate3_precautionary_WARN",
                "GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS": "manufacturing_base_edit_gamma_B_WARN",
                "GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2": "Gate1_and_Gate4_hard_BLOCK",
                "GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION": "Gate2B_and_Gate4_hard_BLOCK",
                "GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT": "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield",
                "GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT": "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth",
            }
            valid_gate_prefixes = ("Gate1", "Gate2A", "Gate2B", "Gate3", "Gate4", "Gate5", "Gate6", "manufacturing_base_edit")
            for row in rows:
                sev_lvl = row.get("severity_level")
                gid = row["graveyard_id"]
                if sev_lvl not in {"BLOCK", "PENALIZE", "WARN", "INFO"} or sev_lvl != expected_severities.get(gid):
                    raise ValueError(f"Invalid or mutated Graveyard severity_level {sev_lvl!r} for {gid}")
                if "evaluator_logit_delta" not in row or "evaluator_target_gate" not in row:
                    raise ValueError("Graveyard rule missing required evaluator_logit_delta or evaluator_target_gate field")
                target_gate = row.get("evaluator_target_gate")
                if (
                    not isinstance(target_gate, str)
                    or not target_gate.strip()
                    or not target_gate.startswith(valid_gate_prefixes)
                    or not re.match(r"^(Gate(1|2A|2B|3|4|5|6)|manufacturing_base_edit)(_|$)", target_gate)
                    or target_gate != expected_target_gates.get(gid)
                ):
                    raise ValueError(f"Graveyard rule {gid} evaluator_target_gate {target_gate!r} diverges from canonical target gate {expected_target_gates.get(gid)!r}")
                logit_delta = row.get("evaluator_logit_delta")
                if sev_lvl == "BLOCK":
                    if logit_delta is not None or not target_gate.endswith("_hard_BLOCK"):
                        raise ValueError(f"Graveyard BLOCK rule {gid} must have evaluator_logit_delta=None and *_hard_BLOCK target gate")
                else:
                    if isinstance(logit_delta, bool) or not isinstance(logit_delta, (int, float)) or not math.isfinite(float(logit_delta)) or float(logit_delta) > 0.0:
                        raise ValueError(f"Graveyard {sev_lvl} rule {gid} requires a finite non-positive evaluator_logit_delta")
                    if gid in expected_logit_deltas and not math.isclose(float(logit_delta), float(expected_logit_deltas[gid]), rel_tol=0.0, abs_tol=1e-12):
                        raise ValueError(f"Graveyard rule {gid} evaluator_logit_delta={logit_delta} diverges from canonical policy {expected_logit_deltas[gid]}")
                pattern = row.get("trigger_pattern")
                modules = pattern.get("requires_any_module") if isinstance(pattern, dict) else None
                if not isinstance(modules, list) or not modules or any(not isinstance(module, str) or module not in known_modules for module in modules):
                    raise ValueError("Graveyard requires a nonempty known-module trigger list")
                if set(modules) != expected_trigger_modules.get(gid):
                    raise ValueError(f"Graveyard rule {gid} trigger modules {modules!r} diverge from canonical trigger set")
            _GRAVEYARD_CACHE.update(identity=identity, mtime=stat.st_mtime_ns, entries={row["graveyard_id"]: row for row in rows})
        return copy.deepcopy(_GRAVEYARD_CACHE["entries"])
    except (OSError, TypeError, KeyError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot load Graveyard policy: {exc}") from exc


def evaluate_graveyard_rules_for_stack(
    allo_k: str, auto_k: str, cyto_k: str, er_k: str, niche_k: str, saf_k: str, total_b: int
) -> list[dict]:
    """
    Executable Idea Graveyard Rule Engine (Audit §6.4 + Audit A08/A09 + Grok P2-2).
    Dynamically reads severity_level, curation_status, and trigger_pattern.requires_any_module from
    data/graveyard/idea_graveyard_v0.5.json so JSON rule updates immediately drive stack blocking and warnings.
    """
    selected = [(allo_k, ALLO_CLOAK_MODULES), (auto_k, AUTOANTIGEN_MODULES),
                (cyto_k, CYTOKINE_GRANZYME_MODULES), (er_k, ER_AMYLOID_MODULES),
                (niche_k, NICHE_VASCULAR_MODULES), (saf_k, SAFETY_SWITCH_MODULES)]
    if any(not isinstance(key, str) or key not in choices for key, choices in selected):
        raise ValueError("Unknown or mistyped stack module")
    minimum_b = sum(choices[key]["base_edits"] for key, choices in selected)
    if type(total_b) is not int or total_b < minimum_b:
        raise ValueError("Base-edit count cannot understate selected module edits")
    gy_map = _load_graveyard_entries_by_id()
    active_modules = {allo_k, auto_k, cyto_k, er_k, niche_k, saf_k}

    def _matches_gy(gid: str, default_modules: set[str], extra_cond: bool = True) -> bool:
        if not extra_cond:
            return False
        entry = gy_map.get(gid)
        req_mods = set((entry.get("trigger_pattern") or {}).get("requires_any_module") or default_modules) if entry else default_modules
        return bool(active_modules & req_mods)

    def _sev(gid: str, default_sev: str) -> str:
        entry = gy_map.get(gid) or {}
        val = str(entry.get("severity_level", default_sev)).strip().upper()
        if val not in {"BLOCK", "PENALIZE", "WARN", "INFO"}:
            raise ValueError(f"Invalid severity_level {val!r} for {gid}")
        return val

    triggered = []
    if _matches_gy("GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER", {"CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER",
            "severity_level": _sev("GRAVEYARD_001_STAT1_JAK1_KO_WITH_GAS_IND_PROMOTER", "BLOCK"),
            "reason": "Clonal STAT1-/- blocks GAS/ISRE promoter induction of PD-L1/SERPINB9 (Sintov 2022 PMID:36055241; Parent 2021 PMID:34407395).",
        })
    if _matches_gy("GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH", {"SAF_G6PC2_ENDOGENOUS_ICASP9"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH",
            "severity_level": _sev("GRAVEYARD_002_LINEAGE_RESTRICTED_G6PC2_OR_INS_ICASP9_SUICIDE_SWITCH", "BLOCK"),
            "reason": "Lineage-restricted (G6PC2::iCasp9) suicide switch fails to ablate dedifferentiated, teratoma, or SLC18A1+ off-target cells (PMID:27626039, PMID:35795779).",
        })
    if _matches_gy("GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK", {"AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR", "SAF_G6PC2_ENDOGENOUS_ICASP9"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK",
            "severity_level": _sev("GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK", "WARN"),
            "evidence_tier": str((gy_map.get("GRAVEYARD_003_G6PC2_NULL_FASTING_HYPOGLYCEMIA_RISK") or {}).get("curation_status", "precautionary_mechanistic_hazard")),
            "reason": "Precautionary mechanistic hazard: Biallelic G6PC2-/- shifts the beta-cell glucose dose-response curve leftward and lowers fasting blood glucose (Bosma 2020 PMID:32569842; Overway 2022 PMID:34954144).",
        })
    if _matches_gy("GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF", {"ER_RNLS_KO", "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF",
            "severity_level": _sev("GRAVEYARD_004_RNLS_KO_HUMAN_IPSC_BETA_GSIS_AND_OCR_TRADEOFF", "WARN"),
            "reason": "Human RNLS-/- iPSC-beta aggregates show elevated basal OCR while poor dynamic GSIS persisted in both WT and RNLS KO in MacDonald et al. 2025 Fig. 6 (PMID:39988068 vs mouse NIT-1 screen Cai et al. 2020 PMID:32719542).",
        })
    if _matches_gy("GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT", {"ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT",
            "severity_level": _sev("GRAVEYARD_005_HLA_C_MINUS21T_LEADER_PEPTIDE_NKG2A_BLINDSPOT", "WARN"),
            "reason": "Requires HLA-C haplotype selection (Xu 2019 PMID:30853558) and empirical verification that a single retained HLA-C allele supplies sufficient leader peptide to stabilize surface HLA-E after HLA-A/HLA-B deletion (Ramsuran 2018 PMID:29302013; Parent 2021 PMID:34407395; Lin 2023 PMID:37264229).",
        })
    if _matches_gy(
        "GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS",
        {
            "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR",
            "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE",
            "AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT",
            "ER_IAPP_PRAMLINTIDE_RNLS_WT",
        },
        extra_cond=(total_b >= 2),
    ):
        triggered.append({
            "graveyard_id": "GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS",
            "severity_level": _sev("GRAVEYARD_006_MULTIPLEX_BASE_PRIME_EDIT_NICKASE_TRANSLOCATIONS", "WARN"),
            "reason": f"Multiplex base/prime editing (B={total_b}) carries precautionary nickase single-strand nick / DSB and deaminase bystander structural-variant risk (extrapolated from Wu 2024 PMID:39529170).",
        })
    if _matches_gy("GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2", {"AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2",
            "severity_level": _sev("GRAVEYARD_007_HIGH_PLI_KNOCKOUTS_PTPN2_TNFAIP3_GAD2", "BLOCK"),
            "reason": "Gross knockout or coding/splice disruption of high-pLI beta-cell brake or paracrine gene (PTPN2-/- [K=3 KO], TNFAIP3-/-, or multiplex GAD2 coding disruption [B=4]) is prohibited (PMID:19336676, PMID:10523611, PMID:31581152, PMID:2550826, PMID:32432213).",
        })
    if _matches_gy("GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION", {"CYTO_CONSTITUTIVE_IL10_VIRAL_OE"}):
        triggered.append({
            "graveyard_id": "GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION",
            "severity_level": _sev("GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION", "BLOCK"),
            "evidence_tier": str((gy_map.get("GRAVEYARD_008_CONSTITUTIVE_IL10_VIRAL_OVEREXPRESSION") or {}).get("curation_status", "conservative_project_exclusion_policy_under_mechanistic_uncertainty")),
            "reason": "Conservative project exclusion policy under mechanistic uncertainty: unbuffered constitutive viral IL-10 overexpression alone is excluded in favor of contact-dependent checkpoints (PD-L1) and intracellular SERPINB9 (Polukort 2016 PMID:27183617 is retained strictly as non-islet mucosal mast-cell context).",
        })
    if _matches_gy("GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT", {"ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO"}):
        gy009 = gy_map.get("GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT") or {}
        triggered.append({
            "graveyard_id": "GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT",
            "severity_level": _sev("GRAVEYARD_009_RESTING_B2M_PROMOTER_HOUR0_INNATE_DENSITY_DEFICIT", "WARN"),
            "evaluator_logit_delta": float(gy009.get("evaluator_logit_delta", -0.20)),
            "evaluator_target_gate": str(gy009.get("evaluator_target_gate", "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield")),
            "evaluator_disposition_rule": str(gy009.get("evaluator_disposition_rule", "Applies -0.20 logit discount in Gate 3 + requires Decision Experiment #1 QuantiBRITE gate")),
            "reason": "Requires QuantiBRITE Hour-0 verification that resting endogenous B2M promoter achieves inhibitory HLA-E and CD47 surface density on Stage 7 beta-cells and iECs (Han 2019 PMID:31040209; Hu 2023 PMID:37156915; Xu 2024 PMID:38508139; Lin 2023 PMID:37264229).",
        })
    if _matches_gy("GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT", {"ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", "ALLO_B2M_CIITA_KO_CD47_LENTI"}):
        gy010 = gy_map.get("GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT") or {}
        triggered.append({
            "graveyard_id": "GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT",
            "severity_level": _sev("GRAVEYARD_010_MHC_NULL_INDIRECT_APC_BYSTANDER_BLINDSPOT", "WARN"),
            "evaluator_logit_delta": float(gy010.get("evaluator_logit_delta", -0.25)),
            "evaluator_target_gate": str(gy010.get("evaluator_target_gate", "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth")),
            "evaluator_disposition_rule": str(gy010.get("evaluator_disposition_rule", "Applies -0.25 logit discount in Gate 2B + requires Decision Experiment #2 indirect APC / SA-FasL gate")),
            "reason": "Donor MHC loss does not block host APC presentation of shed autoantigens and alloantigens. An INS leader/DRiP edit does not remove other epitopes or establish indirect immune protection; the warning remains until exact-product assays qualify the residual hazard (PMID:28263308, PMID:36055241, PMID:35559682).",
        })
    return triggered


def signed_gene_signal(
    target_profiles: dict,
    genes: list[str] | tuple[str, ...],
    scale: float = 0.25,
    is_knockout: bool = False,
    branch_key: str | None = None,
    mode: str | None = None,
) -> float:
    """
    Computes a strictly sign-sensitive, perturbation-direction-aware evidence signal from target_profiles.
    - Dynamically reads `signed_gene_signal_slope`, `signed_gene_signal_clip_bounds`, and (when `branch_key` is
      provided) `gene_signal_scales_by_gate_and_branch` from `COMPLETE_GATE_COEFFICIENT_REGISTRY['auxiliary_evaluator_hyperparameters']` (Audit R10).
    - When `is_knockout=True` (or `mode='ko'`), uses the loss-of-function constraint score (`net_evidence_ko_score`,
      which applies the gnomAD pLI loss-of-function intolerance penalty).
    - When `mode='gof'` (or for gain-of-function transgenes), uses `net_evidence_gof_score` / `net_evidence_and_constraint_score`.
    - When `mode in ('preserve', 'base_edit')`, uses `net_evidence_preserve_score` (or `net_evidence_and_constraint_score`)
      so gnomAD pLI loss-of-function intolerance is never applied to preserved or base-edited genes (Finding 086).
    """
    if not isinstance(target_profiles, dict) or not isinstance(genes, (list, tuple)) or any(not isinstance(g, str) for g in genes):
        raise TypeError("Gene signal requires a target mapping and gene identifiers")
    if type(is_knockout) is not bool or mode not in (None, "ko", "gof", "preserve", "base_edit"):
        raise ValueError("Unknown gene-signal perturbation mode")
    _finite_real(scale, "Gene-signal scale")
    if not target_profiles or not genes:
        return 0.0
    eff_scale = _finite_real(_branch_scale(branch_key, scale) if branch_key is not None else scale, "Branch gene-signal scale")
    slope = _finite_real(_aux_param("signed_gene_signal_slope", 0.42), "Gene-signal slope")
    clip_bounds = _aux_param("signed_gene_signal_clip_bounds", [-4.5, 4.5])
    if not isinstance(clip_bounds, (list, tuple)) or len(clip_bounds) != 2:
        raise ValueError("Gene-signal clipping requires exactly two finite bounds")
    lo_b = _finite_real(clip_bounds[0], "Lower gene-signal bound")
    hi_b = _finite_real(clip_bounds[1], "Upper gene-signal bound")
    if not lo_b <= 0 <= hi_b or lo_b >= hi_b:
        raise ValueError("Gene-signal clipping bounds must be ordered and bracket zero")
    total = 0.0
    for g in genes:
        st = target_profiles.get(g)
        if st is None:
            continue
        if not isinstance(st, dict):
            raise ValueError("Target profiles must contain explicit gene mappings")
        if not st:
            continue
        if (is_knockout or mode == "ko") and "net_evidence_ko_score" in st:
            net_s = _finite_real(st["net_evidence_ko_score"], f"{g} KO evidence score")
        elif mode == "gof" and "net_evidence_gof_score" in st:
            net_s = _finite_real(st["net_evidence_gof_score"], f"{g} GOF evidence score")
        elif mode in ("preserve", "base_edit") and "net_evidence_preserve_score" in st:
            net_s = _finite_real(st["net_evidence_preserve_score"], f"{g} preserve evidence score")
        else:
            net_s = _finite_real(st.get("net_evidence_and_constraint_score", 0.0), f"{g} evidence score")
        term = _finite_real(eff_scale * (slope * net_s), f"{g} scaled gene signal")
        total = _finite_real(total + term, "Accumulated gene signal")
    bound_scale = (hi_b / 4.5) if total >= 0 else (abs(lo_b) / 4.5)
    # Compare before multiplying; overflow must never become a positive score.
    if bound_scale == 0:
        return 0.0
    if total >= 0 and total >= hi_b / bound_scale:
        return hi_b
    if total < 0 and total <= lo_b / bound_scale:
        return lo_b
    return max(lo_b, min(hi_b, total * bound_scale))


# Immutable per-evaluation snapshot. Hot-loop lookups never touch the filesystem.
_WETLAB_PARAMETER_DELTAS = MappingProxyType({})
_WETLAB_CALIBRATION_STATUS = {"artifact_present": False, "applied_to_production_priors": False,
                               "status_reason": "No evaluation context prepared"}


def prepare_wetlab_evaluation_context(report: dict | None = None, path: Path | None = None) -> dict:
    """Validate and freeze one posterior at the evaluation boundary.

    The endpoint-to-logit scales are declared engineering mappings, not learned
    causal calibrations. Unmapped observations remain visible and unapplied.
    """
    global _WETLAB_PARAMETER_DELTAS, _WETLAB_CALIBRATION_STATUS
    _WETLAB_PARAMETER_DELTAS = MappingProxyType({})
    status = {"calibration_artifact_path": "data/extracted/wetlab_bayesian_posterior_calibration.json",
              "artifact_present": False, "applied_to_production_priors": False,
              "raw_endpoint_posterior_deltas": {}, "mapped_nonzero_parameter_deltas": {},
              "unmapped_endpoint_deltas": {}, "retained_unobserved_endpoints": [], "consumed_parameter_deltas": {}, "mapping_estimand": "uncalibrated engineering logit shift"}
    try:
        if report is None:
            source = path or PROJECT_ROOT / status["calibration_artifact_path"]
            if not source.exists():
                status["status_reason"] = "Not_Found"
                _WETLAB_CALIBRATION_STATUS = status
                return dict(status)
            report = json.loads(source.read_text(encoding="utf-8"))
        if not isinstance(report, dict):
            raise ValueError("Posterior report must be a dictionary")
        validator_spec = importlib.util.spec_from_file_location("openislet_wetlab_validator", PROJECT_ROOT / "scripts/ingest_wetlab_results.py")
        validator = importlib.util.module_from_spec(validator_spec)
        validator_spec.loader.exec_module(validator)
        status.update(artifact_present=True, calibration_mode=report.get("calibration_mode"),
                      epistemic_lineage_status=report.get("epistemic_lineage_status"))
        if (report.get("calibration_mode") != "EMPIRICAL_WETLAB_INGESTION"
                or report.get("epistemic_lineage_status") not in validator.SUPPORTED_EMPIRICAL_LINEAGES
                or report.get("contains_synthetic_records", False) is not False):
            status["status_reason"] = "Rejected synthetic, mixed or unsupported lineage"
            _WETLAB_CALIBRATION_STATUS = status
            return dict(status)
        raw_updates = report.get("conjugate_posterior_updates")
        digest = hashlib.sha256(json.dumps(raw_updates, sort_keys=True, allow_nan=False).encode()).hexdigest()
        if report.get("conjugate_posterior_updates_sha256") != digest:
            raise ValueError("Missing or mismatched conjugate posterior SHA-256")
        # Use fresh identifiers because the validator also guards batch replay.
        validator._validate_prior_posterior_report(report, False, "RANKER-EVALUATION-" + digest, "evaluation:" + digest)
        updates = dict(raw_updates)
        updates.update(report.get("optional_pouch_and_36run_screen_updates_v80_070") or {})
        endpoint_map = {
            "stage7_sc_beta_yield_fraction_arm1b": ("stage7_trihormonal_120um_mvf_cao2_niche_prior", 1.0),
            "cd47_igv_hour0_fold_over_baseline_arm1b": ("resting_b2m_promoter_hour0_density_discount", 0.05),
            "serpinb9_p2a_gzmb_log10_kass_M_inv_s_inv": ("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience::CYTO_CXCL10_SWAP_PDL1_SERPINB9", 0.25),
            "aa_vs1_ucoe_icasp9_ap1903_24h_log10_kill": ("Gate5_72h_Vascular_Oxygen_and_80yr_Safety::SAF_AAVS1_CAG_ICASP9_BIALLELIC", 0.25),
            "pouch_min_slab_po2_14d_mmhg": ("stage7_trihormonal_120um_mvf_cao2_niche_prior", 0.01),
            "ifng_48h_hla_dr_dq_negative_pct_arm1b": ("Gate2A_Allo_Adaptive_T_and_DSA_Stealth::ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", 0.01),
        }
        genesis_prior_means = {
            row[0]: float(row[1])
            for row in (*validator.ENDPOINT_SPECS, *validator.EXTENDED_STAGE0_STAGE1_SPECS)
        }
        genesis_prior_sds = {
            row[0]: float(row[2])
            for row in (*validator.ENDPOINT_SPECS, *validator.EXTENDED_STAGE0_STAGE1_SPECS)
        }
        is_multibatch = len(report.get("ingested_batch_ids", [])) > 1
        mapped = {}
        for endpoint, update in updates.items():
            if endpoint not in genesis_prior_means or endpoint not in genesis_prior_sds:
                raise ValueError(f"Unrecognized calibration endpoint {endpoint!r}")
            genesis_m = genesis_prior_means[endpoint]
            genesis_sd = genesis_prior_sds[endpoint]
            if "genesis_prior_mean_exact" in update:
                gm_raw = validator._require_finite_float(update["genesis_prior_mean_exact"], f"{endpoint}.genesis_prior_mean_exact")
                if not math.isclose(gm_raw, genesis_m, rel_tol=1e-12, abs_tol=0.0):
                    raise ValueError(f"Endpoint {endpoint!r} genesis_prior_mean_exact={gm_raw} does not match canonical prior {genesis_m}")
            if "genesis_prior_sd_exact" in update:
                gsd_raw = validator._require_finite_float(update["genesis_prior_sd_exact"], f"{endpoint}.genesis_prior_sd_exact")
                if not math.isclose(gsd_raw, genesis_sd, rel_tol=1e-12, abs_tol=0.0):
                    raise ValueError(f"Endpoint {endpoint!r} genesis_prior_sd_exact={gsd_raw} does not match canonical prior sd {genesis_sd}")
            delta = float(update.get("posterior_mean_exact", update["posterior_mean"])) - genesis_m
            status["raw_endpoint_posterior_deltas"][endpoint] = delta
            if endpoint in endpoint_map:
                prior_sd_cur = float(update.get("prior_sd_exact", update.get("prior_sd", genesis_sd)))
                has_historical_update = is_multibatch and (
                    prior_sd_cur < genesis_sd * (1.0 - 1e-12)
                )
                if update.get("status") == "UPDATED_FROM_OBSERVED_READOUT" or has_historical_update:
                    key, scale = endpoint_map[endpoint]
                    mapped[key] = mapped.get(key, 0.0) + scale * delta
                else:
                    status["retained_unobserved_endpoints"].append(endpoint)
            else:
                status["unmapped_endpoint_deltas"][endpoint] = delta
        mapped = {key: value for key, value in mapped.items() if value != 0.0}
        _WETLAB_PARAMETER_DELTAS = MappingProxyType(mapped)
        status["mapped_nonzero_parameter_deltas"] = dict(mapped)
        status["consumed_parameter_deltas"] = dict(mapped)
        status["active_posterior_logit_deltas"] = dict(mapped)
        status["arm1b_scoped_gate1_only_parameter_deltas"] = {
            k: v for k, v in mapped.items() if k != "stage7_trihormonal_120um_mvf_cao2_niche_prior"
        }
        status["endpoint_scope_qualification_note_astra_r1"] = (
            "stage7_sc_beta_yield_fraction_arm1b (scale 1.0) and pouch_min_slab_po2_14d_mmhg (scale 0.01) both accumulate "
            "into stage7_trihormonal_120um_mvf_cao2_niche_prior in the shared engineering calibration map; arm1b_scoped_gate1_only_parameter_deltas "
            "isolates construct-specific Gate 2A/3/4/5 updates from shared niche priors (AUDIT-REF)."
        )
        status["applied_to_production_priors"] = bool(mapped)
        status["validated_report_sha256"] = hashlib.sha256(json.dumps(report, sort_keys=True, allow_nan=False).encode()).hexdigest()
        status["status_reason"] = "Validated empirical mapped shifts frozen for evaluation" if mapped else "Validated observations; no nonzero consumed parameter shift"
    except (ValueError, TypeError, KeyError, OSError) as exc:
        status["status_reason"] = f"Rejected invalid posterior: {exc}"
    _WETLAB_CALIBRATION_STATUS = status
    return dict(status)


def _load_active_wetlab_posterior_deltas():
    """Compatibility accessor; returns the frozen snapshot, with zero I/O."""
    return _WETLAB_PARAMETER_DELTAS


def _cd47_gate1_gsis_penalty_signal(allo_k: str, target_profiles: dict) -> float:
    """
    Routes CD47's deleterious GSIS counter-evidence (Ghimire et al. 2023 PMID:37820008: full-length CD47
    receptor signaling inhibits beta-cell GSIS via Cdc42; mitigated by signaling-incompetent CD47-IgV in
    Xu et al. 2024 PMID:38508139) into Gate1_Maturation_GSIS_and_Paracrine for CD47-expressing allo-cloak
    modules (Grok P1-3; live-wired to cd47_gsis_penalty_scale per Audit R10).
    """
    if allo_k not in ("ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", "ALLO_B2M_CIITA_KO_CD47_LENTI"):
        return 0.0
    st = (target_profiles or {}).get("CD47") or {}
    del_pen = _finite_real(st.get("weighted_deleterious_penalty_score", 0.0), "CD47 deleterious evidence penalty")
    if del_pen < 0:
        raise ValueError("CD47 deleterious evidence penalty must be nonnegative")
    mult = 1.0 if allo_k == "ALLO_B2M_CIITA_KO_CD47_LENTI" else 0.5
    pen_scale = _finite_real(_aux_param("cd47_gsis_penalty_scale", -0.25), "CD47 GSIS penalty scale")
    return pen_scale * mult * del_pen


def sigmoid(x: float) -> float:
    x = _finite_real(x, "Gate logit")
    return 1.0 / (1.0 + math.exp(-max(-15.0, min(15.0, x))))


_DEFAULT_CLAMPED_SIGMOID = sigmoid


def sigmoid_unclipped(x: float) -> float:
    """Overflow-safe logistic evaluation without [-15, 15] saturation clamping or [-745, -700) truncation (AUDIT-REF L1-01)."""
    x = _finite_real(x, "Gate logit")
    if x >= 0.0:
        if x > 700.0:
            return 1.0
        return 1.0 / (1.0 + math.exp(-x))
    if x >= -700.0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def compute_series_system_utility(
    allo_k: str,
    auto_k: str,
    cyto_k: str,
    er_k: str,
    niche_k: str,
    saf_k: str,
    target_profiles: dict,
    recipient_context: str = "Pediatric_Early_Onset_T1D",
) -> dict:
    """
    Computes a Prior-Anchored Multi-Gate Scenario Index in [0, 1].
    - Returns 0.0 when target_profiles is empty or contains zero extracted empirical studies (AT-20).
    - Dynamically wires `ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE` and `COMPLETE_GATE_COEFFICIENT_REGISTRY`
      (including all 8 `auxiliary_evaluator_hyperparameters`) into gate logits and manufacturing
      penalties (`locus_dsb_penalty_rate_lambda_K`, `base_edit_nickase_penalty_factor_gamma_B`,
      `indirect_apc_presentation_residual_hazard_cap`, `resting_b2m_promoter_hour0_density_discount`,
      `stage7_trihormonal_120um_mvf_cao2_niche_prior`, `tyk2_stage0_differentiation_yield_penalty`).
    - Propagates gene-level posterior standard errors (`posterior_se`) into `lower_10pct_cred_bound` (AT-18).
    """
    if recipient_context not in ("Pediatric_Early_Onset_T1D", "Adult_Onset_LADA"):
        raise ValueError("Unknown recipient context; select an explicit modeled cohort")
    if not isinstance(target_profiles, dict):
        raise TypeError("Target profiles must be a dictionary")
    # Validate all counts before reducing; a generator would stop at the first
    # positive count and leave malformed later profiles unchecked.
    has_studies = bool(target_profiles) and any([
        _profile_has_empirical_studies(st) for st in target_profiles.values()
    ])
    if not has_studies:
        zero_gates = {
            "Gate1_Maturation_GSIS_and_Paracrine": 0.0,
            "Gate2A_Allo_Adaptive_T_and_DSA_Stealth": 0.0,
            "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth": 0.0,
            "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield": 0.0,
            "Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience": 0.0,
            "Gate5_72h_Vascular_Oxygen_and_80yr_Safety": 0.0,
        }
        return {
            "gate_scores": zero_gates,
            "gate_logits_unrounded": zero_gates,
            "unpenalized_gate_product": 0.0,
            "geometric_mean_gate_score": 0.0,
            "geometric_mean_six_subgate_score": 0.0,
            "limiting_bottleneck_gate": "Gate1_Maturation_GSIS_and_Paracrine",
            "limiting_bottleneck_score": 0.0,
            "manufacturing_locus_penalty_factor": 0.0,
            "series_system_joint_survival_probability": 0.0,
            "series_system_joint_survival_probability_unrounded": 0.0,
            "lower_10pct_cred_bound": 0.0,
            "estimand_note": "Zero evidence loaded -> 0.0 relative index.",
        }

    total_cuts = (
        ALLO_CLOAK_MODULES[allo_k]["cuts"]
        + AUTOANTIGEN_MODULES[auto_k]["cuts"]
        + CYTOKINE_GRANZYME_MODULES[cyto_k]["cuts"]
        + ER_AMYLOID_MODULES[er_k]["cuts"]
        + NICHE_VASCULAR_MODULES[niche_k]["cuts"]
        + SAFETY_SWITCH_MODULES[saf_k]["cuts"]
    )
    total_base_edits = (
        ALLO_CLOAK_MODULES[allo_k]["base_edits"]
        + AUTOANTIGEN_MODULES[auto_k]["base_edits"]
        + CYTOKINE_GRANZYME_MODULES[cyto_k]["base_edits"]
        + ER_AMYLOID_MODULES[er_k]["base_edits"]
        + NICHE_VASCULAR_MODULES[niche_k]["base_edits"]
        + SAFETY_SWITCH_MODULES[saf_k]["base_edits"]
    )

    _pw_ent = COMPLETE_GATE_COEFFICIENT_REGISTRY.get("pw_base_scale", 0.08)
    pw = _finite_real(_pw_ent.get("logit_value", 0.08) if isinstance(_pw_ent, dict) else _pw_ent, "Structural logit scale")

    # Gate 1: Lineage Purity, Dynamic GSIS & Low-Glucose Shut-Off
    logit_g1 = (
        pw * 5.0
        + _rdelta("Gate1_Maturation_GSIS_and_Paracrine", "base_intercept", 0.40)
        + _cd47_gate1_gsis_penalty_signal(allo_k, target_profiles)
    )
    if niche_k == "NICHE_TRIHORMONAL_120UM_MVF_CAO2":
        logit_g1 += (
            pw * 5.625
            + _pdelta("stage7_trihormonal_120um_mvf_cao2_niche_prior", 0.45, "Gate1_Maturation_GSIS_and_Paracrine", "NICHE_TRIHORMONAL_120UM_MVF_CAO2", 0.45)
            + _rdelta("Gate1_Maturation_GSIS_and_Paracrine", "NICHE_TRIHORMONAL_120UM_MVF_CAO2", 0.45)
            + signed_gene_signal(target_profiles, ("INS", "GAD2"), scale=0.20, mode="preserve", branch_key="Gate1_NICHE_INS_GAD2")
        )
    elif niche_k == "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION":
        logit_g1 -= pw * 5.625 - _rdelta("Gate1_Maturation_GSIS_and_Paracrine", "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION", -0.45)

    if allo_k == "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO":
        # Clonal TYK2-/- carried through Stage 0->7 hPSC differentiation compromises endocrine progenitor yield by ~40-60% (Chandra et al. 2022 PMID:36289205)
        logit_g1 += (
            signed_gene_signal(target_profiles, ["TYK2"], scale=0.20, is_knockout=True, branch_key="Gate1_ALLO_HLA_C_TYK2_KO")
            - pw * 4.0
            + _pdelta("tyk2_stage0_differentiation_yield_penalty", -0.32, "Gate1_Maturation_GSIS_and_Paracrine", "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO", -0.32)
            + _rdelta("Gate1_Maturation_GSIS_and_Paracrine", "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO", -0.32)
        )

    if auto_k == "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR":
        logit_g1 += (
            signed_gene_signal(target_profiles, ["G6PC2", "SLC30A8", "PTPN2"], scale=0.35, is_knockout=True, branch_key="Gate1_AUTO_FULL_KO")
            + signed_gene_signal(target_profiles, ["GAD2"], scale=0.35, mode="base_edit", branch_key="Gate1_AUTO_FULL_KO")
            - pw * 6.875
            + _rdelta("Gate1_Maturation_GSIS_and_Paracrine", "AUTO_FULL_KO_G6PC2_SLC30A8_PLUS_BASE_INS_GAD2_IAPP_CXADR", -0.55)
        )

    if er_k in ("ER_RNLS_KO", "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE"):
        logit_g1 += (
            signed_gene_signal(target_profiles, ["RNLS"], scale=0.25, is_knockout=True, branch_key="Gate1_ER_RNLS_KO")
            - pw * 4.75
            + _rdelta("Gate1_Maturation_GSIS_and_Paracrine", "ER_RNLS_KO_or_PRAMLINTIDE", -0.38)
        )

    # Gate 2A: Allogeneic CD8+/CD4+ T-Cell & DSA Evasion
    logit_g2a = pw * 1.25 + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "base_intercept", 0.10)
    if allo_k == "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO":
        logit_g2a += (
            pw * 6.875
            + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", 0.55)
            + signed_gene_signal(target_profiles, ["CIITA", "HLA-E", "CD47"], scale=0.30, branch_key="Gate2A_ALLO_B2M_HLAE_CD47")
        )
    elif allo_k == "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO":
        logit_g2a += (
            pw * 5.0
            + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO", 0.40)
            + signed_gene_signal(target_profiles, ["CIITA", "TYK2"], scale=0.28, is_knockout=True, branch_key="Gate2A_ALLO_HLA_C_TYK2")
            + signed_gene_signal(target_profiles, ("HLA-C",), scale=0.28, mode="preserve", branch_key="Gate2A_ALLO_HLA_C_TYK2")
        )
    elif allo_k == "ALLO_B2M_CIITA_KO_CD47_LENTI":
        logit_g2a += (
            pw * 4.75
            + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "ALLO_B2M_CIITA_KO_CD47_LENTI", 0.38)
            + signed_gene_signal(target_profiles, ["CIITA", "CD47"], scale=0.28, branch_key="Gate2A_ALLO_UP421")
        )
    elif allo_k == "ALLO_B2M_KO_PDL1_A20_KI":
        logit_g2a += (
            pw * -1.25
            + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "ALLO_B2M_KO_PDL1_A20_KI", -0.10)
            + signed_gene_signal(target_profiles, ["CD274", "TNFAIP3"], scale=0.25, branch_key="Gate2A_ALLO_VCTX210")
        )
    elif allo_k == "ALLO_NONE_UNEDITED_OR_SYSTEMIC_IMMUNOSUPPRESSION":
        logit_g2a += (
            (pw * -7.5 + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "ALLO_NONE_MACRODEVICE", -0.60))
            if niche_k == "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION"
            else (pw * -22.5 + _rdelta("Gate2A_Allo_Adaptive_T_and_DSA_Stealth", "ALLO_NONE_NAKED", -1.80))
        )

    # Gate 2B: Autoimmune T1D Memory & Indirect APC Presentation Stealth
    logit_g2b = (
        (pw * 2.5 + _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "base_intercept_pediatric", 0.20))
        if recipient_context == "Pediatric_Early_Onset_T1D"
        else (pw * 3.75 + _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "base_intercept_lada", 0.30))
    )
    if allo_k in ("ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", "ALLO_B2M_CIITA_KO_CD47_LENTI", "ALLO_B2M_KO_PDL1_A20_KI"):
        logit_g2b += (
            pw * 6.875
            - pw * 3.125
            + _pdelta("indirect_apc_presentation_residual_hazard_cap", -0.25, "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "B2M_NULL_DIRECT_CLOAK_PLUS_INDIRECT_APC_DISCOUNT", 0.30)
            + _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "B2M_NULL_DIRECT_CLOAK_PLUS_INDIRECT_APC_DISCOUNT", 0.30)
            + signed_gene_signal(target_profiles, ["CIITA", "CD47"], scale=0.22, branch_key="Gate2B_ALLO_B2M_NULL")
        )
    elif allo_k == "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO":
        if auto_k == "AUTO_NONE":
            logit_g2b -= pw * 6.875 - _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "HLA_C_RETAINED_AUTO_NONE", -0.55)
        elif auto_k == "AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT":
            logit_g2b += (
                pw * 3.125
                + _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "HLA_C_RETAINED_AUTO_INS_EDIT", 0.25)
                + signed_gene_signal(target_profiles, ["TYK2"], scale=0.28, is_knockout=True, branch_key="Gate2B_HLA_C_INS_EDIT")
                + signed_gene_signal(target_profiles, ["INS"], scale=0.28, mode="base_edit", branch_key="Gate2B_HLA_C_INS_EDIT")
            )
        else:
            logit_g2b += (
                pw * 2.5
                + _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "HLA_C_RETAINED_AUTO_FULL_KO", 0.20)
                + signed_gene_signal(target_profiles, ["G6PC2", "SLC30A8"], scale=0.25, is_knockout=True, branch_key="Gate2B_HLA_C_FULL_KO")
                + signed_gene_signal(target_profiles, ["INS"], scale=0.25, mode="base_edit", branch_key="Gate2B_HLA_C_FULL_KO")
            )
    elif allo_k == "ALLO_NONE_UNEDITED_OR_SYSTEMIC_IMMUNOSUPPRESSION":
        logit_g2b -= pw * 13.75 - _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "ALLO_NONE", -1.10)

    if cyto_k == "CYTO_CXCL10_SWAP_PDL1_SERPINB9":
        logit_g2b += (
            pw * 4.375
            + _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "CYTO_CXCL10_SWAP_PDL1_SERPINB9", 0.35)
            + signed_gene_signal(target_profiles, ["CD274", "SERPINB9", "CXCL10"], scale=0.25, branch_key="Gate2B_CYTO_CXCL10_SWAP")
        )
    elif cyto_k == "CYTO_CONSTITUTIVE_IL10_VIRAL_OE":
        logit_g2b -= (
            pw * 9.375
            - _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "CYTO_CONSTITUTIVE_IL10_VIRAL_OE", -0.75)
            - signed_gene_signal(target_profiles, ["IL10"], scale=0.25, branch_key="Gate2B_CYTO_IL10_OE")
        )
    elif cyto_k == "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1":
        logit_g2b -= pw * 10.625 - _rdelta("Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth", "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1", -0.85)

    # Gate 3: Innate NK Missing-Self, Macrophage Phagocytosis & IBMIR Shield
    logit_g3 = pw * 2.5 + _rdelta("Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "base_intercept", 0.20)
    if allo_k == "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO":
        logit_g3 += (
            pw * 6.875
            - pw * 2.5
            + _pdelta("resting_b2m_promoter_hour0_density_discount", -0.20, "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", 0.35)
            + _rdelta("Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO", 0.35)
            + signed_gene_signal(target_profiles, ["HLA-E", "CD47"], scale=0.35, branch_key="Gate3_ALLO_B2M_HLAE_CD47")
        )
    elif allo_k == "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO":
        logit_g3 += (
            pw * 4.375
            + _rdelta("Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO", 0.35)
            + signed_gene_signal(target_profiles, ["TYK2"], scale=0.30, is_knockout=True, branch_key="Gate3_ALLO_HLA_C_TYK2")
            + signed_gene_signal(target_profiles, ("HLA-C",), scale=0.30, mode="preserve", branch_key="Gate3_ALLO_HLA_C_TYK2")
        )
    elif allo_k == "ALLO_B2M_CIITA_KO_CD47_LENTI":
        logit_g3 += (
            pw * 3.75
            + _rdelta("Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "ALLO_B2M_CIITA_KO_CD47_LENTI", 0.30)
            + signed_gene_signal(target_profiles, ["CD47"], scale=0.28, branch_key="Gate3_ALLO_UP421")
        )
    elif allo_k == "ALLO_B2M_KO_PDL1_A20_KI":
        logit_g3 -= pw * 9.375 - _rdelta("Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "ALLO_B2M_KO_PDL1_A20_KI", -0.75)
    elif allo_k == "ALLO_NONE_UNEDITED_OR_SYSTEMIC_IMMUNOSUPPRESSION":
        logit_g3 -= pw * 3.125 - _rdelta("Gate3_Innate_NK_Macrophage_and_IBMIR_Shield", "ALLO_NONE", -0.25)

    # Gate 4: Intrinsic Cytokine, Granzyme B, ER Stress & Amyloid Resilience
    logit_g4 = pw * 1.25 + _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "base_intercept", 0.10)
    if cyto_k == "CYTO_CXCL10_SWAP_PDL1_SERPINB9":
        logit_g4 += (
            pw * 5.625
            + _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "CYTO_CXCL10_SWAP_PDL1_SERPINB9", 0.45)
            + signed_gene_signal(target_profiles, ["CD274", "SERPINB9", "CXCL10"], scale=0.32, branch_key="Gate4_CYTO_CXCL10_SWAP")
        )
    elif cyto_k == "CYTO_CONSTITUTIVE_IL10_VIRAL_OE":
        logit_g4 -= (
            pw * 12.5
            - _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "CYTO_CONSTITUTIVE_IL10_VIRAL_OE", -1.00)
            - signed_gene_signal(target_profiles, ["IL10"], scale=0.30, branch_key="Gate4_CYTO_IL10_OE")
        )
    elif cyto_k == "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1":
        logit_g4 -= (
            pw * 16.875
            - _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1", -1.35)
            - signed_gene_signal(target_profiles, ["STAT1"], scale=0.30, is_knockout=True, branch_key="Gate4_CYTO_STAT1_KO")
        )
    elif cyto_k == "CYTO_NONE":
        if allo_k == "ALLO_B2M_KO_PDL1_A20_KI":
            logit_g4 += (
                pw * 1.875
                + _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "CYTO_NONE_WITH_VCTX210_A20", 0.15)
                + signed_gene_signal(target_profiles, ["TNFAIP3", "CD274"], scale=0.25, branch_key="Gate4_CYTO_VCTX210_A20")
            )
        else:
            logit_g4 -= pw * 8.125 - _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "CYTO_NONE_DEFAULT", -0.65)

    if er_k == "ER_RNLS_KO":
        logit_g4 += (
            pw * 2.25
            + _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "ER_RNLS_KO", 0.18)
            + signed_gene_signal(target_profiles, ["RNLS"], scale=0.25, is_knockout=True, branch_key="Gate4_ER_RNLS_KO")
        )
    elif er_k == "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE":
        logit_g4 += (
            pw * 3.125
            + _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "ER_RNLS_KO_PLUS_IAPP_PRAMLINTIDE", 0.25)
            + signed_gene_signal(target_profiles, ["RNLS"], scale=0.25, is_knockout=True, branch_key="Gate4_ER_RNLS_KO")
            + signed_gene_signal(target_profiles, ["IAPP"], scale=0.25, mode="base_edit", branch_key="Gate4_ER_RNLS_KO")
        )
    elif niche_k == "NICHE_TRIHORMONAL_120UM_MVF_CAO2":
        logit_g4 += pw * 3.125 + _rdelta("Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience", "ER_NONE_WITH_STAGE7_120UM_NICHE", 0.25)

    # Gate 5: 72h Vascular/Oxygen Engraftment & 80-Year Pan-Lineage Prodrug Erasability
    logit_g5 = pw * 1.875 + _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "base_intercept", 0.15)
    if niche_k == "NICHE_TRIHORMONAL_120UM_MVF_CAO2":
        logit_g5 += (
            pw * 5.625
            + _pdelta("stage7_trihormonal_120um_mvf_cao2_niche_prior", 0.45, "Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "NICHE_TRIHORMONAL_120UM_MVF_CAO2", 0.45)
            + _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "NICHE_TRIHORMONAL_120UM_MVF_CAO2", 0.45)
        )
    elif niche_k == "NICHE_NAKED_AVASCULAR_TRANSPLANT":
        logit_g5 -= pw * 6.875 - _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "NICHE_NAKED_AVASCULAR_TRANSPLANT", -0.55)
    elif niche_k == "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION":
        logit_g5 -= pw * 11.875 - _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION", -0.95)

    if saf_k == "SAF_AAVS1_CAG_ICASP9_BIALLELIC":
        logit_g5 += (
            pw * 4.75
            + _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "SAF_AAVS1_CAG_ICASP9_BIALLELIC", 0.38)
            + signed_gene_signal(target_profiles, ["CASP9"], scale=0.25, mode="gof", branch_key="Gate5_SAF_ICASP9")
        )
    elif saf_k == "SAF_G6PC2_ENDOGENOUS_ICASP9":
        logit_g5 += (
            pw * -5.625
            + _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "SAF_G6PC2_ENDOGENOUS_ICASP9", -0.45)
            + signed_gene_signal(target_profiles, ["G6PC2"], scale=0.25, is_knockout=True, branch_key="Gate5_SAF_ICASP9")
            + signed_gene_signal(target_profiles, ["CASP9"], scale=0.25, mode="gof", branch_key="Gate5_SAF_ICASP9")
        )
    elif saf_k == "SAF_NONE":
        logit_g5 -= (
            (pw * 5.625 - _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "SAF_NONE_PEDIATRIC", -0.45))
            if recipient_context == "Pediatric_Early_Onset_T1D"
            else (pw * 3.125 - _rdelta("Gate5_72h_Vascular_Oxygen_and_80yr_Safety", "SAF_NONE_LADA", -0.25))
        )

    _gate_sigmoid = sigmoid_unclipped if sigmoid is _DEFAULT_CLAMPED_SIGMOID else sigmoid
    g1 = _gate_sigmoid(logit_g1)
    g2a = _gate_sigmoid(logit_g2a)
    g2b = _gate_sigmoid(logit_g2b)
    g3 = _gate_sigmoid(logit_g3)
    g4 = _gate_sigmoid(logit_g4)
    g5 = _gate_sigmoid(logit_g5)

    g2_combined = g2a * g2b
    unpenalized_product = g1 * g2_combined * g3 * g4 * g5

    # Multiplex DSB translocation penalty + nickase base-edit penalty wired to ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE
    lambda_k = _finite_real(
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.get("locus_dsb_penalty_rate_lambda_K", {}).get("point_estimate", 0.15), "Locus penalty rate"
    )
    gamma_b = _finite_real(
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.get("base_edit_nickase_penalty_factor_gamma_B", {}).get("point_estimate", 0.965), "Base-edit penalty factor"
    )
    mfg_penalty = engineering_complexity_discount(total_cuts, total_base_edits, lambda_k, gamma_b)
    penalized_index = unpenalized_product * mfg_penalty
    geom_mean = unpenalized_product ** (1.0 / 5.0)
    geom_mean_six_subgates = unpenalized_product ** (1.0 / 6.0)

    gates_unrounded = {
        "Gate1_Maturation_GSIS_and_Paracrine": g1,
        "Gate2A_Allo_Adaptive_T_and_DSA_Stealth": g2a,
        "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth": g2b,
        "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield": g3,
        "Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience": g4,
        "Gate5_72h_Vascular_Oxygen_and_80yr_Safety": g5,
    }
    gate_logits_unrounded = {
        "Gate1_Maturation_GSIS_and_Paracrine": logit_g1,
        "Gate2A_Allo_Adaptive_T_and_DSA_Stealth": logit_g2a,
        "Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth": logit_g2b,
        "Gate3_Innate_NK_Macrophage_and_IBMIR_Shield": logit_g3,
        "Gate4_Intrinsic_Cytokine_ER_Granzyme_and_Amyloid_Resilience": logit_g4,
        "Gate5_72h_Vascular_Oxygen_and_80yr_Safety": logit_g5,
    }
    gates = {k: round(v, 4) for k, v in gates_unrounded.items()}
    bottleneck_gate = min(gates_unrounded.items(), key=lambda kv: kv[1])

    # Propagate gene-level posterior SEs into the 10th-percentile lower credible bound (AT-18, live-wired to lower_10pct_cred_bound_hyperparameters per Audit R10)
    lb_params = _aux_param("lower_10pct_cred_bound_hyperparameters", {})
    if not isinstance(lb_params, dict):
        raise ValueError("Heuristic floor parameters must be an explicit mapping")
    z_mult = _finite_real(lb_params.get("z_multiplier", 1.2816), "Floor z multiplier")
    cut_disc = _finite_real(lb_params.get("locus_cut_linear_discount", 0.008), "Floor locus discount")
    def_g_se = _finite_real(lb_params.get("default_gene_se", 0.35), "Default gene dispersion")
    no_g_se = _finite_real(lb_params.get("no_gene_se", 0.45), "Missing gene dispersion")
    if min(z_mult, cut_disc, def_g_se, no_g_se) < 0:
        raise ValueError("Heuristic floor discounts and dispersions must be nonnegative")

    stack_genes = sorted(
        set(
            ALLO_CLOAK_MODULES[allo_k]["genes"]
            + AUTOANTIGEN_MODULES[auto_k]["genes"]
            + CYTOKINE_GRANZYME_MODULES[cyto_k]["genes"]
            + ER_AMYLOID_MODULES[er_k]["genes"]
            + NICHE_VASCULAR_MODULES[niche_k]["genes"]
            + SAFETY_SWITCH_MODULES[saf_k]["genes"]
        )
    )
    if stack_genes:
        errors = [_profile_standard_error(target_profiles.get(g, {}), def_g_se) for g in stack_genes]
        largest = max(errors)
        # Scaled RMS avoids squaring huge finite SEs or underflowing tiny SEs.
        se_combined = largest * math.sqrt(sum((se / largest) ** 2 for se in errors) / len(errors)) if largest else 0.0
    else:
        se_combined = no_g_se
    b_quad_disc = _sec_param("lower_10pct_cred_bound_b_edit_quadratic_discount", 0.005)
    if b_quad_disc < 0:
        raise ValueError("Base-edit floor discount must be nonnegative")
    lower_10pct_bound = round(
        max(0.0, penalized_index * math.exp(-z_mult * se_combined) - cut_disc * total_cuts), 4
    )
    lower_10pct_bound_with_b_quad = round(
        max(0.0, penalized_index * math.exp(-z_mult * se_combined) - cut_disc * total_cuts - b_quad_disc * (total_base_edits ** 2)), 4
    )

    return {
        "gate_scores": gates,
        "gate_logits_unrounded": gate_logits_unrounded,
        "unpenalized_gate_product": round(unpenalized_product, 4),
        "geometric_mean_gate_score": round(geom_mean, 4),
        "geometric_mean_six_subgate_score": round(geom_mean_six_subgates, 4),
        "limiting_bottleneck_gate": bottleneck_gate[0],
        "limiting_bottleneck_score": round(bottleneck_gate[1], 4),
        "manufacturing_locus_penalty_factor": round(mfg_penalty, 4),
        "series_system_joint_survival_probability": round(penalized_index, 4),
        "series_system_joint_survival_probability_unrounded": penalized_index,
        "lower_10pct_cred_bound": lower_10pct_bound,
        "lower_10pct_cred_bound_with_b_edit_quadratic_discount": lower_10pct_bound_with_b_quad,
        "heuristic_conservative_floor_index": lower_10pct_bound,
        "lower_bound_field_note": (
            "Heuristic conservative scenario floor index combining inverse-root evidence sparsity and locus burden "
            "(NOT a calibrated Bayesian 10th-percentile credible interval; retained under key lower_10pct_cred_bound for backward compatibility; "
            "lower_10pct_cred_bound_with_b_edit_quadratic_discount additionally subtracts 0.005 * B^2 for nickase base-edit stacks per Finding 013)."
        ),
        "estimand_note": "Prior-Anchored Multi-Gate Scenario Index (NOT a calibrated human survival probability).",
    }


def _combination_order_key(result: dict) -> tuple:
    return (
        result["is_feasible_candidate"],
        result["composite_series_system_utility_unrounded"],
        result["composite_series_system_utility"],
        result["unblocked_diagnostic_utility"],
        -result["total_genomic_cuts_K"],
        -result["total_base_edits_B"],
        "AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT" in result["combination_id"],
        result["combination_id"],
    )


def run_combinatorial_multi_gate_search(target_profiles: dict) -> list[dict]:
    results = []
    w_ped, w_adult = _composite_utility_weights()
    for allo_k, auto_k, cyto_k, er_k, niche_k, saf_k in product(
        ALLO_CLOAK_MODULES.keys(),
        AUTOANTIGEN_MODULES.keys(),
        CYTOKINE_GRANZYME_MODULES.keys(),
        ER_AMYLOID_MODULES.keys(),
        NICHE_VASCULAR_MODULES.keys(),
        SAFETY_SWITCH_MODULES.keys(),
    ):
        total_cuts = (
            ALLO_CLOAK_MODULES[allo_k]["cuts"]
            + AUTOANTIGEN_MODULES[auto_k]["cuts"]
            + CYTOKINE_GRANZYME_MODULES[cyto_k]["cuts"]
            + ER_AMYLOID_MODULES[er_k]["cuts"]
            + NICHE_VASCULAR_MODULES[niche_k]["cuts"]
            + SAFETY_SWITCH_MODULES[saf_k]["cuts"]
        )
        total_base_edits = (
            ALLO_CLOAK_MODULES[allo_k]["base_edits"]
            + AUTOANTIGEN_MODULES[auto_k]["base_edits"]
            + CYTOKINE_GRANZYME_MODULES[cyto_k]["base_edits"]
            + ER_AMYLOID_MODULES[er_k]["base_edits"]
            + NICHE_VASCULAR_MODULES[niche_k]["base_edits"]
            + SAFETY_SWITCH_MODULES[saf_k]["base_edits"]
        )

        ped = compute_series_system_utility(
            allo_k, auto_k, cyto_k, er_k, niche_k, saf_k, target_profiles, "Pediatric_Early_Onset_T1D"
        )
        lada = compute_series_system_utility(
            allo_k, auto_k, cyto_k, er_k, niche_k, saf_k, target_profiles, "Adult_Onset_LADA"
        )
        raw_composite = w_ped * ped["series_system_joint_survival_probability"] + w_adult * lada["series_system_joint_survival_probability"]
        graveyard_flags = evaluate_graveyard_rules_for_stack(
            allo_k, auto_k, cyto_k, er_k, niche_k, saf_k, total_base_edits
        )
        has_block = any(r.get("severity_level") == "BLOCK" for r in graveyard_flags)
        if has_block:
            ped["unblocked_diagnostic_survival_index"] = ped["series_system_joint_survival_probability"]
            lada["unblocked_diagnostic_survival_index"] = lada["series_system_joint_survival_probability"]
            ped["series_system_joint_survival_probability"] = 0.0
            ped["series_system_joint_survival_probability_unrounded"] = 0.0
            ped["lower_10pct_cred_bound"] = 0.0
            ped["lower_10pct_cred_bound_with_b_edit_quadratic_discount"] = 0.0
            ped["heuristic_conservative_floor_index"] = 0.0
            lada["series_system_joint_survival_probability"] = 0.0
            lada["series_system_joint_survival_probability_unrounded"] = 0.0
            lada["lower_10pct_cred_bound"] = 0.0
            lada["lower_10pct_cred_bound_with_b_edit_quadratic_discount"] = 0.0
            lada["heuristic_conservative_floor_index"] = 0.0
            composite = 0.0
            exact_comp_unrounded = 0.0
        else:
            composite = raw_composite
            exact_comp_unrounded = (
                w_ped * ped["series_system_joint_survival_probability_unrounded"]
                + w_adult * lada["series_system_joint_survival_probability_unrounded"]
            )

        results.append({
            "combination_id": f"{allo_k}__{auto_k}__{cyto_k}__{er_k}__{niche_k}__{saf_k}",
            "modules": {
                "allo_cloak": ALLO_CLOAK_MODULES[allo_k]["label"],
                "autoantigen_stealth": AUTOANTIGEN_MODULES[auto_k]["label"],
                "cytokine_and_granzyme_shield": CYTOKINE_GRANZYME_MODULES[cyto_k]["label"],
                "er_stress_and_amyloid_shield": ER_AMYLOID_MODULES[er_k]["label"],
                "vascular_and_oxygen_niche": NICHE_VASCULAR_MODULES[niche_k]["label"],
                "safety_switch": SAFETY_SWITCH_MODULES[saf_k]["label"],
            },
            "total_genomic_cuts_K": total_cuts,
            "total_base_edits_B": total_base_edits,
            "multiplex_base_edit_plus_dsb_translocation_risk_flag": bool(total_base_edits >= 2 and total_cuts >= 4),
            "unswitched_immune_evasive_clinical_hold": bool(
                allo_k != "ALLO_NONE_UNEDITED_OR_SYSTEMIC_IMMUNOSUPPRESSION" and saf_k == "SAF_NONE"
            ),
            "is_feasible_candidate": not has_block,
            "unblocked_diagnostic_utility": round(raw_composite, 4),
            "utility_pediatric_early_onset_t1d": ped,
            "utility_adult_onset_lada": lada,
            "composite_series_system_utility": round(composite, 4),
            "composite_series_system_utility_unrounded": exact_comp_unrounded,
            "triggered_graveyard_rules": graveyard_flags,
        })

    results.sort(key=_combination_order_key, reverse=True)
    rounded_tie_counts = Counter(
        (r["is_feasible_candidate"], r["composite_series_system_utility"]) for r in results
    )
    for idx, r in enumerate(results, 1):
        r["global_rank_among_all_combinations"] = idx
        r["rank_tie_group_size"] = rounded_tie_counts[
            (r["is_feasible_candidate"], r["composite_series_system_utility"])
        ]
    return results


def _build_odds_composition_and_profile_scope_audit_iip_0003() -> dict:
    """Re-judged four-hurdle conditional scenario specification, score/rank measure labels, and OQ-016 (IIP-0003)."""
    h1_1, h1_3 = 0.88, 0.88
    h2_1, h2_3 = 0.78, 0.82
    h3_1, h3_3 = 0.50, 0.72
    h4_1, h4_3 = 0.60, 0.72
    p1_exact = round(h1_1 * h2_1 * h3_1 * h4_1, 8)
    p3_exact = round(h1_3 * h2_3 * h3_3 * h4_3, 8)
    return {
        "proposal_id": "IIP-0003",
        "title": (
            "The published chance-of-working estimate does not compute the quantity it states, "
            "and omits failure modes the blueprint already documents"
        ),
        "layer": "Layer 1 & Cross-Blueprint Odds Specification",
        "open_question_id": "OQ-016",
        "open_question_title": "Which assumptions and definition of success sit behind the four-hurdle odds estimate?",
        "status": "accepted_in_open_draft",
        "reassessment_date": "2026-11-05",
        "credit": (
            "Linus Schumacher (ORCID 0000-0003-0797-3406, GitHub @ljschumacher, University of Edinburgh), "
            "IIP-0003, https://github.com/openislet/openislet/issues/3."
        ),
        "sequential_conditional_definition": {
            "one_design_label": "One design (Design A)",
            "one_design_rule": (
                "Reassessed and adopted on 2026-11-05 as a sequential conditional judgement: for Design A tested alone, "
                "each step is the judged chance that Design A clears that hurdle given that Design A cleared the earlier hurdles."
            ),
            "three_designs_label": "Three designs in parallel (A, B and X)",
            "three_designs_rule": (
                "Reassessed and adopted on 2026-11-05 as a sequential surviving-portfolio conditional judgement: for Designs A, B and X "
                "advanced in parallel, each step is the judged chance that at least one same design still in the running clears that hurdle "
                "(the same design must clear every hurdle so far)."
            ),
            "composition_rule": (
                "Chance of getting past every hurdle = hurdle 1 x hurdle 2 x hurdle 3 x hurdle 4, "
                "each the chance that a design still in the running clears it."
            ),
        },
        "selected_hurdle_scope_and_excluded_pre_trial_conditions": {
            "selection_basis": (
                "The four hurdles are a deliberate selection covering the four in-body biological steps that have published human or modelled anchors. "
                "They assume this product's cells work as designed and pass the lab safety checks before any trial, so they are not the chance of a cure."
            ),
            "condition_line_short": (
                "Assumes the edited cells work as designed and pass the lab safety and DNA checks before any trial."
            ),
            "assumed_passed_outside_the_four_hurdles": [
                {
                    "id": "GATE1_PRODUCT_FUNCTION",
                    "description": (
                        "Product-specific Gate 1 endocrine function: OpenIslet's own edited 60:30:10 cells retain regulated "
                        "glucose-stimulated insulin secretion and rapid low-glucose shut-off (Hurdle 1 scores the class-level "
                        "human stem-cell-islet precedent from zimislecel PMID:40544428; product-specific function is a stated condition)."
                    ),
                },
                {
                    "id": "OQ-012",
                    "description": (
                        "Pre-trial safety-switch qualification: biallelic AAVS1 iCasp9 clonogenic escape bound <= 1.0e-8 across "
                        "endocrine and vascular stromal lineages under combined 1% O2 hypoxia and inflammatory cytokines."
                    ),
                },
                {
                    "id": "OQ-006",
                    "description": (
                        "Pre-trial edited-line genomic stability at manufacturing scale (rule GY-006): >=63,000x target-enriched "
                        "duplex sequencing for TP53/PPM1D/BCOR, Chr20q11.21 ddPCR, 15-25 kb long-read LOH phasing, and CAST-Seq translocation checks."
                    ),
                },
            ],
        },
        "reassessed_hurdles": [
            {
                "hurdle_number": 1,
                "title": "Do lab-grown islet cells work as a treatment in people? (Gate 1)",
                "one_design_design_a": h1_1,
                "three_designs_parallel": h1_3,
                "decision": "kept_after_2026_11_05_reassessment",
                "reason_for_three_design_figure": (
                    "All three designs use the same lab-grown stem-cell islet cells and 60:30:10 cell mix, "
                    "so trying three designs adds no second route on this shared step."
                ),
            },
            {
                "hurdle_number": 2,
                "title": "Does the genetic cloak stop rejection without anti-rejection drugs? (Gates 2A and 3)",
                "one_design_design_a": h2_1,
                "three_designs_parallel": h2_3,
                "decision": "kept_after_2026_11_05_reassessment",
                "reason_for_three_design_figure": (
                    "Designs A and X share the same immune cloak, so only Design B (which keeps one matched HLA-C gene) "
                    "offers a genuinely different way past rejection if that cloak falls short."
                ),
            },
            {
                "hurdle_number": 3,
                "title": "Can the cells survive the autoimmune attack that caused the diabetes? (Gates 2B and 4)",
                "one_design_design_a": h3_1,
                "three_designs_parallel": h3_3,
                "decision": "kept_after_2026_11_05_reassessment",
                "reason_for_three_design_figure": (
                    "All three designs rely on the same local anti-inflammatory shield at CXCL10, while Designs X and B add an "
                    "insulin-flag edit and Design B also dampens interferon signalling, giving partial backup rather than three separate tries."
                ),
            },
            {
                "hurdle_number": 4,
                "title": "Does the oxygen layer keep the graft alive until blood vessels grow in? (Gate 5 oxygen bridge)",
                "one_design_design_a": h4_1,
                "three_designs_parallel": h4_3,
                "decision": "kept_after_2026_11_05_reassessment",
                "reason_for_three_design_figure": (
                    "All three designs use the same oxygen-releasing sheet and pre-built microvessels, while Design X tests a thinner "
                    "oxygen sheet and smaller clusters if the baseline layout runs short."
                ),
            },
        ],
        "multiplied_conditional_totals": {
            "one_design_exact": p1_exact,
            "one_design_rounded_pct": 21,
            "three_designs_exact": p3_exact,
            "three_designs_rounded_pct": 37,
        },
        "section3_score_and_rank_measure_reconciliation": {
            "overall_50_50_composite_score": {
                "Design_A": 0.1085,
                "Design_X_grid_2": 0.1046,
                "Design_B_grid_25": 0.0673,
                "definition": "Equal 50:50 weight on Pediatric_Early_Onset_T1D and Adult_Onset_LADA (composite_series_system_utility)",
            },
            "childrens_pediatric_t1d_score": {
                "Design_A": 0.1070,
                "Design_X_grid_2": 0.1032,
                "Design_B_grid_25": 0.0664,
                "definition": "Pediatric_Early_Onset_T1D cohort utility (utility_pediatric_early_onset_t1d.series_system_joint_survival_probability)",
            },
            "adults_lada_score": {
                "Design_A": 0.1100,
                "Design_X_grid_2": 0.1061,
                "Design_B_grid_25": 0.0683,
                "definition": "Adult_Onset_LADA cohort utility (utility_adult_onset_lada.series_system_joint_survival_probability)",
            },
            "rank_column_measure": "overall_50_50_composite_score",
            "rank_column_explanation": (
                "global_rank_among_all_combinations (#1 for Design A, #2 for Design X, #25 for Design B out of 2,160) "
                "is sorted on the overall 50:50 composite score (composite_series_system_utility_unrounded / composite_series_system_utility "
                "in _combination_order_key). On the children's score alone (Pediatric_Early_Onset_T1D), Design A ranks #1, Design X ranks #2, and Design B ranks #24."
            ),
        },
        "new_open_question": {
            "id": "OQ-016",
            "title": "Which assumptions and definition of success sit behind the four-hurdle odds estimate?",
            "plain": (
                "The 21% and 37% figures are rough team judgements for clearing four in-body hurdles in sequence, "
                "assuming the edited cells work as designed and pass lab safety checks before any trial. "
                "Moving from one design to three depends on how much the designs overlap at each hurdle, "
                "and the estimate depends on which clinical outcome counts as success."
            ),
            "detail": (
                "Raised in IIP-0003 (issue #3). Each step in the four-hurdle estimate is judged for a design still in the running "
                "after the earlier hurdles, and the three-design figure depends on how strongly Designs A, B and X share failure modes "
                "at each hurdle. Lab release gates (OQ-012 safety-switch escape, OQ-006 DNA rearrangements, and product-specific Gate 1 "
                "insulin and glucagon function) sit outside the four hurdles as pre-trial pass/fail conditions, and a single clinical endpoint "
                "(such as one-year insulin independence without anti-rejection drugs and without severe hypoglycaemia) must be fixed before "
                "replacing starting judgements with data."
            ),
            "settle": (
                "Fix one clinical definition of success, write down the per-design and shared failure assumptions at each hurdle, "
                "and replace the starting judgements with measured pass rates as lab and preclinical results arrive."
            ),
            "tests": ["1", "2A", "2B", "3", "4", "5"],
            "designs": ["A", "B", "X"],
            "linked_iip": "IIP-0003",
        },
    }


def _rank_precision_summary(a_ranks: list[int], b_ranks: list[int], scheme: str, *, max_rank: int | None = None) -> dict:
    """Finite simulation precision only; these are not clinical probabilities."""
    grid_size = math.prod(len(modules) for modules in (
        ALLO_CLOAK_MODULES, AUTOANTIGEN_MODULES, CYTOKINE_GRANZYME_MODULES,
        ER_AMYLOID_MODULES, NICHE_VASCULAR_MODULES, SAFETY_SWITCH_MODULES,
    ))
    max_rank = grid_size if max_rank is None else max_rank
    if type(max_rank) is not int or not 1 <= max_rank <= grid_size:
        raise ValueError("Rank population must be a positive integer within the configured grid")
    if not isinstance(a_ranks, list) or not isinstance(b_ranks, list) or not a_ranks or len(a_ranks) != len(b_ranks) or any(type(rank) is not int or rank < 1 for rank in a_ranks + b_ranks):
        raise ValueError("Rank precision requires aligned nonempty positive integer rank draws")
    if any(rank > max_rank for rank in a_ranks + b_ranks):
        raise ValueError("Rank draws must lie within their declared population")
    def wilson(k, n):
        z = 1.959963984540054
        p = k / n
        den = 1 + z*z/n
        center = (p + z*z/(2*n))/den
        half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
        return [center-half, center+half]
    def quantiles(values):
        v = sorted(values)
        def quantile(q):
            index = (len(v)-1)*q
            lo = int(index)
            return v[lo] + (v[min(lo+1,len(v)-1)]-v[lo])*(index-lo)
        return {str(q): quantile(q) for q in (0.025, 0.25, 0.5, 0.75, 0.975)}
    n = len(a_ranks)
    return {"random_seed": 20261022, "sampling_scheme": scheme,
            "precision_status": "EXPLORATORY_FINITE_MONTE_CARLO_NOT_CLINICAL_PROBABILITY",
            "scenario_a_rank1_frequency_wilson95": wilson(sum(r==1 for r in a_ranks),n),
            "scenario_a_top3_frequency_wilson95": wilson(sum(r<=3 for r in a_ranks),n),
            "scenario_a_rank_quantiles_linear": quantiles(a_ranks),
            "scenario_b_rank_quantiles_linear": quantiles(b_ranks)}


def _numeric_scalar_paths(value, prefix=""):
    if isinstance(value, dict):
        return [p for key, val in value.items() for p in _numeric_scalar_paths(val, f"{prefix}.{key}".strip("."))]
    if isinstance(value, list):
        return [p for idx, val in enumerate(value) for p in _numeric_scalar_paths(val, f"{prefix}[{idx}] pre".replace(" pre", ""))]
    return [prefix] if isinstance(value, (int,float)) and not isinstance(value,bool) else []


def _direct_output_score(combination: dict, gate_noise: list[float]) -> float:
    """Re-evaluate exact production utility under paired output-logit noise.

    Display-rounded manufacturing fields are never numerical inputs. The same
    recipient weights and engineering prior registry serve both estimands.
    """
    ped = combination["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"]
    adult = combination["utility_adult_onset_lada"]["gate_logits_unrounded"]
    if len(gate_noise) != len(ped) or set(ped) != set(adult):
        raise ValueError("Output-logit sensitivity requires aligned gates and one noise term per gate")
    if not combination.get("is_feasible_candidate", True):
        return 0.0
    if combination["utility_pediatric_early_onset_t1d"].get("estimand_note", "").startswith("Zero evidence loaded"):
        return 0.0
    mfg = engineering_complexity_discount(
        combination["total_genomic_cuts_K"], combination["total_base_edits_B"],
        _finite_real(ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["locus_dsb_penalty_rate_lambda_K"]["point_estimate"], "Locus penalty rate"),
        _finite_real(ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE["base_edit_nickase_penalty_factor_gamma_B"]["point_estimate"], "Base-edit penalty factor"),
    )
    _gate_sigmoid = sigmoid_unclipped if sigmoid is _DEFAULT_CLAMPED_SIGMOID else sigmoid
    ped_prod = adult_prod = 1.0
    for gate, noise in zip(ped, gate_noise):
        noise = _finite_real(noise, "Output-logit perturbation")
        ped_prod *= _gate_sigmoid(ped[gate] + noise)
        adult_prod *= _gate_sigmoid(adult[gate] + noise)
    w_ped, w_adult = _composite_utility_weights()
    return mfg * (w_ped * ped_prod + w_adult * adult_prod)


def compute_direct_output_logit_sensitivity_analysis(target_profiles, all_combos, studies=None, omics=None):
    """Seed-replayable evaluation with registry restoration even after errors."""
    import copy
    original_registry = copy.deepcopy(COMPLETE_GATE_COEFFICIENT_REGISTRY)
    original_priors = copy.deepcopy(ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE)
    try:
        return _compute_direct_output_logit_sensitivity_analysis(target_profiles, all_combos, studies, omics)
    finally:
        COMPLETE_GATE_COEFFICIENT_REGISTRY.clear()
        COMPLETE_GATE_COEFFICIENT_REGISTRY.update(original_registry)
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.clear()
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.update(original_priors)


def _compute_direct_output_logit_sensitivity_analysis(
    target_profiles: dict,
    all_combos: list[dict],
    studies: list[dict] | None = None,
    omics: dict | None = None,
) -> dict:
    """
    Computes (1) the exact variance decomposition across the uniform 2,160-stack grid (six pooled
    pediatric gate logits with all three gene-signal channels zeroed vs loaded) and (2) a direct
    output-logit Monte Carlo perturbation analysis adding N(0, sigma_logit^2) noise directly to the
    six output gate logits (resolving Audit A03 + AUDIT-REF full 8/8 auxiliary leaf sampling
    and upstream profile rebuild).
    """
    import random
    import statistics

    zero_profiles = {
        g: {
            **v,
            "net_evidence_and_constraint_score": 0.0,
            "net_evidence_ko_score": 0.0,
            "net_evidence_gof_score": 0.0,
            "net_evidence_preserve_score": 0.0,
            "weighted_deleterious_penalty_score": 0.0,
        }
        for g, v in target_profiles.items()
    }
    corpus_only_profiles = {
        g: {
            **v,
            "net_evidence_ko_score": round(float(v.get("net_evidence_ko_score", 0.0)) - float(v.get("omics_ko_constraint_adjustment", 0.0)), 4),
            "net_evidence_gof_score": round(float(v.get("net_evidence_gof_score", 0.0)) - float(v.get("omics_gof_constraint_adjustment", 0.0)), 4),
            "net_evidence_preserve_score": round(float(v.get("net_evidence_preserve_score", 0.0)) - float(v.get("omics_gof_constraint_adjustment", 0.0)), 4),
            "net_evidence_and_constraint_score": round(
                float(v.get("net_evidence_and_constraint_score", 0.0))
                - (
                    float(v.get("omics_gof_constraint_adjustment", 0.0))
                    if v.get("default_perturbation_direction") == "Gain_of_Function_or_KnockIn"
                    else float(v.get("omics_ko_constraint_adjustment", 0.0))
                ),
                4,
            ),
        }
        for g, v in target_profiles.items()
    }
    omics_only_profiles = {
        g: {
            **v,
            "net_evidence_ko_score": round(float(v.get("omics_ko_constraint_adjustment", 0.0)), 4),
            "net_evidence_gof_score": round(float(v.get("omics_gof_constraint_adjustment", 0.0)), 4),
            "net_evidence_preserve_score": round(float(v.get("omics_gof_constraint_adjustment", 0.0)), 4),
            "net_evidence_and_constraint_score": round(
                float(v.get("omics_gof_constraint_adjustment", 0.0))
                if v.get("default_perturbation_direction") == "Gain_of_Function_or_KnockIn"
                else float(v.get("omics_ko_constraint_adjustment", 0.0)),
                4,
            ),
            "weighted_deleterious_penalty_score": 0.0,
        }
        for g, v in target_profiles.items()
    }
    zero_combos = run_combinatorial_multi_gate_search(zero_profiles)
    zero_by_id = {c["combination_id"]: c for c in zero_combos}
    identifiers = [candidate["combination_id"] for candidate in all_combos]
    if len(identifiers) != len(set(identifiers)) or set(identifiers) != set(zero_by_id):
        raise ValueError("Sensitivity requires each evaluated grid combination exactly once")
    all_combos = sorted(all_combos, key=_combination_order_key, reverse=True)
    corpus_only_by_id = {c["combination_id"]: c for c in run_combinatorial_multi_gate_search(corpus_only_profiles)}
    omics_only_by_id = {c["combination_id"]: c for c in run_combinatorial_multi_gate_search(omics_only_profiles)}

    structural_logits = []
    actual_logits = []
    corpus_only_logits = []
    omics_only_logits = []
    for c in all_combos:
        z_c = zero_by_id[c["combination_id"]]
        co_c = corpus_only_by_id[c["combination_id"]]
        om_c = omics_only_by_id[c["combination_id"]]
        for gate_name, logit in c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"].items():
            actual_logits.append(float(logit))
            structural_logits.append(float(z_c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"][gate_name]))
            corpus_only_logits.append(float(co_c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"][gate_name]))
            omics_only_logits.append(float(om_c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"][gate_name]))

    delta_logits = [a - b for a, b in zip(actual_logits, structural_logits)]
    delta_corpus_only = [a - b for a, b in zip(corpus_only_logits, structural_logits)]
    delta_omics_only = [a - b for a, b in zip(omics_only_logits, structural_logits)]
    vs = statistics.pvariance(structural_logits)
    ve = statistics.pvariance(delta_logits)
    ve_corpus = statistics.pvariance(delta_corpus_only)
    ve_omics = statistics.pvariance(delta_omics_only)
    vt = statistics.pvariance(actual_logits)
    mean_s = statistics.mean(structural_logits)
    mean_e = statistics.mean(delta_logits)
    cov_se = sum((a - mean_s) * (b - mean_e) for a, b in zip(structural_logits, delta_logits)) / len(delta_logits)

    chassis_a_id = (
        "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO__AUTO_NONE__"
        "CYTO_CXCL10_SWAP_PDL1_SERPINB9__ER_NONE__"
        "NICHE_TRIHORMONAL_120UM_MVF_CAO2__SAF_AAVS1_CAG_ICASP9_BIALLELIC"
    )
    chassis_b_id = (
        "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO__AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT__"
        "CYTO_CXCL10_SWAP_PDL1_SERPINB9__ER_NONE__"
        "NICHE_TRIHORMONAL_120UM_MVF_CAO2__SAF_AAVS1_CAG_ICASP9_BIALLELIC"
    )

    feasible_combos = [c for c in all_combos if c.get("is_feasible_candidate", True)]
    mc_results = {}
    for sigma in (0.04, 0.15, 0.25, 0.50):
        rng = random.Random(20261022)
        n_draws = 100
        a_ranks = []
        b_ranks = []
        for _ in range(n_draws):
            scored = []
            for candidate in feasible_combos:
                noises = [rng.gauss(0.0, sigma) for _ in range(6)]
                comp = _direct_output_score(candidate, noises)
                scored.append((comp, -candidate["total_genomic_cuts_K"], -candidate["total_base_edits_B"], candidate["combination_id"]))
            scored.sort(reverse=True)
            rank_map = {item[3]: idx for idx, item in enumerate(scored, 1)}
            a_ranks.append(rank_map[chassis_a_id])
            b_ranks.append(rank_map[chassis_b_id])
        mc_results[f"sigma_logit_{sigma:.2f}"] = {
            **_rank_precision_summary(a_ranks, b_ranks, "independent per-stack six-gate Gaussian noise shared across recipient contexts", max_rank=len(feasible_combos)),
            "effective_output_logit_sd": sigma,
            "n_monte_carlo_draws": n_draws,
            "scenario_a_rank1_count": sum(1 for r in a_ranks if r == 1),
            "scenario_a_top3_count": sum(1 for r in a_ranks if r <= 3),
            "scenario_a_rank_median": float(statistics.median(a_ranks)),
            "scenario_a_rank_min": min(a_ranks),
            "scenario_a_rank_max": max(a_ranks),
            "scenario_b_median_rank": float(statistics.median(b_ranks)),
            "scenario_b_rank_min": min(b_ranks),
            "scenario_b_rank_max": max(b_ranks),
        }

    # Shared-parameter registry perturbation Monte Carlo (Audit R6 + AUDIT-REF + AUDIT-REF):
    # Perturbs all 48 literal gate-offset terms across the 6 gates, all 8/8 auxiliary_evaluator_hyperparameters
    # groups (39 auxiliary scalar leaves), and all 6/6 ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE entries,
    # rebuilding target_profiles inside each draw so upstream evidence/omics hyperparameter draws propagate
    # end-to-end into target_profiles -> gate logits -> stack rankings.
    import copy
    orig_reg = copy.deepcopy(COMPLETE_GATE_COEFFICIENT_REGISTRY)
    orig_priors = copy.deepcopy(ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE)
    shared_param_mc = {
        "sampled_scalar_leaf_registry_ultracode_v80_14": {
            "n_gate_logit_offset_leaves_sampled": sum(1 for name,group in orig_reg.items() if name.startswith("Gate") for entry in group.values() if isinstance(entry,dict) and "logit_value" in entry),
            "n_net_gate_offset_entries_after_composite_collapse": sum(1 for name, group in orig_reg.items() if name.startswith("Gate") for entry in group.values() if isinstance(entry,dict) and "logit_value" in entry),
            "n_auxiliary_hyperparameter_groups_sampled": "8 / 8",
            "n_auxiliary_scalar_leaves_sampled": len(_numeric_scalar_paths({key: val["parameter_value"] for key,val in orig_reg["auxiliary_evaluator_hyperparameters"].items()})),
            "n_engineering_prior_leaves_sampled": len(orig_priors),
            "rebuilds_upstream_target_profiles_per_draw": bool(studies is not None and omics is not None),
            "auxiliary_groups_sampled": [
                "signed_gene_signal_slope (1 scalar)",
                "signed_gene_signal_clip_bounds (2 scalars)",
                "gene_signal_scales_by_gate_and_branch (22 scalars)",
                "cd47_gsis_penalty_scale (1 scalar)",
                "lrr_shrinkage_hyperparameters (2 scalars: prior_weight, per_empirical_study_slope)",
                "effect_direction_weights (4 scalars: protective, deleterious, secretome_or_differentiation_penalty_addon, shrunk_lrr_multiplier)",
                "omics_constraint_weights (3 scalars: gnomad_pli_ko_slope, gnomad_loeuf_ko_slope, opentargets_t1d_slope)",
                "lower_10pct_cred_bound_hyperparameters (4 scalars: z_multiplier, locus_cut_linear_discount, default_gene_se, no_gene_se)",
            ],
        }
    }
    for sigma in (0.04, 0.15, 0.25, 0.50):
        rng_sp = random.Random(20261022)
        n_draws_sp = 50
        sp_a_ranks = []
        sp_b_ranks = []
        for _ in range(n_draws_sp):
            for gk, gdict in orig_reg.items():
                if isinstance(gdict, dict) and gk not in ("pw_base_scale", "auxiliary_evaluator_hyperparameters", "secondary_ranking_tie_break_and_display_constants"):
                    for ek, edict in gdict.items():
                        if isinstance(edict, dict) and "logit_value" in edict:
                            COMPLETE_GATE_COEFFICIENT_REGISTRY[gk][ek]["logit_value"] = (
                                float(edict["logit_value"]) + rng_sp.gauss(0.0, sigma)
                            )
            orig_aux = orig_reg["auxiliary_evaluator_hyperparameters"]
            aux_live = COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"]
            aux_live["signed_gene_signal_slope"]["parameter_value"] = max(
                0.10, float(orig_aux["signed_gene_signal_slope"]["parameter_value"]) * (1.0 + rng_sp.gauss(0.0, 0.25 * sigma))
            )
            hi_clip = max(1.5, float(orig_aux["signed_gene_signal_clip_bounds"]["parameter_value"][1]) * (1.0 + rng_sp.gauss(0.0, 0.15 * sigma)))
            aux_live["signed_gene_signal_clip_bounds"]["parameter_value"] = [-hi_clip, hi_clip]
            aux_live["cd47_gsis_penalty_scale"]["parameter_value"] = min(
                -0.02, float(orig_aux["cd47_gsis_penalty_scale"]["parameter_value"]) * (1.0 + rng_sp.gauss(0.0, 0.25 * sigma))
            )
            for b_key, b_val in orig_aux["gene_signal_scales_by_gate_and_branch"]["parameter_value"].items():
                aux_live["gene_signal_scales_by_gate_and_branch"]["parameter_value"][b_key] = max(
                    0.05, float(b_val) * (1.0 + rng_sp.gauss(0.0, 0.25 * sigma))
                )
            for s_key, s_val in orig_aux["lrr_shrinkage_hyperparameters"]["parameter_value"].items():
                aux_live["lrr_shrinkage_hyperparameters"]["parameter_value"][s_key] = max(
                    0.05, float(s_val) * (1.0 + rng_sp.gauss(0.0, 0.20 * sigma))
                )
            for d_key, d_val in orig_aux["effect_direction_weights"]["parameter_value"].items():
                aux_live["effect_direction_weights"]["parameter_value"][d_key] = max(
                    0.05, float(d_val) * (1.0 + rng_sp.gauss(0.0, 0.20 * sigma))
                )
            for o_key, o_val in orig_aux["omics_constraint_weights"]["parameter_value"].items():
                base_ov = float(o_val)
                perturbed_ov = base_ov * (1.0 + rng_sp.gauss(0.0, 0.20 * sigma))
                aux_live["omics_constraint_weights"]["parameter_value"][o_key] = (
                    min(-0.02, perturbed_ov) if base_ov < 0 else max(0.01, perturbed_ov)
                )
            for l_key, l_val in orig_aux["lower_10pct_cred_bound_hyperparameters"]["parameter_value"].items():
                aux_live["lower_10pct_cred_bound_hyperparameters"]["parameter_value"][l_key] = max(
                    0.001, float(l_val) * (1.0 + rng_sp.gauss(0.0, 0.15 * sigma))
                )
            for p_key, p_meta in orig_priors.items():
                lo_p, hi_p = float(p_meta["plausible_range"][0]), float(p_meta["plausible_range"][1])
                base_p = float(p_meta["point_estimate"])
                sd_p = 0.20 * (hi_p - lo_p) * sigma
                ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE[p_key]["point_estimate"] = min(
                    hi_p, max(lo_p, base_p + rng_sp.gauss(0.0, sd_p))
                )
            draw_profiles = (
                build_target_evidence_profiles(studies, omics)
                if (studies is not None and omics is not None)
                else target_profiles
            )
            draw_combos = run_combinatorial_multi_gate_search(draw_profiles)
            rank_by_cid = {c["combination_id"]: c["global_rank_among_all_combinations"] for c in draw_combos}
            sp_a_ranks.append(rank_by_cid[chassis_a_id])
            sp_b_ranks.append(rank_by_cid[chassis_b_id])
        shared_param_mc[f"sigma_shared_param_{sigma:.2f}"] = {
            **_rank_precision_summary(sp_a_ranks, sp_b_ranks, "one correlated registry/prior draw shared by every stack; upstream profiles rebuilt", max_rank=len(draw_combos)),
            "shared_registry_parameter_sd": sigma,
            "n_monte_carlo_draws": n_draws_sp,
            "scenario_a_rank1_count": sum(1 for r in sp_a_ranks if r == 1),
            "scenario_a_top3_count": sum(1 for r in sp_a_ranks if r <= 3),
            "scenario_a_rank_median": float(statistics.median(sp_a_ranks)),
            "scenario_a_rank_min": min(sp_a_ranks),
            "scenario_a_rank_max": max(sp_a_ranks),
            "scenario_b_median_rank": float(statistics.median(sp_b_ranks)),
            "scenario_b_rank_min": min(sp_b_ranks),
            "scenario_b_rank_max": max(sp_b_ranks),
        }
    COMPLETE_GATE_COEFFICIENT_REGISTRY.clear()
    COMPLETE_GATE_COEFFICIENT_REGISTRY.update(copy.deepcopy(orig_reg))
    ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.clear()
    ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE.update(copy.deepcopy(orig_priors))

    # Feasible-subset (480 non-blocked stacks) and within-top-tier variance decomposition (AUDIT-REF):
    feas_struct_logits = []
    feas_actual_logits = []
    top_tier_struct_logits = []
    top_tier_actual_logits = []
    for c in all_combos:
        if not c.get("is_feasible_candidate", True):
            continue
        z_c = zero_by_id[c["combination_id"]]
        is_top_tier = (
            "NICHE_TRIHORMONAL_120UM_MVF_CAO2" in c["combination_id"]
            and "SAF_AAVS1_CAG_ICASP9_BIALLELIC" in c["combination_id"]
            and (
                "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO" in c["combination_id"]
                or "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO" in c["combination_id"]
                or "ALLO_B2M_CIITA_KO_CD47_LENTI" in c["combination_id"]
            )
        )
        for gate_name, la in c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"].items():
            ls = z_c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"][gate_name]
            feas_actual_logits.append(float(la))
            feas_struct_logits.append(float(ls))
            if is_top_tier:
                top_tier_actual_logits.append(float(la))
                top_tier_struct_logits.append(float(ls))

    feas_delta = [a - b for a, b in zip(feas_actual_logits, feas_struct_logits)]
    vs_feas = statistics.pvariance(feas_struct_logits)
    ve_feas = statistics.pvariance(feas_delta)
    vt_feas = statistics.pvariance(feas_actual_logits)
    top_delta = [a - b for a, b in zip(top_tier_actual_logits, top_tier_struct_logits)]
    vs_top = statistics.pvariance(top_tier_struct_logits)
    ve_top = statistics.pvariance(top_delta)
    vt_top = statistics.pvariance(top_tier_actual_logits)

    k4b0_loaded_list = [
        c for c in all_combos
        if c.get("is_feasible_candidate", True)
        and c.get("total_genomic_cuts_K") == 4
        and c.get("total_base_edits_B") == 0
    ]
    k4b0_zero_list = [
        c for c in zero_combos
        if c.get("is_feasible_candidate", True)
        and c.get("total_genomic_cuts_K") == 4
        and c.get("total_base_edits_B") == 0
    ]
    k4_s_logits = []
    k4_a_logits = []
    for c in k4b0_loaded_list:
        zc = zero_by_id[c["combination_id"]]
        for g_name, logit in c["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"].items():
            k4_a_logits.append(float(logit))
            k4_s_logits.append(float(zc["utility_pediatric_early_onset_t1d"]["gate_logits_unrounded"][g_name]))

    k4_delta = [a - b for a, b in zip(k4_a_logits, k4_s_logits)]
    vs_k4 = statistics.pvariance(k4_s_logits)
    ve_k4 = statistics.pvariance(k4_delta)
    vt_k4 = statistics.pvariance(k4_a_logits)

    feas_zero_list = [c for c in zero_combos if c.get("is_feasible_candidate", True)]
    top25_overlap = len(
        set(c["combination_id"] for c in feasible_combos[:25])
        & set(c["combination_id"] for c in feas_zero_list[:25])
    )
    top50_overlap = len(
        set(c["combination_id"] for c in feasible_combos[:50])
        & set(c["combination_id"] for c in feas_zero_list[:50])
    )
    def _tie_aware_spearman_from_scores(loaded_list: list[dict], zero_list: list[dict]) -> float:
        def _avg_ranks_by_id(items: list[dict]) -> dict[str, float]:
            grouped: dict[float, list[int]] = {}
            for pos, it in enumerate(items, 1):
                sc = float(it.get("composite_series_system_utility_unrounded", it.get("composite_series_system_utility", 0.0)))
                grouped.setdefault(sc, []).append(pos)
            avg_for_sc = {sc: sum(positions) / len(positions) for sc, positions in grouped.items()}
            return {
                it["combination_id"]: avg_for_sc[float(it.get("composite_series_system_utility_unrounded", it.get("composite_series_system_utility", 0.0)))]
                for it in items
            }

        rx = _avg_ranks_by_id(loaded_list)
        ry = _avg_ranks_by_id(zero_list)
        cids = list(rx.keys())
        mx = sum(rx[cid] for cid in cids) / len(cids)
        my = sum(ry[cid] for cid in cids) / len(cids)
        num = sum((rx[cid] - mx) * (ry[cid] - my) for cid in cids)
        den_x = math.sqrt(sum((rx[cid] - mx) ** 2 for cid in cids))
        den_y = math.sqrt(sum((ry[cid] - my) ** 2 for cid in cids))
        return num / (den_x * den_y) if (den_x > 0 and den_y > 0) else 1.0

    feasible_spearman_rho = _tie_aware_spearman_from_scores(feasible_combos, feas_zero_list)
    k4b0_spearman_rho = _tie_aware_spearman_from_scores(k4b0_loaded_list, k4b0_zero_list)

    oat_prior_tornado = {}
    for pk, pmeta in orig_priors.items():
        lo_v, hi_v = float(pmeta["plausible_range"][0]), float(pmeta["plausible_range"][1])
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE[pk]["point_estimate"] = lo_v
        lo_map = {c["combination_id"]: c for c in run_combinatorial_multi_gate_search(target_profiles)}
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE[pk]["point_estimate"] = hi_v
        hi_map = {c["combination_id"]: c for c in run_combinatorial_multi_gate_search(target_profiles)}
        ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE[pk]["point_estimate"] = pmeta["point_estimate"]
        oat_prior_tornado[pk] = {
            "plausible_range": [lo_v, hi_v],
            "point_estimate": float(pmeta["point_estimate"]),
            "scenario_a_rank_at_low": lo_map[chassis_a_id]["global_rank_among_all_combinations"],
            "scenario_a_rank_at_high": hi_map[chassis_a_id]["global_rank_among_all_combinations"],
            "scenario_a_composite_at_low": lo_map[chassis_a_id]["composite_series_system_utility"],
            "scenario_a_composite_at_high": hi_map[chassis_a_id]["composite_series_system_utility"],
            "scenario_a_composite_swing": round(
                hi_map[chassis_a_id]["composite_series_system_utility"] - lo_map[chassis_a_id]["composite_series_system_utility"], 4
            ),
            "scenario_b_rank_at_low": lo_map[chassis_b_id]["global_rank_among_all_combinations"],
            "scenario_b_rank_at_high": hi_map[chassis_b_id]["global_rank_among_all_combinations"],
            "scenario_b_composite_at_low": lo_map[chassis_b_id]["composite_series_system_utility"],
            "scenario_b_composite_at_high": hi_map[chassis_b_id]["composite_series_system_utility"],
            "scenario_b_composite_swing": round(
                hi_map[chassis_b_id]["composite_series_system_utility"] - lo_map[chassis_b_id]["composite_series_system_utility"], 4
            ),
            "scenario_a_remains_rank1_at_both_bounds": (
                lo_map[chassis_a_id]["global_rank_among_all_combinations"] == 1
                and hi_map[chassis_a_id]["global_rank_among_all_combinations"] == 1
            ),
        }
    orig_pw_ent = copy.deepcopy(COMPLETE_GATE_COEFFICIENT_REGISTRY.get("pw_base_scale", {"logit_value": 0.08, "plausible_range": [0.05, 0.12]}))
    orig_pw = float(orig_pw_ent.get("logit_value", 0.08)) if isinstance(orig_pw_ent, dict) else float(orig_pw_ent)
    pw_range = list(orig_pw_ent.get("plausible_range", [0.05, 0.12])) if isinstance(orig_pw_ent, dict) else [0.05, 0.12]
    if isinstance(COMPLETE_GATE_COEFFICIENT_REGISTRY.get("pw_base_scale"), dict):
        COMPLETE_GATE_COEFFICIENT_REGISTRY["pw_base_scale"]["logit_value"] = float(pw_range[0])
    else:
        COMPLETE_GATE_COEFFICIENT_REGISTRY["pw_base_scale"] = float(pw_range[0])
    lo_pw_map = {c["combination_id"]: c for c in run_combinatorial_multi_gate_search(target_profiles)}
    if isinstance(COMPLETE_GATE_COEFFICIENT_REGISTRY.get("pw_base_scale"), dict):
        COMPLETE_GATE_COEFFICIENT_REGISTRY["pw_base_scale"]["logit_value"] = float(pw_range[1])
    else:
        COMPLETE_GATE_COEFFICIENT_REGISTRY["pw_base_scale"] = float(pw_range[1])
    hi_pw_map = {c["combination_id"]: c for c in run_combinatorial_multi_gate_search(target_profiles)}
    COMPLETE_GATE_COEFFICIENT_REGISTRY["pw_base_scale"] = orig_pw_ent
    oat_prior_tornado["pw_base_scale"] = {
        "plausible_range": pw_range,
        "point_estimate": orig_pw,
        "scenario_a_rank_at_low": lo_pw_map[chassis_a_id]["global_rank_among_all_combinations"],
        "scenario_a_rank_at_high": hi_pw_map[chassis_a_id]["global_rank_among_all_combinations"],
        "scenario_a_composite_at_low": lo_pw_map[chassis_a_id]["composite_series_system_utility"],
        "scenario_a_composite_at_high": hi_pw_map[chassis_a_id]["composite_series_system_utility"],
        "scenario_a_composite_swing": round(
            hi_pw_map[chassis_a_id]["composite_series_system_utility"] - lo_pw_map[chassis_a_id]["composite_series_system_utility"], 4
        ),
        "scenario_b_rank_at_low": lo_pw_map[chassis_b_id]["global_rank_among_all_combinations"],
        "scenario_b_rank_at_high": hi_pw_map[chassis_b_id]["global_rank_among_all_combinations"],
        "scenario_b_composite_at_low": lo_pw_map[chassis_b_id]["composite_series_system_utility"],
        "scenario_b_composite_at_high": hi_pw_map[chassis_b_id]["composite_series_system_utility"],
        "scenario_b_composite_swing": round(
            hi_pw_map[chassis_b_id]["composite_series_system_utility"] - lo_pw_map[chassis_b_id]["composite_series_system_utility"], 4
        ),
        "scenario_a_remains_rank1_at_both_bounds": (
            lo_pw_map[chassis_a_id]["global_rank_among_all_combinations"] == 1
            and hi_pw_map[chassis_a_id]["global_rank_among_all_combinations"] == 1
        ),
    }

    top10_loaded = [c["combination_id"] for c in all_combos[:10]]
    top10_zero = [c["combination_id"] for c in zero_combos[:10]]
    top10_overlap = len(set(top10_loaded) & set(top10_zero))
    top10_exact_order_matches = sum(1 for a, b in zip(top10_loaded, top10_zero) if a == b)

    return {
        "methodology_note": (
            "Explicit variance decomposition across the uniform 2,160-stack grid (six pooled pediatric gate logits, "
            "comparing loaded empirical gene-signal channels against all three gene-signal channels zeroed), "
            "independent per-stack output-logit perturbation sensitivity, correlated shared-parameter registry "
            "Monte Carlo sensitivity (perturbing COMPLETE_GATE_COEFFICIENT_REGISTRY and lambda_K/gamma_B once per draw "
            "so all stacks sharing a module shift together), and deterministic One-at-a-Time (OAT) plausible-range prior "
            "tornado diagnostics across all 6 structural priors + pw_base_scale. Also reports exact top-10 membership overlap when all "
            "empirical gene signals are zeroed (F-04: structural priors and Graveyard rules determine macro-tier "
            "membership, whereas empirical gene signals act as local within-tier tie-breakers)."
        ),
        "engineering_priors_plausible_range_tornado_diagnostics": oat_prior_tornado,
        "uniform_grid_pooled_pediatric_gate_logit_variance_decomposition": {
            "design": "Uniform 2,160-stack grid, 6 pooled pediatric gate logits (12,960 gate-logit evaluations), all gene-signal channels zeroed vs loaded",
            "var_structural_priors": round(vs, 6),
            "var_empirical_gene_signal_delta": round(ve, 6),
            "var_empirical_gene_signal_delta_corpus_only": round(ve_corpus, 6),
            "var_empirical_gene_signal_delta_omics_only": round(ve_omics, 6),
            "twice_covariance": round(2.0 * cov_se, 6),
            "var_total": round(vt, 6),
            "structural_share_of_component_variance_sum": round(vs / (vs + ve), 6),
            "empirical_gene_signal_share_of_component_variance_sum": round(ve / (vs + ve), 6),
            "empirical_corpus_only_share_of_component_variance_sum": round(ve_corpus / (vs + ve_corpus), 6),
            "empirical_omics_only_share_of_component_variance_sum": round(ve_omics / (vs + ve_omics), 6),
            "structural_share_of_total_variance": round(vs / vt, 6),
            "empirical_gene_signal_share_of_total_variance": round(ve / vt, 6),
            "twice_covariance_share_of_total_variance": round((2.0 * cov_se) / vt, 6),
            "top10_membership_overlap_loaded_vs_zero_gene_signal": f"{top10_overlap} / 10",
            "top10_exact_rank_matches_loaded_vs_zero_gene_signal": f"{top10_exact_order_matches} / 10",
        },
        "feasible_subset_and_within_top_tier_variance_decomposition_sol_v80_r2_14": {
            "feasible_subset_n_stacks": len(feasible_combos),
            "feasible_var_structural_priors": round(vs_feas, 6),
            "feasible_var_empirical_gene_signal_delta": round(ve_feas, 6),
            "feasible_var_total": round(vt_feas, 6),
            "feasible_structural_share_of_component_variance_sum": round(vs_feas / (vs_feas + ve_feas), 6),
            "feasible_empirical_share_of_component_variance_sum": round(ve_feas / (vs_feas + ve_feas), 6),
            "feasible_structural_share_of_total_variance": round(vs_feas / vt_feas, 6),
            "feasible_spearman_rho_loaded_vs_zero_gene_signal": round(feasible_spearman_rho, 6),
            "top25_membership_overlap_loaded_vs_zero_gene_signal": f"{top25_overlap} / 25",
            "top50_membership_overlap_loaded_vs_zero_gene_signal": f"{top50_overlap} / 50",
            "top_tier_code_defined_n_stacks": len(top_tier_actual_logits) // 6,
            "top_tier_var_structural_priors": round(vs_top, 6),
            "top_tier_var_empirical_gene_signal_delta": round(ve_top, 6),
            "top_tier_twice_covariance": round(vt_top - vs_top - ve_top, 6),
            "top_tier_var_total": round(vt_top, 6),
            "top_tier_structural_share_of_component_variance_sum": round(vs_top / (vs_top + ve_top), 6),
            "top_tier_empirical_share_of_component_variance_sum": round(ve_top / (vs_top + ve_top), 6),
            "top_tier_structural_share_of_total_variance": round(vs_top / vt_top, 6),
            "feasible_k4_b0_n_stacks": len(k4b0_loaded_list),
            "feasible_k4_b0_var_structural_priors": round(vs_k4, 6),
            "feasible_k4_b0_var_empirical_gene_signal_delta": round(ve_k4, 6),
            "feasible_k4_b0_twice_covariance": round(vt_k4 - vs_k4 - ve_k4, 6),
            "feasible_k4_b0_var_total": round(vt_k4, 6),
            "feasible_k4_b0_structural_share_of_component_variance_sum": round(vs_k4 / (vs_k4 + ve_k4), 6),
            "feasible_k4_b0_empirical_share_of_component_variance_sum": round(ve_k4 / (vs_k4 + ve_k4), 6),
            "feasible_k4_b0_structural_share_of_total_variance": round(vs_k4 / vt_k4, 6),
            "feasible_k4_b0_spearman_rho_loaded_vs_zero_gene_signal": round(k4b0_spearman_rho, 6),
        },
        "direct_output_gate_logit_perturbation_monte_carlo": mc_results,
        "shared_parameter_registry_monte_carlo": shared_param_mc,
    }


CANDIDATE_ARCHITECTURES = [
    {
        "architecture_id": "ISLET1_CHASSIS_A_UNIVERSAL_LEAN_4LOCUS",
        "name": "Islet-1 Universal Chassis A ('Lean 4-Locus Swap-In Core + Biallelic AAVS1-[455bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE)]-EEF1A1-core-iCasp9')",
        "matched_combinatorial_search_id": (
            "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO__AUTO_NONE__"
            "CYTO_CXCL10_SWAP_PDL1_SERPINB9__ER_NONE__"
            "NICHE_TRIHORMONAL_120UM_MVF_CAO2__SAF_AAVS1_CAG_ICASP9_BIALLELIC"
        ),
        "genomic_locus_edit_count_K": 4,
        "base_edit_count_B": 0,
        "graveyard_warn_reconciliation_note": (
            "Triggers GRAVEYARD_009 (WARN: resting B2M-promoter Hour-0 surface density deficit) and GRAVEYARD_010 "
            "(WARN: MHC-null indirect APC bystander blindspot). In the evaluator, WARN rules apply pre-registered "
            "quantitative logit discounts (-0.20 in Gate 3 for GRAVEYARD_009; -0.25 in Gate 2B for GRAVEYARD_010) "
            "and attach mandatory Stage-0/Stage-1 Go/No-Go decision gates (Decision Experiment #1 QuantiBRITE flow + "
            "Arm 1B 2A-free bidirectional comparison; Decision Experiment #2 HLA-DQ2/DQ8 DC indirect presentation & "
            "incremental SA-FasL Delta_TUNEL+ <= +5.0 pp screen) without hard-zeroing the stack (Finding 040)."
        ),
        "locus1_and_locus4_engineering_audit_notes": {
            "locus1_arm1a_vs_arm1b_donor_size_and_delivery": (
                "Arm 1A bicistronic donor (4,880 bp donor / 5,170 bp with ITRs) exceeds the 4,700 bp single-AAV6 packaging limit "
                "and requires non-viral linear dsDNA/ceDNA HDR or Bridge-RNA recombinase (PMID:36008610, PMID:38926615, PMID:41642947) "
                "while leaving a 21-aa cytosolic -GSG-P2A tail on HLA-E*01:03; Arm 1B scarless 2A-free bidirectional cassette "
                "[SV40pA-CD47-IgV-GPI(antisense) <- EEF1A1-core | Syn_5xISRE_minP -> HLA-E*01:03(VMAPRTLVL)-bGHpA] is 3,974 bp donor / "
                "4,264 bp with ITRs (+436 bp headroom below the 4,700 bp single-AAV6 limit) and eliminates all 2A peptide tails (Findings 075, 076, 077)."
            ),
            "hla_e_vmaprtlvl_hcmv_ul40_and_kw_motif_note": (
                "VMAPRTLVL matches the HLA-A*02:01 leader and has ~12% sequence overlap with clinical HCMV UL40 isolates (vs ~43% for VMAPRTLIL; "
                "Lin et al. 2023 PMID:37264229), reducing but not eliminating NKG2C+ adaptive NK recognition in HCMV-seropositive recipients; "
                "Stage-0 Decision Experiment #1 compares WT HLA-E*01:03(VMAPRTLVL) against cytoplasmic tail KW-motif retention-signal removal "
                "(He et al. 2025 PMC12582296) across HCMV+ NKG2C+ and NKG2A- KIR2DL1/2/3+ primary NK donors (Finding 078)."
            ),
            "cd47_igv_gpi_vs_dct_ectodomain_note": (
                "CD47-IgV-GPI was evaluated in non-islet systems (Jurkat/macrophage phagocytosis, HUVEC angiogenesis, and hematopoietic stem-cell/immunodeficient-mouse assays; Xu et al. 2024 Cell Rep Med PMID:38508139, PMC10983038) rather than mature beta-cells; "
                "Stage-0 Decision Experiment #1 compares CD47-IgV-GPI against cytoplasmically truncated transmembrane CD47-dCT and full-length "
                "CD47-WT (OA3-323; Ghimire et al. 2023 PMID:37820008 Cdc42/GSIS counter-evidence) for GPI-PLC sensitivity, QuantiBRITE density, "
                "and dynamic perifusion GSIS (Finding 079)."
            ),
            "sa_fasl_sa_thbd_immunosuppression_and_fratricide_note": (
                "Published murine and NHP SA-FasL PEG-microgel studies (Headen et al. 2018 PMID:29867165; Lei et al. 2022 PMID:35559682) included transient short-course "
                "rapamycin, whereas Turan et al. 2024 (PMID:38709159) evaluated direct islet surface biotinylation co-displaying SA-THBD (SA-TM) and SA-CD47 in a murine syngeneic intraportal IBMIR/engraftment model without rapamycin and Kazmi et al. 2026 (PMID:41842939) demonstrated immunosuppression-free murine BALB/c->C57BL/6 intraportal allograft survival with "
                "direct surface SA-THBD + SA-CD47 (without SA-FasL, PEG microgels, or rapamycin); because IL-1beta sensitizes human beta-cells to Fas (CD95) apoptosis (Loweth et al. 1998 PMID:9588443), "
                "Decision Experiment #2 enforces an incremental fratricide kill gate of Delta_TUNEL+ <= +5.0 pp above cytokine-matched SA-THBD-only control AND anti-Fas ZB4 mAb rescue <= 5.0 pp (Finding 080; L3B-01/02/03 [DEC-077])."
            ),
            "locus4_a2ucoe_xiap_and_smac_salvage_note": (
                "Locus 4 uses the 455 bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE alias; insulation function pending) + 250 bp cHS4 insulator driving the 230 bp human EEF1A1 core promoter "
                "(3,457 bp donor / 3,747 bp with ITRs <= 4,700 bp single-AAV6 limit; PMID:25605798, PMID:28860464, PMID:23307570, PMID:26455413). "
                "Because hypoxia/cytokines upregulate XIAP (BIRC4), which inhibits processed Caspase-9 via its BIR3 domain (Srinivasula et al. 2001 PMID:11242052), "
                "Decision Experiment #4 evaluates an unvalidated preclinical salvage hypothesis of co-dosing AP1903 with an IAP antagonist "
                "(IV birinapant [PMC4742409] or oral xevinapant [noting Phase III TrilynX termination, doi:10.1200/JCO-25-00272]) (Findings 081, 086)."
            ),
            "step1_step2_sequential_editing_and_biallelic_hdr_yield_note": (
                "Sequential 2-step editing (<=2 loci/step) with orthogonal nucleases (HiFi-SpCas9 + AsCas12a Ultra) and transient i53/RAD51 mRNA "
                "targets ~25-40% single-locus biallelic HDR (~6-16% dual-locus biallelic clone yield per step across 96-384 screened clones; "
                "Schiroli et al. 2019 PMID:30905619; Kosicki et al. 2018 PMID:30010673; Turchiano et al. 2021 PMID:33626327). "
                "In Step 2, targeting Locus 3 (CXCL10, Chr4q21) and Locus 4 (AAVS1, Chr19q13.42) yields N_effective = 2 simultaneously cleaved "
                "chromosomal loci (1 heterotypic locus pair), or N_total = 4 diploid alleles (6 allele pairs: 4 inter-chromosomal + 2 homologous sister/homolog pairs) (Findings 028, 082)."
            ),
        },
        "edit_stack_specification": [
            {
                "locus": 1,
                "edit": "B2M-ISRE::[EEF1A1-core + chimeric intron]-HLA-E*01:03(VMAPRTLVL)-GSG-P2A-CD47-IgV-GPI-WPREmut6 (or Stage-0 Pre-Pilot Scarless 2A-Free Bidirectional Arm 1B; Endogenous B2M exon 1 knock-in)",
                "role": "Eliminates classical HLA-A/B/C in 1 locus while resolving the C-terminal GPI-signal vs 2A-scar conflict and explicitly accounting for 2A cytosolic tail topology (Audit F-07; Red-Team RT2-02/RT2-03 & RT1-01; Audit R01/R17): for the GPI-anchored construct (Xu et al. 2024 PMID:38508139), places HLA-E*01:03(VMAPRTLVL) 5'-upstream of GSG-P2A (Liu et al. 2017 PMID:28526819) and places CD47-IgV-GPI at the 3' cistron so the 37-aa DAF (CD55) GPI-anchoring signal sequence resides at the true extreme C-terminus without a 21-aa P2A stub blocking ER GPI-transamidase (PIGK) omega-site cleavage (and ER Signal Peptidase SEC11A/C cleaves the N-terminal Pro-SP of CD47-IgV-GPI). Note on membrane topology (Audit R01): because HLA-E*01:03 single-chain trimer is a Type I transmembrane protein whose C-terminal cytoplasmic tail resides in the cytosol (pH ~7.2), whereas Furin (PCSK3) is a trans-Golgi lumenal protease (pH ~6.0, Ca2+-dependent), Furin is omitted and in Arm 1A HLA-E*01:03 retains the 21-aa cytosolic -GSG-P2A peptide tail. To test a completely scarless, read-through-free architecture where neither protein carries a 2A tail, the Stage-0 Construct Pre-Pilot (Decision Experiment #1; GPI-PLC sensitivity + QuantiBRITE flow + anti-2A/anti-HLA-E Western blot) compares (a) Arm 1A bicistronic (4,880 bp donor / 5,170 bp with ITRs, non-viral linear dsDNA/ceDNA HDR per Shy et al. 2023 PMID:36008610 supports a different nonviral editing context; exact ceDNA/Arm 1A delivery in hPSCs remains unvalidated; note that bacterial IS110 Bridge-RNA recombinases in Durrant et al. 2024 PMID:38926615 and Pelea et al. 2026 PMID:41642947 are exploratory and not validated in human iPSCs), (b) scarless 2A-free bidirectional Arm 1B [SV40pA-CD47-IgV-GPI(antisense) <- EEF1A1-core | Syn_5xISRE_minP -> HLA-E*01:03(VMAPRTLVL)-bGHpA] (3,974 bp donor / 4,264 bp with ITRs, single-AAV6 compatible <= 4,700 bp), and (c) transmembrane 5'-CD47-dCT-tPT2A-HLA-E and 5'-CD47-WT(OA3-323)-tPT2A-HLA-E (clinical UP421 benchmark PMID:37156915, PMID:40757665; beta-cell Cdc42/GSIS counter-evidence in Ghimire et al. 2023 PMID:37820008; ~12% VMAPRTLVL overlap with clinical HCMV UL40 isolates in Lin et al. 2023 PMID:37264229; optional cytoplasmic KW-motif per He et al. 2025 PMC12582296).",
                "supporting_pmids": ["37156915", "28504668", "38508139", "33443542", "28526819", "36008610"],
            },
            {
                "locus": 2,
                "edit": "CIITA-/- (Biallelic Knockout)",
                "role": "Eliminates donor HLA-DR/DP/DQ Class II direct CD4+ presentation (gnomAD pLI = 2.0e-14). Paired in Gate 5 with transient biodegradable SA-FasL + SA-THBD microgels (Headen et al. 2018 PMID:29867165; Lei et al. 2022 PMID:35559682; Turan et al. 2024 PMID:38709159; Kazmi et al. 2026 PMID:41842939 — noting that murine and NHP SA-FasL PEG microgels [Headen et al. 2018 PMID:29867165; Lei et al. 2022 PMID:35559682] used transient short-course rapamycin, whereas Turan et al. 2024 [PMID:38709159] tested direct surface SA-THBD + SA-CD47 in a murine syngeneic intraportal IBMIR/engraftment model without rapamycin and Kazmi et al. 2026 [PMID:41842939] tested direct surface SA-THBD + SA-CD47 in mouse BALB/c->C57BL/6 intraportal allografts without SA-FasL, PEG microgels, or systemic immunosuppression, and IL-1beta sensitizes human islet beta-cells to Fas [CD95] apoptosis [Loweth et al. 1998 PMID:9588443], which is gated in Decision Experiment #2 as incremental apoptosis Delta_TUNEL+ <= +5.0 percentage points above cytokine-matched SA-THBD-only control AND anti-Fas ZB4 mAb rescue <= 5.0 pp; L3B-01/02/03 [DEC-077]).",
                "supporting_pmids": ["37156915", "30853558", "29867165", "35559682", "38709159", "41842939"],
            },
            {
                "locus": 3,
                "edit": "CXCL10::SERPINB9-GSG-P2A-PD-L1 (Inverted-Cistron Endogenous CXCL10 Exon 1 Swap-In, 4,090 bp donor / 4,380 bp with ITRs <= 4,700 bp single-AAV6 limit)",
                "role": "Deletes beta-cell CXCL10 chemokine while comparing inverted cistron order SERPINB9-GSG-P2A-PD-L1 against PD-L1-GSG-P2A-SERPINB9 (Red-Team RT2-01, Decision Experiment #3): in the inverted topology, leaderless cytosolic SERPINB9 (PI-9) is translated first on free cytosolic ribosomes and released into the cytosol upon P2A skipping (Liu et al. 2017 PMID:28526819), after which the nascent PD-L1 ER signal peptide (aa 1-18) engages SRP/Sec61 and is cleaved by ER Signal Peptidase (SEC11A/C) after Ala18 to yield mature PD-L1 (aa 19-290) with its native C-terminal cytoplasmic tail (preserving Cys272 S-palmitoylation by DHHC3 [Yao et al. 2019 PMID:30952982] and CMTM6 stabilization [Burr et al. 2017 PMID:28813417], while empirically testing in Decision Experiment #3 whether the 21-aa P2A C-terminal tail on SERPINB9 affects PI-9 Granzyme B inhibition kinetics; Sintov et al. 2022 PMID:36055241; Tanaka et al. 2009 PMID:19641142; Colli et al. 2018 PMID:30269996).",
                "supporting_pmids": ["36055241", "19641142", "32814902", "30269996", "28526819", "30952982", "28813417"],
            },
            {
                "locus": 4,
                "edit": "AAVS1-[455bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE)]-EEF1A1-core-iCasp9 Biallelic Safe-Harbor Suicide Switch (3,457 bp donor / 3,747 bp with ITRs <= 4,700 bp single-AAV6 limit)",
                "role": "Pan-lineage inducible Caspase-9 safety switch insulated with the 455 bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE alias; insulation function pending) + 250 bp cHS4 chromatin boundary (Muller-Kuller et al. 2015 PMID:25605798; Zhang et al. 2017 PMID:28860464; Pfaff et al. 2013 PMID:23307570) driving the 230 bp human EEF1A1 core promoter with the unvalidated aim of reducing Stage 0->7 epigenetic silencing without bidirectional perturbation of PPP1R12C (Ordovas et al. 2015 PMID:26455413; Wunderlich et al. 2022 PMID:35795779; Yagyu et al. 2016 PMID:27626039). Exact AAV6 packaging ledger (Audit R17 / v0.8.0): 2x145 bp ITRs (290 bp) + 2x600 bp AAVS1 HAs (1,200 bp) + 250 bp cHS4 core insulator + 455 bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE alias; insulation function pending) + 230 bp EEF1A1 core promoter + 9 bp Kozak + 1,191 bp FKBP12(F36V)-DeltaCasp9 + 122 bp SV40 late polyA = 3,457 bp donor / 3,747 bp total with ITRs (+953 bp headroom below the 4,700 bp single-AAV6 packaging limit). Evaluates an unvalidated preclinical salvage hypothesis of co-dosing AP1903 (0.4 mg/kg) with a SMAC mimetic / IAP antagonist (IV birinapant [PMC4742409] or oral xevinapant [noting Phase III TrilynX termination, doi:10.1200/JCO-25-00272]) to displace XIAP (BIRC4) BIR3-domain inhibition of dimerized Caspase-9 in hypoxic/cytokine-stressed grafts (Srinivasula et al. 2001 PMID:11242052; Decision Experiment #4).",
                "supporting_pmids": ["35795779", "27626039", "36599351", "25605798", "23307570", "28860464", "26455413", "11242052"],
            },
        ],
        "counter_evidence_and_residual_risks": [
            "Hour-0 Basal B2M Promoter Density & NKG2A- KIR+ NK Subset Gap (Audit §3.1; Red-Team RT1-01/RT2-02): Resting Stage 7 beta-cells have low basal B2M transcription and ~35-45% of human CD56dim NK cells are NKG2A- KIR2DL1/2/3+; addressed by comparing Arm 1A bicistronic B2M-ISRE::[EEF1A1-core+intron]-HLA-E*01:03(VMAPRTLVL)-GSG-P2A-CD47-IgV-GPI-WPREmut6 and Arm 1B scarless 2A-free bidirectional cassette against 5'-CD47(dCT vs WT)-tPT2A-HLA-E and verifying >200,000 CD47 molecules/cell by QuantiBRITE flow on FACS-gated NKG2A- KIR+ NK subsets and HCMV+ NKG2C+ adaptive NK subsets (Decision Experiment #1).",
            "Indirect Allo/Auto-Antigen Presentation, FasL Extrinsic Death & CaO2/Packing Physics (Audit §3.4; Red-Team RT1-02/RT4-01/RT4-02; AUDIT-REF; Findings 051, 056; L3B-01/02/03 [DEC-077/DEC-078]): Addressed without extra genomic cuts by testing transient biodegradable SA-FasL + SA-THBD microgels (SA-FasL PEG microgels [Headen et al. 2018 PMID:29867165; Lei et al. 2022 PMID:35559682, both with transient short-course rapamycin] + SA-THBD [direct islet-surface SA-THBD + SA-CD47 without rapamycin in Turan et al. 2024 PMID:38709159 (syngeneic) and Kazmi et al. 2026 PMID:41842939 (allogeneic)]; gated on cytokine-induced beta-cell Fas fratricide PMID:9588443), hydrophobic PDMS-CaO2 (OxySite) + covalently immobilized Catalase (100-200 U/mg CaO2 = 79,150-158,300 U/pouch input activity, [H2O2] < 5 uM, plus 72-h 5.710 mmol Ca2+ / 11.420 mmol OH- export and CaCO3/Ca(OH)2 pore-patency release gates; Pedraza et al. 2012 PMID:22371586, Montazeri et al. 2016 PMID:26970510, Lenzen et al. 1996 PMID:8720919; R3-SCI-08, R3-SCI-09), isogenic 3:1 iEC:iPericyte iMVFs (de-novo hPSC iMVFs [Wimmer et al. 2019 PMID:30651639] benchmarked against primary adipose MVFs [Aghazadeh et al. 2021 PMID:34480863]; R3-SCI-10), and macro-scale oxygen packing (single-cluster Stage 7 human SC-islet Krogh radius at Delta_pO2 = 15 mmHg and V_max = 0.034 mol/m^3/s is R_crit = 74.6 um [149.2 um diameter] in free saline [alpha*D_free = 2.1e-14 mol/cm/s/mmHg; vs R_crit = 57.8 um / 115.6 um diameter in bulk hydrogel at alpha*D_gel = 1.26e-14 mol/cm/s/mmHg] and R_crit = 67.9 um [135.8 um diameter] under conservative intra-islet diffusivity [alpha*D_islet = 1.742e-14 mol/cm/s/mmHg in Layer 4, where zero-order drop at R = 60 um is 11.71 mmHg], or 92.6 um [185.2 um diameter] at Delta_pO2 = 23 mmHg in free saline vs rat islet 74.6 um at V_max = 0.0524 mol/m^3/s; passive slab <=1,198 IEQ/cm^2 at 2L_slab = 350 um, or H_crit_fixed_sigma <= 139.8 um at fixed 3,000 IEQ/cm^2 [2L_crit_fixed_phi <= 221.2 um at fixed volume fraction phi = 0.1515], or active 72h PDMS-CaO2 inventory of 0.7915 g [791.5 mg / 10.981 mmol] total loaded pure CaO2 [0.4116 g / 52.0% hydrolyzed over 72 h at t_1/2 = 68 h; embedded in 3.166 g / 2.533 g / 1.706 g total PDMS-CaO2 composite mass for Layouts A / B / Design X; delivering 0.4116 g reacted CaO2 at f_72h = 0.45 requires 0.915 g pure CaO2 / 2.033 g composite at 45 wt%, vs 4.635 g pure CaO2 at uncalibrated bulk disk f_72h = 0.0888 or 8.91 g if scaling the full 0.7915 g inventory by 1/0.0888] across 7.00 cm^3 / 200 cm^2).",
            "Why G6PC2-/- and RNLS-/- are excluded from Chassis A Core: G6PC2-/- lowers fasting glucose and risks fasting/exercise hypoglycemia (Bosma 2020 PMID:32569842), while human RNLS-/- iPSC-beta cells exhibited persistent poor dynamic GSIS and higher OCR in MacDonald et al. 2025 Fig. 6 (PMID:39988068).",
        ],
    },
    {
        "architecture_id": "ISLET1_CHASSIS_B_HLA_C_RETAINED_6LOCUS_HYPOTHESIS",
        "name": "Islet-1 Chassis B ('HLA-C-Retained KIR-Matched Superdonor + TYK2-/- + INS Base-Edit Hypothesis', K=6, B=1)",
        "matched_combinatorial_search_id": (
            "ALLO_HLA_C_RETAINED_SUPERDONOR_PLUS_TYK2_KO__AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT__"
            "CYTO_CXCL10_SWAP_PDL1_SERPINB9__ER_NONE__"
            "NICHE_TRIHORMONAL_120UM_MVF_CAO2__SAF_AAVS1_CAG_ICASP9_BIALLELIC"
        ),
        "genomic_locus_edit_count_K": 6,
        "base_edit_count_B": 1,
        "edit_stack_specification": [
            {
                "locus": "1-4 (HLA-A-/-, HLA-B-/-, CIITA-/-, TYK2-/-)",
                "edit": "Biallelic KO of HLA-A, HLA-B, CIITA, and TYK2 (retaining 1 recipient-matched HLA-C allele)",
                "role": "4 distinct chromosomal KOs. Retains 1 recipient-matched HLA-C allele (Xu et al. 2019 PMID:30853558) to engage KIR2DL1/2/3 (C1/C2 epitope groups) while TYK2-/- blocks IFN-alpha-induced MHC-I hyper-expression (Chandra et al. 2022 PMID:36289205). Note on leader-peptide genetics: the canonical position -21M vs -21T dimorphism is primarily an HLA-B polymorphism whereas HLA-A and HLA-C typically encode -21M leader peptides (Ramsuran et al. 2018 PMID:29302013), but because HLA-C has ~10-fold lower basal surface expression than HLA-A/B, deleting HLA-A and HLA-B requires empirical verification that the retained HLA-C allele alone stabilizes sufficient surface HLA-E for CD94/NKG2A+ NK inhibition (Lin et al. 2023 PMID:37264229).",
                "supporting_pmids": ["30853558", "37264229", "36289205", "34407395"],
            },
            {
                "locus": 5,
                "edit": "CXCL10-ISRE/NFkB::[hPGK-min]-SERPINB9-[internal-promoter-or-eIF4G-leader]-PD-L1 (2A-Free Inverted Swap-In for HLA-C-Retained Arm)",
                "role": "Endogenous inflammation-inducible + minimal basal hPGK floor driving SERPINB9 and PD-L1 without viral 2A peptides (Red-Team RT2-04: prevents proteasome/TAP processing of foreign porcine viral 2A peptide neo-epitopes onto retained HLA-C and compensates for TYK2-/- blunting of early IFN-alpha CXCL10 induction; note that because the human eIF4G 5'-leader acts predominantly via cryptic promoter activity rather than a canonical viral IRES in mammalian bicistronic assays [Han & Zhang 2002 PMID:12370285], Northern blot / 5'-RACE transcript integrity or an explicit dual-promoter architecture is required).",
                "supporting_pmids": ["36055241", "19641142", "32814902", "30269996", "12370285"],
            },
            {
                "locus": 6,
                "edit": "Biallelic AAVS1-[455bp HNRNPA2B1-proximal candidate (legacy CBX3-A2UCOE)]-EEF1A1-core-iCasp9 Safe-Harbor Switch (+ 1 INS Leader/DRiP Base Edit, B=1)",
                "role": "Pan-lineage A2UCOE-insulated EEF1A1-core-iCasp9 suicide switch (synchronized with Locus 4 of Chassis A, 3,457 bp donor / 3,747 bp with ITRs; PMID:35795779, PMID:25605798, PMID:23307570, PMID:26455413) + INS leader/DRiP base edit to reduce direct presentation of insulin epitopes on the retained HLA-C allele (PMID:42690930, PMID:28263308).",
                "supporting_pmids": ["35795779", "25605798", "23307570", "26455413", "42690930", "28263308"],
            },
        ],
        "counter_evidence_and_residual_risks": [
            "Higher Multiplex Edit Burden (K=6 DSB loci + B=1 base-edit locus): assigned engineering utility discount 0.7149 versus 1.0000 for K=4 Chassis A; physical translocation rates and clonal yield are unmeasured.",
            "Constitutive TYK2-/- Stage 0->7 Endocrine Yield Loss (Chandra et al. 2022 PMID:36289205 Fig. 2i-j): Constitutive TYK2 knockout compromises Stage 5 endocrine precursor emergence via KRAS dysregulation and reduces Stage 7 monohormonal INS+ beta-cell numbers by ~32% (26.9% vs 41.8% INS+ cells in Fig. 6). Stage-specific post-differentiation editing or inducible TYK2 degradation/allosteric TYK2i (BMS-986165) must be compared against constitutive KO.",
            "HLA-C vs HLA-B Signal Peptide & Expression Context (Ramsuran et al. 2018 PMID:29302013; Lin et al. 2023 PMID:37264229): Whereas HLA-B exhibits the canonical -21M/-21T dimorphism and HLA-C typically encodes -21M, deleting HLA-A and HLA-B removes the dominant HLA-A/-21M HLA-B peptide supply and leaves endogenous HLA-E surface stabilization dependent on lower-expression HLA-C and KIR2DL C1/C2 group matching.",
            "Direct & Indirect Presentation of Non-Insulin Autoantigens: Retained HLA-C can still present GAD65, ZnT8, IGRP, and IA-2 peptides if IFN-gamma/JAK1/JAK2 bypasses TYK2-/-.",
        ],
    },
    {
        "architecture_id": "BENCHMARK_1_SANA_UP421_HU2024_HIP",
        "name": "Clinical Benchmark 1: Sana UP421 / Hu 2024 HIP (B2M-/- + CIITA-/- + Lentiviral CD47 OE, K=3)",
        "matched_combinatorial_search_id": (
            "ALLO_B2M_CIITA_KO_CD47_LENTI__AUTO_NONE__CYTO_NONE__"
            "ER_NONE__NICHE_NAKED_AVASCULAR_TRANSPLANT__SAF_NONE"
        ),
        "genomic_locus_edit_count_K": 3,
        "crispr_dsb_locus_count": 2,
        "lentiviral_integration_locus_count": 1,
        "edit_accounting_note": (
            "Counts 2 targeted CRISPR nuclease DSB knockouts (B2M-/-, CIITA-/-) + 1 semi-random lentiviral CD47 "
            "integration as K=3 total genomic modification loci in the uniform search grid; note that under the "
            "K<=4 DSB translocation penalty threshold exp(-lambda_K * max(0, K - 4)), both K=2 and K=3 yield "
            "manufacturing_locus_penalty_factor = 1.0000, whereas lentiviral insertional mutagenesis risk is tracked "
            "separately in Gate 5 safety (Finding 025)."
        ),
        "base_edit_count_B": 0,
        "edit_stack_specification": [],
        "counter_evidence_and_residual_risks": [
            "Sole Tier-1 Human Clinical Proof-of-Concept (PMID:40757665; 14-month correspondence follow-up in Carlsson et al. 2026 PMID:42430558, indexed as a zero-weight Project_Summary): Demonstrated C-peptide survival in n=1 human T1D intramuscular forearm sentinel implant (sub-curative sentinel dose) using primary cadaveric islets.",
            "Lacks intrinsic cytokine/granzyme shield (PMID:31900242), pre-vascularized oxygen bridge, or pan-lineage suicide switch required for a 0.8B-cell hPSC-derived pediatric curative dose.",
        ],
    },
    {
        "architecture_id": "BENCHMARK_2_CRISPR_VIACYTE_VCTX210",
        "name": "Clinical Benchmark 2: CRISPR/ViaCyte VCTX210 (B2M-/- + PD-L1 KI + TNFAIP3/A20 KI, K=3)",
        "matched_combinatorial_search_id": (
            "ALLO_B2M_KO_PDL1_A20_KI__AUTO_NONE__CYTO_NONE__"
            "ER_NONE__NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION__SAF_NONE"
        ),
        "genomic_locus_edit_count_K": 3,
        "base_edit_count_B": 0,
        "edit_stack_specification": [],
        "counter_evidence_and_residual_risks": [
            "Lacks CIITA-/- (leaving HLA Class II intact) and lacks CD47 or HLA-E (leaving B2M-/- cells exposed to NK missing-self lysis), delivered inside a macro-encapsulation device.",
        ],
    },
    {
        "architecture_id": "BENCHMARK_3_UNEDITED_MACROENCAPSULATED_VX264",
        "name": "Clinical Benchmark 3: Unedited SC-Islets (Vertex VX-264 Macrodevice / zimislecel VX-880 with Systemic Immunosuppression, K=0)",
        "matched_combinatorial_search_id": (
            "ALLO_NONE_UNEDITED_OR_SYSTEMIC_IMMUNOSUPPRESSION__AUTO_NONE__CYTO_NONE__"
            "ER_NONE__NICHE_UNEDITED_MACRODEVICE_ENCAPSULATION__SAF_NONE"
        ),
        "genomic_locus_edit_count_K": 0,
        "base_edit_count_B": 0,
        "clinical_modality_separation_note": (
            "Separates two distinct unedited SC-islet clinical modalities (Finding 025): "
            "(1) Vertex VX-264 (NCT05791201): unedited SC-islets inside a channel-array macro-encapsulation device without systemic immunosuppression "
            "(modeled directly by this search-grid combination; clinical program discontinued March 2025 due to insufficient C-peptide across the diffusion barrier); "
            "(2) Vertex zimislecel / VX-880 (NCT04786262 / NCT06832410; Reichman et al. 2025 NEJM): naked intraportal infusion of unedited SC-islets "
            "achieving insulin independence in 10/12 adult T1D recipients at Day 365, but requiring lifelong systemic calcineurin/mTOR/anti-thymocyte immunosuppression "
            "that is contraindicated for the immunosuppression-free pediatric Target Product Profile."
        ),
        "edit_stack_specification": [],
        "counter_evidence_and_residual_risks": [
            "0 genetic edits (K=0); relies on chronic systemic immunosuppression (zimislecel / VX-880, NCT04786262 / NCT06832410) or a cell-impermeable diffusion membrane (VX-264, NCT05791201, >150 um) that causes central hypoxic necrosis and fibrotic foreign-body reaction.",
        ],
    },
    {
        "architecture_id": "FLAWED_STACK_CLONAL_STAT1_OR_TYK2_KO_WITH_INDUCIBLE_CHECKPOINTS",
        "name": "Flawed Epistatic Stack: Clonal STAT1-/- Combined with GAS/ISRE-Inducible Checkpoints + G6PC2::iCasp9 (K=4)",
        "matched_combinatorial_search_id": (
            "ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO__AUTO_NONE__"
            "CYTO_STAT1_OR_TYK2_CLONAL_KO_WITH_GAS_IND_PDL1__ER_NONE__"
            "NICHE_NAKED_AVASCULAR_TRANSPLANT__SAF_G6PC2_ENDOGENOUS_ICASP9"
        ),
        "genomic_locus_edit_count_K": 4,
        "base_edit_count_B": 0,
        "edit_stack_specification": [],
        "counter_evidence_and_residual_risks": [
            "Triggers Graveyard 001 (BLOCK): Clonal STAT1 knockout abolishes IFN-gamma/STAT1-mediated induction of GAS/ISRE promoters.",
            "Triggers Graveyard 002 (BLOCK): Lineage-restricted G6PC2::iCasp9 fails to ablate dedifferentiated, teratoma, or SLC18A1+ off-target cells.",
        ],
    },
]


# ---------------------------------------------------------------------------------------------------------------------------------
# Exact-rational algebra for the HIP-Core wet-lab design.  Every multiplier / leverage / VIF reported in the protocol JSON is computed
# here from the run table with `fractions.Fraction` (no floating-point arithmetic and no hand-typed numbers), and is re-derived
# independently by tests/test_design_algebra_exact.py.
# ---------------------------------------------------------------------------------------------------------------------------------
HIP_FACTOR_KEYS = ("X1", "X2", "X3", "X4", "X5", "X6")
HIP_2FI_PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]
HIP_TERM_NAMES = (
    ["I"]
    + list(HIP_FACTOR_KEYS)
    + [f"{HIP_FACTOR_KEYS[i]}:{HIP_FACTOR_KEYS[j]}" for i, j in HIP_2FI_PAIRS]
    + ["stage_block_B"]
)
# Current Scenario A lead genotype: X1 = +1 HLA-E single-chain-trimer knock-in, X2 = +1 constitutive hEF1a-core CD47 package (3' of GSG-P2A in Arm 1A or antisense 2A-free in Arm 1B), X3 = +1 CXCL10 Locus-3 swap;
# TYK2 / TXNIP / RNLS wild type (X4..X6 = -1).
HIP_LEAD_CHASSIS_A_CODED = (1, 1, 1, -1, -1, -1)
HIP_LEAD_CHASSIS_BLOCK = "Tier2_Lead_Chassis_A_Confirmatory_Replicate"


def _hip_frac_text(x: "Fraction") -> str:
    return f"{x.numerator}/{x.denominator} ({float(x):.6f})" if x.denominator != 1 else f"{x.numerator} ({float(x):.6f})"


def _hip_frac_plain(x: "Fraction") -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def _hip_model_rows(runs: list[dict], include_block: bool = True) -> list[list[int]]:
    if not isinstance(runs, list) or not runs:
        raise ValueError("HIP design requires nonempty run rows")
    rows = []
    for run in runs:
        if not isinstance(run, dict) or not isinstance(run.get("coded_levels"), dict):
            raise ValueError("Each HIP run requires an explicit coded-level mapping")
        levels = [run["coded_levels"].get(k) for k in HIP_FACTOR_KEYS]
        if any(type(level) is not int or level not in (-1, 1) for level in levels):
            raise ValueError("HIP coded factors must be exact integer -1/+1 values")
        if include_block and (type(run.get("stage_block")) is not int or run["stage_block"] not in (0, 1)):
            raise ValueError("HIP stage block must be exact integer 0/1")
        row = [1] + levels + [levels[i] * levels[j] for i, j in HIP_2FI_PAIRS]
        if include_block:
            row.append(int(run["stage_block"]))
        rows.append(row)
    return rows


def _hip_gram(rows: list[list[int]]) -> list[list["Fraction"]]:
    p = len(rows[0])
    gram = [[0] * p for _ in range(p)]
    for row in rows:
        for a in range(p):
            if row[a]:
                for b in range(p):
                    gram[a][b] += row[a] * row[b]
    return [[Fraction(v) for v in r] for r in gram]


def _hip_inverse(matrix: list[list["Fraction"]]) -> list[list["Fraction"]]:
    n = len(matrix)
    aug = [list(row) + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if aug[i][col] != 0)  # StopIteration => singular; callers check rank first
        aug[col], aug[pivot] = aug[pivot], aug[col]
        inv_p = 1 / aug[col][col]
        aug[col] = [v * inv_p for v in aug[col]]
        for i in range(n):
            if i != col and aug[i][col] != 0:
                f = aug[i][col]
                aug[i] = [u - f * v for u, v in zip(aug[i], aug[col])]
    return [row[n:] for row in aug]


def _hip_pivot_columns(rows: list[list[int]]) -> list[int]:
    a = [[Fraction(v) for v in r] for r in rows]
    pivots: list[int] = []
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c] != 0), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        inv_p = 1 / a[r][c]
        a[r] = [v * inv_p for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c] != 0:
                f = a[i][c]
                a[i] = [u - f * v for u, v in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == len(a):
            break
    return pivots


def _hip_hat_diag(rows: list[list[int]]) -> list["Fraction"]:
    """Diagonal of the column-space projector (valid for rank-deficient matrices through their pivot columns)."""
    pivots = _hip_pivot_columns(rows)
    sub = [[r[c] for c in pivots] for r in rows]
    inv = _hip_inverse(_hip_gram(sub))
    out = []
    for r in sub:
        left = [sum(r[a] * inv[a][b] for a in range(len(r)) if r[a]) for b in range(len(r))]
        out.append(sum(left[b] * r[b] for b in range(len(r))))
    return out


def _hip_design_summary(runs: list[dict]) -> dict:
    """Exact-rational summary of the 22-column genetic model [I, X1..X6, 15 2FIs] and the 23-column model with additive stage block B."""
    n = len(runs)
    x23 = _hip_model_rows(runs, True)
    x22 = [row[:22] for row in x23]
    piv22, piv23 = _hip_pivot_columns(x22), _hip_pivot_columns(x23)
    if len(piv22) != 22 or len(piv23) != 23:
        raise ValueError("Design is not full rank; algebra below would be undefined")
    inv22, inv23 = _hip_inverse(_hip_gram(x22)), _hip_inverse(_hip_gram(x23))

    def centered_vif(rows: list[list[int]], inv: list[list["Fraction"]]) -> list["Fraction"]:
        out = []
        for j in range(len(rows[0])):
            col = [Fraction(r[j]) for r in rows]
            mean = sum(col) / n
            out.append(inv[j][j] * sum((v - mean) ** 2 for v in col))
        return out

    b = [r[22] for r in x23]
    xtb = [sum(x22[i][j] * b[i] for i in range(n)) for j in range(22)]
    ssr = Fraction(sum(v * v for v in b)) - sum(xtb[a] * inv22[a][c] * xtb[c] for a in range(22) for c in range(22))
    mean_b = Fraction(sum(b), n)
    cross = [sum(Fraction(x23[i][j]) * (b[i] - mean_b) for i in range(n)) for j in range(1, 22)]
    tier1_rows = [row[:22] for row, run in zip(x23, runs) if run["block"].startswith("Tier1")]
    hat22 = _hip_hat_diag(x22)
    hat23 = _hip_hat_diag(x23)
    loo_max_main_23 = Fraction(0, 1)
    loo_max_2fi_23 = Fraction(0, 1)
    loo_max_b_23 = Fraction(0, 1)
    for i_row, r_i in enumerate(x23):
        h_i = hat23[i_row]
        if h_i < 1:
            denom = Fraction(1, 1) - h_i
            v_i = [sum(inv23[a][c] * r_i[c] for c in range(23) if r_i[c]) for a in range(23)]
            diag_loo = [inv23[a][a] + (v_i[a] * v_i[a]) / denom for a in range(23)]
            loo_max_main_23 = max(loo_max_main_23, max(diag_loo[1:7]))
            loo_max_2fi_23 = max(loo_max_2fi_23, max(diag_loo[7:22]))
            loo_max_b_23 = max(loo_max_b_23, diag_loo[22])
    loo_min_ssr_b_23 = Fraction(1, 1) / loo_max_b_23 if loo_max_b_23 > 0 else Fraction(0, 1)
    return {
        "n": n,
        "rank22": len(piv22),
        "rank23": len(piv23),
        "rank_tier1_22": len(_hip_pivot_columns(tier1_rows)),
        "rank29_with_stage_by_main": len(_hip_pivot_columns([row + [row[22] * row[k] for k in range(1, 7)] for row in x23])),
        "inv22": inv22,
        "inv23": inv23,
        "diag22": [inv22[i][i] for i in range(22)],
        "diag23": [inv23[i][i] for i in range(23)],
        "vif23": centered_vif(x23, inv23),
        "hat22": hat22,
        "hat23": hat23,
        "ssr": ssr,
        "mean_b": mean_b,
        "cross": cross,
        "loo_max_main_23": loo_max_main_23,
        "loo_max_2fi_23": loo_max_2fi_23,
        "loo_max_b_23": loo_max_b_23,
        "loo_min_ssr_b_23": loo_min_ssr_b_23,
        # c = (X'X)^-1 X'B: least-squares coefficients of B on the 22 genetic columns; Var_23(b_k) = Var_22(b_k) + c_k^2 / SSR
        "proj_coef": [sum(inv22[a][d] * xtb[d] for d in range(22)) for a in range(22)],
    }


def _hip_prediction_variance(inv23: list[list["Fraction"]], levels: tuple, stage: int) -> "Fraction":
    row = [1] + list(levels) + [levels[i] * levels[j] for i, j in HIP_2FI_PAIRS] + [stage]
    return sum(row[a] * inv23[a][c] * row[c] for a in range(23) for c in range(23))


def _hip_leverage_groups(run_numbers: list[int], hat: list["Fraction"]) -> list[dict]:
    groups: dict = {}
    for rn, h in zip(run_numbers, hat):
        groups.setdefault(h, []).append(rn)
    return [
        {"h_ii": _hip_frac_text(h), "h_ii_rational": _hip_frac_plain(h), "h_ii_float": round(float(h), 6), "n_runs": len(rns), "runs": rns}
        for h, rns in sorted(groups.items(), key=lambda kv: -kv[0])
    ]


def _hip_main_effect_block(summary: dict) -> dict:
    out = {}
    for k, key in enumerate(HIP_FACTOR_KEYS, start=1):
        v22, v23 = summary["diag22"][k], summary["diag23"][k]
        out[key] = {
            "var_22col": _hip_frac_text(v22),
            "var_23col": _hip_frac_text(v23),
            "block_inflation_pct": round(float((v23 / v22 - 1) * 100), 4),
            "centered_vif_23col": round(float(summary["vif23"][k]), 6),
        }
    return out


def _hip_protocol_algebra_fields(runs: list[dict]) -> dict:
    """Return every design-algebra field of the HIP-Core protocol dict, computed from `runs` (exact rationals)."""
    n = len(runs)
    run_numbers = [r["run_number"] for r in runs]
    s = _hip_design_summary(runs)
    base_runs = [r for r in runs if r["block"] != HIP_LEAD_CHASSIS_BLOCK]
    s_base = _hip_design_summary(base_runs)
    extra = []
    for rn in (3, 14):
        src = next(r for r in runs if r["run_number"] == rn)
        extra.append({"coded_levels": src["coded_levels"], "stage_block": 1, "block": "Tier2_InterStage_Bridge_Control_Optional", "run_number": 1000 + rn})
    s_up = _hip_design_summary(runs + extra)
    extra_3_13 = [
        {"coded_levels": next(r for r in runs if r["run_number"] == rn)["coded_levels"], "stage_block": 1, "block": "Tier2_InterStage_Bridge_Control_Optional", "run_number": 1100 + rn}
        for rn in (3, 13)
    ]
    s_up_3_13 = _hip_design_summary(runs + extra_3_13)
    extra_5_11 = [
        {"coded_levels": next(r for r in runs if r["run_number"] == rn)["coded_levels"], "stage_block": 1, "block": "Tier2_InterStage_Bridge_Control_Optional", "run_number": 1200 + rn}
        for rn in (5, 11)
    ]
    s_up_5_11 = _hip_design_summary(runs + extra_5_11)
    extra_40 = [
        {"coded_levels": next(r for r in runs if r["run_number"] == rn)["coded_levels"], "stage_block": 1, "block": "Tier2_InterStage_Bridge_Control_Optional", "run_number": 1250 + rn}
        for rn in (2, 3, 4, 5)
    ]
    s_up_40 = _hip_design_summary(runs + extra_40)
    extra_42 = [
        {"coded_levels": next(r for r in runs if r["run_number"] == rn)["coded_levels"], "stage_block": 1, "block": "Tier2_InterStage_Bridge_Control_Optional", "run_number": 1300 + rn}
        for rn in (2, 3, 5, 9, 14, 15)
    ]
    s_up_42 = _hip_design_summary(runs + extra_42)
    n_geno = len({tuple(int(r["coded_levels"][k]) for k in HIP_FACTOR_KEYS) for r in runs})
    run_ids = {r["run_number"]: r["run_id"] for r in runs}

    max_h23, max_h22 = max(s["hat23"]), max(s["hat22"])
    unit_rows = sorted({run_ids[rn] for rn, h in zip(run_numbers, s["hat23"]) if h == 1} | {run_ids[rn] for rn, h in zip(run_numbers, s["hat22"]) if h == 1})
    mains = _hip_main_effect_block(s)
    mains_base = _hip_main_effect_block(s_base)
    infl = [m["block_inflation_pct"] for m in mains.values()]
    infl_base = [m["block_inflation_pct"] for m in mains_base.values()]

    # ---- 2FI multiplier tiers (grouped by the exact 23-column multiplier) ----
    tiers: dict = {}
    for idx in range(7, 22):
        tiers.setdefault(s["diag23"][idx], []).append(idx)
    two_fi = {}
    for v, idxs in sorted(tiers.items()):
        pairs = [HIP_TERM_NAMES[i] for i in idxs]
        key = "_and_".join(p.replace(":", "_") for p in pairs) + f"_{len(pairs)}pair" + ("s" if len(pairs) != 1 else "")
        two_fi[key] = {
            "pairs": pairs,
            "var_23col": _hip_frac_text(v),
            "var_22col_by_pair": {HIP_TERM_NAMES[i]: _hip_frac_text(s["diag22"][i]) for i in idxs},
            "centered_vif_23col_by_pair": {HIP_TERM_NAMES[i]: round(float(s["vif23"][i]), 6) for i in idxs},
        }
    max_2fi, min_2fi = max(s["diag23"][7:22]), min(s["diag23"][7:22])
    mean_main_23 = sum(s["diag23"][1:7]) / 6
    mean_2fi_23 = sum(s["diag23"][7:22]) / 15
    if n == 36:
        assert mean_main_23 == Fraction(137, 2424), f"Expected 6-main-effect mean 137/2424, got {mean_main_23}"
        assert mean_2fi_23 == Fraction(2097, 8080), f"Expected 15-2FI mean 2097/8080, got {mean_2fi_23}"
    max_2fi_base = max(s_base["diag23"][7:22])
    max_2fi_up = max(s_up["diag23"][7:22])
    max_vif_up = max(s_up["vif23"][7:22])

    # ---- stage-block cross-products and the unbridged reference design ----
    main_cross: dict = {}
    for k, key in enumerate(HIP_FACTOR_KEYS):
        main_cross.setdefault(s["cross"][k], []).append(key)
    cross_text = ", ".join(f"{_hip_frac_plain(v)} for {','.join(ks)}" for v, ks in sorted(main_cross.items()))
    fi_cross = sorted({_hip_frac_plain(v) for v in s["cross"][6:]}, key=lambda t: Fraction(t))
    unbridged = [r for r in runs if r["run_number"] <= 24]
    x23u = _hip_model_rows(unbridged, True)
    rank22u = len(_hip_pivot_columns([row[:22] for row in x23u]))
    rank23u = len(_hip_pivot_columns(x23u))
    hat_u = _hip_hat_diag(x23u)
    aug_lev = sorted({_hip_frac_plain(h) for rn, h in zip([r["run_number"] for r in unbridged], hat_u) if 17 <= rn <= 24})

    lead_runs = [r["run_number"] for r in runs if tuple(int(r["coded_levels"][k]) for k in HIP_FACTOR_KEYS) == HIP_LEAD_CHASSIS_A_CODED]
    pv0, pv1 = _hip_prediction_variance(s["inv23"], HIP_LEAD_CHASSIS_A_CODED, 0), _hip_prediction_variance(s["inv23"], HIP_LEAD_CHASSIS_A_CODED, 1)
    pv0_base = _hip_prediction_variance(s_base["inv23"], HIP_LEAD_CHASSIS_A_CODED, 0)
    pv1_base = _hip_prediction_variance(s_base["inv23"], HIP_LEAD_CHASSIS_A_CODED, 1)

    base_c_main = s_base["proj_coef"][1:7]
    if all(v == 0 for v in base_c_main):
        base_statement = (
            f"In the base {len(base_runs)}-run design (Runs 01-{len(base_runs)}) the least-squares projection of B on the 22-column genetic model has coefficient exactly 0 on all six main-effect "
            f"columns (it loads only on the intercept, {_hip_frac_plain(s_base['proj_coef'][0])}, and on the 15 two-factor interaction columns), so by the partitioned-inverse identity "
            f"Var_23(b_k) = Var_22(b_k) + c_k^2 / SSR the six main-effect multipliers are identical with and without B ({_hip_frac_plain(s_base['diag23'][1])} in both models): exactly 0.0% inflation."
        )
    else:
        base_statement = f"The base {len(base_runs)}-run design shows {min(infl_base)}-{max(infl_base)}% main-effect block inflation."
    c_groups: dict = {}
    for key, v in zip(HIP_FACTOR_KEYS, s["proj_coef"][1:7]):
        c_groups.setdefault(v, []).append(key)
    c_text = "; ".join(f"{','.join(ks)}: {_hip_frac_plain(v)}" for v, ks in c_groups.items())
    per_factor_text = ", ".join(f"{k} {v['block_inflation_pct']}%" for k, v in mains.items())
    inflation_note = (
        f"Adding the additive stage-block indicator B to the {n}-run design changes the main-effect variance multipliers by {min(infl)}%-{max(infl)}% "
        f"({per_factor_text}); exact per-factor values are given below. "
        + base_statement
        + (f" The {n - len(base_runs)} lead-chassis confirmatory replicates (genotype {list(HIP_LEAD_CHASSIS_A_CODED)}, stage 1 only) move that projection onto the main effects "
           f"(c_k = {c_text}), which is the sole source of the inflation above; "
           f"in exchange the worst-case 2FI multiplier falls from {_hip_frac_text(max_2fi_base)} to {_hip_frac_text(max_2fi)} and the lead genotype's prediction variance "
           f"(stage 0 / stage 1) from {_hip_frac_plain(pv0_base)} / {_hip_frac_plain(pv1_base)} to {_hip_frac_plain(pv0)} / {_hip_frac_plain(pv1)} (single observation = 1)."
           if n != len(base_runs) else "")
    )

    def lev_profile(hat: list["Fraction"], with_block: bool) -> dict:
        maxh = max(hat)
        prof = {
            "max_leverage_h_ii": float(maxh),
            "max_leverage_rational": _hip_frac_plain(maxh),
            "leverage_groups": _hip_leverage_groups(run_numbers, hat),
            "trace_sum_of_leverages": int(sum(hat)),
            "min_rank_after_deleting_any_single_run": (23 if with_block else 22) - (1 if any(h == 1 for h in hat) else 0),
        }
        if with_block:
            prof["stage_block_projection_residual_l2_norm"] = round(float(s["ssr"]) ** 0.5, 6)
            prof["stage_block_projection_residual_ssr"] = _hip_frac_text(s["ssr"])
            prof["single_run_deletion_35run_worst_case_bounds"] = {
                "max_main_effect_var_multiplier_23col": _hip_frac_text(s["loo_max_main_23"]),
                "max_2fi_var_multiplier_23col": _hip_frac_text(s["loo_max_2fi_23"]),
                "max_stage_block_B_var_multiplier_23col": _hip_frac_text(s["loo_max_b_23"]),
                "min_stage_block_projection_residual_ssr": _hip_frac_text(s["loo_min_ssr_b_23"]),
            }
        return prof

    stage_block_b = {
        "var_23col_rational": _hip_frac_text(s["diag23"][22]),
        "centered_vif_23col": round(float(s["vif23"][22]), 6),
        "projection_residual_ssr": _hip_frac_text(s["ssr"]),
        "projection_coefficient_c_of_B_on_intercept": _hip_frac_plain(s["proj_coef"][0]),
        "projection_coefficient_c_of_B_on_main_effects": {key: _hip_frac_plain(v) for key, v in zip(HIP_FACTOR_KEYS, s["proj_coef"][1:7])},
        "rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols": f"{s['rank29_with_stage_by_main']} / 29",
    }
    optional_upgrade = {
        "description": (
            f"Adding 2 further stage-1 bridge controls (repeating Tier-1 Run 03 {list(next(r for r in runs if r['run_number'] == 3)['coded_levels'].values())} and Run 14 "
            f"{list(next(r for r in runs if r['run_number'] == 14)['coded_levels'].values())} alongside Runs 01 and 16) to the {n}-run design gives N = {n + 2} runs with 4 bridge genotypes: "
            f"Var(beta_hat_stage_block) falls from {_hip_frac_plain(s['diag23'][22])} to {_hip_frac_plain(s_up['diag23'][22])} sigma^2 and the worst-case 2FI multiplier from "
            f"{_hip_frac_plain(max_2fi)} to {_hip_frac_plain(max_2fi_up)}."
        ),
        "n_runs_total": n + 2,
        "bridge_runs_repeated_from_tier1": ["HIP_FACTORIAL_RUN_01", "HIP_FACTORIAL_RUN_16", "HIP_FACTORIAL_RUN_03", "HIP_FACTORIAL_RUN_14"],
        "stage_block_ssr": _hip_frac_text(s_up["ssr"]),
        "stage_block_var_multiplier": _hip_frac_text(s_up["diag23"][22]),
        "stage_block_centered_vif": round(float(s_up["vif23"][22]), 6),
        "rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols": f"{s_up['rank29_with_stage_by_main']} / 29",
        "max_2fi_var_multiplier_rational": _hip_frac_plain(max_2fi_up),
        "max_2fi_var_multiplier": round(float(max_2fi_up), 6),
        "max_2fi_centered_vif": round(float(max_vif_up), 6),
        "max_2fi_reduction_vs_this_design_pct": round(float((1 - max_2fi_up / max_2fi) * 100), 2),
        "alternative_4bridge_pair_3_13_metrics": {
            "bridge_runs_repeated_from_tier1": ["HIP_FACTORIAL_RUN_01", "HIP_FACTORIAL_RUN_16", "HIP_FACTORIAL_RUN_03", "HIP_FACTORIAL_RUN_13"],
            "stage_block_ssr": _hip_frac_text(s_up_3_13["ssr"]),
            "stage_block_var_multiplier": _hip_frac_text(s_up_3_13["diag23"][22]),
            "max_2fi_var_multiplier_rational": _hip_frac_plain(max(s_up_3_13["diag23"][7:22])),
            "rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols": f"{s_up_3_13['rank29_with_stage_by_main']} / 29",
        },
        "alternative_4bridge_pair_5_11_metrics": {
            "bridge_runs_repeated_from_tier1": ["HIP_FACTORIAL_RUN_01", "HIP_FACTORIAL_RUN_16", "HIP_FACTORIAL_RUN_05", "HIP_FACTORIAL_RUN_11"],
            "stage_block_ssr": _hip_frac_text(s_up_5_11["ssr"]),
            "stage_block_var_multiplier": _hip_frac_text(s_up_5_11["diag23"][22]),
            "max_2fi_var_multiplier_rational": _hip_frac_plain(max(s_up_5_11["diag23"][7:22])),
            "rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols": f"{s_up_5_11['rank29_with_stage_by_main']} / 29",
        },
        "minimal_full_rank_29col_40run_6bridge_specification": {
            "n_runs_total": n + 4,
            "bridge_runs_repeated_from_tier1": [
                "HIP_FACTORIAL_RUN_01", "HIP_FACTORIAL_RUN_16", "HIP_FACTORIAL_RUN_02",
                "HIP_FACTORIAL_RUN_03", "HIP_FACTORIAL_RUN_04", "HIP_FACTORIAL_RUN_05",
            ],
            "rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols": f"{s_up_40['rank29_with_stage_by_main']} / 29",
            "note": (
                "Because the 36-run design already contains 2 bridge genotypes (Runs 01 and 16) plus 8 Tier-2 de-aliasing genotypes "
                "and 1 Lead-Chassis-A genotype in Stage 1 (giving rank 25/29, deficiency 4), adding just 4 extra Stage-1 bridge runs "
                "(e.g. repeating Tier-1 Runs 02, 03, 04, and 05 -> 6 total bridge genotypes in a 40-run = 16 Tier-1 + 6 bridge + 18 Tier-2 design) "
                "is the minimal augmentation achieving full rank 29/29 on [I, X1..X6, 15 2FIs, B, B*X1..B*X6] (Finding 075)."
            ),
        },
        "full_rank_29col_42run_8bridge_upgrade_requirement": {
            "n_runs_total": n + 6,
            "bridge_runs_repeated_from_tier1": [
                "HIP_FACTORIAL_RUN_01", "HIP_FACTORIAL_RUN_16", "HIP_FACTORIAL_RUN_02", "HIP_FACTORIAL_RUN_03",
                "HIP_FACTORIAL_RUN_05", "HIP_FACTORIAL_RUN_09", "HIP_FACTORIAL_RUN_14", "HIP_FACTORIAL_RUN_15",
            ],
            "rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols": f"{s_up_42['rank29_with_stage_by_main']} / 29",
            "note": (
                "While 4 extra Stage-1 bridge runs (6 bridges total, 40 runs) is the algebraic minimum for full rank 29/29 on the 29-column model "
                "[I, X1..X6, 15 2FIs, B, B*X1..B*X6], the 8-bridge 42-run design (16 Tier-1 + 8 bridge + 18 Tier-2) provides 2 additional "
                "degrees of freedom for inter-stage interaction conditioning (Findings 075, 084)."
            ),
        },
    }
    dealiasing_note = (
        f"In the 16-run 2^(6-2)_IV screening block alone (Runs 01-16), the 22-column second-order model matrix has rank {s['rank_tier1_22']}/22. "
        f"Adding an unbridged 8-run block in a separate batch (Runs 01-24 with B) gives rank {rank22u}/22 genetic but {rank23u}/23 with the stage block, i.e. B is linearly dependent on the 2FIs, "
        f"and every augmentation row has leverage {', '.join(aug_lev)}. Therefore the {n - 16}-run Tier-2 block (Runs 17-{n}) includes 2 inter-stage bridge repeats of Tier-1 Runs 01 and 16 (Runs 25-26), "
        f"duplicated Tier-2 de-aliasing runs (Runs 17-24 and 27-34) and {n - len(base_runs)} lead-chassis confirmatory replicates (Runs {'-'.join(str(x) for x in (min(lead_runs), max(lead_runs)))}), "
        f"achieving full rank {s['rank22']}/22 on genetic effects (max h_ii = {_hip_frac_plain(max_h22)}, trace = {int(sum(s['hat22']))}), full rank {s['rank23']}/23 with the stage-block column B "
        f"(max h_ii = {_hip_frac_plain(max_h23)}, trace = {int(sum(s['hat23']))}; projection residual SSR = {_hip_frac_plain(s['ssr'])}), and {len(unit_rows)} unit-leverage rows "
        "(leverage groups are listed in exact_rational_leverage_profile_with_stage_block_23cols)."
    )
    protocol_fields = {
        "design_name": (
            f"HIP-Core Two-Stage Bridged & Replicated Wet-Lab Protocol: Tier-1 16-Run 2^(6-2) Resolution IV Screening + Tier-2 {n - 16}-Run Bridged, Replicated & Lead-Chassis-Confirmatory "
            f"Multi-Coset Full-Rank (Rank 22/22 & 23/23 with Stage Block, {len(unit_rows)} Unit-Leverage Rows, Max Leverage {_hip_frac_plain(max_h23)} = {float(max_h23):.6f})"
        ),
        "design_class": "2^(6-2)_IV_plus_2Bridge_plus_2x8_Replicated_MultiCoset_plus_2x_LeadChassisA_Confirmatory_FullRank_Augmentation",
        "n_factors": 6,
        "n_runs_tier1_screening": 16,
        "n_runs_tier2_dealiasing": n - 16,
        "n_runs": n,
        "n_distinct_genotypes": n_geno,
        "model_matrix_rank_intercept_6main_15twofi": f"{s['rank22']} / 22 (Full Rank across {n} combined runs; {s['rank_tier1_22']} / 22 in Tier-1 16-run screening block alone)",
        "model_matrix_rank_with_stage_block_indicator_23cols": (
            f"{s['rank23']} / 23 (Full Rank including additive Stage/Batch Block column B; stage block B is linearly independent of and jointly estimable with all 6 main effects and all 15 "
            f"two-factor interactions, with projection residual SSR = ||(I - H_22) B||_2^2 = {_hip_frac_plain(s['ssr'])}; non-orthogonal centered cross-products X^T (B - mean(B)) with "
            f"mean(B) = {_hip_frac_plain(s['mean_b'])} are {cross_text} for the main effects and {{{', '.join(fi_cross)}}} across the 15 2FIs; note that adding all 6 stage_block * X_j "
            f"interaction columns yields rank {s['rank29_with_stage_by_main']} / 29, so the design identifies an additive stage offset rather than all 6 stage-by-factor interactions)"
        ),
        "unit_leverage_augmentation_rows": unit_rows,
        "max_hat_matrix_leverage_across_all_runs": float(max_h23),
        "exact_rational_leverage_profile_with_stage_block_23cols": lev_profile(s["hat23"], True),
        "exact_rational_leverage_profile_genetic_only_22cols": lev_profile(s["hat22"], False),
        "exact_rational_parameter_variance_and_vif_conditioning": {
            "algebraic_note_main_effect_block_inflation": inflation_note,
            "intercept_var_multiplier": {"var_22col": _hip_frac_text(s["diag22"][0]), "var_23col": _hip_frac_text(s["diag23"][0])},
            "main_effects_var_multiplier_by_factor": mains,
            "main_effect_block_inflation_pct_range": [min(infl), max(infl)],
            "two_factor_interactions_var_multiplier_tiers_23col": two_fi,
            "two_factor_interaction_var_multiplier_range_23col": {"min": _hip_frac_text(min_2fi), "max": _hip_frac_text(max_2fi)},
            "stage_block_B": stage_block_b,
            "lead_chassis_a_prediction_variance": {
                "coded_levels": list(HIP_LEAD_CHASSIS_A_CODED),
                "runs_with_this_genotype": [run_ids[rn] for rn in lead_runs],
                "var_multiplier_stage0": _hip_frac_text(pv0),
                "var_multiplier_stage1": _hip_frac_text(pv1),
                "single_observation_reference": 1.0,
                "base_34_run_design_for_comparison": {"var_multiplier_stage0": _hip_frac_text(pv0_base), "var_multiplier_stage1": _hip_frac_text(pv1_base)},
            },
            "base_34_run_design_reference": {
                "n_runs": len(base_runs),
                "main_effect_var_multiplier_23col": _hip_frac_text(s_base["diag23"][1]),
                "main_effect_block_inflation_pct_range": [min(infl_base), max(infl_base)],
                "min_2fi_var_multiplier_23col": _hip_frac_text(min(s_base["diag23"][7:22])),
                "max_2fi_var_multiplier_23col": _hip_frac_text(max_2fi_base),
                "stage_block_var_23col": _hip_frac_text(s_base["diag23"][22]),
                "max_leverage_h_ii_23col": _hip_frac_text(max(s_base["hat23"])),
                "projection_coefficient_c_of_B_on_intercept": _hip_frac_plain(s_base["proj_coef"][0]),
                "projection_coefficient_c_of_B_on_main_effects": {key: _hip_frac_plain(v) for key, v in zip(HIP_FACTOR_KEYS, s_base["proj_coef"][1:7])},
            },
            "optional_38run_4bridge_upgrade_specification": optional_upgrade,
        },
    }
    # ---- Exact power, MDE, and 108-independent-clone vs 75-shared-genotype-clone covariance algebra (Audit R03 & R04) ----
    x23_rows = _hip_model_rows(runs, True)
    x22_rows = [row[:22] for row in x23_rows]
    # Shared-genotype covariance Sigma_shared = (sigma_clone^2 / n_clones) * Z Z^T + (sigma_eps^2 / (n_clones * n_batches)) * I_n
    sigma_clone_sq = Fraction(64, 1)       # 8.0^2 pp^2
    sigma_eps_sq = Fraction(49, 1)         # 7.0^2 pp^2 (donor-averaged clone-by-batch cell residual; if per-donor well residual is 7.0 pp across n_donors = 3, divide by 3 -> 49/18)
    n_clones_f = Fraction(3, 1)
    n_batches_f = Fraction(2, 1)
    n_donors_f = Fraction(3, 1)
    v_clone_run = sigma_clone_sq / n_clones_f                      # 64/3 = 21.333333 pp^2
    v_eps_run = sigma_eps_sq / (n_clones_f * n_batches_f)          # 49/6 = 8.166667 pp^2
    v_eps_run_3donors = sigma_eps_sq / (n_clones_f * n_batches_f * n_donors_f)  # 49/18 = 2.722222 pp^2
    sigma_run_sq = v_clone_run + v_eps_run                         # 177/6 = 29.50 pp^2
    sigma_run_sq_3donors = v_clone_run + v_eps_run_3donors         # 433/18 = 24.055556 pp^2
    sigma_run_pp = math.sqrt(float(sigma_run_sq))                  # 5.431390 pp
    sigma_run_3donors_pp = math.sqrt(float(sigma_run_sq_3donors))  # 4.904646 pp

    geno_tuples = [tuple(int(r["coded_levels"][k]) for k in HIP_FACTOR_KEYS) for r in runs]
    sigma_shared = [
        [
            (v_clone_run if geno_tuples[i] == geno_tuples[j] else Fraction(0, 1))
            + (v_eps_run if i == j else Fraction(0, 1))
            for j in range(n)
        ]
        for i in range(n)
    ]
    sigma_shared_inv = _hip_inverse(sigma_shared)

    def _sandwich_and_gls_diag(x_rows: list[list[int]], inv_xtx: list[list["Fraction"]]) -> tuple[list["Fraction"], list["Fraction"]]:
        p_cols = len(x_rows[0])
        # Meat M = X^T Sigma_shared X
        xt_sig = [
            [sum(Fraction(x_rows[i][a]) * sigma_shared[i][j] for i in range(n)) for j in range(n)]
            for a in range(p_cols)
        ]
        meat = [
            [sum(xt_sig[a][j] * Fraction(x_rows[j][b]) for j in range(n)) for b in range(p_cols)]
            for a in range(p_cols)
        ]
        # Sandwich = inv_xtx * meat * inv_xtx
        left_m = [
            [sum(inv_xtx[a][k] * meat[k][b] for k in range(p_cols)) for b in range(p_cols)]
            for a in range(p_cols)
        ]
        ols_diag = [sum(left_m[a][k] * inv_xtx[k][a] for k in range(p_cols)) for a in range(p_cols)]
        # GLS = (X^T Sigma_shared^-1 X)^-1
        xt_siginv = [
            [sum(Fraction(x_rows[i][a]) * sigma_shared_inv[i][j] for i in range(n)) for j in range(n)]
            for a in range(p_cols)
        ]
        xt_siginv_x = [
            [sum(xt_siginv[a][j] * Fraction(x_rows[j][b]) for j in range(n)) for b in range(p_cols)]
            for a in range(p_cols)
        ]
        gls_cov = _hip_inverse(xt_siginv_x)
        gls_diag = [gls_cov[a][a] for a in range(p_cols)]
        return ols_diag, gls_diag

    ols23_shared_var, gls23_shared_var = _sandwich_and_gls_diag(x23_rows, s["inv23"])
    ols22_shared_var, gls22_shared_var = _sandwich_and_gls_diag(x22_rows, s["inv22"])

    # Central-t quantile sum approximations & exact noncentral-t power multipliers at 80% power:
    # - Normal two-sided alpha=0.05: z_0.975 + z_0.80 = 1.959964 + 0.841621 = 2.801585
    # - Clone-level df=52 two-sided alpha=0.05: central-t sum = 2.855256; exact noncentral-t ncp = 2.854667
    # - Conservative run-table df=13 two-sided alpha=0.05: central-t sum = 3.030521; exact noncentral-t ncp = 3.031292
    # - BH-FDR q=0.10 across 21 genetic terms (rank-1 top discovery alpha = 0.10/21 = 0.004762, two-sided): z_{1 - 0.05/21} + z_0.80 = 3.664382
    mult_norm_a05 = 2.801585
    mult_t52_a05 = 2.855256
    mult_nct52_exact_a05 = 2.854667
    mult_t13_a05 = 3.030521
    mult_nct13_exact_a05 = 3.031292
    mult_bh21_q10 = 3.664382

    min_main_22, max_main_22 = min(s["diag22"][1:7]), max(s["diag22"][1:7])
    min_main_23, max_main_23 = min(s["diag23"][1:7]), max(s["diag23"][1:7])
    min_2fi_22, max_2fi_22 = min(s["diag22"][7:22]), max(s["diag22"][7:22])
    min_2fi_23, max_2fi_23 = min_2fi, max_2fi

    best_2fi_pairs_22 = [HIP_TERM_NAMES[i] for i in range(7, 22) if s["diag22"][i] == min_2fi_22]
    worst_2fi_pairs_22 = [HIP_TERM_NAMES[i] for i in range(7, 22) if s["diag22"][i] == max_2fi_22]
    best_2fi_pairs_23 = [HIP_TERM_NAMES[i] for i in range(7, 22) if s["diag23"][i] == min_2fi_23]
    worst_2fi_pairs_23 = [HIP_TERM_NAMES[i] for i in range(7, 22) if s["diag23"][i] == max_2fi_23]
    n_distinct_with_stage = len({(tuple(int(r["coded_levels"][k]) for k in HIP_FACTOR_KEYS), 0 if r["run_number"] <= 16 else 1) for r in runs})

    def _se_mde_entry(var_frac: "Fraction", sigma_val: float = sigma_run_pp, pairs: list[str] | None = None) -> dict:
        se_half = sigma_val * math.sqrt(float(var_frac))
        se_full_2b = 2.0 * se_half
        se_dod_4b = 4.0 * se_half
        res = {
            "var_multiplier_rational": _hip_frac_text(var_frac),
            "half_effect_beta_se_pp": round(se_half, 2),
            "full_contrast_2beta_se_pp": round(se_full_2b, 2),
            "difference_of_differences_4beta_se_pp": round(se_dod_4b, 2),
            "mde_80pct_power_2beta_normal_alpha_0_05_pp": round(mult_norm_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_t52_central_approx_alpha_0_05_pp": round(mult_t52_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_t52_alpha_0_05_pp": round(mult_t52_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_nct52_exact_alpha_0_05_pp": round(mult_nct52_exact_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_t13_central_approx_alpha_0_05_pp": round(mult_t13_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_t13_alpha_0_05_pp": round(mult_t13_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_nct13_exact_alpha_0_05_pp": round(mult_nct13_exact_a05 * se_full_2b, 2),
            "mde_80pct_power_2beta_bh_fdr_q_0_10_across_21_terms_pp": round(mult_bh21_q10 * se_full_2b, 2),
            "mde_80pct_power_4beta_dod_normal_alpha_0_05_pp": round(mult_norm_a05 * se_dod_4b, 2),
            "mde_80pct_power_4beta_dod_t52_alpha_0_05_pp": round(mult_t52_a05 * se_dod_4b, 2),
            "mde_80pct_power_4beta_dod_nct52_exact_alpha_0_05_pp": round(mult_nct52_exact_a05 * se_dod_4b, 2),
            "mde_80pct_power_4beta_dod_t13_alpha_0_05_pp": round(mult_t13_a05 * se_dod_4b, 2),
            "mde_80pct_power_4beta_dod_nct13_exact_alpha_0_05_pp": round(mult_nct13_exact_a05 * se_dod_4b, 2),
        }
        if pairs is not None:
            res["interaction_pairs"] = pairs
        return res

    # Stage-by-main interaction expansion rank & differentiation-batch random-effect sensitivity (AUDIT-REF):
    x29_stage_by_main_rows = [
        row + [row[22] * row[j] for j in range(1, 7)]
        for row in x23_rows
    ]
    rank_29col_stage_by_main = len(_hip_pivot_columns(x29_stage_by_main_rows))
    # Stage-1 split-batch incidence: Batch S1-A = Runs 17-24, 25, 35 (0-based indices 16..23, 24, 34);
    # Batch S1-B = Runs 27-34, 26, 36 (0-based indices 26..33, 25, 35) (AUDIT-REF):
    s1_batch_a_indices = {16, 17, 18, 19, 20, 21, 22, 23, 24, 34}
    s1_batch_b_indices = {25, 26, 27, 28, 29, 30, 31, 32, 33, 35}
    batch_sensitivity_sweep = {}
    for b_sd_label, b_var_frac in [
        ("batch_sd_0p0_pp", Fraction(0, 1)),
        ("batch_sd_2p5_pp", Fraction(25, 8)),   # (2.5^2) / 2 batches per stage = 6.25 / 2 = 25/8
        ("batch_sd_5p0_pp", Fraction(25, 2)),   # (5.0^2) / 2 batches per stage = 25 / 2 = 12.5
    ]:
        tau_sq_frac = b_var_frac * 2
        sig_b = [
            [
                sigma_shared[i][j]
                + (b_var_frac if x23_rows[i][22] == x23_rows[j][22] else Fraction(0, 1))
                for j in range(n)
            ]
            for i in range(n)
        ]
        sig_b_inv = _hip_inverse(sig_b)
        xt_sbinv = [
            [sum(Fraction(x23_rows[i][a]) * sig_b_inv[i][j] for i in range(n)) for j in range(n)]
            for a in range(23)
        ]
        xt_sbinv_x = [
            [sum(xt_sbinv[a][j] * Fraction(x23_rows[j][b]) for j in range(n)) for b in range(23)]
            for a in range(23)
        ]
        gls_b_cov = _hip_inverse(xt_sbinv_x)
        gls_b_diag = [gls_b_cov[a][a] for a in range(23)]

        # Exact GLS covariance under literal Stage-1 split-batch allocation (Batch S1-A vs Batch S1-B):
        sig_split_2b = [
            [
                sigma_shared[i][j]
                + (
                    b_var_frac
                    if (x23_rows[i][22] == 0 and x23_rows[j][22] == 0)
                    else (
                        tau_sq_frac
                        if ((i in s1_batch_a_indices and j in s1_batch_a_indices) or (i in s1_batch_b_indices and j in s1_batch_b_indices))
                        else Fraction(0, 1)
                    )
                )
                for j in range(n)
            ]
            for i in range(n)
        ]
        sig_split_2b_inv = _hip_inverse(sig_split_2b)
        xt_s2binv_x = [
            [
                sum(
                    Fraction(x23_rows[i][a]) * sig_split_2b_inv[i][j] * Fraction(x23_rows[j][b])
                    for i in range(n)
                    for j in range(n)
                )
                for b in range(23)
            ]
            for a in range(23)
        ]
        gls_split_2b_diag = [_hip_inverse(xt_s2binv_x)[a][a] for a in range(23)]

        # Literal 1-batch-per-Stage-1-row residual (v_eps = 49/6 on Stage-0 rows 0..15, v_eps = 49/3 on Stage-1 rows 16..35):
        sig_split_1b = [
            [
                sig_split_2b[i][j] + (Fraction(49, 6) if (i == j and x23_rows[i][22] == 1) else Fraction(0, 1))
                for j in range(n)
            ]
            for i in range(n)
        ]
        sig_split_1b_inv = _hip_inverse(sig_split_1b)
        xt_s1binv_x = [
            [
                sum(
                    Fraction(x23_rows[i][a]) * sig_split_1b_inv[i][j] * Fraction(x23_rows[j][b])
                    for i in range(n)
                    for j in range(n)
                )
                for b in range(23)
            ]
            for a in range(23)
        ]
        gls_split_1b_diag = [_hip_inverse(xt_s1binv_x)[a][a] for a in range(23)]

        # Mode 1 (108 independent run-nested clone differentiations, diagonal run-mean variance + literal S1-A vs S1-B batch covariance; AUDIT-REF):
        sig_mode1_split_2b = [
            [
                (Fraction(177, 6) if i == j else Fraction(0, 1))
                + (
                    b_var_frac
                    if (x23_rows[i][22] == 0 and x23_rows[j][22] == 0)
                    else (
                        tau_sq_frac
                        if ((i in s1_batch_a_indices and j in s1_batch_a_indices) or (i in s1_batch_b_indices and j in s1_batch_b_indices))
                        else Fraction(0, 1)
                    )
                )
                for j in range(n)
            ]
            for i in range(n)
        ]
        sig_mode1_split_2b_inv = _hip_inverse(sig_mode1_split_2b)
        xt_m1_2b_x = [
            [
                sum(
                    Fraction(x23_rows[i][a]) * sig_mode1_split_2b_inv[i][j] * Fraction(x23_rows[j][b])
                    for i in range(n)
                    for j in range(n)
                )
                for b in range(23)
            ]
            for a in range(23)
        ]
        gls_mode1_split_2b_diag = [_hip_inverse(xt_m1_2b_x)[a][a] for a in range(23)]

        sig_mode1_split_1b = [
            [
                sig_mode1_split_2b[i][j] + (Fraction(49, 6) if (i == j and x23_rows[i][22] == 1) else Fraction(0, 1))
                for j in range(n)
            ]
            for i in range(n)
        ]
        sig_mode1_split_1b_inv = _hip_inverse(sig_mode1_split_1b)
        xt_m1_1b_x = [
            [
                sum(
                    Fraction(x23_rows[i][a]) * sig_mode1_split_1b_inv[i][j] * Fraction(x23_rows[j][b])
                    for i in range(n)
                    for j in range(n)
                )
                for b in range(23)
            ]
            for a in range(23)
        ]
        gls_mode1_split_1b_diag = [_hip_inverse(xt_m1_1b_x)[a][a] for a in range(23)]

        batch_sensitivity_sweep[b_sd_label] = {
            "assumed_common_differentiation_batch_sd_pp": float(math.sqrt(float(b_var_frac * 2))),
            "batches_per_stage": 2,
            "within_stage_batch_covariance_model": "COMPLETE_RUN_BY_BATCH_CROSSING_STAGE_COMMON_AVERAGE_TAU2_OVER_2",
            "main_effect_half_effect_beta_gls_se_range_pp": [
                round(math.sqrt(float(min(gls_b_diag[1:7]))), 4),
                round(math.sqrt(float(max(gls_b_diag[1:7]))), 4),
            ],
            "two_factor_interaction_half_effect_beta_gls_se_range_pp": [
                round(math.sqrt(float(min(gls_b_diag[7:22]))), 4),
                round(math.sqrt(float(max(gls_b_diag[7:22]))), 4),
            ],
            "stage_block_B_gls_se_pp": round(math.sqrt(float(gls_b_diag[22])), 4),
            "stage1_split_batch_s1a_vs_s1b_incidence_2batch_residual_gls_se_sol61_r3": {
                "covariance_mode": "MODE2_75_SHARED_GENOTYPE_CLONES_ACROSS_25_GENOTYPES",
                "main_effect_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_split_2b_diag[1:7]))), 4),
                    round(math.sqrt(float(max(gls_split_2b_diag[1:7]))), 4),
                ],
                "two_factor_interaction_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_split_2b_diag[7:22]))), 4),
                    round(math.sqrt(float(max(gls_split_2b_diag[7:22]))), 4),
                ],
                "stage_block_B_gls_se_pp": round(math.sqrt(float(gls_split_2b_diag[22])), 4),
            },
            "stage1_split_batch_s1a_vs_s1b_incidence_1batch_stage1_residual_gls_se_sol61_r3": {
                "covariance_mode": "MODE2_75_SHARED_GENOTYPE_CLONES_ACROSS_25_GENOTYPES",
                "main_effect_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_split_1b_diag[1:7]))), 4),
                    round(math.sqrt(float(max(gls_split_1b_diag[1:7]))), 4),
                ],
                "two_factor_interaction_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_split_1b_diag[7:22]))), 4),
                    round(math.sqrt(float(max(gls_split_1b_diag[7:22]))), 4),
                ],
                "stage_block_B_gls_se_pp": round(math.sqrt(float(gls_split_1b_diag[22])), 4),
            },
            "stage1_split_batch_s1a_vs_s1b_mode1_independent_run_nested_clones_2batch_residual_gls_se_kimi31_r3": {
                "covariance_mode": "MODE1_108_INDEPENDENT_RUN_NESTED_CLONE_DIFFERENTIATIONS",
                "main_effect_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_mode1_split_2b_diag[1:7]))), 4),
                    round(math.sqrt(float(max(gls_mode1_split_2b_diag[1:7]))), 4),
                ],
                "two_factor_interaction_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_mode1_split_2b_diag[7:22]))), 4),
                    round(math.sqrt(float(max(gls_mode1_split_2b_diag[7:22]))), 4),
                ],
                "stage_block_B_gls_se_pp": round(math.sqrt(float(gls_mode1_split_2b_diag[22])), 4),
            },
            "stage1_split_batch_s1a_vs_s1b_mode1_independent_run_nested_clones_1batch_stage1_residual_gls_se_kimi31_r3": {
                "covariance_mode": "MODE1_108_INDEPENDENT_RUN_NESTED_CLONE_DIFFERENTIATIONS",
                "main_effect_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_mode1_split_1b_diag[1:7]))), 4),
                    round(math.sqrt(float(max(gls_mode1_split_1b_diag[1:7]))), 4),
                ],
                "two_factor_interaction_half_effect_beta_gls_se_range_pp": [
                    round(math.sqrt(float(min(gls_mode1_split_1b_diag[7:22]))), 4),
                    round(math.sqrt(float(max(gls_mode1_split_1b_diag[7:22]))), 4),
                ],
                "stage_block_B_gls_se_pp": round(math.sqrt(float(gls_mode1_split_1b_diag[22])), 4),
            },
        }

    power_and_covariance_summary = {
        "stage_by_main_interaction_and_batch_covariance_qualification_ultracode_v80_12": {
            "stage_by_main_29col_expansion": {
                "n_rows": n,
                "n_columns_including_6_stage_by_main_interactions": 29,
                "exact_rational_rank": rank_29col_stage_by_main,
                "rank_deficiency": 29 - rank_29col_stage_by_main,
                "tier1_16run_rank": s["rank_tier1_22"],
                "tier1_plus_unreplicated_tier2_26run_rank_23col": s["rank23"],
                "active_tier1_main_effect_columns_note": (
                    "All 6 main effects (X1..X6) are active and clear of 2FIs in the 16-run Tier-1 screening block (rank 14/22); "
                    "adding the 8 Tier-2 de-aliasing genotypes + 2 bridge controls (26 runs) or the full replicated 36-run table "
                    "achieves full rank 22/22 (genetic) and 23/23 (with stage_block_B), whereas the 29-column expansion "
                    "[I, X1..X6, 15 2FIs, B, B*X1..B*X6] has rank 25/29 (deficiency 4) in the 36-run design, rank 26/29 (deficiency 3) "
                    "in the 38-run 4-bridge upgrade, and achieves full rank 29/29 with a minimum of 4 extra Stage-1 bridge runs (6 bridges total, 40 runs; e.g. Runs 01, 16, 02, 03, 04, 05) or in the 8-bridge 42-run upgrade (Findings 021, 075, 084)."
                ),
                "identifiability_note": (
                    "The two inter-stage bridge genotypes (Runs 01/25 and Runs 16/26) together with the Stage-1 de-aliasing and lead runs (Runs 17-24, 27-36) jointly identify the additive stage shift B "
                    f"(where genetic columns X1..X6 and 15 2FIs use {{-1, +1}} coding [full High-minus-Low contrast = 2*beta_j] while stage_block_B uses 0/1 indicator coding [0 = Stage 0 Runs 01-16, 1 = Stage 1 Runs 17-36, so beta_B is already the full Stage-1 minus Stage-0 shift], giving rank([X | B]) = 23/23 and projection residual SSR_B = ||(I - H_22) B||_2^2 = {_hip_frac_plain(s['ssr'])} = {float(s['ssr']):.4f} in 0/1 stage indicator coding [Var(beta_B) = {_hip_frac_plain(s['diag23'][22])} * sigma^2 = 0.688119 * sigma^2, SE(beta_B) = 0.8295 * sigma, MDE(beta_B) = 2.514 * sigma at df=13]; "
                    f"note that the two matched inter-stage bridge pairs alone contribute 1.0000 of the {float(s['ssr']):.4f} SSR_B in 0/1 indicator coding [68.81% of precision, with the remaining 63/139 = 0.4532 from the 18 non-bridge Stage-1 runs (Runs 17-24, 27-36); in {{-1, +1}} stage coding where B_pm = 2*B_01 - 1, SSR_B_pm = 4*SSR_B = 808/139 = 5.8129 with bridge sum 4.0000]), "
                    "whereas adding all 6 stage x main-effect terms B*X1..B*X6 yields 29 columns "
                    "with exact rank 25/29 (identifiability of the 23-column model assumes genetic main effects and 2FIs are invariant "
                    "across Stage 0 and Stage 1 blocks, or requires concurrent single-stage blocking / >=6 bridge genotypes in a >=40-run design; AUDIT-REF & AUDIT-REF)."
                ),
            },
            "differentiation_batch_random_intercept_sensitivity_sweep": batch_sensitivity_sweep,
        },
        "variance_assumptions": {
            "sigma_clone_pp": 8.0,
            "sigma_eps_donor_averaged_batch_cell_pp": 7.0,
            "n_clones_per_genotype_or_run": 3,
            "n_batches": 2,
            "n_donors": 3,
            "sigma_run_donor_averaged_batch_cell_pp": round(sigma_run_pp, 2),
            "sigma_run_if_7pp_is_per_donor_well_residual_divided_by_n_donors_pp": round(sigma_run_3donors_pp, 2),
            "clone_level_df_75_clones": 52,
            "clone_level_df_108_clones": 85,
            "n_distinct_design_rows_with_stage_indicator": n_distinct_with_stage,
            "conservative_run_table_df_23col": n - 23,
            "pure_error_df_within_stage_replicates_23col": n - n_distinct_with_stage,
            "lack_of_fit_df_higher_order_interactions_23col": n_distinct_with_stage - 23,
            "conservative_run_table_df_22col": n - 22,
        },
        "mode1_108_independent_run_nested_clones_exact_iid_run_covariance": {
            "allocation_description": "3 independently derived/edited hPSC clones per run across all 36 runs (108 independent founder clone units total), yielding strictly diagonal Cov(Y_bar) = sigma_run^2 * I_36.",
            "single_stage_22col_genetic_model": {
                "tier1_16run_rank": s["rank_tier1_22"],
                "combined_36run_rank_22col": s["rank22"],
                "main_effect_best_X1_to_X3": {**_se_mde_entry(min_main_22, pairs=["X5"]), "deprecated_key_alias_note": "Legacy key name retained for compatibility; actual min-variance main effect is X5 (see main_effect_best_X5)."},
                "main_effect_best_X5": _se_mde_entry(min_main_22, pairs=["X5"]),
                "main_effect_worst_X4_to_X6": {**_se_mde_entry(max_main_22, pairs=["X1", "X2", "X3"]), "deprecated_key_alias_note": "Legacy key name retained for compatibility; actual max-variance main effects are X1, X2, X3 (see main_effect_worst_X1_to_X3)."},
                "main_effect_worst_X1_to_X3": _se_mde_entry(max_main_22, pairs=["X1", "X2", "X3"]),
                "two_factor_interaction_best_6pairs": _se_mde_entry(min_2fi_22, pairs=best_2fi_pairs_22),
                "two_factor_interaction_worst_2pairs_X1X5_X4X6": {**_se_mde_entry(max_2fi_22, pairs=worst_2fi_pairs_22), "deprecated_key_alias_note": "Legacy key name retained for compatibility; actual 22-column worst-case 2FI is X2:X5 (1295/4448; see two_factor_interaction_worst_1pair_1295_over_4448_X2X5)."},
                "two_factor_interaction_worst_1pair_1295_over_4448_X2X5": _se_mde_entry(max_2fi_22, pairs=worst_2fi_pairs_22),
            },
            "two_stage_23col_model_with_additive_stage_block_B": {
                "combined_36run_rank_23col": s["rank23"],
                "rank_if_augmented_with_6_stage_by_main_29col": rank_29col_stage_by_main,
                "main_effect_best_X1_to_X3": {**_se_mde_entry(min_main_23, pairs=["X5"]), "deprecated_key_alias_note": "Legacy key name retained for compatibility; actual min-variance main effect is X5 (see main_effect_best_X5)."},
                "main_effect_best_X5": _se_mde_entry(min_main_23, pairs=["X5"]),
                "main_effect_worst_X4_to_X6": {**_se_mde_entry(max_main_23, pairs=["X1", "X2", "X3"]), "deprecated_key_alias_note": "Legacy key name retained for compatibility; actual max-variance main effects are X1, X2, X3 (see main_effect_worst_X1_to_X3)."},
                "main_effect_worst_X1_to_X3": _se_mde_entry(max_main_23, pairs=["X1", "X2", "X3"]),
                "two_factor_interaction_best_2pairs_115_over_1616_X4X5_X5X6": _se_mde_entry(min_2fi_23, pairs=best_2fi_pairs_23),
                "two_factor_interaction_worst_1pair_1191_over_1616_X1X3": _se_mde_entry(max_2fi_23, pairs=worst_2fi_pairs_23),
                "all_15_two_factor_interactions_by_contrast_name": {
                    HIP_TERM_NAMES[i]: _se_mde_entry(s["diag23"][i], pairs=[HIP_TERM_NAMES[i]])
                    for i in range(7, 22)
                },
                "stage_block_B": {
                    **_se_mde_entry(s["diag23"][22]),
                    "full_stage_difference_beta_B_se_pp": round(sigma_run_pp * math.sqrt(float(s["diag23"][22])), 4),
                    "full_stage_difference_beta_B_mde_80pct_alpha05_two_sided_df13_central_approx_pp": round(
                        (2.160369 + 0.841621) * sigma_run_pp * math.sqrt(float(s["diag23"][22])), 2
                    ),
                    "full_stage_difference_beta_B_mde_80pct_alpha05_two_sided_df13_pp": round(
                        0.5 * float(_se_mde_entry(s["diag23"][22])["mde_80pct_power_2beta_t13_alpha_0_05_pp"]), 2
                    ),
                    "coding_convention_note": (
                        "stage_block_B (column 23) is coded as a 0/1 indicator (0 = Stage 0 Runs 01-16, 1 = Stage 1 Runs 17-36), "
                        "so beta_B is already the full Stage-1 minus Stage-0 block difference (Var(beta_B) = 139/202 * sigma^2, "
                        "SE(beta_B) = 0.8295 * sigma = 4.5055 pp, MDE(beta_B) = 2.514 * sigma = 13.66 pp at df=13); the 2x multiplier fields in _se_mde_entry "
                        "apply only if stage_block_B is rescaled to {-1, +1} half-effect coding (AUDIT-REF & AUDIT-REF)."
                    ),
                },
            },
        },
        "mode2_75_shared_genotype_clones_across_25_distinct_genotypes_gls_and_ols_sandwich": {
            "allocation_description": (
                "75 independently edited hPSC clones (3 clones per each of the 25 distinct genotypes), where each run cell averages n_batches = 2 independent differentiation batches per run "
                "(giving v_eps_run = sigma_eps^2 / (n_clones * n_batches) = 49/6 pp^2 and diagonal run variance 29.50 pp^2 -> sigma_run = 5.43 pp; "
                "note that if each run were a single differentiation batch n_batches = 1 with bridged/replicated runs 25-26, 27-34, 36 serving as the second batch for those 11 genotypes, "
                "single-batch unbridged runs would have v_eps_run = 49/3 pp^2 -> sigma_run = sqrt(113/3) = 6.14 pp). "
                "Here Cov(Y_bar) = (sigma_clone^2 / n_clones) * Z Z^T + (sigma_eps^2 / (n_clones * n_batches)) * I_36: "
                f"shared clone random effects cancel inside the inter-stage bridge contrasts (cutting stage_block_B SE from "
                f"{round(sigma_run_pp * math.sqrt(float(s['diag23'][22])), 2):.2f} pp [Mode 1 i.i.d.] to "
                f"{round(math.sqrt(float(gls23_shared_var[22])), 2):.2f} pp [Mode 2 GLS] / "
                f"{round(math.sqrt(float(ols23_shared_var[22])), 2):.2f} pp [Mode 2 OLS sandwich]) "
                "while slightly increasing genetic SEs relative to 108 independent clones."
            ),
            "single_stage_22col_ols_sandwich_se_pp": {
                "main_effect_half_effect_beta_se_range_pp": [
                    round(math.sqrt(float(min(ols22_shared_var[1:7]))), 2),
                    round(math.sqrt(float(max(ols22_shared_var[1:7]))), 2),
                ],
                "main_effect_full_contrast_2beta_se_range_pp": [
                    round(2.0 * math.sqrt(float(min(ols22_shared_var[1:7]))), 2),
                    round(2.0 * math.sqrt(float(max(ols22_shared_var[1:7]))), 2),
                ],
                "two_factor_interaction_half_effect_beta_se_range_pp": [
                    round(math.sqrt(float(min(ols22_shared_var[7:22]))), 2),
                    round(math.sqrt(float(max(ols22_shared_var[7:22]))), 2),
                ],
                "two_factor_interaction_full_contrast_2beta_se_range_pp": [
                    round(2.0 * math.sqrt(float(min(ols22_shared_var[7:22]))), 2),
                    round(2.0 * math.sqrt(float(max(ols22_shared_var[7:22]))), 2),
                ],
            },
            "two_stage_23col_ols_sandwich_se_pp": {
                "main_effect_half_effect_beta_se_range_pp": [
                    round(math.sqrt(float(min(ols23_shared_var[1:7]))), 2),
                    round(math.sqrt(float(max(ols23_shared_var[1:7]))), 2),
                ],
                "main_effect_full_contrast_2beta_se_range_pp": [
                    round(2.0 * math.sqrt(float(min(ols23_shared_var[1:7]))), 2),
                    round(2.0 * math.sqrt(float(max(ols23_shared_var[1:7]))), 2),
                ],
                "two_factor_interaction_half_effect_beta_se_range_pp": [
                    round(math.sqrt(float(min(ols23_shared_var[7:22]))), 2),
                    round(math.sqrt(float(max(ols23_shared_var[7:22]))), 2),
                ],
                "two_factor_interaction_full_contrast_2beta_se_range_pp": [
                    round(2.0 * math.sqrt(float(min(ols23_shared_var[7:22]))), 2),
                    round(2.0 * math.sqrt(float(max(ols23_shared_var[7:22]))), 2),
                ],
                "stage_block_B_se_pp": round(math.sqrt(float(ols23_shared_var[22])), 2),
            },
            "two_stage_23col_gls_lmm_se_pp": {
                "main_effect_half_effect_beta_se_range_pp": [
                    round(math.sqrt(float(min(gls23_shared_var[1:7]))), 2),
                    round(math.sqrt(float(max(gls23_shared_var[1:7]))), 2),
                ],
                "main_effect_full_contrast_2beta_se_range_pp": [
                    round(2.0 * math.sqrt(float(min(gls23_shared_var[1:7]))), 2),
                    round(2.0 * math.sqrt(float(max(gls23_shared_var[1:7]))), 2),
                ],
                "two_factor_interaction_half_effect_beta_se_range_pp": [
                    round(math.sqrt(float(min(gls23_shared_var[7:22]))), 2),
                    round(math.sqrt(float(max(gls23_shared_var[7:22]))), 2),
                ],
                "two_factor_interaction_full_contrast_2beta_se_range_pp": [
                    round(2.0 * math.sqrt(float(min(gls23_shared_var[7:22]))), 2),
                    round(2.0 * math.sqrt(float(max(gls23_shared_var[7:22]))), 2),
                ],
                "stage_block_B_se_pp": round(math.sqrt(float(gls23_shared_var[22])), 2),
            },
        },
        # Top-level summary fields retained for direct inspection (explicitly labeling 22-col vs 23-col and best-case vs worst-case 2FI; Audit R03):
        "main_effect_22col_best_half_effect_beta_j_se_pp": round(sigma_run_pp * math.sqrt(float(min_main_22)), 2),
        "main_effect_22col_best_full_contrast_2beta_j_se_pp": round(2.0 * sigma_run_pp * math.sqrt(float(min_main_22)), 2),
        "main_effect_22col_best_full_contrast_mde_80pct_power_t52_alpha_0_05_pp": round(mult_t52_a05 * 2.0 * sigma_run_pp * math.sqrt(float(min_main_22)), 2),
        "main_effect_23col_best_half_effect_beta_j_se_pp": round(sigma_run_pp * math.sqrt(float(min_main_23)), 2),
        "main_effect_23col_best_full_contrast_2beta_j_se_pp": round(2.0 * sigma_run_pp * math.sqrt(float(min_main_23)), 2),
        "main_effect_23col_worst_half_effect_beta_j_se_pp": round(sigma_run_pp * math.sqrt(float(max_main_23)), 2),
        "main_effect_23col_worst_full_contrast_2beta_j_se_pp": round(2.0 * sigma_run_pp * math.sqrt(float(max_main_23)), 2),
        "best_case_2fi_23col_115_over_1616_half_effect_beta_jk_se_pp": round(sigma_run_pp * math.sqrt(float(min_2fi_23)), 2),
        "best_case_2fi_23col_115_over_1616_full_contrast_2beta_jk_se_pp": round(2.0 * sigma_run_pp * math.sqrt(float(min_2fi_23)), 2),
        "best_case_2fi_23col_115_over_1616_full_contrast_mde_80pct_power_t52_alpha_0_05_pp": round(mult_t52_a05 * 2.0 * sigma_run_pp * math.sqrt(float(min_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_half_effect_beta_jk_se_pp": round(sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_jk_se_pp": round(2.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_dod_contrast_4beta_jk_se_pp": round(4.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_normal_alpha_0_05_pp": round(mult_norm_a05 * 2.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_t52_alpha_0_05_pp": round(mult_t52_a05 * 2.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_full_contrast_2beta_mde_80pct_power_t13_alpha_0_05_pp": round(mult_t13_a05 * 2.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_dod_4beta_mde_80pct_power_normal_alpha_0_05_pp": round(mult_norm_a05 * 4.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_dod_4beta_mde_80pct_power_t52_alpha_0_05_pp": round(mult_t52_a05 * 4.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "worst_case_2fi_23col_1191_over_1616_dod_4beta_mde_80pct_power_t13_alpha_0_05_pp": round(mult_t13_a05 * 4.0 * sigma_run_pp * math.sqrt(float(max_2fi_23)), 2),
        "gsis_stimulation_index_22col_best_mde_80pct_power_at_sigma_run_0_45_SI_units": round(mult_t52_a05 * 2.0 * 0.45 * math.sqrt(float(min_main_22)), 2),
    }

    return {
        "protocol_fields": protocol_fields,
        "power_and_covariance_summary": power_and_covariance_summary,
        "n_distinct_genotypes": n_geno,
        "max_leverage_text": f"{_hip_frac_plain(max_h23)} = {float(max_h23):.6f}",
        "ssr_text": _hip_frac_plain(s["ssr"]),
        "all_2fi_estimable": s["rank22"] == 22,
        "full_rank_with_block": s["rank23"] == 23,
        "dealiasing_note": dealiasing_note,
    }


def _locus1_physical_cassette_id(x1: int, x2: int) -> str:
    """Explicit 4-cell X1 x X2 physical construct mapping (Audit R02): eliminates orphan 2A tags when X1 == -1."""
    if x1 == -1 and x2 == -1:
        return "LOC1_X1M_X2M_B2M_CD47WT_MONOCISTRONIC"
    if x1 == -1 and x2 == 1:
        return "LOC1_X1M_X2P_HEF1A_CD47GPI_MONOCISTRONIC"
    if x1 == 1 and x2 == -1:
        return "LOC1_X1P_X2M_B2M_HLAE_TPT2A_CD47WT_BICISTRONIC"
    return "LOC1_X1P_X2P_HEF1A_HLAE_TPT2A_CD47GPI_BICISTRONIC"


LOCUS1_X1_BY_X2_PHYSICAL_CONSTRUCT_MATRIX = {
    "LOC1_X1M_X2M_B2M_CD47WT_MONOCISTRONIC": {
        "coded_X1": -1,
        "coded_X2": -1,
        "cistron_valency": "Monocistronic (HLA-E absent; NO 2A peptide tag)",
        "physical_allele_architecture": "B2M exon 1 KO/KI :: endogenous resting B2M promoter -> Kozak-ATG-CD47-WT(full-length, OA3-323)-bGHpA",
        "purpose": "Tests endogenous resting B2M promoter driving full-length CD47-WT without HLA-E or any 2A peptide tag.",
    },
    "LOC1_X1M_X2P_HEF1A_CD47GPI_MONOCISTRONIC": {
        "coded_X1": -1,
        "coded_X2": 1,
        "cistron_valency": "Monocistronic (HLA-E absent; NO 2A peptide tag)",
        "physical_allele_architecture": "B2M exon 1 KO/KI :: B2M-ISRE-[native human EEF1A1-core (230 bp, 19 CpGs, 60.0% GC) + chimeric intron] -> Kozak-ATG-CD47-IgV-GPI(DAF)-WPREmut6-bGHpA",
        "purpose": "Tests high-floor constitutive/synthetic promoter driving signaling-incompetent GPI-anchored CD47-IgV-GPI without HLA-E or any 2A tag.",
    },
    "LOC1_X1P_X2M_B2M_HLAE_TPT2A_CD47WT_BICISTRONIC": {
        "coded_X1": 1,
        "coded_X2": -1,
        "cistron_valency": "Bicistronic tPT2A (HLA-E 5'-upstream, CD47-WT 3'-downstream)",
        "physical_allele_architecture": "B2M exon 1 KO/KI :: endogenous resting B2M promoter -> HLA-E*01:03(VMAPRTLVL)-GSG-tPT2A-CD47-WT(full-length, OA3-323)-bGHpA",
        "purpose": "Legacy lean chassis Locus 1 (Runs 17 & 27): tests endogenous B2M promoter driving bicistronic HLA-E + full-length CD47-WT.",
    },
    "LOC1_X1P_X2P_HEF1A_HLAE_TPT2A_CD47GPI_BICISTRONIC": {
        "coded_X1": 1,
        "coded_X2": 1,
        "cistron_valency": "Bicistronic tPT2A (HLA-E 5'-upstream, CD47-IgV-GPI 3'-downstream)",
        "physical_allele_architecture": "B2M exon 1 KO/KI :: B2M-ISRE-[native human EEF1A1-core (230 bp, 19 CpGs, 60.0% GC) + chimeric intron] -> HLA-E*01:03(VMAPRTLVL)-GSG-tPT2A-CD47-IgV-GPI(DAF)-WPREmut6-bGHpA",
        "purpose": "Current Scenario A Lead Chassis Locus 1 (Runs 35 & 36): tests high-floor EEF1A1-core+intron driving 5'-HLA-E*01:03 (retaining cytosolic -GSG-P2A tail) and 3'-CD47-IgV-GPI (with unscarred extreme C-terminal DAF GPI signal).",
        "embodiment_a_star_promoted_physical_cassette_id": "LOC1_X1P_X2P_ARM1B_ASTAR_4143BP_BIDIR_2AFREE",
        "embodiment_a_star_promoted_architecture_argon_r4": "B2M exon 1 KO/KI :: antisense native EEF1A1_core_230 (canonical pinned 4,143-bp Embodiment A* sha256 0995f8e35791d219... carries native CTTTTT 5'TOP; optional EEF1A1[5'TOP-CAGAGA] challenger sha256 28b6882c63d7e6fb...)-CD47-IgV-GPI + 154-bp SPA+MAZ4 pause shield / sense Syn_5xISRE_minP-HLA-E*01:03(VMAPRTLVL, 20-aa L1) (4,143 bp donor / 4,433 bp packaged in single AAV6; mandatory for unconfounded beta_{X1 x X4} Syn_5xISRE x TYK2 P1104A interaction estimation; AUDIT-REF & DEC-078)",
    },
    "design_identifiability_note": (
        "Because X1 = -1 cells are monocistronic (expressing CD47 alone at B2M exon 1) whereas X1 = +1 cells are bicistronic "
        "(co-expressing HLA-E*01:03 and CD47 via GSG-tPT2A), the coded X1 main effect and X1:X2 interaction jointly capture HLA-E surface "
        "presence plus any monocistronic-to-bicistronic translational stoichiometry shift on downstream CD47. Therefore, every factorial run "
        "records absolute QuantiBRITE surface molecules/cell for both HLA-E and CD47 as continuous covariates in a secondary mediation analysis, "
        "and Stage-0 Decision Experiment #1 compares the bicistronic tPT2A cassette against a 2A-free bidirectional/dual-promoter cassette (Audit R01/R02)."
    ),
}


def build_champion_challenger_portfolio(ranked_architectures: list[dict], all_combos: list[dict], wetlab_status: dict) -> dict:
    """Keep scored grid candidates separate from unmaterialized product hypotheses."""
    if not isinstance(ranked_architectures, list) or not isinstance(all_combos, list) or not all_combos:
        raise ValueError("Portfolio requires named architectures and a nonempty evaluated grid")
    arch_by_id = {a["architecture_id"]: a for a in ranked_architectures}
    if not {"ISLET1_CHASSIS_A_UNIVERSAL_LEAN_4LOCUS", "ISLET1_CHASSIS_B_HLA_C_RETAINED_6LOCUS_HYPOTHESIS"} <= arch_by_id.keys():
        raise ValueError("Portfolio requires both incumbent architecture identities")
    arch_a = arch_by_id["ISLET1_CHASSIS_A_UNIVERSAL_LEAN_4LOCUS"]
    arch_b = arch_by_id["ISLET1_CHASSIS_B_HLA_C_RETAINED_6LOCUS_HYPOTHESIS"]
    non_ab_feasible = [
        c for c in all_combos
        if c.get("is_feasible_candidate", True)
        and c["combination_id"] not in (arch_a["matched_combinatorial_search_id"], arch_b["matched_combinatorial_search_id"])
    ]
    if not non_ab_feasible:
        raise ValueError("Portfolio has no eligible non-incumbent grid challenger")
    non_ab_feasible.sort(key=lambda c: c["global_rank_among_all_combinations"])
    emergent_x = non_ab_feasible[0]
    zero_warn_x = next((c for c in non_ab_feasible if len(c.get("triggered_graveyard_rules", [])) == 0), None)
    er_stress_x = next((c for c in non_ab_feasible if "__ER_NONE__" not in c["combination_id"]), None)
    def watchlist_reference(candidate):
        if candidate is None:
            return {"status": "NO_QUALIFYING_FEASIBLE_STACK", "combination_id": None, "global_rank": None, "K": None, "B": None, "composite_series_system_utility": None}
        return {"status": "PRIOR_RANKED_WATCHLIST_ONLY", "combination_id": candidate["combination_id"],
                "global_rank": candidate["global_rank_among_all_combinations"], "K": candidate["total_genomic_cuts_K"],
                "B": candidate["total_base_edits_B"], "composite_series_system_utility": candidate["composite_series_system_utility"]}
    w_ped_pf, w_adult_pf = _composite_utility_weights()
    a_ped_u = arch_a["utility_pediatric_early_onset_t1d"]
    a_lada_u = arch_a["utility_adult_onset_lada"]
    a_g3_logit_unrounded = float(a_ped_u["gate_logits_unrounded"]["Gate3_Innate_NK_Macrophage_and_IBMIR_Shield"])
    a_g3_prob_unrounded = sigmoid_unclipped(a_g3_logit_unrounded)
    x_g3_prob_unrounded = sigmoid_unclipped(a_g3_logit_unrounded + 0.20)
    a_comp_unrounded = (
        w_ped_pf * float(a_ped_u["series_system_joint_survival_probability_unrounded"])
        + w_adult_pf * float(a_lada_u["series_system_joint_survival_probability_unrounded"])
    )
    x_comp_unrounded = a_comp_unrounded * (x_g3_prob_unrounded / a_g3_prob_unrounded)
    a_g3_disp = round(a_g3_prob_unrounded, 4)
    x_g3_disp = round(x_g3_prob_unrounded, 4)
    a_comp_disp = round(a_comp_unrounded, 4)
    x_comp_disp = round(a_comp_disp * x_g3_disp / a_g3_disp, 4)
    rel_gain_exact_pct = round(100.0 * (x_g3_prob_unrounded / a_g3_prob_unrounded - 1.0), 2)

    b_ped_u = arch_b["utility_pediatric_early_onset_t1d"]
    b_lada_u = arch_b["utility_adult_onset_lada"]
    b_ped_unrounded = float(b_ped_u["series_system_joint_survival_probability_unrounded"])
    b_lada_unrounded = float(b_lada_u["series_system_joint_survival_probability_unrounded"])
    b_comp_unrounded = w_ped_pf * b_ped_unrounded + w_adult_pf * b_lada_unrounded
    b_g2b_ped_logit = float(b_ped_u["gate_logits_unrounded"]["Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth"])
    b_g2b_lada_logit = float(b_lada_u["gate_logits_unrounded"]["Gate2B_Autoimmune_T1D_Memory_and_Indirect_APC_Stealth"])
    b_g3_logit = float(b_ped_u["gate_logits_unrounded"]["Gate3_Innate_NK_Macrophage_and_IBMIR_Shield"])
    b_ped_gy10 = b_ped_unrounded * (sigmoid_unclipped(b_g2b_ped_logit - 0.25) / sigmoid_unclipped(b_g2b_ped_logit))
    b_lada_gy10 = b_lada_unrounded * (sigmoid_unclipped(b_g2b_lada_logit - 0.25) / sigmoid_unclipped(b_g2b_lada_logit))
    b_comp_gy10 = w_ped_pf * b_ped_gy10 + w_adult_pf * b_lada_gy10
    b_g3_ratio_m015 = sigmoid_unclipped(b_g3_logit - 0.15) / sigmoid_unclipped(b_g3_logit)
    b_g3_ratio_m020 = sigmoid_unclipped(b_g3_logit - 0.20) / sigmoid_unclipped(b_g3_logit)
    b_ped_gy10_m015 = b_ped_gy10 * b_g3_ratio_m015
    b_comp_gy10_m015 = b_comp_gy10 * b_g3_ratio_m015
    b_ped_gy10_m020 = b_ped_gy10 * b_g3_ratio_m020
    b_comp_gy10_m020 = b_comp_gy10 * b_g3_ratio_m020

    champion_challenger_portfolio = {
        "decision_id": "DEC-007",
        "portfolio_model": "3-Slot Champion-Challenger Architecture (Blueprint A, Blueprint B, Blueprint X)",
        "outward_rankings_frozen": False,
        "slots": {
            "A": {
                "slot_id": "A",
                "architecture_id": arch_a["architecture_id"],
                "role": "Lean Universal Off-the-Shelf Champion (K=4, B=0)",
                "matched_combinatorial_search_id": arch_a["matched_combinatorial_search_id"],
                "global_rank": arch_a["global_rank_among_all_combinations"],
                "K": arch_a["genomic_locus_edit_count_K"],
                "B": arch_a["base_edit_count_B"],
                "composite_series_system_utility": arch_a["composite_series_system_utility"],
                "stage0_product_variant_ids": ["ISLET1_A4_ARM1B", "ISLET1_A_STAR_4LOCUS", "ISLET1_A4_ARM1A"],
            },
            "B": {
                "slot_id": "B",
                "architecture_id": arch_b["architecture_id"],
                "role": "HLA-C-Retained KIR-Matched Superdonor Comparator (Scored Grid #25: K=6, B=1 ISLET1_B1_6LOCUS_TYK2_KO_INS [0.0673] + K=5, B=1 ISLET1_B2_5LOCUS_TYK2_PE counterfactual [0.0782 complexity-discount sensitivity with gate logits held fixed; 0.0724 dropping +0.25 Gate-2B INS credit; 0.0569 full AUTO_NONE with TYK2-KO Gate-1 penalty or 0.0629 with TYK2 P1104A WT-yield]; Promoted Unscored Stage-0 Molecular Plan B'-SCT: K=4, B=0 ISLET1_B_SCT_C1_CIITA / ISLET1_B_SCT_C2_CIITA / ISLET1_B_SCT_C1C2_CIITA [counterfactual 0.1085])",
                "matched_combinatorial_search_id": arch_b["matched_combinatorial_search_id"],
                "global_rank": arch_b["global_rank_among_all_combinations"],
                "K": arch_b["genomic_locus_edit_count_K"],
                "B": arch_b["base_edit_count_B"],
                "composite_series_system_utility": arch_b["composite_series_system_utility"],
                "stage0_product_variant_ids": [
                    "ISLET1_B1_6LOCUS_TYK2_KO_INS",
                    "ISLET1_B2_5LOCUS_TYK2_PE",
                    "ISLET1_B_SCT_C1_CIITA",
                    "ISLET1_B_SCT_C2_CIITA",
                    "ISLET1_B_SCT_C1C2_CIITA",
                    "ISLET1_B_SCT_C1_CIITA_TYK2_PE",
                    "ISLET1_B_SCT_C2_CIITA_TYK2_PE",
                    "ISLET1_B_SCT_C1C2_CIITA_TYK2_PE",
                ],
            },
            "X": {
                "slot_id": "X",
                "architecture_id": "ISLET1_CHASSIS_X_EMERGENT_CHALLENGER",
                "grid_alias": "X-grid (ISLET1_CHASSIS_X_CONSERVATIVE_BASE_EDIT_CHALLENGER)",
                "role": "Highest-ranked feasible non-A/non-B grid challenger (K=4, B=1 conservative prime/base-edit challenger; exact utility tie at #2/#3 between #2 AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT + ER_NONE [K=4, B=1, discount = 0.9650] and #3 AUTO_NONE + ER_IAPP_PRAMLINTIDE_RNLS_WT [K=4, B=1, discount = 0.9650] at 0.10467791551016123 exact unrounded composite utility [0.10467792 8-dec; 0.1047 direct 4-dec; 0.1046 displayed via IEEE-754 binary64 cohort-rounded average round((0.1032 + 0.1061)/2, 4) = round(0.10464999999999999, 4) = 0.1046 (coinciding with exact-decimal half-even 0.10465 -> 0.1046, whereas literal binary64 round(0.10465, 4) = 0.1047)], resolved to #2 under reverse=True sort by the Gate-2B T1D autoantigen stealth priority tie-breaker for AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT while #3 ER_IAPP_PRAMLINTIDE_RNLS_WT serves as the Gate-4 ER/amyloid Pareto watchlist challenger, whereas K=5, B=0 AUTO_NONE + ER_RNLS_KO [discount = exp(-0.15) = 0.860708] ranks #9 at 0.08192907 [0.0819]); strictly distinct from off-grid Stage-0 molecular/biophysical Design X (K=4, B=0: ISLET1_A4_X_MOL_RNLS_WT, ISLET1_A4_DXS_ARM1C_L4CD47, ISLET1_A4_L3_BFB_BASALFLOOR, ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500; L1B-F1 & L1B-F2)",
                "evaluation_scope": "COMBINATORIAL_GRID_ONLY",
                "molecular_proposal_included_in_score": False,
                "matched_combinatorial_search_id": emergent_x["combination_id"],
                "global_rank": emergent_x["global_rank_among_all_combinations"],
                "K": emergent_x["total_genomic_cuts_K"],
                "B": emergent_x["total_base_edits_B"],
                "composite_series_system_utility": emergent_x["composite_series_system_utility"],
                "modules": emergent_x["modules"],
                "triggered_graveyard_rules": emergent_x.get("triggered_graveyard_rules", []),
                "stage0_molecular_product_variant_ids": [
                    "ISLET1_A4_X_MOL_RNLS_WT",
                    "ISLET1_A4_DXS_ARM1C_L4CD47",
                    "ISLET1_A4_L3_BFB_BASALFLOOR",
                    "ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500",
                ],
            },
        },
        "watchlist_qualification": "Prior-ranked candidates under declared filters, not experimentally established Pareto superiority",
        "pareto_challenger_watchlist": {
            "zero_graveyard_warn_challenger": watchlist_reference(zero_warn_x),
            "gate4_er_amyloid_resilience_challenger": watchlist_reference(er_stress_x),
        },
        "blueprint_x_stage0_molecular_challenger_spec": {
            "hypothesis_id": "ISLET1_X_STAGE0_MOLECULAR_HYPOTHESIS_K4_B0",
            "distinct_from_grid_challenger": True,
            "status": "UNMATERIALIZED_UNMEASURED_HYPOTHESIS",
            "execution_supported": False,
            "K": 4, "B": 0,
            "engineering_complexity_discount": engineering_complexity_discount(4, 0),
            "physical_manufacturing_yield": None,
            "exact_donor_length_bp": None,
            "packaged_length_bp": None,
            "candidate_arithmetic_donor_length_derivation_bp": 4056,
            "candidate_arithmetic_packaged_length_with_itrs_bp": 4346,
            "fully_loaded_arm1x_donor_length_with_spa_maz4_and_20aa_l1_bp": 4225,
            "fully_loaded_arm1x_packaged_length_with_itrs_bp": 4515,
            "fully_loaded_arm1x_aav6_headroom_bp": 185,
            "arm1x_part_sum_reconciliation_muse_r1": (
                "Reconciles the two Arm 1X length figures across the architecture: (a) 4,056 bp donor (4,346 bp with 2x145 bp ITRs) "
                "is the promoter-swap-only core (3,974 bp Arm 1B - 148 bp Syn_5xISRE_minP + 230 bp Syn_3xISRE_CpGminP_230bp_sense); "
                "(b) 4,225 bp donor (4,515 bp packaged, +185 bp AAV6 headroom below 4,700 bp) is the fully loaded Arm 1X assembly "
                "incorporating the 154-bp non-palindromic SPA+spacer+MAZ4 sense terminator (+154 bp), 20-aa (G4S)4 L1 linker (+15 bp), "
                "and 4-nt S3 PAM+seed shield (+0 bp). In the 8-run (24-clone) 2^3 factorial, SPA+MAZ4 (+154 bp), 20-aa L1 (+15 bp), and "
                "4-nt S3 shield are held fixed as common platform background across all 8 cells so W_P, W_D, and W_C remain orthogonal "
                "(AUDIT-REF & Q2-X1)."
            ),
            "external_de_novo_test_candidate_status": "UNMATERIALIZED_STAGE0_PROPOSAL_AWAITING_DE_NOVO_SYNTHESIS (arithmetic donor derivation: 3,974 bp Arm 1B - 148 bp Syn_5xISRE_minP + 230 bp sense Syn_3xISRE_CpGminP_230bp_sense [or 230 bp EEF1A1_core comparator] = 4,056 bp promoter-only core / 4,225 bp fully loaded with 154-bp SPA+MAZ4 + 15-bp 20-aa L1; exact_donor_length_bp remains None until nucleotide sequence is materialized on disk)",
            "resting_hlae_surface_density_molecules_per_cell": None,
            "resting_cd47_surface_density_molecules_per_cell": None,
            "composite_series_system_utility": None,
            "off_grid_molecular_design_x_locus1_and_lamina_counterfactual_projection_muse_r1": {
                "computation_mode": "LIVE_EVALUATOR_DERIVED_KIMI31_R3",
                "grid_champion_a_gate3_pediatric": a_g3_disp,
                "grid_champion_a_composite_utility": a_comp_disp,
                "off_grid_design_x_constitutive_hlae_gate3_pediatric_projection": x_g3_disp,
                "off_grid_design_x_constitutive_hlae_composite_utility_projection": x_comp_disp,
                "off_grid_design_x_constitutive_hlae_composite_utility_projection_exact_unrounded_inputs": round(x_comp_unrounded, 4),
                "unconditional_composite_utility_if_graveyard_009_retained": a_comp_disp,
                "off_grid_projection_conditionality_status_sonnet55_r2": "CONDITIONAL_ON_HOUR0_HLAE_15000_MOLECULES_PER_CELL_GATE_REMOVING_GRAVEYARD_009_MINUS_0P20_LOGIT",
                "relative_composite_gain_over_champion_a_pct": 7.06,
                "legacy_rounded_relative_composite_gain_pct_reference": 7.07,
                "relative_composite_gain_over_champion_a_exact_unrounded_inputs_pct": rel_gain_exact_pct,
                "relative_composite_gain_from_4dec_gate3_ratio_pct": round(100.0 * (x_g3_disp / a_g3_disp - 1.0), 2),
                "derivation_note": (
                    "While the 2,160-stack grid evaluates fixed locus indicators (where X-grid = ISLET1_CHASSIS_X_CONSERVATIVE_BASE_EDIT_CHALLENGER "
                    "has K=4, B=1, utility = 0.1046 [0.10467791551016123 exact unrounded / 0.10467792 8-dec / 0.1047 direct 4-dec; displayed as 0.1046 via IEEE-754 binary64 cohort-rounded average round((0.1032 + 0.1061)/2, 4) = round(0.10464999999999999, 4) = 0.1046 (coinciding with exact-decimal half-even 0.10465 -> 0.1046, whereas literal binary64 round(0.10465, 4) = 0.1047)]), off-grid molecular Design X (K=4, B=0, discount = 1.0000) replaces Arm 1B's resting-attenuated "
                    "Syn_5xISRE_minP promoter (10,890 molecules/cell at Hour 0 [< 15,000 NKG2A threshold; 0.6050x WT 18,000], incurring the -0.20 resting NK vulnerability penalty on Gate 3) "
                    "with constitutive HLA-E expression (Arm 1C promoter-trapped B2M::HLA-E SCT in DX-S or Arm 1A-lite EEF1A1_core) plus "
                    "AAVS1::iCasp9-GSG-P2A-CD47-IgV-GPI and the 100 um PDMS-CaO2 lamina (eps = 0.42, <=100 um clusters). Removing the "
                    "-0.20 Hour-0 resting HLA-E discount (GRAVEYARD_009; sigma(0.5598637 + 0.20) = 0.681324 vs 0.636421 while GRAVEYARD_010's -0.25 stays) increases Gate 3 pediatric survival from 0.6364 to 0.6813 (+4.49 pp) and lifts "
                    f"projected off-grid K=4, B=0 composite utility from 0.1085 ({a_comp_unrounded:.8f} unrounded) to 0.1162 (0.1085 * 0.6813 / 0.6364 = 0.116155 from displayed 4-decimal inputs; 0.1161 [{x_comp_unrounded:.8f}] from exact unrounded evaluator inputs; exact relative gain is +7.06% [7.0556% unrounded / 7.0553% from 0.6813/0.6364; legacy +7.07% retained in legacy_rounded_relative_composite_gain_pct_reference]; strictly conditional on passing Decision Experiment #1 Hour-0 HLA-E >= 15,000 molecules/cell; unconditional baseline equals 0.1085; AUDIT-REF, D1-02, Q2-X4, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF & L1B-F2 [DEC-077])."
                ),
            },
            "genomic_budget": "Proposed K=4, B=0; utility discount 1 is not measured manufacturing yield",
            "locus1_arm_1x_challenger": "Proposed bidirectional synthetic-promoter HLA-E/CD47 cassette with WT alpha3 Asp227/Thr228 conditional default (W_D = -1 in the 8-run 2^3 Stage-0 factorial) vs optional HLA-E mature-heavy-chain D227K/T228A CD8-null satellite coupon (W_D = +1; activated only if Gate 0A-NK SPR CD8 K_D > 1,000 uM, CD94/NKG2A K_D <= 2.5 uM, 0/6 CMV+ UL40 CTL donors > 5.0% lysis, and >= 80% CD8+ Treg suppression retention all pass; R9-S1-L3A-01 [DEC-077]) and CD47 mature-IgV C15G substitution. No exact Arm 1X donor ledger, structural qualification or expression measurement is available. CD8-null and improved surface-density function remain hypotheses.",
            "seed_line_and_perigraft_shield": "Proposed genotype-selected seed line and spatial SA-THBD/SA-FHL-1 plus conditional SA-FasL/CXCL12 microgels. INS rs3842752 identity, allele direction and functional transfer require independent qualification; seed selection is not an editing event, but does not establish ER-stress or DRiP reduction.",
            "gate5_annular_lamina_topology": "Layout B is a conditional Layer 4 geometry scenario. Its simulated oxygen values do not validate this molecular product or change the grid score; use the regenerated Layer 4 artifact for scenario parameters and numerical outputs.",
            "pending_assay_gates": ["Exact nucleotide and edited-allele ledger", "Lineage-specific expression and HLA-E receptor function", "INS allele/function qualification", "Spatial immune and complement assays", "Matched physical oxygen and viability qualification"],
            "stage0_arm1x_2cubed_deconfounding_subdesign_astra_r1": {
                "design_id": "STAGE0_ARM1X_2CUBED_COMPONENT_DECONFOUNDING_FACTORIAL",
                "runs": 8,
                "clones_per_run": 3,
                "total_clonal_evaluations": 24,
                "external_unmodified_3974bp_arm1b_reference_control_clones": 3,
                "fixed_platform_background_across_all_8_cells": (
                    "154-bp SPA+spacer+MAZ4 sense pause/terminator + 20-aa (G4S)4 L1 linker + 4-nt S3 PAM+seed silent shield "
                    "(anchored against 3 parallel clones of unmodified 3,974 bp canonical Arm 1B so the +169 bp platform shift is directly measured; AUDIT-REF)"
                ),
                "factors_coded_minus1_plus1": {
                    "W_P": "-1 = Syn_5xISRE_minP_148bp (Arm 1B / Embodiment A* baseline, 4,143 bp donor background; represented in the product registry by materialized ISLET1_A4_X_MOL_RNLS_WT) vs +1 = Syn_3xISRE_CpGminP_230bp_sense (unmaterialized Arm 1X hybrid constitutive-floor sense promoter challenger, 4,056 bp promoter-only core / 4,225 bp fully loaded donor; optional constitutive EEF1A1_core_230 sense comparator; DEC-078)",
                    "W_D": "-1 = HLA-E WT alpha3 D227/T228 vs +1 = HLA-E alpha3 D227K/T228A CD8-null substitution",
                    "W_C": "-1 = CD47 WT Cys15 vs +1 = CD47 C15G (C33G precursor) free-thiol substitution",
                },
                "orthogonal_variance_multiplier_in_plus_minus_1_coding": 0.125,
                "full_contrast_2beta_variance_multiplier": 0.500,
                "main_effect_and_interaction_se_pp_at_sigma_run_4p6188": 1.6330,
                "main_effect_and_interaction_mde_80pct_power_two_sided_alpha05_df13_pp": 4.95,
                "full_biological_contrast_2beta_se_pp_at_sigma_run_4p6188": 3.2660,
                "full_biological_contrast_2beta_mde_80pct_power_two_sided_alpha05_df13_pp": 9.90,
                "two_batch_hip_core_matched_sigma_run_pp": 5.4314,
                "main_effect_and_interaction_se_pp_at_two_batch_sigma_run_5p4314": 1.9203,
                "main_effect_and_interaction_mde_80pct_power_two_sided_alpha05_df16_at_5p4314_pp": 5.82,
                "main_effect_mde_exact_df16_noncentral_t_at_5p4314_pp": 5.7312,
                "full_biological_contrast_2beta_se_pp_at_two_batch_sigma_run_5p4314": 3.8406,
                "full_biological_contrast_2beta_mde_80pct_power_two_sided_alpha05_df16_at_5p4314_pp": 11.64,
                "full_biological_contrast_2beta_mde_exact_df16_noncentral_t_at_5p4314_pp": 11.4624,
                "standalone_two_sided_power_for_5pp_full_contrast_at_n8_two_batch": 0.2558,
                "standalone_two_sided_power_exact_df16_noncentral_t_for_5pp_full_contrast": 0.2318,
                "standalone_two_sided_power_approximation_note_muse_r2": {
                    "exact_df16_noncentral_t_power": 0.2318,
                    "normal_approx_power_unrounded_sigma_run_sqrt_177_over_6_5p431390": 0.2558,
                    "normal_approx_power_3dec_z_rounded_minus0p658": 0.2558,
                    "normal_approx_power_rounded_sigma_run_5p43_exact": 0.2559,
                    "normal_approx_power_rounded_sigma_run_5p43": 0.2559,
                    "noncentrality_parameter_delta_5pp_over_3p840573": 1.301889,
                    "clone_only_sigma_eps_zero_sigma_run_sqrt_64_over_3_pp": 4.6188,
                    "single_batch_with_assay_residual_sigma_run_sqrt_113_over_3_pp": 6.1373,
                    "note": (
                        "Reconciles the power and variance-component figures: 0.2318 is the exact df = 16 two-sided noncentral-t power at alpha = 0.05 "
                        "(t_crit = 2.119905, delta = 5.0 / sqrt(177/12) = 1.301889), 0.2558 (0.255779 unrounded; 0.255815 at 3-decimal z-rounding 1.302 - 1.960 = -0.658) is the standard normal approximation at unrounded "
                        "sigma_run = sqrt(177/6) = 5.431390 pp, 0.2559 (0.255905) is the exact normal approximation at 2-decimal rounded sigma_run = 5.43 pp, and sigma_run = 4.6188 pp (sqrt(64/3)) is the "
                        "clone-only (sigma_eps = 0) lower-bound reference whereas a single batch with sigma_eps = 7.0 pp has sigma_run = sqrt(113/3) = 6.1373 pp "
                        "(half-effect SE = 2.1699 pp, full-contrast SE = 4.3397 pp; AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                    ),
                },
                "clones_per_cell_required_for_standalone_80pct_power_on_5pp_full_contrast": 15,
                "clones_per_cell_if_batch_and_assay_residual_omitted_legacy_r9": 9,
                "clones_per_cell_exact_power_at_r9_under_two_batch_variance_model": 0.6030,
                "clones_per_cell_required_for_exact_80pct_power_under_two_batch_variance_model_r15": 15,
                "clones_per_cell_exact_power_at_r15_df112_under_two_batch_variance_model": 0.8228,
                "arm1x_tier_endpoint_role_distinction_sonnet55_r2": {
                    "tier1_r3_molecular_screen_24_clones": (
                        "r = 3 clones/cell (24 factorial clones + 3 unmodified 3,974-bp Arm 1B controls + 3 mandatory 15-aa L1 controls = 30 clones) "
                        "serves as the Stage-0 molecular-endpoint screen for high-SNR biochemical readouts (W_P: QuantiBRITE Hour-0 HLA-E >= 15,000 molecules/cell; "
                        "W_D: recombinant CD8alpha/beta SPR K_D > 1,000 uM; W_C: DTNB/maleimide-PEG2-biotin free-thiol elimination & SIRPa avidity)."
                    ),
                    "tier2_r15_functional_lysis_confirmation": (
                        "r = 15 clones/cell (df = 112, 82.28% exact noncentral-t power at two-batch sigma_run = 5.4314 pp) is required only for "
                        "standalone 5.0-pp functional immune-evasion / lysis contrast confirmation on surviving parsimonious sub-combinations (AUDIT-REF)."
                    ),
                    "external_control_pairwise_two_sample_se_note_kimi31_r3": (
                        "Note on contrast topology (AUDIT-REF & AUDIT-REF): because sigma_run(r)^2 = (64 + 49/2)/r = 88.5 / r already divides by the within-cell clone count r, "
                        "SE(2*beta_j) = sqrt(88.5 / (2*r)) = 3.8406 pp at r = 3 (df = 16, 23.18% power) or 1.7176 pp at r = 15 (df = 112, 82.28% power) "
                        "applies to each orthogonal 4-cell vs 4-cell main/interaction contrast within the 2^3 factorial, whereas pairwise two-sample comparisons between a single "
                        "factorial cell (r_cell clones) and an external reference arm (r_ctrl clones of unmodified 3,974-bp Arm 1B or canonical 15-aa L1) "
                        "have SE_pairwise = sqrt(sigma_run(r_cell)^2 + sigma_run(r_ctrl)^2) = sqrt(88.5 * (1/r_cell + 1/r_ctrl)): at r_cell = r_ctrl = 3, SE_pairwise = sqrt(88.5 * 2/3) = 7.6811 pp "
                        "(two-sample t df = 4, power = 8.04% for 5.0 pp); at r_cell = 15 with r_ctrl = 3, SE_pairwise = sqrt(88.5 * (1/15 + 1/3)) = 5.9498 pp (df = 16, power = 12.41%); "
                        "and if the external reference arm is co-expanded to r_ctrl = 15 alongside r_cell = 15, SE_pairwise = sqrt(88.5 * 2/15) = 3.4351 pp (df = 28, exact noncentral-t power = 29.01% [0.290119; 29.07% is an earlier approximation]; AUDIT-REF)."
                    ),
                },
                "matched_l1_linker_15aa_vs_20aa_satellite_contrast_astra_r2": (
                    "MANDATORY Stage-0 satellite arm (3 clones of 4,128-bp Arm 1B + SPA+MAZ4 + S3 with canonical 15-aa L1=(G4S)3 vs 3 clones of 4,143-bp L1-Run 1 with 20-aa L1=(G4S)4, "
                    "plus 3 clones on the promoted W_P/W_D/W_C background at 4,210 bp vs 4,225 bp) to deconfound the +15-bp 20-aa L1=(G4S)4 linker extension from the +154-bp SPA+MAZ4 "
                    "terminator shield and adopt L1=(G4S)4 as the default platform linker across A*, B-prime, and X unless 15-aa achieves >= 20-aa Hour-0 HLA-E density and NKG2A+ protection "
                    "(AUDIT-REF & AUDIT-REF)."
                ),
                "mandatory_satellite_status": "MANDATORY_STAGE0_LINKER_DECONFOUNDING_CONTROL",
                "estimand_scale_reconciliation_note_sol61_r2": (
                    "In {-1, +1} coding, beta_j = (mu_{+1} - mu_{-1}) / 2 is the regression half-effect (Var(beta_j) = sigma_run^2 / 8 = 0.125 * sigma_run^2, "
                    "SE = 1.6330 pp single-batch / 1.9203 pp two-batch, MDE_80%(beta_j) = 4.95 pp [df=13] / 5.7312 pp exact df=16 noncentral-t [5.82 pp conservative df=13 multiplier]), "
                    "whereas the full biological level difference Delta_j = mu_{+1} - mu_{-1} = 2 * beta_j has Var(2*beta_j) = sigma_run^2 / 2 = 0.500 * sigma_run^2, "
                    "SE(2*beta_j) = 3.2660 pp single-batch / 3.8406 pp two-batch, and MDE_80%(2*beta_j) = 9.90 pp (df=13) / 11.4624 pp exact df=16 noncentral-t (11.64 pp conservative df=13 multiplier). "
                    "A +5.0 pp full contrast corresponds to beta_j = +2.50 pp (23.18% exact df=16 noncentral-t / 25.58% normal-approx standalone two-sided power at r=3 clones/cell [N=8 runs, delta = 5.0 / 3.840573 = 1.301889]; "
                    "under sigma_run^2(r) = 64/r + 49/(2r) with df = 8(r-1), r = 9 clones/cell yields SE(2*beta_j) = 2.2174 pp and 60.30% power, while r = 15 clones/cell "
                    "[SE(2*beta_j) = 1.7176 pp, df = 112] achieves 82.28% >= 80% standalone power; AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
                ),
                "hierarchical_parsimony_promotion_rule": (
                    "Promote the full 3-change Arm 1X (+1, +1, +1) only if each positive factor improves its target "
                    "readout (W_P: resting Hour-0 HLA-E >= 15,000 molecules/cell; W_D: CD8alpha/beta SPR KD > 1,000 uM without ER export loss; "
                    "W_C: free-thiol elimination and SIRPa avidity) without degrading any co-primary gate; otherwise promote the highest-scoring "
                    "parsimonious sub-combination (closing AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
                ),
                "locus3_w_bf_basal_floor_2arm_satellite_muse_r3": {
                    "satellite_id": "STAGE0_LOCUS3_W_BF_BASAL_FLOOR_2ARM_SATELLITE",
                    "n_arms": 2,
                    "clones_per_arm": 3,
                    "total_clones": 6,
                    "arms": {
                        "W_BF_minus1_pure_inducible_dxl3_3799bp": "ISLET1_A4_L3_SCARFREE_3799 (3,799-bp scar-free dual-promoter CXCL10::SERPINB9 + Syn_5xISRE_minP->PD-L1)",
                        "W_BF_plus1_basal_floor_l3_bfb_4078bp": "ISLET1_A4_L3_BFB_BASALFLOOR (4,078-bp donor / 4,368-bp packaged: 49-bp SPA + 230-bp basal-floor promoter -> scar-free SERPINB9 + 82-bp pause + Syn_5xISRE_minP -> PD-L1)",
                    },
                    "decision_experiment_3_split_muse_r3": {
                        "DECISION_EXP_3A_2A_STUB_PROTEIN_QUALITY": "Compares canonical 4,090-bp SERPINB9-GSG-P2A-PD-L1 against scar-free 3,799-bp DX-L3 for SERPINB9 half-life (t1/2 >= 6 h) and Granzyme B k_ass (>= 5.0e5 M^-1 s^-1).",
                        "DECISION_EXP_3B_HOUR0_KINETIC_RACE_AND_BASAL_FLOOR": "Compares pure-inducible DX-L3 (3,799 bp) against basal-floor L3-BF-B (4,078 bp) in a 4-h unstimulated Hour-0 perforin/Granzyme B race (requiring [SERPINB9]_0 >= 0.30 uM, 4-h Annexin-V+ <= 15%, and 24-h IFN-g [SERPINB9]_24 >= 0.80 uM; Gate 0C-kinetic & Rule S1).",
                    },
                },
            },
        },
        "two_step_promotion_gate": {
            "step1_in_silico_qualification": "Evaluate the exact identified challenger with the same declared utility and edit accounting as the incumbent; no BLOCK, no higher-K exemption from the complexity discount, and preregistered endpoint trade-offs. A prior score does not establish physical qualification.",
            "step2_empirical_promotion": "Cannot displace Blueprint A or B without independently reviewed exact-product empirical evidence, all co-primary Stage-0/Stage-1 gates and an explicit promotion decision. Running the ingestor or applying any posterior shift alone is insufficient.",
            "promotion_authorized": False,
            "exact_product_empirical_gate_results": None,
            "current_promotion_status": (
                "CALIBRATION_ACTIVE_PROMOTION_NOT_AUTHORIZED"
                if wetlab_status.get("applied_to_production_priors")
                else "PRIOR_ANCHORED_CHALLENGER_AWAITING_WETLAB_READOUTS"
            ),
        },
        "abxyz_champion_contender_wildcard_tournament_registry_dec020": {
            "decision_ids": ["DEC-020", "DEC-021", "DEC-022", "DEC-023", "DEC-024", "DEC-025", "DEC-026", "DEC-027", "DEC-028", "DEC-029", "DEC-030", "DEC-036"],
            "review_focus_allocation_protocol": "Strictly locked to Stage-0 Execution Triad (A / B / X) and supporting Layers 1-5 (DEC-026..DEC-036 scope lock superseding historical DEC-020 80%/20% allocation; existing Y/Z wildcard specs preserved without expansion)",
            "tier1_stage0_execution_triad_80pct_focus": {
                "A_Champion_Islet1": {
                    "slot": "A",
                    "role": "Incumbent Universal Off-the-Shelf Champion on the throne (K=4, B=0)",
                    "architecture_summary": "Locus 1 B2M::Arm 1B (3,974 bp) vs Arm 1A (4,880 bp) + Locus 2 CIITA Exon 5 KO + Locus 3 CXCL10::SERPINB9-GSG-P2A-PD-L1 (4,090 bp) + Locus 4 AAVS1::cHS4-A2UCOE-EEF1A1-iCasp9 (3,457 bp) + 60:30:10 <=120 um clusters in 350 um PDMS-CaO2/rhCAT pouch",
                    "b2m_null_orphan_hc_erad_status_muse_r3": "B2M_NULL_ORPHAN_HC_ERAD_UNQUALIFIED",
                    "embodiment_a_star_platform_upgrade_and_conditionality_sonnet55_r2": {
                        "canonical_a_rank_conditionality_status": "CONDITIONAL_ON_LAYER3_15AA_L1_GEODESIC_AND_LAYER4_119UM_LAMINA_SURFACE_GSIS_QUALIFICATION",
                        "defended_champion_status": "EMBODIMENT_A_STAR_PROMOTED_AS_DEFENDED_CHAMPION_CONFIGURATION",
                        "embodiment_a_star_locus1_donor_length_bp": 4143,
                        "embodiment_a_star_locus1_packaged_aav6_length_bp": 4433,
                        "embodiment_a_star_definition": (
                            "Embodiment A* is the defended Champion configuration and pairs the canonical Blueprint A molecular stack (K=4, B=0) with the zero-cost platform upgrades shared with X: "
                            "(1) Locus 1 at 4,143 bp donor / 4,433 bp packaged (3,974 bp Arm 1B + 154 bp SPA+spacer+MAZ4 sense read-through terminator shield + 15 bp 20-aa L1=(G4S)4 linker, "
                            "matching L1-Run 1 of the Arm 1X 2^3 sub-design and passing the Layer-3 d < 0.65*Lc surface-route screen at 13.84-22.96 pN; without the +154-bp SPA+MAZ4 shield the linker-only intermediate is 3,989 bp donor / 4,279 bp packaged), "
                            "(2) 100 um thinned central PDMS-CaO2 lamina (a = 100/350, 125 um outer bands, <=100 um micro-islet cluster cap, passing Layer-4 coupled_complete_product_pass with 18.44 mmHg gel-only cluster-surface pO2 >= 15.0 mmHg), and "
                            "(3) unconditional 3-locus recut-immunity PAM-shield bundle (B2M 4-nt S3 shield, CXCL10 chr4:76023451C>G, AAVS1 chr19:55115754C>G), isolating A* vs X strictly to the molecular hypotheses (W_P, W_D, W_C, DX-S, DX-L3; AUDIT-REF, D4-01, D5-01 & AUDIT-REF)."
                        ),
                    },
                },
                "B_Stratified_KIR_SubBank_PlanB_and_Bprime": {
                    "slot": "B",
                    "role": "HLA-C-Retained KIR-Matched Stratified Sub-Bank Comparator (Plan B / Design B-prime; Q2-P1..P3)",
                    "architecture_summary": (
                        "Retains 1 KIR-matched HLA-C allele with functional SP-1C VMAPRTLIL leader (C1: HLA-C*01:02 or C*03:04; "
                        "C2: HLA-C*04:01 or C*06:02; avoiding non-functional SP-3C VMAPRALLL [P6=Ala, T6A with P5=Arg] in HLA-C*07:02 per Lin et al. 2023 PMID:37264229, PMC10690437, Table 1) "
                        "+ TYK2 P1104A (rs34536443) hypomorph + Design B-prime AAVS1::HLA-E*01:03(SCT) + CD47-IgV-GPI + iCasp9 cloak (Q2-P1..P3; "
                        "disclosing both 91.37% conservative unmaterialized-cloak-omitted and 92.84% complexity-discounted survival; AUDIT-REF & AUDIT-REF)"
                    ),
                    "b_prime_4digit_hlac_rearming_and_b_prime_sct_promotion_muse_r3": {
                        "primary_manufacturable_b_prime_status": "B_PRIME_SCT_PROMOTED_AS_PRIMARY_MANUFACTURABLE_SLOT_B_CANDIDATE",
                        "endogenous_hlac_retained_status": "EXPLORATORY_REQUIRES_4DIGIT_HLA_C_MATCH_OR_7_TO_12_LINE_SUPERDONOR_BANK",
                        "primary_product_variant_ids": ["ISLET1_B_SCT_C1_CIITA", "ISLET1_B_SCT_C2_CIITA", "ISLET1_B_SCT_C1C2_CIITA"],
                        "primary_citations": ["PMID:30853558", "PMID:21743387", "PMID:10850706", "PMID:11323700", "PMID:10429675", "PMID:17475856", "PMID:37207206"],
                        "tethered_human_self_peptides_and_f_pocket_clearance_argon_r4": {
                            "C1_HLA_C0304_self_peptide": "GAVDPLLAL (human KPNA1/KPNA2 importin-alpha self-peptide crystallized with HLA-C*03:04:KIR2DL2 in PDB:1EFX, monomeric KIR2DL2 K_D = 7.2 uM; replaces viral Influenza A M1 58-66 GILGFVFTL)",
                            "C2_HLA_C0401_self_peptide": "QYDDAVYKL (crystallographic consensus nonamer designed from HLA-Cw4 pool sequencing and crystallized with HLA-C*04:01:KIR2DL1 in PDB:1IM9, monomeric K_D = 7.2 uM for KIR2DL1 and 30.0 uM for activating KIR2DS1; Fan & Wiley 1999 PMID:10429675, Fan et al. 2001 PMID:11323700, Stewart et al. 2005 PMID:16141329; not an established endogenous TUBB self-peptide, paired with verified endogenous HLA-C*04:01 self-peptide challenger in Gate 0B-SCT; AUDIT-REF, AUDIT-REF & AUDIT-REF)",
                            "l1_and_l2_linkers": "20-aa (G4S)4 L1 and 20-aa (G4S)4 L2 with HLA-C Tyr84Ala lateral F-pocket exit channel evaluated side-by-side against three mechanistically distinct comparators: (1) WT-Y84 with a 25-aa (G4S)5 L1 linker (0 engineered cysteines), (2) DT-1 peptide-linker disulfide trap HLA-C(Y84C, A139-WT) + L1=GCGGS(G4S)3 (2 engineered cysteines; Truscott et al. 2007 PMID:17475856), and (3) DT-2 alpha1-alpha2 F-pocket staple HLA-C(Y84C/A139C) with non-cysteine (G4S)4 L1 (2 engineered cysteines; Finton et al. 2023 PMID:37207206), because native Tyr84:OH/CE1 directly contact KIR2DL Asp183:OD2 at 3.24-3.56 A in PDB:1EFX/1IM9 (Moradi et al. 2021 PMID:33846289; AUDIT-REF & AUDIT-REF)",
                            "ciita_exon5_constitutive_insulated_cassette": "EEF1A1_core_230(5'TOP-CAGAGA) + Chimeric_Intron(133) + B2M_SP-SelfPeptide-(G4S)4-B2M-(G4S)4-HLA-C(Y84A) + bGHpA(225) + SPA+MAZ4(154, insulating adjacent DEXI on chr16p13.13) (3,465 bp donor / 3,755 bp packaged)",
                        },
                        "sol61_r4_heavy_chain_dsa_kir2ds1_and_y84a_contact_qualification": {
                            "evidence_class": "STRUCTURAL_AND_ALLOREACTIVITY_GATE_QUALIFICATION",
                            "primary_citations": ["PMID:24804568", "PMID:16034111", "PMID:11907065", "PMID:10429675", "PMID:16141329", "PMID:23554313", "PMID:33846289", "PMID:17475856", "PMID:37207206"],
                            "polymorphic_heavy_chain_dsa_and_single_pmhc_tcr_recognition": (
                                "Covalent single-peptide tethering eliminates endogenous beta-cell autoantigen presentation (PPI, IGRP, GAD65, ZNT8) and reduces polyclonal peptide-dependent alloreactivity, "
                                "but does not erase solvent-exposed polymorphic alpha1/alpha2 heavy-chain eplets on HLA-C*03:04 or HLA-C*04:01 recognized by pre-formed/de novo anti-HLA-C DSAs (Aubert et al. 2014 PMID:24804568) "
                                "or high-density single-pMHC alloreactive CD8+ TCRs (Whitelegg et al. 2005 PMID:16034111; Yu et al. 2002 PMID:11907065); furthermore, QYDDAVYKL is a synthetic consensus nonamer capable of eliciting peptide-focused CD8+ TCR responses even in HLA-C*04:01-matched hosts. Therefore, B'-SCT requires both (a) a negative recipient-serum anti-HLA-C CDC/ADCC flow crossmatch against the exact B'-SCT product and "
                                "(b) a negative recipient CD8+ MLC (<=10% 4-h lysis, <=5% 72-h CFSE-low) across ALL recipients (including 4-digit matched hosts on QYDDAVYKL arms)."
                            ),
                            "c2_qyddavykl_kir2ds1_activating_receptor_gate": (
                                "HLA-C*04:01:QYDDAVYKL binds activating KIR2DS1 with K_D = 30 uM vs inhibitory KIR2DL1 K_D = 7.2 uM (Stewart et al. 2005 PMID:16141329) and triggers cytolysis by KIR2DS1+/KIR2DL1- NK subsets (~33% of KIR B-haplotype donors; Pittari et al. 2013 PMID:23554313); "
                                "all C2-containing B'-SCT products (ISLET1_B_SCT_C2_CIITA, ISLET1_B_SCT_C2_CIITA_TYK2_PE, ISLET1_B_SCT_C1C2_CIITA, ISLET1_B_SCT_C1C2_CIITA_TYK2_PE) require primary sorted KIR2DS1+/KIR2DL1- NK non-activation (CD107a <=10%, 4-h lysis <=10%) alongside an endogenous HLA-C*04:01 self-peptide challenger."
                            ),
                            "tyr84_kir_asp183_direct_contact_and_sct_comparators": (
                                "In PDB:1EFX (HLA-C*03:04:GAVDPLLAL:KIR2DL2), Tyr84:OH lies 3.558 A from KIR2DL2 Asp183:OD2 (4 heavy-atom pairs < 4.5 A); in PDB:1IM9 (HLA-C*04:01:QYDDAVYKL:KIR2DL1), Tyr84:OH and CE1 lie 3.304 A and 3.242 A from KIR2DL1 Asp183:OD2 (8 heavy-atom pairs < 4.5 A; Moradi et al. 2021 PMID:33846289). "
                                "Gate 0B-SCT therefore tests Y84A side-by-side against WT-Y84 with a 25-aa (G4S)5 L1 linker, DT-1 HLA-C(Y84C, A139-WT) + L1-Cys2 (Truscott et al. 2007 PMID:17475856), and DT-2 HLA-C(Y84C/A139C) + (G4S)4 L1 (Finton et al. 2023 PMID:37207206) via monomeric KIR2DL1/2/3 SPR (requiring K_D <= 10.29 uM, i.e., >=70% of native K_D = 7.2 uM, with <=15.0 uM dimeric KIR-Fc screen)."
                            ),
                        },
                        "kir_group_only_4digit_hlac_mismatch_exposure": {
                            "HLA_C_03_04_C1_af_0p09": {"self_4digit_match_pct": 17.2, "foreign_4digit_mismatch_pct": 82.8},
                            "HLA_C_01_02_C1_af_0p05": {"self_4digit_match_pct": 9.8, "foreign_4digit_mismatch_pct": 90.2},
                            "HLA_C_04_01_C2_af_0p11": {"self_4digit_match_pct": 20.8, "foreign_4digit_mismatch_pct": 79.2},
                        },
                        "immunological_mechanism": (
                            "KIR2DL2/3 (C1) and KIR2DL1 (C2) recognize the dimorphic position-80 Asn80/Lys80 epitope shared across HLA-C alleles, "
                            "whereas host CD8+ T-cell alpha/beta TCRs and anti-HLA-C DSAs recognize polymorphic 4-digit HLA-C heavy chains (Tran et al. 2011 PMID:21743387; Xu et al. 2019 PMID:30853558). "
                            "Thus a 2-line C1/C2 sub-bank matched only by KIR group presents a foreign 4-digit HLA-C alloantigen to 79.2-90.2% of routed recipients, "
                            "and presents endogenous beta-cell autoantigens (PPI, IGRP, GAD65, ZNT8) to HLA-C-restricted autoreactive memory CD8+ T cells in the remaining 9.8-20.8% 4-digit self-matched T1D recipients. "
                            "Promoted B'-SCT (ISLET1_B_SCT_C1_CIITA / ISLET1_B_SCT_C2_CIITA) resolves endogenous autoantigen presentation and avoids the 1.41-Mb chr6p21.3 HLA-A/HLA-B cis-deletion hazard by knocking out B2M (Locus 1: Arm 1B/1X) "
                            "and inserting a synthetic single-chain B2M-(G4S)4-HLA-C*03:04(GAVDPLLAL, Y84A) or C*04:01(QYDDAVYKL, Y84A) SCT loaded with a single covalently tethered peptide (human KPNA2 self-peptide GAVDPLLAL on C1; HLA-Cw4 crystallographic consensus peptide QYDDAVYKL paired with an endogenous HLA-C*04:01 self-peptide challenger on C2) under constitutive EEF1A1(5'TOP-CAGAGA) + SPA+MAZ4 DEXI insulation at CIITA Exon 5 (Locus 2), preserving K_loci=4, B=0, P_passes=2."
                        ),
                        "pass_aware_editing_complexity_tuple_muse_r3": {
                            "A_star_Champion": {"K_loci": 4, "B_base_edits": 0, "P_passes": 2, "concurrent_DSB_max": 2, "complexity_discount": 1.0000},
                            "Design_X_Molecular": {"K_loci": 4, "B_base_edits": 0, "P_passes": 2, "concurrent_DSB_max": 2, "complexity_discount": 1.0000},
                            "B_prime_SCT_primary": {"K_loci": 4, "B_base_edits": 0, "P_passes": 2, "concurrent_DSB_max": 2, "complexity_discount": 1.0000},
                            "B_prime_SCT_plus_TYK2_PE": {"K_loci": 4, "B_base_edits": 1, "P_passes": 3, "concurrent_DSB_max": 2, "complexity_discount": 0.9650},
                            "B2_endogenous_HLAC_5locus": {"K_loci": 5, "B_base_edits": 1, "P_passes": 3, "concurrent_DSB_max": 2, "complexity_discount": round(engineering_complexity_discount(5, 1), 6)},
                        },
                        "design_b_prime_sct_k4_b0_counterfactual_projection_argon_r4": {
                            "k_loci": 4,
                            "b_base_edits": 0,
                            "p_passes": 2,
                            "engineering_complexity_discount": 1.0000,
                            "gate1_differentiation_yield_tyk2_wt_restored": 0.7200,
                            "universal_population_composite_utility_projection": a_comp_disp,
                            "kir_licensed_nkg2a_low_subpopulation_composite_utility_projection": round(a_comp_unrounded * (0.6700 / a_g3_disp), 4),
                            "note": (
                                "Restores TYK2_WT Gate-1 differentiation yield (eliminating B1's -0.320 TYK2_KO penalty) and K=4, B=0 engineering_complexity_discount = 1.0000 "
                                "for ISLET1_B_SCT_C1_CIITA / ISLET1_B_SCT_C2_CIITA (conditional on passing Gate 0B-SCT monomeric KIR2DL1/2/3 SPR K_D <= 10.29 uM [<=15.0 uM dimeric KIR-Fc screen]; AUDIT-REF & AUDIT-REF)."
                            ),
                        },
                        "stage0_gate_0b_allele_pass_criteria": (
                            "Gate 0B-allele & Gate 0B-SCT: monomeric KIR2DL2/3 and KIR2DL1 SPR K_D <= 10.29 uM (>= 70% of native K_D = 7.2 uM; <= 15.0 uM dimeric KIR-Fc screen), "
                            "4-h CD8+ CTL lysis <= 10% AND 72-h MLR CD8+ CFSE-low <= 5% across n >= 6 KIR-ligand-matched PBMC donors (and across all recipients on B'-SCT), "
                            "negative recipient anti-HLA-C CDC/ADCC serum crossmatch against the exact B'-SCT product, "
                            "primary sorted KIR2DS1+/KIR2DL1- NK non-activation (CD107a <= 10%, 4-h lysis <= 10%) on all C2 arms, "
                            "AND <= 10% lysis against HLA-C-matched T1D autoreactive CD8+ T cells (Rule S5: if > 10%, demote endogenous HLA-C retention to exploratory and advance B'-SCT)."
                        ),
                    },
                    "b_prime_stage0_isogenic_2x2_satellite_arm_sonnet55_r2": {
                        "design_id": "STAGE0_BPRIME_ISOGENIC_2X2_LEADER_BY_HLAE_SCT_SATELLITE",
                        "total_clones": 12,
                        "clones_per_cell": 3,
                        "architecture": "Isogenic 2x2 satellite on the Design-A seed line using the chr6-free single-chain B2M-G4S4-HLA-C knock-in at CIITA Exon 5 (avoiding the 1.41 Mb chr6p21.3 cis-deletion hazard)",
                        "factors": {
                            "Factor_1_HLA_C_Leader_Class": "SP-1C functional leader VMAPRTLIL (HLA-C*03:04 C1 / C*04:01 C2; CD94/NKG2A K_D = 7.35 uM) vs SP-3C non-functional leader VMAPRALLL (HLA-C*07:02 C1, P6=Ala)",
                            "Factor_2_Exogenous_HLA_E_SCT": "Exogenous AAVS1::HLA-E*01:03(VMAPRTLVL) SCT cloak Present (Design B-prime) vs Absent (uncloaked Plan B)",
                        },
                        "isogenic_leader_swap_c1_c2_kir_deconfounding_control_kimi31_r3": (
                            "Because comparing native HLA-C*03:04 (SP-1C VMAPRTLIL, C1 Asn80) against native HLA-C*07:02 (SP-3C VMAPRALLL, C1 Asn80) "
                            "or HLA-C*04:01 (SP-1C VMAPRTLIL, C2 Lys80) is an allele+leader composite contrast across heavy-chain backbones and C1/C2 KIR2DL2/3 vs KIR2DL1 "
                            "licensing groups, the satellite includes a sequence-identical heavy-chain leader-swap control on the common HLA-C*03:04 (C1) backbone "
                            "(SP-1C[VMAPRTLIL]-HLA-C*03:04 vs chimeric SP-3C[VMAPRALLL]-HLA-C*03:04, alongside native HLA-C*07:02 [C1] and HLA-C*04:01 [C2] reference controls) "
                            "assayed across FACS-sorted NKG2A+KIR-, NKG2A-KIR2DL2/3+ (C1-licensed), NKG2A-KIR2DL1+ (C2-licensed), and HCMV+ NKG2C+ primary NK subsets "
                            "plus Hour-0 and 48h IFN-g HLA-E 3D12 QuantiBRITE density (AUDIT-REF)."
                        ),
                        "co_primary_readouts_and_kill_gates": [
                            "QuantiBRITE absolute surface HLA-E density (must clear >= 15,000 molecules/cell in SP-1C without SCT; if SP-1C without SCT falls < 15,000 molecules/cell, uncloaked Plan B is terminated)",
                            "Primary human NK lysis across KIR-ligand-typed donors separated into NKG2A+KIR-, C1-licensed KIR2DL2/3+NKG2A-, C2-licensed KIR2DL1+NKG2A-, and HCMV+ NKG2C+ subsets (SP-3C without SCT is expected to exhibit reduced natural signal-peptide surface HLA-E processing [Lin et al. 2023 PMID:37264229 Fig. 1c, even though exogenous peptide pulsing can stabilize surface HLA-E in Fig. 1a] and/or P6-Ala-impaired CD94/NKG2A engagement [Fig. 1-2], losing NKG2A+KIR- protection relative to SP-1C; if SP-3C without SCT matches SP-1C on both QuantiBRITE density and NKG2A+KIR- protection, the leader-class premise is falsified; AUDIT-REF, AUDIT-REF & AUDIT-REF)",
                        ],
                    },
                    "graveyard_005_and_010_counterfactual_sensitivity_sonnet55_r2": {
                        "computation_mode": "LIVE_EVALUATOR_DERIVED_KIMI31_R3",
                        "canonical_chassis_b_index": round(b_comp_unrounded, 4),
                        "canonical_chassis_b_composite_index": round(b_comp_unrounded, 4),
                        "canonical_chassis_b_pediatric_index": round(b_ped_unrounded, 4),
                        "counterfactual_with_graveyard_010_indirect_apc_minus_0p25_gate2b_logit": round(b_ped_gy10, 4),
                        "counterfactual_composite_with_graveyard_010_minus_0p25_gate2b_logit": round(b_comp_gy10, 4),
                        "counterfactual_pediatric_with_graveyard_010_minus_0p25_gate2b_logit": round(b_ped_gy10, 4),
                        "counterfactual_with_both_gy010_minus_0p25_and_uncloaked_sp3c_gy005_minus_0p20_gate3_logit": round(b_ped_gy10_m020, 4),
                        "legacy_linear_approx_counterfactual_with_both_gy010_and_gy005_reference": 0.0574,
                        "counterfactual_linear_approx_with_both_gy010_and_gy005": 0.0574,
                        "exact_composite_with_minus_0p15_gate3_logit": round(b_comp_gy10_m015, 4),
                        "exact_pediatric_with_minus_0p15_gate3_logit": round(b_ped_gy10_m015, 4),
                        "exact_composite_trajectory_50_50_ped_lada": [
                            round(b_comp_unrounded, 4),
                            round(b_comp_gy10, 4),
                            round(b_comp_gy10_m015, 4),
                        ],
                        "exact_pediatric_trajectory": [
                            round(b_ped_unrounded, 4),
                            round(b_ped_gy10, 4),
                            round(b_ped_gy10_m015, 4),
                        ],
                        "exact_composite_with_minus_0p20_gate3_logit": round(b_comp_gy10_m020, 4),
                        "exact_pediatric_with_minus_0p20_gate3_logit": round(b_ped_gy10_m020, 4),
                        "design_b_prime_k5_b1_counterfactual_projection_sonnet55_r3": {
                            "status": "HYPOTHETICAL_STAGE0_CHALLENGER_COUNTERFACTUAL_PROJECTION",
                            "k_cuts": 5,
                            "b_base_edits": 1,
                            "engineering_complexity_discount_k5_b1": round(engineering_complexity_discount(5, 1), 6),
                            "engineering_complexity_discount_k6_b1_reference": round(engineering_complexity_discount(6, 1), 6),
                            "complexity_discount_ratio_k5_over_k6": round(engineering_complexity_discount(5, 1) / engineering_complexity_discount(6, 1), 6),
                            "projected_composite_utility_at_k5_b1": round(b_comp_unrounded * (engineering_complexity_discount(5, 1) / engineering_complexity_discount(6, 1)), 4),
                            "projected_pediatric_utility_at_k5_b1": round(b_ped_unrounded * (engineering_complexity_discount(5, 1) / engineering_complexity_discount(6, 1)), 4),
                            "projected_composite_with_gy010_at_k5_b1": round(b_comp_gy10 * (engineering_complexity_discount(5, 1) / engineering_complexity_discount(6, 1)), 4),
                            "projected_pediatric_with_gy010_at_k5_b1": round(b_ped_gy10 * (engineering_complexity_discount(5, 1) / engineering_complexity_discount(6, 1)), 4),
                            "projected_composite_auto_none_tyk2_ko_gate1_retained_at_k5_b1": 0.0569,
                            "projected_pediatric_auto_none_tyk2_ko_gate1_retained_at_k5_b1": 0.0556,
                            "projected_composite_auto_none_tyk2_p1104a_wt_gate1_yield_at_k5_b1": 0.0629,
                            "projected_pediatric_auto_none_tyk2_p1104a_wt_gate1_yield_at_k5_b1": 0.0615,
                            "note": (
                                "Separates the grid-evaluated K=6, B=1 Chassis B (ISLET1_B1_6LOCUS_TYK2_KO_INS, complexity discount exp(-0.15*2)*0.965 = 0.714890 [0.7148895830 exact], composite = 0.0673, pediatric = 0.0664) "
                                "from the K=5, B=1 Design B-prime edit-budget sensitivity counterfactual (ISLET1_B2_5LOCUS_TYK2_PE [5 DSB loci: cis-linked HLA-A/HLA-B + CIITA + CXCL10 + AAVS1 + 1 TYK2 P1104A prime edit; omitting B1's INS leader/DRiP base edit to isolate a single B=1 edit, whereas retaining both TYK2 P1104A PE and INS base edit would be K=5, B=2 with discount 0.801514] or chr6-free CIITA::B2M-G4S4-HLA-C knock-in [K=4] + 5th DSB locus + TYK2 P1104A PE, "
                                "complexity discount exp(-0.15*1)*0.965 = 0.830583 [0.8305831973 exact]: holding all 6 gate logits fixed yields a pure 1.161834x complexity-discount sensitivity of 0.0782 composite / 0.0771 pediatric; dropping the +0.25 Gate-2B INS-edit credit [equivalent to symmetric GRAVEYARD_010 -0.25 Gate-2B logit shift] yields 0.0724 composite / 0.0711 pediatric; and re-evaluating the full AUTO_NONE Gate-2B branch [-0.65 pediatric / -0.50 LADA Gate-2B logit shift relative to AUTO_DUAL_INS_LEADER_AND_DRIP_EDIT] yields 0.0569 composite / 0.0556 pediatric if the -0.32 TYK2-KO Gate-1 penalty is retained or 0.0629 composite / 0.0615 pediatric if hypomorphic TYK2 P1104A preserves WT Gate-1 differentiation yield; AUDIT-REF, AUDIT-REF, L1B-F2 & L1B-COD-01 [DEC-077/DEC-078])."
                            ),
                        },
                        "note": (
                            "Documents symmetric Graveyard counterfactuals for Slot B (AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF): while GRAVEYARD_010's warning rule triggers on 2 B2M-null modules (ALLO_B2M_HLAE_CD47_SWAP_CIITA_KO, ALLO_B2M_CIITA_KO_CD47_LENTI) and the -0.25 Gate-2B logit discount applies across all 3 B2M-null modules (adding ALLO_B2M_KO_PDL1_A20_KI) in the Gate-2B evaluator branch, host APCs still present shed donor peptides via indirect MHC-II "
                            f"even when HLA-C is retained in CIITA-/- Chassis B; applying GRAVEYARD_010's -0.25 Gate-2B logit discount symmetrically lowers Chassis B from {b_comp_unrounded:.4f} composite ({b_comp_unrounded:.6f} unrounded) / {b_ped_unrounded:.4f} pediatric ({b_ped_unrounded:.6f} unrounded) to {b_comp_gy10:.4f} composite ({b_comp_gy10:.6f} unrounded) / {b_ped_gy10:.4f} pediatric ({b_ped_gy10:.6f} unrounded, -7.8%), "
                            f"and additionally applying an uncloaked SP-3C Gate-3 leader-deficit discount lowers it to {b_comp_gy10_m015:.4f} composite ({b_comp_gy10_m015:.6f} unrounded) / {b_ped_gy10_m015:.4f} pediatric ({b_ped_gy10_m015:.6f} unrounded) at -0.15 logit ({b_comp_gy10_m020:.4f} composite [{b_comp_gy10_m020:.6f} unrounded] / {b_ped_gy10_m020:.4f} pediatric [{b_ped_gy10_m020:.6f} unrounded] at -0.20 logit; 0.0574 linear approximation), leaving rank order A >> B unchanged."
                        ),
                    },
                },
                "X_Locked_Number1_Evolutionary_Contender": {
                    "slot": "X",
                    "role": "Locked #1 Evolutionary Contender built to dethrone Champion A head-to-head in Stage 0 (never overwritten by unvetted speculative ideas)",
                    "architecture_summary": (
                        "K=4, B=0 (RNLS-WT / TXNIP-WT preserving full +3.39 mmHg O2 headroom): tests DX-S Split-Cassette "
                        "(Arm 1C promoter-trapped B2M::HLA-E SCT [3,459 bp donor / 3,749 bp packaged] + relocated AAVS1::cHS4-A2UCOE-EEF1A1-iCasp9-GSG-P2A-CD47-IgV-GPI [4,048 bp donor / 4,338 bp packaged, with CD47-IgV-GPI strictly downstream of P2A so CD55 GPI signal remains C-terminal; Q3-P1 & AUDIT-REF]) "
                        "and Arm 1A-lite (3,891 bp donor / 4,181 bp packaged; Q1-P1) against Arm 1X (4,056 bp promoter-only core / 4,225 bp fully loaded donor with 154-bp non-palindromic SPA+spacer+MAZ4 sense terminator [Q1-P8], "
                        "4-nt S3 AsCas12a PAM+seed shield [Q1-P5], and 20-aa L1=(G4S)4 linker [Q1-P4]), paired with Scar-Free Dual-Promoter Locus 3 (3,717 bp unspaced / 3,799 bp donor with 82-bp inter-unit pause spacer / 4,089 bp packaged) "
                        "and DX-b 100 um thinned PDMS-CaO2 lamina (eps = 0.42) with <=100 um (R <= 50 um) micro-islet cluster diameter cap (Q3-P3, Q4-P4)"
                    ),
                    "b2m_null_orphan_hc_erad_status_muse_r3": "B2M_NULL_ORPHAN_HC_ERAD_UNQUALIFIED",
                    "dxs_single_locus_shield_plus_switch_correlation_pricing_muse_r3": {
                        "product_variant_id": "ISLET1_A4_DXS_ARM1C_L4CD47",
                        "mechanism": (
                            "In DX-S, relocating CD47-IgV-GPI into AAVS1 (Locus 4: cHS4-A2UCOE-EEF1A1 -> iCasp9-GSG-P2A-CD47-IgV-GPI) couples "
                            "macrophage stealth (CD47) and the suicide switch (iCasp9) onto the same single-promoter bicistronic transcript at Chr19q13.42. "
                            "Whereas total 19q13.42 LOH or complete epigenetic silencing of AAVS1 is self-policing under intact host macrophage surveillance "
                            "(CD47-null cells lose SIRPa inhibition and are cleared by macrophages), partial silencing or 5'-selective iCasp9 frameshift/point mutation "
                            "can impair AP1903 killing while retaining sufficient CD47 to survive, or focal AAVS1 silencing during transient macrophage depletion can cause correlated loss of both."
                        ),
                        "lambda_couple_prior_range": [0.0, 0.05],
                        "lambda_couple_central_prior": 0.025,
                        "correlation_multiplier_1_minus_lambda_couple": 0.9750,
                        "unpenalized_dxs_composite_projection": x_comp_disp,
                        "correlation_priced_dxs_composite_projection_at_lambda_0p025": round(x_comp_unrounded * (1.0 - 0.025), 4),
                        "stage0_joint_retention_assay_gate": (
                            "60-day differentiation + 14-day cytokine/5-aza-dC challenge dual-color flow & AP1903 kill assay requiring "
                            "joint (CD47-low AND iCasp9-low) double-loss fraction <= 10^-6 and post-AP1903 surviving CD47-high fraction <= 10^-8 (Rule S7; AUDIT-REF)."
                        ),
                    },
                    "composite_4way_research_lead_specification_muse_r4": {
                        "composite_lead_id": "ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500",
                        "status": "PRE_REGISTERED_STAGE0_COMPOSITE_RESEARCH_LEAD_PENDING_5_GATE_COUPON_QUALIFICATION",
                        "locus_1_b2m": (
                            "Arm 1X (unmaterialized W_P=+1 2^3 factorial challenger: hybrid constitutive-floor Syn_3xISRE_CpGminP_230bp_sense [4,056 bp promoter-only core / 4,225 bp fully loaded donor, with optional constitutive EEF1A1_core_230(5'TOP-CAGAGA) sense comparator] -> HLA-E(WT D227/T228 conditional default vs D227K/T228A) + antisense CD47(C15G)-IgV-GPI "
                            "+ 154-bp SPA+spacer+MAZ4 sense terminator + 4-nt S3 PAM+seed shield + 20-aa (G4S)4 L1 linker; represented in the product registry by the materialized 4,143-bp Embodiment A* baseline under RNLS+/+ lock [ISLET1_A4_X_MOL_RNLS_WT, carrying sense Syn_5xISRE_minP_148bp = W_P=-1] and evaluated head-to-head against defended champion Embodiment A* "
                            "4,143 bp [ISLET1_A4_ARM1B / ISLET1_A_STAR_4LOCUS] across sorted pure CD56dim NKG2A- NKG2C+ adaptive NK subsets [partitioned into KIR2DL1/2/3+ vs KIR2DL1/2/3-; <= 8% 4-h lysis, CD107a <= 10%] "
                            "and HLA-E-restricted CMV UL40 CD8+ CTLs comparing WT HLA-E vs CD8-null HLA-E(D227K/T228A); noting per Kaiser et al. 2008 PNAS PMID:18448674 Table 1 that soluble non-covalent HLA-E*01:03/VMAPRTLFL "
                            "binds CD94/NKG2A with K_D = 0.80 uM and CD94/NKG2C with K_D = 5.20 uM [K_C / K_A = 6.50x; Vales-Gomez et al. 1999 PMID:10428963 Table I reports 2.33 uM vs 10.77 uM = 4.62x; vs Lin et al. 2023 PMID:37264229 VMAPRTLVL 2.10 uM vs 6.52 uM = 3.10x] and is therefore an NKG2A-high-affinity comparator that remains an absolute NKG2C agonist at elevated surface density; AUDIT-REF, AUDIT-REF, R9-S2-L2-03 & L1B-NUM-01 [DEC-044/DEC-045/DEC-073/DEC-074/DEC-078])"
                        ),
                        "locus_2_ciita": "Biallelic CIITA Exon-1/Exon-5 KO (OI-GRNA-02 AsCas12a Ultra in Pass 1 alongside B2M to prevent inter-chromosomal translocations with Pass-2 HiFi-SpCas9; GRAVEYARD_006) in Plan A/X; or CIITA Exon-5::HLA-C*03:04(C1, Asn80, Y84A, GAVDPLLAL) / HLA-C*04:01(C2, Lys80, Y84A, QYDDAVYKL) 3,465-bp SCT knock-in in Plan B'-SCT (ISLET1_B_SCT_C1_CIITA / ISLET1_B_SCT_C2_CIITA; never the GRAVEYARD_005 non-functional SP-3C VMAPRALLL HLA-C*07:02 trap allele; with mandatory 4-digit HLA-C typing, anti-HLA-C DSA, KIR2DS1 counter-screen, and recipient CD8+ MLC gates; AUDIT-REF)",
                        "locus_3_cxcl10": "L3-BF-B 4,078-bp Scar-Free Hybrid Basal-Floor + Inducible Dual-Promoter (ISLET1_A4_L3_BFB_BASALFLOOR: 230-bp Syn_3xISRE_2xNFkB_hPGK1_Sp1_basalfloor driving scar-free 376-aa SERPINB9 [0.30 uM Hour-0 basal floor, 0.89 uM 24-h IFN-gamma surge; 6-pulse autonomous intracellular buffer under f_ISR = 0.45 paired with upstream synaptic degranulation suppression p_degran <= 0.15 for 12-pulse survival; AUDIT-REF] + 82-bp pause spacer + 148-bp Syn_5xISRE_minP driving native-Met1 290-aa PD-L1)",
                        "locus_4_aavs1": "Canonical independent AAVS1::cHS4-A2UCOE-EEF1A1-iCasp9 (3,457 bp donor, preserving locus independence from Chr15 CD47 unless DX-S clears the N_eff >= 299,573,228 effective-clonogenic-cell Rule S7 gate, and paired with primary mechanical containment + en-bloc surgical omental pouch retrievability under Gate 0S-Retrieval; AUDIT-REF, AUDIT-REF & AUDIT-REF)",
                        "scaffold_and_dosing": "DX-B-LowLoad-1500: 100 um thinned central PDMS-CaO2 lamina (eps = 0.42, 46.39 wt% CaO2, >= 32,494 IU localized catalase), D_O2,gel >= 1.50e-5 cm^2/s, <= 100 um (R <= 50 um; <= 80 um / R <= 40 um chronic fibrotic-capsule resilience target for capsule_dz up to 100 um; AUDIT-REF) micro-islet cluster cap, and 1,500 IEQ/cm^2 chronic areal loading across two 166.67 cm^2 portal-draining omental pouches (333.33 cm^2 total for 500,000 IEQ)",
                        "stage0_combined_stress_coupon_gate": "STAGE0_COMBINED_STRESS_KILL_CHAIN_COUPON_MUSE_R4",
                        "indirect_apc_hip_cd4_gate_0d_sol61_r5": {
                            "finding_id": "AUDIT-REF",
                            "decision_id": "DEC-044 / DEC-045 / DEC-046 / DEC-078",
                            "primary_citations": ["PMID:26912858", "PMID:35559682", "PMID:27872147", "PMID:39554513", "PMID:2124002", "PMID:1986378", "PMID:1512546", "PMID:8120026"],
                            "mechanism": (
                                "While B2M-/- and CIITA-/- eliminate direct MHC-I and MHC-II presentation on graft cells, every metabolically active beta-cell graft continuously sheds "
                                "insulin secretory granules, exosomes (Cianciaruso et al. 2017 Diabetes PMID:27872147), and apoptotic bodies containing Hybrid Insulin Peptides (HIPs: insulin C-peptide fragments covalently cross-linked to ChgA WE14 or IAPP2; "
                                "Delong et al. 2016 Science PMID:26912858; HLA-DQ8 presentation characterized by Bhattacharjee et al. 2024 PNAS Nexus PMID:39554513), and transgenic SCT linker-junction neo-peptides (plus potential exogenous serum beta2-microglobulin re-association with orphan endogenous HLA-A/B/C heavy chains in B2M-/- cells, as established for beta2m-free MHC-I heavy chains by Vitiello et al. 1990 Science PMID:2124002, Rock et al. 1991 PNAS PMID:1986378, Bix & Raulet 1992 J Exp Med PMID:1512546, and Demaria et al. 1994 J Biol Chem PMID:8120026) that are endocytosed by recipient HLA-DQ2 (DQA1*05:01/DQB1*02:01), HLA-DQ8 (DQA1*03:01/DQB1*03:02), and HLA-DR4+ dendritic cells and macrophages "
                                "for indirect CD4+ T-cell activation, triggering localized TNF-alpha/IFN-gamma/IL-1beta/NO/ROS bystander toxicity that bypasses MHC-null immune cloaks."
                            ),
                            "stage0_gate_0d_indirect_apc_co_primary_assay": (
                                "Freeze-thawed / cytokine-stressed SC-islet lysate, isolated Stage-7 SC-islet secretory exosomes (Cianciaruso et al. 2017 PMID:27872147), and synthetic HLA-E/HLA-C SCT (G4S)3-4 linker-junction peptides "
                                "(with semi-direct intact-exosome MHC cross-dressing controls and a project-defined Stage-0 Gate 0d 2.0 ug/mL exogenous human beta2m W6/32 + CD8+ CTL challenge grounded in Vitiello et al. 1990 PMID:2124002, Rock et al. 1991 PMID:1986378, Bix & Raulet 1992 PMID:1512546, and Demaria et al. 1994 PMID:8120026, alongside HIP-DQ8 controls per Bhattacharjee et al. 2024 PMID:39554513) co-cultured with primary HLA-DQ2+, HLA-DQ8+, and HLA-DR4+ human monocyte-derived dendritic cells and macrophages "
                                "and insulin-B:9-23 / HIP-reactive and recipient polyclonal CD4+ T-cell responder lines in the presence vs absence of SA-FasL (100-150 ng/uL) + SA-THBD microgels (Lei et al. 2022 Sci Adv PMID:35559682), "
                                "requiring CD4+ T-cell proliferation Stimulation Index SI <= 1.5, IFN-gamma ELISpot/ELISA SI <= 2.0, >= 80% FasL-mediated apoptosis of activated CD25+CD95+ effector CD4+ T cells, "
                                "transwell bystander SC-beta Annexin-V+ apoptosis <= 15%, and full-evidence Stage-0 promotion gating (AUDIT-REF, AUDIT-REF, AUDIT-REF & L2B-NUM-02 [DEC-044/DEC-045/DEC-046/DEC-078])."
                            ),
                        },
                    },
                },
            },
            "tier2_radical_outside_the_box_wildcards_20pct_focus": {
                "Y_Radical_ZeroDSB_and_Active_Indirect_Cloak_Wildcard": {
                    "slot": "Y",
                    "role": "Radical Biology & Genome Engineering Wildcard competing to dethrone X as #1 Contender",
                    "architecture_summary": (
                        "Non-exposed-DSB genome engineering (single-landing-pad PASTE / TwinPE + Bxb1 serine integrase at B2M/AAVS1 + heritable dCas9-KRAB-DNMT3A/3L CRISPRoff "
                        "epigenetic silencing of CIITA/CXCL10; avoiding intended exposed-DSB editing and AAV6 donor delivery while directly enumerating residual nick-associated DSB, "
                        "integrase pseudo-site insertion, and off-target methylation risks; AUDIT-REF) + active indirect-APC & classical/terminal complement "
                        "shield (separate native CD55 + CD59 transcription units or membrane-bound CD55-CD59 chimera + membrane CD200 or sCTLA4-Ig)"
                    ),
                    "stage0_quantitative_kill_gates_muse_r1": [
                        "Y-01 (Epigenetic Mitotic-Silence Retention Gate): CRISPRoff silencing of B2M, CIITA, and CXCL10 must maintain <=1% single-cell HLA-I/HLA-II/CXCL10 functional protein escape and >=99% promoter CpG methylation across 35-day SC-islet differentiation plus 14-day IFN-gamma + TNF-alpha + IL-1beta + hypoxia re-challenge across all product lineages; any locus with >1% escape is rejected",
                        "Y-02 (TwinPE/PASTE Full-Cargo Landing-Pad Fidelity Gate): Bxb1 cargo integration must achieve >=25% seamless integration of the full-length 7-9 kb poly-cargo (not merely attB installation) prior to clonal enrichment, with 100% single-copy sequence verification in the selected clone, zero concatemer/backbone/pseudo-att-site insertions, CAST-Seq translocation frequency <10^-5, and zero TP53/Chr20q selection",
                        "Y-03 (CD39/CD73 Paracrine Adenosine GSIS-Safety Kill-Gate): Because adenosine A1/A3 receptor activation on beta/delta cells inhibits adenylyl cyclase and suppresses GSIS, CD39/CD73 is excluded from the default core and any tested ectonucleotidase comparator that lowers perifusion SI by >10% or elevates basal insulin by >15% is immediately disqualified in favor of separate CD55 + CD59 units plus membrane CD200/sCTLA4-Ig",
                    ],
                },
                "Z_Radical_ZeroCaO2_Microvascular_Ribbon_Wildcard": {
                    "slot": "Z",
                    "role": "Radical Biophysical & Vascular Niche Wildcard competing to dethrone X as #1 Contender",
                    "architecture_summary": (
                        "Zero-CaO2, zero-H2O2, zero-alkali pleated microvascular ribbon (<=120 um total thickness, L_half = 60 um <= passive Krogh limit) "
                        "housing <=80 um tri-hormonal micro-islets with 10% v/v perfluorooctyl bromide (PFOB) passive oxygen-shuttle nanoemulsion "
                        "(20x aqueous O2 solubility without chemical hydrolysis, equilibrium storage capacity alpha_eff/alpha_aq = 2.90x, steady-state partial-pressure Maxwell-Garnett permeability P_eff/P_aq = 1.2836x at D_PFC=D_aq [bounded above by 4/3 = 1.3333x as r->inf vs 6/7 = 0.8571x for impermeable spheres], apparent concentration diffusivity D_app/D_aq = 0.4426x; AUDIT-REF) and pre-inosculated iEC:iPericyte vascular cords (target t_50 <= 24 h)"
                    ),
                    "stage0_quantitative_kill_gates_muse_r1": [
                        "Z-01 (Permeation, Storage Capacity & Emulsion Stability Coupon): Steady-state partial-pressure Maxwell-Garnett permeability P_eff/P_aq in [1.20x, 1.33x] (1.2836x at D_PFC=D_aq, phi=0.10, K_p=20) with equilibrium O2 storage capacity alpha_eff/alpha_aq = 2.90x and DLS droplet Z-average <= 220 nm with <10% Ostwald ripening over 14 days at 37 C (AUDIT-REF)",
                        "Z-02 (72-h Hypoxic Chamber Survival Without CaO2): Under pO2_ext = 25 mmHg and 0 mg CaO2, 80 um micro-islet core pO2 >= 5.0 mmHg, cell viability >= 90%, and perifusion SI >= 2.5",
                        "Z-03 (PFOB Innate Macrophage & Complement Screen): Primary human macrophage IL-1beta/TNF-alpha and serum C3a/C5a generation <=1.5x blank clinical-grade alginate/PEG hydrogel",
                        "Z-04 (Microfluidic Inosculation Kinetics): Perfusable iEC + pericyte microvascular lumen anastomosis within t_50 <= 24-36 h",
                    ],
                },
            },
            "q5_dual_clock_governance_and_screen_confirm_firewall_sonnet55_r1": {
                "d1_d2_screen_confirm_firewall_q5_p1": (
                    "Stage-0 / Stage-1 HIP-Core (Runs 01-36 / 01-40) is strictly a Discovery Screen (Dataset D1) used to select a single locked "
                    "winner genotype G*; promotion to IND-enabling GLP studies requires an independent Confirmatory Cohort (Dataset D2: n >= 6 new "
                    "independent differentiation batches x n >= 6 unseen HLA/KIR/SIRPA-typed human immune donors, exact one-sided binomial/t gate "
                    "with zero post-hoc tuning; AUDIT-REF)."
                ),
                "quantitative_tournament_promotion_and_budget_firewall_rules_muse_r1": {
                    "promotion_posterior_superiority_rule": (
                        "A contender (X displacing A/B, or Y/Z displacing X) must demonstrate a pre-registered Bayesian posterior superiority probability "
                        "P(Delta_primary >= +5.0 pp | D2) > 0.95 on its target bottleneck gate on the unseen Confirmatory Cohort D2 while satisfying "
                        "a strict -3.0 pp non-inferiority margin across all other co-primary gates (GSIS SI, core pO2, HDR yield, genomic safety, and iCasp9 >=99.9% lethality)."
                    ),
                    "quarterly_slot_freeze_cadence": "Wildcard slots Y and Z are frozen on a quarterly review cadence; literature-only proposals cannot overwrite X without wet-lab Stage-0 kill-coupon clearance.",
                    "wildcard_budget_and_fte_cap": "Strict <=20% bench-capacity and experimental-run firewall on Tier-2 wildcards (Y/Z) so >=80% of Stage-0 bench capacity remains locked on the A/B/X execution triad.",
                    "canine_track_stage_gate_dependency": "No in vivo client-owned diabetic dog study proceeds until (1) human D1 winner lock G* is achieved and (2) the in vitro Canine Homology Gate (DLA-88/CD94-NKG2A + canine CD47/SIRPa KD <= 2.0 uM + <=10% canine PBMC killing) passes.",
                },
                "dual_clock_timeline_and_budget_envelope_q5_p2": {
                    "optimistic_parallel_fast_track_clock": "2026-2030 (parallelized Stage 0-1 bench and validation critical path)",
                    "conservative_sequential_risk_adjusted_clock": "54-78 months (sequential Stage 0-1 validation through Phase I/IIa FIH)",
                },
                "canine_dla88_cd47_sirpa_homology_gate_q5_p4": (
                    "Before initiating the client-owned diabetic dog veterinary track, require an in vitro Canine Homology Gate: "
                    "(1) SPR + flow verification that canine DLA-88 leader/SCT engages canine CD94/NKG2A and canine CD47-IgV-GPI binds canine SIRPa "
                    "with KD <= 2.0 uM (and zero human SIRPa xeno-incompatibility), and (2) canine PBMC/NK/macrophage killing <=10%."
                ),
            },
        },
    }

    return champion_challenger_portfolio


def main() -> None:
    graph, studies, omics = load_inputs()
    target_profiles = build_target_evidence_profiles(studies, omics)
    all_combos = run_combinatorial_multi_gate_search(target_profiles)
    combo_by_id = {c["combination_id"]: c for c in all_combos}

    # Build Pareto frontier by chromosomal locus budget K strictly over feasible (non-BLOCKed) stacks, and also track clinically switched B<=1 candidates (AUDIT-REF; DEC-049):
    pareto_by_k = {}
    clinical_switched_pareto_by_k = {}
    for c in all_combos:
        if not c.get("is_feasible_candidate", True):
            continue
        k = c["total_genomic_cuts_K"]
        if k not in pareto_by_k or c["composite_series_system_utility"] > pareto_by_k[k]["composite_series_system_utility"]:
            pareto_by_k[k] = {
                "total_genomic_cuts_K": k,
                "total_base_edits_B": c["total_base_edits_B"],
                "best_combination_id": c["combination_id"],
                "composite_series_system_utility": c["composite_series_system_utility"],
                "global_rank_among_all_combinations": c["global_rank_among_all_combinations"],
                "unswitched_immune_evasive_clinical_hold": c.get("unswitched_immune_evasive_clinical_hold", False),
            }
        if not c.get("unswitched_immune_evasive_clinical_hold", False) and not c.get("multiplex_base_edit_plus_dsb_translocation_risk_flag", False):
            if k not in clinical_switched_pareto_by_k or c["composite_series_system_utility"] > clinical_switched_pareto_by_k[k]["composite_series_system_utility"]:
                clinical_switched_pareto_by_k[k] = c
    for k, p_row in pareto_by_k.items():
        cs = clinical_switched_pareto_by_k.get(k)
        p_row["clinical_switched_best_combination_id"] = cs["combination_id"] if cs else None
        p_row["clinical_switched_composite_series_system_utility"] = cs["composite_series_system_utility"] if cs else None
        p_row["clinical_switched_global_rank"] = cs["global_rank_among_all_combinations"] if cs else None

    ranked_architectures = []
    for arch in CANDIDATE_ARCHITECTURES:
        cid = arch["matched_combinatorial_search_id"]
        c = combo_by_id[cid]
        if c["total_genomic_cuts_K"] != arch["genomic_locus_edit_count_K"]:
            raise ValueError(
            f"Cut count mismatch for {arch['architecture_id']}: "
            f"{c['total_genomic_cuts_K']} vs {arch['genomic_locus_edit_count_K']}"
        )
        if c["total_base_edits_B"] != arch["base_edit_count_B"]:
            raise ValueError(
            f"Base-edit count mismatch for {arch['architecture_id']}: "
            f"{c['total_base_edits_B']} vs {arch['base_edit_count_B']}"
        )
        arch_out = dict(arch)
        product_id = {"ISLET1_CHASSIS_A_UNIVERSAL_LEAN_4LOCUS": "ISLET1_A4_ARM1B",
                      "ISLET1_CHASSIS_B_HLA_C_RETAINED_6LOCUS_HYPOTHESIS": "ISLET1_B1_6LOCUS_TYK2_KO_INS"}.get(arch["architecture_id"])
        if product_id:
            contract = _load_local_helper("product_contract").get_product_definition(product_id, root=PROJECT_ROOT)
            if contract["K"] != c["total_genomic_cuts_K"] or contract["B"] != c["total_base_edits_B"]:
                raise ValueError(f"Product contract K/B mismatch for {product_id}")
            arch_out["product_definition"] = contract
        arch_out["is_feasible_candidate"] = c["is_feasible_candidate"]
        arch_out["unblocked_diagnostic_utility"] = c["unblocked_diagnostic_utility"]
        arch_out["global_rank_among_all_combinations"] = c["global_rank_among_all_combinations"]
        arch_out["utility_pediatric_early_onset_t1d"] = c["utility_pediatric_early_onset_t1d"]
        arch_out["utility_adult_onset_lada"] = c["utility_adult_onset_lada"]
        arch_out["composite_series_system_utility"] = c["composite_series_system_utility"]
        arch_out["search_grid_lower_10pct_cred_bound_pediatric"] = c["utility_pediatric_early_onset_t1d"]["lower_10pct_cred_bound"]
        arch_out["lower_10pct_cred_bound_with_b_edit_quadratic_discount_pediatric"] = c["utility_pediatric_early_onset_t1d"]["lower_10pct_cred_bound_with_b_edit_quadratic_discount"]
        arch_out["triggered_graveyard_rules"] = c["triggered_graveyard_rules"]
        ranked_architectures.append(arch_out)

    ranked_architectures.sort(
        key=lambda a: (
            a["is_feasible_candidate"],
            a["composite_series_system_utility"],
            a["unblocked_diagnostic_utility"],
        ),
        reverse=True,
    )
    for idx, a in enumerate(ranked_architectures, 1):
        a["rank_among_6_named_architectures"] = idx

    # Build the explicit 36-run HIP-Core Wet-Lab Protocol (all design algebra below is computed exactly from this run table):
    # - Tier-1: 16-Run 2^(6-2) Resolution IV Screening Block (Runs 01-16, rank 14/22)
    # - Tier-2: 20-Run Bridge + Replicated De-Aliasing + Lead-Chassis Confirmatory Block (Runs 17-36):
    #   * Runs 17-24: 8 De-Aliasing Runs (Replicate A; Run 17 is the LEGACY lean chassis with the endogenous B2M-promoter bicistronic CD47 of Locus 1, X2 = -1)
    #   * Runs 25-26: 2 Tier-1 Bridge Reference Controls (repeating Tier-1 Run 01 [-1,-1,-1,-1,-1,-1] and Run 16 [+1,+1,+1,+1,+1,+1] inside Tier-2 so rank([X | B]) = 23/23)
    #   * Runs 27-34: 8 Replicated Tier-2 De-Aliasing Runs (Replicate B, eliminating all unit-leverage h_ii = 1.0 rows)
    #   * Runs 35-36: 2 Lead-Chassis-A confirmatory replicates (current Scenario A genotype [+1,+1,+1,-1,-1,-1]: constitutive hEF1a-core CD47 package X2 = +1), which the
    #     34-run design never built (its prediction variance at that genotype was 33/7 to 47/7 sigma^2) and which roughly halve the worst-case 2FI variance multiplier
    factor_defs = {
        "X1": {
            "symbol": "A",
            "name": "Primary_MHC_I_Missing_Self_Element",
            "gate": "Gate_2A_3_Allo_and_NK_Missing_Self",
            "factor_class": "Single_Locus_Allelic_Comparison (see locus1_x1_by_x2_physical_construct_matrix for exact 4-cell monocistronic vs bicistronic physical cassette mapping)",
            "level_minus1": "B2M-/- monocistronic CD47 knock-in without HLA-E (LOC1_X1M_X2M_B2M_CD47WT_MONOCISTRONIC when X2=-1; LOC1_X1M_X2P_HEF1A_CD47GPI_MONOCISTRONIC when X2=+1; NO orphan 2A tag)",
            "level_plus1": "B2M::HLA-E*01:03(VMAPRTLVL)-GSG-tPT2A-CD47 bicistronic exon 1 knock-in (LOC1_X1P_X2M_B2M_HLAE_TPT2A_CD47WT_BICISTRONIC when X2=-1; LOC1_X1P_X2P_HEF1A_HLAE_TPT2A_CD47GPI_BICISTRONIC when X2=+1; PMID:28504668, PMID:39621955; ~12% overlap with clinical HCMV UL40 isolates in PMID:37264229; cytoplasmic KW-motif and 2A-free bidirectional cassette evaluated in Stage-0 Pre-Pilot)",
        },
        "X2": {
            "symbol": "B",
            "name": "CD47_Promoter_and_Ectodomain_Package",
            "gate": "Gate_3_Innate_NK_Macrophage_Hour0_Stoichiometry",
            "factor_class": "Composite_Engineering_Package (strictly locked to one pre-registered construct per (X1, X2) cell in locus1_x1_by_x2_physical_construct_matrix; single-component Promoter-Only vs Ectodomain-Only [CD47-IgV-GPI vs CD47-dCT vs CD47-WT] and 2A-free bidirectional controls run in Stage-0 Decision Experiment #1)",
            "level_minus1": "Endogenous resting B2M-promoter driving full-length CD47-WT (monocistronic Kozak-ATG-CD47-WT when X1=-1; bicistronic -GSG-tPT2A-CD47-WT when X1=+1; tests Graveyard 009 Hour-0 basal density deficit and full-length CD47-WT Cdc42/GSIS signaling; PMID:37156915, PMID:37820008)",
            "level_plus1": "Locked High-Floor Constitutive/Synthetic Promoter [native human EEF1A1-core (230 bp, 19 CpGs, 60.0% GC) + chimeric intron] driving GPI-anchored CD47-IgV-GPI(DAF) (monocistronic Kozak-ATG-CD47-IgV-GPI-WPREmut6 when X1=-1; bicistronic 3'-cistron -GSG-tPT2A-CD47-IgV-GPI-WPREmut6 when X1=+1; pre-registered uniform amendment rule: if Stage-0 GPI-PLC/QuantiBRITE gate fails, all 13 X2=+1 runs across the 36-run table [8 in Tier-1, 1 in Run 26 bridge, 2 in Runs 24/34, 2 in Runs 35/36] uniformly substitute CD47-dCT prior to Stage-1 initiation; Xu et al. 2024 PMID:38508139, PMID:40668892, PMID:30778232)",
        },
        "X3": {
            "symbol": "C",
            "name": "Locus3_Chemokine_KO_and_Inducible_Shield_Package",
            "gate": "Gate_2B_3_Cytokine_and_Granzyme_B_Shield",
            "factor_class": "Composite_Single_Locus_Swap_Package (locked to the single inverted-cistron construct CXCL10::SERPINB9-GSG-tPT2A-PD-L1 in the 36-run factorial; CXCL10-/- empty-KO vs PD-L1-only vs SERPINB9-only and PD-L1-GSG-tPT2A-SERPINB9 orientation controls run in the Stage-0 Pre-Pilot / Decision Experiment #3)",
            "level_minus1": "Wild-type CXCL10 locus (no PD-L1 or SERPINB9 knock-in)",
            "level_plus1": "Locked CXCL10::SERPINB9-GSG-tPT2A-PD-L1 biallelic swap-in (deletes graft CXCL10 + expresses cytosolic SERPINB9 5'-upstream of tPT2A and untagged C-terminal PD-L1 3'-downstream; PMID:36055241, PMID:19641142, PMID:32814902, PMID:30269996, PMID:28526819, PMID:30952982)",
        },
        "X4": {
            "symbol": "D",
            "name": "TYK2_Intracellular_IFN_Signal_Dampening",
            "gate": "Gate_2B_Autoimmune_IFN_Axis_and_Stage0_Yield_Gate",
            "factor_class": "Single_Locus_Knockout (evaluated alongside Stage 6/7 post-endocrine inducible TYK2 degrader / allosteric TYK2i BMS-986165 satellite control to separate developmental KRAS/EP defect from mature beta-cell IFN-alpha protection; PMID:36289205)",
            "level_minus1": "TYK2+/+ (wild-type TYK2 preserving full Stage 0->7 endocrine differentiation yield and IFN/STAT1 promoter induction)",
            "level_plus1": "TYK2-/- (constitutive biallelic knockout dampening IFN-alpha/IL-12/IL-23 signaling while measuring Stage 0->7 monohormonal INS+ yield penalty; Chandra et al. 2022 PMID:36289205)",
        },
        "X5": {
            "symbol": "E",
            "generator_tier1": "X5 = X1 * X2 * X3",
            "name": "TXNIP_Oxidative_Stress_Shield",
            "gate": "Gate_1_Metabolic_Oxidative_Survival",
            "factor_class": "Single_Locus_Knockout",
            "level_minus1": "TXNIP+/+ (wild-type TXNIP)",
            "level_plus1": "TXNIP-/- (biallelic thioredoxin-interacting protein knockout protecting against glucotoxic/oxidative beta-cell apoptosis; Chen et al. 2008 PMID:18552236; Shalev 2014 PMID:24911120)",
        },
        "X6": {
            "symbol": "F",
            "generator_tier1": "X6 = X2 * X3 * X4",
            "name": "RNLS_ER_Stress_vs_GSIS_Shutoff_Test",
            "gate": "Gate_1_Metabolic_GSIS_Shutoff_Graveyard004",
            "factor_class": "Single_Locus_Knockout",
            "level_minus1": "RNLS+/+ (wild-type RNLS)",
            "level_plus1": "RNLS-/- (renalase knockout testing ER-stress protection vs 2.8 mM glucose shut-off delay; PMID:32719542 vs PMID:39988068)",
        },
    }
    factorial_runs = []
    tier1_tuples = []
    for run_idx in range(16):
        x1 = 1 if (run_idx & 1) else -1
        x2 = 1 if (run_idx & 2) else -1
        x3 = 1 if (run_idx & 4) else -1
        x4 = 1 if (run_idx & 8) else -1
        x5 = x1 * x2 * x3
        x6 = x2 * x3 * x4
        tier1_tuples.append((x1, x2, x3, x4, x5, x6))
        factorial_runs.append({
            "run_id": f"HIP_FACTORIAL_RUN_{run_idx + 1:02d}",
            "run_number": run_idx + 1,
            "stage_block": 0,
            "block": "Tier1_16Run_Resolution_IV_Screening",
            "coded_levels": {"X1": x1, "X2": x2, "X3": x3, "X4": x4, "X5": x5, "X6": x6},
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(x1, x2),
            "embodiment_a_star_promoted_locus1_cassette_id": "LOC1_X1P_X2P_ARM1B_ASTAR_4143BP_BIDIR_2AFREE",
            "locus3_basalfloor_promoted_cassette_id": "LOC3_X3P_L3_BFB_4078BP_BASALFLOOR_SCARFREE",
            "decoded_genotype": {
                "X1_MHC_I": "B2M::HLA-E*01:03(VMAPRTLVL)" if x1 == 1 else "B2M-/-",
                "X2_CD47": "Constitutive_HighFloor_3prime_CD47_IgV_GPI_Locked" if x2 == 1 else "Endogenous_B2M_Promoter_CD47_WT",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1" if x3 == 1 else "CXCL10+/+",
                "X4_TYK2": "TYK2-/-" if x4 == 1 else "TYK2+/+",
                "X5_TXNIP": "TXNIP-/-" if x5 == 1 else "TXNIP+/+",
                "X6_RNLS": "RNLS-/-" if x6 == 1 else "RNLS+/+",
            },
        })

    tier2_core_tuples = [
        (1, -1, 1, -1, -1, -1),   # Run 17 (and Replicate Run 27): LEGACY lean chassis (endogenous B2M-promoter bicistronic CD47, X2 = -1); NOT the current Scenario A lead (Runs 35-36)
        (-1, -1, -1, -1, -1, 1),  # Run 18 (and Replicate Run 28): Single-factor RNLS-/- on UP421 baseline
        (-1, -1, -1, -1, 1, -1),  # Run 19 (and Replicate Run 29): Single-factor TXNIP-/- on UP421 baseline
        (-1, -1, -1, -1, 1, 1),   # Run 20 (and Replicate Run 30): Dual metabolic TXNIP-/- + RNLS-/- on UP421 baseline
        (-1, -1, -1, 1, -1, -1),  # Run 21 (and Replicate Run 31): Single-factor TYK2-/- on UP421 baseline
        (-1, -1, -1, 1, 1, -1),   # Run 22 (and Replicate Run 32): TYK2-/- + TXNIP-/- de-aliasing run
        (-1, -1, 1, -1, -1, -1),  # Run 23 (and Replicate Run 33): Single-factor CXCL10::SERPINB9-GSG-tPT2A-PD-L1 on UP421 baseline
        (-1, 1, -1, -1, -1, -1),  # Run 24 (and Replicate Run 34): Single-factor Constitutive_HighFloor_Monocistronic_CD47_IgV_GPI on UP421 baseline
    ]
    tier2_bridge_tuples = [
        tier1_tuples[0],   # Run 25: Bridge Control repeating Tier-1 Run 01 (-1,-1,-1,-1,-1,-1) inside Tier-2 batch
        tier1_tuples[15],  # Run 26: Bridge Control repeating Tier-1 Run 16 (+1,+1,+1,+1,+1,+1) inside Tier-2 batch
    ]
    tier2_all_tuples = (
        [(t, "Tier2_DeAliasing_Replicate_A") for t in tier2_core_tuples]
        + [(t, "Tier2_InterStage_Bridge_Control") for t in tier2_bridge_tuples]
        + [(t, "Tier2_DeAliasing_Replicate_B") for t in tier2_core_tuples]
    )
    for offset, ((x1, x2, x3, x4, x5, x6), sub_block) in enumerate(tier2_all_tuples, 17):
        factorial_runs.append({
            "run_id": f"HIP_FACTORIAL_RUN_{offset:02d}",
            "run_number": offset,
            "stage_block": 1,
            "block": sub_block,
            "coded_levels": {"X1": x1, "X2": x2, "X3": x3, "X4": x4, "X5": x5, "X6": x6},
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(x1, x2),
            "embodiment_a_star_promoted_locus1_cassette_id": "LOC1_X1P_X2P_ARM1B_ASTAR_4143BP_BIDIR_2AFREE",
            "locus3_basalfloor_promoted_cassette_id": "LOC3_X3P_L3_BFB_4078BP_BASALFLOOR_SCARFREE",
            "decoded_genotype": {
                "X1_MHC_I": "B2M::HLA-E*01:03(VMAPRTLVL)" if x1 == 1 else "B2M-/-",
                "X2_CD47": "Constitutive_HighFloor_3prime_CD47_IgV_GPI_Locked" if x2 == 1 else "Endogenous_B2M_Promoter_CD47_WT",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1" if x3 == 1 else "CXCL10+/+",
                "X4_TYK2": "TYK2-/-" if x4 == 1 else "TYK2+/+",
                "X5_TXNIP": "TXNIP-/-" if x5 == 1 else "TXNIP+/+",
                "X6_RNLS": "RNLS-/-" if x6 == 1 else "RNLS+/+",
            },
        })

    for rep_idx, run_number in enumerate((35, 36), start=1):
        x1, x2, x3, x4, x5, x6 = HIP_LEAD_CHASSIS_A_CODED
        factorial_runs.append({
            "run_id": f"HIP_FACTORIAL_RUN_{run_number:02d}",
            "run_number": run_number,
            "stage_block": 1,
            "block": HIP_LEAD_CHASSIS_BLOCK,
            "coded_levels": {"X1": x1, "X2": x2, "X3": x3, "X4": x4, "X5": x5, "X6": x6},
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(x1, x2),
            "embodiment_a_star_promoted_locus1_cassette_id": "LOC1_X1P_X2P_ARM1B_ASTAR_4143BP_BIDIR_2AFREE",
            "locus3_basalfloor_promoted_cassette_id": "LOC3_X3P_L3_BFB_4078BP_BASALFLOOR_SCARFREE",
            "decoded_genotype": {
                "X1_MHC_I": "B2M::HLA-E*01:03(VMAPRTLVL)",
                "X2_CD47": "Constitutive_HighFloor_3prime_CD47_IgV_GPI_Locked",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1",
                "X4_TYK2": "TYK2+/+",
                "X5_TXNIP": "TXNIP+/+",
                "X6_RNLS": "RNLS+/+",
            },
        })
    dalg = _hip_protocol_algebra_fields(factorial_runs)


    whole_product_safety_model = {
        "model_status": "Conditional_Preclinical_Safety_Target_Specification (Stage 3 entry gate; NOT demonstrated clinical safety)",
        "patient_level_target_justification": (
            "Targets <10^-3 per-patient probability of any prodrug-resistant clonogenic escape across a full curative pediatric dose "
            "(lambda = N_prolif * p_escape = 8.0e4 * 1.0e-8 = 8.0e-4 under independent single-cell escape, bounded below by the non-i.i.d. common-mode floor)."
        ),
        "effective_denominator_definition": (
            "Effective detectable clonogenic opportunities (N_prolif = 8.0e4 Ki-67+/SOX9+/OCT4+ residual proliferating cells capable of colony formation "
            "across undifferentiated, ductal/progenitor, and off-target states, plus <10^-6 TRA-1-60+/OCT4+ gating on both endocrine and iEC/iPericyte vascular streams), "
            "NOT total post-mitotic endocrine cells dispensed (8.0e8) or total living cells dispensed (9.335e8 = 8.0e8 endocrine + 1.335e8 iEC/iPericytes; R3-SCI-01, R3-SCI-02)."
        ),
        "total_graft_cells_per_curative_dose": 8.0e8,
        "total_endocrine_cells_per_curative_dose": 8.0e8,
        "total_vascular_stromal_cells_per_curative_dose": 1.335e8,
        "total_living_cells_per_curative_dose": 9.335e8,
        "post_mitotic_endocrine_fraction": 0.9999,
        "lot_release_max_residual_proliferating_fraction_f_prolif": 1.0e-4,
        "max_clonogenic_proliferating_cells_per_dose_N_prolif": 8.0e4,
        "target_per_clonogenic_cell_escape_bound_p_escape": 1.0e-8,
        "lot_release_clonogenic_assay_cells_screened_for_95pct_confidence": 3.0e8,
        "exact_clonogenic_assay_cells_for_95pct_confidence_unrounded": 299573227.3553991,
        "required_physical_cells_plated_at_33p4pct_detection_efficiency": 896925831,
        "required_physical_cells_plated_for_literal_300M_effective_cells_at_33p4pct": 898203593,
        "master_cell_bank_pooled_validation_cells_plated": 3.0e8,
        "lot_release_feasible_sub_lot_aliasing_cells_plated": 1.0e7,
        "expected_surviving_clonogenic_cells_per_dose_lambda": 8.0e-4,
        "per_dose_probability_of_any_clonogenic_escape": 0.0007997,
        "whole_product_proliferating_cells_total_N_prolif": 93350,
        "whole_product_expected_surviving_clonogenic_cells_lambda": 0.0009335,
        "whole_product_per_dose_probability_of_any_clonogenic_escape": 0.0009331,
        "pan_lineage_proliferating_cell_budget_sol61_r2": {
            "endocrine_proliferating_cells_at_f_prolif_1e_4": 80000,
            "vascular_iec_pericyte_proliferating_cells_at_f_prolif_1e_4": 13350,
            "total_pan_lineage_proliferating_cells_N_prolif_total": 93350,
            "expected_surviving_clonogenic_cells_pan_lineage_lambda": 0.0009335,
            "per_dose_probability_of_any_clonogenic_escape_pan_lineage": 0.0009331,
            "maximum_common_mode_q_to_keep_pan_lineage_total_risk_below_1e_3": 0.0000670,
            "maximum_common_mode_q_to_keep_pan_lineage_total_risk_below_1e_3_unrounded": 0.0000669981,
            "note": (
                "Explicitly accounts for both the 8.0e8 endocrine cell stream (80,000 residual proliferating cells at f_prolif <= 10^-4) "
                "and the 1.335e8 vascular iEC/iPericyte stream (13,350 residual proliferating cells if f_prolif <= 10^-4, or 93,350 total "
                "pan-lineage clonogenic opportunities per 9.335e8-cell graft): at p_escape <= 1.0e-8, pan-lineage Poisson hazard is "
                "lambda_total = 93,350 * 1.0e-8 = 9.335e-4 < 1.0e-3 (P(>=1 escape) = 0.0009331), keeping P_escape_total = q + (1 - q)*r_iid <= 1.0e-3 "
                "requires q <= (10^-3 - r_iid)/(1 - r_iid) = 0.0000669981 (~6.70e-5, vs 6.696e-5 under linear union-bound subtraction), and achieving N_eff = 299,573,227.36 "
                "effective clonogenic opportunities at empirical spike-in detection efficiency eta_detection = 0.334 requires plating "
                "N_input = 299,573,227.36 / 0.334 = 896,925,831 (~8.97e8) physical clonogenic cells across cumulative MCB qualification lots "
                "(AUDIT-REF & S55R2-D1-05)."
            ),
        },
        "clonogenic_validation_tier_distinction": {
            "routine_per_lot_sub_lot_screen": {
                "cells_plated_n": 1.0e7,
                "cells_plated_n_is_effective_clonogenic_opportunities_legacy_alias": True,
                "cells_plated_n_unit_semantics_note": "Deprecated effective-opportunity alias (N_eff = 1.0e7); use required_physical_cells_plated_at_33p4pct_detection_efficiency (29,940,120) for physical cell plating.",
                "effective_clonogenic_opportunities_n": 1.0e7,
                "physical_cells_plated_n": 29940120,
                "required_physical_cells_plated_at_33p4pct_detection_efficiency": 29940120,
                "zero_event_upper_95pct_per_cell_bound_p95": 2.9957e-7,
                "standalone_per_dose_escape_risk_at_p95_for_N_prolif_8e4": 0.02368,
                "spiked_ipsc_positive_control_denominator_reconciliation": {
                    "spiked_viable_ipscs_per_well_at_0_1pct_of_1_9356e7_cells": 19356,
                    "single_cell_matrigel_clonogenic_plating_efficiency": 0.20,
                    "effective_colony_forming_ipsc_units_n_eff_per_well": 3871,
                    "zero_colony_upper_95pct_bound_per_viable_spiked_ipsc": 0.0001548,
                    "zero_colony_upper_95pct_bound_per_colony_forming_unit": 0.0007739,
                    "note": (
                        "Reconciles the two spiked-iPSC denominators (Finding 055): spiking 0.1% undifferentiated iPSCs into "
                        "1.9356x10^7 Stage 7 cells introduces 19,356 total viable iPSCs per well; at 20% single-cell Matrigel/Y-27632 "
                        "clonogenic plating efficiency, this corresponds to n_eff = 0.20 * 19,356 = 3,871 effective colony-forming iPSC units "
                        "(or 23,227 colony-forming units pooled across 6 replicate wells)."
                    ),
                },
                "role": "Routine per-lot quality-control release screen across N_eff = 1.0x10^7 effective clonogenic opportunities (29,940,120 physical clonogenic cells plated at eta_detection = 0.334) catching gross batch silencing or contamination (cannot by itself prove p_escape <= 10^-8; AUDIT-REF).",
            },
            "master_cell_bank_pooled_validation": {
                "cells_plated_n": 3.0e8,
                "cells_plated_n_is_effective_clonogenic_opportunities_legacy_alias": True,
                "cells_plated_n_unit_semantics_note": "Deprecated effective-opportunity alias (N_eff = 3.0e8); use physical_cells_plated_n (896,925,831 for exact 10^-8 bound or 898,203,593 for literal 300M) for physical cell plating.",
                "effective_clonogenic_opportunities_n": 3.0e8,
                "effective_clonogenic_opportunities_screened_n_eff": 3.0e8,
                "effective_clonogenic_opportunities_for_1e8_bound_unrounded_n": 299573227.3553991,
                "effective_clonogenic_opportunities_at_896925831_physical_cells_n": 299573227.554,
                "physical_cells_plated_n": 896925831,
                "physical_cells_plated_at_eta_0p334_n_input": 896925831,
                "physical_cells_plated_for_1e8_bound_at_eta_0p334_n": 896925831,
                "physical_cells_plated_for_literal_300M_at_eta_0p334_n": 898203593,
                "physical_cells_plated_for_1e8_bound_at_clopper_pearson_eta_l95_0p269610_n": 1111135446,
                "required_physical_cells_plated_at_clopper_pearson_lower_95pct_efficiency_26p961pct": 1111135446,
                "spike_in_detection_efficiency_point_estimate_50_of_150": 0.333333,
                "spike_in_detection_efficiency_nominal_protocol_eta": 0.334,
                "spike_in_detection_efficiency_clopper_pearson_lower_95pct_eta_l95": 0.269610,
                "zero_event_upper_95pct_per_cell_bound_at_299p57M_p95": 1.0000e-8,
                "zero_event_upper_95pct_per_cell_bound_p95": 9.9858e-9,
                "wunderlich_2022_reference_benchmark_cells_tested": 8.0e8,
                "wunderlich_2022_zero_event_upper_95pct_bound": 3.7447e-9,
                "joint_confidence_and_pan_lineage_scope_qualification_astra_r3": (
                    "Reconciles exact vs rounded denominators and qualification scope (AUDIT-REF & S2-L1-06): "
                    "(1) N_eff = -ln(0.05)/10^-8 = 299,573,227.36 effective clonogenic opportunities yields p_95 = 1.0000e-8 and requires ceil(299,573,227.36 / 0.334) = 896,925,831 physical cells at nominal eta = 0.334 "
                    "(whereas literal rounded N_eff = 300,000,000 yields p_95 = 9.9858e-9 and requires ceil(300,000,000 / 0.334) = 898,203,593 physical cells); "
                    "(2) in the 50/150 = 0.333333 empirical spike-in calibration, the one-sided 95% Clopper-Pearson lower bound is eta_L95 = Beta_inv(0.05; 50, 101) = 0.269610 (requiring ceil(299,573,227.36 / 0.269610) = 1,111,135,446 physical cells if propagating eta_L95, with Bonferroni/union-bound simultaneous confidence >= 90% across two independent 95% bounds or >= 85% across three); "
                    "(3) because StemFlex + AP+ colony counting specifically assays residual hPSCs, pan-lineage N_prolif = 93,350 qualification additionally requires parallel post-AP1903 proliferative recovery screening of the 13,350 vascular iEC/iPericyte stream in EGM-2/pericyte medium with CD31+/CD144+/PDGFRB+/EdU+ readout; "
                    "(4) N_prolif = 93,350 represents the implant-time residual Ki67+ clonogenic pool (requiring pre-implant biallelic iCasp9 verification so post-engraftment vascular/stromal turnover cannot expand an unswitched clone)."
                ),
                "lot_release_vs_mcb_qualification_distinction_sonnet55_r2": (
                    "The 8.97x10^8 physical clonogenic cell screen (N_eff = 2.9957x10^8 [~3.0x10^8] at eta_detection = 0.334; 898,203,593 physical cells for literal 3.0000x10^8) is a one-time pre-clinical Master Cell Bank (MCB) "
                    "qualification pooled across 3-5 qualification differentiation campaigns (~1 clinical-dose equivalent total), whereas routine clinical lot release "
                    "uses the N_eff = 1.0x10^7 sub-lot clonogenic assay (29,940,120 physical cells plated at eta = 0.334) paired with targeted duplex sequencing of iCasp9 and promoter CpG methylation bisulfite sequencing (S55R2-D1-09 & AUDIT-REF)."
                ),
                "role": "Cumulative pre-clinical Master Cell Bank (MCB) qualification required to demonstrate p_escape <= 10^-8 across >=2.9957x10^8 (~3.0x10^8) effective clonogenic opportunities (896,925,831 physical clonogenic cells plated at eta_detection = 0.334 for exact 10^-8 bound, or 898,203,593 physical cells for literal 300M effective opportunities; <=0.08% endocrine-only / <=0.0933% pan-lineage i.i.d. per-dose risk; Wunderlich et al. 2022 PMID:35795779).",
            },
            "one_time_mcb_qualification": {
                "cells_plated_n": 3.0e8,
                "cells_plated_n_is_effective_clonogenic_opportunities_legacy_alias": True,
                "cells_plated_n_unit_semantics_note": "Deprecated effective-opportunity alias (N_eff = 3.0e8); use physical_cells_plated_n (896,925,831 for exact 10^-8 bound or 898,203,593 for literal 300M) for physical cell plating.",
                "effective_clonogenic_opportunities_n": 3.0e8,
                "effective_clonogenic_opportunities_screened_n_eff": 3.0e8,
                "effective_clonogenic_opportunities_for_1e8_bound_unrounded_n": 299573227.3553991,
                "physical_cells_plated_n": 896925831,
                "physical_cells_plated_for_literal_300M_at_eta_0p334_n": 898203593,
                "required_physical_cells_plated_at_clopper_pearson_lower_95pct_efficiency_26p961pct": 1111135446,
                "spike_in_detection_efficiency_nominal_protocol_eta": 0.334,
                "spike_in_detection_efficiency_clopper_pearson_lower_95pct_eta_l95": 0.269610,
                "zero_event_upper_95pct_per_cell_bound_at_299p57M_p95": 1.0000e-8,
                "zero_event_upper_95pct_per_cell_bound_p95": 9.9858e-9,
            },
        },
        "quantitative_common_mode_failure_floor": {
            "exact_conditional_mixture_formula": "P_escape_total = q_common_mode + (1 - q_common_mode) * (1 - (1 - p_iid) ** N_prolif)",
            "conservative_union_bound_approximation": "P_escape_total <= min(1.0, (1 - exp(-N_prolif * p_iid)) + q_common_mode)",
            "formula": "P_escape_total = q_common_mode + (1 - q_common_mode) * r_iid (bounded above by min(1, r_iid + q_common_mode))",
            "single_locus_aa_vs1_only_common_mode_floor_estimate": 1.0e-5,
            "dual_unlinked_loci_aavs1_plus_gapdh_common_mode_floor_target": 1.0e-7,
            "maximum_common_mode_q_to_keep_total_risk_below_1e_3": 0.0002005,
            "maximum_common_mode_q_to_keep_pan_lineage_93350_risk_below_1e_3": 0.0000670,
            "exact_unrounded_common_mode_q_max_for_pan_lineage_1e_3": 0.000066998080,
            "total_pan_lineage_escape_risk_at_rounded_q_6p70e_5": 0.0010000019,
            "nonparametric_zero_failure_batches_required_for_95pct_upper_bound_6p70e_5": 44713,
            "nonparametric_zero_failure_batches_required_for_95pct_upper_bound_exact_qmax": 44713,
            "nonparametric_zero_failure_batches_required_for_95pct_upper_bound_literal_6p70e_5": 44711,
            "scenario_assumption_status": (
                "Illustrative conditional engineering scenario assumptions (NOT empirically measured clinical failure frequencies). "
                "When endocrine-only r_iid = 1 - (1 - 10^-8)^(8e4) = 7.9968e-4, keeping P_escape_total = q + (1 - q)*r_iid <= 1.0e-3 requires q <= 2.005e-4; "
                "when pan-lineage r_iid = 1 - (1 - 10^-8)^(93,350) = 9.33064516e-4 (including 13,350 vascular iEC/iPericyte proliferating cells), keeping P_escape_total <= 1.0e-3 "
                "requires exact q <= (10^-3 - r_iid)/(1 - r_iid) = 0.000066998080 (~6.70e-5, where rounded-up 6.70e-5 gives P_escape_total = 0.0010000019 and linear union-bound subtraction gives 6.6935e-5; "
                "plating N_eff = -ln(0.05)/10^-8 = 299,573,227.36 effective cells at point-estimate eta = 0.334 requires 896,925,831 physical cells [or N_eff / eta_95_lower when propagating the Clopper-Pearson lower 95% confidence bound on the 50/150 spike-in recovery], whereas literal N_eff = 300,000,000 requires 898,203,593 physical cells, and if non-proliferating Stage-7 cells do not form colonies in the clonogenic assay then N_eff denominators count MACS-enriched proliferating/Ki67+ cells plated rather than bulk post-mitotic beta cells). "
                "Crucially, q_max = 6.70e-5 is a conditional preclinical risk-budget parameter (proving q_common <= q_max via non-parametric zero-failure batch testing 1 - 0.05^(1/n_batches) <= q_max would require "
                "ceil(ln(0.05)/ln(1 - 0.000066998080)) = 44,713 independent differentiation batches at exact unrounded q_max, or ceil(ln(0.05)/ln(1 - 6.70e-5)) = 44,711 batches at literal rounded 6.70e-5), and the Stage-0 5-aza-dC / BCL2L1 spike-in assay is a mechanistic positive-control challenge rather than an empirical prevalence bound on q_common (AUDIT-REF, AUDIT-REF & AUDIT-REF)."
            ),
            "note": (
                "Even when i.i.d. per-cell escape is <=10^-8, shared batch-level promoter methylation, pre-differentiation clonal mutation, "
                "or avascular prodrug under-exposure creates a non-i.i.d. common-mode floor unless mitigated by unlinked dual-locus insertion and post-differentiation QC."
            ),
            "common_mode_silencing_and_en_bloc_surgical_retrieval_gate_sol61_r5": {
                "finding_id": "AUDIT-REF",
                "decision_id": "DEC-044 / DEC-046",
                "primary_citations": ["PMID:35795779", "PMID:26455413"],
                "single_cell_correlated_escape_formula": "P_escape_cell = p_common_cell + (1 - p_common_cell) * (p_ind ** 2)",
                "per_graft_vs_per_cell_common_mode_distinction_astra6_r5": (
                    "Layer-1 q_common_mode <= 6.6998e-5 is the per-graft/per-batch catastrophic failure probability in P_escape_total = q_common_mode + (1 - q_common_mode)*(1 - (1 - p_cell)^93350), "
                    "which at q_common_mode = 6.6998e-5 matches the 10^-3 total dose-level ceiling when per-cell escape is at the two-independent-allele floor p_cell = p_ind^2 = 10^-8 (leaving remaining cell-level headroom p_cell <= 1.00000094e-8, i.e. p_common_cell <= 9.4e-16 when q_common_mode = 6.6998e-5 is simultaneously active, whereas under a pure cell-level budget q_common_mode = 0 the 10^-3 dose-level ceiling allows p_cell <= 1.0718e-8, i.e. per-cell p_common_cell <= 7.18e-10; DEC-078); "
                    "if per-cell trans-epigenetic/XIAP resistance is p_common_cell in [10^-6, 10^-5] (Layer 2), expected surviving proliferating cells per 93,350-cell dose are E = [0.0943, 0.9344] and Poisson dose-level escape probability is P(>=1 survivor) = 1 - exp(-E) = [8.9975%, 60.7192%] (90x-607x above 10^-3; AUDIT-REF [DEC-046/DEC-078])."
                ),
                "mcb_effective_clonogenic_trial_floor_n_eff": 299573228,
                "en_bloc_surgical_pouch_retrievability_gate_0s": (
                    "Positions iCasp9 + AP1903 as a high-depth pharmacological debulking switch (>= 4-5 log reduction, tested across N_eff >= 299,573,228 effective clonogenic trials with explicit k = 0 surviving colonies in Stage-0 Gate S9) "
                    "paired with primary mechanical containment and en-bloc surgical retrievability of the dual pre-vascularized omental pouches (Gate 0S-Retrieval), "
                    "so patient safety does not rely solely on extrapolated 10^-8 two-allele independence in the presence of cell-wide trans-acting epigenetic silencing or XIAP/BCL-xL apoptotic resistance (AUDIT-REF & AUDIT-REF [DEC-044/DEC-046])."
                ),
            },
        },
        "non_iid_common_mode_failure_gates": [
            {
                "mode": "Sequential_Second_Hit_Deletion_or_PreExpansion_LOH_at_Chr19q13",
                "description": (
                    "Note on mechanism (Audit A05): Once suicide cassettes are biallelically integrated at both AAVS1 (PPP1R12C Chr19q13.42) alleles, "
                    "ordinary homologous mitotic recombination or loss of a single 19q arm does NOT simultaneously eliminate both cassettes. "
                    "Loss of both copies at a single locus instead requires either (a) LOH/arm loss acting on a subclone that originally received only a monoallelic integration, "
                    "or (b) sequential second-hit focal deletion/frameshift/silencing of the remaining allele during master cell bank expansion."
                ),
                "mitigation": (
                    "Verify 100% biallelic integration in single-cell-derived master cell bank clones by digital droplet PCR (ddPCR) + long-read sequencing across Chr19q13.42; "
                    "if sequential two-hit escape exceeds 10^-8 in expanded MCBs, escalate to the K=5 dual-safe-harbor architecture adding a second unlinked cassette at GAPDH-P2A-iCasp9 (Chr12p13.31)."
                ),
            },
            {
                "mode": "Shared_Promoter_Epigenetic_CpG_Methylation",
                "description": "Correlated de novo DNMT3A/3B methylation during 35-day Stage 0->7 differentiation can silence both copies of a viral CAG promoter at AAVS1 (Ordovas et al. 2015 PMID:26455413).",
                "mitigation": "Replace viral CAG with a 455 bp HNRNPA2B1-proximal candidate (function pending; published 0.7 kb CBX3-A2UCOE is a different element) + native human EEF1A1-core promoter (230 bp, 19 CpGs, 60.0% GC) flanked by cHS4 insulators (Muller-Kuller et al. 2015 PMID:25605798; Pfaff et al. 2013 PMID:23307570; Zhang et al. 2017 PMID:28860464) and require 60-day post-differentiation + inflammatory cytokine (IFN-g/TNF-a/IL-1b) challenge before AP1903 clonogenic escape testing (adjusting the zero-event upper 95% bound by empirical colony-forming plating efficiency).",
            },
            {
                "mode": "Hypoxia_and_Cytokine_Induced_XIAP_Blockade_of_Processed_Caspase9",
                "description": (
                    "Because AP1903 chemically dimerizes iCasp9 in the cytosol, bypassing an upstream requirement for mitochondrial outer membrane permeabilization (MOMP) to initiate procaspase-9 activation (though downstream caspase-3 feedback can still recruit mitochondrial amplification and Smac/DIABLO release; Yuan et al. 2022 PMID:35013114, Guerrero et al. 2012 PMID:23539542), "
                    "upregulation of X-linked inhibitor of apoptosis protein (XIAP / BIRC4) in hypoxic or cytokine-stressed grafts can bind the exposed ATPF tetrapeptide "
                    "of cleaved Caspase-9 via its BIR3 domain and arrest executioner Caspase-3/7 activation (Srinivasula et al. 2001 Nature PMID:11242052; Red-Team RT3-03)."
                ),
                "mitigation": (
                    "Test Decision Experiment #4 clonogenic escape under combined 1% O2 hypoxia + IFN-g/TNF-a/IL-1b (XIAP-high state) and evaluate the unvalidated preclinical salvage hypothesis of co-dosing "
                    "AP1903 (2 x 0.4 mg/kg) with a SMAC mimetic / IAP antagonist (IV birinapant [PMC4742409] or oral xevinapant [noting Phase III TrilynX termination, doi:10.1200/JCO-25-00272]) to displace XIAP from iCasp9."
                ),
            },
            {
                "mode": "Avascular_Core_Prodrug_Perfusion_Deficit",
                "description": (
                    "Krogh oxygen diffusion (R <= 60 um < R_crit = 74.6 um conservative rat islet benchmark [149.2 um diameter at V_max = 0.0524 mol/m^3/s] and < R_crit = 92.6 um Stage 7 human SC-islet benchmark [185.2 um diameter at V_max = 0.034 mol/m^3/s]; Finding 056) prevents central anoxic necrosis but does not by itself prove uniform lipophilic AP1903 dimerizer "
                    "pharmacokinetic penetration into fibrotic or poorly perfused in vivo graft pockets."
                ),
                "mitigation": "Empirically validate 2-dose AP1903 (2 x 0.4 mg/kg, 48h apart) mass-spectrometry tissue penetration and >99.99% bioluminescent/histological graft clearance in vascularized in vivo sentinel grafts.",
            },
            {
                "mode": "Multiplex_DSB_Induced_TP53_and_Chr20q11_21_BCL2L1_Clonal_Selection",
                "description": (
                    "Cas9 DSBs trigger p53-dependent arrest/apoptosis in 80-95% of wild-type hPSCs (Ihry et al. 2018 PMID:29892062; Haapaniemi et al. 2018 PMID:29892067), "
                    "strongly selecting for pre-existing dominant-negative TP53 mutations (Merkle et al. 2017 PMID:28445466) and sub-karyotypic Chr20q11.21 (BCL2L1/BCL-xL) copy-number gains (Avery et al. 2013 PMID:24286026; Red-Team RT3-01)."
                ),
                "mitigation": (
                    "Execute K=4 editing as a 2-step sequential workflow (<=2 loci/step using orthogonal HiFi-SpCas9 + AsCas12a Ultra with transient i53/RAD51 mRNA HDR enhancement; Schiroli et al. 2019 PMID:30905619; "
                    "in each step, N_effective = 2 chromosomal loci are cleaved simultaneously [1 heterotypic locus pair], corresponding to N_total = 4 diploid alleles [4 inter-chromosomal + 2 homologous allele pairs]; Finding 028) "
                    "and enforce 3-timepoint longitudinal QC (Seed -> MCB -> Stage 7 EOPC): error-corrected duplex sequencing of TP53/PPM1D/BCOR (noting that 2,000x duplex depth yields only 1 - exp(-0.2) = 18.1% single-molecule sensitivity at VAF = 10^-4, whereas >=95% single-molecule sensitivity at VAF = 10^-4 requires >=29,956 independent duplex molecules [>=30,000x target-enriched duplex depth, and >=62,958 duplex molecules (~63,000x depth) for >=95% power to observe >=3 mutant duplex molecules, vs 87.5% power at 50,000x; Finding 082]), "
                    "Chr20q11.21 (BCL2L1) ddPCR (CN = 2.00 +/- 0.15 [1.85, 2.15]), functional Nutlin-3a p53 transcriptional response (10 uM Nutlin-3a 8-h p21/CDKN1A and MDM2 mRNA fold induction >= 5.0x; obs_nutlin3a_p21_mdm2_fold_induction >= 5.0), 15-25 kb long-read ONT/PacBio + +/-50 kb flanking SNP LOH phasing (Kosicki et al. 2018 PMID:30010673), and CAST-Seq translocation assay (Turchiano et al. 2021 PMID:33626327)."
                ),
            },
        ],
    }

    col_sums = {k: sum(int(r["coded_levels"][k]) for r in factorial_runs) for k in HIP_FACTOR_KEYS}
    plus_one_counts = {k: sum(1 for r in factorial_runs if int(r["coded_levels"][k]) == 1) for k in HIP_FACTOR_KEYS}

    added_bridge_runs_37_to_40 = [
        {
            "run_id": "HIP_FACTORIAL_RUN_37",
            "run_number": 37,
            "tier": "Tier4_Stage0_Stage1_Bridge_Extension",
            "block": "Tier4_Stage0_Stage1_Bridge_Extension",
            "stage_block": 0,
            "stage_block_B": 0,
            "X1": 1,
            "X2": 1,
            "X3": 1,
            "X4": -1,
            "X5": -1,
            "X6": -1,
            "continuous_Z3_tier": 1.0,
            "continuous_Z3_tier_mode_b_quadratic_q3_option": 0.0,
            "bridged_with_runs": ["HIP_FACTORIAL_RUN_35", "HIP_FACTORIAL_RUN_36"],
            "coded_levels": {
                "X1": 1,
                "X2": 1,
                "X3": 1,
                "X4": -1,
                "X5": -1,
                "X6": -1,
            },
            "named_coded_levels": {
                "X1_HLAE_SCT": 1,
                "X2_CD47": 1,
                "X3_CXCL10_KO_SENTINEL": 1,
                "X4_TYK2_KO": -1,
                "X5_TXNIP_KO": -1,
                "X6_RNLS_KO": -1,
            },
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(1, 1),
            "decoded_genotype": {
                "X1_MHC_I": "B2M::HLA-E*01:03(VMAPRTLVL)",
                "X2_CD47": "Constitutive_HighFloor_3prime_CD47_IgV_GPI_Locked",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1",
                "X4_TYK2": "TYK2+/+",
                "X5_TXNIP": "TXNIP+/+",
                "X6_RNLS": "RNLS+/+",
            },
            "role": "Stage-0 Lead Genotype bridge anchor pairing with Stage-1 Runs 35-36 (bridge pair 3) & Mode-B quadratic Q3=Z3^2-mean(Z3^2) identification anchor",
        },
        {
            "run_id": "HIP_FACTORIAL_RUN_38",
            "run_number": 38,
            "tier": "Tier4_Stage0_Stage1_Bridge_Extension",
            "block": "Tier4_Stage0_Stage1_Bridge_Extension",
            "stage_block": 1,
            "stage_block_B": 1,
            "X1": 1,
            "X2": -1,
            "X3": 1,
            "X4": -1,
            "X5": -1,
            "X6": 1,
            "continuous_Z3_tier": 1.0,
            "continuous_Z3_tier_mode_b_quadratic_q3_option": 0.0,
            "bridged_with_runs": ["HIP_FACTORIAL_RUN_06"],
            "coded_levels": {
                "X1": 1,
                "X2": -1,
                "X3": 1,
                "X4": -1,
                "X5": -1,
                "X6": 1,
            },
            "named_coded_levels": {
                "X1_HLAE_SCT": 1,
                "X2_CD47": -1,
                "X3_CXCL10_KO_SENTINEL": 1,
                "X4_TYK2_KO": -1,
                "X5_TXNIP_KO": -1,
                "X6_RNLS_KO": 1,
            },
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(1, -1),
            "decoded_genotype": {
                "X1_MHC_I": "B2M::HLA-E*01:03(VMAPRTLVL)",
                "X2_CD47": "Endogenous_B2M_Promoter_CD47_WT",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1",
                "X4_TYK2": "TYK2+/+",
                "X5_TXNIP": "TXNIP+/+",
                "X6_RNLS": "RNLS-/-",
            },
            "role": "Stage-1 bridge repeat of Tier-1 Run 06 (bridge pair 4) & Stage x Main (B*X6, B*X2) de-aliasing anchor",
        },
        {
            "run_id": "HIP_FACTORIAL_RUN_39",
            "run_number": 39,
            "tier": "Tier4_Stage0_Stage1_Bridge_Extension",
            "block": "Tier4_Stage0_Stage1_Bridge_Extension",
            "stage_block": 1,
            "stage_block_B": 1,
            "X1": -1,
            "X2": 1,
            "X3": 1,
            "X4": -1,
            "X5": -1,
            "X6": -1,
            "continuous_Z3_tier": 1.0,
            "continuous_Z3_tier_mode_b_quadratic_q3_option": 0.0,
            "bridged_with_runs": ["HIP_FACTORIAL_RUN_07"],
            "coded_levels": {
                "X1": -1,
                "X2": 1,
                "X3": 1,
                "X4": -1,
                "X5": -1,
                "X6": -1,
            },
            "named_coded_levels": {
                "X1_HLAE_SCT": -1,
                "X2_CD47": 1,
                "X3_CXCL10_KO_SENTINEL": 1,
                "X4_TYK2_KO": -1,
                "X5_TXNIP_KO": -1,
                "X6_RNLS_KO": -1,
            },
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(-1, 1),
            "decoded_genotype": {
                "X1_MHC_I": "B2M-/-",
                "X2_CD47": "Constitutive_HighFloor_3prime_CD47_IgV_GPI_Locked",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1",
                "X4_TYK2": "TYK2+/+",
                "X5_TXNIP": "TXNIP+/+",
                "X6_RNLS": "RNLS+/+",
            },
            "role": "Stage-1 bridge repeat of Tier-1 Run 07 (bridge pair 5) & Stage x Main (B*X1, B*X2) de-aliasing anchor",
        },
        {
            "run_id": "HIP_FACTORIAL_RUN_40",
            "run_number": 40,
            "tier": "Tier4_Stage0_Stage1_Bridge_Extension",
            "block": "Tier4_Stage0_Stage1_Bridge_Extension",
            "stage_block": 1,
            "stage_block_B": 1,
            "X1": 1,
            "X2": -1,
            "X3": 1,
            "X4": 1,
            "X5": -1,
            "X6": -1,
            "continuous_Z3_tier": 1.0,
            "continuous_Z3_tier_mode_b_quadratic_q3_option": 0.0,
            "bridged_with_runs": ["HIP_FACTORIAL_RUN_14"],
            "coded_levels": {
                "X1": 1,
                "X2": -1,
                "X3": 1,
                "X4": 1,
                "X5": -1,
                "X6": -1,
            },
            "named_coded_levels": {
                "X1_HLAE_SCT": 1,
                "X2_CD47": -1,
                "X3_CXCL10_KO_SENTINEL": 1,
                "X4_TYK2_KO": 1,
                "X5_TXNIP_KO": -1,
                "X6_RNLS_KO": -1,
            },
            "locus1_physical_cassette_id": _locus1_physical_cassette_id(1, -1),
            "decoded_genotype": {
                "X1_MHC_I": "B2M::HLA-E*01:03(VMAPRTLVL)",
                "X2_CD47": "Endogenous_B2M_Promoter_CD47_WT",
                "X3_Locus3": "CXCL10::SERPINB9-GSG-tPT2A-PD-L1",
                "X4_TYK2": "TYK2-/-",
                "X5_TXNIP": "TXNIP+/+",
                "X6_RNLS": "RNLS+/+",
            },
            "role": "Stage-1 bridge repeat of Tier-1 Run 14 (bridge pair 6) & Stage x Main (B*X4, B*X5) de-aliasing anchor",
        },
    ]
    all_40_bridge_runs = factorial_runs + added_bridge_runs_37_to_40
    s_40_bridge = _hip_design_summary(all_40_bridge_runs)
    x23_40 = _hip_model_rows(all_40_bridge_runs)
    x29_40 = [row + [row[22] * row[k] for k in range(1, 7)] for row in x23_40]
    hat29_40 = _hip_hat_diag(x29_40)
    unit_leverage_runs_29col_40 = [r["run_id"] for r, h in zip(all_40_bridge_runs, hat29_40) if h == 1]
    if (
        s_40_bridge["rank22"] != 22
        or s_40_bridge["rank23"] != 23
        or s_40_bridge["rank29_with_stage_by_main"] != 29
        or s_40_bridge["diag23"][22] != Fraction(2893, 9566)
        or unit_leverage_runs_29col_40 != [
            "HIP_FACTORIAL_RUN_25",
            "HIP_FACTORIAL_RUN_26",
            "HIP_FACTORIAL_RUN_37",
            "HIP_FACTORIAL_RUN_38",
            "HIP_FACTORIAL_RUN_39",
            "HIP_FACTORIAL_RUN_40",
        ]
    ):
        raise ValueError(f"40-run bridge algebraic verification failed: {s_40_bridge['diag23'][22]}, unit={unit_leverage_runs_29col_40}")

    hip_core_protocol = {
        "protocol_id": "OPENISLET_HIP_CORE_2_6_MINUS_2_PLUS_BRIDGED_REPLICATED_v0.7.0",
        "execution_readiness_status": "Pre-Registered_Stage1_Protocol_Specification (Locked 4-cell X1xX2 physical construct matrix + exact 108-independent vs 75-shared-clone GLS/OLS power plan + external UP421 control arm)",
        **dalg["protocol_fields"],
        "column_sums_across_36_runs_by_factor": col_sums,
        "plus_one_counts_across_36_runs_by_factor": plus_one_counts,
        "residual_degrees_of_freedom_single_rep_per_run": len(factorial_runs) - 23,
        "locus1_x1_by_x2_physical_construct_matrix": LOCUS1_X1_BY_X2_PHYSICAL_CONSTRUCT_MATRIX,
        "execution_row_contract": {"n_run_rows": len(factorial_runs), "n_unique_genotypes": dalg["n_distinct_genotypes"],
                                   "execute_repeat_rows": True, "run_ids": [r["run_id"] for r in factorial_runs]},
        "unevaluated_61_clone_alternative": {"status": "UNEVALUATED_NOT_AUTHORIZED_BY_75_OR_108_CLONE_POWER", "n_clones": 61,
                                            "covariance_matrix": None, "standard_errors": None, "minimum_detectable_effects": None},
        "hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1": {
            "design_id": "OPENISLET_HIP_CORE_40RUN_FULL_RANK_29_OF_29_BRIDGE_EXTENSION",
            "total_runs": len(all_40_bridge_runs),
            "distinct_genotypes": 25,
            "distinct_genotypes_mode_b_continuous_midpoint_z3_0": 29,
            "matched_stage0_stage1_bridge_genotype_pairs": 6,
            "matched_stage0_stage1_bridge_genotype_pairs_mode_b_continuous_z3_0": 2,
            "total_clones_mode1_run_nested": 120,
            "total_clones_mode2_shared_genotype": 75,
            "total_clones_mode2_shared_genotype_mode_b_continuous_z3_0": 87,
            "mode_b_genotype_and_clone_sharing_note_astra_s1": (
                "In binary Mode A (Z3 = +1 on Runs 37-40), Runs 37-40 replicate Tier-1 Runs 06, 07, and 14 across stages, yielding 25 distinct genotypes, "
                "75 shared-genotype founder clones (at 3 clones/genotype), and 6 cross-stage bridge pairs (Runs 25-26 + 37-40). In monoallelic midpoint Mode B "
                "(where Runs 37-40 carry physical monoallelic Z3 = 0 knock-ins rather than Z3 = +1), Runs 37-40 constitute 4 distinct midpoint genotypes (3 of which "
                "are distinct from each other, plus Run 38/40 sharing a midpoint genotype if coded identically, or 4 new midpoint genotypes distinct from the 25 binary genotypes = 28-29 total genotypes and 84-87 founder clones), "
                "leaving Runs 25-26 as the 2 binary cross-stage bridges unless the 44-run extension (Runs 01-40 binary + Runs 41-44 midpoint, 29 genotypes / 87 shared clones) is executed (AUDIT-REF)."
            ),
            "native_hip_model_rows_contract_verified": True,
            "added_stage1_bridge_runs_37_to_40": added_bridge_runs_37_to_40,
            "operating_modes_and_condition_number_reconciliation_argon_r2": {
                "mode_a_primary_binary_bridge_Z3_plus1": {
                    "description": "Runs 37-40 executed at binary Z3 = X3 = +1.0 (full bicistronic CXCL10::SERPINB9-P2A-PD-L1 sentinel knock-in, matching Runs 35-36 and Tier-1 Runs 06, 07, 14).",
                    "rank_22col_main_plus_15_2fi": f"{s_40_bridge['rank22']}/22",
                    "condition_number_22col_raw": 11.2549,
                    "condition_number_22col_unscaled_with_intercept": 11.2549,
                    "legacy_centered_subblock_proxy_reference": 3.60,
                    "rank_23col_with_stage_block_B": f"{s_40_bridge['rank23']}/23",
                    "condition_number_23col_raw_B01": 15.21,
                    "condition_number_23col_B01_unit_column_norm_exact": 15.7407,
                    "condition_number_23col_B_pm1_raw": 14.7461,
                    "condition_number_23col_B_pm1_unit_column_norm": 13.4892,
                    "condition_number_23col_B_pm1_unit_column_norm_exact": 14.7461,
                    "condition_number_23col_B_pm1_unit_norm_equality_note": (
                        "When B in {-1, +1}, all 23 columns in {-1, +1}^40 have identical Euclidean norm sqrt(40), so unit-column-norm scaling "
                        "is an exact scalar multiple of the raw matrix and the exact Euclidean 2-norm condition number is identically 14.7461 "
                        "(whereas 15.7407 is the unit-column-norm condition number when B in {0, 1} has norm sqrt(24) and 13.4892 is a legacy centered-block proxy; S55R2-S2-02)."
                    ),
                    "condition_number_23col_scaled_B_pm1": 6.91,
                    "var_beta_B_over_sigma2_23col_exact_fraction": _hip_frac_plain(s_40_bridge["diag23"][22]),
                    "var_beta_B_over_sigma2_23col_decimal": round(float(s_40_bridge["diag23"][22]), 6),
                    "max_hat_leverage_23col_exact_fraction": _hip_frac_plain(max(s_40_bridge["hat23"])),
                    "max_hat_leverage_23col": round(float(max(s_40_bridge["hat23"])), 6),
                    "unit_leverage_runs_23col": [],
                    "rank_29col_with_6_B_times_Xi": f"{s_40_bridge['rank29_with_stage_by_main']}/29",
                    "condition_number_29col_raw_B01": 32.68,
                },
                "mode_b_quadratic_sentinel_dose_curvature_Z3_0": {
                    "description": (
                        "Optional dose-response satellite mode setting Z3 = 0.0 (defined as a single pre-registered intermediate sentinel dose: "
                        "mono-allelic CXCL10::SERPINB9-P2A-PD-L1/+ knock-in on Runs 37-40, avoiding confounding mono-allelic copy number with single-cistron "
                        "composition; AUDIT-REF) so Z3^2 in {0, 1} is non-constant across the 40 runs and Q3 = Z3^2 - 0.9 is linearly independent of "
                        "the intercept (whereas at binary Z3 = +/-1 on all 40 runs, Z3^2 = 1 is constant and Q3 = 0 has rank 23/24; AUDIT-REF)."
                    ),
                    "rank_24col_with_Q3": "24/24",
                    "condition_number_24col_raw_unscaled": 16.4523,
                    "condition_number_24col_euclidean_unit_column_norm_scaled": 16.7269,
                    "condition_number_24col_B_pm1_raw": 15.9055,
                    "condition_number_24col_B_pm1_unit_column_norm": 15.8960,
                    "condition_number_24col_column_scaled": 7.42,
                    "condition_number_15col_2fi_subblock_exact": 7.4225,
                    "max_hat_leverage_24col_mode_b_indicator_exact_fraction": "16180/19049",
                    "max_hat_leverage_24col_mode_b_indicator_decimal": 0.849388,
                    "var_beta_B_over_sigma2_24col_mode_b_indicator_exact_fraction": "6393/19049",
                    "var_beta_B_over_sigma2_24col_mode_b_indicator_decimal": 0.335608,
                    "continuous_dose_polynomial_Z3_matrix_astra_r2": {
                        "description": (
                            "Substituting continuous Z3 in {-1, 0, +1} directly into the main-effect column Z3 and all five Z3*X_k two-factor interaction columns "
                            "on midpoint Runs 37-40 (plus centered Q3 = Z3^2 - 0.9) yields an exact continuous-dose polynomial response-surface matrix where purely linear "
                            "y = Z3 and bilinear y = Z3*X5 fit beta_Q3 = 0 (whereas retaining binary X3 = +1 on Runs 37-40 forms a 3-level categorical indicator matrix "
                            "where pure curvature without interactions is parameterized by beta_Q3 - beta_X3; AUDIT-REF)."
                        ),
                        "rank_24col_continuous_Z3": "24/24",
                        "condition_number_24col_continuous_Z3_raw_B01": 16.2370,
                        "condition_number_24col_continuous_Z3_unit_column_norm_B01": 13.6894,
                        "condition_number_24col_continuous_Z3_raw_B_pm1": 15.6064,
                        "condition_number_24col_continuous_Z3_unit_column_norm_B_pm1": 12.7725,
                        "var_beta_B_over_sigma2_24col_continuous_Z3_exact_fraction": "1233461/3661912",
                        "var_beta_B_over_sigma2_24col_continuous_Z3": 0.336835,
                        "max_hat_leverage_24col_continuous_Z3_exact_fraction": "90371537/106195448",
                        "max_hat_leverage_24col_continuous_Z3_decimal": 0.850993,
                        "optional_44run_satellite_preserving_29col_bridge_rank": {
                            "total_runs": 44,
                            "structure": "40 binary Mode-A bridge runs (Runs 01-40 at Z3 = +/-1) + 4 mono-allelic Z3 = 0 midpoint satellite runs (Runs 41-44)",
                            "distinct_genotypes": 29,
                            "rank_24col_continuous_Z3": "24/24",
                            "residual_df_24col": 20,
                            "condition_number_24col_raw_B01": 14.7745,
                            "rank_30col_with_6_B_times_Xi_and_Q3": "30/30",
                            "residual_df_30col": 14,
                            "condition_number_30col_raw_B01": 27.0250,
                            "max_hat_leverage_30col": 0.95864,
                            "unit_leverage_runs_30col": [],
                        },
                    },
                    "condition_number_scaling_note_sol61_r2": (
                        "16.4523 (raw) and 16.7269 (unit-column-norm) are the exact 2-norm condition numbers of the 40x24 indicator matrix retaining X3=+1 on Runs 37-40 (15.9055 raw / 15.8960 unit-column-norm with B in {-1,+1}); "
                        "16.2370 (raw) and 13.6894 (unit-column-norm) are the exact condition numbers of the 40x24 continuous-dose polynomial matrix substituting Z3=0 into Z3 and all five Z3*X_k columns on Runs 37-40 (15.6064 raw / 12.7725 unit-column-norm with B in {-1,+1}); "
                        "and 7.42 (exact 7.4225) is the condition number of the 15-column two-factor interaction subblock X_2FI retained for backward-compatible schema checks "
                        "(AUDIT-REF, AUDIT-REF, AUDIT-REF, & S55R2-D1-03)."
                    ),
                },
            },
            "twenty_nine_column_interaction_power_diagnostics_argon_r2": {
                "residual_df_29col": 11,
                "pure_error_df": 9,
                "lack_of_fit_df_29col": 2,
                "clone_level_df_120_clones": 52,
                "clone_level_df_120_clones_mode1_iid_120_minus_29": 91,
                "clone_level_df_75_shared_genotype_clones_23col_reference": 52,
                "max_hat_leverage_29col": 1.0,
                "unit_leverage_runs_29col": unit_leverage_runs_29col_40,
                "unit_leverage_diagnostic_and_optional_46run_extension_sol61_r2": {
                    "unit_leverage_explanation": (
                        "In the 29-column exploratory model (23col + 6 B*X_i) on 40 runs, the 6 single-replicate inter-stage bridge rows "
                        "(Runs 25, 26, 37, 38, 39, 40) have hat leverage h_ii = 1.0000 (whereas in the primary 23-column model on 40 runs max h_ii = 8111/9566 = 0.847899 "
                        "[and 16180/19049 = 0.849388 in the 24-column Mode-B indicator matrix or 90371537/106195448 = 0.850993 in the 24-column continuous-Z3 polynomial matrix, vs 175/202 = 0.866337 in the 36-run 23-column baseline; S55R2-D1-04] "
                        "with 0 unit-leverage rows). Because ordinary residuals e_i = (1 - h_ii)*y_i - ... vanish identically on h_ii = 1 rows, "
                        "leave-one-out studentized deletion diagnostics cannot flag an aberrant observation on those 6 rows under the 29-column fit."
                    ),
                    "optional_46run_duplicated_bridge_extension": {
                        "total_runs": 46,
                        "added_duplicate_bridge_runs_41_to_46": "Duplicate each of the 6 bridge rows (Runs 25, 26, 37, 38, 39, 40) once within its respective stage block (0 new founder genotypes in Mode 2; 75 shared-genotype clones total)",
                        "max_hat_leverage_29col_in_46run": 0.875,
                        "max_hat_leverage_29col_exact_fraction": "7/8",
                        "unit_leverage_runs_29col_in_46run": [],
                        "residual_df_29col_in_46run": 17,
                        "pure_error_df_in_46run": 15,
                        "lack_of_fit_df_29col_in_46run": 2,
                        "stage_block_B_var_over_sigma2_29col_in_46run_exact": "3/4",
                        "stage_block_B_var_over_sigma2_29col_in_46run_decimal": 0.75,
                        "legacy_uncentered_stage_block_proxy_fraction": "9/16",
                    },
                },
                "student_t_multiplier_80pct_power_two_sided_alpha05_df11": 3.0765,
                "student_t_multiplier_80pct_power_two_sided_alpha05_df52": 2.8553,
                "per_term_variance_multipliers_se_and_mde_at_sigma_run_5p4314": {
                    "stage_block_B_29col": {"exact_var_over_sigma2": "1/1", "var_over_sigma2": 1.000000, "se_pp": 5.4314, "mde_80pct_df11_pp": 16.71, "mde_80pct_df52_pp": 15.51},
                    "B_times_X2_CD47_best": {"exact_var_over_sigma2": "5/8", "var_over_sigma2": 0.625000, "half_effect_se_pp": 4.2939, "half_effect_mde_80pct_df11_pp": 13.21, "full_interstage_shift_2beta_se_pp": 8.5878, "full_shift_mde_80pct_df11_pp": 26.42},
                    "B_times_X1_HLAE": {"exact_var_over_sigma2": "27/32", "var_over_sigma2": 0.843750, "half_effect_se_pp": 4.9891, "half_effect_mde_80pct_df11_pp": 15.35, "full_interstage_shift_2beta_se_pp": 9.9781, "full_shift_mde_80pct_df11_pp": 30.70},
                    "B_times_X4_TYK2": {"exact_var_over_sigma2": "29/32", "var_over_sigma2": 0.906250, "half_effect_se_pp": 5.1705, "half_effect_mde_80pct_df11_pp": 15.91, "full_interstage_shift_2beta_se_pp": 10.3411, "full_shift_mde_80pct_df11_pp": 31.81},
                    "B_times_X6_RNLS": {"exact_var_over_sigma2": "29/32", "var_over_sigma2": 0.906250, "half_effect_se_pp": 5.1705, "half_effect_mde_80pct_df11_pp": 15.91, "full_interstage_shift_2beta_se_pp": 10.3411, "full_shift_mde_80pct_df11_pp": 31.81},
                    "B_times_X3_CXCL10_SENTINEL": {"exact_var_over_sigma2": "5/4", "var_over_sigma2": 1.250000, "half_effect_se_pp": 6.0725, "half_effect_mde_80pct_df11_pp": 18.68, "full_interstage_shift_2beta_se_pp": 12.1450, "full_shift_mde_80pct_df11_pp": 37.36},
                    "B_times_X5_TXNIP_worst": {"exact_var_over_sigma2": "47/32", "var_over_sigma2": 1.468750, "half_effect_se_pp": 6.5824, "half_effect_mde_80pct_df11_pp": 20.25, "full_interstage_shift_2beta_se_pp": 13.1648, "full_shift_mde_80pct_df11_pp": 40.50},
                },
                "estimator_hierarchy_rule": (
                    "The 23-column additive-block model (df = 17, Var(beta_B)/sigma^2 = 2893/9566 = 0.302425, SE = 2.9870 pp, exact df=17 MDE_80% = 8.88 pp [8.95 pp at conservative df=15]) "
                    "is the primary confirmatory inter-stage bridge estimator; the 29-column 6 B*X_i model (df = 11, half-effect MDE_80% = 13.21-20.25 pp) "
                    "is an exploratory / AICc-stepwise heterogeneity screen for large inter-stage sign flips (AUDIT-REF)."
                ),
                "legacy_scalar_alias_reconciliation_note_astra_r3": (
                    "Backward-compatibility note on scalar keys (AUDIT-REF): clone_level_df_120_clones = 52 retains the Mode-2 75-clone 23-column reference (75 - 23 = 52; see clone_level_df_120_clones_mode1_iid_120_minus_29 = 91 for Mode-1 120-clone 29-column df); "
                    "condition_number_24col_with_q3 = 7.42 retains the 15-column binary 2FI subblock condition number (7.4225; see condition_number_24col_mode_b_indicator_raw = 16.4523 and condition_number_24col_mode_b_indicator_unit_column_norm = 16.7269 for the full 40x24 Mode-B matrix); "
                    "and stage_block_mde_80pct_power_two_sided_alpha05_df17_pp = 8.95 retains the conservative df=15 Student-t multiplier (2.997 * 2.9870 = 8.95 pp; see stage_block_mde_80pct_power_two_sided_alpha05_df17_exact_student_t_pp = 8.88 for exact df=17)."
                ),
            },
            "full_24_column_2fi_plus_quadratic_q3_rank": "24/24",
            "condition_number_24col_with_q3": 7.42,
            "condition_number_15col_binary_2fi_subblock": 7.4225,
            "condition_number_24col_mode_b_indicator_raw": 16.4523,
            "condition_number_24col_mode_b_indicator_unit_column_norm": 16.7269,
            "full_29_column_model_rank": "29/29",
            "condition_number_29col": 32.68,
            "var_beta_B_over_sigma2_23col_in_36run_baseline": 0.688119,
            "var_beta_B_over_sigma2_23col_in_40run_bridge": 0.302425,
            "var_beta_B_over_sigma2_24col_mode_b_indicator_in_40run_bridge": 0.335608,
            "var_beta_B_over_sigma2_24col_continuous_Z3_in_40run_bridge": 0.336835,
            "stage_block_se_pp_23col_at_sigma_run_5p4314": 2.9870,
            "stage_block_mde_80pct_power_two_sided_alpha05_df17_exact_student_t_pp": 8.88,
            "stage_block_mde_80pct_power_two_sided_alpha05_df17_pp": 8.95,
            "stage_block_mde_80pct_power_two_sided_alpha05_df15_conservative_pp": 8.95,
            "stage_block_36run_baseline_mde_80pct_df13_exact_student_t_pp": 13.66,
            "stage_block_36run_baseline_mde_80pct_conservative_pp": 14.01,
            "note": (
                "Adding 4 explicit Stage-0/Stage-1 bridge runs (Runs 37-40 with full 6-factor coded_levels X1..X6: Run 37 = Lead Genotype "
                "[+1,+1,+1,-1,-1,-1] in Stage 0 bridging Runs 35-36, and Runs 38-40 in Stage 1 repeating Tier-1 Runs 06, 07, and 14; "
                "zero new founder genotypes required in Mode 2 [75 clones total]) expands matched Stage-0/Stage-1 bridge pairs from 2 to 6, "
                "cuts the 23-column stage-block variance multiplier Var(beta_B)/sigma^2 by 56.1% from 55/80 = 0.6875 [139/202 = 0.688119 with Q3] "
                "(SE = 4.5055 pp, exact df=13 MDE_80% = 13.66 pp [14.01 pp conservative]) to 2893/9566 = 0.302425 in 23-col Mode A (6393/19049 = 0.335608 "
                "in 24-col Mode B indicator / 1233461/3661912 = 0.336835 in 24-col continuous Z3; SE = 2.9870 pp, exact df=17 "
                "Student-t MDE_80% = 8.88 pp [8.95 pp conservative df=15 multiplier]), achieves rank 24/24 (cond = 7.4225 on the 15-col binary 2FI subblock / "
                "16.4523 raw / 16.7269 unit-column-norm on the full 40x24 Mode B Z3=0.0 indicator matrix) for the "
                "23-column + quadratic sentinel curvature Q3 model, and achieves full column rank 29/29 (cond = 32.68, residual df = 11) for the 29-column "
                "model (intercept + 6 main effects + 15 2FIs + stage_block_B + 6 stage-by-main B*X_i interactions; closing AUDIT-REF, "
                "AUDIT-REF, AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
            ),
        },
        "biological_replication_and_blocking_plan": {
            "biological_replicate_unit": (
                "Mode 1 (108 Independent Run-Nested Founder Clones): 3 independently derived/edited hPSC founder clones per run across all 36 runs (108 independent founder clone units total, giving diagonal Cov(Y_bar) = sigma_run^2 * I_36). "
                "Mode 2 (75 Shared-Genotype Clones Across 25 Distinct Genotypes): n_clones = 3 independently edited, longitudinally QC-verified hPSC clones per distinct genotype "
                "(25 distinct genotypes x 3 clones = 75 edited clones), where the 11 bridged/replicated runs re-differentiate the same 3 genotype clones in a second differentiation batch "
                "(analyzed by exact GLS/REML under Cov(Y_bar) = (sigma_clone^2 / n_clones) * Z Z^T + (sigma_eps^2 / (n_clones * n_batches)) * I_36) "
                "x n_batches = 2 independent Stage 0->7 differentiation batches x n_donors = 3 HLA/KIR/SIRPA-typed primary human immune donors per run"
            ),
            "external_reference_control_arms_run_in_every_batch": [
                {
                    "control_id": "CTRL_UP421_LENTI_CD47",
                    "genotype": "B2M-/- + CIITA-/- + lentiviral CMV/EF1a-CD47-WT (>180k molecules/cell; K=3)",
                    "role": "External clinical benchmark control (Sana UP421 / Hu et al. 2024 PMID:37156915, PMID:40757665) separating lentiviral multi-copy CD47-WT overexpression from single-locus B2M exon-1 knock-in runs (Run 01 is LOC1_X1M_X2M_B2M_CD47WT_MONOCISTRONIC: B2M-/- + CIITA-/- + monocistronic endogenous B2M-promoter CD47-WT without 2A tag).",
                },
                {
                    "control_id": "CTRL_UNEDITED_WT_ISOGENIC",
                    "genotype": "Unedited parental wild-type hPSC Stage 7 SC-islet clusters (K=0)",
                    "role": "Positive lysis / alloreactivity baseline and unedited dynamic perifusion GSIS reference.",
                },
            ],
            "randomization_and_blocking": (
                f"Option A (Preferred Single-Stage Execution, 22-column model): Execute and randomize all {len(factorial_runs)} coded run rows, including all 11 repeat/bridge/lead rows (25 distinct genotypes), within shared differentiation and immune-donor blocks. Do not collapse the 36-row information matrix to 25 rows. "
                "Option B (Sequential Two-Stage Execution, 23-column model with stage_block_B): Include Runs 25-26 (exact repeats of Tier-1 Runs 01 and 16 inside Tier-2), duplicated Tier-2 Runs 17-24 and 27-34, and the two lead-chassis confirmatory replicates Runs 35-36, "
                f"giving rank([X | B]) = 23/23 and 0 unit-leverage rows (max h_ii = {dalg['max_leverage_text']}, projection residual SSR = ||(I - H_22) B||_2^2 = {dalg['ssr_text']}) so inter-stage batch shift B is linearly independent and jointly estimable with all 21 genetic main/interaction effects."
            ),
            "run_to_batch_allocation_and_within_stage_batch_covariance_sonnet55_r2": {
                "stage0_batch_allocation": (
                    "Stage-0 Runs 01-16 (plus Run 37 in the 40-run bridge) are differentiated across n_batches = 2 independent Stage 0->7 differentiation batches "
                    "(Batch S0-A and Batch S0-B), with every Stage-0 run executed in both Batch S0-A and Batch S0-B (complete run-by-batch crossing within Stage 0) "
                    "so all 6 main effects X1..X6 and 15 2FIs are strictly orthogonal to the within-Stage-0 batch contrast (Batch S0-A vs Batch S0-B)."
                ),
                "stage1_batch_allocation": (
                    "Stage-1 Runs 17-36 (plus Runs 38-40 in the 40-run bridge) are differentiated across n_batches = 2 independent Stage 0->7 differentiation batches "
                    "(Batch S1-A executing Tier-2 Block A [Runs 17-24], bridge Run 25, and Lead Run 35; Batch S1-B executing exact Tier-2 duplicate Block B [Runs 27-34], "
                    "bridge Run 26, Lead Run 36, and Runs 38-40), with explicit within-stage batch random-intercept covariance tau_batch^2 * B_batch * B_batch^T "
                    "included in the REML/GLS mixed model (or eliminated for main effects via complete run-by-batch crossing; under literal S1-A vs S1-B split-batch incidence "
                    "B_batch * B_batch^T at tau_batch = 2.5 / 5.0 pp, Mode 2 [75 shared-genotype clones across 25 distinct genotypes, Cov = (64/3)*Z*Z^T + v_eps*I + tau^2*B_batch*B_batch^T] "
                    "gives main-effect GLS SE = [1.2655, 1.3461] / [1.2655, 1.3462] pp, 2FI SE = [1.8054, 4.0872] / [1.8055, 4.0873] pp, and SE(B) = 3.7139 / 5.7046 pp under 2-batch residual v_eps = 49/6 "
                    "(or main [1.2802, 1.3502] / [1.2802, 1.3503] pp, 2FI [1.9109, 4.5345] / [1.9109, 4.5347] pp, and SE(B) = 4.1258 / 5.9810 pp under 1-batch-per-Stage-1-row residual v_eps = 49/3), "
                    "whereas Mode 1 [108 independent run-nested clone differentiations, Cov = diag(64/3 + v_eps) + tau^2*B_batch*B_batch^T] gives main-effect GLS SE = [1.2359, 1.3086] / [1.2363, 1.3093] pp, "
                    "2FI SE = [1.4499, 4.6648] / [1.4502, 4.6655] pp, and SE(B) = 5.1540 / 6.7319 pp under 2-batch residual v_eps = 49/6 (or main [1.2526, 1.3159] / [1.2530, 1.3167] pp, "
                    "2FI [1.5876, 5.0166] / [1.5879, 5.0173] pp, and SE(B) = 5.3910 / 6.9151 pp under 1-batch-per-Stage-1-row residual v_eps = 49/3), as tabulated in differentiation_batch_variance_sensitivity_sweep; "
                    "S55R2-S2-07, AUDIT-REF & AUDIT-REF)."
                ),
            },
        },
        "statistical_analysis_and_power_plan": {
            "linear_mixed_effects_model_specification": (
                "Y_{ijkd} = x_i^T * beta + u_clone(g(i),j) + b_batch(k) + v_donor(d) + epsilon_{ijkd}, where x_i in R^23 is the coded design row "
                "[1, X1..X6, 15 2FIs, stage_block_B], u_clone(g(i),j) ~ N(0, sigma_clone^2) is the random intercept for independently edited hPSC clone j "
                "of genotype g(i) (preventing clonal pseudo-replication), b_batch(k) ~ N(0, sigma_batch^2) is the differentiation batch random effect, "
                "v_donor(d) ~ N(0, sigma_donor^2) is the PBMC/NK/macrophage donor random effect, and epsilon_{ijkd} ~ N(0, sigma_eps^2) is the residual assay error."
            ),
            "effective_run_level_variance_formula": (
                "Defining sigma_eps = 7.0 percentage points (pp) as the donor-averaged clone-by-batch cell residual (averaged across n_donors = 3 immune donors after donor-block centering), "
                "averaging over n_clones = 3 clones and n_batches = 2 differentiation batches per run yields effective run-mean variance "
                "sigma_run^2 = sigma_clone^2 / n_clones + sigma_eps^2 / (n_clones * n_batches) = 8.0^2 / 3 + 7.0^2 / 6 = 177/6 = 29.50 pp^2 -> sigma_run = 5.43 pp "
                "(whereas if 7.0 pp represents the single-donor well residual epsilon_{ijkd}, averaging across n_donors = 3 donors gives "
                "sigma_run,3donors = sqrt(8.0^2 / 3 + 7.0^2 / (3 * 2 * 3)) = sqrt(433/18) = 4.90 pp; and single-clone pseudo-replication n_clones = 1 would inflate sigma_run to sqrt(8.0^2 + 7.0^2/2) = 9.41 pp). "
                "Both Mode 1 (108 independent run-nested clone differentiations, Cov(Y_bar) = sigma_run^2 * I_36) and Mode 2 (75 shared-genotype clones across 25 distinct genotypes, "
                "Cov(Y_bar) = (sigma_clone^2 / n_clones) * Z Z^T + (sigma_eps^2 / (n_clones * n_batches)) * I_36 with exact GLS and OLS sandwich SEs) are reported below."
            ),
            "minimum_detectable_effects_at_80pct_power": dalg["power_and_covariance_summary"],
            "organ_on_chip_stage2_6channel_power_qualification": {
                "n_channels_per_arm": 6,
                "assumed_effect_size_cohens_d": 1.67,
                "two_sided_t_test_df10_ncp_2_887": {
                    "noncentrality_parameter_delta_exact": 2.892525,
                    "power_at_two_sided_alpha_0_05": 0.7413,
                    "power_at_two_sided_alpha_0_10_or_one_sided_alpha_0_05": 0.8518,
                    "power_at_n7_df12_two_sided_alpha_0_05": 0.8178,
                    "power_at_n8_df14_two_sided_alpha_0_05": 0.8736,
                    "min_integer_n_per_arm_for_80pct_power_at_two_sided_alpha_0_05": 7,
                    "n_channels_per_arm_required_for_80pct_power_at_two_sided_alpha_0_05": 8,
                    "organ_on_chip_one_vs_two_sided_power_reconciliation_sonnet55_r2": {
                        "primary_directional_gate_test": "Pre-registered one-sided alpha = 0.05 superiority gate (exact power = 0.8518 at n = 6 channels/arm, df = 10, delta = 2.892525) paired with an automatic adverse-direction kill rule if observed mean shift is adverse (delta_obs <= 0)",
                        "strict_two_sided_alpha_0_05_test": "Two-sided alpha = 0.05 (equivalent to one-sided alpha = 0.025; exact power = 0.7413 at n = 6, 0.8178 at n = 7, and 0.8736 at n = 8 channels/arm on a symmetric 2x4 chip layout; S55R2-S2-08)",
                    },
                    "note": (
                        "In the Stage-2 microfluidic organ-on-a-chip validation (n = 6 independent microfluidic channels per arm, "
                        "two-sample t-test df = 10, exact noncentrality parameter delta = 1.67 * sqrt(3) = 2.892525), n = 6 channels per arm "
                        "achieves 85.18% exact noncentral-t power at a pre-registered one-sided alpha = 0.05 (85.18% at two-sided alpha = 0.10, "
                        "t_crit = 1.81246), whereas two-sided alpha = 0.05 (t_crit = 2.22814) yields 74.13% exact noncentral-t power "
                        "(69.9% under conservative normal-quantile approximation) and requires n = 7 channels per arm (df = 12, exact power = 81.78%) "
                        "or n = 8 channels per arm on symmetric 2x4 chip layouts (df = 14, exact power = 87.36%) to exceed 80% two-sided power "
                        "(Finding 028, AUDIT-REF, & S55R2-S2-08)."
                    ),
                },
            },
            "multiplicity_and_gate_decision_rule": (
                "Co-primary gate endpoints (Gate 1 dynamic GSIS SI >= 3.0 & shut-off <= 3.0 min; Gate 3 NK/macrophage specific lysis <= 10%; Gate 2A/2B autoreactive/allo CTL lysis <= 10%) "
                "are evaluated via hierarchical gatekeeping (Gate 1 viability/maturation gate first, followed by Benjamini-Hochberg FDR control at q = 0.10 across the 6 main effects and 15 two-factor interactions)."
            ),
        },
        "lean_chassis_a_exact_run_id": "HIP_FACTORIAL_RUN_35",
        "lean_chassis_a_replicate_run_id": "HIP_FACTORIAL_RUN_36",
        "lean_chassis_a_coded_levels": dict(zip(HIP_FACTOR_KEYS, HIP_LEAD_CHASSIS_A_CODED)),
        "legacy_lean_chassis_endogenous_promoter_run_ids": ["HIP_FACTORIAL_RUN_17", "HIP_FACTORIAL_RUN_27"],
        "resolution_tier1": "IV",
        "base_factors": ["X1", "X2", "X3", "X4"],
        "defining_generators_tier1": ["X5 = X1 * X2 * X3", "X6 = X2 * X3 * X4"],
        "defining_relation_tier1": "I = X1*X2*X3*X5 = X2*X3*X4*X6 = X1*X4*X5*X6",
        "constant_isogenic_backbone": {
            "cell_source": "Isogenic human Stage 7 tri-hormonal (60% beta, 30% alpha, 10% delta) SC-islet clusters (<= 120 um diameter)",
            "fixed_edits": [
                "CIITA-/- (Locus 2 MHC-II KO)",
                "Biallelic AAVS1-[455bp/0.7kb CBX3-A2UCOE]-EEF1A1-core-iCasp9 (or AAVS1-CAG-iCasp9; Locus 4 pan-lineage safety switch)",
            ],
            "dosing_unit": (
                "Two-Step Lineage Purification & Tri-Hormonal Re-Aggregation Workflow: "
                "Step 1 (Live-Cell Surface Sorting): Stage 7 SC-islet cultures (after Stage 7 week 0 [S7w0] Aurora kinase inhibitor ZM447439 pulse to suppress off-target enterochromaffin proliferation; Balboa et al. 2022 PMID:35241836) "
                "are dissociated and purified without fixation/permeabilization via live-cell surface MACS/FACS into beta-like cells (CD49a+/ITGA1+ live-cell enrichment per Veres et al. 2019 Nature PMID:31068696, which achieves <=80% beta-cell purity with <5% enterochromaffin cells, "
                "co-evaluated with adult beta-cell ectonucleotidase NTPDase3+/ENTPD3+ per Saunders et al. 2019 Cell Metab PMID:30449685 and CD9- negative selection per Li et al. 2020 PMID:33096048), "
                "alpha-like cells (CD26+/DPP4+ / CD49a-, excluding CD49a+ beta and polyhormonal cells, validated by intracellular GCG+ QC aliquot), and delta-like cells (validated by intracellular SST+ QC aliquot). "
                "Step 2 (60:30:10 Pre-Aggregation Input Seeding & Standardized 1,600-Cells/IEQ Volumetric Dosing): Using the unified conversion 1 IEQ (150 um sphere equivalent = 1.77x10^6 um^3) = 1,600 endocrine cells "
                "(matching 500,000 IEQ = 8.0x10^8 cells at curative scale), purified live endocrine subsets are seeded at a pre-aggregation input ratio of 1.2x10^6 beta : 6.0x10^5 alpha : 2.0x10^5 delta cells per 1,000 IEQ batch "
                "(2.0x10^6 total input cells before ~20% AggreWell re-aggregation/handling loss, yielding ~1.6x10^6 viable aggregated endocrine cells per 1,000 IEQ distributed across ~1,950 uniform <=120 um diameter spheroids of ~820 cells/spheroid) and dosed at 1,000 IEQ (1.6x10^6 cells) per assay well. "
                "A matched fixed/permeabilized destructive QC aliquot from each re-aggregated batch is stained intracellularly for NKX6-1+/C-peptide+/INS+ (BD clone R11-560), monohormonal GCG+/C-peptide-, SST+, and intracellular SLC18A1+ (VMAT1+) "
                "so exact viable beta-cell count per 1,000 IEQ is recorded as a continuous covariate alongside per-beta-cell normalized secondary readouts."
            ),
        },
        "factors": factor_defs,
        "alias_structure": {
            "tier1_16run_rank": 14,
            "tier1_16run_main_effects_clear_of_all_2fi": True,
            "tier2_combined_all_main_and_2fi_estimable_full_rank_22_of_22": dalg["all_2fi_estimable"],
            "tier2_combined_with_stage_block_full_rank_23_of_23": dalg["full_rank_with_block"],
            "tier1_main_effect_3fi_aliases": {
                "X1": "X1 + X2*X3*X5 + X4*X5*X6",
                "X2": "X2 + X1*X3*X5 + X3*X4*X6",
                "X3": "X3 + X1*X2*X5 + X2*X4*X6",
                "X4": "X4 + X2*X3*X6 + X1*X5*X6",
                "X5": "X5 + X1*X2*X3 + X1*X4*X6",
                "X6": "X6 + X2*X3*X4 + X1*X4*X5",
            },
            "tier1_two_factor_interaction_alias_chains_resolved_by_tier2": [
                {"chain_id": "2FI_CHAIN_1", "aliased_pairs_in_tier1": "X1*X2 + X3*X5", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "HLA-E x CD47 co-inhibitory NK/macrophage synergy de-aliased from CXCL10::SERPINB9/PD-L1 x TXNIP-/-"},
                {"chain_id": "2FI_CHAIN_2", "aliased_pairs_in_tier1": "X1*X3 + X2*X5", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "HLA-E x CXCL10::SERPINB9/PD-L1 de-aliased from CD47 x TXNIP-/-"},
                {"chain_id": "2FI_CHAIN_3", "aliased_pairs_in_tier1": "X1*X4 + X5*X6", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "Endogenous B2M::HLA-E x TYK2-/- IFN-promoter attenuation de-aliased from TXNIP-/- x RNLS-/- dual metabolic edit"},
                {"chain_id": "2FI_CHAIN_4", "aliased_pairs_in_tier1": "X1*X5 + X2*X3 + X4*X6", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "CD47 x CXCL10::SERPINB9/PD-L1 and TYK2-/- x RNLS-/- de-aliased from X1*X5"},
                {"chain_id": "2FI_CHAIN_5", "aliased_pairs_in_tier1": "X1*X6 + X4*X5", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "HLA-E x RNLS-/- de-aliased from TYK2-/- x TXNIP-/-"},
                {"chain_id": "2FI_CHAIN_6", "aliased_pairs_in_tier1": "X2*X4 + X3*X6", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "B2M-CD47 vs constitutive CD47 x TYK2-/- promoter blunting de-aliased from CXCL10::SERPINB9/PD-L1 x RNLS-/-"},
                {"chain_id": "2FI_CHAIN_7", "aliased_pairs_in_tier1": "X2*X6 + X3*X4", "dealiased_in_combined_design": dalg["all_2fi_estimable"], "biological_interpretation": "Critical epistatic pair CXCL10::SERPINB9/PD-L1 x TYK2-/- de-aliased from CD47 x RNLS-/-"},
            ],
            "co_primary_tyk2_inducible_promoter_interaction_estimands_muse_r3": {
                "co_primary_interactions": ["X1:X4", "X3:X4"],
                "mode1_23col_var_beta_over_sigma2": {
                    "X1:X4": "347/3232 (0.107364)",
                    "X3:X4": "191/808 (0.236386)",
                },
                "mode1_23col_mean_2fi_var_beta_over_sigma2_across_15_pairs": "2097/8080 (0.259530; sum = 12582/3232; vs 6-main-effect mean 137/2424 = 0.056518)",
                "mode1_23col_mean_main_effect_var_beta_over_sigma2_across_6_factors": "137/2424 (0.056518; sum = 137/404 = 1096/3232)",
                "mode1_23col_half_effect_se_pp_at_sigma_run_5p4314": {
                    "X1:X4": 1.7797,
                    "X3:X4": 2.6407,
                    "mean_across_15_2fi_pairs": 2.7669,
                    "arithmetic_mean_across_15_2fi_pairs": 2.5924,
                    "arithmetic_mean_across_15_2fi_pairs_exact": 2.5915,
                    "mean_across_6_main_effects": 1.2912,
                    "arithmetic_mean_across_6_main_effects": 1.2910,
                },
                "mode1_23col_full_contrast_2beta_se_pp_at_sigma_run_5p4314": {
                    "X1:X4": 3.5593,
                    "X3:X4": 5.2814,
                    "mean_across_15_2fi_pairs": 5.5338,
                    "arithmetic_mean_across_15_2fi_pairs": 5.1848,
                    "arithmetic_mean_across_15_2fi_pairs_exact": 5.1831,
                    "mean_across_6_main_effects": 2.5824,
                    "arithmetic_mean_across_6_main_effects": 2.5819,
                },
                "mode1_23col_difference_in_differences_4beta_se_pp_at_sigma_run_5p4314": {
                    "X1:X4": 7.1187,
                    "X3:X4": 10.5629,
                },
                "clustered_covariance_half_effect_se_pp": {
                    "X1:X4": {"mode1_donor_block_se_pp": 1.8872, "mode2_crossed_clone_donor_se_pp": 2.1427},
                    "X3:X4": {"mode1_donor_block_se_pp": 2.8332, "mode2_crossed_clone_donor_se_pp": 2.9847},
                },
                "biological_rationale": (
                    "Elevates beta_{X1:X4} (B2M::HLA-E x TYK2; exact (X^T X)^-1 diagonal 347/3232 = 0.107364, half-effect SE = 1.7797 pp, 2*beta contrast SE = 3.5593 pp, 4*beta DiD SE = 7.1187 pp) "
                    "and beta_{X3:X4} (CXCL10::SERPINB9/PD-L1 x TYK2; exact (X^T X)^-1 diagonal 191/808 = 0.236386, half-effect SE = 2.6407 pp, 2*beta contrast SE = 5.2814 pp, 4*beta DiD SE = 10.5629 pp; "
                    "where 2097/8080 = 0.259530 [sum = 12582/3232, RMS half-effect SE = 2.7669 pp, arithmetic mean = 2.5915 pp exact / 2.5924 pp from 4-decimal rounded diagonals] is the mean variance multiplier across all 15 2FI pairs and 137/2424 = 0.056518 [sum = 137/404 = 1096/3232 across [187/3232, 187/3232, 187/3232, 23/404, 167/3232, 23/404], RMS half-effect SE = 1.2912 pp, arithmetic mean = 1.2910 pp, full-contrast RMS SE = 2.5824 pp; retracting the prior 215/1616 = 0.133045 typo] is the 6-main-effect mean variance multiplier) to co-primary interaction estimands "
                    "in the 36-run HIP-Core regression alongside the 6 main effects, directly testing whether TYK2 loss/hypomorphism blunts ISRE-driven "
                    "HLA-E and SERPINB9/PD-L1 induction under Type-I-IFN-dominant early autoimmune attack (AUDIT-REF, Rule S8, AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                ),
                "isogenic_tyk2_p1104a_by_l3_bfb_hybrid_promoter_satellite_sol61_r4": {
                    "design_id": "STAGE0_TYK2_GENOTYPE_BY_LOCUS3_PROMOTER_3X2_SATELLITE",
                    "factors": {
                        "Factor_1_TYK2_Genotype": ["TYK2_WT", "TYK2_P1104A_HOMOZYGOUS_PRIME_EDIT", "TYK2_KO_BIALLELIC"],
                        "Factor_2_Locus3_Promoter": ["Canonical_Syn_5xISRE_minP", "L3_BF_B_Hybrid_Syn_3xISRE_2xNFkB_hPGK1_Sp1_basalfloor_230"],
                    },
                    "l3_bfb_materialized_construct_binding_kimi31_r4": {
                        "construct_id": "ISLET1_A4_L3_BFB_BASALFLOOR",
                        "donor_length_bp": 4078,
                        "packaged_aav6_length_bp": 4368,
                        "donor_sequence_sha256": "0a5f714d536cfb427039e50966fbc88cc6f7dbec8252f26bf43f32a359017c7e",
                        "promoter_element_composition": "Syn_3xISRE_2xNFkB_enhancer_126bp + hPGK1_Sp1_basal_core_104bp (230 bp driving SERPINB9; paired with 148-bp Syn_5xISRE_minP driving PD-L1)",
                    },
                    "canonical_syn5xisre_comparator_construct_bindings_astra6_r4": {
                        "matched_2afree_comparator": {
                            "comparator_id": "ISLET1_A4_L3_SYN5XISRE_2AFREE_COMPARATOR_3996",
                            "donor_length_bp": 3996,
                            "packaged_aav6_length_bp": 4286,
                            "donor_sequence_sha256": "afae235efdf6bd013f6bb690ac53f7d0ebbf5a6590e441d225debb1bb42f85ba",
                            "duplicated_50mer_windows": 108,
                            "architecture": "Matched 2A-free Locus-3 comparator identical to ISLET1_A4_L3_BFB_BASALFLOOR except the 230-bp L3-BF-B promoter is replaced by the 148-bp Syn_5xISRE_minP element (isolating the pure promoter swap without confounding against tPT2A; carries a 160-bp inter-cistron Syn_5xISRE_minP + Kozak + ATG direct repeat [108 duplicated 50-mers], so paired with the repeat-free 3,799-bp CXCL10p::SERPINB9 + Syn_5xISRE_minP::PD-L1 comparator)",
                        },
                        "canonical_scarfree_3799_comparator": {
                            "construct_id": "ISLET1_A4_L3_SCARFREE_3799",
                            "donor_length_bp": 3799,
                            "packaged_aav6_length_bp": 4089,
                            "donor_sequence_sha256": "3d6bf68e2f0b70911cfea69255e4733db2862ce3baa2b92ae4206561e284bcf7",
                            "pam_shielded_3799bp_donor_sequence_sha256": "13cccea841fbb056c1c5f8ab3f083e6729dd35817436c7788d6b294a6e3508c4",
                            "duplicated_50mer_windows": 4,
                            "inter_unit_promoter_duplications": 0,
                            "architecture": "Canonical Locus-3 scar-free SERPINB9 (endogenous CXCL10 promoter) + Syn_5xISRE_minP -> PD-L1 donor (zero inter-cistron promoter repeat)",
                        },
                    },
                    "replication_and_interaction_estimand_specification_kimi31_r4": {
                        "n_arms": 6,
                        "biological_clones_per_arm_r": 3,
                        "total_independent_clones": 18,
                        "cytokine_regimes_per_clone": 4,
                        "total_clone_by_regime_evaluations": 72,
                        "formal_interaction_estimands": {
                            "DeltaDelta_P1104A_vs_WT_by_Promoter": "((mu[P1104A, L3_BF_B] - mu[P1104A, Syn_5xISRE]) - (mu[WT, L3_BF_B] - mu[WT, Syn_5xISRE])) under Combined_IFN_alpha_1000U_plus_low_IFN_gamma_10U_plus_TNFa_10ng_mL_24h",
                            "DeltaDelta_KO_vs_P1104A_by_Promoter": "((mu[KO, L3_BF_B] - mu[KO, Syn_5xISRE]) - (mu[P1104A, L3_BF_B] - mu[P1104A, Syn_5xISRE])) under IFN_alpha_1000U_mL_24h_early_TypeI_flare",
                            "mechanistic_loss_of_typeI_ifn_control_cell": "TYK2_KO_BIALLELIC x Canonical_Syn_5xISRE_minP serves as the pre-registered negative-induction control for pure Type-I IFN (IFN-alpha) stimulation, isolating NF-kB + basal-floor rescue on L3-BF-B.",
                        },
                    },
                    "cytokine_regimes_tested": [
                        "Unstimulated_Hour0",
                        "IFN_alpha_1000U_mL_24h_early_TypeI_flare",
                        "IFN_gamma_100U_mL_24h_CTL_NK_synapse",
                        "Combined_IFN_alpha_1000U_plus_low_IFN_gamma_10U_plus_TNFa_10ng_mL_24h",
                    ],
                    "pass_criteria": (
                        "Directly resolves the TYK2 P1104A vs KO and L3-BF-B (Syn_3xISRE_2xNFkB + hPGK1_Sp1_basal_core, 4,078-bp donor SHA-256 0a5f714d536cfb427039e50966fbc88cc6f7dbec8252f26bf43f32a359017c7e) vs matched 2A-free Syn_5xISRE (3,996-bp donor SHA-256 afae235efdf6bd013f6bb690ac53f7d0ebbf5a6590e441d225debb1bb42f85ba and repeat-free 3,799-bp CXCL10p::SERPINB9 + Syn_5xISRE_minP::PD-L1 comparator SHA-256 3d6bf68e2f0b70911cfea69255e4733db2862ce3baa2b92ae4206561e284bcf7 / PAM-shielded 13cccea841fbb056c1c5f8ab3f083e6729dd35817436c7788d6b294a6e3508c4) contrast across r = 3 independent clones/arm (18 clones total) outside the 36-run main factorial (which tests TYK2 KO vs WT and canonical Locus 3 vs KO): "
                        "requires Hour-0 [SERPINB9]_0 >= 0.30 uM and >= 70% shield induction under combined IFN-alpha + low IFN-gamma + TNF-alpha, while testing whether human SC-beta TYK2 P1104A retains WT-like IFN-alpha/beta sensitivity at high cytokine concentrations as observed in human primary leukocytes (Boisson-Dupuis et al. 2018 PMID:30578352; Howard et al. 2026 PMID:42268925; AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                    ),
                },
            },
            "endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_gate0a_audit_astra6_r4": {
                "primary_citations": ["PMID:40133488", "PMID:37820008", "PMID:38508139", "PMID:29345283", "PMID:39508879", "PMID:30951669"],
                "gate_0a_tsp1_and_protease_shedding_panel": (
                    "Propagates Layer 3 endogenous 5-TM CD47 TSP-1 sensitivity (Ghimire et al. 2023 PMID:37820008; Kale et al. 2025 PMID:40133488; Xu et al. 2024 PMID:38508139) "
                    "and inflammatory protease ectodomain shedding (Hira-Miyazawa et al. 2018 PMID:29345283) directly into the Stage-0 protocol: "
                    "(1) Decision Experiment #5 adds a recombinant human TSP-1 (1-5 ug/mL) perifusion challenge and an exploratory endogenous CD47 5-TM KO + pure CD47-IgV-GPI replacement satellite arm; "
                    "(2) Gate 0A-surface / Gate 0C-kinetic adds a 4-h neutrophil elastase (100 nM) + recombinant MMP-7/MMP-9/MMP-13 challenge requiring >= 80% surface epitope retention by QuantiBRITE flow and < 10% soluble sPD-L1/sCD47 decoy accumulation by ELISA (AUDIT-REF & AUDIT-REF)."
                ),
                "cell_free_product_effluent_pdl1_functional_screen_iip_0005": {
                    "iip_id": "IIP-0005",
                    "decision_id": "DEC-083",
                    "affects_tests": ["Gate_2B", "Gate_4"],
                    "affects_designs": ["Design_A", "Design_B", "Design_X"],
                    "primary_citations": ["PMID:39508879", "PMID:30951669", "PMID:29345283"],
                    "related_graveyard_rules": ["GY-008", "GY-001", "GY-010"],
                    "scope_and_non_conflation_note": (
                        "Supplements post-protease soluble sPD-L1/sCD47 abundance ELISA (Hira-Miyazawa et al. 2018 PMID:29345283) and Gate 0d graft-antigen/exosome recipient-APC assays "
                        "with a staged, bounded screen of minimally processed cell-free medium conditioned by the product's cells (distinct from Readout 5 CaO2 hydrolysis effluent) on an independently stimulated PD-1-responsive assay "
                        "(Rao et al. 2025 PMID:39508879; Poggio et al. 2019 PMID:30951669). Leaves all construct sequences, ranking weights, Stage-0 COMPLETE/PASS rules, and acceptable shedding percentages unchanged."
                    ),
                    "staged_protocol_clauses": [
                        "1. Dedicated post-washout conditioned-medium harvest from the materialised lead and a named matched parental or package-minus comparator under resting and one planned inflammatory condition, documenting viable input, conditioning time/volume, and maintained product PD-L1 state without reusing protease-challenge supernatant or antigen lysate.",
                        "2. Minimally processed cell-free effluent with a verified cell-free collection boundary and qualified clarification/filtration and matrix recovery for biologically relevant carrier forms so a fixed 0.45-um filter does not silently remove membrane/vesicle carriers.",
                        "3. One validated independently stimulated PD-1-responsive assay with stimulation/medium controls, matched control effluent, responder viability, active anti-PD-L1 blockade vs isotype (measured in blanks, control effluent, and lead effluent), and a form-relevant functional positive reference recovered through the actual processing train, reported as functional effects with uncertainty across independent preparations rather than inferred ng/mL PD-L1 equivalents.",
                        "4. Inflammatory cytokine carryover quantification and lot-matched cytokine-spiked unconditioned blanks (plus neutralisation or source/receiver controls where needed), reporting anti-PD-L1 rescue as PD-L1-pathway-dependent activity until receiver/accessory-cell PD-L1 induction and product-side source attribution are resolved.",
                        "5. Staged follow-up rule: a reproducible, source-attributed positive directs one decision-limiting follow-up (primary-human responder, carrier attribution, or extra-graft concentration-time bound); a sensitive controlled negative deprioritises further work within the tested product, state, processing method, and exposure bound; failed recovery, insensitive responder, carryover, or induction decay remains unresolved and triggers targeted calibration."
                    ],
                },
            },
            "dealiasing_note": dalg["dealiasing_note"],
        },
        "whole_product_clonogenic_safety_model": whole_product_safety_model,
        "primary_readouts": [
            "Stage 0->7 Absolute Endocrine Yield & Lineage Purity: Absolute monohormonal NKX6-1+/C-peptide+ beta-cell yield per 10^6 starting hPSCs, % GCG+ alpha, % SST+ delta, % SLC18A1+ off-target enterochromaffin cells (<5% release gate after FACS re-aggregation into <=120 um 60:30:10 spheroids), and batch differentiation failure rate across constitutive TYK2-/- vs TYK2+/+ vs Stage 6/7 TYK2i (Chandra et al. 2022 PMID:36289205 ~32% INS+ reduction gate; Gate 1)",
            "Hour-0 and +IFN-gamma QuantiBRITE molecules/cell for HLA-E and CD47 (comparing bicistronic tPT2A vs 2A-free bidirectional vs CD47-IgV-GPI / CD47-dCT / CD47-WT) on Stage 7 CD49a+/NTPDase3+ beta-cells and iECs (Gate 3)",
            "Primary human innate immune evasion across two unconfounded sub-panels (3A Myeloid Phagocytosis Panel: CD14+ monocyte-derived macrophages genotyped for SIRPA V1/V1, V1/V2, and V2/V2 [n >= 3 donors per genotype, pHrodo-Red 4h engulfment < 5.0% confirmatory release gate, <= 8.0% Tier-1 36-run screening gate]; 3B Sorted Primary NK Subset Panel: purified CD3- CD56dim NK subsets stratified into HCMV-/NKG2A-high, NKG2A- KIR2DL1+/C2 & KIR2DL2/3+/C1 licensed, and HCMV+/NKG2C-high [>=25% NKG2C+ of CD56dim] donors [n >= 4 donors per subset for confirmatory release, >= 3 Tier-1 screen minimum; 4h specific lysis <= 8.0% in NKG2A+NKG2C- and <= 10.0% across all subsets, CD107a degranulation <= 10.0% confirmatory / <= 12.0% Tier-1 screen]; AUDIT-REF & AUDIT-REF; Gate 3)",
            "Autoreactive T1D CD8+ TCR (PPI/IGRP/GAD65) and allogeneic PBMC survival + HLA-DQ2/DQ8 DC indirect presentation bystander assay + cytokine-stimulated beta-cell Fas (CD95) incremental apoptosis screen (Delta_TUNEL+ <= +5.0 pp above cytokine-matched SA-THBD-only control AND anti-Fas ZB4 neutralizing mAb rescue <= 5.0 pp) under SA-FasL/SA-THBD microgels (Gate 2A/2B)",
            "Dynamic perifusion GSIS stimulation index, 16.7 -> 2.8 mM glucose 90% shut-off time (<3 min gate), and Seahorse basal/maximal OCR (Gate 1)",
            "Satellite 24h and 72h 10 nM AP1903 suicide-switch eradication depth (Two-Tier Gate: Tier-1 36-run bulk potency screen requiring >= 99.9% 24h kill and >= 4-log 72h kill across all 36 runs; Tier-2 lead/MCB clonogenic escape qualification across N_eff = N_clonogenic,input * eta_detection >= 2.996e8 detectable clonogenic opportunities calibrated with 1-, 10-, and 100-cell AP1903-resistant spike-in recovery controls; AUDIT-REF & 012; Gate 5)",
        ],
        "tier1_screening_16_runs": factorial_runs[:16],
        "tier2_dealiasing_augmentation_20_runs": factorial_runs[16:],
        "runs": factorial_runs,
    }
    protocols_dir = PROJECT_ROOT / "data" / "protocols"
    protocols_dir.mkdir(parents=True, exist_ok=True)

    direct_logit_sensitivity = compute_direct_output_logit_sensitivity_analysis(
        target_profiles, all_combos, studies=studies, omics=omics
    )

    wetlab_calibration_link_status = dict(_WETLAB_CALIBRATION_STATUS)

    raw_edit_rows = [edit for study in studies for arm in study.get("arms", []) for edit in arm.get("gate2_3_4_genetic_and_epitope_edits", [])]
    expanded_edit_genes = [gene for edit in raw_edit_rows for gene in expand_gene_symbol_tokens(edit.get("gene_symbol_hgnc", ""))]
    aux_values = {key: value["parameter_value"] for key,value in COMPLETE_GATE_COEFFICIENT_REGISTRY["auxiliary_evaluator_hyperparameters"].items()}
    reg_entries = {f"{gate}.{key}": value for gate,group in COMPLETE_GATE_COEFFICIENT_REGISTRY.items() if gate.startswith("Gate") for key,value in group.items() if isinstance(value,dict) and "logit_value" in value}
    reference = json.loads((PROJECT_ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text())
    promoter = reference["ensembl_verified_homology_arms_and_elements"]["EEF1A1_core_promoter_230bp"]["sequence"].upper()
    reverse_promoter = promoter.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    fasta_records = {}
    active_fasta = PROJECT_ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta"
    for chunk in active_fasta.read_text().split(">"):
        if chunk.strip():
            lines = chunk.splitlines()
            fasta_records[lines[0]] = "".join(lines[1:]).strip().upper()
    promoter_occurrences = [{"construct_id": cid, "forward_occurrences": seq.count(promoter),
                             "reverse_complement_occurrences": seq.count(reverse_promoter)}
                            for cid,seq in fasta_records.items() if promoter in seq or reverse_promoter in seq]
    layer1_audit_ledger = {
        "method": "Computed from this evaluation's source studies, graph, registry, exact design and diagnostics; not an audit closure assertion",
        "source_counts": {"studies": len(studies), "arms": sum(len(study.get("arms", [])) for study in studies),
                          "atomic_effects": sum(len(study.get("atomic_effects", [])) for study in studies),
                          "raw_edit_rows": len(raw_edit_rows), "expanded_edit_gene_occurrences": len(expanded_edit_genes)},
        "corpus_eligibility": graph.get("corpus_epistemic_and_fidelity_diagnostics", {}),
        "profile_counts": {"genes": len(target_profiles), "candidate_target_profiles_with_typed_numeric_lrr": sum(profile.get("signed_lrr_weight_total",0)>0 for profile in target_profiles.values()), "independently_adjudicated_direction_effects": 0, "direction_method": "rule_based_heuristic; no independent adjudication dataset supplied"},
        "native_eef1a1_promoter": {"length_bp": len(promoter), "cpg_dinucleotides": promoter.count("CG"),
                                    "gc_fraction": (promoter.count("G")+promoter.count("C"))/len(promoter),
                                    "sha256": hashlib.sha256(promoter.encode()).hexdigest(),
                                    "construct_occurrences": promoter_occurrences,
                                    "total_exact_or_rc_occurrences": sum(row["forward_occurrences"]+row["reverse_complement_occurrences"] for row in promoter_occurrences)},
        "registry": {"source_literal_offset_terms": sum(len(re.findall(r"(?<![\w.])0\.08\s*\*", value.get("expression", ""))) for value in reg_entries.values()),
                     "net_gate_offset_entries": len(reg_entries), "gate_offset_paths": sorted(reg_entries),
                     "auxiliary_groups": len(aux_values), "auxiliary_scalar_leaves": len(_numeric_scalar_paths(aux_values)),
                     "auxiliary_scalar_paths": _numeric_scalar_paths(aux_values), "engineering_priors": len(ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE)},
        "search": {"all_stacks": len(all_combos), "feasible_stacks": sum(c.get("is_feasible_candidate",True) for c in all_combos),
                   "top_tier_stacks": direct_logit_sensitivity["feasible_subset_and_within_top_tier_variance_decomposition_sol_v80_r2_14"]["top_tier_code_defined_n_stacks"]},
        "design": {"run_rows": len(factorial_runs), "unique_genotypes": dalg["n_distinct_genotypes"],
                   "factor_column_sums": col_sums, "factor_plus_one_counts": plus_one_counts},
        "variance_diagnostics": direct_logit_sensitivity["uniform_grid_pooled_pediatric_gate_logit_variance_decomposition"],
    }

    champion_challenger_portfolio = build_champion_challenger_portfolio(ranked_architectures, all_combos, wetlab_calibration_link_status)

    p_esc_ref = 1e-8
    p_20q_ref = 0.0
    p_20q_mosaic_lod_bound = 0.01
    p_meta_ref = 1e-6
    f_ret_ref = 0.99
    p_harm_ref = 1.0 - (1.0 - p_esc_ref) * (1.0 - p_20q_ref) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    p_esc_8psc_dose = 1.0 - math.exp(-8.0 * p_esc_ref)
    p_harm_8psc_pre_retrieval = 1.0 - (1.0 - p_esc_8psc_dose) * (1.0 - p_20q_ref) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    p_harm_8psc_with_99pct_retrieval = 1.0 - (1.0 - p_esc_8psc_dose * (1.0 - f_ret_ref)) * (1.0 - p_20q_ref) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    # Exact per-residual-PSC escape rate under 1% 20q11.21 trisomic mosaicism (q_esc = 0.01 * 1e-4 + 0.99 * 1e-8 = 1.0099e-6) without double-counting the 1e-8 dual-purge factor already reflected in N_psc_post_purge = 8.0 (AUDIT-REF & MW-R6-S2-01; DEC-052/053):
    q_esc_1pct_mosaic_20q11 = p_20q_mosaic_lod_bound * 1e-4 + (1.0 - p_20q_mosaic_lod_bound) * 1e-8
    p_esc_8psc_if_1pct_mosaic_20q11 = 1.0 - math.exp(-8.0 * q_esc_1pct_mosaic_20q11)
    p_harm_8psc_if_1pct_mosaic_20q11_true_pre_retrieval_no_explant_discount = 1.0 - (1.0 - p_esc_8psc_if_1pct_mosaic_20q11) * (1.0 - p_meta_ref)
    p_harm_8psc_if_1pct_mosaic_20q11_and_1e4_ap1903_resist = 1.0 - (1.0 - p_esc_8psc_if_1pct_mosaic_20q11) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    p_esc_only_8psc_if_1pct_mosaic_20q11_with_99pct_retrieval = p_esc_8psc_if_1pct_mosaic_20q11 * (1.0 - f_ret_ref)
    p_harm_8psc_if_1pct_mosaic_20q11_with_99pct_retrieval = 1.0 - (1.0 - p_esc_8psc_if_1pct_mosaic_20q11 * (1.0 - f_ret_ref)) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    p_harm_8psc_simple_1pct_residual_of_8p0892e_6 = p_harm_8psc_if_1pct_mosaic_20q11_and_1e4_ap1903_resist * (1.0 - f_ret_ref)
    # Whole-product post-purge burden (9.335 PSCs across 9.335e8 total cells at f_PSC,eff = 1e-8) and pre-Tier-A intermediate whole-product burden (93.35 PSCs across 9.335e8 total cells at f_PSC = 1e-7; R9-S1-L1B-01 & AUDIT-REF [DEC-074/075]):
    p_esc_9p335psc_if_1pct_mosaic_20q11 = 1.0 - math.exp(-9.335 * q_esc_1pct_mosaic_20q11)
    p_harm_9p335psc_if_1pct_mosaic_20q11_true_pre_retrieval = 1.0 - (1.0 - p_esc_9p335psc_if_1pct_mosaic_20q11) * (1.0 - p_meta_ref)
    p_harm_9p335psc_if_1pct_mosaic_20q11_with_99pct_retrieval = 1.0 - (1.0 - p_esc_9p335psc_if_1pct_mosaic_20q11 * (1.0 - f_ret_ref)) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    q_esc_1pct_mosaic_10x_resist = p_20q_mosaic_lod_bound * 1e-7 + (1.0 - p_20q_mosaic_lod_bound) * 1e-8
    p_esc_93p35psc_10000x_mosaic = 1.0 - math.exp(-93.35 * q_esc_1pct_mosaic_20q11)
    p_harm_93p35psc_10000x_mosaic_pre_retrieval = 1.0 - (1.0 - p_esc_93p35psc_10000x_mosaic) * (1.0 - p_meta_ref)
    p_harm_93p35psc_10000x_mosaic_with_99pct_retrieval = 1.0 - (1.0 - p_esc_93p35psc_10000x_mosaic * (1.0 - f_ret_ref)) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    p_esc_93p35psc_10x_mosaic = 1.0 - math.exp(-93.35 * q_esc_1pct_mosaic_10x_resist)
    p_harm_93p35psc_10x_mosaic_pre_retrieval = 1.0 - (1.0 - p_esc_93p35psc_10x_mosaic) * (1.0 - p_meta_ref)
    p_harm_93p35psc_10x_mosaic_with_99pct_retrieval = 1.0 - (1.0 - p_esc_93p35psc_10x_mosaic * (1.0 - f_ret_ref)) * (1.0 - p_meta_ref * (1.0 - f_ret_ref))
    whole_product_safety_model["unconditional_patient_level_tumor_harm_and_two_stage_psc_purge_model_sonnet55_r5"] = {
        "finding_ids": ["AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "MW-R5-S1-02", "AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "MW-R6-S2-01", "MW-R6-S1-02", "R9-S1-L1B-01", "AUDIT-REF", "R9-S3-L1-01", "R9-S3-L1-02"],
        "provenance": "Audit Convergence Round 5 Step 5, Audit Full Round 5 Step 6, AUDIT-REF Round 6 Step 2, Audit Round 6 Step 4, Audit Full Round 6 Step 5, GPT AUDIT-REF Audit Round 9 Step 2 & Audit Round 9 Step 3 (DEC-047/048/049/050/052/053/074/075; Tang 2011 PMID:21841799; Avery 2013 PMID:24286026; Veres 2019 PMID:31068696; Liang 2018 PMID:30429614; Straathof 2005 PMID:15728125)",
        "evidence_class": "STAGE0_GATE_CHECKLIST_AND_CONDITIONAL_SCENARIO_BOUND_NOT_RANKER_SCORE_INPUT",
        "model_conditioning_scope_sol61_r6": "conditional_scenario_bound_with_mosaic_20q11_and_untriggered_intragraft_terms",
        "conditional_scenario_bound_not_unconditional_correlated_selection_disclaimer_muse_r6": (
            "This equation is a conditional scenario bound under parameterized post-gate residual fractions; "
            "because 20q11.21/BCL2L1 amplification, partial pluripotency retention, and iCasp9/AP1903 apoptosis resistance "
            "are mechanistically coupled along the mitochondrial MOMP/XIAP axis (Avery et al. 2013 PMID:24286026), "
            "Stage-0 Gate 0S directly measures the coupled joint survival probability empirically via the post-stress 20q11.21 outgrowth assay "
            "(obs_post_stress_20q11_bcl2l1_outgrowth_fold <= 1.0) and the sorted hazard-state proliferative clonogenic assay "
            "(0/299,573,228 colonies under pulsed human-PK AP1903 monotherapy; MW-R6-S2-01 & MW-R6-S1-02 [DEC-053])."
        ),
        "unconditional_patient_level_harm_equation": (
            "P_harm_conditional_post_ablation = 1 - (1 - P_escape_given_AP1903) * (1 - p_20q11_bclxl_resistant_clone) * (1 - p_untriggered_extragraft_metastasis * (1 - f_en_bloc_retrieval_completeness)); "
            "Extended dose-level scenario bound: P_harm_dose = 1 - (1 - P_untriggered_intragraft_teratoma * (1 - P_clinical_detection_and_ablation)) * (1 - P_triggered_post_ablation_escape_with_mosaic_20q11) * (1 - p_untriggered_extragraft_metastasis * (1 - f_en_bloc_retrieval_completeness))"
        ),
        "computed_reference_patient_harm_upper_bound": round(p_harm_ref, 12),
        "computed_dose_level_8psc_patient_harm_pre_retrieval": round(p_harm_8psc_pre_retrieval, 12),
        "computed_dose_level_8psc_patient_harm_with_99pct_en_bloc_retrieval": round(p_harm_8psc_with_99pct_retrieval, 12),
        "computed_dose_level_8psc_escape_if_1pct_mosaic_20q11_pre_retrieval": round(p_esc_8psc_if_1pct_mosaic_20q11, 10),
        "computed_dose_level_8psc_harm_if_1pct_mosaic_20q11_true_pre_retrieval_no_explant_discount": round(p_harm_8psc_if_1pct_mosaic_20q11_true_pre_retrieval_no_explant_discount, 10),
        "computed_dose_level_8psc_harm_with_1pct_mosaic_20q11_pre_retrieval_bound": round(p_harm_8psc_if_1pct_mosaic_20q11_true_pre_retrieval_no_explant_discount, 10),
        "computed_dose_level_8psc_harm_with_1pct_mosaic_20q11_subclone_bound": round(p_harm_8psc_if_1pct_mosaic_20q11_and_1e4_ap1903_resist, 10),
        "computed_dose_level_8psc_escape_only_with_1pct_mosaic_20q11_and_99pct_en_bloc_retrieval": round(p_esc_only_8psc_if_1pct_mosaic_20q11_with_99pct_retrieval, 12),
        "computed_dose_level_8psc_harm_with_1pct_mosaic_20q11_and_99pct_en_bloc_retrieval": round(p_harm_8psc_if_1pct_mosaic_20q11_with_99pct_retrieval, 12),
        "computed_post_retrieval_1pct_residual_8psc_harm_with_1pct_mosaic_20q11_including_1e8_metastasis_floor": round(p_harm_8psc_if_1pct_mosaic_20q11_with_99pct_retrieval, 12),
        "computed_dose_level_8psc_harm_simple_1pct_residual_of_8p0892e_6": round(p_harm_8psc_simple_1pct_residual_of_8p0892e_6, 12),
        "computed_post_retrieval_1pct_residual_8psc_harm_with_1pct_mosaic_20q11_with_retrieval_attenuated_metastasis": round(p_harm_8psc_simple_1pct_residual_of_8p0892e_6, 12),
        "computed_whole_product_9p335psc_harm_with_1pct_mosaic_20q11_true_pre_retrieval": round(p_harm_9p335psc_if_1pct_mosaic_20q11_true_pre_retrieval, 10),
        "computed_whole_product_9p335psc_harm_with_1pct_mosaic_20q11_and_99pct_en_bloc_retrieval": round(p_harm_9p335psc_if_1pct_mosaic_20q11_with_99pct_retrieval, 12),
        "computed_intermediate_93p35psc_10000x_mosaic_harm_pre_retrieval": round(p_harm_93p35psc_10000x_mosaic_pre_retrieval, 10),
        "computed_intermediate_93p35psc_10000x_mosaic_harm_with_99pct_retrieval": round(p_harm_93p35psc_10000x_mosaic_with_99pct_retrieval, 10),
        "computed_intermediate_93p35psc_10x_mosaic_harm_pre_retrieval": round(p_harm_93p35psc_10x_mosaic_pre_retrieval, 10),
        "computed_intermediate_93p35psc_10x_mosaic_harm_with_99pct_retrieval": round(p_harm_93p35psc_10x_mosaic_with_99pct_retrieval, 12),
        "pre_vs_post_retrieval_metastasis_term_reconciliation_muse_r6": (
            "True pre-retrieval dose-level harm without any retrieval discount on extragraft dissemination (under E[N_PSC,endo] = 8.0 and 1% 20q11.21 mosaic at 10,000x AP1903 resistance, q_esc = 0.01 * 1e-4 + 0.99 * 1e-8 = 1.0099e-6) is "
            "1 - (1 - 8.0792e-6) * (1 - 1e-6) = 9.0792e-6; when the 99% en-bloc retrieval factor (1 - f_ret = 0.01) is applied to the "
            "1e-6 extragraft dissemination term (1e-8 residual dissemination), pre-explant intragraft escape + residual dissemination is "
            "1 - (1 - 8.0792e-6) * (1 - 1e-8) = 8.0892e-6, post-99%-explant intragraft escape alone is 8.0792e-6 * 0.01 = 8.0792e-8, "
            "and combined post-99%-explant escape plus residual dissemination is 1 - (1 - 8.0792e-8) * (1 - 1e-8) = 9.0792e-8 "
            "(whereas 0.01 * 8.0892e-6 = 8.0892e-8 would double-discount the 1e-8 dissemination term by (1 - f_ret)^2; at whole-product E[N_PSC,whole] = 9.335, P_pre = 1.0427e-5 and P_post = 1.0427e-7, and at the illustrative pre-Tier-A intermediate whole-product bound f_PSC = 1e-7 across 9.335e8 total cells [E[N_PSC,whole] = 93.35, whereas the 1e7-cell sub-lot screen sample itself has expectation 1.0 PSC], P_pre = 9.5270e-5 / P_post = 9.5270e-7 under 10,000x mosaic resistance [q_esc = 1.0099e-6] or P_pre = 2.0175e-6 / P_post = 2.0175e-8 under 10x mosaic resistance [q_esc = 0.01 * 1e-7 + 0.99 * 1e-8 = 1.09e-8]; MW-R6-S2-01, R9-S1-L1B-01 & AUDIT-REF [DEC-053/DEC-074/DEC-075])."
        ),
        "reference_bound_inputs": {
            "P_escape_given_AP1903_per_clonogenic_unit": p_esc_ref,
            "P_escape_given_AP1903_across_8_residual_pscs": round(p_esc_8psc_dose, 12),
            "p_20q11_bclxl_resistant_clone_post_ddpcr_gate": p_20q_ref,
            "p_20q11_subclonal_mosaic_single_cell_fish_ddpcr_lod_ceiling": p_20q_mosaic_lod_bound,
            "p_untriggered_extragraft_metastasis": p_meta_ref,
            "f_en_bloc_retrieval_completeness": f_ret_ref,
        },
        "two_stage_orthogonal_psc_purge_requirement": {
            "raw_flow_cytometry_lod_floor": 1e-4,
            "naive_residual_psc_count_at_1e4_lod_in_8e8_dose": 80000,
            "stage1_anti_ssea5_tra160_immunomagnetic_depletion_log10": 3.0,
            "stage2_ym155_birc5_survivin_selective_purge_log10": 3.0,
            "theoretical_independent_6log_residual_psc_fraction": 1e-10,
            "theoretical_independent_6log_residual_psc_per_8e8_dose": 0.08,
            "independence_qualification_note": "Theoretical independent multiplicative 6-log depletion (10^-4 * 10^-6 = 10^-10, 0.08 PSCs/dose) is bounded by a conservative 4-log net effective release gate ceiling (10^-4 -> 10^-8, <= 8.0 PSCs/dose) to account for potential partial correlation between surface marker retention and survivin inhibitor sensitivity, verified empirically via the 14-day mTeSR amplification + OCT4/NANOG/TDGF1 ddPCR release gate (obs_residual_psc_effective_fraction <= 1e-8 across obs_residual_psc_effective_trials >= 299,573,228).",
            "orthogonal_combined_residual_psc_fraction_ceiling": 1e-8,
            "expected_residual_psc_per_8e8_dose_after_two_stage_purge": 8.0,
            "lot_release_ddpcr_oct4_nanog_tdgf1_after_14d_mtesr_amplification_ceiling": 1e-8,
            "wetlab_contract_endpoint": "obs_residual_psc_effective_fraction",
            "wetlab_contract_denominator_endpoint": "obs_residual_psc_effective_trials",
        },
        "chr20q11_21_bcl2l1_amplicon_surveillance_gate": {
            "recurrent_psc_hotspot": "chr20q11.21 (ID1 / BCL2L1 [BCL-xL] / HM13; Avery 2013 PMID:24286026)",
            "mechanistic_coupling_to_icasp9": "BCL-xL amplification inhibits mitochondrial MOMP amplification and cooperates with hypoxic XIAP to raise the AP1903 dimerization threshold in residual proliferating cells",
            "mandatory_mcb_and_stage7_ddpcr_cnv_gate": "1.85 <= obs_20q11_bcl2l1_cnv_ratio <= 2.15 (95% CI within [1.85, 2.15]) plus 0/300 G-banded + Optical Genome Mapping structural variants and <= 1% subclonal mosaicism by 1,000-interphase-nuclei 20q11.21/BCL2L1 FISH or single-cell ddPCR",
            "wetlab_contract_endpoint": "obs_20q11_bcl2l1_cnv_ratio",
        },
        "hazard_state_proliferative_survivor_clonogenic_gate": {
            "protocol": (
                "Sort >= 3.0e8 Stage-7 cells for Ki67+ / SOX2+ / OCT4+ / SLC18A1+ hazard-state fractions, challenge under combined 1% O2 + IFN-g/TNF-a/IL-1b + pulsed human-PK AP1903 "
                "(2 x 2h pulses of 10-50 nM free AP1903 monotherapy, requiring obs_ap1903_monotherapy_pulsed_pk_log10_kill >= 4.0 and 0/299,573,228 clonogenic survivors WITHOUT birinapant; "
                "parallel 100 nM birinapant arm run strictly as an exploratory mechanistic XIAP-depletion control), and culture for 28 days in mTeSR1 + Matrigel (AUDIT-REF & AUDIT-REF [DEC-050])"
            ),
            "mandatory_pass_criterion": "obs_hazard_state_prolif_survivor_count_k == 0 colonies across obs_hazard_state_prolif_effective_trials >= 299,573,228 effective clonogenic trials under AP1903 monotherapy (obs_ap1903_monotherapy_pulsed_pk_log10_kill >= 4.0)",
            "wetlab_contract_endpoint": "obs_hazard_state_prolif_survivor_count_k",
            "wetlab_contract_denominator_endpoint": "obs_hazard_state_prolif_effective_trials",
            "wetlab_contract_monotherapy_kill_endpoint": "obs_ap1903_monotherapy_pulsed_pk_log10_kill",
        },
        "en_bloc_retrieval_and_cell_retention_gate": {
            "mandatory_preclinical_retrieval_completeness": "obs_explant_en_bloc_recovery_fraction >= 0.99 (lower 95% CI >= 0.99 across vascularized Day-90 omental explants with obs_vascularized_explant_en_bloc_recovery_verified == True and < 10^-6 extra-capsular Alu/transgene qPCR shedding)",
            "wetlab_contract_endpoint": "obs_explant_en_bloc_recovery_fraction",
        },
        "antiviral_sanctuary_and_cell_intrinsic_viral_fuse_gate_0v_muse_r6": {
            "finding_ids": ["MW-R6-S1-01"],
            "decision_id": "DEC-053",
            "primary_citations": ["PMID:19266182", "PMID:23974921"],
            "lead_product_tyk2_wt_policy": (
                "All lead clinical candidates (Embodiment A* ISLET1_A4_ARM1B, Arm 1X ISLET1_A4_X_MOL_RNLS_WT, and composite challenger ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500) "
                "retain wild-type TYK2 (TYK2-WT) to preserve full IFNAR1/2-STAT1/2 Type-I interferon and IL-12/IL-23 antiviral innate immunity; "
                "optional TYK2 P1104A (rs34536443) satellite variants are strictly quarantined from clinical promotion unless they clear Stage-0 Gate 0V."
            ),
            "stage0_gate_0v_live_viral_challenge_protocol": (
                "Challenge engineered SC-islet organoids at MOI 0.01 and 0.1 with diabetogenic Coxsackievirus B1/B3/B4 (CVB1/CVB3/CVB4, CXADR-tropic), human cytomegalovirus (HCMV), and Epstein-Barr virus (EBV) "
                "for 28 days in the presence of NK cells and autologous/allogeneic PBMCs; require >= 3.0 log10 reduction in viral genome copies by Day 14 "
                "(obs_gate0v_cvb_day14_log10_viral_reduction >= 3.0), 0% persistent VP1+ or IE1+ reservoir cells at Day 28 "
                "(obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct <= 0.0%), intact cytosolic RNA-sensing MDA5 (IFIH1) / RIG-I (DDX58) / MAVS pathway activation, and bystander uninfected beta-cell GSIS SI >= 2.5."
            ),
            "cell_intrinsic_mda5_rigi_mavs_intactness_gate": (
                "Verify intact cell-intrinsic MDA5 (IFIH1) / RIG-I (DDX58) / MAVS -> IRF3/IFN-beta/IFN-lambda1 induction upon poly(I:C) or CVB1/CVB4 challenge so enteroviral abortive clearance operates autonomously upstream of MHC-I/CD8+ CTLs."
            ),
            "cell_intrinsic_viral_fuse_challenger": {
                "challenger_id": "VIRAL_FUSE_PKR_OAS_ABORTIVE_CHALLENGER",
                "mechanism": (
                    "Cell-intrinsic viral-replication-triggered abortive suicide fuse (dsRNA-sensing PKR/OAS1-RNase L or IFIH1/MAVS-inducible caspase-9/GSDME circuit) "
                    "that selectively aborts virally infected cells before virion assembly and paracrine spread without requiring MHC-I/CD8+ CTL recognition."
                ),
            },
            "wetlab_contract_endpoints": [
                "obs_gate0v_cvb_day14_log10_viral_reduction",
                "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct",
            ],
        },
        "coupled_shield_switch_coinheritance_and_post_stress_outgrowth_gate_muse_r6": {
            "finding_ids": ["MW-R6-S1-02"],
            "decision_id": "DEC-053",
            "primary_citations": ["PMID:24286026", "PMID:21841799", "PMID:15728125"],
            "co_inherited_locus1_surface_eraser_tag_challenger": (
                "Co-inheritance challenger coupling an antibody-depletable compact surface tag (RQR8 or tEGFR) or bicistronic CD47-P2A-iCasp9 / orthogonal second suicide switch directly into the "
                "Locus-1 B2M immune-shield bicistron so that promoter-level epigenetic silencing or whole-locus LOH obligately co-silences HLA-E/CD47 "
                "and exposes the cell to immediate innate NK/macrophage clearance, while independent AAVS1 iCasp9 is retained alongside because downstream coding frameshift or P2A ribosome skipping (~5-15% incomplete cleavage) does not abrogate upstream HLA-E/CD47 translation."
            ),
            "dual_orthogonal_eraser_requirement": (
                "Require two mechanistically independent elimination modalities (AAVS1 iCasp9/AP1903 monotherapy + orthogonal cell-surface antibody-mediated CDC/ADCC "
                "or HSV-TK/ganciclovir backup), each independently achieving >= 4-log10 depletion under combined 1% O2 hypoxia + 20q11.21/BCL2L1 trisomy."
            ),
            "post_stress_20q11_outgrowth_assay": (
                "Subject 10^8 differentiated Stage-7 cells to 72-h 1% O2 hypoxia + IFN-g/TNF-a/IL-1b + sub-lethal AP1903 cycling followed by 28-day outgrowth, "
                "requiring zero positive selection or clonal enrichment of 20q11.21/BCL2L1 copy number or TP53/PPM1D/BCOR variants "
                "(obs_post_stress_20q11_bcl2l1_outgrowth_fold <= 1.0 point estimate and upper95 <= 1.20x) and functional p53 Nutlin-3a response >= 5.0x (obs_nutlin3a_p21_mdm2_fold_induction >= 5.0)."
            ),
            "wetlab_contract_endpoint": "obs_post_stress_20q11_bcl2l1_outgrowth_fold",
            "functional_p53_nutlin3a_endpoint": "obs_nutlin3a_p21_mdm2_fold_induction",
            "min_nutlin3a_p21_mdm2_fold_induction": 5.0,
        },
        "wt_iapp_pouch_accumulation_and_pramlintide_lead_promotion_rule_muse_r6": {
            "finding_ids": ["MW-R6-S1-14"],
            "decision_id": "DEC-053",
            "primary_citations": ["PMID:18753660", "PMID:20835230"],
            "pouch_iapp_mass_balance": (
                "40 IU/day insulin (~240 nmol/day) co-secretes 2.4-12.0 nmol/day human WT IAPP (at 1-5% molar ratio); distributed into a 7.0 mL pouch volume "
                "with restricted macromolecular washout and matrix binding, local interstitial WT IAPP generation reaches 0.34-1.71 uM/day, crossing the micromolar "
                "nucleation threshold for toxic beta-sheet oligomerization and NLRP3 inflammasome IL-1beta activation (Westermark et al. 2008 PMID:18753660; Masters et al. 2010 PMID:20835230)."
            ),
            "conditional_lead_promotion_rule": (
                "Stage-0 Gate 0E-amyloid runs a 7-to-14-day confined high-glucose (16.7 mM) + exogenous preformed WT human IAPP fibril-seeding challenge in the 350-um hydrogel coupon; "
                "if WT IAPP exhibits Thioflavin-T (ThT) amyloid area > 1.0% (> 2.0x baseline), A11 toxic oligomer positivity, macrophage NLRP3 IL-1beta release (> 50 pg/mL), or > 3-5% beta-cell "
                "caspase-3/apoptosis, non-amyloidogenic pramlintide-like IAPP(A25P/S28P/S29P) (product candidate ISLET1_A4_IAPP_PRAM / satellite coupon ER_IAPP_PRAMLINTIDE_RNLS_WT, K=0 DSBs, B=1 base edit) is promoted "
                "via sequential DSB-free base editing (or paired only with DSB-only K<=4 stacks so simultaneous base-editor burden remains B<=1 and avoids the GRAVEYARD_006 B>=2 WARN engineering-complexity discount rule; DEC-054/057/058)."
            ),
        },
    }
    hip_core_protocol["layer1_ultracode_audit_reconciliation_ledger"] = layer1_audit_ledger
    with open(protocols_dir / "hip_core_fractional_factorial_v0.6.2.json", "w", encoding="utf-8") as f:
        json.dump(hip_core_protocol, f, indent=2, ensure_ascii=False)
        f.write("\n")

    out_doc = {
        "schema_version": "0.7.0",
        "project": "OpenIslet.org (Project Islet-1) — Sign-Sensitive, Endpoint-Oriented & Locus-Penalized Evaluator",
        "champion_challenger_portfolio": champion_challenger_portfolio,
        "layer1_ultracode_audit_reconciliation_ledger": layer1_audit_ledger,
        "wetlab_bayesian_calibration_integration_status": wetlab_calibration_link_status,
        "estimand_disclaimer": (
            "All scores in this file represent a Relative Multi-Gate Evidence & Feasibility Index in [0, 1] "
            "combining sign-preserved comparative effect sizes, counter-evidence penalties, human genetic "
            "constraints (gnomAD pLI/LOEUF, OpenTargets), explicit engineering priors (with sensitivity ranges), "
            "and multiplex chromosomal/nickase edit penalties. They are NOT calibrated human clinical 1-year "
            "graft survival probabilities, and heuristic_evidence_sparsity_se / heuristic_conservative_floor_index "
            "are deterministic scenario penalties rather than calibrated Bayesian posterior intervals."
        ),
        "based_on_validated_ring1_studies": len(studies),
        "total_atomic_effects_synthesized": graph.get("total_atomic_effects_extracted", 0),
        "total_effect_modifiers_synthesized": graph.get("total_effect_modifier_interactions_extracted", 0),
        "total_multi_hop_causal_chains_synthesized": len(graph.get("multi_hop_cross_study_causal_chains", [])),
        "total_valid_combinatorial_stacks_evaluated": len(all_combos),
        "total_feasible_unblocked_stacks": sum(1 for c in all_combos if c.get("is_feasible_candidate", True)),
        "total_graveyard_blocked_stacks": sum(1 for c in all_combos if not c.get("is_feasible_candidate", True)),
        "engineering_priors_and_sensitivity_table": ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE,
        "complete_gate_coefficient_registry": COMPLETE_GATE_COEFFICIENT_REGISTRY,
        "direct_output_logit_sensitivity_and_variance_decomposition": direct_logit_sensitivity,
        "whole_product_clonogenic_safety_model": whole_product_safety_model,
        "unconditional_patient_level_tumor_harm_and_two_stage_psc_purge_model_sonnet55_r5": whole_product_safety_model["unconditional_patient_level_tumor_harm_and_two_stage_psc_purge_model_sonnet55_r5"],
        "non_hla_neoepitope_and_complement_flow_crossmatch_audit_sonnet55_r5": {
            "finding_ids": ["AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "MW-R5-S2-F05", "AUDIT-REF"],
            "provenance": "Audit Convergence Round 5 Step 5 & Audit Round 5 Step 6 (DEC-047/048/049; Ravindranath 2022 PMID:35323192; Deuse 2021 PMID:33416832; Straathof 2005 PMID:15728125; Moran 1991 PMID:1824699)",
            "chimeric_surface_junction_neoepitope_inventory_count": 306,
            "exact_boundary_constrained_8_to_11mer_window_count": 326,
            "exact_boundary_constrained_8_to_11mer_window_count_with_hlac_l2_alpha1": 360,
            "sliding_8_to_11mer_window_formula": "Nominal single-cut upper-bound formula: 9 chimeric junction entries * (7 + 8 + 9 + 10 = 34) junction-spanning 8..11-mer windows = 306 windows; exact boundary-length-constrained window count is 326 across the 9 entries (33 each for the three 9-mer tethered peptides [VMAPRTLVL, GAVDPLLAL, QYDDAVYKL] and the 9-aa retained CD55 stalk [PNKGSGTTS] where len(flank)=9 caps 11-mers at 9, 34 each for the four >=10-aa-flanked single-cut junctions [L1|B2M, B2M|L2, L2|HLA-E_alpha1, SERPINB9|GSG-P2A], and 58 across the two-boundary 6-aa linker in FKBP12-SGGGSG-dCasp9; or 360 when including the 34-window L2|HLA-C_alpha1 junction GGGGSGGGGSGGGGSGGGGS-GSHSMRYFSTSVSRPGRGEP)",
            "surface_ectodomain_junctions_audited": [
                "VMAPRTLVL-GGGGSGGGGSGGGGSGGGGS (HLA-E SCT 20-aa L1 linker in 4,143-bp Embodiment A*; or VMAPRTLVL-GGGGSGGGGSGGGGS 15-aa L1 in 3,974-bp canonical donor)",
                "GGGGSGGGGSGGGGSGGGGS-IQRTPKIQVY (L1 C-terminus to mature B2M Ile1-Tyr10)",
                "VNHVTLSQPKIVKWDRDM-GGGGSGGGGSGGGGSGGGGS (mature B2M C-terminus to 20-aa L2 linker)",
                "GGGGSGGGGSGGGGSGGGGS-GSHSLKYFHTSVSRPGRGEP (L2 C-terminus to mature HLA-E alpha1 Gly1-Pro20)",
                "ELKYRVVSWFS-PNKGSGTTS (CD47-IgV Phe137-Ser138 to 9-aa retained CD55/DAF stalk Pro139-Ser147 omega-site [native CD55 Pro345-Ser353]; 28-aa hydrophobic propeptide GTTRLLSGHTCFTLTGLLGTLVTMGLLT 148..175 is cleaved at Ser147 upon GPI transamidation)",
                "GAVDPLLAL-GGGGSGGGGSGGGGSGGGGS (Track 2 B'-SCT C1 HLA-C*03:04[Y84A] 20-aa L1 linker; plus L2-to-HLA-C alpha1 junction GGGGSGGGGSGGGGSGGGGS-GSHSMRYFSTSVSRPGRGEP)",
                "QYDDAVYKL-GGGGSGGGGSGGGGSGGGGS (Track 2 B'-SCT C2 HLA-C*04:01[Y84A] 20-aa L1 linker; plus L2-to-HLA-C alpha1 junction GGGGSGGGGSGGGGSGGGGS-GSHSMRYFSTSVSWPGRGEP)",
            ],
            "intracellular_chimeric_junctions_audited_for_indirect_apc_class_ii_presentation": [
                "NSILFCGRFSSP-GSGATNFSLLKQAGDVEENPG (Locus 3 canonical R01 SERPINB9 C-terminus [Asn365-Pro376] to 21-aa retained GSG-P2A stub after Gly20|Pro21 ribosomal skipping [14 junction-spanning 15-mer windows]; eliminated in scar-free L3-BF-B and 3,799-bp comparator)",
                "ATLVFDVELLKLE-SGGGSG-GFGDVGALESLRGNAD (Locus 4 AAVS1 iCasp9 FKBP12[F36V] C-terminus [Ala96-Glu108] + 6-aa SGGGSG linker [109..114] + dCARD-CASP9 [Gly135-Asp150] [18 junction-spanning 15-mer windows; 32 intracellular 15-mer windows total across both intracellular junctions, and 136 boundary-constrained 15-mer windows across all 10 ectodomain + intracellular junctions])",
            ],
            "mandatory_stage0_gate_0d_0g_flow_crossmatch_and_cdc_protocol": {
                "recipient_serum_panel": "n >= 50 human sera (25 high-cPRA >= 80% sensitized allogeneic sera + 25 T1D sera, including anti-HLA-E*01:03 Luminex SAB screening)",
                "assays": [
                    "Flow cytometric IgG and IgM crossmatch on Stage-7 Embodiment A* and B'-SCT clusters (pass: MFI shift <= 1.5x isotype)",
                    "20% normal human serum (NHS) C3b/C4d/C5b-9 MAC deposition and CDC lysis with vs without SA-THBD + SA-FHL-1 microgels (pass: C5b-9 MAC+ <= 5.0%, CDC lysis <= 5.0%)",
                    "Opsonized macrophage ADCP with recipient serum + primary human macrophages (pass: phagocytosis index <= 5.0%)",
                ],
            },
        },
        "ranker_graveyard006_and_unswitched_stealth_audit_sonnet55_r5": {
            "finding_id": "AUDIT-REF",
            "stacks_with_multiplex_base_edit_B_ge_2_and_cuts_K_ge_4": sum(
                1 for c in all_combos if c.get("multiplex_base_edit_plus_dsb_translocation_risk_flag")
            ),
            "stacks_with_unswitched_immune_evasive_clinical_hold": sum(
                1 for c in all_combos if c.get("unswitched_immune_evasive_clinical_hold")
            ),
        },
        "pareto_frontier_by_locus_budget_K": [pareto_by_k[k] for k in sorted(pareto_by_k.keys())],
        "all_evaluated_combinatorial_stacks": all_combos,
        "top_15_combinatorial_search_results": all_combos[:15],
        "bottom_5_combinatorial_traps": all_combos[-5:],
        "ranked_candidate_architectures": ranked_architectures,
        "target_gene_evidence_and_constraint_profiles": target_profiles,
        "hip_core_fractional_factorial_protocol": hip_core_protocol,
        "odds_composition_and_profile_scope_audit_iip_0003": _build_odds_composition_and_profile_scope_audit_iip_0003(),
        "step7_high_value_of_information_gaps": [
            {
                "voi_id": "DECISION_EXP_1_HOUR0_B2M_PROMOTER_SURFACE_DENSITY",
                "title": "Resting Hour-0 QuantiBRITE Surface Density & Head-to-Head Construct Comparison (Bicistronic HLA-E-GSG-tPT2A-CD47-IgV-GPI vs 2A-Free Bidirectional Cassette vs CD47-dCT vs CD47-WT; HLA-E*01:03[VMAPRTLVL] WT vs KW-Motif) Across Stage 7 Beta, Alpha, Delta & iEC Lineages",
                "why_needed": "Mature resting human Stage 7 beta-cells have low basal B2M transcription before IFN-gamma exposure, and Xu et al. 2024 (PMID:38508139) tested GPI-anchored CD47-IgV-GPI in HUVECs/teratomas whereas Ghimire et al. 2023 (PMID:37820008) showed full-length CD47-WT inhibits beta-cell GSIS via Cdc42 and Lin et al. 2023 (PMID:37264229) showed VMAPRTLVL occurs in ~12% of clinical HCMV UL40 isolates. Because Golgi Furin cannot cleave the cytosolic C-terminal -GSG-P2A tail on HLA-E*01:03 in the bicistronic tPT2A cassette (Audit R01), Stage-0 compares bicistronic HLA-E-GSG-tPT2A-CD47-IgV-GPI against a scarless 2A-free bidirectional cassette [SV40pA-CD47-IgV-GPI(antisense) <- EEF1A1-core | Syn_5xISRE_minP -> HLA-E*01:03(VMAPRTLVL)-bGHpA] as well as 5'-CD47-dCT/WT-tPT2A-HLA-E for Hour-0 QuantiBRITE density (>200k CD47, >15k HLA-E), anti-2A read-through Western blot, dynamic GSIS, iEC tube formation, and protection against FACS-gated NKG2A- KIR2DL1/2/3+ and HCMV+ NKG2C+ adaptive NK subsets.",
                "kill_criterion": "If CD47-IgV-GPI fails to match CD47-WT NK/macrophage protection across SIRPA V1/V2 donors, or if bicistronic tPT2A shows >5% unskipped read-through or impaired HLA-E trafficking relative to the 2A-free bidirectional cassette, lock the highest-performing construct uniformly across all 13 X2=+1 factorial runs prior to Stage 1.",
            },
            {
                "voi_id": "DECISION_EXP_2_INDIRECT_APC_PRESENTATION_AND_BYSTANDER_KILLING",
                "title": "Host HLA-DQ2/DQ8 Dendritic Cell Indirect Presentation, Macrophage Bystander Killing & Cytokine-Stimulated Beta-Cell Fas (CD95) Incremental Fratricide Screen Under SA-FasL / SA-THBD Microgels",
                "why_needed": "Deleting donor MHC-I/II eliminates direct presentation but leaves indirect presentation of shed autoantigens intact (Audit §3.4). Moreover, published murine and NHP SA-FasL PEG-microgel studies (Headen et al. 2018 PMID:29867165; Lei et al. 2022 PMID:35559682) used transient short-course rapamycin (whereas Turan et al. 2024 PMID:38709159 evaluated direct islet surface biotinylation co-displaying SA-THBD [SA-TM] and SA-CD47 in a murine syngeneic intraportal IBMIR/engraftment model without rapamycin, and Kazmi et al. 2026 PMID:41842939 evaluated direct surface SA-THBD + SA-CD47 in murine BALB/c->C57BL/6 intraportal allografts without SA-FasL, PEG microgels, or systemic immunosuppression), and human beta-cells upregulate Fas (CD95) apoptosis sensitivity under IL-1b (Loweth et al. 1998 PMID:9588443), creating a potential SA-FasL beta-cell fratricide hazard (Audit F6; Audit R06; AUDIT-REF & L3B-01/02/03 [DEC-076/077]).",
                "kill_criterion": "Because IL-1b + IFN-g pre-treatment alone elevates baseline beta-cell apoptosis above 5%, the kill gate is defined on the incremental SA-FasL treatment contrast (Audit R06): if SA-FasL + SA-THBD microgels increase cleaved-Caspase-8/TUNEL+ beta-cell apoptosis by Delta_TUNEL+ > +5.0 percentage points above the matched cytokine-treated SA-THBD-only vehicle control (or if blocking anti-Fas ZB4 antibody rescues >5.0 pp of apoptosis in the SA-FasL arm), drop SA-FasL from the hydrogel and retain only SA-THBD.",
            },
            {
                "voi_id": "DECISION_EXP_3_CXCL10_INDUCTION_KINETICS_RACE",
                "title": "Head-to-Head Topology & Kinetic Race: CXCL10::SERPINB9-GSG-tPT2A-PD-L1 vs CXCL10::PD-L1-GSG-tPT2A-SERPINB9 vs First-Contact Perforin/Granzyme B Killing",
                "why_needed": "Compares the inverted cistron topology (SERPINB9-GSG-tPT2A-PD-L1, which leaves PD-L1's C-terminal Cys272 palmitoylation/CMTM6 tail untagged [PMID:30952982, PMID:28813417] while appending the 21-aa P2A tail to cytosolic SERPINB9) against PD-L1-GSG-tPT2A-SERPINB9 for cytosolic PI-9 half-life, Granzyme B k_ass, surface PD-L1 recycling, and Hour-0 perforin/Granzyme B race kinetics (Audit §3.3; Audit F4).",
                "kill_criterion": "Select the cistron orientation with superior cytosolic Granzyme B inhibition and surface PD-L1 stability; if Hour-0 CTL/NK killing exceeds 15% before CXCL10 induction, include the minimal basal hPGK core floor.",
            },
            {
                "voi_id": "DECISION_EXP_4_BIALLELIC_AAVS1_ICASP9_ESCAPE_FREQUENCY_BOUND",
                "title": "Whole-Product Clonogenic Escape Target (<10^-8 per Residual Proliferating Cell) Under Hypoxia/Cytokine XIAP Stress & >=63,000x Target-Enriched Duplex Sequencing Gate (>=62,958 Duplex Molecules for >=3 Mutant Molecules at 10^-4 VAF)",
                "why_needed": "Under N_prolif <= 8x10^4 effective clonogenic opportunities per 8x10^8-cell dose, pooled MCB validation (>=3.0x10^8 cells; PMID:35795779) must prove p_escape <= 10^-8 under 1% O2 + IFN-g/TNF-a/IL-1b (evaluating whether XIAP BIR3 binding to cleaved Caspase-9 [PMID:11242052] causes resistance and testing the preclinical AP1903 + SMAC mimetic salvage hypothesis [Audit F8]), while longitudinal TP53/PPM1D/BCOR duplex sequencing must sample >=62,958 independent duplex molecules (~63,000x target-enriched duplex depth) to achieve >=95% Poisson power to observe >=3 mutant duplex molecules at VAF = 10^-4 (vs >=29,956 molecules [~30,000x] for >=1 mutant molecule and 87.5% power for >=3 molecules at 50,000x; Finding 077).",
                "kill_criterion": "No in vivo pediatric translation until residual proliferating fraction is <=10^-4, upper 95% confidence bound on per-clonogenic-cell escape is <=10^-8, and >=63,000x target-enriched duplex sequencing (>=62,958 independent duplex molecules) confirms zero TP53/PPM1D/BCOR driver subclones >=10^-4 VAF with >=3 mutant molecules.",
            },
            {
                "voi_id": "DECISION_EXP_5_FUNCTION_UNDER_PROTECTION_SHUTOFF_PERIFUSION",
                "title": "Head-to-Head Dynamic Perifusion (16.7 mM -> 2.8 mM Low-Glucose Shut-Off Time & Seahorse OCR) Across Isogenic WT, RNLS-/-, G6PC2-/-, and SLC30A8 Allelic Series",
                "why_needed": "Tests RNLS-/- under an OCR precaution; MacDonald 2025 PMID:39988068 Fig. 6 showed poor GSIS in both WT and KO and does not establish a KO-specific GSIS penalty and G6PC2-/- fasting hypoglycemia risk (Bosma 2020 PMID:32569842; Audit §7.5).",
                "kill_criterion": "Any survival edit that delays 2.8 mM glucose insulin shut-off beyond 3 minutes or increases basal insulin leak is permanently excluded.",
            },
            {
                "voi_id": "OQ_016_FOUR_HURDLE_ODDS_ASSUMPTIONS_AND_SUCCESS_ENDPOINT_IIP_0003",
                "open_question_id": "OQ-016",
                "title": "Which assumptions and definition of success sit behind the four-hurdle odds estimate?",
                "why_needed": "Raised by Linus Schumacher (ORCID 0000-0003-0797-3406, GitHub @ljschumacher, University of Edinburgh) in IIP-0003 (https://github.com/openislet/openislet/issues/3). The 21% and 37% figures are rough team judgements for clearing four in-body hurdles in sequence, assuming the edited cells work as designed and pass lab safety checks before any trial. Moving from one design to three depends on how much Designs A, B and X overlap at each hurdle, and the estimate depends on which clinical outcome counts as success.",
                "kill_criterion": "Fix one clinical definition of success, write down the per-design and shared failure assumptions at each hurdle, and replace the starting judgements with measured pass rates as lab and preclinical results arrive.",
            },
        ],
    }

    if not (PROJECT_ROOT / ".openislet_public_bundle").is_file():
        with open(OUT_PATH, "w", encoding="utf-8") as f:
            json.dump(out_doc, f, indent=2, ensure_ascii=False)
            f.write("\n")
    for alias_name in ("islet1_stack_ranking_v0.5.json",):
        alias_path = PROJECT_ROOT / "data" / "extracted" / alias_name
        with open(alias_path, "w", encoding="utf-8") as f:
            json.dump(out_doc, f, indent=2, ensure_ascii=False)
            f.write("\n")

    print(
        f"[ISLET-1 v0.7.0 EVALUATOR COMPLETE] Evaluated {len(all_combos)} combinations across {len(studies)} studies:"
    )
    for a in ranked_architectures:
        ped_u = a["utility_pediatric_early_onset_t1d"]
        lada_u = a["utility_adult_onset_lada"]
        print(
            f"  Named #{a['rank_among_6_named_architectures']} (Global #{a['global_rank_among_all_combinations']}/{len(all_combos)}, feasible={a['is_feasible_candidate']}): "
            f"{a['architecture_id']} (K={a['genomic_locus_edit_count_K']}, B={a['base_edit_count_B']}) | "
            f"Ped Index={ped_u['series_system_joint_survival_probability']:.4f} | "
            f"LADA Index={lada_u['series_system_joint_survival_probability']:.4f} | "
            f"Composite={a['composite_series_system_utility']:.4f} | "
            f"Mfg={ped_u['manufacturing_locus_penalty_factor']:.3f} | "
            f"Limiting Gate={ped_u['limiting_bottleneck_gate']} ({ped_u['limiting_bottleneck_score']:.2f})"
        )


if __name__ == "__main__":
    main()
