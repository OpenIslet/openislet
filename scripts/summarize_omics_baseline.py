#!/usr/bin/env python3
"""
Compiles the live OpenTargets T1D genetic associations (MONDO_0005147) and gnomAD
loss-of-function constraint scores into data/omics/target_prior_table.json for Project Islet-1 (OpenIslet.org).
Covers 33 primary candidate targets across all 5 Gates and all 12 Failure Modes, with explicit
provenance labeling distinguishing the 3 direct OpenTargets MONDO_0005147 top-25 GraphQL hits (INS, PTPN2, TYK2)
from the 17 curated literature prior fallbacks and 13 unindexed/null structural targets (AUDIT-REF).
"""

from __future__ import annotations

import glob
import json
import math
import re
import sys
from pathlib import Path

try:
    from scripts.study_artifact_exports import publish_byte_artifacts
except ModuleNotFoundError as exc:
    if exc.name not in {"scripts", "scripts.study_artifact_exports"}:
        raise
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from study_artifact_exports import publish_byte_artifacts

BASE_DIR = Path(__file__).resolve().parent.parent
OMICS_DIR = BASE_DIR / "data" / "omics"

ENGINEERING_ANNOTATIONS = {
    "B2M": {
        "opentargets_t1d_score": None,
        "gate": "Gate2_Allo_HLA_ClassI_or_II",
        "failure_mode_defended": "FM2_NK_MissingSelf_and_80yr_Viral_Cancer_BlindSpot",
        "recommended_edit": "Endogenous_Exon1_SwapIn_B2M_HLA_E_VMAPRTLVL_P2A_CD47_with_Hour0_MFI_Gate",
        "hpa_islet_expression_status": "Expressed_Basal_Low_IFN_Inducible_High_Descriptive_Context_Only",
        "rationale": "Total MHC-I light chain (pLI=0.948). Knocking HLA-E*01:03 single-chain trimer (VMAPRTLVL leader) + CD47 into B2M exon 1 eliminates HLA-A/B/C in 1 locus, but requires calibrated Hour-0 surface molecule/cell verification in resting Stage 7 beta-cells and iECs (which have low basal B2M transcription) and ablates non-classical B2M partners (FcRn, HFE, CD1d, MR1).",
    },
    "CASP9": {
        "opentargets_t1d_score": 0.0582,
        "gate": "Gate5_80yr_Lineage_Independent_Safety_Switch",
        "failure_mode_defended": "FM2_Oncogenic_Teratoma_or_Viral_Escape_Erasability",
        "recommended_edit": "Biallelic_AAVS1_CAG_or_GAPDH_P2A_iCasp9_KnockIn",
        "hpa_islet_expression_status": "Ubiquitous_Apoptotic_Effector_Descriptive_Context_Only",
        "rationale": "Inducible Caspase-9 (AP1903/rimiducid dimerizer). Must be knocked in biallelically at a pan-lineage housekeeping/safe-harbor locus (AAVS1-CAG or GAPDH-P2A; Wunderlich 2022 PMID:35795779) with post-differentiation anti-proliferative purification so whole-graft clonogenic escape expectation remains <<1 event per 8x10^8-cell patient dose.",
    },
    "CD274": {
        "opentargets_t1d_score": 0.1143,
        "gate": "Gate2_and_Gate4_Inducible_Checkpoint_Shield",
        "failure_mode_defended": "FM5_Bystander_and_Autoreactive_T_Cell_Checkpoints",
        "recommended_edit": "Endogenous_CXCL10_Promoter_SwapIn_SERPINB9_P2A_PDL1",
        "hpa_islet_expression_status": "Inducible_Islet_Beta_Cell_IFN_Response_Descriptive_Context_Only",
        "rationale": "PD-L1 checkpoint upregulated on surviving beta-cells in Joslin 50-Year T1D Medalists (Colli 2018 PMID:30269996) and HILOs (Yoshihara 2020 PMID:32814902). Requires kinetic verification of induction lag vs first-contact CTL killing.",
    },
    "CD47": {
        "opentargets_t1d_score": 0.0468,
        "gate": "Gate3_NK_and_Macrophage_Phagocytosis_Shield",
        "failure_mode_defended": "FM2_NK_and_Macrophage_MissingSelf_Clearance",
        "recommended_edit": "High_Density_CoExpression_Above_SIRPa_Threshold",
        "hpa_islet_expression_status": "Ubiquitous_Membrane_Protein_Requires_SupraPhysiological_Overexpression_Descriptive_Context_Only",
        "rationale": "Binds inhibitory SIRPa on macrophages and NK cells (Deuse 2019 PMID:30778232; Hu 2024 PMID:37156915). Requires verification that resting Hour-0 surface density on Stage 7 beta-cells and iECs exceeds the donor SIRPa polymorphism threshold.",
    },
    "CD55": {
        "opentargets_t1d_score": None,
        "gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "failure_mode_defended": "FM10_IBMIR_ThromboInflammation_and_Complement_Lysis",
        "recommended_edit": "CoExpression_Complement_Decay_Accelerating_Factor",
        "hpa_islet_expression_status": "Endothelial_and_Islet_Surface_Descriptive_Context_Only",
        "rationale": "Complement decay-accelerating factor (DAF); blocks C3 convertase deposition from autoantibodies/IBMIR (pLI=3.9e-13).",
    },
    "CD59": {
        "opentargets_t1d_score": None,
        "gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "failure_mode_defended": "FM10_IBMIR_ThromboInflammation_and_Complement_Lysis",
        "recommended_edit": "Endogenous_Preservation_or_Boost",
        "hpa_islet_expression_status": "High_Endogenous_Beta_Cell_Surface_Descriptive_Context_Only",
        "rationale": "Protectin; blocks Membrane Attack Complex (C5b-9 MAC) pore formation on beta-cells.",
    },
    "CHGA": {
        "opentargets_t1d_score": None,
        "gate": "Gate2_Autoantigen_Stealth_TCR_or_AutoAb",
        "failure_mode_defended": "FM4_PostTranslational_Neoantigens_HIPs_DRiPs_Deamidation",
        "recommended_edit": "WE14_Cleavage_Site_BaseEdit_or_KO",
        "hpa_islet_expression_status": "High_Pan_Endocrine_Secretory_Granule_Descriptive_Context_Only",
        "rationale": "Chromogranin A (WE14 peptide); primary donor fragment for C-peptide-WE14 Hybrid Insulin Peptides (HIPs). LoF tolerant in gnomAD (pLI=4.3e-11, oe_lof=0.754).",
    },
    "CIITA": {
        "opentargets_t1d_score": None,
        "gate": "Gate2_Allo_HLA_ClassI_or_II",
        "failure_mode_defended": "FM1_Alloreactive_T_Cell_and_DSA_Rejection",
        "recommended_edit": "Biallelic_KO",
        "hpa_islet_expression_status": "Cytokine_Inducible_Transactivator_Descriptive_Context_Only",
        "rationale": "Master MHC Class II transactivator (pLI=2.0e-14); biallelic KO is intended to suppress donor HLA-DR/DP/DQ-mediated direct presentation; host APC indirect presentation and B-cell help remain possible and require testing.",
    },
    "CXADR": {
        "opentargets_t1d_score": None,
        "gate": "Gate4_Intrinsic_IFN_JAK_STAT_Block",
        "failure_mode_defended": "FM8_Enterovirus_Reinfection_and_MDA5_Interferon_Storm",
        "recommended_edit": "Beta_Cell_Granule_Isoform_Splicing_BaseEdit",
        "hpa_islet_expression_status": "Islet_Beta_Cell_Secretory_Granule_Isoform_Descriptive_Context_Only",
        "rationale": "Coxsackievirus and Adenovirus Receptor (CAR); beta-cell specific C-terminal PDZ isoform localizes to insulin granules and mediates Coxsackie B viral entry during exocytosis (pLI=0.0061, LOEUF=0.722).",
    },
    "CXCL10": {
        "opentargets_t1d_score": 0.1197,
        "gate": "Gate4_Chemokine_Silence_and_Inducible_Landing_Pad",
        "failure_mode_defended": "FM5_CXCR3_T_Cell_Recruitment_Loop",
        "recommended_edit": "Biallelic_Coding_Replacement_CXCL10_SERPINB9_P2A_PDL1",
        "hpa_islet_expression_status": "Absent_Basal_High_IFN_Inducible_Islet_Chemokine_Descriptive_Context_Only",
        "rationale": "Primary beta-cell CXCR3 chemokine induced by IFN-alpha/gamma and enterovirus (Tanaka 2009 PMID:19641142). Replacing the targeted CXCL10 coding region with SERPINB9-P2A-PD-L1 aims to eliminate donor CXCL10 secretion; exact edited-allele and secretion assays are pending (note: host endothelial/macrophage CXCL9/10/11 remains unedited).",
    },
    "ENTPD1": {
        "opentargets_t1d_score": None,
        "gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "failure_mode_defended": "FM10_IBMIR_ThromboInflammation_and_Complement_Lysis",
        "recommended_edit": "Inducible_Expression_CD39",
        "hpa_islet_expression_status": "Endothelial_and_Vascular_Enriched_Descriptive_Context_Only",
        "rationale": "CD39 ectonucleotidase; hydrolyzes extracellular ATP/ADP to AMP; downstream CD73 activity is required to produce adenosine. Exact-product thrombo-inflammatory protection remains unmeasured.",
    },
    "F3": {
        "opentargets_t1d_score": None,
        "gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "failure_mode_defended": "FM10_IBMIR_ThromboInflammation_and_Complement_Lysis",
        "recommended_edit": "Biallelic_KO_or_Promoter_Attenuation",
        "hpa_islet_expression_status": "Expressed_on_Isolated_and_Stressed_Islets_Descriptive_Context_Only",
        "rationale": "Tissue Factor (CD142); expressed on human islets and triggers Instant Blood-Mediated Inflammatory Reaction (IBMIR) clotting cascade upon blood contact. LoF tolerant in gnomAD (pLI=0.0024, LOEUF=0.796).",
    },
    "G6PC2": {
        "opentargets_t1d_score": 0.245,
        "gate": "Gate2_Autoantigen_Stealth_TCR_or_AutoAb",
        "failure_mode_defended": "FM3_Recurrent_Autoantibodies_Exocytosis_ADCC_and_Epitope_Spreading",
        "recommended_edit": "Defer_KO_Until_2_8mM_Shutoff_Validated_NEVER_Use_For_iCasp9",
        "hpa_islet_expression_status": "Islet_Beta_Cell_Enriched_Mature_Lineage_Restricted_Descriptive_Context_Only",
        "rationale": "IGRP CD8+ autoantigen (pLI=2.4e-16). Removed from Universal Chassis A core in v0.5/v0.6: G6pc2-/- lowers fasting glucose and risks fasting/exercise hypoglycemia (Bosma 2020 PMID:32569842; Overway 2022 PMID:34954144), and lineage-restricted G6PC2::iCasp9 fails to ablate dedifferentiated, teratoma, or SLC18A1+ off-target cells (Yagyu 2016 PMID:27626039; Wunderlich 2022 PMID:35795779).",
    },
    "GAD2": {
        "opentargets_t1d_score": 0.375,
        "gate": "Gate2_Autoantigen_Stealth_TCR_or_AutoAb",
        "failure_mode_defended": "FM3_Recurrent_Autoantibodies_Exocytosis_ADCC_and_Epitope_Spreading",
        "recommended_edit": "BaseEdit_ABE_CBE_TCR_and_Bcell_Epitopes_NO_GROSS_KO",
        "hpa_islet_expression_status": "Islet_Beta_Cell_and_CNS_Enriched_Descriptive_Context_Only",
        "rationale": "GAD65 autoantigen (GADA, dominant in LADA); pLI=1.000; proposed epitope edits must separately demonstrate retained intra-islet GABA function and reduced antibody/TCR recognition; neither follows from locus constraint alone.",
    },
    "GLIS3": {
        "opentargets_t1d_score": 0.412,
        "gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "failure_mode_defended": "FM7_ER_Stress_Unfolded_Proinsulin_and_CTL_Bystander_Death",
        "recommended_edit": "Preserve_Expression_and_Check_Family_Genotype",
        "hpa_islet_expression_status": "Islet_Beta_Cell_Transcription_Factor_Descriptive_Context_Only",
        "rationale": "T1D/neonatal diabetes transcription factor preventing BIM-induced beta-cell apoptosis.",
    },
    "HLA-C": {
        "opentargets_t1d_score": 0.0704,
        "gate": "Gate2_Allo_HLA_ClassI_and_Gate3_NK_KIR",
        "failure_mode_defended": "FM1_and_FM2_HLA_C_Superdonor_Retention",
        "recommended_edit": "Retain_Matched_Allele_Only_in_Chassis_B_Hypothesis",
        "hpa_islet_expression_status": "Basal_Low_IFN_Inducible_High_Descriptive_Context_Only",
        "rationale": "Retaining 1 matched HLA-C allele in HLA-A/B/CIITA-/- iPSCs inhibits KIR2DL+ NK cells (Xu 2019 PMID:30853558), but requires testing whether retained HLA-C signal peptide suffices to stabilize endogenous HLA-E against NKG2A+ NK cells (Lin 2023 PMID:37264229) and presents beta-cell autoantigens.",
    },
    "HLA-E": {
        "opentargets_t1d_score": 0.0275,
        "gate": "Gate3_NK_MissingSelf_Shield",
        "failure_mode_defended": "FM2_NK_MissingSelf_Lysis",
        "recommended_edit": "SingleChain_Trimer_HLA_E_0103_VMAPRTLVL_KnockIn",
        "hpa_islet_expression_status": "Basal_Low_Requires_SingleChain_Trimer_Rescue_upon_B2M_Disruption_Descriptive_Context_Only",
        "rationale": "Non-classical MHC-I binding inhibitory CD94/NKG2A on NK cells (Gornalusse 2017 PMID:28504668). The proposed VMAPRTLVL peptide also occurs in some CMV UL40 isolates; NKG2A/NKG2C receptor discrimination and HCMV-positive donor protection require direct assays (PMID:37264229).",
    },
    "IAPP": {
        "opentargets_t1d_score": None,
        "gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "failure_mode_defended": "FM9_Human_IAPP_Amyloid_Oligomer_Toxicity_and_HIP_Neoantigens",
        "recommended_edit": "Pramlintide_Proline_Triplet_PrimeEdit_A25P_S28P_S29P_or_KO",
        "hpa_islet_expression_status": "High_Islet_Beta_Cell_CoSecreted_Hormone_Descriptive_Context_Only",
        "rationale": "Human Islet Amyloid Polypeptide (Amylin) misfolds into toxic amyloid plaques in transplanted human islets (NLRP3 inflammasome trigger + Ins-IAPP Hybrid Insulin Peptide source). Population LoF constraint is not evidence that IAPP loss is functionally harmless in an engineered graft. The proposed A25P/S28P/S29P substitutions are motivated by pramlintide sequence; complete prevention of amyloid or hybrid-insulin-peptide formation and preserved processing/secretion are unmeasured exact-product hypotheses.",
    },
    "IFIH1": {
        "opentargets_t1d_score": 0.485,
        "gate": "Gate4_Intrinsic_IFN_JAK_STAT_Block",
        "failure_mode_defended": "FM8_Enterovirus_Reinfection_and_MDA5_Interferon_Storm",
        "recommended_edit": "Biallelic_KO_or_Protective_I923V_BaseEdit",
        "hpa_islet_expression_status": "IFN_Inducible_Cytosolic_Sensor_Descriptive_Context_Only",
        "rationale": "Cytosolic viral sensor MDA5; natural human LoF variants reduce T1D risk by ~50% (pLI=6.0e-49, oe_lof=1.164).",
    },
    "INS": {
        "opentargets_t1d_score": 0.7333,
        "gate": "Gate1_Differentiation_and_Gate2_Autoantigen",
        "failure_mode_defended": "FM3_and_FM4_PPI_Leader_and_INS_DRiP_Autoantigen",
        "recommended_edit": "Preserve_Mature_Insulin_Optional_Leader_or_Silent_DRiP_AUG_Edit_in_HLA_Retained_Arms",
        "hpa_islet_expression_status": "Extreme_Islet_Beta_Cell_Core_Hormone_Descriptive_Context_Only",
        "rationale": "#1 OpenTargets T1D locus (score=0.733). Mature insulin is essential; signal-peptide swap (Kobaisi 2026 PMID:42690930) or silent out-of-frame AUG removal (Kracht 2017 PMID:28263308) is only needed in HLA-I-retained arms or if indirect presentation of DRiP peptides is significant.",
    },
    "MANF": {
        "opentargets_t1d_score": None,
        "gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "failure_mode_defended": "FM7_ER_Stress_Unfolded_Proinsulin_and_CTL_Bystander_Death",
        "recommended_edit": "Preserve_or_Secreted_Trophic_Support",
        "hpa_islet_expression_status": "ER_Secretory_Trophic_Factor_Descriptive_Context_Only",
        "rationale": "Mesencephalic astrocyte-derived neurotrophic factor; resolves ER unfolded protein stress and promotes beta-cell survival.",
    },
    "NLRC5": {
        "opentargets_t1d_score": None,
        "gate": "Gate2_Allo_HLA_ClassI_or_II",
        "failure_mode_defended": "FM1_Alloreactive_T_Cell_and_DSA_Rejection",
        "recommended_edit": "Biallelic_KO",
        "hpa_islet_expression_status": "IFN_Gamma_Inducible_MHC_I_Transactivator_Descriptive_Context_Only",
        "rationale": "MHC-I transactivator; ablates cytokine-inducible MHC-I hyper-expression while preserving basal HLA-C/E (pLI=6.4e-17).",
    },
    "PTPN2": {
        "opentargets_t1d_score": 0.5104798651518475,
        "gate": "Gate4_Intrinsic_IFN_JAK_STAT_Block",
        "failure_mode_defended": "FM6_Interferon_MHC_I_Hyperexpression_and_Apoptosis",
        "recommended_edit": "Preserve_or_Stabilize_NEVER_KO",
        "hpa_islet_expression_status": "Ubiquitous_Islet_Protective_Phosphatase_Descriptive_Context_Only",
        "rationale": "Endogenous STAT1 phosphatase brake (pLI=0.990); T1D risk alleles reduce PTPN2.",
    },
    "PTPRN": {
        "opentargets_t1d_score": 0.318,
        "gate": "Gate2_Autoantigen_Stealth_TCR_or_AutoAb",
        "failure_mode_defended": "FM3_Recurrent_Autoantibodies_Exocytosis_ADCC_and_Epitope_Spreading",
        "recommended_edit": "Biallelic_KO_or_Epitope_BaseEdit",
        "hpa_islet_expression_status": "Islet_Secretory_Granule_Membrane_Descriptive_Context_Only",
        "rationale": "IA-2 autoantigen (IA-2A); LoF tolerant (pLI=8.4e-15), subject to dynamic perifusion non-inferiority gate.",
    },
    "RNLS": {
        "opentargets_t1d_score": 0.276,
        "gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "failure_mode_defended": "FM7_ER_Stress_Unfolded_Proinsulin_and_CTL_Bystander_Death",
        "recommended_edit": "Hypothesis_Arm_Only_Requires_Human_iPSC_Beta_GSIS_and_OCR_Validation",
        "hpa_islet_expression_status": "Low_Basal_Islet_Expression_Descriptive_Context_Only",
        "rationale": "CRISPR screen hit protecting mouse NIT-1 (PMID:32719542) and human iPSC-beta cells (PMID:38975343) against ER stress/CTLs (pLI=1.3e-13); MacDonald et al. 2025 (PMID:39988068 Fig. 6) report that poor dynamic GSIS persisted in both WT and RNLS-/- human iPSC-beta aggregates (no KO-induced differential GSIS deficit) alongside higher basal oxygen consumption (OCR) in RNLS-/- aggregates. Excluded from Universal Chassis A core as a metabolic precaution pending head-to-head perifusion validation.",
    },
    "SERPINB9": {
        "opentargets_t1d_score": None,
        "gate": "Gate4_Intrinsic_Granzyme_Fas_Apoptosis_Shield",
        "failure_mode_defended": "FM5_GranzymeB_and_FasL_Perforin_Cytotoxicity",
        "recommended_edit": "Inflammation_Inducible_SwapIn_at_CXCL10",
        "hpa_islet_expression_status": "Low_Endogenous_Beta_Cell_Expression_Requires_Transgenic_Armoring_Descriptive_Context_Only",
        "rationale": "Endogenous Granzyme B inhibitor; neutralizes CTL and NK granule cytotoxicity when co-expressed with PD-L1 under the endogenous CXCL10 promoter.",
    },
    "SLC30A8": {
        "opentargets_t1d_score": 0.289,
        "gate": "Gate2_Autoantigen_Stealth_TCR_or_AutoAb",
        "failure_mode_defended": "FM3_Recurrent_Autoantibodies_Exocytosis_ADCC_and_Epitope_Spreading",
        "recommended_edit": "Optional_Allelic_Series_R325W_or_Haploinsufficiency_Pending_Perifusion_Gate",
        "hpa_islet_expression_status": "High_Islet_Beta_and_Alpha_Cell_Granule_Transporter_Descriptive_Context_Only",
        "rationale": "ZnT8 autoantigen (pLI=1.9e-7). Human LoF protects against T2D (PMID:24584071) and boosts SC-islet GSIS/hypoxia survival (PMID:31676859, PMID:42078397), but complete biallelic KO abolishes insulin dense-core zinc crystallization and alters hepatic insulin clearance. Excluded from lean Chassis A core.",
    },
    "STAT1": {
        "opentargets_t1d_score": 0.0707,
        "gate": "Gate4_Intrinsic_IFN_JAK_STAT_Block",
        "failure_mode_defended": "FM6_Interferon_Signaling_vs_Epistatic_Promoter_Trap",
        "recommended_edit": "Avoid_Clonal_KO_When_Using_GAS_ISRE_Inducible_Shields",
        "hpa_islet_expression_status": "Ubiquitous_IFN_Signal_Transducer_Descriptive_Context_Only",
        "rationale": "Strongly constrained in gnomAD (pLI=1.000, LOEUF=0.302, oe_lof=0.215, oe_mis=0.385). Clonal STAT1-/- blocks IFN-gamma apoptosis (Maestas 2026 PMID:41886504) but creates an epistatic deadlock if paired with GAS/ISRE-inducible PD-L1/SERPINB9 promoters without an NF-kB bypass.",
    },
    "TGM2": {
        "opentargets_t1d_score": None,
        "gate": "Gate2_Autoantigen_Stealth_TCR_or_AutoAb",
        "failure_mode_defended": "FM4_PostTranslational_Neoantigens_HIPs_DRiPs_Deamidation",
        "recommended_edit": "Biallelic_KO",
        "hpa_islet_expression_status": "Expressed_in_Endothelial_and_Stressed_Islet_Cells_Descriptive_Context_Only",
        "rationale": "Tissue transglutaminase 2; deamidates glutamine residues on islet proteins under ER calcium stress to generate high-affinity HLA-DQ2/DQ8 neoantigens. Completely LoF tolerant (pLI=7.0e-13, oe_lof=0.701).",
    },
    "THBD": {
        "opentargets_t1d_score": 0.065,
        "gate": None,
        "failure_mode_defended": None,
        "recommended_edit": "Optional_Surface_CoDisplay_for_Intravascular_IBMIR_Sites",
        "hpa_islet_expression_status": "Endothelial_Enriched_Intronless_Gene_Descriptive_Context_Only",
        "rationale": "Thrombomodulin (CD141); intronless gene (gnomAD LoF metrics null; oe_mis=1.041) that converts thrombin into an activator of anticoagulant Protein C (PMID:38709159).",
    },
    "TNFAIP3": {
        "opentargets_t1d_score": 0.354,
        "gate": "Gate4_Intrinsic_IFN_JAK_STAT_Block",
        "failure_mode_defended": "FM6_Interferon_MHC_I_Hyperexpression_and_Apoptosis",
        "recommended_edit": "Gain_of_Function_Inducible_OE_NEVER_KO",
        "hpa_islet_expression_status": "NFkB_Inducible_Islet_Survival_Brake_Descriptive_Context_Only",
        "rationale": "A20 ubiquitin-editing NF-kB/JNK brake; extreme LoF intolerance (pLI=1.000, LOEUF=0.217); must never be knocked out.",
    },
    "TXNIP": {
        "opentargets_t1d_score": 0.195,
        "gate": "Gate4_Intrinsic_Hypoxia_Metabolic_Shield",
        "failure_mode_defended": "FM11_72h_Ischemic_Hypoxia_and_Necrotic_DAMP_Alarm",
        "recommended_edit": "Biallelic_KO_or_Verapamil_Preconditioning",
        "hpa_islet_expression_status": "Glucose_and_ER_Stress_Inducible_Beta_Cell_ProApoptotic_Factor_Descriptive_Context_Only",
        "rationale": "Mediates glucotoxic and ER-stress mitochondrial apoptosis (Chen et al. 2008 PMID:18552236; Shalev 2014 PMID:24911120); LoF tolerant (pLI=0.0009).",
    },
    "TYK2": {
        "opentargets_t1d_score": 0.583522184815907,
        "gate": "Gate4_Intrinsic_IFN_JAK_STAT_Block",
        "failure_mode_defended": "FM6_Interferon_MHC_I_Hyperexpression_and_Apoptosis",
        "recommended_edit": "Biallelic_KO_or_Selective_Kinase_Domain_BaseEdit",
        "hpa_islet_expression_status": "Ubiquitous_Type_I_IFN_Kinase_Stage3_5_Developmental_Role_Descriptive_Context_Only",
        "rationale": "Blocks IFN-alpha/beta/IL-12/23 signaling and MHC-I hyper-expression (human LoF protects against T1D; OpenTargets #9, pLI=2.8e-7), but complete Stage-0 hPSC knockout reduces Stage 4-5 endocrine precursor yield by ~30-50% (Chandra et al. 2022 PMID:36289205).",
    },
}


def finite_optional_metric(value, label: str, upper: float | None = None, *, required: bool = False):
    if value is None and not required:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite real number")
    try:
        finite = math.isfinite(value)
    except OverflowError:
        finite = False
    if not finite or value < 0 or (upper is not None and value > upper):
        raise ValueError(f"{label} is outside its finite nonnegative domain")
    return value


def main() -> None:
    ot_file = OMICS_DIR / "t1d_associated_targets_top100.json"
    ot_scores = {}
    if ot_file.exists():
        with open(ot_file, encoding="utf-8") as f:
            ot_data = json.load(f)
        if ot_data.get("disease", {}).get("id") != "MONDO_0005147":
            raise ValueError("OpenTargets source must identify MONDO_0005147")
        for row in ot_data.get("disease", {}).get("associatedTargets", {}).get("rows", []):
            sym = row.get("target", {}).get("approvedSymbol")
            if not isinstance(sym, str) or not sym.strip() or sym in ot_scores:
                raise ValueError("OpenTargets symbols must be nonempty and unique")
            ot_scores[sym] = finite_optional_metric(row.get("score"), f"OpenTargets {sym} score", 1.0, required=True)

    records = []
    skipped_gnomad_files = []
    direct_graphql_targets = []
    curated_fallback_targets = []
    null_targets = []
    for p in sorted(glob.glob(str(OMICS_DIR / "gnomad_*.json"))):
        with open(p, encoding="utf-8") as f:
            raw = json.load(f)
        if raw.get("data") is not None and not isinstance(raw["data"], dict):
            raise ValueError("gnomAD data must be an object or null")
        gene = (raw.get("data") or {}).get("gene")
        if gene is not None and not isinstance(gene, dict):
            raise ValueError("gnomAD gene must be an object or null")
        if not gene or not gene.get("symbol"):
            err_msgs = [e.get("message", "") for e in (raw.get("errors") or []) if isinstance(e, dict)]
            skipped_gnomad_files.append({
                "file": Path(p).name,
                "status": "Upstream_GraphQL_Error_Skipped_Non_Stack_Target",
                "error_messages": err_msgs,
                "note": raw.get("openislet_audit_note", "Gene object null in raw gnomAD response; excluded from 27-gene stack prior table."),
            })
            continue
        sym = gene.get("symbol")
        if not isinstance(sym, str) or Path(p).name != f"gnomad_{sym}.json":
            raise ValueError("gnomAD filename and gene symbol must match")
        gid = gene.get("gene_id")
        if not isinstance(gid, str) or re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?", gid) is None:
            raise ValueError("gnomAD source must identify a canonical human Ensembl gene")
        c = gene.get("gnomad_constraint")
        if c is not None and not isinstance(c, dict):
            raise ValueError("gnomAD constraints must be an object or null")
        c = c or {}
        for key in ("pLI", "pli", "oe_lof", "oe_lof_lower", "oe_lof_upper", "oe_mis", "oe_mis_lower", "oe_mis_upper"):
            if key in c:
                finite_optional_metric(c[key], f"{sym} {key}", 1.0 if key in {"pLI", "pli"} else None)
        for prefix in ("oe_lof", "oe_mis"):
            lower, upper = c.get(prefix + "_lower"), c.get(prefix + "_upper")
            if lower is not None and upper is not None and lower > upper:
                raise ValueError(f"{sym} {prefix} interval is reversed")
        ann = ENGINEERING_ANNOTATIONS.get(sym, {})
        finite_optional_metric(ann.get("opentargets_t1d_score"), f"{sym} curated prior", 1.0)
        pli_val = c.get("pLI") if c.get("pLI") is not None else c.get("pli")
        if sym in ot_scores:
            ot_val = ot_scores[sym]
            curated_val = ot_scores[sym]
            ot_prov = "OpenTargets_MONDO_0005147_Direct_GraphQL_Top25"
            direct_graphql_targets.append(sym)
        elif ann.get("opentargets_t1d_score") is not None:
            ot_val = None
            curated_val = ann.get("opentargets_t1d_score")
            ot_prov = "Curated_Literature_Prior_Fallback_Not_In_Top25_GraphQL"
            curated_fallback_targets.append(sym)
        else:
            ot_val = None
            curated_val = None
            ot_prov = "Not_Indexed_Null"
            null_targets.append(sym)
        records.append({
            "gene_symbol": sym,
            "ensembl_gene_id": gid,
            "opentargets_t1d_score": ot_val,
            "curated_literature_prior_t1d_score": curated_val,
            "opentargets_score_source": ot_prov,
            "opentargets_score_provenance": ot_prov,
            "gnomad_pli": pli_val,
            "gnomad_oe_lof": c.get("oe_lof"),
            "gnomad_loeuf_upper": c.get("oe_lof_upper"),
            "gnomad_oe_mis": c.get("oe_mis"),
            "hpa_islet_expression_status": ann.get("hpa_islet_expression_status"),
            "functional_gate": ann.get("gate"),
            "failure_mode_defended": ann.get("failure_mode_defended"),
            "recommended_edit": ann.get("recommended_edit"),
            "engineering_rationale": ann.get("rationale"),
        })

    out_path = OMICS_DIR / "target_prior_table.json"
    prior_payload = {
        "generated_at": "2026-10-24",
        "schema_version": "0.8.1",
        "omics_integration_note": (
            "Quantitative evaluator incorporates gnomAD v2.1.1/v4 loss-of-function constraint "
            "(gnomad_pli, gnomad_loeuf_upper) and OpenTargets/literature T1D genetic association priors. "
            "opentargets_t1d_score is populated ONLY for the 3 direct hits in the top-25 MONDO_0005147 "
            "GraphQL snapshot (INS, PTPN2, TYK2; non-top-25 genes have opentargets_t1d_score=null), while "
            "curated literature T1D association priors are stored separately in "
            "curated_literature_prior_t1d_score with explicit opentargets_score_source / "
            "opentargets_score_provenance tags (AUDIT-REF). Human Protein Atlas / GTEx islet "
            "expression fields (hpa_islet_expression_status) are retained as qualitative biological context "
            "metadata only."
        ),
        "opentargets_provenance_summary": {
            "opentargets_release": "26.09",
            "source_file": "data/omics/t1d_associated_targets_top100.json",
            "source_file_rows_returned": len(ot_scores),
            "direct_graphql_top25_count": len(direct_graphql_targets),
            "direct_graphql_top25_targets": direct_graphql_targets,
            "curated_literature_prior_fallback_count": len(curated_fallback_targets),
            "curated_literature_prior_fallback_targets": curated_fallback_targets,
            "not_indexed_null_count": len(null_targets),
            "not_indexed_null_targets": null_targets,
        },
        "skipped_gnomad_files_audit": skipped_gnomad_files,
        "targets": records,
    }

    extracted_summary_path = BASE_DIR / "data" / "extracted" / "omics_baseline_summary.json"
    summary_payload = {
        "generated_at": "2026-10-24",
        "schema_version": "0.8.1",
        "raw_opentargets_file": "data/omics/t1d_associated_targets_top100.json",
        "raw_opentargets_rows_returned": len(ot_scores),
        "raw_opentargets_symbols": sorted(ot_scores.keys()),
        "direct_graphql_top25_targets_in_prior_table": direct_graphql_targets,
        "curated_literature_prior_fallback_targets": curated_fallback_targets,
        "not_indexed_null_targets": null_targets,
        "skipped_gnomad_files_audit": skipped_gnomad_files,
        "target_count": len(records),
        "targets": records,
    }
    publish_byte_artifacts({
        out_path: (json.dumps(prior_payload, indent=2, allow_nan=False) + "\n").encode("utf-8"),
        extracted_summary_path: (json.dumps(summary_payload, indent=2, allow_nan=False) + "\n").encode("utf-8"),
    }, lock_path=BASE_DIR / "data" / ".omics-publication.lock")
    print(f"Wrote {len(records)} target priors to {out_path} and {extracted_summary_path} (direct GraphQL: {len(direct_graphql_targets)}, curated fallback: {len(curated_fallback_targets)}, null: {len(null_targets)})")


if __name__ == "__main__":
    main()
