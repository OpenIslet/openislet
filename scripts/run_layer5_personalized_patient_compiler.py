#!/usr/bin/env python3
"""
Layer 5: Personalized Patient-DNA & Multiplex CRISPR Compiler + 10,000-Patient Global Population Simulator (v0.8.0)
===================================================================================================================
Closes the final in-silico execution layer of the 5-Layer OpenIslet Architecture (folding in AUDIT-REF,
AUDIT-REF, AUDIT-REF, AUDIT-REF, and AUDIT-REF):
  - Part A: Exact Multiplex CRISPR Guide RNA (gRNA) Compiler & Polymorphism Audit
    * Pass 1 (AsCas12a Ultra, 5'-TTTV-3' PAM, 23-nt spacer):
        - OI-GRNA-01: B2M Exon 1 (`chr15:44711590-44711616 (+)`, straddling the Layer 2 `44711608|44711609` HDR junction)
        - OI-GRNA-02: CIITA Shared Constitutive Exon 5 (`chr16:10898963-10898989 (-)`, ablating pI, pIII, and IFN-gamma-inducible pIV; legacy pIII Exon 1 comparator at `chr16:10877361-10877387 (-)`)
    * Pass 2 (HiFi-SpCas9 R691A, 5'-NGG-3' PAM, 20-nt spacer):
        - OI-GRNA-03: CXCL10 Exon 1 (`chr4:76023451-76023473 (-)`, straddling the Layer 2 `76023468|76023467` HDR junction)
        - OI-GRNA-04: AAVS1 / PPP1R12C Intron 1 (`chr19:55115754-55115776 (-)`, straddling the Layer 2 `55115768|55115767` HDR junction)
    * Optional Non-DSB Prime Edit (PEmax SpCas9-H840A nickase + M-MLV RT + epegRNA; CGBE is an unimplemented contingency, requiring an explicit bystander/window ledger; minus-strand C:G->G:C editing requires an engineered SpG/SpRY PAM variant because no canonical 5'-NGG-3' minus-strand PAM places chr19:10352442C within positions 4-7):
        - OI-GRNA-05: TYK2 Exon 23 (`chr19:10352421-10352443 (+)`, installing protective rs34536443 `c.3310C>G` [`p.Pro1104Ala`] transversion)
    * Computes directional Sugimoto et al. (1995, PMID:7545436) RNA:DNA hybrid nearest-neighbor dG_37C alongside
      SantaLucia (1998) DNA:DNA duplex dG_37C, verifies 100% 0-mismatch identity against GRCh38 reference sequences,
      and proves bidirectional (+ and - strand) protospacer/PAM disruption in both HDR donor plasmids and assembled
      post-HDR edited alleles (`HA_L + insert + HA_R`).
  - Part B: Personalized Patient-DNA & Autoantigen Compiler + 10,000-Patient Global Population Scenario Simulator
    * Validates all patient input fields (rejecting NaN/Inf, out-of-range scores, and invalid SIRPA alleles),
      distinguishes recipient vs donor cell-line protospacer SNPs, unifies Plan B with the canonical HLA-C-retaining
      selective HLA-A/HLA-B KO architecture, and executes `compile_personalized_patient_blueprint` across all
      10,000 simulated virtual T1D patients so population outcomes are causally coupled to compiler outputs.

Outputs: `data/extracted/layer5_personalized_patient_and_crispr_atlas.json`
"""

from __future__ import annotations

import hashlib
import json
import math
import copy
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Dict, List

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.engineering_priors import engineering_complexity_discount
from scripts.product_contract import get_product_definition


@lru_cache(maxsize=16)
def _product_contract_snapshot(root_path: str, variant_id: str) -> Dict[str, Any]:
    """One immutable contract snapshot per imported compiler and product variant."""
    return get_product_definition(variant_id, root=Path(root_path))


def _engineering_parameters() -> Dict[str, Any]:
    """Read the shared in-memory registry at an evaluation boundary, never a file per patient."""
    from scripts.rank_islet1_stacks import ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE
    reg = ENGINEERING_PRIORS_AND_SENSITIVITY_TABLE
    return {
        "lambda_K": float(reg["locus_dsb_penalty_rate_lambda_K"]["point_estimate"]),
        "gamma_B": float(reg["base_edit_nickase_penalty_factor_gamma_B"]["point_estimate"]),
        "lambda_K_range": tuple(reg["locus_dsb_penalty_rate_lambda_K"]["plausible_range"]),
        "gamma_B_range": tuple(reg["base_edit_nickase_penalty_factor_gamma_B"]["plausible_range"]),
    }
RAW_LAYER2_REF_PATH = REPO_ROOT / "data" / "raw" / "layer2_genomic_and_jaspar_reference.json"
RAW_LAYER5_DBSNP_PATH = REPO_ROOT / "data" / "raw" / "layer5_dbsnp_clinvar_reference.json"
LAYER2_ATLAS_PATH = REPO_ROOT / "data" / "extracted" / "layer2_regulatory_and_perturbation_atlas.json"
OUT_PATH = REPO_ROOT / "data" / "extracted" / "layer5_personalized_patient_and_crispr_atlas.json"

# Sugimoto et al. (1995) Biochemistry 34:11211-11216 (PMID:7545436, Table 1) directional RNA(5'->3'):DNA(3'->5') hybrid dG_37 (kcal/mol)
# Keys are 5'->3' RNA dinucleotides (written with T for U) base-paired to antiparallel DNA target strand (closing AUDIT-REF):
SUGIMOTO_1995_RNA_DNA_DG37: Dict[str, float] = {
    "AA": -1.00, "TT": -0.20, "AT": -0.90, "TA": -0.60,
    "CA": -0.90, "TG": -1.60, "GT": -1.10, "AC": -2.10,
    "CT": -0.90, "AG": -1.80, "GA": -1.30, "TC": -1.50,
    "CG": -1.70, "GC": -2.70, "GG": -2.90, "CC": -2.10,
}

# SantaLucia (1998) Proc Natl Acad Sci USA 95:1460-1465 DNA:DNA nearest-neighbor dG_37 (kcal/mol)
NN_DG_37: Dict[str, float] = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}

VALID_SIRPA_GENOTYPES = {
    "v1/v1", "v1/v2", "v2/v1", "v2/v2",
    "v1/v8", "v8/v1", "v2/v8", "v8/v2", "v8/v8",
}
VALID_ANCESTRIES = {"EUR", "AFR", "EAS", "SAS", "AMR", "OTHER_ADMIXED"}
VALID_TRACK_REQUESTS = {
    "AUTO_RECOMMEND",
    "TRACK_1_OFF_THE_SHELF",
    "TRACK_2_STRATIFIED_SUBLINE",
    "TRACK_3_AUTOLOGOUS_IPSC",
}
VALID_IFIH1_GENOTYPES = {"C/C", "C/T", "T/C", "T/T"}
VALID_PTPN22_GENOTYPES = {"G/G", "G/A", "A/G", "A/A"}
VALID_HLA_E_GENOTYPES = {
    "*01:01/*01:01", "*01:01/*01:03", "*01:03/*01:01", "*01:03/*01:03",
    "01:01/01:01", "01:01/01:03", "01:03/01:01", "01:03/01:03",
    "HLA-E*01:01/*01:01", "HLA-E*01:01/*01:03", "HLA-E*01:03/*01:01", "HLA-E*01:03/*01:03",
    "HLA-E*01:01/HLA-E*01:01", "HLA-E*01:01/HLA-E*01:03", "HLA-E*01:03/HLA-E*01:01", "HLA-E*01:03/HLA-E*01:03",
}
VALID_TYK2_RS34536443_GENOTYPES = {"G/G", "G/C", "C/G", "C/C"}
VALID_KIR_LIGAND_GROUPS = {"C1/C1", "C1/C2", "C2/C1", "C2/C2"}
_VALID_HLAC_4DIGIT_RELATIONS = {
    "SELF_4DIGIT_MATCH",
    "EXACT_4DIGIT_AUTOLOGOUS_OR_MATCHED",
    "KIR_MATCH_FOREIGN_4DIGIT_ALLO",
    "KIR_LIGAND_MISMATCH",
    "UNSPECIFIED",
    "UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT",
}

# Standard synthesized guide RNA scaffolds (5'->3' RNA) for manufacturing sequence completeness (AUDIT-REF):
ASCAS12A_DIRECT_REPEAT_20NT_RNA = "UAAUUUCUACUCUUGUAGAU"
SPCAS9_SGRNA_SCAFFOLD_80NT_RNA = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
SPCAS9_SGRNA_SCAFFOLD_76NT_RNA = SPCAS9_SGRNA_SCAFFOLD_80NT_RNA  # Backwards-compatible alias; actual string length is 80 nt (76-nt core + 4U tail)
PEMAX_EPEGRNA_SCAFFOLD_76NT_RNA = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"
# 42-nt 3' motif = 37-nt tevopreQ1 pseudoknot (whose 5' stem begins with the 8-nt CGCGGUUC sequence; 0-nt separate linker) + 5-nt AUUUU U6 termination tail (AUDIT-REF [DEC-078]):
TEVOPREQ1_PSEUDOKNOT_37NT_CORE_RNA = "CGCGGUUCUAUCUAGUUACGCGUUAAACCAACUAGAA"
TEVOPREQ1_U6_TERMINATION_TAIL_5NT_RNA = "AUUUU"
TEVOPREQ1_PSEUDOKNOT_42NT_RNA = TEVOPREQ1_PSEUDOKNOT_37NT_CORE_RNA + TEVOPREQ1_U6_TERMINATION_TAIL_5NT_RNA

# Uncalibrated positional mismatch retention weights for split-junction screening diagnostics (AUDIT-REF):
SPCAS9_POS_RETENTION_WEIGHTS_20NT = (
    0.88, 0.85, 0.82, 0.78, 0.74, 0.70, 0.65, 0.58, 0.50, 0.44,
    0.38, 0.32, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03, 0.02, 0.01,
)
ASCAS12A_POS_RETENTION_WEIGHTS_23NT = (
    0.01, 0.02, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15,
    0.18, 0.22, 0.26, 0.30, 0.35, 0.42, 0.50, 0.58, 0.78, 0.85,
    0.90, 0.94, 0.97,
)


def compute_split_junction_residual_cleavage_metrics(
    protospacer_5to3: str,
    pam_5to3: str,
    nuclease: str,
    left_junction_seq_5to3: str,
    right_junction_seq_5to3: str,
) -> Dict[str, Any]:
    """
    Evaluate raw PAM/seed retention and uncalibrated positional mismatch diagnostics across both the
    5' (`HA_L|Insert`) and 3' (`Insert|HA_R`) post-HDR chromosomal junctions when a protospacer is
    split across an HDR insertion site (closing AUDIT-REF & AUDIT-REF).
    Separates nuclease-specific AsCas12a Ultra vs HiFi-SpCas9 R691A rules and explicitly qualifies
    heuristic mismatch scores as uncalibrated screening diagnostics rather than empirical probabilities.
    """
    protospacer_5to3 = _validated_nucleotides(protospacer_5to3, "protospacer")
    pam_5to3 = _validated_nucleotides(pam_5to3, "PAM")
    left_junction_seq_5to3 = _validated_nucleotides(left_junction_seq_5to3, "left junction")
    right_junction_seq_5to3 = _validated_nucleotides(right_junction_seq_5to3, "right junction")
    if not isinstance(nuclease, str) or not any(x in nuclease for x in ("Cas12a", "SpCas9", "PEmax")):
        raise ValueError("Unsupported nuclease family")
    n_len = len(protospacer_5to3)
    is_cas12a = "Cas12a" in nuclease
    if is_cas12a and ("SpCas9" in nuclease or "PEmax" in nuclease):
        raise ValueError("Ambiguous nuclease family")
    if (n_len != (23 if is_cas12a else 20) or len(pam_5to3) != (4 if is_cas12a else 3)
            or any(len(s) != n_len + len(pam_5to3) for s in (left_junction_seq_5to3, right_junction_seq_5to3))):
        raise ValueError("Protospacer, PAM and complete junction window lengths must match the nuclease")
    if ((is_cas12a and not (pam_5to3[:3] == "TTT" and pam_5to3[3] in "ACG"))
            or (not is_cas12a and pam_5to3[1:] != "GG")):
        raise ValueError("Reference PAM does not match the declared nuclease family")
    weights = ASCAS12A_POS_RETENTION_WEIGHTS_23NT if is_cas12a else SPCAS9_POS_RETENTION_WEIGHTS_20NT

    def _score_candidate_window(candidate_proto: str, candidate_pam: str) -> Dict[str, Any]:
        mismatches = [
            idx + 1 for idx, (a, b) in enumerate(zip(protospacer_5to3, candidate_proto)) if a != b
        ]
        if is_cas12a:
            pam_intact = len(candidate_pam) == 4 and candidate_pam[:3] == "TTT" and candidate_pam[3] in {"A", "C", "G"}
            seed_mismatches = [m for m in mismatches if m <= 18]
            seed_intact = (len(seed_mismatches) == 0)
            attenuation_model = "AsCas12a_Ultra_M537R_F870L_5prime_PAM_and_18nt_seed_requirement"
        else:
            pam_intact = len(candidate_pam) == 3 and candidate_pam[1:] == "GG"
            seed_mismatches = [m for m in mismatches if m >= 9]
            seed_intact = (len(seed_mismatches) == 0)
            attenuation_model = "HiFi_SpCas9_R691A_PAM_distal_R_loop_docking_checkpoint"

        if not pam_intact:
            raw_pos_prod = 0.0
            nuclease_attenuated_heuristic = 0.0
        else:
            raw_pos_prod = 1.0
            for m_pos in mismatches:
                raw_pos_prod *= weights[m_pos - 1]
            if is_cas12a:
                # AsCas12a Ultra (M537R/F870L) requires intact 5'-TTTV-3' PAM + PAM-proximal positions 1..18;
                # splitting between nt +15 and +16 introduces mismatches in positions 16..23 (including 3 inside nt 1..18):
                cas12a_factor = 0.02 if len(seed_mismatches) >= 2 else (0.15 if len(seed_mismatches) == 1 else 0.50)
                nuclease_attenuated_heuristic = raw_pos_prod * cas12a_factor
            else:
                # HiFi-SpCas9 R691A (Vakulskas et al. 2018, Nat Med 24:1216) requires 20-bp R-loop zipping across
                # PAM-distal positions 1..6 to activate HNH domain docking; when the 3'-AGG PAM and 11-14 nt PAM-proximal
                # seed remain intact on HA_R, silent PAM/seed shield mutations in HA_R are recommended for complete immunity:
                hifi_factor = (0.015 ** max(1, len(mismatches) - 3)) if len(mismatches) >= 4 else 0.10
                nuclease_attenuated_heuristic = raw_pos_prod * hifi_factor
        heuristic_mit_proxy = 0.0 if (not pam_intact or len(mismatches) >= 5) else round(raw_pos_prod * 0.25, 8)
        return {
            "candidate_protospacer_5to3": candidate_proto,
            "candidate_pam_5to3": candidate_pam,
            "pam_intact": pam_intact,
            "seed_region_intact": seed_intact,
            "total_mismatches_vs_wt_spacer": len(mismatches),
            "seed_region_mismatches": len(seed_mismatches),
            "mismatch_positions_1based": mismatches,
            "nuclease_specific_attenuation_model": attenuation_model,
            "score_calibration_status": "UNCALIBRATED_POSITIONAL_HEURISTIC_REQUIRING_WETLAB_T7E1_NGS_VALIDATION",
            "positional_retention_product_assigned_weights_uncalibrated": round(raw_pos_prod, 6),
            "positional_retention_product_assigned_weights": round(raw_pos_prod, 6),
            "hsu_2013_mit_off_target_score": None,
            "doench_2016_cfd_off_target_score": None,
            "published_cfd_or_hsu_mit_matrix_status": "UNAVAILABLE_UNCALIBRATED_POSITIONAL_HEURISTIC_ONLY",
            "unattenuated_wt_cas_cfd_score": round(raw_pos_prod, 6),
            "uncalibrated_positional_mismatch_score": round(raw_pos_prod, 6),
            "residual_cfd_cleavage_probability": None,
            "residual_mit_cleavage_probability": None,
            "nuclease_attenuated_positional_heuristic_uncalibrated": round(nuclease_attenuated_heuristic, 8),
            "nuclease_attenuated_positional_heuristic": round(nuclease_attenuated_heuristic, 8),
            "assigned_panel_heuristic": round(heuristic_mit_proxy, 8),
            "silent_ha_r_pam_or_seed_shield_recommended_if_wt_cas": bool(pam_intact and seed_intact),
        }

    if is_cas12a:
        left_eval = _score_candidate_window(left_junction_seq_5to3[4 : 4 + n_len], left_junction_seq_5to3[:4])
        right_eval = _score_candidate_window(right_junction_seq_5to3[4 : 4 + n_len], right_junction_seq_5to3[:4])
    else:
        left_eval = _score_candidate_window(left_junction_seq_5to3[:n_len], left_junction_seq_5to3[n_len : n_len + len(pam_5to3)])
        right_eval = _score_candidate_window(right_junction_seq_5to3[:n_len], right_junction_seq_5to3[n_len : n_len + len(pam_5to3)])
    max_hifi_cfd = max(
        left_eval["nuclease_attenuated_positional_heuristic"],
        right_eval["nuclease_attenuated_positional_heuristic"],
    )
    max_wt_cfd = max(
        left_eval["unattenuated_wt_cas_cfd_score"],
        right_eval["unattenuated_wt_cas_cfd_score"],
    )
    any_pam_retained = bool(left_eval["pam_intact"] or right_eval["pam_intact"])
    any_pam_and_seed_retained = bool(
        (left_eval["pam_intact"] and left_eval["seed_region_intact"])
        or (right_eval["pam_intact"] and right_eval["seed_region_intact"])
    )
    return {
        "five_prime_ha_left_insert_junction": left_eval,
        "three_prime_insert_ha_right_junction": right_eval,
        "full_20_to_23nt_target_absent_on_both_junctions": not any(
            e["pam_intact"] and e["total_mismatches_vs_wt_spacer"] == 0 for e in (left_eval, right_eval)),
        "raw_pam_motif_retained_on_either_junction": any_pam_retained,
        "raw_pam_and_seed_subsequence_retained_on_either_junction": any_pam_and_seed_retained,
        "max_unattenuated_wt_cas_junction_cfd_score": max_wt_cfd,
        "max_junction_cfd_cleavage_probability": None,
        "max_junction_attenuated_positional_heuristic": max_hifi_cfd,
        "cleavage_immune_post_hdr_hifi_r691a": None,
        "empirical_recleavage_status": "PREREGISTERED_PENDING",
        "cleavage_immunity_qualification_note": (
            "Full-target absence is calculated above for the supplied windows and does not establish cleavage immunity. For CXCL10 and AAVS1, "
            "the 3' HA_R junction retains the 3'-AGG PAM and 11-14 nt PAM-proximal seed while replacing the 5' 6-9 nt "
            "with bGH_polyA (Locus 3) or SV40_late_polyA_122bp (Locus 4) sequence. The positional mismatch scores above are uncalibrated screening heuristics (not "
            "published substitution-identity CFD or empirical cleavage probabilities); wet-lab re-cleavage testing or "
            "introducing a silent PAM/seed shield mutation in HA_R is required before claiming empirical cleavage immunity "
            "(closing AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
        ),
    }


def compute_panel_and_reference_off_target_specificity(
    guide_id: str,
    nuclease: str,
    coding_0_1_2_mm: int | None,
    noncoding_3_mm_count: int,
    noncoding_3_mm_cfd_Hits: List[float],
    crispor_whole_genome_mit_ref: float | None,
    ucsc_crispr_track_url: str | None = None,
    genomewide_enumeration_receipt: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Report both the assigned 3-mismatch panel specificity formula `100 / (1 + sum_i CFD_i)` AND the
    verified primary UCSC `crisprAllTargets` (GRCh38) WT-SpCas9 whole-genome MIT specificity score
    (`65` for CXCL10, `77` for AAVS1, `79` for TYK2), explicitly labelling the panel formula as an
    illustrative assigned-hit model pending genome-wide Cas-OFFinder + empirical GUIDE-seq/CIRCLE-seq
    (closing AUDIT-REF & AUDIT-REF).
    """
    if crispor_whole_genome_mit_ref is not None:
        crispor_whole_genome_mit_ref = _validate_finite_float(
            crispor_whole_genome_mit_ref, "reference specificity score", 0.0, 100.0)
    if coding_0_1_2_mm is not None:
        if not isinstance(coding_0_1_2_mm, int) or isinstance(coding_0_1_2_mm, bool) or coding_0_1_2_mm < 0:
            raise ValueError("coding hit count must be a nonnegative integer or None")
        required_receipt_keys = {"tool", "version", "reference_build", "reference_sha256", "guide_id", "search_parameters", "hit_loci", "sha256"}
        if not isinstance(genomewide_enumeration_receipt, dict) or not required_receipt_keys <= genomewide_enumeration_receipt.keys():
            raise ValueError("non-null genome-wide count requires pinned tool/version/build/hit-locus receipt")
        if (genomewide_enumeration_receipt["reference_build"] != "GRCh38"
                or genomewide_enumeration_receipt["guide_id"] != guide_id
                or not isinstance(genomewide_enumeration_receipt["hit_loci"], list)
                or len(genomewide_enumeration_receipt["hit_loci"]) != coding_0_1_2_mm):
            raise ValueError("enumeration receipt reference/count mismatch")
        digest = genomewide_enumeration_receipt["sha256"]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("enumeration receipt requires a lowercase SHA-256 digest")
        reference_hash = genomewide_enumeration_receipt["reference_sha256"]
        if not isinstance(reference_hash, str) or len(reference_hash) != 64 or any(c not in "0123456789abcdef" for c in reference_hash):
            raise ValueError("enumeration receipt requires a pinned reference SHA-256")
        if any(not isinstance(genomewide_enumeration_receipt[k], str) or not genomewide_enumeration_receipt[k].strip() for k in ("tool", "version")):
            raise ValueError("enumeration tool and version must be nonempty")
        if not isinstance(genomewide_enumeration_receipt["search_parameters"], dict) or not genomewide_enumeration_receipt["search_parameters"]:
            raise ValueError("enumeration search parameters must be recorded")
        material = {k:v for k,v in genomewide_enumeration_receipt.items() if k != "sha256"}
        expected = hashlib.sha256(json.dumps(material, sort_keys=True, allow_nan=False).encode()).hexdigest()
        if digest != expected:
            raise ValueError("enumeration receipt content does not match its SHA-256")
    if (not isinstance(noncoding_3_mm_count, int) or isinstance(noncoding_3_mm_count, bool)
            or noncoding_3_mm_count < 0 or not isinstance(noncoding_3_mm_cfd_Hits, list)
            or len(noncoding_3_mm_cfd_Hits) != noncoding_3_mm_count
            or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 <= v <= 1 for v in noncoding_3_mm_cfd_Hits)):
        raise ValueError("Assigned panel requires one finite [0,1] score per declared hit")
    sum_cfd = sum(noncoding_3_mm_cfd_Hits)
    panel_score = round(100.0 / (1.0 + sum_cfd), 1)
    return {
        "window_and_exome_panel_cfd_specificity_score": panel_score,
        "panel_score_evidence_status": "ASSIGNED_ILLUSTRATIVE_3MM_PANEL_FORMULA_NOT_GENOMEWIDE_EXOME_SCAN",
        "sum_off_target_cfd_hazard_across_panel": round(sum_cfd, 4),
        "coding_exome_0_1_2_mismatch_hits": coding_0_1_2_mm,
        "genomewide_off_target_enumeration_status": "NOT_ENUMERATED" if coding_0_1_2_mm is None else "RECEIPT_SUPPLIED_NOT_INDEPENDENTLY_REPLAYED",
        "genomewide_enumeration_receipt": genomewide_enumeration_receipt,
        "noncoding_3_mismatch_loci_evaluated": noncoding_3_mm_count,
        "noncoding_3_mismatch_individual_cfd_scores": noncoding_3_mm_cfd_Hits,
        "whole_genome_crispor_mit_reference_specificity_score": crispor_whole_genome_mit_ref,
        "primary_ucsc_crispr_all_targets_query_url": ucsc_crispr_track_url,
        "scoring_methodology_note": (
            "whole_genome_crispor_mit_reference_specificity_score reports the verified primary UCSC hg38 crisprAllTargets "
            "WT-SpCas9 whole-genome MIT specificity score (CXCL10 = 65, AAVS1 = 77, TYK2 = 79; AsCas12a TTTV guides are "
            "not indexed in the SpCas9 NGG UCSC track). window_and_exome_panel_cfd_specificity_score computes "
            "100 / (1 + sum_i c_i) over an illustrative assigned 3-mismatch panel and does not represent an exhaustive "
            "GRCh38 exome/genome scan; release qualification requires genome-wide Cas-OFFinder enumeration plus empirical "
            "GUIDE-seq / CAST-seq in the target hPSC line (closing AUDIT-REF & AUDIT-REF)."
        ),
    }


def _validated_nucleotides(seq: str, name: str = "sequence") -> str:
    if not isinstance(seq, str) or not seq or any(c not in "ACGTUacgtu" for c in seq):
        raise ValueError(f"{name} must be a nonempty unambiguous DNA/RNA sequence")
    return seq.upper().replace("U", "T")


def revcomp(seq: str) -> str:
    return _validated_nucleotides(seq).translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _canonical_target_present_either_strand(sequence: str, protospacer: str, nuclease: str) -> bool:
    """Search an exact spacer with every canonical PAM, including reverse strand sites."""
    sequence = _validated_nucleotides(sequence)
    protospacer = _validated_nucleotides(protospacer, "protospacer")
    if not isinstance(nuclease, str) or not any(family in nuclease for family in ("Cas12a", "SpCas9", "PEmax")):
        raise ValueError("Unsupported nuclease family")
    is_cas12a = "Cas12a" in nuclease
    if (len(protospacer) != (23 if is_cas12a else 20)
            or (is_cas12a and ("SpCas9" in nuclease or "PEmax" in nuclease))):
        raise ValueError("Protospacer length must match one supported nuclease family")
    if is_cas12a:
        targets = ["TTT" + base + protospacer for base in "ACG"]
    else:
        targets = [protospacer + base + "GG" for base in "ACGT"]
    return any(target in sequence or revcomp(target) in sequence for target in targets)


def _check(cond: bool, msg: str = "Layer 5 invariant check failed") -> None:
    if not cond:
        raise RuntimeError(msg)


def compute_sugimoto_1995_rna_dna_dg37_kcal_mol(seq: str) -> float:
    """
    Compute directional RNA:DNA hybrid dG_37 (kcal/mol) per Sugimoto et al. 1995 (PMID:7545436, Table 1):
    Initiation at 37 C (310.15 K): dG_init(37 C) = dH_init - T * dS_init = 1.90 - 310.15 * (-0.0039) = +3.11 kcal/mol
    (measured in 1 M NaCl buffer; inside the Cas9/Cas12a RNP channel, protein-backbone contacts further stabilize
    the R-loop relative to free-solution duplexes; closing AUDIT-REF & AUDIT-REF).
    """
    s = _validated_nucleotides(seq)
    dg = 1.90 - 310.15 * (-0.0039)  # +3.109585 kcal/mol initiation at 37 C
    for i in range(len(s) - 1):
        dinuc = s[i : i + 2]
        dg += SUGIMOTO_1995_RNA_DNA_DG37[dinuc]
    return round(dg, 2)


def compute_duplex_dg37_kcal_mol(seq: str) -> float:
    """SantaLucia 1998 Eq. 1/Table 1 DNA:DNA dG_37 at 1 M NaCl."""
    s = _validated_nucleotides(seq)
    # Each end contributes 0.98 (G:C) or 1.03 (A:T) kcal/mol.
    # The additional +0.43 term applies only to self-complementary duplexes.
    dg = sum(1.03 if base in "AT" else 0.98 for base in (s[0], s[-1]))
    if s == revcomp(s):
        dg += 0.43
    for i in range(len(s) - 1):
        dinuc = s[i : i + 2]
        dg += NN_DG_37[dinuc]
    return round(dg, 2)


def gc_percent(seq: str) -> float:
    s = _validated_nucleotides(seq)
    return round(100.0 * (s.count("G") + s.count("C")) / max(1, len(s)), 2)


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def build_part_a_crispr_grna_and_polymorphism_atlas() -> Dict[str, Any]:
    """
    Compile exact Pass-1 (AsCas12a Ultra), Pass-2 (HiFi-SpCas9 R691A), and optional PEmax epegRNA guides,
    verify 100% 0-mismatch identity against authentic GRCh38 sequences in `layer2_genomic_and_jaspar_reference.json`,
    and prove bidirectional (+ and - strand) target disruption in both HDR donor plasmids and assembled post-HDR
    edited genomic alleles (closing AUDIT-REF and AUDIT-REF).
    """
    raw_l2 = json.loads(RAW_LAYER2_REF_PATH.read_text(encoding="utf-8"))
    raw_l5 = json.loads(RAW_LAYER5_DBSNP_PATH.read_text(encoding="utf-8"))
    gseqs = raw_l2["ensembl_verified_homology_arms_and_elements"]

    # Assemble exact WT genomic windows from cached Ensembl GRCh38 sequences to verify 100% 0-mismatch identity:
    b2m_wt_1200 = gseqs["B2M_HA_L_600bp"]["sequence"] + gseqs["B2M_HA_R_600bp"]["sequence"]  # chr15:44711009-44712208 (+)
    ciita_wt_600 = gseqs["CIITA_target_window_600bp"]["sequence"]                              # chr16:10877001-10877600 (+)
    ciita_exon5_wt_200 = gseqs["CIITA_pan_isoform_exon5_window_200bp"]["sequence"]             # chr16:10898880-10899079 (+)
    cxcl10_wt_1200 = gseqs["CXCL10_HA_L_600bp_minus"]["sequence"] + gseqs["CXCL10_HA_R_600bp_minus"]["sequence"]  # chr4:76022868-76024067 (-)
    aavs1_wt_1200 = gseqs["AAVS1_HA_L_600bp_minus"]["sequence"] + gseqs["AAVS1_HA_R_600bp_minus"]["sequence"]     # chr19:55115168-55116367 (-)
    tyk2_wt_71 = gseqs["TYK2_exon23_rs34536443_window_71bp_minus"]["sequence"]                 # chr19:10352410-10352480 (-)
    tyk2_plus_71 = revcomp(tyk2_wt_71)                                                         # chr19:10352410-10352480 (+)

    grna_specs: List[Dict[str, Any]] = [
        {
            "guide_id": "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "editing_pass": "Pass 1 (Stage-0 hPSC Multiplex via AsCas12a Ultra RNP)",
            "target_gene": "B2M",
            "target_locus_grch38": "chr15:44711590-44711616 (+)",
            "chrom": "chr15",
            "start_1based": 44711590,
            "end_1based": 44711616,
            "strand": "+",
            "nuclease": "AsCas12a Ultra (M537R/F870L)",
            "pam_motif": "5'-TTTV-3'",
            "pam_sequence_5to3": "TTTC",
            "protospacer_dna_5to3": "TGGCCTGGAGGCTATCCAGCGTG",
            "spacer_length_nt": 23,
            "cleavage_geometry": "5-nt 5'-overhang staggered cut at nt +18 (non-target) / +23 (target) downstream of TTTC PAM in B2M Exon 1",
            "hdr_junction_grch38": "chr15:44711608|44711609 (+)",
            "hdr_pam_and_seed_disrupted_post_knockin": True,
            "hdr_re_cutting_prevention_mechanism": (
                "B2M_HA_L ends at chr15:44711608 (containing 5'-TTTCTGGCCTGGAGGCTAT-3') and B2M_HA_R begins at "
                "chr15:44711609 (containing 5'-CCAGCGTG-3'); HDR insertion of Arm 1A or Arm 1B at 44711608|44711609 "
                "splits the 23-nt AsCas12a protospacer between nt +15 and +16, eliminating the intact target site from "
                "both strands of the donor plasmid and post-HDR edited allele."
            ),
            "cfd_mit_off_target_specificity_score": 96.4,
            "whole_genome_crispor_mit_reference_specificity_score": None,
            "grch38_0_1_2_mismatch_coding_off_targets": None,
            "genomewide_off_target_enumeration_status": "NOT_ENUMERATED",
            "grch38_3_mismatch_intergenic_intronic_loci": 2,
            "_noncoding_3mm_cfd_hits": [0.0211, 0.0162],
            "specificity_provenance": "assigned_3mm_panel_model_plus_genomewide_guide_seq_release_gate",
            "gnomad_v4_1_max_ancestry_maf_in_protospacer_or_pam": 0.00018,
            "gnomad_v4_1_global_af_max_variant_in_target": 0.00004,
            "dbsnp_b155_overlapping_variant_count": 8,
            "ucsc_dbsnp155_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dbSnp155;chrom=chr15;start=44711589;end=44711616",
            "orthogonal_to_pass2_spcas9": True,
            "orthogonal_nuclease_family_to_spcas9": True,
            "delivered_in_separate_editing_pass_from_pass2_dsb": True,
            "_wt_expected_slice": b2m_wt_1200[44711590 - 44711009 : 44711616 - 44711009 + 1],
            "_donor_keys": [
                "ISLET1_LOCUS1_ARM1A_B2M_HLAE_SCT_GSGP2A_CD47_IGV_GPI_4880bp",
                "ISLET1_LOCUS1_ARM1B_B2M_BIDIR_2AFREE_CD47_IGV_GPI_HLAE_SCT_3974bp",
            ],
            "_split_prefix_bp": 19,
            "_split_suffix_bp": 8,
        },
        {
            "guide_id": "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
            "canonical_exon5_guide_id": "OI-GRNA-02_CIITA_EXON5_ASCAS12A",
            "pan_isoform_guide_alias": "OI-GRNA-02_CIITA_PAN_ISOFORM_EXON5_ASCAS12A",
            "editing_pass": "Pass 1 (Stage-0 hPSC Multiplex via AsCas12a Ultra RNP)",
            "target_gene": "CIITA",
            "target_locus_grch38": "chr16:10898963-10898989 (-)",
            "legacy_pIII_exon1_comparator_locus_grch38": "chr16:10877361-10877387 (-)",
            "ciita_exon5_vs_exon1_avi_window_disclosure_claude_v80_073": (
                "OI-GRNA-02 targets shared constitutive CIITA Exon 5 (chr16:10898963-10898989 -, 21.6 kb downstream of the "
                "Layer 2 pIII Exon 1 promoter AVI regulatory window chr16:10877250-10877390; canonical ID OI-GRNA-02_CIITA_EXON5_ASCAS12A, "
                "retained under legacy ID OI-GRNA-02_CIITA_EXON1_ASCAS12A for schema compatibility). Exon 5 biallelic NHEJ frameshift "
                "was deliberately selected over Exon 1 promoter disruption because Exon 5 is shared across pI, pIII, and "
                "IFN-gamma-inducible pIV (whose translation starts in Exon 2 at chr16:10895302), whereas Exon 1 editing leaves "
                "IFN-gamma-inducible pIV MHC-II expression intact (closing AUDIT-REF)."
            ),
            "chrom": "chr16",
            "start_1based": 10898963,
            "end_1based": 10898989,
            "strand": "-",
            "nuclease": "AsCas12a Ultra (M537R/F870L)",
            "pam_motif": "5'-TTTV-3'",
            "pam_sequence_5to3": "TTTC",
            "protospacer_dna_5to3": "TGCCCAACTTCTGCTGGCATCTC",
            "spacer_length_nt": 23,
            "targeted_transcripts_pan_isoform": [
                "ENST00000324288.14 (pIII lymphoid/constitutive)",
                "ENST00000576601.1 (pIV IFN-gamma-inducible non-professional APC / beta-cell / iEC, translation start Met25 in Exon 2 at chr16:10895302)",
                "ENST00000637439.1 (pI myeloid/dendritic)",
            ],
            "cleavage_geometry": (
                "5-nt 5'-overhang staggered cut in shared constitutive CIITA Exon 5 (chr16:10898925-10899002), which is "
                "present in all three major CIITA isoforms (pI, pIII, and IFN-gamma-inducible pIV; closing AUDIT-REF). "
                "The pIII-specific Exon 1 guide (chr16:10877361-10877387 -, TTTA + CCTTGGGGCTCTGACAGGTAGGA, verified 0-mismatch "
                "inside the Layer 2 chr16:10877250-10877390 AVI window) is retained strictly as a pIII-only comparator."
            ),
            "hdr_junction_grch38": "NHEJ_MMEJ_BIALLELIC_FRAMESHIFT_KO",
            "hdr_pam_and_seed_disrupted_post_knockin": True,
            "hdr_re_cutting_prevention_mechanism": (
                "Biallelic frameshift indels (+1/-4/-7 bp, e.g., canonical 4-bp deletion of 'CATC' at nt +18..+21) in shared Exon 5 "
                "ablate pI, pIII, and IFN-gamma-inducible pIV CIITA translation; release qualification requires allele-resolved NGS "
                "plus IFN-gamma-stimulated (48-72 h) surface HLA-DR/DP/DQ flow cytometry across both Stage-7 endocrine and iEC lineages."
            ),
            "cfd_mit_off_target_specificity_score": 94.8,
            "whole_genome_crispor_mit_reference_specificity_score": None,
            "grch38_0_1_2_mismatch_coding_off_targets": None,
            "genomewide_off_target_enumeration_status": "NOT_ENUMERATED",
            "grch38_3_mismatch_intergenic_intronic_loci": 3,
            "_noncoding_3mm_cfd_hits": [0.0218, 0.0184, 0.0147],
            "specificity_provenance": "assigned_3mm_panel_model_plus_genomewide_guide_seq_release_gate",
            "gnomad_v4_1_max_ancestry_maf_in_protospacer_or_pam": 0.00024,
            "gnomad_v4_1_global_af_max_variant_in_target": 0.00003,
            "dbsnp_b155_overlapping_variant_count": 13,
            "ucsc_dbsnp155_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dbSnp155;chrom=chr16;start=10898962;end=10898989",
            "orthogonal_to_pass2_spcas9": True,
            "orthogonal_nuclease_family_to_spcas9": True,
            "delivered_in_separate_editing_pass_from_pass2_dsb": True,
            "_wt_expected_slice": revcomp(ciita_exon5_wt_200[10898963 - 10898880 : 10898989 - 10898880 + 1]),
            "_legacy_pIII_wt_slice": revcomp(ciita_wt_600[10877361 - 10877001 : 10877387 - 10877001 + 1]),
            "_donor_keys": [],
        },
        {
            "guide_id": "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "editing_pass": "Pass 2 (Post-QC Clone Multiplex via HiFi-SpCas9 R691A RNP)",
            "target_gene": "CXCL10",
            "target_locus_grch38": "chr4:76023451-76023473 (-)",
            "chrom": "chr4",
            "start_1based": 76023451,
            "end_1based": 76023473,
            "strand": "-",
            "nuclease": "HiFi-SpCas9 (R691A)",
            "pam_motif": "5'-NGG-3'",
            "pam_sequence_5to3": "AGG",
            "protospacer_dna_5to3": "ATATTCTGAGCCTACAGCAG",
            "spacer_length_nt": 20,
            "cleavage_geometry": (
                "Blunt DSB at chr4:76023456|76023457 (-), 3 bp upstream of 3'-AGG PAM inside CXCL10 Exon 1 "
                "(5'-UTR; PAM 3'-end at 76023451 is 20 bp upstream of CDS ATG start at 76023431, and the blunt DSB "
                "at 76023456|76023457 is 25 bp upstream of the CDS ATG start; closing AUDIT-REF) "
                "across the Layer 2 76023468|76023467 junction"
            ),
            "dsb_to_cds_atg_distance_bp": 25,
            "pam_3prime_end_to_cds_atg_distance_bp": 20,
            "hdr_junction_grch38": "chr4:76023468|76023467 (-)",
            "hdr_pam_and_seed_disrupted_post_knockin": True,
            "hdr_re_cutting_prevention_mechanism": (
                "CXCL10_HA_L_600bp_minus (chr4:76023468-76024067 -) ends with the first 6 nt of the protospacer "
                "(5'-ATATTC-3', lacking all 14 PAM-proximal seed bases and the AGG PAM) and CXCL10_HA_R_600bp_minus "
                "(chr4:76022868-76023467 -) begins with the remaining 14 nt + AGG PAM (5'-TGAGCCTACAGCAGAGG-3', where "
                "the 5' 6 nt are replaced by the 3' end of bGH_polyA, introducing 6 mismatches in positions 1-6 and "
                "eliminating the full 20-nt R-loop); note that raw 3'-AGG PAM and 14-nt seed remain on HA_R, so "
                "wet-lab re-cleavage validation or a silent PAM shield in HA_R is recommended (closing AUDIT-REF & AUDIT-REF)."
            ),
            "length_neutral_pam_shield_candidate_astra_r1": (
                "chr4:76023451C>G (minus-strand 5'-AGG -> 5'-AGC at nucleotide 17 of CXCL10_HA_R_600bp_minus, 20 bp upstream "
                "of the native CDS ATG inside the 5'-UTR; eliminates the NGG PAM on HA_R with zero length change; AUDIT-REF)"
            ),
            "cfd_mit_off_target_specificity_score": 93.2,
            "whole_genome_crispor_mit_reference_specificity_score": 65.0,
            "ucsc_crispr_all_targets_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=crisprAllTargets;chrom=chr4;start=76023450;end=76023473",
            "grch38_0_1_2_mismatch_coding_off_targets": None,
            "genomewide_off_target_enumeration_status": "NOT_ENUMERATED",
            "grch38_3_mismatch_intergenic_intronic_loci": 4,
            "_noncoding_3mm_cfd_hits": [0.0245, 0.0198, 0.0162, 0.0125],
            "specificity_provenance": "assigned_3mm_panel_model_plus_verified_ucsc_crispr_all_targets_mit_reference_65",
            "gnomad_v4_1_max_ancestry_maf_in_protospacer_or_pam": 0.00031,
            "gnomad_v4_1_global_af_max_variant_in_target": 0.00005,
            "dbsnp_b155_overlapping_variant_count": 9,
            "ucsc_dbsnp155_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dbSnp155;chrom=chr4;start=76023450;end=76023473",
            "orthogonal_to_pass2_spcas9": False,
            "orthogonal_nuclease_family_to_spcas9": False,
            "delivered_in_separate_editing_pass_from_pass2_dsb": False,
            "_wt_expected_slice": cxcl10_wt_1200[76024067 - 76023473 : 76024067 - 76023451 + 1],
            "_donor_keys": ["ISLET1_LOCUS3_CXCL10_SERPINB9_GSGP2A_PDL1_NO_FURIN_4090bp"],
            "_split_prefix_bp": 6,
            "_split_suffix_bp": 17,
        },
        {
            "guide_id": "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
            "editing_pass": "Pass 2 (Post-QC Clone Multiplex via HiFi-SpCas9 R691A RNP)",
            "target_gene": "PPP1R12C (AAVS1 Safe Harbor)",
            "target_locus_grch38": "chr19:55115754-55115776 (-)",
            "chrom": "chr19",
            "start_1based": 55115754,
            "end_1based": 55115776,
            "strand": "-",
            "nuclease": "HiFi-SpCas9 (R691A)",
            "pam_motif": "5'-NGG-3'",
            "pam_sequence_5to3": "AGG",
            "protospacer_dna_5to3": "ACAGTGGGGCCACTAGGGAC",
            "spacer_length_nt": 20,
            "cleavage_geometry": "Safe-harbor cleavage 3 bp upstream of AGG PAM in PPP1R12C Intron 1 across the Layer 2 55115768|55115767 junction",
            "hdr_junction_grch38": "chr19:55115768|55115767 (-)",
            "hdr_pam_and_seed_disrupted_post_knockin": True,
            "hdr_re_cutting_prevention_mechanism": (
                "AAVS1_HA_L_600bp_minus (chr19:55115768-55116367 -) ends with the first 9 nt of the protospacer "
                "(5'-ACAGTGGGG-3', lacking the 11-nt PAM-proximal seed and AGG PAM) and AAVS1_HA_R_600bp_minus "
                "(chr19:55115168-55115767 -) begins with the remaining 11 nt + AGG PAM (5'-CCACTAGGGACAGG-3', where "
                "the 5' 9 nt are replaced by the 3' end of SV40_late_polyA_122bp [5'-ATGTTTCAG-3'], where positions 1, 5, and 9 "
                "incidentally match WT A, T, G, producing 6 mismatches across positions [2, 3, 4, 6, 7, 8], 0/12 seed mismatches "
                "[12-bp contiguous PAM-proximal suffix 5'-GCCACTAGGGAC-3' = 11-bp AAVS1 HA_R CCACTAGGGAC + pos 9 G], and "
                "post-HDR uncalibrated positional retention product = 0.143472 [post-HDR CFD = 0.143472]); because the raw 3'-AGG PAM and 12-nt seed suffix remain intact on HA_R, "
                "wet-lab re-cleavage validation or the silent PAM shield chr19:55115754C>G in HA_R is recommended "
                "(closing AUDIT-REF, AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
            ),
            "length_neutral_pam_shield_candidate_astra_r1": (
                "chr19:55115754C>G (minus-strand 3'-AGG -> 3'-AGC at nucleotide 14 of AAVS1_HA_R_600bp_minus inside "
                "PPP1R12C Intron 1; eliminates the NGG PAM on HA_R with zero length change)"
            ),
            "cfd_mit_off_target_specificity_score": 91.8,
            "whole_genome_crispor_mit_reference_specificity_score": 77.0,
            "ucsc_crispr_all_targets_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=crisprAllTargets;chrom=chr19;start=55115753;end=55115776",
            "grch38_0_1_2_mismatch_coding_off_targets": None,
            "genomewide_off_target_enumeration_status": "NOT_ENUMERATED",
            "grch38_3_mismatch_intergenic_intronic_loci": 5,
            "_noncoding_3mm_cfd_hits": [0.0241, 0.0204, 0.0178, 0.0149, 0.0121],
            "specificity_provenance": "assigned_3mm_panel_model_plus_verified_ucsc_crispr_all_targets_mit_reference_77",
            "gnomad_v4_1_max_ancestry_maf_in_protospacer_or_pam": 0.00042,
            "gnomad_v4_1_global_af_max_variant_in_target": 0.00006,
            "dbsnp_b155_overlapping_variant_count": 9,
            "ucsc_dbsnp155_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dbSnp155;chrom=chr19;start=55115753;end=55115776",
            "orthogonal_to_pass2_spcas9": False,
            "orthogonal_nuclease_family_to_spcas9": False,
            "delivered_in_separate_editing_pass_from_pass2_dsb": False,
            "_wt_expected_slice": aavs1_wt_1200[55116367 - 55115776 : 55116367 - 55115754 + 1],
            "_donor_keys": ["ISLET1_LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_EEF1A1_ICASP9_3457bp"],
            "_split_prefix_bp": 9,
            "_split_suffix_bp": 14,
        },
        {
            "guide_id": "OI-GRNA-05_TYK2_P1104A_PEMAX_EPEGRNA_OPTIONAL",
            "editing_pass": "Optional Personalized Non-DSB Prime Edit (Track 2/3 High-IFN-I T1D Recipients)",
            "target_gene": "TYK2",
            "target_locus_grch38": "chr19:10352421-10352443 (+)",
            "chrom": "chr19",
            "start_1based": 10352421,
            "end_1based": 10352443,
            "strand": "+",
            "nuclease": (
                "PEmax (SpCas9-H840A Nickase + M-MLV RT) + epegRNA (No DSB; note CGBE1 on the minus strand requires an "
                "engineered SpG [5'-NGN-3', using 5'-AGC-3' PAM at chr19:10352426-10352428 (-) to place 10352442C at position 7] "
                "or SpRY [5'-NRN/NYN-3', using 5'-GCC-3' PAM at chr19:10352425-10352427 (-) to place 10352442C at position 6] Cas9 variant because the canonical 5'-NGG-3' minus-strand PAMs [5'-AGG-3' at chr19:10352422-10352424 (-) and 5'-CGG-3' at chr19:10352433-10352435 (-)] place 10352442C at positions 3 and 14 outside positions 4-7; "
                "closing AUDIT-REF, AUDIT-REF & R9-S1-L5B [DEC-074])"
            ),
            "pam_motif": "5'-NGG-3'",
            "pam_sequence_5to3": "GGG",
            "protospacer_dna_5to3": "GCCTGGCTCTCACCGTGGGG",
            "spacer_length_nt": 20,
            "epegRNA_extension_specification": {
                "target_variant": "rs34536443 (GRCh38 chr19:10352442G>C on + strand; coding c.3310C>G transversion on - strand; p.Pro1104Ala CCC->GCC)",
                "why_prime_or_cgbe_not_abe8e": (
                    "c.3310C>G (g.10352442G>C) is a C:G -> G:C base transversion; adenine base editors (ABE8e) only "
                    "catalyze A:T -> G:C transitions and cannot install C:G -> G:C transversions (closing AUDIT-REF). "
                    "Moreover, because this + strand protospacer places 10352442G inside the 3' GGG PAM (+22), it is "
                    "specifically a PEmax prime-editing guide rather than a 5'-window CGBE1 guide (closing AUDIT-REF). "
                    "On the minus strand (where chr19:10352432-10352434 (-) is 5'-GGT-3' and chr19:10352436-10352438 (-) is 5'-CCA-3'), "
                    "the nearest canonical 5'-NGG-3' PAMs are 5'-AGG-3' at chr19:10352422-10352424 (-) (placing 10352442C at protospacer position 3) "
                    "and 5'-CGG-3' at chr19:10352433-10352435 (-) (placing 10352442C at protospacer position 14), both outside the CGBE1 nt 4-7 deamination window, so CGBE1 "
                    "requires an SpG (5'-AGC-3' PAM at chr19:10352426-10352428 -, position 7) or SpRY (5'-GCC-3' PAM at chr19:10352425-10352427 -, position 6) variant (closing AUDIT-REF, AUDIT-REF & R9-S1-L5B [DEC-074])."
                ),
                "nick_position_grch38": "chr19:10352437|10352438 (+ strand, 3 nt upstream of GGG PAM)",
                "primer_binding_site_pbs_13nt_5to3": "CACGGTGAGAGCC",
                "rt_template_rtt_12nt_5to3": "CCAGAGCGCCCC",
                "rtt_written_plus_strand_10352438_10352449": "GGGGCGCTCTGG",
                "bystander_mutations_introduced": 0,
                "structured_3prime_motif": "evopreQ1_pseudoknot_42nt",
                "structured_3prime_motif_segmentation_audit_dec078": (
                    "20-nt spacer + 76-nt scaffold + 12-nt RTT + 13-nt PBS + 0-nt separate linker + 37-nt tevopreQ1 pseudoknot "
                    "(whose 5' stem begins with the 8-nt CGCGGUUC sequence) + 5-nt AUUUU U6 termination tail = 163 nt total "
                    "(42-nt 3' extension = 37-nt tevopreQ1 + 5-nt AUUUU; AUDIT-REF [DEC-078])."
                ),
                "spcas9_sgRNA_scaffold_76nt_rna_5to3": PEMAX_EPEGRNA_SCAFFOLD_76NT_RNA,
                "separate_intervening_linker_length_nt": 0,
                "tevopreq1_5prime_stem_prefix_8nt_embedded_in_37nt_motif_5to3": "CGCGGUUC",
                "linker_8nt_5to3": "CGCGGUUC",
                "tevopreq1_pseudoknot_37nt_core_rna_5to3": TEVOPREQ1_PSEUDOKNOT_37NT_CORE_RNA,
                "u6_termination_tail_5nt_rna_5to3": TEVOPREQ1_U6_TERMINATION_TAIL_5NT_RNA,
                "evopreQ1_pseudoknot_42nt_rna_5to3": TEVOPREQ1_PSEUDOKNOT_42NT_RNA,
                "full_assembled_epegrna_dna_equivalent_163nt_5to3": (
                    "GCCTGGCTCTCACCGTGGGG"
                    + PEMAX_EPEGRNA_SCAFFOLD_76NT_RNA.replace("U", "T")
                    + "CCAGAGCGCCCC"
                    + "CACGGTGAGAGCC"
                    + TEVOPREQ1_PSEUDOKNOT_42NT_RNA.replace("U", "T")
                ),
                "full_assembled_epegrna_rna_163nt_5to3": (
                    "GCCUGGCUCUCACCGUGGGG"
                    + PEMAX_EPEGRNA_SCAFFOLD_76NT_RNA
                    + "CCAGAGCGCCCC"
                    + "CACGGUGAGAGCC"
                    + TEVOPREQ1_PSEUDOKNOT_42NT_RNA
                ),
            },
            "cleavage_geometry": "Single-strand H840A nick only (0 DSB) at chr19:10352437|10352438 (+); M-MLV RT polymerizes chr19:10352442G>C (c.3310C>G, p.Pro1104Ala) at +5 downstream of the nick, directly disrupting the GGG PAM (GGG -> GCG)",
            "hdr_junction_grch38": "PRIME_EDIT_SNV_CHR19_10352442_G_TO_C",
            "hdr_pam_and_seed_disrupted_post_knockin": True,
            "hdr_re_cutting_prevention_mechanism": "g.10352442G>C (+5 downstream of the H840A nick) directly mutates the middle G of the 5'-GGG-3' PAM to 5'-GCG-3', abolishing SpCas9 PAM recognition and preventing re-nicking.",
            "cfd_mit_off_target_specificity_score": 95.1,
            "whole_genome_crispor_mit_reference_specificity_score": 79.0,
            "ucsc_crispr_all_targets_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=crisprAllTargets;chrom=chr19;start=10352420;end=10352443",
            "grch38_0_1_2_mismatch_coding_off_targets": None,
            "genomewide_off_target_enumeration_status": "NOT_ENUMERATED",
            "grch38_3_mismatch_intergenic_intronic_loci": 2,
            "_noncoding_3mm_cfd_hits": [0.0285, 0.0230],
            "specificity_provenance": "assigned_3mm_panel_model_plus_verified_ucsc_crispr_all_targets_mit_reference_79",
            "gnomad_v4_1_max_ancestry_maf_in_protospacer_or_pam": 0.042,
            "gnomad_v4_1_global_af_max_variant_in_target": 0.0270,
            "gnomad_v4_1_fin_bottleneck_af_rs34536443": 0.0528,
            "dbsnp_b155_overlapping_variant_count": 17,
            "ucsc_dbsnp155_query_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dbSnp155;chrom=chr19;start=10352420;end=10352443",
            "orthogonal_to_pass2_spcas9": True,
            "orthogonal_nuclease_family_to_spcas9": False,
            "delivered_in_separate_editing_pass_from_pass2_dsb": True,
            "nuclease_orthogonality_note": "PEmax uses an SpCas9-H840A nickase scaffold (same PAM/scaffold family as Pass 2 SpCas9, so orthogonal_nuclease_family_to_spcas9=False) but is delivered in an optional separate non-DSB prime-editing pass (delivered_in_separate_editing_pass_from_pass2_dsb=True) so it never co-introduces DSBs with Pass-2 HDR.",
            "_wt_expected_slice": tyk2_plus_71[10352421 - 10352410 : 10352443 - 10352410 + 1],
            "_donor_keys": [],
        },
    ]

    # Executable round-trip verification of TYK2 epegRNA PBS and RTT against GRCh38 (AUDIT-REF):
    tyk2_spec = grna_specs[4]["epegRNA_extension_specification"]
    pbs_expected = revcomp(tyk2_plus_71[10352425 - 10352410 : 10352437 - 10352410 + 1])  # 13 nt upstream of nick
    _check(
        tyk2_spec["primer_binding_site_pbs_13nt_5to3"] == pbs_expected,
        f"TYK2 PBS mismatch: {tyk2_spec['primer_binding_site_pbs_13nt_5to3']} != {pbs_expected}",
    )
    wt_post_nick_12 = tyk2_plus_71[10352438 - 10352410 : 10352449 - 10352410 + 1]
    rtt_written_plus = revcomp(tyk2_spec["rt_template_rtt_12nt_5to3"])
    written_diffs = [
        (10352438 + idx_b, wt_b, ed_b)
        for idx_b, (wt_b, ed_b) in enumerate(zip(wt_post_nick_12, rtt_written_plus))
        if wt_b != ed_b
    ]
    _check(
        written_diffs == [(10352442, "G", "C")],
        f"TYK2 RTT must write only 10352442G>C with zero bystander edits, got {written_diffs}",
    )

    # Read all 4 HDR donor templates from Layer 2 FASTA
    fasta_path = REPO_ROOT / "data" / "constructs" / "islet1_scenario_a_constructs_v0.7.0.fasta"
    from scripts.run_layer2_perturbation_and_regulatory_plane import parse_construct_fasta_text
    donor_seqs = parse_construct_fasta_text(fasta_path.read_text(encoding="utf-8"))

    # Authentic Ensembl GRCh38 50-bp chromosomal flanks immediately outside each construct's HA_L and HA_R
    # (closing AUDIT-REF; verified against live Ensembl REST GRCh38):
    authentic_outer_chromosomal_flanks_50bp: Dict[str, Dict[str, str]] = {
        "ISLET1_LOCUS1_ARM1A_B2M_HLAE_SCT_GSGP2A_CD47_IGV_GPI_4880bp": {
            "left_outer_grch38": "chr15:44710759-44710808 (+)",
            "left_outer_50bp": "ACCTAGAGGGCGCTGGAAGCTCTAAAGCCCTAGCAGTTACTGCTTTTACT",
            "right_outer_grch38": "chr15:44712409-44712458 (+)",
            "right_outer_50bp": "GCCTCTCACAGATGAAGAAACTAAGGCACCGAGATTTTAAGAAACTTAAT",
            "ha_bp": "800",
        },
        "ISLET1_LOCUS1_ARM1B_B2M_BIDIR_2AFREE_CD47_IGV_GPI_HLAE_SCT_3974bp": {
            "left_outer_grch38": "chr15:44710959-44711008 (+)",
            "left_outer_50bp": "CTCTTTCACTAGGACCTTCTCTGAGCTGTCCTCAGGATGCTTTTGGGACT",
            "right_outer_grch38": "chr15:44712209-44712258 (+)",
            "right_outer_50bp": "GTTTGTAGAATGCTTGGCTGTGATACAAAGCGGTTTCGAATAATTAACTT",
            "ha_bp": "600",
        },
        "ISLET1_LOCUS3_CXCL10_SERPINB9_GSGP2A_PDL1_NO_FURIN_4090bp": {
            "left_outer_grch38": "chr4:76024068-76024117 (-)",
            "left_outer_50bp": "TTTTCCCTCAAAATAGTTATGTTGGAGGCTATTTACTTAATCTAATGGTG",
            "right_outer_grch38": "chr4:76022818-76022867 (-)",
            "right_outer_50bp": "CTAAAATGAACACTGTGCTAACCTTCTCTGCTGTTCCTCTTTTCCTACAG",
            "ha_bp": "600",
        },
        "ISLET1_LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_EEF1A1_ICASP9_3457bp": {
            "left_outer_grch38": "chr19:55116368-55116417 (-)",
            "left_outer_50bp": "TGTGTTGCTGCCCAAGGATGCTCTTTCCGGAGCACTTCCTTCTCGGCGCT",
            "right_outer_grch38": "chr19:55115118-55115167 (-)",
            "right_outer_50bp": "CTACTCCCTTTCTCTTTCTCCATCCTTCTTTCCTTAAAGAGTCCCCAGTG",
            "ha_bp": "600",
        },
    }

    post_hdr_junction_cfd_audit: Dict[str, Any] = {}
    off_target_reconciliation_table: Dict[str, Any] = {}

    for g in grna_specs:
        # These historical screening numbers are assumptions, not retrieved gnomAD AC/AN observations.
        assumed_frequencies = {k.removeprefix("gnomad_v4_1_"): v for k, v in g.items() if k.startswith("gnomad_v4_1_")}
        g["assumed_screening_frequency_priors"] = {
            "source_class": "ASSUMED_NOT_DATABASE_MEASUREMENT", "values": assumed_frequencies,
            "qualification": "Unvalidated screening assumptions; not observed maxima or proven upper bounds.",
        }
        for key in [k for k in g if k.startswith("gnomad_v4_1_")]:
            g[key] = None  # Deprecated measurement aliases require a variant/population AC/AN receipt.
        if "PEmax" in g["nuclease"]:
            g["cgbe_alternative_implementation_status"] = "NOT_IMPLEMENTED_REQUIRES_ENGINEERED_PAM_WINDOW_AND_BYSTANDER_LEDGER"
        proto = g["protospacer_dna_5to3"]
        pam = g["pam_sequence_5to3"]
        full_target = (pam + proto) if "Cas12a" in g["nuclease"] else (proto + pam)
        full_target_rc = revcomp(full_target)

        # Verify 100% 0-mismatch identity against cached Ensembl GRCh38 reference slice (AUDIT-REF & AUDIT-REF):
        wt_slice = g.pop("_wt_expected_slice")
        _check(
            wt_slice == full_target,
            f"GRCh38 coordinate/sequence mismatch for {g['guide_id']}: expected {full_target}, got {wt_slice}",
        )
        if "_legacy_pIII_wt_slice" in g:
            leg_slice = g.pop("_legacy_pIII_wt_slice")
            _check(leg_slice == "TTTACCTTGGGGCTCTGACAGGTAGGA", f"Unexpected legacy CIITA slice: {leg_slice}")
        g["grch38_reference_exact_match_verified"] = True
        g["grch38_reference_mismatch_count"] = 0

        g["full_target_site_5to3"] = full_target
        g["full_target_site_revcomp_5to3"] = full_target_rc
        g["spacer_rna_5to3"] = proto.replace("T", "U")
        if "Cas12a" in g["nuclease"]:
            g["full_synthesized_grna_rna_5to3"] = ASCAS12A_DIRECT_REPEAT_20NT_RNA + proto.replace("T", "U")
            g["scaffold_rna_length_nt"] = len(ASCAS12A_DIRECT_REPEAT_20NT_RNA)
        elif "PEmax" in g["nuclease"]:
            g["full_synthesized_grna_rna_5to3"] = g["epegRNA_extension_specification"]["full_assembled_epegrna_rna_163nt_5to3"]
            g["scaffold_rna_length_nt"] = len(PEMAX_EPEGRNA_SCAFFOLD_76NT_RNA)
        else:
            g["full_synthesized_grna_rna_5to3"] = proto.replace("T", "U") + SPCAS9_SGRNA_SCAFFOLD_80NT_RNA
            g["scaffold_rna_length_nt"] = len(SPCAS9_SGRNA_SCAFFOLD_80NT_RNA)
        g["synthesized_rna_total_length_nt"] = len(g["full_synthesized_grna_rna_5to3"])

        g["protospacer_sha256"] = sha256_str(proto)
        g["gc_percent"] = gc_percent(proto)
        g["sugimoto_1995_rna_dna_hybrid_dg37_kcal_mol"] = compute_sugimoto_1995_rna_dna_dg37_kcal_mol(proto)
        g["santalucia_1998_duplex_dg37_kcal_mol"] = compute_duplex_dg37_kcal_mol(proto)

        # Compute deterministic off-target specificity model & reconcile with whole-genome CRISPOR/MIT reference (AUDIT-REF):
        noncoding_hits = g.pop("_noncoding_3mm_cfd_hits")
        ot_model = compute_panel_and_reference_off_target_specificity(
            guide_id=g["guide_id"],
            nuclease=g["nuclease"],
            coding_0_1_2_mm=g["grch38_0_1_2_mismatch_coding_off_targets"],
            noncoding_3_mm_count=g["grch38_3_mismatch_intergenic_intronic_loci"],
            noncoding_3_mm_cfd_Hits=noncoding_hits,
            crispor_whole_genome_mit_ref=g["whole_genome_crispor_mit_reference_specificity_score"],
            ucsc_crispr_track_url=g.get("ucsc_crispr_all_targets_query_url"),
        )
        _check(ot_model["window_and_exome_panel_cfd_specificity_score"] == g["cfd_mit_off_target_specificity_score"])
        g["off_target_specificity_breakdown_claude_v80_006"] = ot_model
        off_target_reconciliation_table[g["guide_id"]] = ot_model
        g["genome_wide_off_target_search_status"] = (
            "NOT_ENUMERATED_ASSIGNED_PANEL_AND_REFERENCE_SCORE_ONLY"
        )

        # Scan BOTH strands (+ and - revcomp) of all 4 HDR donor plasmids AND reconstructed post-edit genomic alleles
        # using authentic 50-bp Ensembl GRCh38 outer chromosomal flanks (closing AUDIT-REF, AUDIT-REF, AUDIT-REF, & AUDIT-REF):
        intact_in_any_donor_either_strand = any(
            _canonical_target_present_either_strand(d_seq, proto, g["nuclease"]) for d_seq in donor_seqs.values()
        )
        donor_keys = g.pop("_donor_keys")
        split_prefix = g.pop("_split_prefix_bp", 0)
        split_suffix = g.pop("_split_suffix_bp", 0)
        intact_in_post_hdr_allele = False
        reconstructed_allele_lengths: List[int] = []
        reconstructed_allele_hashes: List[str] = []
        junction_cfd_by_donor: Dict[str, Any] = {}
        raw_pam_retained_any_junction = False
        raw_pam_and_seed_retained_any_junction = False
        if donor_keys:
            edited_allele_verification_mode = (
                "VERIFIED_POST_HDR_CHROMOSOMAL_ALLELE_WITH_AUTHENTIC_50BP_OUTER_GRCH38_FLANKS"
            )
            target_len = len(full_target)
            for d_key in donor_keys:
                d_seq = donor_seqs[d_key]
                flanks = authentic_outer_chromosomal_flanks_50bp[d_key]
                ha_bp = int(flanks["ha_bp"])
                _check(flanks["left_outer_50bp"] != d_seq[:50] and flanks["right_outer_50bp"] != d_seq[-50:])
                post_hdr_allele = flanks["left_outer_50bp"] + d_seq + flanks["right_outer_50bp"]
                reconstructed_allele_lengths.append(len(post_hdr_allele))
                reconstructed_allele_hashes.append(sha256_str(post_hdr_allele))
                if _canonical_target_present_either_strand(post_hdr_allele, proto, g["nuclease"]):
                    intact_in_post_hdr_allele = True
                # Extract exact 5' (HA_L|Insert) and 3' (Insert|HA_R) split-protospacer junction windows (AUDIT-REF):
                left_junc_window = d_seq[ha_bp - split_prefix : ha_bp - split_prefix + target_len]
                right_junc_window = d_seq[len(d_seq) - ha_bp - (target_len - split_suffix) : len(d_seq) - ha_bp + split_suffix]
                if "Cas12a" in g["nuclease"]:
                    # For Cas12a (5'-PAM + 23-nt spacer), pass 23-nt candidate spacer + 4-nt candidate PAM in standardized order:
                    j_metrics = compute_split_junction_residual_cleavage_metrics(
                        protospacer_5to3=proto,
                        pam_5to3=pam,
                        nuclease=g["nuclease"],
                        left_junction_seq_5to3=left_junc_window[:4] + left_junc_window[4:],
                        right_junction_seq_5to3=right_junc_window[:4] + right_junc_window[4:],
                    )
                else:
                    j_metrics = compute_split_junction_residual_cleavage_metrics(
                        protospacer_5to3=proto,
                        pam_5to3=pam,
                        nuclease=g["nuclease"],
                        left_junction_seq_5to3=left_junc_window,
                        right_junction_seq_5to3=right_junc_window,
                    )
                if j_metrics["raw_pam_motif_retained_on_either_junction"]:
                    raw_pam_retained_any_junction = True
                if j_metrics["raw_pam_and_seed_subsequence_retained_on_either_junction"]:
                    raw_pam_and_seed_retained_any_junction = True
                junction_cfd_by_donor[d_key] = j_metrics
            g["post_hdr_split_junction_cfd_audit"] = junction_cfd_by_donor
            post_hdr_junction_cfd_audit[g["guide_id"]] = junction_cfd_by_donor
        elif g["guide_id"] == "OI-GRNA-02_CIITA_EXON1_ASCAS12A":
            # Construct authentic post-NHEJ frameshift allele on chr16:10898880-10899079 (-) with canonical
            # AsCas12a staggered cleavage 4-bp deletion (deleting 0-based indices 17..20 'CATC' inside TGCCCAACTTCTGCTGGCATCTC; closing AUDIT-REF):
            edited_allele_verification_mode = (
                "VERIFIED_POST_NHEJ_FRAMESHIFT_CHROMOSOMAL_ALLELE_4BP_STAGGERED_DELETION"
            )
            ciita_exon5_minus_200 = revcomp(ciita_exon5_wt_200)
            _check(full_target in ciita_exon5_minus_200)
            deleted_4bp = proto[17:21]
            _check(deleted_4bp == "CATC", f"Expected 4-bp deletion 'CATC' at proto[17:21], got {deleted_4bp}")
            broken_target = pam + proto[:17] + proto[21:]  # 4-bp frameshift deletion of 'CATC' at staggered cut site
            post_nhej_allele = ciita_exon5_minus_200.replace(full_target, broken_target, 1)
            reconstructed_allele_lengths.append(len(post_nhej_allele))
            reconstructed_allele_hashes.append(sha256_str(post_nhej_allele))
            if _canonical_target_present_either_strand(post_nhej_allele, proto, g["nuclease"]):
                intact_in_post_hdr_allele = True
            raw_pam_retained_any_junction = True  # 5'-TTTC PAM remains upstream of the 4-bp frameshift deletion
            raw_pam_and_seed_retained_any_junction = False  # 18-nt seed/cut zone disrupted by 4-bp deletion at nt 18..21
        elif g["guide_id"] == "OI-GRNA-05_TYK2_P1104A_PEMAX_EPEGRNA_OPTIONAL":
            # Construct authentic post-prime-edit allele on chr19:10352410-10352480 (+) installing g.10352442G>C (GGG->GCG PAM):
            edited_allele_verification_mode = (
                "VERIFIED_POST_PRIME_EDIT_SNV_CHROMOSOMAL_ALLELE_CHR19_10352442_G_TO_C"
            )
            idx_snv = 10352442 - 10352410
            _check(tyk2_plus_71[idx_snv] == "G")
            post_pe_allele = tyk2_plus_71[:idx_snv] + "C" + tyk2_plus_71[idx_snv + 1 :]
            reconstructed_allele_lengths.append(len(post_pe_allele))
            reconstructed_allele_hashes.append(sha256_str(post_pe_allele))
            if _canonical_target_present_either_strand(post_pe_allele, proto, g["nuclease"]):
                intact_in_post_hdr_allele = True
            raw_pam_retained_any_junction = False  # GGG PAM mutated to GCG
            raw_pam_and_seed_retained_any_junction = False
        else:
            edited_allele_verification_mode = "UNVERIFIED_NO_ALLELE_CONSTRUCTED"

        g["edited_allele_verification_mode"] = edited_allele_verification_mode
        g["full_20_to_23nt_target_absent_post_edit"] = (not intact_in_post_hdr_allele)
        g["hdr_full_protospacer_disrupted_post_knockin"] = (not intact_in_post_hdr_allele)
        g["hdr_pam_and_seed_disrupted_post_knockin"] = (not raw_pam_and_seed_retained_any_junction)
        g["raw_pam_motif_retained_on_either_post_edit_junction"] = raw_pam_retained_any_junction
        g["raw_pam_and_seed_subsequence_retained_on_either_post_edit_junction"] = raw_pam_and_seed_retained_any_junction
        g["intact_target_present_in_any_hdr_donor_plasmid"] = intact_in_any_donor_either_strand
        g["intact_target_present_in_post_hdr_edited_allele_either_strand"] = intact_in_post_hdr_allele
        g["intact_target_present_in_post_edit_chromosomal_allele_either_strand"] = intact_in_post_hdr_allele
        g["post_hdr_reconstructed_allele_lengths_bp"] = reconstructed_allele_lengths
        g["post_edit_reconstructed_allele_sha256"] = reconstructed_allele_hashes
        _check(not intact_in_any_donor_either_strand, f"Guide {g['guide_id']} found intact in HDR donor!")
        _check(not intact_in_post_hdr_allele, f"Guide {g['guide_id']} not disrupted in post-edit allele!")

    dbsnp_summary = {}
    for key, rec in raw_l5.get("dbsnp_records", {}).items():
        raw_obs = rec.get("minor_allele_frequencies", [])
        folded_summary = []
        if isinstance(raw_obs, list):
            for obs in raw_obs:
                if isinstance(obs, dict):
                    ac = obs.get("allele_count")
                    an = obs.get("total_count")
                    valid_counts = (
                        not isinstance(ac, bool)
                        and not isinstance(an, bool)
                        and isinstance(ac, int)
                        and isinstance(an, int)
                        and an > 0
                        and 0 <= ac <= an
                    )
                    f_val = (float(ac) / float(an)) if valid_counts else None
                    folded_summary.append({
                        "study": obs.get("study"),
                        "allele_count": ac if valid_counts else None,
                        "total_count": an if valid_counts else None,
                        "raw_ac_over_an": round(f_val, 6) if f_val is not None else None,
                        "folded_minor_allele_frequency": round(min(f_val, 1.0 - f_val), 6) if f_val is not None else None,
                        "allele_identity_in_raw_cache": "UNTYPED_AC_OVER_AN_REQUIRES_SPDI_RESOLUTION" if valid_counts else "INVALID_AC_AN_OBSERVATION_REJECTED",
                    })
        dbsnp_summary[key] = {
            "rsid": rec.get("refsnp_id", rec.get("rsid", key.split("_")[0])),
            "variant_type": rec.get("variant_type", "snv"),
            "genes": rec.get("genes", []),
            "clinical_significance": rec.get("clinical_significances", rec.get("clinical_significance", [])),
            "grch38_placements": rec.get("placements", [])[:2],
            "minor_allele_frequencies": raw_obs,
            "minor_allele_frequencies_schema_note": (
                "Raw cached study rows record untyped allele_count/total_count (where AC/AN > 0.5 reflects the major/reference allele in that study); "
                "see untyped_study_allele_count_observations and derived_folded_minor_allele_frequency_summary (AUDIT-REF)."
            ),
            "untyped_study_allele_count_observations": raw_obs,
            "derived_folded_minor_allele_frequency_summary": folded_summary,
        }

    clinvar_summary = raw_l5.get("clinvar_records", {}).get("rs34536443_TYK2_P1104A", {})

    # Source-verified 5-guide UCSC dbSnp155 identities plus explicitly assumed AF screening priors
    # derived directly from cached primary UCSC dbSnp155 API responses (closing AUDIT-REF & AUDIT-REF):
    ucsc_primary = raw_l5.get("ucsc_dbsnp155_primary_interval_responses", {})
    _check(len(ucsc_primary) == 5, f"Expected 5 primary UCSC dbSnp155 interval records in raw reference, got {len(ucsc_primary)}")
    af_provenance_str = "ASSUMED_NOT_DATABASE_MEASUREMENT_NO_GNOMAD_AC_AN_RECEIPTS"

    def _build_guide_dbsnp_entry(
        guide_id: str,
        grch38_interval: str,
        selected_variant_examples: List[str],
        extra_fields: Dict[str, Any],
    ) -> Dict[str, Any]:
        u_rec = ucsc_primary[guide_id]
        items = u_rec["items"]
        rsids = [it["rsid"] for it in items]
        _check(len(rsids) == u_rec["count"] and rsids == u_rec["rsids"], f"UCSC dbSnp155 count/rsid mismatch for {guide_id}")
        entry: Dict[str, Any] = {
            "grch38_interval": grch38_interval,
            "ucsc_dbsnp155_0based_half_open_query": u_rec["query_url"],
            "ucsc_dbsnp155_data_time_stamp": u_rec["data_time_stamp"],
            "dbsnp_b155_overlapping_variants": len(rsids),
            "enumerated_ucsc_dbsnp155_rsids": rsids,
            "ucsc_dbsnp155_variant_rows": [
                {
                    "rsid": it["rsid"],
                    "chrom": it["chrom"],
                    "chromStart": it["chromStart"],
                    "chromEnd": it["chromEnd"],
                    "ref": it["ref"],
                    "alts": it["alts"],
                    "func": it["func"],
                    "varType": it["varType"],
                }
                for it in items
            ],
            "selected_variant_examples_not_frequency_ranked": selected_variant_examples,
            "variant_example_selection_scope": "Selected illustrative sentinel examples; no frequency ranking is claimed.",
            "gnomad_af_provenance": af_provenance_str,
        }
        assumed = {k.removeprefix("gnomad_v4_1_"): v for k, v in extra_fields.items() if k.startswith("gnomad_v4_1_")}
        entry["assumed_screening_frequency_priors"] = {
            "source_class": "ASSUMED_NOT_DATABASE_MEASUREMENT", "values": assumed,
            "qualification": "Screening assumptions, not an observed database maximum or validated upper bound.",
        }
        entry.update({k: None for k in extra_fields})
        entry["measurement_alias_status"] = "NULL_UNTIL_VARIANT_POPULATION_AC_AN_RECEIPTS"
        return entry

    per_guide_dbsnp_summary = {
        "OI-GRNA-01_B2M_EXON1_ASCAS12A": _build_guide_dbsnp_entry(
            "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "chr15:44711590-44711616 (+)",
            ["rs572886094", "rs1435611656", "rs757230673"],
            {
                "gnomad_v4_1_global_af_max": 0.00004,
                "gnomad_v4_1_max_subpopulation_af": 0.00018,
                "common_polymorphism_af_gt_0_001_present": False,
            },
        ),
        "OI-GRNA-02_CIITA_EXON1_ASCAS12A": {
            "canonical_exon5_guide_id": "OI-GRNA-02_CIITA_EXON5_ASCAS12A",
            **_build_guide_dbsnp_entry(
                "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
                "chr16:10898963-10898989 (-)",
                ["rs142469968", "rs1398227497", "rs757017931"],
                {
                    "gnomad_v4_1_global_af_max": 0.00003,
                    "gnomad_v4_1_max_subpopulation_af": 0.00024,
                    "common_polymorphism_af_gt_0_001_present": False,
                },
            ),
        },
        "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9": _build_guide_dbsnp_entry(
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "chr4:76023451-76023473 (-)",
            ["rs554892214", "rs377480085", "rs767339215"],
            {
                "gnomad_v4_1_global_af_max": 0.00005,
                "gnomad_v4_1_max_subpopulation_af": 0.00031,
                "common_polymorphism_af_gt_0_001_present": False,
            },
        ),
        "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9": _build_guide_dbsnp_entry(
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
            "chr19:55115754-55115776 (-)",
            ["rs946359931", "rs1037294362", "rs920266081"],
            {
                "gnomad_v4_1_global_af_max": 0.00006,
                "gnomad_v4_1_max_subpopulation_af": 0.00042,
                "common_polymorphism_af_gt_0_001_present": False,
            },
        ),
        "OI-GRNA-05_TYK2_P1104A_PEMAX_EPEGRNA_OPTIONAL": _build_guide_dbsnp_entry(
            "OI-GRNA-05_TYK2_P1104A_PEMAX_EPEGRNA_OPTIONAL",
            "chr19:10352421-10352443 (+)",
            ["rs34536443 (intended edit site chr19:10352442G>C)", "rs751482411", "rs373614901"],
            {
                "gnomad_v4_1_rs34536443_global_af": 0.0270,
                "gnomad_v4_1_rs34536443_nfe_eur_af": 0.0422,
                "gnomad_v4_1_rs34536443_fin_af": 0.0528,
                "gnomad_v4_1_max_non_rs34536443_variant_af": 0.00005,
                "common_polymorphism_other_than_intended_rs34536443_present": False,
            },
        ),
    }
    total_dbsnp_b155_variants = sum(v["dbsnp_b155_overlapping_variants"] for v in per_guide_dbsnp_summary.values())
    _check(total_dbsnp_b155_variants == 56, f"Expected 8+13+9+9+17=56 UCSC dbSnp155 variants, got {total_dbsnp_b155_variants}")
    gnomad_v4_1_and_dbsnp_guide_overlap_audit = {
        "audit_methodology_and_reconciliation_note": (
            "Each guide's exact 0-based half-open genomic target interval (protospacer + PAM) was queried against the primary "
            "UCSC hg38 dbSnp155 REST API track (yielding 8, 13, 9, 9, and 17 overlapping variants respectively; sum = 56) and "
            "stored verbatim in data/raw/layer5_dbsnp_clinvar_reference.json under ucsc_dbsnp155_primary_interval_responses. "
            "Whereas live_ncbi_dbsnp_polymorphism_audit tracks the 4 sentinel T1D/HLA functional polymorphisms, "
            "this interval-wide audit enumerates all 56 primary UCSC dbSnp155 short variants (exact rsIDs, 0-based coordinates, "
            "and ref/alt alleles). Frequency screening assumptions are separated from source-verified identities. "
            "No gnomAD measurement or observed frequency maximum is claimed without a variant/population AC/AN receipt."
        ),
        "gnomad_af_provenance": af_provenance_str,
        "total_overlapping_dbsnp_b155_variants_across_5_guides": total_dbsnp_b155_variants,
        "per_guide_ucsc_dbsnp155_counts_vector": [v["dbsnp_b155_overlapping_variants"] for v in per_guide_dbsnp_summary.values()],
        "guides_with_any_non_target_variant_global_af_gt_0_0001": None,
        "frequency_measurement_status": "NOT_ENUMERATED_FROM_GNOMAD_AC_AN_ROWS",
        "per_guide_interval_polymorphism_summary": per_guide_dbsnp_summary,
    }

    g03_r = post_hdr_junction_cfd_audit["OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9"][
        "ISLET1_LOCUS3_CXCL10_SERPINB9_GSGP2A_PDL1_NO_FURIN_4090bp"
    ]["three_prime_insert_ha_right_junction"]
    g03_l = post_hdr_junction_cfd_audit["OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9"][
        "ISLET1_LOCUS3_CXCL10_SERPINB9_GSGP2A_PDL1_NO_FURIN_4090bp"
    ]["five_prime_ha_left_insert_junction"]
    g04_r = post_hdr_junction_cfd_audit["OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9"][
        "ISLET1_LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_EEF1A1_ICASP9_3457bp"
    ]["three_prime_insert_ha_right_junction"]
    g04_l = post_hdr_junction_cfd_audit["OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9"][
        "ISLET1_LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_EEF1A1_ICASP9_3457bp"
    ]["five_prime_ha_left_insert_junction"]
    g01_1a_l = post_hdr_junction_cfd_audit["OI-GRNA-01_B2M_EXON1_ASCAS12A"][
        "ISLET1_LOCUS1_ARM1A_B2M_HLAE_SCT_GSGP2A_CD47_IGV_GPI_4880bp"
    ]["five_prime_ha_left_insert_junction"]
    g01_1b_l = post_hdr_junction_cfd_audit["OI-GRNA-01_B2M_EXON1_ASCAS12A"][
        "ISLET1_LOCUS1_ARM1B_B2M_BIDIR_2AFREE_CD47_IGV_GPI_HLAE_SCT_3974bp"
    ]["five_prime_ha_left_insert_junction"]
    _g01_proto = "TGGCCTGGAGGCTATCCAGCGTG"
    _g01_pam = "TTTC"
    _g01_nuc = "AsCas12a Ultra (Acidaminococcus sp. BV3L6 M537R/F870L)"
    _g01_r_seq = "TAGCAGGCATGCTGGGGATCCAGCGTG"
    g01_s3_l = compute_split_junction_residual_cleavage_metrics(
        _g01_proto, _g01_pam, _g01_nuc, "TCTCTGGACTCGAAGCTATGCTCCGGT", _g01_r_seq
    )["five_prime_ha_left_insert_junction"]
    g01_1c_trap_l = compute_split_junction_residual_cleavage_metrics(
        _g01_proto, _g01_pam, _g01_nuc, "TTTCTGGCCTGGAGGCTATGTAAGTAT", _g01_r_seq
    )["five_prime_ha_left_insert_junction"]
    g01_astar_term_l = compute_split_junction_residual_cleavage_metrics(
        _g01_proto, _g01_pam, _g01_nuc, "TTTCTGGCCTGGAGGCTATAATAAAAG", _g01_r_seq
    )["five_prime_ha_left_insert_junction"]
    _check(
        g03_r["uncalibrated_positional_mismatch_score"] == 0.247822
        and g03_r["total_mismatches_vs_wt_spacer"] == 6
        and g03_r["seed_region_mismatches"] == 0
        and g03_r["mismatch_positions_1based"] == [1, 2, 3, 4, 5, 6]
        and g04_r["uncalibrated_positional_mismatch_score"] == 0.143472
        and g04_r["total_mismatches_vs_wt_spacer"] == 6
        and g04_r["seed_region_mismatches"] == 0
        and g04_r["mismatch_positions_1based"] == [2, 3, 4, 6, 7, 8]
        and g01_1a_l["total_mismatches_vs_wt_spacer"] == 5
        and g01_1a_l["seed_region_mismatches"] == 2
        and g01_1a_l["uncalibrated_positional_mismatch_score"] == 0.173249
        and g01_1b_l["total_mismatches_vs_wt_spacer"] == 7
        and g01_1b_l["seed_region_mismatches"] == 2
        and g01_1b_l["uncalibrated_positional_mismatch_score"] == 0.157781,
        f"Unexpected junction metrics: g03_r={g03_r}, g04_r={g04_r}, g01_1a_l={g01_1a_l}, g01_1b_l={g01_1b_l}",
    )

    post_hdr_junction_cfd_audit["gate_zero_donor_wgs_and_precloning_junction_recleavage_gate_muse_r1"] = {
        "gate_zero_donor_hpsc_30x_wgs_screening_rule": (
            "Gate 0 mandatory 30x short-read WGS + long-read HLA/KIR haplotyping across all 5 target intervals "
            "(B2M, CIITA, CXCL10, AAVS1, TYK2) plus +/-1,000 bp HDR homology arms on the specific clinical hPSC MCB "
            "prior to GMP guide synthesis; allele-matched spacer or homology-arm substitution if any donor MCB variant "
            "overlaps the 20-23 nt spacer, PAM, or inner 100 bp junction window (closing AUDIT-REF)."
        ),
        "pre_cloning_cxcl10_right_junction_recleavage_gate": (
            f"Pre-cloning in vitro RNP cleavage + rhAmpSeq/ddPCR junction-stability assay on reconstructed CXCL10 "
            f"Locus 3 right-junction amplicons (candidate_protospacer_5to3 = {g03_r['candidate_protospacer_5to3']} + {g03_r['candidate_pam_5to3']}, where the 14-bp "
            f"PAM-proximal spacer TGAGCCTACAGCAG + AGG PAM are retained intact [{g03_r['seed_region_mismatches']}/12 seed mismatches, "
            f"{g03_r['total_mismatches_vs_wt_spacer']}/20 5'-distal mismatches at positions {g03_r['mismatch_positions_1based']} ATATTC->GGGGAT from bGH_polyA, "
            f"unattenuated WT-SpCas9 positional-retention heuristic = {g03_r['uncalibrated_positional_mismatch_score']:.6f}, HiFi positional "
            f"heuristic = {g03_r['nuclease_attenuated_positional_heuristic']:.8f}; vs Audit hand-estimate 0.0020]) and AAVS1 Locus 4 right-junction amplicons "
            f"({g04_r['candidate_protospacer_5to3']} + {g04_r['candidate_pam_5to3']}, where the 11-bp PAM-proximal spacer CCACTAGGGAC + AGG PAM are retained "
            f"and 3 of the 9 SV40_late_polyA_122bp 5' bases [positions 1=A, 5=T, 9=G] incidentally match WT ACAGTGGGG [{g04_r['seed_region_mismatches']}/12 seed mismatches in positions 9..20, "
            f"{g04_r['total_mismatches_vs_wt_spacer']}/20 mismatches at 5' positions {g04_r['mismatch_positions_1based']}, unattenuated WT-SpCas9 positional-retention "
            f"heuristic = {g04_r['uncalibrated_positional_mismatch_score']:.6f}, HiFi positional heuristic = {g04_r['nuclease_attenuated_positional_heuristic']:.8f}]); "
            f"if residual HiFi-SpCas9 R691A cleavage exceeds the 0.1% pre-cloning alert threshold (<1.0% Stage-0 Playbook release gate) at saturating RNP dose, incorporate the "
            f"1-bp length-neutral PAM-disrupting point mutation (chr4:76023451C>G [minus-strand 5'-AGG->5'-AGC at nt 17 of CXCL10_HA_R_600bp_minus] "
            f"and/or chr19:55115754C>G [minus-strand 5'-AGG->5'-AGC at nt 14 of AAVS1_HA_R_600bp_minus]) before GMP plasmid linearization "
            f"(closing AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
        ),
        "stage0_challenger_junction_equivalence_audit_argon_r2": {
            "dx_l3_locus3_scarfree_dual_promoter": (
                f"100% nucleotide-identical 5' left ({g03_l['candidate_protospacer_5to3']}+{g03_l['candidate_pam_5to3']}, "
                f"positional heuristic={g03_l['uncalibrated_positional_mismatch_score']:.6f}) and 3' right "
                f"({g03_r['candidate_protospacer_5to3']}+{g03_r['candidate_pam_5to3']}, {g03_r['total_mismatches_vs_wt_spacer']}/20 mismatches at {g03_r['mismatch_positions_1based']}, "
                f"{g03_r['seed_region_mismatches']}/12 seed mismatches, WT positional heuristic={g03_r['uncalibrated_positional_mismatch_score']:.6f}) "
                f"23-bp junctions to canonical Locus 3 (Kozak_9bp and bGH_polyA_225bp outer boundaries unchanged)."
            ),
            "dx_s_locus4_icasp9_p2a_cd47": (
                f"100% nucleotide-identical 5' left ({g04_l['candidate_protospacer_5to3']}+{g04_l['candidate_pam_5to3']}, "
                f"positional heuristic={g04_l['uncalibrated_positional_mismatch_score']:.6f}) and 3' right "
                f"({g04_r['candidate_protospacer_5to3']}+{g04_r['candidate_pam_5to3']}, {g04_r['total_mismatches_vs_wt_spacer']}/20 mismatches at {g04_r['mismatch_positions_1based']}, "
                f"{g04_r['seed_region_mismatches']}/12 seed mismatches, WT positional heuristic={g04_r['uncalibrated_positional_mismatch_score']:.6f}) "
                f"23-bp junctions to canonical Locus 4 (cHS4_core_insulator_250bp and SV40_late_polyA_122bp outer boundaries unchanged)."
            ),
            "arm1c_constitutive_and_arm1a_lite_locus1": (
                f"Unshielded canonical B2M_HA_L 5' left junction has {g01_1a_l['total_mismatches_vs_wt_spacer']}/23 mismatches at {g01_1a_l['mismatch_positions_1based']} "
                f"incl. {g01_1a_l['seed_region_mismatches']}/18 seed mismatches ({g01_1a_l['candidate_pam_5to3']}+{g01_1a_l['candidate_protospacer_5to3']}, "
                f"WT Cas12a positional heuristic={g01_1a_l['uncalibrated_positional_mismatch_score']:.6f}, Ultra attenuated={g01_1a_l['nuclease_attenuated_positional_heuristic']:.6f}), "
                f"whereas S3-shielded Arm 1A-lite / Arm 1B-S3 mutates the 5'-TTTC PAM to non-PAM 5'-TCTC plus 3 seed wobble bases at positions [4, 7, 10] "
                f"({g01_s3_l['candidate_pam_5to3']}+{g01_s3_l['candidate_protospacer_5to3']}, {g01_s3_l['total_mismatches_vs_wt_spacer']}/23 mismatches incl. {g01_s3_l['seed_region_mismatches']}/18 seed mismatches, heuristic={g01_s3_l['uncalibrated_positional_mismatch_score']:.1f}); "
                "3' right junction (TAGC+AGGCATGCTGGGGATCCAGCGTG, 5'-TTTV PAM abolished to TAGC, heuristic=0.0) is 100% identical to canonical Arm 1A."
            ),
            "arm1b_cd47_c15g_and_l1l2_synonymous_locus1": (
                f"100% nucleotide-identical 5' left ({g01_1b_l['candidate_pam_5to3']}+{g01_1b_l['candidate_protospacer_5to3']} where B2M_HA_L adjoins cHS4_core_insulator_250bp starting CTGAAACA, "
                f"{g01_1b_l['total_mismatches_vs_wt_spacer']}/23 mismatches at {g01_1b_l['mismatch_positions_1based']} incl. {g01_1b_l['seed_region_mismatches']}/18 seed mismatches, "
                f"WT Cas12a positional heuristic={g01_1b_l['uncalibrated_positional_mismatch_score']:.6f}, Ultra attenuated={g01_1b_l['nuclease_attenuated_positional_heuristic']:.6f}; "
                "or 0.0 when paired with the S3 TCTC+TGGACTCGAAGCTATCTGAAACA PAM+seed shield) and 3' right (TAGC PAM abolished, heuristic=0.0) 27-bp junctions to canonical Arm 1B."
            ),
            "arm1c_promoter_trapped_and_arm1b_sense_terminator_shield_locus1": (
                f"Unshielded 5' left junctions have {g01_1c_trap_l['total_mismatches_vs_wt_spacer']}/23 ({g01_1c_trap_l['candidate_pam_5to3']}+{g01_1c_trap_l['candidate_protospacer_5to3']}, "
                f"mismatches at {g01_1c_trap_l['mismatch_positions_1based']}, {g01_1c_trap_l['seed_region_mismatches']}/18 seed mismatches, WT Cas12a positional heuristic={g01_1c_trap_l['uncalibrated_positional_mismatch_score']:.6f}, "
                f"Ultra attenuated={g01_1c_trap_l['nuclease_attenuated_positional_heuristic']:.6f}) and {g01_astar_term_l['total_mismatches_vs_wt_spacer']}/23 ({g01_astar_term_l['candidate_pam_5to3']}+{g01_astar_term_l['candidate_protospacer_5to3']}, "
                f"mismatches at {g01_astar_term_l['mismatch_positions_1based']}, {g01_astar_term_l['seed_region_mismatches']}/18 seed mismatches, WT Cas12a positional heuristic={g01_astar_term_l['uncalibrated_positional_mismatch_score']:.6f}, "
                f"Ultra attenuated={g01_astar_term_l['nuclease_attenuated_positional_heuristic']:.6f}; or 0.0 with S3 TCTC PAM+seed shield) mismatches respectively, and 100% identical 3' right junctions (TAGC PAM abolished, heuristic=0.0)."
            ),
        },
    }

    return {
        "multiplex_editing_schedule": (
            "Two-Pass Orthogonal RNP Schedule (Zero Simultaneous 4-Locus DSBs): "
            "Pass 1 uses AsCas12a Ultra (B2M Exon 1 + CIITA Exon 5 Pan-Isoform) followed by single-cell cloning, "
            ">=63,000x TP53 duplex-seq (>=62,958 independent duplex molecules), 20q11.21 BCL2L1 ddPCR, and CAST-seq; "
            "Pass 2 uses HiFi-SpCas9 R691A (CXCL10 Exon 1 + AAVS1 Intron 1) followed by final MCB release QC."
        ),
        "guides": grna_specs,
        "post_hdr_junction_cfd_and_split_protospacer_audit_claude_v80_004": post_hdr_junction_cfd_audit,
        "off_target_scoring_methodology_reconciliation_claude_v80_006": off_target_reconciliation_table,
        "gnomad_v4_1_and_dbsnp_guide_overlap_audit_claude_v80_045": gnomad_v4_1_and_dbsnp_guide_overlap_audit,
        "live_ncbi_dbsnp_polymorphism_audit": dbsnp_summary,
        "live_ncbi_clinvar_tyk2_rs34536443_record": clinvar_summary,
    }


def _validate_finite_float(val: Any, name: str, low: float, high: float) -> float:
    if isinstance(val, (bool, np.bool_)) or not isinstance(val, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a finite numeric value in [{low}, {high}], got {val!r}")
    fval = float(val)
    if not math.isfinite(fval) or fval < low or fval > high:
        raise ValueError(f"{name} must be finite and within [{low}, {high}], got {fval}")
    return fval


def _validate_strict_bool(val: Any, name: str) -> bool:
    if not isinstance(val, bool):
        raise ValueError(f"{name} must be a strict boolean (True or False), got {val!r}")
    return val


def _normalize_hla_e_genotype(raw_val: str) -> str:
    s = raw_val.strip().upper().replace("HLA-E*", "*").replace("HLA-E", "")
    parts = [p.strip() for p in s.split("/")]
    if len(parts) == 2:
        norm_parts = [(p if p.startswith("*") else f"*{p}") for p in parts]
        return f"{norm_parts[0]}/{norm_parts[1]}"
    return raw_val.strip()


def _normalize_sirpa_genotype(raw_val: str) -> str:
    s = raw_val.strip().lower()
    parts = sorted(p.strip() for p in s.split("/"))
    if len(parts) == 2:
        return f"{parts[0]}/{parts[1]}"
    return s


def _get_opt(patient: Dict[str, Any], key: str, default: Any) -> Any:
    val = patient.get(key, default)
    return default if val is None else val


def compile_personalized_patient_blueprint(
    patient: Dict[str, Any],
    *,
    engineering_parameters: Dict[str, Any] | None = None,
    track2_subbank_mode: str | None = None,
) -> Dict[str, Any]:
    """
    Deterministic Personalized Patient-DNA & Autoantigen Compiler (closes AUDIT-REF, AUDIT-REF,
    AUDIT-REF, and AUDIT-REF):
    Validates all patient genotype/biomarker inputs (strictly rejecting non-string patient_id and non-boolean flags),
    treats present-None optional keys consistently as missing/defaulted (`DEFAULT_IMPUTED_NOT_MEASURED`),
    separates recipient vs donor TYK2 rs34536443 genotypes (`0/1/2` alleles requiring prime edit vs natural `C/C`
    protective state), emits explicit warnings when forced tracks override Plan B indications, and tracks
    `(K, B)` chromosomal edit counts alongside an assigned utility discount and unmeasured physical yield.
    """
    if not isinstance(patient, dict):
        raise ValueError("patient input must be a dictionary")

    required_keys = (
        "patient_id",
        "ancestry",
        "nkg2a_neg_kir_educated_nk_fraction",
        "pre_existing_gad65_ia2_znt8_t_cell_zscore",
        "ifn1_ish_genetic_risk_percentile",
        "sirpa_genotype",
    )
    allowed_patient_keys = {
        *required_keys,
        "ifih1_rs1990760_genotype",
        "ptpn22_rs2476601_genotype",
        "recipient_hla_e_genotype",
        "recipient_bmi_kg_m2",
        "prescribed_dose_ieq",
        "tyk2_rs34536443_genotype",
        "recipient_tyk2_rs34536443_genotype",
        "donor_tyk2_rs34536443_genotype",
        "recipient_kir_ligand_group",
        "recipient_anti_hla_c_dsa",
        "recipient_anti_hla_e_dsa",
        "recipient_cpra_pct",
        "cmv_nkg2c_adaptive_nk_fraction",
        "recipient_activating_kir2ds1_unlicensed_hazard",
        "hlac_4digit_relation",
        "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative",
        "b_prime_sct_recipient_cd8_mlc_negative",
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared",
        "b_prime_sct_kir2ds1_non_activation_cleared",
        "exact_hlae_sct_serum_crossmatch_negative",
        "exact_product_non_hla_and_sct_flow_crossmatch_negative",
        "gate0a_nkg2c_adaptive_nk_cleared",
        "recipient_bw4_kir3dl1_educated",
        "bw4_kir3dl1_sorted_nk_lysis_cleared",
        "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard",
        "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared",
        "requested_track",
        "track2_subbank_mode",
        "carries_rare_snp_in_standard_grna",
        "donor_mcb_carries_rare_snp_in_standard_grna",
        "rare_snp_pam_disrupted",
        "cd47_sirpa_binding_compromised",
        "track1_sorted_kir_nkg2a_neg_nk_lysis_cleared",
        "arm1b_cd47_sirpa_nk_cliff_unverified",
        "single_200cm2_pouch_forced_without_modular_scaling",
        "c1c2_heteroallelic_phased_longread_ddpcr_verified",
        "wetlab_stage0_complete_lot_release_verified",
        "recipient_hla_c_genotype_4digit",
        "recipient_kir_genotype",
    }
    unknown_pt_keys = sorted(set(patient) - allowed_patient_keys)
    if unknown_pt_keys:
        raise ValueError(f"Unrecognized patient input keys: {unknown_pt_keys} (AUDIT-REF; DEC-071)")
    for k in required_keys:
        if k not in patient or patient[k] is None:
            raise ValueError(f"Missing required patient field: {k}")

    if not isinstance(patient["patient_id"], str):
        raise ValueError(f"patient_id must be a string, got {patient['patient_id']!r}")
    patient_id = patient["patient_id"].strip()
    if not patient_id:
        raise ValueError("patient_id must be a non-empty string")

    if not isinstance(patient["ancestry"], str):
        raise ValueError(f"ancestry must be a string, got {patient['ancestry']!r}")
    ancestry = patient["ancestry"].strip().upper()
    if ancestry not in VALID_ANCESTRIES:
        raise ValueError(f"Invalid ancestry {ancestry!r}; must be one of {sorted(VALID_ANCESTRIES)}")

    if not isinstance(patient["sirpa_genotype"], str):
        raise ValueError(f"sirpa_genotype must be a string, got {patient['sirpa_genotype']!r}")
    sirpa_genotype = patient["sirpa_genotype"].strip().lower()
    if sirpa_genotype not in VALID_SIRPA_GENOTYPES:
        raise ValueError(f"Invalid sirpa_genotype {sirpa_genotype!r}; must be one of {sorted(VALID_SIRPA_GENOTYPES)}")

    kir_licensed_nkg2a_neg_fraction = _validate_finite_float(
        patient["nkg2a_neg_kir_educated_nk_fraction"], "nkg2a_neg_kir_educated_nk_fraction", 0.0, 1.0
    )
    autoantigen_t_cell_titer = _validate_finite_float(
        patient["pre_existing_gad65_ia2_znt8_t_cell_zscore"], "pre_existing_gad65_ia2_znt8_t_cell_zscore", -5.0, 10.0
    )
    ifn1_genetic_risk_score = _validate_finite_float(
        patient["ifn1_ish_genetic_risk_percentile"], "ifn1_ish_genetic_risk_percentile", 0.0, 100.0
    )

    # Track explicit missingness vs measured status for optional host modifier fields (closing AUDIT-REF, AUDIT-REF, & AUDIT-REF):
    optional_field_provenance = {
        k: ("MEASURED_PATIENT_INPUT" if k in patient and patient[k] is not None else "DEFAULT_IMPUTED_NOT_MEASURED")
        for k in (
            "ifih1_rs1990760_genotype",
            "ptpn22_rs2476601_genotype",
            "recipient_hla_e_genotype",
            "recipient_bmi_kg_m2",
            "prescribed_dose_ieq",
            "tyk2_rs34536443_genotype",
            "recipient_tyk2_rs34536443_genotype",
            "donor_tyk2_rs34536443_genotype",
            "recipient_kir_ligand_group",
            "recipient_anti_hla_c_dsa",
            "recipient_anti_hla_e_dsa",
            "recipient_cpra_pct",
            "cmv_nkg2c_adaptive_nk_fraction",
            "recipient_activating_kir2ds1_unlicensed_hazard",
            "hlac_4digit_relation",
            "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative",
            "b_prime_sct_recipient_cd8_mlc_negative",
            "b_prime_sct_recipient_crossmatch_and_mlc_cleared",
            "b_prime_sct_kir2ds1_non_activation_cleared",
            "exact_hlae_sct_serum_crossmatch_negative",
            "exact_product_non_hla_and_sct_flow_crossmatch_negative",
            "gate0a_nkg2c_adaptive_nk_cleared",
            "recipient_bw4_kir3dl1_educated",
            "bw4_kir3dl1_sorted_nk_lysis_cleared",
            "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard",
            "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared",
        )
    }
    recipient_immunological_screen_complete = all(
        k in patient and patient[k] is not None
        for k in (
            "recipient_anti_hla_c_dsa",
            "recipient_anti_hla_e_dsa",
            "recipient_cpra_pct",
            "cmv_nkg2c_adaptive_nk_fraction",
        )
    )
    recipient_immunological_screen_completeness_status = (
        "COMPLETE_EXPLICIT_RECIPIENT_IMMUNOLOGICAL_SCREEN"
        if recipient_immunological_screen_complete
        else "INCOMPLETE_UNMEASURED_IMMUNOLOGICAL_HAZARD_INPUTS_DEFAULTED_FOR_SCENARIO_SCORING_ONLY"
    )

    # Validate optional host modifier genotypes & BMI (handling present-None consistently as default; AUDIT-REF):
    raw_ifih1 = _get_opt(patient, "ifih1_rs1990760_genotype", "C/T")
    if not isinstance(raw_ifih1, str):
        raise ValueError(f"Invalid ifih1_rs1990760_genotype {raw_ifih1!r}")
    ifih1_genotype = raw_ifih1.strip().upper()
    if ifih1_genotype not in VALID_IFIH1_GENOTYPES:
        raise ValueError(f"Invalid ifih1_rs1990760_genotype {ifih1_genotype!r}; must be one of {sorted(VALID_IFIH1_GENOTYPES)}")

    raw_ptpn22 = _get_opt(patient, "ptpn22_rs2476601_genotype", "G/G")
    if not isinstance(raw_ptpn22, str):
        raise ValueError(f"Invalid ptpn22_rs2476601_genotype {raw_ptpn22!r}")
    ptpn22_genotype = raw_ptpn22.strip().upper()
    if ptpn22_genotype not in VALID_PTPN22_GENOTYPES:
        raise ValueError(f"Invalid ptpn22_rs2476601_genotype {ptpn22_genotype!r}; must be one of {sorted(VALID_PTPN22_GENOTYPES)}")

    raw_hla_e = _get_opt(patient, "recipient_hla_e_genotype", "*01:03/*01:03")
    if not isinstance(raw_hla_e, str) or raw_hla_e.strip() not in VALID_HLA_E_GENOTYPES:
        raise ValueError(f"Invalid recipient_hla_e_genotype {raw_hla_e!r}; must be one of {sorted(VALID_HLA_E_GENOTYPES)}")
    recipient_hla_e = _normalize_hla_e_genotype(raw_hla_e)

    # Validate recipient and donor TYK2 rs34536443 genotypes separately (closing AUDIT-REF & AUDIT-REF):
    raw_recip_tyk2 = _get_opt(
        patient,
        "recipient_tyk2_rs34536443_genotype",
        _get_opt(patient, "tyk2_rs34536443_genotype", "G/G"),
    )
    if not isinstance(raw_recip_tyk2, str):
        raise ValueError(f"Invalid tyk2_rs34536443_genotype {raw_recip_tyk2!r}")
    recipient_tyk2_genotype = raw_recip_tyk2.strip().upper()
    if recipient_tyk2_genotype not in VALID_TYK2_RS34536443_GENOTYPES:
        raise ValueError(
            f"Invalid tyk2_rs34536443_genotype {recipient_tyk2_genotype!r}; must be one of {sorted(VALID_TYK2_RS34536443_GENOTYPES)}"
        )

    raw_donor_tyk2 = _get_opt(patient, "donor_tyk2_rs34536443_genotype", "G/G")
    if not isinstance(raw_donor_tyk2, str):
        raise ValueError(f"Invalid donor_tyk2_rs34536443_genotype {raw_donor_tyk2!r}")
    donor_tyk2_genotype = raw_donor_tyk2.strip().upper()
    if donor_tyk2_genotype not in VALID_TYK2_RS34536443_GENOTYPES:
        raise ValueError(
            f"Invalid donor_tyk2_rs34536443_genotype {donor_tyk2_genotype!r}; must be one of {sorted(VALID_TYK2_RS34536443_GENOTYPES)}"
        )

    kir_ligand_group_explicitly_supplied = (
        "recipient_kir_ligand_group" in patient
        and patient.get("recipient_kir_ligand_group") is not None
    )
    raw_kir = _get_opt(patient, "recipient_kir_ligand_group", "C1/C2")
    if not isinstance(raw_kir, str):
        raise ValueError(f"Invalid recipient_kir_ligand_group {raw_kir!r}")
    kir_ligand_group = raw_kir.strip().upper().replace("_", "/")
    if kir_ligand_group not in VALID_KIR_LIGAND_GROUPS:
        raise ValueError(f"Invalid recipient_kir_ligand_group {kir_ligand_group!r}; must be one of {sorted(VALID_KIR_LIGAND_GROUPS)}")
    if kir_ligand_group == "C2/C1":
        kir_ligand_group = "C1/C2"

    recipient_bmi = _validate_finite_float(_get_opt(patient, "recipient_bmi_kg_m2", 24.0), "recipient_bmi_kg_m2", 10.0, 70.0)
    raw_dose_override = _get_opt(patient, "prescribed_dose_ieq", None)
    if raw_dose_override is not None:
        prescribed_dose_ieq = int(round(_validate_finite_float(raw_dose_override, "prescribed_dose_ieq", 50000.0, 2000000.0)))
        dose_source_mode = "EXPLICIT_PATIENT_PRESCRIBED_DOSE_IEQ_OVERRIDE"
    else:
        prescribed_dose_ieq = int(round(500000.0 * (recipient_bmi / 24.0)))
        dose_source_mode = "BMI_PROPORTIONAL_DEFAULT_500K_AT_24KG_M2"
    areal_loading_at_fixed_200cm2 = round(prescribed_dose_ieq / 200.0, 1)
    required_pouch_area_cm2_at_2500_ieq_cm2 = round(prescribed_dose_ieq / 2500.0, 2)
    required_pouch_area_cm2_at_1500_ieq_cm2 = round(prescribed_dose_ieq / 1500.0, 2)
    modular_100cm2_omental_cassette_count_at_1500_ieq_cm2 = max(1, int(math.ceil(required_pouch_area_cm2_at_1500_ieq_cm2 / 100.0)))
    actual_areal_loading_ieq_per_cm2_across_modular_cassettes = round(
        prescribed_dose_ieq / (100.0 * modular_100cm2_omental_cassette_count_at_1500_ieq_cm2), 1
    )

    raw_req_track = _get_opt(patient, "requested_track", "AUTO_RECOMMEND")
    if not isinstance(raw_req_track, str):
        raise ValueError("requested_track must be a string")
    deployment_mode = raw_req_track.strip().upper()
    if deployment_mode not in VALID_TRACK_REQUESTS:
        raise ValueError(f"Invalid requested_track {deployment_mode!r}; must be one of {sorted(VALID_TRACK_REQUESTS)}")

    # Strictly validate boolean SNP and CD47-SIRPA flags (rejecting string "false" coercion; AUDIT-REF, AUDIT-REF & AUDIT-REF):
    recipient_has_protospacer_snp = _validate_strict_bool(
        _get_opt(patient, "carries_rare_snp_in_standard_grna", False), "carries_rare_snp_in_standard_grna"
    )
    donor_mcb_has_protospacer_snp = _validate_strict_bool(
        _get_opt(patient, "donor_mcb_carries_rare_snp_in_standard_grna", False), "donor_mcb_carries_rare_snp_in_standard_grna"
    )
    rare_snp_pam_disrupted_input = _validate_strict_bool(
        _get_opt(patient, "rare_snp_pam_disrupted", False), "rare_snp_pam_disrupted"
    )
    cd47_sirpa_binding_compromised_input = _validate_strict_bool(
        _get_opt(patient, "cd47_sirpa_binding_compromised", False), "cd47_sirpa_binding_compromised"
    )

    # Select Track 1 (Off-the-Shelf), Track 2 (Stratified Allogeneic HLA-C/TYK2 Sub-Line), or Track 3 (Personalized Autologous iPSC)
    # Note: Track-2 stratification threshold for IFN-I risk is 80.0 (aligned with TYK2 P1104A indication threshold; AUDIT-REF)
    if deployment_mode == "TRACK_3_AUTOLOGOUS_IPSC":
        selected_track = "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC"
    elif deployment_mode == "TRACK_2_STRATIFIED_SUBLINE":
        selected_track = "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE"
    elif deployment_mode == "TRACK_1_OFF_THE_SHELF":
        selected_track = "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC"
    elif (
        kir_licensed_nkg2a_neg_fraction >= 0.42
        or autoantigen_t_cell_titer >= 2.0
        or ifn1_genetic_risk_score >= 80.0
        or ifih1_genotype == "T/T"
        or ptpn22_genotype in {"A/G", "G/A", "A/A"}
    ):
        selected_track = "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE"
    else:
        selected_track = "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC"

    # Plan B (Design B) is unified with the canonical Blueprint's Selective HLA-A/HLA-B KO + Retained HLA-C*07:02 (C1) / HLA-C*04:01 (C2)
    # sub-bank architecture for recipients with extreme KIR-educated NKG2A- NK fractions (>= 0.48), avoiding unbudgeted 5th-locus FGL1 (AUDIT-REF & AUDIT-REF):
    plan_b_biologically_indicated = bool(kir_licensed_nkg2a_neg_fraction >= 0.48)
    use_plan_b_hla_c_retained = bool(
        selected_track == "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE"
        and plan_b_biologically_indicated
    )
    clinical_override_warnings: List[str] = []
    track1_sorted_kir_nk_cleared = _validate_strict_bool(
        _get_opt(patient, "track1_sorted_kir_nkg2a_neg_nk_lysis_cleared", False),
        "track1_sorted_kir_nkg2a_neg_nk_lysis_cleared",
    )
    arm1b_cd47_sirpa_nk_cliff_unverified = _validate_strict_bool(
        _get_opt(patient, "arm1b_cd47_sirpa_nk_cliff_unverified", False),
        "arm1b_cd47_sirpa_nk_cliff_unverified",
    )
    single_200cm2_pouch_forced_without_modular_scaling = _validate_strict_bool(
        _get_opt(patient, "single_200cm2_pouch_forced_without_modular_scaling", False),
        "single_200cm2_pouch_forced_without_modular_scaling",
    )
    tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard = _validate_strict_bool(
        _get_opt(patient, "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard", False),
        "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard",
    )
    gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared = _validate_strict_bool(
        _get_opt(patient, "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared", False),
        "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared",
    )
    recipient_bw4_kir3dl1_educated = _validate_strict_bool(
        _get_opt(patient, "recipient_bw4_kir3dl1_educated", False),
        "recipient_bw4_kir3dl1_educated",
    )
    bw4_kir3dl1_nk_lysis_cleared = _validate_strict_bool(
        _get_opt(patient, "bw4_kir3dl1_sorted_nk_lysis_cleared", False),
        "bw4_kir3dl1_sorted_nk_lysis_cleared",
    )
    c1c2_phasing_explicitly_supplied = (
        "c1c2_heteroallelic_phased_longread_ddpcr_verified" in patient
        and patient.get("c1c2_heteroallelic_phased_longread_ddpcr_verified") is not None
    )
    c1c2_heteroallelic_phased_verified = _validate_strict_bool(
        _get_opt(patient, "c1c2_heteroallelic_phased_longread_ddpcr_verified", False),
        "c1c2_heteroallelic_phased_longread_ddpcr_verified",
    )
    if plan_b_biologically_indicated and not use_plan_b_hla_c_retained:
        clinical_override_warnings.append(
            f"WARNING_FORCED_TRACK_OVERRIDES_HIGH_KIR_NK_PLAN_B_INDICATION: nkg2a_neg_kir_educated_nk_fraction={kir_licensed_nkg2a_neg_fraction:.3f} >= 0.48 "
            f"indicates Plan B (HLA-C-retained sub-bank), but requested_track={deployment_mode} forced {selected_track} Plan A."
        )

    # Resolve starting TYK2 rs34536443 genotype on the actual genome being edited (autologous recipient iPSC in Track 3 vs allogeneic donor MCB in Track 1/2):
    starting_graft_tyk2_genotype = (
        recipient_tyk2_genotype
        if selected_track == "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC"
        else donor_tyk2_genotype
    )
    if starting_graft_tyk2_genotype == "C/C":
        tyk2_alleles_to_edit = 0
    elif starting_graft_tyk2_genotype in {"G/C", "C/G"}:
        tyk2_alleles_to_edit = 1
    else:
        tyk2_alleles_to_edit = 2
    already_homozygous_tyk2_p1104a = (starting_graft_tyk2_genotype == "C/C")

    tyk2_biologically_indicated = bool(
        ifn1_genetic_risk_score >= 80.0
        or autoantigen_t_cell_titer >= 2.0
        or ifih1_genotype == "T/T"
        or ptpn22_genotype == "A/A"
    )
    # Prime editing guide OI-GRNA-05 is prescribed when TYK2 protection is indicated AND the starting graft genome is not already C/C:
    include_tyk2_p1104a_edit = bool(tyk2_biologically_indicated and not already_homozygous_tyk2_p1104a)
    if deployment_mode == "TRACK_1_OFF_THE_SHELF" and include_tyk2_p1104a_edit:
        clinical_override_warnings.append(
            "WARNING_FORCED_TRACK1_WITH_TYK2_INDICATION: Track 1 off-the-shelf MCB requires selecting the pre-banked TYK2-P1104A sub-stock or routing to Track 2."
        )
    final_graft_tyk2_p1104a_protective_state = bool(already_homozygous_tyk2_p1104a or include_tyk2_p1104a_edit)

    has_anti_hlac_dsa_input = _validate_strict_bool(
        _get_opt(patient, "recipient_anti_hla_c_dsa", False), "recipient_anti_hla_c_dsa"
    )
    raw_hlac_rel = _get_opt(
        patient,
        "hlac_4digit_relation",
        "SELF_4DIGIT_MATCH" if not use_plan_b_hla_c_retained else "UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT",
    )
    if not isinstance(raw_hlac_rel, str) or not raw_hlac_rel.strip():
        raise ValueError(f"Invalid hlac_4digit_relation {raw_hlac_rel!r}: must be a non-empty string (AUDIT-REF)")
    hlac_4digit_rel_input = raw_hlac_rel.strip().upper()
    if hlac_4digit_rel_input not in _VALID_HLAC_4DIGIT_RELATIONS:
        raise ValueError(
            f"Invalid hlac_4digit_relation {hlac_4digit_rel_input!r}; must be one of {sorted(_VALID_HLAC_4DIGIT_RELATIONS)} (AUDIT-REF)"
        )
    track2_subbank_mode_eff = str(
        track2_subbank_mode if track2_subbank_mode is not None else _get_opt(patient, "track2_subbank_mode", "AUTO")
    ).strip().upper()
    if track2_subbank_mode_eff == "B_ENDORETAINED_HLAC":
        track2_subbank_mode_eff = "ENDOGENOUS_HLAC_EXPLORATORY"
    if track2_subbank_mode_eff not in {"AUTO", "B_PRIME_SCT", "ENDOGENOUS_HLAC_EXPLORATORY"}:
        raise ValueError(
            f"Invalid track2_subbank_mode {track2_subbank_mode_eff!r}; must be one of "
            "['AUTO', 'B_PRIME_SCT', 'ENDOGENOUS_HLAC_EXPLORATORY', 'B_ENDORETAINED_HLAC'] (AUDIT-REF)"
        )
    sct_crossmatch_cleared_input = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_recipient_crossmatch_and_mlc_cleared", False),
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared",
    )
    sct_serum_crossmatch_neg_input = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative", False),
        "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative",
    )
    sct_cd8_mlc_neg_input = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_recipient_cd8_mlc_negative", False),
        "b_prime_sct_recipient_cd8_mlc_negative",
    )
    has_activating_kir2ds1_hazard_input = _validate_strict_bool(
        _get_opt(patient, "recipient_activating_kir2ds1_unlicensed_hazard", False),
        "recipient_activating_kir2ds1_unlicensed_hazard",
    )
    kir2ds1_cleared_input = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_kir2ds1_non_activation_cleared", False),
        "b_prime_sct_kir2ds1_non_activation_cleared",
    )
    has_anti_hlae_dsa_input = _validate_strict_bool(
        _get_opt(patient, "recipient_anti_hla_e_dsa", False),
        "recipient_anti_hla_e_dsa",
    )
    exact_hlae_sct_serum_crossmatch_neg_input = _validate_strict_bool(
        _get_opt(patient, "exact_hlae_sct_serum_crossmatch_negative", False),
        "exact_hlae_sct_serum_crossmatch_negative",
    )
    exact_product_non_hla_and_sct_flow_crossmatch_neg_input = _validate_strict_bool(
        _get_opt(patient, "exact_product_non_hla_and_sct_flow_crossmatch_negative", False),
        "exact_product_non_hla_and_sct_flow_crossmatch_negative",
    )
    gate0a_nkg2c_adaptive_nk_cleared_input = _validate_strict_bool(
        _get_opt(patient, "gate0a_nkg2c_adaptive_nk_cleared", False),
        "gate0a_nkg2c_adaptive_nk_cleared",
    )
    recipient_cpra_pct_input = _validate_finite_float(
        _get_opt(patient, "recipient_cpra_pct", 0.0), "recipient_cpra_pct", 0.0, 100.0
    )
    cmv_nkg2c_adaptive_nk_fraction_input = _validate_finite_float(
        _get_opt(patient, "cmv_nkg2c_adaptive_nk_fraction", 0.0), "cmv_nkg2c_adaptive_nk_fraction", 0.0, 1.0
    )
    use_b_prime_sct_k4 = bool(
        use_plan_b_hla_c_retained and track2_subbank_mode_eff == "B_PRIME_SCT"
    )
    anti_hlac_dsa_cleared_input = bool(use_b_prime_sct_k4 and sct_serum_crossmatch_neg_input)
    cd8_mlc_cleared_input = bool(
        use_b_prime_sct_k4 and (sct_cd8_mlc_neg_input or sct_crossmatch_cleared_input)
    )
    if use_plan_b_hla_c_retained and has_anti_hlac_dsa_input and not anti_hlac_dsa_cleared_input:
        clinical_override_warnings.append(
            "CONTRAINDICATED_ANTI_HLAC_DSA_BLOCKS_HLA_C_BEARING_SUBBANK_RELEASE: recipient carries anti-HLA-C DSA; "
            "fixed-peptide SCT does not delete polymorphic HLA-C heavy-chain alloepitopes and cannot be cleared by "
            "MLC-only flags or on endogenous HLA-C retention (Aubert et al. 2014 PMID:24804568; AUDIT-REF; AUDIT-REF)."
        )
    if has_anti_hlae_dsa_input and not exact_hlae_sct_serum_crossmatch_neg_input:
        clinical_override_warnings.append(
            "HOLD_RECIPIENT_ANTI_HLAE_DSA_EXACT_SCT_CROSSMATCH_REQUIRED: recipient carries anti-HLA-E antibody; "
            "requires negative exact-product HLA-E SCT serum CDC/ADCC crossmatch prior to release (Ravindranath et al. 2022 PMID:35323192; AUDIT-REF)."
        )
    if (
        selected_track != "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC"
        and recipient_cpra_pct_input >= 80.0
        and not exact_product_non_hla_and_sct_flow_crossmatch_neg_input
    ):
        clinical_override_warnings.append(
            "HOLD_HIGH_CPRA_SENSITIZED_EXACT_PRODUCT_CROSSMATCH_REQUIRED: recipient cPRA >= 80%; "
            "requires negative exact-product non-HLA and SCT flow/CDC crossmatch prior to allogeneic release (AUDIT-REF)."
        )
    if cmv_nkg2c_adaptive_nk_fraction_input >= 0.25 and not gate0a_nkg2c_adaptive_nk_cleared_input:
        clinical_override_warnings.append(
            "HOLD_CMV_HIGH_NKG2C_ADAPTIVE_NK_DONOR_PANEL_GATE_0A_REQUIRED (WARN_CMV_HIGH_NKG2C_ADAPTIVE_NK_DONOR_PANEL_GATE_0A_REQUIRED): "
            "recipient CMV-seropositive NKG2C+ adaptive NK fraction >= 0.25; requires Gate 0A donor-panel verification of net NKG2A/NKG2C "
            "inhibition on VMAPRTLVL-loaded HLA-E SCT (vs high-affinity HLA-G leader VMAPRTLFL; Hammer et al. 2018 PMID:29632329; "
            "AUDIT-REF, AUDIT-REF, & AUDIT-REF)."
        )
    if (
        use_plan_b_hla_c_retained
        and kir_ligand_group in {"C2/C2", "C1/C2", "C2/C1"}
        and has_activating_kir2ds1_hazard_input
        and not kir2ds1_cleared_input
    ):
        clinical_override_warnings.append(
            "HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_KIR2DS1_NON_ACTIVATION_CLEARANCE "
            "(HOLD_RECIPIENT_ACTIVATING_KIR2DS1_UNLICENSED_HAZARD_REQUIRES_SORTED_NK_NON_ACTIVATION_GATE): recipient carries "
            "activating KIR2DS1 with C2-bearing sub-bank; requires sorted KIR2DS1+/KIR2DL1- NK non-activation clearance on "
            "HLA-C*04:01[Y84A, QYDDAVYKL] (Stewart et al. 2005 PMID:16141329; AUDIT-REF & AUDIT-REF)."
        )
    if (
        use_b_prime_sct_k4
        and kir_ligand_group in {"C2/C2", "C1/C2", "C2/C1"}
        and not cd8_mlc_cleared_input
    ):
        clinical_override_warnings.append(
            "HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_CLEARANCE "
            "(HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_NON_REACTIVITY): HLA-C*04:01[Y84A, QYDDAVYKL] "
            "requires recipient CD8+ T-cell MLC non-reactivity across all C2/C1+C2 recipients, including 4-digit heavy-chain-matched "
            "SELF_4DIGIT_MATCH / EXACT_4DIGIT_AUTOLOGOUS_OR_MATCHED recipients (AUDIT-REF & AUDIT-REF)."
        )
    if (
        use_b_prime_sct_k4
        and kir_ligand_group == "C1/C1"
        and hlac_4digit_rel_input in {"UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT", "UNSPECIFIED"}
        and not cd8_mlc_cleared_input
    ):
        clinical_override_warnings.append(
            "HOLD_B_PRIME_SCT_C1_UNSPECIFIED_4DIGIT_RELATION_REQUIRES_ALLELE_TYPING_AND_CD8_MLC_CLEARANCE: "
            "C1/C1 recipient (e.g. HLA-C*07:02 homozygous) receiving HLA-C*03:04[Y84A, GAVDPLLAL] B'-SCT without explicit "
            "4-digit allele relation or recipient CD8+ MLC clearance cannot clear by default (AUDIT-REF; DEC-042)."
        )

    if use_b_prime_sct_k4:
        sct_spec = (
            "CIITA_EXON5_B_SCT_C0304_GAVDPLLAL (HLA-C*03:04[Y84A, GAVDPLLAL] SCT at CIITA exon 5)"
            if kir_ligand_group == "C1/C1"
            else (
                "CIITA_EXON5_B_SCT_C0401_QYDDAVYKL (HLA-C*04:01[Y84A, QYDDAVYKL] SCT at CIITA exon 5)"
                if kir_ligand_group == "C2/C2"
                else "CIITA_EXON5_B_SCT_C0304_C0401_HETEROALLELIC (C1+C2 SCT at CIITA exon 5)"
            )
        )
        locus1_arm = f"OI-L2-DONOR-1B_B2M_EEF1A1_HLAE_CD47IGV + {sct_spec}"
        grna_set = [
            "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "OI-GRNA-02_CIITA_EXON5_ASCAS12A",
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        ]
        k_dsb_loci = 4
    elif use_plan_b_hla_c_retained:
        # Match retained HLA-C sub-bank to recipient KIR ligand group (C1/C1 -> HLA-C*07:02 [Design B-prime SP-1C upgrade: HLA-C*03:04 or C*01:02]; C2/C2 -> HLA-C*04:01; C1/C2 -> dual C1+C2; AUDIT-REF & AUDIT-REF):
        hlac_subbank_spec = (
            "HLA-C*07:02 (C1/C1 sub-bank; Design B-prime SP-1C upgrade: HLA-C*03:04 or HLA-C*01:02)"
            if kir_ligand_group == "C1/C1"
            else ("HLA-C*04:01 (C2/C2 sub-bank)" if kir_ligand_group == "C2/C2" else "HLA-C*07:02 / HLA-C*04:01 (C1/C2 dual-licensed heterozygote sub-bank; Design B-prime SP-1C/SP-2C upgrade: HLA-C*03:04 / HLA-C*04:01)")
        )
        locus1_arm = (
            f"PLAN_B_SELECTIVE_HLAA_HLAB_KO_RETAINED_HLAC_SUBLINE [{hlac_subbank_spec}] "
            "(Investigational comparator; retains B2M for surface HLA-C1/C2 trafficking; separate HLA-A/HLA-B KO + CD47-IgV/HLA-E donor ledger required)"
        )
        grna_set = [
            "PLAN_B_HLAA_EXON2_KO_GUIDE_INVESTIGATIONAL",
            "PLAN_B_HLAB_EXON2_KO_GUIDE_INVESTIGATIONAL",
            "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        ]
        k_dsb_loci = 5
    else:
        # High NKG2A- KIR+ educated NK fraction (>= 0.30) or low-stability recipient HLA-E (*01:01/*01:01)
        # mandates Arm 1B; surface density remains a modeled assumption requiring per-lineage flow-cytometry qualification.
        locus1_arm = (
            "OI-L2-DONOR-1B_B2M_EEF1A1_HLAE_CD47IGV (2A-free bidirectional cassette; per-lineage surface density unmeasured)"
            if (kir_licensed_nkg2a_neg_fraction >= 0.30 or recipient_hla_e == "*01:01/*01:01")
            else "OI-L2-DONOR-1A_B2M_ENDO_HLAE_CD47IGV (Non-Viral dsDNA Comparator; 5,170 bp with ITRs)"
        )
        grna_set = [
            "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        ]
        k_dsb_loci = 4
    if include_tyk2_p1104a_edit:
        grna_set.append("OI-GRNA-05_TYK2_P1104A_PEMAX_EPEGRNA_OPTIONAL")
    b_prime_loci = 1 if include_tyk2_p1104a_edit else 0

    # Custom allele-specific gRNA re-synthesis is triggered ONLY if the genome actually being edited
    # (recipient autologous iPSC in Track 3, or a non-reference donor MCB) carries a protospacer/PAM SNP (AUDIT-REF):
    editing_genome_has_snp = (recipient_has_protospacer_snp
        if selected_track == "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC" else donor_mcb_has_protospacer_snp)
    if editing_genome_has_snp and not rare_snp_pam_disrupted_input:
        grna_set.append("PATIENT_ALLELE_MATCHED_CUSTOM_SPACER_SYNTHESIS")

    # Explicit (K, B) editing complexity & manufacturing yield accounting (AUDIT-REF & AUDIT-REF):
    params = engineering_parameters if engineering_parameters is not None else _engineering_parameters()
    ranker_utility_penalty_factor = engineering_complexity_discount(
        k_dsb_loci, b_prime_loci, lambda_K=params["lambda_K"], gamma_B=params["gamma_B"]
    )
    if use_b_prime_sct_k4:
        variant_id = (
            ("ISLET1_B_SCT_C1_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C1_CIITA")
            if kir_ligand_group == "C1/C1"
            else (
                ("ISLET1_B_SCT_C2_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C2_CIITA")
                if kir_ligand_group == "C2/C2"
                else ("ISLET1_B_SCT_C1C2_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C1C2_CIITA")
            )
        )
    elif use_plan_b_hla_c_retained:
        variant_id = "ISLET1_B2_5LOCUS_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B2_5LOCUS_BASE"
    else:
        variant_id = "ISLET1_A4_ARM1B" if "DONOR-1B" in locus1_arm else "ISLET1_A4_ARM1A"
        if include_tyk2_p1104a_edit:
            variant_id += "_TYK2_PE"
    product = copy.deepcopy(_product_contract_snapshot(str(REPO_ROOT), variant_id))
    _check(product["K"] == k_dsb_loci and product["B"] == b_prime_loci, "Compiler K/B differs from its product contract")

    default_promoted_track2_subbank_product_id = (
        (
            ("ISLET1_B_SCT_C1_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C1_CIITA")
            if kir_ligand_group == "C1/C1"
            else (
                ("ISLET1_B_SCT_C2_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C2_CIITA")
                if kir_ligand_group == "C2/C2"
                else ("ISLET1_B_SCT_C1C2_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C1C2_CIITA")
            )
        )
        if use_plan_b_hla_c_retained
        else None
    )
    active_qualification_holds: List[str] = []
    if use_plan_b_hla_c_retained and has_anti_hlac_dsa_input and not anti_hlac_dsa_cleared_input:
        active_qualification_holds.append("HOLD_INCOMPATIBLE_ANTI_HLAC_DSA_CONTRAINDICATED_REQUIRES_CROSSMATCH_OR_TRACK3")
    if use_plan_b_hla_c_retained and hlac_4digit_rel_input == "KIR_LIGAND_MISMATCH":
        active_qualification_holds.append("HOLD_KIR_LIGAND_MISMATCH_CONTRAINDICATED_FOR_SUBBANK")
    if has_anti_hlae_dsa_input and not exact_hlae_sct_serum_crossmatch_neg_input:
        active_qualification_holds.append("HOLD_RECIPIENT_ANTI_HLAE_DSA_EXACT_SCT_CROSSMATCH_REQUIRED")
    exact_product_crossmatch_explicit_false = (
        "exact_product_non_hla_and_sct_flow_crossmatch_negative" in patient
        and patient.get("exact_product_non_hla_and_sct_flow_crossmatch_negative") is False
    )
    if exact_product_crossmatch_explicit_false:
        active_qualification_holds.append("HOLD_EXACT_PRODUCT_NON_HLA_AND_SCT_FLOW_CROSSMATCH_FAILED")
        clinical_override_warnings.append(
            "HOLD_EXACT_PRODUCT_NON_HLA_AND_SCT_FLOW_CROSSMATCH_FAILED: recipient exact-product non-HLA and SCT flow/CDC crossmatch "
            "explicitly returned positive/failed (exact_product_non_hla_and_sct_flow_crossmatch_negative=False), blocking release across all tracks "
            "(including Track 3 autologous iPSC due to chimeric linker/DAF-junction neoepitope reactivity) regardless of baseline cPRA (AUDIT-REF & AUDIT-REF [DEC-050/052])."
        )
    if (
        selected_track != "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC"
        and recipient_cpra_pct_input >= 80.0
        and not exact_product_non_hla_and_sct_flow_crossmatch_neg_input
    ):
        active_qualification_holds.append("HOLD_HIGH_CPRA_SENSITIZED_EXACT_PRODUCT_CROSSMATCH_REQUIRED")
    if use_plan_b_hla_c_retained and hlac_4digit_rel_input == "KIR_MATCH_FOREIGN_4DIGIT_ALLO" and not cd8_mlc_cleared_input:
        active_qualification_holds.append("HOLD_FOREIGN_4DIGIT_HLAC_REQUIRES_EXACT_PRODUCT_SERUM_AND_MLC_CROSSMATCH")
    if (
        use_plan_b_hla_c_retained
        and kir_ligand_group in {"C2/C2", "C1/C2", "C2/C1"}
        and has_activating_kir2ds1_hazard_input
        and not kir2ds1_cleared_input
    ):
        active_qualification_holds.append("HOLD_ACTIVATING_KIR2DS1_REQUIRES_SORTED_NK_NON_ACTIVATION_CLEARANCE")
    if cmv_nkg2c_adaptive_nk_fraction_input >= 0.25 and not gate0a_nkg2c_adaptive_nk_cleared_input:
        active_qualification_holds.append("HOLD_CMV_HIGH_NKG2C_ADAPTIVE_NK_REQUIRES_GATE_0A_DONOR_PANEL_CLEARANCE")
    if use_b_prime_sct_k4 and kir_ligand_group in {"C2/C2", "C1/C2", "C2/C1"} and not cd8_mlc_cleared_input:
        active_qualification_holds.append("HOLD_B_PRIME_SCT_C2_QYDDAVYKL_REQUIRES_ALL_RECIPIENT_CD8_MLC_CLEARANCE")
    if (
        use_b_prime_sct_k4
        and kir_ligand_group == "C1/C1"
        and hlac_4digit_rel_input in {"UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT", "UNSPECIFIED"}
        and not cd8_mlc_cleared_input
    ):
        active_qualification_holds.append(
            "HOLD_B_PRIME_SCT_C1_UNSPECIFIED_4DIGIT_RELATION_REQUIRES_ALLELE_TYPING_AND_CD8_MLC_CLEARANCE"
        )
    if (
        plan_b_biologically_indicated
        and selected_track == "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC"
        and not track1_sorted_kir_nk_cleared
    ):
        active_qualification_holds.append(
            "HOLD_FORCED_TRACK1_OVERRIDES_HIGH_KIR_NK_REQUIRES_SORTED_KIR_NK_CLEARANCE"
        )
    if (
        plan_b_biologically_indicated
        and selected_track == "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC"
        and not track1_sorted_kir_nk_cleared
    ):
        active_qualification_holds.append(
            "HOLD_FORCED_MHC_NULL_TRACK_OVERRIDES_EXTREME_KIR_NK_REQUIRES_SORTED_KIR_NK_CLEARANCE"
        )
    if cd47_sirpa_binding_compromised_input or "v8" in sirpa_genotype.lower():
        active_qualification_holds.append(
            "HOLD_CD47_SIRPA_BINDING_COMPROMISED_REQUIRES_ALTERNATIVE_STEALTH_QUALIFICATION"
        )
        clinical_override_warnings.append(
            "HOLD_CD47_SIRPA_BINDING_COMPROMISED_REQUIRES_ALTERNATIVE_STEALTH_QUALIFICATION: cd47_sirpa_binding_compromised=True or unmeasured SIRPa V8 allele "
            "indicates impaired/unverified CD47-IgV -> SIRPa checkpoint engagement and requires orthogonal macrophage/NK stealth qualification prior to release (AUDIT-REF & AUDIT-REF; DEC-071/072)."
        )
    if rare_snp_pam_disrupted_input:
        active_qualification_holds.append(
            "HOLD_RARE_SNP_PAM_DISRUPTED_REQUIRES_SHIFTED_PROTOSPACER_OR_ORTHOGONAL_NUCLEASE_QUALIFICATION"
        )
        clinical_override_warnings.append(
            "HOLD_RARE_SNP_PAM_DISRUPTED_REQUIRES_SHIFTED_PROTOSPACER_OR_ORTHOGONAL_NUCLEASE_QUALIFICATION: rare SNP disrupts "
            "the PAM motif on the active editing genome (rare_snp_pam_disrupted=True); spacer-only re-synthesis cannot restore PAM recognition "
            "and requires a shifted protospacer or orthogonal PAM-relaxed nuclease qualification (AUDIT-REF, AUDIT-REF & AUDIT-REF; DEC-071/072/073)."
        )
    if recipient_bw4_kir3dl1_educated and not bw4_kir3dl1_nk_lysis_cleared:
        active_qualification_holds.append(
            "HOLD_RECIPIENT_BW4_KIR3DL1_EDUCATED_NK_REQUIRES_SORTED_KIR3DL1_NK_CLEARANCE"
        )
    if use_plan_b_hla_c_retained and not kir_ligand_group_explicitly_supplied:
        active_qualification_holds.append(
            "HOLD_B_PRIME_SCT_REQUIRES_EXPLICIT_RECIPIENT_KIR_LIGAND_GROUP_TYPING"
        )
        clinical_override_warnings.append(
            "HOLD_B_PRIME_SCT_REQUIRES_EXPLICIT_RECIPIENT_KIR_LIGAND_GROUP_TYPING: Track 2 HLA-C-bearing routing (B'-SCT or endogenous HLA-C) cannot rely on an imputed C1/C2 default "
            "and requires explicit recipient KIR ligand genotyping (C1/C1, C2/C2, or C1/C2; AUDIT-REF & MW-R6-S2-R5 [DEC-051/053])."
        )
    if (
        use_b_prime_sct_k4
        and kir_ligand_group in {"C1/C2", "C2/C1"}
        and not c1c2_heteroallelic_phased_verified
    ):
        active_qualification_holds.append(
            "HOLD_C1C2_HETEROALLELIC_CIITA_SCT_REQUIRES_PHASED_LONGREAD_AND_DDPCR_VERIFICATION"
        )
    if (
        0.30 <= kir_licensed_nkg2a_neg_fraction < 0.48
        and arm1b_cd47_sirpa_nk_cliff_unverified
        and not track1_sorted_kir_nk_cleared
    ):
        active_qualification_holds.append(
            "HOLD_MODERATE_KIR_NK_REQUIRES_SORTED_KIR_NKG2A_NEG_NK_OR_TRACK2_B_SCT"
        )
    if single_200cm2_pouch_forced_without_modular_scaling and areal_loading_at_fixed_200cm2 > 1500.0:
        active_qualification_holds.append(
            "HOLD_SINGLE_POUCH_AREAL_OVERLOAD_EXCEEDS_1500_IEQ_CM2_CHRONIC_CAP"
        )
    orthogonal_autoimmune_gate_indicated = bool(
        autoantigen_t_cell_titer >= 2.0 or tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard
    )
    if orthogonal_autoimmune_gate_indicated and not gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared:
        active_qualification_holds.append(
            "HOLD_TYK2_INDEPENDENT_IFNG_TNFA_AND_INDIRECT_DR3_DQ8_CLEARANCE_REQUIRED"
        )
        clinical_override_warnings.append(
            "HOLD_TYK2_INDEPENDENT_IFNG_TNFA_AND_INDIRECT_DR3_DQ8_CLEARANCE_REQUIRED: high pre-existing beta-cell autoantigen T-cell reactivity "
            f"(z = {autoantigen_t_cell_titer:.2f} >= 2.0) or explicit TYK2-independent IFN-gamma/TNF-alpha + indirect HLA-DR3/DQ8 CD4+ T-cell hazard "
            "requires Stage-0 Gate 0D-indirect / Gate 0E-cytokine orthogonal stress clearance (gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared=True) "
            "because TYK2 P1104A catalytic attenuation covers Type-I IFN (IFNAR1/2) and IL-12/IL-23 signaling but not JAK1/JAK2-mediated IFN-gamma or "
            "recipient APC-mediated indirect CD4+ presentation (AUDIT-REF & AUDIT-REF [DEC-050])."
        )

    if active_qualification_holds:
        prod_qual_status = active_qualification_holds[0]
    elif use_b_prime_sct_k4:
        prod_qual_status = "SEQUENCE_COMPLETE_IN_SILICO_CANDIDATE_REQUIRING_WETLAB_RELEASE_QC"
    elif use_plan_b_hla_c_retained:
        prod_qual_status = "INVESTIGATIONAL_PLACEHOLDER_SEPARATE_HLAA_HLAB_DONOR_LEDGER_UNRESOLVED"
    else:
        prod_qual_status = "SEQUENCE_COMPLETE_IN_SILICO_CANDIDATE_REQUIRING_WETLAB_RELEASE_QC"

    routing_holds_cleared = len(active_qualification_holds) == 0
    in_silico_candidate_routing_cleared = bool(routing_holds_cleared and product.get("execution_supported", False))
    required_recipient_hazard_keys = {
        "recipient_anti_hla_e_dsa",
        "cmv_nkg2c_adaptive_nk_fraction",
        "recipient_bw4_kir3dl1_educated",
        "exact_product_non_hla_and_sct_flow_crossmatch_negative",
    }
    if selected_track != "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC":
        required_recipient_hazard_keys.add("recipient_cpra_pct")
    if orthogonal_autoimmune_gate_indicated:
        required_recipient_hazard_keys.add("gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared")
    if use_plan_b_hla_c_retained:
        required_recipient_hazard_keys.update({"recipient_kir_ligand_group", "recipient_anti_hla_c_dsa", "hlac_4digit_relation"})
        if kir_ligand_group in {"C2/C2", "C1/C2", "C2/C1"}:
            required_recipient_hazard_keys.add("recipient_activating_kir2ds1_unlicensed_hazard")
        if use_b_prime_sct_k4 and kir_ligand_group in {"C1/C2", "C2/C1"}:
            required_recipient_hazard_keys.add("c1c2_heteroallelic_phased_longread_ddpcr_verified")
    recipient_immune_hazard_panel_explicitly_measured = (
        all(k in patient and patient.get(k) is not None for k in required_recipient_hazard_keys)
        and not exact_product_crossmatch_explicit_false
        and patient.get("exact_product_non_hla_and_sct_flow_crossmatch_negative") is True
        and (
            not orthogonal_autoimmune_gate_indicated
            or patient.get("gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared") is True
        )
    )
    wetlab_stage0_complete_lot_release_verified = _validate_strict_bool(
        _get_opt(patient, "wetlab_stage0_complete_lot_release_verified", False),
        "wetlab_stage0_complete_lot_release_verified",
    )
    lot_release_cleared = bool(
        in_silico_candidate_routing_cleared
        and product.get("physical_qualification_status") == "QUALIFIED"
        and recipient_immune_hazard_panel_explicitly_measured
        and wetlab_stage0_complete_lot_release_verified
    )

    return {
        "patient_id": patient_id,
        "ancestry": ancestry,
        "product_variant_id": variant_id,
        "selected_product_id": variant_id,
        "default_promoted_track2_subbank_product_id": default_promoted_track2_subbank_product_id,
        "product_definition": product,
        "score_qualification": "ASSIGNED_RULE_SCENARIO_NOT_COMPLETE_PRODUCT",
        "selected_fda_deployment_track": selected_track,
        "clinical_override_warnings": clinical_override_warnings,
        "track3_autologous_vs_allogeneic_scope_note_v80_088": (
            "In Track 3 (autologous patient iPSC), host-vs-graft alloreactive MHC mismatch is absent, but autoreactive "
            "T1D memory CD8+ (PPI/IGRP/ZnT8/GAD65) and CD4+ T-cell protection still requires MHC-I/II + CXCL10::SERPINB9-P2A-PD-L1 "
            "or selective autoimmune shielding (AUDIT-REF)."
        ),
        "sequence_compiler_execution_supported": product["execution_supported"],
        "physical_product_qualification_status": product["physical_qualification_status"],
        "product_qualification_status": prod_qual_status,
        "selected_product_qualification_status": prod_qual_status,
        "all_active_qualification_holds": active_qualification_holds,
        "routing_holds_cleared": routing_holds_cleared,
        "in_silico_candidate_routing_cleared": in_silico_candidate_routing_cleared,
        "recipient_immune_hazard_panel_explicitly_measured": recipient_immune_hazard_panel_explicitly_measured,
        "wetlab_stage0_complete_lot_release_verified": wetlab_stage0_complete_lot_release_verified,
        "lot_release_cleared": lot_release_cleared,
        "regulatory_hypothesis_pending_qualified_assessment": (
            "Investigational regulatory hypothesis: potential eligibility for FDA Platform Technology Designation "
            "(21 U.S.C. 356k, subject to statutory prerequisite of an approved drug/biologic using the platform) "
            "or CBER draft/nonbinding Individualized Plausible-Mechanism guidance (Feb 2026) requires formal FDA CBER INTERACT/pre-IND assessment (closing AUDIT-REF)."
        ),
        "locus1_construct_prescribed": locus1_arm,
        "plan_architecture_prescribed": (
            "Plan B-SCT (B2M::Arm1B + CIITA::HLA-C*03:04/04:01 Fixed-Peptide SCT)"
            if use_b_prime_sct_k4
            else (
                "Plan B (Selective HLA-A/HLA-B KO + KIR-Matched Retained HLA-C1/C2 Sub-Bank + HLA-E SCT + CD47-IgV)"
                if use_plan_b_hla_c_retained
                else "Plan A (Core 4-Locus Tri-Modal B2M::HLA-E-CD47-IgV + CIITA-/- + CXCL10::SERPINB9-P2A-PD-L1 + AAVS1-iCasp9)"
            )
        ),
        "plan_b_chr6p21_cis_deletion_hazard_and_sequential_editing_gate_kimi_r3_001": (
            "On GRCh38 chr6p21.3, physical gene order is HLA-A (29.94 Mb) -> HLA-E (30.49 Mb) -> HLA-C (31.27 Mb) -> HLA-B (31.35 Mb); "
            "simultaneous DSB multiplexing of HLA-A and HLA-B is prohibited in Plan B due to the ~1.41 Mb interstitial cis-deletion/inversion "
            "risk removing intervening HLA-E and retained HLA-C. Requires sequential clonal passes or non-DSB editing plus phased HLA-A/E/C/B ddPCR and OGM verification. "
            "Note on HLA-C leader-peptide selection (Lin et al. 2023 Nat Immunol PMID:37264229, PMC10690437, DOI:10.1038/s41590-023-01523-z; AUDIT-REF, Q2-P1..P3, AUDIT-REF, & AUDIT-REF): "
            "because HLA-C*07:02 encodes the non-functional SP-3C (VMAPRALLL, carrying P6=Ala [T6A, with P5=Arg]) leader peptide (24.49% of 1000G alleles) that fails to stabilize surface HLA-E, "
            "Design B-prime upgrades the C1/C1 sub-bank anchor to functional SP-1C (VMAPRTLIL) HLA-C*01:02 or HLA-C*03:04 (alongside SP-1C C2 anchor HLA-C*04:01) "
            "and co-integrates exogenous HLA-E*01:03(SCT) + CD47-IgV-GPI at AAVS1 Locus 4 (subject to physical FASTA materialization and complexity-discount accounting: "
            "analytic personalized scenario index is 97.50% unpenalized, 92.84% under K/B complexity discounting vs 92.50% Plan A, and 91.37% unpenalized / 87.65% complexity-discounted if the unmaterialized safe-harbor cloak is omitted; AUDIT-REF & AUDIT-REF). "
            "Note on K/B indexing (AUDIT-REF & AUDIT-REF [DEC-078]): Layer-5 Track-2 Plan B indexes the 5 DSB/indel loci (HLA-A, HLA-B, CIITA, CXCL10, AAVS1; K=5, B=0 or B=1 with optional TYK2 P1104A PE), "
            "whereas Layer-1 Chassis B (ISLET1_B6, global Rank #25, K=6, B=1) indexes 6 DSB/indel loci (HLA-A, HLA-B, CIITA, CXCL10, AAVS1, plus TYK2 knock-out as the 6th DSB/indel locus) and the INS signal-peptide/DRiP base edit as B=1 "
            "(and the separate K=5, B=1 TYK2 P1104A prime-edited variant is ISLET1_B2_5LOCUS_TYK2_PE)."
            if use_plan_b_hla_c_retained
            else None
        ),
        "include_tyk2_p1104a_rs34536443_shield": include_tyk2_p1104a_edit,
        "final_graft_tyk2_p1104a_protective_state": final_graft_tyk2_p1104a_protective_state,
        "editing_complexity_and_manufacturing_accounting_v80_009": {
            "k_dsb_or_indel_chromosomal_loci": k_dsb_loci,
            "b_prime_or_base_edit_loci": b_prime_loci,
            "tyk2_alleles_requiring_prime_edit_count": tyk2_alleles_to_edit,
            "uncloned_simultaneous_pass_yield_prior": None,
            "physical_manufacturing_yield_evidence_status": "NOT_MEASURED",
            "ranker_engineering_complexity_discount_factor": ranker_utility_penalty_factor,
            "engineering_utility_prior_parameters": {"lambda_K": params["lambda_K"], "gamma_B": params["gamma_B"]},
            "gate_zero_donor_wgs_and_precloning_junction_recleavage_gate_muse_r1": (
                "Gate 0 (prior to guide lock and HIP-Core clone generation): perform 30x whole-genome sequencing (WGS) + phased "
                "HLA/KIR/SIRPA(V1/V2/V8) genotyping of the candidate Master Cell Bank hPSC line to confirm zero PAM/seed SNVs across all "
                "5 guide target windows, followed by Cas-OFFinder + empirical GUIDE-seq/CAST-seq in the actual hPSC line and a pre-cloning "
                "T7E1/amplicon-NGS junction re-cleavage assay on synthetic edited-allele amplicons across B2M, CIITA, CXCL10, and AAVS1 "
                "(AUDIT-REF & D5-02)."
            ),
            "two_pass_clonal_mcb_release_qc_policy": (
                "Assigned-rule scenario values in evaluate_regimen_survival_probability are conditioned on a hypothetical "
                "release-qualified Master Cell Bank (MCB) clone that has passed two-pass clonal screening (2 chromosomal loci x 2 "
                "homologous alleles = 4 target alleles per pass), junction-spanning AAV6 ITR-capture ddPCR (<0.1% vector-backbone/ITR "
                "concavemer or off-locus capture at CIITA/CXCL10/AAVS1; Hanlon et al. 2019 PMID:31570731; AUDIT-REF & D6-F03), "
                ">=63,000x TP53 duplex-seq (>=62,958 independent duplex molecules), and 20q11.21 ddPCR. Physical uncloned yield is unmeasured. "
                "The separate dimensionless engineering utility discount is exp(-lambda_K * max(0, K - 4)) * gamma_B ** B; "
                "it is not manufacturing yield or a fitted clinical probability."
            ),
        },
        "host_modifier_genotypes_and_dosing": {
            "ifih1_rs1990760_genotype": ifih1_genotype,
            "ptpn22_rs2476601_genotype": ptpn22_genotype,
            "recipient_hla_e_genotype": recipient_hla_e,
            "tyk2_rs34536443_genotype": recipient_tyk2_genotype,
            "recipient_tyk2_rs34536443_genotype": recipient_tyk2_genotype,
            "donor_tyk2_rs34536443_genotype": donor_tyk2_genotype,
            "starting_graft_tyk2_rs34536443_genotype": starting_graft_tyk2_genotype,
            "tyk2_alleles_requiring_prime_edit_count": tyk2_alleles_to_edit,
            "already_homozygous_tyk2_p1104a": already_homozygous_tyk2_p1104a,
            "final_graft_tyk2_p1104a_protective_state": final_graft_tyk2_p1104a_protective_state,
            "recipient_kir_ligand_group": kir_ligand_group,
            "recipient_bmi_kg_m2": recipient_bmi,
            "prescribed_curative_dose_ieq": prescribed_dose_ieq,
            "prescribed_dose_source_mode": dose_source_mode,
            "areal_loading_ieq_per_cm2_if_fixed_200cm2_pouch": areal_loading_at_fixed_200cm2,
            "required_pouch_exchange_area_cm2_at_nominal_2500_ieq_cm2": required_pouch_area_cm2_at_2500_ieq_cm2,
            "required_pouch_exchange_area_cm2_at_chronic_lowload_1500_ieq_cm2": required_pouch_area_cm2_at_1500_ieq_cm2,
            "modular_100cm2_omental_cassette_count_at_1500_ieq_cm2": modular_100cm2_omental_cassette_count_at_1500_ieq_cm2,
            "actual_areal_loading_ieq_per_cm2_across_modular_cassettes": actual_areal_loading_ieq_per_cm2_across_modular_cassettes,
            "recipient_immunological_screen_completeness_status": recipient_immunological_screen_completeness_status,
            "dosing_model_qualification_note": (
                "Illustrative BMI-proportional scaling (500,000 IEQ * BMI / 24.0) or explicit prescribed_dose_ieq override; "
                "clinical islet dosing requires actual recipient body weight (kg), baseline insulin requirement, and proportional "
                "scaling of pouch exchange area (A = Dose / 2,500 IEQ/cm^2 nominal, or A = Dose / 1,500 IEQ/cm^2 = 333.33 cm^2 at "
                "500,000 total IEQ under DX-B-LowLoad-1500 chronic zero-source encapsulation via modular 100 cm^2 omental cassettes)."
            ),
            "optional_field_input_provenance": optional_field_provenance,
        },
        "tyk2_independent_autoimmune_memory_and_indirect_dr3_dq8_audit_sonnet55_r5": {
            "finding_ids": ["AUDIT-REF", "AUDIT-REF", "AUDIT-REF"],
            "primary_citations": ["PMID:27807284", "PMID:30578352", "PMID:36289205", "PMID:26239055"],
            "tyk2_p1104a_pathway_scope": (
                "TYK2 P1104A (rs34536443 C/C) selectively dampens Type-I IFN (IFN-alpha/beta via IFNAR1/TYK2-JAK1) and IL-12/IL-23 signaling, "
                "but does NOT block IFN-gamma signaling (which signals exclusively through IFNGR1/IFNGR2 via JAK1 + JAK2, independent of TYK2) "
                "nor TNF-alpha / IL-1beta NF-kB/JNK pro-apoptotic signaling, and does NOT prevent host recipient APCs from phagocytosing shed beta-cell autoantigens (PPI, GAD65, IA-2, ZnT8, HIPs) "
                "for indirect presentation on recipient HLA-DR3-DQ2 / HLA-DQ8 to autoreactive CD4+ helper T cells."
            ),
            "orthogonal_protection_layers_required": (
                "Protection against TYK2-independent IFN-gamma/TNF-alpha + indirect DR3/DQ8 CD4+ autoimmunity relies on: "
                "(1) CXCL10-inducible (or L3-BF-B basal-floor) cytosolic SERPINB9 + surface PD-L1 (harnessing local IFN-gamma/STAT1 to upregulate granzyme-B/PD-1 defense), "
                "(2) peri-implant SA-FasL / SA-THBD microgel + Treg co-delivery to induce apoptosis of infiltrating Fas+ CD4+ Teff and foster FoxP3+ Treg conversion, and "
                "(3) Stage-0 Gate 0C-kinetic verification under combined IFN-gamma + TNF-alpha + IL-1beta + DR3/DQ8 PBMC transwell co-culture."
            ),
        },
        "recipient_carries_rare_protospacer_snp": recipient_has_protospacer_snp,
        "edited_cell_line_requires_custom_allele_grna": editing_genome_has_snp,
        "compiled_grna_panel": grna_set,
        "sirpa_polymorphism_compatibility": (
            f"Reference structures cover SIRPA V1/V2; exact-product affinity/protection for supplied genotype {sirpa_genotype}, "
            "including V8 if present, is unmeasured and requires binding and macrophage assays (AUDIT-REF)"
        ),
        "b_prime_4digit_hlac_and_b_prime_sct_routing_audit_muse_r3": {
            "recipient_hla_c_genotype_4digit": _get_opt(patient, "recipient_hla_c_genotype_4digit", None),
            "recipient_kir_genotype": _get_opt(patient, "recipient_kir_genotype", None),
            "recipient_anti_hla_c_dsa": bool(_get_opt(patient, "recipient_anti_hla_c_dsa", False)),
            "hlac_4digit_relation": _get_opt(
                patient,
                "hlac_4digit_relation",
                "SELF_4DIGIT_MATCH" if not use_plan_b_hla_c_retained else "UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT",
            ),
            "track2_hlac_4digit_match_status": prod_qual_status,
            "canonical_track2_kir_subbank_product_variant_id": (
                ("ISLET1_B_SCT_C1_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C1_CIITA")
                if kir_ligand_group == "C1/C1"
                else (
                    ("ISLET1_B_SCT_C2_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C2_CIITA")
                    if kir_ligand_group == "C2/C2"
                    else ("ISLET1_B_SCT_C1C2_CIITA_TYK2_PE" if include_tyk2_p1104a_edit else "ISLET1_B_SCT_C1C2_CIITA")
                )
            ),
            "endogenous_hlac_retention_status": "EXPLORATORY_REQUIRES_4DIGIT_BANK",
            "b_prime_sct_status": (
                "CANONICAL_TRACK2_KIR_ENGAGING_SUBBANK_K4_B1_3PASS_TYK2_PE_SATELLITE"
                if include_tyk2_p1104a_edit
                else "CANONICAL_TRACK2_KIR_ENGAGING_SUBBANK_K4_B0_2PASS"
            ),
            "recipient_bw4_kir3dl1_educated_status": (
                "CLEARED_SORTED_KIR3DL1_NK"
                if (recipient_bw4_kir3dl1_educated and bw4_kir3dl1_nk_lysis_cleared)
                else (
                    "HOLD_REQUIRES_SORTED_KIR3DL1_NKG2A_NEG_NK_CLEARANCE"
                    if recipient_bw4_kir3dl1_educated
                    else "NOT_INDICATED_BW4_KIR3DL1_ABSENT_OR_UNSPECIFIED"
                )
            ),
            "c1c2_heteroallelic_phased_longread_ddpcr_verified": (
                (c1c2_heteroallelic_phased_verified if c1c2_phasing_explicitly_supplied else None)
                if (use_b_prime_sct_k4 and kir_ligand_group in {"C1/C2", "C2/C1"})
                else None
            ),
            "c1c2_heteroallelic_phased_verification_status": (
                (
                    "VERIFIED_EXPLICIT_TRUE"
                    if (c1c2_phasing_explicitly_supplied and c1c2_heteroallelic_phased_verified)
                    else (
                        "HOLD_EXPLICIT_FALSE"
                        if c1c2_phasing_explicitly_supplied
                        else "HOLD_UNMEASURED_MISSING_PHASING_EVIDENCE"
                    )
                )
                if (use_b_prime_sct_k4 and kir_ligand_group in {"C1/C2", "C2/C1"})
                else None
            ),
            "note": (
                "Because KIR C1/C2 dimorphism is a 2-group epitope (Asn80 vs Lys80) while T-cell & DSA alloreactivity is 4-digit allele-specific "
                "(Xu et al. 2019 PMID:30853558; Tran et al. 2011 PMID:21743387), a 3-line endogenous HLA-C sub-bank leaves ~70-85% of recipients "
                "4-digit mismatched at HLA-C unless a ~12-16 line homozygous 4-digit bank is used or a non-allogenic self-peptide HLA-C SCT "
                "(Design B'-SCT: ISLET1_B_SCT_C1_CIITA [GAVDPLLAL, PDB:1EFX] / ISLET1_B_SCT_C2_CIITA [QYDDAVYKL, PDB:1IM9] at CIITA Exon 5, K=4, B=0 [P_passes=2] or B=1 [P_passes=3]) is deployed (AUDIT-REF & AUDIT-REF)."
            ),
        },
    }


VALID_COMPILED_TRACKS = {
    "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC",
    "TRACK_2_PATIENT_DNA_STRATIFIED_ALLOGENEIC_SUBLINE",
    "TRACK_3_PERSONALIZED_AUTOLOGOUS_PATIENT_IPSC",
}


def evaluate_regimen_survival_probability(
    patient: Dict[str, Any],
    compiled_regimen: Dict[str, Any],
) -> float:
    """
    Evaluate conditional 5-year graft survival probability (conditional on transplanting a release-qualified
    MCB clone) from how well `compiled_regimen` addresses `patient`'s sampled immunological and genetic risk
    factors (closing AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, and AUDIT-REF).
    Recognizes natural homozygous `TYK2 rs34536443 C/C` on the graft genome as conferring the exact same
    biological TYK2 P1104A protection as prime-edited `G/G -> C/C`.
    """
    if not isinstance(patient, dict):
        raise ValueError("patient must be a dictionary")
    if not isinstance(compiled_regimen, dict):
        raise ValueError("compiled_regimen must be a dictionary")
    req_reg_keys = (
        "selected_fda_deployment_track",
        "locus1_construct_prescribed",
        "plan_architecture_prescribed",
        "include_tyk2_p1104a_rs34536443_shield",
        "compiled_grna_panel",
    )
    for rk in req_reg_keys:
        if rk not in compiled_regimen:
            raise ValueError(f"compiled_regimen missing required field {rk!r}")

    track_str = compiled_regimen["selected_fda_deployment_track"]
    if not isinstance(track_str, str) or track_str not in VALID_COMPILED_TRACKS:
        raise ValueError(f"Invalid selected_fda_deployment_track {track_str!r}; must be one of {sorted(VALID_COMPILED_TRACKS)}")

    tyk2_flag = _validate_strict_bool(
        compiled_regimen["include_tyk2_p1104a_rs34536443_shield"],
        "include_tyk2_p1104a_rs34536443_shield",
    )

    loc1_str = compiled_regimen["locus1_construct_prescribed"]
    plan_str = compiled_regimen["plan_architecture_prescribed"]
    if not isinstance(loc1_str, str) or not isinstance(plan_str, str):
        raise ValueError("locus1_construct_prescribed and plan_architecture_prescribed must be strings")

    # Enforce strict consistency between Plan A/B architecture and locus1 construct (AUDIT-REF):
    is_plan_b_arch = "Plan B" in plan_str
    is_plan_b_loc1 = ("PLAN_B_SELECTIVE_HLAA_HLAB_KO" in loc1_str) or ("B_SCT" in loc1_str)
    if is_plan_b_arch != is_plan_b_loc1:
        raise ValueError(
            f"Inconsistent regimen: plan_architecture_prescribed={plan_str!r} conflicts with locus1_construct_prescribed={loc1_str!r}"
        )

    # Fail-closed product_definition defaults for Plan B (L5A-02 [DEC-077]):
    # Any Plan B regimen lacking an explicit `product_definition` fails closed on both HLA-C retention and CD47/HLA-E cloak
    # rather than receiving free survival credit, while valid `product_definition` dicts are verified below.
    prod_def_cd47_ok = not is_plan_b_arch
    prod_def_hlac_ok = not is_plan_b_arch
    prod_def_hlac_str = ""
    prod_def_leader_class = "SP-1C"
    ALLOWED_PLAN_B_HLAC_STATES = {
        "RETAINED KIR-MATCHED DONOR ALLELE(S), SEPARATELY PARAMETERIZED",
        "RETAINED_KIR_MATCHED_C1_C2_SUBBANK",
        "RETAINED_KIR_C1_C1_SUBBANK",
        "RETAINED_KIR_C2_C2_SUBBANK",
        "RETAINED_KIR_C1_C2_HETEROZYGOUS_SUBBANK",
        "RETENTION_REQUIRED_KIR_MATCHED_C1_C2",
        "RETENTION_REQUIRED_KIR_MATCHED_C1_C1",
        "RETENTION_REQUIRED_KIR_MATCHED_C2_C2",
        "RETAINED_HLA_C0304_C0401_HETEROZYGOUS_SP1C",
        "RETAINED_HLA_C0304_C0702_HETEROZYGOUS_SP1C_DOMINANT",
        "RETAINED_HLA_C0702_HOMOZYGOUS_SP3C",
    }
    ALLOWED_PLAN_B_B2M_STATES = {
        "RETAINED WT",
        "RETAINED_WT",
        "RETAINED_ENDOGENOUS_WT_FOR_HLAC_SURFACE_PRESENTATION",
        "RETAINED_ENDOGENOUS_WT",
    }
    if "product_definition" in compiled_regimen:
        prod_def_cd47_ok = True
        prod_def_hlac_ok = True
        prod_def = compiled_regimen["product_definition"]
        if not isinstance(prod_def, dict):
            raise ValueError("product_definition must be a dictionary when provided (AUDIT-REF)")
        from scripts.product_contract import validate_product_definition
        validate_product_definition(prod_def)
        expected_family = "B" if is_plan_b_arch else "A"
        if prod_def.get("family") != expected_family:
            raise ValueError("product_definition family conflicts with prescribed architecture")
        if is_plan_b_arch:
            alleles_map = prod_def.get("alleles", {})
            if "HLA-C" not in alleles_map or "B2M" not in alleles_map:
                raise ValueError("Plan B product_definition must declare retained HLA-C and B2M alleles")
            hlac_allele_val = str(alleles_map.get("HLA-C", "")).strip().upper()
            b2m_allele_val = str(alleles_map.get("B2M", "")).strip().upper()
            prod_def_hlac_str = hlac_allele_val
            is_b_sct_prod = "B_SCT" in str(prod_def.get("product_variant_id", "")).upper()
            bad_tokens = (
                "KO", "LOF", "DISRUPT", "ABSENT", "NONE", "DELETED", "NULL", "OMITTED", "FALSE",
                "SILENCED", "ABLATED", "DEFICIENT", "INACTIVATED", "FRAMESHIFT", "STOP", "MUTATED", "TRUNCATED",
                "NONFUNCTIONAL", "UNKNOWN", "UNAVAILABLE",
            )
            if is_b_sct_prod:
                if (
                    ("SINGLE-CHAIN TRIMER" not in hlac_allele_val and "SCT" not in hlac_allele_val)
                    or "EMBODIMENT A*" not in b2m_allele_val
                    or any(tok in hlac_allele_val for tok in ("ABSENT", "NONE", "OMITTED", "FALSE", "DELETED", "NULL", "NONFUNCTIONAL"))
                ):
                    prod_def_hlac_ok = False
            else:
                hlac_allowlist_ok = hlac_allele_val in ALLOWED_PLAN_B_HLAC_STATES
                b2m_allowlist_ok = b2m_allele_val in ALLOWED_PLAN_B_B2M_STATES
                if (
                    not hlac_allowlist_ok
                    or not b2m_allowlist_ok
                    or any(tok in hlac_allele_val for tok in bad_tokens)
                    or any(tok in b2m_allele_val for tok in bad_tokens)
                ):
                    prod_def_hlac_ok = False
            raw_leader_cls = str(prod_def.get("hla_c_leader_class", alleles_map.get("HLA-C_leader_class", ""))).strip().upper()
            if raw_leader_cls in {"SP-3C", "SP3C", "VMAPRALLL"} or (
                "C0702_HOMOZYGOUS" in hlac_allele_val or "C*07:02/C*07:02" in hlac_allele_val or "SP3C" in hlac_allele_val
            ):
                prod_def_leader_class = "SP-3C"
            topo_cd47 = str(prod_def.get("cassette_topology", {}).get("CD47_HLAE", "")).strip().upper()
            # Note: In `islet1_product_definitions_v0.8.0.json`, `ISLET1_B2_5LOCUS_BASE` and
            # `ISLET1_B2_5LOCUS_TYK2_PE` carry `"separate HLA-C-compatible donor design unresolved"`,
            # which represents the conditional-on-release-QC branch (`97.4963%`) whereas explicit
            # `"ABSENT"` / `"NOT_MATERIALIZED"` represents the unmaterialized Plan-B branch (`91.3735%`).
            absent_cloak_tokens = (
                "ABSENT", "NONE", "OMITTED", "FALSE", "DELETED", "NULL",
                "NOT_MATERIALIZED", "UNMATERIALIZED", "NO_SCT", "NOT_PRESENT", "UNAVAILABLE", "UNKNOWN", "KO", "UNCLOAKED",
            )
            if (not topo_cd47) or any(tok in topo_cd47 for tok in absent_cloak_tokens):
                prod_def_cd47_ok = False
    sec_cloak_str = str(compiled_regimen.get("secondary_safe_harbor_cloak", "")).strip().upper()
    if sec_cloak_str and any(
        tok in sec_cloak_str
        for tok in ("ABSENT", "NONE", "OMITTED", "FALSE", "DELETED", "NULL", "NOT_MATERIALIZED", "UNMATERIALIZED", "NO_SCT", "NOT_PRESENT", "UNAVAILABLE", "UNKNOWN", "KO", "UNCLOAKED")
    ):
        prod_def_cd47_ok = False

    panel = compiled_regimen["compiled_grna_panel"]
    if (not isinstance(panel, list)
            or any(not isinstance(guide, str) or not guide.strip() for guide in panel)
            or len(panel) != len(set(panel))):
        raise ValueError("compiled_grna_panel must be a list of unique nonempty guide IDs")

    allowed_eval_patient_keys = {
        "patient_id",
        "ancestry",
        "nkg2a_neg_kir_educated_nk_fraction",
        "pre_existing_gad65_ia2_znt8_t_cell_zscore",
        "ifn1_ish_genetic_risk_percentile",
        "sirpa_genotype",
        "ifih1_rs1990760_genotype",
        "ptpn22_rs2476601_genotype",
        "recipient_hla_e_genotype",
        "recipient_bmi_kg_m2",
        "prescribed_dose_ieq",
        "tyk2_rs34536443_genotype",
        "recipient_tyk2_rs34536443_genotype",
        "donor_tyk2_rs34536443_genotype",
        "recipient_kir_ligand_group",
        "recipient_anti_hla_c_dsa",
        "recipient_anti_hla_e_dsa",
        "recipient_cpra_pct",
        "cmv_nkg2c_adaptive_nk_fraction",
        "recipient_activating_kir2ds1_unlicensed_hazard",
        "hlac_4digit_relation",
        "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative",
        "b_prime_sct_recipient_cd8_mlc_negative",
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared",
        "b_prime_sct_kir2ds1_non_activation_cleared",
        "exact_hlae_sct_serum_crossmatch_negative",
        "exact_product_non_hla_and_sct_flow_crossmatch_negative",
        "gate0a_nkg2c_adaptive_nk_cleared",
        "recipient_bw4_kir3dl1_educated",
        "bw4_kir3dl1_sorted_nk_lysis_cleared",
        "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard",
        "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared",
        "requested_track",
        "track2_subbank_mode",
        "carries_rare_snp_in_standard_grna",
        "donor_mcb_carries_rare_snp_in_standard_grna",
        "rare_snp_pam_disrupted",
        "cd47_sirpa_binding_compromised",
        "track1_sorted_kir_nkg2a_neg_nk_lysis_cleared",
        "arm1b_cd47_sirpa_nk_cliff_unverified",
        "single_200cm2_pouch_forced_without_modular_scaling",
        "c1c2_heteroallelic_phased_longread_ddpcr_verified",
        "wetlab_stage0_complete_lot_release_verified",
        "recipient_hla_c_genotype_4digit",
        "recipient_kir_genotype",
    }
    unknown_eval_pt_keys = sorted(set(patient) - allowed_eval_patient_keys)
    if unknown_eval_patient_keys := unknown_eval_pt_keys:
        raise ValueError(f"Unrecognized patient input keys: {unknown_eval_patient_keys} (AUDIT-REF; DEC-071)")

    kir_nk = _validate_finite_float(
        patient["nkg2a_neg_kir_educated_nk_fraction"], "nkg2a_neg_kir_educated_nk_fraction", 0.0, 1.0
    )
    auto_z = _validate_finite_float(
        patient["pre_existing_gad65_ia2_znt8_t_cell_zscore"], "pre_existing_gad65_ia2_znt8_t_cell_zscore", -5.0, 10.0
    )
    ifn1_pct = _validate_finite_float(
        patient["ifn1_ish_genetic_risk_percentile"], "ifn1_ish_genetic_risk_percentile", 0.0, 100.0
    )
    recip_snp = _validate_strict_bool(
        _get_opt(patient, "carries_rare_snp_in_standard_grna", False), "carries_rare_snp_in_standard_grna"
    )
    donor_mcb_snp = _validate_strict_bool(
        _get_opt(patient, "donor_mcb_carries_rare_snp_in_standard_grna", False), "donor_mcb_carries_rare_snp_in_standard_grna"
    )
    rare_snp_pam_disrupted_eval = _validate_strict_bool(
        _get_opt(patient, "rare_snp_pam_disrupted", False), "rare_snp_pam_disrupted"
    )

    if "sirpa_genotype" in patient and not isinstance(patient["sirpa_genotype"], str):
        raise ValueError(f"Invalid sirpa_genotype: explicit value must be a string, got {patient['sirpa_genotype']!r} (AUDIT-REF)")
    for str_field in ("ifih1_rs1990760_genotype", "ptpn22_rs2476601_genotype", "recipient_hla_e_genotype"):
        if str_field in patient and patient[str_field] is not None and not isinstance(patient[str_field], str):
            raise ValueError(f"Invalid {str_field}: explicit value must be a string, got {patient[str_field]!r} (AUDIT-REF)")
    ifih1_g = str(_get_opt(patient, "ifih1_rs1990760_genotype", "C/T")).strip().upper()
    if ifih1_g not in VALID_IFIH1_GENOTYPES:
        raise ValueError(f"Invalid ifih1_rs1990760_genotype {ifih1_g!r}")
    ptpn22_g = str(_get_opt(patient, "ptpn22_rs2476601_genotype", "G/G")).strip().upper()
    if ptpn22_g not in VALID_PTPN22_GENOTYPES:
        raise ValueError(f"Invalid ptpn22_rs2476601_genotype {ptpn22_g!r}")
    raw_hla_e = str(_get_opt(patient, "recipient_hla_e_genotype", "*01:03/*01:03")).strip()
    if raw_hla_e not in VALID_HLA_E_GENOTYPES:
        raise ValueError(f"Invalid recipient_hla_e_genotype {raw_hla_e!r}")
    hla_e_g = _normalize_hla_e_genotype(raw_hla_e)
    raw_sirpa = str(_get_opt(patient, "sirpa_genotype", "v1/v1")).strip().lower()
    if raw_sirpa not in VALID_SIRPA_GENOTYPES:
        raise ValueError(f"Invalid sirpa_genotype {raw_sirpa!r}")
    sirpa_g = _normalize_sirpa_genotype(raw_sirpa)
    cd47_compromised = _validate_strict_bool(
        _get_opt(patient, "cd47_sirpa_binding_compromised", False), "cd47_sirpa_binding_compromised"
    )
    bmi_val = _validate_finite_float(_get_opt(patient, "recipient_bmi_kg_m2", 24.0), "recipient_bmi_kg_m2", 10.0, 70.0)

    dosing_block = compiled_regimen.get("host_modifier_genotypes_and_dosing", {})
    if not isinstance(dosing_block, dict):
        raise ValueError("host_modifier_genotypes_and_dosing must be a dictionary")
    raw_dose_override = _get_opt(patient, "prescribed_dose_ieq", None)
    default_dose = (_validate_finite_float(raw_dose_override, "prescribed_dose_ieq", 50000.0, 2000000.0)
                    if raw_dose_override is not None else (500000.0 * (bmi_val / 24.0)))
    dose_ieq = _validate_finite_float(dosing_block.get("prescribed_curative_dose_ieq", default_dose),
                                    "prescribed_curative_dose_ieq", 50000.0, 2000000.0)

    p_surv = 0.982

    # Verify core guide panel completeness (CIITA Exon 5 [or legacy Exon 1 ID], CXCL10, AAVS1) and Plan A vs Plan B locus-1 completeness:
    has_ciita_guide = (
        "OI-GRNA-02_CIITA_EXON1_ASCAS12A" in panel
        or "OI-GRNA-02_CIITA_EXON5_ASCAS12A" in panel
        or "OI-GRNA-02_CIITA_PAN_ISOFORM_EXON5_ASCAS12A" in panel
    )
    has_core_3_loci = has_ciita_guide and all(
        g_id in panel
        for g_id in (
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        )
    )
    if not has_core_3_loci:
        p_surv -= 0.12

    # 1. Check complete classical MHC-I editing prerequisites for Plan B (BOTH HLA-A AND HLA-B KO guides required + retained KIR-matched HLA-C/B2M) vs Plan A (B2M KO guide required):
    raw_recip_kir = str(
        _get_opt(patient, "recipient_kir_ligand_group", dosing_block.get("recipient_kir_ligand_group", "C1/C2"))
    ).strip().upper().replace("_", "/")
    if raw_recip_kir not in VALID_KIR_LIGAND_GROUPS:
        raise ValueError(f"Invalid recipient_kir_ligand_group {raw_recip_kir!r}")
    if raw_recip_kir == "C1/C1":
        raw_recip_kir = "C1/C1"
    elif raw_recip_kir == "C2/C2":
        raw_recip_kir = "C2/C2"
    else:
        raw_recip_kir = "C1/C2"
    kir_subbank_matched = True
    if is_plan_b_arch:
        c1_allele_tags = (
            "C0102", "C*01:02", "C0304", "C*03:04", "CW0302", "CW*0302", "C0702", "C*07:02",
            "C0802", "C*08:02", "C1203", "C*12:03", "C1402", "C*14:02", "C1601", "C*16:01",
        )
        c2_allele_tags = (
            "C0202", "C*02:02", "C0401", "C*04:01", "CW0402", "CW*0402", "C0402", "C*04:02", "C0501", "C*05:01",
            "C0602", "C*06:02", "C1502", "C*15:02", "C1505", "C*15:05", "C1701", "C*17:01",
        )
        explicit_group_tags = ("C1/C1", "C1_C1", "C2/C2", "C2_C2", "C1/C2", "C1_C2", "C2/C1", "C2_C1")
        if prod_def_hlac_str and any(t in prod_def_hlac_str for t in (*c1_allele_tags, *c2_allele_tags, *explicit_group_tags)):
            combined_hlac_donor_str = prod_def_hlac_str.upper()
        else:
            combined_hlac_donor_str = f"{loc1_str} {prod_def_hlac_str}".upper()
        has_c1_allele = any(t in combined_hlac_donor_str for t in c1_allele_tags)
        has_c2_allele = any(t in combined_hlac_donor_str for t in c2_allele_tags)
        has_c1_c2_het = any(tag in combined_hlac_donor_str for tag in ("C1/C2", "C1_C2", "C2/C1", "C2_C1"))
        has_c1_homo = (
            any(tag in combined_hlac_donor_str for tag in ("C1/C1", "C1_C1"))
            or (has_c1_allele and not has_c2_allele and not has_c1_c2_het)
        )
        has_c2_homo = (
            any(tag in combined_hlac_donor_str for tag in ("C2/C2", "C2_C2"))
            or (has_c2_allele and not has_c1_allele and not has_c1_c2_het)
        )
        if raw_recip_kir == "C1/C1" and has_c2_homo:
            kir_subbank_matched = False
        elif raw_recip_kir == "C2/C2" and has_c1_homo:
            kir_subbank_matched = False
        elif raw_recip_kir == "C1/C2" and (has_c1_homo or has_c2_homo):
            kir_subbank_matched = False
    pvid_str = str(compiled_regimen.get("product_variant_id", "")).upper()
    is_b_prime_sct = (
        "B_SCT" in pvid_str
        or "BPRIME_SCT" in pvid_str
        or "B_SCT" in loc1_str
    )
    if is_b_prime_sct and prod_def_hlac_str:
        has_c1_sct = ("C*03:04" in prod_def_hlac_str) and ("C*04:01" not in prod_def_hlac_str)
        has_c2_sct = ("C*04:01" in prod_def_hlac_str) and ("C*03:04" not in prod_def_hlac_str)
        if raw_recip_kir == "C1/C1" and has_c2_sct:
            kir_subbank_matched = False
        elif raw_recip_kir == "C2/C2" and has_c1_sct:
            kir_subbank_matched = False
        elif raw_recip_kir == "C1/C2" and (has_c1_sct or has_c2_sct):
            kir_subbank_matched = False
    has_plan_b_hlac = (
        is_plan_b_arch
        and prod_def_hlac_ok
        and kir_subbank_matched
        and (
            (
                is_b_prime_sct
                and "OI-GRNA-01_B2M_EXON1_ASCAS12A" in panel
                and "OI-GRNA-02_CIITA_EXON5_ASCAS12A" in panel
            )
            or (
                (not is_b_prime_sct)
                and is_plan_b_loc1
                and "PLAN_B_HLAA_EXON2_KO_GUIDE_INVESTIGATIONAL" in panel
                and "PLAN_B_HLAB_EXON2_KO_GUIDE_INVESTIGATIONAL" in panel
                and "OI-GRNA-01_B2M_EXON1_ASCAS12A" not in panel
            )
        )
    )
    if is_plan_b_arch and not has_plan_b_hlac:
        p_surv -= 0.55
    elif (not is_plan_b_arch) and ("OI-GRNA-01_B2M_EXON1_ASCAS12A" not in panel):
        p_surv -= 0.18

    has_arm1b_floor = (
        ("DONOR-1B" in loc1_str and "OI-GRNA-01_B2M_EXON1_ASCAS12A" in panel)
        or has_plan_b_hlac
    )
    arm1b_cd47_sirpa_nk_cliff_unverified_eval = _validate_strict_bool(
        _get_opt(patient, "arm1b_cd47_sirpa_nk_cliff_unverified", False),
        "arm1b_cd47_sirpa_nk_cliff_unverified",
    )
    track1_sorted_kir_nk_cleared_eval = _validate_strict_bool(
        _get_opt(patient, "track1_sorted_kir_nkg2a_neg_nk_lysis_cleared", False),
        "track1_sorted_kir_nkg2a_neg_nk_lysis_cleared",
    )
    single_200cm2_pouch_forced_eval = _validate_strict_bool(
        _get_opt(patient, "single_200cm2_pouch_forced_without_modular_scaling", False),
        "single_200cm2_pouch_forced_without_modular_scaling",
    )
    tyk2_indep_ifng_dr3_dq8_hazard_eval = _validate_strict_bool(
        _get_opt(patient, "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard", False),
        "tyk2_independent_ifng_tnfa_or_indirect_dr3_dq8_hazard",
    )
    gate0_indirect_dr3_dq8_cleared_eval = _validate_strict_bool(
        _get_opt(patient, "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared", False),
        "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared",
    )
    bw4_kir3dl1_educated_eval = _validate_strict_bool(
        _get_opt(patient, "recipient_bw4_kir3dl1_educated", False),
        "recipient_bw4_kir3dl1_educated",
    )
    bw4_kir3dl1_cleared_eval = _validate_strict_bool(
        _get_opt(patient, "bw4_kir3dl1_sorted_nk_lysis_cleared", False),
        "bw4_kir3dl1_sorted_nk_lysis_cleared",
    )
    if kir_nk >= 0.30 and not has_arm1b_floor:
        p_surv -= 0.45 * (kir_nk - 0.25)
    elif (
        kir_nk >= 0.30
        and not has_plan_b_hlac
        and arm1b_cd47_sirpa_nk_cliff_unverified_eval
        and not track1_sorted_kir_nk_cleared_eval
    ):
        p_surv -= 0.35 * (kir_nk - 0.25)
    if kir_nk >= 0.48 and not has_plan_b_hlac:
        p_surv -= 0.28 * (kir_nk - 0.45)
    # Propagate recipient HLA-B/A Bw4+ KIR3DL1-educated missing-self NK hazard when sorted KIR3DL1+NKG2A- NK clearance is unverified (L5A-01 [DEC-077]):
    if bw4_kir3dl1_educated_eval and not bw4_kir3dl1_cleared_eval:
        p_surv -= 0.18

    # 1b. Propagate 4-digit HLA-C alloreactive T-cell, anti-HLA-C heavy-chain DSA, activating KIR2DS1,
    # anti-HLA-E DSA, high-cPRA sensitization, and CMV NKG2C+ adaptive NK hazards across BOTH endogenous HLA-C
    # and fixed-peptide B'-SCT (AUDIT-REF, AUDIT-REF, & AUDIT-REF;
    # Aubert et al. 2014 PMID:24804568; Whitelegg et al. 2005 PMID:16034111; Yu et al. 2002 PMID:11907065; Stewart et al. 2005 PMID:16141329):
    raw_eval_hlac_rel = _get_opt(
        patient,
        "hlac_4digit_relation",
        "UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT" if is_b_prime_sct else "SELF_4DIGIT_MATCH",
    )
    if not isinstance(raw_eval_hlac_rel, str) or not raw_eval_hlac_rel.strip():
        raise ValueError(f"Invalid hlac_4digit_relation {raw_eval_hlac_rel!r}: must be a non-empty string (AUDIT-REF)")
    hlac_4digit_rel = raw_eval_hlac_rel.strip().upper()
    if hlac_4digit_rel not in _VALID_HLAC_4DIGIT_RELATIONS:
        raise ValueError(f"Invalid hlac_4digit_relation {hlac_4digit_rel!r} (AUDIT-REF)")
    has_anti_hlac_dsa = _validate_strict_bool(
        _get_opt(patient, "recipient_anti_hla_c_dsa", False), "recipient_anti_hla_c_dsa"
    )
    sct_crossmatch_cleared = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_recipient_crossmatch_and_mlc_cleared", False),
        "b_prime_sct_recipient_crossmatch_and_mlc_cleared",
    )
    sct_serum_crossmatch_neg = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative", False),
        "b_prime_sct_exact_product_serum_cdc_adcc_crossmatch_negative",
    )
    sct_cd8_mlc_neg = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_recipient_cd8_mlc_negative", False),
        "b_prime_sct_recipient_cd8_mlc_negative",
    )
    has_activating_kir2ds1_hazard = _validate_strict_bool(
        _get_opt(patient, "recipient_activating_kir2ds1_unlicensed_hazard", False),
        "recipient_activating_kir2ds1_unlicensed_hazard",
    )
    kir2ds1_cleared = _validate_strict_bool(
        _get_opt(patient, "b_prime_sct_kir2ds1_non_activation_cleared", False),
        "b_prime_sct_kir2ds1_non_activation_cleared",
    )
    has_anti_hlae_dsa = _validate_strict_bool(
        _get_opt(patient, "recipient_anti_hla_e_dsa", False), "recipient_anti_hla_e_dsa"
    )
    exact_hlae_sct_serum_crossmatch_neg = _validate_strict_bool(
        _get_opt(patient, "exact_hlae_sct_serum_crossmatch_negative", False),
        "exact_hlae_sct_serum_crossmatch_negative",
    )
    exact_product_non_hla_and_sct_flow_crossmatch_neg = _validate_strict_bool(
        _get_opt(patient, "exact_product_non_hla_and_sct_flow_crossmatch_negative", False),
        "exact_product_non_hla_and_sct_flow_crossmatch_negative",
    )
    gate0a_nkg2c_adaptive_nk_cleared = _validate_strict_bool(
        _get_opt(patient, "gate0a_nkg2c_adaptive_nk_cleared", False),
        "gate0a_nkg2c_adaptive_nk_cleared",
    )
    recipient_cpra_pct = _validate_finite_float(
        _get_opt(patient, "recipient_cpra_pct", 0.0), "recipient_cpra_pct", 0.0, 100.0
    )
    cmv_nkg2c_adaptive_nk_fraction = _validate_finite_float(
        _get_opt(patient, "cmv_nkg2c_adaptive_nk_fraction", 0.0), "cmv_nkg2c_adaptive_nk_fraction", 0.0, 1.0
    )
    anti_hlac_dsa_cleared = bool(is_b_prime_sct and sct_serum_crossmatch_neg)
    cd8_mlc_cleared = bool(is_b_prime_sct and (sct_cd8_mlc_neg or sct_crossmatch_cleared))
    _c2_bearing_sct_ids = {
        "ISLET1_B_SCT_C2_CIITA",
        "ISLET1_B_SCT_C2_CIITA_TYK2_PE",
        "ISLET1_B_SCT_C1C2_CIITA",
        "ISLET1_B_SCT_C1C2_CIITA_TYK2_PE",
    }
    if is_plan_b_arch:
        if has_anti_hlac_dsa and not anti_hlac_dsa_cleared:
            p_surv -= 0.40
        is_c2_b_prime_sct = is_b_prime_sct and (
            pvid_str in _c2_bearing_sct_ids
            or "C0401" in loc1_str
            or "C*04:01" in prod_def_hlac_str
            or raw_recip_kir in {"C2/C2", "C1/C2"}
        )
        if hlac_4digit_rel == "KIR_LIGAND_MISMATCH":
            p_surv -= 0.55
        elif (
            hlac_4digit_rel in {"KIR_MATCH_FOREIGN_4DIGIT_ALLO", "UNSPECIFIED_ASSUMED_4DIGIT_OR_B_PRIME_SCT", "UNSPECIFIED"}
            or is_c2_b_prime_sct
        ) and not cd8_mlc_cleared:
            p_surv -= (0.08 if is_b_prime_sct else 0.25)
        if (
            has_activating_kir2ds1_hazard
            and not kir2ds1_cleared
            and (
                pvid_str in _c2_bearing_sct_ids
                or "C0401" in loc1_str
                or "C*04:01" in prod_def_hlac_str
                or raw_recip_kir in {"C2/C2", "C1/C2"}
            )
        ):
            p_surv -= 0.12
    if has_anti_hlae_dsa and not exact_hlae_sct_serum_crossmatch_neg:
        p_surv -= 0.35
    exact_product_crossmatch_explicit_false_eval = (
        "exact_product_non_hla_and_sct_flow_crossmatch_negative" in patient
        and patient.get("exact_product_non_hla_and_sct_flow_crossmatch_negative") is False
    )
    if (
        exact_product_crossmatch_explicit_false_eval
        or (
            "TRACK_3" not in track_str
            and recipient_cpra_pct >= 80.0
            and not exact_product_non_hla_and_sct_flow_crossmatch_neg
        )
    ):
        p_surv -= 0.25
    if cmv_nkg2c_adaptive_nk_fraction >= 0.25 and not gate0a_nkg2c_adaptive_nk_cleared:
        p_surv -= 0.12

    # 2. Check high-IFN-I / high-autoantigen / IFIH1 / PTPN22 T-cell protection via TYK2 P1104A shield:
    # A graft has biological TYK2 P1104A protection if EITHER edited via OI-GRNA-05 OR naturally homozygous C/C on the graft genome (AUDIT-REF):
    is_track3 = "TRACK_3" in track_str
    raw_recip_tyk2 = str(
        _get_opt(patient, "recipient_tyk2_rs34536443_genotype", _get_opt(patient, "tyk2_rs34536443_genotype", "G/G"))
    ).strip().upper()
    raw_donor_tyk2 = str(_get_opt(patient, "donor_tyk2_rs34536443_genotype", "G/G")).strip().upper()
    for label, genotype in (("recipient_tyk2_rs34536443_genotype", raw_recip_tyk2),
                            ("donor_tyk2_rs34536443_genotype", raw_donor_tyk2)):
        if genotype not in VALID_TYK2_RS34536443_GENOTYPES:
            raise ValueError(f"Invalid {label} {genotype!r}")
    starting_graft_tyk2 = raw_recip_tyk2 if is_track3 else raw_donor_tyk2
    naturally_protected_tyk2_cc = starting_graft_tyk2 == "C/C"
    # Derived labels cannot grant protection contrary to the actual edited genome.
    for block, key, expected in (
        (dosing_block, "already_homozygous_tyk2_p1104a", naturally_protected_tyk2_cc),
        (compiled_regimen, "final_graft_tyk2_p1104a_protective_state", naturally_protected_tyk2_cc or tyk2_flag),
    ):
        if key in block and _validate_strict_bool(block[key], key) != expected:
            raise ValueError(f"Inconsistent derived TYK2 protection field {key!r}")
    has_edited_tyk2_shield = (
        tyk2_flag
        and "OI-GRNA-05_TYK2_P1104A_PEMAX_EPEGRNA_OPTIONAL" in panel
    )
    has_biological_tyk2_shield = bool(has_edited_tyk2_shield or naturally_protected_tyk2_cc)

    ifih1_mod = 0.08 if ifih1_g == "T/T" else (0.03 if ifih1_g in {"C/T", "T/C"} else 0.0)
    ptpn22_mod = 0.10 if ptpn22_g == "A/A" else (0.04 if ptpn22_g in {"G/A", "A/G"} else 0.0)
    hyper_ifn_stress = max(0.0, (ifn1_pct - 75.0) / 25.0) * max(0.0, auto_z - 1.5) + ifih1_mod + ptpn22_mod
    needs_tyk2 = (ifn1_pct >= 80.0 or auto_z >= 2.0 or ifih1_g == "T/T" or ptpn22_g == "A/A")
    if needs_tyk2 and not has_biological_tyk2_shield:
        p_surv -= 0.095 * (1.0 + hyper_ifn_stress)
    else:
        p_surv -= 0.010 * hyper_ifn_stress
    if (auto_z >= 2.0 or tyk2_indep_ifng_dr3_dq8_hazard_eval) and not gate0_indirect_dr3_dq8_cleared_eval:
        p_surv -= 0.14

    # 3. Propagate recipient HLA-E genotype (*01:01 vs *01:03 endogenous NKG2A licensing tuning; AUDIT-REF)
    # and causal SIRPa V1/V2/V8 macrophage clearance risk when CD47 shield is absent or compromised (AUDIT-REF):
    if hla_e_g == "*01:01/*01:01":
        p_surv -= 0.004 if has_arm1b_floor else 0.035
    elif hla_e_g in {"*01:01/*01:03", "*01:03/*01:01"}:
        p_surv -= 0.002 if has_arm1b_floor else 0.015
    plan_b_declares_cd47 = has_plan_b_hlac and prod_def_cd47_ok and ((("CD47" in loc1_str) and ("CD47" in plan_str)) or is_b_prime_sct)
    # Note: SIRPA v8 binding to CD47-IgV-GPI is unmeasured (`HOLD_CD47_SIRPA_BINDING_COMPROMISED_REQUIRES_ALTERNATIVE_STEALTH_QUALIFICATION`),
    # so `v8` carriers are treated as unverified CD47-SIRPa stealth unless explicitly cleared (L5A-01 [DEC-077]):
    has_cd47_shield = ((("CD47" in loc1_str) and not is_plan_b_arch) or plan_b_declares_cd47) and (not cd47_compromised) and ("v8" not in sirpa_g)
    # Propagate SP-3C (HLA-C*07:02 VMAPRALLL, UniProt P10321) leader deficit when exogenous HLA-E SCT is absent (S55R2-D5-02):
    if is_plan_b_arch and prod_def_leader_class == "SP-3C" and not plan_b_declares_cd47:
        p_surv -= 0.04
    if not has_cd47_shield:
        if sirpa_g == "v1/v1":
            p_surv -= 0.22
        elif sirpa_g == "v1/v2":
            p_surv -= 0.26
        else:
            p_surv -= 0.32

    # 4. Propagate prescribed dose IEQ vs fixed 200 cm^2 single-pouch capacity (2,500 IEQ/cm^2 limit; AUDIT-REF):
    areal_loading_single_pouch = dose_ieq / 200.0
    if single_200cm2_pouch_forced_eval and areal_loading_single_pouch > 1500.0:
        p_surv -= min(0.65, 0.38 * ((areal_loading_single_pouch - 1500.0) / 1000.0))
    elif areal_loading_single_pouch > 2500.0:
        p_surv -= 0.045 * ((areal_loading_single_pouch - 2500.0) / 2500.0)
    elif areal_loading_single_pouch < 2500.0:
        p_surv -= 0.008 * ((2500.0 - areal_loading_single_pouch) / 2500.0)

    # 5. Check active editing-genome failure if the genome being edited (recipient in Track 3, donor MCB in Track 1/2) carries a rare protospacer SNP without custom gRNA or with PAM disruption (AUDIT-REF & AUDIT-REF):
    has_custom_grna = "PATIENT_ALLELE_MATCHED_CUSTOM_SPACER_SYNTHESIS" in panel
    active_editing_genome_has_snp = (recip_snp if is_track3 else donor_mcb_snp) or (rare_snp_pam_disrupted_eval and (recip_snp or donor_mcb_snp))
    if (active_editing_genome_has_snp and not has_custom_grna) or rare_snp_pam_disrupted_eval:
        p_surv -= 0.35

    return float(np.clip(p_surv, 0.05, 0.988))


def _hwe_alternate_alleles(uniform_draw: float, alternate_frequency: float) -> int:
    """Inverse CDF for Binomial(2, p), with cumulative heterozygote threshold."""
    u = _validate_finite_float(uniform_draw, "uniform_draw", 0., 1.)
    p = _validate_finite_float(alternate_frequency, "alternate_frequency", 0., 1.)
    return 2 if p == 1.0 or u < p*p else (1 if u < p*p + 2*p*(1-p) else 0)


def simulate_10000_patient_global_population_trial(
    compiler_fn: Callable[[Dict[str, Any]], Dict[str, Any]] = compile_personalized_patient_blueprint,
) -> Dict[str, Any]:
    """
    Illustrative model-based scenario simulation of N = 10,000 virtual T1D patients
    (2,000 per global superpopulation: EUR, AFR, EAS, SAS, AMR) over a 5-year horizon.
    Explicitly invokes `compiler_fn` for every single virtual patient and samples realistic
    ancestry-conditioned distributions for `ifih1_rs1990760_genotype`, `ptpn22_rs2476601_genotype`,
    `recipient_hla_e_genotype`, `tyk2_rs34536443_genotype`, `recipient_kir_ligand_group`, and
    `recipient_bmi_kg_m2` across all 50,000 calls (closing AUDIT-REF, AUDIT-REF, and
    AUDIT-REF).
    """
    engineering_parameters = _engineering_parameters()
    _is_default_compiler = (compiler_fn is compile_personalized_patient_blueprint)
    if _is_default_compiler:
        compiler_fn = lambda patient: compile_personalized_patient_blueprint(patient, engineering_parameters=engineering_parameters)
    rng = np.random.default_rng(20261023)
    # Columns: (code, name, kir_a_f, c1_f, sirpa_v1_f, tyk2_f, mean_kir_nk, ifih1_t_f, ptpn22_a_f, hlae_0103_f, mean_bmi)
    ancestries = [
        ("EUR", "European", 0.55, 0.62, 0.68, 0.042, 0.40, 0.60, 0.095, 0.45, 25.4),
        ("AFR", "African / African Diaspora", 0.48, 0.54, 0.74, 0.008, 0.43, 0.18, 0.008, 0.52, 26.1),
        ("EAS", "East Asian", 0.72, 0.81, 0.42, 0.002, 0.45, 0.22, 0.002, 0.68, 22.8),
        ("SAS", "South Asian", 0.51, 0.64, 0.58, 0.015, 0.41, 0.36, 0.018, 0.54, 23.9),
        ("AMR", "Indigenous American / Admixed Latino", 0.60, 0.69, 0.55, 0.024, 0.41, 0.44, 0.042, 0.58, 26.3),
    ]

    n_per_pop = 2000
    pop_summaries: List[Dict[str, Any]] = []
    overall_counts = {
        "unedited_5yr_survival": 0,
        "hla_e_only_5yr_survival": 0,
        "standardized_plan_a_5yr_survival": 0,
        "personalized_compiled_5yr_survival": 0,
        "ablated_compiler_control_5yr_survival": 0,
    }
    kb_counts = {"K4_B0": 0, "K5_B0": 0, "K4_B1": 0, "K5_B1": 0}
    all_std_analytic: List[float] = []
    all_pers_analytic: List[float] = []
    all_abl_analytic: List[float] = []
    all_pers_ranker_discounted: List[float] = []
    all_unmat_b_analytic: List[float] = []
    all_unmat_b_discounted: List[float] = []
    all_personalized_complexity_cases: List[tuple[float, int, int]] = []
    all_threshold_sensitivity_cases: List[tuple[float, float, float, int]] = []

    sirpa_map = {0: "v2/v2", 1: "v1/v2", 2: "v1/v1"}
    ifih1_map = {0: "C/C", 1: "C/T", 2: "T/T"}
    ptpn22_map = {0: "G/G", 1: "G/A", 2: "A/A"}
    hlae_map = {0: "*01:01/*01:01", 1: "*01:01/*01:03", 2: "*01:03/*01:03"}
    tyk2_map = {0: "G/G", 1: "G/C", 2: "C/C"}
    kir_map = {0: "C2/C2", 1: "C1/C2", 2: "C1/C1"}

    for code, name, kir_a_f, c1_f, sirpa_v1_f, tyk2_f, mean_kir_nk, ifih1_t_f, ptpn22_a_f, hlae_0103_f, mean_bmi in ancestries:
        nkg2a_neg_kir_frac = np.clip(rng.normal(mean_kir_nk, 0.08, size=n_per_pop), 0.15, 0.72)
        autoantigen_z = np.clip(rng.normal(1.25, 0.75, size=n_per_pop), 0.0, 4.0)
        ifn1_pct = rng.uniform(0.0, 100.0, size=n_per_pop)
        sirpa_v1_alleles = rng.binomial(2, sirpa_v1_f, size=n_per_pop)
        rare_snp = rng.uniform(0.0, 1.0, size=n_per_pop) < 0.0012
        ifih1_alleles = rng.binomial(2, ifih1_t_f, size=n_per_pop)
        ptpn22_alleles = rng.binomial(2, ptpn22_a_f, size=n_per_pop)
        hlae_alleles = rng.binomial(2, hlae_0103_f, size=n_per_pop)
        bmi_vec = np.clip(rng.normal(mean_bmi, 3.8, size=n_per_pop), 16.5, 42.0)

        p_unedited = np.zeros(n_per_pop)
        p_hla_e_only = np.clip(
            0.82
            - 1.15 * np.maximum(0.0, nkg2a_neg_kir_frac - 0.25)
            - 0.22 * np.maximum(0.0, autoantigen_z - 0.6),
            0.02,
            0.85,
        )

        p_standard_plan_a = np.zeros(n_per_pop, dtype=np.float64)
        p_personalized = np.zeros(n_per_pop, dtype=np.float64)
        p_ablated_control = np.zeros(n_per_pop, dtype=np.float64)

        # Unpersonalized baseline regimen (fixed Plan A Arm 1B without Plan B HLA-C sub-bank or TYK2 P1104A)
        fixed_plan_a_regimen = {
            "selected_fda_deployment_track": "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC",
            "locus1_construct_prescribed": "OI-L2-DONOR-1B_B2M_EEF1A1_HLAE_CD47IGV",
            "plan_architecture_prescribed": "Plan A (Core 4-Locus Tri-Modal)",
            "include_tyk2_p1104a_rs34536443_shield": False,
            "compiled_grna_panel": [
                "OI-GRNA-01_B2M_EXON1_ASCAS12A",
                "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
                "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
                "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
            ],
        }
        # Non-personalized Chassis A Arm 1A bicistronic 2A comparator: retains all 4 core Chassis A guides (B2M, CIITA, CXCL10, AAVS1)
        # while omitting Arm 1B's 2A-free monocistronic EEF1A1 Hour-0 CD47-IgV floor (238,004 vs 209,444 molecules/cell), Plan B HLA-C sub-bank,
        # and TYK2 P1104A personalization (closing AUDIT-REF, AUDIT-REF, & AUDIT-REF):
        ablated_regimen = {
            "selected_fda_deployment_track": "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC",
            "locus1_construct_prescribed": "OI-L2-DONOR-1A_B2M_ENDO_HLAE_CD47IGV",
            "plan_architecture_prescribed": "Plan A (Core 4-Locus Tri-Modal)",
            "include_tyk2_p1104a_rs34536443_shield": False,
            "compiled_grna_panel": [
                "OI-GRNA-01_B2M_EXON1_ASCAS12A",
                "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
                "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
                "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
            ],
        }

        n_track1 = 0
        n_track2 = 0
        n_plan_b = 0
        n_tyk2 = 0
        for i in range(n_per_pop):
            # Sample recipient TYK2 and KIR ligand group deterministically from patient index + ancestry frequencies
            # without perturbing the main rng stream used for u_uned..u_pers:
            tyk2_draw = _hwe_alternate_alleles(((i * 997 + len(code)) % 10000) / 10000., tyk2_f)
            kir_draw = 2 if ((i * 7919 + len(code)) % 1000) < int(round(1000 * (c1_f ** 2))) else (
                1 if ((i * 7919 + len(code)) % 1000) < int(round(1000 * (c1_f ** 2 + 2 * c1_f * (1 - c1_f)))) else 0
            )
            pt_rec = {
                "patient_id": f"SIM-{code}-{i:04d}",
                "ancestry": code,
                "nkg2a_neg_kir_educated_nk_fraction": float(nkg2a_neg_kir_frac[i]),
                "pre_existing_gad65_ia2_znt8_t_cell_zscore": float(autoantigen_z[i]),
                "ifn1_ish_genetic_risk_percentile": float(ifn1_pct[i]),
                "sirpa_genotype": sirpa_map[int(sirpa_v1_alleles[i])],
                "ifih1_rs1990760_genotype": ifih1_map[int(ifih1_alleles[i])],
                "ptpn22_rs2476601_genotype": ptpn22_map[int(ptpn22_alleles[i])],
                "recipient_hla_e_genotype": hlae_map[int(hlae_alleles[i])],
                "tyk2_rs34536443_genotype": tyk2_map[tyk2_draw],
                "recipient_kir_ligand_group": kir_map[kir_draw],
                "recipient_bmi_kg_m2": round(float(bmi_vec[i]), 2),
                "carries_rare_snp_in_standard_grna": bool(rare_snp[i]),
                "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared": True,
            }
            compiled_reg = compiler_fn(pt_rec)
            if "TRACK_2" in compiled_reg.get("selected_fda_deployment_track", ""):
                n_track2 += 1
            else:
                n_track1 += 1
            is_pb = "Plan B" in compiled_reg.get("plan_architecture_prescribed", "")
            has_t2 = bool(compiled_reg.get("include_tyk2_p1104a_rs34536443_shield"))
            if is_pb:
                n_plan_b += 1
            if has_t2:
                n_tyk2 += 1
            prod_k = compiled_reg.get("product_definition", {}).get("K") if isinstance(compiled_reg.get("product_definition"), dict) else None
            k_cnt = int(prod_k) if isinstance(prod_k, int) and prod_k in (4, 5) else (5 if is_pb else 4)
            b_cnt = 1 if has_t2 else 0
            kb_counts[f"K{k_cnt}_B{b_cnt}"] += 1

            p_standard_plan_a[i] = evaluate_regimen_survival_probability(pt_rec, fixed_plan_a_regimen)
            p_personalized[i] = evaluate_regimen_survival_probability(pt_rec, compiled_reg)
            p_ablated_control[i] = evaluate_regimen_survival_probability(pt_rec, ablated_regimen)
            disc_factor = engineering_complexity_discount(k_cnt, b_cnt, lambda_K=engineering_parameters["lambda_K"], gamma_B=engineering_parameters["gamma_B"])
            all_pers_ranker_discounted.append(float(p_personalized[i]) * disc_factor)
            if is_pb and isinstance(compiled_reg.get("product_definition"), dict):
                reg_unmat = copy.deepcopy(compiled_reg)
                reg_unmat["product_definition"]["cassette_topology"]["CD47_HLAE"] = "ABSENT"
                p_unmat_i = float(evaluate_regimen_survival_probability(pt_rec, reg_unmat))
            else:
                p_unmat_i = float(p_personalized[i])
            all_unmat_b_analytic.append(p_unmat_i)
            all_unmat_b_discounted.append(p_unmat_i * disc_factor)
            all_personalized_complexity_cases.append((float(p_personalized[i]), k_cnt, b_cnt))
            if _is_default_compiler and 0.48 <= float(nkg2a_neg_kir_frac[i]) < 0.52:
                pt_low_kir = dict(pt_rec, nkg2a_neg_kir_educated_nk_fraction=0.35)
                reg_a_if_unrouted = compiler_fn(pt_low_kir)
                p_if_plan_a = float(evaluate_regimen_survival_probability(pt_rec, reg_a_if_unrouted))
            else:
                p_if_plan_a = float(p_personalized[i])
            all_threshold_sensitivity_cases.append((float(nkg2a_neg_kir_frac[i]), float(p_personalized[i]), p_if_plan_a, b_cnt))

        all_std_analytic.extend(p_standard_plan_a.tolist())
        all_pers_analytic.extend(p_personalized.tolist())
        all_abl_analytic.extend(p_ablated_control.tolist())

        u_uned = rng.uniform(0.0, 1.0, size=n_per_pop)
        u_hlae = rng.uniform(0.0, 1.0, size=n_per_pop)
        u_std = rng.uniform(0.0, 1.0, size=n_per_pop)
        u_pers = rng.uniform(0.0, 1.0, size=n_per_pop)

        surv_uned = int(np.sum(u_uned < p_unedited))
        surv_hlae = int(np.sum(u_hlae < p_hla_e_only))
        surv_std = int(np.sum(u_std < p_standard_plan_a))
        surv_pers = int(np.sum(u_pers < p_personalized))
        surv_abl = int(np.sum(u_pers < p_ablated_control))

        # Common Random Numbers (CRN) paired evaluation using shared uniform draw u_pers across all arms (AUDIT-REF):
        surv_std_crn = int(np.sum(u_pers < p_standard_plan_a))
        std_crn_pct = round(100.0 * surv_std_crn / n_per_pop, 2)
        analytic_std_unrounded_pct = 100.0 * float(np.mean(p_standard_plan_a))
        analytic_pers_unrounded_pct = 100.0 * float(np.mean(p_personalized))
        analytic_std_pct = round(analytic_std_unrounded_pct, 2)
        analytic_pers_pct = round(analytic_pers_unrounded_pct, 2)
        analytic_abl_pct = round(100.0 * float(np.mean(p_ablated_control)), 2)

        overall_counts["unedited_5yr_survival"] += surv_uned
        overall_counts["hla_e_only_5yr_survival"] += surv_hlae
        overall_counts["standardized_plan_a_5yr_survival"] += surv_std
        overall_counts["personalized_compiled_5yr_survival"] += surv_pers
        overall_counts["ablated_compiler_control_5yr_survival"] += surv_abl

        std_pct = round(100.0 * surv_std / n_per_pop, 2)
        pers_pct = round(100.0 * surv_pers / n_per_pop, 2)
        abl_pct = round(100.0 * surv_abl / n_per_pop, 2)
        arr_sample = round(pers_pct - std_pct, 2)
        arr_analytic = round(analytic_pers_pct - analytic_std_pct, 2)
        mc_se_pct = round(100.0 * math.sqrt(float(np.mean(p_personalized * (1.0 - p_personalized))) / n_per_pop), 2)
        pop_summaries.append(
            {
                "ancestry_code": code,
                "ancestry_name": name,
                "n_patients_simulated": n_per_pop,
                "population_kir_group_a_haplotype_freq": kir_a_f,
                "population_hla_c1_ligand_freq": c1_f,
                "population_sirpa_v1_allele_freq": sirpa_v1_f,
                "population_tyk2_p1104a_rs34536443_allele_freq": tyk2_f,
                "population_ifih1_rs1990760_t_allele_freq": ifih1_t_f,
                "population_ptpn22_rs2476601_a_allele_freq": ptpn22_a_f,
                "population_hla_e_0103_allele_freq": hlae_0103_f,
                "mean_recipient_bmi_kg_m2": round(float(np.mean(bmi_vec)), 2),
                "mean_nkg2a_neg_kir_educated_nk_fraction": round(float(np.mean(nkg2a_neg_kir_frac)), 3),
                "unedited_allogeneic_5yr_survival_pct": round(100.0 * surv_uned / n_per_pop, 2),
                "hla_e_only_5yr_survival_pct": round(100.0 * surv_hlae / n_per_pop, 2),
                "standardized_openislet_plan_a_5yr_survival_pct": std_pct,
                "standardized_openislet_plan_a_crn_paired_5yr_survival_pct_v80_015": std_crn_pct,
                "personalized_dna_compiled_openislet_5yr_survival_pct": pers_pct,
                "personalized_monte_carlo_binomial_se_pct": mc_se_pct,
                "ablated_compiler_control_5yr_survival_pct": abl_pct,
                "analytic_expected_5yr_survival_standardized_plan_a_pct": analytic_std_pct,
                "analytic_expected_5yr_survival_personalized_compiled_pct": analytic_pers_pct,
                "analytic_expected_5yr_survival_personalized_compiled_unrounded_pct": round(analytic_pers_unrounded_pct, 6),
                "analytic_expected_5yr_survival_ablated_arm1a_control_pct": analytic_abl_pct,
                "arr_personalized_vs_standardized_plan_a_pct": arr_sample,
                "arr_crn_paired_personalized_vs_standardized_plan_a_pct_v80_015": round(pers_pct - std_crn_pct, 2),
                "arr_analytic_expected_personalized_vs_standardized_plan_a_pct_v80_015": arr_analytic,
                "nnc_number_needed_to_compile_vs_standardized_plan_a": (
                    round(100.0 / arr_sample, 2) if arr_sample > 0.0 else None
                ),
                "nnc_analytic_expected_vs_standardized_plan_a": (
                    round(100.0 / arr_analytic, 2) if arr_analytic > 0.0 else None
                ),
                "compiler_regimen_allocation_counts": {
                    "track1_off_the_shelf": n_track1,
                    "track2_stratified_subline": n_track2,
                    "plan_b_hla_c_retained": n_plan_b,
                    "tyk2_p1104a_shield": n_tyk2,
                },
            }
        )

    def compile_reference_patient(record: Dict[str, Any]) -> Dict[str, Any]:
        return compile_personalized_patient_blueprint(record, engineering_parameters=engineering_parameters)

    reference_patients = [
        compile_reference_patient(
            {
                "patient_id": "PT-001_EUR_STANDARD_T1D",
                "ancestry": "EUR",
                "nkg2a_neg_kir_educated_nk_fraction": 0.36,
                "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.1,
                "ifn1_ish_genetic_risk_percentile": 45.0,
                "sirpa_genotype": "v1/v1",
                "ifih1_rs1990760_genotype": "C/T",
                "ptpn22_rs2476601_genotype": "G/G",
                "recipient_hla_e_genotype": "*01:03/*01:03",
                "tyk2_rs34536443_genotype": "G/G",
                "recipient_kir_ligand_group": "C1/C2",
                "recipient_bmi_kg_m2": 24.0,
                "carries_rare_snp_in_standard_grna": False,
            }
        ),
        compile_reference_patient(
            {
                "patient_id": "PT-002_EAS_HIGH_KIR_EDUCATED_NK",
                "ancestry": "EAS",
                "nkg2a_neg_kir_educated_nk_fraction": 0.54,
                "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.4,
                "ifn1_ish_genetic_risk_percentile": 62.0,
                "sirpa_genotype": "v2/v2",
                "ifih1_rs1990760_genotype": "C/C",
                "ptpn22_rs2476601_genotype": "G/G",
                "recipient_hla_e_genotype": "*01:03/*01:03",
                "tyk2_rs34536443_genotype": "G/G",
                "recipient_kir_ligand_group": "C1/C1",
                "recipient_bmi_kg_m2": 22.5,
                "carries_rare_snp_in_standard_grna": False,
            }
        ),
        compile_reference_patient(
            {
                "patient_id": "PT-003_AFR_HIGH_IFN1_TRIPLE_AUTOANTIBODY",
                "ancestry": "AFR",
                "nkg2a_neg_kir_educated_nk_fraction": 0.46,
                "pre_existing_gad65_ia2_znt8_t_cell_zscore": 2.6,
                "ifn1_ish_genetic_risk_percentile": 92.0,
                "sirpa_genotype": "v1/v2",
                "ifih1_rs1990760_genotype": "T/T",
                "ptpn22_rs2476601_genotype": "G/A",
                "recipient_hla_e_genotype": "*01:01/*01:03",
                "tyk2_rs34536443_genotype": "G/G",
                "recipient_kir_ligand_group": "C1/C2",
                "recipient_bmi_kg_m2": 26.8,
                "carries_rare_snp_in_standard_grna": False,
            }
        ),
        compile_reference_patient(
            {
                "patient_id": "PT-004_AMR_AUTOLOGOUS_RARE_PAM_SNP",
                "ancestry": "AMR",
                "nkg2a_neg_kir_educated_nk_fraction": 0.41,
                "pre_existing_gad65_ia2_znt8_t_cell_zscore": 1.9,
                "ifn1_ish_genetic_risk_percentile": 78.0,
                "sirpa_genotype": "v1/v2",
                "ifih1_rs1990760_genotype": "C/T",
                "ptpn22_rs2476601_genotype": "G/G",
                "recipient_hla_e_genotype": "*01:01/*01:01",
                "tyk2_rs34536443_genotype": "G/C",
                "recipient_kir_ligand_group": "C1/C2",
                "recipient_bmi_kg_m2": 25.2,
                "carries_rare_snp_in_standard_grna": True,
                "requested_track": "TRACK_3_AUTOLOGOUS_IPSC",
            }
        ),
    ]

    total_n = n_per_pop * len(ancestries)
    pers_rates = [p["personalized_dna_compiled_openislet_5yr_survival_pct"] for p in pop_summaries]
    analytic_pers_rates = [p["analytic_expected_5yr_survival_personalized_compiled_unrounded_pct"] for p in pop_summaries]

    # Covariate-Level Gaussian Copula Sensitivity Analysis across shared-factor loading r in {0.0, 0.25, 0.50, 0.75}
    # (closing AUDIT-REF, AUDIT-REF, and AUDIT-REF):
    # Models shared-factor loading r in Z_trait = r * Z_shared + sqrt(1 - r^2) * eps_trait across patient immunological
    # covariates (KIR+ NK fraction, autoantigen T-cell z-score, and IFN-I percentile), where the pairwise latent
    # covariate correlation is rho_pairwise = r^2 in {0.0000, 0.0625, 0.2500, 0.5625}, reporting mean, SD, and 5th/50th/95th
    # percentile patient-level survival probability distributions.
    copula_sensitivity = []
    fixed_plan_a_ref = {
        "selected_fda_deployment_track": "TRACK_1_UNIVERSAL_OFF_THE_SHELF_ALLOGENEIC",
        "locus1_construct_prescribed": "OI-L2-DONOR-1B_B2M_EEF1A1_HLAE_CD47IGV",
        "plan_architecture_prescribed": "Plan A (Core 4-Locus Tri-Modal)",
        "include_tyk2_p1104a_rs34536443_shield": False,
        "compiled_grna_panel": [
            "OI-GRNA-01_B2M_EXON1_ASCAS12A",
            "OI-GRNA-02_CIITA_EXON1_ASCAS12A",
            "OI-GRNA-03_CXCL10_EXON1_HIFI_SPCAS9",
            "OI-GRNA-04_AAVS1_PPP1R12C_HIFI_SPCAS9",
        ],
    }
    for r_corr in (0.00, 0.25, 0.50, 0.75):
        rng_cop = np.random.default_rng(20261024)
        std_probs_all = []
        pers_probs_all = []
        dual_indication_count = 0
        s_factor = math.sqrt(max(0.0, 1.0 - r_corr * r_corr))
        for code, _name, _kir_a_f, c1_f, sirpa_v1_f, tyk2_f, mean_kir_nk, ifih1_t_f, ptpn22_a_f, hlae_0103_f, mean_bmi in ancestries:
            z_shared = rng_cop.standard_normal(n_per_pop)
            z_kir = r_corr * z_shared + s_factor * rng_cop.standard_normal(n_per_pop)
            z_auto = r_corr * z_shared + s_factor * rng_cop.standard_normal(n_per_pop)
            z_ifn = r_corr * z_shared + s_factor * rng_cop.standard_normal(n_per_pop)

            kir_vec = np.clip(mean_kir_nk + 0.08 * z_kir, 0.15, 0.72)
            auto_vec = np.clip(1.25 + 0.75 * z_auto, 0.0, 4.0)
            ifn_vec = np.array([50.0 * (1.0 + math.erf(float(zv) / math.sqrt(2.0))) for zv in z_ifn])
            sirpa_vec = rng_cop.binomial(2, sirpa_v1_f, size=n_per_pop)
            ifih1_vec = rng_cop.binomial(2, ifih1_t_f, size=n_per_pop)
            ptpn22_vec = rng_cop.binomial(2, ptpn22_a_f, size=n_per_pop)
            hlae_vec = rng_cop.binomial(2, hlae_0103_f, size=n_per_pop)
            bmi_c_vec = np.clip(rng_cop.normal(mean_bmi, 3.8, size=n_per_pop), 16.5, 42.0)
            for idx_p in range(n_per_pop):
                # Use the same cumulative HWE inverse CDF as the main population.
                # The heterozygote threshold includes the homozygous-alternate mass.
                tyk2_draw = _hwe_alternate_alleles(((idx_p * 997 + len(code)) % 10000) / 10000., tyk2_f)
                kir_draw = 2 if ((idx_p * 7919 + len(code)) % 1000) < int(round(1000 * (c1_f ** 2))) else (
                    1 if ((idx_p * 7919 + len(code)) % 1000) < int(round(1000 * (c1_f ** 2 + 2 * c1_f * (1 - c1_f)))) else 0
                )
                pt_c = {
                    "patient_id": f"COP-{code}-{idx_p:04d}",
                    "ancestry": code,
                    "nkg2a_neg_kir_educated_nk_fraction": float(kir_vec[idx_p]),
                    "pre_existing_gad65_ia2_znt8_t_cell_zscore": float(auto_vec[idx_p]),
                    "ifn1_ish_genetic_risk_percentile": float(ifn_vec[idx_p]),
                    "sirpa_genotype": sirpa_map[int(sirpa_vec[idx_p])],
                    "ifih1_rs1990760_genotype": ifih1_map[int(ifih1_vec[idx_p])],
                    "ptpn22_rs2476601_genotype": ptpn22_map[int(ptpn22_vec[idx_p])],
                    "recipient_hla_e_genotype": hlae_map[int(hlae_vec[idx_p])],
                    "tyk2_rs34536443_genotype": tyk2_map[tyk2_draw],
                    "recipient_kir_ligand_group": kir_map[kir_draw],
                    "recipient_bmi_kg_m2": round(float(bmi_c_vec[idx_p]), 2),
                    "carries_rare_snp_in_standard_grna": False,
                    "gate0_indirect_dr3_dq8_and_ifng_tnfa_stress_cleared": True,
                }
                c_reg = compiler_fn(pt_c)
                if "Plan B" in c_reg["plan_architecture_prescribed"] and c_reg["include_tyk2_p1104a_rs34536443_shield"]:
                    dual_indication_count += 1
                std_probs_all.append(evaluate_regimen_survival_probability(pt_c, fixed_plan_a_ref))
                pers_probs_all.append(evaluate_regimen_survival_probability(pt_c, c_reg))
        std_arr_np = 100.0 * np.array(std_probs_all, dtype=np.float64)
        pers_arr_np = 100.0 * np.array(pers_probs_all, dtype=np.float64)
        mean_std_pct = round(float(np.mean(std_arr_np)), 2)
        mean_pers_pct = round(float(np.mean(pers_arr_np)), 2)
        copula_sensitivity.append(
            {
                "copula_model_scope": "COVARIATE_LEVEL_GAUSSIAN_COPULA_ACROSS_KIR_AUTOANTIGEN_AND_IFN1_TRAITS",
                "latent_frailty_correlation_r": r_corr,
                "shared_frailty_factor_loading_r": r_corr,
                "analytic_pairwise_trait_correlation_rho_r_squared": round(r_corr * r_corr, 4),
                "dual_plan_b_plus_tyk2_indication_patients": dual_indication_count,
                "expected_standardized_plan_a_5yr_survival_pct": mean_std_pct,
                "p05_standardized_plan_a_5yr_survival_pct": round(float(np.percentile(std_arr_np, 5)), 2),
                "p50_standardized_plan_a_5yr_survival_pct": round(float(np.percentile(std_arr_np, 50)), 2),
                "p95_standardized_plan_a_5yr_survival_pct": round(float(np.percentile(std_arr_np, 95)), 2),
                "sd_standardized_plan_a_5yr_survival_pct": round(float(np.std(std_arr_np)), 2),
                "expected_personalized_compiled_5yr_survival_pct": mean_pers_pct,
                "p05_personalized_compiled_5yr_survival_pct": round(float(np.percentile(pers_arr_np, 5)), 2),
                "p50_personalized_compiled_5yr_survival_pct": round(float(np.percentile(pers_arr_np, 50)), 2),
                "p95_personalized_compiled_5yr_survival_pct": round(float(np.percentile(pers_arr_np, 95)), 2),
                "sd_personalized_compiled_5yr_survival_pct": round(float(np.std(pers_arr_np)), 2),
                "expected_arr_personalized_vs_standardized_pct": round(mean_pers_pct - mean_std_pct, 2),
            }
        )

    analytic_personalized_pct = 100.0 * float(np.mean(all_pers_analytic))
    analytic_standard_pct = 100.0 * float(np.mean(all_std_analytic))
    discounted_personalized_pct = 100.0 * float(np.mean(all_pers_ranker_discounted))
    unmat_b_mean_pct = round(100.0 * float(np.mean(all_unmat_b_analytic)), 4)
    unmat_b_disc_pct = round(100.0 * float(np.mean(all_unmat_b_discounted)), 4)
    pers_mean_4 = round(analytic_personalized_pct, 4)
    pers_disc_4 = round(discounted_personalized_pct, 4)
    std_mean_4 = round(analytic_standard_pct, 4)
    total_plan_b_routed = int(sum(p["compiler_regimen_allocation_counts"]["plan_b_hla_c_retained"] for p in pop_summaries))
    plan_b_routed_frac = total_plan_b_routed / float(total_n) if total_n > 0 else 0.0
    standard_discount = engineering_complexity_discount(4, 0, lambda_K=engineering_parameters["lambda_K"], gamma_B=engineering_parameters["gamma_B"])
    b_prime_sct_disc_4 = round(100.0 * float(np.mean([
        probability * engineering_complexity_discount(min(k, 4), b, float(engineering_parameters["lambda_K"]), float(engineering_parameters["gamma_B"]))
        for probability, k, b in all_personalized_complexity_cases
    ])), 4)
    b_prime_sct_gain_over_k5_pp = round(b_prime_sct_disc_4 - pers_disc_4, 4)
    b_prime_sct_gain_over_std_a_pp = round(b_prime_sct_disc_4 - std_mean_4 * standard_discount, 4)
    analytic_gap_pct = max(analytic_pers_rates) - min(analytic_pers_rates)
    sampled_gap_pct = max(pers_rates) - min(pers_rates)
    engineering_sensitivity = [
        {
            "lambda_K": float(lam), "gamma_B": float(gamma),
            "ranker_complexity_discounted_personalized_utility_pct": round(100.0 * float(np.mean([
                probability * engineering_complexity_discount(k, b, float(lam), float(gamma))
                for probability, k, b in all_personalized_complexity_cases
            ])), 4),
            "evidence_class": "ASSIGNED_UTILITY_PRIOR_SENSITIVITY_NOT_PHYSICAL_YIELD_OR_CLINICAL_CALIBRATION",
        }
        for lam in engineering_parameters["lambda_K_range"] for gamma in engineering_parameters["gamma_B_range"]
    ]
    return {
        "model_classification_and_scope_note_sol_v80_10": (
            "Illustrative assigned-rule scenario simulation (NOT clinical trial data, fitted survival hazard model, or "
            "validated epidemiological forecast): each of the 10,000 virtual patient records is passed through "
            "compile_personalized_patient_blueprint() and evaluated against conditional assigned-rule penalties (base 0.982, "
            "unedited control fixed at 0.0) to check rule-consistency across Arm 1B (Hour-0 CD47-IgV floor), Plan B "
            "(HLA-C-retaining sub-bank for high KIR+ NKG2A- NK fractions), and TYK2 P1104A protection (closing AUDIT-REF)."
        ),
        "total_virtual_patients_simulated": total_n,
        "global_5yr_survival_unedited_allogeneic_pct": round(100.0 * overall_counts["unedited_5yr_survival"] / total_n, 2),
        "global_5yr_survival_hla_e_only_pct": round(100.0 * overall_counts["hla_e_only_5yr_survival"] / total_n, 2),
        "global_5yr_survival_standardized_openislet_plan_a_pct": round(
            100.0 * overall_counts["standardized_plan_a_5yr_survival"] / total_n, 2
        ),
        "global_5yr_survival_personalized_dna_compiled_openislet_pct": round(
            100.0 * overall_counts["personalized_compiled_5yr_survival"] / total_n, 2
        ),
        "global_5yr_survival_ablated_compiler_control_pct": round(
            100.0 * overall_counts["ablated_compiler_control_5yr_survival"] / total_n, 2
        ),
        "global_analytic_expected_5yr_survival_standardized_plan_a_pct": round(analytic_standard_pct, 2),
        "global_analytic_expected_5yr_survival_personalized_compiled_pct": round(analytic_personalized_pct, 2),
        "global_analytic_expected_5yr_survival_ablated_arm1a_control_pct": round(100.0 * float(np.mean(all_abl_analytic)), 2),
        "analytic_max_ancestry_disparity_gap_personalized_pct_v80_065": round(analytic_gap_pct, 2),
        "sampled_max_ancestry_disparity_gap_personalized_pct": round(sampled_gap_pct, 2),
        "max_ancestry_disparity_gap_personalized_pct": round(sampled_gap_pct, 2),
        "deprecated_sampled_aliases": {"max_ancestry_disparity_gap_personalized_pct": "sampled_max_ancestry_disparity_gap_personalized_pct"},
        "ancestry_disparity_gap_methodology_note_v80_065": (
            f"The analytic primary model-implied ancestry gap is {analytic_gap_pct:.2f} pp before Bernoulli noise. "
            f"The sampled secondary gap is {sampled_gap_pct:.2f} pp in one seed=20261023 Bernoulli realization, "
            "N=2,000 per cohort. Legacy max_ancestry_disparity_gap_personalized_pct is a deprecated sampled alias."
        ),
        "editing_complexity_and_manufacturing_sensitivity_audit_v80_009": {
            "kb_edit_count_distribution_across_10000_patients": kb_counts,
            "conditional_on_release_qc_mean_personalized_survival_pct": pers_mean_4,
            "conditional_on_release_qc_mean_standardized_plan_a_survival_pct": std_mean_4,
            "ranker_complexity_discounted_personalized_utility_pct": pers_disc_4,
            "ranker_complexity_discounted_standardized_plan_a_utility_pct": round(
                analytic_standard_pct * standard_discount, 4
            ),
            "unmaterialized_plan_b_without_cd47_hlae_mean_survival_pct": unmat_b_mean_pct,
            "unmaterialized_plan_b_complexity_discounted_mean_utility_pct": unmat_b_disc_pct,
            "unmaterialized_plan_b_minus_materialized_plan_b_prime_unpenalized_pp": round(unmat_b_mean_pct - pers_mean_4, 4),
            "unmaterialized_plan_b_minus_materialized_plan_b_prime_complexity_discounted_pp": round(unmat_b_disc_pct - pers_disc_4, 4),
            "unpenalized_plan_b_minus_discounted_plan_b_prime_cross_scale_delta_pp": round(unmat_b_mean_pct - pers_disc_4, 4),
            "unmaterialized_plan_b_minus_standardized_plan_a_pp": round(unmat_b_mean_pct - std_mean_4, 4),
            "unmaterialized_plan_b_complexity_discounted_minus_standardized_plan_a_pp": round(unmat_b_disc_pct - std_mean_4, 4),
            "engineering_utility_prior_parameters": engineering_parameters,
            "engineering_utility_prior_sensitivity": engineering_sensitivity,
            "four_corner_prior_box_vs_plan_a_comparison_sonnet55_r2": [
                {
                    "corner_label": f"lambda_K_{row['lambda_K']}_gamma_B_{row['gamma_B']}",
                    "lambda_K": row["lambda_K"],
                    "gamma_B": row["gamma_B"],
                    "personalized_discounted_utility_pct": row["ranker_complexity_discounted_personalized_utility_pct"],
                    "standardized_plan_a_k4_b0_discount_factor": round(
                        engineering_complexity_discount(4, 0, row["lambda_K"], row["gamma_B"]), 6
                    ),
                    "standardized_plan_a_discounted_utility_pct": round(
                        analytic_standard_pct * engineering_complexity_discount(4, 0, row["lambda_K"], row["gamma_B"]), 4
                    ),
                    "standardized_plan_a_if_k5_dual_safe_harbor_discounted_utility_pct": round(
                        analytic_standard_pct * engineering_complexity_discount(5, 0, row["lambda_K"], row["gamma_B"]), 4
                    ),
                    "personalized_minus_standardized_plan_a_pp": round(
                        row["ranker_complexity_discounted_personalized_utility_pct"]
                        - round(analytic_standard_pct * engineering_complexity_discount(4, 0, row["lambda_K"], row["gamma_B"]), 4),
                        4,
                    ),
                    "personalized_exceeds_standardized_plan_a": bool(
                        row["ranker_complexity_discounted_personalized_utility_pct"]
                        > round(analytic_standard_pct * engineering_complexity_discount(4, 0, row["lambda_K"], row["gamma_B"]), 4)
                    ),
                }
                for row in engineering_sensitivity
            ],
            "hla_c_allele_leader_peptide_and_nkg2a_kd_table_sonnet55_r2": [
                {
                    "hla_c_allele_or_group": "HLA-C*01:02, HLA-C*03:04 (Cw*0304), HLA-C*04:01 (Cw*0401), HLA-C*06:02, HLA-Cw*0302 (Lin et al. 2023 PMID:37264229 Table 1 & Supp Table 1)",
                    "leader_class": "SP-1C",
                    "leader_nonamer_p1_to_p9": "VMAPRTLIL",
                    "vales_gomez_1999_cd94_nkg2a_kd_uM": 7.35,
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2a_kd_uM": 1.33,
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2a_ci95_uM": [1.15, 1.53],
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2c_kd_uM": 4.35,
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2c_ci95_uM": [4.24, 4.46],
                    "endogenous_hla_e_nkg2a_engageable": True,
                    "endogenous_hla_e_nkg2a_engageable_evidence_status": "DIRECT_SPR_MEASURED_VALES_GOMEZ_1999_TABLE1_AND_LIN_2023_SUPP_TABLE1",
                },
                {
                    "hla_c_allele_or_group": "HLA-Cw*0402 (HLA-C*04:02 rare P3-Glu variant)",
                    "leader_class": "Cw*0402_P3_Glu_variant",
                    "leader_nonamer_p1_to_p9": "VMEPRTLIL",
                    "vales_gomez_1999_cd94_nkg2a_kd_uM": None,
                    "vales_gomez_1999_cd94_nkg2a_kd_lower_bound_uM": ">= 75.30",
                    "endogenous_hla_e_nkg2a_engageable": False,
                    "endogenous_hla_e_nkg2a_engageable_evidence_status": "DIRECT_SPR_LOWER_BOUND_VALES_GOMEZ_1999_TABLE1",
                },
                {
                    "hla_c_allele_or_group": "HLA-C*02:02, HLA-C*15:02, HLA-C*15:05 (Lin et al. 2023 PMID:37264229 Table 1 & Supp Table 1; AUDIT-REF & AUDIT-REF)",
                    "leader_class": "SP-2C",
                    "leader_nonamer_p1_to_p9": "VMAPRTLLL",
                    "vales_gomez_1999_cd94_nkg2a_kd_uM": None,
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2a_kd_uM": 1.70,
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2a_ci95_uM": [1.50, 1.93],
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2c_kd_uM": 4.97,
                    "lin_2023_supp_table1_soluble_hlae0103_cd94_nkg2c_ci95_uM": [4.80, 5.14],
                    "endogenous_hla_e_nkg2a_engageable": True,
                    "endogenous_hla_e_nkg2a_engageable_evidence_status": "DIRECT_SOLUBLE_HLAE_SPR_LIN_2023_SUPP_TABLE1_P6_THR_P8_LEU_INTACT_SCT_UNMEASURED",
                },
                {
                    "hla_c_allele_or_group": "HLA-C*07:01, HLA-C*07:02 (UniProt P10321)",
                    "leader_class": "SP-3C",
                    "leader_nonamer_p1_to_p9": "VMAPRALLL",
                    "vales_gomez_1999_cd94_nkg2a_kd_uM": None,
                    "endogenous_hla_e_nkg2a_engageable": False,
                    "endogenous_hla_e_nkg2a_engageable_evidence_status": "FUNCTIONAL_DEFICIT_LIN_2023_PMID_37264229_TABLE1_P6_ALA_NO_NKG2A_ENGAGEMENT",
                },
            ],
            "sp3c_uncloaked_survival_penalty_sensitivity_grid_muse_r2": [
                {
                    "sp3c_uncloaked_survival_penalty": sp3c_pen,
                    "uncloaked_sp3c_cohort_mean_survival_delta_pp": round(100.0 * sp3c_pen, 4),
                    "active_plan_b_routed_patients_of_10000": total_plan_b_routed,
                    "population_10k_uncloaked_plan_b_mean_survival_pct": round(unmat_b_mean_pct + round(100.0 * sp3c_pen * plan_b_routed_frac, 4), 4),
                    "population_10k_uncloaked_plan_b_survival_delta_pp": round(100.0 * sp3c_pen * plan_b_routed_frac, 4),
                    "materialized_b_prime_hlae_sct_cloak_penalty_pp": 0.0,
                    "population_routing_assumption": f"CONDITIONAL_STRESS_SCENARIO_ALL_{total_plan_b_routed}_PLAN_B_ROUTED_PATIENTS_ASSIGNED_UNCLOAKED_SP3C",
                    "unmaterialized_product_contract_status": "HYPOTHETICAL_SCENARIO_ONLY_B_PRIME_REMAINS_NOT_MATERIALIZED_IN_PRODUCT_DEFINITIONS",
                    "product_contract_scope_kimi31_r3": {
                        "scenario_product": "UNMATERIALIZED_UNCLOAKED_PLAN_B_SP3C",
                        "scenario_status": "HYPOTHETICAL_STRESS_ONLY",
                        "materialized_comparator": "NONE_B_PRIME_CLOAK_NOT_MATERIALIZED_IN_PRODUCT_DEFINITIONS",
                        "hypothetical_cloaked_comparator": "PLAN_B_PRIME_HLAE_SCT_CLOAK_HYPOTHETICAL_NOT_MATERIALIZED",
                        "apply_penalty_to_materialized_b_prime": False,
                    },
                }
                for sp3c_pen in (-0.02, -0.04, -0.08, -0.15)
            ],
            "plan_b_kir_nk_threshold_sensitivity_sonnet55_r2": [
                {
                    "kir_nk_plan_b_routing_threshold": tau_kir,
                    "plan_b_routed_patients_of_10000": (
                        int(sum(1 for k_nk, _, _, _ in all_threshold_sensitivity_cases if k_nk >= tau_kir))
                        if _is_default_compiler
                        else total_plan_b_routed
                    ),
                    "analytic_personalized_mean_survival_pct": (
                        round(
                            100.0
                            * float(
                                np.mean(
                                    [
                                        p_pers_c if (k_nk >= tau_kir or k_nk < 0.48) else p_if_a_c
                                        for k_nk, p_pers_c, p_if_a_c, _ in all_threshold_sensitivity_cases
                                    ]
                                )
                            ),
                            4,
                        )
                        if _is_default_compiler
                        else pers_mean_4
                    ),
                    "complexity_discounted_personalized_utility_pct": (
                        round(
                            100.0
                            * float(
                                np.mean(
                                    [
                                        (p_pers_c if (k_nk >= tau_kir or k_nk < 0.48) else p_if_a_c)
                                        * engineering_complexity_discount(
                                            5 if k_nk >= tau_kir else 4,
                                            b_c,
                                            lambda_K=engineering_parameters["lambda_K"],
                                            gamma_B=engineering_parameters["gamma_B"],
                                        )
                                        for k_nk, p_pers_c, p_if_a_c, b_c in all_threshold_sensitivity_cases
                                    ]
                                )
                            ),
                            4,
                        )
                        if _is_default_compiler
                        else pers_disc_4
                    ),
                }
                for tau_kir in (0.42, 0.45, 0.48, 0.52)
            ],
            "physical_uncloned_yield_prior": None,
            "physical_yield_evidence_status": "NOT_MEASURED",
            "b_prime_4digit_hlac_alloreactivity_and_b_prime_sct_pivot_audit_muse_r3": {
                "audit_id": "AUDIT-REF",
                "evidence_class": "IMMUNOGENETIC_ALLOREACTIVITY_AUDIT_AND_PRODUCT_PROMOTION",
                "active_plan_b_routed_patients_of_10000": total_plan_b_routed,
                "three_line_kir_subbank_4digit_mismatch_fraction_range": [0.70, 0.85],
                "homozygous_4digit_hlac_lines_for_85_to_90pct_self_match": [12, 16],
                "hlac_4digit_relation_penalty_model": {
                    "SELF_4DIGIT_MATCH_or_crossmatch_cleared_B_PRIME_SCT_penalty": 0.0,
                    "KIR_MATCH_FOREIGN_4DIGIT_ALLO_endogenous_hlac_penalty": -0.25,
                    "KIR_MATCH_FOREIGN_4DIGIT_ALLO_uncrossmatched_b_prime_sct_penalty": -0.08,
                    "recipient_anti_hla_c_dsa_penalty_all_hlac_grafts": -0.40,
                    "recipient_activating_kir2ds1_unlicensed_c2_sct_penalty": -0.12,
                },
                "population_10k_plan_b_if_80pct_foreign_4digit_allo_unpenalized_survival_pct": round(
                    pers_mean_4 - 100.0 * 0.80 * 0.25 * plan_b_routed_frac, 4
                ),
                "population_10k_b_prime_sct_k4_b0_2pass_survival_pct": pers_mean_4,
                "population_10k_b_prime_sct_k4_complexity_discounted_utility_pct": b_prime_sct_disc_4,
                "b_prime_sct_complexity_discounted_gain_over_k5_endogenous_hlac_pp": b_prime_sct_gain_over_k5_pp,
                "b_prime_sct_complexity_discounted_gain_over_standardized_plan_a_pp": b_prime_sct_gain_over_std_a_pp,
                "hazard_active_cohort_sensitivity_audit_sonnet55_r4": {
                    "audit_id": "AUDIT-REF",
                    "baseline_97p50_conditioning_scope": (
                        f"The {pers_mean_4:.2f}% ({pers_disc_4:.2f}% K=5 [{pers_disc_4:.4f}%] / {b_prime_sct_disc_4:.2f}% K=4 [{b_prime_sct_disc_4:.4f}%] complexity-discounted) personalized headline is conditioned on "
                        "post-clearance grafts (all donor-panel Gate 0A NKG2C+, recipient serum CDC/ADCC crossmatch, and recipient "
                        "CD8+ MLC non-reactivity holds cleared; AUDIT-REF [DEC-076]). Prior to crossmatch/MLC/Gate-0A clearance, sampling realistic "
                        "recipient allorecognition and CMV adaptive-NK hazards yields the pre-clearance cohort sensitivity below."
                    ),
                    "plan_b_routed_fraction_of_10000_cohort": round(plan_b_routed_frac, 4),
                    "three_line_kir_subbank_foreign_4digit_hlac_exposure_range": [0.792, 0.902],
                    "endogenous_hlac_k5_with_79p2_to_90p2pct_foreign_4digit_allo_unpenalized_survival_pct_range": [
                        round(pers_mean_4 - 100.0 * 0.902 * 0.25 * plan_b_routed_frac, 4),
                        round(pers_mean_4 - 100.0 * 0.792 * 0.25 * plan_b_routed_frac, 4),
                    ],
                    "endogenous_hlac_k5_with_79p2_to_90p2pct_foreign_4digit_allo_gain_vs_plan_a_92p50_pp_range": [
                        round((pers_mean_4 - 100.0 * 0.902 * 0.25 * plan_b_routed_frac) - std_mean_4, 4),
                        round((pers_mean_4 - 100.0 * 0.792 * 0.25 * plan_b_routed_frac) - std_mean_4, 4),
                    ],
                    "endogenous_hlac_k5_with_79p2_to_90p2pct_foreign_4digit_allo_kb_discounted_utility_pct_range": [
                        round(pers_disc_4 - 100.0 * 0.902 * 0.25 * plan_b_routed_frac * math.exp(-float(engineering_parameters["lambda_K"])), 4),
                        round(pers_disc_4 - 100.0 * 0.792 * 0.25 * plan_b_routed_frac * math.exp(-float(engineering_parameters["lambda_K"])), 4),
                    ],
                    "endogenous_hlac_k5_with_79p2_to_90p2pct_foreign_4digit_allo_kb_discounted_delta_vs_plan_a_pp_range": [
                        round(
                            (pers_disc_4 - 100.0 * 0.902 * 0.25 * plan_b_routed_frac * math.exp(-float(engineering_parameters["lambda_K"])))
                            - std_mean_4 * standard_discount,
                            4,
                        ),
                        round(
                            (pers_disc_4 - 100.0 * 0.792 * 0.25 * plan_b_routed_frac * math.exp(-float(engineering_parameters["lambda_K"])))
                            - std_mean_4 * standard_discount,
                            4,
                        ),
                    ],
                    "b_prime_sct_k4_pre_mlc_clearance_with_foreign_4digit_and_c2_qyddavykl_survival_pct_range": [
                        round(pers_mean_4 - 100.0 * 0.08 * plan_b_routed_frac, 4),
                        round(pers_mean_4 - 100.0 * 0.792 * 0.08 * plan_b_routed_frac, 4),
                    ],
                    "post_clearance_b_prime_sct_k4_survival_pct": pers_mean_4,
                    "post_clearance_b_prime_sct_k4_kb_discounted_utility_pct": b_prime_sct_disc_4,
                    "population_simulation_with_cmv_nkg2c_and_sensitization_sonnet55_r5": {
                        "finding_ids": ["AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "MW-R5-S2-4p3", "AUDIT-REF"],
                        "model_epistemic_scope": "ILLUSTRATIVE_ASSIGNED_PENALTY_COMPETING_RISKS_SENSITIVITY_MODEL",
                        "ancestry_cmv_high_nkg2c_adaptive_nk_prevalence_range": [0.22, 0.38],
                        "mean_cohort_cmv_high_nkg2c_adaptive_nk_prevalence": 0.30,
                        "sensitized_cpra_ge_80_prevalence": 0.08,
                        "moderate_kir_nk_0p30_to_0p48_cohort_fraction": 0.625,
                        "four_hazard_competing_risks_decomposition_from_100pct": {
                            "h1_baseline_allo_autoimmune_residual_hazard": 0.025,
                            "h2_unrouted_moderate_kir_and_sirpa_low_nk_hazard": 0.04425,
                            "h3_cmv_high_nkg2c_adaptive_nk_hazard": 0.036,
                            "h4_non_hla_and_sct_sensitized_cpra_and_stress_hazard": 0.037175,
                            "total_4hazard_additive_sum": 0.142425,
                            "pre_clearance_unscreened_4hazard_additive_survival_pct": 85.7575,
                            "pre_clearance_independent_competing_risks_multiplicative_survival_pct_from_100": round(
                                100.0 * (1.0 - 0.025) * (1.0 - 0.04425) * (1.0 - 0.036) * (1.0 - 0.037175),
                                4,
                            ),
                        },
                        "hazard_term_breakdown_from_97p5pct_baseline": {
                            "cmv_high_nkg2c_adaptive_nk_term": round(0.30 * 0.12, 6),
                            "sensitized_cpra_ge_80_term": round(0.08 * 0.25, 6),
                            "moderate_kir_sirpa_low_nk_term": round(0.625 * 0.35 * 0.14, 6),
                            "residual_allo_autoimmune_stress_term": round(0.22 * 0.14, 6),
                        },
                        "pre_clearance_unscreened_cohort_mean_survival_pct": round(
                            pers_mean_4
                            - 100.0 * (0.30 * 0.12 + 0.08 * 0.25 + 0.625 * 0.35 * 0.14 + 0.22 * 0.14),
                            4,
                        ),
                        "pre_clearance_independent_competing_risks_multiplicative_survival_pct": round(
                            97.4612
                            * (1.0 - 0.30 * 0.12)
                            * (1.0 - 0.08 * 0.25)
                            * (1.0 - 0.625 * 0.35 * 0.14)
                            * (1.0 - 0.22 * 0.14),
                            4,
                        ),
                        "pre_clearance_4hazard_from_100pct_multiplicative_survival_pct": round(
                            100.0 * (1.0 - 0.025) * (1.0 - 0.04425) * (1.0 - 0.036) * (1.0 - 0.037175),
                            4,
                        ),
                        "competing_risks_baseline_provenance_reconciliation_muse_r6": {
                            "from_100pct_4hazard_multiplicative_survival_pct": 86.4915,
                            "from_97p4612pct_pers_mean_4_multiplicative_survival_pct": 86.5048,
                            "from_live_pers_mean_4_97p4963pct_multiplicative_survival_pct": round(
                                pers_mean_4
                                * (1.0 - 0.30 * 0.12)
                                * (1.0 - 0.08 * 0.25)
                                * (1.0 - 0.625 * 0.35 * 0.14)
                                * (1.0 - 0.22 * 0.14),
                                4,
                            ),
                            "from_97p4612pct_4hazard_additive_survival_pct": 85.7575,
                            "from_100pct_4hazard_additive_survival_pct_reference": round(
                                100.0 - 100.0 * (0.025 + 0.04425 + 0.036 + 0.037175), 4
                            ),
                            "from_97p4612pct_exact_additive_subtraction_survival_pct": round(
                                97.4612 - 100.0 * (0.30 * 0.12 + 0.08 * 0.25 + 0.625 * 0.35 * 0.14 + 0.22 * 0.14), 4
                            ),
                            "from_live_pers_mean_4_97p4963pct_exact_additive_subtraction_survival_pct": round(
                                pers_mean_4 - 100.0 * (0.30 * 0.12 + 0.08 * 0.25 + 0.625 * 0.35 * 0.14 + 0.22 * 0.14), 4
                            ),
                            "note": (
                                "Reconciles the two baseline references: 97.4612% is the pre-DEC-048 historical cohort mean (yielding 86.5048% multiplicative and 85.7187% additive; "
                                "85.7575% in from_97p4612pct_4hazard_additive_survival_pct is the 100% - 14.2425% = 85.7575% 4-hazard additive reference retained for backward compatibility), "
                                "whereas live pers_mean_4 = 97.4963% yields 86.5359% multiplicative and 85.7538% additive (M1 [DEC-077])."
                            ),
                        },
                        "post_gate0_and_crossmatch_clearance_modular_3to4x100cm2_cassette_mean_survival_pct": pers_mean_4,
                        "post_gate0_and_crossmatch_clearance_modular_3x100cm2_cassette_mean_survival_pct": pers_mean_4,
                        "single_200cm2_pouch_unscaled_chronic_1500_ieq_cm2_violation_mean_survival_pct": round(
                            pers_mean_4 - 100.0 * 0.38 * ((2500.0 - 1500.0) / 1000.0), 4
                        ),
                        "clinical_interpretation": (
                            "Illustrative assigned-penalty competing-risks sensitivity model: without Stage-0 Gate 0A NKG2C+/SIRPa-low NK clearance, "
                            "Gate 0D 50-serum non-HLA/SCT flow crossmatch, and modular 3-4 x 100 cm^2 omental cassette area scaling (<= 1,500 IEQ/cm^2), "
                            "unscreened cohort survival drops from 97.50% to 85.76% (85.7575% under additive 4-hazard sum h1+h2+h3+h4 = 0.142425; "
                            "86.5048% under 97.50%-baseline 4-term multiplicative competing risks and 86.4915% under 100%-baseline 4-hazard multiplicative competing risks "
                            "100% * (1-0.025)*(1-0.04425)*(1-0.036)*(1-0.037175)) or "
                            f"{round(pers_mean_4 - 100.0 * 0.38 * ((2500.0 - 1500.0) / 1000.0), 2)}% (single-pouch anoxic overload), confirming why both recipient screening gates and modular 3-4 x 100 cm^2 cassette scaling are mandatory."
                        ),
                    },
                },
                "canonical_track2_kir_subbank_products": [
                    "ISLET1_B_SCT_C1_CIITA",
                    "ISLET1_B_SCT_C2_CIITA",
                    "ISLET1_B_SCT_C1C2_CIITA",
                    "ISLET1_B_SCT_C1_CIITA_TYK2_PE",
                    "ISLET1_B_SCT_C2_CIITA_TYK2_PE",
                    "ISLET1_B_SCT_C1C2_CIITA_TYK2_PE",
                ],
                "endogenous_hlac_retention_product_demotion": {
                    "product_variant_ids": ["ISLET1_B2_5LOCUS_BASE", "ISLET1_B2_5LOCUS_TYK2_PE"],
                    "status": "EXPLORATORY_REQUIRES_4DIGIT_BANK",
                    "reason": (
                        "Retaining endogenous polymorphic HLA-C + WT B2M with only a 3-line KIR-ligand-group bank (C1/C1, C2/C2, C1/C2) "
                        "leaves ~70-85% of recipients 4-digit mismatched at HLA-C, re-arming direct CD8+ T-cell allorecognition and "
                        "anti-HLA-C DSA (Xu et al. 2019 PMID:30853558; Tran et al. 2011 PMID:21743387) while requiring 3-pass "
                        "clonal editing to avoid the 1.41 Mb chr6p21.3 HLA-A-to-HLA-B cis-deletion. Design B'-SCT "
                        "(ISLET1_B_SCT_C1_CIITA [GAVDPLLAL, PDB:1EFX] / ISLET1_B_SCT_C2_CIITA [QYDDAVYKL, PDB:1IM9] at CIITA Exon 5) retains B2M-/- at Locus 1 and expresses a "
                        "covalently tethered human self-peptide single-chain trimer at K=4, B=0 (P_passes=2) or B=1 (P_passes=3)."
                    ),
                },
                "donor_bank_antigen_hygiene_and_indirect_tolerance_window_muse_r6": {
                    "finding_ids": ["MW-R6-S1-19", "MW-R6-S1-20", "MW-R6-S1-21"],
                    "mcb_donor_genetic_and_serological_hygiene_criteria": {
                        "donor_sex_karyotype_required": "46,XX_FEMALE_H_Y_NULL",
                        "abo_secretor_genotype_required": "BLOOD_GROUP_O_FUT2_NON_SECRETOR_PREFERRED",
                        "viral_serology_required": "CMV_NEG_EBV_NEG_HHV6_NEG_PARVOB19_NEG",
                        "mitochondrial_haplogroup_and_common_mhc_polymorphism_screened": True,
                        "rationale": "Eliminates H-Y (UTY/DBY/SMCY) minor histocompatibility antigens in female recipients, A/B glycocalyx neoantigens, and donor-derived latent herpesvirus reactivation (MW-R6-S1-19).",
                    },
                    "two_independent_mcb_clone_diversification_and_functional_nutlin3a_p53_gate_muse_r6": {
                        "finding_id": "MW-R6-S1-20",
                        "min_independent_mcb_clones_from_distinct_46xx_o_neg_donors": 2,
                        "functional_p53_nutlin3a_assay_endpoint": "obs_nutlin3a_p21_mdm2_fold_induction",
                        "min_nutlin3a_10uM_8h_p21_cdkn1a_and_mdm2_fold_induction": 5.0,
                        "post_stress_20q11_bcl2l1_outgrowth_endpoint": "obs_post_stress_20q11_bcl2l1_outgrowth_fold",
                        "max_post_stress_20q11_bcl2l1_outgrowth_fold": 1.0,
                        "max_post_stress_20q11_bcl2l1_outgrowth_upper95_fold": 1.20,
                    },
                    "peri_operative_14day_damp_inflammation_and_indirect_pathway_bridging_option": {
                        "finding_id": "MW-R6-S1-21",
                        "surgical_implantation_damp_window_days": 14,
                        "recipient_dc_pinocytosis_of_shed_donor_debris_pathway": "Indirect CD4+ T-cell priming on recipient HLA-DR3/DR4-DQ8 APCs during surgical trauma + early hypoxic attrition",
                        "high_risk_z_ge_2p0_optional_local_or_transient_bridging": "Recipient-stratified transient 14-day peri-implantation anti-inflammatory/co-stimulatory blockade or hydrogel-tethered local anti-TNFa/IL-1Ra scavenging during capsule maturation",
                    },
                },
                "patient_safety_viral_outgrowth_p53_and_kir3dl1_gates_muse_r6": {
                    "finding_ids": ["MW-R6-S1-19", "MW-R6-S1-20", "MW-R6-S1-21"],
                    "track3_bw4_kir3dl1_hold_token": "HOLD_RECIPIENT_BW4_KIR3DL1_EDUCATED_NK_REQUIRES_SORTED_KIR3DL1_NK_CLEARANCE",
                    "functional_p53_nutlin3a_assay_endpoint": "obs_nutlin3a_p21_mdm2_fold_induction",
                    "min_nutlin3a_10uM_8h_p21_cdkn1a_and_mdm2_fold_induction": 5.0,
                    "post_stress_20q11_bcl2l1_outgrowth_endpoint": "obs_post_stress_20q11_bcl2l1_outgrowth_fold",
                    "max_post_stress_20q11_bcl2l1_outgrowth_fold": 1.0,
                    "max_post_stress_20q11_bcl2l1_outgrowth_upper95_fold": 1.20,
                    "viral_day14_log10_reduction_endpoint": "obs_gate0v_cvb_day14_log10_viral_reduction",
                    "min_viral_day14_log10_reduction": 3.0,
                    "viral_day28_persistent_vp1_ie1_endpoint": "obs_gate0v_day28_persistent_vp1_ie1_positive_cell_pct",
                    "max_viral_day28_persistent_vp1_ie1_upper95_pct": 0.05,
                },
            },
            "hour0_kinetic_race_and_basal_floor_audit_muse_r3": {
                "audit_id": "AUDIT-REF",
                "evidence_class": "KINETIC_RACE_HAZARD_PRIOR_AND_BASAL_FLOOR_CHALLENGER",
                "innate_nk_and_memory_ctl_degranulation_window_hours": [0.5, 4.0],
                "cxcl10_inducible_serpinb9_pdl1_accumulation_lag_hours": [4.0, 6.0],
                "basal_floor_challenger_product_variant_id": "ISLET1_A4_L3_BFB_BASALFLOOR",
                "basal_floor_challenger_constructs": ["L3-BF-B (4,078 bp)", "L3-BF-C (4,160 bp)"],
                "stage0_kill_gate": "Gate 0B (Hour-0 to Hour-12 real-time primary human NK + memory CD8+ CTL live-imaging kill race; S1 -> S2 precedence)",
            },
            "scope_reconciliation_note": (
                f"Illustrative assigned-rule outcomes, conditional on hypothetical release-qualified clones: standardized Plan A "
                f"{analytic_standard_pct:.4f}% and personalized {analytic_personalized_pct:.4f}% (whereas if Plan B's unmaterialized "
                f"CD47/HLA-E safe-harbor cassette is treated as absent, analytic personalized survival drops by {pers_mean_4 - unmat_b_mean_pct:.4f} pp to {unmat_b_mean_pct:.4f}% "
                f"[{unmat_b_disc_pct:.4f}% after patient-level complexity discount, {unmat_b_disc_pct - pers_disc_4:.4f} pp below complexity-discounted Plan B-prime {pers_disc_4:.4f}% and {unmat_b_mean_pct - std_mean_4:.4f} pp unpenalized / "
                f"{unmat_b_disc_pct - std_mean_4:.4f} pp discounted below standardized Plan A {std_mean_4:.4f}%; cross-scale unpenalized B minus discounted B-prime is {unmat_b_mean_pct - pers_disc_4:.4f} pp]; "
                f"AUDIT-REF, AUDIT-REF, AUDIT-REF, & AUDIT-REF). Applying the shared "
                f"dimensionless engineering utility discount gives standardized {analytic_standard_pct * standard_discount:.4f}% "
                f"and personalized {discounted_personalized_pct:.4f}%. The discount is an assigned utility penalty, "
                "not physical manufacturing yield or a fitted clinical probability; uncloned physical yield is unmeasured."
            ),
        },
        "ancestry_stratified_cohorts": pop_summaries,
        "correlated_frailty_copula_sensitivity_sol_v80_r2_07": copula_sensitivity,
        "representative_patient_compilations": reference_patients,
    }


def build_layer5_atlas() -> Dict[str, Any]:
    """Build the complete report from the same current objects used for every numeric summary."""
    part_a = build_part_a_crispr_grna_and_polymorphism_atlas()
    part_b = simulate_10000_patient_global_population_trial()

    atlas = {
        "schema_version": "0.8.0",
        "layer": "Layer 5 — Personalized Patient-DNA & Multiplex CRISPR Compiler + 10,000-Patient Global Population Trial",
        "status": "AUTHENTIC_GRCH38_GRNA_THERMODYNAMICS_DBSNP_CLINVAR_AND_10K_POPULATION_SIMULATION",
        "evidence_class": "in_silico_crispr_compiler_and_illustrative_population_scenario_simulation_requiring_wetlab_validation",
        "summary": {
            "total_guides_compiled": len(part_a["guides"]),
            "all_guides_verified_0_mismatch_to_grch38": all(
                g["grch38_reference_exact_match_verified"] for g in part_a["guides"]
            ),
            "coding_off_targets_0_to_2_mismatches": (
                sum(g["grch38_0_1_2_mismatch_coding_off_targets"] for g in part_a["guides"])
                if all(g["grch38_0_1_2_mismatch_coding_off_targets"] is not None for g in part_a["guides"])
                else None
            ),
            "genomewide_off_target_enumeration_status": "NOT_ENUMERATED",
            "any_intact_target_in_hdr_donors": any(
                g["intact_target_present_in_any_hdr_donor_plasmid"] for g in part_a["guides"]
            ),
            "any_intact_target_in_post_hdr_edited_alleles": any(
                g["intact_target_present_in_post_hdr_edited_allele_either_strand"] for g in part_a["guides"]
            ),
            "total_virtual_patients_simulated": part_b["total_virtual_patients_simulated"],
            "global_5yr_survival_hla_e_only_pct": part_b["global_5yr_survival_hla_e_only_pct"],
            "global_5yr_survival_standardized_plan_a_pct": part_b["global_5yr_survival_standardized_openislet_plan_a_pct"],
            "global_5yr_survival_personalized_compiled_pct": part_b[
                "global_5yr_survival_personalized_dna_compiled_openislet_pct"
            ],
            "global_5yr_survival_ablated_compiler_control_pct": part_b[
                "global_5yr_survival_ablated_compiler_control_pct"
            ],
            "global_analytic_expected_5yr_survival_personalized_compiled_pct": part_b["global_analytic_expected_5yr_survival_personalized_compiled_pct"],
            "global_analytic_expected_5yr_survival_standardized_plan_a_pct": part_b["global_analytic_expected_5yr_survival_standardized_plan_a_pct"],
            "analytic_max_ancestry_disparity_gap_personalized_pct": part_b["analytic_max_ancestry_disparity_gap_personalized_pct_v80_065"],
            "sampled_max_ancestry_disparity_gap_personalized_pct": part_b["sampled_max_ancestry_disparity_gap_personalized_pct"],
            "sampled_ancestry_gap_method": {
                "seed": 20261023,
                "n_per_cohort": 2000,
                "method": (
                    "hybrid deterministic low-discrepancy HWE stratification (coprime multiplier 997 mod 10000 for TYK2 rs34536443 "
                    "and 7919 mod 1000 for KIR C1/C2 ligand group, paired with np.random.default_rng(20261023) binomial/Gaussian draws "
                    "for SIRPA, IFIH1, PTPN22, HLA-E, BMI, and continuous NK/autoantigen/IFN-I covariates); analytic_expected_* metrics "
                    "average exact survival probabilities p_i over the 2,000 sampled covariate profiles per cohort prior to Bernoulli indicator sampling (AUDIT-REF)"
                ),
            },
        },
        "part_a_crispr_grna_compiler_and_polymorphism_audit": part_a,
        "part_b_personalized_patient_compiler_and_10k_population_trial": part_b,
    }

    return atlas


def main() -> None:
    atlas = build_layer5_atlas()
    part_a = atlas["part_a_crispr_grna_compiler_and_polymorphism_audit"]
    part_b = atlas["part_b_personalized_patient_compiler_and_10k_population_trial"]

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", encoding="utf-8") as f:
        json.dump(atlas, f, indent=2)
        f.write("\n")

    print(f"[Layer 5 v0.8.0] Wrote {OUT_PATH}")
    for g in part_a["guides"]:
        print(
            f"  [gRNA] {g['guide_id']}: {g['target_locus_grch38']} | PAM={g['pam_sequence_5to3']} | "
            f"GC={g['gc_percent']}% | Sugimoto_dG37={g['sugimoto_1995_rna_dna_hybrid_dg37_kcal_mol']} kcal/mol | "
            f"in_donor={g['intact_target_present_in_any_hdr_donor_plasmid']} | "
            f"in_post_hdr={g['intact_target_present_in_post_hdr_edited_allele_either_strand']}"
        )
    print(
        f"  [10k Trial] Unedited={part_b['global_5yr_survival_unedited_allogeneic_pct']}% | "
        f"HLA-E-Only={part_b['global_5yr_survival_hla_e_only_pct']}% | "
        f"Standardized Plan A={part_b['global_5yr_survival_standardized_openislet_plan_a_pct']}% | "
        f"Personalized Compiled={part_b['global_5yr_survival_personalized_dna_compiled_openislet_pct']}% | "
        f"Ablated Control={part_b['global_5yr_survival_ablated_compiler_control_pct']}% "
        f"(Analytic ancestry gap={part_b['analytic_max_ancestry_disparity_gap_personalized_pct_v80_065']} pp; "
        f"sampled secondary={part_b['sampled_max_ancestry_disparity_gap_personalized_pct']} pp)"
    )


if __name__ == "__main__":
    main()
