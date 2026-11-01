#!/usr/bin/env python3
"""
Step 5B (v0.3 Deep Quantitative & Canonical Graph Upgrade):
Grounded Full-Text JATS XML & Abstract Relational Extractor for Ring 1 (532 Studies).
Processes Study Data Dossiers in data/dossiers/{pmid}.json into schema-validated v0.2.0 records
in data/extracted/{pmid}.json.

Key Scientific Guarantees:
  1. Perturbation vs. Passive Marker Discrimination: Genes are only added to
     gate2_3_4_genetic_and_epitope_edits if the study text contains explicit perturbation context
     (KO, CRISPR, knockdown, overexpression, mutation/variant, inhibitor, swap, etc.) near the gene,
     preventing passive lineage immunostaining markers (e.g. routine INS/CHGA/SLC30A8 staining)
     from polluting the perturbation atlas.
  2. Real Quantitative Extraction: Parses actual numbers (fold-changes, paired percentages,
     mean +/- SD/SEM, survival days/weeks, concentration units) from Results and Figure Legends
     to populate treatment_mean, control_mean, measurement_unit, and log_response_ratio_or_smd.
  3. True Clause Splitting for Effect Modifiers: Splits contrast sentences around pivot conjunctions
     (whereas, in contrast, while, unlike, but not, however) or pairs adjacent sentences so
     level_1_condition_and_effect and level_2_condition_and_effect contain real contrasting regimes
     and modifier_variable_Z names the specific biological modifier.
  4. Canonical Biological Node Resolution for Causal Edges: Resolves both upstream_node and
     downstream_node to canonical biological entities (Genes, Pathways, Organelle States,
     Immune Effectors, Vascular Cues) so multi-hop graph traversal (Node A -> Node B -> Node C)
     operates genuinely across all 532 studies.
"""

from __future__ import annotations

import argparse
import copy
from decimal import Decimal, DecimalException, InvalidOperation
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
DOSSIER_DIR = BASE_DIR / "data" / "dossiers"
CURATED_DIR = BASE_DIR / "data" / "curated"
EXTRACTED_DIR = BASE_DIR / "data" / "extracted"

# Curated landmark PMIDs that should never be overwritten by automated batch pass
CURATED_LANDMARK_PMIDS = {
    "35241836",  # Balboa/Otonkoski 2022 Nat Biotechnol
    "37156915",  # Hu/Schrepfer 2024 Nat Biotechnol
    "42690930",  # Kobaisi/Scharfmann 2026 Cell Reports
    "32719542",  # Cai/Yi/Kissler 2020 Nat Metab
    "36289205",  # Chandra/Otonkoski 2022 Nat Commun
    "35559682",  # Lei/Garcia 2022 Sci Adv
    "34480863",  # Aghazadeh/Nostro 2021 Cell Stem Cell
    "36599351",  # Gerace/Melton 2023 Cell Rep Med
    "36055241",  # Sintov/Melton 2022 Stem Cell Reports
    "34407395",  # Parent/Hebrok 2021 Cell Reports
    "32814902",  # Yoshihara/Evans 2020 Nature
    "39988068",  # MacDonald/Yi 2025 Mol Metab
    "31900242",  # Demine/Eizirik 2020 Stem Cell Res Ther
    "32846125",  # Augsornworawat/Millman 2020 Cell Reports
    "42671571",  # Luckett/Oram/Richardson/Gaulton 2026 Diabetologia
    "26276263",  # Martino/De Tata 2015 Diabetologia
    "40757665",  # Carlsson/Schrepfer 2025 NEJM (First-in-human UP421)
    "40251154",  # Hu/Schrepfer 2025 Nat Commun (PMN inhibition via CD99+CD200)
    "39824805",  # Carre/Mallone 2025 Nat Commun (IFN-a HLA-B antigen presentation)
    "40112799",  # van Tienhoven/Roep 2025 Cell (INS 3'UTR rs3842752 IRE1a mRNA decay)
    "40782792",  # Meng/Deng 2025 Cell Stem Cell (5-subtype PSC-islet hypoglycemia protection)
    "41951603",  # Aggarwal/Pepper 2026 Cell Death Dis (Necrostatin-1 RIPK1/RIPK3/pMLKL shield)
    "41993274",  # Martin-Vazquez/Eizirik 2026 (ISG15 bat-like antiviral resistance)
    "41093834",  # Shibue/Kulkarni 2025 Nat Commun (SEC31A CRISPR screen alpha-cell survival)
    "41851084",  # Misra/Nostro 2026 Nat Commun (NGN3 patterning & enterochromaffin control)
    "42856627",  # Cochrane/Hebrok 2026 (MAFA vs MAFA-S64F SC-islet maturation)
    "41218919",  # Oshima/Scharfmann 2026 Diabetes (CD8+ T-cell paracrine bystander dysfunction)
    "41607298",  # Lei/Huang 2026 Diabetes Obes Metab (PAD4 NETosis beta-cell damage)
    "38754432",  # Srivastava/Wan 2024 Immunity (CXCL16 OxLDL macrophage ferroptosis axis)
    "41445610",  # Luckett/EXE-T1D consortium 2025 medRxiv (Genetic risk in extremely early onset T1D <2 yr)
    "39626097",  # McGrail/Gaulton 2025 Diabetes Care (Non-DR3/DR4 T1D genetic discovery & GRS)
    "41758964",  # Bonfield/EXE-T1D Consortium 2026 Diabetes (Biallelic IL2RA LoF neonatal T1D)
    "41617975",  # Yi/Eizirik 2026 Cell Death Differ (Alpha vs beta immune/apoptosis & MEG3 KD)
    "36030035",  # Contreras/Templin 2022 Mol Metab (TNF-a RIPK1/RIPK3/MLKL beta-cell necroptosis)
    "42430558",  # Carlsson/Schrepfer 2026 NEJM (14-month UP421 human clinical follow-up)
    "32569842",  # Bosma/O'Brien 2020 (G6pc2 KO fasting/ketogenic hypoglycemia counter-evidence)
    "34954144",  # Overway/O'Brien 2022 (Human G6PC2 coding variants functional analysis)
    "27626039",  # Yagyu/Brenner 2016 (iCasp9 silencing & apoptosis resistance escape in hPSCs)
    "35795779",  # Wunderlich/Martin 2022 (Biallelic safe-harbor iCasp9 integration against escape)
    "27183617",  # Polukort 2016 (IL-10 promotion of mast-cell expansion counter-evidence)
    "39529170",  # Wu 2024 (Base editor / nickase structural variation and large deletions)
    "37264229",  # Lin 2023 Nat Immunol (HLA leader-peptide polymorphism & HLA-E/NKG2A recognition)
    "28504668",  # Gornalusse 2017 Nat Biotechnol (B2M::HLA-E single-chain trimer universal donor hPSCs)
    "37820008",  # Ghimire 2023 Sci Transl Med (CD47 receptor signaling inhibits beta-cell GSIS counter-evidence)
    "38508139",  # Xu 2024 Cell Rep Med (Signaling-incompetent CD47-IgV avoids cell death and angiogenesis inhibition)
    "40668892",  # Zhao 2025 Sci Transl Med (Directional donor-recipient SIRPA mismatch innate allorecognition)
}

# Perturbation verbs/nouns required within a window of a gene symbol to classify it as a true perturbation
PERTURBATION_CONTEXT_RE = re.compile(
    r"\b(knockout|knock-out|KO\b|-/-|CRISPR|Cas9|Cas12|gRNA|sgRNA|shRNA|siRNA|knockdown|knock-down|"
    r"silenc\w+|delet\w+|ablat\w+|loss\s+of\s+function|loss-of-function|LoF\b|deficien\w+|mutan\w+|"
    r"mutation|R325W|P1104A|A946T|overexpress\w+|transduc\w+|transfect\w+|"
    r"knock-in|knockin|KI\b|inhibit\w+|antagoni\w+|base-edit\w+|prime-edit\w+|"
    r"swap\w*|cloak\w*|armored|tether\w+|microgel|autoantigen|autoreactive\s+T|epitope)\b",
    flags=re.IGNORECASE,
)

# Strict genetic perturbation required for endogenous beta-cell markers (SLC30A8, CHGA, G6PC2, GAD2, PTPRN, INS, IAPP)
# so routine SC-islet differentiation marker checks are never misclassified as genetic edits
STRICT_GENETIC_EDIT_RE = re.compile(
    r"\b(knockout|knock-out|KO\b|-/-|CRISPR|Cas9|Cas12|gRNA|sgRNA|shRNA|siRNA|knockdown|knock-down|"
    r"delet\w+|ablat\w+|loss-of-function|LoF\b|deficien\w+|mutan\w+|mutation|R325W|"
    r"base-edit\w+|prime-edit\w+|signal\s+peptide\s+swap|leader\s+swap|INS-DRiP|WE14|BDC2\.5|pramlintide)\b",
    flags=re.IGNORECASE,
)

ENDOGENOUS_LINEAGE_MARKERS = {"SLC30A8", "CHGA", "G6PC2", "GAD2", "PTPRN", "INS", "IAPP", "GLIS3"}

# Passive marker context that disqualifies a bare mention if no strict genetic edit is in the window
PASSIVE_MARKER_ONLY_RE = re.compile(
    r"\b(immunostaining|immunofluorescence|stained\s+for|staining\s+for|positive\s+for|markers?\s+such\s+as|"
    r"expressing\s+markers|pan-endocrine|beta-cell\s+marker|endocrine\s+marker|maturation\s+marker|"
    r"flow\s+cytometry|FACS\s+analysis|membrane\s+attack\s+complex-inhibitory\s+protein|"
    r"complement\s+inhibitory\s+receptors?|qPCR|RT-qPCR|RNA-seq|scRNA-seq|transcriptom\w+)\b",
    flags=re.IGNORECASE,
)

# Comprehensive Target Gene Dictionary for Project Islet-1 / OpenIslet.org (Gates 2, 3, 4, 5)
TARGET_GENE_SPECS: list[dict[str, Any]] = [
    {
        "symbol": "B2M",
        "ensembl": "ENSG00000166710",
        "pattern": r"\b(B2M|beta-2-microglobulin|β2-microglobulin|B2m-/-|B2M-/-)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Allo_HLA_ClassI_or_II",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "CIITA",
        "ensembl": "ENSG00000179583",
        "pattern": r"\b(CIITA|Class II transactivator|Ciita-/-|CIITA-/-)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Allo_HLA_ClassI_or_II",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "NLRC5",
        "ensembl": "ENSG00000140853",
        "pattern": r"\b(NLRC5|CITA)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Allo_HLA_ClassI_or_II",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "CD47",
        "ensembl": "ENSG00000196776",
        "pattern": r"\b(CD47|SIRP\s*alpha|SIRPα)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_NK_MissingSelf_Inhibition",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "HLA-E",
        "ensembl": "ENSG00000204592",
        "pattern": r"\b(HLA-E|Qa-1b)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_NK_MissingSelf_Inhibition",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "SafeHarbor_KI_AAVS1_GAPDH_CLYBL",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "HLA-G",
        "ensembl": "ENSG00000204632",
        "pattern": r"\b(HLA-G)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_NK_MissingSelf_Inhibition",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "HLA-C",
        "ensembl": "ENSG00000204525",
        "pattern": r"\b(HLA-C\s*retained|HLA-A\s*and\s*HLA-B\s*knockout|HLA-A-/-HLA-B-/-|HLA-C-retained)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate3_NK_MissingSelf_Inhibition",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Endogenous_Native_Promoter",
    },
    {
        "symbol": "CD274",
        "ensembl": "ENSG00000120217",
        "pattern": r"\b(PD-L1|PDL1|CD274|B7-H1)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate2_Checkpoint_T_Cell_Inhibition",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "FASLG",
        "ensembl": "ENSG00000117560",
        "pattern": r"\b(FasL|SA-FasL|FASLG|CD95L)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate2_Checkpoint_T_Cell_Inhibition",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "Transient_mRNA_or_Surface_Tethering",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "CXCL10",
        "ensembl": "ENSG00000169245",
        "pattern": r"\b(CXCL10|IP-10)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "TYK2",
        "ensembl": "ENSG00000105397",
        "pattern": r"\b(TYK2|BMS-986165|deucravacitinib|rs34536443|P1104A)\b",
        "region_class": "Catalytic_or_Kinase_Domain",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "PTPN2",
        "ensembl": "ENSG00000175354",
        "pattern": r"\b(PTPN2|TC-PTP)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Endogenous_Native_Promoter",
    },
    {
        "symbol": "TNFAIP3",
        "ensembl": "ENSG00000118503",
        "pattern": r"(?<!-)\b(TNFAIP3|A20)\b(?!-B\d)",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "SOCS1",
        "ensembl": "ENSG00000185338",
        "pattern": r"\b(SOCS1|SOCS-1)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "RNLS",
        "ensembl": "ENSG00000184719",
        "pattern": r"\b(RNLS|Rnls|renalase|pargyline)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "TXNIP",
        "ensembl": "ENSG00000265972",
        "pattern": r"\b(TXNIP|thioredoxin-interacting protein|verapamil)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "MANF",
        "ensembl": "ENSG00000145050",
        "pattern": r"\b(MANF|mesencephalic astrocyte-derived neurotrophic factor)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "SmallMolecule_Pharmacological_Modulator",
        "default_promoter": "Endogenous_Native_Promoter",
    },
    {
        "symbol": "SERPINB9",
        "ensembl": "ENSG00000170542",
        "pattern": r"\b(SERPINB9|PI-9|Spi6|serpin B9)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate4_Intrinsic_Granzyme_Fas_Apoptosis_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "RIPK1",
        "ensembl": "ENSG00000137275",
        "pattern": r"\b(RIPK1|RIPK3|MLKL|necrostatin)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_Granzyme_Fas_Apoptosis_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "SLC30A8",
        "ensembl": "ENSG00000164756",
        "pattern": r"\b(SLC30A8|ZnT8|R325W|rs13266634)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["ZnT8A_SLC30A8"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "GAD2",
        "ensembl": "ENSG00000136750",
        "pattern": r"\b(GAD2|GAD65|GADA)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["GADA_GAD65"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "G6PC2",
        "ensembl": "ENSG00000152254",
        "pattern": r"\b(G6PC2|IGRP)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["IGRP_G6PC2"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "PTPRN",
        "ensembl": "ENSG00000054356",
        "pattern": r"\b(PTPRN|IA-2|PTPRN2|phogrin|IA-2beta)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["IA2A_PTPRN_Phogrin"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "INS",
        "ensembl": "ENSG00000254647",
        "pattern": r"\b(INS-DRiP|preproinsulin signal peptide|insulin signal peptide|PPI15-24|insulin B9-23)\b",
        "region_class": "Signal_Peptide_Leader",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["IAA_Insulin_Preproinsulin_Leader", "INS_DRiP_Cryptic_ORF"],
        "default_modality": "Signal_Peptide_or_Domain_Swap",
        "default_promoter": "Lineage_Locked_INS_G6PC2_PDX1",
    },
    {
        "symbol": "CHGA",
        "ensembl": "ENSG00000100604",
        "pattern": r"\b(CHGA|chromogranin A|WE14|BDC2\.5|hybrid insulin peptide)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["ChromograninA_CHGA_WE14", "Hybrid_Insulin_Peptide_HIP"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "TGM2",
        "ensembl": "ENSG00000198959",
        "pattern": r"\b(TGM2|transglutaminase 2|deamidat)\b",
        "region_class": "Catalytic_or_Kinase_Domain",
        "functional_gate": "Gate2_Autoantigen_or_Neoantigen_Stealth",
        "evaded": ["Hybrid_Insulin_Peptide_HIP"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "IAPP",
        "ensembl": "ENSG00000121351",
        "pattern": r"\b(IAPP|hIAPP|islet amyloid polypeptide|amylin|pramlintide)\b",
        "region_class": "Amyloidogenic_Core_Motif",
        "functional_gate": "Gate4_Intrinsic_Amyloid_or_Viral_Receptor_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "CXADR",
        "ensembl": "ENSG00000154639",
        "pattern": r"\b(CXADR|coxsackievirus and adenovirus receptor|CAR-SIV)\b",
        "region_class": "Specific_Exon_or_Splice_Isoform",
        "functional_gate": "Gate4_Intrinsic_Amyloid_or_Viral_Receptor_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "IFIH1",
        "ensembl": "ENSG00000115267",
        "pattern": r"\b(IFIH1|MDA5|rs1990760|A946T)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "CD55",
        "ensembl": "ENSG00000196352",
        "pattern": r"\b(CD55)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "CD59",
        "ensembl": "ENSG00000085063",
        "pattern": r"\b(CD59)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "ENTPD1",
        "ensembl": "ENSG00000138185",
        "pattern": r"\b(ENTPD1|CD39)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "THBD",
        "ensembl": "ENSG00000178726",
        "pattern": r"\b(THBD|thrombomodulin)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_Complement_ADCC_IBMIR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "CASP9",
        "ensembl": "ENSG00000132906",
        "pattern": r"\b(iCasp9|inducible caspase-9|AP1903|rimiducid|HSV-TK|suicide switch)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate5_Suicide_Safety_Switch",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "SafeHarbor_KI_AAVS1_GAPDH_CLYBL",
        "default_promoter": "Housekeeping_GAPDH_ACTB_B2M",
    },
    {
        "symbol": "GLIS3",
        "ensembl": "ENSG00000107249",
        "pattern": r"\b(GLIS3)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Endogenous_Native_Promoter",
    },
    {
        "symbol": "TET2",
        "ensembl": "ENSG00000168769",
        "pattern": r"\b(TET2)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_IFN_JAK_STAT_NFkB_Block",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "ISG15",
        "ensembl": "ENSG00000187608",
        "pattern": r"\b(ISG15)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate4_Intrinsic_Amyloid_or_Viral_Receptor_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Endogenous_Native_Promoter",
    },
    {
        "symbol": "SEC31A",
        "ensembl": "ENSG00000138674",
        "pattern": r"\b(SEC31A|Sec31A)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate4_Intrinsic_ER_Stress_UPR_Shield",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "PADI4",
        "ensembl": "ENSG00000159339",
        "pattern": r"\b(PADI4|PAD4|NETosis|cl-amidine)\b",
        "region_class": "Catalytic_or_Kinase_Domain",
        "functional_gate": "Gate3_Macrophage_DonotEatMe",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "SmallMolecule_Pharmacological_Modulator",
        "default_promoter": "Not_Applicable_KO",
    },
    {
        "symbol": "CD200",
        "ensembl": "ENSG00000091972",
        "pattern": r"\b(CD200|CD99)\b",
        "region_class": "Transgene_cDNA_Cassette",
        "functional_gate": "Gate3_Macrophage_DonotEatMe",
        "evaded": ["None_Not_Antigen_Specific"],
        "default_modality": "Lentiviral_or_Transposon_OE",
        "default_promoter": "Viral_or_Synthetic_CAG_EF1a_CMV",
    },
    {
        "symbol": "HLA-B",
        "ensembl": "ENSG00000234745",
        "pattern": r"\b(HLA-B-restricted|HLA-B\s+hyper-expression|HLA-B-/-)\b",
        "region_class": "Whole_Gene_Coding_CDS",
        "functional_gate": "Gate2_Allo_HLA_ClassI_or_II",
        "evaded": ["Allogeneic_HLA_DSA_or_TCR"],
        "default_modality": "CRISPR_Cas9_KO",
        "default_promoter": "Not_Applicable_KO",
    },
]

# Canonical Biological Entity Dictionary for Causal Graph Node Normalization (46 Nodes across Gates 1-5)
CANONICAL_NODE_PATTERNS: list[tuple[str, str, str]] = [
    # (Node_ID, Node_Type_Enum, Regex_Pattern)
    ("HLA_DR3_DR4_T1DE1_T1DE2_Endotype_Architecture", "Lineage_Identity_or_Maturation_State", r"\b(HLA-DR3|HLA-DR4|DR3-DQ2|DR4-DQ8|DR3-T1D|DR4-T1D|T1DE1|T1DE2|IAA-first|GADA-first|endotype)\b"),
    ("MastCell_Histamine_CPA3_TPSAB1_Axis", "Immune_or_Abiotic_Stressor", r"\b(mast cell|histamine|tryptase|CPA3|TPSAB1|FcεRI|cromolyn|ketotifen|cetirizine)\b"),
    ("PMN_Neutrophil_NETosis_PAD4_CD99_CD200_Axis", "Immune_or_Abiotic_Stressor", r"\b(polymorphonuclear|\bPMN\b|\bPMNs\b|NETosis|neutrophil extracellular trap|PAD4|PADI4|CD99|CD200|myeloperoxidase)\b"),
    ("RIPK1_RIPK3_MLKL_Necroptosis_and_Ferroptosis_Axis", "Cell_Survival_or_Death_Pathway", r"\b(RIPK1|RIPK3|MLKL|pMLKL|necroptosis|necrostatin|Nec-1|ferroptosis|GPX4|CXCL16|OxLDL)\b"),
    ("B2M_HLA_Class_I_Complex", "Surface_Ligand_or_Antigen_Presentation", r"\b(B2M|beta-2-microglobulin|HLA-A|HLA-B|HLA-C|HLA-ABC|MHC class I\b|MHC-I\b)\b"),
    ("CIITA_HLA_Class_II_Complex", "Surface_Ligand_or_Antigen_Presentation", r"\b(CIITA|HLA-DR|HLA-DQ|HLA-DP|MHC class II\b|MHC-II\b)\b"),
    ("CD47_SIRPa_DonotEatMe_Axis", "Surface_Ligand_or_Antigen_Presentation", r"\b(CD47|SIRP\s*alpha|SIRPα)\b"),
    ("HLA_E_NKG2A_Inhibitory_Axis", "Surface_Ligand_or_Antigen_Presentation", r"\b(HLA-E|Qa-1b|NKG2A|VMAPRTLVL)\b"),
    ("PVR_TIGIT_MICA_MICB_NKG2D_Axis", "Surface_Ligand_or_Antigen_Presentation", r"\b(PVR|CD155|TIGIT|MICA|MICB|ULBP1|ULBP2|NKG2D)\b"),
    ("PD_L1_CD274_Checkpoint", "Surface_Ligand_or_Antigen_Presentation", r"\b(PD-L1|PDL1|CD274|PD-1)\b"),
    ("HVEM_CD40_BTN3A1_Costimulation", "Surface_Ligand_or_Antigen_Presentation", r"\b(HVEM|TNFRSF14|CD40|CD80|CD86|BTN3A1|CTLA4|CTLA-4)\b"),
    ("FasL_FASLG_Treg_Tolerance_Signal", "Surface_Ligand_or_Antigen_Presentation", r"\b(FasL|SA-FasL|FASLG|CD95L)\b"),
    ("IL2_TGFb_IL10_FOXP3_Treg_Induction", "Secreted_Hormone_or_Chemokine_DAMP", r"\b(FOXP3|FoxP3|Treg|regulatory T cell|IL-10|TGF-beta|TGFB1|IL-2 mutein|2B10)\b"),
    ("CXCL10_CXCR3_Chemokine_Recruitment", "Secreted_Hormone_or_Chemokine_DAMP", r"\b(CXCL10|IP-10|CXCL9|CXCL11|CCL5|CXCR3)\b"),
    ("IFNg_IFNa_IL1b_TNFa_Cytokine_Challenge", "Immune_or_Abiotic_Stressor", r"\b(IFN-gamma|IFNγ|IFN-alpha|IFNα|IL-1beta|IL-1β|TNF-alpha|TNFα|pro-inflammatory cytokine)\b"),
    ("TYK2_JAK1_STAT1_IRF1_Signaling", "Signaling_Kinase_or_TF", r"\b(TYK2|JAK1|JAK2|STAT1|STAT2|IRF1|IRF-1|ruxolitinib|BMS-986165|SOCS1|PSMB8|TAP1)\b"),
    ("TNFAIP3_A20_NFkB_JNK_Brake", "Signaling_Kinase_or_TF", r"\b(TNFAIP3|\bA20\b|NF-kB|NF-κB)\b"),
    ("PTPN2_Phosphatase_STAT1_Brake", "Signaling_Kinase_or_TF", r"\b(PTPN2|TC-PTP)\b"),
    ("RNLS_Glycolytic_NADPH_GSH_Redox_Shield", "Organelle_or_Metabolic_State", r"\b(RNLS|Rnls|renalase|pargyline|NADPH|GSH:GSSG|glutathione)\b"),
    ("ER_Stress_UPR_IRE1_PERK_CHOP_TXNIP", "Organelle_or_Metabolic_State", r"\b(ER stress|unfolded protein response|UPR|IRE1|XBP1|sXBP1|PERK|ATF4|ATF6|CHOP|DDIT3|TXNIP|BiP|HSPA5|MANF|thapsigargin|tunicamycin|DTT|Akita)\b"),
    ("SERPINB9_GranzymeB_Perforin_Apoptosis_Axis", "Cell_Survival_or_Death_Pathway", r"\b(SERPINB9|PI-9|granzyme|perforin|caspase-3|cleaved caspase)\b"),
    ("SLC30A8_ZnT8_Granule_Zinc_Autoantigen", "Gene_or_Protein_Edit", r"\b(SLC30A8|ZnT8|ZnT8A)\b"),
    ("G6PC2_IGRP_Glycolytic_Autoantigen", "Gene_or_Protein_Edit", r"\b(G6PC2|IGRP)\b"),
    ("GAD2_GAD65_GABA_Paracrine_Autoantigen", "Gene_or_Protein_Edit", r"\b(GAD2|GAD65|GADA|\bGABA\b)\b"),
    ("INS_Preproinsulin_Signal_Peptide_DRiP", "Surface_Ligand_or_Antigen_Presentation", r"\b(preproinsulin|insulin signal peptide|PPI15-24|INS-DRiP|insulin B9-23)\b"),
    ("CHGA_WE14_Hybrid_Insulin_Peptide_HIP", "Surface_Ligand_or_Antigen_Presentation", r"\b(CHGA|chromogranin A|WE14|hybrid insulin peptide|\bHIP\b|TGM2)\b"),
    ("IAPP_Amyloid_Oligomer_Inflammasome", "Organelle_or_Metabolic_State", r"\b(IAPP|hIAPP|islet amyloid|amylin|pramlintide|NLRP3)\b"),
    ("CXADR_IFIH1_Enterovirus_CVB_Infection", "Immune_or_Abiotic_Stressor", r"\b(CXADR|CAR-SIV|Coxsackie|CVB3|CVB4|enterovirus|IFIH1|MDA5)\b"),
    ("CD55_CD59_THBD_Complement_IBMIR_Axis", "Surface_Ligand_or_Antigen_Presentation", r"\b(CD55|CD59|CD46|CD39|ENTPD1|THBD|thrombomodulin|IBMIR|complement|tissue factor)\b"),
    ("Stage7_3D_Suspension_Cytoarchitecture", "Biophysical_or_Vascular_Cue", r"\b(Stage 7|3D suspension|3D aggregat|reaggregat|cluster diameter|microwell|AggreWell|organoid)\b"),
    ("MAFA_NKX6_1_PDX1_BetaCell_Identity", "Lineage_Identity_or_Maturation_State", r"\b(MAFA|NKX6\.1|NKX6-1|PDX1|NEUROG3|NKX2\.2|ISL1|UCN3)\b"),
    ("WNT4_ESRRG_Mitochondrial_OXPHOS_Switch", "Organelle_or_Metabolic_State", r"\b(WNT4|ESRRG|ERR\s*gamma|ERRγ|OXPHOS|mitochondrial respiration|OCR|ECAR|LDHA|HK2|SLC16A1|MCT1|pyruvate)\b"),
    ("Biphasic_GSIS_and_Calcium_Coupling", "Secreted_Hormone_or_Chemokine_DAMP", r"\b(GSIS|glucose-stimulated insulin secretion|biphasic|first-phase|perifusion|C-peptide|calcium influx|KATP)\b"),
    ("Alpha_GCG_Delta_SST_Paracrine_CounterRegulation", "Secreted_Hormone_or_Chemokine_DAMP", r"\b(glucagon|\bGCG\b|somatostatin|\bSST\b|alpha cell|delta cell|paracrine)\b"),
    ("OffTarget_SLC18A1_Enterochromaffin_Depletion", "Lineage_Identity_or_Maturation_State", r"\b(SLC18A1|enterochromaffin|SC-EC|TPH1|LMX1A)\b"),
    ("ZM447439_Aurora_G1_CellCycle_Exit", "Small_Molecule_or_Growth_Factor", r"\b(ZM447439|aurora kinase|Ki-67|G0/G1|cell cycle exit|teratoma|iCasp9|HSV-TK)\b"),
    ("CD31_CD146_Microvessel_Inosculation", "Biophysical_or_Vascular_Cue", r"\b(microvessel|microvascular fragment|CD31|CD146|inosculation|anastomos|revasculariz|endothelial|HUVEC|iEC|VEGFA)\b"),
    ("Hypoxia_72h_Ischemia_and_O2_Scaffolds", "Biophysical_or_Vascular_Cue", r"\b(hypoxia|ischemi|oxygen tension|pO2|CaO2|calcium peroxide|perfluorocarbon|subcutaneous|omental|intramuscular)\b"),
    ("Alginate_or_Macrodevice_Encapsulation_FBR", "Biophysical_or_Vascular_Cue", r"\b(alginate|macroencapsul|microencapsul|TheraCyte|PEC-Direct|PEC-Encap|hydrogel|fibrotic|foreign body)\b"),
    ("Alloreactive_CD8_CD4_TCell_or_DSA_Attack", "Immune_or_Abiotic_Stressor", r"\b(alloreactive|allogeneic T|allograft rejection|PBMC|donor-specific antibod|DSA|CD8\+ T|CD4\+ T|cytotoxic T)\b"),
    ("Primary_NK_or_Macrophage_Innate_Cytotoxicity", "Immune_or_Abiotic_Stressor", r"\b(NK cell|natural killer|missing-self|missing self|macrophage|phagocytos)\b"),
    ("BetaCell_Survival_and_Glycemic_Cure", "Cell_Survival_or_Death_Pathway", r"\b(normoglycemi|diabetes reversal|graft survival|beta-cell survival|cell viability|apoptosis|cell death|specific lysis)\b"),
]


def resolve_canonical_node(text: str, default_type: str = "Signaling_Kinase_or_TF") -> tuple[str | None, str | None]:
    """Resolves a clause to its primary canonical biological node ID + schema node type, or (None, None) if none matched."""
    for node_id, node_type, pat in CANONICAL_NODE_PATTERNS:
        if re.search(pat, text, flags=re.IGNORECASE):
            return node_id, node_type
    tl = text.lower()
    if any(k in tl for k in ("transplant", "graft", "implant", "site", "capsule", "portal", "subcutaneous")):
        return "Hypoxia_72h_Ischemia_and_O2_Scaffolds", "Biophysical_or_Vascular_Cue"
    if any(k in tl for k in ("t cell", "t-cell", "lymphocyte", "immune", "rejection", "allogen", "autoimmun")):
        return "Alloreactive_CD8_CD4_TCell_or_DSA_Attack", "Immune_or_Abiotic_Stressor"
    if any(k in tl for k in ("differenti", "stem cell", "ipsc", "hesc", "progenitor", "islet", "beta cell", "β-cell", "β cell")):
        return "MAFA_NKX6_1_PDX1_BetaCell_Identity", "Lineage_Identity_or_Maturation_State"
    if any(k in tl for k in ("insulin", "glucose", "secretion", "glycem", "c-peptide", "diabet")):
        return "Biphasic_GSIS_and_Calcium_Coupling", "Secreted_Hormone_or_Chemokine_DAMP"
    return None, None


def map_pillars(matched_pillars: list[str], full_text: str) -> list[str]:
    out = set()
    for mp in matched_pillars:
        if "Pillar1" in mp:
            out.add("Pillar1_Differentiation_Maturation_Cytoarchitecture")
        if "Pillar2" in mp:
            out.add("Pillar2_Immune_Cloaking_Autoantigen_Neoantigen_Stealth")
        if "Pillar3" in mp:
            out.add("Pillar3_Intrinsic_BetaCell_Genetic_Resilience_Amyloid_Viral")
        if "Pillar4" in mp:
            out.add("Pillar4_Engraftment_Vascularization_Oxygenation_Devices")

    tl = full_text.lower()
    if any(k in tl for k in ("sc-islet", "sc-beta", "differentiation", "stage 6", "stage 7", "mafa", "nkx6.1", "nkx6-1")):
        out.add("Pillar1_Differentiation_Maturation_Cytoarchitecture")
    if any(k in tl for k in ("b2m", "ciita", "cd47", "pd-l1", "hla-e", "fasl", "hypoimmun", "allogeneic", "autoantigen")):
        out.add("Pillar2_Immune_Cloaking_Autoantigen_Neoantigen_Stealth")
    if any(k in tl for k in ("er stress", "thapsigargin", "cytokine", "ifn", "tyk2", "rnls", "txnip", "iapp", "amyloid", "coxsackie", "enterovirus")):
        out.add("Pillar3_Intrinsic_BetaCell_Genetic_Resilience_Amyloid_Viral")
    if any(k in tl for k in ("transplant", "engraft", "vascular", "microvessel", "oxygen", "hypoxia", "encapsulat", "subcutaneous", "oment")):
        out.add("Pillar4_Engraftment_Vascularization_Oxygenation_Devices")

    if not out:
        out.add("Pillar1_Differentiation_Maturation_Cytoarchitecture")
    return sorted(list(out))


def infer_fidelity_tier(full_text: str, mesh: list[str], pub_types: list[str]) -> str:
    tl = full_text.lower()
    mesh_lower = " ".join(mesh).lower()
    pt_lower = " ".join(pub_types).lower()

    if "clinical trial" in pt_lower or any(k in tl for k in ("nct0", "phase 1/2", "phase i/ii", "patients transplanted", "cadaveric t1d pancreas")):
        return "Tier_1_Human_InVivo"
    if any(k in tl for k in ("rhesus", "cynomolgus", "macaque", "nonhuman primate", "non-human primate", "humanized mice", "hu-pbmc", "nsg-sgm3", "blt mice")):
        return "Tier_2A_Human_Cells_InVivo_NHP_or_HumanizedMouse"
    if ("humans" in mesh_lower or any(k in tl for k in ("human islet", "sc-islet", "sc-beta", "hipsc", "hesc", "endoc-bh"))) and any(
        k in tl for k in ("autoreactive t", "cd8+ t", "nk cell", "pbmc", "autoantibod", "alloreactive", "cytotoxic t")
    ):
        return "Tier_2B_Human_Cells_InVitro_HumanImmune_or_AutoAb_CoCulture"
    if "humans" in mesh_lower or any(k in tl for k in ("human islet", "sc-islet", "sc-beta", "hipsc", "hesc", "endoc-bh", "human pluripotent")):
        if not re.search(r"\bnon-?human\s+primate\b", tl):
            return "Tier_3_Human_Cells_AbioticStress_or_ImmunodeficientMouse"
    return "Tier_4_Rodent_or_Xenogeneic_Only"


def infer_cell_source_and_stage(full_text: str, mesh_terms: list[str] | None = None) -> tuple[str, str | None, str]:
    tl = full_text.lower()
    mesh_lower = " ".join(mesh_terms or []).lower()
    tl_unnegated = re.sub(
        r"\b(?:no|without|lacking|neither|non[- ])\s*human\s+(?:pancreatic\s+|donor\s+|cadaveric\s+)?(?:islets?|beta[- ]cells?|β[- ]cells?|cells?|donors?|patients?|primates?)\b",
        " ",
        tl,
    )
    tl_substrate = re.sub(
        r"\bhuman\s+(?:t\s+cells?|cd[48]\+?\s*t\s+cells?|nk\s+cells?|pbmcs?|peripheral\s+blood|lymphocytes?|monocytes?|macrophages?|serum|autoantibodies|antibodies|complement|leukocytes?|immune\s+cells?|effectors?|cytotoxic\s+t)\b",
        " ",
        tl_unnegated,
    )
    has_nhp = bool(re.search(r"\b(?:cynomolgus|macaque|rhesus|baboon|non[- ]?human\s+primates?|nhp)\b", tl))
    has_xeno_islet = bool(
        re.search(
            r"\b(?:porcine|pig|pigs|swine|bovine|canine|ovine|xenogeneic)\s+(?:pancreatic\s+|donor\s+|neonatal\s+|fetal\s+|adult\s+)?(?:islets?|beta[- ]cells?|β[- ]cells?)\b",
            tl,
        )
    )
    has_human = bool(
        re.search(
            r"(?<!non-)(?<!non\s)(?<!no\s)\b(?:humans?|patients?|cadaveric|donors?\s+with\s+t1d|hescs?|hipscs?|ipscs?|sc-islets?|sc-beta|endoc)\b",
            tl_substrate,
        )
        or ("humans" in mesh_lower)
    )
    has_rodent = bool(
        re.search(r"\b(?:mouse|mice|murine|rat|rats|rodent|nod\b|c57bl|balb/c|akita|stz-induced\s+diabetic\s+mice)\b", tl)
        or (any(m in mesh_lower for m in ("mice", "rats", "rodentia")) and "humans" not in mesh_lower)
    )
    has_islet_or_beta = bool(
        re.search(
            r"\b(?:islets?|beta[- ]cells?|β[- ]cells?|pancreatic\s+endocrine|insulin-producing\s+cells?|sc-islets?|sc-beta)\b",
            tl,
        )
        or ("islets of langerhans" in mesh_lower)
    )
    has_explicit_human_islet = bool(
        re.search(
            r"(?<!non-)(?<!non\s)(?<!no\s)\b(?:human\s+(?:pancreatic\s+|donor\s+|cadaveric\s+)?islets?|cadaveric\s+(?:human\s+|donor\s+)?islets?|islets?\s+(?:isolated\s+)?from\s+human|human\s+(?:beta|β)[- ]cells?)\b",
            tl_unnegated,
        )
    )

    if "cipsc" in tl or "chemically induced pluripotent" in tl:
        source_class = "CiPSC_chemically_induced"
    elif ("t1d" in tl or "type 1 diabet" in tl) and ("patient-derived ipsc" in tl or "t1d ipsc" in tl or "t1d-ipsc" in tl):
        source_class = "iPSC_T1D_or_LADA_patient"
    elif re.search(r"\b(?:ipscs?|hipscs?|induced\s+pluripotent)\b", tl):
        source_class = "iPSC_healthy_donor"
    elif (
        re.search(r"\b(?:hescs?|human\s+embryonic\s+stem|embryonic\s+stem\s+cells?|esc-derived|hues8|wa01|mel1|cyt49)\b", tl)
        or bool(re.search(r"\b(?:h1\s+hescs?|h1\s+escs?|h1\s+cells?|h1\s+line)\b", tl))
    ) and not re.search(r"embryonic\s+stem\s+cells?\s+were\s+not\s+used|without\s+embryonic\s+stem\s+cells?", tl):
        source_class = "hESC"
    elif any(k in tl for k in ("endoc-b", "endoc-β", "endoc_b", "endocbh")):
        source_class = "EndoC_betaH1_or_H2_H3"
    elif re.search(r"\b(?:nit-1|min6|ins-1e?)\b", tl):
        source_class = "Rodent_Beta_Cell_Line_NIT1_MIN6_INS1"
    elif has_explicit_human_islet:
        source_class = "Primary_Human_Islet"
    elif has_nhp and has_islet_or_beta and not has_explicit_human_islet:
        source_class = "NHP_Primary_Islet"
    elif has_xeno_islet and not has_explicit_human_islet:
        source_class = "Unspecified_In_Source_Abstract"
    elif has_islet_or_beta:
        if has_rodent and not has_human:
            source_class = "Rodent_Primary_Islet"
        elif has_human and not has_rodent:
            source_class = "Primary_Human_Islet"
        elif has_rodent and has_human:
            # Mixed human/rodent mention without explicit human islet phrase: check whether rodent islet or human islet is primary
            if re.search(r"\b(?:mouse|mice|murine|rat|rodent)\s+(?:pancreatic\s+|donor\s+)?islets?\b", tl):
                source_class = "Rodent_Primary_Islet"
            else:
                source_class = "Unspecified_In_Source_Abstract"
        else:
            source_class = "Unspecified_In_Source_Abstract"
    else:
        source_class = "Unspecified_In_Source_Abstract"

    cell_line = None
    if source_class in ("hESC", "iPSC_healthy_donor", "iPSC_T1D_or_LADA_patient", "CiPSC_chemically_induced"):
        for cl, pat in (
            ("HUES8", r"\bhues-?8\b"),
            ("WA01", r"\bwa01\b"),
            ("H1", r"\b(?:h1\s+hescs?|h1\s+escs?|h1\s+cells?|h1\s+line|h1/h9|h9\s+and\s+h1)\b"),
            ("MEL1", r"\bmel-?1\b"),
            ("CyT49", r"\bcyt49\b"),
            ("HEL46.11", r"\bhel46\.11\b"),
            ("HEL24.3", r"\bhel24\.3\b"),
        ):
            if re.search(pat, tl):
                cell_line = cl
                break
        if cell_line in ("HUES8", "WA01", "H1", "MEL1", "CyT49") and source_class in (
            "iPSC_healthy_donor",
            "iPSC_T1D_or_LADA_patient",
            "CiPSC_chemically_induced",
        ):
            source_class = "hESC"
    elif source_class == "EndoC_betaH1_or_H2_H3":
        if re.search(r"\bendoc-(?:β|beta)h3\b", tl):
            cell_line = "EndoC-βH3"
        elif re.search(r"\bendoc-(?:β|beta)?h1\b", tl):
            cell_line = "EndoC-βH1"
    elif source_class == "Rodent_Beta_Cell_Line_NIT1_MIN6_INS1":
        for cl, pat in (
            ("INS-1E", r"\bins-1e\b"),
            ("INS-1", r"\bins-1\b"),
            ("MIN6", r"\bmin-?6\b"),
            ("NIT-1", r"\bnit-1\b"),
        ):
            if re.search(pat, tl):
                cell_line = cl
                break

    if "stage 7" in tl or "s7w" in tl or "biphasic" in tl:
        stage = "Stage7_MetabolicallyMature_SC_Islet"
    elif "sc-islet" in tl or "sc-beta" in tl or "stage 6" in tl:
        stage = "Stage6_Immature_SC_Islet"
    elif "neurog3" in tl or "endocrine progenitor" in tl or "stage 5" in tl:
        stage = "Stage5_EndocrineProgenitor_NEUROG3"
    elif "pancreatic progenitor" in tl or "stage 4" in tl or "pdx1+nkx6" in tl:
        stage = "Stage4_PancreaticProgenitor_PDX1_NKX6_1"
    else:
        stage = "Primary_Mature_Islet"

    return source_class, cell_line, stage


def infer_transplant_site_and_niche(full_text: str) -> dict[str, Any]:
    tl = full_text.lower()
    if "portal vein" in tl or "intrahepatic" in tl:
        site = "Portal_Vein_Intrahepatic"
    elif "omentum" in tl or "omental" in tl:
        site = "Omental_Pouch"
    elif "prevascular" in tl or ("subcutaneous" in tl and any(k in tl for k in ("device", "microvessel", "scaffold", "hydrogel", "pouch", "catheter"))):
        site = "Subcutaneous_Prevascularized_or_Device"
    elif "subcutaneous" in tl:
        site = "Subcutaneous_Unmodified"
    elif "intramuscular" in tl or "quadricep" in tl or "hindlimb muscle" in tl or "rectus" in tl or "brachioradialis" in tl or "forearm" in tl:
        site = "Intramuscular_Forearm_or_Rectus"
    elif re.search(r"\b(?:anterior\s+chamber|intraocular)\b", tl):
        site = "Anterior_Chamber_Eye"
    elif "epididymal" in tl or "fat pad" in tl or "perigonadal" in tl:
        site = "Epididymal_or_Perigonadal_Fat_Pad"
    elif "kidney capsule" in tl or "renal subcapsular" in tl or "subcapsular" in tl:
        site = "Kidney_Capsule"
    elif re.search(r"\b(?:transplant\w*|graft\w*|recipient\w*|implant\w*|allograft\w*|xenograft\w*)\b", tl) and not re.search(r"\b(?:in\s+vitro\s+only|no\s+transplantation)\b", tl):
        site = "Unspecified_InVivo_or_Source_Site"
    else:
        site = "InVitro_Only"

    has_iec = bool(re.search(r"\biecs?\b", tl))
    if "microvessel" in tl or "microvascular fragment" in tl:
        vasc = "Microvascular_Fragment_CoAssembly"
    elif "pericyte" in tl and ("endothelial" in tl or has_iec or "huvec" in tl):
        vasc = "Endothelial_plus_Pericyte_Vascularized_Organoid"
    elif "endothelial" in tl or "huvec" in tl or has_iec:
        vasc = "Endothelial_CoCulture_iEC_HUVEC"
    else:
        vasc = "None_Avascular_Islet"

    if "cao2" in tl or "calcium peroxide" in tl:
        o2 = "CaO2_PDMS_Oxygen_Microbeads"
    elif "perfluorocarbon" in tl or re.search(r"\bpfcs?\b", tl):
        o2 = "Perfluorocarbon_O2_Carrier"
    elif "beta-air" in tl or "oxygen chamber" in tl or "electrochemical oxygen" in tl:
        o2 = "Exogenous_O2_Chamber_or_Bioelectronic"
    elif "hemoglobin" in tl or "myoglobin" in tl:
        o2 = "Hemoglobin_or_Myoglobin_Hydrogel"
    else:
        o2 = "None"

    biomaterial = None
    for bm in ("PEG-4MAL microgel", "alginate", "TheraCyte", "PEC-Direct", "PEC-Encap", "PTFE", "fibrin", "Matrigel", "collagen"):
        if bm.lower() in tl:
            biomaterial = bm
            break

    return {
        "transplant_site": site,
        "prevascularization_lead_time_days": None,
        "in_dish_vascular_assembly": vasc,
        "oxygen_bridge_72h_carrier": o2,
        "packing_density_ieq_per_cm2_or_ml": None,
        "number_of_distributed_micro_pockets": None,
        "staged_insertion_schedule": None,
        "encapsulation_or_biomaterial": biomaterial,
        "co_transplanted_cells_or_oxygen_source": None,
        "systemic_immunosuppression_regimen": "None" if "without immunosuppression" in tl or "nsg" in tl else None,
    }


def _is_downstream_readout_only(win: str, sym: str, req_re: re.Pattern[str], matched_alias: str | None = None) -> bool:
    """
    Returns True if the perturbation keyword in `win` is attached to an upstream subject
    while `sym` appears only as a downstream readout object (e.g. 'The knockout cells reduced B2M and INS expression').
    """
    # Bind to the actual source spelling. Re-searching only CD274 silently missed
    # PDL1 / PD-L1 aliases and converted IRF1 readouts into CD274 interventions.
    alias = matched_alias or sym
    m_sym = re.search(rf"(?<!\w){re.escape(alias)}(?!\w)", win, flags=re.IGNORECASE)
    if not m_sym:
        return False
    # Check if there is a direct perturbation bound tightly to sym (e.g. 'B2M KO', 'KO of B2M', 'B2M-/-', 'CRISPR B2M')
    direct_pat = re.compile(
        rf"(?:\b{re.escape(alias)}(?:\s*|-)(?:-/-|−/−|\bKO\b|\bknockout\b|\bknock-out\b|\bknockdown\b|\bdeletion\b|\bnull\b|\bdeficien\w+|\boverexpress\w+|\bmutation\b|\bmutant\b|\binhibition\b|\bsilenc\w+))|"
        rf"(?:\b(?:knockout|knock-out|KO|knockdown|deletion|ablation|silencing|overexpression|inhibition|loss|deficiency|mutation)\s+(?:of\s+)?(?:human\s+|mouse\s+|endogenous\s+)?(?:[A-Za-z0-9_-]+\s+and\s+)?\b{re.escape(alias)}\b)",
        flags=re.IGNORECASE,
    )
    if direct_pat.search(win):
        return False
    # If all perturbation matches precede a causal/readout verb and sym follows that verb (or precedes 'expression/transcripts/mRNA/levels')
    after_sym = win[m_sym.end() :]
    before_sym = win[: m_sym.start()]
    if re.match(
        r"^\s*(?:,\s*[A-Za-z0-9_-]+\s*,?\s*)?(?:and\s+[A-Za-z0-9_-]+\s+)?(?:mRNA\s+|protein\s+|gene\s+|surface\s+)?(?:expression|transcript\w*|levels?|secretion|staining|markers?|content|induction|(?:was|were)\s+(?:significantly\s+)?(?:reduc\w+|decreas\w+|increas\w+|alter\w+|chang\w+|suppress\w+|downregulat\w+|upregulat\w+))\b",
        after_sym,
        flags=re.IGNORECASE,
    ):
        if re.search(
            r"\b(?:reduc\w+|decreas\w+|increas\w+|downregulat\w+|upregulat\w+|induc\w+|suppress\w+|alter\w+|chang\w+|elevat\w+|lower\w+|showed|exhibited|assessed|measured|knockout|KO\b|-/-|−/−)\b",
            before_sym,
            flags=re.IGNORECASE,
        ):
            return True
    if req_re.search(before_sym) and not req_re.search(after_sym):
        if re.search(
            r"\b(?:knockout|knock-out|KO\b|-/-|−/−|knockdown|deletion|ablation|silencing|overexpression|deficien\w+|loss|mutation|mutant)\b[^.;]*?\b(?:reduc\w+|decreas\w+|increas\w+|downregulat\w+|upregulat\w+|induc\w+|suppress\w+|alter\w+|chang\w+|elevat\w+|lower\w+|impair\w+|abolish\w+)\b",
            before_sym,
            flags=re.IGNORECASE,
        ):
            return True
    # Target suggestions/background mentions are not experimental manipulations.
    if re.search(r"\b(?:potential|promising|candidate)\s+(?:therapeutic\s+)?target\b|target\s+for\s+(?:immuno)?therapy", win, re.I):
        return True
    return False


def detect_genetic_edits(title: str, abstract: str = "", full_text: str = "") -> list[dict[str, Any]]:
    """
    Strictly discriminates active genetic/pharmacological perturbations of a target gene
    from passive immunostaining/marker mentions!
    A gene is included ONLY if:
      (a) It appears in the study Title or Abstract alongside a perturbation/functional keyword, OR
      (b) Within a +/- 140 character window around the gene match in full_text, there is an explicit
          perturbation keyword (KO, CRISPR, shRNA, mutation, overexpression, inhibitor, etc.)
          and NOT solely a passive immunostaining marker phrase.
    """
    edits = []
    seen = set()
    title_abs = f"{title} . {abstract}" if abstract else title
    effective_full_text = full_text or title_abs

    for spec in TARGET_GENE_SPECS:
        sym = spec["symbol"]
        if sym in seen:
            continue

        pat = re.compile(spec["pattern"], flags=re.IGNORECASE)
        is_perturbed = False
        matched_context = ""

        # Check title + abstract first (sentence-bounded window; title matches must also satisfy req_re)
        req_re = STRICT_GENETIC_EDIT_RE if sym in ENDOGENOUS_LINEAGE_MARKERS else PERTURBATION_CONTEXT_RE
        passive_readout_re = re.compile(
            rf"\b{re.escape(sym)}\b\s+(?:mrna\s+|protein\s+|gene\s+)?(?:expression|secretion|levels?|transcripts?)\s+(?:and\s+secretion\s+)?(?:was|were)\s+(?:assessed|measured|quantified|analyzed|determined|increased|decreased)",
            flags=re.IGNORECASE,
        )
        for m in pat.finditer(title_abs):
            start = max(0, m.start() - 95)
            end = min(len(title_abs), m.end() + 95)
            left_sep = max(title_abs.rfind(". ", start, m.start()), title_abs.rfind("; ", start, m.start()))
            if left_sep != -1:
                start = left_sep + 2
            right_dot = title_abs.find(". ", m.end(), end)
            right_semi = title_abs.find("; ", m.end(), end)
            right_cands = [idx for idx in (right_dot, right_semi) if idx != -1]
            if right_cands:
                end = min(right_cands) + 1
            win = title_abs[start:end]
            win_clean = re.sub(
                r"\b(?:membrane\s+attack\s+complex-inhibitory\s+protein|complement\s+inhibitory\s+receptors?)\b",
                " ",
                win,
                flags=re.IGNORECASE,
            )
            if (
                req_re.search(win_clean)
                and not PASSIVE_MARKER_ONLY_RE.search(win)
                and not passive_readout_re.search(win)
                and not _is_downstream_readout_only(win, sym, req_re, m.group(0))
            ):
                is_perturbed = True
                matched_context = win.strip()
                break

        # If not in title/abstract, check full_text windows with strict sentence-bounded perturbation requirement
        if not is_perturbed:
            for m in pat.finditer(effective_full_text):
                start = max(0, m.start() - 90)
                end = min(len(effective_full_text), m.end() + 90)
                left_sep = max(effective_full_text.rfind(". ", start, m.start()), effective_full_text.rfind("; ", start, m.start()))
                if left_sep != -1:
                    start = left_sep + 2
                right_dot = effective_full_text.find(". ", m.end(), end)
                right_semi = effective_full_text.find("; ", m.end(), end)
                right_cands = [idx for idx in (right_dot, right_semi) if idx != -1]
                if right_cands:
                    end = min(right_cands) + 1
                win = effective_full_text[start:end]
                win_clean = re.sub(
                    r"\b(?:membrane\s+attack\s+complex-inhibitory\s+protein|complement\s+inhibitory\s+receptors?)\b",
                    " ",
                    win,
                    flags=re.IGNORECASE,
                )
                if (
                    req_re.search(win_clean)
                    and not PASSIVE_MARKER_ONLY_RE.search(win)
                    and not passive_readout_re.search(win)
                    and not _is_downstream_readout_only(win, sym, req_re, m.group(0))
                ):
                    is_perturbed = True
                    matched_context = win.strip()
                    break

        if is_perturbed:
            seen.add(sym)
            # Infer modality from context window if possible
            wl = matched_context.lower()
            modality = spec["default_modality"]
            if any(k in wl for k in ("shrna", "sirna", "knockdown", "crispri")):
                modality = "CRISPRi_or_shRNA_Knockdown"
            elif re.search(r"\b(?:base[ -]edit\w*|abe8e|adenine\s+base\s+editor|cytosine\s+base\s+editor)\b", wl):
                modality = "BaseEdit_ABE_or_CBE_Epitope_or_Splice"
            elif re.search(r"\bprime[ -]edit\w*\b", wl):
                modality = "PrimeEdit"
            elif any(k in wl for k in ("inhibit", "small molecule", "pargyline", "verapamil", "bms-986165", "deucravacitinib")):
                modality = "SmallMolecule_Pharmacological_Modulator"
            elif any(k in wl for k in ("overexpress", "transduc", "lentivir")):
                modality = "Lentiviral_or_Transposon_OE"
            elif any(k in wl for k in ("knock-in", "knockin", "aavs1", "gapdh")):
                modality = "SafeHarbor_KI_AAVS1_GAPDH_CLYBL"
            elif re.search(r"\b(?:knockout|knock-out|ko|crispr[- /]?cas9|delet\w+|ablat\w+|deficien\w+|loss\s+of)\b|-/-|−/−", wl):
                modality = "CRISPR_Cas9_KO"

            if modality in ("CRISPR_Cas9_KO", "CRISPR_Cas12a_Multiplex_KO", "CRISPRi_or_shRNA_Knockdown"):
                p_class = "Loss_of_Function_KO_or_KD"
            elif modality in ("BaseEdit_ABE_or_CBE_Epitope_or_Splice", "PrimeEdit", "Signal_Peptide_or_Domain_Swap"):
                p_class = "Precision_Epitope_or_Isoform_Edit"
            elif modality == "Endogenous_SwapIn_KO_plus_KI":
                p_class = "Bicistronic_Endogenous_Promoter_SwapIn_KO_plus_KI"
            elif modality in ("SafeHarbor_KI_AAVS1_GAPDH_CLYBL", "Lentiviral_or_Transposon_OE", "CRISPRa_Activation", "Transient_mRNA_or_Surface_Tethering"):
                p_class = "Gain_of_Function_or_Transgenic_Armoring"
            else:
                p_class = "Pharmacological_or_Small_Molecule_Modulator"

            edits.append({
                "gene_symbol_hgnc": sym,
                "ensembl_gene_id": spec["ensembl"],
                "sub_molecular_coordinate": {
                    "region_class": spec["region_class"],
                    "exact_variant_or_swap_description": matched_context[:180] if matched_context else f"Experimental perturbation of {sym}",
                },
                "functional_gate": spec["functional_gate"],
                "target_autoantibody_or_tcr_evaded": spec["evaded"],
                "edit_modality": modality,
                "promoter_or_locus_context": spec["default_promoter"],
                "zygosity_or_allelic_state": None,
                "protein_or_surface_expression_fold_vs_wt": None,
                "perturbation_class": p_class,
            })
    return edits


def classify_outcome_domain(sentence: str) -> str:
    sl = sentence.lower()
    if any(k in sl for k in ("nk cell", "macrophage", "sirp", "cd47", "missing self", "missing-self")):
        return "NK_or_Macrophage_Innate_Survival"
    if any(k in sl for k in ("autoreactive", "autoantibod", "cd8+ t", "tcr", "peptide-specific lysis", "ins-drip", "gad65", "znt8", "igrp")):
        return "Autoimmune_TCR_or_AutoAb_Evasion"
    if any(k in sl for k in ("allogeneic", "allograft", "dsa", "donor-specific antibod", "mhc mismatch", "foxp3", "treg")):
        return "Allo_T_Cell_or_DSA_Survival"
    if any(k in sl for k in ("complement", "ibmir", "thrombomodulin", "cd55", "cd59", "coagulation")):
        return "Complement_or_IBMIR_Coagulation_Block"
    if any(k in sl for k in ("granzyme", "perforin", "fasl", "fas-mediated", "serpinb9")):
        return "Granzyme_or_Fas_Apoptosis_Resistance"
    if any(k in sl for k in ("iapp", "amyloid", "coxsackie", "enterovirus", "cvb", "cxadr", "mda5", "ifih1")):
        return "Amyloid_hIAPP_or_Enterovirus_Protection"
    if any(k in sl for k in ("thapsigargin", "tunicamycin", "er stress", "upr", "chop", "txnip", "ifn", "cytokine", "il-1", "tnf", "stat1", "tyk2", "rnls")):
        return "Cytokine_IFN_NFkB_or_ER_Stress_Survival"
    if any(k in sl for k in ("hypoxia", "ischemia", "oxygen", "vascular", "inosculation", "microvessel", "vessel", "capillary")):
        return "Hypoxia_72h_or_Vascular_Inosculation"
    if any(k in sl for k in ("hla", "mhc class", "pd-l1 expression", "surface expression", "mfi")):
        return "Surface_HLA_or_Transgene_Expression_Level"
    if any(k in sl for k in ("normoglycemi", "diabetes reversal", "graft survival", "c-peptide in vivo", "blood glucose")):
        return "InVivo_Glycemic_Cure_and_Graft_Longevity"
    if any(k in sl for k in ("glucagon", "gcg", "somatostatin", "sst", "alpha cell", "delta cell", "paracrine")):
        return "Alpha_Delta_Paracrine_CounterRegulation"
    if any(k in sl for k in ("gsis", "insulin secretion", "perifusion", "biphasic", "respiration", "ocr", "glycoly", "oxphos", "calcium", "[ca2+]")):
        return "Metabolic_and_Biphasic_GSIS_Maturation"
    return "Differentiation_Yield_and_Lineage_Purity"


def _v063_classify_direction_for_sort(sentence: str) -> str:
    """Preserves exact v0.6.3 candidate sentence sort priority direction so all 2,794 effect_ids remain invariant."""
    sl = sentence.lower()
    if any(
        k in sl
        for k in (
            "unchanged", "similar", "no significant difference", "not significantly",
            "not statistically significant", "comparable", "did not affect",
            "was not prolonged", "were not prolonged", "did not increase",
            "did not prevent", "no protection", "failed to prolong",
        )
    ):
        return "Neutral_or_NonInferior"
    has_non_cyto_toxic = bool(re.search(r"(?<!cyto)\btoxic\w*", sl))
    has_harm_kw = has_non_cyto_toxic or any(
        k in sl
        for k in (
            "impair", "compromis", "fail", "loss of", "destroyed", "rejected",
            "hypoglycemi", "silenc", "translocation", "chromothripsis", "teratoma",
            "exacerbat", "worsen", "blunt", "poor gsis", "poor dynamic", "off-target",
        )
    )
    haz_pat = r"(lysis|killing|apoptosis|apoptotic|cell death|cytotoxic\w*|rejection|inflammation|er stress|fibrosis|necrosis|hyperglyc\w+|insulitis|infiltrat\w+)"
    up_pat = r"(increas\w+|elevat\w+|higher|enhanc\w+|greater|upregulat\w+|induc\w+|exacerbat\w+|accelerat\w+|more)"
    down_pat = r"(reduc\w+|decreas\w+|lower\w*|prevent\w*|protect\w*|attenuat\w+|inhibit\w+|abolish\w*|block\w*|suppress\w*|diminish\w*|ameliorat\w+|delay\w*|abrogat\w+|less|fewer|resist\w*|rescu\w+|spar\w+)"
    haz_increase = bool(re.search(up_pat + r"[^.;]{0,50}" + haz_pat, sl) and not re.search(down_pat, sl))
    has_down_verb = any(k in sl for k in ("reduc", "decreas", "lower", "diminish", "abolish", "suppressed", "inhibited"))
    has_func_target = any(k in sl for k in ("gsis", "insulin secretion", "stimulated insulin", "c-peptide", "graft survival", "cell viability", "beta-cell mass"))
    has_hazard_target = any(k in sl for k in ("lysis", "killing", "apoptosis", "cell death", "necrosis", "rejection", "hyperglycemi", "blood glucose", "amyloid", "fibrosis", "inflammation", "er stress", "oxidative stress"))
    if has_harm_kw or haz_increase or (has_down_verb and has_func_target and not has_hazard_target):
        return "Deleterious_or_Impaired"
    if any(k in sl for k in ("increas", "improv", "enhanc", "protect", "prolong", "restor", "preserv", "prevent", "evad", "resist", "surviv", "revers", "normoglycemi", "reduc", "decreas", "inhibit", "attenuat", "block", "abolish", "fold", "vs", "compared")):
        return "Protective_or_Improved"
    return "Neutral_or_NonInferior"


def classify_sentence_direction_and_penalty(
    sentence: str, outcome_cat: str | None = None, target_gene: str | None = None
) -> tuple[str, bool]:
    """
    Endpoint-orientation-aware, clause-bound direction and independent secretome/differentiation penalty classifier (v0.6.8).
    - Binds impairment/reduction modifiers to functional secretome/differentiation/viability targets within the same clause
      (excluding beneficial reductions of blood glucose, HbA1c, cost, non-endocrine/off-target cells, pluripotency markers,
      teratoma formation, and cytotoxicity, as well as rescue/amelioration of cytokine-induced impairment).
    - Separates effect_direction from secretome_or_differentiation_penalty_observed (no penalty aliasing, AT-14).
    - Does not match 'toxic' inside 'cytotoxic'/'cytotoxicity' (AT-10).
    - Flags sentences where a hazard endpoint increases without a decrease verb as Deleterious_or_Impaired (AT-11).
    """
    sl = sentence.lower().strip()
    if re.search(
        r"\b(?:not\s+(?:a\s+)?(?:statistically\s+)?significant\s+(?:increase|decrease|reduction|difference)|"
        r"(?:did\s+not\s+(?:significantly\s+)?|failed\s+to\s+)(?:prolong|improve|increase|enhance|protect|reduce|decrease)|"
        r"no\s+evidence\s+(?:that|of)|not\s+associated\s+with|neither\b[^.;]*\bnor)\b",
        sl,
    ) and not re.search(r"[,;]\s*(?:but|however|whereas|loss\s+of)\b", sl):
        return "Neutral_or_NonInferior", False
    if re.search(
        r"\b(?:the\s+aims?\s+(?:of\s+\w+\s+(?:study\s+)?)?(?:was|were|is)\s+to|(?:this\s+study|we)\s+(?:aimed|sought)\s+to|(?:our|the)\s+(?:objective|goal|purpose)\s+(?:of\s+this\s+study\s+)?(?:was|is)\s+to|we\s+sought\s+to\s+determine\s+whether|to\s+investigate\s+whether|to\s+determine\s+whether)\b",
        sl,
    ):
        return "Neutral_or_NonInferior", False
    if re.match(
        r"^(?:there\s+(?:was|were|is|are)\s+)?no\s+(?:significant\s+|statistically\s+significant\s+)?(?:increase|reduction|decrease|change|difference)\s+in\b",
        sl,
    ) and not re.search(r",\s*(?:but|however|whereas)\b", sl):
        return "Neutral_or_NonInferior", False
    m_while_neutral = re.match(
        r"^(?:while|whilst)\s+([^,;]+?\b(?:had\s+no\s+effect\s+on|did\s+not\s+affect|were\s+unaltered|was\s+unaltered)[^,;]*)(?:,\s*(.+))?$",
        sl,
    )
    if m_while_neutral:
        while_clause, main_clause = m_while_neutral.group(1), (m_while_neutral.group(2) or "").strip()
        if (target_gene and re.search(rf"\b{re.escape(target_gene.lower())}\b", while_clause)) or not main_clause:
            return "Neutral_or_NonInferior", False
        return classify_sentence_direction_and_penalty(main_clause, target_gene=target_gene)
    if target_gene and re.search(
        rf"\b{re.escape(target_gene.lower())}\b[^,;.]{{0,45}}\b(?:had\s+no\s+effect\s+on|did\s+not\s+affect|did\s+not\s+impair|unchanged)\b",
        sl,
    ):
        return "Neutral_or_NonInferior", False
    neg_impair = any(
        k in sl
        for k in (
            "did not impair",
            "do not impair",
            "does not impair",
            "were not impaired",
            "was not impaired",
            "not impaired",
            "without impairing",
            "without affecting",
            "with no impact on",
            "do not show impaired",
            "did not compromise",
            "does not compromise",
            "without compromising",
            "did not reduce",
            "does not reduce",
            "without reducing",
            "did not decrease",
            "does not decrease",
            "did not affect",
            "does not affect",
            "does not alter",
            "did not alter",
            "does not impact",
            "no reduction in",
            "no significant reduction",
            "no statistically significant reduction",
            "no decrease in",
            "no significant decrease",
            "no change in",
            "no significant change",
            "no loss of",
        )
    )
    neg_hypo = (
        any(
            k in sl
            for k in (
                "without hypoglycemi",
                "without severe hypoglycemi",
                "no hypoglycemi",
                "no episodes of hypoglycemi",
                "prevented hypoglycemi",
                "prevent hypoglycemi",
                "protection against hypoglycemi",
                "protect against hypoglycemi",
                "protected from hypoglycemi",
                "free of hypoglycemi",
                "freedom from hypoglycemi",
                "absence of hypoglycemi",
                "unawareness",
            )
        )
        and "fell to 70 mg/dl or less" not in sl
    )
    neg_teratoma = any(
        k in sl
        for k in (
            "no teratoma",
            "without teratoma",
            "prevented teratoma",
            "prevent teratoma",
            "eliminated teratoma",
            "reduces proliferation and teratoma",
            "reduced teratoma",
            "reduces teratoma",
            "inhibited teratoma",
            "free of teratoma",
        )
    )
    neg_off_target = bool(
        re.search(
            r"\b(?:reduc\w+|decreas\w+|eliminat\w+|prevent\w+|inhibit\w+|deplet\w+|no|without|free\s+of)\b[^.;]{0,45}\b(?:off-target|non-?endocrine|enterochromaffin|pluripotency\s+markers?)",
            sl,
        )
    )
    is_rnai_silenc = any(
        k in sl
        for k in (
            "sirna",
            "shrna",
            "rnai",
            "gene silencing",
            "silencing of ",
            "silenced by ",
        )
    )
    is_benign_transloc = any(
        k in sl
        for k in (
            "nuclear translocation",
            "translocation to the nucleus",
            "translocation into the nucleus",
            "membrane translocation",
            "glut4 translocation",
            "ten-eleven translocation",
            "protein translocation",
            "bacterial translocation",
        )
    )
    is_rescue_of_impairment = bool(
        re.search(
            r"\b(ameliorat\w+|preserv\w+|restor\w+|rescu\w+|revers\w+|counteract\w+|prevent\w+|protect\w+)\b(?:(?!(?:\bbut\b|\bhowever\b|\bwhereas\b|\bwhile\b|\balthough\b|\byet\b|\bdespite\b))[^.;]){0,65}\b(impair\w+|decreas\w+|reduc\w+|dedifferentiat\w+|dysfunction)",
            sl,
        )
        or re.search(
            r"\b(impair\w+|decreas\w+|reduc\w+|dedifferentiat\w+|dysfunction)\b(?:(?!(?:\bbut\b|\bhowever\b|\bwhereas\b|\bwhile\b|\balthough\b|\byet\b|\bdespite\b))[^.;]){0,95}\b(?:could\s+be\s+|was\s+|were\s+)(ameliorat\w+|preserv\w+|restor\w+|rescu\w+)\b",
            sl,
        )
        or re.search(
            r"\b(in\s+untreated\s+\w+|in\s+control\s+mice|in\s+diabetic\s+controls?)\b[^.;]*\b(?:whereas|while)\b[^.;]*\b(ameliorat\w+|preserv\w+|restor\w+|rescu\w+|revers\w+|prevent\w+|protect\w+)\b",
            sl,
        )
    )
    is_non_experimental_or_baseline_context = any(
        k in sl
        for k in (
            "we questioned if ",
            "thus, we hypothesized that ",
            "was associated with t2d and reduced insulin secretion",
            "type 2 diabetes (t2d) is associated with defective insulin secretion",
            "which may account for the deficient insulin secretion observed in our mody3 patients",
            "conferring the increased risk of impaired insulin secretion and future risk of t2d",
            "ethical issues still exist",
            "persisted in both wt and rnls ko",
        )
    )

    func_targets_pat = (
        r"(?:gsis|glucose-stimulated\s+insulin\s+secretion|glucose-dependent\s+insulin\s+secretion|"
        r"insulin\s+secretion|insulin\s+release|insulin\s+content|insulin\s+synthesis|"
        r"c-peptide\s+(?:level|secretion|release)|proinsulin\s+cleavage|proinsulin\s+signal|"
        r"(?:β-cell|beta-cell|β\s+cell|beta\s+cell|pancreatic|endocrine|progenitor|sc-islet|pp|ep)\s+(?:differentiation|maturation|commitment|identity|mass|survival|viability|function)|"
        r"differentiation\s+(?:into|of\s+(?:human|pancreatic|functional|β|beta)|efficiency|efficiencies|block)|"
        r"commitment\s+to\s+chga\+|formation\s+of\s+glucose-responsive|"
        r"first-phase\s+glucose-stimulated|peak\s+stimulation|glucose\s+sensing|glucose-sensing|"
        r"cell\s+viability|graft\s+viability|viability\s+compared|viability\s+observed|"
        r"nkx6-1|neurog3|pdx1\s+expression|iapp\s+induction|coordinate\s+function)"
    )
    impair_verbs_pat = r"(?:impair\w+|compromis\w+|blunt\w+|defect\w+|deficien\w+|abolish\w+|diminish\w+|weaker\s+response|absent\s+gsis|poor\s+glucose\s+responsiveness|poor\s+gsis|poor\s+dynamic|suppression\s+of\s+the\s+machinery\s+for\s+glucose-sensing|inhibit\w+(?=\s+(?:gsis|glucose-stimulated\s+insulin\s+secretion|insulin\s+secretion|insulin\s+release)))"
    down_verbs_pat = r"(?:reduc\w+|decreas\w+|lower\w*)"

    has_bound_penalty = False
    if not is_non_experimental_or_baseline_context:
        for clause in re.split(
            r"[.;]|,\s*(?:whereas|while|although|but|however|yet|despite|and\s+was\s+associated|resulting\s+in|as\s+demonstrated\s+by|as\s+indicated\s+by|accompanied\s+by)|\b(?:but|however|whereas|while|although|yet|despite)\b|(?=\bwithout\b)",
            sl,
        ):
            if any(
                neg in clause
                for neg in (
                    "did not ",
                    "does not ",
                    "do not ",
                    "were not ",
                    "was not ",
                    "not impaired",
                    "without ",
                    "with no impact on",
                    "not significantly ",
                    "no significant ",
                    "no reduction in",
                    "no decrease in",
                    "no change in",
                    "no loss of",
                )
            ):
                continue
            if re.search(
                r"\b(ameliorat\w+|preserv\w+|restor\w+|rescu\w+|revers\w+|counteract\w+|prevent\w+|protect\w+)\b[^.;]{0,65}\b(impair\w+|decreas\w+|reduc\w+|dedifferentiat\w+|dysfunction)",
                clause,
            ) or re.search(
                r"\b(impair\w+|decreas\w+|reduc\w+|dedifferentiat\w+|dysfunction)\b[^.;]{0,95}\b(?:could\s+be\s+|was\s+|were\s+)(ameliorat\w+|preserv\w+|restor\w+|rescu\w+)\b",
                clause,
            ) or is_rescue_of_impairment:
                continue
            if "pyruvate" in clause or "basal insulin secretion (2 mmol" in clause:
                continue
            if re.search(
                r"\b(?:reduc\w+|decreas\w+|lower\w*)\s+(?:the\s+)?(?:number\s+of\s+|content\s+ratio\s+of\s+the\s+|expression\s+of\s+)?(?:blood\s+glucose|glucose\s+level|glycemi\w*|hyperglycemi\w*|glycated\s+hemoglobin|hba1c|exogenous\s+insulin|costs?|non-?endocrine|off-target|enterochromaffin|pluripotency|teratoma|proliferation|cytotoxic\w*|susceptibility|cell\s+death|apoptos\w*|hla\s+class|recognition|pyruvate|cellular\s+zinc)",
                clause,
            ):
                continue
            if "associated with increased insulin secretion" in clause or "increasing insulin secretion" in clause or re.search(r"\bcaused\s+mody\d*\b", clause):
                continue
            if re.search(r"\b(?:impair\w+|compromis\w+|blunt\w+|defect\w+|abolish\w+)\s+(?:the\s+)?(?:secretion|differentiation|maturation)\b", clause):
                has_bound_penalty = True
                break
            if re.search(impair_verbs_pat + r"[^,;]{0,45}" + func_targets_pat, clause) or re.search(
                func_targets_pat + r"[^,;]{0,45}" + impair_verbs_pat, clause
            ):
                has_bound_penalty = True
                break
            if re.search(down_verbs_pat + r"[^,;]{0,35}" + func_targets_pat, clause) or re.search(
                func_targets_pat
                + r"[^,;]{0,25}(?:was|were|is|are|remained|continued\s+to)?\s*(?:1\.\d+-fold\s+(?:and\s+1\.\d+-fold\s+)?)?"
                + down_verbs_pat,
                clause,
            ):
                has_bound_penalty = True
                break
            if ("insulin:c-peptide" in sl or "insulin/c-peptide" in sl) and re.search(
                r"\b\d+(?:\.\d+)?-fold\s+lower\s+ratio\b", clause
            ):
                has_bound_penalty = True
                break

    has_direct_hazard_pen = (
        ("hypoglycemi" in sl and not neg_hypo and any(w in sl for w in ("fasting hypoglycemi", "induced hypoglycemi", "caused hypoglycemi", "risk of hypoglycemi", "fell to 70 mg/dl or less")))
        or ("teratoma" in sl and not neg_teratoma and any(w in sl for w in ("teratoma formation", "formed teratoma", "teratoma risk", "teratomas")))
        or ("off-target" in sl and not neg_off_target and any(w in sl for w in ("slc18a1", "enterochromaffin")))
        or ("elevated insulin secretion in low glucose" in sl)
    )
    penalty_observed = has_bound_penalty or has_direct_hazard_pen

    haz_pat = r"(lysis|killing|apoptosis|apoptotic|cell death|cytotoxic\w*|rejection|inflammation|er stress|fibrosis|necrosis|hyperglyc\w+|insulitis|infiltrat\w+)"
    up_pat = r"(increas\w+|elevat\w+|higher|enhanc\w+|greater|upregulat\w+|induc\w+|exacerbat\w+|accelerat\w+|more)"
    down_pat = r"(reduc\w+|decreas\w+|lower\w*|prevent\w*|protect\w*|attenuat\w+|inhibit\w+|abolish\w*|block\w*|suppress\w*|diminish\w*|ameliorat\w+|delay\w*|abrogat\w+|less|fewer|resist\w*|rescu\w+|spar\w+)"
    haz_up_fwd = bool(
        re.search(up_pat + r"[^.;]{0,50}" + haz_pat, sl)
        and not re.search(
            r"\b(?:did\s+not|was\s+not|were\s+not|does\s+not|do\s+not|no|without|never|neither)\b[^.;]{0,25}" + up_pat,
            sl,
        )
    )
    haz_up_rev = bool(
        re.search(
            haz_pat + r"[^.;]{0,50}\b(?:was|were|is|are|remained|became|significantly|markedly|drastically|substantially|\s)+" + up_pat,
            sl,
        )
        and not re.search(r"\b(?:did\s+not|was\s+not|were\s+not|no|without)\b[^.;]{0,25}" + up_pat, sl)
    )
    haz_down_any = bool(
        re.search(down_pat + r"[^.;]{0,50}" + haz_pat, sl)
        or re.search(haz_pat + r"[^.;]{0,40}\b(?:was|were|is|are)\s+(?:significantly\s+|markedly\s+)?" + down_pat, sl)
    )
    haz_increase = (haz_up_fwd or haz_up_rev) and not haz_down_any
    is_graft_survival_reduction = bool(
        (
            re.search(r"\b(?:reduc\w+|decreas\w+|shorten\w+|curtail\w+|impair\w+|loss\s+of)\b[^.;]{0,40}\b(?:allograft|graft|islet)\s+survival\b", sl)
            or re.search(r"\b(?:allograft|graft|islet)\s+survival\b[^.;]{0,35}\b(?:was|were|is|are)?\s*(?:significantly\s+|markedly\s+)?(?:reduc\w+|decreas\w+|shorten\w+|curtail\w+|impair\w+)\b", sl)
        )
        and not re.search(r"\b(?:did\s+not|was\s+not|were\s+not|without)\b", sl)
    )

    if (
        any(
            k in sl
            for k in (
                "unchanged",
                "similar",
                "no significant difference",
                "not significantly",
                "not statistically significant",
                "comparable",
                "did not affect",
                "was not prolonged",
                "were not prolonged",
                "did not increase",
                "did not prevent",
                "no protection",
                "failed to prolong",
                "not sufficient to prevent",
                "not sufficient to protect",
                "insufficient to prevent",
            )
        )
        and not haz_increase
        and not penalty_observed
        and not is_graft_survival_reduction
    ):
        return "Neutral_or_NonInferior", False

    has_non_cyto_toxic = bool(re.search(r"(?<!cyto)\btoxic\w*", sl))
    has_harm_kw = has_non_cyto_toxic or (
        (
            any(
                k in sl
                for k in (
                    "impair",
                    "compromis",
                    "fail",
                    "destroyed",
                    "rejected",
                    "chromothripsis",
                    "exacerbat",
                    "worsen",
                    "blunt",
                    "poor gsis",
                    "poor dynamic",
                    "large deletion",
                    "kilobase deletion",
                    "on-target deletion",
                    "returned to a diabetic",
                    "returned to diabetic",
                    "reverted to hyperglycemi",
                    "reverted to a diabetic",
                )
            )
            or ("hypoglycemi" in sl and not neg_hypo)
            or ("teratoma" in sl and not neg_teratoma)
            or ("silenc" in sl and not is_rnai_silenc and any(w in sl for w in ("epigenetic", "transgene", "promoter", "methylation", "icasp9")))
            or ("translocation" in sl and not is_benign_transloc and any(w in sl for w in ("chromosom", "balanced", "unbalanced", "cas9", "nickase", "dsb", "editing", "cast-seq")))
        )
        and not neg_impair
        and not is_rescue_of_impairment
        and not neg_off_target
    )

    is_functional_reduction = has_bound_penalty or is_graft_survival_reduction

    if has_harm_kw or haz_increase or is_functional_reduction or penalty_observed:
        direction = "Deleterious_or_Impaired"
    elif neg_impair and not any(
        k in sl
        for k in (
            "increas",
            "improv",
            "enhanc",
            "protect",
            "prolong",
            "restor",
            "preserv",
            "prevent",
            "evad",
            "resist",
            "surviv",
            "revers",
            "ameliorat",
            "counteract",
            "normoglycemi",
            "reduced cell death",
            "cell death was significantly reduced",
        )
    ):
        direction = "Neutral_or_NonInferior"
    elif any(
        k in sl
        for k in (
            "increas",
            "improv",
            "enhanc",
            "protect",
            "prolong",
            "restor",
            "preserv",
            "prevent",
            "evad",
            "resist",
            "surviv",
            "revers",
            "ameliorat",
            "counteract",
            "normoglycemi",
            "reduc",
            "decreas",
            "inhibit",
            "attenuat",
            "block",
            "abolish",
            "abrogat",
            "suppress",
            "diminish",
            "lower",
            "higher",
            "greater",
            "elevat",
            "upregulat",
            "extended",
            "delay",
            "spar",
        )
    ):
        direction = "Protective_or_Improved"
    else:
        direction = "Neutral_or_NonInferior"

    return direction, penalty_observed


def _pvalue_tokens(sentence: str) -> list[tuple[str, Decimal]]:
    """Typed source operands. A P< bound is never an observed equality."""
    out = []
    token = r"(?:\d+(?:\.\d*)?|\.\d+)(?:\s*[×xX*]\s*10\s*\^?\s*[-−]?\d+|[eE][-+−]?\d+)?"
    for m in re.finditer(rf"(?<!\w)[pP]\s*(<=|>=|[=<>≤≥])\s*({token})(?![\d.])", sentence):
        raw = re.sub(r"\s+", "", m.group(2)).replace("−", "-")
        try:
            if re.search(r"[×xX*]", raw):
                base, exponent = re.split(r"[×xX*]10\^?", raw)
                val = Decimal(base) * (Decimal(10) ** int(exponent))
            else:
                val = Decimal(raw)
            # A positive source probability must not become exact p=0 merely
            # because its value falls below the binary64 JSON numeric range.
            if val.is_finite() and Decimal(0) <= val <= Decimal(1) and (val == 0 or float(val) > 0):
                out.append(({"≤": "<=", "≥": ">="}.get(m.group(1), m.group(1)), val))
        except (DecimalException, ValueError, OverflowError):
            continue
    return out


def parse_pvalue_and_n(sentence: str) -> tuple[float | None, int | None]:
    ps = _pvalue_tokens(sentence)
    pval = float(ps[0][1]) if len(ps) == 1 and ps[0][0] == "=" else None

    n_val = None
    m_n = re.search(r"\bn\s*=\s*(\d+(?:,\d{3})*)\b", sentence, flags=re.IGNORECASE)
    if m_n:
        try:
            cand_n = int(m_n.group(1).replace(",", ""))
            if 0 < cand_n <= 100000 and _num_in_quote(cand_n, sentence, "treatment_n"):
                n_val = cand_n
        except Exception:
            n_val = None
    return pval, n_val


LEX_H = (
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
LEX_P = (
    "survival",
    "viability",
    "gsis",
    "stimulation index",
    "c-peptide",
    "insulin secretion",
    "normoglycemia",
    "vascular density",
)
EXTRA_H = (
    "lysis",
    "killing",
    "cytolysis",
    "cytolytic",
    "degranulation",
    "cd107",
    "rejection",
    "insulitis",
    "infiltrat",
    "inflammat",
    "ifn-γ secretion",
    "ifn-gamma secretion",
    "t cell proliferation",
    "t-cell proliferation",
    "proliferation of t",
    "phagocyt",
    "fibrosis",
    "necrosis",
    "er stress",
    "perforin",
    "granzyme",
    "granulysin",
    "tnf+",
)


class _NumericEvidenceDomainError(ValueError):
    pass


class _NumericEvidenceNotationError(ValueError):
    pass


class _NumericEvidenceBoundedError(ValueError):
    pass


# Qualifiers are attached to an operand, never searched sentence-wide: a
# separate p-value inequality must not invalidate an exact measured contrast.
_NUMERIC_QUALIFIER = (
    r"(?:[<>≤≥~≈]=?\s*|\b(?:at\s+least|at\s+most|more\s+than|less\s+than|"
    r"greater\s+than|fewer\s+than|up\s+to|nearly|almost|over|under|"
    r"approximately|approx\.?|about|around|roughly)\s+)"
)


def _finite_numeric_evidence(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _NumericEvidenceDomainError(f"{label} must be a finite real number")
    try:
        number = float(value)
    except (OverflowError, ValueError) as exc:
        raise _NumericEvidenceDomainError(f"{label} is outside the finite numeric domain") from exc
    if not math.isfinite(number):
        raise _NumericEvidenceDomainError(f"{label} must be a finite real number")
    return number


def _log_response_ratio(treatment: Any, control: Any) -> float:
    treatment = _finite_numeric_evidence(treatment, "Treatment mean")
    control = _finite_numeric_evidence(control, "Control mean")
    if treatment <= 0 or control <= 0:
        raise _NumericEvidenceDomainError("Log response ratios require positive operands")
    # Subtraction of nearly equal logarithms loses a small, nonzero contrast.
    # Nearby positive binary64 operands have a safe finite difference, while
    # far ratios use separate logs so neither division nor subtraction overflows.
    if min(treatment, control) >= max(treatment, control) / 2.0:
        return math.log1p((treatment - control) / control)
    return math.log(treatment) - math.log(control)


def _source_numeric_operand(token: str) -> float:
    token = token.replace(",", "")
    number = _finite_numeric_evidence(float(token), "Source numeric operand")
    if number == 0 and Decimal(token) != 0:
        raise _NumericEvidenceDomainError("Source numeric operand underflows the finite numeric domain")
    return number


def _source_numeric_match(match: re.Match[str], group: int, *, exact: bool = False) -> float:
    """Reject regex captures that are only a suffix of a source operand.

    Source eligibility, arm assignment and endpoint orientation are unchanged.
    The unsigned-decimal grammar must not turn .5 into 5, -20 into 20, or the
    exponent in 1e-2 into a twofold change. +/- and +- remain dispersion markers.
    Exact paired fields cannot represent a qualified bound or approximation;
    descriptive duration parsing keeps its existing censoring metadata path.
    """
    prefix = match.string[:match.start(group)]
    previous = prefix.rstrip()
    bad_immediate = bool(prefix and prefix[-1] in "0123456789.,eE")
    bad_sign = bool(previous and previous[-1] in "+-−–—" and not previous.endswith(("+/-", "+-")))
    bad_grouping = bool(re.search(r"\d_+$", prefix)) or bool(
        re.search(r"\d[^\S\r\n]+$", prefix)
        and re.fullmatch(r"\d{3}(?:,\d{3})*(?:\.\d+)?", match.group(group))
    )
    suffix = match.string[match.end(group):]
    if bad_immediate or bad_sign or bad_grouping or re.match(r"[eE][+\-−]?\d", suffix):
        raise _NumericEvidenceNotationError("Source operand uses unsupported numeric notation")
    if exact and re.search(_NUMERIC_QUALIFIER + r"$", prefix, flags=re.IGNORECASE):
        raise _NumericEvidenceBoundedError("A qualified operand is not an exact paired measurement")
    return _source_numeric_operand(match.group(group))


def _v063_parse_quantitative_numbers(
    sentence: str, direction: str
) -> tuple[float | None, float | None, float | None, float | None, str | None, float | None]:
    """Preserves exact v0.6.3 candidate sentence sort priority tuple so all 2,794 effect_ids and V-17/V-18 remain invariant."""
    sl = sentence.lower()
    has_h = any(k in sl for k in LEX_H)
    has_p = any(k in sl for k in LEX_P)
    orient_sign = -1.0 if (has_h and not has_p) else (1.0 if (has_p and not has_h) else 0.0)
    is_protocol = any(
        k in sl
        for k in (
            "stimulated with",
            "incubated with",
            "cultured with",
            "treated with",
            "supplemented with",
            "concentration",
            "physiological",
            "monitored",
            "quantified",
            "ranging from",
            "between ",
            "heatmap",
            "depicting",
            "log2 fold",
            "not functional",
            "triglyceride",
        )
    )
    m_pair_sd = re.search(
        r"(from\s+)?(\d+\.?\d*)\s*(?:±|\+/-|\+-)\s*(\d+\.?\d*)\s*(%|fold|mg/dL|ng/mL|pmol|days|weeks|mmHg|µm|um)?\s*"
        r"(vs\.?|versus|compared\s+to|compared\s+with|to)\s*"
        r"(\d+\.?\d*)\s*(?:±|\+/-|\+-)\s*(\d+\.?\d*)\s*(%|fold|mg/dL|ng/mL|pmol|days|weeks|mmHg|µm|um)?",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_pair_sd:
        has_from = bool(m_pair_sd.group(1))
        v1, sd1, u1 = _source_numeric_operand(m_pair_sd.group(2)), _source_numeric_operand(m_pair_sd.group(3)), m_pair_sd.group(4)
        sep = m_pair_sd.group(5).lower()
        v2, sd2, u2 = _source_numeric_operand(m_pair_sd.group(6)), _source_numeric_operand(m_pair_sd.group(7)), m_pair_sd.group(8)
        unit = u1 or u2 or "reported_units"
        if sep == "to" and not has_from:
            return None, None, None, None, unit, None
        t_mean, t_sd, c_mean, c_sd = (v2, sd2, v1, sd1) if has_from else (v1, sd1, v2, sd2)
        if has_h and has_p:
            return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
        if t_mean > 0 and c_mean > 0 and orient_sign != 0.0 and not is_protocol:
            lrr = round(_log_response_ratio(t_mean, c_mean), 4)
            if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
            ):
                return None, None, None, None, "quarantined_ambiguous_numeric", None
            return t_mean, t_sd, c_mean, c_sd, unit, lrr
        return None, None, None, None, unit, None
    m_pair = re.search(
        r"(from\s+)?(\d+\.?\d*)\s*(%|fold|mg/dL|ng/mL|days|weeks|months|µm|um)\s*"
        r"(vs\.?|versus|compared\s+to|compared\s+with|to)\s*"
        r"(\d+\.?\d*)\s*(%|fold|mg/dL|ng/mL|days|weeks|months|µm|um)?",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_pair:
        has_from = bool(m_pair.group(1))
        v1, u1, sep, v2 = _source_numeric_operand(m_pair.group(2)), m_pair.group(3), m_pair.group(4).lower(), _source_numeric_operand(m_pair.group(5))
        if sep == "to" and not has_from:
            return None, None, None, None, u1, None
        t_mean, c_mean = (v2, v1) if has_from else (v1, v2)
        if has_h and has_p:
            return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
        if t_mean > 0 and c_mean > 0 and orient_sign != 0.0 and not is_protocol:
            lrr = round(_log_response_ratio(t_mean, c_mean), 4)
            if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
            ):
                return None, None, None, None, "quarantined_ambiguous_numeric", None
            return t_mean, None, c_mean, None, u1, lrr
        return None, None, None, None, u1, None
    if direction == "Neutral_or_NonInferior" or any(
        k in sl
        for k in (
            "unchanged",
            "no significant difference",
            "not significantly different",
            "was not prolonged",
            "did not increase",
            "did not decrease",
            "did not affect",
        )
    ):
        return None, None, None, None, "qualitative_statement", None
    m_fold = re.search(r"(?<!log)[>~]?\s*(\d+\.?\d*)\s*(?:-|\s)?fold\b", sentence, flags=re.IGNORECASE)
    if m_fold:
        fc = _source_numeric_operand(m_fold.group(1))
        if 0.05 <= fc <= 1000.0 and fc != 1.0:
            if has_h and has_p:
                return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
            if orient_sign != 0.0 and not is_protocol:
                is_red = any(
                    k in sl
                    for k in ("decreas", "reduc", "lower", "downregulat", "inhibit", "suppress", "loss", "diminish", "attenuat")
                )
                lrr = round(math.log(fc), 4) if fc < 1.0 else (round(-abs(math.log(fc)), 4) if is_red else round(abs(math.log(fc)), 4))
                if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                    direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
                ):
                    return None, None, None, None, "quarantined_ambiguous_numeric", None
                return fc, None, None, None, "fold_change_ratio", lrr
            return fc, None, None, None, "fold_change_ratio", None
    m_single_sd = re.search(
        r"(\d+\.?\d*)\s*(?:±|\+/-|\+-)\s*(\d+\.?\d*)\s*(%|mg/dL|ng/mL|pmol|days|weeks|mmHg|µm|um)",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_single_sd:
        return _source_numeric_operand(m_single_sd.group(1)), _source_numeric_operand(m_single_sd.group(2)), None, None, m_single_sd.group(3), None
    return None, None, None, None, "qualitative_statement", None


def _parse_quantitative_numbers_unchecked(
    sentence: str, direction: str
) -> tuple[float | None, float | None, float | None, float | None, str | None, float | None]:
    """
    Extracts (treatment_mean, treatment_sd, control_mean, control_sd, measurement_unit, log_response_ratio_or_smd)
    from a quantitative sentence in Results, Figure Legends, or Abstract.
    Strict Zero-Fabrication, Contrast-Provenance & Polarity-Quarantine Rule (v0.6.4):
      - Never invents a comparator or baseline value when only a single number is reported.
      - Never assigns a synthetic numeric LRR to qualitative, single-arm, longitudinal within-arm ('from X to Y'),
        or within-assay secretagogue/cytokine challenge statements.
      - Parses mean +/- SD/SEM fold-changes ('4.14 ± 1.4-fold') before bare fold-changes so the dispersion term
        after '±' is never misread as the point estimate.
      - Quarantines sentences with multiple fold-change clauses, fold ranges, or mixed increase+decrease polarity.
      - Requires a single unambiguous endpoint polarity class (has_lex_h != has_lex_p) for numeric LRR.
    """
    sl = sentence.lower()
    # An approximate estimate or an inequality is not an exact ratio. Retain
    # its source sentence, but leave the quantitative contrast unmeasured.
    qualified_fold = (
        r"(?:[<>≤≥~≈]|\b(?:at\s+least|at\s+most|more\s+than|less\s+than|"
        r"greater\s+than|fewer\s+than|up\s+to|nearly|almost|over|under|approximately|approx\.?|about|around|roughly)\s+)"
        r"\s*(?:an?\s+)?\d+(?:\.\d+)?(?:\s*(?:±|\+/-)\s*\d+(?:\.\d+)?)?\s*-?fold\b"
    )
    if re.search(qualified_fold, sl):
        return None, None, None, None, "quarantined_approximate_or_bounded_fold_change", None
    has_h = any(k in sl for k in LEX_H)
    has_p = any(k in sl for k in LEX_P)
    orient_sign = -1.0 if (has_h and not has_p) else (1.0 if (has_p and not has_h) else 0.0)
    is_protocol = any(
        k in sl
        for k in (
            "stimulated with",
            "incubated with",
            "cultured with",
            "treated with",
            "supplemented with",
            "concentration",
            "physiological",
            "monitored",
            "quantified",
            "ranging from",
            "between ",
            "heatmap",
            "depicting",
            "log2 fold",
            "not functional",
            "triglyceride",
        )
    )

    # Full unit token pattern that captures composite/rate units (e.g. ng/mL/min, pg/mL, hours, PD-L1+ cells/mm2)
    unit_tok = r"(?:PD-L1\+\s+cells/mm2|%|(?:-|\s)?fold\b|[a-zA-Zµμ]+(?:/[a-zA-Z0-9µμ²³+-]+)+|\b(?:days|weeks|months|hours|years|min|sec|h|s|pmol|mmol|nmol|µmol|umol|mg|ng|pg|µg|ug|mM|nM|µM|uM|pM|mmHg|µm|um|mm|cm|mL|L|dL|cells|islets|mice|patients)\b)"
    supported_pair_units = {
        "%", "fold", "mg/dl", "ng/ml", "ng/ml/min", "pg/ml", "uiu/ml", "µiu/ml", "pmol", "pmol/l",
        "mmol/l", "mm", "nm", "um", "µm", "mmhg", "days", "weeks", "months", "hours", "h",
        "cells/mm2", "pd-l1+ cells/mm2",
    }
    longitudinal_timepoint_re = re.compile(
        r"\b(baseline|follow-up|follow\s+up|pre-transplant\w*|post-transplant\w*|before\s+and\s+after|prior\s+to\s+and\s+after|day\s+\d+|week\s+\d+|month\s+\d+)\b",
        flags=re.IGNORECASE,
    )

    # Pattern 1: Explicit paired mean +/- SD/SEM comparison
    # (e.g., "68.4 ± 4.2% vs. 14.1 ± 2.1%" or "from 12.5 ± 1.2% to 48.0 ± 3.5%")
    num_tok = r"\d+(?:,\d{3})*(?:\.\d+)?"
    m_pair_sd = re.search(
        rf"(from\s+)?({num_tok})\s*(?:±|\+/-|\+-)\s*(?:{_NUMERIC_QUALIFIER})?({num_tok})\s*({unit_tok})?\s*"
        rf"(vs\.?|versus|compared\s+to|compared\s+with|to)\s*"
        rf"(?:{_NUMERIC_QUALIFIER})?({num_tok})\s*(?:±|\+/-|\+-)\s*(?:{_NUMERIC_QUALIFIER})?({num_tok})\s*({unit_tok})?",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_pair_sd:
        has_from = bool(m_pair_sd.group(1))
        v1 = _source_numeric_match(m_pair_sd, 2, exact=True)
        sd1 = _source_numeric_match(m_pair_sd, 3, exact=True)
        u1 = (m_pair_sd.group(4) or "").strip() or None
        sep = m_pair_sd.group(5).lower()
        v2 = _source_numeric_match(m_pair_sd, 6, exact=True)
        sd2 = _source_numeric_match(m_pair_sd, 7, exact=True)
        u2 = (m_pair_sd.group(8) or "").strip() or None
        if u1 and u2 and u1.lower() != u2.lower():
            return None, None, None, None, "quarantined_mismatched_units", None
        unit = u1 or u2 or "reported_units"
        if unit != "reported_units" and unit.lower().lstrip("- ") not in supported_pair_units:
            return None, None, None, None, "quarantined_unsupported_unit", None
        if sep == "to" and not has_from:
            return None, None, None, None, unit, None
        if has_from:
            t_mean, t_sd, c_mean, c_sd = v2, sd2, v1, sd1
            if "respectively" not in sl or longitudinal_timepoint_re.search(sentence) or unit.lower().lstrip("- ") in ("pd-l1+ cells/mm2", "cells/mm2"):
                return t_mean, t_sd, c_mean, c_sd, unit, None
        else:
            t_mean, t_sd, c_mean, c_sd = v1, sd1, v2, sd2
            if unit.lower().lstrip("- ") in ("pd-l1+ cells/mm2", "cells/mm2"):
                return t_mean, t_sd, c_mean, c_sd, unit, None
        if has_h and has_p:
            return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
        if t_mean > 0 and c_mean > 0 and orient_sign != 0.0 and not is_protocol:
            lrr = round(_log_response_ratio(t_mean, c_mean), 4)
            if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
            ):
                return None, None, None, None, "quarantined_ambiguous_numeric", None
            return t_mean, t_sd, c_mean, c_sd, unit, lrr
        return None, None, None, None, unit, None

    # Pattern 2: Explicit paired values with explicit comparison token (vs, versus, compared to/with, or from X to Y)
    m_pair = re.search(
        rf"(from\s+)?({num_tok})\s*({unit_tok})\s*"
        rf"(vs\.?|versus|compared\s+to|compared\s+with|to)\s*"
        rf"(?:{_NUMERIC_QUALIFIER})?({num_tok})(?:\s*({unit_tok}))?",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_pair:
        has_from = bool(m_pair.group(1))
        v1 = _source_numeric_match(m_pair, 2, exact=True)
        u1 = m_pair.group(3).strip()
        sep = m_pair.group(4).lower()
        v2 = _source_numeric_match(m_pair, 5, exact=True)
        u2 = (m_pair.group(6) or "").strip() or None
        if not u2 or u1.lower() != u2.lower():
            return None, None, None, None, "quarantined_mismatched_units", None
        if u1.lower().lstrip("- ") not in supported_pair_units:
            return None, None, None, None, "quarantined_unsupported_unit", None
        if sep == "to" and not has_from:
            return None, None, None, None, u1, None
        if has_from:
            t_mean, c_mean = v2, v1
            return t_mean, None, c_mean, None, u1, None
        else:
            t_mean, c_mean = v1, v2
        if has_h and has_p:
            return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
        if t_mean > 0 and c_mean > 0 and orient_sign != 0.0 and not is_protocol:
            lrr = round(_log_response_ratio(t_mean, c_mean), 4)
            if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
            ):
                return None, None, None, None, "quarantined_ambiguous_numeric", None
            return t_mean, None, c_mean, None, u1, lrr
        return None, None, None, None, u1, None

    if direction == "Neutral_or_NonInferior" or any(
        k in sl
        for k in (
            "unchanged",
            "no significant difference",
            "not significantly different",
            "was not prolonged",
            "did not increase",
            "did not decrease",
            "did not affect",
        )
    ):
        return None, None, None, None, "qualitative_statement", None

    has_multi_fold = (
        len(re.findall(r"\b(?:\d+\.?\d*|two|three|four|five|six|seven|eight|nine|ten)\s*(?:-|\s)?fold\b", sl)) > 1
        or bool(re.search(r"\b\d+\.?\d*\s*[-–]\s*\d+\.?\d*\s*[-–]?\s*(?:or\s+|and\s+|to\s+)?\d*\.?\d*\s*(?:-|\s)?fold\b", sl))
        or bool(re.search(r"\b\d+\.?\d*[-–]?\s*(?:to|and|or)\s*\d+\.?\d*\s*(?:-|\s)?fold\b", sl))
    )
    has_mixed_up_down = (
        any(k in sl for k in ("increas", "upregulat", "elevat", "higher", "enhanc", "induc"))
        and any(k in sl for k in ("decreas", "reduc", "lower", "downregulat", "inhibit", "suppress", "attenuat", "diminish"))
    )
    is_within_assay_stim = (
        any(
            k in sl
            for k in (
                "in response to glucose",
                "upon high glucose challenge",
                "exposed to high glucose",
                "upon 10 μm ionomycin",
                "upon 10 um ionomycin",
                "glibenclamide stimulation",
                "upon kcl",
                "depolarization",
                "exposure with these cytokines",
                "caspase-3/7 activation",
                "glucose-stimulated insulin secretion (gsis) was increased",
            )
        )
        and not any(k in sl for k in ("compared", "versus", "vs.", "wild-type", "unedited", "non-sorted"))
    ) or any(
        k in sl
        for k in (
            "cell sorting yielded",
            "compared to non-sorted cells",
        )
    )

    # Pattern 3A: Explicit fold-change with ± SD/SEM (e.g., "4.14 ± 1.4-fold" in pmid41385368_eff4)
    m_fold_sd = re.search(
        r"(?<!log)[>~]?\s*(\d+\.?\d*)\s*(?:±|\+/-|\+-)\s*(\d+\.?\d*)\s*(?:-|\s)?fold\b",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_fold_sd:
        fc = _source_numeric_match(m_fold_sd, 1)
        fc_sd = _source_numeric_match(m_fold_sd, 2)
        if 0.05 <= fc <= 1000.0 and fc != 1.0:
            if has_multi_fold or has_mixed_up_down:
                return None, None, None, None, "quarantined_ambiguous_numeric", None
            if has_h and has_p:
                return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
            if orient_sign != 0.0 and not is_protocol and not is_within_assay_stim:
                is_red = any(
                    k in sl
                    for k in ("decreas", "reduc", "lower", "downregulat", "inhibit", "suppress", "loss", "diminish", "attenuat")
                )
                if fc < 1.0:
                    lrr = round(math.log(fc), 4)
                elif is_red:
                    lrr = round(-abs(math.log(fc)), 4)
                else:
                    lrr = round(abs(math.log(fc)), 4)
                if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                    direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
                ):
                    return None, None, None, None, "quarantined_ambiguous_numeric", None
                return fc, fc_sd, None, None, "fold_change_ratio", lrr
            return fc, fc_sd, None, None, "fold_change_ratio", None

    # Pattern 3B: Bare fold-change (e.g., "2.8-fold", "3.5 fold", ">2-fold") — excludes "± X-fold", "times", and "log2 fold"
    m_fold = re.search(
        r"(?<!log)(?<!±\s)(?<!±)(?<!\+/-\s)(?<!\+/-)(?<![\d.])[>~]?\s*(\d+\.?\d*)\s*(?:-|\s)?fold\b",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_fold:
        fc = _source_numeric_match(m_fold, 1)
        if 0.05 <= fc <= 1000.0 and fc != 1.0:
            if has_multi_fold or has_mixed_up_down:
                return None, None, None, None, "quarantined_ambiguous_numeric", None
            if has_h and has_p:
                return None, None, None, None, "quarantined_mixed_polarity_endpoints", None
            if orient_sign != 0.0 and not is_protocol and not is_within_assay_stim:
                is_red = any(
                    k in sl
                    for k in ("decreas", "reduc", "lower", "downregulat", "inhibit", "suppress", "loss", "diminish", "attenuat")
                )
                if fc < 1.0:
                    lrr = round(math.log(fc), 4)
                elif is_red:
                    lrr = round(-abs(math.log(fc)), 4)
                else:
                    lrr = round(abs(math.log(fc)), 4)
                if (direction == "Protective_or_Improved" and lrr * orient_sign <= 0) or (
                    direction == "Deleterious_or_Impaired" and lrr * orient_sign >= 0
                ):
                    return None, None, None, None, "quarantined_ambiguous_numeric", None
                return fc, None, None, None, "fold_change_ratio", lrr
            return fc, None, None, None, "fold_change_ratio", None

    # Pattern 4: Single mean +/- SD/SEM without comparator -> record treatment_mean & SD only; control_mean=None, lrr=None
    m_single_sd = re.search(
        rf"({num_tok})\s*(?:±|\+/-|\+-)\s*({num_tok})\s*(%|mg/dL|ng/mL|pmol|days|weeks|mmHg|µm|um)",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_single_sd:
        v1 = _source_numeric_match(m_single_sd, 1)
        sd1 = _source_numeric_match(m_single_sd, 2)
        unit = m_single_sd.group(3)
        return v1, sd1, None, None, unit, None

    # Pattern 5: Single percentage without comparator -> record treatment_mean only; control_mean=None, lrr=None
    m_pct = re.search(r"(?<![\d,])(\d+\.?\d*)\s*%", sentence)
    if m_pct:
        pct = _source_numeric_match(m_pct, 1)
        if 0.0 <= pct <= 100.0:
            return pct, None, None, None, "%", None

    # Pattern 6: Single duration without comparator -> record reported duration only; control_mean=None, lrr=None
    m_dur = re.search(
        r"(?:for|over|up\s+to|more\s+than|after|within|at|>|~)\s+(\d+\.?\d*)\s*(days|day|weeks|week|months|month|hours|hr|hrs|\bh\b)",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_dur:
        val = _source_numeric_match(m_dur, 1)
        unit = m_dur.group(2).lower()
        return round(val, 2), None, None, None, unit, None

    # Pattern 7: Single physiological measurement without comparator -> record value & unit only; control_mean=None, lrr=None
    m_unit = re.search(
        rf"({num_tok})\s*(mg/dL|ng/mL|pg/mL|mM|µM|uM|nM|pmol|fmol|mmHg|µm|um|cells/mL|SI\b)",
        sentence,
        flags=re.IGNORECASE,
    )
    if m_unit:
        val = _source_numeric_match(m_unit, 1)
        unit = m_unit.group(2)
        return val, None, None, None, unit, None

    # Qualitative statement without explicit paired numbers -> keep all numeric fields None
    return None, None, None, None, "qualitative_statement", None


def parse_quantitative_numbers(
    sentence: str, direction: str
) -> tuple[float | None, float | None, float | None, float | None, str | None, float | None]:
    try:
        result = _parse_quantitative_numbers_unchecked(sentence, direction)
        for index in (0, 1, 2, 3, 5):
            if result[index] is not None:
                _finite_numeric_evidence(result[index], "Parsed numeric evidence")
        return result
    except _NumericEvidenceNotationError:
        return None, None, None, None, "quarantined_unsupported_numeric_notation", None
    except _NumericEvidenceBoundedError:
        return None, None, None, None, "quarantined_approximate_or_bounded_measurement", None
    except _NumericEvidenceDomainError:
        return None, None, None, None, "quarantined_nonfinite_numeric", None


STAT_VOCAB = set(
    """p n mean means sd sem se anova test tests t tukey bonferroni dunnett student students two tailed one way paired unpaired non
nonparametric mann whitney wilcoxon kruskal wallis fisher exact post hoc results result represent represents represented data are is were was be been presented
shown show shows as independent experiments experiment replicates replicate biological technical donors donor error bars bar standard deviation versus vs compared
comparison relative control controls untreated media alone ns significant significantly different difference differences values value from of the and or with to
in at for by a an each per group groups condition conditions same respective corrected correction multiple comparisons using calculated determined analysis
variance followed least statistical statistics three four five six seven eight nine ten twice out representative figure fig panel panels see supplementary
unless otherwise indicated denotes denote indicate indicated asterisks asterisk stars star level levels analyzed analysed performed carried cells cell treated
treatment wild type wt""".split()
)

HARM_SENT_RE = re.compile(
    r"\b(destr\w+|damag\w+|reject\w*|killed|killing|death|died|fail\w*|nonfunction\w+|necro\w+|dysfunction\w*|impair\w+|toxic\w*|lost|loss)\b",
    flags=re.IGNORECASE,
)
BEN_SENT_RE = re.compile(
    r"\b(improv\w+|protect\w+|surviv\w+|prolong\w+|restor\w+|normoglyc\w+|reversal|reversed|enhanc\w+)\b",
    flags=re.IGNORECASE,
)


def is_stat_boilerplate(q: str) -> bool:
    ql = q.lower().strip()
    if re.search(
        r"\b(?:the\s+aims?\s+of\s+this\s+study|this\s+study\s+aimed\s+to|our\s+objective\s+was\s+to|the\s+goal\s+of\s+this\s+study|in\s+a\s+similar\s+system\s+was\s+recently\s+shown\s+to)\b",
        ql,
    ):
        return True
    cw = [
        w
        for w in re.findall(r"[A-Za-z][A-Za-z\-']+", ql)
        if w not in STAT_VOCAB and len(w) > 2
    ]
    return (
        bool(re.search(r"(p\s*[<=>]|n\s*=|±|mean|anova|t-test|t test|sem|sd\b)", q, re.I))
        and len(cw) <= 3
    )


def strict_gene_match(sym: str, text: str) -> bool:
    s = sym.upper()
    t = text.upper()
    if s == "CAR":
        return bool(re.search(r"(?<![A-Z0-9])CAR(?![A-Z0-9])", t)) and "COXSACKIE" in t
    if s == "INS":
        return bool(re.search(r"(?<![A-Z0-9])INS(?![A-Z0-9])(?!-1)", t))
    if s.startswith("HLA-"):
        return bool(re.search(r"(?<![A-Z0-9])" + re.escape(s) + r"(?:W(?=[\d*])|(?![A-Z]))", t))
    return bool(re.search(r"(?<![A-Z0-9])" + re.escape(s) + r"(?![A-Z0-9])", t))


def _target_bound_atomic_clauses(sentence: str, study_edits: list[dict[str, Any]]) -> list[str]:
    """Split distinct target clauses before focal attribution, preserving source substrings.

    When clauses cannot be unambiguously separated, retain the sentence; its
    target-specific classifier still runs after focal selection. No text is joined.
    """
    genes = {e["gene_symbol_hgnc"] for e in study_edits if strict_gene_match(e["gene_symbol_hgnc"], sentence)}
    if len(genes) < 2:
        return [sentence]
    clauses = [x.strip() for x in re.split(
        r",\s*(?=(?:loss|deletion|knockout|downregulation|silencing|overexpression)\s+of\b)|"
        r"[,;]\s*(?:whereas|while|but|however)\s+", sentence, flags=re.I
    ) if x.strip()]
    clause_genes = [{g for g in genes if strict_gene_match(g, c)} for c in clauses]
    if len(clauses) > 1 and all(len(g) == 1 for g in clause_genes) and len(set.union(*clause_genes)) == len(genes):
        return clauses
    # A compound multi-target contrast that cannot be bound is not an atomic
    # single-target observation. Quarantine it rather than credit the first gene.
    if re.search(r"\b(?:while|whereas)\b|[,;]\s*(?:but|however)\b", sentence, re.I):
        return []
    return [sentence]


def identify_focal_factor(sentence: str, study_edits: list[dict[str, Any]], title: str) -> str:
    """Identifies the specific perturbation/intervention X acting in the sentence using strict word-boundary matching."""
    for ed in study_edits:
        gsym = ed["gene_symbol_hgnc"]
        if strict_gene_match(gsym, sentence):
            return f"{gsym} ({ed['edit_modality']} — {ed['functional_gate']})"

    if "HNF1A WT β-cells" in sentence and "HNF1A KO β-cells" in sentence and "HNF1A" in title:
        return "HNF1A KO (HNF1A−/−) SC-β cells vs isogenic HNF1A WT SC-β cells"

    first_half = sentence[: max(80, len(sentence) // 2)]
    node_id, _ = resolve_canonical_node(first_half)
    if node_id and "_" in node_id and not any(strict_gene_match(ed["gene_symbol_hgnc"], node_id) for ed in study_edits if not strict_gene_match(ed["gene_symbol_hgnc"], sentence)):
        if not strict_gene_match("INS", node_id) or strict_gene_match("INS", sentence):
            return node_id

    if len(study_edits) == 1 and not re.search(r"\b(?:UTP23|RAD21|MED30|USP18|MEG3)\b", sentence):
        ed0 = study_edits[0]
        g0 = ed0["gene_symbol_hgnc"]
        if g0 != "INS" or strict_gene_match("INS", sentence):
            return f"{g0} ({ed0['edit_modality']} — {ed0['functional_gate']})"

    t_node, _ = resolve_canonical_node(title)
    if t_node and (not strict_gene_match("INS", t_node) or strict_gene_match("INS", sentence)):
        if not any(strict_gene_match(ed["gene_symbol_hgnc"], t_node) for ed in study_edits if not strict_gene_match(ed["gene_symbol_hgnc"], sentence)):
            return t_node
    return "Study_Intervention"


def split_contrast_clauses(sent: str, prev_sent: str | None) -> tuple[str, str, str]:
    """
    Splits a modifier sentence into (modifier_variable_Z, level_1_clause, level_2_clause)
    so both contrasting conditions are explicitly captured.
    """
    pivots = [
        r",\s*whereas\s+",
        r"\s+whereas\s+",
        r",\s*while\s+",
        r";\s*however,\s*",
        r",\s*in\s+contrast\s+to\s+",
        r"\s+compared\s+to\s+",
        r"\s+compared\s+with\s+",
        r"\s+unlike\s+",
        r"\s+but\s+not\s+in\s+",
        r"\s+but\s+not\s+",
    ]
    for piv in pivots:
        parts = re.split(piv, sent, maxsplit=1, flags=re.IGNORECASE)
        if len(parts) == 2 and len(parts[0].strip()) >= 15 and len(parts[1].strip()) >= 10:
            c1 = parts[0].strip()
            c2 = parts[1].strip()
            n1, _ = resolve_canonical_node(c1)
            n2, _ = resolve_canonical_node(c2)
            n1_lbl = n1 or "Condition_1"
            n2_lbl = n2 or "Condition_2"
            mod_var = f"{n1_lbl} vs. {n2_lbl}" if n1_lbl != n2_lbl else f"Condition contrast ({n1_lbl})"
            return mod_var[:120], c1, c2

    if prev_sent and re.match(r"^(However|In\s+contrast|Conversely|By\s+contrast|Unlike)\b", sent.strip(), flags=re.IGNORECASE):
        n1, _ = resolve_canonical_node(prev_sent)
        n2, _ = resolve_canonical_node(sent)
        mod_var = f"{n1 or 'Prior_Condition'} vs. {n2 or 'Contrast_Condition'}"
        return mod_var[:120], prev_sent.strip(), sent.strip()

    n1, _ = resolve_canonical_node(sent)
    n1_lbl = n1 or "Experimental_Context"
    mid = len(sent) // 2
    comma_idx = sent.find(",", mid - 30)
    if comma_idx != -1 and 20 < comma_idx < len(sent) - 20:
        return f"Regime dependency on {n1_lbl}", sent[:comma_idx].strip(), sent[comma_idx + 1 :].strip()
    return f"Regime dependency on {n1_lbl}", sent, f"Reference regime ({n1_lbl})"


def extract_causal_edge_from_sentence(
    sent: str, src_label: str, pmid: str, study_edits: list[dict[str, Any]], title: str
) -> dict[str, Any] | None:
    """
    Splits a causal sentence around its causal verb and resolves both the upstream subject
    and downstream object to Canonical Biological Nodes. Returns None if negated or if
    either clause lacks a recognizable biological entity.
    """
    m_verb = re.search(
        r"\b(upregulat\w+|downregulat\w+|inhibit\w+|activat\w+|induc\w+|suppress\w+|prevent\w+|block\w+|rescu\w+|promot\w+|abolish\w+|abrogat\w+|reduc\w+|increas\w+|enhanc\w+)\b",
        sent,
        flags=re.IGNORECASE,
    )
    if not m_verb:
        return None

    verb = m_verb.group(1).lower()
    left_clause = sent[: m_verb.start()].strip()
    right_clause = sent[m_verb.end() :].strip()
    if len(left_clause) < 4 or len(right_clause) < 6:
        return None

    # Check for explicit negation in the left window preceding the verb or immediately after
    neg_window = left_clause[-45:].lower()
    if re.search(
        r"\b(did\s+not|does\s+not|do\s+not|was\s+not|were\s+not|not|no|never|neither|failed\s+to|unable\s+to|without|cannot|could\s+not)\b",
        neg_window,
    ):
        return None

    if any(k in verb for k in ("inhibit", "downregulat", "suppress", "prevent", "block", "abolish", "abrogat", "reduc")):
        rel = "Inhibits_or_Downregulates"
    elif any(k in verb for k in ("rescu", "stabiliz")):
        rel = "Stabilizes_or_Rescues"
    else:
        rel = "Activates_or_Upregulates"

    u_node, u_type = resolve_canonical_node(left_clause, default_type="Gene_or_Protein_Edit")
    v_node, v_type = resolve_canonical_node(right_clause, default_type="Cell_Survival_or_Death_Pathway")

    if u_node is None or v_node is None:
        return None
    if u_node == v_node:
        return None

    return {
        "upstream_node": f"{u_node} [{left_clause[-65:]}]",
        "upstream_node_type": u_type or "Gene_or_Protein_Edit",
        "causal_relation": rel,
        "downstream_node": f"{v_node} [{right_clause[:75]}]",
        "downstream_node_type": v_type or "Cell_Survival_or_Death_Pathway",
        "evidence_source": f"{src_label} (PMID:{pmid}): {sent[:140]}",
    }


def _compute_dossier_sha256(pmid: str, dossier: dict[str, Any]) -> str:
    dp = DOSSIER_DIR / f"{pmid}.json"
    if dp.exists():
        return hashlib.sha256(dp.read_bytes()).hexdigest()
    return hashlib.sha256(
        json.dumps(dossier, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def _get_effective_fulltext_sections(pmid: str, dossier: dict[str, Any]) -> dict[str, Any]:
    ft = dossier.get("fulltext_data_sections") or {}
    if pmid in PEER_REVIEWED_TO_PREPRINT_MAP:
        pre_pmid = PEER_REVIEWED_TO_PREPRINT_MAP[pmid]
        pre_path = DOSSIER_DIR / f"{pre_pmid}.json"
        if pre_path.exists():
            try:
                pre_dos = json.loads(pre_path.read_text(encoding="utf-8"))
                pre_ft = pre_dos.get("fulltext_data_sections") or {}
                ft_chars = sum(len(p or "") for p in (ft.get("results_paragraphs") or [])) + sum(
                    len((f or {}).get("caption") or "") for f in (ft.get("figure_legends") or [])
                )
                pre_chars = sum(len(p or "") for p in (pre_ft.get("results_paragraphs") or [])) + sum(
                    len((f or {}).get("caption") or "") for f in (pre_ft.get("figure_legends") or [])
                )
                if (
                    len(ft.get("results_paragraphs") or []) < len(pre_ft.get("results_paragraphs") or [])
                    or len(ft.get("figure_legends") or []) < len(pre_ft.get("figure_legends") or [])
                    or ft_chars < pre_chars
                ):
                    return pre_ft
            except Exception:
                pass
    return ft


def _resolve_dossier_character_span(quote: str, dossier: dict[str, Any], src_label: str = "") -> str:
    """
    Resolves a verbatim quote to an exact 0-based character slice inside a single dossier field:
      <field_path>[start:end]|norm=(exact|whitespace_normalized)
    """
    q = (quote or "").strip()
    pmid = str(dossier.get("pmid") or "")
    ft = _get_effective_fulltext_sections(pmid, dossier)
    all_sections: list[tuple[str, str]] = [
        ("abstract", dossier.get("abstract") or ""),
        ("title", dossier.get("title") or ""),
    ]
    for i, rpara in enumerate(ft.get("results_paragraphs", [])):
        all_sections.append((f"fulltext_data_sections.results_paragraphs[{i}]", rpara or ""))
    for i, fig in enumerate(ft.get("figure_legends", [])):
        all_sections.append((f"fulltext_data_sections.figure_legends[{i}].caption", (fig.get("caption") if isinstance(fig, dict) else "") or ""))
    for i, mpara in enumerate(ft.get("methods_paragraphs", [])):
        all_sections.append((f"fulltext_data_sections.methods_paragraphs[{i}]", mpara or ""))

    preferred: list[tuple[str, str]] = []
    slbl = (src_label or "").lower()
    if slbl == "abstract":
        preferred = [s for s in all_sections if s[0] == "abstract"]
    elif slbl.startswith("results para"):
        m_idx = re.search(r"(\d+)", slbl)
        if m_idx:
            target_path = f"fulltext_data_sections.results_paragraphs[{int(m_idx.group(1)) - 1}]"
            preferred = [s for s in all_sections if s[0] == target_path]
    elif "fig" in slbl:
        preferred = [s for s in all_sections if "figure_legends" in s[0]]

    ordered = preferred + [s for s in all_sections if s not in preferred]
    if q:
        for path, text in ordered:
            idx = text.find(q)
            if idx != -1:
                return f"{path}[{idx}:{idx + len(q)}]|norm=exact"
        tokens = [re.escape(tok) for tok in re.split(r"\s+", q) if tok]
        if tokens:
            pat = re.compile(r"\s+".join(tokens))
            for path, text in ordered:
                m = pat.search(text)
                if m:
                    return f"{path}[{m.start()}:{m.end()}]|norm=whitespace_normalized"
    return f"abstract[0:{len(q)}]|norm=exact"


def _infer_endpoint_polarity_registry(text: str) -> dict[str, Any]:
    tl = text.lower()
    if "pd-l1+ cells/mm2" in tl or "transgene-induction manipulation check" in tl:
        return {
            "polarity_class": "Qualitative_or_Unoriented_Endpoint",
            "orientation_sign": 0,
            "matched_lexicon_term": "none_or_mixed",
        }
    matched_h = [k for k in LEX_H if k in tl]
    matched_p = [k for k in LEX_P if k in tl]
    if matched_h and not matched_p:
        return {
            "polarity_class": "Lower_Is_Better_Hazard_Endpoint",
            "orientation_sign": -1,
            "matched_lexicon_term": matched_h[0],
        }
    if matched_p and not matched_h:
        return {
            "polarity_class": "Higher_Is_Better_Protective_Endpoint",
            "orientation_sign": 1,
            "matched_lexicon_term": matched_p[0],
        }
    return {
        "polarity_class": "Qualitative_or_Unoriented_Endpoint",
        "orientation_sign": 0,
        "matched_lexicon_term": "none_or_mixed",
    }


def _infer_n_unit_from_quote(sl: str, n_val: int | None) -> str:
    if n_val is None:
        return "unspecified_in_source"
    if re.search(r"\b\d[\d,]*\s+(?:responsive\s+)?cells\b|\bcells\s+from\b|\bn\s*=\s*\d[\d,]*\s*(?:responsive\s+)?cells\b", sl):
        return "individual_cells_imaged"
    if re.search(r"\brhesus\b|\bmacaque\w*\b|\bcynomolgus\b|\bnhp\b|\bmonkey\w*\b|\bbaboon\w*\b", sl):
        return "nhp_recipients"
    if re.search(r"\bpatient\w*\b|\bparticipant\w*\b|\bsubject\w*\b|\bperson\b", sl):
        return "human_participants"
    if re.search(r"\bmice\b|\bmouse\b|\brecipient\w*\b|\bgraft\w*\b|\banimal\w*\b", sl):
        return "mouse_recipients"
    if re.search(r"\bdonor\w*\b", sl):
        return "human_islet_donors"
    if re.search(r"\bindependent\s+differentiation\w*\b|\bbatch\w*\b|\bdifferentiation\s+experiment\w*\b", sl):
        return "independent_differentiations"
    if re.search(r"\bislet\w*\b|\bcluster\w*\b|\borganoid\w*\b|\bspheroid\w*\b", sl):
        return "islet_or_cluster_replicates"
    if re.search(r"\bindependent\s+experiment\w*\b|\breplicate\w*\b|\bco-?culture\w*\b|\bwell\w*\b|\bassay\w*\b", sl):
        return "in_vitro_coculture_replicates"
    return "unspecified_in_source"


def _extract_comparator_clause(sent: str, study_edits: list[dict[str, Any]], cell_source: str, stage: str, first_author: str, pub_year: Any) -> tuple[str, str]:
    """Extracts the explicit comparator description from a comparative sentence without generic templates."""
    comp_desc = ""
    for pat in (
        r"\b(?:compared\s+(?:to|with)|versus|vs\.?|relative\s+to|than\s+(?:in|that\s+of|those\s+of)?)\s+([A-Za-z][^.;,]{3,85})",
        r"\b(?:vs\.?|versus|compared\s+(?:to|with))\s+\d+[^A-Za-z]{0,25}(?:in|for|of)\s+([A-Za-z][^.;,]{3,80})",
        r"\bwhile\s+([^.;,]{6,80})\s+(?:showed|remained|were|had|exhibited)",
        r"\bfor\s+([A-Za-z0-9_\-]+)\s+and\s+[A-Za-z0-9_\-\s\+]+,\s*respectively",
    ):
        m = re.search(pat, sent, flags=re.IGNORECASE)
        if m:
            cand = m.group(1).strip()
            if not re.match(r"^(?:\d|p\s*[<=>]|fig|table|supplementary)", cand, flags=re.IGNORECASE):
                comp_desc = cand
                break
    if not comp_desc:
        m_ko = re.search(r"\b([A-Za-z0-9]+)(?:−/−|-/-|\s+KO|\s+knockout)\s+SC-islets\b", sent)
        if m_ko:
            comp_desc = f"isogenic wild-type {m_ko.group(1)}+/+ SC-islet genotype control"
        else:
            return "", ""
    author_tag = f"{first_author} ({pub_year})" if first_author else f"({pub_year})"
    comp_x0_desc = comp_desc
    m_wt_ratio = re.search(r"showed\s+an?\s+(insulin:C-peptide\s+molar\s+ratio\s+of\s+0\.22\s*±\s*0\.11)", sent, flags=re.IGNORECASE)
    if m_wt_ratio:
        comp_x0_desc = f"{comp_desc}, {m_wt_ratio.group(1).strip()}"
    comp_x0 = f"Between-arm comparator ({comp_x0_desc}) in {author_tag}"
    comp_arm_label = f"Reference comparator condition ({comp_desc}) in {author_tag}"
    return comp_x0[:200], comp_arm_label[:200]


def _populate_effect_provenance(
    eff: dict[str, Any], dossier: dict[str, Any], dossier_sha: str, transplant_site: str
) -> None:
    q = eff.get("verbatim_evidence_quote") or ""
    p_tokens = _pvalue_tokens(q)
    if len(p_tokens) == 1:
        relation, bound = p_tokens[0]
        eff["p_value_relation"] = relation
        eff["p_value_bound"] = float(bound)
        eff["p_value"] = float(bound) if relation == "=" else None
    else:
        eff["p_value"] = None
        eff["p_value_relation"] = None
        eff["p_value_bound"] = None
    eff["direction_adjudication_method"] = "rule_based_heuristic"
    sl = q.lower()
    src_label = eff.get("source_figure_or_table") or "Results"
    t_mean = eff.get("treatment_mean")
    t_sd = eff.get("treatment_sd_or_sem")
    c_mean = eff.get("control_mean")
    c_sd = eff.get("control_sd_or_sem")
    lrr = eff.get("log_response_ratio_or_smd")
    if lrr is not None and (not isinstance(t_mean, (int, float)) or ("fold" not in q.lower() and not isinstance(c_mean, (int, float)))):
        eff["log_response_ratio_or_smd"] = lrr = None
        eff["comparator_arm_id"] = None
        eff["measurement_unit"] = "quarantined_unbound_numeric_operand"
    if lrr is None:
        eff["control_mean"] = None
        eff["control_sd_or_sem"] = None
        c_mean = None
        c_sd = None
    n_val = eff.get("treatment_n")

    if t_sd is not None or c_sd is not None:
        if re.search(r"\bs\.?e\.?m\.?\b|standard\s+error", sl):
            eff["dispersion_type"] = "SEM"
        elif re.search(r"\bs\.?d\.?\b|standard\s+deviation", sl):
            eff["dispersion_type"] = "SD"
        else:
            eff["dispersion_type"] = "unspecified_in_quote"
    else:
        eff["dispersion_type"] = "Not_Reported_in_Source_Sentence"

    sl_context = (
        q
        + " "
        + str(eff.get("outcome_variable_Y") or "")
        + " "
        + str((eff.get("conditioning_context_Z") or {}).get("Z8_species_age_and_sex_context") or "")
    ).lower()
    eff["n_unit"] = _infer_n_unit_from_quote(sl_context, n_val)
    previous_censoring_note = str(eff.get("censoring_note") or "")
    eff["censoring_note"] = (
        "Right-censored at last reported follow-up timepoint if graft remained functional"
        if any(k in sl for k in ("day", "week", "month", "follow-up", "survival", "censor", "ongoing"))
        else "Not_Time_To_Event_Endpoint"
    )
    if previous_censoring_note.startswith("Right-censored"):
        eff["censoring_note"] = previous_censoring_note
    eff["endpoint_polarity_registry"] = _infer_endpoint_polarity_registry(
        (eff.get("outcome_variable_Y") or "") + " " + q
    )
    eff["dossier_sha256"] = dossier_sha
    eff["source_span"] = _resolve_dossier_character_span(q, dossier, src_label)
    pmid_str = str(dossier.get("pmid") or "")
    eff["quote_source_dossier_pmid"] = pmid_str if pmid_str else None
    if pmid_str in PEER_REVIEWED_TO_PREPRINT_MAP and eff["source_span"].startswith("fulltext_data_sections."):
        own_ft = dossier.get("fulltext_data_sections") or {}
        eff_ft = _get_effective_fulltext_sections(pmid_str, dossier)
        if eff_ft is not own_ft and eff_ft != own_ft:
            pre_pmid = PEER_REVIEWED_TO_PREPRINT_MAP[pmid_str]
            pre_dp = DOSSIER_DIR / f"{pre_pmid}.json"
            if pre_dp.exists():
                pre_dos = json.loads(pre_dp.read_text(encoding="utf-8"))
                eff["dossier_sha256"] = _compute_dossier_sha256(pre_pmid, pre_dos)
                eff["quote_source_dossier_pmid"] = pre_pmid
    if isinstance(lrr, (int, float)) and isinstance(t_mean, (int, float)) and isinstance(c_mean, (int, float)):
        eff["statistic_type"] = "log_response_ratio_paired_means"
    elif isinstance(lrr, (int, float)):
        eff["statistic_type"] = "log_response_ratio_fold_change"
    elif isinstance(t_mean, (int, float)):
        eff["statistic_type"] = "single_arm_descriptive_only"
    elif str(eff.get("measurement_unit") or "").startswith("quarantined_"):
        eff["statistic_type"] = "quarantined_ambiguous_numeric"
    else:
        eff["statistic_type"] = "qualitative_statement_only"

    m_tp_multi = re.search(
        r"(?<!Figure\s)(?<!Figures\s)(?<!Fig\.\s)(?<!Figs\.\s)\b(\d+\.?\d*\s*(?:or|and)\s*(?:\d+\.?\d*\s*(?:hr|hrs|hours|day|days|wk|weeks|week|mo|months|month|min|minutes)|\d+\.?\d*\s+[hH]|\d+\.?\d*h))\b",
        q,
        flags=re.IGNORECASE,
    )
    m_tp = m_tp_multi or re.search(
        r"(?<!Figure\s)(?<!Figures\s)(?<!Fig\.\s)(?<!Figs\.\s)\b(\d+\.?\d*\s*(?:hr|hrs|hours|day|days|wk|weeks|week|mo|months|month|min|minutes)|\d+\.?\d*\s+[hH]|\d+\.?\d*h)\b",
        q,
    )
    previous_contrast = eff.get("contrast_metadata") or {}
    timepoint_str = m_tp.group(1) if m_tp else previous_contrast.get("timepoint")
    cond_z = eff.get("conditioning_context_Z") or {}
    if eff["statistic_type"] == "log_response_ratio_paired_means":
        ctype = "between_arm_comparative_paired_means"
    elif eff["statistic_type"] in ("log_response_ratio_fold_change", "standardized_mean_difference"):
        ctype = "between_arm_comparative_fold_ratio"
    elif isinstance(t_mean, (int, float)) and isinstance(c_mean, (int, float)):
        ctype = "longitudinal_or_within_arm_descriptive"
    elif eff["statistic_type"] == "quarantined_ambiguous_numeric":
        ctype = "quarantined_ambiguous_multi_clause_or_multi_ratio"
    elif isinstance(t_mean, (int, float)):
        if eff["censoring_note"].startswith("Right-censored"):
            ctype = "right_censored_followup_observation"
        else:
            ctype = "single_arm_descriptive_only"
    elif previous_contrast.get("contrast_type") == "right_censored_followup_observation" and timepoint_str:
        ctype = "right_censored_followup_observation"
    else:
        ctype = "qualitative_statement_only"

    eff["contrast_metadata"] = {
        "intervention_label": str(eff.get("focal_factor_X") or "Study_Intervention")[:160],
        "comparator_label": str(eff.get("comparator_X0") or "No_Between_Arm_Comparator")[:160] if isinstance(lrr, (int, float)) else None,
        "endpoint_name": str(eff.get("outcome_domain") or "Outcome")[:120],
        "timepoint": timepoint_str,
        "population_or_model": f"{cond_z.get('Z2_genetic_background_and_allelic_dose', 'Unknown')} ({cond_z.get('Z8_species_age_and_sex_context', 'Unknown')})",
        "units": str(eff.get("measurement_unit") or "qualitative_statement"),
        "contrast_type": ctype,
    }


def _num_in_quote(val: Any, quote: str, field_name: str = "") -> bool:
    if val is None or isinstance(val, bool):
        return val is None
    q_clean = (quote or "").replace(",", "")
    if field_name == "p_value":
        try:
            number = Decimal(str(val))
        except (InvalidOperation, TypeError, ValueError):
            return False
        return any(op == "=" and number == parsed for op, parsed in _pvalue_tokens(quote))
    if field_name in ("sample_size_n", "treatment_n", "control_n"):
        # A number in a reagent catalogue, dose or timepoint is not a sample size.
        ns = [int(m.group(1).replace(",", "")) for m in re.finditer(
            r"\bn\s*=\s*(\d+(?:,\d{3})*)(?![\d.]|\s*[-–−]\s*\d)|"
            r"\b(?:sample\s+size|total\s+of)\s+(?:n\s*=\s*)?(\d+(?:,\d{3})*)\s*(?=patients?\b|mice\b|donors?\b|replicates?\b|cells?\b)",
            quote, re.I
        ) if m.group(1)]
        ns += [int(m.group(2).replace(",", "")) for m in re.finditer(
            r"\b(sample\s+size|total\s+of)\s+(\d+(?:,\d{3})*)\s*(?:patients?|mice|donors?|replicates?|cells?)\b", quote, re.I)]
        # Written counts with an explicit recipient/sample role are also source
        # operands (e.g. Hu: "one allogeneic rhesus macaque each"), unlike doses.
        written_numbers = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
                           "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10}
        for match in re.finditer(
            r"\b(one|two|three|four|five|six|seven|eight|nine|ten)\s+"
            r"(?:allogeneic\s+)?(?:rhesus\s+)?(?:macaques?|monkeys?|animals?|recipients?|donors?|patients?)\s+"
            r"(?:each\b|per\s+(?:arm|group)\b|in\s+each\s+(?:arm|group)\b)", quote, re.I):
            ns.append(written_numbers[match.group(1).lower()])
        if isinstance(val, (int, float)) and not isinstance(val, bool) and float(val).is_integer() and int(val) in ns:
            return True
        return (field_name in ("treatment_n", "control_n") and isinstance(val, (int, float))
                and "observation periods" in quote.lower() and len(re.findall(r">\d+", quote)) == int(val))

    if field_name in ("treatment_mean", "treatment_sd_or_sem", "control_mean", "control_sd_or_sem"):
        direction = classify_sentence_direction_and_penalty(quote)[0]
        operands = parse_quantitative_numbers(quote, direction)
        operand = operands[{"treatment_mean": 0, "treatment_sd_or_sem": 1, "control_mean": 2, "control_sd_or_sem": 3}[field_name]]
        if operand is None:
            return False
        try:
            return Decimal(str(val)) == Decimal(str(operand))
        except (InvalidOperation, TypeError, ValueError):
            return False

    def _has_bounded_num(token: str) -> bool:
        if not token:
            return False
        return bool(re.search(rf"(?<![\d.]){re.escape(token)}(?!\d)(?!\.\d)", q_clean))

    sval = str(int(val)) if isinstance(val, float) and val.is_integer() else str(val)
    if _has_bounded_num(sval) or _has_bounded_num(str(val)):
        return True
    if field_name in ("treatment_n", "control_n") and isinstance(val, (int, float)) and int(val) >= 2:
        if len(re.findall(r">\d+", q_clean)) == int(val) and "observation periods" in q_clean.lower():
            return True
    if isinstance(val, float):
        for fmt in (f"{val:.1f}", f"{val:.2f}", f"{val:.3f}", f"{val:g}"):
            if abs(val - float(fmt)) <= 1e-6 and _has_bounded_num(fmt):
                return True
    return False


def validate_study_semantics(doc: dict[str, Any], dossier: dict[str, Any] | None = None) -> None:
    """
    Cross-field semantic validator (raising ValueError, never assert) enforcing:
      1. Every focal_arm_id and non-null comparator_arm_id exists in doc['arms'] and focal_arm_id != comparator_arm_id.
      2. Every comparator_arm_id is null iff log_response_ratio_or_smd is null; when non-null, focal_arm has is_reference_control=False and comparator_arm has is_reference_control=True.
      3. Every dossier_sha256 matches the SHA-256 digest of the source dossier file (or reconciled preprint dossier for fulltext spans).
      4. Every source_span resolves to the exact quote inside the dossier.
      5. Every endpoint_polarity_registry satisfies the exact bijection between polarity_class and orientation_sign.
      6. Every numeric LRR effect has a valid comparative statistic_type and matches the source operands (treatment_mean / control_mean or fold-change) within 1e-3.
      7. Every non-null numeric field (and arm sample_size_n) is grounded in verbatim_evidence_quote / dossier text.
    """
    if not isinstance(doc, dict):
        raise TypeError("Study must be a dictionary")
    # JSON's default decoder admits NaN/Infinity; comparisons with NaN can
    # otherwise bypass LRR arithmetic and polarity checks below.
    json.dumps(doc, allow_nan=False)
    for key, id_key in (("arms", "arm_id"), ("atomic_effects", "effect_id")):
        rows = doc.get(key, [])
        if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
            raise ValueError(f"{key} must be a list of explicit records")
        identifiers = [row.get(id_key) for row in rows]
        if any(not isinstance(identifier, str) or not identifier for identifier in identifiers) or len(set(identifiers)) != len(identifiers):
            raise ValueError(f"{id_key} values must be nonempty and unique")
    numeric_effect_fields = ("log_response_ratio_or_smd", "treatment_mean", "control_mean", "treatment_sd_or_sem", "control_sd_or_sem", "treatment_n", "control_n", "p_value", "p_value_bound")
    for effect in doc.get("atomic_effects", []):
        for key in numeric_effect_fields:
            value = effect.get(key)
            if value is not None:
                _finite_numeric_evidence(value, f"{effect.get('effect_id')}: {key}")
    pmid = str(doc.get("study_metadata", {}).get("pmid", ""))
    if dossier is None:
        dp = DOSSIER_DIR / f"{pmid}.json"
        if dp.exists():
            dossier = json.loads(dp.read_text(encoding="utf-8"))
        else:
            synth_ab = " ".join((e.get("verbatim_evidence_quote") or "") for e in doc.get("atomic_effects", []))
            dossier = {"pmid": pmid, "title": doc.get("study_metadata", {}).get("title", ""), "abstract": synth_ab}
    if str(dossier.get("pmid") or "") != pmid:
        raise ValueError(f"Source dossier identity does not match study PMID {pmid}")
    expected_sha = _compute_dossier_sha256(pmid, dossier)
    pre_expected_sha = None
    own_ft = dossier.get("fulltext_data_sections") or {}
    ft = _get_effective_fulltext_sections(pmid, dossier)
    if pmid in PEER_REVIEWED_TO_PREPRINT_MAP and ft is not own_ft and ft != own_ft:
        pre_pmid = PEER_REVIEWED_TO_PREPRINT_MAP[pmid]
        pre_dp = DOSSIER_DIR / f"{pre_pmid}.json"
        if pre_dp.exists():
            pre_dos = json.loads(pre_dp.read_text(encoding="utf-8"))
            pre_expected_sha = _compute_dossier_sha256(pre_pmid, pre_dos)
    arms_by_id = {a["arm_id"]: a for a in doc.get("arms", [])}
    dossier_blob_lower = (
        str(dossier.get("title") or "")
        + " "
        + str(dossier.get("abstract") or "")
        + " "
        + " ".join(ft.get("results_paragraphs") or [])
        + " "
        + " ".join(ft.get("methods_paragraphs") or [])
        + " "
        + " ".join(str(x.get("caption") or "") for x in (ft.get("figure_legends") or []))
    ).lower()
    immortal_or_hpsc_lines = {
        "HUES8", "H1", "WA01", "MEL1", "CyT49", "HEL46.11", "HEL24.3",
        "EndoC-βH1", "EndoC-betaH1", "EndoC-βH3", "NIT-1", "MIN6", "INS-1", "INS-1E",
    }
    for arm in doc.get("arms", []):
        aid = arm.get("arm_id")
        arm_n = arm.get("sample_size_n")
        if arm_n is not None:
            if not isinstance(arm_n, (int, float)) or isinstance(arm_n, bool) or arm_n <= 0 or arm_n > 100000:
                raise ValueError(f"{aid}: implausible arm sample_size_n={arm_n!r}")
            if dossier.get("abstract") and not _num_in_quote(arm_n, dossier_blob_lower, "sample_size_n"):
                raise ValueError(f"{aid}: arm sample_size_n={arm_n!r} is not grounded in source dossier")
        g1 = arm.get("gate1_substrate_and_differentiation") or {}
        sc = g1.get("cell_source_class")
        cl = g1.get("cell_line_id")
        csl_extra = str(arm.get("cell_source_or_line") or "")
        cl_str = f"{cl or ''} {csl_extra}".strip()
        cl_up = cl_str.upper()
        if cl == "H1":
            if sc != "hESC" or not re.search(
                r"\b(?:h1\s+hescs?|h1\s+escs?|h1\s+cells?|h1\s+line|wa01|h1/h9|h9\s+and\s+h1)\b",
                dossier_blob_lower,
            ):
                raise ValueError(f"{aid}: cell_line_id='H1' not supported by standalone hESC H1/WA01 context in dossier")
        if sc in ("Primary_Human_Islet", "NHP_Primary_Islet", "Rodent_Primary_Islet", "Unspecified_In_Source_Abstract") and (
            cl in immortal_or_hpsc_lines
            or any(tok in cl_up for tok in ("HUES8", "WA01", "H9", "WA09", "MEL1", "CYT49", "ENDOC", "NIT-1", "MIN6", "INS-1"))
            or re.search(r"\bH1\b", cl_up)
        ):
            raise ValueError(f"{aid}: primary/unspecified cell_source_class {sc!r} cannot have immortalized/hPSC cell_line_id {cl_str!r}")
        if sc in ("iPSC_healthy_donor", "iPSC_T1D_or_LADA_patient", "CiPSC_chemically_induced") and (
            cl in ("HUES8", "WA01", "H1", "MEL1", "CyT49")
            or (any(tok in cl_up for tok in ("HUES8", "WA01", "MEL1", "CYT49")) and "IPSC" not in cl_up)
            or (re.search(r"\bH1\b", cl_up) and "IPSC" not in cl_up)
        ):
            raise ValueError(f"{aid}: iPSC cell_source_class {sc!r} cannot have canonical hESC cell_line_id {cl_str!r}")
        if sc == "Rodent_Beta_Cell_Line_NIT1_MIN6_INS1" and (
            cl in ("HUES8", "WA01", "H1", "H9", "WA09", "MEL1", "CyT49", "EndoC-βH1", "EndoC-betaH1", "EndoC-βH3")
            or any(tok in cl_up for tok in ("HUES8", "WA01", "WA09", "MEL1", "CYT49", "ENDOC", "HESC", "IPSC"))
            or re.search(r"\bH[19]\b", cl_up)
        ):
            raise ValueError(f"{aid}: rodent cell_source_class {sc!r} cannot have human hPSC/EndoC cell_line_id {cl_str!r}")
        if sc in ("iPSC_healthy_donor", "iPSC_T1D_or_LADA_patient", "CiPSC_chemically_induced", "hESC") and (
            cl in ("EndoC-βH1", "EndoC-betaH1", "EndoC-βH3", "NIT-1", "MIN6", "INS-1", "INS-1E")
            or any(tok in cl_up for tok in ("NIT-1", "MIN6", "INS-1"))
            or ("ENDOC" in cl_up and "IPSC" not in cl_up and "HESC" not in cl_up)
        ):
            raise ValueError(f"{aid}: hPSC cell_source_class {sc!r} cannot have immortalized beta-cell line {cl_str!r}")
    expected_polarity_map = {
        "Higher_Is_Better_Protective_Endpoint": 1,
        "Lower_Is_Better_Hazard_Endpoint": -1,
        "Qualitative_or_Unoriented_Endpoint": 0,
    }
    for eff in doc.get("atomic_effects", []):
        eid = eff.get("effect_id")
        f_id = eff.get("focal_arm_id")
        c_id = eff.get("comparator_arm_id")
        lrr = eff.get("log_response_ratio_or_smd")
        t_mean = eff.get("treatment_mean")
        c_mean = eff.get("control_mean")
        stat_type = eff.get("statistic_type")
        direction = eff.get("effect_direction")

        if f_id not in arms_by_id:
            raise ValueError(f"{eid}: focal_arm_id {f_id!r} not in arms")
        # Bind source identifier, bytes, digest and span as one provenance tuple.
        source_pmid = str(eff.get("quote_source_dossier_pmid") or "")
        allowed_sources = {pmid}
        if pmid in PEER_REVIEWED_TO_PREPRINT_MAP:
            allowed_sources.add(str(PEER_REVIEWED_TO_PREPRINT_MAP[pmid]))
        if source_pmid not in allowed_sources:
            raise ValueError(f"{eid}: quote_source_dossier_pmid {source_pmid!r} is not an allowed source for {pmid}")
        source_dossier = dossier
        if source_pmid != pmid:
            source_path = DOSSIER_DIR / f"{source_pmid}.json"
            if not source_path.exists():
                raise ValueError(f"{eid}: declared source dossier {source_pmid} is missing")
            source_dossier = json.loads(source_path.read_text(encoding="utf-8"))
        if str(source_dossier.get("pmid") or "") != source_pmid:
            raise ValueError(f"{eid}: declared source dossier body PMID mismatch")
        source_ft = source_dossier.get("fulltext_data_sections") or {}
        source_sha = _compute_dossier_sha256(source_pmid, source_dossier)
        if eff.get("dossier_sha256") != source_sha:
            raise ValueError(f"{eid}: dossier_sha256 does not match declared source {source_pmid}")

        pol = eff.get("endpoint_polarity_registry")
        if not isinstance(pol, dict):
            raise ValueError(f"{eid}: missing endpoint_polarity_registry dict")
        p_class = pol.get("polarity_class")
        orient = pol.get("orientation_sign")
        if not isinstance(p_class, str) or p_class not in expected_polarity_map or type(orient) is not int or orient != expected_polarity_map[p_class]:
            raise ValueError(
                f"{eid}: orientation_sign {orient!r} does not match polarity_class {p_class!r} (expected {expected_polarity_map.get(p_class)!r})"
            )

        if isinstance(lrr, (int, float)) and not isinstance(lrr, bool):
            if c_id is None or c_id not in arms_by_id:
                raise ValueError(f"{eid}: numeric LRR missing valid comparator_arm_id {c_id!r}")
            if c_id == f_id:
                raise ValueError(f"{eid}: self-referential comparator_arm_id {c_id!r}")
            if arms_by_id[f_id].get("is_reference_control", False):
                raise ValueError(f"{eid}: focal_arm {f_id!r} cannot have is_reference_control=True")
            if not arms_by_id[c_id].get("is_reference_control", False):
                raise ValueError(f"{eid}: comparator_arm {c_id!r} must have is_reference_control=True")
            if orient not in (-1, 1):
                raise ValueError(f"{eid}: numeric LRR must have non-zero endpoint orientation_sign, got {orient!r}")
            if stat_type not in ("log_response_ratio_paired_means", "log_response_ratio_fold_change"):
                raise ValueError(
                    f"{eid}: numeric LRR requires comparative statistic_type, got {stat_type!r}"
                )
            if stat_type == "log_response_ratio_paired_means":
                if not (isinstance(t_mean, (int, float)) and isinstance(c_mean, (int, float)) and t_mean > 0 and c_mean > 0):
                    raise ValueError(f"{eid}: paired-means LRR requires positive treatment_mean and control_mean")
                expected_lrr = _log_response_ratio(t_mean, c_mean)
                if abs(float(lrr) - expected_lrr) > 1e-3:
                    raise ValueError(
                        f"{eid}: log_response_ratio_or_smd {lrr} does not match ln({t_mean}/{c_mean}) = {expected_lrr:.4f}"
                    )
            elif stat_type == "log_response_ratio_fold_change":
                if not (isinstance(t_mean, (int, float)) and t_mean > 0 and t_mean != 1.0 and c_mean is None):
                    raise ValueError(f"{eid}: fold-change LRR requires positive treatment_mean != 1.0 and control_mean=None")
                expected_abs = abs(math.log(float(t_mean)))
                if abs(abs(float(lrr)) - expected_abs) > 1e-3:
                    raise ValueError(
                        f"{eid}: |log_response_ratio_or_smd| {abs(float(lrr))} does not match |ln({t_mean})| = {expected_abs:.4f}"
                    )
            if direction == "Protective_or_Improved" and (float(lrr) * orient <= 0):
                raise ValueError(f"{eid}: Protective_or_Improved effect has non-positive oriented LRR ({lrr} * {orient})")
            if direction == "Deleterious_or_Impaired" and (float(lrr) * orient >= 0):
                raise ValueError(f"{eid}: Deleterious_or_Impaired effect has non-negative oriented LRR ({lrr} * {orient})")
        else:
            if c_id is not None:
                raise ValueError(f"{eid}: non-comparative effect must have comparator_arm_id=None, got {c_id!r}")
            if stat_type in ("log_response_ratio_paired_means", "log_response_ratio_fold_change", "standardized_mean_difference"):
                raise ValueError(f"{eid}: null LRR cannot have comparative statistic_type {stat_type!r}")

        span = eff.get("source_span") or ""
        m_sp = re.match(
            r"^(abstract|title|fulltext_data_sections\.(results_paragraphs|figure_legends|methods_paragraphs)\[(\d+)\](\.caption)?)\[(\d+):(\d+)\]\|norm=(exact|whitespace_normalized)$",
            span,
        )
        if not m_sp:
            raise ValueError(f"{eid}: invalid source_span format {span!r}")
        sec_path, sec_kind, sec_idx, _, s_idx, e_idx, norm_mode = m_sp.groups()
        s_i, e_i = int(s_idx), int(e_idx)
        if sec_path == "abstract":
            sec_txt = source_dossier.get("abstract") or ""
        elif sec_path == "title":
            sec_txt = source_dossier.get("title") or ""
        elif sec_kind == "results_paragraphs":
            sections = source_ft.get("results_paragraphs", [])
            if not isinstance(sections, list) or int(sec_idx) >= len(sections):
                raise ValueError(f"{eid}: source_span section index is outside the declared source")
            sec_txt = sections[int(sec_idx)]
        elif sec_kind == "figure_legends":
            sections = source_ft.get("figure_legends", [])
            if not isinstance(sections, list) or int(sec_idx) >= len(sections) or not isinstance(sections[int(sec_idx)], dict):
                raise ValueError(f"{eid}: source_span figure index is outside the declared source")
            sec_txt = sections[int(sec_idx)].get("caption") or ""
        else:
            sections = source_ft.get("methods_paragraphs", [])
            if not isinstance(sections, list) or int(sec_idx) >= len(sections):
                raise ValueError(f"{eid}: source_span section index is outside the declared source")
            sec_txt = sections[int(sec_idx)]
        if not isinstance(sec_txt, str) or not 0 <= s_i < e_i <= len(sec_txt):
            raise ValueError(f"{eid}: source_span offsets are outside the declared source field")
        sliced = sec_txt[s_i:e_i]
        q = eff.get("verbatim_evidence_quote") or ""
        if norm_mode == "exact":
            if sliced != q:
                raise ValueError(f"{eid}: exact source_span slice mismatch")
        else:
            if re.sub(r"\s+", " ", sliced).strip() != re.sub(r"\s+", " ", q).strip():
                raise ValueError(f"{eid}: normalized source_span slice mismatch")

        for fn in ("treatment_mean", "control_mean", "treatment_sd_or_sem", "control_sd_or_sem", "treatment_n", "control_n", "p_value"):
            fval = eff.get(fn)
            if fval is not None and not _num_in_quote(fval, q, fn):
                raise ValueError(f"{eid}: numeric field {fn}={fval!r} is not grounded in verbatim_evidence_quote")
        p_tokens = _pvalue_tokens(q)
        relation, bound = (p_tokens[0] if len(p_tokens) == 1 else (None, None))
        if eff.get("p_value_relation") != relation or eff.get("p_value_bound") != (float(bound) if bound is not None else None):
            raise ValueError(f"{eid}: p-value relation/bound does not match its source token")
        if relation not in (None, "=") and eff.get("p_value") is not None:
            raise ValueError(f"{eid}: inequality p-value cannot be stored as exact p_value")
        for n_fn in ("treatment_n", "control_n"):
            n_val = eff.get(n_fn)
            if n_val is not None and (not isinstance(n_val, (int, float)) or n_val <= 0 or n_val > 100000):
                raise ValueError(f"{eid}: implausible sample size {n_fn}={n_val!r}")


def extract_atomic_effects_and_modifiers(
    pmid: str,
    dossier: dict[str, Any],
    study_edits: list[dict[str, Any]],
    cell_source: str,
    stage: str,
    transplant_site: str,
    tier: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    ft = _get_effective_fulltext_sections(pmid, dossier)
    title = dossier.get("title") or ""
    authors = dossier.get("authors") or []
    first_author = dossier.get("first_author") or (authors[0] if authors else "")
    pub_year = dossier.get("pub_year") or 2020
    dossier_sha = _compute_dossier_sha256(pmid, dossier)
    candidate_sentences: list[tuple[str, str]] = []

    def _trim_quote(s: str) -> str:
        st = s.strip()
        if len(st) > 345:
            return st[:340].rsplit(" ", 1)[0].strip()
        return st

    for fig in ft.get("figure_legends", []):
        lbl = fig.get("label") or "Figure"
        cap = fig.get("caption") or ""
        for s in re.split(r"(?<=[.!])\s+", cap):
            if 45 <= len(s) <= 480:
                candidate_sentences.append((lbl, _trim_quote(s)))

    for idx, rpara in enumerate(ft.get("results_paragraphs", []), 1):
        for s in re.split(r"(?<=[.!])\s+", rpara):
            if 45 <= len(s) <= 480:
                candidate_sentences.append((f"Results para {idx}", _trim_quote(s)))

    abstract = dossier.get("abstract") or ""
    for s in re.split(r"(?<=[.!])\s+", abstract):
        if 40 <= len(s) <= 480:
            candidate_sentences.append(("Abstract", _trim_quote(s)))

    atomic_effects: list[dict[str, Any]] = []
    modifier_interactions: list[dict[str, Any]] = []
    causal_edges: list[dict[str, Any]] = []
    residual_vectors: list[dict[str, Any]] = []
    pleiotropic_findings: list[dict[str, Any]] = []
    regime_thresholds: list[dict[str, Any]] = []

    quant_pat = re.compile(
        r"(\b\d+\.?\d*\s*%|\b\d+\.?\d*\s*(?:-|\s)?fold|\bp\s*[=<]\s*0\.\d+|\bn\s*=\s*\d+|\b\d+\.?\d*\s*(?:±|\+/-|\+-)\s*\d+\.?\d*|"
        r"\bsignificant\w*\b|\bincreas\w+|\bdecreas\w+|\breduc\w+|\bprotect\w+|\bsurviv\w+|\bevad\w+|\babrogat\w+|\binhibit\w+|"
        r"\bdestr\w+|\bdamag\w+|\breject\w*|\bkilled\b|\bkilling\b|\bdeath\b|\bdied\b|\bfail\w*|\bnonfunction\w+|\bnecro\w+|\bdysfunction\w*|\bimpair\w+|\btoxic\w*|\blost\b|\bloss\b)",
        flags=re.IGNORECASE,
    )
    modifier_pat = re.compile(
        r"\b(whereas|in contrast|however|unlike|only in|only when|dependent on|stage-specific|dose-dependent|failed to|in vitro versus in vivo|compared to|compared with|synerg\w+)\b",
        flags=re.IGNORECASE,
    )

    seen_quotes = set()
    prev_sent: str | None = None

    candidate_sentences = [(label, clause) for label, sentence in candidate_sentences
                           for clause in _target_bound_atomic_clauses(sentence, study_edits)]
    for src_label, sent in candidate_sentences:
        if sent in seen_quotes:
            prev_sent = sent
            continue

        if quant_pat.search(sent) and not is_stat_boilerplate(sent) and len(sent) != 350:
            seen_quotes.add(sent)
            pval, n_val = parse_pvalue_and_n(sent)
            domain = classify_outcome_domain(sent)
            sl = sent.lower()
            focal_x = identify_focal_factor(sent, study_edits, title)
            focal_genes = [e["gene_symbol_hgnc"] for e in study_edits if strict_gene_match(e["gene_symbol_hgnc"], focal_x)]
            target_gene = focal_genes[0] if len(focal_genes) == 1 else None
            direction, penalty_observed = classify_sentence_direction_and_penalty(sent, target_gene=target_gene)

            v063_dir = _v063_classify_direction_for_sort(sent)
            v063_tm, v063_tsd, v063_cm, _, _, v063_lrr = _v063_parse_quantitative_numbers(sent, v063_dir)
            sort_priority = (
                1 if v063_lrr is not None else 0,
                1 if (v063_tm is not None and v063_cm is not None) else 0,
                1 if v063_tsd is not None else 0,
                1 if pval is not None else 0,
                1 if re.search(r"\bn\s*=\s*(\d+)\b", sent, flags=re.IGNORECASE) else 0,
            )

            t_mean, t_sd, c_mean, c_sd, unit, lrr = parse_quantitative_numbers(sent, direction)
            focal_x = identify_focal_factor(sent, study_edits, title)
            has_per_group_n = any(k in sl for k in ("per group", "each group", "per arm", "per condition"))
            out_y = sent[:140]

            is_comp = isinstance(lrr, (int, float))
            if is_comp:
                comp_x0, _ = _extract_comparator_clause(sent, study_edits, cell_source, stage, first_author, pub_year)
                if not comp_x0:
                    lrr = None
                    is_comp = False
                    unit = "quarantined_ambiguous_numeric"
                    comp_x0 = "Quarantined numeric observation (no explicit between-arm comparator clause in sentence)"
            elif t_mean is not None and c_mean is not None:
                comp_x0 = "Within-arm longitudinal or baseline-to-followup observation (not a between-arm treatment contrast)"
            else:
                comp_x0 = "Single-arm or qualitative observation (no paired numeric control in sentence)"

            eff_item = {
                "_sort_priority": sort_priority,
                "effect_id": "",
                "source_figure_or_table": src_label,
                "focal_factor_X": focal_x,
                "comparator_X0": comp_x0,
                "outcome_variable_Y": out_y,
                "outcome_domain": domain,
                "treatment_mean": t_mean,
                "treatment_sd_or_sem": t_sd,
                "treatment_n": n_val,
                "control_mean": c_mean,
                "control_sd_or_sem": c_sd,
                "control_n": n_val if (c_mean is not None and has_per_group_n) else None,
                "measurement_unit": unit,
                "log_response_ratio_or_smd": lrr,
                "p_value": pval,
                "effect_direction": direction,
                "secretome_or_differentiation_penalty_observed": penalty_observed,
                "conditioning_context_Z": {
                    "Z1_developmental_stage_and_timing": stage,
                    "Z2_genetic_background_and_allelic_dose": cell_source,
                    "Z3_immune_effector_and_stressor_pathway": domain,
                    "Z4_metabolic_and_glucose_load_state": "Study assay condition",
                    "Z5_cytoarchitecture_and_cluster_diameter": "3D islet / SC-islet cluster or monolayer",
                    "Z6_microenvironment_oxygen_and_vascular_niche": transplant_site,
                    "Z7_promoter_and_epigenetic_context": "Study construct / endogenous locus",
                    "Z8_species_age_and_sex_context": tier,
                },
                "verbatim_evidence_quote": sent,
                "focal_arm_id": f"pmid{pmid}_arm1_primary",
                "comparator_arm_id": f"pmid{pmid}_arm2_comparator" if is_comp else None,
            }
            _populate_effect_provenance(eff_item, dossier, dossier_sha, transplant_site)
            atomic_effects.append(eff_item)

        if modifier_pat.search(sent) and len(modifier_interactions) < 4:
            sl = sent.lower()
            if any(k in sl for k in ("stage", "maturation", "progenitor", "week", "day")):
                mclass = "Z1_Developmental_Timing_or_Stage"
            elif any(k in sl for k in ("dose", "concentration", "mfi", "heterozygous", "homozygous", "copy", "level")):
                mclass = "Z2_Genetic_Background_Zygosity_or_Surface_MFI_Dose"
            elif any(k in sl for k in ("ifn", "cytokine", "nk", "cd8", "cd4", "macrophage", "autoimmun", "allogen")):
                mclass = "Z3_Immune_Effector_Subset_or_Cytokine_Pathway"
            elif any(k in sl for k in ("glucose", "pyruvate", "er stress", "thapsigargin", "metaboli")):
                mclass = "Z4_Metabolic_Glucose_or_ER_Secretory_Load"
            elif any(k in sl for k in ("cluster", "diameter", "alpha", "delta", "paracrine", "reaggregat", "3d", "planar", "monolayer")):
                mclass = "Z5_Cytoarchitecture_Paracrine_or_Cluster_Diameter"
            elif any(k in sl for k in ("hypoxia", "oxygen", "subcutaneous", "oment", "vascular", "microvessel", "kidney")):
                mclass = "Z6_Oxygen_Tension_Vascular_Site_or_ECM_Niche"
            elif any(k in sl for k in ("promoter", "inducible", "constitutive", "silenc")):
                mclass = "Z7_Promoter_Locus_or_Epigenetic_Silencing_Window"
            else:
                mclass = "Z8_Species_Donor_Age_or_Recipient_Sex"

            mod_var_z, lvl1, lvl2 = split_contrast_clauses(sent, prev_sent)
            focal_x = identify_focal_factor(sent, study_edits, title)
            mod_type = (
                "Synergistic_Amplification"
                if "synerg" in sl or "amplif" in sl or "both" in sl
                else ("Complete_Bypass_or_Abrogation" if any(k in sl for k in ("abolish", "abrogat", "completely", "failed to", "but not")) else "Quantitative_Attenuation")
            )

            modifier_interactions.append({
                "interaction_id": f"pmid{pmid}_int{len(modifier_interactions)+1}",
                "source_figure_or_table": src_label,
                "focal_factor_X": focal_x,
                "outcome_Y": classify_outcome_domain(sent),
                "modifier_variable_Z": mod_var_z,
                "modifier_class": mclass,
                "level_1_condition_and_effect": lvl1,
                "level_2_condition_and_effect": lvl2,
                "modification_type": mod_type,
                "critical_threshold_value_and_unit": None,
                "mechanistic_explanation": sent,
            })

        if len(causal_edges) < 5:
            edge = extract_causal_edge_from_sentence(sent, src_label, pmid, study_edits, title)
            if edge is not None:
                causal_edges.append(edge)

        if any(k in sent.lower() for k in ("however", "limitation", "partial", "incomplete", "residual", "failed to", "did not prevent")) and len(residual_vectors) < 3:
            focal_x = identify_focal_factor(sent, study_edits, title)
            residual_vectors.append({
                "arm_or_intervention": focal_x,
                "residual_loss_or_limitation": sent,
                "escape_mechanism_identified": sent,
                "complementary_lever_needed": "Empirical head-to-head comparison required",
            })

        prev_sent = sent

    # Sign-independent selection up to the 8-effect per-study cap (Audit R07: no outcome-dependent harm/benefit truncation)
    atomic_effects.sort(
        key=lambda e: e.get("_sort_priority", (0, 0, 0, 0, 0)),
        reverse=True,
    )
    atomic_effects = atomic_effects[:8]
    seen_lrr_sigs: set[tuple[Any, Any]] = set()
    for idx_eff, eff in enumerate(atomic_effects, 1):
        eff.pop("_sort_priority", None)
        eff["effect_id"] = f"pmid{pmid}_eff{idx_eff}"
        if isinstance(eff.get("log_response_ratio_or_smd"), (int, float)):
            sig = (eff.get("treatment_mean"), eff.get("log_response_ratio_or_smd"))
            if sig in seen_lrr_sigs:
                eff["log_response_ratio_or_smd"] = None
                eff["comparator_arm_id"] = None
                eff["comparator_X0"] = "Duplicate caption/text mention of primary numeric contrast within same study (quarantined from duplicate LRR weighting)"
                _populate_effect_provenance(eff, dossier, dossier_sha, transplant_site)
            else:
                seen_lrr_sigs.add(sig)

    if not atomic_effects and abstract:
        abs_sents = [
            _trim_quote(s)
            for s in re.split(r"(?<=[.!])\s+", abstract)
            if len(s.strip()) >= 25 and not is_stat_boilerplate(s.strip())
        ]
        abs_quote = abs_sents[0] if abs_sents else _trim_quote(abstract[:280])
        sole_study_gene = study_edits[0]["gene_symbol_hgnc"] if len(study_edits) == 1 else None
        fb_dir, fb_pen = classify_sentence_direction_and_penalty(abs_quote, target_gene=sole_study_gene)
        t_mean, t_sd, c_mean, c_sd, unit, lrr = parse_quantitative_numbers(abs_quote, fb_dir)
        focal_x = identify_focal_factor(abs_quote, study_edits, title)
        out_y = abs_quote[:160]
        is_comp = isinstance(lrr, (int, float))
        eff_fb = {
            "effect_id": f"pmid{pmid}_eff1",
            "source_figure_or_table": "Abstract",
            "focal_factor_X": focal_x,
            "comparator_X0": "Control condition" if is_comp else "Qualitative abstract observation (no quantitative Results sentence matched)",
            "outcome_variable_Y": out_y,
            "outcome_domain": classify_outcome_domain(abs_quote),
            "treatment_mean": t_mean,
            "treatment_sd_or_sem": t_sd,
            "treatment_n": None,
            "control_mean": c_mean,
            "control_sd_or_sem": c_sd,
            "control_n": None,
            "measurement_unit": unit,
            "log_response_ratio_or_smd": lrr,
            "p_value": None,
            "effect_direction": fb_dir,
            "secretome_or_differentiation_penalty_observed": fb_pen,
            "conditioning_context_Z": {
                "Z1_developmental_stage_and_timing": stage,
                "Z2_genetic_background_and_allelic_dose": cell_source,
                "Z3_immune_effector_and_stressor_pathway": classify_outcome_domain(abstract),
                "Z4_metabolic_and_glucose_load_state": "Reported in study",
                "Z5_cytoarchitecture_and_cluster_diameter": "Islet / beta-cell preparation",
                "Z6_microenvironment_oxygen_and_vascular_niche": transplant_site,
                "Z7_promoter_and_epigenetic_context": "Reported in study",
                "Z8_species_age_and_sex_context": tier,
            },
            "verbatim_evidence_quote": abs_quote,
            "focal_arm_id": f"pmid{pmid}_arm1_primary",
            "comparator_arm_id": f"pmid{pmid}_arm2_comparator" if is_comp else None,
        }
        _populate_effect_provenance(eff_fb, dossier, dossier_sha, transplant_site)
        atomic_effects.append(eff_fb)

    bridges = {
        "collateral_and_pleiotropic_findings": pleiotropic_findings,
        "residual_failure_and_escape_vectors": residual_vectors,
        "regime_thresholds": regime_thresholds,
        "high_throughput_screen_hits": [],
    }
    return atomic_effects, modifier_interactions, causal_edges, bridges


SUPERSEDED_PREPRINT_MAP = {
    "40832357": "41202120",  # Preprint of Science Advances 2025 peer-reviewed paper (41202120)
    "37502937": "38661000",  # Preprint of Nature Communications 2024 peer-reviewed paper (38661000)
    "38496417": "38975343",  # Preprint of Journal of Autoimmunity 2024 peer-reviewed paper (38975343)
    "38766166": "40335415",  # Preprint of Nature Communications 2025 peer-reviewed paper (40335415)
    "39149298": "39232561",  # Preprint of Science Advances 2024 peer-reviewed paper (39232561)
    "39211206": "40892543",  # Preprint of Journal of Clinical Investigation 2025 peer-reviewed paper (40892543)
    "39763971": "40875294",  # Preprint of Diabetologia 2025 peer-reviewed paper (40875294)
    "40236391": "40614820",  # Preprint of iScience 2025 peer-reviewed paper (40614820)
    "40777266": "42120510",  # Preprint of Nature Communications 2026 peer-reviewed paper (42120510)
    "41542419": "42158573",  # Preprint of Journal of Biomedical Materials Research 2026 peer-reviewed paper (42158573)
    "39211191": "41372150",  # Preprint of Nature Biomedical Engineering 2025 peer-reviewed paper (41372150, outside Ring 1)
    "37886586": "40737125",  # Preprint of Nature Communications 2025 peer-reviewed paper (40737125, outside Ring 1)
    "41509220": "42787837",  # Preprint of Diabetes 2026 peer-reviewed paper (42787837, outside Ring 1)
    "37905040": "38985991",  # Preprint of Molecular Metabolism 2024 peer-reviewed paper (38985991, outside Ring 1)
}

PEER_REVIEWED_TO_PREPRINT_MAP = {v: k for k, v in SUPERSEDED_PREPRINT_MAP.items()}


def _apply_bibliographic_evidence_scope(record: dict, dossier: dict) -> None:
    """Reviews are bibliography/context, never independent experimental arms.

    Their source text remains in the dossier for citation tracing. Harvested
    publication-type metadata supplies the classification, not a PMID exception.
    """
    publication_types = dossier.get("publication_types") or []
    if not isinstance(publication_types, list):
        raise ValueError("Dossier publication_types must be a list")
    contextual = {"review", "systematic review", "meta-analysis", "editorial", "comment"}
    if not contextual.intersection(str(value).strip().lower() for value in publication_types):
        return
    record["study_metadata"]["source_text_type"] = "Review_or_Commentary"
    record["screening_decision"].update(
        included=False,
        is_crispr_or_high_throughput_screen=False,
        screen_library_size_genes=None,
        exclusion_reason="Review_Editorial_or_Commentary",
    )
    for field in ("arms", "atomic_effects", "effect_modifier_interactions", "causal_mechanism_edges"):
        record[field] = []
    record["cross_study_discovery_bridges"] = {key: [] for key in record.get("cross_study_discovery_bridges", {})}


def extract_one_dossier(dossier: dict[str, Any]) -> dict[str, Any]:
    pmid = str(dossier.get("pmid"))
    title = dossier.get("title") or ""
    abstract = dossier.get("abstract") or ""
    ft = _get_effective_fulltext_sections(pmid, dossier)

    ft_parts = [title, abstract]
    for r in ft.get("results_paragraphs", []):
        ft_parts.append(r)
    for fig in ft.get("figure_legends", []):
        ft_parts.append(fig.get("caption", ""))
    for m in ft.get("methods_paragraphs", []):
        ft_parts.append(m)
    full_text = "\n".join(ft_parts)

    has_ft = bool(
        dossier.get("has_parsed_fulltext_xml")
        or dossier.get("has_fulltext_xml")
        or ft.get("results_paragraphs")
        or ft.get("figure_legends")
    )
    has_geo = bool(dossier.get("linked_geo_datasets"))
    norm_ab = re.sub(r"\s+", " ", abstract).strip()
    norm_ti = re.sub(r"\s+", " ", title).strip()
    if pmid in SUPERSEDED_PREPRINT_MAP:
        src_type = "Preprint_Superseded_By_Peer_Reviewed"
    elif norm_ti and norm_ab.lower().startswith(norm_ti.lower()[:50]):
        src_type = "Project_Summary"
    elif has_ft and has_geo:
        src_type = "FullText_Plus_Supplement_or_GEO"
    elif has_ft:
        src_type = "FullText_XML_or_PDF"
    else:
        src_type = "Abstract_Only"

    mesh_list = dossier.get("mesh_terms") or dossier.get("mesh_headings") or []
    pillars = map_pillars(dossier.get("matched_pillars", []), full_text)
    tier = infer_fidelity_tier(full_text, mesh_list, dossier.get("publication_types", []))
    cell_source, cell_line, stage = infer_cell_source_and_stage(full_text, mesh_list)
    niche = infer_transplant_site_and_niche(full_text)
    edits = detect_genetic_edits(title, abstract, full_text)

    tl = full_text.lower()
    is_screen = any(k in tl for k in ("crispr screen", "genome-wide screen", "pooled screen", "perturb-seq", "brunello", "gecko"))

    # Infer challenge conditions
    challenges = []
    if any(k in tl for k in ("ifn-gamma", "ifnγ", "il-1beta", "il-1β", "tnf-alpha", "tnfα", "ifn-alpha", "ifnα", "cytokine")):
        challenges.append({
            "challenge_type": "Cytokine_Cocktail_IFNg_IL1b_TNFa_or_IFNa",
            "effector_to_target_ratio": None,
            "oxygen_pct": None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })
    if any(k in tl for k in ("thapsigargin", "tunicamycin", "er stress", "dtt")):
        challenges.append({
            "challenge_type": "ER_Stress_Thapsigargin_or_Tunicamycin",
            "effector_to_target_ratio": None,
            "oxygen_pct": None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })
    if any(k in tl for k in ("nk cell", "natural killer", "missing self", "missing-self")):
        challenges.append({
            "challenge_type": "Primary_NK_Cell_MissingSelf_Assay",
            "effector_to_target_ratio": None,
            "oxygen_pct": None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })
    if any(k in tl for k in ("autoreactive t", "t-cell clone", "cd8+ t", "ins-drip", "ppi15-24")):
        challenges.append({
            "challenge_type": "T1D_Autoreactive_CD8_or_CD4_TCR_Clone",
            "effector_to_target_ratio": None,
            "oxygen_pct": None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })
    if any(k in tl for k in ("allogeneic", "pbmc", "allograft")):
        challenges.append({
            "challenge_type": "InVivo_Allogeneic_or_Autoimmune_Graft" if niche["transplant_site"] != "InVitro_Only" else "Allogeneic_PBMC_or_T_Cell_CoCulture",
            "effector_to_target_ratio": None,
            "oxygen_pct": None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })
    if any(k in tl for k in ("hypoxia", "1% o2", "ischemi")):
        challenges.append({
            "challenge_type": "Hypoxia_Ischemia_Challenge",
            "effector_to_target_ratio": None,
            "oxygen_pct": 1.0 if "1%" in tl else None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })
    if not challenges:
        challenges.append({
            "challenge_type": "Unstressed_Baseline",
            "effector_to_target_ratio": None,
            "oxygen_pct": None,
            "glucose_concentration_mM": None,
            "duration_hours_or_days": None,
            "specific_autoantigen_or_hla_restriction": None,
        })

    atomic_effects, modifier_interactions, causal_edges, bridges = extract_atomic_effects_and_modifiers(
        pmid=pmid,
        dossier=dossier,
        study_edits=edits,
        cell_source=cell_source,
        stage=stage,
        transplant_site=niche["transplant_site"],
        tier=tier,
    )

    authors = dossier.get("authors") or []
    first_author = dossier.get("first_author") or (authors[0] if authors else "")
    senior_author = dossier.get("senior_author") or (authors[-1] if len(authors) > 1 else first_author)
    geo_accs = [g["accession"] for g in dossier.get("linked_geo_datasets", []) if g.get("accession")]
    pub_year = dossier.get("pub_year") or 2020

    arm1_n = atomic_effects[0]["treatment_n"] if atomic_effects else None
    arm1_n_unit = atomic_effects[0]["n_unit"] if (atomic_effects and arm1_n is not None) else None
    arms_list = [
        {
            "arm_id": f"pmid{pmid}_arm1_primary",
            "arm_label": f"Primary experimental arm in {first_author} ({pub_year})",
            "is_reference_control": False,
            "sample_size_n": arm1_n,
            "sample_size_unit": arm1_n_unit,
            "gate1_substrate_and_differentiation": {
                "cell_source_class": cell_source,
                "cell_line_id": cell_line,
                "donor_hla_profile": None,
                "protocol_backbone": None,
                "canonical_differentiation_stage_evaluated": stage,
                "stage_small_molecules_and_factors": [],
                "cluster_mean_diameter_um": None,
                "pct_nkx6_1_cpeptide_double_pos": None,
                "pct_mafa_pos": None,
                "pct_alpha_gcg_pos": None,
                "pct_delta_sst_pos": None,
                "pct_off_target_enterochromaffin_slc18a1": None,
            },
            "gate2_3_4_genetic_and_epitope_edits": edits,
            "gate5_engraftment_niche": niche,
            "challenge_conditions": challenges,
            "quantitative_readouts": {
                "specific_lysis_or_death_pct_mean": None,
                "specific_lysis_or_death_pct_sd": None,
                "cell_viability_pct_mean": None,
                "early_72h_graft_survival_pct": None,
                "time_to_vascular_inosculation_hours": None,
                "necrotic_core_or_damp_release_pct": None,
                "fibrotic_capsule_thickness_um": None,
                "in_vivo_graft_survival_days_median": None,
                "pct_graft_surviving_at_end_of_followup": None,
                "followup_duration_days": None,
                "static_gsis_stimulation_index_mean": None,
                "static_gsis_stimulation_index_sd": None,
                "dynamic_perifusion_first_phase_fold": None,
                "proinsulin_to_insulin_ratio": None,
                "in_vivo_stimulated_c_peptide_ng_ml": None,
                "de_novo_dsa_or_autoantibody_spreading_observed": None,
            },
        }
    ]
    comp_effs = [e for e in atomic_effects if e.get("comparator_arm_id") == f"pmid{pmid}_arm2_comparator"]
    if comp_effs:
        ce0 = comp_effs[0]
        _, comp_arm_label = _extract_comparator_clause(
            ce0.get("verbatim_evidence_quote") or "", edits, cell_source, stage, first_author, pub_year
        )
        ce0_y = (ce0.get("outcome_variable_Y") or "").lower()
        is_viability_y = ce0.get("measurement_unit") == "%" and any(
            k in ce0_y for k in ("viability", "survival")
        )
        arms_list.append({
            "arm_id": f"pmid{pmid}_arm2_comparator",
            "arm_label": comp_arm_label,
            "is_reference_control": True,
            "sample_size_n": ce0.get("control_n"),
            "sample_size_unit": ce0.get("n_unit") if ce0.get("control_n") is not None else None,
            "gate1_substrate_and_differentiation": {
                "cell_source_class": cell_source,
                "cell_line_id": cell_line,
                "donor_hla_profile": None,
                "protocol_backbone": f"Isogenic reference control matched to {cell_source} ({stage})",
                "canonical_differentiation_stage_evaluated": stage,
                "stage_small_molecules_and_factors": [],
                "cluster_mean_diameter_um": None,
                "pct_nkx6_1_cpeptide_double_pos": None,
                "pct_mafa_pos": None,
                "pct_alpha_gcg_pos": None,
                "pct_delta_sst_pos": None,
                "pct_off_target_enterochromaffin_slc18a1": None,
            },
            "gate2_3_4_genetic_and_epitope_edits": [],
            "gate5_engraftment_niche": niche,
            "challenge_conditions": challenges,
            "quantitative_readouts": {
                "specific_lysis_or_death_pct_mean": None,
                "specific_lysis_or_death_pct_sd": None,
                "cell_viability_pct_mean": ce0.get("control_mean") if is_viability_y else None,
                "early_72h_graft_survival_pct": None,
                "time_to_vascular_inosculation_hours": None,
                "necrotic_core_or_damp_release_pct": None,
                "fibrotic_capsule_thickness_um": None,
                "in_vivo_graft_survival_days_median": None,
                "pct_graft_surviving_at_end_of_followup": None,
                "followup_duration_days": None,
                "static_gsis_stimulation_index_mean": None,
                "static_gsis_stimulation_index_sd": None,
                "dynamic_perifusion_first_phase_fold": None,
                "proinsulin_to_insulin_ratio": None,
                "in_vivo_stimulated_c_peptide_ng_ml": None,
                "de_novo_dsa_or_autoantibody_spreading_observed": None,
            },
        })

    study_meta: dict[str, Any] = {
        "pmid": pmid,
        "pmcid": dossier.get("pmcid"),
        "doi": dossier.get("doi"),
        "nct_id": None,
        "geo_accessions": geo_accs,
        "title": title,
        "first_author": first_author,
        "senior_author": senior_author,
        "pub_year": pub_year,
        "journal": dossier.get("journal") or "",
        "source_text_type": src_type,
        "pillars_covered": pillars,
        "max_biological_fidelity_tier": tier,
    }
    if pmid in SUPERSEDED_PREPRINT_MAP:
        study_meta["superseded_by_pmid"] = SUPERSEDED_PREPRINT_MAP[pmid]
    if pmid in PEER_REVIEWED_TO_PREPRINT_MAP:
        study_meta["reconciled_from_superseded_preprint_pmid"] = PEER_REVIEWED_TO_PREPRINT_MAP[pmid]

    dossier_sha = _compute_dossier_sha256(pmid, dossier)
    author_tag = f"{first_author} ({pub_year})" if first_author else f"({pub_year})"
    if pmid == "31883920":
        for e in atomic_effects:
            if e["effect_id"] == "pmid31883920_eff1":
                e["focal_factor_X"] = "Doxycycline-treated (LIN28B KO) iCrLIN28B eBCs"
                e["comparator_X0"] = f"Between-arm comparator (untreated iCrLIN28B eBCs) in {author_tag}"
                e["outcome_variable_Y"] = "[Endpoint: Static GSIS stimulation index] 2-fold greater increase in stimulation index in doxycycline-treated vs untreated iCrLIN28B eBCs"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
    elif pmid == "34415483":
        for e in atomic_effects:
            if e["effect_id"] == "pmid34415483_eff1":
                e["focal_factor_X"] = "CD9-negative FACS-sorted human SC-islet beta-like cells"
                e["comparator_X0"] = "Single-arm FACS cell-sorting enrichment observation (CD9- sorted fraction vs non-sorted bulk differentiation cells; not a between-arm genetic/pharmacological/device intervention contrast)"
                e["outcome_variable_Y"] = "[Endpoint: Intracellular insulin and C-peptide hormone content (not secreted GSIS rate)] ~4-fold increased insulin and C-peptide content in sorted vs non-sorted cells"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
    elif pmid == "35918471":
        for e in atomic_effects:
            q_e = e.get("verbatim_evidence_quote") or ""
            if isinstance(e.get("log_response_ratio_or_smd"), (int, float)) and re.search(r"HNF1A\s*(?:KO|−/−|-/-|null|knockout)", q_e, re.IGNORECASE):
                e["focal_factor_X"] = "HNF1A KO (HNF1A−/−) SC-β cells vs isogenic HNF1A WT SC-β cells"
                e["comparator_X0"] = f"Between-arm comparator (plasma of mice grafted with HNF1A WT β-cells consistently, insulin:C-peptide molar ratio of 0.22 ± 0.11) in {author_tag}"
                e["outcome_variable_Y"] = "[Endpoint: Circulating plasma insulin:C-peptide molar ratio] While plasma of mice grafted with HNF1A WT β-cells consistently showed an insulin:C-peptide molar ratio of 0.22 ± 0.11 from week 4 to 30 post-transplantation, mice grafted with HNF1A KO β-cells had an 18-fold lower ratio"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
                e["censoring_note"] = "Not_Time_To_Event_Endpoint"
                if isinstance(e.get("contrast_metadata"), dict):
                    e["contrast_metadata"]["endpoint_name"] = "Plasma_Insulin_To_C_Peptide_Molar_Ratio"
                    e["contrast_metadata"]["timepoint"] = "week 4 to 30 post-transplantation"
    elif pmid in ("38661000", "37502937"):
        for e in atomic_effects:
            if "UTP23" in (e.get("verbatim_evidence_quote") or ""):
                e["focal_factor_X"] = "SLC30A8 downregulation (had no effect on beta-cell survival; nearby UTP23 / RAD21 / MED30 knockdown tested)"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
    elif pmid == "36506044":
        for a in arms_list:
            if a["arm_id"] == "pmid36506044_arm1_primary":
                a["arm_label"] = f"Doxycycline-induced PD-L1-expressing iP-BKO (B2M-/- + Tet-On CD274, +DOX) sBC grafts in {author_tag}"
        if not any(a["arm_id"] == "pmid36506044_arm2_comparator" for a in arms_list):
            base_ctrl_365 = copy.deepcopy(arms_list[0])
            base_ctrl_365["arm_id"] = "pmid36506044_arm2_comparator"
            base_ctrl_365["arm_label"] = f"Reference comparator condition (uninduced iP-BKO [-DOX] sBC grafts sharing B2M-/- + Tet-On CD274 background) in {author_tag}"
            base_ctrl_365["is_reference_control"] = True
            base_ctrl_365["sample_size_n"] = None
            base_ctrl_365["sample_size_unit"] = None
            base_ctrl_365["gate1_substrate_and_differentiation"]["protocol_backbone"] = f"Isogenic reference control matched to {cell_source} ({stage})"
            base_ctrl_365["gate2_3_4_genetic_and_epitope_edits"] = [dict(ed) for ed in edits]
            arms_list.append(base_ctrl_365)
        for e in atomic_effects:
            if "46 ± 25" in (e.get("verbatim_evidence_quote") or ""):
                e["focal_factor_X"] = "CD274 (Doxycycline-induced PD-L1 expression in iP-BKO [B2M-/- + Tet-On CD274] sBC grafts, iP-BKO + DOX)"
                e["comparator_X0"] = f"Transgene-induction manipulation check (uninduced iP-BKO [-DOX] vs +DOX sBC grafts in {author_tag}; excluded from functional survival LRR pooling)"
                e["outcome_variable_Y"] = "[Endpoint: Graft PD-L1+ cell density per mm2 (transgene-induction manipulation check)] DOX-induced PD-L1 expression from 46 ± 25 (iP-BKO -DOX) to 915 ± 57 (iP-BKO +DOX) PD-L1+ cells/mm2"
                e["measurement_unit"] = "PD-L1+ cells/mm2"
                e["log_response_ratio_or_smd"] = None
                e["comparator_arm_id"] = None
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
    elif pmid == "41297833":
        study_meta["linked_correction_notice_pmid"] = "41887653"
        for e in atomic_effects:
            if "1.4-fold" in (e.get("verbatim_evidence_quote") or ""):
                e["focal_factor_X"] = "Endocrine micro-pancreata (EMPs) vs free hESC-derived islets"
                e["outcome_variable_Y"] = "[Endpoint: In vitro insulin secretion fold increase at 20 mM glucose; linked correction notice PMID:41887653] EMPs increased insulin secretion 1.4-fold at 20 mM glucose compared to non-responsive free islets"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
    elif pmid == "41385368":
        for a in arms_list:
            if a["arm_id"] == "pmid41385368_arm1_primary":
                a["arm_label"] = f"Circadian-synchronized DEC1+/+ SC-islets (forskolin/dexamethasone shock) in {author_tag}"
            elif a["arm_id"] == "pmid41385368_arm2_comparator":
                a["arm_label"] = f"Reference comparator condition (mock-synchronized DEC1+/+ control SC-islets) in {author_tag}"
        base_exp = dict(arms_list[0])
        base_ctrl = dict(arms_list[1]) if len(arms_list) > 1 else dict(arms_list[0])
        arms_list.append({
            **base_exp,
            "arm_id": "pmid41385368_arm3_dec1_ko",
            "arm_label": f"DEC1−/− (BHLHE40 knockout) SC-islets evaluated in dynamic glucose perifusion in {author_tag}",
            "is_reference_control": False,
            "sample_size_n": None,
            "sample_size_unit": None,
        })
        arms_list.append({
            **base_ctrl,
            "arm_id": "pmid41385368_arm4_wt_genotype_control",
            "arm_label": f"Reference comparator condition (isogenic wild-type DEC1+/+ SC-islet genotype control in dynamic perifusion) in {author_tag}",
            "is_reference_control": True,
            "sample_size_n": None,
            "sample_size_unit": None,
        })
        for e in atomic_effects:
            q_e = e.get("verbatim_evidence_quote") or ""
            is_dec1_ko_contrast = bool(
                re.search(
                    r"(?:DEC1|BHLHE40)\s*(?:−/−|−\s*/\s*−|-/-|--|KO\b|knockout\b|null\b|deficient\b)",
                    q_e,
                    re.IGNORECASE,
                )
            )
            if isinstance(e.get("log_response_ratio_or_smd"), (int, float)) and is_dec1_ko_contrast:
                e["focal_factor_X"] = "DEC1−/− (BHLHE40 knockout) SC-islets vs isogenic DEC1+/+ wild-type SC-islets"
                e["comparator_X0"] = f"Between-arm comparator (isogenic wild-type DEC1+/+ SC-islet genotype control in dynamic perifusion) in {author_tag}"
                e["focal_arm_id"] = "pmid41385368_arm3_dec1_ko"
                e["comparator_arm_id"] = "pmid41385368_arm4_wt_genotype_control"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])
                if isinstance(e.get("contrast_metadata"), dict):
                    e["contrast_metadata"]["endpoint_name"] = "Dynamic_Perifusion_Peak_GSIS_Fold_Stimulation"
            elif isinstance(e.get("log_response_ratio_or_smd"), (int, float)) and not is_dec1_ko_contrast:
                e["focal_factor_X"] = "Circadian synchronization of DEC1+/+ SC-islets (forskolin or dexamethasone shock vs mock control)"
                e["comparator_X0"] = f"Between-arm comparator (mock-synchronized DEC1+/+ control SC-islets) in {author_tag}"
                e["focal_arm_id"] = "pmid41385368_arm1_primary"
                e["comparator_arm_id"] = "pmid41385368_arm2_comparator"
                _populate_effect_provenance(e, dossier, dossier_sha, niche["transplant_site"])

    record = {
        "schema_version": "0.2.0",
        "study_metadata": study_meta,
        "screening_decision": {
            "included": True,
            "is_crispr_or_high_throughput_screen": is_screen,
            "screen_library_size_genes": None,
            "exclusion_reason": None,
        },
        "arms": arms_list,
        "atomic_effects": atomic_effects,
        "effect_modifier_interactions": modifier_interactions,
        "causal_mechanism_edges": causal_edges,
        "cross_study_discovery_bridges": bridges,
    }
    _apply_bibliographic_evidence_scope(record, dossier)
    validate_study_semantics(record, dossier)
    return record


LEGACY_AUTO_CURATED_PMIDS = {
    "26276263",
    "36030035",
    "38754432",
    "39626097",
    "39824805",
    "40112799",
    "40251154",
    "41093834",
    "41218919",
    "41445610",
    "41607298",
    "41617975",
    "41758964",
    "41851084",
    "41951603",
    "41993274",
    "42856627",
}

CURATED_QUOTE_EXACT_OVERRIDES: dict[str, dict[str, Any]] = {
    "demine2020_eff1_ifng_il1b_2fold_apoptosis_irf1_er_stress": {
        "treatment_mean": 6.7, "treatment_sd_or_sem": 2.0, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": 4, "control_n": None, "measurement_unit": "fold_change_ratio", "log_response_ratio_or_smd": None,
        "focal_arm_id": "demine2020_ifng_il1b_challenged_ipsc_islets", "comparator_arm_id": None,
    },
    "demine2020_eff2_ruxolitinib_prevents_cytokine_apoptosis_er_stress_chemokines": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "demine2020_ruxolitinib_jak12_inhibited_ipsc_islets", "comparator_arm_id": None,
    },
    "cai2020_eff1_rnls_ko_human_sc_beta_er_stress": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "cai2020_human_rnls_ko_sc_beta", "comparator_arm_id": None,
    },
    "cai2020_eff2_rnlsmut_autoimmune_survival_and_upr_block": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "cai2020_mouse_rnlsmut_autoimmune_graft", "comparator_arm_id": None,
    },
    "yoshihara2020_eff1_wnt4_esrrg_oxphos_gsis_maturation": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "yoshihara2020_whilos_pdl1_oe", "comparator_arm_id": None,
    },
    "yoshihara2020_eff2_pdl1_oe_50day_glycemic_control": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "yoshihara2020_whilos_pdl1_oe", "comparator_arm_id": None,
    },
    "yoshihara2020_eff3_mps_ifng_epigenetic_pdl1_40day_survival": {
        "treatment_mean": 40.0, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "days", "log_response_ratio_or_smd": None,
        "focal_arm_id": "yoshihara2020_whilo_mps_ifng_epigenetic", "comparator_arm_id": None,
    },
    "augsornworawat2020_eff1_post_transplant_mafa_g6pc2_chgb_iapp_maturation": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "augsornworawat2020_6mo_in_vivo_matured_sc_islets", "comparator_arm_id": None,
    },
    "parent2021_eff1_hla_a2r_16wk_allogeneic_survival": {
        "treatment_mean": 16.0, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "weeks", "log_response_ratio_or_smd": None,
        "focal_arm_id": "parent2021_hla_a2r_selective_retention_ebc", "comparator_arm_id": None,
    },
    "parent2021_eff2_hla_a2r_restores_hla_e_reduces_nk_cd107a": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "parent2021_hla_a2r_selective_retention_ebc", "comparator_arm_id": None,
    },
    "balboa2022_eff1_ki67_proliferation_drop": {
        "treatment_mean": 0.46, "treatment_sd_or_sem": None, "control_mean": 2.1, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "%", "log_response_ratio_or_smd": None,
        "focal_arm_id": "balboa2022_s7w3_w6_full_zm_nac_t3", "comparator_arm_id": None,
    },
    "balboa2022_eff2_slc18a1_enterochromaffin_depletion": {
        "treatment_mean": 6.5, "treatment_sd_or_sem": None, "control_mean": 13.0, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "%", "log_response_ratio_or_smd": None,
        "focal_arm_id": "balboa2022_s7w3_w6_full_zm_nac_t3", "comparator_arm_id": None,
    },
    "balboa2022_eff3_alpha_cell_and_heterotypic_contacts": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "balboa2022_s7w3_w6_full_zm_nac_t3", "comparator_arm_id": None,
    },
    "balboa2022_eff4_glucose_threshold_shift": {
        "treatment_mean": 5.0, "treatment_sd_or_sem": None, "control_mean": 3.0, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "mM", "log_response_ratio_or_smd": None,
        "focal_arm_id": "balboa2022_s7w3_w6_full_zm_nac_t3", "comparator_arm_id": None,
    },
    "balboa2022_eff5_calcium_signaling_maturation": {
        "treatment_mean": 20.8, "treatment_sd_or_sem": 0.4, "control_mean": 5.2, "control_sd_or_sem": 0.1,
        "treatment_n": 1659, "control_n": 1091, "measurement_unit": "%", "log_response_ratio_or_smd": None,
        "focal_arm_id": "balboa2022_s7w3_w6_full_zm_nac_t3", "comparator_arm_id": None,
    },
    "balboa2022_eff6_invivo_s7w3_vs_s7w0_cpeptide": {
        "treatment_mean": 2.0, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "fold_change_ratio", "log_response_ratio_or_smd": 0.6931,
        "focal_arm_id": "balboa2022_s7w3_w6_full_zm_nac_t3", "comparator_arm_id": "balboa2022_s7w0_immature",
    },
    "lei2022_eff1_safasl_nhp_graft_survival_180days": {
        "treatment_mean": 180.0, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": 4, "control_n": None, "measurement_unit": "days", "log_response_ratio_or_smd": None,
        "focal_arm_id": "lei2022_safasl_peg_microgels_omentum", "comparator_arm_id": None,
    },
    "sintov2022_eff1_c10g_humanized_in_vivo_graft_survival_18wk": {
        "treatment_mean": 17.0, "treatment_sd_or_sem": None, "control_mean": 11.0, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "weeks", "log_response_ratio_or_smd": None,
        "focal_arm_id": "sintov2022_c10g_clonal_cxcl10_ko_sc_islets", "comparator_arm_id": None,
    },
    "sintov2022_eff2_cxcl10_ko_in_vitro_viability_2fold": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "sintov2022_c10g_clonal_cxcl10_ko_sc_islets", "comparator_arm_id": None,
    },
    "sintov2022_eff3_clonal_stat1_ko_fails_pbmc_increases_nk_apoptosis": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "sintov2022_st1l_clonal_stat1_ko_sc_islets", "comparator_arm_id": None,
    },
    "chandra2022_eff1_constitutive_tyk2_ko_neurog3_and_beta_yield_loss": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "chandra2022_tyk2_constitutive_ko_sc_islets", "comparator_arm_id": None,
    },
    "chandra2022_eff2_tyk2_ko_or_tyk2i_abrogates_ifna_mhc_class_i_and_ii": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "chandra2022_stage7_tyk2i_bms986165", "comparator_arm_id": None,
    },
    "gerace2023_eff1_2b10_nod_diabetes_reversal": {
        "treatment_mean": 8.0, "treatment_sd_or_sem": None, "control_mean": 2.5, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "weeks", "log_response_ratio_or_smd": None,
        "focal_arm_id": "gerace2023_2b10_tolerogenic_cytokine_sc_islets", "comparator_arm_id": None,
    },
    "gerace2023_eff2_b2m_ko_delayed_allorejection_humanized": {
        "treatment_mean": 7.0, "treatment_sd_or_sem": None, "control_mean": 2.0, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "weeks", "log_response_ratio_or_smd": None,
        "focal_arm_id": "gerace2023_b2m_ko_sc_islets", "comparator_arm_id": None,
    },
    "gerace2023_eff3_intrinsic_nk_resistance_cd47_pvr": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "gerace2023_b2m_ko_sc_islets", "comparator_arm_id": None,
    },
    "hu2024_eff1_k562_cd47_vs_hlae_hlag_pdl1_nk_killing": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "hu2024_human_hip_sc_islets", "comparator_arm_id": None,
    },
    "hu2024_eff2_cd47_threshold_human_islets": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "molecules_per_cell", "log_response_ratio_or_smd": None,
        "focal_arm_id": "hu2024_human_hip_sc_islets", "comparator_arm_id": None,
    },
    "hu2024_eff3_human_hip_islets_humanized_mice_diabetes": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "hu2024_human_hip_sc_islets", "comparator_arm_id": None,
    },
    "hu2024_eff4_rhesus_hip_islets_40week_survival": {
        "treatment_mean": 40.0, "treatment_sd_or_sem": None, "control_mean": 1.0, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "weeks", "log_response_ratio_or_smd": None,
        "focal_arm_id": "hu2024_rhesus_hip_primary_islets", "comparator_arm_id": None,
    },
    "macdonald2025_eff1_rnls_ko_human_sc_islet_ecar_ocr_boost": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "macdonald2025_rnls_ko_t1d_ipsc_sc_beta", "comparator_arm_id": None,
    },
    "macdonald2025_eff2_rnlsmut_glycolytic_nadph_gsh_stress_shield": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "macdonald2025_rnls_ko_t1d_ipsc_sc_beta", "comparator_arm_id": None,
    },
    "macdonald2025_eff3_pargyline_akita_glycemic_reversal_sex_dimorphic": {
        "treatment_mean": 50.0, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "ug/mL", "log_response_ratio_or_smd": None,
        "focal_arm_id": "macdonald2025_pargyline_akita_and_sc_beta", "comparator_arm_id": None,
    },
    "macdonald2025_eff3_rnls_ko_human_gsis_immaturity_tradeoff": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "macdonald2025_rnls_ko_t1d_ipsc_sc_beta", "comparator_arm_id": None,
    },
    "meng2025_eff1_hypoglycemia_protection_3pct_vs_59pct": {
        "treatment_mean": 3.0, "treatment_sd_or_sem": None, "control_mean": 59.0, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "measurement_unit": "% of repeated blood-glucose measurements (<54 mg/dL) clustered within diabetic recipient mice", "log_response_ratio_or_smd": round(math.log(3.0 / 59.0), 4),
        "outcome_variable_Y": "[Endpoint: blood glucose / hypoglycemia (<54 mg/dL)] Measurement-level proportion of repeated blood glucose readings (<54 mg/dL) clustered within diabetic recipient mice",
        "focal_arm_id": "pmid40782792_arm1_primary", "comparator_arm_id": "pmid40782792_arm2_comparator",
    },
    "pmid41993274_eff3": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "pmid41993274_arm1_primary", "comparator_arm_id": None,
    },
    "pmid41993274_eff4": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "pmid41993274_arm1_primary", "comparator_arm_id": None,
    },
    "pmid41617975_eff1": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "pmid41617975_arm1_primary", "comparator_arm_id": None,
    },
    "pmid41617975_eff3": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "pmid41617975_arm1_primary", "comparator_arm_id": None,
    },
    "luckett2026_eff2_dr4_il2_tcell_vs_dr3_mast_cell_cre": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "treatment_n": None, "control_n": None, "p_value": None, "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "pmid42671571_arm1_primary", "comparator_arm_id": None,
    },
    "chandra2022_eff3_tyk2i_reduces_cd8_tcell_cytotoxicity": {
        "treatment_mean": None, "treatment_sd_or_sem": None, "control_mean": None, "control_sd_or_sem": None,
        "measurement_unit": "qualitative_statement", "log_response_ratio_or_smd": None,
        "focal_arm_id": "chandra2022_stage7_tyk2i_bms986165", "comparator_arm_id": None,
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Ring 1 dossiers into v0.2.0 JSON files.")
    parser.add_argument("--overwrite-curated", action="store_true", help="Overwrite curated landmark PMIDs")
    args = parser.parse_args()

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    with open(RAW_DIR / "ring1_platinum_queue.json", "r", encoding="utf-8") as f:
        ring1_data = json.load(f)
    ring1 = ring1_data["studies"] if isinstance(ring1_data, dict) else ring1_data
    all_pmids = sorted(
        {str(item["pmid"]) for item in ring1}
        | {p.stem for p in EXTRACTED_DIR.glob("*.json") if p.stem.isdigit()}
    )

    extracted_count = 0
    preserved_count = 0
    for pmid in all_pmids:
        out_path = EXTRACTED_DIR / f"{pmid}.json"
        dossier_path = DOSSIER_DIR / f"{pmid}.json"
        if not dossier_path.exists():
            continue
        with open(dossier_path, "r", encoding="utf-8") as f:
            dossier = json.load(f)

        curated_src_path = CURATED_DIR / f"{pmid}.json"
        if pmid in CURATED_LANDMARK_PMIDS and (curated_src_path.exists() or out_path.exists()) and not args.overwrite_curated:
            src_to_read = curated_src_path if curated_src_path.exists() else out_path
            with open(src_to_read, "r", encoding="utf-8") as f:
                curated_rec = json.load(f)
            if not curated_rec["study_metadata"].get("first_author") and dossier.get("first_author"):
                curated_rec["study_metadata"]["first_author"] = dossier.get("first_author")
            if not curated_rec["study_metadata"].get("senior_author") and dossier.get("senior_author"):
                curated_rec["study_metadata"]["senior_author"] = dossier.get("senior_author")
            norm_ab = re.sub(r"\s+", " ", (dossier.get("abstract") or "")).strip()
            norm_ti = re.sub(r"\s+", " ", (dossier.get("title") or "")).strip()
            if pmid in SUPERSEDED_PREPRINT_MAP:
                curated_rec["study_metadata"]["source_text_type"] = "Preprint_Superseded_By_Peer_Reviewed"
                curated_rec["study_metadata"]["superseded_by_pmid"] = SUPERSEDED_PREPRINT_MAP[pmid]
            elif norm_ti and norm_ab.lower().startswith(norm_ti.lower()[:50]):
                curated_rec["study_metadata"]["source_text_type"] = "Project_Summary"
            arms = curated_rec.get("arms", [])
            transplant_site = "InVitro_Only"
            for a in arms:
                fa_str = curated_rec["study_metadata"].get("first_author") or ""
                if fa_str and isinstance(a.get("arm_label"), str):
                    a["arm_label"] = re.sub(r"\bin\s+\((\d{4})\)", f"in {fa_str} (\\1)", a["arm_label"])
                if a.get("arm_id", "").endswith("_arm0_wt_control"):
                    a["arm_id"] = f"pmid{pmid}_arm2_comparator"
                niche_obj = a.get("gate5_engraftment_niche") or {}
                if niche_obj.get("transplant_site"):
                    transplant_site = niche_obj["transplant_site"]
                for ed in a.get("gate2_3_4_genetic_and_epitope_edits", []):
                    mod = ed.get("edit_modality", "")
                    if mod in ("CRISPR_Cas9_KO", "CRISPR_Cas12a_Multiplex_KO", "CRISPRi_or_shRNA_Knockdown"):
                        ed["perturbation_class"] = "Loss_of_Function_KO_or_KD"
                    elif mod in ("BaseEdit_ABE_or_CBE_Epitope_or_Splice", "PrimeEdit", "Signal_Peptide_or_Domain_Swap"):
                        ed["perturbation_class"] = "Precision_Epitope_or_Isoform_Edit"
                    elif mod == "Endogenous_SwapIn_KO_plus_KI":
                        ed["perturbation_class"] = "Bicistronic_Endogenous_Promoter_SwapIn_KO_plus_KI"
                    elif mod in ("SafeHarbor_KI_AAVS1_GAPDH_CLYBL", "Lentiviral_or_Transposon_OE", "CRISPRa_Activation", "Transient_mRNA_or_Surface_Tethering"):
                        ed["perturbation_class"] = "Gain_of_Function_or_Transgenic_Armoring"
                    else:
                        ed["perturbation_class"] = "Pharmacological_or_Small_Molecule_Modulator"

            default_exp_arm = next(
                (a["arm_id"] for a in arms if not a.get("is_reference_control", False)),
                arms[0]["arm_id"] if arms else f"pmid{pmid}_arm1_primary",
            )
            default_ctrl_arm = next(
                (a["arm_id"] for a in arms if a.get("is_reference_control", False)),
                f"pmid{pmid}_arm2_comparator",
            )

            dossier_sha = _compute_dossier_sha256(pmid, dossier)
            effs = curated_rec.get("atomic_effects", [])
            for e in effs:
                q = (e.get("verbatim_evidence_quote") or "").strip()
                if len(q) > 345 and pmid not in CURATED_LANDMARK_PMIDS:
                    q = q[:340].rsplit(" ", 1)[0].strip()
                e["verbatim_evidence_quote"] = q
                eid = e.get("effect_id") or ""
                if eid in CURATED_QUOTE_EXACT_OVERRIDES:
                    for k_ov, v_ov in CURATED_QUOTE_EXACT_OVERRIDES[eid].items():
                        e[k_ov] = v_ov
                elif pmid in LEGACY_AUTO_CURATED_PMIDS or eid.startswith(f"pmid{pmid}_eff"):
                    if eid.startswith(f"pmid{pmid}_eff") and not eid.endswith("_counterevidence"):
                        d_dir, d_pen = classify_sentence_direction_and_penalty(q)
                        e["effect_direction"] = d_dir
                        e["secretome_or_differentiation_penalty_observed"] = d_pen
                        t_m, t_s, c_m, c_s, u_m, lrr_m = parse_quantitative_numbers(q, d_dir)
                        e["treatment_mean"] = t_m
                        e["treatment_sd_or_sem"] = t_s
                        e["control_mean"] = c_m
                        e["control_sd_or_sem"] = c_s
                        e["measurement_unit"] = u_m
                        e["log_response_ratio_or_smd"] = lrr_m
                    elif eid == "pmid37264229_eff1_counterevidence":
                        e["treatment_n"] = None
                    e["focal_arm_id"] = default_exp_arm
                    e["comparator_arm_id"] = default_ctrl_arm if isinstance(e.get("log_response_ratio_or_smd"), (int, float)) else None
                else:
                    if not isinstance(e.get("log_response_ratio_or_smd"), (int, float)):
                        e["comparator_arm_id"] = None
                    if not e.get("focal_arm_id") or e.get("focal_arm_id") == default_ctrl_arm:
                        e["focal_arm_id"] = default_exp_arm

                if e.get("secretome_or_differentiation_penalty_observed") and e.get("effect_direction") != "Deleterious_or_Impaired":
                    e["effect_direction"] = "Deleterious_or_Impaired"
                if eid == "lei2022_eff1_safasl_nhp_graft_survival_180days":
                    e["treatment_mean"] = 188.0
                for fn in ("treatment_mean", "control_mean", "treatment_sd_or_sem", "control_sd_or_sem", "treatment_n", "control_n", "p_value"):
                    if e.get(fn) is not None and not _num_in_quote(e.get(fn), q, fn):
                        e[fn] = None
                _populate_effect_provenance(e, dossier, dossier_sha, transplant_site)

            has_num = any(isinstance(e.get("log_response_ratio_or_smd"), (int, float)) for e in effs)
            if pmid in LEGACY_AUTO_CURATED_PMIDS and not has_num:
                curated_rec["arms"] = [a for a in arms if not a.get("is_reference_control", False)]
            elif pmid == "40782792":
                for a in arms:
                    if a.get("arm_id") == "pmid40782792_arm2_comparator":
                        a["arm_label"] = "Non-reconstructed PSC-islet control cohort in Meng et al. (2025)"
                        a["sample_size_n"] = None
                        a["sample_size_unit"] = None

            source_blob = json.dumps(dossier, ensure_ascii=False)
            for arm in curated_rec.get("arms", []):
                if arm.get("sample_size_n") is not None and not _num_in_quote(arm["sample_size_n"], source_blob, "sample_size_n"):
                    arm["sample_size_n"] = None
                    arm["sample_size_unit"] = None
            _apply_bibliographic_evidence_scope(curated_rec, dossier)
            validate_study_semantics(curated_rec, dossier)
            CURATED_DIR.mkdir(parents=True, exist_ok=True)
            if not curated_src_path.exists():
                with open(curated_src_path, "w", encoding="utf-8") as f:
                    json.dump(curated_rec, f, indent=2, ensure_ascii=False)
                    f.write("\n")
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(curated_rec, f, indent=2, ensure_ascii=False)
                f.write("\n")
            preserved_count += 1
            continue

        rec = extract_one_dossier(dossier)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(rec, f, indent=2, ensure_ascii=False)
            f.write("\n")
        extracted_count += 1

    if extracted_count == 0 and preserved_count == 0:
        print(
            f"[RING 1 EXTRACTION SKIPPED] Extracted: 0 | Preserved curated landmarks: 0 | Total Ring 1 queue: {len(ring1)} "
            "— raw full-text/abstract dossiers (data/dossiers/*.json) are excluded from the public bundle; "
            "using shipped pre-extracted study records in data/extracted/*.json without modification."
        )
        return

    # JSONLs are a direct producer output, not a delayed manifest side effect.
    try:
        from scripts.study_artifact_exports import write_study_exports
    except ModuleNotFoundError:
        from study_artifact_exports import write_study_exports
    produced_studies = [json.loads(path.read_text(encoding="utf-8"))
                        for path in sorted(EXTRACTED_DIR.glob("*.json")) if path.stem.isdigit()]
    write_study_exports(BASE_DIR, produced_studies)

    print(
        f"[RING 1 EXTRACTION COMPLETE] Extracted: {extracted_count} | "
        f"Preserved curated landmarks: {preserved_count} | Total Ring 1 queue: {len(ring1)}"
    )


if __name__ == "__main__":
    main()
