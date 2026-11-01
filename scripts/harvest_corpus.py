#!/usr/bin/env python3
"""
OpenIslet.org (Project Islet-1) — Multi-Source Corpus Finder & Harvester (v0.3)
Discovers and deduplicates primary experimental sources across 8 Discovery Channels:
  1. PubMed 4-Pillar Expanded Core Queries (9 sub-pillars covering all 5 Gates & 12 Failure Modes)
  2. PubMed Wide-Net Sensitivity Queries (--mode wide)
  3. Europe PMC Full-Text Body Target Discovery (catches 26 target genes in Results/Tables of SC-islet papers)
  4. Preprint Server Discovery via Europe PMC (SRC:PPR — bioRxiv, medRxiv, Research Square)
  5. 1-Hop Citation-Graph Snowballing (forward citations + backward references of 44 Golden Seed PMIDs)
  6. HPAP (Human Pancreas Analysis Program) & nPOD Primary Human T1D Islet Multi-Omics Sub-Corpus
  7. Patent Literature Discovery via Europe PMC (SRC:PAT — WIPO/EPO/USPTO SC-islet & hypoimmune patents)
  8. NCBI GEO (GSE) Functional Genomics, CRISPR Screen & Single-Cell Datasets + Linked PMIDs
Includes built-in recall verification against the 44-paper Golden Benchmark Seed Set and Multi-Pillar Overlap Matrix.
"""

import argparse
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"

# 44 Verified Primary Experimental Golden Benchmark Seed PMIDs spanning all 12 Failure Modes (FM-1..FM-12) and all 5 Gates
GOLDEN_SEED_PMIDS = {
    # Gate 1 / FM-12: SC-Islet Differentiation, Maturation & Alpha/Delta Paracrine Cytoarchitecture
    "25303535": "Pagliuca/Melton 2014 Cell (First functional human SC-beta cells)",
    "27163171": "Millman/Melton 2016 Nat Commun (T1D patient iPSC-derived SC-beta cells)",
    "32094658": "Hogrebe/Millman 2020 Nat Biotechnol (Latrunculin A cytoskeleton-guided SC-islets)",
    "35241836": "Balboa/Otonkoski 2022 Nat Biotechnol (Stage 7 functional & metabolic maturation of SC-islets)",
    "32846125": "Augsornworawat/Millman 2020 Cell Rep (scRNA-seq maturation of SC-islets post-transplant)",
    "37188763": "Augsornworawat/Millman 2023 Nat Cell Biol (snMulti-omics of SC-islets & enterochromaffin off-targets)",
    "38684699": "Wu/Deng 2024 Cell Discov (Endoderm stem cell-derived islet tissue in human patient)",
    "30412052": "Balboa/Otonkoski 2018 Nat Genet (INS mutations in iPSC beta-cell development)",
    "24658304": "Rodriguez-Diaz/Caicedo 2014 Nat Med (Alpha-cell paracrine cholinergic signaling in human islets)",
    "28325479": "Chambers/Seeley 2017 Cell Metab (Role of pancreatic preproglucagon / intra-islet GLP-1 in glucose homeostasis)",
    "16461897": "Cabrera/Caicedo 2006 PNAS (Unique cytoarchitecture of human vs rodent pancreatic islets)",
    # Gate 2 & 3 / FM-1, FM-2, FM-3, FM-4: Hypoimmune Cloaking, HLA-C Retention, Autoantigens & Neoantigens
    "30853558": "Xu/Yamanaka 2019 Cell Stem Cell (HLA-A/B/CIITA KO + matched HLA-C retained iPSC superdonors)",
    "30778232": "Deuse/Schrepfer 2019 Nat Biotechnol (B2M/CIITA KO + CD47 hypoimmunogenic iPSCs)",
    "37156915": "Hu/Schrepfer 2024 Nat Biotechnol (Hypoimmune iPSC islets survive in allogeneic NHP)",
    "37043559": "Hu/Schrepfer 2023 Cell Stem Cell (Hypoimmune primary human islets evade allogeneic and autoimmune rejection)",
    "38335966": "Hu/Schrepfer 2024 Sci Transl Med (Hypoimmune pseudoislets cure diabetes in allogeneic rhesus macaques)",
    "39758817": "Schrepfer 2025 Cell Stem Cell (Hypoimmune islet transplantation survival and CD47 threshold insights)",
    "34407395": "Parent/Hebrok 2021 Cell Rep (Selective HLA deletion retaining HLA-A2 and HLA-E in human SC-beta cells)",
    "32814902": "Yoshihara/Evans 2020 Nature (WNT4 + inducible PD-L1 immune-evasive human islet-like organoids)",
    "36599351": "Gerace/Melton 2023 Cell Rep Med (Engineering human SC-islets with IL-10/TGFb/IL-2 to promote tolerance)",
    "36055241": "Sintov/Melton 2022 Cell Stem Cell (Whole-genome CRISPR screen in SC-islets identifies JAK/STAT/CXCL10 axis)",
    "29844327": "Osum/Fife 2018 J Autoimmun (IFN-gamma drives PD-L1 on islet beta cells)",
    "30269996": "Colli/Eizirik/Richardson 2018 Diabetologia (PDL1 expressed in human T1D islets via IRF1)",
    "28263308": "Kracht/Roep 2017 Nat Med (Out-of-frame INS-DRiP neoantigen targeted by T1D CD8+ T-cells)",
    "42690930": "2026 Cell Rep (Engineering the insulin signal peptide via CHGA swap to protect human beta cells from autoimmune destruction without ER stress)",
    "26912858": "Delong/Haskins 2016 Science (Pathogenic CD4 T cells in T1D recognize Hybrid Insulin Peptides / HIPs)",
    # Gate 4 / FM-5, FM-6, FM-7, FM-8, FM-9: CRISPR Screens, TYK2, RNLS, SLC30A8, GAD65, IAPP Amyloid, CXADR
    "32719542": "Cai/Kissler 2020 Nat Metab (Genome-scale in vivo CRISPR screen identifies RNLS for beta-cell protection)",
    "39988068": "Kissler 2025 (Renalase inhibition defends against acute and chronic beta-cell stress)",
    "36289205": "Chandra/Otonkoski 2022 Nat Commun (TYK2 regulates beta-cell development and responses to IFN-alpha)",
    "32444635": "Colli/Eizirik 2020 Nat Commun (Multi-omics landscape of IFN-alpha responses in human beta cells)",
    "31900242": "Demine/Eizirik 2020 Cells (Pro-inflammatory cytokines & ER stress in human iPSC-derived beta cells)",
    "24584071": "Flannick 2014 Nat Genet (Loss-of-function mutations in SLC30A8 protect against diabetes)",
    "31676859": "Dwivedi/Groop 2019 Nat Genet (Loss of ZnT8/SLC30A8 function protects against diabetes by enhanced insulin secretion)",
    "32432213": "Menegaz/Phelps 2020 Nat Metab (Pulsatile GABA secretion from cytosolic GAD65 in the human beta cell)",
    "20699420": "Keenan/King 2010 Diabetes (Joslin 50-Year Medalist residual insulin production & beta-cell survival)",
    "18753660": "Westermark 2008 NEJM (Widespread amyloid deposition in transplanted human pancreatic islets)",
    "23086037": "Masters/O'Neill 2013 (IAPP/NLRP3 inflammasome activation in diabetes)",
    "19641142": "Tanaka 2009 Diabetes (Enterovirus infection, CXCL10, and CXCR3 circuit in beta-cell failure)",
    "30074059": "Ifie/Richardson 2018 Diabetologia (Secretory granule CAR-SIV / CXADR isoform in human pancreatic beta cells)",
    "23737756": "Nogueira/Eizirik 2013 PLoS Genet (GLIS3 modulates beta-cell apoptosis via Bim splicing)",
    "19336676": "Moore/Eizirik 2009 Diabetes (PTPN2 modulates IFN-gamma-induced beta-cell apoptosis)",
    "22442301": "Xu/Shalev 2012 Diabetes (Verapamil/calcium channel blockers suppress TXNIP and prevent beta-cell loss)",
    # Gate 5 / FM-10, FM-11: Vascular Inosculation, Pre-vascularization, Oxygen Scaffolds, FasL Microgels, IBMIR
    "33135383": "Rambaughty/Laschke 2021 (Fusion of islet cells with functional microvascular blood vessels)",
    "34480863": "Aghazadeh/Nostro 2021 Sci Transl Med (Microvessels support engraftment of human islets and hESC progenitors)",
    "35559682": "Lei/Garcia 2022 Sci Adv (FasL microgels with allogeneic islets in nonhuman primates without chronic immunosuppression)",
    # v0.4 Blind-Spot & 2025-2026 Frontier Seeds (T1D DR3/DR4 Endotypes, Mast Cells, PMNs, UP421 NEJM, HLA-B, INS 3'UTR RIDD, 5-Subtype PSC-Islets, Necroptosis/NETosis)
    "42671571": "Luckett/Oram/Richardson/Gaulton 2026 Diabetologia (HLA-DR3 mast-cell cRE vs HLA-DR4 T-cell/IL2 cRE stratified T1D GWAS)",
    "26276263": "Martino/De Tata 2015 Diabetologia (Mast cells infiltrate human T1D islets; histamine caspase-independent beta-cell death)",
    "41445610": "Luckett/EXE-T1D 2025 (Genetic risk in extremely early-onset <2 yr T1D across 6,773 cases)",
    "41758964": "Bonfield/EXE-T1D 2026 (Biallelic pathogenic IL2RA variants cause neonatal monogenic autoimmune diabetes)",
    "40757665": "Carlsson/Schrepfer 2025 NEJM (First-in-human UP421 B2M-/- CIITA-/- + CD47 allogeneic islet survival with no immunosuppression)",
    "40251154": "Hu/Schrepfer 2025 Nat Commun (PMN/neutrophil cytotoxicity against hypoimmune cells averted via CD99 + CD200)",
    "39824805": "Carre/Mallone 2025 Nat Commun (IFN-alpha promotes HLA-B-restricted conventional & alternative antigen presentation in human beta cells)",
    "40112799": "van Tienhoven/Roep 2025 Cell (Protective INS 3'UTR SNP rs3842752 creates IRE1a mRNA decay motif, reducing ER stress & INS-DRiP)",
    "40782792": "Meng/Deng 2025 Cell Stem Cell (Reconstructed 5-subtype alpha/beta/delta/epsilon/gamma human PSC-islets protect against hypoglycemia in vivo)",
    "41851084": "Misra/Nostro 2026 Nat Commun (Persistent NGN3 drives enterochromaffin off-targets; patterning maximizes human SC-islet yield)",
    "42856627": "Cochrane/Hebrok 2026 (MAFA vs senescent MAFA-S64F coordinates SC-islet metabolic maturation and stress resilience)",
    "41993274": "Martin-Vazquez/Eizirik 2026 (Bat-like ISG15 expression in alpha vs beta cells protects human islets from CVB1 & IFN-alpha)",
    "41093834": "Shibue/Kulkarni 2025 Nat Commun (Genome-wide CRISPR screen identifies SEC31A as key regulator of alpha-cell survival under ER stress)",
    "41951603": "Aggarwal/Pepper 2026 Cell Death Dis (RIPK1 inhibition with Necrostatin-1 blocks RIPK1/RIPK3/pMLKL necroptosis and improves human islet graft survival)",
    "41607298": "Lei/Huang 2026 Diabetes Obes Metab (PAD4-mediated neutrophil extracellular trap formation aggravates beta-cell damage in T1D)",
    "41218919": "Oshima/Scharfmann 2026 Diabetes (CD8+ T-cell killing of human beta cells via INS15-24 triggers inflammatory paracrine bystander dysfunction)",
}

CORE_PILLARS = {
    "Pillar1A_SC_Islet_Differentiation_Maturation": (
        '(("stem cell-derived islet*"[Title/Abstract] OR "stem-cell-derived islet*"[Title/Abstract] '
        'OR "SC-islet*"[Title/Abstract] OR "SC-beta cell*"[Title/Abstract] OR "SC-beta-cell*"[Title/Abstract] '
        'OR "iPSC-derived beta"[Title/Abstract] OR "hESC-derived beta"[Title/Abstract] '
        'OR "iPSC-derived islet*"[Title/Abstract] OR "hESC-derived islet*"[Title/Abstract] '
        'OR "stem cell-derived beta"[Title/Abstract] OR "stem-cell-derived beta"[Title/Abstract] '
        'OR "endoderm stem cell-derived islet*"[Title/Abstract] OR "CiPSC-islet*"[Title/Abstract] '
        'OR (("pluripotent stem cell*"[Title/Abstract] OR "iPSC*"[Title/Abstract] '
        'OR "hESC*"[Title/Abstract] OR "human embryonic stem cell*"[Title/Abstract]) '
        'AND ("pancreatic beta"[Title/Abstract] OR "pancreatic differentiation"[Title/Abstract] '
        'OR "beta-like cell*"[Title/Abstract] OR "insulin-producing cell*"[Title/Abstract] '
        'OR "islet-like"[Title/Abstract] OR "islet organoid*"[Title/Abstract] '
        'OR "endocrine progenitor*"[Title/Abstract] '
        'OR (("beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract]) '
        'AND ("differentiation"[Title/Abstract] OR "maturation"[Title/Abstract] '
        'OR "cytoskeleton"[Title/Abstract] OR "NKX6-1"[Title/Abstract] OR "NKX6.1"[Title/Abstract] '
        'OR "MAFA"[Title/Abstract] OR "enterochromaffin"[Title/Abstract] OR "SLC18A1"[Title/Abstract] '
        'OR "latrunculin"[Title/Abstract] OR "WNT4"[Title/Abstract] OR "Stage 6"[Title/Abstract] '
        'OR "Stage 7"[Title/Abstract]))))))'
    ),
    "Pillar1B_Islet_Paracrine_Alpha_Delta_Cytoarchitecture": (
        '(("human islet*"[Title/Abstract] OR "pancreatic islet*"[Title/Abstract] '
        'OR "islet*"[Title/Abstract] OR "SC-islet*"[Title/Abstract] '
        'OR "beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract]) '
        'AND ("cytoarchitecture"[Title/Abstract] OR "islet architecture"[Title/Abstract] '
        'OR ("paracrine"[Title/Abstract] AND ("GABA"[Title/Abstract] OR "acetylcholine"[Title/Abstract] OR "somatostatin"[Title/Abstract] OR "glucagon"[Title/Abstract] OR "cholinergic"[Title/Abstract] OR "GLP-1"[Title/Abstract])) '
        'OR ("pulsatile GABA"[Title/Abstract]) '
        'OR ("islet GLP-1"[Title/Abstract] OR ("pancreatic"[Title/Abstract] AND "preproglucagon"[Title/Abstract] AND "GLP-1"[Title/Abstract])) '
        'OR (("connexin-36"[Title/Abstract] OR "Cx36"[Title/Abstract] OR "GJD2"[Title/Abstract]) AND ("islet*"[Title/Abstract] OR "beta cell*"[Title/Abstract])) '
        'OR (("counter-regulat*"[Title/Abstract] OR "counterregulat*"[Title/Abstract]) AND "glucagon"[Title/Abstract] AND "islet*"[Title/Abstract])))'
    ),
    "Pillar2A_Islet_Hypoimmune_Cloaking_Autoantigens_Neoantigens": (
        '(("hypoimmun*"[Title/Abstract] OR "immune evasi*"[Title/Abstract] '
        'OR "immune cloaking"[Title/Abstract] OR "universal donor"[Title/Abstract] '
        'OR "localized immune tolerance"[Title/Abstract] '
        'OR "B2M"[Title/Abstract] OR "beta-2-microglobulin"[Title/Abstract] OR "beta2-microglobulin"[Title/Abstract] '
        'OR "CIITA"[Title/Abstract] OR "NLRC5"[Title/Abstract] '
        'OR "HLA-E"[Title/Abstract] OR "HLA-G"[Title/Abstract] OR "HLA-C"[Title/Abstract] '
        'OR ("HLA"[Title/Abstract] AND "selective"[Title/Abstract] AND "deletion"[Title/Abstract]) '
        'OR "CD47"[Title/Abstract] OR "SIRPalpha"[Title/Abstract] OR "SIRP-alpha"[Title/Abstract] '
        'OR "PD-L1"[Title/Abstract] OR "PDL1"[Title/Abstract] OR "CD274"[Title/Abstract] '
        'OR "FASL"[Title/Abstract] OR "Fas ligand"[Title/Abstract] OR "CTLA4-Ig"[Title/Abstract] '
        'OR "CD55"[Title/Abstract] OR "CD59"[Title/Abstract] OR "CD39"[Title/Abstract] OR "ENTPD1"[Title/Abstract] '
        'OR "hybrid insulin peptide*"[Title/Abstract] OR "peptide fusion"[Title/Abstract] '
        'OR "defective ribosomal"[Title/Abstract] OR "INS-DRiP*"[Title/Abstract] '
        'OR "neoantigen*"[Title/Abstract] OR "neo-epitope*"[Title/Abstract] OR "neoepitope*"[Title/Abstract] '
        'OR "epitope spreading"[Title/Abstract] OR "WE14"[Title/Abstract] '
        'OR (("GAD65"[Title/Abstract] OR "GAD2"[Title/Abstract] OR "ZnT8"[Title/Abstract] '
        'OR "IA-2"[Title/Abstract] OR "PTPRN"[Title/Abstract] OR "IGRP"[Title/Abstract] '
        'OR "preproinsulin"[Title/Abstract] OR "TSPAN7"[Title/Abstract]) '
        'AND ("epitope*"[Title/Abstract] OR "knockout"[Title/Abstract] OR "editing"[Title/Abstract] '
        'OR "autoreactive CD8*"[Title/Abstract] OR "ADCC"[Title/Abstract] OR "exocytosis"[Title/Abstract]))) '
        'AND ("islet*"[Title/Abstract] OR "pseudoislet*"[Title/Abstract] OR "beta cell*"[Title/Abstract] '
        'OR "beta-cell*"[Title/Abstract] OR "SC-islet*"[Title/Abstract] OR "EndoC-betaH1"[Title/Abstract]))'
    ),
    "Pillar2B_Foundational_StemCell_Hypoimmune_HLA_Engineering": (
        '(("pluripotent stem cell*"[Title/Abstract] OR "iPSC*"[Title/Abstract] OR "hESC*"[Title/Abstract]) '
        'AND ("hypoimmun*"[Title/Abstract] OR "universal donor"[Title/Abstract] '
        'OR "immune cloaking"[Title/Abstract] OR "immune compatibility"[Title/Abstract] '
        'OR ("B2M"[Title/Abstract] AND ("CIITA"[Title/Abstract] OR "HLA"[Title/Abstract])) '
        'OR ("HLA-C"[Title/Abstract] AND ("CIITA"[Title/Abstract] OR "HLA-A"[Title/Abstract])) '
        'OR ("HLA-E"[Title/Abstract] AND "B2M"[Title/Abstract]) '
        'OR ("CD47"[Title/Abstract] AND "B2M"[Title/Abstract])))'
    ),
    "Pillar2C_T1D_Endotypes_DR3_DR4_MastCell_PMN_UnconventionalT": (
        '((("type 1 diabetes"[Title/Abstract] OR "autoimmune diabetes"[Title/Abstract] OR "LADA"[Title/Abstract]) '
        'AND ("HLA-DR3"[Title/Abstract] OR "HLA-DR4"[Title/Abstract] OR "DR3-DQ2"[Title/Abstract] OR "DR4-DQ8"[Title/Abstract] '
        'OR "endotype*"[Title/Abstract] OR "T1DE1"[Title/Abstract] OR "T1DE2"[Title/Abstract] OR "EXE-T1D"[Title/Abstract]) '
        'AND ("pathway*"[Title/Abstract] OR "heterogeneity"[Title/Abstract] OR "cis-regulatory"[Title/Abstract] '
        'OR "mast cell*"[Title/Abstract] OR "insulitis"[Title/Abstract] OR "progression"[Title/Abstract])) '
        'OR (("mast cell*"[Title/Abstract] OR "histamine"[Title/Abstract] OR "NETosis"[Title/Abstract] '
        'OR "neutrophil extracellular trap*"[Title/Abstract] OR "PAD4"[Title/Abstract] OR "PADI4"[Title/Abstract] '
        'OR "polymorphonuclear"[Title/Abstract] OR "BTN3A1"[Title/Abstract] OR "MAIT cell*"[Title/Abstract] '
        'OR "CD99"[Title/Abstract] OR "CD200"[Title/Abstract]) '
        'AND ("islet*"[Title/Abstract] OR "beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] '
        'OR "hypoimmun*"[Title/Abstract] OR "insulitis"[Title/Abstract])))'
    ),
    "Pillar3A_Intrinsic_BetaCell_CRISPR_Cytokine_ER_Stress_Targets": (
        '(("beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] '
        'OR "pancreatic islet*"[Title/Abstract] OR "human islet*"[Title/Abstract] '
        'OR "islet*"[Title/Abstract] OR "SC-islet*"[Title/Abstract] OR "EndoC-betaH1"[Title/Abstract]) '
        'AND ("CRISPR"[Title/Abstract] OR "genetic screen*"[Title/Abstract] '
        'OR "genome-wide screen*"[Title/Abstract] OR "genome-scale"[Title/Abstract] '
        'OR "gene editing"[Title/Abstract] OR "base editing"[Title/Abstract] '
        'OR "prime editing"[Title/Abstract] OR "Perturb-seq"[Title/Abstract] '
        'OR "RNLS"[Title/Abstract] OR "renalase"[Title/Abstract] '
        'OR "TYK2"[Title/Abstract] OR "TET2"[Title/Abstract] OR "TNFAIP3"[Title/Abstract] '
        'OR ("A20"[Title/Abstract] AND ("NF-kappaB"[Title/Abstract] OR "NF-kB"[Title/Abstract] OR "apoptosis"[Title/Abstract] OR "transplant*"[Title/Abstract])) '
        'OR "SERPINB9"[Title/Abstract] OR "proteinase inhibitor 9"[Title/Abstract] OR "PI-9"[Title/Abstract] '
        'OR "GLIS3"[Title/Abstract] OR "PTPN2"[Title/Abstract] OR "TXNIP"[Title/Abstract] '
        'OR ("verapamil"[Title/Abstract] AND "diabetes"[Title/Abstract]) '
        'OR ("calcium channel blocker*"[Title/Abstract] AND "beta-cell*"[Title/Abstract]) '
        'OR "MANF"[Title/Abstract] OR "XBP1"[Title/Abstract] '
        'OR "SLC30A8"[Title/Abstract] OR ("ZnT8"[Title/Abstract] AND ("loss of ZnT8"[Title/Abstract] OR "loss-of-function"[Title/Abstract] OR "knockout"[Title/Abstract] OR "secretion"[Title/Abstract] OR "hypoxia"[Title/Abstract])) '
        'OR "G6PC2"[Title/Abstract] OR "TGM2"[Title/Abstract] OR "transglutaminase 2"[Title/Abstract] '
        'OR "IFIH1"[Title/Abstract] OR "MDA5"[Title/Abstract]) '
        'AND ("diabetes"[Title/Abstract] OR "autoimmun*"[Title/Abstract] '
        'OR "cytokine*"[Title/Abstract] OR "interferon*"[Title/Abstract] '
        'OR "ER stress"[Title/Abstract] OR "hypoxia"[Title/Abstract] '
        'OR "survival"[Title/Abstract] OR "apoptosis"[Title/Abstract] OR "islet*"[Title/Abstract]))'
    ),
    "Pillar3B_Amyloid_Enterovirus_Senescence_Medalists": (
        '((("beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] '
        'OR "pancreatic islet*"[Title/Abstract] OR "human islet*"[Title/Abstract] '
        'OR "SC-islet*"[Title/Abstract] OR "islet*"[Title/Abstract] OR "diabetes"[Title/Abstract]) '
        'AND ("CXADR"[Title/Abstract] OR "CAR-SIV"[Title/Abstract] '
        'OR "Coxsackie and adenovirus receptor"[Title/Abstract] '
        'OR ("coxsackie*"[Title/Abstract] AND ("beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] OR "islet*"[Title/Abstract] OR "tropism"[Title/Abstract] OR "persistent"[Title/Abstract] OR "granule*"[Title/Abstract])) '
        'OR ("enterovir*"[Title/Abstract] AND ("CXCL10"[Title/Abstract] OR "CXADR"[Title/Abstract] OR "CAR-SIV"[Title/Abstract] OR "beta-cell failure"[Title/Abstract] OR "secretory granule*"[Title/Abstract])) '
        'OR ("CXCL10"[Title/Abstract] AND "CXCR3"[Title/Abstract] AND ("islet*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] OR "beta cell*"[Title/Abstract])) '
        'OR "senolytic*"[Title/Abstract] OR ("p16INK4a"[Title/Abstract] AND ("islet*"[Title/Abstract] OR "beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract])) '
        'OR ("Medalist*"[Title/Abstract] AND ("Joslin"[Title/Abstract] OR "50 year*"[Title/Abstract] OR "residual insulin"[Title/Abstract])) '
        'OR (("IAPP"[Title/Abstract] OR "islet amyloid polypeptide"[Title/Abstract] OR "islet amyloid"[Title/Abstract] OR "pramlintide"[Title/Abstract]) '
        'AND ("transplant*"[Title/Abstract] OR "graft*"[Title/Abstract] OR "NLRP3"[Title/Abstract] OR "inflammasome"[Title/Abstract] OR "proline"[Title/Abstract])) '
        'OR ("amyloid*"[Title/Abstract] AND "transplanted"[Title/Abstract] AND "islet*"[Title/Abstract]) '
        'OR ("NLRP3"[Title/Abstract] AND "inflammasome"[Title/Abstract] AND ("IAPP"[Title/Abstract] OR "islet*"[Title/Abstract] OR "diabetes"[Title/Abstract])))))'
    ),
    "Pillar3C_HPAP_nPOD_Human_T1D_Islet_Atlas": (
        '(("Human Pancreas Analysis Program"[Title/Abstract] OR "HPAP"[Title/Abstract] '
        'OR "PANC-DB"[Title/Abstract] '
        'OR "Network for Pancreatic Organ Donors with Diabetes"[Title/Abstract] '
        'OR "nPOD"[Title/Abstract]) '
        'AND ("islet*"[Title/Abstract] OR "beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] '
        'OR "pancrea*"[Title/Abstract] OR "type 1 diabetes"[Title/Abstract] OR "autoimmun*"[Title/Abstract]))'
    ),
    "Pillar4A_SC_Islet_and_Allogeneic_Engraftment_Devices": (
        '(("stem cell-derived islet*"[Title/Abstract] OR "stem-cell-derived islet*"[Title/Abstract] '
        'OR "SC-islet*"[Title/Abstract] OR "allogeneic islet*"[Title/Abstract] '
        'OR "islet allograft*"[Title/Abstract] OR "VX-880"[Title/Abstract] '
        'OR "zimislecel"[Title/Abstract] OR "VX-264"[Title/Abstract] '
        'OR "PEC-Direct"[Title/Abstract] OR "PEC-Encap"[Title/Abstract] '
        'OR "VC-01"[Title/Abstract] OR "VC-02"[Title/Abstract] '
        'OR "ViaCyte"[Title/Abstract] OR "UP421"[Title/Abstract] '
        'OR "VCTX210"[Title/Abstract] OR "VCTX211"[Title/Abstract] OR "Cell Pouch"[Title/Abstract]) '
        'AND ("engraftment"[Title/Abstract] OR "vasculariz*"[Title/Abstract] '
        'OR "revasculariz*"[Title/Abstract] OR "inosculation"[Title/Abstract] '
        'OR "encapsulat*"[Title/Abstract] OR "hypoxia"[Title/Abstract] '
        'OR "oxygen*"[Title/Abstract] OR "hydrogel"[Title/Abstract] '
        'OR "scaffold"[Title/Abstract] OR "subcutaneous"[Title/Abstract] '
        'OR "omental"[Title/Abstract] OR "intramuscular"[Title/Abstract] '
        'OR "IBMIR"[Title/Abstract] OR "tissue factor"[Title/Abstract] '
        'OR "FasL"[Title/Abstract] OR "microgel*"[Title/Abstract] '
        'OR "survival"[Title/Abstract] OR "clinical trial"[Title/Abstract] OR "C-peptide"[Title/Abstract]))'
    ),
    "Pillar4B_Islet_Vascularization_Inosculation_Oxygen_IBMIR_SafetySwitch": (
        '(("islet transplant*"[Title/Abstract] OR "transplanted islet*"[Title/Abstract] '
        'OR "islet graft*"[Title/Abstract] OR "islet organoid*"[Title/Abstract] '
        'OR "human islet*"[Title/Abstract] '
        'OR ("islet cells"[Title/Abstract] AND "transplantation"[Title/Abstract])) '
        'AND ("inosculation"[Title/Abstract] OR "prevasculariz*"[Title/Abstract] '
        'OR "pre-vasculariz*"[Title/Abstract] OR "microvascular fragment*"[Title/Abstract] '
        'OR "microvessel*"[Title/Abstract] '
        'OR "functional blood vessel*"[Title/Abstract] OR "endothelial co-cultur*"[Title/Abstract] '
        'OR ("endothelial cell*"[Title/Abstract] AND "pericyte*"[Title/Abstract]) '
        'OR "parathyroid"[Title/Abstract] OR "oxygen-generating"[Title/Abstract] '
        'OR "oxygen-releasing"[Title/Abstract] OR "calcium peroxide"[Title/Abstract] '
        'OR "CaO2"[Title/Abstract] OR "perfluorocarbon"[Title/Abstract] '
        'OR "IBMIR"[Title/Abstract] OR "instant blood-mediated"[Title/Abstract] '
        'OR ("tissue factor"[Title/Abstract] AND ("clotting"[Title/Abstract] OR "thrombo*"[Title/Abstract] OR "coagulation"[Title/Abstract] OR "IBMIR"[Title/Abstract])) '
        'OR "omental pouch"[Title/Abstract] OR "zwitterionic"[Title/Abstract] '
        'OR "iCasp9"[Title/Abstract] OR "inducible caspase"[Title/Abstract] OR "suicide gene"[Title/Abstract]))'
    ),
}

WIDE_QUERIES = {
    "Wide_All_StemCell_Beta_Islet": (
        '(("stem cell*"[Title/Abstract] OR "iPSC*"[Title/Abstract] '
        'OR "hESC*"[Title/Abstract] OR "pluripotent"[Title/Abstract]) '
        'AND ("beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] '
        'OR "pancreatic islet*"[Title/Abstract] OR "insulin-producing"[Title/Abstract] '
        'OR "SC-islet*"[Title/Abstract]))'
    ),
    "Wide_All_Islet_Transplant_Rejection_Survival_Vascular": (
        '(("islet transplant*"[Title/Abstract] OR "islet graft*"[Title/Abstract] '
        'OR "allogeneic islet*"[Title/Abstract] OR "transplanted islet*"[Title/Abstract]) '
        'AND ("rejection"[Title/Abstract] OR "survival"[Title/Abstract] '
        'OR "tolerance"[Title/Abstract] OR "engraftment"[Title/Abstract] '
        'OR "encapsulat*"[Title/Abstract] OR "hypoxia"[Title/Abstract] '
        'OR "vasculariz*"[Title/Abstract] OR "revasculariz*"[Title/Abstract] '
        'OR "oxygen*"[Title/Abstract] OR "IBMIR"[Title/Abstract] '
        'OR "immunosuppress*"[Title/Abstract]))'
    ),
    "Wide_All_T1D_BetaCell_Gene_Perturbation_Protection": (
        '(("type 1 diabetes"[Title/Abstract] OR "autoimmune diabetes"[Title/Abstract] '
        'OR "LADA"[Title/Abstract]) '
        'AND ("beta cell*"[Title/Abstract] OR "beta-cell*"[Title/Abstract] OR "islet*"[Title/Abstract]) '
        'AND ("protect*"[Title/Abstract] OR "survival"[Title/Abstract] '
        'OR "apoptosis"[Title/Abstract] OR "knockout"[Title/Abstract] '
        'OR "overexpress*"[Title/Abstract] OR "CRISPR"[Title/Abstract] '
        'OR "base edit*"[Title/Abstract] OR "prime edit*"[Title/Abstract] '
        'OR "inhibit*"[Title/Abstract] OR "C-peptide"[Title/Abstract]))'
    ),
}

# Europe PMC Supplementary Discovery Queries (Full-Text Body Targets & Preprints)
EUROPEPMC_SUPPLEMENTARY_QUERIES = {
    "EPMC_FullText_SC_Islet_26_Target_Genes": (
        'TITLE_ABS:("SC-islet*" OR "stem cell-derived islet*" OR "stem-cell-derived islet*" '
        'OR "SC-beta" OR "iPSC-derived beta" OR "hESC-derived beta" OR "EndoC-betaH1") '
        'AND ("TYK2" OR "RNLS" OR "IFIH1" OR "CXADR" OR "NLRC5" OR "SLC30A8" OR "G6PC2" '
        'OR "PTPRN" OR "GAD2" OR "IAPP" OR "CHGA" OR "TGM2" OR "SERPINB9" OR "TNFAIP3" '
        'OR "PTPN2" OR "TXNIP" OR "MANF" OR "GLIS3" OR "B2M" OR "CIITA" OR "CD47" '
        'OR "CD274" OR "CD55" OR "CD59" OR "ENTPD1" OR "CXCL10") AND HAS_FT:y'
    ),
    "EPMC_Preprints_SC_Islet_and_Hypoimmune": (
        'SRC:PPR AND ("SC-islet" OR "SC-islets" OR "stem cell-derived islet" '
        'OR "stem-cell-derived islet" OR "stem cell-derived islets" '
        'OR ("hypoimmune" AND "islet") OR ("hypoimmunogenic" AND "islet") '
        'OR ("B2M" AND "islet") OR ("CD47" AND "islet") OR ("HLA-E" AND "islet") '
        'OR ("VX-880" OR "zimislecel" OR "VX-264" OR "UP421"))'
    ),
}

# Europe PMC Patent Discovery Queries (SRC:PAT)
EUROPEPMC_PATENT_QUERIES = {
    "Patent_SC_Islet_and_Hypoimmune_Engineering": (
        'SRC:PAT AND (islet AND (pluripotent OR iPSC OR hESC OR stem '
        'OR hypoimmunogenic OR hypoimmune OR B2M OR CIITA OR CD47 '
        'OR HLA OR encapsulation OR vascularization))'
    ),
}

# NCBI GEO (GDS) Functional Genomics & Single-Cell Dataset Queries
GEO_GDS_QUERIES = {
    "GEO_SC_Islet_Differentiation_and_MultiOmics": (
        '("SC-islet*" OR "stem cell-derived islet*" OR "stem-cell-derived islet*" '
        'OR "SC-beta" OR "iPSC-derived beta" OR "hESC-derived beta" '
        'OR "iPSC-derived islet*" OR "hESC-derived islet*") AND "gse"[ETYP]'
    ),
    "GEO_BetaCell_CRISPR_Screens_and_Gene_Perturbation": (
        '(("CRISPR" OR "Perturb-seq" OR "Cas9" OR "Cas12a" OR "genetic screen" OR "genome-wide screen") '
        'AND ("beta cell*" OR "beta-cell*" OR "pancreatic islet*" OR "human islet*" '
        'OR "EndoC-betaH1" OR "SC-islet*" OR "NIT-1" OR "MIN6")) AND "gse"[ETYP]'
    ),
    "GEO_HPAP_and_nPOD_Human_T1D_Islet_Omics": (
        '(("Human Pancreas Analysis Program" OR "HPAP" OR "nPOD" '
        'OR "Network for Pancreatic Organ Donors") '
        'AND ("islet*" OR "pancrea*" OR "beta cell*" OR "beta-cell*")) AND "gse"[ETYP]'
    ),
    "GEO_Hypoimmune_and_Target_Gene_Islet_Studies": (
        '(("hypoimmune" OR "hypoimmunogenic" OR "B2M" OR "CIITA" OR "CD47" '
        'OR "TYK2" OR "RNLS" OR "SLC30A8" OR "CXADR" OR "TET2" OR "MANF" OR "TNFAIP3") '
        'AND ("pancreatic islet*" OR "human islet*" OR "SC-islet*" OR "beta-cell*" OR "EndoC-betaH1")) AND "gse"[ETYP]'
    ),
}

NON_REVIEW_FILTER = ' NOT ("Review"[Publication Type] OR "Editorial"[Publication Type] OR "Comment"[Publication Type])'


def esearch_pmids(query: str, retmax: int = 10000) -> list[str]:
    if type(retmax) is not int or retmax <= 0:
        raise ValueError("PubMed page size must be a positive integer")
    pmids = []
    seen = set()
    retstart = 0
    expected_total = None
    while True:
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retstart": str(retstart),
            "retmax": str(min(retmax, 5000)),
        }
        api_key = os.environ.get("NCBI_API_KEY")
        if api_key:
            params["api_key"] = api_key
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        res = data["esearchresult"]
        batch = res.get("idlist", [])
        raw_total = res.get("count")
        if not isinstance(raw_total, str) or not raw_total.isdigit():
            raise ValueError("PubMed response has an invalid result count")
        total = int(raw_total)
        if expected_total is None:
            expected_total = total
        elif total != expected_total:
            raise ValueError("PubMed result count changed during pagination; retry a stable query snapshot")
        if not isinstance(batch, list) or any(not isinstance(pid, str) or not pid.isdigit() or int(pid) <= 0 for pid in batch):
            raise ValueError("PubMed response has invalid PMID identifiers")
        if len(set(batch)) != len(batch) or seen.intersection(batch):
            raise ValueError("PubMed pagination repeated PMID identifiers; query snapshot is incomplete")
        if len(batch) > min(retmax, 5000) or retstart + len(batch) > expected_total:
            raise ValueError("PubMed response exceeds its declared result count or page size")
        if not batch and retstart < total:
            raise ValueError(f"PubMed pagination truncated at {retstart} of {total}; partition the query before claiming complete retrieval")
        pmids.extend(batch)
        seen.update(batch)
        retstart += len(batch)
        time.sleep(0.35)
        if retstart >= total:
            break
    return pmids


def _validated_europepmc_page(data: dict, total, fetched: int, seen_cursors: set[str], seen_records: set[tuple[str, str]]):
    """Distinguish a complete result set from an exhausted or repeated cursor."""
    if not isinstance(data, dict) or not isinstance(data.get("resultList"), dict):
        raise ValueError("Europe PMC response has invalid result container")
    count = data.get("hitCount")
    results = data.get("resultList", {}).get("result", [])
    if (type(count) is not int or count < 0 or not isinstance(results, list)
            or any(not isinstance(record, dict) for record in results)):
        raise ValueError("Europe PMC response has invalid result count or records")
    if total is not None and count != total:
        raise ValueError("Europe PMC result count changed during pagination")
    for record in results:
        pmid = record.get("pmid")
        if pmid not in (None, "") and (not isinstance(pmid, str) or not pmid.isascii()
                or not pmid.isdigit() or not 1 <= len(pmid) <= 9 or pmid[0] == "0"):
            raise ValueError("Europe PMC result has invalid PMID identity")
        source = record.get("source", "MED" if record.get("pmid") else "")
        identifier = record.get("id", record.get("pmid"))
        if not isinstance(source, str) or not source or not isinstance(identifier, str) or not identifier:
            raise ValueError("Europe PMC result lacks a valid source identity")
        identity = (source, identifier)
        if identity in seen_records:
            raise ValueError("Europe PMC pagination repeated a record identity")
        seen_records.add(identity)
    fetched += len(results)
    if fetched > count:
        raise ValueError("Europe PMC response exceeds its declared result count")
    if fetched == count:
        return results, count, fetched, None
    next_cursor = data.get("nextCursorMark")
    if not results or not isinstance(next_cursor, str) or not next_cursor or next_cursor in seen_cursors:
        raise ValueError("Europe PMC pagination ended before retrieval completion")
    seen_cursors.add(next_cursor)
    return results, count, fetched, next_cursor


def europepmc_search_ids(query: str, max_pages: int = 10) -> tuple[set[str], list[dict]]:
    """Queries Europe PMC REST API for additional PMIDs and Preprint (PPR/DOI) records."""
    if type(max_pages) is not int or max_pages <= 0:
        raise ValueError("Europe PMC page budget must be a positive integer")
    pmids: set[str] = set()
    preprints: list[dict] = []
    cursor = "*"
    total, fetched, seen_cursors = None, 0, {cursor}
    seen_records: set[tuple[str, str]] = set()
    for _ in range(max_pages):
        params = {
            "query": query,
            "format": "json",
            "pageSize": "1000",
            "resultType": "lite",
            "cursorMark": cursor,
        }
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        results, total, fetched, next_cursor = _validated_europepmc_page(data, total, fetched, seen_cursors, seen_records)
        for r in results:
            pub_type = str(r.get("pubType", "")).lower()
            if "review" in pub_type:
                continue
            pmid = r.get("pmid")
            if pmid:
                pmids.add(str(pmid))
            elif r.get("source") == "PPR":
                preprints.append({
                    "id": r.get("id"),
                    "doi": r.get("doi"),
                    "title": r.get("title"),
                    "pub_year": r.get("pubYear"),
                    "authorString": r.get("authorString"),
                    "journal": r.get("bookOrReportDetails", {}).get("publisher", "Preprint"),
                })
        if next_cursor is None:
            break
        cursor = next_cursor
        time.sleep(0.35)
    else:
        raise ValueError("Europe PMC page budget exhausted before retrieval completion; partition the query or increase the explicit budget")
    return pmids, preprints


def elink_citation_snowball(seed_pmids: list[str]) -> tuple[set[str], set[str]]:
    """Uses NCBI ELink to retrieve all primary PubMed articles that cite OR are cited by the 42 Golden Seeds."""
    params = {
        "dbfrom": "pubmed",
        "db": "pubmed",
        "id": ",".join(seed_pmids),
        "linkname": "pubmed_pubmed_citedin,pubmed_pubmed_refs",
        "retmode": "json",
    }
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=45) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    cited_by: set[str] = set()
    refs: set[str] = set()
    for lset in data.get("linksets", []):
        for db in lset.get("linksetdbs", []):
            lname = db.get("linkname", "")
            links = [str(x) for x in db.get("links", [])]
            if lname == "pubmed_pubmed_citedin":
                cited_by.update(links)
            elif lname == "pubmed_pubmed_refs":
                refs.update(links)
    time.sleep(0.35)
    return cited_by, refs


def europepmc_search_patents(query: str, max_pages: int = 5) -> list[dict]:
    """Queries Europe PMC REST API for patent records (SRC:PAT)."""
    if type(max_pages) is not int or max_pages <= 0:
        raise ValueError("Europe PMC page budget must be a positive integer")
    patents_by_id: dict[str, dict] = {}
    cursor = "*"
    total, fetched, seen_cursors = None, 0, {cursor}
    seen_records: set[tuple[str, str]] = set()
    for _ in range(max_pages):
        params = {
            "query": query,
            "format": "json",
            "pageSize": "1000",
            "resultType": "lite",
            "cursorMark": cursor,
        }
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        results, total, fetched, next_cursor = _validated_europepmc_page(data, total, fetched, seen_cursors, seen_records)
        for r in results:
            if r.get("source") == "PAT":
                pat_id = str(r.get("id", ""))
                if pat_id and pat_id not in patents_by_id:
                    patents_by_id[pat_id] = {
                        "id": pat_id,
                        "title": r.get("title"),
                        "pub_year": r.get("pubYear"),
                        "authorString": r.get("authorString"),
                        "patent_source": r.get("source"),
                    }
        if next_cursor is None:
            break
        cursor = next_cursor
        time.sleep(0.35)
    else:
        raise ValueError("Europe PMC page budget exhausted before retrieval completion; partition the query or increase the explicit budget")
    return list(patents_by_id.values())


def geo_gds_search_datasets(queries: dict[str, str]) -> tuple[list[dict], set[str]]:
    """Queries NCBI GEO DataSets (db=gds) for GSE series, filters out amyloid-beta/TGF-beta false positives, and extracts linked PubMed IDs."""
    gds_id_to_categories: dict[str, list[str]] = {}
    api_key = os.environ.get("NCBI_API_KEY")
    for q_name, q_str in queries.items():
        params = {
            "db": "gds",
            "term": q_str,
            "retmode": "json",
            "retmax": "2000",
        }
        if api_key:
            params["api_key"] = api_key
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        result = data.get("esearchresult", {})
        id_list = result.get("idlist", [])
        count = result.get("count")
        if (not isinstance(count, str) or not count.isascii() or not count.isdigit()
                or not isinstance(id_list, list)
                or any(not isinstance(uid, str) or not uid.isascii() or not uid.isdigit() or int(uid) <= 0 for uid in id_list)
                or len(set(id_list)) != len(id_list)):
            raise ValueError("GEO search response has invalid count or identifiers")
        if len(id_list) != int(count):
            raise ValueError("GEO search retrieval is incomplete; partition the query before claiming complete discovery")
        print(f"[GEO GDS] {q_name}: {len(id_list)} raw GSE series")
        for gid in id_list:
            gds_id_to_categories.setdefault(str(gid), []).append(q_name)
        time.sleep(0.35)

    all_gds_ids = sorted(gds_id_to_categories.keys(), reverse=True)
    datasets: list[dict] = []
    linked_pmids: set[str] = set()
    batch_size = 200
    islet_terms = (
        "islet", "pancrea", "sc-beta", "sc-islet", "beta-cell", "beta cell",
        "β-cell", "β cell", "beta-like", "β-like", "endoc-betah", "endoc-βh",
        "min6", "nit-1", "ins-1", "hpap", "npod", "insulin-producing",
        "insulin-secreting", "endocrine progenitor",
    )
    for i in range(0, len(all_gds_ids), batch_size):
        chunk = all_gds_ids[i : i + batch_size]
        params = {
            "db": "gds",
            "id": ",".join(chunk),
            "retmode": "json",
        }
        if api_key:
            params["api_key"] = api_key
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        data_bytes = urllib.parse.urlencode(params).encode("utf-8")
        req = urllib.request.Request(url, data=data_bytes, method="POST")
        with urllib.request.urlopen(req, timeout=45) as resp:
            sdata = json.loads(resp.read().decode("utf-8"))
        res_map = sdata.get("result", {})
        for gid in chunk:
            item = res_map.get(gid)
            if not isinstance(item, dict) or item.get("error"):
                raise ValueError(f"GEO metadata response is missing or rejected requested identifier {gid}")
            title_str = item.get("title", "") or ""
            summary_str = item.get("summary", "") or ""
            combined_lower = (title_str + " " + summary_str).lower()
            # Exclude NCBI GDS hyphen-split false positives (e.g. amyloid-beta, TGF-beta in brain/liver)
            if not any(term in combined_lower for term in islet_terms):
                continue
            acc = item.get("accession", "")
            pmids_for_gse = [str(p) for p in item.get("pubmedids", []) if p]
            linked_pmids.update(pmids_for_gse)
            datasets.append({
                "gds_uid": gid,
                "accession": acc,
                "title": title_str,
                "summary": summary_str[:400],
                "taxon": item.get("taxon", ""),
                "gds_type": item.get("gdstype", ""),
                "n_samples": item.get("n_samples", 0),
                "pub_date": item.get("pdat", ""),
                "pubmed_ids": pmids_for_gse,
                "matched_geo_queries": gds_id_to_categories.get(gid, []),
            })
        time.sleep(0.35)
    return datasets, linked_pmids


def compute_overlap_breakdown(pillar_map: dict[str, list[str]]) -> dict:
    """Computes multi-pillar and multi-gate overlap statistics across all discovered PMIDs."""
    subpillar_count_hist: dict[int, int] = {}
    major_pillar_count_hist: dict[int, int] = {}
    major_pairwise: dict[str, int] = {}

    def get_major_pillars(tags: list[str]) -> set[str]:
        maj: set[str] = set()
        for t in tags:
            if t.startswith("Pillar1"):
                maj.add("P1_Differentiation_Cytoarchitecture")
            elif t.startswith("Pillar2"):
                maj.add("P2_Allo_Auto_Cloaking")
            elif t.startswith("Pillar3"):
                maj.add("P3_Intrinsic_Stress_Amyloid_Viral_HPAP")
            elif t.startswith("Pillar4"):
                maj.add("P4_Engraftment_Vascularization")
        return maj

    sc_islet_multi_gate = {
        "P1A_SC_Islet_total": 0,
        "P1A_and_P2_Cloaking": 0,
        "P1A_and_P3_Intrinsic_Targets": 0,
        "P1A_and_P4_Engraftment_Vascular": 0,
        "P1A_and_2plus_other_Major_Pillars": 0,
        "All_4_Major_Pillars_Simultaneously": 0,
    }
    high_density_pmids_3plus_subpillars: list[str] = []

    major_labels = [
        "P1_Differentiation_Cytoarchitecture",
        "P2_Allo_Auto_Cloaking",
        "P3_Intrinsic_Stress_Amyloid_Viral_HPAP",
        "P4_Engraftment_Vascularization",
    ]

    for pmid, tags in pillar_map.items():
        core_tags = [t for t in tags if t.startswith("Pillar1") or t.startswith("Pillar2") or t.startswith("Pillar3") or t.startswith("Pillar4")]
        n_sub = len(core_tags)
        subpillar_count_hist[n_sub] = subpillar_count_hist.get(n_sub, 0) + 1
        if n_sub >= 3:
            high_density_pmids_3plus_subpillars.append(pmid)

        maj = get_major_pillars(tags)
        n_maj = len(maj)
        major_pillar_count_hist[n_maj] = major_pillar_count_hist.get(n_maj, 0) + 1
        if n_maj == 4:
            sc_islet_multi_gate["All_4_Major_Pillars_Simultaneously"] += 1

        for i in range(len(major_labels)):
            for j in range(i + 1, len(major_labels)):
                if major_labels[i] in maj and major_labels[j] in maj:
                    pair_key = f"{major_labels[i]}__x__{major_labels[j]}"
                    major_pairwise[pair_key] = major_pairwise.get(pair_key, 0) + 1

        if "Pillar1A_SC_Islet_Differentiation_Maturation" in tags:
            sc_islet_multi_gate["P1A_SC_Islet_total"] += 1
            if "P2_Allo_Auto_Cloaking" in maj:
                sc_islet_multi_gate["P1A_and_P2_Cloaking"] += 1
            if "P3_Intrinsic_Stress_Amyloid_Viral_HPAP" in maj:
                sc_islet_multi_gate["P1A_and_P3_Intrinsic_Targets"] += 1
            if "P4_Engraftment_Vascularization" in maj:
                sc_islet_multi_gate["P1A_and_P4_Engraftment_Vascular"] += 1
            if len(maj - {"P1_Differentiation_Cytoarchitecture"}) >= 2:
                sc_islet_multi_gate["P1A_and_2plus_other_Major_Pillars"] += 1

    return {
        "subpillar_hit_distribution": {str(k): subpillar_count_hist[k] for k in sorted(subpillar_count_hist.keys())},
        "major_pillar_hit_distribution": {str(k): major_pillar_count_hist[k] for k in sorted(major_pillar_count_hist.keys())},
        "major_pillar_pairwise_intersections": major_pairwise,
        "sc_islet_multi_gate_intersections": sc_islet_multi_gate,
        "high_density_3plus_subpillars_count": len(high_density_pmids_3plus_subpillars),
        "high_density_3plus_subpillars_sample_pmids": sorted(high_density_pmids_3plus_subpillars, reverse=True)[:50],
    }


def efetch_metadata(pmids: list[str], pillar_map: dict[str, list[str]], out_path: Path) -> int:
    if not isinstance(pmids, list) or any(not isinstance(pid, str) or not pid.isdigit() or int(pid) <= 0 for pid in pmids) or len(set(pmids)) != len(pmids):
        raise ValueError("Requested corpus PMIDs must be unique positive numeric strings")
    requested_pmids = set(pmids)
    existing_pmids: set[str] = set()
    if out_path.exists():
        with open(out_path, "r", encoding="utf-8") as fin:
            for line_number, line in enumerate(fin, 1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                    pmid = record.get("pmid") if isinstance(record, dict) else None
                    if not isinstance(pmid, str) or not pmid.isdigit() or int(pmid) <= 0 or pmid in existing_pmids:
                        raise ValueError("invalid or duplicate PMID")
                    existing_pmids.add(pmid)
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"Cannot resume malformed corpus {out_path.name} line {line_number}: {exc}") from exc
    if not existing_pmids <= requested_pmids:
        raise ValueError("Cannot resume a corpus containing PMIDs outside the requested query snapshot")
    remaining_pmids = [p for p in pmids if p not in existing_pmids]
    if existing_pmids:
        print(f"[RESUME] {len(existing_pmids)} PMIDs already in {out_path.name}; fetching {len(remaining_pmids)} remaining...", flush=True)
    batch_size = 200
    written = len(existing_pmids)
    needs_record_separator = out_path.exists() and out_path.stat().st_size > 0 and not out_path.read_bytes().endswith(b"\n")
    with open(out_path, "a", encoding="utf-8") as fout:
        for i in range(0, len(remaining_pmids), batch_size):
            chunk = remaining_pmids[i : i + batch_size]
            params = {
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "xml",
            }
            api_key = os.environ.get("NCBI_API_KEY")
            if api_key:
                params["api_key"] = api_key
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            data_bytes = urllib.parse.urlencode(params).encode("utf-8")
            xml_bytes = b""
            for attempt in range(4):
                try:
                    req = urllib.request.Request(url, data=data_bytes, method="POST")
                    with urllib.request.urlopen(req, timeout=60) as resp:
                        xml_bytes = resp.read()
                    break
                except Exception as exc:
                    if attempt == 3:
                        raise
                    print(f"  [WARN] efetch retry {attempt + 1}/3 after {exc}...", flush=True)
                    time.sleep(1.5 * (attempt + 1))
            root = ET.fromstring(xml_bytes)
            batch_records = []
            for article in root.findall(".//PubmedArticle"):
                pmid_el = article.find(".//MedlineCitation/PMID")
                pmid = pmid_el.text if pmid_el is not None else ""
                title_el = article.find(".//Article/ArticleTitle")
                title = "".join(title_el.itertext()).strip() if title_el is not None else ""
                abs_texts = [
                    "".join(el.itertext()).strip()
                    for el in article.findall(".//Article/Abstract/AbstractText")
                ]
                abstract = "\n".join(t for t in abs_texts if t)
                yr_el = article.find(".//Article/Journal/JournalIssue/PubDate/Year")
                if yr_el is None:
                    yr_el = article.find(".//Article/Journal/JournalIssue/PubDate/MedlineDate")
                pub_year = yr_el.text[:4] if (yr_el is not None and yr_el.text) else None
                jour_el = article.find(".//Article/Journal/Title")
                journal = jour_el.text if jour_el is not None else ""
                doi = None
                pmcid = None
                for aid in article.findall(".//PubmedData/ArticleIdList/ArticleId"):
                    id_type = aid.attrib.get("IdType", "")
                    if id_type == "doi":
                        doi = aid.text
                    elif id_type == "pmc":
                        pmcid = aid.text
                mesh_terms = [
                    el.text
                    for el in article.findall(".//MedlineCitation/MeshHeadingList/MeshHeading/DescriptorName")
                    if el.text
                ]
                keywords = [
                    el.text
                    for el in article.findall(".//MedlineCitation/KeywordList/Keyword")
                    if el.text
                ]
                pub_types = [
                    el.text
                    for el in article.findall(".//Article/PublicationTypeList/PublicationType")
                    if el.text
                ]
                authors = article.findall(".//Article/AuthorList/Author")
                first_author = None
                senior_author = None
                if authors:
                    fa_last = authors[0].findtext("LastName") or authors[0].findtext("CollectiveName") or ""
                    fa_init = authors[0].findtext("Initials") or ""
                    first_author = f"{fa_last} {fa_init}".strip()
                    la_last = authors[-1].findtext("LastName") or authors[-1].findtext("CollectiveName") or ""
                    la_init = authors[-1].findtext("Initials") or ""
                    senior_author = f"{la_last} {la_init}".strip()
                rec = {
                    "pmid": pmid,
                    "pmcid": pmcid,
                    "doi": doi,
                    "pub_year": int(pub_year) if (pub_year and pub_year.isdigit()) else None,
                    "journal": journal,
                    "first_author": first_author,
                    "senior_author": senior_author,
                    "title": title,
                    "abstract": abstract,
                    "mesh_terms": mesh_terms,
                    "keywords": keywords,
                    "publication_types": pub_types,
                    "matched_pillars": pillar_map.get(pmid, []),
                    "has_open_pmc_fulltext": pmcid is not None,
                }
                batch_records.append(rec)
            returned_pmids = [rec["pmid"] for rec in batch_records]
            if len(returned_pmids) != len(set(returned_pmids)) or set(returned_pmids) != set(chunk):
                raise ValueError("PubMed efetch batch must contain exactly the requested unique PMIDs; checkpoint not advanced")
            if needs_record_separator:
                fout.write("\n")
                needs_record_separator = False
            for rec in batch_records:
                fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
                written += 1
            fout.flush()
            print(f"Fetched {written}/{len(pmids)} records...", flush=True)
            time.sleep(0.34)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description="Discover and harvest OpenIslet Evidence Corpus sources (v0.3).")
    parser.add_argument("--mode", choices=["core", "wide", "all_sources"], default="core", help="core (9 sub-pillars), wide (broad sensitivity), or all_sources (All 8 Channels)")
    parser.add_argument("--include-reviews", action="store_true", help="Include Review/Editorial publication types")
    parser.add_argument("--dry-run-counts", action="store_true", help="Print live counts, overlap breakdown, and specialty repo manifests without downloading full PubMed XML")
    parser.add_argument("--from-cached-ids", action="store_true", help="Skip esearch discovery and fetch XML metadata directly for PMIDs in discovered_pmids_{mode}.json")
    args = parser.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if args.from_cached_ids:
        id_manifest = RAW_DIR / f"discovered_pmids_{args.mode}.json"
        with open(id_manifest, "r", encoding="utf-8") as fid:
            cached = json.load(fid)
        pillar_map = cached["pmid_to_pillars"]
        all_pmids = sorted(pillar_map.keys(), reverse=True)
        out_file = RAW_DIR / f"corpus_manifest_{args.mode}.jsonl"
        print(f"[CACHED IDS] Loaded {len(all_pmids)} PMIDs from {id_manifest}. Fetching PubMed XML metadata...", flush=True)
        n_written = efetch_metadata(all_pmids, pillar_map, out_file)
        print(f"Done. Wrote {n_written} records to {out_file}", flush=True)
        return
    source_queries = WIDE_QUERIES if args.mode == "wide" else CORE_PILLARS
    pillar_map: dict[str, list[str]] = {}

    for name, q in source_queries.items():
        full_q = q if args.include_reviews else (f"({q}){NON_REVIEW_FILTER}")
        ids = esearch_pmids(full_q)
        print(f"[{args.mode.upper()}] {name}: {len(ids)} PMIDs")
        for pid in ids:
            pillar_map.setdefault(pid, []).append(name)

    pubmed_core_pmids = set(pillar_map.keys())
    print(f"\n[{args.mode.upper()}] Deduplicated PubMed Query PMIDs: {len(pubmed_core_pmids)}")

    # Golden Seed Recall Check
    caught_seeds = [s for s in GOLDEN_SEED_PMIDS if s in pubmed_core_pmids]
    missed_seeds = [s for s in GOLDEN_SEED_PMIDS if s not in pubmed_core_pmids]
    print(f"[RECALL CHECK] Golden Benchmark Seed Recall: {len(caught_seeds)}/{len(GOLDEN_SEED_PMIDS)} ({100*len(caught_seeds)/len(GOLDEN_SEED_PMIDS):.1f}%)")
    if missed_seeds:
        for ms in missed_seeds:
            print(f"  -> Missed Seed PMID:{ms} ({GOLDEN_SEED_PMIDS[ms]})")

    preprint_records: list[dict] = []
    patent_records: list[dict] = []
    geo_datasets: list[dict] = []

    if args.mode == "all_sources":
        # Channel 3 & 4: Europe PMC Full-Text Target Mentions + Preprints (bioRxiv/medRxiv)
        seen_ppr_ids: set[str] = set()
        for ep_name, ep_q in EUROPEPMC_SUPPLEMENTARY_QUERIES.items():
            ep_pmids, ep_pprs = europepmc_search_ids(ep_q)
            new_ep = ep_pmids - set(pillar_map.keys())
            print(f"[EUROPEPMC] {ep_name}: {len(ep_pmids)} PMIDs ({len(new_ep)} newly discovered) + {len(ep_pprs)} Preprints")
            for pid in ep_pmids:
                pillar_map.setdefault(pid, []).append(ep_name)
            for ppr in ep_pprs:
                ppr_id = str(ppr.get("id", ""))
                if ppr_id and ppr_id not in seen_ppr_ids:
                    seen_ppr_ids.add(ppr_id)
                    preprint_records.append(ppr)

        # Channel 5: 1-Hop Citation Snowballing (Cited-By + References of the 44 Golden Seeds)
        cited_by, refs = elink_citation_snowball(list(GOLDEN_SEED_PMIDS.keys()))
        snowball_union = cited_by | refs
        wide_union: set[str] = set()
        for w_name, w_q in WIDE_QUERIES.items():
            w_ids = esearch_pmids(f"({w_q}){NON_REVIEW_FILTER}")
            wide_union.update(w_ids)
        domain_filtered_snowball = snowball_union & wide_union
        new_snowball = domain_filtered_snowball - set(pillar_map.keys())
        print(f"[CITATION SNOWBALL] 44 Golden Seeds -> {len(cited_by)} citing + {len(refs)} referenced = {len(snowball_union)} raw PMIDs")
        print(f"[CITATION SNOWBALL] Filtered to primary islet/beta-cell/T1D domain: {len(domain_filtered_snowball)} PMIDs ({len(new_snowball)} newly discovered beyond Core 9 Sub-Pillars)")
        print(f"[WIDE-NET SENSITIVITY POOL] Total Wide-Net Union: {len(wide_union)} primary non-review PMIDs")
        for pid in domain_filtered_snowball:
            pillar_map.setdefault(pid, []).append("Pillar5_Citation_Snowball_of_44_Golden_Seeds")

        # Channel 7: Europe PMC Patent Literature (SRC:PAT)
        for pat_name, pat_q in EUROPEPMC_PATENT_QUERIES.items():
            pats = europepmc_search_patents(pat_q)
            print(f"[PATENTS] {pat_name}: {len(pats)} deduplicated patent records")
            patent_records.extend(pats)
        if patent_records:
            pat_file = RAW_DIR / "patents_sc_islets_hypoimmune.jsonl"
            with open(pat_file, "w", encoding="utf-8") as fpat:
                for pat in patent_records:
                    fpat.write(json.dumps(pat, ensure_ascii=False) + "\n")
            print(f"[PATENTS] Saved {len(patent_records)} patent records to {pat_file}")

        # Channel 8: NCBI GEO (db=gds) Functional Genomics & Single-Cell Datasets + Linked PMIDs
        geo_datasets, geo_linked_pmids = geo_gds_search_datasets(GEO_GDS_QUERIES)
        new_geo_pmids = geo_linked_pmids - set(pillar_map.keys())
        print(f"[GEO GDS] Total islet-verified GSE series: {len(geo_datasets)} | Linked PMIDs: {len(geo_linked_pmids)} ({len(new_geo_pmids)} newly discovered PMIDs)")
        for pid in geo_linked_pmids:
            pillar_map.setdefault(pid, []).append("Channel8_GEO_Functional_Genomics_Linked_PMID")
        if geo_datasets:
            geo_file = RAW_DIR / "functional_genomics_geo_datasets.json"
            with open(geo_file, "w", encoding="utf-8") as fgeo:
                json.dump({
                    "generated_at": "2026-10-20",
                    "total_gse_series": len(geo_datasets),
                    "total_linked_pmids": len(geo_linked_pmids),
                    "newly_discovered_pmids_from_geo": len(new_geo_pmids),
                    "datasets": geo_datasets,
                }, fgeo, indent=2, ensure_ascii=False)
            print(f"[GEO GDS] Saved {len(geo_datasets)} GSE dataset records to {geo_file}")

        # Save preprint manifest
        if preprint_records:
            ppr_file = RAW_DIR / "preprints_sc_islets_hypoimmune.jsonl"
            with open(ppr_file, "w", encoding="utf-8") as fp:
                for ppr in preprint_records:
                    fp.write(json.dumps(ppr, ensure_ascii=False) + "\n")
            print(f"[PREPRINTS] Saved {len(preprint_records)} deduplicated bioRxiv/medRxiv/PPR records to {ppr_file}")

    all_pmids = sorted(pillar_map.keys(), reverse=True)
    overlap_stats = compute_overlap_breakdown(pillar_map)

    print(f"\n[SUMMARY] Total Deduplicated Primary PubMed PMIDs: {len(all_pmids)}")
    print(f"[OVERLAP] Multi-Sub-Pillar Hit Distribution: {json.dumps(overlap_stats['subpillar_hit_distribution'])}")
    print(f"[OVERLAP] Multi-Major-Pillar Hit Distribution: {json.dumps(overlap_stats['major_pillar_hit_distribution'])}")
    print(f"[OVERLAP] SC-Islet Multi-Gate Intersections: {json.dumps(overlap_stats['sc_islet_multi_gate_intersections'])}")
    print(f"[OVERLAP] Pairwise Major Pillar Intersections: {json.dumps(overlap_stats['major_pillar_pairwise_intersections'])}")

    if args.dry_run_counts:
        id_manifest = RAW_DIR / f"discovered_pmids_{args.mode}.json"
        with open(id_manifest, "w", encoding="utf-8") as fid:
            json.dump({
                "generated_at": "2026-10-20",
                "mode": args.mode,
                "total_deduplicated_pmids": len(all_pmids),
                "total_preprints": len(preprint_records),
                "total_patents": len(patent_records),
                "total_geo_gse_series": len(geo_datasets),
                "golden_seed_recall": f"{len(caught_seeds)}/{len(GOLDEN_SEED_PMIDS)}",
                "overlap_breakdown": overlap_stats,
                "pmid_to_pillars": pillar_map,
            }, fid, indent=2)
        print(f"[MANIFEST] Saved deduplicated PMID & source map ({len(all_pmids)} PMIDs) to {id_manifest}")
        return

    out_file = RAW_DIR / f"corpus_manifest_{args.mode}.jsonl"
    n_written = efetch_metadata(all_pmids, pillar_map, out_file)
    print(f"Done. Wrote {n_written} records to {out_file}")


if __name__ == "__main__":
    main()



