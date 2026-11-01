# OPENISLET.ORG (`PROJECT ISLET-1`) — Ring 1 Cross-Study Discovery Atlas, Whole-Islet Architecture & Dual-Structure Strategy (`v0.7.0`)

<!-- OPENISLET_GENERATED_FACTS_START -->

**Current executable facts (generated from on-disk artifacts; physical gates pending).**

Current corpus: 563 studies, 548 active deduplicated, 547 active empirical, 529 strict peer-reviewed empirical; 3,093 effects, 592 arms, 20 reference controls, 231 edit entries and 236 expanded gene occurrences under the ranker expander. Direction/penalty labels are rule-based heuristics; 43 flagged penalties are not an independently adjudicated accuracy dataset.

Layer 1 evaluates 2,160 prior-anchored stacks (480 feasible, 48 code-defined top tier). Registry: 48 source terms / 46 net offsets / 39 auxiliary scalar leaves / 6 engineering priors. Structural share of component variance is 99.8542% from unrounded latent logits. Current utility-target profiles with typed numeric LRR evidence: 0; the shrinkage branch is checked separately against an actual HNF1A-KO source contrast when this count is zero. Seeded direct and shared-parameter Monte Carlo are exploratory sensitivity simulations with finite-draw Wilson intervals and rank quantiles, not clinical probabilities or fitted posterior intervals. CHK-12 recomputes and compares the complete current replay.

HIP-Core executes all 36 coded rows, including all 11 repeat/bridge/lead rows, across 25 distinct genotypes. The quoted covariance/power applies to 108 independent run-nested clones or 75 shared-genotype clones; 61-clone power is UNEVALUATED with null SE/MDE. Current factor sums are [-10, -10, -6, -10, -6, -10]; +1 counts are [13, 13, 15, 13, 15, 13].

| Construct | Donor bp | Packaged bp, including 2 x 145-bp ITRs |
| :--- | ---: | ---: |
| `LOCUS1_ARM1A_TPT2A_BICISTRONIC` | 4880 | 5170 |
| `LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01` | 3974 | 4264 |
| `LOCUS3_CXCL10_SERPINB9_TPT2A_PDL1_R01` | 4090 | 4380 |
| `LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_ICASP9_R17` | 3457 | 3747 |

Native EEF1A1 core: 230 bp, 19 CpG dinucleotides, 3 exact/RC donor occurrences; these are native CpG-containing promoter sequences.

Release manifest: `1,351` pinned files. Machine facts: `data/extracted/documentation_facts_v0.8.0.json`. These facts supersede stale numeric summaries in explicitly historical rationale sections. Outward rankings remain unfrozen (`outward_rankings_frozen: false`); no exact product has passed physical qualification.

<!-- OPENISLET_GENERATED_FACTS_END -->

> **Earlier-revision rationale:** The body originated in an earlier design revision and retains historical reasoning. Unrendered historical values are not current release claims; the generated facts above and current machine-readable artifacts govern active numerical outputs.

**Date:** 2026-10-23  
**Project Root:** `./` (`OpenIslet` repository root)  
**Primary Graph Asset:** `data/extracted/cross_study_discovery_graph.json`  
**Extraction Schema:** `schemas/extraction_schema_v0.2.json`  

---

## 1. Executive Summary & Current Project State

**Project Islet-1 (`OpenIslet.org`)** is an open-science discovery commons paired with a translational cell-therapy engineering program to design an immunosuppressant-free, hypoimmunogenic ("hypoallergenic") 3D stem-cell-derived islet (`SC-islet`) replacement therapy for Type 1 Diabetes (`T1D` — addressing both early-onset pediatric T1D and adult-onset `LADA`).

### Completed Pipeline Milestones (Steps 1–5)

| Milestone | Validated Metric | File / Directory Location |
| :--- | :---: | :--- |
| **Step 1: Multi-Source Literature & Omics Harvest** | **`7,703` deduplicated manifest records** (`6,859` PubMed + `196` GEO `GSE` series + `146` patents + `70` preprints + `114` clinical trials; `45/45 = 100%` Golden Seed recall) | `scripts/harvest_corpus.py`, `data/raw/corpus_manifest_all_sources.jsonl` |
| **Step 2: Multi-Gate Identification & Stratification** | **`563` Ring 1 Platinum & Counter-Evidence primary studies** (plus `2,175` Ring 2 High-Confidence Multi-Source studies) | `scripts/identify_and_stratify_corpus.py`, `data/raw/corpus_stratification_manifest.json` |
| **Step 3: Relational Extraction Schema (`v0.2.0`)** | 5-Gate arm covariates, atomic effects ($X \to Y \mid \mathbf{Z}$), effect modifiers ($\partial \text{Effect}/\partial Z \neq 0$), causal edges ($A \to B \to C$), and cross-study discovery bridges | `schemas/extraction_schema_v0.2.json` |
| **Step 4: Ring 1 Study Data Dossiers** | **`563` / `563` Ring 1 dossiers built** (`282` Full-Text Europe PMC JATS XML + `281` Abstract+GEO/MeSH, plus `2` auxiliary cited dossiers = `565` total dossier files) | `scripts/build_study_dossiers.py`, `data/dossiers/{pmid}.json` |
| **Step 5: Ring 1 Relational Extraction & Cross-Study Discovery Graph** | **`563` / `563` schema-validated study JSONs** (`0` schema errors; `46` deeply curated landmarks + `517` dossier-grounded extractions across `2,160` combinatorial scenario stacks) | `scripts/extract_ring1_corpus.py`, `scripts/build_cross_study_graph.py`, `data/extracted/{pmid}.json`, `data/extracted/cross_study_discovery_graph.json` |
| **Human Genetics & Safety Priors** | Open Targets (`MONDO_0005147` T1D) + `gnomAD v4` loss-of-function constraint (`pLI`, `LOEUF`) across `33` raw-backed targets | `scripts/summarize_omics_baseline.py`, `data/omics/target_prior_table.json` |

---

## 2. Ring 1 Cross-Study Discovery Graph Statistics (`v0.7.0`)

Across the **`563` validated Ring 1 studies (`548` active deduplicated empirical studies)** in `data/extracted/cross_study_discovery_graph.json`:
- **Experimental & Control Arms ($X \text{ vs. } X_0$):** `592` (`20` explicit reference-control comparator arms; `0` templated `_arm0_wt_control`)
- **Quantitative & Qualitative Atomic Effects ($X \to Y \mid \mathbf{Z}$):** `3,093` (`0` active utility-target comparative numeric LRR profiles; `4` auxiliary `HNF1A` comparative numeric LRR effects and `3` manipulation checks across `7` total typed numeric LRR records; `pmid36506044_eff2` transgene-induction manipulation check quarantined)
- **Effect-Modifier Interactions ($\partial \text{Effect}/\partial Z \neq 0$):** `1,282`
- **Directed Causal Mechanism Edges ($A \to B \to C$):** `1,093`
- **Residual Failure & Escape Vectors ($<100\%$ survival bottlenecks):** `751`
- **Shared Cross-Study Bridge Nodes:** `42` multi-study hubs (capped deterministically at `500` heuristic multi-hop hypothesis paths via `multi_hop_hypothesis_paths_enumeration_policy`)
- **Perturbed/Evaluated Target Gene Hubs (`target_gene_atlas`):** `39` distinct genetic/epitope target hubs.

---

## 3. Why We Engineer the Entire Hypoimmunogenic 3D Islet (`SC-Islet`), Not Just Isolated $\beta$ Cells

A core architectural principle established by our Cross-Study Discovery Graph is that the therapeutic graft must be a **complete, multi-hormone, hypoimmunogenic 3D stem-cell-derived islet (`SC-islet`) mini-organ** (`~120–150 µm` diameter) plus co-assembled microvessels, rather than isolated $\beta$ cells:

### 3.1 Physiological Requirement: Paracrine Counter-Regulation & 3D Bioenergetics
1. **$\beta$ Cells (`~55–65%`, `INS+ / NKX6.1+ / MAFA+`):** Secrete insulin and amylin (`IAPP`) in a biphasic manner to clear glucose.
2. **$\alpha$ Cells (`~15–25%`, `GCG+`):** Secrete glucagon to prevent hypoglycemia and provide local intra-islet `GCG / GLP-1 -> cAMP` tone required to prime $\beta$-cell first-phase insulin release.
3. **$\delta$ Cells (`~5–10%`, `SST+`):** Secrete somatostatin (`SST`), acting as the local paracrine safety brake on $\beta$ and $\alpha$ cells to prevent post-prandial insulin overshoot and life-threatening hypoglycemia (critical in young pediatric recipients).
4. **3D Aggregation Requirement for Mitochondrial OXPHOS:** `MacDonald & Yi 2025` (`PMID:39988068`, `Fig. 6C-F`) and `Balboa & Otonkoski 2022` (`PMID:35241836`) demonstrated that SC-$\beta$ cells in 2D planar culture fail to increase mitochondrial respiration (`OCR`) or synchronize `[Ca2+]i` oscillations; only in **3D aggregated `SC-islets`** (with off-target serotonin-producing `SLC18A1+` enterochromaffin cells depleted via Stage 6 Aurora kinase inhibitor `ZM447439`) do they achieve adult-like mitochondrial coupling and biphasic GSIS.

### 3.2 Immunological Requirement: Why Uncloaked $\alpha$, $\delta$, or Endothelial Cells Destroy the Graft
1. **In Vivo Classical HLA Upregulation Is Highest in SC-$\alpha$ Cells:**
   - Connecting `Augsornworawat & Millman 2020` (`PMID:32846125`, `GSE151117`), `Gerace & Melton 2023` (`PMID:36599351`, `Fig. 1C-E`), and `Demine & Eizirik 2020` (`PMID:31900242`, `Fig. 6a-c`) reveals that after in vivo transplantation and under `IFN-γ + IL-1β` challenge, **`GCG+` SC-$\alpha$ cells upregulate classical `HLA-A/B/C` even more strongly than `INS+` SC-$\beta$ cells**.
   - Consequently, if immune-cloaking transgenes (`CD47`, `PD-L1`, or tolerogenic cytokines) were restricted to an insulin (`INS`) promoter, host alloreactive T cells and antibodies would attack the uncloaked $\alpha$ and $\delta$ cells and destroy the entire islet via bystander `IFN-γ / TNF-α / Granzyme B` release.
2. **Stem-Cell Endothelial Cells (`CD31+`) Express High NKG2D Ligands (`MICA/MICB`, `ULBP1/2`):**
   - While `B2M-/-` SC-$\beta$ cells intrinsically lack `MICA/MICB` and `ULBP1/2` (`Gerace 2023`, `Fig. 3A-E`), **`B2M-/-` SC-endothelial cells express high `MICA/MICB` and `ULBP1/2` (`Fig. S4C-D`)**, making co-assembled blood vessels (`Aghazadeh & Nostro 2021`, `PMID:34480863`) prime targets for primary human NK lysis unless protected by **pan-cellular high-threshold `CD47` overexpression (`Hu & Schrepfer 2024`, `PMID:37156915`)**.
3. **Solution — Stage 0 Master Pluripotent Stem Cell (`hPSC`) Engineering:**
   - By introducing all edits into the **Stage 0 master `hPSC` clone** prior to differentiation—using **pan-cellular housekeeping loci (`GAPDH` or `CAG/EF1α`)** for universal stealth transgenes (`CD47`, `2B10` Treg cytokines) alongside endogenous knockouts (`B2M-/-`, `CIITA-/-`, `CXCL10-/-`, `RNLS-/-`) and endogenous $\beta$-cell-specific edits (`INS-SP^{CHGA}` signal peptide swap, `PMID:42690930`)—**every endocrine cell ($\beta, \alpha, \delta$) and hPSC-derived microvessel in the differentiated 3D `SC-islet` automatically inherits the complete hypoimmunogenic cloak.**

---

## 4. Seven Emergent Multi-Gate Engineering Laws from the Ring 1 Graph

1. **Selective `CXCL10-/-` Beats Clonal `STAT1-/-` in Full Allo-PBMC Settings (`PMID:36055241` $\times$ `PMID:36289205` $\times$ `PMID:31900242`):**
   - Clonal `STAT1-/-` (`Sintov & Melton 2022`, `Fig. 4D-J`) fails in human PBMC co-culture and increases NK apoptosis because deleting `STAT1` strips away `IFN-γ`-inducible protective checkpoints **`PD-L1` (`CD274`), `HLA-E`, and `SOCS1`**. Conversely, clonal **`CXCL10-/-`** breaks the `CXCR3+` T/NK recruitment loop while preserving inducible `PD-L1`, `HLA-E`, and `SOCS1`, achieving 100% graft survival for 18 weeks in humanized mice (`Fig. 5B-D`). If `TYK2` is inhibited (`Chandra & Otonkoski 2022`), checkpoints (`CD47`/`PD-L1`) must be constitutively expressed rather than `IFN`-inducible.
2. **Primary Human `CD56+` NK Cells vs. `NK-92` Artifacts & `HLA-E` Leader Peptide Dependence (`PMID:36599351` $\times$ `PMID:34407395` $\times$ `PMID:37156915`):**
   - Knocking out `HLA-A/B/C` (`HLA-ABCnull`) abolishes surface `HLA-E` even when `B2M` and `HLA-E` are intact because `HLA-E` requires the classical HLA leader nonamer (`VMAPRTLVL`) for ER export (`Parent & Hebrok 2021`, `Fig. 2H-I`). Furthermore, `HLA-E` only inhibits the `NKG2A+` subset (`~40–50%`) of primary human NK cells (`Fig. S3B`), whereas `NK-92 MI` (`96.1% NKG2A+`) overestimates `HLA-E` protection. High-threshold **`CD47` overexpression (`Hu & Schrepfer 2024`)** inhibits both primary NK cells and macrophages across the entire graft.
3. **Non-HLA Co-Stimulatory Escape (`HVEM`, `CD40`, `BTN3A1`) Requires Local Treg or High-Threshold `CD47` Pairing (`PMID:36599351`):**
   - `IFN-γ` upregulates non-HLA T-cell co-activating ligands `HVEM (TNFRSF14)`, `CD40`, and `BTN3A1` on `B2M-/-` SC-$\beta$ cells (`GSE200021`, `Fig. 2E`), causing delayed rejection at 7 weeks unless paired with localized Treg recruitment (`GAPDH::IL-2[N88D]-TGFB1-IL10` `2B10` cassette, `Fig. 4E-H`) or high-threshold `CD47` (`Hu 2024`).
4. **Stage 7 / In Vivo Maturation Simultaneously Upregulates `MAFA` and T1D Autoantigens (`G6PC2`, `CHGB`, `INS` Leader, `IAPP`) (`PMID:35241836` $\times$ `PMID:32846125` $\times$ `PMID:42690930`):**
   - Maturation into `MAFA+` adult-like $\beta$ cells inevitably elevates T1D autoantigens (`G6PC2/IGRP`, `CHGB`, `PPI15-24` leader peptide) and amyloidogenic `IAPP`. Swapping the insulin signal peptide (`INS-SP^{CHGA}`, `Kobaisi & Scharfmann 2026`) completely abolishes `HLA-A2/PPI15-24` CD8+ T-cell recognition without impairing proinsulin folding or GSIS.
5. **Genetic `RNLS-/-` Glycolytic/NADPH Redox Shielding vs. Human iPSC-$\beta$ Dynamic GSIS Trade-Off (`PMID:39988068` $\times$ `PMID:32719542` $\times$ `PMID:32814902` $\times$ `PMID:35241836`):**
   - **RNLS precaution:** MacDonald 2025 (PMID:39988068 Fig. 6) reports poor dynamic GSIS in both WT and RNLS-KO iPSC-beta aggregates. Higher basal OCR in KO motivates a precautionary metabolic penalty, not a demonstrated knockout-specific GSIS decrement; the exact-product dynamic-perifusion comparison is pending.
6. **Multi-Pulse `IFN-γ` Epigenetic Memory (`MPS`) vs. Continuous `IFN-γ + IL-1β` Dedifferentiation (`PMID:32814902` $\times$ `PMID:31900242`):**
   - Continuous 48h `IFN-γ + IL-1β` triggers `sXBP1/CHOP` ER stress and selective `MAFA/NKX2.2/INS` dedifferentiation (`Demine 2020`), whereas brief pulsatile `IFN-γ` conditioning (`MPS`: `3 × 2h at 10 ng/mL` with recovery washes, `Yoshihara 2020`) imprints stable `ATAC-seq` accessibility at `CD274 (PD-L1)` without sustained `STAT1/IRF1` opening or dedifferentiation.
7. **Pre-Assembled `CD31+CD146+` Microvessels + `SA-FasL` Microgels Solve the 72-Hour Subcutaneous Ischemic Bottleneck (`PMID:34480863` $\times$ `PMID:35559682`):**
   - Ready-made `CD31+CD146+` microvessels (`Aghazadeh & Nostro 2021`) anastomose with host vessels in `<3–5 days` in unmodified subcutaneous tissue (whereas single-cell HUVECs fail), and co-delivery with `SA-FasL` PEG-4MAL microgels (`Lei & García 2022`) induces local `FoxP3+` Treg tolerance in primates without chronic calcineurin-inhibitor toxicity.

<!-- OPENISLET_CURRENT_LAYER_SCOPE -->
Layer 2 uses the pinned release pairing Ensembl release 112 / GENCODE v46 (GRCh38.p14, 1-based closed coordinates; ENST00000309268.11 reconciled) and contains 11 literature-anchored stage-specific prior profiles with 1,000 independent prior-scenario draws per endpoint; it is not a fitted dynamic GRN/ODE or held-out validated differentiation predictor. Arm 1B constitutive EEF1A1 drives antisense CD47 and inducible sense 5xISRE drives HLA-E. Computed Arm 1A/1B Hour-0 CD47 is 209,444/238,004 molecules/cell; Arm 1B/WT is 4.96679814x. Resting Arm 1B HLA-E is 10,890 and IFNg24h HLA-E 69,516; the resting endpoint is below the assumed 15,000 threshold. Exact edited beta/alpha/delta/iEC density is unmeasured, with a preregistered Stage-0 acceptance target for Hour-0 density. TYK2-KO Stage-6 prior is 0.364 versus 0.678 WT (-46.3%); this is distinct from the transcript-level endpoint. HNRNPA2B1-proximal 455-bp CpG-core candidate (GRCh38 chr7:26200580-26201034); anti-silencing function requires head-to-head lineage qualification against the published CBX3-side fragment. The local 3,118-row AVI table is quarantined in the private cache: end-to-end batch generation provenance for the full 3,118-row table is incomplete (as a whole it is not a certified model scan), a limited 28 Oct 2026 live re-run against alphagenome v0.10.0 verified 59/59 sampled rows (9 candidate edits plus 50 random rows), AlphaGenome non-commercial output terms conflict with CC-BY-4.0 public redistribution, and zero AlphaGenome outputs are shipped publicly or used as active ranking weights or insertion-safety certifications. Scenario-B P2A profiles are unsupported provisional comparators; the exact 2A-free Scenario-B donor/prediction remains unavailable.

Reference-geometry compatibility screen: 3CDG is a CD94/NKG2A complex. P8-Phe8 contributes 24 heavy-atom contacts to shared CD94 and 0 to NKG2A in copy 1, and 21/0 in copy 2. Removing Phe ring atoms is an atom-truncation diagnostic, not a relaxed Val8 model. This structure cannot discriminate NKG2A from NKG2C. Receptor-specific K_D values are working soluble-ectodomain reference priors; direct SPR and HCMV+ NKG2C+ NK-donor assays remain preregistered. Reference SIRPA V1/V2 sequence and geometry comparison; pan-allelic protection is unverified. Asp101 is inserted in SIRPA V1 (4CMM) relative to V2 (2JJS). Soluble crystal CD47 carries C33G/Cys15-to-Gly; engineered CD47 retains WT Cys15. Reference soluble-C15G structures place CD47 Val115 about 30.93 A from the nearest SIRPA atom. This coarse anchor-distance screen does not model membrane orientation, the engineered WT-Cys15/GPI stalk or V8; binding, membrane and phagocytosis assays remain pending. The 21-residue linker-2 path includes unresolved HLA-E Gly1 between the 20-aa (G4S)4 linker and HLA-E Ser2 N; 73.50 A contour, 0.4234 extension, 9.901 pN are coarse WLC quantities. The 3CDG reference-crystal Arg107 C-alpha distance is 17.02 A, nearest-heavy-atom distance 14.82 A; neither is a measured engineered Gly107 geometry. Foldseek is a cached-hit replay with missing historical database/tool/date provenance, not a reproducible original search. 8H4X was superseded by 8ZCR (Yan, Zhou et al., 2024, Acta Crystallogr F 80:286-293; DOI:10.1107/S2053230X24009439; PMID:39382088). Exact-product function remains unresolved.

Layer 4 is a conditional model, not validated graft survival. The physical source lamina half-width fraction is 0.34; its continuous central volume is derived from that same geometry and control-volume overlap at every mesh. Nominal minimum slab oxygen is 19.92 mmHg at 18 h. The 8.64-mmHg homogenized fixed-surface Dirichlet core is an upper-bound screen; operative uniform-supply / conservative gel-only two-zone cores are 4.08/2.15 mmHg. Physical supply topology remains a Stage-2 assay gate; mesh error is evaluated on unrounded values, with NUMERICALLY_UNRESOLVED whenever an uncertainty interval straddles the threshold. Reference fixed-functional-beta-IEQ ODE curves have pde_coupling_applied=false and cannot inherit a PDE-derived full-product verdict; matched PDE/ODE states are identified separately. Reference ODE outputs: OPENISLET_60_30_10_TRIHORMONAL_G6PC2_WT: minimum 4.59 mM / insulin shutoff t50 4.93 min; BETA_ONLY_MONOCULTURE_NO_ALPHA_DELTA: minimum 3.06 mM / insulin shutoff t50 9.24 min; G6PC2_KO_DISRUPTED_GLYCOLYTIC_RHEOSTAT_GY002: minimum 2.70 mM / insulin shutoff t50 10.44 min. Low-glucose challenge reference-ODE minima: lead 3.9218, composite-reflex OAT1 3.8805, pure-beta-clamp OAT1A 3.9186 mM; these share the specified prior kinetic constants and are not measured clinical protection. Nominal modeled peak peroxide is 0.725 uM and maximum pH 7.405.

Layer 5 is an illustrative assigned-rule population simulator, not a clinical forecast. Personalized analytic expectation is 97.5000%, standardized Plan A analytic 92.5000%; sampled personalized 97.40%. Analytic ancestry gap is 0.40 percentage points; sampled gap 1.80 points, with finite sample/seed scope. Genome-wide off-target enumeration is NOT_ENUMERATED and coding-hit count is null; assigned illustrative three-mismatch panel scores are not exhaustive genome/exome scans or empirical HiFi-SpCas9 specificity. Variant identity is source-verified dbSNP interval identity; rare-variant screening frequencies are explicit assumed priors, with gnomAD row-level AC/AN validation pending. Manufacturing yield is null/NOT_MEASURED; exp(-lambda_K*max(0,K-4))*gamma_B**B is a dimensionless engineering utility discount. The executable optional TYK2 design is PEmax+epegRNA. CGBE remains an unimplemented contingency requiring an engineered SpG/SpRY-compatible PAM, editing-window/bystander ledger and separate validation; no CGBE guide is qualified here.
<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->
