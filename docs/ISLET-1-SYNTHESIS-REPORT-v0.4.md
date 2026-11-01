# OPENISLET.ORG (`PROJECT ISLET-1`) — Synthesis Report, Scenario-Output Freeze & Retirement of Undefined `Chassis C/D` (`v0.4` / `v0.7.0 Wave 0`)

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

**Author:** Chief Scientist (`Gemini 4 Argon (High Effort)`, incorporating audits from `Claude 5.5 Sonnet Max Effort`, `ChatGPT 6 Astra xhigh`, `Kimi K3.1 Pro xhigh`, `ChatGPT 6.1 Pro / Sol 6.1`, `Spark Muse Wonder`, `Grok 4.7 Heavy`, and the 5-Agent Red-Team Swarm)  
**Date:** 2026-10-23  
**Project Root:** `./` (`OpenIslet` repository root)  
**Release Manifest (`P0-8`):** [`manifest.json`](../manifest.json) (`1,351` pinned files + `public_redistribution_allowlist_spec` + `export_public_bundle` with `<=240`-char non-OA per-article excerpt cap)  
**Companion Blueprint & Evaluator Spec:** [`docs/ISLET-1-CANDIDATE-BLUEPRINTS-AND-OPTIMIZER-v0.2.md`](ISLET-1-CANDIDATE-BLUEPRINTS-AND-OPTIMIZER-v0.2.md)

---

## 1. Purpose of This Synthesis Document (`P0-8` & `AT-00` Resolution)

During the Round-3 and Revision-2 adversarial audits (`REVIEW-2026-10-21.md` and `OPENISLET-v0.6.8-AUDIT-REPORT.md`), the external review brief referenced four legacy/brief path names (`data/extracted/ring1_study_dossiers_v0.5.jsonl`, `data/extracted/ring1_atomic_effects_v0.5.jsonl`, `data/extracted/islet1_stack_ranking_v0.5.json`, and `docs/ISLET-1-SYNTHESIS-REPORT-v0.4.md`) as well as two undefined candidate names (`Chassis C` and `Chassis D`).

To eliminate any drift between external review briefs and the repository working tree:
1. **OpenIslet & Islet-1 Framing:** **OpenIslet** is a global community project using AI (`Layer 1: Frontier Reasoning Auditors` + `Layer 2: HuggingFaceBio/Carbon-X 14.8B Local Bilingual DNA Co-Pilot` + `Layer 3: AlphaGenome / Evo 2 / Bridge-RNA / AlphaFold / Mol*`) to harness scientific knowledge and collective wisdom at scale to accelerate a cure for Type 1 Diabetes. Together, we are developing **Islet-1**—an open-source blueprint for a personalised stem-cell islet designed to cure T1D without anti-rejection drugs, spanning **Stage 0** (The `Islet-1` Blueprint & Construct Pre-Pilot), **Stage 1** (Automated Wet-Lab Design & Robotic Screening), **Stage 2** (3D Vascularized Organoid & Oxygen-Bridge Assembly), **Stage 3** (Preclinical In-Vivo & Safety-Switch Target Verification), **Stage 4** (Robotic GMP Manufacturing & Phase 1/2 Trial), and **Stage 5** (Genetically & Immunologically Matched Outpatient Cure stratified across `HLA-C/KIR`, `HLA-DR3/DR4 + Autoantibodies`, `mtDNA` haplogroup, and `ABO Type-O` vascular backbone where biologically justified).
2. All canonical repository files, SHA-256 checksums, and record counts are now pinned in the machine-verifiable root [`manifest.json`](../manifest.json) (`scripts/generate_release_manifest.py`).
  - **Heuristic extraction labels:** Current corpus: 563 studies, 548 active deduplicated, 547 active empirical, 529 strict peer-reviewed empirical; 3,093 effects, 592 arms, 20 reference controls, 231 edit entries and 236 expanded gene occurrences under the ranker expander. Direction/penalty labels are rule-based heuristics; 43 flagged penalties are not an independently adjudicated accuracy dataset.

---

## 2. Formal Retirement of `Chassis C` and `Chassis D`

Neither `Chassis C` nor `Chassis D` was ever defined in the `Project Islet-1` codebase; both originated in an external review prompt template. Based on empirical verification against the `563`-study Ring 1 corpus and primary literature, **both are formally retired**:

1. **Retired `Chassis C` ("`TM4SF4+` vascular-first minimally edited niche"):**
   - **Why It Is Retired (`Marker Confusion`):** In the human pancreatic islet literature (`Baron et al. 2016 PMID:27667365`), `TM4SF4` is a human pancreatic $\alpha$-cell/ductal tetraspanin rather than a vascular endothelial marker (and in our extracted Ring 1 corpus [`PMID:35918471`, `PMID:34035238`], `TM4SF4` appears as a transcript downregulated in `HNF1A`-deficient $\beta$-like/progenitor cells).
   - **What Is Preserved:** The legitimate scientific concept of an **isogenic minimally edited (`HIP` core: `B2M-/- + CIITA-/- + CD47`) + pre-vascularized pouch control arm** is preserved as the baseline control in the **HIP-Core Fractional-Factorial Pilot** (`§3` below).
2. **Retired `Chassis D` ("`CADM1-/- + CXCL10-/- + TNFAIP3 + HLA-A*02:01-/-` T1D-specific evasion"):**
   - **Why It Is Retired (`Category Mismatch & Unvalidated Knockout`):**
     - `HLA-A*02:01-/-` allele-specific deletion is an autologous or single-allele tactic that leaves all remaining classical HLA-A/B/C and HLA-II allotypes intact (and `IFN-α` shifts $\beta$-cell peptidome presentation disproportionately toward `HLA-B`, `PMID:39824805`).
     - `CADM1` mediates islet nerve/$\alpha$-cell gap-junction coupling and has zero hPSC/SC-islet knockout validation.
   - **What Is Preserved:** Orthogonal, corpus-supported immune-evasion candidates—specifically the **`CD155` (`PVR` `rs1058402` Ala67Thr) high-affinity `TIGIT`-agonist variant** (`Brown et al. 2025`, `PMID:41202120`), **`CD155/B7-H3` activating-ligand ablation** (`Chimienti et al. 2022`, `PMID:36170817`), **`CD54/CD58` adhesion ablation** (`Hammer et al. 2024`, `PMID:38981470`), and **`TNFAIP3 (A20)` inducible cytokine resilience**—are promoted into the **HIP-Core Fractional-Factorial Pilot** (`§3`).

---

## 3. Five Joint-PM Governance Rulings (`Wave 0` Active State)

1. **Outward-Facing Ranking Freeze & Scenario-Output Relabelling:**
   - Outward-facing rankings and `"Named #1/#2"` labels are frozen and relabelled as **Prior-Anchored Scenario Outputs (`Scenario Output A`, `Scenario Output B`)** governed by stated constants (`engineering_priors_and_sensitivity_table` and `complete_gate_coefficient_registry`). This freeze lifts only with Kev's explicit approval after Waves 1–2 close.
2. **Retirement of `Chassis C` and `Chassis D`:**
   - Recorded above and enforced across all documentation.
3. **Independent In-Repo Acceptance & Verification Suites (`tests/`):**
   - Verified via `tests/test_adversarial_acceptance.py` (`16/16 PASS`, self-verifying both in the full repo and inside the exported public bundle), `tests/test_citations_resolve.py`, and `tests/test_design_algebra_exact.py`.
4. **Wet-Lab Priority — The HIP-Core Fractional-Factorial Pilot First:**
   - Rather than testing single modules in isolation, wet-lab execution (subject to Stage-0 gate review) prioritizes an isogenic, matched-viable-$\beta$-dose fractional-factorial pilot anchored on the clinically supported **`B2M-/- + CIITA-/- + CD47` (`HIP`) core** (`PMID:40757665`, `PMID:42430558`, `PMID:40812299`), testing whether each incremental edit (`HLA-E` single-chain trimer peptide panel $\times$ $\alpha 3$-domain TCR-attenuated mutant [`PMID:39621955`], `CXCL10::SERPINB9-GSG-P2A-PD-L1` Locus 3, `CD47-IgV-GPI` vs `CD47-dCT` vs full-length `CD47-WT` [`PMID:38508139`] and 2A-free bidirectional cassette vs bicistronic `tPT2A`, and `CD155` engineering [`PMID:41202120`]) delivers statistically significant marginal protection across primary human NK, macrophage (`SIRPA` V1/V2 stratified [`PMID:40668892`]), bystander cytokine, and indirect `HLA-DQ2/DQ8` APC channels while preserving $<3\text{ min}$ dynamic GSIS $2.8\text{ mM}$ shut-off.
5. **Read-Only MCP Gate:**
   - The OpenIslet MCP server exposes strictly **read-only** evidence tables and Graveyard checks (`search_studies`, `get_study`, `get_comparison`, `get_conflicting_evidence`, `check_stack_against_graveyard`, `get_release_manifest`) with **no rankings and no write tools** until `P0-3`, `P0-5`, and `P0-6` pass.

Current encoded element: single 22-aa GSG-P2A element (21-aa upstream stub plus downstream Pro; legacy tPT2A alias). SCT encodes B2M signal MSRSVALAVLALLSLSGLEA (20 aa); ER cleavage and mature N-terminus remain processing assay gates.

<!-- OPENISLET_CURRENT_LAYER_SCOPE -->
Layer 2 uses the pinned release pairing Ensembl release 112 / GENCODE v46 (GRCh38.p14, 1-based closed coordinates; ENST00000309268.11 reconciled) and contains 11 literature-anchored stage-specific prior profiles with 1,000 independent prior-scenario draws per endpoint; it is not a fitted dynamic GRN/ODE or held-out validated differentiation predictor. Arm 1B constitutive EEF1A1 drives antisense CD47 and inducible sense 5xISRE drives HLA-E. Computed Arm 1A/1B Hour-0 CD47 is 209,444/238,004 molecules/cell; Arm 1B/WT is 4.96679814x. Resting Arm 1B HLA-E is 10,890 and IFNg24h HLA-E 69,516; the resting endpoint is below the assumed 15,000 threshold. Exact edited beta/alpha/delta/iEC density is unmeasured, with a preregistered Stage-0 acceptance target for Hour-0 density. TYK2-KO Stage-6 prior is 0.364 versus 0.678 WT (-46.3%); this is distinct from the transcript-level endpoint. HNRNPA2B1-proximal 455-bp CpG-core candidate (GRCh38 chr7:26200580-26201034); anti-silencing function requires head-to-head lineage qualification against the published CBX3-side fragment. The local 3,118-row AVI table is quarantined in the private cache: end-to-end batch generation provenance for the full 3,118-row table is incomplete (as a whole it is not a certified model scan), a limited 28 Oct 2026 live re-run against alphagenome v0.10.0 verified 59/59 sampled rows (9 candidate edits plus 50 random rows), AlphaGenome non-commercial output terms conflict with CC-BY-4.0 public redistribution, and zero AlphaGenome outputs are shipped publicly or used as active ranking weights or insertion-safety certifications. Scenario-B P2A profiles are unsupported provisional comparators; the exact 2A-free Scenario-B donor/prediction remains unavailable.

Reference-geometry compatibility screen: 3CDG is a CD94/NKG2A complex. P8-Phe8 contributes 24 heavy-atom contacts to shared CD94 and 0 to NKG2A in copy 1, and 21/0 in copy 2. Removing Phe ring atoms is an atom-truncation diagnostic, not a relaxed Val8 model. This structure cannot discriminate NKG2A from NKG2C. Receptor-specific K_D values are working soluble-ectodomain reference priors; direct SPR and HCMV+ NKG2C+ NK-donor assays remain preregistered. Reference SIRPA V1/V2 sequence and geometry comparison; pan-allelic protection is unverified. Asp101 is inserted in SIRPA V1 (4CMM) relative to V2 (2JJS). Soluble crystal CD47 carries C33G/Cys15-to-Gly; engineered CD47 retains WT Cys15. Reference soluble-C15G structures place CD47 Val115 about 30.93 A from the nearest SIRPA atom. This coarse anchor-distance screen does not model membrane orientation, the engineered WT-Cys15/GPI stalk or V8; binding, membrane and phagocytosis assays remain pending. The 21-residue linker-2 path includes unresolved HLA-E Gly1 between the 20-aa (G4S)4 linker and HLA-E Ser2 N; 73.50 A contour, 0.4234 extension, 9.901 pN are coarse WLC quantities. The 3CDG reference-crystal Arg107 C-alpha distance is 17.02 A, nearest-heavy-atom distance 14.82 A; neither is a measured engineered Gly107 geometry. Foldseek is a cached-hit replay with missing historical database/tool/date provenance, not a reproducible original search. 8H4X was superseded by 8ZCR (Yan, Zhou et al., 2024, Acta Crystallogr F 80:286-293; DOI:10.1107/S2053230X24009439; PMID:39382088). Exact-product function remains unresolved.

Layer 4 is a conditional model, not validated graft survival. The physical source lamina half-width fraction is 0.34; its continuous central volume is derived from that same geometry and control-volume overlap at every mesh. Nominal minimum slab oxygen is 19.92 mmHg at 18 h. The 8.64-mmHg homogenized fixed-surface Dirichlet core is an upper-bound screen; operative uniform-supply / conservative gel-only two-zone cores are 4.08/2.15 mmHg. Physical supply topology remains a Stage-2 assay gate; mesh error is evaluated on unrounded values, with NUMERICALLY_UNRESOLVED whenever an uncertainty interval straddles the threshold. Reference fixed-functional-beta-IEQ ODE curves have pde_coupling_applied=false and cannot inherit a PDE-derived full-product verdict; matched PDE/ODE states are identified separately. Reference ODE outputs: OPENISLET_60_30_10_TRIHORMONAL_G6PC2_WT: minimum 4.59 mM / insulin shutoff t50 4.93 min; BETA_ONLY_MONOCULTURE_NO_ALPHA_DELTA: minimum 3.06 mM / insulin shutoff t50 9.24 min; G6PC2_KO_DISRUPTED_GLYCOLYTIC_RHEOSTAT_GY002: minimum 2.70 mM / insulin shutoff t50 10.44 min. Low-glucose challenge reference-ODE minima: lead 3.9218, composite-reflex OAT1 3.8805, pure-beta-clamp OAT1A 3.9186 mM; these share the specified prior kinetic constants and are not measured clinical protection. Nominal modeled peak peroxide is 0.725 uM and maximum pH 7.405.

Layer 5 is an illustrative assigned-rule population simulator, not a clinical forecast. Personalized analytic expectation is 97.5000%, standardized Plan A analytic 92.5000%; sampled personalized 97.40%. Analytic ancestry gap is 0.40 percentage points; sampled gap 1.80 points, with finite sample/seed scope. Genome-wide off-target enumeration is NOT_ENUMERATED and coding-hit count is null; assigned illustrative three-mismatch panel scores are not exhaustive genome/exome scans or empirical HiFi-SpCas9 specificity. Variant identity is source-verified dbSNP interval identity; rare-variant screening frequencies are explicit assumed priors, with gnomAD row-level AC/AN validation pending. Manufacturing yield is null/NOT_MEASURED; exp(-lambda_K*max(0,K-4))*gamma_B**B is a dimensionless engineering utility discount. The executable optional TYK2 design is PEmax+epegRNA. CGBE remains an unimplemented contingency requiring an engineered SpG/SpRY-compatible PAM, editing-window/bystander ledger and separate validation; no CGBE guide is qualified here.
<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->
