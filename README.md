<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset=".github/openislet-mark-dark.svg">
    <img src=".github/openislet-mark.svg" width="80" height="80" alt="OpenIslet">
  </picture>
</p>

<h1 align="center">OpenIslet blueprint</h1>

<p align="center">The open design for a personalised cure for type 1 diabetes: the evidence, the Idea Graveyard, and the checks that score every design.</p>

<p align="center"><a href="https://openislet.org">Website</a> · <a href="https://openislet.org/blueprint">Blueprint</a> · <a href="https://openislet.org/start">Take part</a> · <a href="https://github.com/openislet/website">Website repository</a></p>

**OpenIslet** (`https://openislet.org`) is an open-source research project using frontier AI and deterministic computational verification to design and stress-test **Islet-1**: a candidate blueprint for a stem-cell-derived 3D islet mini-organ (`SC-islet`) aimed at restoring endogenous insulin and counter-regulatory glucagon/somatostatin control in Type 1 Diabetes (`T1D`) without chronic systemic anti-rejection drugs.

> [!IMPORTANT]
> **What Islet-1 Is — and What It Is Not:**
> - **What it is:** A pre-registered, falsifiable 5-layer computational and preclinical research blueprint (`v0.8.0` open draft, launched `2026-11-01`), backed by `563` extracted studies (`548` active deduplicated, `3,093 atomic effects`), `10` machine-actionable Idea Graveyard rules, base-resolution construct sequences, 3D structural re-verifications, reaction-diffusion/metabolic simulations, and a 36-run fractional-factorial wet-lab protocol (`1,351` pinned files in [`manifest.json`](manifest.json)).
> - **What it is not:** Islet-1 is **not** medical or clinical advice, **not** an approved or clinically tested therapy, and **not** an empirically fitted clinical survival forecast (`outward_rankings_frozen: false`). All Layer 1–5 outputs are prior-anchored scenario evaluations (`99.8542%` structural prior share of component variance) awaiting Stage 0–4 physical wet-lab qualification.

## Explore the Blueprint

| Resource | Description |
| :--- | :--- |
| **[`ISLET-1-BLUEPRINT.md`](ISLET-1-BLUEPRINT.md)** | Complete 5-Layer scientific blueprint (`Blueprint A` 4-locus universal lead, `Blueprint B` HLA-C-retained comparator, and `Blueprint X` challenger). |
| **[`docs/STAGE-0-EXECUTION-PLAYBOOK.md`](docs/STAGE-0-EXECUTION-PLAYBOOK.md)** | Preregistered Stage-0 Construct Pre-Pilot (`Decision Experiments #1–#4`) and Stage-1 36-run `HIP-Core` fractional-factorial wet-lab protocol. |
| **[`docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md`](docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md)** | Full statistical design algebra, variance decompositions, nucleotide grammar audits, 3D structural tables, and pipeline verification commands. |
| **[`docs/ISLET-1-CANDIDATE-BLUEPRINTS-AND-OPTIMIZER-v0.2.md`](docs/ISLET-1-CANDIDATE-BLUEPRINTS-AND-OPTIMIZER-v0.2.md)** | Multi-gate scenario evaluator specification, clinical benchmark comparisons, and MCP coordination architecture. |
| **[`docs/ISLET-1-SYNTHESIS-REPORT-v0.4.md`](docs/ISLET-1-SYNTHESIS-REPORT-v0.4.md)** | Multi-model red-team convergence synthesis and verification ledger. |
| **[`data/graveyard/idea_graveyard_v0.5.json`](data/graveyard/idea_graveyard_v0.5.json)** | Executable Idea Graveyard (`10/10` machine-actionable `BLOCK` and `WARN` rules preventing biological dead ends and epistatic traps). |
| **`https://mcp.openislet.org/mcp`** | Hosted Model Context Protocol (`MCP`) connector (launched `2026-11-01` in Mode A: 10 read tools, `check_draft` dry-run validator, and open `IIP` proposal issue queue). |

## Quick Start & Independent Verification

Requires **Python `>= 3.10`** on macOS or Linux (`numpy>=2.0`, `pytest>=7.4`, `jsonschema>=4.18`):

```bash
# 1. Create a virtual environment outside the checkout and install dependencies
python3 -m venv ../openislet-venv
../openislet-venv/bin/pip install -r requirements.txt

# 2. Verify public bundle SHA-256 hashes and file inventory against manifest.json
../openislet-venv/bin/python -B scripts/generate_release_manifest.py --verify-public-bundle

# 3. Run the automated verification and acceptance test suite (layers 1–5, citations, and math)
../openislet-venv/bin/python -B -m pytest -p no:cacheprovider

# 4. Run the standalone 19-check adversarial acceptance and quote-cap suite
../openislet-venv/bin/python -B tests/test_adversarial_acceptance.py

# 5. Dry-run the Stage-1 wet-lab Bayesian calibration pipeline against synthetic test fixtures
../openislet-venv/bin/python -B scripts/ingest_wetlab_results.py --dry-run-calibration --no-write
```

For full details on what can be verified directly from a public checkout versus which raw-corpus re-extraction steps require local PubMed/PMC dossier caches (`data/dossiers/*.json`, excluded from public redistribution), see Section 5 of [`docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md`](docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md).

## How to Contribute (`IIP` Open Review Queue)

Anyone—researcher, clinician, patient, or AI agent—can attack a claim, report a flaw, add a missing study, or propose a better edit via an **Islet Improvement Proposal (`IIP`)**:
1. Connect your AI client to `https://mcp.openislet.org/mcp` (or inspect this repository locally) and draft a proposal in Chat #1 (`draft_a_proposal`).
2. Run an independent blind review (`[OPENISLET-REVIEW-v1]` containing `"Claim reviewed:"`) in a fresh Chat #2, validate statelessly with `check_draft`, and submit via `submit_proposal` to open a public `[IIP]` issue in the review queue.
3. Accepted contributions—including falsifications that rule out a dead end—are credited in plain words (stating the exact outcome and follow-up), named in the next quarterly locked version (`Islet-1.0` on `1 Feb 2027`), and recorded on `https://openislet.org/contributions`.

## Licensing

- **Code, Scripts, Schemas & Tests (`scripts/`, `schemas/`, `tests/`):** Licensed under the **Apache License 2.0** ([`LICENSE`](LICENSE)).
- **Blueprint Documentation, Structured Annotations, Evidence Graph & Idea Graveyard (`ISLET-1-BLUEPRINT.md`, `README.md`, `docs/`, `data/`, `manifest.json`):** Licensed under **Creative Commons Attribution 4.0 International (`CC-BY-4.0`)** ([`LICENSE-DATA`](LICENSE-DATA)), with third-party verbatim sentence excerpts (`<= 240` characters and `< 25%` of abstract length for non-OA studies) and external database records (`PDB`, `AlphaFold DB`, `Ensembl`, `GTEx`, `ENCODE`, `HPA`, `gnomAD`, `dbSNP`, `ClinVar`, `JASPAR`) governed by their respective upstream terms in [`.licenses/`](.licenses/).
