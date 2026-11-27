#!/usr/bin/env python3
"""Focused standalone verification suite for IIP-0034 (DEC-113).

Covers:
  1. Design Y-FullStack v2 Mode-B FASTA SHA-256, 6 role-scoped donor records, mandatory policy blocks,
     nonfinite JSON rejection, and strict product_contract whole-record full_record_sha256 pinning.
  2. public_excerpt_policy GENES additions and dual redaction prefix separation (CC-licensed vs non-OA).
  3. Registry-loaded compatibility-constrained 2^(4N) multi-paradigm portfolio artifact, 4,096-state
     incumbent parity, W_B0 modular permission disclosure, W_B9/W_B3 sensitivity & adverse-mixture
     erasure weights, nonfinite H4 mass rejection, and +0.00 pp live invariance.
  4. Standalone cross-runtime planar Route 3P & portfolio verifier (--check byte/SHA-256 parity,
     DEC-113 decision-tree SHA-256 & CommonMark 0.31.2 Playbook s4.1.19 SHA-256, 96-cell typed
     product-by-gate contract & 30-endpoint DEC-102 crosswalk, missing FASTA & registry mutations,
     wetlab-free minimal replay, and CRLF/value tamper rejection).
  5. Ring-1 corpus expansion (589 -> 594 studies, 693 arms, 62 reference-control arms, 260 raw edits /
     262 expanded gene occurrences, 3,158 effects, 78 curated landmarks), Sato 179-char Fig. 3 quote &
     quantitative eff3, Izadifar beta-iPLC line 2 & Endogenous_Locus_KI_Without_KO, Shi 3-way split
     controls & widened eff2 span, Liu 3 arms / HEK293T / Fig. 3f, Levesque shared TRAC on both arms &
     3.3/2.2 SEM, and 39571901 preservation.
  6. Context-Evidence Graph expansion (1,172 -> 1,176 studies, 7,740 arms, 14,167 outcomes, 2,483 edges,
     45 curated ineligible edges, typed intervention/comparator arm ID lists, secondary source pins,
     and M1522:L4 provenance disclosure).
  7. Extraction schema enum additions, case-insensitive/padded MULTI_LOCUS_SAFE_HARBOR filtering,
     and 0/2,160 Layer-1 stack invariance.
  8. Reconciled Stage-0 decision-tree + Playbook s4.1.19 + product registry, mathematical parity with
     incumbent DEC-103 wetlab_contract helper on valid zero-event cutoffs, public-compatible cross-study
     graph build, and 0-byte invariance of wetlab_endpoint_contract_v0.8.0.json (366 / 77 / 122).
  9. CommonMark 0.31.2 Playbook heading enforcement (fenced 4.1.19, multiline Setext, blockquoted/list
     duplicates rejected; literal code/HTML-comment examples ignored), python -O integrity enforcement,
     mutually exclusive CLI flags (--write --check rejected), and contradiction guards across all 3 docs.
 10. Exact 50-digit conservative Decimal R3P-05 boundary enforcement: binary64 ULP straddle (n=5,991,464,
     cap=4.999999206572849e-7 -> HOLD), huge n=10^308 non-zero scientific string, single-dose Blueprint
     93,350 counterexample (3.20909e-5 / 6.41818e-5 -> HOLD), NumPy 1-elem/2-elem array rejection,
     control/zero-width scope character rejection, mandatory f_clonogenic_provenance & 3-way Bonferroni
     alpha allocation, and unconditional k >= 1 FAIL_HARD_STOP with separate attribution_status.
 11. Strict CEG referential ownership (cross-study arm/comparator/edge joins, false outcomes_count, bare
     A999, wrong embedded outcome_id, unresolved candidate_effector_topics), partial-dossier refusal,
     explicit verification_mode reporting, and recursive rolling-40-char + union-shingle _SourceGuard.
 12. FVM geometric/grid guards & actual mass-balance residuals, converged 5,000-mode vs 45-mode Fourier
     thresholds & 0.235 s spatial H200 margin, mechanics critical adhesion & 2-ribbon warning, and
     exhaustive 451-colony integer minimization across 11 qualifying splits.
"""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

EXPECTED_DT_SHA256 = "46609af00da04f998d26b8ffe47f43350587989ef40e9be843ba6849198c5952"
EXPECTED_PB_SHA256 = "7511695ef9dbade1f9a231e3ee522d6c8b749ab8a3c0e6962ac1271a147d1004"
EXPECTED_PARETO_SHA256 = "591ddeb9a706c0bc95fc509765b512aadad8a403287a796e63dd52e7db548575"


def _load_mod(rel: str, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _build_minimal_verifier_tree(dst: Path) -> Path:
    (dst / "scripts").mkdir(parents=True)
    (dst / "data" / "extracted").mkdir(parents=True)
    (dst / "data" / "protocols").mkdir(parents=True)
    (dst / "data" / "constructs").mkdir(parents=True)
    (dst / "docs").mkdir(parents=True)
    for rel_s in (
        "scripts/verify_iip0034_planar_and_portfolio.py",
        "scripts/product_contract.py",
        "scripts/engineering_priors.py",
        "scripts/public_excerpt_policy.py",
        "scripts/rank_islet1_stacks.py",
        "scripts/run_layer4_tissue_and_metabolic_twin.py",
    ):
        shutil.copy2(ROOT / rel_s, dst / rel_s)
    for rel_d in (
        "data/extracted/iip0034_planar_and_portfolio_verification_v0.8.0.json",
        "data/extracted/islet1_multiparadigm_portfolio_v0.8.0.json",
        "data/extracted/islet1_combinatorial_pareto_and_ladder_v0.8.0.json",
        "data/protocols/islet1_product_definitions_v0.8.0.json",
        "data/protocols/stage0_decision_tree_v1.json",
        "docs/STAGE-0-EXECUTION-PLAYBOOK.md",
    ):
        shutil.copy2(ROOT / rel_d, dst / rel_d)
    for fa in (ROOT / "data" / "constructs").glob("*.fasta"):
        shutil.copy2(fa, dst / "data" / "constructs" / fa.name)
    return dst


def test_01_design_y_v2_fasta_sha256_and_strict_product_contract(tmp_path: Path) -> None:
    pc = _load_mod("scripts/product_contract.py", "openislet_pc_iip0034_test")
    fasta_path = ROOT / pc.DESIGN_Y_V2_FASTA_PATH
    assert fasta_path.is_file()
    raw = fasta_path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == pc.DESIGN_Y_V2_FASTA_SHA256

    defs = pc.load_product_contracts(ROOT, require_staged_fastas=True)
    assert len(defs["products"]) == 38
    by_id = {p["product_variant_id"]: p for p in defs["products"]}
    for pvid in pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS:
        assert pvid in by_id
        p = by_id[pvid]
        assert p["execution_supported"] is False
        assert p["physical_route_id"] == "ROUTE_3P_MODULAR_PLANAR_ENVELOPE"

    bk = by_id["ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR"]
    bk_gates = " ".join(bk["pending_assay_gates"])
    assert "11-TU" not in bk_gates and "CD46_WT" not in bk_gates
    assert "12-TU" in bk_gates and "CD55+CFH" in bk_gates
    assert "3 ectopic GPI" in bk_gates and "3x ectopic" not in bk_gates
    assert "activated only if primary CD46_WT fails Gate 0V MCB viral clearance" in bk["execution_scope"]
    assert "activated only if primary CD46_WT fails Gate 0V MCB viral clearance" in bk["alleles"]["B2M"]
    assert "activated only if primary CD46_WT fails Gate 0V MCB viral clearance" in bk["staged_design_y_fullstack_v2_donor_binding"]["stage0_promotion_gates"]

    ctrl = by_id["ISLET1_Y_FULLSTACK_V2_CTRL_EQUAL_SPACING"]
    ctrl_gates = " ".join(ctrl["pending_assay_gates"])
    assert "Y-05" in ctrl_gates and "CD59" in ctrl_gates and "CD47-IgV-GPI" in ctrl_gates and "180,000" in ctrl_gates
    assert "1,523 bp longer" in ctrl_gates

    dec053 = by_id["ISLET1_Y_FULLSTACK_V2_DEC053_CLYBL_CDK1_HTMPK"]
    assert "Not encoded by this standalone chr13 donor" in dec053["hla_retention"]["classical_MHC_I"]

    for excl_id in (
        "ISLET1_Y_FULLSTACK_V2_B2M_HALF_A_ATTB_GT",
        "ISLET1_Y_FULLSTACK_V2_B2M_HALF_B_ATTB_CA",
        "ISLET1_Y_FULLSTACK_V2_DEC053_CLYBL_CDK1_HTMPK",
        "ISLET1_Y_FULLSTACK_V2_CTRL_EQUAL_SPACING",
    ):
        assert "EXCLUDED_FROM_FULLSTACK_CLEARANCE" in by_id[excl_id]["staged_design_y_fullstack_v2_donor_binding"]["stage0_promotion_gates"]

    copy_root = tmp_path / "repo"
    (copy_root / "data" / "protocols").mkdir(parents=True)
    (copy_root / "data" / "constructs").mkdir(parents=True)
    orig_prod_path = ROOT / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json"
    p_file = copy_root / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json"
    shutil.copy2(orig_prod_path, p_file)
    for fa in (ROOT / "data" / "constructs").glob("*.fasta"):
        shutil.copy2(fa, copy_root / "data" / "constructs" / fa.name)
    assert len(pc.load_product_contracts(copy_root, require_staged_fastas=True)["products"]) == 38

    # 1. Tamper FASTA
    t_fa = copy_root / pc.DESIGN_Y_V2_FASTA_PATH
    t_fa.write_bytes(t_fa.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)
    shutil.copy2(fasta_path, t_fa)

    # 2. Flip execution_supported=True
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU":
            p["execution_supported"] = True
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 3. Mutate physical_route_id
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU":
            p["physical_route_id"] = "ROUTE_1_MACROENCAPSULATED"
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="ROUTE_3P_MODULAR_PLANAR_ENVELOPE"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 4. Mutate donor_sha256 -> caught by full_record_sha256
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU":
            p["staged_design_y_fullstack_v2_donor_binding"]["donor_sha256"] = "0" * 64
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="full_record_sha256"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 5. Truncate backup pending_assay_gates
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR":
            p["pending_assay_gates"] = ["R3P-06 alone"]
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="full_record_sha256"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 6. Copy primary pending_assay_gates onto spacing control
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    prim_gates = next(p["pending_assay_gates"] for p in p_doc["products"] if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU")
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_CTRL_EQUAL_SPACING":
            p["pending_assay_gates"] = list(prim_gates)
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="full_record_sha256"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 7. Reverse backup activation trigger in execution_scope & alleles.B2M
    for field_mut in ("execution_scope", "alleles"):
        p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
        for p in p_doc["products"]:
            if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR":
                if field_mut == "execution_scope":
                    p["execution_scope"] = p["execution_scope"].replace("fails Gate 0V", "passes Gate 0V")
                else:
                    p["alleles"]["B2M"] = p["alleles"]["B2M"].replace("fails Gate 0V", "passes Gate 0V")
        p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match="full_record_sha256"):
            pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 8. Mandatory policy block removal and root-level NaN rejected
    for drop_pol in ("mcb_donor_blocking_policy", "design_b_four_state_lifecycle_policy"):
        p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
        p_doc.pop(drop_pol, None)
        p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match=drop_pol):
            pc.load_product_contracts(copy_root, require_staged_fastas=True)

    raw_nan = orig_prod_path.read_text(encoding="utf-8").replace("{", '{"injected_nan": NaN, ', 1)
    p_file.write_text(raw_nan, encoding="utf-8")
    with pytest.raises(ValueError, match="Non-standard JSON constant"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 9. Legal JSON exponent overflowing binary64 (1e400 -> inf) rejected in product_contract AND verifier _load_json_strict
    raw_ovf = orig_prod_path.read_text(encoding="utf-8").replace("{", '{"overflow_test": 1e400, ', 1)
    p_file.write_text(raw_ovf, encoding="utf-8")
    with pytest.raises(ValueError, match="Non-finite JSON float overflow"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)
    vmod_ovf = _load_mod("scripts/verify_iip0034_planar_and_portfolio.py", "openislet_v34_ovf_test")
    with pytest.raises(ValueError, match="Non-finite JSON float overflow"):
        vmod_ovf._load_json_strict(p_file)


def test_02_public_excerpt_policy_genes_and_dual_redaction_prefix() -> None:
    pep = _load_mod("scripts/public_excerpt_policy.py", "openislet_pep_iip0034_test")
    for g in ("VTCN1", "AIFM2", "CDK1", "DTYMK"):
        assert g in pep.GENES
    assert "Endogenous_Locus_KI_Without_KO" in pep.MODALITIES
    assert pep.PUBLISHER_LICENSE_REGISTRY["17375075"] == "CC_BY_NC_ND_4_0"
    assert pep.redaction_token_prefix_for_pmid("17375075") == "[REDACTED_CONSERVATIVE_EXCERPT:"
    assert pep.redaction_token_prefix_for_pmid("42020736") == "[REDACTED_NON_OA_EXCERPT:"
    h64 = "a" * 64
    assert pep.validate_token(f"[REDACTED_CONSERVATIVE_EXCERPT: see PMID:17375075 field=verbatim_evidence_quote quote_sha256={h64}]", expected_pmid="17375075", expected_field="verbatim_evidence_quote")
    assert pep.validate_token(f"[REDACTED_NON_OA_EXCERPT: see PMID:42020736 field=verbatim_evidence_quote quote_sha256={h64}]", expected_pmid="42020736", expected_field="verbatim_evidence_quote")


def test_03_compatibility_constrained_2pow4n_portfolio_and_4096_parity() -> None:
    rank_mod = _load_mod("scripts/rank_islet1_stacks.py", "openislet_rank_iip0034_test")
    mp_path = ROOT / "data" / "extracted" / "islet1_multiparadigm_portfolio_v0.8.0.json"
    assert mp_path.is_file()
    shipped = json.loads(mp_path.read_text(encoding="utf-8"))
    recomputed = rank_mod.compute_multiparadigm_2pow4n_portfolio_iip0034()
    assert shipped == recomputed

    assert shipped["live_headline_delta_pp_today"] == 0.0
    assert "R3P" not in rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["UNHARDENED_A_B_X"]
    assert rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["DESIGN_Y_FULLSTACK_V2"] == frozenset({"R3P"})
    assert "54.9679%" in shipped.get("w_b0_modular_permissions_note", "")

    all_rows = shipped["conservative_zero_h2_h3_lift_4x3_architectural"] + shipped["stage0_gated_conditional_4x3_off_the_shelf"]
    for row in all_rows:
        assert float(row["incumbent_4096_state_max_abs_err"]) <= 1e-12
    by_sc = {r["scenario"]: r for r in all_rows}
    w_b0 = by_sc["W_B0_4Arm_Y@R2_ZeroLift_phiYX0.75 (No Route 3P)"]
    assert w_b0["stage0_qualification_status"] == "HYPOTHETICAL_UNREGISTERED_LEGACY_R2_COUNTERFACTUAL"
    w_b9 = by_sc["W_B9_4Arm_Y@PrePerf3P_Nested(90/92)_pY4=0.68_c4u0.90_phi1.0_50"]
    assert abs(float(w_b9["locked_pct"]) - 57.9434) < 1e-3
    assert abs(float(w_b9["modular_k100_pct"]) - 59.2796) < 1e-3
    sens_phi = shipped["phi_h4_sensitivity_counterexamples"][0]
    assert abs(float(sens_phi["locked_pct"]) - 57.7545) < 1e-3
    sens_marg = shipped["phi_h4_sensitivity_counterexamples"][1]
    assert abs(float(sens_marg["locked_pct"]) - 56.8286) < 1e-3
    assert sens_marg["clears_57p92_locked"] is False
    neg_sens = shipped["cross_hurdle_negative_coupling_sensitivity"]
    assert abs(float(neg_sens[0]["negative_cross_hurdle_locked_pct"]) - 57.4151) < 1e-3 and neg_sens[0]["clears_57p92_locked"] is False
    assert abs(float(neg_sens[0]["adverse_mix_weight_erasing_gain_pct"]) - 2.721) < 0.02
    assert abs(float(neg_sens[1]["negative_cross_hurdle_locked_pct"]) - 57.2375) < 1e-3 and neg_sens[1]["clears_57p92_locked"] is False
    assert abs(float(neg_sens[1]["adverse_mix_weight_erasing_gain_pct"]) - 0.365) < 0.01
    assert shipped["five_arm_bridge_registration_status"]["D1B_AUTO_C_V2"] == "UNREGISTERED_STAGE0_FORWARD_BRIDGE_PLACEHOLDER_NO_IIP0034_PRODUCT_OR_FASTA"

    # Nonfinite H4 probability mass must raise ValueError (never silently return locked_pct=nan)
    bad_dists = copy.deepcopy(list(rank_mod._live_canonical_incumbent_dists_iip0034()))
    if isinstance(bad_dists[3], dict):
        k0 = next(iter(bad_dists[3]))
        bad_dists[3][k0] = float("nan")
    else:
        bad_dists[3] = [float("nan")] + list(bad_dists[3])[1:]
    with pytest.raises(ValueError, match="non-finite|Non-finite|finite|mass"):
        rank_mod._run_dp_and_verify_parity_iip0034(
            "NONFINITE_TEST",
            [(0.88, 0.78, 0.57, 0.69), (0.88, 0.79, 0.64, 0.70), (0.88, 0.82, 0.73, 0.69)],
            bad_dists,
            [0.88, 0.87, 0.89, 0.85],
            0.50778198,
            0.54116275,
            0.57917640,
            stage0_status="SELF_CHECK_NONFINITE",
        )

    pareto_sha = hashlib.sha256(
        (ROOT / "data" / "extracted" / "islet1_combinatorial_pareto_and_ladder_v0.8.0.json").read_bytes()
    ).hexdigest()
    assert pareto_sha == EXPECTED_PARETO_SHA256


def test_04_planar_and_portfolio_verifier_replay_and_adversarial_controls(tmp_path: Path) -> None:
    v_script = ROOT / "scripts" / "verify_iip0034_planar_and_portfolio.py"
    art_path = ROOT / "data" / "extracted" / "iip0034_planar_and_portfolio_verification_v0.8.0.json"
    assert v_script.is_file() and art_path.is_file()

    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"

    res = subprocess.run(
        [sys.executable, "-B", str(v_script), "--check"],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "OK:" in res.stdout

    doc = json.loads(art_path.read_text(encoding="utf-8"))
    assert doc["proposal_id"] == "IIP-0034"
    assert doc["decision_id"] == "DEC-113"
    assert doc["live_headline_odds_delta_pp"] == 0.0
    reg_aud = doc["design_y_v2_product_registry_audit"]
    assert reg_aud["registry_loaded_via_product_contract"] is True
    assert reg_aud["dec113_decision_tree_sha256_verified"] == EXPECTED_DT_SHA256
    assert reg_aud["dec113_playbook_sec4119_sha256_verified"] == EXPECTED_PB_SHA256
    assert reg_aud["resolved_design_y_v2_routes_from_registry"] == ["R3P"]
    typed_aud = reg_aud["typed_product_by_gate_contract_audit"]
    assert typed_aud["products_checked"] == 6 and typed_aud["gates_per_product"] == 16
    assert typed_aud["cells_verified"] == 96 and typed_aud["dec102_endpoints_mapped"] == 30
    assert typed_aud["contradiction_guards_passed"] == 4 and typed_aud["typed_predicates_verified"] == 24
    assert doc["portfolio_compatibility_and_4096_state_audit"]["unhardened_abx_barred_from_route3p"] is True
    assert doc["portfolio_compatibility_and_4096_state_audit"]["design_y_v2_bound_strictly_to_route3p_matching_registry"] is True
    assert doc["portfolio_compatibility_and_4096_state_audit"]["max_4096_state_incumbent_abs_err"] <= 1e-12
    assert abs(float(doc["portfolio_compatibility_and_4096_state_audit"]["w_b3_equalized_y_marginals_sensitivity_counterexample"]["locked_pct"]) - 56.8286) < 1e-3

    # Minimal wetlab-free tree: verifies zero wetlab dependency AND proves BOTH verifier --check AND ranker producer
    # fail on missing registry, missing FASTA, route/alleles.B2M/gate/scope mutation, or missing DEC-113 block.
    min_root = _build_minimal_verifier_tree(tmp_path / "min_verifier_tree")
    assert not (min_root / "scripts" / "wetlab_contract.py").exists()

    r_min = subprocess.run(
        [sys.executable, "-B", str(min_root / "scripts" / "verify_iip0034_planar_and_portfolio.py"), "--check"],
        cwd=str(min_root),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "OK:" in r_min.stdout

    min_prod_path = min_root / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json"
    orig_prod_bytes = min_prod_path.read_bytes()

    def _assert_check_and_ranker_fail() -> None:
        r_v = subprocess.run(
            [sys.executable, "-B", str(min_root / "scripts" / "verify_iip0034_planar_and_portfolio.py"), "--check"],
            cwd=str(min_root), env=env, capture_output=True, text=True,
        )
        assert r_v.returncode != 0
        r_rk = subprocess.run(
            [sys.executable, "-B", "-c", "import sys; sys.path.insert(0, 'scripts'); import rank_islet1_stacks as r; r.compute_multiparadigm_2pow4n_portfolio_iip0034()"],
            cwd=str(min_root), env=env, capture_output=True, text=True,
        )
        assert r_rk.returncode != 0

    for mutate_fn in (
        lambda d: [p.__setitem__("physical_route_id", "ROUTE_1_MACROENCAPSULATED") for p in d["products"] if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU"],
        lambda d: [p.__setitem__("pending_assay_gates", ["R3P-06 alone"]) for p in d["products"] if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR"],
        lambda d: [p.__setitem__("execution_scope", p["execution_scope"].replace("fails Gate 0V", "passes Gate 0V")) for p in d["products"] if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR"],
        lambda d: [p["alleles"].__setitem__("B2M", p["alleles"]["B2M"].replace("fails Gate 0V", "passes Gate 0V")) for p in d["products"] if p["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR"],
    ):
        mut_doc = json.loads(orig_prod_bytes.decode("utf-8"))
        mutate_fn(mut_doc)
        min_prod_path.write_text(json.dumps(mut_doc, indent=2) + "\n", encoding="utf-8")
        _assert_check_and_ranker_fail()
    min_prod_path.write_bytes(orig_prod_bytes)

    # Missing registry & missing FASTA fail both verifier and ranker
    min_prod_path.unlink()
    _assert_check_and_ranker_fail()
    min_prod_path.write_bytes(orig_prod_bytes)

    min_fa_path = min_root / "data" / "constructs" / "islet1_design_y_fullstack_v2_v0.8.0.fasta"
    orig_fa_bytes = min_fa_path.read_bytes()
    min_fa_path.unlink()
    _assert_check_and_ranker_fail()
    min_fa_path.write_bytes(orig_fa_bytes)

    # Strict byte/SHA-256 check rejects CRLF verification JSON, trailing-newline portfolio JSON, and value tampers
    vmod = _load_mod("scripts/verify_iip0034_planar_and_portfolio.py", "openislet_v34_tamper_test")
    crlf_art = tmp_path / "iip0034_crlf.json"
    crlf_art.write_bytes(art_path.read_bytes().replace(b"\n", b"\r\n"))
    old_out, old_p, old_argv = vmod.DEFAULT_OUT_PATH, vmod.PORTFOLIO_PATH, sys.argv
    try:
        vmod.DEFAULT_OUT_PATH = crlf_art
        sys.argv = [str(v_script), "--check"]
        with pytest.raises(SystemExit, match="byte/SHA-256 drift"):
            vmod.main()
    finally:
        vmod.DEFAULT_OUT_PATH = old_out
        sys.argv = old_argv

    padded_mp = tmp_path / "islet1_multiparadigm_portfolio_padded.json"
    padded_mp.write_bytes((ROOT / "data" / "extracted" / "islet1_multiparadigm_portfolio_v0.8.0.json").read_bytes() + b"\n")
    try:
        vmod.PORTFOLIO_PATH = padded_mp
        sys.argv = [str(v_script), "--check"]
        with pytest.raises(SystemExit, match="byte/SHA-256 drift"):
            vmod.main()
    finally:
        vmod.PORTFOLIO_PATH = old_p
        sys.argv = old_argv


def test_05_ring1_corpus_5_new_studies_and_39571901_preservation() -> None:
    ext_files = sorted((ROOT / "data" / "extracted").glob("[0-9]*.json"))
    assert len(ext_files) == 604
    docs = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in ext_files}
    for pmid in ("17375075", "39696686", "42020736", "42056515", "42749803"):
        assert pmid in docs, f"Missing new Ring-1 study {pmid}"
    assert sum(len(d.get("arms", [])) for d in docs.values()) == 708
    assert sum(1 for d in docs.values() for a in d.get("arms", []) if a.get("is_reference_control") is True) == 65
    assert sum(len(a.get("gate2_3_4_genetic_and_epitope_edits", [])) for d in docs.values() for a in d.get("arms", [])) == 261
    assert sum(len(d.get("atomic_effects", [])) for d in docs.values()) == 3174
    assert "39571901" in docs

    # Sato: tight <=350-char Fig. 3 span, quantitative eff3 fold-change, AZT challenge type, null zygosity
    sato_e1, _, sato_e3 = docs["17375075"]["atomic_effects"]
    sato_arms_by_id = {a["arm_id"]: a for a in docs["17375075"]["arms"]}
    assert sato_e1["source_span"] == "fulltext_data_sections.figure_legends[0].caption[363:542]|norm=exact"
    assert len(sato_e1["verbatim_evidence_quote"]) <= 350
    assert sato_e1["treatment_n"] == 3 and sato_e1["control_n"] == 3
    assert sato_arms_by_id[sato_e1["comparator_arm_id"]]["is_reference_control"] is True
    assert sato_arms_by_id[sato_e3["comparator_arm_id"]]["is_reference_control"] is True
    assert len(sato_arms_by_id[sato_e3["comparator_arm_id"]]["gate2_3_4_genetic_and_epitope_edits"]) == 1
    assert abs(float(sato_e3["treatment_mean"]) - 4.0) < 1e-6
    assert abs(float(sato_e3["treatment_sd_or_sem"]) - 0.3) < 1e-6
    assert abs(float(sato_e3["log_response_ratio_or_smd"]) - 1.3863) < 1e-3
    assert sato_e3["statistic_type"] == "log_response_ratio_fold_change"
    for a in docs["17375075"]["arms"]:
        for ed in a.get("gate2_3_4_genetic_and_epitope_edits", []):
            assert ed["zygosity_or_allelic_state"] is None
    assert docs["17375075"]["arms"][0]["challenge_conditions"][0]["challenge_type"] == "SmallMolecule_Prodrug_Ablation_Challenge"

    # Izadifar: mouse beta-iPLC line 2, Endogenous_Locus_KI_Without_KO, null zygosity, clean arm3, descriptive eff1
    iza_arms = docs["39696686"]["arms"]
    assert iza_arms[0]["gate1_substrate_and_differentiation"]["cell_line_id"].startswith("beta-iPLC line 2")
    assert "triple-transgenic" in iza_arms[1]["gate1_substrate_and_differentiation"]["cell_line_id"].lower()
    assert len(iza_arms[1].get("gate2_3_4_genetic_and_epitope_edits", [])) == 0
    if not str(iza_arms[1]["arm_label"]).startswith("[REDACTED_"):
        assert "triple-transgenic" in iza_arms[1]["arm_label"].lower() and "unedited" not in iza_arms[1]["arm_label"].lower()
    assert iza_arms[2]["challenge_conditions"][0]["glucose_concentration_mM"] is None
    for idx_ed in (0, 2, 3, 4):
        eds = iza_arms[idx_ed]["gate2_3_4_genetic_and_epitope_edits"]
        assert len(eds) == 1 and eds[0]["ensembl_gene_id"] == "ENSMUSG00000019942"
        assert eds[0]["edit_modality"] == "Endogenous_Locus_KI_Without_KO"
        assert eds[0]["zygosity_or_allelic_state"] is None
    iza_e1, iza_e2 = docs["39696686"]["atomic_effects"]
    assert iza_e1["effect_direction"] == "Descriptive_or_Unoriented"
    assert iza_e1["endpoint_polarity_registry"]["orientation_sign"] == 0
    assert iza_e2["measurement_unit"] == "ng_per_mL_per_1e6_cells"

    # Shi: 4 arms (3 explicit comparator controls matching CEG) and widened eff2 span
    assert len(docs["42020736"]["arms"]) == 4
    shi_e1, shi_e2 = docs["42020736"]["atomic_effects"]
    assert shi_e1["source_span"] == "abstract[860:1112]|norm=exact"
    assert shi_e2["source_span"] == "abstract[693:859]|norm=exact"
    assert shi_e2["endpoint_polarity_registry"]["orientation_sign"] == 0

    # Liu: 3 arms (unused arm4 removed), HEK293T, Fig. 3f descriptive polarity 0, eff1 SD
    assert len(docs["42056515"]["arms"]) == 3
    assert docs["42056515"]["arms"][0]["gate1_substrate_and_differentiation"]["cell_line_id"] == "HEK293T"
    assert docs["42056515"]["arms"][0]["gate1_substrate_and_differentiation"]["cell_source_class"] == "Human_NonIslet_Immortalized_Cell_Line"
    liu_e1, liu_e2 = docs["42056515"]["atomic_effects"]
    assert liu_e1["source_span"] == "fulltext_data_sections.figure_legends[0].caption[52:355]|norm=exact"
    assert liu_e1["dispersion_type"] == "SD"
    assert "TRAC" not in str(liu_e1.get("verbatim_evidence_quote", ""))
    assert liu_e2["measurement_unit"] == "qualitative_gel_assembly_feasibility"
    assert liu_e2["endpoint_polarity_registry"]["orientation_sign"] == 0
    assert (
        liu_e2["source_figure_or_table"].startswith("Fig. 3f")
        or liu_e2["source_figure_or_table"].startswith("[REDACTED_")
    ) and "3f-g" not in json.dumps(liu_e2)

    # Levesque: shared TRAC PrimeEdit on BOTH resting & activated arms, quantitative 3.3 vs 2.2 SEM
    lev_arms = docs["42749803"]["arms"]
    assert len(lev_arms) == 2
    for arm in lev_arms:
        eds = arm["gate2_3_4_genetic_and_epitope_edits"]
        assert len(eds) == 1 and eds[0]["gene_symbol_hgnc"] == "TRAC" and eds[0]["ensembl_gene_id"] == "ENSG00000277734"
    lev_e1 = docs["42749803"]["atomic_effects"][0]
    assert (
        lev_e1["focal_factor_X"] == "Prime_Assembly_TRAC_ssDNA_Primary_Human"
        or (str(lev_e1["focal_factor_X"]).startswith("[REDACTED_") and "genes=TRAC" in str(lev_e1["focal_factor_X"]))
    )
    assert abs(float(lev_e1["treatment_mean"]) - 3.3) < 1e-6
    assert abs(float(lev_e1["control_mean"]) - 2.2) < 1e-6
    assert lev_e1["dispersion_type"] == "SEM"
    assert lev_e1["effect_direction"] == "Descriptive_or_Unoriented"
    assert lev_e1["endpoint_polarity_registry"]["orientation_sign"] == 0


def test_06_ceg_expansion_1176_and_45_curated_ineligible_edges() -> None:
    ceg = json.loads((ROOT / "data" / "extracted" / "context_evidence_graph_v1.json").read_text(encoding="utf-8"))
    assert len(ceg["studies"]) == 1182
    assert len(ceg["arms"]) == 7749
    assert len(ceg["outcomes"]) == 14175
    assert len(ceg["edges"]) == 2494
    elig = sum(1 for e in ceg["edges"] if e.get("primary_certification_eligible") is True)
    inelig = sum(1 for e in ceg["edges"] if e.get("primary_certification_eligible") is False)
    assert (elig, inelig) == (2273, 221)

    assert "primary mouse neurons" not in json.dumps(ceg["studies"]["PMID:42020736"]).lower()
    assert ceg["outcomes"]["PMID:17375075:O1"]["source_pin"]["span"] == "fulltext_data_sections.figure_legends[0].caption[363:542]|norm=exact"
    assert len(ceg["outcomes"]["PMID:17375075:O1"].get("secondary_source_pins", [])) == 2
    assert len(ceg["outcomes"]["PMID:42056515:O1"].get("secondary_source_pins", [])) >= 1
    assert ceg["outcomes"]["PMID:42749803:O1"]["intervention_value"].startswith("Mean 3.3% resting")
    assert "tier_qualification_note" in ceg["studies"]["PMID:42749803"]

    for oid in (
        "PMID:17375075:O1", "PMID:17375075:O2", "PMID:17375075:O3",
        "PMID:42020736:O1",
        "PMID:42056515:O1", "PMID:42056515:O2",
        "PMID:42749803:O1", "PMID:42749803:O2",
    ):
        o = ceg["outcomes"][oid]
        assert isinstance(o.get("intervention_arm_ids"), list) and len(o["intervention_arm_ids"]) >= 1
        assert isinstance(o.get("comparator_arm_ids"), list)

    edges_by_id = {e["edge_id"]: e for e in ceg["edges"]}
    for eid in ("PMID:17375075:L1", "PMID:42020736:L1", "PMID:42056515:L1", "PMID:42749803:L1", "M1522:L4"):
        assert eid in edges_by_id and edges_by_id[eid]["primary_certification_eligible"] is False
    assert "M1522:O7" in edges_by_id["M1522:L4"]["cited_outcome_tuples"]
    assert "provenance_note" in edges_by_id["M1522:L4"]


def test_07_schema_extensions_and_2160_stack_invariance() -> None:
    schema = json.loads((ROOT / "schemas" / "extraction_schema_v0.2.json").read_text(encoding="utf-8"))
    props = schema["properties"]
    assert "Genome_Integration_and_Cargo_Insertion" in props["atomic_effects"]["items"]["properties"]["outcome_domain"]["enum"]
    assert "in_vitro_monoculture_replicates" in props["arms"]["items"]["properties"]["sample_size_unit"]["enum"]
    assert "in_vitro_monoculture_replicates" in props["atomic_effects"]["items"]["properties"]["n_unit"]["enum"]
    assert "Human_NonIslet_Immortalized_Cell_Line" in props["arms"]["items"]["properties"]["gate1_substrate_and_differentiation"]["properties"]["cell_source_class"]["enum"]
    assert "SmallMolecule_Prodrug_Ablation_Challenge" in props["arms"]["items"]["properties"]["challenge_conditions"]["items"]["properties"]["challenge_type"]["enum"]
    ed_props = props["arms"]["items"]["properties"]["gate2_3_4_genetic_and_epitope_edits"]["items"]["properties"]
    assert "Endogenous_Locus_KI_Without_KO" in ed_props["edit_modality"]["enum"]
    assert "Platform_Genome_Engineering_Delivery_Benchmark" in ed_props["functional_gate"]["enum"]
    assert "Safe_Harbor_or_MultiLocus_Target" in ed_props["promoter_or_locus_context"]["enum"]

    rank_mod = _load_mod("scripts/rank_islet1_stacks.py", "openislet_rank_inv_iip0034_test")
    # Case-insensitive & whitespace-padded sentinel AND redacted gene tokens must produce zero gene tokens
    for sentinel in (
        "MULTI_LOCUS_SAFE_HARBOR",
        "  multi_locus_safe_harbor  ",
        "Multi_Locus_Safe_Harbor",
        "[REDACTED_NON_OA_EXCERPT: see PMID:42056515 field=gene_symbol_hgnc quote_sha256=" + "a" * 64 + "]",
        "[REDACTED_CONSERVATIVE_EXCERPT: see PMID:17375075 field=gene_symbol_hgnc quote_sha256=" + "b" * 64 + "]",
    ):
        assert rank_mod.expand_gene_symbol_tokens(sentinel) == []

    _, studies_all, omics = rank_mod.load_inputs()
    new5 = {"17375075", "39696686", "42020736", "42056515", "42749803"}
    studies_pre = [s for s in studies_all if str(s.get("study_metadata", {}).get("pmid")) not in new5]
    prof_all = rank_mod.build_target_evidence_profiles(studies_all, omics)
    trac_pmids = {str(p) for p in (prof_all.get("TRAC") or {}).get("supporting_pmids", [])}
    assert "42056515" not in trac_pmids
    assert "42749803" in trac_pmids
    for bogus_tok in ("MULTI", "LOCUS", "SAFE", "HARBOR", "MULTI_LOCUS_SAFE_HARBOR"):
        assert bogus_tok not in prof_all, f"Spurious target profile created for {bogus_tok!r}"
    c_post = rank_mod.run_combinatorial_multi_gate_search(prof_all)
    c_pre = rank_mod.run_combinatorial_multi_gate_search(rank_mod.build_target_evidence_profiles(studies_pre, omics))
    assert len(c_post) == len(c_pre) == 2160
    assert sum(1 for a, b in zip(c_post, c_pre) if a != b) == 0


def test_08_stage0_decision_tree_dec103_parity_and_wetlab_invariance(tmp_path: Path) -> None:
    dt = json.loads((ROOT / "data" / "protocols" / "stage0_decision_tree_v1.json").read_text(encoding="utf-8"))
    blk = dt.get("design_y_fullstack_v2_and_route3p_gates_iip_0034")
    assert isinstance(blk, dict)
    assert blk["decision_id"] == "DEC-113" and blk["contribution_id"] == "IIP-0034"
    assert blk["product_by_gate_contract"]["cell_count"] == 96
    assert len(blk["product_by_gate_contract"]["cells"]) == 96
    assert len(blk["dec102_endpoint_crosswalk"]) == 30

    # Verify 1:1 coverage of EVERY DEC-102 check in wetlab_contract.py (when present)
    wc_path = ROOT / "scripts" / "wetlab_contract.py"
    if wc_path.is_file():
        import re as _re
        wc_src = wc_path.read_text(encoding="utf-8")
        m_b = _re.search(r"_dy_bool_keys\s*=\s*\((.*?)\)", wc_src, _re.S)
        m_q = _re.search(r"_dy_quant_keys\s*=\s*\((.*?)\)", wc_src, _re.S)
        assert m_b is not None and m_q is not None
        wc_keys = (
            set(_re.findall(r"'([^']+)'", m_b.group(1)))
            | set(_re.findall(r"'([^']+)'", m_q.group(1)))
            | {"obs_kir_nkg2a_neg_nk_lysis_pct", "obs_gsis_si", "design_y_integrase_modality", "design_y_silencing_modality", "design_y_payload_arm"}
        )
        mapped_keys = {k.strip() for r in blk["dec102_endpoint_crosswalk"] for k in r["wetlab_contract_key"].split(";")}
        assert wc_keys <= mapped_keys, f"Unmapped DEC-102 wetlab_contract keys: {sorted(wc_keys - mapped_keys)}"

    vmod = _load_mod("scripts/verify_iip0034_planar_and_portfolio.py", "openislet_v34_dec103_test")
    wc = _load_mod("scripts/wetlab_contract.py", "openislet_wc_dec103_test")
    valid_scope = {
        "product_variant_id": "ISLET1_Y_FULLSTACK_V2_B2M_11TU",
        "construct_sha256": "e77a358ca5bc55c9e0a3070ab46b115c3314660d3a9be6221ef9e91d698cab4e",
        "donor_seedline_id": "DONOR_SEED_A",
        "mcb_id": "MCB_01",
        "lot_id": "LOT_A",
        "device_assembly_id": "DEV_A",
        "physical_route_id": "ROUTE_3P_MODULAR_PLANAR_ENVELOPE",
    }
    scope_ok = {
        "f_clonogenic": 1.0,
        "f_clonogenic_provenance": "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
        "expected_scope": valid_scope,
        "observed_scope": dict(valid_scope),
        "sampling_independence_model": "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES",
        "clustered_sampling": False,
    }
    for k_v, n_v, eta_v, exp_status in [
        (0, 2_995_730, 1.0, "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        (0, 2_995_731, 1.0, "PASS"),
        (0, 5_991_462, 0.50, "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        (0, 5_991_463, 0.50, "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        (0, 5_991_464, 0.50, "PASS"),
        (1, 5_000_000, 1.0, "FAIL_HARD_STOP"),
    ]:
        r_v = vmod.evaluate_r3p05_dec103_thinning(k=k_v, n_raw=n_v, eta_spike_recovery_l95=eta_v, **scope_ok)
        r_wc = wc.evaluate_rare_event_count_iip0024(
            gate_id="R3P_05_ACUTE_GATE0S_RETRIEVAL",
            k=k_v,
            n=n_v,
            max_u95_fraction=1e-6,
            require_observed_zero=True,
            eta_spike_recovery_l95=eta_v,
            f_clonogenic=1.0,
            require_explicit_calibration=True,
        )
        assert r_v["status"] == r_wc["status"] == exp_status

    # Public-compatible cross-study graph build on disposable copy
    cg_root = tmp_path / "cg_tree"
    (cg_root / "scripts").mkdir(parents=True)
    (cg_root / "schemas").mkdir(parents=True)
    (cg_root / "data" / "extracted").mkdir(parents=True)
    shutil.copy2(ROOT / "scripts" / "build_cross_study_graph.py", cg_root / "scripts" / "build_cross_study_graph.py")
    shutil.copy2(ROOT / "schemas" / "extraction_schema_v0.2.json", cg_root / "schemas" / "extraction_schema_v0.2.json")
    if (ROOT / "scripts" / "public_excerpt_policy.py").is_file():
        shutil.copy2(ROOT / "scripts" / "public_excerpt_policy.py", cg_root / "scripts" / "public_excerpt_policy.py")
    for jp in (ROOT / "data" / "extracted").glob("[0-9]*.json"):
        shutil.copy2(jp, cg_root / "data" / "extracted" / jp.name)
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    res_cg = subprocess.run(
        [sys.executable, "-B", str(cg_root / "scripts" / "build_cross_study_graph.py")],
        cwd=str(cg_root), env=env, capture_output=True, text=True, check=True,
    )
    assert "604 studies" in res_cg.stdout

    # Confirm wetlab_endpoint_contract_v0.8.0.json stays untouched from baseline (366 / 77 / 122)
    ct = json.loads((ROOT / "data" / "protocols" / "wetlab_endpoint_contract_v0.8.0.json").read_text(encoding="utf-8"))
    assert len(ct["stage0_decision_tree_endpoints"]) == 366
    assert len(ct["metadata_keys"]) == 77
    assert len(ct["boolean_policy_keys"]) == 122


def test_09_commonmark_playbook_optimized_cli_and_contradiction_guards(tmp_path: Path) -> None:
    min_root = _build_minimal_verifier_tree(tmp_path / "cm_tree")
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    v_script = min_root / "scripts" / "verify_iip0034_planar_and_portfolio.py"

    # 1. Combined --write --check CLI flags must fail cleanly without rewriting drifted files
    r_both = subprocess.run(
        [sys.executable, "-B", str(v_script), "--write", "--check"],
        cwd=str(min_root), env=env, capture_output=True, text=True,
    )
    assert r_both.returncode != 0

    # 2. Python -O (optimized bytecode mode) must still enforce Pareto SHA-256 and all integrity checks
    pareto_p = min_root / "data" / "extracted" / "islet1_combinatorial_pareto_and_ladder_v0.8.0.json"
    orig_pareto = pareto_p.read_bytes()
    pareto_p.write_bytes(orig_pareto + b"\n")
    for flag in ("--check", "--write"):
        r_opt = subprocess.run(
            [sys.executable, "-B", "-O", str(v_script), flag],
            cwd=str(min_root), env=env, capture_output=True, text=True,
        )
        assert r_opt.returncode != 0, f"python -O {flag} failed to reject tampered Pareto"
    pareto_p.write_bytes(orig_pareto)

    # 3. CommonMark 0.31.2 Playbook attacks:
    pb_p = min_root / "docs" / "STAGE-0-EXECUTION-PLAYBOOK.md"
    orig_pb = pb_p.read_text(encoding="utf-8")

    # 3a. Wrapping real s4.1.19 inside a code fence leaves 0 rendered headings -> REJECTED
    fenced_real = orig_pb.replace("#### 4.1.19", "```markdown\n#### 4.1.19", 1) + "\n```\n"
    pb_p.write_text(fenced_real, encoding="utf-8")
    assert subprocess.run([sys.executable, "-B", str(v_script), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0

    # 3b. Competing rendered headings (multiline Setext, blockquote/list/standard ATX, HTML char refs, backslash escapes) -> REJECTED
    for competing in (
        "\n\n4.1.19 Competing multiline\nSetext heading\n---\n",
        "\n\n> #### 4.1.19 Blockquoted competing heading\n",
        "\n\n- #### 4.1.19 List-contained competing heading\n",
        "\n\n#### 4.1.19 Duplicate ATX heading\n",
        "\n\n### 4.1.19 Level-3 ATX duplicate\n",
        "\n\n#### 4&#46;1&#46;19 HTML decimal char-ref duplicate\n",
        "\n\n#### 4&#x2e;1&#x2e;19 HTML hex char-ref duplicate\n",
        "\n\n#### 4\\.1\\.19 Backslash-escaped dot duplicate\n",
    ):
        pb_p.write_text(orig_pb + competing, encoding="utf-8")
        r_dup = subprocess.run([sys.executable, "-B", str(v_script), "--check"], cwd=str(min_root), env=env, capture_output=True, text=True)
        assert r_dup.returncode != 0, f"Competing heading NOT rejected: {competing!r}"

    # 3c. Non-rendered literal examples outside s4.1.19 (fenced backticks/tildes, indented code, HTML comment) -> IGNORED (PASSES!)
    for ignored_example in (
        "\n\n```markdown\n#### 4.1.19 Literal fenced backtick example\n```\n",
        "\n\n~~~\n#### 4.1.19 Literal fenced tilde example\n~~~\n",
        "\n\n\t#### 4.1.19 Tab-indented code block example\n",
        "\n\n<!--\n#### 4.1.19 HTML comment example\n-->\n",
    ):
        pb_p.write_text(orig_pb + ignored_example, encoding="utf-8")
        r_ign = subprocess.run([sys.executable, "-B", str(v_script), "--check"], cwd=str(min_root), env=env, capture_output=True, text=True)
        assert r_ign.returncode == 0, f"Literal code/comment example falsely rejected: {ignored_example!r}\n{r_ign.stderr}"
    pb_p.write_text(orig_pb, encoding="utf-8")

    # 4. Contradiction guards, duplicate JSON keys, AND re-pinned semantic cell/crosswalk tampers rejected!
    dt_p = min_root / "data" / "protocols" / "stage0_decision_tree_v1.json"
    orig_dt = dt_p.read_text(encoding="utf-8")
    dup_dt = orig_dt.replace('"design_y_fullstack_v2_and_route3p_gates_iip_0034": {', '"design_y_fullstack_v2_and_route3p_gates_iip_0034": {"bogus": 1},\n  "design_y_fullstack_v2_and_route3p_gates_iip_0034": {', 1)
    dt_p.write_text(dup_dt, encoding="utf-8")
    assert subprocess.run([sys.executable, "-B", str(v_script), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0

    for bad_txt in (
        orig_dt.replace("clonogenic escape/survival upper bound U95 <= 1e-8", "iCasp9 kill <= 1e-8"),
        orig_dt.replace("3 ectopic GPI", "3x ectopic GPI"),
    ):
        dt_p.write_text(bad_txt, encoding="utf-8")
        assert subprocess.run([sys.executable, "-B", str(v_script), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0
    dt_p.write_text(orig_dt, encoding="utf-8")

    # Re-pinned semantic tamper tests (simulating a revised document with updated EXPECTED_DEC113_DECISION_TREE_SHA256!)
    vmod_sem = _load_mod("scripts/verify_iip0034_planar_and_portfolio.py", "openislet_v34_sem_test")
    orig_dt_sha = vmod_sem.EXPECTED_DEC113_DECISION_TREE_SHA256
    for sem_label, mutate_blk in (
        (
            "unrecognized_y07_operator",
            lambda b: [c.__setitem__("operator", "BOGUS_OP") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU" and c["gate_id"] == "Y-07"],
        ),
        (
            "placeholder_y07_threshold",
            lambda b: [c.__setitem__("threshold", "TBD") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU" and c["gate_id"] == "Y-07"],
        ),
        (
            "drifted_y07_threshold_0p50",
            lambda b: [c.__setitem__("threshold", c["threshold"].replace("0.85", "0.50")) for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_Y_FULLSTACK_V2_B2M_11TU" and c["gate_id"] == "Y-07"],
        ),
        (
            "duplicate_crosswalk_entry",
            lambda b: b["dec102_endpoint_crosswalk"].append(copy.deepcopy(b["dec102_endpoint_crosswalk"][0])),
        ),
    ):
        dt_mut = json.loads(orig_dt)
        blk_m = dt_mut["design_y_fullstack_v2_and_route3p_gates_iip_0034"]
        mutate_blk(blk_m)
        dt_p.write_text(json.dumps(dt_mut, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        new_sha = hashlib.sha256(
            json.dumps(blk_m, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        try:
            vmod_sem.EXPECTED_DEC113_DECISION_TREE_SHA256 = new_sha
            with pytest.raises(ValueError, match="Unrecognized|placeholder|mismatch|exactly 30|Duplicate"):
                vmod_sem._resolve_and_verify_design_y_v2_registry(min_root)
        finally:
            vmod_sem.EXPECTED_DEC113_DECISION_TREE_SHA256 = orig_dt_sha
    dt_p.write_text(orig_dt, encoding="utf-8")


def test_10_r3p05_exact_decimal_ulp_provenance_alpha_numpy_and_single_dose() -> None:
    vmod = _load_mod("scripts/verify_iip0034_planar_and_portfolio.py", "openislet_v34_r3p05_test")
    valid_scope = {
        "product_variant_id": "ISLET1_Y_FULLSTACK_V2_B2M_11TU",
        "construct_sha256": "e77a358ca5bc55c9e0a3070ab46b115c3314660d3a9be6221ef9e91d698cab4e",
        "donor_seedline_id": "DONOR_SEED_A",
        "mcb_id": "MCB_01",
        "lot_id": "LOT_A",
        "device_assembly_id": "DEV_A",
        "physical_route_id": "ROUTE_3P_MODULAR_PLANAR_ENVELOPE",
    }
    base_ctx = {
        "f_clonogenic": 1.0,
        "f_clonogenic_provenance": "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
        "expected_scope": valid_scope,
        "observed_scope": dict(valid_scope),
        "sampling_independence_model": "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES",
        "clustered_sampling": False,
    }

    # 1. Binary64 ULP straddle MUST HOLD (exact U95 = 4.99999920657284937e-7 > binary64 cap 4.99999920657284866e-7)
    r_ulp = vmod.evaluate_r3p05_dec103_thinning(k=0, n_raw=5_991_464, eta_spike_recovery_l95=1.0, max_u95_fraction=4.999999206572849e-7, **base_ctx)
    assert r_ulp["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"

    # 2. Huge n=10^308 MUST PASS with strictly non-zero scientific bound string
    r_huge = vmod.evaluate_r3p05_dec103_thinning(k=0, n_raw=10**308, eta_spike_recovery_l95=1.0, **base_ctx)
    assert r_huge["status"] == "PASS" and not r_huge["u95_exact_dec103_str"].startswith("0.0")

    # 3. Single-dose Blueprint counterexample (93,350 opportunities -> 3.20909e-5 / 6.41818e-5 -> HOLD)
    r_sd1 = vmod.evaluate_r3p05_dec103_thinning(k=0, n_raw=93_350, eta_spike_recovery_l95=1.0, **base_ctx)
    assert r_sd1["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
    assert abs(float(r_sd1["u95_exact_dec103_str"]) - 3.20909e-5) < 1e-8
    r_sd2 = vmod.evaluate_r3p05_dec103_thinning(k=0, n_raw=93_350, eta_spike_recovery_l95=0.50, **base_ctx)
    assert r_sd2["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
    assert abs(float(r_sd2["u95_exact_dec103_str"]) - 6.41818e-5) < 1e-8

    # 4. NumPy 1-elem and 2-elem arrays for sampling_independence_model -> structured HOLD (no PASS or ValueError)
    for arr in (
        np.array(["INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"]),
        np.array(["INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"]),
    ):
        r_np = vmod.evaluate_r3p05_dec103_thinning(
            k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0,
            f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
            expected_scope=valid_scope, observed_scope=dict(valid_scope),
            sampling_independence_model=arr, clustered_sampling=False,
        )
        assert r_np["status"] == "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE"

    # 5. Mandatory f_clonogenic_provenance and 3-way Bonferroni alpha allocation
    r_noprov = vmod.evaluate_r3p05_dec103_thinning(
        k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0,
        expected_scope=valid_scope, observed_scope=dict(valid_scope),
        sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False,
    )
    assert r_noprov["status"] == "HOLD_MISSING_OR_INVALID_CLONOGENIC_PROVENANCE"
    r_est_no_a = vmod.evaluate_r3p05_dec103_thinning(
        k=0, n_raw=8_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0,
        f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
        expected_scope=valid_scope, observed_scope=dict(valid_scope),
        sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False,
    )
    assert r_est_no_a["status"] == "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION"
    r_est_ok = vmod.evaluate_r3p05_dec103_thinning(
        k=0, n_raw=7_377_760, eta_spike_recovery_l95=0.50, f_clonogenic=1.0,
        f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
        alpha_obs=Decimal("0.025"), alpha_f_clonogenic=Decimal("0.025"), alpha_eta_recovery=Decimal("0.05"),
        expected_scope=valid_scope, observed_scope=dict(valid_scope),
        sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False,
    )
    assert r_est_ok["status"] == "PASS" and r_est_ok["empirical_clearance_status"] == "MATHEMATICAL_COUNT_SUBGATE_ONLY_EMPIRICAL_RETRIEVAL_AND_WHOLE_PRODUCT_PENDING"

    # 5b. Under-allocated L95 recovery alpha (alpha_eta=0.01 < 0.05 while eta=0.50 < 1.0) MUST HOLD!
    r_under_eta = vmod.evaluate_r3p05_dec103_thinning(
        k=0, n_raw=10_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0,
        f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
        alpha_obs=Decimal("0.01"), alpha_eta_recovery=Decimal("0.01"), alpha_f_clonogenic=Decimal("0.08"),
        expected_scope=valid_scope, observed_scope=dict(valid_scope),
        sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False,
    )
    assert r_under_eta["status"] == "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION"
    assert r_under_eta["reason"] == "ETA_SPIKE_RECOVERY_L95_REQUIRES_ALPHA_ETA_AT_LEAST_0P05"

    # 5c. High-precision 120-digit straddle at fixed 1e-6 cap MUST HOLD (immune to half-even Decimal ln/exp rounding!)
    from decimal import localcontext as _localctx
    with _localctx() as _ctx120:
        _ctx120.prec = 120
        _p_true = Decimal(1) - (Decimal("0.05").ln() / Decimal(3_000_000)).exp()
        # Choose f_straddle strictly below _p_true / 1e-6 by 1e-60 so true U95 = 1e-6 + ~1e-66 > 1e-6
        _f_straddle = (_p_true / Decimal("1e-6")) - Decimal("1e-60")
    r_hp_straddle = vmod.evaluate_r3p05_dec103_thinning(
        k=0, n_raw=3_000_000, eta_spike_recovery_l95=Decimal(1), f_clonogenic=_f_straddle,
        f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
        max_u95_fraction=Decimal("1e-6"),
        expected_scope=valid_scope, observed_scope=dict(valid_scope),
        sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False,
    )
    assert r_hp_straddle["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"

    # 6. Control chars, NUL, zero-width chars, and whitespace padding rejected in scope identifiers
    for bad_lot in (" LOT_A ", "LOT\nA", "LOT\x00A", "LOT\u200bA"):
        s_bad = dict(valid_scope, lot_id=bad_lot)
        r_bad = vmod.evaluate_r3p05_dec103_thinning(
            k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0,
            f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
            expected_scope=s_bad, observed_scope=dict(s_bad),
            sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False,
        )
        assert r_bad["status"] == "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING"

    # 7. Unconditional k >= 1 hard veto tracks attribution separately
    r_k1_good = vmod.evaluate_r3p05_dec103_thinning(k=1, n_raw=5_000_000, eta_spike_recovery_l95=1.0, **base_ctx)
    assert r_k1_good["status"] == "FAIL_HARD_STOP" and r_k1_good["attribution_status"] == "EXACT_7KEY_SCOPE_ATTRIBUTED"
    r_k1_unattr = vmod.evaluate_r3p05_dec103_thinning(k=1, n_raw=0, eta_spike_recovery_l95=1.0, max_u95_fraction=10**400)
    assert r_k1_unattr["status"] == "FAIL_HARD_STOP" and r_k1_unattr["attribution_status"].startswith("UNATTRIBUTED_OR_INVALID_AUX_INPUTS_UNCONDITIONAL_STOP")


def test_11_ceg_strict_ownership_source_completeness_and_rolling40_sourceguard(tmp_path: Path) -> None:
    ceg = json.loads((ROOT / "data" / "extracted" / "context_evidence_graph_v1.json").read_text(encoding="utf-8"))
    bceg = _load_mod("scripts/build_context_evidence_graph.py", "openislet_bceg_strict_test")
    assert len(bceg.CURATED_INELIGIBLE_EDGE_IDS) == 66

    det = bceg.verify_ceg_payload_detailed(ceg, repo_root=ROOT)
    has_dossiers = (ROOT / "data" / "dossiers").is_dir()
    assert det["reused_m1522_outcomes_exact_span_verified"] is False
    assert det["izadifar_dossier_contains_figure_legends"] is False
    if has_dossiers:
        assert det["verification_mode"] in ("iip0034_5dossiers_8new_exact_spans_verified", "iip0035_11dossiers_8new_exact_spans_verified")
        assert len(det["iip0034_dossiers_verified"]) == 5
    else:
        assert det["verification_mode"] == "structural_public_only"
        # Requesting require_sources=True without dossiers MUST fail (never silently downgrade)
        with pytest.raises(ValueError, match="Incomplete local dossier set"):
            bceg.verify_ceg_payload(ceg, repo_root=ROOT, require_sources=True)

    # 1. Cross-study comparator_arm_id -> REJECTED
    bad1 = copy.deepcopy(ceg)
    bad1["outcomes"]["PMID:42020736:O1"]["comparator_arm_id"] = "VAR:PMID:17375075:A5_jurkat_htmpk_no_azt_ctrl"
    with pytest.raises(ValueError, match="belongs to"):
        bceg.verify_ceg_payload(bad1, repo_root=ROOT)

    # 2. Arm owner reassigned to another study -> REJECTED
    bad2 = copy.deepcopy(ceg)
    bad2["arms"]["VAR:PMID:42020736:A2_passige_ctrl"]["queue_id"] = "PMID:17375075"
    with pytest.raises(ValueError, match="queue_id|arm_ids"):
        bceg.verify_ceg_payload(bad2, repo_root=ROOT)

    # 3. False outcomes_count -> REJECTED
    bad3 = copy.deepcopy(ceg)
    bad3["studies"]["PMID:42020736"]["outcomes_count"] = 99
    with pytest.raises(ValueError, match="outcomes_count"):
        bceg.verify_ceg_payload(bad3, repo_root=ROOT)

    # 4. Bare A999 unresolved arm token in comparison_arms -> REJECTED
    bad4 = copy.deepcopy(ceg)
    bad4["outcomes"]["PMID:42020736:O1"]["comparison_arms"] = "A1_quadpe vs A999"
    with pytest.raises(ValueError, match="unresolved.*arm"):
        bceg.verify_ceg_payload(bad4, repo_root=ROOT)

    # 5. Wrong embedded outcome_id -> REJECTED
    bad5 = copy.deepcopy(ceg)
    bad5["outcomes"]["PMID:42020736:O1"]["outcome_id"] = "PMID:42020736:O99"
    with pytest.raises(ValueError, match="outcome_id"):
        bceg.verify_ceg_payload(bad5, repo_root=ROOT)

    # 6. Unresolved candidate_effector_topics -> REJECTED
    bad6 = copy.deepcopy(ceg)
    bad6["studies"]["PMID:42020736"]["candidate_effector_topics"] = ["NOT_A_REAL_EFFECTOR"]
    with pytest.raises(ValueError, match="unresolved effector"):
        bceg.verify_ceg_payload(bad6, repo_root=ROOT)

    # 7. Consistently rewired cross-study outcome cited by new IIP-0034 edge -> REJECTED
    bad7 = copy.deepcopy(ceg)
    for e in bad7["edges"]:
        if e["edge_id"] == "PMID:42020736:L1":
            e["cited_outcome_tuples"] = ["PMID:17375075:O1"]
            e["cited_outcome_ids"] = ["PMID:17375075:O1"]
    bad7["outcomes"]["PMID:42020736:O1"]["cited_by_edge_ids"] = []
    bad7["outcomes"]["PMID:17375075:O1"]["cited_by_edge_ids"] = sorted(
        set(bad7["outcomes"]["PMID:17375075:O1"].get("cited_by_edge_ids", [])) | {"PMID:42020736:L1"}
    )
    with pytest.raises(ValueError, match="cross-study outcome"):
        bceg.verify_ceg_payload(bad7, repo_root=ROOT)

    # 8. Typed-arm coherence & nested source_pin.quote_sha256 & renamed M1522 unvalidated_transfer edge -> ALL REJECTED!
    bad8a = copy.deepcopy(ceg)
    bad8a["outcomes"]["PMID:17375075:O1"]["comparator_arm_ids"] = ["VAR:PMID:17375075:A6_primary_tcell_htmpk_no_azt_ctrl"]
    with pytest.raises(ValueError, match="disagree with comparison_arms"):
        bceg.verify_ceg_payload(bad8a, repo_root=ROOT)

    bad8b = copy.deepcopy(ceg)
    bad8b["outcomes"]["PMID:17375075:O1"]["comparator_arm_ids"] = ["VAR:PMID:17375075:A5_jurkat_htmpk_no_azt_ctrl"] * 2
    with pytest.raises(ValueError, match="duplicate entries"):
        bceg.verify_ceg_payload(bad8b, repo_root=ROOT)

    bad8c = copy.deepcopy(ceg)
    bad8c["outcomes"]["PMID:17375075:O1"]["comparator_arm_ids"] = list(bad8c["outcomes"]["PMID:17375075:O1"]["intervention_arm_ids"])
    with pytest.raises(ValueError, match="overlap"):
        bceg.verify_ceg_payload(bad8c, repo_root=ROOT)

    bad8d = copy.deepcopy(ceg)
    bad8d["outcomes"]["PMID:17375075:O1"]["source_pin"]["quote_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="conflicts with top-level quote_sha256"):
        bceg.verify_ceg_payload(bad8d, repo_root=ROOT)

    bad8e = copy.deepcopy(ceg)
    for e in bad8e["edges"]:
        if e["edge_id"] == "M1522:L4":
            e["edge_id"] = "M1522:L99"
            e["primary_certification_eligible"] = True
    for o in bad8e["outcomes"].values():
        if "M1522:L4" in (o.get("cited_by_edge_ids") or []):
            o["cited_by_edge_ids"] = ["M1522:L99" if x == "M1522:L4" else x for x in o["cited_by_edge_ids"]]
    with pytest.raises(ValueError, match="unvalidated_transfer.*cannot be primary_certification_eligible"):
        bceg.verify_ceg_payload(bad8e, repo_root=ROOT)

    # Source-backed tests (when local dossiers are present): partial dossier refusal + rolling-40 + arm edit + reused M1522
    if has_dossiers:
        partial_root = tmp_path / "partial_dossiers"
        (partial_root / "data" / "dossiers").mkdir(parents=True)
        for pm in ("17375075", "39696686", "42020736", "42056515"):  # omit 42749803
            shutil.copy2(ROOT / "data" / "dossiers" / f"{pm}.json", partial_root / "data" / "dossiers" / f"{pm}.json")
        with pytest.raises(ValueError, match="Incomplete local dossier set"):
            bceg.verify_ceg_payload(ceg, repo_root=partial_root)

        # Rolling 40-char substring embedded in unrelated prose -> REJECTED
        shi_abs = json.loads((ROOT / "data" / "dossiers" / "42020736.json").read_text(encoding="utf-8"))["abstract"]
        bad_roll = copy.deepcopy(ceg)
        bad_roll["outcomes"]["PMID:42020736:O1"]["quantitative_summary"] = "Wrapper prefix " + shi_abs[900:955] + " wrapper suffix."
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_roll, repo_root=ROOT)

        # 6-word shingle inside arm edits_or_interventions -> REJECTED
        bad_arm_sh = copy.deepcopy(ceg)
        bad_arm_sh["arms"]["VAR:PMID:42020736:A1_quadpe"]["edits_or_interventions"] = "Test " + shi_abs[860:920]
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_arm_sh, repo_root=ROOT)

        # Reused M1522 outcome tampered with Izadifar shingle -> REJECTED
        iza_abs = json.loads((ROOT / "data" / "dossiers" / "39696686.json").read_text(encoding="utf-8"))["abstract"]
        bad_m1522 = copy.deepcopy(ceg)
        bad_m1522["outcomes"]["M1522:O6"]["quantitative_summary"] = "Injected " + iza_abs[100:170]
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_m1522, repo_root=ROOT)


def test_12_fvm_guards_converged_fourier_451_colonies_and_mechanics() -> None:
    vmod = _load_mod("scripts/verify_iip0034_planar_and_portfolio.py", "openislet_v34_phys_test")
    art = json.loads((ROOT / "data" / "extracted" / "iip0034_planar_and_portfolio_verification_v0.8.0.json").read_text(encoding="utf-8"))

    # 1. FVM input validation guards (negative sleeve, oversized cluster, degenerate grid)
    with pytest.raises(ValueError, match="sleeve_um"):
        vmod.solve_planar_route3p_fvm(label="NEG_SLV", area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=65.0, bed_po2_mmhg=35.0, sleeve_um=-1.0, sleeve_d_factor=0.5, k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097)
    with pytest.raises(ValueError, match="cluster_diam_um"):
        vmod.solve_planar_route3p_fvm(label="OVERSIZE", area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=200.0, bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0, k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097)
    with pytest.raises(ValueError, match="nx|nz|grid"):
        vmod.solve_planar_route3p_fvm(label="TINY_GRID", area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=65.0, bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0, k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097, nx=1, nz=1)

    fvm_by_lbl = {r["label"]: r for r in art["planar_3d_fvm_oxygen_envelopes"]}
    assert fvm_by_lbl["SubARS300_H200_Day0_Clean_65um_PASS"]["passes_15_and_5"] is True
    assert fvm_by_lbl["Hybrid333_H200_Day0_Clean_65um_PASS"]["passes_15_and_5"] is True
    assert float(fvm_by_lbl["Hybrid333_H180_VPIO120_PostAnast_GelOnlyNoIsletPerf_CORE_FAIL"]["worst_core_mmhg"]) < 0.1
    assert fvm_by_lbl["Hybrid333_H180_Day0_10umFibrin_0p5D_65um_FACE_FAIL"]["passes_15_and_5"] is False
    for row in art["planar_3d_fvm_oxygen_envelopes"]:
        assert row["grid_nx_ny_nz"] == [28, 28, 32]
        assert float(row["rel_mass_balance_error_str"]) < 1e-10 and row["rel_mass_balance_error_str"] != "0.0"

    # 2. Converged 5,000-mode vs 45-mode Fourier shutoff & 0.235 s spatial H200 margin
    sh = art["planar_1d_insulin_shutoff_fall_t50"]
    assump = sh["model_assumptions"]
    assert assump["reference_fourier_modes"] == 5000 and assump["comparison_truncated_fourier_modes"] == 45
    assert abs(float(assump["h200_uniform_clean_cutoff_intrinsic_t50_45mode_min"]) - 1.817482) < 1e-4
    assert abs(float(assump["h200_uniform_clean_cutoff_intrinsic_t50_5000mode_min"]) - 1.809045) < 1e-4
    assert abs(float(assump["h200_spatial_65um_clean_cutoff_intrinsic_t50_5000mode_min"]) - 1.75364) < 1e-4
    assert abs(float(assump["h200_spatial_65um_clean_margin_below_3min_seconds"]) - 0.235) < 0.01
    sp_by_case = {r["case"]: r for r in sh["spatially_weighted_1d_projection_cases"]}
    assert abs(float(sp_by_case["Spatial_65um_Bilayer_H180_clean"]["modeled_fall_t50_converged_5000mode_min"]) - 2.846828) < 1e-4
    assert sp_by_case["Spatial_120um_CenteredVPIO_H180_clean_FAIL"]["passes_lt_3p0_min"] is False

    # 3. Exhaustive integer colony minimization (451 across 11 splits)
    col = art["sequential_bxb1_colony_screening_audit"]
    assert col["two_stage_equal_margin_95pct_joint_total"] == 455
    assert col["two_stage_exact_minimum_95pct_joint_total"] == 451
    assert col["two_stage_exact_minimum_95pct_joint_allocation_count"] == 11
    assert col["two_stage_exact_minimum_95pct_joint_allocations"][0] == [261, 190]
    assert col["two_stage_exact_minimum_95pct_joint_allocations"][-1] == [271, 180]

    # 4. Mechanics critical adhesion & concurrent two-ribbon warning
    mech_by_cfg = {m["configuration"]: m for m in art["planar_mechanics_and_clearance"]}
    assert abs(float(mech_by_cfg["ROUTE_3P_SUBARS_300_H180"]["critical_uniform_adhesion_at_15N_kpa"]) - 3.0) < 1e-3
    assert abs(float(mech_by_cfg["ROUTE_3P_HYBRID_333_H180"]["critical_uniform_adhesion_at_15N_kpa"]) - 2.7003) < 1e-3
    assert abs(float(mech_by_cfg["ROUTE_3P_4MOD_SUBARS_260_H180"]["critical_uniform_adhesion_at_15N_kpa"]) - 2.3077) < 1e-3
    assert abs(float(mech_by_cfg["ROUTE_3P_HYBRID_333_H180"]["two_ribbon_concurrent_2kpa_pullout_warning_n"]) - 22.22) < 0.02


def test_13_source_free_public_reexport_preserves_metadata_and_ranker_profiles(tmp_path: Path) -> None:
    """Source-free public re-export must preserve structured metadata, 262 expanded genes, and 51 target profiles."""
    grm_path = ROOT / "scripts" / "generate_release_manifest.py"
    if not grm_path.is_file():
        pytest.skip("generate_release_manifest.py not present in minimal checkout")

    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    reexp_out = tmp_path / "reexported_pub"
    if (ROOT / "data" / "dossiers").is_dir():
        # Full-repo checkout: export a canonical public bundle first, then re-export source-free from it!
        pub_stage = tmp_path / "pub_stage_1"
        subprocess.run(
            [sys.executable, "-B", str(grm_path), "--export-public-bundle", str(pub_stage)],
            cwd=str(ROOT),
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        src_free = pub_stage
    else:
        src_free = ROOT
    assert not (src_free / "data" / "dossiers").exists()

    subprocess.run(
        [sys.executable, "-B", str(src_free / "scripts" / "generate_release_manifest.py"), "--export-public-bundle", str(reexp_out)],
        cwd=str(src_free),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )

    # Verify manifest snapshot integrity on re-exported bundle
    grm_re = _load_mod("scripts/generate_release_manifest.py", "openislet_grm_reexp_test")
    snap = grm_re.verify_public_bundle_snapshot(reexp_out)
    assert snap["passed"] is True

    # Verify Izadifar cell_line_id, Liu MULTI_LOCUS_SAFE_HARBOR & HEK293T, and all 5 new measurement_units survive source-free re-export!
    iza_re = json.loads((reexp_out / "data" / "extracted" / "39696686.json").read_text(encoding="utf-8"))
    assert iza_re["arms"][0]["gate1_substrate_and_differentiation"]["cell_line_id"].startswith("beta-iPLC line 2")
    assert iza_re["atomic_effects"][1]["measurement_unit"] == "ng_per_mL_per_1e6_cells"
    assert iza_re["arms"][0]["gate2_3_4_genetic_and_epitope_edits"][0]["edit_modality"] == "Endogenous_Locus_KI_Without_KO"

    liu_re = json.loads((reexp_out / "data" / "extracted" / "42056515.json").read_text(encoding="utf-8"))
    assert liu_re["arms"][0]["gate1_substrate_and_differentiation"]["cell_line_id"] == "HEK293T"
    for arm in liu_re["arms"]:
        for ed in arm.get("gate2_3_4_genetic_and_epitope_edits", []):
            assert ed["gene_symbol_hgnc"] == "MULTI_LOCUS_SAFE_HARBOR"

    # Verify expanded gene count == 262 and target evidence profiles == 51 (zero spurious redacted-token profiles!)
    rank_mod = _load_mod("scripts/rank_islet1_stacks.py", "openislet_rank_reexp_test")
    re_studies = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((reexp_out / "data" / "extracted").glob("[0-9]*.json"))]
    assert len(re_studies) == 604
    expanded_genes = [
        tok
        for st in re_studies
        for a in st.get("arms", [])
        for ed in a.get("gate2_3_4_genetic_and_epitope_edits", [])
        for tok in rank_mod.expand_gene_symbol_tokens(str(ed.get("gene_symbol_hgnc") or ""))
    ]
    assert len(expanded_genes) == 262, f"Expected 262 expanded genes after source-free re-export, got {len(expanded_genes)}"
    _, _, omics = rank_mod.load_inputs()
    profs = rank_mod.build_target_evidence_profiles(re_studies, omics)
    assert len(profs) == 51, f"Expected 51 target profiles after source-free re-export, got {len(profs)}"
    for k_prof in profs:
        assert not k_prof.startswith("[REDACTED") and k_prof not in {"MULTI", "LOCUS", "SAFE", "HARBOR"}

