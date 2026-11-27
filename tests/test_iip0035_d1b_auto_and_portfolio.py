#!/usr/bin/env python3
"""Focused standalone Option-2 (Zero-LIMS-Touch) verification suite for IIP-0035 (DEC-114).

Covers:
  1. Autologous D1b-Auto-C v2 Mode-B FASTA SHA-256, 5 role-scoped donor records (L1A/L1B full-product
     candidates vs L2/L3/L4 component donors strictly EXCLUDED_FROM_FULLSTACK_CLEARANCE), whole-record
     canonical full_record_sha256 pinning, nonfinite JSON/1e400 rejection, and product count 33 -> 38.
  2. public_excerpt_policy 11-PMID license registry (CC BY 4.0 / CC BY-NC 4.0 / CC BY-NC-ND 4.0 vs
     non-OA dual redaction prefix separation).
  3. New manifest-pinned data/extracted/islet1_5arm_autologous_portfolio_v0.8.0.json replay: 11 Wave-C
     4-arm rows, 12 Wave-B+C 5-arm rows, 2 W_U9 phi_C counterexamples, 7-point p1 sweep, continuous
     affine p1* (0.672989 / 0.746615), 1,048,576-state (2^20) exact unrounded alias verification
     (W_BC1..W_BC5 <-> W_U1, W_U2, W_U8, W_U11, W_U12, max_abs_diff == 0.0 <= 1e-12), NaN mass guard,
     and 0-byte invariance of frozen islet1_multiparadigm_portfolio_v0.8.0.json & Pareto JSON.
  4. Standalone cross-runtime verifier scripts/verify_iip0035_d1b_auto_and_portfolio.py (--check parity,
     ZERO wetlab_contract.py import verified via AST + sys.modules, ZERO assert statements / python -O
     integrity enforcement, mutually exclusive --write --check CLI, CRLF/padding tamper rejection,
     DEC-114 decision-tree SHA-256, CommonMark 0.31.2 Playbook s4.1.20 attacks & contradiction guards,
     and re-pinned semantic cell tamper rejection).
  5. Exact 80-digit sound Decimal Clopper-Pearson rare-event evaluator (.next_minus() on ln & exp,
     alpha_eta >= 0.05 when eta < 1.0, 120-digit straddle immune, huge-n=10^308 non-zero scientific
     string, D1-04 299,573,225 vs 299,573,226 & D1-08 58 [5.0339%] vs 59 [4.9508%], explicit k >= 1
     FAIL_HARD_STOP) + standalone evaluate_stage0_d1b_auto_v2_gates_iip0035 reference matrix (structured
     lot-matched global_gate0v_record with mcb_lot_id, bare-boolean rejection, L2/L3/L4 exclusion,
     SP-3C/SP-4C C*18:01 compound-het exclusion, k>=1 consistency, and "UNAVAILABLE" hold).
  6. Ring-1 corpus expansion (594 -> 604 studies, 707 arms, 65 ref-control arms, 261 raw edits / 262
     expanded gene occurrences, 3,174 effects, 88 curated landmarks; Jetley 41201433 retained in CEG
     M0051 + dossier cache only so its triple HLA-A/B/CIITA edits do not mutate frozen Pareto),
     structured challenge/edit annotations (Amdare 41764741 NonIslet_Immune_Endothelial_or_Hematopoietic_Model /
     Peptide_Restimulation_Islet_Derived_T_Cell_Line_Assay, Uenaka 41690310 MULTI_LOCUS_SAFE_HARBOR,
     Sutherland/Deuse/Lin InVivo_Allogeneic_or_Autoimmune_Graft, Deuse 31427818 2-arm mouse in-vivo +
     human in-vitro Autologous_T_Cell_Neoantigen_InVitro_Assay), Lin 38843830 (Tier_2C_NHP_Allogeneic_or_Autologous_NonIslet
     / NHP_iPSC_Cardiomyocyte_NonIslet), Rouhani 35953586 (Fig. 5 neural differentiation penalty), Cheng 41999746
     (2 arms), and 0-byte preservation of 17375075 & 39571901.
  7. Context-Evidence Graph expansion (1,176 -> 1,182 studies, 7,749 arms, 14,175 outcomes, 2,494 edges,
     65 curated ineligible edges, typed intervention_arm_ids/comparator_arm_ids on all 8 new outcomes,
     explicit provenance_note on all 5 additive cross-links, verify_ceg_payload_detailed reporting
     iip0035_11dossiers_8new_exact_spans_verified & reused_legacy_crosslink_outcomes_exact_span_verified=False,
     partial-dossier refusal, and rolling-40 _SourceGuard).
  8. Extraction schema additions, 0 / 2,160 Layer-1 stack invariance, 0-byte invariance of all 3 protected
     LIMS files (366 / 77 / 122), and source-free public re-export compatibility.
"""
from __future__ import annotations

import ast
import copy
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

EXPECTED_DEC114_DT_SHA256 = "1699fa9bd04dedebd4e4846c75a4b72db10d18ec9f06d80e6004fbf37177299f"
EXPECTED_DEC114_PB_SHA256 = "e3fb5ec701678dfdeb62a77c15866327842b60b6b4c503cb86dfc79885a626a8"
EXPECTED_PARETO_SHA256 = "591ddeb9a706c0bc95fc509765b512aadad8a403287a796e63dd52e7db548575"
EXPECTED_IIP0034_PORTFOLIO_SHA256 = "3ed4d53f1756b8b8615fbbc4e9a15755b9ab4968d206da61160eeb04d833cd45"
EXPECTED_IIP0034_VERIF_SHA256 = "d95970fb9423ff5307e9f481ed7027688fdff7bd558e6c439f91f63c734bd08e"
EXPECTED_LIMS_SHA256 = {
    "scripts/wetlab_contract.py": "331ecd2f2d8602dfceae8b89263d0f360cfadfd0638a46d7d668a376294cfc73",
    "scripts/ingest_wetlab_results.py": "e050593675170515382257ce379449db596c87cc06ffc26d31d128e127630cb1",
    "data/protocols/wetlab_endpoint_contract_v0.8.0.json": "757bc1d9c2722c129cc4fd5962e8995b9f5bc5f7035f8530f5716da9aa2ac26f",
}


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
        "scripts/verify_iip0035_d1b_auto_and_portfolio.py",
        "scripts/product_contract.py",
        "scripts/engineering_priors.py",
        "scripts/public_excerpt_policy.py",
        "scripts/rank_islet1_stacks.py",
        "scripts/wetlab_contract.py",
        "scripts/ingest_wetlab_results.py",
    ):
        shutil.copy2(ROOT / rel_s, dst / rel_s)
    for rel_d in (
        "data/extracted/iip0035_d1b_auto_and_portfolio_verification_v0.8.0.json",
        "data/extracted/islet1_5arm_autologous_portfolio_v0.8.0.json",
        "data/extracted/islet1_multiparadigm_portfolio_v0.8.0.json",
        "data/extracted/iip0034_planar_and_portfolio_verification_v0.8.0.json",
        "data/extracted/islet1_combinatorial_pareto_and_ladder_v0.8.0.json",
        "data/protocols/islet1_product_definitions_v0.8.0.json",
        "data/protocols/stage0_decision_tree_v1.json",
        "data/protocols/wetlab_endpoint_contract_v0.8.0.json",
        "docs/STAGE-0-EXECUTION-PLAYBOOK.md",
    ):
        shutil.copy2(ROOT / rel_d, dst / rel_d)
    for fa in (ROOT / "data" / "constructs").glob("*.fasta"):
        shutil.copy2(fa, dst / "data" / "constructs" / fa.name)
    return dst


def test_01_d1b_auto_v2_fasta_sha256_and_strict_product_contract(tmp_path: Path) -> None:
    pc = _load_mod("scripts/product_contract.py", "openislet_pc_iip0035_test")
    fasta_path = ROOT / pc.D1B_AUTO_V2_FASTA_PATH
    assert fasta_path.is_file()
    raw = fasta_path.read_bytes()
    assert len(raw) == 21516
    assert hashlib.sha256(raw).hexdigest() == pc.D1B_AUTO_V2_FASTA_SHA256

    defs = pc.load_product_contracts(ROOT, require_staged_fastas=True)
    assert len(defs["products"]) == 38
    by_id = {p["product_variant_id"]: p for p in defs["products"]}
    for pvid in pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS:
        assert pvid in by_id
        p = by_id[pvid]
        assert p["execution_supported"] is False
        assert p["lifecycle_state"] == "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"
        assert p["physical_route_id"] == "ROUTE_3P_MODULAR_PLANAR_ENVELOPE"

    for comp_id in (
        "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1",
        "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9",
        "ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y",
    ):
        assert "EXCLUDED_FROM_FULLSTACK_CLEARANCE" in by_id[comp_id]["staged_d1b_auto_v2_donor_binding"]["stage0_promotion_gates"]
        assert "Not encoded by this standalone" in by_id[comp_id]["hla_retention"]["classical_MHC_I"]

    copy_root = tmp_path / "repo"
    (copy_root / "data" / "protocols").mkdir(parents=True)
    (copy_root / "data" / "constructs").mkdir(parents=True)
    orig_prod_path = ROOT / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json"
    p_file = copy_root / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json"
    shutil.copy2(orig_prod_path, p_file)
    for fa in (ROOT / "data" / "constructs").glob("*.fasta"):
        shutil.copy2(fa, copy_root / "data" / "constructs" / fa.name)
    assert len(pc.load_product_contracts(copy_root, require_staged_fastas=True)["products"]) == 38

    # 1. Tamper FASTA -> rejected
    t_fa = copy_root / pc.D1B_AUTO_V2_FASTA_PATH
    t_fa.write_bytes(t_fa.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)
    shutil.copy2(fasta_path, t_fa)

    # 2. Flip execution_supported=True -> rejected
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN":
            p["execution_supported"] = True
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 3. Mutate physical_route_id -> rejected
    p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
    for p in p_doc["products"]:
        if p["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN":
            p["physical_route_id"] = "ROUTE_2_OPEN_PORE_PREVASC"
    p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="ROUTE_3P_MODULAR_PLANAR_ENVELOPE"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 4. Mutate pending_assay_gates / donor_sha256 -> caught by canonical full_record_sha256
    for mut_fn in (
        lambda p: p["staged_d1b_auto_v2_donor_binding"].__setitem__("donor_sha256", "0" * 64),
        lambda p: p.__setitem__("pending_assay_gates", ["D1-09 alone"]),
        lambda p: p.__setitem__("execution_scope", p["execution_scope"] + " TAMPERED"),
    ):
        p_doc = json.loads(orig_prod_path.read_text(encoding="utf-8"))
        for p in p_doc["products"]:
            if p["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR":
                mut_fn(p)
        p_file.write_text(json.dumps(p_doc, indent=2) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match="full_record_sha256"):
            pc.load_product_contracts(copy_root, require_staged_fastas=True)

    # 5. Nonfinite JSON / 1e400 float overflow rejected
    raw_ovf = orig_prod_path.read_text(encoding="utf-8").replace("{", '{"overflow_test": 1e400, ', 1)
    p_file.write_text(raw_ovf, encoding="utf-8")
    with pytest.raises(ValueError, match="Non-finite JSON float overflow"):
        pc.load_product_contracts(copy_root, require_staged_fastas=True)


def test_02_public_excerpt_policy_11_pmids_and_dual_prefix() -> None:
    pep = _load_mod("scripts/public_excerpt_policy.py", "openislet_pep_iip0035_test")
    expected_lic = {
        "2642862": "PUBLISHER_STANDARD_NON_OA",
        "31427818": "PUBLISHER_STANDARD_NON_OA",
        "35953586": "CC_BY_4_0",
        "36944335": "PUBLISHER_STANDARD_NON_OA",
        "38843830": "PMC_AUTHOR_MANUSCRIPT_NON_OA",
        "39753706": "PUBLISHER_STANDARD_NON_OA",
        "41201433": "CC_BY_NC_ND_4_0",
        "41436438": "CC_BY_4_0",
        "41690310": "CC_BY_NC_4_0",
        "41764741": "PUBLISHER_STANDARD_PENDING_PMC_EMBARGO",
        "41999746": "PUBLISHER_STANDARD_NON_OA",
    }
    for pm, lic in expected_lic.items():
        assert pep.PUBLISHER_LICENSE_REGISTRY.get(pm) == lic
    assert pep.redaction_token_prefix_for_pmid("41436438") == "[REDACTED_CONSERVATIVE_EXCERPT:"
    assert pep.redaction_token_prefix_for_pmid("35953586") == "[REDACTED_CONSERVATIVE_EXCERPT:"
    assert pep.redaction_token_prefix_for_pmid("41201433") == "[REDACTED_CONSERVATIVE_EXCERPT:"
    assert pep.redaction_token_prefix_for_pmid("38843830") == "[REDACTED_NON_OA_EXCERPT:"
    assert pep.redaction_token_prefix_for_pmid("2642862") == "[REDACTED_NON_OA_EXCERPT:"


def test_03_5arm_autologous_portfolio_replay_alias_2pow20_and_iip0034_invariance() -> None:
    rank_mod = _load_mod("scripts/rank_islet1_stacks.py", "openislet_rank_iip0035_test")
    p5_path = ROOT / "data" / "extracted" / "islet1_5arm_autologous_portfolio_v0.8.0.json"
    assert p5_path.is_file()
    shipped = json.loads(p5_path.read_text(encoding="utf-8"))
    recomputed = rank_mod.compute_5arm_autologous_portfolio_iip0035(repo_root=ROOT)
    assert shipped == recomputed

    assert shipped["proposal_id"] == "IIP-0035" and shipped["decision_id"] == "DEC-114"
    assert shipped["live_headline_delta_pp_today"] == 0.0
    assert shipped["five_arm_bridge_registration_status"]["D1B_AUTO_C_V2"] == "REGISTERED_AND_FASTA_VERIFIED_IN_IIP0035_MODE_B_NON_EXECUTABLE_ROUTE_3P_ONLY"
    gap = shipped["global_access_positioning"]
    assert "unified_5arm_verdict_zero_h3_lift_conditional" not in gap
    assert "FAIL registered-route locked" in gap["unified_5arm_verdict_registered_route_only_locked"]
    assert "59.15-60.64% at p1=0.62" in gap["unified_5arm_verdict_registered_route_only_locked"]
    assert "W_U8 [phi_YX=0.75, phi_C=0.50]" in gap["unified_5arm_verdict_class_level_route_transfer"]
    assert "W_U9 [phi_YX=1.00, phi_C=0.50]" in gap["unified_5arm_verdict_class_level_route_transfer"]

    wave4 = shipped["stage0_gated_4arm_autologous_d1b_auto_v2"]
    wave5 = shipped["stage0_gated_5arm_unified_waveBC_portfolio"]
    cex = shipped["stage0_gated_5arm_phi_c_sensitivity_counterexamples"]
    assert len(wave4) == 11 and len(wave5) == 12 and len(cex) == 2
    assert len(shipped["autologous_p1_breakeven_sweep"]) == 7

    by4 = {r["scenario"]: r for r in wave4}
    by5 = {r["scenario"]: r for r in wave5}
    w_c0 = by4["W_C0_CTRL_Unedited_Auto_RecurrenceStress@R2 (modeled p3=0.10; motivated by Sutherland1989 & Vendrame2010)"]
    assert w_c0["cell_architecture_classes"][3] == "UNEDITED_AUTO_NEGATIVE_CONTROL"
    assert abs(float(w_c0["locked_pct"]) - 51.5102) < 1e-3
    assert abs(float(by4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["locked_pct"]) - 57.3555) < 1e-3
    assert abs(float(by4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_k100_pct"]) - 59.2886) < 1e-3
    assert abs(float(by4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_registered_route_only_k100_pct"]) - 59.2886) < 1e-3
    assert abs(float(by4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_k100_under_class_level_r1_r2_expansion_pct"]) - 60.0661) < 1e-3
    assert float(by4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_k100_under_r2_only_expansion_pct"]) <= float(by4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_k100_under_class_level_r1_r2_expansion_pct"])
    assert abs(float(by4["W_C5_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.50"]["locked_pct"]) - 56.7069) < 1e-3
    assert abs(float(by5["W_U1_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.70)@3P_c4u0.92_phiYX0.75_25"]["locked_pct"]) - 61.0374) < 1e-3
    assert abs(float(by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["modular_k100_pct"]) - 61.5921) < 1e-3
    assert abs(float(by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["modular_registered_route_only_k100_pct"]) - 61.5921) < 1e-3
    assert abs(float(by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["modular_k100_under_class_level_r1_r2_expansion_pct"]) - 62.3076) < 1e-3
    assert abs(float(by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["class_level_route_transfer_delta_pp"]) - 0.7155) < 1e-3
    # Verify exact scope unification (all unqualified metrics, exclusion mass, & pass flags strictly use registered-route-only scope;
    # R2-only and class-level R1+R2 counterfactual metrics carry explicit suffixes!)
    for row_all in wave4 + wave5 + cex:
        assert float(row_all["locked_pct"]) == float(row_all["locked_registered_route_only_pct"])
        assert float(row_all["modular_k100_pct"]) == float(row_all["modular_registered_route_only_k100_pct"])
        assert float(row_all["modular_k90_pct"]) == float(row_all["modular_registered_route_only_k90_pct"])
        assert float(row_all["modular_registered_route_only_k100_pct"]) <= float(row_all["modular_k100_under_r2_only_expansion_pct"]) + 1e-4
        assert float(row_all["modular_k100_under_r2_only_expansion_pct"]) <= float(row_all["modular_k100_under_class_level_r1_r2_expansion_pct"]) + 1e-4
        assert float(row_all["prohibited_locked_mass_excluded_pp"]) >= float(row_all["prohibited_locked_mass_excluded_under_class_level_expansion_pp"]) - 1e-4
        assert float(row_all["prohibited_modular_mass_excluded_pp"]) >= float(row_all["prohibited_modular_mass_excluded_under_class_level_expansion_pp"]) - 1e-4
        assert abs(
            float(row_all["class_level_route_transfer_delta_pp"])
            - (float(row_all["modular_k100_under_class_level_r1_r2_expansion_pct"]) - float(row_all["modular_registered_route_only_k100_pct"]))
        ) <= 0.00015
        assert row_all["clears_57p92_locked_registered_route_only"] is (float(row_all["locked_registered_route_only_pct"]) > 57.917640)
        assert row_all["clears_57p92_locked"] is row_all["clears_57p92_locked_registered_route_only"]
        assert row_all["clears_57p92_locked_under_hypothetical_r2_transfer"] is (float(row_all["locked_under_hypothetical_r2_transfer_pct"]) > 57.917640)
        assert row_all["clears_qL_14p51_and_qM_8p29"] is row_all["clears_qL_14p51_and_qM_8p29_registered_route_only"]
        assert row_all["clears_baseline_implied_residual_risk_barriers"] is row_all["clears_qL_14p51_and_qM_8p29_registered_route_only"]
        assert abs(float(row_all["baseline_implied_qL_pct"]) - 14.5046) < 1e-4
        assert abs(float(row_all["baseline_implied_qM_pct"]) - 8.2848) < 1e-4
    assert shipped["registered_staged_mode_b_routes_from_product_contract"] == {
        "D1B_AUTO_C_V2": ["R3P"],
        "DESIGN_Y_FULLSTACK_V2": ["R3P"],
    }
    assert shipped["incumbent_and_control_route_rules_inherited_from_iip0034"] == {
        "UNEDITED_AUTO_NEGATIVE_CONTROL": ["R1", "R2"],
        "UNHARDENED_A_B_X": ["R1", "R2"],
    }
    assert "DESIGN_Y_HYPOTHETICAL_R2_COUNTERFACTUAL" not in shipped["incumbent_and_control_route_rules_inherited_from_iip0034"]
    assert "DESIGN_Y_HYPOTHETICAL_R2_COUNTERFACTUAL" in shipped["modular_compatibility_matrix_class_level_hypothetical"]
    for row_hyp in wave4 + wave5 + cex:
        if row_hyp["has_unregistered_h4_assignment"]:
            assert not str(row_hyp["stage0_qualification_status"]).endswith("_PASS"), (
                f"Hypothetical R2 row carries misleading _PASS: {row_hyp['scenario']} -> {row_hyp['stage0_qualification_status']}"
            )
            assert row_hyp["clears_57p92_locked_registered_route_only"] is False
            assert row_hyp["clears_qL_14p51_and_qM_8p29_registered_route_only"] is False
            assert row_hyp["clears_qL_14p51_and_qM_8p29"] is False
            assert float(row_hyp["prohibited_locked_mass_excluded_pp"]) > 0.0
    for _hyp_u in (
        "W_U5_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@R2_c4u0.90_phiYX0.75_50",
        "W_U7_5Arm_HOSTILE_ZeroH3Lift_Y@3P(0.68)_D1(0.62,0.65)@R2_c4u90_50",
        "W_U8_5Arm_EXTREME_HOSTILE_Y@3P(0.68)_D1(0.55,0.60)@R2_c4u90_phi50",
        "W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90",
    ):
        assert by5[_hyp_u]["clears_qL_14p51_and_qM_8p29"] is False
        assert by5[_hyp_u]["clears_qL_14p51_and_qM_8p29_registered_route_only"] is False
        assert by5[_hyp_u]["clears_qL_14p51_and_qM_8p29_under_hypothetical_r2_transfer"] is True
        assert float(by5[_hyp_u]["locked_pct"]) < 57.917640
    assert abs(float(by5["W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90"]["locked_pct"]) - 55.9793) < 1e-3
    assert abs(float(by5["W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90"]["locked_under_hypothetical_r2_transfer_pct"]) - 58.1921) < 1e-3
    assert by5["W_U8_5Arm_EXTREME_HOSTILE_Y@3P(0.68)_D1(0.55,0.60)@R2_c4u90_phi50"]["clears_57p92_locked_registered_route_only"] is False
    assert by5["W_U8_5Arm_EXTREME_HOSTILE_Y@3P(0.68)_D1(0.55,0.60)@R2_c4u90_phi50"]["has_unregistered_h4_assignment"] is True
    assert by5["W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90"]["clears_57p92_locked_registered_route_only"] is False
    assert by5["W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90"]["has_unregistered_h4_assignment"] is True
    assert by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["clears_57p92_locked_registered_route_only"] is True
    assert by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["has_unregistered_h4_assignment"] is False
    assert by5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["has_unregistered_modular_route_sensitivity"] is True
    assert abs(float(cex[0]["locked_under_hypothetical_r2_transfer_pct"]) - 57.8614) < 1e-3 and cex[0]["clears_57p92_locked"] is False and cex[0]["clears_57p92_locked_under_hypothetical_r2_transfer"] is False
    assert abs(float(cex[1]["locked_under_hypothetical_r2_transfer_pct"]) - 57.6628) < 1e-3 and cex[1]["clears_57p92_locked"] is False and cex[1]["clears_57p92_locked_under_hypothetical_r2_transfer"] is False

    # Thread-safety check: concurrent calls to compute_5arm_autologous_portfolio_iip0035() return identical dicts!
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as _pool:
        _conc_res = list(_pool.map(lambda _: rank_mod.compute_5arm_autologous_portfolio_iip0035(repo_root=ROOT), range(4)))
    for _cr in _conc_res:
        assert _cr == shipped

    p1e = shipped["autologous_p1_continuous_breakeven_exact"]
    assert abs(float(p1e["phi_0p25_p1_star"]) - 0.672989) < 1e-6
    assert abs(float(p1e["phi_0p50_p1_star"]) - 0.746615) < 1e-6

    amap = shipped["scenario_alias_deduplication_map"]
    assert float(amap["alias_full_joint_state_max_abs_diff"]) == 0.0
    assert len(amap["repeated_bridge_aliases"]) == 5
    assert amap["distinct_core_scenario_laws"] == amap["total_labeled_core_rows_across_iip0034_and_iip0035"] - 5

    # NaN rejection in _max_abs_full_joint_state_diff_iip0035
    with pytest.raises(ValueError, match="Non-finite"):
        rank_mod._max_abs_full_joint_state_diff_iip0035(
            ((float("nan"), 1.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0)),
            ((0.0, 1.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0)),
        )

    # Exact 0-byte invariance of frozen IIP-0034 & Pareto artifacts
    assert hashlib.sha256((ROOT / "data" / "extracted" / "islet1_multiparadigm_portfolio_v0.8.0.json").read_bytes()).hexdigest() == EXPECTED_IIP0034_PORTFOLIO_SHA256
    assert hashlib.sha256((ROOT / "data" / "extracted" / "iip0034_planar_and_portfolio_verification_v0.8.0.json").read_bytes()).hexdigest() == EXPECTED_IIP0034_VERIF_SHA256
    assert hashlib.sha256((ROOT / "data" / "extracted" / "islet1_combinatorial_pareto_and_ladder_v0.8.0.json").read_bytes()).hexdigest() == EXPECTED_PARETO_SHA256


def test_04_standalone_verifier_zero_lims_import_python_O_and_commonmark_attacks(tmp_path: Path) -> None:
    v_script = ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"
    art_path = ROOT / "data" / "extracted" / "iip0035_d1b_auto_and_portfolio_verification_v0.8.0.json"
    assert v_script.is_file() and art_path.is_file()

    # 1. Static AST audit: ZERO assert statements and ZERO wetlab_contract / ingest_wetlab_results imports!
    tree = ast.parse(v_script.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        assert not isinstance(node, ast.Assert), "Verifier must contain ZERO assert statements (python -O safety)"
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert "wetlab_contract" not in alias.name and "ingest_wetlab_results" not in alias.name
        elif isinstance(node, ast.ImportFrom) and node.module:
            assert "wetlab_contract" not in node.module and "ingest_wetlab_results" not in node.module

    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    for extra_flags in (["-B"], ["-B", "-O"]):
        res = subprocess.run(
            [sys.executable, *extra_flags, str(v_script), "--check"],
            cwd=str(ROOT), env=env, capture_output=True, text=True, check=True,
        )
        assert "OK:" in res.stdout

    doc = json.loads(art_path.read_text(encoding="utf-8"))
    assert doc["proposal_id"] == "IIP-0035" and doc["decision_id"] == "DEC-114"
    assert doc["option_2_zero_lims_touch"] is True
    reg = doc["d1b_auto_v2_product_registry_audit"]
    assert reg["dec114_decision_tree_sha256_verified"] == EXPECTED_DEC114_DT_SHA256
    assert reg["dec114_playbook_sec4120_sha256_verified"] == EXPECTED_DEC114_PB_SHA256
    assert reg["typed_product_by_gate_contract_audit"]["cells_verified"] == 45
    assert reg["typed_product_by_gate_contract_audit"]["shared_predicates_verified_zero_skipped"] == 10
    assert reg["typed_product_by_gate_contract_audit"]["l1b_specific_predicates_verified_zero_skipped"] == 1

    # 2. Minimal isolated tree: CLI mutual exclusion, python -O tamper rejection, CommonMark 0.31.2 attacks, LIMS tamper rejection
    min_root = _build_minimal_verifier_tree(tmp_path / "min_v35_tree")
    min_v = min_root / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"
    assert subprocess.run([sys.executable, "-B", str(min_v), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode == 0
    assert subprocess.run([sys.executable, "-B", str(min_v), "--write", "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0

    # LIMS tamper caught by verifier --check (even under python -O!)
    wc_p = min_root / "scripts" / "wetlab_contract.py"
    orig_wc = wc_p.read_bytes()
    wc_p.write_bytes(orig_wc + b"\n")
    assert subprocess.run([sys.executable, "-B", "-O", str(min_v), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0
    wc_p.write_bytes(orig_wc)

    # CommonMark 0.31.2 Playbook §4.1.20 attacks
    pb_p = min_root / "docs" / "STAGE-0-EXECUTION-PLAYBOOK.md"
    orig_pb = pb_p.read_text(encoding="utf-8")
    pb_p.write_text(orig_pb.replace("#### 4.1.20", "```markdown\n#### 4.1.20", 1) + "\n```\n", encoding="utf-8")
    assert subprocess.run([sys.executable, "-B", str(min_v), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0
    for competing in (
        "\n\n4.1.20 Competing multiline\nSetext heading\n---\n",
        "\n\n> #### 4.1.20 Blockquoted duplicate\n",
        "\n\n- #### 4.1.20 List duplicate\n",
        "\n\n#### 4&#46;1&#46;20 HTML decimal char-ref duplicate\n",
        "\n\n#### 4&#x2e;1&#x2e;20 HTML hex char-ref duplicate\n",
        "\n\n#### 4\\.1\\.20 Backslash-escaped dot duplicate\n",
        "\n\n> ```markdown\n> quoted unclosed fence (CommonMark 0.31.2 example 128)\n\n#### 4.1.20 Escaped second heading after quoted fence\n",
        "\n\n1. First item\n\n   #### 4.1.20 Ordered-list continuation heading\n",
        "\n\n10. Tenth item\n\n    #### 4.1.20 Four-space ordered-list continuation heading\n",
        "\n\n[foo]: /url\n4.1.20 Real Setext heading after link ref def\n---\n",
        "\n\n#### 4.1.<span>20</span> Inline HTML tag split heading\n",
    ):
        pb_p.write_text(orig_pb + competing, encoding="utf-8")
        assert subprocess.run([sys.executable, "-B", str(min_v), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0
    for ignored in (
        "\n\n```markdown\n#### 4.1.20 Literal fenced backtick example\n```\n",
        "\n\n~~~\n#### 4.1.20 Literal fenced tilde example\n~~~\n",
        "\n\n\t#### 4.1.20 Tab-indented code example\n",
        "\n\n<!--\n#### 4.1.20 HTML comment example\n-->\n",
        "\n\n<pre>\n#### 4.1.20 HTML pre block example\n</pre>\n",
        "\n\n<script>\n#### 4.1.20 HTML script block example\n</script>\n",
        "\n\n<div>\n#### 4.1.20 HTML div block example\n</div>\n",
        "\n\n[4.1.20 link ref def]: /url\n---\n",
    ):
        pb_p.write_text(orig_pb + "\n\n#### 4.1.21 Non-binding example section\n" + ignored, encoding="utf-8")
        assert subprocess.run([sys.executable, "-B", str(min_v), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode == 0
    # Visible prose, HTML5 short comment <!-->, or tampered OPENISLET_CURRENT_LAYER_SCOPE footer under §4.1.20 MUST fail closed!
    for pb_smuggle in (
        orig_pb.replace(
            "<!-- END_IIP0035_PLAYBOOK_SECTION_4_1_20 -->",
            "<!-- END_IIP0035_PLAYBOOK_SECTION_4_1_20 -->\n\nOverride D1-06 NK lysis cap to <= 99.0%.\n",
        ),
        orig_pb.replace(
            "<!-- END_IIP0035_PLAYBOOK_SECTION_4_1_20 -->",
            "<!-- END_IIP0035_PLAYBOOK_SECTION_4_1_20 -->\n<!--> Smuggled visible text after HTML5 short comment -->\n",
        ),
        orig_pb.replace(
            "Layer 5 is an illustrative",
            "Gate D1-06 is waived. Layer 5 is an illustrative",
        ),
    ):
        pb_p.write_text(pb_smuggle, encoding="utf-8")
        assert subprocess.run([sys.executable, "-B", str(min_v), "--check"], cwd=str(min_root), env=env, capture_output=True).returncode != 0
    pb_p.write_text(orig_pb, encoding="utf-8")
    # Direct unit test of CommonMark 0.31.2 Example 216 ([foo]: /url\n--- -> 0 headings) & thematic breaks (*** / ___ before === / ---)
    vmod_cm = _load_mod("scripts/verify_iip0035_d1b_auto_and_portfolio.py", "openislet_v35_cm216_test")
    assert vmod_cm._parse_commonmark_headings("[foo]: /url\n---\n") == []
    assert vmod_cm._parse_commonmark_headings("***\n===\n") == []
    assert vmod_cm._parse_commonmark_headings("___\n---\n") == []
    h_ex = vmod_cm._parse_commonmark_headings("[foo]: /url\nReal Setext\n---\n")
    assert len(h_ex) == 1 and h_ex[0]["kind"] == "setext" and h_ex[0]["level"] == 2 and h_ex[0]["rendered_title"] == "Real Setext"

    # Re-pinned semantic cell tamper tests (bad operator, placeholder, drifted threshold, wrong external_ref,
    # missing_data_disposition PASS, extra smuggled threshold, forbidden k==0_unconditional_veto)
    dt_p = min_root / "data" / "protocols" / "stage0_decision_tree_v1.json"
    orig_dt = dt_p.read_text(encoding="utf-8")
    vmod_sem = _load_mod("scripts/verify_iip0035_d1b_auto_and_portfolio.py", "openislet_v35_sem_test")
    orig_dt_sha = vmod_sem.EXPECTED_DEC114_DECISION_TREE_SHA256
    for mutate_blk in (
        lambda b: [c.__setitem__("operator", "BOGUS_OP") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-04"],
        lambda b: [c.__setitem__("threshold", "TBD") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-06"],
        lambda b: [c.__setitem__("threshold", c["threshold"].replace("8.0", "25.0")) for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-06"],
        lambda b: [c.__setitem__("threshold", c["threshold"] + ";override<=99.0") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-06"],
        lambda b: [c.__setitem__("external_ref", "DEC-114:D1-02") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-06"],
        lambda b: [c.__setitem__("missing_data_disposition", "PASS") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-06"],
        lambda b: [c.__setitem__("threshold", "k==0_unconditional_veto") for c in b["product_by_gate_contract"]["cells"] if c["product_variant_id"] == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN" and c["gate_id"] == "D1-04"],
    ):
        dt_mut = json.loads(orig_dt)
        blk_m = dt_mut["autologous_d1b_auto_c_v2_vanguard_gates_iip_0035"]
        mutate_blk(blk_m)
        dt_p.write_text(json.dumps(dt_mut, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        new_sha = hashlib.sha256(json.dumps(blk_m, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()
        try:
            vmod_sem.EXPECTED_DEC114_DECISION_TREE_SHA256 = new_sha
            with pytest.raises(ValueError):
                vmod_sem._resolve_and_verify_d1b_auto_v2_registry(min_root)
        finally:
            vmod_sem.EXPECTED_DEC114_DECISION_TREE_SHA256 = orig_dt_sha
    dt_p.write_text(orig_dt, encoding="utf-8")

    # Reversing D1-09 deadline inequality (<= 41.0 weeks -> >= 41.0 weeks) in decision-tree prose MUST fail predicate check!
    dt_rev = json.loads(orig_dt)
    blk_r = dt_rev["autologous_d1b_auto_c_v2_vanguard_gates_iip_0035"]
    blk_r["rules"]["d1_06_d1_09_kir3dl1_nrd1_hlac0304_and_cmc_turnaround"] = blk_r["rules"][
        "d1_06_d1_09_kir3dl1_nrd1_hlac0304_and_cmc_turnaround"
    ].replace("<= 41.0 weeks", ">= 41.0 weeks")
    dt_p.write_text(json.dumps(dt_rev, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    rev_sha = hashlib.sha256(json.dumps(blk_r, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()
    try:
        vmod_sem.EXPECTED_DEC114_DECISION_TREE_SHA256 = rev_sha
        with pytest.raises(ValueError, match="D1-09_turnaround_41_weeks"):
            vmod_sem._resolve_and_verify_d1b_auto_v2_registry(min_root)
    finally:
        vmod_sem.EXPECTED_DEC114_DECISION_TREE_SHA256 = orig_dt_sha
        dt_p.write_text(orig_dt, encoding="utf-8")

    # CRLF drift rejected on --check
    vmod_t = _load_mod("scripts/verify_iip0035_d1b_auto_and_portfolio.py", "openislet_v35_crlf_test")
    crlf_art = tmp_path / "v35_crlf.json"
    crlf_art.write_bytes(art_path.read_bytes().replace(b"\n", b"\r\n"))
    old_out, old_argv = vmod_t.DEFAULT_OUT_PATH, sys.argv
    try:
        vmod_t.DEFAULT_OUT_PATH = crlf_art
        sys.argv = [str(v_script), "--check"]
        with pytest.raises(SystemExit, match="byte/SHA-256 drift"):
            vmod_t.main()
    finally:
        vmod_t.DEFAULT_OUT_PATH = old_out
        sys.argv = old_argv


def test_05_sound_80digit_decimal_rare_event_and_standalone_gate_matrix() -> None:
    vmod = _load_mod("scripts/verify_iip0035_d1b_auto_and_portfolio.py", "openislet_v35_rare_test")

    # 1. Exact D1-04 (1e-8) boundary: 299,573,225 HOLD vs 299,573,226 PASS
    r_lo = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=299_573_225, max_u95_fraction=Decimal("1e-8")
    )
    r_hi = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=299_573_226, max_u95_fraction=Decimal("1e-8")
    )
    assert r_lo["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
    assert r_hi["status"] == "PASS"
    assert "MATHEMATICAL_COUNT_SUBGATE_ONLY" in r_hi["empirical_clearance_status"]

    # 2. Exact D1-08 (5.0%) boundary: 58 (5.0339%) HOLD vs 59 (4.9508%) PASS
    r_58 = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=58, max_u95_fraction=Decimal("0.05")
    )
    r_59 = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=59, max_u95_fraction=Decimal("0.05")
    )
    assert r_58["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
    assert r_58["u95_exact_sound_decimal_pct_str"].startswith("5.03393")
    assert r_59["status"] == "PASS"
    assert r_59["u95_exact_sound_decimal_pct_str"].startswith("4.95076")

    # 3. 120-digit straddle immune (sound .next_minus() on ln & exp!)
    with localcontext() as ctx120:
        ctx120.prec = 120
        p_true = Decimal(1) - (Decimal("0.05").ln() / Decimal(299_573_226)).exp()
        f_straddle = (p_true / Decimal("1e-8")) - Decimal("1e-60")
    r_strad = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=299_573_226, max_u95_fraction=Decimal("1e-8"),
        eta_recovery_l95=Decimal(1), f_eligible=f_straddle,
    )
    assert r_strad["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"

    # 4. Under-allocated alpha_eta (< 0.05 when eta < 1.0) rejected
    r_under = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=800_000_000, max_u95_fraction=Decimal("1e-8"),
        eta_recovery_l95=Decimal("0.50"), alpha_obs=Decimal("0.05"), alpha_eta_recovery=Decimal("0.02"),
    )
    assert r_under["status"] == "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION"

    # 5. Positive event k >= 1 is ALWAYS FAIL_HARD_STOP
    assert vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=1, n_raw=10**12, max_u95_fraction=Decimal("1e-8")
    )["status"] == "FAIL_HARD_STOP"

    # 6. Exported u95_exact_sound_decimal_str NEVER rounds down: reusing displayed bound as cap MUST PASS!
    for gid, n_val, cap in (
        ("D1_04", 299_573_226, Decimal("1e-8")),
        ("D1_08", 59, Decimal("0.05")),
        ("GATE0V", 5990, Decimal("0.0005")),
    ):
        res_b = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
            gate_id=gid, k=0, n_raw=n_val, max_u95_fraction=cap
        )
        reused = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
            gate_id=gid, k=0, n_raw=n_val, max_u95_fraction=Decimal(res_b["u95_exact_sound_decimal_str"])
        )
        assert reused["status"] == "PASS", f"Displayed u95 rounded below exact bound on {gid} n={n_val}: {res_b}"

    # 7. Bonferroni coverage label MUST floor (ROUND_FLOOR) inside prec=80 directed arithmetic:
    #    - 90.006% -> GE_90.00PCT
    #    - 0.04 + 0.05 + 1e-30 -> exact coverage 90.99999...% MUST floor to GE_90.99PCT, NEVER round up to 91.00%!
    #    - 0.05 + 0.05 + 1e-30 -> exceeds 0.10 MUST HOLD!
    r_floor = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=120, max_u95_fraction=Decimal("0.05"),
        eta_recovery_l95=Decimal("0.5"), alpha_obs=Decimal("0.04994"), alpha_eta_recovery=Decimal("0.05"),
    )
    assert r_floor["coverage_regime"] == "BONFERRONI_JOINT_COVERAGE_GE_90.00PCT_CALIBRATED"
    assert r_floor["alpha_joint_total_str"] == "0.09994"

    for tiny_af in (Decimal("1e-30"), Decimal("1e-500")):
        r_tiny = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
            gate_id="D1_08", k=0, n_raw=120, max_u95_fraction=Decimal("0.05"),
            eta_recovery_l95=Decimal("0.5"), f_eligible=Decimal("0.9"),
            f_eligible_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
            alpha_obs=Decimal("0.04"), alpha_eta_recovery=Decimal("0.05"), alpha_f_eligible=tiny_af,
        )
        assert r_tiny["coverage_regime"] == "BONFERRONI_JOINT_COVERAGE_GE_90.99PCT_CALIBRATED", r_tiny
        r_over = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
            gate_id="D1_08", k=0, n_raw=120, max_u95_fraction=Decimal("0.05"),
            eta_recovery_l95=Decimal("0.5"), f_eligible=Decimal("0.9"),
            f_eligible_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
            alpha_obs=Decimal("0.05"), alpha_eta_recovery=Decimal("0.05"), alpha_f_eligible=tiny_af,
        )
        assert r_over["status"] == "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION", r_over

    # 8. Conflicting nonzero alpha_f_eligible under KNOWN_OR_CONSERVATIVELY_ESTABLISHED MUST HOLD (never silently discard)!
    r_conf_af = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=120, max_u95_fraction=Decimal("0.05"),
        f_eligible_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
        alpha_f_eligible=Decimal("0.04"),
    )
    assert r_conf_af["status"] == "HOLD_CONFLICTING_F_ELIGIBLE_PROVENANCE_AND_ALPHA", r_conf_af

    # 9. Explicitly clustered / non-independent sampling MUST HOLD!
    r_clust = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=500, max_u95_fraction=Decimal("0.05"),
        clustered_sampling=True,
    )
    assert r_clust["status"] == "HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED", r_clust
    r_bad_indep = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=500_000_000, max_u95_fraction=Decimal("1e-8"),
        sampling_independence_model="CLUSTERED_REPLICATES",
    )
    assert r_bad_indep["status"] == "HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED", r_bad_indep

    # 10. Huge finite n_raw=10^100 MUST keep a strictly positive lower enclosure (zero 80-digit cancellation)!
    r_huge100 = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=10**100, max_u95_fraction=Decimal("0.05"),
    )
    assert r_huge100["status"] == "PASS"
    assert Decimal(r_huge100["u95_math_lower_pct_str"]) > Decimal(0), r_huge100

    # 11. V11 exact boundary equality (n=2, alpha=0.04, cap=0.8 -> 1 - 0.04^(1/2) == 0.8 exact -> PASS!)
    #     and exact upper-bound export u95_math_upper_pct_exact_str >= true U95!
    r_eq = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="ARITHMETIC", k=0, n_raw=2, alpha_obs=Decimal("0.04"), max_u95_fraction=Decimal("0.8"),
    )
    assert r_eq["status"] == "PASS", r_eq
    r_10e10 = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="ARITHMETIC", k=0, n_raw=10**10, max_u95_fraction=Decimal("0.05"),
    )
    assert Decimal(r_10e10["u95_math_upper_pct_exact_str"]) >= Decimal(r_10e10["u95_math_lower_pct_str"]) > Decimal(0)


def test_05b_table_driven_conformance_matrix_full_panel_helper_and_docs(tmp_path: Path) -> None:
    """Table-driven conformance test covering D1-05 missing/null/inconsistent inputs, lot conflicts,
    FAIL_HARD_STOP precedence over HOLD/incomplete fields, exact raw-JSON lexeme ingress (--eval-json),
    D1-08 unrounded lower-enclosure floor vs N=58 denominator HOLD, D1-03 head-to-head comparator,
    standalone L1B Gate 0V L95 + 72h early-titer veto & non-LIMS status, and locus-specific component cells.
    """
    vmod = _load_mod("scripts/verify_iip0035_d1b_auto_and_portfolio.py", "openislet_v35_conf_test")
    pc = _load_mod("scripts/product_contract.py", "openislet_pc_conf_test")
    by_id = {p["product_variant_id"]: p for p in pc.load_product_contracts(ROOT, require_staged_fastas=True)["products"]}
    p_l1a = by_id["ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN"]
    p_l1b = by_id["ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR"]

    base_l1a: dict[str, Any] = {
        "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
        "d1_01_4digit_hla_abc_kir_genotyping_verified": True,
        "d1_auto_hlac_sp3c4c_genotype_class": "ELIGIBLE_VL9_COMPETENT_HLAC",
        "d1_02_n1_phased_wgs_offgrid_guard_011_cleared": True,
        "obs_d1_02_min_mismatches_to_patient_hlac_e": 3,
        "obs_d1_02_wt_hlac_allele_retention_pct": 99.5,
        "obs_d1_02_hlac_bystander_ct_pct": 0.1,
        "obs_d1_02_chr6p21_interlocus_del_pct": 0.0,
        "d1_auto_cipsc_source_arm": "ADSC_SUN_PROTECTED_CIPSC",
        "d1_03_head_to_head_adsc_vs_blood_comparator_attested": True,
        "d1_03_adsc_vs_blood_bcor_tp53_ppm1d_wgs_clean": True,
        "obs_d1_03_chrm_numt_amplicon_depth_x": 1500.0,
        "obs_d1_03_chrm_novel_nonsyn_het_pct": 0.2,
        "d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified": True,
        "obs_d1_04_g0_beta_htmpk_leak_pct": 0.05,
        "obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct": 0.05,
        "obs_d1_04_dec053_ap1903_icasp9_survivor_frac": 0.0,
        "obs_d1_04_ap1903_survivor_count_k": 0,
        "obs_d1_04_ap1903_effective_trials_n": 299_573_226,
        "d1_04_azt_kill_measured_after_forced_chr19_loss": True,
        "obs_d1_04_dec053_azt_dividing_log10_kill": 4.5,
        "obs_d1_05_l3_day28_retention_pct": 88.0,
        "obs_d1_05_l4_day28_retention_pct": 86.0,
        "obs_d1_05_l3_l4_day28_retention_pct": 86.0,
        "obs_d1_05_fsp1_ferroptosis_survival_pct": 87.0,
        "obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": 4.2,
        "obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct": 5.0,
        "d1_08_patient_carries_hlac0304": True,
        "d1_08_hlac0304_positive_control_valid": True,
        "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 4.96,
        "obs_d1_08_hlac0304_lysed_events_k": 0,
        "obs_d1_08_hlac0304_scored_targets_n": 59,
        "obs_d1_08_t1d_serum_c5b9_u95_pct": 4.96,
        "obs_d1_08_t1d_serum_c5b9_pos_events_k": 0,
        "obs_d1_08_t1d_serum_c5b9_scored_cells_n": 59,
        "obs_d1_09_autologous_turnaround_weeks": 36.0,
    }
    base_l1b = {**base_l1a, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"}
    g0v_ok = {
        "mcb_lot_id": "LOT_A",
        "cvb_day14_log10_reduction": 3.8,
        "cvb_day14_log10_reduction_l95": 3.2,
        "cvb4_72h_early_titer_log10_pfu_u95": 1.5,
        "day28_vp1_ie1_positive_k": 0,
        "day28_vp1_ie1_cells_scored_n": 5990,
        "gate0v_verdict": "PASS",
    }

    # Table: (case_name, product, obs, lot_arg, g0v_arg, exp_met, exp_status_prefix, exp_exact_reason)
    conformance_table = [
        ("l1a_baseline_pass", p_l1a, base_l1a, "LOT_A", None, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"),
        ("l1b_baseline_pass_5990", p_l1b, base_l1b, "LOT_A", g0v_ok, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"),
        # D1-03 mandatory head-to-head ADSC vs blood comparator attestation
        ("d1_03_missing_head_to_head_HOLD", p_l1a, {k: v for k, v in base_l1a.items() if k != "d1_03_head_to_head_adsc_vs_blood_comparator_attested"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_03_HEAD_TO_HEAD_ADSC_VS_BLOOD_COMPARATOR_REQUIRED"),
        # D1-04 canonical raw_trials_n alias accepted and conflict rejected
        ("d1_04_canonical_raw_trials_n_PASS", p_l1a, {**{k: v for k, v in base_l1a.items() if k != "obs_d1_04_ap1903_effective_trials_n"}, "obs_d1_04_ap1903_raw_trials_n": 299_573_226}, "LOT_A", None, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"),
        ("d1_04_raw_vs_effective_conflict_FAIL", p_l1a, {**base_l1a, "obs_d1_04_ap1903_raw_trials_n": 300_000_000}, "LOT_A", None, False, "FAIL_", "FAIL_D1_04_CONFLICTING_RAW_AND_EFFECTIVE_TRIALS_N"),
        # D1-05 missing, null, UNAVAILABLE, and inconsistent aggregate min
        ("d1_05_absent_l3_l4_only_min86", p_l1a, {k: v for k, v in base_l1a.items() if k not in ("obs_d1_05_l3_day28_retention_pct", "obs_d1_05_l4_day28_retention_pct")}, "LOT_A", None, False, "HOLD_", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"),
        ("d1_05_null_l3_l4_only_min86", p_l1a, {**base_l1a, "obs_d1_05_l3_day28_retention_pct": None, "obs_d1_05_l4_day28_retention_pct": None}, "LOT_A", None, False, "HOLD_", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"),
        ("d1_05_unavailable_l3", p_l1a, {**base_l1a, "obs_d1_05_l3_day28_retention_pct": "UNAVAILABLE"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"),
        ("d1_05_inconsistent_min_99_vs_88_86", p_l1a, {**base_l1a, "obs_d1_05_l3_l4_day28_retention_pct": 99.0}, "LOT_A", None, False, "FAIL_", "FAIL_D1_05_INCONSISTENT_REPORTED_MIN_L3_L4_RETENTION"),
        # Lot-ID conflict binding (explicit arg A + Gate 0V A vs observations lot_id B or mcb_lot_id B)
        ("l1b_conflict_obs_lot_id_B", p_l1b, {**base_l1b, "lot_id": "LOT_B"}, "LOT_A", g0v_ok, False, "HOLD_", "HOLD_D1_CONFLICTING_LOT_ID_BINDING"),
        ("l1b_conflict_obs_mcb_lot_id_B", p_l1b, {**base_l1b, "mcb_lot_id": "LOT_B"}, "LOT_A", g0v_ok, False, "HOLD_", "HOLD_D1_CONFLICTING_LOT_ID_BINDING"),
        ("l1a_conflict_obs_lot_id_B", p_l1a, {**base_l1a, "lot_id": "LOT_B"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_CONFLICTING_LOT_ID_BINDING"),
        # L1B Gate 0V Day-14 L95 and 72h early-titer veto
        ("l1b_missing_day14_l95_HOLD", p_l1b, base_l1b, "LOT_A", {k: v for k, v in g0v_ok.items() if k != "cvb_day14_log10_reduction_l95"}, False, "HOLD_", "HOLD_D1_L1B_GATE0V_INCOMPLETE_RECORD"),
        ("l1b_day14_l95_below_3_FAIL", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "cvb_day14_log10_reduction_l95": "2.99"}, False, "FAIL_", "FAIL_D1_L1B_GATE0V_DAY14_L95_BELOW_3LOG_SWITCH_TO_L1A_IL1RN"),
        ("l1b_72h_early_titer_above_2_FAIL", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "cvb4_72h_early_titer_log10_pfu_u95": "2.01"}, False, "FAIL_", "FAIL_D1_L1B_GATE0V_EARLY_72H_TITER_VETO_SWITCH_TO_L1A_IL1RN"),
        # Known positive event k >= 1 MUST stay FAIL_HARD_STOP even with PENDING verdict or missing/incomplete fields
        ("l1b_k1_with_pending_verdict", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "day28_vp1_ie1_positive_k": 1, "gate0v_verdict": "PENDING"}, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN"),
        ("l1b_k1_with_missing_day14", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "day28_vp1_ie1_positive_k": 1, "cvb_day14_log10_reduction": None}, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN"),
        ("d1_04_k1_with_missing_l3_ret", p_l1a, {**base_l1a, "obs_d1_04_ap1903_survivor_count_k": 1, "obs_d1_05_l3_day28_retention_pct": None}, "LOT_A", None, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED"),
        ("d1_08_k1_with_missing_turnaround", p_l1a, {**base_l1a, "obs_d1_08_t1d_serum_c5b9_pos_events_k": 1, "obs_d1_09_autologous_turnaround_weeks": None}, "LOT_A", None, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_08_T1D_SERUM_C5B9_EVENT_DETECTED"),
        # Zero-event counts on both sides of EVERY boundary: full panel agrees with helper (HOLD below, PASS at)
        ("d1_04_zero_event_299573225_HOLD", p_l1a, {**base_l1a, "obs_d1_04_ap1903_effective_trials_n": 299_573_225}, "LOT_A", None, False, "HOLD_", "HOLD_D1_04_AP1903_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        ("d1_08_hlac03_playbook_5p0339_at_58_HOLD", p_l1a, {**base_l1a, "obs_d1_08_hlac0304_scored_targets_n": 58, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "5.0339"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_08_HLAC0304_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        ("d1_08_hlac03_fraudulent_4p9_at_58_FAIL", p_l1a, {**base_l1a, "obs_d1_08_hlac0304_scored_targets_n": 58, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "4.9"}, "LOT_A", None, False, "FAIL_", "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND"),
        ("d1_08_c5b9_zero_event_58_HOLD", p_l1a, {**base_l1a, "obs_d1_08_t1d_serum_c5b9_scored_cells_n": 58, "obs_d1_08_t1d_serum_c5b9_u95_pct": "5.0339"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_08_T1D_SERUM_C5B9_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        ("l1b_gate0v_zero_event_5989_HOLD", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "day28_vp1_ie1_cells_scored_n": 5989}, False, "HOLD_", "HOLD_D1_L1B_GATE0V_DAY28_INSUFFICIENT_DENOMINATOR"),
        # D1-08 unrounded mathematical bound (4.9507609894% > 4.95076098882269% PASSES; 4.945% FAILS; omitted HOLD)
        ("d1_08_exact_conservative_4p9507609894_PASS", p_l1a, {**base_l1a, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "4.9507609894", "obs_d1_08_t1d_serum_c5b9_u95_pct": "4.9507609894"}, "LOT_A", None, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"),
        ("d1_04_json_float_0p1_equals_str_0p1_PASS", p_l1a, {**base_l1a, "obs_d1_04_g0_beta_htmpk_leak_pct": 0.1, "obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct": "0.1"}, "LOT_A", None, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"),
        ("l1b_day14_exact_str_3p0_PASS", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "cvb_day14_log10_reduction": "3.0", "cvb_day14_log10_reduction_l95": "3.0"}, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"),
        ("l1b_day14_str_2p9999999999999999_FAIL", p_l1b, base_l1b, "LOT_A", {**g0v_ok, "cvb_day14_log10_reduction": "2.9999999999999999", "cvb_day14_log10_reduction_l95": "2.9999999999999999"}, False, "FAIL_", "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN"),
        ("l1b_day14_huge_int_no_overflow_with_d1_04_k1_FAIL_HARD_STOP", p_l1b, {**base_l1b, "obs_d1_04_ap1903_survivor_count_k": 1}, "LOT_A", {**g0v_ok, "cvb_day14_log10_reduction": 10**500}, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED"),
        ("d1_08_reported_4p945_below_4p9508_strict_FAIL", p_l1a, {**base_l1a, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 4.945}, "LOT_A", None, False, "FAIL_", "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND"),
        ("d1_08_omitted_u95_at_k0_n59_HOLD", p_l1a, {k: v for k, v in base_l1a.items() if k != "obs_d1_08_hlac0304_tcell_lysis_u95_pct"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_08_HLAC0304_REPORTED_U95_REQUIRED"),
        ("d1_08_calibrated_eta_0p5_at_n59_HOLD", p_l1a, {**base_l1a, "obs_d1_08_hlac0304_eta_recovery_l95": 0.5, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 9.91}, "LOT_A", None, False, "HOLD_", "HOLD_D1_08_HLAC0304_INSUFFICIENT_INDEPENDENT_DENOMINATOR"),
        ("d1_02_exact_decimal_1e_minus_400_FAIL", p_l1a, {**base_l1a, "obs_d1_02_chr6p21_interlocus_del_pct": "1e-400"}, "LOT_A", None, False, "FAIL_", "FAIL_D1_02_CHR6P21_INTERLOCUS_DELETION_DETECTED"),
        ("d1_06_exact_decimal_8p0000000000000001_FAIL", p_l1a, {**base_l1a, "obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": "8.0000000000000001"}, "LOT_A", None, False, "FAIL_", "FAIL_D1_06_KIR3DL1_MISSING_BW4_NK_LYSIS_EXCEEDS_8PCT"),
        ("l1a_missing_lot_id_HOLD", p_l1a, base_l1a, None, None, False, "HOLD_", "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"),
        ("l1a_tampered_product_def_FAIL", {**p_l1a, "execution_supported": True, "physical_route_id": "ROUTE_2_OPEN_PORE_PREVASC"}, base_l1a, "LOT_A", None, False, "FAIL_", "FAIL_D1_INVALID_OR_TAMPERED_PRODUCT_DEFINITION"),
        ("d1_08_reported_zero_u95_at_k0_n59_FAIL", p_l1a, {**base_l1a, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 0.0, "obs_d1_08_t1d_serum_c5b9_u95_pct": 0.0}, "LOT_A", None, False, "FAIL_", "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND"),
        ("malformed_categorical_array_with_k1_FAIL_HARD_STOP", p_l1a, {**base_l1a, "d1_auto_hlac_sp3c4c_genotype_class": ["BAD_LIST"], "obs_d1_04_ap1903_survivor_count_k": 1}, "LOT_A", None, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED"),
        # D1-08 two-branch allele-conditional conformance + positive-control validity
        ("d1_08_missing_carrier_status_HOLD", p_l1a, {k: v for k, v in base_l1a.items() if k != "d1_08_patient_carries_hlac0304"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED"),
        ("d1_08_c0304_pos_missing_pos_ctrl_HOLD", p_l1a, {k: v for k, v in base_l1a.items() if k != "d1_08_hlac0304_positive_control_valid"}, "LOT_A", None, False, "HOLD_", "HOLD_D1_08_HLAC0304_POSITIVE_CONTROL_VALIDITY_REQUIRED"),
        ("d1_08_c0304_pos_invalid_pos_ctrl_FAIL", p_l1a, {**base_l1a, "d1_08_hlac0304_positive_control_valid": False}, "LOT_A", None, False, "FAIL_", "FAIL_D1_08_HLAC0304_POSITIVE_CONTROL_INVALID"),
        (
            "d1_08_non_c0304_alternate_branch_PASS_no_c0304_counts",
            p_l1a,
            {
                **{k: v for k, v in base_l1a.items() if k not in ("d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")},
                "d1_08_patient_carries_hlac0304": False,
                "d1_08_matched_allele_positive_control_valid": True,
                "d1_08_non_c0304_matched_allele_screen_passed": True,
            },
            "LOT_A", None, True, "REFERENCE_", "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING",
        ),
        (
            "d1_08_non_c0304_missing_matched_allele_screen_HOLD",
            p_l1a,
            {
                **{k: v for k, v in base_l1a.items() if k not in ("d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")},
                "d1_08_patient_carries_hlac0304": False,
                "d1_08_matched_allele_positive_control_valid": True,
            },
            "LOT_A", None, False, "HOLD_", "HOLD_D1_08_MATCHED_ALLELE_CD8_SCREEN_REQUIRED",
        ),
        (
            "d1_08_non_c0304_failed_matched_allele_screen_FAIL",
            p_l1a,
            {
                **{k: v for k, v in base_l1a.items() if k not in ("d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")},
                "d1_08_patient_carries_hlac0304": False,
                "d1_08_matched_allele_positive_control_valid": True,
                "d1_08_non_c0304_matched_allele_screen_passed": False,
            },
            "LOT_A", None, False, "FAIL_", "FAIL_D1_08_MATCHED_ALLELE_CD8_SCREEN_FAILED",
        ),
        (
            "d1_08_non_c0304_with_unexpected_c0304_k1_FAIL_HARD_STOP",
            p_l1a,
            {
                **base_l1a,
                "d1_08_patient_carries_hlac0304": False,
                "d1_08_matched_allele_positive_control_valid": True,
                "d1_08_non_c0304_matched_allele_screen_passed": True,
                "obs_d1_08_hlac0304_lysed_events_k": 1,
            },
            "LOT_A", None, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED",
        ),
        # L1A optional Gate 0V record failing must report consistent status (NOT "NOT_APPLICABLE_L1A_NO_CD46"!)
        ("l1a_optional_gate0v_k1_consistent_FAIL_HARD_STOP", p_l1a, base_l1a, "LOT_A", {**g0v_ok, "day28_vp1_ie1_positive_k": 1}, False, "FAIL_HARD_STOP_", "FAIL_HARD_STOP_D1_L1A_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"),
    ]
    for cname, prod, obs, lot_arg, g0v_arg, exp_met, exp_pref, exp_st in conformance_table:
        out = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
            obs, product_def=prod, evaluated_lot_id=lot_arg, global_gate0v_record=g0v_arg
        )
        assert "gate_pass" not in out and "CLEARED" not in out["status"]
        assert out["clinical_or_batch_release_authorized"] is False
        assert out["controlling_contract_sha256"] == EXPECTED_DEC114_DT_SHA256
        assert out["controlling_playbook_sec4120_sha256"] == EXPECTED_DEC114_PB_SHA256
        assert out["panel_coverage_scope"] == "PER_ASSAY_INDIVIDUAL_COVERAGE_NOT_SIMULTANEOUS_PANEL_FAMILYWISE_GUARANTEE"
        assert out["controlling_contract_cells_loaded"] == 9
        if cname == "l1a_optional_gate0v_k1_consistent_FAIL_HARD_STOP":
            assert out["global_gate0v_standalone_subcheck_status"] != "NOT_APPLICABLE_L1A_NO_CD46"
            assert out["global_gate0v_lims_ledger_verified"] is False
        assert out["reference_criteria_met"] is exp_met, f"{cname}: reference_criteria_met={out['reference_criteria_met']} expected {exp_met} ({out})"
        assert out["status"].startswith(exp_pref), f"{cname}: status={out['status']!r} missing prefix {exp_pref!r}"
        assert out["status"] == exp_st, f"{cname}: status={out['status']!r} expected {exp_st!r} (reasons={out['reasons']})"

    # Verify Gate 0V standalone-vs-LIMS separation and full-panel coverage/alpha propagation!
    out_l1a = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_l1a["global_gate0v_standalone_subcheck_status"] == "NOT_APPLICABLE_L1A_NO_CD46"
    assert out_l1a["global_gate0v_lims_ledger_verified"] is None
    assert "global_gate0v_cleared" not in out_l1a
    assert out_l1a["evaluation_kind"] == "REFERENCE_CONTRACT_COMPUTATION_NOT_EMPIRICAL_RELEASE_DECISION"
    assert out_l1a["coverage_regime"] == "OBSERVATION_LEVEL_95PCT_UNIT_RECOVERY_AND_ELIGIBILITY_ASSUMED"
    assert out_l1a["rare_event_calibration_summary"]["d1_04_ap1903"]["alpha_joint_total_str"] == "0.05"
    out_l1b = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=g0v_ok)
    assert out_l1b["global_gate0v_standalone_subcheck_status"] == "STANDALONE_SUBCHECK_PASSED_LIVE_LIMS_LEDGER_UNVERIFIED"
    assert out_l1b["global_gate0v_lims_ledger_verified"] is False
    out_cal = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_08_hlac0304_eta_recovery_l95": "0.5", "obs_d1_08_hlac0304_scored_targets_n": 120, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "4.95"},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_cal["reference_criteria_met"] is True
    assert out_cal["coverage_regime"].startswith("MIXED_CALIBRATED_AND_UNIT_ASSUMED:")
    assert out_cal["rare_event_calibration_summary"]["d1_08_hlac0304"]["coverage_regime"] == "BONFERRONI_JOINT_COVERAGE_GE_90.00PCT_CALIBRATED"
    assert out_cal["rare_event_calibration_summary"]["d1_08_hlac0304"]["alpha_joint_total_str"] == "0.10"

    # Prove evaluator is genuinely driven by the 45-cell table (changing D1-06 threshold 8.0 -> 3.0 flips 4.2% to FAIL!)
    dt_blk = json.loads((ROOT / "data" / "protocols" / "stage0_decision_tree_v1.json").read_text(encoding="utf-8"))[
        "autologous_d1b_auto_c_v2_vanguard_gates_iip_0035"
    ]
    dt_tight = copy.deepcopy(dt_blk)
    for c in dt_tight["product_by_gate_contract"]["cells"]:
        if c["product_variant_id"] == p_l1a["product_variant_id"] and c["gate_id"] == "D1-06":
            c["threshold"] = "kir3dl1_missing_bw4_autologous_nk_lysis_pct<=3.0"
    out_tight = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", decision_tree_block=dt_tight, require_canonical_sha=False
    )
    assert out_tight["reference_criteria_met"] is False
    assert out_tight["status"] == "FAIL_D1_06_KIR3DL1_MISSING_BW4_NK_LYSIS_EXCEEDS_8PCT"
    assert out_tight["controlling_contract_canonical_sha_verified"] is False

    # Exact raw-JSON lexeme ingress (parse_stage0_observation_json_exact & CLI --eval-json):
    # Unquoted numeric JSON 8.0000000000000001 and 1e-400 MUST NOT round/underflow via binary64!
    raw_d106_json = json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": base_l1a})
    raw_d106_json = raw_d106_json.replace('"obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": 4.2', '"obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": 8.0000000000000001')
    parsed_106 = vmod.parse_stage0_observation_json_exact(raw_d106_json)
    out_106 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(parsed_106["observations"], product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_106["reference_criteria_met"] is False and out_106["status"] == "FAIL_D1_06_KIR3DL1_MISSING_BW4_NK_LYSIS_EXCEEDS_8PCT"

    raw_d102_json = json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": base_l1a})
    raw_d102_json = raw_d102_json.replace('"obs_d1_02_chr6p21_interlocus_del_pct": 0.0', '"obs_d1_02_chr6p21_interlocus_del_pct": 1e-400')
    parsed_102 = vmod.parse_stage0_observation_json_exact(raw_d102_json)
    out_102 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(parsed_102["observations"], product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_102["reference_criteria_met"] is False and out_102["status"] == "FAIL_D1_02_CHR6P21_INTERLOCUS_DELETION_DETECTED"

    ev_path = tmp_path / "eval_d106_exact.json"
    ev_path.write_text(raw_d106_json, encoding="utf-8")
    env_ev = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    res_ev = subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_path)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    )
    ev_out = json.loads(res_ev.stdout)
    assert ev_out["reference_criteria_met"] is False and ev_out["status"] == "FAIL_D1_06_KIR3DL1_MISSING_BW4_NK_LYSIS_EXCEEDS_8PCT"
    assert ev_out["controlling_contract_sha256"] == EXPECTED_DEC114_DT_SHA256
    assert ev_out["clinical_or_batch_release_authorized"] is False

    # --eval-json MUST reject conflicting product identity declarations (req.product_variant_id vs obs/product_def)!
    ev_conf = tmp_path / "eval_conflict_pvid.json"
    ev_conf.write_text(json.dumps({
        "product_variant_id": p_l1a["product_variant_id"],
        "product_def": p_l1b,
        "evaluated_lot_id": "LOT_A",
        "observations": base_l1a,
    }), encoding="utf-8")
    res_conf = subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_conf)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    )
    conf_out = json.loads(res_conf.stdout)
    assert conf_out["reference_criteria_met"] is False
    assert "FAIL_D1_CONFLICTING_PRODUCT_IDENTITY_DECLARATIONS" in conf_out["reasons"]

    # Direct Python evaluator MUST ALSO verify Playbook §4.1.20 SHA-256 (same integrity rule as --eval-json)!
    min_ev = _build_minimal_verifier_tree(tmp_path / "min_eval_tamper")
    pb_ev_p = min_ev / "docs" / "STAGE-0-EXECUTION-PLAYBOOK.md"
    orig_pb_ev = pb_ev_p.read_text(encoding="utf-8")
    pb_ev_p.write_text(orig_pb_ev.replace("<= 8.0%", "<= 99.0%"), encoding="utf-8")
    with pytest.raises(ValueError, match="Playbook section 4.1.20 SHA-256 mismatch"):
        vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", repo_root=min_ev)
    pb_ev_p.write_text(orig_pb_ev, encoding="utf-8")

    # --eval-json MUST fail closed if stage0_decision_tree_v1.json is tampered!
    dt_ev_p = min_ev / "data" / "protocols" / "stage0_decision_tree_v1.json"
    dt_ev_doc = json.loads(dt_ev_p.read_text(encoding="utf-8"))
    dt_ev_doc["autologous_d1b_auto_c_v2_vanguard_gates_iip_0035"]["decision_id"] = "DEC-TAMPERED"
    dt_ev_p.write_text(json.dumps(dt_ev_doc, indent=2) + "\n", encoding="utf-8")
    res_tamper = subprocess.run(
        [sys.executable, "-B", str(min_ev / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_path)],
        cwd=str(min_ev), env=env_ev, capture_output=True, text=True,
    )
    assert res_tamper.returncode != 0

    # V8 Parity & Hardening checks across direct Python and CLI --eval-json:
    # A. Matching-lot L1A Day-28 positive event inside observations["global_gate0v_record"] MUST FAIL_HARD_STOP on BOTH Python and CLI!
    obs_l1a_g0v_k1 = {**base_l1a, "global_gate0v_record": {**g0v_ok, "day28_vp1_ie1_positive_k": 1}}
    out_py_l1a_k1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_l1a_g0v_k1, product_def=p_l1a, evaluated_lot_id="LOT_A")
    ev_l1a_k1 = tmp_path / "eval_l1a_obs_g0v_k1.json"
    ev_l1a_k1.write_text(json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": obs_l1a_g0v_k1}), encoding="utf-8")
    out_cli_l1a_k1 = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_l1a_k1)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert out_py_l1a_k1["status"] == out_cli_l1a_k1["status"] == "FAIL_HARD_STOP_D1_L1A_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"
    assert out_py_l1a_k1["reasons"] == out_cli_l1a_k1["reasons"]

    # B. Passing canonical product_def (with float kinship_pi_hat_lt=0.05 parsed as Decimal by exact parser) MUST PASS on BOTH Python and CLI!
    ev_pdef = tmp_path / "eval_canonical_pdef.json"
    ev_pdef.write_text(json.dumps({"product_def": p_l1a, "evaluated_lot_id": "LOT_A", "observations": base_l1a}), encoding="utf-8")
    out_cli_pdef = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_pdef)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert out_cli_pdef["reference_criteria_met"] is True and out_cli_pdef["status"] == "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"

    # C. Tampered staged FASTA MUST fail closed on BOTH direct Python and CLI!
    min_fa = _build_minimal_verifier_tree(tmp_path / "min_fasta_tamper")
    fa_t = min_fa / "data" / "constructs" / "islet1_d1b_auto_v2_v0.8.0.fasta"
    fa_t.write_bytes(fa_t.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", repo_root=min_fa)
    assert subprocess.run(
        [sys.executable, "-B", str(min_fa / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_pdef)],
        cwd=str(min_fa), env=env_ev, capture_output=True, text=True,
    ).returncode != 0

    # D. Conflicting envelope lot alias (mcb_lot_id/lot_id) in nested CLI --eval-json MUST HOLD!
    ev_env_lot = tmp_path / "eval_envelope_lot_conflict.json"
    ev_env_lot.write_text(json.dumps({
        "product_variant_id": p_l1b["product_variant_id"],
        "evaluated_lot_id": "LOT_A",
        "mcb_lot_id": "LOT_CONFLICT_Z",
        "global_gate0v_record": g0v_ok,
        "observations": base_l1b,
    }), encoding="utf-8")
    out_env_lot = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_env_lot)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert out_env_lot["reference_criteria_met"] is False
    assert "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in out_env_lot["reasons"]

    # E. Clustered sampling rejected on D1-04, D1-08, AND L1B Gate-0V; huge-n=10^100 rejects reported_u95_pct=0!
    out_cl_d104 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_04_ap1903_clustered_sampling": True}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_cl_d104["reference_criteria_met"] is False and any("CLUSTERED_SAMPLING_REJECTED" in r for r in out_cl_d104["reasons"])
    out_cl_g0v = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record={**g0v_ok, "clustered_sampling": True})
    assert out_cl_g0v["reference_criteria_met"] is False and any("CLUSTERED_SAMPLING_REJECTED" in r for r in out_cl_g0v["reasons"])
    out_huge_z = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_08_t1d_serum_c5b9_scored_cells_n": 10**100, "obs_d1_08_t1d_serum_c5b9_u95_pct": 0.0}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_huge_z["reference_criteria_met"] is False and "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND" in out_huge_z["reasons"]

    # F. Clean component-only evaluation (L2, L3, L4): ZERO spurious inapplicable full-product HOLDs!
    p_l2 = by_id["ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1"]
    p_l3 = by_id["ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9"]
    p_l4 = by_id["ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y"]
    out_l3 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            "obs_d1_04_ap1903_survivor_count_k": 0,
            "obs_d1_04_ap1903_raw_trials_n": 299_573_226,
            "obs_d1_05_l3_day28_retention_pct": 88.0,
            "obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct": 5.0,
        },
        product_def=p_l3, evaluated_lot_id="LOT_L3",
    )
    assert out_l3["reference_criteria_met"] is False
    assert out_l3["status"] == "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"
    assert out_l3["component_subgate_evaluation"]["component_subgates_passed"] is True
    assert out_l3["reasons"] == ["HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"]
    assert not any("D1_03" in r or "D1_08" in r for r in out_l3["reasons"])
    out_l2 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            "d1_08_patient_carries_hlac0304": True,
            "d1_08_hlac0304_positive_control_valid": True,
            "obs_d1_08_hlac0304_lysed_events_k": 0,
            "obs_d1_08_hlac0304_scored_targets_n": 59,
            "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 4.96,
        },
        product_def=p_l2, evaluated_lot_id="LOT_L2",
    )
    assert out_l2["component_subgate_evaluation"]["component_subgates_passed"] is True
    assert out_l2["reasons"] == ["HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"]
    out_l4_fail = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            "d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified": True,
            "d1_04_azt_kill_measured_after_forced_chr19_loss": True,
            "obs_d1_04_g0_beta_htmpk_leak_pct": 0.5,
            "obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct": 0.05,
            "obs_d1_04_dec053_azt_dividing_log10_kill": 4.5,
            "obs_d1_05_l4_day28_retention_pct": 86.0,
            "obs_d1_05_fsp1_ferroptosis_survival_pct": 87.0,
        },
        product_def=p_l4, evaluated_lot_id="LOT_L4",
    )
    assert out_l4_fail["status"] == "FAIL_D1_04_G0_BETA_HTMPK_F105Y_LEAK_EXCEEDS_0P1PCT"
    assert not any("D1_03" in r or "D1_08" in r for r in out_l4_fail["reasons"])

    # G. Raw D1-01 patient_hla_c_alleles derivation (Lin 2023 Table 1 allowlist) + D1-02 COMPOUND_ALL_PASS + cross-gate D1-08 consistency + sampling conflict checks!
    out_raw_ok = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            **{k: v for k, v in base_l1a.items() if k != "d1_auto_hlac_sp3c4c_genotype_class"},
            "patient_hla_c_alleles": ["HLA-C*03:04", "HLA-C*07:01"],
            "obs_d1_02_min_mismatches_to_patient_hlac_e": 3,
        },
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_raw_ok["reference_criteria_met"] is True
    assert out_raw_ok["d1_01_allele_derivation_mode"] == "DERIVED_FROM_SUBMITTED_4DIGIT_HLA_C_ALLELES"
    assert out_raw_ok["d1_02_offgrid_guard_evidence_mode"] == "CALLER_ATTESTED_PHASED_WGS_AND_SUBMITTED_MISMATCH_COUNT_COMPOUND_ALL_PASS"
    out_raw_excl = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "patient_hla_c_alleles": ["HLA-C*07:01", "HLA-C*18:01"]},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert "FAIL_D1_01_ALLELE_CLASS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE" in out_raw_excl["reasons"]
    # Unclassified allele outside Lin 2023 Table 1 without a confirmed VL9-competent partner MUST HOLD!
    out_unclass = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            **{k: v for k, v in base_l1a.items() if k not in ("d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")},
            "patient_hla_c_alleles": ["HLA-C*07:01", "HLA-C*99:99"],
            "d1_08_patient_carries_hlac0304": False,
            "d1_08_matched_allele_positive_control_valid": True,
            "d1_08_non_c0304_matched_allele_screen_passed": True,
        },
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_unclass["reference_criteria_met"] is False
    assert "HOLD_D1_01_UNCLASSIFIED_HLA_C_SIGNAL_PEPTIDE_ALLELE_REQUIRES_IPD_IMGT_VERIFICATION" in out_unclass["reasons"]
    # Cross-gate D1-01/D1-08 contradiction: raw HLA-C*03:04 carrier with d1_08_patient_carries_hlac0304=False MUST FAIL on BOTH Python and CLI!
    obs_c0304_skip = {
        **{k: v for k, v in base_l1a.items() if k not in ("d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")},
        "patient_hla_c_alleles": ["HLA-C*04:01", "HLA-C*03:04"],
        "d1_08_patient_carries_hlac0304": False,
        "d1_08_matched_allele_positive_control_valid": True,
        "d1_08_non_c0304_matched_allele_screen_passed": True,
    }
    out_py_c0304_skip = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_c0304_skip, product_def=p_l1a, evaluated_lot_id="LOT_A")
    ev_c0304_skip = tmp_path / "eval_c0304_skip.json"
    ev_c0304_skip.write_text(json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": obs_c0304_skip}), encoding="utf-8")
    out_cli_c0304_skip = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_c0304_skip)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert out_py_c0304_skip["reference_criteria_met"] is False and out_cli_c0304_skip["reference_criteria_met"] is False
    assert "FAIL_D1_08_HLAC0304_CARRIER_STATUS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE" in out_py_c0304_skip["reasons"]
    assert out_py_c0304_skip["reasons"] == out_cli_c0304_skip["reasons"]
    # Reverse cross-gate contradiction (non-C*03:04 alleles claiming d1_08_patient_carries_hlac0304=True) MUST also FAIL!
    out_rev_c0304 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "patient_hla_c_alleles": ["HLA-C*04:01", "HLA-C*06:02"], "d1_08_patient_carries_hlac0304": True},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_rev_c0304["reference_criteria_met"] is False
    assert "FAIL_D1_08_HLAC0304_CARRIER_STATUS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE" in out_rev_c0304["reasons"]
    # L2 component-only path MUST also enforce patient_hla_c_alleles vs d1_08_patient_carries_hlac0304!
    out_l2_c0304_skip = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            "patient_hla_c_alleles": ["HLA-C*04:01", "HLA-C*03:04"],
            "d1_08_patient_carries_hlac0304": False,
            "d1_08_matched_allele_positive_control_valid": True,
            "d1_08_non_c0304_matched_allele_screen_passed": True,
        },
        product_def=p_l2, evaluated_lot_id="LOT_L2",
    )
    assert "FAIL_D1_08_HLAC0304_CARRIER_STATUS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE" in out_l2_c0304_skip["reasons"]
    # Contradictory panel-level clustered_sampling=True vs assay-level False MUST HOLD on BOTH Python & CLI (plus L1B Gate-0V)!
    obs_samp_conf = {
        **base_l1a,
        "clustered_sampling": True,
        "obs_d1_04_ap1903_clustered_sampling": False,
        "obs_d1_08_hlac0304_clustered_sampling": False,
        "obs_d1_08_t1d_serum_c5b9_clustered_sampling": False,
    }
    out_py_samp = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_samp_conf, product_def=p_l1a, evaluated_lot_id="LOT_A")
    ev_samp = tmp_path / "eval_samp_conflict.json"
    ev_samp.write_text(json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": obs_samp_conf}), encoding="utf-8")
    out_cli_samp = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_samp)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert out_py_samp["reference_criteria_met"] is False and out_cli_samp["reference_criteria_met"] is False
    assert any("CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS" in r for r in out_py_samp["reasons"])
    assert out_py_samp["reasons"] == out_cli_samp["reasons"]
    out_l1b_samp = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1b, "clustered_sampling": True},
        product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record={**g0v_ok, "clustered_sampling": False},
    )
    assert out_l1b_samp["reference_criteria_met"] is False
    assert any("CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS" in r for r in out_l1b_samp["reasons"])
    out_raw_mm2 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_02_min_mismatches_to_patient_hlac_e": 2},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert "FAIL_D1_02_N1_PHASED_WGS_OFFGRID_GUARD_011_HLAC_E_MISMATCH_CONFLICT" in out_raw_mm2["reasons"]

    # H. V10 Red-Team Regression Guards:
    # 1. D1-01 raw SP-3C/SP-4C nonfunctional pair MUST FAIL (not HOLD or false mismatch) when class is omitted, UNAVAILABLE, or canonical!
    for cls_val in (None, "UNAVAILABLE", "EXCLUDED_SP3C4C_HLAC07_17_1801_NONFUNCTIONAL_HOMO_OR_COMPOUND_HET_ROUTE_TO_ALLO_Y_OR_D1A"):
        obs_sp34 = {
            **{k: v for k, v in base_l1a.items() if k not in ("d1_auto_hlac_sp3c4c_genotype_class", "d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")},
            "patient_hla_c_alleles": ["HLA-C*07:01", "HLA-C*18:01"],
            "d1_08_patient_carries_hlac0304": False,
            "d1_08_matched_allele_positive_control_valid": True,
            "d1_08_non_c0304_matched_allele_screen_passed": True,
        }
        if cls_val is not None:
            obs_sp34["d1_auto_hlac_sp3c4c_genotype_class"] = cls_val
        out_sp34 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_sp34, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert out_sp34["status"] == "FAIL_D1_01_EXCLUDED_AH8P1_SP3C4C_NONFUNCTIONAL_HLAC_ROUTE_TO_ALLO_Y_OR_D1A", (cls_val, out_sp34)
        assert "FAIL_D1_01_ALLELE_CLASS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE" not in out_sp34["reasons"]

    # 2. sampling_design="independent" MUST PASS (no "depend" substring bug!), whereas "dependent" rejects clustered/non-independent and None rejects malformed!
    out_sd_ok = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": "independent"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_sd_ok["reference_criteria_met"] is True, out_sd_ok
    out_sd_dep = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": "dependent"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_sd_dep["reference_criteria_met"] is False and any("CLUSTERED_SAMPLING_REJECTED" in r for r in out_sd_dep["reasons"]), out_sd_dep
    out_sd_none = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": None}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_sd_none["reference_criteria_met"] is False and any("MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION" in r for r in out_sd_none["reasons"]), out_sd_none

    # 3. D1-05 measured silencing MUST FAIL even when sibling locus or aggregate min is missing/None!
    out_d105_sil = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**{k: v for k, v in base_l1a.items() if k != "obs_d1_05_l3_l4_day28_retention_pct"}, "obs_d1_05_l3_day28_retention_pct": 50.0, "obs_d1_05_l4_day28_retention_pct": None},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_d105_sil["status"] == "FAIL_D1_05_L3_STAGE7_DAY28_RETENTION_BELOW_80PCT", out_d105_sil
    out_d105_min_om = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**{k: v for k, v in base_l1a.items() if k != "obs_d1_05_l3_l4_day28_retention_pct"}, "obs_d1_05_l3_day28_retention_pct": 70.0, "obs_d1_05_l4_day28_retention_pct": 86.0},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_d105_min_om["status"] == "FAIL_D1_05_L3_STAGE7_DAY28_RETENTION_BELOW_80PCT", out_d105_min_om
    assert "FAIL_D1_05_L3_L4_STAGE7_DAY28_RETENTION_BELOW_80PCT" in out_d105_min_om["reasons"], out_d105_min_om

    # 4. global_gate0v_record conflict: 0 != False strictly typed AND k>=1 hard-stop preserved across EITHER conflicting copy!
    out_g0v_bool = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1b, "global_gate0v_record": {**g0v_ok, "day28_vp1_ie1_positive_k": False}},
        product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=g0v_ok,
    )
    assert "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS" in out_g0v_bool["reasons"], out_g0v_bool
    out_g0v_conf_k1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1b, "global_gate0v_record": {**g0v_ok, "day28_vp1_ie1_positive_k": 1}},
        product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=g0v_ok,
    )
    assert out_g0v_conf_k1["status"].startswith("FAIL_HARD_STOP_"), out_g0v_conf_k1

    # 5. Out-of-range AND below-threshold Gate-0V Day-14 fields with missing 72h field still FAIL; non-dict observations fails closed; CLI rejects rogue envelope keys & padded pvid!
    out_g0v_oor = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A",
        global_gate0v_record={k: v for k, v in {**g0v_ok, "cvb_day14_log10_reduction": "-1.0"}.items() if k != "cvb4_72h_early_titer_log10_pfu_u95"},
    )
    assert out_g0v_oor["status"] == "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE", out_g0v_oor
    out_g0v_low = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A",
        global_gate0v_record={k: v for k, v in {**g0v_ok, "cvb_day14_log10_reduction": "2.0", "cvb_day14_log10_reduction_l95": "1.5"}.items() if k != "cvb4_72h_early_titer_log10_pfu_u95"},
    )
    assert out_g0v_low["status"] == "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN", out_g0v_low
    out_nondict = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035("not_a_dict", product_def=p_l1a, evaluated_lot_id="LOT_A")  # type: ignore[arg-type]
    assert out_nondict["reference_criteria_met"] is False and "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION" in out_nondict["reasons"]
    ev_rogue = tmp_path / "eval_rogue_key.json"
    ev_rogue.write_text(json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": base_l1a, "rogue_key": 1}), encoding="utf-8")
    assert subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_rogue)], cwd=str(ROOT), env=env_ev, capture_output=True).returncode != 0
    ev_pad = tmp_path / "eval_pad_pvid.json"
    ev_pad.write_text(json.dumps({"product_variant_id": " " + p_l1a["product_variant_id"] + " ", "evaluated_lot_id": "LOT_A", "observations": base_l1a}), encoding="utf-8")
    out_pad = json.loads(subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(ev_pad)], cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True).stdout)
    assert out_pad["reference_criteria_met"] is False and "HOLD_D1_CANONICAL_PRODUCT_VARIANT_ID_REQUIRED" in out_pad["reasons"]

    # Cross-check Playbook s4.1.20 & typed decision-tree cells (HOLD for 299,573,225 / 58 + locus-specific COMPONENT_ONLY cells)
    pb_txt = (ROOT / "docs" / "STAGE-0-EXECUTION-PLAYBOOK.md").read_text(encoding="utf-8")
    assert "299,573,225` HOLD" in pb_txt and "N = 58` HOLD" in pb_txt
    for cell in dt_blk["product_by_gate_contract"]["cells"]:
        assert "299573225_FAIL" not in cell["threshold"] and "N=58_FAIL" not in cell["threshold"]
        if cell["disposition"] == "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK":
            assert "L3&L4" not in cell["threshold"] and "C5b-9_U95<=5.0%" not in cell["threshold"]

    # I. V11 Blind-Review & Audit Closure Guards:
    # 1. CommonMark 0.31.2 official examples 79, 93, 216, 278 exact heading boundary conformance
    assert vmod._parse_commonmark_headings("> foo\nbar\n===\n") == []
    assert vmod._parse_commonmark_headings("-\n  foo\n-\n  ```\n  bar\n  ```\n-\n      baz\n") == []
    assert vmod._parse_commonmark_headings("[foo]: /url\n===\n") == []
    assert [h["rendered_title"] for h in vmod._parse_commonmark_headings("## \n#\n### ###\n")] == ["", "", ""]

    # 2. Strict canonical sampling_design allowlist & provenance (rejects 'not independent' / 'dependent' as non-independent, None as malformed; accepts 'independent biological replicates')
    for bad_phrase in ("not independent", "not independent biological opportunities", "dependent"):
        out_bp = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": bad_phrase}, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert out_bp["reference_criteria_met"] is False and any("CLUSTERED_SAMPLING_REJECTED" in r for r in out_bp["reasons"]), (bad_phrase, out_bp)
    out_bp_none = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": None}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_bp_none["reference_criteria_met"] is False and any("MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION" in r for r in out_bp_none["reasons"]), out_bp_none
    out_ibr = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": "independent biological replicates"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_ibr["reference_criteria_met"] is True
    assert out_ibr["rare_event_calibration_summary"]["d1_04_ap1903"]["sampling_independence_provenance"] == "EXPLICIT_INDEPENDENT"
    assert out_l1a["rare_event_calibration_summary"]["d1_04_ap1903"]["sampling_independence_provenance"] == "UNSPECIFIED_DEFAULT_ASSUMED_INDEPENDENT"

    # 3. Direct Python auto-resolves product_def by product_variant_id identically to --eval-json CLI; CLI rejects 'observations': null!
    id_only_obs = {"product_variant_id": p_l1a["product_variant_id"]}
    id_only_path = tmp_path / "id_only.json"
    id_only_path.write_text(json.dumps(id_only_obs), encoding="utf-8")
    out_cli_id = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(id_only_path)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    out_py_id = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(id_only_obs)
    assert out_py_id["controlling_contract_cells_loaded"] == out_cli_id["controlling_contract_cells_loaded"] == 9
    assert out_py_id["status"] == out_cli_id["status"] == "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
    out_py_full = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "product_variant_id": p_l1a["product_variant_id"]}, evaluated_lot_id="LOT_A")
    assert out_py_full["reference_criteria_met"] is True
    null_obs_path = tmp_path / "null_obs.json"
    null_obs_path.write_text(json.dumps({"product_variant_id": p_l1a["product_variant_id"], "observations": None}), encoding="utf-8")
    assert subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(null_obs_path)], cwd=str(ROOT), env=env_ev, capture_output=True).returncode != 0
    # On-disk decision tree tampering MUST fail closed even if caller passes a valid decision_tree_block with require_canonical_sha=True!
    with pytest.raises(ValueError, match="On-disk controlling DEC-114 decision-tree SHA-256 mismatch"):
        vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", decision_tree_block=dt_blk, repo_root=min_ev, require_canonical_sha=True)

    # 4. Gate-0V k>=1 hard stop NEVER suppressed by lot-alias conflict or MCB lot mismatch, AND standalone status agrees on dual conflict in BOTH orderings!
    out_lot_k1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1b, "lot_id": "LOT_B"}, product_def=p_l1b, evaluated_lot_id="LOT_A",
        global_gate0v_record={**g0v_ok, "day28_vp1_ie1_positive_k": 1},
    )
    assert out_lot_k1["status"].startswith("FAIL_HARD_STOP_") and "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in out_lot_k1["reasons"]
    out_mcb_k1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A",
        global_gate0v_record={**g0v_ok, "mcb_lot_id": "LOT_OTHER", "day28_vp1_ie1_positive_k": 1},
    )
    assert out_mcb_k1["status"].startswith("FAIL_HARD_STOP_") and "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH" in out_mcb_k1["reasons"]
    for g_kw, g_in in ((g0v_ok, {**g0v_ok, "day28_vp1_ie1_positive_k": 1}), ({**g0v_ok, "day28_vp1_ie1_positive_k": 1}, g0v_ok)):
        out_dual = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1b, "global_gate0v_record": g_in}, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=g_kw)
        assert out_dual["status"].startswith("FAIL_HARD_STOP_")
        assert out_dual["global_gate0v_standalone_subcheck_status"].startswith("FAIL_HARD_STOP_"), out_dual

    # 5. D1-02 mismatch upper bound (10^12 rejected), D1-03 deep depth (>100k accepted), D1-04 forced chr19 loss enforced on L1A & L4, D1-08 sub-bound rejected!
    out_mm_huge = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_02_min_mismatches_to_patient_hlac_e": 10**12}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION" in out_mm_huge["reasons"]
    out_depth_deep = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_03_chrm_numt_amplicon_depth_x": 150000.0}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_depth_deep["reference_criteria_met"] is True
    out_no_chr19_loss = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "d1_04_azt_kill_measured_after_forced_chr19_loss": False}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_no_chr19_loss["status"] == "FAIL_D1_04_AZT_KILL_NOT_MEASURED_AFTER_FORCED_CHR19_LOSS"
    out_l4_no_loss = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({
        "d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified": True,
        "d1_04_azt_kill_measured_after_forced_chr19_loss": False,
        "obs_d1_04_g0_beta_htmpk_leak_pct": 0.05, "obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct": 0.05,
        "obs_d1_04_dec053_azt_dividing_log10_kill": 4.5, "obs_d1_05_l4_day28_retention_pct": 86.0, "obs_d1_05_fsp1_ferroptosis_survival_pct": 87.0,
    }, product_def=p_l4, evaluated_lot_id="LOT_L4")
    assert out_l4_no_loss["status"] == "FAIL_D1_04_AZT_KILL_NOT_MEASURED_AFTER_FORCED_CHR19_LOSS"
    out_u95_sub = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "4.9507609887"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_u95_sub["status"] == "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND"
    # L2:D1-08 is self-contained (tightening L2's own cell fails L2 without touching L1A)
    dt_l2_tight = copy.deepcopy(dt_blk)
    for c in dt_l2_tight["product_by_gate_contract"]["cells"]:
        if c["product_variant_id"] == p_l2["product_variant_id"] and c["gate_id"] == "D1-08":
            c["threshold"] = c["threshold"].replace("U95<=5.0%", "U95<=2.0%")
    out_l2_tight = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {"d1_08_patient_carries_hlac0304": True, "d1_08_hlac0304_positive_control_valid": True, "obs_d1_08_hlac0304_lysed_events_k": 0, "obs_d1_08_hlac0304_scored_targets_n": 59, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 4.96},
        product_def=p_l2, evaluated_lot_id="LOT_L2", decision_tree_block=dt_l2_tight, require_canonical_sha=False,
    )
    assert out_l2_tight["component_subgate_evaluation"]["component_subgates_passed"] is False

    # J. Strict null-derivation blocking, Python/CLI conflict parity, on-disk contract verification accuracy & high-N exact straddle:
    # 1. Explicit "d1_08_patient_carries_hlac0304": null with valid raw alleles MUST block derivation and HOLD on L1A, L1B, and L2 (Python + CLI), while omitted key derives cleanly!
    raw_noncarrier = ["HLA-C*04:01", "HLA-C*06:02"]
    noncarrier_d108 = {
        "patient_hla_c_alleles": raw_noncarrier,
        "d1_08_matched_allele_positive_control_valid": True,
        "d1_08_non_c0304_matched_allele_screen_passed": True,
    }
    _strip_c0304 = ("d1_08_patient_carries_hlac0304", "d1_08_hlac0304_positive_control_valid", "obs_d1_08_hlac0304_tcell_lysis_u95_pct", "obs_d1_08_hlac0304_lysed_events_k", "obs_d1_08_hlac0304_scored_targets_n")
    for _pdef_nc, _base_nc, _g0v_nc, _is_comp_nc in (
        (p_l1a, {k: v for k, v in base_l1a.items() if k not in _strip_c0304}, None, False),
        (p_l1b, {k: v for k, v in base_l1b.items() if k not in _strip_c0304}, g0v_ok, False),
        (p_l2, {}, None, True),
    ):
        obs_omit = {**_base_nc, **noncarrier_d108}
        obs_null = {**obs_omit, "d1_08_patient_carries_hlac0304": None}
        out_omit = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_omit, product_def=_pdef_nc, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_nc)
        out_null = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_null, product_def=_pdef_nc, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_nc)
        if _is_comp_nc:
            assert out_omit["component_subgate_evaluation"]["component_subgates_passed"] is True, out_omit
            assert out_null["component_subgate_evaluation"]["component_subgates_passed"] is False, out_null
            assert "HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED" in out_null["component_subgate_evaluation"]["component_reasons"], out_null
        else:
            assert out_omit["reference_criteria_met"] is True, out_omit
            assert out_null["reference_criteria_met"] is False, out_null
            assert out_null["status"] == "HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED", out_null
        # Verify --eval-json CLI matches Python identically for explicit null!
        cli_null_p = tmp_path / f"cli_null_{_pdef_nc['product_variant_id']}.json"
        cli_env_obj: dict[str, object] = {"product_variant_id": _pdef_nc["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": obs_null}
        if _g0v_nc is not None:
            cli_env_obj["global_gate0v_record"] = _g0v_nc
        cli_null_p.write_text(json.dumps(cli_env_obj), encoding="utf-8")
        out_null_cli = json.loads(subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(cli_null_p)],
            cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
        ).stdout)
        assert out_null_cli["status"] == out_null["status"], (out_null_cli, out_null)
        assert out_null_cli["reasons"] == out_null["reasons"], (out_null_cli, out_null)

    # 2. Conflicting product IDs (L2 envelope/kwarg vs L1A in observations) with D1-04 AP1903 k=1 survivor MUST emit IDENTICAL FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED across Python and --eval-json CLI!
    conflict_obs_k1 = {**base_l1a, "product_variant_id": p_l1a["product_variant_id"], "obs_d1_04_ap1903_survivor_count_k": 1}
    out_conf_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        conflict_obs_k1, product_variant_id=p_l2["product_variant_id"], evaluated_lot_id="LOT_A",
    )
    out_conf_pdef = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        conflict_obs_k1, product_def=p_l2, evaluated_lot_id="LOT_A",
    )
    conf_cli_p = tmp_path / "conflict_l2_l1a_k1.json"
    conf_cli_p.write_text(json.dumps({
        "product_variant_id": p_l2["product_variant_id"],
        "evaluated_lot_id": "LOT_A",
        "observations": conflict_obs_k1,
    }), encoding="utf-8")
    out_conf_cli = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(conf_cli_p)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    for _co in (out_conf_py, out_conf_pdef, out_conf_cli):
        assert _co["status"] == "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED", _co
        assert "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED" in _co["reasons"], _co
        assert "FAIL_D1_CONFLICTING_PRODUCT_IDENTITY_DECLARATIONS" in _co["reasons"], _co
    assert out_conf_py["status"] == out_conf_cli["status"] == out_conf_pdef["status"]
    assert out_conf_py["reasons"] == out_conf_cli["reasons"] == out_conf_pdef["reasons"]

    # 3. require_canonical_sha=False with deleted on-disk stage0_decision_tree_v1.json MUST NOT report controlling_contract_canonical_sha_verified=True!
    missing_dt_root = tmp_path / "missing_dt_root"
    shutil.copytree(min_ev, missing_dt_root)
    (missing_dt_root / "data" / "protocols" / "stage0_decision_tree_v1.json").unlink()
    out_no_disk_dt = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A",
        decision_tree_block=dt_blk, repo_root=missing_dt_root, require_canonical_sha=False,
    )
    assert out_no_disk_dt["controlling_contract_canonical_sha_verified"] is False, out_no_disk_dt
    assert out_no_disk_dt["controlling_contract_on_disk_verified"] is False, out_no_disk_dt
    assert out_no_disk_dt["controlling_contract_in_memory_matches_canonical_sha"] is True, out_no_disk_dt

    # 4. Exact rational straddle equality at n = 1,001 (alpha = 0.95^1001, target = 0.05) MUST decide True on exact equality and False just below!
    import decimal as _dec
    with _dec.localcontext() as _dctx:
        _dctx.prec = 2500
        _alpha_1001 = _dec.Decimal("0.95") ** 1001
        _target_below = _dec.Decimal("0.05") - _dec.Decimal("1e-50")
    assert _target_below < _dec.Decimal("0.05")
    assert vmod._zero_event_u95_le_target(
        _alpha_1001, 1001, _dec.Decimal("1"), _dec.Decimal("1"), _dec.Decimal("0.05"),
        u95_lower=_dec.Decimal("0.049999"), u95_upper=_dec.Decimal("0.050001"),
    ) is True
    assert vmod._zero_event_u95_le_target(
        _alpha_1001, 1001, _dec.Decimal("1"), _dec.Decimal("1"), _target_below,
        u95_lower=_dec.Decimal("0.049999"), u95_upper=_dec.Decimal("0.050001"),
    ) is False

    # 5. Fullwidth Unicode digits & whitespace-padded HLA-C alleles MUST fail closed across L1A/L1B/L2 (Python + CLI), and fullwidth numeric strings MUST fail!
    for _bad_alleles in (["HLA-C*04:01", "HLA-C*\uff10\uff13:\uff10\uff14"], ["HLA-C*04:01", " HLA-C*03:04 "]):
        for _pdef_u, _base_u, _g0v_u, _is_comp_u in (
            (p_l1a, {k: v for k, v in base_l1a.items() if k not in _strip_c0304}, None, False),
            (p_l1b, {k: v for k, v in base_l1b.items() if k not in _strip_c0304}, g0v_ok, False),
            (p_l2, {}, None, True),
        ):
            obs_u = {**_base_u, **noncarrier_d108, "patient_hla_c_alleles": _bad_alleles}
            out_u_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_u, product_def=_pdef_u, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_u)
            if _is_comp_u:
                assert out_u_py["component_subgate_evaluation"]["component_subgates_passed"] is False, out_u_py
            else:
                assert out_u_py["reference_criteria_met"] is False, out_u_py
            assert "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION" in out_u_py["reasons"], out_u_py
            assert "HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED" in out_u_py["reasons"], out_u_py
            cli_u_p = tmp_path / f"cli_unicode_{_pdef_u['product_variant_id']}_{abs(hash(tuple(_bad_alleles)))}.json"
            cli_u_env: dict[str, object] = {"product_variant_id": _pdef_u["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": obs_u}
            if _g0v_u is not None:
                cli_u_env["global_gate0v_record"] = _g0v_u
            cli_u_p.write_text(json.dumps(cli_u_env, ensure_ascii=False), encoding="utf-8")
            out_u_cli = json.loads(subprocess.run(
                [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(cli_u_p)],
                cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
            ).stdout)
            assert out_u_cli["status"] == out_u_py["status"] and out_u_cli["reasons"] == out_u_py["reasons"], (out_u_cli, out_u_py)
    out_fw_num = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_08_t1d_serum_c5b9_u95_pct": "\uff14.\uff19\uff16"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION" in out_fw_num["reasons"], out_fw_num

    # 6. Exact 59-trial reported-U95 floor boundary (80th-decimal exact rational check: value strictly above 1-0.05^(1/59) PASSES, value strictly below FAILS across Python + CLI)!
    rep_59_above = "4.9507609888226947385979516996770924243183522678215692561372644501801537712135560"
    rep_59_below = "4.9507609888226947385979516996770924243183522678215692561372644501801537712135553"
    out_59_ok = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_08_t1d_serum_c5b9_u95_pct": rep_59_above, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": rep_59_above},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_59_ok["reference_criteria_met"] is True, out_59_ok
    out_59_bad = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_08_t1d_serum_c5b9_u95_pct": rep_59_below},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert out_59_bad["status"] == "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND", out_59_bad
    cli_59_p = tmp_path / "cli_59_exact.json"
    cli_59_p.write_text(json.dumps({
        "product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A",
        "observations": {**base_l1a, "obs_d1_08_t1d_serum_c5b9_u95_pct": rep_59_above, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": rep_59_above},
    }), encoding="utf-8")
    out_59_cli = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(cli_59_p)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert out_59_cli["reference_criteria_met"] is True, out_59_cli

    # 7. Exhaustive contract-first 45-cell omission/null/UNAVAILABLE sweep across ALL required fields on L1A, L1B (+Gate-0V), L2, L3, and L4!
    #    Explicitly confirms D1-02 COMPOUND_ALL_PASS (omitting EITHER phased_wgs OR min_mismatches HOLDS on BOTH Python & CLI)
    #    and D1-01 explicit null blocks allele derivation!
    _d101_null_obs = {**base_l1a, "patient_hla_c_alleles": ["HLA-C*03:04", "HLA-C*07:01"], "d1_auto_hlac_sp3c4c_genotype_class": None}
    out_d101_null = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_d101_null_obs, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert out_d101_null["reference_criteria_met"] is False and "HOLD_D1_01_MISSING_GENOTYPING_OR_UNKNOWN_HLAC_CLASS" in out_d101_null["reasons"], out_d101_null
    assert out_d101_null["d1_01_allele_derivation_mode"] == "EXPLICIT_NULL_OR_UNAVAILABLE_BLOCKED_DERIVATION_HOLD"
    _dt_cells_by_key = {
        (c["product_variant_id"], c["gate_id"]): c
        for c in json.loads((ROOT / "data" / "protocols" / "stage0_decision_tree_v1.json").read_text(encoding="utf-8"))[
            "autologous_d1b_auto_c_v2_vanguard_gates_iip_0035"
        ]["product_by_gate_contract"]["cells"]
    }
    _req_full_fields: list[tuple[str, str, str]] = [
        ("d1_01_4digit_hla_abc_kir_genotyping_verified", "D1-01", "HOLD_D1_01_MISSING_GENOTYPING_OR_UNKNOWN_HLAC_CLASS"),
        ("d1_auto_hlac_sp3c4c_genotype_class", "D1-01", "HOLD_D1_01_MISSING_GENOTYPING_OR_UNKNOWN_HLAC_CLASS"),
        ("d1_02_n1_phased_wgs_offgrid_guard_011_cleared", "D1-02", "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC"),
        ("obs_d1_02_min_mismatches_to_patient_hlac_e", "D1-02", "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC"),
        ("obs_d1_02_wt_hlac_allele_retention_pct", "D1-02", "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC"),
        ("obs_d1_02_hlac_bystander_ct_pct", "D1-02", "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC"),
        ("obs_d1_02_chr6p21_interlocus_del_pct", "D1-02", "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC"),
        ("d1_auto_cipsc_source_arm", "D1-03", "HOLD_D1_03_INCOMPLETE_WGS_OR_CHRM_QC"),
        ("d1_03_head_to_head_adsc_vs_blood_comparator_attested", "D1-03", "HOLD_D1_03_HEAD_TO_HEAD_ADSC_VS_BLOOD_COMPARATOR_REQUIRED"),
        ("d1_03_adsc_vs_blood_bcor_tp53_ppm1d_wgs_clean", "D1-03", "HOLD_D1_03_INCOMPLETE_WGS_OR_CHRM_QC"),
        ("obs_d1_03_chrm_numt_amplicon_depth_x", "D1-03", "HOLD_D1_03_INCOMPLETE_WGS_OR_CHRM_QC"),
        ("obs_d1_03_chrm_novel_nonsyn_het_pct", "D1-03", "HOLD_D1_03_INCOMPLETE_WGS_OR_CHRM_QC"),
        ("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified", "D1-04", "HOLD_D1_04_MISSING_SWITCH_OR_GY018_INTEGER_KN"),
        ("d1_04_azt_kill_measured_after_forced_chr19_loss", "D1-04", "HOLD_D1_04_AZT_KILL_FORCED_CHR19_LOSS_CONDITION_REQUIRED"),
        ("obs_d1_04_g0_beta_htmpk_leak_pct", "D1-04", "HOLD_D1_04_MISSING_SWITCH_OR_GY018_INTEGER_KN"),
        ("obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct", "D1-04", "HOLD_D1_04_MISSING_SWITCH_OR_GY018_INTEGER_KN"),
        ("obs_d1_04_ap1903_survivor_count_k", "D1-04", "HOLD_D1_04_AP1903_GY018_INTEGER_KN_REQUIRED"),
        ("obs_d1_04_ap1903_effective_trials_n", "D1-04", "HOLD_D1_04_AP1903_GY018_INTEGER_KN_REQUIRED"),
        ("obs_d1_04_dec053_azt_dividing_log10_kill", "D1-04", "HOLD_D1_04_MISSING_SWITCH_OR_GY018_INTEGER_KN"),
        ("obs_d1_05_l3_day28_retention_pct", "D1-05", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"),
        ("obs_d1_05_l4_day28_retention_pct", "D1-05", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"),
        ("obs_d1_05_fsp1_ferroptosis_survival_pct", "D1-05", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"),
        ("obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct", "D1-06", "HOLD_D1_06_MISSING_KIR3DL1_NK_ASSAY"),
        ("obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct", "D1-07", "HOLD_D1_07_MISSING_VTCN1_SHEDDING_ASSAY"),
        ("d1_08_patient_carries_hlac0304", "D1-08", "HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED"),
        ("d1_08_hlac0304_positive_control_valid", "D1-08", "HOLD_D1_08_HLAC0304_POSITIVE_CONTROL_VALIDITY_REQUIRED"),
        ("obs_d1_08_hlac0304_lysed_events_k", "D1-08", "HOLD_D1_08_HLAC0304_GY018_INTEGER_KN_REQUIRED"),
        ("obs_d1_08_hlac0304_scored_targets_n", "D1-08", "HOLD_D1_08_HLAC0304_GY018_INTEGER_KN_REQUIRED"),
        ("obs_d1_08_hlac0304_tcell_lysis_u95_pct", "D1-08", "HOLD_D1_08_HLAC0304_REPORTED_U95_REQUIRED"),
        ("obs_d1_08_t1d_serum_c5b9_pos_events_k", "D1-08", "HOLD_D1_08_T1D_SERUM_C5B9_GY018_INTEGER_KN_REQUIRED"),
        ("obs_d1_08_t1d_serum_c5b9_scored_cells_n", "D1-08", "HOLD_D1_08_T1D_SERUM_C5B9_GY018_INTEGER_KN_REQUIRED"),
        ("obs_d1_08_t1d_serum_c5b9_u95_pct", "D1-08", "HOLD_D1_08_T1D_SERUM_C5B9_REPORTED_U95_REQUIRED"),
        ("obs_d1_09_autologous_turnaround_weeks", "D1-09", None),
    ]
    for _pdef_sw, _base_sw, _g0v_sw in ((p_l1a, base_l1a, None), (p_l1b, base_l1b, g0v_ok)):
        for _fk, _gid_sw, _exp_hold in _req_full_fields:
            _cell_mdd = _dt_cells_by_key[(_pdef_sw["product_variant_id"], _gid_sw)]["missing_data_disposition"]
            for _mode_sw in ("OMIT", "NULL", "UNAVAILABLE"):
                _obs_sw = {k: v for k, v in _base_sw.items() if k != _fk}
                if _mode_sw == "NULL":
                    _obs_sw[_fk] = None
                elif _mode_sw == "UNAVAILABLE":
                    _obs_sw[_fk] = "UNAVAILABLE"
                _out_sw = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
                    _obs_sw, product_def=_pdef_sw, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_sw,
                )
                assert _out_sw["reference_criteria_met"] is False, (_pdef_sw["product_variant_id"], _fk, _mode_sw, _out_sw)
                assert _out_sw["clinical_or_batch_release_authorized"] is False
                if _exp_hold is not None:
                    assert _exp_hold in _out_sw["reasons"], (_pdef_sw["product_variant_id"], _fk, _mode_sw, _exp_hold, _out_sw["reasons"])
                assert _cell_mdd in _out_sw["reasons"], (_pdef_sw["product_variant_id"], _gid_sw, _fk, _mode_sw, _cell_mdd, _out_sw["reasons"])
                if _fk == "d1_auto_hlac_sp3c4c_genotype_class":
                    assert _out_sw["d1_01_allele_derivation_mode"] == "MISSING_OR_UNAVAILABLE_GENOTYPE_CLASS_HOLD"
    # CLI parity on BOTH D1-02 single-field omissions AND FAIL cases (wgs_ok=False, min_mm=2) on L1A & L1B!
    for _pdef_d2, _base_d2, _g0v_d2 in ((p_l1a, base_l1a, None), (p_l1b, base_l1b, g0v_ok)):
        for _omit_k in ("d1_02_n1_phased_wgs_offgrid_guard_011_cleared", "obs_d1_02_min_mismatches_to_patient_hlac_e"):
            _obs_d2 = {k: v for k, v in _base_d2.items() if k != _omit_k}
            _py_d2 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_obs_d2, product_def=_pdef_d2, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_d2)
            _cli_d2_p = tmp_path / f"cli_d102_{_pdef_d2['product_variant_id']}_{_omit_k}.json"
            _env_d2: dict[str, object] = {"product_variant_id": _pdef_d2["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": _obs_d2}
            if _g0v_d2 is not None:
                _env_d2["global_gate0v_record"] = _g0v_d2
            _cli_d2_p.write_text(json.dumps(_env_d2), encoding="utf-8")
            _cli_d2 = json.loads(subprocess.run(
                [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_cli_d2_p)],
                cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
            ).stdout)
            assert _py_d2["reference_criteria_met"] is False and _cli_d2["reference_criteria_met"] is False
            assert _py_d2["status"] == _cli_d2["status"] == "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC", (_py_d2, _cli_d2)
            assert _py_d2["reasons"] == _cli_d2["reasons"]
            assert _py_d2["d1_02_offgrid_guard_evidence_mode"] == _cli_d2["d1_02_offgrid_guard_evidence_mode"] == "INCOMPLETE_D1_02_PHASED_WGS_OR_MIN_MISMATCH_HOLD"
        for _fail_tag, _fail_ov in (("wgs_false", {"d1_02_n1_phased_wgs_offgrid_guard_011_cleared": False}), ("min_mm_2", {"obs_d1_02_min_mismatches_to_patient_hlac_e": 2})):
            _obs_d2f = {**_base_d2, **_fail_ov}
            _py_d2f = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_obs_d2f, product_def=_pdef_d2, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_d2)
            _cli_d2fp = tmp_path / f"cli_d102f_{_pdef_d2['product_variant_id']}_{_fail_tag}.json"
            _env_d2f: dict[str, object] = {"product_variant_id": _pdef_d2["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": _obs_d2f}
            if _g0v_d2 is not None:
                _env_d2f["global_gate0v_record"] = _g0v_d2
            _cli_d2fp.write_text(json.dumps(_env_d2f), encoding="utf-8")
            _cli_d2f = json.loads(subprocess.run(
                [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_cli_d2fp)],
                cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
            ).stdout)
            assert _py_d2f["reference_criteria_met"] is False and _cli_d2f["reference_criteria_met"] is False
            assert _py_d2f["status"] == _cli_d2f["status"] == "FAIL_D1_02_N1_PHASED_WGS_OFFGRID_GUARD_011_HLAC_E_MISMATCH_CONFLICT", (_py_d2f, _cli_d2f)
            assert _py_d2f["reasons"] == _cli_d2f["reasons"]
            assert _py_d2f["d1_02_offgrid_guard_evidence_mode"] == _cli_d2f["d1_02_offgrid_guard_evidence_mode"] == "FAILED_D1_02_PHASED_WGS_OR_MIN_MISMATCH"
    # L1B Gate-0V every-field omission sweep
    for _gk in ("mcb_lot_id", "cvb_day14_log10_reduction", "cvb_day14_log10_reduction_l95", "cvb4_72h_early_titer_log10_pfu_u95", "day28_vp1_ie1_positive_k", "day28_vp1_ie1_cells_scored_n", "gate0v_verdict"):
        _g_om = {k: v for k, v in g0v_ok.items() if k != _gk}
        _out_gom = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=_g_om)
        assert _out_gom["reference_criteria_met"] is False, (_gk, _out_gom)
    # Component L2/L3/L4 every-field omission sweep
    p_l3 = by_id["ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9"]
    _comp_fixtures = [
        (p_l2, {"d1_08_patient_carries_hlac0304": True, "d1_08_hlac0304_positive_control_valid": True, "obs_d1_08_hlac0304_lysed_events_k": 0, "obs_d1_08_hlac0304_scored_targets_n": 59, "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 4.96}),
        (p_l3, {"obs_d1_04_ap1903_survivor_count_k": 0, "obs_d1_04_ap1903_effective_trials_n": 299_573_226, "obs_d1_05_l3_day28_retention_pct": 88.0, "obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct": 5.0}),
        (p_l4, {"d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified": True, "d1_04_azt_kill_measured_after_forced_chr19_loss": True, "obs_d1_04_g0_beta_htmpk_leak_pct": 0.05, "obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct": 0.05, "obs_d1_04_dec053_azt_dividing_log10_kill": 4.5, "obs_d1_05_l4_day28_retention_pct": 86.0, "obs_d1_05_fsp1_ferroptosis_survival_pct": 87.0}),
    ]
    for _pc_sw, _cfix in _comp_fixtures:
        assert vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_cfix, product_def=_pc_sw, evaluated_lot_id="LOT_C")["component_subgate_evaluation"]["component_subgates_passed"] is True
        for _ck in list(_cfix):
            _c_om = {k: v for k, v in _cfix.items() if k != _ck}
            _out_com = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_c_om, product_def=_pc_sw, evaluated_lot_id="LOT_C")
            assert _out_com["component_subgate_evaluation"]["component_subgates_passed"] is False, (_pc_sw["product_variant_id"], _ck, _out_com)

    # M. V15 Closure Guards:
    # 1. Symmetric order-invariant dual global_gate0v_record permutation matrix across Python AND --eval-json CLI!
    #    Covers int vs Decimal/float type mismatches (k=1 vs 1.0, k=0 vs 0.0, n=5990 vs 5990.0), lot mismatches + k>=1,
    #    dual independent failures, and explicit-null outer/inner conflicts on BOTH L1B and L1A.
    def _cli_eval_g0v(tag: str, env_json_str: str) -> dict[str, object]:
        _p = tmp_path / f"cli_g0v_{tag}.json"
        _p.write_text(env_json_str, encoding="utf-8")
        return json.loads(subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_p)],
            cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
        ).stdout)

    _g0v_k1_int = {**g0v_ok, "day28_vp1_ie1_positive_k": 1}
    _g0v_k1_dec = {**g0v_ok, "day28_vp1_ie1_positive_k": Decimal("1.0")}
    _g0v_k0_dec = {**g0v_ok, "day28_vp1_ie1_positive_k": Decimal("0.0")}
    _g0v_n_dec = {**g0v_ok, "day28_vp1_ie1_cells_scored_n": Decimal("5990.0")}
    _g0v_lot_k1 = {**g0v_ok, "mcb_lot_id": "LOT_OTHER", "day28_vp1_ie1_positive_k": 1}
    _g0v_d14_low = {**g0v_ok, "cvb_day14_log10_reduction": Decimal("2.0"), "cvb_day14_log10_reduction_l95": Decimal("1.5")}
    _g0v_t72_hi = {**g0v_ok, "cvb4_72h_early_titer_log10_pfu_u95": Decimal("2.5")}

    for _arm_tag, _pdef_g, _base_g, _hs_exp in (
        ("L1B", p_l1b, base_l1b, "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN"),
        ("L1A", p_l1a, base_l1a, "FAIL_HARD_STOP_D1_L1A_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"),
    ):
        for _ptag, _rA, _rB, _req_reasons in (
            ("k1_int_vs_dec", _g0v_k1_int, _g0v_k1_dec, (_hs_exp, "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS", "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE")),
            ("k0_int_vs_dec", g0v_ok, _g0v_k0_dec, ("FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS", "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE")),
            ("n_int_vs_dec", g0v_ok, _g0v_n_dec, ("FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS", "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE")),
            ("lot_k1_vs_ok", _g0v_lot_k1, g0v_ok, (_hs_exp, "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS")),
            ("d14_vs_t72", _g0v_d14_low, _g0v_t72_hi, ("FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",)),
        ):
            _fwd = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**_base_g, "global_gate0v_record": _rA}, product_def=_pdef_g, evaluated_lot_id="LOT_A", global_gate0v_record=_rB)
            _rev = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**_base_g, "global_gate0v_record": _rB}, product_def=_pdef_g, evaluated_lot_id="LOT_A", global_gate0v_record=_rA)
            assert _fwd == _rev, (_arm_tag, _ptag, _fwd, _rev)
            assert _fwd["reference_criteria_met"] is False and _fwd["clinical_or_batch_release_authorized"] is False
            for _rr_e in _req_reasons:
                assert _rr_e in _fwd["reasons"], (_arm_tag, _ptag, _rr_e, _fwd["reasons"])
            if _hs_exp in _req_reasons:
                assert _fwd["status"] == _hs_exp and _fwd["global_gate0v_standalone_subcheck_status"] == _hs_exp
            else:
                assert _fwd["global_gate0v_standalone_subcheck_status"] == "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS"
            # Verify raw JSON CLI matches direct Python byte-for-byte on BOTH orderings (preserving 1 vs 1.0 & 0 vs 0.0 lexemes)!
            def _to_json_obj(d: dict[str, object]) -> dict[str, object]:
                return {k: (float(v) if isinstance(v, Decimal) else v) for k, v in d.items()}
            for _ord_tag, _o_rec, _i_rec in (("fwd", _rB, _rA), ("rev", _rA, _rB)):
                _cli_out = _cli_eval_g0v(
                    f"{_arm_tag}_{_ptag}_{_ord_tag}",
                    json.dumps({
                        "product_variant_id": _pdef_g["product_variant_id"],
                        "evaluated_lot_id": "LOT_A",
                        "global_gate0v_record": _to_json_obj(_o_rec),
                        "observations": {**_base_g, "global_gate0v_record": _to_json_obj(_i_rec)},
                    }),
                )
                assert _cli_out == _fwd, (_arm_tag, _ptag, _ord_tag, _cli_out, _fwd)
        # Explicit-null outer vs valid inner AND valid outer vs explicit-null inner both conflict symmetrically!
        _null_in = _cli_eval_g0v(f"{_arm_tag}_null_in", json.dumps({"product_variant_id": _pdef_g["product_variant_id"], "evaluated_lot_id": "LOT_A", "global_gate0v_record": g0v_ok, "observations": {**_base_g, "global_gate0v_record": None}}))
        _null_out = _cli_eval_g0v(f"{_arm_tag}_null_out", json.dumps({"product_variant_id": _pdef_g["product_variant_id"], "evaluated_lot_id": "LOT_A", "global_gate0v_record": None, "observations": {**_base_g, "global_gate0v_record": g0v_ok}}))
        assert _null_in == _null_out and _null_in["reference_criteria_met"] is False and "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS" in _null_in["reasons"]

    # 2. Padded numeric strings & padded/non-ASCII-space sampling_design rejected across Python & CLI!
    for _pad_k, _pad_v in (
        ("obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct", " 5.0 "),
        ("obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct", "5.0\t"),
        ("obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct", " 4.0"),
        ("obs_d1_09_autologous_turnaround_weeks", "35.0 "),
    ):
        _out_padn = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, _pad_k: _pad_v}, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert _out_padn["reference_criteria_met"] is False and "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION" in _out_padn["reasons"], (_pad_k, _pad_v, _out_padn)
    _out_g0v_pad = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record={**g0v_ok, "cvb_day14_log10_reduction": " 3.5 "})
    assert _out_g0v_pad["status"] == "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE"
    for _bad_sd in (" independent ", "independent\u00a0replicates", "\tindependent"):
        _out_bsd = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": _bad_sd}, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert _out_bsd["reference_criteria_met"] is False and any("MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION" in r for r in _out_bsd["reasons"]), (_bad_sd, _out_bsd)

    # 3. Non-carrier D1-08 sampling contradiction parity across full-product (L1A & L1B) AND component L2!
    _base_l1a_nc = {**base_l1a, "d1_08_patient_carries_hlac0304": False, "d1_08_matched_allele_positive_control_valid": True, "d1_08_non_c0304_matched_allele_screen_passed": True}
    _base_l1b_nc = {**base_l1b, "d1_08_patient_carries_hlac0304": False, "d1_08_matched_allele_positive_control_valid": True, "d1_08_non_c0304_matched_allele_screen_passed": True}
    _base_l2_nc = {"d1_08_patient_carries_hlac0304": False, "d1_08_matched_allele_positive_control_valid": True, "d1_08_non_c0304_matched_allele_screen_passed": True}
    for _nc_pdef, _nc_obs, _nc_g0v in ((p_l1a, _base_l1a_nc, None), (p_l1b, _base_l1b_nc, g0v_ok), (p_l2, _base_l2_nc, None)):
        _nc_conf = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**_nc_obs, "sampling_design": "independent", "obs_d1_08_hlac0304_clustered_sampling": True}, product_def=_nc_pdef, evaluated_lot_id="LOT_A", global_gate0v_record=_nc_g0v)
        assert "HOLD_D1_08_HLAC0304_HOLD_CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS" in _nc_conf["reasons"], (_nc_pdef["product_variant_id"], _nc_conf)
        assert _nc_conf["reference_criteria_met"] is False
        _nc_cl = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**_nc_obs, "obs_d1_08_hlac0304_sampling_design": "not independent"}, product_def=_nc_pdef, evaluated_lot_id="LOT_A", global_gate0v_record=_nc_g0v)
        assert "HOLD_D1_08_HLAC0304_HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED" in _nc_cl["reasons"], (_nc_pdef["product_variant_id"], _nc_cl)
        assert _nc_cl["reference_criteria_met"] is False

    # N. V16 Closure Guards:
    # 1. Nested lot alias conflicts inside global_gate0v_record (mcb_lot_id vs lot_id/batch_id/mcb_batch_id/evaluated_lot_id)
    #    HOLD on zero events and preserve FAIL_HARD_STOP on k>=1 across Python & CLI!
    for _alias_k in ("lot_id", "batch_id", "mcb_batch_id", "evaluated_lot_id"):
        _g0v_al_conf = {**g0v_ok, _alias_k: "LOT_CONFLICT"}
        _out_al_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=_g0v_al_conf)
        assert _out_al_py["reference_criteria_met"] is False and "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in _out_al_py["reasons"], (_alias_k, _out_al_py)
        _out_al_cli = _cli_eval_g0v(f"nested_lot_{_alias_k}", json.dumps({"product_variant_id": p_l1b["product_variant_id"], "evaluated_lot_id": "LOT_A", "global_gate0v_record": _g0v_al_conf, "observations": base_l1b}))
        assert _out_al_cli == _out_al_py
        _out_al_k1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record={**_g0v_al_conf, "day28_vp1_ie1_positive_k": 1})
        assert _out_al_k1["status"].startswith("FAIL_HARD_STOP_") and "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in _out_al_k1["reasons"]

    # 2. Panel-level calibration (eta_recovery_l95=0.5, alpha_eta_recovery=0.05) propagates into nested Gate-0V Day-28!
    #    At n=5,990 with eta=0.5, U95 = 0.001 > 0.0005 -> MUST HOLD (not silently default Gate-0V to eta=1.0)!
    _panel_cal_obs = {
        **base_l1b,
        "eta_recovery_l95": Decimal("0.5"),
        "alpha_eta_recovery": Decimal("0.05"),
        "obs_d1_04_ap1903_effective_trials_n": 600_000_000,
        "obs_d1_08_hlac0304_scored_targets_n": 120,
        "obs_d1_08_hlac0304_tcell_lysis_u95_pct": Decimal("4.95"),
        "obs_d1_08_t1d_serum_c5b9_scored_cells_n": 120,
        "obs_d1_08_t1d_serum_c5b9_u95_pct": Decimal("4.95"),
    }
    _out_pcal_hold = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_panel_cal_obs, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record=g0v_ok)
    assert _out_pcal_hold["reference_criteria_met"] is False and "HOLD_D1_L1B_GATE0V_DAY28_INSUFFICIENT_DENOMINATOR" in _out_pcal_hold["reasons"], _out_pcal_hold
    _out_pcal_pass = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_panel_cal_obs, product_def=p_l1b, evaluated_lot_id="LOT_A", global_gate0v_record={**g0v_ok, "day28_vp1_ie1_cells_scored_n": 12_000})
    assert _out_pcal_pass["reference_criteria_met"] is True and _out_pcal_pass["rare_event_calibration_summary"]["l1b_gate0v_day28"]["effective_detection_factor_str"] == "0.5", _out_pcal_pass

    # 3. clustered_sampling: False alone reports EXPLICIT_NON_CLUSTERED_ONLY_INDEPENDENCE_ASSUMED (never overclaiming EXPLICIT_INDEPENDENT)!
    _out_nc_only = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "clustered_sampling": False}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert _out_nc_only["reference_criteria_met"] is True
    assert _out_nc_only["rare_event_calibration_summary"]["d1_04_ap1903"]["sampling_independence_provenance"] == "EXPLICIT_NON_CLUSTERED_ONLY_INDEPENDENCE_ASSUMED"

    # 4. Huge-denominator (n=10^100) 180-digit conservative reported U95 PASSES via Taylor -expm1 (no false FAIL),
    #    while strictly below true bound still FAILS, AND caller-enabled decimal.Inexact trap never breaks evaluation!
    import decimal as _dec_v16
    with _dec_v16.localcontext(_dec_v16.Context(prec=220)) as _c220:
        _x100 = (-_dec_v16.Decimal("0.05").ln()) / _dec_v16.Decimal(10**100)
        _true_u95_pct_100 = (_x100 - (_x100 * _x100) / 2) * 100
        _cons_180_pct = (_true_u95_pct_100 * _dec_v16.Decimal("1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"))
        _below_180_pct = (_true_u95_pct_100 * _dec_v16.Decimal("0.99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999"))
    _out_huge_ok = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_08_t1d_serum_c5b9_scored_cells_n": 10**100, "obs_d1_08_t1d_serum_c5b9_u95_pct": format(_cons_180_pct, "e")},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_huge_ok["reference_criteria_met"] is True, _out_huge_ok
    _out_huge_bad = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_08_t1d_serum_c5b9_scored_cells_n": 10**100, "obs_d1_08_t1d_serum_c5b9_u95_pct": format(_below_180_pct, "e")},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_huge_bad["status"] == "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND", _out_huge_bad
    with _dec_v16.localcontext() as _ctx_trap:
        _ctx_trap.traps[_dec_v16.Inexact] = True
        _ctx_trap.traps[_dec_v16.Rounded] = True
        _out_trapped = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert _out_trapped["reference_criteria_met"] is True, _out_trapped

    # 5. CommonMark 0.31.2 §2.2 space-plus-tab expansion (' \t', '  \t', '   \t' expand to col 4 indented code, NOT Setext heading!)
    for _cm_tab_in in (" \tfoo\n===\n", "  \tfoo\n===\n", "   \tfoo\n===\n", " \t# Not a heading\n"):
        assert vmod._parse_commonmark_headings(_cm_tab_in) == [], (_cm_tab_in, vmod._parse_commonmark_headings(_cm_tab_in))

    # O. Consolidated Evaluator, Math, CommonMark, Contract, & Provenance Closure Guards:
    # 1. Surrogate Unicode rejected in strings/IDs; mixed-type dict keys in extra_lot_bindings rejected cleanly (no TypeError)!
    _out_surr = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_\ud800")
    assert _out_surr["reference_criteria_met"] is False and "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" in _out_surr["reasons"]
    _out_mixed_lot = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A",
        extra_lot_bindings={1: "LOT_A", "mcb_lot_id": "LOT_A"},  # type: ignore[dict-item]
    )
    assert _out_mixed_lot["reference_criteria_met"] is False and "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION" in _out_mixed_lot["reasons"]

    # 2. Symmetric conflicting canonical product_def + product-neutral Gate-0V hard stop + unhidden nested Gate-0V lot errors on L1A/L1B!
    _out_pdef_conf_1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "product_def": p_l1b}, product_def=p_l1a, evaluated_lot_id="LOT_A",
        global_gate0v_record={**g0v_ok, "day28_vp1_ie1_positive_k": 1},
    )
    _out_pdef_conf_2 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "product_def": p_l1a}, product_def=p_l1b, evaluated_lot_id="LOT_A",
        global_gate0v_record={**g0v_ok, "day28_vp1_ie1_positive_k": 1},
    )
    assert _out_pdef_conf_1["status"] == _out_pdef_conf_2["status"] == "FAIL_HARD_STOP_D1_UNRESOLVED_PRODUCT_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"
    assert "FAIL_D1_CONFLICTING_PRODUCT_IDENTITY_DECLARATIONS" in _out_pdef_conf_1["reasons"]
    _out_l1a_g0v_bad_lot = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A",
        global_gate0v_record={**g0v_ok, "mcb_lot_id": "LOT_A", "lot_id": "LOT_OTHER"},
    )
    assert _out_l1a_g0v_bad_lot["reference_criteria_met"] is False and "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in _out_l1a_g0v_bad_lot["reasons"]

    # 3. Unlisted HLA-C allele blocks negative HLA-C*03:04 derivation (23-allele Lin 2023 Table 1 gate)!
    _out_unlisted_c0304 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            **{k: v for k, v in base_l1a.items() if k not in _strip_c0304},
            "patient_hla_c_alleles": ["HLA-C*04:01", "HLA-C*99:99"],
            "d1_08_matched_allele_positive_control_valid": True,
            "d1_08_non_c0304_matched_allele_screen_passed": True,
        },
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_unlisted_c0304["reference_criteria_met"] is False
    assert "HOLD_D1_08_UNLISTED_RAW_HLAC_ALLELE_BLOCKS_NEGATIVE_DERIVATION" in _out_unlisted_c0304["reasons"]

    # 4. Non-canonical diagnostic mode enforces full 45-cell contract & returns NON_CANONICAL_DIAGNOSTIC_CRITERIA_MET_UNPINNED!
    dt_pass_unpin = copy.deepcopy(dt_blk)
    for _c_up in dt_pass_unpin["product_by_gate_contract"]["cells"]:
        if _c_up["product_variant_id"] == p_l1a["product_variant_id"] and _c_up["gate_id"] == "D1-06":
            _c_up["threshold"] = "kir3dl1_missing_bw4_autologous_nk_lysis_pct<=9.0"
    _out_unpin_ok = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", decision_tree_block=dt_pass_unpin, require_canonical_sha=False,
    )
    assert _out_unpin_ok["reference_criteria_met"] is False
    assert _out_unpin_ok["diagnostic_criteria_met"] is True
    assert _out_unpin_ok["status"] == "NON_CANONICAL_DIAGNOSTIC_CRITERIA_MET_UNPINNED"
    dt_corrupt_44 = copy.deepcopy(dt_blk)
    dt_corrupt_44["product_by_gate_contract"]["cells"] = dt_corrupt_44["product_by_gate_contract"]["cells"][:-1]
    with pytest.raises(ValueError, match="exact 45-cell Cartesian contract"):
        vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
            base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", decision_tree_block=dt_corrupt_44, require_canonical_sha=False,
        )

    # 5. Component provenance fields, product-neutral component Gate-0V hard stop, & single-locus L1A/L1B scope disclosures!
    assert out_l1a["assembled_multilocus_manifest_bound"] is False
    assert out_l1a["reference_scenario_scope"] == "SINGLE_LOCUS_L1_DONOR_REFERENCE_SCENARIO_KEY_NOT_ASSEMBLED_MULTILOCUS_MANIFEST"
    assert isinstance(out_l1a["assumption_only_path_flags"], list) and len(out_l1a["assumption_only_path_flags"]) >= 1
    assert "coverage_regime" in out_l3 and "rare_event_calibration_summary" in out_l3 and "sampling_design_provenance_by_assay" in out_l3
    _out_l3_g0v_k1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {"obs_d1_04_ap1903_survivor_count_k": 0, "obs_d1_04_ap1903_raw_trials_n": 299_573_226, "obs_d1_05_l3_day28_retention_pct": 88.0, "obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct": 5.0},
        product_def=p_l3, evaluated_lot_id="LOT_L3", global_gate0v_record={**g0v_ok, "mcb_lot_id": "LOT_L3", "day28_vp1_ie1_positive_k": 1},
    )
    assert _out_l3_g0v_k1["status"] == "FAIL_HARD_STOP_OUT_OF_SCOPE_COMPONENT_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"

    # 6. Explicit-null & conflicting panel-vs-local calibration holds + bound-group split + ASCII sampling_design + non-contradictory dependent!
    _out_null_cal = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "eta_recovery_l95": None}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert _out_null_cal["reference_criteria_met"] is False and any("EXPLICIT_NULL_CALIBRATION_PARAMETER" in r for r in _out_null_cal["reasons"])
    _out_conf_cal = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "eta_recovery_l95": Decimal("0.5"), "alpha_eta_recovery": Decimal("0.05"), "obs_d1_04_ap1903_eta_recovery_l95": Decimal("0.8"), "obs_d1_04_ap1903_alpha_eta_recovery": Decimal("0.05")},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_conf_cal["reference_criteria_met"] is False and any("CONFLICTING_PANEL_AND_LOCAL_CALIBRATION" in r for r in _out_conf_cal["reasons"])
    _out_split_cal = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "eta_recovery_l95": Decimal("0.5"), "alpha_eta_recovery": Decimal("0.05"), "obs_d1_04_ap1903_eta_recovery_l95": Decimal("0.5")},
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_split_cal["reference_criteria_met"] is False and any("PARTIAL_LOCAL_CALIBRATION_SPLITS_BOUND_GROUP" in r for r in _out_split_cal["reasons"])
    _out_turk_sd = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "sampling_design": "\u0131ndependent"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert _out_turk_sd["reference_criteria_met"] is False and any("MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION" in r for r in _out_turk_sd["reasons"])
    _out_dep_nc = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "clustered_sampling": False, "independent_trials": False}, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert _out_dep_nc["reference_criteria_met"] is False and any("NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED" in r for r in _out_dep_nc["reasons"])
    assert not any("CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS" in r for r in _out_dep_nc["reasons"])

    # 7. Hostile Decimal context (prec=1, ROUND_DOWN, Inexact/Rounded traps), a_eta=1e1000000, reported_u95=MIN_ETINY / NaN, & 3-term Bonferroni!
    with _dec_v16.localcontext(_dec_v16.Context(prec=1, rounding=_dec_v16.ROUND_DOWN, traps=[_dec_v16.Inexact, _dec_v16.Rounded])):
        _out_prec1 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert _out_prec1["reference_criteria_met"] is True
        _out_prec1_d105 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035({**base_l1a, "obs_d1_05_l3_l4_day28_retention_pct": "86.0000002"}, product_def=p_l1a, evaluated_lot_id="LOT_A")
        assert _out_prec1_d105["status"] == "FAIL_D1_05_INCONSISTENT_REPORTED_MIN_L3_L4_RETENTION"
    _r_ovf = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=120, max_u95_fraction=Decimal("0.05"),
        eta_recovery_l95=Decimal("0.5"), alpha_eta_recovery=Decimal("1e1000000"),
    )
    assert _r_ovf["status"] == "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION"
    _r_etiny = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=59, max_u95_fraction=Decimal("0.05"),
        reported_u95_pct=Decimal("1e-999999"),
    )
    assert _r_etiny["reported_u95_meets_exact_bound"] is False
    _out_etiny = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {**base_l1a, "obs_d1_08_t1d_serum_c5b9_u95_pct": "1e-999999"}, product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_etiny["status"] == "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND"
    _r_nan = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_08", k=0, n_raw=59, max_u95_fraction=Decimal("0.05"),
        reported_u95_pct=Decimal("NaN"),
    )
    assert _r_nan["status"] == "HOLD_INVALID_REPORTED_U95_PCT"

    # 8. CommonMark 0.31.2 advanced block + inline regressions & required OPENISLET_CURRENT_LAYER_SCOPE Playbook footer!
    for _cm_src, _cm_exp_titles in (
        ("-\n  ```python\n  # not a heading\n  ```\n", []),
        ("1. item\n\n   # heading in list\n", ["heading in list"]),
        ("> foo\nbar\n---\n", []),
        ("<div\n  class=\"x\">\n# not heading\n</div>\n", []),
        ("# `4.1.20` **bold** *em*\n", ["4.1.20 bold em"]),
        ("# [4.1.20](/url \"title\")\n", ["4.1.20"]),
    ):
        assert [h["rendered_title"] for h in vmod._parse_commonmark_headings(_cm_src)] == _cm_exp_titles, _cm_src
    min_ft = _build_minimal_verifier_tree(tmp_path / "min_footer_del")
    pb_ft_p = min_ft / "docs" / "STAGE-0-EXECUTION-PLAYBOOK.md"
    pb_ft_txt = pb_ft_p.read_text(encoding="utf-8")
    pb_ft_p.write_text(pb_ft_txt.split("<!-- OPENISLET_CURRENT_LAYER_SCOPE -->")[0], encoding="utf-8")
    assert subprocess.run(
        [sys.executable, "-B", str(min_ft / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--check"],
        cwd=str(min_ft), env=env_ev, capture_output=True,
    ).returncode != 0

    # 9. Product contract binary64 round-trip precision check, kinship Decimal overflow immunity, & reserved prefix guard!
    pc_t = _load_mod("scripts/product_contract.py", "openislet_pc_v17_test")
    assert pc_t._is_valid_kinship_pi_hat(10**309) is False
    pc_repo = tmp_path / "pc_v17_repo"
    (pc_repo / "data" / "protocols").mkdir(parents=True)
    (pc_repo / "data" / "constructs").mkdir(parents=True)
    for _fa in (ROOT / "data" / "constructs").glob("*.fasta"):
        shutil.copy2(_fa, pc_repo / "data" / "constructs" / _fa.name)
    _orig_p_txt = (ROOT / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json").read_text(encoding="utf-8")
    _p_rt = pc_repo / "data" / "protocols" / "islet1_product_definitions_v0.8.0.json"
    _p_rt.write_text(_orig_p_txt.replace('"kinship_pi_hat_lt": 0.05', '"kinship_pi_hat_lt": 0.050000000000000001', 1), encoding="utf-8")
    with pytest.raises(ValueError, match="Precision-losing or underflowing JSON float"):
        pc_t.load_product_contracts(pc_repo, require_staged_fastas=True)
    _p_doc_fake = json.loads(_orig_p_txt)
    _fake_prod = copy.deepcopy(next(p for p in _p_doc_fake["products"] if p["product_variant_id"] == p_l1a["product_variant_id"]))
    _fake_prod["product_variant_id"] = "ISLET1_D1_FAKE_UNPINNED"
    _fake_prod.pop("staged_d1b_auto_v2_donor_binding", None)
    _p_doc_fake["products"].append(_fake_prod)
    _p_rt.write_text(json.dumps(_p_doc_fake, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Product identifier conflicts with its family"):
        pc_t.load_product_contracts(pc_repo, require_staged_fastas=True)
    _p_doc_al = json.loads(_orig_p_txt)
    _p_doc_al["products"][0]["aliases"] = list(_p_doc_al["products"][0].get("aliases", [])) + ["D1AUTO_V2_FAKE"]
    _p_rt.write_text(json.dumps(_p_doc_al, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="reserved D1b identifier/alias"):
        pc_t.load_product_contracts(pc_repo, require_staged_fastas=True)
    # L4:D1-04 confidence_method is deterministic component audit (not Clopper-Pearson)
    _l4_d104_cell = _dt_cells_by_key[(p_l4["product_variant_id"], "D1-04")]
    assert _l4_d104_cell["confidence_method"] == "deterministic_boolean_percent_and_log10_kill_component_assay_audit"

    # 10. Release manifest refuses --regenerate-repo-manifest on a sealed public bundle & refuses wrong marker bytes!
    grm_path_o = ROOT / "scripts" / "generate_release_manifest.py"
    if grm_path_o.is_file():
        sealed_chk = tmp_path / "sealed_chk"
        (sealed_chk / "scripts").mkdir(parents=True)
        shutil.copy2(grm_path_o, sealed_chk / "scripts" / "generate_release_manifest.py")
        (sealed_chk / ".openislet_public_bundle").write_bytes(b"openislet_public_bundle_v0.8.0\n")
        assert subprocess.run(
            [sys.executable, "-B", str(sealed_chk / "scripts" / "generate_release_manifest.py")],
            cwd=str(sealed_chk), env=env_ev, capture_output=True,
        ).returncode != 0
        assert subprocess.run(
            [sys.executable, "-B", str(sealed_chk / "scripts" / "generate_release_manifest.py"), "--regenerate-repo-manifest"],
            cwd=str(sealed_chk), env=env_ev, capture_output=True,
        ).returncode != 0

    # P. Table-driven conformance guards for diagnostic mode, component lot conflicts, unresolved-product Gate 0V, partial calibration & CommonMark Ex. 518:
    # 1. Missing staged FASTA under require_canonical_sha=False MUST NOT grant canonical REFERENCE_CRITERIA_MET or reference_criteria_met=True!
    no_fa_root = _build_minimal_verifier_tree(tmp_path / "min_no_d1_fasta")
    (no_fa_root / "data" / "constructs" / "islet1_d1b_auto_v2_v0.8.0.fasta").unlink()
    with pytest.raises(ValueError):
        vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", repo_root=no_fa_root, require_canonical_sha=True)
    out_no_fa_diag = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        base_l1a, product_def=p_l1a, evaluated_lot_id="LOT_A", repo_root=no_fa_root, require_canonical_sha=False,
    )
    assert out_no_fa_diag["controlling_staged_registry_verified"] is False, out_no_fa_diag
    assert out_no_fa_diag["evaluation_mode"] == "NON_CANONICAL_DIAGNOSTIC_UNPINNED", out_no_fa_diag
    assert out_no_fa_diag["reference_criteria_met"] is False, out_no_fa_diag
    assert out_no_fa_diag["diagnostic_criteria_met"] is True, out_no_fa_diag
    assert out_no_fa_diag["status"] == "NON_CANONICAL_DIAGNOSTIC_CRITERIA_MET_UNPINNED", out_no_fa_diag
    assert out_no_fa_diag["clinical_or_batch_release_authorized"] is False

    # 2. Component donor L2 under unpinned 20% diagnostic threshold MUST report component_subgates_passed=False, diagnostic_passed=True, & explicit diagnostic status!
    dt_l2_loose = copy.deepcopy(dt_blk)
    for _c in dt_l2_loose["product_by_gate_contract"]["cells"]:
        if _c["product_variant_id"] == p_l2["product_variant_id"] and _c["gate_id"] == "D1-08":
            _c["threshold"] = _c["threshold"].replace("U95<=5.0%", "U95<=20.0%")
    obs_l2_n15 = {
        "d1_08_patient_carries_hlac0304": True,
        "d1_08_hlac0304_positive_control_valid": True,
        "obs_d1_08_hlac0304_lysed_events_k": 0,
        "obs_d1_08_hlac0304_scored_targets_n": 15,
        "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 18.2,
    }
    out_l2_strict_n15 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(obs_l2_n15, product_def=p_l2, evaluated_lot_id="LOT_L2")
    assert out_l2_strict_n15["component_subgate_evaluation"]["component_subgates_passed"] is False
    out_l2_diag_n15 = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        obs_l2_n15, product_def=p_l2, evaluated_lot_id="LOT_L2", decision_tree_block=dt_l2_loose, require_canonical_sha=False,
    )
    assert out_l2_diag_n15["evaluation_mode"] == "NON_CANONICAL_DIAGNOSTIC_UNPINNED", out_l2_diag_n15
    assert out_l2_diag_n15["component_subgate_evaluation"]["component_subgates_passed"] is False, out_l2_diag_n15
    assert out_l2_diag_n15["component_subgate_evaluation"]["component_subgates_diagnostic_passed"] is True, out_l2_diag_n15
    assert out_l2_diag_n15["status"] == "NON_CANONICAL_DIAGNOSTIC_COMPONENT_SUBGATES_MET_UNPINNED_EXCLUDED_FROM_FULLSTACK", out_l2_diag_n15
    assert out_l2_diag_n15["empirical_clearance_status"] == "NON_CANONICAL_DIAGNOSTIC_COMPONENT_UNPINNED_NO_RELEASE", out_l2_diag_n15

    # 3. Diagnostic validation MUST reject operator='NOT_APPLICABLE' on active component cell L3:D1-07!
    dt_bad_l3_op = copy.deepcopy(dt_blk)
    for _c in dt_bad_l3_op["product_by_gate_contract"]["cells"]:
        if _c["product_variant_id"] == p_l3["product_variant_id"] and _c["gate_id"] == "D1-07":
            _c["operator"] = "NOT_APPLICABLE"
    with pytest.raises(ValueError, match="Component cell .* operator/disposition mismatch"):
        vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
            {"obs_d1_04_ap1903_survivor_count_k": 0, "obs_d1_04_ap1903_raw_trials_n": 299_573_226, "obs_d1_05_l3_day28_retention_pct": 88.0, "obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct": 5.0},
            product_def=p_l3, evaluated_lot_id="LOT_L3", decision_tree_block=dt_bad_l3_op, require_canonical_sha=False,
        )

    # 4. Component L2 with nested zero-event Gate-0V lot-alias conflict OR mismatched mcb_lot_id MUST HOLD in BOTH Python and --eval-json CLI!
    obs_l2_pass = {
        "d1_08_patient_carries_hlac0304": True,
        "d1_08_hlac0304_positive_control_valid": True,
        "obs_d1_08_hlac0304_lysed_events_k": 0,
        "obs_d1_08_hlac0304_scored_targets_n": 59,
        "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 4.96,
    }
    for _cname, _bad_g0v in (
        ("alias_conflict", {**g0v_ok, "mcb_lot_id": "LOT_A", "lot_id": "LOT_B"}),
        ("mcb_mismatch", {**g0v_ok, "mcb_lot_id": "LOT_B"}),
    ):
        _obs_c = {**obs_l2_pass, "global_gate0v_record": _bad_g0v}
        _out_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_obs_c, product_def=p_l2, evaluated_lot_id="LOT_A")
        assert _out_py["component_subgate_evaluation"]["component_subgates_passed"] is False, (_cname, _out_py)
        assert _out_py["status"] == "HOLD_D1_CONFLICTING_LOT_ID_BINDING", (_cname, _out_py)
        assert "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in _out_py["reasons"], (_cname, _out_py)
        assert "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE" in _out_py["reasons"], (_cname, _out_py)
        assert _out_py["global_gate0v_standalone_subcheck_status"] == "HOLD_D1_CONFLICTING_LOT_ID_BINDING", (_cname, _out_py)
        _cli_p = tmp_path / f"l2_g0v_{_cname}.json"
        _cli_p.write_text(json.dumps({"product_variant_id": p_l2["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": _obs_c}), encoding="utf-8")
        _out_cli = json.loads(subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_cli_p)],
            cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
        ).stdout)
        assert _out_cli["status"] == _out_py["status"] and _out_cli["reasons"] == _out_py["reasons"], (_out_cli, _out_py)
        assert _out_cli["component_subgate_evaluation"]["component_subgates_passed"] is False

    # 5. Positive Gate-0V event with NO canonical product ID and caller-supplied L1B arm MUST emit product-neutral hard stop (NEVER SWITCH_TO_L1A)!
    _obs_unres_g0v = {
        **base_l1b,
        "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR",
        "global_gate0v_record": {**g0v_ok, "day28_vp1_ie1_positive_k": 1},
    }
    _out_unres_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_obs_unres_g0v, evaluated_lot_id="LOT_A")
    assert _out_unres_py["status"] == "FAIL_HARD_STOP_D1_UNRESOLVED_PRODUCT_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT", _out_unres_py
    assert all("SWITCH_TO_L1A" not in r for r in _out_unres_py["reasons"]), _out_unres_py
    assert "HOLD_D1_CANONICAL_PRODUCT_VARIANT_ID_REQUIRED" in _out_unres_py["reasons"], _out_unres_py
    _cli_unres = tmp_path / "unres_l1b_arm_g0v_k1.json"
    _cli_unres.write_text(json.dumps({"evaluated_lot_id": "LOT_A", "observations": _obs_unres_g0v}), encoding="utf-8")
    _out_unres_cli = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_cli_unres)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert _out_unres_cli["status"] == _out_unres_py["status"] and _out_unres_cli["reasons"] == _out_unres_py["reasons"]

    # 6. Local f_eligible=0.9 with default eta=1.0 & default f_eligible_provenance MUST NOT claim CALIBRATED and MUST surface d1_04_ap1903 default assumptions!
    _out_pf = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
        {
            **base_l1a,
            "sampling_design": "independent",
            "obs_d1_04_ap1903_f_eligible": "0.9",
            "obs_d1_04_ap1903_raw_trials_n": 332_859_141,
            "obs_d1_04_ap1903_effective_trials_n": 332_859_141,
        },
        product_def=p_l1a, evaluated_lot_id="LOT_A",
    )
    assert _out_pf["reference_criteria_met"] is True, _out_pf
    _ap_meta = _out_pf["rare_event_calibration_summary"]["d1_04_ap1903"]
    assert "CALIBRATED" not in str(_ap_meta["coverage_regime"]), _ap_meta
    assert _ap_meta["coverage_regime"] == "OBSERVATION_LEVEL_GE_95.00PCT_PARTIAL_FACTOR_WITH_DEFAULT_ASSUMPTIONS", _ap_meta
    assert "d1_04_ap1903:DEFAULT_UNIT_RECOVERY_ASSUMED(eta=1.0)" in _out_pf["assumption_only_path_flags"], _out_pf["assumption_only_path_flags"]
    assert "d1_04_ap1903:DEFAULT_KNOWN_ELIGIBILITY_PROVENANCE_ASSUMED" in _out_pf["assumption_only_path_flags"], _out_pf["assumption_only_path_flags"]

    # 7. CommonMark 0.31.2 §6.5 nested link (Example 518) and image-in-link / link-in-image parity!
    assert [h["rendered_title"] for h in vmod._parse_commonmark_headings("# [foo [bar](/uri)](/uri)\n")] == ["[foo bar](/uri)"]
    assert [h["rendered_title"] for h in vmod._parse_commonmark_headings("# [![foo](/img)](/uri)\n")] == ["foo"]
    assert [h["rendered_title"] for h in vmod._parse_commonmark_headings("# ![foo [bar](/uri)](/img)\n")] == ["foo bar"]
    assert [h["rendered_title"] for h in vmod._parse_commonmark_headings("# [foo <https://example.org>](/uri)\n")] == ["[foo https://example.org](/uri)"]

    # Q. Table-driven conformance guards for allocation-only / missing-estimate calibration & symmetric panel-vs-nested Gate-0V failure under lot conflict:
    # 1. Missing f_eligible estimate with INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION + alpha_f_eligible MUST HOLD (never substitute 1.0 or award CALIBRATED)!
    _obs_miss_f = {
        **base_l1a,
        "obs_d1_04_ap1903_alpha_obs": "0.05",
        "obs_d1_04_ap1903_alpha_f_eligible": "0.05",
        "obs_d1_04_ap1903_f_eligible_provenance": "INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
    }
    _out_miss_f_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_obs_miss_f, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert _out_miss_f_py["reference_criteria_met"] is False, _out_miss_f_py
    assert "HOLD_D1_04_AP1903_HOLD_MISSING_F_ELIGIBLE_ESTIMATE_FOR_ALPHA_ALLOCATION" in _out_miss_f_py["reasons"], _out_miss_f_py
    assert "CALIBRATED" not in str(_out_miss_f_py["rare_event_calibration_summary"]["d1_04_ap1903"].get("coverage_regime")), _out_miss_f_py
    _cli_mf = tmp_path / "miss_f_elig.json"
    _cli_mf.write_text(json.dumps({"product_variant_id": p_l1a["product_variant_id"], "evaluated_lot_id": "LOT_A", "observations": _obs_miss_f}), encoding="utf-8")
    _out_miss_f_cli = json.loads(subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_cli_mf)],
        cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
    ).stdout)
    assert _out_miss_f_cli["status"] == _out_miss_f_py["status"] and _out_miss_f_cli["reasons"] == _out_miss_f_py["reasons"]

    # 2. Panel-level allocation-only alpha_eta_recovery=0.05 (no eta_recovery_l95) MUST HOLD and NEVER label the panel CALIBRATED!
    _obs_alloc_only = {**base_l1a, "alpha_eta_recovery": "0.05"}
    _out_alloc_py = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(_obs_alloc_only, product_def=p_l1a, evaluated_lot_id="LOT_A")
    assert _out_alloc_py["reference_criteria_met"] is False, _out_alloc_py
    assert "HOLD_D1_04_AP1903_HOLD_MISSING_ETA_RECOVERY_ESTIMATE_FOR_ALPHA_ALLOCATION" in _out_alloc_py["reasons"], _out_alloc_py
    assert "CALIBRATED" not in str(_out_alloc_py["coverage_regime"]), _out_alloc_py
    # Direct low-level CP calls with omitted or unit estimate under alpha allocation MUST ALSO fail closed!
    _cp_unit_eta = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=299_573_226, max_u95_fraction=Decimal("1e-8"),
        alpha_obs=Decimal("0.05"), alpha_eta_recovery=Decimal("0.05"),
    )
    assert _cp_unit_eta["status"] == "HOLD_MISSING_OR_UNIT_ETA_RECOVERY_FOR_ALPHA_ALLOCATION", _cp_unit_eta
    _cp_unit_f = vmod.evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="D1_04", k=0, n_raw=299_573_226, max_u95_fraction=Decimal("1e-8"),
        f_eligible_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
        alpha_obs=Decimal("0.05"), alpha_f_eligible=Decimal("0.05"),
    )
    assert _cp_unit_f["status"] == "HOLD_MISSING_OR_UNIT_F_ELIGIBLE_FOR_ALPHA_ALLOCATION", _cp_unit_f

    # 3. L1B Gate-0V explicit verdict=='FAIL' (and numeric failure) MUST preserve SWITCH_TO_L1A alongside HOLD_D1_CONFLICTING_LOT_ID_BINDING
    #    symmetrically whether the conflicting lot alias is nested in global_gate0v_record OR at panel observation level (Python + CLI)!
    for _g0v_fail_rec, _exp_fail_code in (
        ({**g0v_ok, "gate0v_verdict": "FAIL"}, "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN"),
        ({**g0v_ok, "cvb_day14_log10_reduction": "2.0", "cvb_day14_log10_reduction_l95": "1.8"}, "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN"),
    ):
        _out_nested = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
            base_l1b, product_def=p_l1b, evaluated_lot_id="LOT_A",
            global_gate0v_record={**_g0v_fail_rec, "lot_id": "LOT_B"},
        )
        _out_panel = vmod.evaluate_stage0_d1b_auto_v2_gates_iip0035(
            {**base_l1b, "lot_id": "LOT_B"}, product_def=p_l1b, evaluated_lot_id="LOT_A",
            global_gate0v_record=_g0v_fail_rec,
        )
        for _res_lot in (_out_nested, _out_panel):
            assert _res_lot["status"] == _exp_fail_code, _res_lot
            assert _exp_fail_code in _res_lot["reasons"] and "HOLD_D1_CONFLICTING_LOT_ID_BINDING" in _res_lot["reasons"], _res_lot
            assert _res_lot["clinical_or_batch_release_authorized"] is False
        assert _out_nested["status"] == _out_panel["status"] and _out_nested["reasons"] == _out_panel["reasons"]
        _cli_pl = tmp_path / "l1b_panel_lot_fail.json"
        _cli_pl.write_text(json.dumps({
            "product_variant_id": p_l1b["product_variant_id"], "evaluated_lot_id": "LOT_A",
            "observations": {**base_l1b, "lot_id": "LOT_B"}, "global_gate0v_record": _g0v_fail_rec,
        }), encoding="utf-8")
        _out_pl_cli = json.loads(subprocess.run(
            [sys.executable, "-B", str(ROOT / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--eval-json", str(_cli_pl)],
            cwd=str(ROOT), env=env_ev, capture_output=True, text=True, check=True,
        ).stdout)
        assert _out_pl_cli["status"] == _out_panel["status"] and _out_pl_cli["reasons"] == _out_panel["reasons"]


def test_06_ring1_604_studies_celltype_grounding_and_frozen_studies() -> None:
    ext_files = sorted((ROOT / "data" / "extracted").glob("[0-9]*.json"))
    assert len(ext_files) == 604
    docs = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in ext_files}
    new10 = ("2642862", "31427818", "35953586", "36944335", "38843830", "39753706", "41436438", "41690310", "41764741", "41999746")
    for pm in new10:
        assert pm in docs, f"Missing new Ring-1 study {pm}"
    assert "41201433" not in docs, "Jetley 41201433 stays in CEG M0051 + dossier cache only so Layer-1 Pareto remains 0-byte identical"
    assert sum(len(d.get("arms", [])) for d in docs.values()) == 708
    assert sum(1 for d in docs.values() for a in d.get("arms", []) if a.get("is_reference_control") is True) == 65
    assert sum(len(a.get("gate2_3_4_genetic_and_epitope_edits", [])) for d in docs.values() for a in d.get("arms", [])) == 261
    assert sum(len(d.get("atomic_effects", [])) for d in docs.values()) == 3174

    # Grounding & structured challenge/edit checks (Amdare non-islet T-cell line peptide assay,
    # Deuse 2-arm mouse autologous neoantigen in-vivo vs human in-vitro split, Uenaka non-Stage0_PSC + unredacted curator description,
    # Rouhani 2-arm 696-sample survey vs BCOR-mutant neural subset, Lin autologous vs allogeneic NHP,
    # Sutherland adult whole-pancreas clinical series, all 15 new arms in_dish_vascular_assembly=None):
    for pm in new10:
        for arm in docs[pm]["arms"]:
            assert arm["gate5_engraftment_niche"]["in_dish_vascular_assembly"] is None, (pm, arm["arm_id"])
    assert len(docs["31427818"]["arms"]) == 2 and len(docs["31427818"]["atomic_effects"]) == 2
    assert docs["31427818"]["arms"][0]["challenge_conditions"][0]["challenge_type"] == "InVivo_Autologous_Neoantigen_Graft"
    assert docs["31427818"]["arms"][1]["gate1_substrate_and_differentiation"]["cell_source_class"] == "hPSC_Derived_Unspecified_Line"
    assert docs["31427818"]["arms"][1]["challenge_conditions"][0]["challenge_type"] == "Autologous_T_Cell_Neoantigen_InVitro_Assay"
    amd_arm = docs["41764741"]["arms"][0]
    assert amd_arm["gate1_substrate_and_differentiation"]["cell_source_class"] == "NonIslet_Immune_Endothelial_or_Hematopoietic_Model"
    assert amd_arm["gate1_substrate_and_differentiation"]["canonical_differentiation_stage_evaluated"] is None
    amd_ch = amd_arm["challenge_conditions"][0]
    assert amd_ch["challenge_type"] == "Peptide_Restimulation_Islet_Derived_T_Cell_Line_Assay"
    assert "HLA-C*03:04" in str(amd_ch["specific_autoantigen_or_hla_restriction"])
    assert docs["41690310"]["arms"][0]["gate1_substrate_and_differentiation"]["canonical_differentiation_stage_evaluated"] is None
    assert docs["41690310"]["arms"][0]["gate2_3_4_genetic_and_epitope_edits"][0]["gene_symbol_hgnc"] == "MULTI_LOCUS_SAFE_HARBOR"
    assert "REDACTED" not in str(docs["41690310"]["arms"][0]["gate2_3_4_genetic_and_epitope_edits"][0]["sub_molecular_coordinate"]["exact_variant_or_swap_description"])
    assert docs["2642862"]["arms"][0]["gate1_substrate_and_differentiation"]["canonical_differentiation_stage_evaluated"] is None
    assert docs["2642862"]["arms"][0]["challenge_conditions"][0]["challenge_type"] == "InVivo_Allogeneic_or_Autoimmune_Graft"
    assert docs["38843830"]["study_metadata"]["max_biological_fidelity_tier"] == "Tier_2C_NHP_Allogeneic_or_Autologous_NonIslet"
    assert all(a["gate1_substrate_and_differentiation"]["cell_source_class"] == "NHP_iPSC_Cardiomyocyte_NonIslet" for a in docs["38843830"]["arms"])
    assert docs["38843830"]["arms"][0]["challenge_conditions"][0]["challenge_type"] == "InVivo_Autologous_Graft_No_IS"
    assert docs["38843830"]["arms"][1]["challenge_conditions"][0]["challenge_type"] == "InVivo_Allogeneic_or_Autoimmune_Graft"
    assert len(docs["35953586"]["arms"]) == 2 and docs["35953586"]["arms"][0]["sample_size_n"] == 696 and docs["35953586"]["arms"][1]["sample_size_n"] is None
    assert len(docs["35953586"]["atomic_effects"]) == 2
    assert docs["35953586"]["atomic_effects"][0]["focal_arm_id"] == docs["35953586"]["arms"][0]["arm_id"]
    assert docs["35953586"]["atomic_effects"][1]["focal_arm_id"] == docs["35953586"]["arms"][1]["arm_id"]
    assert docs["35953586"]["atomic_effects"][0]["secretome_or_differentiation_penalty_observed"] is False
    assert docs["35953586"]["atomic_effects"][1]["secretome_or_differentiation_penalty_observed"] is True
    assert len(docs["39753706"]["arms"]) == 2 and docs["39753706"]["arms"][1]["is_reference_control"] is True
    assert len(docs["39753706"]["atomic_effects"]) == 3
    assert docs["39753706"]["atomic_effects"][0]["comparator_arm_id"] == docs["39753706"]["arms"][1]["arm_id"]
    assert docs["39753706"]["atomic_effects"][2]["comparator_arm_id"] is None
    assert len(docs["41436438"]["atomic_effects"]) == 2
    assert len(docs["41999746"]["arms"]) == 2 and docs["41999746"]["arms"][1]["is_reference_control"] is True
    assert "REDACTED" not in str(docs["2642862"]["arms"][0]["arm_label"])
    # Verify ALL 10 new Ring-1 records respect the <= 350-char local excerpt limit when unredacted
    for pm in new10:
        for eff in docs[pm]["atomic_effects"]:
            q = str(eff.get("verbatim_evidence_quote") or "")
            if not q.startswith("[REDACTED_"):
                assert len(q) <= 350, f"PMID {pm} {eff['effect_id']} quote len={len(q)} > 350"

    # Frozen studies 17375075 & 39571901 untouched (strictly in full private repo checkouts, NOT dossier-augmented public bundles!)
    if (ROOT / "data" / "dossiers").is_dir() and not (ROOT / ".openislet_public_bundle").is_file():
        assert hashlib.sha256((ROOT / "data" / "extracted" / "17375075.json").read_bytes()).hexdigest() == "b8decbb5b7ef6807b9ed53eddee4159c98d524200dbbfb8601d0dcf02f5ea3b5"
        assert hashlib.sha256((ROOT / "data" / "extracted" / "39571901.json").read_bytes()).hexdigest() == "73d155458b30d631326a5362403e8fa30d03a92a93dbcf44da960e686a3998f1"


def test_07_ceg_1182_typed_arms_crosslink_provenance_and_sourceguard(tmp_path: Path) -> None:
    ceg = json.loads((ROOT / "data" / "extracted" / "context_evidence_graph_v1.json").read_text(encoding="utf-8"))
    bceg = _load_mod("scripts/build_context_evidence_graph.py", "openislet_bceg_iip0035_test")
    assert len(bceg.CURATED_INELIGIBLE_EDGE_IDS) == 66
    assert len(ceg["studies"]) == 1182
    assert len(ceg["arms"]) == 7749
    assert len(ceg["outcomes"]) == 14175
    assert len(ceg["edges"]) == 2494
    elig = sum(1 for e in ceg["edges"] if e.get("primary_certification_eligible") is True)
    inelig = sum(1 for e in ceg["edges"] if e.get("primary_certification_eligible") is False)
    assert (elig, inelig) == (2273, 221)

    # All 8 new IIP-0035 outcomes carry typed disjoint intervention_arm_ids, comparator_arm_ids, focal_arm_id & singular intervention_arm_id
    for oid, exp_iv, exp_cmp in (
        ("M0189:O1", ["VAR:M0189:A1"], ["VAR:M0189:A2_allogeneic_ctrl"]),
        ("M0189:O2", ["VAR:M0189:A2_allogeneic_ctrl"], []),
        ("M1407:O1", ["VAR:M1407:A1"], []),
        ("PMID:36944335:O1", ["VAR:PMID:36944335:A1"], []),
        ("PMID:39753706:O1", ["VAR:PMID:39753706:A1"], ["VAR:PMID:39753706:A2_original_protocol_ctrl"]),
        ("PMID:39753706:O2", ["VAR:PMID:39753706:A1"], []),
        ("PMID:41436438:O1", ["VAR:PMID:41436438:A1"], []),
        ("PMID:41999746:O1", ["VAR:PMID:41999746:A1"], ["VAR:PMID:41999746:A2_p53_suppressed"]),
    ):
        o = ceg["outcomes"][oid]
        assert o.get("intervention_arm_ids") == exp_iv
        assert o.get("comparator_arm_ids") == exp_cmp
        assert o.get("focal_arm_id") == exp_iv[0]
        assert o.get("intervention_arm_id") == exp_iv[0]

    # Reconciled stat_type & tier on reused outcomes (including M0051:O2/O3/O4/O6, M0241:O2/O3, M0206:O10/O13/O14/O16, M0159:O1/O2/O3!)
    for oid_np in ("M0051:O1", "M0051:O2", "M0051:O3", "M0051:O4", "M0051:O5", "M0051:O6", "M0154:O1", "M0154:O2", "M0154:O3", "M0154:O4", "M0241:O1", "M0241:O2", "M0241:O3", "M0241:O4", "M0206:O16", "M0206:O17"):
        assert ceg["outcomes"][oid_np]["stat_type"] == "qualitative_no_pvalue"
    assert ceg["outcomes"]["M0206:O10"]["stat_type"] == "descriptive_scalar"
    assert ceg["outcomes"]["M0206:O10"]["intervention_value"] == "BCOR mutations present in 11 F-hiPSCs"
    assert ceg["outcomes"]["M0206:O10"]["comparator_value"] == "0 reads with BCOR mutation in founding fibroblasts"
    assert ceg["outcomes"]["M0206:O13"]["stat_type"] == "descriptive_scalar"
    assert ceg["outcomes"]["M0206:O14"]["stat_type"] == "descriptive_scalar"
    for oid_159 in ("M0159:O1", "M0159:O2", "M0159:O3"):
        assert ceg["outcomes"][oid_159]["manuscript_tier"] == "T5_Rodent_or_CellLine"
    edges_by_id = {e["edge_id"]: e for e in ceg["edges"]}
    for leg_inelig in ("M0051:L1", "M0051:L2", "M0051:L3", "M0154:L1", "M0154:L2", "M0159:L1", "M0159:L2", "M0206:L2", "M0241:L1", "M0241:L2"):
        assert edges_by_id[leg_inelig]["primary_certification_eligible"] is False
        assert "legacy_certification_adjudication_note" in edges_by_id[leg_inelig]
    assert edges_by_id["M0206:L1"]["primary_certification_eligible"] is True
    assert "legacy_certification_adjudication_note" in edges_by_id["M0206:L1"]
    assert "human in vivo graft rejection" not in edges_by_id["M0241:L3"]["claim_summary"].lower()
    assert "EFF:Microvascular_Inosculation_Angiogenesis_IBMIR" not in (ceg["studies"]["M0189"].get("candidate_effector_topics") or [])
    assert "EFF:Microvascular_Inosculation_Angiogenesis_IBMIR" not in (ceg["outcomes"]["M0189:O1"].get("candidate_effector_topics") or [])
    assert "EFF:Microvascular_Inosculation_Angiogenesis_IBMIR" not in (edges_by_id["M0189:L1"].get("candidate_effector_topics") or [])
    assert len(ceg["screening_funnel"]["iip0035_record_level_screening_ledger"]) == 11
    assert "M0206:O17" in edges_by_id["M0206:L3"]["cited_outcome_tuples"]
    assert edges_by_id["M0159:L3"]["cited_outcome_tuples"] == ["M0159:O2"]
    assert edges_by_id["M0241:L3"]["cited_outcome_tuples"] == ["M0241:O1", "M0241:O2", "M0241:O4"]
    # Atlas includes dedicated D1b design bucket (not lumped into A/B/X)
    atlas = json.loads((ROOT / "data" / "extracted" / "layer1_context_evidence_graph_atlas.json").read_text(encoding="utf-8"))
    assert "D1b" in atlas["gate_and_design_claim_balance"]["by_design"]
    assert atlas["gate_and_design_claim_balance"]["by_design"]["D1b"]["total_links"] >= 11
    for cl_eid in ("M0051:L4", "M0154:L3", "M0159:L3", "M0206:L3", "M0241:L3"):
        assert edges_by_id[cl_eid]["primary_certification_eligible"] is False
        assert isinstance(edges_by_id[cl_eid].get("provenance_note"), str) and len(edges_by_id[cl_eid]["provenance_note"]) >= 40
    # Sutherland M1407:L1 & Amdare M0159:L3 map to adaptive Gate2B (not innate Gate3)
    assert "Gate2B" in edges_by_id["M1407:L1"]["target_gate_or_module"]
    assert "Gate2B" in edges_by_id["M0159:L3"]["target_gate_or_module"]

    det = bceg.verify_ceg_payload_detailed(ceg, repo_root=ROOT)
    has_dossiers = (ROOT / "data" / "dossiers").is_dir()
    assert det["reused_legacy_crosslink_outcomes_exact_span_verified"] is False
    if has_dossiers:
        assert det["verification_mode"] == "iip0035_11dossiers_8new_exact_spans_verified"
        assert det["iip0035_sourceguard_executed"] is True
        assert len(det["iip0035_sourceguard_scanned_study_ids"]) == 6
        assert len(det["iip0035_sourceguard_scanned_edge_ids"]) == 11
        assert det["iip0035_sourceguard_scan_scope"] == "executed_on_6_new_study_families_11_new_edges_27_touched_legacy_outcomes_and_11_adjudicated_legacy_edges"
        assert len(det["iip0035_sourceguard_scanned_touched_legacy_outcome_ids"]) == 27
        assert len(det["iip0035_sourceguard_scanned_adjudicated_legacy_edge_ids"]) == 11
        assert len(det["iip0035_dossiers_verified"]) == 11
        assert det["iip0035_verified_exact_span_outcome_count"] == 8
        assert det["iip0035_verified_exact_span_outcome_ids"] == [
            "M0189:O1", "M0189:O2", "M1407:O1", "PMID:36944335:O1",
            "PMID:39753706:O1", "PMID:39753706:O2", "PMID:41436438:O1", "PMID:41999746:O1",
        ]
    else:
        assert det["verification_mode"] == "structural_public_only"
        assert det["iip0035_sourceguard_executed"] is False
        assert det["iip0035_sourceguard_scanned_study_ids"] == []
        assert det["iip0035_sourceguard_scanned_edge_ids"] == []
        assert det["iip0035_sourceguard_scan_scope"] == "skipped_public_bundle_no_local_dossiers"
        assert det["iip0035_verified_exact_span_outcome_count"] == 0
        with pytest.raises(ValueError, match="Incomplete local dossier set"):
            bceg.verify_ceg_payload(ceg, repo_root=ROOT, require_sources=True)
    # Verify --verify CLI stdout explicitly distinguishes verification_mode & exact_spans_verified
    env_bceg = dict(os.environ)
    env_bceg["PYTHONDONTWRITEBYTECODE"] = "1"
    bceg_cli_out = subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts" / "build_context_evidence_graph.py"), "--verify"],
        cwd=str(ROOT), env=env_bceg, capture_output=True, text=True, check=True,
    ).stdout
    assert f"mode={det['verification_mode']}" in bceg_cli_out, bceg_cli_out
    assert f"exact_spans_verified={det['iip0035_verified_exact_span_outcome_count']}" in bceg_cli_out, bceg_cli_out

    # Tamper probe: mutating a cross-link edge's outcome set or deleting a required edge MUST fail closed!
    bad_edge_bind = copy.deepcopy(ceg)
    for e in bad_edge_bind["edges"]:
        if e["edge_id"] == "M0159:L3":
            e["cited_outcome_tuples"] = ["M0159:O1"]
    with pytest.raises(ValueError, match="cited_outcome_tuples"):
        bceg.verify_ceg_payload(bad_edge_bind, repo_root=ROOT)

    bad_leg_elig = copy.deepcopy(ceg)
    for e in bad_leg_elig["edges"]:
        if e["edge_id"] == "M0154:L2":
            e["primary_certification_eligible"] = True
    with pytest.raises(ValueError, match="Curated-ineligible edge M0154:L2|Adjudicated legacy qualitative-only edge M0154:L2"):
        bceg.verify_ceg_payload(bad_leg_elig, repo_root=ROOT)

    # Removing required singular intervention_arm_id MUST fail closed!
    bad_iv_missing = copy.deepcopy(ceg)
    bad_iv_missing["outcomes"]["M0189:O1"].pop("intervention_arm_id", None)
    with pytest.raises(ValueError, match="intervention_arm_id"):
        bceg.verify_ceg_payload(bad_iv_missing, repo_root=ROOT)

    # Tamper probes: arm overlap, focal-in-comparator, missing cross-link provenance_note, mismatched owner PMID!
    bad_arm = copy.deepcopy(ceg)
    bad_arm["outcomes"]["M0189:O1"]["comparator_arm_ids"] = ["VAR:M0189:A1"]
    with pytest.raises(ValueError, match="overlap"):
        bceg.verify_ceg_payload(bad_arm, repo_root=ROOT)

    bad_focal = copy.deepcopy(ceg)
    bad_focal["outcomes"]["M0189:O1"]["intervention_arm_ids"] = ["VAR:M0189:A2_allogeneic_ctrl"]
    bad_focal["outcomes"]["M0189:O1"]["intervention_arm_id"] = "VAR:M0189:A2_allogeneic_ctrl"
    bad_focal["outcomes"]["M0189:O1"]["comparator_arm_ids"] = ["VAR:M0189:A1"]
    bad_focal["outcomes"]["M0189:O1"]["comparison_arms"] = "A2_allogeneic_ctrl vs A1"
    with pytest.raises(ValueError, match="focal_arm_id"):
        bceg.verify_ceg_payload(bad_focal, repo_root=ROOT)

    bad_pnote = copy.deepcopy(ceg)
    for e in bad_pnote["edges"]:
        if e["edge_id"] == "M0159:L3":
            e["provenance_note"] = "This filler note is longer than forty characters without disclosing provenance."
    with pytest.raises(ValueError, match="provenance_note"):
        bceg.verify_ceg_payload(bad_pnote, repo_root=ROOT)

    bad_iv_sing = copy.deepcopy(ceg)
    bad_iv_sing["outcomes"]["M0189:O1"]["intervention_arm_id"] = "VAR:M0189:A2_allogeneic_ctrl"
    with pytest.raises(ValueError, match="intervention_arm_id"):
        bceg.verify_ceg_payload(bad_iv_sing, repo_root=ROOT)

    bad_pmid = copy.deepcopy(ceg)
    bad_pmid["studies"]["PMID:41436438"]["pmid"] = "99999999"
    bad_pmid["outcomes"]["PMID:41436438:O1"]["quote_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="owner PMID"):
        bceg.verify_ceg_payload(bad_pmid, repo_root=ROOT)

    # V8 CEG semantic & record-level ledger tamper probes:
    # 1. Duplicate queue_id / unrecognized transition / fractional 0.5 delta in 11-row screening ledger MUST fail closed!
    bad_led_dup = copy.deepcopy(ceg)
    bad_led_dup["screening_funnel"]["iip0035_record_level_screening_ledger"][1]["queue_id"] = bad_led_dup["screening_funnel"]["iip0035_record_level_screening_ledger"][0]["queue_id"]
    with pytest.raises(ValueError, match="unexpected or duplicate queue_id"):
        bceg.verify_ceg_payload(bad_led_dup, repo_root=ROOT)
    bad_led_flt = copy.deepcopy(ceg)
    bad_led_flt["screening_funnel"]["iip0035_record_level_screening_ledger"][0]["d_screened_total"] = 0.5
    with pytest.raises(ValueError, match="strict int"):
        bceg.verify_ceg_payload(bad_led_flt, repo_root=ROOT)
    # 2. Promoting adjudicated ineligible M0206:L2 back to eligible OR stripping M0206:L1 quantitative outcomes OR downgrading O9/O12/O14 to qualitative_statement_only MUST fail closed!
    bad_m206_l2_elig = copy.deepcopy(ceg)
    for e in bad_m206_l2_elig["edges"]:
        if e["edge_id"] == "M0206:L2":
            e["primary_certification_eligible"] = True
    with pytest.raises(ValueError, match="M0206:L2"):
        bceg.verify_ceg_payload(bad_m206_l2_elig, repo_root=ROOT)
    bad_m206 = copy.deepcopy(ceg)
    for e in bad_m206["edges"]:
        if e["edge_id"] == "M0206:L1":
            e["cited_outcome_tuples"] = ["M0206:O16"]
            e["cited_outcome_ids"] = ["M0206:O16"]
    for _oid_m1 in ("M0206:O1", "M0206:O4", "M0206:O5", "M0206:O9", "M0206:O12", "M0206:O14"):
        bad_m206["outcomes"][_oid_m1]["cited_by_edge_ids"] = [
            x for x in (bad_m206["outcomes"][_oid_m1].get("cited_by_edge_ids") or []) if x != "M0206:L1"
        ]
    with pytest.raises(ValueError, match="M0206:L1 lost required quantitative"):
        bceg.verify_ceg_payload(bad_m206, repo_root=ROOT)
    bad_m206_qs = copy.deepcopy(ceg)
    for _oid_m2 in ("M0206:O9", "M0206:O10", "M0206:O12", "M0206:O14"):
        bad_m206_qs["outcomes"][_oid_m2]["stat_type"] = "qualitative_statement_only"
        for _fk in ("endpoint_name", "intervention_value"):
            if isinstance(bad_m206_qs["outcomes"][_oid_m2].get(_fk), str):
                bad_m206_qs["outcomes"][_oid_m2][_fk] = re.sub(
                    r"stat=[^|\]]+", "stat=qualitative_statement_only", bad_m206_qs["outcomes"][_oid_m2][_fk]
                )
        bad_m206_qs["outcomes"][_oid_m2]["quantitative_summary"] = "Qualitative statement without numbers."
    with pytest.raises(ValueError, match="M0206:L1 lost required quantitative"):
        bceg.verify_ceg_payload(bad_m206_qs, repo_root=ROOT)
    # V10 CEG Red-Team Guards: citation-year-only qualitative summary, forged zero-hash stub, M0189:O1 manuscript_tier tamper, M0189:L1 claim_summary tamper, and missing discovery_provenance MUST fail closed!
    bad_m206_yr = copy.deepcopy(ceg)
    bad_m206_yr["outcomes"]["M0206:O9"]["quantitative_summary"] = "Rouhani 2022 qualitative statement without measurements."
    with pytest.raises(ValueError, match="M0206:L[12]"):
        bceg.verify_ceg_payload(bad_m206_yr, repo_root=ROOT)
    bad_m206_stub = copy.deepcopy(ceg)
    bad_m206_stub["outcomes"]["M0206:O9"]["quantitative_summary"] = "[Structured outcome | stat=mean_sd | direction=unfavorable | text_sha256=" + ("0" * 64) + "]"
    with pytest.raises(ValueError, match="M0206:L[12]"):
        bceg.verify_ceg_payload(bad_m206_stub, repo_root=ROOT)
    bad_m189_tier = copy.deepcopy(ceg)
    bad_m189_tier["outcomes"]["M0189:O1"]["manuscript_tier"] = "T1_Human_Clinical"
    with pytest.raises(ValueError, match="full-object.*digest"):
        bceg.verify_ceg_payload(bad_m189_tier, repo_root=ROOT)
    bad_m189_edge = copy.deepcopy(ceg)
    for _ee in bad_m189_edge["edges"]:
        if _ee["edge_id"] == "M0189:L1":
            _ee["claim_summary"] = "Tampered claim summary."
    with pytest.raises(ValueError, match="full-object digest"):
        bceg.verify_ceg_payload(bad_m189_edge, repo_root=ROOT)
    bad_disc_prov = copy.deepcopy(ceg)
    for _lr in bad_disc_prov["screening_funnel"]["iip0035_record_level_screening_ledger"]:
        if _lr["queue_id"] == "PMID:39753706":
            _lr.pop("discovery_provenance", None)
    with pytest.raises(ValueError, match="discovery_provenance"):
        bceg.verify_ceg_payload(bad_disc_prov, repo_root=ROOT)
    # 3. Consistently reversing Wang PMID:39753706:O1 treatment vs control arms, flipping direction, OR mutating O2 10->2 days MUST fail closed!
    bad_wang = copy.deepcopy(ceg)
    wo = bad_wang["outcomes"]["PMID:39753706:O1"]
    wo["intervention_arm_ids"] = ["VAR:PMID:39753706:A2_original_protocol_ctrl"]
    wo["intervention_arm_id"] = "VAR:PMID:39753706:A2_original_protocol_ctrl"
    wo["focal_arm_id"] = "VAR:PMID:39753706:A2_original_protocol_ctrl"
    wo["comparator_arm_ids"] = ["VAR:PMID:39753706:A1"]
    wo["comparison_arms"] = "A2_original_protocol_ctrl vs A1"
    with pytest.raises(ValueError, match="arm-role semantics"):
        bceg.verify_ceg_payload(bad_wang, repo_root=ROOT)
    bad_wang_dir = copy.deepcopy(ceg)
    bad_wang_dir["outcomes"]["PMID:39753706:O1"]["direction"] = "unfavorable_or_toxic"
    with pytest.raises(ValueError, match="direction="):
        bceg.verify_ceg_payload(bad_wang_dir, repo_root=ROOT)
    bad_wang_mag = copy.deepcopy(ceg)
    bad_wang_mag["outcomes"]["PMID:39753706:O2"]["endpoint_name"] = bad_wang_mag["outcomes"]["PMID:39753706:O2"]["endpoint_name"].replace("10-day", "2-day")
    bad_wang_mag["outcomes"]["PMID:39753706:O2"]["intervention_value"] = bad_wang_mag["outcomes"]["PMID:39753706:O2"]["intervention_value"].replace("10 days", "2 days")
    bad_wang_mag["outcomes"]["PMID:39753706:O2"]["quantitative_summary"] = bad_wang_mag["outcomes"]["PMID:39753706:O2"]["quantitative_summary"].replace("10-day", "2-day")
    with pytest.raises(ValueError, match="structured endpoint/magnitude"):
        bceg.verify_ceg_payload(bad_wang_mag, repo_root=ROOT)

    # Full-object CEG tamper probes (arm label, study biological_model, M0206:L2 claim_summary, M0206:O14 n_summary, Fu n_summary, Wang comparator_value, Rouhani M0206:O10 intervention_value/figure_or_table, and legacy M0154:O1)!
    bad_arm_lbl = copy.deepcopy(ceg)
    bad_arm_lbl["arms"]["VAR:M0189:A1"]["label"] = "Allogeneic rat kidney graft under immunosuppression"
    with pytest.raises(ValueError, match="arm VAR:M0189:A1 full-object digest"):
        bceg.verify_ceg_payload(bad_arm_lbl, repo_root=ROOT)
    bad_st_mod = copy.deepcopy(ceg)
    bad_st_mod["studies"]["M0189"]["biological_model"] = "Allogeneic rat kidney graft"
    with pytest.raises(ValueError, match="study M0189 full-object digest"):
        bceg.verify_ceg_payload(bad_st_mod, repo_root=ROOT)
    bad_m206_l1_txt = copy.deepcopy(ceg)
    for _ee in bad_m206_l1_txt["edges"]:
        if _ee["edge_id"] == "M0206:L1":
            _ee["claim_summary"] = "Tampered M0206:L1 claim summary"
    with pytest.raises(ValueError, match="edge M0206:L1 full-object digest"):
        bceg.verify_ceg_payload(bad_m206_l1_txt, repo_root=ROOT)
    bad_m206_l2_txt = copy.deepcopy(ceg)
    for _ee in bad_m206_l2_txt["edges"]:
        if _ee["edge_id"] == "M0206:L2":
            _ee["claim_summary"] = "Reversed claim summary"
    with pytest.raises(ValueError, match="edge M0206:L2 full-object digest"):
        bceg.verify_ceg_payload(bad_m206_l2_txt, repo_root=ROOT)
    bad_m206_nsum = copy.deepcopy(ceg)
    bad_m206_nsum["outcomes"]["M0206:O14"]["n_summary"] = "Tampered denominator"
    with pytest.raises(ValueError, match="M0206:O14 full-object digest"):
        bceg.verify_ceg_payload(bad_m206_nsum, repo_root=ROOT)
    bad_fu_nsum = copy.deepcopy(ceg)
    bad_fu_nsum["outcomes"]["PMID:41436438:O1"]["n_summary"] = "100 donors"
    with pytest.raises(ValueError, match="PMID:41436438:O1 full-object.*digest"):
        bceg.verify_ceg_payload(bad_fu_nsum, repo_root=ROOT)
    bad_wang_cmp = copy.deepcopy(ceg)
    bad_wang_cmp["outcomes"]["PMID:39753706:O2"]["comparator_value"] = "5-fold"
    with pytest.raises(ValueError, match="PMID:39753706:O2 full-object.*digest"):
        bceg.verify_ceg_payload(bad_wang_cmp, repo_root=ROOT)
    bad_rouhani_iv = copy.deepcopy(ceg)
    bad_rouhani_iv["outcomes"]["M0206:O10"]["intervention_value"] = "Tampered value"
    with pytest.raises(ValueError, match="M0206:O10 full-object digest"):
        bceg.verify_ceg_payload(bad_rouhani_iv, repo_root=ROOT)
    bad_rouhani_fig = copy.deepcopy(ceg)
    bad_rouhani_fig["outcomes"]["M0206:O10"]["figure_or_table"] = "Fig 99Z"
    with pytest.raises(ValueError, match="M0206:O10 full-object digest"):
        bceg.verify_ceg_payload(bad_rouhani_fig, repo_root=ROOT)
    bad_legacy_o = copy.deepcopy(ceg)
    bad_legacy_o["outcomes"]["M0154:O1"]["quantitative_summary"] = "Tampered qualitative summary"
    with pytest.raises(ValueError, match="M0154:O1 full-object digest"):
        bceg.verify_ceg_payload(bad_legacy_o, repo_root=ROOT)

    if has_dossiers:
        lin_abs = json.loads((ROOT / "data" / "dossiers" / "38843830.json").read_text(encoding="utf-8"))["abstract"]
        bad_sg = copy.deepcopy(ceg)
        bad_sg["arms"]["VAR:M0189:A1"]["label"] = "Prefix " + lin_abs[60:120] + " suffix."
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_sg, repo_root=ROOT)
        # Casefolded UPPERCASE _SourceGuard probe, legacy edge note & proposal-modified legacy outcome M0154:O1 _SourceGuard probes (MUST trip _SourceGuard BEFORE digest check)!
        bad_sg_upper = copy.deepcopy(ceg)
        bad_sg_upper["studies"]["M0189"]["follow_up"] = ("Prefix " + lin_abs[60:120] + " suffix.").upper()
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_sg_upper, repo_root=ROOT)
        bad_sg_leg = copy.deepcopy(ceg)
        for _el in bad_sg_leg["edges"]:
            if _el["edge_id"] == "M0051:L1":
                _el["legacy_certification_adjudication_note"] = "Reclassified note: " + lin_abs[60:120]
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_sg_leg, repo_root=ROOT)
        bad_sg_leg_out = copy.deepcopy(ceg)
        bad_sg_leg_out["outcomes"]["M0154:O1"]["quantitative_summary"] = "Copied prose: " + lin_abs[60:120]
        with pytest.raises(ValueError, match="_SourceGuard"):
            bceg.verify_ceg_payload(bad_sg_leg_out, repo_root=ROOT)

        # Standalone 11-dossier verification (ONLY the 11 IIP-0035 dossiers without the 5 IIP-0034 dossiers!)
        iso11 = tmp_path / "iso11"
        (iso11 / "data" / "dossiers").mkdir(parents=True)
        for pm11 in ("2642862", "31427818", "35953586", "36944335", "38843830", "39753706", "41201433", "41436438", "41690310", "41764741", "41999746"):
            shutil.copy2(ROOT / "data" / "dossiers" / f"{pm11}.json", iso11 / "data" / "dossiers" / f"{pm11}.json")
        det11 = bceg.verify_ceg_payload_detailed(ceg, repo_root=iso11)
        assert det11["verification_mode"] == "iip0035_11dossiers_8new_exact_spans_verified"
        assert len(det11["iip0035_dossiers_verified"]) == 11
        assert det11["iip0034_dossiers_verified"] == []
        # Replacing a legacy cross-link dossier with a PMID-only stub OR mutating dossier doi MUST fail closed!
        stub_p = iso11 / "data" / "dossiers" / "41201433.json"
        orig_412 = stub_p.read_bytes()
        stub_p.write_text('{"pmid": "41201433"}\n', encoding="utf-8")
        with pytest.raises(ValueError, match="IIP-0035 dossier 41201433.json full-object digest"):
            bceg.verify_ceg_payload_detailed(ceg, repo_root=iso11)
        stub_p.write_bytes(orig_412)
        fu_p = iso11 / "data" / "dossiers" / "41436438.json"
        orig_fu = fu_p.read_bytes()
        fu_tamper = json.loads(orig_fu.decode("utf-8"))
        fu_tamper["doi"] = "10.1038/fake"
        fu_p.write_text(json.dumps(fu_tamper, indent=2) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match="IIP-0035 dossier 41436438.json full-object digest"):
            bceg.verify_ceg_payload_detailed(ceg, repo_root=iso11)
        fu_p.write_bytes(orig_fu)
        # Dossier internal PMID mismatch vs filename MUST fail closed!
        d_tamper = json.loads((iso11 / "data" / "dossiers" / "2642862.json").read_text(encoding="utf-8"))
        d_tamper["pmid"] = "99999999"
        (iso11 / "data" / "dossiers" / "2642862.json").write_text(json.dumps(d_tamper, indent=2) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match="mismatches filename stem"):
            bceg.verify_ceg_payload_detailed(ceg, repo_root=iso11)
        # Partial 10-of-11 dossier set MUST still fail closed
        (iso11 / "data" / "dossiers" / "2642862.json").unlink()
        with pytest.raises(ValueError, match="Incomplete local dossier set"):
            bceg.verify_ceg_payload_detailed(ceg, repo_root=iso11)


def test_08_schema_2160_stack_invariance_lims_0byte_and_public_reexport(tmp_path: Path) -> None:
    schema = json.loads((ROOT / "schemas" / "extraction_schema_v0.2.json").read_text(encoding="utf-8"))
    props = schema["properties"]
    assert "Somatic_Mutation_and_Genomic_Stability" in props["atomic_effects"]["items"]["properties"]["outcome_domain"]["enum"]
    assert "Tier_2C_NHP_Allogeneic_or_Autologous_NonIslet" in props["study_metadata"]["properties"]["max_biological_fidelity_tier"]["enum"]
    assert "NHP_iPSC_Cardiomyocyte_NonIslet" in props["arms"]["items"]["properties"]["gate1_substrate_and_differentiation"]["properties"]["cell_source_class"]["enum"]
    ch_enum = props["arms"]["items"]["properties"]["challenge_conditions"]["items"]["properties"]["challenge_type"]["enum"]
    assert "Peptide_Restimulation_Islet_Derived_T_Cell_Line_Assay" in ch_enum
    assert "Autologous_T_Cell_Neoantigen_InVitro_Assay" in ch_enum
    for u in ("sequenced_hiPSC_and_subclone_samples", "pancreas_transplants"):
        assert u in props["arms"]["items"]["properties"]["sample_size_unit"]["enum"]
        assert u in props["atomic_effects"]["items"]["properties"]["n_unit"]["enum"]

    # HIP-Core audit ledger corpus_eligibility synchronized to 604 / 589 / 570 (zero stale 594/579/560!)
    hip_doc = json.loads((ROOT / "data" / "protocols" / "hip_core_fractional_factorial_v0.6.2.json").read_text(encoding="utf-8"))
    l1_led = hip_doc["layer1_ultracode_audit_reconciliation_ledger"]
    assert l1_led["source_counts"]["studies"] == 604
    elig_rec = l1_led["corpus_eligibility"]["ring1_corpus_count_reconciliation"]
    assert elig_rec["validated_studies_total"] == 604
    assert elig_rec["active_deduplicated_studies_total"] == 589
    assert elig_rec["strict_active_peer_reviewed_empirical_studies_total"] == 570

    # Protected LIMS 0-byte invariance & untouched 366 / 77 / 122 counts
    for rel, exp_sha in EXPECTED_LIMS_SHA256.items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == exp_sha
    ct = json.loads((ROOT / "data" / "protocols" / "wetlab_endpoint_contract_v0.8.0.json").read_text(encoding="utf-8"))
    assert (len(ct["stage0_decision_tree_endpoints"]), len(ct["metadata_keys"]), len(ct["boolean_policy_keys"])) == (366, 77, 122)

    # 0 / 2,160 Layer-1 stack invariance
    rank_mod = _load_mod("scripts/rank_islet1_stacks.py", "openislet_rank_inv_iip0035_test")
    _, studies_all, omics = rank_mod.load_inputs()
    new10 = {"2642862", "31427818", "35953586", "36944335", "38843830", "39753706", "41436438", "41690310", "41764741", "41999746"}
    studies_pre = [s for s in studies_all if str(s.get("study_metadata", {}).get("pmid")) not in new10]
    c_post = rank_mod.run_combinatorial_multi_gate_search(rank_mod.build_target_evidence_profiles(studies_all, omics))
    c_pre = rank_mod.run_combinatorial_multi_gate_search(rank_mod.build_target_evidence_profiles(studies_pre, omics))
    assert len(c_post) == len(c_pre) == 2160
    assert sum(1 for a, b in zip(c_post, c_pre) if a != b) == 0

    # Source-free public re-export preserves all 604 studies, 262 expanded genes, 51 profiles & passes verifier --check!
    grm_path = ROOT / "scripts" / "generate_release_manifest.py"
    if not grm_path.is_file():
        pytest.skip("generate_release_manifest.py not present")
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    reexp_out = tmp_path / "reexp_pub"
    if (ROOT / ".openislet_public_bundle").is_file() and ((ROOT / "data" / "dossiers").exists() or (ROOT / "data" / "curated").exists()):
        # Reviewer unpacked partial private dossiers/curated into a sealed public bundle copy;
        # stage the clean sealed public bundle files so re-export exercises the authentic sealed path.
        clean_pub = tmp_path / "clean_pub"
        for p in sorted(ROOT.rglob("*")):
            if not p.is_file() or "__pycache__" in p.parts:
                continue
            rel_s = p.relative_to(ROOT).as_posix()
            if rel_s.startswith(("data/dossiers/", "data/curated/")):
                continue
            dst = clean_pub / rel_s
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dst)
        src_free = clean_pub
    elif (ROOT / "data" / "dossiers").is_dir():
        pub1 = tmp_path / "pub1"
        subprocess.run([sys.executable, "-B", str(grm_path), "--export-public-bundle", str(pub1)], cwd=str(ROOT), env=env, check=True, capture_output=True)
        src_free = pub1
    else:
        src_free = ROOT
    subprocess.run([sys.executable, "-B", str(src_free / "scripts" / "generate_release_manifest.py"), "--export-public-bundle", str(reexp_out)], cwd=str(src_free), env=env, check=True, capture_output=True)
    grm_re = _load_mod("scripts/generate_release_manifest.py", "openislet_grm_reexp_iip0035")
    assert grm_re.verify_public_bundle_snapshot(reexp_out)["passed"] is True
    assert len(list((reexp_out / "data" / "extracted").glob("[0-9]*.json"))) == 604
    # Idempotent byte-identical preservation across ALL public files in sealed public bundle re-export
    # (including 17375075.json, 39571901.json, and Sutherland 2642862.json factual arm_label!)
    src_files = {p.relative_to(src_free).as_posix(): p.read_bytes() for p in src_free.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    re_files = {p.relative_to(reexp_out).as_posix(): p.read_bytes() for p in reexp_out.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    assert set(src_files) == set(re_files)
    drifted = [k for k in sorted(src_files) if src_files[k] != re_files[k]]
    assert drifted == [], f"Sealed public bundle re-export modified {len(drifted)} file(s): {drifted[:5]}"
    assert "REDACTED" not in str(json.loads((reexp_out / "data" / "extracted" / "2642862.json").read_text(encoding="utf-8"))["arms"][0]["arm_label"])
    res_re_v = subprocess.run(
        [sys.executable, "-B", str(reexp_out / "scripts" / "verify_iip0035_d1b_auto_and_portfolio.py"), "--check"],
        cwd=str(reexp_out), env=env, capture_output=True, text=True, check=True,
    )
    assert "OK:" in res_re_v.stdout
    # Tampering .openislet_public_bundle marker bytes MUST fail verify_public_bundle_snapshot!
    (reexp_out / ".openislet_public_bundle").write_bytes(b"wrong_marker_v9\n")
    with pytest.raises(ValueError, match="Manifest snapshot integrity failed"):
        grm_re.verify_public_bundle_snapshot(reexp_out)
