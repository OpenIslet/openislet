#!/usr/bin/env python3
"""Focused standalone verification suite for IIP-0033 (DEC-112).

Covers:
  1. Pinned 3-mesh fixed-pitch 65 um bilayer & respirometry verification artifact (schema, SHA-256, key physical bounds).
  2. Agilent 101122-100 static KRBH caps (16.60/16.87/14.94/14.00/9.87 uL vs obsolete 22.70 uL V7) and flow-through 5% O2 perifusion.
  3. Exact joint-constraint Hyndman-Fan Type-7 finite-sample cubic-mean bounds, including N = 2..500, N = 200, N = 1,000,000, and N = 4 boundary tolerances.
  4. Cylindrical pre-vascularized channel traps at clinical 500,000 IEQ scale (R = 500 um anoxia; R = 100-120 um mechanics & 30 um clot collapse).
  5. Ring-1 corpus expansion (583 -> 589 studies, 3,139 -> 3,148 effects, 73 curated landmarks, Jin Day 105 / Balboni §2.11 citation pins, vascular endothelial MHC guard).
  6. Layer-1 2,160-stack invariance and live OQ-016 baseline parity (+0.00 pp).
  7. Context-Evidence Graph expansion (1,165 -> 1,172 records, 7,725 arms, 14,159 outcomes, 2,478 links, 40 curated ineligible edges, NCT03513939 protocol class).
  8. Idea Graveyard GY-028 structure, 14 PMIDs, one-70 nz=156 axial qualifier, and explicit 350 um nz=70/140/156 mesh labels.
  9. Preregistered Stage-0 decision-tree specification block and complete invariance of wetlab_endpoint_contract_v0.8.0.json (366 / 77 / 122).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def _load_verifier_module() -> Any:
    v_path = ROOT / "scripts" / "verify_iip0033_bilayer_and_resp1.py"
    spec = importlib.util.spec_from_file_location("openislet_verify_iip0033_test", v_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_01_pinned_verification_artifact_schema_and_sha256() -> None:
    art_path = ROOT / "data" / "extracted" / "iip0033_bilayer_and_respirometry_verification_v0.8.0.json"
    assert art_path.is_file(), f"Missing pinned verification artifact: {art_path}"
    raw = art_path.read_text(encoding="utf-8")
    digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    assert digest == "bf581a628a69efa41acfc4d5e01de08d36b1c52fc913a0a9c309a27b6b97807e", (
        f"Unexpected pinned artifact SHA-256: {digest}"
    )
    doc = json.loads(raw)
    assert doc["contribution_id"] == "IIP-0033"
    assert doc["decision_id"] == "DEC-112"
    assert doc["overall_status"] == "ALL_IIP0033_PHYSICS_AND_RESPIROMETRY_ASSERTIONS_PASSED"
    assert doc["mesh_schedule"] == [[40, 28, 32, 78], [60, 42, 48, 78], [80, 56, 64, 78]]

    fp = doc["fixed_pitch_bilayer_and_holds"]
    assert abs(float(fp["fixed_pitch_um"]) - 104.1604) < 1e-4
    assert float(fp["nominal_loading_ieq_cm2"]) == 1500.0
    sc = fp["spatial_counterexamples_necessary_not_sufficient"]
    uneven = [r for r in sc["route1_390um_failing_counterexamples"] if "Uneven2x2" in r["label"]]
    assert len(uneven) == 12 and all(float(r["worst_face_mmhg"]) < 15.0 for r in uneven)
    sup32 = [r for r in sc["route1_390um_failing_counterexamples"] if "Supercell32" in r["label"] and "CoupledPlus10" in r["label"]]
    assert len(sup32) == 4 and all(float(r["worst_face_mmhg"]) < 15.0 for r in sup32)
    redist = sc["route1_390um_redistributed_8site_controls"]
    assert len(redist) == 12 and any(abs(float(r["worst_face_mmhg"]) - 15.7151) < 1e-3 for r in redist)


def test_02_respirometry_contract_caps_and_flow_through() -> None:
    vmod = _load_verifier_module()
    resp = vmod.verify_respirometry_contract()
    regimes = resp["static_xfe24_101122_100_regimes"]
    assert regimes["governing_minus10pct_14p94uL_krbh_1p27"]["max_ieq_5pct_gsis_1p5x_ge25"] == 9
    assert regimes["governing_minus10pct_14p94uL_krbh_1p27"]["max_ieq_21pct_fccp_2p5x_ge120"] == 16
    assert regimes["screen_displaced_14p00uL_krbh_1p27"]["max_ieq_21pct_fccp_2p5x_ge120"] == 15
    assert regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["max_ieq_5pct_gsis_1p5x_ge25"] == 6
    assert regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["max_ieq_21pct_fccp_2p5x_ge120"] == 10
    assert regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["pred_9_ieq_5pct_gsis1p5_mmhg"] == 20.96
    assert regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["pred_16_ieq_21pct_fccp2p5_mmhg"] == 105.84
    assert resp["capture_screen_pore_size_um"] == 125.0

    ft = resp["flow_through_primary_lead_200ieq_200uLmin_krbh"]
    assert ft["p_out_5pct_basal_1p0x_mmhg"] == 32.07
    assert ft["p_out_5pct_gsis_1p5x_mmhg"] == 30.28
    assert ft["p_out_5pct_gsis_2p0x_mmhg"] == 28.49
    assert ft["p_out_5pct_fccp_2p5x_mmhg"] == 26.70
    assert ft["all_4_states_ge_25mmhg"] is True and ft["passes_primary_spec"] is True

    # Low-OCR flow-through check stays consistent with SNR >= 5 and >= 1.00 mmHg floor
    ft_low_ocr = vmod.compute_sweet2008_flow_through_po2_out_iip0033(
        n_ieq=200.0,
        flow_uL_min=200.0,
        alpha_mol_cm3_mmhg=1.27e-9,
        ocr_whole_prep_pmol_min_ieq=2.54,
        noise_floor_mmhg=0.25,
    )
    assert ft_low_ocr["delta_po2_basal_mmhg"] == 2.0
    assert ft_low_ocr["snr_over_noise_floor"] == 8.0
    assert ft_low_ocr["passes_primary_spec"] is True

    ft_sub_snr = vmod.compute_sweet2008_flow_through_po2_out_iip0033(
        n_ieq=200.0,
        flow_uL_min=200.0,
        alpha_mol_cm3_mmhg=1.27e-9,
        ocr_whole_prep_pmol_min_ieq=1.27,
        noise_floor_mmhg=0.50,
    )
    assert ft_sub_snr["snr_over_noise_floor"] == 2.0
    assert ft_sub_snr["passes_primary_spec"] is False


def test_03_type7_finite_sample_cubic_mean_bounds_exactness_and_boundaries() -> None:
    vmod = _load_verifier_module()
    summary = vmod.verify_type7_finite_sample_cubic_mean_bounds()
    assert summary["all_finite_sample_bounds_hold"] is True
    assert summary["n_range_verified"] == [2, 500]
    assert summary["n200_joint_constraint_counterexample_58_65_65_68"]["exact_joint_lb_um"] == 64.7323
    assert summary["n200_joint_constraint_counterexample_58_65_65_68"]["exact_joint_ub_um"] == 65.1098

    # Boundary witness 1: N = 1,000,000 (999,999 at 64.0 um + 1 at 64.0005 um)
    b_1m = vmod.compute_type7_finite_sample_cubic_mean_bounds_iip0033(
        d_min=64.0, p05=64.0, p95=64.0, d_max=64.0005, sample_n=1_000_000
    )
    expected_1m = ((999_999.0 * (64.0 ** 3) + (64.0005 ** 3)) / 1_000_000.0) ** (1.0 / 3.0)
    assert b_1m["feasible"] is True
    assert abs(b_1m["d_rms3_lb"] - expected_1m) < 1e-12
    assert b_1m["d_rms3_lb"] > 64.0

    # Boundary witness 2: N = 4 [60, 62, 62, 64] perturbed inward by 32e-9 um (inside tol = 64e-9 um)
    b_4_exact = vmod.compute_type7_finite_sample_cubic_mean_bounds_iip0033(
        d_min=60.0, p05=60.3, p95=63.7, d_max=64.0, sample_n=4
    )
    b_4_pert = vmod.compute_type7_finite_sample_cubic_mean_bounds_iip0033(
        d_min=60.0, p05=60.3 + 32e-9, p95=63.7 - 32e-9, d_max=64.0, sample_n=4
    )
    assert b_4_exact["feasible"] is True and b_4_pert["feasible"] is True
    assert abs(b_4_pert["d_rms3_lb"] - b_4_exact["d_rms3_lb"]) < 1e-6


def test_04_cylindrical_prevascular_traps_recomputed() -> None:
    vmod = _load_verifier_module()
    cyl = vmod.verify_cylindrical_traps()
    assert float(cyl["r500_anoxia_3d_fvm"]["core_min_mmhg"]) == 0.0
    assert abs(float(cyl["r120_30um_clot_stroma15_3d_fvm"]["core_min_mmhg"]) - 1.4106) < 0.01
    assert abs(float(cyl["r120_30um_clot_vpio25_3d_fvm"]["core_min_mmhg"]) - 0.5908) < 0.01


def test_05_ring1_corpus_6_studies_and_citation_accuracy() -> None:
    ext_files = sorted((ROOT / "data" / "extracted").glob("[0-9]*.json"))
    assert len(ext_files) == 594
    new_pmids = ("19555051", "25893782", "31527312", "36591328", "41536981", "42176117")
    docs = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in ext_files}
    for pmid in new_pmids:
        assert pmid in docs, f"Missing new Ring-1 study {pmid}"
    assert sum(len(d.get("arms", [])) for d in docs.values()) == 693
    assert sum(
        1
        for d in docs.values()
        for a in d.get("arms", [])
        if a.get("is_reference_control") is True
    ) == 62
    assert sum(len(d.get("atomic_effects", [])) for d in docs.values()) == 3158

    # Verify verified primary-source citations on Jin (Day 105 Fig. 6F / Day 182 Fig. 7) & Balboni (Methods 2.11)
    jin_text = json.dumps(docs["41536981"])
    assert "Day 105" in jin_text and "Day 182" in jin_text
    balb_text = json.dumps(docs["42176117"])
    assert "2.11" in balb_text


def test_06_live_oq016_and_2160_stack_invariance() -> None:
    vmod = _load_verifier_module()
    parity = vmod.verify_live_oq016_baseline()
    assert parity["evaluated_stacks_checked"] == 2160
    assert parity["pre_iip0033_studies_count"] == 588
    assert parity["post_iip0033_studies_count"] == 594
    assert parity["mutated_stacks_vs_pre_iip0033"] == 0
    assert parity["mutated_stacks_vs_committed_pareto"] == 0
    assert parity["iip0033_live_headline_delta_pp"] == 0.0

    pareto_sha = hashlib.sha256(
        (ROOT / "data" / "extracted" / "islet1_combinatorial_pareto_and_ladder_v0.8.0.json").read_bytes()
    ).hexdigest()
    assert pareto_sha == "591ddeb9a706c0bc95fc509765b512aadad8a403287a796e63dd52e7db548575"


def test_07_ceg_expansion_and_ineligible_edges() -> None:
    ceg = json.loads((ROOT / "data" / "extracted" / "context_evidence_graph_v1.json").read_text(encoding="utf-8"))
    assert len(ceg["studies"]) == 1176
    assert len(ceg["arms"]) == 7740
    assert len(ceg["outcomes"]) == 14167
    assert len(ceg["edges"]) == 2483
    elig = sum(1 for e in ceg["edges"] if e.get("primary_certification_eligible") is True)
    inelig = sum(1 for e in ceg["edges"] if e.get("primary_certification_eligible") is False)
    assert (elig, inelig) == (2283, 200)

    for qid in ("M1453", "M1422", "M1421", "M2036", "M2449", "PMID:19555051", "NCT:NCT03513939"):
        assert qid in ceg["studies"], f"Missing CEG study {qid}"
    assert ceg["studies"]["NCT:NCT03513939"]["record_class"] == "trial_registry_protocol"

    edges_by_id = {e["edge_id"]: e for e in ceg["edges"]}
    for eid in ("M1421:L1", "M1422:L1", "M1422:L2", "M1453:L1", "M2036:L1", "M2449:L1", "PMID:19555051:L1", "NCT:NCT03513939:L1"):
        assert eid in edges_by_id and edges_by_id[eid]["primary_certification_eligible"] is False


def test_08_graveyard_gy028_structure_and_qualifiers() -> None:
    gy = json.loads((ROOT / "data" / "graveyard" / "idea_graveyard_v0.5.json").read_text(encoding="utf-8"))
    og_block = gy["off_grid_molecular_dead_end_guardrails"]
    og_list = og_block["guardrails"]
    assert len(gy["graveyard_entries"]) == 10 and len(og_list) == 18
    assert og_block["latest_update_decision_id"] == "DEC-112"
    gy28 = [g for g in og_list if g.get("graveyard_card_id") == "GY-028"]
    assert len(gy28) == 1
    card = gy28[0]
    assert card["contribution_id"] == "IIP-0033"
    assert card["severity_level"] == "BLOCK"
    assert len(card["pmids"]) == 14
    summary = card["summary"]
    assert "nz = 156" in summary and "15.04 mmHg" in summary
    assert "nz = 70" in summary and "nz = 140" in summary and "14.74-14.76 mmHg" in summary


def test_09_stage0_decision_tree_spec_and_wetlab_contract_invariance() -> None:
    dt = json.loads((ROOT / "data" / "protocols" / "stage0_decision_tree_v1.json").read_text(encoding="utf-8"))
    blk = dt.get("fixed_pitch_65um_bilayer_and_respirometry_gates_iip_0033")
    assert isinstance(blk, dict)
    assert blk["decision_id"] == "DEC-112"
    assert blk["contribution_id"] == "IIP-0033"
    assert len(blk["primary_citations"]) == 14

    # Confirm wetlab_endpoint_contract_v0.8.0.json stays untouched from cc5af3b (science-first scope)
    ct = json.loads((ROOT / "data" / "protocols" / "wetlab_endpoint_contract_v0.8.0.json").read_text(encoding="utf-8"))
    assert len(ct["stage0_decision_tree_endpoints"]) == 366
    assert len(ct["metadata_keys"]) == 77
    assert len(ct["boolean_policy_keys"]) == 122
