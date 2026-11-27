#!/usr/bin/env python3
"""
Deterministic Planar Route 3P Physics, Colony-Screening, Exact Decimal DEC-103 R3P-05 Thinning,
CommonMark Playbook & Typed 96-Cell Gate Contract Verifier (IIP-0034 / DEC-113)
===============================================================================================
Self-contained replay script for OpenIslet public and audit checkouts. Depends ONLY on
`scripts/run_layer4_tissue_and_metabolic_twin.py`, `scripts/rank_islet1_stacks.py`,
`scripts/product_contract.py`, `numpy`, and the Python standard library (zero dependency on
`scripts/wetlab_contract.py`).

Run with:
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0034_planar_and_portfolio.py
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0034_planar_and_portfolio.py --write
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0034_planar_and_portfolio.py --check
"""
from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
import hashlib
import html
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import product_contract as pc  # noqa: E402
import rank_islet1_stacks as rank_mod  # noqa: E402
import run_layer4_tissue_and_metabolic_twin as l4  # noqa: E402

DEFAULT_OUT_PATH = ROOT / "data" / "extracted" / "iip0034_planar_and_portfolio_verification_v0.8.0.json"
PORTFOLIO_PATH = ROOT / "data" / "extracted" / "islet1_multiparadigm_portfolio_v0.8.0.json"
PARETO_REL_PATH = "data/extracted/islet1_combinatorial_pareto_and_ladder_v0.8.0.json"
DECISION_TREE_REL_PATH = "data/protocols/stage0_decision_tree_v1.json"
PLAYBOOK_REL_PATH = "docs/STAGE-0-EXECUTION-PLAYBOOK.md"
EXPECTED_PARETO_SHA256 = "591ddeb9a706c0bc95fc509765b512aadad8a403287a796e63dd52e7db548575"
EXPECTED_DEC113_DECISION_TREE_SHA256 = "46609af00da04f998d26b8ffe47f43350587989ef40e9be843ba6849198c5952"
EXPECTED_DEC113_PLAYBOOK_SEC4119_SHA256 = "7511695ef9dbade1f9a231e3ee522d6c8b749ab8a3c0e6962ac1271a147d1004"
R3P05_MAX_U95_CEILING_DEC = Decimal("1e-6")
R3P05_MAX_U95_CEILING = 1e-6
ELIGIBLE_R3P05_FULLSTACK_ROLES = frozenset({"FULL_CHR15_11TU_REFERENCE", "GATE0V_CD55_CFH_BACKUP_COMPARATOR"})
ASSEMBLED_DEVICE_LOT_SCOPE_KEYS = (
    "product_variant_id",
    "construct_sha256",
    "donor_seedline_id",
    "mcb_id",
    "lot_id",
    "device_assembly_id",
    "physical_route_id",
)
ASSEMBLED_DEVICE_LOT_SCOPE_KEY_SET = frozenset(ASSEMBLED_DEVICE_LOT_SCOPE_KEYS)
VALID_F_CLONOGENIC_PROVENANCE = frozenset({
    "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
    "INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
})
EXPECTED_96_GATES = (
    "Y-01", "Y-02", "Y-03", "Y-04", "Y-05", "Y-06", "Y-07", "Y-08", "Y-09",
    "R3P-01", "R3P-02", "R3P-03", "R3P-04", "R3P-05", "R3P-06",
    "STAGE0_ENDPOINT_CROSSWALK",
)
EXPECTED_DEC102_ENDPOINT_IDS = frozenset({
    "DEC102_MODALITY_POLARITY_AND_EXPLICIT_SE",
    "DEC102_BOOL_ATTLR_RACE_VERIFIED",
    "DEC102_BOOL_ZERO_2A_SCAR_PROTEOMICS",
    "DEC102_BOOL_DPCPX_ANTAGONIST_CONTROL",
    "DEC102_BOOL_ANTI_CD200_BLOCKADE_REVERSAL",
    "DEC102_BOOL_Y04_FULL_COCULTURE_CLONOGENIC",
    "Y-01_SAMPLING_AGE_60D",
    "Y-01_ALL4_LINEAGES_SCORED",
    "Y-01_LINEAGE_CELLS_SCORED_10K",
    "Y-01_MAX_LINEAGE_ESCAPE_1PCT",
    "Y-01_TADCBED_CBE_OFFTARGET_0P10PCT",
    "Y-02_COUPON_KI_3PCT",
    "Y-02_CLONAL_DDPCR_CN_RANGE",
    "Y-02_ONT_INTACT_HOMOLOG_99PCT",
    "Y-02_WT_B2M_RESIDUAL_0P001",
    "Y-02_BACKBONE_KANR_VAF_1E4",
    "Y-02_WGS_PSEUDOSITE_K0",
    "Y-03_IEC_BETA_SELECTIVITY_50TO1",
    "Y-03_BETA_MIR375_RETENTION_90PCT",
    "Y-03_COCULTURE_GSIS_SI_2P5",
    "Y-03_GSIS_RETENTION_90PCT",
    "Y-03_SUPERNATANT_ADENOSINE_25NM",
    "Y-03_DPCPX_GSIS_DELTA_10PCT",
    "Y-04_COUPON_RETENTION_95PCT",
    "Y-04_COMPLEMENT_C5B9_5PCT",
    "Y-04_INDIRECT_APC_CD4_SI_1P5",
    "Y-04_KIR_NKG2A_NEG_NK_8PCT",
    "Y-04_ICASP9_CLONOGENIC_ESCAPE_1E8",
    "Y-04_STANDALONE_GSIS_SI_3P0",
    "DEC102_WILDCARD_20PCT_FIREWALL",
})
REQUIRED_XWALK_KEYS = frozenset({
    "dec102_endpoint_id",
    "wetlab_contract_key",
    "disposition",
    "replacement_or_binding",
    "operator",
    "threshold",
    "unit",
    "threshold_and_basis",
})
VALID_XWALK_OPERATORS = frozenset({"<=", ">=", "==", "IN_RANGE"})
VALID_XWALK_DISPOSITIONS = frozenset({
    "RETAINED",
    "SUPERSEDED_BY_DEC113_Y05",
    "RETAINED_WITH_DEC113_Y06_SPLIT",
})
VALID_CELL_OPERATORS = frozenset({
    "<=",
    ">=",
    "CONJUNCTIVE_PASS",
    "PROTOCOL_DESIGN_RULE",
    "COMPARATOR_REFERENCE_ONLY",
    "EXCLUDED_FROM_FULLSTACK_CLEARANCE",
    "NOT_APPLICABLE",
})
VALID_CELL_DISPOSITIONS = frozenset({
    "REQUIRED_EXACT_PRODUCT",
    "REQUIRED_AFTER_PRIMARY_GATE0V_FAIL",
    "REQUIRED_AFTER_PRIMARY_GATE0V_FAIL_AND_DELIVERY_HOLD",
    "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK",
    "CONTROL_BENCHMARK_ONLY_EXCLUDED_FROM_FULLSTACK",
    "NOT_APPLICABLE_EXCLUDED_FROM_FULLSTACK",
})
EXPECTED_FULLSTACK_GATE_OPERATORS: Dict[str, str] = {
    "Y-01": "<=",
    "Y-02": "CONJUNCTIVE_PASS",
    "Y-03": "CONJUNCTIVE_PASS",
    "Y-04": "<=",
    "Y-05": "CONJUNCTIVE_PASS",
    "Y-06": "<=",
    "Y-07": "CONJUNCTIVE_PASS",
    "Y-08": "CONJUNCTIVE_PASS",
    "Y-09": "CONJUNCTIVE_PASS",
    "R3P-01": "CONJUNCTIVE_PASS",
    "R3P-02": "<=",
    "R3P-03": ">=",
    "R3P-04": "CONJUNCTIVE_PASS",
    "R3P-05": "<=",
    "R3P-06": "PROTOCOL_DESIGN_RULE",
    "STAGE0_ENDPOINT_CROSSWALK": "CONJUNCTIVE_PASS",
}
_PLACEHOLDER_TOKEN_RE = re.compile(r"(?i)\b(?:TBD|TODO|PLACEHOLDER|FIXME|XXX|PENDING_VALUE|DUMMY)\b|\?\?\?")
REQUIRED_CELL_KEYS = frozenset({
    "product_variant_id",
    "gate_id",
    "disposition",
    "external_ref",
    "measure",
    "operator",
    "threshold",
    "unit",
    "confidence_method",
    "sample_unit",
    "scope",
    "precondition",
    "missing_data_disposition",
})
_FORBIDDEN_ID_CHARS_RE = re.compile(r"[\x00-\x1f\x7f\u200b-\u200f\u2028-\u202f\u2060-\u206f\ufeff]")

PHYSICAL_ROUTE_TO_PORTFOLIO_CODE: Dict[str, str] = {
    "ROUTE_1_MACROENCAPSULATED": "R1",
    "ROUTE_2_OPEN_PORE_PREVASC": "R2",
    "ROUTE_3P_MODULAR_PLANAR_ENVELOPE": "R3P",
}

# Canonical OpenIslet Layer-4 physical constants
ALPHA_LOW_SALT = float(l4.ALPHA_BUNSEN_MOL_CM3_MMHG)  # 1.34e-9 mol/(cm^3 * mmHg)
ALPHA_PHYS_BUFFER = 1.27e-9                           # physiological-buffer sensitivity solubility mol/(cm^3 * mmHg)
D_GEL_CHRON = 1.50e-5                                 # cm^2/s (macroporous fibrin-free slot)
D_ISLET = float(l4.D_O2_ISLET_TISSUE_CM2_S)           # 1.30e-5 cm^2/s
VMAX_PURE = float(l4.VMAX_ISLET_PURE_MOL_CM3_S)       # 3.4e-8 mol/(cm^3 * s)
VOL_IEQ = float(l4.VOL_1_IEQ_CM3)                     # 1.7671458676442586e-6 cm^3
DOSE_IEQ = 500_000.0


def _reject_duplicate_json_keys(pairs: List[tuple[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise ValueError(f"Duplicate JSON key rejected in controlling document: {k!r}")
        out[k] = v
    return out


def _reject_nonfinite_json(token: str) -> Any:
    raise ValueError(f"Non-standard JSON constant rejected: {token!r}")


def _parse_finite_json_float(token: str) -> float:
    val = float(token)
    if not math.isfinite(val):
        raise ValueError(f"Non-finite JSON float overflow rejected: {token!r}")
    return val


def _load_json_strict(path: Path) -> Any:
    return json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_json_keys,
        parse_constant=_reject_nonfinite_json,
        parse_float=_parse_finite_json_float,
    )


def _safe_exact_decimal(val: Any) -> Decimal | None:
    """Convert non-bool int/float/Decimal to exact finite Decimal without raising OverflowError."""
    if isinstance(val, bool):
        return None
    if isinstance(val, Decimal):
        return val if val.is_finite() else None
    if isinstance(val, int):
        return Decimal(val)
    if isinstance(val, float):
        if not math.isfinite(val):
            return None
        return Decimal.from_float(val)
    return None


def _is_canonical_identifier(val: Any) -> bool:
    return (
        isinstance(val, str)
        and bool(val)
        and val == val.strip()
        and _FORBIDDEN_ID_CHARS_RE.search(val) is None
    )


def _check_canonical_scope_integrity(
    expected_scope: Any,
    observed_scope: Any,
) -> tuple[bool, str, str | None]:
    """Validate canonical exact 7-key DEC-103 scope, eligible full-stack product binding, and identifier syntax."""
    if not isinstance(expected_scope, dict) or not isinstance(observed_scope, dict):
        return False, "MISSING_FULL_7KEY_ASSEMBLED_DEVICE_LOT_SCOPE", None
    if (
        set(expected_scope.keys()) != ASSEMBLED_DEVICE_LOT_SCOPE_KEY_SET
        or set(observed_scope.keys()) != ASSEMBLED_DEVICE_LOT_SCOPE_KEY_SET
    ):
        return False, "NON_EXACT_7KEY_SCOPE_MISSING_OR_UNEXPECTED_EXTRA_KEYS", None
    for key in ASSEMBLED_DEVICE_LOT_SCOPE_KEYS:
        ev = expected_scope.get(key)
        ov = observed_scope.get(key)
        if not _is_canonical_identifier(ev) or not _is_canonical_identifier(ov) or ev != ov:
            return False, f"SCOPE_KEY_MISSING_NONCANONICAL_OR_MISMATCHED:{key}", None
    if expected_scope["physical_route_id"] != "ROUTE_3P_MODULAR_PLANAR_ENVELOPE":
        return False, "PHYSICAL_ROUTE_NOT_ROUTE_3P_MODULAR_PLANAR_ENVELOPE", None
    eligible_products = {
        pid: str(spec["donor_sha256"])
        for pid, spec in pc.DESIGN_Y_V2_DONOR_BINDINGS.items()
        if str(spec.get("role") or "") in ELIGIBLE_R3P05_FULLSTACK_ROLES
    }
    pvid = expected_scope["product_variant_id"]
    if pvid not in eligible_products:
        return False, f"PRODUCT_NOT_ELIGIBLE_FULLSTACK_DESIGN_Y_V2:{pvid}", None
    if expected_scope["construct_sha256"] != eligible_products[pvid]:
        return False, "CONSTRUCT_SHA256_MISMATCH_OR_NON_DIGEST", None
    return True, "OK", pvid


def _parse_commonmark_headings(text: str) -> List[Dict[str, Any]]:
    """Stdlib CommonMark 0.31.2 block parser extracting all rendered ATX and multiline Setext headings.

    Properly ignores fenced code blocks (backtick/tilde >= 3 chars), indented code blocks (>= 4 spaces or tab),
    and HTML comments (<!-- ... -->), while recursively stripping blockquote ('>') and list-item container
    prefixes so headings inside blockquotes or list items are still detected as rendered headings.
    """
    norm = text.replace("\r\n", "\n").replace("\r", "\n")
    raw_lines = norm.split("\n")
    line_offsets: List[int] = []
    pos = 0
    for ln in raw_lines:
        line_offsets.append(pos)
        pos += len(ln) + 1

    headings: List[Dict[str, Any]] = []
    in_fence = False
    fence_char = ""
    fence_len = 0
    in_html_comment = False
    para_lines: List[tuple[int, str, int]] = []  # (line_idx, stripped_text, container_depth)

    atx_re = re.compile(r"^ {0,3}(#{1,6})(?:[ \t]+(.*?))?[ \t]*$")
    fence_open_re = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
    setext_re = re.compile(r"^ {0,3}(=+|-+)[ \t]*$")
    bq_re = re.compile(r"^ {0,3}> ?(.*)$")
    ul_re = re.compile(r"^ {0,3}[-+*][ \t]+(.*)$")
    ol_re = re.compile(r"^ {0,3}\d{1,9}[.)][ \t]+(.*)$")

    def _strip_containers(line: str) -> tuple[str, int]:
        depth = 0
        cur = line
        while True:
            mb = bq_re.match(cur)
            if mb:
                depth += 1
                cur = mb.group(1)
                continue
            mu = ul_re.match(cur)
            if mu:
                depth += 1
                cur = mu.group(1)
                continue
            mo = ol_re.match(cur)
            if mo:
                depth += 1
                cur = mo.group(1)
                continue
            break
        return cur, depth

    def _clean_atx_title(raw_title: str) -> str:
        t = raw_title.strip()
        # Strip optional CommonMark closing '#' sequence preceded by space/tab
        t = re.sub(r"(?:^|[ \t])#+[ \t]*$", "", t).strip()
        return t

    def _render_cm_inline_text(raw_title: str) -> str:
        # Decode CommonMark backslash-escaped ASCII punctuation and HTML character references (&#46;, &#x2e;, etc.)
        unesc = re.sub(r'\\([!"#$%&\'()*+,\-./:;<=>?@\[\\\]^_`{|}~])', r"\1", raw_title)
        return html.unescape(unesc)

    for idx, raw_line in enumerate(raw_lines):
        # 1. Active HTML comment block (CommonMark Type 2)
        if in_html_comment:
            if "-->" in raw_line:
                in_html_comment = False
            para_lines.clear()
            continue

        # 2. Active fenced code block
        if in_fence:
            close_pat = rf"^ {{0,3}}{re.escape(fence_char)}{{{fence_len},}}[ \t]*$"
            if re.match(close_pat, raw_line):
                in_fence = False
            para_lines.clear()
            continue

        # Strip blockquote/list containers to inspect leaf block
        leaf, c_depth = _strip_containers(raw_line)

        # 3. Blank line resets paragraph buffer
        if not leaf.strip():
            para_lines.clear()
            continue

        # 4. Opening HTML comment block (0-3 spaces)
        if re.match(r"^ {0,3}<!--", leaf):
            para_lines.clear()
            if "-->" not in leaf[leaf.find("<!--") + 4 :]:
                in_html_comment = True
            continue

        # 5. Opening fenced code block (0-3 spaces; backtick info string may not contain backticks)
        mf = fence_open_re.match(leaf)
        if mf:
            f_tok, info = mf.group(1), mf.group(2)
            if f_tok[0] == "~" or "`" not in info:
                in_fence = True
                fence_char = f_tok[0]
                fence_len = len(f_tok)
                para_lines.clear()
                continue

        # 6. Indented code block (>= 4 spaces or leading tab, only when NOT interrupting a paragraph)
        if not para_lines and (leaf.startswith("    ") or leaf.startswith("\t")):
            continue

        # 7. Multiline Setext underline immediately following non-empty paragraph buffer at same container depth
        ms = setext_re.match(leaf)
        if ms and para_lines and para_lines[-1][2] == c_depth:
            level = 1 if ms.group(1)[0] == "=" else 2
            title = " ".join(x[1] for x in para_lines).strip()
            start_idx = para_lines[0][0]
            headings.append({
                "kind": "setext",
                "level": level,
                "title": title,
                "rendered_title": _render_cm_inline_text(title),
                "raw_first_line": raw_lines[start_idx],
                "start_line": start_idx,
                "end_line": idx,
                "start_offset": line_offsets[start_idx],
                "end_offset": line_offsets[idx] + len(raw_line),
                "container_depth": c_depth,
            })
            para_lines.clear()
            continue

        # 8. ATX heading (1-6 '#', 0-3 leading spaces)
        ma = atx_re.match(leaf)
        if ma:
            para_lines.clear()
            level = len(ma.group(1))
            title = _clean_atx_title(ma.group(2) or "")
            headings.append({
                "kind": "atx",
                "level": level,
                "title": title,
                "rendered_title": _render_cm_inline_text(title),
                "raw_first_line": raw_line,
                "start_line": idx,
                "end_line": idx,
                "start_offset": line_offsets[idx],
                "end_offset": line_offsets[idx] + len(raw_line),
                "container_depth": c_depth,
            })
            continue

        # 9. Otherwise accumulate paragraph line (for potential multiline Setext heading)
        if para_lines and para_lines[-1][2] != c_depth:
            para_lines.clear()
        para_lines.append((idx, leaf.strip(), c_depth))

    return headings


def evaluate_r3p05_dec103_thinning(
    *,
    k: Any,
    n_raw: Any,
    eta_spike_recovery_l95: Any,
    f_clonogenic: Any = None,
    f_clonogenic_provenance: Any = None,
    alpha_obs: Any = None,
    alpha_f_clonogenic: Any = None,
    alpha_eta_recovery: Any = None,
    max_u95_fraction: Any = R3P05_MAX_U95_CEILING,
    expected_scope: Dict[str, Any] | None = None,
    observed_scope: Dict[str, Any] | None = None,
    sampling_independence_model: Any = None,
    clustered_sampling: Any = None,
    backup_gate0v_primary_failed: Any = None,
) -> Dict[str, Any]:
    """Fail-closed exact 50-digit conservative Decimal DEC-103 independent detection-thinning subgate for R3P-05.

    Enforces:
      - Unconditional k >= 1 observed-zero hard veto FIRST (with separate attribution_status field).
      - Strict hard ceiling max_u95_fraction in (0, 1e-6] evaluated in exact 50-digit Decimal.
      - Canonical exact 7-key DEC-103 ASSEMBLED_DEVICE_LOT_SCOPE (exact key set, control/zero-width/whitespace
        rejection, eligible full-stack Design Y v2 chr15 donor_sha256 binding, ROUTE_3P_MODULAR_PLANAR_ENVELOPE,
        and explicit backup_gate0v_primary_failed is True for backup CD55+CFH).
      - Strict isinstance(sampling_independence_model, str) check BEFORE equality (blocks NumPy 1-elem/2-elem arrays).
      - Mandatory explicit f_clonogenic_provenance:
          * 'KNOWN_OR_CONSERVATIVELY_ESTABLISHED': uses alpha_obs = 0.05 (unless stricter explicit alpha_obs given),
            yielding >= 90% Bonferroni joint coverage with 95% LCB recovery.
          * 'INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION': requires explicit alpha_obs and alpha_f_clonogenic
            such that conservative ROUND_CEILING sum alpha_obs + alpha_eta_recovery + alpha_f_clonogenic <= 0.10.
      - Exact 50-digit conservative Decimal upper bound (ROUND_FLOOR on ln(alpha)/n, exp, and f*eta; ROUND_CEILING on
        1 - exp and division) comparing exact integer n_raw without binary64 ULP false PASS or 0.0 underflow display.
      - Mathematical count subgate ONLY: empirical_clearance_status always reports retrieval/whole-product pending.
    """
    is_int_k = isinstance(k, int) and not isinstance(k, bool)
    is_int_n = isinstance(n_raw, int) and not isinstance(n_raw, bool)
    cap_d = _safe_exact_decimal(max_u95_fraction)
    scope_ok, scope_reason, pvid = _check_canonical_scope_integrity(expected_scope, observed_scope)

    # 1. Observed-zero hard veto latches FIRST with separate attribution status
    if is_int_k and k > 0:
        cap_out = float(cap_d) if (cap_d is not None and Decimal(0) < cap_d <= R3P05_MAX_U95_CEILING_DEC) else None
        valid_count_rel = is_int_n and n_raw > 0 and k <= n_raw
        attr_status = (
            "EXACT_7KEY_SCOPE_ATTRIBUTED"
            if (scope_ok and valid_count_rel)
            else f"UNATTRIBUTED_OR_INVALID_AUX_INPUTS_UNCONDITIONAL_STOP:{scope_reason if not scope_ok else 'INVALID_N_OR_K_GT_N'}"
        )
        return {
            "status": "FAIL_HARD_STOP",
            "reason": "OBSERVED_NONZERO_VIABLE_ESCAPE_HARD_VETO",
            "attribution_status": attr_status,
            "empirical_clearance_status": "FAIL_HARD_STOP_OBSERVED_ESCAPE",
            "k": k,
            "n_raw": n_raw if is_int_n else None,
            "cap_fraction": cap_out,
        }

    # 2. Strict hard ceiling cap validation (exact Decimal <= 1e-6)
    if cap_d is None or not (Decimal(0) < cap_d <= R3P05_MAX_U95_CEILING_DEC):
        return {
            "status": "HOLD_RELAXED_OR_INVALID_CAP_FRACTION",
            "k": k if is_int_k else None,
            "n_raw": n_raw if is_int_n else None,
            "preregistered_max_cap_fraction": R3P05_MAX_U95_CEILING,
        }
    cap_float = float(cap_d)

    # 3. Canonical exact 7-key DEC-103 ASSEMBLED_DEVICE_LOT_SCOPE + registered eligible full-stack product/digest binding
    if not scope_ok:
        st_code = (
            "HOLD_INELIGIBLE_OR_UNREGISTERED_PRODUCT_SCOPE"
            if scope_reason.startswith(("PRODUCT_NOT_ELIGIBLE", "CONSTRUCT_SHA256_MISMATCH"))
            else "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING"
        )
        return {
            "status": st_code,
            "reason": scope_reason,
            "k": k if is_int_k else None,
            "n_raw": n_raw if is_int_n else None,
            "cap_fraction": cap_float,
        }
    if pvid == "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR" and backup_gate0v_primary_failed is not True:
        return {
            "status": "HOLD_BACKUP_NOT_ACTIVATED_PRIMARY_GATE0V_NOT_FAILED",
            "reason": "BACKUP_CD55_CFH_REQUIRES_EXPLICIT_PRIMARY_CD46_WT_GATE0V_FAILURE",
            "k": k if is_int_k else None,
            "n_raw": n_raw if is_int_n else None,
            "cap_fraction": cap_float,
        }

    # 4. Strict integer & bounds validation (exact arbitrary-precision integer n_raw supported!)
    if k is None or n_raw is None:
        return {"status": "HOLD_MISSING_ATOMIC_INTEGER_K_OR_N", "cap_fraction": cap_float}
    if not is_int_k or not is_int_n:
        return {"status": "HOLD_NON_INTEGER_COUNT_OR_DENOMINATOR_REJECTED", "cap_fraction": cap_float}
    if n_raw <= 0 or k < 0 or k > n_raw:
        return {"status": "HOLD_INVALID_COUNT_BOUNDS", "cap_fraction": cap_float}

    # 5. Explicit fail-closed sampling independence guard (strict isinstance(..., str) BEFORE equality!)
    if (
        not isinstance(sampling_independence_model, str)
        or sampling_independence_model != "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"
        or clustered_sampling is not False
    ):
        return {
            "status": "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE",
            "cap_fraction": cap_float,
        }

    # 6. Explicit calibration & provenance check (NO implicit 1.0 default; mandatory f_clonogenic_provenance!)
    if eta_spike_recovery_l95 is None or f_clonogenic is None:
        return {"status": "HOLD_MISSING_RECOVERY_OR_CLONOGENIC_FACTOR", "cap_fraction": cap_float}
    eta_d = _safe_exact_decimal(eta_spike_recovery_l95)
    f_d = _safe_exact_decimal(f_clonogenic)
    if eta_d is None or f_d is None:
        return {"status": "HOLD_NONFINITE_RECOVERY_OR_CLONOGENIC_FACTOR", "cap_fraction": cap_float}
    if not (Decimal(0) < eta_d <= Decimal(1)) or not (Decimal(0) < f_d <= Decimal(1)):
        return {"status": "HOLD_OUT_OF_RANGE_RECOVERY_OR_CLONOGENIC_FACTOR", "cap_fraction": cap_float}
    if (eta_d < Decimal(1) or f_d < Decimal(1)) and eta_d < Decimal("0.30"):
        return {"status": "HOLD_SPIKE_RECOVERY_BELOW_0P30_QUALIFICATION_FLOOR", "cap_fraction": cap_float}

    if not isinstance(f_clonogenic_provenance, str) or f_clonogenic_provenance not in VALID_F_CLONOGENIC_PROVENANCE:
        return {
            "status": "HOLD_MISSING_OR_INVALID_CLONOGENIC_PROVENANCE",
            "reason": "REQUIRES_KNOWN_OR_CONSERVATIVELY_ESTABLISHED_OR_INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
            "cap_fraction": cap_float,
        }

    # Resolve Bonferroni alpha budget conservatively (joint error <= 0.10 => >= 90% simultaneous coverage).
    # Because eta_spike_recovery_l95 is an L95 (one-sided 95% LCB, alpha = 0.05) bound whenever eta < 1.0,
    # its Bonferroni error share MUST be >= 0.05 whenever eta < 1.0 (preventing under-allocated L95 recovery claims).
    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_CEILING
        if f_clonogenic_provenance == "KNOWN_OR_CONSERVATIVELY_ESTABLISHED":
            a_obs_d = Decimal("0.05") if alpha_obs is None else _safe_exact_decimal(alpha_obs)
            a_eta_d = (
                (Decimal("0.05") if eta_d < Decimal(1) else Decimal("0"))
                if alpha_eta_recovery is None
                else _safe_exact_decimal(alpha_eta_recovery)
            )
            a_f_d = Decimal("0")
            if a_obs_d is None or a_eta_d is None or not (Decimal(0) < a_obs_d <= Decimal("0.05")) or a_eta_d < Decimal(0):
                return {"status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION", "cap_fraction": cap_float}
            if eta_d < Decimal(1) and a_eta_d < Decimal("0.05"):
                return {
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "reason": "ETA_SPIKE_RECOVERY_L95_REQUIRES_ALPHA_ETA_AT_LEAST_0P05",
                    "cap_fraction": cap_float,
                }
            if (a_obs_d + a_eta_d) > Decimal("0.10"):
                return {
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "reason": "BONFERRONI_ALPHA_SUM_EXCEEDS_0P10_FOR_90PCT_JOINT_COVERAGE",
                    "cap_fraction": cap_float,
                }
        else:
            a_obs_d = _safe_exact_decimal(alpha_obs)
            a_f_d = _safe_exact_decimal(alpha_f_clonogenic)
            a_eta_d = (
                (Decimal("0.05") if eta_d < Decimal(1) else Decimal("0"))
                if alpha_eta_recovery is None
                else _safe_exact_decimal(alpha_eta_recovery)
            )
            if (
                a_obs_d is None
                or a_f_d is None
                or a_eta_d is None
                or not (Decimal(0) < a_obs_d <= Decimal("0.05"))
                or not (Decimal(0) < a_f_d < Decimal("0.10"))
                or a_eta_d < Decimal(0)
            ):
                return {"status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION", "cap_fraction": cap_float}
            if eta_d < Decimal(1) and a_eta_d < Decimal("0.05"):
                return {
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "reason": "ETA_SPIKE_RECOVERY_L95_REQUIRES_ALPHA_ETA_AT_LEAST_0P05",
                    "cap_fraction": cap_float,
                }
            if (a_obs_d + a_eta_d + a_f_d) > Decimal("0.10"):
                return {
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "reason": "BONFERRONI_ALPHA_SUM_EXCEEDS_0P10_FOR_90PCT_JOINT_COVERAGE",
                    "cap_fraction": cap_float,
                }

    # 7. Exact 80-digit strictly sound conservative Decimal upper bound:
    # Decimal.ln() and Decimal.exp() always round half-even (< 1 ULP error at prec=80), so applying .next_minus()
    # to BOTH a_obs_d.ln() and q_lo.exp() strictly encloses the true real transcendental values from below,
    # guaranteeing u95_upper > exact_real_U95 for EVERY input (immune to binary64 and high-precision straddles!).
    with localcontext() as ctx:
        ctx.prec = 80
        ctx.rounding = ROUND_FLOOR
        eff_lower = f_d * eta_d
        if eff_lower <= Decimal(0):
            return {"status": "HOLD_UNDERFLOW_OR_INVALID_DETECTION_FACTOR", "cap_fraction": cap_float}
        ln_alpha_lo = a_obs_d.ln().next_minus()  # strictly < exact ln(alpha_obs) < 0
        n_d = Decimal(n_raw)
        q_lo = ln_alpha_lo / n_d  # ROUND_FLOOR on negative numerator -> strictly < exact ln(alpha_obs) / n_raw
        if q_lo > Decimal("-1e-60"):
            # Since 1 - exp(q) <= -q for all q <= 0, -q_lo is a strict conservative positive upper bound even at n=10^308!
            ctx.rounding = ROUND_CEILING
            p_raw_upper = -q_lo
        else:
            exp_lo = q_lo.exp().next_minus()  # strictly < exact exp(q_lo) < exact alpha_obs**(1/n_raw)
            ctx.rounding = ROUND_CEILING
            p_raw_upper = Decimal(1) - exp_lo
        ctx.rounding = ROUND_CEILING
        u95_upper = p_raw_upper / eff_lower
        if u95_upper > Decimal(1):
            u95_upper = Decimal(1)

    u95_str = format(u95_upper, ".16e")
    eff_str = format(eff_lower, ".12g")
    u95_display = float(u95_upper) if u95_upper >= Decimal("1e-300") else None

    if u95_upper > cap_d:
        return {
            "status": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_AND_EMPIRICAL_RETRIEVAL_PENDING",
            "k": 0,
            "n_raw": n_raw,
            "f_clonogenic_provenance": f_clonogenic_provenance,
            "alpha_obs_str": format(a_obs_d, ".6g"),
            "effective_detection_factor_str": eff_str,
            "u95_exact_dec103_str": u95_str,
            "u95_exact_dec103_display": round(u95_display, 16) if u95_display is not None else None,
            "cap_fraction_str": format(cap_d, ".16e"),
            "cap_fraction": cap_float,
        }
    return {
        "status": "PASS",
        "empirical_clearance_status": "MATHEMATICAL_COUNT_SUBGATE_ONLY_EMPIRICAL_RETRIEVAL_AND_WHOLE_PRODUCT_PENDING",
        "k": 0,
        "n_raw": n_raw,
        "f_clonogenic_provenance": f_clonogenic_provenance,
        "alpha_obs_str": format(a_obs_d, ".6g"),
        "effective_detection_factor_str": eff_str,
        "u95_exact_dec103_str": u95_str,
        "u95_exact_dec103_display": round(u95_display, 16) if u95_display is not None else None,
        "cap_fraction_str": format(cap_d, ".16e"),
        "cap_fraction": cap_float,
    }


def compute_r3p05_dec103_thinning_boundary_audit() -> Dict[str, Any]:
    """Deterministic fail-closed R3P-05 DEC-103 exact 50-digit Decimal boundary & single-dose feasibility audit."""
    valid_scope = {
        "product_variant_id": "ISLET1_Y_FULLSTACK_V2_B2M_11TU",
        "construct_sha256": "e77a358ca5bc55c9e0a3070ab46b115c3314660d3a9be6221ef9e91d698cab4e",
        "donor_seedline_id": "DONOR_SEED_A",
        "mcb_id": "MCB_01",
        "lot_id": "LOT_A",
        "device_assembly_id": "DEV_A",
        "physical_route_id": "ROUTE_3P_MODULAR_PLANAR_ENVELOPE",
    }
    valid_ctx = {
        "f_clonogenic": 1.0,
        "f_clonogenic_provenance": "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
        "expected_scope": valid_scope,
        "observed_scope": dict(valid_scope),
        "sampling_independence_model": "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES",
        "clustered_sampling": False,
    }
    lot_bad = dict(valid_scope, lot_id="LOT_B")
    dev_bad = dict(valid_scope, device_assembly_id="DEV_B")
    padded_lot = dict(valid_scope, lot_id=" LOT_A ")
    ctrl_char_lot = dict(valid_scope, lot_id="LOT\x00A")
    zwsp_lot = dict(valid_scope, lot_id="LOT\u200bA")
    ctrl_scope = dict(
        valid_scope,
        product_variant_id="ISLET1_Y_FULLSTACK_V2_CTRL_EQUAL_SPACING",
        construct_sha256="d5f69bd613cc0d598b846ba061154dc522b8ec6fb064610cfad01562505d2925",
    )
    backup_scope = dict(
        valid_scope,
        product_variant_id="ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR",
        construct_sha256="d7291ca6e610f8251ab04bbf6fc62d7747e7baaf0e680300f46b895920aab9d8",
    )
    extra_key_scope_a = dict(valid_scope, extra_unverified_tag="A")
    extra_key_scope_b = dict(valid_scope, extra_unverified_tag="B")
    nondigest_scope = dict(valid_scope, construct_sha256="not-a-sha256")
    partial_2key = {"lot_id": "LOT_A", "device_assembly_id": "DEV_A"}

    cases = [
        # Single-dose Blueprint feasibility counterexamples (93,350 clonogenic opportunities -> ~32.1x / 64.2x above 1e-6!)
        ("blueprint_single_dose_93350_eta1p0_HOLD_unachievable_in_single_dose", evaluate_r3p05_dec103_thinning(k=0, n_raw=93_350, eta_spike_recovery_l95=1.0, **valid_ctx)),
        ("blueprint_single_dose_93350_eta0p5_HOLD_unachievable_in_single_dose", evaluate_r3p05_dec103_thinning(k=0, n_raw=93_350, eta_spike_recovery_l95=0.50, **valid_ctx)),
        ("blueprint_single_dose_933p5M_raw_f1e4_eta1p0_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=933_500_000, eta_spike_recovery_l95=1.0, f_clonogenic=1e-4, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        # Standard exact boundaries
        ("eta1p0_n2995730_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=2_995_730, eta_spike_recovery_l95=1.0, **valid_ctx)),
        ("eta1p0_n2995731_PASS", evaluate_r3p05_dec103_thinning(k=0, n_raw=2_995_731, eta_spike_recovery_l95=1.0, **valid_ctx)),
        ("eta0p5_n5991462_HOLD_rejects_floor_approximation", evaluate_r3p05_dec103_thinning(k=0, n_raw=5_991_462, eta_spike_recovery_l95=0.50, **valid_ctx)),
        ("eta0p5_n5991463_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=5_991_463, eta_spike_recovery_l95=0.50, **valid_ctx)),
        ("eta0p5_n5991464_PASS", evaluate_r3p05_dec103_thinning(k=0, n_raw=5_991_464, eta_spike_recovery_l95=0.50, **valid_ctx)),
        # Binary64 ULP tie counterexample & huge-n non-zero scientific display
        ("binary64_ulp_straddle_n5991464_cap4p999999206572849e7_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=5_991_464, eta_spike_recovery_l95=1.0, max_u95_fraction=4.999999206572849e-7, **valid_ctx)),
        ("huge_n_1e308_nonzero_sci_display_PASS", evaluate_r3p05_dec103_thinning(k=0, n_raw=10**308, eta_spike_recovery_l95=1.0, **valid_ctx)),
        # Provenance & 3-way Bonferroni alpha allocation
        ("missing_f_clonogenic_provenance_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("estimated_f_clono_without_alpha_allocation_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=6_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION", expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("estimated_f_clono_alpha_sum_exceeds_0p10_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=8_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION", alpha_obs=0.04, alpha_f_clonogenic=0.04, expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("estimated_f_clono_underallocated_l95_eta_0p01_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=10_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION", alpha_obs=Decimal("0.01"), alpha_eta_recovery=Decimal("0.01"), alpha_f_clonogenic=Decimal("0.08"), expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("estimated_f_clono_valid_3way_bonferroni_PASS", evaluate_r3p05_dec103_thinning(k=0, n_raw=7_377_760, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION", alpha_obs=Decimal("0.025"), alpha_f_clonogenic=Decimal("0.025"), alpha_eta_recovery=Decimal("0.05"), expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        # Observed-zero hard veto & attribution status
        ("k1_n5000000_FAIL_observed_zero_veto", evaluate_r3p05_dec103_thinning(k=1, n_raw=5_000_000, eta_spike_recovery_l95=1.0, **valid_ctx)),
        ("k1_overflow_cap_1e400_FAIL_HARD_STOP", evaluate_r3p05_dec103_thinning(k=1, n_raw=5_000_000, eta_spike_recovery_l95=1.0, max_u95_fraction=10**400, **valid_ctx)),
        ("k0_overflow_cap_1e400_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, max_u95_fraction=10**400, **valid_ctx)),
        ("relaxed_cap_0p01_n600_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=600, eta_spike_recovery_l95=0.50, max_u95_fraction=0.01, **valid_ctx)),
        ("nan_cap_fraction_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=1, eta_spike_recovery_l95=1.0, max_u95_fraction=float("nan"), **valid_ctx)),
        # Scope, backup, identifier hygiene & NumPy array resistance
        ("excluded_control_scope_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=ctrl_scope, observed_scope=dict(ctrl_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("inactive_backup_without_gate0v_failure_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=backup_scope, observed_scope=dict(backup_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("activated_backup_with_gate0v_failure_CONDITIONAL_ARITHMETIC_PASS", evaluate_r3p05_dec103_thinning(k=0, n_raw=2_995_731, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=backup_scope, observed_scope=dict(backup_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False, backup_gate0v_primary_failed=True)),
        ("extra_scope_key_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=extra_key_scope_a, observed_scope=extra_key_scope_b, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("non_digest_construct_sha256_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=nondigest_scope, observed_scope=dict(nondigest_scope), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("whitespace_padded_lot_id_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=padded_lot, observed_scope=dict(padded_lot), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("control_char_lot_id_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=ctrl_char_lot, observed_scope=dict(ctrl_char_lot), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("zero_width_char_lot_id_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=zwsp_lot, observed_scope=dict(zwsp_lot), sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("numpy_1elem_sampling_model_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model=np.array(["INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"]), clustered_sampling=False)),
        ("numpy_2elem_sampling_model_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=dict(valid_scope), sampling_independence_model=np.array(["INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"]), clustered_sampling=False)),
        ("omitted_scope_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("partial_2key_scope_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=partial_2key, observed_scope=partial_2key, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("omitted_independence_model_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=valid_scope, clustered_sampling=False)),
        ("null_clustered_sampling_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=valid_scope, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=None)),
        ("clustered_sampling_true_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=6_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=valid_scope, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=True)),
        ("omitted_f_clonogenic_HOLD", evaluate_r3p05_dec103_thinning(k=0, n_raw=3_000_000, eta_spike_recovery_l95=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=valid_scope, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("lot_mismatch_HOLD_gy018_scope", evaluate_r3p05_dec103_thinning(k=0, n_raw=6_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=lot_bad, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
        ("device_mismatch_HOLD_gy018_scope", evaluate_r3p05_dec103_thinning(k=0, n_raw=6_000_000, eta_spike_recovery_l95=0.50, f_clonogenic=1.0, f_clonogenic_provenance="KNOWN_OR_CONSERVATIVELY_ESTABLISHED", expected_scope=valid_scope, observed_scope=dev_bad, sampling_independence_model="INDEPENDENT_BIOLOGICAL_OPPORTUNITIES", clustered_sampling=False)),
    ]
    expected_statuses = {
        "blueprint_single_dose_93350_eta1p0_HOLD_unachievable_in_single_dose": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "blueprint_single_dose_93350_eta0p5_HOLD_unachievable_in_single_dose": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "blueprint_single_dose_933p5M_raw_f1e4_eta1p0_HOLD": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "eta1p0_n2995730_HOLD": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "eta1p0_n2995731_PASS": "PASS",
        "eta0p5_n5991462_HOLD_rejects_floor_approximation": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "eta0p5_n5991463_HOLD": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "eta0p5_n5991464_PASS": "PASS",
        "binary64_ulp_straddle_n5991464_cap4p999999206572849e7_HOLD": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        "huge_n_1e308_nonzero_sci_display_PASS": "PASS",
        "missing_f_clonogenic_provenance_HOLD": "HOLD_MISSING_OR_INVALID_CLONOGENIC_PROVENANCE",
        "estimated_f_clono_without_alpha_allocation_HOLD": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
        "estimated_f_clono_alpha_sum_exceeds_0p10_HOLD": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
        "estimated_f_clono_underallocated_l95_eta_0p01_HOLD": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
        "estimated_f_clono_valid_3way_bonferroni_PASS": "PASS",
        "k1_n5000000_FAIL_observed_zero_veto": "FAIL_HARD_STOP",
        "k1_overflow_cap_1e400_FAIL_HARD_STOP": "FAIL_HARD_STOP",
        "k0_overflow_cap_1e400_HOLD": "HOLD_RELAXED_OR_INVALID_CAP_FRACTION",
        "relaxed_cap_0p01_n600_HOLD": "HOLD_RELAXED_OR_INVALID_CAP_FRACTION",
        "nan_cap_fraction_HOLD": "HOLD_RELAXED_OR_INVALID_CAP_FRACTION",
        "excluded_control_scope_HOLD": "HOLD_INELIGIBLE_OR_UNREGISTERED_PRODUCT_SCOPE",
        "inactive_backup_without_gate0v_failure_HOLD": "HOLD_BACKUP_NOT_ACTIVATED_PRIMARY_GATE0V_NOT_FAILED",
        "activated_backup_with_gate0v_failure_CONDITIONAL_ARITHMETIC_PASS": "PASS",
        "extra_scope_key_HOLD": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "non_digest_construct_sha256_HOLD": "HOLD_INELIGIBLE_OR_UNREGISTERED_PRODUCT_SCOPE",
        "whitespace_padded_lot_id_HOLD": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "control_char_lot_id_HOLD": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "zero_width_char_lot_id_HOLD": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "numpy_1elem_sampling_model_HOLD": "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE",
        "numpy_2elem_sampling_model_HOLD": "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE",
        "omitted_scope_HOLD": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "partial_2key_scope_HOLD": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "omitted_independence_model_HOLD": "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE",
        "null_clustered_sampling_HOLD": "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE",
        "clustered_sampling_true_HOLD": "HOLD_MISSING_OR_CLUSTERED_INDEPENDENCE_EVIDENCE",
        "omitted_f_clonogenic_HOLD": "HOLD_MISSING_RECOVERY_OR_CLONOGENIC_FACTOR",
        "lot_mismatch_HOLD_gy018_scope": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
        "device_mismatch_HOLD_gy018_scope": "HOLD_SCOPE_MISMATCH_NO_CROSS_LOT_OR_DEVICE_POOLING",
    }
    by_name = {name: res for name, res in cases}
    for cname, exp_st in expected_statuses.items():
        got_st = by_name[cname]["status"]
        if got_st != exp_st:
            raise ValueError(f"R3P-05 audit case {cname} returned {got_st!r}, expected {exp_st!r}")
    if by_name["huge_n_1e308_nonzero_sci_display_PASS"]["u95_exact_dec103_str"].startswith("0.0"):
        raise ValueError("Huge-n exact bound must not print as zero")

    return {
        "formula": "k == 0 AND exact_80digit_strictly_sound_Decimal(min(1, (1 - exp_lo(ln_lo(alpha_obs)/n_raw)) / (f_clonogenic * eta_spike_recovery_l95))) <= cap <= 1e-6",
        "arithmetic_precision": "decimal.Decimal prec=80 strictly sound directed bounds (.next_minus() on ln() and exp(), ROUND_FLOOR on eff_prob & q_lo, ROUND_CEILING on 1-exp and division)",
        "required_scope_keys": list(ASSEMBLED_DEVICE_LOT_SCOPE_KEYS),
        "eligible_fullstack_product_ids": sorted(
            pid for pid, spec in pc.DESIGN_Y_V2_DONOR_BINDINGS.items()
            if str(spec.get("role") or "") in ELIGIBLE_R3P05_FULLSTACK_ROLES
        ),
        "preregistered_hard_cap_ceiling": R3P05_MAX_U95_CEILING,
        "min_n_raw_at_eta_1p0": 2_995_731,
        "min_n_raw_at_eta_0p5": 5_991_464,
        "single_dose_blueprint_feasibility_counterexample": {
            "reference_dose_total_living_cells": 933_500_000,
            "blueprint_residual_clonogenic_fraction_upper_assumption": 1e-4,
            "max_clonogenic_opportunities_per_dose": 93_350,
            "u95_at_eta_1p0_str": by_name["blueprint_single_dose_93350_eta1p0_HOLD_unachievable_in_single_dose"]["u95_exact_dec103_str"],
            "u95_at_eta_0p5_str": by_name["blueprint_single_dose_93350_eta0p5_HOLD_unachievable_in_single_dose"]["u95_exact_dec103_str"],
            "ratio_above_1e6_ceiling_eta_1p0": 32.09,
            "ratio_above_1e6_ceiling_eta_0p5": 64.18,
            "disposition": "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
            "interpretation_note": (
                "Under the Blueprint upper inventory of 93,350 clonogenic opportunities per reference dose (80,000 endocrine "
                "+ 13,350 vascular at f_prolif <= 1e-4), observing k == 0 in one assembled dose gives U95 = 3.20909e-5 at "
                "eta=1.0 (32.09x above 1e-6) and 6.41818e-5 at eta=0.50 (64.18x above 1e-6), so single-dose quantitative "
                "R3P-05 clearance holds at HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR; pooling across lots or "
                "assemblies remains prohibited under GY-018, and whole-product empirical retrieval/hemostasis clearance stays separate."
            ),
        },
        "floor_approximation_n_eff_forbidden": True,
        "explicit_7key_scope_independence_and_calibration_provenance_required": True,
        "joint_coverage_note": (
            "Combining a 95% LCB recovery factor with a 95% Clopper-Pearson count inversion yields >= 90% Bonferroni "
            "joint coverage ONLY when f_clonogenic_provenance == 'KNOWN_OR_CONSERVATIVELY_ESTABLISHED'; when "
            "f_clonogenic is independently estimated ('INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION'), explicit "
            "alpha_obs + alpha_eta_recovery + alpha_f_clonogenic <= 0.10 is enforced conservatively."
        ),
        "cases": [{"case": name, **res} for name, res in cases],
    }


def compute_planar_mechanics_table() -> List[Dict[str, Any]]:
    """Dimensional geometry, single-ribbon uniform 2-kPa pull-out, critical adhesion, and Poiseuille shear."""
    mu_pa_s = 0.0012  # cell-suspension vehicle viscosity (1.2 mPa*s)
    tau_adh_pa = 2000.0  # 2.0 kPa uniform adhesion assumption
    configs = [
        ("ROUTE_3P_SUBARS_300_H180", 300.0, 180.0, 6, 2, 3.125, 6.0),
        ("ROUTE_3P_SUBARS_300_H200", 300.0, 200.0, 6, 2, 3.125, 6.0),
        ("ROUTE_3P_HYBRID_333_H180", 333.3, 180.0, 6, 2, 3.472, 8.0),
        ("ROUTE_3P_HYBRID_333_H200", 333.3, 200.0, 6, 2, 3.472, 8.0),
        ("ROUTE_3P_4MOD_SUBARS_260_H180", 260.0, 180.0, 4, 2, 3.250, 6.0),
        ("NAIVE_147CM2_ANOXIC_BARRED", 147.0, 180.0, 4, 2, 2.450, 6.0),
    ]
    rows: List[Dict[str, Any]] = []
    for name, area_cm2, h_um, n_env, ribbons_per_env, width_cm, prescribed_load_min in configs:
        n_ribbons = n_env * ribbons_per_env
        areal_ieq_cm2 = DOSE_IEQ / area_cm2
        vol_ml = area_cm2 * (h_um * 1e-4)
        ribbon_one_sided_m2 = (area_cm2 / n_ribbons) * 1e-4
        pullout_n = 2.0 * ribbon_one_sided_m2 * tau_adh_pa
        crit_adhesion_kpa = (15.0 / (2.0 * ribbon_one_sided_m2)) / 1000.0
        q_m3_s = (vol_ml * 1e-6) / (prescribed_load_min * 60.0)
        w_m = width_cm * 1e-2
        h_m = h_um * 1e-6
        tau_w_pa = 6.0 * mu_pa_s * q_m3_s / (w_m * (h_m ** 2))
        rows.append({
            "configuration": name,
            "active_one_sided_area_cm2": area_cm2,
            "dual_sided_area_cm2": round(2.0 * area_cm2, 1),
            "slot_height_um": h_um,
            "envelopes": n_env,
            "ribbons_total": n_ribbons,
            "areal_loading_ieq_cm2": round(areal_ieq_cm2, 2),
            "ideal_slot_volume_ml": round(vol_ml, 3),
            "prescribed_loading_duration_input_min": prescribed_load_min,
            "pullout_unit": "single_individually_extracted_ribbon",
            "uniform_2kpa_single_ribbon_pullout_n": round(pullout_n, 3),
            "critical_uniform_adhesion_at_15N_kpa": round(crit_adhesion_kpa, 3),
            "two_ribbon_concurrent_2kpa_pullout_warning_n": round(2.0 * pullout_n, 3),
            "poiseuille_loading_wall_shear_pa": round(tau_w_pa, 4),
            "vpio_120um_lateral_width_ratio": round((width_cm * 1e4) / 120.0, 1),
            "vpio_120um_slot_height_ratio": round(h_um / 120.0, 3),
            "height_ratio_caveat": (
                "Slot height is 1.50-1.67x 120 um VPIO diameter; loading duration is a prescribed input to the "
                "steady Poiseuille shear estimate, and anatomical fit / peak entrance shear / atraumatic hemostatic "
                "extraction remain pending Stage-0 R3P-01/R3P-02."
            ),
        })
    return rows


def solve_planar_route3p_fvm(
    *,
    label: str,
    area_cm2: float,
    h_slot_um: float = 180.0,
    cluster_diam_um: float = 65.0,
    bed_po2_mmhg: float = 35.0,
    sleeve_um: float = 0.0,
    sleeve_d_factor: float = 1.0,
    k_perf_gel_s_inv: float = 0.0,
    k_perf_islet_s_inv: float = 0.0,
    ocr_mult: float = 1.097,
    alpha_mol_cm3_mmhg: float = ALPHA_LOW_SALT,
    nx: int = 28,
    nz: int = 32,
) -> Dict[str, Any]:
    """Deterministic 3D FVM of a periodic planar Route 3P envelope unit cell using live Layer-4 solver."""
    for pname, pval in (
        ("area_cm2", area_cm2),
        ("h_slot_um", h_slot_um),
        ("cluster_diam_um", cluster_diam_um),
        ("bed_po2_mmhg", bed_po2_mmhg),
        ("sleeve_d_factor", sleeve_d_factor),
        ("ocr_mult", ocr_mult),
        ("alpha_mol_cm3_mmhg", alpha_mol_cm3_mmhg),
    ):
        if isinstance(pval, bool) or not isinstance(pval, (int, float)) or not math.isfinite(float(pval)) or float(pval) <= 0.0:
            raise ValueError(f"Invalid positive FVM parameter {pname}={pval!r} in {label}")
    for pname, pval in (("sleeve_um", sleeve_um), ("k_perf_gel_s_inv", k_perf_gel_s_inv), ("k_perf_islet_s_inv", k_perf_islet_s_inv)):
        if isinstance(pval, bool) or not isinstance(pval, (int, float)) or not math.isfinite(float(pval)) or float(pval) < 0.0:
            raise ValueError(f"Invalid non-negative FVM parameter {pname}={pval!r} in {label}")
    if cluster_diam_um >= h_slot_um:
        raise ValueError(f"Invalid geometry in {label}: cluster_diam_um ({cluster_diam_um}) must be < h_slot_um ({h_slot_um})")
    if isinstance(nx, bool) or isinstance(nz, bool) or not isinstance(nx, int) or not isinstance(nz, int) or nx < 4 or nz < 4:
        raise ValueError(f"Invalid FVM grid resolution in {label}: nx={nx}, nz={nz} (minimum 4x4x4 required)")
    if abs(alpha_mol_cm3_mmhg - ALPHA_LOW_SALT) > 1e-18 and (k_perf_gel_s_inv > 0.0 or k_perf_islet_s_inv > 0.0):
        raise ValueError(
            f"Solubility override alpha={alpha_mol_cm3_mmhg} cannot be combined with non-zero perfusion in {label} "
            "because the underlying Layer-4 perfusion term uses the global reference solubility."
        )

    loading_ieq_cm2 = DOSE_IEQ / area_cm2
    r_cm = 0.5 * cluster_diam_um * 1e-4
    v_cl = (4.0 / 3.0) * math.pi * (r_cm ** 3)
    bilayer = cluster_diam_um <= 80.0
    n_layer = (0.5 if bilayer else 1.0) * (loading_ieq_cm2 * VOL_IEQ / v_cl)
    pitch_cm = math.sqrt(1.0 / n_layer)
    h_cm = h_slot_um * 1e-4
    dx = pitch_cm / nx
    dz = h_cm / nz
    xs = (np.arange(nx) + 0.5) * dx
    zs = (np.arange(nz) + 0.5) * dz - 0.5 * h_cm
    X, Y, Z = np.meshgrid(xs, xs, zs, indexing="ij")
    dv = dx * dx * dz

    def pdist(a: np.ndarray, c: float) -> np.ndarray:
        d = np.abs(a - c)
        return np.minimum(d, pitch_cm - d)

    in_any = np.zeros_like(X, dtype=bool)
    vmax = np.zeros_like(X, dtype=np.float64)
    masks: List[np.ndarray] = []
    if bilayer:
        z_c = 0.5 * h_cm - r_cm - 15.0e-4
        if z_c < 0.0:
            raise ValueError(f"Bilayer standoff geometry invalid in {label}: z_c={z_c}")
        centers = [(0.5 * pitch_cm, 0.5 * pitch_cm, +z_c), (0.0, 0.0, -z_c)]
    else:
        centers = [(0.5 * pitch_cm, 0.5 * pitch_cm, 0.0)]

    for cx, cy, cz in centers:
        msk = (pdist(X, cx) ** 2 + pdist(Y, cy) ** 2 + (Z - cz) ** 2) <= (r_cm ** 2)
        v_vox = float(np.sum(msk)) * dv
        if v_vox <= 0.0:
            raise ValueError(f"Cluster voxel mask empty at grid ({nx},{nx},{nz}) in {label}")
        in_any |= msk
        vmax[msk] = (VMAX_PURE * ocr_mult) * (v_cl / v_vox)
        masks.append(msk)

    in_gel = ~in_any
    gel_frac = float(np.mean(in_gel))
    if gel_frac <= 0.0:
        raise ValueError(f"Gel domain empty in {label}")
    phi_tot = (loading_ieq_cm2 * VOL_IEQ) / h_cm
    vmax[in_gel] = (phi_tot * VMAX_PURE * 0.15 * ocr_mult) / gel_frac
    perm = np.where(in_any, D_ISLET * alpha_mol_cm3_mmhg, D_GEL_CHRON * alpha_mol_cm3_mmhg)
    kp = np.where(in_any, k_perf_islet_s_inv, k_perf_gel_s_inv)

    res = l4._solve_3d_oxygen_strict_iip_0019(
        perm, vmax, np.zeros_like(X), kp, dx, dx, dz,
        p_bed=bed_po2_mmhg, standoff_z_cm=sleeve_um * 1e-4, p_ext=D_GEL_CHRON * sleeve_d_factor * alpha_mol_cm3_mmhg,
    )
    p3 = res["p_3d"]
    worst_face = min(
        float(l4._compute_exact_fv_interface_stats_iip_0019(p3, perm, m)["fv_harmonic_face_min_mmhg"])
        for m in masks
    )
    worst_core = min(float(np.min(p3[m])) for m in masks)
    raw_mbe = float(res["rel_mass_balance_error"])
    if not math.isfinite(raw_mbe) or raw_mbe >= 1e-10:
        raise ValueError(f"Mass balance residual exceeded 1e-10 in {label}: {raw_mbe}")
    return {
        "label": label,
        "grid_nx_ny_nz": [nx, nx, nz],
        "area_cm2": area_cm2,
        "loading_ieq_cm2": round(loading_ieq_cm2, 2),
        "slot_height_um": h_slot_um,
        "cluster_diam_um": cluster_diam_um,
        "bed_po2_mmhg": bed_po2_mmhg,
        "sleeve_um": sleeve_um,
        "sleeve_d_factor": round(sleeve_d_factor, 4),
        "k_perf_gel_s_inv": k_perf_gel_s_inv,
        "k_perf_islet_s_inv": k_perf_islet_s_inv,
        "ocr_mult": round(ocr_mult, 4),
        "alpha_mol_cm3_mmhg": alpha_mol_cm3_mmhg,
        "pitch_um": round(pitch_cm * 1e4, 3),
        "worst_face_mmhg": round(worst_face, 4),
        "worst_core_mmhg": round(worst_core, 4),
        "passes_15_and_5": bool(worst_face >= 15.0 and worst_core >= 5.0),
        "passes_core_5mmhg_floor": bool(worst_core >= 5.0),
        "rel_mass_balance_error_str": f"{raw_mbe:.3e}",
        "mass_balance_bound_lt_1e10": True,
    }


def _solve_1d_slab_fall_t50_min(modes: List[tuple[float, float]], k_int: float) -> float:
    def _frac(t_s: float) -> float:
        tot = 0.0
        for cn, lam in modes:
            if abs(lam - k_int) < 1e-12:
                term = (1.0 + k_int * t_s) * math.exp(-k_int * t_s)
            else:
                term = (lam * math.exp(-k_int * t_s) - k_int * math.exp(-lam * t_s)) / (lam - k_int)
            tot += cn * term
        return tot

    lo, hi = 0.0, 1200.0
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        if _frac(mid) > 0.5:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) / 60.0


def _uniform_modes(l_eff_cm: float, d_ins_cm2_s: float, n_modes: int) -> List[tuple[float, float]]:
    modes: List[tuple[float, float]] = []
    pi_2l = math.pi / (2.0 * l_eff_cm)
    for n in range(n_modes):
        odd = 2 * n + 1
        a = odd * pi_2l
        modes.append((8.0 / ((odd * math.pi) ** 2), a * a * d_ins_cm2_s))
    return modes


def _spatial_modes(l_cm: float, r_cm: float, z_c: float, d_ins_cm2_s: float, n_modes: int) -> List[tuple[float, float]]:
    modes: List[tuple[float, float]] = []
    pi_2l = math.pi / (2.0 * l_cm)
    for n in range(n_modes):
        odd = 2 * n + 1
        a = odd * pi_2l
        x = a * r_cm
        sph_factor = 3.0 * (math.sin(x) - x * math.cos(x)) / (x ** 3)
        cn = (4.0 / (odd * math.pi)) * ((-1.0) ** n) * math.cos(a * z_c) * sph_factor
        modes.append((cn, a * a * d_ins_cm2_s))
    return modes


def _bisect_intrinsic_t50_for_target(modes: List[tuple[float, float]], target_fall_t50_min: float = 3.0) -> float:
    lo, hi = 0.5, 3.0
    for _ in range(52):
        mid = 0.5 * (lo + hi)
        k_int = math.log(2.0) / (mid * 60.0)
        if _solve_1d_slab_fall_t50_min(modes, k_int) < target_fall_t50_min:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def compute_planar_insulin_shutoff_table() -> Dict[str, Any]:
    """Exact 1D bilateral slab Fourier-series insulin shutoff reporting BOTH 45-mode & converged 5,000-mode reference."""
    d_ins_cm2_s = 5.0e-7  # hydrated hydrogel insulin diffusivity (~1/3 free-water 1.5e-6 cm^2/s)
    t50_intrinsic_min = 1.75  # hardcoded proposal modeling assumption (NOT a reported Nair 2019 fitted half-time)
    k_int = math.log(2.0) / (t50_intrinsic_min * 60.0)
    n_modes_ref = 5000
    n_modes_45 = 45
    cases = [
        ("H180_clean", 180.0, 0.0, 1.0),
        ("H180_plus_15um_fibrin", 180.0, 15.0, 0.5),
        ("H200_clean", 200.0, 0.0, 1.0),
        ("H220_clean_BARRED", 220.0, 0.0, 1.0),
    ]
    rows: List[Dict[str, Any]] = []
    h200_unif_modes_45: List[tuple[float, float]] = []
    h200_unif_modes_ref: List[tuple[float, float]] = []
    for lbl, h_um, df_um, d_fib_factor in cases:
        l_eff_cm = (0.5 * h_um + (df_um / d_fib_factor)) * 1e-4
        m45 = _uniform_modes(l_eff_cm, d_ins_cm2_s, n_modes_45)
        mref = _uniform_modes(l_eff_cm, d_ins_cm2_s, n_modes_ref)
        if lbl == "H200_clean":
            h200_unif_modes_45 = m45
            h200_unif_modes_ref = mref
        t50_45 = _solve_1d_slab_fall_t50_min(m45, k_int)
        t50_ref = _solve_1d_slab_fall_t50_min(mref, k_int)
        rows.append({
            "case": lbl,
            "slot_height_um": h_um,
            "fibrin_sleeve_um": df_um,
            "effective_half_thickness_um": round(l_eff_cm * 1e4, 2),
            "modeled_fall_t50_converged_5000mode_min": round(t50_ref, 5),
            "modeled_fall_t50_45mode_min": round(t50_45, 5),
            "modeled_fall_t50_min": round(t50_ref, 3),
            "passes_lt_3p0_min": bool(t50_ref < 3.0),
        })

    spatial_rows: List[Dict[str, Any]] = []
    h200_sp_modes_ref: List[tuple[float, float]] = []
    h200_sp_t50_ref = 0.0
    for lbl, h_um, diam_um in [
        ("Spatial_65um_Bilayer_H180_clean", 180.0, 65.0),
        ("Spatial_65um_Bilayer_H200_clean", 200.0, 65.0),
        ("Spatial_120um_CenteredVPIO_H180_clean_FAIL", 180.0, 120.0),
        ("Spatial_120um_CenteredVPIO_H200_clean_FAIL", 200.0, 120.0),
    ]:
        l_cm = 0.5 * h_um * 1e-4
        r_cm = 0.5 * diam_um * 1e-4
        z_c = (l_cm - r_cm - 15.0e-4) if diam_um <= 80.0 else 0.0
        m45 = _spatial_modes(l_cm, r_cm, z_c, d_ins_cm2_s, n_modes_45)
        mref = _spatial_modes(l_cm, r_cm, z_c, d_ins_cm2_s, n_modes_ref)
        t50_45 = _solve_1d_slab_fall_t50_min(m45, k_int)
        t50_ref = _solve_1d_slab_fall_t50_min(mref, k_int)
        if lbl == "Spatial_65um_Bilayer_H200_clean":
            h200_sp_modes_ref = mref
            h200_sp_t50_ref = t50_ref
        spatial_rows.append({
            "case": lbl,
            "slot_height_um": h_um,
            "cluster_diam_um": diam_um,
            "sphere_center_offset_um": round(z_c * 1e4, 2),
            "modeled_fall_t50_converged_5000mode_min": round(t50_ref, 6),
            "modeled_fall_t50_45mode_min": round(t50_45, 6),
            "modeled_fall_t50_min": round(t50_ref, 3),
            "passes_lt_3p0_min": bool(t50_ref < 3.0),
        })

    cutoff_unif_45 = _bisect_intrinsic_t50_for_target(h200_unif_modes_45, 3.0)
    cutoff_unif_ref = _bisect_intrinsic_t50_for_target(h200_unif_modes_ref, 3.0)
    cutoff_sp_ref = _bisect_intrinsic_t50_for_target(h200_sp_modes_ref, 3.0)
    margin_sp_h200_s = (3.0 - h200_sp_t50_ref) * 60.0

    return {
        "model_assumptions": {
            "d_insulin_hydrogel_cm2_s": d_ins_cm2_s,
            "unencapsulated_intrinsic_t50_min": t50_intrinsic_min,
            "intrinsic_t50_parameter_provenance": "proposal_modeling_assumption_not_reported_nair_fit",
            "reference_fourier_modes": n_modes_ref,
            "comparison_truncated_fourier_modes": n_modes_45,
            "h200_uniform_clean_cutoff_intrinsic_t50_5000mode_min": round(cutoff_unif_ref, 6),
            "h200_uniform_clean_cutoff_intrinsic_t50_45mode_min": round(cutoff_unif_45, 6),
            "h200_spatial_65um_clean_cutoff_intrinsic_t50_5000mode_min": round(cutoff_sp_ref, 6),
            "h200_spatial_65um_clean_margin_below_3min_seconds": round(margin_sp_h200_s, 3),
            "fibrin_sleeve_model": "effective_homogeneous_transit_thickness_approximation",
            "caveat": (
                "Continuum 1D bilateral slab Fourier model using a 1.75-min intrinsic shutoff modeling assumption "
                "(motivated by Nair 2019 S3d step-down context, NOT a directly reported Nair fitted half-time). "
                f"Converged {n_modes_ref:,}-mode reference shifts uniform slabs by +0.008-0.009 min vs {n_modes_45} modes "
                f"(uniform H200 3.0-min intrinsic cutoff: {cutoff_unif_ref:.6f} min at {n_modes_ref:,} modes "
                f"[1.808982 min at 30,000 modes] vs {cutoff_unif_45:.6f} min at {n_modes_45} modes) and spatial projections "
                f"by <= 0.000014 min (spatial 65 um H200 cutoff: {cutoff_sp_ref:.6f} min, nominal margin {margin_sp_h200_s:.3f} s; "
                "centered 120 um VPIOs fail clean at 3.070 / 3.385 min even where Hybrid H180 nominal oxygen passes). "
                "Borderline H200 remains an unqualified scenario pending empirical chamber-corrected loaded-envelope perifusion (R3P-04)."
            ),
        },
        "cases": rows,
        "spatially_weighted_1d_projection_cases": spatial_rows,
    }


def compute_sequential_colony_screening_audit() -> Dict[str, Any]:
    """Forward binomial colony screening math from explicitly declared planning rates (no back-solving)."""
    p_single_17kb = 0.00219   # 0.219% assumed single-step 17.2 kb efficiency
    p_step1_10p5kb = 0.01294  # 1.294% assumed Step-1 10.45 kb attP-GT efficiency
    p_step2_6p9kb = 0.02117   # 2.117% assumed Step-2 6.91 kb attP-CA efficiency

    def _n_for_conf(p: float, target_conf: float) -> int:
        return int(math.ceil(math.log(1.0 - target_conf) / math.log(1.0 - p)))

    def _joint_prob(n1: int, n2: int) -> float:
        return (1.0 - ((1.0 - p_step1_10p5kb) ** n1)) * (1.0 - ((1.0 - p_step2_6p9kb) ** n2))

    per_stage_for_95_joint = math.sqrt(0.95)
    n_single_95 = _n_for_conf(p_single_17kb, 0.95)
    n1_95_stage = _n_for_conf(p_step1_10p5kb, 0.95)
    n2_95_stage = _n_for_conf(p_step2_6p9kb, 0.95)
    n1_95_joint = _n_for_conf(p_step1_10p5kb, per_stage_for_95_joint)
    n2_95_joint = _n_for_conf(p_step2_6p9kb, per_stage_for_95_joint)
    exact_joint_per_stage_95 = _joint_prob(n1_95_stage, n2_95_stage)
    exact_joint_true_95 = _joint_prob(n1_95_joint, n2_95_joint)

    # Exact exhaustive integer minimization for >= 95% joint probability of >= 1 colony at each step
    best_total = n1_95_joint + n2_95_joint
    min_allocs: List[List[int]] = []
    for n1 in range(1, best_total + 1):
        p1 = 1.0 - ((1.0 - p_step1_10p5kb) ** n1)
        if p1 < 0.95:
            continue
        rem = 1.0 - (0.95 / p1)
        n2 = int(math.ceil(math.log(rem) / math.log(1.0 - p_step2_6p9kb)))
        if _joint_prob(n1, n2) < 0.95:
            n2 += 1
        tot = n1 + n2
        if tot < best_total:
            best_total = tot
            min_allocs = [[n1, n2]]
        elif tot == best_total:
            min_allocs.append([n1, n2])

    return {
        "assumed_forward_per_colony_rates": {
            "single_step_17kb": p_single_17kb,
            "step1_10p5kb_attp_gt": p_step1_10p5kb,
            "step2_6p9kb_attp_ca": p_step2_6p9kb,
        },
        "transferred_planning_rates_note": (
            "Forward planning calculation of P(>=1 acceptable colony at each sequential step) from explicitly declared "
            "independent-draw planning rates (not back-solved from colony counts); Pandey 2024/2025 (PMID:38858586) "
            "provides no exact-product yield for Design Y, and direct hiPSC integration rates remain unmeasured pending Stage-0 Y-05."
        ),
        "single_step_17kb_n95": n_single_95,
        "two_stage_per_stage_95pct_colonies": [n1_95_stage, n2_95_stage],
        "two_stage_per_stage_95pct_total": n1_95_stage + n2_95_stage,
        "two_stage_per_stage_95pct_continuous_lower_bound_pct": round(100.0 * (0.95 * 0.95), 2),
        "two_stage_per_stage_95pct_joint_success_pct": round(100.0 * exact_joint_per_stage_95, 4),
        "per_stage_confidence_for_true_95pct_joint_pct": round(100.0 * per_stage_for_95_joint, 4),
        "two_stage_equal_margin_95pct_joint_colonies": [n1_95_joint, n2_95_joint],
        "two_stage_equal_margin_95pct_joint_total": n1_95_joint + n2_95_joint,
        "two_stage_equal_margin_95pct_joint_exact_success_pct": round(100.0 * exact_joint_true_95, 4),
        "two_stage_exact_minimum_95pct_joint_total": best_total,
        "two_stage_exact_minimum_95pct_joint_allocation_count": len(min_allocs),
        "two_stage_exact_minimum_95pct_joint_allocations": min_allocs,
        "two_stage_equal_margin_455_at_1p0pct_floor_joint_pct": round(
            100.0 * (1.0 - (0.99 ** n1_95_joint)) * (1.0 - (0.99 ** n2_95_joint)), 4
        ),
        "allocation_uniqueness_note": (
            f"{n1_95_joint} + {n2_95_joint} = {n1_95_joint + n2_95_joint} is the equal-margin (sqrt(0.95) per-stage) "
            f"allocation, NOT the minimum total: exact integer minimization gives {best_total} colonies across "
            f"{len(min_allocs)} qualifying splits ({min_allocs[0][0]}+{min_allocs[0][1]} through "
            f"{min_allocs[-1][0]}+{min_allocs[-1][1]})."
        ),
    }


def _resolve_and_verify_design_y_v2_registry(repo_root: Path) -> Dict[str, Any]:
    """Load product registry, staged FASTAs, DEC-113 96-cell decision-tree contract, and CommonMark Playbook s4.1.19."""
    p_defs = pc.load_product_contracts(repo_root, require_staged_fastas=True)
    by_id = {p["product_variant_id"]: p for p in p_defs["products"]}

    # 1. Verify DEC-113 block in stage0_decision_tree_v1.json (strict no-duplicate-keys, no-NaN, canonical SHA-256)
    dt_path = repo_root / DECISION_TREE_REL_PATH
    if not dt_path.is_file():
        raise ValueError(f"Missing Stage-0 decision tree file: {dt_path}")
    dt_doc = _load_json_strict(dt_path)
    dec113_blk = dt_doc.get("design_y_fullstack_v2_and_route3p_gates_iip_0034")
    if not isinstance(dec113_blk, dict):
        raise ValueError("Missing DEC-113 block design_y_fullstack_v2_and_route3p_gates_iip_0034 in stage0_decision_tree_v1.json")
    dec113_sha = hashlib.sha256(
        json.dumps(dec113_blk, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if dec113_sha != EXPECTED_DEC113_DECISION_TREE_SHA256:
        raise ValueError(
            f"DEC-113 decision-tree block SHA-256 mismatch: found {dec113_sha}, expected {EXPECTED_DEC113_DECISION_TREE_SHA256}"
        )

    # 2. Parse Playbook with real CommonMark 0.31.2 block parser & require EXACTLY 1 rendered 4.1.19 heading!
    pb_path = repo_root / PLAYBOOK_REL_PATH
    if not pb_path.is_file():
        raise ValueError(f"Missing Stage-0 execution playbook: {pb_path}")
    pb_raw = pb_path.read_text(encoding="utf-8")
    pb_norm = pb_raw.replace("\r\n", "\n").replace("\r", "\n")
    all_headings = _parse_commonmark_headings(pb_norm)
    matched_4119 = [
        h for h in all_headings
        if re.search(r"(?<![0-9A-Za-z.])4\.1\.19(?![0-9A-Za-z.])", re.sub(r"[`*_]", "", h["rendered_title"]))
    ]
    pb_expected_prefix = "#### 4.1.19 Preregistered `Design Y-FullStack v2`"
    if len(matched_4119) != 1:
        raise ValueError(
            f"CommonMark Playbook audit requires exactly 1 rendered 4.1.19 heading, found {len(matched_4119)}: "
            f"{[(h['kind'], h['level'], h['container_depth'], h['start_line'] + 1) for h in matched_4119]}"
        )
    h4119 = matched_4119[0]
    if (
        h4119["kind"] != "atx"
        or h4119["level"] != 4
        or h4119["container_depth"] != 0
        or not h4119["raw_first_line"].startswith(pb_expected_prefix)
    ):
        raise ValueError(f"Rendered 4.1.19 heading has unexpected CommonMark structure: {h4119}")

    # Determine section end strictly from parser offsets (next top-level heading of level <= 4 OR provenance boundary)
    start_off = h4119["start_offset"]
    next_heading_off = len(pb_norm)
    for h in all_headings:
        if h["start_offset"] > start_off and h["container_depth"] == 0 and h["level"] <= 4:
            next_heading_off = h["start_offset"]
            break
    prov_marker = "\n\nCurrent cHS4 provenance:"
    prov_idx = pb_norm.find(prov_marker, start_off, next_heading_off)
    if prov_idx == -1:
        raise ValueError("Missing expected post-4.1.19 provenance boundary after rendered 4.1.19 heading")
    pb_sec_full = pb_norm[start_off:prov_idx].strip()
    pb_sha = hashlib.sha256(pb_sec_full.encode("utf-8")).hexdigest()
    if pb_sha != EXPECTED_DEC113_PLAYBOOK_SEC4119_SHA256:
        raise ValueError(
            f"DEC-113 Playbook section 4.1.19 SHA-256 mismatch: found {pb_sha}, expected {EXPECTED_DEC113_PLAYBOOK_SEC4119_SHA256}"
        )

    # 3. Verify product roles, non-executable state, and 1:1 route binding
    expected_fullstack_roles = frozenset({"FULL_CHR15_11TU_REFERENCE", "GATE0V_CD55_CFH_BACKUP_COMPARATOR"})
    resolved_fullstack_ids: List[str] = []
    resolved_excluded_ids: List[str] = []
    resolved_fullstack_route_codes: set[str] = set()
    resolved_all_route_codes: set[str] = set()

    for pvid in sorted(pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS):
        if pvid not in by_id:
            raise ValueError(f"Loaded product registry missing staged Design Y-FullStack v2 record {pvid}")
        prod = by_id[pvid]
        if prod.get("execution_supported") is not False:
            raise ValueError(f"Loaded product {pvid} must keep execution_supported=False")
        phys_route = prod.get("physical_route_id")
        mapped_code = PHYSICAL_ROUTE_TO_PORTFOLIO_CODE.get(str(phys_route))
        if mapped_code is None:
            raise ValueError(f"Unmapped physical_route_id {phys_route!r} on {pvid}")
        resolved_all_route_codes.add(mapped_code)
        binding = prod.get("staged_design_y_fullstack_v2_donor_binding") or {}
        role = str(binding.get("role") or "")
        gates_str = str(binding.get("stage0_promotion_gates") or "")
        if role in expected_fullstack_roles:
            if "EXCLUDED_FROM_FULLSTACK_CLEARANCE" in gates_str:
                raise ValueError(f"Full-stack candidate {pvid} must not be marked excluded from full-stack clearance")
            resolved_fullstack_ids.append(pvid)
            resolved_fullstack_route_codes.add(mapped_code)
        else:
            if "EXCLUDED_FROM_FULLSTACK_CLEARANCE" not in gates_str:
                raise ValueError(f"Non-full-stack record {pvid} ({role}) must be marked EXCLUDED_FROM_FULLSTACK_CLEARANCE")
            resolved_excluded_ids.append(pvid)

    portfolio_allowed = rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["DESIGN_Y_FULLSTACK_V2"]
    if resolved_fullstack_route_codes != portfolio_allowed or resolved_all_route_codes != frozenset({"R3P"}):
        raise ValueError(
            f"Loaded product registry routes ({sorted(resolved_fullstack_route_codes)}) disagree with "
            f"ALLOWED_CELL_ROUTE_PAIRS_IIP0034['DESIGN_Y_FULLSTACK_V2'] ({sorted(portfolio_allowed)})"
        )

    # 4. Validate authoritative typed 96-cell product_by_gate_contract & 30-endpoint DEC-102 crosswalk
    #    with strict closed operator/disposition allowlists, duplicate/placeholder rejection, AND
    #    direct cell-to-prose semantic concordance across ALL 3 controlling documents!
    xwalk_rows = dec113_blk.get("dec102_endpoint_crosswalk")
    if not isinstance(xwalk_rows, list) or len(xwalk_rows) != len(EXPECTED_DEC102_ENDPOINT_IDS):
        raise ValueError(
            f"DEC-113 dec102_endpoint_crosswalk must contain exactly {len(EXPECTED_DEC102_ENDPOINT_IDS)} unique entries "
            f"(got {len(xwalk_rows) if isinstance(xwalk_rows, list) else type(xwalk_rows).__name__})"
        )
    seen_xwalk_ids: set[str] = set()
    for xrow in xwalk_rows:
        if not isinstance(xrow, dict) or set(xrow.keys()) != REQUIRED_XWALK_KEYS:
            raise ValueError(f"Malformed dec102_endpoint_crosswalk row: {xrow!r}")
        if any(not isinstance(v, str) or not v.strip() or _PLACEHOLDER_TOKEN_RE.search(v) for v in xrow.values()):
            raise ValueError(f"Empty, non-string, or placeholder field in dec102_endpoint_crosswalk row {xrow.get('dec102_endpoint_id')!r}")
        if xrow["operator"] not in VALID_XWALK_OPERATORS:
            raise ValueError(f"Unrecognized crosswalk operator {xrow['operator']!r} on {xrow['dec102_endpoint_id']!r}")
        if xrow["disposition"] not in VALID_XWALK_DISPOSITIONS:
            raise ValueError(f"Unrecognized crosswalk disposition {xrow['disposition']!r} on {xrow['dec102_endpoint_id']!r}")
        eid = xrow["dec102_endpoint_id"]
        if eid in seen_xwalk_ids:
            raise ValueError(f"Duplicate dec102_endpoint_id rejected in crosswalk: {eid!r}")
        seen_xwalk_ids.add(eid)
    if seen_xwalk_ids != EXPECTED_DEC102_ENDPOINT_IDS:
        raise ValueError(
            f"DEC-113 dec102_endpoint_crosswalk ID set mismatch: missing={sorted(EXPECTED_DEC102_ENDPOINT_IDS - seen_xwalk_ids)}, "
            f"extra={sorted(seen_xwalk_ids - EXPECTED_DEC102_ENDPOINT_IDS)}"
        )

    pbg = dec113_blk.get("product_by_gate_contract")
    if not isinstance(pbg, dict) or pbg.get("cell_count") != 96 or not isinstance(pbg.get("cells"), list) or len(pbg["cells"]) != 96:
        raise ValueError("DEC-113 product_by_gate_contract must define exactly 96 typed product-by-gate cells")
    cell_map: Dict[tuple[str, str], Dict[str, Any]] = {}
    for cell in pbg["cells"]:
        if not isinstance(cell, dict) or set(cell.keys()) != REQUIRED_CELL_KEYS:
            raise ValueError(f"Malformed typed cell in product_by_gate_contract: {cell!r}")
        if any(not isinstance(v, str) or not v.strip() or _PLACEHOLDER_TOKEN_RE.search(v) for v in cell.values()):
            raise ValueError(f"Empty, non-string, or placeholder field in typed cell {cell.get('product_variant_id')}:{cell.get('gate_id')}")
        if cell["operator"] not in VALID_CELL_OPERATORS:
            raise ValueError(f"Unrecognized typed cell operator {cell['operator']!r} on {cell['product_variant_id']}:{cell['gate_id']}")
        if cell["disposition"] not in VALID_CELL_DISPOSITIONS:
            raise ValueError(f"Unrecognized typed cell disposition {cell['disposition']!r} on {cell['product_variant_id']}:{cell['gate_id']}")
        key = (cell["product_variant_id"], cell["gate_id"])
        if key in cell_map:
            raise ValueError(f"Duplicate product-by-gate cell: {key}")
        cell_map[key] = cell
    expected_pairs = {(p, g) for p in pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS for g in EXPECTED_96_GATES}
    if set(cell_map.keys()) != expected_pairs:
        raise ValueError("Typed 96-cell product_by_gate_contract does not match exact 6 products x 16 gates Cartesian set")

    # Verify exact per-(role, gate) operator semantics on all 96 typed cells!
    for pvid in sorted(pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS):
        role = str(by_id[pvid]["staged_design_y_fullstack_v2_donor_binding"]["role"])
        for gid in EXPECTED_96_GATES:
            cell = cell_map[(pvid, gid)]
            if role in expected_fullstack_roles:
                exp_op = EXPECTED_FULLSTACK_GATE_OPERATORS[gid]
                exp_disp = (
                    "REQUIRED_EXACT_PRODUCT"
                    if role == "FULL_CHR15_11TU_REFERENCE"
                    else ("REQUIRED_AFTER_PRIMARY_GATE0V_FAIL_AND_DELIVERY_HOLD" if gid == "Y-05" else "REQUIRED_AFTER_PRIMARY_GATE0V_FAIL")
                )
            elif role in ("SEQUENTIAL_STEP1_ATTB_GT_DONOR", "SEQUENTIAL_STEP2_ATTB_CA_DONOR") and gid == "Y-05":
                exp_op, exp_disp = "CONJUNCTIVE_PASS", "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK"
            elif role == "DEC053_CHR13_CLYBL_SAFETY_SWITCH" and gid == "Y-09":
                exp_op, exp_disp = "CONJUNCTIVE_PASS", "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK"
            elif role == "EQUAL_SPACING_ATTRIBUTION_CONTROL" and gid in ("Y-05", "Y-07"):
                exp_op, exp_disp = "COMPARATOR_REFERENCE_ONLY", "CONTROL_BENCHMARK_ONLY_EXCLUDED_FROM_FULLSTACK"
            elif gid == "STAGE0_ENDPOINT_CROSSWALK":
                exp_op = "EXCLUDED_FROM_FULLSTACK_CLEARANCE"
                exp_disp = (
                    "CONTROL_BENCHMARK_ONLY_EXCLUDED_FROM_FULLSTACK"
                    if role == "EQUAL_SPACING_ATTRIBUTION_CONTROL"
                    else "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK"
                )
            else:
                exp_op, exp_disp = "NOT_APPLICABLE", "NOT_APPLICABLE_EXCLUDED_FROM_FULLSTACK"
            if cell["operator"] != exp_op or cell["disposition"] != exp_disp:
                raise ValueError(
                    f"Typed cell operator/disposition mismatch on {pvid}:{gid}: "
                    f"got ({cell['operator']!r}, {cell['disposition']!r}), expected ({exp_op!r}, {exp_disp!r})"
                )

    dt_y_rule = str((dec113_blk.get("rules") or {}).get("y05_y09_design_y_fullstack_v2_gates") or "")
    dt_r3p_rule = str((dec113_blk.get("rules") or {}).get("r3p01_r3p06_modular_flat_ribbon_planar_envelope_gates") or "")
    dt_xwalk = str((dec113_blk.get("rules") or {}).get("stage0_endpoint_crosswalk") or "")
    pb_clean = pb_sec_full.replace("`", "")

    # Contradiction guards (reject sign reversals, unpinned floor approximations, or broken triggers in ANY doc/cell)
    forbidden_patterns = (
        (r"(?i)iCasp9\s+kill\s*\(?<=\s*1e-8", "Y-04 sign-reversed 'kill <= 1e-8'"),
        (r"(?i)passes\s+Gate\s+0V", "reversed backup trigger 'passes Gate 0V'"),
        (r"n_eff\s*=\s*floor\(", "forbidden floor(n_raw*eta) approximation"),
        (r"(?i)3x\s+ectopic[- ]GPI", "misleading '3x' GPI wording instead of '3 ectopic GPI'"),
    )
    all_checked_texts: List[tuple[str, str]] = [
        ("decision_tree.y05_y09", dt_y_rule),
        ("decision_tree.r3p01_r3p06", dt_r3p_rule),
        ("decision_tree.crosswalk", dt_xwalk),
        ("playbook.sec4119", pb_clean),
    ]
    for (cpvid, cgid), cobj in sorted(cell_map.items()):
        all_checked_texts.append((f"typed_cell.{cpvid}:{cgid}", json.dumps(cobj, ensure_ascii=False)))
    for pvid in sorted(pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS):
        all_checked_texts.append((f"registry.{pvid}", json.dumps(by_id[pvid], ensure_ascii=False)))
    for dlabel, dtxt in all_checked_texts:
        for pat, desc in forbidden_patterns:
            if re.search(pat, dtxt):
                raise ValueError(f"Contradiction guard tripped ({desc}) in {dlabel}")

    # Direct typed-cell AND cross-document predicate concordance (gate_id, cell_threshold_regex, prose_regex)!
    # Checks BOTH the actual typed cell thresholds in product_by_gate_contract AND all 3 controlling prose views!
    predicate_specs: List[tuple[str, str, str, str]] = [
        ("Y-01_cbe_offtarget_0p10", "Y-01", r"CBE_offtarget_U95<=0\.10%", r"CBE\s+off-target\s+U95\s*<=\s*0\.10%"),
        ("Y-01_5_controls_explicit_se", "Y-01", r"5_supporting_controls==True;EXPLICIT_SE_SUPPLIED", r"EXPLICIT_SE_SUPPLIED"),
        ("Y-03_paracrine_gsis_adenosine_dpcpx", "Y-03", r"coculture_GSIS_SI_L95>=2\.5;GSIS_retention_L95>=90%;supernatant_adenosine_U95<=25\.0nM;DPCPX_GSIS_delta_U95<=10\.0%", r"2\.5[^;\n]{0,120}90%[^;\n]{0,120}25\.0\s*nM[^;\n]{0,120}10\.0%"),
        ("Y-04_escape_bound", "Y-04", r"clonogenic_escape_survival_U95<=1e-8\(k==0,N_eff>=299573228\)", r"clonogenic\s+escape/survival[^;\n]{0,60}<=\s*1(?:\.0)?e-8"),
        ("Y-04_standalone_gsis_3p0", "Y-04", r"obs_gsis_si_L95>=3\.0", r"obs_gsis_si[^;\n]{0,20}>=\s*3\.0"),
        ("Y-05_hold_interval", "Y-05", r"0\.8%.*1\.0%.*HOLD", r"0\.8%[^;\n]{0,30}1\.0%[^;\n]{0,25}HOLD"),
        ("Y-05_backup_hold", "Y-05", r"HOLD_PENDING_DEDICATED_12TU_BACKUP_DELIVERY_DONOR_SPEC", r"HOLD_PENDING_DEDICATED_12TU_BACKUP_DELIVERY_DONOR_SPEC"),
        ("Y-05_pigk_bound", "Y-05", r"CD59_and_CD47_IgV_GPI_ratio_L95>=0\.80", r"HOLD_PIGK_GPI_COMPETITION"),
        ("Y-05_cd47_floor", "Y-05", r"CD47_L95>=180000_molecules/cell", r"180,000"),
        ("Y-06_adcc_u95", "Y-06", r"ADCC_total_and_paired_IgG_U95<=5\.0%", r"ADCC\s+U95\s*<=\s*5\.0%"),
        ("Y-06_viral_u95", "Y-06", r"Gate0V_viral_k==0_U95<=0\.0005\(n>=5990\)", r"k\s*==\s*0[^;\n]{0,70}0\.0005[^;\n]{0,45}5,990"),
        ("Y-06_c5b9_u95", "Y-06", r"C5b9_U95<=5\.0%", r"C5b-9\s+U95\s*<=\s*5\.0%"),
        ("Y-07_il1rn_window", "Y-07", r"IL1Ra:IL1beta_in_\[10x,100x\]", r"10[-–]100x"),
        ("Y-07_conjunctive_si_auc", "Y-07", r"pulse1_and_pulse2_SI_and_AUC", r"SI\s+AND\s+first-phase\s+AUC"),
        ("Y-07_lcb_0p85", "Y-07", r"L95>=0\.85", r"0\.85"),
        ("Y-08_vtcn1_fsp1", "Y-08", r"VTCN1_shed_U95_t<=10\.0%;FSP1_AIFM2_viable_L95_t>=80\.0%", r"10%[^;\n]{0,90}80%"),
        ("Y-09_leak_and_4log", "Y-09", r"G0_Cpep\+_Ki67-_leak_U95_t<=0\.1%;AZT_dividing_log10_reduction_L95_t>=4\.0.*ZM447439<=0\.005uM", r"0\.1%[^;\n]{0,120}4-log10[^;\n]{0,200}0\.005\s*uM"),
        ("R3P-01_clearance", "R3P-01", r"d_cluster<H_slot;d_cluster\+2\*fibrin<=H_slot", r"d_cluster\s*<\s*H_slot\s+and\s+d_cluster\s*\+\s*2\*fibrin\s*<=\s*H_slot"),
        ("R3P-02_per_ribbon", "R3P-02", r"F_pull<=15\.0N_per_individually_extracted_ribbon", r"15\.0\s*N[^;\n]{0,45}ribbon"),
        ("R3P-03_ceus_30mmhg", "R3P-03", r"per_envelope_pO2>=30\.0mmHg_AND_CEUS_Doppler", r"CEUS/Doppler[^;\n]{0,60}30\.0\s*mmHg"),
        ("R3P-04_oxygen_shutoff", "R3P-04", r"Day0_worst_face>=15\.0mmHg;Day0_core>=5\.0mmHg;post_anastomosis_2p5x_GSIS_core>=5\.0mmHg;empirical_chamber_corrected_fall_t50<3\.0min", r"15\.0[^\n]{0,80}5\.0\s*mmHg[^\n]{0,180}3\.0\s*min"),
        ("R3P-05_single_dose_93350", "R3P-05", r"single_dose_93350_holds_at_3\.20909e-5", r"93,350[^;\n]{0,120}3\.20909e-5"),
        ("R3P-05_provenance_alpha", "R3P-05", r"f_clonogenic_provenance_and_alpha_allocation\(alpha_eta>=0\.05_when_eta<1", r"alpha_eta\s*>=\s*0\.05"),
        ("R3P-06_ada_targets", "R3P-06", r"140_180_vs_noncritical_perioperative_100_180_mg_dL", r"140[-–]180[^;\n]{0,80}100[-–]180\s*mg/dL"),
    ]
    dt_combined = dt_y_rule + "\n" + dt_r3p_rule
    for doc_name, blob in (
        ("decision_tree.rules", dt_combined),
        ("decision_tree.crosswalk", dt_xwalk),
        ("playbook.sec4119", pb_clean),
    ):
        for pred_name, _gid, _cell_pat, prose_pat in predicate_specs:
            if pred_name == "R3P-06_ada_targets" and doc_name == "decision_tree.crosswalk":
                continue
            if pred_name == "R3P-05_single_dose_93350" and doc_name == "decision_tree.crosswalk":
                if "93,350" not in blob:
                    raise ValueError(f"Missing single-dose 93,350 disclosure in {doc_name}")
                continue
            if not re.search(prose_pat, blob):
                raise ValueError(f"Typed predicate parity failure ({pred_name}) in {doc_name}")

    # Verify EVERY individual product record in registry AND its ACTUAL 16 typed cells match!
    product_summaries: List[Dict[str, Any]] = []
    for pvid in sorted(pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS):
        prod = by_id[pvid]
        bind = prod["staged_design_y_fullstack_v2_donor_binding"]
        role = str(bind["role"])
        gates_txt = " ".join(prod.get("pending_assay_gates") or [])
        if role == "FULL_CHR15_11TU_REFERENCE":
            for pred_name, gid, cell_pat, prose_pat in predicate_specs:
                if pred_name in ("Y-05_backup_hold", "Y-05_pigk_bound", "Y-05_cd47_floor"):
                    continue
                if not re.search(cell_pat, cell_map[(pvid, gid)]["threshold"]):
                    raise ValueError(f"Typed cell threshold semantic mismatch ({pred_name}) on {pvid}:{gid}")
                if not re.search(prose_pat, gates_txt):
                    raise ValueError(f"Registry predicate mismatch ({pred_name}) on {pvid}")
        elif role == "GATE0V_CD55_CFH_BACKUP_COMPARATOR":
            for pred_name, gid, cell_pat, prose_pat in predicate_specs:
                if not re.search(cell_pat, cell_map[(pvid, gid)]["threshold"]):
                    raise ValueError(f"Typed cell threshold semantic mismatch ({pred_name}) on {pvid}:{gid}")
                if not re.search(prose_pat, gates_txt):
                    raise ValueError(f"Registry predicate mismatch ({pred_name}) on {pvid}")
        elif role == "EQUAL_SPACING_ATTRIBUTION_CONTROL":
            if "0.80" not in cell_map[(pvid, "Y-05")]["threshold"] or "0.85" not in cell_map[(pvid, "Y-07")]["threshold"]:
                raise ValueError(f"Control typed cells Y-05/Y-07 missing ratio thresholds on {pvid}")
            for req in ("Y-05", "CD59", "CD47-IgV-GPI", "0.80", "180,000", "Y-07", "SI AND first-phase AUC", "0.85"):
                if req not in gates_txt:
                    raise ValueError(f"Control registry {pvid} missing Y-05/Y-07 comparator clause {req!r}")
        elif role in ("SEQUENTIAL_STEP1_ATTB_GT_DONOR", "SEQUENTIAL_STEP2_ATTB_CA_DONOR"):
            if "0.8%-1.0%_HOLD" not in cell_map[(pvid, "Y-05")]["threshold"] or "384" not in cell_map[(pvid, "Y-05")]["threshold"]:
                raise ValueError(f"Component typed cell Y-05 missing HOLD interval / colony cap on {pvid}")
            if "0.8% <= rate < 1.0% -> HOLD" not in gates_txt or "384" not in gates_txt:
                raise ValueError(f"Component registry {pvid} missing Y-05 HOLD interval / colony cap")
        elif role == "DEC053_CHR13_CLYBL_SAFETY_SWITCH":
            if "0.1%" not in cell_map[(pvid, "Y-09")]["threshold"] or "4.0_log10" not in cell_map[(pvid, "Y-09")]["threshold"]:
                raise ValueError(f"Component typed cell Y-09 missing thresholds on {pvid}")
            if "0.1%" not in gates_txt or "4-log10" not in gates_txt or "0.005 uM" not in gates_txt:
                raise ValueError(f"Component registry {pvid} missing Y-09 thresholds")
            if "Not encoded by this standalone chr13 donor" not in str((prod.get("hla_retention") or {}).get("classical_MHC_I") or ""):
                raise ValueError(f"Standalone DEC053 record {pvid} must distinguish encoded payload from intended host genotype")

        dispositions_by_gate = {g: cell_map[(pvid, g)]["disposition"] for g in EXPECTED_96_GATES}
        product_summaries.append({
            "product_variant_id": pvid,
            "role": role,
            "fullstack_clearance_eligible": role in expected_fullstack_roles,
            "gate_dispositions": dispositions_by_gate,
            "typed_cell_operators_and_thresholds_verified": True,
            "registry_predicates_verified_against_typed_contract": True,
        })

    return {
        "registry_loaded_via_product_contract": True,
        "dec113_decision_tree_sha256_verified": dec113_sha,
        "dec113_playbook_sec4119_sha256_verified": pb_sha,
        "playbook_commonmark_parser": {
            "rendered_4119_headings_found": len(matched_4119),
            "heading_kind": h4119["kind"],
            "heading_level": h4119["level"],
            "container_depth": h4119["container_depth"],
            "total_document_headings_parsed": sum(1 for h in all_headings if h["start_offset"] <= h4119["start_offset"]),
        },
        "typed_product_by_gate_contract_audit": {
            "products_checked": len(product_summaries),
            "gates_per_product": len(EXPECTED_96_GATES),
            "cells_verified": len(cell_map),
            "dec102_endpoints_mapped": len(xwalk_rows),
            "contradiction_guards_passed": len(forbidden_patterns),
            "typed_predicates_verified": len(predicate_specs),
            "product_gate_matrix": product_summaries,
        },
        "staged_fasta_sha256_verified": pc.DESIGN_Y_V2_FASTA_SHA256,
        "staged_design_y_v2_record_count": len(pc.STAGED_DESIGN_Y_V2_PRODUCT_IDS),
        "resolved_fullstack_candidate_ids": resolved_fullstack_ids,
        "resolved_excluded_component_and_control_ids": resolved_excluded_ids,
        "resolved_design_y_v2_routes_from_registry": sorted(resolved_fullstack_route_codes),
        "design_y_v2_bound_strictly_to_route3p_matching_registry": True,
        "scope_note": (
            "Verifies exact FASTA file/record lengths and SHA-256s, whole-record canonical product-entry SHA-256s, "
            "canonical JSON SHA-256 of the DEC-113 stage0_decision_tree_v1.json block (including its 96-cell typed "
            "product_by_gate_contract and 30-endpoint DEC-102 crosswalk matching wetlab_contract.py:5526-5634), "
            "LF-normalized trimmed CommonMark-parsed Playbook section 4.1.19 SHA-256 (with HTML character-reference "
            "decoding), closed cell/crosswalk operator allowlists, placeholder/duplicate rejection, direct typed-cell "
            "threshold concordance against all 3 controlling documents, execution_supported=False, role-scoped eligibility, "
            "and 1:1 physical_route_id alignment against ALLOWED_CELL_ROUTE_PAIRS_IIP0034; does not claim verbatim prose "
            "equality across differently structured views and does not independently recompute full sequence-grammar motif scans."
        ),
    }


def compute_full_verification(repo_root: Path | None = None) -> Dict[str, Any]:
    active_root = ROOT if repo_root is None else Path(repo_root)
    registry_audit = _resolve_and_verify_design_y_v2_registry(active_root)
    r3p05_audit = compute_r3p05_dec103_thinning_boundary_audit()
    mechanics = compute_planar_mechanics_table()
    fvm_cases = [
        solve_planar_route3p_fvm(
            label="Naive_147cm2_H180_Day0_Clean_FAIL",
            area_cm2=147.0, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H180_Day0_Clean_65um_PASS",
            area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H200_Day0_Clean_65um_PASS",
            area_cm2=300.0, h_slot_um=200.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_Clean_65um_PASS",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H200_Day0_Clean_65um_PASS",
            area_cm2=333.3, h_slot_um=200.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_10umFibrin_0p5D_65um_FACE_FAIL",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=10.0, sleeve_d_factor=0.5,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_15umFibrin_0p5D_65um",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=15.0, sleeve_d_factor=0.5,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_30umDenseThrombus_0p5D_CORE_FAIL",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=30.0, sleeve_d_factor=0.5,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H180_Day0_30umDenseThrombus_0p5D_CORE_FAIL",
            area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=35.0, sleeve_um=30.0, sleeve_d_factor=0.5,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H180_Day0_Clean_65um_Bed30mmHg_FACE_FAIL",
            area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=30.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H200_Day0_Clean_65um_Bed30mmHg_FACE_FAIL",
            area_cm2=300.0, h_slot_um=200.0, cluster_diam_um=65.0,
            bed_po2_mmhg=30.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_Clean_65um_Bed30mmHg_PASS",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=30.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H200_Day0_Clean_65um_Bed30mmHg_FACE_FAIL",
            area_cm2=333.3, h_slot_um=200.0, cluster_diam_um=65.0,
            bed_po2_mmhg=30.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="FourMod260_H180_Day0_Clean_65um_Bed30mmHg_FACE_FAIL",
            area_cm2=260.0, h_slot_um=180.0, cluster_diam_um=65.0,
            bed_po2_mmhg=30.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H180_Day0_Avascular_VPIO120_Clean_FACE_FAIL",
            area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_Avascular_VPIO120_Clean_PASS",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H180_Day0_Avascular_VPIO120_14p9umFibrin_0p5D_FAIL",
            area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=14.9, sleeve_d_factor=0.5,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_Day0_Avascular_VPIO120_14p9umFibrin_0p5D_FAIL",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=14.9, sleeve_d_factor=0.5,
            k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        ),
        solve_planar_route3p_fvm(
            label="SubARS300_H180_VPIO120_PostAnast_2p5xGSIS_45umFBR_PASS",
            area_cm2=300.0, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=45.0, sleeve_d_factor=1.0 / 3.0,
            k_perf_gel_s_inv=5.0, k_perf_islet_s_inv=3.0, ocr_mult=1.097 * 2.5,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_VPIO120_PostAnast_2p5xGSIS_45umFBR_PASS",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=45.0, sleeve_d_factor=1.0 / 3.0,
            k_perf_gel_s_inv=5.0, k_perf_islet_s_inv=3.0, ocr_mult=1.097 * 2.5,
        ),
        solve_planar_route3p_fvm(
            label="Hybrid333_H180_VPIO120_PostAnast_GelOnlyNoIsletPerf_CORE_FAIL",
            area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=120.0,
            bed_po2_mmhg=35.0, sleeve_um=45.0, sleeve_d_factor=1.0 / 3.0,
            k_perf_gel_s_inv=5.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097 * 2.5,
        ),
    ]
    shutoff = compute_planar_insulin_shutoff_table()
    colonies = compute_sequential_colony_screening_audit()
    solubility_sensitivity_30mmhg = solve_planar_route3p_fvm(
        label="Hybrid333_H180_Day0_Clean_65um_Bed30mmHg_Alpha1p27e9_FACE_FAIL",
        area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=65.0,
        bed_po2_mmhg=30.0, sleeve_um=0.0, sleeve_d_factor=1.0,
        k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        alpha_mol_cm3_mmhg=ALPHA_PHYS_BUFFER,
    )
    solubility_sensitivity_120um_bed35 = solve_planar_route3p_fvm(
        label="Hybrid333_H180_Day0_Clean_VPIO120_Bed35mmHg_Alpha1p27e9_FAIL",
        area_cm2=333.3, h_slot_um=180.0, cluster_diam_um=120.0,
        bed_po2_mmhg=35.0, sleeve_um=0.0, sleeve_d_factor=1.0,
        k_perf_gel_s_inv=0.0, k_perf_islet_s_inv=0.0, ocr_mult=1.097,
        alpha_mol_cm3_mmhg=ALPHA_PHYS_BUFFER,
    )

    core_fn = getattr(rank_mod, "_compute_multiparadigm_2pow4n_portfolio_iip0034_core", None)
    if core_fn is None:
        core_fn = rank_mod.compute_multiparadigm_2pow4n_portfolio_iip0034
    mp = core_fn()
    all_rows = (
        list(mp.get("conservative_zero_h2_h3_lift_4x3_architectural", []))
        + list(mp.get("stage0_gated_conditional_4x3_off_the_shelf", []))
        + list(mp.get("stage0_gated_5arm_unified_bridge", []))
    )
    max_4096_err = max(float(r["incumbent_4096_state_max_abs_err"]) for r in all_rows)
    if not math.isfinite(max_4096_err) or max_4096_err > 1e-12:
        raise ValueError(f"4096-state parity violation: {max_4096_err}")

    pa = (0.88, 0.78, 0.57, 0.69)
    pb = (0.88, 0.79, 0.64, 0.70)
    px = (0.88, 0.82, 0.73, 0.69)
    c_base = (0.88, 0.87, 0.89, 0.85)
    live = rank_mod.compute_common_shock_2pow12_portfolio_tree(pa, pb, px, c_base)
    inadm_diag = rank_mod._run_dp_and_verify_parity_iip0034(
        "INADMISSIBLE_DIAGONAL_SELF_CHECK",
        [pa, pb, px],
        rank_mod._live_canonical_incumbent_dists_iip0034(),
        list(c_base),
        float(live["locked_2route_shared_ax_pouch_exact"]),
        float(live["modular_2route_shared_ax_pouch_exact"]),
        math.prod(c_base),
        stage0_status="SELF_CHECK_INADMISSIBLE_DIAGONAL",
        cell_classes=("UNHARDENED_A_B_X", "UNHARDENED_A_B_X", "UNHARDENED_A_B_X"),
        route_ids=("R3P", "R3P", "R3P"),
    )
    if not (inadm_diag["locked_pct"] == 0.0 and inadm_diag["modular_k100_pct"] == 0.0 and inadm_diag["modular_k90_pct"] == 0.0):
        raise ValueError(f"Inadmissible diagonal leaked non-zero survival: {inadm_diag}")

    by_sc = {r["scenario"]: r for r in all_rows}
    pareto_path = active_root / PARETO_REL_PATH
    if not pareto_path.is_file():
        raise ValueError(f"Missing Pareto artifact: {pareto_path}")
    pareto_sha = hashlib.sha256(pareto_path.read_bytes()).hexdigest()
    if pareto_sha != EXPECTED_PARETO_SHA256:
        raise ValueError(f"Pareto SHA-256 drift: {pareto_sha} != {EXPECTED_PARETO_SHA256}")

    return {
        "schema_version": "openislet_iip0034_planar_and_portfolio_verification_v0.8.0",
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "license": "CC-BY-4.0",
        "epistemic_status": "DETERMINISTIC_CONTINUUM_AND_PORTFOLIO_REPLAY_STAGE0_PENDING",
        "live_headline_odds_delta_pp": 0.0,
        "pareto_sha256_verified": pareto_sha,
        "design_y_v2_product_registry_audit": registry_audit,
        "r3p05_dec103_thinning_boundary_audit": r3p05_audit,
        "planar_mechanics_and_clearance": mechanics,
        "planar_3d_fvm_oxygen_envelopes": fvm_cases,
        "fvm_model_assumptions_and_convergence_note": (
            "Periodic unit-cell 3D FVM solved on 28x28x32 meshes (grid refinement to 56x56x64-80x80x96 shifts "
            "marginal 120 um face/core values by tenths of mmHg without altering any combined face/core gate class). "
            "Post-anastomosis 2.5x GSIS OCR multiplier scales both islet and stromal gel sinks conservatively; "
            "direct intra-islet perfusion is a modeled scenario pending loaded-envelope measurement."
        ),
        "planar_30mmhg_buffer_solubility_sensitivity": {
            "alpha_nominal_mol_cm3_mmhg": ALPHA_LOW_SALT,
            "alpha_phys_buffer_mol_cm3_mmhg": ALPHA_PHYS_BUFFER,
            "hybrid333_h180_bed30_clean_phys_buffer": solubility_sensitivity_30mmhg,
            "hybrid333_h180_bed35_120um_clean_phys_buffer": solubility_sensitivity_120um_bed35,
            "caveat": (
                "Switching Bunsen solubility from nominal 1.34e-9 to physiological-buffer 1.27e-9 mol/(cm^3*mmHg) "
                "drops Hybrid-333 H180 clean 65 um at a 30 mmHg bed from 15.2152 to 14.4198 mmHg (< 15.0 mmHg face floor) "
                "and drops Day-0 avascular 120 um Hybrid-333 H180 at a 35 mmHg bed from 15.3775/5.9609 to "
                "14.3801/4.5487 mmHg (failing both face and core floors); additionally 10 um of 0.5x D fibrin at 35 mmHg "
                "drops 65 um Hybrid-333 H180 face to 14.72 mmHg (< 15.0 mmHg) despite being inside the 30 um thickness cap."
            ),
        },
        "planar_1d_insulin_shutoff_fall_t50": shutoff,
        "sequential_bxb1_colony_screening_audit": colonies,
        "portfolio_compatibility_and_4096_state_audit": {
            "allowed_cell_route_pairs": {
                k: sorted(v) for k, v in sorted(rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034.items())
            },
            "unhardened_abx_barred_from_route3p": (
                "R3P" not in rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["UNHARDENED_A_B_X"]
                and "R3P" not in rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["UNEDITED_AUTO_NEGATIVE_CONTROL"]
                and "R3P" in rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["DESIGN_Y_FULLSTACK_V2"]
                and "R3P" in rank_mod.ALLOWED_CELL_ROUTE_PAIRS_IIP0034["D1B_AUTO_C_V2"]
            ),
            "design_y_v2_bound_strictly_to_route3p_matching_registry": registry_audit["design_y_v2_bound_strictly_to_route3p_matching_registry"],
            "locked_path_compatibility_enforced": True,
            "inadmissible_diagonal_locked_pct": inadm_diag["locked_pct"],
            "inadmissible_diagonal_modular_k100_pct": inadm_diag["modular_k100_pct"],
            "inadmissible_diagonal_modular_k90_pct": inadm_diag["modular_k90_pct"],
            "ordinary_scenario_count_checked": len(all_rows),
            "total_evaluated_rows_including_sensitivities": (
                len(all_rows)
                + len(mp.get("phi_h4_sensitivity_counterexamples", []))
                + len(mp.get("cross_hurdle_negative_coupling_sensitivity", []))
            ),
            "max_4096_state_incumbent_abs_err": max_4096_err,
            "w_b1_modular_compatible_pct": by_sc["W_B1_4Arm_Y@PrePerf3P_ZeroLift_c4u0.92_phiYX0.75_phiH4=0.25"]["modular_k100_pct"],
            "w_b1_prohibited_modular_mass_excluded_pp": by_sc["W_B1_4Arm_Y@PrePerf3P_ZeroLift_c4u0.92_phiYX0.75_phiH4=0.25"]["prohibited_modular_mass_excluded_pp"],
            "w_b5_modular_compatible_pct": by_sc["W_B5_4Arm_Y@PrePerf3P_Nested(90/92)_c4u0.92_phiYX0.75_phiH4=0.25"]["modular_k100_pct"],
            "w_b5_prohibited_modular_mass_excluded_pp": by_sc["W_B5_4Arm_Y@PrePerf3P_Nested(90/92)_c4u0.92_phiYX0.75_phiH4=0.25"]["prohibited_modular_mass_excluded_pp"],
            "w_b9_locked_phiH4_0p50_pct": by_sc["W_B9_4Arm_Y@PrePerf3P_Nested(90/92)_pY4=0.68_c4u0.90_phi1.0_50"]["locked_pct"],
            "w_b10_prohibited_modular_mass_excluded_pp": by_sc["W_B10_4Arm_Y@Hybrid3P_Nested(90/92)_c4u0.92_phiYX0.75_phiH4=0.25"]["prohibited_modular_mass_excluded_pp"],
            "w_b9_phiH4_0p60_sensitivity_counterexample": mp["phi_h4_sensitivity_counterexamples"][0],
            "w_b3_equalized_y_marginals_sensitivity_counterexample": mp["phi_h4_sensitivity_counterexamples"][1],
            "cross_hurdle_negative_coupling_sensitivity": mp["cross_hurdle_negative_coupling_sensitivity"],
            "five_arm_bridge_registration_status": mp["five_arm_bridge_registration_status"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Deterministic planar Route 3P & portfolio verifier (IIP-0034 / DEC-113).")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help=f"Write canonical JSON artifact to {DEFAULT_OUT_PATH}")
    mode.add_argument("--check", action="store_true", help=f"Verify byte-identical match against {DEFAULT_OUT_PATH}")
    args = parser.parse_args()

    payload = compute_full_verification(ROOT)
    serialized = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    expected_bytes = serialized.encode("utf-8")
    expected_sha = hashlib.sha256(expected_bytes).hexdigest()
    if args.write:
        DEFAULT_OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        DEFAULT_OUT_PATH.write_bytes(expected_bytes)
        print(f"Wrote {DEFAULT_OUT_PATH} ({len(expected_bytes)} bytes, SHA-256={expected_sha})")
    elif args.check:
        if not DEFAULT_OUT_PATH.is_file():
            raise SystemExit(f"Missing verification artifact: {DEFAULT_OUT_PATH}")
        existing_bytes = DEFAULT_OUT_PATH.read_bytes()
        existing_sha = hashlib.sha256(existing_bytes).hexdigest()
        if existing_bytes != expected_bytes or existing_sha != expected_sha:
            raise SystemExit(
                f"Verification artifact byte/SHA-256 drift detected in {DEFAULT_OUT_PATH}: "
                f"found {existing_sha}, expected {expected_sha}"
            )
        if not PORTFOLIO_PATH.is_file():
            raise SystemExit(f"Missing shipped portfolio artifact: {PORTFOLIO_PATH}")
        shipped_portfolio_bytes = PORTFOLIO_PATH.read_bytes()
        expected_portfolio = rank_mod.compute_multiparadigm_2pow4n_portfolio_iip0034()
        expected_portfolio_bytes = (
            json.dumps(expected_portfolio, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        ).encode("utf-8")
        shipped_p_sha = hashlib.sha256(shipped_portfolio_bytes).hexdigest()
        expected_p_sha = hashlib.sha256(expected_portfolio_bytes).hexdigest()
        if shipped_portfolio_bytes != expected_portfolio_bytes or shipped_p_sha != expected_p_sha:
            raise SystemExit(
                f"Shipped multiparadigm portfolio byte/SHA-256 drift detected in {PORTFOLIO_PATH}: "
                f"found {shipped_p_sha}, expected {expected_p_sha}"
            )
        print(
            f"OK: {DEFAULT_OUT_PATH} (SHA-256={expected_sha}) and "
            f"{PORTFOLIO_PATH} (SHA-256={expected_p_sha}) match deterministic replay byte-for-byte."
        )
    else:
        print(serialized)


if __name__ == "__main__":
    main()
