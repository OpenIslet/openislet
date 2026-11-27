#!/usr/bin/env python3
"""
Deterministic Autologous D1b-Auto-C v2 Mode-B Registry, 45-Cell Stage-0 Contract,
80-Digit Sound Decimal Rare-Event Evaluator & Unified 5-Arm Portfolio Verifier
(IIP-0035 / DEC-114 — Option-2 Zero-LIMS-Touch)
=================================================================================
Self-contained cross-runtime replay verifier for OpenIslet full-repo and public-bundle
checkouts. Depends ONLY on `scripts/product_contract.py`, `scripts/rank_islet1_stacks.py`,
and the Python standard library (ZERO import of `scripts/wetlab_contract.py` or
`scripts/ingest_wetlab_results.py`; ZERO `assert` statements so `python -O` is 100% safe).

Run with:
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0035_d1b_auto_and_portfolio.py
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0035_d1b_auto_and_portfolio.py --write
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0035_d1b_auto_and_portfolio.py --check
"""
from __future__ import annotations

import argparse
import decimal
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, localcontext
from fractions import Fraction


def _clean_dec_ctx(prec: int, rounding: str = ROUND_CEILING):
    """Fresh Decimal Context explicitly specifying bounds/clamp/flags/traps so it is 100% immune
    to both active-context AND mutable decimal.DefaultContext customizations."""
    return localcontext(
        Context(
            prec=prec,
            rounding=rounding,
            Emin=-999999,
            Emax=999999,
            capitals=1,
            clamp=0,
            flags=[],
            traps=[],
        )
    )


_MISSING_ARG = object()
_LOT_ALIAS_KEYS = ("evaluated_lot_id", "lot_id", "batch_id", "mcb_lot_id", "mcb_batch_id")
_LOT_ALIAS_KEY_SET = frozenset(_LOT_ALIAS_KEYS)
import hashlib
import html
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import product_contract as pc  # noqa: E402
import rank_islet1_stacks as rank_mod  # noqa: E402

DEFAULT_OUT_PATH = ROOT / "data" / "extracted" / "iip0035_d1b_auto_and_portfolio_verification_v0.8.0.json"
PORTFOLIO_5ARM_REL_PATH = "data/extracted/islet1_5arm_autologous_portfolio_v0.8.0.json"
PORTFOLIO_5ARM_PATH = ROOT / PORTFOLIO_5ARM_REL_PATH
FROZEN_IIP0034_PORTFOLIO_REL_PATH = "data/extracted/islet1_multiparadigm_portfolio_v0.8.0.json"
FROZEN_IIP0034_VERIF_REL_PATH = "data/extracted/iip0034_planar_and_portfolio_verification_v0.8.0.json"
PARETO_REL_PATH = "data/extracted/islet1_combinatorial_pareto_and_ladder_v0.8.0.json"
DECISION_TREE_REL_PATH = "data/protocols/stage0_decision_tree_v1.json"
PLAYBOOK_REL_PATH = "docs/STAGE-0-EXECUTION-PLAYBOOK.md"

EXPECTED_PARETO_SHA256 = "591ddeb9a706c0bc95fc509765b512aadad8a403287a796e63dd52e7db548575"
EXPECTED_IIP0034_PORTFOLIO_SHA256 = "3ed4d53f1756b8b8615fbbc4e9a15755b9ab4968d206da61160eeb04d833cd45"
EXPECTED_IIP0034_VERIF_SHA256 = "d95970fb9423ff5307e9f481ed7027688fdff7bd558e6c439f91f63c734bd08e"
EXPECTED_DEC114_DECISION_TREE_SHA256 = "1699fa9bd04dedebd4e4846c75a4b72db10d18ec9f06d80e6004fbf37177299f"
EXPECTED_DEC114_PLAYBOOK_SEC4120_SHA256 = "e3fb5ec701678dfdeb62a77c15866327842b60b6b4c503cb86dfc79885a626a8"

PROTECTED_LIMS_SHA256: Dict[str, str] = {
    "scripts/wetlab_contract.py": "331ecd2f2d8602dfceae8b89263d0f360cfadfd0638a46d7d668a376294cfc73",
    "scripts/ingest_wetlab_results.py": "e050593675170515382257ce379449db596c87cc06ffc26d31d128e127630cb1",
    "data/protocols/wetlab_endpoint_contract_v0.8.0.json": "757bc1d9c2722c129cc4fd5962e8995b9f5bc5f7035f8530f5716da9aa2ac26f",
}

D1B_FULLSTACK_CANDIDATE_ROLES = frozenset({
    "PRIMARY_LOCUS1_CIITA_IL1RN_CD64T_GATE0V_CLEAN",
    "SECONDARY_GATE0V_LOCUS1_CIITA_CD46WT_CD64T_COMPARATOR",
})
D1B_FULLSTACK_CANDIDATE_IDS = frozenset({
    "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN",
    "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR",
})
D1B_COMPONENT_DONOR_IDS = frozenset({
    "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1",
    "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9",
    "ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y",
})
EXPECTED_45_GATES = (
    "D1-01", "D1-02", "D1-03", "D1-04", "D1-05", "D1-06", "D1-07", "D1-08", "D1-09",
)
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
VALID_CELL_OPERATORS = frozenset({
    "<=",
    ">=",
    "==",
    "<",
    "COMPOUND_ALL_PASS",
    "NOT_APPLICABLE",
})
VALID_CELL_DISPOSITIONS = frozenset({
    "REQUIRED_EXACT_PRODUCT",
    "REQUIRED_EXACT_PRODUCT_WITH_LOT_MATCHED_GATE0V",
    "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK",
    "NOT_APPLICABLE_EXCLUDED_FROM_FULLSTACK",
})
EXPECTED_L1A_GATE_OPERATORS: Dict[str, str] = {
    "D1-01": "COMPOUND_ALL_PASS",
    "D1-02": "COMPOUND_ALL_PASS",
    "D1-03": "COMPOUND_ALL_PASS",
    "D1-04": "COMPOUND_ALL_PASS",
    "D1-05": "COMPOUND_ALL_PASS",
    "D1-06": "<=",
    "D1-07": "<=",
    "D1-08": "COMPOUND_ALL_PASS",
    "D1-09": "<=",
}
EXPECTED_L1B_GATE_OPERATORS: Dict[str, str] = {
    **EXPECTED_L1A_GATE_OPERATORS,
    "D1-09": "COMPOUND_ALL_PASS",
}
COMPONENT_APPLICABLE_GATES: Dict[str, frozenset[str]] = {
    "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1": frozenset({"D1-08"}),
    "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9": frozenset({"D1-04", "D1-05", "D1-07"}),
    "ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y": frozenset({"D1-04", "D1-05"}),
}
EXPECTED_COMPONENT_GATE_OPERATORS: Dict[str, Dict[str, str]] = {
    "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1": {"D1-08": "COMPOUND_ALL_PASS"},
    "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9": {
        "D1-04": "COMPOUND_ALL_PASS",
        "D1-05": ">=",
        "D1-07": "<=",
    },
    "ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y": {
        "D1-04": "COMPOUND_ALL_PASS",
        "D1-05": "COMPOUND_ALL_PASS",
    },
}
EXCLUDED_HLAC_SP3C4C_CLASSES = frozenset({
    "EXCLUDED_SP3C4C_HLAC07_17_1801_NONFUNCTIONAL_HOMO_OR_COMPOUND_HET_ROUTE_TO_ALLO_Y_OR_D1A",
    "EXCLUDED_SP3C4C_HLAC07_17_HOMOZYGOTE_ROUTE_TO_ALLO_Y_OR_D1A",
})
ALLOWED_HLAC_SP3C4C_CLASSES = frozenset({
    "ELIGIBLE_VL9_COMPETENT_HLAC",
    *EXCLUDED_HLAC_SP3C4C_CLASSES,
})
ALLOWED_CIPSC_SOURCES = frozenset({
    "ADSC_SUN_PROTECTED_CIPSC",
    "BLOOD_ERYTHROID_CIPSC_HEAD_TO_HEAD",
    "SUN_EXPOSED_DERMAL_FIBROBLAST_DISCOURAGED",
})
ALLOWED_LOCUS1_ARMS = frozenset({
    "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
    "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR",
})
VALID_F_ELIGIBLE_PROVENANCE = frozenset({
    "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
    "INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION",
})
PHYSICAL_ROUTE_TO_PORTFOLIO_CODE: Dict[str, str] = {
    "ROUTE_1_MACROENCAPSULATED": "R1",
    "ROUTE_2_OPEN_PORE_PREVASC": "R2",
    "ROUTE_3P_MODULAR_PLANAR_ENVELOPE": "R3P",
}
_PLACEHOLDER_TOKEN_RE = re.compile(r"(?i)\b(?:TBD|TODO|PLACEHOLDER|FIXME|XXX|PENDING_VALUE|DUMMY)\b|\?\?\?")
_FORBIDDEN_ID_CHARS_RE = re.compile(r"[\x00-\x1f\x7f\u200b-\u200f\u2028-\u202f\u2060-\u206f\ufeff]")
_DECIMAL_STR_RE = re.compile(r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$", re.ASCII)

_D1_ENDPOINT_BOUNDS: Dict[str, tuple[float, float | None]] = {
    "obs_d1_02_wt_hlac_allele_retention_pct": (0.0, 100.0),
    "obs_d1_02_hlac_bystander_ct_pct": (0.0, 100.0),
    "obs_d1_02_chr6p21_interlocus_del_pct": (0.0, 100.0),
    "obs_d1_03_chrm_numt_amplicon_depth_x": (0.0, None),
    "obs_d1_03_chrm_novel_nonsyn_het_pct": (0.0, 100.0),
    "obs_d1_04_g0_beta_htmpk_leak_pct": (0.0, 100.0),
    "obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct": (0.0, 100.0),
    "obs_d1_04_dec053_ap1903_icasp9_survivor_frac": (0.0, 1.0),
    "obs_d1_04_dec053_azt_dividing_log10_kill": (0.0, 10.0),
    "obs_d1_05_l3_day28_retention_pct": (0.0, 100.0),
    "obs_d1_05_l4_day28_retention_pct": (0.0, 100.0),
    "obs_d1_05_l3_l4_day28_retention_pct": (0.0, 100.0),
    "obs_d1_05_fsp1_ferroptosis_survival_pct": (0.0, 100.0),
    "obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": (0.0, 100.0),
    "obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct": (0.0, 100.0),
    "obs_d1_08_hlac0304_tcell_lysis_u95_pct": (0.0, 100.0),
    "obs_d1_08_t1d_serum_c5b9_u95_pct": (0.0, 100.0),
    "obs_d1_09_autologous_turnaround_weeks": (1.0, 104.0),
}


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
    """Parse controlling JSON float via exact Decimal first, rejecting nonfinite, overflow, underflow, or binary64 precision loss."""
    try:
        d = Decimal(token)
    except Exception as exc:
        raise ValueError(f"Malformed JSON decimal token rejected: {token!r}") from exc
    if not d.is_finite():
        raise ValueError(f"Non-finite JSON decimal token rejected: {token!r}")
    val = float(token)
    if not math.isfinite(val):
        raise ValueError(f"Non-finite JSON float overflow rejected: {token!r}")
    if Decimal(str(val)) != d:
        raise ValueError(
            f"JSON float loses precision or underflows under binary64 (exact={token!r}, float={val!r}); rejected before canonical hashing"
        )
    return val


def _parse_exact_decimal_json_float(token: str) -> Decimal:
    """Parse JSON float token directly to exact finite Decimal without binary64 rounding/underflow."""
    try:
        d = Decimal(token)
    except Exception as exc:
        raise ValueError(f"Malformed JSON decimal token rejected: {token!r}") from exc
    if not d.is_finite():
        raise ValueError(f"Non-finite JSON decimal token rejected: {token!r}")
    return d


def parse_stage0_observation_json_exact(raw_json: str) -> Dict[str, Any]:
    """Decode raw Stage-0 observation JSON preserving exact decimal numeric lexemes (zero binary64 drift)."""
    obj = json.loads(
        raw_json,
        object_pairs_hook=_reject_duplicate_json_keys,
        parse_constant=_reject_nonfinite_json,
        parse_float=_parse_exact_decimal_json_float,
    )
    if not isinstance(obj, dict):
        raise ValueError("Stage-0 observation JSON must decode to a top-level JSON object")
    return obj


def _load_json_strict(path: Path) -> Any:
    return json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_json_keys,
        parse_constant=_reject_nonfinite_json,
        parse_float=_parse_finite_json_float,
    )


def _safe_exact_decimal(val: Any) -> Decimal | None:
    """Convert non-bool int/float/Decimal/numeric-str to exact finite Decimal without binary64 boundary drift.

    Uses Decimal(str(val)) for finite Python floats (shortest round-trip decimal representation via dtoa)
    so ordinary decoded JSON numbers like 0.1 compare equal to exact decimal string "0.1" at inclusive
    thresholds, while numeric strings parse directly via Decimal(s) with zero binary64 underflow/rounding.
    """
    if isinstance(val, bool):
        return None
    if isinstance(val, Decimal):
        return val if val.is_finite() else None
    if isinstance(val, int):
        return Decimal(val)
    if isinstance(val, float):
        if not math.isfinite(val):
            return None
        try:
            d = Decimal(str(val))
        except Exception:
            return None
        return d if d.is_finite() else None
    if isinstance(val, str):
        if val != val.strip() or not _DECIMAL_STR_RE.fullmatch(val):
            return None
        try:
            d = Decimal(val)
        except Exception:
            return None
        return d if d.is_finite() else None
    return None


def _has_surrogate(val: str) -> bool:
    return any(0xD800 <= ord(ch) <= 0xDFFF for ch in val)


def _is_canonical_identifier(val: Any) -> bool:
    return (
        isinstance(val, str)
        and bool(val)
        and val == val.strip()
        and not _has_surrogate(val)
        and _FORBIDDEN_ID_CHARS_RE.search(val) is None
    )


def _cm_is_blank(s: str) -> bool:
    """CommonMark 0.31.2 §2.2: a blank line contains ONLY ASCII spaces (U+0020) and tabs (U+0009)."""
    return all(ch in (" ", "\t") for ch in s)


def _cm_strip_ascii_ws(s: str) -> str:
    return s.strip(" \t")


def _cm_normalize_link_label(label: str) -> str:
    return re.sub(r"[ \t\n]+", " ", label.strip(" \t\n")).casefold()


def _parse_cm_link_ref_defs_from_text(text: str) -> tuple[set[str], int]:
    """Consume zero or more leading CommonMark 0.31.2 §4.7 link reference definitions from `text`.

    Returns (normalized_labels_found, number_of_leading_lines_consumed).
    Handles balanced nested parentheses in destinations (depth <= 32), multi-line quoted/parenthesized
    titles, and rejects whitespace-only labels like `[ ]` per CommonMark §4.7.
    """
    n = len(text)
    pos = 0
    labels: set[str] = set()
    lines_consumed = 0

    def _skip_spaces_tabs(p: int) -> int:
        while p < n and text[p] in (" ", "\t"):
            p += 1
        return p

    def _skip_opt_ws_and_one_nl(p: int) -> int:
        p = _skip_spaces_tabs(p)
        if p < n and text[p] == "\n":
            p += 1
            p = _skip_spaces_tabs(p)
        return p

    while pos < n:
        save_pos = pos
        # Up to 3 leading spaces
        sp = 0
        while pos < n and text[pos] == " " and sp < 3:
            pos += 1
            sp += 1
        if pos >= n or text[pos] != "[":
            pos = save_pos
            break
        pos += 1
        lbl_start = pos
        lbl_chars = 0
        found_close = False
        while pos < n and lbl_chars <= 999:
            ch = text[pos]
            if ch == "\\" and pos + 1 < n:
                pos += 2
                lbl_chars += 2
                continue
            if ch == "[":
                break
            if ch == "]":
                found_close = True
                break
            pos += 1
            lbl_chars += 1
        if not found_close:
            pos = save_pos
            break
        raw_label = text[lbl_start:pos]
        # §4.7: label must contain at least one non-whitespace character (ASCII space/tab/newline)
        if not re.search(r"[^ \t\n]", raw_label):
            pos = save_pos
            break
        pos += 1  # consume ']'
        if pos >= n or text[pos] != ":":
            pos = save_pos
            break
        pos += 1  # consume ':'
        pos = _skip_opt_ws_and_one_nl(pos)
        if pos >= n:
            pos = save_pos
            break
        # Parse link destination (§4.7 / §6.3)
        if text[pos] == "<":
            pos += 1
            dest_ok = False
            while pos < n:
                ch = text[pos]
                if ch == "\n" or ch == "<":
                    break
                if ch == "\\" and pos + 1 < n and text[pos + 1] != "\n":
                    pos += 2
                    continue
                if ch == ">":
                    pos += 1
                    dest_ok = True
                    break
                pos += 1
            if not dest_ok:
                pos = save_pos
                break
        else:
            dest_start = pos
            paren_depth = 0
            dest_bad = False
            while pos < n:
                ch = text[pos]
                oc = ord(ch)
                if ch in (" ", "\t", "\n") or oc < 0x20 or oc == 0x7F:
                    break
                if ch == "\\" and pos + 1 < n and text[pos + 1] not in (" ", "\t", "\n"):
                    pos += 2
                    continue
                if ch == "(":
                    paren_depth += 1
                    if paren_depth > 32:
                        dest_bad = True
                        break
                    pos += 1
                    continue
                if ch == ")":
                    if paren_depth == 0:
                        dest_bad = True
                        break
                    paren_depth -= 1
                    pos += 1
                    continue
                pos += 1
            if dest_bad or paren_depth != 0 or pos == dest_start:
                pos = save_pos
                break

        # After destination, check optional title or end-of-line
        after_dest_pos = pos
        ws_after_dest_start = pos
        p_ws = _skip_spaces_tabs(pos)
        had_nl_after_dest = False
        if p_ws < n and text[p_ws] == "\n":
            had_nl_after_dest = True
            p_title_start = _skip_spaces_tabs(p_ws + 1)
        else:
            p_title_start = p_ws
        has_ws_before_title = (p_title_start > ws_after_dest_start)

        def _try_parse_title(tp: int) -> int | None:
            if tp >= n or text[tp] not in ('"', "'", "("):
                return None
            delim = ")" if text[tp] == "(" else text[tp]
            tp += 1
            while tp < n:
                ch = text[tp]
                if ch == "\\" and tp + 1 < n:
                    tp += 2
                    continue
                if delim == ")" and ch == "(":
                    return None
                if ch == "\n":
                    # Blank line inside title is forbidden
                    nl_next = _skip_spaces_tabs(tp + 1)
                    if nl_next < n and text[nl_next] == "\n":
                        return None
                if ch == delim:
                    tp += 1
                    tp = _skip_spaces_tabs(tp)
                    if tp == n:
                        return tp
                    if text[tp] == "\n":
                        return tp + 1
                    return None
                tp += 1
            return None

        title_end = _try_parse_title(p_title_start) if has_ws_before_title else None
        if title_end is not None:
            pos = title_end
            labels.add(_cm_normalize_link_label(raw_label))
            lines_consumed = text[:pos].count("\n") if pos < n else (text.count("\n") + 1)
            continue
        # No valid title; destination must end at EOL (or EOF)
        p_eol = _skip_spaces_tabs(after_dest_pos)
        if p_eol == n:
            pos = n
            labels.add(_cm_normalize_link_label(raw_label))
            lines_consumed = text.count("\n") + 1
            break
        if text[p_eol] == "\n":
            pos = p_eol + 1
            labels.add(_cm_normalize_link_label(raw_label))
            lines_consumed = text[:pos].count("\n")
            continue
        # Neither valid title nor EOL after destination -> not a valid link ref def
        pos = save_pos
        break

    return labels, lines_consumed


_CM_ENTITY_RE = re.compile(r"&(?:#[0-9]{1,7}|#[xX][0-9A-Fa-f]{1,6}|[A-Za-z][A-Za-z0-9]{1,31});", re.ASCII)
_CM_ASCII_PUNCT = frozenset(r"""!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~""")


def _cm_decode_entity(match: re.Match[str]) -> str:
    raw = match.group(0)
    if raw.startswith("&#"):
        try:
            cp = int(raw[3:-1], 16) if raw[2] in ("x", "X") else int(raw[2:-1], 10)
            if cp == 0 or cp > 0x10FFFF or (0xD800 <= cp <= 0xDFFF):
                return "\uFFFD"
            return chr(cp)
        except Exception:
            return raw
    decoded = html.unescape(raw)
    return decoded


def _cm_render_inline_with_meta(raw: str, link_refs: set[str] | None = None) -> tuple[str, bool]:
    """Render visible plain text of a CommonMark 0.31.2 inline string and report whether any link/autolink formed.

    Enforces CommonMark §6.5 (links may not contain other links or autolinks at any level of nesting,
    e.g. Example 518 `[foo [bar](/uri)](/uri)` -> `[foo bar](/uri)`, while images may contain links
    and links may contain images).
    """
    refs = link_refs or set()
    s = raw
    n = len(s)
    # Tokenize into protected literal spans vs plain text so emphasis resolution runs cleanly on plain tokens
    tokens: List[tuple[str, str]] = []  # ("lit" | "text" | "delim", val)
    has_link = False
    i = 0

    def _find_matching_bracket(start_idx: int) -> int | None:
        depth = 0
        j = start_idx
        while j < n:
            ch = s[j]
            if ch == "\\" and j + 1 < n and s[j + 1] in _CM_ASCII_PUNCT:
                j += 2
                continue
            if ch == "`":
                bt = j
                while bt < n and s[bt] == "`":
                    bt += 1
                run = bt - j
                k = bt
                found_bt = False
                while k < n:
                    if s[k] == "`":
                        k2 = k
                        while k2 < n and s[k2] == "`":
                            k2 += 1
                        if k2 - k == run:
                            j = k2
                            found_bt = True
                            break
                        k = k2
                    else:
                        k += 1
                if found_bt:
                    continue
                j = bt
                continue
            if ch == "[":
                depth += 1
            elif ch == "]":
                if depth == 0:
                    return j
                depth -= 1
            j += 1
        return None

    def _try_consume_inline_dest(p: int) -> int | None:
        if p >= n or s[p] != "(":
            return None
        p += 1
        while p < n and s[p] in (" ", "\t", "\n"):
            p += 1
        if p < n and s[p] == ")":
            return p + 1
        if p < n and s[p] == "<":
            p += 1
            while p < n and s[p] != ">":
                if s[p] == "\n":
                    return None
                if s[p] == "\\" and p + 1 < n:
                    p += 2
                    continue
                p += 1
            if p >= n:
                return None
            p += 1
        else:
            pdepth = 0
            while p < n:
                ch = s[p]
                if ch in (" ", "\t", "\n"):
                    break
                if ch == "\\" and p + 1 < n and s[p + 1] in _CM_ASCII_PUNCT:
                    p += 2
                    continue
                if ch == "(":
                    pdepth += 1
                elif ch == ")":
                    if pdepth == 0:
                        break
                    pdepth -= 1
                p += 1
            if pdepth != 0:
                return None
        while p < n and s[p] in (" ", "\t", "\n"):
            p += 1
        if p < n and s[p] in ('"', "'", "("):
            delim = ")" if s[p] == "(" else s[p]
            p += 1
            while p < n and s[p] != delim:
                if s[p] == "\\" and p + 1 < n:
                    p += 2
                    continue
                p += 1
            if p >= n:
                return None
            p += 1
            while p < n and s[p] in (" ", "\t", "\n"):
                p += 1
        if p < n and s[p] == ")":
            return p + 1
        return None

    while i < n:
        ch = s[i]
        # 1. Code span
        if ch == "`":
            j = i
            while j < n and s[j] == "`":
                j += 1
            run = j - i
            k = j
            matched_end = None
            while k < n:
                if s[k] == "`":
                    k2 = k
                    while k2 < n and s[k2] == "`":
                        k2 += 1
                    if k2 - k == run:
                        matched_end = k
                        break
                    k = k2
                else:
                    k += 1
            if matched_end is not None:
                code_body = s[j:matched_end].replace("\n", " ")
                if len(code_body) >= 2 and code_body[0] == " " and code_body[-1] == " " and any(c != " " for c in code_body):
                    code_body = code_body[1:-1]
                tokens.append(("lit", code_body))
                i = matched_end + run
                continue
            tokens.append(("lit", s[i:j]))
            i = j
            continue
        # 2. Backslash escape (CommonMark §2.4: BEFORE HTML entities!)
        if ch == "\\" and i + 1 < n and s[i + 1] in _CM_ASCII_PUNCT:
            tokens.append(("lit", s[i + 1]))
            i += 2
            continue
        # 3. Semicolon-required HTML character reference (CommonMark §2.5)
        if ch == "&":
            me = _CM_ENTITY_RE.match(s, i)
            if me:
                tokens.append(("lit", _cm_decode_entity(me)))
                i = me.end()
                continue
            tokens.append(("text", "&"))
            i += 1
            continue
        # 4. Autolink or raw inline HTML tag / comment
        if ch == "<":
            m_auto = re.match(r"<([A-Za-z][A-Za-z0-9+.-]{1,31}:[^\s<>]*|[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*)>", s[i:])
            if m_auto:
                tokens.append(("lit", m_auto.group(1)))
                has_link = True
                i += m_auto.end()
                continue
            m_html = re.match(
                r"(?:<!--.*?-->|<[A-Za-z][A-Za-z0-9-]*(?:[ \t\n]+[A-Za-z_:][A-Za-z0-9_.:-]*(?:[ \t\n]*=[ \t\n]*(?:[^ \t\n\"'=<>`]+|'[^']*'|\"[^\"]*\"))?)*[ \t\n]*/?>|</[A-Za-z][A-Za-z0-9-]*[ \t\n]*>)",
                s[i:],
                re.S,
            )
            if m_html:
                # Strip raw inline HTML tags from visible text
                i += m_html.end()
                continue
            tokens.append(("text", "<"))
            i += 1
            continue
        # 5. Links & Images: [label](dest), [label][ref], [label][], [label] (shortcut ref)
        if ch == "[" or (ch == "!" and i + 1 < n and s[i + 1] == "["):
            is_img = (ch == "!")
            open_b = i + 1 if is_img else i
            close_b = _find_matching_bracket(open_b + 1)
            if close_b is not None:
                inner_lbl = s[open_b + 1 : close_b]
                next_p = close_b + 1
                end_pos: int | None = None
                # Inline link/image
                end_inline = _try_consume_inline_dest(next_p)
                if end_inline is not None:
                    end_pos = end_inline
                # Full or collapsed reference link [text][ref] / [text][]
                elif next_p < n and s[next_p] == "[":
                    close_r = _find_matching_bracket(next_p + 1)
                    if close_r is not None:
                        ref_raw = s[next_p + 1 : close_r]
                        ref_key = _cm_normalize_link_label(ref_raw if ref_raw else inner_lbl)
                        if ref_key in refs:
                            end_pos = close_r + 1
                # Shortcut reference link [text]
                elif _cm_normalize_link_label(inner_lbl) in refs:
                    end_pos = close_b + 1
                if end_pos is not None:
                    inner_txt, inner_has_link = _cm_render_inline_with_meta(inner_lbl, refs)
                    # CommonMark §6.5: links may not contain other links/autolinks; images may.
                    if is_img or not inner_has_link:
                        tokens.append(("lit", inner_txt))
                        has_link = (has_link or inner_has_link) if is_img else True
                        i = end_pos
                        continue
            tokens.append(("text", ch))
            i += 1
            continue
        # 6. Emphasis delimiters (* or _)
        if ch in ("*", "_"):
            j = i
            while j < n and s[j] == ch:
                j += 1
            tokens.append(("delim", s[i:j]))
            i = j
            continue
        tokens.append(("text", ch))
        i += 1

    # Resolve matched emphasis delimiter runs (* and _) per CommonMark §6.2 flanking rules
    def _prev_char(idx: int) -> str:
        for k in range(idx - 1, -1, -1):
            if tokens[k][1]:
                return tokens[k][1][-1]
        return " "

    def _next_char(idx: int) -> str:
        for k in range(idx + 1, len(tokens)):
            if tokens[k][1]:
                return tokens[k][1][0]
        return " "

    def _is_ws_ch(c: str) -> bool:
        return c.isspace()

    def _is_punct_ch(c: str) -> bool:
        return c in _CM_ASCII_PUNCT

    delims: List[Dict[str, Any]] = []
    for idx, (tk, tv) in enumerate(tokens):
        if tk != "delim":
            continue
        dch = tv[0]
        pc_c = _prev_char(idx)
        nc_c = _next_char(idx)
        left_flank = (not _is_ws_ch(nc_c)) and (not _is_punct_ch(nc_c) or _is_ws_ch(pc_c) or _is_punct_ch(pc_c))
        right_flank = (not _is_ws_ch(pc_c)) and (not _is_punct_ch(pc_c) or _is_ws_ch(nc_c) or _is_punct_ch(nc_c))
        can_open = left_flank if dch == "*" else (left_flank and (not right_flank or _is_punct_ch(pc_c)))
        can_close = right_flank if dch == "*" else (right_flank and (not left_flank or _is_punct_ch(nc_c)))
        delims.append({
            "idx": idx,
            "ch": dch,
            "rem": len(tv),
            "orig": len(tv),
            "can_open": can_open,
            "can_close": can_close,
        })

    # Match closers against nearest valid openers
    ci = 0
    while ci < len(delims):
        closer = delims[ci]
        if closer["rem"] <= 0 or not closer["can_close"]:
            ci += 1
            continue
        matched_opener = None
        for oi in range(ci - 1, -1, -1):
            opener = delims[oi]
            if opener["rem"] <= 0 or not opener["can_open"] or opener["ch"] != closer["ch"]:
                continue
            if (opener["can_close"] or closer["can_open"]) and (opener["orig"] + closer["orig"]) % 3 == 0:
                if opener["orig"] % 3 != 0 or closer["orig"] % 3 != 0:
                    continue
            matched_opener = oi
            break
        if matched_opener is not None:
            opener = delims[matched_opener]
            use = 2 if (opener["rem"] >= 2 and closer["rem"] >= 2) else 1
            opener["rem"] -= use
            closer["rem"] -= use
            # Clear any intermediate unmatched delimiters between opener and closer
            for mi in range(matched_opener + 1, ci):
                delims[mi]["can_open"] = False
                delims[mi]["can_close"] = False
            if closer["rem"] == 0:
                ci += 1
        else:
            ci += 1

    delim_rem = {d["idx"]: d["rem"] for d in delims}
    out_parts: List[str] = []
    for idx, (tk, tv) in enumerate(tokens):
        if tk == "delim":
            rem = delim_rem.get(idx, len(tv))
            if rem > 0:
                out_parts.append(tv[:rem])
        else:
            out_parts.append(tv)
    return "".join(out_parts), has_link


def _cm_render_inline_visible_text(raw: str, link_refs: set[str] | None = None) -> str:
    """Render visible plain text of a CommonMark 0.31.2 inline heading string."""
    txt, _ = _cm_render_inline_with_meta(raw, link_refs)
    return txt


def _parse_commonmark_headings(text: str) -> List[Dict[str, Any]]:
    """Stdlib CommonMark 0.31.2 block + inline heading parser extracting ATX and multiline Setext headings.

    Tracks exact 4-column absolute tab stops across nested blockquotes and list items (§2.2/§5.1/§5.2),
    strictly checks open fenced/HTML leaf blocks before opening new containers, uses ASCII-only blank
    lines, parses full §4.7 link-reference definitions (nested parens, multi-line titles, non-blank labels),
    and renders visible inline heading text (code spans, backslash escapes before semicolon-required HTML5
    entities, links/images, and emphasis).
    """
    norm = text.replace("\r\n", "\n").replace("\r", "\n")
    raw_lines = norm.split("\n")
    line_offsets: List[int] = []
    pos = 0
    for ln in raw_lines:
        line_offsets.append(pos)
        pos += len(ln) + 1

    headings: List[Dict[str, Any]] = []
    doc_link_refs: set[str] = set()
    in_fence = False
    fence_char = ""
    fence_len = 0
    fence_depth = 0
    in_html_block_type = 0
    html_block_depth = 0
    para_lines: List[tuple[int, str, int]] = []

    atx_re = re.compile(r"^ {0,3}(#{1,6})(?:[ \t]+(.*?))?[ \t]*$")
    fence_open_re = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
    setext_re = re.compile(r"^ {0,3}(=+|-+)[ \t]*$")
    thematic_re = re.compile(r"^[ ]{0,3}(?:(?:\*[ \t]*){3,}|(?:-[ \t]*){3,}|(?:_[ \t]*){3,})$")
    ul_re = re.compile(r"^( {0,3})([-+*])(?:([ \t]+)(.*)|([ \t]*))$")
    ol_re = re.compile(r"^( {0,3})([0-9]{1,9}[.)])(?:([ \t]+)(.*)|([ \t]*))$", re.ASCII)
    # CommonMark 0.31.2 §4.6 HTML block start/end rules (Types 1-7)
    html_t1_start_re = re.compile(r"^ {0,3}<(?:pre|script|style|textarea)(?:[ \t>]|$)", re.IGNORECASE)
    html_t1_end_re = re.compile(r"</(?:pre|script|style|textarea)>", re.IGNORECASE)
    html_t2_start_re = re.compile(r"^ {0,3}<!--")
    html_t3_start_re = re.compile(r"^ {0,3}<\?")
    html_t4_start_re = re.compile(r"^ {0,3}<![A-Za-z]")
    html_t5_start_re = re.compile(r"^ {0,3}<!\[CDATA\[")
    html_t6_start_re = re.compile(
        r"^ {0,3}</?(?:address|article|aside|base|basefont|blockquote|body|caption|center|col|colgroup|"
        r"dd|details|dialog|dir|div|dl|dt|fieldset|figcaption|figure|footer|form|frame|frameset|"
        r"h1|h2|h3|h4|h5|h6|head|header|hr|html|iframe|legend|li|link|main|menu|menuitem|nav|"
        r"noframes|ol|optgroup|option|p|param|search|section|summary|table|tbody|td|tfoot|th|"
        r"thead|title|tr|track|ul)(?:[ \t>]|/>|$)",
        re.IGNORECASE,
    )
    html_t7_start_re = re.compile(
        r"^ {0,3}(?:<[A-Za-z][A-Za-z0-9-]*(?:[ \t]+[A-Za-z_:][A-Za-z0-9_.:-]*(?:[ \t]*=[ \t]*(?:[^ \t\"'=<>`]+|'[^']*'|\"[^\"]*\"))?)*[ \t]*/?>|"
        r"</[A-Za-z][A-Za-z0-9-]*[ \t]*>)[ \t]*$"
    )

    container_stack: List[tuple[str, int]] = []

    def _flush_para_for_refs() -> None:
        if not para_lines:
            return
        p_text = "\n".join(p[1] for p in para_lines)
        found_refs, _ = _parse_cm_link_ref_defs_from_text(p_text)
        if found_refs:
            doc_link_refs.update(found_refs)
        para_lines.clear()

    def _count_leading_indent(s: str, start_col: int) -> tuple[int, int]:
        """Return (indent_cols, char_index_of_first_non_ws) with 4-col tab stops from absolute `start_col`."""
        col = start_col
        i = 0
        n_s = len(s)
        while i < n_s:
            ch = s[i]
            if ch == " ":
                col += 1
                i += 1
            elif ch == "\t":
                col += 4 - (col % 4)
                i += 1
            else:
                break
        return col - start_col, i

    def _consume_cols_abs(s: str, cols: int, start_col: int) -> tuple[bool, str, int]:
        """Consume `cols` display columns starting at absolute column `start_col` (CommonMark §2.2)."""
        col = start_col
        target = start_col + cols
        i = 0
        n_s = len(s)
        while i < n_s and col < target:
            ch = s[i]
            if ch == " ":
                col += 1
                i += 1
            elif ch == "\t":
                step = 4 - (col % 4)
                if col + step > target:
                    rem = (col + step) - target
                    return True, (" " * rem) + s[i + 1 :], target
                col += step
                i += 1
            else:
                return False, s, start_col
        if col >= target:
            return True, s[i:], col
        return False, s, start_col

    def _try_match_bq(s: str, start_col: int) -> tuple[bool, str, int]:
        """Match CommonMark §5.1 blockquote marker ` {0,3}>[optional 1-col space/tab]` at absolute `start_col`."""
        ind_cols, idx_ns = _count_leading_indent(s, start_col)
        if ind_cols > 3 or idx_ns >= len(s) or s[idx_ns] != ">":
            return False, s, start_col
        after_gt_col = start_col + ind_cols + 1
        rest = s[idx_ns + 1 :]
        if not rest:
            return True, "", after_gt_col
        if rest[0] == " ":
            return True, rest[1:], after_gt_col + 1
        if rest[0] == "\t":
            # Optional space after '>' consumes 1 column of the tab; remaining columns become spaces!
            tab_w = 4 - (after_gt_col % 4)
            rem_sp = tab_w - 1
            return True, (" " * rem_sp) + rest[1:], after_gt_col + 1
        return True, rest, after_gt_col

    def _try_match_list(s: str, start_col: int, *, interrupting_same_depth_para: bool) -> tuple[bool, int, str, int]:
        """Match CommonMark §5.2 list item marker at absolute `start_col`."""
        ind_cols, idx_ns = _count_leading_indent(s, start_col)
        if ind_cols > 3 or idx_ns >= len(s):
            return False, 0, s, start_col
        sub = s[idx_ns:]
        if thematic_re.match((" " * ind_cols) + sub):
            return False, 0, s, start_col
        m_ul = ul_re.match(sub)
        m_ol = None if m_ul else ol_re.match(sub)
        if not m_ul and not m_ol:
            return False, 0, s, start_col
        m = m_ul or m_ol
        mark_s = m.group(2)
        sep_s = m.group(3)
        rest_s = m.group(4) if sep_s is not None else ""
        mark_end_col = start_col + ind_cols + len(mark_s)
        if sep_s is None:
            if interrupting_same_depth_para:
                return False, 0, s, start_col
            content_indent = ind_cols + len(mark_s) + 1
            return True, content_indent, "", mark_end_col + 1
        is_empty_content = _cm_is_blank(rest_s)
        if interrupting_same_depth_para:
            if is_empty_content:
                return False, 0, s, start_col
            if m_ol is not None and mark_s not in ("1.", "1)"):
                return False, 0, s, start_col
        # Compute exact display width of `sep_s` from `mark_end_col`
        c_cur = mark_end_col
        for ch in sep_s:
            c_cur += (4 - (c_cur % 4)) if ch == "\t" else 1
        sp_cols = c_cur - mark_end_col
        if is_empty_content:
            content_indent = ind_cols + len(mark_s) + 1
            return True, content_indent, "", mark_end_col + 1
        if sp_cols >= 5:
            # Indented code inside list item: marker consumes 1 space, remaining (sp_cols - 1) spaces stay!
            content_indent = ind_cols + len(mark_s) + 1
            return True, content_indent, (" " * (sp_cols - 1)) + rest_s, mark_end_col + 1
        content_indent = ind_cols + len(mark_s) + sp_cols
        return True, content_indent, rest_s, c_cur

    def _interrupts_paragraph(s: str, start_col: int) -> bool:
        ind_cols, idx_ns = _count_leading_indent(s, start_col)
        if ind_cols >= 4:
            return False
        norm_s = (" " * ind_cols) + s[idx_ns:]
        if atx_re.match(norm_s) or thematic_re.match(norm_s):
            return True
        mf_chk = fence_open_re.match(norm_s)
        if mf_chk and not (mf_chk.group(1)[0] == "`" and "`" in mf_chk.group(2)):
            return True
        return bool(
            html_t1_start_re.match(norm_s)
            or html_t2_start_re.match(norm_s)
            or html_t3_start_re.match(norm_s)
            or html_t4_start_re.match(norm_s)
            or html_t5_start_re.match(norm_s)
            or html_t6_start_re.match(norm_s)
        )

    for idx, raw_line in enumerate(raw_lines):
        has_open_para = bool(para_lines)
        open_para_depth = para_lines[-1][2] if para_lines else 0

        # PHASE 1: Match existing open containers in `container_stack` with exact absolute column tracking
        prev_stack = list(container_stack)
        cur = raw_line
        cur_col = 0
        matched_stack: List[tuple[str, int]] = []
        is_blank_line = _cm_is_blank(raw_line)
        for kind, req_indent in container_stack:
            if kind == "bq":
                ok_bq, rest_bq, next_col = _try_match_bq(cur, cur_col)
                if ok_bq:
                    matched_stack.append((kind, req_indent))
                    cur = rest_bq
                    cur_col = next_col
                    continue
                break
            else:
                if _cm_is_blank(cur):
                    matched_stack.append((kind, req_indent))
                    cur = ""
                    cur_col += req_indent
                    continue
                ok_li, rest_li, next_col = _consume_cols_abs(cur, req_indent, cur_col)
                if ok_li:
                    matched_stack.append((kind, req_indent))
                    cur = rest_li
                    cur_col = next_col
                    continue
                break

        matched_depth = len(matched_stack)

        # PHASE 2: Handle open fenced code block or HTML leaf block BEFORE opening ANY new containers!
        if in_html_block_type:
            if matched_depth < html_block_depth:
                in_html_block_type = 0
                html_block_depth = 0
            else:
                container_stack = matched_stack
                if in_html_block_type == 1 and html_t1_end_re.search(cur):
                    in_html_block_type = 0
                    html_block_depth = 0
                elif in_html_block_type == 2 and "-->" in cur:
                    in_html_block_type = 0
                    html_block_depth = 0
                elif in_html_block_type == 3 and "?>" in cur:
                    in_html_block_type = 0
                    html_block_depth = 0
                elif in_html_block_type == 4 and ">" in cur:
                    in_html_block_type = 0
                    html_block_depth = 0
                elif in_html_block_type == 5 and "]]>" in cur:
                    in_html_block_type = 0
                    html_block_depth = 0
                elif in_html_block_type in (6, 7) and _cm_is_blank(cur):
                    in_html_block_type = 0
                    html_block_depth = 0
                _flush_para_for_refs()
                continue

        if in_fence:
            if matched_depth < fence_depth:
                in_fence = False
                fence_char = ""
                fence_len = 0
                fence_depth = 0
            else:
                container_stack = matched_stack
                ind_c, ns_c = _count_leading_indent(cur, cur_col)
                norm_close = (" " * ind_c) + cur[ns_c:] if ind_c <= 3 else ""
                close_m = re.match(rf"^ {{0,3}}{re.escape(fence_char)}{{{fence_len},}}[ \t]*$", norm_close) if norm_close else None
                if close_m:
                    in_fence = False
                    fence_char = ""
                    fence_len = 0
                    fence_depth = 0
                _flush_para_for_refs()
                continue

        if is_blank_line:
            while container_stack and container_stack[-1][0] == "bq" and len(container_stack) > matched_depth:
                container_stack.pop()
            _flush_para_for_refs()
            continue

        # PHASE 3: Lazy paragraph continuation check (CommonMark §5.1 / §4.3)
        if has_open_para and matched_depth < open_para_depth:
            ok_bq_chk, _, _ = _try_match_bq(cur, cur_col)
            opens_new_container = ok_bq_chk
            if not opens_new_container:
                # At a lower depth than the open paragraph, a list marker starts a NEW container
                # (any ordered start number with non-empty content, or empty list item at depth 0 per Ex. 278)
                ok_li_chk, _, rest_li_chk, _ = _try_match_list(cur, cur_col, interrupting_same_depth_para=False)
                if ok_li_chk:
                    if not _cm_is_blank(rest_li_chk) or matched_depth == 0:
                        opens_new_container = True
            if not opens_new_container and not _interrupts_paragraph(cur, cur_col):
                container_stack = prev_stack
                para_lines.append((idx, _cm_strip_ascii_ws(cur), open_para_depth))
                continue

        container_stack = matched_stack

        # PHASE 4: Open any new nested blockquotes or list items on this line
        while True:
            ok_bq, rest_bq, next_col = _try_match_bq(cur, cur_col)
            if ok_bq:
                container_stack.append(("bq", 0))
                cur = rest_bq
                cur_col = next_col
                continue
            ind_c, ns_c = _count_leading_indent(cur, cur_col)
            norm_lead = (" " * ind_c) + cur[ns_c:] if ind_c <= 3 else cur
            if has_open_para and len(container_stack) == open_para_depth and ind_c <= 3 and setext_re.match(norm_lead):
                break
            if ind_c <= 3 and thematic_re.match(norm_lead):
                break
            same_depth_para = bool(has_open_para and len(container_stack) == open_para_depth)
            ok_li, req_ind, rest_li, next_col_li = _try_match_list(cur, cur_col, interrupting_same_depth_para=same_depth_para)
            if ok_li:
                container_stack.append(("list", req_ind))
                cur = rest_li
                cur_col = next_col_li
                continue
            break

        depth = len(container_stack)
        if _cm_is_blank(cur):
            _flush_para_for_refs()
            continue

        ind_cols, idx_ns = _count_leading_indent(cur, cur_col)
        if not para_lines and ind_cols >= 4:
            continue
        inner = (" " * ind_cols) + cur[idx_ns:] if ind_cols <= 3 else cur[idx_ns:]

        if ind_cols <= 3:
            if html_t1_start_re.match(inner):
                if not html_t1_end_re.search(inner):
                    in_html_block_type = 1
                    html_block_depth = depth
                _flush_para_for_refs()
                continue
            if html_t2_start_re.match(inner):
                if "-->" not in inner[inner.find("<!--") + 4 :]:
                    in_html_block_type = 2
                    html_block_depth = depth
                _flush_para_for_refs()
                continue
            if html_t3_start_re.match(inner):
                if "?>" not in inner[inner.find("<?") + 2 :]:
                    in_html_block_type = 3
                    html_block_depth = depth
                _flush_para_for_refs()
                continue
            if html_t4_start_re.match(inner):
                if ">" not in inner[inner.find("<!") + 2 :]:
                    in_html_block_type = 4
                    html_block_depth = depth
                _flush_para_for_refs()
                continue
            if html_t5_start_re.match(inner):
                if "]]>" not in inner[inner.find("<![CDATA[") + 9 :]:
                    in_html_block_type = 5
                    html_block_depth = depth
                _flush_para_for_refs()
                continue
            if html_t6_start_re.match(inner):
                in_html_block_type = 6
                html_block_depth = depth
                _flush_para_for_refs()
                continue
            if not para_lines and html_t7_start_re.match(inner):
                in_html_block_type = 7
                html_block_depth = depth
                _flush_para_for_refs()
                continue

            mf = fence_open_re.match(inner)
            if mf:
                f_tok, info = mf.group(1), mf.group(2)
                if not (f_tok[0] == "`" and "`" in info):
                    in_fence = True
                    fence_char = f_tok[0]
                    fence_len = len(f_tok)
                    fence_depth = depth
                    _flush_para_for_refs()
                    continue

            if para_lines and para_lines[-1][2] == depth:
                ms = setext_re.match(inner)
                if ms:
                    p_text = "\n".join(p[1] for p in para_lines)
                    found_refs, consumed_n = _parse_cm_link_ref_defs_from_text(p_text)
                    if found_refs:
                        doc_link_refs.update(found_refs)
                    _rem_para = para_lines[consumed_n:]
                    if _rem_para:
                        level = 1 if ms.group(1)[0] == "=" else 2
                        first_idx = _rem_para[0][0]
                        raw_title_joined = "\n".join(p[1] for p in _rem_para)
                        headings.append({
                            "kind": "setext",
                            "level": level,
                            "container_depth": depth,
                            "start_line": first_idx,
                            "end_line": idx,
                            "start_offset": line_offsets[first_idx],
                            "raw_first_line": raw_lines[first_idx],
                            "raw_inline_title": raw_title_joined,
                        })
                    para_lines.clear()
                    continue

            ma = atx_re.match(inner)
            if ma:
                _flush_para_for_refs()
                level = len(ma.group(1))
                title_raw = _cm_strip_ascii_ws(ma.group(2) or "")
                title_raw = _cm_strip_ascii_ws(re.sub(r"(?:^|[ \t]+)#+[ \t]*$", "", title_raw))
                headings.append({
                    "kind": "atx",
                    "level": level,
                    "container_depth": depth,
                    "start_line": idx,
                    "end_line": idx,
                    "start_offset": line_offsets[idx],
                    "raw_first_line": raw_line,
                    "raw_inline_title": title_raw,
                })
                continue

            if thematic_re.match(inner):
                _flush_para_for_refs()
                continue

        if para_lines and para_lines[-1][2] != depth:
            _flush_para_for_refs()
        para_lines.append((idx, _cm_strip_ascii_ws(cur), depth))

    _flush_para_for_refs()

    # Second pass: render visible inline heading titles with the complete document link-reference map
    for h in headings:
        raw_t = h.pop("raw_inline_title")
        h["rendered_title"] = _cm_render_inline_visible_text(raw_t, doc_link_refs)
    return headings


def _exact_pct_to_frac(pct: Decimal) -> Decimal | None:
    """Divide a finite Decimal percentage by 100 with zero context-precision rounding; None if exponent underflows Decimal."""
    if not isinstance(pct, Decimal) or not pct.is_finite():
        return None
    tup = pct.as_tuple()
    new_exp = int(tup.exponent) - 2
    if new_exp < decimal.MIN_ETINY:
        return None
    try:
        return Decimal((tup.sign, tup.digits, new_exp))
    except Exception:
        return None


def _dec_frac_bit_estimate(d: Decimal) -> int:
    """Upper-bound estimate of exact rational bit-length for Fraction(d) BEFORE constructing Fraction."""
    tup = d.as_tuple()
    exp = int(tup.exponent)
    return len(tup.digits) * 4 + abs(exp) * 4 + 16


def _abs_diff_exceeds_tol(a: Decimal, b: Decimal, tol: Decimal) -> bool | None:
    """Context-isolated exact comparison for |a - b| > tol (immune to caller Decimal precision/rounding/traps).

    Uses exact Fraction arithmetic whenever the aligned decimal exponent span is within a bounded budget
    (<= 10,000 digits), falling back to a 120-digit directed outward interval enclosure (returning None
    if the mismatch straddles tol beyond the exact budget).
    """
    if not (a.is_finite() and b.is_finite() and tol.is_finite()):
        return None
    if a == b:
        return False
    ea = int(a.as_tuple().exponent)
    eb = int(b.as_tuple().exponent)
    et = int(tol.as_tuple().exponent)
    span = max(abs(ea), abs(eb), abs(et)) + len(a.as_tuple().digits) + len(b.as_tuple().digits) + len(tol.as_tuple().digits)
    if span <= 10_000:
        return abs(Fraction(a) - Fraction(b)) > Fraction(tol)
    with _clean_dec_ctx(120, ROUND_FLOOR) as ctx_d:
        hi_v, lo_v = (a, b) if a >= b else (b, a)
        diff_lo = hi_v - lo_v
        ctx_d.rounding = ROUND_CEILING
        diff_hi = hi_v - lo_v
    if diff_lo > tol:
        return True
    if diff_hi <= tol:
        return False
    return None


def _sound_one_minus_exp_neg(x_low: Decimal, x_high: Decimal, prec: int) -> tuple[Decimal, Decimal]:
    """Provable outward interval [p_lo, p_hi] enclosing 1 - exp(-x) for any x in [x_low, x_high] > 0.

    For x_high < 1e-6, iterates strictly in paired terms (2m-1, 2m) so the lower bound ALWAYS terminates
    at an even partial sum S_{2m}(x_low) and the upper bound ALWAYS terminates at an odd partial sum
    S_{2m+1}(x_high), guaranteeing p_lo <= 1 - exp(-x) <= p_hi at every exit and any precision.
    """
    wprec = max(prec + 40, 140)
    if x_high < Decimal("1e-6"):
        max_pairs = max(64, prec)
        with _clean_dec_ctx(wprec, ROUND_FLOOR) as ctx_t:
            # 1. Provable lower bound on 1 - exp(-x_low) via even truncation S_{2m}(x_low)
            thresh_lo = x_low * Decimal(1).scaleb(-(prec + 8))
            fact = 1
            pow_lo = x_low
            pow_hi = x_low
            s_lo = Decimal(0)
            for pair_idx in range(1, max_pairs + 1):
                # Odd term k = 2*pair_idx - 1 (added with ROUND_FLOOR)
                ctx_t.rounding = ROUND_FLOOR
                term_odd_lo = pow_lo / Decimal(fact)
                s_lo = s_lo + term_odd_lo
                # Advance to even term k = 2*pair_idx
                k_even = 2 * pair_idx
                fact *= k_even
                ctx_t.rounding = ROUND_FLOOR
                pow_lo = pow_lo * x_low
                ctx_t.rounding = ROUND_CEILING
                pow_hi = pow_hi * x_low
                # Even term subtracted with ROUND_CEILING on magnitude & ROUND_FLOOR on sum
                ctx_t.rounding = ROUND_CEILING
                term_even_hi = pow_hi / Decimal(fact)
                ctx_t.rounding = ROUND_FLOOR
                s_lo = s_lo - term_even_hi
                if term_even_hi <= thresh_lo:
                    break
                fact *= k_even + 1
                ctx_t.rounding = ROUND_FLOOR
                pow_lo = pow_lo * x_low
                ctx_t.rounding = ROUND_CEILING
                pow_hi = pow_hi * x_low
            if s_lo < Decimal(0):
                s_lo = Decimal(0)

            # 2. Provable upper bound on 1 - exp(-x_high) via odd truncation S_{2m+1}(x_high)
            thresh_hi = x_high * Decimal(1).scaleb(-(prec + 8))
            fact = 1
            pow_lo = x_high
            pow_hi = x_high
            # Start with k = 1 (odd upper bound S_1 = x_high)
            ctx_t.rounding = ROUND_CEILING
            s_hi = x_high
            for pair_idx in range(1, max_pairs + 1):
                # Even term k = 2*pair_idx (subtracted with ROUND_FLOOR on magnitude & ROUND_CEILING on sum)
                k_even = 2 * pair_idx
                fact *= k_even
                ctx_t.rounding = ROUND_FLOOR
                pow_lo = pow_lo * x_high
                ctx_t.rounding = ROUND_CEILING
                pow_hi = pow_hi * x_high
                ctx_t.rounding = ROUND_FLOOR
                term_even_lo = pow_lo / Decimal(fact)
                ctx_t.rounding = ROUND_CEILING
                s_hi = s_hi - term_even_lo
                # Odd term k = 2*pair_idx + 1 (added with ROUND_CEILING -> guaranteed odd upper bound S_{2m+1})
                k_odd = k_even + 1
                fact *= k_odd
                ctx_t.rounding = ROUND_FLOOR
                pow_lo = pow_lo * x_high
                ctx_t.rounding = ROUND_CEILING
                pow_hi = pow_hi * x_high
                ctx_t.rounding = ROUND_CEILING
                term_odd_hi = pow_hi / Decimal(fact)
                s_hi = s_hi + term_odd_hi
                if term_odd_hi <= thresh_hi:
                    break
            if s_hi > Decimal(1):
                s_hi = Decimal(1)
        with _clean_dec_ctx(prec, ROUND_FLOOR):
            p_lo_out = +s_lo
        with _clean_dec_ctx(prec, ROUND_CEILING):
            p_hi_out = +s_hi
        return p_lo_out, p_hi_out

    with _clean_dec_ctx(wprec, ROUND_FLOOR) as ctx_e:
        exp_lo = (-x_high).exp().next_minus()
        ctx_e.rounding = ROUND_CEILING
        p_hi_raw = Decimal(1) - exp_lo
        exp_hi = (-x_low).exp().next_plus()
        ctx_e.rounding = ROUND_FLOOR
        p_lo_raw = Decimal(1) - exp_hi
        if p_lo_raw < Decimal(0):
            p_lo_raw = Decimal(0)
        if p_hi_raw > Decimal(1):
            p_hi_raw = Decimal(1)
    with _clean_dec_ctx(prec, ROUND_FLOOR):
        p_lo_out = +p_lo_raw
    with _clean_dec_ctx(prec, ROUND_CEILING):
        p_hi_out = +p_hi_raw
    return p_lo_out, p_hi_out


def _zero_event_u95_compare(
    alpha_d: Decimal,
    n_raw: int,
    f_d: Decimal,
    eta_d: Decimal,
    target_frac: Decimal,
    *,
    u95_lower: Decimal,
    u95_upper: Decimal,
) -> bool | None:
    """Sound three-valued decider for true U95 = (1 - alpha^(1/n)) / (f * eta) <= target_frac.

    Returns:
      - True:  provably true_u95 <= target_frac
      - False: provably target_frac < true_u95
      - None:  beyond the 2,000,000-bit exact rational budget AND target_frac lies strictly inside
               the ~1e-258-relative-width 260-digit Taylor/exp outward enclosure [u_lo_260, u_hi_260].
    """
    if target_frac <= Decimal(0):
        return False
    if u95_upper <= target_frac:
        return True
    if target_frac < u95_lower:
        return False
    # Budget prospective Fraction conversions BEFORE constructing any Fraction (zero DoS on extreme exponents),
    # and compute 1 - f*eta*target directly in exact rational arithmetic (never rounding near-one bases!).
    max_in_bits = max(
        _dec_frac_bit_estimate(alpha_d),
        _dec_frac_bit_estimate(f_d),
        _dec_frac_bit_estimate(eta_d),
        _dec_frac_bit_estimate(target_frac),
    )
    if max_in_bits <= 250_000:
        prod_frac = Fraction(f_d) * Fraction(eta_d) * Fraction(target_frac)
        b_frac = Fraction(1) - prod_frac
        if b_frac <= 0:
            return True
        a_frac = Fraction(alpha_d)
        max_base_bits = max(b_frac.numerator.bit_length(), b_frac.denominator.bit_length())
        if n_raw * max_base_bits <= 2_000_000:
            lhs_num = pow(b_frac.numerator, n_raw) * a_frac.denominator
            rhs_num = a_frac.numerator * pow(b_frac.denominator, n_raw)
            return lhs_num <= rhs_num
    with _clean_dec_ctx(280, ROUND_FLOOR) as ctx_hi:
        eff_lo_hi = f_d * eta_d
        if eff_lo_hi <= Decimal(0):
            return False
        ctx_hi.rounding = ROUND_CEILING
        eff_hi_hi = f_d * eta_d
        ctx_hi.rounding = ROUND_FLOOR
        ln_a_lo = alpha_d.ln().next_minus()
        ctx_hi.rounding = ROUND_CEILING
        ln_a_hi = alpha_d.ln().next_plus()
        n_dec = Decimal(n_raw)
        ctx_hi.rounding = ROUND_FLOOR
        x_lo_260 = (-ln_a_hi) / n_dec
        ctx_hi.rounding = ROUND_CEILING
        x_hi_260 = (-ln_a_lo) / n_dec
        p_lo_260, p_hi_260 = _sound_one_minus_exp_neg(x_lo_260, x_hi_260, 260)
        ctx_hi.rounding = ROUND_FLOOR
        u_lo_260 = p_lo_260 / eff_hi_hi
        ctx_hi.rounding = ROUND_CEILING
        u_hi_260 = p_hi_260 / eff_lo_hi
        if u_lo_260 > Decimal(1):
            u_lo_260 = Decimal(1)
        if u_hi_260 > Decimal(1):
            u_hi_260 = Decimal(1)
    if u_hi_260 <= target_frac:
        return True
    if target_frac < u_lo_260:
        return False
    return None


def _zero_event_u95_le_target(
    alpha_d: Decimal,
    n_raw: int,
    f_d: Decimal,
    eta_d: Decimal,
    target_frac: Decimal,
    *,
    u95_lower: Decimal,
    u95_upper: Decimal,
) -> bool:
    """Fail-closed conservative decider for true U95 <= target_frac (returns True iff provably <= target_frac)."""
    return _zero_event_u95_compare(
        alpha_d, n_raw, f_d, eta_d, target_frac, u95_lower=u95_lower, u95_upper=u95_upper
    ) is True


def evaluate_sound_decimal_clopper_pearson_zero_event(
    *,
    gate_id: str,
    k: Any,
    n_raw: Any,
    max_u95_fraction: Any,
    eta_recovery_l95: Any = Decimal(1),
    f_eligible: Any = Decimal(1),
    f_eligible_provenance: str = "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
    alpha_obs: Any = None,
    alpha_eta_recovery: Any = None,
    alpha_f_eligible: Any = None,
    sampling_independence_model: Any = "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES",
    clustered_sampling: Any = False,
    sampling_independence_conflict: bool = False,
    sampling_malformed: bool = False,
    sampling_independence_provenance: str = "UNSPECIFIED_DEFAULT_ASSUMED_INDEPENDENT",
    calibration_hold_reason: str | None = None,
    calibration_sources: Dict[str, str] | None = None,
    reported_u95_pct: Decimal | None = None,
) -> Dict[str, Any]:
    """Strictly sound outward Decimal rare-event evaluator for zero-event Clopper-Pearson bounds.

    Semantics (pre-empting ambiguous 'k==0 veto' shorthand):
      - Any integer k >= 1 IMMEDIATELY returns FAIL_HARD_STOP (unconditional hard stop).
      - When k == 0, bounds the zero-event upper limit by directed Decimal arithmetic using logarithm
        and exponential enclosures at ordinary arguments and a paired alternating Taylor enclosure at
        small arguments (x < 1e-6). The final upper bound is >= the exact capped upper limit (clamped
        inclusively at 1). Comparisons that remain unresolved within the 2,000,000-bit exact rational
        and 260-digit outward enclosure proof budget HOLD explicitly rather than misclassifying.
      - Whenever eta_recovery_l95 < 1.0, requires alpha_eta_recovery >= 0.05 and Bonferroni sum <= 0.10.
      - Always stamps empirical_clearance_status as mathematical count subgate arithmetic ONLY,
        conditional on validated recovery/eligible-fraction calibration, supplying no empirical clearance.
    """
    cap_d = _safe_exact_decimal(max_u95_fraction)
    cap_valid = cap_d is not None and Decimal(0) < cap_d <= Decimal(1)
    with _clean_dec_ctx(80, ROUND_HALF_EVEN):
        cap_str = format(cap_d, ".16e") if cap_valid and cap_d is not None else None

    is_int_k = isinstance(k, int) and not isinstance(k, bool)
    is_int_n = isinstance(n_raw, int) and not isinstance(n_raw, bool)

    # 1. Explicit k >= 1 is ALWAYS FAIL_HARD_STOP (independent of denominator, calibration, or cap validity)
    if is_int_k and k >= 1:
        return {
            "gate_id": gate_id,
            "status": "FAIL_HARD_STOP",
            "reason": "POSITIVE_EVENT_COUNT_K_GE_1_HARD_STOP",
            "empirical_clearance_status": "FAIL_HARD_STOP_POSITIVE_EVENT_DETECTED",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "k": k,
            "n_raw": n_raw if is_int_n else None,
            "cap_fraction_str": cap_str,
        }

    # 2. Validate cap, integer k == 0, and n_raw >= 1
    if not cap_valid or cap_d is None:
        return {
            "gate_id": gate_id,
            "status": "HOLD_INVALID_CAP_FRACTION",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if k is None or n_raw is None:
        return {
            "gate_id": gate_id,
            "status": "HOLD_MISSING_ATOMIC_INTEGER_K_OR_N",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INCOMPLETE",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if not is_int_k or not is_int_n:
        return {
            "gate_id": gate_id,
            "status": "HOLD_NON_INTEGER_COUNT_OR_DENOMINATOR_REJECTED",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if n_raw <= 0 or k < 0 or k > n_raw:
        return {
            "gate_id": gate_id,
            "status": "HOLD_INVALID_COUNT_BOUNDS",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }

    # 2b. Reject contradictory, malformed, or explicitly clustered/non-independent sampling units (GY-018 / DEC-103 parity)
    if sampling_independence_conflict:
        return {
            "gate_id": gate_id,
            "status": "HOLD_CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS",
            "reason": "CONTRADICTORY_PANEL_AND_ASSAY_SAMPLING_INDEPENDENCE_DECLARATIONS",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if sampling_malformed:
        return {
            "gate_id": gate_id,
            "status": "HOLD_MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION",
            "reason": "MALFORMED_OR_NON_ASCII_SAMPLING_INDEPENDENCE_DECLARATION",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if (
        not isinstance(sampling_independence_model, str)
        or sampling_independence_model != "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"
        or clustered_sampling is not False
    ):
        return {
            "gate_id": gate_id,
            "status": "HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED",
            "reason": "REQUIRES_INDEPENDENT_BIOLOGICAL_OPPORTUNITIES_AND_CLUSTERED_SAMPLING_FALSE",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }

    # 2c. Reject cross-scope calibration conflicts, explicit-null overrides, or split bound groups
    if calibration_hold_reason is not None:
        return {
            "gate_id": gate_id,
            "status": calibration_hold_reason,
            "reason": calibration_hold_reason,
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }

    # 3. Validate calibration factors & provenance (alpha_eta >= 0.05 whenever eta < 1.0)
    eta_d = _safe_exact_decimal(eta_recovery_l95)
    f_d = _safe_exact_decimal(f_eligible)
    if eta_d is None or f_d is None:
        return {
            "gate_id": gate_id,
            "status": "HOLD_NONFINITE_RECOVERY_OR_ELIGIBLE_FACTOR",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if not (Decimal(0) < eta_d <= Decimal(1)) or not (Decimal(0) < f_d <= Decimal(1)):
        return {
            "gate_id": gate_id,
            "status": "HOLD_OUT_OF_RANGE_RECOVERY_OR_ELIGIBLE_FACTOR",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if (eta_d < Decimal(1) or f_d < Decimal(1)) and eta_d < Decimal("0.30"):
        return {
            "gate_id": gate_id,
            "status": "HOLD_RECOVERY_BELOW_0P30_QUALIFICATION_FLOOR",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }
    if not isinstance(f_eligible_provenance, str) or f_eligible_provenance not in VALID_F_ELIGIBLE_PROVENANCE:
        return {
            "gate_id": gate_id,
            "status": "HOLD_MISSING_OR_INVALID_ELIGIBLE_PROVENANCE",
            "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
            "sampling_independence_provenance": sampling_independence_provenance,
            "calibration_sources": calibration_sources,
            "cap_fraction_str": cap_str,
        }

    # Validate reported_u95_pct finiteness BEFORE arithmetic if supplied
    rep_pct_d: Decimal | None = None
    if reported_u95_pct is not None:
        rep_pct_d = _safe_exact_decimal(reported_u95_pct)
        if rep_pct_d is None:
            return {
                "gate_id": gate_id,
                "status": "HOLD_INVALID_REPORTED_U95_PCT",
                "reason": "NONFINITE_OR_MALFORMED_REPORTED_U95_PCT",
                "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                "sampling_independence_provenance": sampling_independence_provenance,
                "calibration_sources": calibration_sources,
                "cap_fraction_str": cap_str,
            }

    exact_joint_alpha_frac: Fraction | None = None
    with _clean_dec_ctx(80, ROUND_CEILING) as ctx:
        if f_eligible_provenance == "KNOWN_OR_CONSERVATIVELY_ESTABLISHED":
            a_obs_d = Decimal("0.05") if alpha_obs is None else _safe_exact_decimal(alpha_obs)
            a_eta_d = (
                (Decimal("0.05") if eta_d < Decimal(1) else Decimal("0"))
                if alpha_eta_recovery is None
                else _safe_exact_decimal(alpha_eta_recovery)
            )
            if alpha_f_eligible is not None:
                _af_chk = _safe_exact_decimal(alpha_f_eligible)
                if _af_chk is None or _af_chk != Decimal("0"):
                    return {
                        "gate_id": gate_id,
                        "status": "HOLD_CONFLICTING_F_ELIGIBLE_PROVENANCE_AND_ALPHA",
                        "reason": "NONZERO_OR_INVALID_ALPHA_F_ELIGIBLE_INCOMPATIBLE_WITH_KNOWN_OR_CONSERVATIVELY_ESTABLISHED_PROVENANCE",
                        "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                        "sampling_independence_provenance": sampling_independence_provenance,
                        "calibration_sources": calibration_sources,
                        "cap_fraction_str": cap_str,
                    }
            a_f_d = Decimal("0")
            if (
                a_obs_d is None
                or a_eta_d is None
                or not (Decimal(0) < a_obs_d <= Decimal("0.05"))
                or not (Decimal(0) <= a_eta_d <= Decimal("0.10"))
            ):
                return {
                    "gate_id": gate_id,
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                    "sampling_independence_provenance": sampling_independence_provenance,
                    "calibration_sources": calibration_sources,
                    "cap_fraction_str": cap_str,
                }
            if eta_d < Decimal(1) and a_eta_d < Decimal("0.05"):
                return {
                    "gate_id": gate_id,
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "reason": "ETA_RECOVERY_L95_REQUIRES_ALPHA_ETA_AT_LEAST_0P05",
                    "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                    "sampling_independence_provenance": sampling_independence_provenance,
                    "calibration_sources": calibration_sources,
                    "cap_fraction_str": cap_str,
                }
        else:
            a_obs_d = _safe_exact_decimal(alpha_obs)
            a_f_d = _safe_exact_decimal(alpha_f_eligible)
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
                or not (Decimal(0) <= a_eta_d <= Decimal("0.10"))
            ):
                return {
                    "gate_id": gate_id,
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                    "sampling_independence_provenance": sampling_independence_provenance,
                    "calibration_sources": calibration_sources,
                    "cap_fraction_str": cap_str,
                }
            if eta_d < Decimal(1) and a_eta_d < Decimal("0.05"):
                return {
                    "gate_id": gate_id,
                    "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                    "reason": "ETA_RECOVERY_L95_REQUIRES_ALPHA_ETA_AT_LEAST_0P05",
                    "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                    "sampling_independence_provenance": sampling_independence_provenance,
                    "calibration_sources": calibration_sources,
                    "cap_fraction_str": cap_str,
                }

        # Never award a calibrated bound or substitute unit default 1.0 when caller allocates alpha or declares an estimated bound!
        if a_eta_d > Decimal(0) and (
            eta_d >= Decimal(1)
            or (isinstance(calibration_sources, dict) and calibration_sources.get("eta_recovery_l95") == "default")
        ):
            return {
                "gate_id": gate_id,
                "status": "HOLD_MISSING_OR_UNIT_ETA_RECOVERY_FOR_ALPHA_ALLOCATION",
                "reason": "ALPHA_ETA_RECOVERY_ALLOCATION_REQUIRES_EXPLICIT_ESTIMATE_LT_1",
                "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                "sampling_independence_provenance": sampling_independence_provenance,
                "calibration_sources": calibration_sources,
                "cap_fraction_str": cap_str,
            }
        if (a_f_d > Decimal(0) or f_eligible_provenance == "INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION") and (
            f_d >= Decimal(1)
            or (isinstance(calibration_sources, dict) and calibration_sources.get("f_eligible") == "default")
        ):
            return {
                "gate_id": gate_id,
                "status": "HOLD_MISSING_OR_UNIT_F_ELIGIBLE_FOR_ALPHA_ALLOCATION",
                "reason": "ESTIMATED_F_ELIGIBLE_PROVENANCE_OR_ALPHA_REQUIRES_EXPLICIT_ESTIMATE_LT_1",
                "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                "sampling_independence_provenance": sampling_independence_provenance,
                "calibration_sources": calibration_sources,
                "cap_fraction_str": cap_str,
            }

        # Directed interval + resource-bounded exact aligned-integer Bonferroni sum check (zero false-reject on exact 0.10!)
        with _clean_dec_ctx(80, ROUND_FLOOR):
            _alpha_sum_lower = a_obs_d + a_eta_d + a_f_d
        with _clean_dec_ctx(80, ROUND_CEILING):
            _alpha_sum_upper = a_obs_d + a_eta_d + a_f_d
        if _alpha_sum_lower > Decimal("0.10"):
            return {
                "gate_id": gate_id,
                "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                "reason": "BONFERRONI_ALPHA_SUM_EXCEEDS_0P10_FOR_90PCT_JOINT_COVERAGE",
                "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                "sampling_independence_provenance": sampling_independence_provenance,
                "calibration_sources": calibration_sources,
                "cap_fraction_str": cap_str,
            }
        if _alpha_sum_upper > Decimal("0.10"):
            # Bracket straddles 0.10: resolve via exact bounded Fraction sum before deciding
            _alpha_bits = max(
                _dec_frac_bit_estimate(a_obs_d),
                _dec_frac_bit_estimate(a_eta_d),
                _dec_frac_bit_estimate(a_f_d),
            )
            if _alpha_bits <= 40_000:
                exact_joint_alpha_frac = Fraction(a_obs_d) + Fraction(a_eta_d) + Fraction(a_f_d)
                if exact_joint_alpha_frac > Fraction(1, 10):
                    return {
                        "gate_id": gate_id,
                        "status": "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
                        "reason": "BONFERRONI_ALPHA_SUM_EXCEEDS_0P10_FOR_90PCT_JOINT_COVERAGE",
                        "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                        "sampling_independence_provenance": sampling_independence_provenance,
                        "calibration_sources": calibration_sources,
                        "cap_fraction_str": cap_str,
                    }
            else:
                return {
                    "gate_id": gate_id,
                    "status": "HOLD_BONFERRONI_ALPHA_SUM_UNRESOLVED_WITHIN_BUDGET",
                    "reason": "BONFERRONI_ALPHA_SUM_STRADDLES_0P10_BEYOND_EXACT_BUDGET",
                    "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                    "sampling_independence_provenance": sampling_independence_provenance,
                    "calibration_sources": calibration_sources,
                    "cap_fraction_str": cap_str,
                }

    # 4. Exact 80-digit sound outward Decimal upper & lower bounds (fresh trap-free Context; paired Taylor -expm1 for tiny x)
    with _clean_dec_ctx(80, ROUND_FLOOR) as ctx:
        eff_lower = f_d * eta_d
        if eff_lower <= Decimal(0):
            return {
                "gate_id": gate_id,
                "status": "HOLD_UNDERFLOW_OR_INVALID_DETECTION_FACTOR",
                "empirical_clearance_status": "HOLD_MATHEMATICAL_SUBGATE_INVALID_INPUT",
                "sampling_independence_provenance": sampling_independence_provenance,
                "calibration_sources": calibration_sources,
                "cap_fraction_str": cap_str,
            }
        ctx.rounding = ROUND_CEILING
        eff_upper = f_d * eta_d
        ctx.rounding = ROUND_FLOOR
        ln_alpha_lo = a_obs_d.ln().next_minus()
        ctx.rounding = ROUND_CEILING
        ln_alpha_hi = a_obs_d.ln().next_plus()
        n_d = Decimal(n_raw)
        ctx.rounding = ROUND_FLOOR
        x_lo = (-ln_alpha_hi) / n_d
        ctx.rounding = ROUND_CEILING
        x_hi = (-ln_alpha_lo) / n_d
        p_raw_lower, p_raw_upper = _sound_one_minus_exp_neg(x_lo, x_hi, 80)
        ctx.rounding = ROUND_CEILING
        u95_upper = p_raw_upper / eff_lower
        if u95_upper > Decimal(1):
            u95_upper = Decimal(1)
        ctx.rounding = ROUND_FLOOR
        u95_lower = p_raw_lower / eff_upper
        if u95_lower > Decimal(1):
            u95_lower = Decimal(1)
        u95_lower_pct = u95_lower * Decimal(100)
        u95_math_lower_pct_str = format(u95_lower_pct, "e") if (u95_lower_pct > Decimal(0) and u95_lower_pct < Decimal("1e-6")) else format(u95_lower_pct, "f")
        ctx.rounding = ROUND_CEILING
        u95_upper_pct_exact = u95_upper * Decimal(100)
        u95_math_upper_pct_exact_str = format(u95_upper_pct_exact, "e") if (u95_upper_pct_exact > Decimal(0) and u95_upper_pct_exact < Decimal("1e-6")) else format(u95_upper_pct_exact, "f")

        # Sound upward-rounded (ROUND_CEILING) serialization so exported strings NEVER round below u95_upper
        if u95_upper == Decimal(0):
            u95_str = "0.0000000000000000e+00"
        else:
            exp_adj = u95_upper.adjusted()
            mant_ceil = u95_upper.scaleb(-exp_adj).quantize(Decimal("1e-16"), rounding=ROUND_CEILING)
            if mant_ceil >= Decimal(10):
                mant_ceil = (mant_ceil / Decimal(10)).quantize(Decimal("1e-16"), rounding=ROUND_CEILING)
                exp_adj += 1
            u95_str = f"{mant_ceil:.16f}e{exp_adj:+03d}"
        u95_pct_ceil = (u95_upper * Decimal(100)).quantize(Decimal("0.00000001"), rounding=ROUND_CEILING)
        u95_pct_str = f"{u95_pct_ceil:.8f}"

    cap_cmp = _zero_event_u95_compare(
        a_obs_d, n_raw, f_d, eta_d, cap_d, u95_lower=u95_lower, u95_upper=u95_upper
    )
    passed = (cap_cmp is True)
    rep_meets_exact: bool | None = None
    if rep_pct_d is not None:
        if rep_pct_d <= Decimal(0):
            rep_meets_exact = False
        elif rep_pct_d < u95_lower_pct:
            # Compare in percent units FIRST so extreme subnormal exponents like MIN_ETINY never throw!
            rep_meets_exact = False
        elif u95_upper_pct_exact <= rep_pct_d:
            rep_meets_exact = True
        else:
            rep_frac_d = _exact_pct_to_frac(rep_pct_d)
            if rep_frac_d is None:
                rep_meets_exact = None
            else:
                rep_meets_exact = _zero_event_u95_compare(
                    a_obs_d,
                    n_raw,
                    f_d,
                    eta_d,
                    rep_frac_d,
                    u95_lower=u95_lower,
                    u95_upper=u95_upper,
                )
    has_bonferroni_alpha = (a_eta_d > Decimal(0)) or (a_f_d > Decimal(0))
    has_nonunit_factor = (eta_d < Decimal(1)) or (f_d < Decimal(1))
    _cal_src = calibration_sources if isinstance(calibration_sources, dict) else {}
    default_cal_assumptions: List[str] = []
    if eta_d == Decimal(1):
        default_cal_assumptions.append(
            "DEFAULT_UNIT_RECOVERY_ASSUMED(eta=1.0)"
            if _cal_src.get("eta_recovery_l95", "default") == "default"
            else "UNIT_RECOVERY_ASSUMED(eta=1.0)"
        )
    if f_d == Decimal(1):
        default_cal_assumptions.append(
            "DEFAULT_UNIT_ELIGIBILITY_ASSUMED(f_eligible=1.0)"
            if _cal_src.get("f_eligible", "default") == "default"
            else "UNIT_ELIGIBILITY_ASSUMED(f_eligible=1.0)"
        )
    elif _cal_src.get("f_eligible_provenance") == "default" and a_f_d == Decimal(0):
        default_cal_assumptions.append("DEFAULT_KNOWN_ELIGIBILITY_PROVENANCE_ASSUMED")
    # Directed 80-digit sound Bonferroni coverage bound (exact when resolved via Fraction, else directed ROUND_CEILING)
    with _clean_dec_ctx(80, ROUND_CEILING) as ctx_cov:
        if exact_joint_alpha_frac is not None:
            joint_alpha_upper = Decimal(exact_joint_alpha_frac.numerator) / Decimal(exact_joint_alpha_frac.denominator)
        else:
            joint_alpha_upper = a_obs_d + a_eta_d + a_f_d
        ctx_cov.rounding = ROUND_FLOOR
        cov_pct_lower = (Decimal(1) - joint_alpha_upper) * Decimal(100)
        cov_floor = cov_pct_lower.quantize(Decimal("0.01"), rounding=ROUND_FLOOR)
        ctx_cov.rounding = ROUND_CEILING
        ja_ceil_6 = joint_alpha_upper.quantize(Decimal("0.000001"), rounding=ROUND_CEILING)
        with _clean_dec_ctx(80, ROUND_HALF_EVEN):
            joint_alpha_str = format(joint_alpha_upper, ".6g") if ja_ceil_6 == joint_alpha_upper else f"{ja_ceil_6:.6f}"
            alpha_obs_fmt = format(a_obs_d, ".6g")
            eff_det_fmt = format(eff_lower, ".12g")
    if has_bonferroni_alpha:
        cov_regime = f"BONFERRONI_JOINT_COVERAGE_GE_{cov_floor:.2f}PCT_CALIBRATED"
        cov_cat = "CALIBRATED"
    elif has_nonunit_factor:
        cov_regime = (
            f"OBSERVATION_LEVEL_GE_{cov_floor:.2f}PCT_PARTIAL_FACTOR_WITH_DEFAULT_ASSUMPTIONS"
            if default_cal_assumptions
            else f"OBSERVATION_LEVEL_GE_{cov_floor:.2f}PCT_KNOWN_FACTORS_OBSERVATION_ALPHA_ONLY"
        )
        cov_cat = "PARTIAL_FACTOR_WITH_DEFAULTS" if default_cal_assumptions else "KNOWN_FACTORS_NO_BONFERRONI"
    else:
        cov_regime = (
            "OBSERVATION_LEVEL_95PCT_UNIT_RECOVERY_AND_ELIGIBILITY_ASSUMED"
            if a_obs_d == Decimal("0.05")
            else f"OBSERVATION_LEVEL_GE_{cov_floor:.2f}PCT_UNIT_RECOVERY_AND_ELIGIBILITY_ASSUMED"
        )
        cov_cat = "UNIT_ASSUMED"
    if cap_cmp is True:
        cp_status = "PASS"
    elif cap_cmp is False:
        cp_status = "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
    else:
        cp_status = "HOLD_ZERO_EVENT_CAP_COMPARISON_UNRESOLVED_WITHIN_BUDGET"
    res_cp: Dict[str, Any] = {
        "gate_id": gate_id,
        "status": cp_status,
        "empirical_clearance_status": (
            "MATHEMATICAL_COUNT_SUBGATE_ONLY_CONDITIONAL_ON_CALIBRATION_NO_EMPIRICAL_SAFETY_CLEARANCE"
            if passed
            else "HOLD_MATHEMATICAL_COUNT_SUBGATE_AND_EMPIRICAL_SAFETY_PENDING"
        ),
        "coverage_regime": cov_regime,
        "coverage_category": cov_cat,
        "default_calibration_assumptions": default_cal_assumptions,
        "sampling_independence_provenance": sampling_independence_provenance,
        "calibration_sources": calibration_sources,
        "k": 0,
        "n_raw": n_raw,
        "alpha_obs_str": alpha_obs_fmt,
        "alpha_joint_total_str": joint_alpha_str,
        "effective_detection_factor_str": eff_det_fmt,
        "u95_exact_sound_decimal_str": u95_str,
        "u95_exact_sound_decimal_pct_str": u95_pct_str,
        "u95_math_lower_pct_str": u95_math_lower_pct_str,
        "u95_math_upper_pct_exact_str": u95_math_upper_pct_exact_str,
        "cap_fraction_str": cap_str,
    }
    if rep_meets_exact is not None:
        res_cp["reported_u95_meets_exact_bound"] = rep_meets_exact
    return res_cp


_SP1C2C_VL9_COMPETENT_HLAC_ALLELES = frozenset({
    # Lin et al. 2023 (PMID:37264229, PMC10690437 Table 1): SP-1C (15 alleles) & SP-2C (3 alleles)
    "HLA-C*01:02", "HLA-C*03:02", "HLA-C*03:03", "HLA-C*03:04", "HLA-C*04:01",
    "HLA-C*04:03", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*08:01", "HLA-C*08:02",
    "HLA-C*12:02", "HLA-C*12:03", "HLA-C*14:02", "HLA-C*16:01", "HLA-C*16:02",
    "HLA-C*02:02", "HLA-C*15:02", "HLA-C*15:05",
})
_SP3C4C_NONFUNC_HLAC_ALLELES = frozenset({
    # Lin et al. 2023 (PMID:37264229, PMC10690437 Table 1): SP-3C (4 alleles) & SP-4C (1 allele)
    "HLA-C*07:01", "HLA-C*07:02", "HLA-C*07:04", "HLA-C*17:01", "HLA-C*18:01",
})
_PINNED_KNOWN_HLAC_TWO_FIELD_INVENTORY = _SP1C2C_VL9_COMPETENT_HLAC_ALLELES | _SP3C4C_NONFUNC_HLAC_ALLELES


CANONICAL_INDEPENDENT_SAMPLING_DESIGNS = frozenset({
    "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES",
    "INDEPENDENT",
    "IID",
    "INDEPENDENT_TRIALS",
    "INDEPENDENT_REPLICATES",
    "INDEPENDENT_BIOLOGICAL_REPLICATES",
})


def _reconcile_sampling_independence_flags(
    primary_map: Dict[str, Any],
    *,
    prefix: str | None = None,
    panel_map: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Fail-closed reconciliation of sampling-independence keys across assay and panel scopes.

    Separates clustering assertions (`clustered_sampling`/`is_clustered`) from independence assertions
    (`independent_trials`/`independent_sample_units`/`sampling_independence_model`/`sampling_design`),
    rejects non-ASCII spellings BEFORE normalization, and reports malformed values distinctly.
    """
    scopes: List[tuple[Any, str | None]] = [(primary_map, prefix)]
    if panel_map is not None:
        scopes.append((panel_map, None))
    elif prefix is not None:
        scopes.append((primary_map, None))

    clustered_true = False
    clustered_false = False
    indep_true = False
    indep_false = False
    malformed_flag = False

    for src, pfx in scopes:
        if not isinstance(src, dict):
            continue
        for key in ("clustered_sampling", "is_clustered", "independent_trials", "independent_sample_units", "sampling_independence_model", "sampling_design"):
            full_k = f"{pfx}_{key}" if pfx else key
            if full_k not in src:
                continue
            val = src[full_k]
            if key in ("clustered_sampling", "is_clustered"):
                if val is True:
                    clustered_true = True
                elif val is False:
                    clustered_false = True
                else:
                    malformed_flag = True
            elif key in ("independent_trials", "independent_sample_units"):
                if val is False:
                    indep_false = True
                elif val is True:
                    indep_true = True
                else:
                    malformed_flag = True
            elif key == "sampling_independence_model":
                if isinstance(val, str) and val.isascii() and val == "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES":
                    indep_true = True
                elif isinstance(val, str) and val.isascii() and val == val.strip() and val:
                    indep_false = True
                else:
                    malformed_flag = True
            elif key == "sampling_design":
                if isinstance(val, str) and val.isascii() and val == val.strip() and val:
                    norm_sd = re.sub(r"[ -]+", "_", val, flags=re.ASCII).upper()
                    if norm_sd in CANONICAL_INDEPENDENT_SAMPLING_DESIGNS:
                        indep_true = True
                    else:
                        indep_false = True
                else:
                    malformed_flag = True

    conflict = bool(
        (clustered_true and clustered_false)
        or (indep_true and indep_false)
        or (clustered_true and indep_true)
    )
    non_independent = bool(clustered_true or indep_false)
    if conflict:
        prov = "CONFLICTING_INDEPENDENT_AND_NON_INDEPENDENT_DECLARATIONS"
    elif malformed_flag:
        prov = "MALFORMED_OR_NON_ASCII_SAMPLING_DECLARATION"
    elif clustered_true:
        prov = "EXPLICIT_CLUSTERED_SAMPLING"
    elif indep_false and clustered_false:
        prov = "EXPLICIT_UNCLUSTERED_BUT_DEPENDENT_SAMPLING"
    elif indep_false:
        prov = "EXPLICIT_NON_INDEPENDENT_OR_UNRECOGNIZED"
    elif indep_true:
        prov = "EXPLICIT_INDEPENDENT"
    elif clustered_false:
        prov = "EXPLICIT_NON_CLUSTERED_ONLY_INDEPENDENCE_ASSUMED"
    else:
        prov = "UNSPECIFIED_DEFAULT_ASSUMED_INDEPENDENT"
    return {
        "sampling_independence_model": (
            "NON_INDEPENDENT_OR_CLUSTERED"
            if (non_independent or malformed_flag)
            else "INDEPENDENT_BIOLOGICAL_OPPORTUNITIES"
        ),
        "clustered_sampling": bool(non_independent or malformed_flag),
        "sampling_independence_conflict": conflict,
        "sampling_malformed": bool(malformed_flag and not conflict),
        "sampling_independence_provenance": prov,
    }


_CAL_PARAM_DEFAULTS: Dict[str, Any] = {
    "eta_recovery_l95": Decimal(1),
    "f_eligible": Decimal(1),
    "f_eligible_provenance": "KNOWN_OR_CONSERVATIVELY_ESTABLISHED",
    "alpha_obs": None,
    "alpha_eta_recovery": None,
    "alpha_f_eligible": None,
}


def _resolve_assay_calibration(
    local_map: Dict[str, Any],
    *,
    local_prefix: str | None = None,
    panel_map: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Shared typed calibration resolver across D1-04, D1-08, L2/L3 components, and Gate-0V Day-28.

    Uses explicit key presence distinct from JSON null:
      - Omitted local key inherits the panel value (or built-in default if panel omitted).
      - Explicit JSON null (`key in src and src[key] is None`) is rejected as malformed/invalid override (HOLD).
      - Conflicting values across panel and local scopes for the SAME parameter are rejected (HOLD).
      - Bound parameter groups ({eta_recovery_l95, alpha_eta_recovery} and {f_eligible, f_eligible_provenance,
        alpha_f_eligible}) cannot be split partially across local and panel scopes (HOLD).
    """
    p_src = panel_map if isinstance(panel_map, dict) else (local_map if local_prefix is not None else None)
    resolved: Dict[str, Any] = {}
    sources: Dict[str, str] = {}
    local_present_keys: set[str] = set()
    panel_present_keys: set[str] = set()
    hold_reason: str | None = None

    for param, default_val in _CAL_PARAM_DEFAULTS.items():
        lk = f"{local_prefix}_{param}" if local_prefix else param
        has_local = isinstance(local_map, dict) and (lk in local_map)
        has_panel = isinstance(p_src, dict) and (param in p_src) and (local_prefix is not None or panel_map is not None)
        if has_local:
            local_present_keys.add(param)
        if has_panel:
            panel_present_keys.add(param)

        l_val = local_map.get(lk) if has_local else _MISSING_ARG
        p_val = p_src.get(param) if has_panel else _MISSING_ARG

        if has_local and l_val is None:
            hold_reason = hold_reason or "HOLD_EXPLICIT_NULL_CALIBRATION_PARAMETER"
            resolved[param] = default_val
            sources[param] = "local_explicit_null_rejected"
            continue
        if has_panel and p_val is None:
            hold_reason = hold_reason or "HOLD_EXPLICIT_NULL_CALIBRATION_PARAMETER"
            resolved[param] = default_val
            sources[param] = "panel_explicit_null_rejected"
            continue

        if has_local and has_panel:
            # Check strict value equivalence (using Decimal for numeric values, exact str for provenance)
            if param == "f_eligible_provenance":
                same = isinstance(l_val, str) and isinstance(p_val, str) and l_val == p_val
            else:
                dl = _safe_exact_decimal(l_val)
                dp = _safe_exact_decimal(p_val)
                same = (dl is not None and dp is not None and dl == dp)
            if not same:
                hold_reason = hold_reason or "HOLD_CONFLICTING_PANEL_AND_LOCAL_CALIBRATION"
            resolved[param] = l_val
            sources[param] = "local_and_panel_matched" if same else "local_panel_conflict"
        elif has_local:
            resolved[param] = l_val
            sources[param] = "local"
        elif has_panel:
            resolved[param] = p_val
            sources[param] = "panel"
        else:
            resolved[param] = default_val
            sources[param] = "default"

    # Enforce bound parameter groups: never combine a local factor with a panel-scoped sibling alpha/provenance
    for grp in (
        ("eta_recovery_l95", "alpha_eta_recovery"),
        ("f_eligible", "f_eligible_provenance", "alpha_f_eligible"),
    ):
        g_loc = {k for k in grp if k in local_present_keys}
        g_pan = {k for k in grp if k in panel_present_keys}
        if g_loc and g_pan and (g_loc != g_pan):
            hold_reason = hold_reason or "HOLD_PARTIAL_LOCAL_CALIBRATION_SPLITS_BOUND_GROUP"

    # Fail closed if caller allocates Bonferroni alpha or declares an estimated bound while omitting the estimate itself!
    _a_eta_res = _safe_exact_decimal(resolved.get("alpha_eta_recovery")) if resolved.get("alpha_eta_recovery") is not None else None
    if sources.get("alpha_eta_recovery") != "default" and sources.get("eta_recovery_l95") == "default":
        if _a_eta_res is None or _a_eta_res > Decimal(0):
            hold_reason = hold_reason or "HOLD_MISSING_ETA_RECOVERY_ESTIMATE_FOR_ALPHA_ALLOCATION"
    if (
        resolved.get("f_eligible_provenance") == "INDEPENDENTLY_ESTIMATED_WITH_ALPHA_ALLOCATION"
        and sources.get("f_eligible") == "default"
    ):
        hold_reason = hold_reason or "HOLD_MISSING_F_ELIGIBLE_ESTIMATE_FOR_ALPHA_ALLOCATION"

    samp_flags = _reconcile_sampling_independence_flags(local_map, prefix=local_prefix, panel_map=panel_map)
    return {
        **resolved,
        **samp_flags,
        "calibration_hold_reason": hold_reason,
        "calibration_sources": sources,
    }


def _has_positive_gate0v_event(record: Any) -> bool:
    if not isinstance(record, dict):
        return False
    d28_k = record.get("day28_vp1_ie1_positive_k")
    return isinstance(d28_k, int) and not isinstance(d28_k, bool) and d28_k >= 1


def _inspect_adverse_gate0v_failure(
    record: Any,
    *,
    day14_min: Decimal = Decimal("3.0"),
    day14_l95_min: Decimal = Decimal("3.0"),
    early_72h_max: Decimal = Decimal("2.0"),
) -> str | None:
    """Detect explicit verdict=='FAIL' or numeric viral-clearance threshold failures on a Gate-0V record
    so adverse Gate-0V failures are preserved symmetrically regardless of whether a lot conflict occurs
    at the panel observation level or inside the nested Gate-0V record."""
    if not isinstance(record, dict):
        return None
    if record.get("gate0v_verdict") == "FAIL":
        return "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN"
    d14 = record.get("cvb_day14_log10_reduction")
    if d14 is not None and not (isinstance(d14, str) and d14 == "UNAVAILABLE"):
        d14_d = _safe_exact_decimal(d14)
        if d14_d is not None and Decimal("0") <= d14_d < day14_min:
            return "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN"
    d14_l95 = record.get("cvb_day14_log10_reduction_l95")
    if d14_l95 is not None and not (isinstance(d14_l95, str) and d14_l95 == "UNAVAILABLE"):
        d14_l95_d = _safe_exact_decimal(d14_l95)
        if d14_l95_d is not None and Decimal("0") <= d14_l95_d < day14_l95_min:
            return "FAIL_D1_L1B_GATE0V_DAY14_L95_BELOW_3LOG_SWITCH_TO_L1A_IL1RN"
    t72 = record.get("cvb4_72h_early_titer_log10_pfu_u95")
    if t72 is not None and not (isinstance(t72, str) and t72 == "UNAVAILABLE"):
        t72_d = _safe_exact_decimal(t72)
        if t72_d is not None and early_72h_max < t72_d <= Decimal("12.0"):
            return "FAIL_D1_L1B_GATE0V_EARLY_72H_TITER_VETO_SWITCH_TO_L1A_IL1RN"
    return None


def _evaluate_structured_lot_matched_gate0v_record(
    record: Any,
    *,
    evaluated_lot_id: Any,
    day14_min: Decimal = Decimal("3.0"),
    day14_l95_min: Decimal = Decimal("3.0"),
    early_72h_max: Decimal = Decimal("2.0"),
    day28_u95_max: Decimal = Decimal("0.0005"),
    day28_n_min: int = 5990,
    panel_observations: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Validate structured lot-matched global_gate0v_record (NEVER accepts a bare caller boolean).

    Preserves unconditional FAIL_HARD_STOP precedence for ANY recorded positive Day-28 VP1/IE1 event (k >= 1)
    even when lot metadata is missing, malformed, or mismatched (recording both the hard stop and lot hold!).
    """
    pos_event = _has_positive_gate0v_event(record)
    lot_hold: str | None = None
    if not _is_canonical_identifier(evaluated_lot_id):
        lot_hold = "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
    elif record is None:
        return {
            "cleared": False,
            "failed": False,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "HOLD_D1_L1B_STRUCTURED_LOT_MATCHED_GLOBAL_GATE0V_RECORD_REQUIRED",
        }
    elif isinstance(record, bool) or not isinstance(record, dict):
        return {
            "cleared": False,
            "failed": False,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "HOLD_D1_L1B_BARE_BOOLEAN_OR_NON_DICT_GATE0V_REJECTED",
        }
    else:
        mcb_lot_id = record.get("mcb_lot_id")
        if not _is_canonical_identifier(mcb_lot_id):
            lot_hold = "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID"
        else:
            _nested_lot_vals: List[str] = []
            _nested_lot_bad = False
            for _nlk in _LOT_ALIAS_KEYS:
                if _nlk not in record:
                    continue
                _nlv = record[_nlk]
                if _nlv is None:
                    _nested_lot_bad = True
                elif not _is_canonical_identifier(_nlv):
                    _nested_lot_bad = True
                elif _nlv not in _nested_lot_vals:
                    _nested_lot_vals.append(_nlv)
            if _nested_lot_bad:
                lot_hold = "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID"
            elif len(_nested_lot_vals) > 1:
                lot_hold = "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
            elif mcb_lot_id != evaluated_lot_id:
                lot_hold = "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"

    if lot_hold is not None:
        if pos_event:
            return {
                "cleared": False,
                "failed": True,
                "hard_stop": True,
                "lot_binding_invalid": True,
                "reason": "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
                "lot_hold_reason": lot_hold,
                "adverse_event_pending_lot_attribution": True,
            }
        # Preserve explicit verdict=='FAIL' or numeric viral-clearance failure alongside nested lot-attribution hold!
        _adv_fail_reason = _inspect_adverse_gate0v_failure(
            record, day14_min=day14_min, day14_l95_min=day14_l95_min, early_72h_max=early_72h_max,
        )
        return {
            "cleared": False,
            "failed": bool(_adv_fail_reason is not None),
            "hard_stop": False,
            "lot_binding_invalid": True,
            "reason": _adv_fail_reason if _adv_fail_reason is not None else lot_hold,
            "lot_hold_reason": lot_hold,
        }
    mcb_lot_id = record.get("mcb_lot_id")

    d14 = record.get("cvb_day14_log10_reduction")
    d14_l95 = record.get("cvb_day14_log10_reduction_l95")
    t72_u95 = record.get("cvb4_72h_early_titer_log10_pfu_u95")
    d28_k = record.get("day28_vp1_ie1_positive_k")
    d28_n = record.get("day28_vp1_ie1_cells_scored_n")
    verdict = record.get("gate0v_verdict")

    # 1. Unconditional FAIL_HARD_STOP precedence for ANY recorded positive Day-28 VP1/IE1 event (k >= 1),
    #    even if Day-14 reduction, L95, 72h early titer, n, or gate0v_verdict is missing/PENDING.
    if isinstance(d28_k, int) and not isinstance(d28_k, bool) and d28_k >= 1:
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": True,
            "lot_binding_invalid": False,
            "reason": "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
        }
    if verdict == "FAIL":
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN",
        }
    d14_present = not (d14 is None or (isinstance(d14, str) and d14 == "UNAVAILABLE"))
    d14_l95_present = not (d14_l95 is None or (isinstance(d14_l95, str) and d14_l95 == "UNAVAILABLE"))
    t72_present = not (t72_u95 is None or (isinstance(t72_u95, str) and t72_u95 == "UNAVAILABLE"))
    d28_k_present = not (d28_k is None or (isinstance(d28_k, str) and d28_k == "UNAVAILABLE"))
    d28_n_present = not (d28_n is None or (isinstance(d28_n, str) and d28_n == "UNAVAILABLE"))
    d14_d = _safe_exact_decimal(d14) if d14_present else None
    d14_l95_d = _safe_exact_decimal(d14_l95) if d14_l95_present else None
    t72_u95_d = _safe_exact_decimal(t72_u95) if t72_present else None
    # Reject invalid / out-of-range numeric fields BEFORE returning incomplete HOLD for missing sibling fields!
    if (
        (d14_present and (d14_d is None or not (Decimal("0") <= d14_d <= Decimal("12.0"))))
        or (d14_l95_present and (d14_l95_d is None or not (Decimal("0") <= d14_l95_d <= Decimal("12.0"))))
        or (d14_d is not None and d14_l95_d is not None and d14_l95_d > d14_d)
        or (t72_present and (t72_u95_d is None or not (Decimal("0") <= t72_u95_d <= Decimal("12.0"))))
        or (d28_k_present and (isinstance(d28_k, bool) or not isinstance(d28_k, int) or d28_k < 0))
        or (d28_n_present and (isinstance(d28_n, bool) or not isinstance(d28_n, int) or d28_n < 1))
        or (d28_k_present and d28_n_present and d28_k > d28_n)
    ):
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE",
        }
    if d14_d is not None and Decimal("0") <= d14_d < day14_min:
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN",
        }
    if d14_l95_d is not None and Decimal("0") <= d14_l95_d < day14_l95_min:
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_L1B_GATE0V_DAY14_L95_BELOW_3LOG_SWITCH_TO_L1A_IL1RN",
        }
    if t72_u95_d is not None and early_72h_max < t72_u95_d <= Decimal("12.0"):
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_L1B_GATE0V_EARLY_72H_TITER_VETO_SWITCH_TO_L1A_IL1RN",
        }

    if any(
        v is None or (isinstance(v, str) and v == "UNAVAILABLE")
        for v in (d14, d14_l95, t72_u95, d28_k, d28_n, verdict)
    ):
        return {
            "cleared": False,
            "failed": False,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "HOLD_D1_L1B_GATE0V_INCOMPLETE_RECORD",
        }
    if (
        d14_d is None or not (Decimal("0") <= d14_d <= Decimal("12.0"))
        or d14_l95_d is None or not (Decimal("0") <= d14_l95_d <= Decimal("12.0"))
        or d14_l95_d > d14_d
        or t72_u95_d is None or not (Decimal("0") <= t72_u95_d <= Decimal("12.0"))
    ):
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE",
        }
    if not (isinstance(d28_k, int) and not isinstance(d28_k, bool)) or not (isinstance(d28_n, int) and not isinstance(d28_n, bool)):
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE",
        }
    if d28_k < 0 or d28_n < 1 or d28_k > d28_n:
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE",
        }
    if verdict != "PASS":
        return {
            "cleared": False,
            "failed": False,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "HOLD_D1_L1B_GATE0V_VERDICT_NOT_PASS",
        }

    _g0v_cal_kw = _resolve_assay_calibration(record, local_prefix=None, panel_map=panel_observations)
    cp_v = evaluate_sound_decimal_clopper_pearson_zero_event(
        gate_id="GATE0V_DAY28_VP1_IE1",
        k=d28_k,
        n_raw=d28_n,
        max_u95_fraction=day28_u95_max,
        **_g0v_cal_kw,
    )
    if cp_v["status"] == "FAIL_HARD_STOP":
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": True,
            "lot_binding_invalid": False,
            "reason": "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
            "day28_u95_str": cp_v.get("u95_exact_sound_decimal_str"),
            "cp_result": cp_v,
        }
    if d14_d < day14_min or d14_l95_d < day14_l95_min:
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_L1B_GATE0V_CD46WT_VIRAL_CLEARANCE_SWITCH_TO_L1A_IL1RN",
            "day28_u95_str": cp_v.get("u95_exact_sound_decimal_str"),
            "cp_result": cp_v,
        }
    if t72_u95_d > early_72h_max:
        return {
            "cleared": False,
            "failed": True,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": "FAIL_D1_L1B_GATE0V_EARLY_72H_TITER_VETO_SWITCH_TO_L1A_IL1RN",
            "day28_u95_str": cp_v.get("u95_exact_sound_decimal_str"),
            "cp_result": cp_v,
        }
    if cp_v["status"] != "PASS" or d28_n < day28_n_min:
        g0v_hold_reason = (
            "HOLD_D1_L1B_GATE0V_DAY28_INSUFFICIENT_DENOMINATOR"
            if cp_v["status"] in ("PASS", "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR")
            else f"HOLD_D1_L1B_GATE0V_DAY28_{cp_v['status']}"
        )
        return {
            "cleared": False,
            "failed": False,
            "hard_stop": False,
            "lot_binding_invalid": False,
            "reason": g0v_hold_reason,
            "day28_u95_str": cp_v.get("u95_exact_sound_decimal_str"),
            "cp_result": cp_v,
        }
    return {
        "cleared": True,
        "failed": False,
        "hard_stop": False,
        "lot_binding_invalid": False,
        "reason": "STANDALONE_SUBCHECK_PASSED_LIVE_LIMS_LEDGER_UNVERIFIED",
        "mcb_lot_id": mcb_lot_id,
        "day28_u95_str": cp_v.get("u95_exact_sound_decimal_str"),
        "cp_result": cp_v,
    }


def _extract_dec_from_cell_threshold(threshold_str: str, pattern: str, label: str) -> Decimal:
    m = re.search(pattern, threshold_str)
    if not m:
        raise ValueError(f"Controlling 45-cell contract missing required threshold clause {label!r} in {threshold_str!r}")
    dv = _safe_exact_decimal(m.group(1))
    if dv is None:
        raise ValueError(f"Controlling 45-cell contract threshold {label!r} is not a valid Decimal: {m.group(1)!r}")
    return dv


EXPECTED_PLAYBOOK_LAYER_SCOPE_BODY_SHA256 = "da7f53e8102cbedbb23e24afcee1b188bb175c46a98a0d1aca67311c98475afa"


def _strict_typed_equal(a: Any, b: Any) -> bool:
    """Strict type-aware JSON equality where bool != int AND int != float/Decimal (while float/Decimal compare exact)."""
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, int) or isinstance(b, int):
        return type(a) is int and type(b) is int and a == b
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_strict_typed_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_strict_typed_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, (float, Decimal)) and isinstance(b, (float, Decimal)):
        da = _safe_exact_decimal(a)
        db = _safe_exact_decimal(b)
        return da is not None and db is not None and da == db
    return type(a) is type(b) and a == b


def _verify_playbook_sec4120(
    repo_root: Path,
    *,
    require_canonical_sha: bool = True,
) -> Dict[str, Any]:
    """Parse docs/STAGE-0-EXECUTION-PLAYBOOK.md §4.1.20 via CommonMark heading scanner & verify canonical SHA-256."""
    pb_path = repo_root / PLAYBOOK_REL_PATH
    if not pb_path.is_file():
        raise ValueError(f"Missing Stage-0 execution playbook: {pb_path}")
    pb_norm = pb_path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
    all_headings = _parse_commonmark_headings(pb_norm)
    matched_4120 = [
        h for h in all_headings
        if re.search(r"(?<![0-9A-Za-z.])4\.1\.20(?![0-9A-Za-z.])", h["rendered_title"])
    ]
    pb_expected_prefix = "#### 4.1.20 Preregistered Autologous `D1b-Auto-C v2`"
    if len(matched_4120) != 1:
        raise ValueError(
            f"CommonMark Playbook audit requires exactly 1 rendered 4.1.20 heading, found {len(matched_4120)}: "
            f"{[(h['kind'], h['level'], h['container_depth'], h['start_line'] + 1) for h in matched_4120]}"
        )
    h4120 = matched_4120[0]
    if (
        h4120["kind"] != "atx"
        or h4120["level"] != 4
        or h4120["container_depth"] != 0
        or not h4120["raw_first_line"].startswith(pb_expected_prefix)
    ):
        raise ValueError(f"Rendered 4.1.20 heading has unexpected CommonMark structure: {h4120}")

    start_off = h4120["start_offset"]
    next_heading_off = len(pb_norm)
    for h in all_headings:
        if h["start_offset"] > start_off and h["container_depth"] == 0 and h["level"] <= 4:
            next_heading_off = h["start_offset"]
            break
    end_marker = "<!-- END_IIP0035_PLAYBOOK_SECTION_4_1_20 -->"
    end_idx = pb_norm.find(end_marker, start_off, next_heading_off)
    if end_idx == -1:
        raise ValueError(f"Missing required post-4.1.20 end marker {end_marker!r} before next top-level heading")
    # Require EXACTLY ONE literal standalone OPENISLET_CURRENT_LAYER_SCOPE footer pair in canonical mode
    # (rejecting complete deletion, duplication, incomplete markers, or any residual prose/comments before next heading).
    trailing_under_4120 = pb_norm[end_idx + len(end_marker):next_heading_off]
    _ls_matches = list(
        re.finditer(
            r"(?m)^<!-- OPENISLET_CURRENT_LAYER_SCOPE -->\n(.*?)\n<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->[ \t]*$",
            trailing_under_4120,
            re.S,
        )
    )
    footer_present = len(_ls_matches) == 1
    footer_verified = False
    if len(_ls_matches) > 1:
        raise ValueError(f"Duplicate OPENISLET_CURRENT_LAYER_SCOPE footers ({len(_ls_matches)}) under Playbook §4.1.20")
    if len(_ls_matches) == 1:
        _ls_m = _ls_matches[0]
        _ls_body = _ls_m.group(1)
        if re.search(r"(?m)^[ ]{0,3}#{1,6}\s", _ls_body):
            raise ValueError("Unexpected Markdown heading smuggled inside OPENISLET_CURRENT_LAYER_SCOPE footer under Playbook §4.1.20")
        if re.search(r"(?i)\b(?:D1-0[1-9]|D1b|DEC-114|IIP-0035|KIR3DL1|VTCN1|HLA-C\*03:04|AP1903|Gate\s*0V|C5b-9|overrid\w*|waiv\w*|exempt\w*|bypass\w*)\b", _ls_body):
            raise ValueError("Unexpected D1b-Auto-C / DEC-114 or override prose inside OPENISLET_CURRENT_LAYER_SCOPE footer under Playbook §4.1.20")
        _ls_sha = hashlib.sha256(_ls_body.strip().encode("utf-8")).hexdigest()
        footer_verified = (_ls_sha == EXPECTED_PLAYBOOK_LAYER_SCOPE_BODY_SHA256)
        if require_canonical_sha and not footer_verified:
            raise ValueError(
                f"OPENISLET_CURRENT_LAYER_SCOPE footer under Playbook §4.1.20 SHA-256 mismatch: found {_ls_sha}, expected {EXPECTED_PLAYBOOK_LAYER_SCOPE_BODY_SHA256}"
            )
        trailing_residual = trailing_under_4120[:_ls_m.start()] + trailing_under_4120[_ls_m.end():]
    else:
        if require_canonical_sha:
            raise ValueError(
                "Missing required canonical OPENISLET_CURRENT_LAYER_SCOPE generated footer under Playbook §4.1.20"
            )
        trailing_residual = trailing_under_4120
    if trailing_residual.strip():
        raise ValueError(
            "Unexpected content (visible text or unpinned comment) appended under Playbook §4.1.20 after end marker before next rendered heading: "
            f"{trailing_residual.strip()[:120]!r}"
        )
    pb_sec_full = pb_norm[start_off:end_idx].strip()
    pb_sha = hashlib.sha256(pb_sec_full.encode("utf-8")).hexdigest()
    if require_canonical_sha and pb_sha != EXPECTED_DEC114_PLAYBOOK_SEC4120_SHA256:
        raise ValueError(
            f"DEC-114 Playbook section 4.1.20 SHA-256 mismatch: found {pb_sha}, expected {EXPECTED_DEC114_PLAYBOOK_SEC4120_SHA256}"
        )
    return {
        "pb_sec_full": pb_sec_full,
        "pb_sha": pb_sha,
        "footer_present": footer_present,
        "footer_verified": footer_verified,
        "canonical_sha_verified": bool(pb_sha == EXPECTED_DEC114_PLAYBOOK_SEC4120_SHA256 and footer_verified),
        "all_headings": all_headings,
        "matched_4120": matched_4120,
        "h4120": h4120,
    }


def _canonicalize_json_for_digest(obj: Any) -> Any:
    """Convert exact finite Decimal scalars matching canonical float JSON representation so Decimal-parsed product_defs hash identically."""
    if isinstance(obj, dict):
        return {k: _canonicalize_json_for_digest(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_canonicalize_json_for_digest(v) for v in obj]
    if isinstance(obj, Decimal):
        if not obj.is_finite():
            raise ValueError(f"Non-finite Decimal in product_def: {obj!r}")
        fv = float(obj)
        if not math.isfinite(fv) or Decimal(str(fv)) != obj:
            raise ValueError(f"Non-canonical Decimal value in product_def: {obj!r}")
        return fv
    return obj


def _load_controlling_contract_for_product(
    pvid: str,
    *,
    decision_tree_block: Dict[str, Any] | None = None,
    repo_root: Path = ROOT,
    require_canonical_sha: bool = True,
) -> Dict[str, Any]:
    """Load, verify (45-cell DEC-114 table, Playbook §4.1.20, AND staged product registry/FASTAs), and parse clauses."""
    pb_info = _verify_playbook_sec4120(repo_root, require_canonical_sha=require_canonical_sha)
    staged_registry_verified = False
    if require_canonical_sha:
        # Enforce identical staged-FASTA & product-registry cryptographic verification across BOTH direct Python and CLI!
        pc.load_product_contracts(repo_root, require_staged_fastas=True)
        staged_registry_verified = True
    disk_blk: Dict[str, Any] | None = None
    disk_sha: str | None = None
    dt_path = repo_root / DECISION_TREE_REL_PATH
    if require_canonical_sha or decision_tree_block is None:
        if not dt_path.is_file():
            raise ValueError(f"Missing controlling Stage-0 decision tree: {dt_path}")
        dt_doc = _load_json_strict(dt_path)
        disk_blk = dt_doc.get("autologous_d1b_auto_c_v2_vanguard_gates_iip_0035")
        if not isinstance(disk_blk, dict):
            raise ValueError("Missing DEC-114 block autologous_d1b_auto_c_v2_vanguard_gates_iip_0035 in stage0_decision_tree_v1.json")
        disk_sha = hashlib.sha256(
            json.dumps(disk_blk, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        if require_canonical_sha and disk_sha != EXPECTED_DEC114_DECISION_TREE_SHA256:
            raise ValueError(
                f"On-disk controlling DEC-114 decision-tree SHA-256 mismatch: found {disk_sha}, expected {EXPECTED_DEC114_DECISION_TREE_SHA256}"
            )
    elif dt_path.is_file():
        # Best-effort non-raising on-disk check when require_canonical_sha=False so status fields never
        # misreport an unverified/absent installed decision tree as canonical-verified!
        try:
            _dt_doc_opt = _load_json_strict(dt_path)
            _db_opt = _dt_doc_opt.get("autologous_d1b_auto_c_v2_vanguard_gates_iip_0035")
            if isinstance(_db_opt, dict):
                disk_blk = _db_opt
                disk_sha = hashlib.sha256(
                    json.dumps(disk_blk, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
                ).hexdigest()
        except Exception:
            disk_blk = None
            disk_sha = None
    if decision_tree_block is None:
        if disk_blk is None:
            raise ValueError("Missing on-disk DEC-114 decision-tree block")
        blk = disk_blk
    else:
        blk = decision_tree_block
    dt_sha = hashlib.sha256(
        json.dumps(blk, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    on_disk_dt_verified = bool(
        disk_blk is not None
        and disk_sha == EXPECTED_DEC114_DECISION_TREE_SHA256
        and _strict_typed_equal(blk, disk_blk)
    )
    if require_canonical_sha and (dt_sha != EXPECTED_DEC114_DECISION_TREE_SHA256 or not on_disk_dt_verified):
        raise ValueError(
            f"Controlling DEC-114 decision-tree SHA-256 mismatch: found {dt_sha}, expected {EXPECTED_DEC114_DECISION_TREE_SHA256}"
        )
    pbg = blk.get("product_by_gate_contract")
    if not isinstance(pbg, dict) or pbg.get("cell_count") != 45 or not isinstance(pbg.get("cells"), list) or len(pbg["cells"]) != 45:
        raise ValueError("Controlling DEC-114 product_by_gate_contract must define the exact 45-cell Cartesian contract")
    all_cell_map: Dict[tuple[str, str], Dict[str, Any]] = {}
    for c in pbg["cells"]:
        if not isinstance(c, dict) or set(c.keys()) != REQUIRED_CELL_KEYS:
            raise ValueError(f"Malformed typed cell in controlling product_by_gate_contract: {c!r}")
        if any(not isinstance(v, str) or not v.strip() or _PLACEHOLDER_TOKEN_RE.search(v) for v in c.values()):
            raise ValueError(f"Empty, non-string, or placeholder field in typed cell {c.get('product_variant_id')}:{c.get('gate_id')}")
        if c["operator"] not in VALID_CELL_OPERATORS:
            raise ValueError(f"Unrecognized typed cell operator {c['operator']!r} on {c['product_variant_id']}:{c['gate_id']}")
        if c["disposition"] not in VALID_CELL_DISPOSITIONS:
            raise ValueError(f"Unrecognized typed cell disposition {c['disposition']!r} on {c['product_variant_id']}:{c['gate_id']}")
        mdd_c = c["missing_data_disposition"]
        if not mdd_c.startswith("HOLD_") or re.search(r"\bPASS\b", mdd_c, re.IGNORECASE):
            raise ValueError(f"Typed cell missing_data_disposition must be fail-closed HOLD_* on {c['product_variant_id']}:{c['gate_id']}: {mdd_c!r}")
        ckey = (c["product_variant_id"], c["gate_id"])
        if ckey in all_cell_map:
            raise ValueError(f"Duplicate product-by-gate cell in controlling contract: {ckey}")
        all_cell_map[ckey] = c
    expected_cart = {(p, g) for p in pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS for g in EXPECTED_45_GATES}
    if set(all_cell_map.keys()) != expected_cart:
        raise ValueError("Controlling product_by_gate_contract does not match the exact 5 products x 9 gates Cartesian set")

    pvid_cells: Dict[str, Dict[str, Any]] = {
        gid: all_cell_map[(pvid, gid)]
        for gid in EXPECTED_45_GATES
        if (pvid, gid) in all_cell_map
    }
    # Verify expected operators & dispositions across ALL 45 cells (all 5 products x 9 gates) even in diagnostic mode!
    for check_pvid in pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS:
        if check_pvid in D1B_FULLSTACK_CANDIDATE_IDS:
            exp_ops = EXPECTED_L1B_GATE_OPERATORS if check_pvid == "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR" else EXPECTED_L1A_GATE_OPERATORS
            for gid in EXPECTED_45_GATES:
                cell = all_cell_map[(check_pvid, gid)]
                exp_disp = (
                    "REQUIRED_EXACT_PRODUCT_WITH_LOT_MATCHED_GATE0V"
                    if (check_pvid == "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR" and gid == "D1-09")
                    else "REQUIRED_EXACT_PRODUCT"
                )
                if cell["operator"] != exp_ops[gid] or cell["disposition"] != exp_disp:
                    raise ValueError(
                        f"Controlling contract operator/disposition mismatch on {check_pvid}:{gid}: "
                        f"got ({cell['operator']!r}, {cell['disposition']!r}), expected ({exp_ops[gid]!r}, {exp_disp!r})"
                    )
        elif check_pvid in D1B_COMPONENT_DONOR_IDS:
            comp_app = COMPONENT_APPLICABLE_GATES[check_pvid]
            comp_exp_ops = EXPECTED_COMPONENT_GATE_OPERATORS[check_pvid]
            for gid in EXPECTED_45_GATES:
                cell = all_cell_map[(check_pvid, gid)]
                if gid in comp_app:
                    if cell["operator"] != comp_exp_ops[gid] or cell["disposition"] != "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK":
                        raise ValueError(
                            f"Component cell {check_pvid}:{gid} operator/disposition mismatch: "
                            f"got ({cell['operator']!r}, {cell['disposition']!r}), "
                            f"expected ({comp_exp_ops[gid]!r}, 'COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK')"
                        )
                else:
                    if cell["operator"] != "NOT_APPLICABLE" or cell["disposition"] != "NOT_APPLICABLE_EXCLUDED_FROM_FULLSTACK":
                        raise ValueError(f"Inapplicable component cell {check_pvid}:{gid} must be NOT_APPLICABLE / NOT_APPLICABLE_EXCLUDED_FROM_FULLSTACK")

    cell_missing_dispositions: Dict[str, str] = {
        gid: str(pvid_cells[gid]["missing_data_disposition"])
        for gid in pvid_cells
    }

    is_component_excluded = bool(
        len(pvid_cells) == 9
        and all("EXCLUDED_FROM_FULLSTACK" in str(c.get("disposition") or "") for c in pvid_cells.values())
    )
    thresholds: Dict[str, Any] = {}
    if len(pvid_cells) == 9 and not is_component_excluded:
        t02 = str(pvid_cells["D1-02"]["threshold"])
        t03 = str(pvid_cells["D1-03"]["threshold"])
        t04 = str(pvid_cells["D1-04"]["threshold"])
        t05 = str(pvid_cells["D1-05"]["threshold"])
        t06 = str(pvid_cells["D1-06"]["threshold"])
        t07 = str(pvid_cells["D1-07"]["threshold"])
        t08 = str(pvid_cells["D1-08"]["threshold"])
        t09 = str(pvid_cells["D1-09"]["threshold"])
        thresholds = {
            "min_mm_reject_le": int(_extract_dec_from_cell_threshold(t02, r"reject_le_([0-9]+)mm", "D1-02:min_mm_reject_le")),
            "min_mm_max": int(_extract_dec_from_cell_threshold(t02, r"min_mm_to_HLAC_E_valid_range_0_to_([0-9]+)", "D1-02:min_mm_max")),
            "wt_hlac_min": _extract_dec_from_cell_threshold(t02, r"wt_hlac_retention_pct>=([0-9.eE+-]+)", "D1-02:wt_hlac"),
            "bystander_lt": _extract_dec_from_cell_threshold(t02, r"bystander_ct_pct<([0-9.eE+-]+)", "D1-02:bystander"),
            "chr6_del_eq": _extract_dec_from_cell_threshold(t02, r"chr6p21_interlocus_del_pct==([0-9.eE+-]+)", "D1-02:chr6_del"),
            "chrm_depth_min": _extract_dec_from_cell_threshold(t03, r"chrm_numt_amplicon_depth_x>=([0-9.eE+-]+)", "D1-03:chrm_depth"),
            "chrm_het_lt": _extract_dec_from_cell_threshold(t03, r"chrm_novel_nonsyn_het_pct<([0-9.eE+-]+)", "D1-03:chrm_het"),
            "htmpk_leak_max": _extract_dec_from_cell_threshold(t04, r"G0_beta_htmpk_leak_pct<=([0-9.eE+-]+)", "D1-04:htmpk_leak"),
            "elispot_max": _extract_dec_from_cell_threshold(t04, r"hlac07_ryafsgvay_elispot_above_bg_pct<=([0-9.eE+-]+)", "D1-04:elispot"),
            "ap1903_u95_max": _extract_dec_from_cell_threshold(t04, r"AP1903_k==0_exact_80digit_sound_Decimal_U95<=([0-9.eE+-]+)", "D1-04:ap1903_u95"),
            "azt_forced_chr19_loss_required": ("azt_kill_measured_after_forced_chr19_loss==True" in t04),
            "azt_kill_min": _extract_dec_from_cell_threshold(t04, r"AZT_dividing_log10_kill>=([0-9.eE+-]+)", "D1-04:azt_kill"),
            "l3_ret_min": _extract_dec_from_cell_threshold(t05, r"L3_stage7_day28_retention_pct>=([0-9.eE+-]+)", "D1-05:l3_ret"),
            "l4_ret_min": _extract_dec_from_cell_threshold(t05, r"L4_stage7_day28_retention_pct>=([0-9.eE+-]+)", "D1-05:l4_ret"),
            "min_l3_l4_ret_min": _extract_dec_from_cell_threshold(t05, r"min_L3_L4_stage7_day28_retention_pct>=([0-9.eE+-]+)", "D1-05:min_l3_l4"),
            "fsp1_surv_min": _extract_dec_from_cell_threshold(t05, r"RSL3_fsp1_ferroptosis_survival_pct>=([0-9.eE+-]+)", "D1-05:fsp1"),
            "kir_nk_max": _extract_dec_from_cell_threshold(t06, r"kir3dl1_missing_bw4_autologous_nk_lysis_pct<=([0-9.eE+-]+)", "D1-06:kir_nk"),
            "vtcn1_shed_max": _extract_dec_from_cell_threshold(t07, r"vtcn1_nrd1_soluble_sb7h4_shedding_pct<=([0-9.eE+-]+)", "D1-07:vtcn1"),
            "hlac03_u95_pct_max": _extract_dec_from_cell_threshold(t08, r"HLAC0304_k==0_exact_80digit_sound_Decimal_U95<=([0-9.eE+-]+)%", "D1-08:hlac03_u95"),
            "t1d_c5b9_u95_pct_max": _extract_dec_from_cell_threshold(t08, r"T1D_serum_C5b9_k==0_exact_80digit_sound_Decimal_U95<=([0-9.eE+-]+)%", "D1-08:c5b9_u95"),
            "d1_08_two_branch_enforced": (
                "d1_08_patient_carries_hlac0304_bool_or_derived_from_patient_hla_c_alleles" in t08
                and "hlac0304_positive_control_valid==True" in t08
                and "matched_allele_positive_control_valid==True" in t08
                and "matched_allele_cd8_screen_passed==True" in t08
            ),
            "turnaround_max": _extract_dec_from_cell_threshold(t09, r"autologous_sequential_6locus_gmp_turnaround_weeks<=([0-9.eE+-]+)", "D1-09:turnaround"),
        }
        if not thresholds["azt_forced_chr19_loss_required"]:
            raise ValueError(f"Controlling 45-cell contract D1-04 missing required forced chr19 loss clause: {t04!r}")
        if not thresholds["d1_08_two_branch_enforced"]:
            raise ValueError(f"Controlling 45-cell contract D1-08 missing required two-branch & positive-control clauses: {t08!r}")
        if "REQUIRED_EXACT_PRODUCT_WITH_LOT_MATCHED_GATE0V" == pvid_cells["D1-09"].get("disposition"):
            thresholds["g0v_day14_min"] = _extract_dec_from_cell_threshold(t09, r"Gate0V_day14_log10_reduction>=([0-9.eE+-]+)", "D1-09:g0v_day14")
            thresholds["g0v_day14_l95_min"] = _extract_dec_from_cell_threshold(t09, r"Gate0V_day14_log10_reduction_l95>=([0-9.eE+-]+)", "D1-09:g0v_day14_l95")
            thresholds["g0v_72h_max"] = _extract_dec_from_cell_threshold(t09, r"Gate0V_cvb4_72h_early_titer_log10_pfu_u95<=([0-9.eE+-]+)", "D1-09:g0v_72h")
            thresholds["g0v_day28_u95_max"] = _extract_dec_from_cell_threshold(t09, r"Gate0V_day28_VP1_IE1_k==0_U95<=([0-9.eE+-]+)", "D1-09:g0v_u95")
            m_n = re.search(r"Gate0V_day28_VP1_IE1_k==0_U95<=[0-9.eE+-]+\(n>=([0-9]+)\)", t09, re.ASCII)
            if not m_n:
                raise ValueError(f"Controlling 45-cell contract missing Gate0V n>= threshold in {t09!r}")
            thresholds["g0v_day28_n_min"] = int(m_n.group(1))
    elif len(pvid_cells) == 9 and is_component_excluded:
        if pvid == "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9":
            thresholds = {
                "ap1903_u95_max": _extract_dec_from_cell_threshold(str(pvid_cells["D1-04"]["threshold"]), r"L3_chr19_AAVS1_AP1903_k==0_exact_80digit_sound_Decimal_U95<=([0-9.eE+-]+)", "L3:D1-04:ap1903_u95"),
                "l3_ret_min": _extract_dec_from_cell_threshold(str(pvid_cells["D1-05"]["threshold"]), r"L3_stage7_day28_retention_pct>=([0-9.eE+-]+)", "L3:D1-05:l3_ret"),
                "vtcn1_shed_max": _extract_dec_from_cell_threshold(str(pvid_cells["D1-07"]["threshold"]), r"L3_vtcn1_nrd1_soluble_sb7h4_shedding_pct<=([0-9.eE+-]+)", "L3:D1-07:vtcn1"),
            }
        elif pvid == "ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y":
            t04_c = str(pvid_cells["D1-04"]["threshold"])
            t05_c = str(pvid_cells["D1-05"]["threshold"])
            thresholds = {
                "htmpk_leak_max": _extract_dec_from_cell_threshold(t04_c, r"L4_G0_beta_htmpk_leak_pct<=([0-9.eE+-]+)", "L4:D1-04:htmpk_leak"),
                "elispot_max": _extract_dec_from_cell_threshold(t04_c, r"L4_hlac07_ryafsgvay_elispot_above_bg_pct<=([0-9.eE+-]+)", "L4:D1-04:elispot"),
                "azt_forced_chr19_loss_required": ("L4_azt_kill_measured_after_forced_chr19_loss==True" in t04_c),
                "azt_kill_min": _extract_dec_from_cell_threshold(t04_c, r"L4_AZT_dividing_log10_kill>=([0-9.eE+-]+)", "L4:D1-04:azt_kill"),
                "l4_ret_min": _extract_dec_from_cell_threshold(t05_c, r"L4_stage7_day28_retention_pct>=([0-9.eE+-]+)", "L4:D1-05:l4_ret"),
                "fsp1_surv_min": _extract_dec_from_cell_threshold(t05_c, r"L4_RSL3_fsp1_ferroptosis_survival_pct>=([0-9.eE+-]+)", "L4:D1-05:fsp1"),
            }
            if not thresholds["azt_forced_chr19_loss_required"]:
                raise ValueError(f"Controlling 45-cell contract L4:D1-04 missing required forced chr19 loss clause: {t04_c!r}")
        elif pvid == "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1":
            # Parse L2's complete normative D1-08 component subgate STRICTLY from L2's own cell (ZERO L1A borrowing or numeric fallback!)
            t08_l2 = str(pvid_cells["D1-08"]["threshold"])
            l2_two_branch_ok = (
                "L2_d1_08_patient_carries_hlac0304_bool_or_derived_from_patient_hla_c_alleles" in t08_l2
                and "hlac0304_positive_control_valid==True" in t08_l2
                and "matched_allele_positive_control_valid==True" in t08_l2
                and "matched_allele_cd8_screen_passed==True" in t08_l2
                and "optional_adverse_T1D_serum_C5b9_k>=1_triggers_FAIL_HARD_STOP" in t08_l2
            )
            if not l2_two_branch_ok:
                raise ValueError(f"Controlling 45-cell contract L2:D1-08 missing required self-contained branch/adverse clauses: {t08_l2!r}")
            thresholds = {
                "hlac03_u95_pct_max": _extract_dec_from_cell_threshold(t08_l2, r"L2_HLAC0304_k==0_exact_80digit_sound_Decimal_U95<=([0-9.eE+-]+)%", "L2:D1-08:hlac03_u95"),
                "l2_d1_08_self_contained_enforced": True,
            }

    return {
        "gate_id": str(blk.get("gate_id") or ""),
        "schema_version": str(pbg.get("schema_version") or ""),
        "sha256": dt_sha,
        "staged_registry_verified": staged_registry_verified,
        "canonical_sha_verified": bool(dt_sha == EXPECTED_DEC114_DECISION_TREE_SHA256 and on_disk_dt_verified),
        "on_disk_decision_tree_verified": on_disk_dt_verified,
        "in_memory_block_matches_canonical_sha": bool(dt_sha == EXPECTED_DEC114_DECISION_TREE_SHA256),
        "playbook_sec4120_sha256": pb_info["pb_sha"],
        "playbook_footer_present": pb_info["footer_present"],
        "playbook_footer_verified": pb_info["footer_verified"],
        "playbook_canonical_sha_verified": pb_info["canonical_sha_verified"],
        "cells_loaded": len(pvid_cells),
        "external_refs": [str(pvid_cells[g].get("external_ref") or "") for g in EXPECTED_45_GATES if g in pvid_cells],
        "cell_missing_dispositions": cell_missing_dispositions,
        "is_component_excluded": is_component_excluded,
        "thresholds": thresholds,
    }


def _aggregate_panel_coverage_regime(cal_summary: Dict[str, Any]) -> str:
    """Aggregate per-assay coverage categories independently of distinct Bonferroni percentages."""
    active_entries = [m for m in cal_summary.values() if isinstance(m, dict) and isinstance(m.get("coverage_regime"), str)]
    if not active_entries:
        return "UNAVAILABLE"
    regimes = sorted({str(m["coverage_regime"]) for m in active_entries})
    if len(regimes) == 1:
        return regimes[0]
    cats = {str(m.get("coverage_category") or ("CALIBRATED" if "BONFERRONI" in str(m["coverage_regime"]) else "UNIT_ASSUMED")) for m in active_entries}
    if cats == {"CALIBRATED"}:
        return "ALL_CALIBRATED_VARIED_BONFERRONI_COVERAGE:" + "|".join(regimes)
    if cats == {"UNIT_ASSUMED"}:
        return "ALL_UNIT_ASSUMED_VARIED_OBSERVATION_ALPHA:" + "|".join(regimes)
    if cats == {"CALIBRATED", "UNIT_ASSUMED"}:
        return "MIXED_CALIBRATED_AND_UNIT_ASSUMED:" + "|".join(regimes)
    if "CALIBRATED" not in cats:
        return "MIXED_PARTIAL_AND_UNIT_ASSUMED:" + "|".join(regimes)
    return "MIXED_CALIBRATED_AND_ASSUMPTION_DEPENDENT:" + "|".join(regimes)


def _collect_assumption_only_flags(cal_summary: Dict[str, Any]) -> List[str]:
    """Surface EVERY default or uncalibrated sampling/recovery/eligibility assumption visibly at headline output."""
    flags: List[str] = []
    for _ak, _am in cal_summary.items():
        if not isinstance(_am, dict):
            continue
        _sp = _am.get("sampling_independence_provenance")
        if _sp in ("UNSPECIFIED_DEFAULT_ASSUMED_INDEPENDENT", "EXPLICIT_NON_CLUSTERED_ONLY_INDEPENDENCE_ASSUMED"):
            flags.append(f"{_ak}:SAMPLING_INDEPENDENCE_ASSUMED({_sp})")
        _cat = _am.get("coverage_category")
        if _cat == "UNIT_ASSUMED":
            flags.append(f"{_ak}:UNIT_RECOVERY_AND_ELIGIBILITY_ASSUMED")
        else:
            _da = _am.get("default_calibration_assumptions")
            if isinstance(_da, list):
                for _item in _da:
                    flags.append(f"{_ak}:{_item}")
    return flags


def evaluate_stage0_d1b_auto_v2_gates_iip0035(
    observations: Dict[str, Any],
    *,
    product_def: Dict[str, Any] | None = None,
    product_variant_id: str | None = None,
    declared_product_variant_id: str | None = None,
    evaluated_lot_id: str | None = None,
    extra_lot_bindings: Dict[str, Any] | None = None,
    global_gate0v_record: Any = _MISSING_ARG,
    outer_gate0v_explicit_null: bool = False,
    decision_tree_block: Dict[str, Any] | None = None,
    repo_root: Path = ROOT,
    require_canonical_sha: bool = True,
) -> Dict[str, Any]:
    """Standalone fail-closed Stage-0 reference evaluator driven directly by the 45-cell DEC-114 contract (D1-01..D1-09).

    Zero wetlab_contract.py dependency. Enforces unified integrity across equivalent observation envelopes
    under default strict mode (verifying both stage0_decision_tree_v1.json and CommonMark Playbook §4.1.20
    SHA-256s), parses numeric/branch rules directly from the 45-cell contract, reconciles conflicting
    product identity declarations symmetrically, and explicitly distinguishes reference-scenario arithmetic
    (keyed to a single-locus L1 donor record) from assembled multi-locus manifests or clinical/batch release
    authorization (`clinical_or_batch_release_authorized` is always False).
    """
    invalid_or_oob_keys: List[str] = []
    explicitly_unavailable: List[str] = []
    if not isinstance(observations, dict):
        return {
            "gate_id": "STAGE0_GATE_AUTOLOGOUS_D1B_AUTO_C_V2_D101_D109_IIP_0035",
            "decision_id": "DEC-114",
            "contribution_id": "IIP-0035",
            "evaluation_kind": "REFERENCE_CONTRACT_COMPUTATION_NOT_EMPIRICAL_RELEASE_DECISION",
            "clinical_or_batch_release_authorized": False,
            "assembled_multilocus_manifest_bound": False,
            "empirical_clearance_status": "HOLD_INVALID_OBSERVATION_PAYLOAD",
            "reference_criteria_met": False,
            "status": "FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION",
            "reasons": ["FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION"],
        }

    def _dec(key: str) -> Decimal | None:
        raw = observations.get(key)
        if isinstance(raw, str) and raw == "UNAVAILABLE":
            explicitly_unavailable.append(key)
            return None
        if raw is None:
            return None
        dv = _safe_exact_decimal(raw)
        if dv is None:
            invalid_or_oob_keys.append(key)
            return None
        bounds = _D1_ENDPOINT_BOUNDS.get(key)
        if bounds is not None:
            if dv < Decimal(str(bounds[0])) or (bounds[1] is not None and dv > Decimal(str(bounds[1]))):
                invalid_or_oob_keys.append(key)
                return None
        return dv

    def _strict_int(key: str) -> int | None:
        raw = observations.get(key)
        if isinstance(raw, str) and raw == "UNAVAILABLE":
            explicitly_unavailable.append(key)
            return None
        if raw is None:
            return None
        if isinstance(raw, bool) or not isinstance(raw, int):
            invalid_or_oob_keys.append(key)
            return None
        if raw < 0:
            invalid_or_oob_keys.append(key)
            return None
        return raw

    reasons: List[str] = []
    obs_pdef_inner = observations.get("product_def")
    has_pdef_obj_conflict = False
    if product_def is not None and not isinstance(product_def, dict):
        reasons.append("FAIL_D1_INVALID_OR_TAMPERED_PRODUCT_DEFINITION")
    if obs_pdef_inner is not None and not isinstance(obs_pdef_inner, dict):
        reasons.append("FAIL_D1_CONFLICTING_PRODUCT_IDENTITY_DECLARATIONS")
        has_pdef_obj_conflict = True
    elif isinstance(product_def, dict) and isinstance(obs_pdef_inner, dict) and not _strict_typed_equal(obs_pdef_inner, product_def):
        reasons.append("FAIL_D1_CONFLICTING_PRODUCT_IDENTITY_DECLARATIONS")
        has_pdef_obj_conflict = True

    # Feed embedded product_variant_id from EVERY supplied product_def (both outer and inner) symmetrically!
    outer_pdef_id = product_def.get("product_variant_id") if isinstance(product_def, dict) else None
    inner_pdef_id = obs_pdef_inner.get("product_variant_id") if isinstance(obs_pdef_inner, dict) else None
    declared_pvids: List[str] = []
    pvid_invalid = False
    for _cand_pvid in (outer_pdef_id, inner_pdef_id, product_variant_id, declared_product_variant_id, observations.get("product_variant_id")):
        if _cand_pvid is None:
            continue
        if _is_canonical_identifier(_cand_pvid):
            if _cand_pvid not in declared_pvids:
                declared_pvids.append(_cand_pvid)
        else:
            pvid_invalid = True
            invalid_or_oob_keys.append("product_variant_id")
    has_pvid_conflict = len(declared_pvids) > 1 or has_pdef_obj_conflict
    if len(declared_pvids) > 1:
        reasons.append("FAIL_D1_CONFLICTING_PRODUCT_IDENTITY_DECLARATIONS")

    eff_product_def = product_def if isinstance(product_def, dict) else (obs_pdef_inner if isinstance(obs_pdef_inner, dict) else None)
    # Under ANY identity or object conflict, keep identity symmetrically unresolved (order-independent!).
    if has_pvid_conflict:
        pvid = ""
        eff_product_def = None
    elif eff_product_def is None and len(declared_pvids) == 1 and not pvid_invalid:
        _loaded_pc = pc.load_product_contracts(repo_root, require_staged_fastas=require_canonical_sha)
        _pc_by_id = {p["product_variant_id"]: p for p in _loaded_pc["products"]}
        _single_pvid = declared_pvids[0]
        if _single_pvid in _pc_by_id:
            eff_product_def = _pc_by_id[_single_pvid]
            pvid = _single_pvid
        else:
            pvid = _single_pvid
    else:
        pvid = declared_pvids[0] if (len(declared_pvids) == 1 and not pvid_invalid) else ""
    product_def = eff_product_def
    contract = _load_controlling_contract_for_product(
        pvid,
        decision_tree_block=decision_tree_block,
        repo_root=repo_root,
        require_canonical_sha=require_canonical_sha,
    )
    thr = contract["thresholds"]
    cell_mdd = contract["cell_missing_dispositions"]
    is_fullstack = pvid in D1B_FULLSTACK_CANDIDATE_IDS
    is_component = (pvid in D1B_COMPONENT_DONOR_IDS) or contract["is_component_excluded"]

    if not pvid:
        if not has_pvid_conflict:
            reasons.append("HOLD_D1_CANONICAL_PRODUCT_VARIANT_ID_REQUIRED")
    elif not (is_fullstack or is_component):
        reasons.append("FAIL_D1_CANONICAL_PRODUCT_NOT_D1B_AUTO_V2")
    elif not isinstance(product_def, dict):
        reasons.append("FAIL_D1_INVALID_OR_TAMPERED_PRODUCT_DEFINITION")
    else:
        if is_component:
            reasons.append("HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE")
        exp_rec = pc.D1B_AUTO_V2_DONOR_BINDINGS.get(pvid) or {}
        try:
            canon_pdef = _canonicalize_json_for_digest(product_def)
            got_p_sha = hashlib.sha256(
                json.dumps(canon_pdef, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
            ).hexdigest()
        except Exception:
            got_p_sha = ""
        if (
            product_def.get("execution_supported") is not False
            or product_def.get("physical_route_id") != "ROUTE_3P_MODULAR_PLANAR_ENVELOPE"
            or product_def.get("lifecycle_state") != "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"
            or not exp_rec.get("full_record_sha256")
            or got_p_sha != exp_rec["full_record_sha256"]
        ):
            reasons.append("FAIL_D1_INVALID_OR_TAMPERED_PRODUCT_DEFINITION")

    canon_l1_arm: str | None = None
    if pvid == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN":
        canon_l1_arm = "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"
    elif pvid == "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR":
        canon_l1_arm = "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"

    l1_arm = observations.get("d1_auto_locus1_donor_arm")
    if isinstance(l1_arm, str) and l1_arm == "UNAVAILABLE":
        explicitly_unavailable.append("d1_auto_locus1_donor_arm")
        l1_arm = None
    elif l1_arm is not None and not isinstance(l1_arm, str):
        invalid_or_oob_keys.append("d1_auto_locus1_donor_arm")
        l1_arm = None
    elif l1_arm is not None and l1_arm not in ALLOWED_LOCUS1_ARMS:
        reasons.append("HOLD_D1_01_UNKNOWN_LOCUS1_DONOR_ARM")
    elif canon_l1_arm is not None and l1_arm is not None and l1_arm != canon_l1_arm:
        reasons.append("FAIL_D1_01_CANONICAL_PRODUCT_LOCUS1_ARM_MISMATCH")
    effective_l1_arm = None if has_pvid_conflict else (canon_l1_arm if canon_l1_arm is not None else l1_arm)

    # Strict lot-ID binding using the single shared _LOT_ALIAS_KEY_SET across Python and CLI!
    lot_candidates: List[str] = []
    lot_invalid = False
    lot_pairs: List[tuple[str, Any]] = [
        ("evaluated_lot_id_arg", evaluated_lot_id),
        ("mcb_lot_id", observations.get("mcb_lot_id")),
        ("lot_id", observations.get("lot_id")),
        ("evaluated_lot_id", observations.get("evaluated_lot_id")),
        ("batch_id", observations.get("batch_id")),
        ("mcb_batch_id", observations.get("mcb_batch_id")),
    ]
    if extra_lot_bindings is not None:
        if not isinstance(extra_lot_bindings, dict):
            lot_invalid = True
        else:
            # Validate key types & membership in _LOT_ALIAS_KEY_SET BEFORE sorting (prevents TypeError on mixed keys & rejects unknown alias keys!)
            valid_ek_pairs: List[tuple[str, Any]] = []
            for ek, ev in extra_lot_bindings.items():
                if not isinstance(ek, str) or ek not in _LOT_ALIAS_KEY_SET:
                    lot_invalid = True
                    invalid_or_oob_keys.append("extra_lot_bindings")
                else:
                    valid_ek_pairs.append((ek, ev))
            for ek, ev in sorted(valid_ek_pairs, key=lambda kv: kv[0]):
                lot_pairs.append((f"envelope.{ek}", ev))
    for lot_src_key, lot_raw in lot_pairs:
        if lot_raw is None:
            continue
        if isinstance(lot_raw, str) and lot_raw == "UNAVAILABLE":
            explicitly_unavailable.append(lot_src_key)
            continue
        if not _is_canonical_identifier(lot_raw):
            lot_invalid = True
            continue
        if lot_raw not in lot_candidates:
            lot_candidates.append(lot_raw)
    if lot_invalid:
        reasons.append("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID")
    lot_conflict = len(lot_candidates) > 1
    if lot_conflict:
        reasons.append("HOLD_D1_CONFLICTING_LOT_ID_BINDING")
    lot_id_eff = lot_candidates[0] if (len(lot_candidates) == 1 and not lot_invalid) else None
    if lot_id_eff is None and not lot_conflict and not lot_invalid:
        reasons.append("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID")

    # Resolve global_gate0v_record symmetrically across Python (private _MISSING_ARG sentinel + explicit None == JSON null!) and CLI
    obs_g0v_present = "global_gate0v_record" in observations
    obs_g0v_rec = observations.get("global_gate0v_record")
    outer_g0v_present = (global_gate0v_record is not _MISSING_ARG) or bool(outer_gate0v_explicit_null)
    outer_g0v_val = None if global_gate0v_record is _MISSING_ARG else global_gate0v_record
    g0v_supplied_records: List[Any] = []
    if outer_g0v_present:
        g0v_supplied_records.append(outer_g0v_val)
    if obs_g0v_present:
        g0v_supplied_records.append(obs_g0v_rec)
    g0v_conflict = bool(
        outer_g0v_present
        and obs_g0v_present
        and not _strict_typed_equal(outer_g0v_val, obs_g0v_rec)
    )
    g0v_any_positive_k = any(_has_positive_gate0v_event(c) for c in g0v_supplied_records)
    if g0v_conflict:
        reasons.append("FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS")
        if g0v_any_positive_k:
            if has_pvid_conflict or not is_fullstack:
                reasons.append("FAIL_HARD_STOP_D1_UNRESOLVED_PRODUCT_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT")
            elif effective_l1_arm == "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR":
                reasons.append("FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN")
            else:
                reasons.append("FAIL_HARD_STOP_D1_L1A_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT")
    eff_g0v_record = (
        outer_g0v_val
        if (outer_g0v_present and outer_g0v_val is not None)
        else (obs_g0v_rec if (obs_g0v_present and obs_g0v_rec is not None) else None)
    )

    def _cal_kwargs(prefix: str) -> Dict[str, Any]:
        """Resolve assay-local and panel-level calibration & sampling via shared typed resolver."""
        return _resolve_assay_calibration(observations, local_prefix=prefix)

    # Upfront raw 2-field patient_hla_c_alleles parsing shared across D1-01, full-product D1-08, AND L2 component path
    raw_hlac_present = "patient_hla_c_alleles" in observations
    raw_hlac_alleles = observations.get("patient_hla_c_alleles")
    parsed_hlac_alleles: List[str] | None = None
    derived_carries_hlac0304: bool | None = None
    unlisted_hlac_blocks_d1_08_negative = False
    if not raw_hlac_present:
        d1_01_raw_hla_c_input_status = "NOT_SUPPLIED"
    elif raw_hlac_alleles is None:
        d1_01_raw_hla_c_input_status = "EXPLICIT_NULL"
    elif isinstance(raw_hlac_alleles, str) and raw_hlac_alleles == "UNAVAILABLE":
        explicitly_unavailable.append("patient_hla_c_alleles")
        d1_01_raw_hla_c_input_status = "EXPLICIT_UNAVAILABLE"
    else:
        if (
            not isinstance(raw_hlac_alleles, list)
            or len(raw_hlac_alleles) != 2
            or not all(isinstance(a, str) and a == a.strip() and re.fullmatch(r"HLA-C\*[0-9]{2}:[0-9]{2}", a, re.ASCII) for a in raw_hlac_alleles)
        ):
            invalid_or_oob_keys.append("patient_hla_c_alleles")
            d1_01_raw_hla_c_input_status = "MALFORMED"
        else:
            parsed_hlac_alleles = list(raw_hlac_alleles)
            d1_01_raw_hla_c_input_status = "VALID_PAIR"
            if "HLA-C*03:04" in parsed_hlac_alleles:
                # Positive HLA-C*03:04 evidence is ALWAYS retained regardless of partner allele status
                derived_carries_hlac0304 = True
            elif all(a in _PINNED_KNOWN_HLAC_TWO_FIELD_INVENTORY for a in parsed_hlac_alleles):
                # Both alleles are verified against the version-pinned Lin 2023 Table 1 HLA-C inventory
                derived_carries_hlac0304 = False
            else:
                # Unlisted two-field token (e.g. HLA-C*00:00 / HLA-C*99:99) MUST HOLD negative carrier inference!
                derived_carries_hlac0304 = None
                unlisted_hlac_blocks_d1_08_negative = True

    def _cal_meta_from_cp(cp_obj: Dict[str, Any] | None, fallback_samp: Dict[str, Any] | None = None) -> Dict[str, Any] | None:
        if isinstance(cp_obj, dict):
            out_m: Dict[str, Any] = {
                "coverage_regime": cp_obj.get("coverage_regime"),
                "coverage_category": cp_obj.get("coverage_category"),
                "default_calibration_assumptions": cp_obj.get("default_calibration_assumptions"),
                "sampling_independence_provenance": cp_obj.get("sampling_independence_provenance"),
                "alpha_obs_str": cp_obj.get("alpha_obs_str"),
                "alpha_joint_total_str": cp_obj.get("alpha_joint_total_str"),
                "effective_detection_factor_str": cp_obj.get("effective_detection_factor_str"),
                "calibration_sources": cp_obj.get("calibration_sources"),
                "status": cp_obj.get("status"),
            }
            return out_m
        if isinstance(fallback_samp, dict):
            return {
                "coverage_regime": None,
                "coverage_category": None,
                "default_calibration_assumptions": None,
                "sampling_independence_provenance": fallback_samp.get("sampling_independence_provenance"),
                "alpha_obs_str": None,
                "alpha_joint_total_str": None,
                "effective_detection_factor_str": None,
                "calibration_sources": None,
                "status": "ARITHMETIC_WITHHELD",
            }
        return None

    # Clean locus-specific component-only evaluation path (carries full per-assay bound provenance & product-neutral Gate-0V adverse status!)
    if is_component:
        comp_gates = sorted(COMPONENT_APPLICABLE_GATES.get(pvid, ()))
        comp_reasons: List[str] = []
        comp_ap_cp: Dict[str, Any] | None = None
        comp_hlac_cp: Dict[str, Any] | None = None
        comp_t1d_cp: Dict[str, Any] | None = None
        comp_samp_meta: Dict[str, Any] = {}
        comp_d1_08_branch = "COMPONENT_ONLY_NOT_FULLSTACK"
        if pvid == "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9":
            ap_u95_cap_c = thr["ap1903_u95_max"]
            ap_surv_c = _dec("obs_d1_04_dec053_ap1903_icasp9_survivor_frac")
            ap_k_c = _strict_int("obs_d1_04_ap1903_survivor_count_k")
            ap_n_r = _strict_int("obs_d1_04_ap1903_raw_trials_n")
            ap_n_e = _strict_int("obs_d1_04_ap1903_effective_trials_n")
            if ap_n_r is not None and ap_n_e is not None and ap_n_r != ap_n_e:
                comp_reasons.append("FAIL_D1_04_CONFLICTING_RAW_AND_EFFECTIVE_TRIALS_N")
            ap_n_c = ap_n_r if ap_n_r is not None else ap_n_e
            if ap_surv_c is not None and ap_surv_c > ap_u95_cap_c:
                comp_reasons.append("FAIL_D1_04_DEC053_AP1903_ICASP9_SURVIVOR_FRAC_EXCEEDS_1E_8")
            _l3_ap_cal = _cal_kwargs("obs_d1_04_ap1903")
            comp_samp_meta["d1_04_ap1903"] = _l3_ap_cal
            if ap_k_c is None or ap_n_c is None or ap_n_c < 1 or ap_k_c > ap_n_c:
                if ap_k_c is not None and ap_k_c >= 1:
                    comp_reasons.append("FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED")
                else:
                    comp_reasons.append("HOLD_D1_04_AP1903_GY018_INTEGER_KN_REQUIRED")
            else:
                comp_ap_cp = evaluate_sound_decimal_clopper_pearson_zero_event(
                    gate_id="D1_04_AP1903_ICASP9", k=ap_k_c, n_raw=ap_n_c,
                    max_u95_fraction=ap_u95_cap_c, **_l3_ap_cal,
                )
                if comp_ap_cp["status"] == "FAIL_HARD_STOP":
                    comp_reasons.append("FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED")
                elif comp_ap_cp["status"] != "PASS":
                    comp_reasons.append(
                        "HOLD_D1_04_AP1903_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                        if comp_ap_cp["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                        else f"HOLD_D1_04_AP1903_{comp_ap_cp['status']}"
                    )
                if ap_surv_c is not None:
                    with _clean_dec_ctx(40):
                        ap_exp_c = Decimal(ap_k_c) / Decimal(ap_n_c)
                        ap_tol_c = max(Decimal("0.05") * ap_exp_c, Decimal("1e-15"))
                        if (
                            (ap_k_c == 0 and ap_surv_c != Decimal(0))
                            or (ap_k_c > 0 and ap_surv_c <= Decimal(0))
                            or (ap_k_c > 0 and abs(ap_surv_c - ap_exp_c) > ap_tol_c)
                        ):
                            comp_reasons.append("FAIL_D1_04_AP1903_SURVIVOR_FRAC_INCONSISTENT_WITH_KN")
            l3_c = _dec("obs_d1_05_l3_day28_retention_pct")
            if l3_c is None:
                comp_reasons.append("HOLD_D1_05_L3_STAGE7_DAY28_RETENTION_REQUIRED")
            elif l3_c < thr["l3_ret_min"]:
                comp_reasons.append("FAIL_D1_05_L3_STAGE7_DAY28_RETENTION_BELOW_80PCT")
            vtcn_c = _dec("obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct")
            if vtcn_c is None:
                comp_reasons.append("HOLD_D1_07_MISSING_VTCN1_SHEDDING_ASSAY")
            elif vtcn_c > thr["vtcn1_shed_max"]:
                comp_reasons.append("FAIL_D1_07_VTCN1_NRD1_SB7H4_SHEDDING_EXCEEDS_10PCT_MUTAGENESIS_REQUIRED")
        elif pvid == "ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y":
            chr13_c = observations.get("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified")
            if isinstance(chr13_c, str) and chr13_c == "UNAVAILABLE":
                explicitly_unavailable.append("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified")
                chr13_c = None
            elif chr13_c is not None and not isinstance(chr13_c, bool):
                invalid_or_oob_keys.append("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified")
                chr13_c = None
            elif chr13_c is False:
                comp_reasons.append("FAIL_D1_04_CHR13_CLYBL_TANDEM_INDEPENDENCE_FROM_CHR19_UNVERIFIED")
            if chr13_c is None:
                comp_reasons.append("HOLD_D1_04_L4_CHR13_INDEPENDENCE_REQUIRED")
            azt_chr19_c = observations.get("d1_04_azt_kill_measured_after_forced_chr19_loss")
            if isinstance(azt_chr19_c, str) and azt_chr19_c == "UNAVAILABLE":
                explicitly_unavailable.append("d1_04_azt_kill_measured_after_forced_chr19_loss")
                azt_chr19_c = None
            elif azt_chr19_c is not None and not isinstance(azt_chr19_c, bool):
                invalid_or_oob_keys.append("d1_04_azt_kill_measured_after_forced_chr19_loss")
                azt_chr19_c = None
            elif azt_chr19_c is False:
                comp_reasons.append("FAIL_D1_04_AZT_KILL_NOT_MEASURED_AFTER_FORCED_CHR19_LOSS")
            if azt_chr19_c is None:
                comp_reasons.append("HOLD_D1_04_AZT_KILL_FORCED_CHR19_LOSS_CONDITION_REQUIRED")
            leak_c = _dec("obs_d1_04_g0_beta_htmpk_leak_pct")
            elis_c = _dec("obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct")
            azt_c = _dec("obs_d1_04_dec053_azt_dividing_log10_kill")
            l4_c = _dec("obs_d1_05_l4_day28_retention_pct")
            fsp_c = _dec("obs_d1_05_fsp1_ferroptosis_survival_pct")
            if leak_c is not None and leak_c > thr["htmpk_leak_max"]:
                comp_reasons.append("FAIL_D1_04_G0_BETA_HTMPK_F105Y_LEAK_EXCEEDS_0P1PCT")
            if elis_c is not None and elis_c > thr["elispot_max"]:
                comp_reasons.append("FAIL_D1_04_HLAC07_RYAFSGVAY_ELISPOT_EXCEEDS_0P1PCT_G0_GATE")
            if azt_c is not None and azt_c < thr["azt_kill_min"]:
                comp_reasons.append("FAIL_D1_04_DEC053_AZT_DIVIDING_KILL_BELOW_4LOG_AFTER_CHR19_LOSS")
            if l4_c is not None and l4_c < thr["l4_ret_min"]:
                comp_reasons.append("FAIL_D1_05_L4_STAGE7_DAY28_RETENTION_BELOW_80PCT")
            if fsp_c is not None and fsp_c < thr["fsp1_surv_min"]:
                comp_reasons.append("FAIL_D1_05_FSP1_FERROPTOSIS_SURVIVAL_BELOW_80PCT")
            if any(x is None for x in (leak_c, elis_c, azt_c, l4_c, fsp_c)) or azt_chr19_c is not True:
                comp_reasons.append("HOLD_D1_L4_COMPONENT_SUBGATES_INCOMPLETE")
        elif pvid == "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1":
            hlac03_cap_c = thr["hlac03_u95_pct_max"]
            carries_c_present = "d1_08_patient_carries_hlac0304" in observations
            carries_c = observations.get("d1_08_patient_carries_hlac0304")
            carries_c_raw_explicit_unavail_or_invalid = False
            if carries_c_present and carries_c is None:
                carries_c_raw_explicit_unavail_or_invalid = True
            elif isinstance(carries_c, str) and carries_c == "UNAVAILABLE":
                explicitly_unavailable.append("d1_08_patient_carries_hlac0304")
                carries_c = None
                carries_c_raw_explicit_unavail_or_invalid = True
            elif carries_c is not None and not isinstance(carries_c, bool):
                invalid_or_oob_keys.append("d1_08_patient_carries_hlac0304")
                carries_c = None
                carries_c_raw_explicit_unavail_or_invalid = True
            if derived_carries_hlac0304 is not None:
                if isinstance(carries_c, bool) and carries_c != derived_carries_hlac0304:
                    comp_reasons.append("FAIL_D1_08_HLAC0304_CARRIER_STATUS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE")
                    if derived_carries_hlac0304 is True:
                        carries_c = True
                elif carries_c is None and not carries_c_raw_explicit_unavail_or_invalid:
                    carries_c = derived_carries_hlac0304
            elif unlisted_hlac_blocks_d1_08_negative and carries_c is None and not carries_c_raw_explicit_unavail_or_invalid:
                comp_reasons.append("HOLD_D1_08_UNLISTED_RAW_HLAC_ALLELE_BLOCKS_NEGATIVE_DERIVATION")
            hlac03_t_c = _dec("obs_d1_08_hlac0304_tcell_lysis_u95_pct")
            hlac03_k_c = _strict_int("obs_d1_08_hlac0304_lysed_events_k")
            hlac03_n_c = _strict_int("obs_d1_08_hlac0304_scored_targets_n")
            t1d_k_c_extra = _strict_int("obs_d1_08_t1d_serum_c5b9_pos_events_k")
            if t1d_k_c_extra is not None and t1d_k_c_extra >= 1:
                comp_reasons.append("FAIL_HARD_STOP_D1_08_T1D_SERUM_C5B9_EVENT_DETECTED")
            _l2_hlac_cal = _cal_kwargs("obs_d1_08_hlac0304")
            comp_samp_meta["d1_08_hlac0304"] = _l2_hlac_cal
            if carries_c is True:
                comp_d1_08_branch = "L2_COMPONENT_HLAC0304_CARRIER_CD8_SUBGATE"
                c0304_pc_c = observations.get("d1_08_hlac0304_positive_control_valid")
                if isinstance(c0304_pc_c, str) and c0304_pc_c == "UNAVAILABLE":
                    explicitly_unavailable.append("d1_08_hlac0304_positive_control_valid")
                    c0304_pc_c = None
                elif c0304_pc_c is not None and not isinstance(c0304_pc_c, bool):
                    invalid_or_oob_keys.append("d1_08_hlac0304_positive_control_valid")
                    c0304_pc_c = None
                elif c0304_pc_c is False:
                    comp_reasons.append("FAIL_D1_08_HLAC0304_POSITIVE_CONTROL_INVALID")
                if c0304_pc_c is None:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_POSITIVE_CONTROL_VALIDITY_REQUIRED")
                if hlac03_k_c is None or hlac03_n_c is None or hlac03_n_c < 1 or hlac03_k_c > hlac03_n_c:
                    if hlac03_k_c is not None and hlac03_k_c >= 1:
                        comp_reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
                    else:
                        comp_reasons.append("HOLD_D1_08_HLAC0304_GY018_INTEGER_KN_REQUIRED")
                        if hlac03_t_c is None:
                            comp_reasons.append("HOLD_D1_08_HLAC0304_REPORTED_U95_REQUIRED")
                else:
                    comp_hlac_cp = evaluate_sound_decimal_clopper_pearson_zero_event(
                        gate_id="D1_08_HLAC0304_TCELL_LYSIS", k=hlac03_k_c, n_raw=hlac03_n_c,
                        max_u95_fraction=_exact_pct_to_frac(hlac03_cap_c), reported_u95_pct=hlac03_t_c,
                        **_l2_hlac_cal,
                    )
                    if comp_hlac_cp["status"] == "FAIL_HARD_STOP":
                        comp_reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
                    elif comp_hlac_cp["status"] != "PASS":
                        comp_reasons.append(
                            "HOLD_D1_08_HLAC0304_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                            if comp_hlac_cp["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                            else f"HOLD_D1_08_HLAC0304_{comp_hlac_cp['status']}"
                        )
                        if (
                            hlac03_t_c is not None
                            and hlac03_k_c == 0
                            and hlac03_t_c <= hlac03_cap_c
                            and comp_hlac_cp.get("reported_u95_meets_exact_bound") is False
                        ):
                            comp_reasons.append("FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND")
                    else:
                        if hlac03_t_c is None:
                            comp_reasons.append("HOLD_D1_08_HLAC0304_REPORTED_U95_REQUIRED")
                        elif hlac03_t_c > hlac03_cap_c:
                            comp_reasons.append("HOLD_D1_08_HLAC0304_REPORTED_U95_ABOVE_GATE_LIMIT")
                        elif hlac03_t_c <= Decimal(0) or comp_hlac_cp.get("reported_u95_meets_exact_bound") is False:
                            comp_reasons.append("FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND")
                        elif comp_hlac_cp.get("reported_u95_meets_exact_bound") is None:
                            comp_reasons.append("HOLD_D1_08_REPORTED_U95_EXACT_FLOOR_INCONCLUSIVE")
            elif carries_c is False:
                comp_d1_08_branch = "L2_COMPONENT_NON_HLAC0304_MATCHED_ALLELE_CD8_SUBGATE"
                ma_pc_c = observations.get("d1_08_matched_allele_positive_control_valid")
                ma_ok_c = observations.get("d1_08_non_c0304_matched_allele_screen_passed")
                if isinstance(ma_pc_c, str) and ma_pc_c == "UNAVAILABLE":
                    explicitly_unavailable.append("d1_08_matched_allele_positive_control_valid")
                    ma_pc_c = None
                elif ma_pc_c is not None and not isinstance(ma_pc_c, bool):
                    invalid_or_oob_keys.append("d1_08_matched_allele_positive_control_valid")
                    ma_pc_c = None
                elif ma_pc_c is False:
                    comp_reasons.append("FAIL_D1_08_MATCHED_ALLELE_POSITIVE_CONTROL_INVALID")
                if ma_pc_c is None:
                    comp_reasons.append("HOLD_D1_08_MATCHED_ALLELE_POSITIVE_CONTROL_VALIDITY_REQUIRED")
                if isinstance(ma_ok_c, str) and ma_ok_c == "UNAVAILABLE":
                    explicitly_unavailable.append("d1_08_non_c0304_matched_allele_screen_passed")
                    ma_ok_c = None
                elif ma_ok_c is not None and not isinstance(ma_ok_c, bool):
                    invalid_or_oob_keys.append("d1_08_non_c0304_matched_allele_screen_passed")
                    ma_ok_c = None
                elif ma_ok_c is False:
                    comp_reasons.append("FAIL_D1_08_MATCHED_ALLELE_CD8_SCREEN_FAILED")
                if ma_ok_c is None:
                    comp_reasons.append("HOLD_D1_08_MATCHED_ALLELE_CD8_SCREEN_REQUIRED")
                if hlac03_k_c is not None and hlac03_k_c >= 1:
                    comp_reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
                if _l2_hlac_cal["sampling_independence_conflict"]:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_HOLD_CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS")
                elif _l2_hlac_cal["sampling_malformed"]:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_HOLD_MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION")
                elif _l2_hlac_cal["clustered_sampling"]:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED")
            else:
                comp_d1_08_branch = "L2_COMPONENT_UNSPECIFIED_HLAC0304_CARRIER_STATUS_HOLD"
                comp_reasons.append("HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED")
                if hlac03_k_c is not None and hlac03_k_c >= 1:
                    comp_reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
                if _l2_hlac_cal["sampling_independence_conflict"]:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_HOLD_CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS")
                elif _l2_hlac_cal["sampling_malformed"]:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_HOLD_MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION")
                elif _l2_hlac_cal["clustered_sampling"]:
                    comp_reasons.append("HOLD_D1_08_HLAC0304_HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED")
        # Never drop any reported cross-gate adverse safety hard stop (k >= 1) on any component donor path!
        _comp_ap_k_adv = _strict_int("obs_d1_04_ap1903_survivor_count_k")
        if _comp_ap_k_adv is not None and _comp_ap_k_adv >= 1:
            comp_reasons.append("FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED")
        _comp_hlac_k_adv = _strict_int("obs_d1_08_hlac0304_lysed_events_k")
        if _comp_hlac_k_adv is not None and _comp_hlac_k_adv >= 1:
            comp_reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
        _comp_c5b9_k_adv = _strict_int("obs_d1_08_t1d_serum_c5b9_pos_events_k")
        if _comp_c5b9_k_adv is not None and _comp_c5b9_k_adv >= 1:
            comp_reasons.append("FAIL_HARD_STOP_D1_08_T1D_SERUM_C5B9_EVENT_DETECTED")
        comp_out_of_scope_g0v_status: str | None = None
        comp_g0v_standalone_status = "NOT_APPLICABLE_COMPONENT_DONOR"
        if g0v_any_positive_k:
            comp_reasons.append("FAIL_HARD_STOP_OUT_OF_SCOPE_COMPONENT_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT")
            comp_out_of_scope_g0v_status = "FAIL_HARD_STOP_OUT_OF_SCOPE_ADVERSE_GATE0V_RECORD_RECEIVED"
            comp_g0v_standalone_status = "FAIL_HARD_STOP_OUT_OF_SCOPE_ADVERSE_GATE0V_RECORD_RECEIVED"
        if g0v_conflict:
            comp_reasons.append("FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS")
            if comp_out_of_scope_g0v_status is None:
                comp_out_of_scope_g0v_status = "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS"
                comp_g0v_standalone_status = "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS"
        # Inspect ALL supplied nested Gate-0V records on component donors for lot-alias conflicts,
        # mismatches against evaluated_lot_id, malformed lot aliases, non-dict records, or explicit FAIL verdicts!
        for _comp_g0v_rec in g0v_supplied_records:
            if _comp_g0v_rec is None:
                continue
            if isinstance(_comp_g0v_rec, bool) or not isinstance(_comp_g0v_rec, dict):
                comp_reasons.append("HOLD_D1_COMPONENT_OUT_OF_SCOPE_MALFORMED_GATE0V_RECORD")
                if comp_out_of_scope_g0v_status is None:
                    comp_out_of_scope_g0v_status = "HOLD_D1_COMPONENT_OUT_OF_SCOPE_MALFORMED_GATE0V_RECORD"
                    comp_g0v_standalone_status = "HOLD_D1_COMPONENT_OUT_OF_SCOPE_MALFORMED_GATE0V_RECORD"
                continue
            if _inspect_adverse_gate0v_failure(_comp_g0v_rec) is not None:
                comp_reasons.append("FAIL_OUT_OF_SCOPE_COMPONENT_GLOBAL_GATE0V_VERDICT_FAIL")
                if comp_out_of_scope_g0v_status is None:
                    comp_out_of_scope_g0v_status = "FAIL_OUT_OF_SCOPE_COMPONENT_GLOBAL_GATE0V_VERDICT_FAIL"
                    comp_g0v_standalone_status = "FAIL_OUT_OF_SCOPE_COMPONENT_GLOBAL_GATE0V_VERDICT_FAIL"
            _cn_lot_vals: List[str] = []
            _cn_lot_bad = False
            for _cnlk in _LOT_ALIAS_KEYS:
                if _cnlk not in _comp_g0v_rec:
                    continue
                _cnlv = _comp_g0v_rec[_cnlk]
                if _cnlv is None or not _is_canonical_identifier(_cnlv):
                    _cn_lot_bad = True
                elif _cnlv not in _cn_lot_vals:
                    _cn_lot_vals.append(_cnlv)
            if _cn_lot_bad:
                comp_reasons.append("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID")
                if comp_out_of_scope_g0v_status is None:
                    comp_out_of_scope_g0v_status = "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
                    comp_g0v_standalone_status = "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
            if len(_cn_lot_vals) > 1 or (lot_id_eff is not None and any(_v != lot_id_eff for _v in _cn_lot_vals)):
                comp_reasons.append("HOLD_D1_CONFLICTING_LOT_ID_BINDING")
                if comp_out_of_scope_g0v_status is None:
                    comp_out_of_scope_g0v_status = "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                    comp_g0v_standalone_status = "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
        if invalid_or_oob_keys:
            comp_reasons.append("FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION")
        if explicitly_unavailable:
            comp_reasons.append("HOLD_D1_OBSERVATION_EXPLICITLY_UNAVAILABLE")
        all_c = list(dict.fromkeys(reasons + comp_reasons))
        non_excl_reasons = [r for r in all_c if r != "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"]
        comp_is_canonical_eval = bool(
            require_canonical_sha
            and contract["staged_registry_verified"]
            and contract["canonical_sha_verified"]
            and contract["playbook_canonical_sha_verified"]
        )
        comp_raw_subgates_passed = len(non_excl_reasons) == 0
        comp_subgates_passed = bool(comp_raw_subgates_passed and comp_is_canonical_eval)
        if not comp_is_canonical_eval:
            all_c.append("HOLD_D1_NON_CANONICAL_DIAGNOSTIC_CONTRACT_UNPINNED")
        ord_c = (
            [r for r in all_c if r.startswith("FAIL_HARD_STOP_")]
            + [r for r in all_c if r.startswith("FAIL_") and not r.startswith("FAIL_HARD_STOP_")]
            + [r for r in all_c if not r.startswith("FAIL_") and r != "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"]
            + [r for r in all_c if r == "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"]
        )
        comp_cal_summary: Dict[str, Any] = {}
        comp_sampling_prov: Dict[str, Any] = {}
        if pvid == "ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9":
            comp_cal_summary["d1_04_ap1903"] = _cal_meta_from_cp(comp_ap_cp, comp_samp_meta.get("d1_04_ap1903"))
            comp_sampling_prov["d1_04_ap1903"] = (comp_samp_meta.get("d1_04_ap1903") or {}).get("sampling_independence_provenance")
        elif pvid == "ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1":
            comp_cal_summary["d1_08_hlac0304"] = _cal_meta_from_cp(comp_hlac_cp, comp_samp_meta.get("d1_08_hlac0304"))
            comp_sampling_prov["d1_08_hlac0304"] = (comp_samp_meta.get("d1_08_hlac0304") or {}).get("sampling_independence_provenance")
        comp_cov_regime = _aggregate_panel_coverage_regime(comp_cal_summary)
        comp_assumption_flags = _collect_assumption_only_flags(comp_cal_summary)
        comp_top_status = (
            "NON_CANONICAL_DIAGNOSTIC_COMPONENT_SUBGATES_MET_UNPINNED_EXCLUDED_FROM_FULLSTACK"
            if (not comp_is_canonical_eval and comp_raw_subgates_passed)
            else ord_c[0]
        )
        comp_emp_status = (
            "HOLD_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE"
            if comp_is_canonical_eval
            else "NON_CANONICAL_DIAGNOSTIC_COMPONENT_UNPINNED_NO_RELEASE"
        )
        return {
            "gate_id": "STAGE0_GATE_AUTOLOGOUS_D1B_AUTO_C_V2_D101_D109_IIP_0035",
            "decision_id": "DEC-114",
            "contribution_id": "IIP-0035",
            "controlling_contract_block_id": contract["gate_id"],
            "controlling_contract_schema_version": contract["schema_version"],
            "controlling_contract_sha256": contract["sha256"],
            "controlling_staged_registry_verified": contract["staged_registry_verified"],
            "controlling_contract_canonical_sha_verified": contract["canonical_sha_verified"],
            "controlling_contract_on_disk_verified": contract["on_disk_decision_tree_verified"],
            "controlling_contract_in_memory_matches_canonical_sha": contract["in_memory_block_matches_canonical_sha"],
            "controlling_playbook_sec4120_sha256": contract["playbook_sec4120_sha256"],
            "controlling_playbook_footer_verified": contract["playbook_footer_verified"],
            "controlling_playbook_canonical_sha_verified": contract["playbook_canonical_sha_verified"],
            "controlling_contract_cells_loaded": contract["cells_loaded"],
            "controlling_contract_external_refs": contract["external_refs"],
            "evaluation_mode": "CANONICAL_PINNED" if comp_is_canonical_eval else "NON_CANONICAL_DIAGNOSTIC_UNPINNED",
            "evaluation_kind": "COMPONENT_SUBGATE_REFERENCE_COMPUTATION_EXCLUDED_FROM_FULLSTACK_CLEARANCE",
            "reference_scenario_scope": "SINGLE_LOCUS_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK",
            "assembled_multilocus_manifest_bound": False,
            "clinical_or_batch_release_authorized": False,
            "canonical_product_variant_id": pvid,
            "effective_locus1_donor_arm": None,
            "evaluated_lot_id": lot_id_eff if _is_canonical_identifier(lot_id_eff) else None,
            "d1_08_allele_branch": comp_d1_08_branch,
            "d1_08_non_c0304_screen_attestation_scope": (
                "EXTERNALLY_ADJUDICATED_CALLER_ATTESTATION_PENDING_NORMATIVE_PROTOCOL"
                if comp_d1_08_branch == "L2_COMPONENT_NON_HLAC0304_MATCHED_ALLELE_CD8_SUBGATE"
                else None
            ),
            "global_gate0v_standalone_subcheck_status": comp_g0v_standalone_status,
            "global_gate0v_out_of_scope_adverse_record_status": comp_out_of_scope_g0v_status,
            "global_gate0v_lims_ledger_verified": None,
            "panel_coverage_scope": "PER_ASSAY_INDIVIDUAL_COVERAGE_NOT_SIMULTANEOUS_PANEL_FAMILYWISE_GUARANTEE",
            "coverage_regime": comp_cov_regime,
            "rare_event_calibration_summary": comp_cal_summary,
            "sampling_design_provenance_by_assay": comp_sampling_prov,
            "assumption_only_path_flags": comp_assumption_flags,
            "component_subgate_evaluation": {
                "applicable_gates": comp_gates,
                "evaluation_mode": "CANONICAL_PINNED" if comp_is_canonical_eval else "NON_CANONICAL_DIAGNOSTIC_UNPINNED",
                "component_subgates_passed": comp_subgates_passed,
                "component_subgates_diagnostic_passed": comp_raw_subgates_passed,
                "component_reasons": comp_reasons,
            },
            "d1_04_ap1903_u95_str": comp_ap_cp.get("u95_exact_sound_decimal_str") if comp_ap_cp else None,
            "d1_08_hlac0304_u95_pct_str": comp_hlac_cp.get("u95_exact_sound_decimal_pct_str") if comp_hlac_cp else None,
            "d1_08_t1d_serum_c5b9_u95_pct_str": comp_t1d_cp.get("u95_exact_sound_decimal_pct_str") if comp_t1d_cp else None,
            "empirical_clearance_status": comp_emp_status,
            "reference_criteria_met": False,
            "status": comp_top_status,
            "reasons": ord_c,
        }

    # Full-panel / unresolved-product evaluation path
    g0v_kw = {
        "day14_min": thr.get("g0v_day14_min", Decimal("3.0")),
        "day14_l95_min": thr.get("g0v_day14_l95_min", Decimal("3.0")),
        "early_72h_max": thr.get("g0v_72h_max", Decimal("2.0")),
        "day28_u95_max": thr.get("g0v_day28_u95_max", Decimal("0.0005")),
        "day28_n_min": int(thr.get("g0v_day28_n_min", 5990)),
        "panel_observations": observations,
    }
    _hs_g0v_l1a_or_unres = (
        "FAIL_HARD_STOP_D1_UNRESOLVED_PRODUCT_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"
        if (has_pvid_conflict or not is_fullstack)
        else "FAIL_HARD_STOP_D1_L1A_GLOBAL_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT"
    )
    _fail_g0v_l1a_or_unres = (
        "FAIL_D1_UNRESOLVED_PRODUCT_GLOBAL_GATE0V_VIRAL_SANCTUARY_FAILED"
        if (has_pvid_conflict or not is_fullstack)
        else "FAIL_D1_L1A_GLOBAL_GATE0V_VIRAL_SANCTUARY_FAILED"
    )
    _g0v_adverse_fail_reasons = sorted({
        r_f
        for rec in g0v_supplied_records
        if (r_f := _inspect_adverse_gate0v_failure(
            rec,
            day14_min=g0v_kw["day14_min"],
            day14_l95_min=g0v_kw["day14_l95_min"],
            early_72h_max=g0v_kw["early_72h_max"],
        )) is not None
    })
    if is_fullstack and not has_pvid_conflict and effective_l1_arm == "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR":
        if lot_conflict or lot_invalid:
            _lot_block_reason = "HOLD_D1_CONFLICTING_LOT_ID_BINDING" if lot_conflict else "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
            if g0v_any_positive_k:
                reasons.append("FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN")
                g0v = {
                    "cleared": False,
                    "failed": True,
                    "hard_stop": True,
                    "lot_binding_invalid": True,
                    "reason": "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
                    "lot_hold_reason": _lot_block_reason,
                    "adverse_event_pending_lot_attribution": True,
                }
            elif _g0v_adverse_fail_reasons:
                # Preserve Gate-0V viral-clearance failure & switch-to-L1A disposition symmetrically with nested lot conflicts!
                for _af_r in _g0v_adverse_fail_reasons:
                    reasons.append(_af_r)
                g0v = {
                    "cleared": False,
                    "failed": True,
                    "hard_stop": False,
                    "lot_binding_invalid": True,
                    "reason": _g0v_adverse_fail_reasons[0],
                    "lot_hold_reason": _lot_block_reason,
                }
            elif g0v_conflict:
                g0v = {
                    "cleared": False,
                    "failed": True,
                    "hard_stop": False,
                    "lot_binding_invalid": True,
                    "reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                }
            else:
                g0v = {
                    "cleared": False,
                    "failed": False,
                    "hard_stop": False,
                    "lot_binding_invalid": True,
                    "reason": _lot_block_reason,
                }
        else:
            g0v_evals = [
                _evaluate_structured_lot_matched_gate0v_record(rec, evaluated_lot_id=lot_id_eff, **g0v_kw)
                for rec in (g0v_supplied_records if g0v_supplied_records else [None])
            ]
            if not g0v_conflict:
                g0v = g0v_evals[0]
                if not g0v["cleared"]:
                    reasons.append(g0v["reason"])
                    if g0v.get("lot_hold_reason"):
                        reasons.append(g0v["lot_hold_reason"])
                    if g0v["reason"] == "HOLD_D1_L1B_STRUCTURED_LOT_MATCHED_GLOBAL_GATE0V_RECORD_REQUIRED" and "D1-09" in cell_mdd:
                        reasons.append(cell_mdd["D1-09"])
            else:
                _per_rec_reasons = sorted({
                    str(ev["reason"]) for ev in g0v_evals if not ev.get("cleared") and ev.get("reason")
                })
                _per_rec_lot_holds = sorted({
                    str(ev["lot_hold_reason"]) for ev in g0v_evals if ev.get("lot_hold_reason")
                })
                for _rr in _per_rec_reasons:
                    reasons.append(_rr)
                for _lh in _per_rec_lot_holds:
                    reasons.append(_lh)
                if g0v_any_positive_k or any(ev.get("hard_stop") for ev in g0v_evals):
                    reasons.append("FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN")
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": True,
                        "reason": "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
                        "conflict_reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                    }
                else:
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": False,
                        "reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                    }
    else:
        g0v = {"cleared": None, "failed": False, "hard_stop": False, "reason": "NOT_APPLICABLE_L1A_NO_CD46"}
        if g0v_conflict:
            if lot_conflict or lot_invalid:
                _lot_blk = "HOLD_D1_CONFLICTING_LOT_ID_BINDING" if lot_conflict else "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
                if g0v_any_positive_k:
                    reasons.append(_hs_g0v_l1a_or_unres)
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": True,
                        "reason": _hs_g0v_l1a_or_unres,
                        "conflict_reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                        "lot_hold_reason": _lot_blk,
                    }
                elif _g0v_adverse_fail_reasons:
                    reasons.append(_fail_g0v_l1a_or_unres)
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": False,
                        "lot_binding_invalid": True,
                        "reason": _fail_g0v_l1a_or_unres,
                        "conflict_reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                        "lot_hold_reason": _lot_blk,
                    }
                else:
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": False,
                        "reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                    }
            else:
                _l1a_evals = [
                    _evaluate_structured_lot_matched_gate0v_record(rec, evaluated_lot_id=lot_id_eff, **g0v_kw)
                    for rec in g0v_supplied_records
                ]
                _l1a_mapped: set[str] = set()
                _l1a_hs = g0v_any_positive_k
                for _ev in _l1a_evals:
                    if _ev.get("hard_stop"):
                        _l1a_hs = True
                        _l1a_mapped.add(_hs_g0v_l1a_or_unres)
                        if _ev.get("lot_hold_reason"):
                            _lhr_hs0 = str(_ev["lot_hold_reason"])
                            _l1a_mapped.add(
                                "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                                if _lhr_hs0 == "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"
                                else ("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" if _lhr_hs0 == "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID" else _lhr_hs0)
                            )
                    elif _ev.get("failed"):
                        _l1a_mapped.add(_fail_g0v_l1a_or_unres)
                        if _ev.get("lot_hold_reason"):
                            _lhr0 = str(_ev["lot_hold_reason"])
                            _l1a_mapped.add(
                                "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                                if _lhr0 == "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"
                                else ("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" if _lhr0 == "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID" else _lhr0)
                            )
                        if _ev.get("reason") == "FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE":
                            _l1a_mapped.add("FAIL_D1_GLOBAL_GATE0V_INVALID_OR_OUT_OF_RANGE")
                    elif _ev.get("lot_binding_invalid"):
                        _lh_r = str(_ev.get("lot_hold_reason") or _ev.get("reason") or "HOLD_D1_CONFLICTING_LOT_ID_BINDING")
                        _l1a_mapped.add(
                            "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                            if _lh_r == "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"
                            else ("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" if _lh_r == "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID" else _lh_r)
                        )
                for _mr in sorted(_l1a_mapped):
                    reasons.append(_mr)
                if _l1a_hs:
                    reasons.append(_hs_g0v_l1a_or_unres)
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": True,
                        "reason": _hs_g0v_l1a_or_unres,
                        "conflict_reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                    }
                else:
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": False,
                        "reason": "FAIL_D1_CONFLICTING_GLOBAL_GATE0V_DECLARATIONS",
                    }
        elif eff_g0v_record is not None:
            if lot_conflict or lot_invalid:
                _lot_blk = "HOLD_D1_CONFLICTING_LOT_ID_BINDING" if lot_conflict else "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID"
                if g0v_any_positive_k:
                    reasons.append(_hs_g0v_l1a_or_unres)
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": True,
                        "reason": _hs_g0v_l1a_or_unres,
                        "lot_hold_reason": _lot_blk,
                        "adverse_event_pending_lot_attribution": True,
                    }
                elif _g0v_adverse_fail_reasons:
                    reasons.append(_fail_g0v_l1a_or_unres)
                    g0v = {
                        "cleared": False,
                        "failed": True,
                        "hard_stop": False,
                        "lot_binding_invalid": True,
                        "reason": _fail_g0v_l1a_or_unres,
                        "lot_hold_reason": _lot_blk,
                    }
                else:
                    g0v = {"cleared": False, "failed": False, "hard_stop": False, "lot_binding_invalid": True, "reason": _lot_blk}
            else:
                g0v_chk = _evaluate_structured_lot_matched_gate0v_record(eff_g0v_record, evaluated_lot_id=lot_id_eff, **g0v_kw)
                if g0v_chk.get("hard_stop"):
                    reasons.append(_hs_g0v_l1a_or_unres)
                    if g0v_chk.get("lot_hold_reason"):
                        _lhr_hs1 = str(g0v_chk["lot_hold_reason"])
                        reasons.append(
                            "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                            if _lhr_hs1 == "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"
                            else ("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" if _lhr_hs1 == "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID" else _lhr_hs1)
                        )
                    g0v = {
                        **g0v_chk,
                        "cleared": False,
                        "reason": _hs_g0v_l1a_or_unres,
                    }
                elif g0v_chk.get("failed"):
                    reasons.append(_fail_g0v_l1a_or_unres)
                    if g0v_chk.get("lot_hold_reason"):
                        _lhr1 = str(g0v_chk["lot_hold_reason"])
                        reasons.append(
                            "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                            if _lhr1 == "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"
                            else ("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" if _lhr1 == "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID" else _lhr1)
                        )
                    g0v = {
                        **g0v_chk,
                        "cleared": False,
                        "reason": _fail_g0v_l1a_or_unres,
                    }
                elif g0v_chk.get("lot_binding_invalid"):
                    # Any malformed or conflicting nested lot alias in an optional L1A/unresolved Gate-0V record MUST hold attribution!
                    _lhr2 = str(g0v_chk.get("lot_hold_reason") or g0v_chk.get("reason") or "HOLD_D1_CONFLICTING_LOT_ID_BINDING")
                    _mapped_lh = (
                        "HOLD_D1_CONFLICTING_LOT_ID_BINDING"
                        if _lhr2 == "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH"
                        else ("HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID" if _lhr2 == "HOLD_D1_L1B_GATE0V_MISSING_OR_NONCANONICAL_MCB_LOT_ID" else _lhr2)
                    )
                    reasons.append(_mapped_lh)
                    g0v = {
                        **g0v_chk,
                        "cleared": False,
                        "reason": _mapped_lh,
                    }
                elif g0v_chk.get("cleared") is True:
                    g0v = {
                        **g0v_chk,
                        "reason": "OPTIONAL_GATE0V_RECORD_PASSED_NOT_REQUIRED_FOR_L1A_LIVE_LIMS_UNVERIFIED",
                    }
                else:
                    g0v = {
                        **g0v_chk,
                        "reason": f"OPTIONAL_GATE0V_RECORD_NONPASSING_NOT_REQUIRED_FOR_L1A:{g0v_chk.get('reason')}",
                    }

    # D1-01 (supports direct derivation from raw 2-field patient_hla_c_alleles against Lin 2023 PMC10690437 Table 1 allowlist)
    geno_ok = observations.get("d1_01_4digit_hla_abc_kir_genotyping_verified")
    hlac_cls_present = "d1_auto_hlac_sp3c4c_genotype_class" in observations
    hlac_cls = observations.get("d1_auto_hlac_sp3c4c_genotype_class")
    d1_01_allele_derivation_mode = "CALLER_ATTESTED_GENOTYPE_CLASS_RAW_ALLELES_NOT_SUPPLIED"
    if d1_01_raw_hla_c_input_status == "MALFORMED":
        d1_01_allele_derivation_mode = "MALFORMED_RAW_HLA_C_ALLELES_REJECTED"
    elif d1_01_raw_hla_c_input_status == "EXPLICIT_UNAVAILABLE":
        d1_01_allele_derivation_mode = "RAW_HLA_C_ALLELES_EXPLICITLY_UNAVAILABLE_CALLER_ATTESTED"
    elif d1_01_raw_hla_c_input_status == "EXPLICIT_NULL":
        d1_01_allele_derivation_mode = "RAW_HLA_C_ALLELES_EXPLICIT_NULL_CALLER_ATTESTED"
    if isinstance(geno_ok, str) and geno_ok == "UNAVAILABLE":
        explicitly_unavailable.append("d1_01_4digit_hla_abc_kir_genotyping_verified")
        geno_ok = None
    elif geno_ok is not None and not isinstance(geno_ok, bool):
        invalid_or_oob_keys.append("d1_01_4digit_hla_abc_kir_genotyping_verified")
        geno_ok = None
    elif geno_ok is False:
        reasons.append("FAIL_D1_01_4DIGIT_HLA_ABC_KIR_GENOTYPING_MISSING")
    d1_01_exclusion_detected_from_raw = False
    if parsed_hlac_alleles is not None:
        norm_a = parsed_hlac_alleles
        d1_01_allele_derivation_mode = "DERIVED_FROM_SUBMITTED_4DIGIT_HLA_C_ALLELES"
        both_nonfunc = all(a in _SP3C4C_NONFUNC_HLAC_ALLELES for a in norm_a)
        any_confirmed_vl9 = any(a in _SP1C2C_VL9_COMPETENT_HLAC_ALLELES for a in norm_a)
        if both_nonfunc:
            d1_01_exclusion_detected_from_raw = True
            is_homo = (norm_a[0] == norm_a[1])
            derived_cls = (
                "EXCLUDED_SP3C4C_HLAC07_17_HOMOZYGOTE_ROUTE_TO_ALLO_Y_OR_D1A"
                if (is_homo and norm_a[0] != "HLA-C*18:01")
                else "EXCLUDED_SP3C4C_HLAC07_17_1801_NONFUNCTIONAL_HOMO_OR_COMPOUND_HET_ROUTE_TO_ALLO_Y_OR_D1A"
            )
            matching_excl = (
                {"EXCLUDED_SP3C4C_HLAC07_17_HOMOZYGOTE_ROUTE_TO_ALLO_Y_OR_D1A", "EXCLUDED_SP3C4C_HLAC07_17_1801_NONFUNCTIONAL_HOMO_OR_COMPOUND_HET_ROUTE_TO_ALLO_Y_OR_D1A"}
                if (is_homo and norm_a[0] != "HLA-C*18:01")
                else {"EXCLUDED_SP3C4C_HLAC07_17_1801_NONFUNCTIONAL_HOMO_OR_COMPOUND_HET_ROUTE_TO_ALLO_Y_OR_D1A"}
            )
            # Raw SP-3C/SP-4C nonfunctional pair ALWAYS enforces exclusion FAIL (even if hlac_cls is None or "UNAVAILABLE")
            reasons.append("FAIL_D1_01_EXCLUDED_AH8P1_SP3C4C_NONFUNCTIONAL_HLAC_ROUTE_TO_ALLO_Y_OR_D1A")
            if not hlac_cls_present:
                hlac_cls = derived_cls
            elif hlac_cls is None or (isinstance(hlac_cls, str) and hlac_cls == "UNAVAILABLE"):
                d1_01_allele_derivation_mode = "RAW_EXCLUSION_FAIL_WITH_EXPLICIT_NULL_OR_UNAVAILABLE_CLASS"
            elif isinstance(hlac_cls, str) and hlac_cls not in matching_excl:
                reasons.append("FAIL_D1_01_ALLELE_CLASS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE")
        elif any_confirmed_vl9:
            derived_cls = "ELIGIBLE_VL9_COMPETENT_HLAC"
            if not hlac_cls_present:
                hlac_cls = derived_cls
            elif isinstance(hlac_cls, str) and hlac_cls != "UNAVAILABLE" and hlac_cls != derived_cls:
                reasons.append("FAIL_D1_01_ALLELE_CLASS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE")
        else:
            d1_01_allele_derivation_mode = "UNCLASSIFIED_HLA_C_ALLELE_OUTSIDE_LIN2023_TABLE1_HOLD"
            reasons.append("HOLD_D1_01_UNCLASSIFIED_HLA_C_SIGNAL_PEPTIDE_ALLELE_REQUIRES_IPD_IMGT_VERIFICATION")
            if isinstance(hlac_cls, str) and hlac_cls in EXCLUDED_HLAC_SP3C4C_CLASSES:
                reasons.append("FAIL_D1_01_ALLELE_CLASS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE")
    d1_01_unknown_cls = False
    if isinstance(hlac_cls, str) and hlac_cls == "UNAVAILABLE":
        explicitly_unavailable.append("d1_auto_hlac_sp3c4c_genotype_class")
        hlac_cls = None
    elif hlac_cls is not None and not isinstance(hlac_cls, str):
        invalid_or_oob_keys.append("d1_auto_hlac_sp3c4c_genotype_class")
        hlac_cls = None
    elif hlac_cls is not None and hlac_cls not in ALLOWED_HLAC_SP3C4C_CLASSES:
        d1_01_unknown_cls = True
        reasons.append("HOLD_D1_01_UNKNOWN_HLAC_SP3C4C_GENOTYPE_CLASS")
    elif hlac_cls in EXCLUDED_HLAC_SP3C4C_CLASSES:
        reasons.append("FAIL_D1_01_EXCLUDED_AH8P1_SP3C4C_NONFUNCTIONAL_HLAC_ROUTE_TO_ALLO_Y_OR_D1A")
    if geno_ok is None or hlac_cls is None or d1_01_unknown_cls:
        reasons.append(cell_mdd.get("D1-01", "HOLD_D1_01_MISSING_GENOTYPING_OR_UNKNOWN_HLAC_CLASS"))
    if parsed_hlac_alleles is None and hlac_cls is None and d1_01_raw_hla_c_input_status != "MALFORMED":
        d1_01_allele_derivation_mode = "MISSING_OR_UNAVAILABLE_GENOTYPE_CLASS_HOLD"
    elif (
        parsed_hlac_alleles is not None
        and hlac_cls is None
        and not d1_01_exclusion_detected_from_raw
        and d1_01_allele_derivation_mode == "DERIVED_FROM_SUBMITTED_4DIGIT_HLA_C_ALLELES"
    ):
        d1_01_allele_derivation_mode = "EXPLICIT_NULL_OR_UNAVAILABLE_BLOCKED_DERIVATION_HOLD"

    # D1-02 (COMPOUND_ALL_PASS: requires BOTH caller-attested phased WGS == True AND submitted min mismatches in 3..25 plus all 3 numeric WGS readouts)
    wgs_ok = observations.get("d1_02_n1_phased_wgs_offgrid_guard_011_cleared")
    min_mm_hlac_e = _strict_int("obs_d1_02_min_mismatches_to_patient_hlac_e")
    _min_mm_max = int(thr.get("min_mm_max", 25))
    if min_mm_hlac_e is not None and min_mm_hlac_e > _min_mm_max:
        invalid_or_oob_keys.append("obs_d1_02_min_mismatches_to_patient_hlac_e")
        min_mm_hlac_e = None
    if isinstance(wgs_ok, str) and wgs_ok == "UNAVAILABLE":
        explicitly_unavailable.append("d1_02_n1_phased_wgs_offgrid_guard_011_cleared")
        wgs_ok = None
    elif wgs_ok is not None and not isinstance(wgs_ok, bool):
        invalid_or_oob_keys.append("d1_02_n1_phased_wgs_offgrid_guard_011_cleared")
        wgs_ok = None
    elif wgs_ok is False:
        reasons.append("FAIL_D1_02_N1_PHASED_WGS_OFFGRID_GUARD_011_HLAC_E_MISMATCH_CONFLICT")
    if min_mm_hlac_e is not None and min_mm_hlac_e <= thr.get("min_mm_reject_le", 2):
        reasons.append("FAIL_D1_02_N1_PHASED_WGS_OFFGRID_GUARD_011_HLAC_E_MISMATCH_CONFLICT")
    wt_hlac = _dec("obs_d1_02_wt_hlac_allele_retention_pct")
    bystander = _dec("obs_d1_02_hlac_bystander_ct_pct")
    chr6_del = _dec("obs_d1_02_chr6p21_interlocus_del_pct")
    if wt_hlac is not None and wt_hlac < thr.get("wt_hlac_min", Decimal("99.0")):
        reasons.append("FAIL_D1_02_WT_HLAC_ALLELE_RETENTION_BELOW_99PCT")
    if bystander is not None and bystander >= thr.get("bystander_lt", Decimal("0.5")):
        reasons.append("FAIL_D1_02_HLAC_BYSTANDER_CT_DEAMINATION_GE_0P5PCT")
    if chr6_del is not None and chr6_del > thr.get("chr6_del_eq", Decimal("0")):
        reasons.append("FAIL_D1_02_CHR6P21_INTERLOCUS_DELETION_DETECTED")
    wgs_ok_effective = (
        wgs_ok is True
        and min_mm_hlac_e is not None
        and min_mm_hlac_e > thr.get("min_mm_reject_le", 2)
    )
    d1_02_offgrid_guard_evidence_mode = (
        "CALLER_ATTESTED_PHASED_WGS_AND_SUBMITTED_MISMATCH_COUNT_COMPOUND_ALL_PASS"
        if wgs_ok_effective
        else (
            "FAILED_D1_02_PHASED_WGS_OR_MIN_MISMATCH"
            if (wgs_ok is False or (min_mm_hlac_e is not None and min_mm_hlac_e <= thr.get("min_mm_reject_le", 2)))
            else "INCOMPLETE_D1_02_PHASED_WGS_OR_MIN_MISMATCH_HOLD"
        )
    )
    if any(x is None for x in (wgs_ok, min_mm_hlac_e, wt_hlac, bystander, chr6_del)):
        reasons.append(cell_mdd.get("D1-02", "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC"))

    # D1-03 (requires head-to-head ADSC vs blood CiPSC comparator attestation matching Playbook §4.1.20)
    cipsc_src = observations.get("d1_auto_cipsc_source_arm")
    h2h_ok = observations.get("d1_03_head_to_head_adsc_vs_blood_comparator_attested")
    bcor_ok = observations.get("d1_03_adsc_vs_blood_bcor_tp53_ppm1d_wgs_clean")
    d1_03_unknown_src = False
    if isinstance(cipsc_src, str) and cipsc_src == "UNAVAILABLE":
        explicitly_unavailable.append("d1_auto_cipsc_source_arm")
        cipsc_src = None
    elif cipsc_src is not None and not isinstance(cipsc_src, str):
        invalid_or_oob_keys.append("d1_auto_cipsc_source_arm")
        cipsc_src = None
    elif cipsc_src is not None and cipsc_src not in ALLOWED_CIPSC_SOURCES:
        d1_03_unknown_src = True
        reasons.append("HOLD_D1_03_UNKNOWN_CIPSC_SOURCE_ARM")
    elif cipsc_src == "SUN_EXPOSED_DERMAL_FIBROBLAST_DISCOURAGED":
        reasons.append("FAIL_D1_03_SUN_EXPOSED_FIBROBLAST_UV_BURDEN_BARRED")
    if isinstance(h2h_ok, str) and h2h_ok == "UNAVAILABLE":
        explicitly_unavailable.append("d1_03_head_to_head_adsc_vs_blood_comparator_attested")
        h2h_ok = None
    elif h2h_ok is not None and not isinstance(h2h_ok, bool):
        invalid_or_oob_keys.append("d1_03_head_to_head_adsc_vs_blood_comparator_attested")
        h2h_ok = None
    elif h2h_ok is False:
        reasons.append("FAIL_D1_03_HEAD_TO_HEAD_ADSC_VS_BLOOD_COMPARATOR_NOT_MET")
    if h2h_ok is None:
        reasons.append("HOLD_D1_03_HEAD_TO_HEAD_ADSC_VS_BLOOD_COMPARATOR_REQUIRED")
    if isinstance(bcor_ok, str) and bcor_ok == "UNAVAILABLE":
        explicitly_unavailable.append("d1_03_adsc_vs_blood_bcor_tp53_ppm1d_wgs_clean")
        bcor_ok = None
    elif bcor_ok is not None and not isinstance(bcor_ok, bool):
        invalid_or_oob_keys.append("d1_03_adsc_vs_blood_bcor_tp53_ppm1d_wgs_clean")
        bcor_ok = None
    elif bcor_ok is False:
        reasons.append("FAIL_D1_03_BCOR_TP53_PPM1D_OR_20Q11_1Q32_CLONAL_HIT_DETECTED")
    chrm_depth = _dec("obs_d1_03_chrm_numt_amplicon_depth_x")
    chrm_het = _dec("obs_d1_03_chrm_novel_nonsyn_het_pct")
    if chrm_depth is not None and chrm_depth < thr.get("chrm_depth_min", Decimal("1000.0")):
        reasons.append("FAIL_D1_03_CHRM_NUMT_AMPLICON_DEPTH_BELOW_1000X")
    if chrm_het is not None and chrm_het >= thr.get("chrm_het_lt", Decimal("1.0")):
        reasons.append("FAIL_D1_03_CHRM_NOVEL_NONSYN_HETEROPLASMY_GE_1PCT")
    if any(x is None for x in (cipsc_src, h2h_ok, bcor_ok, chrm_depth, chrm_het)) or d1_03_unknown_src:
        reasons.append(cell_mdd.get("D1-03", "HOLD_D1_03_INCOMPLETE_WGS_OR_CHRM_QC"))

    # D1-04
    chr13_ok = observations.get("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified")
    if isinstance(chr13_ok, str) and chr13_ok == "UNAVAILABLE":
        explicitly_unavailable.append("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified")
        chr13_ok = None
    elif chr13_ok is not None and not isinstance(chr13_ok, bool):
        invalid_or_oob_keys.append("d1_04_chr13_clybl_tandem_404bp_independent_from_chr19_verified")
        chr13_ok = None
    elif chr13_ok is False:
        reasons.append("FAIL_D1_04_CHR13_CLYBL_TANDEM_INDEPENDENCE_FROM_CHR19_UNVERIFIED")
    azt_chr19_loss_ok = observations.get("d1_04_azt_kill_measured_after_forced_chr19_loss")
    if isinstance(azt_chr19_loss_ok, str) and azt_chr19_loss_ok == "UNAVAILABLE":
        explicitly_unavailable.append("d1_04_azt_kill_measured_after_forced_chr19_loss")
        azt_chr19_loss_ok = None
    elif azt_chr19_loss_ok is not None and not isinstance(azt_chr19_loss_ok, bool):
        invalid_or_oob_keys.append("d1_04_azt_kill_measured_after_forced_chr19_loss")
        azt_chr19_loss_ok = None
    elif azt_chr19_loss_ok is False:
        reasons.append("FAIL_D1_04_AZT_KILL_NOT_MEASURED_AFTER_FORCED_CHR19_LOSS")
    if azt_chr19_loss_ok is None:
        reasons.append("HOLD_D1_04_AZT_KILL_FORCED_CHR19_LOSS_CONDITION_REQUIRED")
    htmpk_leak = _dec("obs_d1_04_g0_beta_htmpk_leak_pct")
    elispot = _dec("obs_d1_04_hlac07_ryafsgvay_elispot_above_bg_pct")
    ap_surv = _dec("obs_d1_04_dec053_ap1903_icasp9_survivor_frac")
    ap_k = _strict_int("obs_d1_04_ap1903_survivor_count_k")
    ap_n_raw_key = _strict_int("obs_d1_04_ap1903_raw_trials_n")
    ap_n_eff_key = _strict_int("obs_d1_04_ap1903_effective_trials_n")
    if ap_n_raw_key is not None and ap_n_eff_key is not None and ap_n_raw_key != ap_n_eff_key:
        reasons.append("FAIL_D1_04_CONFLICTING_RAW_AND_EFFECTIVE_TRIALS_N")
    ap_n = ap_n_raw_key if ap_n_raw_key is not None else ap_n_eff_key
    azt_kill = _dec("obs_d1_04_dec053_azt_dividing_log10_kill")
    ap_u95_cap = thr.get("ap1903_u95_max", Decimal("1e-8"))
    if htmpk_leak is not None and htmpk_leak > thr.get("htmpk_leak_max", Decimal("0.1")):
        reasons.append("FAIL_D1_04_G0_BETA_HTMPK_F105Y_LEAK_EXCEEDS_0P1PCT")
    if elispot is not None and elispot > thr.get("elispot_max", Decimal("0.1")):
        reasons.append("FAIL_D1_04_HLAC07_RYAFSGVAY_ELISPOT_EXCEEDS_0P1PCT_G0_GATE")
    if ap_surv is not None and ap_surv > ap_u95_cap:
        reasons.append("FAIL_D1_04_DEC053_AP1903_ICASP9_SURVIVOR_FRAC_EXCEEDS_1E_8")

    _ap_cal_kw = _cal_kwargs("obs_d1_04_ap1903")
    ap_cp: Dict[str, Any] | None = None
    if ap_k is None or ap_n is None or ap_n < 1 or ap_k > ap_n:
        if ap_k is not None and ap_k >= 1:
            reasons.append("FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED")
        else:
            reasons.append("HOLD_D1_04_AP1903_GY018_INTEGER_KN_REQUIRED")
    else:
        ap_cp = evaluate_sound_decimal_clopper_pearson_zero_event(
            gate_id="D1_04_AP1903_ICASP9",
            k=ap_k,
            n_raw=ap_n,
            max_u95_fraction=ap_u95_cap,
            **_ap_cal_kw,
        )
        if ap_cp["status"] == "FAIL_HARD_STOP":
            reasons.append("FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED")
        elif ap_cp["status"] != "PASS":
            reasons.append(
                "HOLD_D1_04_AP1903_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                if ap_cp["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                else f"HOLD_D1_04_AP1903_{ap_cp['status']}"
            )
        if ap_surv is not None:
            with _clean_dec_ctx(40):
                ap_expected = Decimal(ap_k) / Decimal(ap_n)
                ap_tol = max(Decimal("0.05") * ap_expected, Decimal("1e-15"))
                if (
                    (ap_k == 0 and ap_surv != Decimal(0))
                    or (ap_k > 0 and ap_surv <= Decimal(0))
                    or (ap_k > 0 and abs(ap_surv - ap_expected) > ap_tol)
                ):
                    reasons.append("FAIL_D1_04_AP1903_SURVIVOR_FRAC_INCONSISTENT_WITH_KN")
    if azt_kill is not None and azt_kill < thr.get("azt_kill_min", Decimal("4.0")):
        reasons.append("FAIL_D1_04_DEC053_AZT_DIVIDING_KILL_BELOW_4LOG_AFTER_CHR19_LOSS")
    if any(x is None for x in (chr13_ok, azt_chr19_loss_ok, htmpk_leak, elispot, azt_kill, ap_k, ap_n)):
        reasons.append(cell_mdd.get("D1-04", "HOLD_D1_04_MISSING_SWITCH_OR_GY018_INTEGER_KN"))

    # D1-05 (MANDATORY per-locus L3 AND L4 Stage-7 Day-28 retention; context-isolated exact min comparison!)
    l3_ret = _dec("obs_d1_05_l3_day28_retention_pct")
    l4_ret = _dec("obs_d1_05_l4_day28_retention_pct")
    day28_ret = _dec("obs_d1_05_l3_l4_day28_retention_pct")
    fsp1_surv = _dec("obs_d1_05_fsp1_ferroptosis_survival_pct")
    if l3_ret is not None and l3_ret < thr.get("l3_ret_min", Decimal("80.0")):
        reasons.append("FAIL_D1_05_L3_STAGE7_DAY28_RETENTION_BELOW_80PCT")
    if l4_ret is not None and l4_ret < thr.get("l4_ret_min", Decimal("80.0")):
        reasons.append("FAIL_D1_05_L4_STAGE7_DAY28_RETENTION_BELOW_80PCT")
    if l3_ret is None or l4_ret is None or fsp1_surv is None:
        reasons.append(cell_mdd.get("D1-05", "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1"))
    if l3_ret is not None and l4_ret is not None:
        actual_min_l3_l4 = min(l3_ret, l4_ret)
        if actual_min_l3_l4 < thr.get("min_l3_l4_ret_min", Decimal("80.0")):
            reasons.append("FAIL_D1_05_L3_L4_STAGE7_DAY28_RETENTION_BELOW_80PCT")
        if day28_ret is not None:
            _min_diff_exceeds = _abs_diff_exceeds_tol(day28_ret, actual_min_l3_l4, Decimal("1e-9"))
            if _min_diff_exceeds is True:
                reasons.append("FAIL_D1_05_INCONSISTENT_REPORTED_MIN_L3_L4_RETENTION")
            elif _min_diff_exceeds is None:
                reasons.append("HOLD_D1_05_REPORTED_MIN_L3_L4_COMPARISON_INCONCLUSIVE")
    if day28_ret is not None and day28_ret < thr.get("min_l3_l4_ret_min", Decimal("80.0")):
        reasons.append("FAIL_D1_05_L3_L4_STAGE7_DAY28_RETENTION_BELOW_80PCT")
    if fsp1_surv is not None and fsp1_surv < thr.get("fsp1_surv_min", Decimal("80.0")):
        reasons.append("FAIL_D1_05_FSP1_FERROPTOSIS_SURVIVAL_BELOW_80PCT")

    # D1-06 & D1-07
    kir_nk = _dec("obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct")
    vtcn1_shed = _dec("obs_d1_07_vtcn1_nrd1_sb7h4_shedding_pct")
    if kir_nk is None:
        reasons.append(cell_mdd.get("D1-06", "HOLD_D1_06_MISSING_KIR3DL1_NK_ASSAY"))
    elif kir_nk > thr.get("kir_nk_max", Decimal("8.0")):
        reasons.append("FAIL_D1_06_KIR3DL1_MISSING_BW4_NK_LYSIS_EXCEEDS_8PCT")
    if vtcn1_shed is None:
        reasons.append(cell_mdd.get("D1-07", "HOLD_D1_07_MISSING_VTCN1_SHEDDING_ASSAY"))
    elif vtcn1_shed > thr.get("vtcn1_shed_max", Decimal("10.0")):
        reasons.append("FAIL_D1_07_VTCN1_NRD1_SB7H4_SHEDDING_EXCEEDS_10PCT_MUTAGENESIS_REQUIRED")

    # D1-08 (Two-branch allele-conditional residual CD8+ T-cell screen & intact T1D serum C5b-9)
    hlac03_u95_pct_cap = thr.get("hlac03_u95_pct_max", Decimal("5.0"))
    t1d_c5b9_u95_pct_cap = thr.get("t1d_c5b9_u95_pct_max", Decimal("5.0"))

    carries_c0304_present = "d1_08_patient_carries_hlac0304" in observations
    carries_c0304 = observations.get("d1_08_patient_carries_hlac0304")
    carries_c0304_raw_explicit_unavail_or_invalid = False
    if carries_c0304_present and carries_c0304 is None:
        carries_c0304_raw_explicit_unavail_or_invalid = True
    elif isinstance(carries_c0304, str) and carries_c0304 == "UNAVAILABLE":
        explicitly_unavailable.append("d1_08_patient_carries_hlac0304")
        carries_c0304 = None
        carries_c0304_raw_explicit_unavail_or_invalid = True
    elif carries_c0304 is not None and not isinstance(carries_c0304, bool):
        invalid_or_oob_keys.append("d1_08_patient_carries_hlac0304")
        carries_c0304 = None
        carries_c0304_raw_explicit_unavail_or_invalid = True
    if derived_carries_hlac0304 is not None:
        if isinstance(carries_c0304, bool) and carries_c0304 != derived_carries_hlac0304:
            reasons.append("FAIL_D1_08_HLAC0304_CARRIER_STATUS_MISMATCH_WITH_RAW_HLA_C_GENOTYPE")
            if derived_carries_hlac0304 is True:
                carries_c0304 = True
        elif carries_c0304 is None and not carries_c0304_raw_explicit_unavail_or_invalid:
            carries_c0304 = derived_carries_hlac0304
    elif unlisted_hlac_blocks_d1_08_negative and carries_c0304 is None and not carries_c0304_raw_explicit_unavail_or_invalid:
        reasons.append("HOLD_D1_08_UNLISTED_RAW_HLAC_ALLELE_BLOCKS_NEGATIVE_DERIVATION")

    hlac03_t = _dec("obs_d1_08_hlac0304_tcell_lysis_u95_pct")
    hlac03_k = _strict_int("obs_d1_08_hlac0304_lysed_events_k")
    hlac03_n = _strict_int("obs_d1_08_hlac0304_scored_targets_n")
    t1d_c5b9 = _dec("obs_d1_08_t1d_serum_c5b9_u95_pct")
    t1d_k = _strict_int("obs_d1_08_t1d_serum_c5b9_pos_events_k")
    t1d_n = _strict_int("obs_d1_08_t1d_serum_c5b9_scored_cells_n")

    _hlac03_cal_kw = _cal_kwargs("obs_d1_08_hlac0304")
    _t1d_cal_kw = _cal_kwargs("obs_d1_08_t1d_serum_c5b9")

    hlac03_cp: Dict[str, Any] | None = None
    d1_08_branch_ok = False
    c0304_pc: Any = None
    ma_pc: Any = None
    ma_ok: Any = None
    if carries_c0304 is True:
        d1_08_allele_branch = "HLAC0304_CARRIER_INSULIN_ACHAIN_AND_C5B9"
        c0304_pc = observations.get("d1_08_hlac0304_positive_control_valid")
        if isinstance(c0304_pc, str) and c0304_pc == "UNAVAILABLE":
            explicitly_unavailable.append("d1_08_hlac0304_positive_control_valid")
            c0304_pc = None
        elif c0304_pc is not None and not isinstance(c0304_pc, bool):
            invalid_or_oob_keys.append("d1_08_hlac0304_positive_control_valid")
            c0304_pc = None
        elif c0304_pc is False:
            reasons.append("FAIL_D1_08_HLAC0304_POSITIVE_CONTROL_INVALID")
        if c0304_pc is None:
            reasons.append("HOLD_D1_08_HLAC0304_POSITIVE_CONTROL_VALIDITY_REQUIRED")
        else:
            d1_08_branch_ok = bool(c0304_pc is True)

        if hlac03_k is None or hlac03_n is None or hlac03_n < 1 or hlac03_k > hlac03_n:
            if hlac03_k is not None and hlac03_k >= 1:
                reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
            else:
                reasons.append("HOLD_D1_08_HLAC0304_GY018_INTEGER_KN_REQUIRED")
                if hlac03_t is None:
                    reasons.append("HOLD_D1_08_HLAC0304_REPORTED_U95_REQUIRED")
        else:
            hlac03_cp = evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_08_HLAC0304_TCELL_LYSIS",
                k=hlac03_k,
                n_raw=hlac03_n,
                max_u95_fraction=_exact_pct_to_frac(hlac03_u95_pct_cap),
                reported_u95_pct=hlac03_t,
                **_hlac03_cal_kw,
            )
            if hlac03_cp["status"] == "FAIL_HARD_STOP":
                reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
                if hlac03_t is not None and hlac03_t <= Decimal(0):
                    reasons.append("FAIL_D1_08_INCONSISTENT_ZERO_UPPER_BOUND_WITH_EVENTS")
            elif hlac03_cp["status"] != "PASS":
                reasons.append(
                    "HOLD_D1_08_HLAC0304_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                    if hlac03_cp["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                    else f"HOLD_D1_08_HLAC0304_{hlac03_cp['status']}"
                )
                if (
                    hlac03_t is not None
                    and hlac03_k == 0
                    and hlac03_t <= hlac03_u95_pct_cap
                    and hlac03_cp.get("reported_u95_meets_exact_bound") is False
                ):
                    reasons.append("FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND")
            else:
                if hlac03_t is None:
                    reasons.append("HOLD_D1_08_HLAC0304_REPORTED_U95_REQUIRED")
                elif hlac03_t > hlac03_u95_pct_cap:
                    reasons.append("HOLD_D1_08_HLAC0304_REPORTED_U95_ABOVE_GATE_LIMIT")
                elif hlac03_t <= Decimal(0) or hlac03_cp.get("reported_u95_meets_exact_bound") is False:
                    reasons.append("FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND")
                elif hlac03_cp.get("reported_u95_meets_exact_bound") is None:
                    reasons.append("HOLD_D1_08_REPORTED_U95_EXACT_FLOOR_INCONCLUSIVE")
    elif carries_c0304 is False:
        d1_08_allele_branch = "NON_HLAC0304_MATCHED_ALLELE_SCREEN_AND_C5B9"
        ma_pc = observations.get("d1_08_matched_allele_positive_control_valid")
        ma_ok = observations.get("d1_08_non_c0304_matched_allele_screen_passed")
        if isinstance(ma_pc, str) and ma_pc == "UNAVAILABLE":
            explicitly_unavailable.append("d1_08_matched_allele_positive_control_valid")
            ma_pc = None
        elif ma_pc is not None and not isinstance(ma_pc, bool):
            invalid_or_oob_keys.append("d1_08_matched_allele_positive_control_valid")
            ma_pc = None
        elif ma_pc is False:
            reasons.append("FAIL_D1_08_MATCHED_ALLELE_POSITIVE_CONTROL_INVALID")
        if ma_pc is None:
            reasons.append("HOLD_D1_08_MATCHED_ALLELE_POSITIVE_CONTROL_VALIDITY_REQUIRED")

        if isinstance(ma_ok, str) and ma_ok == "UNAVAILABLE":
            explicitly_unavailable.append("d1_08_non_c0304_matched_allele_screen_passed")
            ma_ok = None
        elif ma_ok is not None and not isinstance(ma_ok, bool):
            invalid_or_oob_keys.append("d1_08_non_c0304_matched_allele_screen_passed")
            ma_ok = None
        elif ma_ok is False:
            reasons.append("FAIL_D1_08_MATCHED_ALLELE_CD8_SCREEN_FAILED")
        if ma_ok is None:
            reasons.append("HOLD_D1_08_MATCHED_ALLELE_CD8_SCREEN_REQUIRED")

        d1_08_branch_ok = bool(ma_pc is True and ma_ok is True)
        if hlac03_k is not None and hlac03_k >= 1:
            reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
        if _hlac03_cal_kw["sampling_independence_conflict"]:
            reasons.append("HOLD_D1_08_HLAC0304_HOLD_CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS")
        elif _hlac03_cal_kw["sampling_malformed"]:
            reasons.append("HOLD_D1_08_HLAC0304_HOLD_MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION")
        elif _hlac03_cal_kw["clustered_sampling"]:
            reasons.append("HOLD_D1_08_HLAC0304_HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED")
    else:
        d1_08_allele_branch = "UNSPECIFIED_HLAC0304_CARRIER_STATUS_HOLD"
        reasons.append("HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED")
        if hlac03_k is not None and hlac03_k >= 1:
            reasons.append("FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED")
        if _hlac03_cal_kw["sampling_independence_conflict"]:
            reasons.append("HOLD_D1_08_HLAC0304_HOLD_CONFLICTING_SAMPLING_INDEPENDENCE_DECLARATIONS")
        elif _hlac03_cal_kw["sampling_malformed"]:
            reasons.append("HOLD_D1_08_HLAC0304_HOLD_MALFORMED_SAMPLING_INDEPENDENCE_DECLARATION")
        elif _hlac03_cal_kw["clustered_sampling"]:
            reasons.append("HOLD_D1_08_HLAC0304_HOLD_NON_INDEPENDENT_OR_CLUSTERED_SAMPLING_REJECTED")

    t1d_cp: Dict[str, Any] | None = None
    if t1d_k is None or t1d_n is None or t1d_n < 1 or t1d_k > t1d_n:
        if t1d_k is not None and t1d_k >= 1:
            reasons.append("FAIL_HARD_STOP_D1_08_T1D_SERUM_C5B9_EVENT_DETECTED")
        else:
            reasons.append("HOLD_D1_08_T1D_SERUM_C5B9_GY018_INTEGER_KN_REQUIRED")
            if t1d_c5b9 is None:
                reasons.append("HOLD_D1_08_T1D_SERUM_C5B9_REPORTED_U95_REQUIRED")
    else:
        t1d_cp = evaluate_sound_decimal_clopper_pearson_zero_event(
            gate_id="D1_08_T1D_SERUM_C5B9",
            k=t1d_k,
            n_raw=t1d_n,
            max_u95_fraction=_exact_pct_to_frac(t1d_c5b9_u95_pct_cap),
            reported_u95_pct=t1d_c5b9,
            **_t1d_cal_kw,
        )
        if t1d_cp["status"] == "FAIL_HARD_STOP":
            reasons.append("FAIL_HARD_STOP_D1_08_T1D_SERUM_C5B9_EVENT_DETECTED")
            if t1d_c5b9 is not None and t1d_c5b9 <= Decimal(0):
                reasons.append("FAIL_D1_08_INCONSISTENT_ZERO_UPPER_BOUND_WITH_EVENTS")
        elif t1d_cp["status"] != "PASS":
            reasons.append(
                "HOLD_D1_08_T1D_SERUM_C5B9_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                if t1d_cp["status"] == "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR"
                else f"HOLD_D1_08_T1D_SERUM_C5B9_{t1d_cp['status']}"
            )
            if (
                t1d_c5b9 is not None
                and t1d_k == 0
                and t1d_c5b9 <= t1d_c5b9_u95_pct_cap
                and t1d_cp.get("reported_u95_meets_exact_bound") is False
            ):
                reasons.append("FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND")
        else:
            if t1d_c5b9 is None:
                reasons.append("HOLD_D1_08_T1D_SERUM_C5B9_REPORTED_U95_REQUIRED")
            elif t1d_c5b9 > t1d_c5b9_u95_pct_cap:
                reasons.append("HOLD_D1_08_T1D_SERUM_C5B9_REPORTED_U95_ABOVE_GATE_LIMIT")
            elif t1d_c5b9 <= Decimal(0) or t1d_cp.get("reported_u95_meets_exact_bound") is False:
                reasons.append("FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND")
            elif t1d_cp.get("reported_u95_meets_exact_bound") is None:
                reasons.append("HOLD_D1_08_REPORTED_U95_EXACT_FLOOR_INCONCLUSIVE")
    if (
        carries_c0304 is None
        or (carries_c0304 is True and (c0304_pc is None or hlac03_k is None or hlac03_n is None or hlac03_t is None))
        or (carries_c0304 is False and (ma_pc is None or ma_ok is None))
        or t1d_k is None
        or t1d_n is None
        or t1d_c5b9 is None
    ):
        reasons.append(cell_mdd.get("D1-08", "HOLD_D1_08_MISSING_HLAC0304_OR_C5B9_GY018_INTEGER_KN"))

    # D1-09
    turnaround = _dec("obs_d1_09_autologous_turnaround_weeks")
    if turnaround is None:
        reasons.append(cell_mdd.get("D1-09", "HOLD_D1_09_MISSING_TURNAROUND_OR_GATE0V_PROVENANCE"))
        if effective_l1_arm == "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR":
            reasons.append("HOLD_D1_09_MISSING_TURNAROUND_WEEKS")
    elif turnaround > thr.get("turnaround_max", Decimal("41.0")):
        reasons.append("FAIL_D1_09_AUTOLOGOUS_SEQUENTIAL_TURNAROUND_EXCEEDS_41WKS")

    if invalid_or_oob_keys:
        reasons.append("FAIL_D1_INVALID_OR_OUT_OF_RANGE_OBSERVATION")
    if explicitly_unavailable:
        reasons.append("HOLD_D1_OBSERVATION_EXPLICITLY_UNAVAILABLE")

    required_nums = (
        wt_hlac, bystander, chr6_del, chrm_depth, chrm_het,
        htmpk_leak, elispot, azt_kill,
        l3_ret, l4_ret, fsp1_surv, kir_nk, vtcn1_shed, turnaround,
    )
    if (
        any(v is None for v in required_nums)
        or hlac_cls is None
        or cipsc_src is None
        or effective_l1_arm is None
        or geno_ok is not True
        or not wgs_ok_effective
        or h2h_ok is not True
        or bcor_ok is not True
        or chr13_ok is not True
        or azt_chr19_loss_ok is not True
        or not d1_08_branch_ok
    ):
        reasons.append("HOLD_D1B_AUTO_V2_INCOMPLETE_D101_D109_PANEL")

    dedup: List[str] = []
    for r in reasons:
        if r not in dedup:
            dedup.append(r)
    ordered = (
        [r for r in dedup if r.startswith("FAIL_HARD_STOP_")]
        + [r for r in dedup if r.startswith("FAIL_") and not r.startswith("FAIL_HARD_STOP_")]
        + [r for r in dedup if not r.startswith("FAIL_")]
    )
    reference_criteria_met = len(ordered) == 0

    g0v_cp = g0v.get("cp_result") if isinstance(g0v, dict) else None
    cal_summary: Dict[str, Any] = {
        "d1_04_ap1903": _cal_meta_from_cp(ap_cp, _ap_cal_kw),
        "d1_08_hlac0304": _cal_meta_from_cp(hlac03_cp, _hlac03_cal_kw),
        "d1_08_t1d_serum_c5b9": _cal_meta_from_cp(t1d_cp, _t1d_cal_kw),
        "l1b_gate0v_day28": _cal_meta_from_cp(g0v_cp),
    }
    sampling_prov_by_assay: Dict[str, Any] = {
        "d1_04_ap1903": _ap_cal_kw.get("sampling_independence_provenance"),
        "d1_08_hlac0304": _hlac03_cal_kw.get("sampling_independence_provenance"),
        "d1_08_t1d_serum_c5b9": _t1d_cal_kw.get("sampling_independence_provenance"),
        "l1b_gate0v_day28": g0v_cp.get("sampling_independence_provenance") if isinstance(g0v_cp, dict) else None,
    }
    panel_cov_regime = _aggregate_panel_coverage_regime(cal_summary)

    # Surface ALL assumption-only / default-calibration paths visibly at headline output
    assumption_only_flags: List[str] = _collect_assumption_only_flags(cal_summary)

    is_l1b_arm = bool(is_fullstack and not has_pvid_conflict and effective_l1_arm == "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR")
    if not is_l1b_arm:
        g0v_standalone_status = str(g0v.get("reason") or ("NOT_APPLICABLE_L1A_NO_CD46" if is_fullstack else "UNRESOLVED_PRODUCT_GATE0V_NOT_EVALUATED"))
        g0v_lims_verified = None if (eff_g0v_record is None and not g0v_conflict) else False
    else:
        if g0v.get("cleared") is True and (lot_invalid or lot_conflict):
            g0v_standalone_status = "INVALID_LOT_BINDING_PREVENTS_GATE0V_MATCH"
        elif g0v.get("cleared") is True:
            g0v_standalone_status = "STANDALONE_SUBCHECK_PASSED_LIVE_LIMS_LEDGER_UNVERIFIED"
        else:
            g0v_standalone_status = str(g0v.get("reason") or "HOLD_D1_L1B_GATE0V_UNVERIFIED")
        g0v_lims_verified = False

    is_canonical_eval = bool(
        require_canonical_sha
        and contract["staged_registry_verified"]
        and contract["canonical_sha_verified"]
        and contract["playbook_canonical_sha_verified"]
    )
    if reference_criteria_met:
        top_status = (
            "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING"
            if is_canonical_eval
            else "NON_CANONICAL_DIAGNOSTIC_CRITERIA_MET_UNPINNED"
        )
    else:
        top_status = ordered[0]

    return {
        "gate_id": "STAGE0_GATE_AUTOLOGOUS_D1B_AUTO_C_V2_D101_D109_IIP_0035",
        "decision_id": "DEC-114",
        "contribution_id": "IIP-0035",
        "controlling_contract_block_id": contract["gate_id"],
        "controlling_contract_schema_version": contract["schema_version"],
        "controlling_contract_sha256": contract["sha256"],
        "controlling_staged_registry_verified": contract["staged_registry_verified"],
        "controlling_contract_canonical_sha_verified": contract["canonical_sha_verified"],
        "controlling_contract_on_disk_verified": contract["on_disk_decision_tree_verified"],
        "controlling_contract_in_memory_matches_canonical_sha": contract["in_memory_block_matches_canonical_sha"],
        "controlling_playbook_sec4120_sha256": contract["playbook_sec4120_sha256"],
        "controlling_playbook_footer_verified": contract["playbook_footer_verified"],
        "controlling_playbook_canonical_sha_verified": contract["playbook_canonical_sha_verified"],
        "controlling_contract_cells_loaded": contract["cells_loaded"],
        "controlling_contract_external_refs": contract["external_refs"],
        "evaluation_mode": "CANONICAL_PINNED" if is_canonical_eval else "NON_CANONICAL_DIAGNOSTIC_UNPINNED",
        "evaluation_kind": "REFERENCE_CONTRACT_COMPUTATION_NOT_EMPIRICAL_RELEASE_DECISION",
        "reference_scenario_scope": (
            "SINGLE_LOCUS_L1_DONOR_REFERENCE_SCENARIO_KEY_NOT_ASSEMBLED_MULTILOCUS_MANIFEST"
            if is_fullstack
            else "UNRESOLVED_OR_NON_CANONICAL_PRODUCT_SCOPE"
        ),
        "assembled_multilocus_manifest_bound": False,
        "clinical_or_batch_release_authorized": False,
        "canonical_product_variant_id": pvid or None,
        "effective_locus1_donor_arm": effective_l1_arm,
        "evaluated_lot_id": lot_id_eff if _is_canonical_identifier(lot_id_eff) else None,
        "d1_01_raw_hla_c_input_status": d1_01_raw_hla_c_input_status,
        "d1_01_allele_derivation_mode": d1_01_allele_derivation_mode,
        "d1_02_offgrid_guard_evidence_mode": d1_02_offgrid_guard_evidence_mode,
        "d1_08_allele_branch": d1_08_allele_branch,
        "d1_08_non_c0304_screen_attestation_scope": (
            "EXTERNALLY_ADJUDICATED_CALLER_ATTESTATION_PENDING_NORMATIVE_PROTOCOL"
            if d1_08_allele_branch == "NON_HLAC0304_MATCHED_ALLELE_SCREEN_AND_C5B9"
            else None
        ),
        "d1_08_assay_coverage_limit_note": (
            "Evaluates specified CD8+ peptide and T1D serum C5b-9 complement readouts only; "
            "GY-023 cellular ADCC/NK-ADCC and longitudinal proteostasis remain open requirements."
        ),
        "global_gate0v_standalone_subcheck_status": g0v_standalone_status,
        "global_gate0v_lims_ledger_verified": g0v_lims_verified,
        "panel_coverage_scope": "PER_ASSAY_INDIVIDUAL_COVERAGE_NOT_SIMULTANEOUS_PANEL_FAMILYWISE_GUARANTEE",
        "coverage_regime": panel_cov_regime,
        "rare_event_calibration_summary": cal_summary,
        "sampling_design_provenance_by_assay": sampling_prov_by_assay,
        "assumption_only_path_flags": assumption_only_flags,
        "d1_04_ap1903_u95_str": ap_cp.get("u95_exact_sound_decimal_str") if ap_cp else None,
        "d1_08_hlac0304_u95_pct_str": hlac03_cp.get("u95_exact_sound_decimal_pct_str") if hlac03_cp else None,
        "d1_08_t1d_serum_c5b9_u95_pct_str": t1d_cp.get("u95_exact_sound_decimal_pct_str") if t1d_cp else None,
        "empirical_clearance_status": (
            "MATHEMATICAL_SUBGATE_REFERENCE_CHECK_ONLY_EMPIRICAL_STAGE0_PENDING"
            if (reference_criteria_met and is_canonical_eval)
            else (
                "NON_CANONICAL_DIAGNOSTIC_UNPINNED_NO_RELEASE"
                if reference_criteria_met
                else "HOLD_OR_FAIL_STAGE0_GATE"
            )
        ),
        "reference_criteria_met": reference_criteria_met if is_canonical_eval else False,
        "diagnostic_criteria_met": reference_criteria_met,
        "status": top_status,
        "reasons": ordered,
    }


def compute_rare_event_and_gate_matrix_audit(
    by_id: Dict[str, Dict[str, Any]],
    *,
    repo_root: Path = ROOT,
) -> Dict[str, Any]:
    """Deterministic audit of 80-digit sound Decimal Clopper-Pearson boundaries & standalone D1-01..D1-09 matrix."""
    # 1. Exact 80-digit sound Decimal Clopper-Pearson boundaries
    cp_cases = [
        (
            "d1_04_ap1903_k0_n299573225_HOLD_above_1e8",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_04_AP1903", k=0, n_raw=299_573_225, max_u95_fraction=Decimal("1e-8")
            ),
            "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        ),
        (
            "d1_04_ap1903_k0_n299573226_PASS_at_1e8",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_04_AP1903", k=0, n_raw=299_573_226, max_u95_fraction=Decimal("1e-8")
            ),
            "PASS",
        ),
        (
            "d1_04_ap1903_k1_FAIL_HARD_STOP",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_04_AP1903", k=1, n_raw=300_000_000, max_u95_fraction=Decimal("1e-8")
            ),
            "FAIL_HARD_STOP",
        ),
        (
            "d1_04_ap1903_underallocated_alpha_eta_0p01_HOLD",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_04_AP1903", k=0, n_raw=700_000_000, max_u95_fraction=Decimal("1e-8"),
                eta_recovery_l95=Decimal("0.50"), alpha_obs=Decimal("0.05"), alpha_eta_recovery=Decimal("0.01"),
            ),
            "HOLD_INSUFFICIENT_BONFERRONI_ALPHA_ALLOCATION",
        ),
        (
            "d1_04_ap1903_huge_n_1e308_nonzero_sci_PASS",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_04_AP1903", k=0, n_raw=10**308, max_u95_fraction=Decimal("1e-8")
            ),
            "PASS",
        ),
        (
            "d1_08_k0_n58_HOLD_5p0340pct",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_08_HLAC0304_AND_C5B9", k=0, n_raw=58, max_u95_fraction=Decimal("0.05")
            ),
            "HOLD_INCONCLUSIVE_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        ),
        (
            "d1_08_k0_n59_PASS_4p9508pct",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_08_HLAC0304_AND_C5B9", k=0, n_raw=59, max_u95_fraction=Decimal("0.05")
            ),
            "PASS",
        ),
        (
            "d1_08_k1_n1000_FAIL_HARD_STOP",
            evaluate_sound_decimal_clopper_pearson_zero_event(
                gate_id="D1_08_HLAC0304_AND_C5B9", k=1, n_raw=1000, max_u95_fraction=Decimal("0.05")
            ),
            "FAIL_HARD_STOP",
        ),
    ]
    cp_results: List[Dict[str, Any]] = []
    for cname, res, exp_st in cp_cases:
        if res["status"] != exp_st:
            raise ValueError(f"Rare-event audit case {cname} returned {res['status']!r}, expected {exp_st!r}")
        cp_results.append({"case": cname, **res})
    huge_res = next(r for r in cp_results if r["case"] == "d1_04_ap1903_huge_n_1e308_nonzero_sci_PASS")
    if huge_res["u95_exact_sound_decimal_str"].startswith("0.0"):
        raise ValueError("Huge-n sound Decimal bound must not print as zero")

    # 2. Standalone D1-01..D1-09 reference matrix
    base_obs: Dict[str, Any] = {
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
    valid_g0v = {
        "mcb_lot_id": "MCB_AUTO_LOT_01",
        "cvb_day14_log10_reduction": 3.8,
        "cvb_day14_log10_reduction_l95": 3.2,
        "cvb4_72h_early_titer_log10_pfu_u95": 1.5,
        "day28_vp1_ie1_positive_k": 0,
        "day28_vp1_ie1_cells_scored_n": 6000,
        "gate0v_verdict": "PASS",
    }
    p_l1a = by_id["ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN"]
    p_l1b = by_id["ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR"]
    p_l2 = by_id["ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1"]
    p_l3 = by_id["ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9"]
    p_l4 = by_id["ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y"]

    gate_cases = [
        (
            "l1a_full_panel_pass",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING",
        ),
        (
            "l1b_structured_lot_matched_gate0v_pass",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01", global_gate0v_record=valid_g0v, repo_root=repo_root,
            ),
            "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING",
        ),
        (
            "l1b_bare_boolean_gate0v_rejected_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01", global_gate0v_record=True, repo_root=repo_root,
            ),
            "HOLD_D1_L1B_BARE_BOOLEAN_OR_NON_DICT_GATE0V_REJECTED",
        ),
        (
            "l1b_mismatched_mcb_lot_id_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01",
                global_gate0v_record={**valid_g0v, "mcb_lot_id": "MCB_OTHER_LOT_99"}, repo_root=repo_root,
            ),
            "HOLD_D1_L1B_GATE0V_MCB_LOT_ID_MISMATCH",
        ),
        (
            "l1b_gate0v_viral_failure_routes_to_l1a_FAIL_HARD_STOP",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01",
                global_gate0v_record={**valid_g0v, "day28_vp1_ie1_positive_k": 1}, repo_root=repo_root,
            ),
            "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
        ),
        (
            "l1b_gate0v_k1_with_pending_verdict_still_FAIL_HARD_STOP",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01",
                global_gate0v_record={**valid_g0v, "day28_vp1_ie1_positive_k": 1, "gate0v_verdict": "PENDING"}, repo_root=repo_root,
            ),
            "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
        ),
        (
            "l1b_gate0v_k1_with_missing_day14_still_FAIL_HARD_STOP",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01",
                global_gate0v_record={**valid_g0v, "day28_vp1_ie1_positive_k": 1, "cvb_day14_log10_reduction": None}, repo_root=repo_root,
            ),
            "FAIL_HARD_STOP_D1_L1B_GATE0V_DAY28_VP1_IE1_POSITIVE_EVENT_SWITCH_TO_L1A_IL1RN",
        ),
        (
            "l1b_conflicting_obs_lot_id_vs_explicit_arg_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR", "lot_id": "MCB_OTHER_LOT_B"},
                product_def=p_l1b, evaluated_lot_id="MCB_AUTO_LOT_01",
                global_gate0v_record=valid_g0v, repo_root=repo_root,
            ),
            "HOLD_D1_CONFLICTING_LOT_ID_BINDING",
        ),
        (
            "d1_05_missing_l3_l4_with_only_aggregate_min_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    k: v for k, v in {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"}.items()
                    if k not in ("obs_d1_05_l3_day28_retention_pct", "obs_d1_05_l4_day28_retention_pct")
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_05_MISSING_PER_LOCUS_DURABILITY_OR_FSP1",
        ),
        (
            "d1_05_inconsistent_reported_min_99_vs_88_86_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_05_l3_l4_day28_retention_pct": 99.0},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_05_INCONSISTENT_REPORTED_MIN_L3_L4_RETENTION",
        ),
        (
            "d1_04_k0_n299573225_full_panel_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_04_ap1903_effective_trials_n": 299_573_225},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_04_AP1903_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        ),
        (
            "d1_08_k0_n58_displayed_5p0339_full_panel_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **base_obs,
                    "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
                    "obs_d1_08_hlac0304_scored_targets_n": 58,
                    "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "5.0339",
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_08_HLAC0304_INSUFFICIENT_INDEPENDENT_DENOMINATOR",
        ),
        (
            "d1_08_k0_n59_conservative_unrounded_4p9507609894_PASS",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **base_obs,
                    "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
                    "obs_d1_08_hlac0304_scored_targets_n": 59,
                    "obs_d1_08_hlac0304_tcell_lysis_u95_pct": "4.9507609894",
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING",
        ),
        (
            "l1b_arm_mismatch_with_l1a_product_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1B_CD46WT_CD64T_GATE0V_COMPARATOR"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_01_CANONICAL_PRODUCT_LOCUS1_ARM_MISMATCH",
        ),
        (
            "component_l2_excluded_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(base_obs, product_def=p_l2, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root),
            "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE",
        ),
        (
            "component_l3_excluded_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(base_obs, product_def=p_l3, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root),
            "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE",
        ),
        (
            "component_l4_excluded_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(base_obs, product_def=p_l4, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root),
            "HOLD_D1_COMPONENT_DONOR_EXCLUDED_FROM_FULLSTACK_CLEARANCE",
        ),
        (
            "sp3c_sp4c_c1801_compound_het_excluded_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **base_obs,
                    "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
                    "d1_auto_hlac_sp3c4c_genotype_class": "EXCLUDED_SP3C4C_HLAC07_17_1801_NONFUNCTIONAL_HOMO_OR_COMPOUND_HET_ROUTE_TO_ALLO_Y_OR_D1A",
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_01_EXCLUDED_AH8P1_SP3C4C_NONFUNCTIONAL_HLAC_ROUTE_TO_ALLO_Y_OR_D1A",
        ),
        (
            "d1_04_k1_ap1903_hard_stop_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_04_ap1903_survivor_count_k": 1},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED",
        ),
        (
            "d1_08_k1_with_zero_reported_u95_inconsistent_FAIL_HARD_STOP",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **base_obs,
                    "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
                    "obs_d1_08_hlac0304_lysed_events_k": 1,
                    "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 0.0,
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_HARD_STOP_D1_08_HLAC0304_LYSIS_EVENT_DETECTED",
        ),
        (
            "explicit_unavailable_entry_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": "UNAVAILABLE"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_06_MISSING_KIR3DL1_NK_ASSAY",
        ),
        (
            "d1_02_missing_phased_wgs_with_min_mm3_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {k: v for k, v in {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"}.items() if k != "d1_02_n1_phased_wgs_offgrid_guard_011_cleared"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC",
        ),
        (
            "d1_02_missing_min_mm_with_phased_wgs_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {k: v for k, v in {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"}.items() if k != "obs_d1_02_min_mismatches_to_patient_hlac_e"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_02_INCOMPLETE_PHASED_WGS_QC",
        ),
        (
            "d1_02_exact_decimal_1e_minus_400_del_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_02_chr6p21_interlocus_del_pct": "1e-400"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_02_CHR6P21_INTERLOCUS_DELETION_DETECTED",
        ),
        (
            "d1_06_exact_decimal_8p0000000000000001_nk_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_06_kir3dl1_missing_bw4_nk_lysis_pct": "8.0000000000000001"},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_06_KIR3DL1_MISSING_BW4_NK_LYSIS_EXCEEDS_8PCT",
        ),
        (
            "l1a_missing_lot_id_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"},
                product_def=p_l1a, evaluated_lot_id=None, repo_root=repo_root,
            ),
            "HOLD_D1_INVALID_OR_MISSING_EVALUATED_LOT_ID",
        ),
        (
            "l1a_tampered_product_def_route_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"},
                product_def={**p_l1a, "execution_supported": True, "physical_route_id": "ROUTE_2_OPEN_PORE_PREVASC"},
                evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_INVALID_OR_TAMPERED_PRODUCT_DEFINITION",
        ),
        (
            "d1_08_reported_zero_u95_at_k0_n59_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY", "obs_d1_08_hlac0304_tcell_lysis_u95_pct": 0.0},
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_08_REPORTED_U95_BELOW_EXACT_CLOPPER_PEARSON_BOUND",
        ),
        (
            "malformed_categorical_array_with_k1_still_FAIL_HARD_STOP",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **base_obs,
                    "d1_auto_locus1_donor_arm": ["UNHASHABLE_ARRAY"],
                    "obs_d1_04_ap1903_survivor_count_k": 1,
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_HARD_STOP_D1_04_AP1903_SURVIVOR_DETECTED",
        ),
        (
            "d1_08_c0304_pos_missing_positive_control_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    k: v for k, v in {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"}.items()
                    if k != "d1_08_hlac0304_positive_control_valid"
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_08_HLAC0304_POSITIVE_CONTROL_VALIDITY_REQUIRED",
        ),
        (
            "d1_08_c0304_pos_invalid_positive_control_FAIL",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **base_obs,
                    "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
                    "d1_08_hlac0304_positive_control_valid": False,
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "FAIL_D1_08_HLAC0304_POSITIVE_CONTROL_INVALID",
        ),
        (
            "d1_08_non_c0304_alternate_matched_allele_branch_PASS",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    **{
                        k: v for k, v in base_obs.items()
                        if k not in (
                            "d1_08_hlac0304_positive_control_valid",
                            "obs_d1_08_hlac0304_tcell_lysis_u95_pct",
                            "obs_d1_08_hlac0304_lysed_events_k",
                            "obs_d1_08_hlac0304_scored_targets_n",
                        )
                    },
                    "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY",
                    "d1_08_patient_carries_hlac0304": False,
                    "d1_08_matched_allele_positive_control_valid": True,
                    "d1_08_non_c0304_matched_allele_screen_passed": True,
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "REFERENCE_CRITERIA_MET_EMPIRICAL_RELEASE_PENDING",
        ),
        (
            "d1_08_missing_carrier_status_HOLD",
            evaluate_stage0_d1b_auto_v2_gates_iip0035(
                {
                    k: v for k, v in {**base_obs, "d1_auto_locus1_donor_arm": "D1AUTO_V2_L1A_IL1RN_CD64T_PRIMARY"}.items()
                    if k != "d1_08_patient_carries_hlac0304"
                },
                product_def=p_l1a, evaluated_lot_id="MCB_AUTO_LOT_01", repo_root=repo_root,
            ),
            "HOLD_D1_08_HLAC0304_CARRIER_STATUS_REQUIRED",
        ),
    ]
    gate_results: List[Dict[str, Any]] = []
    for cname, res, exp_st in gate_cases:
        if res["status"] != exp_st:
            raise ValueError(f"Gate matrix case {cname} returned {res['status']!r}, expected {exp_st!r}")
        if res.get("clinical_or_batch_release_authorized") is not False:
            raise ValueError(f"Gate matrix case {cname} must never authorize clinical or batch release")
        if res.get("controlling_contract_sha256") != EXPECTED_DEC114_DECISION_TREE_SHA256:
            raise ValueError(f"Gate matrix case {cname} missing controlling DEC-114 SHA-256")
        if res.get("controlling_playbook_sec4120_sha256") != EXPECTED_DEC114_PLAYBOOK_SEC4120_SHA256:
            raise ValueError(f"Gate matrix case {cname} missing controlling Playbook §4.1.20 SHA-256")
        gate_results.append({"case": cname, **res})

    return {
        "arithmetic_precision": (
            "decimal.Decimal prec=80 strictly sound directed outward bounds (ln/exp at ordinary arguments "
            "and paired alternating Taylor -expm1 enclosure at small x < 1e-6; bounded 2M-bit exact rational "
            "decider + 260-digit outward enclosure; alpha_eta >= 0.05 when eta < 1.0)"
        ),
        "count_subgate_framing": (
            "All zero-event rare-event thresholds (D1-04 AP1903 U95 <= 1.0e-8 at N_eff = 299,573,226 PASS vs 299,573,225 HOLD; "
            "D1-08 HLA-C*03:04 & C5b-9 U95 <= 5.0% at N = 59 PASS [4.9508% display] vs N = 58 HOLD [5.0339% display / 5.0340% conservative ceil]) "
            "represent mathematical count subgate arithmetic only, conditional on validated recovery/eligible-fraction calibration, "
            "and supply no empirical safety clearance."
        ),
        "clopper_pearson_boundary_cases": cp_results,
        "standalone_stage0_gate_reference_matrix": gate_results,
    }


def _resolve_and_verify_d1b_auto_v2_registry(repo_root: Path) -> tuple[Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """Load product_contract, verify FASTA & 5 full_record_sha256s, 45-cell typed contract & CommonMark Playbook §4.1.20."""
    contracts = pc.load_product_contracts(repo_root, require_staged_fastas=True)
    by_id: Dict[str, Dict[str, Any]] = {p["product_variant_id"]: p for p in contracts["products"]}

    # 1. Verify DEC-114 stage0_decision_tree_v1.json canonical SHA-256
    dt_path = repo_root / DECISION_TREE_REL_PATH
    if not dt_path.is_file():
        raise ValueError(f"Missing Stage-0 decision tree: {dt_path}")
    dt_doc = _load_json_strict(dt_path)
    dec114_blk = dt_doc.get("autologous_d1b_auto_c_v2_vanguard_gates_iip_0035")
    if not isinstance(dec114_blk, dict):
        raise ValueError("Missing DEC-114 block autologous_d1b_auto_c_v2_vanguard_gates_iip_0035 in stage0_decision_tree_v1.json")
    dec114_sha = hashlib.sha256(
        json.dumps(dec114_blk, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if dec114_sha != EXPECTED_DEC114_DECISION_TREE_SHA256:
        raise ValueError(
            f"DEC-114 decision-tree block SHA-256 mismatch: found {dec114_sha}, expected {EXPECTED_DEC114_DECISION_TREE_SHA256}"
        )

    # 2. Parse Playbook with CommonMark 0.31.2 parser & require EXACTLY 1 rendered 4.1.20 heading + end-marker
    pb_info = _verify_playbook_sec4120(repo_root, require_canonical_sha=True)
    all_headings = pb_info["all_headings"]
    matched_4120 = pb_info["matched_4120"]
    h4120 = pb_info["h4120"]
    pb_sec_full = pb_info["pb_sec_full"]
    pb_sha = pb_info["pb_sha"]

    # 3. Verify EXACT admitted ISLET1_D1_* namespace (no unpinned synthetic entries!) + roles, non-executable state & ROUTE_3P binding
    d1_namespace_ids = {p for p in by_id if p.startswith("ISLET1_D1_")}
    if d1_namespace_ids != set(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS):
        raise ValueError(
            f"ISLET1_D1_* namespace must contain ONLY the exact 5 pinned IIP-0035 IDs: "
            f"found {sorted(d1_namespace_ids)}, expected {sorted(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS)}"
        )
    resolved_fullstack_ids: List[str] = []
    resolved_excluded_ids: List[str] = []
    resolved_route_codes: set[str] = set()
    for pvid in sorted(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS):
        if pvid not in by_id:
            raise ValueError(f"Loaded product registry missing staged D1b-Auto-C v2 record {pvid}")
        prod = by_id[pvid]
        if prod.get("execution_supported") is not False:
            raise ValueError(f"Loaded product {pvid} must keep execution_supported=False")
        phys_route = prod.get("physical_route_id")
        mapped_code = PHYSICAL_ROUTE_TO_PORTFOLIO_CODE.get(str(phys_route))
        if mapped_code != "R3P":
            raise ValueError(f"Staged D1b-Auto-C v2 record {pvid} must bind strictly to ROUTE_3P_MODULAR_PLANAR_ENVELOPE, got {phys_route!r}")
        resolved_route_codes.add(mapped_code)
        binding = prod.get("staged_d1b_auto_v2_donor_binding") or {}
        role = str(binding.get("role") or "")
        gates_str = str(binding.get("stage0_promotion_gates") or "")
        if role in D1B_FULLSTACK_CANDIDATE_ROLES:
            if "EXCLUDED_FROM_FULLSTACK_CLEARANCE" in gates_str:
                raise ValueError(f"Full-product candidate {pvid} must not be marked excluded from full-stack clearance")
            resolved_fullstack_ids.append(pvid)
        else:
            if "EXCLUDED_FROM_FULLSTACK_CLEARANCE" not in gates_str:
                raise ValueError(f"Component donor {pvid} ({role}) must be marked EXCLUDED_FROM_FULLSTACK_CLEARANCE")
            resolved_excluded_ids.append(pvid)

    # 4. Validate authoritative 45-cell typed product_by_gate_contract (ALL 13 fields checked semantically!)
    pbg = dec114_blk.get("product_by_gate_contract")
    if not isinstance(pbg, dict) or pbg.get("cell_count") != 45 or not isinstance(pbg.get("cells"), list) or len(pbg["cells"]) != 45:
        raise ValueError("DEC-114 product_by_gate_contract must define exactly 45 typed product-by-gate cells")
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
    expected_pairs = {(p, g) for p in pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS for g in EXPECTED_45_GATES}
    if set(cell_map.keys()) != expected_pairs:
        raise ValueError("Typed 45-cell product_by_gate_contract does not match exact 5 products x 9 gates Cartesian set")

    expected_l1a_exact_thresholds: Dict[str, str] = {
        "D1-01": "two_field_HLA_ABC_plus_KIR_genotyping==True;hlac_sp3c4c_class==ELIGIBLE_VL9_COMPETENT_HLAC;exclude_SP3C_SP4C_C07_C17_C1801_homozygotes_and_compound_heterozygotes",
        "D1-02": "OFFGRID_GUARD_011_reject_le_2mm_to_HLAC_E==True;min_mm_to_HLAC_E_valid_range_0_to_25;wt_hlac_retention_pct>=99.0;bystander_ct_pct<0.5;chr6p21_interlocus_del_pct==0.0",
        "D1-03": "cipsc_source_in_{ADSC_SUN_PROTECTED_CIPSC,BLOOD_ERYTHROID_CIPSC_HEAD_TO_HEAD};head_to_head_adsc_vs_blood_comparator_attested==True;TP53_BCOR_PPM1D_20q11.21_1q32.1_MDM4_hits==0;chrm_numt_amplicon_depth_x>=1000.0;chrm_novel_nonsyn_het_pct<1.0",
        "D1-04": "chr13_clybl_independent_from_chr19==True(distinct_chromosomes_joint_failure_unmeasured);G0_beta_htmpk_leak_pct<=0.1;hlac07_ryafsgvay_elispot_above_bg_pct<=0.1;AP1903_k>=1_is_FAIL_HARD_STOP;AP1903_k==0_exact_80digit_sound_Decimal_U95<=1.0e-8(math_count_subgate_only:exact_binomial_N_eff=299573226_PASS_vs_299573225_HOLD;Poisson_ceiling=299573228);azt_kill_measured_after_forced_chr19_loss==True;AZT_dividing_log10_kill>=4.0(admissible_range_0_to_10)",
        "D1-05": "L3_stage7_day28_retention_pct>=80.0;L4_stage7_day28_retention_pct>=80.0;min_L3_L4_stage7_day28_retention_pct>=80.0;RSL3_fsp1_ferroptosis_survival_pct>=80.0",
        "D1-06": "kir3dl1_missing_bw4_autologous_nk_lysis_pct<=8.0",
        "D1-07": "vtcn1_nrd1_soluble_sb7h4_shedding_pct<=10.0",
        "D1-08": "d1_08_patient_carries_hlac0304_bool_or_derived_from_patient_hla_c_alleles(omitted_derives_from_raw_alleles_requiring_both_alleles_in_23_allele_Lin2023_Table1_inventory_for_negative_derivation;contradiction_FAIL;null_or_UNAVAILABLE_or_unlisted_allele_blocks_negative_derivation_HOLD);if_C0304_pos:hlac0304_positive_control_valid==True_and_HLAC0304_assay;if_C0304_neg:matched_allele_positive_control_valid==True_and_matched_allele_cd8_screen_passed==True;HLAC0304_k>=1_is_FAIL_HARD_STOP;HLAC0304_k==0_exact_80digit_sound_Decimal_U95<=5.0%;T1D_serum_C5b9_k>=1_is_FAIL_HARD_STOP;T1D_serum_C5b9_k==0_exact_80digit_sound_Decimal_U95<=5.0%(math_count_subgate_only:N=59_PASS_4.9508%_vs_N=58_HOLD_5.0340%_conservative_ceiling);reported_u95_pct_mandatory_ge_unrounded_exact_bound_and_le_5.0%(missing_reported_u95_HOLD;n<59_returns_denominator_HOLD_unless_underreported_le_5.0%)",
        "D1-09": "autologous_sequential_6locus_gmp_turnaround_weeks<=41.0(admissible_range_1_to_104_weeks);modeled_realistic_range_29_to_41_weeks;L1A_primary_Gate0V_viral_receptor_clean",
    }
    expected_l1b_exact_thresholds: Dict[str, str] = {
        **expected_l1a_exact_thresholds,
        "D1-09": "autologous_sequential_6locus_gmp_turnaround_weeks<=41.0(admissible_range_1_to_104_weeks);structured_global_gate0v_record_required_with_mcb_lot_id==evaluated_lot_id;Gate0V_day14_log10_reduction>=3.0;Gate0V_day14_log10_reduction_l95>=3.0;Gate0V_cvb4_72h_early_titer_log10_pfu_u95<=2.0;Gate0V_day28_VP1_IE1_k>=1_is_FAIL_HARD_STOP;Gate0V_day28_VP1_IE1_k==0_U95<=0.0005(n>=5990);Gate0V_verdict==PASS;never_bare_caller_boolean;standalone_subcheck_not_live_LIMS_ledger_clearance;any_Gate0V_fail_routes_to_L1A_IL1RN",
    }
    expected_component_exact_specs: Dict[tuple[str, str], tuple[str, str]] = {
        ("ISLET1_D1_AUTO_V2_L2_CXCL10_SERPINB9_PDL1", "D1-08"): (
            "COMPOUND_ALL_PASS",
            "EXCLUDED_FROM_FULLSTACK_CLEARANCE;L2_d1_08_patient_carries_hlac0304_bool_or_derived_from_patient_hla_c_alleles(omitted_derives_from_raw_alleles;contradiction_FAIL;null_or_UNAVAILABLE_blocks_derivation_HOLD);if_C0304_pos:hlac0304_positive_control_valid==True_and_L2_HLAC0304_k>=1_is_FAIL_HARD_STOP_and_L2_HLAC0304_k==0_exact_80digit_sound_Decimal_U95<=5.0%_and_reported_u95_pct_ge_unrounded_exact_bound_and_le_5.0%;if_C0304_neg:matched_allele_positive_control_valid==True_and_matched_allele_cd8_screen_passed==True_and_any_HLAC0304_k>=1_is_FAIL_HARD_STOP;optional_adverse_T1D_serum_C5b9_k>=1_triggers_FAIL_HARD_STOP_without_requiring_C5b9_interval_in_L2_subgate;full_product_D1-08_requires_assembled_L1A_or_L1B_stack",
        ),
        ("ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9", "D1-04"): (
            "COMPOUND_ALL_PASS",
            "EXCLUDED_FROM_FULLSTACK_CLEARANCE;L3_chr19_AAVS1_AP1903_k>=1_is_FAIL_HARD_STOP;L3_chr19_AAVS1_AP1903_k==0_exact_80digit_sound_Decimal_U95<=1.0e-8(math_count_subgate_only:exact_binomial_N_eff=299573226_PASS_vs_299573225_HOLD;Poisson_ceiling=299573228);L4_AZT_switch_requires_L4_and_assembled_stack",
        ),
        ("ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9", "D1-05"): (
            ">=",
            "EXCLUDED_FROM_FULLSTACK_CLEARANCE;L3_stage7_day28_retention_pct>=80.0;L4_retention_and_RSL3_require_L4_and_assembled_stack",
        ),
        ("ISLET1_D1_AUTO_V2_L3_AAVS1_VTCN1_UCOE2S_ICASP9", "D1-07"): (
            "<=",
            "EXCLUDED_FROM_FULLSTACK_CLEARANCE;L3_vtcn1_nrd1_soluble_sb7h4_shedding_pct<=10.0",
        ),
        ("ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y", "D1-04"): (
            "COMPOUND_ALL_PASS",
            "EXCLUDED_FROM_FULLSTACK_CLEARANCE;L4_chr13_clybl_independent_from_chr19==True(distinct_chromosomes_joint_failure_unmeasured);L4_G0_beta_htmpk_leak_pct<=0.1;L4_hlac07_ryafsgvay_elispot_above_bg_pct<=0.1;L4_azt_kill_measured_after_forced_chr19_loss==True;L4_AZT_dividing_log10_kill>=4.0(admissible_range_0_to_10);L3_AP1903_switch_requires_L3_and_assembled_stack",
        ),
        ("ISLET1_D1_AUTO_V2_L4_CLYBL_TANDEM_TUBA1B_FSP1_CHS4_CDK1_HTMPK_F105Y", "D1-05"): (
            "COMPOUND_ALL_PASS",
            "EXCLUDED_FROM_FULLSTACK_CLEARANCE;L4_stage7_day28_retention_pct>=80.0;L4_RSL3_fsp1_ferroptosis_survival_pct>=80.0;L3_retention_requires_L3_and_assembled_stack",
        ),
    }

    for pvid in sorted(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS):
        for gid in EXPECTED_45_GATES:
            cell = cell_map[(pvid, gid)]
            if pvid == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN":
                exp_op, exp_disp = EXPECTED_L1A_GATE_OPERATORS[gid], "REQUIRED_EXACT_PRODUCT"
                exp_ext_ref = f"DEC-114:{gid}"
                exp_thr = expected_l1a_exact_thresholds[gid]
            elif pvid == "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR":
                exp_op = EXPECTED_L1B_GATE_OPERATORS[gid]
                exp_disp = "REQUIRED_EXACT_PRODUCT_WITH_LOT_MATCHED_GATE0V" if gid == "D1-09" else "REQUIRED_EXACT_PRODUCT"
                exp_ext_ref = "DEC-114:D1-09+DEC-103:GATE0V" if gid == "D1-09" else f"DEC-114:{gid}"
                exp_thr = expected_l1b_exact_thresholds[gid]
            else:
                if gid in COMPONENT_APPLICABLE_GATES[pvid]:
                    exp_op, exp_thr = expected_component_exact_specs[(pvid, gid)]
                    exp_disp = "COMPONENT_EVIDENCE_ONLY_EXCLUDED_FROM_FULLSTACK"
                    exp_ext_ref = f"DEC-114:{gid}:COMPONENT_ONLY"
                else:
                    exp_op, exp_disp = "NOT_APPLICABLE", "NOT_APPLICABLE_EXCLUDED_FROM_FULLSTACK"
                    exp_ext_ref = f"DEC-114:{gid}:COMPONENT_NOT_APPLICABLE"
                    exp_thr = "NOT_APPLICABLE;EXCLUDED_FROM_FULLSTACK_CLEARANCE"
            if cell["operator"] != exp_op or cell["disposition"] != exp_disp:
                raise ValueError(
                    f"Typed cell operator/disposition mismatch on {pvid}:{gid}: "
                    f"got ({cell['operator']!r}, {cell['disposition']!r}), expected ({exp_op!r}, {exp_disp!r})"
                )
            # Enforce exact gate-scoped external_ref (never points to a different gate!)
            if cell["external_ref"] != exp_ext_ref:
                raise ValueError(
                    f"Typed cell external_ref mismatch on {pvid}:{gid}: got {cell['external_ref']!r}, expected {exp_ext_ref!r}"
                )
            # Enforce fail-closed missing_data_disposition (MUST start with HOLD_ on all 45 cells, NEVER PASS!)
            mdd = cell["missing_data_disposition"]
            if not mdd.startswith("HOLD_") or re.search(r"\bPASS\b", mdd, re.IGNORECASE):
                raise ValueError(
                    f"Typed cell missing_data_disposition must be fail-closed HOLD_* (never PASS) on {pvid}:{gid}: {mdd!r}"
                )
            # Enforce exact full-string threshold equality on ALL 45 cells (no extra smuggled conflicting thresholds!)
            if cell["threshold"] != exp_thr:
                raise ValueError(
                    f"Typed cell threshold exact-equality violation on {pvid}:{gid}: "
                    f"got {cell['threshold']!r}, expected {exp_thr!r}"
                )

    # 5. Contradiction guards across ALL controlling texts, typed cells, and registry records
    rules_dict = dec114_blk.get("rules") or {}
    dt_rules_blob = "\n".join(str(rules_dict.get(k) or "") for k in sorted(rules_dict))
    pb_clean = pb_sec_full.replace("`", "")

    forbidden_patterns = (
        (r"k==0_unconditional_veto", "ambiguous k==0_unconditional_veto shorthand"),
        (r"3,701 bp", "CHK-19 forbidden '3,701 bp' comma literal"),
        (r"n_eff\s*=\s*floor\(", "forbidden floor(n_raw*eta) approximation"),
        (r"(?i)bare\s+boolean\s+is\s+sufficient", "prohibited bare-boolean Gate 0V claim"),
    )
    all_checked_texts: List[tuple[str, str]] = [
        ("decision_tree.rules", dt_rules_blob),
        ("playbook.sec4120", pb_clean),
    ]
    for (cpvid, cgid), cobj in sorted(cell_map.items()):
        all_checked_texts.append((f"typed_cell.{cpvid}:{cgid}", json.dumps(cobj, ensure_ascii=False)))
    for pvid in sorted(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS):
        all_checked_texts.append((f"registry.{pvid}", json.dumps(by_id[pvid], ensure_ascii=False)))
    for dlabel, dtxt in all_checked_texts:
        for pat, desc in forbidden_patterns:
            if re.search(pat, dtxt):
                raise ValueError(f"Contradiction guard tripped ({desc}) in {dlabel}")

    # 6. Direct semantic predicate concordance checked across BOTH typed cells AND decision-tree rules AND Playbook s4.1.20 AND full-product registry!
    #    ZERO skipped predicates (pre-empting job-18 Caveat #2)!
    shared_predicates: List[tuple[str, str, str, str]] = [
        ("D1-01_sp3c4c_exclusion", "D1-01", r"exclude_SP3C_SP4C_C07_C17_C1801_homozygotes_and_compound_heterozygotes", r"SP-3C/SP-4C[^\n]{0,160}C\*18:01[^\n]{0,160}compound\s+heterozygotes"),
        ("D1-02_offgrid_guard_011_thresholds", "D1-02", r"min_mm_to_HLAC_E_valid_range_0_to_25;wt_hlac_retention_pct>=99\.0;bystander_ct_pct<0\.5;chr6p21_interlocus_del_pct==0\.0", r"OFFGRID_GUARD_011[^\n]{0,200}>=\s*99\.0%[^\n]{0,120}<\s*0\.5%[^\n]{0,120}0\.0%"),
        ("D1-03_chrm_1000x_and_1pct", "D1-03", r"chrm_numt_amplicon_depth_x>=1000\.0;chrm_novel_nonsyn_het_pct<1\.0", r">=\s*1,000x\s+NUMT-filtered\s+chrM[^\n]{0,160}<\s*1\.0%"),
        ("D1-04_htmpk_leak_elispot_azt", "D1-04", r"G0_beta_htmpk_leak_pct<=0\.1;hlac07_ryafsgvay_elispot_above_bg_pct<=0\.1;.*azt_kill_measured_after_forced_chr19_loss==True;AZT_dividing_log10_kill>=4\.0", r"<=\s*0\.1%[^\n]{0,140}RYAFSGVAY[^\n]{0,40}<=\s*0\.1%[^\n]{0,160}>=\s*4\.0-log10"),
        ("D1-04_ap1903_explicit_k_ge_1_and_299573226", "D1-04", r"AP1903_k>=1_is_FAIL_HARD_STOP;AP1903_k==0_exact_80digit_sound_Decimal_U95<=1\.0e-8\(math_count_subgate_only:exact_binomial_N_eff=299573226_PASS_vs_299573225_HOLD;Poisson_ceiling=299573228\)", r"k\s*>=\s*1\s+is\s+FAIL_HARD_STOP[^\n]{0,160}<=\s*1\.0e-8[^\n]{0,480}299,573,226[^\n]{0,60}299,573,225\s+HOLD"),
        ("D1-05_per_locus_l3_l4_and_fsp1", "D1-05", r"L3_stage7_day28_retention_pct>=80\.0;L4_stage7_day28_retention_pct>=80\.0;min_L3_L4_stage7_day28_retention_pct>=80\.0;RSL3_fsp1_ferroptosis_survival_pct>=80\.0", r">=\s*80\.0%[^\n]{0,140}min\(L3,\s*L4\)\s*>=\s*80\.0%[^\n]{0,120}RSL3[^\n]{0,100}>=\s*80\.0%"),
        ("D1-06_kir3dl1_nk_8pct", "D1-06", r"kir3dl1_missing_bw4_autologous_nk_lysis_pct<=8\.0", r"KIR3DL1\+[^\n]{0,100}<=\s*8\.0%"),
        ("D1-07_vtcn1_shedding_10pct", "D1-07", r"vtcn1_nrd1_soluble_sb7h4_shedding_pct<=10\.0", r"VTCN1/NRD1[^\n]{0,100}<=\s*10\.0%"),
        ("D1-08_hlac0304_c5b9_explicit_k_ge_1_and_59_vs_58", "D1-08", r"d1_08_patient_carries_hlac0304_bool_or_derived_from_patient_hla_c_alleles\(omitted_derives_from_raw_alleles_requiring_both_alleles_in_23_allele_Lin2023_Table1_inventory_for_negative_derivation;contradiction_FAIL;null_or_UNAVAILABLE_or_unlisted_allele_blocks_negative_derivation_HOLD\);if_C0304_pos:hlac0304_positive_control_valid==True_and_HLAC0304_assay;if_C0304_neg:matched_allele_positive_control_valid==True_and_matched_allele_cd8_screen_passed==True;HLAC0304_k>=1_is_FAIL_HARD_STOP;HLAC0304_k==0_exact_80digit_sound_Decimal_U95<=5\.0%;T1D_serum_C5b9_k>=1_is_FAIL_HARD_STOP;T1D_serum_C5b9_k==0_exact_80digit_sound_Decimal_U95<=5\.0%\(math_count_subgate_only:N=59_PASS_4\.9508%_vs_N=58_HOLD_5\.0340%_conservative_ceiling\)", r"HLA-C\*03:04[^\n]{0,520}C5b-9[^\n]{0,180}k\s*>=\s*1\s+is\s+FAIL_HARD_STOP[^\n]{0,200}<=\s*5\.0%[^\n]{0,280}N\s*=\s*59[^\n]{0,80}4\.9508%[^\n]{0,60}N\s*=\s*58\s+HOLD[^\n]{0,60}5\.0340%"),
        ("D1-09_turnaround_41_weeks", "D1-09", r"autologous_sequential_6locus_gmp_turnaround_weeks<=41\.0", r"<=\s*41\.0\s*weeks"),
    ]
    l1b_specific_predicates: List[tuple[str, str, str, str]] = [
        ("D1-09_l1b_lot_matched_gate0v_record", "D1-09", r"structured_global_gate0v_record_required_with_mcb_lot_id==evaluated_lot_id;Gate0V_day14_log10_reduction>=3\.0;Gate0V_day14_log10_reduction_l95>=3\.0;Gate0V_cvb4_72h_early_titer_log10_pfu_u95<=2\.0;Gate0V_day28_VP1_IE1_k>=1_is_FAIL_HARD_STOP;Gate0V_day28_VP1_IE1_k==0_U95<=0\.0005\(n>=5990\);Gate0V_verdict==PASS;never_bare_caller_boolean", r"global_gate0v_record[^\n]{0,140}mcb_lot_id[^\n]{0,220}>=\s*3\.0-log10[^\n]{0,340}<=\s*0\.0005[^\n]{0,60}>=\s*5,990[^\n]{0,80}gate0v_verdict\s*==\s*\"?PASS\"?[^\n]{0,120}bare\s+caller\s+boolean"),
    ]

    for doc_name, blob in (("decision_tree.rules", dt_rules_blob), ("playbook.sec4120", pb_clean)):
        for pred_name, _gid, _cell_pat, prose_pat in shared_predicates + l1b_specific_predicates:
            if not re.search(prose_pat, blob):
                raise ValueError(f"Typed predicate parity failure ({pred_name}) in {doc_name}")

    product_summaries: List[Dict[str, Any]] = []
    for pvid in sorted(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS):
        prod = by_id[pvid]
        bind = prod["staged_d1b_auto_v2_donor_binding"]
        role = str(bind["role"])
        gates_txt = " ".join(prod.get("pending_assay_gates") or [])
        if pvid == "ISLET1_D1_AUTO_V2_L1A_CIITA_CD64T_IL1RN":
            for pred_name, gid, cell_pat, prose_pat in shared_predicates:
                if not re.search(cell_pat, cell_map[(pvid, gid)]["threshold"]):
                    raise ValueError(f"Typed cell threshold semantic mismatch ({pred_name}) on {pvid}:{gid}")
                if not re.search(prose_pat, gates_txt):
                    raise ValueError(f"Registry predicate mismatch ({pred_name}) on {pvid}")
        elif pvid == "ISLET1_D1_AUTO_V2_L1B_CIITA_CD64T_CD46WT_COMPARATOR":
            for pred_name, gid, cell_pat, prose_pat in shared_predicates + l1b_specific_predicates:
                if not re.search(cell_pat, cell_map[(pvid, gid)]["threshold"]):
                    raise ValueError(f"Typed cell threshold semantic mismatch ({pred_name}) on {pvid}:{gid}")
                if not re.search(prose_pat, gates_txt):
                    raise ValueError(f"Registry predicate mismatch ({pred_name}) on {pvid}")
        else:
            if "COMPONENT_ONLY_EXCLUDED_FROM_FULLSTACK_CLEARANCE" not in gates_txt:
                raise ValueError(f"Component donor {pvid} pending_assay_gates must declare COMPONENT_ONLY_EXCLUDED_FROM_FULLSTACK_CLEARANCE")
            if "Not encoded by this standalone" not in str((prod.get("hla_retention") or {}).get("classical_MHC_I") or ""):
                raise ValueError(f"Component donor {pvid} must distinguish encoded payload from intended host genotype")
            for gid in sorted(COMPONENT_APPLICABLE_GATES[pvid]):
                if "EXCLUDED_FROM_FULLSTACK_CLEARANCE" not in cell_map[(pvid, gid)]["threshold"]:
                    raise ValueError(f"Component cell {pvid}:{gid} missing EXCLUDED_FROM_FULLSTACK_CLEARANCE marker")

        product_summaries.append({
            "product_variant_id": pvid,
            "role": role,
            "fullstack_clearance_eligible": pvid in D1B_FULLSTACK_CANDIDATE_IDS,
            "gate_dispositions": {g: cell_map[(pvid, g)]["disposition"] for g in EXPECTED_45_GATES},
            "typed_cell_operators_and_thresholds_verified": True,
            "registry_predicates_verified_against_typed_contract": True,
        })

    audit = {
        "registry_loaded_via_product_contract": True,
        "dec114_decision_tree_sha256_verified": dec114_sha,
        "dec114_playbook_sec4120_sha256_verified": pb_sha,
        "playbook_commonmark_parser": {
            "rendered_4120_headings_found": len(matched_4120),
            "heading_kind": h4120["kind"],
            "heading_level": h4120["level"],
            "container_depth": h4120["container_depth"],
            "headings_up_to_sec4120": sum(1 for h in all_headings if h["start_offset"] <= h4120["start_offset"]),
        },
        "typed_product_by_gate_contract_audit": {
            "products_checked": len(product_summaries),
            "gates_per_product": len(EXPECTED_45_GATES),
            "cells_verified": len(cell_map),
            "contradiction_guards_passed": len(forbidden_patterns),
            "shared_predicates_verified_zero_skipped": len(shared_predicates),
            "l1b_specific_predicates_verified_zero_skipped": len(l1b_specific_predicates),
            "product_gate_matrix": product_summaries,
        },
        "staged_fasta_sha256_verified": pc.D1B_AUTO_V2_FASTA_SHA256,
        "staged_d1b_auto_v2_record_count": len(pc.STAGED_D1B_AUTO_V2_PRODUCT_IDS),
        "resolved_fullstack_candidate_ids": resolved_fullstack_ids,
        "resolved_excluded_component_ids": resolved_excluded_ids,
        "resolved_d1b_auto_v2_routes_from_registry": sorted(resolved_route_codes),
        "author_declared_sequence_motif_metadata_note": (
            "Sequence motif tallies (0 introns, 0x 2A, 0 uATG, 0 foreign peptides) are author-declared staging "
            "metadata on the product definitions and are not independently recomputed by this verifier."
        ),
    }
    return audit, by_id


def _verify_protected_lims_and_frozen_artifacts(repo_root: Path) -> Dict[str, Any]:
    """Verify exact 0-byte-diff SHA-256 invariance of the 3 protected LIMS files & frozen IIP-0034/Pareto artifacts."""
    lims_verified: Dict[str, Dict[str, Any]] = {}
    for rel, exp_sha in sorted(PROTECTED_LIMS_SHA256.items()):
        fp = repo_root / rel
        if not fp.is_file():
            raise ValueError(f"Missing protected LIMS file for 0-byte invariance verification: {fp}")
        raw = fp.read_bytes()
        got_sha = hashlib.sha256(raw).hexdigest()
        if got_sha != exp_sha:
            raise ValueError(f"Protected LIMS file {rel} mutated: found SHA-256 {got_sha}, expected {exp_sha}")
        lims_verified[rel] = {"sha256": got_sha, "bytes": len(raw), "zero_byte_diff_verified": True}

    ct = _load_json_strict(repo_root / "data" / "protocols" / "wetlab_endpoint_contract_v0.8.0.json")
    lims_counts = {
        "stage0_decision_tree_endpoints": len(ct.get("stage0_decision_tree_endpoints", [])),
        "metadata_keys": len(ct.get("metadata_keys", [])),
        "boolean_policy_keys": len(ct.get("boolean_policy_keys", [])),
    }
    if (lims_counts["stage0_decision_tree_endpoints"], lims_counts["metadata_keys"], lims_counts["boolean_policy_keys"]) != (366, 77, 122):
        raise ValueError(f"Protected LIMS contract counts drifted from 366/77/122: {lims_counts}")

    frozen_checks: Dict[str, str] = {}
    for rel, exp_sha in (
        (PARETO_REL_PATH, EXPECTED_PARETO_SHA256),
        (FROZEN_IIP0034_PORTFOLIO_REL_PATH, EXPECTED_IIP0034_PORTFOLIO_SHA256),
        (FROZEN_IIP0034_VERIF_REL_PATH, EXPECTED_IIP0034_VERIF_SHA256),
    ):
        fp = repo_root / rel
        if not fp.is_file():
            raise ValueError(f"Missing frozen artifact for invariance verification: {fp}")
        got_sha = hashlib.sha256(fp.read_bytes()).hexdigest()
        if got_sha != exp_sha:
            raise ValueError(f"Frozen artifact {rel} drifted: found SHA-256 {got_sha}, expected {exp_sha}")
        frozen_checks[rel] = got_sha

    return {
        "protected_lims_files_zero_byte_diff": lims_verified,
        "protected_lims_contract_counts_unchanged": lims_counts,
        "frozen_upstream_artifacts_sha256_verified": frozen_checks,
    }


def compute_full_verification(repo_root: Path | None = None) -> Dict[str, Any]:
    active_root = ROOT if repo_root is None else Path(repo_root)
    registry_audit, by_id = _resolve_and_verify_d1b_auto_v2_registry(active_root)
    rare_event_and_gate_audit = compute_rare_event_and_gate_matrix_audit(by_id, repo_root=active_root)
    invariance_audit = _verify_protected_lims_and_frozen_artifacts(active_root)

    # Independent replay of the 5-arm autologous portfolio artifact via rank_islet1_stacks
    p5 = rank_mod.compute_5arm_autologous_portfolio_iip0035(repo_root=active_root)
    wavec_4 = list(p5.get("stage0_gated_4arm_autologous_d1b_auto_v2", []))
    wavebc_5 = list(p5.get("stage0_gated_5arm_unified_waveBC_portfolio", []))
    phi_c_cex = list(p5.get("stage0_gated_5arm_phi_c_sensitivity_counterexamples", []))
    if len(wavec_4) != 11 or len(wavebc_5) != 12 or len(phi_c_cex) != 2:
        raise ValueError(
            f"Unexpected 5-arm portfolio row counts: wavec_4={len(wavec_4)}, wavebc_5={len(wavebc_5)}, phi_c_cex={len(phi_c_cex)}"
        )
    all_new_rows = wavec_4 + wavebc_5 + phi_c_cex
    max_4096_err = max(float(r["incumbent_4096_state_max_abs_err"]) for r in all_new_rows)
    if not math.isfinite(max_4096_err) or max_4096_err > 1e-12:
        raise ValueError(f"4,096-state incumbent parity violation: {max_4096_err}")

    alias_map = p5.get("scenario_alias_deduplication_map") or {}
    if float(alias_map.get("alias_full_joint_state_max_abs_diff", 1.0)) > 1e-12:
        raise ValueError(f"1,048,576-state alias verification exceeded 1e-12: {alias_map}")

    # Inadmissible diagonal self-check (unhardened A/B/X barred from R3P -> 0.0% survival)
    pa = (0.88, 0.78, 0.57, 0.69)
    pb = (0.88, 0.79, 0.64, 0.70)
    px = (0.88, 0.82, 0.73, 0.69)
    c_base = (0.88, 0.87, 0.89, 0.85)
    live = rank_mod.compute_common_shock_2pow12_portfolio_tree(pa, pb, px, c_base)
    inadm_diag = rank_mod._run_dp_and_verify_parity_iip0034(
        "IIP0035_INADMISSIBLE_DIAGONAL_SELF_CHECK",
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

    by_4 = {r["scenario"]: r for r in wavec_4}
    by_5 = {r["scenario"]: r for r in wavebc_5}
    p1_exact = p5["autologous_p1_continuous_breakeven_exact"]
    if abs(float(p1_exact["phi_0p25_p1_star"]) - 0.672989) > 1e-6 or abs(float(p1_exact["phi_0p50_p1_star"]) - 0.746615) > 1e-6:
        raise ValueError(f"Continuous p1* break-even drifted: {p1_exact}")

    return {
        "schema_version": "openislet_iip0035_d1b_auto_and_portfolio_verification_v0.8.0",
        "proposal_id": "IIP-0035",
        "decision_id": "DEC-114",
        "license": "CC-BY-4.0",
        "option_2_zero_lims_touch": True,
        "epistemic_status": "DETERMINISTIC_STAGE0_CONTRACT_RARE_EVENT_AND_5ARM_PORTFOLIO_REPLAY_STAGE0_PENDING",
        "live_headline_odds_delta_pp": 0.0,
        "d1b_auto_v2_product_registry_audit": registry_audit,
        "rare_event_and_stage0_gate_matrix_audit": rare_event_and_gate_audit,
        "protected_lims_and_frozen_artifact_invariance": invariance_audit,
        "five_arm_autologous_portfolio_audit": {
            "wave_c_4arm_count": len(wavec_4),
            "wave_bc_5arm_count": len(wavebc_5),
            "phi_c_counterexamples_count": len(phi_c_cex),
            "p1_sweep_points_count": len(p5.get("autologous_p1_breakeven_sweep", [])),
            "max_4096_state_incumbent_abs_err": max_4096_err,
            "inadmissible_diagonal_locked_pct": inadm_diag["locked_pct"],
            "inadmissible_diagonal_modular_k100_pct": inadm_diag["modular_k100_pct"],
            "five_arm_bridge_registration_status": p5["five_arm_bridge_registration_status"],
            "w_c0_negative_control_locked_pct": by_4["W_C0_CTRL_Unedited_Auto_RecurrenceStress@R2 (modeled p3=0.10; motivated by Sutherland1989 & Vendrame2010)"]["locked_pct"],
            "w_c4_p1_0p62_phi25_locked_pct": by_4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["locked_pct"],
            "w_c4_p1_0p62_phi25_registered_route_modular_pct": by_4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_registered_route_only_k100_pct"],
            "w_c4_p1_0p62_phi25_r2_only_modular_pct": by_4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_k100_under_r2_only_expansion_pct"],
            "w_c4_p1_0p62_phi25_class_level_r1_r2_modular_pct": by_4["W_C4_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.25"]["modular_k100_under_class_level_r1_r2_expansion_pct"],
            "w_c5_p1_0p62_phi50_locked_pct": by_4["W_C5_4Arm_D1b@R3P_S2_Screened_p1=0.62_c4u0.92_phi=0.50"]["locked_pct"],
            "w_u1_5arm_zero_h3_phi25_locked_pct": by_5["W_U1_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.70)@3P_c4u0.92_phiYX0.75_25"]["locked_pct"],
            "w_u3_5arm_p1_0p62_registered_route_modular_pct": by_5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["modular_registered_route_only_k100_pct"],
            "w_u3_5arm_p1_0p62_r2_only_modular_pct": by_5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["modular_k100_under_r2_only_expansion_pct"],
            "w_u3_5arm_p1_0p62_class_level_r1_r2_modular_pct": by_5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["modular_k100_under_class_level_r1_r2_expansion_pct"],
            "w_u3_class_level_route_transfer_delta_pp": by_5["W_U3_5Arm_ZeroH3Lift_Y@3P_D1(p1=0.62)@3P_c4u0.92_phiYX0.75_25"]["class_level_route_transfer_delta_pp"],
            "w_u8_5arm_extreme_hostile_phi50_registered_locked_pct": by_5["W_U8_5Arm_EXTREME_HOSTILE_Y@3P(0.68)_D1(0.55,0.60)@R2_c4u90_phi50"]["locked_pct"],
            "w_u8_5arm_extreme_hostile_phi50_hypothetical_r2_locked_pct": by_5["W_U8_5Arm_EXTREME_HOSTILE_Y@3P(0.68)_D1(0.55,0.60)@R2_c4u90_phi50"]["locked_under_hypothetical_r2_transfer_pct"],
            "w_u9_5arm_ultra_hostile_phi50_registered_locked_pct": by_5["W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90"]["locked_pct"],
            "w_u9_5arm_ultra_hostile_phi50_hypothetical_r2_locked_pct": by_5["W_U9_5Arm_ULTRA_HOSTILE_phiYX1.0_Y@3P(0.68)_D1(0.55,0.60)@R2_c90"]["locked_under_hypothetical_r2_transfer_pct"],
            "w_u9_phi_c_0p625_registered_locked_pct": phi_c_cex[0]["locked_pct"],
            "w_u9_phi_c_0p625_hypothetical_r2_locked_pct": phi_c_cex[0]["locked_under_hypothetical_r2_transfer_pct"],
            "w_u9_phi_c_0p70_registered_locked_pct": phi_c_cex[1]["locked_pct"],
            "w_u9_phi_c_0p70_hypothetical_r2_locked_pct": phi_c_cex[1]["locked_under_hypothetical_r2_transfer_pct"],
            "w_u11_5arm_nested_upside_phi25_locked_pct": by_5["W_U11_5Arm_NestedY(90/92)_Y@Hyb3P_D1(p1=0.70)@3P_c4u0.92_phi25"]["locked_pct"],
            "autologous_p1_continuous_breakeven_exact": p1_exact,
            "scenario_alias_deduplication_map": alias_map,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Deterministic autologous D1b-Auto-C v2, 45-cell contract, rare-event & 5-arm portfolio verifier (IIP-0035 / DEC-114)."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true", help=f"Write canonical JSON artifact to {DEFAULT_OUT_PATH}")
    mode.add_argument("--check", action="store_true", help=f"Verify byte-identical match against {DEFAULT_OUT_PATH} and {PORTFOLIO_5ARM_PATH}")
    mode.add_argument(
        "--eval-json",
        type=Path,
        default=None,
        help="Evaluate a raw Stage-0 observation JSON file preserving exact decimal numeric lexemes (parse_float=Decimal)",
    )
    args = parser.parse_args()

    if args.eval_json is not None:
        # Cryptographically verify product registry, 45-cell DEC-114 contract, and Playbook §4.1.20 BEFORE evaluating!
        _reg_audit, by_id = _resolve_and_verify_d1b_auto_v2_registry(ROOT)
        dt_doc = _load_json_strict(ROOT / DECISION_TREE_REL_PATH)
        dec114_blk = dt_doc["autologous_d1b_auto_c_v2_vanguard_gates_iip_0035"]
        req = parse_stage0_observation_json_exact(args.eval_json.read_text(encoding="utf-8"))
        if "observations" in req and not isinstance(req["observations"], dict):
            raise SystemExit("Invalid --eval-json payload: 'observations' must be a JSON object when provided (null/non-object rejected).")
        obs = req["observations"] if isinstance(req.get("observations"), dict) else req
        if obs is not req:
            _allowed_env_keys = frozenset({
                "observations", "product_variant_id", "product_def",
                "evaluated_lot_id", "mcb_lot_id", "lot_id", "batch_id", "mcb_batch_id",
                "global_gate0v_record",
            })
            _unknown_env_keys = sorted(k for k in req.keys() if k not in _allowed_env_keys)
            if _unknown_env_keys:
                raise SystemExit(
                    f"Invalid --eval-json envelope: unexpected top-level key(s) outside 'observations': {_unknown_env_keys}"
                )
        # Pass envelope fields directly into evaluate_stage0_d1b_auto_v2_gates_iip0035() WITHOUT pre-selecting
        # a product_def from distinct_pvids, guaranteeing 100% branch-for-branch parity with direct Python calls!
        env_pdef = req.get("product_def") if (obs is not req and "product_def" in req) else None
        env_pvid = req.get("product_variant_id") if obs is not req else None
        env_lots = (
            {
                k: req[k]
                for k in ("mcb_lot_id", "lot_id", "batch_id", "mcb_batch_id")
                if k in req and req[k] is not None
            }
            if obs is not req
            else None
        )
        res = evaluate_stage0_d1b_auto_v2_gates_iip0035(
            obs,
            product_def=env_pdef,
            product_variant_id=env_pvid,
            evaluated_lot_id=(req.get("evaluated_lot_id") if obs is not req else None),
            extra_lot_bindings=env_lots,
            global_gate0v_record=(
                req["global_gate0v_record"]
                if (obs is not req and "global_gate0v_record" in req)
                else _MISSING_ARG
            ),
            decision_tree_block=dec114_blk,
            repo_root=ROOT,
            require_canonical_sha=True,
        )
        print(json.dumps(res, indent=2, ensure_ascii=False, allow_nan=False))
        return

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
        if not PORTFOLIO_5ARM_PATH.is_file():
            raise SystemExit(f"Missing shipped 5-arm autologous portfolio artifact: {PORTFOLIO_5ARM_PATH}")
        shipped_p5_bytes = PORTFOLIO_5ARM_PATH.read_bytes()
        expected_p5 = rank_mod.compute_5arm_autologous_portfolio_iip0035(repo_root=ROOT)
        expected_p5_bytes = (
            json.dumps(expected_p5, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        ).encode("utf-8")
        shipped_p5_sha = hashlib.sha256(shipped_p5_bytes).hexdigest()
        expected_p5_sha = hashlib.sha256(expected_p5_bytes).hexdigest()
        if shipped_p5_bytes != expected_p5_bytes or shipped_p5_sha != expected_p5_sha:
            raise SystemExit(
                f"Shipped 5-arm autologous portfolio byte/SHA-256 drift detected in {PORTFOLIO_5ARM_PATH}: "
                f"found {shipped_p5_sha}, expected {expected_p5_sha}"
            )
        print(
            f"OK: {DEFAULT_OUT_PATH} (SHA-256={expected_sha}) and "
            f"{PORTFOLIO_5ARM_PATH} (SHA-256={expected_p5_sha}) match deterministic replay byte-for-byte."
        )
    else:
        print(serialized)


if __name__ == "__main__":
    main()
