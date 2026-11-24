#!/usr/bin/env python3
"""
Deterministic Stage-0 Physics, Multi-Mesh 3D FVM & Respirometry Verifier (IIP-0033 / DEC-112)
==============================================================================================
Self-contained replay script for OpenIslet public and audit checkouts. Depends ONLY on
`scripts/run_layer4_tissue_and_metabolic_twin.py`, `scripts/rank_islet1_stacks.py`,
`numpy`, and the Python standard library.

Run with:
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0033_bilayer_and_resp1.py
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0033_bilayer_and_resp1.py --write
  PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/verify_iip0033_bilayer_and_resp1.py --check

Verifies deterministically:
  1. Static Agilent 101122-100 islet capture microplate O2 mass balance across ALL THREE
     published Agilent chamber volumes (16.60 uL FluxPak 103518-100, 16.87 uL brochure
     5991-8657EN p. 9, and 9.87 uL XF24 support FAQ -> 6/10 IEQ caps where 9/16 IEQ drops to
     20.96/105.84 mmHg), plus 14.94 uL (-10% off 16.60 uL -> 9/16 IEQ with +0.61 mmHg air-FCCP
     margin), 14.00 uL screen-displaced sensitivity (9/15 IEQ), and obsolete 22.70 uL XF24 V7
     (Gerencser 2009 PMID:19555051, mouse-liver mitochondria n=4) under KRBH alpha = 1.27e-9
     mol/(cm^3*mmHg) vs low-salt/BSA (1.34e-9) and pure water (1.39e-9), alongside 125 um stock
     capture-screen pore retention calibration requirements.
  2. Primary flow-through 5% O2 perifusion modeled on Sweet et al. (2008 PMID:18021279; 750
     handpicked islets, glucose/carbachol) methodological precedent: Stage-0 modeled 200 IEQ
     @ 200 uL/min -> 32.07 / 30.28 / 28.49 / 26.70 mmHg across 1x-2.5x (SNR = 17.9x) with exact
     unrounded float threshold comparisons.
  3. Fixed-pitch (p = 104.1604 um, 1,500.0 IEQ/cm^2 nominal) 2x2 mixed-lot 65 um staggered
     bilayer multi-mesh 3D FVM (nx = 40, 60, 80) calling live OpenIslet l4 Newton-PCG solver:
       - PASS envelopes: Route 2 140 um staggered Design Z (Chr30 No-FBR & Chr35 15 um 3x FBR),
         Route 1 390 um CaO2 (nominal & coupled +10% OCR + 16.5 um standoff stress), Route 1 350 um
         nominal & pure +10% OCR at fixed 15.0 um standoff (15.03-15.06 mmHg), and p95 = 70 um
         on Route 1 390 um under pure +10% OCR at fixed 15.0 um standoff (15.21-15.24 mmHg).
       - HOLD/FAIL proofs: Route 1 350 um under coupled +10% OCR + 16.5 um standoff (14.62-14.64 mmHg
         < 15 mmHg), p95 = 70 um on Route 1 390 um under coupled +10% OCR + 16.5 um standoff
         (14.89-14.93 mmHg < 15 mmHg), 160 um aligned backup (exact 1,500, fixed 104.16 um mold, and
         3.1464 cap), and unconstrained d_rms3 = 68 um mold overload (1,717.4 IEQ/cm^2).
       - Explicit Day-0 Avascular (k_perf = 0, D_gel_acute) vs Chronic Perfused separation.
       - GY-027 monolayer fixed-pitch vs re-derived 75 um packing gap check (-1.35 vs +16.29 um).
  4. GY-028 cylindrical pre-vascularized channel traps at clinical 500,000 IEQ scale:
       - R = 500 um (1 mm ID) Day-0 core anoxia (0.0000 mmHg).
       - R = 120 um 30 um dense thrombus sleeve collapse (1.4106 / 0.5908 mmHg).
       - R = 100-120 um track length, surgical footprint, pull-out force, and VPIO jamming.
  5. Live OQ-016 2^12 baseline parity check (+0.00 pp live odds invariance).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import rank_islet1_stacks as rank_mod  # noqa: E402
import run_layer4_tissue_and_metabolic_twin as l4  # noqa: E402

DEFAULT_OUT_PATH = ROOT / "data" / "extracted" / "iip0033_bilayer_and_respirometry_verification_v0.8.0.json"

# Canonical OpenIslet Layer-4 physical constants
ALPHA_LOW_SALT = float(l4.ALPHA_BUNSEN_MOL_CM3_MMHG)  # 1.34e-9 mol/(cm^3 * mmHg)
ALPHA_PURE_WATER_37C = 1.39e-9                        # 1.39e-9 mol/(cm^3 * mmHg)
ALPHA_KRBH = 1.27e-9                                  # 1.27e-9 mol/(cm^3 * mmHg)
D_GEL_ACUTE = float(l4.D_O2_HYDROGEL_CM2_S)           # 0.940254e-5 cm^2/s
D_GEL_CHRON = 1.50e-5                                 # cm^2/s
D_ISLET = float(l4.D_O2_ISLET_TISSUE_CM2_S)           # 1.30e-5 cm^2/s
VMAX_PURE = float(l4.VMAX_ISLET_PURE_MOL_CM3_S)       # 3.4e-8 mol/(cm^3 * s)
VOL_IEQ = float(l4.VOL_1_IEQ_CM3)                     # 1.7671458676442586e-6 cm^3
K_GEL_CHRON = D_GEL_CHRON * ALPHA_LOW_SALT
K_ISLET = D_ISLET * ALPHA_LOW_SALT
K_PERF_CHRON = 0.14 * 0.97                            # 0.1358 s^-1
DOSE_IEQ = 500_000.0
VOL_500K_CM3 = DOSE_IEQ * VOL_IEQ                     # 0.883573 cm^3


def _stable_mass_balance_err(raw_err: float) -> float:
    """Assert strict < 1e-10 mass-balance conservation while clamping sub-1e-10 BLAS/NumPy ULP noise for byte-exact cross-runtime --check."""
    err = abs(float(raw_err))
    assert math.isfinite(err) and err < 1e-10, f"Mass balance residual exceeded 1e-10: {raw_err!r}"
    return 0.0 if err < 1e-10 else round(err, 10)


def _is_valid_finite_number_iip0033(val: Any) -> bool:
    return not isinstance(val, bool) and isinstance(val, (int, float)) and math.isfinite(float(val))


def _require_positive_finite_iip0033(val: Any, name: str) -> float:
    if not _is_valid_finite_number_iip0033(val) or float(val) <= 0.0:
        raise ValueError(f"{name} must be a strictly positive finite number, got {val!r}")
    return float(val)


def compute_static_xfe24_101122_100_cap_iip0033(
    *,
    volume_uL: float = 14.94,
    alpha_mol_cm3_mmhg: float = 1.27e-9,
    ocr_whole_prep_pmol_min_ieq: float = 4.5479,
    read_min: float = 3.0,
) -> Dict[str, Any]:
    """Compute exact Agilent 101122-100 static microchamber O2 depletion and max IEQ/well caps (unrounded threshold comparison; IIP-0033 / DEC-112)."""
    volume_uL = _require_positive_finite_iip0033(volume_uL, "volume_uL")
    alpha_mol_cm3_mmhg = _require_positive_finite_iip0033(alpha_mol_cm3_mmhg, "alpha_mol_cm3_mmhg")
    ocr_whole_prep_pmol_min_ieq = _require_positive_finite_iip0033(ocr_whole_prep_pmol_min_ieq, "ocr_whole_prep_pmol_min_ieq")
    read_min = _require_positive_finite_iip0033(read_min, "read_min")
    p_5pct = 0.05 * (760.0 - 47.0)       # 35.6500 mmHg
    p_21pct = 0.2095 * (760.0 - 47.0)    # 149.3735 mmHg
    vol_cm3 = volume_uL * 1e-3
    ocr_mol_min_ieq = ocr_whole_prep_pmol_min_ieq * 1e-12
    dp_basal_per_ieq = (ocr_mol_min_ieq * read_min) / (vol_cm3 * alpha_mol_cm3_mmhg)
    return {
        "volume_uL": float(volume_uL),
        "alpha_mol_cm3_mmhg": float(alpha_mol_cm3_mmhg),
        "ocr_whole_prep_pmol_min_ieq": float(ocr_whole_prep_pmol_min_ieq),
        "read_min": float(read_min),
        "dp_basal_3min_mmhg_per_ieq": round(dp_basal_per_ieq, 4),
        "dp_basal_3min_mmhg_per_ieq_unrounded": float(dp_basal_per_ieq),
        "max_ieq_5pct_basal_ge25": int(math.floor((p_5pct - 25.0) / dp_basal_per_ieq)),
        "max_ieq_5pct_gsis_1p5x_ge25": int(math.floor((p_5pct - 25.0) / (1.5 * dp_basal_per_ieq))),
        "max_ieq_5pct_gsis_2p0x_ge25": int(math.floor((p_5pct - 25.0) / (2.0 * dp_basal_per_ieq))),
        "max_ieq_5pct_fccp_2p5x_ge25": int(math.floor((p_5pct - 25.0) / (2.5 * dp_basal_per_ieq))),
        "max_ieq_21pct_fccp_2p5x_ge120": int(math.floor((p_21pct - 120.0) / (2.5 * dp_basal_per_ieq))),
    }


def compute_sweet2008_flow_through_po2_out_iip0033(
    *,
    n_ieq: float = 200.0,
    flow_uL_min: float = 200.0,
    alpha_mol_cm3_mmhg: float = 1.27e-9,
    ocr_whole_prep_pmol_min_ieq: float = 4.5479,
    noise_floor_mmhg: float = 0.20,
) -> Dict[str, Any]:
    """Compute flow-through 5% O2 perifusion steady-state outflow pO2 (Sweet 2008 methodological precedent; exact unrounded threshold comparison; IIP-0033 / DEC-112)."""
    n_ieq = _require_positive_finite_iip0033(n_ieq, "n_ieq")
    flow_uL_min = _require_positive_finite_iip0033(flow_uL_min, "flow_uL_min")
    alpha_mol_cm3_mmhg = _require_positive_finite_iip0033(alpha_mol_cm3_mmhg, "alpha_mol_cm3_mmhg")
    ocr_whole_prep_pmol_min_ieq = _require_positive_finite_iip0033(ocr_whole_prep_pmol_min_ieq, "ocr_whole_prep_pmol_min_ieq")
    noise_floor_mmhg = _require_positive_finite_iip0033(noise_floor_mmhg, "noise_floor_mmhg")
    p_5pct = 0.05 * (760.0 - 47.0)  # 35.6500 mmHg
    q_cm3_min = flow_uL_min * 1e-3
    dp_basal = (n_ieq * ocr_whole_prep_pmol_min_ieq * 1e-12) / (q_cm3_min * alpha_mol_cm3_mmhg)
    p_b_raw = p_5pct - dp_basal
    p_15_raw = p_5pct - 1.5 * dp_basal
    p_20_raw = p_5pct - 2.0 * dp_basal
    p_25_raw = p_5pct - 2.5 * dp_basal
    return {
        "n_ieq": float(n_ieq),
        "flow_uL_min": float(flow_uL_min),
        "alpha_mol_cm3_mmhg": float(alpha_mol_cm3_mmhg),
        "delta_po2_basal_mmhg": round(dp_basal, 2),
        "snr_over_noise_floor": round(dp_basal / noise_floor_mmhg, 2),
        "p_out_5pct_basal_1p0x_mmhg": round(p_b_raw, 2),
        "p_out_5pct_gsis_1p5x_mmhg": round(p_15_raw, 2),
        "p_out_5pct_gsis_2p0x_mmhg": round(p_20_raw, 2),
        "p_out_5pct_fccp_2p5x_mmhg": round(p_25_raw, 2),
        "p_out_5pct_fccp_2p5x_unrounded_mmhg": float(p_25_raw),
        "all_4_states_ge_25mmhg": bool(p_25_raw >= 25.0),
        "passes_primary_spec": bool(p_25_raw >= 25.0 and dp_basal >= max(1.00, 5.0 * noise_floor_mmhg)),
    }


def compute_type7_finite_sample_cubic_mean_bounds_iip0033(
    *,
    d_min: float,
    p05: float,
    p95: float,
    d_max: float,
    sample_n: int | None = None,
) -> Dict[str, Any]:
    """Exact joint-constraint Hyndman-Fan Type-7 cubic-mean bounds over sorted d_min = x_0 <= ... <= x_{N-1} = d_max simultaneously satisfying both Type-7 linear interpolation equalities for any N >= 2 (IIP-0033 / DEC-112)."""
    d_min_f = _require_positive_finite_iip0033(d_min, "d_min")
    p05_f = _require_positive_finite_iip0033(p05, "p05")
    p95_f = _require_positive_finite_iip0033(p95, "p95")
    d_max_f = _require_positive_finite_iip0033(d_max, "d_max")
    if not (d_min_f <= p05_f <= p95_f <= d_max_f):
        raise ValueError(f"Non-monotonic quantiles: d_min={d_min_f}, p05={p05_f}, p95={p95_f}, d_max={d_max_f}")
    if sample_n is None:
        asymp_lb = (0.05 * (d_min_f ** 3) + 0.90 * (p05_f ** 3) + 0.05 * (p95_f ** 3)) ** (1.0 / 3.0)
        asymp_ub = (0.05 * (p05_f ** 3) + 0.90 * (p95_f ** 3) + 0.05 * (d_max_f ** 3)) ** (1.0 / 3.0)
        ext_lb_10k = ((9999.0 * (d_min_f ** 3) + 1.0 * (d_max_f ** 3)) / 10000.0) ** (1.0 / 3.0)
        ext_ub_10k = ((1.0 * (d_min_f ** 3) + 9999.0 * (d_max_f ** 3)) / 10000.0) ** (1.0 / 3.0)
        return {
            "sample_n": None,
            "feasible": True,
            "d_rms3_lb": float(max(asymp_lb, ext_lb_10k)),
            "d_rms3_ub": float(min(asymp_ub, ext_ub_10k)),
        }
    if isinstance(sample_n, bool) or not isinstance(sample_n, int) or sample_n < 2:
        raise ValueError(f"sample_n must be an integer >= 2 when provided, got {sample_n!r}")
    n = int(sample_n)
    tol = 1e-9 * max(1.0, abs(d_max_f))
    i05, r05 = divmod(n - 1, 20)
    i95, r95 = divmod(19 * (n - 1), 20)
    f05 = r05 / 20.0
    f95 = r95 / 20.0

    # Case 1: N == 2 (0 free variables: x_0 = d_min, x_1 = d_max)
    if n == 2:
        q05_exp = (1.0 - f05) * d_min_f + f05 * d_max_f
        q95_exp = (1.0 - f95) * d_min_f + f95 * d_max_f
        if abs(q05_exp - p05_f) > tol or abs(q95_exp - p95_f) > tol:
            return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
        val = ((d_min_f ** 3 + d_max_f ** 3) / 2.0) ** (1.0 / 3.0)
        return {"sample_n": n, "feasible": True, "d_rms3_lb": float(val), "d_rms3_ub": float(val)}

    # Case 2: N == 3 (0 free variables: x_1 constrained by BOTH equations)
    if n == 3:
        x1_a = (p05_f - (1.0 - f05) * d_min_f) / f05
        x1_b = (p95_f - f95 * d_max_f) / (1.0 - f95)
        anchor_tol_3 = tol / min(f05, 1.0 - f95)
        if abs(x1_a - x1_b) > 2.0 * anchor_tol_3:
            return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
        x1 = 0.5 * (x1_a + x1_b)
        if x1 < d_min_f - anchor_tol_3 or x1 > d_max_f + anchor_tol_3:
            return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
        x1 = min(max(x1, d_min_f), d_max_f)
        if (
            abs((1.0 - f05) * d_min_f + f05 * x1 - p05_f) > tol
            or abs((1.0 - f95) * x1 + f95 * d_max_f - p95_f) > tol
        ):
            return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
        val = ((d_min_f ** 3 + x1 ** 3 + d_max_f ** 3) / 3.0) ** (1.0 / 3.0)
        return {"sample_n": n, "feasible": True, "d_rms3_lb": float(val), "d_rms3_ub": float(val)}

    # Case 3: 4 <= N <= 20 (i05 == 0, i95 == n - 2; both interior anchors x_1 and x_{n-2} uniquely fixed;
    # propagate quantile tolerance tol through 1/f05 and 1/(1-f95) before final exact quantile residual check)
    if n <= 20:
        y1 = (p05_f - (1.0 - f05) * d_min_f) / f05
        y2 = (p95_f - f95 * d_max_f) / (1.0 - f95)
        anchor_tol = tol / min(f05, 1.0 - f95)
        if y1 < d_min_f - anchor_tol or y2 > d_max_f + anchor_tol or y1 > y2 + 2.0 * anchor_tol:
            return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
        y1 = min(max(y1, d_min_f), d_max_f)
        y2 = min(max(y2, d_min_f), d_max_f)
        if y1 > y2:
            y1 = y2 = 0.5 * (y1 + y2)
        if (
            abs((1.0 - f05) * d_min_f + f05 * y1 - p05_f) > tol
            or abs((1.0 - f95) * y2 + f95 * d_max_f - p95_f) > tol
        ):
            return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
        s_min = (d_min_f ** 3) + (n - 3) * (y1 ** 3) + (y2 ** 3) + (d_max_f ** 3)
        s_max = (d_min_f ** 3) + (y1 ** 3) + (n - 3) * (y2 ** 3) + (d_max_f ** 3)
        return {
            "sample_n": n,
            "feasible": True,
            "d_rms3_lb": float((s_min / float(n)) ** (1.0 / 3.0)),
            "d_rms3_ub": float((s_max / float(n)) ** (1.0 / 3.0)),
        }

    # Case 4: N >= 21 with r05 == 0 (exact integer quantile indices: f05 == f95 == 0; 4 fixed anchors)
    if r05 == 0:
        s_min = i05 * (d_min_f ** 3) + (i95 - i05) * (p05_f ** 3) + (n - 1 - i95) * (p95_f ** 3) + (d_max_f ** 3)
        s_max = (d_min_f ** 3) + i05 * (p05_f ** 3) + (i95 - i05) * (p95_f ** 3) + (n - 1 - i95) * (d_max_f ** 3)
        return {
            "sample_n": n,
            "feasible": True,
            "d_rms3_lb": float((s_min / float(n)) ** (1.0 / 3.0)),
            "d_rms3_ub": float((s_max / float(n)) ** (1.0 / 3.0)),
        }

    # Case 5: N >= 22 with r05 > 0 (6 distinct anchors [0, i05, i05+1, i95, i95+1, n-1], 2 free variables (u, v) = (x[i05+1], x[i95]))
    u_lo = p05_f
    u_hi = min(p95_f, (p05_f - (1.0 - f05) * d_min_f) / f05)
    v_lo = max(p05_f, (p95_f - f95 * d_max_f) / (1.0 - f95))
    v_hi = p95_f
    if u_lo > u_hi + tol or v_lo > v_hi + tol or u_lo > v_hi + tol:
        return {"sample_n": n, "feasible": False, "d_rms3_lb": None, "d_rms3_ub": None}
    u_hi = max(u_lo, u_hi)
    v_lo = min(v_lo, v_hi)

    w0_min, w1_min, w2_min, w3_min, w4_min, w5_min = (
        i05, 1, i95 - (i05 + 1), 1, (n - 1) - (i95 + 1), 1
    )
    w0_max, w1_max, w2_max, w3_max, w4_max, w5_max = (
        1, i05, 1, i95 - (i05 + 1), 1, (n - 1) - (i95 + 1)
    )

    def _y1(u_val: float) -> float:
        return (p05_f - f05 * u_val) / (1.0 - f05)

    def _y4(v_val: float) -> float:
        return (p95_f - (1.0 - f95) * v_val) / f95

    def _eval_min(u_val: float, v_val: float) -> float:
        return (
            w0_min * (d_min_f ** 3)
            + w1_min * (_y1(u_val) ** 3)
            + w2_min * (u_val ** 3)
            + w3_min * (v_val ** 3)
            + w4_min * (_y4(v_val) ** 3)
            + w5_min * (d_max_f ** 3)
        )

    def _eval_max(u_val: float, v_val: float) -> float:
        return (
            w0_max * (d_min_f ** 3)
            + w1_max * (_y1(u_val) ** 3)
            + w2_max * (u_val ** 3)
            + w3_max * (v_val ** 3)
            + w4_max * (_y4(v_val) ** 3)
            + w5_max * (d_max_f ** 3)
        )

    verts: List[Tuple[float, float]] = []
    for u_cand, v_cand in (
        (u_lo, v_lo),
        (u_lo, v_hi),
        (u_hi, v_hi),
        (u_hi, v_lo),
        (u_hi, u_hi),
        (v_lo, v_lo),
        (u_lo, u_lo),
        (v_hi, v_hi),
    ):
        if (
            u_lo - tol <= u_cand <= u_hi + tol
            and v_lo - tol <= v_cand <= v_hi + tol
            and u_cand <= v_cand + tol
        ):
            uc = min(max(u_cand, u_lo), u_hi)
            vc = min(max(v_cand, v_lo), v_hi)
            if uc > vc:
                uc = vc = 0.5 * (uc + vc)
            verts.append((uc, vc))

    s_max = max(_eval_max(u_val, v_val) for u_val, v_val in verts)

    rho_u = math.sqrt((w2_min * (1.0 - f05)) / (w1_min * f05))
    u_star = p05_f / (f05 + (1.0 - f05) * rho_u)
    rho_v = math.sqrt((w3_min * f95) / (w4_min * (1.0 - f95)))
    v_star = p95_f / ((1.0 - f95) + f95 * rho_v)

    u_c = min(max(u_star, u_lo), u_hi)
    v_c = min(max(v_star, v_lo), v_hi)

    min_cands: List[Tuple[float, float]] = list(verts)
    if u_c <= v_c + tol:
        uc_f, vc_f = u_c, v_c
        if uc_f > vc_f:
            uc_f = vc_f = 0.5 * (uc_f + vc_f)
        min_cands.append((uc_f, vc_f))
    for uc_e, vc_e in ((u_c, v_lo), (u_c, v_hi), (u_lo, v_c), (u_hi, v_c)):
        if (
            u_lo - tol <= uc_e <= u_hi + tol
            and v_lo - tol <= vc_e <= v_hi + tol
            and uc_e <= vc_e + tol
        ):
            uc_cl = min(max(uc_e, u_lo), u_hi)
            vc_cl = max(min(max(vc_e, v_lo), v_hi), uc_cl)
            min_cands.append((uc_cl, vc_cl))
    t_lo = max(u_lo, v_lo)
    t_hi = min(u_hi, v_hi)
    if t_lo <= t_hi + tol:
        t_lo_c = min(t_lo, t_hi)
        t_hi_c = max(t_lo, t_hi)
        min_cands.append((t_lo_c, t_lo_c))
        min_cands.append((t_hi_c, t_hi_c))
        a1 = w1_min * f05 / ((1.0 - f05) ** 3)
        a4 = w4_min * (1.0 - f95) / (f95 ** 3)
        qa = (w2_min + w3_min) - a1 * (f05 ** 2) - a4 * ((1.0 - f95) ** 2)
        qb = 2.0 * a1 * p05_f * f05 + 2.0 * a4 * p95_f * (1.0 - f95)
        qc = -a1 * (p05_f ** 2) - a4 * (p95_f ** 2)
        roots: List[float] = []
        if abs(qa) <= 1e-15 * max(1.0, abs(qb), abs(qc)):
            if abs(qb) > 1e-15:
                roots.append(-qc / qb)
        else:
            disc = qb * qb - 4.0 * qa * qc
            if disc >= 0.0:
                sd = math.sqrt(disc)
                q_term = -0.5 * (qb + math.copysign(sd, qb if qb != 0.0 else 1.0))
                if abs(q_term) > 1e-30:
                    roots.append(q_term / qa)
                    roots.append(qc / q_term)
                else:
                    roots.append(-qb / (2.0 * qa))
        for rt in roots:
            if t_lo_c - tol <= rt <= t_hi_c + tol:
                rtc = min(max(rt, t_lo_c), t_hi_c)
                min_cands.append((rtc, rtc))

    s_min = min(_eval_min(u_val, v_val) for u_val, v_val in min_cands)
    return {
        "sample_n": n,
        "feasible": True,
        "d_rms3_lb": float((s_min / float(n)) ** (1.0 / 3.0)),
        "d_rms3_ub": float((s_max / float(n)) ** (1.0 / 3.0)),
    }


def _type7_quantile(vals: Sequence[float], p: float) -> float:
    """Exact pure-Python Hyndman-Fan Type 7 linear quantile (matches numpy.quantile(..., method='linear'))."""
    s = sorted(float(x) for x in vals)
    n = len(s)
    assert n >= 1
    if n == 1:
        return s[0]
    idx = p * (n - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    frac = idx - lo
    return s[lo] + frac * (s[hi] - s[lo])


def verify_type7_finite_sample_cubic_mean_bounds() -> Dict[str, Any]:
    """Verify exact joint-constraint finite-sample Hyndman-Fan Type-7 cubic-mean bounds over N = 2..500 including N=18 (17x10+1x300) and N=200 (58/65/65/68)."""
    checked_cases = 0
    for n in range(2, 501):
        # Stress extreme right-skew (n-1 at 10, 1 at 300), left-skew (1 at 10, n-1 at 300),
        # two-tail extreme, and smooth ramp samples across every finite N = 2..500
        samples = [
            [10.0] * (n - 1) + [300.0],
            [10.0] + [300.0] * (n - 1),
            [10.0] + [65.0] * max(0, n - 2) + [300.0],
            [58.0 + 12.0 * (i / (n - 1)) for i in range(n)],
        ]
        for arr in samples:
            d_min = min(arr)
            p05 = _type7_quantile(arr, 0.05)
            p95 = _type7_quantile(arr, 0.95)
            d_max = max(arr)
            d_rms3 = (sum(x ** 3 for x in arr) / n) ** (1.0 / 3.0)
            b = compute_type7_finite_sample_cubic_mean_bounds_iip0033(
                d_min=d_min, p05=p05, p95=p95, d_max=d_max, sample_n=n
            )
            assert b["feasible"] is True
            assert b["d_rms3_lb"] - 1e-9 <= d_rms3 <= b["d_rms3_ub"] + 1e-9, (
                f"Finite-sample Type-7 bound failed at N={n}: {b['d_rms3_lb']} <= {d_rms3} <= {b['d_rms3_ub']}"
            )
            checked_cases += 1

    # Explicit N = 18 counterexample (17x10 + 1x300): exact joint constraints force unique vector (lb == ub == 114.4954 um)
    arr_18 = [10.0] * 17 + [300.0]
    p05_18 = _type7_quantile(arr_18, 0.05)
    p95_18 = _type7_quantile(arr_18, 0.95)
    d_rms3_18 = (sum(x ** 3 for x in arr_18) / 18.0) ** (1.0 / 3.0)
    b_18 = compute_type7_finite_sample_cubic_mean_bounds_iip0033(
        d_min=10.0, p05=p05_18, p95=p95_18, d_max=300.0, sample_n=18
    )
    asymp_ub_18 = (0.05 * (p05_18 ** 3) + 0.90 * (p95_18 ** 3) + 0.05 * (300.0 ** 3)) ** (1.0 / 3.0)
    assert abs(p05_18 - 10.0) < 1e-12 and abs(p95_18 - 53.5) < 1e-12
    assert d_rms3_18 > asymp_ub_18  # proves why N-aware bound is required for small finite N
    assert b_18["feasible"] is True
    assert abs(b_18["d_rms3_lb"] - d_rms3_18) < 1e-9 and abs(b_18["d_rms3_ub"] - d_rms3_18) < 1e-9

    # Explicit N = 200 pinned-extremum counterexample (min=p05=p95=58, max=70, claimed d_rms3=58.001 < exact lb=58.0732)
    b_200_pin = compute_type7_finite_sample_cubic_mean_bounds_iip0033(
        d_min=58.0, p05=58.0, p95=58.0, d_max=70.0, sample_n=200
    )
    exact_200_lb = ((199.0 * (58.0 ** 3) + 1.0 * (70.0 ** 3)) / 200.0) ** (1.0 / 3.0)
    assert b_200_pin["feasible"] is True
    assert abs(b_200_pin["d_rms3_lb"] - exact_200_lb) < 1e-9
    assert b_200_pin["d_rms3_lb"] > 58.001

    # Explicit N = 200 joint-constraint counterexample (min=58, p05=p95=65, max=68 -> x[9..190]=65, exact lb=64.7323, ub=65.1098)
    b_200_joint = compute_type7_finite_sample_cubic_mean_bounds_iip0033(
        d_min=58.0, p05=65.0, p95=65.0, d_max=68.0, sample_n=200
    )
    exact_200_joint_lb = ((9.0 * (58.0 ** 3) + 190.0 * (65.0 ** 3) + (68.0 ** 3)) / 200.0) ** (1.0 / 3.0)
    exact_200_joint_ub = (((58.0 ** 3) + 190.0 * (65.0 ** 3) + 9.0 * (68.0 ** 3)) / 200.0) ** (1.0 / 3.0)
    assert b_200_joint["feasible"] is True
    assert abs(b_200_joint["d_rms3_lb"] - exact_200_joint_lb) < 1e-9
    assert abs(b_200_joint["d_rms3_ub"] - exact_200_joint_ub) < 1e-9
    assert 64.71 < b_200_joint["d_rms3_lb"] - 1e-9 and b_200_joint["d_rms3_lb"] <= 64.74 <= b_200_joint["d_rms3_ub"]

    return {
        "n_range_verified": [2, 500],
        "stress_samples_verified_total": checked_cases,
        "all_finite_sample_bounds_hold": True,
        "n18_counterexample_17x10_1x300": {
            "sample_n": 18,
            "d_min_um": 10.0,
            "type7_p05_um": round(p05_18, 4),
            "type7_p95_um": round(p95_18, 4),
            "d_max_um": 300.0,
            "true_d_rms3_um": round(d_rms3_18, 4),
            "asymptotic_0p05_ub_um": round(asymp_ub_18, 4),
            "exact_finite_sample_lb_um": round(b_18["d_rms3_lb"], 4),
            "exact_finite_sample_ub_um": round(b_18["d_rms3_ub"], 4),
            "finite_sample_bound_holds": True,
        },
        "n200_pinned_extremum_counterexample_58_58_58_70": {
            "sample_n": 200,
            "d_min_um": 58.0,
            "type7_p05_um": 58.0,
            "type7_p95_um": 58.0,
            "d_max_um": 70.0,
            "claimed_impossible_d_rms3_um": 58.001,
            "exact_pinned_extremum_lb_um": round(b_200_pin["d_rms3_lb"], 4),
            "rejected_by_pinned_extremum_lb": bool(58.001 < b_200_pin["d_rms3_lb"] - 1e-9),
        },
        "n200_joint_constraint_counterexample_58_65_65_68": {
            "sample_n": 200,
            "d_min_um": 58.0,
            "type7_p05_um": 65.0,
            "type7_p95_um": 65.0,
            "d_max_um": 68.0,
            "claimed_impossible_d_rms3_um": 64.71,
            "accepted_feasible_d_rms3_um": 64.74,
            "exact_joint_lb_um": round(b_200_joint["d_rms3_lb"], 4),
            "exact_joint_ub_um": round(b_200_joint["d_rms3_ub"], 4),
            "rejected_64p71_by_exact_joint_lb": bool(64.71 < b_200_joint["d_rms3_lb"] - 1e-9),
            "accepted_64p74_within_exact_joint_bounds": bool(
                b_200_joint["d_rms3_lb"] - 1e-9 <= 64.74 <= b_200_joint["d_rms3_ub"] + 1e-9
            ),
        },
        "formula_note": (
            "For any finite sample of size N >= 2 with Hyndman-Fan Type-7 linear interpolation constraints "
            "(1-f05)*x[i05] + f05*x[i05+1] = p05 and (1-f95)*x[i95] + f95*x[i95+1] = p95 over sorted "
            "d_min = x[0] <= ... <= x[N-1] = d_max, compute_type7_finite_sample_cubic_mean_bounds_iip0033 "
            "checks exact feasibility and computes the exact global min/max of ((1/N)*sum(x_k^3))^(1/3) "
            "(0-variable exact evaluation for N <= 20 and N == 1 mod 20; exact 2D separable strictly-convex "
            "polytope minimization and vertex maximization for N >= 22 with fractional indices; 1e-9 tolerance)."
        ),
    }


def solve_design_z_supercell_fixed_pitch(
    *,
    label: str,
    top_diams_um: Sequence[float],
    bot_diams_um: Sequence[float],
    slab_h_um: float = 140.0,
    z_c_um: float = 27.5,
    staggered: bool = True,
    bed_po2_mmhg: float = 30.0,
    standoff_um: float = 15.0,
    standoff_d_factor: float = 1.0,
    d_gel_cm2_s: float = D_GEL_CHRON,
    k_perf_s_inv: float = K_PERF_CHRON,
    ocr_mult: float = 1.097,
    stromal_ocr_mult: float | None = None,
    nominal_d_um: float = 65.0,
    nominal_loading_ieq_cm2: float = 1500.0,
    nx_tot: int = 40,
    nz: int = 28,
) -> Dict[str, Any]:
    """Fixed-pitch (p = 104.1604 um for 65 um @ 1,500 IEQ/cm^2) 2x2 or 1x1 supercell 3D FVM."""
    n_side = int(round(math.sqrt(len(top_diams_um))))
    assert len(top_diams_um) == n_side * n_side and len(bot_diams_um) == n_side * n_side
    r_nom_cm = 0.5 * nominal_d_um * 1e-4
    v_nom = (4.0 / 3.0) * math.pi * (r_nom_cm ** 3)
    n_layer_nom = 0.5 * ((nominal_loading_ieq_cm2 * VOL_IEQ) / v_nom)
    pitch_cm = math.sqrt(1.0 / n_layer_nom)
    L_cm = n_side * pitch_cm
    slab_h_cm = slab_h_um * 1e-4
    z_c_cm = z_c_um * 1e-4

    dx = L_cm / nx_tot
    dz = slab_h_cm / nz
    xs = (np.arange(nx_tot) + 0.5) * dx
    zs = (np.arange(nz) + 0.5) * dz - 0.5 * slab_h_cm
    X, Y, Z = np.meshgrid(xs, xs, zs, indexing="ij")
    dv = dx * dx * dz

    def pdist(a: np.ndarray, c: float) -> np.ndarray:
        d = np.abs(a - c)
        return np.minimum(d, L_cm - d)

    in_any = np.zeros_like(X, dtype=bool)
    vmax = np.zeros_like(X, dtype=np.float64)
    sphere_masks: List[Tuple[str, float, np.ndarray]] = []
    tot_vol = 0.0
    all_diams = list(top_diams_um) + list(bot_diams_um)

    for layer_idx, d_list in enumerate((top_diams_um, bot_diams_um)):
        z_sign = +1.0 if layer_idx == 0 else -1.0
        for idx, d_um in enumerate(d_list):
            ix, iy = divmod(idx, n_side)
            if layer_idx == 1 and staggered:
                cx = ((ix + 1.0) * pitch_cm) % L_cm
                cy = ((iy + 1.0) * pitch_cm) % L_cm
            else:
                cx = (ix + 0.5) * pitch_cm
                cy = (iy + 0.5) * pitch_cm
            r_cm = 0.5 * d_um * 1e-4
            v_true = (4.0 / 3.0) * math.pi * (r_cm ** 3)
            tot_vol += v_true
            msk = (pdist(X, cx) ** 2 + pdist(Y, cy) ** 2 + (Z - z_sign * z_c_cm) ** 2) <= (r_cm ** 2)
            assert not np.any(in_any & msk), f"Overlap in {label}"
            in_any |= msk
            v_vox = float(np.sum(msk)) * dv
            vmax[msk] = (VMAX_PURE * ocr_mult) * (v_true / v_vox)
            sphere_masks.append((f"{'top' if layer_idx == 0 else 'bot'}_{idx}_{d_um:.2f}um", d_um, msk))

    in_gel = ~in_any
    s_mult = ocr_mult if stromal_ocr_mult is None else stromal_ocr_mult
    t_phi_nom = (nominal_loading_ieq_cm2 * VOL_IEQ) / slab_h_cm
    vmax[in_gel] = (t_phi_nom * VMAX_PURE * 0.15 * s_mult) / float(np.mean(in_gel))
    perm = np.where(in_any, K_ISLET, d_gel_cm2_s * ALPHA_LOW_SALT)
    kp = np.where(in_gel, k_perf_s_inv, 0.0)

    res = l4._solve_3d_oxygen_strict_iip_0019(
        perm, vmax, np.zeros_like(X), kp, dx, dx, dz,
        p_bed=bed_po2_mmhg, standoff_z_cm=standoff_um * 1e-4, p_ext=D_GEL_CHRON * standoff_d_factor * ALPHA_LOW_SALT,
    )
    p3 = res["p_3d"]
    worst_face, worst_core = 1e9, 1e9
    for _, _, msk in sphere_masks:
        fv = l4._compute_exact_fv_interface_stats_iip_0019(p3, perm, msk)
        worst_face = min(worst_face, float(fv["fv_harmonic_face_min_mmhg"]))
        worst_core = min(worst_core, float(np.min(p3[msk])))

    d_rms3 = (sum(d ** 3 for d in all_diams) / len(all_diams)) ** (1.0 / 3.0)
    achieved_load = (tot_vol / VOL_IEQ) / (L_cm * L_cm)
    min_wall_clear_um = (0.5 * slab_h_um) - z_c_um - 0.5 * max(all_diams)
    return {
        "label": label,
        "staggered_lower_layer": bool(staggered),
        "fixed_pitch_um": round(pitch_cm * 1e4, 4),
        "nx_tot": nx_tot,
        "nz": nz,
        "cell_site_count": len(all_diams),
        "d_rms3_um": round(d_rms3, 3),
        "cell_type7_p95_um": round(_type7_quantile(all_diams, 0.95), 4),
        "d_max_um": round(max(all_diams), 2),
        "min_wall_clearance_um": round(min_wall_clear_um, 3),
        "achieved_loading_ieq_cm2": round(achieved_load, 2),
        "k_perf_s_inv": k_perf_s_inv,
        "worst_face_mmhg": round(worst_face, 4),
        "worst_core_mmhg": round(worst_core, 4),
        "paired_face_and_core_mmhg": [round(worst_face, 4), round(worst_core, 4)],
        "passes_15_and_5": bool(worst_face >= 15.0 and worst_core >= 5.0),
        "mass_balance_err": _stable_mass_balance_err(res["rel_mass_balance_error"]),
        "mass_balance_bound_lt_1e10": True,
    }


def solve_route1_supercell_fixed_pitch(
    *,
    label: str,
    top_diams_um: Sequence[float],
    bot_diams_um: Sequence[float],
    slab_um: float = 390.0,
    lam_um: float = 100.0,
    buf_um: float = 20.0,
    staggered: bool = False,
    bed_po2_mmhg: float = 35.0,
    ocr_mult: float = 1.097,
    standoff_phys_um: float = 15.0,
    nominal_d_um: float = 65.0,
    nominal_loading_ieq_cm2: float = 1500.0,
    nx_tot: int = 40,
    nz: int = 78,
) -> Dict[str, Any]:
    """Fixed-pitch (p = 104.1604 um) supercell 3D FVM for chronic vascularized Route 1 pouch (k_perf=0.1358 s^-1, zero active CaO2 source)."""
    n_side = int(round(math.sqrt(len(top_diams_um))))
    assert len(top_diams_um) == n_side * n_side and len(bot_diams_um) == n_side * n_side
    r_nom_cm = 0.5 * nominal_d_um * 1e-4
    v_nom = (4.0 / 3.0) * math.pi * (r_nom_cm ** 3)
    n_face_nom = 0.5 * ((nominal_loading_ieq_cm2 * VOL_IEQ) / v_nom)
    pitch_cm = math.sqrt(1.0 / n_face_nom)
    L_cm = n_side * pitch_cm
    slab_cm = slab_um * 1e-4
    lam_half_cm = 0.5 * lam_um * 1e-4
    band_inner_z = lam_half_cm + buf_um * 1e-4
    band_outer_z = 0.5 * slab_cm
    z_c_cm = 0.5 * (band_inner_z + band_outer_z)

    dx = L_cm / nx_tot
    dz = slab_cm / nz
    xs = (np.arange(nx_tot) + 0.5) * dx
    zs = (np.arange(nz) + 0.5) * dz - 0.5 * slab_cm
    X, Y, Z = np.meshgrid(xs, xs, zs, indexing="ij")
    dv = dx * dx * dz

    def pdist(a: np.ndarray, c: float) -> np.ndarray:
        d = np.abs(a - c)
        return np.minimum(d, L_cm - d)

    in_sph = np.zeros_like(X, dtype=bool)
    vmax = np.zeros_like(X, dtype=np.float64)
    sphere_masks: List[Tuple[str, float, np.ndarray]] = []
    tot_vol = 0.0
    all_diams = list(top_diams_um) + list(bot_diams_um)

    for z_sign, d_list in ((+1.0, top_diams_um), (-1.0, bot_diams_um)):
        for idx, d_um in enumerate(d_list):
            ix, iy = divmod(idx, n_side)
            if z_sign < 0.0 and staggered:
                cx = ((ix + 1.0) * pitch_cm) % L_cm
                cy = ((iy + 1.0) * pitch_cm) % L_cm
            else:
                cx = (ix + 0.5) * pitch_cm
                cy = (iy + 0.5) * pitch_cm
            r_cm = 0.5 * d_um * 1e-4
            v_true = (4.0 / 3.0) * math.pi * (r_cm ** 3)
            tot_vol += v_true
            msk = (pdist(X, cx) ** 2 + pdist(Y, cy) ** 2 + (Z - z_sign * z_c_cm) ** 2) <= (r_cm ** 2)
            assert not np.any(in_sph & msk), f"Overlap in {label}"
            in_sph |= msk
            v_vox = float(np.sum(msk)) * dv
            vmax[msk] = (VMAX_PURE * ocr_mult) * (v_true / v_vox)
            sphere_masks.append((f"{'top' if z_sign > 0 else 'bot'}_{idx}_{d_um:.2f}um", d_um, msk))

    in_lam = np.abs(Z) <= lam_half_cm
    in_buf = (np.abs(Z) > lam_half_cm) & (np.abs(Z) <= band_inner_z)
    perm = np.where(in_lam, K_GEL_CHRON * (0.42 ** 1.5), K_GEL_CHRON)
    perm[in_sph] = K_ISLET
    in_gel = (~in_sph) & (~in_lam) & (~in_buf)
    t_phi_nom = (nominal_loading_ieq_cm2 * VOL_IEQ) / slab_cm
    vmax[in_gel] = (t_phi_nom * VMAX_PURE * 0.15 * ocr_mult) / float(np.mean(in_gel))
    kp = np.where(in_gel, K_PERF_CHRON, 0.0)

    res = l4._solve_3d_oxygen_strict_iip_0019(
        perm, vmax, np.zeros_like(X), kp, dx, dx, dz,
        p_bed=bed_po2_mmhg, standoff_z_cm=standoff_phys_um * 1e-4, p_ext=K_GEL_CHRON / 3.0,
    )
    p3 = res["p_3d"]
    worst_face, worst_core = 1e9, 1e9
    for _, _, msk in sphere_masks:
        fv = l4._compute_exact_fv_interface_stats_iip_0019(p3, perm, msk)
        worst_face = min(worst_face, float(fv["fv_harmonic_face_min_mmhg"]))
        worst_core = min(worst_core, float(np.min(p3[msk])))

    d_rms3 = (sum(d ** 3 for d in all_diams) / len(all_diams)) ** (1.0 / 3.0)
    achieved_load = (tot_vol / VOL_IEQ) / (L_cm * L_cm)
    stress_regime = (
        "NOMINAL_1P097X_OCR_15P0UM_STANDOFF"
        if abs(ocr_mult - 1.097) < 1e-6 and abs(standoff_phys_um - 15.0) < 1e-6
        else (
            "PURE_PLUS10PCT_OCR_15P0UM_STANDOFF"
            if abs(standoff_phys_um - 15.0) < 1e-6
            else "COUPLED_PLUS10PCT_OCR_AND_16P5UM_STANDOFF"
        )
    )
    return {
        "label": label,
        "model_regime": "CHRONIC_VASCULARIZED_KPERF_0P1358_ZERO_ACTIVE_CAO2",
        "stress_regime": stress_regime,
        "staggered_lower_layer": bool(staggered),
        "slab_um": slab_um,
        "fixed_pitch_um": round(pitch_cm * 1e4, 4),
        "nx_tot": nx_tot,
        "nz": nz,
        "cell_site_count": len(all_diams),
        "d_rms3_um": round(d_rms3, 3),
        "cell_type7_p95_um": round(_type7_quantile(all_diams, 0.95), 4),
        "d_max_um": round(max(all_diams), 2),
        "achieved_loading_ieq_cm2": round(achieved_load, 2),
        "ocr_mult": round(ocr_mult, 5),
        "standoff_phys_um": standoff_phys_um,
        "worst_face_mmhg": round(worst_face, 4),
        "worst_core_mmhg": round(worst_core, 4),
        "passes_15_and_5": bool(worst_face >= 15.0 and worst_core >= 5.0),
        "mass_balance_err": _stable_mass_balance_err(res["rel_mass_balance_error"]),
        "mass_balance_bound_lt_1e10": True,
    }


def solve_3d_preperf_microtrack_fvm(
    *,
    label: str,
    r_lumen_um: float,
    tissue_bridge_um: float,
    phi_lumen: float = 0.15,
    cluster_diam_um: float = 65.0,
    bed_po2_mmhg: float = 35.0,
    sleeve_um: float = 0.0,
    sleeve_d_factor: float = 1.0,
    k_perf_wall_s_inv: float = 25.0,
    k_perf_lumen_s_inv: float = 0.0,
    ocr_mult: float = 1.097,
    stromal_frac: float = 0.15,
    gsis_mult: float = 1.0,
    nx_cross: int = 40,
) -> Dict[str, Any]:
    """Exact 3D FVM of periodic parallel cylindrical pre-vascularized micro-tracks."""
    r_lum_cm = r_lumen_um * 1e-4
    r_cl_cm = 0.5 * cluster_diam_um * 1e-4
    v_cl = (4.0 / 3.0) * math.pi * (r_cl_cm ** 3)
    pitch_y_cm = (2.0 * r_lumen_um + tissue_bridge_um) * 1e-4
    h_z_cm = pitch_y_cm

    if r_lumen_um <= 165.0:
        r_off_cm = max(0.0, r_lum_cm - r_cl_cm - 5.0e-4)
        m_ring = 2 if r_lumen_um <= 105.0 else (4 if r_lumen_um <= 130.0 else (6 if r_lumen_um <= 145.0 else 8))
        n_cl_unit = 1 + m_ring
        l_x_cm = (n_cl_unit * v_cl) / (phi_lumen * math.pi * (r_lum_cm ** 2))
        nx_ax = max(28, int(round(nx_cross * (l_x_cm / pitch_y_cm))))
        ny = nx_cross
        nz = nx_cross
    else:
        n_cl_unit = 0
        l_x_cm = 0.0200
        nx_ax = 12
        ny = nx_cross
        nz = nx_cross

    dx = l_x_cm / nx_ax
    dy = pitch_y_cm / ny
    dz = h_z_cm / nz
    xs = (np.arange(nx_ax) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dy - 0.5 * pitch_y_cm
    zs = (np.arange(nz) + 0.5) * dz - 0.5 * h_z_cm
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    R_cyl = np.sqrt(Y ** 2 + Z ** 2)
    dv = dx * dy * dz

    in_lumen = R_cyl <= r_lum_cm
    sleeve_cm = sleeve_um * 1e-4
    in_sleeve = (R_cyl > r_lum_cm) & (R_cyl <= (r_lum_cm + sleeve_cm))
    in_wall = R_cyl > (r_lum_cm + sleeve_cm)

    def pdist_x(a: np.ndarray, c: float) -> np.ndarray:
        d = np.abs(a - c)
        return np.minimum(d, l_x_cm - d)

    d_gel_use = D_GEL_ACUTE if k_perf_lumen_s_inv == 0.0 else D_GEL_CHRON
    k_gel_use = d_gel_use * ALPHA_LOW_SALT
    k_sleeve_use = D_GEL_CHRON * sleeve_d_factor * ALPHA_LOW_SALT

    perm = np.where(in_wall, K_GEL_CHRON, k_gel_use)
    if sleeve_um > 0.0:
        perm[in_sleeve] = k_sleeve_use

    vmax = np.zeros_like(X, dtype=np.float64)
    sphere_masks: List[Tuple[str, np.ndarray]] = []
    in_any_cl = np.zeros_like(X, dtype=bool)

    if n_cl_unit > 0:
        centers = [("center_axis", 0.25 * l_x_cm, 0.0, 0.0)]
        for k in range(m_ring):
            ang = (2.0 * math.pi * k) / m_ring
            centers.append((f"ring_{k}", 0.75 * l_x_cm, r_off_cm * math.cos(ang), r_off_cm * math.sin(ang)))
        for tag, xc, yc, zc in centers:
            msk = (pdist_x(X, xc) ** 2 + (Y - yc) ** 2 + (Z - zc) ** 2) <= (r_cl_cm ** 2)
            assert not np.any(in_any_cl & msk), f"Strict non-overlap violation in {label} ({tag})"
            assert np.all(in_lumen[msk]), f"Sphere protruded outside lumen in {label} ({tag})"
            in_any_cl |= msk
            v_vox = float(np.sum(msk)) * dv
            vmax[msk] = (VMAX_PURE * ocr_mult * gsis_mult) * (v_cl / v_vox)
            sphere_masks.append((tag, msk))
        perm[in_any_cl] = K_ISLET
        in_lumen_fluid = in_lumen & (~in_any_cl)
        vol_lumen_true = math.pi * (r_lum_cm ** 2) * l_x_cm
        total_stroma_rate = (phi_lumen * vol_lumen_true) * (VMAX_PURE * stromal_frac * ocr_mult)
        vmax[in_lumen_fluid] = total_stroma_rate / (float(np.sum(in_lumen_fluid)) * dv)
    else:
        perm[in_lumen] = K_ISLET * phi_lumen + k_gel_use * (1.0 - phi_lumen)
        vmax[in_lumen] = VMAX_PURE * ocr_mult * (gsis_mult + stromal_frac) * phi_lumen

    kp = np.zeros_like(X, dtype=np.float64)
    kp[in_wall] = k_perf_wall_s_inv
    if k_perf_lumen_s_inv > 0.0:
        kp[in_lumen & (~in_any_cl) if n_cl_unit > 0 else in_lumen] = k_perf_lumen_s_inv

    res = l4._solve_3d_oxygen_strict_iip_0019(
        perm, vmax, np.zeros_like(X), kp, dx, dy, dz,
        p_bed=bed_po2_mmhg, standoff_z_cm=15.0e-4, p_ext=K_GEL_CHRON,
    )
    p3 = res["p_3d"]
    if n_cl_unit > 0:
        worst_face = min(
            float(l4._compute_exact_fv_interface_stats_iip_0019(p3, perm, msk)["fv_harmonic_face_min_mmhg"])
            for _, msk in sphere_masks
        )
        core_min = float(np.min(p3[in_any_cl]))
    else:
        worst_face = float(np.min(p3[in_lumen]))
        core_min = worst_face

    return {
        "label": label,
        "r_lumen_um": r_lumen_um,
        "nx_cross": nx_cross,
        "grid_shape": [nx_ax, ny, nz],
        "sleeve_um": sleeve_um,
        "stromal_frac": stromal_frac,
        "worst_face_mmhg": round(worst_face, 4),
        "core_min_mmhg": round(core_min, 4),
        "passes_15_and_5": bool(worst_face >= 15.0 and core_min >= 5.0),
        "mass_balance_err": _stable_mass_balance_err(res["rel_mass_balance_error"]),
        "mass_balance_bound_lt_1e10": True,
    }


def verify_respirometry_contract() -> Dict[str, Any]:
    """Verify Agilent 16.60 / 16.87 / 14.94 / 14.00 / 9.87 uL KRBH caps & Sweet 2008 flow-through."""
    p_5pct = 0.05 * (760.0 - 47.0)       # 35.6500 mmHg
    p_21pct = 0.2095 * (760.0 - 47.0)    # 149.3735 mmHg
    ocr_endo_1p097 = VOL_IEQ * VMAX_PURE * 1.097 * 60.0
    ocr_stroma_15 = 0.15 * ocr_endo_1p097
    ocr_total_1p097 = ocr_endo_1p097 + ocr_stroma_15  # 4.547919e-12 mol/(min*IEQ)
    ocr_whole_pmol = round(ocr_total_1p097 * 1e12, 4) # 4.5479 pmol/(min*IEQ)

    regime_specs = [
        ("agilent_fluxpak_16p60uL_krbh_1p27", 16.60, ALPHA_KRBH, "Agilent FluxPak 103518-100 published 16.60 uL + KRBH"),
        ("agilent_brochure_16p87uL_krbh_1p27", 16.87, ALPHA_KRBH, "Agilent brochure 5991-8657EN p.9 published 16.87 uL (101122-100) + KRBH"),
        ("governing_minus10pct_14p94uL_krbh_1p27", 14.94, ALPHA_KRBH, "Conservative -10% margin off 16.60 uL (14.94 uL) + KRBH — 9/16 IEQ (+0.61 mmHg 16-IEQ air-FCCP margin)"),
        ("screen_displaced_14p00uL_krbh_1p27", 14.00, ALPHA_KRBH, "Screen-displaced 14.00 uL sensitivity + KRBH — 9/15 IEQ"),
        ("agilent_xf24_faq_9p87uL_krbh_1p27", 9.87, ALPHA_KRBH, "Agilent XF24 support FAQ 9.87 uL volume + KRBH — 6/10 IEQ (9/16 IEQ predicts 20.96/105.84 mmHg FAIL)"),
        ("obsolete_v7_22p70uL_lowsalt_1p34", 22.70, ALPHA_LOW_SALT, "OBSOLETE Gerencser 2009 V7 22.70 uL (mouse-liver mitochondria n=4) + 1.34e-9 comparator (false-pass 25.56/121.34 mmHg)"),
        ("obsolete_v7_22p70uL_purewater_1p39", 22.70, ALPHA_PURE_WATER_37C, "OBSOLETE Gerencser 2009 V7 22.70 uL + 1.39e-9 pure water comparator (false-pass 25.92/122.35 mmHg)"),
    ]
    static_regimes: Dict[str, Any] = {}
    for key, vol_uL, alpha, desc in regime_specs:
        cap = compute_static_xfe24_101122_100_cap_iip0033(
            volume_uL=vol_uL,
            alpha_mol_cm3_mmhg=alpha,
            ocr_whole_prep_pmol_min_ieq=ocr_whole_pmol,
            read_min=3.0,
        )
        dp_exact = (ocr_whole_pmol * 1e-12 * 3.0) / ((vol_uL * 1e-3) * alpha)
        static_regimes[key] = {
            "description": desc,
            "volume_uL": vol_uL,
            "alpha_mol_cm3_mmhg": alpha,
            "dp_basal_3min_mmhg_per_ieq": cap["dp_basal_3min_mmhg_per_ieq"],
            "max_ieq_5pct_basal_ge25": cap["max_ieq_5pct_basal_ge25"],
            "max_ieq_5pct_gsis_1p5x_ge25": cap["max_ieq_5pct_gsis_1p5x_ge25"],
            "max_ieq_21pct_fccp_2p5x_ge120": cap["max_ieq_21pct_fccp_2p5x_ge120"],
            "pred_9_ieq_5pct_gsis1p5_mmhg": round(p_5pct - 9.0 * 1.5 * dp_exact, 2),
            "pred_16_ieq_21pct_fccp2p5_mmhg": round(p_21pct - 16.0 * 2.5 * dp_exact, 2),
            "pred_15_ieq_5pct_gsis1p5_mmhg": round(p_5pct - 15.0 * 1.5 * dp_exact, 2),
            "pred_25_ieq_21pct_fccp2p5_mmhg": round(p_21pct - 25.0 * 2.5 * dp_exact, 2),
        }

    # Exact assertions locking all published Agilent hardware specifications
    assert static_regimes["governing_minus10pct_14p94uL_krbh_1p27"]["max_ieq_5pct_gsis_1p5x_ge25"] == 9
    assert static_regimes["governing_minus10pct_14p94uL_krbh_1p27"]["max_ieq_21pct_fccp_2p5x_ge120"] == 16
    assert abs(static_regimes["governing_minus10pct_14p94uL_krbh_1p27"]["pred_16_ieq_21pct_fccp2p5_mmhg"] - 120.61) <= 0.01
    assert static_regimes["screen_displaced_14p00uL_krbh_1p27"]["max_ieq_21pct_fccp_2p5x_ge120"] == 15
    assert static_regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["max_ieq_5pct_gsis_1p5x_ge25"] == 6
    assert static_regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["max_ieq_21pct_fccp_2p5x_ge120"] == 10
    assert static_regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["pred_9_ieq_5pct_gsis1p5_mmhg"] == 20.96
    assert static_regimes["agilent_xf24_faq_9p87uL_krbh_1p27"]["pred_16_ieq_21pct_fccp2p5_mmhg"] == 105.84

    # Primary flow-through perifusion modeled on Sweet 2008 methodological precedent
    ft_lead = compute_sweet2008_flow_through_po2_out_iip0033(
        n_ieq=200.0, flow_uL_min=200.0, alpha_mol_cm3_mmhg=ALPHA_KRBH, ocr_whole_prep_pmol_min_ieq=ocr_whole_pmol
    )
    assert ft_lead["p_out_5pct_basal_1p0x_mmhg"] == 32.07
    assert ft_lead["p_out_5pct_gsis_1p5x_mmhg"] == 30.28
    assert ft_lead["p_out_5pct_gsis_2p0x_mmhg"] == 28.49
    assert ft_lead["p_out_5pct_fccp_2p5x_mmhg"] == 26.70
    assert ft_lead["all_4_states_ge_25mmhg"] is True and ft_lead["passes_primary_spec"] is True

    return {
        "capture_screen_pore_size_um": 125.0,
        "requires_calibrated_volume_and_cluster_retention_for_static_caps": True,
        "sweet2008_methodological_precedent_note": (
            "Sweet et al. (2008 PMID:18021279; 750 handpicked islets, glucose/carbachol) establishes the flow-through "
            "differential-OCR methodological precedent; 200 IEQ @ 200 uL/min in 5% O2 KRBH and outlet values "
            "(32.07/30.28/28.49/26.70 mmHg) are modeled Stage-0 assumptions and predictions."
        ),
        "static_xfe24_101122_100_regimes": static_regimes,
        "flow_through_primary_lead_200ieq_200uLmin_krbh": ft_lead,
    }


def verify_fixed_pitch_bilayer_and_holds(meshes: Sequence[Tuple[int, int, int, int]]) -> Dict[str, Any]:
    """Run fixed-pitch 2x2 mixed-lot 65 um bilayer 3D FVM across meshes and prove all HOLD cases fail."""
    d_comp_68 = (2.0 * (65.0 ** 3) - (68.0 ** 3)) ** (1.0 / 3.0)
    d_comp_70 = (2.0 * (65.0 ** 3) - (70.0 ** 3)) ** (1.0 / 3.0)
    fixed_pitch_um = math.sqrt(1.0 / (0.5 * 1500.0 * VOL_IEQ / ((4.0 / 3.0) * math.pi * (32.5e-4 ** 3)))) * 1e4
    assert abs(fixed_pitch_um - 104.1604) < 1e-4

    # Verify GY-027 monolayer fixed-pitch vs re-derived 75 um packing gap
    mono_p65_um = math.sqrt(1.0 / (1500.0 * VOL_IEQ / ((4.0 / 3.0) * math.pi * (32.5e-4 ** 3)))) * 1e4
    mono_p75_um = math.sqrt(1.0 / (1500.0 * VOL_IEQ / ((4.0 / 3.0) * math.pi * (37.5e-4 ** 3)))) * 1e4
    gy027_packing_check = {
        "mono_65um_pitch_um": round(mono_p65_um, 4),
        "mono_65um_ordered_gap_um": round(mono_p65_um - 65.0, 2),
        "mono_75um_under_fixed_65um_pitch_gap_um": round(mono_p65_um - 75.0, 2),
        "mono_75um_rederived_pitch_um": round(mono_p75_um, 4),
        "mono_75um_rederived_gap_um": round(mono_p75_um - 75.0, 2),
    }
    assert gy027_packing_check["mono_75um_under_fixed_65um_pitch_gap_um"] == -1.35
    assert gy027_packing_check["mono_75um_rederived_gap_um"] == 16.29

    pass_runs: List[Dict[str, Any]] = []
    hold_runs: List[Dict[str, Any]] = []
    avascular_runs: List[Dict[str, Any]] = []

    for nx_t, nz_dz, nz_160, nz_r1 in meshes:
        r_z30 = solve_design_z_supercell_fixed_pitch(
            label=f"PASS_R2_140Stag_Chr30_NoFBR_kperf0p1358_p95_68_2xAdj_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            slab_h_um=140.0, z_c_um=27.5, staggered=True,
            bed_po2_mmhg=30.0, standoff_um=15.0, standoff_d_factor=1.0,
            d_gel_cm2_s=D_GEL_CHRON, k_perf_s_inv=K_PERF_CHRON, nx_tot=nx_t, nz=nz_dz,
        )
        assert r_z30["passes_15_and_5"]
        pass_runs.append(r_z30)

        r_z35 = solve_design_z_supercell_fixed_pitch(
            label=f"PASS_R2_140Stag_Chr35_15FBR3x_kperf0p1358_p95_68_2xAdj_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            slab_h_um=140.0, z_c_um=27.5, staggered=True,
            bed_po2_mmhg=35.0, standoff_um=15.0, standoff_d_factor=1.0 / 3.0,
            d_gel_cm2_s=D_GEL_CHRON, k_perf_s_inv=K_PERF_CHRON, nx_tot=nx_t, nz=nz_dz,
        )
        assert r_z35["passes_15_and_5"]
        pass_runs.append(r_z35)

        r_r1n = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_390_Nominal1p097_p95_68_2xAdj_nx{nx_t}",
            top_diams_um=[68.0, 68.0, d_comp_68, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, ocr_mult=1.097, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        assert r_r1n["passes_15_and_5"]
        pass_runs.append(r_r1n)

        r_r1s2 = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_390_CoupledPlus10OCR_16p5umStandoff_p95_68_2xAdj_nx{nx_t}",
            top_diams_um=[68.0, 68.0, d_comp_68, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        assert r_r1s2["passes_15_and_5"]
        pass_runs.append(r_r1s2)

        r_r1n_stag = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_390_Staggered_Nominal1p097_p95_68_2xAdj_nx{nx_t}",
            top_diams_um=[68.0, 68.0, d_comp_68, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=True, ocr_mult=1.097, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        assert r_r1n_stag["passes_15_and_5"]
        pass_runs.append(r_r1n_stag)

        r_r1s2_stag = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_390_Staggered_CoupledPlus10OCR_16p5umStandoff_p95_68_2xAdj_nx{nx_t}",
            top_diams_um=[68.0, 68.0, d_comp_68, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        assert r_r1s2_stag["passes_15_and_5"]
        pass_runs.append(r_r1s2_stag)

        r_350n = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_350_Nominal1p097_max68_t7p95_66p95_1x68_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=350.0, staggered=False, ocr_mult=1.097, standoff_phys_um=15.0, nx_tot=nx_t, nz=70,
        )
        assert r_350n["passes_15_and_5"]
        pass_runs.append(r_350n)

        r_350n_stag = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_350_Staggered_Nominal1p097_max68_t7p95_66p95_1x68_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=350.0, staggered=True, ocr_mult=1.097, standoff_phys_um=15.0, nx_tot=nx_t, nz=70,
        )
        assert r_350n_stag["passes_15_and_5"]
        pass_runs.append(r_350n_stag)

        r_350_pure = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_350_PurePlus10OCR_15umStandoff_max68_t7p95_66p95_1x68_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=350.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=70,
        )
        assert r_350_pure["passes_15_and_5"]
        pass_runs.append(r_350_pure)

        r_350_pure_stag = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_350_Staggered_PurePlus10OCR_15umStandoff_max68_t7p95_66p95_1x68_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=350.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=70,
        )
        assert r_350_pure_stag["passes_15_and_5"]
        pass_runs.append(r_350_pure_stag)

        r_p70_pure = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_390_PurePlus10OCR_15umStandoff_max70_t7p95_68p25_1x70_nx{nx_t}",
            top_diams_um=[70.0, 65.0, 65.0, d_comp_70],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        assert r_p70_pure["passes_15_and_5"]
        pass_runs.append(r_p70_pure)

        r_p70_pure_stag = solve_route1_supercell_fixed_pitch(
            label=f"PASS_R1_390_Staggered_PurePlus10OCR_15umStandoff_max70_t7p95_68p25_1x70_nx{nx_t}",
            top_diams_um=[70.0, 65.0, 65.0, d_comp_70],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        assert r_p70_pure_stag["passes_15_and_5"]
        pass_runs.append(r_p70_pure_stag)

        h_350s = solve_route1_supercell_fixed_pitch(
            label=f"HOLD_R1_350_CoupledPlus10OCR_16p5umStandoff_max68_t7p95_66p95_1x68_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=350.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=70,
        )
        h_350s_stag = solve_route1_supercell_fixed_pitch(
            label=f"HOLD_R1_350_Staggered_CoupledPlus10OCR_16p5umStandoff_max68_t7p95_66p95_1x68_nx{nx_t}",
            top_diams_um=[68.0, 65.0, 65.0, d_comp_68],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=350.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=70,
        )
        assert not h_350s["passes_15_and_5"] and not h_350s_stag["passes_15_and_5"]
        hold_runs.extend([h_350s, h_350s_stag])

        h_p70 = solve_route1_supercell_fixed_pitch(
            label=f"HOLD_R1_390_CoupledPlus10OCR_16p5umStandoff_max70_t7p95_68p25_1x70_nx{nx_t}",
            top_diams_um=[70.0, 65.0, 65.0, d_comp_70],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        h_p70_stag = solve_route1_supercell_fixed_pitch(
            label=f"HOLD_R1_390_Staggered_CoupledPlus10OCR_16p5umStandoff_max70_t7p95_68p25_1x70_nx{nx_t}",
            top_diams_um=[70.0, 65.0, 65.0, d_comp_70],
            bot_diams_um=[65.0, 65.0, 65.0, 65.0],
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        assert not h_p70["passes_15_and_5"] and not h_p70_stag["passes_15_and_5"]
        hold_runs.extend([h_p70, h_p70_stag])

        h_160e = solve_design_z_supercell_fixed_pitch(
            label=f"HOLD_160Aligned_70um_Exact1500_Chr35_15FBR3x_nx{nx_t}",
            top_diams_um=[70.0], bot_diams_um=[70.0],
            slab_h_um=160.0, z_c_um=37.5, staggered=False,
            bed_po2_mmhg=35.0, standoff_um=15.0, standoff_d_factor=1.0 / 3.0,
            d_gel_cm2_s=D_GEL_CHRON, k_perf_s_inv=K_PERF_CHRON, nominal_d_um=70.0,
            nx_tot=nx_t, nz=nz_160,
        )
        assert not h_160e["passes_15_and_5"]
        hold_runs.append(h_160e)

        h_160m = solve_design_z_supercell_fixed_pitch(
            label=f"HOLD_160Aligned_70um_Fixed104Mold_Chr35_15FBR3x_nx{nx_t}",
            top_diams_um=[70.0], bot_diams_um=[70.0],
            slab_h_um=160.0, z_c_um=37.5, staggered=False,
            bed_po2_mmhg=35.0, standoff_um=15.0, standoff_d_factor=1.0 / 3.0,
            d_gel_cm2_s=D_GEL_CHRON, k_perf_s_inv=K_PERF_CHRON, nominal_d_um=65.0,
            nx_tot=nx_t, nz=nz_160,
        )
        assert not h_160m["passes_15_and_5"]
        hold_runs.append(h_160m)

        h_160c = solve_design_z_supercell_fixed_pitch(
            label=f"HOLD_160Aligned_65um_At3p1464Cap_Chr30_15FBR3x_nx{nx_t}",
            top_diams_um=[65.0], bot_diams_um=[65.0],
            slab_h_um=160.0, z_c_um=37.5, staggered=False,
            bed_po2_mmhg=30.0, standoff_um=15.0, standoff_d_factor=1.0 / 3.0,
            d_gel_cm2_s=D_GEL_CHRON, k_perf_s_inv=K_PERF_CHRON,
            ocr_mult=3.1464 / (VOL_IEQ * VMAX_PURE * 60.0 * 1e12), stromal_ocr_mult=1.097,
            nx_tot=nx_t, nz=nz_160,
        )
        assert not h_160c["passes_15_and_5"]
        hold_runs.append(h_160c)

    # Mold overload without d_rms3 <= 65.0 um cap (68 um uniform = 1,717.4 IEQ/cm^2)
    h_ov_z = solve_design_z_supercell_fixed_pitch(
        label="HOLD_MoldOverload_drms3_68um_R2_140Stag_Chr30_nx60",
        top_diams_um=[68.0], bot_diams_um=[68.0], slab_h_um=140.0, z_c_um=27.5, staggered=True,
        bed_po2_mmhg=30.0, standoff_um=15.0, standoff_d_factor=1.0,
        d_gel_cm2_s=D_GEL_CHRON, k_perf_s_inv=K_PERF_CHRON, ocr_mult=1.097, nx_tot=60, nz=42,
    )
    h_ov_r1 = solve_route1_supercell_fixed_pitch(
        label="HOLD_MoldOverload_drms3_68um_R1_390_Plus10_nx60",
        top_diams_um=[68.0], bot_diams_um=[68.0], slab_um=390.0,
        ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=60, nz=78,
    )
    assert not h_ov_z["passes_15_and_5"] and not h_ov_r1["passes_15_and_5"]
    hold_runs.extend([h_ov_z, h_ov_r1])

    # Explicit Day-0 Avascular (k_perf = 0, D_gel_acute) separation (40 mmHg is a conservative Stage-0 policy floor; 38 mmHg passes)
    for lbl, bed, ocr_m in [
        ("Day0_Avascular_kperf0_140Stag_bed30_1p00x_nx60", 30.0, 1.00),
        ("Day0_Avascular_kperf0_140Stag_bed38_1p097x_nx60", 38.0, 1.097),
        ("Day0_Avascular_kperf0_140Stag_bed40_1p00x_nx60", 40.0, 1.00),
        ("Day0_Avascular_kperf0_140Stag_bed40_1p097x_nx60", 40.0, 1.097),
    ]:
        avascular_runs.append(solve_design_z_supercell_fixed_pitch(
            label=lbl, top_diams_um=[65.0], bot_diams_um=[65.0],
            slab_h_um=140.0, z_c_um=27.5, staggered=True,
            bed_po2_mmhg=bed, standoff_um=15.0, standoff_d_factor=1.0,
            d_gel_cm2_s=D_GEL_ACUTE, k_perf_s_inv=0.0, ocr_mult=ocr_m, nx_tot=60, nz=42,
        ))
    assert not avascular_runs[0]["passes_15_and_5"]
    assert avascular_runs[1]["passes_15_and_5"] and avascular_runs[2]["passes_15_and_5"] and avascular_runs[3]["passes_15_and_5"]

    # Spatial counterexamples proving aggregate lot histograms (p05>=58, d_rms3<=65, p95<=68, d_max<=70, <=1500 IEQ/cm^2)
    # are NECESSARY-NOT-SUFFICIENT for Route 1 chronic vascularized +10% OCR qualification, mandating worst-location stress pO2 >= 15.0 mmHg
    # on BOTH Route 1 350 um and Route 1 390 um (plus direct loaded-coupon verification across routes).
    d_comp_68_uneven = (2.0 * (64.999 ** 3) - (68.0 ** 3)) ** (1.0 / 3.0)
    top_32 = [65.0] * 16
    top_32[0] = 70.0
    top_32[10] = 59.08328752
    bot_32 = [65.0] * 16

    spatial_counterexamples: List[Dict[str, Any]] = []
    redistributed_8site_controls: List[Dict[str, Any]] = []
    for nx_t, _, _, nz_r1 in meshes:
        # Uneven 2x2 lot: 4x68 top + 4x61.692 bottom (d_rms3=64.999, p95=d_max=68, 1499.93 IEQ/cm^2)
        c_un_aligned_pure = solve_route1_supercell_fixed_pitch(
            label=f"COUNTEREXAMPLE_Uneven2x2_R1_390_Aligned_PurePlus10_15p0um_nx{nx_t}",
            top_diams_um=[68.0] * 4, bot_diams_um=[d_comp_68_uneven] * 4,
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        c_un_coupled = solve_route1_supercell_fixed_pitch(
            label=f"COUNTEREXAMPLE_Uneven2x2_R1_390_Aligned_CoupledPlus10_16p5um_nx{nx_t}",
            top_diams_um=[68.0] * 4, bot_diams_um=[d_comp_68_uneven] * 4,
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        c_un_stag_pure = solve_route1_supercell_fixed_pitch(
            label=f"COUNTEREXAMPLE_Uneven2x2_R1_390_Staggered_PurePlus10_15p0um_nx{nx_t}",
            top_diams_um=[68.0] * 4, bot_diams_um=[d_comp_68_uneven] * 4,
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        c_un_stag_coupled = solve_route1_supercell_fixed_pitch(
            label=f"COUNTEREXAMPLE_Uneven2x2_R1_390_Staggered_CoupledPlus10_16p5um_nx{nx_t}",
            top_diams_um=[68.0] * 4, bot_diams_um=[d_comp_68_uneven] * 4,
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        assert not c_un_aligned_pure["passes_15_and_5"]
        assert not c_un_coupled["passes_15_and_5"]
        assert not c_un_stag_pure["passes_15_and_5"]
        assert not c_un_stag_coupled["passes_15_and_5"]
        spatial_counterexamples.extend([c_un_aligned_pure, c_un_coupled, c_un_stag_pure, c_un_stag_coupled])

        # Identical-histogram redistributed 8-site controls (2x68 + 2x61.692 diagonally spaced per layer -> PASS >= 15.0 mmHg!)
        c_red_pure_aligned = solve_route1_supercell_fixed_pitch(
            label=f"CONTROL_Redistributed2x2_R1_390_Aligned_PurePlus10_15p0um_nx{nx_t}",
            top_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            bot_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        c_red_pure_stag = solve_route1_supercell_fixed_pitch(
            label=f"CONTROL_Redistributed2x2_R1_390_Staggered_PurePlus10_15p0um_nx{nx_t}",
            top_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            bot_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=nx_t, nz=nz_r1,
        )
        c_red_aligned = solve_route1_supercell_fixed_pitch(
            label=f"CONTROL_Redistributed2x2_R1_390_Aligned_CoupledPlus10_16p5um_nx{nx_t}",
            top_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            bot_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        c_red_stag = solve_route1_supercell_fixed_pitch(
            label=f"CONTROL_Redistributed2x2_R1_390_Staggered_CoupledPlus10_16p5um_nx{nx_t}",
            top_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            bot_diams_um=[68.0, d_comp_68_uneven, d_comp_68_uneven, 68.0],
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        assert c_red_pure_aligned["passes_15_and_5"] and c_red_pure_stag["passes_15_and_5"]
        assert c_red_aligned["passes_15_and_5"] and c_red_stag["passes_15_and_5"]
        redistributed_8site_controls.extend([c_red_pure_aligned, c_red_pure_stag, c_red_aligned, c_red_stag])

    # 32-site (4x4 per layer) supercell across staggered meshes + nx=80 aligned/pure checks:
    for nx_t, _, _, nz_r1 in meshes:
        c_32_stag_m = solve_route1_supercell_fixed_pitch(
            label=f"COUNTEREXAMPLE_Supercell32_R1_390_Staggered_CoupledPlus10_16p5um_nx{nx_t}",
            top_diams_um=top_32, bot_diams_um=bot_32,
            slab_um=390.0, staggered=True, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=nx_t, nz=nz_r1,
        )
        assert not c_32_stag_m["passes_15_and_5"]
        spatial_counterexamples.append(c_32_stag_m)

    c_32_aligned = solve_route1_supercell_fixed_pitch(
        label="COUNTEREXAMPLE_Supercell32_R1_390_Aligned_CoupledPlus10_16p5um_nx80",
        top_diams_um=top_32, bot_diams_um=bot_32,
        slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=16.5, nx_tot=80, nz=78,
    )
    c_32_pure = solve_route1_supercell_fixed_pitch(
        label="COUNTEREXAMPLE_Supercell32_R1_390_Aligned_PurePlus10_15p0um_nx80",
        top_diams_um=top_32, bot_diams_um=bot_32,
        slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=80, nz=78,
    )
    assert not c_32_aligned["passes_15_and_5"]
    assert not c_32_pure["passes_15_and_5"]
    spatial_counterexamples.extend([c_32_aligned, c_32_pure])

    # Axial-grid refinement audit for the marginal 32-site pure +10% OCR case (nz = 78, 117, 156 at nx = 80):
    c_32_pure_nz117 = solve_route1_supercell_fixed_pitch(
        label="AXIAL_REFINEMENT_Supercell32_R1_390_Aligned_PurePlus10_15p0um_nx80_nz117",
        top_diams_um=top_32, bot_diams_um=bot_32,
        slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=80, nz=117,
    )
    c_32_pure_nz156 = solve_route1_supercell_fixed_pitch(
        label="AXIAL_REFINEMENT_Supercell32_R1_390_Aligned_PurePlus10_15p0um_nx80_nz156",
        top_diams_um=top_32, bot_diams_um=bot_32,
        slab_um=390.0, staggered=False, ocr_mult=1.097 * 1.10, standoff_phys_um=15.0, nx_tot=80, nz=156,
    )
    assert c_32_pure_nz117["passes_15_and_5"] and c_32_pure_nz156["passes_15_and_5"]
    axial_refinement_32site_pure_plus10 = [c_32_pure, c_32_pure_nz117, c_32_pure_nz156]

    # Route 2 Design Z 140 um staggered check on both uneven counterexample distributions
    c_z_un = solve_design_z_supercell_fixed_pitch(
        label="CHECK_Uneven2x2_R2_140Stag_Chr30_NoFBR_nx40",
        top_diams_um=[68.0] * 4, bot_diams_um=[d_comp_68_uneven] * 4,
        slab_h_um=140.0, z_c_um=27.5, staggered=True,
        bed_po2_mmhg=30.0, standoff_um=15.0, standoff_d_factor=1.0, nx_tot=40, nz=28,
    )
    c_z_32 = solve_design_z_supercell_fixed_pitch(
        label="CHECK_Supercell32_R2_140Stag_Chr30_NoFBR_nx80",
        top_diams_um=top_32, bot_diams_um=bot_32,
        slab_h_um=140.0, z_c_um=27.5, staggered=True,
        bed_po2_mmhg=30.0, standoff_um=15.0, standoff_d_factor=1.0, nx_tot=80, nz=28,
    )
    assert c_z_un["passes_15_and_5"] and c_z_32["passes_15_and_5"]

    _seg_pure_nx40 = next(
        r for r in spatial_counterexamples if r["label"] == "COUNTEREXAMPLE_Uneven2x2_R1_390_Aligned_PurePlus10_15p0um_nx40"
    )
    _red_pure_nx40 = next(
        r for r in redistributed_8site_controls if r["label"] == "CONTROL_Redistributed2x2_R1_390_Aligned_PurePlus10_15p0um_nx40"
    )
    _seg_ref = next(
        (r for r in spatial_counterexamples if r["label"].startswith("COUNTEREXAMPLE_Uneven2x2_R1_390_Staggered_PurePlus10_15p0um_nx80")),
        spatial_counterexamples[2],
    )
    _red_ref = next(
        (r for r in redistributed_8site_controls if r["label"].endswith("_Staggered_PurePlus10_15p0um_nx80")),
        redistributed_8site_controls[1],
    )
    assert abs(_seg_pure_nx40["worst_face_mmhg"] - 14.2990) < 1e-3
    assert abs(_red_pure_nx40["worst_face_mmhg"] - 15.7151) < 1e-3
    _d0_30 = avascular_runs[0]
    _d0_38 = avascular_runs[1]
    _d0_40_1x = avascular_runs[2]
    _d0_40_1097 = avascular_runs[3]

    return {
        "fixed_pitch_um": round(fixed_pitch_um, 4),
        "nominal_loading_ieq_cm2": 1500.0,
        "route1_fvm_model_scope_note": (
            "Route 1 3D FVM runs model chronic vascularized hydrogel perfusion (k_perf = 0.1358 s^-1, zero active CaO2 source) "
            "across both aligned and staggered 2-layer registries on all 3 meshes, and explicitly separate pure +10% OCR (15.0 um standoff) "
            "from coupled +10% OCR + 16.5 um standoff. On 8-site supercells, one-68 and one-70 tails carry cell maxima of 68.0/70.0 um and "
            "exact 8-site Type-7 p95 of 66.95/68.25 um (both serialized per run). Aggregate lot histogram tolerances (p05 >= 58, d_rms3 <= 65, "
            f"p95 <= 68, d_max <= 70, <= 1,500 IEQ/cm^2) are necessary-not-sufficient: segregated 8-site placement "
            f"({_seg_pure_nx40['worst_face_mmhg']:.4f} vs computed redistributed control {_red_pure_nx40['worst_face_mmhg']:.4f} mmHg at nx=40, "
            f"and {_seg_ref['worst_face_mmhg']:.4f} vs {_red_ref['worst_face_mmhg']:.4f} mmHg staggered at nx={_seg_ref['nx_tot']}, all at 15.0 um standoff "
            f"and identical histogram) and coupled 32-site stress "
            f"({c_32_aligned['worst_face_mmhg']:.4f} mmHg) robustly fail < 15.0 mmHg on Route 1 390 um (while the marginal 32-site pure +10% OCR "
            f"case is axially mesh-sensitive: {c_32_pure['worst_face_mmhg']:.4f} / {c_32_pure_nz117['worst_face_mmhg']:.4f} / "
            f"{c_32_pure_nz156['worst_face_mmhg']:.4f} mmHg at nz = 78 / 117 / 156), mandating explicit worst-location stress pO2 >= 15.0 mmHg "
            "across BOTH Route 1 350 um and Route 1 390 um plus direct loaded-coupon verification."
        ),
        "day0_bed_40mmhg_policy_floor_note": (
            "40.0 mmHg Day-0 avascular bed pO2 is a conservative Stage-0 policy floor rather than the physical minimum: "
            f"at 40.0 mmHg bed pO2, paired (face, core) minima are {_d0_40_1x['worst_face_mmhg']:.2f} / {_d0_40_1x['worst_core_mmhg']:.2f} mmHg "
            f"at 1.00x OCR and {_d0_40_1097['worst_face_mmhg']:.2f} / {_d0_40_1097['worst_core_mmhg']:.2f} mmHg at 1.097x OCR; "
            f"at 1.097x demand, 38.0 mmHg bed pO2 passes at {_d0_38['worst_face_mmhg']:.4f} mmHg face / {_d0_38['worst_core_mmhg']:.4f} mmHg core, "
            f"whereas 30.0 mmHg fails at {_d0_30['worst_face_mmhg']:.2f} / {_d0_30['worst_core_mmhg']:.2f} mmHg."
        ),
        "gy027_monolayer_packing_check": gy027_packing_check,
        "pass_runs_multimesh": pass_runs,
        "hold_proofs_all_fail": hold_runs,
        "spatial_counterexamples_necessary_not_sufficient": {
            "route1_390um_failing_counterexamples": spatial_counterexamples,
            "route1_390um_redistributed_8site_controls": redistributed_8site_controls,
            "axial_refinement_sensitivity_32site_pure_plus10": axial_refinement_32site_pure_plus10,
            "axial_refinement_note": (
                "For 32-site Route 1 390 um pure +10% OCR at nx = 80, axial refinement nz = 78 -> 117 -> 156 shifts worst face "
                f"from {c_32_pure['worst_face_mmhg']:.4f} mmHg (below 15.0) to {c_32_pure_nz117['worst_face_mmhg']:.4f} and "
                f"{c_32_pure_nz156['worst_face_mmhg']:.4f} mmHg (above 15.0). That single marginal case is mesh-sensitive; "
                f"the segregated 8-site placement counterexample ({_seg_pure_nx40['worst_face_mmhg']:.4f} vs computed redistributed control "
                f"{_red_pure_nx40['worst_face_mmhg']:.4f} mmHg at nx = 40, 15.0 um standoff) and coupled +10% OCR + 16.5 um standoff stress "
                f"({c_32_aligned['worst_face_mmhg']:.4f} mmHg) carry the robust failure signal."
            ),
            "route2_140um_design_z_uneven_checks": [c_z_un, c_z_32],
        },
        "day0_avascular_split": avascular_runs,
    }


def verify_cylindrical_traps() -> Dict[str, Any]:
    """Re-verify GY-028 cylindrical pre-vascularized channel traps at 500,000 IEQ scale."""
    r500 = solve_3d_preperf_microtrack_fvm(
        label="GY028_Cylindrical_R500um_1mmID_Day0_Anoxia",
        r_lumen_um=500.0, tissue_bridge_um=500.0, phi_lumen=0.15,
        bed_po2_mmhg=35.0, k_perf_wall_s_inv=25.0, nx_cross=40,
    )
    assert r500["core_min_mmhg"] == 0.0 and r500["worst_face_mmhg"] == 0.0

    r120_15 = solve_3d_preperf_microtrack_fvm(
        label="GY028_Cylindrical_R120um_30umDenseThrombus_Stroma15",
        r_lumen_um=120.0, tissue_bridge_um=160.0, phi_lumen=0.15,
        bed_po2_mmhg=35.0, sleeve_um=30.0, sleeve_d_factor=0.5,
        k_perf_wall_s_inv=25.0, k_perf_lumen_s_inv=0.0, stromal_frac=0.15, nx_cross=40,
    )
    r120_25 = solve_3d_preperf_microtrack_fvm(
        label="GY028_Cylindrical_R120um_30umDenseThrombus_VPIOStroma25",
        r_lumen_um=120.0, tissue_bridge_um=160.0, phi_lumen=0.15,
        bed_po2_mmhg=35.0, sleeve_um=30.0, sleeve_d_factor=0.5,
        k_perf_wall_s_inv=25.0, k_perf_lumen_s_inv=0.0, stromal_frac=0.25, nx_cross=40,
    )
    assert abs(r120_15["core_min_mmhg"] - 1.4106) < 0.01
    assert abs(r120_25["core_min_mmhg"] - 0.5908) < 0.01

    mech = []
    for r_um, bridge_um, phi in [(100.0, 150.0, 0.15), (120.0, 160.0, 0.15)]:
        d_um = 2.0 * r_um
        r_m = r_um * 1e-6
        d_m = d_um * 1e-6
        v_lum_m3 = (VOL_500K_CM3 / phi) * 1e-6
        l_tot_m = v_lum_m3 / (math.pi * r_m ** 2)
        pitch_um = d_um + bridge_um
        area_cm2 = (l_tot_m * 100.0) * (pitch_um * 1e-4)
        pad_side_cm = math.sqrt(0.5 * area_cm2)
        n_tracks = 2 * int(round((pad_side_cm * 1e4) / pitch_um))
        l_single_m = pad_side_cm / 100.0
        f_pull_all_N = n_tracks * (math.pi * d_m * l_single_m * 2000.0)
        mech.append({
            "r_um": r_um,
            "d_um": d_um,
            "phi": phi,
            "total_track_length_m": round(l_tot_m, 1),
            "total_area_cm2": round(area_cm2, 1),
            "two_square_pads_side_cm": round(pad_side_cm, 2),
            "n_tracks_total": n_tracks,
            "aspect_ratio_L_over_D": round(l_single_m / d_m, 1),
            "pull_force_all_rods_2kPa_N": round(f_pull_all_N, 1),
            "vpio_jamming_d_over_D_100um": round(100.0 / d_um, 3),
            "vpio_jamming_d_over_D_150um": round(150.0 / d_um, 3),
        })

    return {
        "r500_anoxia_3d_fvm": r500,
        "r120_30um_clot_stroma15_3d_fvm": r120_15,
        "r120_30um_clot_vpio25_3d_fvm": r120_25,
        "r100_r120_mechanics_and_anatomy": mech,
    }


def verify_live_oq016_baseline() -> Dict[str, Any]:
    """Programmatically verify 2,160-stack Layer-1 parity vs committed Pareto JSON & pre-IIP-0033 subset and OQ-016 live tree."""
    oq016_spec = rank_mod._build_odds_composition_and_profile_scope_audit_iip_0003()
    ref_notes = oq016_spec["multiplied_conditional_totals"]["iip_0020_sensitivity_and_calibration_notes"]
    ref_tree = ref_notes["programmatic_2pow12_generator_audit_iip_0022"]
    lck = float(ref_tree["locked_2route_shared_ax_pouch_exact"])
    mod = float(ref_tree["modular_2route_shared_ax_pouch_exact"])

    # Also verify independent recomputation from reassessed_hurdles and live rank_islet1_stacks helper
    hurdles = oq016_spec["reassessed_hurdles"]
    pa0 = [float(h["one_design_design_a"]) for h in hurdles]
    pb0 = [float(h["one_design_design_b"]) for h in hurdles]
    px0 = [float(h["one_design_design_x"]) for h in hurdles]
    assert pa0 == ref_tree["recovered_marginals_check"]["design_a"]
    assert pb0 == ref_tree["recovered_marginals_check"]["design_b"]
    assert px0 == ref_tree["recovered_marginals_check"]["design_x"]
    c0 = [0.88, 0.87, 0.89, 0.85]
    recomputed_tree = rank_mod.compute_common_shock_2pow12_portfolio_tree(pa0, pb0, px0, c0)
    recomputed_lck = float(recomputed_tree["locked_2route_shared_ax_pouch_exact"])
    recomputed_mod = float(recomputed_tree["modular_2route_shared_ax_pouch_exact"])
    ceil = math.prod(c0)

    pareto_path = ROOT / "data" / "extracted" / "islet1_combinatorial_pareto_and_ladder_v0.8.0.json"
    pareto_doc = json.loads(pareto_path.read_text(encoding="utf-8"))
    committed_stacks = pareto_doc["evaluated_stacks"]

    _, studies_all, omics = rank_mod.load_inputs()
    iip0033_pmids = {"19555051", "25893782", "31527312", "36591328", "41536981", "42176117"}
    studies_pre = [s for s in studies_all if str(s.get("study_metadata", {}).get("pmid")) not in iip0033_pmids]
    combos_post = rank_mod.run_combinatorial_multi_gate_search(rank_mod.build_target_evidence_profiles(studies_all, omics))
    combos_pre = rank_mod.run_combinatorial_multi_gate_search(rank_mod.build_target_evidence_profiles(studies_pre, omics))

    assert len(combos_post) == len(combos_pre) == len(committed_stacks) == 2160
    mutated_vs_pre = sum(1 for a, b in zip(combos_post, combos_pre) if a != b)
    mutated_vs_committed = sum(
        1
        for a, c in zip(combos_post, committed_stacks)
        if a["combination_id"] != c["combination_id"]
        or a["global_rank_among_all_combinations"] != c["rank"]
        or round(float(a["composite_series_system_utility"]), 4) != float(c["composite_score"])
    )
    max_abs_score_delta = max(
        abs(float(a["composite_series_system_utility"]) - float(b["composite_series_system_utility"]))
        for a, b in zip(combos_post, combos_pre)
    )
    max_abs_committed_delta = max(
        abs(round(float(a["composite_series_system_utility"]), 4) - float(c["composite_score"]))
        for a, c in zip(combos_post, committed_stacks)
    )
    computed_delta_pp = round(
        100.0 * abs(recomputed_lck - lck)
        + 100.0 * abs(recomputed_mod - mod)
        + 100.0 * max_abs_score_delta
        + 100.0 * max_abs_committed_delta,
        6,
    )

    assert mutated_vs_pre == 0 and mutated_vs_committed == 0
    assert abs(100.0 * lck - 50.778198) < 1e-4
    assert abs(100.0 * mod - 54.116275) < 1e-4
    assert abs(100.0 * ceil - 57.917640) < 1e-4
    assert computed_delta_pp == 0.0
    return {
        "live_locked_pct": round(100.0 * lck, 6),
        "live_modular_pct": round(100.0 * mod, 6),
        "live_shared_ceiling_pct": round(100.0 * ceil, 6),
        "evaluated_stacks_checked": len(combos_post),
        "pre_iip0033_studies_count": len(studies_pre),
        "post_iip0033_studies_count": len(studies_all),
        "mutated_stacks_vs_pre_iip0033": mutated_vs_pre,
        "mutated_stacks_vs_committed_pareto": mutated_vs_committed,
        "max_abs_composite_score_delta": round(max_abs_score_delta, 12),
        "max_abs_committed_rounded_delta": round(max_abs_committed_delta, 12),
        "iip0033_live_headline_delta_pp": computed_delta_pp,
    }


def build_verification_artifact(*, fast: bool = False) -> Dict[str, Any]:
    meshes = [(40, 28, 32, 78)] if fast else [(40, 28, 32, 78), (60, 42, 48, 78), (80, 56, 64, 78)]
    return {
        "schema_version": "openislet_iip0033_bilayer_respirometry_verification_v0.8.0",
        "contribution_id": "IIP-0033",
        "decision_id": "DEC-112",
        "mesh_schedule": [list(m) for m in meshes],
        "respirometry_contract": verify_respirometry_contract(),
        "type7_finite_sample_cubic_mean_bounds": verify_type7_finite_sample_cubic_mean_bounds(),
        "fixed_pitch_bilayer_and_holds": verify_fixed_pitch_bilayer_and_holds(meshes),
        "gy028_cylindrical_trap_anchors": verify_cylindrical_traps(),
        "oq016_live_baseline_parity": verify_live_oq016_baseline(),
        "overall_status": "ALL_IIP0033_PHYSICS_AND_RESPIROMETRY_ASSERTIONS_PASSED",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Deterministic IIP-0033 3-mesh bilayer, respirometry & GY-028 verifier.")
    parser.add_argument("--write", action="store_true", help=f"Write verified artifact to {DEFAULT_OUT_PATH}")
    parser.add_argument("--check", action="store_true", help=f"Verify exact match against {DEFAULT_OUT_PATH}")
    parser.add_argument("--fast", action="store_true", help="Run single nx=40 mesh only (smoke check)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_PATH, help="Output JSON path")
    args = parser.parse_args()

    payload = build_verification_artifact(fast=args.fast)
    serialized = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    if args.write:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(serialized, encoding="utf-8")
        print(f"Wrote {args.out} (SHA-256: {digest}) — ALL ASSERTIONS PASSED.")
    elif args.check:
        if not args.out.is_file():
            raise SystemExit(f"Missing verification artifact: {args.out}")
        existing = args.out.read_text(encoding="utf-8")
        if existing != serialized:
            raise SystemExit(f"Artifact drift detected vs {args.out} (computed SHA-256: {digest})")
        print(f"Verified exact match against {args.out} (SHA-256: {digest}) — PASS.")
    else:
        print(f"IIP-0033 deterministic replay passed (SHA-256: {digest}; meshes={payload['mesh_schedule']}).")


if __name__ == "__main__":
    main()
