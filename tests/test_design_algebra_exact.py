#!/usr/bin/env python3
"""
OpenIslet (Project Islet-1) Exact-Algebra Gate for the HIP-Core wet-lab design.

Re-derives, from the run table alone and with `fractions.Fraction` only (no floating-point linear algebra, no import from
scripts/rank_islet1_stacks.py), every design-algebra figure published in data/protocols/hip_core_fractional_factorial_v0.6.2.json:

  * ranks of the 22-column genetic model [I, X1..X6, 15 two-factor interactions] and of the 23-column model with the additive stage block B
  * hat-matrix leverages (all runs, both models), unit-leverage rows, minimum rank after deleting any single run (brute force)
  * variance multipliers diag((X'X)^-1) for every main effect / 2FI / intercept / B, centered VIFs, main-effect block inflation
  * projection of B on the genetic model (SSR, coefficients c) and the partitioned-inverse identity Var_23(b_k) = Var_22(b_k) + c_k^2 / SSR
  * prediction variance at the lead (Scenario A) genotype, the base 34-run reference design and the optional 38-run 4-bridge upgrade

It also checks the run-table structure (generators, bridge/replicate/lead blocks, decoded genotypes) and that the headline figures
are quoted exactly in README.md and CLAUDE.md, with no stale unqualified 34-run claims in any repository document.

    python3 tests/test_design_algebra_exact.py          # CLI summary, non-zero exit on any failure
    python3 -m pytest tests/test_design_algebra_exact.py
"""

from __future__ import annotations

import json
import re
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL_PATH = ROOT / "data" / "protocols" / "hip_core_fractional_factorial_v0.6.2.json"
MANIFEST_PATH = ROOT / "manifest.json"

FACTORS = ("X1", "X2", "X3", "X4", "X5", "X6")
PAIRS = [(i, j) for i in range(6) for j in range(i + 1, 6)]
TERMS = ["I"] + list(FACTORS) + ["%s:%s" % (FACTORS[i], FACTORS[j]) for i, j in PAIRS] + ["B"]
LEAD_LEVELS = (1, 1, 1, -1, -1, -1)  # Scenario A lead: HLA-E KI, constitutive hEF1a-core CD47 package, Locus-3 swap; TYK2/TXNIP/RNLS wild type
LEGACY_LEVELS = (1, -1, 1, -1, -1, -1)  # legacy lean chassis (endogenous B2M-promoter bicistronic CD47)
LEAD_BLOCK = "Tier2_Lead_Chassis_A_Confirmatory_Replicate"
EXPECTED_BLOCKS = (
    (range(1, 17), "Tier1_16Run_Resolution_IV_Screening", 0),
    (range(17, 25), "Tier2_DeAliasing_Replicate_A", 1),
    (range(25, 27), "Tier2_InterStage_Bridge_Control", 1),
    (range(27, 35), "Tier2_DeAliasing_Replicate_B", 1),
    (range(35, 37), LEAD_BLOCK, 1),
)
# Numbers that were true only of the earlier hand-typed 34-run design / optional 36-run upgrade; they must not reappear unqualified in prose.
STALE_TOKENS = ("773/896", "0.862723", "24/29", "1.7274", "1.9459", "additive stage-block variance inflation")
QUALIFIER_RE = re.compile(r"base|unreplicated|never built|legacy|superseded|previous|earlier", re.I)
_PRIME = (1 << 61) - 1


# ---------------------------------------------------------------------------------------------------------------------------------
# exact linear algebra (independent re-implementation: max-magnitude pivoting, inverse verified against the Gram matrix)
# ---------------------------------------------------------------------------------------------------------------------------------
def model_rows(runs, with_block):
    rows = []
    for run in runs:
        lv = [int(run["coded_levels"][k]) for k in FACTORS]
        row = [1] + lv + [lv[i] * lv[j] for i, j in PAIRS]
        if with_block:
            row.append(int(run["stage_block"]))
        rows.append(row)
    return rows


def gram(rows):
    p = len(rows[0])
    return [[sum(r[a] * r[b] for r in rows) for b in range(p)] for a in range(p)]


def invert(matrix):
    n = len(matrix)
    a = [[Fraction(v) for v in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda i: abs(a[i][col]))
        if a[pivot][col] == 0:
            raise ValueError("singular matrix: the design cannot estimate every model term")
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [v / scale for v in a[col]]
        for i in range(n):
            if i != col and a[i][col] != 0:
                f = a[i][col]
                a[i] = [u - f * v for u, v in zip(a[i], a[col])]
    return [row[n:] for row in a]


def quad(row, inv):
    p = len(row)
    return sum(row[a] * inv[a][b] * row[b] for a in range(p) if row[a] for b in range(p) if row[b])


def rank_exact(rows):
    a = [[Fraction(v) for v in r] for r in rows]
    r = 0
    for c in range(len(a[0])):
        piv = next((i for i in range(r, len(a)) if a[i][c] != 0), None)
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        for i in range(r + 1, len(a)):
            if a[i][c] != 0:
                f = a[i][c] / a[r][c]
                a[i] = [u - f * v for u, v in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def rank_mod_p(rows):
    """Rank over GF(p); never exceeds the rational rank, so a full-rank answer is a rigorous certificate (used for the 36 row deletions)."""
    a = [[v % _PRIME for v in r] for r in rows]
    r = 0
    for c in range(len(a[0])):
        piv = next((i for i in range(r, len(a)) if a[i][c]), None)
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        inv = pow(a[r][c], -1, _PRIME)
        for i in range(r + 1, len(a)):
            if a[i][c]:
                f = a[i][c] * inv % _PRIME
                a[i] = [(u - f * v) % _PRIME for u, v in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def min_rank_after_single_deletion(rows):
    p = len(rows[0])
    lowest = p
    for i in range(len(rows)):
        sub = rows[:i] + rows[i + 1:]
        rk = rank_mod_p(sub)
        if rk < p:
            rk = rank_exact(sub)
        lowest = min(lowest, rk)
    return lowest


class Design(object):
    """Exact-rational summary of a run table."""

    def __init__(self, runs):
        self.runs = runs
        self.n = len(runs)
        self.x22 = model_rows(runs, False)
        self.x23 = model_rows(runs, True)
        self.inv22, self.inv23 = invert(gram(self.x22)), invert(gram(self.x23))
        for rows, inv in ((self.x22, self.inv22), (self.x23, self.inv23)):
            g, p = gram(rows), len(rows[0])
            if not all(sum(g[a][k] * inv[k][b] for k in range(p)) == int(a == b) for a in range(p) for b in range(p)):
                raise ArithmeticError("exact inverse failed its identity check")
        self.diag22 = [self.inv22[i][i] for i in range(22)]
        self.diag23 = [self.inv23[i][i] for i in range(23)]
        self.hat22 = [quad(r, self.inv22) for r in self.x22]
        self.hat23 = [quad(r, self.inv23) for r in self.x23]
        b = [r[22] for r in self.x23]
        xtb = [sum(r[j] * bi for r, bi in zip(self.x22, b)) for j in range(22)]
        self.proj = [sum(self.inv22[a][d] * xtb[d] for d in range(22)) for a in range(22)]
        self.ssr = sum(v * v for v in b) - sum(xtb[a] * self.proj[a] for a in range(22))
        self.mean_b = Fraction(sum(b), self.n)
        self.cross = [sum(r[j] * (bi - self.mean_b) for r, bi in zip(self.x23, b)) for j in range(1, 22)]
        self.vif23 = [None] + [self._vif(j) for j in range(1, 23)]
        self.rank22, self.rank23 = rank_exact(self.x22), rank_exact(self.x23)
        self.rank29 = rank_exact([r + [r[22] * r[k] for k in range(1, 7)] for r in self.x23])
        self.rank_tier1 = rank_exact([r for r, run in zip(self.x22, runs) if int(run["stage_block"]) == 0])

    def _vif(self, j):
        col = [Fraction(r[j]) for r in self.x23]
        mean = sum(col) / self.n
        return self.inv23[j][j] * sum((v - mean) ** 2 for v in col)

    def prediction_variance(self, levels, stage):
        row = [1] + list(levels) + [levels[i] * levels[j] for i, j in PAIRS] + [stage]
        return quad(row, self.inv23)


# ---------------------------------------------------------------------------------------------------------------------------------
# helpers for comparing published text with exact values
# ---------------------------------------------------------------------------------------------------------------------------------
def frac_of(text):
    """'175/202 (0.866337)' -> Fraction(175, 202); '1 (1.000000)' -> Fraction(1); '3/278' -> Fraction(3, 278)."""
    return Fraction(str(text).split()[0])


def plain(x):
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def pct_inflation(v22, v23):
    return round(float((v23 / v22 - 1) * 100), 4)


def load_protocol():
    return json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------------------------------------------------
# protocol verification
# ---------------------------------------------------------------------------------------------------------------------------------
def verify_protocol(proto):
    """Return a list of human-readable failures (empty list = every published figure is the exact value)."""
    problems = []

    def fail(msg):
        problems.append(msg)

    def same(label, published, exact):
        try:
            ok = frac_of(published) == exact
        except (ValueError, ZeroDivisionError, TypeError):
            ok = False
        if not ok:
            fail("%s: published %r != exact %s" % (label, published, plain(exact)))
            return
        m = re.search(r"\((-?\d+\.\d+)\)", str(published))
        if m and abs(float(m.group(1)) - float(exact)) > 5.1e-7:
            fail("%s: published decimal %s disagrees with exact %s = %.6f" % (label, m.group(1), plain(exact), float(exact)))

    def near(label, published, exact, tol=1e-6):
        if published is None or abs(float(published) - float(exact)) > tol:
            fail("%s: published %r != exact %.9f" % (label, published, float(exact)))

    runs = proto["runs"]
    n = len(runs)
    levels = [tuple(int(r["coded_levels"][k]) for k in FACTORS) for r in runs]

    execution = proto.get("execution_row_contract", {})
    if (execution.get("n_run_rows") != n or execution.get("n_unique_genotypes") != len(set(levels))
            or execution.get("execute_repeat_rows") is not True
            or execution.get("run_ids") != [row["run_id"] for row in runs]):
        fail("Execution contract must retain every coded row, including repeat/bridge rows used in the information matrix")
    alternate = proto.get("unevaluated_61_clone_alternative", {})
    if (alternate.get("status") != "UNEVALUATED_NOT_AUTHORIZED_BY_75_OR_108_CLONE_POWER"
            or any(alternate.get(key) is not None for key in ("covariance_matrix", "standard_errors", "minimum_detectable_effects"))):
        fail("61-clone shortcut must remain explicitly unevaluated unless its own covariance and power are derived")

    # ---- run-table structure -------------------------------------------------------------------------------------------------
    if [r["run_number"] for r in runs] != list(range(1, n + 1)):
        fail("run_number is not 1..%d in order" % n)
    if [r["run_id"] for r in runs] != ["HIP_FACTORIAL_RUN_%02d" % i for i in range(1, n + 1)]:
        fail("run_id is not HIP_FACTORIAL_RUN_01..%02d in order" % n)
    if any(v not in (-1, 1) for lv in levels for v in lv):
        fail("coded levels must be +/-1")
    if n != 36:
        fail("expected the 36-run design (16 screening + 8 + 2 + 8 + 2), found %d runs" % n)
    for rng, block, stage in EXPECTED_BLOCKS:
        for rn in rng:
            if rn <= n and (runs[rn - 1]["block"] != block or int(runs[rn - 1]["stage_block"]) != stage):
                fail("Run %02d should be block %s / stage %d, found %s / %s" % (rn, block, stage, runs[rn - 1]["block"], runs[rn - 1]["stage_block"]))
    tier1 = levels[:16]
    if len(set(tier1)) != 16:
        fail("Tier-1 screening runs are not 16 distinct genotypes")
    if any(lv[4] != lv[0] * lv[1] * lv[2] or lv[5] != lv[1] * lv[2] * lv[3] for lv in tier1):
        fail("Tier-1 runs violate the generators X5 = X1*X2*X3, X6 = X2*X3*X4")
    if any(sum(lv[k] for lv in tier1) != 0 for k in range(6)):
        fail("Tier-1 columns are not balanced")
    if n >= 36:
        if levels[24] != tier1[0] or levels[25] != tier1[15]:
            fail("Runs 25-26 must repeat Tier-1 Run 01 (all -1) and Run 16 (all +1)")
        if levels[26:34] != levels[16:24]:
            fail("Runs 27-34 must repeat the levels of Runs 17-24")
        if len(set(levels[16:24])) != 8:
            fail("Runs 17-24 are not 8 distinct de-aliasing genotypes")
        if levels[34:36] != [LEAD_LEVELS, LEAD_LEVELS]:
            fail("Runs 35-36 must both carry the Scenario A lead genotype %s" % (LEAD_LEVELS,))
    if len(set(levels)) != proto.get("n_distinct_genotypes"):
        fail("n_distinct_genotypes published %r != %d" % (proto.get("n_distinct_genotypes"), len(set(levels))))
    if proto.get("n_runs") != n or proto.get("n_runs_tier1_screening") != 16 or proto.get("n_runs_tier2_dealiasing") != n - 16:
        fail("n_runs / n_runs_tier1_screening / n_runs_tier2_dealiasing disagree with the run table")
    if proto.get("residual_degrees_of_freedom_single_rep_per_run") != n - 23:
        fail("residual_degrees_of_freedom_single_rep_per_run published %r != %d" % (proto.get("residual_degrees_of_freedom_single_rep_per_run"), n - 23))
    tier2_keys = [k for k in proto if re.fullmatch(r"tier2_dealiasing_augmentation_(\d+)_runs", k)]
    if len(tier2_keys) != 1 or len(proto[tier2_keys[0]]) != n - 16 or int(re.fullmatch(r"tier2_dealiasing_augmentation_(\d+)_runs", tier2_keys[0]).group(1)) != n - 16:
        fail("tier2_dealiasing_augmentation_<N>_runs key/length does not match the %d Tier-2 runs" % (n - 16))
    if len(proto.get("tier1_screening_16_runs", [])) != 16:
        fail("tier1_screening_16_runs must list 16 runs")
    if proto.get("tier1_screening_16_runs") != runs[:16]:
        fail("tier1_screening_16_runs does not match runs[:16]")
    if proto.get("tier2_dealiasing_augmentation_20_runs") != runs[16:]:
        fail("tier2_dealiasing_augmentation_20_runs does not match runs[16:]")
    expected_col_sums = {k: sum(int(r["coded_levels"][k]) for r in runs) for k in FACTORS}
    if proto.get("column_sums_across_36_runs_by_factor") != expected_col_sums:
        fail("column_sums_across_36_runs_by_factor %r != %r" % (proto.get("column_sums_across_36_runs_by_factor"), expected_col_sums))

    # decoded genotype strings are a pure function of the coded levels, and the lean-chassis pointers hit the right genotype
    seen = {}
    for r in runs:
        for key, text in r["decoded_genotype"].items():
            fk = key.split("_")[0]
            seen.setdefault((fk, int(r["coded_levels"][fk])), set()).add(text)
    for (fk, lv), texts in sorted(seen.items()):
        if len(texts) != 1:
            fail("decoded_genotype for %s at level %+d is not unique across runs: %s" % (fk, lv, sorted(texts)))
    # R02: 4-cell physical construct matrix for (X1, X2) must eliminate any orphan 2A tag when X1 = -1
    loc1_mat = proto.get("locus1_x1_by_x2_physical_construct_matrix", {})
    expected_loc1_ids = {
        (-1, -1): "LOC1_X1M_X2M_B2M_CD47WT_MONOCISTRONIC",
        (-1, 1): "LOC1_X1M_X2P_HEF1A_CD47GPI_MONOCISTRONIC",
        (1, -1): "LOC1_X1P_X2M_B2M_HLAE_TPT2A_CD47WT_BICISTRONIC",
        (1, 1): "LOC1_X1P_X2P_HEF1A_HLAE_TPT2A_CD47GPI_BICISTRONIC",
    }
    for (x1_lv, x2_lv), cid in expected_loc1_ids.items():
        cell_obj = loc1_mat.get(cid, {})
        if cell_obj.get("coded_X1") != x1_lv or cell_obj.get("coded_X2") != x2_lv:
            fail("locus1_x1_by_x2_physical_construct_matrix cell %s coded levels %r,%r != %d,%d" % (cid, cell_obj.get("coded_X1"), cell_obj.get("coded_X2"), x1_lv, x2_lv))
        if x1_lv == -1 and "NO 2A" not in str(cell_obj.get("cistron_valency", "")):
            fail("locus1_x1_by_x2_physical_construct_matrix cell %s must specify NO 2A peptide tag" % cid)
    for r in runs:
        x1_lv, x2_lv = int(r["coded_levels"]["X1"]), int(r["coded_levels"]["X2"])
        if r.get("locus1_physical_cassette_id") != expected_loc1_ids[(x1_lv, x2_lv)]:
            fail("%s locus1_physical_cassette_id %r != %r" % (r["run_id"], r.get("locus1_physical_cassette_id"), expected_loc1_ids[(x1_lv, x2_lv)]))
        if x1_lv == -1 and "P2A" in r["decoded_genotype"]["X2_CD47"]:
            fail("%s has X1=-1 (no HLA-E) but decoded_genotype.X2_CD47 still mentions an orphan P2A: %r" % (r["run_id"], r["decoded_genotype"]["X2_CD47"]))
    by_id = {r["run_id"]: r for r in runs}
    lead_pointer_ids = [proto.get("lean_chassis_a_exact_run_id"), proto.get("lean_chassis_a_replicate_run_id")]
    for rid in lead_pointer_ids:
        run = by_id.get(rid)
        if run is None or tuple(int(run["coded_levels"][k]) for k in FACTORS) != LEAD_LEVELS:
            fail("lean_chassis_a pointer %r does not reference a run at the lead genotype %s" % (rid, LEAD_LEVELS))
            continue
        g = run["decoded_genotype"]
        if not ("HLA-E" in g["X1_MHC_I"] and "Constitutive" in g["X2_CD47"] and "CXCL10::SERPINB9" in g["X3_Locus3"] and all("+/+" in g[k] for k in ("X4_TYK2", "X5_TXNIP", "X6_RNLS"))):
            fail("lean_chassis_a pointer %r decodes to %s, not the lean chassis (HLA-E KI + constitutive CD47 + Locus-3 swap, TYK2/TXNIP/RNLS +/+)" % (rid, g))
    if dict(proto.get("lean_chassis_a_coded_levels", {})) != dict(zip(FACTORS, LEAD_LEVELS)):
        fail("lean_chassis_a_coded_levels is not %s" % (LEAD_LEVELS,))
    for rid in proto.get("legacy_lean_chassis_endogenous_promoter_run_ids", []):
        if tuple(int(by_id[rid]["coded_levels"][k]) for k in FACTORS) != LEGACY_LEVELS:
            fail("legacy lean chassis pointer %s is not at %s" % (rid, LEGACY_LEVELS))

    # ---- exact algebra, 36-run design ---------------------------------------------------------------------------------------
    try:
        d = Design(runs)
    except (ValueError, ArithmeticError) as exc:
        return problems + ["exact algebra failed on the published run table: %s" % exc]
    if (d.rank22, d.rank23) != (22, 23):
        fail("rank(X) = %d/22 and rank([X|B]) = %d/23 (expected full rank)" % (d.rank22, d.rank23))
    if d.rank_tier1 != 14:
        fail("Tier-1 screening block rank %d != 14" % d.rank_tier1)
    if proto["alias_structure"].get("tier1_16run_rank") != d.rank_tier1:
        fail("alias_structure.tier1_16run_rank published %r != %d" % (proto["alias_structure"].get("tier1_16run_rank"), d.rank_tier1))
    if proto["alias_structure"].get("tier2_combined_all_main_and_2fi_estimable_full_rank_22_of_22") is not (d.rank22 == 22):
        fail("alias_structure full-rank-22 flag disagrees with the computed rank")
    if proto["alias_structure"].get("tier2_combined_with_stage_block_full_rank_23_of_23") is not (d.rank23 == 23):
        fail("alias_structure full-rank-23 flag disagrees with the computed rank")
    for chain in proto["alias_structure"].get("tier1_two_factor_interaction_alias_chains_resolved_by_tier2", []):
        if chain.get("dealiased_in_combined_design") is not (d.rank22 == 22):
            fail("%s dealiased flag disagrees with the computed rank" % chain.get("chain_id"))
    if not str(proto["model_matrix_rank_intercept_6main_15twofi"]).startswith("%d / 22" % d.rank22) or ("%d / 22 in Tier-1" % d.rank_tier1) not in proto["model_matrix_rank_intercept_6main_15twofi"]:
        fail("model_matrix_rank_intercept_6main_15twofi text disagrees with ranks %d/22 (Tier-1 %d/22)" % (d.rank22, d.rank_tier1))
    rank23_text = str(proto["model_matrix_rank_with_stage_block_indicator_23cols"])
    if not rank23_text.startswith("%d / 23" % d.rank23) or plain(d.ssr) not in rank23_text or ("rank %d / 29" % d.rank29) not in rank23_text:
        fail("model_matrix_rank_with_stage_block_indicator_23cols text must state rank %d/23, SSR %s and rank %d/29" % (d.rank23, plain(d.ssr), d.rank29))
    if ("mean(B) = %s" % plain(d.mean_b)) not in rank23_text:
        fail("rank text must quote mean(B) = %s" % plain(d.mean_b))
    for tag, group in (("23-col", d.hat23), ("22-col", d.hat22)):
        if sum(group) != (23 if tag == "23-col" else 22):
            fail("trace of the %s hat matrix is %s, not its rank" % (tag, sum(group)))
        if any(h <= 0 or h > 1 for h in group):
            fail("%s leverage outside (0, 1]" % tag)

    # leverage profiles (both models) and unit-leverage rows
    unit_rows = sorted({runs[i]["run_id"] for i, h in enumerate(d.hat23) if h == 1} | {runs[i]["run_id"] for i, h in enumerate(d.hat22) if h == 1})
    if sorted(proto["unit_leverage_augmentation_rows"]) != unit_rows:
        fail("unit_leverage_augmentation_rows %r != %r" % (proto["unit_leverage_augmentation_rows"], unit_rows))
    near("max_hat_matrix_leverage_across_all_runs", proto["max_hat_matrix_leverage_across_all_runs"], max(d.hat23), 1e-12)
    for key, hat, rows, with_block in (
        ("exact_rational_leverage_profile_with_stage_block_23cols", d.hat23, d.x23, True),
        ("exact_rational_leverage_profile_genetic_only_22cols", d.hat22, d.x22, False),
    ):
        prof = proto[key]
        same(key + ".max_leverage_rational", prof["max_leverage_rational"], max(hat))
        near(key + ".max_leverage_h_ii", prof["max_leverage_h_ii"], max(hat), 1e-12)
        exact_groups = {}
        for i, h in enumerate(hat):
            exact_groups.setdefault(h, []).append(runs[i]["run_number"])
        published_groups = {frac_of(g["h_ii_rational"]): sorted(g["runs"]) for g in prof["leverage_groups"]}
        if published_groups != exact_groups:
            fail("%s.leverage_groups differ from the exact leverage groups" % key)
        for g in prof["leverage_groups"]:
            same(key + ".group " + g["h_ii_rational"], g["h_ii"], frac_of(g["h_ii_rational"]))
            if g["n_runs"] != len(g["runs"]):
                fail("%s group %s n_runs != len(runs)" % (key, g["h_ii_rational"]))
        if prof["trace_sum_of_leverages"] != sum(hat):
            fail("%s.trace_sum_of_leverages %r != %s" % (key, prof["trace_sum_of_leverages"], sum(hat)))
        brute = min_rank_after_single_deletion(rows)
        if prof["min_rank_after_deleting_any_single_run"] != brute:
            fail("%s.min_rank_after_deleting_any_single_run published %r != brute-force %d" % (key, prof["min_rank_after_deleting_any_single_run"], brute))
        if (brute == len(rows[0])) != (max(hat) < 1):
            fail("%s: brute-force deletion rank contradicts the leverage criterion h_ii < 1" % key)
        if with_block:
            same(key + ".stage_block_projection_residual_ssr", prof["stage_block_projection_residual_ssr"], d.ssr)
            near(key + ".stage_block_projection_residual_l2_norm", prof["stage_block_projection_residual_l2_norm"], float(d.ssr) ** 0.5)

    # variance multipliers, VIFs, block inflation
    vv = proto["exact_rational_parameter_variance_and_vif_conditioning"]
    same("intercept var_22col", vv["intercept_var_multiplier"]["var_22col"], d.diag22[0])
    same("intercept var_23col", vv["intercept_var_multiplier"]["var_23col"], d.diag23[0])
    inflations = []
    for k, key in enumerate(FACTORS, start=1):
        pub = vv["main_effects_var_multiplier_by_factor"][key]
        same(key + " var_22col", pub["var_22col"], d.diag22[k])
        same(key + " var_23col", pub["var_23col"], d.diag23[k])
        infl = pct_inflation(d.diag22[k], d.diag23[k])
        inflations.append(infl)
        near(key + " block_inflation_pct", pub["block_inflation_pct"], infl, 1e-9)
        near(key + " centered_vif_23col", pub["centered_vif_23col"], d.vif23[k])
        # partitioned-inverse identity ties the independently inverted 22- and 23-column models together
        if d.diag23[k] != d.diag22[k] + d.proj[k] ** 2 / d.ssr:
            fail("%s violates Var_23 = Var_22 + c^2/SSR" % key)
    if d.inv23[22][22] != 1 / d.ssr:
        fail("Var(B) != 1/SSR in the 23-column model")
    for k in range(22):
        if d.diag23[k] != d.diag22[k] + d.proj[k] ** 2 / d.ssr:
            fail("%s violates Var_23 = Var_22 + c^2/SSR" % TERMS[k])
            break
    if vv["main_effect_block_inflation_pct_range"] != [min(inflations), max(inflations)]:
        fail("main_effect_block_inflation_pct_range %r != %r" % (vv["main_effect_block_inflation_pct_range"], [min(inflations), max(inflations)]))
    note = vv["algebraic_note_main_effect_block_inflation"]
    for k, key in enumerate(FACTORS):
        if ("%s %s%%" % (key, inflations[k])) not in note:
            fail("algebraic_note_main_effect_block_inflation must quote '%s %s%%'" % (key, inflations[k]))
    seen_pairs = []
    for tier in vv["two_factor_interactions_var_multiplier_tiers_23col"].values():
        same("2FI tier " + ",".join(tier["pairs"]), tier["var_23col"], d.diag23[TERMS.index(tier["pairs"][0])])
        for pair in tier["pairs"]:
            idx = TERMS.index(pair)
            seen_pairs.append(pair)
            if d.diag23[idx] != frac_of(tier["var_23col"]):
                fail("2FI %s is not at its tier multiplier %s" % (pair, tier["var_23col"]))
            same("2FI %s var_22col" % pair, tier["var_22col_by_pair"][pair], d.diag22[idx])
            near("2FI %s centered_vif_23col" % pair, tier["centered_vif_23col_by_pair"][pair], d.vif23[idx])
    if sorted(seen_pairs) != sorted(TERMS[7:22]):
        fail("2FI tiers must cover each of the 15 interactions exactly once")
    same("2FI range min", vv["two_factor_interaction_var_multiplier_range_23col"]["min"], min(d.diag23[7:22]))
    same("2FI range max", vv["two_factor_interaction_var_multiplier_range_23col"]["max"], max(d.diag23[7:22]))

    sb = vv["stage_block_B"]
    same("stage_block_B.var_23col_rational", sb["var_23col_rational"], d.diag23[22])
    near("stage_block_B.centered_vif_23col", sb["centered_vif_23col"], d.vif23[22])
    same("stage_block_B.projection_residual_ssr", sb["projection_residual_ssr"], d.ssr)
    same("stage_block_B.c_intercept", sb["projection_coefficient_c_of_B_on_intercept"], d.proj[0])
    for k, key in enumerate(FACTORS, start=1):
        same("stage_block_B.c[%s]" % key, sb["projection_coefficient_c_of_B_on_main_effects"][key], d.proj[k])
    if d.rank29 != 25:
        fail("rank([X_23 | X_stage_x_main]) = %d != 25 for 36-run design" % d.rank29)
    if sb["rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols"] != "%d / 29" % d.rank29:
        fail("rank with 6 stage x main interactions published %r != %d / 29" % (sb["rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols"], d.rank29))

    # lead-genotype prediction variance and the base 34-run reference design
    lead = vv["lead_chassis_a_prediction_variance"]
    if tuple(lead["coded_levels"]) != LEAD_LEVELS or lead["runs_with_this_genotype"] != [r["run_id"] for r, lv in zip(runs, levels) if lv == LEAD_LEVELS]:
        fail("lead_chassis_a_prediction_variance genotype / run list is wrong")
    same("lead prediction variance stage 0", lead["var_multiplier_stage0"], d.prediction_variance(LEAD_LEVELS, 0))
    same("lead prediction variance stage 1", lead["var_multiplier_stage1"], d.prediction_variance(LEAD_LEVELS, 1))
    base_runs = [r for r in runs if r["block"] != LEAD_BLOCK]
    base = Design(base_runs)
    base_ref = vv["base_34_run_design_reference"]
    if base_ref["n_runs"] != len(base_runs) or len(base_runs) != 34 or (base.rank22, base.rank23) != (22, 23):
        fail("base 34-run reference must be the 34 non-lead runs at full rank")
    same("base main-effect multiplier (23-col)", base_ref["main_effect_var_multiplier_23col"], base.diag23[1])
    base_infl = [pct_inflation(base.diag22[k], base.diag23[k]) for k in range(1, 7)]
    if base_ref["main_effect_block_inflation_pct_range"] != [min(base_infl), max(base_infl)] or any(v != 0 for v in base_infl):
        fail("base 34-run main-effect block inflation must be exactly 0.0%% (computed %r, published %r)" % (base_infl, base_ref["main_effect_block_inflation_pct_range"]))
    if any(base.proj[k] != 0 for k in range(1, 7)):
        fail("base 34-run projection of B must have zero coefficient on every main effect")
    same("base max 2FI multiplier", base_ref["max_2fi_var_multiplier_23col"], max(base.diag23[7:22]))
    same("base min 2FI multiplier", base_ref["min_2fi_var_multiplier_23col"], min(base.diag23[7:22]))
    same("base Var(B)", base_ref["stage_block_var_23col"], base.diag23[22])
    same("base max leverage (23-col)", base_ref["max_leverage_h_ii_23col"], max(base.hat23))
    same("base projection coefficient of B on the intercept", base_ref["projection_coefficient_c_of_B_on_intercept"], base.proj[0])
    for k, key in enumerate(FACTORS, start=1):
        same("base projection coefficient of B on %s" % key, base_ref["projection_coefficient_c_of_B_on_main_effects"][key], base.proj[k])
    same("base lead prediction variance stage 0", lead["base_34_run_design_for_comparison"]["var_multiplier_stage0"], base.prediction_variance(LEAD_LEVELS, 0))
    same("base lead prediction variance stage 1", lead["base_34_run_design_for_comparison"]["var_multiplier_stage1"], base.prediction_variance(LEAD_LEVELS, 1))

    # optional 38-run upgrade: repeat Tier-1 Runs 03 and 14 in stage 1 next to the Run 01 / Run 16 bridge repeats
    up = vv["optional_38run_4bridge_upgrade_specification"]
    extra = [{"coded_levels": runs[2]["coded_levels"], "stage_block": 1, "block": "Tier2_Optional_Bridge"}, {"coded_levels": runs[13]["coded_levels"], "stage_block": 1, "block": "Tier2_Optional_Bridge"}]
    up_design = Design(runs + extra)
    if up["n_runs_total"] != n + 2 or up["bridge_runs_repeated_from_tier1"] != ["HIP_FACTORIAL_RUN_01", "HIP_FACTORIAL_RUN_16", "HIP_FACTORIAL_RUN_03", "HIP_FACTORIAL_RUN_14"]:
        fail("optional upgrade run count / bridge list wrong")
    same("upgrade Var(B)", up["stage_block_var_multiplier"], up_design.diag23[22])
    same("upgrade SSR", up["stage_block_ssr"], up_design.ssr)
    near("upgrade centered VIF of B", up["stage_block_centered_vif"], up_design.vif23[22])
    same("upgrade max 2FI multiplier", up["max_2fi_var_multiplier_rational"], max(up_design.diag23[7:22]))
    near("upgrade max 2FI multiplier (float)", up["max_2fi_var_multiplier"], max(up_design.diag23[7:22]))
    near("upgrade max 2FI centered VIF", up["max_2fi_centered_vif"], max(up_design.vif23[7:22]))
    near("upgrade max 2FI reduction pct", up["max_2fi_reduction_vs_this_design_pct"], (1 - max(up_design.diag23[7:22]) / max(d.diag23[7:22])) * 100, 0.006)
    if up_design.rank29 != 26:
        fail("upgrade rank([X_23 | X_stage_x_main]) = %d != 26 for 38-run design" % up_design.rank29)
    if up["rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols"] != "%d / 29" % up_design.rank29:
        fail("upgrade rank with 6 stage x main interactions published %r != %d / 29" % (up["rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols"], up_design.rank29))
    for token in (plain(d.diag23[22]), plain(up_design.diag23[22]), plain(max(d.diag23[7:22])), plain(max(up_design.diag23[7:22]))):
        if token not in up["description"]:
            fail("optional upgrade description must quote %s" % token)

    # narrative fields quote the exact numbers
    if ("Max Leverage %s = %.6f" % (plain(max(d.hat23)), float(max(d.hat23)))) not in proto["design_name"] or ("%d Unit-Leverage Rows" % len(unit_rows)) not in proto["design_name"]:
        fail("design_name must quote 'Max Leverage %s = %.6f' and '%d Unit-Leverage Rows'" % (plain(max(d.hat23)), float(max(d.hat23)), len(unit_rows)))
    dealias = proto["alias_structure"]["dealiasing_note"]
    for token in ("max h_ii = %s" % plain(max(d.hat22)), "max h_ii = %s" % plain(max(d.hat23)), "SSR = %s" % plain(d.ssr)):
        if token not in dealias:
            fail("alias_structure.dealiasing_note must quote '%s'" % token)
    plan = proto["biological_replication_and_blocking_plan"]["randomization_and_blocking"]
    for token in (plain(max(d.hat23)), plain(d.ssr), "%d distinct genotypes" % len(set(levels))):
        if token not in plan:
            fail("biological_replication_and_blocking_plan must quote '%s'" % token)

    # R03 & R04: verify statistical_analysis_and_power_plan MDE and 108-independent vs 75-shared-clone covariance table
    mde = proto.get("statistical_analysis_and_power_plan", {}).get("minimum_detectable_effects_at_80pct_power", {})
    m1 = mde.get("mode1_108_independent_run_nested_clones_exact_iid_run_covariance", {})
    m1_23 = m1.get("two_stage_23col_model_with_additive_stage_block_B", {})
    m1_22 = m1.get("single_stage_22col_genetic_model", {})
    same("power 23col best 2FI multiplier", m1_23.get("two_factor_interaction_best_2pairs_115_over_1616_X4X5_X5X6", {}).get("var_multiplier_rational"), min(d.diag23[7:22]))
    same("power 23col worst 2FI multiplier", m1_23.get("two_factor_interaction_worst_1pair_1191_over_1616_X1X3", {}).get("var_multiplier_rational"), max(d.diag23[7:22]))
    same("power 22col best 2FI multiplier", m1_22.get("two_factor_interaction_best_6pairs", {}).get("var_multiplier_rational"), min(d.diag22[7:22]))
    same("power 22col worst 2FI multiplier", m1_22.get("two_factor_interaction_worst_2pairs_X1X5_X4X6", {}).get("var_multiplier_rational"), max(d.diag22[7:22]))
    sigma_run_sq = Fraction(177, 6)
    worst_se_beta_23 = (float(sigma_run_sq * max(d.diag23[7:22])) ** 0.5)
    near("power 23col worst 2FI SE(beta)", m1_23.get("two_factor_interaction_worst_1pair_1191_over_1616_X1X3", {}).get("half_effect_beta_se_pp"), round(worst_se_beta_23, 2), 0.01)
    near("power 23col worst 2FI SE(2*beta)", m1_23.get("two_factor_interaction_worst_1pair_1191_over_1616_X1X3", {}).get("full_contrast_2beta_se_pp"), round(2.0 * worst_se_beta_23, 2), 0.01)
    near("power 23col worst 2FI SE(4*beta)", m1_23.get("two_factor_interaction_worst_1pair_1191_over_1616_X1X3", {}).get("difference_of_differences_4beta_se_pp"), round(4.0 * worst_se_beta_23, 2), 0.01)
    m2 = mde.get("mode2_75_shared_genotype_clones_across_25_distinct_genotypes_gls_and_ols_sandwich", {})
    if "two_stage_23col_gls_lmm_se_pp" not in m2 or "two_stage_23col_ols_sandwich_se_pp" not in m2:
        fail("mode2 shared-genotype clone covariance summary missing GLS or OLS sandwich SE entries")
    else:
        # Exact rational verification of Mode 2 shared-genotype covariance Sigma_shared
        sigma_clone_sq = Fraction(64, 1)
        sigma_eps_sq = Fraction(49, 1)
        n_clones_m2 = 3
        n_batches_m2 = 2
        cov_shared = [[Fraction(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if levels[i] == levels[j]:
                    cov_shared[i][j] += sigma_clone_sq / n_clones_m2
                if i == j:
                    cov_shared[i][j] += sigma_eps_sq / (n_clones_m2 * n_batches_m2)
        cov_inv = invert(cov_shared)
        xt_sinv = [[sum(d.x23[r][a] * cov_inv[r][c] for r in range(n)) for c in range(n)] for a in range(23)]
        xt_sinv_x = [[sum(xt_sinv[a][r] * d.x23[r][b] for r in range(n)) for b in range(23)] for a in range(23)]
        gls_cov23 = invert(xt_sinv_x)
        gls_main_ses = [float(gls_cov23[k][k]) ** 0.5 for k in range(1, 7)]
        gls_2fi_ses = [float(gls_cov23[k][k]) ** 0.5 for k in range(7, 22)]
        gls_stage_se = float(gls_cov23[22][22]) ** 0.5
        gls_pub = m2["two_stage_23col_gls_lmm_se_pp"]
        near("mode2 GLS stage_block_B SE", gls_pub.get("stage_block_B_se_pp"), round(gls_stage_se, 2), 0.01)
        near("mode2 GLS min main SE", gls_pub.get("main_effect_half_effect_beta_se_range_pp", [None, None])[0], round(min(gls_main_ses), 2), 0.01)
        near("mode2 GLS max main SE", gls_pub.get("main_effect_half_effect_beta_se_range_pp", [None, None])[1], round(max(gls_main_ses), 2), 0.01)
        near("mode2 GLS min 2FI SE", gls_pub.get("two_factor_interaction_half_effect_beta_se_range_pp", [None, None])[0], round(min(gls_2fi_ses), 2), 0.01)
        near("mode2 GLS max 2FI SE", gls_pub.get("two_factor_interaction_half_effect_beta_se_range_pp", [None, None])[1], round(max(gls_2fi_ses), 2), 0.01)

    va = mde.get("variance_assumptions", {})
    near("variance_assumptions sigma_run", va.get("sigma_run_donor_averaged_batch_cell_pp"), round(float(sigma_run_sq) ** 0.5, 2), 0.01)
    near("variance_assumptions sigma_run 3donors", va.get("sigma_run_if_7pp_is_per_donor_well_residual_divided_by_n_donors_pp"), round((433.0 / 18.0) ** 0.5, 2), 0.01)
    if va.get("clone_level_df_75_clones") != 52 or va.get("clone_level_df_108_clones") != 85:
        fail("variance_assumptions clone_level_df_75_clones / 108_clones mismatch")
    n_distinct_stage = len({(lv, int(r["stage_block"])) for lv, r in zip(levels, runs)})
    if va.get("conservative_run_table_df_23col") != n - 23 or va.get("pure_error_df_within_stage_replicates_23col") != n - n_distinct_stage or va.get("lack_of_fit_df_higher_order_interactions_23col") != n_distinct_stage - 23:
        fail("variance_assumptions 13 df = 9 pure-error + 4 lack-of-fit mismatch")

    sb_qual = mde.get("stage_by_main_interaction_and_batch_covariance_qualification_ultracode_v80_12", {}).get("stage_by_main_29col_expansion", {})
    if sb_qual.get("exact_rational_rank") != d.rank29:
        fail("stage_by_main_29col_expansion.exact_rational_rank != %d" % d.rank29)
    if sb_qual.get("rank_deficiency") != 29 - d.rank29:
        fail("stage_by_main_29col_expansion.rank_deficiency != %d" % (29 - d.rank29))

    # S55R2-D1-03, D1-04, S2-02: exact rational verification of 40-run bridge 23col, 24col Mode-B indicator, and 24col continuous-Z3 matrices
    b40 = proto.get("hip_core_40run_full_rank_2fi_plus_quadratic_q3_bridge_option_astra_r1", {})
    added_37_40 = b40.get("added_stage1_bridge_runs_37_to_40", [])
    if len(added_37_40) == 4:
        runs40 = runs + added_37_40
        d40 = Design(runs40)
        modes40 = b40.get("operating_modes_and_condition_number_reconciliation_argon_r2", {})
        mode_a = modes40.get("mode_a_primary_binary_bridge_Z3_plus1", {})
        mode_b = modes40.get("mode_b_quadratic_sentinel_dose_curvature_Z3_0", {})
        cont_z3 = mode_b.get("continuous_dose_polynomial_Z3_matrix_astra_r2", {})
        same("40-run 23col Var(B)", mode_a.get("var_beta_B_over_sigma2_23col_exact_fraction"), d40.diag23[22])
        same("40-run 23col max leverage", mode_a.get("max_hat_leverage_23col_exact_fraction"), max(d40.hat23))
        if mode_a.get("condition_number_23col_B_pm1_unit_column_norm_exact") != mode_a.get("condition_number_23col_B_pm1_raw"):
            fail("40-run 23col B_pm1 unit-column-norm exact condition number must equal raw condition number 14.7461")
        # 24-col Mode-B indicator matrix (retaining binary X3=+1 on Runs 37-40, adding Q3 = Z3^2 - 9/10)
        q3_col = [Fraction(1, 10) if i < 36 else Fraction(-9, 10) for i in range(40)]
        x24_ind = [list(r) + [q3_col[i]] for i, r in enumerate(d40.x23)]
        inv24_ind = invert(gram(x24_ind))
        hat24_ind = [quad(r, inv24_ind) for r in x24_ind]
        same("40-run 24col Mode-B indicator max leverage", mode_b.get("max_hat_leverage_24col_mode_b_indicator_exact_fraction"), max(hat24_ind))
        same("40-run 24col Mode-B indicator Var(B)", mode_b.get("var_beta_B_over_sigma2_24col_mode_b_indicator_exact_fraction"), inv24_ind[22][22])
        # 24-col continuous-Z3 polynomial matrix (substituting Z3=0 into X3 and all 5 X3*X_k columns on Runs 37-40)
        x24_cont = []
        for i, run in enumerate(runs40):
            lv = [Fraction(int(run["coded_levels"][k])) for k in FACTORS]
            if i >= 36:
                lv[2] = Fraction(0)
            row = [Fraction(1)] + lv + [lv[a] * lv[b] for a, b in PAIRS] + [Fraction(int(run["stage_block"])), q3_col[i]]
            x24_cont.append(row)
        inv24_cont = invert(gram(x24_cont))
        hat24_cont = [quad(r, inv24_cont) for r in x24_cont]
        same("40-run 24col continuous-Z3 max leverage", cont_z3.get("max_hat_leverage_24col_continuous_Z3_exact_fraction"), max(hat24_cont))
        same("40-run 24col continuous-Z3 Var(B)", cont_z3.get("var_beta_B_over_sigma2_24col_continuous_Z3_exact_fraction"), inv24_cont[22][22])
        # AUDIT-REF & ALG-04: verify Runs 37-40 decoded_genotype & locus1_physical_cassette_id match bridged source runs 35, 06, 07, 14
        by_num = {r["run_number"]: r for r in runs}
        for r_add, src_num in zip(added_37_40, (35, 6, 7, 14)):
            src_r = by_num[src_num]
            if r_add.get("locus1_physical_cassette_id") != src_r.get("locus1_physical_cassette_id"):
                fail("Run %d locus1_physical_cassette_id %r != source Run %d %r" % (r_add["run_number"], r_add.get("locus1_physical_cassette_id"), src_num, src_r.get("locus1_physical_cassette_id")))
            if r_add.get("decoded_genotype") != src_r.get("decoded_genotype"):
                fail("Run %d decoded_genotype != source Run %d decoded_genotype" % (r_add["run_number"], src_num))
        # AUDIT-REF: verify 46-run duplicated-bridge 29-column leverage (7/8) and Var(beta_B)/sigma^2 (3/4)
        runs46 = runs40 + [by_num[25], by_num[26]] + list(added_37_40)
        x23_46 = model_rows(runs46, True)
        x29_46 = [r + [r[22] * r[k] for k in range(1, 7)] for r in x23_46]
        inv29_46 = invert(gram(x29_46))
        hat29_46 = [quad(r, inv29_46) for r in x29_46]
        if max(hat29_46) != Fraction(7, 8) or inv29_46[22][22] != Fraction(3, 4):
            fail("46-run 29col max leverage %s != 7/8 or Var(B) %s != 3/4" % (plain(max(hat29_46)), plain(inv29_46[22][22])))
    # AUDIT-REF: verify split-batch S1-A vs S1-B 2-batch residual GLS SEs under Mode 1 (108 clones) and Mode 2 (75 clones)
    b_sweep = mde.get("stage_by_main_interaction_and_batch_covariance_qualification_ultracode_v80_12", {}).get("differentiation_batch_random_intercept_sensitivity_sweep", {}).get("batch_sd_2p5_pp", {})
    m2_2b = b_sweep.get("stage1_split_batch_s1a_vs_s1b_incidence_2batch_residual_gls_se_sol61_r3", {})
    m1_2b = b_sweep.get("stage1_split_batch_s1a_vs_s1b_mode1_independent_run_nested_clones_2batch_residual_gls_se_kimi31_r3", {})
    if m1_2b.get("main_effect_half_effect_beta_gls_se_range_pp") != [1.2359, 1.3086] or m1_2b.get("stage_block_B_gls_se_pp") != 5.154:
        fail("Mode-1 split-batch S1-A/S1-B 2-batch residual GLS SE mismatch")
    if m2_2b.get("main_effect_half_effect_beta_gls_se_range_pp") != [1.2655, 1.3461] or m2_2b.get("stage_block_B_gls_se_pp") != 3.7139:
        fail("Mode-2 split-batch S1-A/S1-B 2-batch residual GLS SE mismatch")
    return problems


# ---------------------------------------------------------------------------------------------------------------------------------
# documentation parity
# ---------------------------------------------------------------------------------------------------------------------------------
def headline_tokens(proto):
    vv = proto["exact_rational_parameter_variance_and_vif_conditioning"]
    up = vv["optional_38run_4bridge_upgrade_specification"]
    base = vv["base_34_run_design_reference"]
    lead = vv["lead_chassis_a_prediction_variance"]
    inflation = {"main-effect block inflation %s%%" % ("%.2f" % m["block_inflation_pct"]): "%.2f%%" % m["block_inflation_pct"] for m in vv["main_effects_var_multiplier_by_factor"].values()}
    return dict(inflation, **{
        "stage-block variance Var(B)": plain(frac_of(vv["stage_block_B"]["var_23col_rational"])),
        "stage-block projection SSR": plain(frac_of(vv["stage_block_B"]["projection_residual_ssr"])),
        "max leverage": proto["exact_rational_leverage_profile_with_stage_block_23cols"]["max_leverage_rational"],
        "max 2FI multiplier": plain(frac_of(vv["two_factor_interaction_var_multiplier_range_23col"]["max"])),
        "min 2FI multiplier": plain(frac_of(vv["two_factor_interaction_var_multiplier_range_23col"]["min"])),
        "lead prediction variance (stage 0)": plain(frac_of(lead["var_multiplier_stage0"])),
        "lead prediction variance (stage 1)": plain(frac_of(lead["var_multiplier_stage1"])),
        "38-run Var(B)": plain(frac_of(up["stage_block_var_multiplier"])),
        "38-run max 2FI multiplier": up["max_2fi_var_multiplier_rational"],
        "rank with stage x main interactions": vv["stage_block_B"]["rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols"].replace(" ", ""),
        "38-run rank with stage x main interactions": up["rank_if_augmented_with_all_6_stage_by_main_effect_interactions_29cols"].replace(" ", ""),
        "base 34-run max 2FI multiplier": plain(frac_of(base["max_2fi_var_multiplier_23col"])),
        "base 34-run min 2FI multiplier": plain(frac_of(base["min_2fi_var_multiplier_23col"])),
    })


def verify_docs(proto):
    problems = []
    n = len(proto["runs"])
    docs = {p.relative_to(ROOT).as_posix(): p.read_text(encoding="utf-8") for p in [ROOT / "README.md", ROOT / "CLAUDE.md"] + sorted((ROOT / "docs").glob("*.md")) if p.exists()}
    for rel in [r for r in ("docs/ISLET-1-DERIVATIONS-AND-AUDIT-TABLES.md", "CLAUDE.md") if r in docs]:
        text = docs[rel].replace(" ", "")
        for label, token in headline_tokens(proto).items():
            if token.replace(" ", "") not in text:
                problems.append("%s does not quote the %s %s from the protocol JSON" % (rel, label, token))
        if ("%d-run" % n) not in docs[rel]:
            problems.append("%s does not describe the protocol as the %d-run design" % (rel, n))
    for rel, text in docs.items():
        for lineno, line in enumerate(text.splitlines(), 1):
            for token in STALE_TOKENS:
                if token in line:
                    problems.append("%s:%d carries stale design claim %r" % (rel, lineno, token))
            for m in re.finditer(r"\b34[- ]run", line, re.I):
                window = line[max(0, m.start() - 120): m.end() + 120]
                if not QUALIFIER_RE.search(window):
                    problems.append("%s:%d mentions a '34-run' design without saying it is the base/legacy design" % (rel, lineno))
    if MANIFEST_PATH.exists():
        gate = json.loads(MANIFEST_PATH.read_text(encoding="utf-8")).get("publication_freeze", {}).get("pilot_experiment_gate", "")
        if ("%d-run" % n) not in gate:
            problems.append("manifest.json publication_freeze.pilot_experiment_gate does not say %d-run: %r" % (n, gate))
    return problems


# ---------------------------------------------------------------------------------------------------------------------------------
# pytest entry points and CLI
# ---------------------------------------------------------------------------------------------------------------------------------
def test_published_design_algebra_is_exact():
    problems = verify_protocol(load_protocol())
    if problems:
        raise AssertionError("\n".join(problems))


def test_docs_quote_the_computed_design_figures():
    problems = verify_docs(load_protocol())
    if problems:
        raise AssertionError("\n".join(problems))


def summarize(proto):
    d = Design(proto["runs"])
    return {
        "n_runs": d.n,
        "rank_X": d.rank22,
        "rank_X_plus_stage_block_B": d.rank23,
        "unit_leverage_rows": sorted({proto["runs"][i]["run_id"] for i, h in enumerate(d.hat23) if h == 1}),
        "max_leverage_23col": plain(max(d.hat23)),
        "var_B_multiplier": plain(d.diag23[22]),
        "max_2fi_multiplier": plain(max(d.diag23[7:22])),
        "min_2fi_multiplier": plain(min(d.diag23[7:22])),
        "main_effect_block_inflation_pct": [pct_inflation(d.diag22[k], d.diag23[k]) for k in range(1, 7)],
        "lead_genotype_prediction_variance": [plain(d.prediction_variance(LEAD_LEVELS, 0)), plain(d.prediction_variance(LEAD_LEVELS, 1))],
    }


def main():
    proto = load_protocol()
    problems = verify_protocol(proto) + verify_docs(proto)
    print("HIP-Core exact-algebra gate: %s" % json.dumps(summarize(proto), sort_keys=True))
    if problems:
        print("FAIL (%d)" % len(problems))
        for p in problems:
            print("  -", p)
        return 1
    print("PASS: every published design figure equals its exact-rational recomputation, and the docs quote the same figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
