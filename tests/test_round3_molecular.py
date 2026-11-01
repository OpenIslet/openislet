"""Round-3 molecular numerical oracles and temporary malformed-input probes."""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import run_layer2_perturbation_and_regulatory_plane as l2
from scripts import run_layer3_structural_binder_plane as l3


def require(condition, message="Round-3 independent molecular oracle failed"):
    if not condition:
        raise AssertionError(message)


@pytest.fixture(scope="module")
def source():
    atlas = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text())
    reference = json.loads((ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text())
    return atlas["construct_sequence_and_grammar_verification"], reference


def retile(parts):
    cursor = 1
    for part in parts:
        part["start_1based"] = cursor
        part["end_1based"] = cursor + part["length_bp"] - 1
        cursor += part["length_bp"]


@pytest.mark.parametrize("index", range(4))
@pytest.mark.parametrize("mutation", ["left_flip", "right_flip", "internal", "different_locus", "coordinate_suffix"])
def test_homology_contract_requires_terminal_correct_locus_and_donor_orientation(source, index, mutation):
    construct = copy.deepcopy(source[0])
    selected = construct["constructs"][index]
    parts = selected["part_ledger"]
    if mutation == "left_flip":
        parts[0]["orientation"] = "-"
    elif mutation == "right_flip":
        parts[-1]["orientation"] = "-"
    elif mutation == "internal":
        parts[-1], parts[-2] = parts[-2], parts[-1]
        retile(parts)
    elif mutation == "coordinate_suffix":
        parts[0]["source_accession_or_coordinate"] += "wrong"
    else:
        other = construct["constructs"][(index + 2) % 4]["part_ledger"]
        for destination, replacement in [(parts[0], other[0]), (parts[-1], other[-1])]:
            for key in ["part_name", "raw_reference_key", "source_accession_or_coordinate", "sha256"]:
                destination[key] = replacement[key]
    result = l2.verify_construct_part_grammar(construct)
    require(not result["all_construct_grammars_verified"], "malformed homology contract was certified")


@pytest.mark.parametrize("token,expected", [("+", "+1"), ("+1", "+1"), ("1", "+1"), ("-", "-1"), ("-1", "-1")])
def test_grch38_strand_tokens_normalize_without_partial_matching(token, expected):
    require(l2._parse_grch38_coord(f"GRCh38:chr15:1-2:{token}") == ("chr15", 1, 2, expected))
    require(l2._parse_grch38_coord(f"GRCh38:chr15:1-2:{token}extra") is None)


def test_unchanged_reference_constructs_bind_to_current_dna_and_consistent_alias_lengths(source):
    construct = copy.deepcopy(source[0])
    for entry in construct["processed_polypeptide_ledger"].values():
        entry["mature_length_aa"] = entry["processed_length_aa"] = len(entry["sequence"])
    l3.validate_processed_polypeptide_ledger(construct["processed_polypeptide_ledger"])
    result = l2.validate_construct_fasta_binding(construct, ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta")
    require(result["processed_product_count"] == 6)


@pytest.mark.parametrize("product_index", range(6))
def test_all_present_processed_length_labels_must_agree_with_dna(source, product_index):
    construct = copy.deepcopy(source[0])
    entry = list(construct["processed_polypeptide_ledger"].values())[product_index]
    missing_key = "processed_length_aa" if "mature_length_aa" in entry else "mature_length_aa"
    entry[missing_key] = 1
    with pytest.raises(ValueError, match="length"):
        l3.validate_processed_polypeptide_ledger(construct["processed_polypeptide_ledger"])
    with pytest.raises(ValueError, match="expected product"):
        l2.validate_construct_fasta_binding(construct, ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta")


def small_pwm():
    return {"matrix_id": "ROUND3_ENUMERATED", "tf_name": "synthetic_numerical_fixture", "length": 3,
            "log2_odds_pwm": {"A": [1., -.5, .25], "C": [-1., 1., -.25],
                              "G": [0., .5, 1.25], "T": [.5, -1., -.5]}}


def enumerated_tail(sequence, pwm, word, strand):
    counts = {base: Fraction(sequence.count(base), len(sequence)) for base in "ACGT"}
    if strand == "-":
        counts = {base: counts[l2.revcomp(base)] for base in "ACGT"}
    target = sum(Fraction(str(pwm[base][j])) for j, base in enumerate(word))
    total = Fraction(0)
    for candidate in itertools.product("ACGT", repeat=len(word)):
        score = sum(Fraction(str(pwm[base][j])) for j, base in enumerate(candidate))
        if score >= target:
            probability = Fraction(1)
            for base in candidate:
                probability *= counts[base]
            total += probability
    return float(total)


def test_pwm_exact_tails_and_multiplicity_match_complete_word_enumeration():
    sequence, entry = "AAAACGTT", small_pwm()
    result = l2.scan_jaspar_pwm(sequence, entry, top_k=100, n_matrices_in_scan=7)
    require(len(result["top_hits"]) == 12)
    for hit in result["top_hits"]:
        expected = enumerated_tail(sequence, entry["log2_odds_pwm"], hit["matched_kmer"], hit["strand"])
        require(math.isclose(hit["exact_markov0_p_value_unadjusted"], expected, rel_tol=2e-11))
        require(math.isclose(hit["bonferroni_region_across_all_matrices_p_value"], min(1., expected * 12 * 7), rel_tol=2e-10))
        require(0. <= hit["relative_score_0_to_1"] <= 1.)
        start, end = hit["start_1based_in_window"], hit["end_1based_in_window"]
        expected_word = sequence[start - 1:end]
        if hit["strand"] == "-":
            expected_word = l2.revcomp(expected_word)
        require(hit["matched_kmer"] == expected_word)


def test_pwm_exact_reverse_complement_scan_is_coordinate_and_tail_equivariant():
    sequence, entry = "AAAACGTT", small_pwm()
    forward = l2.scan_jaspar_pwm(sequence, entry, top_k=100)
    reverse = l2.scan_jaspar_pwm(l2.revcomp(sequence), entry, top_k=100)
    lookup = {(hit["strand"], hit["start_1based_in_window"]): hit for hit in reverse["top_hits"]}
    for hit in forward["top_hits"]:
        opposite = lookup[("-" if hit["strand"] == "+" else "+", len(sequence) - hit["end_1based_in_window"] + 1)]
        for key in ["matched_kmer", "log2_odds_score", "exact_markov0_p_value_unadjusted",
                    "bonferroni_region_across_all_matrices_p_value"]:
            require(hit[key] == opposite[key], key)


@pytest.mark.parametrize("sequence,test_count,skipped", [("AC", 0, 0), ("NNN", 2, 2), ("ACGTN", 6, 2)])
def test_pwm_window_boundaries_and_ambiguous_windows_are_explicit(sequence, test_count, skipped):
    result = l2.scan_jaspar_pwm(sequence, small_pwm(), top_k=100)
    require(result["total_positions_tested_both_strands"] == test_count)
    require(result["skipped_ambiguous_windows_count"] == skipped)
    require(len(result["top_hits"]) == test_count - skipped)


def test_off_lattice_pwm_is_rejected_before_inconsistent_probability_calibration():
    entry = {"matrix_id": "ROUND3_OFF_GRID", "tf_name": "synthetic", "length": 1,
             "log2_odds_pwm": {"A": [.00006], "C": [-.00004], "G": [-.00004], "T": [-.00004]}}
    with pytest.raises(ValueError, match="0.0001-bit lattice"):
        l2.scan_jaspar_pwm("AACGT", entry)


def test_pfm_rounding_cannot_certify_a_nine_unit_pwm_mismatch(source):
    entry = copy.deepcopy(source[1]["jaspar_pwms"]["MA0079.5"])
    entry["log2_odds_pwm"]["A"][0] = round(entry["log2_odds_pwm"]["A"][0] + .0009, 4)
    with pytest.raises(ValueError, match="PFM-to-PWM mismatch"):
        l2.scan_jaspar_pwm("ACGTACGTACGT", entry)


@pytest.mark.parametrize("weight", [300., 1e308])
def test_pwm_dense_lattice_limit_rejects_before_input_controlled_allocation(weight):
    entry = {"matrix_id": "ROUND3_OVERSIZE", "tf_name": "synthetic", "length": 1,
             "log2_odds_pwm": {"A": [weight], "C": [0.], "G": [0.], "T": [0.]}}
    with pytest.raises(ValueError, match="lattice"):
        l2.scan_jaspar_pwm("ACGT", entry)


def cif_rows(rows, columns=None, wrap=False):
    default = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "auth_comp_id",
               "auth_asym_id", "auth_seq_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
               "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_PDB_model_num"]
    columns = columns or default
    tokens = [[str(row.get(column, {"group_PDB": "ATOM", "id": i + 1, "type_symbol": "C",
               "label_atom_id": "CA", "label_alt_id": ".", "auth_comp_id": "ALA", "auth_asym_id": "A",
               "auth_seq_id": i + 1, "label_seq_id": i + 1, "pdbx_PDB_ins_code": "?", "Cartn_x": i,
               "Cartn_y": 0, "Cartn_z": 0, "occupancy": 1, "B_iso_or_equiv": 90,
               "pdbx_PDB_model_num": 1}[column])) for column in columns] for i, row in enumerate(rows)]
    if wrap:
        values = list(itertools.chain.from_iterable(tokens))
        lines = [" ".join(values[i:i + 7]) for i in range(0, len(values), 7)]
    else:
        lines = [" ".join(row) for row in tokens]
    return "data_round3\nloop_\n" + "\n".join("_atom_site." + column for column in columns) + "\n" + "\n# internal comment\n".join(lines) + "\n#\n"


@pytest.mark.parametrize("wrap", [False, True])
def test_mmcif_column_order_row_wrapping_and_comments_do_not_change_atom_identity(tmp_path, wrap):
    standard = cif_rows([{}, {}])
    columns = ["id", "group_PDB", "Cartn_z", "Cartn_x", "label_alt_id", "Cartn_y", "type_symbol",
               "auth_comp_id", "label_atom_id", "auth_asym_id", "auth_seq_id", "label_seq_id",
               "pdbx_PDB_ins_code", "B_iso_or_equiv", "occupancy", "pdbx_PDB_model_num"]
    p = tmp_path / "coordinates.cif"
    p.write_text(standard)
    expected = l3.parse_mmcif_atoms(p)
    p.write_text(cif_rows([{}, {}], columns, wrap))
    require(l3.parse_mmcif_atoms(p) == expected)


def test_mmcif_preserves_unquoted_prime_atom_name_and_model_conformer_occupancy_filters(tmp_path):
    rows = [{"label_atom_id": "C1'", "auth_seq_id": 1},
            {"auth_seq_id": 2, "pdbx_PDB_model_num": 2},
            {"auth_seq_id": 3, "occupancy": 0},
            {"auth_seq_id": 4, "label_alt_id": "B"},
            {"auth_seq_id": 5, "label_alt_id": "A"}]
    p = tmp_path / "selection.cif"
    p.write_text(cif_rows(rows))
    atoms = l3.parse_mmcif_atoms(p)
    require([(atom["seq_id"], atom["atom_id"]) for atom in atoms] == [(1, "C1'"), (5, "CA")])
    with pytest.raises(ValueError, match="one named atom conformer"):
        l3.parse_mmcif_atoms(p, allowed_alt_ids=(".", "A", "B"))


def test_mmcif_incomplete_and_duplicate_rows_fail_closed(tmp_path):
    p = tmp_path / "invalid.cif"
    p.write_text(cif_rows([{}]).replace(" 90 1\n", " 90\n"))
    with pytest.raises(ValueError, match="Incomplete atom row"):
        l3.parse_mmcif_atoms(p)
    p.write_text(cif_rows([{"auth_seq_id": 1}, {"auth_seq_id": 1}]))
    with pytest.raises(ValueError, match="Duplicate selected atom identity"):
        l3.parse_mmcif_atoms(p)


def test_ca_polymer_order_uses_label_sequence_with_author_insertion_codes(tmp_path):
    rows = [{"auth_seq_id": 185, "label_seq_id": 3},
            {"auth_seq_id": 184, "pdbx_PDB_ins_code": "A", "label_seq_id": 2},
            {"auth_seq_id": 184, "label_seq_id": 1}]
    p = tmp_path / "insertion.cif"
    p.write_text(cif_rows(rows))
    atoms = l3.parse_mmcif_atoms(p, keep_insertion_codes=True)
    ordered = l3.get_ca_ordered_list_with_ins_codes(atoms, "A")
    require([(a["seq_id"], a["ins_code"]) for a in ordered] == [(184, ""), (184, "A"), (185, "")])
    require(l3.get_ca_ordered_list_with_ins_codes(list(reversed(atoms)), "A") == ordered)
    for atom in atoms:
        atom.pop("label_seq_id")
    require([(a["seq_id"], a["ins_code"]) for a in l3.get_ca_ordered_list_with_ins_codes(atoms, "A")] == [(184, ""), (184, "A"), (185, "")])


def test_real_granzyme_ca_sequence_is_independent_of_coordinate_row_order():
    atoms = l3.parse_mmcif_atoms(ROOT / "data/structures/pdb_00001iau.cif", keep_insertion_codes=True)
    expected = l3.get_ca_ordered_list_with_ins_codes(atoms, "A")
    require(len(expected) == 226)
    require(l3.get_ca_ordered_list_with_ins_codes(list(reversed(atoms)), "A") == expected)


@pytest.mark.parametrize("scale", [1e-300, 1., 1e300])
def test_kabsch_reflection_has_closed_form_nonzero_proper_rotation_residual(scale):
    tetrahedron = np.array([[1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]]) * scale
    reflection = tetrahedron * [-1., 1., 1.]
    rmsd, _, _ = l3.kabsch_rmsd(tetrahedron, reflection)
    require(math.isclose(rmsd / scale, 2., rel_tol=2e-14), "reflection was fitted as an improper rotation")
    rotation = tetrahedron[:, [1, 2, 0]]
    require(l3.kabsch_rmsd(tetrahedron, rotation)[0] / scale < 2e-14)


@pytest.mark.parametrize("scale", [1e-300, 1., 1e300])
def test_kabsch_degenerate_collapsed_cloud_has_analytic_target_variance(scale):
    p = np.zeros((3, 3))
    q = np.array([[-1., 0., 0.], [0., 0., 0.], [1., 0., 0.]]) * scale
    require(math.isclose(l3.kabsch_rmsd(p, q)[0] / scale, math.sqrt(2. / 3.), rel_tol=2e-14))


def test_kabsch_identical_cloud_at_maximum_finite_endpoint_is_exact():
    p = np.array([[1.7976931348623157e308, 7.04933019228406e307, -1.2272271573865657e307],
                  [-9.583840270745158e307, 7.295485509215793e307, -4.778109374709363e307],
                  [-9.821353204128491e306, 8.398382431041934e307, -8.230348949619342e307],
                  [7.365816733740578e307, -6.348038704285288e307, 7.634748063713624e307],
                  [7.858149653431424e307, 4.852899740514347e307, -2.6111649057717124e306]])
    rmsd, aligned, target = l3.kabsch_rmsd(p, p)
    require(rmsd == 0. and np.array_equal(aligned, p) and np.array_equal(target, p))


@pytest.mark.parametrize("span,contour,persistence", [(1e-16, 1e308, math.ulp(0.)), (1e-15, 1e308, math.ulp(0.)),
                                                       (30., 73.5, 4.), (1e-8, 73.5, 8.)])
def test_wlc_force_matches_600_digit_original_equation_through_subnormal_ratios(span, contour, persistence):
    with localcontext() as context:
        context.prec = 600
        fraction = Decimal.from_float(span) / Decimal.from_float(contour)
        force = (Decimal(1) / (4 * (1 - fraction) ** 2) - Decimal(".25") + fraction)
        expected = force * Decimal.from_float(l3.KBT_310K_PN_ANGSTROM) / Decimal.from_float(persistence)
    result = l3.worm_like_chain_metrics(1, span, l0=contour, lp=persistence)
    require(result["marko_siggia_entropic_tension_pN_at_selected_lp"] == round(float(expected), 3))
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("maximum", [10 ** 400, True, float("nan"), -.1])
def test_pae_maximum_domain_rejects_overflow_and_invalid_values(maximum):
    with pytest.raises(ValueError, match="declared maximum"):
        l3.validated_pae_matrix({"predicted_aligned_error": [[0., 1.], [2., 0.]], "max_predicted_aligned_error": maximum}, 2, "round3")


def test_pae_asymmetric_matrix_remains_unmodified_under_real_maximum_bound():
    matrix = [[0., .25], [12., 0.]]
    result = l3.validated_pae_matrix({"predicted_aligned_error": matrix, "max_predicted_aligned_error": 12.}, 2, "round3")
    require(np.array_equal(result, matrix))


@pytest.mark.parametrize("matrix", [np.full((3, 3), sys.float_info.max),
                                  np.array([[sys.float_info.max, 0.], [0., sys.float_info.max]])])
def test_pae_mean_matches_decimal_convex_mean_at_finite_endpoint(matrix):
    with localcontext() as context:
        context.prec = 400
        expected = sum(Decimal.from_float(float(value)) for value in matrix.flat) / matrix.size
    require(l3._pae_mean(matrix) == float(expected))


def test_pae_mean_preserves_exact_zero_and_rejects_empty_input():
    require(l3._pae_mean(np.zeros((3, 3))) == 0.)
    with pytest.raises(ValueError, match="nonempty"):
        l3._pae_mean(np.empty((0, 0)))


def test_prior_quantiles_have_finite_convex_extreme_value_oracle():
    maximum = sys.float_info.max
    samples = [-maximum, -maximum / 2, maximum / 2, maximum]
    result = l2.summarize_prior_draws(samples)
    with localcontext() as context:
        context.prec = 100
        for name, probability in [("p10", Decimal(".1")), ("median", Decimal(".5")), ("p90", Decimal(".9"))]:
            index = probability * (len(samples) - 1)
            lo, fraction = int(index), index - int(index)
            expected = (1 - fraction) * Decimal.from_float(samples[lo]) + fraction * Decimal.from_float(samples[min(lo + 1, len(samples) - 1)])
            require(math.isfinite(result[name]))
            require(math.isclose(result[name], float(expected), rel_tol=3e-16, abs_tol=.001))


@pytest.mark.parametrize("num_a", [7, 8, 9])
def test_pwm_relative_threshold_uses_exact_four_fifths_integer_boundary(num_a):
    entry = {"matrix_id": "ROUND3_THRESHOLD", "tf_name": "synthetic", "length": 10,
             "log2_odds_pwm": {base: [.0001 if base == "A" else 0.] * 10 for base in "ACGT"}}
    result = l2.scan_jaspar_pwm("A" * num_a + "C" * (10 - num_a) + "T" * 90, entry, top_k=200)
    hit = next(row for row in result["top_hits"] if row["strand"] == "+" and row["start_1based_in_window"] == 1)
    require(hit["relative_score_0_to_1"] == num_a / 10)
    for flag, probability in [("passes_rel_score_0p80_and_p_le_0p05", "bonferroni_window_p_value"),
                              ("passes_region_wide_bonferroni_0p05", "bonferroni_region_across_all_matrices_p_value")]:
        require(hit[flag] == (num_a >= 8 and hit[probability] <= .05), "decision disagrees with exact integer-lattice boundary")


@pytest.mark.parametrize("value", [40., 1e308])
def test_interface_b_factor_constant_moments_are_finite_and_serialize(value):
    def atom(chain, seq, x):
        return dict(group="ATOM", atom_id="CA", element="C", comp_id="ALA", chain=chain,
                    seq_id=seq, ins_code="", res_key=str(seq), x=x, y=0., z=0.,
                    b_iso=value, alt_id=".", occupancy=1.)
    result = l3.compute_interface_metrics([atom("A", 1, 0.), atom("A", 2, 1.)],
                                          [atom("B", 1, 1.5), atom("B", 2, 2.5)])
    moments = result["interface_crystallographic_b_factor_profile_A2"]
    for partner in [1, 2]:
        require(moments[f"partner_{partner}_contact_mean_b_iso_A2"] == value)
        require(moments[f"partner_{partner}_contact_std_b_iso_A2"] == 0.)
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("values", [[-sys.float_info.max, sys.float_info.max],
                                   [1e308, math.nextafter(1e308, 0.)], [0., 1., 2.]])
def test_population_moments_match_decimal_for_extreme_spread_and_local_offsets(values):
    with localcontext() as context:
        context.prec = 900
        numbers = [Decimal.from_float(value) for value in values]
        mean = sum(numbers) / len(numbers)
        sd = (sum((value - mean) ** 2 for value in numbers) / len(numbers)).sqrt()
    observed_mean, observed_sd = l3._finite_population_moments(values)
    require(math.isclose(observed_mean, float(mean), rel_tol=3e-15, abs_tol=1e-300))
    require(math.isclose(observed_sd, float(sd), rel_tol=3e-15, abs_tol=1e-300))


@pytest.mark.parametrize("num_residues,l0", [(20, 3.5), (2, .5), (1, 1e308)])
def test_wlc_near_contour_gap_matches_original_900_digit_equation(num_residues, l0):
    contour = num_residues * l0
    span = math.nextafter(contour, 0.)
    with localcontext() as context:
        context.prec = 900
        extension = Decimal.from_float(span) / Decimal.from_float(contour)
        expected = (1 / (4 * (1 - extension) ** 2) - Decimal(".25") + extension)
        expected *= Decimal.from_float(l3.KBT_310K_PN_ANGSTROM) / Decimal(4)
    observed = l3.worm_like_chain_metrics(num_residues, span, l0=l0)["marko_siggia_entropic_tension_pN_at_selected_lp"]
    require(math.isclose(observed, float(expected), rel_tol=3e-15), "near-contour subtraction lost force accuracy")
