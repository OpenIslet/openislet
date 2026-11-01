"""Independent three-pass molecular regression and metamorphic oracles.

Mutant inputs are temporary; no assay observations or sequence designs are invented.
"""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import random
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import run_layer2_perturbation_and_regulatory_plane as l2
from scripts import run_layer3_structural_binder_plane as l3


def require(condition, message="independent molecular oracle failed"):
    if not condition:
        raise AssertionError(message)


@pytest.fixture(scope="module")
def source():
    atlas = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text())
    fasta = ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta"
    return atlas, fasta


@pytest.mark.parametrize("mutation", ["duplicate", "preamble", "unknown", "ambiguity", "lowercase", "length", "hash", "missing"])
def test_reference_fasta_rejects_silent_parser_loss(source, mutation):
    _, fasta = source
    text = fasta.read_text()
    first_header = text.splitlines()[0]
    if mutation == "duplicate":
        text += first_header + "\n" + "A" * 4880 + "\n"
    elif mutation == "preamble":
        text = "ACGT\n" + text
    elif mutation == "unknown":
        text = text.replace("ISLET1_LOCUS1_ARM1A", "ISLET1_LOCUS1_ARM1X", 1)
    elif mutation in {"ambiguity", "lowercase"}:
        lines = text.splitlines(); lines[1] = ("N" if mutation == "ambiguity" else "a") + lines[1][1:]
        text = "\n".join(lines)
    elif mutation == "length":
        text = text.replace("length=4880", "length=1", 1)
    elif mutation == "hash":
        text = text.replace(first_header, first_header.split("sha256=")[0] + "sha256=" + "0" * 64, 1)
    else:
        text = text[:text.rfind("\n>")]
    with pytest.raises(ValueError):
        l2.parse_construct_fasta_text(text)


def test_fasta_record_order_is_not_a_hidden_identity_key(source):
    _, fasta = source
    text = fasta.read_text()
    shuffled = "\n".join(">" + row.strip() for row in reversed(text.split(">")[1:])) + "\n"
    require(l2.parse_construct_fasta_text(shuffled) == l2.parse_construct_fasta_text(text))


@pytest.mark.parametrize("product", [
    "HLAE_0103_VMAPRTLVL_SCT_Arm1B_mature_surface",
    "HLAE_0103_VMAPRTLVL_SCT_Arm1A_mature_with_cytosolic_P2A_stub",
    "CD47_IgV_GPI_mature_cell_surface", "SERPINB9_cytosolic_with_21aa_GSG_P2A_stub",
    "PDL1_mature_cell_surface", "iCasp9_cytosolic_monomer",
])
def test_self_resealed_polypeptide_must_still_match_dna(source, product):
    atlas, fasta = source
    cv = copy.deepcopy(atlas["construct_sequence_and_grammar_verification"])
    entry = cv["processed_polypeptide_ledger"][product]
    entry["sequence"] = ("A" if entry["sequence"][0] != "A" else "G") + entry["sequence"][1:]
    entry["sha256"] = hashlib.sha256(entry["sequence"].encode()).hexdigest()
    if product == "CD47_IgV_GPI_mature_cell_surface":
        entry["cd47_ectodomain_120aa_sha256"] = hashlib.sha256(entry["sequence"][:120].encode()).hexdigest()
    # The old self-hash gate accepts this; the DNA-derived binding must reject it.
    l3.validate_processed_polypeptide_ledger(cv["processed_polypeptide_ledger"])
    with pytest.raises(ValueError, match="not the expected product"):
        l2.validate_construct_fasta_binding(cv, fasta)


def test_production_layer3_rejects_self_resealed_wrong_product_before_output(source, tmp_path, monkeypatch):
    atlas, _ = source
    mutant = copy.deepcopy(atlas)
    entry = mutant["construct_sequence_and_grammar_verification"]["processed_polypeptide_ledger"]["PDL1_mature_cell_surface"]
    entry["sequence"] = "A" * len(entry["sequence"])
    entry["sha256"] = hashlib.sha256(entry["sequence"].encode()).hexdigest()
    src, out = tmp_path / "atlas.json", tmp_path / "output.json"
    src.write_text(json.dumps(mutant)); out.write_text("previous output\n")
    monkeypatch.setattr(l3, "LAYER2_ATLAS_PATH", src)
    monkeypatch.setattr(l3, "OUT_PATH", out)
    with pytest.raises(ValueError, match="not the expected product"):
        l3.main()
    require(out.read_text() == "previous output\n", "invalid product overwrote a prior output")


@pytest.mark.parametrize("index", range(4))
@pytest.mark.parametrize("mutation", ["kozak_hash", "kozak_name", "total_length", "zero_part", "duplicate_part"])
def test_grammar_binds_initiation_modules_and_complete_positive_tiling(source, index, mutation):
    atlas, _ = source
    cv = copy.deepcopy(atlas["construct_sequence_and_grammar_verification"])
    construct = cv["constructs"][index]
    if mutation == "total_length":
        construct["total_length_bp"] = 1
    elif mutation == "zero_part":
        construct["part_ledger"][1]["length_bp"] = 0
    elif mutation == "duplicate_part":
        construct["part_ledger"][1]["part_name"] = construct["part_ledger"][0]["part_name"]
    else:
        kozak = next(p for p in construct["part_ledger"] if p["part_name"].startswith("Kozak_"))
        if mutation == "kozak_hash": kozak["sha256"] = "0" * 64
        else: kozak["part_name"] = "unknown_initiation_site"
    with pytest.raises(ValueError):
        l2.verify_construct_part_grammar(cv)


@pytest.mark.parametrize("length", [1, 2, 3, 4, 999, 1000])
def test_prior_quantiles_match_independent_numpy_convention(length):
    rng = random.Random(8472)
    draws = [rng.uniform(-25.0, 31.0) for _ in range(length)]
    before = list(draws)
    result = l2.summarize_prior_draws(draws)
    expected = np.quantile(draws, [.5, .1, .9], method="linear")
    require(list(result.values()) == [round(float(x), 3) for x in expected])
    require(draws == before, "summary mutated source draws")


def test_even_sample_median_is_midpoint_not_upper_order_statistic():
    require(l2.summarize_prior_draws([0.0, 10.0]) == {"median": 5.0, "p10": 1.0, "p90": 9.0})


@pytest.mark.parametrize("value", [True, "1", complex(1, 0), float("nan"), float("inf"), 10**400])
def test_draw_summaries_reject_coercions(value):
    with pytest.raises(ValueError): l2.summarize_prior_draws([0, value])


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_scaled_kabsch_keeps_reflection_error_and_rigid_rotation_invariance(scale):
    points = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])
    mirror = points * np.array([-1., 1., 1.])
    require(math.isclose(l3.kabsch_rmsd(points*scale, mirror*scale)[0]/scale, .5, abs_tol=1e-14),
            "reflected tetrahedron error must remain 0.5 times scale")
    rotation = np.array([[0.,1.,0.], [-1.,0.,0.], [0.,0.,1.]])
    target = (points @ rotation + np.array([2.,3.,4.])) * scale
    result, aligned, original = l3.kabsch_rmsd(points*scale, target)
    require(result/scale < 1e-13)
    require(np.max(np.abs(aligned/scale-target/scale)) < 1e-13)
    require(np.array_equal(original, target))


@pytest.mark.parametrize("value", [True, "1", complex(1, 0), complex(1, 1), float("nan"), float("inf")])
def test_geometry_does_not_silently_coerce_coordinates(value):
    p = [[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]]; bad = copy.deepcopy(p); bad[1][0] = value
    with pytest.raises(ValueError): l3.kabsch_rmsd(bad, p)
    with pytest.raises(ValueError): l3.euclidean_distance([value,0,0], [0,0,0])


@pytest.mark.parametrize("distance", [1e-200, 1e200])
def test_euclidean_norm_preserves_finite_nonzero_scale(distance):
    require(l3.euclidean_distance([distance, 0, 0], [0,0,0]) == distance)


def test_unrepresentable_distance_rejects_infinity():
    with pytest.raises(ValueError): l3.euclidean_distance([1e308,0,0], [-1e308,0,0])


def test_kabsch_large_translation_preserves_small_representable_cloud():
    points = 1e16 + np.eye(3) * 4
    target = points - points[0]
    rmsd, aligned, _ = l3.kabsch_rmsd(points, target)
    require(rmsd < 1e-13)
    require(np.max(np.abs(aligned-target)) < 1e-13)


def test_kabsch_collapsed_huge_cloud_does_not_erase_tiny_target_spread():
    points = np.full((3,3), 1e200)
    target = np.eye(3) * 1e-200
    rmsd, aligned, _ = l3.kabsch_rmsd(points, target)
    require(math.isclose(rmsd / 1e-200, math.sqrt(2/3), rel_tol=1e-14))
    require(np.max(np.abs(aligned / 1e-200 - np.full((3,3), 1/3))) < 1e-14)


def test_wlc_small_extension_matches_80_digit_limit():
    with localcontext() as ctx:
        ctx.prec = 80
        x, lp = Decimal("1e-17"), Decimal("1e-17")
        expected = (1 / (4 * (1-x)**2) - Decimal(".25") + x) * Decimal("42.8005") / lp
    result = l3.worm_like_chain_metrics(1, 1e-17, l0=1., lp=1e-17)
    require(result["marko_siggia_entropic_tension_pN_at_selected_lp"] == round(float(expected), 3))


@pytest.mark.parametrize("maximum", [True, -1., float("inf"), "31.75", 1.])
def test_pae_declared_maximum_is_a_real_bound(maximum):
    payload = {"predicted_aligned_error": [[0.,2.], [1.,0.]], "max_predicted_aligned_error": maximum}
    with pytest.raises(ValueError): l3.validated_pae_matrix(payload, 2, "synthetic")


def test_pae_is_allowed_to_be_directionally_asymmetric():
    matrix = [[0., 2.], [1., 0.]]
    require(np.array_equal(l3.validated_pae_matrix({"predicted_aligned_error":matrix,
        "max_predicted_aligned_error":31.75}, 2, "asymmetric"), matrix))


def test_duplicate_selected_atom_rejected_and_zero_occupancy_not_resolved(tmp_path):
    path = ROOT / "data/structures/pdb_00001mhe.cif"
    atoms = l3.parse_mmcif_atoms(path)
    require(all(a["occupancy"] > 0 for a in atoms))
    require(not any(a["atom_id"] == "CA" and a["chain"] == "A" and a["seq_id"] == 223 for a in atoms), "zero-occupancy GLY223 became resolved")
    lines = path.read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("ATOM "))
    lines.insert(index + 1, lines[index])
    mutant = tmp_path / "duplicate.cif"; mutant.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="Duplicate selected atom"):
        l3.parse_mmcif_atoms(mutant)


@pytest.mark.parametrize("mutation", ["missing_residue", "identity", "plddt", "metadata_missing", "metadata_bool"])
def test_alphafold_complete_reference_and_confidence_are_checked(mutation):
    atoms = l3.parse_mmcif_atoms(ROOT / "data/structures/AF-P50453-F1-model_v6.cif")
    metadata = json.loads((ROOT / "data/structures/AF-P50453-F1-metadata.json").read_text())
    record = metadata[0] if isinstance(metadata, list) else metadata
    if mutation == "missing_residue": atoms = [a for a in atoms if a["seq_id"] != 1]
    elif mutation == "identity":
        for atom in atoms:
            if atom["seq_id"] == 1: atom["comp_id"] = "ALA"
    elif mutation == "plddt": atoms[0]["b_iso"] = 101.
    elif mutation == "metadata_missing": del record["globalMetricValue"]
    else: record["globalMetricValue"] = True
    with pytest.raises(ValueError): l3.validate_alphafold_reference(atoms, metadata, "P50453", 376)
