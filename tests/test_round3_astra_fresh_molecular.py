"""Fresh molecular source-identity probes and independent numerical oracles.

Run with --basetemp beneath the audit report directory; all mutations are copies.
No construct generation, atlas export, or source-data rewriting is performed.
"""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import run_layer2_perturbation_and_regulatory_plane as l2
from scripts import run_layer3_structural_binder_plane as l3


def require(condition, message="Fresh molecular oracle failed"):
    if not condition:
        raise AssertionError(message)


@pytest.fixture(scope="module")
def source():
    atlas = json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text())
    ref = json.loads((ROOT / "data/raw/layer2_genomic_and_jaspar_reference.json").read_text())
    return atlas["construct_sequence_and_grammar_verification"], ref


def resealed_precursor_mutant(cv, prefix, aa_position, replacement, path):
    """Independently retranslate, rehash and reslice a copied donor after one change."""
    mutant = copy.deepcopy(cv)
    records = l2.parse_construct_fasta_text(l2.FASTA_OUT.read_text())
    coding = {}
    by_id = {entry["construct_id"]: entry for entry in mutant["constructs"]}
    for record, cid in l2.CANONICAL_FASTA_CONSTRUCT_IDS.items():
        con, dna = by_id[cid], records[record]
        coding[cid] = {}
        for part in con["part_ledger"]:
            if part["part_type"] != "CDS":
                continue
            lo, hi = part["start_1based"] - 1, part["end_1based"]
            segment = dna[lo:hi]
            oriented = segment if part["orientation"] == "+" else l2.revcomp(segment)
            if part["part_name"].startswith(prefix):
                offset = 3 * (aa_position - 1)
                require(l2.translate_dna(oriented)[aa_position - 1] != replacement)
                oriented = oriented[:offset] + l2.CODON_TABLE[replacement] + oriented[offset + 3:]
                segment = oriented if part["orientation"] == "+" else l2.revcomp(oriented)
                dna = dna[:lo] + segment + dna[hi:]
                part["sha256"] = l2.sha256_str(segment)
                part["translated_aa_sequence"] = l2.translate_dna(oriented)
            coding[cid][part["part_name"]] = l2.translate_dna(oriented)
        con["donor_sha256"] = l2.sha256_str(dna)
        records[record] = dna
    a, b, c, d = [coding[cid] for cid in l2.CANONICAL_FASTA_CONSTRUCT_IDS.values()]
    def one(parts, head):
        return next(value for name, value in parts.items() if name.startswith(head))
    p2a = one(a, "GSG_P2A_")
    expected = {
        "HLAE_0103_VMAPRTLVL_SCT_Arm1B_mature_surface": one(b, "HLAE_")[20:-1],
        "HLAE_0103_VMAPRTLVL_SCT_Arm1A_mature_with_cytosolic_P2A_stub": one(a, "HLAE_")[20:] + p2a[:-1],
        "CD47_IgV_GPI_mature_cell_surface": one(a, "CD47_")[18:147],
        "SERPINB9_cytosolic_with_21aa_GSG_P2A_stub": one(c, "SERPINB9_") + p2a[:-1],
        "PDL1_mature_cell_surface": one(c, "PDL1_")[18:-1],
        "iCasp9_cytosolic_monomer": one(d, "iCasp9_")[:-1],
    }
    for key, aa in expected.items():
        entry = mutant["processed_polypeptide_ledger"][key]
        entry["sequence"], entry["sha256"] = aa, l2.sha256_str(aa)
        if key.startswith("CD47_"):
            entry["cd47_ectodomain_120aa_sha256"] = l2.sha256_str(aa[:120])
    text = "".join(f">{key} length={len(dna)} sha256={l2.sha256_str(dna)}\n{dna}\n" for key, dna in records.items())
    path.write_text(text)
    mutant["fasta_sha256"] = l2.sha256_str(text)
    return mutant


def resealed_noncoding_mutant(cv, name, path):
    mutant = copy.deepcopy(cv)
    records = l2.parse_construct_fasta_text(l2.FASTA_OUT.read_text())
    by_id = {entry["construct_id"]: entry for entry in mutant["constructs"]}
    for record, cid in l2.CANONICAL_FASTA_CONSTRUCT_IDS.items():
        con, dna = by_id[cid], records[record]
        for part in con["part_ledger"]:
            if part["part_name"] == name:
                lo, hi = part["start_1based"] - 1, part["end_1based"]
                replacement = "A" if dna[lo] != "A" else "T"
                dna = dna[:lo] + replacement + dna[lo + 1:]
                part["sha256"] = l2.sha256_str(dna[lo:hi])
        con["donor_sha256"] = l2.sha256_str(dna)
        records[record] = dna
    text = "".join(f">{key} length={len(dna)} sha256={l2.sha256_str(dna)}\n{dna}\n" for key, dna in records.items())
    path.write_text(text)
    mutant["fasta_sha256"] = l2.sha256_str(text)
    return mutant


@pytest.mark.parametrize("prefix,position,replacement", [
    ("HLAE_", 2, "P"),                 # B2M-derived signal peptide
    ("HLAE_", 270, "R"),               # HLA-E mature Gly107 -> Arg107
    ("CD47_", 2, "P"),                 # CD47 signal peptide
    ("CD47_", 175, "P"),               # discarded GPI propeptide
    ("PDL1_", 2, "P"),                 # PD-L1 signal peptide
    ("SERPINB9_", 340, "D"),            # reactive-center P1 identity
    ("iCasp9_", 37, "F"),              # revert engineered FKBP F36V
])
def test_resealed_donor_cannot_inherit_wrong_reference_precursor_processing(source, tmp_path, prefix, position, replacement):
    path = tmp_path / "resealed.fasta"
    mutant = resealed_precursor_mutant(source[0], prefix, position, replacement, path)
    l3.validate_processed_polypeptide_ledger(mutant["processed_polypeptide_ledger"])
    with pytest.raises(ValueError, match="named reference precursor"):
        l2.validate_construct_fasta_binding(mutant, path)


@pytest.mark.parametrize("accession", sorted(l2.REFERENCE_UNIPROT_SEQUENCE_SHA256))
def test_native_source_identity_cannot_be_changed_under_existing_accession(source, accession):
    ref = copy.deepcopy(source[1])
    entry = ref["uniprot_protein_sequences"][accession]
    entry["sequence"] = "A" + entry["sequence"][1:]
    with pytest.raises(ValueError, match="frozen reference identity"):
        l2.validated_reference_coding_modules(ref)


@pytest.mark.parametrize("mutation", ["sequence", "coordinate"])
def test_homology_registry_is_not_authenticated_by_its_own_resealed_input(source, mutation):
    ref = copy.deepcopy(source[1])
    key = "B2M_HA_L_800bp"
    entry = ref["ensembl_verified_homology_arms_and_elements"][key]
    pin = ref["validated_homology_arm_source_registry"][key]
    if mutation == "sequence":
        entry["sequence"] = ("A" if entry["sequence"][0] != "A" else "T") + entry["sequence"][1:]
        entry["sha256"] = pin["sequence_sha256"] = l2.sha256_str(entry["sequence"])
    else:
        for field in ("start_1based", "end_1based"):
            entry[field] -= 1
            pin[field] -= 1
        entry["region"] = f"GRCh38:{entry['chrom']}:{entry['start_1based']}-{entry['end_1based']}:strand=+1"
    with pytest.raises(ValueError, match="frozen GRCh38 reference"):
        l2.validate_homology_arm_sources(ref)
    cv = copy.deepcopy(source[0])
    cv["validated_homology_arm_source_registry"] = ref["validated_homology_arm_source_registry"]
    require(l2.verify_construct_part_grammar(cv)["all_construct_grammars_verified"] is False)


@pytest.mark.parametrize("name", ["cHS4_core_insulator_250bp", "Chimeric_Intron_133bp", "WPREmut6_589bp"])
@pytest.mark.parametrize("mutation", ["relabel", "flip", "misclassify"])
def test_reference_grammar_requires_all_regulatory_parts(source, name, mutation):
    mutant = copy.deepcopy(source[0])
    for con in mutant["constructs"]:
        for part in con["part_ledger"]:
            if part["part_name"] == name:
                if mutation == "relabel":
                    part["part_name"] = "unrelated_filler"
                elif mutation == "flip":
                    part["orientation"] = "-"
                else:
                    part["part_type"] = "PROMOTER"
    require(l2.verify_construct_part_grammar(mutant)["all_construct_grammars_verified"] is False)


@pytest.mark.parametrize("name", ["cHS4_core_insulator_250bp", "Chimeric_Intron_133bp", "WPREmut6_589bp",
    "EEF1A1_core_promoter_230bp", "Syn_5xISRE_minP_148bp_sense", "bGH_polyA_225bp", "SV40_late_polyA_122bp_antisense"])
def test_resealed_noncoding_bytes_cannot_reuse_reference_part_labels(source, tmp_path, name):
    path = tmp_path / "noncoding.fasta"
    mutant = resealed_noncoding_mutant(source[0], name, path)
    with pytest.raises(ValueError, match="grammar is invalid"):
        l2.validate_construct_fasta_binding(mutant, path)


def test_unchanged_four_donors_bind_and_reference_modules_are_complete(source):
    modules = l2.validated_reference_coding_modules(source[1])
    require({key: len(value) for key, value in modules.items()} == {
        "HLAE_": 500, "CD47_": 175, "SERPINB9_": 376, "PDL1_": 290, "iCasp9_": 396, "GSG_P2A_": 22})
    require({key: l2.sha256_str(value) for key, value in modules.items()} == l2.REFERENCE_CODING_MODULE_SHA256)
    result = l2.validate_construct_fasta_binding(source[0], l2.FASTA_OUT)
    require(result["processing_status"] == "EXPECTED_NOT_EMPIRICALLY_MEASURED")


@pytest.mark.parametrize("name", sorted(l3.STRUCTURAL_REFERENCE_SHA256))
def test_reference_source_bytes_cannot_be_replaced_under_fixed_structure_labels(tmp_path, monkeypatch, name):
    original = (l3.STRUCT_DIR / name).read_bytes()
    (tmp_path / name).write_bytes(original + b"\n")
    monkeypatch.setattr(l3, "STRUCT_DIR", tmp_path)
    with pytest.raises(ValueError, match="source identity mismatch"):
        l3.verify_structural_reference_inputs(name)


def test_all_retained_coordinate_and_pae_inputs_match_independent_release_pins():
    l3.verify_structural_reference_inputs(*l3.STRUCTURAL_REFERENCE_SHA256)


def test_named_complex_rejects_finite_coordinate_mutation_before_geometry(source, tmp_path, monkeypatch):
    text = (l3.STRUCT_DIR / "pdb_00004zqk.cif").read_text()
    # Change exactly one atom coordinate, preserving valid mmCIF syntax.
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("ATOM "):
            fields = line.split()
            fields[10] = str(float(fields[10]) + 0.001)
            lines[i] = " ".join(fields)
            break
    (tmp_path / "pdb_00004zqk.cif").write_text("\n".join(lines) + "\n")
    monkeypatch.setattr(l3, "STRUCT_DIR", tmp_path)
    with pytest.raises(ValueError, match="source identity mismatch"):
        l3.analyze_complex_4_pdl1_pd1(source[0]["processed_polypeptide_ledger"])


@pytest.fixture(scope="module")
def af():
    atoms = l3.parse_mmcif_atoms(l3.STRUCT_DIR / "AF-P50453-F1-model_v6.cif")
    metadata = json.loads((l3.STRUCT_DIR / "AF-P50453-F1-metadata.json").read_text())
    return atoms, metadata


@pytest.mark.parametrize("field,value", [("modelEntityId", "AF-WRONG-F1"), ("entryId", "AF-WRONG-F1"),
    ("latestVersion", 5), ("latestVersion", True), ("chainId", "B"),
    ("cifUrl", "https://alphafold.ebi.ac.uk/files/AF-P50453-F1-model_v5.cif"),
    ("paeDocUrl", "https://alphafold.ebi.ac.uk/files/AF-P50453-F1-predicted_aligned_error_v5.json")])
def test_alphafold_model_version_identity_is_bound(af, field, value):
    metadata = dict(af[1], **{field: value})
    with pytest.raises(ValueError, match="model/version identity mismatch"):
        l3.validate_alphafold_reference(af[0], metadata, "P50453", 376)


def test_alphafold_confidence_is_residue_weighted_with_metadata_difference_preserved(af):
    original = copy.deepcopy(af[1])
    result = l3.validate_alphafold_reference(af[0], af[1], "P50453", 376)
    ca = [atom["b_iso"] for atom in af[0] if atom["atom_id"] == "CA"]
    expected = float(sum(Decimal(str(value)) for value in ca) / Decimal(376))
    require(math.isclose(result["global_plddt_audit"]["coordinate_mean_unrounded"], expected, abs_tol=1e-12))
    require(result["coordinate_global_plddt"] == 91.16)
    require(result["global_plddt_audit"]["metadata_supplied_global_plddt"] == 91.19)
    require(result["global_plddt_audit"]["agreement_at_two_decimal_places"] is False)
    require(af[1] == original, "Validator must not rewrite supplied metadata")


def test_alphafold_atom_confidence_and_identity_cannot_disagree_with_its_residue(af):
    atoms = copy.deepcopy(af[0])
    next(atom for atom in atoms if atom["atom_id"] != "CA")["b_iso"] = 99.99
    with pytest.raises(ValueError, match="disagrees within its declared residue"):
        l3.validate_alphafold_reference(atoms, af[1], "P50453", 376)


def test_coordinate_confidence_reaches_complex_output(source):
    result = l3.analyze_complex_3_serpinb9_gzmb(source[0]["processed_polypeptide_ledger"])
    models = result["structures_analyzed"][:2]
    require([model["global_plddt"] for model in models] == [91.16, 92.41])
    require(all(model["global_plddt_source_audit"]["agreement_at_two_decimal_places"] is False for model in models))


@pytest.fixture(scope="module")
def reference_complex_outputs(source):
    ledger = source[0]["processed_polypeptide_ledger"]
    return {name: getattr(l3, name)(ledger) for name in (
        "analyze_complex_1_hlae_nkg2a", "analyze_complex_2_cd47_sirpa",
        "analyze_complex_3_serpinb9_gzmb", "analyze_complex_4_pdl1_pd1")}


@pytest.mark.parametrize("analysis_name,name,pae_label", [
    ("analyze_complex_1_hlae_nkg2a", "pdb_00003cdg.cif", None),
    ("analyze_complex_2_cd47_sirpa", "pdb_00004cmm.cif", None),
    ("analyze_complex_3_serpinb9_gzmb", "pdb_00001iau.cif", None),
    ("analyze_complex_3_serpinb9_gzmb", "AF-P50453-F1-model_v6.cif", None),
    ("analyze_complex_3_serpinb9_gzmb", "AF-P50453-F1-predicted_aligned_error_v6.json", "SERPINB9"),
    ("analyze_complex_3_serpinb9_gzmb", "AF-P10144-F1-predicted_aligned_error_v6.json", "Granzyme B"),
    ("analyze_complex_4_pdl1_pd1", "pdb_00004zqk.cif", None),
])
@pytest.mark.parametrize("restore_after_consume", [False, True])
def test_named_analysis_consumes_verified_snapshot_during_file_replacement(
        source, reference_complex_outputs, tmp_path, monkeypatch,
        analysis_name, name, pae_label, restore_after_consume):
    """Change real disk bytes after verification, possibly restoring after parse.

    A preflight plus final recheck is insufficient: the restore case must still
    consume the verified original, with identical A/B conformers and provenance.
    """
    original_dir = l3.STRUCT_DIR
    original_bytes = (original_dir / name).read_bytes()
    if pae_label is not None:
        payload = json.loads(original_bytes)
        record = payload[0] if isinstance(payload, list) else payload
        record["predicted_aligned_error"] = [
            [0.0] * len(row) for row in record["predicted_aligned_error"]]
        altered_bytes = json.dumps(payload).encode()
    else:
        lines = original_bytes.decode().splitlines()
        for i, line in enumerate(lines):
            if line.startswith("ATOM "):
                fields = line.split()
                fields[10] = str(float(fields[10]) + 0.125)
                lines[i] = " ".join(fields)
                break
        altered_bytes = ("\n".join(lines) + "\n").encode()
    require(original_bytes != altered_bytes)
    shutil.copytree(original_dir, tmp_path / "structures")
    structures = tmp_path / "structures"
    target = structures / name
    original_verify = l3.verify_structural_reference_inputs
    original_parse = l3.parse_mmcif_atoms
    original_pae = l3.validated_pae_matrix
    observed = []

    def replace_after_verified_read(*names):
        snapshots = original_verify(*names)
        require(name in names)
        target.write_bytes(altered_bytes)
        return snapshots

    def observe_parse(path, *args, **kwargs):
        atoms = original_parse(path, *args, **kwargs)
        if Path(path).name == name:
            expected_kwargs = {key: value for key, value in kwargs.items() if key != "source_bytes"}
            expected = original_parse(original_dir / name, *args, **expected_kwargs)
            observed.append(atoms == expected)
            if restore_after_consume:
                target.write_bytes(original_bytes)
        return atoms

    def observe_pae(payload, expected_length, label):
        matrix = original_pae(payload, expected_length, label)
        if label == pae_label:
            expected = original_pae(json.loads(original_bytes), expected_length, label)
            observed.append(np.array_equal(matrix, expected))
            if restore_after_consume:
                target.write_bytes(original_bytes)
        return matrix

    monkeypatch.setattr(l3, "STRUCT_DIR", structures)
    monkeypatch.setattr(l3, "verify_structural_reference_inputs", replace_after_verified_read)
    monkeypatch.setattr(l3, "parse_mmcif_atoms", observe_parse)
    monkeypatch.setattr(l3, "validated_pae_matrix", observe_pae)
    result = getattr(l3, analysis_name)(source[0]["processed_polypeptide_ledger"])
    require(observed and all(observed), "An unverified replacement reached the coordinate/PAE consumer")
    require(result == reference_complex_outputs[analysis_name], "Geometry or provenance came from a different read")
    require(target.read_bytes() == (original_bytes if restore_after_consume else altered_bytes))
    if pae_label is not None:
        model = result["structures_analyzed"][0 if pae_label == "SERPINB9" else 1]
        require(model["pae_sha256"] == hashlib.sha256(original_bytes).hexdigest())
    else:
        structure_rows = result.get("structures_analyzed", result.get("experimental_structures_analyzed", []))
        require(any(row["sha256"] == hashlib.sha256(original_bytes).hexdigest()
                    for row in structure_rows))


def test_verified_structure_snapshot_is_immutable_and_generic_parser_keeps_file_api(tmp_path, monkeypatch):
    name = "pdb_00004zqk.cif"
    snapshots = l3.verify_structural_reference_inputs(name)
    with pytest.raises(TypeError):
        snapshots[name] = b"replacement"
    with pytest.raises(TypeError):
        snapshots[name][0] = 0
    expected = l3.parse_mmcif_atoms(l3.STRUCT_DIR / name)
    missing = tmp_path / name
    require(l3.parse_mmcif_atoms(missing, source_bytes=snapshots[name]) == expected)
    with pytest.raises(FileNotFoundError):
        l3.parse_mmcif_atoms(missing)


def test_fresh_pwm_dp_matches_full_probability_enumeration_on_each_strand():
    seq = "AACGT"
    weights = {"A": [0.1, 0.7, -0.1], "C": [0.2, -0.2, 0.3],
               "G": [-0.3, 0.1, 0.4], "T": [0.4, -0.1, 0.5]}
    entry = {"matrix_id": "FRESH_ENUM", "tf_name": "SYNTHETIC", "length": 3, "log2_odds_pwm": weights}
    result = l2.scan_jaspar_pwm(seq, entry, top_k=6, n_matrices_in_scan=1)
    for hit in result["top_hits"]:
        background = seq if hit["strand"] == "+" else l2.revcomp(seq)
        tail = Fraction(0)
        threshold = round(hit["log2_odds_score"] * 10000)
        for word in itertools.product("ACGT", repeat=3):
            if sum(round(weights[base][j] * 10000) for j, base in enumerate(word)) >= threshold:
                prob = Fraction(1)
                for base in word:
                    prob *= Fraction(background.count(base), len(background))
                tail += prob
        require(math.isclose(hit["exact_markov0_p_value_unadjusted"], float(tail), abs_tol=1e-12))


def test_fresh_kabsch_reflection_retains_closed_form_residual():
    tetrahedron = np.array([[1., 1., 1.], [-1., -1., 1.], [-1., 1., -1.], [1., -1., -1.]])
    reflected = tetrahedron * [1., 1., -1.]
    rmsd, aligned, target = l3.kabsch_rmsd(tetrahedron, reflected)
    require(math.isclose(rmsd, 2., abs_tol=1e-14))
    require(math.isclose(float(np.sqrt(np.mean(np.sum((aligned - target) ** 2, axis=1)))), 2., abs_tol=1e-14))


def test_fresh_wlc_force_and_second_moment_match_independent_decimal_equations():
    with localcontext() as context:
        context.prec = 70
        contour, span, persistence = Decimal("24.5"), Decimal("12"), Decimal("6")
        x = span / contour
        force = Decimal("42.8005") / persistence * (1 / (4 * (1 - x) ** 2) - Decimal("0.25") + x)
        rms = (2 * persistence * contour - 2 * persistence ** 2 * (1 - (-contour / persistence).exp())).sqrt()
    result = l3.worm_like_chain_metrics(7, 12., lp=6.)
    require(result["marko_siggia_entropic_tension_pN_at_selected_lp"] == round(float(force), 3))
    require(result["wlc_rms_end_to_end_angstrom"] == round(float(rms), 3))
