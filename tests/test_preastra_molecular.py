"""Independent source/negative controls for PREASTRA-L2-A01..09 and L3-A01..06.

All mutations and structural replays are temporary; no canonical artifact is written.
"""
from __future__ import annotations
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import run_layer2_perturbation_and_regulatory_plane as l2
from scripts import run_layer3_structural_binder_plane as l3


def require(value, message="independent molecular assertion failed"):
    if not value:
        raise AssertionError(message)


def load(relative):
    return json.loads((ROOT / relative).read_text())


@pytest.fixture
def ref():
    return load("data/raw/layer2_genomic_and_jaspar_reference.json")


@pytest.fixture
def atlas():
    return load("data/extracted/layer2_regulatory_and_perturbation_atlas.json")


@pytest.fixture
def cv(atlas):
    return copy.deepcopy(atlas["construct_sequence_and_grammar_verification"])


def test_source_registry_baseline(ref):
    require(len(l2.validate_homology_arm_sources(ref)) == 8)


@pytest.mark.parametrize("mutation", ["overlap", "reversed", "wrong_region", "wrong_strand", "wrong_chrom", "wrong_sequence", "wrong_length", "wrong_assembly", "missing_registry"])
def test_source_registry_rejects_raw_mutants(ref, mutation):
    arms = ref["ensembl_verified_homology_arms_and_elements"]
    entry = arms["B2M_HA_R_600bp"]
    if mutation == "overlap":
        entry.update(start_1based=44711009, end_1based=44711608, region="GRCh38:chr15:44711009-44711608:strand=+1")
    elif mutation == "reversed":
        entry["start_1based"], entry["end_1based"] = entry["end_1based"], entry["start_1based"]
    elif mutation == "wrong_region": entry["region"] = "GRCh38:chr15:1-600:strand=+1"
    elif mutation == "wrong_strand": entry["strand"] = -1
    elif mutation == "wrong_chrom": entry["chrom"] = "chr16"
    elif mutation == "wrong_sequence": entry["sequence"] = "A" * 600
    elif mutation == "wrong_length": entry["length_bp"] = 599
    elif mutation == "wrong_assembly": entry["assembly"] = "GRCh37"
    else: del ref["validated_homology_arm_source_registry"]
    with pytest.raises(ValueError): l2.validate_homology_arm_sources(ref)


@pytest.mark.parametrize("mutation", ["missing_locus4", "duplicate_locus4", "missing_locus4_audit", "missing_icasp9", "wrong_promoter", "destroyed_p2a", "reversed_ha_interval"])
def test_incomplete_or_corrupt_grammar_rejects(cv, ref, atlas, mutation):
    if mutation == "missing_locus4": cv["constructs"] = cv["constructs"][:-1]
    elif mutation == "duplicate_locus4": cv["constructs"].append(copy.deepcopy(cv["constructs"][-1]))
    elif mutation == "missing_locus4_audit": del cv["regulatory_element_audits"]["CBX3_A2UCOE_core_455bp_R17"]
    else:
        target = cv["constructs"][-1] if mutation in ("missing_icasp9", "wrong_promoter") else cv["constructs"][0]
        for part in target["part_ledger"]:
            if mutation == "missing_icasp9" and "iCasp9" in part["part_name"]: part["part_name"] = "unidentified_CDS"
            elif mutation == "wrong_promoter" and part["part_type"] == "PROMOTER": part["orientation"] = "-"
            elif mutation == "destroyed_p2a" and "P2A" in part["part_name"]: part["translated_aa_sequence"] = "A" * 22
            elif mutation == "reversed_ha_interval" and "_HA_R_" in part["part_name"]: part["source_accession_or_coordinate"] = "GRCh38:chr15:44712408-44711609:+1"
    try:
        result = l2.verify_construct_part_grammar(cv)
    except ValueError:
        result = {"all_construct_grammars_verified": False}
    require(result["all_construct_grammars_verified"] is False, mutation)
    with pytest.raises((ValueError, AssertionError)):
        l2.simulate_stage0_to_stage7_perturbation_plane(ref, cv, atlas["jaspar_pwm_regulatory_motif_scans"])


def test_candidate_identity_and_native_chs4_are_computed(ref, cv):
    candidate = cv["regulatory_element_audits"]["CBX3_A2UCOE_core_455bp_R17"]
    require(candidate["canonical_display_name"] == "HNRNPA2B1-proximal 455-bp CpG-core candidate (GRCh38 chr7:26200580-26201034)")
    require(candidate["anti_silencing_functional_status"] == "PREREGISTERED_PENDING")
    elements = ref["ensembl_verified_homology_arms_and_elements"]
    insert = elements["cHS4_core_insulator_250bp"]["sequence"]
    result = l2.compute_chs4_alignment(insert, ref)
    require([result[k] for k in ["longest_exact_contiguous_match_bp", "insert_start_1based", "reference_start_1based", "left_flank_outside_match_bp", "right_flank_outside_match_bp"]] == [235, 3, 10, 2, 13])
    require(result == cv["regulatory_element_audits"]["cHS4_core_insulator_250bp_ULTRACODE_V80_13"]["computed_native_accession_alignment"])
    mutated = insert[:40] + ("A" if insert[40] != "A" else "C") + insert[41:]
    require(l2.compute_chs4_alignment(mutated, ref) != result)
    ref["primary_accession_references"]["U78775.2"]["sequence"] += "A"
    with pytest.raises(ValueError): l2.compute_chs4_alignment(insert, ref)


def test_avi_quarantine_is_not_a_model_and_public_mode_needs_no_archive(atlas, monkeypatch, tmp_path):
    monkeypatch.setattr(l2, "AVI_ARCHIVE_PATH", tmp_path / "missing.json")
    monkeypatch.setattr(l2, "AVI_SEED_PATH", tmp_path / "missing-mirror.json")
    summary = l2.recompute_avi_summaries_from_full_archive()
    require(summary == atlas["alphagenome_avi_saturation_mutagenesis_and_variants"])
    require(summary["source_origin"] == "UNRESOLVED_ARCHIVAL_TABLE")
    require(summary["origin_and_redistribution_rights"] == "UNVERIFIED")
    require(summary["total_snvs_scored_and_archived"] is None)
    require(all(row["max_avi_phred"] is None for row in summary["window_saturation_mutagenesis_summaries"].values()))


def test_avi_private_mirror_checks_integrity_not_generation(monkeypatch, tmp_path):
    archive_path = ROOT / "data/raw/alphagenome_avi_3118_snvs_full_archive.json"
    mirror_path = ROOT / "data/raw/alphagenome_avi_26_seed_anchors.json"
    if not archive_path.exists() or not mirror_path.exists():
        pytest.skip("Private unresolved-origin archive intentionally excluded from this distribution")
    archive, mirror = json.loads(archive_path.read_text()), json.loads(mirror_path.read_text())
    require(l2.verify_avi_integrity_mirror(archive)["supplied_rows_count"] == 3118)
    window = next(iter(archive["windows"]))
    # Non-highlighted identical edits in both supplied tables legitimately pass mirror integrity.
    archive["windows"][window]["variants"][-1]["avi_phred"] += 0.001
    ledger = mirror["deterministic_window_reconstruction_ledger"][window]
    rows = archive["windows"][window]["variants"]
    ledger["compact_rows"] = rows
    ledger["row_schema"] = "full_variant_dict"
    digest = hashlib.sha256(json.dumps(rows, indent=2).encode()).hexdigest()
    ledger["serialized_variants_array_sha256"] = digest
    archive["windows"][window]["serialized_variants_array_sha256"] = digest
    path = tmp_path / "mirror.json"; path.write_text(json.dumps(mirror)); monkeypatch.setattr(l2, "AVI_SEED_PATH", path)
    result = l2.verify_avi_integrity_mirror(archive)
    require(result["model_predictions_provided"] is False)
    require(result["integrity_status"] == "SUPPLIED_COPIES_MATCH_NOT_INDEPENDENT_GENERATION")
    archive["windows"][window]["variants"][-1]["avi_phred"] += 1
    with pytest.raises(ValueError): l2.verify_avi_integrity_mirror(archive)


def test_profiles_have_typed_products_prior_bands_and_derived_notes(atlas):
    plane = atlas["stage0_to_stage7_perturbation_plane"]
    profiles = {p["genotype_id"]: p for p in plane["perturbation_profiles"]}
    require(plane["profile_count"] == len(profiles) == 11)
    a, b, wt = [profiles[x] for x in ["SCENARIO_A_ARM_1A_BICISTRONIC_P2A_COMPARATOR", "SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS", "WT_UNEDITED_CONTROL"]]
    require(abs(a["hour0_basal_cd47_molecules_per_cell"] / b["hour0_basal_cd47_molecules_per_cell"] - .88) < 2 / b["hour0_basal_cd47_molecules_per_cell"])
    note = plane["hour0_hlae_inducible_vs_constitutive_kinetics_note"]
    require(abs(note["arm1b_cd47_to_wt_fold"] - b["hour0_basal_cd47_molecules_per_cell"] / wt["hour0_basal_cd47_molecules_per_cell"]) < 1e-12)
    for p in profiles.values():
        require(p["hour0_basal_hlae_molecules_per_cell"] == p["hour0_resting_unstimulated_hlae_molecules_per_cell"])
        require(p["product_definition_id"])
        if p["pdl1_downstream_2a_skip_factor_prior"] is not None:
            require(p["pdl1_downstream_2a_skip_factor_prior"] == .88)
            require(p["ifng_24h_pdl1_molecules_per_cell"] == round(p["ifng_24h_pdl1_before_2a_skip_prior_molecules_per_cell"] * .88))
        require("stage7_sc_beta_yield_prior_scenario_interval_p10_p90" in p)
        require("no posterior fitting" in p["uncertainty_semantics"])
        for key, meta in p["deprecated_uncertainty_aliases"].items(): require(p[key] == p[meta["deprecated_alias_of"]])
        if p["genotype_id"].startswith("SCENARIO_B"):
            require(p["exact_product_defined"] is False and p["source_construct"] is None)
            require(p["no_scenario_b_product_benefit_claim"] is True and p["product_definition"]["execution_supported"] is False)
    require(a["product_definition_id"] == "ISLET1_A4_ARM1A" and b["product_definition_id"] == "ISLET1_A4_ARM1B")
    ko = profiles["HIP_CORE_X4_TYK2_KO_PERTURBATION"]
    changes = note["tyk2_relative_yield_comparisons"]
    require(abs(changes["stage6_relative_change_pct"] - 100 * (ko["stage6_neurog3_ep_fraction"] / wt["stage6_neurog3_ep_fraction"] - 1)) < 1e-12)
    require(abs(changes["stage7_relative_change_pct"] - 100 * (ko["stage7_sc_beta_yield_fraction"] / wt["stage7_sc_beta_yield_fraction"] - 1)) < 1e-12)


def test_exact_encoded_signal_and_single_p2a(cv, ref):
    require(ref["ensembl_gencode_release_pairing"].startswith("Ensembl release 112 / GENCODE v46"))
    for c in cv["constructs"]:
        for part in c["part_ledger"]:
            if "HLAE" in part["part_name"] and part["part_type"] == "CDS": require(part["translated_aa_sequence"][:20] == "MSRSVALAVLALLSLSGLEA")
            if "P2A" in part["part_name"]: require(part["length_bp"] == 66 and part["translated_aa_sequence"] == "GSGATNFSLLKQAGDVEENPGP")


def make_structural_fixture(tmp_path):
    for rel in ["scripts/run_layer3_structural_binder_plane.py", "scripts/run_layer2_perturbation_and_regulatory_plane.py", "scripts/product_contract.py", "tests/test_layer2_and_layer3_planes.py", "data/extracted/layer2_regulatory_and_perturbation_atlas.json", "data/extracted/layer3_structural_verification_atlas.json", "data/extracted/ring1_atomic_effects_v0.5.jsonl", "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta", "data/raw/layer2_genomic_and_jaspar_reference.json", "data/protocols/islet1_product_definitions_v0.8.0.json"]:
        p = tmp_path / rel; p.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(ROOT / rel, p)
    shutil.copytree(ROOT / "data/structures", tmp_path / "data/structures")
    return tmp_path


def test_layer3_raw_mutation_rejects_stale_pass_normal_and_optimized(tmp_path):
    fixture = make_structural_fixture(tmp_path)
    spec = importlib.util.spec_from_file_location("molecular_suite", ROOT / "tests/test_layer2_and_layer3_planes.py")
    suite = importlib.util.module_from_spec(spec); spec.loader.exec_module(suite)
    committed = load("data/extracted/layer3_structural_verification_atlas.json")
    fresh = suite.verify_fresh_layer3_parity(committed)
    require(fresh["summary"]["all_reference_coordinate_geometric_screens_passed"] is True)
    af = fixture / "data/structures/AF-P10144-F1-model_v6.cif"
    cols, out = [], []
    for line in af.read_text().splitlines():
        if line.startswith("_atom_site."): cols.append(line.split(".", 1)[1].strip())
        if line.startswith(("ATOM ", "HETATM ")):
            fields = line.split(); x, residue = cols.index("Cartn_x"), cols.index("label_seq_id")
            fields[x] = f"{float(fields[x]) + 15 * math.sin(1.7 * int(fields[residue])):.3f}"; line = " ".join(fields)
        out.append(line)
    af.write_text("\n".join(out) + "\n")
    # Fixed-source analyses now reject substituted coordinates before creating
    # any differently labelled atlas. Both producer and gate must fail closed.
    atlas_path = fixture / "data/extracted/layer3_structural_verification_atlas.json"
    original_atlas_bytes = atlas_path.read_bytes()
    for optimized in [False, True]:
        for relative in ["scripts/run_layer3_structural_binder_plane.py", "tests/test_layer2_and_layer3_planes.py"]:
            command = [sys.executable, "-B"] + (["-O"] if optimized else []) + [str(fixture / relative)]
            result = subprocess.run(command, cwd=fixture, text=True, capture_output=True)
            require(result.returncode != 0 and "Structural reference source identity mismatch" in result.stderr,
                    f"Fixed-source analysis accepted raw mutant under optimized={optimized}: {result.stderr[-1200:]}")
            require(atlas_path.read_bytes() == original_atlas_bytes, "Rejected source must preserve the previous atlas")


def test_layer3_reference_qualification_hash_alias_and_foldseek_provenance():
    atlas = load("data/extracted/layer3_structural_verification_atlas.json")
    c1, c2, c3, _ = atlas["complexes"]
    contacts = c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis"]["p8_contacts_by_receptor_subunit_in_3CDG"]
    require([contacts[x] for x in ["copy1_P_to_CD94_chainJ_heavy_contacts", "copy1_P_to_NKG2A_chainK_heavy_contacts", "copy2_Q_to_CD94_chainE_heavy_contacts", "copy2_Q_to_NKG2A_chainF_heavy_contacts"]] == [24, 0, 21, 0])
    require(c1["vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis"]["p8_coordinates_discriminate_nkg2a_over_nkg2c"] is False)
    hashes = c2["layer2_processed_polypeptide_hashes"]
    require(hashes["CD47_crystal_2JJS_chainC_resolved_116aa_C15G_sha256"] == hashes["CD47_crystal_2JJS_chainC_resolved_116aa_with_PCA1_C15G_sha256"])
    require(hashes["CD47_crystal_2JJS_chainC_resolved_115aa_ATOM_CA_C15G_sha256"] != hashes["CD47_crystal_2JJS_chainC_resolved_116aa_C15G_sha256"])
    provenance = c3["foldseek_structural_homolog_verification"]["search_provenance"]
    require(provenance["qualification_status"] == "HISTORICAL_CACHED_HITS_DATABASE_SNAPSHOT_NOT_REPRODUCIBLE")
    require(provenance["tool_or_api_version"] is None and provenance["query_date"] is None)
    require(all(v is None for v in provenance["database_snapshot_versions"].values()))
    require(provenance["source_artifact_sha256"] == hashlib.sha256((ROOT / "data/structures/serpinb9_foldseek_results.md").read_bytes()).hexdigest())
    require(atlas["summary"]["all_exact_product_functional_constraints_passed"] is False)


def test_sirpa_insertion_and_wlc_identity_are_reference_geometry():
    atlas = load("data/extracted/layer3_structural_verification_atlas.json")
    c1, c2 = atlas["complexes"][:2]
    audit = c2["gpi_anchor_and_sirpa_polymorphism_verification"]
    require(audit["sirpa_v1_v2_polymorphic_positions_tested"] == [14, 20, 22, 24, 27, 45, 65, 66, 70, 77, 79, 101])
    insertions = [x for x in audit["sirpa_v1_v2_insertion_aware_alignment_differences"] if x["difference_type"] == "V1_single_residue_insertion"]
    require(len(insertions) == 1 and insertions[0]["sirpa_v1_4CMM_seq_id"] == 101 and insertions[0]["sirpa_v1_4CMM_aa"] == "D")
    wlc = c1["single_chain_trimer_linker_physics"]["linker_2_B2M_Met99_to_HLAE_Ser2"]
    require(wlc["unresolved_n_terminal_hlae_residue_identity"] == "Gly1")
    require(wlc["linker_length_aa"] == 21 and wlc["design_linker_length_aa"] == 20)
    require(abs(wlc["contour_length_angstrom"] - 73.5) < 1e-10)
    require(abs(wlc["marko_siggia_entropic_tension_pN_at_310K_lp4A"] - 9.901) < .001)
    metrics = c1["hlae_0103_gly107_stability_verification"]
    require(metrics["crystal_residue_107_identity_in_3CDG"] == "ARG")
    require(metrics["crystal_arg107_ca_distance_to_cd94_nkg2a_interface_angstrom"] == 17.02)
    require(metrics["crystal_arg107_min_heavy_atom_distance_to_cd94_nkg2a_interface_angstrom"] == 14.82)


def test_cd47_hashes_reconstruct_from_raw_cif_without_atlas_aliases():
    import shlex
    three = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL PCA".split(), "ARNDCQEGHILKMFPSTWYVQ"))
    cols, residues, pca = [], {}, False
    for line in (ROOT / "data/structures/pdb_00002jjs.cif").read_text().splitlines():
        if line.startswith("_atom_site."): cols.append(line.split(".", 1)[1].strip())
        elif line.startswith(("ATOM ", "HETATM ")):
            row = dict(zip(cols, shlex.split(line)))
            if row["auth_asym_id"] != "C" or row["label_atom_id"] != "CA" or row["label_alt_id"] not in (".", "?", "A"): continue
            if row["label_comp_id"] == "PCA":
                require(int(row["auth_seq_id"]) == 1); pca = True
            elif row["group_PDB"] == "ATOM": residues[int(row["auth_seq_id"])] = three[row["label_comp_id"]]
    seq115 = "".join(residues[k] for k in sorted(residues))
    require(pca and len(seq115) == 115)
    seq116 = "Q" + seq115
    hashes = load("data/extracted/layer3_structural_verification_atlas.json")["complexes"][1]["layer2_processed_polypeptide_hashes"]
    for key in hashes:
        if "116aa" in key: require(hashes[key] == hashlib.sha256(seq116.encode()).hexdigest())
        if "115aa" in key: require(hashes[key] == hashlib.sha256(seq115.encode()).hexdigest())


def test_foldseek_parser_live_and_missing_hit_rejects(monkeypatch, tmp_path, cv):
    structures = tmp_path / "structures"
    shutil.copytree(ROOT / "data/structures", structures)
    monkeypatch.setattr(l3, "STRUCT_DIR", structures)
    md = structures / "serpinb9_foldseek_results.md"
    text = md.read_text()
    require("Yan, Zhou et al. 2024 Acta Crystallogr Sect F Struct Biol Commun 80:286-293" in text)
    require("Li et al. 2024 J. Biol. Chem." not in text)
    md.write_text(text.replace("| 8h4x-assembly1.cif.gz_A Structure of cleaved Serpin B9 | 85.9%", "| 8h4x-assembly1.cif.gz_A Structure of cleaved Serpin B9 | 84.9%"))
    result = l3.analyze_complex_3_serpinb9_gzmb(cv["processed_polypeptide_ledger"])
    require(result["foldseek_structural_homolog_verification"]["same_protein_human_serpinb9_pdb_hit_and_supersession"]["q_cov_pct"] == 84.9)
    md.write_text("\n".join(line for line in text.splitlines() if not line.startswith("| 8h4x-assembly1.cif.gz_A")))
    with pytest.raises(ValueError, match="Required Foldseek hit"):
        l3.analyze_complex_3_serpinb9_gzmb(cv["processed_polypeptide_ledger"])


def raw_cif_atom_rows(path):
    """Independent reader for the geometric test oracles."""
    import shlex
    columns, rows = [], []
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."): columns.append(line.split(".", 1)[1].strip())
        elif line.startswith(("ATOM ", "HETATM ")):
            r = dict(zip(columns, shlex.split(line)))
            if r.get("pdbx_PDB_model_num", "1") != "1" or r["label_alt_id"] not in (".", "?", "A") or r["type_symbol"].upper() == "H": continue
            rows.append({"chain": r["auth_asym_id"], "residue": int(r["auth_seq_id"]), "atom": r["label_atom_id"], "xyz": tuple(float(r[k]) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))})
    return rows


def test_wlc_force_recomputed_from_raw_endpoints_and_boltzmann_constant():
    atoms = raw_cif_atom_rows(ROOT / "data/structures/pdb_00003cdg.cif")
    def point(chain, residue, atom):
        hits = [r["xyz"] for r in atoms if (r["chain"],r["residue"],r["atom"]) == (chain,residue,atom)]
        require(len(hits) == 1)
        return hits[0]
    distance = math.sqrt(sum((a-b)**2 for a,b in zip(point("B",99,"C"),point("A",2,"N"))))
    wlc = load("data/extracted/layer3_structural_verification_atlas.json")["complexes"][0]["single_chain_trimer_linker_physics"]["linker_2_B2M_Met99_to_HLAE_Ser2"]
    require(abs(wlc["spanned_euclidean_distance_angstrom"]-distance) < .001)
    extension = distance / (21 * 3.5)
    force = (1.380649e-23 * 310 * 1e22 / 4) * (1/(4*(1-extension)**2)-.25+extension)
    require(abs(wlc["marko_siggia_entropic_tension_pN_at_310K_lp4A"]-force) < .001)


def test_shared_cd94_p8_contacts_recomputed_independently():
    atoms = raw_cif_atom_rows(ROOT / "data/structures/pdb_00003cdg.cif")
    def count(p,r):
        return sum(sum((x-y)**2 for x,y in zip(a["xyz"],b["xyz"])) <= 4.5**2 for a in atoms if a["chain"] == p and a["residue"] == 8 for b in atoms if b["chain"] == r)
    require([count("P","J"),count("P","K"),count("Q","E"),count("Q","F")] == [24,0,21,0])


def test_kinetics_note_tracks_changed_expression_inputs(ref,cv,atlas,tmp_path):
    changed = copy.deepcopy(ref)
    changed["gtex_v8_v10_human_tissue_expression"]["CD47"]["pancreas_median_tpm"] *= 1.5
    # Freeze the shared corpus once: isolate expression changes from concurrent corpus probes.
    effects = tmp_path / "effects.jsonl"
    effects.write_bytes(l2.RING1_EFFECTS_PATH.read_bytes())
    base = l2.simulate_stage0_to_stage7_perturbation_plane(ref,cv,atlas["jaspar_pwm_regulatory_motif_scans"],effects_path=effects)
    result = l2.simulate_stage0_to_stage7_perturbation_plane(changed,cv,atlas["jaspar_pwm_regulatory_motif_scans"],effects_path=effects)
    profiles = {p["genotype_id"]:p for p in result["perturbation_profiles"]}
    original = {p["genotype_id"]:p for p in base["perturbation_profiles"]}
    note = result["hour0_hlae_inducible_vs_constitutive_kinetics_note"]
    for entry in note["profiles"]:
        p = profiles[entry["genotype_id"]]
        for label,key in [("resting_hlae_molecules_per_cell","hour0_resting_unstimulated_hlae_molecules_per_cell"),("ifng24h_hlae_molecules_per_cell","ifng_24h_induced_hlae_molecules_per_cell"),("hour0_cd47_molecules_per_cell","hour0_basal_cd47_molecules_per_cell")]: require(entry[label] == p[key])
        require(abs(p["hour0_basal_cd47_molecules_per_cell"]-1.5*original[p["genotype_id"]]["hour0_basal_cd47_molecules_per_cell"]) <= 1)
    require(note["arm1b_cd47_to_wt_fold"] == profiles["SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS"]["hour0_basal_cd47_molecules_per_cell"]/profiles["WT_UNEDITED_CONTROL"]["hour0_basal_cd47_molecules_per_cell"])


def test_native_chs4_adaptation_ledger_identifies_real_bases(ref,cv):
    ledger = cv["regulatory_element_audits"]["cHS4_core_insulator_250bp_ULTRACODE_V80_13"]["computed_native_accession_alignment"]["anchored_native_adjacent_adaptation_ledger"]
    native = ref["primary_accession_references"]["U78775.2"]["sequence"]
    insert = ref["ensembl_verified_homology_arms_and_elements"]["cHS4_core_insulator_250bp"]["sequence"]
    require(ledger["native_bases_omitted_before_following_exact_match"] == native[244:245] == "T")
    require(ledger["native_omission_interval_1based_closed"] == [245,245])
    require(ledger["following_exact_native_sequence"] == native[245:251] == insert[237:243] == "AGGCTG")
    require(ledger["remaining_engineered_tail_sequence"] == insert[243:] == "CCGCGGG")
    require(ledger["remaining_tail_contains_sacii_recognition_motif"] is True)


def test_layer2_prior_plane_replays_current_source_inputs(ref,cv,atlas):
    fresh = l2.simulate_stage0_to_stage7_perturbation_plane(ref,cv,atlas["jaspar_pwm_regulatory_motif_scans"])
    require(fresh == atlas["stage0_to_stage7_perturbation_plane"],
            "Layer2 prior plane is stale: current raw expression/effects/grammar replay differs from committed atlas")
    provenance = fresh["input_provenance"]["ring1_atomic_effects"]
    require(provenance["path"] == "data/extracted/ring1_atomic_effects_v0.5.jsonl")
    require(provenance["sha256"] == hashlib.sha256(l2.RING1_EFFECTS_PATH.read_bytes()).hexdigest())
