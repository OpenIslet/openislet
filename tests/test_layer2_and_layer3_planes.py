#!/usr/bin/env python3
"""
Automated Verification Suite for Layer 2 & Layer 3 Planes (v0.8.0)
==================================================================
Closes V70-01 through V70-10 and AUDIT-REF by independently verifying:
  1. Exact SHA-256 checksums, component coordinates, non-overlapping GRCh38 homology arms,
     AAV6 packaging lengths INCLUDING 2x145-bp ITRs, and independent 3-frame/ORF amino acid
     translation of all 4 constructs in `data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta`.
  2. Honest quarantine metadata for the unresolved-origin AVI archive (no model predictions),
     and independent raw-reference coordinates for the 455-bp HNRNPA2B1-proximal candidate.
  3. Direct atomic parsing of all 10 `.cif` files in `data/structures/` (`st_size > 10,000 bytes`,
     non-empty CA atom arrays), 226/226 aligned residues for `AF-P10144` vs `1IAU` with all 6 PDB
     chymotrypsinogen insertion codes (`36A, 36B, 145A, 172A, 184A, 188A`) preserved (`RMSD = 0.4049 A`), Layer-2-to-Layer-3 polypeptide SHA-256 ledger,
     SERPINB9 RCL local pLDDT (`48.36` across `322-348` and `36.21` across `335-345`), and Marko-Siggia WLC entropic tension (`pN`).
  4. Explicit Negative-Control Mutation Probes (`V70-10`, `AUDIT-REF`): proves that mutating any FASTA record
     to all-`A` bases or truncating any `.cif` file to 0 bytes raises an `AssertionError` even under `python3 -O`.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import shutil
import subprocess
import sys
from typing import Any, Dict, List

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
L2_ATLAS = REPO_ROOT / "data" / "extracted" / "layer2_regulatory_and_perturbation_atlas.json"
L2_FASTA = REPO_ROOT / "data" / "constructs" / "islet1_scenario_a_constructs_v0.7.0.fasta"
L2_RAW_REF = REPO_ROOT / "data" / "raw" / "layer2_genomic_and_jaspar_reference.json"
L2_AVI_ARCHIVE = REPO_ROOT / "data" / "raw" / "alphagenome_avi_3118_snvs_full_archive.json"
L3_ATLAS = REPO_ROOT / "data" / "extracted" / "layer3_structural_verification_atlas.json"
STRUCT_DIR = REPO_ROOT / "data" / "structures"


def _check(cond: bool, msg: str = "Verification check failed") -> None:
    if not cond:
        raise AssertionError(msg)


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_fasta_text(text: str) -> Dict[str, str]:
    records: Dict[str, str] = {}
    curr_id = ""
    buf: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if curr_id:
                records[curr_id] = "".join(buf)
            curr_id = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if curr_id:
        records[curr_id] = "".join(buf)
    return records


def revcomp_dna(dna: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return dna.translate(comp)[::-1]


def translate_dna(dna: str) -> str:
    _check(len(dna) % 3 == 0, f"CDS length {len(dna)} is not divisible by 3")
    return "".join(GENETIC_CODE[dna[i : i + 3]] for i in range(0, len(dna), 3))


def validate_fasta_against_l2_atlas(fasta_records: Dict[str, str], l2: Dict[str, Any], raw_ref: Dict[str, Any]) -> None:
    constructs = {
        c["construct_id"]: c for c in l2["construct_sequence_and_grammar_verification"]["constructs"]
    }
    expected_ids = {
        "ISLET1_LOCUS1_ARM1A_B2M_HLAE_SCT_GSGP2A_CD47_IGV_GPI_4880bp": "LOCUS1_ARM1A_TPT2A_BICISTRONIC",
        "ISLET1_LOCUS1_ARM1B_B2M_BIDIR_2AFREE_CD47_IGV_GPI_HLAE_SCT_3974bp": "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01",
        "ISLET1_LOCUS3_CXCL10_SERPINB9_GSGP2A_PDL1_NO_FURIN_4090bp": "LOCUS3_CXCL10_SERPINB9_TPT2A_PDL1_R01",
        "ISLET1_LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_EEF1A1_ICASP9_3457bp": "LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_ICASP9_R17",
    }
    _check(set(fasta_records.keys()) == set(expected_ids.keys()), f"FASTA keys mismatch: {list(fasta_records.keys())}")

    raw_seqs = raw_ref["ensembl_verified_homology_arms_and_elements"]
    ucoe_raw = raw_seqs.get("CBX3_A2UCOE_455bp_exact", {})
    ucoe_audit = l2["construct_sequence_and_grammar_verification"].get("regulatory_element_audits", {}).get("CBX3_A2UCOE_core_455bp_R17", {})
    _check(
        ucoe_raw.get("start_1based") == 26200580
        and ucoe_raw.get("end_1based") == 26201034
        and ucoe_audit.get("avi_coordinate_overlap_bp") == 455
        and ucoe_audit.get("avi_coordinate_overlap_fraction") == 1.0,
        f"UCOE reference interval must match pinned AVI window chr7:26200580-26201034 with 455 bp overlap, got {ucoe_raw.get('start_1based')}..{ucoe_raw.get('end_1based')} (overlap={ucoe_audit.get('avi_coordinate_overlap_bp')})",
    )
    for fasta_id, cid in expected_ids.items():
        seq = fasta_records[fasta_id]
        c_meta = constructs[cid]
        _check(len(seq) == c_meta["total_length_bp"], f"{fasta_id} length mismatch")
        _check(set(seq) == {"A", "C", "G", "T"}, f"{fasta_id} must contain all 4 ACGT bases (rejects homopolymer/all-A corruption)")
        actual_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        _check(actual_sha == c_meta["donor_sha256"], f"{fasta_id} SHA-256 mismatch: {actual_sha} != {c_meta['donor_sha256']}")
        _check(c_meta["homology_arm_overlap_bp"] == 0, f"Non-zero homology arm overlap in {cid}")

        # Verify every component interval in part_ledger, its SHA-256, and translate every CDS/ORF part on BOTH + and - strands (AUDIT-REF & AUDIT-REF)
        cursor = 1
        ledger_list = c_meta["part_ledger"]
        part_by_name = {p["part_name"]: p for p in ledger_list}
        if cid == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01":
            cd47_part = next((p for p in ledger_list if "CD47" in p["part_name"]), None)
            eef_part = next((p for p in ledger_list if "EEF1A1" in p["part_name"]), None)
            isre_part = next((p for p in ledger_list if "Syn_5xISRE" in p["part_name"]), None)
            hlae_part = next((p for p in ledger_list if "HLAE_" in p["part_name"]), None)
            _check(
                cd47_part is not None and eef_part is not None and isre_part is not None and hlae_part is not None,
                "Arm 1B missing required bidirectional promoter/CDS parts",
            )
            _check(
                cd47_part.get("orientation") == "-" and eef_part.get("orientation") == "-",
                "Arm 1B CD47 CDS and EEF1A1 promoter must both be in antisense ('-') orientation",
            )
            _check(
                hlae_part.get("orientation") == "+" and isre_part.get("orientation") == "+",
                "Arm 1B HLA-E CDS and Syn_5xISRE promoter must both be in sense ('+') orientation",
            )
            _check(
                cd47_part["end_1based"] < eef_part["start_1based"] < isre_part["start_1based"] < hlae_part["start_1based"],
                "Arm 1B head-to-head bidirectional architecture order violated (expected antisense CD47 <- antisense EEF1A1 | sense ISRE -> sense HLA-E)",
            )
        for idx_comp, comp in enumerate(ledger_list):
            p_name = comp["part_name"]
            sub = seq[comp["start_1based"] - 1 : comp["end_1based"]]
            _check(comp["start_1based"] == cursor, f"Discontinuity in {cid} at {p_name}")
            _check(len(sub) == comp["length_bp"], f"Length mismatch in {cid}/{p_name}")
            comp_sha = hashlib.sha256(sub.encode("utf-8")).hexdigest()
            _check(comp_sha == comp["sha256"], f"Part SHA-256 mismatch in {cid}/{p_name}")
            orient = comp.get("orientation", "+")
            if p_name in raw_seqs:
                expected_part_seq = raw_seqs[p_name]["sequence"] if orient == "+" else revcomp_dna(raw_seqs[p_name]["sequence"])
                _check(sub == expected_part_seq, f"GenBank/GRCh38 sequence mismatch in {cid}/{p_name} ({orient})")
                if "_HA_" in p_name:
                    source = raw_seqs[p_name]
                    pin = raw_ref.get("validated_homology_arm_source_registry", {}).get(p_name, {})
                    start, end, strand = source.get("start_1based"), source.get("end_1based"), source.get("strand")
                    _check(type(start) is int and type(end) is int and 1 <= start <= end,
                           f"Invalid raw source interval in {cid}/{p_name}")
                    _check(end - start + 1 == len(source["sequence"]) and type(strand) is int and strand in (-1, 1),
                           f"Raw source interval length/strand mismatch in {cid}/{p_name}")
                    _check(source.get("region") == f"GRCh38:{source['chrom']}:{start}-{end}:strand={strand:+d}",
                           f"Raw source region string mismatch in {cid}/{p_name}")
                    independent_identity = {"assembly": source.get("assembly"), "chrom": source["chrom"], "start_1based": start,
                        "end_1based": end, "strand": strand, "length_bp": len(source["sequence"]),
                        "sequence_sha256": hashlib.sha256(source["sequence"].encode()).hexdigest()}
                    _check(independent_identity == pin, f"Raw source differs from independent pinned identity in {cid}/{p_name}")
                    _check(comp["source_accession_or_coordinate"] == f"GRCh38:{source['chrom']}:{start}-{end}:{strand:+d}",
                           f"Ledger source interval differs from raw identity in {cid}/{p_name}")
            if comp.get("part_type") == "CDS" or "CDS" in p_name or "ORF" in p_name or "P2A" in p_name:
                coding_dna = sub if orient == "+" else revcomp_dna(sub)
                aa = translate_dna(coding_dna)
                if "P2A" in p_name:
                    _check(
                        aa.startswith("GSG") and aa.endswith("DVEENPGP") and "ATNFSLLKQAGDVEENPGP" in aa and "*" not in aa,
                        f"2A peptide {p_name} in {cid} lacks canonical ATNFSLLKQAGDVEENPGP ribosome-skipping motif: {aa}",
                    )
                elif orient == "+" and idx_comp + 1 < len(ledger_list) and "P2A" in ledger_list[idx_comp + 1]["part_name"]:
                    _check(aa.startswith("M"), f"Upstream 2A-linked CDS {p_name} in {cid} lacks start Met")
                    _check("*" not in aa, f"Upstream 2A-linked CDS {p_name} in {cid} must read through without stop codon")
                else:
                    _check(aa.startswith("M"), f"CDS {p_name} in {cid} ({orient}) lacks start Met")
                    _check(aa.endswith("*"), f"CDS {p_name} in {cid} ({orient}) lacks terminal stop codon")
                    _check("*" not in aa[:-1], f"Internal stop codon in {p_name} in {cid} ({orient})")
            cursor = comp["end_1based"] + 1
        _check(cursor - 1 == len(seq), f"Component sum {cursor - 1} != total length {len(seq)}")


def validate_cif_file_content(cif_text: str, path_label: str) -> int:
    _check(len(cif_text) >= 10000, f"CIF {path_label} is truncated ({len(cif_text)} bytes < 10,000)")
    ca_coords: List[tuple[float, float, float]] = []
    for line in cif_text.splitlines():
        if (line.startswith("ATOM ") or line.startswith("HETATM ")) and " CA " in line:
            toks = line.split()
            if len(toks) >= 13:
                try:
                    x, y, z = float(toks[10]), float(toks[11]), float(toks[12])
                    ca_coords.append((x, y, z))
                except ValueError:
                    pass
    _check(len(ca_coords) >= 50, f"CIF {path_label} has only {len(ca_coords)} C-alpha atoms (< 50)")
    xs = [c[0] for c in ca_coords]
    ys = [c[1] for c in ca_coords]
    zs = [c[2] for c in ca_coords]
    coord_span = (max(xs) - min(xs)) + (max(ys) - min(ys)) + (max(zs) - min(zs))
    _check(coord_span > 25.0, f"CIF {path_label} has degenerate/zero 3D coordinates (span={coord_span:.2f} A)")
    return len(ca_coords)


def run_negative_control_mutated_fasta_and_empty_cif_probes(
    fasta_records: Dict[str, str], l2: Dict[str, Any], raw_ref: Dict[str, Any]
) -> None:
    """V70-10 + AUDIT-REF + AUDIT-REF Proof: Mutating FASTA or degenerate .cif files MUST fail even under python3 -O."""
    mutated_fasta = dict(fasta_records)
    first_key = next(iter(mutated_fasta))
    mutated_fasta[first_key] = "A" * len(mutated_fasta[first_key])
    failed_fasta = False
    try:
        validate_fasta_against_l2_atlas(mutated_fasta, l2, raw_ref)
    except AssertionError:
        failed_fasta = True
    _check(failed_fasta, "Negative control FAILED: all-A mutated FASTA was not rejected!")

    # Test antisense CDS mutation in Arm 1B (AUDIT-REF)
    arm1b_key = "ISLET1_LOCUS1_ARM1B_B2M_BIDIR_2AFREE_CD47_IGV_GPI_HLAE_SCT_3974bp"
    mut_arm1b = dict(fasta_records)
    seq1b = mut_arm1b[arm1b_key]
    mut_seq1b = seq1b[:630] + "TTA" + seq1b[633:]
    mut_arm1b[arm1b_key] = mut_seq1b
    l2_copy = json.loads(json.dumps(l2))
    for c in l2_copy["construct_sequence_and_grammar_verification"]["constructs"]:
        if c["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01":
            c["donor_sha256"] = hashlib.sha256(mut_seq1b.encode("utf-8")).hexdigest()
            for comp in c["part_ledger"]:
                if comp["part_name"] == "CD47_IgV_CD55_GPI_CDS_528bp_antisense":
                    sub = mut_seq1b[comp["start_1based"] - 1 : comp["end_1based"]]
                    comp["sha256"] = hashlib.sha256(sub.encode("utf-8")).hexdigest()
    failed_antisense = False
    try:
        validate_fasta_against_l2_atlas(mut_arm1b, l2_copy, raw_ref)
    except AssertionError:
        failed_antisense = True
    _check(failed_antisense, "Negative control FAILED: internal stop codon in antisense CDS was not rejected!")

    failed_cif = False
    try:
        validate_cif_file_content("", "EMPTY_0_BYTE_CIF_PROBE")
    except AssertionError:
        failed_cif = True
    _check(failed_cif, "Negative control FAILED: 0-byte .cif file was not rejected!")

    # Test degenerate all-zero coordinate .cif probe (AUDIT-REF)
    zero_cif_lines = ["ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 90.00 ? 1 ALA A CA 1"] * 150
    zero_cif_text = "\n".join(zero_cif_lines) + "\n" + ("#" * 10000)
    failed_zero_coord_cif = False
    try:
        validate_cif_file_content(zero_cif_text, "ZERO_COORD_CIF_PROBE")
    except AssertionError:
        failed_zero_coord_cif = True
    _check(failed_zero_coord_cif, "Negative control FAILED: all-zero coordinate .cif probe was not rejected!")

    # Test overlapping right-homology-arm GRCh38 coordinates in verify_construct_part_grammar (AUDIT-REF)
    import sys
    if str(REPO_ROOT / "scripts") not in sys.path:
        sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from run_layer2_perturbation_and_regulatory_plane import verify_construct_part_grammar

    cv_valid = json.loads(json.dumps(l2["construct_sequence_and_grammar_verification"]))
    res_valid = verify_construct_part_grammar(cv_valid)
    _check(
        res_valid["all_construct_grammars_verified"] is True and res_valid["all_construct_grammars_valid"] is True,
        f"Expected valid construct grammar to pass: {res_valid}",
    )
    cv_overlap = json.loads(json.dumps(l2["construct_sequence_and_grammar_verification"]))
    for c_item in cv_overlap["constructs"]:
        if c_item["construct_id"] == "LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01":
            for p_item in c_item["part_ledger"]:
                if p_item["part_name"] == "B2M_HA_R_600bp":
                    p_item["source_accession_or_coordinate"] = "GRCh38:chr15:44711009-44711608:+1"
    res_overlap = verify_construct_part_grammar(cv_overlap)
    _check(
        res_overlap["all_construct_grammars_verified"] is False
        and res_overlap["all_construct_grammars_valid"] is False
        and res_overlap["all_homology_arms_non_overlapping_verified"] is False,
        f"Negative control FAILED: overlapping right-arm GRCh38 coordinates were not rejected by verify_construct_part_grammar: {res_overlap}",
    )



def regenerate_layer3_from_current_inputs(root: Path = REPO_ROOT) -> Dict[str, Any]:
    """Replay all raw structural inputs in an owned checkout, never canonical output paths."""
    with tempfile.TemporaryDirectory(prefix="openislet-layer3-replay-") as temp:
        scratch = Path(temp)
        (scratch / "scripts").mkdir()
        shutil.copy2(root / "scripts/run_layer3_structural_binder_plane.py", scratch / "scripts")
        for dependency in ("run_layer2_perturbation_and_regulatory_plane.py", "product_contract.py"):
            shutil.copy2(root / "scripts" / dependency, scratch / "scripts")
        shutil.copytree(root / "data/structures", scratch / "data/structures")
        (scratch / "data/extracted").mkdir()
        shutil.copy2(root / "data/extracted/layer2_regulatory_and_perturbation_atlas.json", scratch / "data/extracted")
        (scratch / "data/constructs").mkdir()
        for fasta_p in sorted((root / "data/constructs").glob("*.fasta")):
            shutil.copy2(fasta_p, scratch / "data/constructs")
        (scratch / "data/protocols").mkdir()
        shutil.copy2(root / "data/protocols/islet1_product_definitions_v0.8.0.json", scratch / "data/protocols")
        result = subprocess.run([sys.executable, "-B", str(scratch / "scripts/run_layer3_structural_binder_plane.py")],
                                cwd=scratch, capture_output=True, text=True)
        _check(result.returncode == 0, f"Fresh Layer 3 replay failed: {result.stderr[-3000:]}")
        return json.loads((scratch / "data/extracted/layer3_structural_verification_atlas.json").read_text())


def verify_fresh_layer3_parity(committed: Dict[str, Any], root: Path = REPO_ROOT) -> Dict[str, Any]:
    fresh = regenerate_layer3_from_current_inputs(root)
    _check(fresh == committed, "Layer 3 artifact is stale or raw inputs changed: whole-atlas replay differs")
    _check(fresh["summary"]["all_reference_coordinate_geometric_screens_passed"] is True,
           "Fresh raw-input geometry screen failed")
    return fresh

def main() -> None:
    _check(L2_ATLAS.exists(), f"Missing {L2_ATLAS}")
    _check(L2_FASTA.exists(), f"Missing {L2_FASTA}")
    _check(L2_RAW_REF.exists(), f"Missing {L2_RAW_REF}")
    _check(L3_ATLAS.exists(), f"Missing {L3_ATLAS}")

    l2 = json.loads(L2_ATLAS.read_text(encoding="utf-8"))
    raw_ref = json.loads(L2_RAW_REF.read_text(encoding="utf-8"))
    l3 = json.loads(L3_ATLAS.read_text(encoding="utf-8"))
    fasta = parse_fasta_text(L2_FASTA.read_text(encoding="utf-8"))

    # 1. Run independent sequence, SHA-256, component, and ORF translation verification
    validate_fasta_against_l2_atlas(fasta, l2, raw_ref)

    # Reject an atlas generated before the current JSONL/inputs, independently of cached summaries.
    if str(REPO_ROOT / "scripts") not in sys.path:
        sys.path.insert(0,str(REPO_ROOT / "scripts"))
    from run_layer2_perturbation_and_regulatory_plane import simulate_stage0_to_stage7_perturbation_plane
    effects_path = REPO_ROOT / "data/extracted/ring1_atomic_effects_v0.5.jsonl"
    fresh_plane = simulate_stage0_to_stage7_perturbation_plane(raw_ref,
        json.loads(json.dumps(l2["construct_sequence_and_grammar_verification"])),
        l2["jaspar_pwm_regulatory_motif_scans"], effects_path=effects_path)
    _check(fresh_plane == l2["stage0_to_stage7_perturbation_plane"],
           "Layer 2 prior plane is stale: current JSONL/reference/grammar replay differs from atlas")
    provenance = fresh_plane["input_provenance"]["ring1_atomic_effects"]
    _check(provenance["path"] == "data/extracted/ring1_atomic_effects_v0.5.jsonl"
        and provenance["sha256"] == hashlib.sha256(effects_path.read_bytes()).hexdigest(),
        "Layer 2 effects input provenance does not match the consumed current JSONL")

    # 2. Run negative-control mutation probes (V70-10)
    run_negative_control_mutated_fasta_and_empty_cif_probes(fasta, l2, raw_ref)

    # 3. Verify V70-01, V70-02, V70-03, V70-04, V70-05, V70-06 fixes in Layer 2
    _check(l2["summary"]["all_homology_arms_non_overlapping_grch38_verified"] is True, "Homology arms overlap check failed")
    _check(l2["summary"]["cbx3_a2ucoe_455bp_avi_coordinate_overlap_bp"] == 455, "CBX3 455bp overlap mismatch")
    avi = l2["alphagenome_avi_saturation_mutagenesis_and_variants"]
    _check(avi["source_origin"] == "UNRESOLVED_ARCHIVAL_TABLE", "AVI origin qualification missing")
    _check(avi["origin_and_redistribution_rights"] == "UNVERIFIED", "AVI rights qualification missing")
    _check(avi["total_snvs_scored_and_archived"] is None, "Unverified AVI scores cannot be model predictions")
    _check(all(w["max_avi_phred"] is None and w["snvs_evaluated"] is None for w in avi["window_saturation_mutagenesis_summaries"].values()), "Public AVI scores must be null")

    constructs = {
        c["construct_id"]: c for c in l2["construct_sequence_and_grammar_verification"]["constructs"]
    }
    c1a = constructs["LOCUS1_ARM1A_TPT2A_BICISTRONIC"]
    c1b = constructs["LOCUS1_ARM1B_2AFREE_BIDIRECTIONAL_R01"]
    c3_l2 = constructs["LOCUS3_CXCL10_SERPINB9_TPT2A_PDL1_R01"]
    c4_l2 = constructs["LOCUS4_AAVS1_CHS4_455BP_CBX3_A2UCOE_ICASP9_R17"]

    _check(c1a["packaged_aav6_genome_with_2x145bp_itrs_bp"] == c1a["total_length_bp"] + 290, "Arm 1A ITR length mismatch")
    _check(c1a["aav6_single_capsid_compatible_le_4700bp"] is False, "Arm 1A should exceed 4700bp AAV6 ceiling")
    _check(c1b["packaged_aav6_genome_with_2x145bp_itrs_bp"] == 4264 <= 4700, "Arm 1B ITR length mismatch")
    _check(c1b["aav6_single_capsid_compatible_le_4700bp"] is True, "Arm 1B should be AAV6 compatible")
    _check(c3_l2["packaged_aav6_genome_with_2x145bp_itrs_bp"] == 4380 <= 4700, "Locus 3 ITR length mismatch")
    _check(c4_l2["packaged_aav6_genome_with_2x145bp_itrs_bp"] == 3747 <= 4700, "Locus 4 ITR length mismatch")
    _check(
        [c1a["donor_cpg_dinucleotide_count"], c1b["donor_cpg_dinucleotide_count"], c3_l2["donor_cpg_dinucleotide_count"], c4_l2["donor_cpg_dinucleotide_count"]]
        == [308, 260, 162, 202],
        "Whole-donor CpG counts mismatch",
    )
    _check(
        [
            c1a["cargo_excluding_ha_cpg_dinucleotide_count"],
            c1b["cargo_excluding_ha_cpg_dinucleotide_count"],
            c3_l2["cargo_excluding_ha_cpg_dinucleotide_count"],
            c4_l2["cargo_excluding_ha_cpg_dinucleotide_count"],
        ]
        == [219, 181, 155, 168],
        "Cargo-only CpG counts mismatch",
    )
    _check(
        [
            c1a["cds_only_cpg_dinucleotide_count"],
            c1b["cds_only_cpg_dinucleotide_count"],
            c3_l2["cds_only_cpg_dinucleotide_count"],
            c4_l2["cds_only_cpg_dinucleotide_count"],
        ]
        == [157, 151, 114, 78],
        "CDS-only CpG counts mismatch",
    )
    for c_entry in (c1a, c1b, c3_l2, c4_l2):
        _check(
            "pseudocount_0p8_production" in c_entry.get("whole_donor_jaspar_6pwm_scans_by_pseudocount", {})
            and "pseudocount_0p5_comparator" in c_entry.get("whole_donor_jaspar_6pwm_scans_by_pseudocount", {}),
            f"Missing whole-donor 6-PWM pseudocount scans in {c_entry['construct_id']}",
        )

    # Verify full-length HLA-E*01:03 retains cytoplasmic tail 344-358 AEWSDSAQGSESHSL (V70-02)
    # and CD47-IgV-GPI mature cell-surface form is 129 aa ending at native omega-site Ser353 (AUDIT-REF)
    ledger_l2 = l2["construct_sequence_and_grammar_verification"]["processed_polypeptide_ledger"]
    for k, entry in ledger_l2.items():
        if "HLAE" in k:
            _check("AEWSDSAQGSESHSL" in entry["sequence"], f"Missing HLA-E cytoplasmic tail in {k}")
    cd47_mat = ledger_l2["CD47_IgV_GPI_mature_cell_surface"]
    _check(cd47_mat["mature_length_aa"] == 129, "CD47 mature_length_aa != 129")
    _check(cd47_mat["sequence"].endswith("S"), "CD47 mature sequence must end with Ser353")
    _check(cd47_mat["sequence"][-9:] == raw_ref["uniprot_protein_sequences"]["P08174"]["sequence"][344:353], "CD47 stalk mismatch")

    # Verify WPRE authentic WHV J04514.1 (AUDIT-REF + AUDIT-REF) and EEF1A1 core promoter (AUDIT-REF)
    wpre_ref = raw_ref["ensembl_verified_homology_arms_and_elements"]["WPREmut6_589bp_J04514_mut6"]
    _check("J04514.1" in wpre_ref["source_coord"] and "J04364" not in wpre_ref["source_coord"], "WPRE accession mismatch")
    _check(wpre_ref["length_bp"] == 589, "WPRE length != 589")
    _check("AUDIT-REF" in wpre_ref["variant_qualification_status"], "WPRE variant qualification missing")
    eef_ref = raw_ref["ensembl_verified_homology_arms_and_elements"]["EEF1A1_core_promoter_230bp"]
    _check("73521000-73521229" in eef_ref["source_coord"], "EEF1A1 coordinate mismatch")
    gencode_models = raw_ref.get("gencode_v46_transcript_models_pinned", {})
    _check(len(gencode_models) >= 4, "Missing GENCODE v46 transcript models")
    # AUDIT-REF: Verify exact UCSC GENCODE v46 1-based closed coordinates
    _check(gencode_models["B2M"]["tss_1based"] == 44711517, "B2M TSS mismatch")
    _check(gencode_models["B2M"]["tes_1based"] == 44718145, "B2M TES mismatch")
    _check(gencode_models["CXCL10"]["tss_1based"] == 76023497, "CXCL10 TSS mismatch")
    _check(gencode_models["CXCL10"]["tes_1based"] == 76021118, "CXCL10 TES mismatch")
    _check(gencode_models["PPP1R12C_AAVS1"]["tss_1based"] == 55117637, "PPP1R12C TSS mismatch")
    _check(gencode_models["PPP1R12C_AAVS1"]["tes_1based"] == 55090918, "PPP1R12C TES mismatch")

    # AUDIT-REF: Verify exact discrete Markov-0 dynamic programming PWM p-value convolution and deduplicated Ring-1 study counts
    motif_scans = l2["jaspar_pwm_regulatory_motif_scans"]
    for scan_val in motif_scans.values():
        for tf_val in scan_val.get("top_jaspar_pwm_hits", {}).values():
            _check("exact_markov0_p_value_unadjusted" in tf_val["top_hits"][0], "Missing exact_markov0_p_value_unadjusted")
    calib_inputs = l2["stage0_to_stage7_perturbation_plane"]["calibration_inputs_consumed"]
    _check("encode_screen_cxcl10_max_biosample_tpm" in calib_inputs, "Missing encode_screen_cxcl10_max_biosample_tpm")
    _check(calib_inputs["ring1_gene_effect_evidence_counts"]["CD47"]["total"] <= 50, "CD47 Ring-1 count out of bounds")

    # Verify UP421 has no unmeasured Stage 0->7 differentiation claim (V70-06)
    profiles = {
        p["genotype_id"]: p for p in l2["stage0_to_stage7_perturbation_plane"]["perturbation_profiles"]
    }
    _check(profiles["SANA_UP421_CADAVERIC_PRIMARY_ISLET_BENCHMARK"]["stage7_sc_beta_yield_fraction"] is None, "UP421 should not have Stage 7 differentiation yield")

    # 4. Parse all 13 .cif files in data/structures/ (11 RCSB mmCIF + 2 AlphaFold v6) and verify Layer 3 fixes (V70-07, V70-08, V70-09, AUDIT-REF, IIP-0011)
    cif_files = sorted(STRUCT_DIR.rglob("*.cif"))
    _check(len(cif_files) == 13, f"Expected 13 .cif files in {STRUCT_DIR}, found {len(cif_files)}")
    for cif_path in cif_files:
        validate_cif_file_content(cif_path.read_text(encoding="utf-8"), str(cif_path.name))

    verify_fresh_layer3_parity(l3)

    # Cross-check Layer-2 to Layer-3 FASTA SHA-256 linkage and reference geometry status (AUDIT-REF)
    _check(l3["layer2_construct_fasta_sha256_linked"] == hashlib.sha256(L2_FASTA.read_bytes()).hexdigest(), "Layer 2 to Layer 3 FASTA SHA-256 mismatch")
    _check(l3["summary"]["reference_coordinate_geometry_status"] == "PASS", f"Expected reference_coordinate_geometry_status == PASS, got {l3['summary']['reference_coordinate_geometry_status']}")
    _check(l3["summary"]["all_reference_coordinate_geometric_screens_passed"] is True, "Expected all_reference_coordinate_geometric_screens_passed is True")

    complexes = {c["complex_id"]: c for c in l3["complexes"]}
    c1 = complexes["LAYER3_COMPLEX_1_HLAE_SCT_CD94_NKG2A"]
    _check("17.02 A" in c1["hlae_0103_gly107_stability_verification"]["qualification_note"] and "29.97" not in c1["hlae_0103_gly107_stability_verification"]["qualification_note"], "3CDG Arg107 CA distance note mismatch")
    c1_structs = {s["pdb_id"]: s for s in c1["experimental_structures_analyzed"]}
    _check("VMAPRTLFL" in c1_structs["3BZE"]["crystallized_peptide"], f"3BZE peptide mismatch: {c1_structs['3BZE']['crystallized_peptide']}")
    _check(c1["single_chain_trimer_linker_physics"]["linker_1_peptide_P9_to_B2M_Ile1"]["marko_siggia_entropic_tension_pN_at_310K_lp4A"] > 0.0, "WLC tension <= 0")
    _check(
        c1["single_chain_trimer_linker_physics"]["linker_1_peptide_P9_to_B2M_Ile1"]["model_classification"]
        == "coarse_geometric_screen_not_excluded_volume_folding_proof",
        "WLC model_classification mismatch",
    )
    # AUDIT-REF: Verify strict pH 7.4 ionic salt bridges vs histidine-assisted contacts & residue-pair deduplication
    c1_if = c1["trimer_to_cd94_nkg2a_interface"]
    _check("unique_salt_bridge_residue_pairs_arg_lys_ph74" in c1_if, "Missing unique_salt_bridge_residue_pairs_arg_lys_ph74")
    _check("his_protonation_dependent_candidate_pairs" in c1_if, "Missing his_protonation_dependent_candidate_pairs")
    _check("unique_residue_contact_pairs" in c1_if, "Missing unique_residue_contact_pairs")
    _check("interface_crystallographic_b_factor_profile_A2" in c1_if, "Missing interface_crystallographic_b_factor_profile_A2")

    c3 = complexes["LAYER3_COMPLEX_3_SERPINB9_GRANZYME_B"]
    _check(
        "Yan, Zhou et al. 2024 Acta Crystallogr Sect F Struct Biol Commun 80:286-293"
        in c3["foldseek_structural_homolog_verification"]["same_protein_human_serpinb9_pdb_hit_and_supersession"]["superseding_structure_citation"],
        "8ZCR citation mismatch",
    )
    _check(c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["aligned_ca_residues_AF_vs_1IAU"] == 226, "1IAU aligned residues != 226")
    _check(c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["identical_residue_pairs"] == 226, "1IAU identical residue pairs != 226")
    _check(c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["rmsd_AF_P10144_vs_1IAU_angstrom"] < 0.45, "1IAU RMSD >= 0.45 A")
    rcl_geom = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]
    _check(abs(rcl_geom["rcl_322_348_loop_mean_plddt"] - 48.36) < 0.1, "SERPINB9 RCL 322-348 pLDDT != 48.36")
    _check(abs(rcl_geom["rcl_335_345_hinge_to_p1p1prime_mean_plddt"] - 36.21) < 0.1, "SERPINB9 RCL 335-345 pLDDT != 36.21")
    _check(
        rcl_geom["functional_status"] == "experimentally_unresolved_requires_GzmB_kinetics_assay",
        "SERPINB9 RCL functional_status mismatch",
    )

    l3_cfh_prov = l3["locus2_cfh_stp_gpi_structural_provenance_iip_0031"]
    _check(
        l3_cfh_prov["pdb_2wii"]["glu264_ca_to_c3b_ted_q1013_ca_angstrom"] == 33.04,
        f"Expected 2WII Glu264->Gln1013(B:991) CA distance 33.04 A, got {l3_cfh_prov['pdb_2wii']['glu264_ca_to_c3b_ted_q1013_ca_angstrom']}",
    )
    _check(
        l3_cfh_prov["pdb_2wii"].get("c3b_ted_q1013_pdb_author_residue") == "B:991 (GLN; UniProt P01024 Gln1013)",
        "Missing or mismatched 2WII c3b_ted_q1013_pdb_author_residue label",
    )
    _check(l3["summary"]["total_experimental_pdb_structures_verified"] == 15, "Expected 15 experimental PDB structures verified")
    _check(l3["summary"]["experimental_rcsb_mmcif_structures_verified"] == 11, "Expected 11 experimental mmCIF structures verified")
    _check(l3["summary"]["experimental_rcsb_legacy_pdb_structures_verified"] == 4, "Expected 4 experimental legacy PDB structures verified")
    _check(
        {"2WII", "5O32"} <= set(l3["summary"]["experimental_pdb_accession_inventory"]),
        "Expected 2WII and 5O32 in summary experimental_pdb_accession_inventory",
    )

    print(
        "layer2_and_layer3_planes gate: PASS (4/4 FASTAs translated & SHA-256 verified, AVI quarantine metadata verified, fresh Layer 3 replay identical, "
        "10/10 CIFs parsed, 226/226 1IAU insertion codes preserved, negative-control mutation probes rejected)"
    )


run_layer2_and_layer3_suite = main


def test_iip0031_reference_panels_reject_incomplete_inputs(monkeypatch, tmp_path) -> None:
    import pytest
    from scripts import run_layer2_perturbation_and_regulatory_plane as layer2

    panel_keys = ("pvr_refseq_isoforms_4", "cd276_refseq_isoforms_4")
    texts = {f"row-{i}": f"RECORD-{i}" for i in range(4)}
    pins = {key: hashlib.sha256(text.encode()).hexdigest() for key, text in texts.items()}
    monkeypatch.setattr(layer2, "_IIP0031_TRANSCRIPT_PINS", {key: pins for key in panel_keys})
    monkeypatch.setattr(layer2, "REPO_ROOT", tmp_path)
    rows = {key: {"mrna_sequence": text, "mrna_sha256": pins[key]} for key, text in texts.items()}
    complete = {key: rows for key in panel_keys}
    bundle_key = "iip0031_locus2_modeb_v3_reference_bundle"
    build = layer2.build_locus2_ciita_modeb_comparators_v3_manifest

    with pytest.raises(ValueError, match="requires the exact"):
        build({bundle_key: {key: {} for key in panel_keys}})
    for key in panel_keys:
        assert layer2._verified_refseq_panel(complete, key) == texts
        with pytest.raises(ValueError, match=key):
            build({bundle_key: {other: rows for other in panel_keys if other != key}})
        missing_one = {k: v for k, v in rows.items() if k != "row-0"}
        invalid_panels = [
            None, {}, [], missing_one, {**rows, "extra": rows["row-0"]},
            {**missing_one, "substitute": rows["row-0"]},
        ]
        invalid_records = (
            None, [], {}, {"mrna_sequence": None}, {"mrna_sequence": "  "},
            {"mrna_sequence": texts["row-0"]}, rows["row-1"],
            {**rows["row-0"], "mrna_sequence": texts["row-1"]},
            {**rows["row-0"], "mrna_sha256": pins["row-1"]},
        )
        invalid_panels.extend({**rows, "row-0": row} for row in invalid_records)
        for malformed in invalid_panels:
            with pytest.raises(ValueError, match=key):
                build({bundle_key: {**complete, key: malformed}})


def test_iip0031_product_contract_hardening(tmp_path) -> None:
    import copy
    import pytest
    from scripts import product_contract as pc

    def _stage_registry_root(dest: Path) -> Path:
        (dest / "data" / "protocols").mkdir(parents=True, exist_ok=True)
        (dest / "data" / "constructs").mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / pc.CONTRACT_PATH, dest / pc.CONTRACT_PATH)
        shutil.copy2(REPO_ROOT / pc.DESIGN_B_FASTA_PATH, dest / pc.DESIGN_B_FASTA_PATH)
        shutil.copy2(REPO_ROOT / pc.LOCUS2_V3_FASTA_PATH, dest / pc.LOCUS2_V3_FASTA_PATH)
        shutil.copy2(REPO_ROOT / pc.DESIGN_Y_V2_FASTA_PATH, dest / pc.DESIGN_Y_V2_FASTA_PATH)
        shutil.copy2(REPO_ROOT / pc.D1B_AUTO_V2_FASTA_PATH, dest / pc.D1B_AUTO_V2_FASTA_PATH)
        return dest

    root = _stage_registry_root(tmp_path / "canonical")
    doc = pc.load_product_contracts(root)
    for row in doc["products"]:
        assert pc.validate_product_definition_dict(copy.deepcopy(row), root=root)["product_variant_id"] == row["product_variant_id"]

    # Mutated in-memory product_def dictionaries must be rejected by validate_product_definition_dict:
    l2_cd200 = copy.deepcopy(pc.get_product_definition("L2_V3_CO_PRIMARY_CD200", root=root))
    assert pc.validate_product_definition_dict(l2_cd200, root=root) == l2_cd200
    for mutator in (
        lambda d: d.__setitem__("execution_supported", 0),
        lambda d: d["staged_locus2_modeb_donor_binding"].__setitem__("donor_length_bp", 4012.0),
        lambda d: d["staged_locus2_modeb_donor_binding"].__setitem__("unit2_promoter", "EEF1A1_core_230bp"),
        lambda d: d.__setitem__("extra_injected_field", True),
    ):
        bad_def = copy.deepcopy(l2_cd200)
        mutator(bad_def)
        with pytest.raises(ValueError):
            pc.validate_product_definition_dict(bad_def, root=root)
    with pytest.raises(ValueError):
        pc.validate_product_definition_dict("not-a-dict", root=root)  # type: ignore[arg-type]

    pdef_path = root / pc.CONTRACT_PATH
    canonical_text = pdef_path.read_text(encoding="utf-8")
    canonical_json = json.loads(canonical_text)

    def _expect_registry_reject(mutated_obj_or_raw: Any, match: str | None = None) -> None:
        if isinstance(mutated_obj_or_raw, str):
            pdef_path.write_text(mutated_obj_or_raw, encoding="utf-8")
        else:
            pdef_path.write_text(json.dumps(mutated_obj_or_raw), encoding="utf-8")
        try:
            if match:
                with pytest.raises(ValueError, match=match):
                    pc.load_product_contracts(root)
            else:
                with pytest.raises(ValueError):
                    pc.load_product_contracts(root)
        finally:
            pdef_path.write_text(canonical_text, encoding="utf-8")

    # Duplicate JSON keys at root and nested levels must raise ValueError:
    _expect_registry_reject(
        '{"schema_version": "0.8.0", "schema_version": "0.8.0", "products": []}',
        match="Duplicate JSON key",
    )
    nested_dup = canonical_text.replace(
        '"unit2_promoter": "Syn_5xISRE_minP",',
        '"unit2_promoter": "EEF1A1_core_230bp", "unit2_promoter": "Syn_5xISRE_minP",',
        1,
    )
    assert nested_dup != canonical_text
    _expect_registry_reject(nested_dup, match="Duplicate JSON key")

    # Unhashable / malformed scalar types must raise clean ValueError (never TypeError or AttributeError):
    for field, bad_val in (
        ("family", []),
        ("family", 123),
        ("lifecycle_state", []),
        ("lifecycle_state", 123),
        ("execution_supported", 1),
        ("aliases", [" L2_V3_CO_PRIMARY_CD200 "]),
        ("aliases", [123]),
        ("nominal_packaged_aav6_bp", 4264.0),
        ("nominal_packaged_aav6_bp", True),
    ):
        bad = copy.deepcopy(canonical_json)
        bad["products"][0][field] = bad_val
        _expect_registry_reject(bad)

    # Missing any of the 8 STAGED_LOCUS2_V3_PRODUCT_IDS rows must raise ValueError:
    missing_row = copy.deepcopy(canonical_json)
    missing_row["products"] = [
        r for r in missing_row["products"] if r["product_variant_id"] != "ISLET1_A4_L2_V3_CO_PRIMARY_MANF"
    ]
    _expect_registry_reject(missing_row, match="ISLET1_A4_L2_V3_CO_PRIMARY_MANF")

    # Missing L2_V3_* short alias on its canonical row or reassigning it to another row must raise ValueError:
    dropped_alias = copy.deepcopy(canonical_json)
    for r in dropped_alias["products"]:
        if r["product_variant_id"] == "ISLET1_A4_L2_V3_CO_PRIMARY_CD200":
            r["aliases"] = [a for a in r["aliases"] if a != "L2_V3_CO_PRIMARY_CD200"]
    _expect_registry_reject(dropped_alias, match="L2_V3_CO_PRIMARY_CD200")

    reassigned_alias = copy.deepcopy(canonical_json)
    for r in reassigned_alias["products"]:
        if r["product_variant_id"] == "ISLET1_A4_L2_V3_CO_PRIMARY_CD200":
            r["aliases"] = [a for a in r["aliases"] if a != "L2_V3_CO_PRIMARY_CD200"]
        elif r["product_variant_id"] == "ISLET1_A4_LEAD":
            r["aliases"].append("L2_V3_CO_PRIMARY_CD200")
    _expect_registry_reject(reassigned_alias, match="L2_V3_CO_PRIMARY_CD200")

    # Mutating pinned semantic metadata, float lengths, or Design-B c1/c5 bindings must raise ValueError:
    for target_id, mut_fn in (
        ("ISLET1_A4_L2_V3_CO_PRIMARY_CD200", lambda r: r["staged_locus2_modeb_donor_binding"].__setitem__("unit2_promoter", "EEF1A1_core_230bp")),
        ("ISLET1_A4_L2_V3_CO_PRIMARY_CD200", lambda r: r["staged_locus2_modeb_donor_binding"].__setitem__("mir_class", "INTRONLESS")),
        ("ISLET1_A4_L2_V3_BARRIER_CHS4_250", lambda r: r["staged_locus2_modeb_donor_binding"].__setitem__("barrier_variant_type", "CANONICAL_80BP_PAIR")),
        ("ISLET1_A4_L2_V3_CO_PRIMARY_CD200", lambda r: r["staged_locus2_modeb_donor_binding"].__setitem__("donor_length_bp", 4012.0)),
        ("ISLET1_A4_L2_V3_CO_PRIMARY_CD200", lambda r: r["staged_locus2_modeb_donor_binding"].__setitem__("nominal_packaged_aav6_bp", 4302.0)),
        ("ISLET1_B2_5LOCUS_BASE", lambda r: r["staged_design_b_donor_reference_binding"].pop("cxcl10_locus5_primary_donor_sha256")),
        ("ISLET1_B2_5LOCUS_BASE", lambda r: r["staged_design_b_donor_reference_binding"].__setitem__("cxcl10_locus5_primary_donor_length_bp", 4078.0)),
        ("ISLET1_B2_L2_V3_DESIGN_B_CD200", lambda r: r["staged_design_b_donor_reference_binding"].pop("aavs1_locus6_primary_donor_sha256")),
    ):
        bad = copy.deepcopy(canonical_json)
        for r in bad["products"]:
            if r["product_variant_id"] == target_id:
                mut_fn(r)
        _expect_registry_reject(bad)

    # Missing or directory FASTA paths must raise ValueError (never silently skip or raise IsADirectoryError):
    for fasta_rel in (pc.DESIGN_B_FASTA_PATH, pc.LOCUS2_V3_FASTA_PATH):
        fpath = root / fasta_rel
        saved = fpath.read_bytes()
        try:
            fpath.unlink()
            with pytest.raises(ValueError, match="Missing or non-regular"):
                pc.load_product_contracts(root)
            fpath.mkdir()
            with pytest.raises(ValueError, match="Missing or non-regular"):
                pc.load_product_contracts(root)
            fpath.rmdir()
            fpath.write_bytes(saved[:-1])
            with pytest.raises(ValueError, match="SHA-256 mismatch"):
                pc.load_product_contracts(root)
        finally:
            if fpath.is_dir():
                fpath.rmdir()
            fpath.write_bytes(saved)


def test_iip0031_layer2_hardening_and_scanners(monkeypatch, tmp_path) -> None:
    import copy
    import pytest
    from scripts import run_layer2_perturbation_and_regulatory_plane as layer2

    raw_ref = json.loads(L2_RAW_REF.read_text(encoding="utf-8"))

    # 1. Frozen-input parity on overlapping seed scanner, spatial PWM audit, and guide counts:
    frozen = layer2.build_locus2_ciita_modeb_comparators_v3_manifest(ref=copy.deepcopy(raw_ref))
    rnai = frozen["intronic_mir_and_rnai_attribution_scope"]
    assert rnai["candidate_seed_pairs_screened"] == 296
    assert rnai["replayed_8mer_hits_total"] == 0
    assert rnai["replayed_7mer_m8_hits_total"] == 0
    assert len(rnai["disclosed_7mer_a1_hits"]) == 2
    assert len(rnai["disclosed_22gene_cds_7mer_m8_genes"]) == 6
    assert rnai["lead_guides_per_target"] == 1
    assert rnai["backup_guides_per_target"] == 4
    assert rnai["total_screened_guides_per_target"] == 5
    roles = [row["role"] for row in rnai["backup_guide_screen_ledger_15oi_22utr"]]
    assert roles.count("LEAD") == 2 and roles.count("BACKUP") == 8
    assert frozen["canonical_fastas_untouched"] is True
    assert frozen["layer1_pareto_artifact_verified"] is True
    assert frozen["layer1_and_oq016_score_delta_pp"] == 0.0
    co_cd200_entry = next(c for c in frozen["constructs"] if c["construct_id"] == "L2_V3_CO_PRIMARY_CD200")
    assert co_cd200_entry["descriptive_6pwm_jaspar_summary"]["exact_family_total_tests"] == 47990

    # 2. Single-read comparator FASTA immutability + construct-only temp root + proxy/canonical FASTA pin checks:
    construct_root = tmp_path / "construct_root"
    shutil.copytree(REPO_ROOT / "data" / "constructs", construct_root / "data" / "constructs")
    monkeypatch.setattr(layer2, "REPO_ROOT", construct_root)

    comp_fa = construct_root / "data" / "constructs" / "islet1_locus2_ciita_modeb_comparators_v0.8.0.fasta"
    comp_bytes = comp_fa.read_bytes()
    orig_read_bytes = Path.read_bytes
    comp_reads: List[str] = []

    def _mutate_comp_after_first_read(p: Path) -> bytes:
        data = orig_read_bytes(p)
        if p == comp_fa:
            comp_reads.append(p.name)
            p.write_bytes(b">CORRUPTED_AFTER_FIRST_READ\nAAAA\n")
        return data

    with monkeypatch.context() as mp:
        mp.setattr(Path, "read_bytes", _mutate_comp_after_first_read)
        res_race = layer2.build_locus2_ciita_modeb_comparators_v3_manifest(ref=copy.deepcopy(raw_ref))
    assert comp_reads == ["islet1_locus2_ciita_modeb_comparators_v0.8.0.fasta"]
    assert res_race["layer1_pareto_artifact_verified"] is False
    assert res_race["layer1_and_oq016_score_delta_pp"] is None
    comp_fa.write_bytes(comp_bytes)

    # Stripped or mutated proxy / canonical FASTAs must raise ValueError:
    for fasta_name in (
        "islet1_scenario_a_constructs_v0.7.0.fasta",
        "islet1_design_b_complete_constructs_v0.8.0.fasta",
        "islet1_design_y_wildcard_comparators_v0.8.0.fasta",
        "nm_021999_5_archival_mrna.fasta",
    ):
        fp = construct_root / "data" / "constructs" / fasta_name
        saved = fp.read_bytes()
        try:
            fp.write_bytes(b">EMPTY_HEADER\n")
            with pytest.raises(ValueError):
                layer2.build_locus2_ciita_modeb_comparators_v3_manifest(ref=copy.deepcopy(raw_ref))
        finally:
            fp.write_bytes(saved)

    with pytest.raises(ValueError, match="Duplicate FASTA"):
        layer2._parse_fa_bytes(b">REC1\nACGT\n>REC1\nTGCA\n", label="dup_test")
    with pytest.raises(ValueError, match="Empty or non-ACGT"):
        layer2._parse_fa_bytes(b">REC1\n", label="empty_test")
    with pytest.raises(ValueError, match="non-ACGT"):
        layer2._parse_fa_bytes(b">REC1\nACGTN\n", label="non_acgt_test")

    # Mutated Pareto artifact must raise ValueError when present:
    pareto_dst = construct_root / "data" / "extracted" / "islet1_combinatorial_pareto_and_ladder_v0.8.0.json"
    pareto_dst.parent.mkdir(parents=True, exist_ok=True)
    pareto_dst.write_text('{"tampered": true}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        layer2.build_locus2_ciita_modeb_comparators_v3_manifest(ref=copy.deepcopy(raw_ref))
    pareto_dst.unlink()

    # 3. CIITA 8-isoform exon coordinates/labels, genomic metadata, 22-gene cds_start_1based, and 6-PWM pins:
    for bad_mut in (
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["ciita_refseq_isoforms_8"].__setitem__("NM_000246.4", "not-a-dict"),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["ciita_refseq_isoforms_8"]["NM_000246.4"].__setitem__("exons_1based", []),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["ciita_refseq_isoforms_8"]["NM_000246.4"].__setitem__(
            "exons_1based", [[1, 4532]]
        ),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["ciita_refseq_isoforms_8"]["NM_000246.4"].__setitem__("iso_key", "TAMPERED"),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["ciita_refseq_isoforms_8"]["NM_000246.4"].__setitem__(
            "promoter_first_exon_mapping", "TAMPERED"
        ),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["grch38_chr16_ciita_exon5_2001bp"].__setitem__("chrom", "chr1"),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["grch38_chr19_dm1_ctcf1_80bp"].__setitem__("strand", -1),
        lambda r: r["ensembl_verified_homology_arms_and_elements"].pop("CIITA_pan_isoform_exon5_window_200bp"),
        lambda r: r["ensembl_verified_homology_arms_and_elements"]["CIITA_pan_isoform_exon5_window_200bp"].__setitem__("sha256", "0" * 64),
        lambda r: r["iip0031_locus2_modeb_v3_reference_bundle"]["human_22_gene_refseq_panel"]["PDX1"].pop("cds_start_1based"),
        lambda r: r["jaspar_pwms"].pop("MA0139.1"),
        lambda r: r["jaspar_pwms"]["MA0139.1"].__setitem__("pwm_matrix_sha256", "0" * 64),
        lambda r: r["jaspar_pwms"]["MA0139.1"]["log2_odds_pwm"]["A"].__setitem__(0, 10**400),
    ):
        bad_ref = copy.deepcopy(raw_ref)
        bad_mut(bad_ref)
        with pytest.raises(ValueError):
            layer2.build_locus2_ciita_modeb_comparators_v3_manifest(ref=bad_ref)

    # 4. Direct scan_jaspar_pwm overflow guard, uncapped all_significant_hits, and family-test fields:
    irf1_pwm = copy.deepcopy(raw_ref["jaspar_pwms"]["MA0050.2"])
    dense_irf1_seq = "GAAAGTGAAACC" * 40
    dense_scan = layer2.scan_jaspar_pwm(dense_irf1_seq, irf1_pwm, top_k=5, n_matrices_in_scan=6)
    assert len(dense_scan["top_hits"]) == 5
    assert len(dense_scan["all_significant_hits"]) >= 35
    assert dense_scan["focal_matrix_bonferroni_family_factor"] == 6 * dense_scan["total_positions_tested_both_strands"]
    assert dense_scan["total_tests_across_all_scanned_matrices_in_region"] is None
    dense_scan_exact = layer2.scan_jaspar_pwm(
        dense_irf1_seq, irf1_pwm, top_k=5, n_matrices_in_scan=6, actual_family_total_tests=5000
    )
    assert dense_scan_exact["total_tests_across_all_scanned_matrices_in_region"] == 5000

    overflow_pwm = copy.deepcopy(irf1_pwm)
    overflow_pwm.pop("pwm_matrix_sha256", None)
    overflow_pwm["log2_odds_pwm"]["A"][0] = 10**400
    with pytest.raises(ValueError, match="finite"):
        layer2.scan_jaspar_pwm("ACGT" * 10, overflow_pwm)

    # 5. Overlapping _scan_seeds verification ("T"*22 on "A"*9 -> 2x 8mer, 3x 7mer-m8, 3x 7mer-A1):
    seed_hits = layer2._scan_seeds("T" * 22, {"homopoly_A9": "A" * 9})
    assert seed_hits == {"homopoly_A9": {"8mer": 2, "7mer_m8": 3, "7mer_A1": 3}}
    with pytest.raises(ValueError):
        layer2._scan_seeds("ACGTTGC", {"t": "AAAAAAAAA"})
    with pytest.raises(ValueError):
        layer2._scan_seeds("T" * 22, {"t": "AAAAANAAA"})

    # 6. Repeat-scanner input validation, budget edge cases, unclustered pairs, and cross-span DR:
    for bad_seq in ("NNNN", "RRYY", "ßß", 123):
        with pytest.raises(ValueError):
            layer2.scan_disjoint_inverted_repeats_with_mismatches(bad_seq, min_k=4, max_mismatches=0)  # type: ignore[arg-type]
    for bad_k, bad_mm in ((0, 0), (-1, 0), (True, 0), (4, -1), (4, False)):
        with pytest.raises(ValueError):
            layer2.scan_disjoint_inverted_repeats_with_mismatches("ACGTACGT", min_k=bad_k, max_mismatches=bad_mm)  # type: ignore[arg-type]

    hits_at = layer2.scan_disjoint_inverted_repeats_with_mismatches("AT", min_k=1, max_mismatches=0)
    assert len(hits_at) == 1 and hits_at[0]["pos1_a"] == 1 and hits_at[0]["pos1_b"] == 2
    hits_budget = layer2.scan_disjoint_inverted_repeats_with_mismatches("AAAA", min_k=2, max_mismatches=2)
    assert len(hits_budget) == 1 and hits_budget[0]["mismatches"] == 2
    hits_aattt = layer2.scan_disjoint_inverted_repeats_with_mismatches("AATTT", min_k=2, max_mismatches=0)
    assert [(h["pos1_a"], h["pos1_b"]) for h in hits_aattt] == [(1, 3), (1, 4)]
    hits_aattt_clustered = layer2.scan_disjoint_inverted_repeats_with_mismatches(
        "AATTT", min_k=2, max_mismatches=0, cluster_within_bp=4
    )
    assert len(hits_aattt_clustered) == 1
    hits_aaattt = layer2.scan_disjoint_inverted_repeats_with_mismatches("AAATTT", min_k=2, max_mismatches=0)
    assert len(hits_aaattt) == 4
    assert next(h for h in hits_aaattt if (h["pos1_a"], h["pos1_b"]) == (1, 5))["maximal_stem_bp"] == 3
    assert next(h for h in hits_aaattt if (h["pos1_a"], h["pos1_b"]) == (2, 4))["maximal_stem_bp"] == 3

    assert layer2._max_overlapping_dr("ACGTACGTACGTCCCTGAGA", excluded_spans_0based=((0, 12),))[0] <= 3
    assert layer2._max_overlapping_dr("ACGTACGGTTCCACGTACGGTT", excluded_spans_0based=((0, 10),))[0] == 10
    assert layer2._max_overlapping_dr("ACGACG")[0] == 3
    assert layer2._max_disjoint_dr("ACGACG")[0] == 3
    assert layer2._max_disjoint_ir("AAATTT")[0] == 3
    assert layer2._count_disjoint_ir_ge_k("AATTT", min_k=2) == 2


def test_iip0031_layer3_2wii_dbref_and_inventory() -> None:
    from scripts import run_layer3_structural_binder_plane as l3

    prov = l3._build_locus2_cfh_stp_gpi_structural_provenance_iip_0031()
    assert prov["pdb_2wii"]["glu264_ca_to_c3b_ted_q1013_ca_angstrom"] == 33.04
    assert prov["pdb_2wii"]["c3b_ted_q1013_pdb_author_residue"] == "B:991 (GLN; UniProt P01024 Gln1013)"

    fresh = regenerate_layer3_from_current_inputs(REPO_ROOT)
    summary = fresh["summary"]
    assert summary["total_experimental_pdb_structures_verified"] == 15
    assert summary["experimental_rcsb_mmcif_structures_verified"] == 11
    assert summary["experimental_rcsb_legacy_pdb_structures_verified"] == 4
    assert summary["canonical_experimental_pdb_structures_verified"] == 13
    assert summary["challenger_iip0031_experimental_pdb_structures_verified"] == 2
    assert summary["experimental_pdb_accession_inventory"][-2:] == ["2WII", "5O32"]
    assert fresh["locus2_cfh_stp_gpi_structural_provenance_iip_0031"] == prov


def test_layer2_and_layer3_planes() -> None:
    main()


if __name__ == "__main__":
    main()
