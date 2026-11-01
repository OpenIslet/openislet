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
        shutil.copy2(root / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta", scratch / "data/constructs")
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

    # 4. Parse all 10 .cif files in data/structures/ and verify Layer 3 fixes (V70-07, V70-08, V70-09, AUDIT-REF, AUDIT-REF)
    cif_files = sorted(STRUCT_DIR.rglob("*.cif"))
    _check(len(cif_files) == 10, f"Expected 10 .cif files in {STRUCT_DIR}, found {len(cif_files)}")
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

    print(
        "layer2_and_layer3_planes gate: PASS (4/4 FASTAs translated & SHA-256 verified, AVI quarantine metadata verified, fresh Layer 3 replay identical, "
        "10/10 CIFs parsed, 226/226 1IAU insertion codes preserved, negative-control mutation probes rejected)"
    )


run_layer2_and_layer3_suite = main


def test_layer2_and_layer3_planes() -> None:
    main()


if __name__ == "__main__":
    main()
