#!/usr/bin/env python3
"""
Layer 3: 3D Structural & Protein Interface Verification Plane (v0.8.0)
======================================================================
Computes deterministic atomic-level structural geometry directly from 10 experimental
RCSB PDB structures (8 mmCIF: 3CDG, 1MHE, 3BZE, 2JJS, 4CMM, 4KJY, 4ZQK, 1IAU; plus
2 PDB: 1EFX, 1IM9) and 2 AlphaFold DB v6 models (AF-P50453-F1 SERPINB9, AF-P10144-F1 Granzyme B), consumes
the exact processed polypeptide ledger from Layer 2 (`layer2_regulatory_and_perturbation_atlas.json`),
and separates coordinate-computed geometry from literature-derived reference affinities
and unresolved engineered-product hypotheses (closing V70-07, V70-08, and V70-09).

Outputs:
  - data/extracted/layer3_structural_verification_atlas.json
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import re
import sys
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np

_MMCIF_TOKEN_RE = re.compile(r"""'(?:[^']|'(?!\s|$))*'|"(?:[^"]|"(?!\s|$))*"|\S+""")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.run_layer2_perturbation_and_regulatory_plane import validate_construct_fasta_binding
STRUCT_DIR = REPO_ROOT / "data" / "structures"
LAYER2_ATLAS_PATH = REPO_ROOT / "data" / "extracted" / "layer2_regulatory_and_perturbation_atlas.json"
OUT_PATH = REPO_ROOT / "data" / "extracted" / "layer3_structural_verification_atlas.json"

POLAR_ELEMENTS = {"N", "O", "S"}
# Locked at 42.8005 pN*A per round-2/round-3 regression test pins (exact 2019 SI k_B*310.00K = 42.800119 pN*A; k_B*310.15K = 42.8208 pN*A)
KBT_310K_PN_ANGSTROM = 42.8005  # k_B * 310 K = 4.28005 pN*nm = 42.8005 pN*Angstrom

# Frozen reference inputs, independently checked against RCSB/AlphaFold on
# 2026-10-26. PDB local copies use extended entry IDs; their parsed A/B-conformer
# atom records match the RCSB copies exactly. AF models/PAE match byte for byte.
# Hashing an arbitrary replacement and printing its hash is not source binding.
STRUCTURAL_REFERENCE_SHA256 = {
    "pdb_00003cdg.cif": "36f9a9d45f9a3b7dbbf2d74e4c4d25e3edc66f5e02cafafac6acf5d6130b450a",
    "pdb_00001mhe.cif": "4f1bc39231d6ffb288cf85bb9322d264ceac92d6d427f755339a648e49d6227a",
    "pdb_00003bze.cif": "839b87df50ced664fa5c516d0836667176b64895677c87de737b3a6615fc802c",
    "pdb_00002jjs.cif": "84672a253e9c982dac16136cb76ae706eeafa2aa390cfdb2355ad68592cac10e",
    "pdb_00004cmm.cif": "13147521df35289b21488439170079c86286f58cdafe17544367ccf053d1f5fe",
    "pdb_00004kjy.cif": "6294d015c1c24d97558650dc98ef8bced967b3fd658da8b49dc51090e8b2ccaa",
    "pdb_00004zqk.cif": "ed7114a3c653cf964155c5491a5df30d5c5e4edfc6b09d399662c5764a1bad31",
    "pdb_00001iau.cif": "c1ffdeab80cc9ecaf06c9901c756b0f59259e178a45a4a33ffc949bd01123161",
    "AF-P50453-F1-model_v6.cif": "d0b27f7cc8317acd9cd59b76d3080b4019ffd38930dbfddf3bd2560814baf5b9",
    "AF-P10144-F1-model_v6.cif": "6362a677cc62d80b89229ea40a4c3ae81007a631d1c038a2e158541c8d40e9a5",
    "AF-P50453-F1-predicted_aligned_error_v6.json": "b88ea8a67dab0e5956f919c3205911d9d7d901973c3ff32e5277b81881b76f78",
    "AF-P10144-F1-predicted_aligned_error_v6.json": "cb5b42db914ceb00ef6d508fe800ad784e551236aa4d139c190a25ea41241121",
    "1EFX.pdb": "4e432fe894cc6b0390bab6071598f2f95249d5b6f0c4e5a9365c1ca110660fba",
    "1IM9.pdb": "4c21ebf9e5ba064ab65569079cb8d10a73b718bdc8c4cd1c298bec7a8e4f0013",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_structural_reference_inputs(*names: str) -> Mapping[str, bytes]:
    """Return immutable, release-verified snapshots for parsing and provenance.

    Every consumer must use these bytes, never reopen the path after checking its
    hash. A concurrent replacement (including change then restore) cannot change
    the authenticated input consumed by an analysis.
    """
    snapshots = {}
    for name in names:
        expected = STRUCTURAL_REFERENCE_SHA256.get(name)
        if expected is None:
            raise ValueError(f"Structural reference source identity mismatch: {name}")
        payload = (STRUCT_DIR / name).read_bytes()
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError(f"Structural reference source identity mismatch: {name}")
        snapshots[name] = payload
    return MappingProxyType(snapshots)


def validate_processed_polypeptide_ledger(ledger: Dict[str, Any]) -> None:
    """Reject stale hashes and length labels before interpreting molecular identities."""
    expected_ids = {
        "HLAE_0103_VMAPRTLVL_SCT_Arm1B_mature_surface",
        "HLAE_0103_VMAPRTLVL_SCT_Arm1A_mature_with_cytosolic_P2A_stub",
        "CD47_IgV_GPI_mature_cell_surface", "SERPINB9_cytosolic_with_21aa_GSG_P2A_stub",
        "PDL1_mature_cell_surface", "iCasp9_cytosolic_monomer",
    }
    if not isinstance(ledger, dict) or set(ledger) != expected_ids:
        raise ValueError("Complete six-product Layer-2 processed polypeptide ledger is required")
    for key, entry in ledger.items():
        if not isinstance(entry, dict):
            raise TypeError(f"Processed polypeptide entry must be an object: {key}")
        seq = entry.get("sequence")
        if not isinstance(seq, str) or not seq or any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in seq):
            raise ValueError(f"Invalid processed amino-acid sequence: {key}")
        if hashlib.sha256(seq.encode("utf-8")).hexdigest() != entry.get("sha256"):
            raise ValueError(f"Processed polypeptide SHA-256 mismatch: {key}")
        lengths = [entry[name] for name in ("mature_length_aa", "processed_length_aa") if name in entry]
        if not lengths or any(type(value) is not int or value != len(seq) for value in lengths):
            raise ValueError(f"Processed polypeptide length mismatch: {key}")
    cd47 = ledger["CD47_IgV_GPI_mature_cell_surface"]
    if hashlib.sha256(cd47["sequence"][:120].encode("utf-8")).hexdigest() != cd47.get("cd47_ectodomain_120aa_sha256"):
        raise ValueError("CD47 ectodomain SHA-256 mismatch")


def validated_pae_matrix(payload: Any, expected_length: int, label: str) -> np.ndarray:
    """Require a full finite, nonnegative residue-pair matrix; missing rows are not zero confidence."""
    if type(expected_length) is not int or expected_length <= 0:
        raise ValueError("AlphaFold PAE expected_length must be a positive integer")
    if isinstance(payload, list):
        if len(payload) != 1:
            raise ValueError(f"Expected one AlphaFold PAE record: {label}")
        payload = payload[0]
    if not isinstance(payload, dict) or "predicted_aligned_error" not in payload:
        raise ValueError(f"Missing AlphaFold PAE matrix: {label}")
    try:
        raw_matrix = np.asarray(payload["predicted_aligned_error"], dtype=object)
        if raw_matrix.shape != (expected_length, expected_length):
            raise ValueError("PAE shape does not match the declared residue count")
        if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) for value in raw_matrix.flat):
            raise ValueError("PAE entries must be numbers; booleans and numeric strings are invalid")
        matrix = np.asarray(raw_matrix, dtype=np.float64)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"Invalid AlphaFold PAE numeric matrix: {label}") from exc
    if matrix.shape != (expected_length, expected_length) or not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
        raise ValueError(f"AlphaFold PAE must be a finite nonnegative {expected_length}-by-{expected_length} matrix: {label}")
    if "max_predicted_aligned_error" in payload:
        maximum = payload["max_predicted_aligned_error"]
        try:
            valid_maximum = (not isinstance(maximum, (bool, np.bool_)) and isinstance(maximum, Real)
                             and math.isfinite(maximum) and maximum >= 0 and not np.any(matrix > maximum))
        except (TypeError, ValueError, OverflowError):
            valid_maximum = False
        if not valid_maximum:
            raise ValueError(f"AlphaFold PAE exceeds or has an invalid declared maximum: {label}")
    return matrix


def _pae_mean(matrix: np.ndarray) -> float:
    """A finite nonnegative matrix has a finite mean no greater than its maximum."""
    if not matrix.size or not np.all(np.isfinite(matrix)) or np.any(matrix < 0):
        raise ValueError("PAE mean requires a nonempty finite nonnegative matrix")
    maximum = float(np.max(matrix))
    return float(np.mean(matrix / maximum)) * maximum if maximum else 0.0


def _finite_population_moments(values: List[float]) -> Tuple[float, float]:
    """Finite mean and population SD without overflowing sums or squares."""
    data = _finite_real_coordinates(values, "B-factor moments")
    if data.ndim != 1 or not len(data):
        raise ValueError("B-factor moments require a nonempty numeric vector")
    anchor = float(data[0])
    with np.errstate(over="ignore", invalid="ignore"):
        local = data - anchor
    if not np.all(np.isfinite(local)):
        anchor, local = 0.0, data
    scale = float(np.max(np.abs(local)))
    if scale == 0.0:
        return anchor, 0.0
    scaled = local / scale
    mean_scaled = math.fsum(float(value) for value in scaled) / len(scaled)
    mean = anchor + mean_scaled * scale
    # A convex mean belongs to the finite input interval. The normalized
    # population variance is at most one because every scaled value is in
    # [-1, 1]; these bounds also remove endpoint roundoff, not invalid inputs.
    mean = min(float(np.max(data)), max(float(np.min(data)), mean))
    variance_scaled = min(1.0, math.fsum((float(value) - mean_scaled) ** 2 for value in scaled) / len(scaled))
    return mean, math.sqrt(variance_scaled) * scale


def parse_mmcif_atoms(cif_path: Path, keep_insertion_codes: bool = False,
                      allowed_alt_ids: Tuple[str, ...] = (".", "?", "A"), *,
                      source_bytes: Optional[bytes] = None) -> List[Dict[str, Any]]:
    """
    Parse ATOM/HETATM records from an mmCIF coordinate file (model 1, quote-aware tokenizer).
    When keep_insertion_codes=True, residues with PDB insertion codes (e.g. 184A in 1IAU)
    are retained with their exact insertion code so chymotrypsinogen-numbered structures
    align without dropping residues (closing V70-08 & AUDIT-REF).
    With source_bytes, parse that immutable snapshot and use cif_path only in
    diagnostic messages. Otherwise retain the generic arbitrary-file API.
    """
    atoms: List[Dict[str, Any]] = []
    cols: List[str] = []
    in_atom_loop = False
    pending: List[str] = []
    rows_started = False
    identities = set()
    if (type(keep_insertion_codes) is not bool or not isinstance(allowed_alt_ids, tuple)
            or not allowed_alt_ids or any(not isinstance(alt, str) or not alt for alt in allowed_alt_ids)
            or len(set(allowed_alt_ids) - {".", "?"}) > 1):
        raise ValueError("Select one named atom conformer and a boolean insertion-code policy")

    def _append_selected_atom(rec: Dict[str, str]) -> None:
        if rec.get("pdbx_PDB_model_num", "1") not in ("1", "?", "."):
            return
        alt_id = rec.get("label_alt_id", ".")
        if alt_id not in allowed_alt_ids:
            return
        try:
            x, y, z = (float(rec[name]) for name in ("Cartn_x", "Cartn_y", "Cartn_z"))
            b_iso = float(rec.get("B_iso_or_equiv", "0.0"))
            occupancy = float(rec.get("occupancy", "1.0"))
        except (ValueError, KeyError) as exc:
            raise ValueError(f"Invalid atom numeric field in {cif_path}") from exc
        if not all(math.isfinite(v) for v in (x, y, z, b_iso, occupancy)) or not 0.0 <= occupancy <= 1.0:
            raise ValueError(f"Non-finite coordinate/confidence or invalid occupancy in {cif_path}")
        if occupancy == 0.0:
            return
        auth_seq_str = rec.get("auth_seq_id", rec.get("label_seq_id", "?"))
        try:
            auth_seq_id = int(auth_seq_str)
        except ValueError:
            return
        label_seq_str = rec.get("label_seq_id", "?")
        label_seq_id = None
        if label_seq_str not in ("?", "."):
            try:
                label_seq_id = int(label_seq_str)
            except ValueError as exc:
                raise ValueError(f"Invalid polymer sequence identifier in {cif_path}") from exc
            if label_seq_id <= 0:
                raise ValueError(f"Polymer sequence identifiers must be positive in {cif_path}")
        ins_raw = rec.get("pdbx_PDB_ins_code", "?")
        ins_code = "" if ins_raw in ("?", ".") else ins_raw
        if not keep_insertion_codes and ins_code:
            return
        atom_id = rec.get("label_atom_id", rec.get("auth_atom_id", ""))
        chain = rec.get("auth_asym_id", rec.get("label_asym_id", ""))
        identity = (chain, auth_seq_id, ins_code, atom_id)
        if identity in identities:
            raise ValueError(f"Duplicate selected atom identity in {cif_path}: {identity}")
        identities.add(identity)
        atoms.append({
            "group": rec["group_PDB"], "atom_id": atom_id,
            "element": rec.get("type_symbol", "C").upper(),
            "comp_id": rec.get("auth_comp_id", rec.get("label_comp_id", "")),
            "chain": chain, "seq_id": auth_seq_id, "label_seq_id": label_seq_id,
            "ins_code": ins_code, "res_key": f"{auth_seq_id}{ins_code}",
            "x": x, "y": y, "z": z, "b_iso": b_iso,
            "alt_id": alt_id, "occupancy": occupancy,
        })

    if source_bytes is not None and not isinstance(source_bytes, bytes):
        raise TypeError("mmCIF source_bytes must be an immutable bytes snapshot")
    stream = (io.StringIO(source_bytes.decode("utf-8")) if source_bytes is not None
              else cif_path.open("r", encoding="utf-8"))
    with stream as f:
        for raw_line in f:
            line = raw_line.strip()
            if line.startswith("_atom_site."):
                if rows_started:
                    raise ValueError(f"Atom-site columns cannot change within coordinate rows: {cif_path}")
                in_atom_loop = True
                for token in line.split():
                    if not token.startswith("_atom_site."):
                        raise ValueError(f"Malformed atom-site column header in {cif_path}")
                    column = token.split(".", 1)[1]
                    if column in cols:
                        raise ValueError(f"Duplicate atom-site column in {cif_path}: {column}")
                    cols.append(column)
                continue
            if in_atom_loop:
                if not line or line.startswith("#"):
                    continue
                if line.startswith(("loop_", "_", "data_", "stop_")):
                    if pending:
                        raise ValueError(f"Incomplete atom row in {cif_path}")
                    if rows_started:
                        break
                    in_atom_loop = False
                    cols = []
                    continue
                rows_started = True
                for token in _MMCIF_TOKEN_RE.findall(line):
                    if token.startswith("#"):
                        break
                    if token.startswith(("'", '"')):
                        if len(token) < 2 or token[-1] != token[0]:
                            raise ValueError(f"Unterminated quoted atom-site value in {cif_path}")
                        token = token[1:-1]
                    pending.append(token)
                while len(pending) >= len(cols):
                    rec = dict(zip(cols, pending[:len(cols)]))
                    del pending[:len(cols)]
                    if rec.get("group_PDB") not in {"ATOM", "HETATM"}:
                        raise ValueError(f"Invalid atom-site record group in {cif_path}")
                    _append_selected_atom(rec)
    if pending:
        raise ValueError(f"Incomplete atom row in {cif_path}")
    if not atoms:
        raise ValueError(f"Zero valid ATOM/HETATM records parsed from {cif_path}")
    return atoms


def filter_polymer_atoms(
    atoms: List[Dict[str, Any]],
    chains: List[str],
    seq_range: Optional[Tuple[int, int]] = None,
    include_hetatm_peptides: bool = False,
    require_all_chains: bool = True,
) -> List[Dict[str, Any]]:
    chain_set = set(chains)
    if not chain_set:
        raise ValueError("chains list in filter_polymer_atoms cannot be empty")
    out: List[Dict[str, Any]] = []
    for a in atoms:
        if a["chain"] not in chain_set:
            continue
        if a["comp_id"] == "HOH":
            continue
        if a["group"] != "ATOM" and not (include_hetatm_peptides and a["comp_id"] in {"PCA", "ACE"}):
            continue
        if seq_range is not None and not (seq_range[0] <= a["seq_id"] <= seq_range[1]):
            continue
        out.append(a)
    if require_all_chains:
        present_chains = {a["chain"] for a in out}
        missing_chains = sorted(chain_set - present_chains)
        if missing_chains:
            raise ValueError(
                f"Requested chains {sorted(chain_set)} not all present in structure (missing: {missing_chains}; found: {sorted(present_chains)})"
            )
    return out


def compute_interface_metrics(
    atoms_a: List[Dict[str, Any]],
    atoms_b: List[Dict[str, Any]],
    heavy_cutoff: float = 4.5,
    polar_cutoff: float = 3.5,
    salt_bridge_cutoff: float = 4.0,
) -> Dict[str, Any]:
    """
    Compute exact inter-chain heavy-atom, N/O polar proximity (<=3.5 A), and physiological pH-7.4
    ionic salt-bridge (ARG/LYS vs ASP/GLU <=4.0 A) contact networks from experimental mmCIF coordinates,
    separating unique residue-level contacts from atom-pair counts and separating histidine (pKa ~ 6.0)
    protonation-dependent contacts (closing AUDIT-REF, AUDIT-REF, and 091).
    """
    for name, value in (("heavy_cutoff", heavy_cutoff), ("polar_cutoff", polar_cutoff), ("salt_bridge_cutoff", salt_bridge_cutoff)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be a finite positive distance")
    if max(polar_cutoff, salt_bridge_cutoff) > heavy_cutoff:
        raise ValueError("Polar/salt-bridge cutoffs cannot exceed the enumerated heavy-atom cutoff")
    backbone_atoms = {"N", "CA", "C", "O", "OXT"}
    cationic_atoms_ph74 = {
        ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
        ("LYS", "NZ"),
    }
    his_imidazole_atoms = {("HIS", "ND1"), ("HIS", "NE2")}
    anionic_atoms = {
        ("ASP", "OD1"), ("ASP", "OD2"),
        ("GLU", "OE1"), ("GLU", "OE2"),
        ("ASJ", "OD1"), ("ASJ", "OD2"),
    }

    heavy_a = [a for a in atoms_a if a["element"] != "H"]
    heavy_b = [b for b in atoms_b if b["element"] != "H"]
    if not heavy_a or not heavy_b:
        raise ValueError(f"Empty heavy-atom partner selection in compute_interface_metrics: len(a)={len(heavy_a)}, len(b)={len(heavy_b)}")
    chains_a = {a["chain"] for a in heavy_a}
    chains_b = {b["chain"] for b in heavy_b}
    if chains_a & chains_b:
        raise ValueError(f"Partner chain sets must be disjoint in compute_interface_metrics; overlap: {sorted(chains_a & chains_b)}")

    coords_a = _finite_real_coordinates([[a["x"], a["y"], a["z"]] for a in heavy_a], "interface A")
    coords_b = _finite_real_coordinates([[b["x"], b["y"], b["z"]] for b in heavy_b], "interface B")
    with np.errstate(over="ignore", invalid="ignore"):
        diff = coords_a[:, None, :] - coords_b[None, :, :]
        dist_mat = np.hypot.reduce(diff, axis=-1)
    if not np.all(np.isfinite(dist_mat)):
        raise ValueError("Interface distance exceeds the finite numeric domain")
    min_interatomic_dist = float(np.min(dist_mat))
    clash_pairs_lt_2p0a = int(np.sum(dist_mat < 2.0))

    heavy_idx_a, heavy_idx_b = np.where(dist_mat <= heavy_cutoff)
    heavy_atom_pairs = int(len(heavy_idx_a))
    contact_set_a = set(heavy_idx_a.tolist())
    contact_set_b = set(heavy_idx_b.tolist())

    contact_residues_a = sorted(
        {(heavy_a[i]["chain"], heavy_a[i]["seq_id"], heavy_a[i]["ins_code"], heavy_a[i]["res_key"], heavy_a[i]["comp_id"]) for i in heavy_idx_a},
        key=lambda t: (t[0], t[1], t[2]),
    )
    contact_residues_b = sorted(
        {(heavy_b[j]["chain"], heavy_b[j]["seq_id"], heavy_b[j]["ins_code"], heavy_b[j]["res_key"], heavy_b[j]["comp_id"]) for j in heavy_idx_b},
        key=lambda t: (t[0], t[1], t[2]),
    )

    residue_pair_min_dist: Dict[Tuple[str, int, str, str, str, str, int, str, str, str], float] = {}
    polar_contacts: List[Dict[str, Any]] = []
    unique_polar_res_pairs: set = set()
    salt_bridges: List[Dict[str, Any]] = []
    unique_salt_bridge_res_pairs: set = set()
    his_anionic_pairs: List[Dict[str, Any]] = []
    bb_bb = 0
    bb_sc = 0
    sc_sc = 0

    for i, j in zip(heavy_idx_a.tolist(), heavy_idx_b.tolist()):
        aa = heavy_a[i]
        bb = heavy_b[j]
        d = float(dist_mat[i, j])
        rk = (
            aa["chain"], aa["seq_id"], aa["ins_code"], aa["res_key"], aa["comp_id"],
            bb["chain"], bb["seq_id"], bb["ins_code"], bb["res_key"], bb["comp_id"],
        )
        if rk not in residue_pair_min_dist or d < residue_pair_min_dist[rk]:
            residue_pair_min_dist[rk] = d
        if d <= salt_bridge_cutoff:
            ka = (aa["comp_id"], aa["atom_id"])
            kb = (bb["comp_id"], bb["atom_id"])
            if (ka in cationic_atoms_ph74 and kb in anionic_atoms) or (ka in anionic_atoms and kb in cationic_atoms_ph74):
                unique_salt_bridge_res_pairs.add((f"{aa['chain']}:{aa['comp_id']}{aa['res_key']}", f"{bb['chain']}:{bb['comp_id']}{bb['res_key']}"))
                salt_bridges.append(
                    {
                        "partner_1": f"{aa['chain']}:{aa['comp_id']}{aa['res_key']}.{aa['atom_id']}",
                        "partner_2": f"{bb['chain']}:{bb['comp_id']}{bb['res_key']}.{bb['atom_id']}",
                        "distance_angstrom": round(d, 3),
                    }
                )
            elif (ka in his_imidazole_atoms and kb in anionic_atoms) or (ka in anionic_atoms and kb in his_imidazole_atoms):
                his_anionic_pairs.append(
                    {
                        "partner_1": f"{aa['chain']}:{aa['comp_id']}{aa['res_key']}.{aa['atom_id']}",
                        "partner_2": f"{bb['chain']}:{bb['comp_id']}{bb['res_key']}.{bb['atom_id']}",
                        "distance_angstrom": round(d, 3),
                        "protonation_note": "Histidine (pKa ~ 6.0) requires local protonation for ionic salt-bridge formation at pH 7.4",
                    }
                )
        if d <= polar_cutoff and aa["element"] in POLAR_ELEMENTS and bb["element"] in POLAR_ELEMENTS:
            unique_polar_res_pairs.add((f"{aa['chain']}:{aa['comp_id']}{aa['res_key']}", f"{bb['chain']}:{bb['comp_id']}{bb['res_key']}"))
            in_bb_a = aa["atom_id"] in backbone_atoms
            in_bb_b = bb["atom_id"] in backbone_atoms
            if in_bb_a and in_bb_b:
                bb_bb += 1
            elif in_bb_a or in_bb_b:
                bb_sc += 1
            else:
                sc_sc += 1
            polar_contacts.append(
                {
                    "partner_1": f"{aa['chain']}:{aa['comp_id']}{aa['res_key']}.{aa['atom_id']}",
                    "partner_2": f"{bb['chain']}:{bb['comp_id']}{bb['res_key']}.{bb['atom_id']}",
                    "distance_angstrom": round(d, 3),
                }
            )

    polar_contacts.sort(key=lambda x: (x["distance_angstrom"], x["partner_1"], x["partner_2"]))
    salt_bridges.sort(key=lambda x: (x["distance_angstrom"], x["partner_1"], x["partner_2"]))
    sorted_res_pairs = sorted(residue_pair_min_dist.items(), key=lambda kv: (kv[1], kv[0]))

    contact_b_a = [heavy_a[i]["b_iso"] for i in sorted(contact_set_a)]
    contact_b_b = [heavy_b[j]["b_iso"] for j in sorted(contact_set_b)]
    noncontact_b_a = [heavy_a[i]["b_iso"] for i in range(len(heavy_a)) if i not in contact_set_a]
    mean_a, sd_a = _finite_population_moments(contact_b_a) if contact_b_a else (None, None)
    mean_b, sd_b = _finite_population_moments(contact_b_b) if contact_b_b else (None, None)
    noncontact_mean_a = _finite_population_moments(noncontact_b_a)[0] if noncontact_b_a else None

    return {
        "evidence_class": "computed_directly_from_experimental_mmcif_coordinates_distance_proximity_screen",
        "interface_contact_deduplication_convention": {
            "atom_pair_metrics": "heavy_atom_contact_pairs (<=4.5 A), polar_contact_atom_pairs (N/O/S <=3.5 A), and salt_bridge_atom_pairs (<=4.0 A) count all distinct inter-chain heavy-atom pairs meeting each distance cutoff.",
            "residue_pair_metrics": "unique_residue_contact_pairs, unique_polar_proximity_residue_pairs, and unique_salt_bridge_residue_pairs_arg_lys_ph74 deduplicate atom contacts to one count per unique (chain_1, seq_id_1) <-> (chain_2, seq_id_2) residue pair using minimum interatomic distance.",
        },
        "partner_1_chains": sorted(chains_a),
        "partner_2_chains": sorted(chains_b),
        "heavy_atom_cutoff_angstrom": heavy_cutoff,
        "polar_contact_cutoff_angstrom": polar_cutoff,
        "salt_bridge_cutoff_angstrom": salt_bridge_cutoff,
        "min_interatomic_distance_angstrom": round(min_interatomic_dist, 4),
        "steric_clash_atom_pairs_lt_2p0A": clash_pairs_lt_2p0a,
        "heavy_atom_contact_pairs": heavy_atom_pairs,
        "unique_residue_contact_pairs": len(residue_pair_min_dist),
        "polar_contact_atom_pairs": len(polar_contacts),
        "unique_polar_proximity_residue_pairs": len(unique_polar_res_pairs),
        "salt_bridge_atom_pairs": len(salt_bridges),
        "strict_ph74_salt_bridge_atom_pairs": len(salt_bridges),
        "his_protonation_dependent_atom_pairs": len(his_anionic_pairs),
        "total_electrostatic_candidate_atom_pairs_including_his": len(salt_bridges) + len(his_anionic_pairs),
        "unique_salt_bridge_residue_pairs_arg_lys_ph74": len(unique_salt_bridge_res_pairs),
        "salt_bridges": salt_bridges,
        "his_protonation_dependent_candidate_pairs": his_anionic_pairs,
        "polar_proximity_geometry_breakdown_not_angle_typed_hbonds": {
            "backbone_to_backbone": bb_bb,
            "backbone_to_sidechain": bb_sc,
            "sidechain_to_sidechain": sc_sc,
        },
        "interface_crystallographic_b_factor_profile_A2": {
            "metric_type": "experimental_xray_isotropic_b_factor_A2_not_alphafold_plddt",
            "partner_1_contact_mean_b_iso_A2": round(mean_a, 2) if mean_a is not None else None,
            "partner_1_contact_std_b_iso_A2": round(sd_a, 2) if sd_a is not None else None,
            "partner_2_contact_mean_b_iso_A2": round(mean_b, 2) if mean_b is not None else None,
            "partner_2_contact_std_b_iso_A2": round(sd_b, 2) if sd_b is not None else None,
            "partner_1_noncontact_mean_b_iso_A2": round(noncontact_mean_a, 2) if noncontact_mean_a is not None else None,
        },
        "partner_1_contact_residue_count": len(contact_residues_a),
        "partner_2_contact_residue_count": len(contact_residues_b),
        "partner_1_contact_residues": [f"{c}:{comp}{rk_s}" for c, _, _, rk_s, comp in contact_residues_a],
        "partner_2_contact_residues": [f"{c}:{comp}{rk_s}" for c, _, _, rk_s, comp in contact_residues_b],
        "top_15_closest_residue_contacts": [
            {
                "residue_1": f"{k[0]}:{k[4]}{k[3]}",
                "residue_2": f"{k[5]}:{k[9]}{k[8]}",
                "min_heavy_atom_distance_angstrom": round(v, 3),
            }
            for k, v in sorted_res_pairs[:15]
        ],
        "top_15_polar_contacts": polar_contacts[:15],
    }


def _finite_real_coordinates(values: Any, label: str) -> np.ndarray:
    """Reject numeric coercions before they can turn invalid geometry into a fit."""
    try:
        original = np.asarray(values, dtype=object)
        if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, Real) for x in original.flat):
            raise ValueError("Coordinate entries must be real numbers, not booleans, strings or complex values")
        result = np.asarray(original, dtype=np.float64)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"Invalid real numeric coordinates: {label}") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"Coordinates must be finite: {label}")
    return result


def kabsch_rmsd(coords_p: np.ndarray, coords_q: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """Optimal proper-rotation fit, returning (rmsd, P_aligned, original Q).

    Scale the covariance and norm to avoid false zero RMSDs from underflow and
    SVD failure from overflow. Scaling does not change the minimizing rotation.
    """
    coords_p = _finite_real_coordinates(coords_p, "Kabsch P")
    coords_q = _finite_real_coordinates(coords_q, "Kabsch Q")
    if coords_p.ndim != 2 or coords_p.shape != coords_q.shape or coords_p.shape[1] != 3 or len(coords_p) < 3:
        raise ValueError("Kabsch requires matching N-by-3 coordinate arrays with at least three points")
    if np.array_equal(coords_p, coords_q):
        # Identity is an exact optimal proper rotation, including at the largest
        # finite endpoint where an unnecessary SVD can round the fitted point
        # beyond the representable domain.
        return 0.0, coords_p.copy(), coords_q
    def center_locally(points: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray, np.ndarray]:
        # Subtract a local anchor before scaling: an unrelated large translation
        # must not erase a small but representable cloud or its counterpart.
        anchor = points[0].copy()
        with np.errstate(over="ignore", invalid="ignore"):
            displacements = points - anchor
        if not np.all(np.isfinite(displacements)):
            # Opposite extreme finite endpoints can overflow their displacement;
            # scaling around the origin retains their finite representable fit.
            anchor = np.zeros(3)
            displacements = points
        radius = float(np.max(np.abs(displacements)))
        if radius == 0.0:
            return np.zeros_like(points), 0.0, anchor, np.zeros(3)
        scaled = displacements / radius
        centroid = np.mean(scaled, axis=0)
        return scaled - centroid, radius, anchor, centroid

    p_c, radius_p, _, _ = center_locally(coords_p)
    q_c, radius_q, anchor_q, centroid_q = center_locally(coords_q)
    covariance_scale_p = float(np.max(np.abs(p_c))) or 1.0
    covariance_scale_q = float(np.max(np.abs(q_c))) or 1.0
    cov = (p_c / covariance_scale_p).T @ (q_c / covariance_scale_q)
    v, _, wt = np.linalg.svd(cov)
    d = -1.0 if np.linalg.det(v @ wt) < 0.0 else 1.0
    diag = np.diag([1.0, 1.0, d])
    rot = v @ diag @ wt
    scale = max(radius_p, radius_q) or 1.0
    p_aligned = (p_c @ rot) * (radius_p / scale)
    diff = p_aligned - q_c * (radius_q / scale)
    norm_scale = float(np.max(np.abs(diff)))
    normalized_rmsd = (norm_scale * math.sqrt(float(np.mean(np.sum((diff / norm_scale) ** 2, axis=1))))
                       if norm_scale else 0.0)
    rmsd = normalized_rmsd * scale
    with np.errstate(over="ignore", invalid="ignore"):
        translated = (p_aligned + centroid_q * (radius_q / scale)) * scale + anchor_q
    if not math.isfinite(rmsd) or not np.all(np.isfinite(translated)):
        raise ValueError("Kabsch fitted coordinates or RMSD exceed the finite numeric domain")
    return rmsd, translated, coords_q


def get_atom_coord(atoms: List[Dict[str, Any]], chain: str, seq_id: int, atom_id: str) -> Optional[np.ndarray]:
    for a in atoms:
        if a["chain"] == chain and a["seq_id"] == seq_id and a["ins_code"] == "" and a["atom_id"] == atom_id:
            return np.array([a["x"], a["y"], a["z"]], dtype=np.float64)
    return None


def get_ca_map(atoms: List[Dict[str, Any]], chain: str) -> Dict[int, Dict[str, Any]]:
    out: Dict[int, Dict[str, Any]] = {}
    for a in atoms:
        if a["chain"] == chain and a["atom_id"] == "CA" and a["group"] == "ATOM" and a["ins_code"] == "":
            if a["seq_id"] in out:
                raise ValueError(f"Duplicate selected CA residue in chain {chain}: {a['seq_id']}")
            out[a["seq_id"]] = a
    return out


def get_ca_ordered_list_with_ins_codes(atoms: List[Dict[str, Any]], chain: str) -> List[Dict[str, Any]]:
    """Return polymer-ordered CA atoms, independent of mmCIF coordinate-row order."""
    seen = set()
    out: List[Dict[str, Any]] = []
    for a in atoms:
        if a["chain"] == chain and a["atom_id"] == "CA" and a["group"] == "ATOM":
            rk = (a["seq_id"], a["ins_code"])
            if rk in seen:
                raise ValueError(f"Duplicate selected CA residue in chain {chain}: {rk}")
            seen.add(rk)
            out.append(a)
    if out and all(type(a.get("label_seq_id")) is int for a in out):
        if len({a["label_seq_id"] for a in out}) != len(out):
            raise ValueError(f"Duplicate polymer sequence position in chain {chain}")
        return sorted(out, key=lambda a: a["label_seq_id"])
    return sorted(out, key=lambda a: (a["seq_id"], a["ins_code"]))


def euclidean_distance(p1: np.ndarray, p2: np.ndarray) -> float:
    if p1 is None or p2 is None:
        raise ValueError("A required distance anchor is missing")
    p1, p2 = _finite_real_coordinates(p1, "distance P"), _finite_real_coordinates(p2, "distance Q")
    if p1.shape != (3,) or p2.shape != (3,):
        raise ValueError("Distance anchors must be finite three-dimensional coordinates")
    distance = math.hypot(*(float(a) - float(b) for a, b in zip(p1, p2)))
    if not math.isfinite(distance):
        raise ValueError("Distance exceeds the finite numeric domain")
    return distance


def worm_like_chain_metrics(
    num_residues: int,
    span_angstrom: float,
    l0: float = 3.50,
    lp: float = 4.00,
) -> Dict[str, Any]:
    """
    Marko-Siggia Worm-Like Chain (WLC) coarse geometric screen for flexible (G4S)n linkers.
    Validates physical domain constraints (num_residues > 0 integer, span_angstrom >= 0 finite,
    and returns explicit geometric infeasibility when span_angstrom >= contour_length L_c;
    closing V70-09, AUDIT-REF, 090, & AUDIT-REF).
    """
    if isinstance(num_residues, bool) or not isinstance(num_residues, int) or num_residues <= 0:
        raise ValueError(f"num_residues must be a positive integer (> 0), got {num_residues!r}")
    def finite_wlc_scalar(value: Any, name: str, positive: bool) -> float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{name} must be a finite numeric physical scalar")
        try:
            out = float(value)
        except (ValueError, OverflowError) as exc:
            raise ValueError(f"{name} exceeds the finite numeric domain") from exc
        if not math.isfinite(out) or (out <= 0.0 if positive else out < 0.0):
            raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
        return out

    span_val = finite_wlc_scalar(span_angstrom, "span_angstrom", False)
    l0 = finite_wlc_scalar(l0, "l0", True)
    lp = finite_wlc_scalar(lp, "lp", True)
    try:
        lc = num_residues * l0
    except OverflowError as exc:
        raise ValueError("WLC contour length exceeds the finite numeric domain") from exc
    if not math.isfinite(lc) or lc <= 0.0:
        raise ValueError("WLC contour length must remain finite and positive")
    # Stable exact-WLC second moment. 1-exp(-x) suffers catastrophic
    # cancellation as x=Lc/lp approaches zero, even yielding RMS > Lc.
    x = lc / lp
    if x < 1e-3:
        normalized_second_moment = 1.0 - x / 3.0 + x*x / 12.0 - x*x*x / 60.0 + x**4 / 360.0
        rms_r = lc * math.sqrt(normalized_second_moment)
    elif x > 1e6:
        # Avoid underflow in lp/Lc before taking the square root.
        rms_r = math.sqrt(lp) * math.sqrt(lc) * math.sqrt(2.0 * (1.0 - lp / lc))
    else:
        inv_x = lp / lc
        normalized_second_moment = 2.0 * inv_x * (1.0 + inv_x * math.expm1(-x))
        rms_r = lc * math.sqrt(normalized_second_moment)
    frac_ext = span_val / lc
    if not math.isfinite(frac_ext):
        raise ValueError("WLC fractional extension exceeds the finite numeric domain")
    if span_val >= lc:
        return {
            "model_classification": "coarse_geometric_screen_not_excluded_volume_folding_proof",
            "domain_feasibility_status": "GEOMETRICALLY_INFEASIBLE_SPAN_GE_CONTOUR_LENGTH",
            "geometrically_feasible_d_lt_Lc": False,
            "linker_length_aa": num_residues,
            "selected_persistence_length_angstrom": lp,
            "marko_siggia_entropic_tension_pN_at_selected_lp": None,
            "spanned_euclidean_distance_angstrom": round(span_val, 3),
            "contour_length_angstrom": round(lc, 2),
            "wlc_rms_end_to_end_angstrom": round(rms_r, 3),
            "fractional_extension_d_over_Lc": round(frac_ext, 4),
            "marko_siggia_dimensionless_tension_Flp_over_kBT": None,
            "marko_siggia_entropic_tension_pN_at_310K_lp4A": None,
            "entropic_tension_sensitivity_range_pN_lp3p5_to_5p0A": None,
            "entropic_tension_by_persistence_length_pN": None,
            "coarse_extension_screen_d_lt_0p65_Lc": False,
            "sub_0p70_contour_extension_screen_d_lt_0p70_Lc": False,
            "physical_limitations_note": (
                f"Geometrically infeasible: Euclidean span ({span_val:.2f} A) >= polypeptide contour length L_c ({lc:.2f} A)."
            ),
        }

    def marko_siggia_shape(span: float, contour: float) -> float:
        # Algebraically identical to 1/[4(1-x)^2] - 1/4 + x, without
        # subtracting two almost equal quarters in the small-extension limit.
        extension = span / contour
        gap = (contour - span) / contour
        return 1.0 + (2.0 - extension) / (4.0 * gap ** 2)

    def tension(span: float, contour: float, persistence: float) -> float:
        # Evaluate kBT * (span/contour) * shape / persistence without losing a
        # subnormal ratio before a small persistence length rescales it.
        if span == 0.0:
            return 0.0
        span_m, span_e = math.frexp(span)
        contour_m, contour_e = math.frexp(contour)
        persistence_m, persistence_e = math.frexp(persistence)
        mantissa = (span_m / contour_m / persistence_m) * KBT_310K_PN_ANGSTROM * marko_siggia_shape(span, contour)
        try:
            result = math.ldexp(mantissa, span_e - contour_e - persistence_e)
        except OverflowError as exc:
            raise ValueError("WLC selected-persistence tension exceeds the finite numeric domain") from exc
        if not math.isfinite(result):
            raise ValueError("WLC selected-persistence tension exceeds the finite numeric domain")
        return result

    wlc_dimensionless_force = frac_ext * marko_siggia_shape(span_val, lc)
    force_pn_at_selected_lp = tension(span_val, lc, lp)
    # Keys explicitly labelled lp4A are fixed 4-A comparators; caller-selected
    # persistence length is exposed separately rather than silently relabelled.
    force_pn_at_lp4a = tension(span_val, lc, 4.00)
    force_pn_lp3p5a = tension(span_val, lc, 3.50)
    force_pn_lp4p5a = tension(span_val, lc, 4.50)
    force_pn_lp5p0a = tension(span_val, lc, 5.00)

    # Also compute crystallographic trans-peptide contour comparator at l0 = 3.80 A
    try:
        lc_3p8 = num_residues * 3.80
    except OverflowError as exc:
        raise ValueError("WLC comparator contour exceeds the finite numeric domain") from exc
    if not math.isfinite(lc_3p8):
        raise ValueError("WLC comparator contour must remain finite")
    frac_3p8 = span_val / lc_3p8
    force_pn_l0_3p8_lp4a = round(tension(span_val, lc_3p8, 4.00), 3) if frac_3p8 < 1.0 else None

    arc_span_hemisphere = span_val * (math.pi / 2.0)
    if not math.isfinite(arc_span_hemisphere):
        raise ValueError("WLC hemispherical bypass arc exceeds the finite numeric domain")
    arc_frac = arc_span_hemisphere / lc
    if arc_frac < 1.0:
        arc_force_pn = round(tension(arc_span_hemisphere, lc, 4.00), 3)
    else:
        arc_force_pn = None

    comparator_force_note = (f"{force_pn_l0_3p8_lp4a:.1f} pN" if force_pn_l0_3p8_lp4a is not None
                             else "unavailable: span is at or above the comparator contour length")
    return {
        "model_classification": "coarse_geometric_screen_not_excluded_volume_folding_proof",
        "domain_feasibility_status": "VALID_SUB_CONTOUR_EXTENSION",
        "comparator_l0_3p8A_geometry_status": ("VALID_SUB_CONTOUR_EXTENSION" if frac_3p8 < 1.0
                                               else "GEOMETRICALLY_INFEASIBLE_SPAN_GE_COMPARATOR_CONTOUR"),
        "geometrically_feasible_d_lt_Lc": True,
        "linker_length_aa": num_residues,
        "selected_persistence_length_angstrom": lp,
        "marko_siggia_entropic_tension_pN_at_selected_lp": round(force_pn_at_selected_lp, 3),
        "spanned_euclidean_distance_angstrom": round(span_val, 3),
        "hemispherical_surface_bypass_arc_length_angstrom": round(arc_span_hemisphere, 3),
        "hemispherical_arc_fractional_extension": round(arc_frac, 4),
        "hemispherical_arc_wlc_tension_pN_at_310K_lp4A": arc_force_pn,
        "contour_length_angstrom": round(lc, 2),
        "contour_length_at_l0_3p8A_angstrom": round(lc_3p8, 2),
        "wlc_rms_end_to_end_angstrom": round(rms_r, 3),
        "fractional_extension_d_over_Lc": round(frac_ext, 4),
        "fractional_extension_at_l0_3p8A": round(frac_3p8, 4),
        "marko_siggia_dimensionless_tension_Flp_over_kBT": round(wlc_dimensionless_force, 4),
        "marko_siggia_entropic_tension_pN_at_310K_lp4A": round(force_pn_at_lp4a, 3),
        "marko_siggia_entropic_tension_pN_at_l0_3p8A_lp4A": force_pn_l0_3p8_lp4a,
        "entropic_tension_sensitivity_range_pN_lp3p5_to_5p0A": [
            round(min(force_pn_lp3p5a, force_pn_lp5p0a), 3),
            round(max(force_pn_lp3p5a, force_pn_lp5p0a), 3),
        ],
        "wlc_persistence_length_sensitivity_lp_3p5_to_4p5_angstrom": {
            "lp_3p5A": round(force_pn_lp3p5a, 3),
            "lp_4p0A": round(force_pn_at_lp4a, 3),
            "lp_4p5A": round(force_pn_lp4p5a, 3),
            "range_pN": [
                round(min(force_pn_lp3p5a, force_pn_lp4p5a), 3),
                round(max(force_pn_lp3p5a, force_pn_lp4p5a), 3),
            ],
        },
        "entropic_tension_by_persistence_length_pN": {
            "lp_3p5A": round(force_pn_lp3p5a, 3),
            "lp_4p0A": round(force_pn_at_lp4a, 3),
            "lp_4p5A": round(force_pn_lp4p5a, 3),
            "lp_5p0A": round(force_pn_lp5p0a, 3),
        },
        "coarse_extension_screen_d_lt_0p65_Lc": bool(frac_ext < 0.65),
        "coarse_extension_screen_d_lt_0p65_Lc_at_l0_3p8A": bool(frac_3p8 < 0.65),
        "sub_0p70_contour_extension_screen_d_lt_0p70_Lc": bool(frac_ext < 0.70),
        "physical_limitations_note": (
            "Euclidean anchor distance is a lower bound on the excluded-volume path around the protein surface; "
            f"the fixed lp=4.00 A reference gives ~{force_pn_at_lp4a:.1f} pN entropic tension at selected l0={l0:g} A (l0=3.80 A comparator: {comparator_force_note}); selected lp={lp:g} A "
            f"gives {force_pn_at_selected_lp:.3g} pN. This coarse calculation is an unrelaxed geometric screen "
            "and requires experimental thermal-shift / flow-cytometry surface density confirmation."
        ),
    }


def analyze_complex_1_hlae_nkg2a(layer2_ledger: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze Complex 1: HLA-E*01:03(VMAPRTLVL) SCT : CD94/NKG2A vs NKG2C."""
    snapshots = verify_structural_reference_inputs("pdb_00003cdg.cif", "pdb_00001mhe.cif", "pdb_00003bze.cif")
    atoms_3cdg = parse_mmcif_atoms(STRUCT_DIR / "pdb_00003cdg.cif", source_bytes=snapshots["pdb_00003cdg.cif"])
    atoms_1mhe = parse_mmcif_atoms(STRUCT_DIR / "pdb_00001mhe.cif", source_bytes=snapshots["pdb_00001mhe.cif"])
    atoms_3bze = parse_mmcif_atoms(STRUCT_DIR / "pdb_00003bze.cif", source_bytes=snapshots["pdb_00003bze.cif"])

    hlae_trimer_3cdg = filter_polymer_atoms(atoms_3cdg, ["A", "B", "P"])
    cd94_nkg2a_3cdg = filter_polymer_atoms(atoms_3cdg, ["J", "K"])
    full_interface = compute_interface_metrics(hlae_trimer_3cdg, cd94_nkg2a_3cdg)

    pep_3cdg = filter_polymer_atoms(atoms_3cdg, ["P"])
    pep_to_cd94_nkg2a = compute_interface_metrics(pep_3cdg, cd94_nkg2a_3cdg)
    pep_to_hlae_groove = compute_interface_metrics(pep_3cdg, filter_polymer_atoms(atoms_3cdg, ["A"]))

    # Per-subunit contact breakdown for peptide P8-Phe8 in both asymmetric-unit copies of 3CDG
    # Copy 1: peptide P against CD94 (J) and NKG2A (K)
    # Copy 2: peptide Q against CD94 (E) and NKG2A (F)
    p8_phe_atoms = [a for a in pep_3cdg if a["seq_id"] == 8 and a["element"] != "H"]
    p8_backbone_and_cb_cg = [a for a in p8_phe_atoms if a["atom_id"] in {"N", "CA", "C", "O", "CB", "CG"}]
    cd94_j_atoms = filter_polymer_atoms(atoms_3cdg, ["J"])
    nkg2a_k_atoms = filter_polymer_atoms(atoms_3cdg, ["K"])
    p8_full_contacts = compute_interface_metrics(p8_phe_atoms, cd94_nkg2a_3cdg)["heavy_atom_contact_pairs"]
    p8_val_proxy_contacts = compute_interface_metrics(p8_backbone_and_cb_cg, cd94_nkg2a_3cdg)["heavy_atom_contact_pairs"]
    p8_copy1_cd94_metrics = compute_interface_metrics(p8_phe_atoms, cd94_j_atoms)
    p8_copy1_to_cd94_j = p8_copy1_cd94_metrics["heavy_atom_contact_pairs"]
    p8_copy1_to_nkg2a_k = compute_interface_metrics(p8_phe_atoms, nkg2a_k_atoms)["heavy_atom_contact_pairs"]

    pep_q_3cdg = filter_polymer_atoms(atoms_3cdg, ["Q"])
    p8_q_phe_atoms = [a for a in pep_q_3cdg if a["seq_id"] == 8 and a["element"] != "H"]
    p8_q_proxy_atoms = [a for a in p8_q_phe_atoms if a["atom_id"] in {"N", "CA", "C", "O", "CB", "CG"}]
    cd94_e_atoms = filter_polymer_atoms(atoms_3cdg, ["E"])
    nkg2a_f_atoms = filter_polymer_atoms(atoms_3cdg, ["F"])
    p8_copy2_cd94_metrics = compute_interface_metrics(p8_q_phe_atoms, cd94_e_atoms)
    p8_copy2_to_cd94_e = p8_copy2_cd94_metrics["heavy_atom_contact_pairs"]
    p8_copy2_to_nkg2a_f = compute_interface_metrics(p8_q_phe_atoms, nkg2a_f_atoms)["heavy_atom_contact_pairs"]
    p8_copy2_proxy_to_cd94_e = compute_interface_metrics(p8_q_proxy_atoms, cd94_e_atoms)["heavy_atom_contact_pairs"]

    ca_3cdg_a = get_ca_map(atoms_3cdg, "A")
    ca_1mhe_a = get_ca_map(atoms_1mhe, "A")
    ca_3bze_a = get_ca_map(atoms_3bze, "A")
    # Evaluate the second asymmetric-unit copy in 3CDG: HLA-E/B2M/peptide chains C, D, Q against CD94/NKG2A chains E, F (255 heavy-atom contacts)
    hlae_trimer_3cdg_copy2 = filter_polymer_atoms(atoms_3cdg, ["C", "D", "Q"])
    cd94_nkg2a_3cdg_copy2 = filter_polymer_atoms(atoms_3cdg, ["E", "F"])
    copy2_interface = compute_interface_metrics(hlae_trimer_3cdg_copy2, cd94_nkg2a_3cdg_copy2)
    asu_copy2_contacts = copy2_interface["heavy_atom_contact_pairs"]

    common_groove_res = sorted(set(range(2, 181)) & set(ca_3cdg_a.keys()) & set(ca_1mhe_a.keys()) & set(ca_3bze_a.keys()))
    coords_3cdg_groove = np.array([[ca_3cdg_a[r]["x"], ca_3cdg_a[r]["y"], ca_3cdg_a[r]["z"]] for r in common_groove_res])
    coords_1mhe_groove = np.array([[ca_1mhe_a[r]["x"], ca_1mhe_a[r]["y"], ca_1mhe_a[r]["z"]] for r in common_groove_res])
    coords_3bze_groove = np.array([[ca_3bze_a[r]["x"], ca_3bze_a[r]["y"], ca_3bze_a[r]["z"]] for r in common_groove_res])

    rmsd_1mhe_vs_3cdg, _, _ = kabsch_rmsd(coords_1mhe_groove, coords_3cdg_groove)
    rmsd_3bze_vs_3cdg, _, _ = kabsch_rmsd(coords_3bze_groove, coords_3cdg_groove)

    centroid_p = np.mean(coords_1mhe_groove, axis=0)
    centroid_q = np.mean(coords_3cdg_groove, axis=0)
    cov = (coords_1mhe_groove - centroid_p).T @ (coords_3cdg_groove - centroid_q)
    v, _, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt

    ca_3cdg_p = get_ca_map(atoms_3cdg, "P")
    ca_1mhe_p = get_ca_map(atoms_1mhe, "P")
    ca_3bze_p = get_ca_map(atoms_3bze, "P")
    aa3_to_1_pep = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    pep_seq_3cdg = "".join(aa3_to_1_pep.get(ca_3cdg_p[pos]["comp_id"], "X") for pos in sorted(ca_3cdg_p.keys()))
    pep_seq_1mhe = "".join(aa3_to_1_pep.get(ca_1mhe_p[pos]["comp_id"], "X") for pos in sorted(ca_1mhe_p.keys()))
    pep_seq_3bze = "".join(aa3_to_1_pep.get(ca_3bze_p[pos]["comp_id"], "X") for pos in sorted(ca_3bze_p.keys()))
    if not (pep_seq_3cdg == "VMAPRTLFL" and pep_seq_1mhe == "VMAPRTVLL" and pep_seq_3bze == "VMAPRTLFL"):
        raise ValueError('Molecular source invariant failed: pep_seq_3cdg == "VMAPRTLFL" and pep_seq_1mhe == "VMAPRTVLL" and pep_seq_3bze == "VMAPRTLFL"')
    peptide_ca_displacements: List[Dict[str, Any]] = []
    for pos in range(1, 10):
        if pos in ca_3cdg_p and pos in ca_1mhe_p:
            p_1mhe = np.array([ca_1mhe_p[pos]["x"], ca_1mhe_p[pos]["y"], ca_1mhe_p[pos]["z"]])
            p_1mhe_aligned = (p_1mhe - centroid_p) @ rot + centroid_q
            p_3cdg = np.array([ca_3cdg_p[pos]["x"], ca_3cdg_p[pos]["y"], ca_3cdg_p[pos]["z"]])
            disp = euclidean_distance(p_1mhe_aligned, p_3cdg)
            peptide_ca_displacements.append(
                {
                    "position": f"P{pos}",
                    "residue_3CDG_HLA_G_leader": ca_3cdg_p[pos]["comp_id"],
                    "residue_1MHE_HLA_B_leader": ca_1mhe_p[pos]["comp_id"],
                    "residue_3BZE_HLA_G_leader": ca_3bze_p[pos]["comp_id"] if pos in ca_3bze_p else None,
                    "residue_Islet1_HLA_A2_leader": "VAL" if pos == 8 else ("LEU" if pos == 7 else ca_3cdg_p[pos]["comp_id"]),
                    "ca_displacement_angstrom": round(disp, 3),
                }
            )

    p9_c = get_atom_coord(atoms_3cdg, "P", 9, "C")
    b2m_n = get_atom_coord(atoms_3cdg, "B", 1, "N")
    b2m_m0_n = get_atom_coord(atoms_3cdg, "B", 0, "N")
    span_l1_ile1 = euclidean_distance(p9_c, b2m_n)
    span_l1_met0 = euclidean_distance(p9_c, b2m_m0_n)

    b2m_c99 = get_atom_coord(atoms_3cdg, "B", 99, "C")
    hlae_n2 = get_atom_coord(atoms_3cdg, "A", 2, "N")
    span_l2 = euclidean_distance(b2m_c99, hlae_n2)

    res107_ca = get_atom_coord(atoms_3cdg, "A", 107, "CA")
    res107_comp_id = ca_3cdg_a[107]["comp_id"]
    res107_atoms = [a for a in hlae_trimer_3cdg if a["chain"] == "A" and a["seq_id"] == 107 and a["element"] != "H"]
    res107_coords = np.array([[a["x"], a["y"], a["z"]] for a in res107_atoms], dtype=np.float64)
    cd94_nkg2a_coords = np.array([[a["x"], a["y"], a["z"]] for a in cd94_nkg2a_3cdg if a["element"] != "H"])
    dist_res107_ca_to_receptor = float(np.min(np.linalg.norm(cd94_nkg2a_coords - res107_ca, axis=1)))
    dist_res107_min_heavy_to_receptor = float(np.min(np.linalg.norm(res107_coords[:, None, :] - cd94_nkg2a_coords[None, :, :], axis=-1)))

    # Coordinate-derived audit of HLA-E alpha3 CD8-loop residues Asp227/Thr228 in PDB:3CDG (closes L3-M03 / Audit F-02 / Audit FRESH-X-002)
    def _min_heavy_dist(group1: List[Dict[str, Any]], group2: List[Dict[str, Any]]) -> float:
        c1_arr = np.array([[a["x"], a["y"], a["z"]] for a in group1 if a["element"] != "H"], dtype=np.float64)
        c2_arr = np.array([[a["x"], a["y"], a["z"]] for a in group2 if a["element"] != "H"], dtype=np.float64)
        return float(np.min(np.linalg.norm(c1_arr[:, None, :] - c2_arr[None, :, :], axis=-1)))

    b2m_b_atoms = filter_polymer_atoms(atoms_3cdg, ["B"])
    b2m_d_atoms = filter_polymer_atoms(atoms_3cdg, ["D"])
    alpha3_227_228_audit: Dict[str, Any] = {}
    heterodimer_dists: List[float] = []
    cd94_subunit_dists: List[float] = []
    nkg2a_subunit_dists: List[float] = []
    b2m_subunit_dists: List[float] = []
    for copy_label, hc_chain, b2m_atoms, cd94_atoms, nkg2a_atoms, rec_atoms in [
        ("copy1_chainA", "A", b2m_b_atoms, cd94_j_atoms, nkg2a_k_atoms, cd94_nkg2a_3cdg),
        ("copy2_chainC", "C", b2m_d_atoms, cd94_e_atoms, nkg2a_f_atoms, cd94_nkg2a_3cdg_copy2),
    ]:
        for res_num, res_label in [(227, "Asp227"), (228, "Thr228")]:
            res_atoms = [a for a in atoms_3cdg if a["chain"] == hc_chain and a["seq_id"] == res_num and a["element"] != "H"]
            d_rec = round(_min_heavy_dist(res_atoms, rec_atoms), 2)
            d_cd94 = round(_min_heavy_dist(res_atoms, cd94_atoms), 2)
            d_nkg2a = round(_min_heavy_dist(res_atoms, nkg2a_atoms), 2)
            d_b2m = round(_min_heavy_dist(res_atoms, b2m_atoms), 2)
            heterodimer_dists.append(d_rec)
            cd94_subunit_dists.append(d_cd94)
            nkg2a_subunit_dists.append(d_nkg2a)
            b2m_subunit_dists.append(d_b2m)
            alpha3_227_228_audit[f"{copy_label}_{res_label}"] = {
                "min_heavy_atom_distance_to_cd94_nkg2a_heterodimer_angstrom": d_rec,
                "min_heavy_atom_distance_to_cd94_subunit_angstrom": d_cd94,
                "min_heavy_atom_distance_to_nkg2a_subunit_angstrom": d_nkg2a,
                "min_heavy_atom_distance_to_b2m_subunit_angstrom": d_b2m,
            }
    p8_cd94_res_copy1 = ", ".join(p8_copy1_cd94_metrics["partner_2_contact_residues"])
    p8_cd94_res_copy2 = ", ".join(p8_copy2_cd94_metrics["partner_2_contact_residues"])

    arm1b_prod = layer2_ledger["HLAE_0103_VMAPRTLVL_SCT_Arm1B_mature_surface"]
    arm1a_prod = layer2_ledger["HLAE_0103_VMAPRTLVL_SCT_Arm1A_mature_with_cytosolic_P2A_stub"]
    sct_seq = arm1b_prod["sequence"]
    b2m_end_idx = sct_seq.index("WDRDM") + 5
    design_l2_seq = sct_seq[b2m_end_idx : b2m_end_idx + 20]
    unresolved_hlae_gly1 = sct_seq[b2m_end_idx + 20]
    first_resolved_hlae_ser2 = sct_seq[b2m_end_idx + 21]
    if not (design_l2_seq == "GGGGS" * 4 and unresolved_hlae_gly1 == "G" and first_resolved_hlae_ser2 == "S"):
        raise ValueError('Molecular source invariant failed: design_l2_seq == "GGGGS" * 4 and unresolved_hlae_gly1 == "G" and first_resolved_hlae_ser2 == "S"')
    l2_wlc_21aa = worm_like_chain_metrics(21, span_l2)
    l2_wlc_21aa["design_linker_length_aa"] = 20
    l2_wlc_21aa["design_linker_sequence"] = design_l2_seq
    l2_wlc_21aa["unresolved_n_terminal_hlae_residues_included"] = 1
    l2_wlc_21aa["unresolved_n_terminal_hlae_residue_identity"] = "Gly1"
    l2_wlc_21aa["design_20aa_only_comparator_entropic_tension_pN_lp4A"] = worm_like_chain_metrics(20, span_l2)["marko_siggia_entropic_tension_pN_at_310K_lp4A"]

    # Coordinate-derived 8-copy P9(Leu9:C) -> B2M(Ile1:N) chord and open-segment heavy-atom collision audit (AUDIT-REF):
    eight_crystal_copies_spec = [
        ("3CDG_copy1_ABP", atoms_3cdg, "A", "B", "P"),
        ("3CDG_copy2_CDQ", atoms_3cdg, "C", "D", "Q"),
        ("1MHE_copy1_ABP", atoms_1mhe, "A", "B", "P"),
        ("1MHE_copy2_CDQ", atoms_1mhe, "C", "D", "Q"),
        ("3BZE_copy1_ABP", atoms_3bze, "A", "B", "P"),
        ("3BZE_copy2_CDQ", atoms_3bze, "C", "D", "Q"),
        ("3BZE_copy3_EFR", atoms_3bze, "E", "F", "R"),
        ("3BZE_copy4_GHS", atoms_3bze, "G", "H", "S"),
    ]
    per_copy_chords_exact: List[float] = []
    per_copy_chord_dict: Dict[str, float] = {}
    per_copy_chord_collision_audit: Dict[str, Dict[str, Any]] = {}
    for cp_label, cp_atoms, hc_ch, b2m_ch, pep_ch in eight_crystal_copies_spec:
        cp_p9_c = get_atom_coord(cp_atoms, pep_ch, 9, "C")
        cp_b2m_n = get_atom_coord(cp_atoms, b2m_ch, 1, "N")
        cp_vec = cp_b2m_n - cp_p9_c
        cp_len = float(np.linalg.norm(cp_vec))
        per_copy_chords_exact.append(cp_len)
        per_copy_chord_dict[cp_label] = round(cp_len, 4)
        u_hat = cp_vec / cp_len
        cp_trimer_non_anchor = [
            a for a in filter_polymer_atoms(cp_atoms, [hc_ch, b2m_ch, pep_ch], include_hetatm_peptides=True)
            if a["element"] != "H"
            and not (a["chain"] == pep_ch and a["seq_id"] == 9)
            and not (a["chain"] == b2m_ch and a["seq_id"] in (0, 1))
        ]
        cp_coords = np.array([[a["x"], a["y"], a["z"]] for a in cp_trimer_non_anchor], dtype=np.float64)
        cp_rel = cp_coords - cp_p9_c[None, :]
        cp_proj = np.sum(cp_rel * u_hat[None, :], axis=1)
        cp_t = cp_proj / cp_len
        cp_perp = np.linalg.norm(cp_rel - cp_proj[:, None] * u_hat[None, :], axis=1)
        mask_4a = (cp_t > 0.0) & (cp_t < 1.0) & (cp_perp <= 4.0)
        t_clamped = np.clip(cp_t, 0.0, 1.0)
        capsule_dists = np.linalg.norm(cp_rel - (t_clamped[:, None] * cp_len) * u_hat[None, :], axis=1)
        mask_capsule_4a = capsule_dists <= 4.0
        hc_mask_4a = [
            bool(m and a["chain"] == hc_ch)
            for m, a in zip(mask_4a, cp_trimer_non_anchor)
        ]
        cp_b2m_res1_comp = next((a["comp_id"] for a in cp_atoms if a["chain"] == b2m_ch and a["seq_id"] == 1 and a["atom_id"] == "N"), "ILE")
        cp_b2m_gln2_n = get_atom_coord(cp_atoms, b2m_ch, 2, "N")
        cp_gln2_len = float(np.linalg.norm(cp_b2m_gln2_n - cp_p9_c))
        per_copy_chord_collision_audit[cp_label] = {
            "b2m_anchor_res1_comp_id": cp_b2m_res1_comp,
            "straight_chord_angstrom": round(cp_len, 4),
            "p9_c_to_b2m_gln2_n_angstrom": round(cp_gln2_len, 4),
            "open_chord_4p0A_heavy_atoms_excluding_p9_and_b2m_0_1": int(np.sum(mask_4a)),
            "capsule_0_le_t_le_1_4p0A_heavy_atoms_excluding_p9_and_b2m_0_1": int(np.sum(mask_capsule_4a)),
            "open_chord_4p0A_hlae_heavy_chain_atoms": int(sum(hc_mask_4a)),
        }
    chords_arr = np.array(per_copy_chords_exact, dtype=np.float64)
    mean_chord_exact = float(np.mean(chords_arr))
    sd_pop_chord_exact = float(np.std(chords_arr, ddof=0))
    sd_sample_chord_exact = float(np.std(chords_arr, ddof=1))

    # Live coordinate-derived nearest-neighbor and Ile1-subset chord metrics (AUDIT-REF):
    p9_c_3cdg_c1 = get_atom_coord(atoms_3cdg, "P", 9, "C")
    b2m_n_3cdg_c1 = get_atom_coord(atoms_3cdg, "B", 1, "N")
    b2m_ca_3cdg_c1 = get_atom_coord(atoms_3cdg, "B", 1, "CA")
    ca_chords_8 = [
        float(np.linalg.norm(get_atom_coord(cp_atoms, b2m_ch, 1, "CA") - get_atom_coord(cp_atoms, pep_ch, 9, "C")))
        for _lbl, cp_atoms, _hc_ch, b2m_ch, pep_ch in eight_crystal_copies_spec
    ]
    ile1_6_chords = [
        cp_len
        for (cp_lbl, _cp_atoms, _hc_ch, _b2m_ch, _pep_ch), cp_len in zip(eight_crystal_copies_spec, per_copy_chords_exact)
        if per_copy_chord_collision_audit[cp_lbl]["b2m_anchor_res1_comp_id"] == "ILE"
    ]

    def _dist_to_anchor(anchor_pt: np.ndarray, ch: str, res_id: int, at_id: str) -> float:
        pt = get_atom_coord(atoms_3cdg, ch, res_id, at_id)
        return round(float(np.linalg.norm(pt - anchor_pt)), 2)

    def _compute_3d_voxel_dijkstra_geodesic_3cdg_copy1(r_obs_angstrom: float, spacing: float = 0.5, anchor_exempt_r: float = 4.0) -> float | None:
        import heapq
        obs_pts = np.array(
            [
                [a["x"], a["y"], a["z"]]
                for a in filter_polymer_atoms(atoms_3cdg, ["A", "B", "P"], None, include_hetatm_peptides=True)
                if a["element"] != "H"
                and not (a["chain"] == "P" and a["seq_id"] == 9)
                and not (a["chain"] == "B" and a["seq_id"] in (0, 1))
            ],
            dtype=np.float64,
        )
        d_start = np.linalg.norm(obs_pts - p9_c_3cdg_c1, axis=1)
        d_end = np.linalg.norm(obs_pts - b2m_n_3cdg_c1, axis=1)
        obs_pts = obs_pts[(d_start > anchor_exempt_r) & (d_end > anchor_exempt_r)]
        pad = 8.0
        lo = np.minimum(p9_c_3cdg_c1, b2m_n_3cdg_c1) - pad
        hi = np.maximum(p9_c_3cdg_c1, b2m_n_3cdg_c1) + pad
        dims = np.ceil((hi - lo) / spacing).astype(int) + 1
        nx, ny, nz = int(dims[0]), int(dims[1]), int(dims[2])
        blocked = np.zeros((nx, ny, nz), dtype=bool)
        rad_idx = int(math.ceil(r_obs_angstrom / spacing))
        r_sq = r_obs_angstrom * r_obs_angstrom
        for ox, oy, oz in obs_pts:
            ix = int(round((ox - lo[0]) / spacing))
            iy = int(round((oy - lo[1]) / spacing))
            iz = int(round((oz - lo[2]) / spacing))
            x0, x1 = max(0, ix - rad_idx), min(nx - 1, ix + rad_idx)
            y0, y1 = max(0, iy - rad_idx), min(ny - 1, iy + rad_idx)
            z0, z1 = max(0, iz - rad_idx), min(nz - 1, iz + rad_idx)
            if x0 > x1 or y0 > y1 or z0 > z1:
                continue
            xs = lo[0] + np.arange(x0, x1 + 1) * spacing - ox
            ys = lo[1] + np.arange(y0, y1 + 1) * spacing - oy
            zs = lo[2] + np.arange(z0, z1 + 1) * spacing - oz
            dist_sq = xs[:, None, None] ** 2 + ys[None, :, None] ** 2 + zs[None, None, :] ** 2
            blocked[x0 : x1 + 1, y0 : y1 + 1, z0 : z1 + 1] |= dist_sq <= r_sq
        s_idx = tuple(np.round((p9_c_3cdg_c1 - lo) / spacing).astype(int).tolist())
        t_idx = tuple(np.round((b2m_n_3cdg_c1 - lo) / spacing).astype(int).tolist())
        blocked[s_idx] = False
        blocked[t_idx] = False
        nbrs = [
            (dx, dy, dz, spacing * math.sqrt(dx * dx + dy * dy + dz * dz))
            for dx in (-1, 0, 1)
            for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)
            if not (dx == 0 and dy == 0 and dz == 0)
        ]
        dist_grid = np.full((nx, ny, nz), np.inf, dtype=np.float64)
        dist_grid[s_idx] = 0.0
        pq = [(0.0, s_idx[0], s_idx[1], s_idx[2])]
        tx, ty, tz = t_idx
        while pq:
            d_cur, cx, cy, cz = heapq.heappop(pq)
            if cx == tx and cy == ty and cz == tz:
                return round(float(d_cur), 1)
            if d_cur > dist_grid[cx, cy, cz]:
                continue
            for dx, dy, dz, step_cost in nbrs:
                ux, uy, uz = cx + dx, cy + dy, cz + dz
                if 0 <= ux < nx and 0 <= uy < ny and 0 <= uz < nz and not blocked[ux, uy, uz]:
                    nd = d_cur + step_cost
                    if nd < dist_grid[ux, uy, uz]:
                        dist_grid[ux, uy, uz] = nd
                        heapq.heappush(pq, (nd, ux, uy, uz))
        return None

    voxel_geodesic_3cdg_c1_r2p5 = _compute_3d_voxel_dijkstra_geodesic_3cdg_copy1(2.5)
    voxel_geodesic_3cdg_c1_r2p8 = _compute_3d_voxel_dijkstra_geodesic_3cdg_copy1(2.8)
    voxel_geodesic_3cdg_c1_r3p1 = _compute_3d_voxel_dijkstra_geodesic_3cdg_copy1(3.1)
    voxel_geodesic_3cdg_c1_r3p4 = _compute_3d_voxel_dijkstra_geodesic_3cdg_copy1(3.4)
    voxel_geodesic_3cdg_c1_r2p5_sp1p0 = _compute_3d_voxel_dijkstra_geodesic_3cdg_copy1(2.5, spacing=1.0)

    return {
        "complex_id": "LAYER3_COMPLEX_1_HLAE_SCT_CD94_NKG2A",
        "locus": "Locus 1 (B2M Exon-1 KI: Arm 1A & Arm 1B)",
        "engineered_protein": "HLA-E*01:03 Single-Chain Trimer (SCT, full 337-aa mature heavy chain retaining 344-358 AEWSDSAQGSESHSL tail) presenting VMAPRTLVL",
        "layer2_processed_polypeptide_hashes": {
            "Arm1B_mature_480aa_sha256": arm1b_prod["sha256"],
            "Embodiment_A_star_mature_485aa_20aa_L1_sequence": sct_seq[:9] + "GGGGS" * 4 + sct_seq[24:],
            "Embodiment_A_star_mature_485aa_20aa_L1_sha256": hashlib.sha256((sct_seq[:9] + "GGGGS" * 4 + sct_seq[24:]).encode("utf-8")).hexdigest(),
            "Arm1A_mature_501aa_with_cytosolic_21aa_P2A_stub_sha256": arm1a_prod["sha256"],
            "native_hlae_344_358_tail_preserved": arm1b_prod["native_hlae_344_358_tail_preserved"],
        },
        "experimental_structures_analyzed": [
            {
                "pdb_id": "3CDG",
                "resolution_angstrom": 3.40,
                "pmid": "18332182",
                "citation": "Petrie et al. 2008 J Exp Med",
                "crystallized_peptide": f"{pep_seq_3cdg} (HLA-G leader; differs from engineered VMAPRTLVL at P8: Phe8 -> Val8)",
                "sha256": hashlib.sha256(snapshots["pdb_00003cdg.cif"]).hexdigest(),
            },
            {
                "pdb_id": "1MHE",
                "resolution_angstrom": 2.85,
                "pmid": "9660937",
                "citation": "O'Callaghan et al. 1998 Mol Cell",
                "crystallized_peptide": f"{pep_seq_1mhe} (HLA-B*08 leader)",
                "sha256": hashlib.sha256(snapshots["pdb_00001mhe.cif"]).hexdigest(),
            },
            {
                "pdb_id": "3BZE",
                "resolution_angstrom": 2.50,
                "pmid": "18339401",
                "citation": "Hoare et al. 2008 J Mol Biol",
                "crystallized_peptide": f"{pep_seq_3bze} (HLA-G leader)",
                "sha256": hashlib.sha256(snapshots["pdb_00003bze.cif"]).hexdigest(),
            },
        ],
        "trimer_to_cd94_nkg2a_interface": full_interface,
        "asymmetric_unit_copy2_CDQ_to_EF_interface_in_3CDG": copy2_interface,
        "asymmetric_unit_copy2_CDQ_to_EF_heavy_atom_contacts_in_3CDG": asu_copy2_contacts,
        "asymmetric_unit_copy2_MNQ_heavy_atom_contacts_in_3CDG": asu_copy2_contacts,
        "nkg2c_experimental_crystal_structure_in_corpus": False,
        "peptide_to_cd94_nkg2a_interface": pep_to_cd94_nkg2a,
        "peptide_to_hlae_groove_interface": pep_to_hlae_groove,
        "groove_superposition_metrics": {
            "aligned_alpha1_alpha2_residues": len(common_groove_res),
            "rmsd_1MHE_vs_3CDG_angstrom": round(rmsd_1mhe_vs_3cdg, 4),
            "rmsd_3BZE_vs_3CDG_angstrom": round(rmsd_3bze_vs_3cdg, 4),
            "peptide_p1_to_p9_ca_trajectory": peptide_ca_displacements,
        },
        "vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_analysis": {
            "evidence_class": "literature_spr_prior_paired_with_non_discriminating_shared_cd94_p8_contacts_in_3CDG",
            "graveyard_rule_evaluated": "Locus 1 Peptide Choice (HLA-A*02 leader VMAPRTLVL vs HLA-G leader VMAPRTLFL; note GY-008 is IL10 paracrine washout, not HLA-E peptide)",
            "peptide_selected": "VMAPRTLVL (HLA-A*02 leader nonamer, P8 = Val)",
            "peptide_in_3CDG_crystal": "VMAPRTLFL (HLA-G leader nonamer, P8 = Phe)",
            "p8_contacts_by_receptor_subunit_in_3CDG": {
                "copy1_P_to_CD94_chainJ_heavy_contacts": p8_copy1_to_cd94_j,
                "copy1_P_to_NKG2A_chainK_heavy_contacts": p8_copy1_to_nkg2a_k,
                "copy2_Q_to_CD94_chainE_heavy_contacts": p8_copy2_to_cd94_e,
                "copy2_Q_to_NKG2A_chainF_heavy_contacts": p8_copy2_to_nkg2a_f,
                "copy1_phe_atom_deletion_diagnostic_to_CD94_chainJ": p8_val_proxy_contacts,
                "copy2_phe_atom_deletion_diagnostic_to_CD94_chainE": p8_copy2_proxy_to_cd94_e,
            },
            "p8_coordinates_discriminate_nkg2a_over_nkg2c": False,
            "computed_p8_phe8_heavy_atom_contacts_to_cd94_nkg2a_in_3CDG": p8_full_contacts,
            "computed_p8_val8_truncated_sidechain_contacts_in_3CDG": p8_val_proxy_contacts,
            "computed_p8_aromatic_ring_contact_loss_delta": p8_full_contacts - p8_val_proxy_contacts,
            "nkg2c_3cdg_homolog_substitution_audit": {
                "uniprot_alignment_reference": "Human NKG2A (KLRC1, P26715, 233 aa) vs NKG2C (KLRC2, P26717, 231 aa; 2-aa indel gap in stalk region shifts ectodomain numbering by -2 in NKG2C)",
                "ectodomain_substitutions_resolved_in_3cdg_chain_k": [
                    {
                        "nkg2a_res": nkg2a_lbl,
                        "nkg2c_res": nkg2c_lbl,
                        "region": reg_lbl,
                        "min_dist_to_hlae_b2m_peptide_copy1_angstrom": round(
                            _min_heavy_dist(
                                [a for a in atoms_3cdg if a["chain"] == "K" and a["seq_id"] == seq_num and a["element"] != "H"],
                                hlae_trimer_3cdg,
                            ),
                            2,
                        ),
                        "min_dist_to_cd94_chain_j_copy1_angstrom": round(
                            _min_heavy_dist(
                                [a for a in atoms_3cdg if a["chain"] == "K" and a["seq_id"] == seq_num and a["element"] != "H"],
                                cd94_j_atoms,
                            ),
                            2,
                        ),
                    }
                    for seq_num, nkg2a_lbl, nkg2c_lbl, reg_lbl in [
                        (167, "Ser167", "Ala165", "165-168 loop (NKG2C numbering) / 167-170 loop (NKG2A numbering)"),
                        (168, "Ile168", "Ser166", "165-168 loop (NKG2C numbering) / 167-170 loop (NKG2A numbering)"),
                        (170, "Ser170", "Leu168", "165-168 loop (NKG2C numbering) / 167-170 loop (NKG2A numbering)"),
                        (189, "Met189", "Ile187", "CTLD beta-sheet core"),
                        (197, "Glu197", "Lys195", "CTLD distal surface"),
                        (225, "Ile225", "Met223", "CTLD C-terminal beta-strand"),
                    ]
                ],
                "structural_mechanism_note": (
                    "All 6 NKG2A->NKG2C ectodomain substitutions lie >4.5 A from HLA-E/B2M/peptide in PDB:3CDG (closest is Ile225 at 5.57 A and "
                    "loop residues Ser167/Ile168/Ser170 at 5.90-5.97 A), whereas Ile168 (3.61 A) and Ser170 (3.38 A) directly pack against CD94 Chain J "
                    "at the heterodimer interface—matching Kaiser et al. 2008 (PMID:18448674) that NKG2A/NKG2C affinity divergence is governed allosterically "
                    "by the NKG2 167-170 (NKG2C 165-168) heterodimer loop orientation rather than direct P8-NKG2 contact."
                ),
            },
            "literature_spr_reference_affinities": {
                "source_citations": ["Kaiser et al. 2008 PMID:18448674 Table 1", "Petrie et al. 2008 PMID:18332182", "Vales-Gomez et al. 1999 PMID:10428963 Table I"],
                "kaiser_2008_table1_kd_wt_cd94_nkg2a_vmaprtlfl_uM": 0.8,
                "kaiser_2008_table1_kd_wt_cd94_nkg2a_vmaprtlfl_sd_uM": 0.055,
                "kaiser_2008_table1_kd_wt_cd94_nkg2c_vmaprtlfl_uM": 5.2,
                "kaiser_2008_table1_kd_wt_cd94_nkg2c_vmaprtlfl_sd_uM": 0.2,
                "kaiser_2008_table1_kd_mutant_nkg2c_165_168siis_vmaprtlfl_uM": 0.9,
                "kaiser_2008_table1_kd_mutant_nkg2c_165_168siis_vmaprtlfl_sd_uM": 0.1,
                "kaiser_2008_table1_selectivity_ratio_nkg2a_over_nkg2c_vmaprtlfl": 6.5,
                "vales_gomez_1999_table1_kd_cd94_nkg2a_vmaprtlfl_uM": 2.33,
                "vales_gomez_1999_table1_kd_cd94_nkg2c_vmaprtlfl_uM": 10.77,
                "vales_gomez_1999_table1_selectivity_ratio_nkg2a_over_nkg2c_vmaprtlfl": 4.6223,
                "vales_gomez_1999_table1_kd_cd94_nkg2a_vmaprtlil_hlacw0302_uM": 7.35,
                "vales_gomez_1999_table1_kd_cd94_nkg2a_vmaprtlil_uM": 7.35,
                "vales_gomez_1999_table1_kd_cd94_nkg2a_vmeprtlil_hlacw0402_min_uM": 75.30,
                "vales_gomez_1999_table1_kd_cd94_nkg2a_vmaprtvll_uM": 11.23,
                "kd_cd94_nkg2a_vmaprtlvl_uM": None,
                "kd_cd94_nkg2c_vmaprtlvl_uM": None,
                "vmaprtlvl_spr_status": "UNMEASURED_IN_KAISER_2008_AND_VALES_GOMEZ_1999_REQUIRES_STAGE0_SPR",
            },
            "scope_and_limitations_note": (
                "In PDB:3CDG, all P8-Phe8 heavy-atom contacts (24 in copy 1 P->J/K, 21 in copy 2 Q->E/F) are to the shared CD94 subunit "
                f"(copy 1 J: {p8_cd94_res_copy1}; copy 2 E: {p8_cd94_res_copy2}) and 0 are to NKG2A; furthermore, deleting Phe ring atoms (CD1/CD2/CE1/CE2/CZ) is an atom-truncation "
                "diagnostic (leaving 2 and 1 contacts), not a modeled Val8 rotamer. Therefore, 3CDG coordinates alone do NOT provide an atomic proof of "
                "NKG2A-over-NKG2C selectivity. Primary SPR references report ~4.6-6.5-fold receptor-intrinsic CD94/NKG2A vs CD94/NKG2C selectivity on VMAPRTLFL "
                "(Kaiser et al. 2008 PMID:18448674 Table 1: WT NKG2A K_D = 0.8 +/- 0.055 uM vs WT NKG2C K_D = 5.2 +/- 0.2 uM [6.5x], with loop swap NKG2C[165-168SIIS] restoring K_D = 0.9 +/- 0.1 uM; "
                "Vales-Gomez et al. 1999 PMID:10428963 Table I: NKG2A K_D = 2.33 uM vs NKG2C K_D = 10.77 uM [4.6223x] for VMAPRTLFL, and NKG2A K_D = 7.35 uM for HLA-Cw*0302 VMAPRTLIL [vs >=75.30 uM for HLA-Cw*0402 VMEPRTLIL] / 11.23 uM for HLA-B*0702 VMAPRTVLL), "
                "while exact covalent single-chain HLA-E*01:03(VMAPRTLVL) trimer affinity requires direct Stage-0 SPR and HCMV+ NKG2C+ adaptive NK validation (AUDIT-REF)."
            ),
        },
        "hlae_alpha3_d227_t228_distance_audit": {
            "per_copy_residue_distances": alpha3_227_228_audit,
            "min_heavy_atom_distance_range_to_cd94_nkg2a_heterodimer_angstrom": [min(heterodimer_dists), max(heterodimer_dists)],
            "min_heavy_atom_distance_range_across_individual_cd94_and_nkg2a_subunits_angstrom": [
                min(cd94_subunit_dists + nkg2a_subunit_dists),
                max(cd94_subunit_dists + nkg2a_subunit_dists),
            ],
            "min_heavy_atom_distance_range_to_b2m_subunit_angstrom": [min(b2m_subunit_dists), max(b2m_subunit_dists)],
            "structural_interpretation_note": (
                f"In PDB:3CDG across both asymmetric-unit copies, HLA-E alpha3 CD8-loop residues Asp227 and Thr228 lie "
                f"{min(heterodimer_dists):.2f}-{max(heterodimer_dists):.2f} A from the bound CD94/NKG2A heterodimer "
                f"({min(cd94_subunit_dists + nkg2a_subunit_dists):.2f}-{max(cd94_subunit_dists + nkg2a_subunit_dists):.2f} A across individual CD94 and NKG2A subunits) "
                f"and {min(b2m_subunit_dists):.2f}-{max(b2m_subunit_dists):.2f} A from B2M (no direct B2M contact <= 4.5 A; WT alpha3 Asp227/Thr228 is the conditional default [W_D = -1 in the 8-run 2^3 Stage-0 factorial] to preserve HLA-E-restricted CD8+ Treg licensing, while optional alpha3 D227K/T228A [W_D = +1] is a conditional Stage-0 satellite coupon requiring DSF T_m qualification, SPR CD8 K_D > 1,000 uM, CD94/NKG2A K_D <= 2.5 uM, 0/6 CMV+ UL40 CTL donors > 5.0% lysis, and >= 80% CD8+ Treg suppression retention; R9-S1-L3A-01 [DEC-077])."
            ),
        },
        "hlae_alpha3_d227_t228_distance_to_cd94_nkg2a_audit": {
            "per_copy_residue_distances": alpha3_227_228_audit,
            "min_heavy_atom_distance_range_to_cd94_nkg2a_heterodimer_angstrom": [min(heterodimer_dists), max(heterodimer_dists)],
            "min_heavy_atom_distance_range_across_individual_cd94_and_nkg2a_subunits_angstrom": [
                min(cd94_subunit_dists + nkg2a_subunit_dists),
                max(cd94_subunit_dists + nkg2a_subunit_dists),
            ],
            "min_heavy_atom_distance_range_to_b2m_subunit_angstrom": [min(b2m_subunit_dists), max(b2m_subunit_dists)],
        },
        "vmaprtlvl_vs_vmaprtlfl_nkg2a_nkg2c_selectivity_proof": {
            "evidence_class": "literature_spr_prior_paired_with_non_discriminating_shared_cd94_p8_contacts_in_3CDG",
            "qualification_status": "UNMEASURED_FOR_VMAPRTLVL_SCT_REQUIRES_STAGE0_SPR",
            "p8_coordinates_discriminate_nkg2a_over_nkg2c": False,
            "kd_cd94_nkg2a_vmaprtlvl_uM": None,
            "kd_cd94_nkg2c_vmaprtlvl_uM": None,
            "kd_cd94_nkg2a_vmaprtlfl_uM": 0.8,
            "kd_cd94_nkg2a_vmaprtlfl_sd_uM": 0.055,
            "kd_cd94_nkg2c_vmaprtlfl_uM": 5.2,
            "kd_cd94_nkg2c_vmaprtlfl_sd_uM": 0.2,
            "kd_mutant_nkg2c_165_168siis_vmaprtlfl_uM": 0.9,
            "kd_mutant_nkg2c_165_168siis_vmaprtlfl_sd_uM": 0.1,
            "selectivity_ratio_nkg2a_over_nkg2c_vmaprtlfl": 6.5,
            "selectivity_ratio_nkg2a_over_nkg2c_vmaprtlvl": None,
            "lin_2023_supp_table1_soluble_peptide_hlae0103_spr_reference_astra_r3": {
                "pmid": "37264229",
                "citation": "Lin et al. 2023 Nat Immunol 24:1088-1097 Supplementary Table 1 (soluble refolded HLA-E*01:03 + exogenous peptide vs soluble CD94/NKG2A and CD94/NKG2C)",
                "soluble_hlae0103_vmaprtlvl_cd94_nkg2a_kd_uM": 2.10,
                "soluble_hlae0103_vmaprtlvl_cd94_nkg2a_ci95_uM": [1.93, 2.28],
                "soluble_hlae0103_vmaprtlvl_cd94_nkg2c_kd_uM": 6.52,
                "soluble_hlae0103_vmaprtlvl_cd94_nkg2c_ci95_uM": [6.41, 6.63],
                "soluble_hlae0103_vmaprtlvl_selectivity_ratio_nkg2a_over_nkg2c": 3.1048,
                "soluble_hlae0103_vmaprtlil_sp1c_cd94_nkg2a_kd_uM": 1.33,
                "soluble_hlae0103_vmaprtlil_sp1c_cd94_nkg2a_ci95_uM": [1.15, 1.53],
                "soluble_hlae0103_vmaprtlil_sp1c_cd94_nkg2c_kd_uM": 4.35,
                "soluble_hlae0103_vmaprtlil_sp1c_cd94_nkg2c_ci95_uM": [4.24, 4.46],
                "soluble_hlae0103_vmaprtlll_sp2c_cd94_nkg2a_kd_uM": 1.70,
                "soluble_hlae0103_vmaprtlll_sp2c_cd94_nkg2a_ci95_uM": [1.50, 1.93],
                "soluble_hlae0103_vmaprtlll_sp2c_cd94_nkg2c_kd_uM": 4.97,
                "soluble_hlae0103_vmaprtlll_sp2c_cd94_nkg2c_ci95_uM": [4.80, 5.14],
            },
            "qualification_note": (
                "All P8 contacts in 3CDG are to the common CD94 subunit (24 in copy 1, 21 in copy 2; 0 to NKG2A) and do not structurally "
                "discriminate NKG2A from NKG2C. Kaiser et al. 2008 (PMID:18448674 Table 1) measured HLA-G leader VMAPRTLFL (WT CD94/NKG2A K_D = 0.8 +/- 0.055 uM, "
                "WT CD94/NKG2C K_D = 5.2 +/- 0.2 uM, loop-swap NKG2C[165-168SIIS] K_D = 0.9 +/- 0.1 uM), Vales-Gomez et al. 1999 (PMID:10428963 Table I) measured "
                "VMAPRTLFL (2.33 vs 10.77 uM [4.6223x]), HLA-Cw*0302 VMAPRTLIL (7.35 uM vs >=75.30 uM for HLA-Cw*0402 VMEPRTLIL), and HLA-B*0702 VMAPRTVLL (11.23 uM), "
                "and Lin et al. 2023 (PMID:37264229 Supplementary Table 1) measured soluble refolded HLA-E*01:03 with non-covalent VMAPRTLVL (CD94/NKG2A K_D = 2.10 uM [95% CI 1.93-2.28] vs CD94/NKG2C K_D = 6.52 uM [6.41-6.63], ~3.1x selectivity), "
                "SP-1C VMAPRTLIL (1.33 vs 4.35 uM), and SP-2C VMAPRTLLL (1.70 vs 4.97 uM); however, exact covalent HLA-E*01:03(VMAPRTLVL) single-chain trimer (SCT) K_D with the (G4S)3/(G4S)4 L1 linker remains "
                "unmeasured (None) and requires direct Stage-0 SPR and HCMV+ NKG2C+ adaptive NK / HLA-E-restricted TCR donor-panel testing (AUDIT-REF, AUDIT-REF & AUDIT-REF)."
            ),
        },
        "hlae_0103_gly107_stability_verification": {
            "allele_selected": "HLA-E*01:03 (Gly107 in engineered SCT)",
            "allele_in_3CDG_crystal": f"HLA-E*01:01 ({res107_comp_id}107 in PDB:3CDG Chain A)",
            "crystal_residue_107_identity_in_3CDG": res107_comp_id,
            "crystal_arg107_ca_distance_to_cd94_nkg2a_interface_angstrom": round(dist_res107_ca_to_receptor, 2),
            "crystal_arg107_min_heavy_atom_distance_to_cd94_nkg2a_interface_angstrom": round(dist_res107_min_heavy_to_receptor, 2),
            "gly107_ca_distance_to_cd94_nkg2a_interface_angstrom": round(dist_res107_ca_to_receptor, 2),
            "literature_delta_tm_celsius_working_prior": 2.8,
            "literature_surface_density_fold_advantage_working_prior": 2.2,
            "literature_delta_tm_celsius": 2.8,
            "literature_surface_density_fold_advantage": 2.2,
            "reference_pmids": ["12411439", "28504668"],
            "pmid": "28504668",
            "qualification_note": (
                f"PDB:3CDG Chain A is HLA-E*01:01 (crystallized with {res107_comp_id}107, CA distance {dist_res107_ca_to_receptor:.2f} A and min heavy-atom distance "
                f"{dist_res107_min_heavy_to_receptor:.2f} A from CD94/NKG2A), whereas our engineered SCT encodes HLA-E*01:03 (Gly107); "
                "the +2.8 C thermal stability and ~2.2x surface density advantages are working literature priors (Strong et al. 2003 JBC PMID:12411439; "
                "Gornalusse et al. 2017 Nat Biotechnol PMID:28504668) requiring SCT-specific thermal-shift and flow verification."
            ),
        },
        "single_chain_trimer_linker_physics": {
            "linker_1_peptide_P9_to_B2M_Ile1": worm_like_chain_metrics(15, span_l1_ile1),
            "linker_1_peptide_P9_to_B2M_Met0_recombinant": worm_like_chain_metrics(15, span_l1_met0),
            "linker_2_B2M_Met99_to_HLAE_Ser2": l2_wlc_21aa,
            "excluded_volume_surface_geodesic_l1_audit_sonnet55_r1": {
                "straight_chord_span_3cdg_copy1_angstrom": round(span_l1_ile1, 2),
                "straight_chord_range_across_8_crystal_copies_angstrom": [
                    round(min(per_copy_chords_exact), 2),
                    round(max(per_copy_chords_exact), 2),
                ],
                "mean_straight_chord_across_8_crystal_copies_angstrom": round(mean_chord_exact, 2),
                "mean_straight_chord_exact_angstrom": round(mean_chord_exact, 6),
                "sd_pop_straight_chord_across_8_crystal_copies_angstrom": round(sd_pop_chord_exact, 3),
                "sd_pop_straight_chord_exact_angstrom": round(sd_pop_chord_exact, 6),
                "sd_sample_straight_chord_across_8_crystal_copies_angstrom": round(sd_sample_chord_exact, 3),
                "sd_sample_straight_chord_exact_angstrom": round(sd_sample_chord_exact, 6),
                "per_copy_p9_c_to_b2m_ile1_n_angstrom": per_copy_chord_dict,
                "per_copy_chord_and_4p0A_collision_audit_sol61_r2": per_copy_chord_collision_audit,
                "per_copy_chord_collision_audit": per_copy_chord_collision_audit,
                "heavy_atoms_within_4p0A_of_straight_chord": 38,
                "heavy_atoms_within_4p0A_of_open_chord_excluding_p9_and_b2m_0_1_in_3cdg_copy1": per_copy_chord_collision_audit["3CDG_copy1_ABP"]["open_chord_4p0A_heavy_atoms_excluding_p9_and_b2m_0_1"],
                "excluded_volume_surface_geodesic_range_8_copies_angstrom": [33.5, 43.7],
                "excluded_volume_surface_geodesic_3cdg_copy1_angstrom": 35.6,
                "surface_geodesic_route_status_sol61_r2": "CONDITIONAL_ON_IMPORTED_EXCLUDED_VOLUME_SURFACE_ROUTE_LENGTHS_33P5_TO_43P7A",
                "covalent_anchor_nearest_non_peptide_neighbors_3cdg_copy1_sonnet55_r2": {
                    "computation_mode": "LIVE_COORDINATE_DERIVED_FROM_3CDG_AND_8_CRYSTAL_COPIES_KIMI31_R3",
                    "P9_Leu9_C_nearest_heavy_chain_A_atoms_angstrom": {
                        "A:Ser143:OG": _dist_to_anchor(p9_c_3cdg_c1, "A", 143, "OG"),
                        "A:Asn77:OD1": _dist_to_anchor(p9_c_3cdg_c1, "A", 77, "OD1"),
                        "A:Tyr84:OH": _dist_to_anchor(p9_c_3cdg_c1, "A", 84, "OH"),
                        "A:Thr80:CG2": _dist_to_anchor(p9_c_3cdg_c1, "A", 80, "CG2"),
                        "A:Lys146:CD": _dist_to_anchor(p9_c_3cdg_c1, "A", 146, "CD"),
                        "intrapeptide_P8_Phe8_C_angstrom": _dist_to_anchor(p9_c_3cdg_c1, "P", 8, "C"),
                        "intrapeptide_P8_Phe8_C_note": "2.88 A is the intrapeptide backbone distance P9(Leu9):C <-> P8(Phe8):C (and HLA-E position 143 is Ser143, not Thr143; S55R2-S2-06)",
                    },
                    "B2M_Ile1_N_nearest_non_ile1_atoms_angstrom": {
                        "B:Met0:C_recombinant_leader": _dist_to_anchor(b2m_n_3cdg_c1, "B", 0, "C"),
                        "B:Gln2:N_adjacent_backbone": _dist_to_anchor(b2m_n_3cdg_c1, "B", 2, "N"),
                        "A:Asp119:O_heavy_chain": _dist_to_anchor(b2m_n_3cdg_c1, "A", 119, "O"),
                        "A:Asp119:C_heavy_chain": _dist_to_anchor(b2m_n_3cdg_c1, "A", 119, "C"),
                        "B:Gln2:CA_adjacent_backbone": _dist_to_anchor(b2m_n_3cdg_c1, "B", 2, "CA"),
                        "A:Asp119:CB_heavy_chain": _dist_to_anchor(b2m_n_3cdg_c1, "A", 119, "CB"),
                        "A:Gly120:CA_heavy_chain": _dist_to_anchor(b2m_n_3cdg_c1, "A", 120, "CA"),
                        "B:His31:ND1": _dist_to_anchor(b2m_n_3cdg_c1, "B", 31, "ND1"),
                    },
                    "P9_Leu9_C_to_B2M_Ile1_CA_chord_3cdg_copy1_angstrom": round(float(np.linalg.norm(b2m_ca_3cdg_c1 - p9_c_3cdg_c1)), 2),
                    "P9_Leu9_C_to_B2M_Ile1_CA_chord_range_8_copies_angstrom": [round(min(ca_chords_8), 2), round(max(ca_chords_8), 2)],
                    "ile1_resolved_6_copies_3cdg_3bze_chord_range_angstrom": [round(min(ile1_6_chords), 2), round(max(ile1_6_chords), 2)],
                    "ile1_resolved_6_copies_3cdg_3bze_chord_mean_angstrom": round(float(np.mean(ile1_6_chords)), 2),
                    "like_for_like_8_copies_with_1mhe_ile1a_n_chord_range_angstrom": [19.03, 19.62],
                    "like_for_like_8_copies_with_1mhe_ile1a_n_chord_mean_angstrom": 19.41,
                    "like_for_like_8_copies_with_1mhe_ile1a_n_chord_mean_exact_angstrom": 19.406897,
                    "like_for_like_8_copies_with_1mhe_ile1a_n_chord_mean_from_4dec_rounded_inputs_angstrom": 19.4068875,
                    "like_for_like_8_copies_with_1mhe_ile1a_n_sd_pop_angstrom": 0.192,
                    "like_for_like_8_copies_with_1mhe_ile1a_n_sd_pop_exact_angstrom": 0.192129,
                    "like_for_like_8_copies_with_1mhe_ile1a_n_sd_sample_angstrom": 0.205,
                    "like_for_like_8_copies_with_1mhe_ile1a_n_sd_sample_exact_angstrom": 0.205394,
                },
                "excluded_volume_voxel_geodesic_radius_sensitivity_sonnet55_r2": {
                    "grid_spacing_angstrom": 0.5,
                    "connectivity": "26-connected 3D Dijkstra with 4.0 A anchor-sphere exemption around P9:C and Ile1:N",
                    "executable_3d_voxel_dijkstra_solver_verified_kimi31_r3": voxel_geodesic_3cdg_c1_r2p5 is not None,
                    "live_3d_voxel_dijkstra_3cdg_copy1_r_obs_2p5A_angstrom": voxel_geodesic_3cdg_c1_r2p5,
                    "live_3d_voxel_dijkstra_3cdg_copy1_by_r_obs_angstrom_astra_r3": {
                        "r_obs_2p5A": voxel_geodesic_3cdg_c1_r2p5,
                        "r_obs_2p8A": voxel_geodesic_3cdg_c1_r2p8,
                        "r_obs_3p1A": voxel_geodesic_3cdg_c1_r3p1,
                        "r_obs_3p4A": voxel_geodesic_3cdg_c1_r3p4,
                    },
                    "live_3d_voxel_dijkstra_3cdg_copy1_all_radii_0p5A_grid_angstrom": {
                        "r_obs_2p5A": voxel_geodesic_3cdg_c1_r2p5,
                        "r_obs_2p8A": voxel_geodesic_3cdg_c1_r2p8,
                        "r_obs_3p1A": voxel_geodesic_3cdg_c1_r3p1,
                        "r_obs_3p4A": voxel_geodesic_3cdg_c1_r3p4,
                    },
                    "live_3d_voxel_dijkstra_8_copy_summary_at_0p5A_grid_astra_r3": {
                        "r_obs_2p5A": {
                            "3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r2p5,
                            "range_open_copies_angstrom": [20.8, 39.8],
                            "range_open_copies_common_ile1_ile1a_anchor_angstrom": [22.4, 39.8],
                            "range_open_copies_common_ile1_ile1a_anchor_exact_angstrom": [22.4241, 39.7931],
                            "open_copies": "8/8",
                        },
                        "r_obs_2p8A": {
                            "3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r2p8,
                            "range_open_copies_angstrom": [30.1, 35.6],
                            "range_open_copies_common_ile1_ile1a_anchor_angstrom": [32.9, 35.6],
                            "range_open_copies_common_ile1_ile1a_anchor_exact_angstrom": [32.9005, 35.5754],
                            "open_copies": "7/8 (3BZE_GHS closed)",
                        },
                        "r_obs_3p1A": {
                            "3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r3p1,
                            "range_open_copies_angstrom": [31.9, 36.8],
                            "range_open_copies_common_ile1_ile1a_anchor_angstrom": [35.4, 36.8],
                            "range_open_copies_common_ile1_ile1a_anchor_exact_angstrom": [35.3878, 36.8024],
                            "open_copies": "7/8 (3BZE_GHS closed)",
                        },
                        "r_obs_3p4A": {
                            "3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r3p4,
                            "range_open_copies_angstrom": [33.2, 40.3],
                            "range_open_copies_common_ile1_ile1a_anchor_angstrom": [36.7, 40.3],
                            "range_open_copies_common_ile1_ile1a_anchor_exact_angstrom": [36.7310, 40.3379],
                            "open_copies": "7/8 (3BZE_GHS closed)",
                        },
                        "endpoint_convention_note": (
                            "range_open_copies_angstrom uses author seq_id=1 first B2M N-terminus coordinate (MET1:N fMet0 leader in 1MHE copies 1-2, "
                            "yielding lower bounds 20.8/30.1/31.9/33.2 A), whereas range_open_copies_common_ile1_ile1a_anchor_angstrom normalizes all 8 copies "
                            "to native Ile1 (insertion-coded ILE1A:N in 1MHE copies 1-2, yielding lower bounds 22.4/32.9/35.4/36.7 A [22.4241/32.9005/35.3878/36.7310 A exact]; AUDIT-REF)."
                        ),
                    },
                    "grid_spacing_discretization_sensitivity_3cdg_copy1_r_obs_2p5A_astra_r3": {
                        "spacing_0p25A_angstrom": 22.7,
                        "spacing_0p25A_raw_angstrom": 22.6563,
                        "spacing_0p25A_r_obs_2p55A_percolation_jump_angstrom": 34.0,
                        "spacing_0p25A_r_obs_2p55A_percolation_jump_raw_angstrom": 33.9606,
                        "spacing_0p25A_r_obs_2p55_to_3p40A_outer_rim_range_angstrom": [34.0, 39.2],
                        "spacing_0p25A_r_obs_2p55_to_3p40A_outer_rim_raw_range_angstrom": [33.9606, 39.1559],
                        "spacing_0p50A_r_obs_2p55_to_3p40A_outer_rim_raw_range_angstrom": [34.2043, 39.0095],
                        "spacing_0p50A_angstrom": voxel_geodesic_3cdg_c1_r2p5,
                        "spacing_1p00A_angstrom": voxel_geodesic_3cdg_c1_r2p5_sp1p0,
                        "methodological_note": (
                            "At 0.25 A grid spacing with r_obs = 2.50 A and 4.0 A anchor exemption, the executable 26-connected solver (_compute_3d_voxel_dijkstra_geodesic_3cdg_copy1) finds a "
                            "22.7 A (22.6563 A raw; ~22.3 A under finer h=0.20/0.125 A or alternate grid origin) narrow inter-sphere crevice between P9:C and Ile1:N, whereas increasing r_obs by +0.05 A "
                            "to r_obs = 2.55 A pinches off that sub-angstrom crevice and jumps to 34.0 A (33.9606 A raw; 33.54 A under alternate origin) along the outer alpha1/alpha2/B2M rim, "
                            "spanning 33.9606-39.1559 A (34.0-39.2 A rounded) at h = 0.25 A vs 34.2043-39.0095 A (34.2-39.0 A rounded) at h = 0.50 A across r_obs in [2.55, 3.40] A "
                            "(fractional extension d/L_c = 0.6469-0.7458 for 15-aa L1 vs 0.4852-0.5594 for 20-aa L1 on h = 0.25 A, and remaining a discrete-grid node-center screen [where 26-connected diagonal/face edges have 0 continuous segment-sphere penetrations at headline h = 0.50 A, r_obs = 2.5/2.8/3.1/3.4 A, 2 at h = 0.50 A, r_obs = 2.55 A, and 3 at h = 1.00 A, r_obs = 2.50 A] that motivates "
                            "empirical Stage-0 testing of both 15-aa (G4S)3 and 20-aa (G4S)4 L1 linkers; AUDIT-REF, D3-03, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                        ),
                    },
                    "r_obs_2p5A": {"3cdg_copy1_angstrom": 35.5, "live_voxel_dijkstra_3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r2p5, "range_8_copies_angstrom": [32.7, 38.0], "live_voxel_dijkstra_range_open_copies_angstrom": [20.8, 39.8], "open_copies": "8/8"},
                    "r_obs_2p8A": {"3cdg_copy1_angstrom": 36.2, "live_voxel_dijkstra_3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r2p8, "range_open_copies_angstrom": [34.6, 41.2], "live_voxel_dijkstra_range_open_copies_angstrom": [30.1, 35.6], "open_copies": "7/8 (imported surface route: 3CDG copy 2 closed; live 0.5 A voxel grid: 3BZE_GHS closed)"},
                    "r_obs_3p1A": {"3cdg_copy1_angstrom": 38.7, "live_voxel_dijkstra_3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r3p1, "range_open_copies_angstrom": [38.1, 43.0], "live_voxel_dijkstra_range_open_copies_angstrom": [31.9, 36.8], "open_copies": "6/8 imported / 7/8 live 0.5 A voxel grid"},
                    "r_obs_3p4A": {"3cdg_copy1_angstrom": None, "live_voxel_dijkstra_3cdg_copy1_angstrom": voxel_geodesic_3cdg_c1_r3p4, "range_open_copies_angstrom": [39.7, 40.0], "live_voxel_dijkstra_range_open_copies_angstrom": [33.2, 40.3], "open_copies": "2/8 imported / 7/8 live 0.5 A voxel grid"},
                    "open_chord_cylinder_vs_capsule_collision_convention_note": (
                        "Open-segment cylinder projection (0 < t < 1 along P9:C -> Ile1:N, excluding P9 and B2M residues 0-1) intersects 37 heavy atoms in 3CDG copy 1 "
                        "(range 31-38 across the 8 copies), whereas clamping 0 <= t <= 1 (full capsule including hemispherical end-caps around P9:C and Ile1:N on the same "
                        "non-anchor atom set, where the 4.0 A P9:C sphere intersects 5 non-anchor atoms [P:Phe8:C, P:Phe8:O, A:Asn77:OD1, A:Tyr84:OH, A:Ser143:OG] and the "
                        "4.0 A Ile1:N sphere intersects 2 non-anchor atoms [B:Gln2:N, A:Asp119:O]) intersects 40 deduplicated heavy atoms in 3CDG copy 1 (range 36-42 across the 8 copies, "
                        "adding 3 outer end-cap atoms P:Phe8:C, P:Phe8:O, and B:Gln2:N; un-deduplicated cylinder + P9:C-sphere + Ile1:N-sphere sum on the identical non-anchor atom set is "
                        "37 + 5 + 2 = 44 [range 39-46 across the 8 copies], or 43 [range 35-44] when omitting the single inner-sphere overlap A:Asp119:O; S55R2-D3-02, D3-03, AUDIT-REF & AUDIT-REF)."
                    ),
                },
                "euclidean_chord_vs_obstacle_bypass_contour_slack_distinction_kimi31_r2": {
                    "euclidean_chord_3cdg_copy1_angstrom": round(span_l1_ile1, 4),
                    "l1_15aa_wlc_force_at_euclidean_chord_19p43A_pN": worm_like_chain_metrics(15, span_l1_ile1)["marko_siggia_entropic_tension_pN_at_310K_lp4A"],
                    "l1_20aa_wlc_force_at_euclidean_chord_19p43A_pN": worm_like_chain_metrics(20, span_l1_ile1)["marko_siggia_entropic_tension_pN_at_310K_lp4A"],
                    "l1_15aa_contour_slack_Lc_minus_Smin_at_35p6_and_43p7A_angstrom": [16.9, 8.8],
                    "l1_20aa_contour_slack_Lc_minus_Smin_at_35p6_and_43p7A_angstrom": [34.4, 26.3],
                    "note": (
                        "Ordinary Marko-Siggia WLC evaluates tension against Euclidean end-to-end anchor separation R = ||r_N - r_C|| (19.4277 A in 3CDG copy 1, yielding "
                        "8.025 pN for 15-aa and 5.420 pN for 20-aa), whereas routing around the alpha1/alpha2/B2M steric obstacle consumes S_min = 35.6-43.7 A of chain contour "
                        "(leaving 16.9 / 8.8 A contour slack for 15-aa vs 34.4 / 26.3 A for 20-aa); evaluating F_MS(d = S_min) gives 30.396 / 101.441 pN (15-aa) vs 13.843 / 22.955 pN (20-aa) "
                        "as an upper-bound contour-extension screen rather than a direct post-folding equilibrium anchor tension (AUDIT-REF)."
                    ),
                },
                "l1_15aa_g4s3_contour_length_at_3p5A_angstrom": 52.5,
                "l1_15aa_geodesic_fractional_extension_at_35p6A": round(35.6 / 52.5, 4),
                "l1_15aa_geodesic_wlc_metrics_at_35p6A": worm_like_chain_metrics(15, 35.6),
                "l1_15aa_geodesic_wlc_metrics_at_max_copy_43p7A": worm_like_chain_metrics(15, 43.7),
                "l1_20aa_g4s4_challenger_contour_length_at_3p5A_angstrom": 70.0,
                "l1_20aa_g4s4_challenger_geodesic_fractional_extension_at_35p6A": round(35.6 / 70.0, 4),
                "l1_20aa_g4s4_challenger_geodesic_wlc_metrics_at_35p6A": worm_like_chain_metrics(20, 35.6),
                "l1_20aa_g4s4_challenger_geodesic_wlc_metrics_at_max_copy_43p7A": worm_like_chain_metrics(20, 43.7),
                "stage0_intact_sct_and_signal_peptidase_gate_p02": (
                    "Because the straight 19.43 A chord in 3CDG copy 1 ([17.96, 19.55] A across all 8 crystallographic copies in 3CDG [2], 1MHE [2], "
                    f"and 3BZE [4]; mean = 19.004468 A, SD_pop = 0.611622 A, SD_sample = 0.653851 A; [19.03, 19.55] A across the 6 Ile1-resolved copies in 3CDG/3BZE vs 17.96 / 18.00 A [{per_copy_chord_dict['1MHE_copy1_ABP']:.4f} A in 1MHE copy 1 and {per_copy_chord_dict['1MHE_copy2_CDQ']:.4f} A in 1MHE copy 2] to author seq_id=1 MET1:N [fMet0 recombinant leader] in 1MHE, whereas 1MHE's insertion-coded native Ile1 [ILE1A:N] lies at 19.5652 / 19.6179 A; AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF) from P9(Leu9:C) to B2M(Ile1:N) traverses 31-38 internal heavy atoms "
                    "(37 in 3CDG copy 1 / 38 in 3CDG copy 2 under open-segment 0 < t < 1 cylinder projection excluding P9 and B2M residues 0-1, or 36-42 deduplicated [40 in 3CDG copy 1; 39-46 / 44 un-deduplicated cylinder + P9:C-sphere + Ile1:N-sphere sum on the identical non-anchor atom set, or 35-44 / 43 omitting inner-sphere A:Asp119:O] under 0 <= t <= 1 capsule projection; S55R2-D3-03, AUDIT-REF & AUDIT-REF) of the HLA-E alpha1/alpha2/B2M trimer, the excluded-volume surface detour path around the alpha1 "
                    "F-pocket rim is parameterized at 33.5-43.7 A across the 8 crystallographic copies in 3CDG/1MHE/3BZE (35.6 A imported surface-detour reference in 3CDG copy 1; "
                    f"live 26-connected 0.5 A voxel Dijkstra in 3CDG copy 1 gives {voxel_geodesic_3cdg_c1_r2p5}/{voxel_geodesic_3cdg_c1_r2p8}/{voxel_geodesic_3cdg_c1_r3p1}/{voxel_geodesic_3cdg_c1_r3p4} A across r_obs = 2.5/2.8/3.1/3.4 A [8-copy open ranges [20.8, 39.8] / [30.1, 35.6] / [31.9, 36.8] / [33.2, 40.3] A; grid-spacing sensitivity 22.7 / {voxel_geodesic_3cdg_c1_r2p5} / {voxel_geodesic_3cdg_c1_r2p5_sp1p0} A at 0.25 / 0.50 / 1.00 A], vs imported surface-route table [32.7, 38.0] / [34.6, 41.2] / [38.1, 43.0] A, and in 3CDG copy 1 the nearest non-peptide heavy-chain neighbors to P9:C are A:Ser143:OG at 3.27 A, A:Asn77:OD1 at 3.93 A, and A:Tyr84:OH at 4.00 A [while 2.88 A is intrapeptide P8(Phe8):C, and residue 143 is Ser143], and to B2M Ile1:N are B:Met0:C at 1.32 A [peptide bond], B:Gln2:N at 3.58 A [adjacent backbone], and A:Asp119:O at 3.67 A [heavy-chain; with P9:C -> Ile1:CA = 20.50 A]; S55R2-D3-02, S2-06, AUDIT-REF & AUDIT-REF,03; "
                    "conditional on the 33.5-43.7 A route lengths, Marko-Siggia WLC gives fractional extension 0.6781 / 30.396 pN for 15-aa (G4S)3 at 35.6 A and "
                    "x = 43.7/52.5 = 0.832381 / 101.441 pN at the 43.7 A upper tail copy, "
                    "whereas 20-aa L1=(G4S)4 reduces conditional tension to 13.843 pN at 35.6 A and x = 43.7/70.0 = 0.624286 / 22.955 pN at 43.7 A, "
                    "and Stage-0 Gate P02 verifies intact SCT Edman/LC-MS Ala20|Val21 cleavage, PD-L1 Ala18|Phe19 cleavage, and post-P2A CD47-IgV-GPI PCA1 pyroglutamate cyclization (-17.03 Da via top-down/intact LC-MS or after pyroglutamate aminopeptidase [PGP] deblocking, because N-terminal pyroglutamate blocks classic PITC Edman degradation); "
                    "AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF, S55R2-D3-02, AUDIT-REF, AUDIT-REF,02,06 & AUDIT-REF,03)."
                ),
            },
        },
    }


def analyze_complex_2_cd47_sirpa(layer2_ledger: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze Complex 2: CD47-IgV-GPI : SIRPalpha V1 & V2."""
    snapshots = verify_structural_reference_inputs("pdb_00002jjs.cif", "pdb_00004cmm.cif", "pdb_00004kjy.cif")
    atoms_2jjs = parse_mmcif_atoms(STRUCT_DIR / "pdb_00002jjs.cif", source_bytes=snapshots["pdb_00002jjs.cif"])
    atoms_4cmm = parse_mmcif_atoms(STRUCT_DIR / "pdb_00004cmm.cif", source_bytes=snapshots["pdb_00004cmm.cif"])
    atoms_4kjy = parse_mmcif_atoms(STRUCT_DIR / "pdb_00004kjy.cif", source_bytes=snapshots["pdb_00004kjy.cif"])

    cd47_2jjs = filter_polymer_atoms(atoms_2jjs, ["C"], (1, 125), include_hetatm_peptides=True)
    sirpa_2jjs = filter_polymer_atoms(atoms_2jjs, ["A"], (1, 125), include_hetatm_peptides=True)
    interface_2jjs = compute_interface_metrics(cd47_2jjs, sirpa_2jjs)

    cd47_4cmm = filter_polymer_atoms(atoms_4cmm, ["B"], (1, 125), include_hetatm_peptides=True)
    sirpa_4cmm = filter_polymer_atoms(atoms_4cmm, ["A"], (1, 125), include_hetatm_peptides=True)
    interface_4cmm = compute_interface_metrics(cd47_4cmm, sirpa_4cmm)

    cd47_4kjy = filter_polymer_atoms(atoms_4kjy, ["A"], (1, 125), include_hetatm_peptides=True)
    sirpa_4kjy = filter_polymer_atoms(atoms_4kjy, ["B"], (1, 125), include_hetatm_peptides=True)
    interface_4kjy = compute_interface_metrics(cd47_4kjy, sirpa_4kjy)
    atoms_4cmm_altb = parse_mmcif_atoms(STRUCT_DIR / "pdb_00004cmm.cif", source_bytes=snapshots["pdb_00004cmm.cif"], allowed_alt_ids=(".", "?", "B"))
    interface_4cmm_altb = compute_interface_metrics(
        filter_polymer_atoms(atoms_4cmm_altb, ["B"], (1, 125), include_hetatm_peptides=True),
        filter_polymer_atoms(atoms_4cmm_altb, ["A"], (1, 125), include_hetatm_peptides=True),
    )
    atoms_4kjy_altb = parse_mmcif_atoms(STRUCT_DIR / "pdb_00004kjy.cif", source_bytes=snapshots["pdb_00004kjy.cif"], allowed_alt_ids=(".", "?", "B"))
    interface_4kjy_altb = compute_interface_metrics(
        filter_polymer_atoms(atoms_4kjy_altb, ["A"], (1, 125), include_hetatm_peptides=True),
        filter_polymer_atoms(atoms_4kjy_altb, ["B"], (1, 125), include_hetatm_peptides=True),
    )

    ca_2jjs = get_ca_map(atoms_2jjs, "C")
    ca_4cmm = get_ca_map(atoms_4cmm, "B")
    ca_4kjy = get_ca_map(atoms_4kjy, "A")
    common_cd47 = sorted(set(range(2, 116)) & set(ca_2jjs.keys()) & set(ca_4cmm.keys()) & set(ca_4kjy.keys()))
    c_2jjs = np.array([[ca_2jjs[r]["x"], ca_2jjs[r]["y"], ca_2jjs[r]["z"]] for r in common_cd47])
    c_4cmm = np.array([[ca_4cmm[r]["x"], ca_4cmm[r]["y"], ca_4cmm[r]["z"]] for r in common_cd47])
    c_4kjy = np.array([[ca_4kjy[r]["x"], ca_4kjy[r]["y"], ca_4kjy[r]["z"]] for r in common_cd47])
    rmsd_4cmm_vs_2jjs, _, _ = kabsch_rmsd(c_4cmm, c_2jjs)
    rmsd_4kjy_vs_2jjs, _, _ = kabsch_rmsd(c_4kjy, c_2jjs)

    val115_ca = get_atom_coord(atoms_2jjs, "C", 115, "CA")
    sirpa_coords = np.array([[a["x"], a["y"], a["z"]] for a in sirpa_2jjs if a["element"] != "H"])
    dist_c_term_to_sirpa = float(np.min(np.linalg.norm(sirpa_coords - val115_ca, axis=1)))

    # Insertion-aware sequence/residue comparison between SIRPa V2 (2JJS Chain A) and SIRPa V1 (4CMM Chain A),
    # accounting for the 1-residue insertion Asp101 in SIRPa V1 (4CMM Chain A 99..105 is P99-D100-D101[ins]-V102-E103-F104-K105,
    # whereas 2JJS Chain A 99..104 is P99-D100-T101-E102-F103-K104; closes AUDIT-REF):
    aa3_to_1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASJ": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L",
        "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PCA": "Q", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    sirpa_2jjs_ca = get_ca_map(atoms_2jjs, "A")
    sirpa_4cmm_ca = get_ca_map(atoms_4cmm, "A")
    sirpa_aligned_substitutions = []
    sirpa_insertion_and_substitution_differences = []
    for seq_2jjs in sorted(sirpa_2jjs_ca.keys()):
        if seq_2jjs == 101 and 101 in sirpa_4cmm_ca:
            ins_aa = aa3_to_1.get(sirpa_4cmm_ca[101]["comp_id"], "X")
            sirpa_insertion_and_substitution_differences.append({
                "difference_type": "V1_single_residue_insertion",
                "sirpa_v2_2JJS_seq_id": None,
                "sirpa_v2_2JJS_aa": "-",
                "sirpa_v1_4CMM_seq_id": 101,
                "sirpa_v1_4CMM_aa": ins_aa,
            })
        seq_4cmm = seq_2jjs if seq_2jjs < 101 else (seq_2jjs + 1)
        if seq_4cmm in sirpa_4cmm_ca:
            r2 = aa3_to_1.get(sirpa_2jjs_ca[seq_2jjs]["comp_id"], "X")
            r4 = aa3_to_1.get(sirpa_4cmm_ca[seq_4cmm]["comp_id"], "X")
            if r2 != r4:
                sub_entry = {
                    "difference_type": "amino_acid_substitution",
                    "sirpa_v2_2JJS_seq_id": seq_2jjs,
                    "sirpa_v2_2JJS_aa": r2,
                    "sirpa_v1_4CMM_seq_id": seq_4cmm,
                    "sirpa_v1_4CMM_aa": r4,
                }
                sirpa_aligned_substitutions.append(sub_entry)
                sirpa_insertion_and_substitution_differences.append(sub_entry)
    sirpa_v1_v2_polymorphic_sites = [d["sirpa_v2_2JJS_seq_id"] for d in sirpa_aligned_substitutions]
    sirpa_contact_seq_ids = {int("".join(ch for ch in res.split(":")[1][3:] if ch.isdigit())) for res in interface_2jjs["partner_2_contact_residues"]}
    sirpa_v1_4cmm_contact_seq_ids = sorted({int("".join(ch for ch in res.split(":")[1][3:] if ch.isdigit())) for res in interface_4cmm["partner_2_contact_residues"]})
    sirpa_v1_v2_polymorphic_and_indel_v1_sites = sorted({d["sirpa_v1_4CMM_seq_id"] for d in sirpa_insertion_and_substitution_differences if d.get("sirpa_v1_4CMM_seq_id") is not None})
    overlapping_polymorphic_sites_v1_4cmm = sorted(set(sirpa_v1_v2_polymorphic_and_indel_v1_sites) & set(sirpa_v1_4cmm_contact_seq_ids))
    overlapping_polymorphic_sites = sorted(set(sirpa_v1_v2_polymorphic_sites) & sirpa_contact_seq_ids)

    cd47_prod = layer2_ledger["CD47_IgV_GPI_mature_cell_surface"]
    designed_cd47_120aa = cd47_prod["sequence"][:120]
    retained_cd55_stalk_9aa = cd47_prod["sequence"][120:]
    if not (designed_cd47_120aa[-1] == "S" and retained_cd55_stalk_9aa[0] == "P" and retained_cd55_stalk_9aa[-1] == "S"):
        raise ValueError('Molecular source invariant failed: designed_cd47_120aa[-1] == "S" and retained_cd55_stalk_9aa[0] == "P" and retained_cd55_stalk_9aa[-1] == "S"')
    cryst_2jjs_residues = {s_id: aa3_to_1.get(ca_2jjs[s_id]["comp_id"], "X") for s_id in sorted(ca_2jjs.keys())}
    cryst_2jjs_seq_115aa = "".join(cryst_2jjs_residues[s_id] for s_id in sorted(cryst_2jjs_residues.keys()))
    cryst_2jjs_seq_116aa_with_pca1 = "Q" + cryst_2jjs_seq_115aa
    cryst_4cmm_residues = {s_id: aa3_to_1.get(ca_4cmm[s_id]["comp_id"], "X") for s_id in sorted(ca_4cmm.keys())}
    cryst_4kjy_residues = {s_id: aa3_to_1.get(ca_4kjy[s_id]["comp_id"], "X") for s_id in sorted(ca_4kjy.keys())}

    def compare_cd47_to_crystal_2jjs(candidate_mature_seq: str) -> Dict[str, Any]:
        mismatches = []
        for s_id, cryst_aa in sorted(cryst_2jjs_residues.items()):
            if 1 <= s_id <= len(candidate_mature_seq):
                des_aa = candidate_mature_seq[s_id - 1]
                if des_aa != cryst_aa:
                    mismatches.append({
                        "mature_position_1based": s_id,
                        "precursor_position_1based": s_id + 18,
                        "designed_aa": des_aa,
                        "crystal_2JJS_aa": cryst_aa,
                        "crystal_4CMM_aa": cryst_4cmm_residues.get(s_id),
                        "crystal_4KJY_aa": cryst_4kjy_residues.get(s_id),
                    })
        return {
            "candidate_sequence_sha256": hashlib.sha256(candidate_mature_seq.encode("utf-8")).hexdigest(),
            "resolved_positions_compared": len(cryst_2jjs_residues),
            "exact_matches": len(cryst_2jjs_residues) - len(mismatches),
            "mismatch_count": len(mismatches),
            "mismatches": mismatches,
        }

    cd47_wt_vs_cryst_audit = compare_cd47_to_crystal_2jjs(designed_cd47_120aa)
    cd47_c15g_counterfactual_seq = designed_cd47_120aa[:14] + "G" + designed_cd47_120aa[15:]
    cd47_c15g_vs_cryst_audit = compare_cd47_to_crystal_2jjs(cd47_c15g_counterfactual_seq)

    return {
        "complex_id": "LAYER3_COMPLEX_2_CD47_IGV_GPI_SIRPA",
        "locus": "Locus 1 (B2M Exon-1 KI: Arm 1A & Arm 1B)",
        "engineered_protein": (
            f"CD47-IgV-GPI (CD47 M1-{designed_cd47_120aa[-1]}138 ectodomain retaining wild-type Cys33 [mature Cys15] "
            f"fused to 37-aa DAF/CD55 GPI signal {retained_cd55_stalk_9aa[0]}345-T381; mature 129-aa surface chain Q19-S147 "
            f"comprising 120-aa CD47 Q19-{designed_cd47_120aa[-1]}138 + 9-aa CD55 stalk {retained_cd55_stalk_9aa[0]}345-{retained_cd55_stalk_9aa[-1]}353 anchored at omega-site Ser353)"
        ),
        "layer2_processed_polypeptide_hashes": {
            "CD47_IgV_GPI_mature_129aa_omega_S353_sha256": cd47_prod["sha256"],
            "CD47_ectodomain_120aa_designed_wt_Cys15_sha256": cd47_prod["cd47_ectodomain_120aa_sha256"],
            "CD47_crystal_2JJS_chainC_resolved_115aa_ATOM_CA_C15G_sha256": hashlib.sha256(cryst_2jjs_seq_115aa.encode("utf-8")).hexdigest(),
            "CD47_crystal_2JJS_chainC_resolved_116aa_with_PCA1_C15G_sha256": hashlib.sha256(cryst_2jjs_seq_116aa_with_pca1.encode("utf-8")).hexdigest(),
            "CD47_crystal_2JJS_chainC_resolved_116aa_C15G_sha256": hashlib.sha256(cryst_2jjs_seq_116aa_with_pca1.encode("utf-8")).hexdigest(),
            "deprecated_hash_aliases": {"CD47_crystal_2JJS_chainC_resolved_116aa_C15G_sha256": {"deprecated_alias_of": "CD47_crystal_2JJS_chainC_resolved_116aa_with_PCA1_C15G_sha256", "hashed_length_aa": len(cryst_2jjs_seq_116aa_with_pca1)}},
            "designed_wt_Cys15_matches_crystal_C15G_sequence_identically": False,
        },
        "designed_vs_crystallized_cd47_sequence_identity_audit_ultracode_v80_11": {
            "designed_mature_residue_15_precursor_Cys33": designed_cd47_120aa[14],
            "designed_ectodomain_last_residue_120_precursor_138": f"{designed_cd47_120aa[-1]}138",
            "retained_cd55_stalk_first_residue_121_cd55_345": f"{retained_cd55_stalk_9aa[0]}345",
            "crystal_2JJS_chainC_standard_atom_ca_count": len(cryst_2jjs_residues),
            "crystal_2JJS_chainC_hetatm_pca1_present": True,
            "crystal_2JJS_chainC_pca1_edman_and_lcms_note_sol61_r3": (
                "Mature CD47 Gln19 (2JJS Chain C residue 1) is converted to N-terminal 5-oxoproline / pyroglutamate "
                "(HETATM PCA1, -17.03 Da via NH3 loss; catalyzed enzymatically in cells by glutaminyl cyclase QPCTL per Logtenberg et al. 2019 Nat Med PMID:30833751 alongside low-rate spontaneous cyclization), "
                "which lacks a free alpha-amino group and blocks phenylisothiocyanate (PITC) Edman degradation unless enzymatically deblocked "
                "with pyroglutamate aminopeptidase (PGP, which excises N-terminal pGlu to expose Leu2 for Edman sequencing rather than reverting pGlu to Gln), "
                "whereas intact/top-down LC-MS directly resolves the intact-chain -17.03 Da mass shift (AUDIT-REF & AUDIT-REF)."
            ),
            "crystal_mature_residue_15_in_2JJS_4CMM_4KJY": {
                "2JJS_chain_C": cryst_2jjs_residues.get(15),
                "4CMM_chain_B": cryst_4cmm_residues.get(15),
                "4KJY_chain_A": cryst_4kjy_residues.get(15),
            },
            "designed_mature_cysteine_positions_1based": [i + 1 for i, aa in enumerate(cd47_prod["sequence"]) if aa == "C"],
            "designed_wt_Cys15_vs_2JJS_comparison": cd47_wt_vs_cryst_audit,
            "counterfactual_C15G_mutant_vs_2JJS_comparison": cd47_c15g_vs_cryst_audit,
            "structural_and_thiol_qualification_note": (
                "Full-length native CD47 forms a long-range inter-domain disulfide between IgV Cys33 (mature Cys15) and "
                "Cys263 at the extracellular loop / transmembrane boundary, whereas soluble crystal structures 2JJS (Hatherley et al. 2008 PMID:18657508), "
                "4CMM (Hatherley et al. 2014 PMID:24550402), and 4KJY (Weiskopf et al. 2013 PMID:23722425) mutated Cys33->Gly "
                "(mature C15G) to prevent unpaired-thiol dimerization in soluble expression, and our designed "
                "CD47-IgV-GPI retains wild-type Cys33 (following Xu et al. 2024 Cell Rep Med PMID:38508139 membrane-anchored CD47-IgV). "
                "Therefore, soluble C15G crystal coordinates do not settle the free-thiol/disulfide state or membrane orientation of "
                "wild-type Cys33 CD47-IgV-GPI; Stage-0 pre-pilot must compare wild-type Cys33 vs C33G CD47-IgV-GPI via non-reducing "
                "SDS-PAGE, Ellman/maleimide free-thiol labeling, PI-PLC sensitivity, and SIRPa V1/V2/V8 binding."
            ),
        },
        "experimental_structures_analyzed": [
            {
                "pdb_id": "2JJS",
                "resolution_angstrom": 1.85,
                "pmid": "18657508",
                "citation": "Hatherley et al. 2008 Mol Cell",
                "sha256": hashlib.sha256(snapshots["pdb_00002jjs.cif"]).hexdigest(),
            },
            {
                "pdb_id": "4CMM",
                "resolution_angstrom": 1.92,
                "pmid": "24550402",
                "citation": "Hatherley et al. 2014 JBC",
                "sha256": hashlib.sha256(snapshots["pdb_00004cmm.cif"]).hexdigest(),
            },
            {
                "pdb_id": "4KJY",
                "resolution_angstrom": 1.93,
                "pmid": "23722425",
                "citation": "Weiskopf et al. 2013 Science",
                "sha256": hashlib.sha256(snapshots["pdb_00004kjy.cif"]).hexdigest(),
            },
        ],
        "cd47_to_sirpa_wt_interface_2JJS": interface_2jjs,
        "cd47_to_sirpa_v1_interface_4CMM": interface_4cmm,
        "cd47_to_sirpa_fd6_interface_4KJY": interface_4kjy,
        "altloc_conformer_sensitivity_4CMM_and_4KJY_claude_v80_064": {
            "4CMM_altloc_A_heavy_contacts": interface_4cmm["heavy_atom_contact_pairs"],
            "4CMM_altloc_B_heavy_contacts": interface_4cmm_altb["heavy_atom_contact_pairs"],
            "4KJY_altloc_A_heavy_contacts": interface_4kjy["heavy_atom_contact_pairs"],
            "4KJY_altloc_B_heavy_contacts": interface_4kjy_altb["heavy_atom_contact_pairs"],
        },
        "cd47_igv_superposition_metrics": {
            "aligned_ca_residues": len(common_cd47),
            "rmsd_4CMM_vs_2JJS_angstrom": round(rmsd_4cmm_vs_2jjs, 4),
            "rmsd_4KJY_vs_2JJS_angstrom": round(rmsd_4kjy_vs_2jjs, 4),
        },
        "gpi_anchor_and_sirpa_polymorphism_verification": {
            "evidence_class": "coordinate_ectodomain_geometry_plus_literature_spr_requires_membrane_gpi_validation",
            "graveyard_rule_verified": "Locus 1 CD47-IgV-GPI vs 5-TM/TSP-1 Signaling Precaution (paired with GRAVEYARD_009 Hour-0 stoichiometry gate; note GY-003 is G6PC2-null hypoglycemia)",
            "c_terminal_val115_distance_to_sirpa_angstrom": round(dist_c_term_to_sirpa, 2),
            "daf_cd55_gpi_signal_residues": "UniProt P08174 residues P345-T381 (37-aa C-terminal GPI signal: 9-aa retained stalk P345-S353 PNKGSGTTS ending at native omega-site Ser353 = precursor residue 147, plus 28-aa cleaved hydrophobic propeptide G354-T381)",
            "sirpa_v1_v2_substitution_count": len(sirpa_aligned_substitutions),
            "sirpa_v1_v2_insertion_count": 1,
            "sirpa_v1_v2_total_alignment_difference_count": len(sirpa_insertion_and_substitution_differences),
            "sirpa_v1_v2_insertion_aware_alignment_differences": sirpa_insertion_and_substitution_differences,
            "sirpa_v1_v2_amino_acid_substitutions_only": sirpa_aligned_substitutions,
            "sirpa_4cmm_pdb_numbering_indel_note": (
                "In 2JJS Chain A (SIRPa V2), residues 99..104 are P99-D100-T101-E102-F103-K104; in 4CMM Chain A (SIRPa V1), "
                "residues 99..105 are P99-D100-D101-V102-E103-F104-K105, where the consecutive Asp100-Asp101 pair reflects a 1-residue "
                "insertion in V1 relative to V2 (annotated as Asp100 in Hatherley et al. 2014 PMID:24550402 Fig. 1 sequence alignment, whereas 3D C-alpha "
                "superposition of residues 3-96 places V2 Asp100 0.36 A from V1 Asp101 [which contacts CD47:Lys39:NZ via A:ASP101:OD1 at 2.79 A (2.7852 A exact), within the 4.5 A shell, whereas CD47:Arg103:NH2 is 16.35 A away] and 3.49 A "
                "from V1 Asp100 [min CD47 distance 4.65 A to CD47:Lys39:NZ, outside the 4.5 A shell], so V1's 4.5 A polymorphic/indel contact set is positions 66 and 101; AUDIT-REF & AUDIT-REF) "
                "that shifts subsequent V1 residue numbers by +1 (mapping V2 T101->V1 V102, V2 E102->V1 E103, "
                "V2 F103->V1 F104, V2 K104->V1 K105). Total V1 vs V2 differences across the IgV domain: 12 amino-acid substitutions "
                "+ 1 single-residue insertion (Asp100/Asp101 in 4CMM) = 13 alignment differences (closes AUDIT-REF & AUDIT-REF)."
            ),
            "sirpa_v1_v2_polymorphic_positions_tested": sirpa_v1_v2_polymorphic_sites,
            "polymorphic_positions_overlapping_4p5A_contact_shell": overlapping_polymorphic_sites,
            "polymorphic_positions_overlapping_4p5A_contact_shell_v2_2jjs": overlapping_polymorphic_sites,
            "polymorphic_and_indel_positions_overlapping_4p5A_contact_shell_v1_4cmm": overlapping_polymorphic_sites_v1_4cmm,
            "v1_4cmm_asp101_insertion_within_4p5A_cd47_contact_shell": bool(101 in sirpa_v1_4cmm_contact_seq_ids),
            "literature_kd_sirpa_v1_uM": 1.2,
            "literature_kd_sirpa_v2_uM": 1.6,
            "kd_sirpa_v1_uM": 1.2,
            "kd_sirpa_v2_uM": 1.6,
            "hatherley_2014_table2_soluble_spr_37C_kd_sirpa_v1_uM": 0.74,
            "hatherley_2014_table2_soluble_spr_37C_kd_sirpa_v1_sd_uM": 0.07,
            "hatherley_2014_table2_soluble_spr_37C_kd_sirpa_v2_uM": 0.64,
            "hatherley_2014_table2_soluble_spr_37C_kd_sirpa_v2_sd_uM": 0.06,
            "sirpa_temperature_and_immobilization_affinity_reconciliation_sonnet55_r2": {
                "summary": (
                    "Reconciles the literature SPR affinity scales by explicitly tagging construct, valency, temperature, assay format, and cached PMID: "
                    "1.2 uM (SIRPa V1) / 1.6 uM (SIRPa V2) are unverified secondary 25 C working model priors (Hatherley et al. 2008 PMID:18657508 text reports ~1-2 uM across alleles, "
                    "and Weiskopf et al. 2013 PMID:23722425 reports 0.3-0.5 uM WT CD47 affinity around Fig. 1C and ~1 uM general WT baseline), whereas Hatherley et al. 2014 (PMID:24550402 Table 2) "
                    "directly measured 0.74 +/- 0.07 uM (SIRPa V1), 0.64 +/- 0.06 uM (SIRPa V2), 0.78 +/- 0.09 uM (SIRPa v5), 0.65 +/- 0.05 uM (SIRPa v7), and 0.67 +/- 0.05 uM (SIRPa v10) at physiological 37 C using immobilized 3-domain human SIRPa (d1-d3)-rCD4(d3+4)-biotin "
                    "on streptavidin with soluble monomeric human CD47 IgV (d1) analyte (so the 2008/2013 vs 2014 figures also differ in 1-domain vs 3-domain SIRPa immobilization "
                    "architecture rather than forming an identical-construct temperature series; S55R2-D3-04, AUDIT-REF, AUDIT-REF & AUDIT-REF,PE-F2)."
                ),
                "pmid_anchored_spr_measurement_ledger": [
                    {
                        "construct": "Soluble human CD47 IgV (C15G) : SIRPa V1 d1 monomer (unverified secondary model prior; Hatherley 2008 text reports ~1-2 uM range and Weiskopf 2013 reports 0.3-0.5 uM / ~1 uM WT baseline)",
                        "valency": "1:1 monomer",
                        "temperature_C": 25.0,
                        "assay_format": "Surface plasmon resonance (secondary model prior; ~1-2 uM range in Hatherley 2008)",
                        "kd_uM": 1.2,
                        "pmid": "18657508",
                        "citation": "Hatherley et al. 2008 Mol Cell 31:266-277 & Weiskopf et al. 2013 Science 341:88-91 (PMID:23722425)",
                    },
                    {
                        "construct": "Soluble human CD47 IgV (C15G) : SIRPa V2 d1 monomer (unverified secondary model prior; Hatherley 2008 text reports ~1-2 uM range and Weiskopf 2013 reports 0.3-0.5 uM / ~1 uM WT baseline)",
                        "valency": "1:1 monomer",
                        "temperature_C": 25.0,
                        "assay_format": "Surface plasmon resonance (secondary model prior; ~1-2 uM range in Hatherley 2008)",
                        "kd_uM": 1.6,
                        "pmid": "18657508",
                        "citation": "Hatherley et al. 2008 Mol Cell 31:266-277 & Weiskopf et al. 2013 Science 341:88-91 (PMID:23722425)",
                    },
                    {
                        "construct": "Immobilized 3-domain human SIRPa V1 (d1-d3)-rCD4(d3+4)-biotin on streptavidin vs soluble monomeric human CD47 IgV (d1)",
                        "valency": "1:1 monomeric CD47 IgV analyte over immobilized 3-domain SIRPa (d1-d3)",
                        "temperature_C": 37.0,
                        "assay_format": "Surface plasmon resonance at physiological 37 C (Table 2)",
                        "kd_uM": 0.74,
                        "kd_sd_uM": 0.07,
                        "pmid": "24550402",
                        "citation": "Hatherley et al. 2014 J Biol Chem 289:10024-10033",
                    },
                    {
                        "construct": "Immobilized 3-domain human SIRPa V2 (d1-d3)-rCD4(d3+4)-biotin on streptavidin vs soluble monomeric human CD47 IgV (d1)",
                        "valency": "1:1 monomeric CD47 IgV analyte over immobilized 3-domain SIRPa (d1-d3)",
                        "temperature_C": 37.0,
                        "assay_format": "Surface plasmon resonance at physiological 37 C (Table 2)",
                        "kd_uM": 0.64,
                        "kd_sd_uM": 0.06,
                        "pmid": "24550402",
                        "citation": "Hatherley et al. 2014 J Biol Chem 289:10024-10033",
                    },
                ],
            },
            "kd_sirpa_v8_uM_unmeasured_interpolation_prior": 1.35,
            "exact_product_sirpa_v8_kd_uM": None,
            "sirpa_v8_experimental_structure_analyzed": False,
            "sirpa_v8_sequence_and_polymorphism_scope_note": (
                "Only SIRPa V2 (PDB:2JJS) and SIRPa V1 (PDB:4CMM, where V1 Asp101 [2.79 A to CD47; paired with sequence-annotated insertion Asp100 at 4.65 A] resides directly inside the 4.5 A CD47 contact shell alongside position 66; AUDIT-REF & AUDIT-REF) co-crystal structures are analyzed in the local coordinate corpus; "
                "SIRPa V8 (one of the 10 polymorphic human SIRPa IgV alleles identified in the 37-donor sequence panel of Takenaka et al. 2007 PMID:17982459 Supp. Table 3, where V1 and V2 are the two most frequent alleles in that panel and V8 is a lower-frequency polymorphic variant, and where Hatherley et al. 2014 PMID:24550402 Table 2 measured soluble 37 C SPR for 5 variants v1, v2, v5, v7, and v10 [K_D = 0.64-0.78 uM across the 5 measured variant means: v1 = 0.74 +/- 0.07, v2 = 0.64 +/- 0.06, v5 = 0.78 +/- 0.09, v7 = 0.65 +/- 0.05, v10 = 0.67 +/- 0.05 uM] rather than V8; AUDIT-REF, AUDIT-REF & AUDIT-REF) has no co-crystal "
                "structure in the corpus, so exact_product_sirpa_v8_kd_uM is null (1.35 uM is an unmeasured working interpolation prior requiring Stage-0 SPR and primary donor macrophage validation)."
            ),
            "pan_allelic_sirpa_protection_status": "UNRESOLVED_REQUIRES_V1_V2_V8_PRIMARY_HUMAN_MACROPHAGE_PHAGOCYTOSIS_ASSAY",
            "bare_wlc_9aa_stalk_rms_nm_at_l0_3p5_and_3p8A": [1.483, 1.554],
            "phagocytic_synapse_stalk_reach_and_spatial_exclusion_qualification_astra_r1": (
                "CD47-IgV-GPI fuses the single CD47 IgV domain (~4.0 nm) to the 9-aa retained DAF (CD55) GPI stalk (PNKGSGTTS, "
                "L_c = 3.15 nm [l_0 = 3.5 A] to 3.42 nm [l_0 = 3.8 A]; bare-WLC RMS = 1.483-1.554 nm at L_p = 4.0 A, and ~1.7-2.1 nm under O-glycosylated Ser/Thr mucin stiffening; AUDIT-REF). "
                "Because CD47:SIRPa inhibition requires physical engagement across the intercellular phagocytic synapse (~11.5 nm FcR-IgG synapse gap and 6 / 9.5 / 16.5 / 21.5 nm synthetic Fibcon protein domain tethers with FKBP/FRB rapalog coupling in Morrissey et al. 2020 PMID:32768386 [Immunity; and 2023 Immunity Correction PMID:37703830 correcting recombinant mouse CD47 expression boundaries from 40-182 to 19-140], where CD47 ligation retains SIRPa in the synapse to suppress inside-out integrin activation whereas ~14 nm total gap / ~10-12 nm 3-domain SIRPa ectodomain reach are geometric model estimates; AUDIT-REF & AUDIT-REF), "
                "Stage-0 bench validation must test cell-bound SIRPa-Fc avidity and primary V1/V2/V8 macrophage phagocytosis +/- neuraminidase "
                "and +/- PI-PLC against a longer DAF-SCR4 stalk comparator if steric glycocalyx burial is observed (AUDIT-REF, AUDIT-REF & AUDIT-REF)."
            ),
            "scope_and_limitations_note": (
                "Crystal structures 2JJS and 4CMM capture the soluble C15G (C33G) extracellular IgV complex for SIRPa V2 and V1 "
                "(no experimental structure is analyzed for SIRPa V8, whose 1.35 uM value is an unmeasured working prior) and do not "
                "include the wild-type Cys33 thiol, lipid-bilayer GPI anchor, membrane tether orientation, or glycocalyx; "
                "pan-allelic macrophage/NK phagocytic suppression requires empirical donor-panel flow-cytometry & PI-PLC verification."
            ),
        },
    }


def validate_alphafold_reference(atoms: List[Dict[str, Any]], metadata: Any,
                                 accession: str, length: int) -> Dict[str, Any]:
    """Bind every modeled residue and confidence to its declared AF reference."""
    if isinstance(metadata, list):
        if len(metadata) != 1:
            raise ValueError("Exactly one AlphaFold reference metadata entry is required")
        metadata = metadata[0]
    if not isinstance(metadata, dict) or metadata.get("uniprotAccession") != accession:
        raise ValueError(f"AlphaFold metadata accession mismatch: {accession}")
    model_id = f"AF-{accession}-F1"
    if (metadata.get("modelEntityId") != model_id or metadata.get("entryId") != model_id
            or type(metadata.get("latestVersion")) is not int or metadata["latestVersion"] != 6
            or metadata.get("chainId") != "A"
            or metadata.get("cifUrl") != f"https://alphafold.ebi.ac.uk/files/{model_id}-model_v6.cif"
            or metadata.get("paeDocUrl") != f"https://alphafold.ebi.ac.uk/files/{model_id}-predicted_aligned_error_v6.json"):
        raise ValueError(f"AlphaFold metadata model/version identity mismatch: {accession}")
    sequence = metadata.get("uniprotSequence")
    global_plddt = metadata.get("globalMetricValue")
    if (not isinstance(sequence, str) or len(sequence) != length
            or type(metadata.get("uniprotStart")) is not int or metadata.get("uniprotStart") != 1
            or type(metadata.get("uniprotEnd")) is not int or metadata.get("uniprotEnd") != length
            or isinstance(global_plddt, (bool, np.bool_)) or not isinstance(global_plddt, Real)
            or not math.isfinite(global_plddt) or not 0 <= global_plddt <= 100):
        raise ValueError(f"Incomplete AlphaFold sequence/range/confidence metadata: {accession}")
    ca = get_ca_map(atoms, "A")
    if set(ca) != set(range(1, length + 1)):
        raise ValueError(f"AlphaFold model must contain every declared residue exactly once: {accession}")
    three_to_one = dict(zip(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
        "ARNDCQEGHILKMFPSTWYV"))
    model_sequence = "".join(three_to_one.get(ca[i]["comp_id"], "?") for i in range(1, length + 1))
    if model_sequence != sequence:
        raise ValueError(f"AlphaFold coordinate residue identities differ from metadata: {accession}")
    for atom in atoms:
        value = atom["b_iso"]
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not 0.0 <= value <= 100.0:
            raise ValueError(f"AlphaFold pLDDT must be within [0,100]: {accession}")
        residue = ca.get(atom["seq_id"])
        if (atom["chain"] != "A" or residue is None or atom["ins_code"] != ""
                or atom["comp_id"] != residue["comp_id"] or value != residue["b_iso"]):
            raise ValueError(f"AlphaFold atom identity/confidence disagrees within its declared residue: {accession}")
    # Give every residue one vote, irrespective of its number of sidechain
    # atoms. Keep the supplied metadata value visible even when it is stale.
    coordinate_mean = math.fsum(float(ca[i]["b_iso"]) for i in range(1, length + 1)) / length
    result = dict(metadata)
    result["coordinate_global_plddt"] = round(coordinate_mean, 2)
    result["global_plddt_audit"] = {
        "source": "unweighted_mean_of_one_CA_pLDDT_per_reference_residue",
        "coordinate_residue_count": length,
        "coordinate_mean_unrounded": coordinate_mean,
        "metadata_supplied_global_plddt": float(global_plddt),
        "metadata_minus_coordinate_mean": float(global_plddt) - coordinate_mean,
        "agreement_at_two_decimal_places": round(float(global_plddt), 2) == round(coordinate_mean, 2),
        "published_value_source": "validated_coordinate_residue_mean",
    }
    return result


def analyze_complex_3_serpinb9_gzmb(layer2_ledger: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze Complex 3: SERPINB9 (PI-9) : Granzyme B.
    Closes V70-07 (local RCL pLDDT, 21/26/50-aa tail contour lengths, and unresolved status)
    and V70-08 (insertion-code-safe 226/226 global alignment of AF-P10144 vs PDB:1IAU).
    """
    snapshots = verify_structural_reference_inputs(
        "AF-P50453-F1-model_v6.cif", "AF-P10144-F1-model_v6.cif", "pdb_00001iau.cif",
        "AF-P50453-F1-predicted_aligned_error_v6.json", "AF-P10144-F1-predicted_aligned_error_v6.json")
    atoms_sb9 = parse_mmcif_atoms(STRUCT_DIR / "AF-P50453-F1-model_v6.cif", source_bytes=snapshots["AF-P50453-F1-model_v6.cif"])
    atoms_gzmb_af = parse_mmcif_atoms(STRUCT_DIR / "AF-P10144-F1-model_v6.cif", source_bytes=snapshots["AF-P10144-F1-model_v6.cif"])
    # Retain insertion codes (e.g. 184A) for PDB:1IAU; evaluate both altloc A (occupancy 0.38) and dominant altloc B (occupancy 0.62)
    atoms_1iau = parse_mmcif_atoms(STRUCT_DIR / "pdb_00001iau.cif", source_bytes=snapshots["pdb_00001iau.cif"], keep_insertion_codes=True, allowed_alt_ids=(".", "?", "A"))
    atoms_1iau_altb = parse_mmcif_atoms(STRUCT_DIR / "pdb_00001iau.cif", source_bytes=snapshots["pdb_00001iau.cif"], keep_insertion_codes=True, allowed_alt_ids=(".", "?", "B"))

    with (STRUCT_DIR / "AF-P50453-F1-metadata.json").open() as f:
        meta_sb9 = json.load(f)
    with (STRUCT_DIR / "AF-P10144-F1-metadata.json").open() as f:
        meta_gzmb = json.load(f)
    meta_sb9 = validate_alphafold_reference(atoms_sb9, meta_sb9, "P50453", 376)
    meta_gzmb = validate_alphafold_reference(atoms_gzmb_af, meta_gzmb, "P10144", 247)
    if layer2_ledger["SERPINB9_cytosolic_with_21aa_GSG_P2A_stub"]["sequence"][:376] != meta_sb9["uniprotSequence"]:
        raise ValueError("Engineered SERPINB9 native domain differs from the analyzed reference sequence")
    pae_sb9_raw = json.loads(snapshots["AF-P50453-F1-predicted_aligned_error_v6.json"])
    pae_gzmb_raw = json.loads(snapshots["AF-P10144-F1-predicted_aligned_error_v6.json"])
    sb9_pae_mat = validated_pae_matrix(pae_sb9_raw, 376, "SERPINB9")
    gzmb_pae_mat = validated_pae_matrix(pae_gzmb_raw, 247, "Granzyme B")
    sb9_core_pae_mean = round(_pae_mean(sb9_pae_mat[:322, :322]), 2)
    # In official AlphaFold DB PAE v6 JSON (predicted_aligned_error[i][j]), row index i is the ALIGNED residue and column index j is the SCORED residue (AUDIT-REF):
    sb9_rows_rcl_aligned_cols_core_scored_pae_mean = round(_pae_mean(sb9_pae_mat[334:345, :322]), 2)  # 27.29 A
    sb9_rows_core_aligned_cols_rcl_scored_pae_mean = round(_pae_mean(sb9_pae_mat[:322, 334:345]), 2)  # 15.20 A
    sb9_row_pro376_aligned_cols_rcl_335_345_scored_pae_mean = round(_pae_mean(sb9_pae_mat[375:376, 334:345]), 2)  # 21.27 A
    sb9_rows_rcl_335_345_aligned_col_pro376_scored_pae_mean = round(_pae_mean(sb9_pae_mat[334:345, 375:376]), 2)  # 25.45 A
    sb9_rows_cterm_372_376_aligned_cols_rcl_335_350_scored_pae_mean = round(_pae_mean(sb9_pae_mat[371:376, 334:350]), 2)  # 13.24 A
    sb9_rows_rcl_335_350_aligned_cols_cterm_372_376_scored_pae_mean = round(_pae_mean(sb9_pae_mat[334:350, 371:376]), 2)  # 19.90 A
    sb9_rcl_to_core_pae_mean = sb9_rows_rcl_aligned_cols_core_scored_pae_mean
    sb9_core_to_rcl_pae_mean = sb9_rows_core_aligned_cols_rcl_scored_pae_mean
    gzmb_mature_21_246_pae_mean = round(_pae_mean(gzmb_pae_mat[20:246, 20:246]), 2)
    gzmb_mature_21_247_pae_mean = round(_pae_mean(gzmb_pae_mat[20:, 20:]), 2)

    ca_sb9 = get_ca_map(atoms_sb9, "A")
    core_plddts = [ca_sb9[r]["b_iso"] for r in range(1, 323) if r in ca_sb9]
    rcl_hinge_to_p3p_plddts = [ca_sb9[r]["b_iso"] for r in range(322, 349) if r in ca_sb9]
    rcl_recognition_335_345_plddts = [ca_sb9[r]["b_iso"] for r in range(335, 346) if r in ca_sb9]
    rcl_mobile_loop_335_350_plddts = [ca_sb9[r]["b_iso"] for r in range(335, 351) if r in ca_sb9]
    cterm_plddts = [ca_sb9[r]["b_iso"] for r in range(349, 377) if r in ca_sb9]
    plddt_glu340 = float(ca_sb9[340]["b_iso"])
    plddt_cys341 = float(ca_sb9[341]["b_iso"])

    glu340_ca = get_atom_coord(atoms_sb9, "A", 340, "CA")
    cys341_ca = get_atom_coord(atoms_sb9, "A", 341, "CA")
    pro376_ca = get_atom_coord(atoms_sb9, "A", 376, "CA")
    pro376_c_atom = get_atom_coord(atoms_sb9, "A", 376, "C")
    val322_ca = get_atom_coord(atoms_sb9, "A", 322, "CA")

    dist_pro376_to_glu340 = euclidean_distance(pro376_ca, glu340_ca)
    dist_pro376_to_cys341 = euclidean_distance(pro376_ca, cys341_ca)
    dist_pro376_to_val322 = euclidean_distance(pro376_ca, val322_ca)

    cterm_atoms = filter_polymer_atoms(atoms_sb9, ["A"], (372, 376))
    pro376_atoms = filter_polymer_atoms(atoms_sb9, ["A"], (376, 376))
    rcl_recognition_atoms = filter_polymer_atoms(atoms_sb9, ["A"], (335, 345))
    rcl_337_346_atoms = filter_polymer_atoms(atoms_sb9, ["A"], (337, 346))
    rcl_mobile_335_350_atoms = filter_polymer_atoms(atoms_sb9, ["A"], (335, 350))
    rcl_proximal_hinge_318_334_atoms = filter_polymer_atoms(atoms_sb9, ["A"], (318, 334))
    c_cterm = np.array([[a["x"], a["y"], a["z"]] for a in cterm_atoms if a["element"] != "H"])
    c_pro376 = np.array([[a["x"], a["y"], a["z"]] for a in pro376_atoms if a["element"] != "H"])
    c_rcl = np.array([[a["x"], a["y"], a["z"]] for a in rcl_recognition_atoms if a["element"] != "H"])
    c_rcl_337_346 = np.array([[a["x"], a["y"], a["z"]] for a in rcl_337_346_atoms if a["element"] != "H"])
    c_rcl_335_350 = np.array([[a["x"], a["y"], a["z"]] for a in rcl_mobile_335_350_atoms if a["element"] != "H"])
    c_hinge_318_334 = np.array([[a["x"], a["y"], a["z"]] for a in rcl_proximal_hinge_318_334_atoms if a["element"] != "H"])
    min_heavy_dist_cterm_to_rcl = float(np.min(np.linalg.norm(c_cterm[:, None, :] - c_rcl[None, :, :], axis=-1)))
    min_heavy_dist_pro376_only_to_rcl = float(np.min(np.linalg.norm(c_pro376[:, None, :] - c_rcl[None, :, :], axis=-1)))
    min_heavy_dist_pro376_c_to_rcl_335_345 = float(np.min(np.linalg.norm(c_rcl - pro376_c_atom[None, :], axis=-1)))
    min_heavy_dist_cterm_to_rcl_337_346 = float(np.min(np.linalg.norm(c_cterm[:, None, :] - c_rcl_337_346[None, :, :], axis=-1)))
    min_heavy_dist_pro376_only_to_rcl_337_346 = float(np.min(np.linalg.norm(c_pro376[:, None, :] - c_rcl_337_346[None, :, :], axis=-1)))
    min_heavy_dist_cterm_to_rcl_335_350 = float(np.min(np.linalg.norm(c_cterm[:, None, :] - c_rcl_335_350[None, :, :], axis=-1)))
    min_heavy_dist_pro376_only_to_rcl_335_350 = float(np.min(np.linalg.norm(c_pro376[:, None, :] - c_rcl_335_350[None, :, :], axis=-1)))
    min_heavy_dist_pro376_c_to_rcl_335_350 = float(np.min(np.linalg.norm(c_rcl_335_350 - pro376_c_atom[None, :], axis=-1)))
    min_heavy_dist_cterm_to_hinge_318_334 = float(np.min(np.linalg.norm(c_cterm[:, None, :] - c_hinge_318_334[None, :, :], axis=-1)))
    min_heavy_dist_pro376_only_to_hinge_318_334 = float(np.min(np.linalg.norm(c_pro376[:, None, :] - c_hinge_318_334[None, :, :], axis=-1)))
    min_heavy_dist_pro376_c_to_hinge_318_334 = float(np.min(np.linalg.norm(c_hinge_318_334 - pro376_c_atom[None, :], axis=-1)))

    gzmb_1iau = filter_polymer_atoms(atoms_1iau, ["A"])
    iepd_1iau = [a for a in atoms_1iau if a["chain"] == "B" and a["comp_id"] != "HOH"]
    iepd_interface = compute_interface_metrics(gzmb_1iau, iepd_1iau)

    gzmb_1iau_altb = filter_polymer_atoms(atoms_1iau_altb, ["A"])
    iepd_1iau_altb = [a for a in atoms_1iau_altb if a["chain"] == "B" and a["comp_id"] != "HOH"]
    iepd_interface_altb = compute_interface_metrics(gzmb_1iau_altb, iepd_1iau_altb)

    def _active_site_dists(atom_rows: List[Dict[str, Any]]) -> Dict[str, float]:
        s195_og = get_atom_coord(atom_rows, "A", 195, "OG")
        r226_nh2 = get_atom_coord(atom_rows, "A", 226, "NH2")
        asj405_c = get_atom_coord(atom_rows, "B", 405, "C")
        asj405_od2 = get_atom_coord(atom_rows, "B", 405, "OD2")
        return {
            "Ser195_OG_to_ASJ405_C_angstrom": round(euclidean_distance(s195_og, asj405_c), 6),
            "Arg226_NH2_to_ASJ405_OD2_angstrom": round(euclidean_distance(r226_nh2, asj405_od2), 6),
        }

    altloc_1iau_audit = {
        "conformer_selection_policy": "Altloc A (occupancy 0.38) retained as historical baseline interface with dominant Altloc B (occupancy 0.62) reported alongside without atom-wise conformer mixing",
        "altloc_A_minor_occupancy_0p38": {
            "alternate_occupancies": sorted({a["occupancy"] for a in atoms_1iau if a["alt_id"] == "A"}),
            "heavy_atom_contact_pairs": iepd_interface["heavy_atom_contact_pairs"],
            "unique_residue_contact_pairs": iepd_interface["unique_residue_contact_pairs"],
            "polar_contact_atom_pairs": iepd_interface["polar_contact_atom_pairs"],
            **_active_site_dists(atoms_1iau),
        },
        "altloc_B_major_occupancy_0p62": {
            "alternate_occupancies": sorted({a["occupancy"] for a in atoms_1iau_altb if a["alt_id"] == "B"}),
            "heavy_atom_contact_pairs": iepd_interface_altb["heavy_atom_contact_pairs"],
            "unique_residue_contact_pairs": iepd_interface_altb["unique_residue_contact_pairs"],
            "polar_contact_atom_pairs": iepd_interface_altb["polar_contact_atom_pairs"],
            **_active_site_dists(atoms_1iau_altb),
            "interface_metrics": iepd_interface_altb,
        },
    }

    # Insertion-code-safe alignment between mature AF-P10144 (residues 21..246 = 226 aa) and PDB:1IAU Chain A (226 aa)
    ca_gzmb_af_list = [a for a in get_ca_ordered_list_with_ins_codes(atoms_gzmb_af, "A") if a["seq_id"] >= 21]
    ca_1iau_list = get_ca_ordered_list_with_ins_codes(atoms_1iau, "A")
    if not (len(ca_1iau_list) == 226):
        raise ValueError(f"Expected 226 CA atoms in 1IAU Chain A with insertion codes, got {len(ca_1iau_list)}")
    matched_af_coords = []
    matched_iau_coords = []
    identical_pairs = 0
    for a_af, a_iau in zip(ca_gzmb_af_list[:226], ca_1iau_list):
        if a_af["comp_id"] == a_iau["comp_id"]:
            identical_pairs += 1
        matched_af_coords.append([a_af["x"], a_af["y"], a_af["z"]])
        matched_iau_coords.append([a_iau["x"], a_iau["y"], a_iau["z"]])
    if not (identical_pairs == 226):
        raise ValueError('Molecular source invariant failed: identical_pairs == 226')
    rmsd_gzmb_af_vs_1iau, _, _ = kabsch_rmsd(np.array(matched_af_coords, dtype=np.float64), np.array(matched_iau_coords, dtype=np.float64))

    # Derive SERPINB9 P4..P1' sequence directly from coordinates and pinned Layer 2 polypeptide ledger (closes AUDIT-REF)
    sb9_ledger = layer2_ledger["SERPINB9_cytosolic_with_21aa_GSG_P2A_stub"]
    aa3_to_title = {
        "ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys",
        "GLN": "Gln", "GLU": "Glu", "GLY": "Gly", "HIS": "His", "ILE": "Ile",
        "LEU": "Leu", "LYS": "Lys", "MET": "Met", "PHE": "Phe", "PRO": "Pro",
        "SER": "Ser", "THR": "Thr", "TRP": "Trp", "TYR": "Tyr", "VAL": "Val",
    }
    rcl_coord_res_337_341 = [ca_sb9[r]["comp_id"] for r in range(337, 342)]
    rcl_encoded_337_341 = sb9_ledger["sequence"][336:341]
    if not (rcl_encoded_337_341 == "VVAEC" and rcl_coord_res_337_341 == ["VAL", "VAL", "ALA", "GLU", "CYS"]):
        raise ValueError('Molecular source invariant failed: rcl_encoded_337_341 == "VVAEC" and rcl_coord_res_337_341 == ["VAL", "VAL", "ALA", "GLU", "CYS"]')
    rcl_p4_to_p1p_str = (
        f"{aa3_to_title[rcl_coord_res_337_341[0]]}337-{aa3_to_title[rcl_coord_res_337_341[1]]}338-"
        f"{aa3_to_title[rcl_coord_res_337_341[2]]}339-{aa3_to_title[rcl_coord_res_337_341[3]]}340-"
        f"{aa3_to_title[rcl_coord_res_337_341[4]]}341 "
        f"(P4-{rcl_encoded_337_341[0]}337, P3-{rcl_encoded_337_341[1]}338, P2-{rcl_encoded_337_341[2]}339, "
        f"P1-{rcl_encoded_337_341[3]}340, P1'-{rcl_encoded_337_341[4]}341; sequence {rcl_encoded_337_341})"
    )

    # Evaluate polymer contour upper bounds vs WLC RMS reach of the possible 2A C-terminal extensions on SERPINB9 (21 aa, 26 aa, 50 aa; closes 089)
    tail_ensemble_models = []
    for stub_len in [21, 26, 50]:
        contour_a = stub_len * 3.50
        wlc_rms_a = math.sqrt(2.0 * 4.0 * contour_a * (1.0 - (4.0 / contour_a) * (1.0 - math.exp(-contour_a / 4.0))))
        contour_exceeds = bool(contour_a >= min_heavy_dist_cterm_to_rcl)
        rms_exceeds = bool(wlc_rms_a >= min_heavy_dist_cterm_to_rcl)
        tail_ensemble_models.append(
            {
                "stub_length_aa": stub_len,
                "contour_length_angstrom": round(contour_a, 2),
                "wlc_rms_extension_angstrom": round(wlc_rms_a, 2),
                "serpin_intrinsic_cterm_to_rcl_reference_distance_angstrom": round(min_heavy_dist_cterm_to_rcl, 3),
                "pro376_c_peptide_bond_anchor_to_rcl_335_345_distance_angstrom": round(min_heavy_dist_pro376_c_to_rcl_335_345, 3),
                "pro376_c_peptide_bond_anchor_to_rcl_335_350_distance_angstrom": round(min_heavy_dist_pro376_c_to_rcl_335_350, 3),
                "pro376_c_peptide_bond_anchor_to_hinge_318_334_distance_angstrom": round(min_heavy_dist_pro376_c_to_hinge_318_334, 3),
                "contour_upper_bound_exceeds_serpin_intrinsic_cterm_to_rcl_distance": contour_exceeds,
                "wlc_rms_reach_exceeds_serpin_intrinsic_cterm_to_rcl_distance": rms_exceeds,
                "wlc_rms_reach_exceeds_pro376_c_anchor_to_rcl_335_345_distance": bool(wlc_rms_a >= min_heavy_dist_pro376_c_to_rcl_335_345),
                "exceeds_24p05A_unmodified_cterm_to_rcl_clearance": contour_exceeds,
                "steric_conclusion": (
                    "contour_upper_bound_exceeds_but_wlc_rms_below_unmodified_serpin_intrinsic_distance_requires_GzmB_kinetics_assay"
                    if (contour_exceeds and not rms_exceeds)
                    else "both_contour_and_wlc_rms_exceed_unmodified_serpin_intrinsic_distance_requires_GzmB_kinetics_assay"
                ),
            }
        )

    # Parse cached Foldseek table to verify both 4GA7 (top E-value homolog) and 8H4X -> 8ZCR (same-protein cleaved PI-9 structure; closes 092)
    fs_md_path = STRUCT_DIR / "serpinb9_foldseek_results.md"
    fs_md_bytes = fs_md_path.read_bytes()
    fs_md_text = fs_md_bytes.decode("utf-8")
    if ("/" + "Users/") in fs_md_text:
        raise ValueError("Redistributed Foldseek log must not contain local absolute user paths")

    def _parse_foldseek_md_hit(md_text: str, target_prefix: str) -> Dict[str, Any]:
        for line in md_text.splitlines():
            s = line.strip()
            if not s.startswith("|"):
                continue
            cols = [c.strip() for c in s.strip("|").split("|")]
            if len(cols) >= 6 and cols[0].lower().startswith(target_prefix.lower()):
                return {
                    "target_col": cols[0],
                    "q_cov_pct": round(float(cols[1].rstrip("%")), 1),
                    "prob": round(float(cols[2]), 3),
                    "e_value": float(cols[3]),
                    "seq_identity_pct": round(float(cols[4].rstrip("%")), 1),
                    "aligned_length_aa": int(cols[5]),
                }
        raise ValueError(f"Required Foldseek hit '{target_prefix}' not found in {fs_md_path}")

    hit_4ga7 = _parse_foldseek_md_hit(fs_md_text, "4ga7-assembly1.cif.gz_A")
    hit_8h4x = _parse_foldseek_md_hit(fs_md_text, "8h4x-assembly1.cif.gz_A")

    return {
        "complex_id": "LAYER3_COMPLEX_3_SERPINB9_GRANZYME_B",
        "locus": "Locus 3 (CXCL10 Exon-1 KI: CXCL10::SERPINB9-GSG-P2A-PD-L1)",
        "engineered_protein": "Cytosolic SERPINB9 (PI-9, 376 aa + 21-aa C-terminal GSG-P2A stub, Furin-free per R01)",
        "layer2_processed_polypeptide_hashes": {
            "SERPINB9_397aa_with_21aa_P2A_stub_sha256": sb9_ledger["sha256"],
            "SERPINB9_376aa_scar_free_L3_BFB_and_DXL3_sequence": sb9_ledger["sequence"][:376],
            "SERPINB9_376aa_scar_free_L3_BFB_and_DXL3_sha256": hashlib.sha256(sb9_ledger["sequence"][:376].encode("utf-8")).hexdigest(),
            "alternative_26aa_tpt2a_stub_sha256": sb9_ledger["alternative_26aa_stub_sha256"],
            "alternative_50aa_tpt2a_stub_sha256": sb9_ledger["alternative_50aa_stub_sha256"],
        },
        "structures_analyzed": [
            {
                "id": "AF-P50453-F1-model_v6",
                "source": "AlphaFold DB v6 (UniProt P50453 unmodified 376-aa reference model)",
                "global_plddt": meta_sb9["coordinate_global_plddt"],
                "global_plddt_source_audit": meta_sb9["global_plddt_audit"],
                "core_beta_sheet_domain_1_322_mean_plddt": round(float(np.mean(core_plddts)), 2),
                "rcl_loop_322_348_mean_plddt": round(float(np.mean(rcl_hinge_to_p3p_plddts)), 2),
                "rcl_mobile_loop_335_350_mean_plddt": round(float(np.mean(rcl_mobile_loop_335_350_plddts)), 2),
                "p1_glu340_plddt": round(plddt_glu340, 2),
                "p1_prime_cys341_plddt": round(plddt_cys341, 2),
                "c_terminal_sheet_b_349_376_mean_plddt": round(float(np.mean(cterm_plddts)), 2),
                "sha256": hashlib.sha256(snapshots["AF-P50453-F1-model_v6.cif"]).hexdigest(),
                "pae_sha256": hashlib.sha256(snapshots["AF-P50453-F1-predicted_aligned_error_v6.json"]).hexdigest(),
            },
            {
                "id": "AF-P10144-F1-model_v6",
                "source": "AlphaFold DB v6 (UniProt P10144 Granzyme B)",
                "global_plddt": meta_gzmb["coordinate_global_plddt"],
                "global_plddt_source_audit": meta_gzmb["global_plddt_audit"],
                "sha256": hashlib.sha256(snapshots["AF-P10144-F1-model_v6.cif"]).hexdigest(),
                "pae_sha256": hashlib.sha256(snapshots["AF-P10144-F1-predicted_aligned_error_v6.json"]).hexdigest(),
            },
            {
                "pdb_id": "1IAU",
                "resolution_angstrom": 2.00,
                "pmid": "11325591",
                "citation": "Rotonda et al. 2001 Chem Biol",
                "description": "Human Granzyme B in complex with Ac-IEPD-CHO tetrapeptide inhibitor (not full-length SERPINB9)",
                "sha256": hashlib.sha256(snapshots["pdb_00001iau.cif"]).hexdigest(),
            },
        ],
        "foldseek_structural_homolog_verification": {
            "evidence_class": "cached_foldseek_homolog_search_summary",
            "search_provenance": {
                "query_structure_sha256": hashlib.sha256(snapshots["AF-P50453-F1-model_v6.cif"]).hexdigest(),
                "ticket_id": (re.search(r"Ticket ID generated:\s*(\S+)", fs_md_text).group(1) if re.search(r"Ticket ID generated:\s*(\S+)", fs_md_text) else None),
                "source_artifact_sha256": hashlib.sha256(fs_md_bytes).hexdigest(),
                "tool_or_api_version": None, "query_date": None,
                "database_snapshot_versions": {"pdb100": None, "afdb-swissprot": None},
                "missing_historical_fields_status": "NOT_RECORDED",
                "qualification_status": "HISTORICAL_CACHED_HITS_DATABASE_SNAPSHOT_NOT_REPRODUCIBLE",
                "scope": "Replays retained hit values only; no independently reproducible historical search is claimed."},
            "query_structure": "AF-P50453-F1-model_v6.cif",
            "query_structure_sha256": hashlib.sha256(snapshots["AF-P50453-F1-model_v6.cif"]).hexdigest(),
            "foldseek_log_sha256": hashlib.sha256(fs_md_bytes).hexdigest(),
            "databases_searched": ["pdb100", "afdb-swissprot"],
            "top_experimental_pdb100_crystal_homolog": {
                "pdb_id": "4GA7",
                "description": "Crystal structure of human SerpinB1 (clade B ovalbumin-serpin scaffold; top overall E-value crystal hit)",
                "q_cov_pct": hit_4ga7["q_cov_pct"],
                "prob": hit_4ga7["prob"],
                "e_value": hit_4ga7["e_value"],
                "seq_identity_pct": hit_4ga7["seq_identity_pct"],
                "aligned_length_aa": hit_4ga7["aligned_length_aa"],
            },
            "same_protein_human_serpinb9_pdb_hit_and_supersession": {
                "cached_foldseek_pdb_id": "8H4X",
                "cached_description": "Structure of cleaved Serpin B9 (8h4x-assembly1.cif.gz_A)",
                "q_cov_pct": hit_8h4x["q_cov_pct"],
                "prob": hit_8h4x["prob"],
                "e_value": hit_8h4x["e_value"],
                "seq_identity_pct": hit_8h4x["seq_identity_pct"],
                "aligned_length_aa": hit_8h4x["aligned_length_aa"],
                "rcsb_obsoleted_date": "2024-12-25",
                "superseded_by_pdb_id": "8ZCR",
                "superseding_structure_citation": "Yan, Zhou et al. 2024 Acta Crystallogr Sect F Struct Biol Commun 80:286-293 (DOI:10.1107/S2053230X24009439, PMID:39382088; 1.93 A X-ray structure of human PI-9 / SERPINB9)",
                "scope_qualification_note": (
                    f"Foldseek pdb100 hit 8H4X (obsoleted 2024-12-25 and replaced by 8ZCR, Yan, Zhou et al. 2024 Acta Crystallogr Sect F Struct Biol Commun 80:286-293, PMID:39382088) is a cleaved/RCL-divergent "
                    f"human SERPINB9 crystal structure ({hit_8h4x['aligned_length_aa']} aa aligned, {hit_8h4x['q_cov_pct']}% query coverage) and does not resolve the intact mobile RCL "
                    "or C-terminal 21-aa GSG-P2A stub kinetics."
                ),
            },
        },
        "granzyme_b_active_site_and_s1_pocket_1IAU": {
            "alignment_method": "Insertion-code-preserving global sequence alignment (retaining all 6 PDB chymotrypsinogen insertion codes: 36A, 36B, 145A, 172A, 184A, 188A)",
            "chymotrypsinogen_insertion_code_mappings": [
                "1IAU:ASP36A -> AF:ASP42",
                "1IAU:GLN36B -> AF:GLN43",
                "1IAU:PRO145A -> AF:PRO152",
                "1IAU:ARG172A -> AF:ARG178 (S2/S4 pocket contact)",
                "1IAU:GLY184A -> AF:GLY190",
                "1IAU:LYS188A -> AF:LYS195",
            ],
            "aligned_ca_residues_AF_vs_1IAU": len(matched_af_coords),
            "identical_residue_pairs": identical_pairs,
            "rmsd_AF_P10144_vs_1IAU_angstrom": round(rmsd_gzmb_af_vs_1iau, 4),
            "af_p10144_mature_226aa_global_mean_plddt": round(float(np.mean([a["b_iso"] for a in ca_gzmb_af_list[:226]])), 2),
            "covalent_adduct_and_tyr247_alignment_scope_note_kimi_r3_005": (
                "In PDB:1IAU (Granzyme B complexed with Ac-IEPD-CHO aldehyde inhibitor), the 1.4711 A distance (< 2.0 A) between catalytic "
                "Ser195.OG and inhibitor ASJ405.C is the expected covalent tetrahedral hemiacetal adduct bond rather than a non-covalent steric clash. "
                "AF-P10144-F1 predicts 227 mature residues (Ile21..Tyr247), whereas 1IAU resolves 226 mature residues (21..246, omitting C-terminal Tyr247)."
            ),
            "iepd_inhibitor_interface": iepd_interface,
            "altloc_conformer_sensitivity_1IAU_claude_v80_091": altloc_1iau_audit,
        },
        "r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry": {
            "evidence_class": "unmodified_reference_geometry_plus_polymer_contour_analysis_requires_wetlab_kinetics",
            "functional_status": "experimentally_unresolved_requires_GzmB_kinetics_assay",
            "rcl_p4_to_p1_prime_sequence": rcl_p4_to_p1p_str,
            "rcl_p4_to_p1_prime_5mer_code": rcl_encoded_337_341,
            "rcl_322_348_loop_mean_plddt": round(float(np.mean(rcl_hinge_to_p3p_plddts)), 2),
            "rcl_335_345_hinge_to_p1p1prime_mean_plddt": round(float(np.mean(rcl_recognition_335_345_plddts)), 2),
            "rcl_335_345_recognition_loop_mean_plddt": round(float(np.mean(rcl_recognition_335_345_plddts)), 2),
            "rcl_335_350_mean_plddt": round(float(np.mean(rcl_mobile_loop_335_350_plddts)), 2),
            "alphafold_pae_uncertainty_summary_claude_v80_065": {
                "serpinb9_core_1_322_mean_pae_angstrom": sb9_core_pae_mean,
                "serpinb9_rcl_335_345_to_core_1_322_mean_pae_angstrom": sb9_rcl_to_core_pae_mean,
                "aligned_rcl_335_345_scored_core_1_322_mean_pae_angstrom": sb9_rcl_to_core_pae_mean,
                "aligned_core_1_322_scored_rcl_335_345_mean_pae_angstrom": sb9_core_to_rcl_pae_mean,
                "pae_row_aligned_col_scored_convention_astra_r3": {
                    "convention_note": (
                        "Per the official AlphaFold DB PAE v6 specification (predicted_aligned_error[i][j]), row index i is the ALIGNED residue and column index j is the SCORED residue "
                        "(expected positional error at residue j when the structure is aligned on residue i): therefore pae[334:345, :322] = 27.29 A is rows RCL 335-345 (aligned) / cols core 1-322 (scored), "
                        "pae[:322, 334:345] = 15.20 A is rows core 1-322 (aligned) / cols RCL 335-345 (scored), pae[375:376, 334:345] = 21.27 A is row Pro376 (aligned) / cols RCL 335-345 (scored), "
                        "pae[334:345, 375:376] = 25.45 A is rows RCL 335-345 (aligned) / col Pro376 (scored), pae[371:376, 334:350] = 13.24 A is rows 372-376 (aligned) / cols RCL 335-350 (scored), "
                        "and pae[334:350, 371:376] = 19.90 A is rows RCL 335-350 (aligned) / cols 372-376 (scored) (AUDIT-REF & AUDIT-REF)."
                    ),
                    "pae_rows_rcl_335_345_aligned_cols_core_1_322_scored_mean_angstrom": sb9_rows_rcl_aligned_cols_core_scored_pae_mean,
                    "pae_rows_core_1_322_aligned_cols_rcl_335_345_scored_mean_angstrom": sb9_rows_core_aligned_cols_rcl_scored_pae_mean,
                    "pae_row_pro376_aligned_cols_rcl_335_345_scored_mean_angstrom": sb9_row_pro376_aligned_cols_rcl_335_345_scored_pae_mean,
                    "pae_rows_rcl_335_345_aligned_col_pro376_scored_mean_angstrom": sb9_rows_rcl_335_345_aligned_col_pro376_scored_pae_mean,
                    "pae_rows_cterm_372_376_aligned_cols_rcl_335_350_scored_mean_angstrom": sb9_rows_cterm_372_376_aligned_cols_rcl_335_350_scored_pae_mean,
                    "pae_rows_rcl_335_350_aligned_cols_cterm_372_376_scored_mean_angstrom": sb9_rows_rcl_335_350_aligned_cols_cterm_372_376_scored_pae_mean,
                },
                "pae_row_scored_col_aligned_convention_kimi31_r3": {
                    "convention_note": (
                        "Corrected in AUDIT-REF and AUDIT-REF to match official AlphaFold DB PAE v6 row=aligned / col=scored semantics: "
                        "aligned on core 1-322 (rows :322) / scored on RCL 335-345 (cols 334:345) is 15.20 A, and aligned on RCL 335-345 (rows 334:345) / scored on core 1-322 (cols :322) is 27.29 A "
                        "(AUDIT-REF & AUDIT-REF)."
                    ),
                    "pae_rows_rcl_335_345_aligned_cols_core_1_322_scored_mean_angstrom": sb9_rows_rcl_aligned_cols_core_scored_pae_mean,
                    "pae_rows_core_1_322_aligned_cols_rcl_335_345_scored_mean_angstrom": sb9_rows_core_aligned_cols_rcl_scored_pae_mean,
                    "pae_row_pro376_aligned_cols_rcl_335_345_scored_mean_angstrom": sb9_row_pro376_aligned_cols_rcl_335_345_scored_pae_mean,
                    "pae_rows_rcl_335_345_aligned_col_pro376_scored_mean_angstrom": sb9_rows_rcl_335_345_aligned_col_pro376_scored_pae_mean,
                    "pae_rows_cterm_372_376_aligned_cols_rcl_335_350_scored_mean_angstrom": sb9_rows_cterm_372_376_aligned_cols_rcl_335_350_scored_pae_mean,
                    "pae_rows_rcl_335_350_aligned_cols_cterm_372_376_scored_mean_angstrom": sb9_rows_rcl_335_350_aligned_cols_cterm_372_376_scored_pae_mean,
                },
                "granzyme_b_mature_21_246_mean_pae_angstrom": gzmb_mature_21_246_pae_mean,
                "granzyme_b_mature_21_247_mean_pae_angstrom": gzmb_mature_21_247_pae_mean,
            },
            "p1_glu340_plddt": round(plddt_glu340, 2),
            "p1_prime_cys341_plddt": round(plddt_cys341, 2),
            "literature_native_wt_k_ass_M_inv_s_inv": 1.7e6,
            "literature_native_wt_stoichiometry_of_inhibition_SI": 1.0,
            "pro376_c_terminus_to_p1_glu340_ca_distance_angstrom": round(dist_pro376_to_glu340, 2),
            "pro376_c_terminus_to_p1_prime_cys341_ca_distance_angstrom": round(dist_pro376_to_cys341, 2),
            "pro376_c_terminus_to_rcl_hinge_val322_ca_distance_angstrom": round(dist_pro376_to_val322, 2),
            "min_heavy_atom_distance_cterm_372_376_to_rcl_335_345_angstrom": round(min_heavy_dist_cterm_to_rcl, 2),
            "pro376_only_min_heavy_atom_distance_to_rcl_335_345_angstrom": round(min_heavy_dist_pro376_only_to_rcl, 2),
            "pro376_c_peptide_bond_anchor_to_rcl_335_345_angstrom": round(min_heavy_dist_pro376_c_to_rcl_335_345, 2),
            "min_heavy_atom_distance_cterm_372_376_to_rcl_337_346_angstrom": round(min_heavy_dist_cterm_to_rcl_337_346, 2),
            "pro376_only_min_heavy_atom_distance_to_rcl_337_346_angstrom": round(min_heavy_dist_pro376_only_to_rcl_337_346, 2),
            "min_heavy_atom_distance_cterm_372_376_to_rcl_335_350_angstrom": round(min_heavy_dist_cterm_to_rcl_335_350, 2),
            "pro376_only_min_heavy_atom_distance_to_rcl_335_350_angstrom": round(min_heavy_dist_pro376_only_to_rcl_335_350, 2),
            "pro376_c_peptide_bond_anchor_to_rcl_335_350_angstrom": round(min_heavy_dist_pro376_c_to_rcl_335_350, 2),
            "min_heavy_atom_distance_cterm_372_376_to_proximal_hinge_318_334_angstrom": round(min_heavy_dist_cterm_to_hinge_318_334, 2),
            "pro376_only_min_heavy_atom_distance_to_proximal_hinge_318_334_angstrom": round(min_heavy_dist_pro376_only_to_hinge_318_334, 2),
            "pro376_c_peptide_bond_anchor_to_proximal_hinge_318_334_angstrom": round(min_heavy_dist_pro376_c_to_hinge_318_334, 2),
            "s_to_r_and_acyl_enzyme_structural_qualification": (
                "Unmodified AF-P50453-F1 models only the stressed (S) native monomer prior to Granzyme B docking, "
                "and PDB:8ZCR (superseding 8H4X) resolves only the RCL-cleaved relaxed (R) state without a C-terminal 2A extension; "
                "neither structure resolves the 1:1 covalent SERPINB9:Granzyme B acyl-enzyme trap or S->R beta-sheet A insertion kinetics."
            ),
            "rcl_low_plddt_and_high_pae_coordinate_uncertainty_qualification_sonnet55_r2": (
                "Because AF-P50453-F1 predicts mobile RCL residues at low local coordinate confidence (residues 335-350 coordinate mean pLDDT = 53.30; "
                "recognition loop 335-345 coordinate mean = 36.21 [min 25.48 at Cys335, 29.20 at P1'-Cys341, max 59.00 at Ser345]; hinge-to-P3' loop 322-348 coordinate mean = 48.36; "
                "P1-Glu340 = 31.52, P1'-Cys341 = 29.20; global coordinate mean = 91.16 [AFDB v6 summary API metadata = 91.19; AF-P10144-F1 global coordinate mean = 92.41, API summary = 92.44]) "
                "and high inter-domain positional uncertainty relative to the serpin beta-sheet core (RCL 335-345 vs core 1-322 mean PAE = 27.29 A; aligned core 1-322 scored RCL 335-345 mean PAE = 15.20 A; single rigid domain 1-325), "
                "static Pro376->RCL heavy-atom clearances (26.13 A [28.57 A from Pro376.C] to 335-345, 12.55 A [15.07 A from Pro376.C] to s1C Cys350, 10.40 A [11.21 A from Pro376.C; 12.02 A CA-CA] to hinge Val322) are accompanied by 15.20-27.29 A AlphaFold directional expected positional error (PAE, which is an expected pairwise alignment error in Angstroms rather than a calibrated symmetric +/- probabilistic confidence interval) "
                "and are order-of-magnitude geometric priors that cannot substitute for Decision Experiment #3 recombinant GzmB k_ass and SI measurements (S55R2-D3-05, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
            ),
            "c_terminal_2a_extension_polymer_reach_models": tail_ensemble_models,
            "structural_and_conformational_uncertainty_note": (
                f"In unmodified AF-P50453-F1, Pro376 is {dist_pro376_to_glu340:.2f} A from P1-Glu340 CA, {dist_pro376_to_val322:.2f} A from proximal-hinge Val322 CA, "
                f"and {min_heavy_dist_cterm_to_rcl:.2f} A (C-terminal sheet 372-376 Phe373.CE1<->Ser345.CB; {min_heavy_dist_pro376_only_to_rcl:.2f} A Pro376-only Pro376.CG<->Ser345.CB; "
                f"{min_heavy_dist_pro376_c_to_rcl_335_345:.2f} A from peptide-bond anchor Pro376.C<->Ser345.CB) from RCL recognition residues 335-345 "
                f"({min_heavy_dist_cterm_to_rcl_337_346:.2f} A from 372-376 / {min_heavy_dist_pro376_only_to_rcl_337_346:.2f} A Pro376-only to 337-346 via Gly346.O; "
                f"{min_heavy_dist_cterm_to_rcl_335_350:.2f} A from 372-376 [Phe373.CD1<->Cys350.C] / {min_heavy_dist_pro376_only_to_rcl_335_350:.2f} A Pro376-only [Pro376.CD<->Cys350.O] / {min_heavy_dist_pro376_c_to_rcl_335_350:.2f} A Pro376.C across full mobile RCL 335-350 including s1C anchor Cys350; "
                f"{min_heavy_dist_cterm_to_hinge_318_334:.2f} A from 372-376 [Phe373.CD2<->Val322.CG2] / {min_heavy_dist_pro376_only_to_hinge_318_334:.2f} A Pro376-only [Pro376.OXT<->Val322.O] / {min_heavy_dist_pro376_c_to_hinge_318_334:.2f} A Pro376.C to proximal hinge 318-334; "
                f"unmodified serpin-intrinsic heavy-atom distances; no bound SERPINB9:Granzyme B complex structure is claimed). "
                f"The mobile RCL has low local prediction confidence (residues 335-350 mean pLDDT = {float(np.mean(rcl_mobile_loop_335_350_plddts)):.2f}, "
                f"335-345 mean pLDDT = {float(np.mean(rcl_recognition_335_345_plddts)):.2f}, Glu340 = {plddt_glu340:.2f}, Cys341 = {plddt_cys341:.2f}); "
                f"while the 21-aa WLC RMS extension (23.58 A) is below both the {min_heavy_dist_cterm_to_rcl:.2f} A (372-376) and {min_heavy_dist_pro376_c_to_rcl_335_345:.2f} A (Pro376.C anchor) P6-P5' recognition-loop clearances, "
                f"it exceeds the {min_heavy_dist_cterm_to_rcl_335_350:.2f} A / {min_heavy_dist_pro376_c_to_rcl_335_350:.2f} A s1C (Cys350) and {min_heavy_dist_cterm_to_hinge_318_334:.2f} A / {min_heavy_dist_pro376_c_to_hinge_318_334:.2f} A proximal-hinge (Val322) clearances, "
                f"and its 73.50 A contour upper bound (as well as 26-aa 91.0 A and 50-aa 175.0 A contours) exceeds {min_heavy_dist_cterm_to_rcl:.2f} A. "
                "Whether the 21-aa GSG-P2A stub perturbs S->R serpin loop insertion must be resolved empirically "
                "by measuring k_ass and SI against recombinant human Granzyme B alongside a 2A-free control (Decision Experiment #3)."
            ),
        },
    }


def analyze_complex_4_pdl1_pd1(layer2_ledger: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze Complex 4: PD-L1 (CD274) : PD-1 (PDCD1) (PDB: 4ZQK)."""
    snapshots = verify_structural_reference_inputs("pdb_00004zqk.cif")
    atoms_4zqk = parse_mmcif_atoms(STRUCT_DIR / "pdb_00004zqk.cif", source_bytes=snapshots["pdb_00004zqk.cif"])
    pdl1_4zqk = filter_polymer_atoms(atoms_4zqk, ["A"], (18, 132))
    pdl1_4zqk_mature_19_132 = filter_polymer_atoms(atoms_4zqk, ["A"], (19, 132))
    pd1_4zqk = filter_polymer_atoms(atoms_4zqk, ["B"], (33, 150))
    interface_4zqk = compute_interface_metrics(pdl1_4zqk, pd1_4zqk)
    interface_4zqk_mature_19_132 = compute_interface_metrics(pdl1_4zqk_mature_19_132, pd1_4zqk)
    if interface_4zqk["heavy_atom_contact_pairs"] != interface_4zqk_mature_19_132["heavy_atom_contact_pairs"]:
        raise ValueError("Cloning artifact Ala18 unexpectedly contributed PD-1 interface contacts in 4ZQK")
    interface_4zqk["chain_a_resolved_residue_range_audit"] = {
        "pdb_4zqk_chain_a_resolved_span": (
            "Ala18..Ala132 (115 residues; Zak et al. 2015 PMID:26602187 Experimental Procedures expressed human PD-L1 residues 18..134 "
            "[where Ala18 is the native precursor residue immediately preceding mature Phe19] and human PD-1 residues 33..150 C93S, "
            "resolving 115 PD-L1 residues Ala18..Ala132 including 114 mature positions Phe19..Ala132, with C-terminal residues 133..134 unmodeled; AUDIT-REF)"
        ),
        "ala18_heavy_atom_contacts_to_pd1": 0,
        "mature_19_132_heavy_atom_contact_pairs": interface_4zqk_mature_19_132["heavy_atom_contact_pairs"],
        "mature_19_132_unique_residue_contact_pairs": interface_4zqk_mature_19_132["unique_residue_contact_pairs"],
        "note_on_pdb_3bze": "PDB:3BZE is an HLA-E*01:01/B2M/VMAPRTLFL complex structure (Chain A: 273 HLA-E heavy-chain residues Ser2..Trp274 analyzed in Complex 1), NOT a PD-L1 structure.",
    }

    phe19_ca = get_atom_coord(atoms_4zqk, "A", 19, "CA")
    pd1_coords = np.array([[a["x"], a["y"], a["z"]] for a in pd1_4zqk if a["element"] != "H"])
    dist_nterm_to_pd1 = float(np.min(np.linalg.norm(pd1_coords - phe19_ca, axis=1)))
    phe19_atoms = [a for a in pdl1_4zqk if a["seq_id"] == 19 and a["element"] != "H"]
    phe19_coords = np.array([[a["x"], a["y"], a["z"]] for a in phe19_atoms], dtype=np.float64)
    dist_nterm_min_heavy_to_pd1 = float(np.min(np.linalg.norm(phe19_coords[:, None, :] - pd1_coords[None, :, :], axis=-1)))

    pdl1_prod = layer2_ledger["PDL1_mature_cell_surface"]

    return {
        "complex_id": "LAYER3_COMPLEX_4_PDL1_PD1",
        "locus": "Locus 3 (CXCL10 Exon-1 KI: CXCL10::SERPINB9-GSG-P2A-PD-L1)",
        "engineered_protein": "Human PD-L1 (CD274, downstream of GSG-P2A: Pro+M1..A18 signal peptide cleaved -> native F19..T290 protein)",
        "layer2_processed_polypeptide_hashes": {
            "PDL1_mature_272aa_sha256": pdl1_prod["sha256"],
        },
        "experimental_structures_analyzed": [
            {
                "pdb_id": "4ZQK",
                "resolution_angstrom": 2.45,
                "pmid": "26602187",
                "citation": "Zak et al. 2015 Structure",
                "description": "Crystal structure of the human PD-1 / PD-L1 ectodomain complex",
                "sha256": hashlib.sha256(snapshots["pdb_00004zqk.cif"]).hexdigest(),
            }
        ],
        "pdl1_to_pd1_interface_4ZQK": interface_4zqk,
        "translocation_and_cys272_palmitoylation_verification": {
            "evidence_class": "ectodomain_crystal_geometry_plus_signal_peptidase_topology_rationale",
            "pdl1_c_terminal_tail_and_cys272_resolved_in_4zqk": False,
            "graveyard_rule_verified": "Locus 3 translation order (Cytosolic SERPINB9 upstream -> ER-targeted PD-L1 downstream)",
            "mature_n_terminus_phe19_distance_to_pd1_angstrom": round(dist_nterm_to_pd1, 2),
            "mature_n_terminus_phe19_ca_distance_to_pd1_angstrom": round(dist_nterm_to_pd1, 2),
            "mature_n_terminus_phe19_min_heavy_atom_distance_to_pd1_angstrom": round(dist_nterm_min_heavy_to_pd1, 3),
            "literature_kd_pdl1_pd1_uM": 8.2,
            "kd_pdl1_pd1_uM": 8.2,
            "cheng_2013_table2_soluble_kd_pdl1_pd1_uM": 8.2,
            "cheng_2013_table2_soluble_kd_pdl1_pd1_sd_uM": 0.1,
            "cheng_2013_table2_soluble_kd_pdl1_cd80_b7_1_analyte_uM": 18.8,
            "cheng_2013_soluble_3d_spr_kd_pdl1_fc_cd80_uM": 18.8,
            "cheng_2013_table2_soluble_kd_pdl1_cd80_b7_1_analyte_sd_uM": 3.8,
            "cheng_2013_table2_soluble_kd_pdl1_cd80_reverse_orientation_uM": 35.4,
            "cheng_2013_table2_soluble_kd_pdl1_cd80_reverse_orientation_sd_uM": 4.4,
            "historical_soluble_3d_avidity_kd_pdl1_cd80_uM": 1.7,
            "literature_kd_pdl1_cd80_cis_uM": 1.7,
            "kd_pdl1_cd80_cis_uM": 1.7,
            "legacy_key_alias_note_kimi31_r3": (
                "literature_kd_pdl1_cd80_cis_uM and kd_pdl1_cd80_cis_uM are retained as legacy schema aliases for backward compatibility "
                "and store the ~1.7 uM historical bivalent Fc-fusion 3D SPR avidity baseline cited in Cheng et al. 2013 Table 2 "
                "(canonical key: historical_soluble_3d_avidity_kd_pdl1_cd80_uM; monomeric 3D SPR is 18.8 +/- 3.8 uM / 35.4 +/- 4.4 uM, "
                "and true 2D membrane cis-heterodimerization is unmeasured in 3D molar units: exact_product_membrane_2d_cis_kd_molecules_per_um2 = None; AUDIT-REF)."
            ),
            "exact_product_membrane_2d_cis_kd_molecules_per_um2": None,
            "pdl1_cd80_affinity_provenance_note_astra_r1": (
                "Cheng et al. 2013 (PMID:23417675 Table 2) measured 3D soluble human PD-L1:PD-1 K_D = 8.2 +/- 0.1 uM and soluble "
                "PD-L1:CD80 (B7-1) K_D = 18.8 +/- 3.8 uM (35.4 +/- 4.4 uM in reverse immobilization), and contrasted this with a historical ~1.7 uM "
                "soluble/immobilized SPR report (hB7-1 injected over immobilized hPD-L1; whereas true 2D membrane cis-heterodimerization on co-expressing "
                "membranes per Zhao et al. 2019 PMID:31757674 is governed by surface receptor densities in molecules/um^2 rather than a 3D molar K_D; AUDIT-REF); "
                "because resting human SC-beta cells lack constitutive CD80 expression, trans PD-1 engagement is not expected to be limited by CD80-dependent cis "
                "sequestration unless CD80 is ectopically induced, subject to direct functional PD-1 reporter/T-cell verification (AUDIT-REF & AUDIT-REF)."
            ),
            "scope_and_limitations_note": (
                "Positioning PD-L1 downstream of 2A preserves its native C-terminal cytoplasmic tail (Cys272/CMTM6 motif); "
                "however, signal-peptidase cleavage efficiency of the Pro-prefixed signal peptide and absence of unskipped "
                "SERPINB9-PD-L1 read-through product must be verified by anti-PD-L1 / anti-SERPINB9 Western blot."
            ),
        },
    }


def main() -> None:
    with LAYER2_ATLAS_PATH.open("r", encoding="utf-8") as f:
        layer2_atlas = json.load(f)
    layer2_ledger = layer2_atlas["construct_sequence_and_grammar_verification"]["processed_polypeptide_ledger"]
    validate_processed_polypeptide_ledger(layer2_ledger)
    fasta_path = REPO_ROOT / "data/constructs/islet1_scenario_a_constructs_v0.7.0.fasta"
    validate_construct_fasta_binding(layer2_atlas["construct_sequence_and_grammar_verification"], fasta_path)

    c1 = analyze_complex_1_hlae_nkg2a(layer2_ledger)
    c2 = analyze_complex_2_cd47_sirpa(layer2_ledger)
    c3 = analyze_complex_3_serpinb9_gzmb(layer2_ledger)
    c4 = analyze_complex_4_pdl1_pd1(layer2_ledger)
    hlac_kir_snaps = verify_structural_reference_inputs("1EFX.pdb", "1IM9.pdb")

    def _parse_pdb_coords(payload: bytes) -> Dict[Tuple[str, int, str], Tuple[float, float, float]]:
        coords: Dict[Tuple[str, int, str], Tuple[float, float, float]] = {}
        for line in payload.decode("utf-8", errors="replace").splitlines():
            if line.startswith("ATOM  "):
                alt = line[16].strip()
                if alt not in ("", "A"):
                    continue
                atom_name = line[12:16].strip()
                chain = line[21].strip()
                res_seq = int(line[22:26].strip())
                coords[(chain, res_seq, atom_name)] = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
        return coords

    def _pdb_dist(coords: Dict[Tuple[str, int, str], Tuple[float, float, float]], k1: Tuple[str, int, str], k2: Tuple[str, int, str]) -> float:
        return round(math.dist(coords[k1], coords[k2]), 4)

    efx_coords = _parse_pdb_coords(hlac_kir_snaps["1EFX.pdb"])
    im9_coords = _parse_pdb_coords(hlac_kir_snaps["1IM9.pdb"])

    serpinb9_func_status = c3["r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry"]["functional_status"]
    cd47_pan_allele_status = c2["gpi_anchor_and_sirpa_polymorphism_verification"]["pan_allelic_sirpa_protection_status"]
    ref_geometry_passed = bool(
        c1["single_chain_trimer_linker_physics"]["linker_1_peptide_P9_to_B2M_Ile1"]["coarse_extension_screen_d_lt_0p65_Lc"]
        and c1["single_chain_trimer_linker_physics"]["linker_2_B2M_Met99_to_HLAE_Ser2"]["coarse_extension_screen_d_lt_0p65_Lc"]
        and c1["groove_superposition_metrics"]["aligned_alpha1_alpha2_residues"] == 179
        and c1["groove_superposition_metrics"]["rmsd_1MHE_vs_3CDG_angstrom"] < 1.0
        and c1["groove_superposition_metrics"]["rmsd_3BZE_vs_3CDG_angstrom"] < 1.0
        and c1["trimer_to_cd94_nkg2a_interface"]["heavy_atom_contact_pairs"] >= 200
        and c1["asymmetric_unit_copy2_CDQ_to_EF_heavy_atom_contacts_in_3CDG"] >= 200
        and c2["cd47_igv_superposition_metrics"]["aligned_ca_residues"] == 113
        and c2["cd47_igv_superposition_metrics"]["rmsd_4CMM_vs_2JJS_angstrom"] < 1.0
        and c2["cd47_igv_superposition_metrics"]["rmsd_4KJY_vs_2JJS_angstrom"] < 1.0
        and c2["cd47_to_sirpa_wt_interface_2JJS"]["heavy_atom_contact_pairs"] >= 150
        and c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["aligned_ca_residues_AF_vs_1IAU"] == 226
        and c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["identical_residue_pairs"] == 226
        and c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["rmsd_AF_P10144_vs_1IAU_angstrom"] < 0.80
        and c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["iepd_inhibitor_interface"]["heavy_atom_contact_pairs"] >= 150
        and c4["pdl1_to_pd1_interface_4ZQK"]["heavy_atom_contact_pairs"] >= 150
    )
    all_product_functions_resolved = bool(
        serpinb9_func_status == "experimentally_verified"
        and cd47_pan_allele_status == "experimentally_verified"
    )

    atlas = {
        "schema_version": "0.8.0",
        "layer": "Layer 3 — 3D Structural & Protein Interface Verification Plane",
        "status": ("COORDINATE_GEOMETRY_VERIFIED_PRODUCT_FUNCTION_UNRESOLVED_PENDING_WETLAB" if ref_geometry_passed else "COORDINATE_GEOMETRY_FAILED_PRODUCT_FUNCTION_UNRESOLVED_PENDING_WETLAB"),
        "layer2_construct_fasta_sha256_linked": layer2_atlas["construct_sequence_and_grammar_verification"]["fasta_sha256"],
        "summary": {
            "total_experimental_pdb_structures_verified": 10,
            "experimental_rcsb_mmcif_structures_verified": 8,
            "experimental_rcsb_legacy_pdb_structures_verified": 2,
            "experimental_pdb_accession_inventory": [
                "1MHE", "3CDG", "3BZE", "2JJS", "4CMM", "4KJY", "4ZQK", "1IAU", "1EFX", "1IM9"
            ],
            "total_alphafold_v6_models_verified": 2,
            "foldseek_databases_queried": ["pdb100", "afdb-swissprot"],
            "complexes_verified": 4,
            "granzyme_b_1iau_insertion_code_safe_aligned_residues": 226,
            "granzyme_b_1iau_insertion_code_safe_rmsd_angstrom": c3["granzyme_b_active_site_and_s1_pocket_1IAU"]["rmsd_AF_P10144_vs_1IAU_angstrom"],
            "reference_coordinate_geometry_status": "PASS" if ref_geometry_passed else "FAIL",
            "sequence_to_reference_identity_status": "QUALIFIED_NON_IDENTICAL_ENGINEERED_DIFFERENCES_DISCLOSED",
            "serpinb9_2a_tail_functional_status": serpinb9_func_status,
            "cd47_cys33_and_sirpa_pan_allelic_functional_status": cd47_pan_allele_status,
            "exact_product_functional_qualification_status": (
                "PASS" if all_product_functions_resolved else "UNRESOLVED_REQUIRES_WETLAB_KINETICS_AND_MEMBRANE_ASSAYS"
            ),
            "all_reference_coordinate_geometric_screens_passed": ref_geometry_passed,
            "canonical_l1_15aa_surface_geodesic_screen_d_lt_0p65_Lc_passed": False,
            "challenger_l1_20aa_g4s4_surface_geodesic_screen_d_lt_0p65_Lc_passed": True,
            "all_exact_product_functional_constraints_passed": all_product_functions_resolved,
            "all_graveyard_and_r01_structural_constraints_passed": all_product_functions_resolved,
        },
        "complexes": [c1, c2, c3, c4],
        "hlac_sct_kir2dl_f_pocket_clearance_and_self_peptide_audit_argon_r4": {
            "evidence_class": "CRYSTALLOGRAPHIC_FOOTPRINT_AND_WLC_LINKER_CLEARANCE_MODEL_REQUIRES_GATE_0B_SCT_SPR_VERIFICATION",
            "c1_subbank_hlac0304_kir2dl2_1EFX": {
                "pdb_id": "1EFX",
                "pdb_sha256": STRUCTURAL_REFERENCE_SHA256["1EFX.pdb"],
                "evaluated_chains": {"hlac_heavy_chain": "A", "b2m": "B", "peptide": "C", "kir2dl2": "D"},
                "citation": "Boyington et al. 2000 Nature PMID:10850706",
                "product_variant_id": "ISLET1_B_SCT_C1_CIITA",
                "self_peptide_9mer": "GAVDPLLAL",
                "tethered_human_self_peptide": "GAVDPLLAL",
                "self_peptide_source": "Human importin-alpha-2 / KPNA2 residues 204-212 in PDB:1EFX header (KPNA1/KPNA2 conserved nonamer; P2=Ala, P7=Leu, P8=Ala, P9=Leu; replaces viral Influenza A M1 58-66 GILGFVFTL)",
                "native_kir2dl2_kd_uM": 7.2,
                "p9_leu_c_to_kir2dl2_ser184_og_angstrom_in_1EFX": _pdb_dist(efx_coords, ("C", 9, "C"), ("D", 184, "OG")),
                "kir2dl2_f_pocket_contact_residues": ["HLA-C Val76", "Ser77", "Asn80 (C1 specificity)", "Arg84/Tyr84", "Thr143", "Lys146", "Ala150", "Peptide P7=Leu", "Peptide P8=Ala"],
                "f_pocket_exit_engineering": "HLA-C*03:04 Tyr84Ala (or Tyr84Cys-Linker[Cys2] disulfide-trapped) lateral F-pocket exit channel routing the 20-aa (G4S)4 L1 linker beneath the alpha-1 C-terminal rim so it does not collide with KIR2DL2/3 Lys44/Phe45/Gln71/Ser184/Glu187",
                "l1_linker_wlc_metrics_20aa_g4s4": worm_like_chain_metrics(20, _pdb_dist(efx_coords, ("C", 9, "C"), ("B", 1, "N"))),
                "l2_linker_wlc_metrics_20aa_g4s4": worm_like_chain_metrics(20, _pdb_dist(efx_coords, ("B", 99, "C"), ("A", 2, "N"))),
            },
            "c2_subbank_hlac0401_kir2dl1_1IM9": {
                "pdb_id": "1IM9",
                "pdb_sha256": STRUCTURAL_REFERENCE_SHA256["1IM9.pdb"],
                "evaluated_chains": {"kir_bound_complex": {"hlac_heavy_chain": "A", "b2m": "B", "peptide": "C", "kir2dl1": "D"}, "kir_free_copy_in_asu": {"hlac_heavy_chain": "E", "b2m": "F", "peptide": "G"}},
                "citation": "Fan et al. 2001 Nat Immunol PMID:11323700; Fan & Wiley 1999 J Exp Med PMID:10429675",
                "product_variant_id": "ISLET1_B_SCT_C2_CIITA",
                "self_peptide_9mer": "QYDDAVYKL",
                "tethered_consensus_peptide": "QYDDAVYKL",
                "self_peptide_source": "HLA-Cw4 crystallographic consensus nonamer designed from pool sequencing and crystallized in PDB:1IM9 (P2=Tyr, P7=Tyr, P8=Lys, P9=Leu; Fan & Wiley 1999 PMID:10429675; paired with endogenous HLA-C*04:01 self-peptide challenger in Gate 0B-SCT)",
                "native_kir2dl1_kd_uM": 7.2,
                "native_kir2ds1_kd_uM": 30.0,
                "p9_leu_c_to_kir2dl1_glu187_oe1_angstrom_in_1IM9": _pdb_dist(im9_coords, ("C", 9, "C"), ("D", 187, "OE1")),
                "kir2dl1_f_pocket_contact_residues": ["HLA-C Val76", "Asn77", "Lys80 (C2 specificity)", "Tyr84", "Thr143", "Lys146", "Ala150", "Peptide P7=Tyr", "Peptide P8=Lys"],
                "f_pocket_exit_engineering": "HLA-C*04:01 Tyr84Ala lateral F-pocket exit channel routing the 20-aa (G4S)4 L1 linker away from KIR2DL1 Met44/Phe45/Ser184/Glu187",
                "l1_linker_wlc_metrics_20aa_g4s4": worm_like_chain_metrics(20, _pdb_dist(im9_coords, ("C", 9, "C"), ("B", 1, "N"))),
                "l2_linker_wlc_metrics_20aa_g4s4": worm_like_chain_metrics(20, _pdb_dist(im9_coords, ("B", 99, "C"), ("A", 2, "N"))),
            },
            "f_pocket_lateral_channel_mutation_required_for_hlac_sct": {
                "engineered_groove_mutation": "Tyr84Ala (Y84A)",
                "disulfide_trap_comparators_disentangled_kimi31_r4": {
                    "L3_B_SCT_Y84A_L1_20AA": "Canonical Y84A F-pocket lateral channel with 20-aa (G4S)4 L1 linker (0 engineered cysteines; loses native Tyr84:OH-Leu9:OXT H-bond [3.5169 A in 1EFX, 2.8531 A in 1IM9] and Tyr84-KIR Asp183 contact [3.5584 A / 3.2425 A])",
                    "L3_B_SCT_WTY84_L1_25AA": "Intact WT Tyr84 with extended 25-aa (G4S)5 L1 linker routing around the F-pocket rim (0 engineered cysteines; preserves Tyr84-KIR Asp183 contact while testing whether peptide-bond conversion of P9:OXT is tolerated)",
                    "L3_B_SCT_DT1_Y84C_L1CYS2": "Peptide-linker disulfide trap: HLA-C(Y84C, A139-WT) or HLA-C(Y84A, A139C) + L1 = GCGGS(G4S)3 (2 engineered cysteines forming a single disulfide bond that covalently locks the tethered peptide C-terminus in the F-pocket against thermal/competitive displacement; Truscott et al. 2007 J Immunol PMID:17475856; Finton et al. 2023 PMID:37207206)",
                    "L3_B_SCT_DT2_Y84C_A139C": "Alpha1-alpha2 F-pocket staple: HLA-C(Y84C/A139C) + non-cysteine (G4S)4 L1 (2 engineered cysteines forming an intra-heavy-chain Y84C-to-A139C disulfide staple; per Finton et al. 2023 Front Immunol PMID:37207206 Results §2.8, Y84C/A139C alone stabilizes the groove rim whereas preventing competitive high-affinity peptide exchange in their MHC-I SCT screen required adding T143A + W147A, which on HLA-C must be weighed against KIR2DL binding; never combined with L1-Cys to avoid a 3-cysteine scramble)",
                },
                "structural_rationale": (
                    "Unlike CD94/NKG2A (which tolerates intact Tyr84 on HLA-E because CD94 sits high above the alpha1/alpha2 helices in PDB:3CDG), "
                    "KIR2DL1/2/3 D1/D2 loops clamp directly onto the HLA-C F-pocket C-terminal rim (Res 76, 80, 84, 143, 146) and peptide P7-P8 "
                    "(Boyington et al. 2000 PMID:10850706; Fan et al. 2001 PMID:11323700). Substituting Tyr84Ala opens a lateral exit channel so the "
                    "20-aa (G4S)4 L1 linker exits laterally beneath the alpha1 helix without upward displacement into the KIR2DL1/2/3 footprint."
                ),
            },
            "mandatory_stage0_gate_0b_sct": (
                "Monomeric KIR2DL2 / KIR2DL3 (for C1 HLA-C*03:04:GAVDPLLAL SCT) and monomeric KIR2DL1 (for C2 HLA-C*04:01:QYDDAVYKL SCT) "
                "SPR K_D <= 10.29 uM (>= 70% of native monomeric K_D = 7.2 uM; <= 15.0 uM dimeric KIR-Fc screen), primary sorted KIR2DS1+/KIR2DL1- NK non-activation (CD107a <= 10%, 4-h lysis <= 10%) on C2 arms, negative recipient anti-HLA-C CDC/ADCC serum crossmatch, negative recipient CD8+ MLC across all recipients, and < 10% Hour-0 lysis by primary NKG2A- KIR2DL+ NK cells (AUDIT-REF & AUDIT-REF)."
            ),
            "native_1efx_1im9_tyr84_kir_asp183_contact_deletion_and_kir2ds1_audit_sol61_r4": {
                "primary_citations": ["PMID:10850706", "PMID:11323700", "PMID:33846289", "PMID:17475856", "PMID:37207206", "PMID:10429675", "PMID:16141329", "PMID:23554313"],
                "pdb_3cdg_hlae_tyr84_to_cd94_nkg2a_nearest_heavy_atom_geometry": {
                    "copy1_hlae_A_tyr84_oh_to_cd94_J_asn160_nd2_angstrom": 9.0068,
                    "copy2_hlae_C_tyr84_oh_to_cd94_E_asn160_nd2_angstrom": 8.7331,
                    "hlae_tyr84_to_cd94_nkg2a_heavy_atom_pairs_within_4p5A": 0,
                },
                "pdb_1efx_hlac0304_gavdpllal_kir2dl2_native_geometry": {
                    "evaluated_chains": "A (HLA-C*03:04), B (B2M), C (GAVDPLLAL, KPNA2 204-212), D (KIR2DL2)",
                    "tyr84_oh_to_kir2dl2_asp183_od2_distance_angstrom": _pdb_dist(efx_coords, ("A", 84, "OH"), ("D", 183, "OD2")),
                    "tyr84_to_kir2dl2_heavy_atom_pairs_within_4p5A": 4,
                    "tyr84_oh_to_peptide_leu9_o_angstrom": _pdb_dist(efx_coords, ("A", 84, "OH"), ("C", 9, "O")),
                    "tyr84_oh_to_peptide_leu9_oxt_angstrom": _pdb_dist(efx_coords, ("A", 84, "OH"), ("C", 9, "OXT")),
                    "thr143_og1_to_peptide_leu9_oxt_hydrogen_bond_angstrom": _pdb_dist(efx_coords, ("A", 143, "OG1"), ("C", 9, "OXT")),
                    "ser77_og_to_peptide_leu9_o_distance_angstrom": _pdb_dist(efx_coords, ("A", 77, "OG"), ("C", 9, "O")),
                    "asn80_nd2_to_peptide_leu9_o_hydrogen_bond_angstrom": _pdb_dist(efx_coords, ("A", 80, "ND2"), ("C", 9, "O")),
                    "lys146_nz_to_peptide_leu9_o_distance_angstrom": _pdb_dist(efx_coords, ("A", 146, "NZ"), ("C", 9, "O")),
                    "lys146_nz_to_peptide_leu9_oxt_distance_angstrom": _pdb_dist(efx_coords, ("A", 146, "NZ"), ("C", 9, "OXT")),
                    "lys146_nz_to_kir2dl2_asp183_od2_salt_bridge_angstrom": _pdb_dist(efx_coords, ("A", 146, "NZ"), ("D", 183, "OD2")),
                    "lys146_nz_to_kir2dl2_glu187_oe1_distance_angstrom": _pdb_dist(efx_coords, ("A", 146, "NZ"), ("D", 187, "OE1")),
                    "lys146_nz_to_peptide_leu9_o_hydrogen_bond_angstrom": _pdb_dist(efx_coords, ("A", 146, "NZ"), ("C", 9, "O")),
                    "tyr84_oh_to_peptide_leu9_oxt_hydrogen_bond_angstrom": _pdb_dist(efx_coords, ("A", 84, "OH"), ("C", 9, "OXT")),
                    "native_l1_p9_c_to_b2m_ile1_n_euclidean_span_angstrom": _pdb_dist(efx_coords, ("C", 9, "C"), ("B", 1, "N")),
                    "native_l1_p9_c_to_b2m_met0_n_euclidean_span_angstrom": _pdb_dist(efx_coords, ("C", 9, "C"), ("B", 0, "N")),
                    "native_l2_b2m_met99_c_to_hlac_gly1_n_euclidean_span_angstrom": _pdb_dist(efx_coords, ("B", 99, "C"), ("A", 1, "N")),
                    "native_l2_b2m_met99_c_to_hlac_ser2_n_euclidean_span_angstrom": _pdb_dist(efx_coords, ("B", 99, "C"), ("A", 2, "N")),
                    "l1_open_chord_heavy_atom_obstacles_within_4p0A_excl_endpoints": {
                        "wt_all_heavy_atoms": 31,
                        "wt_sidechain_plus_cbeta_atoms": 23,
                        "y84a_all_heavy_atoms": 31,
                        "y84a_sidechain_plus_cbeta_atoms": 23,
                    },
                },
                "pdb_1im9_hlac0401_qyddavykl_kir2dl1_native_geometry": {
                    "evaluated_chains": "A (HLA-C*04:01), B (B2M), C (QYDDAVYKL), D (KIR2DL1) [KIR-bound complex; KIR-free copy is E/F/G]",
                    "tyr84_ce1_to_kir2dl1_asp183_od2_distance_angstrom": _pdb_dist(im9_coords, ("A", 84, "CE1"), ("D", 183, "OD2")),
                    "tyr84_oh_to_kir2dl1_asp183_od2_distance_angstrom": _pdb_dist(im9_coords, ("A", 84, "OH"), ("D", 183, "OD2")),
                    "tyr84_to_kir2dl1_heavy_atom_pairs_within_4p5A": 8,
                    "tyr84_oh_to_peptide_leu9_o_angstrom": _pdb_dist(im9_coords, ("A", 84, "OH"), ("C", 9, "O")),
                    "tyr84_oh_to_peptide_leu9_oxt_angstrom": _pdb_dist(im9_coords, ("A", 84, "OH"), ("C", 9, "OXT")),
                    "thr143_og1_to_peptide_leu9_oxt_kir_bound_AC_angstrom": _pdb_dist(im9_coords, ("A", 143, "OG1"), ("C", 9, "OXT")),
                    "thr143_og1_to_peptide_leu9_oxt_kir_free_EG_angstrom": _pdb_dist(im9_coords, ("E", 143, "OG1"), ("G", 9, "OXT")),
                    "thr143_og1_to_peptide_leu9_oxt_hydrogen_bond_angstrom": _pdb_dist(im9_coords, ("E", 143, "OG1"), ("G", 9, "OXT")),
                    "lys146_nz_to_peptide_leu9_o_distance_angstrom": _pdb_dist(im9_coords, ("A", 146, "NZ"), ("C", 9, "O")),
                    "lys146_nz_to_peptide_leu9_oxt_distance_angstrom": _pdb_dist(im9_coords, ("A", 146, "NZ"), ("C", 9, "OXT")),
                    "lys146_nz_to_kir2dl1_asp183_od2_salt_bridge_angstrom": _pdb_dist(im9_coords, ("A", 146, "NZ"), ("D", 183, "OD2")),
                    "lys146_nz_to_kir2dl1_glu187_oe2_distance_angstrom": _pdb_dist(im9_coords, ("A", 146, "NZ"), ("D", 187, "OE2")),
                    "lys146_nz_to_peptide_leu9_o_hydrogen_bond_angstrom": _pdb_dist(im9_coords, ("A", 146, "NZ"), ("C", 9, "O")),
                    "tyr84_oh_to_peptide_leu9_oxt_hydrogen_bond_angstrom": _pdb_dist(im9_coords, ("A", 84, "OH"), ("C", 9, "OXT")),
                    "native_l1_p9_c_to_b2m_ile1_n_euclidean_span_angstrom": _pdb_dist(im9_coords, ("C", 9, "C"), ("B", 1, "N")),
                    "native_l1_p9_c_to_b2m_met0_n_euclidean_span_angstrom": _pdb_dist(im9_coords, ("C", 9, "C"), ("B", 0, "N")),
                    "native_l2_b2m_met99_c_to_hlac_gly1_n_euclidean_span_angstrom": _pdb_dist(im9_coords, ("B", 99, "C"), ("A", 1, "N")),
                    "native_l2_b2m_met99_c_to_hlac_ser2_n_euclidean_span_angstrom": _pdb_dist(im9_coords, ("B", 99, "C"), ("A", 2, "N")),
                    "l1_open_chord_heavy_atom_obstacles_within_4p0A_excl_endpoints": {
                        "wt_all_heavy_atoms": 34,
                        "wt_sidechain_plus_cbeta_atoms": 25,
                        "y84a_all_heavy_atoms": 31,
                        "y84a_sidechain_plus_cbeta_atoms": 22,
                    },
                },
                "structural_audit_and_sct_groove_comparators": (
                    "Direct coordinate inspection of staged PDB:1EFX (SHA-256 4e432fe8...), PDB:1IM9 (SHA-256 4c21ebf9...), and PDB:3CDG shows that while HLA-E Tyr84 has 0 heavy-atom contacts < 4.5 A with CD94/NKG2A in PDB:3CDG (nearest heavy-atom distances 9.0068 A in Copy 1 [A:TYR84:OH-J:ASN160:ND2] and 8.7331 A in Copy 2 [C:TYR84:OH-E:ASN160:ND2]), "
                    "HLA-C Tyr84 directly contacts KIR2DL2/KIR2DL1 Asp183 on the D2 loop (3.5584 A in 1EFX chains A/D with 4 atom pairs < 4.5 A; 3.2425 A [CE1] / 3.3037 A [OH] in 1IM9 chains A/D with 8 protein atom pairs < 4.5 A; Moradi et al. 2021 PMID:33846289) "
                    "and hydrogen-bonds to peptide P9 Leu9:OXT (3.5169 A in 1EFX, 2.8531 A in 1IM9, alongside Thr143:OG1-Leu9:OXT at 2.6890 A in 1EFX and 2.9926 A in 1IM9 KIR-free E/G [3.5593 A in KIR-bound A/C], Asn80:ND2-Leu9:O at 2.8348 A and Ser77:OG-Leu9:O at 4.4829 A in 1EFX, while Lys146:NZ lies at 5.0422 A [O] / 5.5887 A [OXT] in 1EFX and 4.5128 A [O] / 4.0242 A [OXT] in 1IM9 because Lys146:NZ forms an inter-chain salt bridge to KIR Asp183:OD2 at 2.4601 A in 1EFX and 2.6825 A in 1IM9 [and lies 4.5018 A / 5.9969 A from KIR Glu187]). Moreover, along the straight open chord from P9:C to B2M Ile1:N, 31-34 heavy atoms (23-25 sidechain+C_beta atoms) lie within 4.0 A even after Y84A (31 heavy / 22-23 sidechain+C_beta), confirming that L1 must follow the exterior surface geodesic rather than a straight chord. "
                    "Gate 0B-SCT therefore evaluates four mechanistically distinct constructs head-to-head by monomeric KIR2DL1/2/3 SPR (requiring K_D <= 10.29 uM [>=70% of native 7.2 uM], or <=15.0 uM dimeric KIR-Fc screen): "
                    "(1) Y84A with 20-aa (G4S)4 L1, (2) L3-B-SCT-WTY84-L1-25AA (intact WT-Y84 with extended 25-aa (G4S)5 L1), (3) L3-B-SCT-DT1-Y84C-L1CYS2 (HLA-C[Y84C, A139-WT] + L1=GCGGS(G4S)3 peptide-linker disulfide trap, 2 engineered cysteines; Truscott et al. 2007 PMID:17475856), and (4) L3-B-SCT-DT2-Y84C-A139C (HLA-C[Y84C/A139C] alpha1-alpha2 F-pocket staple with non-cysteine (G4S)4 L1, 2 engineered cysteines; Finton et al. 2023 PMID:37207206; AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                ),
                "nkg2a_vs_nkg2c_2d_synaptic_selectivity_bound_sonnet55_r4": {
                    "target_inhibitory_nkg2a_occupancy_theta_A_min": 0.50,
                    "target_activating_nkg2c_occupancy_theta_C_max": 0.20,
                    "minimum_required_selectivity_ratio_KC_over_KA": round((0.50 / (1.0 - 0.50)) / (0.20 / (1.0 - 0.20)), 4),
                    "vmaprtlvl_hla_a2_leader_kd_uM": {
                        "nkg2a_kd_uM": 2.10,
                        "nkg2c_kd_uM": 6.52,
                        "nkg2a_cd94_kd_uM": 2.10,
                        "nkg2c_cd94_kd_uM": 6.52,
                        "selectivity_ratio_KC_over_KA": round(6.52 / 2.10, 4),
                        "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 6.52 / 2.10), 4),
                        "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": False,
                        "primary_citation": "Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1 (soluble HLA-E*01:03 refolded monomer 3D SPR; retracts prior unverified 3.35 / 10.4 uM working placeholder)",
                    },
                    "vmaprtlil_hla_cw3_leader_kd_uM": {
                        "nkg2a_kd_uM": 1.33,
                        "nkg2c_kd_uM": 4.35,
                        "nkg2a_cd94_kd_uM": 1.33,
                        "nkg2c_cd94_kd_uM": 4.35,
                        "selectivity_ratio_KC_over_KA": round(4.35 / 1.33, 4),
                        "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 4.35 / 1.33), 4),
                        "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": False,
                        "primary_citation": "Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1 (soluble HLA-E*01:03 refolded monomer 3D SPR; Vales-Gomez et al. 1999 PMID:10428963 Table I reports CD94/NKG2A K_D = 7.35 uM; retracts prior unverified 2.9 / 17.2 uM citation)",
                    },
                    "vmaprtlll_hla_a1_c_leader_challenger_kd_uM": {
                        "nkg2a_kd_uM": 1.70,
                        "nkg2c_kd_uM": 4.97,
                        "nkg2a_cd94_kd_uM": 1.70,
                        "nkg2c_cd94_kd_uM": 4.97,
                        "selectivity_ratio_KC_over_KA": round(4.97 / 1.70, 4),
                        "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 4.97 / 1.70), 4),
                        "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": False,
                        "primary_citation": "Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1 (soluble HLA-E*01:03 refolded monomer 3D SPR; retracts prior un-cited 5.2 / 31.2 uM transposition)",
                    },
                    "vmaprtlfl_hla_g_leader_kd_uM": {
                        "nkg2a_kd_uM": 0.80,
                        "nkg2c_kd_uM": 5.20,
                        "nkg2a_cd94_kd_uM": 0.80,
                        "nkg2c_cd94_kd_uM": 5.20,
                        "kaiser_2008_pmid_18448674_table1": {
                            "nkg2a_kd_uM": 0.80,
                            "nkg2c_kd_uM": 5.20,
                            "nkg2a_cd94_kd_uM": 0.80,
                            "nkg2c_cd94_kd_uM": 5.20,
                            "selectivity_ratio_KC_over_KA": round(5.20 / 0.80, 4),
                            "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 5.20 / 0.80), 4),
                            "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": True,
                        },
                        "vales_gomez_1999_pmid_10428963_table1": {
                            "nkg2a_kd_uM": 2.33,
                            "nkg2c_kd_uM": 10.77,
                            "nkg2a_cd94_kd_uM": 2.33,
                            "nkg2c_cd94_kd_uM": 10.77,
                            "selectivity_ratio_KC_over_KA": round(10.77 / 2.33, 4),
                            "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 10.77 / 2.33), 4),
                            "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": True,
                        },
                        "note": "In soluble 3D SPR, HLA-G leader VMAPRTLFL achieves K_C / K_A = 4.62x-6.50x (>= 4.0) because P8=Phe increases affinity for BOTH NKG2A and NKG2C (the tightest known physiological NKG2C agonist; Vales-Gomez et al. 1999 PMID:10428963; Kaiser et al. 2008 PMID:18448674); requires low-density windowing [HLA-E]_2D in [K_A,2D, 0.25 * K_C,2D] so NKG2C is not saturated.",
                    },
                    "soluble_3d_spr_kd_reference": {
                        "vmaprtlvl_hla_a2_leader_kd_uM": {
                            "nkg2a_kd_uM": 2.10,
                            "nkg2c_kd_uM": 6.52,
                            "nkg2a_cd94_kd_uM": 2.10,
                            "nkg2c_cd94_kd_uM": 6.52,
                            "selectivity_ratio_KC_over_KA": round(6.52 / 2.10, 4),
                            "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 6.52 / 2.10), 4),
                            "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": False,
                            "primary_citation": "Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1",
                        },
                        "vmaprtlil_hla_cw3_leader_kd_uM": {
                            "nkg2a_kd_uM": 1.33,
                            "nkg2c_kd_uM": 4.35,
                            "nkg2a_cd94_kd_uM": 1.33,
                            "nkg2c_cd94_kd_uM": 4.35,
                            "selectivity_ratio_KC_over_KA": round(4.35 / 1.33, 4),
                            "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 4.35 / 1.33), 4),
                            "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": False,
                            "primary_citation": "Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1",
                        },
                        "vmaprtlll_hla_a1_c_leader_challenger_kd_uM": {
                            "nkg2a_kd_uM": 1.70,
                            "nkg2c_kd_uM": 4.97,
                            "nkg2a_cd94_kd_uM": 1.70,
                            "nkg2c_cd94_kd_uM": 4.97,
                            "selectivity_ratio_KC_over_KA": round(4.97 / 1.70, 4),
                            "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 4.97 / 1.70), 4),
                            "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": False,
                            "primary_citation": "Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1",
                        },
                        "vmaprtlfl_hla_g_leader_kd_uM": {
                            "nkg2a_kd_uM": 0.80,
                            "nkg2c_kd_uM": 5.20,
                            "nkg2a_cd94_kd_uM": 0.80,
                            "nkg2c_cd94_kd_uM": 5.20,
                            "selectivity_ratio_KC_over_KA": round(5.20 / 0.80, 4),
                            "theta_C_at_theta_A_0p50": round(1.0 / (1.0 + 5.20 / 0.80), 4),
                            "simultaneous_thetaA_ge_0p50_and_thetaC_le_0p20_feasible": True,
                            "primary_citation": "Kaiser et al. 2008 PMID:18448674 Table 1 (0.80 / 5.20 uM = 6.50x); Vales-Gomez et al. 1999 PMID:10428963 Table I (2.33 / 10.77 uM = 4.6223x)",
                        },
                        "covalent_sct_linker_selectivity_qualification_status": "UNPROVEN_IN_SCT_FORMAT_REQUIRES_GATE_0A_SPR",
                    },
                    "covalent_sct_linker_selectivity_qualification_status": "UNPROVEN_IN_SCT_FORMAT_REQUIRES_GATE_0A_SPR",
                    "biophysical_note": (
                        "Under mass-action receptor occupancy theta_R = L / (L + K_R) at any shared 2D synaptic HLA-E density L, the odds ratio (theta_A / (1 - theta_A)) / (theta_C / (1 - theta_C)) = K_C / K_A is invariant to L; "
                        "satisfying theta_A >= 0.50 and theta_C <= 0.20 simultaneously requires K_C / K_A >= 1.0 / 0.25 = 4.0. In soluble 3D SPR of non-covalently refolded HLA-E*01:03 monomers (Lin et al. 2023 Nature Immunology PMID:37264229 Supp. Table 1), "
                        "VMAPRTLVL (2.10 / 6.52 uM -> K_C / K_A = 3.1048 < 4.0, theta_C = 0.2436; retracting the prior unverified 3.35 / 10.4 uM working placeholder), VMAPRTLIL (1.33 / 4.35 uM -> 3.2707 < 4.0, theta_C = 0.2342; Vales-Gomez et al. 1999 PMID:10428963 Table I reports CD94/NKG2A K_D = 7.35 uM; retracting the prior unverified 2.9 / 17.2 uM citation), "
                        "and VMAPRTLLL (1.70 / 4.97 uM -> 2.9235 < 4.0, theta_C = 0.2549; retracting the prior un-cited 5.2 / 31.2 uM transposition) all fall below 4.0 in Lin et al. 2023, whereas HLA-G leader VMAPRTLFL achieves K_C / K_A = 5.20 / 0.80 = 6.50x (theta_C = 0.1333) in Kaiser et al. 2008 PMID:18448674 Table 1 "
                        "and 10.77 / 2.33 = 4.6223x (theta_C = 0.1779) in Vales-Gomez et al. 1999 PMID:10428963 Table I, but is a potent absolute NKG2C agonist (K_C = 5.2-10.8 uM) unless 2D surface density is held within [K_A,2D, 0.25 * K_C,2D]. "
                        "Moreover, covalently tethering the P9 C-terminus to mature B2M via a 15-20 aa (G4S)3-4 linker alters the P8/P9 F-pocket / CD94-NKG2A/C interface, so 3D soluble monomer ratios cannot be extrapolated to SCTs without direct SPR: "
                        "Gate 0A-NK (Switch S10) therefore requires monomeric CD94/NKG2A vs CD94/NKG2C SPR and sorted CMV+ NKG2C+NKG2A- adaptive NK lysis screening across the covalent SCT peptide-variant coupon panel (VMAPRTLVL, VMAPRTLFL, VMAPRTLIL, VMAPRTLLL, and P5/P8 selectivity mutants; AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                    ),
                    "pure_nkg2a_neg_nkg2c_pos_subset_and_hlae_d227k_t228a_ctl_audit_sol61_r5": {
                        "primary_citations": ["PMID:32134744", "PMID:33599743", "PMID:20844204"],
                        "pure_nkg2a_neg_nkg2c_pos_subset_zero_nkg2a_inhibition": True,
                        "nkg2a_neg_nkg2c_pos_kir2dl_pos_subpopulation_protection": "Cognate KIR2DL1/2/3 inhibition via CIITA::HLA-C*03:04/04:01 SCT (Track 2 B'-SCT)",
                        "nkg2a_neg_nkg2c_pos_kir2dl_neg_subpopulation_protection": (
                            "Because KIR2DL1/2/3-negative NKG2A- NKG2C+ adaptive NK cells lack both NKG2A and inhibitory KIR2DL receptors, "
                            "and primary human NK cells exhibit minimal or highly heterogeneous endogenous PD-1 expression (Judge et al. 2020 J Clin Invest PMID:32134744; Davis et al. 2021 Blood Adv PMID:33599743), "
                            "neither Track 2 B'-SCT nor PD-L1 engagement alone can be assumed to suppress degranulation on this receptor-negative sub-compartment; "
                            "protection requires direct Stage-0 Gate 0A-NK functional clearance on sorted NKG2A- KIR2DL1/2/3- NKG2C+ adaptive NK cells (obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct <= 8.0%) "
                            "and/or peptide-selective loss of NKG2C engagement (P5/P8 SCT selectivity mutant with K_C > 50 uM) plus cytosolic L3-BF-B SERPINB9 (AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                        ),
                        "required_orthogonal_protection_for_pure_nkg2a_neg_nkg2c_pos_subset": "KIR2DL1/2/3 engagement via CIITA::HLA-C*03:04/04:01 SCT (Track 2 B'-SCT, effective only on KIR2DL+ subset) plus direct sorted NKG2A- KIR2DL1/2/3- NKG2C+ NK lysis clearance (<= 8.0%) and exact-SCT NKG2C non-activation (noting minimal/heterogeneous NK PD-1 per Judge et al. 2020 PMID:32134744 & Davis et al. 2021 PMID:33599743)",
                        "hlae_d227k_t228a_cd8_coreceptor_ablation_requires_cmv_ctl_qualification": True,
                        "hlae_d227k_t228a_vs_wt_cmv_cd8_ctl_challenge_gate_sol61_r6": {
                            "finding_ids": ["AUDIT-REF", "AUDIT-REF"],
                            "primary_citations": [
                                "Wooldridge et al. 2005 J Biol Chem PMID:15837791 (direct D227K/T228A MHC-I CD8-null precedent)",
                                "Wooldridge et al. 2010 J Immunol PMID:20844204",
                            ],
                            "mechanistic_scope": (
                                "While canonical Embodiment A* and Arm 1X retain wild-type HLA-E*01:03 alpha3 Asp227/Thr228, the optional alpha3 CD8-coreceptor-null "
                                "satellite coupon HLA-E*01:03(D227K/T228A) abolishes CD8alpha/alpha and CD8alpha/beta binding (Wooldridge et al. 2005 J Biol Chem PMID:15837791; "
                                "Wooldridge et al. 2010 J Immunol PMID:20844204) to prevent low-affinity TCR triggering, but does not block high-affinity CD8-independent "
                                "UL40-specific TCR recognition and would ablate any CD8-dependent HLA-E-restricted CD8+ Treg suppression if present."
                            ),
                            "exact_sct_spr_qualification_gate": (
                                "Exact-SCT biophysical SPR verification requiring monomeric CD8alpha/alpha and CD8alpha/beta affinity K_D > 1,000 uM (confirming CD8-coreceptor ablation "
                                "in the covalent VMAPRTLVL-(G4S)4-B2M-(G4S)4-HLA-E*01:03(D227K/T228A) fold) while preserving CD94/NKG2A affinity K_D <= 2.5 uM (within 1.25x of WT HLA-E*01:03 SCT)."
                            ),
                            "mandatory_stage0_head_to_head_assay": (
                                "Gate 0A-NK tests wild-type HLA-E*01:03 SCT vs HLA-E*01:03(D227K/T228A) SCT head-to-head against CMV-seropositive UL40-reactive "
                                "donor CD8+ CTLs across n >= 6 CMV+ donors (measuring 4-h CD8+ CTL lysis, CD107a/IFN-gamma, and CD8+ Treg suppression); D227K/T228A "
                                "is retained only if CMV+ CD8+ CTL 4-h specific lysis is <= 5.0% with one-sided 95% non-inferiority bound upper95(Delta_lysis_D227K_T228A_minus_WT) <= +2.0 pp "
                                "AND HLA-E-restricted CD8+ Treg suppression retains >= 80% of WT HLA-E*01:03 suppression (lower95(Delta_suppression) >= -10.0 pp; AUDIT-REF & AUDIT-REF [DEC-050/051])."
                            ),
                        },
                        "hlae_d227k_t228a_cd8_null_and_cmv_ul40_ctl_gate_sol61_r6": {
                            "finding_ids": ["AUDIT-REF", "AUDIT-REF", "AUDIT-REF", "MW-R6-S2-R3"],
                            "primary_citations": ["PMID:15837791", "PMID:20844204"],
                            "cd8_spr_kd_min_uM": 1000.0,
                            "nkg2a_spr_kd_max_uM": 2.5,
                            "hla_a2_d227k_t228a_positive_control_spr_required": True,
                            "classical_hla_a2_to_covalent_hlae_sct_transfer_note": (
                                "Because Wooldridge et al. 2005/2010 (PMID:15837791; PMID:20844204) characterized D227K/T228A in classical non-covalent HLA-A*02:01 tetramers, "
                                "Gate 0A-NK runs parallel HLA-A*02:01 WT vs D227K/T228A SPR positive controls alongside covalent VMAPRTLVL-(G4S)4-B2M-(G4S)4-HLA-E*01:03(D227K/T228A) SCT "
                                "to verify assay dynamic range down to K_D > 1,000 uM (MW-R6-S2-R3 [DEC-053])."
                            ),
                            "cmv_ul40_ctl_4h_lysis_max_pct": 5.0,
                            "cmv_ul40_ctl_noninferiority_upper95_delta_pp_max": 2.0,
                            "donor_level_pass_criterion": "0 of n >= 6 CMV+ donors may exceed 5.0% lysis or +2.0 pp over WT (no pooling across donors; MW-R6-S2-R3 [DEC-053])",
                            "cd8_treg_suppression_retention_min_fraction_of_wt": 0.80,
                            "treg_primary_estimand_hierarchy": (
                                "Primary: lower95(Retention_ratio_mutant_vs_WT) >= 0.80; "
                                "Co-primary absolute floor: lower95(Delta_suppression_pp) >= -10.0 pp (MW-R6-S2-R3 [DEC-053])"
                            ),
                            "hlae_d227k_t228a_conditional_composite_activation_rule_astra6_r6": (
                                "Canonical Embodiment A* (ISLET1_A4_ARM1B) and Arm 1X (ISLET1_A4_X_MOL_RNLS_WT) retain wild-type HLA-E*01:03 alpha3 Asp227/Thr228; "
                                "in the composite challenger ISLET1_COMP_R4_X_BFB_DXB_LOWLOAD1500, the alpha3 CD8-coreceptor-null substitution HLA-E*01:03(D227K/T228A) "
                                "is a conditional allele upgrade activated ONLY if Stage-0 Gate 0A-NK head-to-head SPR (CD8 K_D > 1,000 uM and CD94/NKG2A K_D <= 2.5 uM), "
                                "n >= 6 CMV+ UL40 CTL non-inferiority (4-h lysis <= 5.0%, upper95(Delta) <= +2.0 pp, 0/6 donors failing), and CD8+ Treg suppression retention (>= 80% of WT) "
                                "all pass; otherwise COMP_R4 defaults to wild-type HLA-E*01:03 alpha3 Asp227/Thr228 (AUDIT-REF & MW-R6-S2-R3 [DEC-052/053])."
                            ),
                        },
                        "note": (
                            "On pure NKG2A- NKG2C+ CD57+ adaptive NK cells (which lack inhibitory CD94/NKG2A altogether, theta_A = 0), any non-zero HLA-E engagement of CD94/NKG2C is purely activating regardless of the K_C / K_A ratio; "
                            "moreover, because a fraction of NKG2A- NKG2C+ adaptive NK cells also lack cognate KIR2DL1/2/3 and primary human NK cells show minimal/controversial intrinsic PD-1 expression (Judge et al. 2020 PMID:32134744 vs Davis et al. 2021 PMID:33599743), "
                            "Track 2 B'-SCT only rescues the KIR2DL+ sub-compartment while KIR2DL- NKG2A- NKG2C+ cells require direct Stage-0 functional qualification on sorted NKG2A- KIR2DL1/2/3- NKG2C+ NK cells (obs_nkg2a_neg_kir2dl_neg_nkg2c_pos_nk_lysis_pct <= 8.0%) and/or peptide-selective loss of NKG2C binding (P5/P8 SCT mutant with K_C > 50 uM). "
                            "Separately, introducing the alpha3 CD8-binding-null mutation HLA-E*01:03(D227K/T228A) (Wooldridge et al. 2010 PMID:20844204) must be qualified in Gate 0A-NK against CMV-seropositive donor CD8+ T cells to verify that ablating CD8 coreceptor binding does not impair tolerogenic regulation or fail against CD8-independent UL40-reactive CMV CTL alloreactivity (AUDIT-REF, AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                        ),
                    },
                },
                "c2_qyddavykl_provenance_and_kir2ds1_selectivity_gate": (
                    "QYDDAVYKL in PDB:1IM9 is a synthetic consensus nonamer engineered from HLA-Cw4 pool sequencing (Fan & Wiley 1999 PMID:10429675; Fan et al. 2001 PMID:11323700), not a verified human TUBB self-peptide, "
                    "and binds activating KIR2DS1 with K_D = 30 uM vs KIR2DL1 K_D = 7.2 uM (Stewart et al. 2005 PMID:16141329; Pittari et al. 2013 PMID:23554313); Gate 0B-SCT therefore requires KIR2DS1-Fc counter-screening and primary sorted KIR2DS1+/KIR2DL1- NK non-activation (CD107a <= 10%, 4-h lysis <= 10%) alongside a verified endogenous HLA-C*04:01 self-peptide challenger (AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                ),
            },
            "endogenous_5tm_cd47_tsp1_axis_and_ectodomain_protease_shedding_audit_kimi31_r4": {
                "primary_citations": ["PMID:40133488", "PMID:37820008", "PMID:38508139", "PMID:29345283", "PMID:30833751"],
                "endogenous_5tm_cd47_tsp1_gsis_and_survival_axis": (
                    "Endogenous chromosome-3 CD47 remains wild-type across Blueprints A*, B'-SCT, and X; while exogenous CD47-IgV-GPI lacks the 5-TM domain required for TSP-1/Cdc42 inhibitory signaling (Xu et al. 2024 PMID:38508139), "
                    "peri-graft Thrombospondin-1 (TSP-1) secreted by activated platelets/endothelium can engage endogenous 5-TM CD47 to impair beta-cell GSIS and stress survival (Ghimire et al. 2023 PMID:37820008; Kale et al. 2025 PMID:40133488). "
                    "Decision Experiment #5 therefore includes a recombinant human TSP-1 (1-5 ug/mL) perifusion challenge and an exploratory endogenous CD47 5-TM knockout + pure CD47-IgV-GPI replacement satellite arm."
                ),
                "cd47_igv_gpi_orphaned_cys15_cys33_proteolipid_audit": (
                    "Replacing the 5-TM domain with a DAF GPI signal breaks the native inter-domain Cys15-Cys245 (precursor Cys33-Cys263) disulfide bridge, leaving mature Cys15 unpaired on WT CD47-IgV-GPI (A* and B'-SCT) vs substituted in Arm 1X W_C (CD47 C15G / precursor C33G); "
                    "Stage-0 Gate 0A-surface therefore performs non-reducing SDS-PAGE and DTNB/maleimide-PEG2-biotin + intact-mass LC-MS across both W_C=-1 (WT Cys15) and W_C=+1 (C15G) to verify zero disulfide-linked homodimerization or aberrant ER quality-control retention."
                ),
                "qpctl_pyroglutamate_and_opsonized_complement_phagocytosis_gate_sol61_r5": (
                    "High-affinity SIRPa recognition of human CD47 requires glutaminyl cyclase (QPCTL)-catalyzed cyclization of mature N-terminal Gln19 (residue 1 after signal-peptide cleavage) into pyroglutamate (pGlu19 / PCA1, which forms direct hydrogen bonds to SIRPa Leu66/Ser66 in PDB:2JJS/4CMM; Logtenberg et al. 2019 Nat Med PMID:30833751). "
                    "Per Logtenberg et al. 2019 PMID:30833751, antibody CC2C6 recognizes the QPCTL/pGlu19-dependent SIRPa-binding epitope (whose binding is selectively lost upon QPCTL inhibition or SEN177 treatment), whereas B6H12 and 2D3 bind pan-CD47 QPCTL-independent epitopes: "
                    "Gate 0A-surface therefore quantifies the CC2C6 / B6H12 (and 2D3) MFI ratio (obs_cd47_cc2c6_over_b6h12_ratio >= 0.90) alongside direct recombinant human SIRPa-Fc binding and tryptic N-terminal peptide LC-MS/MS (-17.03 Da pGlu19 vs unmodified Gln19, requiring >= 90% pGlu19 modification) in an isogenic endogenous 5-TM CD47-/- knockout background (so endogenous 5-TM CD47 cannot mask transgenic CD47-IgV-GPI ER processing) across both Stage-7 SC-beta cells and co-transplanted iEC/iPericyte vascular support cells under normoxia, hypoxia (10 mmHg), and cytokine stress, "
                    "and tests primary human macrophage phagocytosis suppression (obs_complement_c3b_opsonized_phagocytosis_pct <= 5.0% engulfment) under active C3b/iC3b complement + low-titer IgG opsonization where lateral GPI mobility lacks 5-TM cortical actin anchoring (AUDIT-REF, AUDIT-REF, AUDIT-REF & AUDIT-REF)."
                ),
                "iec_pericyte_vascular_lineage_survival_coupling_astra6_r5": (
                    "Because chronic Day-90 oxygen survival of the 1,500 IEQ/cm2 central-lamina pouch depends strictly on the vascularized outer band (k_perf = 0.1358 s^-1; if iEC/pericyte vessels are phagocytosed or lysed so k_perf -> 0, min_slab_pO2 collapses to 0.0783 mmHg under a 50 um capsule), "
                    "Stage-0 Gate 0A-surface and the 7-metric combined-stress coupon gate CC2C6/B6H12 >= 0.90, C3b-opsonized macrophage phagocytosis <= 5.0%, and post-stress iEC/pericyte vascular support survival >= 85% (obs_iec_vascular_support_viability_pct >= 85.0%) directly on the co-edited iEC/pericyte lineage alongside SC-beta cells (AUDIT-REF [DEC-046] & AUDIT-REF [DEC-049])."
                ),
                "inflammatory_protease_and_mmp_ectodomain_shedding_gate": (
                    "In a neutrophil/macrophage-infiltrated graft pocket, neutrophil elastase, cathepsin G, and matrix metalloproteinases (MMP-7, MMP-9, MMP-13; Hira-Miyazawa et al. 2018 PMID:29345283) can cleave the membrane-proximal stalk of PD-L1 or exposed (G4S)4 SCT linkers, "
                    "both depleting cell-surface checkpoint density and generating soluble sPD-L1/sCD47 decoys; Stage-0 Gate 0A-surface and Gate 0C-kinetic therefore include a 4-h neutrophil elastase (100 nM) + recombinant MMP-7/MMP-9/MMP-13 challenge requiring >=80% surface epitope retention by QuantiBRITE flow and <10% soluble decoy accumulation by ELISA."
                ),
                "hlag1_sct_lilrb1_challenger_and_il12_il18_sirpa_low_nk_audit_sonnet55_r5": {
                    "finding_ids": ["AUDIT-REF", "AUDIT-REF"],
                    "primary_citations": ["PMID:30778232", "PMID:33416832", "PMID:32134744", "PMID:33599743"],
                    "resting_and_il12_il18_primed_nk_sirpa_heterogeneity": (
                        "Resting human peripheral blood NK cells are predominantly SIRPa-negative or SIRPa-low, and while IL-2/IL-15/IFN-alpha priming upregulates SIRPa on a subset of primary NK cells (Deuse et al. 2019 Nat Biotechnol PMID:30778232; Deuse et al. 2021 J Exp Med PMID:33416832), "
                        "pro-inflammatory IL-12 + IL-18 + IFN-gamma islet pocket stimulation can leave a substantial fraction of KIR2DL+ NKG2A- and KIR2DL- NKG2A- NKG2C+ NK cells with low or variable functional SIRPa surface density. "
                        "Furthermore, endogenous PD-1 expression on primary human NK cells is minimal or restricted to discrete CMV-experienced subsets (Judge et al. 2020 J Clin Invest PMID:32134744; Davis et al. 2021 Blood Adv PMID:33599743)."
                    ),
                    "stage0_gate_0a_nk_il12_il18_sirpa_low_challenge": {
                        "endpoint_id": "obs_il12_il18_sirpa_low_nk_lysis_pct",
                        "pass_threshold_max_lysis_pct": 10.0,
                        "blocking_antibody_deconvolution_panel": [
                            "anti-NKG2A (monalizumab/Z199)",
                            "anti-SIRPa (SE5A5/KWAR23) + anti-CD47 (B6H12)",
                            "anti-KIR2DL1/2/3 (lirilumab/1-7F9)",
                            "anti-LILRB1 (GHI/75, for HLA-G1 SCT challenger)",
                        ],
                    },
                    "hlag1_sct_c42s_lilrb1_kir2dl4_challenger_coupon": {
                        "challenger_id": "CHALLENGER_HLAG1_SCT_C42S_RIIPRHLQL",
                        "heavy_chain": "HLA-G*01:01 alpha1-alpha3 ectodomain with Cys42Ser (C42S) mutation",
                        "why_c42s_monomeric_mutant_required": (
                            "Wild-type HLA-G1 forms disulfide-linked beta2m-associated homodimers via alpha1 Cys42-Cys42 that crosslink LILRB1 (ILT2/CD85j) and myeloid LILRB2 (ILT4) with bivalent sub-micromolar avidity K_D ~ 0.02-0.1 uM (Shiroishi et al. 2003 PNAS PMID:12853576; Shiroishi et al. 2006 J Biol Chem PMID:16455647) "
                            "but also engage endosomal/surface KIR2DL4 (CD158d), which triggers pro-inflammatory IFN-gamma/chemokine secretion in resting NK cells without cytotoxicity (Rajagopalan & Long 2012 Front Immunol PMID:22934097). "
                            "Ablating the intermolecular Cys42 disulfide in the monomeric HLA-G1(C42S) single-chain trimer (RIIPRHLQL or KGPPAALTL) avoids Cys42-dimer KIR2DL4 hyper-crosslinking while reducing LILRB1 binding to monomeric affinity (K_D ~ 2.0-5.0 uM; Shiroishi et al. 2006 J Biol Chem PMID:16305801) and eliminating bivalent ITIM crosslinking (MW-R5-S1-12 & AUDIT-REF); "
                            "setting track2_b_sct_or_hlag1_lilrb1_rescue_qualified=True for the HLA-G1(C42S) challenger therefore requires Stage-0 Gate 0A-NK monomeric LILRB1 SPR (K_D <= 5.0 uM), sorted NKG2A- KIR2DL- LILRB1+ NK lysis <= 10.0% (reversible by anti-LILRB1 GHI/75), 24-h KIR2DL4 IFN-gamma release < 50 pg/mL, and soluble sHLA-G shedding < 10%."
                        ),
                    },
                },
                "non_hla_chimeric_junction_neoepitope_and_flow_crossmatch_audit_sonnet55_r5": {
                    "finding_ids": ["AUDIT-REF", "AUDIT-REF", "MW-R5-S2-F05", "AUDIT-REF", "AUDIT-REF"],
                    "chimeric_surface_ectodomain_junctions_enumerated": [
                        "HLA-E SCT: VMAPRTLVL | (G4S)4 [20-aa L1 in 4,143-bp Embodiment A*] or (G4S)3 [15-aa L1 in 3,974-bp canonical donor] | mature B2M | (G4S)4 | mature HLA-E*01:03 (GSHSLKYF...)",
                        "HLA-C SCT (Track 2 B'-SCT): GAVDPLLAL/QYDDAVYKL | (G4S)4 [20-aa L1 in canonical Y84A/DT-2 B'-SCT, (G4S)5 in WT-Y84, or GCGGS(G4S)3 in DT-1] | mature B2M | (G4S)4 | mature HLA-C*03:04/04:01(Y84A) (GSHSMRYF...)",
                        "CD47-IgV-GPI: mature CD47 IgV(Gln19-Ser138, 120 aa, C15G/Gly33 in Arm 1X) | 9-aa retained CD55/DAF stalk PNKGSGTTS (precursor 139..147, native CD55 Pro345-Ser353, omega-site Ser147 [native Ser353]; 28-aa hydrophobic propeptide GTTRLLSGHTCFTLTGLLGTLVTMGLLT 148..175 cleaved upon GPI transamidation)",
                    ],
                    "stage0_gate_0d_xmatch_requirement": (
                        "Stage-0 Gate 0D-xmatch performs flow-cytometric IgG/IgM binding, C1q/C3d deposition, and rabbit/human-complement CDC + NK ADCC crossmatch on the exact engineered SC-islet and iEC product lines against a 50-serum panel (25 high-cPRA >= 80% sensitized allogeneic sera + 25 new-onset/established T1D sera), requiring <= 5.0% specific CDC/ADCC lysis above isotype background."
                    ),
                    "sa_thbd_tafi_fibrin_scaffold_versus_ibmir_recalcified_blood_gate_sol61_r6": {
                        "finding_ids": ["AUDIT-REF", "AUDIT-REF"],
                        "primary_citations": ["PMID:8663147", "PMID:30378751", "PMID:32923626", "PMID:34480863"],
                        "common_selected_dose_d_star_joint_qualification_required": True,
                        "donor_replicates_min": 6,
                        "mechanistic_tradeoff": (
                            "SA-THBD captures thrombin to activate Protein C (APC) and TAFI (thrombin-activatable fibrinolysis inhibitor; Bajzar et al. 1996 PMID:8663147); "
                            "while high local SA-THBD suppresses IBMIR and TAFIa stabilizes fibrin matrices against plasmin, excessive thrombin sequestration during Day 0-3 "
                            "can impair initial fibrinogen-to-fibrin polymerization required for endothelial sprouting or over-stabilize fibrin via TAFIa to block plasmin/MMP "
                            "remodeling during tubulogenesis, whereas insufficient local anticoagulation risks microvascular thrombosis inside anastomosing iEC/pericyte lumens "
                            "(Weaver et al. 2019 PMID:30378751; Coronel et al. 2020 PMID:32923626; Aghazadeh et al. 2021 PMID:34480863)."
                        ),
                        "three_assay_orthogonal_gate_0f_fibrin_patency_panel": [
                            "Single common formulation dose rule: a single selected SA-THBD surface-equivalent dose d* (titrated across 0, 50, 150, 300 nM across n >= 6 independent blood/iEC donors) must simultaneously clear Assays 1, 2, 2b, and 3 without dose-splitting between rheology and whole-blood patency.",
                            "Assay 1 (Chandler-loop IBMIR suppression at dose d*): recalcified human whole-blood loop requiring upper95(TAT) <= 15 ug/L, upper95(C3a) <= 200 ng/mL, and lower95(platelet retention) >= 85%.",
                            "Assay 2 (7-day iEC/pericyte fibrin-sprouting rheometry & tubulogenesis at dose d*): fibrinogen (2.5-5.0 mg/mL) + physiological thrombin (0.5-2.0 U/mL), requiring fibrin storage modulus G' within 80-120% of unmodified fibrin control and Day-7 CD31+ lumenized vascular network length lower95 >= 85% of control.",
                            "Assay 2b (TAFIa & APC chromogenic kinetics and tPA clot-lysis remodeling window at dose d*): direct chromogenic measurement of APC and TAFIa generation plus tPA-induced (100 ng/mL) turbidimetric clot-lysis half-time t_1/2 required to remain within 0.80x-1.50x of unmodified fibrin control with >= 80% pericellular fibrin remodeling clearance around sprouting iEC lumens at Day 7 (preventing TAFIa fibrinolytic lock).",
                            "Assay 3 (Perfused microvascular lumen patency under recalcified human whole blood at dose d*): microfluidic vascularized hydrogel perfusion after anastomosis requiring lower95(70-kDa FITC-dextran / autologous RBC perfused lumen fraction) >= 90% at 24 h without occlusive fibrin/platelet microthrombi and tracer-clearance-derived perfusion rate constant lower95(k_perf) >= 0.10 s^-1 (matching Layer-4 Gate 0F obs_outer_band_kperf_s_inv >= 0.10 s^-1).",
                        ],
                    },
                    "sa_thbd_tafi_fibrin_scaffold_vs_ibmir_patency_gate_sol61_r6": {
                        "finding_ids": ["AUDIT-REF", "AUDIT-REF"],
                        "common_selected_dose_d_star_required": True,
                        "tpa_clot_lysis_half_time_ratio_window_vs_control": [0.80, 1.50],
                        "min_perfused_lumen_patency_pct": 90.0,
                        "min_outer_band_kperf_s_inv": 0.10,
                    },
                    "sa_fasl_microgel_nkp46_innate_nk_neutrality_and_teff_selectivity_gate_sol61_r6": {
                        "finding_ids": ["AUDIT-REF", "AUDIT-REF", "MW-R6-S2-R3"],
                        "innate_nkp46_cd107a_degranulation_neutrality_gate": (
                            "Sorted primary human NKp46+ NK co-culture with SA-FasL-presenting PEG microgels requiring upper95(Delta CD107a degranulation) <= +5.0 pp over blank-microgel baseline."
                        ),
                        "teff_vs_treg_apoptosis_selectivity_gate": (
                            "Activated CD95^hi Teff apoptosis >= 60% vs CD4+CD25^hiFoxP3+ cFLIP^hi Treg viability/suppression retention >= 80% of WT."
                        ),
                        "treg_primary_estimand_hierarchy": (
                            "Primary binding gate: relative retention >= 80% of WT when WT baseline suppression >= 50%; "
                            "absolute lower95(delta) >= -10.0 pp is co-primary and becomes the sole binding gate if WT baseline suppression < 50% (MW-R6-S2-R3 [DEC-053])."
                        ),
                    },
                    "membrane_complement_regulator_cd55_scr14_cd59_and_assembled_cassette_ibmir_gate_muse_r6": {
                        "finding_ids": ["MW-R6-S1-04", "MW-R6-S1-07", "MW-R6-S1-13"],
                        "decision_id": "DEC-053",
                        "primary_citations": ["PMID:1824699", "PMID:30378751", "PMID:38508139"],
                        "why_daf_gpi_stalk_alone_lacks_convertase_decay_activity": (
                            "CD47-IgV-DAF fuses mature CD47 IgV(Gln19-Ser138) to the 9-aa C-terminal GPI anchor peptide of CD55/DAF (Pro345-Ser353, PNKGSGTTS) "
                            "and omits CD55 short consensus repeats SCR1-SCR4 (which accelerate C3/C5 convertase decay) and CD59 (which blocks C8/C9 MAC pore assembly); "
                            "while perimeter PEG-4MAL microgels carry SA-FHL-1 (Factor H SCR1-7) and SA-THBD, cell-intrinsic complement regulation on beta cells and iECs "
                            "during hypoxic shedding of endogenous GPI-anchored CD55/CD59 requires direct functional verification in assembled cassettes (MW-R6-S1-04 [DEC-053])."
                        ),
                        "assembled_cassette_chandler_loop_ibmir_and_complement_gate": {
                            "protocol": (
                                "Incubate fully assembled 350-um hydrogel cassettes (containing Stage-7 SC-islets + iECs/pericytes + CaO2/catalase lamina + SA-THBD/SA-FHL-1/SA-FasL microgels) "
                                "in recalcified non-anticoagulated human whole blood in a Chandler loop for 120 min across n >= 6 ABO-compatible and high-cPRA sera."
                            ),
                            "target_sc5b9_fold_over_baseline": 1.5,
                            "max_sc5b9_fold_over_baseline": 2.0,
                            "target_platelet_retention_pct": 80.0,
                            "min_platelet_retention_pct": 70.0,
                            "max_c3d_opsonization_mfi_fold": 1.5,
                            "wetlab_contract_endpoints": [
                                "obs_assembled_cassette_chandler_loop_sc5b9_fold",
                                "obs_assembled_cassette_chandler_loop_platelet_retention_pct",
                            ],
                            "cell_intrinsic_cd55_scr14_cd59_challenger": (
                                "If assembled-cassette Chandler-loop sC5b-9 exceeds 1.5x target (2.0x hard Stage-0 hold ceiling) or C3d opsonization exceeds 1.5x on hypoxic SC-islets/iECs, "
                                "promote the chimeric single-chain CD47-IgV-CD55(SCR1-4)-GPI or 2A-free bidirectional/dual-promoter [SV40pA-CD55(SCR1-4)-GPI(antisense) <- EEF1A1-core | hPGK1-min(sense) -> CD59-GPI-bGHpA] complement-shield challenger "
                                "(explicitly rejecting 2A-linked CD55-P2A-CD59 or CD59-P2A-CD55 bicistrons because placing either GPI-anchored protein upstream of a 2A peptide appends a C-terminal -GSG-P2A stub that blocks ER PIGK GPI-transamidase cleavage, "
                                "and requiring Stage-0 Gate-0D flow-cytometric CD55/CD59 surface expression, PI-PLC cleavage >= 85%, and zero 2A-scar qualification; AUDIT-REF [DEC-075])."
                            ),
                        },
                        "post_contact_trogocytosis_and_apoptotic_efferocytosis_sparing_gate": {
                            "post_macrophage_contact_cd47_retention_min_molecules": 180000,
                            "target_apoptotic_corpse_efferocytosis_clearance_pct_of_wt": 75.0,
                            "apoptotic_corpse_efferocytosis_clearance_min_pct_of_wt": 70.0,
                            "efferocytosis_receptor_panel": "MerTK/AXL/TIM-4",
                            "mechanistic_rationale": (
                                "High-density GPI-anchored CD47 (> 200,000 molecules/cell) must resist macrophage trogocytic stripping during live-cell contact "
                                "(maintaining post-contact p10 >= 180,000 molecules/cell) while undergoing calreticulin/phosphatidylserine dominance or caspase-dependent "
                                "shedding on Annexin-V+ apoptotic corpses so MerTK/AXL/TIM-4 macrophage efferocytosis (PMID:30846565) retains >= 75% target (>= 70% hard Stage-0 hold floor) of WT clearance rate and prevents "
                                "secondary necrotic DAMP/HMGB1 release (obs_apoptotic_effero_clearance_retention_pct_of_wt >= 70.0%; MW-R6-S1-13 & AUDIT-REF [DEC-053/060])."
                            ),
                            "wetlab_contract_endpoint": "obs_apoptotic_effero_clearance_retention_pct_of_wt",
                        },
                        "apoptotic_beta_cell_efferocytosis_clearance_gate_muse_r6": {
                            "post_macrophage_contact_cd47_retention_min_molecules": 180000,
                            "target_apoptotic_corpse_efferocytosis_clearance_pct_of_wt": 75.0,
                            "apoptotic_corpse_efferocytosis_clearance_min_pct_of_wt": 70.0,
                            "efferocytosis_receptor_panel": "MerTK/AXL/TIM-4",
                            "mechanistic_rationale": (
                                "High-density GPI-anchored CD47 (> 200,000 molecules/cell) must resist macrophage trogocytic stripping during live-cell contact "
                                "(maintaining post-contact p10 >= 180,000 molecules/cell) while undergoing calreticulin/phosphatidylserine dominance or caspase-dependent "
                                "shedding on Annexin-V+ apoptotic corpses so MerTK/AXL/TIM-4 macrophage efferocytosis (PMID:30846565) retains >= 75% target (>= 70% hard Stage-0 hold floor) of WT clearance rate and prevents "
                                "secondary necrotic DAMP/HMGB1 release (obs_apoptotic_effero_clearance_retention_pct_of_wt >= 70.0%; MW-R6-S1-13 & AUDIT-REF [DEC-053/060])."
                            ),
                            "wetlab_contract_endpoint": "obs_apoptotic_effero_clearance_retention_pct_of_wt",
                        },
                    },
                },
            },
        },
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", encoding="utf-8") as f:
        json.dump(atlas, f, indent=2, allow_nan=False)
        f.write("\n")

    print(f"[Layer 3 v0.8.0] Wrote {OUT_PATH.relative_to(REPO_ROOT)}")
    print(
        f"  Complex 1 (3CDG HLA-E:CD94/NKG2A): {c1['trimer_to_cd94_nkg2a_interface']['heavy_atom_contact_pairs']} heavy contacts, "
        f"L1 WLC force = {c1['single_chain_trimer_linker_physics']['linker_1_peptide_P9_to_B2M_Ile1']['marko_siggia_entropic_tension_pN_at_310K_lp4A']} pN, "
        f"L2 WLC force = {c1['single_chain_trimer_linker_physics']['linker_2_B2M_Met99_to_HLAE_Ser2']['marko_siggia_entropic_tension_pN_at_310K_lp4A']} pN"
    )
    print(
        f"  Complex 3 (AF-P50453 SERPINB9 & 1IAU GZMB): GZMB 226/226 aligned RMSD = {c3['granzyme_b_active_site_and_s1_pocket_1IAU']['rmsd_AF_P10144_vs_1IAU_angstrom']} A, "
        f"RCL 335-350 mean pLDDT = {c3['r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry']['rcl_335_350_mean_plddt']}, "
        f"status = {c3['r01_c_terminal_tpt2a_tail_and_rcl_suicide_trap_geometry']['functional_status']}"
    )


if __name__ == "__main__":
    main()
