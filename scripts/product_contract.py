"""Versioned exact-product identifiers shared by the ranking and scenario planes."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import re
from pathlib import Path

CONTRACT_PATH = Path("data/protocols/islet1_product_definitions_v0.8.0.json")
DESIGN_B_FASTA_PATH = Path("data/constructs/islet1_design_b_complete_constructs_v0.8.0.fasta")
DESIGN_B_FASTA_SHA256 = "4bb20389773c76b2705031ba675487825dfbd94d920a692cde75b7d56af22aa4"
DESIGN_B_C1_ID = "DESIGNB_L5A_CXCL10_BFB_SYN6_ACC_4078bp"
DESIGN_B_C1_SHA256 = "3392af83d0b37b63c7db2094e88ae3b7034811e3a7b64f678f45ffd56194bf56"
DESIGN_B_C5_ID = "DESIGNB_L6_AAVS1_CD47_ICASP9_2AFREE_HUMAN_ACC_4172bp"
DESIGN_B_C5_SHA256 = "ba264b14b5cf81b191a8e2126cf5f2c9949559a3782f6ff69bb206690d705882"
DESIGN_B_BOUND_C1_C5_PAIR_SHA256 = "82c351af082e0f97361e715c36e2ab5872b8d097efa1c9f59dfda53ea60496ea"
STAGED_DESIGN_B2_PRODUCT_IDS = frozenset({
    "ISLET1_B2_5LOCUS_BASE",
    "ISLET1_B2_5LOCUS_TYK2_PE",
})
ALLOWED_DESIGN_B_LIFECYCLE_STATES = frozenset({
    "STAGED_REFERENCE_FASTA_ARTIFACT",
    "REGISTRY_DONOR_REFERENCE_BOUND",
    "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC",
    "EXECUTION_SUPPORTED_QUALIFIED_PRODUCT",
})


def _parse_fasta_map(fasta_text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    cur_id: str | None = None
    buf: list[str] = []
    for raw in fasta_text.splitlines():
        ln = raw.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            if cur_id is not None:
                records[cur_id] = "".join(buf).upper()
            cur_id = ln[1:].split()[0]
            buf = []
        else:
            buf.append(ln)
    if cur_id is not None:
        records[cur_id] = "".join(buf).upper()
    return records


def load_product_contracts(root: Path | str | None = None) -> dict:
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    payload = json.loads((root / CONTRACT_PATH).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("products"), list):
        raise ValueError("Product registry must contain a products list")
    if "mcb_donor_blocking_policy" in payload:
        pol = payload["mcb_donor_blocking_policy"]
        if (
            not isinstance(pol, dict)
            or type(pol.get("min_independent_abo_oo_seedlines")) is not int
            or pol["min_independent_abo_oo_seedlines"] < 2
            or not isinstance(pol.get("kinship_pi_hat_lt"), (int, float))
            or isinstance(pol.get("kinship_pi_hat_lt"), bool)
            or not (0.0 < float(pol["kinship_pi_hat_lt"]) <= 0.05)
            or type(pol.get("donor1_paired_lot_blocks_min")) is not int
            or pol["donor1_paired_lot_blocks_min"] < 6
            or type(pol.get("confirmatory_welch_tost_lots_per_bank")) is not int
            or pol["confirmatory_welch_tost_lots_per_bank"] < 6
        ):
            raise ValueError("Invalid mcb_donor_blocking_policy in product registry (IIP-0024; DEC-103)")
    if "design_b_four_state_lifecycle_policy" in payload:
        lc_pol = payload["design_b_four_state_lifecycle_policy"]
        states = lc_pol.get("states") if isinstance(lc_pol, dict) else None
        required_state_keys = {
            "STATE_1_STAGED_REFERENCE_FASTA_ARTIFACT",
            "STATE_2_REGISTRY_DONOR_REFERENCE_BOUND",
            "STATE_3_STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC",
            "STATE_4_EXECUTION_SUPPORTED_QUALIFIED_PRODUCT",
        }
        if (
            not isinstance(states, dict)
            or set(states.keys()) != required_state_keys
            or states["STATE_1_STAGED_REFERENCE_FASTA_ARTIFACT"].get("fasta_sha256") != DESIGN_B_FASTA_SHA256
            or states["STATE_2_REGISTRY_DONOR_REFERENCE_BOUND"].get("bound_c1_c5_pair_sha256") != DESIGN_B_BOUND_C1_C5_PAIR_SHA256
            or states["STATE_3_STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"].get("execution_supported") is not False
            or states["STATE_4_EXECUTION_SUPPORTED_QUALIFIED_PRODUCT"].get("currently_achieved_family_b_count") != 0
        ):
            raise ValueError("Invalid design_b_four_state_lifecycle_policy in product registry (IIP-0026; DEC-105)")
    fasta_file = root / DESIGN_B_FASTA_PATH
    if fasta_file.is_file():
        raw_bytes = fasta_file.read_bytes()
        if hashlib.sha256(raw_bytes).hexdigest() != DESIGN_B_FASTA_SHA256:
            raise ValueError("Staged Design B FASTA SHA-256 mismatch against registry pin (IIP-0026; DEC-105)")
        seq_map = _parse_fasta_map(raw_bytes.decode("utf-8"))
        c1_seq = seq_map.get(DESIGN_B_C1_ID, "")
        c5_seq = seq_map.get(DESIGN_B_C5_ID, "")
        pair_sha = hashlib.sha256((c1_seq + "|" + c5_seq).encode("utf-8")).hexdigest()
        if (
            hashlib.sha256(c1_seq.encode("utf-8")).hexdigest() != DESIGN_B_C1_SHA256
            or hashlib.sha256(c5_seq.encode("utf-8")).hexdigest() != DESIGN_B_C5_SHA256
            or pair_sha != DESIGN_B_BOUND_C1_C5_PAIR_SHA256
        ):
            raise ValueError("Staged Design B c1+c5 pair digest mismatch against registry pin (IIP-0026; DEC-105)")
    products = payload["products"]
    for product in products:
        validate_product_definition(product)
    if not products or len({p["product_variant_id"] for p in products}) != len(products):
        raise ValueError("Product variants must be nonempty and unique")
    identifiers = set()
    for product in products:
        validate_product_definition(product)
        for identifier in [product["product_variant_id"], *product.get("aliases", [])]:
            if identifier in identifiers:
                raise ValueError("Product identifiers and aliases must be globally unique")
            identifiers.add(identifier)
    return payload


def validate_product_definition(product: dict) -> dict:
    if not isinstance(product, dict):
        raise TypeError("Product definition must be a dictionary")
    required = {"product_variant_id", "family", "dsb_or_indel_loci", "base_or_prime_edit_loci",
                "alleles", "cassette_topology", "K", "B", "sequence_materialization_status",
                "physical_qualification_status", "execution_supported", "donor_lineage", "hla_retention", "pending_assay_gates", "physical_manufacturing_yield"}
    if required - product.keys():
        raise ValueError(f"Incomplete product definition: {sorted(required - product.keys())}")
    if not isinstance(product["product_variant_id"], str) or not re.fullmatch(r"ISLET1_[A-Z0-9_]+", product["product_variant_id"]):
        raise ValueError("Invalid canonical product identifier")
    if product["family"] not in {"A", "B"}:
        raise ValueError("Unknown product family")
    if not (
        product["product_variant_id"].startswith("ISLET1_" + product["family"])
        or (product["family"] == "A" and product["product_variant_id"].startswith("ISLET1_COMP_"))
        or (
            product["family"] == "A"
            and product["product_variant_id"].startswith("ISLET1_Z_")
            and product["execution_supported"] is False
        )
    ):
        raise ValueError("Product identifier conflicts with its family")
    if product["physical_manufacturing_yield"] is not None:
        raise ValueError("Physical manufacturing yield remains unmeasured in this release")
    if type(product["execution_supported"]) is not bool:
        raise ValueError("execution_supported must be boolean")
    aliases = product.get("aliases", [])
    if not isinstance(aliases, list) or any(not isinstance(a, str) or not re.fullmatch(r"[A-Z0-9_]+", a) for a in aliases):
        raise ValueError("Aliases must be bounded identifiers")
    if len(set(aliases)) != len(aliases) or product["product_variant_id"] in aliases:
        raise ValueError("Product aliases must be unique and distinct from the canonical identifier")
    for key in ("alleles", "cassette_topology"):
        value = product[key]
        if not isinstance(value, dict) or not value or any(not isinstance(k, str) or not k.strip() or not isinstance(v, str) or not v.strip() for k,v in value.items()):
            raise ValueError(f"{key} must be a nonempty explicit mapping")
    for key in ("donor_lineage", "hla_retention"):
        if not isinstance(product[key], dict) or not product[key]:
            raise ValueError(f"{key} must explicitly define donor/HLA scope")
    if "donor_seedline_id_required" in product["donor_lineage"] and type(product["donor_lineage"]["donor_seedline_id_required"]) is not bool:
        raise ValueError("donor_seedline_id_required must be boolean (IIP-0024; DEC-103)")
    if "physical_route_id" in product and (not isinstance(product["physical_route_id"], str) or not re.fullmatch(r"[A-Z0-9_]+", product["physical_route_id"])):
        raise ValueError("physical_route_id must be a bounded uppercase identifier (IIP-0024; DEC-103)")
    gates = product["pending_assay_gates"]
    if not isinstance(gates, list) or not gates or any(not isinstance(g, str) or not g.strip() for g in gates):
        raise ValueError("Explicit pending assay gates are required")
    known_loci = {"B2M", "CIITA", "CXCL10", "AAVS1", "HLA-A", "HLA-B", "TYK2", "INS"}
    for loci_key, count_key in [("dsb_or_indel_loci", "K"), ("base_or_prime_edit_loci", "B")]:
        loci = product[loci_key]
        if not isinstance(loci, list) or any(not isinstance(x, str) or x not in known_loci for x in loci) or len(loci) != len(set(loci)):
            raise ValueError(f"{loci_key} must list unique loci")
        if type(product[count_key]) is not int or product[count_key] != len(loci):
            raise ValueError(f"{count_key} does not match the explicit locus list")
    if not (set(product["dsb_or_indel_loci"]) | set(product["base_or_prime_edit_loci"])) <= product["alleles"].keys():
        raise ValueError("Every edited locus needs an explicit allele definition")
    # Enforce 4-state Design B lifecycle while preserving Family A reference-compiler semantics (IIP-0026; DEC-105):
    if product["family"] == "A":
        if product["sequence_materialization_status"] != "REFERENCE_DONORS_AVAILABLE_REQUIRES_EDITED_ALLELE_QC":
            raise ValueError("Sequence materialization status conflicts with the product family")
    elif product["product_variant_id"] in STAGED_DESIGN_B2_PRODUCT_IDS:
        if product["sequence_materialization_status"] != "STAGED_RESEARCH_COMPARATOR_DONORS_BOUND_PENDING_ALLELE_QC":
            raise ValueError("Staged Design B2 variants require STAGED_RESEARCH_COMPARATOR_DONORS_BOUND_PENDING_ALLELE_QC (IIP-0026; DEC-105)")
        if product.get("lifecycle_state") != "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC":
            raise ValueError("Staged Design B2 variants must remain at State 3 STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC (IIP-0026; DEC-105)")
        binding = product.get("staged_design_b_donor_reference_binding")
        if (
            not isinstance(binding, dict)
            or binding.get("fasta_path") != str(DESIGN_B_FASTA_PATH)
            or binding.get("fasta_sha256") != DESIGN_B_FASTA_SHA256
            or binding.get("cxcl10_locus5_primary_donor_id") != DESIGN_B_C1_ID
            or binding.get("cxcl10_locus5_primary_donor_sha256") != DESIGN_B_C1_SHA256
            or binding.get("aavs1_locus6_primary_donor_id") != DESIGN_B_C5_ID
            or binding.get("aavs1_locus6_primary_donor_sha256") != DESIGN_B_C5_SHA256
            or binding.get("bound_c1_c5_pair_sha256") != DESIGN_B_BOUND_C1_C5_PAIR_SHA256
        ):
            raise ValueError("Invalid or mismatched staged_design_b_donor_reference_binding (IIP-0026; DEC-105)")
    else:
        if product["sequence_materialization_status"] != "NOT_MATERIALIZED_SEPARATE_DONOR_LEDGER_REQUIRED":
            raise ValueError("Sequence materialization status conflicts with the product family")
    if "lifecycle_state" in product and product["lifecycle_state"] not in ALLOWED_DESIGN_B_LIFECYCLE_STATES:
        raise ValueError("Invalid lifecycle_state in product definition (IIP-0026; DEC-105)")
    if product["physical_qualification_status"] != "PREREGISTERED_PENDING":
        raise ValueError("This release contains no physically qualified exact product")
    if product["family"] == "B" and product["execution_supported"]:
        raise ValueError("Scenario B research comparators must keep execution_supported=False pending State-4 physical qualification")
    if "locus1_donor_variants_metadata" in product:
        l1_meta = product["locus1_donor_variants_metadata"]
        if not isinstance(l1_meta, dict) or not l1_meta:
            raise ValueError("locus1_donor_variants_metadata must be a non-empty dictionary")
        expected_variants = {
            "canonical_arm1b_3974bp": (
                3974,
                4264,
                "f74728767987c18bdcb8c9d49773ee8dbe10e6006e57aa622638e4fb07062b0b",
            ),
            "embodiment_a_star_4143bp": (
                4143,
                4433,
                "0995f8e35791d219dcc905df4bd50145c9b2a43ddb84efdb502b41dc1f2ae6a8",
            ),
        }
        for v_key, (exp_len, exp_pkg_len, exp_sha) in expected_variants.items():
            if v_key not in l1_meta or not isinstance(l1_meta[v_key], dict):
                raise ValueError(f"Missing required locus1 donor variant metadata: {v_key}")
            v_rec = l1_meta[v_key]
            if v_rec.get("donor_length_bp") != exp_len or v_rec.get("packaged_aav6_length_bp") != exp_pkg_len:
                raise ValueError(f"Invalid donor/packaged length in {v_key}")
            sha_val = v_rec.get("donor_sha256")
            if not isinstance(sha_val, str) or not re.fullmatch(r"[0-9a-f]{64}", sha_val) or sha_val != exp_sha:
                raise ValueError(f"Invalid or mismatched donor_sha256 in {v_key}: {sha_val}")
    return product


def get_product_definition(product_id: str, root: Path | str | None = None) -> dict:
    if not isinstance(product_id, str) or not product_id or product_id != product_id.strip():
        raise ValueError(f"Unknown product variant: {product_id!r}")
    for product in load_product_contracts(root)["products"]:
        if product_id == product["product_variant_id"] or product_id in product.get("aliases", []):
            return deepcopy(product)
    raise ValueError(f"Unknown product variant: {product_id}")
