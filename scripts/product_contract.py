"""Versioned exact-product identifiers shared by the ranking and scenario planes."""
from __future__ import annotations

from copy import deepcopy
import json
import re
from pathlib import Path

CONTRACT_PATH = Path("data/protocols/islet1_product_definitions_v0.8.0.json")


def load_product_contracts(root: Path | str | None = None) -> dict:
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    payload = json.loads((root / CONTRACT_PATH).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("products"), list):
        raise ValueError("Product registry must contain a products list")
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
    if not (product["product_variant_id"].startswith("ISLET1_" + product["family"]) or (product["family"] == "A" and product["product_variant_id"].startswith("ISLET1_COMP_"))):
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
    expected_sequence_status = "REFERENCE_DONORS_AVAILABLE_REQUIRES_EDITED_ALLELE_QC" if product["family"] == "A" else "NOT_MATERIALIZED_SEPARATE_DONOR_LEDGER_REQUIRED"
    if product["sequence_materialization_status"] != expected_sequence_status:
        raise ValueError("Sequence materialization status conflicts with the product family")
    if product["physical_qualification_status"] != "PREREGISTERED_PENDING":
        raise ValueError("This release contains no physically qualified exact product")
    if product["family"] == "B" and product["execution_supported"]:
        raise ValueError("Unmaterialized Scenario B designs must not be executable products")
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
    for product in load_product_contracts(root)["products"]:
        if product_id == product["product_variant_id"] or product_id in product.get("aliases", []):
            return deepcopy(product)
    raise ValueError(f"Unknown product variant: {product_id}")
