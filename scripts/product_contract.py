"""Versioned exact-product identifiers shared by the ranking and scenario planes."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import re
from pathlib import Path

CONTRACT_PATH = Path("data/protocols/islet1_product_definitions_v0.8.0.json")
DESIGN_B_FASTA_PATH = Path("data/constructs/islet1_design_b_complete_constructs_v0.8.0.fasta")
DESIGN_B_FASTA_SHA256 = "4bb20389773c76b2705031ba675487825dfbd94d920a692cde75b7d56af22aa4"
DESIGN_B_C1_ID = "DESIGNB_L5A_CXCL10_BFB_SYN6_ACC_4078bp"
DESIGN_B_C1_LENGTH_BP = 4078
DESIGN_B_C1_SHA256 = "3392af83d0b37b63c7db2094e88ae3b7034811e3a7b64f678f45ffd56194bf56"
DESIGN_B_C5_ID = "DESIGNB_L6_AAVS1_CD47_ICASP9_2AFREE_HUMAN_ACC_4172bp"
DESIGN_B_C5_LENGTH_BP = 4172
DESIGN_B_C5_SHA256 = "ba264b14b5cf81b191a8e2126cf5f2c9949559a3782f6ff69bb206690d705882"
DESIGN_B_BOUND_C1_C5_PAIR_SHA256 = "82c351af082e0f97361e715c36e2ab5872b8d097efa1c9f59dfda53ea60496ea"
STAGED_DESIGN_B2_PRODUCT_IDS = frozenset({
    "ISLET1_B2_5LOCUS_BASE",
    "ISLET1_B2_5LOCUS_TYK2_PE",
})
LOCUS2_V3_FASTA_PATH = Path("data/constructs/islet1_locus2_ciita_modeb_comparators_v0.8.0.fasta")
LOCUS2_V3_FASTA_SHA256 = "25ffcb10e2bad77a4b34725c7eff0532990ebf7b80864cf7217d33f5ccbe5b6d"
RULE_S11_GATE_ID = "S11_Locus2_CIITA_ModeB_Vascular_Complement_and_Efferocytosis_Gate"
LOCUS2_V3_SCRAMBLE_CONTROL_ID = "L2_V3_CTRL_NT_SCRAMBLE"
LOCUS2_V3_SCRAMBLE_CONTROL_SHA256 = "d1c8d1bbe01e89dac2256d48b51d2dac03b7cba41d076a0e9bbed6cc249c43ff"
LOCUS2_V3_PINNED_B2MKO_INTRONLESS_CONTROL_IDS: dict[str, str] = {
    "CD200": "L2_V3_B2MKO_INTRONLESS_CD200_CTRL",
    "MANF": "L2_V3_B2MKO_INTRONLESS_MANF_CTRL",
    "CCL22": "L2_V3_B2MKO_INTRONLESS_CCL22_CTRL",
}
LOCUS2_V3_BARRIER_MATCHED_80BP_COMPARATOR_ID = "L2_V3_CO_PRIMARY_CD200"
LOCUS2_V3_DONOR_BINDINGS: dict[str, dict[str, object]] = {
    "ISLET1_A4_L2_V3_CO_PRIMARY_CD200": {
        "locus2_donor_id": "L2_V3_CO_PRIMARY_CD200",
        "locus2_donor_sha256": "723e052172adfe5ae6c74e35e44cf727202d2716010772041d8c12abd76529c8",
        "donor_length_bp": 4012,
        "nominal_packaged_aav6_bp": 4302,
        "family": "A",
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "has_intronic_mir": True,
        "mir_class": "ACTIVE_CD276_PVR",
        "unit2_cargo": "CD200",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
    "ISLET1_A4_L2_V3_CO_PRIMARY_MANF": {
        "locus2_donor_id": "L2_V3_CO_PRIMARY_MANF",
        "locus2_donor_sha256": "e4d9ba0c1019ecc04930f38dd484d433f3a476e0af886765fa4779a3225a014a",
        "donor_length_bp": 3724,
        "nominal_packaged_aav6_bp": 4014,
        "family": "A",
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "has_intronic_mir": True,
        "mir_class": "ACTIVE_CD276_PVR",
        "unit2_cargo": "MANF",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
    "ISLET1_B2_L2_V3_DESIGN_B_CD200": {
        "locus2_donor_id": "L2_V3_DESIGN_B_CD200",
        "locus2_donor_sha256": "4b2b5713e29a190bcca115507a876fd48ce9ef92b0e6f793a08cd925a33334e1",
        "donor_length_bp": 3692,
        "nominal_packaged_aav6_bp": 3982,
        "family": "B",
        "background_id": "DESIGN_B_ENDOGENOUS_HLAC_INTACT",
        "has_intronic_mir": False,
        "mir_class": "INTRONLESS_NO_MIR",
        "unit2_cargo": "CD200",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
    "ISLET1_B2_L2_V3_DESIGN_B_MANF": {
        "locus2_donor_id": "L2_V3_DESIGN_B_MANF",
        "locus2_donor_sha256": "841a52222c703426c44d4f18fe799b646db71479b008d4b8fdbdcb7c2c2038c8",
        "donor_length_bp": 3404,
        "nominal_packaged_aav6_bp": 3694,
        "family": "B",
        "background_id": "DESIGN_B_ENDOGENOUS_HLAC_INTACT",
        "has_intronic_mir": False,
        "mir_class": "INTRONLESS_NO_MIR",
        "unit2_cargo": "MANF",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
    "ISLET1_A4_L2_V3_ALT2_CCL22": {
        "locus2_donor_id": "L2_V3_ALT2_CCL22",
        "locus2_donor_sha256": "ba84d57577ff5143b3a2b813c4b12aab5139570eb60903a3c7aebae6fdaad30e",
        "donor_length_bp": 3457,
        "nominal_packaged_aav6_bp": 3747,
        "family": "A",
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "has_intronic_mir": True,
        "mir_class": "ACTIVE_CD276_PVR",
        "unit2_cargo": "CCL22",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
    "ISLET1_A4_L2_V3_CTRL_NT_SCRAMBLE": {
        "locus2_donor_id": "L2_V3_CTRL_NT_SCRAMBLE",
        "locus2_donor_sha256": "d1c8d1bbe01e89dac2256d48b51d2dac03b7cba41d076a0e9bbed6cc249c43ff",
        "donor_length_bp": 4012,
        "nominal_packaged_aav6_bp": 4302,
        "family": "A",
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "has_intronic_mir": True,
        "mir_class": "NON_TARGETING_SCRAMBLE",
        "unit2_cargo": "CD200",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
    "ISLET1_A4_L2_V3_BARRIER_CHS4_250": {
        "locus2_donor_id": "L2_V3_BARRIER_CHS4_250",
        "locus2_donor_sha256": "f58995d9396ed09432e14a2d231a995aecd76122b852ee07ddd629cb07b6af66",
        "donor_length_bp": 4182,
        "nominal_packaged_aav6_bp": 4472,
        "family": "A",
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "has_intronic_mir": True,
        "mir_class": "ACTIVE_CD276_PVR",
        "unit2_cargo": "CD200",
        "unit2_promoter": "Syn_5xISRE_minP",
        "barrier_variant_type": "cHS4_250bp_C79T_C80A_COMBINED_LENGTH_AND_SEQUENCE",
    },
    "ISLET1_A4_L2_V3_CD274P_ZERO_TATATA": {
        "locus2_donor_id": "L2_V3_CD274P_ZERO_TATATA",
        "locus2_donor_sha256": "ecd6df0d48f0fc6f897482d14e4ec2c92e1f4740e58e11b77c40fe9760de5577",
        "donor_length_bp": 4114,
        "nominal_packaged_aav6_bp": 4404,
        "family": "A",
        "background_id": "DESIGN_A_OR_TRACK1_X_B2M_NULL",
        "has_intronic_mir": True,
        "mir_class": "ACTIVE_CD276_PVR",
        "unit2_cargo": "CD200",
        "unit2_promoter": "human_CD274_prom_250bp",
        "barrier_variant_type": "cHS4_1xFII_80bp",
    },
}
STAGED_LOCUS2_V3_PRODUCT_IDS = frozenset(LOCUS2_V3_DONOR_BINDINGS.keys())
STAGED_LOCUS2_V3_DONOR_ALIASES = frozenset(str(v["locus2_donor_id"]) for v in LOCUS2_V3_DONOR_BINDINGS.values())
DESIGN_Y_V2_FASTA_PATH = Path("data/constructs/islet1_design_y_fullstack_v2_v0.8.0.fasta")
DESIGN_Y_V2_FASTA_SHA256 = "9826f7917c73f56b62f38db9df0fd61d2ffbf90da3f454392c9f415c16de336c"
DESIGN_Y_V2_DONOR_BINDINGS: dict[str, dict[str, object]] = {
    "ISLET1_Y_FULLSTACK_V2_B2M_11TU": {
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "fasta_record_id": "DESIGN_Y_FULLSTACK_v2_B2M_11TU",
        "donor_sha256": "e77a358ca5bc55c9e0a3070ab46b115c3314660d3a9be6221ef9e91d698cab4e",
        "donor_length_bp": 17237,
        "tu_count": 11,
        "pigk_ectopic_gpi_count": 1,
        "encoded_2a_motifs_6frame": 0,
        "role": "FULL_CHR15_11TU_REFERENCE",
        "stage0_promotion_gates": "Y-01..Y-09 + R3P-01..R3P-06 (primary 11-TU CD46_WT full-stack candidate)",
        "full_record_sha256": "3aa9ed8d81156d22c293dfc849c6ef246aa9589a2d5dddb824668a46e153a584",
    },
    "ISLET1_Y_FULLSTACK_V2_B2M_HALF_A_ATTB_GT": {
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "fasta_record_id": "DESIGN_Y_FULLSTACK_v2_B2M_HALF_A_ATTB_GT",
        "donor_sha256": "bd8dd98ae5a81f1ac1157523ac7195fe409af34e52ae310528071c896e0708df",
        "donor_length_bp": 10450,
        "tu_count": 7,
        "pigk_ectopic_gpi_count": 1,
        "encoded_2a_motifs_6frame": 0,
        "role": "SEQUENTIAL_STEP1_ATTB_GT_DONOR",
        "stage0_promotion_gates": "COMPONENT_ONLY_EXCLUDED_FROM_FULLSTACK_CLEARANCE (Y-05 Step-1 attP-GT knock-in and attP-CA landing-site fidelity only)",
        "full_record_sha256": "372af0bc127fd9ffc5cbc2753d8e3d14d0493bdf7e660ddcde3674edebf8ba7a",
    },
    "ISLET1_Y_FULLSTACK_V2_B2M_HALF_B_ATTB_CA": {
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "fasta_record_id": "DESIGN_Y_FULLSTACK_v2_B2M_HALF_B_ATTB_CA",
        "donor_sha256": "097f8b46bd8be54f6305372ddcba2b7e863aa62c3736b912808bd2ec3206fa2b",
        "donor_length_bp": 6909,
        "tu_count": 4,
        "pigk_ectopic_gpi_count": 0,
        "encoded_2a_motifs_6frame": 0,
        "role": "SEQUENTIAL_STEP2_ATTB_CA_DONOR",
        "stage0_promotion_gates": "COMPONENT_ONLY_EXCLUDED_FROM_FULLSTACK_CLEARANCE (Y-05 Step-2 attP-CA knock-in into HALF_A clone only)",
        "full_record_sha256": "61587ac04bf28422c0aa8dd72be1680b4c4a5c7d6f6e5052727566155e4b2613",
    },
    "ISLET1_Y_FULLSTACK_V2_DEC053_CLYBL_CDK1_HTMPK": {
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "fasta_record_id": "DESIGN_Y_FULLSTACK_v2_DEC053_CHR13_CLYBL_CDK1_HTMPK",
        "donor_sha256": "f631592ed4c08fd433a03d2b085c04196b4b5b812b6f7fd49916c06e865d1323",
        "donor_length_bp": 2297,
        "tu_count": 1,
        "pigk_ectopic_gpi_count": 0,
        "encoded_2a_motifs_6frame": 0,
        "role": "DEC053_CHR13_CLYBL_SAFETY_SWITCH",
        "stage0_promotion_gates": "COMPONENT_ONLY_EXCLUDED_FROM_FULLSTACK_CLEARANCE (Y-09 independent chr13 CLYBL::CDK1p::hTMPK-F105Y module only)",
        "full_record_sha256": "ee224e637491625ded5434cc4f50628b2dabcaa111f6c8e523b0d5378d1463f7",
    },
    "ISLET1_Y_FULLSTACK_V2_CD55_CFH_COMPARATOR": {
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "fasta_record_id": "DESIGN_Y_FULLSTACK_v2_CD55_CFH_COMPARATOR",
        "donor_sha256": "d7291ca6e610f8251ab04bbf6fc62d7747e7baaf0e680300f46b895920aab9d8",
        "donor_length_bp": 18760,
        "tu_count": 12,
        "pigk_ectopic_gpi_count": 3,
        "encoded_2a_motifs_6frame": 0,
        "role": "GATE0V_CD55_CFH_BACKUP_COMPARATOR",
        "stage0_promotion_gates": "Y-01..Y-09 + R3P-01..R3P-06 (12-TU CD55+CFH Gate 0V backup full-stack candidate; activated only if primary CD46_WT fails Gate 0V MCB viral clearance)",
        "full_record_sha256": "ffeb49b2055cda9039ddf514ee87403b2b1cf7d91adee3d151448c918d43bb33",
    },
    "ISLET1_Y_FULLSTACK_V2_CTRL_EQUAL_SPACING": {
        "proposal_id": "IIP-0034",
        "decision_id": "DEC-113",
        "fasta_record_id": "DESIGN_Y_FULLSTACK_v2_CTRL_EQUAL_SPACING",
        "donor_sha256": "d5f69bd613cc0d598b846ba061154dc522b8ec6fb064610cfad01562505d2925",
        "donor_length_bp": 17237,
        "tu_count": 11,
        "pigk_ectopic_gpi_count": 1,
        "encoded_2a_motifs_6frame": 0,
        "role": "EQUAL_SPACING_ATTRIBUTION_CONTROL",
        "stage0_promotion_gates": "CONTROL_ONLY_EXCLUDED_FROM_FULLSTACK_CLEARANCE (pooled 5-TU v2-vs-v1 spacing/retention, Y-05 conditional backup CD59/CD47-IgV-GPI surface-density ratio, and Y-07 paired 2-pulse GSIS comparator only)",
        "full_record_sha256": "bddf38d1aed8747f45e6ad112a02d4baced270418ec3f6d135dfb632958a4e99",
    },
}
STAGED_DESIGN_Y_V2_PRODUCT_IDS = frozenset(DESIGN_Y_V2_DONOR_BINDINGS.keys())


def _reject_nonfinite_product_json(token: str) -> Any:
    raise ValueError(f"Non-standard JSON constant rejected in product registry: {token!r}")


def _parse_finite_product_json_float(token: str) -> float:
    val = float(token)
    if not math.isfinite(val):
        raise ValueError(f"Non-finite JSON float overflow rejected in product registry: {token!r}")
    return val
ALLOWED_DESIGN_B_LIFECYCLE_STATES = frozenset({
    "STAGED_REFERENCE_FASTA_ARTIFACT",
    "REGISTRY_DONOR_REFERENCE_BOUND",
    "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC",
    "EXECUTION_SUPPORTED_QUALIFIED_PRODUCT",
})


def _reject_duplicate_json_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    out: dict[str, object] = {}
    for key, value in pairs:
        if key in out:
            raise ValueError(f"Duplicate JSON key in product registry: {key!r}")
        out[key] = value
    return out


def _strict_json_equal(a: object, b: object) -> bool:
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_strict_json_equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(_strict_json_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return math.isnan(a) == math.isnan(b) and a == b
    return a == b


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
                raw_seq = "".join(buf)
                if not raw_seq or not all(ch in "ACGTacgt" for ch in raw_seq):
                    raise ValueError(f"Empty or non-ACGT sequence in FASTA record {cur_id!r}")
                records[cur_id] = raw_seq.upper()
            tokens = ln[1:].split()
            if not tokens:
                raise ValueError("Empty FASTA record header")
            cur_id = tokens[0]
            if cur_id in records:
                raise ValueError(f"Duplicate FASTA record identifier: {cur_id!r}")
            buf = []
        else:
            if cur_id is None:
                raise ValueError("FASTA sequence data encountered before first header")
            buf.append(ln)
    if cur_id is not None:
        raw_seq = "".join(buf)
        if not raw_seq or not all(ch in "ACGTacgt" for ch in raw_seq):
            raise ValueError(f"Empty or non-ACGT sequence in FASTA record {cur_id!r}")
        records[cur_id] = raw_seq.upper()
    if not records:
        raise ValueError("Empty FASTA file")
    return records


def _validate_staged_design_b_binding(binding: object, product_variant_id: str) -> None:
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
        raise ValueError(
            f"Invalid or mismatched staged Design B donor binding for {product_variant_id} (IIP-0026; DEC-105)"
        )
    if "cxcl10_locus5_primary_donor_length_bp" in binding and (
        type(binding["cxcl10_locus5_primary_donor_length_bp"]) is not int
        or binding["cxcl10_locus5_primary_donor_length_bp"] != DESIGN_B_C1_LENGTH_BP
    ):
        raise ValueError(f"Invalid cxcl10_locus5_primary_donor_length_bp for {product_variant_id}")
    if "aavs1_locus6_primary_donor_length_bp" in binding and (
        type(binding["aavs1_locus6_primary_donor_length_bp"]) is not int
        or binding["aavs1_locus6_primary_donor_length_bp"] != DESIGN_B_C5_LENGTH_BP
    ):
        raise ValueError(f"Invalid aavs1_locus6_primary_donor_length_bp for {product_variant_id}")


def load_product_contracts(
    root: Path | str | None = None,
    *,
    require_staged_fastas: bool = True,
) -> dict:
    if type(require_staged_fastas) is not bool:
        raise ValueError("require_staged_fastas must be a bool")
    root = Path(root) if root is not None else Path(__file__).resolve().parents[1]
    contract_file = root / CONTRACT_PATH
    if not contract_file.is_file():
        raise ValueError(f"Missing or non-regular product registry file: {contract_file}")
    try:
        raw_json = contract_file.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"Product registry is not valid UTF-8: {contract_file}") from exc
    payload = json.loads(raw_json, object_pairs_hook=_reject_duplicate_json_keys, parse_constant=_reject_nonfinite_product_json, parse_float=_parse_finite_product_json_float)
    if not isinstance(payload, dict) or not isinstance(payload.get("products"), list):
        raise ValueError("Product registry must contain a products list")
    if require_staged_fastas and ("mcb_donor_blocking_policy" not in payload or "design_b_four_state_lifecycle_policy" not in payload):
        raise ValueError("Product registry missing required mcb_donor_blocking_policy or design_b_four_state_lifecycle_policy")
    if "mcb_donor_blocking_policy" in payload:
        pol = payload["mcb_donor_blocking_policy"]
        if (
            not isinstance(pol, dict)
            or type(pol.get("min_independent_abo_oo_seedlines")) is not int
            or pol["min_independent_abo_oo_seedlines"] < 2
            or not isinstance(pol.get("kinship_pi_hat_lt"), (int, float))
            or isinstance(pol.get("kinship_pi_hat_lt"), bool)
            or not math.isfinite(float(pol["kinship_pi_hat_lt"]))
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
            or any(not isinstance(v, dict) for v in states.values())
            or states["STATE_1_STAGED_REFERENCE_FASTA_ARTIFACT"].get("fasta_sha256") != DESIGN_B_FASTA_SHA256
            or (
                "construct_count" in states["STATE_1_STAGED_REFERENCE_FASTA_ARTIFACT"]
                and (
                    type(states["STATE_1_STAGED_REFERENCE_FASTA_ARTIFACT"]["construct_count"]) is not int
                    or states["STATE_1_STAGED_REFERENCE_FASTA_ARTIFACT"]["construct_count"] != 6
                )
            )
            or states["STATE_2_REGISTRY_DONOR_REFERENCE_BOUND"].get("bound_c1_c5_pair_sha256") != DESIGN_B_BOUND_C1_C5_PAIR_SHA256
            or states["STATE_3_STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC"].get("execution_supported") is not False
            or type(states["STATE_4_EXECUTION_SUPPORTED_QUALIFIED_PRODUCT"].get("currently_achieved_family_b_count")) is not int
            or states["STATE_4_EXECUTION_SUPPORTED_QUALIFIED_PRODUCT"].get("currently_achieved_family_b_count") != 0
        ):
            raise ValueError("Invalid design_b_four_state_lifecycle_policy in product registry (IIP-0026; DEC-105)")

    products = payload["products"]
    if not products or any(not isinstance(p, dict) for p in products):
        raise ValueError("Product registry products must be a nonempty list of dictionaries")
    for product in products:
        validate_product_definition(product)
    if len({p["product_variant_id"] for p in products}) != len(products):
        raise ValueError("Product variants must be nonempty and unique")
    identifiers: set[str] = set()
    by_id: dict[str, dict] = {}
    for product in products:
        by_id[product["product_variant_id"]] = product
        for identifier in [product["product_variant_id"], *product.get("aliases", [])]:
            if identifier in identifiers:
                raise ValueError("Product identifiers and aliases must be globally unique")
            identifiers.add(identifier)

    missing_b2 = STAGED_DESIGN_B2_PRODUCT_IDS - by_id.keys()
    if missing_b2:
        raise ValueError(f"Missing required staged Design B2 product definitions: {sorted(missing_b2)} (IIP-0026; DEC-105)")
    missing_l2_v3 = STAGED_LOCUS2_V3_PRODUCT_IDS - by_id.keys()
    if missing_l2_v3:
        raise ValueError(f"Missing required staged Locus-2 v3 product definitions: {sorted(missing_l2_v3)} (IIP-0031; DEC-110)")
    missing_dy2 = STAGED_DESIGN_Y_V2_PRODUCT_IDS - by_id.keys()
    if missing_dy2:
        raise ValueError(f"Missing required staged Design Y-FullStack v2 product definitions: {sorted(missing_dy2)} (IIP-0034; DEC-113)")

    fasta_file = root / DESIGN_B_FASTA_PATH
    if require_staged_fastas or fasta_file.exists():
        if not fasta_file.is_file():
            raise ValueError(f"Missing or non-regular staged Design B FASTA: {fasta_file} (IIP-0026; DEC-105)")
        raw_bytes = fasta_file.read_bytes()
        if hashlib.sha256(raw_bytes).hexdigest() != DESIGN_B_FASTA_SHA256:
            raise ValueError("Staged Design B FASTA SHA-256 mismatch against registry pin (IIP-0026; DEC-105)")
        try:
            seq_map = _parse_fasta_map(raw_bytes.decode("utf-8"))
        except UnicodeDecodeError as exc:
            raise ValueError("Staged Design B FASTA is not valid UTF-8 (IIP-0026; DEC-105)") from exc
        c1_seq = seq_map.get(DESIGN_B_C1_ID, "")
        c5_seq = seq_map.get(DESIGN_B_C5_ID, "")
        pair_sha = hashlib.sha256((c1_seq + "|" + c5_seq).encode("utf-8")).hexdigest()
        if (
            len(c1_seq) != DESIGN_B_C1_LENGTH_BP
            or len(c5_seq) != DESIGN_B_C5_LENGTH_BP
            or hashlib.sha256(c1_seq.encode("utf-8")).hexdigest() != DESIGN_B_C1_SHA256
            or hashlib.sha256(c5_seq.encode("utf-8")).hexdigest() != DESIGN_B_C5_SHA256
            or pair_sha != DESIGN_B_BOUND_C1_C5_PAIR_SHA256
        ):
            raise ValueError("Staged Design B c1+c5 pair digest or length mismatch against registry pin (IIP-0026; DEC-105)")

    l2_fasta_file = root / LOCUS2_V3_FASTA_PATH
    if require_staged_fastas or l2_fasta_file.exists():
        if not l2_fasta_file.is_file():
            raise ValueError(f"Missing or non-regular staged Locus-2 CIITA Mode-B v3 FASTA: {l2_fasta_file} (IIP-0031; DEC-110)")
        l2_raw = l2_fasta_file.read_bytes()
        if hashlib.sha256(l2_raw).hexdigest() != LOCUS2_V3_FASTA_SHA256:
            raise ValueError("Staged Locus-2 CIITA Mode-B v3 FASTA SHA-256 mismatch against registry pin (IIP-0031; DEC-110)")
        try:
            l2_map = _parse_fasta_map(l2_raw.decode("utf-8"))
        except UnicodeDecodeError as exc:
            raise ValueError("Staged Locus-2 CIITA Mode-B v3 FASTA is not valid UTF-8 (IIP-0031; DEC-110)") from exc
        for _pvid, _spec in LOCUS2_V3_DONOR_BINDINGS.items():
            _did = str(_spec["locus2_donor_id"])
            _seq = l2_map.get(_did, "")
            _exp_len = _spec["donor_length_bp"]
            if (
                type(_exp_len) is not int
                or len(_seq) != _exp_len
                or hashlib.sha256(_seq.encode("utf-8")).hexdigest() != str(_spec["locus2_donor_sha256"])
            ):
                raise ValueError(f"Staged Locus-2 donor digest or length mismatch for {_did} (IIP-0031; DEC-110)")

    dy2_fasta_file = root / DESIGN_Y_V2_FASTA_PATH
    if require_staged_fastas or dy2_fasta_file.exists():
        if not dy2_fasta_file.is_file():
            raise ValueError(f"Missing or non-regular staged Design Y-FullStack v2 FASTA: {dy2_fasta_file} (IIP-0034; DEC-113)")
        dy2_raw = dy2_fasta_file.read_bytes()
        if hashlib.sha256(dy2_raw).hexdigest() != DESIGN_Y_V2_FASTA_SHA256:
            raise ValueError("Staged Design Y-FullStack v2 FASTA SHA-256 mismatch against registry pin (IIP-0034; DEC-113)")
        try:
            dy2_map = _parse_fasta_map(dy2_raw.decode("utf-8"))
        except UnicodeDecodeError as exc:
            raise ValueError("Staged Design Y-FullStack v2 FASTA is not valid UTF-8 (IIP-0034; DEC-113)") from exc
        for _pvid, _spec in DESIGN_Y_V2_DONOR_BINDINGS.items():
            _rid = str(_spec["fasta_record_id"])
            _seq = dy2_map.get(_rid, "")
            if len(_seq) != int(_spec["donor_length_bp"]) or hashlib.sha256(_seq.encode("utf-8")).hexdigest() != str(_spec["donor_sha256"]):
                raise ValueError(f"Staged Design Y-FullStack v2 record digest or length mismatch for {_rid} (IIP-0034; DEC-113)")

    return payload


def validate_product_definition(product: dict) -> dict:
    if not isinstance(product, dict):
        raise ValueError("Product definition must be a dictionary")
    required = {
        "product_variant_id", "family", "dsb_or_indel_loci", "base_or_prime_edit_loci",
        "alleles", "cassette_topology", "K", "B", "sequence_materialization_status",
        "physical_qualification_status", "execution_supported", "donor_lineage",
        "hla_retention", "pending_assay_gates", "physical_manufacturing_yield",
    }
    if required - product.keys():
        raise ValueError(f"Incomplete product definition: {sorted(required - product.keys())}")
    pvid = product["product_variant_id"]
    if not isinstance(pvid, str) or pvid != pvid.strip() or not re.fullmatch(r"ISLET1_[A-Z0-9_]+", pvid):
        raise ValueError("Invalid canonical product identifier")
    fam = product["family"]
    if not isinstance(fam, str) or fam not in {"A", "B"}:
        raise ValueError("Unknown product family")
    if type(product["execution_supported"]) is not bool:
        raise ValueError("execution_supported must be boolean")
    if not isinstance(product["sequence_materialization_status"], str):
        raise ValueError("sequence_materialization_status must be a string")
    if not isinstance(product["physical_qualification_status"], str):
        raise ValueError("physical_qualification_status must be a string")
    if "lifecycle_state" in product and (
        not isinstance(product["lifecycle_state"], str)
        or product["lifecycle_state"] not in ALLOWED_DESIGN_B_LIFECYCLE_STATES
    ):
        raise ValueError("Invalid lifecycle_state in product definition (IIP-0026; DEC-105)")
    if not (
        pvid.startswith("ISLET1_" + fam)
        or (fam == "A" and pvid.startswith("ISLET1_COMP_"))
        or (fam == "A" and pvid.startswith("ISLET1_Z_") and product["execution_supported"] is False)
        or (fam == "A" and pvid.startswith("ISLET1_Y_") and product["execution_supported"] is False)
    ):
        raise ValueError("Product identifier conflicts with its family")
    if product["physical_manufacturing_yield"] is not None:
        raise ValueError("Physical manufacturing yield remains unmeasured in this release")
    if "aliases" in product and not isinstance(product["aliases"], list):
        raise ValueError("Aliases must be a list of bounded identifiers")
    aliases = product.get("aliases", [])
    if any(not isinstance(a, str) or not a or a != a.strip() or not re.fullmatch(r"[A-Z0-9_]+", a) for a in aliases):
        raise ValueError("Aliases must be bounded identifiers")
    if len(set(aliases)) != len(aliases) or pvid in aliases:
        raise ValueError("Product aliases must be unique and distinct from the canonical identifier")
    for key in ("alleles", "cassette_topology"):
        value = product[key]
        if not isinstance(value, dict) or not value or any(not isinstance(k, str) or not k.strip() or not isinstance(v, str) or not v.strip() for k, v in value.items()):
            raise ValueError(f"{key} must be a nonempty explicit mapping")
    for key in ("donor_lineage", "hla_retention"):
        if not isinstance(product[key], dict) or not product[key]:
            raise ValueError(f"{key} must explicitly define donor/HLA scope")
    if "donor_seedline_id_required" in product["donor_lineage"] and type(product["donor_lineage"]["donor_seedline_id_required"]) is not bool:
        raise ValueError("donor_seedline_id_required must be boolean (IIP-0024; DEC-103)")
    if "kinship_pi_hat_lt" in product["donor_lineage"] and (
        isinstance(product["donor_lineage"]["kinship_pi_hat_lt"], bool)
        or not isinstance(product["donor_lineage"]["kinship_pi_hat_lt"], (int, float))
        or not math.isfinite(float(product["donor_lineage"]["kinship_pi_hat_lt"]))
        or not (0.0 < float(product["donor_lineage"]["kinship_pi_hat_lt"]) <= 0.05)
    ):
        raise ValueError("kinship_pi_hat_lt must be a finite number in (0, 0.05] (IIP-0024; DEC-103)")
    if "physical_route_id" in product and (not isinstance(product["physical_route_id"], str) or not re.fullmatch(r"[A-Z0-9_]+", product["physical_route_id"])):
        raise ValueError("physical_route_id must be a bounded uppercase identifier (IIP-0024; DEC-103)")
    for k_bp, v_bp in product.items():
        if isinstance(k_bp, str) and k_bp.endswith("_bp") and v_bp is not None:
            if type(v_bp) is not int or v_bp <= 0:
                raise ValueError(f"{k_bp} must be a positive integer")
    gates = product["pending_assay_gates"]
    if not isinstance(gates, list) or not gates or any(not isinstance(g, str) or not g.strip() for g in gates):
        raise ValueError("Explicit pending assay gates are required")
    known_loci = {"B2M", "CIITA", "CXCL10", "AAVS1", "HLA-A", "HLA-B", "TYK2", "INS", "CLYBL"}
    for loci_key, count_key in [("dsb_or_indel_loci", "K"), ("base_or_prime_edit_loci", "B")]:
        loci = product[loci_key]
        if not isinstance(loci, list) or any(not isinstance(x, str) or x not in known_loci for x in loci) or len(loci) != len(set(loci)):
            raise ValueError(f"{loci_key} must list unique loci")
        if type(product[count_key]) is not int or product[count_key] != len(loci):
            raise ValueError(f"{count_key} does not match the explicit locus list")
    if not (set(product["dsb_or_indel_loci"]) | set(product["base_or_prime_edit_loci"])) <= product["alleles"].keys():
        raise ValueError("Every edited locus needs an explicit allele definition")

    # Enforce 4-state Design B lifecycle while preserving Family A reference-compiler semantics (IIP-0026; DEC-105)
    # and binding all 8 unscored Locus-2 CIITA Mode-B v3 research comparators (IIP-0031; DEC-110):
    if pvid in STAGED_LOCUS2_V3_PRODUCT_IDS:
        l2_exp = LOCUS2_V3_DONOR_BINDINGS[pvid]
        expected_short_alias = str(l2_exp["locus2_donor_id"])
        if expected_short_alias not in aliases:
            raise ValueError(f"Locus-2 v3 comparator {pvid} must carry pinned donor alias {expected_short_alias} (IIP-0031; DEC-110)")
        for a in aliases:
            if a in STAGED_LOCUS2_V3_PRODUCT_IDS or (a in STAGED_LOCUS2_V3_DONOR_ALIASES and a != expected_short_alias):
                raise ValueError(f"Locus-2 v3 comparator {pvid} claims reserved staged identifier/alias {a!r} (IIP-0031; DEC-110)")
        if fam != l2_exp["family"]:
            raise ValueError("Locus-2 v3 research comparator family conflicts with registry binding (IIP-0031; DEC-110)")
        if product["execution_supported"] is not False:
            raise ValueError("Locus-2 v3 research comparators must keep execution_supported=False (IIP-0031; DEC-110)")
        if product.get("lifecycle_state") != "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC":
            raise ValueError("Locus-2 v3 research comparators must stay at STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC (IIP-0031; DEC-110)")
        if fam == "A":
            if product["sequence_materialization_status"] != "REFERENCE_DONORS_AVAILABLE_REQUIRES_EDITED_ALLELE_QC":
                raise ValueError("Family A Locus-2 v3 comparator requires REFERENCE_DONORS_AVAILABLE_REQUIRES_EDITED_ALLELE_QC (IIP-0031; DEC-110)")
        else:
            if product["sequence_materialization_status"] != "STAGED_RESEARCH_COMPARATOR_DONORS_BOUND_PENDING_ALLELE_QC":
                raise ValueError("Family B Locus-2 v3 comparator requires STAGED_RESEARCH_COMPARATOR_DONORS_BOUND_PENDING_ALLELE_QC (IIP-0031; DEC-110)")
            _validate_staged_design_b_binding(product.get("staged_design_b_donor_reference_binding"), pvid)
        l2_bind = product.get("staged_locus2_modeb_donor_binding")
        if (
            not isinstance(l2_bind, dict)
            or l2_bind.get("fasta_path") != str(LOCUS2_V3_FASTA_PATH)
            or l2_bind.get("fasta_sha256") != LOCUS2_V3_FASTA_SHA256
            or l2_bind.get("locus2_donor_id") != l2_exp["locus2_donor_id"]
            or l2_bind.get("locus2_donor_sha256") != l2_exp["locus2_donor_sha256"]
            or type(l2_bind.get("donor_length_bp")) is not int
            or l2_bind.get("donor_length_bp") != l2_exp["donor_length_bp"]
            or type(l2_bind.get("nominal_packaged_aav6_bp")) is not int
            or l2_bind.get("nominal_packaged_aav6_bp") != l2_exp["nominal_packaged_aav6_bp"]
            or l2_bind.get("background_id") != l2_exp["background_id"]
            or type(l2_bind.get("has_intronic_mir")) is not bool
            or l2_bind.get("has_intronic_mir") is not l2_exp["has_intronic_mir"]
            or l2_bind.get("mir_class") != l2_exp["mir_class"]
            or l2_bind.get("unit2_cargo") != l2_exp["unit2_cargo"]
            or l2_bind.get("unit2_promoter") != l2_exp["unit2_promoter"]
            or l2_bind.get("barrier_variant_type") != l2_exp["barrier_variant_type"]
            or l2_bind.get("rule_s11_gate_id") != RULE_S11_GATE_ID
        ):
            raise ValueError(f"Invalid or mismatched staged_locus2_modeb_donor_binding for {pvid} (IIP-0031; DEC-110)")
    else:
        if "staged_locus2_modeb_donor_binding" in product:
            raise ValueError(f"Non-staged product {pvid} must not carry staged_locus2_modeb_donor_binding (IIP-0031; DEC-110)")
        for a in aliases:
            if a in STAGED_LOCUS2_V3_PRODUCT_IDS or a in STAGED_LOCUS2_V3_DONOR_ALIASES or a.startswith("L2_V3_"):
                raise ValueError(f"Non-staged product {pvid} must not claim reserved Locus-2 v3 identifier/alias {a!r} (IIP-0031; DEC-110)")
        if fam == "A":
            if product["sequence_materialization_status"] != "REFERENCE_DONORS_AVAILABLE_REQUIRES_EDITED_ALLELE_QC":
                raise ValueError("Sequence materialization status conflicts with the product family")
        elif pvid in STAGED_DESIGN_B2_PRODUCT_IDS:
            if product["sequence_materialization_status"] != "STAGED_RESEARCH_COMPARATOR_DONORS_BOUND_PENDING_ALLELE_QC":
                raise ValueError("Staged Design B2 variants require STAGED_RESEARCH_COMPARATOR_DONORS_BOUND_PENDING_ALLELE_QC (IIP-0026; DEC-105)")
            if product.get("lifecycle_state") != "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC":
                raise ValueError("Staged Design B2 variants must remain at State 3 STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC (IIP-0026; DEC-105)")
            _validate_staged_design_b_binding(product.get("staged_design_b_donor_reference_binding"), pvid)
        else:
            if product["sequence_materialization_status"] != "NOT_MATERIALIZED_SEPARATE_DONOR_LEDGER_REQUIRED":
                raise ValueError("Sequence materialization status conflicts with the product family")

    if "staged_design_b_donor_pair_binding" in product:
        _validate_staged_design_b_binding(product["staged_design_b_donor_pair_binding"], pvid)
    if pvid in STAGED_DESIGN_Y_V2_PRODUCT_IDS:
        _dy2_exp = DESIGN_Y_V2_DONOR_BINDINGS[pvid]
        if product["execution_supported"] is not False:
            raise ValueError(f"Staged Design Y-FullStack v2 comparator {pvid} must keep execution_supported=False (IIP-0034; DEC-113)")
        if product.get("lifecycle_state") != "STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC":
            raise ValueError(f"Staged Design Y-FullStack v2 comparator {pvid} must stay at STAGE0_RESEARCH_COMPARATOR_PENDING_ALLELE_QC (IIP-0034; DEC-113)")
        if product.get("physical_route_id") != "ROUTE_3P_MODULAR_PLANAR_ENVELOPE":
            raise ValueError(f"Staged Design Y-FullStack v2 comparator {pvid} must declare physical_route_id=ROUTE_3P_MODULAR_PLANAR_ENVELOPE (IIP-0034; DEC-113)")
        _dy2_bind = product.get("staged_design_y_fullstack_v2_donor_binding")
        _rec_sha = hashlib.sha256(json.dumps(product, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()
        if _rec_sha != _dy2_exp["full_record_sha256"]:
            raise ValueError(f"Invalid or mutated Design Y-FullStack v2 product record {pvid} (full_record_sha256 mismatch; IIP-0034; DEC-113)")
        if (
            not isinstance(_dy2_bind, dict)
            or _dy2_bind.get("proposal_id") != _dy2_exp["proposal_id"]
            or _dy2_bind.get("decision_id") != _dy2_exp["decision_id"]
            or _dy2_bind.get("fasta_path") != str(DESIGN_Y_V2_FASTA_PATH)
            or _dy2_bind.get("fasta_sha256") != DESIGN_Y_V2_FASTA_SHA256
            or _dy2_bind.get("fasta_record_id") != _dy2_exp["fasta_record_id"]
            or _dy2_bind.get("donor_sha256") != _dy2_exp["donor_sha256"]
            or type(_dy2_bind.get("donor_length_bp")) is not int
            or _dy2_bind.get("donor_length_bp") != _dy2_exp["donor_length_bp"]
            or type(_dy2_bind.get("tu_count")) is not int
            or _dy2_bind.get("tu_count") != _dy2_exp["tu_count"]
            or type(_dy2_bind.get("pigk_ectopic_gpi_count")) is not int
            or _dy2_bind.get("pigk_ectopic_gpi_count") != _dy2_exp["pigk_ectopic_gpi_count"]
            or type(_dy2_bind.get("encoded_2a_motifs_6frame")) is not int
            or _dy2_bind.get("encoded_2a_motifs_6frame") != _dy2_exp["encoded_2a_motifs_6frame"]
            or _dy2_bind.get("role") != _dy2_exp["role"]
            or _dy2_bind.get("stage0_promotion_gates") != _dy2_exp["stage0_promotion_gates"]
        ):
            raise ValueError(f"Invalid or missing staged_design_y_fullstack_v2_donor_binding for {pvid} (IIP-0034; DEC-113)")
    elif "staged_design_y_fullstack_v2_donor_binding" in product:
        raise ValueError(f"Non-Design-Y-v2 product {pvid} must not carry staged_design_y_fullstack_v2_donor_binding (IIP-0034; DEC-113)")
    if product["physical_qualification_status"] != "PREREGISTERED_PENDING":
        raise ValueError("This release contains no physically qualified exact product")
    if fam == "B" and product["execution_supported"]:
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
            if (
                type(v_rec.get("donor_length_bp")) is not int
                or v_rec.get("donor_length_bp") != exp_len
                or type(v_rec.get("packaged_aav6_length_bp")) is not int
                or v_rec.get("packaged_aav6_length_bp") != exp_pkg_len
            ):
                raise ValueError(f"Invalid donor/packaged length in {v_key}")
            if "l1_linker_aa" in v_rec and type(v_rec["l1_linker_aa"]) is not int:
                raise ValueError(f"Invalid l1_linker_aa type in {v_key}")
            sha_val = v_rec.get("donor_sha256")
            if not isinstance(sha_val, str) or not re.fullmatch(r"[0-9a-f]{64}", sha_val) or sha_val != exp_sha:
                raise ValueError(f"Invalid or mismatched donor_sha256 in {v_key}: {sha_val}")
    return product


def get_product_definition(
    product_id: str,
    root: Path | str | None = None,
    *,
    require_staged_fastas: bool = True,
) -> dict:
    if not isinstance(product_id, str) or not product_id or product_id != product_id.strip():
        raise ValueError(f"Unknown product variant: {product_id!r}")
    for product in load_product_contracts(root, require_staged_fastas=require_staged_fastas)["products"]:
        if product_id == product["product_variant_id"] or product_id in product.get("aliases", []):
            return deepcopy(product)
    raise ValueError(f"Unknown product variant: {product_id}")


def validate_product_definition_dict(pdef: dict, root: Path | str | None = None) -> dict:
    if not isinstance(pdef, dict):
        raise ValueError("Product definition must be a dictionary")
    validate_product_definition(pdef)
    canonical = get_product_definition(pdef["product_variant_id"], root=root)
    if not _strict_json_equal(pdef, canonical):
        raise ValueError(
            f"Supplied product_def for {pdef['product_variant_id']!r} diverges from authoritative registry definition"
        )
    return deepcopy(canonical)
