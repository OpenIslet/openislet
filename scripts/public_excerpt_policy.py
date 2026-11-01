"""Bounded opaque redaction references, never a free-text exemption channel."""
from __future__ import annotations

import hashlib
import re

FIELDS = frozenset({
    'verbatim_evidence_quote', 'focal_factor_X', 'outcome_variable_Y', 'comparator_X0',
    'exact_variant_or_swap_description', 'level_1_condition_and_effect',
    'level_2_condition_and_effect', 'mechanistic_explanation', 'evidence_source',
    'observation', 'verbatim_quote', 'statement', 'hop1_evidence', 'hop2_evidence',
    'source_span', 'field', 'omitted_text', 'description', 'note', 'notes',
    'failure_mode_description', 'mechanism', 'condition', 'outcome', 'rationale',
    'quoted_threshold_or_constraint', 'threshold_description', 'effect_description',
})
GENES = frozenset('B2M CASP9 CD200 CD274 CD47 CD55 CD59 THBD CHGA CIITA CXADR CXCL10 ENTPD1 FASLG G6PC2 GAD2 GLIS3 HLA-A HLA-B HLA-C HLA-E HLA-G IAPP IFIH1 IL10 TGFB1 IL2 INS ISG15 JAK1 JAK2 MANF NLRC5 PTPN2 PTPRN RIPK1 RNLS SEC31A SLC30A8 SOCS1 STAT1 TET2 TGM2 TNFAIP3 TXNIP TYK2 SERPINB9 PD-L1 ICASP9 DDIT3 HYOU1 PCSK1 PTPN22 MAFA PDX1 NKX6-1 GCG SST UTP23 HNF1A'.split())
MODALITIES = frozenset('CRISPR_Cas12a_Multiplex_KO BaseEdit_ABE_or_CBE_Epitope_or_Splice PrimeEdit Signal_Peptide_or_Domain_Swap Endogenous_SwapIn_KO_plus_KI SafeHarbor_KI_AAVS1_GAPDH_CLYBL Lentiviral_or_Transposon_OE Transient_mRNA_or_Surface_Tethering CRISPRi_or_shRNA_Knockdown CRISPRa_Activation SmallMolecule_Pharmacological_Modulator none NONE CRISPR_Cas9_KO CRISPR_Cas12a_KO CRISPRi_KD CRISPRa_OE Lentiviral_OE Transposon_OE HDR_KI Base_Editing Prime_Editing Pharmacologic_Inhibitor Pharmacologic_Activation Antibody_Blockade RNAi_KD KO KI'.split())


def is_bounded_metadata(value: object, field: str, *, bibliographic: bool = False) -> bool:
    """Exempt validated identifier values, never an arbitrary metadata key.

    Bibliographic prose has a separate bounded scope inside study_metadata;
    graph properties with those names do not inherit that scope.
    """
    if not isinstance(value,str):return False
    if field=='pmid':return re.fullmatch(r'[1-9][0-9]{0,8}',value) is not None
    if field=='pmcid':return re.fullmatch(r'PMC[1-9][0-9]{0,9}',value) is not None
    if field in {'study','hop1_study','hop2_study'}:
        return len(value)<=160 and re.fullmatch(r'PMID:[1-9][0-9]{0,8}(?: \([^\W\d_][^\x00-\x1f\x7f()]{0,120} [12][0-9]{3}\))?',value) is not None
    if field=='schema_version':return re.fullmatch(r'[0-9]{1,3}(?:\.[0-9]{1,3}){1,2}',value) is not None
    if field=='doi':return len(value)<=256 and re.fullmatch(r'10\.[0-9]{4,9}/[^\s\x00-\x1f\x7f]+',value) is not None
    if not bibliographic:return False
    if field=='citation_short':
        return len(value)<=160 and re.fullmatch(r'[^\W\d_][^\x00-\x1f\x7f()]{0,120} [12][0-9]{3}',value) is not None
    limits={'title':1000,'journal':200}
    return field in limits and 0<len(value)<=limits[field] and re.search(r'[\x00-\x1f\x7f]',value) is None


def validate_token(text: str, *, expected_pmid: str | None = None,
                   source_quote: str | None = None,
                   expected_field: str | None = None) -> bool:
    """Validate fixed metadata and, when available, the actual source digest.

    A missing private source permits an opaque reference, never a claim that its
    acquisition or scientific meaning has been independently authenticated.
    """
    if not isinstance(text, str) or len(text) > 1500:
        return False
    m = re.fullmatch(r'\[REDACTED_NON_OA_EXCERPT: see PMID:([0-9]{7,8}) ([A-Za-z0-9_=, .-]+)\]', text)
    if m is None or (expected_pmid is not None and m[1] != str(expected_pmid)):
        return False
    values = {}
    for part in m[2].split():
        if part.count('=') != 1:
            return False
        key, value = part.split('=', 1)
        if key in values or key not in {'field', 'mod', 'genes', 'pert', 'quote_sha256'} or not value:
            return False
        values[key] = value
    if not {'field', 'quote_sha256'} <= values.keys() or values['field'] not in FIELDS:
        return False
    if expected_field is not None and values['field'] != safe_field(expected_field):
        return False
    digest = values['quote_sha256']
    if not re.fullmatch(r'[0-9a-f]{64}', digest) or digest == '0' * 64:
        return False
    if 'genes' in values and values['genes'] != 'NONE':
        genes = values['genes'].split(',')
        if len(genes) > 32 or len(genes) != len(set(genes)) or any(g not in GENES for g in genes):
            return False
    if values.get('pert', 'NONE') not in {'GOF_OE', 'KO_LOF', 'BOTH', 'NONE'}:
        return False
    if 'mod' in values and values['mod'] not in MODALITIES:
        return False
    return source_quote is None or (isinstance(source_quote, str) and hashlib.sha256(source_quote.encode('utf-8')).hexdigest() == digest)


def safe_field(name: str) -> str:
    return name if name in FIELDS else 'omitted_text'
