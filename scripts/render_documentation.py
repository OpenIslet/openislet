"""Render current scientific facts from artifacts while retaining design rationale.

No numeric value in this module is a surrogate assay result. Historical details
remain in their source documents; the generated block identifies current scope.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import math
import os
import stat

START = '<!-- OPENISLET_GENERATED_FACTS_START -->'
END = '<!-- OPENISLET_GENERATED_FACTS_END -->'
CANDIDATE = 'HNRNPA2B1-proximal 455-bp CpG-core candidate (GRCh38 chr7:26200580-26201034)'
RNLS = ('MacDonald 2025 (PMID:39988068 Fig. 6) reports poor dynamic GSIS in both WT and RNLS-KO iPSC-beta aggregates. '
        'Higher basal OCR in KO motivates a precautionary metabolic penalty, not a demonstrated knockout-specific GSIS decrement; '
        'the exact-product dynamic-perifusion comparison is pending.')
DELIVERY = ('Shy 2023 (PMID:36008610) supports a nonviral ssDNA/dsDNA editing strategy in its tested cell context. '
            'Exact 4,880-bp hPSC donor delivery, ceDNA and mammalian bridge-recombinase transfer are experimental delivery options requiring Stage-0 qualification.')
STATUTE = ('Section 356k(e) concerns actions to expedite development and review of an application using a designated platform; '
           'designation requests and determinations are governed by sections 356k(c)-(d). Section 356k(f) concerns leveraging data, '
           'subject to the same-sponsor or granted right-of-reference conditions. This is a citation correction, not a finding of product eligibility or approval. '
           '[Official statute](https://www.govinfo.gov/link/uscode/21/356k).')
CHS4 = ('engineered 250-bp cHS4-derived variant: 235-bp exact contiguous match at insert 3..237 to U78775.2 positions 10..244, '
        'with 2 bp before and 13 bp after that block; the alignment ledger discloses adjacent native T245 omission and native/adapted tail classification. '
        'Those 15 outside-block bases are not all established synthetic adapters; insulation requires bench qualification')
P2A = 'single 22-aa GSG-P2A element (21-aa upstream stub plus downstream Pro; legacy tPT2A alias)'
STRUCTURE1 = ('Reference-geometry compatibility screen: 3CDG is a CD94/NKG2A complex. P8-Phe8 contributes 24 heavy-atom contacts to shared CD94 and 0 to NKG2A in copy 1, and 21/0 in copy 2. '
              'Removing Phe ring atoms is an atom-truncation diagnostic, not a relaxed Val8 model. This structure cannot discriminate NKG2A from NKG2C. '
              'Receptor-specific K_D values are working soluble-ectodomain reference priors; direct SPR and HCMV+ NKG2C+ NK-donor assays remain preregistered.')
STRUCTURE2 = ('Reference SIRPA V1/V2 sequence and geometry comparison; pan-allelic protection is unverified. '
              'Asp101 is inserted in SIRPA V1 (4CMM) relative to V2 (2JJS). Soluble crystal CD47 carries C33G/Cys15-to-Gly; engineered CD47 retains WT Cys15. '
              'Reference soluble-C15G structures place CD47 Val115 about 30.93 A from the nearest SIRPA atom. '
              'This coarse anchor-distance screen does not model membrane orientation, the engineered WT-Cys15/GPI stalk or V8; binding, membrane and phagocytosis assays remain pending.')


def _safe_location(path):
    for parent in [path.absolute(), *path.absolute().parents]:
        if parent.is_symlink():
            raise ValueError('Documentation paths and ancestors may not be symbolic links')


def _read_regular(path):
    _safe_location(path)
    info=path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError('Documentation inputs must be single-link regular files')
    fd=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
    try:
        opened=os.fstat(fd)
        if (opened.st_dev,opened.st_ino)!=(info.st_dev,info.st_ino):
            raise ValueError('Documentation input changed identity while opening')
        with os.fdopen(fd,'rb',closefd=False) as stream:data=stream.read()
        after=os.fstat(fd)
        if (after.st_size,after.st_mtime_ns,after.st_ctime_ns)!=(len(data),opened.st_mtime_ns,opened.st_ctime_ns):
            raise ValueError('Documentation input changed while reading')
        return data
    finally:os.close(fd)


def read_json(root, relative, default=None):
    path = root / relative
    return json.loads(_read_regular(path)) if path.exists() else ({} if default is None else default)


def strip_trajectories(value):
    if isinstance(value, dict):
        return {k: strip_trajectories(v) for k,v in value.items() if not any(t in k for t in ('trajectory', 'hourly_grid', 'sampled_hourly'))}
    if isinstance(value, list): return [strip_trajectories(v) for v in value]
    return value


def derive_facts(root: Path, manifest_pin_count: int) -> dict:
    source_bytes={}
    def load(relative, default=None):
        key=Path(relative).as_posix();path=root/key
        if not path.exists():return {} if default is None else default
        if key not in source_bytes:source_bytes[key]=_read_regular(path)
        return json.loads(source_bytes[key])
    discovery_path=root/'data/raw/corpus_manifest_all_sources.jsonl'
    discovery_key=discovery_path.relative_to(root).as_posix()
    source_bytes[discovery_key]=_read_regular(discovery_path)
    discovery_rows=[json.loads(line) for line in source_bytes[discovery_key].decode('utf-8').splitlines() if line.strip()]
    discovery={'jsonl_rows':len(discovery_rows),'unique_pubmed_pmids':len({str(row['pmid']) for row in discovery_rows if row.get('pmid')}),'inline_abstract_strings':sum(bool(row.get('abstract')) for row in discovery_rows),'scope':'PubMed PMID records only; patent, preprint, GEO and trial inventories are separate source files, not additive components of this JSONL.'}
    paths = {'graph':'data/extracted/cross_study_discovery_graph.json',
             'rankings':'data/extracted/islet1_candidate_stack_rankings.json',
             'protocol':'data/protocols/hip_core_fractional_factorial_v0.6.2.json',
             'layer2':'data/extracted/layer2_regulatory_and_perturbation_atlas.json',
             'layer3':'data/extracted/layer3_structural_verification_atlas.json',
             'layer4':'data/extracted/layer4_tissue_and_metabolic_twin_atlas.json',
             'layer5':'data/extracted/layer5_personalized_patient_and_crispr_atlas.json'}
    data = {name:load(path) for name,path in paths.items()}
    studies=[load(p.relative_to(root)) for p in sorted((root/'data/extracted').glob('*.json')) if p.stem.isdigit()]
    selected_gencode_models=load('data/raw/layer2_genomic_and_jaspar_reference.json').get('gencode_v46_transcript_models_pinned',{})
    raw_edits=[edit for study in studies for arm in study.get('arms',[]) for edit in arm.get('gate2_3_4_genetic_and_epitope_edits',[])]
    ledger=data['rankings'].get('layer1_ultracode_audit_reconciliation_ledger',{})
    stage=data['layer2'].get('stage0_to_stage7_perturbation_plane',{})
    profiles=stage.get('perturbation_profiles',[])
    construct=data['layer2'].get('construct_sequence_and_grammar_verification',{})
    pde=data['layer4'].get('part_a_336h_spatial_reaction_diffusion_pouch_pde',{})
    regimes=pde.get('regimes_simulated',[])
    lead=next((r for r in regimes if r.get('regime_id')=='LEAD_CENTRAL_LAMINA_PDMS_CAO2_F045_PLUS_CATALASE'),{})
    ode=data['layer4'].get('part_b_24h_trihormonal_paracrine_and_whole_body_metabolic_ode',{})
    patient=data['layer5'].get('part_b_personalized_patient_compiler_and_10k_population_trial',{})
    return {'schema_version':'0.8.0','epistemic_status':'GENERATED_COMPUTATIONAL_FACTS_PHYSICAL_QUALIFICATION_PENDING',
            'manifest_pin_count':manifest_pin_count,'discovery_manifest':discovery,
            'source_sha256':{path:hashlib.sha256(data).hexdigest() for path,data in source_bytes.items()},
            'corpus':{'studies':len(studies),'arms':sum(len(s.get('arms',[])) for s in studies),
                      'reference_control_arms':sum(bool(a.get('is_reference_control')) for s in studies for a in s.get('arms',[])),
                      'atomic_effects':sum(len(s.get('atomic_effects',[])) for s in studies),
                      'raw_edit_entries':len(raw_edits),'total_numeric_lrr_effects':sum(isinstance(e.get('log_response_ratio_or_smd'),(int,float)) for st in studies for e in st.get('atomic_effects',[])),'expanded_gene_occurrences':ledger.get('source_counts',{}).get('expanded_edit_gene_occurrences'),
                      'heuristic_penalty_effects':sum(bool(e.get('secretome_or_differentiation_penalty_observed')) for s in studies for e in s.get('atomic_effects',[])),
                      **data['graph'].get('corpus_epistemic_and_fidelity_diagnostics',{})},
            'layer1':{'audit_ledger':ledger,'ranked_architectures':data['rankings'].get('ranked_candidate_architectures',[]),
                      'sensitivity':data['rankings'].get('direct_output_logit_sensitivity_and_variance_decomposition',{}),
                      'protocol':strip_trajectories(data['protocol'])},
            'layer2':{'summary':data['layer2'].get('summary',{}),'release_pairing':data['layer2'].get('ensembl_gencode_release_pairing'),
                      'constructs':construct.get('constructs',[]),'selected_gencode_models':selected_gencode_models,'regulatory_audits':construct.get('regulatory_element_audits',{}),
                      'profile_count':len(profiles),'profiles':profiles,'kinetics':stage.get('hour0_hlae_inducible_vs_constitutive_kinetics_note',{}),
                      'model_type':stage.get('model_type')},
            'layer3':{'summary':data['layer3'].get('summary',{}),'complexes':data['layer3'].get('complexes',[])},
            'layer4':{'summary':data['layer4'].get('summary',{}),'lead':strip_trajectories(lead),
                      'geometry':data['layer4'].get('part_a_exact_oxygen_boundary_value_and_volume_algebra',{}),
                      'shrinking_core':strip_trajectories(pde.get('shrinking_core_coupled_bvp_safety_comparison_ultracode_v80_09',{})),
                      'grid_convergence':strip_trajectories(pde.get('spatial_grid_convergence_51_101_201_401_preastra_a094',{})),
                      'joint_vascular_iip_0002':strip_trajectories(pde.get('joint_vascular_timing_and_exchange_sensitivity_iip_0002',{})),
                      'challenge':ode.get('low_glucose_challenge_below_4p2mM_oat1_verification_ultracode_v80_15',{}),
                      'reference_architectures':strip_trajectories(ode.get('graft_architectures_compared',[])),
                      'matched_pde_ode':strip_trajectories(ode.get('matched_pde_ode_qualification_scenarios',[]))},
            'layer5':{'summary':data['layer5'].get('summary',{}),'patient':strip_trajectories(patient)}}


def f(value, digits=2):
    if isinstance(value, bool) or (isinstance(value, (int, float)) and not math.isfinite(value)):
        raise ValueError('Rendered numeric facts must be finite numbers, never booleans')
    return 'unavailable' if value is None else f'{value:,.{digits}f}' if isinstance(value,(int,float)) else str(value)


def narratives(facts):
    c=facts['corpus']; l1=facts['layer1']; l2=facts['layer2']; l4=facts['layer4']; l5=facts['layer5']
    a=l1.get('audit_ledger',{}); registry=a.get('registry',{}); search=a.get('search',{})
    vd=l1['sensitivity'].get('uniform_grid_pooled_pediatric_gate_logit_variance_decomposition',{})
    structural_share=vd.get('structural_share_of_component_variance_sum')
    structural_percent=100*structural_share if type(structural_share) in (int,float) else structural_share
    profiles={p['genotype_id']:p for p in l2['profiles']}
    wt=profiles.get('WT_UNEDITED_CONTROL',{}); arm1a=profiles.get('SCENARIO_A_ARM_1A_BICISTRONIC_P2A_COMPARATOR',{}); arm1b=profiles.get('SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS',{})
    cd47=arm1b.get('hour0_basal_cd47_molecules_per_cell'); wt_cd47=wt.get('hour0_basal_cd47_molecules_per_cell')
    resting_hlae=arm1b.get('hour0_resting_unstimulated_hlae_molecules_per_cell')
    resting_relation='unavailable relative to' if resting_hlae is None else ('below' if resting_hlae < 15000 else 'at or above')
    ratio=cd47/wt_cd47 if cd47 is not None and wt_cd47 else None
    stage6=next((p.get('stage6_neurog3_ep_fraction') for p in l2['profiles'] if 'TYK2' in p.get('genotype_id','') and 'KO' in p.get('genotype_id','')),None)
    wt6=wt.get('stage6_neurog3_ep_fraction'); change=(stage6/wt6-1)*100 if stage6 is not None and wt6 else None
    shell=l4['lead'].get('wigner_seitz_120um_cluster_pericellular_shell_audit_v80_001',{})
    p=l5['patient']; s=l5['summary']; lead=l4['lead']
    shrinking=l4['shrinking_core'].get('coupled_shrinking_core_regimes',[])
    shrink_text='; '.join(str(row.get('scenario_id',row.get('regime_id','shrinking-core scenario')))+': '+f(row.get('min_slab_po2_across_14_days_mmhg'),2)+' mmHg at '+f(row.get('time_of_min_slab_po2_hours'),0)+' h' for row in shrinking)
    challenge=l4['challenge']
    reference_ode_text='; '.join(str(row.get('graft_id'))+': minimum '+f(row.get('min_glucose_mM'),2)+' mM / insulin shutoff t50 '+f(row.get('insulin_secretion_shutoff_half_time_min'),2)+' min' for row in l4['reference_architectures'])
    challenge_text='Low-glucose challenge reference-ODE minima: lead '+f(challenge.get('low_glucose_challenge_lead_with_lt_4p2mM_clamp_min_glucose_0p25min_mM'),4)+', composite-reflex OAT1 '+f(challenge.get('low_glucose_challenge_oat1_without_lt_4p2mM_clamp_min_glucose_0p25min_mM'),4)+', pure-beta-clamp OAT1A '+f(challenge.get('low_glucose_challenge_oat1a_beta_clamp_only_min_glucose_0p25min_mM'),4)+' mM; these share the specified prior kinetic constants and are not measured clinical protection.'
    return {
      'reference_ode':reference_ode_text,
      'challenge':challenge_text,
      'shrinking_core':shrink_text,
      'corpus':f"Current corpus: {c['studies']:,} studies, {c.get('active_deduplicated_studies_count','unavailable')} active deduplicated, {c.get('active_empirical_studies_count','unavailable')} active empirical, {c.get('ring1_corpus_count_reconciliation',{}).get('strict_active_peer_reviewed_empirical_studies_total','unavailable')} strict peer-reviewed empirical; {c['atomic_effects']:,} effects, {c['arms']} arms, {c['reference_control_arms']} reference controls, {c['raw_edit_entries']} edit entries and {c.get('expanded_gene_occurrences')} expanded gene occurrences under the ranker expander. Direction/penalty labels are rule-based heuristics; {c['heuristic_penalty_effects']} flagged penalties are not an independently adjudicated accuracy dataset.",
      'layer1':f"Layer 1 evaluates {f(search.get('all_stacks'),0)} prior-anchored stacks ({f(search.get('feasible_stacks'),0)} feasible, {f(search.get('top_tier_stacks'),0)} code-defined top tier). Registry: {f(registry.get('source_literal_offset_terms'),0)} source terms / {f(registry.get('net_gate_offset_entries'),0)} net offsets / {f(registry.get('auxiliary_scalar_leaves'),0)} auxiliary scalar leaves / {f(registry.get('engineering_priors'),0)} engineering priors. Structural share of component variance is {f(structural_percent,4)}% from unrounded latent logits. Current utility-target profiles with typed numeric LRR evidence: {f(a.get('profile_counts',{}).get('candidate_target_profiles_with_typed_numeric_lrr'),0)}; the shrinkage branch is checked separately against an actual HNF1A-KO source contrast when this count is zero. Seeded direct and shared-parameter Monte Carlo are exploratory sensitivity simulations with finite-draw Wilson intervals and rank quantiles, not clinical probabilities or fitted posterior intervals. CHK-12 recomputes and compares the complete current replay.",
      'layer2':f"Layer 2 uses the pinned release pairing {l2.get('release_pairing')} and contains {l2['profile_count']} literature-anchored stage-specific prior profiles with 1,000 independent prior-scenario draws per endpoint; it is not a fitted dynamic GRN/ODE or held-out validated differentiation predictor. Arm 1B constitutive EEF1A1 drives antisense CD47 and inducible sense 5xISRE drives HLA-E. Computed Arm 1A/1B Hour-0 CD47 is {f(arm1a.get('hour0_basal_cd47_molecules_per_cell'),0)}/{f(cd47,0)} molecules/cell; Arm 1B/WT is {f(ratio,8)}x. Resting Arm 1B HLA-E is {f(arm1b.get('hour0_resting_unstimulated_hlae_molecules_per_cell'),0)} and IFNg24h HLA-E {f(arm1b.get('ifng_24h_induced_hlae_molecules_per_cell'),0)}; the resting endpoint is {resting_relation} the assumed 15,000 threshold. Exact edited beta/alpha/delta/iEC density is unmeasured, with a preregistered Stage-0 acceptance target for Hour-0 density. TYK2-KO Stage-6 prior is {f(stage6,3)} versus {f(wt6,3)} WT ({f(change,1)}%); this is distinct from the transcript-level endpoint. {CANDIDATE}; anti-silencing function requires head-to-head lineage qualification against the published CBX3-side fragment. The local 3,118-row AVI table is quarantined in the private cache: end-to-end batch generation provenance for the full 3,118-row table is incomplete (as a whole it is not a certified model scan), a limited 28 Oct 2026 live re-run against alphagenome v0.10.0 verified 59/59 sampled rows (9 candidate edits plus 50 random rows), AlphaGenome non-commercial output terms conflict with CC-BY-4.0 public redistribution, and zero AlphaGenome outputs are shipped publicly or used as active ranking weights or insertion-safety certifications. Scenario-B P2A profiles are unsupported provisional comparators; the exact 2A-free Scenario-B donor/prediction remains unavailable.",
      'layer3':STRUCTURE1+' '+STRUCTURE2+' The 21-residue linker-2 path includes unresolved HLA-E Gly1 between the 20-aa (G4S)4 linker and HLA-E Ser2 N; 73.50 A contour, 0.4234 extension, 9.901 pN are coarse WLC quantities. The 3CDG reference-crystal Arg107 C-alpha distance is 17.02 A, nearest-heavy-atom distance 14.82 A; neither is a measured engineered Gly107 geometry. Foldseek is a cached-hit replay with missing historical database/tool/date provenance, not a reproducible original search. 8H4X was superseded by 8ZCR (Yan, Zhou et al., 2024, Acta Crystallogr F 80:286-293; DOI:10.1107/S2053230X24009439; PMID:39382088). Exact-product function remains unresolved.',
      'layer4':f"Layer 4 is a conditional model, not validated graft survival. The physical source lamina half-width fraction is {f(l4['geometry'].get('selected_source_lamina_half_width_fraction'),2)}; its continuous central volume is derived from that same geometry and control-volume overlap at every mesh. Nominal minimum slab oxygen is {f(lead.get('min_slab_po2_across_14_days_mmhg'),2)} mmHg at {f(lead.get('time_of_min_slab_po2_hours'),0)} h. The {f(lead.get('min_120um_cluster_core_po2_mmhg'),2)}-mmHg homogenized fixed-surface Dirichlet core is an upper-bound screen; operative uniform-supply / conservative gel-only two-zone cores are {f(shell.get('two_zone_fv_uniform_supply_120um_core_po2_mmhg'),2)}/{f(shell.get('two_zone_fv_gel_only_supply_120um_core_po2_mmhg'),2)} mmHg. Physical supply topology remains a Stage-2 assay gate; mesh error is evaluated on unrounded values, with NUMERICALLY_UNRESOLVED whenever an uncertainty interval straddles the threshold. Reference fixed-functional-beta-IEQ ODE curves have pde_coupling_applied=false and cannot inherit a PDE-derived full-product verdict; matched PDE/ODE states are identified separately. Reference ODE outputs: {reference_ode_text}. {challenge_text} Nominal modeled peak peroxide is {f(lead.get('max_interstitial_h2o2_uM'),3)} uM and maximum pH {f(lead.get('max_interstitial_ph'),3)}.",
      'layer5':f"Layer 5 is an illustrative assigned-rule population simulator, not a clinical forecast. Personalized analytic expectation is {f(p.get('global_analytic_expected_5yr_survival_personalized_compiled_pct'),4)}%, standardized Plan A analytic {f(p.get('global_analytic_expected_5yr_survival_standardized_plan_a_pct'),4)}%; sampled personalized {f(p.get('global_5yr_survival_personalized_dna_compiled_openislet_pct'),2)}%. Analytic ancestry gap is {f(s.get('analytic_max_ancestry_disparity_gap_personalized_pct',p.get('analytic_max_ancestry_disparity_gap_personalized_pct_v80_065')),2)} percentage points; sampled gap {f(s.get('sampled_max_ancestry_disparity_gap_personalized_pct',p.get('max_ancestry_disparity_gap_personalized_pct')),2)} points, with finite sample/seed scope. Genome-wide off-target enumeration is NOT_ENUMERATED and coding-hit count is null; assigned illustrative three-mismatch panel scores are not exhaustive genome/exome scans or empirical HiFi-SpCas9 specificity. Variant identity is source-verified dbSNP interval identity; rare-variant screening frequencies are explicit assumed priors, with gnomAD row-level AC/AN validation pending. Manufacturing yield is null/NOT_MEASURED; exp(-lambda_K*max(0,K-4))*gamma_B**B is a dimensionless engineering utility discount. The executable optional TYK2 design is PEmax+epegRNA. CGBE remains an unimplemented contingency requiring an engineered SpG/SpRY-compatible PAM, editing-window/bystander ledger and separate validation; no CGBE guide is qualified here.",
      'design':f"HIP-Core executes all {len(l1['protocol'].get('runs',[]))} coded rows, including all 11 repeat/bridge/lead rows, across 25 distinct genotypes. The quoted covariance/power applies to 108 independent run-nested clones or 75 shared-genotype clones; 61-clone power is UNEVALUATED with null SE/MDE. Current factor sums are {list(l1['protocol'].get('column_sums_across_36_runs_by_factor',{}).values())}; +1 counts are {list(l1['protocol'].get('plus_one_counts_across_36_runs_by_factor',{}).values())}."
    }


def facts_block(facts):
    n=narratives(facts); lines=[START,'','**Current executable facts (generated from on-disk artifacts; physical gates pending).**','',n['corpus'],'',n['layer1'],'',n['design'],'','| Construct | Donor bp | Packaged bp, including 2 x 145-bp ITRs |','| :--- | ---: | ---: |']
    for construct in facts['layer2']['constructs']:
        lines.append(f"| `{construct['construct_id']}` | {construct.get('total_length_bp')} | {construct.get('packaged_aav6_genome_with_2x145bp_itrs_bp')} |")
    promoter=facts['layer1'].get('audit_ledger',{}).get('native_eef1a1_promoter',{})
    lines += ['',f"Native EEF1A1 core: {promoter.get('length_bp')} bp, {promoter.get('cpg_dinucleotides')} CpG dinucleotides, {promoter.get('total_exact_or_rc_occurrences')} exact/RC donor occurrences; these are native CpG-containing promoter sequences.",'',f"Release manifest: `{facts['manifest_pin_count']:,}` pinned files. Machine facts: `data/extracted/documentation_facts_v0.8.0.json`. These facts supersede stale numeric summaries in explicitly historical rationale sections. Outward rankings remain unfrozen (`outward_rankings_frozen: false`); no exact product has passed physical qualification.",'',END]
    return '\n'.join(lines)


def _refine_table_cell(cell, key, current):
    """Update a bounded generated qualifier while retaining authored cell content."""
    start = '<!-- OPENISLET_CELL_' + key + '_START -->'
    end = '<!-- OPENISLET_CELL_' + key + '_END -->'
    note = start + ' Current generated scope: ' + current + ' ' + end
    _validate_markers(cell, start, end)
    if start in cell:
        return re.sub(re.escape(start) + r'.*?' + re.escape(end), lambda _m: note, cell)
    if current.strip() in cell:
        return cell
    # Earlier generated cells had no marker. If only numbers changed, update
    # those values in place rather than appending a second, contradictory copy.
    numeric_shape=lambda value: re.sub(r'(?<![A-Za-z])[-+]?\d+(?:[.,]\d+)*(?:[eE][-+]?\d+)?','#',value.strip())
    if numeric_shape(cell)==numeric_shape(current):
        return cell[:len(cell)-len(cell.lstrip())]+current.strip()+cell[len(cell.rstrip()):]
    return cell.rstrip() + '<br>' + note + ' '


def _score_cell(cell, value, floor=None):
    number = '`' + f(value, 4) + '`'
    if re.search(r'`[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)`', cell):
        cell = re.sub(r'`[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)`', lambda _m: number, cell, count=1)
    else:
        cell = ' **' + number + '** ' + cell.strip()
    if floor is not None:
        cell = re.sub(r'(Floor:\s*`)[0-9.]+(`)', lambda m: m.group(1) + f(floor,4) + m.group(2), cell)
    return cell


def _ranking_header(cells):
    patterns = {'pediatric':r'pediatric','adult':r'adult','composite':r'composite',
                'rank':r'grid\s*position|\brank\b','discount':r'mfg|manufacturing|engineering.*discount',
                'limiting':r'limiting|bottleneck','specification':r'specification'}
    matches = {key: [i for i,cell in enumerate(cells) if re.search(pattern,cell,re.I)] for key,pattern in patterns.items()}
    required = ('pediatric','adult')
    intent=any(re.search(r'index|score|composite|grid\s*position|architecture\s*id',cell,re.I) for cell in cells)
    if not intent: return None
    if any(matches[key] for key in required) and (any(len(matches[key]) != 1 for key in required) or len({matches[key][0] for key in required}) != 2 or len(matches['composite']) > 1):
        raise ValueError('Ranking table requires distinct, unique pediatric/adult columns and at most one composite column')
    mapping = {key: indices[0] if indices else None for key,indices in matches.items()}
    return mapping if all(mapping[key] is not None for key in ('pediatric','adult')) else None


def _validate_markers(text, start, end):
    if text.count(start) != text.count(end) or text.count(start) > 1 or (start in text and text.index(start) > text.index(end)):
        raise ValueError('Malformed or duplicate generated-document markers')


def _table_cells(line):
    # A pipe preceded by an even number of backslashes is a delimiter.
    cells=[];start=0;slashes=0
    for index,char in enumerate(line):
        if char=='|' and slashes%2==0:
            cells.append(line[start:index]);start=index+1
        slashes=slashes+1 if char=='\\' else 0
    cells.append(line[start:])
    return cells


def _protect_link_destinations(text, prefix):
    """Shield link destinations during prose substitutions without a parser dependency.

    Balanced parentheses and backslash escapes are retained byte-for-byte;
    reference definitions and autolinks are protected as complete literals.
    The surrounding authored label remains available to ordinary rendering.
    """
    protected = []
    def hold(value):
        key = f'{prefix}LINK_{len(protected)}\x00'
        protected.append((key, value))
        return key
    # A definition's destination and title may each begin on the next line.
    # Scan escapes and balanced destination parentheses before accepting it;
    # an ordinary following paragraph must remain available to rendering.
    definition = re.compile(r'^ {0,3}\[(?:\\[^\r\n]|[^\[\]\\\r\n])+\]:[ \t]*', re.M)
    def next_content(position):
        return position + len(re.match(r'[ \t]*(?:\r?\n[ \t]*)?', text[position:])[0])
    spans = []
    for match in definition.finditer(text):
        if spans and match.start() < spans[-1][1]:
            continue
        start = end = next_content(match.end())
        if end >= len(text):
            continue
        angle = text[end] == '<'
        if angle:
            end += 1
        depth = 0
        while end < len(text):
            char = text[end]
            if char == '\\' and end+1 < len(text):
                end += 2; continue
            if angle:
                if char == '>':
                    end += 1; break
                if char in '<\r\n':
                    break
            else:
                if char.isspace() or char in '<>':
                    break
                if char == '(':
                    depth += 1
                elif char == ')':
                    depth -= 1
                    if depth < 0: break
            end += 1
        if end == start or depth or (angle and text[end-1:end] != '>'):
            continue
        destination_end = end
        title_start = next_content(end)
        if title_start > end and text[title_start:title_start+1] in ('"', "'", '('):
            closing = ')' if text[title_start] == '(' else text[title_start]
            end = title_start+1
            while end < len(text):
                if text[end] == '\\' and end+1 < len(text):
                    end += 2; continue
                if text[end] == closing:
                    end += 1; break
                if text[end] == '\n' and re.match(r'[ \t]*\r?\n', text[end+1:]):
                    break
                end += 1
            if text[end-1:end] != closing:
                end = destination_end
        else:
            end = destination_end
        end += len(re.match(r'[ \t]*', text[end:])[0])
        if end < len(text) and text[end] not in '\r\n':
            # An invalid title on the next line does not invalidate an otherwise
            # complete destination on the preceding line.
            end = destination_end + len(re.match(r'[ \t]*', text[destination_end:])[0])
            if end < len(text) and text[end] not in '\r\n':
                continue
        spans.append((match.start(), end))
    for start, end in reversed(spans):
        text = text[:start] + hold(text[start:end]) + text[end:]
    # Preserve lexical HTML, including quoted attributes and comments. Managed
    # table delimiters remain visible to their existing scope-update parser.
    html = re.compile(r'''<!--.*?-->|<\?.*?\?>|<!\[CDATA\[.*?\]\]>|<![A-Z][^>]*>|</?[A-Za-z][A-Za-z0-9-]*(?:\s+[A-Za-z_:][A-Za-z0-9_.:-]*(?:\s*=\s*(?:"[^"]*"|'[^']*'|[^\s"'=<>`]+))?)*\s*/?>''', re.S)
    def hold_html(match):
        value = match[0]
        if re.fullmatch(r'<!-- OPENISLET_CELL_[A-Z0-9_]+_(?:START|END) -->', value):
            return value
        return hold(value)
    text = html.sub(hold_html, text)
    pieces = []; start = 0; i = 0
    while i < len(text)-1:
        if text[i:i+2] != '](':
            i += 1; continue
        end = i+2; depth = 1; angle = False; quoted = None
        while end < len(text) and depth:
            char = text[end]
            if char == '\\':
                end += 2; continue
            if quoted:
                if char == quoted: quoted = None
            elif char in ('"', "'") and (end == i+2 or text[end-1].isspace()): quoted = char
            elif char == '<': angle = True
            elif char == '>': angle = False
            elif not angle and char == '(': depth += 1
            elif not angle and char == ')': depth -= 1
            end += 1
        if depth:
            i += 2; continue
        pieces.append(text[start:i+1]); pieces.append(hold(text[i+1:end]))
        start = end; i = end
    pieces.append(text[start:])
    text = ''.join(pieces)
    text = re.sub(r'<(?:[A-Za-z][A-Za-z0-9+.-]*:[^<>\n]*|[^<>\s]+@[^<>\s]+)>', lambda m: hold(m[0]), text)
    return text, protected


def render_text(text, facts, filename):
    # Fenced examples and diagrams are authored literals, not replacement inputs.
    fences=[];unclosed_key=None
    literal_prefix='\x00OPENISLET_LITERAL_FENCE_'
    while literal_prefix in text:
        literal_prefix = literal_prefix + '_'
    def literal_key(index):return f'{literal_prefix}{index}\x00'
    lines=text.splitlines(keepends=True); held=[]; i=0
    while i < len(lines):
        # Also shield fenced blocks inside block quotes and ordinary indented
        # code. Conservative shielding retains authored literals byte-for-byte.
        prefix=re.match(r'^(?: {0,3}>[ \t]?| {0,3}(?:[-+*]|[0-9]{1,9}[.)])[ \t]+)*',lines[i]).group(0)
        candidate=lines[i][len(prefix):]
        match=re.match(r'^ {0,3}(`{3,}|~{3,})([^\n]*)',candidate)
        raw_html=re.match(r'^ {0,3}<(script|pre|style|textarea)(?:\s|>|$)',candidate,re.I)
        if raw_html:
            start=i;closing_html=re.compile(r'</'+raw_html[1]+r'\s*>',re.I)
            while i<len(lines):
                line=lines[i];i+=1
                if closing_html.search(line):break
            key=literal_key(len(fences))
            fences.append((key,''.join(lines[start:i]).rstrip('\n')));held.append(key+'\n')
            if not closing_html.search(line):unclosed_key=key
            continue
        if not match and candidate.strip() and (candidate.startswith('    ') or candidate.startswith('\t')):
            start=i;i+=1
            while i<len(lines) and (not lines[i].strip() or re.sub(r'^(?: {0,3}>[ \t]?)*','',lines[i]).startswith(('    ','\t'))):i+=1
            # Separator blank lines belong to surrounding Markdown. Keeping
            # them outside the literal avoids consuming one on each re-render.
            while i>start+1 and not lines[i-1].strip():i-=1
            key=literal_key(len(fences))
            fences.append((key,''.join(lines[start:i]).rstrip('\n')));held.append(key+'\n')
            continue
        if not match or (match.group(1)[0]=='`' and '`' in match.group(2)):
            held.append(lines[i]); i+=1; continue
        fence=match.group(1); start=i; i+=1
        closed=False
        # Container order matters: '- >' continues as '  >', whereas '> -'
        # continues as '>  '. Stripping all quotes first conflates those forms.
        containers=[]
        for item in re.finditer(r' {0,3}>[ \t]?| {0,3}(?:[-+*]|[0-9]{1,9}[.)])[ \t]+',prefix):
            containers.append(('quote' if '>' in item[0] else 'list',len(item[0].expandtabs(4))))
        closing=re.compile(r'^ {0,3}'+re.escape(fence[0])+'{'+str(len(fence))+r',}[ \t]*(?:\r?\n)?$')
        while i<len(lines):
            inner=lines[i]
            for kind,width in containers:
                if kind=='quote':
                    quote=re.match(r'^ {0,3}>[ \t]?',inner)
                    if quote is None:
                        inner=None;break
                    inner=inner[len(quote[0]):]
                else:
                    expanded=inner.expandtabs(4)
                    if inner.strip() and not expanded.startswith(' '*width):
                        inner=None;break
                    inner=expanded[width:]
            if inner is None:
                # The enclosing list/quote ended, even without a fence closer.
                closed=True;break
            i+=1
            if closing.fullmatch(inner):closed=True;break
        key=literal_key(len(fences))
        fences.append((key,''.join(lines[start:i]).rstrip('\n')))
        if not closed:unclosed_key=key
        held.append(key+'\n')
    text=''.join(held)
    _validate_markers(text,START,END)
    _validate_markers(text,'<!-- OPENISLET_CURRENT_LAYER_SCOPE -->','<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->')
    n=narratives(facts)
    text=re.sub(re.escape(START)+r'.*?'+re.escape(END)+r'\n*','',text,flags=re.S)
    text=re.sub(r'<!-- OPENISLET_CURRENT_LAYER_SCOPE -->.*?<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->\s*','',text,flags=re.S)
    text,protected_links=_protect_link_destinations(text,literal_prefix)
    replacements={
      'WAVE-0 PUBLICATION & RANKING FREEZE (Active':'WAVE-0 PRIOR-ANCHORED PUBLICATION (outward_rankings_frozen: false; physical assay gates pending',
      'while guaranteeing Hour-0 density':'with a preregistered Stage-0 acceptance target for Hour-0 density',
      'GSG-tPT2A':'GSG-P2A',
      '2A-free design preventing':'proposed 2A-free topology (exact donor and benefit unresolved), intended to reduce',
      '100% native Phe19 N-terminus':'predicted native Phe19 N-terminus; processing requires direct qualification',
      'including `Ile22`':'including the unresolved HLA-E Gly1 between the 20-aa (G4S)4 linker and HLA-E Ser2 N',
      '§ 356k(e)` governs review actions on the designation request':'§ 356k(e)` concerns actions to expedite development and review of an application using a designated platform; requests and determinations are governed by § 356k(c)-(d)',
      'Pan-Allelic `SIRPα` Proof':'Reference SIRPA V1/V2 sequence and geometry comparison; pan-allelic protection unverified',
      'Pan-Allelic SIRPa Proof':'Reference SIRPA V1/V2 sequence and geometry comparison; pan-allelic protection unverified',
      '1` `Asp101` insertion in `V2`':'1` `Asp101` insertion in `V1` (4CMM), relative to V2 (2JJS)',
      '[14, 20, 22, 24, 27, 45, 65, 66, 70, 77, 79, 100]':'[14, 20, 22, 24, 27, 45, 65, 66, 70, 77, 79, 101] (2JJS/V2 numbering)',
      'not CpG-free':'CpG-containing', 'CpG-free':'CpG-containing candidate; anti-silencing pending',
      'To guarantee zero false negatives':'To benchmark retrieval rather than guarantee zero false negatives',
      'live NCBI `dbSNP`/`ClinVar`/`gnomAD v4.1` polymorphism audit':'source-verified dbSNP/ClinVar interval identities plus explicitly assumed rare-variant screening priors; gnomAD row-level frequency validation pending',
      'v113 / GENCODE v46':'112 / GENCODE v46',
      'canonical Ensembl v113':'pinned Ensembl/GENCODE (current versions in documentation_facts_v0.8.0.json)',
      'source-adjudicated':'source-bound heuristic', 'Source-adjudicated':'Source-bound heuristic',
      'MVDGTLLLL':'MSRSVALAVLALLSLSGLEA (20 aa encoded B2M signal; processing qualification pending)',
    }
    for old,new in replacements.items(): text=text.replace(old,new)
    # Never consume code-span delimiters while replacing an alias. Earlier
    # permissive backtick matching left valid inline code with an open span.
    text=re.sub(r'(?:minimal[ *]+)?455[ *]*bp[ *]*(?:minimal[ *]+|core[ *]+)?(?:[\[(]or[^\])]+[\])][ *]*)?CBX3-A2UCOE(?:\s+core)?',CANDIDATE+' (legacy CBX3-A2UCOE alias; anti-silencing function pending)',text)
    text=text.replace('AAVS1-[455bp CBX3-A2UCOE]','AAVS1-[HNRNPA2B1-proximal 455bp candidate; legacy CBX3-A2UCOE]')
    text=text.replace('mPGK \\| EEF1A1-core -> HLA-E','EEF1A1-core \\| Syn_5xISRE_minP -> HLA-E')
    text=text.replace('mPGK | EEF1A1-core -> HLA-E','EEF1A1-core | Syn_5xISRE_minP -> HLA-E')
    output=[]
    ranking_columns=None
    source_lines=text.splitlines()
    table_width=None
    for row_index,line in enumerate(source_lines):
        original=line
        if line.startswith('|'):
            next_line=source_lines[row_index+1] if row_index+1 < len(source_lines) else ''
            if next_line.startswith('|') and re.fullmatch(r'[| :\-\t]+',next_line):
                header_cells=_table_cells(line)
                ranking_columns=_ranking_header(header_cells)
                table_width=len(header_cells)
        else: ranking_columns=None;table_width=None
        # Current count assertions must follow regenerated records, including
        # review reclassification; do not freeze values from a prior audit.
        for label, key in [('total arms', 'arms'), ('arms', 'arms'), ('atomic effects', 'atomic_effects'), ('atomic quantitative & qualitative effects', 'atomic_effects')]:
            line=re.sub(r'`[\d,]+`(?= '+label+r'\b)', '`'+f"{facts['corpus'][key]:,}"+'`', line)
        line=re.sub(r'`[\d,]+ / [\d,]+ = 100\.0%`(?= exact)', '`'+f"{facts['corpus']['atomic_effects']:,} / {facts['corpus']['atomic_effects']:,} = 100.0%"+'`', line)
        if 'single-field substring match' in line and 'effects' in line:
            line=re.sub(r'`[\d,]+ / [\d,]+`', '`'+f"{facts['corpus']['atomic_effects']:,} / {facts['corpus']['atomic_effects']:,}"+'`', line)
        strict_empirical=facts['corpus'].get('ring1_corpus_count_reconciliation',{}).get('strict_active_peer_reviewed_empirical_studies_total')
        if strict_empirical is not None:
            line=re.sub(r'`[\d,]+`(?= strict peer-reviewed empirical)', '`'+str(strict_empirical)+'`', line)
        if ('PD-L1' in line or 'Phe19' in line) and '100%' in line and ('cleav' in line or 'N-termin' in line):
            line=line.replace('100%', 'predicted (direct processing qualification pending)')
        if ('HNRNPA2B1' in line or '455' in line or 'AAVS1' in line) and ('iCasp9' in line or 'suicide switch' in line):
            line=line.replace('and tested with', 'and proposed for testing with').replace('insulated by', 'proposed insulation using')
        # Refresh generated labels too: idempotence must not freeze a prior
        # artifact snapshot when the scientific source objects change.
        if line.startswith('- **Validated Ring-1 records:**'):
            line='- **Validated Ring-1 records:** '+n['corpus']+' Private dossier caches are excluded from public redistribution.'
        elif line.lstrip().startswith('- Current corpus:'):
            line=re.match(r'^\s*',line).group(0)+'- '+n['corpus']
        elif line.lstrip().startswith('- **Heuristic extraction labels:**'):
            line='  - **Heuristic extraction labels:** '+n['corpus']
        elif line.startswith('> 1. **Prior-anchored scenario outputs:**'):
            line='> 1. **Prior-anchored scenario outputs:** '+n['layer1']
        elif line.startswith('> **Estimand & prior sensitivity:**'):
            line='> **Estimand & prior sensitivity:** '+n['layer1']
        elif line.startswith('  - Layer 4 is a conditional model'):
            line='  - '+n['layer4']+' Shrinking-core sensitivity: '+n['shrinking_core']+'. Current chemistry extrema are listed in the atlas.'
        elif line.startswith('  - Low-glucose challenge reference-ODE minima:'):
            line='  - '+n['challenge']+' Reference fixed-functional-beta-IEQ ODE curves have pde_coupling_applied=false and cannot inherit a PDE-derived full-product verdict; matched PDE/ODE states are identified separately. Reference ODE outputs: '+n['reference_ode']+'.'
        elif line.startswith('  - Layer 5 is an illustrative assigned-rule population simulator'):
            line='  - '+n['layer5']
        # Replace entire stale assertions, preserving the surrounding document's structure.
        if ('**Layer 1' in line or '1. **Layer 1' in line) and len(line)>400: line='1. **Layer 1:** '+n['corpus']+' '+n['layer1']
        elif ('**Layer 2' in line or '2. **Layer 2' in line) and len(line)>500: line='2. **Layer 2:** '+n['layer2']+' '+CHS4+'. '+P2A+'. '+DELIVERY
        elif ('**Layer 3' in line or '3. **Layer 3' in line) and len(line)>500: line='3. **Layer 3:** '+n['layer3']
        elif ('**Layer 4' in line or '4. **Layer 4' in line) and len(line)>500: line='4. **Layer 4:** '+n['layer4']
        elif ('**Layer 5' in line or '5. **Layer 5' in line) and len(line)>500: line='5. **Layer 5:** '+n['layer5']
        elif '**Estimand & Prior-Dominance Note:**' in line or ('99.' in line and ('variance decomposition' in line or 'structural prior variance' in line)): line='> **Estimand & prior sensitivity:** '+n['layer1']
        elif line.startswith('> 1. **Prior-Anchored Scenario Outputs'): line='> 1. **Prior-anchored scenario outputs:** '+n['layer1']
        elif '**Structural Compatibility Screen for `GY-008`' in line: line='  - '+STRUCTURE1
        elif '**Structural Compatibility Screen for `GY-003`' in line: line='  - '+STRUCTURE2
        elif line.startswith('  - Across the 336-hour'): line='  - '+n['layer4']+' Shrinking-core sensitivity: '+n['shrinking_core']+'. Current chemistry extrema are listed in the atlas.'
        elif line.startswith('  - Executes `compile_personalized_patient_blueprint`'): line='  - '+n['layer5']
        elif ('3,118' in line or '3118' in line) and ('3,092' in line or '26` rows are primary' in line or 'unresolved-origin' in line or 'quarantined in the private cache' in line): line='  - The 3,118-row AVI table is quarantined in the private cache: end-to-end batch generation provenance for the full 3,118-row table is incomplete (as a whole it is not a certified model scan), a limited 28 Oct 2026 live re-run against alphagenome v0.10.0 verified 59/59 sampled rows (9 candidate edits plus 50 random rows), AlphaGenome non-commercial output terms conflict with CC-BY-4.0 public redistribution, and zero AlphaGenome outputs are shipped publicly or used as active ranking weights or insertion-safety certifications. Ensembl release 112 / GENCODE v46 selected transcript identities and coordinates are read from the pinned raw models.'
        elif ('1,000' in line or '1000' in line) and ('credible' in line or 'GRN' in line) and len(line)>200 and not line.startswith('- **Stage 0→7 prior-scenario plane:**'): line='  - '+n['layer2']
        if 'EEF1A1' in line and re.search(r'`?20`?\s+(?:native\s+)?CpG',line):
            count=facts['layer1'].get('audit_ledger',{}).get('native_eef1a1_promoter',{}).get('cpg_dinucleotides')
            line=re.sub(r'`?20`?\s+(native\s+)?CpG(?:\s+dinucleotides|s)?',lambda m: f"`{count}` {m.group(1) or ''}CpG dinucleotides",line)
        if ('`44` adjudicated' in line and 'independently' not in line) or ('false positives' in line and ('`0` false positives' in line or '0` false positives' in line or '39` strict true' in line)) or ('44/44' in line and 'penalt' in line):
            if line.startswith('|') and not line.startswith('| **Stage 0**'):
                cells=_table_cells(line)
                if len(cells)>=5:
                    cells[2]=_refine_table_cell(cells[2],'HEURISTIC_COUNT',str(facts['corpus']['heuristic_penalty_effects'])+' heuristic penalty flags; independent accuracy unmeasured')
                    cells[3]=_refine_table_cell(cells[3],'HEURISTIC_SCOPE','Rule-based null/clause/subject negative controls are regression checks, not an independently adjudicated sensitivity or false-positive estimate.')
                    line='|'.join(cells)
            else: line='  - **Heuristic extraction labels:** '+n['corpus']
        if line.startswith(('- **Stage 0→7 Perturbation Trajectory Plane:**','- **Stage 0→7 prior-scenario plane:**')):
            line='- **Stage 0→7 prior-scenario plane:** '+n['layer2']
        if '3,092' in line and ('effects' in line or 'Atomic Effects' in line): line=line.replace('3,092',f"{facts['corpus']['atomic_effects']:,}")
        if '561` active' in line: line=line.replace('561` active',str(facts['corpus'].get('active_deduplicated_studies_count'))+'` active')
        if '233/233' in line: line=line.replace('233/233',str(facts['corpus']['raw_edit_entries'])+'/'+str(facts['corpus']['raw_edit_entries']))
        if '312` gene tokens' in line: line=line.replace('312` gene tokens',str(facts['corpus'].get('expanded_gene_occurrences'))+'` expanded gene occurrences')
        if line.startswith('| **Prior-Constant & Auxiliary'):
            cells=_table_cells(line);cells[2]=_refine_table_cell(cells[2],'PRIOR_VALUES',n['layer1']);cells[3]=_refine_table_cell(cells[3],'PRIOR_SCOPE','All active registry ranges and source/omics channels are replayed; independent source validity remains a separate question.');line='|'.join(cells)
        if '7` functional LRR pairs' in line: line=line.replace('7` functional LRR pairs',str(facts['corpus'].get('source_bound_functional_intervention_lrr_pairs_count', 4))+'` functional intervention LRR pairs ('+str(facts['corpus']['total_numeric_lrr_effects'])+'` total typed numeric LRR effects)')
        if line.startswith('  - Holding peripheral systemic insulin-action'):
            line='  - '+n['challenge']+' Reference fixed-functional-beta-IEQ ODE curves have pde_coupling_applied=false and cannot inherit a PDE-derived full-product verdict; matched PDE/ODE states are identified separately. Reference ODE outputs: '+n['reference_ode']+'.'
        if 'RNLS' in line and ('poor dynamic GSIS' in line or 'Fig. 6' in line or 'GSIS/OCR tradeoff' in line):
            if line.startswith('| **`GRAVEYARD_004`'):
                cells=_table_cells(line)
                cells[-2]=_refine_table_cell(cells[-2],'RNLS_SCOPE',RNLS)
                line='|'.join(cells)
            elif line.startswith('|'):
                cells=_table_cells(line)
                cells=[(' '+RNLS+' ') if ('poor dynamic GSIS' in cell or ('MacDonald' in cell and '39988068' in cell)) else cell for cell in cells]
                line='|'.join(cells)
            else:
                prefix='> 3. ' if line.startswith('> 3.') else re.match(r'^\s*',line).group(0)+'- '
                line=prefix+'**RNLS precaution:** '+RNLS
        if 'Worm-Like Chain (`WLC`)' in line:
            line=re.sub(r'`HLA-E\*01:03 Gly107` sits.*$', 'The 3CDG reference crystal Arg107 C-alpha is 17.02 A from the receptor; nearest-heavy-atom clearance is 14.82 A. Engineered Gly107 geometry is unmeasured; stability and surface-density values are working priors.',line)
        if line.startswith(('- **Multi-Repository Discovery Manifest','- **PubMed Discovery Manifest')):
            discovery=facts['discovery_manifest']
            line=f"- **PubMed Discovery Manifest (`data/raw/corpus_manifest_all_sources.jsonl`):** {discovery['jsonl_rows']:,} JSONL rows / {discovery['unique_pubmed_pmids']:,} unique PubMed PMIDs, with {discovery['inline_abstract_strings']} inline abstract strings. Separate preprint, patent, GEO and ClinicalTrials.gov inventories are retained in their own raw files and are not additive components of this PMID manifest."
        if line.startswith('| **Raw Corpus Manifest'):
            discovery=facts['discovery_manifest'];cells=_table_cells(line)
            cells[2]=_refine_table_cell(cells[2],'DISCOVERY_COUNT',f"{discovery['jsonl_rows']:,} JSONL rows / {discovery['unique_pubmed_pmids']:,} unique PubMed PMIDs")
            cells[3]=_refine_table_cell(cells[3],'DISCOVERY_SCOPE','PubMed records only; separate preprint/patent/GEO/trial files are independent inventories, not additional rows in this manifest.')
            line='|'.join(cells)
        if line.startswith('All extraction, graph-construction, and combinatorial scenario pipelines'):
            line=f"All extraction, graph-construction and scenario pipelines are hash-pinned in [`manifest.json`](manifest.json) (`{facts['manifest_pin_count']:,}` pinned files). Private research dossiers remain excluded from public bundles; the rights-aware exporter applies the <=240-character / <25%-abstract policy and preserves source/license qualification. Counts by scientific scope are generated above rather than inferred by adding separate inventories."
        if line.startswith('| **Stage 0**'):
            cells=_table_cells(line)
            if len(cells)>=5:
                cells[3]=_refine_table_cell(cells[3],'STAGE0',n['corpus']+' Stage-0 construct, density, differentiation and immune/GSIS acceptance experiments remain preregistered pending.')
                line='|'.join(cells)
        if line.lstrip().startswith('- **`594` experimental') or ('genetic edit entries' in line and '233' in line):
            line='  - '+n['corpus']
        if ('3,092' in line or '3092' in line) and 'atomic' in line and 'effect' in line:
            line='  - Current extracted effects: '+str(facts['corpus']['atomic_effects'])+'; source identifier, digest, typed numeric role and direct character slice are verified together. Inequality P values retain relation/bound with exact p_value=null.'
        if 'schema-validated primary study records' in line:
            line='- **Validated Ring-1 records:** '+n['corpus']+' Private dossier caches are excluded from public redistribution.'
        if '217-bp' in line and 'U78775.1' in line and not line.startswith('Current cHS4'):
            if 'LOCUS4_' in line:
                loc4=next((c for c in facts['layer2']['constructs'] if c['construct_id'].startswith('LOCUS4')), {})
                line='  4. **Locus 4:** '+str(loc4.get('total_length_bp'))+' bp donor / '+str(loc4.get('packaged_aav6_genome_with_2x145bp_itrs_bp'))+' bp including ITRs; '+CHS4+'; '+CANDIDATE+'; 230-bp native EEF1A1 core, iCasp9 and polyA. Exact sequence and source intervals are in the construct ledger.'
            else:
                line='  - **Current cHS4 source:** '+CHS4+'. The earlier 217-bp U78775.1 match is historical.'
        if line.startswith('  - **Pass 1') or line.startswith('  - **Pass 2'):
            line=line.split(' All specificity figures')[0]
            line += ' All specificity figures on this line are assigned illustrative three-mismatch panel scores, not an exhaustive genome/exome scan or HiFi-SpCas9 empirical specificity estimate; genome-wide enumeration is NOT_ENUMERATED.'
        if filename=='README.md' and line.lstrip().startswith('- **Locus 1:**'):
            line='   - **Locus 1:** Arm 1A bicistronic HLA-E-SCT-GSG-P2A-CD47-IgV-GPI is the nonviral comparator; Arm 1B has antisense constitutive EEF1A1-to-CD47 and sense inducible 5xISRE-to-HLA-E. '+DELIVERY+' Exact edited beta/alpha/delta/iEC densities remain unmeasured; QuantiBRITE, matched immune killing and GSIS are preregistered Stage-0 acceptance targets.'
        # Bind named ranking summaries and table values to the current named architecture.
        for arch in facts['layer1']['ranked_architectures']:
            aid=arch['architecture_id']
            if aid not in original: continue
            if original.startswith('|') and ranking_columns is not None:
                columns=_table_cells(line)
                if len(columns) != table_width:
                    raise ValueError('Ranking row width differs from its header')
                pediatric=arch['utility_pediatric_early_onset_t1d'];adult=arch['utility_adult_onset_lada']
                values={'pediatric':pediatric['series_system_joint_survival_probability'],
                        'adult':adult['series_system_joint_survival_probability'],
                        'composite':arch.get('composite_series_system_utility')}
                for key,value in values.items():
                    index=ranking_columns[key]
                    if index is None: continue
                    floor=pediatric.get('heuristic_conservative_floor_index') if key=='pediatric' else None
                    columns[index]=_score_cell(columns[index],value,floor)
                if ranking_columns['rank'] is not None:
                    columns[ranking_columns['rank']]=f" **`#{arch['global_rank_among_all_combinations']} / {facts['layer1']['audit_ledger'].get('search',{}).get('all_stacks')}`** "
                if ranking_columns['discount'] is not None:
                    index=ranking_columns['discount'];columns[index]=_score_cell(columns[index],pediatric.get('manufacturing_locus_penalty_factor'))
                if ranking_columns['limiting'] is not None:
                    columns[ranking_columns['limiting']]=str(pediatric.get('limiting_bottleneck_gate'))+' ('+f(pediatric.get('limiting_bottleneck_score'),4)+')'
                if aid.startswith('ISLET1_CHASSIS_B') and ranking_columns['specification'] is not None:
                    index=ranking_columns['specification']
                    columns[index]=columns[index].split(' Scenario-B scope:')[0]+' Scenario-B scope: unmaterialized six-locus research comparator; exact 2A-free donor/function unavailable. '
                line='|'.join(columns)
            elif ('Scenario Output' in original or 'Reference Profile' in original or 'current global rank #' in original) and original.lstrip().startswith(('1.','2.','-')):
                line=f"- **`{aid}`**: current global rank #{arch['global_rank_among_all_combinations']}, K={arch['genomic_locus_edit_count_K']}, B={arch['base_edit_count_B']}; composite prior-anchored index {f(arch.get('composite_series_system_utility'),4)}. Physical qualification remains pending."
        if original.startswith('|') and (not line.startswith('|') or len(_table_cells(line))!=len(_table_cells(original))):
            raise ValueError('A document substitution would corrupt an authored table row')
        output.append(line)
    text='\n'.join(output).strip()+'\n'
    # Keep old study-design details, but correct active row-allocation and stale counts.
    text=re.sub(r'column sums \$\+2\$ on.*?\$-6\$ on `X5`', 'column sums from the current machine run table '+str(list(facts['layer1']['protocol'].get('column_sums_across_36_runs_by_factor',{}).values())),text)
    text=re.sub(r'`[\d,]+`(?= pinned files)', '`'+f"{facts['manifest_pin_count']:,}"+'`',text)
    text=re.sub(r'`(?:233|312)`(?=\s+(?:edit|expanded|gene-symbol))',lambda m:'`'+str(facts['corpus']['raw_edit_entries'] if m.group(0)=='`233`' else facts['corpus'].get('expanded_gene_occurrences'))+'`',text)
    # A generated footer may itself contain a protected source link. Inspect a
    # read-only expanded view for presence, without exposing URLs to edits.
    presence_text=text
    for key,value in reversed(protected_links):
        presence_text=presence_text.replace(key,value)
    if 'Shy' in text and DELIVERY not in text: text+='\n'+DELIVERY+'\n'
    if 'RNLS' in text and RNLS not in text: text+='\n'+RNLS+'\n'
    if 'cHS4' in text and CHS4 not in text: text+='\nCurrent cHS4 provenance: '+CHS4+'. Earlier 217-bp U78775.1 descriptions are superseded accession matches.\n'
    if '356k' in text and STATUTE not in presence_text: text+='\n'+STATUTE+'\n'
    if 'single 22-aa GSG-P2A' not in text and ('P2A' in text or 'tPT2A' in text): text+='\nCurrent encoded element: '+P2A+'. SCT encodes B2M signal MSRSVALAVLALLSLSGLEA (20 aa); ER cleavage and mature N-terminus remain processing assay gates.\n'
    if 'CGBE1' in text: text=text.replace('or `CGBE1` for', '(CGBE1 is an unimplemented contingency, with no qualified guide) for')
    # A single generated scope paragraph is present in each current companion document.
    # Only an ATX heading can be split safely after one line. In particular a
    # setext underline belongs to its preceding paragraph, including multiple
    # heading lines; prepend generated blocks when there is no ATX heading.
    no_header=re.match(r'^ {0,3}#{1,6}(?:[ \t]+|$)',text.partition('\n')[0]) is None
    scope='\n<!-- OPENISLET_CURRENT_LAYER_SCOPE -->\n'+n['layer2']+'\n\n'+n['layer3']+'\n\n'+n['layer4']+'\n\n'+n['layer5']+'\n<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->\n'
    # An unclosed authored fence consumes the rest of a Markdown document.
    # Put generated scope before it so re-renders never trap generated blocks
    # inside an authored literal and append a second copy on each invocation.
    if unclosed_key is not None:text=text.replace(unclosed_key,scope+'\n'+unclosed_key,1)
    else:text+=scope
    if no_header:first,rest='',text
    else:first,sep,rest=text.partition('\n')
    historical_names = {'ISLET-1-MASTER-SPEC-v0.1.md',
                        'ISLET-1-CANDIDATE-BLUEPRINTS-AND-OPTIMIZER-v0.2.md',
                        'ISLET-1-SYNTHESIS-REPORT-v0.4.md',
                        'PROJECT-ISLET-1-RING1-DISCOVERY-ATLAS-AND-STRATEGY.md'}
    historical_notice = '> **Earlier-revision rationale:** The body originated in an earlier design revision and retains historical reasoning. Unrendered historical values are not current release claims; the generated facts above and current machine-readable artifacts govern active numerical outputs.'
    rest=rest.replace(historical_notice+'\n\n','')
    if filename in historical_names:
        rest=historical_notice+'\n\n'+rest.lstrip()
    rendered=(first+'\n\n' if first else '')+facts_block(facts)+'\n\n'+rest.lstrip()
    if fences or protected_links:
        # Restore after every prose/count substitution. Later protected syntax
        # may contain an earlier token (e.g. HTML inside a link title), so unwind
        # in reverse creation order. The prefix was absent from authored input.
        for key,value in reversed(fences+protected_links):
            if rendered.count(key)!=1:
                raise ValueError('Authored literal rendering lost unique placeholder provenance')
            rendered=rendered.replace(key,value)
        if literal_prefix in rendered:raise ValueError('Unrestored authored literal rendering marker')
    return rendered


def _publish_document_batch(prepared,rollback_journal=None):
    try:from .study_artifact_exports import publish_byte_artifacts
    except ImportError:from study_artifact_exports import publish_byte_artifacts
    publish_byte_artifacts({path:text.encode('utf-8') for path,text in prepared.items()},rollback_journal=rollback_journal)


def write_documentation(root: Path | str, manifest_pin_count: int, render: bool = True, *, rollback_journal=None) -> dict:
    root=Path(root)
    _safe_location(root)
    if isinstance(manifest_pin_count,bool) or not isinstance(manifest_pin_count,int) or manifest_pin_count<0: raise ValueError('manifest_pin_count must be a nonnegative integer')
    if type(render) is not bool:raise ValueError('render must be boolean')
    try:from .study_artifact_exports import PublicationLease
    except ImportError:from study_artifact_exports import PublicationLease
    with PublicationLease(root/'.openislet-documentation.lock') as lease:
        return _write_documentation_locked(root,manifest_pin_count,render,lease,rollback_journal)


def _write_documentation_locked(root,manifest_pin_count,render,lease=None,rollback_journal=None):
    facts=derive_facts(root,manifest_pin_count)
    target=root/'data/extracted/documentation_facts_v0.8.0.json'
    prepared={target:json.dumps(facts,indent=2,ensure_ascii=False,allow_nan=False)+'\n'}
    if render:
        doc_candidates = [root/'README.md', root/'CLAUDE.md', root/'ISLET-1-BLUEPRINT.md', *sorted((root/'docs').glob('*.md'))]
        for path in [p for p in doc_candidates if p.exists()]:
            if path.name in ('README.md', 'PROJECT-DECISIONS-LOG.md', 'OPENISLET-2026-2030-MASTER-TIMELINE.md'):
                text=_read_regular(path).decode('utf-8')
                text=re.sub(r'[\d,]+(?= atomic effects)',f"{facts['corpus']['atomic_effects']:,}",text)
                text=re.sub(r'`[\d,]+`(?= pinned files)', '`'+f"{manifest_pin_count:,}"+'`', text)
                prepared[path]=text
                continue
            text=_read_regular(path).decode('utf-8')
            prepared[path]=render_text(text,facts,path.name)
    for relative,digest in facts.get('source_sha256',{}).items():
        if hashlib.sha256(_read_regular(root/relative)).hexdigest()!=digest:
            raise ValueError('A source changed during documentation rendering')
    if lease is not None:lease.require_owned()
    _publish_document_batch(prepared,rollback_journal=rollback_journal)
    if lease is not None:lease.require_owned()
    return facts


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);parser.add_argument('--manifest-pin-count',type=int,required=True);parser.add_argument('--facts-only',action='store_true');args=parser.parse_args()
    write_documentation(args.root,args.manifest_pin_count,not args.facts_only)
