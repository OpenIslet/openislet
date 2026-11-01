"""Round-2 molecular regression oracles; all assay/source mutations are synthetic and temporary."""
from __future__ import annotations
import copy
import json
from decimal import Decimal, localcontext
from pathlib import Path
import sys
import math
import pytest
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import run_layer2_perturbation_and_regulatory_plane as l2
from scripts import run_layer3_structural_binder_plane as l3


def require(ok, message):
    if not ok:raise AssertionError(message)

@pytest.fixture(scope='module')
def source():
    a=json.loads((ROOT/'data/extracted/layer2_regulatory_and_perturbation_atlas.json').read_text())
    return {'ref':json.loads((ROOT/'data/raw/layer2_genomic_and_jaspar_reference.json').read_text()),
            'constructs':a['construct_sequence_and_grammar_verification'],'motifs':a['jaspar_pwm_regulatory_motif_scans']}

def effect(**patch):
    return {'pmid':'1','effect_id':'pmid1_eff1','source_text_type':'Abstract_Only','focal_factor_X':'CD47',
            'effect_direction':'Protective_or_Improved','secretome_or_differentiation_penalty_observed':False,**patch}

def replay(source,tmp_path,rows=None,ref=None,motifs=None):
    p=tmp_path/'effects.jsonl';p.write_text(''.join(json.dumps(e)+'\n' for e in (rows if rows is not None else [effect()])))
    return l2.simulate_stage0_to_stage7_perturbation_plane(
        copy.deepcopy(ref if ref is not None else source['ref']),copy.deepcopy(source['constructs']),
        copy.deepcopy(motifs if motifs is not None else source['motifs']),p)

@pytest.mark.parametrize('value',[True,False,'0',None,10**400])
def test_pae_rejects_values_that_only_coerce_to_confidence(value):
    with pytest.raises(ValueError):l3.validated_pae_matrix({'predicted_aligned_error':[[0.,value],[1.,0.]]},2,'synthetic')

@pytest.mark.parametrize('length',[True,False,0,-1,2.0])
def test_pae_requires_exact_positive_residue_count(length):
    with pytest.raises(ValueError):l3.validated_pae_matrix({'predicted_aligned_error':[[0.,1.],[1.,0.]]},length,'synthetic')

def test_pae_accepts_actual_numeric_values_without_coercion():
    p=np.array([[0,2],[4,0]],dtype=np.int64)
    require(np.array_equal(l3.validated_pae_matrix({'predicted_aligned_error':p},2,'numeric'),p),'numeric PAE changed')

def test_processed_entry_must_be_object(source):
    p=copy.deepcopy(source['constructs']['processed_polypeptide_ledger']);p['PDL1_mature_cell_surface']=[]
    with pytest.raises(TypeError):l3.validate_processed_polypeptide_ledger(p)

@pytest.mark.parametrize('lp',[3.5,5.0,8.0])
def test_wlc_selected_parameter_does_not_relabel_fixed_reference(lp):
    default=l3.worm_like_chain_metrics(21,30.)
    actual=l3.worm_like_chain_metrics(21,30.,lp=lp)
    for key in ['marko_siggia_entropic_tension_pN_at_310K_lp4A',
                'marko_siggia_entropic_tension_pN_at_l0_3p8A_lp4A',
                'hemispherical_arc_wlc_tension_pN_at_310K_lp4A',
                'entropic_tension_by_persistence_length_pN']:
        require(actual[key]==default[key],key+' silently changed with a different selected lp')
    d=30./73.5;force=(Decimal('0.25')/(Decimal(1)-Decimal(str(d)))**2-Decimal('.25')+Decimal(str(d)))*Decimal(str(l3.KBT_310K_PN_ANGSTROM))/Decimal(str(lp))
    require(actual['marko_siggia_entropic_tension_pN_at_selected_lp']==round(float(force),3),'selected force disagrees with independent formula')
    require(actual['selected_persistence_length_angstrom']==lp,'selected parameter provenance missing')

@pytest.mark.parametrize('lp',[1e2,1e4,1e8,1e12])
def test_wlc_rod_limit_matches_high_precision_second_moment(lp):
    with localcontext() as ctx:
        ctx.prec=80;contour=Decimal('73.5');p=Decimal(str(lp))
        expected=(2*p*contour*(1-(p/contour)*(1-(-contour/p).exp()))).sqrt()
    result=l3.worm_like_chain_metrics(21,30.,lp=lp)
    require(result['wlc_rms_end_to_end_angstrom']==round(float(expected),3),'WLC RMS lost accuracy')
    require(result['wlc_rms_end_to_end_angstrom']<=result['contour_length_angstrom'],'RMS exceeds contour')

@pytest.mark.parametrize('kwargs',[{'num_residues':10**400,'span_angstrom':0.},
                                 {'num_residues':21,'span_angstrom':1.,'l0':1e308},
                                 {'num_residues':21,'span_angstrom':30.,'lp':1e-308},
                                 {'num_residues':21,'span_angstrom':10**400},
                                 {'num_residues':21,'span_angstrom':1.,'lp':10**400}])
def test_wlc_rejects_nonfinite_derived_results(kwargs):
    with pytest.raises(ValueError):l3.worm_like_chain_metrics(**kwargs)

def test_wlc_zero_force_remains_finite_at_extreme_flexibility():
    result=l3.worm_like_chain_metrics(21,0.,lp=1e-308)
    require(result['marko_siggia_entropic_tension_pN_at_selected_lp']==0.,'zero-span force should vanish')
    json.dumps(result,allow_nan=False)

@pytest.mark.parametrize('value',[True,False,'38.4',float('nan'),float('inf'),-.1,1_000_001.,10**400])
def test_consumed_tpm_cannot_be_coerced_clamped_or_overflowed(source,tmp_path,value):
    r=copy.deepcopy(source['ref']);r['gtex_v8_v10_human_tissue_expression']['CD47']['pancreas_median_tpm']=value
    with pytest.raises(ValueError,match='GTEx CD47'):replay(source,tmp_path,ref=r)

def test_zero_cd47_tpm_retains_undefined_ratio_instead_of_crashing(source,tmp_path):
    r=copy.deepcopy(source['ref']);r['gtex_v8_v10_human_tissue_expression']['CD47']['pancreas_median_tpm']=0.
    out=replay(source,tmp_path,ref=r);note=out['hour0_hlae_inducible_vs_constitutive_kinetics_note']
    require(note['arm1b_cd47_to_wt_fold'] is None,'undefined ratio fabricated')
    require(note['arm1b_cd47_to_wt_fold_status']=='UNDEFINED_ZERO_REFERENCE_CD47_DENSITY','missing ratio status')
    require(all(p['hour0_basal_cd47_molecules_per_cell']==0 for p in out['perturbation_profiles']),'zero expression was replaced')
    json.dumps(out,allow_nan=False)

@pytest.mark.parametrize('region,matrix', [('CXCL10_promoter_chr4_76023197_76023797','MA0137.3_STAT1'),
                                        ('B2M_promoter_chr15_44711217_44711817','MA0050.2_IRF1')])
@pytest.mark.parametrize('value',[True,'7.5',float('nan'),float('inf')])
def test_consumed_pwm_scores_must_be_real_finite_values(source,tmp_path,region,matrix,value):
    motifs=copy.deepcopy(source['motifs']);motifs[region]['top_jaspar_pwm_hits'][matrix]['top_hits'][0]['log2_odds_score']=value
    with pytest.raises(ValueError,match='PWM score'):replay(source,tmp_path,motifs=motifs)

def test_missing_consumed_reference_does_not_silently_reuse_default(source,tmp_path):
    r=copy.deepcopy(source['ref']);del r['gtex_v8_v10_human_tissue_expression']['CD47']['pancreas_median_tpm']
    with pytest.raises(ValueError,match='Missing or malformed'):replay(source,tmp_path,ref=r)

@pytest.mark.parametrize('token,canonical',[('PD-L1','CD274'),('PDL1','CD274'),('CD274','CD274'),('HLA-E','HLAE'),('HLAE','HLAE'),('ICASP9','CASP9'),('CASP9','CASP9')])
def test_aliases_count_one_canonical_directional_link(source,tmp_path,token,canonical):
    out=replay(source,tmp_path,[effect(focal_factor_X=token)]);c=out['calibration_inputs_consumed']
    require(c['ring1_directional_study_gene_links_count']==1,'one alias became multiple evidence links')
    require(c['ring1_canonical_gene_effect_evidence_counts'][canonical]['protective']==1,'canonical target lost')
    require(c['legacy_display_counts_must_not_be_summed'] is True,'legacy duplication not qualified')
    require(c['directional_link_normalization_denominator_prior']==148.,'normalization prior retuned')

def test_same_study_alias_rows_do_not_inflate_prior_and_exact_replay_deduplicates(source,tmp_path):
    one=effect(focal_factor_X='PD-L1');alias=effect(effect_id='pmid1_eff2',focal_factor_X='CD274')
    a=replay(source,tmp_path,[one]);b=replay(source,tmp_path,[one,alias,one])
    require(a['perturbation_profiles']==b['perturbation_profiles'],'equivalent aliases changed prior profiles')
    require(b['calibration_inputs_consumed']['ring1_directional_study_gene_links_count']==1,'alias rows inflated evidence')
    require(b['evidence_eligibility']['duplicate_exact_effect_records_skipped']==1,'exact repeat not recorded')

@pytest.mark.parametrize('patch',[{'pmid':'2'},{'effect_direction':'Deleterious_or_Impaired'},{'focal_factor_X':'TYK2'}])
def test_conflicting_duplicate_effect_identity_fails_closed(source,tmp_path,patch):
    with pytest.raises(ValueError):replay(source,tmp_path,[effect(),effect(**patch)])

@pytest.mark.parametrize('patch',[{'pmid':True},{'pmid':''},{'effect_id':''},{'superseded_by_pmid':True},
                                  {'secretome_or_differentiation_penalty_observed':'false'},
                                  {'effect_direction':'Not_Protective'},{'focal_factor_X':{'gene':'CD47'}},
                                  {'arbitrary_nonfinite_field':float('nan')}])
def test_atomic_effect_identity_direction_and_flag_types_are_explicit(source,tmp_path,patch):
    with pytest.raises(ValueError):replay(source,tmp_path,[effect(**patch)])


def test_wlc_infeasible_comparator_is_null_without_crashing_valid_selected_geometry():
    out=l3.worm_like_chain_metrics(21,80.,l0=4.)
    require(out['geometrically_feasible_d_lt_Lc'] is True,'selected geometry incorrectly rejected')
    require(out['marko_siggia_entropic_tension_pN_at_l0_3p8A_lp4A'] is None,'infeasible comparator fabricated force')
    require(out['comparator_l0_3p8A_geometry_status']=='GEOMETRICALLY_INFEASIBLE_SPAN_GE_COMPARATOR_CONTOUR','missing comparator infeasibility status')
    require('unavailable' in out['physical_limitations_note'],'infeasible comparator absent from explanation')
    json.dumps(out,allow_nan=False)


@pytest.mark.parametrize('kwargs',[
    {'num_residues':100000000,'span_angstrom':1.5e308,'l0':1.6e300},
    {'num_residues':4*10**307,'span_angstrom':1.3e308,'l0':4.},
])
def test_wlc_nonfinite_surface_bypass_geometry_fails_explicitly(kwargs):
    with pytest.raises(ValueError,match='bypass arc'):
        l3.worm_like_chain_metrics(**kwargs)


def test_wlc_infeasible_branch_cannot_emit_nonfinite_fractional_extension():
    with pytest.raises(ValueError,match='fractional extension'):
        l3.worm_like_chain_metrics(1,1e308,l0=5e-324)
