"""Two-pass master-audit molecular negative controls; all fixtures stay temporary."""
import ast
import copy
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from scripts import run_layer2_perturbation_and_regulatory_plane as l2
from scripts import run_layer3_structural_binder_plane as l3

def require(value, message="master molecular validation failed"):
    if not value: raise AssertionError(message)

def atlas():
    return json.loads((ROOT / "data/extracted/layer2_regulatory_and_perturbation_atlas.json").read_text())

@pytest.mark.parametrize("file", ["run_layer2_perturbation_and_regulatory_plane.py", "run_layer3_structural_binder_plane.py"])
def test_production_validation_survives_optimization(file):
    tree = ast.parse((ROOT / "scripts" / file).read_text())
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), file)

@pytest.mark.parametrize("dna", ["ATGA", "ATN", "", "atg", None])
def test_invalid_coding_dna_rejected(dna):
    with pytest.raises((ValueError, TypeError)): l2.translate_dna(dna)

@pytest.mark.parametrize("opt", [False, True])
def test_source_invariants_are_rejected_in_real_python_modes(opt):
    code = """from scripts.run_layer2_perturbation_and_regulatory_plane import translate_dna
from scripts.run_layer3_structural_binder_plane import validated_pae_matrix, euclidean_distance
for fn in (lambda:translate_dna('ATGA'),lambda:validated_pae_matrix([{'predicted_aligned_error':[[0.0]]}],376,'truncated'),lambda:euclidean_distance(None,[0,0,0])):
    try: fn()
    except (ValueError, TypeError): pass
    else: raise RuntimeError('invalid input survived')
print('three invalid inputs rejected')
"""
    r = subprocess.run([sys.executable] + (["-O"] if opt else []) + ["-c", code], cwd=ROOT, text=True, capture_output=True)
    require(r.returncode == 0 and "three invalid inputs rejected" in r.stdout, r.stderr)

@pytest.mark.parametrize("mutation", ["shape", "nan", "negative", "infinite"])
def test_pae_requires_complete_finite_matrix(mutation):
    matrix = [[0.0, 1.0], [2.0, 3.0]]
    if mutation == "shape": matrix.pop()
    else: matrix[0][0] = {"nan":float("nan"),"negative":-1.0,"infinite":float("inf")}[mutation]
    with pytest.raises(ValueError): l3.validated_pae_matrix([{"predicted_aligned_error":matrix}], 2, "fixture")
    require(l3.validated_pae_matrix([{"predicted_aligned_error":[[0.0,1.0],[2.0,3.0]]}],2,"valid").shape == (2,2))

def test_missing_reference_anchor_cannot_become_zero_tension(monkeypatch):
    original = l3.parse_mmcif_atoms
    def missing_anchor(path,*args,**kwargs):
        return [a for a in original(path,*args,**kwargs) if not (Path(path).name == "pdb_00003cdg.cif" and a["chain"] == "B" and a["seq_id"] == 99 and a["atom_id"] == "C")]
    monkeypatch.setattr(l3,"parse_mmcif_atoms",missing_anchor)
    with pytest.raises(ValueError, match="anchor is missing"):
        l3.analyze_complex_1_hlae_nkg2a(atlas()["construct_sequence_and_grammar_verification"]["processed_polypeptide_ledger"])

@pytest.mark.parametrize("mutation", ["hash", "length", "sequence", "missing"])
def test_stale_layer2_polypeptide_ledger_rejected(mutation):
    ledger = atlas()["construct_sequence_and_grammar_verification"]["processed_polypeptide_ledger"]
    l3.validate_processed_polypeptide_ledger(ledger)
    key = "PDL1_mature_cell_surface"
    if mutation == "missing": del ledger[key]
    elif mutation == "hash": ledger[key]["sha256"] = "0"*64
    elif mutation == "length": ledger[key]["mature_length_aa"] += 1
    else: ledger[key]["sequence"] = "X" + ledger[key]["sequence"][1:]
    with pytest.raises(ValueError): l3.validate_processed_polypeptide_ledger(ledger)

@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, 0, -1])
def test_invalid_physical_cutoffs_rejected(value):
    with pytest.raises(ValueError): l3.compute_interface_metrics([], [], heavy_cutoff=value)
    with pytest.raises(ValueError): l3.worm_like_chain_metrics(20, 10.0, lp=value)

@pytest.mark.parametrize("points", [[], [[0.,0.,0.]], [[float("nan"),0.,0.]]*3])
def test_kabsch_rejects_incomplete_or_nonfinite_points(points):
    with pytest.raises(ValueError): l3.kabsch_rmsd(np.asarray(points),np.asarray(points))

def test_kabsch_proper_rotation_does_not_fit_a_reflection():
    p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    q=p.copy();q[:,0]*=-1
    rmsd,_,_=l3.kabsch_rmsd(p,q)
    require(rmsd > .1)
    require(l3.kabsch_rmsd(p,p+7)[0] < 1e-12)

@pytest.mark.parametrize("value", [float("nan"),float("inf"),-1.0,0.0,True])
def test_pwm_rejects_invalid_pseudocount(value):
    with pytest.raises(ValueError): l2.scan_jaspar_pwm("ACGT", {}, pseudocount=value)

def test_mmcif_nonfinite_coordinate_fails_closed(tmp_path):
    src=(ROOT/"data/structures/pdb_00003cdg.cif").read_text().splitlines()
    cols=[];changed=False
    for i,line in enumerate(src):
        if line.startswith("_atom_site."): cols.append(line.split(".",1)[1].strip())
        if line.startswith("ATOM ") and not changed:
            parts=line.split();parts[cols.index("Cartn_x")]="NaN";src[i]=" ".join(parts);changed=True
    require(changed)
    p=tmp_path/"bad.cif";p.write_text("\n".join(src)+"\n")
    with pytest.raises(ValueError,match="Non-finite"): l3.parse_mmcif_atoms(p)

def test_l2_excludes_summaries_reviews_and_superseded_evidence(tmp_path):
    a=atlas();ref=json.loads((ROOT/"data/raw/layer2_genomic_and_jaspar_reference.json").read_text())
    data=[{"pmid":"1","effect_id":"active","source_text_type":"Abstract_Only","focal_factor_X":"CD47","effect_direction":"Protective_or_Improved"}]
    p=tmp_path/"effects.jsonl"
    def replay(rows):
        p.write_text("\n".join(json.dumps(x) for x in rows)+"\n")
        return l2.simulate_stage0_to_stage7_perturbation_plane(ref,copy.deepcopy(a["construct_sequence_and_grammar_verification"]),a["jaspar_pwm_regulatory_motif_scans"],p)
    baseline=replay(data)
    for i,source in enumerate(["Project_Summary","Preprint_Superseded_By_Peer_Reviewed","Review_or_Commentary"]):
        data.append({"pmid":str(i+2),"effect_id":source,"source_text_type":source,"focal_factor_X":"CD47","effect_direction":"Protective_or_Improved"})
    filtered=replay(data)
    require(baseline["perturbation_profiles"]==filtered["perturbation_profiles"])
    require(filtered["evidence_eligibility"]["excluded_effect_records"]==3)
    data.append({"pmid":"9","effect_id":"unknown","focal_factor_X":"CD47","effect_direction":"Protective_or_Improved"})
    with pytest.raises(ValueError,match="source_text_type"): replay(data)
