"""Independent documentation arithmetic and authored-structure preservation controls."""
import copy
import importlib.util
import json
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
MODULE=ROOT/'scripts/render_documentation.py'
spec=importlib.util.spec_from_file_location('round2_renderer',MODULE);render=importlib.util.module_from_spec(spec);spec.loader.exec_module(render)
def facts():return copy.deepcopy(json.loads((ROOT/'data/extracted/documentation_facts_v0.8.0.json').read_text()))
def required(ok,message):
    if not ok:raise AssertionError(message)

def test_named_score_columns_and_floor_follow_their_own_estimands():
    f=facts();a=f['layer1']['ranked_architectures'][0];a['composite_series_system_utility']=.4
    a['utility_pediatric_early_onset_t1d']['series_system_joint_survival_probability']=.2
    a['utility_pediatric_early_onset_t1d']['heuristic_conservative_floor_index']=.0123
    a['utility_pediatric_early_onset_t1d']['manufacturing_locus_penalty_factor']=.75
    a['utility_adult_onset_lada']['series_system_joint_survival_probability']=.8
    text='# Test\n| Profile | Grid Position | Architecture ID | Loci | Specification | Pediatric T1D Index | Adult LADA Index | Composite Index | Mfg Factor | Limiting Gate |\n|---|---|---|---|---|---|---|---|---|---|\n'
    text+='| Scenario A | #999 | '+a['architecture_id']+' | K4 | Authored dose 250 IU | **`0.9999`**<br>*(Floor: `0.9999`)* | **`0.8888`** | **`0.7777`** | `0.6666` | Other |\n'
    row=next(x for x in render.render_text(text,f,'ISLET-1-BLUEPRINT.md').splitlines() if a['architecture_id'] in x and x.startswith('|'))
    cells=row.split('|')
    required('0.2000' in cells[6] and '0.0123' in cells[6],'pediatric index or its floor assigned to wrong column')
    required('0.8000' in cells[7],'adult index assigned to wrong column')
    required('0.4000' in cells[8],'composite index assigned to wrong column')
    required('0.7500' in cells[9],'engineering discount assigned to wrong column')
    required('Authored dose 250 IU' in cells[5],'authored specification erased')

@pytest.mark.parametrize('row',[
    '| **`Complex 1`** | Target | 3CDG | 123 contacts; n=17 | 0.1234 A | Authored linker 27 residues |',
    '| **`Complex 2`** | Target | 2JJS | 456 contacts; n=19 | 0.5678 A | Authored linker 31 residues |',
    '| **Stage 0** | Assay | Batch n=17; GSIS >=3; cutoff 25,000 molecules |',
    '| **`GRAVEYARD_004`** | **WARN** | RNLS | poor dynamic GSIS; paired WT n=17; OCR dose=250 |',
    '| **Raw Corpus Manifest** | Authored count n=17 | Authored split 30/70 |',
    '| **Prior-Constant & Auxiliary** | Authored variance17 | Authored per-batch n=19 |',
])
def test_renderer_retains_authored_table_parameters_and_mermaid(row):
    diagram='```mermaid\nflowchart LR\n A[38-45 days] --> B[16-22 weeks]\n```'
    output=render.render_text('# Test\n'+row+'\n\n'+diagram+'\n',facts(),'ISLET-1-BLUEPRINT.md')
    for marker in ['n=17','n=19','0.1234 A','0.5678 A','27 residues','31 residues','GSIS >=3','25,000 molecules','dose=250','count n=17','split 30/70','variance17','per-batch n=19']:
        if marker in row:required(marker in output,'authored engineering parameter erased: '+marker)
    required(diagram in output,'authored Mermaid diagram changed')
