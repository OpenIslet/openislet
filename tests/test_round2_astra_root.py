"""Independent release-boundary and renderer mutation controls; no source writes."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import os
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import generate_release_manifest as export
from scripts import render_documentation as render

def require(ok,why):
    if not ok:raise AssertionError(why)

def facts():return copy.deepcopy(json.loads((ROOT/'data/extracted/documentation_facts_v0.8.0.json').read_text()))

@pytest.mark.parametrize('rel',['../secret','a/../secret','/absolute','a//b','a/./b','a\\b','a\x00b','a\nb','ａｂ.json'])
def test_release_paths_have_one_canonical_identity(rel):
    require(not export._canonical_relative_path(rel),'ambiguous release path accepted')

@pytest.mark.parametrize('kind',['symlink','hardlink','fifo','parent_symlink','directory'])
def test_release_reader_rejects_non_regular_single_identity(tmp_path,kind):
    (tmp_path/'original').write_text('sensitive fixture')
    rel='payload'
    if kind=='symlink':(tmp_path/rel).symlink_to(tmp_path/'original')
    elif kind=='hardlink':os.link(tmp_path/'original',tmp_path/rel)
    elif kind=='fifo':os.mkfifo(tmp_path/rel)
    elif kind=='parent_symlink':
        (tmp_path/'real').mkdir();(tmp_path/'real/file').write_text('fixture');(tmp_path/'alias').symlink_to(tmp_path/'real');rel='alias/file'
    else:(tmp_path/rel).mkdir()
    with pytest.raises(ValueError):export._regular_file_bytes(tmp_path,rel)


def test_pins_check_content_and_size_and_regular_positive(tmp_path):
    p=tmp_path/'payload';p.write_bytes(b'original')
    pin={'sha256':hashlib.sha256(b'original').hexdigest(),'bytes':8};mf={'artifacts':{'payload':pin}}
    export._verify_pinned_source(tmp_path,'payload',mf)
    for field,value in [('bytes',True),('bytes',9),('sha256','0'*64)]:
        bad=copy.deepcopy(mf);bad['artifacts']['payload'][field]=value
        with pytest.raises(ValueError):export._verify_pinned_source(tmp_path,'payload',bad)


def minimal_export(root,monkeypatch):
    root.mkdir();(root/'LICENSE').write_bytes(b'license');(root/'data/dossiers').mkdir(parents=True)
    (root/'data/dossiers/12345678.json').write_bytes(b'{"abstract":"original"}')
    pins={rel:{'sha256':hashlib.sha256((root/rel).read_bytes()).hexdigest(),'bytes':(root/rel).stat().st_size} for rel in ['LICENSE','data/dossiers/12345678.json']}
    mf={'artifacts':pins,'public_redistribution_allowlist_spec':{'public_package_included_paths':['LICENSE','manifest.json'],'local_research_cache_excluded_from_public_bundle':['data/dossiers/*.json']}}
    (root/'manifest.json').write_text(json.dumps(mf));monkeypatch.setattr(export,'ROOT',root);monkeypatch.setattr(export,'DOSSIERS_DIR',root/'data/dossiers')

@pytest.mark.parametrize('mutation',['source','missing','bytes'])
def test_export_preflight_fails_before_publishing_changed_private_source_or_missing_pins(tmp_path,monkeypatch,mutation):
    root=tmp_path/'core';minimal_export(root,monkeypatch)
    if mutation=='source':(root/'data/dossiers/12345678.json').write_text('{"abstract":"tampered"}')
    elif mutation=='missing':(root/'LICENSE').unlink()
    else:
        mf=json.loads((root/'manifest.json').read_text());mf['artifacts']['LICENSE']['bytes']=9000;(root/'manifest.json').write_text(json.dumps(mf))
    dest=tmp_path/'public'
    with pytest.raises((ValueError,FileNotFoundError)):export.export_public_bundle(dest)
    require(not dest.exists(),'invalid release was published')
    require(not list(tmp_path.glob('*.publish-lock')),'publication lock leaked')


def test_export_refuses_repository_descendant_without_creating_it(tmp_path,monkeypatch):
    root=tmp_path/'core';root.mkdir();monkeypatch.setattr(export,'ROOT',root)
    with pytest.raises(ValueError):export.export_public_bundle(root/'new/child')
    require(not (root/'new').exists(),'invalid output validation created source directories')


def test_publication_lock_never_steals_current_owner(tmp_path):
    path=tmp_path/'lock'
    with export._publication_lock(path):
        with pytest.raises(ValueError):
            with export._publication_lock(path):pass
        require(path.is_dir(),'contender removed owner lock')
    require(not path.exists(),'lock remained after successful owner exit')

@pytest.mark.parametrize('text',['α β γ δ ε ζ η θ','Ａ Ｂ Ｃ Ｄ Ｅ Ｆ','source\u200b words remain the same here'])
def test_unicode_source_overlap_counts_real_characters(text):
    source=text.replace('\u200b','')
    require(export._max_verbatim_shingle_chars({'note':text},source)==len(text),'Unicode text escaped aggregate quote accounting')

@pytest.mark.parametrize('text',[
 '# Test\n'+render.START+'\nAUTHORED', '# Test\n'+render.END,
 '# Test\n'+render.END+'\n'+render.START,
 '# Test\n'+render.START+render.END+render.START+render.END,
 '# Test\n<!-- OPENISLET_CURRENT_LAYER_SCOPE -->\nAUTHORED',
])
def test_malformed_generated_blocks_fail_closed(text):
    with pytest.raises(ValueError):render.render_text(text,facts(),'ISLET-1-BLUEPRINT.md')


def test_literal_mermaid_and_code_are_never_subject_to_scientific_replacements():
    diagram='```mermaid\nflowchart LR\n A["GSG-tPT2A"] --> B["MVDGTLLLL"]\n```'
    out=render.render_text('# Test\n'+diagram+'\n',facts(),'ISLET-1-BLUEPRINT.md')
    require(diagram in out,'authored Mermaid changed')
    require(render.render_text(out,facts(),'ISLET-1-BLUEPRINT.md')==out,'renderer not idempotent')

@pytest.mark.parametrize('header',['| ID | Pediatric | Pediatric | Adult | Composite |','| Architecture ID | Pediatric | Unlabelled | Composite |','| ID | Pediatric Adult | Composite |'])
def test_ambiguous_ranking_headers_rejected(header):
    text='# Test\n'+header+'\n|'+'---|'*(len(header.split('|'))-2)+'\n'
    with pytest.raises(ValueError):render.render_text(text,facts(),'ISLET-1-BLUEPRINT.md')

@pytest.mark.parametrize('value',[14999,15000,15001,None])
def test_resting_density_prose_recomputes_threshold_relation(value):
    f=facts();p=next(x for x in f['layer2']['profiles'] if x['genotype_id']=='SCENARIO_A_ARM_1B_UNIVERSAL_LEAN_4LOCUS');p['hour0_resting_unstimulated_hlae_molecules_per_cell']=value
    expected='unavailable relative to' if value is None else 'below' if value<15000 else 'at or above'
    require('the resting endpoint is '+expected+' the assumed 15,000' in render.narratives(f)['layer2'],'stale threshold verdict')

@pytest.mark.parametrize('value',[True,False,float('nan'),float('inf'),float('-inf')])
def test_nonfinite_and_boolean_doc_numbers_are_not_printed(value):
    with pytest.raises(ValueError):render.f(value)

@pytest.mark.parametrize('indent,opening,closing',[('', '```', '````'),('   ','```','```'),('  ','~~~','~~~~'),('', '```','')])
def test_commonmark_literal_fences_keep_authored_bytes(indent,opening,closing):
    diagram=indent+opening+'mermaid\n'+indent+'A["GSG-tPT2A"] --> B["MVDGTLLLL"]\n'+(indent+closing if closing else '')
    out=render.render_text('# Test\n'+diagram,facts(),'ISLET-1-BLUEPRINT.md')
    require(diagram.rstrip('\n') in out,'valid CommonMark literal was rewritten')

@pytest.mark.parametrize('mutation',['none','bytes','boolean_bytes','hash','symlink','hardlink','fifo','missing','traversal'])
def test_acceptance_manifest_boundary_rejects_mutations(tmp_path,mutation):
    p=tmp_path/'file';p.write_bytes(b'fixture')
    mf={'artifacts':{'file':{'sha256':hashlib.sha256(b'fixture').hexdigest(),'bytes':7}}}
    if mutation=='bytes':mf['artifacts']['file']['bytes']=8
    elif mutation=='boolean_bytes':mf['artifacts']['file']['bytes']=True
    elif mutation=='hash':mf['artifacts']['file']['sha256']='0'*64
    elif mutation=='symlink':p.rename(tmp_path/'real');p.symlink_to(tmp_path/'real')
    elif mutation=='hardlink':os.link(p,tmp_path/'alias')
    elif mutation=='fifo':p.unlink();os.mkfifo(p)
    elif mutation=='missing':p.unlink()
    elif mutation=='traversal':mf['artifacts']['../file']=mf['artifacts'].pop('file')
    require(bool(export.manifest_integrity_errors(tmp_path,mf))==(mutation!='none'),'acceptance manifest did not reject exact mutation')

def test_public_manifest_inventory_cannot_omit_or_add_payload(tmp_path):
    (tmp_path/'file').write_bytes(b'x');(tmp_path/'.openislet_public_bundle').write_text('fixture')
    mf={'artifacts':{'file':{'sha256':hashlib.sha256(b'x').hexdigest(),'bytes':1}},'public_export_file_inventory':['file']}
    require(not export.manifest_integrity_errors(tmp_path,mf,True),'clean fixture failed')
    (tmp_path/'extra').write_text('unlisted')
    require('unexpected:extra' in export.manifest_integrity_errors(tmp_path,mf,True),'unlisted payload passed')
    mf['public_export_file_inventory']=['file','file']
    require(export.manifest_integrity_errors(tmp_path,mf,True)==['invalid-public-inventory'],'duplicate inventory passed')


def test_document_batch_validates_before_any_authored_file_is_replaced(tmp_path,monkeypatch):
    (tmp_path/'docs').mkdir()
    for name in ['README.md','CLAUDE.md','ISLET-1-BLUEPRINT.md']:(tmp_path/name).write_text('# Original\n')
    marker='<!-- OPENISLET_CURRENT_LAYER_SCOPE -->x<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->'
    (tmp_path/'docs/bad.md').write_text('# Bad\n'+marker+marker)
    before={p:p.read_bytes() for p in tmp_path.rglob('*.md')}
    monkeypatch.setattr(render,'derive_facts',lambda *args: facts())
    with pytest.raises(ValueError):render.write_documentation(tmp_path,0)
    require(all(p.read_bytes()==b for p,b in before.items()),'failed document batch partially overwrote prose')


def test_unmarked_generated_cell_updates_numbers_without_adding_scope_prose():
    old=' Current corpus: 593 arms, 21 reference controls. '
    current='Current corpus: 592 arms, 20 reference controls.'
    require(render._refine_table_cell(old,'STAGE0',current)==' '+current+' ','numeric update added prose or changed table spacing')


def test_public_cache_exception_does_not_exempt_symlinks(tmp_path):
    (tmp_path/'file').write_bytes(b'x');(tmp_path/'.openislet_public_bundle').write_text('fixture')
    mf={'artifacts':{'file':{'sha256':hashlib.sha256(b'x').hexdigest(),'bytes':1}},'public_export_file_inventory':['file']}
    (tmp_path/'__pycache__').symlink_to(tmp_path.parent,target_is_directory=True)
    require('unsafe:__pycache__' in export.manifest_integrity_errors(tmp_path,mf,True),'cache symlink bypassed CHK07')


def test_public_cache_exception_does_not_exempt_hardlinks(tmp_path):
    (tmp_path/'file').write_bytes(b'x');(tmp_path/'.openislet_public_bundle').write_text('fixture')
    mf={'artifacts':{'file':{'sha256':hashlib.sha256(b'x').hexdigest(),'bytes':1}},'public_export_file_inventory':['file']}
    (tmp_path/'__pycache__').mkdir();os.link(tmp_path/'file',tmp_path/'__pycache__/cached.pyc')
    require('unsafe:__pycache__/cached.pyc' in export.manifest_integrity_errors(tmp_path,mf,True),'cache hardlink bypassed CHK07')


def test_source_normalization_cache_is_content_keyed_and_bounded():
    export._normalized_source.cache_clear()
    phrase='alpha beta gamma delta epsilon zeta'
    require(export._max_verbatim_shingle_chars(phrase,phrase)==len(phrase),'positive source overlap lost')
    require(export._max_verbatim_shingle_chars(phrase,'one two three four five six')==0,'changed content reused an old source overlap')
    for i in range(100):export._normalized_source(str(i)+' independent source')
    require(export._normalized_source.cache_info().currsize<=16,'source cache grew without bound')
