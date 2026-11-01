"""Independent numerical, publication, literal and producer-boundary controls."""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from scripts import engineering_priors as engineering
from scripts import extract_ring1_corpus as extraction
from scripts import generate_release_manifest as export
from scripts import public_excerpt_policy as excerpts
from scripts import render_documentation as renderer
from scripts import study_artifact_exports as materialize

def require(condition,message):
    if not condition:raise AssertionError(message)


def facts():
    return json.loads((ROOT/'data/extracted/documentation_facts_v0.8.0.json').read_text())


@pytest.mark.parametrize('treatment,control',[
    (1e308,math.nextafter(1e308,0.)),(math.nextafter(1e308,0.),1e308),
    (1.,math.nextafter(1.,0.)),(1e-300,math.nextafter(1e-300,0.)),
    (math.ulp(0.),2*math.ulp(0.)),(1e308,math.ulp(0.)),
    (math.ulp(0.),1e308),(.5,1.),(1.,1.),
])
def test_log_response_ratio_against_independent_decimal_oracle(treatment,control):
    with localcontext() as context:
        context.prec=120
        expected=float((Decimal.from_float(treatment)/Decimal.from_float(control)).ln())
    actual=extraction._log_response_ratio(treatment,control)
    require(math.isclose(actual,expected,rel_tol=2e-13,abs_tol=0.),'binary64 log ratio lost a real contrast or overflowed')


def test_engineering_utility_arbitrary_integer_counts_have_mathematical_limits():
    require(engineering.engineering_complexity_discount(10**400,10**400,0.,1.)==1.,'zero penalties must remain exactly neutral')
    require(engineering.engineering_complexity_discount(4,10**400,0.,.965)==0.,'large base count must have its representable zero limit')
    with localcontext() as context:
        context.prec=100
        expected=float((-(Decimal.from_float(math.ulp(0.))*Decimal(10**320-4))).exp())
    actual=engineering.engineering_complexity_discount(10**320,0,math.ulp(0.),1.)
    require(math.isclose(actual,expected,rel_tol=2e-15),'small finite count/penalty product was mistaken for infinity')


@pytest.mark.parametrize('parameter',[10**400,True,float('inf'),float('nan'),-1.])
def test_engineering_parameters_fail_with_explicit_domain_error(parameter):
    with pytest.raises(ValueError):engineering.engineering_complexity_discount(4,0,parameter,.965)


def study(pmid='12345678',**effect):
    return {'study_metadata':{'pmid':pmid,'citation_short':'Synthetic test'},'atomic_effects':[{'effect_id':'test-effect',**effect}]}


@pytest.mark.parametrize('mutation',['duplicate','override','nan','wrong-container'])
def test_jsonl_validation_cannot_publish_either_half(tmp_path,mutation):
    tmp_path=tmp_path.resolve();folder=tmp_path/'data/extracted';folder.mkdir(parents=True)
    paths=[folder/'ring1_study_dossiers_v0.5.jsonl',folder/'ring1_atomic_effects_v0.5.jsonl']
    for path in paths:path.write_bytes(b'previous run\n')
    records=[study()]
    if mutation=='duplicate':records.append(study())
    elif mutation=='override':records[0]['atomic_effects'][0]['pmid']='87654321'
    elif mutation=='nan':records[0]['atomic_effects'][0]['value']=float('nan')
    else:records[0]['atomic_effects']={}
    with pytest.raises(ValueError):materialize.write_study_exports(tmp_path,records)
    require(all(path.read_bytes()==b'previous run\n' for path in paths),'validation left mixed-run or truncated aliases')


def test_jsonl_second_publication_failure_restores_first_and_releases_lock(tmp_path,monkeypatch):
    tmp_path=tmp_path.resolve();folder=tmp_path/'data/extracted';folder.mkdir(parents=True)
    paths=[folder/'ring1_study_dossiers_v0.5.jsonl',folder/'ring1_atomic_effects_v0.5.jsonl']
    for path in paths:path.write_bytes(b'previous run\n')
    original=os.replace
    def fail_effect(source,destination):
        if Path(destination)==paths[1]:raise OSError('injected effect publication failure')
        return original(source,destination)
    monkeypatch.setattr(materialize.os,'replace',fail_effect)
    with pytest.raises(OSError):materialize.write_study_exports(tmp_path,[study()])
    require(all(path.read_bytes()==b'previous run\n' for path in paths),'pair rollback failed')
    require(not (folder/'.study-exports.lock').exists(),'producer lock leaked')


def test_jsonl_writer_does_not_steal_an_existing_lock(tmp_path):
    tmp_path=tmp_path.resolve();folder=tmp_path/'data/extracted';folder.mkdir(parents=True)
    lock=folder/'.study-exports.lock';lock.mkdir()
    with pytest.raises(ValueError):materialize.write_study_exports(tmp_path,[study()])
    require(lock.is_dir(),'contending writer removed the owner lock')


def test_jsonl_ordinary_bytes_match_declared_row_order_and_identity(tmp_path):
    records=[study(value=.125)];a,b,n=materialize.write_study_exports(tmp_path.resolve(),records)
    row={'pmid':'12345678','citation_short':'Synthetic test','source_text_type':None,'superseded_by_pmid':None,'effect_id':'test-effect','value':.125}
    require(a.read_text()==json.dumps(records[0],ensure_ascii=False)+'\n','study bytes changed format')
    require(b.read_text()==json.dumps(row,ensure_ascii=False)+'\n' and n==1,'effect bytes or source identity changed')


def test_shared_publication_restores_files_in_different_parent_directories(tmp_path,monkeypatch):
    tmp_path=tmp_path.resolve();a=tmp_path/'omics/prior.json';b=tmp_path/'extracted/summary.json'
    for path in [a,b]:path.parent.mkdir();path.write_bytes(b'previous')
    original=os.replace
    def fail_second(source,destination):
        if Path(destination)==b:raise OSError('injected second-parent failure')
        return original(source,destination)
    monkeypatch.setattr(materialize.os,'replace',fail_second)
    with pytest.raises(OSError):materialize.publish_byte_artifacts({a:b'new prior',b:b'new summary'},lock_path=tmp_path/'.shared.lock')
    require(a.read_bytes()==b.read_bytes()==b'previous','rollback did not span both output parents')
    require(not (tmp_path/'.shared.lock').exists(),'shared byte publication leaked its writer lock')


@pytest.mark.parametrize('mutation',['bad-bytes','overlap','parent-alias','dotdot'])
def test_shared_publication_rejects_invalid_paths_before_writing(tmp_path,mutation):
    tmp_path=tmp_path.resolve();a=tmp_path/'a';a.write_bytes(b'previous')
    if mutation=='bad-bytes':payloads={a:'not bytes'}
    elif mutation=='overlap':payloads={tmp_path/'new':b'x',tmp_path/'new/child':b'y'}
    elif mutation=='parent-alias':
        alias=tmp_path/'alias';alias.symlink_to(tmp_path,target_is_directory=True);payloads={alias/'a':b'changed'}
    else:payloads={tmp_path/'new/../a':b'changed'}
    with pytest.raises(ValueError):materialize.publish_byte_artifacts(payloads)
    require(a.read_bytes()==b'previous' and not (tmp_path/'new').exists(),'rejected publication touched payload state')


@pytest.mark.parametrize('lock_name,target_name',[
    ('lease','lease/payload'),('LEASE','lease/payload'),
    ('caf\u00e9','cafe\u0301/payload'),('payload/lease','payload'),('payload','payload'),
])
def test_shared_publication_rejects_symmetric_portable_lock_overlap_before_mkdir(tmp_path,lock_name,target_name):
    root=tmp_path.resolve();lock=root/lock_name;target=root/target_name
    with pytest.raises(ValueError,match='overlap'):
        materialize.publish_byte_artifacts({target:b'payload'},lock_path=lock)
    require(not list(root.iterdir()),'overlapping lock and destination created payload, lease or parent state')


@pytest.mark.parametrize('existing',[False,True])
@pytest.mark.parametrize('names',[('a.json','A.JSON'),('café.json','cafe\u0301.json')])
def test_shared_publication_rejects_case_and_unicode_aliases_before_replacement(tmp_path,existing,names):
    tmp_path=tmp_path.resolve();a,b=[tmp_path/name for name in names]
    if existing:a.write_bytes(b'previous')
    with pytest.raises(ValueError):materialize.publish_byte_artifacts({a:b'first',b:b'second'})
    require(a.read_bytes()==b'previous' if existing else not a.exists() and not b.exists(),'aliased destination publication touched bytes')


@pytest.mark.parametrize('intervening_bytes',[False,True])
def test_shared_publication_revocation_stops_commits_and_preserves_new_owner(tmp_path,monkeypatch,intervening_bytes):
    tmp_path=tmp_path.resolve();a=tmp_path/'a';b=tmp_path/'b';lock=tmp_path/'.lease'
    a.write_bytes(b'old a');b.write_bytes(b'old b');original=os.replace
    def revoke_after_first(source,destination):
        result=original(source,destination)
        if Path(destination)==a:
            shutil.rmtree(lock);lock.mkdir();(lock/'owner').write_bytes(b'new owner')
            if intervening_bytes:a.write_bytes(b'new owner bytes')
        return result
    monkeypatch.setattr(materialize.os,'replace',revoke_after_first)
    with pytest.raises(ValueError,match='revoked'):materialize.publish_byte_artifacts({a:b'first',b:b'second'},lock_path=lock)
    require(b.read_bytes()==b'old b','revoked writer continued to the second output')
    require(a.read_bytes()==(b'new owner bytes' if intervening_bytes else b'first'),'revoked writer overwrote intervening bytes or rolled back without authority')
    require((lock/'owner').read_bytes()==b'new owner','old writer removed the newer lease')


def test_shared_publication_rollback_never_overwrites_intervening_bytes(tmp_path,monkeypatch):
    tmp_path=tmp_path.resolve();a=tmp_path/'a';b=tmp_path/'b';a.write_bytes(b'old a');b.write_bytes(b'old b');original=os.replace
    def external_write_before_second_failure(source,destination):
        if Path(destination)==b:
            a.write_bytes(b'external bytes');raise OSError('injected second-file failure')
        return original(source,destination)
    monkeypatch.setattr(materialize.os,'replace',external_write_before_second_failure)
    with pytest.raises(ValueError,match='Rollback refused'):materialize.publish_byte_artifacts({a:b'first',b:b'second'})
    require(a.read_bytes()==b'external bytes' and b.read_bytes()==b'old b','rollback destroyed intervening output')


def test_publication_lease_detects_revocation_and_preserves_new_lock(tmp_path):
    path=tmp_path.resolve()/'lease'
    with export._publication_lock(path) as lease:
        shutil.rmtree(path);path.mkdir();(path/'owner').write_bytes(b'new owner')
        with pytest.raises(ValueError,match='revoked'):lease.require_owned()
    require((path/'owner').read_bytes()==b'new owner','lease cleanup removed the new owner')


@pytest.mark.parametrize('cache',['__pycache__/private-dossier.json','.pytest_cache/private-dossier.json','nested/__pycache__/source.txt'])
def test_public_inventory_has_no_cache_directory_payload_exemption(tmp_path,cache):
    (tmp_path/'payload').write_bytes(b'x');(tmp_path/'.openislet_public_bundle').write_bytes(b'marker')
    manifest={'artifacts':{'payload':{'sha256':hashlib.sha256(b'x').hexdigest(),'bytes':1}},'public_export_file_inventory':['payload']}
    require(export.manifest_integrity_errors(tmp_path,manifest,True)==[],'ordinary inventory failed')
    path=tmp_path/cache;path.parent.mkdir(parents=True);path.write_bytes(b'private source')
    require('unexpected:'+cache in export.manifest_integrity_errors(tmp_path,manifest,True),'undeclared cache payload passed exact inventory')


def opaque(quote,field='verbatim_evidence_quote',pmid='12345678'):
    return f'[REDACTED_NON_OA_EXCERPT: see PMID:{pmid} field={field} quote_sha256={hashlib.sha256(quote.encode()).hexdigest()}]'


def test_opaque_tokens_bind_source_field_and_graph_hop_identity():
    quote='alpha beta gamma delta epsilon zeta'
    good={'source_span':f'abstract[0:{len(quote)}]|norm=exact','verbatim_evidence_quote':opaque(quote)}
    export._validate_imported_redaction_tokens(good,'12345678',{'abstract':quote})
    bad=copy.deepcopy(good);bad['verbatim_evidence_quote']=opaque(quote,'focal_factor_X')
    with pytest.raises(ValueError):export._validate_imported_redaction_tokens(bad,'12345678',{'abstract':quote})
    graph={'hop1_study':'PMID:12345678','hop1_evidence':opaque(quote,'hop1_evidence','87654321')}
    require(bool(export.public_token_errors(graph)),'wrong hop PMID passed opaque reference validation')
    require(bool(export.public_token_errors({'note':'prefix '+opaque(quote,'note')})),'embedded marker prefix passed as an opaque token')


@pytest.mark.parametrize('text',['α β γ δ ε ζ','Ａ Ｂ Ｃ Ｄ Ｅ Ｆ','a\u034flpha beta gamma delta epsilon zeta','a\ufe0flpha beta gamma delta epsilon zeta'])
def test_unicode_source_overlap_counts_original_characters(text):
    source='α β γ δ ε ζ A B C D E F alpha beta gamma delta epsilon zeta'
    require(export._max_verbatim_shingle_chars({'note':text},source)==len(text),'invisible/compatibility Unicode bypassed source accounting')


def test_title_copied_into_free_text_still_consumes_article_budget():
    quote='alpha beta gamma delta epsilon zeta'
    require(export._max_verbatim_shingle_chars({'study_metadata':{'title':quote},'note':quote},quote)==len(quote),'free-text title replay was mistaken for citation metadata')


@pytest.mark.parametrize('literal',[
    '> ```mermaid\n> A["GSG-tPT2A"] --> B["MVDGTLLLL"]\n> ```',
    '  > > ~~~text\n> > GSG-tPT2A MVDGTLLLL\n> > ~~~~',
    '    GSG-tPT2A MVDGTLLLL',
    '\tGSG-tPT2A MVDGTLLLL',
])
def test_nested_and_indented_literals_survive_repeat_render(literal):
    rendered=renderer.render_text('# Test\n\n'+literal+'\n',facts(),'README.md')
    require(literal in rendered,'authored literal bytes were changed')
    require(renderer.render_text(rendered,facts(),'README.md')==rendered,'literal shielding is not idempotent')


@pytest.mark.parametrize('heading',['# Test\n',''])
def test_unclosed_fences_do_not_capture_new_generated_blocks(heading):
    literal='```mermaid\nA["GSG-tPT2A"] --> B["MVDGTLLLL"]\n'
    rendered=renderer.render_text(heading+literal,facts(),'README.md')
    require(literal.rstrip('\n') in rendered,'unclosed authored literal changed')
    require(renderer.render_text(rendered,facts(),'README.md')==rendered,'unclosed fence accumulated generated scope on repeat render')


def test_indented_literal_separator_blank_lines_survive_repeat_render():
    text='# Test\n\n    GSG-tPT2A MVDGTLLLL\n\n\n### Next section\n'
    rendered=renderer.render_text(text,facts(),'README.md')
    require('    GSG-tPT2A MVDGTLLLL\n\n\n### Next section' in rendered,'literal separator blank lines changed')
    require(renderer.render_text(rendered,facts(),'README.md')==rendered,'indented literals consumed a separator on every render')


def test_current_authored_documents_are_idempotent_in_memory():
    current_facts=facts()
    for path in [ROOT/'README.md',ROOT/'CLAUDE.md',ROOT/'ISLET-1-BLUEPRINT.md',*sorted((ROOT/'docs').glob('*.md'))]:
        if not path.exists() or path.name in {'README.md','PROJECT-DECISIONS-LOG.md','OPENISLET-2026-2030-MASTER-TIMELINE.md'}:continue
        first=renderer.render_text(path.read_text(),current_facts,path.name)
        require(renderer.render_text(first,current_facts,path.name)==first,'current document not idempotent: '+path.name)


def test_literal_placeholder_identity_survives_two_digit_indexes():
    literals=[f'```text\nunique literal {index}: GSG-tPT2A MVDGTLLLL\n```' for index in range(24)]
    text='# Literal provenance\n\n'+'\n\n'.join(literals)+'\n'
    rendered=renderer.render_text(text,facts(),'README.md')
    require(all(rendered.count(literal)==1 for literal in literals),'one literal placeholder captured another indexed literal')


def test_placeholder_like_authored_text_is_preserved_without_recursive_expansion():
    literal='```text\nOPENISLET_LITERAL_FENCE_0 OPENISLET_LITERAL_FENCE_10\n```'
    rendered=renderer.render_text('# Literal names\n\n'+literal+'\n',facts(),'README.md')
    require(rendered.count(literal)==1,'literal marker-like authored bytes were rejected or recursively expanded')


def test_duplicate_internal_literal_marker_fails_before_restoration(monkeypatch):
    original=renderer.facts_block
    monkeypatch.setattr(renderer,'facts_block',lambda value:original(value)+'\x00OPENISLET_LITERAL_FENCE_0\x00')
    with pytest.raises(ValueError,match='placeholder provenance'):
        renderer.render_text('# Literal provenance\n\n```text\nretained literal\n```\n',facts(),'README.md')


def test_real_atlas_authored_body_bytes_survive_rendering():
    path=ROOT/'docs/ISLET-1-MASTER-SPEC-v0.1.md'
    def authored_digest(text):
        regions=[(renderer.START,renderer.END),('<!-- OPENISLET_CURRENT_LAYER_SCOPE -->','<!-- OPENISLET_CURRENT_LAYER_SCOPE_END -->')]
        for start,end in regions:text=re.sub(re.escape(start)+'.*?'+re.escape(end),'[GENERATED REGION]',text,flags=re.S)
        return hashlib.sha256(text.encode('utf-8')).hexdigest()
    # Independent entry-baseline oracle: authored rules, equations, gene values
    # and the fenced directory tree must not become trusted corrupted input.
    expected='d56a70021d9dcf39405b5cc4594d09779c947e04584c8cd01456502a4f2412ec'
    source=path.read_text()
    require(authored_digest(source)==expected,'real atlas input no longer matches its protected authored-byte oracle')
    rendered=renderer.render_text(source,facts(),path.name)
    require(authored_digest(rendered)==expected,'rendering changed protected atlas body bytes')


def test_table_pipe_escape_parity_and_corruption_guard():
    require(renderer._table_cells(r'| one\\| two |')==['',' one\\\\',' two ',''],'even backslashes incorrectly hid a table delimiter')
    require(renderer._table_cells(r'| one\| two |')==['',' one\\| two ',''],'escaped pipe became a delimiter')
    text='# Test\n| **Layer 2:** '+'x'*600+' | authored cell |\n'
    with pytest.raises(ValueError):renderer.render_text(text,facts(),'README.md')


def test_document_batch_second_replace_restores_original_files(tmp_path,monkeypatch):
    a=tmp_path/'a.txt';b=tmp_path/'b.txt';a.write_bytes(b'old a');b.write_bytes(b'old b')
    original=os.replace
    def fail_second(source,destination):
        if Path(destination)==b:raise OSError('injected second-document failure')
        return original(source,destination)
    monkeypatch.setattr(renderer.os,'replace',fail_second)
    with pytest.raises(OSError):renderer._publish_document_batch({a:'new a',b:'new b'})
    require(a.read_bytes()==b'old a' and b.read_bytes()==b'old b','failed batch left a partially updated document set')


def test_document_targets_reject_ancestor_symlinks(tmp_path):
    outside=tmp_path/'outside';outside.mkdir();(outside/'a.txt').write_text('original')
    alias=tmp_path/'alias';alias.symlink_to(outside,target_is_directory=True)
    with pytest.raises(ValueError):renderer._publish_document_batch({alias/'a.txt':'changed'})
    require((outside/'a.txt').read_text()=='original','renderer followed a parent alias into another tree')


def test_document_facts_hash_the_bytes_actually_parsed_and_reject_later_change(tmp_path,monkeypatch):
    if not (ROOT/'data/extracted/islet1_candidate_stack_rankings.json').exists(): pytest.skip('Precomputed ranking cache intentionally excluded from public export')
    root=tmp_path.resolve()/'source';root.mkdir()
    files=['data/raw/corpus_manifest_all_sources.jsonl','data/raw/layer2_genomic_and_jaspar_reference.json',
           'data/extracted/cross_study_discovery_graph.json','data/extracted/islet1_candidate_stack_rankings.json',
           'data/protocols/hip_core_fractional_factorial_v0.6.2.json',
           *[f'data/extracted/layer{i}_{name}.json' for i,name in [(2,'regulatory_and_perturbation_atlas'),(3,'structural_verification_atlas'),(4,'tissue_and_metabolic_twin_atlas'),(5,'personalized_patient_and_crispr_atlas')]]]
    for relative in files:
        path=root/relative;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/relative,path)
    graph=root/'data/extracted/cross_study_discovery_graph.json';original=graph.read_bytes()
    read_regular=renderer._read_regular
    def mutate_after_read(path):
        data=read_regular(path)
        if path==graph:graph.write_bytes(b'{"changed_after_parse":true}\n')
        return data
    monkeypatch.setattr(renderer,'_read_regular',mutate_after_read)
    derived=renderer.derive_facts(root,1)
    require(derived['source_sha256'][graph.relative_to(root).as_posix()]==hashlib.sha256(original).hexdigest(),'facts digest authenticates a later read rather than parsed source bytes')
    graph.write_bytes(original)
    with pytest.raises(ValueError,match='source changed'):renderer.write_documentation(root,1,False)
    require(not (root/'data/extracted/documentation_facts_v0.8.0.json').exists(),'changed source facts reached publication')


def test_manifest_generation_failure_restores_previous_outputs(tmp_path,monkeypatch):
    root=tmp_path.resolve()/'source';root.mkdir();(root/'README.md').write_bytes(b'original prose')
    os.chmod(root/'README.md',0o640);monkeypatch.setattr(export,'ROOT',root)
    def fail_after_writes():
        path=root/'data/extracted/ring1_atomic_effects_v0.5.jsonl'
        export._publish_generation_artifacts({root/'README.md':b'partially rendered',path:b'new effects'},root/'.test-generation.lock')
        raise RuntimeError('injected manifest generation failure')
    monkeypatch.setattr(export,'_generate_release_manifest',fail_after_writes)
    with pytest.raises(RuntimeError):export.main()
    require((root/'README.md').read_bytes()==b'original prose' and (root/'README.md').stat().st_mode & 0o777==0o640,'failed generator did not restore bytes and permissions')
    require(not (root/'data/extracted/ring1_atomic_effects_v0.5.jsonl').exists() and not (root/'.openislet-manifest.lock').exists(),'failed generator retained new payload or lock')


@pytest.mark.parametrize('previously_owned',[False,True])
def test_manifest_generation_failure_preserves_separate_cooperating_publication(tmp_path,monkeypatch,previously_owned):
    root=tmp_path.resolve();readme=root/'README.md';readme.write_bytes(b'initial README')
    own=root/'own.json';own.write_bytes(b'initial own');monkeypatch.setattr(export,'ROOT',root)
    def fail_after_independent_publication():
        payloads={own:b'generated own'}
        if previously_owned:payloads[readme]=b'generated README'
        export._publish_generation_artifacts(payloads,root/'.test-generation.lock')
        materialize.publish_byte_artifacts({readme:b'independent rendered README'},lock_path=root/'.independent-documentation.lock')
        raise RuntimeError('injected later generation failure')
    monkeypatch.setattr(export,'_generate_release_manifest',fail_after_independent_publication)
    with pytest.raises(ValueError if previously_owned else RuntimeError):export.main()
    require(readme.read_bytes()==b'independent rendered README','macro rollback overwrote independent publication')
    require(own.read_bytes()==b'initial own','macro rollback failed to restore its unchanged owned output')
    require(not list(root.glob('.*lock')),'macro rollback leaked a cooperating writer lease')


def test_manifest_journal_tracks_repeated_writes_and_child_rollback(tmp_path,monkeypatch):
    root=tmp_path.resolve();a=root/'a';b=root/'b';a.write_bytes(b'initial');b.write_bytes(b'initial b')
    monkeypatch.setattr(export,'ROOT',root);original=os.replace
    def fail_after_repeated_writes():
        export._publish_generation_artifacts({a:b'first'},root/'.child.lock')
        export._publish_generation_artifacts({a:b'second'},root/'.child.lock')
        def fail_b(source,destination):
            if Path(destination)==b:raise OSError('injected child second commit failure')
            return original(source,destination)
        with monkeypatch.context() as patch:
            patch.setattr(materialize.os,'replace',fail_b)
            export._publish_generation_artifacts({a:b'third',b:b'changed b'},root/'.child.lock')
    monkeypatch.setattr(export,'_generate_release_manifest',fail_after_repeated_writes)
    with pytest.raises(OSError):export.main()
    require(a.read_bytes()==b'initial' and b.read_bytes()==b'initial b','journal lost original preimage after a later child rollback')


def test_manifest_generation_revocation_preserves_new_owner_bytes(tmp_path,monkeypatch):
    root=tmp_path.resolve()/'source';root.mkdir();(root/'README.md').write_bytes(b'original')
    monkeypatch.setattr(export,'ROOT',root)
    def revoke_during_generation():
        lock=root/'.openislet-manifest.lock';shutil.rmtree(lock);lock.mkdir();(lock/'owner').write_bytes(b'new owner')
        (root/'README.md').write_bytes(b'new owner bytes')
    monkeypatch.setattr(export,'_generate_release_manifest',revoke_during_generation)
    with pytest.raises(ValueError,match='revoked'):export.main()
    require((root/'README.md').read_bytes()==b'new owner bytes' and (root/'.openislet-manifest.lock/owner').read_bytes()==b'new owner','generator rollback or cleanup overwrote newer ownership')


@pytest.fixture
def export_source(tmp_path,monkeypatch):
    if not (ROOT/'data/dossiers/41385368.json').exists(): pytest.skip('Private source dossier intentionally excluded from public export')
    root=tmp_path.resolve()/'source';root.mkdir()
    files=['scripts/rank_islet1_stacks.py','scripts/engineering_priors.py','scripts/public_excerpt_policy.py','schemas/extraction_schema_v0.2.json','data/extracted/41385368.json','data/dossiers/41385368.json','data/omics/target_prior_table.json']
    for relative in files:
        target=root/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/relative,target)
    for relative,payload in [('LICENSE','license'),('data/extracted/cross_study_discovery_graph.json','{}'),('data/extracted/ring1_study_dossiers_v0.5.jsonl',''),('data/extracted/ring1_atomic_effects_v0.5.jsonl','')]:
        (root/relative).write_text(payload);files.append(relative)
    priv_pat = root / 'data/raw/private_export_guard_patterns.json'
    priv_pat.parent.mkdir(parents=True, exist_ok=True)
    priv_pat.write_text('{"regex_patterns":{"synthetic_guard":"SYNTHETIC_GUARD_TOKEN"}}')
    manifest={'artifacts':{},'public_redistribution_allowlist_spec':{'public_package_included_paths':[p for p in files if not p.startswith('data/dossiers/')]+['manifest.json'],'local_research_cache_excluded_from_public_bundle':['data/dossiers/*.json']}}
    def seal():
        for relative in files:
            data=(root/relative).read_bytes();manifest['artifacts'][relative]={'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}
        (root/'manifest.json').write_text(json.dumps(manifest))
    seal();monkeypatch.setattr(export,'ROOT',root)
    return root,manifest,files,seal


def test_public_export_nominal_has_exact_inventory_and_aliases(export_source,tmp_path):
    root,_,_,_=export_source;destination=tmp_path.resolve()/'public';export.export_public_bundle(destination)
    manifest=json.loads((destination/'manifest.json').read_text())
    require(not export.manifest_integrity_errors(destination,manifest,True),'export shipped undeclared or changed files')
    record=json.loads((destination/'data/extracted/41385368.json').read_text())
    require(json.loads((destination/'data/extracted/ring1_study_dossiers_v0.5.jsonl').read_text())==record,'public study alias diverges')
    require(not (destination/'data/dossiers').exists(),'private source directory shipped')


@pytest.mark.parametrize('caller_poison',['absent','cwd-file','preloaded-module'])
def test_public_export_file_spec_from_external_cwd_uses_own_sibling_code(export_source,tmp_path,caller_poison):
    root,_,_,_=export_source
    for relative in ['scripts/generate_release_manifest.py','scripts/study_artifact_exports.py']:
        shutil.copy2(ROOT/relative,root/relative)
    caller=tmp_path.resolve()/'unrelated-caller';caller.mkdir()
    if caller_poison!='absent':(caller/'study_artifact_exports.py').write_text("raise RuntimeError('caller-supplied sibling executed')\n")
    destination=tmp_path.resolve()/'external-public'
    code="""import importlib.util,json,pathlib,sys,types
before=list(sys.path)
source=pathlib.Path(sys.argv[1]);destination=pathlib.Path(sys.argv[2])
poison=None
if sys.argv[3]=='preloaded-module':
    poison=types.ModuleType('study_artifact_exports')
    def rejected(*args,**kwargs):raise RuntimeError('caller-supplied preloaded sibling used')
    poison.PublicationLease=poison.PublicationJournal=poison.publish_byte_artifacts=rejected
    sys.modules['study_artifact_exports']=poison
spec=importlib.util.spec_from_file_location('externally_loaded_release',source)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
files=module.export_public_bundle(destination)
if sys.path!=before:raise RuntimeError('exporter changed caller import search paths')
if poison is not None and sys.modules['study_artifact_exports'] is not poison:raise RuntimeError('exporter replaced caller module registry')
manifest=json.loads((destination/'manifest.json').read_text())
if module.manifest_integrity_errors(destination,manifest,True):raise RuntimeError('file-spec export failed exact inventory verification')
if (destination/'data/dossiers').exists():raise RuntimeError('file-spec export shipped private dossiers')
print(json.dumps({'files':len(files),'identity':json.loads((destination/'data/extracted/41385368.json').read_text())['study_metadata']['pmid']}))
"""
    flags=['-O'] if sys.flags.optimize else []
    environment={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'};environment.pop('PYTHONPATH',None)
    result=subprocess.run([sys.executable,'-B',*flags,'-c',code,str(root/'scripts/generate_release_manifest.py'),str(destination),caller_poison],cwd=caller,env=environment,capture_output=True,text=True)
    require(result.returncode==0,'external file-spec export failed: '+result.stderr)
    receipt=json.loads(result.stdout)
    require(receipt['files']==11 and receipt['identity']=='41385368','external export did not produce the expected minimal valid bundle')


def test_public_export_canonicalizes_approved_system_alias_before_lease(export_source,tmp_path):
    destination=tmp_path.resolve()/'public'
    # /var/../.. resolves to / on macOS as well as on ordinary Unix layouts;
    # all created files remain under this test's supplied evidence directory.
    alias=Path('/var/../..')/destination.relative_to('/')
    require(alias.resolve()==destination,'system alias fixture did not reach its intended destination')
    export.export_public_bundle(alias)
    require((destination/'.openislet_public_bundle').exists(),'approved system alias was rejected by the canonical-path lease')


def test_public_export_special_allowlist_match_fails_before_publication(export_source,tmp_path):
    root,manifest,_,_=export_source;os.mkfifo(root/'scripts/injected.py')
    manifest['public_redistribution_allowlist_spec']['public_package_included_paths'].append('scripts/*.py');(root/'manifest.json').write_text(json.dumps(manifest))
    destination=tmp_path.resolve()/'public'
    with pytest.raises(ValueError):export.export_public_bundle(destination)
    require(not destination.exists(),'special allowlist file reached publication')


def test_public_export_policy_cannot_reenable_private_dossiers(export_source,tmp_path):
    root,manifest,_,_=export_source
    manifest['public_redistribution_allowlist_spec']['public_package_included_paths'].append('data/dossiers/*.json')
    manifest['public_redistribution_allowlist_spec']['local_research_cache_excluded_from_public_bundle']=[]
    (root/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError):export.export_public_bundle(tmp_path.resolve()/'public')


def test_public_export_article_identity_must_match_filename(export_source,tmp_path):
    root,_,_,seal=export_source;path=root/'data/extracted/41385368.json';record=json.loads(path.read_text());record['study_metadata']['pmid']='12345678';path.write_text(json.dumps(record));seal()
    with pytest.raises(ValueError):export.export_public_bundle(tmp_path.resolve()/'public')


def test_public_export_reads_policy_inputs_from_pinned_staging_snapshot(export_source,tmp_path,monkeypatch):
    root,_,_,_=export_source
    rank_path=root/'scripts/rank_islet1_stacks.py';omics_path=root/'data/omics/target_prior_table.json'
    originals={path:path.read_bytes() for path in [rank_path,omics_path]}
    original_spec=importlib.util.spec_from_file_location
    original_validate=export._validate_imported_redaction_tokens
    changed=False;restored=False
    def spec_from_snapshot(name,path,*args,**kwargs):
        nonlocal changed
        if name=='openislet_rank_helper':
            changed=True;rank_path.write_text("raise RuntimeError('mutable original script executed')\n");omics_path.write_text('not JSON')
        return original_spec(name,path,*args,**kwargs)
    def restore_before_validation(*args,**kwargs):
        nonlocal restored
        if changed and not restored:
            for path,data in originals.items():path.write_bytes(data)
            restored=True
        return original_validate(*args,**kwargs)
    monkeypatch.setattr(importlib.util,'spec_from_file_location',spec_from_snapshot)
    monkeypatch.setattr(export,'_validate_imported_redaction_tokens',restore_before_validation)
    try:export.export_public_bundle(tmp_path.resolve()/'public')
    finally:
        for path,data in originals.items():path.write_bytes(data)
    require(changed and restored,'snapshot counterexample did not exercise the mutable-source interval')


@pytest.mark.parametrize('field',['unattributed_note','title','citation_short','journal','doi','pmid','pmcid','study','hop1_study','hop2_study','schema_version'])
@pytest.mark.parametrize('attributed',[False,True])
def test_public_export_redacts_unattributed_graph_source_text(export_source,tmp_path,field,attributed):
    root,_,_,seal=export_source
    source=json.loads((root/'data/dossiers/41385368.json').read_text())['abstract']
    graph={field:source}
    if attributed:graph={'study':'PMID:41385368','payload':graph}
    totals=export.public_article_verbatim_totals([],graph,{'41385368':(source,240)})
    require(totals['41385368']==len(source),'independent recount exempted source prose by its metadata key')
    (root/'data/extracted/cross_study_discovery_graph.json').write_text(json.dumps(graph));seal()
    destination=tmp_path.resolve()/'public';export.export_public_bundle(destination)
    output=json.loads((destination/'data/extracted/cross_study_discovery_graph.json').read_text())
    if attributed:output=output['payload']
    require(output[field]!=source and excerpts.validate_token(output[field],expected_pmid='41385368'),'unattributed graph prose bypassed article accounting')


def test_public_export_preserves_bounded_graph_metadata_identifiers(export_source,tmp_path):
    root,_,_,seal=export_source
    graph={'pmid':'41385368','pmcid':'PMC1234567','study':'PMID:41385368 (Synthetic A 2025)',
           'hop1_study':'PMID:41385368','hop2_study':'PMID:12345678 (Other B 2024)','schema_version':'0.2.0'}
    for key,value in graph.items():require(excerpts.is_bounded_metadata(value,key),'legitimate identifier failed shared policy validation')
    (root/'data/extracted/cross_study_discovery_graph.json').write_text(json.dumps(graph));seal()
    destination=tmp_path.resolve()/'public';export.export_public_bundle(destination)
    output=json.loads((destination/'data/extracted/cross_study_discovery_graph.json').read_text())
    require(all(output[key]==value for key,value in graph.items()),'bounded graph identifier was changed by rights enforcement')


@pytest.mark.parametrize('field,value',[
    ('pmid',' 41385368'),('pmid','41385368 trailing prose'),('pmid','0'),('pmcid','PMC1234567\n'),
    ('pmcid','PMC'+'1'*11),('study','PMID:41385368 (arbitrary prose)'),
    ('hop1_study','PMID:41385368 '+('x'*170)),('hop2_study','PMID:41385368 (Author 2025)\n'),
    ('schema_version','0.2.0 arbitrary prose'),('schema_version','1234.2.0'),
    ('doi','10.1234/'+('x'*250)),('doi','10.1234/private abstract'),
])
def test_metadata_exemptions_require_bounded_whole_value_syntax(field,value):
    require(not excerpts.is_bounded_metadata(value,field),'malformed or unbounded metadata received a prose exemption')


def test_public_export_cannot_hide_source_prose_in_graph_keys(export_source,tmp_path):
    root,_,_,seal=export_source;source=json.loads((root/'data/dossiers/41385368.json').read_text())['abstract']
    (root/'data/extracted/cross_study_discovery_graph.json').write_text(json.dumps({source:0}));seal()
    with pytest.raises(ValueError,match='field name'):export.export_public_bundle(tmp_path.resolve()/'public')


def test_public_export_discovery_metadata_cannot_ship_new_fulltext_field():
    source={'pmid':'12345678','title':'Synthetic title','abstract':'private abstract','fulltext':'private article','attachment':{'text':'private supplement'}}
    output=export._public_discovery_record(source)
    require(output['title']=='Synthetic title' and not {'abstract','fulltext','attachment'}&output.keys(),'future source-text channel escaped discovery metadata allowlist')


def test_combined_article_budget_counts_study_and_both_graph_hops():
    quote='alpha beta gamma delta epsilon zeta';source=quote+' '+quote
    first_pmid='12345678';second_pmid='87654321'
    studies=[{'study_metadata':{'pmid':first_pmid},'note':quote}]
    graph=[{'study':f'PMID:{first_pmid}','note':quote},{'hop1_study':f'PMID:{first_pmid}','hop1_evidence':quote,'hop2_study':f'PMID:{second_pmid}','hop2_evidence':quote}]
    totals=export.public_article_verbatim_totals(studies,graph,{first_pmid:(source,240),second_pmid:(source,240)})
    require(totals=={first_pmid:3*len(quote),second_pmid:len(quote)},'per-article accounting reset at study/graph or hop boundaries')


def test_unbound_graph_recount_normalizes_each_source_once_without_changing_totals(monkeypatch):
    quote='alpha beta gamma delta epsilon zeta'
    sources={str(10000000+i):(f'unique{i} '+quote,10000) for i in range(24)}
    source_texts={source for source,_ in sources.values()};calls={source:0 for source in source_texts};original=export._source_words
    def counted_words(text):
        if text in calls:calls[text]+=1
        return original(text)
    monkeypatch.setattr(export,'_source_words',counted_words)
    totals=export.public_article_verbatim_totals([],[{'schema_version':quote} for _ in range(50)],sources)
    require(all(total==50*len(quote) for total in totals.values()),'cached unbound recount lost an article or repeated occurrence')
    require(all(count==1 for count in calls.values()),'corpus recount repeatedly normalizes sources beyond the process LRU capacity')


def test_public_publication_rename_failure_restores_old_bundle(export_source,tmp_path,monkeypatch):
    root,_,_,_=export_source;destination=tmp_path.resolve()/'public';destination.mkdir();(destination/'.openislet_public_bundle').write_text('marker');(destination/'old.txt').write_bytes(b'previous bundle')
    original=Path.rename
    def fail_publish(path,target):
        if '.staging-' in path.name and Path(target)==destination:raise OSError('injected publication failure')
        return original(path,target)
    monkeypatch.setattr(Path,'rename',fail_publish)
    with pytest.raises(OSError):export.export_public_bundle(destination)
    require((destination/'old.txt').read_bytes()==b'previous bundle','failed publication destroyed the previous bundle')
    require(not list(destination.parent.glob('.public.staging-*')) and not (destination.parent/'.public.publish-lock').exists(),'publication leaked staging or lock state')


def test_public_export_revocation_during_staging_preserves_destination_and_new_lock(export_source,tmp_path,monkeypatch):
    root,_,_,_=export_source;destination=tmp_path.resolve()/'public';destination.mkdir();(destination/'.openislet_public_bundle').write_text('marker');(destination/'old').write_bytes(b'old bundle')
    lock=destination.parent/'.public.publish-lock';original=shutil.copy2
    def revoke_during_copy(source,target,*args,**kwargs):
        result=original(source,target,*args,**kwargs)
        if Path(source)==root/'LICENSE':
            shutil.rmtree(lock);lock.mkdir();(lock/'owner').write_bytes(b'new owner')
        return result
    monkeypatch.setattr(shutil,'copy2',revoke_during_copy)
    with pytest.raises(ValueError,match='revoked'):export.export_public_bundle(destination)
    require((destination/'old').read_bytes()==b'old bundle' and (lock/'owner').read_bytes()==b'new owner','revoked export changed destination or removed newer lease')
    require(not list(destination.parent.glob('.public.staging-*')),'revoked staging directory leaked')
