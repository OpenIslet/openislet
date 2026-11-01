"""Canonical study/effect JSONL materialization at the producer boundary."""
from __future__ import annotations
import json
import os
from pathlib import Path
import re
import stat
import tempfile
import unicodedata
import uuid


def _portable_path_key(path):
    # APFS/HFS can treat case and canonical Unicode spellings as one path even
    # when Path.resolve() preserves distinct strings. Requiring distinct keys
    # also keeps a multi-file publication portable to case-insensitive volumes.
    return Path(unicodedata.normalize('NFC',str(path)).casefold())


def _safe_destination(path):
    if '..' in path.parts:raise ValueError('Artifact destinations must be canonical paths')
    for parent in [path,*path.parents]:
        if parent.is_symlink():
            raise ValueError('Artifact destinations and ancestors may not be symbolic links')
    if path.exists() and (not stat.S_ISREG(path.lstat().st_mode) or path.stat().st_nlink!=1):
        raise ValueError('Artifact destinations must be single-link regular files')


def _file_snapshot(path):
    _safe_destination(path)
    try:info=path.lstat()
    except FileNotFoundError:return None
    fd=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
    try:
        opened=os.fstat(fd)
        if not stat.S_ISREG(opened.st_mode) or opened.st_nlink!=1 or (opened.st_dev,opened.st_ino)!=(info.st_dev,info.st_ino):
            raise ValueError('Artifact identity changed while opening')
        with os.fdopen(fd,'rb',closefd=False) as stream:data=stream.read()
        after=os.fstat(fd)
        identity=(opened.st_dev,opened.st_ino,opened.st_size,opened.st_mtime_ns,opened.st_ctime_ns)
        if (after.st_dev,after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)!=identity or len(data)!=opened.st_size:
            raise ValueError('Artifact changed while reading')
        return data,opened.st_mode & 0o777,identity
    finally:os.close(fd)


def _stage_bytes(path,data,mode):
    fd,name=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'wb') as stream:
            stream.write(data);stream.flush();os.fsync(stream.fileno())
        os.chmod(name,mode)
        return name
    except BaseException:
        if os.path.exists(name):os.unlink(name)
        raise


class PublicationLease:
    """Cooperative filesystem writer lease with ownership-aware cleanup."""
    def __init__(self,path):
        self.path=Path(path).absolute()
        if '..' in self.path.parts:raise ValueError('Publication lock must be canonical')
        for parent in [self.path,*self.path.parents]:
            if parent.is_symlink():raise ValueError('Publication lock and ancestors may not be symbolic links')
        self.path.parent.mkdir(parents=True,exist_ok=True)
        try:self.path.mkdir()
        except FileExistsError as error:raise ValueError('Another writer owns artifact publication') from error
        info=self.path.lstat();self.identity=(info.st_dev,info.st_ino)
        self.owner=self.path/'owner';self.nonce=uuid.uuid4().hex.encode('ascii')
        try:
            with self.owner.open('xb') as stream:
                stream.write(self.nonce);stream.flush();os.fsync(stream.fileno())
            info=self.owner.lstat();self.owner_identity=(info.st_dev,info.st_ino)
        except BaseException:
            # An incomplete owner record is never authority to remove an
            # unknown lock. Leave ambiguous interrupted state for inspection.
            if self.owner.exists() and self.owner.read_bytes()==self.nonce:self.owner.unlink()
            current=self.path.lstat() if self.path.exists() else None
            if current is not None and (current.st_dev,current.st_ino)==self.identity and not any(self.path.iterdir()):self.path.rmdir()
            raise

    def owns_lock(self):
        try:
            info=self.path.lstat();record=_file_snapshot(self.owner)
            return (stat.S_ISDIR(info.st_mode) and (info.st_dev,info.st_ino)==self.identity and record is not None
                    and record[0]==self.nonce and record[2][:2]==self.owner_identity)
        except (OSError,ValueError):return False

    def require_owned(self):
        if not self.owns_lock():raise ValueError('Artifact publication lease was revoked; refusing further writes')

    def close(self):
        if self.owns_lock():
            self.owner.unlink()
            current=self.path.lstat() if self.path.exists() else None
            if current is not None and (current.st_dev,current.st_ino)==self.identity:self.path.rmdir()

    def __enter__(self):return self

    def __exit__(self,*_):self.close()


class PublicationJournal:
    """Guarded rollback receipts for an explicitly enrolled sequence of writes.

    Independent publishers are not enrolled implicitly. A receipt authorizes
    rollback only while the enclosing authority and original writer lease hold
    and the file still has the exact snapshot produced by that publication.
    """
    def __init__(self):self.entries={}

    def record(self,saved,committed,lock):
        for path,current in committed.items():
            prior=self.entries.get(path)
            original=prior[0] if prior is not None and prior[1]==saved[path] else saved[path]
            self.entries[path]=(original,current,lock)

    def restored(self,path,previous,current):
        prior=self.entries.get(path)
        if prior is not None and prior[1]==previous:
            self.entries[path]=(prior[0],current,prior[2])

    def rollback(self,require_authority):
        conflicts=[]
        groups={}
        for path,entry in reversed(self.entries.items()):groups.setdefault(entry[2],[]).append((path,entry))
        for lock,entries in groups.items():
            require_authority()
            with PublicationLease(lock) as lease:
                for path,(previous,owned,_) in entries:
                    require_authority();lease.require_owned()
                    if _file_snapshot(path)!=owned:
                        conflicts.append(str(path));continue
                    if previous is None:
                        path.unlink()
                    else:
                        name=_stage_bytes(path,previous[0],previous[1])
                        try:
                            require_authority();lease.require_owned()
                            if _file_snapshot(path)!=owned:
                                conflicts.append(str(path));continue
                            os.replace(name,path)
                        finally:
                            if os.path.exists(name):os.unlink(name)
                    require_authority();lease.require_owned()
        if conflicts:raise ValueError('Generation rollback refused for externally changed artifact(s): '+', '.join(conflicts))


def publish_byte_artifacts(payloads: dict[Path, bytes], *, lock_path: Path | None = None,
                           rollback_journal: PublicationJournal | None = None) -> None:
    """Publish prepared bytes with writer locking and rollback on caught failure.

    Parents may differ; every temporary file is staged beside its destination.
    Detected lease revocation stops writes and preserves the newer owner, even
    if an earlier replacement already occurred. Rollback requires the same
    lease and unchanged bytes/identity of each committed output. This does not
    give atomic visibility to readers or survive a process/host crash between
    replacements; immutable generations with one reader pointer are required.
    """
    if not isinstance(payloads,dict) or not payloads or any(not isinstance(data,bytes) for data in payloads.values()):
        raise ValueError('Artifact payloads must be a nonempty mapping of paths to bytes')
    targets={Path(path).absolute():data for path,data in payloads.items()}
    keys={_portable_path_key(path) for path in targets}
    if len(targets)!=len(payloads) or len(keys)!=len(targets):raise ValueError('Artifact destinations must have unique portable identities')
    lock=Path(lock_path).absolute() if lock_path is not None else Path(os.path.commonpath([path.parent for path in targets]))/'.artifact-publication.lock'
    for path in targets:_safe_destination(path)
    lock_key=_portable_path_key(lock)
    if any(key==lock_key or key in lock_key.parents or lock_key in key.parents or any(key in other.parents for other in keys if other!=key) for key in keys):
        raise ValueError('Artifact destinations and lock must not overlap')
    if '..' in lock.parts:raise ValueError('Artifact lock must be canonical')
    for parent in [lock,*lock.parents]:
        if parent.is_symlink():raise ValueError('Artifact lock and ancestors may not be symbolic links')
    existing_identities=set()
    for path in targets:
        if path.exists():
            info=path.stat();identity=(info.st_dev,info.st_ino)
            if identity in existing_identities:raise ValueError('Artifact destinations alias one existing file')
            existing_identities.add(identity)
    lease=PublicationLease(lock)
    require_lock=lease.require_owned
    owns_lock=lease.owns_lock
    saved={};staged={};committed={}
    try:
        for path,data in targets.items():
            require_lock()
            previous=_file_snapshot(path)
            mode=previous[1] if previous is not None else 0o644
            saved[path]=previous
            path.parent.mkdir(parents=True,exist_ok=True)
            staged[path]=_stage_bytes(path,data,mode)
        for path,name in staged.items():
            require_lock()
            if _file_snapshot(path)!=saved[path]:raise ValueError('Artifact changed before publication; refusing to overwrite intervening bytes')
            os.replace(name,path)
            require_lock()
            observed=_file_snapshot(path)
            if observed is None or observed[0]!=targets[path]:raise ValueError('Artifact changed immediately after publication')
            committed[path]=observed
        if rollback_journal is not None:rollback_journal.record(saved,committed,lock)
    except BaseException as failure:
        conflicts=[]
        # A revoked lease no longer authorizes restoring earlier files. Another
        # owner may already have published intervening output bytes.
        if owns_lock():
            for path in reversed(committed):
                require_lock()
                if _file_snapshot(path)!=committed[path]:
                    conflicts.append(str(path));continue
                previous=saved[path]
                if previous is None:path.unlink()
                else:
                    name=_stage_bytes(path,previous[0],previous[1])
                    try:
                        require_lock();os.replace(name,path)
                    finally:
                        if os.path.exists(name):os.unlink(name)
                if rollback_journal is not None:rollback_journal.restored(path,previous,_file_snapshot(path))
        if conflicts:raise ValueError('Rollback refused for externally changed artifact(s): '+', '.join(conflicts)) from failure
        raise
    finally:
        for name in staged.values():
            if os.path.exists(name):os.unlink(name)
        lease.close()


def write_study_exports(root: Path, studies: list[dict], *, rollback_journal: PublicationJournal | None = None) -> tuple[Path,Path,int]:
    """Write the same exact row format consumed by downstream planes and exports.

    The extractor calls this after producing all study JSONs; the manifest uses
    this helper again. Thus Layer 2 cannot observe a previous run's effect cache.
    """
    root=Path(root).absolute()
    if not isinstance(studies, list) or any(not isinstance(study, dict) for study in studies):
        raise ValueError('Studies must be a list of JSON objects')
    effects=[]
    seen_pmids=set()
    for study in studies:
        meta=study.get('study_metadata')
        if not isinstance(meta, dict):
            raise ValueError('Study metadata must identify a study')
        pmid=str(meta.get('pmid',''))
        if re.fullmatch(r'[1-9][0-9]{0,8}',pmid) is None or pmid in seen_pmids:
            raise ValueError('Study PMID must be canonical and unique')
        seen_pmids.add(pmid)
        citation=meta.get('citation_short','')
        atomic_effects=study.get('atomic_effects',[])
        if not isinstance(atomic_effects,list) or any(not isinstance(effect,dict) for effect in atomic_effects):
            raise ValueError('Atomic effects must be a list of JSON objects')
        bound={'pmid':pmid,'citation_short':citation,
               'source_text_type':meta.get('source_text_type'),
               'superseded_by_pmid':meta.get('superseded_by_pmid')}
        for effect in atomic_effects:
            if any(key in effect and effect[key] != value for key,value in bound.items()):
                raise ValueError('Atomic effect may not replace its enclosing study identity')
            effects.append({**bound,**effect})
    studies_path=root/'data/extracted/ring1_study_dossiers_v0.5.jsonl'
    effects_path=root/'data/extracted/ring1_atomic_effects_v0.5.jsonl'
    payloads={path:''.join(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n' for row in rows).encode('utf-8')
              for path,rows in [(studies_path,studies),(effects_path,effects)]}
    publish_byte_artifacts(payloads,lock_path=studies_path.parent/'.study-exports.lock',rollback_journal=rollback_journal)
    return studies_path,effects_path,len(effects)
