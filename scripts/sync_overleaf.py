#!/usr/bin/env python3
"""Copy the reviewed R update to an existing Overleaf-linked local folder.

Requires an explicit target. Refuses to overwrite a diverged manuscript base.
No remote compilation or synchronization is asserted by this local operation.
"""
import argparse, hashlib, json, shutil
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--target',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[1]
target=args.target.resolve()
original=json.loads((root/'verification/original_sha256.json').read_text())
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
current=target/'manuscript.tex'
allowed={original['manuscript.tex'],sha(root/'manuscript/manuscript.tex')}
if not current.is_file() or sha(current) not in allowed:
    raise SystemExit('Target manuscript differs from both the preserved base and reviewed update. Reconcile changes before copying.')
# Back up the original source outside the live target before overwriting it.
backup=root/'verification/local/overleaf_before_R'
backup.mkdir(parents=True,exist_ok=True)
for name in ('manuscript.tex','README.md'):
    p=target/name
    if p.exists() and not (backup/name).exists(): shutil.copy2(p,backup/name)
files={root/'manuscript/manuscript.tex':'manuscript.tex',
       root/'manuscript/manuscript.pdf':'R_Updated_Manuscript.pdf',
       root/'overleaf/README.md':'README.md',
       root/'overleaf/R_Implementation_Update.md':'R_Implementation_Update.md'}
for directory in ('R','scripts','data','results','figures','legacy'):
    for p in (root/directory).rglob('*'):
        if p.is_file(): files[p]=str(p.relative_to(root))
for p in (root/'manuscript/generated').glob('*.tex'): files[p]=str(Path('generated')/p.name)
for p in (root/'verification').glob('*'):
    if p.is_file() and p.suffix in ('.json','.md'): files[p]=str(Path('verification')/p.name)
record={}
for src,rel in files.items():
    dst=target/rel; dst.parent.mkdir(parents=True,exist_ok=True)
    if not dst.exists() or sha(src)!=sha(dst): shutil.copy2(src,dst)
    if sha(src)!=sha(dst): raise RuntimeError(f'Copy mismatch: {rel}')
    record[rel]=sha(dst)
(root/'verification/local/overleaf_sync_receipt.json').write_text(json.dumps({
    'operation':'Local copy to existing Overleaf-linked folder',
    'files_verified':len(record),'sha256':record,
    'live_Overleaf_compilation_verified':False},indent=2)+'\n')
print(f'Copied and hash-verified {len(record)} files. Online compilation is a separate check.')
