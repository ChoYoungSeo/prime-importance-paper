#!/usr/bin/env python3
"""Import the hash-pinned official FIRST supplement into an ignored local cache."""
import argparse, hashlib, zipfile
from pathlib import Path
EXPECTED='759c110ef16e4cce1fcebf17e1d1ed03035a90c92878e02715ff038f1afe1ee4'
p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path,help='Official code.zip or extracted code/lib/first.R');a=p.parse_args()
if zipfile.is_zipfile(a.source):
    with zipfile.ZipFile(a.source) as z:
        matches=[n for n in z.namelist() if n.endswith('/lib/first.R') or n=='lib/first.R']
        if len(matches)!=1:raise SystemExit('Expected exactly one lib/first.R in the official archive')
        data=z.read(matches[0])
else:data=a.source.read_bytes()
if hashlib.sha256(data).hexdigest()!=EXPECTED:raise SystemExit('FIRST source hash mismatch; source was not imported')
root=Path(__file__).resolve().parents[2];target=root/'.verified/first_reference.R';target.parent.mkdir(exist_ok=True);target.write_bytes(data)
print(f'Verified original FIRST reference: {target}')
