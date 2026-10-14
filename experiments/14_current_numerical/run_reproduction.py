#!/usr/bin/env python3
"""Run unchanged current estimators on frozen inputs; default is two fixtures."""
import argparse, importlib.util, os, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('scope',choices=['smoke','all'],nargs='?',default='smoke')
a=p.parse_args()
if not (ROOT/'.verified/first_reference.R').is_file():
    p.error('Import the official FIRST code.zip first: python3 experiments/14_current_numerical/import_first.py /path/to/code.zip')
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
v=module('current',Path(__file__).with_name('validate.py'));v.validate()
w=module('verified',ROOT/'experiments/99_validate/verified_workflow.py')
m,h=w.verify_release(ROOT);w.package_snapshot(ROOT,m,h)
env=w.r_environment(ROOT)
for k in ['PRIME_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS']:env[k]='1'
subprocess.run(['Rscript','--vanilla',str(Path(__file__).with_name('run_reproduction.R')),a.scope],cwd=ROOT,env=env,check=True)
