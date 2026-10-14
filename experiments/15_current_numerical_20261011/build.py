#!/usr/bin/env python3
"""Validate saved outputs and render the approved R figures; no fitting."""
import argparse
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'build/current_numerical_20261011')
    parser.add_argument('--check-only',action='store_true',help='Check R box statistics and layout without rendering images.')
    args=parser.parse_args()
    output=args.output_dir.resolve()
    for protected in ('manuscript','results','software','experiments'):
        path=ROOT/protected
        if output==path or path in output.parents:parser.error('Output must be outside scientific source and saved-result directories.')
    subprocess.run([sys.executable,str(HERE/'validate.py')],check=True)
    for dataset in ('friedman','concrete'):
        command=['Rscript','--vanilla',str(HERE/'build_figures.R'),'--dataset',dataset,'--output-dir',str(output)]
        if args.check_only:command.append('--check')
        subprocess.run(command,check=True,cwd=ROOT)
    print('Validated output:',output)

if __name__=='__main__':main()
