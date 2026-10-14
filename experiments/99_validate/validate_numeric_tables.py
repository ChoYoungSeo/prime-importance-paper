#!/usr/bin/env python3
"""Independently compare mutable table entries with numeric_claims.py output.
No table/source writes. Design n/p fields are verified by input manifests,
not recomputed here. Matching displayed rounding does not certify timing provenance.
"""
import argparse,json,re,sys
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--snapshot',required=True);p.add_argument('--candidate',required=True);p.add_argument('--out',required=True);a=p.parse_args()
m=json.loads(Path(a.snapshot).read_text())['metrics'];c=Path(a.candidate);checks=[];errors=[]
def rows(file):
 out=[]
 for line in (c/'tables'/file).read_text().splitlines():
  if '&' in line and '\\\\' in line:
   out.append([x.strip().removesuffix('\\\\').strip() for x in line.split('&')])
 return out

def check(file,label,expected,actual):
 ok=len(expected)==len(actual)
 if ok:
  for e,v in zip(expected,actual):
   if e is None:continue
   ok=ok and e==v
 record={'table':file,'row':label,'expected':expected,'actual':actual,'match':ok};checks.append(record)
 if not ok:errors.append(record)

def f(x,d=3):return 'unavailable' if x is None else f'{x:.{d}f}'
def s(z,c,stat='median'):return z[c][stat]
def interval(v):return f"{f(v['median'])} [{f(v['q05'])}, {f(v['q95'])}]"
def unique(file,name):
 selected=[r for r in rows(file) if r[0]==name]
 if len(selected)!=1:raise ValueError(f'{file}: expected unique row {name!r}, found {len(selected)}')
 return selected[0][1:]

lowmap={'PRIME':'PRIME','Full local-linear PME':'Full-coordinate local-linear PME','PME-kNN':'PME-kNN','NCMD-PME':'NCMD-PME','Oracle local-linear PME':'Oracle-coordinate LL-PME'}
for method,label in lowmap.items():
 z=m['low']['methods'][method];rho=m['low']['subgroups'][method+'|0.9']
 e=[f(s(z,'half_l1','mean')),f(s(z,'half_l1')),f(s(rho,'half_l1','mean')),f(s(z,'elapsed_sec','mean'),2)]
 check('pme_estimator_diagnostics.tex',label,e,unique('pme_estimator_diagnostics.tex',label))

short=lambda x:x.replace(' screen + PME-kNN',' + PME-kNN')
for method,z in m['correlated']['methods'].items():
 label=short(method);e=[f(s(z,'half_l1','mean')),f(z['active_inclusion']['proportion']),f(s(z,'proxy_mass','mean'),4),f"{z['completed']}/{z['attempted']}"]
 check('pme_correlated_pipeline.tex',label,e,unique('pme_correlated_pipeline.tex',label))

lookup={};response=None
for r in rows('pme_scaling_pmax.tex'):
 if r[0] in ('OTL','Piston','Borehole'):response=r[0].lower()
 if response and len(r)==6:lookup[(response,r[1])]=r[2:]
for response in ('otl','piston','borehole'):
 for method in m['scaling']['methods']:
  z=m['scaling']['subgroups'][response+'|500|'+method];e=[f(s(z,'half_l1')),f(s(z,'elapsed_sec'),2),f(s(z,'selected_size'),1),f"{z['completed']}/{z['attempted']}"]
  check('pme_scaling_pmax.tex',response+' '+method,e,lookup.get((response,short(method)),[]))

for label,module,dataset in [('Concrete','concrete','engineering'),('Airfoil','observed','airfoil'),('Auto MPG','observed','auto_mpg')]:
 z=m[module]['groups'][dataset+'|PRIME'];full=m[module]['groups'][dataset+'|Full LL-PME'];cv=m[module]['cv'][dataset]
 e=[None,None,f(s(z,'selected_size'),0),f(s(z,'value_calls'),0),f(s(full,'value_calls'),0),f(cv['median_paired_ratio'])]
 check('observed_application_main.tex',label,e,unique('observed_application_main.tex',label))

rawmap={'cement':'Cement','slag':'Blast-furnace slag','ash':'Fly ash','water':'Water','superplastic':'Superplasticizer','coarseagg':'Coarse aggregate','fineagg':'Fine aggregate','age':'Age'}
for var,label in rawmap.items():
 v=m['concrete']['coordinates']['raw|PRIME'][var];check('concrete_raw_supplement.tex',label,[f(v['selection_frequency']),interval(v)],unique('concrete_raw_supplement.tex',label))
obsmap={'airfoil':{'Angle':'Angle of attack','Chord_Length':'Chord length','Displacement':'Displacement thickness','Frequency':'Frequency','Velocity':'Velocity'},'auto_mpg':{'acceleration':'Acceleration','cylinders':'Cylinders','displacement':'Displacement','horsepower':'Horsepower','model_year':'Model year','origin_2':'Origin: Europe','origin_3':'Origin: Japan','weight':'Weight'}}
for dataset,mapping in obsmap.items():
 for var,label in mapping.items():
  v=m['observed']['coordinates'][dataset+'|PRIME'][var];full=m['observed']['coordinates'][dataset+'|Full LL-PME'][var]
  check('observed_allocation_supplement.tex',dataset+' '+label,[f(v['selection_frequency'],2),interval(v),interval(full)],unique('observed_allocation_supplement.tex',label))

for ratio,z in m['twinning'].items():
 e=[f(z['query_count']['median'],0),f(z['selected_jaccard']['median']),f(z['active_recall']['median']),f"{f(z['speedup']['median'],2)} [{f(z['speedup']['q25'],2)}, {f(z['speedup']['q75'],2)}]"]
 check('twinning_diagnostic.tex',ratio,e,unique('twinning_diagnostic.tex',ratio))

out={'status':'PASS' if not errors else 'FAIL','rows_checked':len(checks),'checks':checks,'errors':errors,'note':'Rounding/entry parity only. Timing source validity and final-input status remain separate.'}
Path(a.out).parent.mkdir(parents=True,exist_ok=True);Path(a.out).write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ('checks','errors')},indent=2));sys.exit(bool(errors))
