#!/usr/bin/env python3
"""Read-only source-to-prose arithmetic and explicit draft/validation contract.

Never writes inputs or the candidate. --emit-proposal writes two review copies
under --out. Release validation additionally requires explicit final-input,
observed-stage and timing attestations; CSV arithmetic cannot certify methods.
"""
import argparse, csv, difflib, hashlib, json, math, re, statistics, sys
from collections import Counter, defaultdict
from pathlib import Path

NA = {'', 'NA', 'NaN', 'nan', 'NULL'}
TOL = 2e-8
MODULES = {
 'link': ('selection_to_pme_v052_main','selection_to_pme.csv','selection_to_pme_phi.csv',('link_id','method'),100),
 'low': ('pme_v052_main','pme_runs.csv','pme_phi.csv',('task_id','method'),None),
 'detail': ('friedman_detail_v052_main','friedman_detail_runs.csv','friedman_detail_phi.csv',('detail_id','method'),100),
 'correlated': ('correlated_pipelines_v052_main','correlated_pipeline_runs.csv','correlated_pipeline_phi.csv',('task_id','method'),None),
 'scaling': ('scaling_v052_main','scaling_runs.csv','scaling_phi.csv',('task_id','method'),None),
 'concrete': ('concrete_v052_main','bootstrap_runs.csv','bootstrap_phi.csv',('representation','replication','method'),8),
 'observed': ('observed_application_v052_main','bootstrap_runs.csv','bootstrap_phi.csv',('dataset','replication','method'),None),
}
PIPE = ['PRIME','FIRST screen + PME-kNN','NNVS screen + PME-kNN','LASSO screen + PME-kNN','FOCI screen + PME-kNN','FORD screen + PME-kNN']
LINK = ['PRIME screen','FIRST','FOCI','FORD','LASSO']
LOW = ['PRIME','Full local-linear PME','PME-kNN','NCMD-PME','Oracle local-linear PME']
SHORT = {x:x.replace(' screen + PME-kNN','') for x in PIPE}
SHORT.update({'PRIME screen':'PRIME screen','Full local-linear PME':'Full-coordinate LL-PME','Oracle local-linear PME':'Oracle-coordinate LL-PME'})

def number(v):
 if v is None or str(v).strip() in NA:return None
 x=float(v)
 return x if math.isfinite(x) else None

def boolean(v):
 if str(v).upper() in ('TRUE','1'):return True
 if str(v).upper() in ('FALSE','0'):return False
 raise ValueError('Missing/invalid logical value: '+repr(v))

def avg(v):return statistics.mean(v) if v else None
def median(v):return statistics.median(v) if v else None
def fmt(v,d=3):return 'unavailable' if v is None else f'{v:.{d}f}'
def quantile(v,p):
 if not v:return None
 v=sorted(v);h=(len(v)-1)*p;i=math.floor(h)
 return v[i]+(h-i)*(v[min(i+1,len(v)-1)]-v[i])
def key(r,cols):return tuple(r[c] for c in cols)
def grouped(rows,cols):
 out=defaultdict(list)
 for r in rows:out[key(r,cols)].append(r)
 return out

def raw_root(p):
 p=Path(p).resolve()
 return p/'results/raw' if (p/'results/raw').is_dir() else p

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def collapse(s):return re.sub(r'\s+',' ',s).strip()

def human_list(items):
 return items[0] if len(items)==1 else (' and '.join(items) if len(items)==2 else ', '.join(items[:-1])+', and '+items[-1])

def target_token(r,keys):
 t=number(r.get('truth'))
 payload=(key(r,keys),r['variable'],None if t is None else format(t,'.12g'),r.get('active'),r.get('proxy'),r.get('role'))
 return int.from_bytes(hashlib.sha256(repr(payload).encode()).digest(),'big')

class Audit:
 def __init__(self,args):
  self.a=args;self.root=raw_root(args.data_root);self.base=raw_root(args.baseline_root)
  self.fixed=raw_root(args.fixed_root or args.baseline_root)
  self.sources={};self.rows={};self.metrics={};self.errors=[];self.warnings=[];self.boot={};self.blocks=[]
 def demand(self,ok,msg):
  if not ok:self.errors.append(msg)
 def read(self,p,source_name,stream=False):
  if not p.is_file():raise ValueError(f'Missing required input: {p}')
  before=(p.stat().st_size,p.stat().st_mtime_ns,p.resolve())
  def iterator():
   with p.open(newline='') as f:
    yield from csv.DictReader(f)
   after=(p.stat().st_size,p.stat().st_mtime_ns,p.resolve())
   self.demand(before==after,'Input changed during read: '+str(p))
   self.sources[source_name]={'path':str(p),'resolved':str(p.resolve()),'sha256':sha(p),'bytes':p.stat().st_size}
  return iterator() if stream else list(iterator())
 def role_path(self,name,filename):
  folder=MODULES[name][0]
  p=self.root/folder/filename
  if p.is_file():return p
  # CV may legitimately be unchanged, but its inherited source is explicit.
  if filename=='nested_cv.csv':
   self.warnings.append(f'{name} CV inherited from fixed root; final observed-stage attestation required')
   return self.fixed/folder/filename
  return p
 def load_module(self,name):
  folder,runfile,phifile,keys,default_dim=MODULES[name]
  rows=self.read(self.role_path(name,runfile),name+'.runs')
  base=self.read(self.base/folder/runfile,'baseline.'+name+'.runs')
  by={key(r,keys):r for r in rows};base_by={key(r,keys):r for r in base}
  self.demand(len(by)==len(rows),name+': duplicate run keys')
  self.demand(set(by)==set(base_by),name+': missing/extra run keys versus frozen protocol baseline')
  immutable=['seed','benchmark','n','p','rho','design','replication','representation','dataset','method']
  for k,r in by.items():
   if k in base_by:
    for c in immutable:
     if c in r and c in base_by[k]:self.demand(r[c]==base_by[k][c],f'{name} {k}: changed protocol field {c}')
   r['_key']=k;r['_dim']=int(r.get('p') or default_dim)
   r['_ok']=boolean(r['feasible']) if 'feasible' in r else number(r.get('half_l1')) is not None
   r['_loss']=number(r.get('half_l1'));r['_time']=number(r.get('elapsed_sec'))
   if r['_ok']:
    if 'half_l1' in r:self.demand(r['_loss'] is not None,name+': completed run lacks finite loss '+str(k))
   elif 'half_l1' in r:self.demand(r['_loss'] is None,name+': incomplete run has a statistical loss '+str(k))
   if r.get('selected') is not None:
    sel=[int(x) for x in r['selected'].split(',') if x.strip()]
    self.demand(len(set(sel))==len(sel)==int(r['selected_size']),name+': selected-size mismatch '+str(k))
    r['_selected']=set(sel)
   if 'elapsed_sec' in r and r['_time'] is None:self.warnings.append(name+': missing time '+str(k))
  accum={k:{'count':0,'mask':0,'finite':0,'sum':0.,'truthsum':0.,'loss':0.,'inactive':0.,'proxy':0.,'has_truth':False} for k in by}
  baseline_target=0;current_target=0
  if name not in ('concrete','observed'):
   for r in self.read(self.base/folder/phifile,'baseline.'+name+'.phi',stream=True):baseline_target ^= target_token(r,keys)
  for r in self.read(self.role_path(name,phifile),name+'.phi',stream=True):
   if name not in ('concrete','observed'):current_target ^= target_token(r,keys)
   k=key(r,keys)
   if k not in by:
    self.errors.append(name+': orphan allocation row '+str(k));continue
   z=accum[k];v=int(r['variable']);dim=by[k]['_dim']
   self.demand(1<=v<=dim,name+': coordinate out of bounds '+str(k))
   bit=1 << (v-1)
   self.demand(not z['mask']&bit,name+': duplicate allocation coordinate '+str(k))
   z['mask']|=bit;z['count']+=1;x=number(r.get('estimate'));truth=number(r.get('truth'))
   if 'truth' in r:
    self.demand(truth is not None,name+': missing population truth '+str(k));z['has_truth']=True
    if truth is not None:z['truthsum']+=truth
   if x is not None:
    z['finite']+=1;z['sum']+=x
    if by[k]['_ok']:self.demand(-TOL<=x<=1+TOL,name+': allocation outside simplex '+str(k))
    if truth is not None:z['loss']+=abs(x-truth)/2
    active=boolean(r['active']) if 'active' in r else r.get('role')=='active'
    proxy=boolean(r['proxy']) if 'proxy' in r else r.get('role') in ('proxy','correlated proxy')
    if not active:z['inactive']+=x
    if proxy:z['proxy']+=x
   if name in ('concrete','observed'):
    self.boot.setdefault((name,k),{})[v]={'estimate':x,'name':r.get('variable_name',str(v)),'selected':boolean(r['selected'])}
    self.demand(boolean(r['selected'])==(v in by[k].get('_selected',set())),name+': phi/run selection mismatch '+str(k))
  if name not in ('concrete','observed'):self.demand(current_target==baseline_target,name+': changed population target or coordinate roles versus frozen baseline')
  for k,z in accum.items():
   r=by[k];self.demand(z['count']==r['_dim'] and z['mask'].bit_count()==r['_dim'],name+': incomplete coordinate coverage '+str(k))
   if z['has_truth']:self.demand(abs(z['truthsum']-1)<TOL,name+': truth not normalized '+str(k))
   if r['_ok']:
    if self.a.historical_summary and name=='link' and r['method']=='LASSO' and z['finite']==0:
     self.warnings.append(name+': HISTORICAL ONLY: finite run loss but missing entire coordinate vector '+str(k));continue
    self.demand(z['finite']==r['_dim'] and abs(z['sum']-1)<TOL,name+': completed allocation not finite/normalized '+str(k))
    if z['has_truth']:
     self.demand(r['_loss'] is not None and abs(r['_loss']-z['loss'])<TOL,name+': loss does not match coordinate allocations '+str(k))
     for col,field in [('inactive_mass','inactive'),('proxy_mass','proxy')]:
      if col in r and number(r[col]) is not None:self.demand(abs(float(r[col])-z[field])<TOL,name+': '+col+' mismatch '+str(k))
   elif z['finite']:
    self.warnings.append(name+': incomplete run retains finite allocation entries; excluded by explicit status '+str(k))
  self.rows[name]=rows
  self.metrics[name]={'methods':{k[0]:self.summary(v) for k,v in grouped(rows,['method']).items()}}
 def summary(self,rows):
  ok=[r for r in rows if r['_ok']];res={'attempted':len(rows),'completed':len(ok),'failed':len(rows)-len(ok),'reasons':dict(Counter(r.get('error') or r.get('reason') or 'unrecorded reason' for r in rows if not r['_ok']))}
  for c in ('half_l1','inactive_mass','proxy_mass','value_calls','selected_size','elapsed_sec','selection_elapsed','allocation_elapsed'):
   if not any(c in r for r in rows):continue
   selected=ok if c in ('half_l1','inactive_mass','proxy_mass') else rows
   v=[number(r.get(c)) for r in selected];v=[x for x in v if x is not None]
   res[c]={'n_finite':len(v),'basis':'completed' if selected is ok else 'all attempts with finite field','mean':avg(v),'median':median(v),'q05':quantile(v,.05),'q25':quantile(v,.25),'q75':quantile(v,.75),'q95':quantile(v,.95)}
  for c in ('active_inclusion','exact_active_set','exact'):
   if rows and c in rows[0]:res[c]={'count':sum(boolean(r[c]) for r in rows),'denominator':len(rows),'proportion':sum(boolean(r[c]) for r in rows)/len(rows)}
  return res
 def ranked_cells(self,name,cols):
  cells={};wins=Counter();ties=Counter()
  for k,rs in grouped(self.rows[name],cols).items():
   values={m[0]:avg([r['_loss'] for r in v if r['_ok']]) for m,v in grouped(rs,['method']).items()}
   finite={m:v for m,v in values.items() if v is not None}
   winners=[m for m,v in finite.items() if abs(v-min(finite.values()))<=1e-12] if finite else []
   for w in winners:(wins if len(winners)==1 else ties)[w]+=1
   cells['|'.join(k)]={'mean_completed_loss':values,'winners':winners,'completion':{m[0]:[sum(r['_ok'] for r in v),len(v)] for m,v in grouped(rs,['method']).items()}}
  self.metrics[name]['cells']=cells;self.metrics[name]['strict_wins']=dict(wins);self.metrics[name]['tied_wins']=dict(ties)
 def allocation_details(self):
  for name,cols in [('link',['benchmark','design']),('low',['benchmark','n','rho']),('correlated',['benchmark','n']),('scaling',['benchmark','p'])]:self.ranked_cells(name,cols)
  for name in ('low','scaling','correlated'):
   cols=['method','rho'] if name=='low' else (['benchmark','p','method'] if name=='scaling' else ['benchmark','n','method'])
   self.metrics[name]['subgroups']={'|'.join(k):self.summary(v) for k,v in grouped(self.rows[name],cols).items()}
  self.metrics['scaling']['prime_calls_by_p']={k[0]:self.summary(v)['value_calls'] for k,v in grouped([r for r in self.rows['scaling'] if r['method']=='PRIME'],['p']).items()}
  self.metrics['scaling']['prime_coverage_by_response']={k[0]:self.summary(v)['active_inclusion'] for k,v in grouped([r for r in self.rows['scaling'] if r['method']=='PRIME'],['benchmark']).items()}
  for name in ('concrete','observed'):
   rows=self.rows[name];field='representation' if name=='concrete' else 'dataset';groups=grouped(rows,[field,'method'])
   self.metrics[name]['groups']={'|'.join(k):self.summary(v) for k,v in groups.items()};coordinates={}
   for k,rs in groups.items():
    dim=rs[0]['_dim'];coord={}
    for j in range(1,dim+1):
     v=[self.boot[(name,r['_key'])][j] for r in rs];est=[x['estimate'] for x,r in zip(v,rs) if r['_ok'] and x['estimate'] is not None]
     coord[v[0]['name']]={'median':median(est),'q05':quantile(est,.05),'q95':quantile(est,.95),'completed':len(est),'attempted':len(rs),'selection_count':sum(x['selected'] for x in v),'selection_frequency':sum(x['selected'] for x in v)/len(v)}
    coordinates['|'.join(k)]=coord
   self.metrics[name]['coordinates']=coordinates
   pairs={}
   for k,rs in grouped(rows,[field,'replication']).items():
    if name=='concrete' and k[0]=='raw':continue
    methods={r['method']:r for r in rs};self.demand('PRIME' in methods and 'Full LL-PME' in methods,name+': missing paired method '+str(k))
    if 'PRIME' not in methods or 'Full LL-PME' not in methods:continue
    a,b=methods['PRIME'],methods['Full LL-PME'];complete=a['_ok'] and b['_ok']
    vals=pairs.setdefault(k[0],{'attempted':0,'complete':0,'values':[]});vals['attempted']+=1
    if complete:
     av=self.boot[(name,a['_key'])];bv=self.boot[(name,b['_key'])]
     vals['complete']+=1;vals['values'].append(sum(abs(av[j]['estimate']-bv[j]['estimate']) for j in av)/2)
   self.metrics[name]['paired_distances']={k:{'attempted':v['attempted'],'complete':v['complete'],'median':median(v['values']),'mean':avg(v['values'])} for k,v in pairs.items()}
   self.cv(name)
 def cv(self,name):
  rows=self.read(self.role_path(name,'nested_cv.csv'),name+'.cv');cols=(['dataset'] if name=='observed' else [])+['repeat_id','fold','feature_set']
  base=self.read(self.fixed/MODULES[name][0]/'nested_cv.csv','baseline.'+name+'.cv')
  self.demand(len({key(r,cols) for r in rows})==len(rows),name+': duplicated CV key')
  self.demand({key(r,cols) for r in rows}=={key(r,cols) for r in base},name+': missing/extra CV folds or feature sets')
  cohorts=grouped(rows,['dataset']) if name=='observed' else {('engineering',):rows};out={}
  for dataset,rs in cohorts.items():
   med={k[0]:median([number(r['rmse']) for r in v if number(r['rmse']) is not None]) for k,v in grouped(rs,['feature_set']).items()};ratios=[]
   for k,v in grouped(rs,['repeat_id','fold']).items():
    d={r['feature_set']:number(r['rmse']) for r in v};self.demand(all(x is not None and x>0 for x in d.values()),name+': missing/nonpositive held-out RMSE '+str(k))
    if d.get('All coordinates') and d.get('Selected coordinates'):ratios.append(d['Selected coordinates']/d['All coordinates'])
   out[dataset[0]]={'rmse_medians':med,'paired_folds':len(ratios),'median_paired_ratio':median(ratios)}
  self.metrics[name]['cv']=out
 def fixed_endpoints(self):
  selection_path=self.root/'selection_v052_runs.csv'
  selection=self.read(selection_path if selection_path.is_file() else self.fixed/'selection_v052_runs.csv','selection')
  self.demand(len({key(r,('benchmark','n','design','replication','method')) for r in selection})==len(selection),'Duplicate screening records')
  def screening(rs):
   return {'attempted':len(rs),'coverage':avg([int(boolean(r['inclusion'])) for r in rs]),'exact':avg([int(boolean(r['exact'])) for r in rs]),'tpr':avg([float(r['tpr']) for r in rs])}
  self.metrics['screening']={'methods':{k[0]:screening(v) for k,v in grouped(selection,['method']).items()},'prime_by_n':{k[0]:screening(v) for k,v in grouped([r for r in selection if r['method']=='PRIME screen'],['n']).items()}}
  twin_path=self.root/'twinning_v052_main/twinning_runs.csv'
  twin=self.read(twin_path if twin_path.is_file() else self.fixed/'twinning_v052_main/twinning_runs.csv','twinning')
  self.metrics['twinning']={}
  for k,rs in grouped([r for r in twin if r['twin_r']!='all'],['twin_r']).items():
   z={'attempted':len(rs)}
   for c in ('selected_jaccard','active_recall','speedup','query_count'):
    v=[number(r[c]) for r in rs];self.demand(all(x is not None for x in v),'Twinning missing metric '+c)
    v=[x for x in v if x is not None];z[c]={'n_finite':len(v),'median':median(v),'q25':quantile(v,.25),'q75':quantile(v,.75)}
   self.metrics['twinning'][k[0]]=z
  variance=self.read(self.fixed.parent/'provenance/inherited/conditional_variance_summary.csv','fixed.variance_summary')
  ratios=[];cells={}
  for k,rs in grouped(variance,['benchmark','n','rho']).items():
   d={r['method']:number(r['rmse_normalized_tau']) for r in rs}
   self.demand(len(rs)==2 and set(d)=={'Local linear','Local constant'},'Missing/duplicated variance-calibration cell '+str(k))
   if d.get('Local constant') and d.get('Local linear') is not None:
    ratio=d['Local linear']/d['Local constant'];ratios.append(ratio);cells['|'.join(k)]={'ratio':ratio,'reps':[r['reps'] for r in rs]}
  self.metrics['variance_calibration']={'cells':cells,'n_cells':len(cells),'below_diagonal':sum(x<1 for x in ratios),'ratio_median':median(ratios),'ratio_min':min(ratios) if ratios else None,'ratio_max':max(ratios) if ratios else None,'evidence_scope':'inherited summary; not regenerated diagnostic raw observations'}
 def block(self,id,file,start,end,text,deps):
  self.blocks.append({'id':id,'file':file,'start':start,'end':end,'replacement':text.strip()+'\n\n','depends_on':deps})
 def prose(self):
  m=self.metrics;L=m['link']['methods'];p=L['PRIME screen'];f=L['FIRST'];o=L['FOCI'];d=L['FORD']
  def value(z,c='half_l1',stat='mean'):return z.get(c,{}).get(stat)
  def fraction(z):return f"{z['completed']}/{z['attempted']}"
  self.block('main.selector_link','sections/06_experiments.tex','PME based on PRIME-screened coordinates','\\begin{figure}',f"PME based on PRIME-screened coordinates has pooled allocation loss {fmt(value(p))}.\nThe pooled losses over completed runs are {fmt(value(f))} for FIRST,\n{fmt(value(d))} for FORD, and {fmt(value(o))} for FOCI.  Exact allocation\nis completed in {fraction(p)} PRIME-screen runs, with mean inactive-coordinate\nallocation {fmt(value(p,'inactive_mass'))}.  Completion counts are\n{fraction(f)} for FIRST, {fraction(o)} for FOCI, and {fraction(d)} for FORD.",['link'])
  D=m['detail']['methods'];q=D['PRIME screen'];N=q['attempted'];coverage=q['active_inclusion']['count'];exact=q['exact']['count']
  others=['FIRST','FORD','FOCI'];ex=human_list([f"{D[x]['exact']['count']} for {x}" for x in others]);loss=human_list([f"{fmt(value(D[x]))} for {x}" for x in others]);comp=human_list([f"{fraction(D[x])} for {x}" for x in ['PRIME screen','FIRST','FOCI','FORD']])
  self.block('main.friedman_detail','sections/06_experiments.tex','The PRIME screen covers all five active coordinates','\\subsection{Scaling in Ambient Dimension}',f"The PRIME screen covers all five active coordinates in {coverage} of {N}\nreplications and recovers the exact set in {exact}.  Exact-recovery counts\nare {ex}.  Mean PME errors among completed runs\nare {loss}, compared with {fmt(value(q))} for\nthe PRIME-screened calculation.  Completion counts are {comp}.\nThe supplement reports the Branin--Bratley comparison across three sample sizes.",['detail'])
  S=m['scaling'];P=S['methods']['PRIME'];cells=len(S['cells']);wins=S['strict_wins'].get('PRIME',0);ties=S['tied_wins'].get('PRIME',0)
  rank=(f"PRIME has the lowest mean error in all {cells} response--dimension cells." if wins==cells else f"PRIME has the lowest mean error in {wins} of the {cells} response--dimension cells.")
  if ties:rank+=f" It ties for lowest error in {ties} additional cells."
  order=['FORD','FIRST','FOCI','NNVS','LASSO'];loss=human_list([f"{fmt(value(S['methods'][x+' screen + PME-kNN']))} for {x}" for x in order])
  comp=human_list([f"{fraction(S['methods'][x])} for {SHORT[x]}" for x in PIPE])
  calls=S['prime_calls_by_p'];unique={v['median'] for v in calls.values()};calltext=(f"The median number of PRIME subset-value evaluations is {next(iter(unique)):g}\nat every ambient dimension." if len(unique)==1 else 'Median PRIME subset-value counts are '+', '.join(f"{v['median']:g} at p={k}" for k,v in sorted(calls.items(),key=lambda x:int(x[0])))+'.')
  self.block('main.scaling','sections/06_experiments.tex','PRIME has the lowest mean error','These small allocation errors',f"{rank}\nIts pooled mean error over completed runs is {fmt(value(P))}, compared with\n{loss}, each followed by PME-kNN.\nCompletion counts are {comp}.\n{calltext} Final enumeration depends on the selected dimension.",['scaling'])
  eng={k:S['prime_coverage_by_response'][k] for k in ('otl','piston','borehole')};covtext=('In every OTL, piston, and borehole run, the screen omits at least one active\ncoordinate.' if all(x['count']==0 for x in eng.values()) else 'Complete active-set coverage is '+', '.join(f"{x['count']}/{x['denominator']} for {k}" for k,x in eng.items())+'.')
  self.block('main.engineering_coverage','sections/06_experiments.tex','These small allocation errors','\\subsection{Observed-Data Applications}',f"{covtext}\nAllocation accuracy and complete active-set coverage are therefore reported\nas distinct empirical endpoints.",['scaling'])
  C=m['concrete'];E=C['coordinates']['engineering|PRIME'];cv=C['cv']['engineering'];rm=cv['rmse_medians']
  self.block('main.concrete','sections/06_experiments.tex','PRIME and full-coordinate local-linear PME (LL-PME)','\\begin{figure}',f"For concrete, we compare PRIME with full-coordinate local-linear PME (LL-PME).\nAge and the water-to-binder ratio receive median PRIME\nallocations {fmt(E['age']['median'])} and {fmt(E['water_binder_ratio']['median'])}, respectively.\nThe prespecified pair consisting of these two coordinates has median test\nroot mean squared error (RMSE) {fmt(rm['Ratio and age'],2)}, compared with\n{fmt(rm['All coordinates'],2)} for all eight coordinates.",['concrete'])
  datasets=[('Concrete',C,'engineering'),('Airfoil',m['observed'],'airfoil'),('Auto MPG',m['observed'],'auto_mpg')]
  dist=human_list([f"{fmt(z['paired_distances'][k]['median'])} for {label}" for label,z,k in datasets])
  paired='; '.join(f"{label}: {z['paired_distances'][k]['complete']}/{z['paired_distances'][k]['attempted']}" for label,z,k in datasets)
  calls='; '.join(f"{label}: {value(z['groups'][k+'|Full LL-PME'],'value_calls','median'):g} to {value(z['groups'][k+'|PRIME'],'value_calls','median'):g}" for label,z,k in datasets)
  ratios=human_list([f"{fmt(z['cv'][k]['median_paired_ratio'])} for {label}" for label,z,k in datasets])
  complete_pairs=[z['paired_distances'][k] for label,z,k in datasets]
  pairtext=(f"All {complete_pairs[0]['attempted']} bootstrap pairs are complete for each dataset." if all(x['complete']==x['attempted']==complete_pairs[0]['attempted'] for x in complete_pairs) else f"The distances use paired completed allocations ({paired}).")
  reductions=[(label,value(z['groups'][k+'|Full LL-PME'],'value_calls','median'),value(z['groups'][k+'|PRIME'],'value_calls','median')) for label,z,k in datasets]
  calltext=('Screening reduces median subset-value counts '+human_list([f"from {full:g} to {selected:g} for {label}" for label,full,selected in reductions])+'.' if all(selected<full for label,full,selected in reductions) else 'Median full-coordinate and selected-set subset-value counts are '+calls+'.')
  self.block('main.observed_summary','sections/06_experiments.tex','Median screened-to-full allocation distances',None,f"Median screened-to-full allocation distances are {dist}.\n{pairtext}  {calltext}\nThe median paired selected-to-full test-RMSE ratios are {ratios}.\nCoordinate-level bootstrap summaries are reported in the supplement.",['concrete','observed'])
  R=m['correlated'];ranks=(f"PRIME has the lowest mean loss in all {len(R['cells'])} response--sample-size cells." if R['strict_wins'].get('PRIME',0)==len(R['cells']) else f"PRIME has the lowest mean loss in {R['strict_wins'].get('PRIME',0)} of {len(R['cells'])} cells.")
  if R['tied_wins'].get('PRIME',0):ranks+=f" It ties for lowest loss in {R['tied_wins']['PRIME']} additional cells."
  coverage={x:v['active_inclusion']['proportion'] for x,v in R['methods'].items()};proxy={x:value(v,'proxy_mass') for x,v in R['methods'].items()};bestcov=max(coverage.values());bestproxy=min(v for v in proxy.values() if v is not None)
  names=lambda d,v:', '.join(SHORT[x] for x,a in d.items() if a is not None and abs(a-v)<1e-12)
  self.block('supp.correlated','sections/supplement_results.tex','Active-set coverage is the proportion','\\subsection{Ambient dimension}',f"Active-set coverage is the proportion of all runs whose selected set contains\nevery active coordinate.  Proxy allocation is the total PME assigned to the\npaired inactive proxies.  Mean loss and proxy allocation use completed\nallocations, while the table reports completion counts over all attempts.\n{ranks}  The highest pooled active-set coverage is\n{fmt(bestcov)} ({names(coverage,bestcov)}), and the smallest pooled proxy\nallocation is {fmt(bestproxy,4)} ({names(proxy,bestproxy)}).",['correlated'])
  raw=C['coordinates']['raw|PRIME'];namesraw=['cement','water','age','slag'];freq=human_list([f"{fmt(raw[x]['selection_frequency'])} for {x}" for x in namesraw]);vals=human_list([fmt(raw[x]['median']) for x in namesraw])
  freqtext=(f"cement, water, age, and slag are selected in all {raw['cement']['attempted']} bootstrap samples" if all(raw[x]['selection_frequency']==1 and raw[x]['attempted']==raw['cement']['attempted'] for x in namesraw) else 'selection frequencies are '+freq)
  self.block('supp.concrete_raw','sections/supplement_results.tex','In the original representation','\\input{tables/observed_allocation_supplement}',f"In the original representation, {freqtext}.\nTheir median PME values are {vals}, respectively, among\ncompleted bootstrap allocations.",['concrete'])
  self.block('supp.link_lasso_completion','sections/supplement_experiments.tex','The selector-to-PME figure displays four selectors','\\subsection{Representative-query acceleration}',f"The selector-to-PME figure displays four selectors; LASSO is also evaluated\nand completes {L['LASSO']['completed']} of its {L['LASSO']['attempted']} allocations.",['link'])
 def unchanged_prose_checks(self,candidate):
  main=collapse((candidate/'sections/06_experiments.tex').read_text());supp=collapse((candidate/'sections/supplement_results.tex').read_text());m=self.metrics;checks={}
  try:
   n=int(re.search(r'Across the ([0-9,]+) paired datasets',main).group(1).replace(',',''))
   block=main.split('the PRIME screen has active-set coverage',1)[1].split('We next measure',1)[0]
   values=re.findall(r'(?<![0-9])([0-9]+\.[0-9]+)',block)
   p=m['screening']['methods']['PRIME screen'];f=m['screening']['methods']['FIRST'];by_n=m['screening']['prime_by_n']
   expected=[fmt(p[x]) for x in ('coverage','exact','tpr')]+[fmt(f[x]) for x in ('coverage','exact','tpr')]+[fmt(by_n[x]['coverage']) for x in ('100','200')]
   checks['main.recovery']=n==p['attempted']==f['attempted'] and values==expected and all(v['coverage']==1 for k,v in by_n.items() if int(k)>200) and 'and one at every larger sample size' in block
  except (AttributeError,IndexError,KeyError,ValueError):checks['main.recovery']=False
  try:
   block=supp.split('All ',1)[1].split('\\subsection{Coordinate screening}',1)[0];n=int(re.match(r'([0-9]+)',block).group(1));values=re.findall(r'(?<![0-9])([0-9]+\.[0-9]+)',block);v=m['variance_calibration']
   checks['supp.variance_calibration']=n==v['n_cells']==v['below_diagonal'] and values==[fmt(v[x]) for x in ('ratio_median','ratio_min','ratio_max')]
  except (AttributeError,IndexError,KeyError,ValueError):checks['supp.variance_calibration']=False
  try:
   block=supp.split('The median selected-set Jaccard index',1)[1];match=re.search(r'Median runtime speedup increases from ([0-9.]+) at ratio ([0-9]+) to ([0-9.]+) at ratio ([0-9]+)',block);a,ra,b,rb=match.groups();tw=m['twinning'];ordered=[tw[k]['speedup']['median'] for k in sorted(tw,key=int)]
   checks['supp.twinning']=all(x['selected_jaccard']['median']==1 for x in tw.values()) and 'is one for every split ratio' in block and a==fmt(tw[ra]['speedup']['median'],2) and b==fmt(tw[rb]['speedup']['median'],2) and all(x<y for x,y in zip(ordered,ordered[1:]))
  except (AttributeError,IndexError,KeyError,ValueError):checks['supp.twinning']=False
  return checks
 def write(self):
  out=Path(self.a.out).resolve();candidate=Path(self.a.candidate).resolve()
  if out==candidate or candidate in out.parents:raise ValueError('--out must not be inside candidate')
  out.mkdir(parents=True,exist_ok=True)
  paths={b['file'] for b in self.blocks};texts={f:(candidate/f).read_text() for f in paths};checks={}
  for b in self.blocks:checks[b['id']]=collapse(b['replacement']) in collapse(texts[b['file']])
  checks.update(self.unchanged_prose_checks(candidate))
  finalized=self.a.finalized_inputs and self.a.observed_current and self.a.timings_current and not self.a.historical_summary
  report={'stage':'finalized by caller' if finalized else 'provisional; not release validation','attestations':{'inputs_finalized':self.a.finalized_inputs,'observed_neighbor_stage_current':self.a.observed_current,'reported_timing_columns_current':self.a.timings_current},'errors':self.errors,'warnings':self.warnings,'sources':self.sources,'metrics':self.metrics,'prose_contract_matches':checks,'release_validation_pass':bool(finalized and not self.errors and all(checks.values())),'limits':['Caller attestations identify final method/timing stages; arithmetic does not establish source/estimator parity.','Exact generated-block matching validates the prose contract, not arbitrary alternative wording.','Qualitative global conclusions and nonnumeric method claims still require review.']}
  (out/'numeric_claim_snapshot.json').write_text(json.dumps(report,indent=2,allow_nan=False))
  (out/'proposed_blocks.json').write_text(json.dumps(self.blocks,indent=2))
  if self.a.compare_snapshot:
   old=json.loads(Path(self.a.compare_snapshot).read_text());changes=[]
   def walk(a,b,path=''):
    if isinstance(a,dict) and isinstance(b,dict):
     for k in sorted(set(a)|set(b)):walk(a.get(k),b.get(k),path+'.'+k if path else k)
    elif a!=b:changes.append({'claim':path,'before':a,'after':b})
   walk(old['metrics'],self.metrics)
   (out/'numeric_impact.json').write_text(json.dumps(changes,indent=2))
   with (out/'numeric_impact.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['claim','before','after','delta']);w.writeheader()
    for change in changes:
     a,b=change['before'],change['after'];delta=b-a if isinstance(a,(int,float)) and isinstance(b,(int,float)) else ''
     w.writerow(dict(change,delta=delta))
  if self.a.emit_proposal and not self.errors:
   for b in self.blocks:
    t=texts[b['file']];start=b['start']
    beginmark='% numeric-claim '+b['id']+' begin';endmark='% numeric-claim '+b['id']+' end'
    if beginmark in t:
     self.demand(t.count(beginmark)==t.count(endmark)==1,'Ambiguous numeric markers '+b['id'])
     if t.count(beginmark)!=1 or t.count(endmark)!=1:continue
     a=t.index(beginmark);end=t.index(endmark,a)+len(endmark)
     replacement=beginmark+'\n'+b['replacement'].rstrip()+'\n'+endmark
    else:
     self.demand(t.count(start)==1,'Ambiguous/missing prose update anchor '+b['id'])
     if t.count(start)!=1:continue
     a=t.index(start);end=t.find(b['end'],a+len(start)) if b['end'] else len(t)
     replacement=beginmark+'\n'+b['replacement'].rstrip()+'\n'+endmark+'\n\n'
    self.demand(end>=a,'Missing end anchor '+b['id'])
    if end>=a:texts[b['file']]=t[:a]+replacement+t[end:]
   if not self.errors:
    for f,t in texts.items():
     p=out/'proposed'/f;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(t)
     old=(candidate/f).read_text();(out/(Path(f).name+'.diff')).write_text(''.join(difflib.unified_diff(old.splitlines(True),t.splitlines(True),fromfile='candidate/'+f,tofile='proposed/'+f)))
  # Re-save any strict anchor failures rather than reporting a partial success.
  report['errors']=self.errors;report['release_validation_pass']=bool(finalized and not self.errors and all(checks.values()))
  (out/'numeric_claim_snapshot.json').write_text(json.dumps(report,indent=2,allow_nan=False))
  print(json.dumps({'errors':len(self.errors),'warnings':len(self.warnings),'blocks':len(self.blocks),'matching_current_prose_blocks':sum(checks.values()),'release_validation_pass':report['release_validation_pass'],'output':str(out)},indent=2))
  return 1 if self.errors or (self.a.validate_prose and not report['release_validation_pass']) else 0

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--data-root',required=True);p.add_argument('--baseline-root',required=True);p.add_argument('--fixed-root');p.add_argument('--candidate',required=True);p.add_argument('--out',required=True)
 p.add_argument('--compare-snapshot');p.add_argument('--emit-proposal',action='store_true');p.add_argument('--validate-prose',action='store_true')
 p.add_argument('--historical-summary',action='store_true',help='Explicitly retain 19 known missing-vector LASSO historical run summaries for baseline comparison; always disables release validation')
 p.add_argument('--finalized-inputs',action='store_true');p.add_argument('--observed-current',action='store_true');p.add_argument('--timings-current',action='store_true')
 a=p.parse_args();audit=Audit(a)
 try:
  for name in MODULES:audit.load_module(name)
  if audit.errors:
   Path(a.out).mkdir(parents=True,exist_ok=True);(Path(a.out)/'numeric_input_errors.json').write_text(json.dumps(audit.errors,indent=2));print(json.dumps({'errors':audit.errors[:20],'total':len(audit.errors)}));return 1
  audit.allocation_details();audit.fixed_endpoints();audit.prose();return audit.write()
 except (OSError,ValueError,KeyError,TypeError) as e:
  print('Numeric claim audit failed: '+str(e),file=sys.stderr);return 2
if __name__=='__main__':sys.exit(main())
