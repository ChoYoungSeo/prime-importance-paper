#!/usr/bin/env python3
"""Validate the dated numerical overlay using only bundled CSVs (no fitting)."""
import argparse, csv, hashlib, json, math, re, statistics
from collections import defaultdict, Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'results/numerical_20261001'
def rows(name):
    with (DATA/name).open(newline='') as f:return list(csv.DictReader(f))
def close(a,b,tol=1e-10):
    assert math.isfinite(float(a)) and math.isfinite(float(b))
    assert abs(float(a)-float(b))<=tol*max(1,abs(float(a)),abs(float(b))), (a,b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def groups(rr,keys):
    out=defaultdict(list)
    for r in rr:out[tuple(r[k] for k in keys)].append(r)
    return out

def validate(historical_data_only=False):
    manifest=json.loads((DATA/'MANIFEST.json').read_text())
    scope=manifest['public_export_scope']
    assert all(scope[k]=='NOT_IN_PUBLIC_EXPORT' for k in ('manuscript_sources','paper_prose','typeset_table'))
    delegated=[]
    current=json.loads((ROOT/'results/numerical_20261011/manifest.json').read_text()) if historical_data_only else None
    for path,value in manifest['files'].items():
        presentation=(path.startswith('manuscript/') or path.startswith('docs/assets/') or path=='docs/numerical-evidence.md')
        maintenance=path=='experiments/14_current_numerical/validate.py'
        if historical_data_only and (presentation or maintenance):
            # Updated presentation and this validator are pinned by the new overlay.
            # Every numerical input and original scientific reproduction script
            # still has to match the untouched historical manifest below.
            assert path in current['files'] and sha(ROOT/path)==current['files'][path],path
            delegated.append(path)
        else:assert sha(ROOT/path)==value, path
    assert sha(ROOT/'results/current_release.json')==manifest['historical_release_sha256']
    gamma=23.82646377138105/(23.82646377138105+1)
    phi=rows('friedman/native_phi.csv'); runs=rows('friedman/runs.csv')
    assert len(phi)==8100 and len(runs)==810
    assert all(r['status']=='complete' for r in phi+runs)
    assert set(Counter((r['n'],r['rho'],r['method']) for r in runs).values())=={30}
    assert len(Counter((r['n'],r['rho'],r['method']) for r in runs))==27
    runmap={(r['task_id'],r['method']):r for r in runs}
    inp={r['task_id']:r for r in rows('friedman/input_manifest.csv')}
    assert len(inp)==270 and len(runmap)==810
    losses=defaultdict(list)
    for key,rr in groups(phi,['task_id','method']).items():
        r=runmap[key];i=inp[key[0]]
        assert len(rr)==10 and {int(x['variable']) for x in rr}==set(range(1,11))
        for x in rr:
            assert math.isfinite(float(x['estimate']))
            assert x['standardized_X_hash']==r['standardized_X_hash']==i['standardized_X_hash']
            assert x['input_hash']==r['input_hash']==i['input_hash']
        close(sum(float(x['estimate']) for x in rr),r['effect_sum'])
        if key[1]=='FIRST-SHAP':
            assert all(x['truth_native']=='NA' for x in rr)
            assert r['target_relative_half_l1']=='NA'
            continue
        scale=gamma if key[1]=='PME-kNN' else 1
        for x in rr:close(x['truth_native'],scale*float(x['truth_unit_signal']))
        e=sum(abs(float(x['estimate'])-float(x['truth_native'])) for x in rr)/(2*scale)
        close(e,r['target_relative_half_l1'])
        losses[(int(r['n']),float(r['rho']),key[1])].append(e)
    llruns={r['task_id']:r for r in rows('friedman/full_ll_runs.csv')}
    ll=[r for r in rows('friedman/full_ll_allocations.csv') if r['rule']=='PME']
    assert len(ll)==2700 and len(llruns)==270
    for key,rr in groups(ll,['task_id']).items():
        r=llruns[key[0]];assert r['input_md5']==inp[key[0]]['input_file_md5']
        assert len(rr)==10 and {int(x['variable_index']) for x in rr}==set(range(1,11))
        truth={int(x['variable']):float(x['truth_unit_signal']) for x in phi if x['task_id']==key[0] and x['method']=='PRIME'}
        for x in rr:close(x['population_target'],truth[int(x['variable_index'])])
        e=sum(abs(float(x['estimate'])-float(x['population_target'])) for x in rr)/2
        close(e,r['pme_half_l1']);losses[(int(r['n']),float(r['rho']),'Full LL')].append(e)
    table=[]
    for n in (100,200,500):
        for rho in (0,.5,.9):
            vals={m:statistics.mean(losses[n,rho,m]) for m in ('PRIME','Full LL','PME-kNN')}
            assert all(len(losses[n,rho,m])==30 for m in vals)
            assert vals['PRIME']<min(vals['Full LL'],vals['PME-kNN'])
            table.append({'n':n,'rho':rho,**vals})
    # Public scope: validate the saved numeric reference, not withheld TeX/typesetting.
    golden=rows('friedman/earlier_four_method_table.csv')
    assert len(golden)==9
    for raw,r in zip(golden,table):
        assert [float(raw['n']),float(raw['rho'])]==[r['n'],r['rho']]
        for method in ('PRIME','Full LL','PME-kNN'):
            close(raw[method],r[method])
    concrete=rows('concrete/paired_phi.csv'); refs=rows('concrete/full_data_reference.csv')
    assert len(concrete)==2400 and len(refs)==24
    cstatus=rows('concrete/all_status.csv'); assert len(cstatus)==303 and all(r['status']=='complete' for r in cstatus)
    cr={r['id']:r for r in rows('concrete/standardization_receipts.csv')}; assert len(cr)==101
    for key,rr in groups(concrete,['replication','method']).items():
        assert len(rr)==8 and {int(x['variable']) for x in rr}==set(range(1,9))
        assert all(x['n']=='824' and x['status']=='complete' and math.isfinite(float(x['estimate'])) for x in rr)
        assert all(x['standardized_hash']==cr[x['id']]['standardized_hash'] for x in rr)
    idx=groups(rows('concrete/sample_indices.csv'),['replication'])
    for key,rr in idx.items():
        if key==('0',):continue
        assert len(rr)==824 and len({x['source_row'] for x in rr})==824
        assert all(1<=int(x['source_row'])<=1030 for x in rr)
    figure_stats={}
    for num,source,p,nrep in [(4,[x for x in phi if x['n']=='100' and x['rho']=='0.9'],10,30),(5,concrete,8,100)]:
        a={(int(x['replication']),x['method'],int(x['variable'])):float(x['estimate']) for x in source}
        b={}
        for x in rows(f'figures/figure{num}_native.csv'):
            if x['method']=='PME-kNN (noise = FALSE)':continue
            m=x['method'].replace('PME-kNN (noise = TRUE)','PME-kNN')
            b[int(x['replication']),m,int(x['variable'])]=float(x['estimate'])
        assert a.keys()==b.keys() and len(a)==3*p*nrep
        for k in a:close(a[k],b[k],1e-13)
        outside=sum(v<-.30 or v>1.05 for v in b.values()); assert outside==({4:13,5:29}[num])
        figure_stats[str(num)]={'coordinates':len(b),'boxes':3*p,'replications':nrep,'outside_detail':outside}
    tt=rows('timing/timings.csv'); assert len(tt)==120
    assert all(x['status']=='complete' and x['exact_full_state']=='TRUE' and x['standardization_exact']=='TRUE' and float(x['effect_max_abs_gap'])==0 for x in tt)
    med={}
    for key,rr in groups(tt,['method']).items():
        assert len(rr)==30
        assert {x['task_id'] for x in rr}=={x['task_id'] for x in phi if x['n']=='100' and x['rho']=='0.9'}
        for x in rr:assert x['input_sha256']==inp[x['task_id']]['input_file_sha256'] and x['standardized_X_hash']==inp[x['task_id']]['standardized_X_hash']
        med[key[0]]=statistics.median(float(x['elapsed_sec']) for x in rr)
    for r in rows('timing/summary.csv'):close(med[r['method']],r['median_sec'])
    assert [f'{med[m]:.3f}' for m in ('PRIME','PME_kNN_TRUE','FIRST_SHAP')]==['0.019','1.706','0.029']
    assert sum(float(x['estimate'])<0 for x in concrete if x['method']=='PME-kNN')==314
    # Supplementary K6 FALSE/TRUE native summaries use the same saved figure vectors.
    supplemental={}
    for num,pcount in [(4,10),(5,8)]:
        native=rows(f'figures/figure{num}_native.csv')
        by=groups(native,['method'])
        for noise in ('FALSE','TRUE'):
            m=f'PME-kNN (noise = {noise})';vv=[float(x['estimate']) for x in by[(m,)]]
            expected={(4,'FALSE'):(0,.030,.319,0),(4,'TRUE'):(116,-2.618,3.771,13),
                      (5,'FALSE'):(0,.021,.575,0),(5,'TRUE'):(314,-3.583,3.080,29)}[num,noise]
            actual=(sum(v<0 for v in vv),float(f'{min(vv):.3f}'),float(f'{max(vv):.3f}'),sum(v<-.3 or v>1.05 for v in vv))
            assert actual==expected,(num,noise,actual)
            supplemental[f'figure{num}_{noise}']={'negative':actual[0],'minimum':actual[1],'maximum':actual[2],'outside':actual[3]}
    native=rows('figures/figure5_native.csv'); ref=rows('figures/figure5_references.csv')
    def quantile(v,q):
        v=sorted(v);t=(len(v)-1)*q;i=int(t);return v[i]+(t-i)*(v[min(i+1,len(v)-1)]-v[i])
    iqr={}
    for (m,),rr in groups(native,['method']).items():
        byvar=groups(rr,['variable']);assert len(byvar)==8 and all(len(v)==100 for v in byvar.values())
        iqr[m]=statistics.mean(quantile([float(x['estimate']) for x in v],.75)-quantile([float(x['estimate']) for x in v],.25) for v in byvar.values())
    assert [f'{iqr[m]:.3f}' for m in ['PRIME','PME-kNN (noise = FALSE)','PME-kNN (noise = TRUE)','FIRST-SHAP']]==['0.017','0.036','0.096','0.021']
    water=[float(x['estimate']) for x in native if x['method']=='PRIME' and x['variable']=='4']
    assert [f'{quantile(water,q):.3f}' for q in [.25,.5,.75]]==['0.438','0.455','0.478']
    refmap={(x['method'],x['variable']):float(x['estimate']) for x in ref}
    assert f"{refmap['PRIME','4']:.3f}"=='0.385' and f"{refmap['PRIME','8']:.3f}"=='0.513'
    true=[x for x in native if x['method']=='PME-kNN (noise = TRUE)']
    vectors=groups(true,['replication']);assert len(vectors)==100 and all(any(float(x['estimate'])<0 for x in v) for v in vectors.values())
    totals=[sum(float(x['estimate']) for x in v) for v in vectors.values()]
    assert f'{statistics.median(totals):.3f}'=='0.708'
    assert f"{sum(v for (m,j),v in refmap.items() if m=='PME-kNN (noise = TRUE)'):.3f}"=='0.729'
    return {'status':'PASS','public_document_scope':'NOT_IN_PUBLIC_EXPORT','numeric_accuracy_reference_settings':9,'numeric_accuracy_reference_methods':3,'scope':'Bundled saved-output numerical validation; no estimator reruns','historical_data_only':historical_data_only,'historical_hashes_verified':len(manifest['files'])-len(delegated),'current_overlay_hashes_verified':len(delegated),'presentation_and_validator_delegated_to_current_overlay':delegated,'frozen_files_verified':len(manifest['files']),'friedman_runs':810,'friedman_coordinates':8100,'full_LL_coordinates':2700,'table_cells':9,'table_values':table,'concrete_runs':303,'concrete_coordinates':2400,'figures':figure_stats,'serial_timed_calls':120,'serial_median_seconds':med,'gamma':gamma,'legacy_release_unchanged':True,'supplement_native_summaries':supplemental}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical-data-only',action='store_true',help='Keep historical data/scientific-code hashes; require updated presentation and validator hashes in the 2026-10-11 overlay.')
    print(json.dumps(validate(parser.parse_args().historical_data_only),indent=2))
