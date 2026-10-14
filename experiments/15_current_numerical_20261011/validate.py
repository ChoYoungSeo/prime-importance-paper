#!/usr/bin/env python3
"""Validate the portable four-method overlay using the standard library only."""
import csv
from collections import Counter,defaultdict
from decimal import Decimal,ROUND_HALF_UP
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'results/numerical_20261011'
OLD=ROOT/'results/numerical_20261001'
METHODS=('PRIME','PME-NA','FIRST-SHAP','SHAP')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def close(a,b,tolerance=1e-12):
    a,b=float(a),float(b)
    assert math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tolerance*max(1,abs(a),abs(b)),(a,b)
def group(data,keys):
    result=defaultdict(list)
    for row in data:result[tuple(row[k] for k in keys)].append(row)
    return result
def input_path(relative):
    path=(DATA/relative).resolve()
    assert DATA.resolve() in path.parents and not Path(relative).is_absolute(),relative
    return path
def label(dataset,value):
    if dataset=='friedman':return f'{value:.3f} s'
    return str(Decimal(str(value)).quantize(Decimal('.000000001'),rounding=ROUND_HALF_UP).quantize(Decimal('.001'),rounding=ROUND_HALF_UP))+' s'

def validate():
    manifest=read(DATA/'manifest.json')
    scope=manifest['public_export_scope']
    assert not (ROOT/'manuscript').exists(), 'Manuscript directory is outside the public export scope'
    assert all(scope[k]=='NOT_IN_PUBLIC_EXPORT' for k in ('manuscript_sources','paper_prose','typeset_table'))
    assert manifest['schema_version']==1
    for relative,digest in manifest['files'].items():assert sha(ROOT/relative)==digest,relative
    assert sha(ROOT/'results/current_release.json')==manifest['historical_release_sha256']
    assert sha(OLD/'MANIFEST.json')==manifest['historical_overlay_manifest_sha256']
    protocol=read(DATA/'protocol.json');cases=read(DATA/'cases.json');saved=read(DATA/'results/shap_cases.json')
    expected=[f'friedman_n100_rho0.9_r{r:03d}' for r in range(1,31)]+[f'concrete_n0824_r{r:02d}' for r in range(1,101)]
    assert protocol['case_ids']==expected==[c['case_id'] for c in cases]==[c['case_id'] for c in saved]
    case_map={c['case_id']:c for c in cases};saved_map={c['case_id']:c for c in saved}
    assert len(case_map)==len(saved_map)==130
    assert protocol['analysis']=='full_sample_in_sample'
    assert protocol['shap']==dict(data=None,feature_perturbation='tree_path_dependent',model_output='raw',approximate=False,aggregation='np.abs(values).mean(axis=0)',normalization='none')
    assert protocol['candidates']==[dict(n_estimators=n,max_depth=d) for n in (200,500,1000) for d in (2,4)]
    assert protocol['fixed_parameters']==dict(objective='reg:squarederror',learning_rate=.05,reg_lambda=1,min_child_weight=1,subsample=1,colsample_bytree=1,colsample_bylevel=1,colsample_bynode=1,tree_method='hist',n_jobs=1,verbosity=0)
    old_f={r['task_id']:r for r in rows(OLD/'friedman/input_manifest.csv')}
    old_c={r['id']:r for r in rows(OLD/'concrete/paired_phi.csv')}
    old_indices=group(rows(OLD/'concrete/sample_indices.csv'),['replication'])
    bundle_hashes={}
    for c in cases:
        cid=c['case_id'];z=saved_map[cid];n,p=(100,10) if c['dataset']=='friedman' else (824,8)
        assert (c['n'],c['p'])==(n,p)==(z['n'],z['p'])
        raw=gzip.decompress(input_path(c['input_path']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==c['csv_sha256']==z['input_sha256']
        parsed=list(csv.reader(io.StringIO(raw.decode('utf-8'))))
        assert parsed[0]==[f'X{j}' for j in range(1,p+1)]+['Y'] and len(parsed)==n+1
        assert all(len(row)==p+1 and all(math.isfinite(float(v)) for v in row) for row in parsed[1:])
        bundle=c['source_bundle_path']
        assert not Path(bundle).is_absolute()
        if bundle not in bundle_hashes:bundle_hashes[bundle]=sha(ROOT/bundle)
        assert bundle_hashes[bundle]==c['source_bundle_sha256']
        if c['dataset']=='friedman':assert c['source_input_hash']==old_f[cid]['input_hash']
        else:
            assert c['source_input_hash']==old_c[c['source_case_id']]['input_hash']
            indices=c['source_row_indices'];assert len(indices)==len(set(indices))==824 and all(1<=i<=1030 for i in indices)
            assert indices==[int(r['source_row']) for r in old_indices[(str(c['replication']),)]]
        seed=(20261008+int(hashlib.sha256(cid.encode()).hexdigest()[:7],16))%2**32
        assert c['case_seed']==z['case_seed']==seed
        assert c['inner_fold_seed']==z['inner_fold_seed']==(seed+50021)%2**32
        assert c['inner_fit_seeds']==z['inner_fit_seeds']==[(seed+1000003*(j+1)+7919)%2**32 for j in range(3)]
        assert c['fit_seed']==z['fit_seed']==(seed+((1<<p)-1)*7919)%2**32
        assert z['in_sample'] and z['same_training_and_evaluation_rows'] and z['n_train']==z['n_evaluation']==n
        assert z['training_scope']=='full_sample' and z['evaluation_scope']=='same_full_sample_in_sample'
        assert z['feature_perturbation']=='tree_path_dependent' and z['model_output']=='raw' and not z['approximate'] and z['n_background']==0 and z['normalization']=='none'
        assert c['feature_names']==z['feature_names'] and len(z['phi'])==p and all(math.isfinite(v) and v>=0 for v in z['phi'])
        assert all(z['scientific_validation'].values())
        folds=z['inner_validation_indices'];assert len(folds)==3 and sorted(j for fold in folds for j in fold)==list(range(n))
        scores=z['candidate_scores'];assert [s['params'] for s in scores]==protocol['candidates']
        for s in scores:
            assert len(s['inner_fold_r2'])==3 and all(math.isfinite(v) for v in s['inner_fold_r2'])
            close(s['mean_inner_r2'],statistics.mean(s['inner_fold_r2']))
            close(s['se_inner_r2'],statistics.stdev(s['inner_fold_r2'])/math.sqrt(3))
        key=lambda s:(s['params']['n_estimators'],s['params']['max_depth'])
        best=min(scores,key=lambda s:(-s['mean_inner_r2'],*key(s)))
        cutoff=best['mean_inner_r2']-best['se_inner_r2'];close(cutoff,z['one_se_cutoff'])
        chosen=min((s for s in scores if s['mean_inner_r2']>=cutoff-1e-14),key=key)
        assert chosen['params']==z['selected_params']
    values=rows(DATA/'results/plot_values.csv')
    assert len(values)==4400
    vmap={(r['case_id'],r['method'],int(r['variable'])):float(r['estimate']) for r in values};assert len(vmap)==4400
    for c in cases:
        for method in METHODS:
            assert all((c['case_id'],method,j) in vmap for j in range(1,c['p']+1))
        assert all(vmap[c['case_id'],'SHAP',j]==saved_map[c['case_id']]['phi'][j-1] for j in range(1,c['p']+1))
    for dataset,num in [('friedman',4),('concrete',5)]:
        for r in rows(OLD/f'figures/figure{num}_native.csv'):
            method={'PRIME':'PRIME','PME-kNN (noise = TRUE)':'PME-NA','FIRST-SHAP':'FIRST-SHAP'}.get(r['method'])
            if not method:continue
            rep=int(r['replication']);cid=f'friedman_n100_rho0.9_r{rep:03d}' if dataset=='friedman' else f'concrete_n0824_r{rep:02d}'
            assert vmap[cid,method,int(r['variable'])]==float(r['estimate'])
    refs=rows(DATA/'results/friedman_references.csv')
    assert len(refs)==20 and {r['method'] for r in refs}==set(METHODS[:2])
    ref_map={(r['method'],r['variable']):float(r['estimate']) for r in refs}
    for r in rows(OLD/'figures/figure4_references.csv'):
        method={'PRIME':'PRIME','PME-kNN (noise = TRUE)':'PME-NA','FIRST-SHAP':'FIRST-SHAP'}.get(r['method'])
        if method:assert ref_map[method,r['variable']]==float(r['estimate'])
    timings=rows(DATA/'results/timings.csv');assert len(timings)==520
    tmap={(r['case_id'],r['method']):r for r in timings};assert len(tmap)==520
    assert all(r['status']=='PASS' and r['parity_pass']=='True' and math.isfinite(float(r['elapsed_sec'])) and float(r['elapsed_sec'])>0 for r in timings)
    timing_summary=rows(DATA/'results/timing_summary.csv');assert len(timing_summary)==8
    labels={}
    for (dataset,method),rr in group(timings,['dataset','method']).items():
        count=30 if dataset=='friedman' else 100;assert len(rr)==count
        assert {int(r['replication']) for r in rr}==set(range(1,count+1))
        summary=next(s for s in timing_summary if (s['dataset'],s['method'])==(dataset,method))
        median=statistics.median(float(r['elapsed_sec']) for r in rr)
        close(median,summary['median_sec']);assert summary['display']==label(dataset,median)
        labels[dataset,method]=summary['display']
    assert [labels['friedman',m] for m in METHODS]==['0.019 s','1.706 s','0.029 s','0.514 s']
    assert [labels['concrete',m] for m in METHODS]==['0.450 s','13.730 s','0.218 s','2.081 s']
    schedule=rows(DATA/'results/concrete_schedule.csv');assert len(schedule)==400
    assert Counter((r['method'],r['position']) for r in schedule)==Counter({(m,str(p)):25 for m in METHODS for p in range(1,5)})
    patterns=[('PRIME','PME-NA','SHAP','FIRST-SHAP'),('PME-NA','FIRST-SHAP','PRIME','SHAP'),('FIRST-SHAP','SHAP','PME-NA','PRIME'),('SHAP','PRIME','FIRST-SHAP','PME-NA')]
    for r in schedule:
        assert r['method']==patterns[(int(r['replication'])-1)%4][int(r['position'])-1]
        assert tmap[r['case_id'],r['method']]['position']==r['position']
    parity=read(DATA/'results/concrete_timing_checks.json');assert len(parity)==400
    for r in parity:
        assert r['status']=='PASS' and r['parity_pass'] and all(r['checks'].values())
        assert r['input_hash']==case_map[r['case_id']]['source_input_hash']
        close(r['elapsed_sec'],tmap[r['case_id'],r['method']]['elapsed_sec'])
        assert r['source_receipt_sha256']==tmap[r['case_id'],r['method']]['source_receipt_sha256']
        phi=r['phi'] if r['method']=='SHAP' else r['effects']
        for j,value in enumerate(phi,1):close(value,vmap[r['case_id'],r['method'],j],1e-13)
        if r['method']=='SHAP':assert r['selected_params']==saved_map[r['case_id']]['selected_params'] and r['max_mean_abs_shap_error']==0
        else:assert r['max_effect_error']==0
    approved=read(DATA/'provenance/approved_figure_statistics.json')
    for dataset in ('friedman','concrete'):
        stem='PRIME_Friedman' if dataset=='friedman' else 'PRIME_Concrete'
        # Withheld manuscript duplicate is not exported; approved research graphics remain hash-pinned above.
        assert sha(ROOT/f'docs/assets/{dataset}-current.png')==sha(DATA/f'approved_figures/{stem}.png')
        a=approved[dataset]['alignment'];close(a['zero_coordinate_both_panel_npc'],1/22)
        assert a['bottom_zero_whitespace_restored'] and a['first_SHAP_y_range']==[-.05,1.05]
        for method in METHODS[:3]:
            old_method='PME-kNN (noise = TRUE)' if method=='PME-NA' else method
            assert approved[dataset]['baseline_preservation'][old_method]['timing_label']==labels[dataset,method]
    return dict(status='PASS',scope='Saved four-method outputs, portable inputs and retained research graphics; no fitting',
                public_document_scope='NOT_IN_PUBLIC_EXPORT',
                cases=130,plot_coordinates=4400,shap_aggregates=1100,measured_timing_rows=520,
                concrete_paired_timed_calls=400,all_four_methods_positions_balanced=True,
                full_sample_tree_path_shap=True,normalized_shap=False,small_negative_native_estimates_preserved=True,
                historical_release_unchanged=True,labels={f'{d}/{m}':labels[d,m] for d,m in labels},
                manifest_files=len(manifest['files']))

if __name__=='__main__':print(json.dumps(validate(),indent=2))
