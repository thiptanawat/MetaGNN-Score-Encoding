"""Post-evaluation check using origin adjacency matrices, without evaluator imports.

This validation script is not a new analysis arm. It reconstructs the same frozen
estimand from retained observations and predictions after the official evaluation.
"""
import argparse,csv,hashlib,itertools,json
from collections import defaultdict
from pathlib import Path
import numpy as np


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--results-dir',type=Path,required=True);ap.add_argument('--scenario',default='primary');a=ap.parse_args()
    root=a.root;run=a.results_dir
    official_path=run/f'evaluation_{a.scenario}.json'
    official=json.loads(official_path.read_text());ref=official['primary_fixed_panel']
    assert official['partition']=='test' and ref['status']=='complete_predictions'
    lock=json.loads((root/'protocol/PROTOCOL_LOCK.json').read_text())
    for group in ['input_sha256','source_sha256']:
        for path,h in lock[group].items():assert sha(root/path)==h,path
    plan=json.loads((root/'protocol/analysis_plan.json').read_text())
    context={r['context_id']:r['origin_group'] for r in read(root/'manifests/contexts.tsv') if r['partition']=='test'}
    origins=sorted(set(context.values()));oi={s:i for i,s in enumerate(origins)}
    targets=plan['primary_targets'];mids=[t['metabolite_id'] for t in targets]
    arms=['magnitude_g0.50','magnitude_g0.80','magnitude_g1.00','magnitude_g1.25','magnitude_g2.00','ordinal','uniform_pfba']
    raw=defaultdict(list)
    for r in read(root/'data/outcomes_long.tsv'):
        if r['context_id'] in context and r['metabolite_id'] in mids:raw[(r['context_id'],r['metabolite_id'])].append(float(r['value']))
    means={k:sum(v)/len(v) for k,v in raw.items()}
    pred={(r['context_id'],r['arm'],r['exchange_id']):r for r in read(run/'predictions.tsv') if r['scenario']==a.scenario}
    count=np.zeros((len(targets),len(origins),len(origins)))
    value=np.zeros((len(arms),len(targets),len(origins),len(origins)))
    for m,t in enumerate(targets):
        mid=t['metabolite_id'];ex=t['exchange_id'];factor=float(t['sign_factor']);tol=plan['assay_tolerances'][mid]
        for left,right in itertools.combinations(sorted(context),2):
            if context[left]==context[right]:continue
            obs=means[(left,mid)]-means[(right,mid)]
            if abs(obs)<=tol:continue
            i,j=oi[context[left]],oi[context[right]];count[m,i,j]+=1
            for armidx,arm in enumerate(arms):
                lp=pred[(('__reference__' if arm=='uniform_pfba' else left),arm,ex)]
                rp=pred[(('__reference__' if arm=='uniform_pfba' else right),arm,ex)]
                assert lp['status']==rp['status']=='optimal'
                delta=factor*(float(lp['point'])-float(rp['point']))
                assert np.isfinite(delta)
                value[armidx,m,i,j]+=.5 if abs(delta)<=plan['point_tolerance'] else float(delta*obs>0)
    denom=count.sum(axis=(1,2));assert (denom>0).all()
    c=(value.sum(axis=(2,3))/denom).mean(axis=1)
    contrasts=np.array([c[5]-c[:5].mean(),c[5]-c[2]])
    rng=np.random.default_rng(plan['test_seed'])
    mult=np.array([np.bincount(rng.choice(len(origins),len(origins),replace=True),minlength=len(origins)) for _ in range(plan['bootstrap_resamples'])])
    denominator=np.einsum('bi,mij,bj->bm',mult,count,mult,optimize=True)
    valid=(denominator>0).all(axis=1);denominator=denominator[valid];vm=mult[valid]
    samples=np.stack([(np.einsum('bi,mij,bj->bm',vm,v,vm,optimize=True)/denominator).mean(axis=1) for v in value],axis=1)
    contrast_samples=np.stack((samples[:,5]-samples[:,:5].mean(axis=1),samples[:,5]-samples[:,2]),axis=1)
    ci=np.percentile(contrast_samples,[2.5,97.5],axis=0).T
    abs_ci=np.percentile(samples,[2.5,97.5],axis=0).T
    expectedc=np.array([next(r['macro_point_operational'] for r in ref['arms'] if r['arm']==arm) for arm in arms])
    expci=np.array([ref['bootstrap']['delta_grid_interval'],ref['bootstrap']['delta_identity_interval']])
    expabs=np.array([ref['bootstrap']['absolute_C_intervals'][arm] for arm in arms])
    errors={'absolute_C':float(np.max(abs(c-expectedc))),'contrasts':float(np.max(abs(contrasts-np.array([ref['primary']['delta_grid'],ref['primary']['delta_identity']])))),
      'contrast_intervals':float(np.max(abs(ci-expci))),'absolute_intervals':float(np.max(abs(abs_ci-expabs)))}
    assert int(denom.sum())==ref['observed_eligible_pairs']
    assert int((~valid).sum())==ref['bootstrap']['invalid_draws']
    assert max(errors.values())<1e-12,errors
    result={'status':'PASS','scenario':a.scenario,'method':'Independent pair construction and origin adjacency-matrix quadratic forms; no imports from production evaluator.',
      'observed_pairs':int(denom.sum()),'targets':len(targets),'origins':len(origins),'valid_draws':int(valid.sum()),'invalid_draws':int((~valid).sum()),
      'absolute_C':dict(zip(arms,c.tolist())),'delta_grid':float(contrasts[0]),'delta_identity':float(contrasts[1]),
      'delta_grid_interval':ci[0].tolist(),'delta_identity_interval':ci[1].tolist(),'maximum_absolute_differences':errors,
      'official_evaluation_sha256':sha(official_path),'prediction_sha256':sha(run/'predictions.tsv'),'validation_script_sha256':sha(Path(__file__)),
      'scope':'Computational cross-check of the same frozen estimand and paired draws, not independent biological replication.'}
    out=run/f'independent_concordance_{a.scenario}.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
