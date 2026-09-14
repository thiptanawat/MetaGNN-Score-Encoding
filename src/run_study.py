"""CPU-only optimization runner. Reads score/context metadata, never CORE outcomes."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import csv
import json
import os
import time
import traceback
import hashlib
import numpy as np
from src.cpu_model import load_config, build_model, score_weights, digest
from src.cpu_readout import select_flux

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with open(path) as f: return list(csv.DictReader(f, delimiter='\t'))


def input_fingerprint(root, config):
    paths=['configs/study.json','data/raw/Recon3D.json','data/scores_conservative.npz',
        'manifests/contexts.tsv','manifests/metabolite_map.tsv','src/cpu_model.py',
        'src/cpu_readout.py','src/run_study.py','src/medium_v2.py','environment.lock']
    manifest={p:digest(root/p) for p in paths}
    return hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest(),manifest


def one_arm(job):
    root=Path(job['root']);outdir=Path(job['output']);config=job['config']
    context,arm,scenario=job['context'],job['arm'],job['scenario']
    token=hashlib.sha256(f'{context}|{arm}|{scenario}'.encode()).hexdigest()[:16]
    dest=outdir/'arms'/f'{token}.json';dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():
        saved=json.loads(dest.read_text())
        flux_path=outdir/saved.get('flux_file','missing')
        if (saved.get('input_fingerprint')==job['input_fingerprint'] and saved.get('status')=='optimal'
            and flux_path.is_file() and digest(flux_path)==saved.get('flux_sha256')):
            return saved
    record={'run_id':job['run_id'],'context_id':context,'arm':arm,'scenario':scenario,
            'input_fingerprint':job['input_fingerprint']}
    t0=time.monotonic()
    try:
        d=np.load(root/'data/scores_conservative.npz',allow_pickle=False)
        rxns=d['rxn'].astype(str).tolist();contexts=d['contexts'].astype(str).tolist()
        ci=contexts.index(job['reference_context'] if context=='__reference__' else context)
        model,meta=build_model(root,config,scenario)
        w=score_weights(model,rxns,d['A'][:,ci],job['k'],arm,config['epsilon'])
        out=select_flux(model,w,job['panel'],config)
        fpath=dest.with_suffix('.npz')
        np.savez_compressed(fpath,flux=out.pop('flux'),reaction_ids=np.array(out.pop('reaction_ids')),
                            reaction_costs=np.array([w[r.id] for r in model.reactions]))
        record.update(out);record.update({'status':'optimal','model':meta,'flux_file':str(fpath.relative_to(outdir)),
                                          'flux_sha256':digest(fpath)})
    except Exception as exc:
        record.update({'status':'failed','error':repr(exc),'traceback':traceback.format_exc()})
    record['seconds']=time.monotonic()-t0
    temp=dest.with_suffix('.tmp');temp.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');temp.replace(dest)
    return record


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partition',choices=['development','test'],default='development')
    ap.add_argument('--scenarios',nargs='+',default=['primary']);ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--max-contexts',type=int);ap.add_argument('--arms',nargs='+');ap.add_argument('--output',required=True)
    args=ap.parse_args();config=load_config(ROOT)
    if args.partition=='test':
        from src.verify_inputs import verify_protocol_lock
        verify_protocol_lock(ROOT)
    for scenario in args.scenarios:
        if scenario not in config['scenarios']:raise ValueError(scenario)
    metadata=rows(ROOT/'manifests/contexts.tsv')
    contexts=sorted(r['context_id'] for r in metadata if r['partition']==args.partition)
    dev=sorted(r['context_id'] for r in metadata if r['partition']=='development')
    if args.max_contexts:contexts=contexts[:args.max_contexts]
    if not contexts or not dev:raise ValueError('Empty requested partition or development scale set')
    mapped=rows(ROOT/'manifests/metabolite_map.tsv')
    panel=sorted({r['exchange_id'] for r in mapped if str(r['primary_eligible']).lower() in ['1','true','yes']})
    if not panel:raise ValueError('Empty chemistry panel')
    d=np.load(ROOT/'data/scores_conservative.npz',allow_pickle=False)
    names=d['contexts'].astype(str).tolist();a=d['A'][:,[names.index(c) for c in dev]]
    k=float(np.median(a[a>0]))
    fingerprint,manifest=input_fingerprint(ROOT,config)
    output=Path(args.output).resolve();output.mkdir(parents=True,exist_ok=True)
    run_id=os.environ.get('RUN_ID',time.strftime('cpu-%Y%m%dT%H%M%SZ',time.gmtime()))
    arms=args.arms or [f'magnitude_g{g:.2f}' for g in config['gammas']]+['ordinal']
    jobs=[]
    for scenario in args.scenarios:
        for context in contexts:
            for arm in arms:
                jobs.append({'root':str(ROOT),'output':str(output),'config':config,'context':context,
                    'arm':arm,'scenario':scenario,'panel':panel,'k':k,'reference_context':dev[0],
                    'input_fingerprint':fingerprint,'run_id':run_id})
        jobs.append({**jobs[-1],'context':'__reference__','arm':'uniform_pfba'})
    print(json.dumps({'run_id':run_id,'partition':args.partition,'contexts':len(contexts),'targets':len(panel),
                      'jobs':len(jobs),'k':k,'workers':args.workers}),flush=True)
    t0=time.monotonic();results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(one_arm,j) for j in jobs]
        for f in as_completed(futures):
            rec=f.result();results.append(rec)
            print(f"{len(results)}/{len(jobs)} {rec['scenario']} {rec['context_id']} {rec['arm']}: {rec['status']} {rec['seconds']:.1f}s",flush=True)
    results.sort(key=lambda r:(r['scenario'],r['context_id'],r['arm']))
    with open(output/'predictions.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['run_id','context_id','arm','exchange_id','point','range_lo','range_hi','scenario','status'])
        for rec in results:
            for rid in panel:
                pt=rec.get('point',{}).get(rid,'');ran=rec.get('ranges',{}).get(rid,['',''])
                w.writerow([rec['run_id'],rec['context_id'],rec['arm'],rid,pt,*ran,rec['scenario'],rec['status']])
    summary={'run_id':run_id,'partition':args.partition,'config':config,'input_fingerprint':fingerprint,
        'inputs':manifest,'k':k,'contexts':contexts,'chemistry_exchange_ids':panel,'jobs':len(jobs),
        'completed':len(results),'failed':sum(r['status']!='optimal' for r in results),
        'seconds':time.monotonic()-t0,'arms':[{k:v for k,v in r.items() if k not in ['point','ranges','solve_ledger','coordinate_order','traceback']} for r in results]}
    (output/'run_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('FINISHED',summary['failed'],'failures',round(summary['seconds'],1),'seconds',flush=True)
    if summary['failed']:raise SystemExit(1)


if __name__=='__main__':main()
