"""Bounded outcome-free B0/B1 computation; reuse verified saved B2 only."""
import sys
sys.dont_write_bytecode = True
import argparse, csv, hashlib, json, math, time, traceback
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.sparse import coo_matrix

EPS=1e-5
TOL=1e-6

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_tsv(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def dump(path, value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def normalized_range(lo,hi):
    lo,hi=float(lo),float(hi)
    if not math.isfinite(lo+hi):raise ValueError('Nonfinite range')
    if lo>hi:
        if lo-hi>TOL:raise ValueError('Materially inverted range')
        lo=hi=(lo+hi)/2
    return [lo,hi]

def nested_violation(outer,inner):return max(0.,outer[0]-inner[0],inner[1]-outer[1])

def stage_label(b0,b1,b2):
    fixed=[b[1]-b[0]<=EPS for b in [b0,b1,b2]]
    if fixed[0]:return 'already_fixed_B0'
    if fixed[1]:return 'first_fixed_B1'
    if fixed[2]:return 'first_fixed_B2'
    return 'still_variable_B2'

def worker(job):
    root,out=Path(job['root']),Path(job['out']);dest=out/(job['token']+'.json')
    record={'stage':job['stage'],'context':job.get('context'),'arm':job.get('arm'),
            'status':'running','ranges':{},'solve_ledger':[],'audits':[],
            'started_at_utc':datetime.now(timezone.utc).isoformat()}
    t0=time.monotonic();vectors=[];names=[]
    try:
        if time.monotonic()>=job['deadline']:raise TimeoutError('Global budget exhausted before job')
        sys.path.insert(0,str(root))
        from src.cpu_model import build_model,score_weights,set_absolute_cost
        config=json.loads((root/'configs/study.json').read_text())
        model,meta=build_model(root,config,job.get('scenario','primary'))
        record['task']={k:meta[k] for k in ['g_max','task_flux','objective','scenario']}
        ids=[r.id for r in model.reactions]; lookup={r:i for i,r in enumerate(ids)}
        mets={m.id:i for i,m in enumerate(model.metabolites)};ii=[];jj=[];vv=[]
        for j,r in enumerate(model.reactions):
            for m,v in r.metabolites.items():ii.append(mets[m.id]);jj.append(j);vv.append(v)
        S=coo_matrix((vv,(ii,jj)),shape=(len(mets),len(ids))).tocsr()
        lower=np.array([r.lower_bound for r in model.reactions]);upper=np.array([r.upper_bound for r in model.reactions])
        weights=None;primary_cap=None

        def vector_audit(flux):
            if not np.isfinite(flux).all():raise ValueError('Nonfinite full flux')
            result={'mass_balance_residual':float(np.max(np.abs(S@flux))),
                    'bound_violation':float(max(0,np.max(lower-flux),np.max(flux-upper)))}
            if weights is not None:
                result['primary_cost']=float(np.dot(weights,np.abs(flux)))
                result['primary_cap_violation']=max(0.,result['primary_cost']-primary_cap)
            bad={k:v for k,v in result.items() if (k.endswith('residual') or k.endswith('violation')) and v>TOL}
            if bad:raise ValueError(f'Vector validation failed {bad}')
            return result

        def solve(stage,rid='',keep=False):
            for attempt in range(2):
                remaining=job['deadline']-time.monotonic()
                if remaining<=0:raise TimeoutError('Global runtime budget exhausted')
                model.solver.configuration.timeout=max(1,min(60,int(remaining)))
                started=time.monotonic();value=model.slim_optimize(error_value=float('nan'))
                status=model.solver.status
                record['solve_ledger'].append({'stage':stage,'exchange':rid,'attempt':attempt,
                    'status':status,'objective':float(value) if np.isfinite(value) else None,
                    'seconds':time.monotonic()-started})
                if status=='optimal' and np.isfinite(value):
                    flux=np.array([r.forward_variable.primal-r.reverse_variable.primal for r in model.reactions])
                    if keep:
                        aud=vector_audit(flux);record['audits'].append({'stage':stage,'exchange':rid,**aud})
                        names.append(f'{stage}:{rid}');vectors.append(flux)
                    return float(value)
                if attempt==0:
                    from swiglpk import glp_adv_basis
                    glp_adv_basis(model.solver.problem,0)
            raise RuntimeError(f'Nonoptimal {stage}/{rid}: {status}')

        if job['stage']=='B1':
            saved=json.loads(Path(job['saved_json']).read_text())
            with np.load(job['saved_npz'],allow_pickle=False) as d:
                saved_ids=d['reaction_ids'].astype(str).tolist();saved_flux=d['flux'].copy();saved_weights=d['reaction_costs'].copy()
            if ids!=saved_ids:raise ValueError('Saved/model reaction order mismatch')
            with np.load(root/'data/scores_conservative.npz',allow_pickle=False) as d:
                ci=d['contexts'].astype(str).tolist().index(job['context'])
                w=score_weights(model,d['rxn'].astype(str).tolist(),d['A'][:,ci],job['k'],job['arm'],config['epsilon'])
            weights=np.array([w[r] for r in ids]);difference=float(np.max(np.abs(weights-saved_weights)))
            if difference>1e-12:raise ValueError(f'Cost mismatch {difference}')
            for k in ['g_max','task_flux']:
                if abs(meta[k]-saved['model'][k])>1e-9:raise ValueError(f'Task/model mismatch {k}')
            if saved['model']['objective']!=meta['objective'] or saved['model']['scenario']!=job.get('scenario','primary'):raise ValueError('Wrong saved model')
            set_absolute_cost(model,w);primary_opt=solve('primary_reoptimization')
            difference_opt=abs(primary_opt-saved['primary_optimum'])
            if difference_opt>TOL:raise ValueError(f'Primary optimum mismatch {difference_opt}')
            primary_cap=float(saved['primary_cap'])
            cap=model.problem.Constraint(model.objective.expression,ub=primary_cap,name='pilot_primary_cap')
            model.add_cons_vars([cap]);model.solver.update()
            aud=vector_audit(saved_flux)
            secondary=float(np.abs(saved_flux).sum())
            if secondary>saved['secondary_cap']+TOL:raise ValueError('Saved B2 point fails secondary cap')
            coordinate_point_error=max(abs(saved['point'][rid]-saved_flux[lookup[rid]]) for rid in job['panel'])
            if coordinate_point_error>TOL:raise ValueError('Saved B2 full vector and point mismatch')
            record['comparability']={'cost_max_difference':difference,'primary_optimum_difference':difference_opt,
                'saved_primary_optimum':saved['primary_optimum'],'new_primary_optimum':primary_opt,
                'exact_reused_primary_cap':primary_cap,'saved_secondary_cap':saved['secondary_cap'],
                'saved_B2_vector_audit':aud,'saved_B2_secondary_cost':secondary,
                'saved_B2_point_difference':coordinate_point_error,'saved_B2_reoptimized':False}
            record['saved_B2_ranges']={r:normalized_range(*saved['ranges'][r]) for r in job['panel']}
            record['saved_B2_points']={r:float(saved['point'][r]) for r in job['panel']}

        for rid in job['panel']:
            expression=model.reactions.get_by_id(rid).flux_expression
            model.objective=model.problem.Objective(expression,direction='min');lo=solve('range_min',rid,True)
            model.objective=model.problem.Objective(expression,direction='max');hi=solve('range_max',rid,True)
            record['ranges'][rid]=normalized_range(lo,hi)
        record['status']='complete'
    except Exception as exc:
        record['status']='failed_or_incomplete';record['error']=repr(exc);record['traceback']=traceback.format_exc()
    if vectors:
        fluxpath=out/(job['token']+'_extrema.npz')
        np.savez_compressed(fluxpath,reaction_ids=np.array(ids),endpoint=np.array(names),flux=np.stack(vectors))
        record['extrema_file']=fluxpath.name;record['extrema_sha256']=sha(fluxpath)
    record['seconds']=time.monotonic()-t0;record['finished_at_utc']=datetime.now(timezone.utc).isoformat();dump(dest,record)
    return record

def main(root,out):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'run_manifest.json').exists():raise FileExistsError('Do not overwrite existing pilot')
    config=json.loads((root/'configs/study.json').read_text())
    plan=json.loads((root/'protocol/analysis_plan.json').read_text())
    summary=json.loads((root/'results/development_v3/run_summary.json').read_text())
    sources=['data/raw/Recon3D.json','configs/study.json','environment.lock','data/scores_conservative.npz',
        'manifests/contexts.tsv','manifests/metabolite_map.tsv','protocol/analysis_plan.json',
        'src/cpu_model.py','src/cpu_readout.py','src/medium_v2.py','src/run_study.py']
    hashes={str(root/p):sha(root/p) for p in sources}
    for name,h in summary['inputs'].items():
        if name not in sources:raise ValueError(f'Unexpected required original fingerprint file {name}')
        if hashes[str(root/name)]!=h:raise ValueError(f'Original run input differs: {name}')
    if config!=summary['config']:raise ValueError('Saved config mismatch')
    if sha(root/'manifests/contexts.tsv')!=plan['input_sha256']['manifests/contexts.tsv']:raise ValueError('Plan contexts mismatch')
    if plan['point_tolerance']!=EPS or plan['numerical_tolerance']!=TOL:raise ValueError('Tolerance mismatch')
    panel=sorted(t['exchange_id'] for t in plan['primary_targets'])
    if len(panel)!=52 or len(set(panel))!=52:raise ValueError('Unexpected panel')
    contexts=sorted(r['context_id'] for r in read_tsv(root/'manifests/contexts.tsv') if r['partition']=='development')
    if len(contexts)!=11:raise ValueError('Unexpected development allocation')
    with np.load(root/'data/scores_conservative.npz',allow_pickle=False) as d:
        names=d['contexts'].astype(str).tolist();a=d['A'][:,[names.index(c) for c in contexts]];k=float(np.median(a[a>0]))
    if k!=summary['k']:raise ValueError('Scale mismatch')
    selections=[]
    for context in contexts:
        for arm in ['magnitude_g1.00','ordinal']:
            token=hashlib.sha256(f'{context}|{arm}|primary'.encode()).hexdigest()[:16]
            saved_json=root/'results/development_v3/arms'/f'{token}.json';saved=json.loads(saved_json.read_text())
            saved_npz=root/'results/development_v3'/saved['flux_file']
            if (saved['status']!='optimal' or saved['context_id']!=context or saved['arm']!=arm or saved['scenario']!='primary'
                or saved['input_fingerprint']!=summary['input_fingerprint'] or sha(saved_npz)!=saved['flux_sha256']):
                raise ValueError('Saved B2 record identity/hash failure')
            final={(r['stage'],r['exchange_id']):r for r in saved['solve_ledger']}
            required=[('primary',''),('secondary',''),('final_feasibility','')]+[(s,r) for r in panel for s in ['range_min','range_max']]
            if any(final.get(x,{}).get('status')!='optimal' for x in required):raise ValueError('Missing successful original logical stage')
            for stage in ['primary','secondary']:
                prefix='cost' if stage=='primary' else 'secondary'
                expected=saved[f'{stage}_optimum']+config[f'{prefix}_absolute_allowance']+abs(saved[f'{stage}_optimum'])*config[f'{prefix}_relative_allowance']
                if abs(expected-saved[f'{stage}_cap'])>1e-12:raise ValueError('Saved cap formula mismatch')
            hashes[str(saved_json)]=sha(saved_json);hashes[str(saved_npz)]=sha(saved_npz)
            selections.append({'token':token,'context':context,'arm':arm,'stage':'B1',
                               'saved_json':str(saved_json),'saved_npz':str(saved_npz),'k':k})
    hashes[str(root/'results/development_v3/run_summary.json')]=sha(root/'results/development_v3/run_summary.json')
    hashes[str(Path(__file__).resolve())]=sha(__file__)
    hashes[str(out/'AMENDMENT.md')]=sha(out/'AMENDMENT.md')
    start=time.monotonic();deadline=start+480
    common={'root':str(root),'out':str(out),'panel':panel,'deadline':deadline}
    manifest={'started_at_utc':datetime.now(timezone.utc).isoformat(),'stage':'preflight_passed_before_LPs',
        'inputs_sha256':hashes,'expected_B0_jobs':1,'expected_B1_jobs':22,'workers':4,'runtime_budget_seconds':480,
        'profiles':contexts,'targets':panel,'outcomes_read':False,'original_files_mutated':False,
        'python':sys.version,'python_executable':sys.executable,'saved_B2_extrema_reoptimized':False}
    dump(out/'run_manifest.json',manifest)
    b0=worker({**common,'token':'B0','stage':'B0'});print('B0',b0['status'],round(b0['seconds'],2),flush=True)
    results=[]
    if b0['status']=='complete':
        with ProcessPoolExecutor(max_workers=4) as pool:
            jobs=[pool.submit(worker,{**common,**s}) for s in selections]
            for future in as_completed(jobs):
                r=future.result();results.append(r);print(len(results),r['context'],r['arm'],r['status'],round(r['seconds'],2),flush=True)
    rows=[];issues=[]
    for r in results:
        if r['status']!='complete':continue
        for rid in panel:
            z=b0['ranges'][rid];one=r['ranges'][rid];two=r['saved_B2_ranges'][rid];point=r['saved_B2_points'][rid]
            v01=nested_violation(z,one);v12=nested_violation(one,two)
            pv=max(0.,z[0]-point,point-z[1],one[0]-point,point-one[1],two[0]-point,point-two[1])
            if max(v01,v12,pv)>TOL:issues.append({'context':r['context'],'arm':r['arm'],'exchange_id':rid,'B0_B1':v01,'B1_B2':v12,'point':pv})
            w0,w1,w2=[x[1]-x[0] for x in [z,one,two]]
            rows.append({'context':r['context'],'arm':r['arm'],'exchange_id':rid,
                'B0_low':z[0],'B0_high':z[1],'B1_low':one[0],'B1_high':one[1],'B2_low':two[0],'B2_high':two[1],
                'B0_width':w0,'B1_width':w1,'B2_width':w2,'B0_fixed':w0<=EPS,'B1_fixed':w1<=EPS,'B2_fixed':w2<=EPS,
                'B0_to_B1_width_contraction':w0-w1,'B1_to_B2_width_contraction':w1-w2,
                'first_fixed_stage':stage_label(z,one,two),'B0_B1_nesting_violation':v01,
                'B1_B2_nesting_violation':v12,'saved_point_containment_violation':pv})
    if rows:
        with (out/'range_stage_ledger.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    grouped=[]
    for context in contexts:
        for arm in ['magnitude_g1.00','ordinal']:
            subset=[r for r in rows if r['context']==context and r['arm']==arm]
            if not subset:continue
            entry={'context':context,'arm':arm,'targets':len(subset)}
            for stage in ['B0','B1','B2']:entry[stage+'_fixed']=sum(r[stage+'_fixed'] for r in subset)
            for label in ['already_fixed_B0','first_fixed_B1','first_fixed_B2','still_variable_B2']:
                entry[label]=sum(r['first_fixed_stage']==label for r in subset)
            for key in ['B0_to_B1_width_contraction','B1_to_B2_width_contraction']:
                entry[key+'_above_threshold']=sum(r[key]>EPS for r in subset)
            grouped.append(entry)
    hashes_unchanged=all(sha(path)==value for path,value in hashes.items())
    completed=sum(r['status']=='complete' for r in results)
    result={'status':'complete_verified' if completed==22 and not issues and hashes_unchanged else 'partial_or_validation_failure',
        'completed_at_utc':datetime.now(timezone.utc).isoformat(),'seconds':time.monotonic()-start,
        'B0_status':b0['status'],'B1_completed':completed,'B1_failed_or_incomplete':sum(r['status']!='complete' for r in results),
        'B1_unstarted':22-len(results),'reused_B2_records':completed,'target_context_arm_rows':len(rows),
        'validation_issues':issues,'all_input_hashes_unchanged':hashes_unchanged,'summary_by_profile_arm':grouped,
        'outcomes_read':False,'point_tolerance':EPS,'nesting_tolerance':TOL,
        'interpretation':'Conditional within-profile coordinate restriction, not between-profile constancy or a biological bottleneck.'}
    dump(out/'summary.json',result);print(json.dumps({k:v for k,v in result.items() if k!='summary_by_profile_arm'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.root.resolve(),a.output.resolve())
