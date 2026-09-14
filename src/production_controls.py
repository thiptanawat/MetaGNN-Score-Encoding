"""Outcome-blind checks through the complete production CPU point readout."""
from pathlib import Path
import argparse
import csv
import json
import time
import numpy as np
import cobra
from src.cpu_model import build_model, load_config, score_weights, digest
from src.cpu_readout import select_flux

ROOT=Path(__file__).resolve().parents[1]


def read_rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(ROOT));args=ap.parse_args()
    root=Path(args.root);cfg=load_config(root);t0=time.monotonic()
    ctx=read_rows(root/'manifests/contexts.tsv');dev=sorted(r['context_id'] for r in ctx if r['partition']=='development')
    mp=read_rows(root/'manifests/metabolite_map.tsv')
    panel=sorted({r['exchange_id'] for r in mp if str(r['primary_eligible']).lower() in ['1','true','yes']})
    d=np.load(root/'data/scores_conservative.npz',allow_pickle=False);names=d['contexts'].astype(str).tolist();rxns=d['rxn'].astype(str).tolist()
    aa=d['A'][:,[names.index(c) for c in dev]];k=float(np.median(aa[aa>0]));a=d['A'][:,names.index('CO:HT29')]
    results={'context':'CO:HT29','scale_k':k,'panel_size':len(panel),'config':cfg,
        'inputs':{p:digest(root/p) for p in ['configs/study.json','data/scores_conservative.npz','manifests/metabolite_map.tsv','src/cpu_model.py','src/cpu_readout.py','src/production_controls.py','src/medium_v2.py','manifests/contexts.tsv','data/raw/Recon3D.json','environment.lock']},'checks':{}}
    def fresh():return build_model(root,cfg,'primary')[0]
    def delta(p,q):return max(abs(p[r]-q[r]) for r in panel)
    def solve(m,w,factors=None,sec=None):return select_flux(m,w,panel,cfg,coordinate_factors=factors,secondary_costs=sec)
    m=fresh();w=score_weights(m,rxns,a,k,'magnitude_g1.00',cfg['epsilon']);base=solve(m,w)
    print('baseline final readout',round(base['seconds'],2),'s',flush=True)
    rep=solve(fresh(),w)
    change=delta(base['point'],rep['point']);results['checks']['fresh_model_repeatability']={'max_point_change':change,'pass':change<=cfg['prediction_tie_tolerance']}
    print('fresh repeat',change,flush=True)
    z=a/(a+k)
    for gamma in [.5,2.0]:
        zr=(z**gamma)**(1/gamma);ar=k*zr/(1-zr)
        new=fresh();wr=score_weights(new,rxns,ar,k,'magnitude_g1.00',cfg['epsilon']);out=solve(new,wr)
        diff=delta(base['point'],out['point']);max_input=float(np.max(abs(zr-z)))
        results['checks'][f'inverse_gamma_{gamma}']={'max_input_change':max_input,'max_point_change':diff,'pass':max_input<1e-12 and diff<=cfg['prediction_tie_tolerance']}
        print('inverse',gamma,diff,flush=True)
    # Active coordinates ensure the control is not a trivial zero-flux flip.
    active=sorted(rid for rid,v in zip(base['reaction_ids'],base['flux']) if abs(v)>1e-4 and not m.reactions.get_by_id(rid).boundary and rid!=cfg['objective'])
    exchange=sorted([rid for rid in panel if abs(base['point'][rid])>1e-4])
    if not active or not exchange:raise RuntimeError('No active coordinates available for nontrivial controls')
    for rid in [active[0],exchange[0]]:
        new=fresh();r=new.reactions.get_by_id(rid);stoich=dict(r.metabolites);lo,hi=r.bounds
        r.add_metabolites({met:-2*v for met,v in stoich.items()});r.bounds=(-hi,-lo)
        factors={rid:-1} if rid in panel else None
        out=solve(new,w,factors=factors);diff=delta(base['point'],out['point'])
        results['checks'][f'orientation_{rid}']={'max_point_change':diff,'pass':diff<=cfg['prediction_tie_tolerance']}
        print('orientation',rid,diff,flush=True)
    for rid in [active[0],exchange[0]]:
        new=fresh();r=new.reactions.get_by_id(rid);stoich=dict(r.metabolites);lo,hi=r.bounds;scale=10.
        r.add_metabolites({met:v*(1/scale-1) for met,v in stoich.items()});r.bounds=(lo*scale,hi*scale)
        ww=dict(w);ww[rid]/=scale;sec={r.id:1.0 for r in new.reactions};sec[rid]/=scale
        factors={rid:1/scale} if rid in panel else None
        out=solve(new,ww,factors=factors,sec=sec);diff=delta(base['point'],out['point'])
        results['checks'][f'flux_unit_{rid}']={'factor':scale,'max_point_change':diff,'pass':diff<=cfg['prediction_tie_tolerance']}
        print('flux units',rid,diff,flush=True)
    # Closed-medium growth and ATP demand checks are model QC, not CORE prediction.
    new=fresh();new.reactions.get_by_id(cfg['objective']).bounds=(0,1000)
    for r in new.reactions:
        if r.id.startswith('EX_'):r.lower_bound=0
    new.objective=cfg['objective'];growth=float(new.slim_optimize());gs=new.solver.status
    r=cobra.Reaction('CPU_ATP_CHECK');r.add_metabolites({new.metabolites.get_by_id(k):v for k,v in {'atp_c':-1,'h2o_c':-1,'adp_c':1,'pi_c':1,'h_c':1}.items()});r.bounds=(0,1000);new.add_reactions([r]);new.objective=r
    atp=float(new.slim_optimize());ats=new.solver.status
    results['checks']['closed_medium']={'max_task':growth,'max_ATP_demand':atp,'status_task':gs,'status_ATP':ats,'pass':gs=='optimal' and ats=='optimal' and abs(growth)<1e-7 and abs(atp)<1e-7}
    results['baseline_audit']=base['audit'];results['seconds']=time.monotonic()-t0
    results['all_pass']=all(c['pass'] for c in results['checks'].values())
    (root/'results/production_controls.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    print('ALL PASS',results['all_pass'],round(results['seconds'],1),'s',flush=True)
    if not results['all_pass']:raise SystemExit(1)


if __name__=='__main__':main()
