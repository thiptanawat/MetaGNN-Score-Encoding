"""Descriptive prediction sensitivity to fixed-rank score spacing; no CORE values."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from src.cpu_model import digest
from src.evaluation import MAGNITUDE_ARMS,read_tsv


def range_diagnostics(predictions,scenario,contexts,panel,point_tolerance,interval_tolerance):
    """Separate within-context width from lack of variation between contexts."""
    rows=[r for r in read_tsv(predictions) if r['scenario']==scenario]
    index={(r['context_id'],r['arm'],r['exchange_id']):r for r in rows}
    arms=sorted({r['arm'] for r in rows})
    output={}
    for arm in arms:
        points=[];widths=[]
        for c in contexts:
            values=[]
            for target in panel:
                rid=target['exchange_id']
                r=index.get((c,arm,rid)) or index.get(('__reference__',arm,rid))
                if r is None or r['status']!='optimal':raise ValueError(f'Missing successful range {c}/{arm}/{rid}')
                point,lo,hi=(float(r[k]) for k in ('point','range_lo','range_hi'))
                if not np.isfinite([point,lo,hi]).all() or hi<lo-interval_tolerance:
                    raise ValueError(f'Invalid range {c}/{arm}/{rid}')
                values.append(point);widths.append(max(0.,hi-lo))
            points.append(values)
        arr=np.asarray(points);widths=np.asarray(widths)
        spans=np.ptp(arr,axis=0)
        output[arm]={'context_target_cells':len(widths),
          'within_context_width_above_tolerance_count':int((widths>interval_tolerance).sum()),
          'within_context_width_quantiles':{str(q):float(np.quantile(widths,q)) for q in [0,.5,.9,.95,1.]},
          'constant_prediction_targets':int((spans<=point_tolerance).sum()),
          'all_zero_prediction_targets':int((np.max(np.abs(arr),axis=0)<=point_tolerance).sum()),
          'total_targets':len(panel)}
    return {'point_tolerance':point_tolerance,'interval_tolerance':interval_tolerance,
      'interpretation':'Descriptive numerical diagnostics, not certificates of exact uniqueness. Across-context interval overlap can reflect constant narrow predictions, wide within-context ranges, or both.',
      'per_arm':output}


def summarize(predictions,scenario,contexts,panel,tolerance=1e-6):
    rows=[r for r in read_tsv(predictions) if r['scenario']==scenario and r['context_id'] in contexts and r['arm'] in MAGNITUDE_ARMS]
    index={(r['context_id'],r['arm'],r['exchange_id']):r for r in rows}
    spans=[];detail=[];full=[]
    for c in contexts:
        for target in panel:
            rid=target['exchange_id'];vals=[]
            for arm in MAGNITUDE_ARMS:
                r=index.get((c,arm,rid))
                if r is None or r['status']!='optimal':raise ValueError(f'Missing successful prediction {c}/{arm}/{rid}')
                vals.append(float(r['point']))
            span=float(np.ptp(vals));spans.append(span);full.append(vals)
            detail.append({'context_id':c,'metabolite_id':target['metabolite_id'],'exchange_id':rid,'prediction_span':span,'changed_above_tolerance':span>tolerance})
    spans=np.array(spans);arr=np.array(full).reshape(len(contexts),len(panel),len(MAGNITUDE_ARMS))
    strict=0;tiechange=0;total=0
    for i,c in enumerate(contexts):
        for j,d in enumerate(contexts):
            if j<=i or contexts[c]['origin_group']==contexts[d]['origin_group']:continue
            differences=arr[i]-arr[j]
            signs=np.where(np.abs(differences)<=tolerance,0,np.sign(differences))
            strict+=int(((signs.min(axis=1)<0)&(signs.max(axis=1)>0)).sum())
            tiechange+=int(((signs.min(axis=1)!=signs.max(axis=1))&~((signs.min(axis=1)<0)&(signs.max(axis=1)>0))).sum())
            total+=len(panel)
    return {'contexts':len(contexts),'targets':len(panel),'context_target_cells':len(spans),'prediction_tolerance':tolerance,
      'changed_cells':int((spans>tolerance).sum()),'changed_fraction':float((spans>tolerance).mean()),
      'prediction_span_quantiles':{str(q):float(np.quantile(spans,q)) for q in [0,.5,.9,.95,1.]},
      'all_distinct_origin_context_target_pairs':total,'strict_order_reversals_across_grid':strict,
      'strict_order_reversal_fraction':strict/total if total else None,
      'tie_transition_without_strict_reversal_pairs':tiechange,
      'interpretation':'Descriptive across all different-origin context pairs, irrespective of observed CORE separability; not accuracy or biological validation.',
      'detail':detail}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--predictions',type=Path,required=True);ap.add_argument('--partition',choices=['development','test'],default='development');ap.add_argument('--scenario',default='primary');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if a.partition=='test':
        from src.verify_inputs import verify_protocol_lock
        verify_protocol_lock(a.root)
    plan=json.loads((a.root/'protocol/analysis_plan.json').read_text())
    contexts={r['context_id']:r for r in read_tsv(a.root/'manifests/contexts.tsv') if r['partition']==a.partition}
    contexts=dict(sorted(contexts.items()))
    result={'scenario':a.scenario,'partition':a.partition,'predictions_sha256':digest(a.predictions),'analysis_plan_sha256':digest(a.root/'protocol/analysis_plan.json'),
      'primary_fixed_panel':summarize(a.predictions,a.scenario,contexts,plan['primary_targets'],plan['point_tolerance']),
      'broader_chemistry':summarize(a.predictions,a.scenario,contexts,plan['chemistry_targets'],plan['point_tolerance'])}
    result['range_diagnostics_primary']=range_diagnostics(a.predictions,a.scenario,contexts,plan['primary_targets'],plan['point_tolerance'],plan['interval_tolerance'])
    result['range_diagnostics_broader']=range_diagnostics(a.predictions,a.scenario,contexts,plan['chemistry_targets'],plan['point_tolerance'],plan['interval_tolerance'])
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result['primary_fixed_panel'].items() if k!='detail'},indent=2))

if __name__=='__main__':main()
