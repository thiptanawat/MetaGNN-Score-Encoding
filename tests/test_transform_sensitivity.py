import csv
from src.transform_sensitivity import summarize
from src.transform_sensitivity import range_diagnostics
from src.evaluation import MAGNITUDE_ARMS


def test_strict_reversal_is_distinct_from_tie_transition(tmp_path):
    p=tmp_path/'pred.tsv'; rows=[]
    contexts={'a':{'origin_group':'a'},'b':{'origin_group':'b'},'same_a':{'origin_group':'a'}}
    panel=[{'exchange_id':'x','metabolite_id':'x'},{'exchange_id':'y','metabolite_id':'y'}]
    for c in contexts:
        for i,arm in enumerate(MAGNITUDE_ARMS):
            for rid in ['x','y']:
                v=0. if c!='b' else ([-1.,-1.,0.,1.,1.][i] if rid=='x' else [0.,0.,1.,1.,1.][i])
                rows.append(dict(context_id=c,arm=arm,exchange_id=rid,point=v,scenario='s',status='optimal'))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    r=summarize(p,'s',contexts,panel)
    assert r['all_distinct_origin_context_target_pairs']==4
    assert r['strict_order_reversals_across_grid']==2
    assert r['tie_transition_without_strict_reversal_pairs']==2
    assert r['changed_cells']==2


def test_constant_narrow_and_wide_ranges_are_distinct(tmp_path):
    p=tmp_path/'pred.tsv';rows=[]
    contexts={'a':{'origin_group':'a'},'b':{'origin_group':'b'}}
    panel=[{'exchange_id':'x'},{'exchange_id':'y'}]
    for c in contexts:
        for rid in ['x','y']:
            rows.append(dict(context_id=c,arm='ordinal',exchange_id=rid,point=0.,range_lo=0.,
                             range_hi=1e-8 if rid=='x' else 1.,scenario='s',status='optimal'))
    for rid in ['x','y']:
        rows.append(dict(context_id='__reference__',arm='uniform_pfba',exchange_id=rid,point=2.,
                         range_lo=2.,range_hi=2.,scenario='s',status='optimal'))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    r=range_diagnostics(p,'s',contexts,panel,1e-5,1e-5)['per_arm']
    assert r['ordinal']['within_context_width_above_tolerance_count']==2
    assert r['ordinal']['constant_prediction_targets']==2
    assert r['ordinal']['all_zero_prediction_targets']==2
    assert r['uniform_pfba']['context_target_cells']==4
    assert r['uniform_pfba']['within_context_width_above_tolerance_count']==0
    assert r['uniform_pfba']['all_zero_prediction_targets']==0
