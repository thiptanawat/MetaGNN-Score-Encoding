"""New outcome-free saved-prediction audit. Reads no CORE observations or solver."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import argparse, csv, gzip, hashlib, io, json, math, time

SCENARIOS = ('primary', 'half_serum', 'lower_task')
POWER = ('magnitude_g0.50','magnitude_g0.80','magnitude_g1.00','magnitude_g1.25','magnitude_g2.00')
ARM_SETS = {'five_power': POWER, 'six_including_ordinal': POWER + ('ordinal',)}
EPSILON = 1e-5
NUMERICAL_TOLERANCE = 1e-6

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with Path(path).open() as f: return list(csv.DictReader(f, delimiter='\t'))

def oriented(point, lo, hi, factor):
    point, lo, hi, factor = map(float, (point, lo, hi, factor))
    if not all(map(math.isfinite, (point, lo, hi, factor))) or factor not in (-1.0, 1.0):
        raise ValueError('Nonfinite value or invalid orientation')
    inverted = lo > hi
    if inverted:
        if lo-hi > NUMERICAL_TOLERANCE: raise ValueError('Material interval inversion')
        lo = hi = (lo+hi)/2
    if point < lo-NUMERICAL_TOLERANCE or point > hi+NUMERICAL_TOLERANCE:
        raise ValueError('Point outside saved range')
    return ((point, lo, hi) if factor == 1 else (-point, -hi, -lo)), inverted

def profile_pairs(contexts):
    keys = sorted(contexts)
    return [(a,b) for i,a in enumerate(keys) for b in keys[i+1:]
            if contexts[a] != contexts[b]]

def decisions(first, second):
    if len(first) != len(second) or not first: raise ValueError('Arm count mismatch')
    diffs = [a[0]-b[0] for a,b in zip(first,second)]
    signs = [1 if d>EPSILON else -1 if d < -EPSILON else 0 for d in diffs]
    if all(s==0 for s in signs): category='unanimous_ties'
    elif all(s==1 for s in signs): category='unanimous_positive'
    elif all(s==-1 for s in signs): category='unanimous_negative'
    elif 1 in signs and -1 in signs: category='opposing_nonzero_signs'
    else: category='mixed_ties_one_direction'
    matched_pos = min(a[1]-b[2] for a,b in zip(first,second))
    matched_neg = min(b[1]-a[2] for a,b in zip(first,second))
    independent_pos = min(a[1] for a in first)-max(b[2] for b in second)
    independent_neg = min(b[1] for b in second)-max(a[2] for a in first)
    matched = 1 if matched_pos>EPSILON else -1 if matched_neg>EPSILON else 0
    independent = 1 if independent_pos>EPSILON else -1 if independent_neg>EPSILON else 0
    for call in (matched, independent):
        if call and any(s != call for s in signs):
            raise ValueError('Interval call disagrees with a strict point sign')
    if independent and independent != matched:
        raise ValueError('Independent-encoding support is not a subset of matched support')
    return {'point_signs': ';'.join(map(str,signs)), 'point_category':category,
            'matched_interval_sign':matched, 'independent_interval_sign':independent,
            'matched_positive_margin':matched_pos,'matched_negative_margin':matched_neg,
            'independent_positive_margin':independent_pos,'independent_negative_margin':independent_neg}

def counts_to_row(counts):
    n = counts['pairs']
    out = dict(counts)
    for key in ['unanimous_positive','unanimous_negative','unanimous_ties','opposing_nonzero_signs',
                'mixed_ties_one_direction','matched_supported','independent_supported']:
        out.setdefault(key,0)
    out['point_unanimous_nontied'] = out['unanimous_positive']+out['unanimous_negative']
    for key in ['point_unanimous_nontied','unanimous_ties','matched_supported','independent_supported']:
        out[key+'_coverage'] = out[key]/n if n else None
    return out

def run(root, output):
    started = datetime.now(timezone.utc).isoformat(); t0=time.monotonic()
    here = Path(__file__).resolve().parent
    sources = {'plan':root/'protocol/analysis_plan.json', 'contexts':root/'manifests/contexts.tsv',
               'predictions':root/'results/reserved_v3/predictions.tsv',
               'amendment':here/'AMENDMENT.md','audit_script':Path(__file__).resolve(),
               'tests':here/'test_agreement_audit.py'}
    input_hash = {str(p):sha(p) for p in sources.values()}
    plan = json.loads(sources['plan'].read_text())
    if (plan['point_tolerance'] != EPSILON or plan['interval_tolerance'] != EPSILON or
        plan['numerical_tolerance'] != NUMERICAL_TOLERANCE):
        raise ValueError('Recorded thresholds differ from amendment')
    if sha(sources['contexts']) != plan['input_sha256']['manifests/contexts.tsv']:
        raise ValueError('Context input no longer matches recorded plan')
    panel = plan['primary_targets']; targets = {r['exchange_id']:r for r in panel}
    if len(panel) != 52 or len(targets) != 52: raise ValueError('Unexpected fixed panel')
    metadata = [r for r in read_tsv(sources['contexts']) if r['partition']=='test']
    contexts={r['context_id']:r['origin_group'] for r in metadata}
    if len(contexts)!=47 or len(set(contexts.values()))!=44 or len(metadata)!=47:
        raise ValueError('Unexpected profile/origin allocation')
    pairs = profile_pairs(contexts); records={}; inversions=0
    for r in read_tsv(sources['predictions']):
        if (r['scenario'] not in SCENARIOS or r['arm'] not in ARM_SETS['six_including_ordinal']
            or r['context_id'] not in contexts or r['exchange_id'] not in targets): continue
        key=(r['scenario'],r['context_id'],r['exchange_id'],r['arm'])
        if key in records or r['status']!='optimal': raise ValueError('Duplicate or failed prediction')
        records[key],inv=oriented(r['point'],r['range_lo'],r['range_hi'],targets[r['exchange_id']]['sign_factor'])
        inversions+=inv
    expected={(s,c,e,a) for s in SCENARIOS for c in contexts for e in targets for a in ARM_SETS['six_including_ordinal']}
    if set(records)!=expected: raise ValueError('Incomplete selected grid')
    output.mkdir(parents=True,exist_ok=True)
    ledger = output/'agreement_pair_ledger.tsv.gz'
    if ledger.exists(): raise FileExistsError('Refuse to overwrite a completed audit')
    details=[]; totals=[]; examples={}
    keys=['scenario','arm_set','metabolite_id','exchange_id','first_context','second_context',
          'first_origin','second_origin','point_signs','point_category','matched_interval_sign',
          'independent_interval_sign','matched_positive_margin','matched_negative_margin',
          'independent_positive_margin','independent_negative_margin']
    with ledger.open('wb') as raw, gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz, io.TextIOWrapper(gz,newline='') as f:
        writer=csv.DictWriter(f,keys,delimiter='\t');writer.writeheader()
        for scenario in SCENARIOS:
            for arm_set,arms in ARM_SETS.items():
                total=Counter()
                for target in panel:
                    counts=Counter();e=target['exchange_id']
                    for first,second in pairs:
                        one=[records[(scenario,first,e,a)] for a in arms]
                        two=[records[(scenario,second,e,a)] for a in arms]
                        d=decisions(one,two)
                        row={'scenario':scenario,'arm_set':arm_set,'metabolite_id':target['metabolite_id'],
                             'exchange_id':e,'first_context':first,'second_context':second,
                             'first_origin':contexts[first],'second_origin':contexts[second],**d}
                        writer.writerow(row)
                        counts['pairs']+=1;counts[d['point_category']]+=1
                        counts['matched_supported']+=bool(d['matched_interval_sign'])
                        counts['independent_supported']+=bool(d['independent_interval_sign'])
                        kind='independent_supported' if d['independent_interval_sign'] else ('matched_only' if d['matched_interval_sign'] else d['point_category'])
                        exkey=f'{scenario}/{arm_set}/{kind}'
                        if exkey not in examples:
                            examples[exkey]={'ledger_row':row,'arms':list(arms),'first_point_lo_hi':one,'second_point_lo_hi':two}
                    total.update(counts)
                    details.append({'scenario':scenario,'arm_set':arm_set,'metabolite_id':target['metabolite_id'],
                                    'exchange_id':e,**counts_to_row(counts)})
                totals.append({'scenario':scenario,'arm_set':arm_set,**counts_to_row(total)})
    for name,rows in [('agreement_summary.tsv',totals),('agreement_per_target.tsv',details)]:
        with (output/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    result={'status':'completed_outcome_free_post_analysis_audit','started_at_utc':started,
            'finished_at_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-t0,
            'new_lp_solves':0,'core_outcomes_read':False,'accuracy_calculated':False,
            'threshold':EPSILON,'numerical_tolerance':NUMERICAL_TOLERANCE,
            'profiles':len(contexts),'origins':len(set(contexts.values())),'targets':len(panel),
            'distinct_origin_profile_pairs':len(pairs),'selected_prediction_records':len(records),
            'tolerance_level_inverted_intervals_normalized':inversions,
            'interval_point_sign_inconsistencies':0,'support_subset_inconsistencies':0,
            'input_sha256':input_hash,'summary':totals,'per_target':details,'examples':examples,
            'limitations':['Agreement is not accuracy or biological validation.',
                'Unanimous ties are abstention, not informative coverage.',
                'Numerical ranges are not statistical confidence intervals or exact certificates.',
                'Independent-encoding support allows a different encoding for each profile.',
                'Dependent profile pairs and scenarios are not independent sample sizes.',
                '52 targets were originally selected using development outcomes; this audit reads no new outcomes.']}
    (output/'agreement_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (output/'SHA256SUMS.txt').open('w') as f:
        for p in sorted(output.iterdir()):
            if p.is_file() and p.name!='SHA256SUMS.txt':f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps({'seconds':result['elapsed_seconds'],'summary':totals},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.root.resolve(),a.output.resolve())
