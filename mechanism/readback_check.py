"""Read back completed pilot artifacts; no source-model or solver execution."""
from pathlib import Path
import csv,hashlib,json
import numpy as np

def main():
    root=Path(__file__).resolve().parent
    result=json.loads((root/'summary.json').read_text());manifest=json.loads((root/'run_manifest.json').read_text())
    records=[];worst=0.;full_vectors=0
    for p in sorted(root.glob('*.json')):
        record=json.loads(p.read_text())
        if record.get('stage') not in ['B0','B1'] or record.get('status')!='complete':continue
        fp=root/record['extrema_file']
        assert hashlib.sha256(fp.read_bytes()).hexdigest()==record['extrema_sha256']
        with np.load(fp,allow_pickle=False) as d:
            ids={v:i for i,v in enumerate(d['reaction_ids'].astype(str))}
            labels=d['endpoint'].astype(str).tolist();v=d['flux']
            assert v.shape[0]==len(labels)==104 and np.isfinite(v).all()
            for i,label in enumerate(labels):
                stage,rid=label.split(':',1);expected=record['ranges'][rid][0 if stage=='range_min' else 1]
                worst=max(worst,abs(v[i,ids[rid]]-expected))
            full_vectors+=len(labels)
        records.append(p.name)
    assert worst<=1e-6
    for path,h in manifest['inputs_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h
    with (root/'range_stage_ledger.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    summaries=[];nonmonotone=[]
    for arm in ['magnitude_g1.00','ordinal']:
        selected=[r for r in rows if r['arm']==arm]
        s={'arm':arm,'profile_target_rows':len(selected)}
        for stage in ['B0','B1','B2']:s[stage+'_fixed']=sum(r[stage+'_fixed']=='True' for r in selected)
        for label in ['already_fixed_B0','first_fixed_B1','first_fixed_B2','still_variable_B2']:
            s[label]=sum(r['first_fixed_stage']==label for r in selected)
        for key in ['B0_to_B1_width_contraction','B1_to_B2_width_contraction']:
            values=[float(r[key]) for r in selected]
            s[key+'_above_threshold']=sum(v>1e-5 for v in values)
        summaries.append(s)
    for r in rows:
        ranges=[[float(r[stage+'_low']),float(r[stage+'_high'])] for stage in ['B0','B1','B2']]
        for outer,inner in zip(ranges,ranges[1:]):assert inner[0]>=outer[0]-1e-6 and inner[1]<=outer[1]+1e-6
        f=[r[stage+'_fixed']=='True' for stage in ['B0','B1','B2']]
        if any(a and not b for a,b in zip(f,f[1:])):nonmonotone.append({k:r[k] for k in ['context','arm','exchange_id']})
    receipt={'readback_status':'passed','complete_records':len(records),'new_full_extremum_vectors':full_vectors,
             'max_extremum_coordinate_difference':worst,'all_source_hashes_match':True,
             'ledger_rows':len(rows),'nonmonotone_fixedness_threshold_cases':nonmonotone,
             'arm_aggregates':summaries,'B0_fixed_once':sum(hi-lo<=1e-5 for lo,hi in json.loads((root/'B0.json').read_text())['ranges'].values())}
    (root/'readback_check.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
