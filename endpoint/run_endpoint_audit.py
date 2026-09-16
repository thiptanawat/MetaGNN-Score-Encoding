from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib,math,statistics,xlrd
import os
B=Path(__file__).resolve().parent
S=B.parent  # the repository is the study root
OUT=Path(os.environ.get('ENDPOINT_OUT',S/'results/endpoint'));OUT.mkdir(parents=True,exist_ok=True)
RAW=S/'data/raw/core/NIHMS419088-supplement-Database_S1.xls'
assert hashlib.sha256((B/'AUDIT_AMENDMENT.json').read_bytes()).hexdigest()==(B/'AUDIT_AMENDMENT.sha256').read_text().split()[0]
b=xlrd.open_workbook(RAW);s=b.sheet_by_name('Concentrations');hdr=s.row_values(0)
bases=[i for i,x in enumerate(hdr) if str(x).startswith('Baseline')];cultures=[i for i in range(3,s.ncols) if i not in bases]
assert len(bases)==10 and len(cultures)==120
contexts=list(csv.DictReader((S/'manifests/contexts.tsv').open(),delimiter='\t'))
bylabel={r['core_label']:r for r in contexts if r['has_core']=='True'}
assert len(bylabel)==60
human_path=B/'Source_Quality_Provenance.json'
human=set(json.loads(human_path.read_text())['source_cell_lines'])
counts={};assay=[];cult=[];seen={};profiles=[]
for col in cultures:
 label=hdr[col];seen[label]=seen.get(label,0)+1;rep=seen[label];c=bylabel[label]
 volume=50 if label in ['SR','MOLT-4','HL-60(TB)','K562','RPMI 8226','CCRF-CEM'] else 25 if label in ['NCI-H460','HCC-2998','SW620'] else 35
 rowstatus={name:[] for name in ['all_reported_assays','calibrated_assays']}
 for r in range(1,s.nrows):
  vals=[s.cell_value(r,j) for j in bases];finite=all(isinstance(x,(int,float)) and math.isfinite(x) for x in vals)
  fresh=statistics.mean(vals) if finite else None;end=s.cell_value(r,col);valid=fresh is not None and fresh>0 and isinstance(end,(int,float)) and math.isfinite(end) and end>=0
  ratio=end/fresh if valid else None;flag=valid and end<0.1*fresh;status='screen_positive' if flag else 'screen_negative' if valid else 'unresolved'
  cal=s.cell_value(r,2)==1
  d={'core_label':label,'replicate':rep,'context_id':c['context_id'],'origin_group':c['origin_group'],'partition':c['partition'],'source_excel_row':r+1,'source_excel_column':col+1,'method':s.cell_value(r,0),'metabolite':s.cell_value(r,1),'calibrated':cal,'baseline_mean':fresh,'endpoint_value':end,'endpoint_fresh_ratio':ratio,'status':status,'value_interpretation':'uM per source calibrated flag' if cal else 'uncalibrated signal; relative ratio only'}
  assay.append(d);rowstatus['all_reported_assays'].append(d)
  if cal:rowstatus['calibrated_assays'].append(d)
 for name,rows in rowstatus.items():
  flagged=[x for x in rows if x['status']=='screen_positive'];unknown=[x for x in rows if x['status']=='unresolved']
  status='screen_positive' if flagged else 'unresolved' if unknown else 'screen_negative'
  cult.append({'variant':name,'core_label':label,'replicate':rep,'context_id':c['context_id'],'origin_group':c['origin_group'],'partition':c['partition'],'culture_volume_mL':volume,'duration_days':'4_or_5_not_resolved_per_culture','growth_time_series_available':False,'assays_screened':len(rows),'flagged_assay_count':len(flagged),'unresolved_assay_count':len(unknown),'status':status,'flagged_assays':'; '.join(x['metabolite'] for x in flagged),'in_human1_11_line_subset':label in human})
for name in ['all_reported_assays','calibrated_assays']:
 for label in seen:
  rr=[x for x in cult if x['variant']==name and x['core_label']==label];assert len(rr)==2
  status='screen_positive' if any(x['status']=='screen_positive' for x in rr) else 'unresolved' if any(x['status']=='unresolved' for x in rr) else 'screen_negative'
  profiles.append({'variant':name,'core_label':label,'context_id':rr[0]['context_id'],'origin_group':rr[0]['origin_group'],'partition':rr[0]['partition'],'replicate1_status':rr[0]['status'],'replicate2_status':rr[1]['status'],'both_replicates_screen_negative':status=='screen_negative','status':status,'in_human1_11_line_subset':label in human})
 def count(rows):return {k:sum(x['status']==k for x in rows) for k in ['screen_positive','screen_negative','unresolved']}
 all_c=[x for x in cult if x['variant']==name];all_p=[x for x in profiles if x['variant']==name];negative=set(x['core_label'] for x in all_p if x['status']=='screen_negative')
 counts[name]={'source_assays_screened_per_culture':all_c[0]['assays_screened'],'source_cultures':count(all_c),'source_profiles':count(all_p),'matches_Nilsson_84_positive_36_negative':count(all_c)=={'screen_positive':84,'screen_negative':36,'unresolved':0},'profiles_negative_not_Human1':sorted(negative-human),'Human1_not_profiles_negative':sorted(human-negative),'by_allocation':{part:{'cultures':count([x for x in all_c if x['partition']==part]),'profiles':count([x for x in all_p if x['partition']==part]),'screen_negative_origins':len(set(x['origin_group'] for x in all_p if x['partition']==part and x['status']=='screen_negative'))} for part in sorted(set(x['partition'] for x in all_p))}}
for fn,rows in [('culture_assay_ledger.tsv',assay),('culture_quality_ledger.tsv',cult),('profile_quality_ledger.tsv',profiles)]:
 with (OUT/fn).open('w') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
# Independently recompute the screening counts by metabolite-first loops.
for name in counts:
 flags=[False]*120;unresolved=[False]*120
 for r in range(1,s.nrows):
  if name=='calibrated_assays' and s.cell_value(r,2)!=1:continue
  vals=[s.cell_value(r,i) for i in bases]
  fresh=sum(vals)/len(bases) if all(isinstance(v,(int,float)) and math.isfinite(v) for v in vals) else None
  for k,col in enumerate(cultures):
   val=s.cell_value(r,col)
   if fresh is None or fresh<=0 or not isinstance(val,(int,float)) or not math.isfinite(val) or val<0:unresolved[k]=True
   elif val<fresh/10:flags[k]=True
 check={'screen_positive':sum(flags),'screen_negative':sum(not f and not u for f,u in zip(flags,unresolved)),'unresolved':sum(not f and u for f,u in zip(flags,unresolved))}
 assert check==counts[name]['source_cultures']
assert len(assay)==16800 and len(cult)==240 and len(profiles)==120
report={'executed_at_utc':datetime.now(timezone.utc).isoformat(),'status':'Post-analysis reconstruction of the published concentration screen, not a physiological quality certificate','source_cultures':120,'source_profiles':60,'matched_profiles':sum(c['partition'] in ['development','test'] for c in contexts),'rule':'Any endpoint < 0.1 times mean of ten fresh concentrations; profile rule requires both cultures screen-negative','variants':counts,'independent_count_check':'passed','input_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [RAW,S/'manifests/contexts.tsv',human_path,B/'AUDIT_AMENDMENT.json']},'limits':['Original Nilsson Supplemental Table S1 has not been retrieved; count/label comparisons do not verify its unpublished filtering details','No per-culture growth trajectories, exact duration or culture-size normalization correction recovered','A screen-negative result does not establish biological validity','No model outcome performance evaluated and no original results changed']}
(OUT/'coverage_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['executed_at_utc','matched_profiles','variants','independent_count_check']},indent=2))
