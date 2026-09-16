"""Compare only after the isolated arm has been produced and saved."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np

report = Path(__file__).resolve().parent
root = report.parents[1]
isolated = report / 'isolated'
fresh_path = isolated / 'fresh_arm.json'
if not fresh_path.is_file():
    raise RuntimeError('Isolated result must exist before reading any main result')
fresh = json.loads(fresh_path.read_text())
if fresh['status'] != 'optimal' or not fresh.get('completed_at'):
    raise RuntimeError('Isolated result is incomplete')
comparison_started = datetime.now(timezone.utc).isoformat()
context, arm, scenario = 'CO:HT29', 'magnitude_g1.00', 'primary'
token = hashlib.sha256(f'{context}|{arm}|{scenario}'.encode()).hexdigest()[:16]
main_path = root / 'results/development_v3/arms' / f'{token}.json'
main = json.loads(main_path.read_text())  # First main-result read occurs here.
assert main['status'] == 'optimal'
assert (main['context_id'], main['arm'], main['scenario']) == (context, arm, scenario)
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
fresh_npz = isolated / fresh['flux_file']
main_npz = root / 'results/development_v3' / main['flux_file']
assert digest(fresh_npz) == fresh['flux_sha256']
assert digest(main_npz) == main['flux_sha256']
with np.load(fresh_npz, allow_pickle=False) as f, np.load(main_npz, allow_pickle=False) as m:
    ids = f['reaction_ids'].astype(str).tolist()
    assert ids == m['reaction_ids'].astype(str).tolist()
    flux = f['flux'].copy()
    costs = f['reaction_costs'].copy()
    flux_difference = float(np.max(np.abs(flux-m['flux'])))
    costs_difference = float(np.max(np.abs(costs-m['reaction_costs'])))
panel = sorted(fresh['point'])
assert panel == sorted(main['point']) and len(panel) == 96
point_difference = max(abs(fresh['point'][rid]-main['point'][rid]) for rid in panel)
range_difference = max(abs(a-b) for rid in panel for a,b in zip(fresh['ranges'][rid],main['ranges'][rid]))
objective_fields = ['primary_optimum','primary_cap','secondary_optimum','secondary_cap']
objective_differences = {key: abs(fresh[key]-main[key]) for key in objective_fields}
metadata_difference = {key: abs(fresh['model'][key]-main['model'][key]) for key in ('g_max','task_flux')}
copy_manifest = json.loads((report/'copy_manifest.json').read_text())
copied_files = []
for row in copy_manifest['files']:
    relative = row['relative_path']
    copied_files.append({**row, 'isolated_unchanged': digest(isolated/relative)==row['sha256'],
                         'matches_current_main_file': digest(root/relative)==row['sha256']})
paths = ['configs/study.json','data/raw/Recon3D.json','data/scores_conservative.npz','manifests/contexts.tsv',
         'manifests/metabolite_map.tsv','src/cpu_model.py','src/cpu_readout.py','src/run_study.py','src/medium_v2.py','environment.lock']
main_inputs = {path:digest(root/path) for path in paths}
main_fingerprint = hashlib.sha256(json.dumps(main_inputs,sort_keys=True).encode()).hexdigest()
source_matches = main_fingerprint == main['input_fingerprint']
# Independent algebraic audit reads the isolated model/evidence, without another solve.
sys.path.insert(0, str(root))
from src.audit_results import reconstruct_geometry, load_recipe, reconstruct_costs, audit_arm
raw = json.loads((isolated/'data/raw/Recon3D.json').read_text())
geometry = reconstruct_geometry(raw, fresh['config'], scenario, load_recipe(isolated/'src/medium_v2.py'))
with np.load(isolated/'data/scores_conservative.npz',allow_pickle=False) as data:
    supported = data['rxn'].astype(str).tolist()
    names = data['contexts'].astype(str).tolist()
    evidence = data['A'][:,names.index(context)].copy()
expected = reconstruct_costs(geometry, supported, evidence, fresh['scale_k'], arm, fresh['config']['epsilon'])
independent = audit_arm(fresh, flux, ids, costs, geometry, expected, fresh['config'], panel)
tolerance = fresh['config']['validation_tolerance']
all_differences = [flux_difference,costs_difference,point_difference,range_difference,*objective_differences.values(),*metadata_difference.values()]
passed = (independent['pass'] and source_matches and max(all_differences)<=tolerance and
          all(row['isolated_unchanged'] and row['matches_current_main_file'] for row in copied_files))
comparison = {'status':'passed' if passed else 'failed','pass':passed,
    'context_id':context,'arm':arm,'scenario':scenario,'reported_exchanges':len(panel),'full_vector_reactions':len(ids),
    'fresh_completed_at':fresh['completed_at'],'comparison_started_at':comparison_started,
    'fresh_elapsed_seconds':fresh['elapsed_seconds'],'main_result_path':str(main_path),
    'fresh_json_sha256':digest(fresh_path),'main_json_sha256':digest(main_path),
    'fresh_npz_sha256':digest(fresh_npz),'main_npz_sha256':digest(main_npz),'npz_byte_identical':digest(fresh_npz)==digest(main_npz),
    'max_abs_full_flux_difference':flux_difference,'max_abs_cost_weight_difference':costs_difference,
    'max_abs_exchange_point_difference':point_difference,'max_abs_exchange_range_difference':range_difference,
    'objective_differences':objective_differences,'model_metadata_differences':metadata_difference,
    'comparison_tolerance':tolerance,'main_source_fingerprint_current':source_matches,
    'copied_files':copied_files,'independent_fresh_vector_audit':independent,
    'shared_python_environment':copy_manifest['shared_python_environment'],
    'scope':'One actual fresh-directory arm using independently rebuilt data and copied final model sources. Shared Python environment; not a new-environment install or whole-study clean rerun. No CORE performance evaluated.'}
(report/'comparison.json').write_text(json.dumps(comparison,indent=2,allow_nan=False)+'\n')
summary = f'''# Fresh-directory single-arm reproduction

Status: **{'PASS' if passed else 'FAIL'}**. Run on 13 September 2026.

The isolated run reproduced **CO:HT29 / magnitude_g1.00 / primary**, reporting **{len(panel)} exchanges** and a full **{len(ids):,}-reaction flux vector**. It used the fresh data build under `work/reproduction_clean`, with only the final `cpu_model.py`, `cpu_readout.py`, `medium_v2.py`, configuration and environment lock copied as production implementation inputs. A reviewer-written driver invoked the model/readout directly; it did not invoke the full-study runner.

The isolated result was saved at **{fresh['completed_at']}**. The main result was first read after this, at **{comparison_started}**. The isolated calculation took **{fresh['elapsed_seconds']:.2f} seconds**. No CORE outcome comparison was performed.

| Comparison | Maximum absolute difference |
|---|---:|
| Full feasible flux vector | {flux_difference:.12g} |
| Reaction cost weights | {costs_difference:.12g} |
| 96 reported exchange points | {point_difference:.12g} |
| Exchange range endpoints | {range_difference:.12g} |
| Primary/secondary optima and caps | {max(objective_differences.values()):.12g} |
| Saved task maximum and fixed demand | {max(metadata_difference.values()):.12g} |

All copied input/source files {'match' if all(r['matches_current_main_file'] for r in copied_files) else 'do not match'} the current main files. The main arm's source fingerprint is {'current' if source_matches else 'mismatched'}. The compressed flux archives are {'byte-identical' if comparison['npz_byte_identical'] else 'not byte-identical'}.

An independent saved-vector audit {'passed' if independent['pass'] else 'failed'} mass balance, reconstructed bounds, fixed task, both cost caps, point/range correspondence, independently reconstructed weights and final logical LP statuses. Its detailed residuals and any recovered attempts are in `comparison.json`.

**Scope:** both runs used the same existing `.venv` Python environment. This checks one fresh-directory model execution with fresh-derived inputs; it is not a separate environment installation, full cohort reproduction or biological validation. Objective/FVA extrema were compared between two solves, while the algebraic audit itself did not independently optimize them.

Evidence: `copy_manifest.json`, `isolated/fresh_driver.py`, `isolated/fresh_execution.log`, `isolated/fresh_arm.json`, `isolated/fresh_arm.npz`, `compare_after_production.py`, and `comparison.json`.
'''
(report/'README.md').write_text(summary)
manifest = []
for path in sorted(report.rglob('*')):
    if path.is_file() and '__pycache__' not in path.parts and path.name != 'artifact_hashes.json':
        manifest.append({'path':str(path.relative_to(report)),'sha256':digest(path),'bytes':path.stat().st_size})
(report/'artifact_hashes.json').write_text(json.dumps({'files':manifest},indent=2)+'\n')
print(json.dumps({key:comparison[key] for key in ['status','full_vector_reactions','max_abs_full_flux_difference','max_abs_exchange_point_difference','max_abs_exchange_range_difference','npz_byte_identical','fresh_elapsed_seconds']},indent=2))
raise SystemExit(0 if passed else 1)
