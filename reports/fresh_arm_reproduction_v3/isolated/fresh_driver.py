"""Reviewer-written single-arm reproduction driver; no main-run result access."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
import numpy as np
import swiglpk
from src import cpu_model, cpu_readout, medium_v2

root = Path(__file__).resolve().parent
started = datetime.now(timezone.utc).isoformat()
clock = time.monotonic()
for module in (cpu_model, cpu_readout, medium_v2):
    if not Path(module.__file__).resolve().is_relative_to(root):
        raise RuntimeError(f"Production module escaped isolated source directory: {module.__file__}")
with (root / 'manifests/contexts.tsv').open() as f:
    metadata = list(csv.DictReader(f, delimiter='\t'))
with (root / 'manifests/metabolite_map.tsv').open() as f:
    mapping = list(csv.DictReader(f, delimiter='\t'))
development = sorted(row['context_id'] for row in metadata if row['partition'] == 'development')
panel = sorted({row['exchange_id'] for row in mapping if row['primary_eligible'].lower() in ('true','1','yes')})
assert len(panel) == 96
with np.load(root / 'data/scores_conservative.npz', allow_pickle=False) as scores:
    reaction_ids = scores['rxn'].astype(str).tolist()
    contexts = scores['contexts'].astype(str).tolist()
    evidence = scores['A'].copy()
development_values = evidence[:, [contexts.index(context) for context in development]]
k = float(np.median(development_values[development_values > 0]))
config = cpu_model.load_config(root)
context, arm, scenario = 'CO:HT29', 'magnitude_g1.00', 'primary'
model, model_metadata = cpu_model.build_model(root, config, scenario)
weights = cpu_model.score_weights(model, reaction_ids, evidence[:, contexts.index(context)], k, arm, config['epsilon'])
print(json.dumps({'started_at': started, 'context': context, 'arm': arm, 'scenario': scenario, 'exchanges': len(panel), 'scale_k': k}), flush=True)
result = cpu_readout.select_flux(model, weights, panel, config)
output = root / 'fresh_arm.npz'
np.savez_compressed(output, flux=result.pop('flux'), reaction_ids=np.array(result.pop('reaction_ids')),
                    reaction_costs=np.array([weights[r.id] for r in model.reactions]))
result.update({'status': 'optimal', 'context_id': context, 'arm': arm, 'scenario': scenario,
               'model': model_metadata, 'config': config, 'scale_k': k, 'development_contexts': development,
               'flux_file': output.name, 'flux_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
               'started_at': started, 'completed_at': datetime.now(timezone.utc).isoformat(),
               'elapsed_seconds': time.monotonic()-clock,
               'python_executable': sys.executable, 'python_version': platform.python_version(),
               'glpk_version': swiglpk.glp_version(),
               'package_versions': {name: importlib.metadata.version(name) for name in ('cobra','optlang','numpy','scipy','swiglpk')},
               'production_module_paths': {m.__name__: str(Path(m.__file__).resolve()) for m in (cpu_model,cpu_readout,medium_v2)},
               'scope': 'Single-arm model reproduction with freshly regenerated data and shared environment; no CORE performance calculation and no main-result access before saving.'})
(root / 'fresh_arm.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps({'status': result['status'], 'elapsed_seconds': result['elapsed_seconds'], 'flux_sha256': result['flux_sha256'], 'completed_at': result['completed_at']}), flush=True)
