#!/usr/bin/env python3
"""Second reconstruction: targets constant across profiles and their relation to the uniform solution.

    python3 second/summarize.py --root SECOND_ROOT --partition test --predictions results/reserved_v3/predictions.tsv --output results/second_summary.json

Reads the predictions of the arms, the profile table and the target list of the locked analysis plan.
No measurement is read. Fixed in the plan of the second reconstruction, for the 52 primary targets and
the primary scenario, with the point tolerance of the locked plan (1e-5):

  constant under an encoding   the largest minus the smallest point of the target over the profiles of
                               the partition is at most the tolerance;
  at the uniform value         the target's point differs from the point of the uniform cost by at most
                               the tolerance, in every profile;
  moved by the transcript cost for one profile and one encoding, the number of targets whose point
                               differs from the point of the uniform cost by more than the tolerance.
"""
import argparse, csv, json, statistics
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument('--root', type=Path, required=True); ap.add_argument('--partition', choices=['development', 'test'], required=True)
ap.add_argument('--predictions', required=True); ap.add_argument('--output', required=True); ap.add_argument('--scenario', default='primary')
a = ap.parse_args(); root = a.root.resolve()
plan = json.loads((root / 'protocol/analysis_plan.json').read_text()); tol = plan['point_tolerance']
targets = sorted(t['exchange_id'] for t in plan['primary_targets'])
contexts = sorted(r['context_id'] for r in csv.DictReader(open(root / 'manifests/contexts.tsv'), delimiter='\t') if r['partition'] == a.partition)
arms = plan['comparison_arms']
P, U = {}, {}
for r in csv.DictReader(open(root / a.predictions), delimiter='\t'):
    if r['scenario'] != a.scenario or r['exchange_id'] not in targets or r['status'] != 'optimal': continue
    if r['arm'] == 'uniform_pfba': U[r['exchange_id']] = float(r['point'])
    else: P[(r['context_id'], r['arm'], r['exchange_id'])] = float(r['point'])
missing = [(c, m, t) for c in contexts for m in arms for t in targets if (c, m, t) not in P]
if missing or len(U) != len(targets):
    raise SystemExit(f'predictions incomplete: {len(missing)} cells of the encodings missing, {len(targets) - len(U)} of the uniform cost')
const = {m: [t for t in targets if max(P[(c, m, t)] for c in contexts) - min(P[(c, m, t)] for c in contexts) <= tol] for m in arms}
const_all = [t for t in targets if all(t in const[m] for m in arms)]
at_uniform = {m: [t for t in const[m] if all(abs(P[(c, m, t)] - U[t]) <= tol for c in contexts)] for m in arms}
at_uniform_all = [t for t in const_all if all(abs(P[(c, m, t)] - U[t]) <= tol for c in contexts for m in arms)]
moved = {m: [sum(abs(P[(c, m, t)] - U[t]) > tol for t in targets) for c in contexts] for m in arms}
out = dict(note='second reconstruction: targets constant across profiles and their relation to the uniform solution; no measurement read',
           partition=a.partition, scenario=a.scenario, profiles=len(contexts), targets=len(targets), tolerance=tol,
           constant_across_profiles={m: len(const[m]) for m in arms}, constant_under_every_encoding=len(const_all),
           constant_and_at_uniform_value={m: len(at_uniform[m]) for m in arms}, constant_under_every_encoding_and_at_uniform_value=len(at_uniform_all),
           constant_under_every_encoding_not_at_uniform_value=sorted(set(const_all) - set(at_uniform_all)),
           targets_moved_from_uniform_value_per_profile={m: dict(median=statistics.median(v), min=min(v), max=max(v)) for m, v in moved.items()},
           varying_targets_under_every_encoding=sorted(set(targets) - set().union(*[set(const[m]) for m in arms])))
Path(root / a.output).write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: v for k, v in out.items() if k not in ('note', 'varying_targets_under_every_encoding')}, indent=1))
