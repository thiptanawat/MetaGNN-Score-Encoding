#!/usr/bin/env python3
"""Second reconstruction: the reported quantities without the profiles that the second solver completed.

    python3 second/without_second_solver.py --root SECOND_ROOT --check-only
    python3 second/without_second_solver.py --root SECOND_ROOT --output results/without_second_solver.json

Plan: second/WITHOUT_SECOND_SOLVER_PLAN_2026-11-15.md, fixed with this script before the second command was
run. In the reserved run of the second reconstruction the first solver stopped on an infeasible problem in 4
of the 283 arms, and the second solver completed them (results/reserved_v3/second_solver_arms.json). The points
of those 4 arms entered every comparison beside the points of the first solver. This script removes the 4
profiles that own such an arm, under every encoding, and computes the reported quantities again from the saved
records. No linear program is solved and no file of the reserved run is changed.

Recomputed for the remaining profiles, with the functions of the main analysis where they exist:
  encoding     src.evaluation.evaluate_panel: concordance of every encoding on the 52 primary targets, the two
               contrasts and their 95% bootstrap intervals (the seed and the number of resamples of the locked plan);
               src.transform_sensitivity.summarize: predictions changed by the encoding and orderings reversed;
  attribution  the rows of results/mechanism/reserved_range/summary.json: coordinates the network leaves open and
               coordinates first fixed by the primary cost cap; the rule of second/summarize.py: targets constant
               across profiles and how many of them take the value of the uniform solution;
  allowance    the jobs of results/closure/cap_sweep: mean number of fixed coordinates under each cap; the rule of
               closure/cap_concordance.py: orderings certified by the interval rule under each cap.
The readings of Section 6 of second/SECOND_RECONSTRUCTION_PLAN_2026-11-14.md are applied to the remaining
profiles as they were applied to all of them.

Before anything is removed, the same code computes every quantity for all profiles and compares it with the
deposited value; the script stops if one differs. --check-only does this and nothing else. The agreement of the
two solvers on the arms that both solved (step F of that plan) is not recomputed: it removes no profile.
"""
import argparse, csv, glob, json, statistics, sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--root', type=Path, required=True)
ap.add_argument('--output', default='results/without_second_solver.json', help='relative to the root')
ap.add_argument('--check-only', action='store_true', help='reproduce the deposited values for all profiles and stop; nothing is removed and nothing is written')
A = ap.parse_args(); ROOT = A.root.resolve()
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'closure')); sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_arms import verify_second_lock
import src.evaluation as E
import src.transform_sensitivity as TS
from stat_repairs import load as sr_load, build_pairs as sr_build_pairs, TIE

LOCK = verify_second_lock(ROOT)
PLAN = json.loads((ROOT / 'protocol/analysis_plan.json').read_text())
for _rel, _digest in PLAN['input_sha256'].items():
    if E.sha256(ROOT / _rel) != _digest: sys.exit(f'input differs from the analysis plan: {_rel}')
PRED = ROOT / 'results/reserved_v3/predictions.tsv'; SCEN = 'primary'
ENC = PLAN['comparison_arms']; IDENT, MID = 'magnitude_g1.00', 'ordinal'; CAPS = ['exact', 'rel_1pct', 'rel_5pct', 'rel_20pct']
ALL = E.context_table(E.read_tsv(ROOT / 'manifests/contexts.tsv'), 'test')
J = lambda rel: json.loads((ROOT / rel).read_text())


def compute(keep):
    """Every reported quantity of the second reconstruction for the profiles in `keep`."""
    ctx = dict(sorted((c, r) for c, r in ALL.items() if c in keep)); out = dict(profiles=len(ctx), origins=len({r['origin_group'] for r in ctx.values()}))
    # encoding: the locked evaluation on the kept profiles
    ev = E.evaluate_panel(E.prediction_index(E.read_tsv(PRED), SCEN), E.mean_observations(E.read_tsv(ROOT / 'data/outcomes_long.tsv'), ctx), ctx,
                          PLAN['primary_targets'], PLAN['assay_tolerances'], PLAN, resamples=PLAN['bootstrap_resamples'], seed=PLAN['test_seed'])
    b = ev['bootstrap'] or {}
    out['evaluation'] = dict(status=ev['status'], eligible_pairs=ev['observed_eligible_pairs'], zero_pair_targets=ev['zero_observed_pair_targets'],
                             primary=ev['primary'], delta_grid_interval=b.get('delta_grid_interval'), delta_identity_interval=b.get('delta_identity_interval'),
                             valid_draws=b.get('valid_draws'), invalid_draws=b.get('invalid_draws'), origin_groups=b.get('origin_groups'),
                             seed=b.get('seed'), resamples=b.get('resamples'), interval_status=b.get('interval_status'))
    tr = TS.summarize(PRED, SCEN, ctx, PLAN['primary_targets'], PLAN['point_tolerance'])
    out['transform'] = {k: tr[k] for k in ('context_target_cells', 'changed_cells', 'changed_fraction', 'all_distinct_origin_context_target_pairs',
                                           'strict_order_reversals_across_grid', 'strict_order_reversal_fraction', 'tie_transition_without_strict_reversal_pairs')}
    # attribution, part 1: the range stages (rows of the deposited summary, one per profile and encoding)
    arm = defaultdict(lambda: dict(open=0, first_fixed_by_cost_cap=0, fixed_by_network=0, rows=0))
    for r in J('results/mechanism/reserved_range/summary.json')['summary_by_profile_arm']:
        if r['context'] not in ctx: continue
        a = arm[r['arm']]; a['open'] += r['targets'] - r['B0_fixed']; a['first_fixed_by_cost_cap'] += r['first_fixed_B1']; a['fixed_by_network'] += r['B0_fixed']; a['rows'] += 1
    for a in arm.values(): a['share_first_fixed_by_cost_cap'] = a['first_fixed_by_cost_cap'] / a['open']
    out['range_stages'] = {k: arm[k] for k in sorted(arm)}
    # attribution, part 2: targets constant across profiles (the rule of second/summarize.py)
    tol = PLAN['point_tolerance']; targets = sorted(t['exchange_id'] for t in PLAN['primary_targets']); P, U = {}, {}
    for r in csv.DictReader(open(PRED), delimiter='\t'):
        if r['scenario'] != SCEN or r['exchange_id'] not in targets or r['status'] != 'optimal': continue
        if r['arm'] == 'uniform_pfba': U[r['exchange_id']] = float(r['point'])
        elif r['context_id'] in ctx: P[(r['context_id'], r['arm'], r['exchange_id'])] = float(r['point'])
    if any((c, m, t) not in P for c in ctx for m in ENC for t in targets) or len(U) != len(targets): sys.exit('predictions incomplete for the kept profiles')
    const = {m: [t for t in targets if max(P[(c, m, t)] for c in ctx) - min(P[(c, m, t)] for c in ctx) <= tol] for m in ENC}
    const_all = [t for t in targets if all(t in const[m] for m in ENC)]
    at_uni = [t for t in const_all if all(abs(P[(c, m, t)] - U[t]) <= tol for c in ctx for m in ENC)]
    moved = {m: [sum(abs(P[(c, m, t)] - U[t]) > tol for t in targets) for c in ctx] for m in ENC}
    out['constant_targets'] = dict(targets=len(targets), constant_under_every_encoding=len(const_all), constant_under_every_encoding_and_at_uniform_value=len(at_uni),
                                   targets_moved_from_uniform_value_per_profile={m: dict(median=statistics.median(v), min=min(v), max=max(v)) for m, v in moved.items()})
    # allowance: the cap sweep (identity and midrank encodings, one job per profile) and the interval rule under each cap
    jobs = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(ROOT / 'results/closure/cap_sweep/jobs/*.json')))]
    jobs = [d for d in jobs if d['status'] == 'complete']
    test, origin, tg, atol, obs, obs_mean, PP = sr_load(ROOT, SCEN); ex_of = dict(tg)
    arms7 = ['magnitude_g0.50', 'magnitude_g0.80', 'magnitude_g1.00', 'magnitude_g1.25', 'magnitude_g2.00', 'ordinal', 'uniform_pfba']
    ref = {ex: v for (c, ex), v in PP['uniform_pfba'].items() if c == '__reference__'}; PP['uniform_pfba'] = {(c, ex): v for c in test for ex, v in ref.items()}
    kept = [c for c in test if c in ctx]; pairs = sr_build_pairs(kept, origin, tg, atol, obs_mean, PP, arms7); n_pairs = int(sum(len(v) for v in pairs.values()))
    if n_pairs != ev['observed_eligible_pairs']: sys.exit('the two pair builders disagree on the number of eligible pairs')
    cap = {}
    for a in (IDENT, MID):
        rr = [d for d in jobs if d['arm'] == a and d['context'] in ctx]; cap[a] = dict(profiles=len(rr))
        if len(rr) != len(ctx): sys.exit(f'cap sweep: {len(rr)} complete jobs for {a}, {len(ctx)} profiles kept')
        for name in CAPS:
            R = {(d['context'], ex): v for d in rr for ex, v in d['caps'][name]['ranges'].items()}; res = 0
            for m, arr in pairs.items():
                for x, y, truth in arr:
                    p, q = R[(kept[int(x)], ex_of[m])], R[(kept[int(y)], ex_of[m])]
                    res += int(p[0] > q[1] + TIE or q[0] > p[1] + TIE)
            cap[a][name] = dict(mean_fixed_of_52=float(np.mean([d['caps'][name]['fixed'] for d in rr])), resolved_pairs=res, resolved_share=res / n_pairs)
    out['cap_sweep'] = cap
    # the readings of Section 6 of the plan of the second reconstruction
    e = out['evaluation']; inc = lambda iv: iv is not None and iv[0] < 0 < iv[1]
    out['readings'] = dict(
        encoding=bool(e['status'] == 'complete_predictions' and inc(e['delta_grid_interval']) and inc(e['delta_identity_interval'])),
        attribution=bool(all(a['first_fixed_by_cost_cap'] > 0.5 * a['open'] for a in out['range_stages'].values())
                         and out['constant_targets']['constant_under_every_encoding_and_at_uniform_value'] > 0.5 * out['constant_targets']['constant_under_every_encoding']),
        allowance=bool(all(cap[a]['rel_1pct']['resolved_share'] < cap[a]['exact']['resolved_share'] and cap[a]['rel_1pct']['mean_fixed_of_52'] < cap[a]['exact']['mean_fixed_of_52']
                           for a in (IDENT, MID))))
    return out


def same(a, b): return a == b or (isinstance(a, float) and isinstance(b, float) and abs(a - b) <= 1e-12 * max(1.0, abs(a), abs(b)))


def check_against_deposit(full):
    """The all-profile values of this script against the deposited files of the reserved run."""
    dep = J('results/reserved_v3/evaluation_primary.json')['primary_fixed_panel']; e = full['evaluation']; ok = {}
    ok['evaluation'] = (e['eligible_pairs'] == dep['observed_eligible_pairs'] and all(same(e['primary'][k], dep['primary'][k]) for k in ('delta_grid', 'delta_identity'))
                        and all(same(e['primary']['absolute_macro_C'][k], v) for k, v in dep['primary']['absolute_macro_C'].items())
                        and all(same(x, y) for k in ('delta_grid_interval', 'delta_identity_interval') for x, y in zip(e[k], dep['bootstrap'][k]))
                        and e['valid_draws'] == dep['bootstrap']['valid_draws'])
    dt = J('results/reserved_v3/transform_primary.json')['primary_fixed_panel']; ok['transform'] = all(same(full['transform'][k], dt[k]) for k in full['transform'])
    ds = J('results/second_summary.json'); c = full['constant_targets']
    ok['constant_targets'] = (all(c[k] == ds[k] for k in ('constant_under_every_encoding', 'constant_under_every_encoding_and_at_uniform_value'))
                              and c['targets_moved_from_uniform_value_per_profile'] == ds['targets_moved_from_uniform_value_per_profile'])
    dr = J('results/mechanism/reserved_range/summary.json'); ok['range_stages'] = sum(a['rows'] for a in full['range_stages'].values()) == len(dr['summary_by_profile_arm']) == dr['B1_completed']
    dc = J('results/closure/cap_sweep/summary.json'); dk = J('results/closure/cap_sweep/cap_concordance.json')
    ok['cap_sweep'] = all(same(full['cap_sweep'][a][n]['mean_fixed_of_52'], dc[a][n]['mean_fixed_of_52']) and full['cap_sweep'][a][n]['resolved_pairs'] == dk[a][n]['resolved_pairs']
                          and same(full['cap_sweep'][a][n]['resolved_share'], dk[a][n]['resolved_share']) for a in (IDENT, MID) for n in CAPS)
    return ok


def main():
    passes = J('results/reserved_v3/second_solver_arms.json')
    second = sorted((x['context_id'], x['arm']) for x in passes['solved_in_pass2']); removed = sorted({c for c, _ in second})
    if any(c not in ALL for c in removed): sys.exit('an arm of the second solver names a profile outside the reserved partition')
    full = compute(set(ALL)); ok = check_against_deposit(full)
    print('all profiles, against the deposited values:', json.dumps(ok))
    if not all(ok.values()): sys.exit('a value for all profiles differs from the deposited value; nothing was removed')
    if A.check_only: return
    less = compute(set(ALL) - set(removed))
    if less['profiles'] != len(ALL) - len(removed): sys.exit('the number of kept profiles is not what the removal gives')
    out = dict(note='second reconstruction: the reported quantities for all reserved profiles and without the profiles that own an arm completed by the second solver; '
                    'computed from the saved records, no linear program solved',
               plan='second/WITHOUT_SECOND_SOLVER_PLAN_2026-11-15.md', second_reconstruction_plan_sha256=LOCK['second_reconstruction_plan_sha256'],
               predictions_sha256=E.sha256(PRED), arms_completed_by_the_second_solver=[dict(context_id=c, arm=a) for c, a in second], profiles_removed=removed,
               all_profiles_reproduce_the_deposited_values=ok, all_profiles=full, without_removed_profiles=less)
    (ROOT / A.output).write_text(json.dumps(out, indent=1, allow_nan=False) + '\n')
    e = less['evaluation']
    print(json.dumps(dict(output=str(ROOT / A.output), profiles=less['profiles'], origins=less['origins'], eligible_pairs=e['eligible_pairs'], primary=e['primary'],
                          delta_grid_interval=e['delta_grid_interval'], delta_identity_interval=e['delta_identity_interval'], readings=less['readings']), indent=1))


if __name__ == '__main__':
    main()
