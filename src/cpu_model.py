"""CPU model interface. No outcome values are read here."""
from pathlib import Path
import hashlib
import json
import numpy as np
import cobra
from src.medium_v2 import apply_medium


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_config(root):
    return json.loads((Path(root)/'configs/study.json').read_text())


def build_model(root, config, scenario):
    model = cobra.io.load_json_model(str(Path(root)/'data/raw/Recon3D.json'))
    model.solver = config['solver']
    model.solver.configuration.tolerances.feasibility = config['solver_feasibility_tolerance']
    # The pinned optlang GLPK interface exposes primal tolerance but not tol_dj.
    # Set the GLPK simplex dual tolerance on its documented smcp structure.
    model.solver.configuration._smcp.tol_dj = config['solver_optimality_tolerance']
    model.solver.configuration.timeout = config['solver_timeout_seconds']
    sc = config['scenarios'][scenario]
    opened, missing, changes = apply_medium(model, {'assumed_serum': sc['serum_uptake']})
    model.objective = config['objective']
    gmax = model.slim_optimize()
    if model.solver.status != 'optimal' or not np.isfinite(gmax) or gmax <= 0:
        raise RuntimeError(f'Invalid task maximum: {gmax}, {model.solver.status}')
    task = model.reactions.get_by_id(config['objective'])
    g0 = float(gmax) * sc['task_fraction']
    task.bounds = (g0, g0)
    return model, {'g_max': float(gmax), 'task_flux': g0, 'objective': task.id,
                   'objective_name': task.name, 'scenario': scenario,
                   'opened': opened, 'missing': missing, 'bound_changes': changes}


def score_weights(model, reaction_ids, evidence, k, arm, epsilon):
    """Return costs indexed by reaction; boundary costs remain zero."""
    z = np.asarray(evidence, dtype=float)/(np.asarray(evidence, dtype=float)+k)
    if not np.isfinite(z).all() or (z < 0).any() or (z >= 1).any():
        raise ValueError('Invalid finite nonnegative reaction evidence')
    if arm == 'uniform_pfba':
        q = np.zeros_like(z)
    elif arm == 'ordinal':
        from scipy.stats import rankdata
        q = (rankdata(z, method='average')-0.5)/len(z)
    elif arm.startswith('magnitude_g'):
        q = z ** float(arm.removeprefix('magnitude_g'))
    else:
        raise ValueError(arm)
    supported = epsilon+1-q
    supported /= supported.mean()
    w = {r.id: (0.0 if r.boundary else 1.0) for r in model.reactions}
    for rid, val in zip(reaction_ids, supported):
        if rid not in model.reactions or model.reactions.get_by_id(rid).boundary:
            raise ValueError(f'Invalid supported internal reaction: {rid}')
        w[rid] = float(val)
    return w


def set_absolute_cost(model, weights):
    model.objective = model.problem.Objective(0, direction='min')
    model.objective.set_linear_coefficients({v: float(weights[r.id])
        for r in model.reactions for v in (r.forward_variable, r.reverse_variable)
        if weights.get(r.id, 0) != 0})


def audit_vector(model, flux, costs, primary_cap, secondary_costs, secondary_cap):
    ids = [r.id for r in model.reactions]
    flux = np.asarray(flux, dtype=float)
    if not np.isfinite(flux).all():
        raise ValueError('Nonfinite flux')
    balance = {m.id: 0.0 for m in model.metabolites}
    worst_bound = 0.0
    for r, v in zip(model.reactions, flux):
        for met, coef in r.metabolites.items(): balance[met.id] += coef*v
        worst_bound = max(worst_bound, r.lower_bound-v, v-r.upper_bound)
    cost = float(sum(costs[rid]*abs(v) for rid, v in zip(ids, flux)))
    secondary = float(sum(secondary_costs[rid]*abs(v) for rid,v in zip(ids,flux)))
    return {'finite': True, 'max_mass_balance_residual': max(abs(v) for v in balance.values()),
        'max_bound_violation': max(0.0, float(worst_bound)),
        'net_primary_cost': cost, 'primary_cap_violation': max(0.0, cost-primary_cap),
        'net_secondary_cost': secondary, 'secondary_cap_violation': max(0.0, secondary-secondary_cap)}
