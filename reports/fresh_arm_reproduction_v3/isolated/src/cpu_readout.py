"""Outcome-blind sequential LP point selection and pre-selection ranges."""
import time
import numpy as np
from src.cpu_model import set_absolute_cost, audit_vector


def select_flux(model, costs, reported, config, coordinate_factors=None, secondary_costs=None):
    """Canonical coordinate = factor * model flux. Return full jointly feasible flux.

    Ranges are computed before the ordered exchange coordinates are fixed. The
    coordinate factors make the convention invariant to representation changes.
    All constraints added in this call are removed when it returns or fails.
    """
    t0 = time.monotonic()
    factors = {rid: float((coordinate_factors or {}).get(rid, 1.0)) for rid in reported}
    secondary_costs = secondary_costs or {r.id: 1.0 for r in model.reactions}
    report = sorted(reported)
    statuses, added = [], []
    validation_tol = config['validation_tolerance']

    def optimize(stage, rid=''):
        for attempt in range(2):
            start = time.monotonic()
            value = model.slim_optimize(error_value=float('nan'))
            status = model.solver.status
            finite = bool(np.isfinite(value))
            statuses.append({'stage': stage, 'exchange_id': rid, 'attempt': attempt,
                             'status': status, 'objective': float(value) if finite else None,
                             'seconds': time.monotonic()-start})
            if status == 'optimal' and finite:
                return float(value)
            if attempt == 0:
                # Rebuild the numerical basis without changing any LP coefficients,
                # constraints, objective, bounds, or tolerance settings.
                from swiglpk import glp_adv_basis
                glp_adv_basis(model.solver.problem, 0)
        raise RuntimeError(f'{stage} {rid}: {status}, {value}; identical-LP retry exhausted')

    try:
        set_absolute_cost(model, costs)
        primary_opt = optimize('primary')
        primary_expr = model.objective.expression
        primary_cap = primary_opt + config['cost_absolute_allowance'] + abs(primary_opt)*config['cost_relative_allowance']
        cap = model.problem.Constraint(primary_expr, ub=primary_cap, name='cpu_primary_cap')
        model.add_cons_vars([cap]); added.append(cap)
        set_absolute_cost(model, secondary_costs)
        secondary_opt = optimize('secondary')
        secondary_expr = model.objective.expression
        secondary_cap = secondary_opt + config['secondary_absolute_allowance'] + abs(secondary_opt)*config['secondary_relative_allowance']
        cap2 = model.problem.Constraint(secondary_expr, ub=secondary_cap, name='cpu_secondary_cap')
        model.add_cons_vars([cap2]); added.append(cap2)
        model.solver.update()
        ranges = {}
        for rid in report:
            rx = model.reactions.get_by_id(rid)
            expr = factors[rid]*rx.flux_expression
            model.objective = model.problem.Objective(expr, direction='min')
            lo = optimize('range_min', rid)
            model.objective = model.problem.Objective(expr, direction='max')
            hi = optimize('range_max', rid)
            if lo-hi > validation_tol: raise RuntimeError(f'Reversed range {rid}: {lo}, {hi}')
            ranges[rid] = [lo, hi]
        fixes = {}
        for index, rid in enumerate(report):
            expr = factors[rid]*model.reactions.get_by_id(rid).flux_expression
            model.objective = model.problem.Objective(expr, direction='min')
            selected = optimize('coordinate_min', rid)
            delta = config['coordinate_fix_allowance']
            con = model.problem.Constraint(expr, lb=selected-delta, ub=selected+delta, name=f'cpu_coordinate_{index}')
            model.add_cons_vars([con]); added.append(con)
            model.solver.update()
            fixes[rid] = selected
        model.objective = model.problem.Objective(0, direction='min')
        optimize('final_feasibility')
        flux = np.array([r.forward_variable.primal-r.reverse_variable.primal for r in model.reactions], dtype=float)
        point = {rid: float(factors[rid]*flux[model.reactions.index(rid)]) for rid in report}
        audit = audit_vector(model, flux, costs, primary_cap, secondary_costs, secondary_cap)
        audit['max_coordinate_fix_violation'] = max([max(0,abs(point[rid]-fixes[rid])-config['coordinate_fix_allowance']) for rid in report],default=0)
        audit['max_point_range_violation'] = max([max(0,ranges[rid][0]-point[rid],point[rid]-ranges[rid][1]) for rid in report],default=0)
        failures = {k:v for k,v in audit.items() if (k.endswith('violation') or k.endswith('residual')) and v > validation_tol}
        if failures: raise RuntimeError(f'Final flux audit failed: {failures}')
        return {'point':point, 'ranges':ranges, 'flux':flux, 'reaction_ids':[r.id for r in model.reactions],
            'primary_optimum':primary_opt,'primary_cap':primary_cap,'secondary_optimum':secondary_opt,
            'secondary_cap':secondary_cap,'audit':audit,'solve_ledger':statuses,
            'seconds':time.monotonic()-t0, 'coordinate_order':report}
    finally:
        if added: model.remove_cons_vars(added); model.solver.update()
