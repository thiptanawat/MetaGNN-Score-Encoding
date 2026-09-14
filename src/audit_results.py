"""Independently audit saved CPU flux vectors without solving LPs or reading CORE values.

This verifies feasibility, encoding, ledger consistency and artifact identity.
It cannot independently certify optimality, g_max or FVA extrema without solving.
The CLI writes JSON to stdout only, with a nonzero exit status on any failure.
"""
from __future__ import annotations

import argparse
import ast
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from .verify_inputs import RAW_SOURCES, VerificationError, check_hash_mapping, digest, inside, read_json

FINGERPRINT_FILES = (
    "configs/study.json", "data/raw/Recon3D.json", "data/scores_conservative.npz",
    "manifests/contexts.tsv", "manifests/metabolite_map.tsv", "src/cpu_model.py",
    "src/cpu_readout.py", "src/run_study.py", "src/medium_v2.py", "environment.lock",
)


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def truthy(value):
    return str(value).lower() in ("true", "1", "yes")


def load_recipe(path):
    """Read literal recipe data, not the production medium-application function."""
    needed = {"BASAL", "CULTURE", "SERUM", "UPTAKE"}
    recipe = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in needed:
                    recipe[target.id] = ast.literal_eval(node.value)
    if set(recipe) != needed:
        raise VerificationError("Cannot independently read the four literal medium recipe definitions")
    return recipe


@dataclass
class Geometry:
    ids: list
    matrix: object
    lower: np.ndarray
    upper: np.ndarray
    boundary: np.ndarray
    task_index: int


def reconstruct_geometry(raw, config, scenario, recipe):
    """Build sparse S and bounds directly from pinned JSON, without COBRA/solver calls."""
    reactions = raw["reactions"]
    ids = [r["id"] for r in reactions]
    if len(ids) != len(set(ids)):
        raise VerificationError("Duplicate model reaction identifiers")
    metabolites = [m["id"] for m in raw["metabolites"]]
    if len(metabolites) != len(set(metabolites)):
        raise VerificationError("Duplicate model metabolite identifiers")
    mi = {name: i for i, name in enumerate(metabolites)}
    rr, cc, values = [], [], []
    for j, reaction in enumerate(reactions):
        for name, value in reaction["metabolites"].items():
            rr.append(mi[name]); cc.append(j); values.append(float(value))
    matrix = coo_matrix((values, (rr, cc)), shape=(len(metabolites), len(ids))).tocsr()
    lower = np.array([r["lower_bound"] for r in reactions], dtype=float)
    upper = np.array([r["upper_bound"] for r in reactions], dtype=float)
    for j, rid in enumerate(ids):
        if rid.startswith(("DM_", "SK_")):
            lower[j] = max(0, lower[j])
        elif rid.startswith("EX_"):
            lower[j] = 0
            upper[j] = max(1000, upper[j])
    index = {rid: j for j, rid in enumerate(ids)}
    capacities = dict(recipe["UPTAKE"])
    capacities["assumed_serum"] = config["scenarios"][scenario]["serum_uptake"]
    for stem, category in recipe["BASAL"] + recipe["CULTURE"] + recipe["SERUM"]:
        rid = f"EX_{stem}_e"
        if category == "unrepresented_source_b12":
            if rid in index:
                raise VerificationError("Unresolved vitamin-B12 placeholder unexpectedly has a model exchange")
            continue
        if rid in index:
            lower[index[rid]] = -float(capacities[category])
    # Match the documented single-metabolite boundary definition, without calling it.
    boundary = np.array([len(r["metabolites"]) == 1 for r in reactions], dtype=bool)
    return Geometry(ids, matrix, lower, upper, boundary, index[config["objective"]])


def average_ranks(values):
    """Independent average-rank implementation, avoiding the producer's rankdata call."""
    values = np.asarray(values)
    order = np.argsort(values, kind="stable")
    result = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        result[order[start:end]] = (start + 1 + end) / 2
        start = end
    return result


def reconstruct_costs(geometry, supported_ids, evidence, k, arm, epsilon):
    values = np.asarray(evidence, dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0) or not np.isfinite(k) or k <= 0:
        raise VerificationError("Invalid reaction evidence or development scale")
    if len(values) != len(supported_ids) or len(set(supported_ids)) != len(supported_ids):
        raise VerificationError("Invalid supported reaction axis")
    z = values / (values + k)
    if arm == "uniform_pfba":
        q = np.zeros_like(z)
    elif arm == "ordinal":
        q = (average_ranks(z) - 0.5) / len(z)
    elif arm.startswith("magnitude_g"):
        q = z ** float(arm[len("magnitude_g"):])
    else:
        raise VerificationError(f"Unknown arm: {arm}")
    supported_cost = epsilon + 1 - q
    supported_cost /= supported_cost.mean()
    expected = np.where(geometry.boundary, 0.0, 1.0)
    index = {rid: j for j, rid in enumerate(geometry.ids)}
    for rid, cost in zip(supported_ids, supported_cost):
        if rid not in index or geometry.boundary[index[rid]]:
            raise VerificationError(f"Supported evidence refers to an invalid internal reaction: {rid}")
        expected[index[rid]] = cost
    return expected


def expected_stages(panel):
    return ([('primary', ''), ('secondary', '')] +
            [(stage, rid) for rid in panel for stage in ('range_min', 'range_max')] +
            [('coordinate_min', rid) for rid in panel] + [('final_feasibility', '')])


def audit_arm(record, flux, reaction_ids, stored_costs, geometry, expected_costs, config, panel):
    """Pure vector/record audit, independently testable on small known-answer models."""
    failures, metrics = [], {}
    def require(condition, name):
        if not condition: failures.append(name)
    tol = float(config["validation_tolerance"])
    flux = np.asarray(flux, dtype=float)
    stored_costs = np.asarray(stored_costs, dtype=float)
    require(record.get("status") == "optimal", "arm_status")
    require(list(reaction_ids) == geometry.ids, "reaction_axis")
    if flux.shape != (len(geometry.ids),) or stored_costs.shape != flux.shape:
        return {"pass": False, "failures": failures + ["array_shape"], "metrics": metrics}
    if not np.isfinite(flux).all() or not np.isfinite(stored_costs).all():
        return {"pass": False, "failures": failures + ["nonfinite_array"], "metrics": metrics}
    difference = float(np.max(np.abs(stored_costs - expected_costs)))
    metrics["max_cost_weight_difference"] = difference
    require(difference <= 1e-12, "cost_weight_encoding")
    meta = record["model"]
    require(meta.get("objective") == config["objective"], "objective_identity")
    require(meta.get("scenario") == record["scenario"], "scenario_identity")
    gmax = float(meta["g_max"])
    task = gmax * config["scenarios"][record["scenario"]]["task_fraction"]
    require(np.isfinite(gmax) and gmax > 0, "positive_finite_saved_gmax")
    require(abs(float(meta["task_flux"]) - task) <= tol, "task_fraction")
    lower, upper = geometry.lower.copy(), geometry.upper.copy()
    lower[geometry.task_index] = upper[geometry.task_index] = task
    metrics["max_mass_balance_residual"] = float(np.max(np.abs(geometry.matrix @ flux), initial=0))
    metrics["max_bound_violation"] = float(max(0, np.max(lower-flux), np.max(flux-upper)))
    metrics["fixed_task_violation"] = float(abs(flux[geometry.task_index]-task))
    primary_cost = float(np.dot(expected_costs, np.abs(flux)))
    secondary_cost = float(np.abs(flux).sum())
    for name, value in (("primary", primary_cost), ("secondary", secondary_cost)):
        opt = float(record[f"{name}_optimum"])
        cap = float(record[f"{name}_cap"])
        prefix = "cost" if name == "primary" else "secondary"
        formula = opt + config[f"{prefix}_absolute_allowance"] + abs(opt)*config[f"{prefix}_relative_allowance"]
        require(np.isfinite(opt) and np.isfinite(cap), f"{name}_finite_objective")
        require(abs(cap-formula) <= 1e-10*max(1, abs(formula)), f"{name}_cap_formula")
        metrics[f"net_{name}_cost"] = value
        metrics[f"{name}_cap_violation"] = max(0, value-cap)
        require(value >= opt-tol, f"{name}_value_below_claimed_optimum")
    require(record.get("coordinate_order") == panel, "canonical_coordinate_order")
    require(set(record.get("point", {})) == set(panel), "point_panel")
    require(set(record.get("ranges", {})) == set(panel), "range_panel")
    # Retain recovered failures: only the last attempt of each logical LP must be optimal.
    groups = []
    for entry in record.get("solve_ledger", []):
        key = (entry.get("stage"), entry.get("exchange_id", ""))
        if not groups or groups[-1][0] != key:
            groups.append((key, []))
        groups[-1][1].append(entry)
    require([key for key, _ in groups] == expected_stages(panel), "logical_LP_sequence")
    finals = {}
    recovered = 0
    for key, attempts in groups:
        last = attempts[-1]
        value = last.get("objective")
        successful = last.get("status") == "optimal" and isinstance(value, (int, float)) and np.isfinite(value)
        require(successful, f"final_LP_status:{key}")
        if successful: finals[key] = float(value)
        recovered += int(any(a.get("status") != "optimal" for a in attempts[:-1]))
    metrics["logical_LPs"] = len(groups)
    metrics["saved_solver_attempts"] = len(record.get("solve_ledger", []))
    metrics["recovered_LP_groups"] = recovered
    for name in ("primary", "secondary"):
        require(abs(finals.get((name, ''), float('inf')) - record[f"{name}_optimum"]) <= tol, f"{name}_ledger_objective")
    require(abs(finals.get(('final_feasibility', ''), float('inf'))) <= tol, "final_zero_objective")
    index = {rid: j for j, rid in enumerate(geometry.ids)}
    point_diff, range_violation, fix_violation = 0.0, 0.0, 0.0
    for rid in panel:
        if rid not in record.get("point", {}) or rid not in record.get("ranges", {}):
            continue
        point = float(record["point"][rid]); lo, hi = map(float, record["ranges"][rid])
        if not np.isfinite([point, lo, hi]).all():
            require(False, f"nonfinite_point_or_range:{rid}"); continue
        require(lo-hi <= tol, f"reversed_range:{rid}")
        point_diff = max(point_diff, abs(point-flux[index[rid]]))
        range_violation = max(range_violation, lo-point, point-hi)
        require(abs(finals.get(('range_min', rid), float('inf'))-lo) <= tol, f"range_min_ledger:{rid}")
        require(abs(finals.get(('range_max', rid), float('inf'))-hi) <= tol, f"range_max_ledger:{rid}")
        selected = finals.get(('coordinate_min', rid), float('inf'))
        fix_violation = max(fix_violation, abs(point-selected)-config["coordinate_fix_allowance"])
    metrics["max_saved_point_difference"] = point_diff
    metrics["max_point_range_violation"] = max(0, range_violation)
    metrics["max_coordinate_fix_violation"] = max(0, fix_violation)
    for name, value in metrics.items():
        if name.endswith(("violation", "residual")) or name == "max_saved_point_difference":
            require(value <= tol, name)
    # Preserve JSON validity if a missing ledger value produced infinite diagnostic error.
    metrics = {name: value if np.isfinite(value) else None for name, value in metrics.items()}
    return {"pass": not failures, "failures": failures, "metrics": metrics}


def audit_run(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    config = read_json(root / "configs/study.json")
    summary = read_json(output / "run_summary.json")
    if digest(root / "data/raw/Recon3D.json") != RAW_SOURCES["data/raw/Recon3D.json"][0]:
        raise VerificationError("Model differs from the pinned raw Recon3D source")
    check_hash_mapping(root, summary.get("inputs"), FINGERPRINT_FILES)
    fingerprint = hashlib.sha256(json.dumps(summary["inputs"], sort_keys=True).encode()).hexdigest()
    if fingerprint != summary.get("input_fingerprint"):
        raise VerificationError("Run summary fingerprint does not represent its declared inputs")
    if summary.get("config") != config:
        raise VerificationError("Run summary configuration differs from current configuration")
    context_rows = rows(root / "manifests/contexts.tsv")
    development = sorted(r["context_id"] for r in context_rows if r["partition"] == "development")
    mapping = [r for r in rows(root / "manifests/metabolite_map.tsv") if truthy(r["primary_eligible"])]
    panel = sorted(r["exchange_id"] for r in mapping)
    if len(panel) != len(set(panel)) or panel != summary.get("chemistry_exchange_ids"):
        raise VerificationError("Chemical exchange panel is duplicated or differs from the run summary")
    raw = read_json(root / "data/raw/Recon3D.json")
    reactions = {r["id"]: r for r in raw["reactions"]}
    for row in mapping:
        reaction = reactions[row["exchange_id"]]
        if not reaction["id"].startswith("EX_") or len(reaction["metabolites"]) != 1:
            raise VerificationError("Mapped target is not a single-metabolite extracellular exchange")
        met, coefficient = next(iter(reaction["metabolites"].items()))
        if not met.endswith("_e") or float(row["sign_factor"]) != -float(coefficient):
            raise VerificationError(f"Map exchange orientation differs from stoichiometry: {row['exchange_id']}")
    with np.load(root / "data/scores_conservative.npz", allow_pickle=False) as data:
        supported = data["rxn"].astype(str).tolist()
        contexts = data["contexts"].astype(str).tolist()
        evidence = data["A"].copy()
    development_values = evidence[:, [contexts.index(c) for c in development]]
    k = float(np.median(development_values[development_values > 0]))
    if not np.isclose(k, summary.get("k"), atol=1e-12, rtol=0):
        raise VerificationError("Saved scale differs from independently recomputed development scale")
    recipe = load_recipe(root / "src/medium_v2.py")
    geometry_cache, gmax_by_scenario = {}, {}
    results, expected_rows, seen_arms = [], {}, set()
    for declared in summary.get("arms", []):
        identity = (declared["scenario"], declared["context_id"], declared["arm"])
        if identity in seen_arms:
            raise VerificationError(f"Duplicate summary arm: {identity}")
        seen_arms.add(identity)
        scenario, context, arm = identity
        token = hashlib.sha256(f"{context}|{arm}|{scenario}".encode()).hexdigest()[:16]
        record = read_json(output / "arms" / f"{token}.json")
        if (record.get("scenario"), record.get("context_id"), record.get("arm")) != identity:
            raise VerificationError(f"Arm identity mismatch: {identity}")
        if record.get("input_fingerprint") != fingerprint:
            raise VerificationError(f"Arm fingerprint mismatch: {identity}")
        if record.get("status") != "optimal":
            results.append({"identity": identity, "pass": False, "failures": ["arm_status_not_optimal"]})
            continue
        if scenario not in geometry_cache:
            geometry_cache[scenario] = reconstruct_geometry(raw, config, scenario, recipe)
            gmax_by_scenario[scenario] = record["model"]["g_max"]
        if record["model"]["g_max"] != gmax_by_scenario[scenario]:
            raise VerificationError(f"Scenario g_max differs among saved arms: {scenario}")
        geometry = geometry_cache[scenario]
        ci = contexts.index(development[0] if context == "__reference__" else context)
        expected = reconstruct_costs(geometry, supported, evidence[:, ci], k, arm, config["epsilon"])
        fpath = inside(output, record["flux_file"])
        if digest(fpath) != record.get("flux_sha256") or declared.get("flux_sha256") != record.get("flux_sha256"):
            raise VerificationError(f"Flux archive hash mismatch: {identity}")
        with np.load(fpath, allow_pickle=False) as saved:
            result = audit_arm(record, saved["flux"], saved["reaction_ids"].astype(str).tolist(),
                               saved["reaction_costs"], geometry, expected, config, panel)
        results.append({"identity": identity, **result})
        for rid in panel:
            expected_rows[identity + (rid,)] = (record["run_id"], record["point"][rid], *record["ranges"][rid])
    seen_rows = set()
    table_failures = []
    for row in rows(output / "predictions.tsv"):
        identity = (row["scenario"], row["context_id"], row["arm"], row["exchange_id"])
        if identity in seen_rows or identity not in expected_rows:
            table_failures.append({"identity": identity, "failure": "unexpected_duplicate_or_failed_prediction"})
            continue
        seen_rows.add(identity)
        current = (row["run_id"], float(row["point"]), float(row["range_lo"]), float(row["range_hi"]))
        if current != expected_rows[identity] or row["status"] != "optimal":
            table_failures.append({"identity": identity, "failure": "table_differs_from_saved_arm"})
    if seen_rows != set(expected_rows):
        table_failures.append({"failure": "missing_prediction_rows"})
    complete = (summary.get("failed") == 0 and len(results) == summary.get("jobs") == summary.get("completed") and bool(results))
    passed = complete and not table_failures and all(r["pass"] for r in results)
    return {"schema_version": 1, "status": "passed" if passed else "failed", "pass": passed,
            "run_directory": str(output), "partition": summary.get("partition"),
            "run_summary_sha256": digest(output / "run_summary.json"), "prediction_sha256": digest(output / "predictions.tsv"),
            "source_fingerprint": fingerprint, "auditor_sha256": digest(Path(__file__)),
            "arms_checked": len(results), "complete_job_accounting": complete,
            "scenario_geometries_built": sorted(geometry_cache), "development_scale_k": k,
            "table_failures": table_failures, "arms": results,
            "limits": ["No LP, FVA, parameter fitting or outcome evaluation was rerun.",
                       "g_max and objective extrema are checked for internal consistency, not independently re-optimized.",
                       "Final logical LP statuses are required optimal; visible failed attempts followed by recovery are reported.",
                       "No CORE measurement values are read; within-metabolite predictive validity is evaluated separately."]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--run-directory", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        report = audit_run(args.root, args.run_directory)
    except (ValueError, OSError, KeyError, TypeError, IndexError) as exc:
        report = {"status": "failed", "pass": False, "error": str(exc), "exception": type(exc).__name__}
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
