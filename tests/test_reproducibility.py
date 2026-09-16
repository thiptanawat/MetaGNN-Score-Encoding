"""Integrity/gate tests with temporary fixtures; never run cohort optimization."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from copy import deepcopy

import pytest
from src import verify_inputs as v


def put(root, relative, text="fixture\n"):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


@pytest.fixture
def study(tmp_path, monkeypatch):
    root = tmp_path / "study with spaces"
    root.mkdir()
    raw = put(root, "data/raw/source.txt", "original raw data\n")
    monkeypatch.setattr(v, "RAW_SOURCES", {"data/raw/source.txt": (v.digest(raw), "https://example.invalid/original")})
    monkeypatch.setattr(v, "REFERENCE_INPUTS", {})
    monkeypatch.setattr(v, "DATA_OUTPUTS", ("data/derived.tsv", "manifests/contexts.tsv", "manifests/metabolite_map.tsv", "data/outcomes_long.tsv", "data/assay_noise_development.tsv"))
    monkeypatch.setattr(v, "DATA_SOURCE_FILES", ("src/prepare_data.py", "environment.lock"))
    monkeypatch.setattr(v, "SOURCE_FILES", ("src/prepare_data.py", "src/analysis.py"))
    monkeypatch.setattr(v, "LOCK_INPUTS", ("configs/study.json", "environment.lock", v.PLAN_PATH, v.DATA_RECEIPT) + v.DATA_OUTPUTS)
    for path in set(v.DATA_OUTPUTS + v.DATA_SOURCE_FILES + v.SOURCE_FILES + ("configs/study.json",)):
        put(root, path)
    put(root, "configs/study.json", json.dumps({"prediction_tie_tolerance": 1e-5,
                                               "interval_resolution_tolerance": 1e-5,
                                               "validation_tolerance": 1e-6}))
    v.seal_data(root)
    required = ("manifests/contexts.tsv", "manifests/metabolite_map.tsv", "data/outcomes_long.tsv", "data/assay_noise_development.tsv")
    put(root, v.PLAN_PATH, json.dumps({"input_sha256": v.hashes(root, required), "primary_targets": [{"metabolite_id": "one"}], "point_tolerance": 1e-5, "interval_tolerance": 1e-5, "numerical_tolerance": 1e-6}))
    put(root, "manifests/analysis_panel.tsv")
    put(root, "manifests/development_observed_pairs.tsv")
    return root


def lock(study):
    record = {"schema_version": 1, "status": "locked", "reserved_evaluation_authorized": True,
              "locked_at": "2026-09-13T00:00:00+00:00", **v.lock_requirements(study)}
    v.write_json(study / "protocol/PROTOCOL_LOCK.json", record)
    return record


def test_raw_missing_and_modified_inputs_fail(study):
    path = study / "data/raw/source.txt"
    path.write_text("wrong")
    with pytest.raises(v.VerificationError, match="mismatch"):
        v.verify_raw(study)
    path.unlink()
    with pytest.raises(v.VerificationError, match="missing"):
        v.verify_raw(study)


@pytest.mark.parametrize("path", ["data/derived.tsv", "src/prepare_data.py"])
def test_prepared_receipt_detects_output_and_builder_drift(study, path):
    v.verify_prepared(study)
    put(study, path, "changed\n")
    with pytest.raises(v.VerificationError, match="mismatch"):
        v.verify_prepared(study)


def test_lock_is_never_implicitly_created(study):
    v.lock_requirements(study)
    assert not (study / "protocol/PROTOCOL_LOCK.json").exists()
    with pytest.raises(v.VerificationError, match="Cannot read JSON"):
        v.verify_protocol_lock(study)


def test_exact_matching_lock_passes_and_prevents_data_reseal(study):
    record = lock(study)
    assert v.verify_protocol_lock(study) == record
    with pytest.raises(v.VerificationError, match="Refusing"):
        v.seal_data(study)


@pytest.mark.parametrize("field,value", [("reserved_evaluation_authorized", False), ("status", "draft"), ("locked_at", ""), ("schema_version", 99)])
def test_unauthorized_or_incomplete_lock_rejected(study, field, value):
    record = lock(study)
    record[field] = value
    v.write_json(study / "protocol/PROTOCOL_LOCK.json", record)
    with pytest.raises(v.VerificationError):
        v.verify_protocol_lock(study)


@pytest.mark.parametrize("path", ["configs/study.json", "src/analysis.py", v.PLAN_PATH])
def test_lock_rejects_configuration_source_or_plan_change(study, path):
    lock(study)
    if path in (v.PLAN_PATH, "configs/study.json"):
        plan = v.read_json(study / path)
        plan["changed"] = True
        v.write_json(study / path, plan)
    else:
        put(study, path, "changed\n")
    with pytest.raises(v.VerificationError, match="mismatch|does not match"):
        v.verify_protocol_lock(study)


def test_lock_cannot_omit_a_required_source(study):
    record = lock(study)
    record["source_sha256"].pop("src/analysis.py")
    v.write_json(study / "protocol/PROTOCOL_LOCK.json", record)
    with pytest.raises(v.VerificationError, match="omits"):
        v.verify_protocol_lock(study)


def test_lock_path_and_manifest_cannot_escape_study(study, tmp_path):
    record = lock(study)
    other = tmp_path / "other.json"
    v.write_json(other, record)
    with pytest.raises(v.VerificationError, match="canonical"):
        v.verify_protocol_lock(study, other)
    outside = put(tmp_path, "outside.txt")
    for relative in ("../outside.txt", str(outside)):
        with pytest.raises(v.VerificationError, match="Unsafe"):
            v.check_hash_mapping(study, {relative: v.digest(outside)})
    (study / "escape.txt").symlink_to(outside)
    with pytest.raises(v.VerificationError, match="escapes"):
        v.check_hash_mapping(study, {"escape.txt": v.digest(outside)})


@pytest.fixture
def shell_study(tmp_path):
    root = tmp_path / "runner with spaces"
    root.mkdir()
    shutil.copyfile(Path(__file__).resolve().parents[1] / "reproduce.sh", root / "reproduce.sh")
    put(root, "src/production_controls.py", "# placeholder for CLI-routing fixture only\n")
    fake = put(root, "fake-python", f"#!{sys.executable}\n" + r'''
import json, os, sys
args=sys.argv[1:]
with open(os.environ['CALLS'], 'a') as f: f.write(json.dumps(args)+'\n')
module=args[1] if len(args)>1 and args[0]=='-m' else ''
if os.environ.get('FAIL_MODULE')==module: sys.exit(9)
if module=='src.verify_inputs' and '--check' in args:
    check=args[args.index('--check')+1]
    if os.environ.get('FAIL_CHECK')==check: sys.exit(8)
''')
    fake.chmod(0o755)
    calls = root / "calls.jsonl"
    environment = {**os.environ, "PYTHON": str(fake), "CALLS": str(calls), "RUN_ID": "fixture"}
    return root, calls, environment


def invoke(shell_study, args, extra=None):
    root, calls, environment = shell_study
    result = subprocess.run(["bash", str(root / "reproduce.sh"), *args], cwd=root.parent,
                            env={**environment, **(extra or {})}, text=True, capture_output=True)
    records = [json.loads(line) for line in calls.read_text().splitlines()] if calls.exists() else []
    return result, records


def test_default_preflight_is_read_only_stage(shell_study):
    result, calls = invoke(shell_study, [])
    assert result.returncode == 0, result.stderr
    assert len(calls) == 1 and calls[0][-1] == "preflight"
    assert not (shell_study[0] / "results").exists()


def test_smoke_stage_order_uses_new_builders_and_never_evaluates(shell_study):
    result, calls = invoke(shell_study, ["--stage", "smoke"])
    assert result.returncode == 0, result.stderr
    modules = [args[1] for args in calls if len(args)>1 and args[0]=='-m' and args[1]!='src.verify_inputs']
    assert modules == ["src.prepare_data", "src.prepare_analysis", "pytest", "src.production_controls", "src.run_study", "src.audit_results"]
    run = next(args for args in calls if args[1] == "src.run_study")
    assert run[run.index("--max-contexts")+1] == "1"
    assert run[run.index("--partition")+1] == "development"
    assert "--smoke" in calls[-1]


def test_full_development_routes_each_scenario_to_both_summaries(shell_study):
    result, calls = invoke(shell_study, ["--stage", "all"])
    assert result.returncode == 0, result.stderr
    evaluation = [args for args in calls if args[1] == "src.evaluation"]
    transform = [args for args in calls if args[1] == "src.transform_sensitivity"]
    for commands in (evaluation, transform):
        assert [args[args.index("--scenario")+1] for args in commands] == ["primary", "half_serum", "lower_task"]
    audit_index = next(i for i,args in enumerate(calls) if args[1] == "src.audit_results")
    evaluation_index = next(i for i,args in enumerate(calls) if args[1] == "src.evaluation")
    assert audit_index < evaluation_index


def test_cached_array_corruption_is_not_reused(study, monkeypatch):
    from src import run_study, cpu_model
    monkeypatch.setattr(run_study, "input_fingerprint", lambda root, config: ("a"*64, {}))
    monkeypatch.setattr(cpu_model, "load_config", lambda root: {})
    output = study / "results/fixture"
    array = put(output, "arms/one.npz", "saved array bytes")
    record = {"input_fingerprint": "a"*64, "status": "optimal", "scenario": "primary", "context_id": "one",
              "arm": "ordinal", "flux_file": "arms/one.npz", "flux_sha256": v.digest(array)}
    v.write_json(output / "arms/one.json", record)
    assert v.verify_cache(study, output)["verified_successful_cached_arms"] == 1
    array.write_text("corrupted array")
    with pytest.raises(v.VerificationError, match="hash mismatch"):
        v.verify_cache(study, output)


def test_cached_source_drift_requires_new_output(study, monkeypatch):
    from src import run_study, cpu_model
    monkeypatch.setattr(run_study, "input_fingerprint", lambda root, config: ("b"*64, {}))
    monkeypatch.setattr(cpu_model, "load_config", lambda root: {})
    output = study / "results/fixture"
    v.write_json(output / "arms/one.json", {"input_fingerprint": "a"*64, "status": "failed"})
    with pytest.raises(v.VerificationError, match="Cached inputs differ"):
        v.verify_cache(study, output)


def test_failed_stage_stops_before_downstream_work(shell_study):
    result, calls = invoke(shell_study, ["--stage", "all"], {"FAIL_MODULE": "src.prepare_data"})
    assert result.returncode == 9
    assert not any(args[1] in ("src.prepare_analysis", "src.run_study", "src.evaluation") for args in calls)


def test_reserved_gate_failure_precedes_any_build_or_run(shell_study):
    result, calls = invoke(shell_study, ["--stage", "all", "--partition", "test"], {"FAIL_CHECK": "lock"})
    assert result.returncode == 8
    assert all(args[1] == "src.verify_inputs" for args in calls)
    assert calls[-1][-1] == "lock"


def test_reserved_subset_is_rejected_before_python(shell_study):
    result, calls = invoke(shell_study, ["--stage", "run", "--partition", "test", "--max-contexts", "1"])
    assert result.returncode == 2 and not calls


def test_existing_output_requires_explicit_resume(shell_study):
    root, _, _ = shell_study
    (root / "results/development/arms").mkdir(parents=True)
    result, calls = invoke(shell_study, ["--stage", "run"])
    assert result.returncode == 2
    assert not any(args[1] == "src.run_study" for args in calls)


def test_resume_verifies_cache_before_run(shell_study):
    root, _, _ = shell_study
    (root / "results/development/arms").mkdir(parents=True)
    result, calls = invoke(shell_study, ["--stage", "run", "--resume"])
    assert result.returncode == 0, result.stderr
    cache = next(i for i,args in enumerate(calls) if "cache" in args)
    run = next(i for i,args in enumerate(calls) if args[1] == "src.run_study")
    assert cache < run


@pytest.fixture
def saved_toy():
    import numpy as np
    from scipy.sparse import csr_matrix
    from src.audit_results import Geometry
    geometry = Geometry(["IN", "R", "EX_t_e"], csr_matrix([[1., -1., 0.], [0., 1., -1.]]),
                        np.array([1., 0., 0.]), np.array([1., 10., 10.]),
                        np.array([True, False, True]), 1)
    config = {"objective": "R", "scenarios": {"primary": {"task_fraction": 1.0}},
              "validation_tolerance": 1e-6, "cost_absolute_allowance": 1e-7,
              "cost_relative_allowance": 1e-9, "secondary_absolute_allowance": 1e-7,
              "secondary_relative_allowance": 1e-9, "coordinate_fix_allowance": 1e-8}
    record = {"status": "optimal", "scenario": "primary", "model": {"objective": "R", "scenario": "primary", "g_max": 1., "task_flux": 1.},
              "primary_optimum": 1., "primary_cap": 1.+1e-7+1e-9,
              "secondary_optimum": 3., "secondary_cap": 3.+1e-7+3e-9,
              "point": {"EX_t_e": 1.}, "ranges": {"EX_t_e": [1., 1.]}, "coordinate_order": ["EX_t_e"],
              "solve_ledger": [{"stage": name, "exchange_id": rid, "status": "optimal", "objective": value}
                               for name,rid,value in [("primary","",1.),("secondary","",3.),("range_min","EX_t_e",1.),
                                                      ("range_max","EX_t_e",1.),("coordinate_min","EX_t_e",1.),("final_feasibility","",0.)]]}
    return record, np.ones(3), geometry.ids, np.array([0.,1.,0.]), geometry, np.array([0.,1.,0.]), config, ["EX_t_e"]


def test_independent_vector_audit_accepts_known_feasible_solution(saved_toy):
    from src.audit_results import audit_arm
    result = audit_arm(*saved_toy)
    assert result["pass"], result
    assert result["metrics"]["max_mass_balance_residual"] == 0


@pytest.mark.parametrize("tamper,expected", [("flux", "max_mass_balance_residual"), ("weights", "cost_weight_encoding"),
    ("cap", "primary_cap_formula"), ("point", "max_saved_point_difference"), ("axis", "reaction_axis"),
    ("range", "reversed_range:EX_t_e"), ("task", "task_fraction")])
def test_independent_vector_audit_catches_tampering(saved_toy, tamper, expected):
    from src.audit_results import audit_arm
    values = list(deepcopy(saved_toy))
    if tamper == "flux": values[1][0] += .1
    elif tamper == "weights": values[3][1] = 2.
    elif tamper == "cap": values[0]["primary_cap"] += 1.
    elif tamper == "point": values[0]["point"]["EX_t_e"] = 2.
    elif tamper == "axis": values[2] = list(reversed(values[2]))
    elif tamper == "range": values[0]["ranges"]["EX_t_e"] = [1., .5]
    elif tamper == "task": values[0]["model"]["task_flux"] = .5
    result = audit_arm(*values)
    assert not result["pass"] and expected in result["failures"], result


def test_independent_audit_distinguishes_recovered_attempt_from_final_failure(saved_toy):
    from src.audit_results import audit_arm
    values = list(deepcopy(saved_toy))
    values[0]["solve_ledger"].insert(0, {"stage": "primary", "exchange_id": "", "status": "infeasible", "objective": None})
    result = audit_arm(*values)
    assert result["pass"] and result["metrics"]["recovered_LP_groups"] == 1
    values[0]["solve_ledger"][-1]["status"] = "infeasible"
    values[0]["solve_ledger"][-1]["objective"] = None
    result = audit_arm(*values)
    assert not result["pass"] and any(f.startswith("final_LP_status") for f in result["failures"])


def test_independent_ranking_keeps_true_ties():
    from src.audit_results import average_ranks
    assert average_ranks([4., 1., 1., 8.]).tolist() == [3., 1.5, 1.5, 4.]


def test_analysis_rejects_model_scoring_tolerance_mismatch(study):
    cfg=v.read_json(study/'configs/study.json')
    cfg['prediction_tie_tolerance']=1e-4
    v.write_json(study/'configs/study.json',cfg)
    with pytest.raises(v.VerificationError,match='tolerance mismatch'):
        v.verify_analysis(study)
