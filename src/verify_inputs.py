"""Read-only checks and explicit reproduction receipts; never authorizes a test run.

No network access, downloads, model fitting or optimization occurs in this module.
The ``seal-data`` action records a completed data rebuild; it is not a protocol lock.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
RAW_SOURCES = {
    "data/raw/Recon3D.json": ("aba925f17547a42f9fdb4c1f685d89364cbf4979bbe7862e9f793af7169b26d5", "http://bigg.ucsd.edu/static/models/Recon3D.json"),
    "data/raw/cellminer.zip": ("334dca10d425eaf9fa3c7fc386f8d413ef32d8408eb4fc5508ca75e6cf6b734e", "https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__Affy_HG_U133%28A_B%29_GCRMA.zip"),
    "data/raw/core/NIHMS419088-supplement-Database_S1.xls": ("9a089755b35e16fd17d61b50865ca695bbb33ec8bc8616b285cd5e1fdbfa2ce4", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3526189/"),
    "data/raw/core/NIHMS419088-supplement.pdf": ("1bf01aa2715af80924aa3dbc74c00507259f9517877f23e3157344966cc0f187", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3526189/"),
}
REFERENCE_INPUTS = {"protocol/chemical_identity_reference.tsv": "6512193df6635edba02ec46d4ecc758a633b2650048c060aec7dc552f199e752"}
DATA_OUTPUTS = (
    "data/scores_conservative.npz", "data/gene_expression.npz",
    "data/expression_summary.json", "data/preparation_summary.json",
    "data/outcomes_long.tsv", "data/assay_noise_development.tsv",
    "data/source_supplement_layout.txt", "manifests/contexts.tsv",
    "manifests/metabolite_map.tsv", "manifests/source_table_s1.tsv",
    "manifests/probe_gene_map.tsv", "manifests/model_gene_map.tsv",
    "manifests/gpr_coverage.tsv", "manifests/expression_archive_members.tsv",
    "manifests/raw_sources.tsv",
    "manifests/chemical_identity_review.tsv", "manifests/chemical_reference_provenance.json",
)
DATA_SOURCE_FILES = ("src/prepare_data.py", "src/gpr.py", "environment.lock", "protocol/chemical_identity_reference.tsv")
PLAN_PATH = "protocol/analysis_plan.json"
DATA_RECEIPT = "manifests/reproduction_data.json"
SOURCE_FILES = (
    "src/__init__.py", "src/prepare_data.py", "src/gpr.py",
    "src/prepare_analysis.py", "src/evaluation.py", "src/cpu_model.py",
    "src/cpu_readout.py", "src/medium_v2.py", "src/run_study.py",
    "src/production_controls.py", "src/verify_inputs.py", "reproduce.sh",
    "src/audit_results.py",
    "src/transform_sensitivity.py",
)
LOCK_INPUTS = (
    "configs/study.json", "environment.lock", PLAN_PATH, DATA_RECEIPT, "protocol/chemical_identity_reference.tsv",
    "protocol/CPU_DESIGN.md",
    "manifests/analysis_panel.tsv", "manifests/development_observed_pairs.tsv",
) + DATA_OUTPUTS


class VerificationError(ValueError):
    pass


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inside(root, relative):
    """All manifests use relative files within this study, including symlink targets."""
    root = Path(root).resolve()
    rel = Path(relative)
    if rel.is_absolute() or ".." in rel.parts:
        raise VerificationError(f"Unsafe manifest path: {relative}")
    path = (root / rel).resolve()
    if not path.is_relative_to(root):
        raise VerificationError(f"Manifest path escapes study: {relative}")
    if not path.is_file():
        raise VerificationError(f"Required file is missing: {relative}")
    return path


def hashes(root, paths):
    return {str(path): digest(inside(root, path)) for path in sorted(set(paths))}


def check_hash_mapping(root, mapping, required=()):
    if not isinstance(mapping, dict) or not mapping:
        raise VerificationError("Expected a nonempty relative-path SHA256 mapping")
    missing = sorted(set(required) - set(mapping))
    if missing:
        raise VerificationError(f"Hash mapping omits required files: {', '.join(missing)}")
    for path, expected in mapping.items():
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise VerificationError(f"Invalid SHA256 for {path}")
        if digest(inside(root, path)) != expected:
            raise VerificationError(f"SHA256 mismatch: {path}")


def read_json(path):
    try:
        obj = json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise VerificationError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(obj, dict):
        raise VerificationError(f"Expected JSON object: {path}")
    return obj


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def verify_raw(root):
    errors, found = [], {}
    for relative, (expected, url) in RAW_SOURCES.items():
        try:
            path = inside(root, relative)
            actual = digest(path)
            if actual != expected:
                raise VerificationError(f"SHA256 mismatch: {relative}")
            found[relative] = {"sha256": actual, "bytes": path.stat().st_size}
        except VerificationError as exc:
            errors.append(f"{exc}. Obtain the exact original from {url}")
    if errors:
        raise VerificationError("\n".join(errors))
    return found


def verify_references(root):
    if REFERENCE_INPUTS:
        check_hash_mapping(root, REFERENCE_INPUTS)
    return REFERENCE_INPUTS


def verify_environment(root):
    problems = []
    for line in inside(root, "environment.lock").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        if "==" not in line:
            raise VerificationError(f"Unsupported environment requirement: {line}")
        name, version = line.split("==", 1)
        try:
            actual = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            actual = "missing"
        if actual != version:
            problems.append(f"{name}: expected {version}, found {actual}")
    if shutil.which("pdftotext") is None:
        problems.append("pdftotext is missing; install Poppler and add its bin directory to PATH")
    if problems:
        raise VerificationError("Environment does not match environment.lock:\n" + "\n".join(problems))
    return {"python": sys.executable, "python_version": sys.version.split()[0],
            "pdftotext": shutil.which("pdftotext"), "environment_sha256": digest(root / "environment.lock")}


def seal_data(root):
    for name in ("PROTOCOL_LOCK.json", "PROTOCOL_LOCK.md"):
        if (root / "protocol" / name).exists():
            raise VerificationError("Refusing to reseal data after a protocol lock exists")
    verify_raw(root)
    verify_references(root)
    receipt = {"schema_version": 1,
               "inputs": hashes(root, tuple(RAW_SOURCES) + DATA_SOURCE_FILES),
               "outputs": hashes(root, DATA_OUTPUTS),
               "scope": "Hash receipt of a completed prepare_data build; not biological validation or evaluation authorization"}
    write_json(root / DATA_RECEIPT, receipt)
    return receipt


def verify_prepared(root):
    verify_raw(root)
    verify_references(root)
    receipt = read_json(inside(root, DATA_RECEIPT))
    check_hash_mapping(root, receipt.get("inputs"), tuple(RAW_SOURCES) + DATA_SOURCE_FILES)
    check_hash_mapping(root, receipt.get("outputs"), DATA_OUTPUTS)
    return {"data_receipt_sha256": digest(root / DATA_RECEIPT), "outputs": len(DATA_OUTPUTS)}


def verify_analysis(root):
    verify_prepared(root)
    plan = read_json(inside(root, PLAN_PATH))
    required = ("manifests/contexts.tsv", "manifests/metabolite_map.tsv",
                "data/outcomes_long.tsv", "data/assay_noise_development.tsv")
    check_hash_mapping(root, plan.get("input_sha256"), required)
    for name in ("manifests/analysis_panel.tsv", "manifests/development_observed_pairs.tsv"):
        inside(root, name)
    config = read_json(inside(root, "configs/study.json"))
    for plan_key, config_key in (("point_tolerance", "prediction_tie_tolerance"),
                                  ("interval_tolerance", "interval_resolution_tolerance"),
                                  ("numerical_tolerance", "validation_tolerance")):
        if config_key in config and plan.get(plan_key) != config[config_key]:
            raise VerificationError(f"Analysis/configuration tolerance mismatch: {plan_key}")
    if not plan.get("primary_targets"):
        raise VerificationError("Selected primary panel is empty")
    return plan


def lock_requirements(root):
    """Read-only hashes for the investigator's explicit lock decision; creates no lock."""
    root = Path(root).resolve()
    verify_analysis(root)
    return {"analysis_plan_sha256": digest(root / PLAN_PATH),
            "input_sha256": hashes(root, tuple(RAW_SOURCES) + LOCK_INPUTS),
            "source_sha256": hashes(root, SOURCE_FILES)}


def verify_protocol_lock(root, lock_path=None):
    root = Path(root).resolve()
    canonical = root / "protocol/PROTOCOL_LOCK.json"
    if lock_path is not None:
        candidate = Path(lock_path)
        if not candidate.is_absolute():
            candidate = root / candidate
        if candidate.resolve() != canonical.resolve():
            raise VerificationError("Reserved evaluation requires this study's canonical protocol/PROTOCOL_LOCK.json")
    lock = read_json(canonical)
    if lock.get("schema_version") != 1 or lock.get("status") != "locked":
        raise VerificationError("Protocol lock must have schema_version=1 and status='locked'")
    if lock.get("reserved_evaluation_authorized") is not True:
        raise VerificationError("Protocol lock does not explicitly authorize reserved evaluation")
    if not isinstance(lock.get("locked_at"), str) or not lock["locked_at"].strip():
        raise VerificationError("Protocol lock requires its recorded decision time")
    verify_analysis(root)
    if lock.get("analysis_plan_sha256") != digest(root / PLAN_PATH):
        raise VerificationError("Protocol lock does not match the current analysis plan")
    check_hash_mapping(root, lock.get("input_sha256"), tuple(RAW_SOURCES) + LOCK_INPUTS)
    check_hash_mapping(root, lock.get("source_sha256"), SOURCE_FILES)
    return lock


def verify_controls(root):
    verify_prepared(root)
    result = read_json(inside(root, "results/production_controls.json"))
    if result.get("all_pass") is not True or not result.get("checks"):
        raise VerificationError("Production controls have not all passed")
    if any(check.get("pass") is not True for check in result["checks"].values()):
        raise VerificationError("A production control has failed")
    check_hash_mapping(root, result.get("inputs"),
                       ("configs/study.json", "data/scores_conservative.npz",
                        "manifests/metabolite_map.tsv", "src/cpu_model.py",
                        "src/cpu_readout.py", "src/production_controls.py",
                        "src/medium_v2.py", "manifests/contexts.tsv",
                        "data/raw/Recon3D.json", "environment.lock"))
    return result


def verify_cache(root, output, require_complete=False):
    """Check hashes before native per-arm resume; corrupt arrays are not skipped."""
    from src.run_study import input_fingerprint
    from src.cpu_model import load_config
    output = Path(output).resolve()
    expected, _ = input_fingerprint(root, load_config(root))
    checked, records = 0, {}
    for path in sorted((output / "arms").glob("*.json")):
        record = read_json(path)
        if record.get("input_fingerprint") != expected:
            raise VerificationError(f"Cached inputs differ: {path}. Use a new output directory.")
        if record.get("status") == "optimal":
            relative = record.get("flux_file")
            if not isinstance(relative, str):
                raise VerificationError(f"Successful cached arm has no flux file: {path}")
            if digest(inside(output, relative)) != record.get("flux_sha256"):
                raise VerificationError(f"Cached flux array hash mismatch: {path}")
            checked += 1
            key = (record.get("scenario"), record.get("context_id"), record.get("arm"))
            if key in records:
                raise VerificationError(f"Duplicate cached arm identity: {key}")
            records[key] = record
    if require_complete:
        summary = read_json(output / "run_summary.json")
        if summary.get("input_fingerprint") != expected:
            raise VerificationError("Run summary inputs differ from the current source/data")
        if summary.get("failed") != 0 or summary.get("completed") != summary.get("jobs"):
            raise VerificationError("Optimization run is incomplete or failed")
        if not summary.get("jobs") or checked < summary["jobs"]:
            raise VerificationError("Run summary is missing successful arm artifacts")
        selected = {}
        for arm in summary.get("arms", []):
            key = (arm.get("scenario"), arm.get("context_id"), arm.get("arm"))
            if key in selected or key not in records:
                raise VerificationError(f"Run summary arm is duplicated or missing: {key}")
            record = records[key]
            if arm.get("flux_sha256") != record.get("flux_sha256"):
                raise VerificationError(f"Run summary flux hash differs: {key}")
            selected[key] = record
        if len(selected) != summary["jobs"]:
            raise VerificationError("Run summary arm count differs from completed jobs")
        expected_rows = {(key + (rid,)): (record, rid) for key, record in selected.items()
                         for rid in record.get("point", {})}
        seen = set()
        with inside(output, "predictions.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["scenario"], row["context_id"], row["arm"], row["exchange_id"])
                if key in seen or key not in expected_rows:
                    raise VerificationError(f"Unexpected or duplicate prediction row: {key}")
                seen.add(key)
                record, rid = expected_rows[key]
                numbers = (float(row["point"]), float(row["range_lo"]), float(row["range_hi"]))
                if numbers != (record["point"][rid], *record["ranges"][rid]) or row["status"] != "optimal":
                    raise VerificationError(f"Prediction row differs from its arm artifact: {key}")
        if seen != set(expected_rows):
            raise VerificationError("Prediction table omits completed arm/exchange rows")
        return summary
    return {"verified_successful_cached_arms": checked}


def make_summary(root, output, partition, scenarios, smoke=False):
    verify_analysis(root)
    if partition == "test":
        verify_protocol_lock(root)
    run = verify_cache(root, output, require_complete=True)
    verify_results_audit(root, output)
    if run.get("partition") != partition:
        raise VerificationError("Run partition does not match the requested summary")
    summary = {"created_at": datetime.now(timezone.utc).isoformat(), "partition": partition,
               "smoke_only": bool(smoke), "contexts": run["contexts"], "jobs": run["jobs"],
               "input_fingerprint": run["input_fingerprint"], "evaluations": {}}
    lines = ["# Reproduction execution summary", "",
             "One-context development smoke only; no cohort evaluation or scientific contrast was computed." if smoke else
             ("Development feasibility analysis; this is not a reserved evaluation." if partition == "development" else "Evaluation against the verified, explicitly authorized protocol lock."), ""]
    if not smoke:
        for scenario in scenarios:
            path = output / f"evaluation_{scenario}.json"
            result = read_json(path)
            if result.get("scenario") != scenario or result.get("partition") != partition:
                raise VerificationError(f"Evaluation identity mismatch: {path}")
            if result.get("analysis_plan_sha256") != digest(root / PLAN_PATH):
                raise VerificationError(f"Evaluation plan hash mismatch: {path}")
            if result.get("prediction_sha256") != digest(output / "predictions.tsv"):
                raise VerificationError(f"Evaluation prediction hash mismatch: {path}")
            primary = result.get("primary_fixed_panel", {})
            if primary.get("status") != "complete_predictions":
                raise VerificationError(f"Incomplete scientific evaluation: {path}")
            summary["evaluations"][scenario] = {"file": path.name, "sha256": digest(path),
                "status": primary["status"], "primary": primary.get("primary"),
                "bootstrap": primary.get("bootstrap")}
            descriptive_path = output / f"transform_{scenario}.json"
            descriptive = read_json(descriptive_path)
            if descriptive.get("scenario") != scenario or descriptive.get("partition") != partition:
                raise VerificationError(f"Transform summary identity mismatch: {descriptive_path}")
            if descriptive.get("analysis_plan_sha256") != digest(root / PLAN_PATH) or descriptive.get("predictions_sha256") != digest(output / "predictions.tsv"):
                raise VerificationError(f"Transform summary inputs are stale: {descriptive_path}")
            summary["evaluations"][scenario]["transform_summary"] = {"file": descriptive_path.name, "sha256": digest(descriptive_path)}
            lines += [f"## {scenario}", "", f"Evaluation: [{path.name}]({path.name}).", "",
                      "```json", json.dumps(primary.get("primary"), indent=2, allow_nan=False), "```", ""]
            lines += [f"Descriptive transformation sensitivity: [{descriptive_path.name}]({descriptive_path.name}); this is not accuracy.", ""]
    lines += ["These records copy verified analysis outputs. They do not establish manuscript readiness, equivalence, biological causation or journal acceptance.", ""]
    write_json(output / "reproduction_summary.json", summary)
    (output / "reproduction_summary.md").write_text("\n".join(lines))
    return summary


def verify_results_audit(root, output):
    output = Path(output)
    run = verify_cache(root, output, require_complete=True)
    report = read_json(output / "independent_audit.json")
    if report.get("pass") is not True or report.get("status") != "passed":
        raise VerificationError("Independent saved-vector audit did not pass")
    expected = {"run_summary_sha256": digest(output / "run_summary.json"),
                "prediction_sha256": digest(output / "predictions.tsv"),
                "source_fingerprint": run["input_fingerprint"],
                "auditor_sha256": digest(root / "src/audit_results.py")}
    for field, value in expected.items():
        if report.get(field) != value:
            raise VerificationError(f"Independent audit is stale: {field}")
    return {"audit_sha256": digest(output / "independent_audit.json"), "status": "passed"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--check", required=True, choices=["raw", "environment", "preflight", "seal-data", "prepared", "analysis", "lock", "lock-requirements", "controls", "cache", "complete", "audit", "summary"])
    parser.add_argument("--lock-file", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--partition", choices=["development", "test"], default="development")
    parser.add_argument("--scenarios", nargs="+", default=["primary", "half_serum", "lower_task"])
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        if args.check == "raw": result = verify_raw(root)
        elif args.check == "environment": result = verify_environment(root)
        elif args.check == "preflight": result = {"raw": verify_raw(root), "references": verify_references(root), "environment": verify_environment(root)}
        elif args.check == "seal-data": result = seal_data(root)
        elif args.check == "prepared": result = verify_prepared(root)
        elif args.check == "analysis": result = verify_analysis(root)
        elif args.check == "lock": result = verify_protocol_lock(root, args.lock_file)
        elif args.check == "lock-requirements": result = lock_requirements(root)
        elif args.check == "controls": result = verify_controls(root)
        else:
            if args.output is None:
                raise VerificationError(f"--output is required for {args.check}")
            output = args.output.resolve()
            if args.check in ("cache", "complete"):
                result = verify_cache(root, output, require_complete=args.check == "complete")
            elif args.check == "audit":
                result = verify_results_audit(root, output)
            else:
                result = make_summary(root, output, args.partition, args.scenarios, args.smoke)
        print(json.dumps({"check": args.check, "status": "verified", "result": result}, indent=2, allow_nan=False))
        return 0
    except (VerificationError, OSError, ValueError) as exc:
        print(f"Verification failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
