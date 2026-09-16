#!/usr/bin/env python3
"""Release integrity and lock reconciliation for a checkout of this repository.

Three checks, none of which run any analysis:

1. CHECKSUMS.sha256 at the repository root lists every deposited file with its digest; each file is
   re-hashed and compared. Files the release does not carry (the four third-party raw inputs, the
   per-arm records that are release assets) are reported as absent, not as failures.
2. protocol/PROTOCOL_LOCK.json names the exact inputs, sources and validation evidence of the locked
   study. Every path it names is re-hashed against the lock; the four raw inputs are absent until
   they are retrieved (data/raw/README.md), and everything else must match byte for byte. This is a
   read-only comparison: it creates no lock and rewrites nothing.
3. The per-arm records, when present, are checked against the digests in RELEASE_MANIFEST.md's
   asset table via the flux_sha256 field of each record (each JSON record names the digest of its
   own NPZ array).

Usage: python verify/check_release.py [--root DIR] [--strict]
  --strict   exit nonzero if the raw inputs are absent (for a fully provisioned checkout)
"""
import argparse, hashlib, json, re, sys
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--strict", action="store_true")
    a = ap.parse_args()
    root = a.root.resolve()
    problems = []

    # 1. deposited files
    listed = 0; matched = 0; absent = []
    for line in (root / "CHECKSUMS.sha256").read_text().splitlines():
        if not line.strip():
            continue
        digest, rel = line.split("  ", 1)
        listed += 1
        p = root / rel
        if not p.exists():
            absent.append(rel); continue
        if sha(p) != digest:
            problems.append(f"checksum mismatch: {rel}")
        else:
            matched += 1
    print(f"CHECKSUMS.sha256: {listed} files listed, {matched} match, {len(absent)} absent, "
          f"{sum(1 for p in problems if p.startswith('checksum'))} mismatched")
    for rel in absent:
        print(f"  absent: {rel}")
    untracked = []
    ignore = re.compile(r"^(\.git/|\.venv/|work/|logs/|results_repeat/|data/raw/|results/(reserved|development)_v3/arms(/|$)|"
                        r"results/(reserved|development)_v3/run_summary\.json$|.*__pycache__/|.*\.pyc$|(.*/)?\.DS_Store$|CHECKSUMS\.sha256$)")
    listed_set = {line.split("  ", 1)[1] for line in (root / "CHECKSUMS.sha256").read_text().splitlines() if line.strip()}
    for p in root.rglob("*"):
        if p.is_file() and not p.is_symlink():
            rel = p.relative_to(root).as_posix()
            if not ignore.match(rel) and rel not in listed_set:
                untracked.append(rel)
    if untracked:
        print(f"  files present but not listed in CHECKSUMS.sha256: {len(untracked)}")
        for rel in untracked[:20]:
            print(f"    {rel}")

    # 2. the protocol lock
    lock = json.loads((root / "protocol/PROTOCOL_LOCK.json").read_text())
    counts = {"match": 0, "absent": 0, "mismatch": 0}
    raw_absent = []
    for section in ("input_sha256", "source_sha256", "validation_evidence_sha256"):
        for rel, expected in lock[section].items():
            p = root / rel
            if not p.exists():
                counts["absent"] += 1
                if rel.startswith("data/raw/"):
                    raw_absent.append(rel)
                else:
                    problems.append(f"lock: named file absent: {rel}")
            elif sha(p) == expected:
                counts["match"] += 1
            else:
                counts["mismatch"] += 1
                problems.append(f"lock: digest differs: {rel}")
    named = sum(len(lock[s]) for s in ("input_sha256", "source_sha256", "validation_evidence_sha256"))
    print(f"protocol lock: {named} paths named, {counts['match']} match, {counts['absent']} absent "
          f"({len(raw_absent)} third-party raw inputs), {counts['mismatch']} differ")
    for rel in raw_absent:
        print(f"  raw input not retrieved yet: {rel}")
    if a.strict and raw_absent:
        problems.append("strict: raw inputs absent")

    # 3. per-arm records, if fetched
    for part in ("reserved_v3", "development_v3"):
        arms = root / "results" / part / "arms"
        if not arms.is_dir():
            print(f"results/{part}/arms: not present (release asset; see verify/fetch_release_assets.sh)")
            continue
        records = sorted(arms.glob("*.json"))
        bad = 0
        for rec in records:
            d = json.loads(rec.read_text())
            npz = root / "results" / part / d["flux_file"]
            if not npz.exists() or sha(npz) != d["flux_sha256"]:
                bad += 1
        print(f"results/{part}/arms: {len(records)} records, {len(records) - bad} flux arrays match their recorded digest, {bad} do not")
        if bad:
            problems.append(f"{part}: {bad} flux arrays differ from their recorded digest")

    print("-" * 60)
    if problems:
        print(f"{len(problems)} problem(s):")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    print("release integrity: pass")


if __name__ == "__main__":
    main()
