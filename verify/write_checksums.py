#!/usr/bin/env python3
"""Write CHECKSUMS.sha256 at the repository root: one line per deposited file, `<sha256>  <path>`.

Excluded, because they are not part of the deposited tree: the git metadata, virtual environments,
working directories (work/, logs/, results_repeat/), the third-party raw inputs under data/raw/
(except its README), the per-arm records and run summaries that are release assets, byte-compiled
files and the checksum file itself. verify/check_release.py reads the file this script writes.
"""
import hashlib, re, sys
from pathlib import Path

IGNORE = re.compile(r"^(\.git/|\.venv/|work/|logs/|results_repeat/|data/raw/(?!README\.md)|"
                    r"results/(reserved|development)_v3/arms(/|$)|results/(reserved|development)_v3/run_summary\.json$|"
                    r".*__pycache__/|.*\.pyc$|(.*/)?\.DS_Store$|CHECKSUMS\.sha256$)")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    lines = []
    for p in sorted(root.rglob("*")):
        if p.is_symlink() or not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if IGNORE.match(rel):
            continue
        lines.append(f"{sha(p)}  {rel}")
    (root / "CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")
    print(f"{len(lines)} files listed in CHECKSUMS.sha256")


if __name__ == "__main__":
    main()
