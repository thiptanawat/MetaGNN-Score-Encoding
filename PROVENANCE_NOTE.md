# Normalization applied when the working tree was packaged

Nothing in the analysis changed when these files were prepared for release, but some working-tree
detail was normalized and is recorded here rather than left unstated. `RELEASE_MANIFEST.md` lists
the digest of every normalized file before and after.

**Absolute paths.** Several provenance records stored the full filesystem path of each hashed
input, including the account names of the machines the analyses ran on. Those prefixes were
removed so that each record names its input by a repository-relative path. The SHA-256 values
themselves were not touched, so every recorded hash still verifies against the file it names.
Three records named in the protocol lock (`reports/analysis_reproduction_clean_v3.json`,
`reports/fresh_arm_reproduction_v3/comparison.json`, `results/development_v3/independent_audit.json`)
also contain such paths; they are deposited unchanged, because the lock names their digests.

**A working directory name.** One earlier working directory carried a name describing an internal
review stage. References to it were replaced with the neutral name `provenance`. The file it
points to, and its recorded hash, are unchanged.

**One amendment field.** The endpoint reconstruction's amendment record carried a field describing
who asked for the analysis. It was replaced with a description of what the record is. The plan,
the rules and the results in that record are unchanged, and `endpoint/AUDIT_AMENDMENT.sha256`
names the deposited copy (the original digest is in `RELEASE_MANIFEST.md`).

The files affected in v1.0.0 were `endpoint/README.md`, `endpoint/run_endpoint_audit.py`,
`endpoint/AUDIT_AMENDMENT.json`, `results/coverage/agreement_summary.json`,
`results/endpoint/coverage_summary.json` and `results/mechanism/reserved_range/run_manifest.json`
(under their v1.0.0 paths).

**Version 1.1.0.** The two scenario-ladder run manifests under `results/closure/ladders/` recorded
the absolute path of the Python interpreter that ran them; the prefix was removed so that the field
names the interpreter by a relative path. Nothing else in those records changed.

**Version 1.2.0.** The repository was rearranged into the study layout that the protocol lock and
every script expect (`RELEASE_MANIFEST.md`, sections 1 and 3). Seven further provenance records
were normalized in the same way as above (section 5 of the manifest), and the path constants of
four scripts that still named workstation folders were changed to name the repository (section 4).
No result file, protocol record, locked source or locked input was modified.

**Version 1.2.1.** No provenance record changed. The release replaces the figure script for
figures 4 to 6 and two rendered figures (`RELEASE_MANIFEST.md`, section 8).

**Version 1.3.0.** The release adds the second reconstruction (`second/`, `results/second/`;
`RELEASE_MANIFEST.md`, section 9). Every result record under `results/second/` was scanned for
absolute workstation paths, account names and network addresses. Two run manifests of the range
stages carried absolute paths (`results/second/mechanism/reserved_range/run_manifest.json` and
`results/second/mechanism/bottleneck_pilot/run_manifest.json`); the prefix of the second root was
removed from the names of their inputs and the prefix of the interpreter from its path, no digest
value inside them was changed, and their digests before and after are in `RELEASE_MANIFEST.md`,
section 9. Every other result record is deposited byte for byte. No provenance record of an earlier
release changed. Three files of the identifier audit (`second/audit/`) named the working directory of
the workstation, and the test file of the second reconstruction named the location of the script it
tests. They were changed as follows; no number, class or digest value inside any of them was
touched, and the regenerated report text is the deposited text:

| File | Before | After | Change |
|---|---|---|---|
| `second/audit/README.md` | `4f6d03b32d910b439b980758f4b6292c1d13321ff617fccf5a9e4513d3d62d4e` | `f710e599b5c2f3e9a1cc2e78ea9d15cf74bc227393c5c540762f376277dab8c1` | the Method paragraph names the input files by repository path and by role instead of by workstation path |
| `second/audit/scripts/audit_lib.py` | `17eeca3e7afa7b00b5dcf18bf57de07a02bc985b55643ffe39a85654ca7045fb` | `cf71dc3c6dda6bfb3107ee6ba0e456c43af4fd0b1a7854935b649110b7642c44` | paths resolve from the repository root; the Recon3D JSON of the main analysis, the converted Human-GEM model and the Human-GEM `model/` directory are named by `MAIN_MODEL_JSON`, `SECOND_ROOT_JSON` and `HUMANGEM_MODEL_DIR` |
| `second/audit/scripts/make_report.py` | `6016166df530c7130ce6f0418b41a7e96f266110b55056d83c222241979b34c2` | `4f54ddb6eb14e5c6dc45d506c02012f68ef4e4106db45e5c1d05f1e9745d3c87` | the same Method sentence as in the report |
| `second/tests/test_second_summarize.py` | `103e645d51dc5ac223da658e55f631a468028226ab3164e48d0d46aacc60ea26` | `981c32979fe974f4a370bbcdb3cb3ff340362feaf162f9396555f3122a589d19` | the test file moved from `tests/` to `second/tests/`, and the path to `summarize.py` follows |

The plan (`second/SECOND_RECONSTRUCTION_PLAN_2026-11-14.md`), the four Python scripts and
`run_ranges.sh` of the second analysis, the four build records and the six result tables of the
identifier audit are deposited byte for byte, and the copy of the plan under `results/second/protocol/` has the same
digest as the plan.
