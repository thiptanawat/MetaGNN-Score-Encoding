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
