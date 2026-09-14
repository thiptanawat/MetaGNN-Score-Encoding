# Normalization applied when the working tree was packaged

Nothing in the analysis changed when these files were prepared for release, but two kinds of
working-tree detail were normalized and are recorded here rather than left unstated.

**Absolute paths.** Several provenance records stored the full filesystem path of each hashed
input, including the account names of the machines the analyses ran on. Those prefixes were
removed so that each record names its input by a stable relative path. The SHA-256 values
themselves were not touched, so every recorded hash still verifies against the file it names.

**A working directory name.** One earlier working directory carried a name describing an
internal review stage. References to it were replaced with the neutral name `provenance`. The
file it points to, and its recorded hash, are unchanged.

**One amendment field.** The endpoint reconstruction's amendment record carried a field
describing who asked for the analysis. It was replaced with a description of what the record is.
The plan, the rules and the results in that record are unchanged.

The files affected were `endpoint/README.md`, `endpoint/run_endpoint_audit.py`,
`endpoint/AUDIT_AMENDMENT.json`, `results/coverage/agreement_summary.json`,
`results/endpoint/coverage_summary.json` and `results/reserved_range/run_manifest.json`.
