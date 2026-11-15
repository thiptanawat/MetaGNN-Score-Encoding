# Release manifest

Release tag: `v1.3.1`
Archived record: see `SOURCES.md` (the identifier of this release is minted when the release is archived, so the copy of that file inside the release names the identifiers of the releases before it; the repository's `main` branch records the new identifier in the commit that follows the tag).

Release v1.3.1 is release v1.3.0 with one recomputation added (section 10): the reported
quantities of the second reconstruction without the four profiles in which a second solver
completed an arm. It adds the plan of that recomputation, its script, its output and eight tests,
and edits six existing files. Section 10 lists them. No locked file, lock, result or output of the
main analysis or of the second reconstruction changed, and the three release assets are those of
v1.3.0. Release v1.3.0 was release v1.2.1 with the second reconstruction added (section 9).

Release v1.2.1 was release v1.2.0 with the figure script for figures 4 to 6 and two rendered
figures replaced (section 8).

This file states what the release carries, how it relates to the locked study directory the
analyses ran in, what regenerates from what, and what was changed when the working tree was
packaged. It is a packaging record; `protocol/PROTOCOL_LOCK.json` is the scientific record and
was not modified.

## 1. The repository is the study root

The locked study ran in a directory laid out as `configs/`, `data/`, `manifests/`, `protocol/`,
`src/`, `tests/`, `reports/`, `results/`, driven by `reproduce.sh`. Release v1.1.0 deposited the
sources but rearranged the derived data and results, so the driver's commands and the lock's paths
did not resolve from a checkout. Release v1.2.0 restores the study layout: every path named in
`protocol/PROTOCOL_LOCK.json` resolves in the checkout to a file with the recorded digest, except
the four third-party raw inputs, which are not redistributed and are documented in
`data/raw/README.md`.

Lock reconciliation of this release (`python verify/check_release.py`):

| Section of the lock | Paths named | Match byte for byte | Absent |
|---|---:|---:|---|
| `input_sha256` | 30 | 26 | 4, all under `data/raw/` (third-party inputs) |
| `source_sha256` | 14 | 14 | none; `reproduce.sh` is the locked driver, digest `5e7817b817024aabf2185be97b843cd4d6ed391d3d58d84ae6ddb336c84b3652` |
| `validation_evidence_sha256` | 7 | 7 | none |

The lock was not rewritten, re-signed or supplemented. No new lock was created for the post-lock
analyses; those analyses record their own input digests in their run manifests.

## 2. What the release carries, what regenerates, what is verified

| Item | In the release | Regenerates from | Verified by |
|---|---|---|---|
| Locked sources (`src/`, `reproduce.sh`), configuration, protocol, manifests, environment lock | yes, byte for byte | not applicable | lock reconciliation (`verify/check_release.py`) |
| Derived inputs under `data/` (7 files, 3.7 MB) | yes, byte for byte | `reproduce.sh --stage data` from the four raw inputs; `manifests/reproduction_data.json` is the build receipt | lock reconciliation; `python -m src.verify_inputs --root . --check prepared` once the raw inputs are present |
| Completed run summaries, predictions, evaluations, transform sensitivities, independent audits, concordance reconstructions (`results/development_v3/`, `results/reserved_v3/`) | yes | `reproduce.sh --stage all --partition development` and `--partition test` into new directories | `CHECKSUMS.sha256`; three files are lock evidence (see section 1) |
| Per-arm records of the two runs (`results/*/arms/`, 1,698 and 402 files) and the two `run_summary.json` files | release assets (table below), not in the tree | the same two commands | each JSON record names the digest of its NPZ array; `verify/check_release.py` checks every pair when the assets are unpacked; `reports/fresh_arm_reproduction_v3/comparison.json` records a fresh arm that matched its saved record with zero difference in every flux |
| Validation evidence named in the lock (`reports/`, `results/production_controls.json`) | yes, byte for byte | `reproduce.sh --stage controls`, `--stage tests`; the fresh-arm and clean-directory reproductions were one-off records | lock reconciliation |
| Post-lock analyses: code and outputs (`coverage/`, `mechanism/`, `closure/`, `copeland/`, `endpoint/`, `figures/` and their `results/` directories) | yes | `reproduce_extensions.sh <stage>` into `results_repeat/<timestamp>/`; needs the release assets for the mechanism and closure stages, the packages in `environment_extensions.txt` for three stages, and `copeland/fetch_inputs.sh` for the independent evaluation | `CHECKSUMS.sha256`; the run manifests record the digests of their inputs; `verify/verify_manuscript_numbers.py` traces every reported number to these files |
| Second reconstruction: plan, code, driver, build records, identifier audit, tests (`second/`) | yes, byte for byte for the plan, the four Python scripts, `run_ranges.sh` and the build records | `second/reproduce_second.sh` stage by stage from a Human-GEM 2.1.1 checkout, which is not redistributed; its `build` stage prints `identical` or `DIFFERS` for each of the four build records | `CHECKSUMS.sha256`; the plan's SHA-256 `0de8cbf4a737829b4c005b33122568d1abb342166162917f5d5481d4b7cd9a81` is written in the lock of the second root and checked by the driver; `python -m pytest -q second/tests` (5 tests; the locked study's 74 stay under `tests/`) |
| Second reconstruction: outputs (`results/second/`) | yes (the lock and its plan copy, predictions, evaluations, the ledgers, summaries and run manifests of the range stages and of their partial development pilot, the cap sweep and the solver check); the per-arm records and the two `run_summary.json` files are a release asset, not in the tree | the same driver, stages `development` to `evaluate` | `CHECKSUMS.sha256`; `results/second/protocol/PROTOCOL_LOCK.json` names the plan digest; `verify/verify_manuscript_numbers.py` traces every reported number of the second reconstruction to these files; `verify/check_release.py` checks every record of the asset against its flux array when the asset is unpacked |
| Third-party raw inputs (Recon3D, CellMiner archive, the Jain et al. workbook and supplementary PDF) | no | retrieved once as `data/raw/README.md` describes | `python -m src.verify_inputs --root . --check raw` |
| Independent-evaluation inputs (two public git repositories) | no | `copeland/fetch_inputs.sh` clones them and refuses any commit other than the pinned ones | digests of the two `.rda` files checked by that script |
| Human-GEM 2.1.1 (the second reconstruction; public repository, CC BY 4.0) | no | `second/build_root.py` rebuilds the second root from the pinned commit; `second/reproduce_second.sh build` refuses any other commit and checks the five file digests of `second/build_report.json` | the same five digests, listed in `SOURCES.md` |

### Release assets

| Asset | Bytes | SHA-256 | Contents |
|---|---:|---|---|
| `reserved_v3_arms.tar.gz` | 108149951 | `083524997b82bb2e459c2ca078a2c72000566f9ad960b108f11fce23ca91a210` | `arms/` (849 JSON records and 849 NPZ flux arrays of the reserved run: 47 profiles, 6 encodings and 3 scenarios, plus the 3 context-independent uniform-cost records) and `run_summary.json`; unpack into `results/reserved_v3/` |
| `development_v3_arms.tar.gz` | 25567045 | `94709ee80f65a81fed390845c7626622a080ac4c1d3a3634e2ac31dfbc17f351` | `arms/` (201 JSON records and 201 NPZ arrays of the development run) and `run_summary.json`; unpack into `results/development_v3/` |
| `second_reconstruction_arms.tar.gz` | 28732431 | `de9dc74000a43b274afb91694efd1935490daafc3c32890e8c4da694f7bac401` | The per-arm records of the second reconstruction: `results/second/reserved_v3/arms/` (283 JSON records and 283 NPZ flux arrays of the reserved run: 47 profiles under 6 encodings in the primary scenario, plus the uniform-cost record), `results/second/development_v3/arms/` (67 JSON records and 67 NPZ arrays of the development run) and the two `run_summary.json` files; the paths in the archive start at `results/`, so unpack it at the repository root |

`bash verify/fetch_release_assets.sh` downloads the three from the release that matches this
checkout, verifies these digests and unpacks them (`reserved`, `development` or `second` as its
argument fetches one). The first two archives, attached to releases v1.2.1 and v1.3.0, are the
archives of release v1.2.0, unchanged. The third archive is new in v1.3.0. No stage of
`reproduce_extensions.sh` and no check of the reported numbers needs it; when it is unpacked,
`verify/check_release.py` compares each of its flux arrays with the digest in its record.

## 3. Path map from release v1.1.0

Every file of v1.1.0 is present in v1.2.0 with the same content unless listed as replaced.

| v1.1.0 | v1.2.0 |
|---|---|
| `data/inputs/scores_conservative.npz` | `data/scores_conservative.npz` (the lock's path) |
| `data/reserved_v3/*` | `results/reserved_v3/*` |
| `results/evaluation_*.json`, `results/evaluation_*_summary.tsv`, `results/independent_concordance_*.json`, `results/audit_compact.json`, `results/reproduction_summary.*` | `results/reserved_v3/` (the same names) |
| `data/mechanism/reserved_range_stage_ledger.tsv`, `results/reserved_range/*` | `results/mechanism/reserved_range/range_stage_ledger.tsv`, `run_manifest.json`, `summary.json` |
| `data/mechanism/uniform_control_range_stage_ledger.tsv`, `results/uniform_control/summary_uniform.json` | `results/mechanism/uniform_control/range_stage_ledger_uniform.tsv`, `summary_uniform.json` |
| `data/mechanism/cross_cost_*.tsv`, `results/cross_cost.json`, `data/mechanism/range_location_*.tsv`, `results/range_location.json` | `results/mechanism/` (the same names) |
| `data/bottleneck_pilot/development_range_stage_ledger.tsv`, `results/bottleneck_pilot/*` | `results/mechanism/bottleneck_pilot/range_stage_ledger.tsv`, `summary.json`, `readback_check.json` |
| `data/coverage/*`, `results/coverage/*` | `results/coverage/` |
| `data/endpoint/*.tsv`, `results/endpoint/coverage_summary.json` | `results/endpoint/` |
| `data/closure/<run>/predictions.tsv`, `data/closure/ladders/<scenario>/range_stage_ledger.tsv` | `results/closure/<run>/predictions.tsv`, `results/closure/ladders/<scenario>/range_stage_ledger.tsv`, beside their summaries and run manifests |
| `data/copeland/run_predictions.tsv`, `data/copeland/crossed_predictions.tsv` | `results/copeland/run/predictions.tsv`, `results/copeland/crossed/crossed_predictions.tsv` |
| `results/copeland/meantask/` | `results/copeland/run_meantask/` (predictions, run manifest, run summary) and `results/copeland/evaluation_meantask/` (decisions, evaluation) |
| `reproduce.sh` (a hand-written driver whose commands did not match the modules' interfaces) | replaced by the locked driver, byte for byte; the post-lock stages moved to `reproduce_extensions.sh` with the interfaces the scripts have |
| `verify/verify_manuscript_numbers.py` (read the manuscript from a private path and the outputs from workstation folders) | replaced: resolves every artefact inside the repository, takes the manuscript and supplement as optional arguments, and runs values-only without them |
| `closure/solver_check.py` | replaced by the version with the `--panel` selection (chemistry list of 96 exchanges, or the 52-target panel); the first run's outputs stay under `results/closure/solver_check/`, the production-list run is `results/closure/solver_check_96/` |
| `copeland/copeland_rule.py`, `copeland_evaluate.py`, `test_copeland_rule.py` | replaced by the versions that report the planned shared-encoding rule beside the pooled rule (amendment A1 in `copeland/copeland_amendments.json`); the earlier outputs stay under `results/copeland/evaluation/` and `evaluation_meantask/`, the current ones are `evaluation_r5/` and `evaluation_meantask_r5/` |
| `figures/make_figures.py` | replaced by the version that draws figures 4, 5 and 6 from the repository paths |

Added in v1.2.0: `configs/study.json`; the seven derived inputs under `data/`;
`results/production_controls.json`; `results/development_v3/` (evaluations, audits, predictions,
transform sensitivities, reproduction summaries, `run_summary_scalars.json`);
`results/reserved_v3/independent_audit.json`, `interpretation_compact.json` and the three
post-lock range diagnostics; `reports/` (the lock's validation evidence, the hardware record, the
fresh-arm reproduction with its isolated driver and outputs, the validation scripts and their
records); `tests/` (the locked study's seven test modules, 74 tests); the bottleneck pilot's
per-arm records and extrema arrays, `B0.json`, run manifest, execution log and readme; the coverage audit's execution
records; the logs of the mechanism and independent-evaluation runs; `results/closure/solver_check_96/`
and `results/closure/stat_repairs3/` with `closure/stat_repairs3.py`; `results/copeland/prepared/scores_copeland.npz`,
`results/copeland/crossed/run_manifest.json` and the r5 evaluations; `endpoint/Source_Quality_Provenance.json`;
the rendered figures under `results/figures/`; `copeland/fetch_inputs.sh`, `verify/check_release.py`,
`verify/fetch_release_assets.sh`, `verify/write_checksums.py` (the script that writes `CHECKSUMS.sha256`),
`environment_extensions.txt`, `data/raw/README.md` and this file.

## 4. Path constants changed in four scripts

Four scripts carried the paths of the workstation folders the analyses ran in. Their path
constants now name the repository; the computation is unchanged and the outputs they produced are
deposited beside them. Digests before and after (the "before" digest is that of the copy which
produced the deposited outputs; for the first three it is the copy deposited in v1.1.0):

| Script | Before | After | Change |
|---|---|---|---|
| `mechanism/range_location.py` | `10da3883d130b80f616c1edb9e2359b93efda2d977fd1712aab3ed7408381b03` | `e741b836e62d0689d4ea64381f23d1abb9d43116e9d69d634b6342e588ae67b3` | inputs resolved under the repository (`STUDY_ROOT`, `PILOT_LEDGER`, `RESERVED_LEDGER` overrides) |
| `mechanism/cross_cost.py` | `af0fe819b364e7827ad4d6f7bfc50679e652e6ad08836ecbc3bb25ffd8e4a1e8` | `275cac488cdf64b2f807b870644cde1b23bf68c2a871291fe055651caf7c1560` | per-arm records resolved under the repository (`STUDY_ROOT` override) |
| `endpoint/run_endpoint_audit.py` | `352cdce86472d8d68044847c82c414781d98136372e0c9283dfa930f97cc8741` | see `CHECKSUMS.sha256` | study root and provenance file resolved under the repository; outputs written to `results/endpoint/` (`ENDPOINT_OUT` override) |
| `figures/make_figures.py` | `273e56798daa5de54382f0c3fd46afbee4a160cdbbbe95bc8daed73a08d63159` | `0eb6712b02f0c2215b8172e4c3109a3e72a1991a70c3c5d55e969f67c5090cca` | inputs resolved under `results/`, output directory `results/figures/manuscript` (`FIGURE_OUT` override) |
| `mechanism/readback_check.py` | see `CHECKSUMS.sha256` of v1.1.0 | see `CHECKSUMS.sha256` | takes the pilot directory as an optional argument instead of reading its own directory |

`results/mechanism/range_location.json` records its inputs under the keys the script used at the
time (`inputs_sha256`, workstation-relative names); a repeat records repository-relative names with
the same digests.

### Repeats run on the workstation while packaging

With the deposited ledgers and per-arm records as inputs, the patched scripts reproduced the
deposited outputs: the coverage audit's six scenario and arm-set rows agree on every count;
`range_location.json` is identical apart from the names of its inputs; `cross_cost.json` is
identical to 1e-12 in every value apart from the elapsed time, and its three ledgers differ only
in the last digits of floating-point sums (the repeat used the locked NumPy 2.4.4, the original
run a different NumPy build). The endpoint reconstruction flagged the same 119 of 120 cultures.
Figures 4 to 6 rendered from the deposited outputs.

## 5. Records normalized when packaging

As in the earlier releases (`PROVENANCE_NOTE.md`), provenance records that stored absolute
workstation paths were rewritten to repository-relative paths. No digest value inside any record
was changed, and no file named in the protocol lock was touched. The original digests are kept
here so that the archived working copies can be matched to the deposited files.

| File | Original SHA-256 | Deposited SHA-256 |
|---|---|---|
| `results/reserved_v3/independent_audit.json` | `af33baaba4079fbac45ddb3072d180879b6e5f9c6619f4fdf0632d4595618fce` | `c60b1486d7bd908de2d2d859cebfcd4ec4102090a5b80d3689bb8ab00f8363a7` |
| `results/reserved_v3/postlock_range_sensitivity_primary.json` | `0098f8ee2b5ff7dbcb32fd233e2dbcffe428e273eff7f8425548149d61474101` | `389dd65e4f025e5e66bdc1215ec79063f3aec0961be8e635ffcd5db2e2c052cc` |
| `results/reserved_v3/postlock_range_sensitivity_half_serum.json` | `7ec8ddee291c40a8878d2f7c5cc87a72fc977f5f7698bf2ab7320a5fde83253d` | `4b5e415f88d2bad3e0c32c9105ba1d1297766c91bd3b67f41ef444f4bac1cb35` |
| `results/reserved_v3/postlock_range_sensitivity_lower_task.json` | `4810667616bcfe5c49d16c54bfc4c77f3c97106ae713ac047548e40ae10adea9` | `53b1973f279c9bc8d28d263425bebf73a51165c3d182f0765d966437ebdd0f38` |
| `reports/fresh_arm_reproduction_v3/copy_manifest.json` | `0370ffa3dd5ef18da54681f8a3fdad2fc483c4ea7333525714e736e027a412ee` | `6e40ca3d944d7880bb4e922c60352cbe99dd72ba7740fb6f45a0ace5e22e4c79` |
| `results/figures/reserved_v3/figure_manifest.json` | `04661520dffb75f9a668103881095279f3bf0b73ecb2c83469510fe07511637e` | `136c5c04545584530c8d4c5e178bcecbe0226dd0b4353c5813658b5a40cff913` |
| `results/mechanism/bottleneck_pilot/run_manifest.json` | `895b9adb138e90df8f0fde807cef3614f3a3b7fb9e41050bd6ddddd0819614fa` | `06c22771bc5d728831d686a5650b63ace5ce5be6b9b3180f4cd67a7368d48677` |

Three lock-hashed records (`reports/analysis_reproduction_clean_v3.json`,
`reports/fresh_arm_reproduction_v3/comparison.json`, `results/development_v3/independent_audit.json`)
also contain absolute workstation paths, as does `reports/fresh_arm_reproduction_v3/isolated/fresh_arm.json`,
whose digest the comparison record and `artifact_hashes.json` name. They are deposited byte for byte,
because rewriting them would break the records they are named in; the paths they contain identify
the study directory on the workstation and nothing else.

The endpoint reconstruction's amendment record was normalized in v1.0.0 (`PROVENANCE_NOTE.md`), so
`endpoint/AUDIT_AMENDMENT.sha256` names the deposited copy; the original record's digest was
`f7b879fb4051b521a55bfda076b8fac822710f28e0c1fde449f2b5850632c044`.

## 6. Not carried

The two `run_summary.json` files of the study runs (58 MB and 247 MB) travel with the release
assets rather than the tree; the console logs of the study wrapper (332 MB) and the two separate
reserved batch directories that `results/reserved_v3/` was assembled from (`reports/BATCH_ASSEMBLY.json`
in the study directory records the assembly) are not deposited. The narrative review documents
written while the study was assembled are not deposited; the machine records they cite are.

## 7. Clean-checkout verification

The record of the clean-checkout verification is `verify/clean_checkout_log_v1.2.0.md`. It was
performed on release v1.2.0. Release v1.2.1 differs from it only in the four files of section 8,
none of which is an input to an analysis or to the number checker; the checks repeated on a fresh
clone of v1.2.1 are listed in that section. Release v1.3.0 adds files and edits the files listed in
section 9, which also lists the checks repeated on a fresh clone of v1.3.0. Release v1.3.1 adds the
files of section 10 and edits those it lists, with the checks repeated on a fresh clone of v1.3.1.

## 8. Changes in v1.2.1

Release v1.2.1 replaces four files of v1.2.0. No analysis was rerun, and no result table,
protocol record, locked source or locked input was modified.

| File | Digest in v1.2.0 | Digest in v1.2.1 | Change |
|---|---|---|---|
| `figures/make_figures.py` | `0eb6712b02f0c2215b8172e4c3109a3e72a1991a70c3c5d55e969f67c5090cca` | `ea2429f41fc7033f85a8efa97d27bf9e65e805cdb032dfe73cd83ec91fe6a027` | The midrank encoding is labeled "Midrank" in figures 4 and 6 (it was "Ordinal"; its identifier in the code and in the result files stays `ordinal`). The resolved-pair table of figure 4C, the bracket labels and the lactate title of figure 6A, and the note and the inset of figure 6C are placed so that they stay inside the frame and cover no plotted line, and the script now stops if the note or the inset of figure 6C would cover one. The files read and the values drawn are the same |
| `results/figures/manuscript/Figure_4_Constraint_Stage.png` | `91531b4ab33a597dfc3c4da7e709b9fa6de5b7744b6eefbccf0a84ccbc4ec681` | `8bc209380775b6d54f236810c2a6d0dcce8c700d3b305029bdadd37d762de9ae` | rendered with the script above |
| `results/figures/manuscript/Figure_6_Independent_Evaluation.png` | `267030b0ad07f505667368c0da52c17ef11539d27a23ca54ce7d86eec144ba87` | `9a6f0db33dde7104ccbac57e8ad4f32c94b5d40239de6bf00fa0f43ae1361937` | rendered with the script above |
| `results/figures/manuscript/SHA256SUMS.txt` | `bb3965cd310f0b6685f336ca52f295da1a887e32fd8815239f367aa1e75815ba` | `1fc66779a6ef6e30747442efd82b4be82a7f016ef251132761fab45c9a115da8` | names the three figures as the script writes them, so that `sha256sum -c` runs in that directory |

`results/figures/manuscript/Figure_5_Comparator_Formulations.png` is byte-identical to v1.2.0.
The script prints the values it draws; run from v1.2.0 and from v1.2.1 on the same checkout, the
two reports agree on all 260 values. The digest of `figures/make_figures.py` in section 4 is that
of v1.2.0.

Checks repeated on a fresh clone of v1.2.1 before the tag was pushed: `python verify/check_release.py`
(every file listed in `CHECKSUMS.sha256` matches, and the lock reconciliation is that of section
1), `FIGURE_OUT=<new directory> python figures/make_figures.py` (figures 4 to 6 render from the
deposited outputs), and `python verify/verify_manuscript_numbers.py` against the manuscript and
its supplementary article (284 checks, 0 problems).

## 9. Changes in v1.3.0

Release v1.3.0 adds the second reconstruction and edits eleven existing files. No analysis of the
main study was rerun. No locked file, no lock, no result and no output of the main analysis
changed: `protocol/PROTOCOL_LOCK.json` and every file it names are as in v1.2.1, and the lock
reconciliation of section 1 is unchanged. The lock of the second root
(`results/second/protocol/PROTOCOL_LOCK.json`) is a separate record in a separate directory, and
the main lock was not rewritten.

Files added (155):

| Directory | Files | Contents |
|---|---:|---|
| `second/` | 13 | `README.md`, the locked plan `SECOND_RECONSTRUCTION_PLAN_2026-11-14.md` (SHA-256 `0de8cbf4a737829b4c005b33122568d1abb342166162917f5d5481d4b7cd9a81`), `build_root.py`, `run_arms.py`, `in_second_root.py`, `run_ranges.sh`, `summarize.py`, `reproduce_second.sh`, `build_report.json`, `exchange_map.tsv`, `target_map.tsv`, `medium_additions.tsv`, `environment_second.txt` |
| `second/audit/` | 12 | `README.md`, `target_audit.tsv`, `medium_audit.tsv`, `alternatives.tsv`, `medium_alternatives.tsv`, `medium_missing_candidates.tsv`, `yml_crosscheck.txt`, and five scripts under `scripts/` |
| `second/tests/` | 1 | `test_second_summarize.py` (5 tests; kept outside `tests/` so that the locked study's 74 tests stay 74) |
| `results/second/` | 129 | `protocol/` (the lock of the second root and the plan copy it names), `second_summary.json`, `reserved_v3/` (7 files), `development_v3/` (3), `mechanism/` (16: the four files of the interval-position analysis, `reserved_range/` with 3, `uniform_control/` with 4 and `bottleneck_pilot/` with 5), `closure/cap_sweep/` (3 files and 95 job records) and `closure/solver_check_96/` (2) |

Not carried, as for the main runs (section 6): the per-arm records and flux arrays of the second
root and its two `run_summary.json` files, which are the release asset named in section 2; the
per-job records and extrema arrays of the range stages (171 MB for the reserved stage; the ledgers
built from them are deposited); and the logs of the runs. The development pilot of the range
stages on the second root is partial, as the plan states (8 of 22 calculations inside the time
budget fixed in the script); its ledger, summary, run manifest and network-only record are
deposited under `results/second/mechanism/bottleneck_pilot/` because the interval-position analysis
reads the ledger, and its per-job records are not. The plan copy under `results/second/protocol/`
has the SHA-256 of the plan above. The second root's input digests are in its lock, and the
Human-GEM digests are in `second/build_report.json`.

Two run manifests of the second root stored absolute workstation paths and were normalized as in
section 5: the prefix of the second root was removed from the names of the inputs and the prefix
of the interpreter from its path. No digest value inside either record was changed.

| File | Original SHA-256 | Deposited SHA-256 |
|---|---|---|
| `results/second/mechanism/reserved_range/run_manifest.json` | `0df0c6b1b1f2fc066c17c2a77d5694c6539d75ed7f1abf998f93735641d0ccf2` | `6272876193e6d770f975a92d7a4bc58f6d1bf34b347dfaddaa6ce011a1b202df` |
| `results/second/mechanism/bottleneck_pilot/run_manifest.json` | `ce8b1a9dc3509f45401309163c5c001099daca16897d2458b1aaf3a15ba222e1` | `6befeb8805bffe42c7acf5a36b87d679ea8e2bd5380f83a9e2742cce1d426675` |

Files changed (11). The digests before are those of the main branch at commit `cbfa9bb`, that is,
release v1.2.1 followed by the commit that recorded its archived identifier in `CITATION.cff`,
`SOURCES.md` and this file; for the other eight files they are also the digests in release v1.2.1.

| File | Digest before | Digest in v1.3.0 | Change |
|---|---|---|---|
| `README.md` | `a32d889128e7b206d690f73768f5cd873013fb5ba8047299442117a18a7c771a` | `d928d9e0a209e996de07a3c28e9f99663a8d6958152f3162e7ac0e899b2a14d3` | a section on the second reconstruction, one Layout row each for `second/` and `results/second/`, and the driver `second/reproduce_second.sh` under Reproducing; no statement about the main analysis changed |
| `RELEASE_MANIFEST.md` | `639796e7ca0d133389b2d4c5d3645bcc435e091c37ab18ea9db80c40cf02540e` | see `CHECKSUMS.sha256` | this file: tag, opening paragraph, section 2 rows, the third release asset, section 7 and this section |
| `.gitignore` | `ea17fe4a7f6ec2a23be201593767f24423a64dd15d5624314582f33d61c9eabc` | `dacf876115a967aca4e10967d0f929f550806b9b03273b3abe283514b174851b` | four lines: the `arms/` directories and the `run_summary.json` files of `results/second/reserved_v3/` and `results/second/development_v3/`, which the third release asset unpacks |
| `verify/write_checksums.py` | `386bb28e33464ded186a2743be0fa99f0713096c0a7184219de4d4e4e280b8b8` | `a95091b49640de12aa4c8e36be019f422b9d763402c55eda1ca687a8a3e089c9` | the exclusion of the per-arm records and run summaries also covers those under `results/second/` |
| `verify/check_release.py` | `03eb39d92e9b24ab6079fadd4f6060f37afd50be7923adf9d1cab7c930ca905d` | `259315a18bb38180961739c13dc0a7a9129a44f8902daccace7147f57b28cae9` | the same exclusion, and the comparison of each per-arm record with its flux array also runs over the two runs of the second reconstruction when the third asset is unpacked; the checks of the main analysis are unchanged |
| `verify/fetch_release_assets.sh` | `2156caf23db55678e6d54e4ec0701176ca3f85ecb9566017e1cce8c193499fa7` | see `CHECKSUMS.sha256` | fetches, verifies and unpacks the third asset (`second`, and as part of `all`); the first two are handled as before |
| `CITATION.cff` | `ec08e9d7bcc24689cea5264b0500d6f1f6da846dd50097347356c171461ff298` | `4958775ed452b80c5985c2f8793637b8698419964419985411b27b1aff744303` | `version` 1.3.0 and `date-released` 2026-11-15; no new identifier, because the identifier of this release does not exist before it is archived |
| `PROVENANCE_NOTE.md` | `2de0084f376b087415b493ceaa3fa2a2d803717a1e4927719b92054f4b4b18a0` | `b68687ed539cd08df5eb8b776f4bfb3ff3f1b220cb8ad6e7a3e3e1639ecdb53a` | a "Version 1.3.0" paragraph: two run manifests of the second root were normalized (digests above); four added files were edited, with digests before and after |
| `SOURCES.md` | `e44ab0d594c5a51c9006cf6f301ece6886ef0e090f66c7845e49d1257f4560fe` | `7a9b63a43e0bd012a1f84d33cb6efe78cdd06a5e564b9b5c2c649872ea7bda72` | the entry for Human-GEM 2.1.1 (commit, license, five digests, citation, not redistributed) and a sentence on `results/second/`; the archived-record section still names v1.2.1 and earlier |
| `reproduce_extensions.sh` | `4c6e8eddeb0fa3506612b262cbbcf1e242e27cb243423a60611856302ae0a6de` | `155e5581933e0b741626351bd445ee86c93a7d70943b6f341f5f2978d9c26b4e` | the `tests` stage also runs `pytest -q second/tests`; two header lines name `second/reproduce_second.sh`; the usage text prints one more line |
| `verify/verify_manuscript_numbers.py` | `a13bd3279685a597873615cde48b99f7f8c604fd3e11ece066aa4b2738cde6ae` | `0259da88311cda53e2f48256fe950e04b58e69b42f37a20185d62d55e8a862d3` | the 45 checks of Supplementary Methods S22 are added before the report, reading `second/` and `results/second/` and, for the Recon3D values beside them, the cap-sweep files the script already reads; the earlier 284 checks are unchanged |

`CHECKSUMS.sha256` is regenerated and lists the 155 added files.

Checks repeated on a fresh clone of v1.3.0 before the tag was pushed: `python verify/check_release.py`
(every file listed in `CHECKSUMS.sha256` matches, and the lock reconciliation is that of section
1; with the third asset unpacked, the 283 and 67 flux arrays of the second reconstruction match the
digests in their records), `python verify/verify_manuscript_numbers.py` against the manuscript and
its supplementary article (329 checks, 0 problems), `python -m pytest -q second/tests` (5 tests),
`bash second/reproduce_second.sh build` from a Human-GEM checkout at the pinned commit (the four
build records it writes are identical to the deposited copies), `bash second/reproduce_second.sh
evaluate` in the second root that the build wrote, after the lock, the plan copy and the deposited
predictions, ledgers and cap-sweep records were placed in it (the lock accepts the rebuilt inputs,
and the ten files the stage writes are identical to the deposited copies), and the scripts of
`second/audit/scripts/` (the six tables and the report they write are identical to the deposited
copies). The arm runs and the range stages were not repeated for this check.

## 10. Changes in v1.3.1

Release v1.3.1 adds one recomputation to the second reconstruction and edits six existing files. No
analysis was rerun and no linear program was solved. No locked file, no lock, no result and no
output of the main analysis or of the second reconstruction changed, and the three release assets
are those of v1.3.0, with the digests of section 2.

In the reserved run on Human-GEM the second solver completed 4 of the 283 arms (section 9;
`results/second/reserved_v3/second_solver_arms.json`). The recomputation removes the 4 profiles that
own such an arm, under every encoding, and computes the reported quantities again for the remaining
43 profiles from the saved records. Its plan and its script were committed to `main` on 15 November
2026 (commit `18b7168ef8f197877cad0934e3139946e7bdf0b1`), when the script had been run only with
`--check-only`, which reproduces the deposited values for all 47 profiles and removes nothing. The
output was computed after that commit was public.

Files added (4):

| File | Contents |
|---|---|
| `second/WITHOUT_SECOND_SOLVER_PLAN_2026-11-15.md` | the plan of the recomputation; in `main` since commit `18b7168e` (SHA-256 `abcca74bc22efd402ddf2440ce24788e4f6e652fa63ea4f0554ef522484caee1`) |
| `second/without_second_solver.py` | the script; in `main` since commit `18b7168e` (SHA-256 `48e93d62c2bc87e29a96cd953915954d639ee690037b300066ed5c2f724067d2`) |
| `results/second/without_second_solver.json` | its output: every reported quantity for all 47 reserved profiles and for the 43 that remain, with the three readings of the plan of the second reconstruction (SHA-256 `ed674cedc836179c4146207fd31ae319a13718bdb3758cb5272e343e6ccdb740`) |
| `second/tests/test_without_second_solver.py` | 8 checks of that output against the deposited records; `second/tests/` now holds 13 tests |

Files changed (6). The digests before are those of the main branch at commit `7cf0b92`, that is,
release v1.3.0 followed by the commit that recorded its archived identifier.

| File | Digest before | Digest in v1.3.1 | Change |
|---|---|---|---|
| `README.md` | `d928d9e0a209e996de07a3c28e9f99663a8d6958152f3162e7ac0e899b2a14d3` | `da98b88654170b8cdeb2d35ce968a10157a3c5ab78cd93c436a3116e0bacf5c4` | one sentence in the section on the second reconstruction and the two Layout rows of `second/` and `results/second/` |
| `second/README.md` | `bbff286dd97258914db89f409d0aaf7d25ac9a519e62e14469026aa4b3bd1b1a` | `3f82e4f29e2acbbadd4831a74d1c330ff9670f6c0834acc5166814f8dd1227d8` | two rows in the file table, one row in the results table, one paragraph under the commands and one under Results |
| `RELEASE_MANIFEST.md` | `ee850361eb7dc81b1537b445aa9d0e63bd88811d1e55fa8ae665d00b1896bde5` | see `CHECKSUMS.sha256` | this file: tag, archived-record line, opening paragraph, section 7 and this section |
| `CITATION.cff` | `4779a4d161edbe92f7b09dad8e410d598126460f542e954a5313b269138efaf8` | `bad98372bef9deaea474c878639fb10b573827f25d1f182fbe5f0d41ad2ad1c5` | `version` 1.3.1; no new identifier, because the identifier of this release does not exist before it is archived |
| `PROVENANCE_NOTE.md` | `b68687ed539cd08df5eb8b776f4bfb3ff3f1b220cb8ad6e7a3e3e1639ecdb53a` | `fee0bd4a1195af0697c91caae258fc4ce7c97bd1de6c17fabb2f1c7c1b6aa827` | a "Version 1.3.1" paragraph |
| `verify/verify_manuscript_numbers.py` | `0259da88311cda53e2f48256fe950e04b58e69b42f37a20185d62d55e8a862d3` | `20a6076644e8225ed874d50a9e9c80ba431b4288d306dc6367adc5d391ff3cfb` | 20 checks added, 349 in all: 17 for the recomputation and 3 for sentences of the article on the second reconstruction; the claim string of one check follows the article's abstract |

Checks repeated on a fresh clone of the commit to be tagged, before the tag was pushed:
`verify/check_release.py` (every listed file matches; the lock reconciliation as in section 1); the
third asset through `verify/fetch_release_assets.sh second`; the number checker with the article's
text (349 checks, 0 problems); `pytest second/tests` (13 passed); the `build` and `evaluate` stages
of `second/reproduce_second.sh` in a rebuilt second root (the build records and the ten outputs
byte-identical to the deposited files); `second/without_second_solver.py` in that rebuilt root (its
output byte-identical to `results/second/without_second_solver.json`); and the identifier audit.
