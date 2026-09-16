# Release manifest

Release tag: `v1.2.0`
Archived record: see `SOURCES.md` (the identifier of this release is minted when the release is archived, so the copy of that file inside the release names the identifiers of the releases before it; the repository's `main` branch records the new identifier in the commit that follows the tag).

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
| Third-party raw inputs (Recon3D, CellMiner archive, the Jain et al. workbook and supplementary PDF) | no | retrieved once as `data/raw/README.md` describes | `python -m src.verify_inputs --root . --check raw` |
| Independent-evaluation inputs (two public git repositories) | no | `copeland/fetch_inputs.sh` clones them and refuses any commit other than the pinned ones | digests of the two `.rda` files checked by that script |

### Release assets

| Asset | Bytes | SHA-256 | Contents |
|---|---:|---|---|
| `reserved_v3_arms.tar.gz` | 108149951 | `083524997b82bb2e459c2ca078a2c72000566f9ad960b108f11fce23ca91a210` | `arms/` (849 JSON records and 849 NPZ flux arrays of the reserved run: 47 profiles, 6 encodings and 3 scenarios, plus the 3 context-independent uniform-cost records) and `run_summary.json`; unpack into `results/reserved_v3/` |
| `development_v3_arms.tar.gz` | 25567045 | `94709ee80f65a81fed390845c7626622a080ac4c1d3a3634e2ac31dfbc17f351` | `arms/` (201 JSON records and 201 NPZ arrays of the development run) and `run_summary.json`; unpack into `results/development_v3/` |

`bash verify/fetch_release_assets.sh` downloads both from the release that matches this checkout,
verifies these digests and unpacks them.

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

The record of the clean-checkout verification performed on this release is
`verify/clean_checkout_log_v1.2.0.md`.
