# Second reconstruction: three analyses repeated on Human-GEM 2.1.1

The main analysis of this repository ran on one human reconstruction, Recon3D. This directory
repeats three of its analyses on a second reconstruction, Human-GEM 2.1.1, to ask whether three
statements depend on that choice: that changing the spacing of reaction scores while keeping their
ranking changes some exchange predictions without moving concordance, that the near-exact cost cap
fixes most coordinates the network leaves open, and that relaxing the cap widens the ranges and
reduces the orderings the interval rule certifies. The design, the held settings and the rule for
reading each statement were fixed in `SECOND_RECONSTRUCTION_PLAN_2026-11-14.md` (SHA-256
`0de8cbf4a737829b4c005b33122568d1abb342166162917f5d5481d4b7cd9a81`). The plan is not edited here.

The plan was hash-locked before the reserved runs. Its digest was written into the lock of the
second root (`results/second/protocol/PROTOCOL_LOCK.json`, field `locked_at_utc`
`2026-11-14T16:07:03Z`), and the driver does not run the reserved partition without that lock. The
plan was not public before the runs. It is released with the results, so a reader can check the
hash and cannot check the order from a public record. The main analysis, its lock and its results
are not changed by anything in this directory.

## What is in `second/`

| Path | Contents |
|---|---|
| `SECOND_RECONSTRUCTION_PLAN_2026-11-14.md` | The locked plan, byte for byte |
| `build_root.py` | Builds the second root from this repository and a Human-GEM checkout at the pinned commit |
| `run_arms.py`, `in_second_root.py`, `run_ranges.sh`, `summarize.py` | The arm runs with the second-solver pass, the wrapper that runs a main-analysis module inside the second root, the range stages with the cap sweep and the solver check, and the summary of constant targets |
| `reproduce_second.sh` | A driver written for this release that chains the calls below: stages `build`, `development`, `pilot`, `lock`, `reserved`, `ranges` and `evaluate` |
| `build_report.json`, `exchange_map.tsv`, `target_map.tsv`, `medium_additions.tsv` | The build records: model size, the SHA-256 of the five Human-GEM files read, the exchange and target identifier maps, and the six growth-task inputs added to the medium |
| `environment_second.txt` | The package list of the environment the runs were made in |
| `tests/test_second_summarize.py` | Five known-answer tests of `summarize.py` on a synthetic root. They are kept here so that `tests/` stays the locked study's 74 tests; `bash reproduce_extensions.sh tests` runs them |
| `WITHOUT_SECOND_SOLVER_PLAN_2026-11-15.md`, `without_second_solver.py` | Added in release v1.3.1: the plan and the script of the recomputation without the 4 profiles in which the second solver completed an arm. Both were committed to `main` before the script was run without `--check-only` |
| `tests/test_without_second_solver.py` | Eight checks of `results/second/without_second_solver.json` against the deposited records. They need a checkout only and run with the tests above |
| `audit/` | The comparison of the Human-GEM identifiers with those of Recon3D by formula, charge and InChIKey, made after the runs (`audit/README.md`). It was not part of the plan and changed nothing in the analysis. The classes are 88 of 96 target exchanges in full agreement, 3 differing in protonation only, 3 with no stereochemistry on one side, 2 with a different stereochemical block and none a different compound. The scripts need RDKit, which no other part of the repository uses; run from this checkout with RDKit 2026.09.2, they wrote the six tables and the report byte for byte |

## What is in `results/second/`

| Path | Contents |
|---|---|
| `protocol/PROTOCOL_LOCK.json`, `protocol/SECOND_RECONSTRUCTION_PLAN.md` | The lock of the second root and the plan copy it names, so the lock can be checked from a checkout |
| `reserved_v3/` | Predictions of the 47 reserved profiles, the evaluation (`evaluation_primary.json`, `evaluation_primary_summary.tsv`, `evaluation_primary_per_metabolite.tsv`), the transform sensitivity (`transform_primary.json`), the record of the second-solver pass (`second_solver_arms.json`) and the run scalars |
| `development_v3/` | Predictions of the 11 development profiles, the record of the second-solver pass and the run scalars |
| `second_summary.json` | Targets constant across profiles and their relation to the uniform solution |
| `without_second_solver.json` | Release v1.3.1: every reported quantity for all 47 reserved profiles and without the 4 profiles in which the second solver completed an arm, with the three readings of the plan |
| `mechanism/` | The range stages: the reserved cohort (`reserved_range/`: ledger, summary, run manifest), the uniform-cost control (`uniform_control/`), the interval-position analysis (`range_location*`) and the partial development pilot that the plan describes (`bottleneck_pilot/`: ledger, summary, run manifest and the network-only record of the 8 calculations that fit the script's time budget) |
| `closure/cap_sweep/` | The cap sweep (`summary.json`, `cap_concordance.json`, `run_manifest.json`, one job record per profile and arm in `jobs/`) |
| `closure/solver_check_96/` | The check of GLPK against HiGHS over the production list of 96 exchanges |

The per-arm records and flux arrays of the second root and its two `run_summary.json` files are
one release asset (`second_reconstruction_arms.tar.gz`; its digest and contents are in
`RELEASE_MANIFEST.md`, and `bash verify/fetch_release_assets.sh second` downloads, verifies and
unpacks it). They are not in the tree. Human-GEM is not redistributed.

## Inputs to fetch

- Human-GEM at tag v2.1.1, commit `6812be2a5c1504a89ea9e84d770786a84f5f3ce2`
  (`SysBioChalmers/Human-GEM`, CC BY 4.0). The SHA-256 of the five files read are in
  `build_report.json` and in `SOURCES.md`.
- The main repository's own inputs, as its `README.md` describes. `build_root.py` reads only
  derived files in the checkout and no file under `data/raw/`. The two measurement tables are copied
  into the second root only by the `evaluate` stage.

## Environment

The runs were made on a Linux x86-64 workstation (Ubuntu 24.04) under Python 3.12.3 with COBRApy
0.32.1, optlang 1.9.1, GLPK 5.0 (swiglpk 5.0.13), HiGHS (highspy 1.15.1), NumPy 2.4.4 and SciPy
1.17.1. The full package list is `environment_second.txt`. To rebuild an environment, start from
`environment.lock` (Python 3.12) and add `highspy` and `osqp` at the versions in
`environment_extensions.txt`. HiGHS is used for the second solver pass and the solver check. No GPU
is used.

## Commands, in order

`reproduce_second.sh` chains them in the order of the recorded runs (stages `build`, `development`,
`pilot`, `lock`, `reserved`, `ranges`, `evaluate`). It stops if the output it would write exists.
From the repository root:

```
bash second/reproduce_second.sh build       --humangem HUMAN_GEM_CHECKOUT --root SECOND_ROOT
bash second/reproduce_second.sh development --root SECOND_ROOT --workers N
bash second/reproduce_second.sh pilot       --root SECOND_ROOT --workers N
bash second/reproduce_second.sh lock        --root SECOND_ROOT
bash second/reproduce_second.sh reserved    --root SECOND_ROOT --workers N
bash second/reproduce_second.sh ranges      --root SECOND_ROOT --workers N
bash second/reproduce_second.sh evaluate    --root SECOND_ROOT
python -m pytest -q second/tests
```

The pilot is the development pilot of the range stages of the main analysis. Its ledger is read by
the interval-position step of `evaluate` and by nothing else, and no statement of the plan rests on
it. It stops at the time budget fixed in `mechanism/nested_range_pilot.py` (480 seconds), which was
set for Recon3D. The number of calculations it completes on Human-GEM therefore depends on the
machine: the recorded pilot completed 8 of 22, and a repeat on other hardware can complete more or
fewer, which changes `range_location_pilot.tsv` and the development part of `range_location.json`
and nothing else.

The calls underneath, for a reader who wants them one at a time:

```
python second/build_root.py --main . --humangem HUMAN_GEM --out SECOND_ROOT
python second/run_arms.py --root SECOND_ROOT --partition development --output results/development_v3 --workers N
bash   second/run_ranges.sh SECOND_ROOT PYTHON N pilot
python second/run_arms.py --root SECOND_ROOT --write-lock second/SECOND_RECONSTRUCTION_PLAN_2026-11-14.md
python second/run_arms.py --root SECOND_ROOT --partition test --output results/reserved_v3 --workers N
bash   second/run_ranges.sh SECOND_ROOT PYTHON N reserved
python second/in_second_root.py SECOND_ROOT src.evaluation --root SECOND_ROOT --partition test ...
python second/in_second_root.py SECOND_ROOT src.transform_sensitivity --root SECOND_ROOT --partition test ...
python second/summarize.py --root SECOND_ROOT --partition test --predictions results/reserved_v3/predictions.tsv --output results/second_summary.json
python mechanism/range_location.py        # run with STUDY_ROOT=SECOND_ROOT inside SECOND_ROOT/results/mechanism
python closure/cap_concordance.py --root SECOND_ROOT --sweep SECOND_ROOT/results/closure/cap_sweep --out SECOND_ROOT/results/closure/cap_sweep
```

The `build` stage compares its four build records with the copies in this directory and prints
`identical` or `DIFFERS` for each. The deposited numbers are checked against the files of
`results/second/` by `python verify/verify_manuscript_numbers.py`.

The recomputation without the profiles that the second solver completed is not a stage of the
driver. After `evaluate`, `python second/without_second_solver.py --root SECOND_ROOT --check-only`
computes every reported quantity for all 47 profiles, compares it with the outputs of `evaluate`
and writes nothing. The same command without `--check-only` then removes the 4 profiles and writes
`SECOND_ROOT/results/without_second_solver.json`, which is deposited as
`results/second/without_second_solver.json`.

## What is held and what differs

| | Same as the main analysis | Differs, by the rule in the plan |
|---|---|---|
| Data | 59 expression profiles (11 development, 47 reserved from 44 origins), 96 chemistry exchanges, 52 primary targets, measurements, pair policy, tolerances, bootstrap | none |
| Code | `src/`, `mechanism/`, `closure/`, run unchanged in a second root | only the medium module gains six growth-task inputs, appended to a copy |
| Settings | `configs/study.json`, the scenario, the six encodings, the uniform cost | the objective identifier and the version label |
| Model | none | Human-GEM 2.1.1 in place of Recon3D; exchanges renamed from Human-GEM's own BiGG cross-references; its biomass reaction (`MAR13082`) as the objective |
| Evidence | the gene-level matrix, the gene-rule policy | mapped to Human-GEM genes through its gene table |
| Solver | GLPK at the tolerances of the main analysis | an arm whose point selection fails is solved again with HiGHS (`run_arms.py`) |

## Run times

Wall times of the reserved runs that produced `results/second/`. The arms, the range stage and the
cap sweep ran with 16 worker processes, the second-solver check with 8 (the value `run_ranges.sh`
passes), and the uniform-cost control as one process.

| Stage | Workers | Wall time | Recorded in (under `results/second/`) |
|---|---:|---:|---|
| reserved arms, first pass (GLPK) | 16 | 837 s | the run log, which is not deposited |
| reserved arms, second pass (the 4 arms that GLPK did not solve, on HiGHS) | 16 | 674 s | `reserved_v3/run_summary_scalars.json` (`seconds`) |
| range stage | 16 | 3,858 s | `mechanism/reserved_range/summary.json` (`seconds`) |
| uniform-cost control | 1 | 272 s | `mechanism/uniform_control/summary_uniform.json` (`seconds`) |
| cap sweep | 16 | 357 s | `closure/cap_sweep/summary.json` (`_run.seconds`) |
| second-solver check | 8 | 1,348 s | `closure/solver_check_96/solver_check.json` (`seconds`) |

The worker counts are the arguments the driver was given; the `workers` entry of
`configs/study.json` (4), which the run summaries copy, is the default those arguments replace.
Each worker process solves whole jobs on its own copy of the model, so the count sets the wall time
and not the values.

## Results

Each statement is read by the rule fixed in the plan. An interval that includes zero is not
evidence of equivalence, and the plan sets no equivalence margin. All three statements hold on
Human-GEM. Table S6 of the supplementary article gives every quantity beside its value on Recon3D;
the Recon3D column below is copied from it.

**1. Encoding.** Rule: the statement holds if the 95 percent origin-bootstrap interval of each of
the two concordance contrasts includes zero.

| Quantity | Recon3D | Human-GEM 2.1.1 | File under `results/second/` |
|---|---|---|---|
| Predictions changed by the power encodings | 246 of 2,444 (10.07 percent) | 249 of 2,444 (10.19 percent) | `reserved_v3/transform_primary.json` |
| Between-origin orderings reversed | 1,884 of 56,056 (3.36 percent) | 2,092 of 56,056 (3.73 percent) | `reserved_v3/transform_primary.json` |
| Concordance, power encodings | 0.5045 to 0.5085 | 0.4973 to 0.5018 | `reserved_v3/evaluation_primary.json`, `evaluation_primary_summary.tsv` |
| Concordance, midrank encoding | 0.5054 | 0.5013 | `reserved_v3/evaluation_primary.json` |
| Midrank minus mean of the power encodings | -0.0013 (-0.0087 to 0.0069) | 0.0017 (-0.0031 to 0.0063) | `reserved_v3/evaluation_primary.json` |
| Midrank minus identity | -0.0023 (-0.0094 to 0.0057) | 0.0005 (-0.0054 to 0.0059) | `reserved_v3/evaluation_primary.json` |

Outcome: both intervals include zero, with 2,000 of 2,000 origin resamples valid, over 46,447
eligible pairs.

**2. Attribution.** Rule: the statement holds if, under every encoding, the cost cap fixes more
than half of the coordinates that the network leaves open, and if more than half of the targets
that are constant across profiles take the value of the uniform solution.

| Quantity | Recon3D | Human-GEM 2.1.1 | File under `results/second/` |
|---|---|---|---|
| Targets fixed by the network alone | 3 of 52 | 2 of 52 | `mechanism/reserved_range/summary.json`, `range_stage_ledger.tsv` |
| Open cells first fixed by the cost cap, range over the six encodings | 2,077 to 2,110 of 2,303 (90.2 to 91.6 percent) | 2,153 to 2,223 of 2,350 (91.6 to 94.6 percent) | `mechanism/reserved_range/range_stage_ledger.tsv` |
| Free targets first fixed by the cap under the uniform cost | 43 of 49 (87.8 percent) | 45 of 50 (90.0 percent) | `mechanism/uniform_control/summary_uniform.json` |
| Targets constant across profiles under every encoding | 38 of 52 | 46 of 52 | `second_summary.json` |
| Of those, at the value of the uniform solution | 38 of 38 | 46 of 46 | `second_summary.json` |
| Targets moved from the uniform value in a profile, median by encoding (range) | 5 to 8 (3 to 10) | 5 to 6 (4 to 6) | `second_summary.json` |

Outcome: the cap first fixed more than 90 percent of the open cells under every encoding, and
every constant target took the uniform value. The per-encoding counts of Table S6 are in the
ledger and the summary.

**3. Allowance.** Rule: the statement holds if, for the identity and the midrank encoding, the
share of between-origin orderings certified by the interval rule and the share of fixed
coordinates are both lower under the cap of 1 percent than under the near-exact cap.

| Quantity (identity and midrank) | Recon3D | Human-GEM 2.1.1 | File under `results/second/` |
|---|---|---|---|
| Coordinates fixed of 52, near-exact cap | 47.89 and 47.19 | 48.23 and 49.30 | `closure/cap_sweep/summary.json` |
| Coordinates fixed of 52, cap of 1 percent | 3.00 and 3.00 | 2.00 and 2.00 | `closure/cap_sweep/summary.json` |
| Orderings certified, near-exact cap | 4,725 and 6,395 (10.17 and 13.77 percent) | 3,906 and 3,506 (8.41 and 7.55 percent) | `closure/cap_sweep/cap_concordance.json` |
| Orderings certified, cap of 1 percent | none | none | `closure/cap_sweep/cap_concordance.json` |

Outcome: both shares are lower under the cap of 1 percent for both encodings. At caps of 1, 5 and
20 percent the only fixed coordinates are the 2 that the network fixes, and no ordering is
certified.

**Solver passes.** GLPK's point selection returned an infeasible problem in 4 of the 283 reserved
arms. Each was solved again on HiGHS, and no arm failed on both passes
(`reserved_v3/second_solver_arms.json`; Recon3D: none of 283). The second-solver check over the
production list of 96 exchanges (`closure/solver_check_96/solver_check.json` and
`solver_check_rows.tsv`) covers 16 arms and 1,536 coordinates. Table S6 reports it without the one
arm whose saved point is itself from HiGHS: the fixed classification agreed on 1,440 of 1,440
coordinates (15 arms), the largest absolute difference between points was 1.0e-6, and the
between-profile signs agreed in 4,704 of 4,704 comparisons. `verify/verify_manuscript_numbers.py`
recomputes these three figures from the rows file.

**Without the profiles that the second solver completed (release v1.3.1).** The 4 arms that only
HiGHS solved belong to 4 profiles (OV:OVCAR-4, CO:HCT-15, ME:MDA-MB-435, OV:IGROV1). Under a plan
committed before the computation (`WITHOUT_SECOND_SOLVER_PLAN_2026-11-15.md`), those profiles were
removed under every encoding, and the quantities above were computed again for the remaining 43
profiles (41 origins, 38,971 eligible pairs) from the saved records
(`without_second_solver.json`); no linear program was solved. Midrank minus the mean of the power
encodings was 0.0011 (-0.0041 to 0.0063), against 0.0017 (-0.0031 to 0.0063) for all
profiles, and midrank minus identity -0.0002 (-0.0063 to 0.0059), against 0.0005 (-0.0054 to 0.0059).
The cost cap first fixed 1,973 to 2,032 of the 2,150 open cells per encoding
(91.8 to 94.5 percent). 46 of the 52 targets were constant across the profiles and all
46 lay at the uniform value. The near-exact cap fixed 48.26 and 49.26 coordinates and certified
3,308 and 2,974 of the 38,971 orderings (8.49 and 7.63 percent), and the looser caps
none. By the three rules of the plan each statement holds without these profiles. The intervals are
those of a bootstrap over 41 origins in place of 44.

## Limits and what was not repeated

The reserved profiles and the measured exchanges are those of the main analysis, and the target
panel and the scenario logic were fixed on the development profiles under Recon3D. The comparison
on Human-GEM is therefore the same comparison through another network, and the limitation of the
outcome reference described in the main text applies to it unchanged. The medium is the list of the
main analysis, where Human-GEM represents it, completed by six inputs taken from Human-GEM's own
growth task; it is not a second independent description of the culture medium. The bound-based
interface, the pruning-and-sampling workflow, the two sensitivity scenarios and the independent
fibroblast evaluation were not repeated on Human-GEM and remain results on one reconstruction.
No arm was rerun with changed settings, and there was no deviation from the plan.

Before the plan was fixed, the 11 development profiles were run to make the pipeline work, and a
development pilot of the range stages was started; no measured exchange was read for any profile on
Human-GEM, and no concordance was computed, before the plan was fixed (Supplementary Methods S22).
