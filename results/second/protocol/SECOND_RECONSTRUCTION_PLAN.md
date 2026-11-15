# Plan: the encoding, attribution and allowance analyses on a second reconstruction (Human-GEM)

Date: 14 November 2026. This plan was fixed, and its SHA-256 recorded in the lock of the second root (`protocol/PROTOCOL_LOCK.json`, written by `second/run_arms.py --write-lock`), before any reserved profile was run on the second reconstruction and before any measured exchange was read in this analysis. Nothing below was changed afterwards; deviations, if any, are listed at the end under a dated heading.

## 1. Question

The main analysis uses one human reconstruction, Recon3D. Three of its statements could depend on that choice:

1. **Encoding.** Changing the spacing of the reaction scores while keeping their ranking changes some exchange predictions, and no contrast in concordance between encodings is distinguishable from zero.
2. **Attribution.** The near-exact cap on the transcript-weighted cost fixes most of the exchange coordinates that the network leaves open, often at the value that a cost carrying no transcript information also returns.
3. **Allowance.** Relaxing the cost cap widens the admissible ranges and reduces the between-origin orderings that the interval rule can certify.

This plan asks whether each statement holds when the reconstruction is Human-GEM and everything else is held.

## 2. What is held

| Item | As in the main analysis |
|---|---|
| Profiles | the 59 retained expression profiles; the same allocation: 11 development profiles, 47 reserved profiles from 44 origins (`manifests/contexts.tsv`) |
| Targets | the 96 chemistry exchanges handed to the readout and the 52 primary targets of the locked analysis plan (`protocol/analysis_plan.json`) |
| Measurements, pairs, tolerances | the published exchange measurements, the pair policy, the assay tolerances, the tie and interval tolerances and the bootstrap (2,000 origin resamples, the seed of the locked plan), through `src/evaluation.py` |
| Encodings and costs | five power encodings and the midrank encoding of the reaction evidence, the cost construction and the uniform cost of `src/cpu_model.py`; the scale constant is the median of the positive development evidence, recomputed for this reconstruction from the development profiles only |
| Readout | `src/cpu_readout.py`, unchanged: the primary cost cap, the parsimony cap, the ranges before any coordinate is fixed, the lexicographic point |
| Allowances and tolerances | `configs/study.json`, unchanged apart from the objective identifier and the version label |
| Scenario | the primary scenario (task fraction 0.9, serum capacity 1.0) |
| Code | `src/`, `mechanism/` and `closure/` of the main analysis, run in a second root; `second/build_root.py` reports which files are byte-identical |

## 3. What differs, and by which rule

| Item | Rule |
|---|---|
| Model | Human-GEM 2.1.1 (`SysBioChalmers/Human-GEM`, commit `6812be2a5c1504a89ea9e84d770786a84f5f3ce2`), read from `model/Human-GEM.xml` (SHA-256 `dc0efd8009a8b680bde7…`, full value in `second/build_report.json`). Stoichiometry, bounds, gene rules and the biomass reaction are untouched. |
| Exchange identifiers | The pipeline finds its targets and its medium by identifier. Every boundary reaction of Human-GEM is named `EX_<metBiGGID>_e`, where `metBiGGID` is the BiGG identifier that Human-GEM's own metabolite table gives for the exchanged metabolite; where the metabolite has none, the reaction's own BiGG identifier; a boundary reaction takes such an identifier only when no other boundary reaction would, and `EX_<its Human-GEM identifier>` otherwise. All 96 chemistry targets and all 52 primary targets are found this way. For 95 of the 96 the reaction-level cross-reference gives the same identifier (ornithine is `EX_orn__L_e` there). For glyceraldehyde the KEGG identifier of Human-GEM's metabolite (C02154) differs from that of the measured analyte (C00577); it is kept, because Human-GEM itself cross-references the metabolite to the BiGG identifier of the main model. |
| Medium | `src/medium_v2.py`, unchanged, opens the list of the main analysis. Five of its stems have no such exchange in Human-GEM and are recorded as missing, as that module does: vitamin B12 as an unrepresented source (as in the main model), 4-aminobenzoate, nitrate, zinc and manganese. The biomass reaction of Human-GEM cannot carry flux on that list. The inputs of Human-GEM's own growth task ("Growth on Ham's media") that the list does not open were closed one at a time in the order of the task file, and an input stayed closed if biomass could still be produced. Six had to stay open and are added at a capacity of 1.0, the capacity of the assumed serum constituents: aquacob(III)alamin, lipoic acid, linolenate, retinoate, alpha-tocopherol and gamma-tocopherol. Alanine, cysteine, hypoxanthine and thymidine were not needed and stay closed. |
| Objective | Human-GEM's generic human cell biomass reaction (`MAR13082`). |
| Reaction evidence | The gene-level matrix of the main analysis (Entrez identifiers) mapped to Human-GEM's genes through Human-GEM's gene table; the gene rules evaluated by `src/gpr.py` with the conservative policy; the supported set defined as in the main analysis. 2,187 of 2,870 genes are measured and 7,107 of 7,746 internal reactions with a rule are supported. |
| Solver, second pass | The solver of the main analysis (GLPK, tolerance 1e-9) is used for every arm and every range. In 2 of the 67 development arms on Human-GEM its lexicographic point selection returned an infeasible problem; the ranges were computed in every arm. An arm whose readout fails is solved again, by the same worker, on the second solver of the main analysis's cross-solver check (HiGHS through optlang's hybrid interface) at the same tolerances (`second/run_arms.py`). The arms solved this way are listed with the results. |

## 4. What was looked at before this plan was fixed

The 11 development profiles were run on the second root to make the pipeline work: the growth check of Section 3, the failures of the point selection and their repair by the second pass, the time per arm, and the widths of the admissible ranges (in the 65 development arms that GLPK solved, 9 of 6,240 coordinate cells were wider than 1e-5). On one development arm that both solvers solved, the points differed by at most 9.3e-8 and the fixed classification agreed on all 96 coordinates.

The development pilot of the range stages of the main analysis (`mechanism/nested_range_pilot.py`) was also started on the second root. Its time budget of 480 seconds is fixed in the script and was set for the main model; on Human-GEM it allowed 8 of its 22 jobs (4 development profiles under two encodings). In those 8 jobs the network alone fixed 2 of the 52 targets, and the primary cost cap fixed 377 of the 400 coordinate cells that the network leaves open. The pilot reads no measurement. Its ledger is an input of `mechanism/range_location.py` and is reported as partial; no statement of Section 6 rests on it.

No measured exchange was read for any profile in this analysis, development or reserved, before this plan was fixed, and no concordance was computed on Human-GEM. No reserved profile was run.

## 5. Analyses

All scripts are those of the main analysis, run with the second root as their root.

| Step | Script | Output in the second root |
|---|---|---|
| A | `second/run_arms.py --partition test` (the job list and worker of `src/run_study.py`) | `results/reserved_v3/`: 47 profiles by six encodings and the uniform cost; points and ranges of the 96 exchanges |
| B | `src/transform_sensitivity.py` | predictions changed by the encoding, between-origin orderings reversed |
| C | `src/evaluation.py`, primary scenario, reserved partition | concordance of every encoding on the 52 primary targets; the contrasts midrank minus the mean of the power encodings and midrank minus identity, with their 95% bootstrap intervals |
| D | `mechanism/reserved_range_stage.py` (with a time budget of 14,400 seconds in place of the default, an argument of the script), `mechanism/uniform_control.py`, `mechanism/range_location.py`, `second/summarize.py` | coordinates fixed by the network alone, after the primary cost cap and after the parsimony cap; the same under the uniform cost; targets constant across profiles (largest minus smallest point over the 47 profiles at most 1e-5), how many of them lie at the uniform solution's value (within 1e-5 in every profile), and how many targets the transcript cost moves from that value in a profile |
| E | `closure/cap_sweep.py` (caps of 1%, 5% and 20% beside the near-exact cap; identity, midrank and uniform costs) and `closure/cap_concordance.py` | share of fixed coordinates and of orderings certified by the interval rule under each cap |
| F | `closure/solver_check.py --solver hybrid --panel chemistry` | agreement of the two solvers on 16 reserved arms of this reconstruction |

Steps A, B, D, the ranges of E and F read no measurement. Step C and the concordance of E read the measurements, after step A is complete.

## 6. Reading

Every quantity of Section 5 is reported beside its value in the main analysis, whatever it is. For each statement of Section 1 the following is fixed in advance.

1. **Encoding.** The statement holds on Human-GEM if the 95% bootstrap interval of each of the two concordance contrasts includes zero. If an interval excludes zero, a difference between encodings is reported for this reconstruction, with its sign.
2. **Attribution.** The statement holds on Human-GEM if, under every encoding, the primary cost cap fixes more than half of the coordinates that the network leaves open, and if more than half of the targets that are constant across profiles take the value of the uniform solution.
3. **Allowance.** The statement holds on Human-GEM if, for the identity and the midrank encoding, the share of between-origin orderings certified by the interval rule is lower under the 1% cap than under the near-exact cap, and the share of fixed coordinates is lower.

The manuscript reports for each statement whether it held, and reports the numbers in either case.

## 7. What is not repeated

The bound-based interface, the published pruning-and-sampling workflow, the two sensitivity scenarios and the independent fibroblast data are not repeated on Human-GEM. They remain results on one reconstruction, and the manuscript says so.

## 8. Failures

An arm that still fails after the second pass is reported as failed. If a primary arm is missing, the concordance of the primary panel is not estimable under the rule of the main analysis, and that is what is reported. No arm is rerun with changed settings, and no setting of Section 2 is changed after this plan is fixed.

## 9. Deviations

None at the time the plan was fixed.
