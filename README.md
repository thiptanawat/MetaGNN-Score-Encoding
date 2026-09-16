# Objective allowances and growth assumptions shape transcript-informed exchange predictions in a human metabolic model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22745338.svg)](https://doi.org/10.5281/zenodo.22745338)

Code, protocol records, derived data and result tables for a controlled study of the interface
between transcript-derived reaction scores and constraint-based prediction of metabolic exchange.
The repository is laid out as the study root that every script expects, so the locked pipeline and
the analyses added after the lock run from a checkout without any path mapping.

The question is narrow on purpose. Two encodings can rank every supported reaction identically and
still assign different distances between their scores. An optimizer using those numbers as relative
costs can then select different fluxes even though the network, the growth task and the evidence
ordering are unchanged. The study measures how often that reaches the reported exchanges and their
agreement with measurement, and then asks which assumptions restrict the reported exchange
coordinates: the network, the objectives with their allowances, or the transcript costs.

## What the study found

Five power encodings and a midrank encoding were compared in a fixed Recon3D weighted-parsimony
formulation, holding the supported reaction set and its computed ordering constant, and scored
against the published NCI-60 consumption and release measurements for 52 metabolites in 47
reserved profiles from 44 origins. The power encodings changed 246 of 2,444 predictions and
reversed 1,884 of 56,056 between-origin orderings, but no encoding contrast in concordance was
distinguishable from zero (midrank minus grid -0.0013, interval -0.0087 to 0.0069), and 38 targets
were constant across profiles under every encoding.

Nested admissible sets separate what the network, medium and task determine from what the
objectives determine. The network fixed 3 targets. The near-exact cost cap then fixed 90 to 92
percent of the remaining coordinates, a uniform cost with no transcript information fixed 87.8
percent of the same coordinates, and every target that was constant across profiles took the
uniform solution's value; in a typical profile the transcript costs moved 5 to 8 of the 52 targets
away from it, and those few carry a substantial share of the boundary flux. Relaxing the cost cap
alone widened the ranges in proportion to the allowance over the tested grid and cut the orderings
the interval rule certified from 10 to 14 percent of pairs to none at 1 percent of the optimum;
with the parsimony cap kept, the complete two-cap readout still returned a determinate point
(concordance 0.498 to 0.505). Rescoring each saved flux vector under every other encoding's cost
shows that the encodings select different interior solutions, differing in a median of 117 of
10,600 reactions, and that the difference mostly does not survive the projection onto the reported
exchanges.

Two formulations that carry the evidence differently were run through the same evaluation.
Writing the evidence into reaction bounds let more targets vary between profiles when the growth
task followed each profile, and matched the capped cost when the task was common to all profiles
(0.508 and 0.498); a published pruning-and-sampling workflow, run as published, scored 0.494. A
second linear-programming solver reproduced the fixed classification of all 1,536 reported
coordinates it was given, with every point within the classification threshold. On independent
lung-fibroblast condition-response data, a rule that reports a direction only when every encoding
agrees abstained on all 8 declared contrasts under the original condition-independent growth task,
on 7 of 8 under a condition-independent task at the mean measured growth rate, and reported all 8
once measured growth rates were supplied, ordered exactly as those rates.

The conclusion the repository supports is a joint one: near-exact objective caps restricted most
reported exchange coordinates, often to values a transcript-free cost also returns; score spacing
still altered some exchanges and the internal flux distributions; relaxing the allowances widened
the ranges and reduced certified orderings; and on independent data the predictions followed the
imposed growth task. Objective allowances, admissible ranges, discriminatory coverage and measured
agreement should be reported together with selected exchange predictions.

## Layout

The repository is the study root: `reproduce.sh`, `configs/`, `data/`, `manifests/`, `protocol/`,
`src/`, `tests/`, `reports/` and `results/development_v3`, `results/reserved_v3` are the locked
study exactly as `protocol/PROTOCOL_LOCK.json` names them (see `RELEASE_MANIFEST.md` for the
reconciliation). The directories added after the lock sit beside them.

| Path | Contents |
|---|---|
| `reproduce.sh` | The locked study driver, byte for byte the file named in the protocol lock: preflight, data, analysis, tests, controls, run, audit, evaluate, summary |
| `reproduce_extensions.sh` | The analyses added after the lock, in dependency order: coverage, mechanism, closure, independent, figures, verify |
| `configs/study.json` | The frozen study configuration (encodings, solver tolerances, objective allowances, thresholds, scenarios) |
| `protocol/` | The frozen analysis plan, the protocol lock, the design record and the exposure ledger |
| `manifests/` | Context allocation, chemical mapping, gene maps, coverage tables, source records and the data-build receipt |
| `data/` | The derived inputs the lock names: the reaction-evidence matrix, the expression archive extract, the observations, the noise estimates and the supplement layout; `data/raw/README.md` names the four third-party inputs that are not redistributed |
| `src/` | The study pipeline: model construction, medium, gene-protein-reaction evaluation, readout, evaluation, controls, audit and figures 1 to 3 |
| `tests/` | The locked study's 74-test suite, run by `reproduce.sh --stage tests` |
| `results/development_v3/`, `results/reserved_v3/` | The completed development and reserved runs: predictions, evaluations, transform sensitivities, independent audits, concordance reconstructions and reproduction summaries (per-arm records are release assets; see below) |
| `results/production_controls.json`, `reports/` | The validation evidence the lock names, the hardware record, the fresh-arm reproduction and the validation scripts |
| `coverage/`, `results/coverage/` | The encoding-unanimity coverage audit, its tests and outputs |
| `mechanism/`, `results/mechanism/` | Nested admissible sets on the development pilot, the reserved cohort and the uniform-cost control; the interval-position analysis; the cross-encoding cost comparison; with tests |
| `closure/`, `results/closure/` | The objective-allowance sweeps (single cap, fine grid, complete two-cap readout) and their interval scoring, the supporting statistical analyses, the bound-based interface with its common-task controls, the RIPTiDe run, the second-solver checks and the sensitivity-scenario ladders |
| `copeland/`, `results/copeland/` | The independent condition-response evaluation: the frozen plan and its amendment record, the input acquisition script, the culture medium, the preprocessing interface, the runner, the crossed control, the condition-independent mean task and the reporting rule under both of its versions, with tests |
| `endpoint/`, `results/endpoint/` | A reconstruction attempt of a published depletion screen, reported as unresolved |
| `figures/`, `results/figures/` | The figure script for figures 4 to 6 and the rendered figures 1 to 6 |
| `verify/` | The release checker, the release-asset fetcher and the checker that traces every reported number to its artefact |
| `RELEASE_MANIFEST.md` | What this release carries, what regenerates, what was normalized, and the digests of the release assets |

## Reproducing

Everything runs on a CPU with the GLPK solver; no accelerator is used at any stage. The locked
environment is `environment.lock` (Python 3.12); three post-lock stages need the extra packages in
`environment_extensions.txt`.

```
python3.12 -m venv .venv && .venv/bin/python -m pip install -r environment.lock
python verify/check_release.py               # digests of every deposited file and the protocol lock
bash reproduce.sh --stage tests              # the locked study's 74 tests
bash reproduce_extensions.sh tests           # the post-lock analyses' tests
python verify/verify_manuscript_numbers.py   # every reported number against the deposited artefacts
```

To go beyond the deposited outputs, place the four raw inputs as `data/raw/README.md` describes
and run `bash reproduce.sh --stage preflight`; `python -m src.verify_inputs --root . --check lock`
then confirms the complete lock. `bash reproduce.sh --stage all --partition development` and
`--partition test` regenerate the two study runs into new directories (the reserved run took a
summed 49.87 minutes of batch wall time on 12 workers of the recorded CPU). The per-arm records
of the completed runs, which the post-lock mechanism and closure stages read, are release assets
(`verify/fetch_release_assets.sh`), and `reports/fresh_arm_reproduction_v3/comparison.json`
records that a fresh arm reproduced its saved record byte for byte. `reproduce_extensions.sh`
writes each repeat into `results_repeat/<timestamp>/` so that the deposited outputs stay in place
for comparison; the scripts that keep run manifests refuse to overwrite an existing run.

Run the tests first. They are known-answer checks written before the outputs they guard were
interpreted, and a sign convention read the wrong way round fails there rather than reaching a
conclusion.

## Scope

Feasible ranges here are numerical model objects evaluated under declared objective caps. They
are not confidence intervals, no coverage guarantee is attached to them, and the reporting rule
built on them is not a probability of biological correctness. The constraint-stage attribution is
conditional on the declared objective sequence and on this reconstruction, medium and task: a
width contraction localizes a mathematical restriction, not a regulatory bottleneck. The
independent evaluation covers four conditions, two metabolites and four biological replicates in
one cell type, and two measurement-eligible contrasts cannot establish an accuracy rate.

The bound-based interface follows the idea of E-Flux without reproducing that method, and RIPTiDe
was run as published without tuning; no comparative performance claim is made about either
published method. The weighted-parsimony objective of the main study is not a replication of any
published method.

## License

Code is released under the MIT License. Result tables, derived data and protocol records in this
repository are released under CC-BY-4.0. Third-party source data retain the terms of their
original resources and are not redistributed here.
