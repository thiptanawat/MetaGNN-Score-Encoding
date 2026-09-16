# An optimality tolerance, not the spacing of transcript-derived scores, selects the predicted metabolic exchanges

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22745338.svg)](https://doi.org/10.5281/zenodo.22745338)

Code, protocol records and result tables for a controlled study of the interface between
transcript-derived reaction scores and constraint-based metabolic prediction.

The question is narrow on purpose. Two encodings can rank every supported reaction identically
and still assign different distances between their scores. An optimizer using those numbers as
relative costs can then select different fluxes, even though the network, the task and the
evidence ordering are unchanged. This repository holds the experiment that measures how often
that happens, the diagnosis of which modeling stage decides the answer, a sweep of the
optimality tolerance that stage depends on, two formulations that carry the evidence
differently, and an independent test of the reporting rule that follows.

## What the study found

The encoding intervention changes a minority of predictions and none of the agreement. Rank
preserving power encodings changed 246 of 2,444 reserved primary predictions and reversed 1,884
of 56,056 possible between-origin orderings, yet every encoding contrast against measured
extracellular exchange included zero, and the design had less than 50 percent power for any
contrast below 0.045 on the 9 targets that varied between profiles.

The reason is locatable, and it is the optimality tolerance. Nested admissible sets separate
what the network, medium and growth task determine from what the objectives determine. The
network left 49 of 52 exchange targets free, with a mean admissible width of 35 canonical
units. The cost cap at the study's near-exact tolerance first fixed 90 to 92 percent of those
coordinates, at positions that coincided across profiles for 41 or 42 of the 49 targets. A
uniform cost with no transcript information fixed 87.8 percent of the same targets, and every
target that was constant across profiles took the uniform solution's value: in a typical
profile the readout differs from the uniform boundary on only 5 to 8 of the 52 targets under
any cost. Loosening the cap released the coordinates in proportion to the tolerance, tenfold in
range width per decade of relative slack, while the share of measured pairs the interval rule
could resolve fell from 10 to 14 percent to none at 1 percent of the optimum.

The encodings were not agreeing because they agreed. Rescoring each saved flux vector under
every other encoding's cost shows they select different interior solutions, differing in a
median of 117 of 10,600 reactions, with 97.4 percent of them violating another encoding's cost
cap. That difference mostly does not survive the projection onto the reported exchanges, where
a median of 6 of 96 coordinates differ.

Two formulations that do not carry the evidence as a capped cost were run through the same
evaluation. Writing the evidence into reaction bounds with the growth task set per profile let
more targets vary between profiles and resolved about twice as many measured pairs, at an
accuracy of 0.48 to 0.50; with the task held at one value across profiles the same interface
pinned 38 of the 49 free targets and scored 0.508 and 0.498, within 0.007 of the cost-based
arms. RIPTiDe, run in its native form, pruned most reported exchanges to zero and scored 0.494.
The complete two-stage readout at loosened tolerances still returned a determinate point, with
concordance between 0.498 and 0.505 at every tolerance tried. A second linear-programming
solver reproduced the fixed classification of every coordinate it was given, and both
sensitivity scenarios reproduced the attribution.

An independent condition-response evaluation on primary human lung fibroblasts closes the loop.
With the growth task specified as in the original study, identical across conditions, or set to
a condition-independent value at the mean measured growth rate, the reporting rule declined all
eight declared decisions. Supplying measured growth rates produced unanimous directions for all
eight, ordered exactly as the growth rates were, with predicted lactate secretion a near-constant
multiple of the imposed rate.

## Layout

| Path | Contents |
|---|---|
| `protocol/` | The frozen analysis plan, the protocol lock, the design record and the exposure ledger |
| `src/` | The study pipeline: model construction, medium, gene-protein-reaction evaluation, readout, evaluation and controls |
| `manifests/` | Context allocation, chemical mapping, gene maps, coverage tables and source records |
| `mechanism/` | Constraint-stage attribution, the interval-position analysis and the cross-encoding cost comparison, with their tests |
| `closure/` | The optimality-tolerance sweeps (single-stage, fine grid and complete two-stage readout) and their interval scoring, the supporting statistical analyses, the bound-based encoding with its fixed-task control, the RIPTiDe run, the second-solver check, the sensitivity-scenario ladders and the comparator evaluator |
| `copeland/` | The independent condition-response evaluation: the frozen plan, the culture medium, the preprocessing interface, the runner, the crossed control, the condition-independent mean task and the reporting rule, with its tests |
| `coverage/` | The encoding-unanimity coverage audit and its tests |
| `endpoint/` | A reconstruction attempt of a published depletion screen, reported as unresolved |
| `results/` | Summary-level outputs: evaluations, contrasts, stage attribution, position analysis, cross-encoding comparison, and under `results/closure/` the tolerance sweeps, the statistical analyses, the comparator formulations, the solver check and the scenario ladders; `results/copeland/` holds the independent evaluation under all three tasks |
| `figures/` | The figure script |
| `data/` | Derived outputs: per-arm predictions and ranges, the agreement and nested-range ledgers, the cross-encoding tables, the reaction-evidence matrix, and under `data/closure/` the per-arm predictions and ledgers of the comparator, two-stage and scenario-ladder runs |
| `verify/` | A checker that traces every numeric claim in the manuscript back to the artefact that produced it |

Larger derived outputs are under `data/`. See `SOURCES.md` for the three third-party inputs
that are not redistributed.

## Reproducing

`reproduce.sh` names the stages and the order. Three third-party inputs are not redistributed
here and must be retrieved once; `SOURCES.md` gives the retrieval route and the expected
SHA-256 for each. The environment is pinned in `environment.lock`; the solver is GLPK and no
accelerator is required at any stage. The `closure` stage additionally needs the `riptide`
package for the RIPTiDe run and `highspy` for the second-solver check; the versions used are
recorded in the run manifests under `results/closure/`.

Run the tests first. They are known-answer checks written before the outputs they guard were
interpreted, and a sign convention read the wrong way round fails there rather than reaching a
conclusion.

```
python -m unittest discover -s mechanism -p 'test_*.py'
python -m unittest discover -s copeland  -p 'test_*.py'
python -m unittest discover -s coverage  -p 'test_*.py'
```

## Scope

Feasible ranges here are numerical model objects evaluated under declared objective caps. They
are not confidence intervals, no coverage guarantee is attached to them, and the reporting rule
built on them is not a probability of biological correctness. The constraint-stage attribution
is conditional on the declared objective sequence and on this reconstruction, medium and task:
a width contraction localizes a mathematical restriction, not a regulatory bottleneck. The
independent evaluation covers four conditions, two metabolites and four biological replicates
in one cell type, and two measurement-eligible contrasts cannot establish an accuracy rate.

The bound-based encoding follows the idea of E-Flux without reproducing that method, and
RIPTiDe was run as published without tuning; no comparative performance claim is made about
either published method. The weighted-parsimony objective used in the main study is not a
replication of any published method.

## License

Code is released under the MIT License. Result tables and protocol records in this repository
are released under CC-BY-4.0. Third-party source data retain the terms of their original
resources and are not redistributed here.
