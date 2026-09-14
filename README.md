# Where transcript evidence stops influencing predicted metabolic exchange

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22745338.svg)](https://doi.org/10.5281/zenodo.22745338)

Code, protocol records and result tables for a controlled study of the interface between
transcript-derived reaction scores and constraint-based metabolic prediction.

The question is narrow on purpose. Two encodings can rank every supported reaction identically
and still assign different distances between their scores. An optimizer using those numbers as
relative costs can then select different fluxes, even though the network, the task and the
evidence ordering are unchanged. This repository holds the experiment that measures how often
that happens, the diagnosis of which modeling stage decides the answer, and an independent test
of the reporting rule that follows from it.

## What the study found

The encoding intervention changes a minority of predictions and none of the agreement. Rank
preserving power encodings changed 246 of 2,444 reserved primary predictions and reversed 1,884
of 56,056 possible between-origin orderings, yet every encoding contrast against measured
extracellular exchange included zero and every arm sat within about 0.009 of a context
independent reference.

The reason is locatable. Nested admissible sets separate what the network, medium and growth
task determine from what the objectives determine. The network left 49 of 52 exchange targets
free, with a mean admissible width of 35 canonical units. The transcript-weighted cost then
first fixed 90 to 92 percent of those coordinates, and it fixed them at positions that coincided
across profiles: for a median target the profiles moved across about 4 parts in 10 billion of
the width the network had allowed. A uniform cost with no transcript information fixed 87.8
percent of the same targets, so the optimality cap rather than the evidence carried in it
accounts for nearly all of the restriction.

The encodings were not agreeing because they agreed. Rescoring each saved flux vector under
every other encoding's cost shows they select genuinely different interior solutions, differing
in a median of 117 of 10,600 reactions, with 97.4 percent of them violating another encoding's
cost cap. That difference mostly does not survive the projection onto the reported exchanges,
where a median of 6 of 96 coordinates differ.

An independent condition-response evaluation on primary human lung fibroblasts closes the loop.
With the growth task specified as in the original study, identical across conditions, the
reporting rule declined all eight declared decisions. Supplying measured growth rates produced
unanimous directions for all eight, ordered exactly as the growth rates were, correct for one
and incorrect for the other of the two contrasts the measurements resolved.

## Layout

| Path | Contents |
|---|---|
| `protocol/` | The frozen analysis plan, the protocol lock, the design record and the exposure ledger |
| `src/` | The study pipeline: model construction, medium, gene-protein-reaction evaluation, readout, evaluation and controls |
| `manifests/` | Context allocation, chemical mapping, gene maps, coverage tables and source records |
| `mechanism/` | Constraint-stage attribution, the interval-position analysis and the cross-encoding cost comparison, with their tests |
| `copeland/` | The independent condition-response evaluation: the frozen plan, the culture medium, the preprocessing interface, the runner, the crossed control and the reporting rule, with its tests |
| `coverage/` | The encoding-unanimity coverage audit and its tests |
| `endpoint/` | A reconstruction attempt of a published depletion screen, reported as unresolved |
| `results/` | Summary-level outputs: evaluations, contrasts, stage attribution, position analysis, cross-encoding comparison and the independent evaluation |
| `figures/` | The figure script |
| `data/` | Derived outputs: per-arm predictions and ranges, the agreement and nested-range ledgers, the cross-encoding tables and the reaction-evidence matrix |
| `verify/` | A checker that traces every numeric claim in the manuscript back to the artefact that produced it |

Larger derived outputs are under `data/`: per-arm predictions, the full pair ledgers, the
nested-range ledgers, the cross-encoding regret and distance tables and the reaction-evidence
matrix. See `SOURCES.md`.

## Reproducing

`reproduce.sh` names the stages and the order. Three third-party inputs are not redistributed
here and must be retrieved once; `SOURCES.md` gives the retrieval route and the expected
SHA-256 for each. The environment is pinned in `environment.lock`; the solver is GLPK and no
accelerator is required at any stage.

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

No published transcript-to-flux workflow was reimplemented in its native formulation, and no
comparative performance claim is made. The weighted-parsimony objective used here is not a
replication of any published method.

## License

Code is released under the MIT License. Result tables and protocol records in this repository
are released under CC-BY-4.0. Third-party source data retain the terms of their original
resources and are not redistributed here.
