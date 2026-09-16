# Fresh-directory single-arm reproduction

Status: **PASS**. Run on 13 September 2026.

The isolated run reproduced **CO:HT29 / magnitude_g1.00 / primary**, reporting **96 exchanges** and a full **10,600-reaction flux vector**. It used the fresh data build under `work/reproduction_clean`, with only the final `cpu_model.py`, `cpu_readout.py`, `medium_v2.py`, configuration and environment lock copied as production implementation inputs. A reviewer-written driver invoked the model/readout directly; it did not invoke the full-study runner.

The isolated result was saved at **2026-09-13T16:23:20.973089+00:00**. The main result was first read after this, at **2026-09-13T16:23:36.627055+00:00**. The isolated calculation took **33.63 seconds**. No CORE outcome comparison was performed.

| Comparison | Maximum absolute difference |
|---|---:|
| Full feasible flux vector | 0 |
| Reaction cost weights | 0 |
| 96 reported exchange points | 0 |
| Exchange range endpoints | 0 |
| Primary/secondary optima and caps | 0 |
| Saved task maximum and fixed demand | 0 |

All copied input/source files match the current main files. The main arm's source fingerprint is current. The compressed flux archives are byte-identical.

An independent saved-vector audit passed mass balance, reconstructed bounds, fixed task, both cost caps, point/range correspondence, independently reconstructed weights and final logical LP statuses. Its detailed residuals and any recovered attempts are in `comparison.json`.

**Scope:** both runs used the same existing `.venv` Python environment. This checks one fresh-directory model execution with fresh-derived inputs; it is not a separate environment installation, full cohort reproduction or biological validation. Objective/FVA extrema were compared between two solves, while the algebraic audit itself did not independently optimize them.

Evidence: `copy_manifest.json`, `isolated/fresh_driver.py`, `isolated/fresh_execution.log`, `isolated/fresh_arm.json`, `isolated/fresh_arm.npz`, `compare_after_production.py`, and `comparison.json`.
