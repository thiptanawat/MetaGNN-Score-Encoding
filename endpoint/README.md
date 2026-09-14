# Historical depletion-screen reconstruction

**Status: completed reconstruction attempt; source-rule reconciliation remains unresolved.** This is not a completed physiological quality assessment or permission to filter the manuscript's endpoint data.

The [amendment](AUDIT_AMENDMENT.json) was saved and hashed before the new calculations. The [script](run_endpoint_audit.py) reads the original corrected Jain workbook's `Concentrations` sheet, the frozen profile allocation and the previously verified Human1 eleven-line list. It does not read model predictions or calculate accuracy. The original inputs and results remain unchanged.

Nilsson et al. (2018), [Metabolite Depletion Affects Flux Profiling of Cell Lines](https://doi.org/10.1016/j.tibs.2018.03.009), Figure 1B reports 84 depleted and 36 nondepleted cultures under a below-10%-of-fresh criterion. This attempt used the arithmetic mean of ten drift-corrected fresh columns and evaluated the criterion for every available positive-baseline assay. It retained two predeclared variants: all 140 assay rows, and the 115 source-calibrated rows. Uncalibrated rows lack concentration entries in this workbook, so adding them does not add usable quantitative depletion evidence. Three calibrated rows have zero fresh baselines and cannot supply a positive-baseline depletion ratio.

Both variants flagged 119 of 120 cultures, with one unresolved. Every source line had at least one flagged culture. Neither variant reproduces the published 84/36 counts or the historical eleven-line subset. This does **not** show that all 60 lines are biologically invalid. It shows that this literal reconstruction, including its assay eligibility and fresh-baseline interpretation, is insufficient to recover the published screen. The exact source supplementary table/implementation has not been recovered; a general paper statement does not establish every implementation choice.

Do not choose a new assay list or threshold because it recreates a desired count or improves a model result. Determine the original rule from source evidence, then rerun with a separately documented amendment. Matching counts must be followed by matching individual culture identities.

## Files

- [Coverage summary](coverage_summary.json): 120 cultures, 60 source profiles, 58 matched model profiles, both variants and allocation-specific counts; input hashes and limitations.
- [Culture–assay ledger](culture_assay_ledger.tsv): 16,800 source cells, baseline mean, endpoint value, ratio, calibration, spreadsheet coordinates and assessment/missingness.
- [Culture ledger](culture_quality_ledger.tsv): 240 rows, one per culture and reconstruction variant; flags, unresolved-assay counts and source-derived volume.
- [Profile ledger](profile_quality_ledger.tsv): 120 rows, one per source line and variant; both-culture requirement and fixed development/reserved/excluded assignment.
- [Amendment checksum](AUDIT_AMENDMENT.sha256).

The source description supplies culture volumes and a 4- or 5-day window; exact culture durations and complete per-culture growth trajectories were not recovered. No growth-phase correction, uncertainty calibration or quality-filtered prediction evaluation was performed. An independent metabolite-first calculation checked the culture counts against the culture-first ledger.

## Reproduction

Place this entire `strengthening_2026-09-14` directory beside the original `cpu_study_2026-09-13` and `provenance` directories. Use the original study's pinned environment with xlrd. The script reads those source directories and writes only this audit directory. Preserve a copy of this completed audit before rerunning because its output tables are regenerated. Verify hashes in the summary before interpreting a reproduction. A fresh script execution is not a fresh independent biological experiment.
