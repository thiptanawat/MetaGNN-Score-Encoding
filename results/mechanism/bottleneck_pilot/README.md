# Nested-range pilot: which optimization stage narrows exchanges?

**Completed 14 September 2026 Bangkok.** All requested work finished in **475.39 seconds (7 minutes 55 seconds)**, within the eight-minute budget: one shared B0 job, all 22 B1 jobs, and comparison with 22 matching saved B2 records. No failed or unstarted jobs remain. At most four CPU workers were used. No CORE observation table or performance result was read, and no original scientific file or lock was changed.

This is a **post-analysis, outcome-free computational follow-up**, specified in the [dated amendment](AMENDMENT.md) before the new LPs. The original paper reported ranges under two objective caps; it did not establish which optimization stage narrowed those ranges. This pilot addresses that specific mechanistic question within the model. It does not establish a physiological bottleneck or a new predictive-performance result.

## Comparison and results

- **B0:** mass balance, original primary medium/bounds and fixed primary biomass task, with neither cost cap. Shared across profiles, computed once for all 52 fixed targets.
- **B1:** B0 plus the exact saved primary transcript-weighted cost cap. Computed for identity and ordinal encodings in all 11 development profiles.
- **B2:** saved ranges under both primary and secondary absolute-flux caps, before coordinate fixing. Reused only after comparability checks; not reoptimized here.

A range is *numerically fixed* when its width is at most **10⁻⁵ canonical model units**. This is the original numerical convention, not exact uniqueness. The absolute task was 3.2852230666462234 model units, matching the original primary configuration.

| Encoding | Profile–target combinations | Fixed at B0* | Fixed at B1 | Fixed at B2 | First fixed at B1 | First fixed at B2 | Still variable at B2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Identity magnitude | 572 | 33 | 523 | 569 | 490 | 46 | 3 |
| Ordinal | 572 | 33 | 514 | 565 | 481 | 51 | 7 |

*B0 fixes **3 of 52 unique targets**, not 33 independently estimated targets. The table repeats that common result over 11 profiles per arm. The fixed exchanges are `EX_cmp_e`, `EX_hom__L_e` and `EX_sbt__D_e`, each with computed B0 interval [0, 0].*

Across individual profiles, B1 fixes 45–50 identity targets and 44–50 ordinal targets; B2 fixes 49–52 and 51–52, respectively. Adding E1 shrinks all 539 initially non-fixed profile–target widths per arm by more than 10⁻⁵. Adding E2 shrinks 44 identity and 52 ordinal widths by more than that threshold. Width contraction and crossing the fixedness threshold are different events, so those counts need not equal the “first fixed at B2” counts.

**Interpretation:** the primary cost cap accounts for most first numerical fixation in this development pilot; the secondary cap makes a smaller additional contribution. This concerns **within-profile ranges**. Narrow intervals in different profiles can have different centers. The experiment therefore does not, by itself, explain the previously reported **between-profile constancy**, establish preservation/loss of biological information, or identify a causal reaction bottleneck. The development results cannot automatically be generalized to every reserved profile or scenario. The underlying ranges and constraints are numerical model objects, not statistical confidence intervals or measured physiological bounds.

## Comparability and verification

The run manifest pins raw Recon3D, config, environment, score inputs, context/target metadata, inspected production sources and every reused B2 JSON/NPZ. Saved reaction order and reconstructed weights matched exactly. Fresh primary optima also matched saved optima exactly, and the saved cap formulas were verified. Reused full B2 vectors passed mass balance, bounds, task, both cost caps and reported-point checks.

- **2,414 recorded new LP attempts** were optimal on their first attempt: 104 B0 extrema plus 22 × (one primary reoptimization + 104 B1 extrema). Model construction also performed 23 checked task-maximum solves; these are setup solves outside that recorded-attempt subtotal.
- Every new extremum's full vector was audited and saved. [Readback](readback_check.json) checked **2,392 full extremum vectors**, with zero discrepancy between their exchange coordinates and recorded extrema.
- Maximum range-nesting deviations: **1.58 × 10⁻¹²** for B0/B1 and **7.32 × 10⁻¹⁰** for B1/B2; maximum saved-point containment deviation **3.68 × 10⁻¹¹**. All are below the fixed 10⁻⁶ validation tolerance.
- Worst new-vector mass-balance residual: **2.08 × 10⁻¹¹**; bound violation: **7.78 × 10⁻¹⁰**; primary-cap violation: **2.39 × 10⁻⁹**. No validation failures or nonmonotone fixedness classifications occurred.
- Three [toy logic tests](toy_test_record.txt) passed. All original input/source hashes remained unchanged.

Saved B2 extrema were not independently reoptimized, and numerical checks do not constitute exact mathematical certificates. Ranges retain the original reaction-coordinate convention; widths and nesting are unaffected by the exchange sign orientation. No new biological accuracy or uncertainty estimate is supplied.

## Files and rerunning

Start with [summary.json](summary.json), [run manifest](run_manifest.json), [per-profile/target range ledger](range_stage_ledger.tsv), and [readback verification](readback_check.json). Per-job JSONs retain ranges, cap comparisons, vector audits and solve ledgers; matching `*_extrema.npz` files contain complete new extremum vectors. [Execution log](execution.log) and [file hashes](SHA256SUMS.txt) retain completion evidence.

The [pilot script](nested_range_pilot.py) uses the original study's matching environment and inspected model helpers. The [readback script](readback_check.py) performs no solver work. To repeat, copy this pilot directory's scripts and amendment into a **new** output directory, then run with the existing study root and matching Python environment. The runner refuses to overwrite an existing run manifest. Keep bytecode writing disabled with `PYTHONDONTWRITEBYTECODE=1`. No further solver work is needed for the present bounded result.
