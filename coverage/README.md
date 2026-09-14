# Prediction agreement and abstention coverage audit

Completed 14 September 2026 Bangkok. This is a **new post-analysis, outcome-free descriptive calculation** from saved predictions. No CORE outcome table was read, no accuracy or error risk was calculated, and no model was reoptimized. The original study and its lock remain unchanged.

The dated [amendment](AMENDMENT.md) was written before this calculation at 18:39:51 UTC; its SHA256 is `7f97bdc7c1b1e889d8150fc12ea5f62014730e8d5b8ef20304926742f5fcbfda`. Execution began at 18:42:59 UTC and took 3.99 seconds. Prior sensitivity results had already been inspected; this is not a new confirmatory analysis.

## Main finding

Most agreement is agreement on **ties**, rather than a common informative ordering. With the five power encodings, the primary scenario has 49,442 unanimous ties among 56,056 target–profile-pair comparisons (88.20%), but only 3,214 unanimously non-tied orderings (5.73%). Including ordinal reduces the latter to 1,443 (2.57%). These quantities describe computational coverage, not biological correctness.

Every row below has the same denominator: **56,056 comparisons = 52 fixed targets × 1,078 different-origin profile pairs**. These are dependent comparisons from 47 profiles and 44 origins, not 56,056 independent observations.

| Scenario / arm set | Unanimous non-tied points | Unanimous ties | Matched-encoding interval support | Independently varying encoding support |
|---|---:|---:|---:|---:|
| Primary / five powers | 3,214 (5.734%) | 49,442 (88.201%) | 3,204 (5.716%) | 2,636 (4.702%) |
| Half serum / five powers | 1,717 (3.063%) | 47,019 (83.879%) | 1,717 (3.063%) | 845 (1.507%) |
| Lower task / five powers | 2,377 (4.240%) | 46,256 (82.517%) | 2,364 (4.217%) | 1,608 (2.869%) |
| Primary / powers + ordinal | 1,443 (2.574%) | 46,829 (83.540%) | 1,438 (2.565%) | 301 (0.537%) |
| Half serum / powers + ordinal | 542 (0.967%) | 45,908 (81.897%) | 542 (0.967%) | 4 (0.007%) |
| Lower task / powers + ordinal | 721 (1.286%) | 45,253 (80.728%) | 721 (1.286%) | 44 (0.078%) |

The point categories partition the denominator: unanimous positive, unanimous negative, unanimous ties, mixed ties with one strict direction, and opposing strict directions. Interval support overlaps the unanimous non-tied categories; it must not be added to them.

**Matched support** requires the same non-tied ordering throughout each arm's saved ranges, comparing the same arm in the two profiles. **Independently varying support** additionally permits a different arm in each profile, comparing the union-range envelopes. It is deliberately stronger than the original common-encoding experiment. The threshold is the frozen **10⁻⁵ canonical model units**. Neither type of range is a statistical confidence interval.

## Verification

All 11 [toy tests](toy_test_record.txt) passed. They include all-zero abstention, opposite signs, mixed ties, unanimous points with overlapping ranges, strict-threshold boundaries, signed orientation and invalid inputs. All real interval-supported calls agree with every arm's strict point sign; all independently varying calls are also matched-supported. There were zero inconsistencies. Twenty-two tolerance-level inverted input intervals were normalized to their midpoint using the original evaluation convention; material inversions would stop execution.

The compressed ledger was read back independently: **336,336 rows**, six scenario/arm-set groups, with category counts matching the JSON summaries. All recorded input hashes still match. See [ledger readback](ledger_readback.json).

Manual endpoint checks:

- Primary alanine, BT-549 versus MCF7: minimum MCF7 lower endpoint across powers is 11.92443253299; maximum BT-549 upper endpoint is 11.73025091846. The negative-direction gap is 0.19418161453, exceeding 10⁻⁵, so both interval rules support the call.
- Primary alanine, BT-549 versus NCI/ADR-RES: the minimum matched negative margin is 0.00003034744, but the independently varying margin is −0.10424641537. Matched support therefore passes and the stricter rule abstains.
- Primary 4-aminobutyrate, BT-549 versus MCF7: every point is zero; all five signs are zero, and both support rules abstain.

## Files and reproduction

- [JSON summary, per-target counts, examples and exact input hashes](agreement_summary.json)
- [Summary TSV](agreement_summary.tsv) and [per-target TSV](agreement_per_target.tsv)
- [Full compressed per-pair ledger](agreement_pair_ledger.tsv.gz)
- [Independent audit script](agreement_audit.py), [toy tests](test_agreement_audit.py), [execution record](execution_record.json), [file hashes](SHA256SUMS.txt)

Uses Python's standard library. From this directory:

```bash
python3 -m unittest -v test_agreement_audit.py
python3 agreement_audit.py --root ../../cpu_study_2026-09-13 --output ./repeat_audit
```

The script refuses to overwrite an existing ledger. Originals are read-only inputs. The panel was selected using development observations in the original study; “outcome-free” describes this new calculation, not the entire study design. Low coverage may be useful for reporting when predictions are uninformative, but no calibrated abstention policy, clinical safety, superior predictor or biological validation follows from these results. The endpoint-quality issue remains unresolved.
