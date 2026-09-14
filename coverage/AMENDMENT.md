# Post-analysis prediction agreement and abstention audit

Recorded **2026-09-13 18:39:51 UTC / 2026-09-14 01:39:51 Asia/Bangkok**, before calculating the new agreement statistics. This is a new, outcome-free descriptive audit of previously completed predictions. It is not preregistration of the original study, a replacement of its protocol lock, or a predictive-performance analysis. Earlier published-in-package sensitivity/range counts have already been inspected.

## Fixed scope

- Use the original 52 primary targets from `protocol/analysis_plan.json`, the recorded test partition and origin groups from `manifests/contexts.tsv`, and the saved `results/reserved_v3/predictions.tsv`.
- Include primary, half_serum and lower_task scenarios separately. Do not select a favorable scenario, target or encoding.
- Primary arm set: all five power encodings, gamma 0.50, 0.80, 1.00, 1.25 and 2.00. Sensitivity arm set: those five plus ordinal. Exclude the uniform reference from agreement because it supplies no context-specific prediction.
- Enumerate every unordered pair of distinct profiles from different origin groups for every fixed target. Do not read CORE observations, apply an observed-separability filter, or calculate agreement with outcomes.
- Orient points and range endpoints using the frozen target sign factor. Use the original numerical point/order threshold **epsilon = 0.00001 canonical model units**. Normalize a tolerance-level inverted saved interval to its midpoint only as in the recorded evaluation, with its original numerical tolerance of 0.000001. Reject material inversion or point-containment failure.

## Fixed decisions and outputs

For each scenario, target, profile pair and arm set:

1. Record each arm's point-order sign: +1 if the oriented first-minus-second point exceeds epsilon, -1 if below negative epsilon, otherwise 0.
2. Report unanimous non-tied ordering separately from unanimous ties. Partition the remainder into mixed ties with one nonzero direction, or opposing nonzero directions.
3. **Matched-encoding interval support:** call +1 only if every arm satisfies lower(first, arm) > upper(second, arm) + epsilon; call -1 for the reverse inequality. Otherwise abstain.
4. **Independently varying encoding support:** call +1 only if min_arm lower(first) > max_arm upper(second) + epsilon; call -1 for the reverse inequality. Otherwise abstain. This stricter definition allows each profile a different encoding; it is a robustness envelope, not the primary common-encoding experiment.
5. Verify that both interval-supported decisions have the same strict sign as every point decision, and that the independently varying rule is a subset of matched support. Report and stop on inconsistencies rather than hiding them.

Save a per-pair ledger, per-target counts and whole-panel counts/coverage, exact input and source SHA256 digests, and execution timing. Coverage denominators retain all fixed targets/pairs, including constant and zero predictions. No p-values, confidence intervals, risk calibration, accuracy claims or thresholds selected by results will be added. Profile pairs share origins and are not independent sample sizes; arm sets and scenarios are not independent experiments.

## Verification and interpretation

Before the real-data calculation, test all-zero predictions, opposing signs, mixed ties, positive and negative unanimous orderings, matched support without independent-encoding support, sign orientation, numerical boundaries, origin exclusion and invalid inputs. Inspect an actual supported ledger case against its saved endpoints after execution.

High agreement can arise from all predictions being constant, so it must not be called informative coverage. Narrow within-arm ranges and unanimous signs do not establish biological correctness, statistical confidence, calibration or a probability of error. CORE endpoint quality remains an independent unresolved issue. All originals and solver outputs remain unchanged, with no optimization rerun.
