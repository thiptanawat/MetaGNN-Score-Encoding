# Post-analysis nested-range bottleneck pilot

Recorded **2026-09-13 18:49:26 UTC / 2026-09-14 01:49:26 Bangkok**, before this pilot's new LPs. Original outcomes and all completed prediction/sensitivity summaries have previously been inspected. This is an outcome-free computational follow-up, not preregistration, biological validation or an alteration of the frozen study.

## Question and fixed scope

Which stage of the declared optimization narrows each exchange's within-profile feasible range: base constraints, primary transcript-weighted cost, or secondary absolute-flux parsimony? The existing comparison only reports ranges under both caps and cannot attribute their narrowness to one stage.

Use the original primary medium and absolute biomass task obtained by the original primary build; all **52 fixed primary targets**; all **11 development profiles**; **identity magnitude and ordinal** encodings. No reserved-profile expansion or other scenarios. No CORE observations, errors, accuracy or performance-based selection.

- **B0:** original mass balance, reaction bounds, primary medium and fixed biomass, without either cost cap. Compute once because these constraints are shared across profiles.
- **B1:** B0 plus the exact matching saved primary-cost cap. Recompute and compare the primary optimum first; reconstruct and match all reaction costs to the saved NPZ.
- **B2:** reuse saved ranges under both original caps, only after raw model, config, source fingerprints, reaction order, weights, cap formula, metadata, statuses and saved full-vector feasibility are verified. B2 is not reoptimized in this pilot.

Keep GLPK feasibility/optimality tolerances 10⁻⁹, validation/nesting tolerance 10⁻⁶, and range-resolution threshold 10⁻⁵ canonical model units. Width <=10⁻⁵ is numerically fixed for this report; it is not proof of exact uniqueness. Tolerance-level reversed saved intervals follow the original midpoint convention. No cap relaxation, coefficient changes or success-by-dropping-failures.

For every new extremum, audit the full returned vector for mass balance, bounds and the applicable cost cap. Check each reused B2 point against B0/B1 bounds, costs and ranges; check nested ranges B2 ⊆ B1 ⊆ B0 within the fixed tolerance. Preserve all solver attempts and failures; allow the same one advanced-basis retry with unchanged LP. Store full extremum vectors, ranges, audits and exact input/code hashes.

## Reporting and resource bound

Report numerical fixedness and width contraction separately at each stage, including transitions already-fixed-in-B0, first-fixed-at-B1, first-fixed-at-B2 and still-variable. Keep per-profile/target counts and incomplete coverage visible. Within-profile fixedness is not between-profile constancy. Width shrinkage attributes a computational restriction to the added constraint set, not a physiological bottleneck or a uniquely causal reaction.

Use at most four CPU workers and an approximately eight-minute total LP runtime budget. A per-LP time limit is at most the original 60 seconds and may be shortened only to respect the remaining global time budget; timed-out or unstarted work remains partial. No automatic expansion. Original source, results, protocol lock and raw inputs remain untouched. Python bytecode writing is disabled for reused source imports.
