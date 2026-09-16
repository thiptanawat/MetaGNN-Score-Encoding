# Final reserved figure QA

All three final PNGs are covered by individual full-resolution visual inspection. Figure 1 was inspected again after the only presentation repair. Figures 2 and 3 were inspected individually before that repair and verified byte-identical after regeneration.

- Figure 1: PASS. Panel A retains the full 0–1 concordance scale; panel B shows a labelled 0.47–0.53 detail view of the identical saved values. The detail makes the short reserved confidence intervals visible. Both panels preserve the same arm order and 0.5 reference. Labels, uncertainty note, 47-context/44-origin/52-target scope and 0/2000 invalid-draw count fit without overlap.
- Figure 2: PASS. All six interval/point pairs, zero reference, scenario labels, direction statement and legend are visible. Marginal conditional uncertainty and nonindependence of scenario results remain explicit.
- Figure 3: PASS. All four panels and captions fit. A/B distinguish prediction-change denominators from C's observed-eligible pairs. D's 2444-cell and 52-target counts match the primary summary and distinguish width, constancy and zero predictions. Numerical thresholds are 1e-5 throughout; no exact-uniqueness claim is made.

The initial full-scale-only Figure 1 obscured narrow intervals behind its markers; this was corrected by adding the labelled detail panel, without changing any saved number. The plotting source is outside the scientific protocol lock. Both development and reserved figure sets and manifests were regenerated after the repair. The revised development Figure 1 was also individually inspected at full resolution and passed; its Figures 2/3 remained byte-identical to their earlier inspected versions.

Six plotting validation tests passed. Every source, input and output SHA-256 in both final manifests was checked against the corresponding file. Every exported PDF is structurally valid and contains one page; the PNGs are the directly inspected renders of the same Matplotlib figures. This record does not claim a separate raster inspection of each PDF export. Scientific source/configuration files and primary-phase evidence were not changed.

## Inspected final reserved PNG identities

- `figure1_absolute_concordance.png`: `59a34afbc011e4edec49471a99cf688ba33620745e045a9cc04cf2ead0e8a28a`
- `figure2_paired_contrasts.png`: `e6b34fcd0a22f86a86c2d28b885ae530c165dacdbfaa861d056512583e5a1fd3`
- `figure3_descriptive_sensitivity.png`: `783e513c05e50d03af688de3055177eb668076cb2e2a75e50c0f85954d0b28ea`
