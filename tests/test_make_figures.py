"""Summary-provenance and missing-data guards; synthetic values stay in tests."""
import pytest
from src.make_figures import ARM_ORDER, prepare_plot_data, interval


@pytest.fixture
def summaries():
    rows = [{"arm": a, "failed_point_pairs": 0, "failed_range_pairs": 0,
             "macro_point_operational": .5, "eligible_pairs": 10,
             "resolved_interval_pairs": 2} for a in ARM_ORDER]
    shared = {"scenario": "primary", "partition": "development", "analysis_plan_sha256": "plan"}
    evaluation = {**shared, "prediction_sha256": "predictions", "primary_fixed_panel": {
        "status": "complete_predictions", "contributing_target_panel": ["m1", "m2"],
        "contexts": 3, "origins": 3, "arms": rows,
        "primary": {"absolute_macro_C": {a: .5 for a in ARM_ORDER if a != "uniform_pfba"},
                    "delta_grid": 0., "delta_identity": 0.},
        "bootstrap": {"interval_status": "conditional_on_estimability", "invalid_draws": 0,
                      "resamples": 2000, "absolute_C_intervals": {a: [.4, .6] for a in ARM_ORDER},
                      "delta_grid_interval": [-.1, .1], "delta_identity_interval": [-.1, .1]}}}
    sensitivity = {**shared, "predictions_sha256": "predictions", "primary_fixed_panel": {
        "targets": 2, "contexts": 3, "all_distinct_origin_context_target_pairs": 6,
        "context_target_cells": 6, "changed_fraction": .5, "strict_order_reversal_fraction": 1/6,
        "tie_transition_without_strict_reversal_pairs": 2, "prediction_tolerance": 1e-5}}
    sensitivity["range_diagnostics_primary"] = {"point_tolerance": 1e-5, "interval_tolerance": 1e-5,
        "per_arm": {a: {"context_target_cells": 6, "total_targets": 2,
                        "within_context_width_above_tolerance_count": 0,
                        "constant_prediction_targets": 2, "all_zero_prediction_targets": 1} for a in ARM_ORDER}}
    return evaluation, sensitivity, {"point_tolerance": 1e-5, "interval_tolerance": 1e-5}


def test_saved_values_and_distinct_denominators(summaries):
    plotted = prepare_plot_data(*summaries)
    assert plotted["contrasts"]["delta_grid"]["estimate"] == 0
    assert plotted["coverage"][0]["fraction"] == .2
    assert plotted["tie_transition_fraction"] == 2/6
    assert plotted["context_target_cells"] == 6


def test_hash_mismatch_cannot_mix_summary_runs(summaries):
    e, s, p = summaries
    s["predictions_sha256"] = "another run"
    with pytest.raises(ValueError, match="hashes differ"):
        prepare_plot_data(e, s, p)


def test_unreliable_ci_never_becomes_point_only_primary(summaries):
    e, s, p = summaries
    e["primary_fixed_panel"]["bootstrap"]["interval_status"] = "unreliable_invalid_draw_rate"
    with pytest.raises(ValueError, match="not reportable"):
        prepare_plot_data(e, s, p)


def test_reference_failure_is_omitted_not_plotted_at_chance(summaries):
    e, s, p = summaries
    e["primary_fixed_panel"]["arms"][0]["failed_point_pairs"] = 10
    plotted = prepare_plot_data(e, s, p)
    assert "uniform_pfba" not in {r["arm"] for r in plotted["absolute"]}
    e["primary_fixed_panel"]["arms"][1]["failed_point_pairs"] = 1
    with pytest.raises(ValueError, match="Missing successful primary arm"):
        prepare_plot_data(e, s, p)


def test_transform_threshold_cannot_silently_drift_from_plan(summaries):
    e, s, p = summaries
    s["primary_fixed_panel"]["prediction_tolerance"] = 1e-6
    with pytest.raises(ValueError, match="threshold differs"):
        prepare_plot_data(e, s, p)


def test_interval_uses_saved_endpoints_without_forcing_point_inside():
    # Percentile intervals need not contain the sample estimate. Drawing uses
    # hlines plus separate points, not negative error-bar lengths or altered CIs.
    assert interval([.6, .7], "example") == (.6, .7)
    with pytest.raises(ValueError, match="Reversed interval"):
        interval([.7, .6], "example")
