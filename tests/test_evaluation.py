import itertools
import json
from collections import Counter

import numpy as np
import pytest

from src.evaluation import (
    COMPARISON_ARMS, Pair, aggregate_fixed, build_pairs, check_test_lock,
    evaluate_panel, interval_order, mean_observations, paired_group_bootstrap,
    point_score, prediction_index, prediction_values, sha256, write_tsv,
)
from src.prepare_analysis import estimability_risks, prepare, support_distribution


def pair(met="m", a="a", b="b", truth=1):
    return Pair(met, f"EX_{met}", a, b, a, b, truth)


def test_production_comparators_known_order_and_negative_rates():
    for a, b, truth, answer in [(6, 2, 1, 1), (6, 2, -1, 0), (-1, -9, 1, 1),
                                (-9, -1, 1, 0), (0, 0, -1, .5)]:
        assert point_score(a, b, truth, 1e-6) == answer
        assert point_score(b, a, -truth, 1e-6) == answer
    assert interval_order((5, 6), (1, 2), 1e-6) == 1
    assert interval_order((-9, -8), (-1, -.5), 1e-6) == -1
    assert interval_order((1, 5), (5, 9), 1e-6) == 0
    assert point_score(1, 1 + 1e-7, -1, 1e-6) == .5


def test_macro_is_not_pooled_and_single_pair_is_not_dropped():
    pairs = [pair("m1")] + [pair("m2", str(i), str(i + 1)) for i in range(9)]
    scores = np.array([1] + [0] * 9)
    result = aggregate_fixed(scores, pairs, ["m1", "m2"])
    assert result["valid"]
    assert result["macro"][0] == .5
    assert np.mean(scores) == .1


def test_fixed_panel_zero_denominator_invalid_not_chance_or_drop():
    result = aggregate_fixed([1, 0], [pair("m1"), pair("m2", "b", "c")],
                             ["m1", "m2"], weights=[3, 0])
    assert not result["valid"]
    assert result["missing_targets"] == ["m2"]
    assert np.isnan(result["macro"][0])
    # Multiplicity changes weights; it never changes the panel or minimum count.
    valid = aggregate_fixed([1, 0], [pair("m1"), pair("m2", "b", "c")],
                            ["m1", "m2"], weights=[9, 1])
    assert valid["macro"][0] == .5


def test_bootstrap_exact_support_distribution_against_exhaustive_draws():
    supports, mass = support_distribution(3)
    empirical = Counter()
    for draw in itertools.product(range(3), repeat=3):
        empirical[sum(1 << i for i in set(draw))] += 1
    for mask, probability in zip(supports, mass):
        assert probability == pytest.approx(empirical[int(mask)] / 27)
    risks, bad, probability = estimability_risks([pair()], ["m"], ["a", "b", "c"])
    expected = sum(count for mask, count in empirical.items() if not (mask & 1 and mask & 2)) / 27
    assert risks["m"]["invalid_probability"] == pytest.approx(expected)
    assert probability[bad["m"]].sum() == pytest.approx(expected)


def test_bootstrap_rejects_missing_target_draws_and_is_seeded():
    pairs = [pair("m1"), pair("m2", "b", "c")]
    scores = np.full((len(COMPARISON_ARMS), len(pairs)), .5)
    scores[-1] = [1, 0]
    one = paired_group_bootstrap(scores, pairs, ["m1", "m2"], COMPARISON_ARMS,
                                 ["a", "b", "c"], resamples=300, seed=11)
    two = paired_group_bootstrap(scores, pairs, ["m1", "m2"], COMPARISON_ARMS,
                                 ["a", "b", "c"], resamples=300, seed=11)
    assert one == two
    assert one["invalid_draws"] > 0
    assert one["valid_draws"] + one["invalid_draws"] == 300
    assert one["interval_status"] == "unreliable_invalid_draw_rate"
    assert one["delta_grid_interval"] is None
    assert one["fixed_panel_size"] == 2


def test_within_origin_exclusion_is_explicit():
    contexts = {"a": {"origin_group": "same"}, "b": {"origin_group": "same"},
                "c": {"origin_group": "other"}}
    pairs = build_pairs({("m", "a"): 1, ("m", "b"): 2, ("m", "c"): 8},
                        contexts, [{"metabolite_id": "m", "exchange_id": "EX_m"}], {"m": .1})
    assert len(pairs) == 2
    assert all(p.origin_i != p.origin_j for p in pairs)
    with pytest.raises(ValueError, match="within-origin"):
        same = Pair("m", "EX_m", "a", "b", "same", "same", 1)
        paired_group_bootstrap(np.ones((6, 1)), [same], ["m"], COMPARISON_ARMS, ["same", "other"])


def test_prediction_sign_conversion_and_failed_ranges():
    target = {"metabolite_id": "m", "exchange_id": "EX_m", "sign_factor": -1}
    idx = {("a", "ordinal", "EX_m"): {"status": "optimal", "point": 3, "range_lo": 2, "range_hi": 4}}
    point, bounds, note = prediction_values(idx, "a", "ordinal", target, 1e-6)
    assert point == -3 and bounds == (-4, -2) and note is None
    idx[("a", "ordinal", "EX_m")]["range_lo"] = 5
    point, bounds, note = prediction_values(idx, "a", "ordinal", target, 1e-6)
    assert point == -3 and bounds is None and note == "material_range_inversion"
    idx[("a", "ordinal", "EX_m")]["point"] = "nan"
    assert prediction_values(idx, "a", "ordinal", target, 1e-6)[0] is None


def test_prediction_duplicates_fail_closed():
    row = {"scenario": "main", "context_id": "a", "arm": "ordinal", "exchange_id": "EX_m"}
    with pytest.raises(ValueError, match="Duplicate prediction"):
        prediction_index([row, row], "main")


def test_missing_predictions_do_not_become_a_primary_chance_result():
    contexts = {c: {"origin_group": c} for c in "ab"}
    targets = [{"metabolite_id": "m", "exchange_id": "EX_m", "sign_factor": 1}]
    plan = {"point_tolerance": 1e-6, "interval_tolerance": 1e-6,
            "numerical_tolerance": 1e-6, "maximum_invalid_fraction": .01}
    result = evaluate_panel({}, {("m", "a"): 2, ("m", "b"): 1}, contexts, targets,
                            {"m": .1}, plan, resamples=20, seed=1)
    assert result["primary"] is None
    assert result["bootstrap"] is None
    assert result["status"] == "primary_not_estimable_or_prediction_failure"
    assert result["arms"][0]["macro_point_operational"] == .5
    assert result["arms"][0]["failed_point_pairs"] == 1


def test_constant_predictions_keep_chance_and_known_contrast():
    contexts = {c: {"origin_group": c} for c in "abc"}
    target = {"metabolite_id": "m", "exchange_id": "EX_m", "sign_factor": 1}
    idx = {}
    for arm in COMPARISON_ARMS:
        for value, context in enumerate("abc"):
            pt = value if arm == "ordinal" else 0
            idx[(context, arm, "EX_m")] = {"status": "optimal", "point": pt, "range_lo": pt, "range_hi": pt}
    plan = {"point_tolerance": 1e-6, "interval_tolerance": 1e-6,
            "numerical_tolerance": 1e-6, "maximum_invalid_fraction": .01}
    result = evaluate_panel(idx, {("m", c): i for i, c in enumerate("abc")}, contexts,
                            [target], {"m": .1}, plan, resamples=0, seed=1)
    assert result["primary"]["delta_grid"] == .5
    assert result["primary"]["delta_identity"] == .5
    assert result["primary"]["absolute_macro_C"]["magnitude_g1.00"] == .5


def test_development_filter_precedes_outcome_numeric_parsing():
    rows = [{"context_id": "dev", "metabolite_id": "m", "replicate": "1", "value": "3"},
            {"context_id": "test", "metabolite_id": "m", "replicate": "1", "value": "DO_NOT_PARSE_RESERVED"}]
    assert mean_observations(rows, {"dev"}) == {("m", "dev"): 3}
    with pytest.raises(ValueError, match="Duplicate outcome replicate"):
        mean_observations(rows + [rows[0]], {"dev"})


def test_reserved_evaluation_delegates_to_shared_strict_gate(tmp_path, monkeypatch):
    import sys
    import types
    plan = tmp_path / "protocol/analysis_plan.json"
    plan.parent.mkdir()
    plan.write_text("{}")
    check_test_lock("development", plan, None)
    calls = []
    def strict_gate(root, lock):
        calls.append((root, lock))
        raise ValueError("No authorized locked protocol")
    monkeypatch.setitem(sys.modules, "src.verify_inputs", types.SimpleNamespace(verify_protocol_lock=strict_gate))
    with pytest.raises(ValueError, match="No authorized"):
        check_test_lock("test", plan, None, tmp_path)
    assert calls == [(tmp_path, None)]
    with pytest.raises(ValueError, match="canonical"):
        check_test_lock("test", tmp_path / "other.json", None, tmp_path)


def test_prepare_uses_only_development_observations_and_exact_risk(tmp_path):
    contexts = [{"context_id": f"d{i}", "origin_group": f"d{i}", "partition": "development"} for i in range(5)]
    contexts.append({"context_id": "heldout", "origin_group": "heldout", "partition": "test"})
    mapping = [{"metabolite_id": "good", "exchange_id": "EX_good", "status": "mapped", "primary_eligible": "true", "sign_factor": 1},
               {"metabolite_id": "constant", "exchange_id": "EX_constant", "status": "mapped", "primary_eligible": "true", "sign_factor": 1}]
    outcomes = [{"context_id": f"d{i}", "metabolite_id": m, "replicate": "1", "value": i if m == "good" else 0}
                for i in range(5) for m in ("good", "constant")]
    outcomes.append({"context_id": "heldout", "metabolite_id": "good", "replicate": "1", "value": "POISON_RESERVED"})
    noise = [{"metabolite_id": m, "tolerance": .1} for m in ("good", "constant")]
    for relative, rows in [("manifests/contexts.tsv", contexts), ("manifests/metabolite_map.tsv", mapping),
                           ("data/outcomes_long.tsv", outcomes), ("data/assay_noise_development.tsv", noise)]:
        write_tsv(tmp_path / relative, rows, list(rows[0]))
    plan, ledger, pairs = prepare(tmp_path)
    assert [t["metabolite_id"] for t in plan["primary_targets"]] == ["good"]
    assert len(plan["chemistry_targets"]) == 2
    assert plan["development_exact_joint_invalid_probability"] == pytest.approx(5 / 5 ** 5)
    assert plan["development_union_bound_invalid_probability"] <= .01
    assert all("heldout" not in (p.context_i, p.context_j) for p in pairs)


def test_absolute_intervals_share_the_paired_draws_for_constant_known_scores():
    pairs=[pair('m','a','b'),pair('m','a','c'),pair('m','b','c')]
    scores=np.full((len(COMPARISON_ARMS),3),.5);scores[-1]=1.
    result=paired_group_bootstrap(scores,pairs,['m'],COMPARISON_ARMS,['a','b','c'],resamples=300,seed=11,maximum_invalid_fraction=.5)
    assert result['absolute_C_intervals']['ordinal']==[1.,1.]
    assert result['absolute_C_intervals']['magnitude_g1.00']==[.5,.5]
    assert result['delta_grid_interval']==[.5,.5]
    assert len(result['diagnostic_valid_draw_quantiles'])==2
