"""Separate toy validation for the post-lock descriptive range audit."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

MODULE = Path(__file__).with_name("cross_encoding_ranges.py")
spec = importlib.util.spec_from_file_location("range_audit", MODULE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class RangeAuditToys(unittest.TestCase):
    def test_disjoint_intervals_have_positive_certificate(self):
        w = audit.disjoint_witness({"a": (0, -1, 1), "b": (3, 2, 4)}, 1e-5)
        self.assertEqual((w["higher_arm"], w["lower_arm"], w["gap"]), ("b", "a", 1))

    def test_point_changes_do_not_certify_overlapping_ranges(self):
        self.assertIsNone(audit.disjoint_witness({"a": (0, 0, 3), "b": (3, 0, 3)}, 1e-5))

    def test_tolerance_boundary_is_not_strictly_separated(self):
        self.assertIsNone(audit.disjoint_witness({"a": (0, 0, 0), "b": (1e-5, 1e-5, 1e-5)}, 1e-5))

    def test_robust_reversal_needs_two_opposite_orderings(self):
        first = {"a": (3, 2, 4), "b": (0, -1, 1)}
        second = {"a": (0, -1, 1), "b": (3, 2, 4)}
        w = audit.reversal_witness(first, second, 1e-5)
        self.assertEqual(w["first_above_second_arm"], "a")
        self.assertEqual(w["second_above_first_arm"], "b")
        self.assertEqual(w["minimum_gap"], 1)

    def test_overlapping_intervals_do_not_certify_point_reversal(self):
        first = {"a": (0, 0, 3), "b": (3, 0, 3)}
        second = {"a": (3, 0, 3), "b": (0, 0, 3)}
        self.assertIsNone(audit.reversal_witness(first, second, 1e-5))

    def test_derivatives_are_excluded_only_within_origin(self):
        contexts = ["a", "b", "c", "d"]
        pairs = audit.distinct_origin_pairs(contexts, {"a": "x", "b": "x", "c": "y", "d": "z"})
        self.assertEqual(len(pairs), 5)
        self.assertNotIn(("a", "b"), pairs)
        self.assertIn(("a", "c"), pairs)
        self.assertIn(("b", "c"), pairs)

    def test_small_inversion_midpoint_and_orientation(self):
        values, inverted = audit.normalize_record(5e-8, 1e-7, 0, 1e-6, -1.0)
        self.assertTrue(inverted)
        self.assertEqual(values, (-5e-8, -5e-8, -5e-8))

    def test_material_and_nonfinite_records_fail(self):
        for record in [(0, 1e-3, 0), (float("nan"), 0, 1), (0, 0, float("inf")), (2, 0, 1)]:
            with self.subTest(record=record), self.assertRaises(ValueError):
                audit.normalize_record(*record, 1e-6)

    def test_end_to_end_without_an_outcome_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "protocol").mkdir()
            (root / "manifests").mkdir()
            contexts = root / "manifests/contexts.tsv"
            contexts.write_text("context_id\torigin_group\tpartition\na\tx\ttest\nb\ty\ttest\n")
            arms = [f"magnitude_g{k}" for k in range(5)]
            targets = [{"metabolite_id": "m", "exchange_id": "e", "sign_factor": 1.0}]
            plan = {"point_tolerance": 1e-5, "interval_tolerance": 1e-5,
                    "numerical_tolerance": 1e-6, "comparison_arms": arms + ["ordinal"],
                    "chemistry_targets": targets, "primary_targets": targets,
                    "input_sha256": {"manifests/contexts.tsv": audit.sha256(contexts)}}
            plan_path = root / "protocol/analysis_plan.json"
            plan_path.write_text(json.dumps(plan))
            (root / "protocol/PROTOCOL_LOCK.json").write_text(json.dumps({
                "status": "locked", "reserved_evaluation_authorized": True,
                "analysis_plan_sha256": audit.sha256(plan_path)}))
            predictions = root / "predictions.tsv"
            with predictions.open("w") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
                    "context_id", "arm", "exchange_id", "point", "range_lo", "range_hi", "scenario", "status"])
                writer.writeheader()
                for i, arm in enumerate(arms):
                    for context, point in [("a", i), ("b", 4-i)]:
                        writer.writerow(dict(context_id=context, arm=arm, exchange_id="e",
                                             point=point, range_lo=point-.1, range_hi=point+.1,
                                             scenario="primary", status="optimal"))
            result = audit.analyze(root, predictions, "primary", "test")
            p = result["primary_fixed_panel"]
            self.assertEqual(p["disjoint_interval_cells"], 2)
            self.assertEqual(p["robust_interval_reversal_pairs"], 1)
            self.assertEqual(p["distinct_origin_context_target_pairs"], 1)
            self.assertFalse(result["core_outcome_table_read"])
            self.assertEqual(result["new_lp_solves"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
