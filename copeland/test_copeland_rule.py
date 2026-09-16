"""Known-answer checks for the Copeland reporting rule, written before any prediction was read."""
import unittest
from copeland_rule import (EPS, measured_direction, interval_direction, shared_encoding_direction,
                           point_category, agreement)


class Measured(unittest.TestCase):
    def test_complete_separation_upwards(self):
        self.assertEqual(measured_direction([5, 6, 7, 8], [1, 2, 3, 4]), 1)

    def test_complete_separation_downwards(self):
        self.assertEqual(measured_direction([1, 2, 3, 4], [5, 6, 7, 8]), -1)

    def test_one_overlapping_replicate_blocks_a_direction(self):
        self.assertEqual(measured_direction([5, 6, 7, 4], [1, 2, 3, 4]), 0)

    def test_uptake_values_are_ordered_as_numbers_not_magnitudes(self):
        # more negative means MORE uptake but a SMALLER number; the rule must order numerically
        self.assertEqual(measured_direction([-100, -110], [-300, -320]), 1)

    def test_lactate_21_percent_is_separable_upwards(self):
        dmso = [708.6, 1114.5, 1127.4, 950.1]
        bay = [1193.0, 1393.4, 1539.5, 1255.2]
        self.assertEqual(measured_direction(bay, dmso), 1)

    def test_glucose_21_percent_is_not_separable(self):
        dmso = [-411.5, -485.4, -572.2, -418.0]
        bay = [-377.9, -513.5, -800.3, -688.7]
        self.assertEqual(measured_direction(bay, dmso), 0)


class Intervals(unittest.TestCase):
    def test_disjoint_intervals_report_a_direction(self):
        self.assertEqual(interval_direction([(10.0, 11.0)], [(1.0, 2.0)]), 1)
        self.assertEqual(interval_direction([(1.0, 2.0)], [(10.0, 11.0)]), -1)

    def test_touching_intervals_abstain(self):
        self.assertEqual(interval_direction([(2.0, 3.0)], [(1.0, 2.0)]), 0)

    def test_a_single_overlapping_encoding_forces_abstention(self):
        a = [(10.0, 11.0), (1.5, 2.5)]          # the second encoding overlaps B
        self.assertEqual(interval_direction(a, [(1.0, 2.0)]), 0)

    def test_negative_intervals_are_ordered_numerically(self):
        self.assertEqual(interval_direction([(-2.0, -1.0)], [(-10.0, -9.0)]), 1)

    def test_a_gap_smaller_than_epsilon_abstains(self):
        self.assertEqual(interval_direction([(1.0 + EPS / 2, 2.0)], [(0.0, 1.0)]), 0)

    def test_a_gap_larger_than_epsilon_reports(self):
        self.assertEqual(interval_direction([(1.0 + EPS * 10, 2.0)], [(0.0, 1.0)]), 1)

    def test_empty_input_abstains(self):
        self.assertEqual(interval_direction([], [(1.0, 2.0)]), 0)


class SharedEncodingRule(unittest.TestCase):
    """The planned rule compares intervals within an encoding and then asks every encoding to agree."""

    def test_every_encoding_supports_the_same_direction(self):
        a = {"g1": [(10.0, 11.0), (10.5, 12.0)], "g2": [(5.0, 6.0)]}
        b = {"g1": [(1.0, 2.0)], "g2": [(1.0, 2.0), (0.5, 1.5)]}
        self.assertEqual(shared_encoding_direction(a, b), 1)
        self.assertEqual(shared_encoding_direction(b, a), -1)

    def test_one_encoding_overlapping_forces_abstention(self):
        a = {"g1": [(10.0, 11.0)], "g2": [(1.5, 2.5)]}
        b = {"g1": [(1.0, 2.0)], "g2": [(1.0, 2.0)]}
        self.assertEqual(shared_encoding_direction(a, b), 0)

    def test_encodings_disagreeing_in_direction_abstain(self):
        a = {"g1": [(10.0, 11.0)], "g2": [(0.0, 0.5)]}
        b = {"g1": [(1.0, 2.0)], "g2": [(1.0, 2.0)]}
        self.assertEqual(shared_encoding_direction(a, b), 0)

    def test_encodings_are_not_compared_with_each_other(self):
        # each encoding clears within itself, but encoding g2 of A sits below encoding g1 of B;
        # the pooled rule abstains, the shared rule reports
        a = {"g1": [(10.0, 11.0)], "g2": [(3.0, 3.5)]}
        b = {"g1": [(4.0, 5.0)], "g2": [(1.0, 2.0)]}
        self.assertEqual(shared_encoding_direction(a, b), 1)
        pooled_a = a["g1"] + a["g2"]; pooled_b = b["g1"] + b["g2"]
        self.assertEqual(interval_direction(pooled_a, pooled_b), 0)

    def test_missing_encoding_on_one_side_abstains(self):
        self.assertEqual(shared_encoding_direction({"g1": [(10.0, 11.0)]}, {"g1": [(1.0, 2.0)], "g2": [(0.0, 1.0)]}), 0)


class Points(unittest.TestCase):
    def test_unanimous_positive(self):
        self.assertEqual(point_category([5.0, 6.0], [1.0])[0], "unanimous_positive")
        self.assertEqual(point_category([5.0, 6.0], [1.0])[1], 1)

    def test_unanimous_negative(self):
        self.assertEqual(point_category([1.0], [5.0, 6.0])[1], -1)

    def test_identical_predictions_are_a_tie_not_a_disagreement(self):
        self.assertEqual(point_category([3.0, 3.0], [3.0])[0], "unanimous_tie")
        self.assertEqual(point_category([3.0, 3.0], [3.0])[1], 0)

    def test_opposing_signs_are_conflicting(self):
        self.assertEqual(point_category([5.0, 0.0], [1.0])[0], "conflicting")

    def test_a_direction_mixed_with_ties_is_neither_unanimous_nor_conflicting(self):
        self.assertEqual(point_category([5.0, 1.0], [1.0])[0], "mixed_with_ties")

    def test_a_tie_reports_no_direction(self):
        self.assertEqual(point_category([3.0], [3.0])[1], 0)


class Agreement(unittest.TestCase):
    def test_matching_directions_agree(self):
        self.assertIs(agreement(1, 1), True)
        self.assertIs(agreement(-1, -1), True)

    def test_opposite_directions_disagree(self):
        self.assertIs(agreement(1, -1), False)

    def test_abstention_is_not_a_disagreement(self):
        self.assertIsNone(agreement(0, 1))

    def test_an_unresolved_measurement_yields_no_judgement(self):
        self.assertIsNone(agreement(1, 0))


class SignConvention(unittest.TestCase):
    """Both sides use positive for secretion, so a correct call must not need a flip."""

    def test_secretion_increase_is_a_positive_direction_on_both_sides(self):
        measured = measured_direction([1200.0, 1300.0], [700.0, 900.0])     # more lactate out
        predicted = interval_direction([(12.0, 12.1)], [(7.0, 7.1)])        # model secretes more
        self.assertEqual(measured, 1)
        self.assertEqual(predicted, 1)
        self.assertIs(agreement(predicted, measured), True)

    def test_uptake_increase_is_a_negative_direction_on_both_sides(self):
        measured = measured_direction([-800.0, -700.0], [-400.0, -300.0])   # more glucose in
        predicted = interval_direction([(-8.0, -7.0)], [(-4.0, -3.0)])
        self.assertEqual(measured, -1)
        self.assertEqual(predicted, -1)
        self.assertIs(agreement(predicted, measured), True)


if __name__ == "__main__":
    unittest.main(verbosity=1)
