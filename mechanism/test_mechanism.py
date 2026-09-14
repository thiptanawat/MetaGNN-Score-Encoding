"""Known-answer checks for mechanism.py, written before any study output was read.

Each case states an answer that can be verified by hand. A reversed comparator or a sign
convention read the wrong way round would fail here rather than reaching a conclusion.
"""
import math, unittest
import numpy as np
from mechanism import (WIDTH_TOL, width, location, is_fixed, spread, exploration,
                       linear_cost, regret, relative_l1)


class Geometry(unittest.TestCase):
    def test_width_is_high_minus_low(self):
        self.assertEqual(width(2.0, 5.0), 3.0)

    def test_width_of_a_consumption_interval_is_still_positive(self):
        # negative flux means uptake; the interval [-5, -2] is 3 units wide, not -3
        self.assertEqual(width(-5.0, -2.0), 3.0)

    def test_width_of_an_interval_spanning_zero(self):
        self.assertEqual(width(-1.5, 2.5), 4.0)

    def test_location_is_the_midpoint(self):
        self.assertEqual(location(2.0, 6.0), 4.0)
        self.assertEqual(location(-6.0, -2.0), -4.0)
        self.assertEqual(location(-3.0, 3.0), 0.0)

    def test_fixedness_threshold_is_inclusive(self):
        self.assertTrue(is_fixed(0.0, WIDTH_TOL))          # exactly at tolerance counts as fixed
        self.assertFalse(is_fixed(0.0, WIDTH_TOL * 1.001))
        self.assertTrue(is_fixed(-4.0, -4.0))              # degenerate negative interval

    def test_a_wide_interval_at_a_negative_location_is_not_fixed(self):
        self.assertFalse(is_fixed(-10.0, -1.0))


class Spread(unittest.TestCase):
    def test_spread_is_max_minus_min(self):
        self.assertEqual(spread([1.0, 4.0, -2.0]), 6.0)

    def test_identical_profiles_have_zero_spread(self):
        self.assertEqual(spread([7.0, 7.0, 7.0]), 0.0)

    def test_single_profile_has_zero_spread(self):
        self.assertEqual(spread([3.0]), 0.0)

    def test_empty_has_no_spread(self):
        self.assertTrue(math.isnan(spread([])))


class Exploration(unittest.TestCase):
    def test_profiles_pinned_to_one_point_explore_nothing(self):
        # every profile sits at 4.0 inside a shared admissible range 0..10
        self.assertEqual(exploration([4.0, 4.0, 4.0], 0.0, 10.0), 0.0)

    def test_profiles_spanning_the_whole_admissible_range_explore_all_of_it(self):
        self.assertEqual(exploration([0.0, 10.0], 0.0, 10.0), 1.0)

    def test_half_the_range_is_one_half(self):
        self.assertEqual(exploration([2.0, 7.0], 0.0, 10.0), 0.5)

    def test_undefined_when_the_shared_range_is_itself_fixed(self):
        self.assertTrue(math.isnan(exploration([0.0, 0.0], 3.0, 3.0)))

    def test_works_on_a_consumption_range(self):
        self.assertEqual(exploration([-8.0, -2.0], -10.0, 0.0), 0.6)


class Cost(unittest.TestCase):
    def test_cost_uses_absolute_flux(self):
        # uptake and secretion of the same magnitude cost the same
        self.assertEqual(linear_cost([2.0, -2.0], [1.0, 1.0]), 4.0)

    def test_cost_is_weighted(self):
        self.assertEqual(linear_cost([2.0, 3.0], [0.5, 2.0]), 1.0 + 6.0)

    def test_zero_flux_costs_nothing(self):
        self.assertEqual(linear_cost([0.0, 0.0], [5.0, 7.0]), 0.0)


class Regret(unittest.TestCase):
    def test_the_optimal_vector_has_zero_regret(self):
        self.assertEqual(regret([2.0, 0.0], [1.0, 1.0], 2.0), 0.0)

    def test_a_worse_vector_has_positive_regret(self):
        # cost 4 against an optimum of 2 is 100 percent excess
        self.assertAlmostEqual(regret([4.0, 0.0], [1.0, 1.0], 2.0), 1.0)

    def test_a_vector_cheaper_than_the_stated_optimum_reports_negative(self):
        # this must not be silently clipped: it would mean the stored optimum is wrong
        self.assertAlmostEqual(regret([1.0, 0.0], [1.0, 1.0], 2.0), -0.5)

    def test_numerical_noise_at_the_optimum_reads_as_exactly_zero(self):
        self.assertEqual(regret([2.0 + 1e-12, 0.0], [1.0, 1.0], 2.0), 0.0)

    def test_undefined_against_a_zero_optimum(self):
        self.assertTrue(math.isnan(regret([1.0], [1.0], 0.0)))


class Distance(unittest.TestCase):
    def test_identical_vectors_are_at_zero_distance(self):
        self.assertEqual(relative_l1([1.0, -2.0, 3.0], [1.0, -2.0, 3.0]), 0.0)

    def test_disjoint_supports_are_at_the_maximum_distance(self):
        self.assertEqual(relative_l1([1.0, 0.0], [0.0, 1.0]), 2.0)

    def test_sign_flip_is_the_maximum_distance(self):
        self.assertEqual(relative_l1([3.0], [-3.0]), 2.0)

    def test_scale_matters_but_stays_bounded(self):
        self.assertAlmostEqual(relative_l1([1.0], [3.0]), 2 * 2.0 / 4.0)

    def test_two_zero_vectors_have_no_defined_distance(self):
        self.assertTrue(math.isnan(relative_l1([0.0, 0.0], [0.0, 0.0])))


class ConstancyLogic(unittest.TestCase):
    """The inference the analysis will draw, exercised on constructed profiles."""

    def test_narrow_intervals_at_one_location_read_as_shared_pinning(self):
        intervals = [(4.0, 4.0 + 1e-9)] * 5          # five profiles, same place
        locs = [location(lo, hi) for lo, hi in intervals]
        self.assertTrue(all(is_fixed(lo, hi) for lo, hi in intervals))
        self.assertLess(spread(locs), WIDTH_TOL)
        self.assertEqual(exploration(locs, 0.0, 10.0), 0.0)

    def test_narrow_intervals_at_different_locations_do_not_read_as_pinning(self):
        intervals = [(1.0, 1.0 + 1e-9), (8.0, 8.0 + 1e-9)]
        locs = [location(lo, hi) for lo, hi in intervals]
        self.assertTrue(all(is_fixed(lo, hi) for lo, hi in intervals))
        self.assertGreater(spread(locs), WIDTH_TOL)      # narrow but NOT constant
        self.assertAlmostEqual(exploration(locs, 0.0, 10.0), 0.7)


if __name__ == "__main__":
    unittest.main(verbosity=2)
