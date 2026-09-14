import unittest
from nested_range_pilot import normalized_range,nested_violation,stage_label

class NestedLogic(unittest.TestCase):
    def test_nested_stage_attribution(self):
        self.assertEqual(stage_label([0,10],[5,10],[10,10]),'first_fixed_B2')
        self.assertEqual(stage_label([0,10],[5,5],[5,5]),'first_fixed_B1')
        self.assertEqual(stage_label([0,0],[0,0],[0,0]),'already_fixed_B0')
        self.assertEqual(stage_label([0,10],[2,8],[3,7]),'still_variable_B2')
    def test_nesting_catches_outside_endpoints(self):
        self.assertEqual(nested_violation([0,10],[2,8]),0)
        self.assertEqual(nested_violation([0,10],[-1,8]),1)
        self.assertEqual(nested_violation([0,10],[2,12]),2)
    def test_midpoint_is_only_numerical(self):
        self.assertEqual(normalized_range(1e-8,-1e-8),[0,0])
        with self.assertRaises(ValueError):normalized_range(1,-1)

if __name__=='__main__':unittest.main()
