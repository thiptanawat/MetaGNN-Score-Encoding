import unittest
from agreement_audit import decisions, oriented, profile_pairs, EPSILON

def points(values):return [(v,v,v) for v in values]

class AgreementTests(unittest.TestCase):
    def test_zero_is_abstention(self):
        r=decisions(points([0]*5),points([0]*5))
        self.assertEqual(r['point_category'],'unanimous_ties')
        self.assertEqual((r['matched_interval_sign'],r['independent_interval_sign']),(0,0))
    def test_opposing_arms(self):
        r=decisions(points([1,-1]),points([0,0]))
        self.assertEqual(r['point_category'],'opposing_nonzero_signs')
        self.assertEqual(r['matched_interval_sign'],0)
    def test_matched_is_weaker_than_independent(self):
        r=decisions(points([10,0]),points([9,-1]))
        self.assertEqual(r['point_category'],'unanimous_positive')
        self.assertEqual((r['matched_interval_sign'],r['independent_interval_sign']),(1,0))
    def test_both_support_and_reversal(self):
        a=[(2,1.5,2.5),(3,2.5,3.5)];b=[(0,-.5,.5),(0,-.1,.1)]
        for first,second,sign in [(a,b,1),(b,a,-1)]:
            r=decisions(first,second)
            self.assertEqual((r['matched_interval_sign'],r['independent_interval_sign']),(sign,sign))
    def test_point_unanimity_does_not_prove_range_support(self):
        r=decisions([(1,-1,2)]*5,[(0,-2,1)]*5)
        self.assertEqual(r['point_category'],'unanimous_positive')
        self.assertEqual(r['matched_interval_sign'],0)
    def test_mixed_ties(self):
        self.assertEqual(decisions(points([0,1]),points([0,0]))['point_category'],'mixed_ties_one_direction')
    def test_exact_threshold_is_tie(self):
        self.assertEqual(decisions(points([EPSILON]),points([0]))['point_category'],'unanimous_ties')
    def test_orientation_and_small_inversion(self):
        self.assertEqual(oriented(2,1,3,-1),((-2,-3,-1),False))
        x,inverted=oriented(0,1e-8,-1e-8,1)
        self.assertTrue(inverted);self.assertEqual(x,(0,0,0))
    def test_reject_invalid_range_or_point(self):
        for args in [(0,.1,-.1,1),(1,0,0,1),(0,0,0,0),(float('nan'),0,0,1)]:
            with self.assertRaises(ValueError):oriented(*args)
    def test_dependent_origins_excluded(self):
        self.assertEqual(profile_pairs({'a':'o1','b':'o1','c':'o2'}),[('a','c'),('b','c')])
    def test_consistency_failure_is_not_hidden(self):
        with self.assertRaises(ValueError):decisions([(0,1,1)],[(0,0,0)])

if __name__=='__main__':unittest.main()
