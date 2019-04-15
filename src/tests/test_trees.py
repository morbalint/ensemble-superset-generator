import unittest
import json
import numpy as np
import sys
from src.app.tree_algorithms import ordered_last_search, safe_sum

class TestSafeSum(unittest.TestCase):
    def test_summing_ones(self):
        self.assertEqual(safe_sum(np.ones(1024)), 1024)
        self.assertEqual(safe_sum(np.ones(77)), 77)

    def test_my_bin_search(self):
        arr = range(128)
        i = ordered_last_search(arr, lambda x: x < 42)
        self.assertEqual(i, 41)

    def test_sum_probabilites(self):
        with open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8') as f:
            data = json.load(f)
            probs = data['probability']['ALA']['ALA']
            one = safe_sum([safe_sum(row) for row in probs])
            flat = [val for row in probs for val in row]
            one_again = safe_sum(flat)
            # ASK: numerical precision is off even in the reference database
            self.assertAlmostEqual(one, 1.0, 7)
            self.assertAlmostEqual(one_again, 1.0, 7)

if __name__ == '__main__':
    unittest.main()