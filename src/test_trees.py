import unittest
import json
import numpy as np
import sys
from tree_algorithms import ordered_last_search, safe_sum

class TestSafeSum(unittest.TestCase):
    def test_summing_ones(self):
        self.assertEqual(safe_sum(np.ones(1024)), 1024)
        self.assertEqual(safe_sum(np.ones(77)), 77)

    def test_my_bin_search(self):
        arr = range(128)
        i = ordered_last_search(arr, lambda x: x < 42)
        self.assertEqual(i, 41)

if __name__ == '__main__':
    unittest.main()