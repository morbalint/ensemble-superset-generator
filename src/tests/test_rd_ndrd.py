import unittest
import json
import numpy as np
import sys
import src.app.rd as rd

class TestRDMeta(unittest.TestCase):
    def __setup(self):
        with open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8') as f:
            self.data = json.load(f)
            self.meta = rd.RDMeta(
                self.data['phi_resolution'],
                self.data['psi_resolution'],
                self.data['right_neighboor'])
            return self.meta

    def test_meta(self):
        """yes this is not really a good test... :("""
        meta = self.__setup()
        self.assertEqual(meta.NX, 72,
            'NX should be exactly 360°/5.0° = 72 (full circle / phi resolution)')
        self.assertEqual(meta.NY, 72,
            'NY should be exaclty 360°/5.0° = 72 (full circle / psi resolution)')
        self.assertIsInstance(meta.isRightNeighbor, bool,
            'neighbor type should be boolean')

class TestRDist(unittest.TestCase):
    def __setup(self):
        with open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8') as f:
            self.data = json.load(f)
            self.meta = rd.RDMeta(
                self.data['phi_resolution'],
                self.data['psi_resolution'],
                self.data['right_neighboor'])
            self.ALA_ALA = rd.RDist(self.meta, self.data['probability']['ALA']['ALA'])
            return self.ALA_ALA

    def test_cum_sum(self):
        dist = self.__setup()
        self.assertEqual(dist.cumulative_sums[0], 0.0, \
            'the fist element of a matrix of cumulative sums should be 0')
        self.assertAlmostEqual(dist.cumulative_sums[-1] + dist.probabilities[-1][-1], 1.0,
            msg='sum of probailities should be 1')
        for i in range(1,len(dist.cumulative_sums)):
            self.assertGreaterEqual(dist.cumulative_sums[i], dist.cumulative_sums[i-1],
                'Cumulative sums should form a monotonly increasing series! i='+str(i))

    def test_random_draw(self):
        dist = self.__setup()
        (x,y) = dist.draw(isDebug=True)
        self.assertLess(x, dist.meta.NX)
        self.assertLess(y, dist.meta.NY)
        self.assertGreaterEqual(x, 0)
        self.assertGreaterEqual(y, 0)

    def test_draw_zero(self):
        dist = self.__setup()
        (x,y) = dist.draw(rnd=0.0, isDebug=False)
        self.assertEqual(x, 0)
        self.assertEqual(y, 0)

    def test_draw_one(self):
        dist = self.__setup()
        (x,y) = dist.draw(rnd=1.0, isDebug=False)
        self.assertEqual(x, self.meta.NX-1)
        self.assertEqual(y, self.meta.NY-1)

if __name__ == '__main__':
    unittest.main()