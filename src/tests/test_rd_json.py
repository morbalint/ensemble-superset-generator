import unittest
import json
import numpy as np

class TestJsonLoad(unittest.TestCase):
    
    def __setup(self):
        with open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8') as f:
            x = json.load(f)
            return x

    def test_file_load(self):
        f = open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8')
        x = json.load(f)
        self.assertIsNotNone(x, 'Should not be None')
        pass

    def test_angle_resolutions(self):
        x = self.__setup()
        self.assertAlmostEqual(x['phi_resolution'], 5.0, msg='Phi resolution should be 5')
        self.assertAlmostEqual(x['psi_resolution'], 5.0, msg='Psi resolution should be 5')
        pass

    def test_number_of_aa(self):
        x = self.__setup()
        self.assertGreaterEqual(len(x['probability']), 20, \
            'There should be at least 20 amino acids in the probability matrix')
        self.assertGreaterEqual(len(x['probability']['ALA']), 20, \
            'There should be at least 2 types of neighboring amino acids')
        aa_num = len(x['probability'])
        for key, val in x['probability'].items():
            self.assertEqual(len(val),aa_num, 'Neighbor matrix should be square. It is not for '+key)

    def test_probability_lengths(self):
        x = self.__setup()
        NX = int(360 / x['phi_resolution'])
        NY = int(360 / x['psi_resolution'])
        for neighbors in x['probability'].values():
            for probability_matrix in neighbors.values():
                self.assertEqual(len(probability_matrix), NY)
                for row in probability_matrix:
                    self.assertEqual(len(row), NX)

    def test_probability_sums(self):
        x = self.__setup()
        for key, val in x['probability'].items():
            for nkey, nval in val.items():
                # somehow these were not even part of the original dataset, but we json deserialized it
                # to build a correct square neighborhood matrix
                if (key != 'ALA' and nkey != 'SEC') and \
                    (nkey != 'ALA' and key != 'SEC'):
                    self.assertAlmostEqual(np.sum(np.sum(nval)), 1.0, \
                        msg='Sum of probabilities should be 1. It is not for ' + key + ' and ' + nkey)

if __name__ == '__main__':
    unittest.main()