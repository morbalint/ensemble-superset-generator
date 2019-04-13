import unittest
import json

class TestJsonLoad(unittest.TestCase):
    
    def __setup(self):
        f = open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8')
        x = json.load(f)
        return x

    def test_file_load(self):
        f = open('samples/NDRD_R_TCBIG_pretty.json', 'r', encoding='utf8')
        x = json.load(f)
        self.assertIsNotNone(x, 'Should not be None')
        pass

    def test__angle__resolutions(self):
        x = self.__setup()
        self.assertAlmostEqual(x['phi_resolution'], 5.0, msg='Phi resolution should be 5')
        self.assertAlmostEqual(x['psi_resolution'], 5.0, msg='Psi resolution should be 5')
        pass

if __name__ == '__main__':
    unittest.main()