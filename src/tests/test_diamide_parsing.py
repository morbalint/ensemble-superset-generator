import unittest
import json
import numpy as np
import src.app.load_diamides as ld

class TestParisng(unittest.TestCase):

    def test_atom_parsing(self):
        atom = ld.Atom.parse_line(
            'ATOM    157  CA  ALA A  23      50.923  23.451   8.670  1.00 83.18           C ')
        self.assertEqual(atom.id, 157)
        self.assertEqual(atom.name, 'CA')
        self.assertEqual(atom.x, 50.923)
        self.assertEqual(atom.y, 23.451)
        self.assertEqual(atom.z, 8.670)
        self.assertEqual(atom.atom_type, 'C')

    def test_residue_parsing(self):
        res = ld.AAResidue.parse_line('REMARK PDB ID 1c9k residue_number 22 chain A residue_name D phi -88.5 psi 35.0 secondary_structure S')
        self.assertEqual(res.pdb_id, '1c9k')
        self.assertEqual(res.res_num, 22)
        self.assertEqual(res.chain, 'A')
        self.assertEqual(res.aa_type, 'D')
        self.assertEqual(res.phi, -88.5)
        self.assertEqual(res.psi, 35.0)
        self.assertEqual(res.sec_structure, 'S')
    
    def test_file_parsing(self):
        diamide = ld.Diamide.parse_file('samples/DAP-1c9k-A-22_23_24.pdb')
        self.assertEqual(diamide.pdb_id, '1c9k')
        self.assertEqual(diamide.chain, 'A')
        self.assertEqual(diamide.left_aa.aa_type, 'D')
        self.assertEqual(diamide.central_aa.aa_type, 'A')
        self.assertEqual(diamide.right_aa.aa_type, 'P')
        # TODO: test everything. (low prob, high risk)
        self.assertEqual(len(diamide.left_aa.atoms), 3)
        self.assertEqual(len(diamide.central_aa.atoms), 5)
        self.assertEqual(len(diamide.right_aa.atoms), 2)

if __name__ == '__main__':
    unittest.main()