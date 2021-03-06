import sys
import os
import numpy as np
import src.app.rd as rd
import prody as pd

class Atom:
    def __init__(self, id, name, x, y, z, atom_type=None, aa=None):
        self.id = id
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.atom_type = atom_type or name[0]
        self.aa = aa  # link to amino acid

    @staticmethod
    def parse_line(line, aa=None):
        return Atom(
            int(line[4:12]),
            line[12:16].strip(),
            float(line[31:38]),
            float(line[39:46]),
            float(line[47:55]),
            line[77:].strip(),
            aa
        )

class AAResidue:
    """ Amino acid residue (short notation: AA) """

    def __init__(self, pdb_id, res_num, chain, aa_type, phi, psi, sec_structure=None):
        self.pdb_id = pdb_id
        self.chain = chain
        self.res_num = res_num
        self.aa_type = aa_type
        self.phi = phi
        self.psi = psi
        self.sec_structure = sec_structure
        self.atoms = []

    @staticmethod
    def parse_line(line):
        """
        @summary builds an AAResidue instance from a single line of text. example line: 'REMARK PDB ID 1c9k
        residue_number 22 chain A residue_name D phi -88.5 psi 35.0 secondary_structure S'
        """
        vals = line.split()
        pdb_id = vals[3]
        res_num = int(vals[5])
        chain = vals[7]
        res_name = vals[9]
        phi = float(vals[11])
        psi = float(vals[13])
        sec_structure = None
        if len(vals) > 15:
            sec_structure = vals[15]
        return AAResidue(pdb_id, res_num, chain, res_name, phi, psi, sec_structure)

class Diamide:
    """A triplet of amino acids with just enough data to calculate dihedral angles of the central amino acid."""

    def __init__(self, left, central, right,filePath):
        self.left_aa = left
        self.central_aa = central
        self.right_aa = right
        self.pdb_id = central.pdb_id
        self.chain = central.chain
        self.filePath = filePath

    @staticmethod
    def parse_file(filePath):
        """
        Parses a file containing a diamide.
        Input parameter is the file path.
        Output is an instance of the diamide class.
        (File format given by Zita Harmati)"""
        lines = []
        with open(filePath) as f:
            lines = f.readlines()
            
        left_aa = AAResidue.parse_line(lines[0])
        central_aa = AAResidue.parse_line(lines[1])
        right_aa = AAResidue.parse_line(lines[2])
        for line in lines[3:]:
            seq_num = int(line[22:27])
            aa = None
            if left_aa.res_num == seq_num:
                aa = left_aa
            elif central_aa.res_num == seq_num:
                aa = central_aa
            elif right_aa.res_num == seq_num:
                aa = right_aa
            else:
                print('unkown sequence: %i' % seq_num)
            atom = Atom.parse_line(line, aa)
            aa.atoms.append(atom)
        return Diamide(left_aa, central_aa, right_aa, filePath)
