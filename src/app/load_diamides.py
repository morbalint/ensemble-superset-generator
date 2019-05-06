import sys
import os

class Atom:
    def __init__(self, id, name, x, y, z, atom_type=None, aa=None):
        self.id = id
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.atom_type = atom_type or name[0]
        self.aa = aa # link to amino acid

    @staticmethod
    def parse_line(line, aa=None):
        return Atom(
            int(line[4:12]),
            line[12:17].strip(),
            float(line[31:39]),
            float(line[39:47]),
            float(line[47:55]),
            line[77:].strip(),
            aa
        )

class AAResidue:
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
        """example line: 'REMARK PDB ID 1c9k residue_number 22 chain A residue_name D phi -88.5 psi 35.0 secondary_structure S'"""
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
    def __init__(self, left, central, right):
        self.left_aa = left
        self.central_aa = central
        self.right_aa = right
        self.pdb_id = central.pdb_id
        self.chain = central.chain

    @staticmethod
    def parse_file(filePath):
        lines = []
        with open(filePath) as f:
            lines = f.readlines()
        left_aa = AAResidue.parse_line(lines[0])
        central_aa = AAResidue.parse_line(lines[1])
        right_aa = AAResidue.parse_line(lines[2])
        print(type(right_aa))
        for line in lines[3:]:
            seq_num = int(line[24:27])
            aa = None
            if left_aa.res_num == seq_num:
                aa = left_aa
            elif central_aa.res_num == seq_num:
                aa = central_aa
            elif right_aa.res_num == seq_num:
                aa = right_aa
            else
                print(seq_num)
            atom = Atom.parse_line(line, aa)
            aa.atoms.append(atom)
        return Diamide(left_aa, central_aa, right_aa)

class DiamidesDb:
    def __init__(self, listFile, folder=None):
        self.db = []
        folder = folder or os.path.dirname(listFile)
        with open(listFile) as fl:
            for line in fl:
                fpath = os.path.join(folder, line)
                self.db.append(Diamide.parse_file(fpath.rstrip()))
