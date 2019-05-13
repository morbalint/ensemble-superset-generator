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
            line[12:17].strip(),
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
        
        
    def get_angles(self):
        
        return [self.aa_type,self.phi,self.psi]
        
        
    def check_backbone(self):
        
        if len(self.atoms) < 1 :
            
            return 0;

        else:
            return 1;
        

        
    def AA2At_group(self,ID):
        
       group = pd.AtomGroup(self.aa_type)
       coordinates = []
       atom_types = []
       Resname = []
       ResIDs = []
       
       for atom in self.atoms:
           
           coordinates.append([atom.x,atom.y,atom.z])
           atom_types.append(atom.name)
           Resname.append(self.aa_type)
           ResIDs.append(ID)
           
           
           
       group.setCoords(np.asarray(coordinates))
       group.setNames(atom_types)
       group.setResnames(Resname)
       group.setResnums(ResIDs)
       
       return group;




class AAResidue_db:

    def __init__(self, DiamideDb, resol_num):  

        self.dict = {'A': "ALA", 'C': "CYS", 'D': "ASP", 'E': "GLU", 'F': "PHE",
                     'G': "GLY", 'H': "HIS", 'I': "ILE", 'K': "LYS", 'L': "LEU",
                     'M': "MET", 'N': "ASN", 'P': "PRO", 'Q': "GLN", 'R': "ARG",
                     'S': "SER", 'T': "THR", 'V': "VAL", 'W': "TRP", 'Y': "TYR",
                     'U': "SEC", 'O': "PYL"}

        self.db = []
        self.DiamideDb = DiamideDb
        self.resol_num = resol_num
        self.bins = np.arange(180, -180 + (-360 / resol_num), -360 / resol_num)

        for aa in DiamideDb.db:  

            Dangles = aa.get_DiamideAngles()

            # for ind in range(3):

            l = Dangles[1]
            l[1] = self.bins[np.digitize(l[1], self.bins)]
            l[2] = self.bins[np.digitize(l[2], self.bins)]
            # Dangles[1] = l

            # self.db = self.db + l
            self.db.append(l)

    def find_AA(self, neve, phi,
                psi):  

        counter = 0
        

        aa = None
        diamide = None

        for line in self.db:

            if line[0] == neve and line[1] == phi and line[2] == psi:

                diamide = self.DiamideDb.db[counter]
                aa = diamide.central_aa
                
                break

            counter += 1

        return aa,diamide

    def query(self, filename, RDtype, trying_num,ID, aa_type1=None, aa_type2=None):  # Mafa a kompakt lekerdezes.

        aa = None
        diamide = None

        if aa_type1 != None and aa_type2 == None:

            RD = rd.get_RD(filename, RDtype)
            daa = RD.distributions[self.dict[aa_type1]]

        elif aa_type1 != None and aa_type2 != None:

            RD = rd.get_RD(filename, RDtype)
            daa = RD.distributions[aa_type1].distributions[aa_type2]

        for counter in np.arange(trying_num):

            if aa == None:

                [x, y] = daa.draw()
                print(self.bins[np.digitize(x, self.bins)], self.bins[np.digitize(y, self.bins)])
                aa,diamide = self.find_AA(aa_type1, self.bins[np.digitize(x, self.bins)], self.bins[np.digitize(y, self.bins)])
            else:

                return aa,diamide.diamide2AtGroup(ID)

        return aa,diamide

class Diamide:
    """A triplet of amino acids with just enough data to calculate dihedral angles of the central amino acid."""

    def __init__(self, left, central, right,filePath):
        self.left_aa = left
        self.central_aa = central
        self.right_aa = right
        self.pdb_id = central.pdb_id
        self.chain = central.chain
        self.filePath = filePath


    def get_DiamideAngles(self):
        
        l = self.left_aa.get_angles()
        c = self.central_aa.get_angles()
        r = self.right_aa.get_angles()
        
        return [l,c,r];

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


    def check_Diamide_backbone(self):
        
        if self.left_aa.check_backbone() and self.central_aa.check_backbone() and self.right_aa.check_backbone() :
            
            return 1;
        else:
            return 0;
        
        

    def diamide2AtGroup(self,ID):
        
        group  =  self.left_aa.AA2At_group(ID-1)
        group += self.central_aa.AA2At_group(ID)
        group += self.right_aa.AA2At_group(ID+1)
        
        
        return group;
        
        


class DiamidesDb:
    """A class containg a large data set of diamides 
    in some internally optimized format (not yet determined, private to outside users)
    and helper methods to execite simple search queries"""

    def __init__(self, listFile, folder=None):
        """Instantiate a DiamideDb class with the given list of files. 
        Actual input parameter: a file path containing line separated diamide files paths
        """
        self.db = []
        folder = folder or os.path.dirname(listFile)
        with open(listFile) as fl:
            for fpath in [fpath for fpath in [os.path.join(folder, line).rstrip() for line in fl] if fpath.endswith('.pdb')]:
                diamide = Diamide.parse_file(fpath)
                if diamide.check_Diamide_backbone():
                    self.db.append(Diamide.parse_file(fpath))