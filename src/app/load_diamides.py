import sys
import os
import numpy as np
import rd

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

    def get_angles(self):
        return [self.aa_type,self.phi,self.psi];

    @staticmethod
    def parse_line(line):
        """
        @summary builds an AAResidue instance from a single line of text.
        example line: 'REMARK PDB ID 1c9k residue_number 22 chain A residue_name D phi -88.5 psi 35.0 secondary_structure S'
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


class AAResidue_db:
    
    def __init__(self, DiamideDb, resol_num):
        """
        Létrehozom az aminosav adatbázis osztályt
        """

        self.dict = {'A':"ALA",'C':"CYS",'D':"ASP",'E':"GLU",'F':"PHE",
                     'G':"GLY",'H':"HIS",'I':"ILE",'K':"LYS",'L':"LEU",
                     'M':"MET",'N':"ASN",'P':"PRO",'Q':"GLN",'R':"ARG",
                     'S':"SER",'T':"THR",'V':"VAL",'W':"TRP",'Y':"TYR",
                     'U':"SEC",'O':"PYL"}

        self.db = []
        self.DiamideDb = DiamideDb
        self.resol_num = resol_num
        self.bins = np.arange(180,-180+(-360/resol_num),-360/resol_num)

        for aa in DiamideDb.db: #Beparsolom a szükséges adatokat a Diamide adatbázisból 

            Dangles = aa.get_DiamideAngles()

           # for ind in range(3):

            l = Dangles[1]
            l[1] = self.bins[np.digitize(l[1], self.bins)]
            l[2] = self.bins[np.digitize(l[2], self.bins)]
            #Dangles[1] = l

            #self.db = self.db + l
            self.db.append(l)
            
    def find_AA(self, neve, phi, psi):
        """
        Megadsz neki egy aminosav betűt és a két szöget és ő keres neked egyet az adatbázisban
        Ha nem talál,akkor visszatér None értékkel
        """

        counter = 0
        residue = 0
        integer = 0
        aa = None

        for line in self.db:
            if line[0] == neve and line[1] == phi and line[2] == psi:
                    residue = counter%3
                    integer = int(counter/3)
                    Diamide = self.DiamideDb.db[integer]
                    if residue == 0:
                        aa = Diamide.left_aa
                    elif residue == 1:
                        aa = Diamide.central_aa
                    elif residue == 2:
                        aa = Diamide.right_aa
                    break
            counter += 1

        return aa
        
        
    def query(self,filename,RDtype,trying_num,aa_type1 = None, aa_type2 = None):
        """Mafa a kompakt lekérdezés."""

        aa = None

        if aa_type1 != None and aa_type2 == None:
            RD = rd.get_RD(filename,RDtype)
            daa = RD.distributions[self.dict[aa_type1]]
        elif aa_type1 != None and aa_type2 != None:
            RD = rd.get_RD(filename,RDtype)
            daa = RD.distributions[aa_type1].distributions[aa_type2]
        else
            return aa

        for _ in np.arange(trying_num):
            if aa == None :
                
                [x,y] = daa.draw()
                print(self.bins[np.digitize(x, self.bins)],self.bins[np.digitize(y, self.bins)])
                aa = self.find_AA(aa_type1,self.bins[np.digitize(x, self.bins)],self.bins[np.digitize(y, self.bins)])
            else:
                return aa
        return aa
        
class Diamide:
    """A triplet of amino acids with just enough data to calculate dihedral angles of the central amino acid."""
    def __init__(self, left, central, right):
        self.left_aa = left
        self.central_aa = central
        self.right_aa = right
        self.pdb_id = central.pdb_id
        self.chain = central.chain

    def get_DiamideAngles(self):
        l = self.left_aa.get_angles()
        c = self.central_aa.get_angles()
        r = self.right_aa.get_angles()
        return [l,c,r]

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
        return Diamide(left_aa, central_aa, right_aa)

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
            for line in fl:
                fpath = os.path.join(folder, line)
                #print(line)
                self.db.append(Diamide.parse_file(fpath.rstrip()))
