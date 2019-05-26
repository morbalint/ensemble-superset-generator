import os
import numpy as np
import prody

import src.app.rd as rd

aa_map = {  'A': "ALA", 'C': "CYS", 'D': "ASP", 'E': "GLU", 'F': "PHE",
            'G': "GLY", 'H': "HIS", 'I': "ILE", 'K': "LYS", 'L': "LEU",
            'M': "MET", 'N': "ASN", 'P': "PRO", 'Q': "GLN", 'R': "ARG",
            'S': "SER", 'T': "THR", 'V': "VAL", 'W': "TRP", 'Y': "TYR",
            'U': "SEC", 'O': "PYL"  }

class DB:
    def __init__(self, tdrdPath, dataFolder, ndrdPath=None):
        self.tdrd = rd.get_RD(tdrdPath, False)
        self.dataFolder = dataFolder
        if (ndrdPath != None):
            self.isNeighborDependent = True
            self.ndrd = rd.get_RD(ndrdPath, True)
        else:
            self.isNeighborDependent = False


    def query(self, aaType, neighborType=None, maxTries = 100):
        dist = self._get_rd(aaType, neighborType)
        options = []
        counter = 0
        phi = 0
        psi = 0
        while(len(options) < 1 and counter < maxTries):
            phi, psi = dist.draw()
            listFilePath = os.path.join(self.dataFolder, aaType, str(phi), str(psi), 'list.txt')
            with open(listFilePath) as f:
                options = f.readlines()
        if (len(options) < 1):
            raise Exception('can not draw a valid aa')
        selectedIdx = np.random.randint(0, len(options))
        selectedPath = os.path.join(self.dataFolder, aaType, str(phi), str(psi), options[selectedIdx].strip())
        return prody.parsePDB(selectedPath)

    def _get_rd(self, aaType, neighborType=None):
        aa = aa_map[aaType]
        if (self.isNeighborDependent and neighborType != None):
            return self.ndrd.distributions[aa].distributions[aa_map[neighborType]]
        else:
            return self.tdrd.distributions[aa]