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
    def __init__(self, rdFilePath, isNeighboorDependent, dataFolder):
        self.isNeighboorDependent = isNeighboorDependent
        self.rd = rd.get_RD(rdFilePath, isNeighboorDependent)
        self.dataFolder = dataFolder

    def query(self, aaType, neighboorType=None, maxTries = 100):
        dist = None
        aa = aa_map[aaType]
        if (self.isNeighboorDependent):
            if(neighboorType == None):
                raise Exception('neighboor aa type missing')
            else:
                dist = self.rd.distributions[aa].distributions[aa_map[neighboorType]]
        else:
            dist = self.rd.distributions[aa]
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
