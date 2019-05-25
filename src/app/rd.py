import json
import sys
import numpy as np
from src.app.tree_algorithms import ordered_last_search, safe_sum
#from tree_algorithms import ordered_last_search, safe_sum

class RDist:
    """discrete 2 dimensional probability distributrion, with extra data precalculated for speed (memory sacrifice)"""
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.probabilities = probabilities
        flat_prob = [val for row in probabilities for val in row]
        # safer but much slower method
        # self.cumulative_sums = [0.0] + [ safe_sum(flat_prob[:i]) for i in range(1,len(flat_prob)) ]
        self.cumulative_sums = [0.0] + np.cumsum(flat_prob).tolist()[:-1]
        # testing differences. results are less than 1e-25 squared 'error'
        # step1 = [0.0] + [safe_sum(flat_prob[:i]) for i in range(1, len(flat_prob))]
        # step2 = np.subtract(step1, self.cumulative_sums)
        # print(np.dot(step2, step2))
        pass

    def draw(self, rnd=None, isDebug=False):
        """Draws a random value from the discrete distribution"""
        if rnd == None:
            rnd = np.random.uniform()
        idx = ordered_last_search(self.cumulative_sums, lambda x: x <= rnd)
        x = idx % self.meta.NX
        y = int(idx / self.meta.NX)
        if isDebug:
            print(rnd)
            print(self.cumulative_sums[idx])
            print(self.cumulative_sums[idx+1])
        return x,y

class RDMeta:
    """ ramachandran distribibution meta data """
    def __init__(self, phiResolution, psiResolution, isRightNeighboor):
        self.phiResolution = phiResolution
        self.psiResolution = psiResolution
        self.NX = int(360/phiResolution)
        self.NY = int(360/psiResolution)
        self.isRightNeighboor = isRightNeighboor

class TDRD:
    """ Type dependent ramachandran distribution """
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.distributions = {}
        for key, val in probabilities.items():
            self.distributions[key] = RDist(meta, val)
        pass

class NDRD:
    """ Neighboor dependent ramachandran distribution """ 
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.distributions = {}
        for key, val in probabilities.items():
            self.distributions[key] = TDRD(meta, val)

def get_RD(filename,isNeighboorDependent):
    
    if isNeighboorDependent :
        data = json.load(open(filename, 'r', encoding='utf8'))
        meta = RDMeta(data['phi_resolution'], data['psi_resolution'], data['right_neighboor'])
        return NDRD(meta, data['probability'])
    else:
        data = json.load(open(filename))
        meta = RDMeta(data['phi_resolution'], data['psi_resolution'], data['right_neighboor'])
        return TDRD(meta, data['probability'])
