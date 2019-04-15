import json
import sys
import numpy as np
from tree_algorithms import ordered_last_search, safe_sum

class RDist:
    """discrete 2 dimensional probability distributrion, with extra data precalculated for speed (memory sacrifice)"""
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.probabilities = probabilities
        flat_prob = [val for row in probabilities for val in row]
        self.cumulative_sums = [0.0] + [ safe_sum(flat_prob[:i]) for i in range(1,len(flat_prob)) ]
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
    def __init__(self, phiResolution, psiResolution, isRightNeighboor):
        self.phiResolution = phiResolution
        self.psiResolution = psiResolution
        self.NX = int(360/phiResolution)
        self.NY = int(360/psiResolution)
        self.isRightNeighboor = isRightNeighboor

class TDRD:
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.distributions = {}
        for key, val in probabilities.items():
            print('calculating cumulative sums of ' + key)
            self.distributions[key] = RDist(meta, val)
        pass

class NDRD:
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.distributions = {}
        for key, val in probabilities.items():
            print('calculating neighboors of ' + key)
            self.distributions[key] = TDRD(meta, val)

def json_test():
    f = open('./data/NDRD/NDRD_R_TCBIG.json', 'r', encoding='utf8')
    x = json.load(f)
    print('phi resolution: ' + str(x['phi_resolution']))
    print('psi resolution: ' + str(x['psi_resolution']))
    print('number of main aa:' + str(len(x['probability'])))
    print('number of neighboring aa:' + str(len(x['probability']['ALA'])))
    print('probability matrix len x:' + str(len(x['probability']['ALA']['ALA'])))
    print('probability matrix len y:' + str(len(x['probability']['ALA']['ALA'][0])))
    print('probability of ALA-ALA 5,5:' + str(x['probability']['ALA']['ALA'][36][35]))
    pass

def TDRD_test():
    data = json.load(open('samples/TDRD_R_TCBIG.json'))
    meta = RDMeta(data['phi_resolution'], data['psi_resolution'], data['right_neighboor'])
    TDRD_R_TCBIG = TDRD(meta, data['probability'])
    print(meta.NX)  # test computed property
    ALA = TDRD_R_TCBIG.distributions['ALA']
    print(ALA.cumulative_sums[-1] + ALA.probabilities[meta.NX-1][meta.NY-1]) # cummulative sum should be close to 1
    print(ALA.draw(isDebug=True))  # test draw function

def NDRD_test():
    data = json.load(open('data/NDRD/NDRD_R_TCBIG.json', 'r', encoding='utf8'))
    meta = RDMeta(data['phi_resolution'], data['psi_resolution'], data['right_neighboor'])
    NDRD_R_TCBIG = NDRD(meta, data['probability'])
    print(meta.NX)  # test computed property
    ALA = NDRD_R_TCBIG.distributions['ALA'].distributions['ALA']
    print(ALA.cumulative_sums[-1] + ALA.probabilities[meta.NX-1][meta.NY-1]) # cummulative sum should be close to 1
    print(ALA.draw(isDebug=True))  # test draw function
