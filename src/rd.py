import json
import sys
from collections import namedtuple
from typing import NamedTuple
import numpy as np

class RDist:
    """discrete 2 dimensional probability distributrion, with extra data precalculated for speed (memory sacrifice)"""
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.probabilities = probabilities
        self.cum_sum = []
        self.cum_row_sums = []
        sum = 0.0
        # TODO: we need better precision here, the cum_sum is off by .00002985118158 (~2e-5)
        # which is more than comparable to a single probaility (~1e-7)
        for row in probabilities:
            row_sum = 0.0
            cum_row = []
            self.cum_row_sums.append(sum)
            for value in row:
                row_sum += value
                cum_row.append(sum + row_sum)
            sum += row_sum
            self.cum_sum.append(cum_row)
        pass

    def draw(self, isDebug=False):
        """Draws a random value from the discrete distribution"""
        # TODO: how to create binary tree with custom serch algoritm in python ?
        # WISH: optimize algorithm if we have more time.
        rnd = np.random.uniform()
        row_idx = next(idx for idx, val in enumerate(self.cum_row_sums) if val >= rnd) - 1
        col_idx = next(idx for idx, val in enumerate(self.cum_sum[row_idx]) if val >= rnd) - 1
        if isDebug:
            print(rnd)
            print(self.cum_row_sums[row_idx])
            print(self.cum_row_sums[row_idx+1])
            print(self.cum_sum[row_idx][col_idx])
            print(self.cum_sum[row_idx][col_idx+1])
        return row_idx,col_idx

class RDMeta:
    def __init__(self, phiResolution, psiResolution, isRightNeighboor):
        self.phiResolution = phiResolution
        self.psiResolution = psiResolution
        self.NX = int(360/phiResolution)
        self.NY = int(360/psiResolution)
        self.isRightNeighboor = isRightNeighboor

class TDRD(RDMeta):
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.distributions = {}
        for key,val in probabilities.items():
            self.distributions[key] = RDist(meta, val)
        pass

class NDRD:
    def __init__(self, meta, probabilities):
        self.meta = meta
        self.distributions = {}
        for key,val in probabilities.items():
            # TODO: this is probably wrong
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
    print(ALA.cum_sum[meta.NX-1][meta.NY-1] + ALA.probabilities[meta.NX-1][meta.NY-1]) # cummulative sum should be close to 1
    print(ALA.draw())  # test draw function

def NDRD_test():
    data = json.load(open('data/NDRD/NDRD_R_TCBIG.json', 'r', encoding='utf8'))
    meta = RDMeta(data['phi_resolution'], data['psi_resolution'], data['right_neighboor'])
    NDRD_R_TCBIG = NDRD(meta, data['probability'])
    print(meta.NX)  # test computed property
    ALA = NDRD_R_TCBIG.distributions['ALA'].distributions['ALA']
    print(ALA.cum_sum[meta.NX-1][meta.NY-1] + ALA.probabilities[meta.NX-1][meta.NY-1]) # cummulative sum should be close to 1
    print(ALA.draw())  # test draw function

if __name__ == '__main__':
    # json_test()
    TDRD_test()
    NDRD_test()