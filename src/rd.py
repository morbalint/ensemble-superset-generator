import json
import sys
from collections import namedtuple
from typing import NamedTuple
import numpy as np

def load():
    f = open('./data/NDRD/NDRD_R_TCBIG.json', 'r', encoding='utf8')
    return json.load(f)

class RDist:
	def __init__(self, phi_res, psi_res, probabilities):
		self.phi_res = phi_res
        self.psi_res = psi_res
        self.probabilities = probabilities
        self.NX = 360/phi_res
        self.NY = 360/psi_res
        self.cum_sum = []
        self.cum_row_sums = []
        sum = 0.0
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
        
    def draw():
    	rnd = np.random.uniform()
        row_idx = next(idx for idx,val in enumerate(self.cum_row_sums) if val >= rnd)
        col_idx = next(idx for idx,val in enumerate(self.cum_sum[row_idx]) if val >= rnd)
        return row_idx,col_idx
        
class TDRD:
	def __init__(self, phi_res, psi_res, probabilities):
		self.phi_res = phi_res
		self.psi_res = psi_res
		self.distributions = {}
		for key,val in probabilities.iteritems():
			self.distributions[key] = RDist(phi_res, psi_res, val)
		pass

      
def test():
    x = load()
    print('phi resolution: ' + str(x['phi_resolution']))
    print('psi resolution: ' + str(x['psi_resolution']))
    print('number of main aa:' + str(len(x['probability'])))
    print('number of neighboring aa:' + str(len(x['probability']['ALA'])))
    print('probability matrix len x:' + str(len(x['probability']['ALA']['ALA'])))
    print('probability matrix len y:' + str(len(x['probability']['ALA']['ALA'][0])))
    print('probability of ALA-ALA 5,5:' + str(x['probability']['ALA']['ALA'][36][35]))
    # this does not work yet :(
    # y = NamedTuple(x['phi_resolution'], x['psi_resolution'], x['probability'])
    # return y
    data = json.load(open('samples/TDRD_R_TCBIG.json'))
    TDRD_R_TCBIG = TDRD(data['phi_resolution'], data['psi_resolution'], data['probability'])
    print(TDRD_R_TCBIG.distributions['ALA'].NX) # test computed property
    print(TDRD_R_TCBIG.distributions['ALA'].draw())
    pass
    
    
    
    