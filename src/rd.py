import json
import sys
from collections import namedtuple
from typing import NamedTuple

def load():
    f = open('./data/NDRD/NDRD_R_TCBIG.json', 'r', encoding='utf8')
    return json.load(f)

class RDist(NamedTuple):
    phi_res: int
    psi_res: int
    probability: list

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