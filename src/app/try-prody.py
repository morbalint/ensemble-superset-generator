from prody import *
from pylab import *


if __name__ == '__main__':
    print('hello')
    dap = parsePDB('samples/DAP-1c9k-A-22_23_24.pdb')
    print(list(dap['A'])[1].getCoords())