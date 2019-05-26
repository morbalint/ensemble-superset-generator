import os
import sys
sys.path.append(os.getcwd())  # this is fuking ridiculus

import prody

import src.app.chainer as chainer
import src.app.db as db


def structure_builder(sequence, thedb, isNeighborDependent=False, isRightNeighbor=False):
    chain = thedb.query(sequence[0])
    # format chain starting block
    temp = list(chain.iterResidues())
    temp[0].setResnum(0)
    temp[1].setResnum(1)
    temp[0].setChids('A')
    temp[1].setChids('A')
    # build chain
    for i in range(len(sequence) - 1):
        neighborAA = None
        if (isNeighborDependent):
            if (isRightNeighbor):
                if(i+2 < len(sequence)):
                    neighborAA = sequence[i + 2]
            else:
                neighborAA = sequence[i]
        aa = sequence[i+1]
        diamid = thedb.query(aa, neighborAA)
        chain = chainer.appendDiamid2Chain(chain, diamid, i + 2)
    chain = chain.select('not resnum 0').copy()
    chain = chain.select('not resnum ' + str(len(sequence)+1)).copy()
    return chain

if __name__ == '__main__':
    prody.confProDy(verbosity='warning')
    prody.fetchPDB('1d3z')
    pdb = prody.parsePDB('1d3z.pdb.gz')
    sequence = pdb.select('name CA').getSequence()
    thedb = db.DB('samples/TDRD_R_TCBIG.json', 'data/diamides', 'samples/NDRD_R_TCBIG_pretty.json')
    N = 10000
    lenN = str(len(str(N)))
    with open('1d3z_test_ndrd_out.pdb', 'a+') as f:
        f.write('NUMMDL\t%d\n' % N)
        for i in range(N):
            structure = structure_builder(sequence, thedb, True, True)
            f.write('MODEL\t%d\n' % i)
            prody.writePDBStream(f, structure)
            f.write('ENDMDL\n')
            print(('structure {:'+lenN+'d} / {:d} generated').format(i, N))
