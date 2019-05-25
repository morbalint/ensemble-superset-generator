import os
import sys
sys.path.append(os.getcwd())  # this is fuking ridiculus

import prody

import src.app.chainer as chainer
import src.app.db as db


def structure_builder(sequence, thedb, isNeighboorDependent):
    chain = thedb.query(sequence[0])
    for i in range(len(sequence) - 1):
        aa = sequence[i+1]
        diamid = thedb.query(aa)
        chain = chainer.appendDiamid2Chain(chain, diamid, i + 2)
    chain = chain.select('not resnum 0')
    chain = chain.select('not resnum ' + str(len(sequence)))
    return chain

if __name__ == '__main__':
    prody.fetchPDB('1d3z')
    pdb = prody.parsePDB('1d3z.pdb.gz')
    sequence = pdb.select('name CA').getSequence()
    thedb = db.DB('samples/TDRD_R_TCBIG.json', False, 'data/diamides')
    structure = structure_builder(sequence, thedb, False)
    prody.writePDB('1d3z_test.pdb', structure)
    print(structure)
    print(structure.getCoords())
    #print([aa.getResname() for aa in structure.iterResidues()])
