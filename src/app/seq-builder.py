import os
import sys
sys.path.append(os.getcwd())  # this is fuking ridiculus

import src.app.chainer as chainer
import src.app.load_diamides as ld
import prody

def structure_builder(sequence, db, fileName, RDType):
    _, chain = db.query(fileName, RDType, 10, 1, sequence[0])
    for i in range(len(sequence) - 1):
        aa =sequence[i+1]
        _, diamid = db.query(fileName, RDType, 1, i+2, aa)
        chain = chainer.appendDiamid2Chain(chain, diamid, i + 2)
    chain = chain.select('not resnum 0')
    chain = chain.select('not resnum ' + str(len(sequence)))
    return chain

if __name__ == '__main__':
    db = ld.AAResidue_db(ld.DiamidesDb('data/list_of_diamides.txt', 'data/Database-of-diamides-2018-12/'), 5)
    prody.fetchPDB('1d3z')
    pdb = prody.parsePDB('1d3z.pdb')
    structure = structure_builder(pdb.getSequence(), db, 'samples/TDRD_R_TCBIG.json', 'TDRD')
    prody.writePDB('1d3z_test.pdb', structure)
    print(structure)
    print(structure.getCoords())
    #print([aa.getResname() for aa in structure.iterResidues()])
