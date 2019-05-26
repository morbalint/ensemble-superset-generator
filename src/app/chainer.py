import prody
import pylab

def appendDiamid2Chain(chain, nextDiamid, idx):
    chainResidues = list(chain.iterResidues())
    last_center = chainResidues[-2]
    last2align = chainResidues[-1]
    last_center_filtered = last_center.select('name CA C O')
    last2align_filtered = last2align.select('name N CA')
    chain_coords_2B_aligned = (last_center_filtered + last2align_filtered).getCoords()
    nextResidues = list(nextDiamid.iterResidues())
    next_prev = nextResidues[0].select('name CA C O')
    next_center = nextResidues[1]
    next_center.setResnum(idx)
    next_center.setChids('A')
    next_center_atoms_2B_aligned = next_center.select('name N CA')
    next_next = nextResidues[2]
    next_next.setResnum(idx + 1)
    next_next.setChids('A')
    next_coords_2B_aligned = (next_prev + next_center_atoms_2B_aligned).getCoords()
    if (len(next_coords_2B_aligned) != 5 or len(chain_coords_2B_aligned) != 5):
        print('insuficent atoms in one of these: ')
        print('next:')
        print(next_coords_2B_aligned)
        print(nextDiamid)
        print('chain:')
        print(chain_coords_2B_aligned)
        raise Exception('insuficent atoms')

    tran = prody.calcTransformation(next_coords_2B_aligned, chain_coords_2B_aligned)
    next_2B_attached = next_center + next_next
    next_aligned = prody.applyTransformation(tran, next_2B_attached)
    chain_without_last = chain.select('not resnum ' + str(last2align.getResnum())).copy()
    return chain_without_last + next_aligned.copy()

if __name__ == '__main__':
    dap = prody.parsePDB('samples/DAP-1c9k-A-22_23_24.pdb')
    aaa = prody.parsePDB('samples/diamides/AAA-1af7-A-122_123_124.pdb')
    chain = appendDiamid2Chain(dap, aaa, 2)
    chain = appendDiamid2Chain(chain, dap, 3)
    chain = appendDiamid2Chain(chain, aaa, 4)
    chain = appendDiamid2Chain(chain, dap, 5)
    print(chain)
    print(chain.getCoords())
    print([aa.getResname() for aa in chain.iterResidues()])