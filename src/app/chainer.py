from prody import *
from pylab import *
import itertools

def nth_item(iterator, n):
    """gets the nth item from an iterator"""
    if n < 0:
        return list(iterator)[n]
        # n = len(iterator) + n
    return next(itertools.islice(iterator, n, n + 1))

def select_chain_coords_2B_aligned(chain):
    last_center = nth_item(chain.iterResidues(), -2)
    last2align = nth_item(chain.iterResidues(), -1)
    return (last_center.select('name CA C O') + last2align.select('name N CA')).getCoords()

def select_diamid_coords_2B_aligned(diamid):
    return 1

def appendDiamid2Chain(chain, nextDiamid, idx):
    last_center = nth_item(chain.iterResidues(), -2)
    last2align = nth_item(chain.iterResidues(), -1)
    last_center_filtered = last_center.select('name CA C O')
    last2align_filtered = last2align.select('name N CA')
    chain_coords_2B_aligned = (last_center_filtered + last2align_filtered).getCoords()
    next_prev = nth_item(nextDiamid.iterResidues(), 0).select('name CA C O')
    next_center = nth_item(nextDiamid.iterResidues(), 1)
    next_center.setResnum(idx)
    next_center_atoms_2B_aligned = next_center.select('name N CA')
    next_next = nth_item(nextDiamid.iterResidues(), 2)
    next_next.setResnum(idx+1)
    next_coords_2B_aligned = (next_prev + next_center_atoms_2B_aligned).getCoords()
    if (len(next_coords_2B_aligned) != 5 or len(chain_coords_2B_aligned) != 5):
        print('insuficent atoms in one of these: ')
        print('next:')
        print(next_coords_2B_aligned)
        print(nextDiamid)
        print('chain:')
        print(chain_coords_2B_aligned)
        raise error('insuficent atoms')

    tran = calcTransformation(next_coords_2B_aligned, chain_coords_2B_aligned)
    next_2B_attached = next_center + next_next
    next_aligned = applyTransformation(tran, next_2B_attached)
    chain_without_last = chain.select('not resnum ' + str(last2align.getResnum())).copy()
    return chain_without_last + next_aligned.copy()

def run_simple_experiment():
    print('hello')
    dap = parsePDB('samples/DAP-1c9k-A-22_23_24.pdb')
    aaa = parsePDB('samples/diamides/AAA-1af7-A-122_123_124.pdb')
    print('==========================')

    print('1c9k center:')
    _1c9k_center = nth_item(dap.iterResidues(), 1)
    _1c9k_center_coords = _1c9k_center.getCoords()[-3:]
    print(_1c9k_center)
    print(_1c9k_center_coords)
    print('1c9k next:')
    _1c9k_next = nth_item(dap.iterResidues(), 2)
    _1c9k_next_coords = _1c9k_next.getCoords()
    print(_1c9k_next)
    print(_1c9k_next_coords)
    print('1c9k align:')
    _1c9k_align = _1c9k_center + _1c9k_next
    _1c9k_align_coords = _1c9k_align.getCoords()[-5:]
    print(_1c9k_align)
    print(_1c9k_align_coords)
    print('==========================')
    
    print('1af7 prev:')
    _1af7_prev = nth_item(aaa.iterResidues(), 0)
    _1af7_prev_coords = _1af7_prev.getCoords()
    print(_1af7_prev)
    print(_1af7_prev_coords)
    print('1af7 center:')
    _1af7_center = nth_item(aaa.iterResidues(), 1)
    _1af7_center_coords = _1af7_prev.getCoords()[0:2]
    print(_1af7_center)
    print(_1af7_center_coords)
    print('1af7 align:')
    _1af7_align = _1af7_prev + _1af7_center
    _1af7_align_coords = _1af7_align.getCoords()[0:5]
    print(_1af7_align)
    print(_1af7_align_coords)
    print('==========================')

    # calulate alignment
    trans = calcTransformation(_1af7_align_coords, _1c9k_align_coords)
    print(trans)

    # apply alignment
    _1af7_center_aligned = applyTransformation(trans, _1af7_center)
    _1af7_center_aligned_coords = _1af7_center_aligned.getCoords()
    print('1af7 center aligned:')
    print(_1af7_center_aligned_coords)
    print('combined:')
    newChain = _1c9k_center.copy() + _1af7_center_aligned.copy()
    print(newChain)
    print(newChain.getCoords())
    print('==========================')

if __name__ == '__main__':
    dap = parsePDB('samples/DAP-1c9k-A-22_23_24.pdb')
    aaa = parsePDB('samples/diamides/AAA-1af7-A-122_123_124.pdb')
    chain = appendDiamid2Chain(dap, aaa, 2)
    chain = appendDiamid2Chain(chain, dap, 3)
    chain = appendDiamid2Chain(chain, aaa, 4)
    chain = appendDiamid2Chain(chain, dap, 5)
    print(chain)
    print(chain.getCoords())
    print([aa.getResname() for aa in chain.iterResidues()])