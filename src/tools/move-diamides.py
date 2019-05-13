import os
import sys
sys.path.append(os.getcwd()) # this is fuking ridiculus

from itertools import groupby
from shutil import copy2, move

import src.app.load_diamides as ld

# NOTES: 
# list group by: https://stackoverflow.com/questions/37568763/python-group-a-list-into-sublists-by-a-equality-of-projected-value/37568870
# 


def select_center(filePath):
    if len(filePath) > 3:
        return filePath[1]
    else:
        return None

def move_diamides_by_center_type(listFilePath, folder, outFolder):
    with open(listFilePath) as f:
        lines_sorted = sorted([line.strip() for line in f.readlines() if line.strip().endswith('.pdb')], key=select_center)
        groups = [list(it) for k, it in groupby(lines_sorted, select_center)]
        groupNames = [ grp[0][1]  for grp in groups if len(grp) > 0 and select_center(grp[0]) != None]
        for grp in groupNames:
            os.makedirs(os.path.join(outFolder, grp), exist_ok=True)
        for grp in groups:
            name = select_center(grp[0])
            folder_2B_placed = os.path.join(outFolder, name)
            for filePath in grp:
                copy2(os.path.join(folder, filePath), os.path.join(folder_2B_placed, filePath))
    pass

def move_diamides_by_angle(listFile, folder, angle='phi', resolution=5):
    N = int(360/resolution)
    for i in range(N):
        os.makedirs(os.path.join(folder, str(i)), exist_ok=True)
    names = [ [] for _ in range(N) ]
    with open(listFile) as fl:
        for line in [l.strip() for l in fl.readlines() if len(l.strip()) > 0]:
            filePath = os.path.join(folder, line)
            file2chech = os.path.join(os.getcwd(), filePath)
            if os.path.isfile(file2chech):
                dia = ld.Diamide.parse_file(filePath)
                if len(dia.central_aa.atoms) > 3:
                    idx = -1
                    if angle == 'phi':
                        idx = int(int(dia.central_aa.phi + 180.0) / resolution) % N
                    elif angle == 'psi':
                        idx = int(int(dia.central_aa.psi + 180.0) / resolution) % N
                    else:
                        print('wrong angle: ' + str(angle))
                    if idx >= 0:
                        move(filePath, os.path.join(os.path.join(folder, str(idx)), line))
                        names[idx].append(line)
                else:
                    # wrong file generated by Zita
                    os.remove(filePath)
    for i, lines in enumerate(names):
        with open(os.path.join(os.path.join(folder, str(i)), 'list.txt'), 'w') as f:
            f.writelines(lines)

def move_all_diamides_by_angle(list_of_lists, folder, resolution=5):
    N = int(360/resolution)
    with open(list_of_lists) as f:
        for lst in [ line.strip() for line in f.readlines() if len(line.strip()) > 0 ]:
            AA = lst[-5]
            currentFolder = os.path.join(folder, AA)
            move_diamides_by_angle(os.path.join(folder, lst), currentFolder, 'phi', resolution)
            for i in range(N):
                phiFolder = os.path.join(currentFolder, str(i))
                move_diamides_by_angle(os.path.join(phiFolder, 'list.txt'), phiFolder, 'psi', resolution)
            print('done with: ' + AA)

if __name__ == "__main__":
    print(select_center('DAP-1c9k-A-22_23_24.pdb'))
    move_diamides_by_center_type('data/list_of_diamides.txt', 'data/Database-of-diamides-2018-12', 'data/diamides')
    move_all_diamides_by_angle('data/diamides/list_lists.txt', 'data/diamides')
    print('hello world')
