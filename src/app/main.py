import os
import sys

import numpy as np
sys.path.append(os.getcwd()) # this is fuking ridiculus
import load_diamides as ld
import rd as rd
# import rd,load_diamides

def get_pwd_full_data():
    pwd = os.getcwd()
    pwd_home = pwd[0:pwd.find("src")]                     #Az a mappa ami tartalmazza a project fileokat
    pwd_data = "samples/diamides"+os.sep                    #Az a mappa ami tartalmazza a pdb fileokat.
    pwd_full_data = os.path.join(pwd_home,pwd_data)         
    return pwd_full_data                                   

def create_list_file(pwd_full_data):
    list_file = os.listdir(pwd_full_data)                   
    with open( pwd_full_data+'list_file.txt', 'w') as f:   
        for item in list_file:
            if item != 'list_file.txt':
                f.write("%s\n" % item)
        f.close()
def get_pwd_jsons():
    pwd = os.getcwd()
    pwd_home = pwd[0:pwd.find("src")]  
    return pwd_home+"samples/"; 
    
if __name__ == '__main__':
    
    print("hello world!")
   
    pwd = os.getcwd()
    pwd_home = pwd[0:pwd.find("src")]                       #Az a mappa ami tartalmazza a project fileokat
    pwd_data = "samples/diamides"+pwd[pwd.find("src")-1]    #Az a mappa ami tartalmazza a pdb fileokat.
    pwd_full_data = pwd_home+pwd_data                       
    pwd_jsons = pwd_home+"samples/"  
    pwd_list_data = pwd+'/list_file.txt'                                             

    pwd_full_data = get_pwd_full_data()
    create_list_file(pwd_full_data)
    pwd_list_data = pwd_full_data+'list_file.txt'
    pwd_jsons = get_pwd_jsons()
    #load_diamides.DiamidesDb DATAbase
    DATAbase = ld.DiamidesDb(pwd_list_data)
    #(list_file)
    # we must always write our program between hello world printing to remember where we started from.
    resol_num = 5 
    filename1 = pwd_jsons+"NDRD_R_TCBIG_pretty.json"
    filename2 = pwd_jsons+"TDRD_R_TCBIG.json"
    
    
    AA_db = ld.AAResidue_db(DATAbase,resol_num)

    aa,group = AA_db.query(filename2,"TDRD",10,23,'A')  #Returns with the aminoacid with our format and its diamide in atomgroup format 
    #group = diamide.diamide2AtGroup(23)

    print("hello world!")