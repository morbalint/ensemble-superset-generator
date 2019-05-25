import os
import sys

import numpy as np
sys.path.append(os.getcwd()) # this is fuking ridiculus
import src.app.load_diamides as ld
import src.app.rd as rd
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