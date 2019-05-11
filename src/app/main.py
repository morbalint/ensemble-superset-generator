import os
import sys
import numpy as np
sys.path.append(os.getcwd()) # this is fuking ridiculus
#import src.app.rd as rd
import rd,load_diamides


def create_ListFile(pwd_full_data,pwd_filename):
    
    list_file = os.listdir(pwd_full_data)
    f = open( pwd+'/list_file.txt', 'w') 
                  # A lista tartalmát egy fileba írom soronként
    for item in list_file:
        
        f.write("%s\n" % item)
        
    f.close()


#def join_AAs()

if __name__ == '__main__':
    
    
    #rd.json_test()
    #rd.TDRD_test()
    #rd.NDRD_test()
   
    pwd = os.getcwd()
    pwd_home = pwd[0:pwd.find("src")]                       #Az a mappa ami tartalmazza a project fileokat
    pwd_data = "samples/diamides"+pwd[pwd.find("src")-1]    #Az a mappa ami tartalmazza a pdb fileokat.
    pwd_full_data = pwd_home+pwd_data                       
    pwd_jsons = pwd_home+"samples/"  
    pwd_list_data = pwd+'/list_file.txt'                                             

    list_file = os.listdir(pwd_full_data)                  
    create_ListFile(pwd_full_data, pwd_list_data)   # Egy listába írom a pdb fileneveket

    
    DATAbase = load_diamides.DiamidesDb(pwd_list_data,pwd_full_data)
    
    
    
    resol_num = 5 ### FIGYELEM EZ AZ ADAT BELE VAN ÉGETVE A JSON FILEBA! EZT KEZELNI KELL
    filename1 = pwd_jsons+"NDRD_R_TCBIG_pretty.json"
    filename2 = pwd_jsons+"TDRD_R_TCBIG.json"
    
    
    AA_db = load_diamides.AAResidue_db(DATAbase,resol_num) #Itt készítem el az aminosav adatbázist

    alma = AA_db.query(filename2,"TDRD",10,'A')  #Itt kérek egy aminosavat az adatbázisból
    
    print("hello world!")