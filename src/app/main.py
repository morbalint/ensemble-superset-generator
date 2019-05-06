import os
import sys
sys.path.append(os.getcwd()) # this is fuking ridiculus
#import src.app.rd as rd
import rd,load_diamides

if __name__ == '__main__':
    
    print("hello world!")
    #rd.json_test()
    #rd.TDRD_test()
    #rd.NDRD_test()
    #print("hello world!")
    pwd = os.getcwd()
    pwd_home = pwd[0:pwd.find("src")]                       #Az a mappa ami tartalmazza a project fileokat
    pwd_data = "samples/diamides"+pwd[pwd.find("src")-1]    #Az a mappa ami tartalmazza a pdb fileokat.
    pwd_full_data = pwd_home+pwd_data                       #A pdb fileok teljes elérési útvonala, igyekeztem ezt dinamikussá tenni
                                                            #hogy mindkettőnk gépén lefusson

    list_file = os.listdir(pwd_full_data)                   # Egy listába írom a pdb fileneveket
    


    with open( pwd_full_data+'list_file.txt', 'w') as f:                   # A lista tartalmát egy fileba írom soronként
        for item in list_file:
            f.write("%s\n" % item)
        f.close()
    
    
    pwd_list_data = pwd_full_data+'list_file.txt'
    #load_diamides.DiamidesDb DATAbase
    
    DATAbase = load_diamides.DiamidesDb(pwd_list_data)
    #(list_file)