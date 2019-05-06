import os
import sys
sys.path.append(os.getcwd()) # this is fuking ridiculus
import src.app.load_diamides as ld
import src.app.rd as rd
# import rd,load_diamides

def get_pwd_full_data():
    pwd = os.getcwd()
    pwd_home = pwd[0:pwd.find("src")+1]                     #Az a mappa ami tartalmazza a project fileokat
    pwd_data = "samples/diamides"+os.sep                    #Az a mappa ami tartalmazza a pdb fileokat.
    pwd_full_data = os.path.join(pwd_home,pwd_data)         #A pdb fileok teljes elérési útvonala, igyekeztem ezt dinamikussá tenni
    return pwd_full_data                                    #hogy mindkettőnk gépén lefusson

def create_list_file(pwd_full_data):
    list_file = os.listdir(pwd_full_data)                   # Egy listába írom a pdb fileneveket
    with open( pwd_full_data+'list_file.txt', 'w') as f:    # A lista tartalmát egy fileba írom soronként
        for item in list_file:
            if item != 'list_file.txt':
                f.write("%s\n" % item)
        f.close()

if __name__ == '__main__':
    # we must always write our program between hello world printing to remember where we started from.
    print("hello world!")

    pwd_full_data = get_pwd_full_data()
    create_list_file(pwd_full_data)
    pwd_list_data = pwd_full_data+'list_file.txt'
    #load_diamides.DiamidesDb DATAbase
    DATAbase = ld.DiamidesDb(pwd_list_data)
    #(list_file)
    # we must always write our program between hello world printing to remember where we started from.
    print("hello world!")