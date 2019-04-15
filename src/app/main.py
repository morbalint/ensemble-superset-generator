import os
import sys
sys.path.append(os.getcwd()) # this is fuking ridiculus
import src.app.rd as rd

if __name__ == '__main__':
    print("hello world!")
    rd.json_test()
    rd.TDRD_test()
    rd.NDRD_test()
    print("hello world!")