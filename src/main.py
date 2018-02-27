from dat import Dat
import glob
from natsort import natsorted
from makeText import createCSV
from analyzer import Analyzer

def main():
    all_dat = {}
    file_list = glob.glob("./dat/*.dat")
    file_list = natsorted(file_list)
    for file_name in file_list:
        print(file_name)
        all_dat[file_name] = Dat(file_name)
    column = createCSV(all_dat, file_list)
    analyzer = Analyzer(column)


if __name__ == '__main__':
    main()
