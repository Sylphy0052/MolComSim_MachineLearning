import os
import pickle
from dat import Dat
import glob
from natsort import natsorted
from makeText import createCSV
from analyzer import Analyzer

def create_pickle_file(file_list):
    data_dict = {}
    if os.path.isfile('pickle_file_data.txt'):
        with open('pickle_file_data.txt', 'rb') as f:
            data_dict = pickle.load(f)
    else:
        count = 1
        for file_name in file_list:
            print("{}/{} - {} file reading...".format(count, len(file_list), file_name.split("/")[2]))
            data_dict[file_name] = Dat(file_name)
            count += 1
        with open('pickle_file_data.txt', 'wb') as f:
            pickle.dump(data_dict, f)

    return data_dict

def main():
    file_list = natsorted(glob.glob("./dat/*.dat"))
    data_dict = create_pickle_file(file_list)
    column = createCSV(data_dict, file_list)
    analyzer = Analyzer(column)

if __name__ == '__main__':
    main()
