import csv
import os

def createCSV(all_data, file_list):
    if os.path.isfile('output.csv'):
        return all_data[file_list[0]].getColumn()

    f = open('output.csv', 'w')
    writer = csv.writer(f, lineterminator='\n')

    for file_name in file_list:
        data = all_data[file_name]
        allInfo = data.getAllInfo()
        csvArray = allInfo
        csvArray.append(-1)
        for step in data.input_data.steps:
            csvArray[-1] = step
            writer.writerow(csvArray)

    f.close()
    return all_data[file_list[0]].getColumn()
