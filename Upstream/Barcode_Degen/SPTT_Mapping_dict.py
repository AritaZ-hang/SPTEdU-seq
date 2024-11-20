import pickle
import csv
import sys
import os


def csv_to_dict(filename):
    data = {}
    with open(filename, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header row
        for row in reader:
            key = row[0]
            value = row[1]
            data[key] = value
    return data

def main():
    mapping_file = sys.argv[1]
    workdir = sys.argv[2]
    mapping_dict = csv_to_dict(mapping_file)

    filename, file_extension = os.path.splitext(mapping_file)

    with open(filename + ".pkl", "wb") as f:
        pickle.dump(mapping_dict, f)
    
if __name__ == "__main__":
    main()