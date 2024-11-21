import os
import pickle
import itertools
from Levenshtein import distance
import numpy as np
from collections import defaultdict
import pandas as pd

def combination_seqs(input_char, output_length):
    all_combine = ["".join(seq) for seq in itertools.product(input_char, repeat = output_length)]
    return all_combine

def preprocess_bc_list(bc_list):
    # Preprocess the list of barcodes into a dictionary for fast lookup
    bc_dict = defaultdict(list)
    for bc in bc_list:
        bc_dict[len(bc)].append(bc)
    return bc_dict

def find_bc(seq, bc_dict, mismatch_rate=1):
    result = "nobc"
    seq_len = len(seq)
    for i in range(seq_len - mismatch_rate, seq_len + mismatch_rate + 1):
        for bc in bc_dict[i]:
            mismatch = distance(bc, seq)
            if mismatch <= mismatch_rate:
                result = bc
                return result  # Exit early if a match is found
    return result

def generate_bc_dictionary(all_sequences, bc_dict, bc_raw):
    all_seq_match = [find_bc(seq, bc_dict) for seq in all_sequences]
    bc_process = pd.DataFrame({"bc_process": all_seq_match})
    bc_use = pd.merge(bc_process, bc_raw, on="bc_process", how="left")
    all_seq_match = list(bc_use["bc_raw"])
    all_seq_dic = dict(zip(all_sequences, all_seq_match))
    return all_seq_dic

# X & Y are 8nt barcodes.

all_seq_8bp = combination_seqs(["A", "C", "G", "T", "N"], 8)

bc_raw_1_file = "./lib_barcode_X.txt"
bc_raw_2_file = "./lib_barcode_Y.txt"

bc_raw_1 = open(bc_raw_1_file)
bc_raw_2 = open(bc_raw_2_file)

bc1_bc = list(map(lambda x: x.strip(), bc_raw_1.readlines()))
bc2_bc =  list(map(lambda x: x.strip(), bc_raw_2.readlines()))

print(bc1_bc[0:10])
print(len(bc1_bc))
print(bc2_bc[0:10])
print(len(bc2_bc))

bc_raw1 = pd.read_csv("./lib_barcode_X_df.txt", sep = " ")
bc_raw2 = pd.read_csv("./lib_barcode_Y_df.txt", sep = " ")

preprocessed_bc_list1 = preprocess_bc_list(bc1_bc)
print("Generating the dictionary for bc1 barcodes...")
all_seq_dic1 = generate_bc_dictionary(all_seq_8bp, preprocessed_bc_list1, bc_raw1)

preprocessed_bc_list2 = preprocess_bc_list(bc2_bc)
print("Generating the dictionary for bc2 barcodes...")
all_seq_dic2 = generate_bc_dictionary(all_seq_8bp, preprocessed_bc_list2, bc_raw2)

##### save the dictionary into a file
pickle_out1 = open("./lib_barcodeX.pickle", "wb")
pickle.dump(all_seq_dic1, pickle_out1)
pickle_out1.close()
print("Generating the dictionary for ligation barcodes...")

pickle_in1 = open("./lib_barcodeX.pickle", "rb")
bc_dic1 = pickle.load(pickle_in1)
bc_filter1 =  {k: v for k, v in bc_dic1.items() if v != "nobc"}
pickle_out1 = open("./lib_barcodeX.pickle2", "wb")
pickle.dump(bc_filter1, pickle_out1, 2)
pickle_in1.close()
pickle_out1.close()

pickle_out2 = open("./lib_barcodeY.pickle", "wb")
pickle.dump(all_seq_dic2, pickle_out2)
pickle_out2.close()
print("Generating the dictionary for ligation barcodes...")

pickle_in2 = open("./lib_barcodeY.pickle", "rb")
bc_dic2 = pickle.load(pickle_in2)
bc_filter2 =  {k: v for k, v in bc_dic2.items() if v != "nobc"}
pickle_out2 = open("./lib_barcodeY.pickle2", "wb")
pickle.dump(bc_filter2, pickle_out2, 2)
pickle_in2.close()
pickle_out2.close()
