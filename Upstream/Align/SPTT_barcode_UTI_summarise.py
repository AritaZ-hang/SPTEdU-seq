import argparse
import pysam
import pandas as pd
from tqdm import tqdm
from xopen import xopen
import pickle

def args_parse():
    parser = argparse.ArgumentParser(description="Summarise UTI tag. Beads Barcodes is the 1-12nt in R1, UMI tag is the 1-10nt in R2.")
    parser.add_argument("-in_fq_1", type = str, default=None, help = "The filename of input fastq1.")
    parser.add_argument("-in_fq_2", type = str, default=None, help = "The filename of input fastq2.")
    parser.add_argument("-workdir", type = str, default = None, help = "The working directory.")
    return parser.parse_args()

def writeDictToCSV(dict_to_summarize):
    target_csv = pd.DataFrame.from_dict(dict_to_summarize, orient = "index")
    target_csv = target_csv.fillna(0)
    return target_csv

def uti_summarise(in_fq_1: str, in_fq_2: str, workdir: str):
    beads_umi_dict = dict()
    with pysam.FastxFile(in_fq_1) as fastq_file1, pysam.FastxFile(in_fq_2) as fastq_file2:
        for read1, read2 in tqdm(zip(fastq_file1, fastq_file2)):
            fastq_read1 = {"name":read1.name, 
                           "sequence": read1.sequence, 
                           "comment": read1.comment if read1.comment else "", 
                           "quality": read1.quality if read1.quality else ""}
            
            fastq_read2 = {"name":read2.name, 
                           "sequence": read2.sequence, 
                           "comment": read2.comment if read2.comment else "", 
                           "quality": read2.quality if read2.quality else ""}
            barcodes = fastq_read1["sequence"][0:12]
            umi = fastq_read2["sequence"][0:10]

            if barcodes not in beads_umi_dict.keys():
                beads_umi_dict[barcodes] = {}
            if umi not in beads_umi_dict[barcodes].keys():
                beads_umi_dict[barcodes][umi] = 1
            else:
                beads_umi_dict[barcodes][umi] += 1
            
    with open(workdir + "/beads_umi_dict.pkl", "wb") as f:
        pickle.dump(beads_umi_dict, f)
    
    utis = {}
    for group in beads_umi_dict.keys():
        if group not in utis.keys():
            utis[group] = len(beads_umi_dict[group].keys())
    
    with open(workdir + "utis.pkl", "wb") as f:
        pickle.dump(utis, f)

    utis_df = writeDictToCSV(utis)
    utis_df.columns = ["reads"]
    utis_df.to_csv(workdir + '/utis.csv')

    return

def main():
    args = args_parse()
    uti_summarise(in_fq_1 = args.in_fq_1, 
                  in_fq_2 = args.in_fq_2, 
                  workdir = args.workdir)
    return

if __name__ == "__main__":
    main()