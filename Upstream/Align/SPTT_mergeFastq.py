import argparse
import pysam
from tqdm import tqdm
from xopen import xopen

def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'

def args_parse():
    parser = argparse.ArgumentParser(description="Merge BC1 & BC2 & BC3 Fastqs to make a integrate R1.")
    parser.add_argument("-bc1", type = str, default = None, help = "The filename of input bc1.")
    parser.add_argument("-bc2", type = str, default = None, help = "The filename of input bc2.")
    parser.add_argument("-bc3", type = str, default=None, help = "The filename of input bc3.")
    parser.add_argument("-r2", type = str, default=None, help = "The filename of input R2.")
    parser.add_argument("-out_fq_1", type = str, default = None, help = "The output file name of R1. It should have a .gz suffix.")
    parser.add_argument("-out_fq_2", type = str, default = None, help = "The output file name of R2. It should have a .gz suffix.")

    return parser.parse_args()

class Prepare:
    def __init__(self, args):
        self.bc1 = args.bc1
        self.bc2 = args.bc2
        self.bc3 = args.bc3
        self.r2 = args.r2
        self.out_fq_1 = args.out_fq_1
        self.out_fq_2 = args.out_fq_2
    def process_barcode(self):
        '''
        merge bc1 & bc2 & bc3
        '''
        out_h = xopen(self.out_fq_1, "w")
        out_h2 = xopen(self.out_fq_2, "w")
        with pysam.FastxFile(self.bc1) as bc1_file, pysam.FastxFile(self.bc2) as bc2_file, pysam.FastxFile(self.bc3) as bc3_file, pysam.FastxFile(self.r2) as r2_file:
            for read1, read2, read3, read4 in tqdm(zip(bc1_file, bc2_file, bc3_file, r2_file)):
                fastq_read1 = {"name":read1.name, 
                               "sequence": read1.sequence, 
                               "comment":read1.comment if read1.comment else "", 
                               "quality":read1.quality if read1.quality else ""}
                fastq_read2 = {"name":read2.name, 
                               "sequence": read2.sequence, 
                               "comment":read2.comment if read2.comment else "", 
                               "quality": read2.quality if read2.quality else ""}
                fastq_read3 = {"name":read3.name, 
                               "sequence": read3.sequence, 
                               "comment": read3.comment if read3.comment else "", 
                               "quality": read3.quality if read3.quality else ""}
                fastq_read4 = {"name":read4.name, 
                               "sequence": read4.sequence, 
                               "comment": read4.comment if read4.comment else "", 
                               "quality": read4.quality if read4.quality else ""}
                
                if len(fastq_read1["sequence"]) >=4 and len(fastq_read2["sequence"]) >=4 and len(fastq_read3["sequence"]) >=4:
                    bc1 = fastq_read1["sequence"][0:4]
                    bc2 = fastq_read2["sequence"][0:4]
                    bc3 = fastq_read3["sequence"][0:4]
                    if len(fastq_read3["sequence"]) >= 14:
                        umi = fastq_read3["sequence"][4:14]
                    else:
                        umi = fastq_read3["sequence"][4:len(fastq_read3["sequence"])] + "T"*(14-len(fastq_read3["sequence"]))
                        
                    if bc1 and bc2 and bc3 and umi:
                        header = fastq_read1["name"]
                        seq = bc1 + bc2 + bc3 + umi
                        qual = fastq_read1["quality"][0:4] + fastq_read2["quality"][0:4] + fastq_read3["quality"][0:4] + "I"*len(umi)

                        header2=fastq_read4["name"]
                        seq2 = fastq_read4["sequence"]
                        qual2 = fastq_read4["quality"]

                        out_h.write(fastq_line(header, seq, qual))
                        out_h2.write(fastq_line(header2, seq2, qual2))
        
        out_h.close()
        out_h2.close()

def main():
    args = args_parse()
    preprocessing = Prepare(args)
    preprocessing.process_barcode()

if __name__ == "__main__":
    main()