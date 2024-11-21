import argparse
import pickle
import pysam

def arg_parse():
    parser = argparse.ArgumentParser(description = "add CB barcode tag and UB UMI tag to the bam file.")
    parser.add_argument("-technique", type = str, default=None, help = "Options: sptt, stereo, dbit, decoder, slide.")
    parser.add_argument("-in_bam", type = str, default = None, help = "Please input the absolute path of the input bam.")
    parser.add_argument("-out_bam", type = str, default = None, help = "Please input the ouput bam path.")
    return parser.parse_args()

class preprocess:
    def __init__(self, args):
        self.technique = args.technique
        self.in_bam = args.in_bam
        self.out_bam = args.out_bam

        if self.technique == "sptt":
            self.barcode_len = 12
            self.umi_len = 10
        if self.technique == "stereo":
            self.barcode_len = 25
            self.umi_len = 10
        if self.technique == "decoder":
            self.barcode_len = 16
            self.umi_len = 12
        if self.technique == "dbit":
            self.barcode_len = 16
            self.umi_len = 10
        if self.technique == "slide":
            self.barcode_len = 15
            self.umi_len = 8
    
    def process(self):        
        '''
        Dbit-seq & Stereo-seq: read name
        Decoder-seq & SPTT: BC tag
        '''

        bf = pysam.AlignmentFile(self.in_bam, "rb", check_sq = False)
        bf_head_dict = dict(bf.header)

        with pysam.AlignmentFile(self.out_bam, "wb", header = bf_head_dict) as outf:
            for r in bf:
                if self.technique == "dbit" or self.technique == "stereo":
                    to_grep = r.qname.split("_")
                    barcode = to_grep[1]
                    umi = to_grep[2]

                if self.technique == "sptt" or self.technique == "slide":
                    barcode = r.get_tag("XC")
                    umi = r.get_tag("XM")
        
                if self.technique == "decoder":
                    to_grep = r.get_tag("BC")
                    barcode = to_grep[0:self.barcode_len]
                    umi = to_grep[self.barcode_len: self.barcode_len + self.umi_len]
                
                r.set_tag("CB", barcode)
                r.set_tag("UB", umi)
                outf.write(r)
        
        outf.close()
        return
    
def main():
    args = arg_parse()
    to_preprocess = preprocess(args)
    to_preprocess.process()

if __name__ == "__main__":
    main()
