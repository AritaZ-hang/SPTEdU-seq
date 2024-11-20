import argparse
import pysam

'''
BBV withlinker & without linker
'''
def str2bool(v):
    if v.lower() in ("false", "f", "no", "0"):
        return False
    elif v.lower() in ("true", "t", "yes", "1"):
        return True
    else:
        raise "Please input valid options!"

def args_parse():
    parser = argparse.ArgumentParser(description = "A tiny script to cope with the barcode and UMI tag.")
    parser.add_argument("-in_bam", type = str, default = None, help = "The filename of input bam.")
    parser.add_argument("-out_bam", type = str, default = None, help = "The filename of output bam.")
    parser.add_argument("-linker", type = str, default = None, help = "Linker or not.")
    return parser.parse_args()

def add_barcode_and_umi_tags(in_bam:str, out_bam:str, linker_or_not:bool):
    # read in bam file and get header
    bf = pysam.AlignmentFile(in_bam, "rb", check_sq = False)
    bf_head_dict = dict(bf.header)

    # add barcode and umi tags
    with pysam.AlignmentFile(out_bam, "wb", header = bf_head_dict) as outf:
        for r in bf:
            if linker_or_not == True:
                barcode = r.get_tag("BC")
                umi = r.get_tag("UM")
            else:
                barcode = r.get_tag("BC")[0:12]
                umi = r.get_tag("BC")[12:22]
            r.set_tag("XC", barcode)
            r.set_tag("XM", umi)
            outf.write(r)
    outf.close()
    return

def main():
    args = args_parse()
    add_barcode_and_umi_tags(in_bam = args.in_bam, 
                             out_bam = args.out_bam, 
                             linker_or_not=str2bool(args.linker))

if __name__ == "__main__":
    main()