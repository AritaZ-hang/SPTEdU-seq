import argparse
import pickle
import pysam
from tqdm import tqdm

def load_bcdict(dict_path):
    if dict_path is not None:
        dict = pickle.load(open(dict_path, "rb"))
        return dict
    else:
        return None

def args_parse():
    parser = argparse.ArgumentParser(description="correct beads barcodes at the bam level for SPTT-seq.", add_help = True)
    parser.add_argument("-barcode_dict", type = str, default = None, help = "Please input the absolute path of the barcode mapping dict.")
    parser.add_argument("-in_bam", type = str, default = None, help = "Please input the absolute path of the input bam.")
    parser.add_argument("-out_bam", type = str, default = None, help = "Please input the output bam path.")
    parser.add_argument("-barcode_tag", type = str, default = None, help = "Please input the barcode tag in bam file.")
    return parser.parse_args()

class Correct:
    def __init__(self, args):
        self.barcode_dict = load_bcdict(args.barcode_dict)
        self.in_bam = args.in_bam
        self.out_bam = args.out_bam
        self.barcode_tag = args.barcode_tag

    def mapping_list_preprocessing(self):
        if self.barcode_dict is not None:
            all_mapped_bead = list(self.barcode_dict.keys())
        else:
            all_mapped_bead = None
        return all_mapped_bead

    def correct_barcodes(self):
        in_bam = self.in_bam
        out_bam = self.out_bam
        barcode_tag = self.barcode_tag

        all_mapped_bead = set(self.mapping_list_preprocessing())
        # read in bamfile and get header
        bf = pysam.AlignmentFile(in_bam, "rb", check_sq=False)
        bf_head_dict = dict(bf.header)

        # count total and filtered lines
        total_line = 0
        filtered_line = 0

        with pysam.AlignmentFile(out_bam, "wb", header = bf_head_dict) as outf:
            for r in bf:
                total_line += 1
                barcode = r.get_tag(barcode_tag)
                if barcode in all_mapped_bead:
                    new_barcode = self.barcode_dict[barcode]
                    r.set_tag(barcode_tag, new_barcode)
                    outf.write(r)

        outf.close()

        # message filter rate
        print("total line: %f, filtered line: %f, filter rate: %f"%(total_line, filtered_line, float(filtered_line) / float(total_line)))
        com_message = '''~~~~~~~~~~~~~~~correct barcodes done~~~~~~~~~~~~~~~~~~'''
        print(com_message)

def main():
    args = args_parse()
    to_correct = Correct(args)
    to_correct.correct_barcodes()

if __name__ == "__main__":
    main()