import argparse
import pickle
import pysam

def load_bcdict(dict_path):
    if dict_path is not None:
        dict = pickle.load(open(dict_path, "rb"))
        return dict
    else:
        return None
    
def args_parse():
    parser = argparse.ArgumentParser(description="correct bead barcodes for decoder-seq.", add_help=True)
    parser.add_argument("-barcode_xdict", type = str, default=None, help="Barcode X correction dictionary.")
    parser.add_argument("-barcode_ydict", type = str, default = None, help = "Barcode Y correction dictionary.")
    parser.add_argument("-input_bam", type = str, default = None, help = "The bam to correct barcodes.")
    parser.add_argument("-output_bam", type = str, default = None, help = "The output bam.")

    return parser.parse_args()

def correct_barcodes(barcode_xdict, barcode_ydict, input_bam, output_bam):

    if barcode_xdict is None or barcode_ydict is None:
        raise "Please input correct path of barcode correction dictionary."

    barcode_x_len = 8
    barcode_y_len = 8
    umi_len = 12

    bf = pysam.AlignmentFile(input_bam, "rb", check_sq = False)
    bf_head_dict = dict(bf.header)

    # count total and filtered lines
    total_line = 0
    filtered_line = 0

    with pysam.AlignmentFile(output_bam, "wb", header = bf_head_dict) as outf:
        for r in bf:
            total_line += 1
            tags = r.tags
            tags_only = []
            for i in range(len(tags)):
                tags_only.append(tags[i][0])
            if "BC" in tags_only:
                barcode = r.get_tag("BC")
                barcode_x = barcode[0:barcode_x_len]
                barcode_y = barcode[barcode_x_len:barcode_x_len + barcode_y_len]
                umi = barcode[barcode_x_len + barcode_y_len: barcode_x_len + barcode_y_len + umi_len]

                if barcode_x in barcode_xdict and barcode_y in barcode_ydict:
                    filtered_line += 1
                    barcode_corrected = barcode_xdict[barcode_x] + barcode_ydict[barcode_y]

                    r.set_tag("XC", barcode_corrected)
                    r.set_tag("XM", umi)
                    outf.write(r)

        outf.close()
    
    print("total line: %f, filtered line: %f, filter rate: %f"%(total_line, filtered_line, float(filtered_line) / float(total_line)))
    com_message = '''~~~~~~~~~~~~~~~correct barcodes done~~~~~~~~~~~~~~~~~~'''
    print(com_message)

if __name__ == "__main__":
    args = args_parse()
    correct_barcodes(barcode_xdict = load_bcdict(args.barcode_xdict), barcode_ydict = load_bcdict(args.barcode_ydict), input_bam = args.input_bam, output_bam = args.output_bam)

