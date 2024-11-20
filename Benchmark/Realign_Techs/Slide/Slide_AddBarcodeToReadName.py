import sys
import pysam

def addqname(in_bam, out_bam, barcode_tag, umi_tag):
    bf = pysam.AlignmentFile(in_bam, "rb", check_sq = False)
    bf_head_dict = dict(bf.header)
    with pysam.AlignmentFile(out_bam, 'wb', header=bf_head_dict) as outf:
        for r in bf:
            q_name = r.qname
            cell_barcode = r.get_tag(barcode_tag)
            umi = r.get_tag(umi_tag)
            new_q_name = q_name + "_" + cell_barcode + "_" + umi
            r.qname = new_q_name
            outf.write(r)
        outf.close()
    bf.close()
    return

if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    barcode_tag = sys.argv[3]
    umi_tag = sys.argv[4]
    addqname(in_bam, out_bam, barcode_tag, umi_tag)