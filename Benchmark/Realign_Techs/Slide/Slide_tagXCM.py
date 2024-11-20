import sys
import pickle
import pysam


def addtag(in_bam, out_bam):
	# read in bamfile and get header
	bf = pysam.AlignmentFile(in_bam, 'rb', check_sq=False)
	bf_head_dict = dict(bf.header)


	with pysam.AlignmentFile(out_bam, "wb", header=bf_head_dict) as outf:
		for r in bf:  
			barcode1 = r.qname.split("_")[1]
			r.set_tag('XC', barcode1)
			barcode2 = r.qname.split("_")[2]
			r.set_tag('XM', barcode2)
			outf.write(r)
		outf.close()


if __name__ == "__main__":
	in_bam = sys.argv[1]
	out_bam = sys.argv[2]
	addtag(in_bam, out_bam)
