import sys
import os
import pysam as ps
import mappy as mp
from collections import defaultdict as dd


def get_spot_to_clu(spot_to_clu_fn, clu_to_name_tsv):
    spot_to_clu = dd(lambda: 'NotClustered')
    clu_to_name = dict()
    with open(spot_to_clu_fn) as spot_fp, open(clu_to_name_tsv) as clu_fp:
        for line in clu_fp:
            ele = line.strip().rsplit()
            if len(ele) < 2:
                continue
            clu, name = ele[0], ele[1]
            clu_to_name[clu] = name
        for line in spot_fp:
            ele = line.strip().rsplit()
            if len(ele) < 2:
                continue
            spot, clu = ele[0], ele[1]
            if clu not in clu_to_name:
                name = 'NoClusterName'
                # sys.stderr.write('Cluster not in clu_to_name.tsv: {}\n'.format(clu))
                # continue
            else:
                name = clu_to_name[clu]
            spot_to_clu[spot] = name
    return spot_to_clu

def get_bam_fps(bam_fn, spot_to_clu, split_bam_dir, collect_non_clu = False):
    bam_fps = dd(lambda:None)
    all_clus = set(spot_to_clu.values())
    for clu in all_clus:
        if not collect_non_clu:
            split_bam_fn = split_bam_dir + "/clu" + clu + ".bam"
            split_bam_fp = ps.AlignmentFile(split_bam_fn, "wb", template=ps.AlignmentFile(bam_fn))
        else:
            split_bam_fn = split_bam_dir + "/non_clu" + clu + ".bam"
            split_bam_fp = ps.AlignmentFile(split_bam_fn, "wb", template=ps.AlignmentFile(bam_fn))
        bam_fps[clu] = split_bam_fp
    return bam_fps

def get_spot_from_name(qname):
    return qname

def split_bam(bam_fn, bam_fps, spot_to_clu, bc_in_name=True, collect_non_clu=False):
    if not collect_non_clu:
        with ps.AlignmentFile(bam_fn) as bam_fp:
            for read in bam_fp:
                if read.is_unmapped:
                    continue
                if bc_in_name:
                    spot = get_spot_from_name(read.qname)
                    if spot not in spot_to_clu:
                        sys.stderr.write("Unclustered\t{}\n".format(read.qname))
                        continue
                    clu = spot_to_clu[spot]
                bam_fps[clu].write(read)
        return
    
if __name__ == "__main__":
    bam_fn, spot_to_clu_fn, clu_to_name_tsv, split_bam_dir, collect_non_bam = sys.argv[1:]
    bc_in_name=True
    collect_non_bam=False
    spot_to_clu = get_spot_to_clu(spot_to_clu_fn, clu_to_name_tsv)
    if not os.path.exists(split_bam_dir):
        os.mkdir(split_bam_dir)
    bam_fps = get_bam_fps(bam_fn, spot_to_clu, split_bam_dir, collect_non_bam)
    split_bam(bam_fn, bam_fps, spot_to_clu, bc_in_name, collect_non_bam)
