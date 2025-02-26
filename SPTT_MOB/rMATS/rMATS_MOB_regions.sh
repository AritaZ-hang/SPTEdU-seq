#!/bin/bash
rmats="/path/to/rmats/rmats.py"
split_py="/path/to/scripts/bam_by_clu.py"
in_bam="/path/to/SPTT_MOB/SPTT_MOB.bam"
read_len=150
pixel_to_clu_ID_tsv="/path/to/data/pixel_to_clu.txt"
clu_ID_to_region_names="/path/to/data/clu_ID_to_region_names.txt"
N_regions=5
bam_by_clu_dir="bam_by_clu_name"
threads=1
gtf="/path/to/genome_reference/Mus_musculus.GRCm38.88.gtf"
all_clu_bams=($(ls $bam_by_clu_dir/*.bam))
rmats_out_dir="rmats_out_dir_cstat0.05"

for((clu_i=0; clu<$N_regions-1; ++clu_i))
do
    for((clu_j=clu_i+1;clu_j<$N_regions;++clu_j))
    do
        clu_i_bam=${all_clu_bams[$clu_i]}
        clu_j_bam=${all_clu_bams[$clu_j]}
        basename_i=$(basename ${clu_i_bam%.bam})
        basename_j=$(basename ${clu_j_bam%.bam})
        od=${rmats_out_dir}/${basename_i}_${basename_j}
        mkdir -p $od
        echo "$clu_i_bam" > $od/b1.txt
        echo "$clu_j_bam" > $od/b2.txt
        python $rmats --gtf $gtf --b1 $od/b1.txt --b2 $od/b2.txt --od $od --tmp $od/tmp --nthread $threads -t single --readLength $read_len --allow-clipping --variable-read-length --cstat 0.05

    done
done
