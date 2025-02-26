#!/bin/bash

cd `pwd`

rmats="/path/to/rmats_turbo/rmats.py"
read_len=150
threads=10
gtf="/path/to/genome_reference/Mus_musculus.GRCm38.88.gtf"
rmats_out_dir="vs10X_cstat0.05"

mkdir -p $rmats_out_dir

cd $rmats_out_dir

INPUT_SPTT=/path/to/SPTT_MOB/SPTT_MOB.bam
INPUT_VISIUM=/path/to/Visium_MOB/Visium_MOB.bam

od=${rmats_out_dir}/SPTT_10X
mkdir -p $od

echo $INPUT_SPTT > $od/b1.txt
echo $INPUT_VISIUM > $od/b2.txt

python $rmats --gtf $gtf --b1 $od/b1.txt --b2 $od/b2.txt --od $od --tmp $od/tmp --nthread $threads -t single --allow-clipping --readLength $read_len --variable-read-length --cstat 0.05
