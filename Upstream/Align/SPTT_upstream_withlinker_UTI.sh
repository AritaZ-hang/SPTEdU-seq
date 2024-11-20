#!/bin/bash

cd `pwd`
tmpDir=tmp
mkdir -p $tmpDir
sample_name=$(basename `pwd`)
workdir=`pwd`
scripts_dir=/your/scripts/root/


summariseScript=${scripts_dir}/SPTT_upstream_barcode_UTI_summarise.py
linker_seq=TCTTCAGCGTTCCCGAGAGCAGATGCA

###########################
###### filter linker ######
###########################

trimDir=Trim
mkdir -p $trimDir

cutadapt \
--minimum-length 15 \
-G $linker_seq \
--pair-filter=any \
-o H_R1_filtered.fastq.gz \
-p H_R2_filtered.fastq.gz \
--cores 12 H_R1_withlinker.fastq.gz H_R2_withlinker.fastq.gz > ${trimDir}/cutadapt_linker_report.txt

############################
###### UTI Summarise #######
############################
python -u $summariseScript -in_fq_1 H_R1_filtered.fastq.gz -in_fq_2 H_R2_filtered.fastq.gz -workdir $workdir