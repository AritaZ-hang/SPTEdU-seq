#!/bin/bash

cd `pwd`
tmpdir=tmp
mkdir -p $tmpdir
sample_name=$(basename `pwd`)
bbmap_root=/your/bbmap/root/
scripts_dir=/your/scripts/root/

workdir=`pwd`
echo $workdir

linker1_seq=CGACTCACTACACGT
linker2_seq=TCGCTGACACGATCG

#################################################
####### Filter correctly structured reads #######
#################################################

### filter R1 linker ###
${bbmap_root}/bbduk2.sh in=H_R1.fastq.gz in2=H_R2.fastq.gz outm=H_R1.linker1.fastq outm2=H_R2.linker1.fastq fliteral=$linker1_seq k=15 hdist=2

${bbmap_root}/bbduk2.sh in=H_R1.linker1.fastq in2=H_R2.linker1.fastq outm=H_R1.linker2.fastq outm2=H_R2.linker2.fastq fliteral=$linker2_seq k=15 hdist=2 && rm H_R1.linker1.fastq H_R2.linker1.fastq

gzip H_R1.linker2.fastq
gzip H_R2.linker2.fastq

# filter polyT
cutadapt \
--minimum-length 30 \
-a T{10} \
--pair-filter=any \
-o H_R1.use.fastq.gz \
-p H_R2.use.fastq.gz \
--cores 12 H_R1.linker2.fastq.gz H_R2.linker2.fastq.gz > polyT_trim_report.txt && rm H_R1.linker2.fastq.gz H_R2.linker2.fastq.gz

######################################
###### R1 barcode prerpocessing ######
######################################

trimDir=Trim
mkdir -p $trimDir

cutadapt \
--minimum-length 0 \
-a $linker1_seq \
-e 0.2 \
-o H_R1.bc1.fastq.gz \
--cores 12 H_R1.use.fastq.gz > ${trimDir}/cutadapt_bc1_report.txt

cutadapt \
--minimum-length 0 \
-g $linker1_seq \
-e 0.2 \
-o H_R1.bc2.tmp.fastq.gz \
--cores 12 H_R1.use.fastq.gz > ${trimDir}/cutadapt_bc2_tmp_report.txt

cutadapt \
--minimum-length 0 \
-a $linker2_seq \
-e 0.2 \
-o H_R1.bc2.fastq.gz \
--cores 12 H_R1.bc2.tmp.fastq.gz > ${trimDir}/cutadapt_bc2_report.txt 

cutadapt \
--minimum-length 0 \
-g $linker2_seq \
-e 0.2 \
-o H_R1.bc3UMI.fastq.gz \
--cores 12 H_R1.use.fastq.gz > ${trimDir}/cutadapt_bc3UMI_report.txt

##########################################
######## bc1 & bc2 & bc3UMI merge ########
##########################################
script=${scripts_dir}/SPTT_mergeFastq.py

python -u $script -bc1 H_R1.bc1.fastq.gz -bc2 H_R1.bc2.fastq.gz -bc3 H_R1.bc3UMI.fastq.gz -r2 H_R2.use.fastq.gz -out_fq_1 H_R1.final.fastq.gz -out_fq_2 H_R2.final.fastq.gz && rm H_R1.bc1.fastq.gz H_R1.bc2.fastq.gz H_R1.bc2.tmp.fastq.gz H_R1.bc3UMI.fastq.gz H_R1.use.fastq.gz H_R2.use.fastq.gz


cp ${scripts_dir}/SPTT_upstream_clearStruct.sh .
bash SPTT_upstream_clearStruct.sh

