#!/bin/bash

cd `pwd`
tmpdir=tmp
mkdir -p $tmpdir
sample_name=$(basename `pwd`)
bbmap_root=/your/bbmap/root/
scripts_dir=/your/scripts/root/

workdir=`pwd`
echo $workdir

linker3_seq=TCTTCAGCGTTCCCGAGAGCAGATGCA

in_fq_1=H_R1.final.fastq.gz
in_fq_2=H_R2.final.fastq.gz

##### filter linker #####
${bbmap_root}/bbduk2.sh in=$in_fq_2 in2=$in_fq_1 outm=H_R2_withlinker.fastq outm2=H_R1_withlinker.fastq out=H_R2_withoutlinker.fastq out2=H_R1_withoutlinker.fastq fliteral=$linker3_seq k=27 skipr2=t hdist=2

withlinkerDir=WithLinker
mkdir -p $withlinkerDir

withoutlinkerDir=WithoutLinker
mkdir -p $withoutlinkerDir

gzip H_R1_withlinker.fastq H_R2_withlinker.fastq
gzip H_R1_withoutlinker.fastq H_R2_withoutlinker.fastq

mv H_R1_withlinker.fastq.gz H_R2_withlinker.fastq.gz $withlinkerDir
mv H_R1_withoutlinker.fastq.gz H_R2_withoutlinker.fastq.gz $withoutlinkerDir

cd ${workdir}/${withlinkerDir}
cp ${scripts_dir}/SPTT_upstream_withlinker_UTI.sh .

cd ${workdir}/${withoutlinkerDir}
cp ${scripts_dir}/SPTT_upstream_withoutlinker.sh .
