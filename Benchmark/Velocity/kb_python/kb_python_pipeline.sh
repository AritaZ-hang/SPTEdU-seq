#!/bin/bash


dataDir=/your/data/dir/ 
RefDir=/your/ref/dir/
indexDir=${RefDir}/kb_index # kb_index

mkdir -p $indexDir
cd $indexDir

# build ref
kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno ${RefDir}/Mus_musculus.GRCm38.88.fasta ${RefDir}/Mus_musculus.GRCm38.88.gtf

# count unspliced and spliced reads
# R1: 12 spatial barcodes & 10 UMI, R2: reads

workdir=/your/work/dir/
cd $workdir

kb count --h5ad -i ${indexDir}/index.idx -g ${indexDir}/t2g.txt -x 0,0,12:0,12,22:1,0,0 -o SPTT_MOB -c1 ${indexDir}/cdna_t2c.txt -c2 ${indexDir}/intron_t2c.txt --workflow lamanno --filter bustools -t 2 ${dataDir}/H_R1.final.fastq.gz ${dataDir}/H_R2.final.fastq.gz
