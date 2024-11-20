#!/bin/bash

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)
dropseq_root=/your/dropseq/root/
bbmap_root=/your/bbmap/root/
STAR_root=/your/star/root/
script_dir=/your/scripts/root/
genome_dir=/your/genome/dir/

addTagScript=${script_dir}/techs_addTagScript.py

technique=dbit

##### extract beads barcode & UMI from R1 #####
umi_tools extract --stdin H_R1.fastq.gz \
    --extract-method=regex \
    --bc-pattern="(?P<cell_1>.{16})(?P<umi_1>.{10}).*" \
    --stdout H_R1_extracted.fq.gz --read2-in H_R2.fastq.gz --read2-out H_R2_extracted.fastq.gz \
   -L extract.log

$STAR_root \
--runThreadN 12  --runMode alignReads \
 --genomeDir ${genome_dir} \
 --readFilesIn H_R2_extracted.fastq.gz \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 78628592870 \
 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --alignIntronMax 1 

samtools view -b -q 250 Aligned.sortedByCoord.out.bam > dbit.bam
samtools sort -@ 12 -o dbit.sort.bam dbit.bam

python -u $addTagScript -technique $technique -in_bam $in_bam -out_bam $out_bam

${dropseq_root}/TagReadWithGeneFunction I=$out_bam \
O=star_gene_exon_tagged.bam \
ANNOTATIONS_FILE=${genome_dir}/Mus_musculus.GRCm38.88.gtf && rm $out_bam

${dropseq_root}/DigitalExpression -m 8g \
I=star_gene_exon_tagged.bam \
CELL_BARCODE_TAG=CB \
MOLECULAR_BARCODE_TAG=UB \
O=dbit_dge.txt.gz \
SUMMARY=dbit_dge.summary.txt \
NUM_CORE_BARCODES=10000 \
LOCUS_FUNCTION_LIST=INTRONIC \
STRAND_STRATEGY=BOTH \
TMP_DIR=.
