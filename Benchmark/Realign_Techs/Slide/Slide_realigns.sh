#!/bin/bash


workdir=`pwd`

dropseq_root=/your/dropseq/root/
genome_dir=/your/genome/dir/
gtf_file=${genome_dir}/Mus_musculus.GRCm38.88.gtf
star_root=/your/star/root/
script_dir=/your/scripts/dir/

# add barcode and UMI to read names

IN_BAM=Puck_200127_15.bam
OUT_BAM=Slide_raw.bam

python -u ${script_dir}/Slide_AddBarcodeToReadName.py $IN_BAM $OUT_BAM XC XM

# convert bam to fastq
samtools bam2fq $OUT_BAM > slide.fq
gzip slide.fq 

# Star realign
${star_root}/STAR --genomeDir ${genome_dir} --runMode alignReads --readFilesIn slide.fq.gz --readFilesCommand zcat --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --alignIntronMax 1 --limitOutSJcollapsed 5000000 --outSAMtype BAM SortedByCoordinate --runThreadN 12 --limitBAMsortRAM 41143265264

samtools view -b -q 250 Aligned.sortedByCoord.out.bam > star.bam && rm Aligned.sortedByCoord.out.bam

# add cell barcode & UMI information
python -u ${script_dir}/Slide_tagXCM.py star.bam mouse.bam
${dropseq_root}/TagReadWithGeneFunction I=mouse.bam O=star_gene_exon_tagged.bam ANNOTATIONS_FILE=$gtf_file

# sort the bam
samtools sort -o star_gene_exon_tagged.sort.bam star_gene_exon_tagged.bam -@ 16
samtools index star_gene_exon_tagged.sort.bam -@ 16

# DGE
${dropseq_root}/DigitalExpression -m 8g \
I=star_gene_exon_tagged.sort.bam \
CELL_BARCODE_TAG=XC \
MOLECULAR_BARCODE_TAG=XM \
O=slide_dge.txt.gz \
SUMMARY=slide_dge.summary.txt \
NUM_CORE_BARCODES=10000 \
LOCUS_FUNCTION_LIST=INTRONIC \
STRAND_STRATEGY=BOTH \
TMP_DIR=.

# all reads
${dropseq_root}/BamTagHistogram I=star_gene_exon_tagged.sort.bam O=out_cell_readcounts.txt TAG=XC 
