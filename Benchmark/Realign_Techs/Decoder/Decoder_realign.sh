#!/bin/bash

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)
dropseq_root=/your/dropseq/root/
star_root=/your/star/root/
genome_dir=/your/genome/dir/
script_dir=/your/scripts/dir/
decoder_correctBC_script=${script_dir}/Decoder_correctBC.py

############## fastq to sam ############
java -jar ${dropseq_root}/3rdParty/picard/picard.jar FastqToSam F1=H_R1.fastq.gz F2=H_R2.fastq.gz O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name

############## Grep X & Y cell barcodes #############
#### 1-8 represents barcodeX, 39-46 represents barcodeY, and 47-58 represents UMIs.
${dropseq_root}/TagBamWithReadSequenceExtended \ BASE_RANGE=1-8:39-46:47-58 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=true TAG_NAME=BC NUM_BASES_BELOW_QUALITY=4 \
INPUT=H.bam OUTPUT=$tmpdir/H1.bam COMPRESSION_LEVEL=5

#### Filter BAM
${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=$tmpdir/H1.bam OUTPUT=$tmpdir/H2.bam

#### correct barcodeX & barcodeY
python -u $correctBC -barcode_xdict lib_barcodeX.pickle2 -barcode_ydict lib_barcodeY.pickle2 -input_bam $tmpdir/H2.bam -output_bam filtered.bam

java -Xmx100g -jar ${dropseq_root}/3rdParty/picard/picard.jar \
SamToFastq INPUT=filtered.bam FASTQ=R2.fastq

####  Trim
cutadapt -a A{10} -j 0 -O 10 --minimum-length=15 -o R2_trim.fastq R2.fastq

${star_root}/STAR \
--genomeDir ${genome_dir} \
--readFilesIn R2.fastq \
--outFileNamePrefix star \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outFilterMismatchNoverLmax 0.05 \
--outFilterMatchNmin 16 \
--alignIntronMax 1 \
--limitOutSJcollapsed 5000000 \
--outSAMtype BAM Unsorted \
--runThreadN 12 \
--limitBAMsortRAM 41143265264

## MergeBamAlignment
java -Xmx100g -jar ${dropseq_root}/3rdParty/picard/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=${genome_dir}/Mus_musculus.GRCm38.88.fasta \
UNMAPPED_BAM=filtered.bam \
ALIGNED_BAM=starAligned.out.bam \
OUTPUT=merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

${dropseq_root}/TagReadWithGeneFunction I=merged.bam \
O=star_gene_exon_tagged.bam \
ANNOTATIONS_FILE=${genome_dir}/Mus_musculus.GRCm38.88.gtf
