#!/bin/bash

cd `pwd`
tmpDir=tmp
mkdir -p $tmpDir
sample_name=$(basename `pwd`)
dropseq_root=/your/dropseq/root/
STAR_root=/your/star/root/
seqtk_root=/your/seqtk/root/
scripts_dir=/your/scripts/root/
genome_dir=/your/path/to/genome/
addTagScript=${scripts_dir}/SPTT_barcode_and_UMI_deal.py

java -jar ${dropseq_root}/3rdParty/picard/picard.jar FastqToSam F1=H_R1_withoutlinker.fastq.gz F2=H_R2_withoutlinker.fastq.gz O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name

###############################
###### Tag Barcode & UMI ######
###############################
${dropseq_root}/TagBamWithReadSequenceExtended \ BASE_RANGE=1-22 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=true TAG_NAME=BC NUM_BASES_BELOW_QUALITY=4 INPUT=H.bam OUTPUT=${tmpDir}/H1.bam COMPRESSION_LEVEL=5

###############################
######### BAM FILTER ##########
###############################
${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=${tmpDir}/H1.bam OUTPUT=${tmpDir}/H2.bam && rm ${tmpDir}/H1.bam

java -Xmx100g -jar ${dropseq_root}/3rdParty/picard/picard.jar SamToFastq INPUT=${tmpDir}/H2.bam FASTQ=R2.fastq

##############################
###### trim PolyA & rev ######
##############################

trimDir=trim
mkdir -p $trimDir

cutadapt \
-a A{10} \
-j 0 \
--minimum-length=15 \
-o R2_trim.fastq \
--cores 12 R2.fastq > ${trimDir}/cutadapt_polyA_trim_report.txt && rm R2.fastq

$seqtk_root seq -Ar R2_trim.fastq > R2_trim_reverse.fastq && rm R2_trim.fastq

##############################
###### Align and DGE #########
##############################

dgeDir=dge
mkdir -p $dgeDir

$STAR_root \
--genomeDir ${genome_dir} \
--readFilesIn R2_trim_reverse.fastq \
--outFileNamePrefix star \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outFilterMismatchNoverLmax 0.05 \
--outFilterMatchNmin 16 \
--alignIntronMax 1 \
--limitOutSJcollapsed 5000000 \
--outSAMtype BAM Unsorted \
--runThreadN 12 \
--limitBAMsortRAM 41143265264 && rm R2_trim_reverse.fastq

java -Xmx100g -jar ${dropseq_root}/3rdParty/picard/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=${genome_dir}/Mus_musculus.GRCm38.88.fasta \
UNMAPPED_BAM=${tmpDir}/H2.bam \
ALIGNED_BAM=starAligned.out.bam \
OUTPUT=merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false && rm starAligned.out.bam ${tmpDir}/H2.bam

${dropseq_root}/TagReadWithGeneFunction I=merged.bam \
O=${dgeDir}/star_gene_exon_tagged.bam \
ANNOTATIONS_FILE=${genome_dir}/Mus_musculus.GRCm38.88.gtf && rm merged.bam

#################################################
####### Deal with Barcode Tag & UMI tag #########
#################################################
python -u $addTagScript -in_bam ${dgeDir}/star_gene_exon_tagged.bam -out_bam ${dgeDir}/star_gene_exon_tagged.use.bam -linker no && rm ${dgeDir}/star_gene_exon_tagged.bam

${dropseq_root}/DigitalExpression -m 8g \
I=${dgeDir}/star_gene_exon_tagged.use.bam \
CELL_BARCODE_TAG=XC \
MOLECULAR_BARCODE_TAG=XM \
O=${dgeDir}/_dge.txt.gz \
SUMMARY=${dgeDir}/_dge.summary.txt \
NUM_CORE_BARCODES=10000 \
LOCUS_FUNCTION_LIST=INTRONIC \
STRAND_STRATEGY=BOTH \
TMP_DIR=.

${dropseq_root}/BamTagHistogram I=${dgeDir}/star_gene_exon_tagged.use.bam O=${dgeDir}/out_cell_readcounts.txt TAG=XC && gzip ${dgeDir}/out_cell_readcounts.txt

rm -r ${tmpDir}
echo "Done"