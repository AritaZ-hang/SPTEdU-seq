#!/bin/bash


sample=SRRNumber # The SRR Number

dropseq_root=/your/dropseq/root/
bbmap_root=/your/bbmap/root/
star_root=/your/star/root/
genome_dir=/your/genome/dir/
workdir=`pwd`${sample}/
script_dir=/your/scripts/dir/
BarcodeScript=${script_dir}/Velocyto_preprocess.py

cd $workdir

R1=${workdir}/H_R1.fastq.gz
R2=${workdir}/H_R2.fastq.gz

#####################################
####### Trimming R2  && fastqc#######
#####################################

pre_fastqcDir=preTrim_fastqC
mkdir -p $pre_fastqcDir
post_fastqcDir=postTrim_fastqC
mkdir -p $post_fastqcDir

fastqc \
--outdir $pre_fastqcDir \
--threads 12 \
$R2

trimDir=Trim
mkdir -p $trimDir

# TSO & polyA trimming
cutadapt \
--minimum-length 10 \
-A A{100} \
-G CCCATGTACTCTGCGTTGATACCACTGCTT \
--pair-filter=any \
-o H_R1_Atrimmed.fastq.gz \
-p H_R2_Atrimmed.fastq.gz \
--cores 12 H_R1.fastq.gz H_R2.fastq.gz > ${trimDir}/cutadapt_polyA_report.txt

# TSO & polyG trimming
cutadapt \
--minimum-length 10 \
-A G{100} \
-G AAGCAGTGGTATCAACGCAGAGTACATGGG \
--pair-filter=any \
-o H_R1_final.fastq.gz \
-p H_R2_final.fastq.gz \
--cores 12 H_R1_Atrimmed.fastq.gz H_R2_Atrimmed.fastq.gz > ${trimDir}/cutadapt_polyG_report.txt

mv H_R1_final.fastq.gz ${trimDir}
mv H_R2_final.fastq.gz ${trimDir}

fastqc \
--output $post_fastqcDir \
--threads 12 \
${trimDir}/H_R2_final.fastq.gz


##################################
######## Extract barcodes ########
##################################
extractDir=BC_Extract
mkdir -p $extractDir

umi_tools extract --stdin ${trimDir}/H_R1_final.fastq.gz --extract-method=regex --bc-pattern="(?P<cell_1>.{16})(?P<umi_1>.{12}).*" --stdout ${extractDir}/H_R1_extract.fastq.gz --read2-in ${trimDir}/H_R2_final.fastq.gz --read2-out ${extractDir}/H_R2_extract.fastq.gz -L ${extractDir}/extract.log


#######################
######## Align ########
#######################
alignDir=R2_align
mkdir -p $alignDir

${star_root}/STAR \
--genomeDir ${genome_dir} \
--readFilesIn ${extractDir}/H_R2_extract.fastq.gz \
--readFilesCommand zcat \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outFilterMismatchNoverLmax 0.05 \
--outFilterMatchNmin 16 \
--limitOutSJcollapsed 5000000 \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 12 \
--limitBAMsortRAM 41143265264

mv Aligned.sortedByCoord.out.bam ${alignDir}
samtools view -b -q 250 ${alignDir}/Aligned.sortedByCoord.out.bam > STRS.bam

#############################
####### grep barcodes #######
#############################
sed -i 's/^\([^ \t]\+\)/CB:Z:\1/' visium-v1.txt
CB_WHITELIST=visium-v1.txt

python -u $BarcodeScript -technique strs -in_bam STRS.bam -out_bam STRS.processed.bam
mv STRS.processed.bam ${alignDir}
samtools sort -@ 12 -o ${alignDir}/STRS.processed.sort.bam ${alignDir}/STRS.processed.bam
samtools index -@ 12 ${alignDir}/STRS.processed.sort.bam

BAM_FILE=${alignDir}/STRS.processed.sort.bam
mkdir temp
samtools view -@ 16 -H $BAM_FILE > SAM_header
samtools view -@ 16 $BAM_FILE | LC_ALL=C grep -F -f $CB_WHITELIST > filtered_SAM_body
cat SAM_header filtered_SAM_body > filtered.sam
samtools view -@ 16 -b filtered.sam > ${alignDir}/STRS.final.bam && rm filtered.sam filtered_SAM_body

###############################
######## generate DGE #########
###############################
dgeDir=dge
mkdir -p $dgeDir

${dropseq_root}/TagReadWithGeneFunction I=${alignDir}/STRS.final.bam O=${dgeDir}/STRS.gene.tagged.bam \
ANNOTATIONS_FILE=${genome_dir}/Mus_musculus.GRCm38.88.gtf

${dropseq_root}/DigitalExpression -m 8g \
I=${dgeDir}/STRS.gene.tagged.bam \
CELL_BARCODE_TAG=CB \
MOLECULAR_BARCODE_TAG=UB \
O=STRS_top1w_dge.txt.gz \
SUMMARY=STRS_top1w_dge.summary.txt \
NUM_CORE_BARCODES=10000 \
LOCUS_FUNCTION_LIST=INTRONIC \
STRAND_STRATEGY=BOTH \
TMP_DIR=.
