#!/bin/bash

workdir=`pwd`
dropseq_root=/your/dropseq/root/
scripts_dir=/your/scripts/dir/
split_bam=${scripts_dir}/splitbam.py

trimReadsDir=trimreads_10K
mkdir -p $trimReadsDir
cd $trimReadsDir

tmpDir=tmp
mkdir -p $tmpDir

BAM_FILE=${workdir}/use.bam
SAM_HEADER=${workdir}/SAM_header
BARCODE_TAG=XC
UMI_TAG=XM

samtools sort -t $BARCODE_TAG $BAM_FILE -o use.sorted.bam
python -u $split_bam $BARCODE_TAG use.sorted.bam ${tmpDir}/

file3='_sorted.bam'
file4='_sample.sam'
file5='_sample.bam'
file5='_out.bam'
file6='_dge.txt.gz'
file7='_out_cell_readcounts.txt'
file8='_dge.summary.txt'

for k in 1000 2000 4000 6000 8000 10000
do
    cd $tmpDir
    mkdir tmp2
    for j in *.bam
    do
        samtools view -@ 16 $j | shuf -n $k > tmp2/filtered_SAM_body
        cat $SAM_HEADER tmp2/filtered_SAM_body > tmp2/$j$file4
        samtools view -@ 16 -b tmp2/$j$file4 > tmp2/$j$file5 &&rm tmp2/filtered_SAM_body tmp2/$j$file4
        samtools sort -@ 16 tmp2/$j$file5 > tmp2/$j$file3 && rm tmp2/$j$file5
        samtools index tmp2/$j$file3
    done

    samtools merge ../_trim$k$file5 tmp2/*_sorted.bam && rm -r tmp2
    cd ../
    
    ${dropseq_root}/DigitalExpression -m 16g \
    I=_trim$k$file5 \
    CELL_BARCODE_TAG=$BARCODE_TAG \
    MOLECULAR_BARCODE_TAG=$UMI_TAG \
    O=_trim$k$file6 \
    SUMMARY=_trim$k$file8 \
    NUM_CORE_BARCODES=1000 \
    LOCUS_FUNCTION_LIST=INTRONIC \
    TMP_DIR=.
done
