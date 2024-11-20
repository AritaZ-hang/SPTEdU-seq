#!/bin/bash

workdir=`pwd`
spaceranger=/your/path/to/spaceranger/
dropseq_root=/your/dropseq/root/

fastq_dir=Visium_Mouse_Olfactory_Bulb_fastqs
ref_dir=/your/spaceranger/ref/dir/

$spaceranger count --id="10X_MOB" \
--description="Adult Mouse Olfactory Bulb" \
--transcriptome=${ref_dir} \
--fastqs=${fastq_dir} \
--image=Visium_Mouse_Olfactory_Bulb_image.tif \
--slide=V10N30-322 \
--area=A1 \
--localcores=16 \
--localmem=128 \
--create-bam=true \
--output-dir=10X_MOB_Remapped


input_bam=${workdir}/10X_MOB_Remapped/outs/possorted_genome_bam.bam

${dropseq_root}/DigitalExpression -m 8g \
I=$input_bam \
CELL_BARCODE_TAG=CB \
MOLECULAR_BARCODE_TAG=UB \
O=10X_Visium_dge.txt.gz \
SUMMARY=10X_Visium_dge.summary.txt \
NUM_CORE_BARCODES=1000 \
LOCUS_FUNCTION_LIST=INTRONIC \
STRAND_STRATEGY=BOTH \
TMP_DIR=.


