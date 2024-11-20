#!/bin/bash

workdir=`pwd`
ref_bed=mm10_protein_coding_genes.bed

input_bam=${workdir}/gene.subsample.100k.sort.bam
geneBody_coverage.py -r $ref_bed -i $input_bam -o genebody_coverage
