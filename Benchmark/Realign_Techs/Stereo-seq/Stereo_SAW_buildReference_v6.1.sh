#!/bin/bash


referenceDir=/your/reference/dir/
saw_dir=/your/saw/dir/
mkdir $referenceDir/STAR_SJ100_ver6.1
export SINGULARITY_BIND=$referenceDir

singularity exec ${saw_dir}/SAW_6.1.sif mapping \
	--runMode genomeGenerate \
	--genomeDir ${referenceDir}/STAR_SJ100_ver6.1 \
	--genomeFastaFiles ${referenceDir}/Mus_musculus.GRCm38.88.fasta \
	--sjdbGTFfile ${referenceDir}/Mus_musculus.GRCm38.88.gtf \
	--sjdbOverhang 99 \
	--runThreadN 12
