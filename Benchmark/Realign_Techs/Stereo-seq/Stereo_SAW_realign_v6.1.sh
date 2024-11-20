#!/bin/bash

sample=stereo_MOB # Get tissue parameter from command line argument
species=Mus_musculus # Get reference parameter from command line argumentt
tissueType=MOB # Tissue name (it just a label)
casno=SAW_image_v6.1 # Output saw-tools pipeline folder
  
ulimit -n 4096
ulimit -v 33170449147

dataDir=/your/data/Stereo_MOB/
rawData=/your/data/Stereo_MOB/
scriptDir=/SAW_ver6.1/
saw_dir=/your/saw/dir/
outDir=${dataDir}/$casno
referenceDir=/your/reference/dir/

NUMBA_CACHE_DIR=$scriptDir/temp
export SINGULARITY_BIND=$dataDir,$outDir,$rawData

bash ${scriptDir}/stereo-pipeline_v6.0.sh \
-genomeSize 5 \
-splitCount 16 \
-maskFile ${rawData}/FP200009107_E414.barcodeToPos.h5 \
-fq1 ${rawData}/FP200009107_E414_R1.fq.gz \
-fq2 ${rawData}/FP200009107_E414_R2.fq.gz \
-refIndex ${referenceDir}/STAR_SJ100_ver6.1/ \
-speciesName $species \
-tissueType $tissueType \
-annotationFile ${referenceDir}/Mus_musculus.GRCm38.88.gtf \
-outDir $outDir \
-doCellBin N \
-rRNAremove N \
-threads 16 \
-sif ${saw_dir}/SAW_6.1.sif \
-imageRecordFile ${rawData}/FP200009107_E414_SC_20230614_193616_2.1.0.ipr \
-imageCompressedFile ${rawData}/FP200009107_E414_SC_20230614_193616_2.1.0.tar.gz
