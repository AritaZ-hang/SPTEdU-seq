#!/bin/bash

workdir=`pwd`
datadir=/your/data/dir/
RefDir=/your/ref/dir/

velocyto run -m ${RefDir}/modified_mm10_rmsk.gtf --outputfolder Velocyto_Rmsk --verbose ${datadir}/input.bam ${RefDir}/Mus_musculus.GRCm38.88.gtf
