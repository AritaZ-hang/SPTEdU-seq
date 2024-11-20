#!/bin/bash

cd `pwd`
mkdir tmp
tmpdir=tmp
sample_name=$(basename `pwd`)

dropseq_root=/your/dropseq/root/
script_dir=/your/scripts/root/

##### customed scripts #####

barcode_manipulate_script=${script_dir}/SPTT_Barcode_Manipulate.R
degen_barcode_script=${script_dir}/SPTT_degen_barcode.py
mapping_dict_script=${script_dir}/SPTT_Mapping_dict.py
bin_making_script=${script_dir}/SPTT_Bin_Making.R

workdir=`pwd`

top_num=200000
bc_loc_file=bc_loc.csv

###################################
########## grep barcodes ##########
###################################
Rscript $barcode_manipulate_script grep out_cell_readcounts.txt $top_num $workdir
###################################
###### barcode degen ##############
###################################
python -u $degen_barcode_script -bc_loc $bc_loc_file -seq_barcodes "Top"$top_num"_barcodes.txt" -workdir $workdir
Rscript $barcode_manipulate_script mapping ${workdir}/barcode_mapping.txt.gz $workdir
python -u $mapping_dict_script ${workdir}/barcode_mapping.csv $workdir
###################################
####### make bin level coords #####
###################################
Rscript $bin_making_script ${workdir}/barcode_mapping.txt.gz $workdir
python -u $mapping_dict_script ${workdir}/bin10_coordinates.csv $workdir
python -u $mapping_dict_script ${workdir}/bin50_coordinates.csv $workdir
python -u $mapping_dict_script ${workdir}/bin100_coordinates.csv $workdir