#!/usr/bin/Rscript

rm(list = ls()); gc()
library(data.table)
argv = commandArgs(trailingOnly=T)

command = argv[1]
if(command == "grep")
{
    readsfile=as.data.frame(fread(file = argv[2]))[, 1:2]
    colnames(readsfile) = c("reads", "barcode")
    barcodes_num = as.integer(argv[3])
    savedir = argv[4]

    selected_barcodes = readsfile[["barcode"]][1:barcodes_num]

    write.table(selected_barcodes, file = paste0(savedir, "/", "Top", as.character(barcodes_num), "_barcodes.txt"), row.names = F, col.names =F , quote=F)
}

if(command == "mapping")
{
    barcode_mapping = as.data.frame(fread(file = argv[2]))
    savedir = argv[3]

    barcode_mapping_pre_and_now = barcode_mapping[, 1:2]
    write.csv(barcode_mapping_pre_and_now, file = paste0(savedir, "/", "barcode_mapping.csv"), row.names = F)
}