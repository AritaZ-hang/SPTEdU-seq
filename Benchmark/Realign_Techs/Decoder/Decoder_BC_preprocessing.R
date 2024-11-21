rm(list = ls()); gc()
setwd("/workdir/")
library(data.table)
library(dplyr)
library(stringr)

list.files()

lib_barcode_X = as.data.frame(fread(file = "./lib_barcode_X.txt", header = F))
lib_barcode_Y = as.data.frame(fread(file = "./lib_barcode_Y.txt", header = F))

lib_barcode_X = lib_barcode_X$V1
lib_barcode_Y = lib_barcode_Y$V1

lib_barcode_X.df = data.frame(bc_raw = lib_barcode_X, 
                            bc_process = lib_barcode_X)
lib_barcode_X.df = rbind(lib_barcode_X.df, c("nobc", "nobc"))
write.table(lib_barcode_X.df, file = "./lib_barcode_X_df.txt", quote = F, row.names = F)

lib_barcode_Y.df = data.frame(bc_raw = lib_barcode_Y, 
                              bc_process = lib_barcode_Y)
lib_barcode_Y.df = rbind(lib_barcode_Y.df, c("nobc", "nobc"))
write.table(lib_barcode_Y.df, file = "./lib_barcode_Y_df.txt", quote = F, row.names = F)