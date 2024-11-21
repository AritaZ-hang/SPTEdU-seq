rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(mclust)
library(tidyverse)
library(RColorBrewer)
library(BuenColors)
library(scales)
library(CytoTRACE2)
library(ggsignif)
library(ggpubr)
library(ggnewscale)

seob = readRDS(file = "./seob.RDS")

cytotrace2_res = cytotrace2(seob, 
                                  is_seurat = TRUE, 
                                  slot_type = "counts", 
                                  species = "mouse")

meta = cytotrace2_res@meta.data

saveRDS(cytotrace2_res, file = "./CytoTrace2.RDS")
