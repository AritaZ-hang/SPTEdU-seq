rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(mclust)
library(tidyverse)
library(AUCell)
library(GSEABase)

if(!dir.exists("./Total_UMI_Filtered"))
{
  dir.create("./Total_UMI_Filtered")
}

injure=readRDS(file="./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_nPC30_Res0.8.RDS")
meta=readRDS(file="./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_meta_nPC30_Res0.8.RDS")

#
exprMatrix=GetAssayData(injure, slot = "counts")

# Gene sets' activity
injury_genes=(c("Adamts1", "Atf3", "Ccl2", "Ccnd1", "Cd68", "Cebpd", "Cyba", "Fn1", "Gal", "Gap43", "Hmox1", "Hspb1", "Igfbp2", "Jun", "Junb", "Fos", "Lgals1", "Neat1", "Socs3", "Tnc", "S100a10", "Timp1"))
proliferation_genes=(c("Mki67", "Ccnd1", "E2f1", "Gab2", "Mcm6", "Pcna", "Top2a"))
migration_genes=(unique(as.data.frame(fread(file="./migration_score_genes.csv", header=F))$V1))

injury_genes=intersect(injury_genes, rownames(injure))
proliferation_genes=intersect(proliferation_genes, rownames(injure))
migration_genes=intersect(migration_genes, rownames(injure))

injury_gene_sets=GeneSet(injury_genes, setName="injury")
proliferation_gene_sets=GeneSet(proliferation_genes, setName="proliferation")
migration_gene_sets=GeneSet(migration_genes, setName="migration")

##
cells_rankings=AUCell_buildRankings(exprMatrix,plotStats = FALSE)
cells_AUC_injury=AUCell_calcAUC(geneSets = injury_gene_sets, cells_rankings)
cells_AUC_proliferation=AUCell_calcAUC(geneSets=proliferation_gene_sets, cells_rankings)
cells_AUC_migration=AUCell::AUCell_calcAUC(geneSets = migration_gene_sets, cells_rankings)

save(cells_AUC_injury, cells_AUC_proliferation, cells_AUC_migration, file="./Total_UMI_Filtered/AUCs.RData")

#### check thresholds
cells_assignment_injury = AUCell_exploreThresholds(cells_AUC_injury, plotHist = TRUE, assign=TRUE)
cells_assignment_proliferation = AUCell_exploreThresholds(cells_AUC_proliferation, plotHist = TRUE, assign=TRUE)
cells_assignment_migration = AUCell_exploreThresholds(cells_AUC_migration, plotHist = TRUE, assign=TRUE)
save(cells_assignment_injury, cells_assignment_proliferation, cells_assignment_migration, file="./Total_UMI_Filtered/AUCs_assignments.RData")

####

AUC_injury=as.data.frame(t(getAUC(cells_AUC_injury)))
AUC_proliferation=as.data.frame(t(getAUC(cells_AUC_proliferation)))
AUC_migration=as.data.frame(t(getAUC(cells_AUC_migration)))

AUC_injury = AUC_injury[meta$barcodes,]
AUC_proliferation = AUC_proliferation[meta$barcodes,]
AUC_migration = AUC_migration[meta$barcodes,]

meta$AUC_injury=AUC_injury
meta$AUC_proliferation=AUC_proliferation
meta$AUC_migration=AUC_migration

saveRDS(meta, file="./Total_UMI_Filtered/AUC_meta_raw.rds")
