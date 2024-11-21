rm(list = ls())
setwd("/your/workdir/")
library(data.table)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
source('plot_utils.R')

# Barcode Level
rep1_dge = as.data.frame(fread(file = "./Rep1_MOB_top1w_dge.txt.gz"))
rep1_dge = column_to_rownames(rep1_dge, var = "GENE")

rep2_dge = as.data.frame(fread(file = "./Rep2_MOB_top1w_dge.txt.gz"))
rep2_dge = column_to_rownames(rep2_dge, var = "GENE")

# intersect genes
intersect_genes = intersect(rownames(rep1_dge), rownames(rep2_dge))

rep1_dge = rep1_dge[intersect_genes,]
rep2_dge = rep2_dge[intersect_genes,]

# TPM normalization
norm_rep1 = rep1_dge %>% rowSums() %>% {1e6 *. / sum(.)}
norm_rep2 = rep2_dge %>% rowSums() %>% {1e6 *. / sum(.)}

# rna tpm
rna_tpm = tibble(
  "Replicate1" = norm_rep1, 
  "Replicate2" = norm_rep2
)%>%
  filter(Replicate1 != 0, Replicate2 !=0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

P1 = ggplot(rna_tpm, aes(x = Replicate1, y = Replicate2)) + geom_point(size = 0.3) + theme_classic() + xlab("MOB Replicate 1 Barcode level") + ylab("MOB Replicate 2 Barcode level")+ coord_equal() + stat_scatter_density(size = 0.3)+ scale_color_viridis_c(option = 'G') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman", alternative = "two.sided", cor.coef.name = "rho")+ geom_smooth(method = "lm", se = TRUE, color = "#47A1A2", formula = y~x)  + labs(color = "Density")
P1

ggsave(P1, filename="./TPM_correlation_with_MOB_replicates_barcode_level.pdf", width = 4, height = 4)
