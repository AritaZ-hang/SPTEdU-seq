rm(list = ls())
setwd("/your/workdir/")
library(data.table)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
source("plot_utils.R")

spatial = readRDS(file = "./Barcode_level_0109-01MOB.rds")
scRNA = readRDS(file = "./scRNA_MOB_reference.rds")

spatial_counts = GetAssayData(spatial, slot = "counts", assay = "SCT")
scRNA_counts = GetAssayData(scRNA, slot = "counts", assay = "RNA")

table(scRNA$celltype)

# TPM normalization
norm_spatial = spatial_counts %>% rowSums() %>% {1e6 *. / sum(.)}
norm_scRNA = scRNA_counts %>% rowSums() %>% {1e6 *. / sum(.)}

intersect_genes = intersect(rownames(spatial_counts), rownames(scRNA_counts))

rna_tpm = tibble(
  "Spatial" = norm_spatial[intersect_genes], 
  "scRNA" = norm_scRNA[intersect_genes]
)%>%
  filter(Spatial != 0, scRNA !=0) %>%
  mutate(across(everything(), ~log10(1+1e6 * .x / sum(.x))))

P1 = ggplot(rna_tpm, aes(x = Spatial, y = scRNA)) + geom_point(size = 0.3) + theme_classic() + xlab("SPTT(all spatial clusters)") + ylab("scRNA Reference(all clusters)")+ coord_equal() + stat_scatter_density(size = 0.3) + scale_color_viridis_c(option = 'G') + 
  ggpubr::stat_cor(aes(label= paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman", alternative = "two.sided", cor.coef.name = "rho")+ geom_smooth(method = "lm", se = TRUE, color = "#47A1A2", formula = y~x)  + labs(color = "Density")
P1
ggsave(P1, filename = "./TPM_correlation_with_scRNA-seq_Ref.pdf", width = 4, height = 4)
