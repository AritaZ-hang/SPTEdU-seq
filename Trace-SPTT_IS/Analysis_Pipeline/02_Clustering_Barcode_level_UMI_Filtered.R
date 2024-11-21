rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(mclust)
library(tidyverse)
library(BuenColors)
library(scales)

if(!dir.exists("./Total_UMI_Filtered"))
{
  dir.create("./Total_UMI_Filtered")
}

meta=readRDS(file="./Total_map_meta.RDS")
dge=as.data.frame(fread(file="../total_dge/allmapped_dge.txt.gz"))
dge=column_to_rownames(dge, "GENE")
dge = dge[, meta$barcodes]

injure = CreateSeuratObject(counts = dge, 
                            assay = "Spatial")
coords.df = data.frame(x = meta$x_coords, 
                       y = -meta$y_coords, 
                       stringAsFactors = F)
rownames(coords.df) = meta$barcodes

injure@images$image = new(Class = "SlideSeq", assay = "Spatial", key = "image_", coordinates = coords.df)
injure$log_nCount_Spatial = log10(injure$nCount_Spatial)
injure$UTIs = meta$UTIs
injure$log_UTIs = log10(injure$UTIs)
injure$ratio = meta$ratio
injure$ratio_adjusted=meta$ratio_adjusted

#
injure = subset(injure, `nCount_Spatial` >= 240 & `nFeature_Spatial` >= 125)
# injure = subset(injure, `nCount_Spatial` >= 200 & `nFeature_Spatial` >= 100)


rm(dge); gc()

injure=readRDS(file="./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain.RDS")

#

injure = SCTransform(injure, assay = "Spatial", ncells = 300, verbose = TRUE)
injure = RunPCA(injure)
ElbowPlot(injure, ndims = 50)
injure = RunUMAP(injure, dims = 1:30)
injure = FindNeighbors(injure, dims = 1:30)
injure = FindClusters(injure, resolution = 0.8, verbose = T)

cluster_colors = colorRampPalette(c("#14517C", "#2F7FC1", "#E7EFFA", "#96C37D", "#F3D266", "#D8383A", "#F7E1ED", "#F8F3F9", "#C497B2", "#A9B8C6"))(length(unique(injure$seurat_clusters)))

#
pdf(file="./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_Clusters_nPC30_Res0.8.pdf", width = 8, height = 8)
DimPlot(injure, reduction = "umap", label = TRUE, cols = cluster_colors, raster = F, pt.size= 0.1) + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + coord_fixed()
dev.off()

rownames(meta) = meta$barcodes
meta = meta[rownames(injure@meta.data),]
meta$seurat_clusters=injure$seurat_clusters
meta$nCount_SCT=injure$nCount_SCT
meta$nFeature_SCT=injure$nFeature_SCT

pdf(file = "./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_spatialClusters_nPC30_Res0.8.pdf", width = 8, height =8)
ggplot(meta, aes(x = -x_coords, y = -y_coords, color = seurat_clusters)) + geom_point(size = 0.1) + theme_minimal() + scale_color_manual(values = cluster_colors) + xlab("Dim 1") + ylab("Dim 2") + theme(axis.text = element_blank(), legend.position = "none", panel.grid = element_blank()) + coord_fixed() + annotate("segment", x = -9000, xend = -(9000 + 200/0.72), y = -2500, yend = -2500, linewidth = 2) 
dev.off()

#
markers = FindAllMarkers(injure, logfc.threshold = 0.2, verbose = T, only.pos = T, min.pct = 0.2)
markers.filtered = markers %>% filter(p_val_adj <= 0.05)
table(markers.filtered$cluster)
write.csv(markers.filtered, file = "./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_markers_filtered_nPC30_Res0.8.csv", row.names = F)

saveRDS(injure, file = "./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_nPC30_Res0.8.RDS")
saveRDS(meta, file = "./Total_UMI_Filtered/Total_UMI_Filtered_0401-01InjuryBrain_meta_nPC30_Res0.8.RDS")