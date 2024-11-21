rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(RColorBrewer)
library(BuenColors)

meta = readRDS(file = "./Barcode_level_meta_filter.rds")
dge = readRDS(file = "./Barcode_level_dge_filter.rds")

mob = CreateSeuratObject(counts = dge, 
                         assay = "Spatial")
coords.df = data.frame(x = -meta$new_x, 
                       y =-meta$new_y, 
                       stringAsFactors=FALSE)
rownames(coords.df) = meta$degen_barcodes
mob@images$image = new(
  Class = "SlideSeq", 
  assay = "Spatial", 
  key = "image_", 
  coordinates=coords.df
)
mob$log_nCount_Spatial = log10(mob$nCount_Spatial)

mito.genes = grep(pattern = "^mt-", x = rownames(x = mob@assays$RNA), value = T)
mob[['percent.mt']] = PercentageFeatureSet(mob, pattern = "^mt-")
summary(mob$percent.mt)

#### scTransform
mob = SCTransform(mob, assay = "Spatial", ncells = 3000, verbose = TRUE)
mob= RunPCA(mob)
mob = RunUMAP(mob, dims = 1:30)
mob = FindNeighbors(mob, dims = 1:30)
mob = FindClusters(mob, resolution = 0.6, verbose = TRUE)

cluster_colors = colorRampPalette(c("#14517C", "#2F7FC1", "#E7EFFA", "#96C37D", "#F3D266", "#D8383A", "#F7E1ED", "#F8F3F9", "#C497B2", "#A9B8C6"))(length(unique(mob$seurat_clusters)))

cluster_colors2 = c("#14517C", "#E58450", "#F7EEF6", "#B8C874", "#D2E3DA", "#EBA1A9", "#468DC8", "#CAA2BA", "#A9B8C6" )

pdf(file = "./Barcode_Level_0109-01MOB_CellType.pdf", width = 6, height = 6)
DimPlot(mob, reduction = "umap", label = F, cols = cluster_colors2, group.by = "cell_type") + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + coord_fixed() + ggtitle("")
dev.off()

rownames(meta) = meta$degen_barcodes
meta = meta[rownames(mob@meta.data),]
meta$seurat_clusters = mob$seurat_clusters

pdf(file = "./Barcode_Level_0109-01MOB_spatialClusters.pdf", width = 6, height = 6)
ggplot(meta, aes(x = new_x, y = new_y, color = cell_type)) + geom_point(size = 0.01) + theme_minimal() + scale_color_manual(values = cluster_colors2) + xlab("Dim 1") + ylab("Dim 2") + theme(axis.text = element_blank(), legend.position = "none", panel.grid = element_blank()) + coord_fixed() + annotate("segment", x = 6000, xend = (6000 + 200/0.72), y = 2500, yend = 2500, linewidth = 2) 
dev.off()

 

## markers
markers = FindAllMarkers(mob, logfc.threshold = 0.25, verbose = T, only.pos = T, min.pct = 0.25)
markers.filtered = markers %>% filter(p_val_adj <= 0.05)

write.csv(markers.filtered, file = "./Barcode_Level_0109-01MOB_markers_filtered.csv", row.names = F)

mob = FindSpatiallyVariableFeatures(mob, assay = "SCT", slot = "scale.data", features=VariableFeatures(mob)[1:2000], 
                                    selection.method="moransi", x.cuts=100, y.cuts=100)

sp_features=SpatiallyVariableFeatures(mob, selection.method = "moransi")
write.csv(sp_features, file = "./Barcode_level_0109-01SpatiallyVariable_genes.csv")
