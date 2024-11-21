rm(list = ls()); gc()
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(RColorBrewer)
library(BuenColors)
library(mclust)

meta = readRDS(file = "./Barcode_level_meta_filter.rds")
dge = readRDS(file = "./Barcode_level_dge_filter.rds")

hippocampus = CreateSeuratObject(counts = dge, 
                         assay = "Spatial")
coords.df = data.frame(x = -meta$new_x, 
                       y = meta$new_y, 
                       stringAsFactors = FALSE)
rownames(coords.df) = meta$degen_barcodes
hippocampus@images$image = new(
  Class = "SlideSeq", 
  assay = "Spatial", 
  key = "image_", 
  coordinates = coords.df
)

mito.genes = grep(pattern = "^mt-", x = rownames(x = hippocampus@assays$RNA), value = T)
hippocampus[['percent.mt']] = PercentageFeatureSet(hippocampus, pattern = "^mt-")
summary(hippocampus$percent.mt)

hippocampus$log_nCount_Spatial = log10(hippocampus$nCount_Spatial)

### scTransform
hippocampus = SCTransform(hippocampus, assay = "Spatial", ncells = 3000, verbose = TRUE)
hippocampus = RunPCA(hippocampus)
hippocampus = RunUMAP(hippocampus, dims = 1:30)
hippocampus = FindNeighbors(hippocampus, dims = 1:30)
hippocampus = FindClusters(hippocampus, resolution = 0.4, verbose = TRUE)

###

cluster_colors = colorRampPalette(c("#14517C", "#2F7FC1", "#E7EFFA", "#96C37D", "#F3D266", "#D8383A", "#F7E1ED", "#F8F3F9", "#C497B2", "#A9B8C6"))(length(unique(hippocampus$seurat_clusters)))

pdf(file = "./Barcode_Level_0222-01Hippocampus_clusters.pdf", width = 6, height = 6)
DimPlot(hippocampus, reduction = "umap", label = TRUE, cols = cluster_colors) + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + coord_fixed()
dev.off()

rownames(meta) = meta$degen_barcodes
meta = meta[rownames(hippocampus@meta.data),]
meta$seurat_clusters = hippocampus$seurat_clusters

pdf(file = "./Barcode_Level_0222-01Hippocampus_spatialClusters.pdf", width = 6, height = 6)
ggplot(meta, aes(x = -new_x, y = -new_y, color = seurat_clusters)) + geom_point(size = 0.01) + theme_minimal() + scale_color_manual(values = cluster_colors) + xlab("Dim 1") + ylab("Dim 2") + theme(axis.text = element_blank(), legend.position = "none", panel.grid = element_blank()) + coord_fixed() + annotate("segment", x = -8000, xend = -(8000 + 200/0.72), y = -2500, yend = -2500, linewidth = 2) 
dev.off()

# markers
markers = FindAllMarkers(hippocampus, logfc.threshold = 0.25, verbose = T, only.pos  =T, min.pct = 0.25)
markers.filtered = markers %>% filter(p_val_adj <= 0.05)

write.csv(markers.filtered, file = "./Barcode_Level_0222-01_Hippocampus_markers_filtered.csv", row.names = F)

## spatial variable genes
hippocampus = FindSpatiallyVariableFeatures(hippocampus, assay = "SCT", slot = "scale.data", features=VariableFeatures(hippocampus)[1:2000], 
                                    selection.method="moransi", x.cuts=100, y.cuts=100)

sp_features=SpatiallyVariableFeatures(hippocampus, selection.method = "moransi")
write.csv(sp_features, file = "./Barcode_Level_0222-01_Hippocampus_Spatially_Variable_Features.csv")

saveRDS(hippocampus, file = "./Barcode_Level_0222-01Hippocampus.rds")
saveRDS(meta, file = "./Barcode_Level_0222-01Hippocampus_meta.rds")
