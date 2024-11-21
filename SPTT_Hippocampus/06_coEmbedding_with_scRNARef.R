rm(list = ls())
setwd("/your/workdir/")
library(data.table)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

list.files()

hippocampus = readRDS(file = "./Barcode_Level_0222-01Hippocampus.rds")
hippocampus_ref = readRDS(file="./scRNA_hippocampus_reference.RDS")

hippocampus_ref = SCTransform(hippocampus_ref, assay = "RNA", ncells = 3000, verbose = TRUE)
hippocampus_ref = RunPCA(hippocampus_ref)
ElbowPlot(hippocampus_ref, dims = 1:30)
hippocampus_ref = RunUMAP(hippocampus_ref, dims = 1:30)
hippocampus_ref = FindNeighbors(hippocampus_ref, dims = 1:30)
hippocampus_ref = FindClusters(hippocampus_ref, resolution = 0.6, verbose = T)

obj_list = list(hippocampus, hippocampus_ref)

# intergrate scRNA-seq & Spatial
features = SelectIntegrationFeatures(obj_list)
integration.anchors=FindIntegrationAnchors(obj_list, anchor.features = features)
combined = IntegrateData(anchorset = integration.anchors)

DefaultAssay(combined) = "integrated"

combined = ScaleData(combined, verbose = F)
combined = RunPCA(combined, verbose = T)
combined = RunUMAP(combined, reduction="pca", dims=1:30)
combined = FindNeighbors(combined, reduction = "pca")
combined = FindClusters(combined, resolution=0.5)

DimPlot(combined, reduction = "umap", group.by = "celltype", raster=F) + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + ggtitle("") + guides(color = guide_legend(override.aes = list(size=5)))
ggsave(filename="./hippocampus_coembed_with_scRNARef.pdf", width = 7, height = 6)
ggsave(filename='./hippocampus_coembed_with_scRNARef.png', width = 7, height = 6)

saveRDS(combined, file = "./combined_with_scRNA_Ref.RDS")

