rm(list = ls())
setwd("/your/workdir/")
library(data.table)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

list.files()

mob = readRDS(file = "./Barcode_level_0109-01MOB.rds")
mob_ref = readRDS(file = "./scRNA_MOB_reference.rds") 
# previously analyzed by lognormalized, now SCTransform
mob_ref = SCTransform(mob_ref, assay = "RNA", ncells = 3000, verbose=TRUE)
mob_ref = RunPCA(mob_ref)
ElbowPlot(mob_ref, ndims = 50)
mob_ref = RunUMAP(mob_ref, dims = 1:30)
mob_ref = FindNeighbors(mob_ref, dims = 1:30)
mob_ref = FindClusters(mob_ref, resolution = 0.6, verbose = TRUE)


obj_list = list(mob, mob_ref)

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

combined$technique="scRNA"
combined$technique[which(combined$orig.ident == "SeuratProject")] = "Spatial"
combined$simplified_celltype = gsub("_[0-9]+", "", combined$celltype)

colors = brewer.pal(n = 12,name="Paired")

#DimPlot(combined, reduction="umap", group.by = "technique")
DimPlot(combined, reduction = "umap", group.by = "simplified_celltype", raster=F) + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + ggtitle("") + guides(color = guide_legend(override.aes = list(size=5)))
ggsave(filename="./coembed_with_scRNARef.pdf", width = 7, height = 6)
ggsave(filename='./coembed_with_scRNARef.png', width = 7, height = 6)

# correlation heatmap
var_genes = VariableFeatures(combined)
spatial_clusters_avg = AverageExpression(mob, slot = "counts", group.by = "cell_type")[[1]]
scRNA_clusters_avg = AverageExpression(mob_ref, slot = "counts", group.by = "celltype")[[1]]
spatial_clusters_avg = spatial_clusters_avg[var_genes,]
scRNA_clusters_avg = scRNA_clusters_avg[var_genes,]
cor_mat = matrix(nrow = dim(spatial_clusters_avg)[2], ncol = dim(scRNA_clusters_avg)[2])
for(i in 1:dim(spatial_clusters_avg)[2])
{
  for(j in 1:dim(scRNA_clusters_avg)[2])
  {
    cor_mat[i,j] = cor(spatial_clusters_avg[,i], scRNA_clusters_avg[,j], method = "spearman")
  }
}
cor_mat = as.data.frame(cor_mat)
colnames(cor_mat) = colnames(scRNA_clusters_avg)
rownames(cor_mat) = colnames(spatial_clusters_avg)
heatmap_colors = jdb_palette("brewer_yes", type = "continuous")
p1 = pheatmap::pheatmap(cor_mat, scale="none", border = "white", color = heatmap_colors, cluster_rows = T, cluster_cols = T, clustering_method = "ward.D2", )
cor_mat_plot = cor_mat
cor_mat_plot = cor_mat_plot[p1$tree_row$order, p1$tree_col$order]
pheatmap::pheatmap(cor_mat_plot, scale="none", border = "white", color = heatmap_colors, cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", filename = "./combined_integrated_spearman_correlation_clustered_no_tree.pdf", width = 10, height = 4)

saveRDS(combined, file = "./combined_with_scRNA_Ref.RDS")

