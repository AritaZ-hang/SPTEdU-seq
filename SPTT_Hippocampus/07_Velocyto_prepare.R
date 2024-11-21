rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(velocyto.R)
library(tidyverse)
library(pagoda2)
library(SeuratDisk)

list.files()
#############
hippocampus = readRDS(file = "./Barcode_Level_0222-01Hippocampus.rds")
hippocampus.meta = readRDS(file = "./Barcode_level_meta_filter.rds")
rownames(hippocampus.meta) = hippocampus.meta$degen_barcodes
hippocampus.meta = hippocampus.meta[colnames(hippocampus),]
hippocampus_embedding = data.frame(
  barcodes = colnames(hippocampus), 
  x_UMAP = hippocampus@reductions$umap@cell.embeddings[,1],
  y_UMAP = hippocampus@reductions$umap@cell.embeddings[,2], 
  x_coords = hippocampus.meta$new_x, 
  y_coords = hippocampus.meta$new_y, 
  seurat_clusters = hippocampus$seurat_clusters
)
rownames(hippocampus_embedding) = hippocampus_embedding$barcodes
saveRDS(hippocampus_embedding, file = "./Hippocampus_embedding.rds")

# rename the barcodes

hippocampus$x_coords = hippocampus.meta$new_x
hippocampus$y_coords = hippocampus.meta$new_y

hippocampus.rename=RenameCells(hippocampus, 
                               new.names = paste0("Hippocampus_Top5W_degen_V3CE0:", gsub("-1", "", colnames(hippocampus))))

saveRDS(hippocampus.rename, file = "./Barcode_Level_0222-01Hippocampus_Velocyto.rds")

ldat = read.loom.matrices("./Velocyto/Hippocampus_Top5W_degen_V3CE0.loom")
emat = ldat$spliced; nmat = ldat$unspliced
length(intersect(colnames(emat), colnames(hippocampus.rename)))

# export the data
if(!dir.exists("./Velocyto/Velocyto_data/"))
{
  dir.create("./Velocyto/Velocyto_data/")
}

# umap embedding
hippocampus_umap_embed=Embeddings(hippocampus.rename, reduction = "umap") 
# spatial coords
hippocampus_spatial_coords=data.frame(x_coords=-hippocampus.rename$x_coords,
                                      y_coords=-hippocampus.rename$y_coords)
rownames(hippocampus_spatial_coords)=colnames(hippocampus.rename)
hippocampus_spatial_coords=as.matrix(hippocampus_spatial_coords)
# barcodes
hippocampus_barcode=Cells(hippocampus.rename)
# cluster
hippocampus_cluster=as.data.frame(hippocampus.rename@meta.data[, "seurat_clusters"])
colnames(hippocampus_cluster)=c("cluster")
rownames(hippocampus_cluster)=colnames(hippocampus.rename)

write.csv(hippocampus_umap_embed, file = "./Velocyto/Velocyto_data/cell_embeddings.csv")
write.csv(hippocampus_spatial_coords, file = "./Velocyto/Velocyto_data/spatial_coords.csv")
write.csv(hippocampus_barcode,file="./Velocyto/Velocyto_data/cellID_obs.csv",row.names=F)
write.csv(hippocampus_cluster,file="./Velocyto/Velocyto_data/cell_clusters.csv")
