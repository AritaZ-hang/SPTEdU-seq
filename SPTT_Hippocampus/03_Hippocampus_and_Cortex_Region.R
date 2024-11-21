rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(RColorBrewer)
library(BuenColors)
library(mclust)

hippocampus = readRDS(file = "./Barcode_Level_0222-01Hippocampus.rds")
meta = readRDS(file = "./Barcode_Level_0222-01Hippocampus_meta.rds")

###### only hippocampus region
hippocampus_region = meta[which(meta$seurat_clusters %in% c("5", "7", "9")),]
sp.loc = hippocampus_region[, c("new_x", "new_y")]
sp.loc$new_x = as.integer(sp.loc[["new_x"]])
sp.loc$new_y = as.integer(sp.loc[["new_y"]])
nnDist = FNN::get.knnx(sp.loc[, c("new_x", "new_y")], sp.loc[, c("new_x", "new_y")], k = 20)
loc.mat = cbind(nnDist$nn.index, round(nnDist$nn.dist, 2))
average_dist = function(x, loc.mat)
{ 
  return(mean(loc.mat[x, 22:dim(loc.mat)[2]]))
}
nrows = 1:nrow(loc.mat)
beads_average_dist = lapply(nrows, average_dist, loc.mat)
beads_average_dist = do.call(c, beads_average_dist)
hippocampus_region$hippocampus_dist = beads_average_dist
GMM_dist_model = Mclust(hippocampus_region$hippocampus_dist, G = 3, verbose = T)
hippocampus_region$hippocampus_region_gmm = as.character(GMM_dist_model$classification)
hippocampus_region = hippocampus_region %>% filter(hippocampus_region_gmm != "3")

saveRDS(hippocampus_region, file = "./Hippocampus_region.RDS")

### cortex region
cortex_region = meta[which(meta$seurat_clusters %in% c("1", "2", "3")), ]
sp.loc = cortex_region[, c("new_x", "new_y")]
sp.loc$new_x = as.integer(sp.loc[["new_x"]]); sp.loc$new_y = as.integer(sp.loc[["new_y"]])
nnDist = FNN::get.knnx(sp.loc[, c("new_x", "new_y")], sp.loc[, c("new_x", "new_y")], k = 20)
loc.mat = cbind(nnDist$nn.index, round(nnDist$nn.dist, 2))
average_dist = function(x, loc.mat)
{ 
  return(mean(loc.mat[x, 22:dim(loc.mat)[2]]))
}
nrows = 1:nrow(loc.mat)
beads_average_dist = lapply(nrows, average_dist, loc.mat)
beads_average_dist = do.call(c, beads_average_dist)
cortex_region$cortex_dist = beads_average_dist
GMM_dist_model = Mclust(cortex_region$cortex_dist, G = 3, verbose = T)
cortex_region$cortex_region_gmm = as.character(GMM_dist_model$classification)
cortex_region = cortex_region %>% filter(cortex_region_gmm != "3")
saveRDS(cortex_region, file = "./Cortex_Region.RDS")

# for the velocyto
cortex_barcodes=cortex_region$degen_barcodes
cortex_barcodes=paste0("Hippocampus_Top5W_degen_V3CE0:", gsub("-1", "", cortex_barcodes))
write.csv(cortex_barcodes, file = "./Velocyto/cortex_barcodes_obs.csv", row.names=F)