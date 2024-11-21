rm(list = ls())
setwd("/your/workdir/")
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(mclust)

dge_summary = as.data.frame(fread(file = "./total_dge/allmapped_dge.summary.txt")); colnames(dge_summary) = c("barcodes", "genic_reads", "UMIs", "genes")
coords = as.data.frame(fread(file = "./coords/ST110098_A1_FilterBarcodes.csv")); colnames(coords) = c("x_coords", "y_coords", "barcodes")
utis = as.data.frame(fread(file = "./utis/utis.csv", sep=",")); colnames(utis) = c("barcodes", "UTIs")

intersect_beads = intersect(dge_summary$barcodes, coords$barcodes) %>% intersect(utis$barcodes) 

dge_summary.sub = dge_summary[which(dge_summary$barcodes %in% intersect_beads), ]
coords.sub = coords[which(coords$barcodes %in% intersect_beads),]
utis.sub = utis[which(utis$barcodes %in% intersect_beads),]

data = inner_join(dge_summary.sub, coords.sub, by = "barcodes") %>% inner_join(utis.sub)
data = data %>% filter(UMIs>0)

data$ratio = data$UTIs / data$UMIs
summary(data$ratio)
summary(data$UMIs) 
data$ratio_adjusted = log1p(data$ratio)
summary(data$ratio_adjusted) 
data$log2ratio=log2(data$ratio)
summary(data$log2ratio) 

# data UMI Threshold
rank_UMI=data.frame(UMIs=data$UMIs[order(data$UMIs, decreasing = T)])
rank_UMI$rank=1:dim(rank_UMI)[1]

GMM_umi_model=Mclust(log10(rank_UMI$UMIs), G=3, verbose=T)
rank_UMI$umi_gmm = as.character(GMM_umi_model$classification)
umi_threshold=max(rank_UMI$UMIs[which(rank_UMI$umi_gmm=="1")])
data.filter=data%>% filter(UMIs>=umi_threshold)

# dist GMM

sp.loc = data.filter[, c("x_coords", "y_coords")]
sp.loc$new_x = as.integer(sp.loc[["x_coords"]]); 
sp.loc$new_y = as.integer(sp.loc[["y_coords"]])
nnDist = FNN::get.knnx(sp.loc[, c("x_coords", "y_coords")], sp.loc[, c("x_coords", "y_coords")], k = 20)
loc.mat = cbind(nnDist$nn.index, round(nnDist$nn.dist, 2))
average_dist = function(x, loc.mat)
{
  return(mean(loc.mat[x, 22:dim(loc.mat)[2]]))
}
nrows = 1:nrow(loc.mat)
beads_average_dist = lapply(nrows, average_dist, loc.mat)
beads_average_dist = do.call(c, beads_average_dist)
data.filter$average_dist = beads_average_dist
GMM_dist_model = Mclust(data.filter$average_dist, G = 4, verbose = T)
data.filter$dist_gmm = as.character(GMM_dist_model$classification)
data.filter.final=data.filter %>% filter(dist_gmm !="4")

saveRDS(data.filter.final, file="./Total_map/Total_map_meta.RDS")
save(data, data.filter, data.filter.final, file = "./Total_map/Total_map_meta_all.RData")
