rm(list = ls())
setwd("/your/workdir/")

set.seed(1)
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(mclust)
library(factoextra)
library(cluster)
library(dbscan)
library(CytoTRACE2)
library(ggsignif)
library(ggpubr)
library(ggnewscale)

#
data.filter.final = readRDS(file="./Total_map/Total_map_meta.RDS")
data.filter.final = data.filter.final[order(-data.filter.final$ratio_adjusted),]


ratio_rank = data.filter.final$ratio_adjusted
ratio_rank_plot = data.frame(ratio_adjusted = ratio_rank, 
                             rank = 1:length(ratio_rank))
ratio_rank_plot$logRank=log10(ratio_rank_plot$rank)

kd = density(ratio_rank_plot$ratio_adjusted, bw = 0.1, kernel="gaussian", n = 10000, from = min(ratio_rank_plot$ratio_adjusted), to = max(ratio_rank_plot$ratio_adjusted))
descending = diff(c(Inf, kd$y)) < 0L
descending_indices = cumsum(rle(descending)$lengths)
min_threshold = kd$x[descending_indices[[2]]]

filter=data.filter.final %>% filter(data.filter.final$ratio_adjusted > min_threshold)

cluster=kmeans(filter[, c("ratio_adjusted")], 3)

filter$kmeans_cluster=as.character(cluster$cluster)

filter_kmeans_rank=data.frame(kmeans_cluster=filter$kmeans_cluster, 
                              ratio_adjusted=filter$ratio_adjusted)
filter_kmeans_rank=filter_kmeans_rank[order(filter_kmeans_rank$ratio_adjusted, decreasing=T),]
filter_kmeans_rank$rank = 1:dim(filter_kmeans_rank)[1]

mean_ratio_adjusted=aggregate(filter$ratio_adjusted, list(filter$kmeans_cluster), mean)
median_group = mean_ratio_adjusted$Group.1[which(mean_ratio_adjusted$x == median(mean_ratio_adjusted$x))]
threshold_kmeans = max(filter$ratio_adjusted[which(filter$kmeans_cluster == median_group)])
filter2=filter %>% filter(ratio_adjusted >= threshold_kmeans)

# gmm

GMM_dist_model = Mclust(filter2[, c("x_coords", "y_coords", "UTIs", "UMIs", "ratio_adjusted")], verbose = T, G=2)
filter2$gmm_class = as.character(GMM_dist_model$classification)
gmm_max_ratio = aggregate(filter2$ratio_adjusted, list(filter2$gmm_class), max)
gmm_min_ratio = aggregate(filter2$ratio_adjusted, list(filter2$gmm_class), min)
gmm_mean_ratio = aggregate(filter2$ratio_adjusted, list(filter2$gmm_class), mean)

select_group = gmm_mean_ratio$Group.1[which.max(gmm_mean_ratio$x)]
non_select_group = setdiff(gmm_mean_ratio$Group.1, select_group)

final_threshold=max(gmm_max_ratio$x[which(gmm_max_ratio$Group.1 == non_select_group)], gmm_min_ratio$x[which(gmm_min_ratio$Group.1 == select_group)]) # 0.7

# 
filter3 = filter2 %>% filter(gmm_class == select_group & ratio_adjusted >= final_threshold) 

saveRDS(filter3, file = "./Total_map/UTI_threshold_filtered.RDS")
