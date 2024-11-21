rm(list = ls())
setwd("/your/workdir/")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(mclust)
library(Matrix)
library(RColorBrewer)
library(FNN)
library(Seurat)
library(data.table)
library(spacexr)
library(cowplot)
library(BuenColors)

myRCTD = readRDS(file = "./RCTD/hippocampus_RCTD.rds")
hippocampus_deconv = readRDS(file = "./RCTD/hippocampus_deconv_df.rds")
all_colors_use = readRDS(file = "./RCTD/hippocampus_all_colors_use.rds")

# hippocampus region

hippocampus_region  = readRDS(file = "./Hippocampus_region.RDS")
hippocampus_deconv.selected = hippocampus_deconv[hippocampus_region$degen_barcodes, ]

# embedding
hippocampus_deconv.selected.filtered = hippocampus_deconv.selected %>% filter(spot_class != "reject")
celltype_colors = all_colors_use$colors[which(all_colors_use$celltype %in% hippocampus_deconv.selected.filtered$first_type)]

pdf(file = "./RCTD/hippocampus_first_type_CAonly.pdf", width = 8, height = 8)
ggplot(hippocampus_deconv.selected.filtered, aes(x = x_coords, y = y_coords, color = first_type)) + geom_point(size = 0.4) + theme_minimal() + scale_color_manual(values = celltype_colors) + xlab("Dim 1") + ylab("Dim 2") + coord_fixed() + theme(panel.grid = element_blank()) + annotate("segment", x = -7500, xend = -(7500 + 200/0.72), y = -4500, yend = -4500, linewidth = 2) + theme(axis.text = element_blank()) + theme(legend.position = "right") + labs(color = "") + guides(color = guide_legend(override.aes = list(size=5)))
dev.off()

########### first type ratio ##########
first_type_ratio = as.data.frame(table(hippocampus_deconv.selected.filtered$first_type))
first_type_ratio = first_type_ratio %>% filter(first_type_ratio$Freq !=0)
first_type_ratio$Prop = round(first_type_ratio$Freq / sum(first_type_ratio$Freq)*100, 2)

first_type_ratio = first_type_ratio %>% 
  arrange(desc(Prop)) %>% 
  mutate(lab.ypos = cumsum(Prop) - 0.5*Prop)

ggplot(first_type_ratio, aes(x = "", y = Prop, fill = Var1)) + 
  geom_bar(width = 1, stat = "identity", color="white") + 
  coord_polar("y", start = 0) + theme_void() + 
  scale_fill_manual(values=celltype_colors) + 
  labs(fill = "") + 
  guides(fill = guide_legend(override.aes = list(size=4)))
ggsave(filename="./RCTD/RCTD_first_type_ratio_CAonly.pdf", width = 8, height=6)

############# spot class ratio #############
spot_class.freq = as.data.frame(table(hippocampus_deconv.selected.filtered$spot_class))
spot_class.freq = spot_class.freq %>% filter(spot_class.freq$Freq !=0)
spot_class_colors = c("#FFBE7A", "#FA7F6F", "#82B0D2")

ggplot(spot_class.freq, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + scale_fill_manual(values = spot_class_colors) + scale_color_manual(values = spot_class_colors) + theme_classic() + theme(legend.position = "none") + xlab("") + ylab("") + theme(axis.text = element_text(size=12, color="black"))
ggsave(filename="./RCTD/RCTD_spot_class_CAonly.pdf", width = 6, height = 6)


########### seurat clusters vs RCTD ##############
hippocampus_region.filter = hippocampus_region[rownames(hippocampus_deconv.selected.filtered),]
hippocampus_region.filter$RCTD_first_type = hippocampus_deconv.selected.filtered$first_type

predictions = table(hippocampus_region.filter$cell_type, hippocampus_region.filter$RCTD_first_type)
predictions = predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions = as.data.frame(predictions)

fraction_colors = jdb_palette("brewer_violet", type = "continuous")

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile()  + xlab("Manual Annotated Cell Types") + ylab("RCTD Deconvolution Cell Types") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_gradientn(colors=fraction_colors) + labs(fill = "Fraction of Barcodes")
p1
ggsave(p1, filename = "./RCTD/hippocampus_seurat_clusters_RCTD.pdf", width = 8, height = 6)
