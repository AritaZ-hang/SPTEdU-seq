rm(list = ls()); gc()
setwd("/your/workdir/")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(mclust)
library(Matrix)
library(RColorBrewer)
library(FNN)
library(igraph)
library(data.table)
library(spacexr)
library(Seurat)
library(scales)
library(BuenColors)
library(ggrepel)

if(!dir.exists("RCTD"))
  dir.create("RCTD")

mob = readRDS(file = "./Barcode_level_0109-01MOB.rds")
meta = readRDS(file = "./Barcode_level_meta_filter.rds")

mob.counts = GetAssayData(mob, slot = "counts")
coords.input = data.frame(new_x = -meta$new_x, 
                          new_y = -meta$new_y)
rownames(coords.input) = meta$degen_barcodes
nUMI = colSums(mob.counts)
use_beads = colnames(mob.counts)

ref = readRDS(file = "./scRNA_MOB_reference_RCTD.rds")

## create Spatial RNA object ##
puck = SpatialRNA(coords.input, mob.counts, nUMI)

## examine spatialRNA obejct ##
print(dim(puck@counts)) 
pdf(file = "./RCTD/puck_histogram_UMI.pdf", width = 6, height = 6)
hist(log(puck@nUMI, 2))
dev.off()

## running RCTD ##
myRCTD = create.RCTD(puck, ref, max_cores = 12)
myRCTD = run.RCTD(myRCTD, doublet_mode = "doublet")

## RCTD results ##
results = myRCTD@results
norm_weights = normalize_weights(results$weights)
cell_type_names = myRCTD@cell_type_info$info[[2]]
spatialRNA = myRCTD@spatialRNA
resultdir = "./RCTD/RCTD_Plots"
dir.create(resultdir)

## make the plots ##
plot_weights(cell_type_names, spatialRNA, resultdir, norm_weights)
plot_weights_unthreshold(cell_type_names, spatialRNA, resultdir, norm_weights)
plot_weights_doublet(cell_type_names, spatialRNA, resultdir, results$weights_doublet, results$results_df)
plot_cond_occur(cell_type_names, resultdir, norm_weights, spatialRNA)
plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultdir)

saveRDS(myRCTD, file = "./RCTD/mob_RCTD.rds")


#########
mob_deconv = data.frame(x_coords = myRCTD@spatialRNA@coords$x, 
                        y_coords = myRCTD@spatialRNA@coords$y, 
                        spot_class = myRCTD@results$results_df$spot_class, 
                        first_type = myRCTD@results$results_df$first_type,
                        second_type = myRCTD@results$results_df$second_type) 
rownames(mob_deconv) = rownames(myRCTD@results$results_df)

mob_deconv.filter = mob_deconv %>% filter(spot_class != "reject")

colors=c("#E2B8A6", 
         "#B9B3D8", 
         "#ffe8a8", "#ffe59c", "#f3cb8e", "#e7b280", "#da9871", "#ce7e63",
         "#ff8787", 
         "#9cd4c7", "#83caba", 
         "#bfdda3", "#b1d691", "#92c467", 
         "#EDD594", "#E1c379", 
         "#f5bdbc", "#e1a2a4", 
         "#bfe3ef", "#aedcf4")

ggplot(mob_deconv.filter, aes(x = -x_coords, y = -y_coords, color = first_type)) + geom_point(size = 0.3) + theme_minimal() +  xlab("Dim 1") + ylab("Dim 2") + coord_fixed() + theme(panel.grid = element_blank(), axis.text = element_blank(), legend.key.size = unit(6, "mm")) + annotate("segment", x = 6000, xend = (6000 + 200/0.72), y = 2500, yend = 2500, linewidth = 2) + guides(color=guide_legend(override.aes = list(size=5))) + labs(color = "") + scale_color_manual(values = colors)
ggsave(filename="./RCTD_first_type.pdf", width = 8, height = 6)

first_type_ratio = as.data.frame(table(mob_deconv.filter$first_type))
first_type_ratio = first_type_ratio %>% filter(first_type_ratio$Freq != 0)
first_type_ratio$Prop = round(first_type_ratio$Freq / sum(first_type_ratio$Freq)*100, 2)

first_type_ratio = first_type_ratio %>% 
  arrange(desc(Prop)) %>% 
  mutate(lab.ypos = cumsum(Prop) - 0.5*Prop)

ggplot(first_type_ratio, aes(x = "", y = Prop, fill = Var1)) + 
  geom_bar(width = 1, stat = "identity", color="white") + 
  coord_polar("y", start = 0) + theme_void() + 
  scale_fill_manual(values=colors) + 
  geom_text_repel(aes(y = lab.ypos, label=Prop), color="black", size=6) + labs(fill = "") + guides(fill = guide_legend(override.aes = list(size=4)))
ggsave(filename="./RCTD_first_type_ratio.pdf", width = 8, height=6)

### spot class num
spot_class_colors=c("#FFBE7A", "#FA7F6F", "#82B0D2")
spot_class_freq = as.data.frame(table(mob_deconv.filter$spot_class))
spot_class_freq = spot_class_freq %>% filter(spot_class_freq$Freq !=0)

ggplot(spot_class_freq, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + scale_fill_manual(values=spot_class_colors) + scale_color_manual(values=spot_class_colors) + theme_classic() + theme(legend.position = "none") + xlab("") + ylab("") + theme(axis.text = element_text(size=12, color="black"))
ggsave(filename="./RCTD_spot_class.pdf", width= 6, height=6)

saveRDS(mob_deconv, file = "./RCTD/mob_deconv_df.rds")
saveRDS(mob_deconv.filter, file = "./RCTD/mob_deconv.filter.rds")
