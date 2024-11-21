rm(list = ls()); gc()
setwd("/media/ggj/Guo-4T-G/ST_benchmark/Barcode_UMI/0222-01Hippocampus/")
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
library(cowplot)

# make RCTD scRNA-seq reference
ref = readRDS(file = "./scRNA_hippocampus_reference.RDS")

ref_counts = ref@assays$RNA@counts
rownames(ref_counts)
ref_meta = ref@meta.data
ref_cell_types = ref$celltype
names(ref_cell_types) = rownames(ref_meta)
ref_cell_types = as.factor(ref_cell_types)
ref_nUMI = ref@meta.data$nCount_RNA
names(ref_nUMI) = rownames(ref_meta)

#### create the reference object
ref_RCTD = Reference(ref_counts, ref_cell_types, ref_nUMI)
saveRDS(ref_RCTD, file = "./scRNA_hippocampus_reference_RCTD.rds")

rm(list=ls())
gc()

# Deconvolution by RCTD

if(!dir.exists("RCTD"))
  dir.create("RCTD")

hippocampus = readRDS(file = "./Barcode_Level_0222-01Hippocampus.rds")
meta = readRDS(file = "./Barcode_level_meta_filter.rds")

hippocampus.counts = GetAssayData(hippocampus, slot = "counts")
coords.input = data.frame(new_x = -meta$new_x, 
                          new_y = -meta$new_y)
rownames(coords.input) = meta$degen_barcodes
nUMI = colSums(hippocampus.counts)
use_beads = colnames(hippocampus.counts)

ref = readRDS(file = "./scRNA_hippocampus_reference_RCTD.rds")

## create Spatial RNA object ##
puck = SpatialRNA(coords.input, hippocampus.counts, nUMI)

## examine spatialRNA obejct ##
print(dim(puck@counts))
pdf(file = "./RCTD/puck_histogram_UMI.pdf", width = 6, height = 6)
hist(log(puck@nUMI, 2))
dev.off()

## on the plot ##
pdf(file = "./RCTD/plot_puck_continuous.pdf", width = 6, height = 6)
plot_puck_continuous(puck, use_beads, puck@nUMI, ylimit = c(0, round(quantile(puck@nUMI, 0.9))), title = "plot of nUMI")
dev.off()

## running RCTD ##
myRCTD = create.RCTD(puck, ref, max_cores = 5)
myRCTD = run.RCTD(myRCTD, doublet_mode = "doublet")

## RCTD results ##
results = myRCTD@results
norm_weights = normalize_weights(results$weights)
cell_type_names = myRCTD@cell_type_info$info[[2]]
spatialRNA = myRCTD@spatialRNA
resultdir = "./RCTD/RCTD_Plots"
dir.create(resultdir)

saveRDS(myRCTD, file = "./RCTD/hippocampus_RCTD.rds")


## plot
hippocampus_deconv = data.frame(x_coords = myRCTD@spatialRNA@coords$x, 
                                y_coords = myRCTD@spatialRNA@coords$y, 
                                spot_class = myRCTD@results$results_df$spot_class, 
                                first_type = myRCTD@results$results_df$first_type, 
                                second_type = myRCTD@results$results_df$second_type)
rownames(hippocampus_deconv) = rownames(myRCTD@results$results_df)

hippocampus_deconv.clear = hippocampus_deconv %>% filter(spot_class != "reject")

colors=c("#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#6fb3a8", "#b3e19b", "#50aa4b", "#ff9d9f", "#f36569", '#3581b7', "#cdb6da", "#704ba3", "#9a7fbd", "#dba9a8", "#e43030")

all_colors_use = data.frame(celltype=c("Astrocyte", "Endothelial", "Ependyma", "Fibroblast-Like", "Interneuron", "Microglia_Macrophage", "Mural", "Neurogenesis", "Neuron_CA1", "Neuron_CA2CA3", "Neuron_Dentate", "Neuron_Subiculum", "Neuron_Subiculumn_Entorhinal", "Oligodendrocyte", "Polydendrocyte"), 
                            colors = colors)

hippocampus_deconv.clear$first_type = as.character(hippocampus_deconv.clear$first_type)

pdf(file = "./RCTD/hippocampus_first_type.pdf", width = 8, height = 8)
ggplot(hippocampus_deconv.clear, aes(x = x_coords, y = y_coords, color = first_type)) + geom_point(size = 0.3) + theme_minimal() + xlab("Dim 1") + ylab("Dim 2") + theme(axis.text = element_blank()) + coord_fixed() + theme(panel.grid = element_blank(), axis.text = element_blank()) + annotate("segment", x = -8000, xend = -(8000 + 200/0.72), y = -2500, yend = -2500, linewidth = 2) + guides(color=guide_legend(override.aes = list(size=5))) + labs(color = "") + scale_color_manual(values=all_colors_use$colors[which(all_colors_use$celltype %in% unique(hippocampus_deconv.clear$first_type))])
dev.off()

########## spot class ###########
spot_class.freq = as.data.frame(table(hippocampus_deconv.clear$spot_class))
spot_class.freq = spot_class.freq %>% filter(spot_class.freq$Freq !=0)
spot_class_colors = c("#FFBE7A", "#FA7F6F", "#82B0D2")

ggplot(spot_class.freq, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + scale_fill_manual(values=spot_class_colors) + scale_color_manual(values=spot_class_colors) + theme_classic() + theme(legend.position = "none") + xlab("") + ylab("") + theme(axis.text = element_text(size=12, color="black"))
ggsave(filename="./RCTD/RCTD_spot_class.pdf", width= 6, height=6)

saveRDS(hippocampus_deconv, file = "./RCTD/hippocampus_deconv_df.rds")
saveRDS(all_colors_use, file = "./RCTD/hippocampus_all_colors_use.rds")

########### first ############
first_type_ratio = as.data.frame(table(hippocampus_deconv.clear$first_type))
first_type_ratio = first_type_ratio %>% filter(first_type_ratio$Freq !=0)
first_type_ratio$Prop = round(first_type_ratio$Freq / sum(first_type_ratio$Freq)*100, 2)

first_type_ratio = first_type_ratio %>% 
  arrange(desc(Prop)) %>% 
  mutate(lab.ypos = cumsum(Prop) - 0.5*Prop)

ggplot(first_type_ratio, aes(x = "", y = Prop, fill = Var1)) + 
  geom_bar(width = 1, stat = "identity", color="white") + 
  coord_polar("y", start = 0) + theme_void() + 
  scale_fill_manual(values=colors) + 
  labs(fill = "") + 
  guides(fill = guide_legend(override.aes = list(size=4)))
ggsave(filename="./RCTD_first_type_ratio.pdf", width = 8, height=6)