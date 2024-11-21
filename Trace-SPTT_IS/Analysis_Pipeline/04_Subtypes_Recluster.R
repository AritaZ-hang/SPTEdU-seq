# for example, for astrocytes


rm(list = ls())
setwd("/your/workdir/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(mclust)
library(tidyverse)
library(harmony)
library(reshape2)
library(ggalluvial)
library(ggh4x)
library(patchwork)
library(BuenColors)
library(scales)
library(CytoTRACE2)

seob = readRDS(file="./Astro_raw.RDS")
#
seob = SCTransform(seob, assay = "RNA", ncells = 3000, verbose=T, vars.to.regress = c("nCount_RNA"))
seob = RunPCA(seob)
ElbowPlot(seob, ndims = 50)

seob = RunUMAP(seob, dims = 1:20)
seob = FindNeighbors(seob, dims = 1:20)
seob = FindClusters(seob, resolution = 0.3, verbose=T)

markers=FindAllMarkers(seob, logfc.threshold = 0.1, verbose = T, only.pos = T, min.pct = 0.1)
markers.filtered = markers %>% filter(p_val_adj <= 0.05)

write.csv(markers.filtered, "./Markers_filtered.csv", row.names=F)

seob$annotations=""
seob$annotations[which(seob$seurat_clusters=="0")] = "Transit astrocytes"
seob$annotations[which(seob$seurat_clusters=="1")] = "Excitatory neurons"
seob$annotations[which(seob$seurat_clusters=="2")] = "Inhibitory neurons"
seob$annotations[which(seob$seurat_clusters=="3")] = "Reactive astrocyte"
seob$annotations[which(seob$seurat_clusters=="4")] = "Fibroblasts"
seob$annotations[which(seob$seurat_clusters=="5")] = "Restorative astrocyte"
seob$annotations[which(seob$seurat_clusters=="6")] = "Non-reactive astrocyte"
seob$annotations[which(seob$seurat_clusters=="7")] = "Macrophage"


seob$celltype=""
seob$celltype[which(seob$seurat_clusters=="0")] = "Transit astrocytes"
seob$celltype[which(seob$seurat_clusters=="1")] = "Excitatory neurons"
seob$celltype[which(seob$seurat_clusters=="2")] = "Inhibitory neurons"
seob$celltype[which(seob$seurat_clusters=="3")] = "Reactive astrocyte"
seob$celltype[which(seob$seurat_clusters=="4")] = "Fibroblasts"
seob$celltype[which(seob$seurat_clusters=="5")] = "Restorative astrocyte"
seob$celltype[which(seob$seurat_clusters=="6")] = "Non-reactive astrocyte"
seob$celltype[which(seob$seurat_clusters=="7")] = "Macrophage"


# cluster_colors, 8
cluster_colors = c("#bd338f", "#eb8252", "#f5dc83", "#cdd4dc", "#8098A2", "#8FA33f", "#5f7929", "#014820")

DimPlot(seob, reduction = "umap", label = F, cols = cluster_colors, raster=F, group.by = "annotations") + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + ggtitle("") + guides(color=guide_legend(override.aes = list(size=5)))
ggsave(filename="./Astrocytes_subtype_clusters.pdf", width = 7, height = 6)


###
stage_colors = c("#ff847c", "#f03861", "#45171d")
#stage_colors = c("#fecea8", "#ff847c", "#f03861", "#45171d") uninjured, 3dpi, 7dpi, 14dpi
seob$stage = factor(seob$stage, levels = c("3dpi", "7dpi", "14dpi"))

DimPlot(seob, reduction = "umap", label = F, cols = stage_colors, raster=F, group.by = "stage") + theme_minimal() + xlab("UMAP 1") + ylab("UMAP 2") + theme(axis.text = element_blank(), panel.grid = element_blank()) + ggtitle("") + guides(color=guide_legend(override.aes = list(size=5)))
ggsave(filename="./Astrocytes_subtype_stages.pdf", width = 7, height = 6)

# 
stage_annotations=as.data.frame(table(seob$stage, seob$annotations))
stage_annotations$ratio = 0
for(i in stage_annotations$Var1)
{
  stage_annotations$ratio[which(stage_annotations$Var1 ==i)] = stage_annotations$Freq[which(stage_annotations$Var1==i)] / sum(stage_annotations$Freq[which(stage_annotations$Var1 == i)])
}
stage_annotations$Var1 = factor(stage_annotations$Var1, levels = c("uninjured", "3dpi", "7dpi", "14dpi"))

ggplot(stage_annotations, aes(x = Var1, y = ratio*100, fill = Var2, alluvium = Var2, stratum = Var2)) + 
  geom_col(position = "stack", width = 0.5) + 
  geom_alluvium(aes(stratum = Var2), width = 0.5, alpha = 0.4, color = "white", linewidth = 0.8, curve_type="linear") + 
  geom_stratum(width = 0.5, alpha=0.8, color="white") + 
  scale_fill_manual(values=cluster_colors) + 
  labs(x="", y = "Ratio%", title="", fill="Annotation") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle=45, hjust = 1))
ggsave(filename="./Astrocytes_subtype_stage_annotations_trend.pdf", width = 6, height = 7)

###
annotations_stage = as.data.frame(table(seob$annotations, seob$stage))
annotations_stage$ratio = 0
for(i in annotations_stage$Var1)
{
  annotations_stage$ratio[which(annotations_stage$Var1 == i)] = annotations_stage$Freq[which(annotations_stage$Var1 == i)] / sum(annotations_stage$Freq[which(annotations_stage$Var1 == i)])
}
annotations_stage$Var2=factor(annotations_stage$Var2, levels = c("uninjured", "3dpi", "7dpi", "14dpi"))
ggplot(annotations_stage, aes(x = Var1, y = ratio*100, fill = Var2)) + geom_bar(stat = "identity") + theme_minimal() + scale_fill_manual(values = stage_colors) + xlab("") + ylab("Ratio%") + labs(fill = "Stage") + theme(axis.text.x = element_text(angle=45, hjust = 1))
ggsave(filename='./Astrocytes_subtype_annotations_stage_trend.pdf', width = 10, height = 7)

# 
genes = c("Igfbp5", "Meg3", "Slc1a2", "Jarid2", "Camkmt", "Gfap")
p1 = DotPlot(seob, features=genes, group.by = "stage")
genes_df = p1$data
genes_df$id = factor(genes_df$id, levels = c("3dpi", "7dpi", "14dpi"))
ggplot(genes_df, aes(x = id, y = avg.exp , color = features.plot, group = features.plot)) + geom_point(size=3) + theme_classic() + xlab("Stage") + ylab("Average Expression") + geom_line(linewidth = 1, linetype = "dashed") + scale_color_manual(values=c("#FAB9AC", "#7BBC53", "#DE6736", "#3f474b", "#E6B90D")) + labs(color = "Genes") + theme(axis.text = element_text(size=12), axis.title=element_text(size=16))
ggsave(filename="./Astrocytes_genes_Average_Expression.pdf", width = 6, height = 4)

#
colors_map = data.frame(colors=cluster_colors, 
                        celltype=unique(seob$annotations)[order(unique(seob$annotations))])

counts=GetAssayData(seob, assay="SCT", slot="counts")
tmp = t(as.matrix(counts[c(genes),]))
tmp.df = cbind(seob@meta.data, tmp)
tmp.df$stage = factor(tmp.df$stage, levels = c("uninjured", "3dpi", "7dpi", "14dpi"))
tmp.df$timepoint=0
tmp.df$timepoint[which(tmp.df$stage=="uninjured")] = 0
tmp.df$timepoint[which(tmp.df$stage=="3dpi")] = 3
tmp.df$timepoint[which(tmp.df$stage=="7dpi")] = 7
tmp.df$timepoint[which(tmp.df$stage=="14dpi")]=14

# inspired by STRS
for(FEAT in genes)
{
  ggplot(data=tmp.df, 
         mapping = aes_string(x = "timepoint", 
                              y = FEAT, 
                              color="stage")) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Excitatory neurons",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Excitatory neurons")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Excitatory neurons",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Excitatory neurons")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Excitatory neurons",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Excitatory neurons")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Fibroblasts",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Fibroblasts")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Fibroblasts",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Fibroblasts")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Fibroblasts",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Fibroblasts")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Inhibitory neurons",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Inhibitory neurons")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Inhibitory neurons",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Inhibitory neurons")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Inhibitory neurons",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Inhibitory neurons")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Macrophage",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Macrophage")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Macrophage",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Macrophage")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Macrophage",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Macrophage")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Non-reactive astrocyte",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Non-reactive astrocyte")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Non-reactive astrocyte",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Non-reactive astrocyte")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Non-reactive astrocyte",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Non-reactive astrocyte")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Reactive astrocyte",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Reactive astrocyte")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Reactive astrocyte",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Reactive astrocyte")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Reactive astrocyte",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Reactive astrocyte")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) +
    
    stat_summary(data = tmp.df[tmp.df$celltype == "Restorative astrocyte",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Restorative astrocyte")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Restorative astrocyte",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Restorative astrocyte")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Restorative astrocyte",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Restorative astrocyte")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) +
    stat_summary(data = tmp.df[tmp.df$celltype == "Transit astrocytes",], 
                 fun = mean, 
                 color = colors_map$colors[which(colors_map$celltype == "Transit astrocytes")], 
                 geom="line", 
                 size=1) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Transit astrocytes",], 
                 fun = mean,
                 color = colors_map$colors[which(colors_map$celltype == "Transit astrocytes")], 
                 geom="point", 
                 size=3) + 
    stat_summary(data = tmp.df[tmp.df$celltype == "Transit astrocytes",], 
                 fun.data=mean_sdl, 
                 color=colors_map$colors[which(colors_map$celltype == "Transit astrocytes")], 
                 geom = "errorbar", 
                 alpha=0.7, width = 0.5) +
    
    theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(face="italic")) + guides(color=guide_legend(override.aes = list(size=2))) + scale_x_continuous(breaks = c(0, 3, 7, 14))
  ggsave(filename=paste0("./", FEAT, "_lineplot.pdf"), width = 4, height=2)
}

saveRDS(seob, file = "./Astrocyte_anno.RDS")