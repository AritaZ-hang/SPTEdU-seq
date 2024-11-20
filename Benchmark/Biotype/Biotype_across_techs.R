rm(list = ls()); gc()
setwd("/workdir/")
library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

sptt_seq = as.data.frame(fread(file = "SPTT_genes.txt", header = F, sep = ""))
visium = as.data.frame(fread(file = "Visium_genes.txt", header = F, sep = ""))
stereo = as.data.frame(fread(file = "Stereo_genes.txt", header = F, sep = ""))
decoder = as.data.frame(fread(file = "Decoder_genes.txt", header = F, sep = ""))
dbit = as.data.frame(fread(file = "Dbit_genes.txt", header = F, sep = ""))
strs = as.data.frame(fread(file = "STRS_genes.txt", header = F, sep = ""))
slide = as.data.frame(fread(file = "Slide_genes.txt", header = F, sep  =""))

sptt_seq$V1 = gsub("gn:Z:", "", sptt_seq$V1)
visium$V1 = gsub("gn:Z:", "", visium$V1)
stereo$V1 = gsub("gn:Z:", "", stereo$V1)
decoder$V1 = gsub("gn:Z:", "", decoder$V1)
dbit$V1 = gsub("gn:Z:", "", dbit$V1)
strs$V1 = gsub("gn:Z:", "", strs$V1)
slide$V1 = gsub("gn:Z:", "", slide$V1)

# biotype
mouse_use_biotype = readRDS("./mouse_use_biotype.rds")

colnames(sptt_seq) = colnames(visium) = colnames(stereo) = colnames(decoder) = colnames(dbit) = colnames(strs) = colnames(slide)= "gene"

# 

# inner join
sptt.full = left_join(sptt_seq, mouse_use_biotype, by = "gene")
visium.full = left_join(visium, mouse_use_biotype, by = "gene")
stereo.full = left_join(stereo, mouse_use_biotype, by = "gene")
decoder.full = left_join(decoder, mouse_use_biotype, by = "gene")
dbit.full = left_join(dbit, mouse_use_biotype, by = "gene")
strs.full = left_join(strs, mouse_use_biotype, by = "gene")
slide.full = left_join(slide, mouse_use_biotype, by = "gene")

sptt.full[grep(",", sptt.full$gene),]$biotype = "Multi"
visium.full[grep(",", visium.full$gene),]$biotype = "Multi"
stereo.full[grep(",", stereo.full$gene),]$biotype = "Multi"
decoder.full[grep(",", decoder.full$gene),]$biotype = "Multi"
dbit.full[grep(",", dbit.full$gene),]$biotype = "Multi"
strs.full[grep(",", strs.full$gene),]$biotype = "Multi"
slide.full[grep(",", slide.full$gene),]$biotype = "Multi"


# Gm42418, AY036118, Rn18s-rs5
sptt.full[grep("Gm42418", sptt.full$gene),]$biotype = "rRNA"; sptt.full[grep("AY036118", sptt.full$gene),]$biotype = "rRNA"; sptt.full[grep("Rn18s-rs5", sptt.full$gene),]$biotype = "rRNA"
visium.full[grep("Gm42418", visium.full$gene),]$biotype = "rRNA"; visium.full[grep("AY036118", visium.full$gene),]$biotype = "rRNA"; visium.full[grep("Rn18s-rs5", visium.full$gene),]$biotype = "rRNA"
stereo.full[grep("Gm42418", stereo.full$gene),]$biotype = "rRNA"; stereo.full[grep("AY036118", stereo.full$gene),]$biotype = "rRNA"; stereo.full[grep("Rn18s-rs5", stereo.full$gene),]$biotype = "rRNA"
decoder.full[grep("Gm42418", decoder.full$gene),]$biotype = "rRNA"; decoder.full[grep("AY036118", decoder.full$gene),]$biotype = "rRNA"; decoder.full[grep("Rn18s-rs5", decoder.full$gene),]$biotype = "rRNA"
dbit.full[grep("Gm42418", dbit.full$gene),]$biotype = "rRNA"; dbit.full[grep("AY036118", dbit.full$gene),]$biotype = "rRNA"; dbit.full[grep("Rn18s-rs5", dbit.full$gene),]$biotype = "rRNA"
strs.full[grep("Gm42418", strs.full$gene),]$biotype = "rRNA"; strs.full[grep("AY036118", strs.full$gene),]$biotype = "rRNA"; strs.full[grep("Rn18s-rs5", strs.full$gene),]$biotype = "rRNA"
slide.full[grep("Gm42418", slide.full$gene),]$biotype = "rRNA"; slide.full[grep("AY036118", slide.full$gene),]$biotype = "rRNA"; slide.full[grep("Rn18s-rs5", slide.full$gene),]$biotype = "rRNA"


deal_multi = function(full, biotype)
{

  multi = full[which(full$biotype == "Multi"),]
  full = full %>% filter(biotype != "Multi")
  multi_genes = stringr::str_split(multi$gene, pattern = ",")
  multi_target = (lapply(multi_genes, FUN = function(x){
    t = length(unique(x))
    if(t == 1)
      return(x)
  }))
  
  for(i in 1:dim(multi)[1])
  {
    if(!is.null(multi_target[[i]]))
    {
      target_biotype=biotype$biotype[which(biotype$gene == unique(multi_target[[i]]))]
      multi[i, "biotype"] = target_biotype[1]
    }
  }
  full = rbind(full, multi)
  
}

sptt.full = deal_multi(sptt.full, mouse_use_biotype)
visium.full = deal_multi(visium.full, mouse_use_biotype)
stereo.full = deal_multi(stereo.full, mouse_use_biotype)
decoder.full = deal_multi(decoder.full, mouse_use_biotype)
dbit.full = deal_multi(dbit.full, mouse_use_biotype)
strs.full = deal_multi(strs.full, mouse_use_biotype)
slide.full = deal_multi(slide.full, mouse_use_biotype)

table(sptt.full$biotype)
table(visium.full$biotype)
table(stereo.full$biotype)
table(decoder.full$biotype)
table(dbit.full$biotype)
table(strs.full$biotype)
table(slide.full$biotype)

# rRNA & Mt-rRNA

p2.data = data.frame()
p2.data = rbind(c(sum(sptt.full$biotype %in% "rRNA"), sum(sptt.full$biotype %in% "Mt_rRNA")), 
                c(sum(visium.full$biotype %in% "rRNA"), sum(visium.full$biotype %in% "Mt_rRNA")), 
                c(sum(stereo.full$biotype %in% "rRNA"), sum(stereo.full$biotype %in% "Mt_rRNA")), 
                c(sum(decoder.full$biotype %in% "rRNA"), sum(decoder.full$biotype %in% "Mt_rRNA")), 
                c(sum(dbit.full$biotype %in% "rRNA"), sum(dbit.full$biotype %in% "Mt_rRNA")), 
                c(sum(strs.full$biotype %in% "rRNA"), sum(strs.full$biotype %in% "Mt_rRNA")), 
                c(sum(slide.full$biotype %in% "rRNA"), sum(slide.full$biotype %in% "Mt_rRNA")))
p2.data = as.data.frame(p2.data)
colnames(p2.data) = c("rRNA", "Mt_rRNA")

p2.data[1,] = p2.data[1,] / dim(sptt.full)[1]
p2.data[2,] = p2.data[2,] / dim(visium.full)[1]
p2.data[3,] = p2.data[3,] / dim(stereo.full)[1]
p2.data[4,] = p2.data[4,] / dim(decoder.full)[1]
p2.data[5,] = p2.data[5,] / dim(dbit.full)[1]
p2.data[6,] = p2.data[6,] / dim(strs.full)[1]
p2.data[7,] = p2.data[7,] / dim(slide.full)[1]

p2.data$expr = c("SPTT_seq", "10X_Visium", "Stereo_seq", "Decoder_seq", "Dbit_seq", "STRS", "Slide-seqV2")
rRNA_color = c("#BEB8DC", "#E7DAD2")
rRNA_color2 = c("#F1B656", "#397FC7")

p2.melt = reshape2::melt(p2.data)
p2.melt$expr = factor(p2.melt$expr, levels = c("SPTT_seq", "10X_Visium", "Stereo_seq", "Decoder_seq", "Dbit_seq", "STRS", "Slide-seqV2"))
p2 = ggplot(p2.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') + scale_fill_manual(values = rRNA_color2) + theme_classic() + xlab('') + ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_text(color = 'black', size = 9))
p2
ggsave(p2, filename = "./rRNA_mt_rRNA_ratio.pdf", width = 8, height = 4)


################ p1
unique(mouse_use_biotype$biotype)
mouse_pseudogenes_types = unique(mouse_use_biotype$biotype)[grep("pseudogene", unique(mouse_use_biotype$biotype))]
mouse_data.matrix = data.frame()
mouse_target = c("protein_coding", "lincRNA", "misc_RNA", "rRNA", "Mt_rRNA", "TEC", mouse_pseudogenes_types)
mouse_rest = setdiff(unique(mouse_use_biotype$biotype), mouse_target)

data.matrix = rbind(c(sum(sptt.full$biotype %in% "protein_coding"), sum(sptt.full$biotype %in% "lincRNA"), sum(sptt.full$biotype %in% "misc_RNA"), sum(sptt.full$biotype %in% mouse_pseudogenes_types), sum(sptt.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(sptt.full$biotype %in% "TEC"), sum(sptt.full$biotype %in% "Multi"), sum(sptt.full$biotype %in% mouse_rest)), 
                    c(sum(visium.full$biotype %in% "protein_coding"), sum(visium.full$biotype %in% "lincRNA"), sum(visium.full$biotype %in% "misc_RNA"), sum(visium.full$biotype %in% mouse_pseudogenes_types), sum(visium.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(visium.full$biotype %in% "TEC"), sum(visium.full$biotype %in% "Multi"), sum(visium.full$biotype %in% mouse_rest)),
                    c(sum(stereo.full$biotype %in% "protein_coding"), sum(stereo.full$biotype %in% "lincRNA"), sum(stereo.full$biotype %in% "misc_RNA"), sum(stereo.full$biotype %in% mouse_pseudogenes_types), sum(stereo.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(stereo.full$biotype %in% "TEC"), sum(stereo.full$biotype %in% "Multi"), sum(stereo.full$biotype %in% mouse_rest)), 
                    c(sum(decoder.full$biotype %in% "protein_coding"), sum(decoder.full$biotype %in% "lincRNA"), sum(decoder.full$biotype %in% "misc_RNA"), sum(decoder.full$biotype %in% mouse_pseudogenes_types), sum(decoder.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(decoder.full$biotype %in% "TEC"), sum(decoder.full$biotype %in% "Multi"), sum(decoder.full$biotype %in% mouse_rest)), 
                    c(sum(dbit.full$biotype %in% "protein_coding"), sum(dbit.full$biotype %in% "lincRNA"), sum(dbit.full$biotype %in% "misc_RNA"), sum(dbit.full$biotype %in% mouse_pseudogenes_types), sum(dbit.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(dbit.full$biotype %in% "TEC"), sum(dbit.full$biotype %in% "Multi"), sum(dbit.full$biotype %in% mouse_rest)), 
                    c(sum(strs.full$biotype %in% "protein_coding"), sum(strs.full$biotype %in% "lincRNA"), sum(strs.full$biotype %in% "misc_RNA"), sum(strs.full$biotype %in% mouse_pseudogenes_types), sum(strs.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(strs.full$biotype %in% "TEC"), sum(strs.full$biotype %in% "Multi"), sum(strs.full$biotype %in% mouse_rest)), 
                    c(sum(slide.full$biotype %in% "protein_coding"), sum(slide.full$biotype %in% "lincRNA"), sum(slide.full$biotype %in% "misc_RNA"), sum(slide.full$biotype %in% mouse_pseudogenes_types), sum(slide.full$biotype %in% c("rRNA", "Mt_rRNA")), sum(slide.full$biotype %in% "TEC"), sum(slide.full$biotype %in% "Multi"), sum(slide.full$biotype %in% mouse_rest)))

data.matrix = as.data.frame(data.matrix)
colnames(data.matrix) = c("Protein coding", "lncRNA/lincRNA", "miscRNA", "Pseudogenes", "rRNA & MT-rRNA", "TEC", "Multi Annotated", "Other")
data.matrix = data.matrix / rowSums(data.matrix)
data.matrix$expr = c("SPTT_seq", "10X_Visium", "Stereo_seq", "Decoder_seq", "Dbit_seq", "STRS", "Slide-seqV2")
data.melt = reshape2::melt(data.matrix)

biotype_colors = colorRampPalette(c("#FC757B", "#F88455", "#FDCA93", "#FFE59B", "#76CBB4", "#3C9BC9"))(8)
data.melt$expr = factor(data.melt$expr, levels = c("SPTT_seq", "10X_Visium", "Stereo_seq", "Decoder_seq", "Dbit_seq", "STRS", "Slide-seqV2"))

p1 = ggplot(data.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') + scale_fill_manual(values = biotype_colors) + theme_classic() + xlab('') + ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_text(color = 'black', size = 9))

p1

ggsave(p1, filename = "./biotype_total.pdf", width = 8, height = 4)

#
p3.data = data.frame()
p3.data = rbind(c(sum(sptt.full$biotype == "snoRNA"), sum(sptt.full$biotype == "miRNA"), sum(sptt.full$biotype == "scaRNA"), sum(sptt.full$biotype %in% c("misc_RNA")), sum(sptt.full$biotype == "Mt_tRNA"), sum(sptt.full$biotype == "scRNA"), sum(sptt.full$biotype == "snRNA")), 
                c(sum(visium.full$biotype == "snoRNA"), sum(visium.full$biotype == "miRNA"), sum(visium.full$biotype == "scaRNA"), sum(visium.full$biotype %in% c("misc_RNA")), sum(visium.full$biotype == "Mt_tRNA"), sum(visium.full$biotype == "scRNA"), sum(visium.full$biotype == "snRNA")), 
                c(sum(stereo.full$biotype == "snoRNA"), sum(stereo.full$biotype == "miRNA"), sum(stereo.full$biotype == "scaRNA"), sum(stereo.full$biotype %in% c("misc_RNA")), sum(stereo.full$biotype == "Mt_tRNA"), sum(stereo.full$biotype == "scRNA"), sum(stereo.full$biotype == "snRNA")), 
                c(sum(decoder.full$biotype == "snoRNA"), sum(decoder.full$biotype == "miRNA"), sum(decoder.full$biotype == "scaRNA"), sum(decoder.full$biotype %in% c("misc_RNA")), sum(decoder.full$biotype == "Mt_tRNA"), sum(decoder.full$biotype == "scRNA"), sum(decoder.full$biotype == "snRNA")),
                c(sum(dbit.full$biotype == "snoRNA"), sum(dbit.full$biotype == "miRNA"), sum(dbit.full$biotype == "scaRNA"), sum(dbit.full$biotype %in% c("misc_RNA")), sum(dbit.full$biotype == "Mt_tRNA"), sum(dbit.full$biotype == "scRNA"), sum(dbit.full$biotype == "snRNA")), 
                c(sum(strs.full$biotype == "snoRNA"), sum(strs.full$biotype == "miRNA"), sum(strs.full$biotype == "scaRNA"), sum(strs.full$biotype %in% c("misc_RNA")), sum(strs.full$biotype == "Mt_tRNA"), sum(strs.full$biotype == "scRNA"), sum(strs.full$biotype == "snRNA")), 
                c(sum(slide.full$biotype == "snoRNA"), sum(slide.full$biotype == "miRNA"), sum(slide.full$biotype == "scaRNA"), sum(slide.full$biotype %in% c("misc_RNA")), sum(slide.full$biotype == "Mt_tRNA"), sum(slide.full$biotype == "scRNA"), sum(slide.full$biotype == "snRNA")))

p3.data = as.data.frame(p3.data)
colnames(p3.data) = c("snoRNA", "miRNA", "scaRNA", "misc-RNA", "MT-tRNA", "scRNA", "snRNA")

p3.data[1,] = p3.data[1,] / dim(sptt.full)[1]
p3.data[2,] = p3.data[2,] / dim(visium.full)[1]
p3.data[3,] = p3.data[3,] / dim(stereo.full)[1]
p3.data[4,] = p3.data[4,] / dim(decoder.full)[1]
p3.data[5,] = p3.data[5,] / dim(dbit.full)[1]
p3.data[6,] = p3.data[6,] / dim(strs.full)[1]
p3.data[7,] = p3.data[7,] / dim(slide.full)[1]

p3.data$expr = c("SPTT_seq", "10X_Visium", "Stereo_seq", "Decoder_seq", "Dbit_seq", "STRS", "Slide-seqV2")

p3.melt = reshape2::melt(p3.data)
p3.melt$expr = factor(p3.melt$expr, levels = c("SPTT_seq", "10X_Visium", "Stereo_seq", "Decoder_seq", "Dbit_seq", "STRS", "Slide-seqV2"))

rare_rnas_colors = c("#979998", "#C69287", "#E79A90", "#EFBC91", "#E4CD87", "#FAE5B8", "#DDDDDF")
rare_rnas_colors2 = c("#F57C6E", "#F2B56F", "#FAE69E", "#84C3B7", "#71B7ED", "#B8AEEB", "#F2A7DA")

p3 = ggplot(p3.melt, aes(x = expr, y = value * 100, fill = variable)) + geom_bar(stat = 'identity') + scale_fill_manual(values = rare_rnas_colors2) + theme_classic() + xlab('') + ylab('% of annotated reads') + coord_flip() + 
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 9),
        axis.text.y = element_text(color = 'black', size = 9))
p3
ggsave(p3, filename = "./non_coding_RNAs.pdf", width = 8, height = 4)
