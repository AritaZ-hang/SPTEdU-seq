library(Seurat)
library(velocyto.R)
library(tidyverse)
library(data.table)

setwd("/your/workdir/")

# SPTT for example

sptt = read.loom.matrices("SPTT.loom")
sptt_unspliced_reads = as.data.frame(Matrix::colSums(sptt$unspliced))
sptt_spliced_reads = as.data.frame(Matrix::colSums(sptt$spliced))
sptt_ambiguous_reads = as.data.frame(Matrix::colSums(sptt$ambiguous))

sptt_data = data.frame(unspliced = sptt_unspliced_reads$`Matrix::colSums(sptt$unspliced)`, 
                          spliced = sptt_spliced_reads$`Matrix::colSums(sptt$spliced)`, 
                          ambiguous = sptt_ambiguous_reads$`Matrix::colSums(sptt$ambiguous)`,
                          total = sptt_spliced_reads$`Matrix::colSums(sptt$spliced)` + sptt_unspliced_reads$`Matrix::colSums(sptt$unspliced)` + sptt_ambiguous_reads$`Matrix::colSums(sptt$ambiguous)`)
sptt_data$ratio = sptt_data$unspliced / sptt_data$total
