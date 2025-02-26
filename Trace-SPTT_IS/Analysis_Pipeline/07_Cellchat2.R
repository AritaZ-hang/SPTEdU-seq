ptm = Sys.time()
rm(list = ls()); gc()
setwd("YOURWORKDIR")
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(CellChat)
library(patchwork)
options(stringAsFactor=FALSE)
options(future.globals.maxSize = 1000 * 1024^2)

# configures
conversion.factor=0.72
spot.size=10
spatial.factors=data.frame(ratio=conversion.factor,
                           tol=spot.size/2)

data.input = readRDS(file="./counts.RDS")
meta = readRDS(file = "./meta.RDS")
spatial.locs = readRDS(file="./coords.RDS")
colnames(spatial.locs) = c("x", "y")

data.input = NormalizeData(data.input) # Since the data we input is raw counts, it needs to be normalized.

spatial.locs = spatial.locs[rownames(meta),]


cellchat = createCellChat(object = data.input, 
                          meta = meta, 
                          group.by = "labels", 
                          datatype = "spatial", 
                          coordinates=spatial.locs, 
                          spatial.factors = spatial.factors)

spatialDimPlot(cellchat, group.by = "labels", point.size = 0.1)

CellChatDB = CellChatDB.mouse
showDatabaseCategory(CellChatDB)

CellChatDB.use = subsetDB(CellChatDB)

cellchat@DB = CellChatDB.use

cellchat = subsetData(cellchat)
future::plan("multisession", workers = 12)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat, variable.both=F)

cellchat = computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range=250, scale.distance=1, 
                             contact.dependent=TRUE, contact.range=10)

cellchat = filterCommunication(cellchat, min.cells=5)

cellchat = computeCommunProbPathway(cellchat)

cellchat = aggregateNet(cellchat)
df.net = subsetCommunication(cellchat)

saveRDS(cellchat, file="./cellchat.RDS")
saveRDS(df.net, file = "./cellchat_df_net.RDS")
