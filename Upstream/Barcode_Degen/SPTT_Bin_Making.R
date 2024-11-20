#!/usr/bin/Rscript
rm(list = ls()); gc()
library(data.table)
library(tidyverse)

argv = commandArgs(trailingOnly=T)

degen_coords_file=argv[1]
workdir = argv[2]


degen_coords = as.data.frame(fread(file = degen_coords_file))
colnames(degen_coords) = c("raw_barcode", "degen_barcodes", "xcoord", "ycoord")

# Step 1

scale = 0.72
pixel_per_bin_10 = round(10 / scale, 2) # 13.89
pixel_per_bin_25 = round(25 / scale, 2) # 34.72
pixel_per_bin_50 = round(50 / scale, 2) # 69.44
pixel_per_bin_100 = round(100 / scale, 2) # 138.89
pixel_per_bin_110 = round(110 / scale, 2) # 152.78

# Step2

cal_bin_edges = function(pixel_per_bin, coords)
{
  # modified from slide-seq v2
  min_xcoord = min(coords[["xcoord"]])
  max_xcoord = max(coords[["xcoord"]])
  min_ycoord = min(coords[["ycoord"]])
  max_ycoord = max(coords[["ycoord"]])
  
  # calculate bin edges
  xbins = seq(min_xcoord, max_xcoord, by = pixel_per_bin)
  ybins = seq(min_ycoord, max_ycoord, by = pixel_per_bin)
  
  # get centers for new coordinates
  xcoords = xbins + pixel_per_bin / 2
  ycoords = ybins + pixel_per_bin / 2
  
  # add final coordinates so the next step works
  xbins = c(xbins, max_xcoord)
  ybins = c(ybins, max_ycoord)
  
  # calculate new coordinates
  coords[["new_x"]] = sapply(coords[["xcoord"]], function(x) xcoords[which.max(xbins > x) - 1])
  coords[["new_y"]] = sapply(coords[["ycoord"]], function(y) ycoords[which.max(ybins > y) - 1])
  coords[["new_x"]] = as.numeric(coords[["new_x"]])
  coords[["new_y"]] = as.numeric(coords[["new_y"]])
  # remove points with na coordinates
  coords = na.omit(coords)
  
  # the result
  return(coords)
  
}

# Bin10

coords_bin10 = cal_bin_edges(pixel_per_bin = pixel_per_bin_10, degen_coords)
coords_bin10_use = coords_bin10[, c("new_x", "new_y", "raw_barcode", "degen_barcodes")]
save(coords_bin10, coords_bin10_use, file = paste0(workdir, "/", "bin10_coordinates.RData"))
write.csv(coords_bin10_use, file = paste0(workdir, "/", "bin10_coordinates.csv"), row.names = F)

# Bin 50
coords_bin50 = cal_bin_edges(pixel_per_bin = pixel_per_bin_50, degen_coords)
coords_bin50_use = coords_bin50[, c("new_x", "new_y", "raw_barcode", "degen_barcodes")]
save(coords_bin50, coords_bin50_use, file = paste0(workdir, "/", "bin50_coordinates.RData"))
write.csv(coords_bin50_use, file = paste0(workdir, "/", "bin50_coordinates.csv"), row.names = F)

# Bin 100
coords_bin100 = cal_bin_edges(pixel_per_bin = pixel_per_bin_100, degen_coords)
coords_bin100_use = coords_bin100[, c("new_x", "new_y", "raw_barcode", "degen_barcodes")]
save(coords_bin100, coords_bin100_use, file = paste0(workdir, "/", "bin100_coordinates.RData"))
write.csv(coords_bin100_use, file = paste0(workdir, "/", "bin100_coordinates.csv"), row.names = F)