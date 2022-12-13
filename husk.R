library(Biobase)
library(GEOquery)
library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)

load(url("https://ndownloader.figshare.com/files/25950872"))

data.input = data_humanSkin$data
meta = data_humanSkin$meta
husk <- CreateSeuratObject(data.input, meta.data = meta)

husk

#husk.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
husk.data[ , 1:30]
husk
