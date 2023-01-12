library(Biobase)
library(GEOquery)
library(Seurat)

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

load(url("https://ndownloader.figshare.com/files/25950872"))

data.input = data_humanSkin$data
meta = data_humanSkin$meta
seurat_obj <- CreateSeuratObject(data.input, meta.data = meta)