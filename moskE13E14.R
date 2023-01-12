library(Biobase)
library(GEOquery)
library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)


cellchat.E13 <- readRDS(url("https://ndownloader.figshare.com/files/25957094"))
cellchat.E13 <- updateCellChat(cellchat.E13)
cellchat.E14 <- readRDS(url("https://ndownloader.figshare.com/files/25957634"))
cellchat.E14 <- updateCellChat(cellchat.E14)
 

slotNames(cellchat.E13) 
#output:[1] "data.raw"       "data"           "data.signaling" "data.scale"     "data.project"  
#[6] "images"         "net"            "netP"           "meta"           "idents"        
#[11] "DB" 
#syntax:object@slotname

moskE13.data = cellchat.E13@data
moskE13.meta = cellchat.E13@meta
moskE13 <- CreateSeuratObject(counts = moskE13.data, meta.data = moskE13.meta)
moskE13

moskE14.data = cellchat.E14@data
moskE14.meta = cellchat.E14@meta
moskE14 <- CreateSeuratObject(counts = moskE14.data, meta.data = moskE14.meta)
moskE14
