#human data annotated for communication analysis using LIANA
library(base)
library(datasets)
library(Biobase)
library(GEOquery)
library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(magrittr)
library(liana)
library(ggplot2)

load(url("https://ndownloader.figshare.com/files/25950872"))
#husk_data <- load(url("https://ndownloader.figshare.com/files/25950872"))
#husk <- CreateSeuratObject(husk_data)
#husk_url <- "https://ndownloader.figshare.com/files/25950872"

#husk_data <- load(url("https://ndownloader.figshare.com/files/25950872"))
#husk_data <- readRDS("/Users/xiaoyizheng/Downloads/Praktikum/data_humanSkin_CellChat.rda")
#husk <- CreateSeuratObject(husk_data)


husk.data = data_humanSkin$data
husk.meta = data_humanSkin$meta
husk <- CreateSeuratObject(counts = husk.data, meta.data = husk.meta)

husk
husk.data
husk.meta

#husk <-RunUMAP(husk, reduction = "umap", dims = 1:10)

#DimPlot(husk, reduction = "pca", group.by = "cluster")


#Install LIANA
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

#remotes::install_github('saezlab/liana')
remotes::install_github('saezlab/liana', force = TRUE)

#liana takes Seurat and SingleCellExperiment objects as input, containing processed counts and clustered cells.
liana_path <- system.file(package = "liana")
#testdata <- readRDS("/Users/xiaoyizheng/Downloads/Praktikum/cellchat_embryonic_E13.rds")
liana_test <- liana_wrap(husk_an, resource = c("Consensus"))



#testdata <- readRDS(file.path(liana_path , "testdataE13", "input", "cellchat_embryonic_E13.rds"))
#testdata <- readRDS(data.dir = "/Users/xiaoyizheng/Downloads/Praktikum/cellchat_embryonic_E13.rds")

# Liana returns a list of results, each element of which corresponds to a method
liana_test %>% dplyr::glimpse()

liana_test <- liana_test %>% liana_aggregate()
dplyr::glimpse(liana_test)

##scoring: to prioritize interaction in a single sample system
##exception of SingleCellSignalR’s LRscore, take the specificity of the cluster pair into account

liana_test %>%
  liana_dotplot(source_groups = c("FIB"),
                target_groups = c("LC", "TC", "NKT"),
                ntop = 20)
#"APOE+ FIB", "FBN1+ FIB", "COL11A1+ FIB", "Inflam. FIB"