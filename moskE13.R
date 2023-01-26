#I. Setup
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
library(ggplot2)
library(devtools)



#loading data and create Seurat object
E13 <- readRDS(url("https://ndownloader.figshare.com/files/25957094"))
E13 <- updateCellChat(E13)
slotNames(E13) 

E13.data = E13@data
E13.meta = E13@meta
E13 <- CreateSeuratObject(counts = E13.data, meta.data = E13.meta)
E13

#II. Pre-processing

#quality control
E13[["percent.mt"]] <- PercentageFeatureSet(E13, pattern = "^MT-") #mitoc. genes -> Zeichen fuer schlechte Qualitaet
VlnPlot(E13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #All cells have the same value of percent.mt.(maybe bc already filtered?)

head(E13@meta.data, 5) #Show QC metrics for the first 5 cells

E13.data[c("Igf2", "Cck", "Pf4"), 1:30]

plot1 <- FeatureScatter(E13, feature1 = "nCount_RNA", feature2 = "percent.mt") #since all have same mt quality...
plot2 <- FeatureScatter(E13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #positively correlated
plot1 + plot2

#filter "E13"
E13 <- subset(E13, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #hmmm but mt for filtering unneccessary

#III.normalization
E13 <- NormalizeData(E13, normalization.method = "LogNormalize", scale.factor = 10000)

#IV.feature selection
E13 <- FindVariableFeatures(E13, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(E13), 10) # Identify the 10 most highly variable genes
top10

plot1 <- VariableFeaturePlot(E13)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#V. Scaling
all.genes <- rownames(E13)
E13 <- ScaleData(E13, features = all.genes)

#VI. Linear Dimension Reduction
E13 <- RunPCA(E13, features = VariableFeatures(object = E13))
print(E13[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E13, dims = 1:2, reduction = "pca")
DimPlot(E13, reduction = "pca")
#DimPlot(E13, reduction = "pca",label = TRUE)
DimHeatmap(E13, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
DimHeatmap(E13, dims = 1:15, cells = 500, balanced = TRUE)

#VII. Determination of Dimension
E13 <- JackStraw(E13, num.replicate = 100)
E13 <- ScoreJackStraw(E13, dims = 1:20)
JackStrawPlot(E13, dims = 1:15) 
ElbowPlot(E13)

#VIII. Clustering Cells
E13 <- FindNeighbors(E13, dims = 1:10)
E13 <- FindClusters(E13, resolution = 0.35) # 11 clusters aus paper 
head(Idents(E13), 5)

#IX. Non-linear Dimension Reduction
E13 <- RunUMAP(E13, dims = 1:10)
DimPlot(E13, reduction = "umap", label = TRUE)
##saveRDS(E13, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/E13Output.rds") 

#X. Cluster Biomarkers
#retrived from CellChat paper summplm.




VlnPlot(E13, features = c("Col1a1", "Lum", "Crabp1", "Lox",
                          "Sox2",
                          "Bmp4", "Acta2",
                          "Rgs5"), sort = 'decreasing')
#8 = MELA
#9 = Immune




###method to retrieve markers according to Seurat tutorial:
#E13.markers <- FindAllMarkers(E13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
#E13.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)

#feature plots for validation
#all in one:
#FeaturePlot(E13, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))

FeaturePlot(E13, features = c("", "", ...)

#heatmap with top 20 markers found by alt. method above, not from CellChat paper
#E13.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10
#DoHeatmap(E13, features = top10$gene) 
#+ NoLegend()

#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("FIB-1", "FIB-3", "MYL", "FIB-2", "BC", "LYME", "FIB-4", "RBC", "ENDO", "TC", "FIB-5", "SCH", "DEN")
names(new.cluster.ids) <- levels(E13)
E13 <- RenameIdents(E13, new.cluster.ids)
DimPlot(E13, reduction = "umap", label = TRUE, pt.size = 0.5)

#dot plot
markers.to.plot <- c("Crabp1", 
                     "Hba-a2", "Hbb-bs",
                     "Eln", "Ogn",
                     "Ccr7", "H2-DMb1",
                     "C1qb", "Pf4",
                     "Des", "Rgs5",
                     "Cd93", "Pecam1",
                     "Ccl21a", "Lyve1",
                     "Birc5", "Ccnb2",
                     "Nkg7",
                     "Hdc",  "G0s2",
                     "Acp5", "Mmp9")


DotPlot(
  E13,
  assay = NULL,
  features = markers.to.plot,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + RotatedAxis()                   


#XII. Identify significant interaction using LIANA

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github('saezlab/liana', force = TRUE)
#devtools::install_github('saezlab/liana', force = TRUE)
#Install LIANA
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#if (!requireNamespace("remotes", quietly = TRUE))
#  install.packages("remotes")
library(liana)
#remotes::install_github('saezlab/liana')
#remotes::install_github('saezlab/liana', force = TRUE)
liana_path <- system.file(package = "liana")
liana_test <- liana_wrap(E13, resource = c("Consensus"))

# Liana returns a list of results, each element of which corresponds to a method => we need result through CellPhone
liana_test %>% dplyr::glimpse()
#=> list of 5: results through natmi, connectome, logfc, sca, cellphonedb
#extracting results only those via CellPhoneDB
liana_test_cpdb <- liana_test$cellphonedb
liana_test_cpdb
CrossTalkR_input <- select(liana_test_cpdb, source, ligand, target, receptor, lr.mean)
write.csv(CrossTalkR_input, file="/Users/xiaoyizheng/Downloads/Praktikum/E1312_liana_res.csv")


#we probably dont want to aggregate since we only need for now from CellPhoneDB right?
#liana_test <- liana_test %>% liana_aggregate()
#dplyr::glimpse(liana_test)
#liana_test %>%
#  liana_dotplot(source_groups = c("FIB-1"),
#                target_groups = c("FIB-2", "FIB-3", "FIB-4"),
#                ntop = 20)

#XIII. CrossTalkR
#install_github("hadley/devtools")
#install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/CrossTalkeR", force = TRUE, build_vignettes = TRUE)
require(CrossTalkeR)
#library(CrossTalkeR)

suppressPackageStartupMessages({require(CrossTalkeR)})
suppressPackageStartupMessages({require(igraph)})
suppressPackageStartupMessages({require(ggraph)})
suppressPackageStartupMessages({require(ggplot2)})

#vignette('CrossTalkeR')

paths <- c('E1312' = system.file("/Users/xiaoyizheng/Downloads/Praktikum", "E1312_liana_res.csv", package = "CrossTalkeR"))
#paths <- system.file("E1312_liana_res.csv",package = "CrossTalkeR")
#paths
testdata <- read.csv("/Users/xiaoyizheng/Downloads/Praktikum/E1312_liana_res.csv")
testdata
genes <- c('Hba-a2')
data <- generate_report("/Users/xiaoyizheng/Downloads/Praktikum/E1312_liana_res.csv",
                        genes,
                        out_path=paste0(getwd(),'/html/'),
                        threshold=0,
                        out_file = 'E1312_test.html',
                        output_fmt = "html_document",
                        report = TRUE)
