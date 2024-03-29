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
#library(liana)
library(ggplot2)
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github('saezlab/liana', force = TRUE)
#devtools::install_github('saezlab/liana', force = TRUE)

#loading data and create Seurat object
mosk.data <- Read10X(data.dir = "/Users/xiaoyizheng/Downloads/Praktikum/Bioinfo_Einarbeiten/CellChatSelf/mouse/GSE113854_RAW/")
mosk.data

mosk <- CreateSeuratObject(counts = mosk.data, project = "MoSk", min.cells = 3, min.features = 200)
mosk

#II. Pre-processing

#quality control
mosk[["percent.mt"]] <- PercentageFeatureSet(mosk, pattern = "^MT-") #mitoc. genes -> Zeichen fuer schlechte Qualitaet
VlnPlot(mosk, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #All cells have the same value of percent.mt.(maybe bc already filtered?)

head(mosk@meta.data, 5) #Show QC metrics for the first 5 cells

mosk.data[c("Igf2", "Cck", "Pf4"), 1:30]

plot1 <- FeatureScatter(mosk, feature1 = "nCount_RNA", feature2 = "percent.mt") #since all have same mt quality...
plot2 <- FeatureScatter(mosk, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #positively correlated
plot1 + plot2

#filter "mosk"
mosk <- subset(mosk, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #hmmm but mt for filtering unneccessary

#III.normalization
mosk <- NormalizeData(mosk, normalization.method = "LogNormalize", scale.factor = 10000)

#IV.feature selection
mosk <- FindVariableFeatures(mosk, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mosk), 10) # Identify the 10 most highly variable genes
top10

plot1 <- VariableFeaturePlot(mosk)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#V. Scaling
all.genes <- rownames(mosk)
mosk <- ScaleData(mosk, features = all.genes)

#VI. Linear Dimension Reduction
mosk <- RunPCA(mosk, features = VariableFeatures(object = mosk))
print(mosk[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mosk, dims = 1:2, reduction = "pca")
DimPlot(mosk, reduction = "pca")
#DimPlot(mosk, reduction = "pca",label = TRUE)
DimHeatmap(mosk, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
DimHeatmap(mosk, dims = 1:15, cells = 500, balanced = TRUE)

#VII. Determination of Dimension
mosk <- JackStraw(mosk, num.replicate = 100)
mosk <- ScoreJackStraw(mosk, dims = 1:20)
JackStrawPlot(mosk, dims = 1:15) 
ElbowPlot(mosk)

#VIII. Clustering Cells
mosk <- FindNeighbors(mosk, dims = 1:10)
mosk <- FindClusters(mosk, resolution = 0.3) #13 clusters aus paper #res = 0.3 cluster 0->12
head(Idents(mosk), 5)

#IX. Non-linear Dimension Reduction
mosk <- RunUMAP(mosk, dims = 1:10)
DimPlot(mosk, reduction = "umap", label = TRUE)
##saveRDS(mosk, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/moskOutput.rds") 

#X. Cluster Biomarkers
#retrived from CellChat paper summplm.

#"Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"

#FIB2:
VlnPlot(mosk, features = c("Des", "Rgs5"), sort = 'decreasing')

#MYL:
VlnPlot(mosk, features = c("C1qb", "Pf4"), sort = 'decreasing')

#FIB-3:
VlnPlot(mosk, features = c("Eln", "Ogn"), sort = 'decreasing')

#ENDO:
VlnPlot(mosk, features = c("Cd93", "Pecam1"), sort = 'decreasing')

#FIB-4:
VlnPlot(mosk, features = c("Birc5", "Ccnb2"), sort = 'decreasing')

#TC:
VlnPlot(mosk, features = c("Icos", "Nkg7"), sort = 'decreasing')

#BC:
VlnPlot(mosk, features = c("Ccr7", "H2-DMb1"), sort = 'decreasing')

#FIB-5:
VlnPlot(mosk, features = c("Hdc",  "G0s2"), sort = 'decreasing')

#SCH:
VlnPlot(mosk, features = c("Cadm4", "Itih5"), sort = 'decreasing') + labs(title = "SCH")

#RBC:
VlnPlot(mosk, features = c("Hba-a2", "Hbb-bs"), sort = 'decreasing')

#DEN:
VlnPlot(mosk, features = c("Acp5", "Mmp9"), sort = 'decreasing')

#LYME:
VlnPlot(mosk, features = c("Ccl21a", "Lyve1"), sort = 'decreasing')

###method to retrieve markers according to Seurat tutorial:
#mosk.markers <- FindAllMarkers(mosk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
#mosk.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)

#feature plots for validation
#all in one:
#FeaturePlot(mosk, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))
FeaturePlot(mosk, features = "Crabp1")
FeaturePlot(mosk, features = c("Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7"))
FeaturePlot(mosk, features = c("Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))

#heatmap with top 20 markers found by alt. method above, not from CellChat paper
#mosk.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10
#DoHeatmap(mosk, features = top10$gene) 
#+ NoLegend()

#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("FIB-1", "FIB-3", "MYL", "FIB-2", "BC", "LYME", "FIB-4", "RBC", "ENDO", "TC", "FIB-5", "SCH", "DEN")
names(new.cluster.ids) <- levels(mosk)
mosk <- RenameIdents(mosk, new.cluster.ids)
DimPlot(mosk, reduction = "umap", label = TRUE, pt.size = 0.5)

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
  mosk,
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
#Install LIANA
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#if (!requireNamespace("remotes", quietly = TRUE))
#  install.packages("remotes")
library(liana)
#remotes::install_github('saezlab/liana')
#remotes::install_github('saezlab/liana', force = TRUE)
liana_path <- system.file(package = "liana")
liana_test <- liana_wrap(mosk, resource = c("Consensus"))

# Liana returns a list of results, each element of which corresponds to a method => we need result through CellPhone
liana_test %>% dplyr::glimpse()
#=> list of 5: results through natmi, connectome, logfc, sca, cellphonedb
#extracting results only those via CellPhoneDB
liana_test_cpdb <- liana_test$cellphonedb
liana_test_cpdb
CrossTalkR_input <- select(liana_test_cpdb, source, ligand, target, receptor, lr.mean)
write.csv(CrossTalkR_input, file="/Users/xiaoyizheng/Downloads/Praktikum/mosk12_liana_res.csv")


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

paths <- c('mosk12' = system.file("/Users/xiaoyizheng/Downloads/Praktikum", "mosk12_liana_res.csv", package = "CrossTalkeR"))
#paths <- system.file("mosk12_liana_res.csv",package = "CrossTalkeR")
#paths
testdata <- read.csv("/Users/xiaoyizheng/Downloads/Praktikum/mosk12_liana_res.csv")
testdata
genes <- c('Hba-a2')
data <- generate_report("/Users/xiaoyizheng/Downloads/Praktikum/mosk12_liana_res.csv",
                        genes,
                        out_path=paste0(getwd(),'/html/'),
                        threshold=0,
                        out_file = 'mosk12_test.html',
                        output_fmt = "html_document",
                        report = TRUE)
