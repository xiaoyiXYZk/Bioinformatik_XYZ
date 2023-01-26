# I. Setup
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
load(url("https://ndownloader.figshare.com/files/25950872"))

husk.data = data_humanSkin$data
husk.meta = data_humanSkin$meta
husk <- CreateSeuratObject(counts = husk.data, meta.data = husk.meta)

husk
husk.data
husk.data[c("A1BG", "A4GALT", "A4GALT"), 1:30]

#II. Pre-processing

#quality control
husk[["percent.mt"]] <- PercentageFeatureSet(husk, pattern = "^MT-") #mitoc. genes -> Zeichen fuer schlechte Qualitaet
VlnPlot(husk, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #All cells have the same value of percent.mt.(maybe bc already filtered?)

head(husk@meta.data, 5) #Show QC metrics for the first 5 cells



plot1 <- FeatureScatter(husk, feature1 = "nCount_RNA", feature2 = "percent.mt") #since all have same mt quality...
plot2 <- FeatureScatter(husk, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #positively correlated
plot1 + plot2

#filter "husk"
husk <- subset(husk, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #hmmm but mt for filtering unneccessary

#III.normalization
husk <- NormalizeData(husk, normalization.method = "LogNormalize", scale.factor = 10000)

#IV.feature selection
husk <- FindVariableFeatures(husk, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(husk), 10) # Identify the 10 most highly variable genes
top10

plot1 <- VariableFeaturePlot(husk)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#V. Scaling
all.genes <- rownames(husk)
husk <- ScaleData(husk, features = all.genes)

#VI. Linear Dimension Reduction
husk <- RunPCA(husk, features = VariableFeatures(object = husk))
print(husk[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(husk, dims = 1:2, reduction = "pca")
DimPlot(husk, reduction = "pca")
#DimPlot(husk, reduction = "pca",label = TRUE)
DimHeatmap(husk, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
DimHeatmap(husk, dims = 1:15, cells = 500, balanced = TRUE)

#VII. Determination of Dimension
husk <- JackStraw(husk, num.replicate = 100)
husk <- ScoreJackStraw(husk, dims = 1:20)
JackStrawPlot(husk, dims = 1:15) 
ElbowPlot(husk)

#VIII. Clustering Cells
husk <- FindNeighbors(husk, dims = 1:10)
husk <- FindClusters(husk, resolution = 0.7) #in paper: 10 clusters
head(Idents(husk), 5)

#IX. Non-linear Dimension Reduction
husk <- RunUMAP(husk, dims = 1:10)
DimPlot(husk, reduction = "umap", label = TRUE)
##saveRDS(husk, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/huskOutput.rds") 

#X. Cluster Biomarkers
#retrived from CellChat paper summplm.

#"Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"

#FIB2:
VlnPlot(husk, features = c("Des", "Rgs5"), sort = 'decreasing')

#MYL:
VlnPlot(husk, features = c("C1qb", "Pf4"), sort = 'decreasing')

#FIB-3:
VlnPlot(husk, features = c("Eln", "Ogn"), sort = 'decreasing')

#ENDO:
VlnPlot(husk, features = c("Cd93", "Pecam1"), sort = 'decreasing')

#FIB-4:
VlnPlot(husk, features = c("Birc5", "Ccnb2"), sort = 'decreasing')

#TC:
VlnPlot(husk, features = c("Icos", "Nkg7"), sort = 'decreasing')

#BC:
VlnPlot(husk, features = c("Ccr7", "H2-DMb1"), sort = 'decreasing')

#FIB-5:
VlnPlot(husk, features = c("Hdc",  "G0s2"), sort = 'decreasing')

#SCH:
VlnPlot(husk, features = c("Cadm4", "Itih5"), sort = 'decreasing') + labs(title = "SCH")

#RBC:
VlnPlot(husk, features = c("Hba-a2", "Hbb-bs"), sort = 'decreasing')

#DEN:
VlnPlot(husk, features = c("Acp5", "Mmp9"), sort = 'decreasing')

#LYME:
VlnPlot(husk, features = c("Ccl21a", "Lyve1"), sort = 'decreasing')

###method to retrieve markers according to Seurat tutorial:
#husk.markers <- FindAllMarkers(husk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
#husk.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)

#feature plots for validation
#all in one:
#FeaturePlot(husk, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))
FeaturePlot(husk, features = "Crabp1")
FeaturePlot(husk, features = c("Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7"))
FeaturePlot(husk, features = c("Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))

#heatmap with top 20 markers found by alt. method above, not from CellChat paper
#husk.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10
#DoHeatmap(husk, features = top10$gene) 
#+ NoLegend()

#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("FIB-1", "FIB-3", "MYL", "FIB-2", "BC", "LYME", "FIB-4", "RBC", "ENDO", "TC", "FIB-5", "SCH", "DEN")
names(new.cluster.ids) <- levels(husk)
husk <- RenameIdents(husk, new.cluster.ids)
DimPlot(husk, reduction = "umap", label = TRUE, pt.size = 0.5)

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
  husk,
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
liana_test <- liana_wrap(husk, resource = c("Consensus"))

# Liana returns a list of results, each element of which corresponds to a method => we need result through CellPhone
liana_test %>% dplyr::glimpse()
#=> list of 5: results through natmi, connectome, logfc, sca, cellphonedb
#extracting results only those via CellPhoneDB
liana_test_cpdb <- liana_test$cellphonedb
liana_test_cpdb
CrossTalkR_input <- select(liana_test_cpdb, source, ligand, target, receptor, lr.mean)
write.csv(CrossTalkR_input, file="/Users/xiaoyizheng/Downloads/Praktikum/husk12_liana_res.csv")


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

paths <- c('husk12' = system.file("/Users/xiaoyizheng/Downloads/Praktikum", "husk12_liana_res.csv", package = "CrossTalkeR"))
#paths <- system.file("husk12_liana_res.csv",package = "CrossTalkeR")
#paths
testdata <- read.csv("/Users/xiaoyizheng/Downloads/Praktikum/husk12_liana_res.csv")
testdata
genes <- c('Hba-a2')
data <- generate_report("/Users/xiaoyizheng/Downloads/Praktikum/husk12_liana_res.csv",
                        genes,
                        out_path=paste0(getwd(),'/html/'),
                        threshold=0,
                        out_file = 'husk12_test.html',
                        output_fmt = "html_document",
                        report = TRUE)
