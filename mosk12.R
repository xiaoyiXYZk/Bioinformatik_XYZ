#I.Setup
library(dplyr)
library(Seurat)
library(patchwork)
#PC
#mosk.data <- Read10X(data.dir = "/home/xiaoyi/Documents/Praktikum/Bioinformatik/GSE113854_RAW")
#MAC
mosk.data <- Read10X(data.dir = "/Users/xiaoyizheng/Downloads/Praktikum/Bioinfo_Einarbeiten/CellChatSelf/mouse/GSE113854_RAW/")
mosk.data

mosk <- CreateSeuratObject(counts = mosk.data, project = "MoSk", min.cells = 3, min.features = 200)
mosk

#II. Pre-processing
mosk[["percent.mt"]] <- PercentageFeatureSet(mosk, pattern = "^MT-") #mitoc. genes -> Zeichen fuer schlechte Qualitaet
VlnPlot(mosk, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #All cells have the same value of percent.mt.(maybe bc already filtered?)

head(mosk@meta.data, 5) #Show QC metrics for the first 5 cells

mosk.data[c("Igf2", "Cck", "Pf4"), 1:30] #a few genes in the first 30 cells #picked "Igf2", "Cck", "Pf4" according to web

#visualize feature-feature relationships
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

top20 <- head(VariableFeatures(mosk), 20) # Identify the 20 most highly variable genes
top20

plot1 <- VariableFeaturePlot(mosk)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#V. Scaling
all.genes <- rownames(mosk)
mosk <- ScaleData(mosk, features = all.genes)

#optional: mosk <- ScaleData(mosk)

#VI. Linear Dimension Reduction
mosk <- RunPCA(mosk, features = VariableFeatures(object = mosk))
print(mosk[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mosk, dims = 1:2, reduction = "pca")
DimPlot(mosk, reduction = "pca")
DimPlot(mosk, reduction = "pca",label = TRUE)
DimHeatmap(mosk, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
DimHeatmap(mosk, dims = 1:15, cells = 500, balanced = TRUE)
#metafeature: combines information across a correlated feature set

#VII. Determination of Dimension
mosk <- JackStraw(mosk, num.replicate = 100)
mosk <- ScoreJackStraw(mosk, dims = 1:20)
JackStrawPlot(mosk, dims = 1:15) 
ElbowPlot(mosk) #ranking of principle components #how to read it?(the slope = significance?)
#first 8 are good

#VIII. Clustering Cells
mosk <- FindNeighbors(mosk, dims = 1:10)
mosk <- FindClusters(mosk, resolution = 0.3) #13 clusters aus paper #res = 0.3 cluster 0->12
head(Idents(mosk), 5)

#IX. Non-linear Dimension Reduction
mosk <- RunUMAP(mosk, dims = 1:10)
#DimPlot(mosk, reduction = "umap")
DimPlot(mosk, reduction = "umap", label = TRUE)
##saveRDS(mosk, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/moskOutput.rds") 
#CellChatSelf/mosk/output/
##infoRDS("/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/moskOutput.rds")

#X. Cluster Biomarkers
# cluster0.markers <- FindMarkers(mosk, ident.1 = 0, min.pct = 0.25)
# head(cluster0.markers, n = 5)
# VlnPlot(mosk, features = c("Col3a1", "Col5a2", "Col6a3", "Tnn", "Col5a1")) 
# 
# cluster1.markers <- FindMarkers(mosk, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)
# VlnPlot(mosk, features = c("Rpl7", "Rpl31", "Col3a1", "Col5a2", "Spats2l"))
# 
# cluster2.markers <- FindMarkers(mosk, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# VlnPlot(mosk, features = c("Col3a1", "Des", "Serpine2", "Rgs5", "Rgs4"))
# 
# 
# cluster3.markers <- FindMarkers(mosk, ident.1 = 3, min.pct = 0.25)
# head(cluster3.markers, n = 5)
# VlnPlot(mosk, features = c("Col3a1", "Col5a2", "Selp", "F11r", "Cd34")) #"Selp", "F11r", "Cd34" very significant

cluster4.markers <- FindMarkers(mosk, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
VlnPlot(mosk, features = c("Fcer1g", "Il1b", "Tyrobp", "Lyz2", "Fth1", "Ifi27l2a", "Cd74", "Ctss", "Wfdc17", "Cd52"))

# cluster5.markers <- FindMarkers(mosk, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
# VlnPlot(mosk, features = c("Fcer1g", "Il1b", "Ctss", "Laptm5", "Cd52"))

# cluster6.markers <- FindMarkers(mosk, ident.1 = 6, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster6.markers, n = 6)
# VlnPlot(mosk, features = c("", "", "", "", ""))

cluster7.markers <- FindMarkers(mosk, ident.1 = 7, min.pct = 0.25)
markers_cluster7 = head(cluster7.markers, n = 10)
VlnPlot(mosk, features = c("Rgs5", "Gm13889", "Notch3", "Col4a1", "Ndufa4l2", "Col4a2", "Il6", "Des", "Serpine2", "Col18a1"))

# 
# cluster8.markers <- FindMarkers(mosk, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mosk, features = c("", "", "", "", ""))
# 
cluster9.markers <- FindMarkers(mosk, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 10)
VlnPlot(mosk, features = c("Icos", "Cd3d", "AW112010", "Hcst", "Cd52", "Rac2", "Ccl5", "Vps37b", "Ftl1", "Rora"))

cluster10.markers <- FindMarkers(mosk, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 10)
VlnPlot(mosk, features = c("Rgs5", "Il6", "Col4a1", "Col4a2", "Gm13889", "Fth1", "Serpine2", "Higd1b", "Tyrobp", "Ndufa4l2"),sort = 'decreasing')

# 
# cluster11.markers <- FindMarkers(mosk, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mosk, features = c("", "", "", "", ""))
# 
# cluster12.markers <- FindMarkers(mosk, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mosk, features = c("", "", "", "", ""))

VlnPlot(mosk, features = c("Hdc", "G0s2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hba-a2", "Hbb-bs"), sort = 'decreasing')



#aus Paper supplm

library(ggplot2)

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

VlnPlot(mosk, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1"), sort = 'decreasing')
VlnPlot(mosk, features = c("Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4"))
VlnPlot(mosk, features = c("Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))


#cluster2.markers <- FindMarkers(mosk, ident.1 = 2, min.pct = 0.25)
#head(cluster2.markers, n = 5)
#cluster5.markers <- FindMarkers(mosk, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)
#mosk.markers <- FindAllMarkers(mosk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #No features pass logfc.threshold threshold, maybe threshold=0.1, 0.15, 0.2?
mosk.markers <- FindAllMarkers(mosk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
mosk.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#cluster0.markers <- FindMarkers(mosk, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#cluster0.markers <- FindMarkers(mosk, ident.1 = 0, logfc.threshold = 0.2, test.use = "roc", only.pos = TRUE)
#VlnPlot(mosk, features = c("Rgs5", "Tyrobp"))
#VlnPlot(mosk, features = c("Crabp1", "Itm2a"))
#VlnPlot(mosk, features = c("Ptn", "H19"))
#VlnPlot(mosk, features = c("Saa3", "Crabp1"))
#VlnPlot(mosk, features = c("Col1a1", "Col1a2"), slot = "counts", log = TRUE) #plot raw counts
FeaturePlot(mosk, features = c("Sparc", "Tyrobp", "Igfbp7", "Pecam1", "Cd34", "Mfap4", "Ndufa4l2", "Il6", "Birc5")) #found after PCA
FeaturePlot(mosk, features = c("Crabp1", "Itm2a", "Ptn", "H19", "Saa3", "Crabp1", "Rgs5", "Col4a1",
                                "Ctla2a", "Col4a1"))
FeaturePlot(mosk, features = c("Col5a1", "Tnc"))

#features from paper supplm. 2.b
FeaturePlot(mosk, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1"))
FeaturePlot(mosk, features = c("Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4"))
FeaturePlot(mosk, features = c("Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))



mosk.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(mosk, features = top10$gene) 
#+ NoLegend()

#cluster0.markers <- FindMarkers(mosk, ident.1 = 2, min.pct = 0.25)
#head(cluster0.markers, n = 2)
#FeaturePlot(mosk, features = c("Col5a1 ", "Tnc"))


#alternative way to assign IDs
# mosk <- RenameIdents(mosk, `0` = "0_FIB-1", `1` = "1_FIB-3", `2` = "2_MYL",
#                       `3` = "3_FIB-2", `4` = "4_undecided", `5` = "5_LYME", `6` = "6_FIB-4", `7` = "7_undecided", `8` = "8_ENDO", `9` = "9_undecided",
#                       `10` = "10_undecided", `11` = "11_SCH", `12` = "12_DEN")
# DimPlot(mosk, label = TRUE)

#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("FIB-1", "FIB-3", "MYL", "FIB-2", "BC", "LYME", "FIB-4", "RBC", "ENDO", "TC", "FIB-5", "SCH", "DEN")
names(new.cluster.ids) <- levels(mosk)
mosk <- RenameIdents(mosk, new.cluster.ids)
DimPlot(mosk, reduction = "umap", label = TRUE, pt.size = 0.5)
#DimPlot(mosk, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



#FindConservedMarkers(mosk, idents.1= [nr.] )
#Idents(mosk)
#RenameIdents(mosk, '[nr.]'= [celltype] )


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
                     



#install.packages("devtools")
#devtools::install_github("https://github.com/CostaLab/CrossTalkeR", build_vignettes = TRUE)                     
require(CrossTalkeR)  
#CrossTalkeR cookbook:
#vignette('CrossTalkeR')

# paths <- c('CTR' = system.file("extdata",
#                                "ctr_nils_bm_human_newformat.csv",
#                                package = "CrossTalkeR"),
#            'EXP' = system.file("extdata",
#                                "exp_nils_bm_human_newformat.csv",
#                                package = "CrossTalkeR"))

#genes <- c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1")

#output <- system.file("extdata", package = "CrossTalkeR")

# data <- generate_report(paths,
#                         genes,
#                         out_path=paste0(output,'/'),
#                         threshold=0,
#                         out_file = 'vignettes_example.html',
#                         output_fmt = "html_document",
#                         report = TRUE)


#suppressPackageStartupMessages({require(CrossTalkeR)})



require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
library(biomaRt)

allgenes <- rownames(mosk)
matrix1 <- as.data.frame(mosk@assays$RNA@data)
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

#mosk12 uses mouse data so conversion to orthological human gene names
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(mosk@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrix1 <- matrix1[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),]
matrix1$gene <- genesV2$Gene.stable.ID



