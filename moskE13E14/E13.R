library(Biobase)
library(GEOquery)
library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)


cellchat.E13 <- readRDS(url("https://ndownloader.figshare.com/files/25957094"))
cellchat.E13 <- updateCellChat(cellchat.E13)
slotNames(cellchat.E13) 
#output:[1] "data.raw"       "data"           "data.signaling" "data.scale"     "data.project"  
#[6] "images"         "net"            "netP"           "meta"           "idents"        
#[11] "DB" 
#syntax:object@slotname

moskE13.data = cellchat.E13@data
moskE13.meta = cellchat.E13@meta
moskE13 <- CreateSeuratObject(counts = moskE13.data, meta.data = moskE13.meta)
moskE13

#II. Pre-processing
moskE13[["percent.mt"]] <- PercentageFeatureSet(moskE13, pattern = "^MT-") #mitoc. genes -> Zeichen fuer schlechte Qualitaet
VlnPlot(moskE13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #All cells have the same value of percent.mt.(maybe bc already filtered?)

head(moskE13@meta.data, 5) #Show QC metrics for the first 5 cells

moskE13.data[c("Igf2", "Cck", "Pf4"), 1:30] #a few genes in the first 30 cells #picked "Igf2", "Cck", "Pf4" according to web

#visualize feature-feature relationships
plot1 <- FeatureScatter(moskE13, feature1 = "nCount_RNA", feature2 = "percent.mt") #since all have same mt quality...
plot2 <- FeatureScatter(moskE13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #positively correlated
plot1 + plot2

#filter "moskE13"
moskE13 <- subset(moskE13, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #hmmm but mt for filtering unneccessary

#III.normalization
moskE13 <- NormalizeData(moskE13, normalization.method = "LogNormalize", scale.factor = 10000)

#IV.feature selection
moskE13 <- FindVariableFeatures(moskE13, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(moskE13), 10) # Identify the 10 most highly variable genes
top10

top20 <- head(VariableFeatures(moskE13), 20) # Identify the 20 most highly variable genes
top20

plot1 <- VariableFeaturePlot(moskE13)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#V. Scaling
all.genes <- rownames(moskE13)
moskE13 <- ScaleData(moskE13, features = all.genes)

#optional: moskE13 <- ScaleData(moskE13)

#VI. Linear Dimension Reduction
moskE13 <- RunPCA(moskE13, features = VariableFeatures(object = moskE13))
print(moskE13[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(moskE13, dims = 1:2, reduction = "pca")
DimPlot(moskE13, reduction = "pca")
DimPlot(moskE13, reduction = "pca",label = TRUE)
DimHeatmap(moskE13, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
DimHeatmap(moskE13, dims = 1:15, cells = 500, balanced = TRUE)
#metafeature: combines information across a correlated feature set

#VII. Determination of Dimension
moskE13 <- JackStraw(moskE13, num.replicate = 100)
moskE13 <- ScoreJackStraw(moskE13, dims = 1:20)
JackStrawPlot(moskE13, dims = 1:15) 
ElbowPlot(moskE13) #ranking of principle components #how to read it?(the slope = significance?)
#first 8 are good

#VIII. Clustering Cells
moskE13 <- FindNeighbors(moskE13, dims = 1:10)
moskE13 <- FindClusters(moskE13, resolution = 0.3) #13 clusters aus paper #res = 0.3 cluster 0->12
head(Idents(moskE13), 5)

#IX. Non-linear Dimension Reduction
moskE13 <- RunUMAP(moskE13, dims = 1:10)
#DimPlot(moskE13, reduction = "umap")
DimPlot(moskE13, reduction = "umap", label = TRUE)
##saveRDS(moskE13, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/moskE13Output.rds") 
#CellChatSelf/moskE13/output/
##infoRDS("/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/moskE13Output.rds")

#X. Cluster Biomarkers
# cluster0.markers <- FindMarkers(moskE13, ident.1 = 0, min.pct = 0.25)
# head(cluster0.markers, n = 5)
# VlnPlot(moskE13, features = c("Col3a1", "Col5a2", "Col6a3", "Tnn", "Col5a1")) 
# 
# cluster1.markers <- FindMarkers(moskE13, ident.1 = 1, min.pct = 0.25)
# head(cluster1.markers, n = 5)
# VlnPlot(moskE13, features = c("Rpl7", "Rpl31", "Col3a1", "Col5a2", "Spats2l"))
# 
# cluster2.markers <- FindMarkers(moskE13, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# VlnPlot(moskE13, features = c("Col3a1", "Des", "Serpine2", "Rgs5", "Rgs4"))
# 
# 
# cluster3.markers <- FindMarkers(moskE13, ident.1 = 3, min.pct = 0.25)
# head(cluster3.markers, n = 5)
# VlnPlot(moskE13, features = c("Col3a1", "Col5a2", "Selp", "F11r", "Cd34")) #"Selp", "F11r", "Cd34" very significant

cluster4.markers <- FindMarkers(moskE13, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
VlnPlot(moskE13, features = c("Fcer1g", "Il1b", "Tyrobp", "Lyz2", "Fth1", "Ifi27l2a", "Cd74", "Ctss", "Wfdc17", "Cd52"))

# cluster5.markers <- FindMarkers(moskE13, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
# VlnPlot(moskE13, features = c("Fcer1g", "Il1b", "Ctss", "Laptm5", "Cd52"))

# cluster6.markers <- FindMarkers(moskE13, ident.1 = 6, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster6.markers, n = 6)
# VlnPlot(moskE13, features = c("", "", "", "", ""))

cluster7.markers <- FindMarkers(moskE13, ident.1 = 7, min.pct = 0.25)
markers_cluster7 = head(cluster7.markers, n = 10)
VlnPlot(moskE13, features = c("Rgs5", "Gm13889", "Notch3", "Col4a1", "Ndufa4l2", "Col4a2", "Il6", "Des", "Serpine2", "Col18a1"))

# 
# cluster8.markers <- FindMarkers(moskE13, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(moskE13, features = c("", "", "", "", ""))
# 
cluster9.markers <- FindMarkers(moskE13, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 10)
VlnPlot(moskE13, features = c("Icos", "Cd3d", "AW112010", "Hcst", "Cd52", "Rac2", "Ccl5", "Vps37b", "Ftl1", "Rora"))

cluster10.markers <- FindMarkers(moskE13, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 10)
VlnPlot(moskE13, features = c("Rgs5", "Il6", "Col4a1", "Col4a2", "Gm13889", "Fth1", "Serpine2", "Higd1b", "Tyrobp", "Ndufa4l2"),sort = 'decreasing')

# 
# cluster11.markers <- FindMarkers(moskE13, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(moskE13, features = c("", "", "", "", ""))
# 
# cluster12.markers <- FindMarkers(moskE13, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(moskE13, features = c("", "", "", "", ""))

VlnPlot(moskE13, features = c("Hdc", "G0s2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hba-a2", "Hbb-bs"), sort = 'decreasing')



#aus Paper supplm

library(ggplot2)

#FIB2:
VlnPlot(moskE13, features = c("Des", "Rgs5"), sort = 'decreasing')

#MYL:
VlnPlot(moskE13, features = c("C1qb", "Pf4"), sort = 'decreasing')

#FIB-3:
VlnPlot(moskE13, features = c("Eln", "Ogn"), sort = 'decreasing')

#ENDO:
VlnPlot(moskE13, features = c("Cd93", "Pecam1"), sort = 'decreasing')

#FIB-4:
VlnPlot(moskE13, features = c("Birc5", "Ccnb2"), sort = 'decreasing')

#TC:
VlnPlot(moskE13, features = c("Icos", "Nkg7"), sort = 'decreasing')

#BC:
VlnPlot(moskE13, features = c("Ccr7", "H2-DMb1"), sort = 'decreasing')

#FIB-5:
VlnPlot(moskE13, features = c("Hdc",  "G0s2"), sort = 'decreasing')

#SCH:
VlnPlot(moskE13, features = c("Cadm4", "Itih5"), sort = 'decreasing') + labs(title = "SCH")

#RBC:
VlnPlot(moskE13, features = c("Hba-a2", "Hbb-bs"), sort = 'decreasing')

#DEN:
VlnPlot(moskE13, features = c("Acp5", "Mmp9"), sort = 'decreasing')

#LYME:
VlnPlot(moskE13, features = c("Ccl21a", "Lyve1"), sort = 'decreasing')

VlnPlot(moskE13, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1"), sort = 'decreasing')
VlnPlot(moskE13, features = c("Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4"))
VlnPlot(moskE13, features = c("Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))


#cluster2.markers <- FindMarkers(moskE13, ident.1 = 2, min.pct = 0.25)
#head(cluster2.markers, n = 5)
#cluster5.markers <- FindMarkers(moskE13, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)
#moskE13.markers <- FindAllMarkers(moskE13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #No features pass logfc.threshold threshold, maybe threshold=0.1, 0.15, 0.2?
moskE13.markers <- FindAllMarkers(moskE13, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
moskE13.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#cluster0.markers <- FindMarkers(moskE13, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#cluster0.markers <- FindMarkers(moskE13, ident.1 = 0, logfc.threshold = 0.2, test.use = "roc", only.pos = TRUE)
#VlnPlot(moskE13, features = c("Rgs5", "Tyrobp"))
#VlnPlot(moskE13, features = c("Crabp1", "Itm2a"))
#VlnPlot(moskE13, features = c("Ptn", "H19"))
#VlnPlot(moskE13, features = c("Saa3", "Crabp1"))
#VlnPlot(moskE13, features = c("Col1a1", "Col1a2"), slot = "counts", log = TRUE) #plot raw counts
FeaturePlot(moskE13, features = c("Sparc", "Tyrobp", "Igfbp7", "Pecam1", "Cd34", "Mfap4", "Ndufa4l2", "Il6", "Birc5")) #found after PCA
FeaturePlot(moskE13, features = c("Crabp1", "Itm2a", "Ptn", "H19", "Saa3", "Crabp1", "Rgs5", "Col4a1",
                               "Ctla2a", "Col4a1"))
FeaturePlot(moskE13, features = c("Col5a1", "Tnc"))

#features from paper supplm. 2.b
FeaturePlot(moskE13, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1"))
FeaturePlot(moskE13, features = c("Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4"))
FeaturePlot(moskE13, features = c("Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))



moskE13.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(moskE13, features = top10$gene) 
#+ NoLegend()

#cluster0.markers <- FindMarkers(moskE13, ident.1 = 2, min.pct = 0.25)
#head(cluster0.markers, n = 2)
#FeaturePlot(moskE13, features = c("Col5a1 ", "Tnc"))


#alternative way to assign IDs
# moskE13 <- RenameIdents(moskE13, `0` = "0_FIB-1", `1` = "1_FIB-3", `2` = "2_MYL",
#                       `3` = "3_FIB-2", `4` = "4_undecided", `5` = "5_LYME", `6` = "6_FIB-4", `7` = "7_undecided", `8` = "8_ENDO", `9` = "9_undecided",
#                       `10` = "10_undecided", `11` = "11_SCH", `12` = "12_DEN")
# DimPlot(moskE13, label = TRUE)

#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("FIB-1", "FIB-3", "MYL", "FIB-2", "BC", "LYME", "FIB-4", "RBC", "ENDO", "TC", "FIB-5", "SCH", "DEN")
names(new.cluster.ids) <- levels(moskE13)
moskE13 <- RenameIdents(moskE13, new.cluster.ids)
DimPlot(moskE13, reduction = "umap", label = TRUE, pt.size = 0.5)
#DimPlot(moskE13, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



#FindConservedMarkers(moskE13, idents.1= [nr.] )
#Idents(moskE13)
#RenameIdents(moskE13, '[nr.]'= [celltype] )


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
  moskE13,
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

allgenes <- rownames(moskE13)
matrix1 <- as.data.frame(moskE13@assays$RNA@data)
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

#moskE1312 uses mouse data so conversion to orthological human gene names
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(moskE13@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrix1 <- matrix1[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),]
matrix1$gene <- genesV2$Gene.stable.ID



