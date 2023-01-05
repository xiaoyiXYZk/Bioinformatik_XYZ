#I.Setup
library(dplyr)
library(Seurat)
library(patchwork)
#mouse.data <- Read10X(data.dir = "/home/xiaoyi/Documents/Praktikum/Bioinformatik/GSE113854_RAW")
#mouse.data <- Read10X(data.dir = "/Users/xiaoyizheng/Downloads/Praktikum/Bioinfo_Einarbeiten/CellChatSelf/mouse/GSE113854_RAW/")
mouse.data

mouse <- CreateSeuratObject(counts = mouse.data, project = "MoSk", min.cells = 3, min.features = 200)
mouse

#II. Pre-processing
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^MT-") #mitoc. genes -> Zeichen fuer schlechte Qualitaet
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #All cells have the same value of percent.mt.(maybe bc already filtered?)

head(mouse@meta.data, 5) #Show QC metrics for the first 5 cells

mouse.data[c("Igf2", "Cck", "Pf4"), 1:30] #a few genes in the first 30 cells #picked "Igf2", "Cck", "Pf4" according to web

#visualize feature-feature relationships
plot1 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "percent.mt") #since all have same mt quality...
plot2 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #positively correlated
plot1 + plot2

#filter "mouse"
mouse <- subset(mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #hmmm but mt for filtering unneccessary

#III.normalization
mouse <- NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)

#IV.feature selection
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mouse), 10) # Identify the 10 most highly variable genes
top10

top20 <- head(VariableFeatures(mouse), 20) # Identify the 20 most highly variable genes
top20

plot1 <- VariableFeaturePlot(mouse)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#V. Scaling
all.genes <- rownames(mouse)
mouse <- ScaleData(mouse, features = all.genes)

#optional: mouse <- ScaleData(mouse)

#VI. Linear Dimension Reduction
mouse <- RunPCA(mouse, features = VariableFeatures(object = mouse))
print(mouse[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouse, dims = 1:2, reduction = "pca")
DimPlot(mouse, reduction = "pca")
DimPlot(mouse, reduction = "pca",label = TRUE)
DimHeatmap(mouse, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
DimHeatmap(mouse, dims = 1:15, cells = 500, balanced = TRUE)
#metafeature: combines information across a correlated feature set

#VII. Determination of Dimension
mouse <- JackStraw(mouse, num.replicate = 100)
mouse <- ScoreJackStraw(mouse, dims = 1:20)
JackStrawPlot(mouse, dims = 1:15) 
ElbowPlot(mouse) #ranking of principle components #how to read it?(the slope = significance?)
#first 8 are good

#VIII. Clustering Cells
mouse <- FindNeighbors(mouse, dims = 1:10)
mouse <- FindClusters(mouse, resolution = 0.3) #13 clusters aus paper #res = 0.3 cluster 0->12
head(Idents(mouse), 5)

#IX. Non-linear Dimension Reduction
mouse <- RunUMAP(mouse, dims = 1:10)
#DimPlot(mouse, reduction = "umap")
DimPlot(mouse, reduction = "umap", label = TRUE)
#saveRDS(mouse, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/CellChatSelf/mouse/output/mouseOutput.rds")


#X. Cluster Biomarkers
cluster0.markers <- FindMarkers(mouse, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
VlnPlot(mouse, features = c("Col3a1", "Col5a2", "Col6a3", "Tnn", "Col5a1")) 

cluster1.markers <- FindMarkers(mouse, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
VlnPlot(mouse, features = c("Rpl7", "Rpl31", "Col3a1", "Col5a2", "Spats2l"))

cluster2.markers <- FindMarkers(mouse, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
VlnPlot(mouse, features = c("Col3a1", "Des", "Serpine2", "Rgs5", "Rgs4"))


cluster3.markers <- FindMarkers(mouse, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
VlnPlot(mouse, features = c("Col3a1", "Col5a2", "Selp", "F11r", "Cd34")) #"Selp", "F11r", "Cd34" very significant

cluster4.markers <- FindMarkers(mouse, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)
VlnPlot(mouse, features = c("Col3a1", "Mdk", "Tmsb4x", "Eln", "Col1a2"))

cluster5.markers <- FindMarkers(mouse, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
VlnPlot(mouse, features = c("Fcer1g", "Il1b", "Ctss", "Laptm5", "Cd52"))

# cluster6.markers <- FindMarkers(mouse, ident.1 = 6, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster6.markers, n = 6)
# VlnPlot(mouse, features = c("", "", "", "", ""))
# 
# cluster7.markers <- FindMarkers(mouse, ident.1 = 7, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = 7)
# VlnPlot(mouse, features = c("", "", "", "", ""))
# 
# cluster8.markers <- FindMarkers(mouse, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mouse, features = c("", "", "", "", ""))
# 
# cluster9.markers <- FindMarkers(mouse, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mouse, features = c("", "", "", "", ""))
# 
# cluster10.markers <- FindMarkers(mouse, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mouse, features = c("", "", "", "", ""))
# 
# cluster11.markers <- FindMarkers(mouse, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mouse, features = c("", "", "", "", ""))
# 
# cluster12.markers <- FindMarkers(mouse, ident.1 = , ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster.markers, n = )
# VlnPlot(mouse, features = c("", "", "", "", ""))





#aus Paper supplm
VlnPlot(mouse, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1"))
VlnPlot(mouse, features = c("Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4"))
VlnPlot(mouse, features = c("Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))


#cluster2.markers <- FindMarkers(mouse, ident.1 = 2, min.pct = 0.25)
#head(cluster2.markers, n = 5)
#cluster5.markers <- FindMarkers(mouse, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)
#mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #No features pass logfc.threshold threshold, maybe threshold=0.1, 0.15, 0.2?
mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
mouse.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#cluster0.markers <- FindMarkers(mouse, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#cluster0.markers <- FindMarkers(mouse, ident.1 = 0, logfc.threshold = 0.2, test.use = "roc", only.pos = TRUE)
#VlnPlot(mouse, features = c("Rgs5", "Tyrobp"))
#VlnPlot(mouse, features = c("Crabp1", "Itm2a"))
#VlnPlot(mouse, features = c("Ptn", "H19"))
#VlnPlot(mouse, features = c("Saa3", "Crabp1"))
#VlnPlot(mouse, features = c("Col1a1", "Col1a2"), slot = "counts", log = TRUE) #plot raw counts
FeaturePlot(mouse, features = c("Sparc", "Tyrobp", "Igfbp7", "Pecam1", "Cd34", "Mfap4", "Ndufa4l2", "Il6", "Birc5")) #found after PCA
FeaturePlot(mouse, features = c("Crabp1", "Itm2a", "Ptn", "H19", "Saa3", "Crabp1", "Rgs5", "Col4a1",
                                "Ctla2a", "Col4a1"))
FeaturePlot(mouse, features = c("Col5a1", "Tnc"))

#features from paper supplm. 2.b
FeaturePlot(mouse, features = c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1"))
FeaturePlot(mouse, features = c("Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4"))
FeaturePlot(mouse, features = c("Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1"))



mouse.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(mouse, features = top10$gene) 
#+ NoLegend()

#cluster0.markers <- FindMarkers(mouse, ident.1 = 2, min.pct = 0.25)
#head(cluster0.markers, n = 2)
#FeaturePlot(mouse, features = c("Col5a1 ", "Tnc"))


#alternative way to assign IDs
mouse <- RenameIdents(mouse, `0` = "0_FIB-1", `1` = "1_FIB-3", `2` = "2_MYL",
                      `3` = "3_FIB-2", `4` = "4_undecided", `5` = "5_LYME", `6` = "6_FIB-4", `7` = "7_undecided", `8` = "8_ENDO", `9` = "9_undecided",
                      `10` = "10_undecided", `11` = "11_SCH", `12` = "12_DEN")
DimPlot(mouse, label = TRUE)

#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("FIB-1", "FIB-3", "MYL", "FIB-2", "4", "ENDO or LYME", "FIB-4", "7", "ENDO", "9", "10", "SCH", "DEN")
names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5)
#DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



#FindConservedMarkers(mouse, idents.1= [nr.] )
#Idents(mouse)
#RenameIdents(mouse, '[nr.]'= [celltype] )


#dot plot
markers.to.plot <- c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1","Birc5", "Ccnb2", "Icos", "Nkg7", "Ccr7", "H2-DMb1", "Hdc",  "G0s2", "Cadm4", "Itih5", "Hba-a2", "Hbb-bs", "Acp5", "Mmp9", "Ccl21a", "Lyve1")
DotPlot(
  mouse,
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
                     



install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/CrossTalkeR", build_vignettes = TRUE)                     
require(CrossTalkeR)  
vignette('CrossTalkeR')

# paths <- c('CTR' = system.file("extdata",
#                                "ctr_nils_bm_human_newformat.csv",
#                                package = "CrossTalkeR"),
#            'EXP' = system.file("extdata",
#                                "exp_nils_bm_human_newformat.csv",
#                                package = "CrossTalkeR"))

genes <- c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1")

#output <- system.file("extdata", package = "CrossTalkeR")

# data <- generate_report(paths,
#                         genes,
#                         out_path=paste0(output,'/'),
#                         threshold=0,
#                         out_file = 'vignettes_example.html',
#                         output_fmt = "html_document",
#                         report = TRUE)


suppressPackageStartupMessages({require(CrossTalkeR)})
›