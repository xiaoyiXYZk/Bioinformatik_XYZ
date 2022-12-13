#I.Setup
library(dplyr)
library(Seurat)
library(patchwork)

mouse.data <- Read10X(data.dir = "/Users/xiaoyizheng/Downloads/Praktikum/Bioinfo_Einarbeiten/CellChatSelf/mouse/GSE113854_RAW/")
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

#VIII. Clustering Cells
mouse <- FindNeighbors(mouse, dims = 1:10)
mouse <- FindClusters(mouse, resolution = 0.5)
head(Idents(mouse), 5)

#IX. Non-linear Dimension Reduction
mouse <- RunUMAP(mouse, dims = 1:10)
#DimPlot(mouse, reduction = "umap")
DimPlot(mouse, reduction = "umap", label = TRUE)
saveRDS(mouse, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/CellChatSelf/mouse/output/mouseOutput.rds")


#X. Cluster Biomarkers

cluster2.markers <- FindMarkers(mouse, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
cluster5.markers <- FindMarkers(mouse, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
#mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #No features pass logfc.threshold threshold, maybe threshold=0.1, 0.15, 0.2?
mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
mouse.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#cluster0.markers <- FindMarkers(mouse, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster0.markers <- FindMarkers(mouse, ident.1 = 0, logfc.threshold = 0.2, test.use = "roc", only.pos = TRUE)
VlnPlot(mouse, features = c("Rgs5", "Tyrobp"))
VlnPlot(mouse, features = c("Col1a1", "Col1a2"), slot = "counts", log = TRUE) #plot raw counts
#FeaturePlot(mouse, features = c("Sparc", "Tyrobp", "Igfbp7", "Pecam1", "Cd34", "Mfap4", "Ndufa4l2", "Il6",
                               "Birc5"))
FeaturePlot(mouse, features = c("Crabp1", "Itm2a", "Ptn", "H19", "Saa3", "Crabp1", "Rgs5", "Col4a1",
                                "Ctla2a", "Col4a1"))
FeaturePlot(mouse, features = c("Col5a1", "Tnc"))

mouse.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(mouse, features = top10$gene) + NoLegend()


#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")
names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



