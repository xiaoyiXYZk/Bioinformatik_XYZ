library(Biobase)
library(GEOquery)
library(Seurat)
library(CellChat)
library(dplyr)
library(patchwork)


load(url("https://ndownloader.figshare.com/files/25950872"))

husk.data = data_humanSkin$data
husk.meta = data_humanSkin$meta
husk <- CreateSeuratObject(counts = husk.data, meta.data = husk.meta)

husk
husk.data
husk.data[c("A1BG", "A4GALT", "A4GALT"), 1:30]


#II. Pre-processing
husk[["percent.mt"]] <- PercentageFeatureSet(husk, pattern = "^MT-") 
VlnPlot(husk, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #probably already good quality

head(husk@meta.data, 5) #Show QC metrics for the first 5 cells
head(husk@meta.data, 10)

husk.data[c("CCL2", "CCL19", "CCR7"), 1:30] #a few genes in the first 30 cells #genes picked from GEO description GSE147424
husk.data[c("CCL2", "CCL19", "CCR7"), 2000:2020]
#visualize feature-feature relationships
plot1 <- FeatureScatter(husk, feature1 = "nCount_RNA", feature2 = "percent.mt") #since all have same mt quality...
plot2 <- FeatureScatter(husk, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #positively correlated
plot1 + plot2

#filter "husk"
husk <- subset(husk, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 

#III Normalization
husk <- NormalizeData(husk)

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

#optional: husk <- ScaleData(husk)

#Between V. and VI.: spliting the data into two groups according to whether NL or LS
#two Seurat objects for comparative analysis later on in CrossTalkR
huskNL <- husk[which(husk@meta.data$condition == "NL"),]
huskLS <- husk[which(husk@meta.data$condition == "LS"),]

huskNL <- ScaleData(huskNL)
huskLS <- ScaleData(huskLS)

huskNL <- FindVariableFeatures(huskNL, selection.method = "vst", nfeatures = 2000)
top10NL <- head(VariableFeatures(huskNL), 10) # Identify the 10 most highly variable genes
top10NL

huskLS <- FindVariableFeatures(huskLS, selection.method = "vst", nfeatures = 2000)
top10LS <- head(VariableFeatures(huskLS), 10) # Identify the 10 most highly variable genes
top10LS

top10_NL = c("CCL1", "CCL19", "CCL17", "CD74", "IL22",  "IL9",  "IL32", "DCN",  "CCL5", "CD52")
top10_LS = c("COCH", "HLA-DRA", "CFD", "APOD", "GNLY", "CXCL14", "HLA-DPA1", "HLA-DPB1", "FABP7", "JCHAIN")
common_top10 = intersect(top10_NL, top10_LS)
common_top10

#VI. Linear Dimension Reduction

husk <- RunPCA(husk, features = VariableFeatures(object = husk))
print(husk[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(husk, dims = 1:2, reduction = "pca")
DimPlot(husk, reduction = "pca")

huskNL <- RunPCA(huskNL, features = VariableFeatures(object = huskNL))
print(huskNL[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(huskNL, dims = 1:2, reduction = "pca")
DimPlot(huskNL, reduction = "pca")
#DimPlot(huskNL, reduction = "pca",label = TRUE)

huskLS <- RunPCA(huskLS, features = VariableFeatures(object = huskLS))
print(huskLS[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(huskLS, dims = 1:2, reduction = "pca")
DimPlot(huskLS, reduction = "pca")

##DimHeatmap(husk, dims = 1, cells = 500, balanced = TRUE) #heatmap: exploration of the primary sources of heterogeneity
##DimHeatmap(husk, dims = 1:15, cells = 500, balanced = TRUE)
#metafeature: combines information across a correlated feature set

#VII. Determination of Dimension
husk <- JackStraw(husk, num.replicate = 100)
husk <- ScoreJackStraw(husk, dims = 1:20)
JackStrawPlot(husk, dims = 1:15) #accroding to graph: good p-values
ElbowPlot(husk) #if taking the first signf. slope change it would remain only 3PCs -> taking 10 or 15 instead?

# husk <- JackStraw(husk, num.replicate = 100)
# husk <- ScoreJackStraw(husk, dims = 1:20)
# JackStrawPlot(husk, dims = 1:15) #accroding to graph: good p-values
# ElbowPlot(husk) #if taking the first signf. slope change it would remain only 3PCs -> taking 10 or 15 instead?

#VIII. Clustering Cells
husk <- FindNeighbors(husk, dims = 1:10)
husk <- FindClusters(husk, resolution = 0.7) #0.5, 0.4, 0.3
#husk <- FindClusters(husk, resolution = 0.4)
#husk <- FindClusters(husk, resolution = 0.3)
head(Idents(husk), 5)



#IX. Non-linear Dimension Reduction
husk <- RunUMAP(husk, dims = 1:10)
#DimPlot(husk, reduction = "umap")
DimPlot(husk, reduction = "umap", label = TRUE)
#saveRDS(husk, file = "/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/CellChatSelf/husk/output/huskOutput.rds")


#X. Cluster Biomarkers
cluster0.markers <- FindMarkers(husk, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
VlnPlot(husk, features = c("CLEC3B", "FBN1", "FSTL1", "IGFBP5", "MFAP5"))

cluster1.markers <- FindMarkers(husk, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
VlnPlot(husk, features = c("APOE", "CXCL12", "IGFBP7", "APOC1", "MGST1"))

cluster2.markers <- FindMarkers(husk, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
VlnPlot(husk, features = c("B2M", "CD3D", "CD52", "GSN", "IL32"))


cluster3.markers <- FindMarkers(husk, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
VlnPlot(husk, features = c("CFD", "APOD", "GSN", "DCN", "CXCL14"))

cluster4.markers <- FindMarkers(husk, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)
VlnPlot(husk, features = c("CD69", "FOS", "JUNB", "CD52", "DNAJB1"))

cluster5.markers <- FindMarkers(husk, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
VlnPlot(husk, features = c("COL6A5", "TNC", "APCDD1", "COL18A1", "CFD"))

#husk.markers <- FindAllMarkers(husk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)  #No features pass logfc.threshold threshold, maybe threshold=0.1, 0.15, 0.2?
husk.markers <- FindAllMarkers(husk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2) #0.2 is okki
husk.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#cluster0.markers <- FindMarkers(husk, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#cluster0.markers <- FindMarkers(husk, ident.1 = 0, logfc.threshold = 0.2, test.use = "roc", only.pos = TRUE)


#VlnPlot(husk, features = c("Col1a1", "Col1a2"), slot = "counts", log = TRUE) #plot raw counts

#FeaturePlot(husk, features = c("Sparc", "Tyrobp", "Igfbp7", "Pecam1", "Cd34", "Mfap4", "Ndufa4l2", "Il6",
#"Birc5"))

#FeaturePlot(husk, features = c("SLPI", "MFAP5", "APOE", "APOC1", "PTPRC", "CD3D", "APOD", "APOE",
#                                "CD69", "FOS", "TNC", "COL6A5", "LYZ", "HLA-DRA", "CCL1", "GZMB", "COL11A1","ASPN"))
FeaturePlot(husk, features = c("CFD", "IBSP", "COCH", "CCL1", "HLA-DRA", "APOD", "MMP12", "CXCL14", "CCL19", "CD74"))

husk.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(husk, features = top10$gene) + NoLegend()

#cluster0.markers <- FindMarkers(husk, ident.1 = 0, min.pct = 0.25)
#head(cluster0.markers, n = 2)
#FeaturePlot(husk, features = c("CLEC3B", "FBN1"))



#XI. Assigning Cell Type Identity to Clusters
new.cluster.ids <- c("0","1","FIB","3","TC","5","6","7","8","9","10")
names(new.cluster.ids) <- levels(husk)
husk <- RenameIdents(husk, new.cluster.ids)
DimPlot(husk, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#FindConservedMarkers(husk, idents.1= 0 )
#Idents(husk)
#RenameIdents(husk, '[nr.]'= [celltype] )

