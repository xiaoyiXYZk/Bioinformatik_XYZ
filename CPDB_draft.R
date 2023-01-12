#using CellPhoneDB to analyze mosk pw 12d

#1 - Step Extract data from Seurat Object
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
library(dplyr)
library(Seurat)
library(patchwork)
# Basic function to convert mouse to human gene names
data_ <- readRDS("/Users/xiaoyizheng/Downloads/Bioinfo_Einarbeiten/moskOutput.rds")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
allgenes <- rownames(data)
matrix1 <- as.data.frame(data@assays$RNA@data)
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

### If you are using a mouse data, then its needed to convert the gene names to human orthologs
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(alldata@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrix1 <- matrix1[match(genesV2$MGI.symbol,rownames(alldata),nomatch=F),]
matrix1$gene <- genesV2$Gene.stable.ID
#rownames(matrix1) <-  genesV2$Gene.stable.ID

### Subseting the matrix
s1 <- grepl('state1',alldata@meta.data$cond)
s2 <- grepl('state2',alldata@meta.data$cond)
s1[match('gene',colnames(matrix1))] <- TRUE
s2[match('gene',colnames(matrix1))] <- TRUE

## Checking the dimensions
print(dim(matrix1[,s1]))
print(dim(matrix1[,s2]))


## If the cluster names are categorical, you will need to convert it to numerical
alldata@meta.data$sub_cell_type <- as.factor(alldata@meta.data$sub_cell_type)
print(levels(alldata@meta.data$sub_cell_type))
levels(alldata@meta.data$sub_cell_type) <- 1:length(levels(alldata@meta.data$sub_cell_type))
print(1:length(levels(alldata@meta.data$sub_cell_type)))
alldata@meta.data$sub_cell_type <- as.numeric(alldata@meta.data$sub_cell_type)


write.table(matrix1[,s1], 's1_filtered_hcount.csv',row.names=T,sep=',')
write.table(matrix1[,s2], 's2_filtered_hcount.csv',row.names=T,sep=',')
metadata <- data.frame(cells=rownames(alldata@meta.data[grepl('state1',alldata@meta.data$stim),]),cluster=alldata@meta.data$sub_cell_type[grepl('state1',alldata@meta.data$stim)])
metadata_s2 <- data.frame(cells=rownames(alldata@meta.data[!grepl('state1',alldata@meta.data$stim),]),cluster=alldata@meta.data$sub_cell_type[!grepl('state1',alldata@meta.data$stim)]) ## Just negate grepl('state1',alldata@meta.data$stim),]
print('Writing Metadata')
write.csv(metadata, 's1_filtered_meta.csv', row.names=FALSE)
write.csv(metadata_tac, 's2_filtered_meta.csv', row.names=FALSE)
```