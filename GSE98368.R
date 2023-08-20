library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
setwd()
matrix1 <- read.delim('GSE98368_AdGFP_AdMAF_RPKM_table.txt')
colnames(matrix1)[1]<-'Gene'
matrix2<-read.delim('GSE98368_Resting_IFNg_RPKM_table.txt')
matrix<-merge(matrix1,matrix2,by='Gene')
colnames(matrix)<-c('Gene','resting_1','resting_2','IFNG+_1','IFNG+_2','resting_3','resting_4','IFNG+_3','IFNG+_4','resting_5','IFNG+_5','resting_6','IFNG+_6')
row.names(matrix) <- make.names(matrix$Gene,unique = TRUE)
matrix <- matrix[,-1]
data <- CreateSeuratObject(counts = matrix)
data@active.ident
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data,npcs = 11, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:3, cells = 500, balanced = T)
ElbowPlot(data)
data
DimPlot(data, reduction = "pca", pt.size = 7)
#------------Figure 2h
VlnPlot(data, features = c('MRC1','DDX58','IFIH1','DHX58','SEC14L1','CD274','PDCD1LG2','IL1B','IL6','IL6ST','NFKB1'), cols = c('red','grey'))
FindMarkers(data, ident.1 = 'IFNG+', ident.2 = 'resting', features = c('MRC1','DDX58','IFIH1','DHX58','SEC14L1','CD274','PDCD1LG2','IL1B','IL6','IL6ST','NFKB1'))
break
