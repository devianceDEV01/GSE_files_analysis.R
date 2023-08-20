library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
setwd()
matrix<-read.delim('GSE66261_read_count_table.txt.gz')
rownames(matrix)<-make.names(matrix$Name, unique = TRUE)
matrix<-matrix[,-c(1:2,5,8,11,14)]
colnames(matrix)<-c('IL12+_1','IL4+_1',
                     'IL12+_2','IL4+_2',
                     'IL12+_3','IL4+_3',
                     'IL12+_4','IL4+_4')
data <- CreateSeuratObject(counts = matrix)
head(data@meta.data)
gc()
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
gc()
data <- NormalizeData(data)
gc()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
gc()
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
gc()
all.genes <- rownames(data)
all.genes
gc()
data <- ScaleData(data, features = all.genes)
dim(data)
gc()
data <- RunPCA(data, npcs=7,features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:7, cells = 500, balanced = T)
ElbowPlot(data)
DimPlot(data, reduction = "pca", pt.size = 7, label = FALSE)
gc()
#Figure 2e
VlnPlot(data,features = c('TBX21','PDCD1','IFNG','TNF','DDX58','IFIH1','DHX58','SEC14L1','NFKB1'),idents = c(),cols = c('orange','grey'))
FindMarkers(data, ident.1 = 'IL12+', ident.2 = 'IL4+', features = c('TBX21','PDCD1','IFNG','TNF','DDX58','IFIH1','DHX58','SEC14L1','NFKB1'))
break 
