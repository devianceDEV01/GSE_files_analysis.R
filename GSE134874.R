library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
#-------------------------------------------------------------------------------
setwd('/home/amp_prog/Desktop/TAM_manuscript/datasets/GSE134874_PRR_stims')
matrix <- read.delim('GSE134874_Human_umi_counts.tsv')
row.names(matrix) <- make.names(matrix$X,unique = TRUE)
matrix <-matrix[,-c(1)]
data <- CreateSeuratObject(counts = matrix)
#---------------------------------------------------------------------
data@meta.data
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top200 <- head(VariableFeatures(data), 200)
top200
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data,npcs=50,features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
head(data@meta.data)
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
#Figure 2c
VlnPlot(data,features = c('IFNA1','IFNB1','IL12A','IL12B','DDX58','IFIH1','DHX58','NFKB1','MILR1'),idents = c('Ctrl','LS'),cols = c('grey','skyblue'))
FindMarkers(data, ident.1 = 'LS', ident.2 = 'Ctrl', features = c('IFNA1','IFNB1','IL12A','IL12B','DDX58','IFIH1','DHX58','NFKB1','MILR1'))

