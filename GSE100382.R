library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
setwd()
matrix <- read.delim("GSE100382_RNAseq_gene_expression.txt")
dim(matrix)
colnames(matrix) <- c('gene','IFNa_1','IFNa_2','IFNa_3','IFNa+TNF+LPS_1','IFNa+TNF+LPS_2','IFNa+TNF+LPS_3',
                      'IFNa+TNF_1','IFNa+TNF_2','IFNa+TNF_3','IFNa+LPS_1','IFNa+LPS_2','IFNa+LPS_3',
                      'TNF+LPS_1','TNF+LPS_2','TNF+LPS_3','control_1','control_2','control_3',
                      'TNF_1','TNF_2','TNF_3','LPS_1','LPS_2','LPS_3')
row.names(matrix) <- make.names(matrix$gene, unique = TRUE)
matrix <- matrix[,-1]
data <- CreateSeuratObject(counts = matrix)
data@meta.data
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
data <- RunPCA(data,npcs = 23, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
data
table(data@meta.data$orig.ident)
break 
#---- Figure 3c 
VlnPlot(data,features = c('MILR1'),idents = c(),cols = c('black','grey','grey','grey','grey','grey','red','grey'))
FindMarkers(data, ident.1 = 'TNF', ident.2 = 'control', features = c('MILR1'))
#------ Figure 3d
VlnPlot(data,features = c('MILR1','MRC1','DDX58','IFIH1','DHX58','SEC14L1'),idents = c(),cols = c('black','grey','grey','grey','grey','grey','red','grey'))
FindMarkers(data, ident.1 = 'TNF', ident.2 = 'control', features = c('MILR1','MRC1','DDX58','IFIH1','DHX58','SEC14L1'))
#------ Figure 3e
VlnPlot(data,features = c('TNFAIP3','TNFSF10','IFNA1','IFNB1','IL1B','IL4I1','IL6','IL8','IL10'),idents = c(),cols = c('black','grey','grey','grey','grey','grey','red','grey'))
FindMarkers(data, ident.1 = 'TNF', ident.2 = 'control', features = c('TNFAIP3','TNFSF10','IFNA1','IFNB1','IL1B','IL4I1','IL6','IL8','IL10'))
#-----Figure 3f
VlnPlot(data,features = c('IL18','TNFSF10','CSF1','CSF2','CSF3'),idents = c(),cols = c('black','grey','grey','grey','grey','grey','red','grey'))
FindMarkers(data, ident.1 = 'TNF', ident.2 = 'control', features = c('IL18','TNFSF10','CSF1','CSF2','CSF3'))
#-----figure 3g
VlnPlot(data,features = c('BCL2','BCL2L1','MYC','VEGFA','TGFB2','FAS','CASP8','CASP10'),idents = c(),cols = c('black','grey','grey','grey','grey','grey','red','grey'))
FindMarkers(data, ident.1 = 'TNF', ident.2 = 'control', features = c('BCL2','BCL2L1','MYC','VEGFA','TGFB2','FAS','CASP8','CASP10'))
