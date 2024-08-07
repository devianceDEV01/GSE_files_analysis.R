library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.2A/R2_N')
features_path <- 'genes.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'peripheral')
summary(x@active.ident)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.2A/R2_T')
features_path <- 'genes.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'GBM')
summary(y@active.ident)
data<-merge(x,y=c(y),project='pt2.trainer')
table(data@meta.data$orig.ident)
head(data@active.ident)
rm(x,y)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt <25)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
gc()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
table(data@meta.data$orig.ident)
rm(matrix)
gc()
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
table(data@meta.data$orig.ident)
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
data <- subset(data, subset = CD19.groups != "CD19.pos")
gc()
table(data@meta.data$orig.ident)
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
table(data@meta.data$orig.ident)
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
data <- subset(data, subset = CD14.groups != "CD14.neg")
gc()
table(data@meta.data$orig.ident)
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
table(data@meta.data$orig.ident)
VlnPlot(data, features = c('MILR1','PTPRC','CD19','TRAC','CD14'),cols = c())
reg<-FetchData(data,vars = c('ident','CD163','LYZ','S100A8'),layer = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='GBM', 1, 0)
table(reg$ident)
edit<-reg[-c(8351:9062),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
train<-reg
set.seed(13)
model<-glm(ident~CD163+LYZ+S100A8,data = train, family = binomial)
CD163_estimate<-0.67263
CD163_StdError<-0.02027
zvalue<-CD163_estimate/CD163_StdError
wald_value<-pnorm(zvalue)
CD163_pvalue<-(1-wald_value)*2
CD163_pvalue
LYZ_estimate<-0.16445
LYZ_StdError<-0.01241
zvalue<-LYZ_estimate/LYZ_StdError
wald_value<-pnorm(zvalue)
LYZ_pvalue<-(1-wald_value)*2
LYZ_pvalue
S100A8_estimate<-0.16004
S100A8_StdError<-0.01728
zvalue<-S100A8_estimate/S100A8_StdError
wald_value<-pnorm(zvalue)
S100A8_pvalue<-(1-wald_value)*2
S100A8_pvalue
summary(model)
logLik(model)
vif(model)
varImp(model)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/models')
write_rds(model,file = 'GSE162631_gbm_pt2_training.rda')
