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
setwd('R4_N')
features_path <- 'genes.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'peripheral')
summary(x@active.ident)
#---------------------------------------------------------------------------
setwd('R4_T')
features_path <- 'genes.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'GBM')
summary(y@active.ident)
#---------------------------------------------------------------------------
data<-merge(x,y=c(y),project='')
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
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
rm(matrix)
gc()
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD45.groups')
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19-  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD19.groups')
head(data@meta.data)
data <- subset(data, subset = CD19.groups != "CD19.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
head(data@meta.data)
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- CD14+     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD14.groups')
head(data@meta.data)
data <- subset(data, subset = CD14.groups != "CD14.neg")
gc()
table(data@meta.data$orig.ident)
#----------Isolate MILR1+  -------------------------------------
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
DimPlot(data, reduction = 'pca',split.by = 'MILR1.groups')
head(data@meta.data)
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
table(data@meta.data$orig.ident)
VlnPlot(data, features = c('MILR1'),cols = c())
#----split data
set.seed(13)
reg<-FetchData(data,vars = c('ident','CD163','LYZ','S100A8'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='GBM', 1, 0)
table(reg$ident)
#------Even out group numbers and shuffle
edit<-reg[-c(1:316),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
test<-reg
#----------run model
setwd()
model<-read_rds('gbm_pt2_training.rda')
summary(model)
logLik(model)
#------Displaying variance inflation factors
vif(model)
#------Displaying variable importance factors
varImp(model)
newdata = test
summary(newdata)
dim(newdata)
dim(test)
summary(predict(model, newdata, type = 'response'))
predicted<-predict(model, newdata, type = 'response')
pred_factor <- predicted
pred_factor<- round(pred_factor)
pred_factor<-as.factor(pred_factor)
summary(pred_factor)
table(test$ident)
actual<-as.factor(test$ident)
confusion_matrix<-caret::confusionMatrix(pred_factor, actual,positive='1')
confusion_matrix
table(actual)
table(pred_factor)
actuals<-test$ident
df<-cbind(actuals,predicted)
df<-data.frame(df,check.names = FALSE)
optCutoff<-optimalCutoff(actuals = actuals,
                         predictedScores = predicted,
                         optimiseFor = "Ones",
                         returnDiagnostics = TRUE)
head(optCutoff)
auc(actuals,predicted)
ROC_pred<-prediction(df$predicted,df$actuals)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot.roc(actuals, predicted, percent = TRUE, main = 'Validation_ROC', add =  FALSE, asp = NA, print.auc = TRUE)
summary(model)
logLik(model)
confusion_matrix
table(data@meta.data$orig.ident)
table(test$ident)
break 
#----------- Figure 4g
VlnPlot(data, features = c('CD163','LYZ','S100A8'),cols = c('red','grey'))
FindMarkers(data, ident.1 = 'GBM', ident.2 = 'peripheral', features = c('CD163','LYZ','S100A8'))
#----------- Figure 4h
fourfoldplot(as.table(confusion_matrix),color = c('grey','red'),main='Peripheral=0 GBM=1')
plot.roc(actuals, predicted, percent = TRUE, main = 'GBM_pt4_ROC', add =  FALSE, asp = NA, print.auc = TRUE)
#---------- Figure 4i
VlnPlot(data, features = c('DDX58','IFIH1','NFKB1'),cols = c('red','grey'))
FindMarkers(data, ident.1 = 'GBM', ident.2 = 'peripheral', features = c('DDX58','IFIH1','NFKB1'))

