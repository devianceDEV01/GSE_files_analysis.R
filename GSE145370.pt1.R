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
setwd('P1.A')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
a1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(a1@meta.data$orig.ident)
#---------------------------------------------------------------------------
setwd('P2.A')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
a2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(a2@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('P3.A')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
a3 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(a3@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('/home/amp_prog/Desktop/drive_sept2023/models/GSE145370_ESCC/P4.A')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
a4 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(a4@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('P5.A')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
a5 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(a5@meta.data$orig.ident)
#---------------------------------------------------------------------------
setwd('P1.T')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
e1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'ESCC')
table(e1@meta.data$orig.ident)
#-------------------------------------------------------------------------------
data<-merge(a1,y=c(a2,a3,a4,a5,e1),project='ESCC')
rm(a1,a2,a3,a4,a5,e1,matrix,features_path,barcodes_path,matrix_path)
head(data@active.ident)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt <20)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
VlnPlot(data, features = c('PTPRC','CD19','TRAC','CD14','MILR1'),pt.size=0.1)
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
#--------------------------------------------------------------
set.seed(10)
#----split data
reg<-FetchData(data,vars = c('ident','C1QC','PMP22','CLEC5A'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='ESCC', 1, 0)
table(reg$ident)
#------Even out group numbers and shuffle
edit<-reg[-c(433:730),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.5,0.5))
train  <- reg[sample, ]
table(train$ident)
test   <- reg[!sample, ]
table(test$ident)
#-------------------------------------------------------------------------
model<-glm(ident~C1QC+PMP22+CLEC5A,data = train, family = binomial)
summary(model)
logLik(model)
#------Displaying variance inflation factors
vif(model)
#------Displaying variable importances
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
plot.roc(actuals, predicted, percent = TRUE, main = 'Validation_ROC', add =  FALSE, asp = NA, print.auc = TRUE)
summary(model)
logLik(model)
confusion_matrix
table(data@meta.data$orig.ident)
table(test$ident)
break 
#----Figure 6a
VlnPlot(data, features = c('C1QC','PMP22','CLEC5A'),cols = c('grey','yellow'))
FindMarkers(data, ident.1 = 'ESCC', ident.2 = 'adjacent', features = c('C1QC','PMP22','CLEC5A'))
#----Figure 6b
fourfoldplot(as.table(confusion_matrix),color = c('grey','yellow'),main='adjacent=0 ESCC=1')
plot.roc(actuals, predicted, percent = TRUE, main = 'ESCC_validation_ROC', add =  FALSE, asp = NA, print.auc = TRUE)
