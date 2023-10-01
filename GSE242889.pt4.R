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
setwd('4NT_P29')
features_path <- 'genes.tsv'
barcodes_path <- 'barcodes.tsv'
matrix_path <- 'matrix.mtx'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
summary(x@active.ident)
#---------------------------------------------------------------------------
setwd('4T_C29')
features_path <- 'genes.tsv'
barcodes_path <- 'barcodes.tsv'
matrix_path <- 'matrix.mtx'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'HCC')
summary(y@active.ident)
#---------------------------------------------------------------------------
data<-merge(x,y=c(y),project='')
table(data@meta.data$orig.ident)
head(data@active.ident)
rm(x,y)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <25)
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
VlnPlot(data, features = c('PTPRC','CD19','TRAC','CD14','MILR1'),cols = c('red','grey'))
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19-  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
data <- subset(data, subset = CD19.groups != "CD19.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- CD14+     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
data <- subset(data, subset = CD14.groups != "CD14.neg")
gc()
table(data@meta.data$orig.ident)
#----------Isolate MILR1+  -------------------------------------
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
table(data@meta.data$orig.ident)
VlnPlot(data, features = c('MILR1'),cols = c())
table(data@meta.data$orig.ident)
#----split data
reg<-FetchData(data,vars = c('ident','HP','GPNMB','HRG'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='HCC', 1, 0)
table(reg$ident)
#------Even out group numbers and shuffle
edit<-reg[-c(1:212),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
test<-reg
#--------------------------------------------------
set.seed(14)
#----------run model
setwd()
model<-read_rds('hcc.pt5_model.rda')
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
plot.roc(actuals, predicted, percent = TRUE, main = 'Validation_ROC', add =  FALSE, asp = NA, print.auc = TRUE)
summary(model)
logLik(model)
confusion_matrix
table(data@meta.data$orig.ident)
table(test$ident)
break 
#-------Figure 3g
VlnPlot(data, features = c('HP','GPNMB','HRG'),cols = c('grey','skyblue'))
FindMarkers(data, ident.1 = 'HCC', ident.2 = 'adjacent', features = c('HP','GPNMB','HRG'))
#-------Figure 3h
fourfoldplot(as.table(confusion_matrix),color = c('grey','skyblue'),main='Adjacent=0 HCC=1')
plot.roc(actuals, predicted, percent = TRUE, main = 'HCC_pt4__ROC', add =  FALSE, asp = NA, print.auc = TRUE)
