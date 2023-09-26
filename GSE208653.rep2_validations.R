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
setwd()
#---------------------------------------------------------------------------
features_path <- 'GSM6360683_N_2.features.tsv.gz'
barcodes_path <- 'GSM6360683_N_2.barcodes.tsv.gz'
matrix_path <- 'GSM6360683_N_2.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
hpv <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpv')
gc()
#---------------------------------------------------------------------------
features_path <- 'GSM6360685_HSIL_2.features.tsv.gz'
barcodes_path <- 'GSM6360685_HSIL_2.barcodes.tsv.gz'
matrix_path <- 'GSM6360685_HSIL_2.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
lesions <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lesions')
gc()
#---------------------------------------------------------------------------
features_path <- 'GSM6360687_SCC_5.features.tsv.gz'
barcodes_path <- 'GSM6360687_SCC_5.barcodes.tsv.gz'
matrix_path <- 'GSM6360687_SCC_5.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
cervical_cancer <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'cervical_cancer')
gc()
#-------------------------------------------------------------------------------
data<-merge(hpv,y=c(lesions,cervical_cancer),project='hpv_cervical_cancer')
table(data@meta.data$orig.ident)
rm(hpv,lesions,cervical_cancer,matrix)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt <20)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
gc()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
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
table(data@meta.data$orig.ident)
VlnPlot(data, features = c('PTPRC','TRAC','CD19','CD14','MILR1'),cols = c())
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
VlnPlot(data, features = c('PTPRC','TRAC','CD19','CD14','MILR1'),cols = c())
#----------------------------------------------------------------
reg<-FetchData(data,vars = c('ident','CXCL9','IFI6','SLC40A1'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='cervical_cancer', 1, 0)
table(reg$ident)
edit<-reg[-c(1:21),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
test<-reg
#----------run model
setwd()
model<-read_rds('cc_model.rda')
summary(model)
#----------------------------
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
#------Figure ?
VlnPlot(data, features = c('CXCL9','IFI6','SLC40A1'),cols = c('green','grey','grey'))
FindMarkers(data, ident.1 = 'cervical_cancer', ident.2 = 'hpv', features = c('CXCL9','IFI6','SLC40A1'))
#------Figure ?
fourfoldplot(as.table(confusion_matrix),color = c('grey','green'),main='non-tumor=0 cervical_cancer=1')
plot.roc(actuals, predicted, percent = TRUE, main = 'Validation_ROC', add =  FALSE, asp = NA, print.auc = TRUE)
#------Figure ?
VlnPlot(data, features = c('DDX58','IFIH1','NFKB1'),cols = c('green','grey','grey'))
FindMarkers(data, ident.1 = 'cervical_cancer', ident.2 = 'hpv', features = c('DDX58','IFIH1','NFKB1'))
