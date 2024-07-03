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
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.4/GSE155698_RAW/AdjNorm_TISSUE_1/filtered_feature_bc_matrix')
features_path <- 'features.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
summary(x@active.ident)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.4/GSE155698_RAW/PDAC_TISSUE_1/filtered_feature_bc_matrix')
features_path <- 'features.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'PDAC')
summary(y@active.ident)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.4/GSE155698_RAW/AdjNorm_TISSUE_2/filtered_feature_bc_matrix')
features_path <- 'features.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
summary(x2@active.ident)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.4/GSE155698_RAW/PDAC_TISSUE_2/filtered_feature_bc_matrix')
features_path <- 'features.tsv.gz'
barcodes_path <- 'barcodes.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'PDAC')
summary(y2@active.ident)
data<-merge(x,y=c(x2,y,y2),project='PDA')
table(data@meta.data$orig.ident)
head(data@active.ident)
rm(x,y,x2,y2)
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
VlnPlot(data, features = c('PTPRC','TRAC','CD19','CD14','MILR1'),cols = c())
data<-JoinLayers(data)
reg<-FetchData(data,vars = c('ident','CD36','HLA-DRB5','INHBA'),layer = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='PDAC', 1, 0)
table(reg$ident)
edit<-reg[-c(1:61),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
set.seed(14)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.5,0.5))
train  <- reg[sample, ]
table(reg$ident)
test   <- reg[!sample, ]
table(reg$ident)
model<-glm(ident~CD36+`HLA-DRB5`+INHBA,data = train, family = binomial)
summary(model)
logLik(model)
vif(model)
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
#----Figure 4c
VlnPlot(data, features = c('CD36','HLA-DRB5','INHBA'),cols = c('grey','green'))
FindMarkers(data, ident.1 = 'PDAC', ident.2 = 'adjacent', features = c('CD36','HLA-DRB5','INHBA'),
            test.use='bimod')
#----Figure 4d
fourfoldplot(as.table(confusion_matrix),
             color = c('yellow','green'),
             std='all.max',
             main='Adjacent=0 PDA=1')
plot.roc(actuals, predicted,
         percent = TRUE,
         main = 'GBM patient 3 Area Under the ROC curve',
         add =  FALSE,
         asp = NA,
         print.auc = TRUE,
         print.auc.col='red',
         print.auc.cex=1,
         grid=TRUE,
         grid.col='skyblue',
         identity.col='blue',
         identity.lty=8,
         col='red',
         print.thres=TRUE,
         print.thres.pch=5,
         print.thres.col='black',
         ci=TRUE)
