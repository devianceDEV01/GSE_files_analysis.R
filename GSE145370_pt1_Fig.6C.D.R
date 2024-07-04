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
#---------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S133B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent1@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S134B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent2@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S135B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent3 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent3@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S149B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent4 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent4@meta.data$orig.ident)
#---------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S150B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent5 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent5@meta.data$orig.ident)
#---------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S158B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent6 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent6@meta.data$orig.ident)
#---------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S159B_adjacent/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
adjacent7 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'adjacent')
table(adjacent7@meta.data$orig.ident)
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.6/GSE145370_RAW/S133A_ESCC/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
escc1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'ESCC')
table(escc1@meta.data$orig.ident)
#-------------------------------------------------------------------------------
data<-merge(adjacent1,y=c(adjacent2,adjacent3,adjacent4,adjacent5,
                          adjacent6,adjacent7,escc1),project='ESCC')
summary(data@active.ident)
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
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19-  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
head(data@meta.data)
data <- subset(data, subset = CD19.groups != "CD19.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
head(data@meta.data)
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- CD14+     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
head(data@meta.data)
data <- subset(data, subset = CD14.groups != "CD14.neg")
gc()
table(data@meta.data$orig.ident)
#----------Isolate MILR1+  -------------------------------------
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
head(data@meta.data)
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
table(data@meta.data$orig.ident)
VlnPlot(data, features = c('PTPRC','CD19','TRAC','CD14','MILR1'),pt.size=0.1)
data<-JoinLayers(data)
#--------------------------------------------------------------
set.seed(10)
#----split data
reg<-FetchData(data,vars = c('ident','HSPA1B','MMP12',
                             'MT-ATP8','CTSC',
                             'IER5','FTL',
                             'ACP5','MMP9',
                             'ID3','GSN','CLEC11A','RNASE1',
                             'CTSD','C15orf48'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='ESCC', 1, 0)
table(reg$ident)
#------Even out group numbers and shuffle
edit<-reg[-c(575:801),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.5,0.5))
train  <- reg[sample, ]
table(train$ident)
test   <- reg[!sample, ]
table(test$ident)
#-------------------------------------------------------------------------
model<-glm(ident~HSPA1B+MMP12+
             `MT-ATP8`+CTSC+
             IER5+FTL+
             ACP5+MMP9+
             ID3+GSN+CLEC11A+RNASE1+
             CTSD+C15orf48,
           data = train, family = binomial)
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
#----Figure 6C
VlnPlot(data, features = c('HSPA1B','MMP12',
                           'MT-ATP8','CTSC',
                           'IER5','FTL',
                           'ACP5','MMP9',
                           'ID3','GSN','CLEC11A','RNASE1',
                           'CTSD','C15orf48'),cols = c('grey','yellow'))
#----Figure 6D
fourfoldplot(as.table(confusion_matrix),
             color = c('skyblue','yellow'),
             std='all.max',
             main='Adjacent=0 ESCC=1')
plot.roc(actuals, predicted,
         percent = TRUE,
         main = 'ESCC patient 1 Area Under the ROC curve',
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
