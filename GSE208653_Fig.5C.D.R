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
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.5/GSE208653_RAW')
features_path <- 'GSM6360682_N_1.features.tsv.gz'
barcodes_path <- 'GSM6360682_N_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360682_N_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
hpv1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpvNorm')
gc()
features_path <- 'GSM6360686_SCC_4.features.tsv.gz'
barcodes_path <- 'GSM6360686_SCC_4.barcodes.tsv.gz'
matrix_path <- 'GSM6360686_SCC_4.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
cervical_cancer1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpvCC')
gc()
features_path <- 'GSM6360683_N_2.features.tsv.gz'
barcodes_path <- 'GSM6360683_N_2.barcodes.tsv.gz'
matrix_path <- 'GSM6360683_N_2.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
hpv2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpvNorm')
gc()
features_path <- 'GSM6360687_SCC_5.features.tsv.gz'
barcodes_path <- 'GSM6360687_SCC_5.barcodes.tsv.gz'
matrix_path <- 'GSM6360687_SCC_5.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
cervical_cancer2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpvCC')
gc()
data<-merge(hpv1,y=c(cervical_cancer1,hpv2,cervical_cancer2),project='hpv-CC')
table(data@meta.data$orig.ident)
rm(hpv1,cervical_cancer1,hpv2,cervical_cancer2,matrix)
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
table(data@meta.data$orig.ident)
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
VlnPlot(data, features = c('PTPRC','TRAC','CD19','CD14','MILR1'),cols = c())
data<-JoinLayers(data)
set.seed(14)
reg<-FetchData(data,vars = c('ident','CXCL10','FN1','SDC3'),layer = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='hpvCC', 1, 0)
table(reg$ident)
edit<-reg[-c(504:603),]
edit<-edit[-c(228:300),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.5,0.5))
train  <- reg[sample, ]
table(train$ident)
test   <- reg[!sample, ]
table(test$ident)
model<-glm(ident~CXCL10+FN1+SDC3,data = train, family = binomial)
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
#----Figure 5c
VlnPlot(data, features = c('FN1','SDC3','CXCL10'),cols = c('grey','orange'))
FindMarkers(data, ident.1 = 'hpvCC', ident.2 = 'hpvNorm', features = c('FN1','SDC3','CXCL10'),
            test.use='bimod')
#----Figure 5d
fourfoldplot(as.table(confusion_matrix),
             color = c('brown','orange'),
             std='all.max',
             main='hpvNorm=0 hpvCC=1')
plot.roc(actuals, predicted,
         percent = TRUE,
         main = 'hpvNorm v hpvCC Area Under the ROC curve',
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

