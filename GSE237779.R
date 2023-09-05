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
setwd('/home/amp_prog/Desktop/TAM_manuscript/datasets/GSE237779_immune.cells_GBM')
data<-read_rds('GSE237779_integrated_seurat_object_CD45_scRNAseq.rds')
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
VlnPlot(data,features = 'PTPRC')
#-----regression modeling
reg<-FetchData(data,vars = c('ident','CD163','FPR3','S100A8','CLEC5A','PLAU','IFIT3','MILR1','CD14'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='Mono/Macro', 1, 0)
table(reg$ident)
#------Even out group numbers and shuffle
edit<-reg[-c(1:16356),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
set.seed(13)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.5,0.5))
train  <- reg[sample, ]
table(train$ident)
test   <- reg[!sample, ]
table(test$ident)
#----------run model
model<-glm(ident~FPR3+CLEC5A+PLAU+IFIT3+MILR1+CD14,
           data = train, family = binomial)
summary(model)
logLik(model)
#----------Mcfadden's pseudo R squared
null<-model$null.deviance/-2
resdDEV<-model$deviance/-2
pR2<-(null-resdDEV)/null
print(pR2)
#------Displaying variance inflation factors
vif(model)
#------Displaying variable importance factors
varImp(model)
newdata = data.frame(test)
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
caret::confusionMatrix(pred_factor, actual,positive='1')
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
ROC_pred<-prediction(df$predicted,df$actual)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
summary(model)
logLik(model)
caret::confusionMatrix(pred_factor, actual,positive='1')
table(data@meta.data$integrated_celltype)
table(test$ident)
cat('Gene panel is highest in mono/macrophage subsets')
break 
Idents(data)<-data@meta.data$integrated_celltype
VlnPlot(data, features = c('MILR1'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey','red','grey','grey','grey','grey','grey','grey','grey','grey','grey','red'))
FindMarkers(data, ident.1 = 'Macrophage', ident.2 = 'CD4+ Tcm', features = c('MILR1'))

