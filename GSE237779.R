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
data<-read_rds('GSE237779_integrated_seurat_object_CD45_scRNAseq.rds')
#----------Isolate CD45+  -------------------------------------
VlnPlot(data, features = c('PTPRC'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey'))

data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
#-----Figure 3a
VlnPlot(data, features = c('MILR1','MRC1','TNF','IL1B','IL6ST','CXCL8'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey'))
FindMarkers(data, ident.1 = 'Mono/Macro', ident.2 = 'CD8+ T', features = c('MILR1','MRC1','TNF','IL1B','IL6ST','CXCL8'))
#----Figure 3b
VlnPlot(data, features = c('IL18','BCL2','BCL2L1'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey'))
FindMarkers(data, ident.1 = 'Mono/Macro', ident.2 = 'CD8+ T', features = c('IL18','BCL2','BCL2L1'))
table(data@meta.data$integrated_celltype)
break 
