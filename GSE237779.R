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
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
VlnPlot(data,features = 'PTPRC')
Idents(data)<-data@meta.data$integrated_celltype
#----  Figure 1a
VlnPlot(data, features = c('MILR1'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey','red','grey','grey','grey','grey','grey','grey','grey','grey','grey','red'))
FindMarkers(data, ident.1 = 'Macrophage', ident.2 = 'CD4+ Tcm', features = c('MILR1'))
#----- Figure 1b
VlnPlot(data, features = c('MRC1'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey','red','grey','grey','grey','grey','grey','grey','grey','grey','grey','red'))
FindMarkers(data, ident.1 = 'Macrophage', ident.2 = 'CD4+ Tcm', features = c('MRC1'))

