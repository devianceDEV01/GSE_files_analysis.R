library(Seurat)
library(readr)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.A.B')
data<-read_rds('integrated_seurat_object_CD45_scRNAseq.rds')
VlnPlot(data,features = 'PTPRC')
#----  Figure 1a
VlnPlot(data, features = c('MILR1'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey'))
FindMarkers(data, ident.1 = 'Mono/Macro', ident.2 = 'CD4+ T', features = c('MILR1'),
            test.use='negbinom')
#----- Figure 1b
VlnPlot(data, features = c('MRC1'),pt.size=0.1,cols = c('red','grey','grey','grey','grey','grey','grey','grey'))
FindMarkers(data, ident.1 = 'Mono/Macro', ident.2 = 'CD4+ T', features = c('MRC1'),
            test.use='negbinom')

