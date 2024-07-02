library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(Matrix)
library(DESeq2)
library(EnhancedVolcano)
library(BiocParallel)
library(Seurat)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
param<-MulticoreParam(workers = 6,
                      progressbar = TRUE)
register(param)
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.2A/R2_N')
counts <- readMM("matrix.mtx.gz")
genes <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
head(gene_ids)
cell_ids <- read_tsv("barcodes.tsv.gz", col_names = FALSE)$X1
head(cell_ids)
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
matrix<-data.frame(counts)
names<-c(paste0("norm_", 1:11883))
colnames(matrix)<-names
setwd('/home/em_b/Desktop/scRNAseq_manuscript/Fig.2A/R2_T')
counts2 <- readMM("matrix.mtx.gz")
genes2 <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids2 <- genes2$X2
head(gene_ids2)
cell_ids2 <- read_tsv("barcodes.tsv.gz", col_names = FALSE)$X1
head(cell_ids2)
rownames(counts2) <- gene_ids2
colnames(counts2) <- cell_ids2
matrix2<-data.frame(counts2)
names2<-c(paste0("tumor_", 1:15094))
colnames(matrix2)<-names2
matrix<-matrix[,c(1:1000)]
matrix2<-matrix2[,c(1:1000)]
Matrix<-cbind(matrix,matrix2)
Matrix<-Matrix+1
rm(counts,counts2,genes,genes2,matrix,matrix2)
gc()
names<-colnames(Matrix)
condition<-c(replicate(1000,'normal'),replicate(1000,'tumor'))
str(condition)
dim(Matrix)
type<-replicate(2000, 'paired')
gc()
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
head(coldata)
colnames(Matrix)[1:5]
row.names(coldata)[1:5]
gc()
deseq<-DESeqDataSetFromMatrix(countData = Matrix,
                              colData = coldata,
                              design = ~ condition)
gc()
deseq
deseq<-estimateSizeFactors(deseq)
DE<-DESeq(deseq,
          test = c('LRT'),
          useT = TRUE,
          fitType = 'glmGamPoi',
          minmu = 1e-6,
          minReplicatesForReplace = Inf,
          reduced = ~ 1,
          parallel = TRUE,
          BPPARAM = param)
resLRT <- results(DE)
resLRT<-data.frame(resLRT)
resLRT
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
results_df<-data.frame(resLRT)
results_df<-results_df[order(results_df$padj, decreasing=FALSE),]
results_df<-subset(results_df,padj< 0.05)
summary(results_df)
#---Figure 2A
EnhancedVolcano(resLRT,
                lab = row.names(resLRT),
                x='log2FoldChange',
                y='pvalue',
                title = 'Differential gene expressions from glioblastoma tumor core vs paired surrounding peripheral tissue',
                subtitle = 'Gamma-poisson quasi-likelihood ratio test | n=2000 single cells',
                legendLabels = NULL,
                legendIconSize = -1,
                legendPosition = 'bottom',
                pCutoff = 0.001,
                FCcutoff = 0.2,
                shape = 1,
                pointSize = 1,
                borderColour = 'skyblue',
                hlineCol ='blue',
                borderWidth = 1,
                ylim = c(0,200),
                xlim = c(-1.5,2.5))
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(resLRT)
head(geneid)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
pathways<-c('GO:0002262','GO:0002274','GO:0097529','GO:0061515','GO:0002573',
            'GO:0030099','GO:0002275','GO:0002761') 
go@result = go@result[go@result$ID %in% pathways,]
go@result
GO_plot <- pairwise_termsim(go)
#---Figure 2B
dotplot(GO_plot,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        label_format=20,
        title='Gene Ontology Biological Process enrichment terms')
