library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(Matrix)
library(stringr)
library(ggsignif)
library(ggpubr)
library(viridis)
mytheme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent",colour=NA),
    plot.background = element_rect(fill = "transparent",colour=NA),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6,colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.background = element_rect(fill = "transparent",colour=NA),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
#===============================================================================

file_list <- list.files('GSE165816/')
file_list <- file_list[-1]
file_list2 <- file_list[c(23:26,29:31,34,41,51)]
disease1 <- c('Diabetic subject without DFU','Healthy non-diabetic subject',
            'subject with non-healing DFU','subject with healing DFU',
            'subject with healing DFU','Healthy non-diabetic subject',
            'subject with healing DFU','Healthy non-diabetic subject',
            'subject with non-healing DFU','Diabetic subject without DFU')

for(i in 1:10){
  tmp1 <- read.csv(paste0('GSE165816/',file_list2[i]))
  tmp1 <- CreateSeuratObject(counts = tmp1)
  tmp1$sample <- paste0('sample',i)
  tmp1$disease <- disease1[i]
  tmp2 <- paste0('sample',i)
  assign(tmp2,tmp1)
}


obj3 <- merge(sample1,sample2)
meta1 <- obj3@meta.data
for(j in 3:10){
  objname <- paste0('sample',j)
  obj3 <- merge(obj3,get(objname))
}

rm(list=paste0('sample',c(1:10)))
obj3 <- SCTransform(obj3)
obj3 <- RunPCA(obj3)
obj3 <- RunUMAP(obj3,dims=1:15)
obj3 <- FindNeighbors(obj3,dims = 1:10)
obj3 <- FindClusters(obj3,resolution = 0.1)
DimPlot(obj3)
FeaturePlot(obj3,features = c('CD3D','CD3E'))

T_obj <- subset(obj3,seurat_clusters %in% c(0,2,4))
T_obj <- RunPCA(T_obj)
T_obj <- RunUMAP(T_obj,dims = 1:15)
T_obj <- FindNeighbors(T_obj,dims = 1:15)
T_obj <- FindClusters(T_obj,resolution = 0.1)
DimPlot(T_obj)


FeaturePlot(T_obj,features = c('CD4','CD8A','CD44'))

marker1 <- FindMarkers(T_obj,ident.1 = 0,min.pct = 0,logfc.threshold = 0)

MR_gene_list <- read.csv('T2D gene list.csv')

marker2 <- marker1[tmp1,]





CD4_obj <- subset(T_obj,seurat_clusters == 0)
CD4_obj <- RunPCA(CD4_obj)
CD4_obj <- RunUMAP(CD4_obj,dims = 1:15)
CD4_obj$group <- 'Diabetic'
CD4_obj$group[CD4_obj$disease == 'Healthy non-diabetic subject'] <- 'Health'
table(CD4_obj$group)
Marker_CD4_1 <- FindMarkers(CD4_obj,group.by = 'group',ident.1 = 'Diabetic',min.pct = 0,logfc.threshold = 0)
marker2 <- Marker_CD4_1[tmp1,]
marker3 <- subset(marker2,p_val < 0.05)




#===============================================================================



DimPlot(CD4_obj,group.by = 'group')
CD4_obj <- subset(CD4_obj,UMAP_2 > -5)
CD4_obj <- RunPCA(CD4_obj)
CD4_obj <- RunUMAP(CD4_obj,dims = 1:20)
CD4_obj <- FindNeighbors(CD4_obj,dims = 1:15)
CD4_obj <- FindClusters(CD4_obj,resolution = 0.5)

DimPlot(CD4_obj)

FeaturePlot(CD4_obj,features = c('IL7R','CD74','GZMB','FOXP3','LGALS1','VIM','IFITM2','GNLY','NKG7'))#'CCR7','IL7R','CD74','ANXA1','PTPRC'
FeaturePlot(CD4_obj,features = c('IL7R','CD52','GIMAP7','SARAF','BTG1','CXCR4'))
FeaturePlot(CD4_obj,features = c('LEF1','TCF7','SELL','CCR7','GNLY'))
FeaturePlot(CD4_obj,features = c('IL7R','CD40LG','ANXA1','FOS','JUN'))
FeaturePlot(CD4_obj,features = c('FOXP3','SAT1','IL2RA','CTLA4'))
DimPlot(CD4_obj,group.by = 'group')




CD4_naive <- subset(CD4_obj,seurat_clusters %in% c(0,1,5))
CD4_mem <- subset(CD4_obj,seurat_clusters %in% c(3))



marker_5 <- FindMarkers(CD4_mem,group.by = 'group',ident.1 = 'Diabetic',logfc.threshold = 0,min.pct = 0.1)
marker_6 <- FindMarkers(CD4_naive,group.by = 'group',ident.1 = 'Diabetic',logfc.threshold = 0,min.pct = 0.1)






