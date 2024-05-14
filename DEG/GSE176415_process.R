library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
# library(clusterp)
library(pheatmap)
library(Matrix)
library(stringr)
library(glmGamPoi)

filelist <- list.files('GSE176415_RAW')
filelist2 <- data.frame(str_split_fixed(filelist,'_',4))
filelist2 <- paste(filelist2[,1],filelist2[,2],filelist2[,3],'',sep='_')
filelist2 <- unique(filelist2)

sample_list <- str_split_fixed(filelist,'_',4)[,3]
sample_list <- unique(sample_list)


num1 <- 1
for(i in filelist2){
  tmp1 <- myRead10X(data.dir = 'GSE176415_RAW/',prefixed = i)
  tmp1 <- CreateSeuratObject(counts = tmp1)
  tmp1$sample <- sample_list[num1]
  tmp2 <- paste0('sample',num1)
  assign(tmp2,tmp1)
  num1 <- num1+1
}
obj2 <- merge(sample1,sample2)

for(j in 3:7){
  objname <- paste0('sample',j)
  obj2 <- merge(obj2,get(objname))
}


obj2 <- SCTransform(obj2)
meta1 <- obj2@meta.data

obj2 <- subset(obj2,nFeature_RNA >= 150)
obj2 <- RunPCA(obj2)
obj2 <- RunUMAP(obj2,dims=1:10)
obj2 <- FindNeighbors(obj2,dims=1:10)
obj2 <- FindClusters(obj2,resolution = 0.1)
FeaturePlot(obj2,features = c('CD3D','CD3E'))
DimPlot(obj2)


T_obj <- subset(obj2,seurat_clusters==3)
DimPlot(T_obj)

T_obj <- RunPCA(T_obj)
T_obj <- RunUMAP(T_obj,dims=1:10)
T_obj <- FindNeighbors(T_obj,dims=1:10)
T_obj <- FindClusters(T_obj,resolution = 0.1)
DimPlot(T_obj,group.by = 'sample')
DimPlot(T_obj)
FeaturePlot(T_obj,features = c('CD4','CD8A','CD8B','CD44','NKG7','GZMA','CCR7'))

CD4_obj <- subset(T_obj,seurat_clusters %in% c(0,3,4))

CD4_obj <- RunPCA(CD4_obj)
CD4_obj <- RunUMAP(CD4_obj,dims=1:10)
CD4_obj <- FindNeighbors(CD4_obj,dims=1:10)
CD4_obj <- FindClusters(CD4_obj,resolution = 0.5)
DimPlot(CD4_obj,group.by = 'sample')
DimPlot(CD4_obj)
FeaturePlot(CD4_obj,features = c('S100A4','S100A6','CTLA4','LEF1','TCF7','SELL','CCR7','FOXP3'))
FeaturePlot(CD4_obj,features = c('IL7R','CD74','GZMB','FOXP3','LGALS1','VIM','IFITM2','GNLY','NKG7','DDX21','CD69','HSPA8'))

FeaturePlot(CD4_obj,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))

FeaturePlot(CD4_obj,features = c('IL7R','CD52','GIMAP7','SARAF','BTG1','CXCR4'))
FeaturePlot(CD4_obj,features = c('LEF1','TCF7','SELL','CCR7'))
FeaturePlot(CD4_obj,features = c('IL7R','CD40LG','ANXA1','FOS','JUN'))
FeaturePlot(CD4_obj,features = c('FOXP3','SAT1','IL2RA','CTLA4'))
FeaturePlot(CD4_obj,features = c('GNLY','IL2','NKG7','VIM','S100A4','IFITM2','CD74','HSP90AA1'))
FeaturePlot(CD4_obj,features = c('GNLY','CD69','EOMES','GZMA','CXCL13','CTLA4'))

CD4_obj$Diabetes <- 'Healthy'
CD4_obj$Diabetes[cd4_act$sample %in% c('W1','W2','W3')] <- 'Diabetes'


cd4_naive <- subset(CD4_obj,seurat_clusters %in% c(1))
cd4_mem <- subset(CD4_obj,seurat_clusters %in% c(0))



marker2 <- FindMarkers(cd4_mem,ident.1 = 'Diabetes',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)
marker3 <- FindMarkers(cd4_naive,ident.1 = 'Diabetes',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)



