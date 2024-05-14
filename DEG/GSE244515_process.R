library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(Matrix)
library(stringr)
library(glmGamPoi)
file1 <- list.files('GSE244515_RAW')

file_list <- data.frame(str_split_fixed(file1,'_',3))
file_list_1 <- paste(file_list[,1],file_list[,2],'',sep='_')
file_list_1 <- unique(file_list_1)

file_list_1 <- file_list_1[c(1:11,22:27)]

num1 <- 1
for(i in file_list_1){
  tmp1 <- myRead10X(data.dir = 'GSE244515_RAW/',prefixed = i)
  tmp1 <- CreateSeuratObject(counts = tmp1)
  tmp1$sample <- i
  tmp2 <- paste0('sample',num1)
  assign(tmp2,tmp1)
  num1 <- num1+1
}
obj2 <- merge(sample1,sample2)

for(j in 3:17){
  objname <- paste0('sample',j)
  obj2 <- merge(obj2,get(objname))
}



rm(list = paste0('sample',c(1:27)))
obj2 <- SCTransform(obj2)

meta1 <- obj2@meta.data
obj2$sample_2 <- str_split_fixed(obj2$sample,'_',3)[,2]
obj2 <- subset(obj2,nFeature_RNA >= 150)
obj2 <- RunPCA(obj2)
obj2 <- RunUMAP(obj2,dims=1:10)
obj2 <- FindNeighbors(obj2,dims=1:10)
obj2 <- FindClusters(obj2,resolution = 0.1)
FeaturePlot(obj2,features = c('CD3D','CD3E'))
DimPlot(obj2)

rm(obj2)
obj3 <- subset(obj3,seurat_clusters %in% c(2,5,6,7,8))
obj3 <- SCTransform(obj3)
FeaturePlot(obj3,features = c('CD4','CD8A','CD44'))
DimPlot(obj3)
obj3 <- RunPCA(obj3)
obj3 <- RunUMAP(obj3,dims=1:15)
obj3 <- FindNeighbors(obj3,1:10,reduction = "pca")
obj3 <- FindClusters(obj3,resolution = 1)


obj4 <- subset(obj3,seurat_clusters %in% c(seq(0,2),5,seq(7,9),seq(11,15)))
DimPlot(obj4)
FeaturePlot(obj4,features = c('CD4','CD8A','CD44'))
obj4 <- SCTransform(obj4)
obj4 <- RunPCA(obj4)
obj4 <- RunUMAP(obj4,dims=1:15)
obj4 <- FindNeighbors(obj4,dims=1:10)
obj4 <- FindClusters(obj4,resolution = 0.5)

DimPlot(obj4)
FeaturePlot(obj4,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))

FeaturePlot(obj4,features = c('FOXP3','CD8B','TCF7','SELL'))
FeaturePlot(obj4,features = c('IL7R','CD74','GZMB','FOXP3','LGALS1','VIM','IFITM2','GNLY','NKG7'))
FeaturePlot(obj4,features = c('S100A4','S100A6','CTLA4','LEF1','TCF7','SELL','CCR7'))
obj4 <- subset(obj4,seurat_clusters %in% c(0,3,5,6))

cd4_act <- subset(obj4,seurat_clusters %in% c(0,3,5,6))






FeaturePlot(obj4,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))
DimPlot(obj4)
FeaturePlot(obj4,features = c('IL7R','CD52','GIMAP7','SARAF','BTG1','CXCR4'))
FeaturePlot(obj4,features = c('LEF1','TCF7','SELL','CCR7'))
FeaturePlot(obj4,features = c('IL7R','CD40LG','ANXA1','FOS','JUN'))
FeaturePlot(obj4,features = c('FOXP3','SAT1','IL2RA','CTLA4'))
FeaturePlot(obj4,features = c('GNLY','IL2','NKG7','VIM','S100A4','IFITM2','CD74','HSP90AA1'))
FeaturePlot(obj4,features = c('GNLY','CD69','EOMES','GZMA','CXCL13','CTLA4'))
DimPlot(obj4)
obj4$Diabetes <- 'Healthy'
obj4$Diabetes[obj4$sample_2 %in% c('PDDM1', 'PDDM2', 'PDDM3', 'PDDM4', 'PDDM5', 'PDDM6')] <- 'T2D'

cd4_act <- subset(obj4,seurat_clusters %in% c(3,5,6))
cd4_mem <- subset(obj4,seurat_clusters %in% c(0))
cd4_naive <- subset(obj4,seurat_clusters %in% c(1,2,4))









marker2 <- FindMarkers(cd4_mem,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)
marker3 <- FindMarkers(cd4_naive,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)






obj4$Cell_type <- 'CD4_act'
obj4$Cell_type[obj4$seurat_clusters  %in% c(0) ] <- 'CD4_mem'
obj4$Cell_type[obj4$seurat_clusters  %in% c(1,2,4) ] <- 'CD4_naive'


meta2 <- obj4@meta.data
tmp1 <- obj4@assays$SCT@scale.data

tmp2 <- data.frame(meta2$Cell_type,tmp1['GIMAP4',])
colnames(tmp2) <- c('Cell_type','Gene')










