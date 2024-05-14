library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)

data1 <- readRDS('149c33da-2b25-4072-aecf-0ccea9a8f9e6.rds')
data2 <- readRDS('a7445ac6-93ca-43af-810a-3809c9e4d82e.rds')
meta1 <- obj1@meta.data


obj1 <- subset(data1,subclass.l1=='Immune')
obj1 <- SCTransform(obj1)
obj1 <- RunPCA(obj1)
obj1 <- RunUMAP(obj1,dims=1:20)
obj1 <- FindNeighbors(obj1,dims=1:15)
obj1 <- FindClusters(obj1,resolution = 0.1)
DimPlot(obj1)
FeaturePlot(obj1,features = c('ENSG00000010610','ENSG00000153563'))

gene_name1 <- rownames(obj1)
gene.df <- bitr(gene_name1,fromType = 'ENSEMBL',
                toType = c('SYMBOL'),
                OrgDb = org.Hs.eg.db,drop = T)

gene.df <- gene.df[!duplicated(gene.df$SYMBOL),]
gene.df <- gene.df[!duplicated(gene.df$ENSEMBL),]
obj2 <- obj1[rownames(obj1) %in% gene.df$ENSEMBL,]
identical(rownames(obj2),gene.df$ENSEMBL)


T_data1 <- obj2@assays$RNA$counts
identical(rownames(T_data1),gene.df$ENSEMBL)
rownames(T_data1) <- gene.df$SYMBOL
obj3 <- CreateSeuratObject(counts = T_data1,meta.data = obj2@meta.data)
obj3 <- SCTransform(obj3)
obj3 <- RunPCA(obj3)
obj3 <- RunUMAP(obj3,dims=1:20)
obj3 <- FindNeighbors(obj3,dims=1:15)
obj3 <- FindClusters(obj3,resolution = 0.1)
DimPlot(obj3)
FeaturePlot(obj3,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))

obj4 <- subset(obj3,seurat_clusters %in% c(0,1,2,4,6,7))
DimPlot(obj4)
obj4 <- SCTransform(obj4)
obj4 <- RunPCA(obj4)
obj4 <- RunUMAP(obj4,dims=1:20)
obj4 <- FindNeighbors(obj4,dims=1:15)
obj4 <- FindClusters(obj4,resolution = 0.1)
DimPlot(obj4)
FeaturePlot(obj4,features = c('CD4','CD8A','CD8B','CD44','NKG7','GZMA','CCR7'))

obj4 <- subset(obj4,seurat_clusters %in% c(0,1,3,4,5,6,7))
obj4 <- FindClusters(obj4,resolution = 0.5)
DimPlot(obj4)
obj4 <- subset(obj4,seurat_clusters %in% c(0,1,2,3,4,5,6,8,9,10,11,12,13))

obj4 <- FindClusters(obj4,resolution = 0.5)
DimPlot(obj4)
FeaturePlot(obj4,features = c('CD74','GZMB','FOXP3','LGALS1','VIM','IFITM2','GNLY','NKG7'))
FeaturePlot(obj4,features = c('IL7R','S100A4','S100A6','CTLA4','LEF1','TCF7','SELL','CCR7'))
CD4_act <- subset(obj4,seurat_clusters %in% c(2,3,4,5,7,9,11,12,13))
FeaturePlot(CD4_act,features = c('CD74','GZMK','GZMB','FOXP3','LGALS1','VIM','IFITM2','GNLY','NKG7'))
FeaturePlot(CD4_act,features = c('CD4','CD8A','CD8B','CD44','NKG7','GZMA','CCR7'))
DimPlot(CD4_act,group.by = 'sampletype')
CD4_act <- subset(CD4_act,sampletype %in% c('DKD','LD'))

save(CD4_act,file = 'E:/ZJlab/project/WXYsproject/GSE183279/process_data/ScCD4_act.Rdata')
save(obj4,file = 'E:/ZJlab/project/WXYsproject/GSE183279/process_data/ScCD4.Rdata')


load('E:/ZJlab/project/WXYsproject/GSE183279/process_data/ScCD4.Rdata')
FeaturePlot(obj4,features = c('CD4','CD8A','CD8B','CD44','NKG7','GZMA','CCR7'))


FeaturePlot(obj4,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))

FeaturePlot(obj4,features = c('IL7R','CD52','GIMAP7','SARAF','BTG1','CXCR4'))
FeaturePlot(obj4,features = c('LEF1','TCF7','SELL','CCR7','GNLY'))
FeaturePlot(obj4,features = c('IL7R','CD40LG','ANXA1','FOS','JUN'))
FeaturePlot(obj4,features = c('FOXP3','SAT1','IL2RA','CTLA4'))
DimPlot(obj4)
obj4 <- subset(obj4,sampletype %in% c('DKD','LD'))
obj4$Diabetes <- 'Healthy'
obj4$Diabetes[obj4$diabetes_history == 'Yes'] <- 'T2D'
CD4_naive <- subset(obj4,seurat_clusters %in% c(0,1,6))
CD4_mem <- subset(obj4,seurat_clusters %in% c(2,5,11))



marker2 <- FindMarkers(CD4_mem,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)
marker3 <- FindMarkers(CD4_naive,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)


