library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)

data1 <- readRDS('918c5641-444f-4b58-b6ea-131b9785d885.rds')

meta1 <- data1@meta.data
tmp1 <- data.frame(table(meta1$author_cell_type))
obj1 <- subset(data1,author_cell_type %in% c('CD4+ KLRB1 T cells','Gd T cells','NKT','Naive/CM CD4+ T cells'))
dim(obj1)


DimPlot(obj1,group.by = 'author_cell_type')

gene_name1 <- rownames(obj1)
gene.df <- bitr(gene_name1,fromType = 'ENSEMBL',
                toType = c('SYMBOL'),
                OrgDb = org.Hs.eg.db,drop = T)

gene.df <- gene.df[!duplicated(gene.df$SYMBOL),]
gene.df <- gene.df[!duplicated(gene.df$ENSEMBL),]
obj2 <- obj1[rownames(obj1) %in% gene.df$ENSEMBL,]
identical(rownames(obj2),gene.df$ENSEMBL)


T_data1 <- obj2@assays$RNA@counts
identical(rownames(T_data1),gene.df$ENSEMBL)
rownames(T_data1) <- gene.df$SYMBOL
obj3 <- CreateSeuratObject(counts = T_data1,meta.data = obj2@meta.data)
obj3 <- SCTransform(obj3)
obj3 <- RunPCA(obj3)
obj3 <- RunUMAP(obj3,dims=1:20)
obj3 <- FindNeighbors(obj3,dims=1:15)
obj3 <- FindClusters(obj3,resolution = 0.1)
DimPlot(obj3)
FeaturePlot(obj3,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))#

DimPlot(obj3,group.by = 'author_cell_type')
FeaturePlot(obj3,features = c('IL7R','CD74','GZMB','FOXP3','LGALS1','VIM','IFITM2','GNLY','NKG7','DDX21','CD69','HSPA8'))

DimPlot(obj3,group.by ='donor_id')


obj3 <- subset(obj3,donor_id %in% c('H04','H10','H11','H13','H14','H16','H18',
                                       'H21','H22','H23','H25','H37','H38'))
DimPlot(obj3,group.by ='donor_id')

obj3$Diabetes <- 'Healthy'
obj3$Diabetes[obj3$donor_id %in% c('H04','H13','H21','H37')] <- 'T2D'

cd4 <- subset(obj3,author_cell_type %in% c('Naive/CM CD4+ T cells'))
FeaturePlot(cd4,features = c('LEF1','TCF7','SELL','CCR7'))
DimPlot(cd4)

cd4 <- RunPCA(cd4)
cd4 <- RunUMAP(cd4,dims=1:20)
cd4 <- FindNeighbors(cd4,dims=1:20)
cd4 <- FindClusters(cd4,resolution = 0.2)
cd4_act <- subset(obj3,author_cell_type %in% c('CD4+ KLRB1 T cells'))
cd4_naive <- subset(cd4,seurat_clusters %in% c(2,4))
cd4_mem <- subset(cd4,seurat_clusters %in% c(0,1,3))





marker2 <- FindMarkers(cd4_mem,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0)
marker3 <- FindMarkers(cd4_naive,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0)











