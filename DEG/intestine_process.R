library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
data1 <- readRDS('f6f2d5f6-8290-42b2-9f98-e19905fa0271.rds')


tmp1 <- data.frame(table(data1$cell_type))
tmp2 <- data.frame(table(data1$author_cell_type))
obj1 <- subset(data1,cell_type %in% c('activated CD4-positive, alpha-beta T cell','CD4-positive, alpha-beta T cell',
                                      'T follicular helper cell','T-helper 1 cell','T-helper 17 cell'))
dim(obj1)


DimPlot(obj1,group.by = 'author_cell_type')

DimPlot(obj1,group.by = 'cell_type')

gene_name1 <- rownames(obj1)
gene.df <- bitr(gene_name1,fromType = 'ENSEMBL',
                toType = c('SYMBOL'),
                OrgDb = org.Hs.eg.db,drop = T)

gene.df <- gene.df[!duplicated(gene.df$SYMBOL),]
gene.df <- gene.df[!duplicated(gene.df$ENSEMBL),]
obj2 <- obj1[rownames(obj1) %in% gene.df$ENSEMBL,]
identical(rownames(obj2),gene.df$ENSEMBL)


T_data1 <- obj2@assays$RNA$data
identical(rownames(T_data1),gene.df$ENSEMBL)
rownames(T_data1) <- gene.df$SYMBOL
obj3 <- CreateSeuratObject(counts = T_data1,meta.data = obj2@meta.data)
obj3 <- SCTransform(obj3)
obj3 <- RunPCA(obj3)
obj3 <- RunUMAP(obj3,dims=1:20)
obj3 <- FindNeighbors(obj3,dims=1:15)
obj3 <- FindClusters(obj3,resolution = 0.1)
DimPlot(obj3)
DimPlot(obj3,group.by = 'cell_type')
DimPlot(obj3,group.by = 'author_cell_type')
DimPlot(obj3,group.by = 'donor_id')

obj3 <- subset(obj3,donor_id %in% c('A26 (386C)','A30 (398B)','A32 (411C)',
                                          'A34 (417C)','A38 (432C)','A39 (440C)'))

FeaturePlot(obj3,features = c('CD4','CD8A','CD44','NKG7','GZMA','CCR7'))#,'','CD44','NKG7','GZMA','CCR7'


FeaturePlot(obj3,features = c('IL7R','CD52','GIMAP7','SARAF','BTG1','CXCR4'))
FeaturePlot(obj3,features = c('LEF1','TCF7','SELL','CCR7'))
FeaturePlot(obj3,features = c('IL7R','CD40LG','ANXA1','FOS','JUN'))
FeaturePlot(obj3,features = c('FOXP3','SAT1','IL2RA','CTLA4'))
FeaturePlot(obj3,features = c('GNLY','IL2','NKG7','VIM','S100A4','IFITM2','CD74','HSP90AA1'))
FeaturePlot(obj3,features = c('GNLY','CD69','EOMES','GZMA','CXCL13','CTLA4'))


obj3$Diabetes <- 'Healthy'
obj3$Diabetes[obj3$donor_id %in% c('A26 (386C)')] <- 'T2D'
CD4_act <- subset(obj3,author_cell_type %in% c('Activated CD4 T'))
CD4_naive <- subset(obj3,author_cell_type %in% c('SELL+ CD4 T','Th1','Th17'))
CD4_mem <- obj3


marker2 <- FindMarkers(CD4_mem,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)
marker3 <- FindMarkers(CD4_naive,ident.1 = 'T2D',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0,recorrect_umi = FALSE)







