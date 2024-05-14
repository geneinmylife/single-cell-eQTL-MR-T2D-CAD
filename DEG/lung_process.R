library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(SeuratDisk)
library(SingleCellExperiment)
library(anndata)
library(reticulate)
library(rhdf5)

use_condaenv(condaenv = "", required = TRUE)

ad <- import("anndata")
h5ls("lung_5loc_sc_sn_TNK_cellxgene_updnames_allgenes_05122022.h5ad")

final_ad <- ad$read_h5ad("lung_5loc_sc_sn_TNK_cellxgene_updnames_allgenes_05122022.h5ad")

colnames(final_ad$X) <- final_ad$var$features
counts <- final_ad$raw$X
colnames(counts) <- final_ad$raw$var_names
rownames(counts) <- final_ad$raw$obs_names
counts <- t(counts)
meta1 <- final_ad$obs

tmp1 <- final_ad$raw$obs_names



seurat_obj <- CreateSeuratObject(counts = counts, assay = "RNA",
                                         meta.data = meta1)
table(seurat_obj$Celltypes)
meta1 <- seurat_obj@meta.data

CD4_obj <- subset(seurat_obj,Celltypes %in% c('CD4_EM/Effector'))
CD4_obj <- SCTransform(CD4_obj)
CD4_obj <- RunPCA(CD4_obj)
CD4_obj <- RunUMAP(CD4_obj,dims = 1:15)

DimPlot(CD4_obj,group.by = 'Celltypes')
DimPlot(CD4_obj,group.by = 'Donor')

CD4_obj$Diabetes <- 'Healthy'
CD4_obj$Diabetes[CD4_obj$Donor %in% c('A26')] <- 'Diabetes'
table(CD4_obj$Diabetes)

cd4_act <- CD4_obj

marker1 <- FindMarkers(CD4_obj,ident.1 = 'Diabetes',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0)



CD4_Mem <- subset(seurat_obj,Celltypes %in% c('CD4_TRM'))
CD4_Mem <- SCTransform(CD4_Mem)
CD4_Mem <- RunPCA(CD4_Mem)
CD4_Mem<- RunUMAP(CD4_Mem,dims = 1:15)



CD4_Mem$Diabetes <- 'Healthy'
CD4_Mem$Diabetes[CD4_Mem$Donor %in% c('A26')] <- 'Diabetes'
table(CD4_Mem$Diabetes)



marker2 <- FindMarkers(CD4_Mem,ident.1 = 'Diabetes',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0)


CD4_Naive <- subset(seurat_obj,Celltypes %in% c('CD4_naive/CM'))
CD4_Naive <- SCTransform(CD4_Naive)
CD4_Naive <- RunPCA(CD4_Naive)
CD4_Naive <- RunUMAP(CD4_Naive,dims = 1:15)



CD4_Naive$Diabetes <- 'Healthy'
CD4_Naive$Diabetes[CD4_Naive$Donor %in% c('A26')] <- 'Diabetes'
table(CD4_Naive$Diabetes)


marker3 <- FindMarkers(CD4_Naive,ident.1 = 'Diabetes',group.by = 'Diabetes',min.pct = 0.1,logfc.threshold = 0)













