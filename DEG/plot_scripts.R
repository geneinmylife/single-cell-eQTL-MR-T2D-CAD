library(readxl)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(ggheatmap)
library(dplyr)

mytheme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent",colour=NA),
    plot.background = element_rect(fill = "transparent",colour=NA),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 9),
    axis.line = element_line(colour = "black"),
    legend.background = element_rect(fill = "transparent",colour=NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

gene_list <- list()
num1 <- 1
sample_list <- c('lung','liver','intestine','Kidney','Wound-edge','PBMC1','PBMC2')
for(i in c('lung','liver','intestine','GSE183279','GSE176415',
           'GSE244515','GSE165816')){
  tmp1 <- read.csv(paste0(i,'/process_data/DEG_2.csv'))
  tmp1$sample <- sample_list[num1]
  gene_list[[num1]] <- tmp1
  num1 <- num1+1
}
gene_list_2 <- do.call(rbind,gene_list)
table(gene_list_2$sample)


T2D_gene_list <- data.frame(read_xlsx('Gene_list_name_last.xlsx',sheet=1))

T2D_gene_list_2 <- data.frame(read_xlsx('results/combined T2D new results.xlsx'))
T2D_gene_list_3 <- data.frame(read_xlsx('results/T2D revised wald ratio.xlsx',sheet = 4))
T2D_gene_list_2$ENSGID <- str_split_fixed(T2D_gene_list_2$phe,'_',2)[,1]

T2D_gene_list_2$exposure <- str_split_fixed(T2D_gene_list_2$phe,'_',2)[,2]
T2D_gene_list_2$Cell_type <- 'Memory'
T2D_gene_list_2$Cell_type[grepl('TN',T2D_gene_list_2$exposure,ignore.case = T) | grepl('Naive',T2D_gene_list_2$exposure,ignore.case = T)] <- 'Naive'
T2D_gene_list_2$Cell_type[grepl('nTreg',T2D_gene_list_2$exposure,ignore.case = T)] <- 'Treg'
T2D_gene_list_2$Cell_type[grepl('T_ER-stress',T2D_gene_list_2$exposure,ignore.case = T)] <- 'T_ER-stress'
T2D_gene_list_2 <- na.omit(T2D_gene_list_2)
T2D_gene_list_2 <- subset(T2D_gene_list_2,Cell_type %in% c('Naive','Memory'))
length(unique(c(T2D_gene_list_2$ENSGID,T2D_gene_list_3$Gene)))

T2D_gene_list_2_sum <- T2D_gene_list_2 %>% group_by(ENSGID) %>% summarise(mean(beta))
colnames(T2D_gene_list_2_sum) <- c('ENSGID','Mean_b')


T2D_gene_list_3$Cell_type <- 'Memory'
T2D_gene_list_3$Cell_type[grepl('TN',T2D_gene_list_3$Cell_type_time_point,ignore.case = T) | grepl('Naive',T2D_gene_list_3$Cell_type_time_point,ignore.case = T)] <- 'Naive'
T2D_gene_list_3$Cell_type[grepl('nTreg',T2D_gene_list_3$Cell_type_time_point,ignore.case = T)] <- 'Treg'


T2D_gene_list_3_sum <- T2D_gene_list_3 %>% group_by(Gene) %>% summarise(mean(b))
colnames(T2D_gene_list_3_sum) <- c('ENSGID','Mean_b')
T2D_gene_list_sum <- data.frame(rbind(T2D_gene_list_2_sum,T2D_gene_list_3_sum))



mem_DEG_list <- c(T2D_gene_list_2$Gene[T2D_gene_list_2$Cell_type=='Memory'],T2D_gene_list_3$Gene_name[T2D_gene_list_3$Cell_type=='Memory'])
Naive_DEG_list <- c(T2D_gene_list_2$Gene[T2D_gene_list_2$Cell_type=='Naive'],T2D_gene_list_3$Gene_name[T2D_gene_list_3$Cell_type=='Naive'])


check_direction <- function(x,gene_data=T2D_gene_list_sum){
  gene_name <- x[1]
  beta <- T2D_gene_list_sum[gene_name,'Mean_b']
  logFC <- as.numeric(x[3])
  if(beta*logFC>=0){
    return(T)
  } else{
    return(F)
  }
}



act_gene_merge <- data.frame(gene_list_2[(gene_list_2$X %in% T2D_gene_list_sum$Gene_name),])
act_gene_merge$check <- apply(act_gene_merge,1,check_direction)
act_gene_merge <- subset(act_gene_merge,check)
act_gene_merge <- subset(act_gene_merge,p_val<0.05)


result_2 <- act_gene_merge %>% group_by(X) %>% summarise(mean(avg_log2FC),length(X))
colnames(result_2) <- c('Gene','avg_log2FC','num')
result_2$Cell_type <- 'CD4_mem'


mem_DEG <- read.csv("E:/ZJlab/project/WXYsproject/results/merge_genes_with_DEG_of_mem_CD4.csv")
mem_DEG$verifiedinMR <- mem_DEG$Gene %in% mem_DEG_list


Naive_DEG <- read.csv("E:/ZJlab/project/WXYsproject/results/merge_genes_with_DEG_of_naive_CD4.csv")
Naive_DEG$verifiedinMR <- Naive_DEG$Gene %in% Naive_DEG_list




result_1_1 <- subset(result_1,num>=3)
result_2_1 <- subset(result_2,num>=3)


result_1_1 <- result_1_1[order(result_1_1$num,decreasing=T),]
result_2_1 <- result_2_1[order(result_2_1$num,decreasing=T),]


result_merge <- rbind(result_1,result_2)
order <- rev(unique(c(result_2_1$Gene,result_1_1$Gene)))

result_merge$Cell_type <- factor(result_merge$Cell_type,levels = c('CD4_naive','CD4_mem'))
result_merge_1 <- subset(result_merge,Gene %in% c(result_1_1$Gene,result_2_1$Gene))




colnames(result_merge_1) <- c("Gene","avg_log2FC","Confidence","Cell_type") 


order1 <- subset(result_merge_1,Cell_type=='CD4_naive'& avg_log2FC>0)
order2 <- subset(result_merge_1,Cell_type=='CD4_naive'& avg_log2FC<0)
order1 <- order1[order(order1$Confidence,decreasing=T),]
order2 <- order2[order(order2$Confidence,decreasing=F),]
order3 <- rbind(order1,order2)

result_merge_1$Gene <- factor(result_merge_1$Gene,levels = rev(order3$Gene))

p1 <- ggplot()+
  geom_point(data=result_merge_1,aes(x=Cell_type,y=Gene,size=Confidence,color=avg_log2FC))+
  scale_color_gradientn(colors=rev(colorRampPalette(c('#B54152','#D17F6F','#EEAF93','#F2D6C6','white','#C7D9E2','#8BBCD7','#119CD0','#3479AD'))(100)),
                        limits=c(-0.5,0.5),oob = scales::squish)+
  scale_size_area(breaks=c(1,3,5,7))+
  mytheme+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 8))+
  guides(color = guide_colorbar(barwidth = 0.5,barheight=5,title = 'avg_log2FC'))
p1





act_gene_merge <- data.frame(gene_list_2[(gene_list_2$X %in% T2D_gene_list_sum$Gene_name),])


act_gene_merge$check <- apply(act_gene_merge,1,check_direction)
act_gene_merge_2 <- subset(act_gene_merge,check)
act_gene_merge_2 <- subset(act_gene_merge_2,p_val<0.05)



p1 <- ggplot()+
  geom_bar(data = act_gene_merge,aes(x=sample),fill='#FDCFCD')+
  geom_bar(data = act_gene_merge_2,aes(x=sample),fill='#921A7D')+
  mytheme+
  # scale_y_reverse(expand=c(0,0))+
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank()
  )+
  scale_x_discrete()+
  scale_y_continuous(expand=c(0,0))+
  coord_flip()+
  theme(axis.text.x = element_text(size = 8))+
  ylab('')+xlab('')
p1

p2 <- ggplot()+
  geom_bar(data = act_gene_merge,aes(x=sample,fill=check))+
  scale_fill_manual(values = c('#FDCFCD','#921A7D'))+
  # geom_bar(data = act_gene_merge_2,aes(x=sample),fill='#921A7D')+
  mytheme
  theme(axis.text.x = element_text(size = 5,angle=45,hjust = 1))+
  ylab('Detected Gene counts')+xlab('Datasets')
p2






DimPlot(cd4_mem,group.by = 'Diabetes')




mytheme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent",colour=NA),
    plot.background = element_rect(fill = "transparent",colour=NA),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 9),
    axis.line = element_line(colour = "black"),
    legend.background = element_rect(fill = "transparent",colour=NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


plotDF1 <- data.frame(cd4_mem@reductions$umap@cell.embeddings,cd4_mem@assays$SCT@counts['CD52',],cd4_mem$Diabetes)
colnames(plotDF1) <- c('UMAP_1','UMAP_2','Gene','Diabetes')
p3 <- ggplot()+
  geom_point(data = plotDF1,aes(x=UMAP_1,y=UMAP_2,color=Gene),size=0.6)+
  scale_color_gradient(low = '#FDCFCD',high = '#921A7D',limits=c(0,10),oob = scales::squish)+
  mytheme
p3



p4 <- ggplot()+
  geom_point(data = plotDF1,aes(x=UMAP_1,y=UMAP_2,color=Diabetes),size=0.6)+
  scale_color_manual(values = c('#AFAED8','#ADD0B2'))+
  mytheme+guides(color = guide_legend(override.aes = list(size = 5)))
p4







