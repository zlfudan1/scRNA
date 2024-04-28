library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(monocle3)
library(clustree)
library(ggplot2)
library(harmony)
library(devtools)
library(ggrepel)
library(ggrepel)
#------------------------------------------
setwd("./Bcell")
load(file="./scRNA-umap-celltype.Rdata")

sce=subset(scRNA,ident=c("B_cells"))
table(sce@meta.data$orig.ident)
sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)

sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                          "percent.mt", "nCount_RNA"), verbose = T)
sce <- RunPCA(sce, features = VariableFeatures(object = sce),npcs = 50)
# 
library(harmony)
sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:25, reduction = "harmony")
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 
set.resolutions = c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1)
sce  <- FindClusters(object = sce , resolution = set.resolutions, verbose = FALSE) 
# 
for (res in c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1)) {
  sce=FindClusters(sce, #graph.name = "CCA_snn", 
                   resolution = res, algorithm = 1)
}
clustree(sce)
colnames(sce@meta.data)

sce  <- FindClusters(object = sce , resolution = 0.05, verbose = FALSE) 
Bcell <- RunUMAP(sce, reduction = "harmony", dims = 1:25)
Bcell <- RunTSNE(sce, reduction = "harmony", dims =1:25)#

scRNA_harmony <- Bcell
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "4harmonyBcell-0.05_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "top50_diff_genes_wilcox.csv", row.names = F)
#
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 25   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'4harmonyBcell-0.05_TopMarkersTOP25.csv', row.names=F)
#-------------------------------------------------------------------------------
Idents(Bcell)="seurat_clusters"
Bcell2 <- subset(Bcell,ident=c('0', '1', '2','3', "4",'6','7'))
sce <- Bcell2

sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)
#
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                          "percent.mt", "nCount_RNA"), verbose = T)
sce <- RunPCA(sce, features = VariableFeatures(object = sce), npcs = 50)
# 
sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:25, reduction = "harmony")
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 
set.resolutions = c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1)
sce  <- FindClusters(object = sce , resolution = set.resolutions, verbose = FALSE) 
for (res in c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1)) {
  sce=FindClusters(sce,resolution = res, algorithm = 1)
}
clustree(sce)
colnames(sce@meta.data)

sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25)
sce  <- FindClusters(object = sce , resolution = 0.05, verbose = FALSE) 
Bcell2 <- RunUMAP(sce, reduction = "harmony", dims = 1:25)
Bcell2 <- RunTSNE(Bcell2, reduction = "harmony", dims =1:25)#

DimPlot(Bcell2, reduction = "umap",label = T) 
DimPlot(Bcell2, reduction = "tsne",label = T) 
DimPlot(Bcell2, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE,split.by = "group")
DimPlot(Bcell2, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE,split.by = "spatialdistr")
DimPlot(Bcell2, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE,split.by = "metastasis")

scRNA_harmony <- Bcell2
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "4harmonyBcell2-0.05_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "top50_diff_genes_wilcox.csv", row.names = F)
#
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 25  
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'4harmonyBcell2-0.05_TopMarkersTOP25.csv', row.names=F)
#-------------------------------------------------------------------------------
load(file="./4celltypeharmonyBcell-0.1.Rdata")
DimPlot(Bcell, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(Bcell, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(Bcell, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE,split.by = "group")
DimPlot(Bcell, reduction = "umap", group.by='celltype',label = T,raster=FALSE)
VlnPlot(Bcell, 
        features = c('SDC1', 'MS4A1','IGHG4', 'IGHG1', 'TNFRSF17', 'CD38', 'MZB1'),
        pt.size = 0,
        ncol = 2)

#B naive
VlnPlot(Bcell, 
        features = c('TCL1A', 'IGHD', 'IL4R', 'FCER2', 'CD19','CD22','CD83'),
        pt.size = 0,
        ncol = 2)
#B memory
VlnPlot(Bcell, 
        features = c('AIM2', 'TNFRSF13B', 'CD27'),
        pt.size = 0,
        ncol = 2)
#B plasma
VlnPlot(Bcell, 
        features = c('MZB1', 'IGHG3', 'JCHAIN'),
        pt.size = 0,
        ncol = 2)
#B germinal center
VlnPlot(Bcell, 
        features = c('S1PI2', 'LRMP', 'SUGCT', 'MME', 'MKI67', 'AICDA'),
        pt.size = 0,
        ncol = 2)

#
df_Bcell <- Bcell@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = Bcell@meta.data$celltype)%>%
  cbind(seurat_clusters = Bcell@meta.data$seurat_clusters)
colnames(df_Bcell)
#
cluster_order <- c('0','1',"2","3")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D'),
                        c('0','1',"2","3"))
#
cell <- df_Bcell %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:3)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_Bcell)
ggplot(df_Bcell,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
  geom_point(size = 1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.6,'cm'),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5, arrow = arrow(length = unit(0.2,"cm"), type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), type = "closed")) +
  annotate("text", x = min(df_Bcell$UMAP_1) +1, y = min(df_Bcell$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_Bcell$UMAP_1) -0.5, y = min(df_Bcell$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 4,angle=90)+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("Bcell_UMAP.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)

#
B <- A
B$x <- 1
B$lab <- c(3:0)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>% 
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=seurat_clusters)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 4)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)
ggsave("Bcell_umap_legend.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
###-----------------------------------------------------------------------------
#
Bcell_CAYA = Bcell[,Bcell@meta.data$group%in%c('CAYA')]
df_Bcell_CAYA <- Bcell_CAYA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(Bcell_CAYA = Bcell_CAYA@meta.data$celltype)%>%
  cbind(seurat_clusters = Bcell_CAYA@meta.data$seurat_clusters)
colnames(df_Bcell_CAYA)
#
cluster_order <- c('0','1',"2","3")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D'),
                        c('0','1',"2","3"))
#
cell_CAYA <- df_Bcell_CAYA %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell_CAYA) <-cell_CAYA$seurat_clusters
A <- cell_CAYA[cluster_order,]
a <- c(0:3)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_Bcell_CAYA)
ggplot(df_Bcell_CAYA,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
  geom_point(size = 1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.6,'cm'),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5, arrow = arrow(length = unit(0.2,"cm"), type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), type = "closed")) +
  annotate("text", x = min(df_Bcell_CAYA$UMAP_1) +1, y = min(df_Bcell_CAYA$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_Bcell_CAYA$UMAP_1) -0.5, y = min(df_Bcell_CAYA$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 4,angle=90)+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("Bcell_CAYA_UMAP.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
###-----------------------------------------------------------------------------
#
Bcell_ADULT = Bcell[,Bcell@meta.data$group%in%c('ADULT')]
df_Bcell_ADULT <- Bcell_ADULT@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(Bcell_ADULT = Bcell_ADULT@meta.data$celltype)%>%
  cbind(seurat_clusters = Bcell_ADULT@meta.data$seurat_clusters)
colnames(df_Bcell_ADULT)
#
cluster_order <- c('0','1',"2","3")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D'),
                        c('0','1',"2","3"))
#
cell_ADULT <- df_Bcell_ADULT %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell_ADULT) <-cell_ADULT$seurat_clusters
A <- cell_ADULT[cluster_order,]
a <- c(0:3)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_Bcell_ADULT)
ggplot(df_Bcell_ADULT,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
  geom_point(size = 1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.6,'cm'),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5, arrow = arrow(length = unit(0.2,"cm"), type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), type = "closed")) +
  annotate("text", x = min(df_Bcell_ADULT$UMAP_1) +1, y = min(df_Bcell_ADULT$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_Bcell_ADULT$UMAP_1) -0.5, y = min(df_Bcell_ADULT$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 4,angle=90)+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("Bcell_ADULT_UMAP.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
###-----------------------------------------------------------------------------
#
table(Bcell$HT)
prop.table(table(Bcell$celltype))
Cellratio_HT <- prop.table(table(Bcell$celltype, Bcell$HT), margin = 2)
Cellratio_HT <- as.data.frame(Cellratio_HT)
write.csv(Cellratio_HT, "Cellratio_HT.csv", row.names = T)
#
table(Bcell$group)
prop.table(table(Bcell$celltype))
Cellratio_group <- prop.table(table(Bcell$celltype, Bcell$group), margin = 2)
Cellratio_group <- as.data.frame(Cellratio_group)
write.csv(Cellratio_group, "Cellratio_group.csv", row.names = T)
#
table(Bcell$metastasis)
prop.table(table(Bcell$celltype))
Cellratio_metastasis <- prop.table(table(Bcell$celltype, Bcell$metastasis), margin = 2)
Cellratio_metastasis <- as.data.frame(Cellratio_metastasis)
write.csv(Cellratio_metastasis, "Cellratio_metastasis.csv", row.names = T)
#
table(Bcell$spatialdistr)
prop.table(table(Bcell$celltype))
Cellratio_spatialdistr <- prop.table(table(Bcell$celltype, Bcell$spatialdistr), margin = 2)
Cellratio_spatialdistr <- as.data.frame(Cellratio_spatialdistr)
write.csv(Cellratio_spatialdistr, "Cellratio_spatialdistr.csv", row.names = T)

#
Idents(Bcell)="metastasis"
table(Bcell@meta.data$metastasis)
Bcell=RenameIdents(Bcell,c('N0'='N0','N1a'='N1','N1b'='N1'))
Bcell@meta.data$metastasis2=Bcell@active.ident
Idents(Bcell)="metastasis2"
table(Bcell@meta.data$metastasis2)
#
table(Bcell$metastasis2)
prop.table(table(Bcell$celltype))
Cellratio_metastasis2 <- prop.table(table(Bcell$celltype, Bcell$metastasis2), margin = 2)
Cellratio_metastasis2 <- as.data.frame(Cellratio_metastasis2)
write.csv(Cellratio_metastasis2, "Cellratio_metastasis2.csv", row.names = T)

data <- read.csv("cellratio_Bcell.csv", header = T)
data$var <- factor(data$var, levels=c('B_naive','B_GC','B_memory','B_plasma'))
ggplot(data,aes(x = flag, y= Frequence, fill = var))+ 
  geom_col()+ labs(x= NULL, y="Frequence", order = T)+
  geom_bar(stat = 'identity')+ 
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+ 
  theme(axis.text = element_text(colour = 'black'))+
  scale_x_discrete(limits = c('CAYA', 'ADULT', 'Tumor', 'Normal', 'HT', 'NHT', 'N0', 'N1'))+
  theme_classic()+
  theme(axis.title.x = element_text(family = 'Times New Roman', size = 12),
        axis.title.y = element_text(family = 'Times New Roman', size = 12),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size = 8, colour = 'black',
                                   family = 'Times New Roman',vjust = 1, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 11, colour = 'black',
                                   family = 'Times New Roman',vjust = 0.5, hjust = 1))+
  scale_x_discrete(limits = c('CAYA', 'ADULT', 'Tumor', 'Normal', 'HT', 'NHT', 'N0', 'N1'))+
  scale_fill_manual(values = c('#349573', '#C7591A', '#6F6BA5','#CE2B7D'))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 8.5), linetype= "dashed", linewidth = 0.5)
ggsave("cellratio_Bcell.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)