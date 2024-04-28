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
rm(list = ls())
load(file="./scRNA-umap-celltype.Rdata")
DimPlot(scRNA, reduction = "tsne", label=T) 
DimPlot(scRNA, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
FeaturePlot(scRNA,"CD8B")
FeaturePlot(scRNA,"CD8A")
FeaturePlot(scRNA,"CD4")
load(file="./NEWCD4Tcell_scRNA.anchors.Rdata")
setwd("./subcelltype/CD4T")
###################################################################
table(scRNA@meta.data$seurat_clusters)
Idents(scRNA)="seurat_clusters"
CD4sc.Tcell=subset(scRNA,ident=c('4', '6', '8', '15', '30'))
table(CD4sc.Tcell@meta.data$orig.ident)

sce.list2 <- SplitObject(CD4sc.Tcell, split.by = "orig.ident")
for (i in 1:length(sce.list2)) {
  sce.list2[[i]] <- NormalizeData(sce.list2[[i]],verbose=FALSE)
  sce.list2[[i]] <- FindVariableFeatures(sce.list2[[i]], selection.method = "vst",nfeatures = 2000,verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = sce.list2)
scRNA.anchors <- FindIntegrationAnchors(object.list = sce.list2, anchor.features = features)
save(scRNA.anchors,file="NEWCD4Tcell_scRNA.anchors.Rdata")
load(file="./NEWCD4Tcell_scRNA.anchors.Rdata")
CD4T <- IntegrateData(anchorset = scRNA.anchors)

DefaultAssay(CD4T) <- "RNA"
#
CD4T <- NormalizeData(CD4T)
CD4T <- FindVariableFeatures(CD4T, selection.method = "vst", nfeatures = 2000)
#
DefaultAssay(CD4T) <- "integrated"
CD4T<- ScaleData(CD4T, features = VariableFeatures(CD4T))
#
CD4T <- RunPCA(CD4T, features = VariableFeatures(CD4T))
ElbowPlot(CD4T)
CD4T <- FindNeighbors(CD4T, dims = 1:20)
save(scRNA.anchors,file="NEWCD4Tcell_FindNeighbors.Rdata")
###
res.used <- seq(0.1,1,by=0.1)
res.used
#
for(i in res.used){
  Tcell <- FindClusters(object = CD4T, verbose = T, resolution = res.used)
}
#
library(clustree)
clus.tree.out <- clustree(CD4T) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")

clus.tree.out
##
CD4T<- FindClusters(CD4T, resolution = 0.5) 
pc.num=1:20 
#tSNE
CD4T = RunTSNE(CD4T, dims = pc.num)
CD4T <- RunUMAP(CD4T, dims = pc.num)
#
DimPlot(CD4T, reduction = "tsne", label=T) 
DimPlot(CD4T, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
###---------------------------------------
DefaultAssay(CD4T) <- "RNA"

diff.wilcox = FindAllMarkers(CD4T,logfc.threshold = 0.25,only.pos = T,test.use = "MAST")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "CD4T_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "CD4T_top50_diff_genes_wilcox.csv", row.names = F)
###
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50 
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP50<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP50,'NEWCD4T_R0.4_TopMarkersTOP50.csv', row.names=F)

###
Tcell <- FindNeighbors(Tcell, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.4)
Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:20)
Tcell <- RunTSNE(Tcell, reduction = "harmony", dims = 1:20)
DimPlot(Tcell, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(Tcell, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
FeaturePlot(Tcell,"CD8B")

#combinate
scRNA_harmony <- Tcell
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "harmonyTcell-0.1_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "top50_diff_genes_wilcox.csv", row.names = F)
### 
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'harmonyTcell-0.4_TopMarkersTOP50.csv', row.names=F)

scRNA_harmony <- Tcell
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "harmonyTcell-0.1_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "top50_diff_genes_wilcox.csv", row.names = F)
###
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50 
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'harmonyTcell-0.3_TopMarkersTOP50.csv', row.names=F)

Tcell <-scRNA_harmony 
save(Tcell,all.markers,file="harmonyTcell-0.1.Rdata") 
load(file="./harmonyTcell-0.1.Rdata")

Idents(CD4T)="seurat_clusters"
CD4T2 <- subset(CD4T,ident=c('0','1','2',"3", '4', '7', '8', '9', '11'))
CD4T2 <- NormalizeData(CD4T2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
DefaultAssay(CD4T2) <- "RNA"
#
CD4T2 <- NormalizeData(CD4T2)
CD4T2 <- FindVariableFeatures(CD4T2, selection.method = "vst", nfeatures = 2000)
#
DefaultAssay(CD4T2) <- "integrated"
CD4T2<- ScaleData(CD4T2, features = VariableFeatures(CD4T2))
#
CD4T2 <- RunPCA(CD4T2, features = VariableFeatures(CD4T2))
ElbowPlot(CD4T2)
CD4T2 <- FindNeighbors(CD4T2, dims = 1:20)
CD4T2<- FindClusters(CD4T2, resolution = 0.4) 
pc.num=1:20
#tSNE
CD4T2 = RunTSNE(CD4T2, dims = pc.num)
CD4T2 <- RunUMAP(CD4T2, dims = pc.num)
#group_by
DimPlot(CD4T2, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(CD4T2, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)

DefaultAssay(CD4T2) <- "RNA"
#####
diff.wilcox = FindAllMarkers(CD4T2,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)

write.csv(all.markers, "CD8T_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "CD8T_R0.2_top50_diff_genes_wilcox.csv", row.names = F)
###
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP50<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP50,'NEW_CD4T2_R0.4_deleted_TopMarkersTOP50.csv', row.names=F)
save(CD4T2,file = './final_CD4T2.Rdata')
###-----------------------------------------------------------------------------
Idents(CD4T2)="seurat_clusters"
CD4T3 <- subset(CD4T2,ident=c('0','1','2',"3", '4', '5', '8'))
CD4T3 <- NormalizeData(CD4T3) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
DefaultAssay(CD4T3) <- "RNA"
#
CD4T3 <- NormalizeData(CD4T3)
CD4T3 <- FindVariableFeatures(CD4T3, selection.method = "vst", nfeatures = 2000)
#
DefaultAssay(CD4T3) <- "integrated"
CD4T3<- ScaleData(CD4T3, features = VariableFeatures(CD4T3))
#
CD4T3 <- RunPCA(CD4T3, features = VariableFeatures(CD4T3))
ElbowPlot(CD4T3)
CD4T3 <- FindNeighbors(CD4T3, dims = 1:20)
CD4T3<- FindClusters(CD4T3, resolution = 0.4) 
pc.num=1:20 
#tSNE
CD4T3 = RunTSNE(CD4T3, dims = pc.num)
DimPlot(CD4T3, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)

DefaultAssay(CD4T3) <- "RNA"
#####
diff.wilcox = FindAllMarkers(CD4T3,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)

write.csv(all.markers, "CD8T_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "CD8T_R0.2_top50_diff_genes_wilcox.csv", row.names = F)
###
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50  
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP50<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP50,'NEW_CD4T3_R0.5_deleted_2_TopMarkersTOP50.csv', row.names=F)
save(CD4T3,file = './final_CD4T3.Rdata')

M <- CD4T3
M@meta.data <- M@meta.data[,-(8:15)]
M@meta.data <- M@meta.data[,-(9:11)]
CD4T3@meta.data <- M@meta.data
FeaturePlot(CD4T3,features = "CD4")

Idents(CD4T)="celltype"
CD4T3=RenameIdents(CD4T3,c('CD4T2_Tfh_CXCL13'='CD4T_Tfh_CXCL13'))
CD4T3@meta.data$celltype=CD4T3@active.ident
Idents(CD4T3)="celltype"
table(CD4T3@meta.data$celltype)
DimPlot(CD4T3, reduction = "tsne", group.by='celltype',label = T,raster=FALSE,split.by = "group")
CD4T <- CD4T3
save(CD4T,all.markers,file="final_CD4T.Rdata") 
##------------------------------------------------------------------------------
Idents(CD4T)="celltype"
CD4T=RenameIdents(CD4T,c('CD4T2_Tfh_CXCL13'='CD4T_Tfh_CXCL13'))
CD4T@meta.data$celltype=CD4T@active.ident
Idents(CD4T)="celltype"
table(CD4T@meta.data$celltype)
DimPlot(CD4T, reduction = "tsne", group.by='celltype',label = T,raster=FALSE,split.by = "group")
CD4T <- CD4T3

##
setwd("./CD4T")
load(file="./final_CD4T.Rdata")
DimPlot(CD4T, reduction = "tsne", group.by='celltype',label = T,raster=FALSE)
df_CD4T <- CD4T@reductions$tsne@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = CD4T@meta.data$celltype)%>%
  cbind(seurat_clusters = CD4T@meta.data$seurat_clusters)
colnames(df_CD4T)
###
cluster_order <- c('0','1',"2","3", "4", "5", "6")
###
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                          '#D8A326', '#9E732B'),
                        c('0','1',"2","3", "4","5","6"))
###
cell <- df_CD4T %>%group_by(seurat_clusters) %>%
  summarise(tsne_1 = median(tSNE_1),
            tsne_2 = median(tSNE_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:6)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","tSNE_1", "tSNE_2" ,"ID")
colnames(df_CD4T)
ggplot(df_CD4T,aes(x= tSNE_1 , y = tSNE_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
  geom_segment(aes(x = min(tSNE_1) , y = min(tSNE_2),xend = min(tSNE_1)+7, yend = min(tSNE_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(tSNE_1), y = min(tSNE_2),xend = min(tSNE_1),yend = min(tSNE_2)+8),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), 
                                                        type = "closed")) +
  annotate("text", x = min(df_CD4T$tSNE_1) +3, y = min(df_CD4T$tSNE_2) -0.8, label = "tSNE_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_CD4T$tSNE_1) -0.8, y = min(df_CD4T$tSNE_2) + 4, label = "tSNE_2",
           color="black",size = 4,angle=90) +
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
ggsave("CD4T_tsne.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)

#
B <- A
B$x <- 1
B$lab <- c(6:0)
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
ggsave("CD4T_umap_legend.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
###-----------------------------------------------------------------------------
#
table(CD4T$HT)
prop.table(table(CD4T$celltype))
Cellratio_HT <- prop.table(table(CD4T$celltype, CD4T$HT), margin = 2)
Cellratio_HT <- as.data.frame(Cellratio_HT)
write.csv(Cellratio_HT, "Cellratio_HT.csv", row.names = T)
#
table(CD4T$group)
prop.table(table(CD4T$celltype))
Cellratio_group <- prop.table(table(CD4T$celltype, CD4T$group), margin = 2)
Cellratio_group <- as.data.frame(Cellratio_group)
write.csv(Cellratio_group, "Cellratio_group.csv", row.names = T)
#
table(CD4T$metastasis)
prop.table(table(CD4T$celltype))
Cellratio_metastasis <- prop.table(table(CD4T$celltype, CD4T$metastasis), margin = 2)
Cellratio_metastasis <- as.data.frame(Cellratio_metastasis)
write.csv(Cellratio_metastasis, "Cellratio_metastasis.csv", row.names = T)
#
table(CD4T$spatialdistr)
prop.table(table(CD4T$celltype))
Cellratio_spatialdistr <- prop.table(table(CD4T$celltype, CD4T$spatialdistr), margin = 2)
Cellratio_spatialdistr <- as.data.frame(Cellratio_spatialdistr)
write.csv(Cellratio_spatialdistr, "Cellratio_spatialdistr.csv", row.names = T)

setwd("./CD4T")
data <- read.csv("cellratio_CD4T.csv", header = T)
View(data)
data$var <- factor(data$var, levels=c('CD4T_Tn/Tm_ANXA1','CD4T_Tfh_CXCL13','CD4T_Treg_RTKN2','CD4T_Th17_CCL5','CD4T_Treg_TNFRSF4','CD4T_Tfh_TOX2','CD4T_ISG+T'))
ggplot(data,aes(x = flag, y= Frequence, fill = var))+ 
  geom_col()+ labs(x= NULL, y="Frequence", order = T)+
  geom_bar(stat = 'identity')+ 
  labs(x=NULL)+ 
  theme_bw(base_size = 18)+ 
  theme(axis.text = element_text(colour = 'black'))+
  scale_x_discrete(limits = c('CAYA', 'ADULT', 'Tumor', 'Normal', 'HT', 'NHT', 'N0', 'N1a', 'N1b'))+
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
  scale_x_discrete(limits = c('CAYA', 'ADULT', 'Tumor', 'Normal', 'HT', 'NHT', 'N0', 'N1a', 'N1b'))+
  scale_fill_manual(values = c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                               '#D8A326', '#9E732B'))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 9.5), linetype= "dashed", linewidth = 0.5)
ggsave("cellratio.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
#
Idents(CD4T)="metastasis"
table(CD4T@meta.data$metastasis)
CD4T=RenameIdents(CD4T,c('N0'='N0','N1a'='N1','N1b'='N1'))
CD4T@meta.data$metastasis2=CD4T@active.ident
Idents(CD4T)="metastasis2"
table(CD4T@meta.data$metastasis2)
#
table(CD4T$metastasis2)
prop.table(table(CD4T$celltype))
Cellratio_metastasis2 <- prop.table(table(CD4T$celltype, CD4T$metastasis2), margin = 2)
Cellratio_metastasis2 <- as.data.frame(Cellratio_metastasis2)
write.csv(Cellratio_metastasis2, "Cellratio_metastasis2.csv", row.names = T)
data <- read.csv("cellratio_CD4T_2.csv", header = T)
data$var <- factor(data$var, levels=c('CD4T_Tn/Tm_ANXA1','CD4T_Tfh_CXCL13','CD4T_Treg_RTKN2','CD4T_Th17_CCL5','CD4T_Treg_TNFRSF4','CD4T_Tfh_TOX2','CD4T_ISG+T'))
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
  scale_fill_manual(values = c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                               '#D8A326', '#9E732B','#656565','#A4C7DA', '#3172A5', 
                               '#ADD28A','#40973E','#EA9393'))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 8.5), linetype= "dashed", linewidth = 0.5)
ggsave("cellratio_2.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
###-----------------------------------------------------------------------------
##monocle3
setwd("./CD4T")
load(file="./final_CD4T.Rdata")
pc.num <- 1:20
CD4T <- RunUMAP(CD4T, dims = pc.num)
data <- GetAssayData(CD4T, assay = 'RNA', slot = 'counts')
cell_metadata <- CD4T@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#
cds <- preprocess_cds(cds, num_dim = 50)
#
cds <- reduce_dimension(cds,reduction_method='UMAP',preprocess_method = 'PCA')
cds <- reduce_dimension(cds,reduction_method='tSNE',preprocess_method = 'PCA')

plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") + ggtitle('cds.tsne')
##
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(CD4T, reduction = "UMAP")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds)
plot_cells(cds, reduction_method="UMAP", label_groups_by_cluster = F, label_leaves = F, 
           label_branch_points = F)
plot_cells(cds, reduction_method="tSNE", label_groups_by_cluster = F, label_leaves = F, 
           label_branch_points = F)
plot_cells(cds,reduction_method="UMAP",
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=T,
           label_branch_points=T,
           group_label_size=3,cell_size=1.5)
cds <- order_cells(cds,reduction_method = 'UMAP')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
top10 <- subset_pr_test_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

#
plot_genes_in_pseudotime(cds[top10,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)
###-----------------------------------------------------------------------------
#
cds <- reduce_dimension(cds,reduction_method='UMAP',
                        preprocess_method = 'PCA')
cds <- cluster_cells(cds)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(CD4T, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
#
plot_cells(cds,reduction_method="UMAP",
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=T,
           label_branch_points=T,
           group_label_size=3,cell_size=1.5)
#
mycds <- cds
mycds <- learn_graph(mycds,
                     verbose=T,
                     learn_graph_control=list(minimal_branch_len=20,
                                              euclidean_distance_ratio=10))
#
plot_cells(mycds, reduction_method="UMAP",
           color_cells_by = "celltype", 
           label_groups_by_cluster=FALSE,
           label_leaves=T, 
           label_branch_points=TRUE,
           graph_label_size=3,group_label_size=3,cell_size=0.5)
#===============================================================================
mycds1 <- mycds
mycds1 <- order_cells(mycds1)
#
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 4, 
           cell_size=0.5, 
           trajectory_graph_segment_size = 1)

#===============================================================================
#
pd <- pseudotime(mycds1, reduction_method = 'UMAP')
pd <- pseudotime(mycds1, reduction_method = 'tSNE')
CD4T <- AddMetaData(CD4T,metadata = pd,col.name = 'pseudotime')
#
mycds1_res <- graph_test(mycds1, 
                         neighbor_graph="principal_graph", cores=4)
#-------------------------------------------------------------------------------
mycds2 <- mycds1
mycds2_res <- subset_pr_test_res
res_ids <- row.names(subset(mycds2_res, q_value < 0.01))
gene_module_df <- find_gene_modules(mycds2[res_ids,], 
                                    resolution=c(10^seq(-6,-1)))
write.csv(gene_module_df, file = "gene_module_df.csv")
cell_group_df <- tibble::tibble(cell=row.names(colData(mycds2)), 
                                cell_group=colData(mycds2)$celltype)
agg_mat <- aggregate_gene_expression(mycds2, 
                                     gene_module_df, 
                                     cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

#-------------------------------------------------------------------------------
#
genes_sig <- mycds2_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#
plot_genes_in_pseudotime(mycds2[genes_sig,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)
#
plot_genes_in_pseudotime(mycds2[genes_sig,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2)
#
plot_cells(mycds2, genes=genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

#-------------------------------------------------------------------------------
library(viridis)
plot_cells(mycds2,
           genes=c("CCR7","GZMK",'TOX',"TOX2"), 
           reduction_method="tSNE",
           label_cell_groups=F,
           show_trajectory_graph=T, 
           cell_size=1, trajectory_graph_color = "black", 
           label_branch_points = F, 
           label_roots = F, label_leaves = F)+
  scale_color_viridis(option="inferno")

#-------------------------------------------------------------------------------
#
ica_space_df <- t(mycds2@principal_graph_aux$UMAP$dp_mst) %>% as.data.frame() %>% dplyr::mutate(sample_name = rownames(.),sample_state = rownames(.))
pgraph <- principal_graph(mycds2)$UMAP
edge_df <- pgraph %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="UMAP_1",
                       source_prin_graph_dim_2="UMAP_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="UMAP_1",
                       target_prin_graph_dim_2="UMAP_2"),
                   by = "target")

#
trajectory_graph_segment_size = 1
trajectory_graph_color = "red"
#
custom.plot.data <- data.frame(reducedDims(mycds2)[["UMAP"]], 
                               condition = factor(mycds2@colData$group, levels = c("CAYA", "ADULT")),
                               celltype = mycds2 @colData$celltype,
                               pseudotime = pseudotime(mycds2, reduction_method = "UMAP"))
colnames(custom.plot.data)[1:2] <- c("UMAP1", "UMAP2")
#
library(ggplot2)
ggplot(custom.plot.data, aes(x = UMAP1, y = UMAP2, colour = condition)) +
  geom_point(size = 0.5) +
  theme_classic() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
               color = I(trajectory_graph_color), linetype = "solid", na.rm = TRUE, data = edge_df) +
  guides(color = guide_legend(title = 'group', 
                              override.aes = list(size = 4)))+
  scale_colour_manual(values = c('#CC3333',"#006699"))
#=================================================================================
#
library(Hmisc)
library(dplyr)
mycds_order <- mycds2
pData(mycds_order)$pseudotime <- pseudotime(mycds_order)
mycds_order <- mycds_order[,is.finite(pseudotime(mycds_order))]
a <- as.numeric(mycds_order$pseudotime)

mycds_order$pseudotime_bin <-  cut2(as.numeric(mycds_order$pseudotime), 
                                    seq(min(a),max(a), by = 0.15)) 
table(mycds_order@colData@listData[["pseudotime_bin"]])
table(mycds_order@colData@listData[["orig.ident"]],
      mycds_order@colData@listData[["pseudotime_bin"]])
bins <- levels(mycds_order$pseudotime_bin)
proportion_df <- data.frame()
for(i in 1:length(bins)){
  bin = bins[i]
  meta <- subset(as.data.frame(mycds_order@colData@listData), pseudotime_bin == bin)
  if(nrow(meta) == 0) {
    df <- data.frame(group = c("CAYA", "ADULT"),
                     proportion = c(0,0),
                     bin_num = c(i, i),
                     bin= c(bins[i], bins[i]))
    proportion_df <- rbind(proportion_df, df)
  }
  
  else{
    df <- data.frame(table(meta$group) / nrow(meta))
    df <- df %>% dplyr::rename(c(group = Var1, proportion = Freq))
    df$bin_num <- i
    df$bin <- bin
    proportion_df <- rbind(proportion_df, df)
  }
  
}

#
proportion_cor <- cor.test(
  subset(proportion_df[,], group == 'CAYA') %>% .$bin_num,
  subset(proportion_df[,], group == 'CAYA') %>% .$proportion,
  method='pearson')

#
library(ggpubr)
p <- ggscatter(
  subset(proportion_df, group == 'CAYA'),
  color = "#D43D4F",
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of CAYA cells \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')
p
##------------------------------------------------------------------------------
proportion_cor <- cor.test(
  subset(proportion_df[,], group == 'ADULT') %>% .$bin_num,
  subset(proportion_df[,], group == 'ADULT') %>% .$proportion,
  method='pearson')

ggscatter(
  subset(proportion_df, group == 'ADULT'),
  color = "#D43D4F",
  x = 'bin_num', y = 'proportion',
  add ='reg.line',
  add.params = list(color = '#5c4b4d', fill = 'lightgray'),
  conf.int=TRUE,
  cor.coef=TRUE,
  cor.coeff.args = list(method = 'pearson', label.sep = '\n'),
  alpha = 0.8,
  size = 2.5, 
  title = "Propotion of ADULT cells \n along the pseudotime") + 
  xlab('pseudotime') +
  ylab('proportion')

#===============================================================================
pseudotime_table <- subset(mycds2_res, q_value < 0.001)
rownames(pseudotime_table) <- pseudotime_table$gene_short_name
pseudotime_genes <- as.character(rownames(pseudotime_table))
CD4T_sec <- CD4T[, rownames(mycds_order@colData)]
CD4T_sec$pseudotime <- mycds_order@colData$pseudotime
CD4T_sec$pseudotime_bin <- mycds_order@colData$pseudotime_bin
#
Idents(CD4T_sec) <- CD4T_sec$pseudotime_bin
average_exp <- AverageExpression(CD4T_sec, slot='data', assay='RNA', features=pseudotime_genes)
average_exp <- average_exp$RNA
#
library(scales)
library(dplyr)
scaled_exp <- sapply(1:nrow(average_exp), function(i) rescale(as.numeric(average_exp[i,]), to=c(-1,1))) %>% t %>% as.data.frame
rownames(scaled_exp) <- rownames(average_exp)
average_exp <- scaled_exp
write.csv(average_exp, file = 'pseudotime_bin_average_exp.csv')

library(future.apply)

range01 <- function(x){
  cur <- average_exp[x,]
  (cur-min(cur))/(max(cur)-min(cur))
}
scaled <- lapply(rownames(average_exp), range01)
scaled <- do.call(rbind, scaled)
rownames(scaled) <- rownames(average_exp)

ordering <- future_lapply(1:nrow(scaled), function(i){
  match(names(scaled[i,])[scaled[i,] >= 0.99][1], colnames(scaled))
})
ordered <- scaled[rownames(scaled)[order(unlist(ordering))],]
ordered <- ordered[rownames(ordered),]

library(viridis)
install_github("jokergoo/ComplexHeatmap")
install.packages('magick')
library(magick)
library(ComplexHeatmap)
Heatmap(as.matrix(ordered),
        show_column_names = F, 
        show_row_names=F,
        col = viridis(256),
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        use_raster = TRUE)
#============================================================================================

cutColumn_Means <- function(data_exp,
                            cut){
  plot_matrix_combin <- list()
  nums <- ncol(data_exp)/cut
  if (nums-round(nums, 0)==0){
    
    for (i in 1:length(seq(1, ncol(data_exp), cut))){
      num <- seq(1, ncol(data_exp), cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
      
    }
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
    
  }else{
    
    for (i in 1:length(seq(1, ncol(data_exp)-cut, cut))){
      num <- seq(1, ncol(data_exp)-cut, cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
    }
    
    plot_matrix_combin[[length(seq(1, ncol(data_exp)-cut, cut))+1]] <- as.data.frame(rowMeans(data_exp[,(max(num)+cut):ncol(data_exp)]))                       
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
  }
  
  
}

plot_test <- cutColumn_Means(plot_matrix,cut = 25)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


p1 <- pheatmap::pheatmap(plot_matrix_combin, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=4,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


p2 <- pheatmap::pheatmap(plot_test, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=4,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


#
annotation_row <- data.frame(Cluster=factor(cutree(p2$tree_row, 4)))
row.names(annotation_row) <- rownames(plot_test)

rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
names(rowcolor) <- c("1","2","3","4") 
#
ann_colors <- list(Cluster=rowcolor)

p3 <- pheatmap::pheatmap(plot_test, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         cutree_rows=4,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")

#
gene <- c('ISG15', 'GZMK', 'GZMA', 'CCR7', 'CCL4', 'CCL5', 'EOMES', 'TOX', 'TOX2',
          'CXCL13', 'RORC', 'RORA', ' CD40LG', 'IL2RB', 'XCL1', 'XCL5')
#
p4 <- pheatmap::pheatmap(plot_test, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=T, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")
source('C:/Users/guoka/Desktop/add.flag.R')
add.flag(p4,kept.labels = gene,repel.degree = 0.2)

#=======================================================================================
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

module_gene <- as.data.frame(cutree(p3$tree_row, k=4))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
#
Module_GO=data.frame()

for (i in unique(module_gene$Module)) {
  
  data=filter(module_gene,module_gene$Module==i)
  df=bitr(data$gene, 
          fromType="SYMBOL",
          toType=c("ENTREZID"), 
          OrgDb="org.Hs.eg.db")
  
  go <- enrichGO(gene= unique(df$ENTREZID),
                 OrgDb= org.Mm.eg.db,
                 keyType= 'ENTREZID',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}

#
Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
Module_GO <- Module_GO[,c("ID","Description","qvalue","cluster")]
write.csv(Module_GO, file = 'Module_GO.csv')

#=======================================================================================

peuGene_exp <- normalized_counts(mycds2, norm_method = "log")
peuGene_exp <- as.data.frame(peuGene_exp)
dim(peuGene_exp)

pseudotime <- as.data.frame(pseudotime(mycds2))
colnames(pseudotime) <- "pseudotime"
pseudotime <- pseudotime[rowSums(abs(pseudotime)==Inf)==0,,drop=F]
peuGene_exp <- peuGene_exp[,colnames(peuGene_exp) %in% rownames(pseudotime)]
dim(peuGene_exp)

peuGene_exp_CAYA <- peuGene_exp[,grep(pattern="^CAYA",colnames(peuGene_exp))]
pseudotime_CAYA <- pseudotime[colnames(peuGene_exp_CAYA),]

peuGene_exp_ADULT <- peuGene_exp[,grep(pattern="^ADULT",colnames(peuGene_exp))]
pseudotime_ADULT <- pseudotime[colnames(peuGene_exp_ADULT),]

CAYA_anxa1 <- peuGene_exp_CAYA['ANXA1',]
CAYA_anxa1 <- as.data.frame(t(CAYA_anxa1))
CAYA_anxa1$pseudotime <- pseudotime_CAYA
CAYA_anxa1$group <- 'CAYA'

ADULT_anxa1 <- peuGene_exp_ADULT['ANXA1',]
ADULT_anxa1 <- as.data.frame(t(ADULT_anxa1))
ADULT_anxa1$pseudotime <- pseudotime_ADULT
ADULT_anxa1$group <- 'ADULT'

peu_trend <- rbind(CAYA_anxa1, ADULT_anxa1)
#ggplot
ggplot(peu_trend, aes(x=pseudotime, y=ANXA1, color=group))+
  geom_smooth(aes(fill=group))+ 
  xlab('pseudotime') + 
  ylab('Relative expression') +
  ggtitle('Anxa1')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14))+
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))

CAYA_tox2 <- peuGene_exp_CAYA['TOX2',]
CAYA_tox2 <- as.data.frame(t(CAYA_tox2))
CAYA_tox2$pseudotime <- pseudotime_CAYA
CAYA_tox2$group <- 'CAYA'

ADULT_tox2 <- peuGene_exp_ADULT['TOX2',]
ADULT_tox2 <- as.data.frame(t(ADULT_tox2))
ADULT_tox2$pseudotime <- pseudotime_ADULT
ADULT_tox2$group <- 'ADULT'
peu_trend <- rbind(CAYA_tox2, ADULT_tox2)
#ggplot
ggplot(peu_trend, aes(x=pseudotime, y=TOX2, color=group))+
  geom_smooth(aes(fill=group))+
  xlab('pseudotime') + 
  ylab('Relative expression') +
  ggtitle('tox2')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14))+
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))