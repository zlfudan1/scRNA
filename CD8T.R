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
setwd("./CD8+T")
load(file="./scRNA-umap-celltype.Rdata")
table(scRNA@meta.data$seurat_clusters)
Idents(scRNA)="seurat_clusters"
sc.CD8Tcell=subset(scRNA,ident=c("1","2","5","24"))
table(sc.CD8Tcell@meta.data$orig.ident)
sce.list2 <- SplitObject(sc.CD8Tcell, split.by = "orig.ident")
for (i in 1:length(sce.list2)) {
  sce.list2[[i]] <- NormalizeData(sce.list2[[i]],verbose=FALSE)
  sce.list2[[i]] <- FindVariableFeatures(sce.list2[[i]], selection.method = "vst",nfeatures = 2000,verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = sce.list2)
scRNA.anchors <- FindIntegrationAnchors(object.list = sce.list2, anchor.features = features)
CD8T <- IntegrateData(anchorset = scRNA.anchors)
DefaultAssay(CD8T) <- "RNA"
#
CD8T <- NormalizeData(CD8T)
CD8T <- FindVariableFeatures(CD8T, selection.method = "vst", nfeatures = 2000)
#
DefaultAssay(CD8T) <- "integrated"
CD8T<- ScaleData(CD8T, features = VariableFeatures(CD8T))
#
CD8T <- RunPCA(CD8T, features = VariableFeatures(CD8T))
ElbowPlot(CD8T)
CD8T <- FindNeighbors(CD8T, dims = 1:20)
save(CD8T,file = './222CD8T.FindNeighbors.Rdata')
load(file = './222CD8T.FindNeighbors.Rdata')
CD8T<- FindClusters(CD8T, resolution = 0.4) 
pc.num=1:20  
#tSNE
CD8T = RunTSNE(CD8T, dims = pc.num)
CD8T <- RunUMAP(CD8T, dims = pc.num)
#group_by
DimPlot(CD8T, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(CD8T, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
DefaultAssay(CD8T) <- "RNA"
#
diff.wilcox = FindAllMarkers(CD8T,logfc.threshold = 0.25,only.pos = T,test.use = "MAST")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "CD8T_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "CD8T_R0.2_top50_diff_genes_wilcox.csv", row.names = F)
#
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50 
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP50<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP50,'NEW_CD8T_R0.4_TopMarkersTOP50.csv', row.names=F)
#
Idents(CD8T)="seurat_clusters"
CD8T2 <- subset(CD8T,ident=c('0','1','2',"3", '4', '5', '6', '7', '9'))
DefaultAssay(CD8T2) <- "RNA"
#
CD8T2 <- NormalizeData(CD8T2)
CD8T2 <- FindVariableFeatures(CD8T2, selection.method = "vst", nfeatures = 2000)
#
DefaultAssay(CD8T2) <- "integrated"
CD8T2<- ScaleData(CD8T2, features = VariableFeatures(CD8T2))
#
CD8T2 <- RunPCA(CD8T2, features = VariableFeatures(CD8T2))
ElbowPlot(CD8T2)
CD8T2 <- FindNeighbors(CD8T2, dims = 1:20)
CD8T2<- FindClusters(CD8T2, resolution = 0.4) 
pc.num=1:20   
#tSNE
CD8T2 = RunTSNE(CD8T2, dims = pc.num)
CD8T2 <- RunUMAP(CD8T2, dims = pc.num)
#group_by
DimPlot(CD8T2, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(CD8T2, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
DefaultAssay(CD8T2) <- "RNA"
#
diff.wilcox = FindAllMarkers(CD8T2,logfc.threshold = 0.25,only.pos = T,test.use = "MAST")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "CD8T_diff_genes_wilcox.csv", row.names = F)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "CD8T_R0.2_top50_diff_genes_wilcox.csv", row.names = F)
#
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP50<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP50,'NEW_CD8T2_R0.4_deleted_TopMarkersTOP50.csv', row.names=F)
save(CD8T2,file = './final_CD8T2.Rdata')
load(file = './final_CD8T2.Rdata')
Idents(CD8T2)="seurat_clusters"
CD8T2=RenameIdents(CD8T2,c('0'='CD8T_Tn_CCR7','1'='CD8T_Tem_GZMK','2'='CD8T_Tm_CCL5','3'='CD8T_Trm_XCL1','4'='CD8T_Tex_CXCL13', '5'='CD8T_Tn_CD55', '6'='CD8T_Tn_PLCG2', '7'='CD8T_Tc17_KLRB1', '8'='CD8T_Tex_TOX2', '9'='CD8T_ISG+T'))
CD8T2@meta.data$celltype=CD8T2@active.ident
Idents(CD8T2)="celltype"
table(CD8T2@meta.data$celltype)
DimPlot(CD8T2, reduction = "tsne", group.by='celltype',label = T,raster=FALSE)
CD8T <- CD8T2
M <-CD8T
M@meta.data <- M@meta.data[,-(8:15)]
M@meta.data <- M@meta.data[,-(9:10)]
M@meta.data <- M@meta.data[,-10]
CD8T@meta.data <- M@meta.data
DimPlot(CD8T, reduction = "tsne", group.by='celltype',label = T,raster=FALSE)
save(CD8T,file = './final_CD8T.Rdata')
###-----------------------------------------------------------------------------
#
setwd("./CD8T")
load(file = './final_CD8T.Rdata')
DimPlot(CD8T, reduction = "umap", group.by='celltype',label = T,raster=FALSE)
df_CD8T <- CD8T@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = CD8T@meta.data$celltype)%>%
  cbind(seurat_clusters = CD8T@meta.data$seurat_clusters)
colnames(df_CD8T)
#
cluster_order <- c('0','1',"2","3", "4", "5", "6", "7", "8", "9")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                          '#D8A326', '#9E732B', '#656565','#A4C7DA', '#3172A5'),
                        c('0','1',"2","3", "4","5","6", "7", "8", "9"))
#
cell <- df_CD8T %>%group_by(seurat_clusters) %>%
  summarise(umap_1 = median(UMAP_1),
            umap_2 = median(UMAP_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:9)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","umap_1", "umap_2" ,"ID")
colnames(df_CD8T)
ggplot(df_CD8T,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+7, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+8),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), 
                                                        type = "closed")) +
  annotate("text", x = min(df_CD8T$UMAP_1) +3, y = min(df_CD8T$UMAP_2) -0.8, label = "UMAP_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_CD8T$UMAP_1) -0.8, y = min(df_CD8T$UMAP_2) + 4, label = "UMAP_2",
           color="black",size = 4,angle=90)+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 
#-------------------------------------------------------------------------------
cell <- df_CD8T %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:9)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_CD8T)
ggplot(df_CD8T,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+1.5, yend = min(UMAP_2)),
               colour = "black", size=0.5,linewidth = 1,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,linewidth = 1, arrow = arrow(length = unit(0.2,"cm"), 
                                                                       type = "closed")) +
  annotate("text", x = min(df_CD8T$UMAP_1) + 0.7, y = min(df_CD8T$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 3, fontface = 'bold', family = 'serif') + 
  annotate("text", x = min(df_CD8T$UMAP_1) -0.5, y = min(df_CD8T$UMAP_2) + 0.7, label = "UMAP_2",
           color="black",size = 3, fontface = 'bold', family = 'serif', angle = 90) +
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
ggsave("CD8T_umap.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)

#
B <- A
B$x <- 1
B$lab <- c(9:0)
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
ggsave("CD8T_umap_legend.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
#-------------------------------------------------------------------------------
#
table(CD8T$HT)
prop.table(table(CD8T$celltype))
Cellratio_HT <- prop.table(table(CD8T$celltype, CD8T$HT), margin = 2)
Cellratio_HT <- as.data.frame(Cellratio_HT)
write.csv(Cellratio_HT, "Cellratio_HT.csv", row.names = T)
#
table(CD8T$group)
prop.table(table(CD8T$celltype))
Cellratio_group <- prop.table(table(CD8T$celltype, CD8T$group), margin = 2)
Cellratio_group <- as.data.frame(Cellratio_group)
write.csv(Cellratio_group, "Cellratio_group.csv", row.names = T)
#
table(CD8T$metastasis)
prop.table(table(CD8T$celltype))
Cellratio_metastasis <- prop.table(table(CD8T$celltype, CD8T$metastasis), margin = 2)
Cellratio_metastasis <- as.data.frame(Cellratio_metastasis)
write.csv(Cellratio_metastasis, "Cellratio_metastasis.csv", row.names = T)
#
table(CD8T$spatialdistr)
prop.table(table(CD8T$celltype))
Cellratio_spatialdistr <- prop.table(table(CD8T$celltype, CD8T$spatialdistr), margin = 2)
Cellratio_spatialdistr <- as.data.frame(Cellratio_spatialdistr)
write.csv(Cellratio_spatialdistr, "Cellratio_spatialdistr.csv", row.names = T)

setwd("./CD8T")
data <- read.csv("cellratio_CD8T.csv", header = T)
View(data)
data$var <- factor(data$var, levels=c('CD8T_Tn_CCR7','CD8T_Tem_GZMK','CD8T_Tm_CCL5','CD8T_Trm_XCL1','CD8T_Tex_CXCL13', 'CD8T_Tn_CD55', 'CD8T_Tn_PLCG2', 'CD8T_Tc17_KLRB1', 'CD8T_Tex_TOX2', 'CD8T_ISG+T'))
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
                               '#D8A326', '#9E732B', '#656565','#A4C7DA', '#3172A5' ))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 9.5), linetype= "dashed", linewidth = 0.5)
ggsave("cellratio.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
#
Idents(CD8T)="metastasis"
table(CD8T@meta.data$metastasis)
CD8T=RenameIdents(CD8T,c('N0'='N0','N1a'='N1','N1b'='N1'))
CD8T@meta.data$metastasis2=CD8T@active.ident
Idents(CD8T)="metastasis2"
table(CD8T@meta.data$metastasis2)
#
table(CD8T$metastasis2)
prop.table(table(CD8T$celltype))
Cellratio_metastasis2 <- prop.table(table(CD8T$celltype, CD8T$metastasis2), margin = 2)
Cellratio_metastasis2 <- as.data.frame(Cellratio_metastasis2)
write.csv(Cellratio_metastasis2, "Cellratio_metastasis2.csv", row.names = T)
data <- read.csv("cellratio_CD8T_2.csv", header = T)
data$var <- factor(data$var, levels=c('CD8T_Tn_CCR7','CD8T_Tem_GZMK','CD8T_Tm_CCL5','CD8T_Trm_XCL1','CD8T_Tex_CXCL13', 'CD8T_Tn_CD55', 'CD8T_Tn_PLCG2', 'CD8T_Tc17_KLRB1', 'CD8T_Tex_TOX2', 'CD8T_ISG+T'))
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
                               '#D8A326', '#9E732B', '#656565','#A4C7DA', '#3172A5'))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 8.5), linetype= "dashed", linewidth = 0.5)
ggsave("cellratio_2.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
#-------------------------------------------------------------------------------
#Ucell
library(UCell)
library(stringr)
library(ggplot2)
library(viridis)

markers <- list()
markers$CD8T_Tn_CCR7 <- c('LEF1',	'TCF7',	'KLF2',	'LTB',	'CCR7',	'IL7R',	'SELL',	'MAL',	'EEF1A1',	'TRABD2A',	'TPT1',	'EEF1B2',	'NOSIP',	'PABPC1')
markers$CD8T_Tem_GZMK <- c('EOMES',	'GZMK',	'CD74',	'DUSP2',	'CMC1',	'CST7',	'SH2D1A',	'CRTAM',	'GZMA',	'CCL5',	'IL32',	'CCL4',	'HLA-DRB1',	'HLA-DPA1',	'COTL1',	'HLA-DPB1',	'HLA-DQA1',	'HLA-DQB1',	'HLA-DRB5',	'ITM2C',	'APOBEC3G',	'GZMH',	'GZMM',	'NKG7',	'PLEK')
markers$CD8T_Tm_CCL5 <- c('ZFP36L2',	'ZFP36',	'ANXA1',	'LMNA',	'FTH1',	'RGCC',	'S100A4')
markers$CD8T_Trm_XCL1 <- c('ZNF683',	'HOPX',	'ID2',	'ZFP36L2',	'CKLF',	'IL32',	'XCL1',	'CCL5',	'XCL2',	'CXCR3',	'CAPG',	'S100A4',	'LGALS1',	'ITGA1',	'LDLRAD4')
markers$CD8T_Tex_CXCL13 <- c('TOX',	'CXCL13',	'GZMA',	'CCL3',	'CCL5',	'CCL4',	'IFNG',	'CD74',	'CXCR6',	'HAVCR2',	'CD27',	'LYST',	'PDCD1',	'DUSP4',	'CTLA4',	'TNFRSF9',	'HLA-DQA1',	'HLA-DRB1',	'RBPJ',	'FAM3C',	'GZMB',	'KRT86',	'TIGIT',	'PHLDA1')
markers$CD8T_Tn_CD55 <- c('LEF1',	'TCF7',	'KLF2',	'TXK',	'BACH2',	'CCR7',	'IL7R',	'SELL',	'EEF1A1',	'ACTN1',	'TRABD2A',	'EEF1B2',	'NELL2',	'NOSIP',	'PABPC1')
markers$CD8T_Tn_PLCG2 <- c('LEF1',	'TCF7',	'TXK',	'LTB',	'CCR7',	'IL7R',	'SELL',	'MAL',	'ACTN1',	'TRABD2A',	'EEF1B2',	'NOSIP')
markers$CD8T_Tc17_KLRB1 <- c('ZBTB16',	'CEBPD',	'RORA',	'TNF',	'CCR6',	'IL7R',	'IL18RAP',	'KLRB1',	'NCR3')
markers$CD8T_Tex_TOX2 <- c('TOX',	'PRDM1',	'GZMK',	'CXCL13',	'CD74',	'CD27',	'PDCD1',	'DUSP4',	'CTLA4',	'TOX2',	'GEM',	'PDCD1',	'GNG4',	'TNFRSF4',	'NMB')
markers$CD8T_ISG_T <- c('PLSCR1',	'STAT1',	'IRF7',	'SP100',	'TNFSF10',	'IFIT1',	'RSAD2',	'IFIT3',	'IFI44L',	'MX1',	'IFI6',	'OAS1',	'CMPK2',	'ISG15',	'OAS3')
#
marker_score <- AddModuleScore_UCell(CD8T,
                                     features=markers)
#
a1 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tn_CCR7_UCell")
p1 <- FeaturePlot(marker_score,features = a1, cols = viridis(256))+
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tn_CCR7')+
  theme(plot.title = element_text(hjust = 0.5))

a2 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tem_GZMK_UCell")
p2 <- FeaturePlot(marker_score,features = a2, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tem_GZMK')+
  theme(plot.title = element_text(hjust = 0.5))

a3 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tm_CCL5_UCell")
p3 <- FeaturePlot(marker_score,features = a3, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tm_CCL5')+
  theme(plot.title = element_text(hjust = 0.5))

a4 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Trm_XCL1_UCell")
p4 <- FeaturePlot(marker_score,features = a4, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Trm_XCL1')+
  theme(plot.title = element_text(hjust = 0.5))

a5 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tex_CXCL13_UCell")
p5 <- FeaturePlot(marker_score,features = a5, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tex_CXCL13')+
  theme(plot.title = element_text(hjust = 0.5))

a6 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tn_CD55_UCell")
p6 <- FeaturePlot(marker_score,features = a6, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tn_CD55')+
  theme(plot.title = element_text(hjust = 0.5))

a7 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tn_PLCG2_UCell")
p7 <- FeaturePlot(marker_score,features = a7, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tn_PLCG2')+
  theme(plot.title = element_text(hjust = 0.5))

a8 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tc17_KLRB1_UCell")
p8 <- FeaturePlot(marker_score,features = a8, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tc17_KLRB1')+
  theme(plot.title = element_text(hjust = 0.5))

a9 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_Tex_TOX2_UCell")
p9 <- FeaturePlot(marker_score,features = a9, cols = viridis(256))+ 
  theme_void()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_Tex_TOX2')+
  theme(plot.title = element_text(hjust = 0.5))

a10 <- colnames(marker_score@meta.data) %>% str_subset("CD8T_ISG_T_UCell")
p10 <- FeaturePlot(marker_score,features = a10, cols = viridis(256))+ 
  theme_void()+
  theme(plot.title = element_text(size = 10))+
  ggtitle('CD8T_ISG_T')+
  theme(plot.title = element_text(hjust = 0.5))

install.packages('gridExtra')
library(gridExtra)
p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2)
p
ggsave("CD8T_Ucell.pdf",device = cairo_pdf, width = 18, height = 7, dpi = 300, p)
###-----------------------------------------------------------------------------
marker <- c('PDCD1', 'TOX', 'CXCL13', 'TIGIT', 'CTLA4', 'TNFRSF9', 'HAVCR2', 'LAG3',
            'CST7', 'GZMK', 'GZMA', 'NKG7', 'IFNG', 'PRF1', 'GZMB', 'GNLY', 
            'CD28', 'EOMES', 'CCR5', 'CCL4', 'CCL5', 'CD200R1', 'TMEM155', 'CD27', 
            'TOX2', 'TSHZ2', 'BATF', 'GEM', 'CD200', 'SNX9', 'ENTPD1', 'LAYN', 'ETV1', 'PRDM1', 'CSF1', 'MYO7A', 'MYO1E', 'GNG4', 
            'BTLA', 'FOXP3', 'TNFRSF4', 'TCF7', 'EEF1A1', 'SELL', 'CCR7', 'IL6R', 'IGFBP4', 'IGFL2')
marker <- c('TOX', 'PRDM1', 'GZMK', 'CXCL13', 'GZMA', 'CCL3', 'CCL5', 'CCL4', 'IFNG', 'CD74', 'CXCR6', 'HAVCR2', 'CD27', 'LYST', 'PDCD1', 'DUSP4', 'CTLA4', 'TNFRSF9', 'HLA-DQA1', 'HLA-DRB1', 'RBPJ', 'TOX2', 'GZMB', 'IL2RB', 'IL2RG', 'LAYN', 'ENTPD1', 'KRT86', 'TNFRSF18', 'GEM', 'TIGIT')
marker <- c('TOX', 'TOX2', 'CD74', 'CD27', 'PDCD1', 'DUSP4', 'CTLA4', 'TNFRSF18', 'TIGIT', 
            'BATF', 'SNX9', 'PRDM1', 
            'CXCL13', 'GZMA', 'CCL3', 'CCL4', 'CCL5', 'IFNG', 'CXCR6', 'HAVCR2', 'LYST',  'TNFRSF9', 'HLA-DQA1', 'HLA-DRB1')
CD8T_Tex = CD8T[, Idents(CD8T) %in% c('CD8T_Tex_CXCL13', 'CD8T_Tex_TOX2')]
DotPlot(CD8T_Tex, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
#
p <- DotPlot(CD8T_Tex, features = marker)
exp <- p$data
library(ggplot2)
p1 <- ggplot(exp,aes(x=features.plot,y=id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust = 1,vjust=0.5))+
  scale_color_gradient(low="lightgrey",high="blue")+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
#
library(ggtree)
df1<-reshape2::dcast(exp,id~features.plot,value.var = "avg.exp.scaled")
rownames(df1)<-df1[,1] 
df1 <- df1[,-1]
df1 <- t(df1)
head(df1)
p4<-ggtree(hclust(dist(df1)))
p4+geom_tiplab()
#
p3 <- ggtree(hclust(dist(t(df1))))+layout_dendrogram()
p3
#
library(aplot)
p1%>%insert_left(p4)
p1%>%insert_top(p3,height = 0.15)%>%insert_left(p4,width = 0.3)
library(patchwork)
#-------------------------------------------------------------------------------
install.packages("devtools")
devtools::install_github("Simon-Leonard/FlexDotPlot")
library(FlexDotPlot)
library(Seurat)
library(ggplot2)
dp = DotPlot(CD8T_Tex, features = marker) + RotatedAxis()
dot_plot(dp$data[,c(3,4,1,2,5)],
         size_var = "pct.exp", 
         col_var = "avg.exp.scaled",
         size_legend = "Percent Expressed", 
         col_legend = "Average Expression",
         x.lab.pos = "bottom", 
         y.lab.pos = "right",
         display_max_sizes = F, 
         shape.scale=8,
         hclust_method = "ward.D2",
         dend_x_var = c("pct.exp", "avg.exp.scaled"),
         dend_y_var = c("pct.exp", "avg.exp.scaled"),
         text.size = 0.5,
         text.vjust = 0.5,
         size.breaks.number=6,
         color.breaks.number=4,
         x.lab.size.factor = 0.8,
         cols.use = c('#330066','#336699','#66CC66','#FFCC33'),
         y.lab.size.factor = 1)

#
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
p <- DotPlot(CD8T_Tex, features = marker,
             cols = my36colors, group.by = "celltype")+coord_flip()
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)


p2 <- ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#6699CC','#FFFF99','#CC3333'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  geom_vline(xintercept=c(1.5), linetype="dotted",size=1)+
  geom_rect(aes(xmin = 1 - 0.5, xmax = 1 + 0.5, ymin = 13 - 0.5, ymax = 24 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 1 + 0.5, xmax = 2 + 0.5, ymin = 1 - 0.5, ymax =13 - 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(legend.direction = "horizontal",legend.position = "bottom")

library(aplot)
p2%>%
  insert_left(p4,width = 0.3)%>%insert_top(p3,height = 0.3)
####----------------------------------------------------------------------------
#monocle3
data <- GetAssayData(CD8T, assay = 'RNA', slot = 'counts')
cell_metadata <- CD8T@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#
cds <- preprocess_cds(cds, num_dim = 50)
#umap
cds <- reduce_dimension(cds,reduction_method='UMAP',preprocess_method = 'PCA')
cds <- reduce_dimension(cds,reduction_method='tSNE',preprocess_method = 'PCA')

plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") + ggtitle('cds.tsne')
#
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(CD8T, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds,use_partition = FALSE)
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

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5)
p2 <- plot_cells(cds,
                 color_cells_by = "celltype",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5)
p1 | p2
plot_cells(cds,reduction_method = "UMAP",
           color_cells_by = "pseudotime",cell_size = 0.6,
           trajectory_graph_segment_size = 1,
           label_branch_points = F,
           label_roots = T,
           label_leaves = F,
           alpha = 2,
           cell_stroke = T,
           trajectory_graph_color = "black") + 
  theme_classic()

cds@colData$celltype <- factor(cds@colData$celltype)
#
col <- ggsci::scale_color_lancet()
col <- col$palette(length(unique(CD8T$celltype)))
plot_cells(cds,reduction_method = "UMAP",
           color_cells_by = "celltype",
           cell_size = 0.6,
           trajectory_graph_segment_size = 1,
           label_branch_points = F,label_roots = T,
           label_leaves = F,alpha = 2,cell_stroke = T,
           trajectory_graph_color = "black",
           label_cell_groups = T,group_label_size = 4)  + 
  theme_bw() + theme_classic() + scale_color_manual(values = col)
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                          '#D8A326', '#9E732B', '#656565','#A4C7DA', '#3172A5'),
                        c('CD8T_Tn_CCR7','CD8T_Tem_GZMK','CD8T_Tm_CCL5','CD8T_Trm_XCL1',
                          'CD8T_Tex_CXCL13', 'CD8T_Tn_CD55', 'CD8T_Tn_PLCG2', 'CD8T_Tc17_KLRB1', 
                          'CD8T_Tex_TOX2', 'CD8T_ISG+T'))

p1 <- plot_cells(cds,reduction_method = "UMAP",
                 color_cells_by = "celltype",
                 cell_size = 0.6,
                 trajectory_graph_segment_size = 1,
                 label_branch_points = F,label_roots = T,
                 label_leaves = F,alpha = 3,cell_stroke = T,
                 trajectory_graph_color = "black",
                 label_cell_groups = F, group_label_size = 3, graph_label_size = 3)  +
  theme_classic() + 
  theme(axis.text = element_blank(),  axis.title = element_blank(), axis.ticks = element_blank()) + 
  theme(axis.line = element_blank()) +
  theme(legend.position = "none")+
  scale_color_manual(values = col_cluster) + 
  ggtitle("") + #colored by celltype
  theme(plot.title = element_text(hjust = 0.5)) 

p2 <- plot_cells(cds,reduction_method = "UMAP",
                 color_cells_by = "pseudotime",
                 cell_size = 0.6,
                 trajectory_graph_segment_size = 1,
                 label_branch_points = F,label_roots = T,
                 label_leaves = F,alpha = 3,cell_stroke = T,
                 trajectory_graph_color = "black",
                 label_cell_groups = F,graph_label_size = 3)  +
  theme_classic() + 
  theme(axis.text = element_blank(),  axis.title = element_blank(), axis.ticks = element_blank()) + 
  theme(axis.line = element_blank()) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))
p1 | p2
p3 <- p1 + p2
ggsave("monocle3_pseudotime_CD8T.pdf",device = cairo_pdf, width = 11, height = 10, dpi = 300, p2)
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
#
plot_cells(cds,reduction_method="tSNE", genes=top10, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
plot_cells(cds,reduction_method="UMAP", genes=top10, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
###-----------------------------------------------------------------------------
mycds <- cds
plot_cells(mycds, reduction_method="UMAP",
           color_cells_by = "celltype", 
           label_groups_by_cluster=FALSE,
           label_leaves=T, 
           label_branch_points=TRUE,
           graph_label_size=3,group_label_size=3,cell_size=0.5)

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
#-------------------------------------------------------------------------------
#
pd <- pseudotime(mycds1, reduction_method = 'UMAP')
CD8T <- AddMetaData(CD8T,metadata = pd,col.name = 'pseudotime')
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
p1 <- pheatmap::pheatmap(agg_mat,
                         scale="column", clustering_method="ward.D2")
ggsave("CD8T_Module.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p1)
#-------------------------------------------------------------------------------
genes_sig <- mycds2_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(mycds2[genes_sig,], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)
plot_genes_in_pseudotime(mycds2[genes_sig,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2)

plot_cells(mycds2, genes=genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

library(viridis)
plot_cells(mycds2,
           genes=c("CCR7","GZMK",'TOX',"TOX2"), 
           reduction_method="UMAP",
           label_cell_groups=F,
           show_trajectory_graph=T, 
           cell_size=1, trajectory_graph_color = "black", 
           label_branch_points = F, 
           label_roots = F, label_leaves = F)+
  scale_color_viridis(option="inferno")

#-------------------------------------------------------------------------------
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

trajectory_graph_segment_size = 1
trajectory_graph_color = "red"
#提取数据
custom.plot.data <- data.frame(reducedDims(mycds2)[["UMAP"]], 
                               condition = factor(mycds2@colData$group, levels = c("CAYA", "ADULT")),
                               celltype = mycds2 @colData$celltype,
                               pseudotime = pseudotime(mycds2, reduction_method = "UMAP")
)
colnames(custom.plot.data)[1:2] <- c("UMAP1", "UMAP2")

library(ggplot2)
p1 <- ggplot(custom.plot.data, aes(x = UMAP1, y = UMAP2, colour = condition)) +
  geom_point(size = 0.5) +
  theme_classic() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
               color = I(trajectory_graph_color), linetype = "solid", na.rm = TRUE, data = edge_df) +
  guides(color = guide_legend(title = 'group', 
                              override.aes = list(size = 4)))+
  scale_colour_manual(values = c('#CC3333',"#006699")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))
#-------------------------------------------------------------------------------

library(Hmisc)
library(dplyr)
mycds_order <- mycds2
pData(mycds_order)$pseudotime <- pseudotime(mycds_order)
mycds_order <- mycds_order[,is.finite(pseudotime(mycds_order))]
a <- as.numeric(mycds_order$pseudotime)#拟时时间轴

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

proportion_cor <- cor.test(
  subset(proportion_df[,], group == 'CAYA') %>% .$bin_num,
  subset(proportion_df[,], group == 'CAYA') %>% .$proportion,
  method='pearson')

library(ggpubr)
p1 <- ggscatter(
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
p1
##-------------
proportion_cor <- cor.test(
  subset(proportion_df[,], group == 'ADULT') %>% .$bin_num,
  subset(proportion_df[,], group == 'ADULT') %>% .$proportion,
  method='pearson')
p2 <- ggscatter(
  subset(proportion_df, group == 'ADULT'),
  color = "#006699",
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
p1 | p2
p3 <- p1 + p2
ggsave("Propotion_cells_along_pseudotime_CD8T.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p3)
#-------------------------------------------------------------------------------

markers$CD8T_Tn_CCR7 <- c('LEF1',	'TCF7',	'KLF2',	'LTB',	'CCR7',	'IL7R',	'SELL',	'MAL',	'EEF1A1',	'TRABD2A',	'TPT1',	'EEF1B2',	'NOSIP',	'PABPC1')
markers$CD8T_Tem_GZMK <- c('EOMES',	'GZMK',	'CD74',	'DUSP2',	'CMC1',	'CST7',	'SH2D1A',	'CRTAM',	'GZMA',	'CCL5',	'IL32',	'CCL4',	'HLA-DRB1',	'HLA-DPA1',	'COTL1',	'HLA-DPB1',	'HLA-DQA1',	'HLA-DQB1',	'HLA-DRB5',	'ITM2C',	'APOBEC3G',	'GZMH',	'GZMM',	'NKG7',	'PLEK')
markers$CD8T_Tm_CCL5 <- c('ZFP36L2',	'ZFP36',	'ANXA1',	'LMNA',	'FTH1',	'RGCC',	'S100A4')
markers$CD8T_Trm_XCL1 <- c('ZNF683',	'HOPX',	'ID2',	'ZFP36L2',	'CKLF',	'IL32',	'XCL1',	'CCL5',	'XCL2',	'CXCR3',	'CAPG',	'S100A4',	'LGALS1',	'ITGA1',	'LDLRAD4')
markers$CD8T_Tex_CXCL13 <- c('TOX',	'CXCL13',	'GZMA',	'CCL3',	'CCL5',	'CCL4',	'IFNG',	'CD74',	'CXCR6',	'HAVCR2',	'CD27',	'LYST',	'PDCD1',	'DUSP4',	'CTLA4',	'TNFRSF9',	'HLA-DQA1',	'HLA-DRB1',	'RBPJ',	'FAM3C',	'GZMB',	'KRT86',	'TIGIT',	'PHLDA1')
markers$CD8T_Tn_CD55 <- c('LEF1',	'TCF7',	'KLF2',	'TXK',	'BACH2',	'CCR7',	'IL7R',	'SELL',	'EEF1A1',	'ACTN1',	'TRABD2A',	'EEF1B2',	'NELL2',	'NOSIP',	'PABPC1')
markers$CD8T_Tn_PLCG2 <- c('LEF1',	'TCF7',	'TXK',	'LTB',	'CCR7',	'IL7R',	'SELL',	'MAL',	'ACTN1',	'TRABD2A',	'EEF1B2',	'NOSIP')
markers$CD8T_Tc17_KLRB1 <- c('ZBTB16',	'CEBPD',	'RORA',	'TNF',	'CCR6',	'IL7R',	'IL18RAP',	'KLRB1',	'NCR3')
markers$CD8T_Tex_TOX2 <- c('TOX',	'PRDM1',	'GZMK',	'CXCL13',	'CD74',	'CD27',	'PDCD1',	'DUSP4',	'CTLA4',	'TOX2',	'GEM',	'BATF',	'GNG4',	'TNFRSF4',	'NMB')
markers$CD8T_ISG_T <- c('PLSCR1',	'STAT1',	'IRF7',	'SP100',	'TNFSF10',	'IFIT1',	'RSAD2',	'IFIT3',	'IFI44L',	'MX1',	'IFI6',	'OAS1',	'CMPK2',	'ISG15',	'OAS3')

#CD8T_Tn_CCR7
PMN_state_gene <- c('LEF1',	'TCF7',	'KLF2',	'LTB',	'CCR7',	'IL7R',	'SELL',	'MAL',	'EEF1A1',	'TRABD2A',	'TPT1',	'EEF1B2',	'NOSIP',	'PABPC1')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p1 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tn_CCR7')

#---------------
#CD8T_Tem_GZMK
PMN_state_gene <- c('EOMES',	'GZMK',	'CD74',	'DUSP2',	'CMC1',	'CST7',	'SH2D1A',	'CRTAM',	'GZMA',	'CCL5',	'IL32',	'CCL4',	'HLA-DRB1',	'HLA-DPA1',	'COTL1',	'HLA-DPB1',	'HLA-DQA1',	'HLA-DQB1',	'HLA-DRB5',	'ITM2C',	'APOBEC3G',	'GZMH',	'GZMM',	'NKG7',	'PLEK')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p2 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tem_GZMK')
#---------------
#CD8T_Tm_CCL5
PMN_state_gene <- c('ZFP36L2',	'ZFP36',	'ANXA1',	'LMNA',	'FTH1',	'RGCC',	'S100A4')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p3 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tm_CCL5')
#---------------
#CD8T_Trm_XCL1
PMN_state_gene <- c('ZNF683',	'HOPX',	'ID2',	'ZFP36L2',	'CKLF',	'IL32',	'XCL1',	'CCL5',	'XCL2',	'CXCR3',	'CAPG',	'S100A4',	'LGALS1',	'ITGA1',	'LDLRAD4')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p4 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Trm_XCL1')
#---------------
#CD8T_Tex_CXCL13
PMN_state_gene <- c('TOX',	'CXCL13',	'GZMA',	'CCL3',	'CCL5',	'CCL4',	'IFNG',	'CD74',	'CXCR6',	'HAVCR2',	'CD27',	'LYST',	'PDCD1',	'DUSP4',	'CTLA4',	'TNFRSF9',	'HLA-DQA1',	'HLA-DRB1',	'RBPJ',	'FAM3C',	'GZMB',	'KRT86',	'TIGIT',	'PHLDA1')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p5 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tex_CXCL13')
#---------------
#CD8T_Tn_CD55
PMN_state_gene <- c('LEF1',	'TCF7',	'KLF2',	'TXK',	'BACH2',	'CCR7',	'IL7R',	'SELL',	'EEF1A1',	'ACTN1',	'TRABD2A',	'EEF1B2',	'NELL2',	'NOSIP',	'PABPC1')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p6 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tn_CD55')
#---------------
#CD8T_Tn_PLCG2
PMN_state_gene <- c('LEF1',	'TCF7',	'TXK',	'LTB',	'CCR7',	'IL7R',	'SELL',	'MAL',	'ACTN1',	'TRABD2A',	'EEF1B2',	'NOSIP')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p7 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tn_PLCG2')
#---------------
#CD8T_Tc17_KLRB1
PMN_state_gene <- c('ZBTB16',	'CEBPD',	'RORA',	'TNF',	'CCR6',	'IL7R',	'IL18RAP',	'KLRB1',	'NCR3')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p8 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tc17_KLRB1')
#---------------
#CD8T_Tex_TOX2
PMN_state_gene <- c('TOX',	'PRDM1',	'GZMK',	'CXCL13',	'CD74',	'CD27',	'PDCD1',	'DUSP4',	'CTLA4',	'TOX2',	'GEM',	'BATF',	'GNG4',	'TNFRSF4',	'NMB')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p9 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_Tex_TOX2')
#---------------
#CD8T_ISG+T
PMN_state_gene <- c('PLSCR1',	'STAT1',	'IRF7',	'SP100',	'TNFSF10',	'IFIT1',	'RSAD2',	'IFIT3',	'IFI44L',	'MX1',	'IFI6',	'OAS1',	'CMPK2',	'ISG15',	'OAS3')
CD8T_sc <- AddModuleScore(CD8T_sc,
                          features=list('PI' = PMN_state_gene),
                          pool = rownames(CD8T_sc), k=F, nbin=24,
                          name = c('PI'))
modules <- select(CD8T_sc@meta.data, c(pseudotime_bin, PI1))
modules$pseudotime_bin_num <- as.numeric(CD8T_sc$pseudotime_bin)
features <- c('PI')
tmp <- lapply(1:max(modules$pseudotime_bin_num), function(i){
  df <- modules %>% subset(pseudotime_bin_num == i)
  data.frame(
    value = as.numeric(sum(df[,'PI1']) / nrow(df)),
    bin_num = i,
    feature = features)
}) 
plot_df <- Reduce(rbind, tmp)
plot_df <- na.omit(plot_df)
p10 <- ggplot(plot_df, aes(bin_num, value)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', alpha = 0.75) +
  geom_point(size = 2) +
  geom_smooth() +
  xlab('pseudotime') + ylab('module score') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ggtitle('CD8T_ISG+T')

library(gridExtra)
p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2)
p
ggsave("clusters_along_the_pseudotime_CD8T.pdf",device = cairo_pdf, width = 11, height = 5, dpi = 300, p)

#-------------------------------------------------------------------------------
pseudotime_table <- subset(mycds2_res, q_value < 0.001)
rownames(pseudotime_table) <- pseudotime_table$gene_short_name
pseudotime_genes <- as.character(rownames(pseudotime_table))

CD8T_sec <- CD8T[, rownames(mycds_order@colData)]
CD8T_sec$pseudotime <- mycds_order@colData$pseudotime
CD8T_sec$pseudotime_bin <- mycds_order@colData$pseudotime_bin
#
Idents(CD8T_sec) <- CD8T_sec$pseudotime_bin
average_exp <- AverageExpression(CD8T_sec, slot='data', assay='RNA', features=pseudotime_genes)
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

#-------------------------------------------------------------------------------
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(monocle3)
genes <- row.names(subset(mycds2_res, q_value< 0.01 & morans_I > 0.1))

genes <- genes[!grepl("^RP[SL]", genes, ignore.case = F)]
genes <- genes[!grepl("^MT-", genes, ignore.case = F)]
plot_matrix <- exprs(mycds2)[match(genes,
                                   rownames(rowData(mycds2))),
                             order(pseudotime(mycds2))]

plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes;
dim(plot_matrix)


plot_matrix_combin <- list()
for (i in 1:length(seq(1, 34050-50, 50))){
  num <- seq(1, 34050-50, 50)
  A <- plot_matrix[,num[i]:(50+num[i]-1)]
  a <- rowMeans(A)
  a <- as.data.frame(a)
  a <- a$a
  plot_matrix_combin[[i]] <- a
}

length(plot_matrix_combin)
plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
rownames(plot_matrix_combin) <- rownames(plot_matrix)
write.csv(plot_matrix_combin, 'plot_matrix_combin.csv')
#============================================================================================

cutColumn_Means <- function(data_exp,
                            cut
){
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
write.csv(plot_test, 'plot_test.csv')

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
                         cutree_rows=3,
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
                         cutree_rows=3,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)


#
annotation_row <- data.frame(Cluster=factor(cutree(p2$tree_row, 3)))
row.names(annotation_row) <- rownames(plot_test)

rowcolor <- c("#85B22E","#E29827","#922927") 
names(rowcolor) <- c("1","2","3") 
#
ann_colors <- list(Cluster=rowcolor)

p3 <- pheatmap::pheatmap(plot_test, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         cutree_rows=3,
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
gene <- c('LEF1',	'TCF7',	'KLF2',	'CCR7',	'IL7R',	'EOMES',	'GZMK',	'GZMA',	'CCL5',	'IL32',	'GZMH',	'GZMM',	'NKG7',	'ANXA1',	'FTH1',	'S100A4',	'HOPX',	'XCL1',	'XCL2',	'CXCR3', 'TOX',	'CXCL13',	'IFNG',	'CXCR6',	'HAVCR2',	'PDCD1',	'CTLA4',	'TIGIT',	'NCR3',	'TOX2',	'BATF',	'IRF7',	'IFIT1',	'IFI6',	'OAS1')

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
p5 <- add.flag(p4,kept.labels = gene,repel.degree = 0.2)
ggsave("gene_clusters_along_the_pseudotime_CD8T.pdf",device = cairo_pdf, width = 8, height = 10, dpi = 300, p5)


#-------------------------------------------------------------------------------
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db", force = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)

module_gene <- as.data.frame(cutree(p3$tree_row, k=3))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)

Module_GO=data.frame()

for (i in unique(module_gene$Module)) {
  
  data=filter(module_gene,module_gene$Module==i)
  df=bitr(data$gene, 
          fromType="SYMBOL",
          toType=c("ENTREZID"), 
          OrgDb="org.Hs.eg.db")
  
  go <- enrichGO(gene= unique(df$ENTREZID),
                 OrgDb= org.Hs.eg.db,
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

Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
write.csv(Module_GO, file = 'Module_GO_all.csv')
Module_GO <- Module_GO[,c("ID","Description","qvalue","cluster")]
write.csv(Module_GO, file = 'Module_GO.csv')
#-------------------------------------------------------------------------------
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

#PD-1 good
CAYA_PDCD1 <- peuGene_exp_CAYA['PDCD1',]
CAYA_PDCD1 <- as.data.frame(t(CAYA_PDCD1))
CAYA_PDCD1$pseudotime <- pseudotime_CAYA
CAYA_PDCD1$group <- 'CAYA'
ADULT_PDCD1 <- peuGene_exp_ADULT['PDCD1',]
ADULT_PDCD1 <- as.data.frame(t(ADULT_PDCD1))
ADULT_PDCD1$pseudotime <- pseudotime_ADULT
ADULT_PDCD1$group <- 'ADULT'
peu_trend <- rbind(CAYA_PDCD1, ADULT_PDCD1)
p1 <- ggplot(peu_trend, aes(x=pseudotime, y=PDCD1, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('PDCD1')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))

#LAG3 good
CAYA_LAG3 <- peuGene_exp_CAYA['LAG3',]
CAYA_LAG3 <- as.data.frame(t(CAYA_LAG3))
CAYA_LAG3$pseudotime <- pseudotime_CAYA
CAYA_LAG3$group <- 'CAYA'
ADULT_LAG3 <- peuGene_exp_ADULT['LAG3',]
ADULT_LAG3 <- as.data.frame(t(ADULT_LAG3))
ADULT_LAG3$pseudotime <- pseudotime_ADULT
ADULT_LAG3$group <- 'ADULT'
peu_trend <- rbind(CAYA_LAG3, ADULT_LAG3)
p2 <- ggplot(peu_trend, aes(x=pseudotime, y=LAG3, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('LAG3')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))

#NFATC1(NFAT2) good
CAYA_NFATC1 <- peuGene_exp_CAYA['NFATC1',]
CAYA_NFATC1 <- as.data.frame(t(CAYA_NFATC1))
CAYA_NFATC1$pseudotime <- pseudotime_CAYA
CAYA_NFATC1$group <- 'CAYA'
ADULT_NFATC1 <- peuGene_exp_ADULT['NFATC1',]
ADULT_NFATC1 <- as.data.frame(t(ADULT_NFATC1))
ADULT_NFATC1$pseudotime <- pseudotime_ADULT
ADULT_NFATC1$group <- 'ADULT'
peu_trend <- rbind(CAYA_NFATC1, ADULT_NFATC1)
p3 <- ggplot(peu_trend, aes(x=pseudotime, y=NFATC1, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('NFATC1')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))

#CXCR5 good
CAYA_CXCR5 <- peuGene_exp_CAYA['CXCR5',]
CAYA_CXCR5 <- as.data.frame(t(CAYA_CXCR5))
CAYA_CXCR5$pseudotime <- pseudotime_CAYA
CAYA_CXCR5$group <- 'CAYA'
ADULT_CXCR5 <- peuGene_exp_ADULT['CXCR5',]
ADULT_CXCR5 <- as.data.frame(t(ADULT_CXCR5))
ADULT_CXCR5$pseudotime <- pseudotime_ADULT
ADULT_CXCR5$group <- 'ADULT'
peu_trend <- rbind(CAYA_CXCR5, ADULT_CXCR5)
p4 <- ggplot(peu_trend, aes(x=pseudotime, y=CXCR5, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('CXCR5')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))

#IL-2 good
CAYA_IL2 <- peuGene_exp_CAYA['IL2',]
CAYA_IL2 <- as.data.frame(t(CAYA_IL2))
CAYA_IL2$pseudotime <- pseudotime_CAYA
CAYA_IL2$group <- 'CAYA'
ADULT_IL2 <- peuGene_exp_ADULT['IL2',]
ADULT_IL2 <- as.data.frame(t(ADULT_IL2))
ADULT_IL2$pseudotime <- pseudotime_ADULT
ADULT_IL2$group <- 'ADULT'
peu_trend <- rbind(CAYA_IL2, ADULT_IL2)
p5 <- ggplot(peu_trend, aes(x=pseudotime, y=IL2, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('IL2')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))

#IFNG (IFNγ) good
CAYA_IFNG <- peuGene_exp_CAYA['IFNG',]
CAYA_IFNG <- as.data.frame(t(CAYA_IFNG))
CAYA_IFNG$pseudotime <- pseudotime_CAYA
CAYA_IFNG$group <- 'CAYA'
ADULT_IFNG <- peuGene_exp_ADULT['IFNG',]
ADULT_IFNG <- as.data.frame(t(ADULT_IFNG))
ADULT_IFNG$pseudotime <- pseudotime_ADULT
ADULT_IFNG$group <- 'ADULT'
peu_trend <- rbind(CAYA_IFNG, ADULT_IFNG)
p6 <- ggplot(peu_trend, aes(x=pseudotime, y=IFNG, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('IFNG')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))

#HAVCR2 (TIM3) good
CAYA_HAVCR2 <- peuGene_exp_CAYA['HAVCR2',]
CAYA_HAVCR2 <- as.data.frame(t(CAYA_HAVCR2))
CAYA_HAVCR2$pseudotime <- pseudotime_CAYA
CAYA_HAVCR2$group <- 'CAYA'
ADULT_HAVCR2 <- peuGene_exp_ADULT['HAVCR2',]
ADULT_HAVCR2 <- as.data.frame(t(ADULT_HAVCR2))
ADULT_HAVCR2$pseudotime <- pseudotime_ADULT
ADULT_HAVCR2$group <- 'ADULT'
peu_trend <- rbind(CAYA_HAVCR2, ADULT_HAVCR2)
p7 <- ggplot(peu_trend, aes(x=pseudotime, y=HAVCR2, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('HAVCR2')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))
#TOX2 good
CAYA_TOX2 <- peuGene_exp_CAYA['TOX2',]
CAYA_TOX2 <- as.data.frame(t(CAYA_TOX2))
CAYA_TOX2$pseudotime <- pseudotime_CAYA
CAYA_TOX2$group <- 'CAYA'
ADULT_TOX2 <- peuGene_exp_ADULT['TOX2',]
ADULT_TOX2 <- as.data.frame(t(ADULT_TOX2))
ADULT_TOX2$pseudotime <- pseudotime_ADULT
ADULT_TOX2$group <- 'ADULT'
peu_trend <- rbind(CAYA_TOX2, ADULT_TOX2)
p8 <- ggplot(peu_trend, aes(x=pseudotime, y=TOX2, color=group))+
  geom_smooth(aes(fill=group))+
  ggtitle('TOX2')+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.title = element_blank(),
        legend.position = 'none')+
  scale_color_manual(name=NULL, values = c("#006699","#CC3333"))+
  scale_fill_manual(name=NULL, values = c("#006699","#CC3333"))
p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
p
ggsave("Tex_keyfactors.pdf",device = cairo_pdf, width = 10, height = 11, dpi = 300, p)
###-----------------------------------------------------------------------------

devtools::install_github("Japrin/sscVis",dependencies = T)
devtools::install_github("Japrin/STARTRAC")
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
setwd("./Startrac")
read_tcr <- function(tcrfile){
  p3_n <- read.csv(tcrfile)
  p3_n$file <- sub('.filtered_contig_annotations.csv','',sub('^.*/','',tcrfile))
  return(p3_n)
}

tcrfiles <- list.files('./','.filtered_contig_annotations.csv',full.names = T)
tcrfiles

if (all(file.exists(tcrfiles))){
  tcr_list = list()
  for (i in 1:length(tcrfiles)){
    print(i)
    tcr_list[[i]] = read_tcr(tcrfile = tcrfiles[i])
  }
}
lapply(tcr_list,  dim)

vdj <- do.call(rbind, tcr_list);dim(vdj)
head(vdj,2)
table(vdj$file)

vdj <- vdj %>% 
  dplyr::filter(high_confidence =="true" & 
                  chain %in% c("TRA","TRB") &
                  productive =="true")
vdj$Cell_name <- paste0(vdj$file,'_',vdj$barcode)
head(vdj,2)
#
vdj_a <- vdj %>% filter(chain =="TRA") %>% dplyr::arrange(desc(umis), desc(reads)) 
vdj_b <- vdj %>% filter(chain =="TRB") %>% dplyr::arrange(desc(umis), desc(reads)) 

test <- vdj_a %>% 
  dplyr::group_by(Cell_name) %>% 
  dplyr::summarise(reads=max(reads), umis=max(umis)) 
head(test)
vdj_a <- data.frame(inner_join(vdj_a, test)) 

dim(vdj_a)

test <- vdj_b %>% group_by(Cell_name) %>%
  dplyr::summarise(reads = max(reads), umis=max(umis) )
vdj_b <- data.frame(inner_join(vdj_b, test))
dim(vdj_b)

final_vdj = dplyr::full_join(x = vdj_a, y=vdj_b, by = c("Cell_name"), suffix = c(".TRA",".TRB"))
dim(final_vdj)
head(final_vdj,2)

load(file="./final_CD8T.Rdata")
Idents(CD8T)="orig.ident"
table(CD8T@meta.data[["orig.ident"]])
TCRCD8T=subset(CD8T,ident=c("ADULT1P","ADULT1T","CAYA1P","CAYA1T","CAYA2P","CAYA2T","CAYA3P","CAYA3T",
                            "CAYA4P","CAYA4T","CAYA5P","CAYA5T"))

subT <- TCRCD8T
subT@meta.data <- subT@meta.data %>% 
  mutate(Cell_name = rownames(subT@meta.data)) %>% 
  inner_join(final_vdj, by = "Cell_name")

head(subT@meta.data)

subT@meta.data$Clone_AA = paste(subT@meta.data$cdr3.TRA, subT@meta.data$cdr3.TRB, sep="_")

subT@meta.data = subset(subT@meta.data, productive.TRA == "true" & productive.TRB == "true"  ) ; dim(subT@meta.data)
subT@meta.data = subT@meta.data %>% arrange(., Clone_AA)

tmp = subT@meta.data %>% 
  group_by(Clone_AA) %>%
  summarize(Clone_NUM = n()) %>%
  mutate(Clone_ID = paste0("Clone_",rownames(.)))
head(tmp)
subT@meta.data = merge.data.frame(subT@meta.data, tmp) 
head(subT@meta.data,2)

subT.meta <- subT@meta.data %>% 
  select(Cell_name, Clone_ID, orig.ident,Clone_NUM, group,celltype, seurat_clusters, spatialdistr)
head(subT.meta)

subT.meta$clone.status <- ifelse(subT.meta$Clone_NUM >1 ,"Clonal","NoClonal")

subT.meta <- subT.meta %>% 
  select(Cell_name, Clone_ID, clone.status, orig.ident,Clone_NUM, group,celltype, seurat_clusters, spatialdistr)
head(subT.meta)
colnames(subT.meta) <- c('Cell_Name', 'clone.id', 'clone.status', 'patient', 'Clone_NUM', 'group', 'majorCluster', 'seurat_cluster', 'loc')
head(subT.meta)
View(subT.meta)

write.table (subT.meta, file ="subTmeta.txt", sep ="\t", row.names =F, col.names =TRUE, quote =TRUE)

in.dat <- read.table("subTmeta.txt", header=TRUE, stringsAsFactors = F)
in.dat
source("./func.R")
tic('Startrac.run')
out <- Startrac.run(in.dat, proj="CRC", cores=NULL,verbose=F)
View(out@pIndex.migr)

plot(out,index.type="cluster.all",byPatient=F)

plot(out,index.type="cluster.all",byPatient=T)

plot(out,index.type="pairwise.migr",byPatient=F)

plot(out,index.type="pairwise.tran",byPatient=T)

#-------------------------------------------------------------------------------
#miloR
devtools::install_github("MarioniLab/miloR", ref="devel") 
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(ggbeeswarm)
library(scater)
library(scales)
library(forcats)
library(data.table)
library(stringr)
library(dplyr)

DimPlot(CD8T,label = T)
unique(CD8T$group)

CD8T$group<- factor(CD8T$group, levels=c("CAYA","ADULT"))

CD8T <- as.SingleCellExperiment(CD8T)
CD8T_milo <- miloR::Milo(CD8T)
CD8T_milo <- miloR::buildGraph(CD8T_milo, k = 10, d = 30)

CD8T_milo <- makeNhoods(CD8T_milo, 
                        prop = 0.2, 
                        k = 30, 
                        d=50, 
                        refined = TRUE)

CD8T_milo <- countCells(CD8T_milo, meta.data = data.frame(colData(CD8T_milo)), 
                        sample="orig.ident")

traj_design <- data.frame(colData(CD8T_milo))[,c("orig.ident", "group")]
traj_design$orig.ident <- as.factor(traj_design$orig.ident)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident

CD8T_milo <- calcNhoodDistance(CD8T_milo, d=50)

da_results <- testNhoods(CD8T_milo, 
                         design = ~ group, 
                         design.df = traj_design)
CD8T_milo <- buildNhoodGraph(CD8T_milo)

#-------------------------------------------------------------------------------
#UMAP
plotUMAP(CD8T_milo, colour_by = "celltype")

pdf('milo_PMN.pdf', width=8, height=8)
plotNhoodGraphDA(CD8T_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) 
#-------------------------------------------------------------------------------
#
da_results <- annotateNhoods(CD8T_milo, da_results, coldata_col = "celltype")
plotDAbeeswarm(da_results, group.by = "celltype") +
  scale_color_gradient2(low="#070091",
                        mid="lightgrey",
                        high="#910000",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black'))
#-------------------------------------------------------------------------------
#scenic
library(ggheatmap)
library(reshape2)
library(ComplexHeatmap)
rss_data <- read.csv('./CD8T_rssdata.csv', header = T, row.names = 1)
View(rss_data)
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)

library(data.table)
rss_data <- fread('CD8T_rssdata(1).csv', header = T)
rss_data <- as.data.frame(rss_data)
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]

df <- data.frame(colnames(rss_data))
colnames(df) <- 'celltype'

top_anno = HeatmapAnnotation(df = df,
                             border = F,
                             show_annotation_name = F,
                             gp = gpar(col = 'white'),
                             col = list(celltype = c("CD8T_Tn_CCR7"="#349573",
                                                     "CD8T_Tem_GZMK"="#C7591A",
                                                     "CD8T_Tm_CCL5"="#6F6BA5",
                                                     "CD8T_Trm_XCL1"="#CE2B7D",
                                                     "CD8T_Tex_CXCL13"="#6A9D3A",
                                                     "CD8T_Tn_CD55"="#D8A326",
                                                     "CD8T_Tn_PLCG2"="#9E732B",
                                                     "CD8T_Tc17_KLRB1"="#656565",
                                                     "CD8T_Tex_TOX2"="#A4C7DA",
                                                     "CD8T_ISG+T"="#3172A5")))#颜色设置

marker_exp <- t(scale(t(rss_data),scale = T,center = T))

TF <- read.csv('CD8T_rssdata.csv', header = T, row.names = 1)
TF <- TF$x
pdf('CD4T_scenic.pdf', width = 9, height = 8)
p <- Heatmap(marker_exp,
             cluster_rows = T,
             cluster_columns = T,
             show_column_names = F,
             show_row_names = T,
             show_row_dend= F,
             show_column_dend= F,
             row_order = TF,
             column_title = NULL,
             heatmap_legend_param = list(
               title=' '),
             col = colorRampPalette(c("#1a5592","white","#b83d3d"))(100),
             border = 'white',
             rect_gp = gpar(col = "white", lwd = 1),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 10),
             top_annotation = top_anno)
p
dev.off()
ggsave("CD4T_scenic.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)