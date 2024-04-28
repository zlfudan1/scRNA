load(file="./NK-celltype.Rdata")
DimPlot(NK, reduction = "tsne", label=T) 
DimPlot(NK, reduction = "umap", group.by='celltype',label = T,raster=FALSE)
DimPlot(NK, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
#-------------------------------------------------------------------------------
#
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 50   
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP50<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP50,'NEW_NK_R0.3_TopMarkersTOP50.csv', row.names=F)
#-------------------------------------------------------------------------------
setwd("./NK")
load(file="./4celltypeharmonyNK-0.1.Rdata")
DimPlot(NK, reduction = "tsne", label=T) 
DimPlot(NK, reduction = "umap", group.by='celltype',label = T,raster=FALSE, split.by = 'group' )
DimPlot(NK, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
#
df_NK <- NK@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = NK@meta.data$celltype)%>%
  cbind(seurat_clusters = NK@meta.data$seurat_clusters)
colnames(df_NK)
#
cluster_order <- c('0','1',"2","3", "4", "5")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', "#D8A326"),
                        c('0','1',"2","3", "4", "5"))
#
cell <- df_NK %>%group_by(seurat_clusters) %>%
  summarise(tsne_1 = median(UMAP_1),
            tsne_2 = median(UMAP_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:5)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_NK)
ggplot(df_NK,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+2, yend = min(UMAP_2)),
               colour = "black", size=0.5, arrow = arrow(length = unit(0.2,"cm"), type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+2),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), type = "closed")) +
  annotate("text", x = min(df_NK$UMAP_1) +1, y = min(df_NK$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_NK$UMAP_1) -0.5, y = min(df_NK$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 4,angle=90)+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("NK_UMAP.pdf",device = cairo_pdf, width = 6, height = 4.5, dpi = 300)

#
B <- A
B$x <- 1
B$lab <- c(4:0)
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
ggsave("NK_umap_legend.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
#-------------------------------------------------------------------------------

table(NK$HT)
prop.table(table(NK$celltype))
Cellratio_HT <- prop.table(table(NK$celltype, NK$HT), margin = 2)
Cellratio_HT <- as.data.frame(Cellratio_HT)
write.csv(Cellratio_HT, "Cellratio_HT_NK.csv", row.names = T)
#
table(NK$group)
prop.table(table(NK$celltype))
Cellratio_group <- prop.table(table(NK$celltype, NK$group), margin = 2)
Cellratio_group <- as.data.frame(Cellratio_group)
write.csv(Cellratio_group, "Cellratio_group_NK.csv", row.names = T)
#
table(NK$spatialdistr)
prop.table(table(NK$celltype))
Cellratio_spatialdistr <- prop.table(table(NK$celltype, NK$spatialdistr), margin = 2)
Cellratio_spatialdistr <- as.data.frame(Cellratio_spatialdistr)
write.csv(Cellratio_spatialdistr, "Cellratio_spatialdistr_NK.csv", row.names = T)
#
Idents(NK)="metastasis"
table(NK@meta.data$metastasis)
NK=RenameIdents(NK,c('N0'='N0','N1a'='N1','N1b'='N1'))
NK@meta.data$metastasis2=NK@active.ident
Idents(NK)="metastasis2"
table(NK@meta.data$metastasis2)
#
table(NK$metastasis2)
prop.table(table(NK$celltype))
Cellratio_metastasis2 <- prop.table(table(NK$celltype, NK$metastasis2), margin = 2)
Cellratio_metastasis2 <- as.data.frame(Cellratio_metastasis2)
write.csv(Cellratio_metastasis2, "Cellratio_metastasis2_NK.csv", row.names = T)

data <- read.csv("cellratio_NK.csv", header = T)
data$var <- factor(data$var, levels=c('NKT_GNB2L1', 'CD16NK_GZMB', 'CD16NK_XCL1', 'NKT_RACK1', 'CD56NK'))
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
  scale_fill_manual(values = c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', "#D8A326"))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 8.5), linetype= "dashed", linewidth = 0.5)
ggsave("cellratio_NK.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)