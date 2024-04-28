load(file="scRNA-umap-celltype.Rdata")
sce=subset(scRNA,ident=c("Fibroblasts"))
table(sce@meta.data$orig.ident)
sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                          "percent.mt", "nCount_RNA"), verbose = T)
sce <- RunPCA(sce, features = VariableFeatures(object = sce),, npcs = 50)
library(harmony)
sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:25, reduction = "harmony")
#DimPlot(sce,reduction = "umap",label=T,split.by  = 'orig.ident') 
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 
set.resolutions = c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1,1.2)
sce  <- FindClusters(object = sce , resolution = set.resolutions, verbose = FALSE) 
for (res in c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1,1.2)) {
  sce=FindClusters(sce, #graph.name = "CCA_snn", 
                   resolution = res, algorithm = 1)
}
clustree(sce)
colnames(sce@meta.data)

sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25)
sce  <- FindClusters(object = sce , resolution = 0.2, verbose = FALSE) 
Fibroblasts <- RunUMAP(sce, reduction = "harmony", dims = 1:25)
Fibroblasts <- RunTSNE(Fibroblasts, reduction = "harmony", dims =1:25)#

DimPlot(Fibroblasts, reduction = "umap",label = T) 

scRNA_harmony <- Fibroblasts
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "4harmonyFibroblasts-0.2_diff_genes_wilcox.csv", row.names = F)
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 25  
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'4harmonyFibroblasts-0.2_TopMarkersTOP25.csv', row.names=F)
####################################################
#############################################################
Idents(Fibroblasts)="seurat_clusters"
Fibroblasts3 <- subset(Fibroblasts2,ident=c('1','0','2',"3"))

sce <- Fibroblasts3

sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                          "percent.mt", "nCount_RNA"), verbose = T)
sce <- RunPCA(sce, features = VariableFeatures(object = sce),, npcs = 50)
sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:25, reduction = "harmony")
#DimPlot(sce,reduction = "umap",label=T,split.by  = 'orig.ident') 
sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 
set.resolutions = c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1,1.2)
sce  <- FindClusters(object = sce , resolution = set.resolutions, verbose = FALSE) 
# 
for (res in c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1,1.2)) {
  sce=FindClusters(sce, #graph.name = "CCA_snn", 
                   resolution = res, algorithm = 1)
}
clustree(sce)
colnames(sce@meta.data)

sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25)
sce  <- FindClusters(object = sce , resolution = 0.1, verbose = FALSE) 
Fibroblasts3 <- RunUMAP(sce, reduction = "harmony", dims = 1:25)
Fibroblasts3 <- RunTSNE(Fibroblasts3, reduction = "harmony", dims =1:25)#
DimPlot(Fibroblasts3, reduction = "umap",label = T) 
DimPlot(Fibroblasts3, reduction = "tsne",label = T) 
DimPlot(Fibroblasts3, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE,split.by = "group")

scRNA_harmony <- Fibroblasts3
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "4harmonyFibroblasts3-0.1_diff_genes_wilcox.csv", row.names = F)
ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 25  
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'4harmonyFibroblasts3-0.1_TopMarkersTOP25.csv', row.names=F)

M <- Fibroblasts3
M@meta.data <- M@meta.data[,-(11:15)]
M@meta.data <- M@meta.data[,-(12:31)]
Fibroblasts3@meta.data <- M@meta.data

VlnPlot(Fibroblasts, 
        features = c('CFD',"C3",'FBLN1','LAMP5','MYH11','CD36','FAP',"SMMHC"),
        pt.size = 0,
        ncol = 2)
FeaturePlot(Fibroblasts3,features = c('CFD',"C3",'FBLN1','LAMP5','MYH11','CD36','FAP',"POSTN"),
            pt.size = 0,
            ncol = 2)

Idents(Fibroblasts3)="seurat_clusters"
table(Fibroblasts3@meta.data$seurat_clusters)
table(Fibroblasts3@meta.data$orig.ident)
Fibroblasts3=RenameIdents(Fibroblasts3,c('0'='lpmCAF_CD36','1'='iCAF_CFD','2'='myoCAF_MYH11','3'='emCAF_LAMP5'))
Fibroblasts3@meta.data$celltype=Fibroblasts3@active.ident
Idents(Fibroblasts3)="celltype"
table(Fibroblasts3@meta.data$celltype)
DimPlot(Fibroblasts3, reduction = "umap", group.by='celltype',label = T,raster=FALSE,split.by = "group")
Fibroblasts <- Fibroblasts3

Fibroblasts2=Fibroblasts
Fibroblasts2$group2<-Fibroblasts2$orig.ident
table(Fibroblasts2$orig.ident)
Fibroblasts2$group2<-as.character(Fibroblasts2$group2)
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA1P'] <- 'CATA-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA1T'] <- 'CAYA-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA2P'] <- "CATA-P"
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA2T'] <- 'CAYA-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA3P'] <- 'CATA-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA3T'] <- 'CAYA-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA4P'] <- 'CATA-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA4T'] <- 'CAYA-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA5P'] <- 'CATA-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'CAYA5T'] <- 'CAYA-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT1P'] <- 'ADUL-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT1T'] <- 'ADUL-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT2T'] <- 'ADUL-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT2P'] <- 'ADUL-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT3T'] <- 'ADUL-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT3P'] <- 'ADUL-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT4T'] <- 'ADUL-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT4P'] <- 'ADUL-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT5T'] <- 'ADUL-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT5P'] <- 'ADUL-P'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT6T'] <- 'ADUL-T'
Fibroblasts2$group2[Fibroblasts2$group2 == 'ADULT6P'] <- 'ADUL-P'
Fibroblasts=Fibroblasts2
save(Fibroblasts,all.markers,file="4celltypeharmonyFibroblasts-0.1.Rdata") 
#########################################################################
library(ggplot2)
library(dplyr)
library(tidyverse)
library(Seurat)
library(devtools)
#devtools::install_github("sajuukLyu/ggunchull", type = "source")
library(ggunchull)
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
library(ggrepel)
DimPlot(Fibroblasts, reduction = "umap",label = T) 
DimPlot(Fibroblasts, reduction = "tsne",label = T) 
DimPlot(Fibroblasts, reduction = "umap", group.by='celltype',label = T,raster=FALSE,split.by = "group")
levels = c("mCAF_LAMP5","mCAF_CD36","mCAF_STEAP4", "iCAF_CFD","mCAF_MYH11"))
c("mCAF_LAMP5","mCAF_CD36","mCAF_STEAP4", "iCAF_CFD","mCAF_MYH11"))

levels(Fibroblasts@meta.data$celltype)
Fibroblasts$celltype = factor(Fibroblasts$celltype,
                              levels = c("emCAF_LAMP5","lpmCAF_CD36","iCAF_CFD","myoCAF_MYH11"))
nature_col = c("#c3a2b3", "#ccba88", "#7d95b1","#7fb4a4")
nature_col = c("#ad6955", "#60715f", "#b6976e","#899790")
nature_col = c("#D9534F", "#EDE574", "#98DF8A", "#0099CC","#9467BD")
nature_col = c("#f19138", "#b27da3", "#547eaa","#79b8b3")
cluster_order <- c("emCAF_LAMP5","lpmCAF_CD36","iCAF_CFD","myoCAF_MYH11")
scales::show_col(nature_col)
col_cluster <- setNames(c("#f19138", "#b27da3", "#547eaa","#79b8b3"),
                        c("emCAF_LAMP5","lpmCAF_CD36","iCAF_CFD","myoCAF_MYH11"))

sce <- Fibroblasts
df <- data.frame(sce@meta.data, 
                 sce@reductions$umap@cell.embeddings[,1:2])
#sce@reductions$tsne@cell.embeddings[,1:2])
colnames(df)[14:15] <- c("UMAP1","UMAP2")
#colnames(df)[10:11] <- c("tSNE1","tSNE2")
colnames(df)
ggplot(df,aes(x= UMAP1 , y = UMAP2 ,col=factor(celltype, levels = cluster_order))) + 
  geom_point(size = 1.5, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.8,'cm'),
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))

cell <- df %>%group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP1),
            UMAP2 = median(UMAP2))

rownames(cell) <- cell$celltype
A <- cell[cluster_order,]
A

a <- A
a$UMAP1 <-c( A$UMAP1[1]+1,A$UMAP1[2]-1.5,A$UMAP1[3]+1,A$UMAP1[4]+1)
a$UMAP2 <- c( A$UMAP2[1]+2.8,A$UMAP2[2]+4.55,A$UMAP2[3]-4.2,A$UMAP2[4]-3.7)
a

ggplot(df,aes(x= UMAP1 , y = UMAP2 , fill = factor(celltype, levels = cluster_order),color=factor(celltype, levels = cluster_order))) +
  stat_unchull(alpha = 0.5, size = 0.5, lty = 2) +
  geom_point(size = 2, shape=16)+
  scale_fill_manual("",values = col_cluster)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.6,'cm'),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_text_repel(data=a, aes(label=celltype), size=5,min.segment.length = Inf, point.padding = 0.3)+
  geom_segment(aes(x = min(UMAP1)-1.2 , y = min(UMAP2)-1.6,xend = min(UMAP1)+1.8, yend = min(UMAP2)-1.6),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP1)-1.2, y = min(UMAP2)-1.6,xend = min(UMAP1)-1.2,yend = min(UMAP2)+1.8),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), 
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP1)+0.3, y = min(df$UMAP2)-2, label = "UMAP1",
           color="black",size = 7) + 
  annotate("text", x = min(df$UMAP1)-1.6, y = min(df$UMAP2)-0.2, label = "UMAP2",
           color="black",size = 7,angle=90) 
#######################################################################################################
nature_col = c("#7fce7f", "#36638a", "#fd9734", "#cb3398","#b14949","#18A7B2")
cluster_order <- c("Tumor","Normal")
scales::show_col(nature_col)
col_cluster <- setNames(c("#36638a", "#7fce7f"),
                        #c("#dd7457", "#afbace"),
                        c("Tumor","Normal"))
df <- data.frame(sce@meta.data, 
                 sce@reductions$umap@cell.embeddings[,1:2])
#sce@reductions$tsne@cell.embeddings[,1:2])
colnames(df)[14:15] <- c("UMAP1","UMAP2")
#colnames(df)[10:11] <- c("tSNE1","tSNE2")
colnames(df)
table(df$spatialdistr)
ggplot(df,aes(x= UMAP1 , y = UMAP2 ,col=factor(spatialdistr, levels = cluster_order))) + 
  geom_point(size = 1.5, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(1.2,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))
######################################################################################
cluster_order <- c("CAYA","ADULT")
scales::show_col(nature_col)
col_cluster <- setNames(c("#888fba","#d7914e"),
                        c("CAYA","ADULT"))
col_cluster <- setNames(c('#C0000A','#50A5CC'),
                        c("CAYA","ADULT"))
df <- data.frame(sce@meta.data, 
                 sce@reductions$umap@cell.embeddings[,1:2])
#sce@reductions$tsne@cell.embeddings[,1:2])
colnames(df)[15:16] <- c("UMAP1","UMAP2")
#colnames(df)[10:11] <- c("tSNE1","tSNE2")
colnames(df)
table(df$group)
ggplot(df,aes(x= UMAP1 , y = UMAP2 ,col=factor(group, levels = cluster_order))) + 
  geom_point(size = 1.5, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(1.2,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))
###########################################################################################################
Dotplot_anno <- function(object,
                         single_type=T,
                         features=NULL,
                         bulk_feature=NULL,
                         bulk_samples=NULL,
                         group,
                         groups,
                         color,
                         order,
                         celltype_color,
                         level,
                         legend.postion=NULL,
                         heatmap=T){
  if (single_type==T){
    
    data.usage <- DotPlot(object,features = markers)$data
    data.anno <- data.frame(
      features.plot = unique(data.usage$features.plot),
      label = group)
    
    df.plot <- plyr::join(data.usage,data.anno)
    
    
    if(order==T){
      df.plot$id <- factor(df.plot$id,levels = sort(levels(df.plot$id),decreasing = T))
      df.plot$label <- factor(df.plot$label, levels = groups) 
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  as.numeric(id),
                              size = pct.exp, 
                              color = avg.exp.scaled))+
        geom_point() + 
        scale_size(range = c(0,6)) + 
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "Average\nexpression") +
        ylab("") + xlab("") +
        scale_y_continuous(breaks = 1:length(levels(df.plot$id)),
                           labels = levels(df.plot$id),
                           sec.axis = dup_axis())+ 
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(colour = 'black',linewidth = 0.5), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12,color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.position = legend.postion)
      
      
    }else{
      df.plot$id <- factor(df.plot$id,levels = rev(level))
      df.plot$label <- factor(df.plot$label, levels = groups)  
      
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  as.numeric(id),
                              size = pct.exp, 
                              color = avg.exp.scaled))+
        geom_point() + 
        scale_size(range = c(0,6)) + 
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "Average\nexpression") +
        ylab("") + xlab("") +
        scale_y_continuous(breaks = 1:length(levels(df.plot$id)),
                           labels = levels(df.plot$id),
                           sec.axis = dup_axis())+ 
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(colour = 'black',linewidth = 0.5), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12, color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.position = legend.postion)
      
    }
  }else{
    
    data.usage <- object
    data.anno <- data.frame(
      features.plot = unique(data.usage[, colnames(data.usage)%in% bulk_feature]),
      label = group)
    
    df.plot <- plyr::join(data.usage,data.anno)
    df.plot[, colnames(df.plot)%in% bulk_samples]<- as.factor(df.plot[, colnames(df.plot)%in% bulk_samples])
    
    if(order==T){
      df.plot[, colnames(df.plot)%in% bulk_samples] <- factor(df.plot[, colnames(df.plot)%in% bulk_samples],levels = level)
    }else{
      df.plot = df.plot
    }
    
    
    if('ratio' %in% colnames(df.plot)){
      p <- ggplot(df.plot,aes(x=features.plot,
                              y =  id,
                              size=ratio,
                              color = exp))+
        geom_point() + 
        scale_size(range = c(0,5))+
        scale_color_gradientn(colours = color,
                              guide = guide_colorbar(ticks.colour = "black",
                                                     frame.colour = "black"),
                              name = "-Log10(P)") +
        ylab("") + xlab("") +
        facet_grid(~label, scales="free_x",space = "free")+
        theme_classic() +
        theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
              axis.text.y = element_text(size=12, color="black"),
              axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
              
              axis.ticks.y = element_blank(),
              axis.text.y.right = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(colour = 'black',linewidth = 0.5), 
              
              panel.spacing=unit(0, "mm"), 
              strip.text.x = element_text(size=12, color = "black",
                                          vjust = 0.5,margin = margin(b = 3,t=3)),
              strip.background = element_rect(colour="black",size = 1),
              plot.margin=unit(c(1, 1, 1, 1),'cm'),
              legend.margin = margin(-0.2,-0.2,0,0,'cm'),
              legend.key.height = unit(0.4,'cm'),
              legend.key.width = unit(0.4,'cm'),
              legend.position = legend.postion)+
        guides(size=guide_legend(title="Ratio"))
    }else{
      
      if (heatmap==T){
        
        p <- ggplot(df.plot,aes(x=features.plot,
                                y =  id))+
          geom_tile(aes(fill = exp), colour = "black", size = 0.5) + 
          scale_fill_gradientn(colours = color,
                               guide = guide_colorbar(ticks.colour = "black",
                                                      frame.colour = "black"),
                               name = "Relative\nexpression") +
          ylab("") + xlab("") +
          facet_grid(~label, scales="free_x",space = "free")+
          theme_classic()+
          theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
                axis.text.y = element_text(size=12, color="black"),
                axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
                
                axis.ticks.y = element_blank(),
                axis.text.y.right = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line = element_blank(), 
                
                panel.spacing=unit(0, "mm"), 
                strip.text.x = element_text(size=12, color = "black",
                                            vjust = 0.5,margin = margin(b = 2,t=2)),
                strip.background = element_rect(colour="black",size = 0.5),
                plot.margin=unit(c(1, 1, 1, 1),'cm'),
                legend.position = legend.postion)+
          scale_x_discrete(expand = c(0,0))+
          scale_y_discrete(expand = c(0,0))
      }else{
        
        p <- ggplot(df.plot,aes(x=features.plot,
                                y =  id,
                                size=5,
                                color = exp))+
          geom_point() + 
          guides(size="none")+
          scale_color_gradientn(colours = color,
                                guide = guide_colorbar(ticks.colour = "black",
                                                       frame.colour = "black"),
                                name = "Relative\nexpression") +
          ylab("") + xlab("") +
          facet_grid(~label, scales="free_x",space = "free")+
          theme_classic() +
          theme(axis.text.x = element_text(size=12, angle=90, hjust=1, color="black"),
                axis.text.y = element_text(size=12, color="black"),
                axis.title.x = element_text(size=12,colour = 'black',hjust = 1),
                
                axis.ticks.y = element_blank(),
                axis.text.y.right = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line = element_line(colour = 'black',linewidth = 0.5), 
                
                panel.spacing=unit(0, "mm"), 
                strip.text.x = element_text(size=12, color = "black",
                                            vjust = 0.5,margin = margin(b = 3,t=3)),
                strip.background = element_rect(colour="black",size = 1),
                plot.margin=unit(c(1, 1, 1, 1),'cm'),
                legend.position = legend.postion)
      }
      
      
    }
    
    
  }
  
  
  g <- ggplot_gtable(ggplot_build(p))
  strips <- which(grepl('strip-', g$layout$name))
  fills <- celltype_color
  k = 1
  
  for (i in strips) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  plot(g)
  
}
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)
library(grid)
markers <- c( "LAMP5","CYP1B1","AGT","F2R","ENC1","FAP",#lamp5caf
              "CD36","STEAP4","FABP5","C20orf27","GJA4","ARHGDIB",#cd36mcaf
              "CFD","FBLN1","EFEMP1","DCN","C3",	"S100A10",#icaf
              "MYH11","ACTA2","SORBS2","NET1","CASQ2","WFDC1")#myh11caf

# DotPlot(human_data, features = markers, split.by = 'orig.ident')
library(dittoSeq)

human_data <- Fibroblasts

DotPlot(human_data, 
        features = markers, 
        split.by = 'orig.ident',
        cols = dittoColors())

#=====================================================================================
human_data@meta.data$newgroup <- paste0(human_data$celltype,"_",human_data$orig.ident)
Idents(human_data) <- 'newgroup'

unique(human_data$celltype)
Idents(human_data) <- factor(Idents(human_data), levels = c(paste0("emCAF_LAMP5","_",unique(human_data$orig.ident)),
                                                            paste0("lpmCAF_CD36","_",unique(human_data$orig.ident)),
                                                            paste0("iCAF_CFD","_",unique(human_data$orig.ident)),
                                                            paste0("myoCAF_MYH11","_",unique(human_data$orig.ident))))        

DotPlot(human_data, features = markers, cols = c("lightgrey", "red"))

table(human_data@meta.data$newgroup)

#=====================================================================================

myGrob <- grobTree(rectGrob(gp=gpar(fill="blanchedalmond", alpha=0.3)),
                   gTree(x0=0, x1=1, y0=0, y1=1, default.units="npc"))

p<- DotPlot(human_data, features = markers, cols = c("lightgrey", "red"))+
  geom_point(mapping = aes_string(size = 'pct.exp'),shape=21, color='black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size =15,face="bold"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.height = unit(0.4,'cm'),
        panel.grid = element_blank())+
  geom_vline(xintercept = 6.5, linetype=1, cex=0.5)+
  geom_vline(xintercept = 12.5, linetype=1, cex=0.5)+
  geom_vline(xintercept = 18.5, linetype=1, cex=0.5)+
  geom_hline(yintercept = 14.5, linetype=1, cex=0.5)+
  geom_hline(yintercept = 28.5, linetype=1, cex=0.5)+
  geom_hline(yintercept = 47.5, linetype=1, cex=0.5)+
  theme(plot.margin=unit(c(1, 1, 1, 1),'cm'))+
  annotation_custom(myGrob, xmin=-Inf,xmax=6.5,ymin=0,ymax=Inf)+
  annotation_custom(myGrob, xmin=12.5,xmax=18.5,ymin=0,ymax=Inf)
p
sample1 <- p[["data"]][["id"]]
levels(sample1)
group <- data.frame(sample=levels(sample1),
                    group=c(rep("ADULT",4),rep("CAYA",10),rep("ADULT",4),rep("CAYA",10),rep("ADULT",9),rep("CAYA",10),
                            rep("ADULT",4),rep("CAYA",10)),
                    celltype = c(rep("emCAF_LAMP5",14),rep("lpmCAF_CD36",14),
                                 rep("iCAF_CFD",19),rep("myoCAF_MYH11",14)),
                    tissue=c("T","P","T","T","T","P","T","P","T","P","T","P","T","P",
                             "T","P","T","P","T","P","T","P","T","P","T","P","T","P",
                             "T","P","P","P","T","P","T","P","T","T","P","T","P","T","P","T","P","T","P",
                             "T","P","T","T","T","P","T","P","T","P","T","P","T","P"))
group$y <- rownames(group)
group$y <- factor(group$y, levels = rownames(group))
cell_levels = c("emCAF_LAMP5","lpmCAF_CD36","iCAF_CFD","myoCAF_MYH11")
group$celltype <- factor(group$celltype, levels = cell_levels)

p1 <- ggplot(group, aes(x=1,y=y, fill=group))+
  geom_tile()+ 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values = c('#50A5CC','#C0000A'))
p1
p2 <- ggplot(group, aes(x=1,y=y, fill=celltype))+
  geom_tile()+ 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values = c("#f19138", "#b27da3", "#547eaa","#79b8b3"))
p2
p3 <- ggplot(group, aes(x=1,y=y, fill=tissue))+
  geom_tile()+ 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values = c("#7fce7f","#36638a"))
p3

library(patchwork)
plot =p2+p1+p+plot_layout(ncol = 3, widths  = c(0.1,0.1, 4),guides = 'collect')
plot =p2+p3+p1+p+plot_layout(ncol = 4, widths  = c(0.1,0.1,0.1, 4),guides = 'collect')
plot & theme(plot.margin = margin(0.5,0.5,0.5,0.5))
##########################################################################################################
KS_plot_density <- function(obj,
                            marker,
                            dim=c("TSNE","UMAP"),
                            size,
                            ncol=NULL
){
  require(ggplot2)
  require(ggrastr)
  require(Seurat)
  
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  if(dim=="TSNE"){
    
    xtitle = "tSNE1"
    ytitle = "tSNE2"
    
  }
  
  if(dim=="UMAP"){
    
    xtitle = "UMAP1"
    ytitle = "UMAP2"
  }
  
  
  if(length(marker)==1){
    
    plot <- FeaturePlot(obj, features = marker)
    data <- plot$data
    
    
    if(dim=="TSNE"){
      
      colnames(data)<- c("x","y","ident","gene")
      
    }
    
    if(dim=="UMAP"){
      
      colnames(data)<- c("x","y","ident","gene")
    }
    
    
    p <- ggplot(data, aes(x, y)) +
      geom_point_rast(shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = size) +
      geom_density_2d(data=data[data$gene>0,], 
                      aes(x=x, y=y), 
                      bins = 5, colour="black") +
      scale_fill_gradientn(colours = mypalette)+
      scale_colour_gradientn(colours = mypalette)+
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=30, face="bold.italic", hjust = 0.5),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    return(p)
    
  }else{
    
    gene_list <- list()
    
    
    
    for (i in 1:length(marker)) {
      plot <- FeaturePlot(obj, features = marker[i])
      data <- plot$data
      
      
      if(dim=="TSNE"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      if(dim=="UMAP"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      gene_list[[i]] <- data
      names(gene_list) <- marker[i]
    }
    
    plot_list <- list()
    
    
    for (i in 1:length(marker)) {
      
      p <- ggplot(gene_list[[i]], aes(x, y)) +
        geom_point_rast(shape = 21, stroke=0.25,
                        aes(colour=gene, 
                            fill=gene), size = size) +
        geom_density_2d(data=gene_list[[i]][gene_list[[i]]$gene>0,], 
                        aes(x=x, y=y), 
                        bins = 5, colour="black") +
        scale_fill_gradientn(colours = mypalette)+
        scale_colour_gradientn(colours = mypalette)+
        theme_bw()+ggtitle(marker[i])+
        labs(x=xtitle, y=ytitle)+
        theme(
          plot.title = element_text(size=30, face="bold.italic", hjust = 0.5),
          axis.text=element_text(size=8, colour = "black"),
          axis.title=element_text(size=12),
          legend.text = element_text(size =10),
          legend.title=element_blank(),
          aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      plot_list[[i]] <- p
    }
    
    
    Seurat::CombinePlots(plot_list, ncol = ncol)
    
    
  }
  
  
}


human_data <- Fibroblasts
KS_plot_density(obj=human_data, 
                marker=c("LAMP5","CD36","CFD","MYH11"), 
                dim = "UMAP", size =1.5, ncol =4)
ggsave("Plots/UAMP密度图(all).pdf", device = "pdf", width = 20, height = 10, dpi = 600)
KS_plot_density(obj=human_data, 
                marker=c("CD36"), 
                dim = "UMAP", size =2, ncol =4)
ggsave("Plots/UAMP密度图(CD36).pdf", device = "pdf", width = 10, height = 8, dpi = 600)
KS_plot_density(obj=human_data, 
                marker=c("CFD"), 
                dim = "UMAP", size =2, ncol =4)
ggsave("Plots/UAMP密度图(CFD).pdf", device = "pdf", width = 10, height = 8, dpi = 600)
KS_plot_density(obj=human_data, 
                marker=c("MYH11"), 
                dim = "UMAP", size =2, ncol =4)
ggsave("Plots/UAMP密度图(MYH11).pdf", device = "pdf", width = 10, height = 8, dpi = 600)
