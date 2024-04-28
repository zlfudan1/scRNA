library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(harmony)
library(devtools)
library(ggrepel)
library(ggridges)
library(RColorBrewer)
library(scales)
library(tidytree)
library(viridis)

setwd("./Thyrocytes")
load(file="./4celltypeharmonyThyrocytes-0.1.Rdata")
library(monocle)
table(scRNA$celltype)
sce=subset(scRNA,ident=c("Thyrocytes"))
table(sce@meta.data$orig.ident)
sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)
#
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo",
                                          "percent.mt", "nCount_RNA"), verbose = T)
sce <- RunPCA(sce, features = VariableFeatures(object = sce), npcs = 50)
# 
library(harmony)
sce <- RunHarmony(sce, "orig.ident")
names(sce@reductions)
sce <- RunUMAP(sce,  dims = 1:25, reduction = "harmony")

sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25) 

set.resolutions = c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1,1.2)
sce  <- FindClusters(object = sce , resolution = set.resolutions, verbose = FALSE) 
# 
for (res in c(0.05,0.1, 0.2, 0.3, 0.5, 0.8,1,1.2)) {
  sce=FindClusters(sce, resolution = res, algorithm = 1)
}
clustree(sce)
colnames(sce@meta.data)

sce <- FindNeighbors(sce, reduction = "harmony",dims = 1:25)
sce  <- FindClusters(object = sce , resolution = 0.1, verbose = FALSE) 
Thyrocytes <- RunUMAP(sce, reduction = "harmony", dims = 1:25)
Thyrocytes <- RunTSNE(Thyrocytes, reduction = "harmony", dims =1:25)#

DimPlot(Thyrocytes, reduction = "umap",label = T) 
DimPlot(Thyrocytes, reduction = "umap",label = T,split.by = "group")
DimPlot(Thyrocytes, reduction = "umap",label = T,split.by = "spatialdistr")

VlnPlot(Thyrocytes, 
        features = c('EPCAM','TG',"CD24","CLU"),
        pt.size = 0,
        ncol = 2)

scRNA_harmony <- Thyrocytes
diff.wilcox = FindAllMarkers(scRNA_harmony,logfc.threshold = 0.25,only.pos = T,test.use = "wilcox")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "4harmonyThyrocytes-0.1_diff_genes_wilcox.csv", row.names = F)
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
write.csv(ClusterMarkerTOP25,'4harmonyThyrocytes-0.1_TopMarkersTOP25.csv', row.names=F)

M <- Thyrocytes
M@meta.data <- M@meta.data[,-(11:15)]
M@meta.data <- M@meta.data[,-(12)]
M@meta.data <- M@meta.data[,-(14:15)]
Thyrocytes@meta.data <- M@meta.data

Thyrocytes2=Thyrocytes
Thyrocytes2$group2<-Thyrocytes2$orig.ident
table(Thyrocytes2$orig.ident)
Thyrocytes2$group2<-as.character(Thyrocytes2$group2)
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA1P'] <- 'CATA-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA1T'] <- 'CAYA-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA2P'] <- "CATA-P"
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA2T'] <- 'CAYA-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA3P'] <- 'CATA-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA3T'] <- 'CAYA-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA4P'] <- 'CATA-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA4T'] <- 'CAYA-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA5P'] <- 'CATA-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'CAYA5T'] <- 'CAYA-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT1P'] <- 'ADUL-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT1T'] <- 'ADUL-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT2T'] <- 'ADUL-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT2P'] <- 'ADUL-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT3T'] <- 'ADUL-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT3P'] <- 'ADUL-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT4T'] <- 'ADUL-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT4P'] <- 'ADUL-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT5T'] <- 'ADUL-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT5P'] <- 'ADUL-P'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT6T'] <- 'ADUL-T'
Thyrocytes2$group2[Thyrocytes2$group2 == 'ADULT6P'] <- 'ADUL-P'
Thyrocytes=Thyrocytes2
#-------------------------------------------------------------------------------
save(Thyrocytes,all.markers,file="4celltypeharmonyThyrocytes-0.1.Rdata") 
load(file="/home/user-y3/subcelltype/Thyrocytes/Thyrocytes-Score.Rdata")
DimPlot(Thyrocytes, reduction = "umap",label = T, split.by = 'group') 

#
df_Thyrocytes <- Thyrocytes@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = Thyrocytes@meta.data$celltype)%>%
  cbind(seurat_clusters = Thyrocytes@meta.data$seurat_clusters)
colnames(df_Thyrocytes)
#
cluster_order <- c('0','1',"2","3","4","5")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', '#D8A326'),
                        c('0','1',"2","3","4","5"))
#
cell <- df_Thyrocytes %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:5)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_Thyrocytes)
p <- ggplot(df_Thyrocytes,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
  annotate("text", x = min(df_Thyrocytes$UMAP_1) +1, y = min(df_Thyrocytes$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 4) + 
  annotate("text", x = min(df_Thyrocytes$UMAP_1) -0.5, y = min(df_Thyrocytes$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 4,angle=90)+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("Thyrocytes_UMAP.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p)

#
B <- A
B$x <- 1
B$lab <- c(5:0)
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
ggsave("Thyrocytes_umap_legend.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, leg)
###-----------------------------------------------------------------------------
#
Thyrocytes_CAYA = Thyrocytes[,Thyrocytes@meta.data$group%in%c('CAYA')]
df_Thyrocytes_CAYA <- Thyrocytes_CAYA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(Thyrocytes_CAYA = Thyrocytes_CAYA@meta.data$celltype)%>%
  cbind(seurat_clusters = Thyrocytes_CAYA@meta.data$seurat_clusters)
colnames(df_Thyrocytes_CAYA)
#
cluster_order <- c('0','1',"2","3","4","5")
#
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', '#D8A326'),
                        c('0','1',"2","3","4","5"))
#
cell_CAYA <- df_Thyrocytes_CAYA %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell_CAYA) <-cell_CAYA$seurat_clusters
A <- cell_CAYA[cluster_order,]
a <- c(0:5)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_Thyrocytes_CAYA)
p <- ggplot(df_Thyrocytes_CAYA,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 2, fill = NA))+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("Thyrocytes_CAYA_UMAP.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p)
###-----------------------------------------------------------------------------
#
Thyrocytes_ADULT = Thyrocytes[,Thyrocytes@meta.data$group%in%c('ADULT')]
df_Thyrocytes_ADULT <- Thyrocytes_ADULT@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(Thyrocytes_ADULT = Thyrocytes_ADULT@meta.data$celltype)%>%
  cbind(seurat_clusters = Thyrocytes_ADULT@meta.data$seurat_clusters)
colnames(df_Thyrocytes_ADULT)

cluster_order <- c('0','1',"2","3","4","5")

col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', '#D8A326'),
                        c('0','1',"2","3","4","5"))

cell_ADULT <- df_Thyrocytes_ADULT %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell_ADULT) <-cell_ADULT$seurat_clusters
A <- cell_ADULT[cluster_order,]
a <- c(0:5)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df_Thyrocytes_ADULT)
p <- ggplot(df_Thyrocytes_ADULT,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", size = 2, fill = NA))+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3) 

ggsave("Thyrocytes_ADULT_UMAP.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p)
#===============================================================================
#monocle2
library(Seurat)
install.packages("./monocle_2.26.0.tar.gz", repos = NULL, type = "source")
library(monocle)
devtools::load_all("./monocle")
library(dplyr)
#===============================================================================
seurat_to_monocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- GetAssayData(otherCDS, assay = assay, slot = slot)
    data <- data[rowSums(as.matrix(data)) != 0,]
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      } else {
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    } 
  } 
  return(monocle_cds)
}
thyrocyte_data <- Thyrocytes
thyrocyte_monocle_fun <- seurat_to_monocle(thyrocyte_data, assay = "RNA", slot = "counts")

#===============================================================================

#
thyrocyte_monocle <- estimateSizeFactors(thyrocyte_monocle_fun)
thyrocyte_monocle <- estimateDispersions(thyrocyte_monocle)

#
thyrocyte_monocle <- detectGenes(thyrocyte_monocle, min_expr = 0.1)
print(head(fData(thyrocyte_monocle)))
print(head(pData(thyrocyte_monocle)))
expressed_genes <- row.names(subset(fData(thyrocyte_monocle), num_cells_expressed >= 10))
pData(thyrocyte_monocle)$Total_mRNAs <- Matrix::colSums(exprs(thyrocyte_monocle))
thyrocyte_monocle <- thyrocyte_monocle[,pData(thyrocyte_monocle)$Total_mRNAs < 1e6]

#differentialGeneTest
cds_DGT <- thyrocyte_monocle
diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
plot_ordering_genes(cds_DGT)
save(cds_DGT, file = "cds_DGT_Thyrocyte.RData")
load('./2output_of_cds.Rdata')
cds_DGT <- cds
#################################################################################
cds_DGT <- reduceDimension(cds_DGT, max_components = 2,reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT, root_state = NULL, num_paths = NULL, reverse = T) 
save(cds_DGT, file = "cds_orderCells_Thyrocyte.RData")
#-------------------------------------------------------------------------------
#
plot_cell_trajectory(cds_DGT, cell_size = 1.0, color_by = "orig.ident") +
  facet_wrap(~group, nrow = 1)
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime") +facet_wrap(~group, nrow = 1)
plot_cell_trajectory(cds_DGT, color_by = "orig.ident",cell_link_size = 1.5)
plot_cell_trajectory(cds_DGT, color_by = "State")
plot_cell_trajectory(cds_DGT, color_by = 'group')
plot_cell_trajectory(cds_DGT, color_by = 'Cluster') +
  facet_wrap(~group, nrow = 1)
#-------------------------------------------------------------------------------
data_df <- t(reducedDimS(cds_DGT)) %>% as.data.frame() %>%
  select_(Component_1 = 1, Component_2 = 2) %>% 
  rownames_to_column("cells") %>% 
  mutate(pData(cds_DGT)$State) %>% 
  mutate(pData(cds_DGT)$Pseudotime, 
         pData(cds_DGT)$orig.ident, 
         pData(cds_DGT)$celltype,
         pData(cds_DGT)$seurat_clusters,
         pData(cds_DGT)$group)

colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","celltype",'seurat_clusters','group')

View(data_df)
data_df <- as.data.frame(data_df)
data_df$newState[data_df$State=="1"] = c("State3")
data_df$newState[data_df$State=="2"] = c("State1")
data_df$newState[data_df$State=="3"] = c("State2")

data_df$State <- data_df$newState
data_df <- data_df[,-10]
#-------------------------------------------------------------------------------
#
dp_mst <- minSpanningTree(cds_DGT)
reduced_dim_coords <- reducedDimK(cds_DGT)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))

edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% select_(source = "sample_name", 
                                     source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
  left_join(ica_space_df %>% select_(target = "sample_name", 
                                     target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")
#===============================================================================
#
Cellratio <- prop.table(table(data_df$State, data_df$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"orig.ident","Freq")
#==============================================================================
#
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)
g <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color =Pseudotime)) + 
  scale_color_viridis()+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_arc(arrow = arrow(length = unit(0.15, "inches"), 
                         type = "closed",angle=30),
           aes(x0=0,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)+
  geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",
               aes(x0=-15,y0=0,r0=0,r=2.5,amount=Freq,fill=orig.ident))+
  geom_arc_bar(data=subset(Cellratio,State=='2'),stat = "pie",
               aes(x0=2,y0=9,r0=0,r=2.5,amount=Freq,fill=orig.ident))+
  geom_arc_bar(data=subset(Cellratio,State=='3'),stat = "pie",
               aes(x0=10,y0=-8,r0=0,r=2.5,amount=Freq,fill=orig.ident))
ggsave("Thyrocytes_monocle2.pdf",device = cairo_pdf, width = 11, height = 7, dpi = 300, g)
#-------------------------------------------------------------------------
col_cluster <- setNames(c('#9e0142', '#5e4fa2', '#ffd700'),
                        c('State1', 'State2', 'State3'))
g <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color = State)) + 
  scale_color_manual("",values = col_cluster)+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_arc(arrow = arrow(length = unit(0.15, "inches"), 
                         type = "closed",angle=30),
           aes(x0=-2,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)

ggsave("Thyrocytes_State_monocle2.pdf",device = cairo_pdf, width = 11, height = 7, dpi = 300, g)
#-------------------------------------------------------------------------
#
data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% 
  select_(Component_1 = 1, Component_2 = 2) %>% 
  rownames_to_column("cells") %>%
  mutate(pData(cds)$State) %>% 
  mutate(pData(cds)$Pseudotime, 
         pData(cds)$orig.ident, 
         pData(cds)$celltype,
         pData(cds)$seurat_clusters,
         pData(cds)$group)
cds_DGT_CAYA = subset(data_df,data_df$`pData(cds)$group`=='CAYA')
data_df_CAYA <- cds_DGT_CAYA

colnames(data_df_CAYA) <- c("cells","Component_1","Component_2","State",
                            "Pseudotime","orig.ident","celltype",'seurat_clusters','group')

View(data_df_CAYA)
data_df_CAYA <- as.data.frame(data_df_CAYA)
data_df_CAYA$newState[data_df_CAYA$State=="1"] = c("State3")
data_df_CAYA$newState[data_df_CAYA$State=="2"] = c("State1")
data_df_CAYA$newState[data_df_CAYA$State=="3"] = c("State2")

data_df_CAYA$State <- data_df_CAYA$newState
data_df_CAYA <- data_df_CAYA[,-10]

dp_mst <- minSpanningTree(cds_DGT)
reduced_dim_coords <- reducedDimK(cds_DGT)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))

edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% select_(source = "sample_name", 
                                     source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
  left_join(ica_space_df %>% select_(target = "sample_name", 
                                     target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

Cellratio <- prop.table(table(data_df_CAYA$State, data_df_CAYA$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"orig.ident","Freq")

g1 <- ggplot() + 
  geom_point_rast(data = data_df_CAYA, aes(x = Component_1, 
                                           y = Component_2,
                                           color =State)) + 
  scale_color_manual("",values = col_cluster)+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("Thyrocytes_monocle2_CAYA.pdf",device = cairo_pdf, width = 11, height = 7, dpi = 300, g1)
#-------------------------------------------------------------------------------
#
data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% 
  select_(Component_1 = 1, Component_2 = 2) %>% #
  rownames_to_column("cells") %>% 
  mutate(pData(cds)$State) %>%
  mutate(pData(cds)$Pseudotime, 
         pData(cds)$orig.ident, 
         pData(cds)$celltype,
         pData(cds)$seurat_clusters,
         pData(cds)$group)
cds_DGT_ADULT = subset(data_df,data_df$`pData(cds)$group`=='ADULT')
data_df_ADULT <- cds_DGT_ADULT

colnames(data_df_ADULT) <- c("cells","Component_1","Component_2","State",
                             "Pseudotime","orig.ident","celltype",'seurat_clusters','group')

View(data_df_ADULT)
data_df_ADULT <- as.data.frame(data_df_ADULT)
data_df_ADULT$newState[data_df_ADULT$State=="1"] = c("State3")
data_df_ADULT$newState[data_df_ADULT$State=="2"] = c("State1")
data_df_ADULT$newState[data_df_ADULT$State=="3"] = c("State2")

data_df_ADULT$State <- data_df_ADULT$newState
data_df_ADULT <- data_df_ADULT[,-10]

dp_mst <- minSpanningTree(cds_DGT)
reduced_dim_coords <- reducedDimK(cds_DGT)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))

edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% select_(source = "sample_name", 
                                     source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
  left_join(ica_space_df %>% select_(target = "sample_name", 
                                     target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

Cellratio <- prop.table(table(data_df_ADULT$State, data_df_ADULT$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"orig.ident","Freq")

g2 <- ggplot() + 
  geom_point_rast(data = data_df_ADULT, aes(x = Component_1, 
                                            y = Component_2,
                                            color =State)) + 
  scale_color_manual("",values = col_cluster)+#密度色
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("Thyrocytes_monocle2_ADULT.pdf",device = cairo_pdf, width = 11, height = 7, dpi = 300, g2)
#===============================================================================
load("./newThyrocytes-Score.Rdata")
setwd("./thyrocytes/TDS")
sce <- Thyrocytes
colnames(sce@meta.data)
sce@meta.data <- sce@meta.data[,-c(19:20)]

#
BRAF_gene<-readxl::read_xlsx("BRAF1.xlsx")
#
BRAF_gene <- BRAF_gene[1:45,]
gene<-as.list(BRAF_gene)
sce<-AddModuleScore(object=sce,features = gene,ctrl=100,name="BRAF")
colnames(sce@meta.data)
colnames(sce@meta.data)[20]<-"BRAF"

#
RAS_gene<-readxl::read_xlsx("RAS.xlsx")
#
gene<-as.list(RAS_gene)
sce<-AddModuleScore(object=sce,features = gene,ctrl=100,name="RAS")
colnames(sce@meta.data)
colnames(sce@meta.data)[21]<-"RAS"

Thyrocytes <- sce

setwd("./Thyrocytes")
save(Thyrocytes,file="newThyrocytes-Score.Rdata")
load('newThyrocytes-Score.Rdata')
#===============================================================================
#
TDS_gene<-readxl::read_xlsx("TDS.xlsx")
View(TDS_gene)
gene<-as.list(TDS_gene)
sce<-AddModuleScore(object=cds_DGT,features = gene,ctrl=100,name="TDS")
colnames(sce@meta.data)

DimPlot(Thyrocytes, reduction = "umap",label = T,split.by = "group")
VlnPlot(Thyrocytes,features = "TDS",pt.size = 0,group.by = "group")
VlnPlot(Thyrocytes,features = "TDS",pt.size = 0,group.by = "spatialdistr")
VlnPlot(Thyrocytes,features = "TDS",pt.size = 0,group.by = "seurat_clusters")
VlnPlot(Thyrocytes,features = "TDS",pt.size = 0,group.by = "State")
#
pTDS <- ggboxplot(Thyrocytes@meta.data, x="State", y="TDS", width = 0.6, 
                  fill="State",
                  xlab = F, 
                  bxp.errorbar=T,
                  bxp.errorbar.width=0.5, 
                  size=1, 
                  outlier.shape=NA, 
                  legend = "right")

library(rstatix)
Thyrocytes$State = factor(Thyrocytes$State,
                          levels = c("State1","State2","State3"))
datastate <- as.data.frame(Thyrocytes@meta.data)
colnames(datastate)
compare_means(TDS ~ State,  
              data = datastate, 
              method = "wilcox.test")
comparisons <- list( c("State1", "State2"),
                     c("State2", "State3"),
                     c("State1", "State3") )
col_cluster <- setNames(c( '#5e4fa2', '#ffd700','#9e0142'),
                        c('State1', 'State2', 'State3'))

pTDS <- ggplot(data=datastate,aes(x=State,y=TDS,fill=State))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(1.8,1.4,2.2),size=6)+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="State",y="TDS score")
ggsave("Thyrocytes_TDSscore.pdf",device = cairo_pdf, width = 9, height = 7, dpi = 300, pTDS)

pFUSION <- ggplot(data=datastate,aes(x=State,y=FUSION,fill=State))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.3,0.4,0.5),size=6)+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="State",y="FUSION score")
ggsave("Thyrocytes_FUSIONscore.pdf",device = cairo_pdf, width = 9, height = 7, dpi = 300, pFUSION)

pBRAF <- ggplot(data=datastate,aes(x=State,y=BRAF,fill=State))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.6,0.7,0.8),size=6,method = "wilcox.test")+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="State",y="BRAF score")
ggsave("Thyrocytes_BRAFscore.pdf",device = cairo_pdf, width = 9, height = 7, dpi = 300, pBRAF)

pRAS <- ggplot(data=datastate,aes(x=State,y=RAS,fill=State))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.7,0.6,0.8),size=6,method = "t.test")+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="State",y="RAS score")
ggsave("Thyrocytes_RASscore.pdf",device = cairo_pdf, width = 9, height = 7, dpi = 300, pRAS)

#-------------------------------------------------------------------------------
Thyrocytes$group = factor(Thyrocytes$group,
                          levels = c("CAYA","ADULT"))
datagroup <- as.data.frame(Thyrocytes@meta.data)
colnames(datagroup)
compare_means(TDS ~ group,  
              data = datagroup, 
              method = "wilcox.test")
comparisons <- list( c("ADULT","CAYA") )
col_cluster <- setNames(c('#C0000A','#50A5CC'),
                        c("CAYA","ADULT"))

ggplot(data=datagroup,aes(x=group,y=TDS,fill=group))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(1.8,1.4,2.2),size=6)+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="group",y="TDS score")

ggplot(data=datagroup,aes(x=group,y=FUSION,fill=group))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.3,0.4,0.5),size=6)+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="group",y="FUSION score")

ggplot(data=datagroup,aes(x=group,y=BRAF,fill=group))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.6,0.7,0.8),size=6,method = "wilcox.test")+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="group",y="BRAF score")

ggplot(data=datagroup,aes(x=group,y=RAS,fill=group))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.7,0.6,0.8),size=6,method = "t.test")+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="group",y="RAS score")
#-------------------------------------------------------------------------------
Idents(Thyrocytes)="State"
sce=subset(Thyrocytes,ident=c("State1"))
sce$group = factor(sce$group,
                   levels = c("CAYA","ADULT"))
datagroup <- as.data.frame(sce@meta.data)
colnames(datagroup)
compare_means(FUSION ~ group,  
              data = datagroup, 
              method = "wilcox.test")
comparisons <- list( c("ADULT","CAYA") )
col_cluster <- setNames(c('#C0000A','#50A5CC'),
                        c("CAYA","ADULT"))

ggplot(data=datagroup,aes(x=group,y=FUSION,fill=group))+
  geom_violin(scale = "width",color=NA,alpha=0.7)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_boxplot(width=0.3,outlier.shape = NA,position=position_dodge(0.9))+
  stat_compare_means(comparisons = comparisons,aes(label = ..p.signif..),label.y =c(0.3,0.4,0.5),size=6)+
  scale_fill_manual("",values = col_cluster)+
  theme_classic()+ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.15, "cm"),
        #axis.text.x.bottom = element_text(size = 14,color = "black",angle = 45,hjust = 1),
        axis.text.x.bottom = element_text(),
        axis.text.y.left = element_text(size = 14,color = "black"),
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_text(size = 16)) +
  labs(x="group",y="FUSION score")
#-------------------------------------------------------------------------------
#
library(ggplot2)
mydata<- FetchData(Thyrocytes,vars = c("UMAP_1","UMAP_2","TDS"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2, color = TDS))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

#-------------------------------------------------------------------------------
#
colnames(sce@meta.data)
DimPlot(Thyrocytes, reduction = "umap",label = T,split.by = "group")
VlnPlot(Thyrocytes,features = "BRAF",pt.size = 0,group.by = "group")
VlnPlot(Thyrocytes,features = "BRAF",pt.size = 0,group.by = "spatialdistr")
VlnPlot(Thyrocytes,features = "BRAF",pt.size = 0,group.by = "seurat_clusters")
VlnPlot(Thyrocytes,features = "BRAF",pt.size = 0,group.by = "State")

pBRAF <- ggboxplot(Thyrocytes@meta.data, x="State", y="BRAF", width = 0.6, 
                   fill="State",
                   xlab = F, 
                   bxp.errorbar=T,
                   bxp.errorbar.width=0.5, 
                   size=1, 
                   outlier.shape=NA, 
                   legend = "right")

library(ggplot2)
mydata<- FetchData(Thyrocytes,vars = c("UMAP_1","UMAP_2","BRAF"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2, color = BRAF))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
#-------------------------------------------------------------------------------
#
colnames(sce@meta.data)
DimPlot(Thyrocytes, reduction = "umap",label = T,split.by = "group")
VlnPlot(Thyrocytes,features = "RAS",pt.size = 0,group.by = "group")
VlnPlot(Thyrocytes,features = "RAS",pt.size = 0,group.by = "spatialdistr")
VlnPlot(Thyrocytes,features = "RAS",pt.size = 0,group.by = "seurat_clusters")
VlnPlot(Thyrocytes,features = "RAS",pt.size = 0,group.by = "State")

pRAS <- ggboxplot(Thyrocytes@meta.data, x="State", y="RAS", width = 0.6, 
                  fill="State",
                  xlab = F, 
                  bxp.errorbar=T,
                  bxp.errorbar.width=0.5, 
                  size=1, 
                  outlier.shape=NA, 
                  legend = "right")

mydata<- FetchData(Thyrocytes,vars = c("UMAP_1","UMAP_2","RAS"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2, color = RAS))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
#-------------------------------------------------------------------------------
#
FUSION_gene<-readxl::read_xlsx("FUSION.xlsx")
View(FUSION_gene)
gene<-as.list(FUSION_gene)
sce<-AddModuleScore(object=cds_DGT,features = FUSION,ctrl=100,name="FUSION")
colnames(sce@meta.data)
DimPlot(Thyrocytes, reduction = "umap",label = T,split.by = "group")
VlnPlot(Thyrocytes,features = "FUSION",pt.size = 0,group.by = "group")
VlnPlot(Thyrocytes,features = "FUSION",pt.size = 0,group.by = "spatialdistr")
VlnPlot(Thyrocytes,features = "FUSION",pt.size = 0,group.by = "seurat_clusters")
VlnPlot(Thyrocytes,features = "FUSION",pt.size = 0,group.by = "State")

pFUSION <- ggboxplot(Thyrocytes@meta.data, x="State", y="FUSION", width = 0.6, 
                     fill="State",
                     xlab = F,
                     bxp.errorbar=T,
                     bxp.errorbar.width=0.5,
                     size=1,
                     outlier.shape=NA, 
                     legend = "right")

mydata<- FetchData(Thyrocytes,vars = c("UMAP_1","UMAP_2","FUSION"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2, color = FUSION))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
library(gridExtra)
p <- grid.arrange(pTDS, pBRAF, pRAS, pFUSION, nrow = 2)
p
ggsave("TDS_BRAF_RAS_FUSION.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p)

#-------------------------------------------------------------------------------
load(file="./harmonythyrocytes-0.1.Rdata")
#
cds_DGT_pseudotimegenes <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~sm.ns(Pseudotime)")

cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)
marker <- FindAllMarkers(Thyrocytes, only.pos = T,logfc.threshold = 0.5)
marker <- marker[which(marker$p_val_adj<0.05),]
top15 <- marker %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top15_ordergene <- cds_DGT_pseudotimegenes_sig[top15$gene, ]
Time_genes <- top15_ordergene %>% pull(gene_short_name) %>% as.character()
Time_genes <- unique(Time_genes)
write.csv(Time_genes, file = 'Time_genes.csv')
write.csv(cds_DGT@featureData@data$gene_short_name, file = 'cds_DGT@featureData@data$gene_short_name.csv')

Time_genes <- c('C2orf40',	'H3F3A',	'HIST1H4C',	'H3F3B',	'NDUFC2',	'NGFRAP1',	
                'ATP5L',	'MTRNR2L12',	'PLCG2',	'IGFBP5',	'SELM',	'ATP5E',	
                'APOC1',	'APOE',	'ARMCX3',	'CD74',	'CXCL14',	'RACK1',	'HLA-DPA1',	
                'H3-3A',	'S100A2',	'HLA-DRA',	'ECM1',	'CCN1',	'GDF15',	'HLA-DRB5',	
                'CXCL2',	'H3-3B',	'SELENOM',	'CCN2',	'FUS',	'NFKBIZ',	'CCNL1',	
                'SLC38A2',	'PPP1R10',	'GLS',	'VMP1',	'IFRD1',	'SQSTM1',	'INTS6',	
                'HSPH1',	'NEAT1',	'HMGA2',	'LRRK2',	'XIST',	'NUPR1',	'PPDPF',	
                'GABARAP',	'ATP5F1E',	'ATP5MC2',	'SNHG29',	'TFF3',	'MT1F',	'BEX3',	
                'MT1G',	'SELENOW',	'ATP5MG',	'ATP5PF',	'CXCR4',	'SRGN',	'CD69',
                'CD52',	'IL7R',	'IL32',	'BTG1',	'CD3E',	'IGHG4',	'CCL4',	'IGKC',	'RPS29',
                'MT2A',	'RPS17',	'MT1X',	'SPARCL1',	'BGN',	'COL3A1',	'IFI27',	'IGFBP7',
                'IFITM1',	'TAGLN',	'RGS5',	'ACTA2',	'SPARC',	'COL1A2',	'RGS5.1',	'COL1A1')
plot_pseudotime_heatmap(cds_DGT[cds_DGT@featureData@data$gene_short_name %in% Time_genes,], 
                        num_cluster = 4, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

plot_pseudotime_heatmap(cds_DGT[cds_DGT@featureData@data$gene_short_name %in% Time_genes,], 
                        num_cluster = 4, 
                        show_rownames = T, 
                        return_heatmap = T,
                        hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))

plot_pseudotime_heatmap(cds_DGT[cds_DGT@featureData@data$gene_short_name %in% Time_genes,], 
                        num_cluster = 4, 
                        show_rownames = T, 
                        return_heatmap = F,
                        hmcols = viridis(256))
#-------------------------------------------------------------------------------
#
cds_subset <- cds_DGT
newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                       max(pData(cds_subset)$Pseudotime),
                                       length.out = 100)) 
m <- genSmoothCurves(cds_DGT[Time_genes,], 
                     trend_formula = '~sm.ns(Pseudotime, df=3)',  
                     relative_expr = T, new_data = newdata)
m=m[!apply(m,1,sum)==0,]
m <- log10(m+1) 
#
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m[m > 3] = 3
m[m <- 3] = -3

#
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

#
p1 <- pheatmap(m, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_method = "ward.D2",
               cutree_rows=4,
               filename=NA,
               border_color = NA,
               fontsize_row = 8,
               color=colorRampPalette(c("navy","white","firebrick3"))(100),
               clustering_callback = callback)


#
annotation_col = data.frame(
  pseudotime = rescale(newdata$Pseudotime, to = c(-1, 1)))
row.names(annotation_col) <- colnames(m)

#
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))
row.names(annotation_row) <- rownames(m)

rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
names(rowcolor) <- c("1","2","3","4") 

#
ann_colors <- list(pseudotime=viridis(100),
                   Cluster=rowcolor) 

p2 <- pheatmap(m, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_method = "ward.D2",
               cutree_rows=4,
               filename=NA,
               border_color = NA,
               fontsize_row = 8,
               color=viridis(100),
               annotation_col = annotation_col,
               annotation_colors=ann_colors,
               annotation_row = annotation_row,
               clustering_callback = callback,
               annotation_names_col = F,
               annotation_names_row = F,
               main="Pseudotime")
#-------------------------------------------------------------------------------
library(reshape)
heat_gg <- m
heat_gg <- as.data.frame(heat_gg)
heat_gg <- heat_gg%>% mutate(gene=row.names(.)) %>% melt()

p3 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white",
                       name='Pseudotime')+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',
                                   size = 8))+
  scale_y_discrete(position = "right")

#
library(ggtree)
library(tidytree)
d <- m%>% as.data.frame()
hc <- hclust(dist(d), method = "ward.D2")
clus <-cutree(hc, 4)
d1 = data.frame(label=names(clus),member=factor(clus))
ph <- ggtree(as.phylo(hc))

#
cluster <- d1 %>%
  ggplot(aes(x=1, y=label,fill=member))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  labs(fill = "cluster")+
  scale_fill_manual(values = c("#708090",'#68A180','#F3B1A0', '#D6E7A3'))


#
group <- colnames(m) %>% as.data.frame() %>% 
  mutate(group=newdata$Pseudotime) %>%
  mutate(p="group") %>%
  ggplot(aes(.,y=1,fill=group))+
  geom_tile() + 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_gradientn(colours = c("#85B22E","#E29827",'#57C3F3',"#922927"))+
  scale_x_discrete(expand = c(0,0))+
  geom_segment(aes(x = 5, y = 1, xend = 95, yend = 1),
               arrow = arrow(length = unit(0.1, "inches"),
                             type = 'closed'))+
  theme(plot.margin = margin(0,0,0,0))

library(aplot)
p4 <- p3 %>%insert_left(cluster, width = 0.04)%>%
  insert_top(group, height =0.02)%>%
  insert_left(ph,width=0.1)

p4
ggsave("heatmap_monocle2.pdf",device = cairo_pdf, width = 8, height = 11, dpi = 300,p4)

#-------------------------------------------------------------------------------
#
library(ggridges)
plotdf=pData(cds_DGT)
plotdf$celltype <- factor(plotdf$State, 
                          levels = c('1', '2', '3'))

p5 <- ggplot(plotdf, aes(x=Pseudotime,y=State,fill=State))+
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8))+
  scale_x_continuous(position = 'top')
ggsave("rage_monocle2.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300,p5)

p6 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white",
                       name='Pseudotime')+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',
                                   size = 8))+
  scale_y_discrete(position = "right")

p7 <- p6 %>% insert_left(cluster, width = 0.04)%>%
  insert_top(group, height =0.02)%>%
  insert_left(ph,width=0.1)

library(patchwork)
p5+p7+plot_layout(ncol = 1,heights=c(1,5),guides = 'collect')
p5 | p7

#-------------------------------------------------------------------------------
genes <- c('APOC1',	'APOE',	'ARMCX3',	'CD74',	'CXCL14',	'RACK1')
genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(cds_DGT)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#
pData(cds_DGT) = cbind(pData(cds_DGT), gene_exp)

#
data <- pData(cds_DGT)
colnames(data)
#
data<-data %>% select("celltype","group","orig.ident","seurat_clusters","Pseudotime",
                      'APOC1',	'APOE',	'ARMCX3',	'CD74',	'CXCL14',	'RACK1')

features <- c('APOC1',	'APOE',	'ARMCX3',	'CD74',	'CXCL14',	'RACK1')
plist <- list()
#ggplot
for (i in 1:length(features)){
  df <- data[, colnames(data)%in% c("orig.ident","Pseudotime",features[i])]
  colnames(df) <- c("orig.ident","Pseudotime",'gene')
  p <- ggplot(df, aes(x=Pseudotime, 
                      y=gene)) + 
    geom_point(aes(color=orig.ident), size=0.8)+
    geom_smooth(method = "loess",level = 0.95,
                formula = y~x, color='black',se=F)+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          axis.title = element_blank(),
          axis.text  = element_text(size = 15, colour = 'black'),
          legend.position = 'none')+
    ggtitle(features[i])
  plist[[i]] <- p
}

CombinePlots(plist, ncol = 2)

#---------------------------------------------------------------------------

cds_DGT_CAYA = cds_DGT[,cds_DGT@phenoData@data$group %in% c('CAYA')]
plot_cell_trajectory(cds_DGT_CAYA, cell_size = 2.2, color_by = "group") +
  facet_wrap(~celltype, nrow = 2)
plot_cell_trajectory(cds_DGT_CAYA, color_by = "Pseudotime")
plot_cell_trajectory(cds_DGT_CAYA, color_by = "orig.ident",cell_link_size = 1.5)
plot_cell_trajectory(cds_DGT_CAYA, color_by = "State")
plot_cell_trajectory(cds_DGT_CAYA, color_by = 'group')
#ADULT
cds_DGT_ADULT = cds_DGT[,cds_DGT@phenoData@data$group %in% c('ADULT')]
plot_cell_trajectory(cds_DGT_ADULT, cell_size = 2.2, color_by = "group") +
  facet_wrap(~celltype, nrow = 2)
plot_cell_trajectory(cds_DGT_ADULT, color_by = "Pseudotime")
plot_cell_trajectory(cds_DGT_ADULT, color_by = "orig.ident",cell_link_size = 1.5)
plot_cell_trajectory(cds_DGT_ADULT, color_by = "State")
plot_cell_trajectory(cds_DGT_ADULT, color_by = 'group')
#基因表达拟时图,这种方法是5中基因单独展示
markers <- c('EZH2', 'TSHR', 'IL17', 'IL17R', 'CXCR6')
plot_cell_trajectory(cds_DGT, markers = markers, use_color_gradient=T,cell_size = 1,cell_link_size = 1.0)

#各个分组拟时图
plot_cell_trajectory(cds_DGT, color_by = "orig.ident") + facet_wrap(~group, nrow = 1)

#################################################################################
#
cds_DGT_pseudotimegenes <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~sm.ns(Pseudotime)")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)
write.csv(cds_DGT_pseudotimegenes_sig,"cds_DGT_pseudotimegenes_sig.csv")

#
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 1e-320)
Time_genes <- cds_DGT_pseudotimegenes_sig %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds_DGT[Time_genes,], 
                             num_cluster = FALSE, 
                             show_rownames = T, 
                             return_heatmap = T)

#################################################################################

devtools::load_all("./monocle")
load("./ordering_genes.Rdata")
BEAM_res <- BEAM(cds_DGT[ordering_genes,], branch_point = 1, cores = 10)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

my_branched_heatmap = plot_genes_branched_heatmap(
  cds_DGT[row.names(subset(BEAM_res, qval<1e-4)),], 
  branch_point=1,
  num_clusters=3, cores=15, 
  use_gene_short_name=T,
  show_rownames=T)

plot_genes_branched_heatmap(cds_DGT[row.names(BEAM_res)[1:100]], branch_point = 1, num_clusters = 3, cores=15, use_gene_short_name=TRUE, show_rownames=TRUE)

BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_genes_branched_heatmap(cds_DGT[BEAM_genes,], branch_point = 1, num_clusters = 3, cores=15, use_gene_short_name=TRUE, show_rownames=TRUE)
p
save(cds_DGT,diff_test_res,
     file = '2output_of_cds.Rdata')
load(file = './2output_of_cds.Rdata')

methods <- c('duplicate', 'expression', 'cluster')
BEAM_res <- lapply(methods, function(method) {
  BEAM(cds, branch_point = 1, cores = 8, progenitor_method = method)
})
BEAM_res <- BEAM(cds_DGT, branch_point = 1,cores = 8,  progenitor_method = 'duplicate')
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write.csv(BEAM_res,file = "BEAM_res.csv")

plot_genes_branched_heatmap(cds_DGT[row.names(subset(BEAM_res, qval < 1e-120)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 8,
                            use_gene_short_name = F,
                            show_rownames = F)
#-------------------------------------------------------------------------------
library(org.Hs.eg.db)
library(clusterProfiler)
#
deg = all.markers %>% group_by(cluster) %>% top_n(n = 300, wt = avg_log2FC)
table(deg$cluster)
deg1 <- deg[deg$cluster=="0",]
degs.list=deg1$gene
#
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP", 
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP,showCategory = 15 )
erich.go.BP=erich.go.BP@result
write.table(erich.go.BP,"new(mCAF_CD36)erich.go.BP.deg.con.hf.txt",sep = "\t",col.names = NA)

erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.CC,showCategory = 15)


erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.MF,showCategory = 8)

library(dplyr)
markers <- all.markers |> group_by(cluster) |>
  filter(p_val_adj < 0.001) |>
  ungroup()

library(clusterProfiler)
gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(markers, gid, by=c('gene' = 'SYMBOL'))

gcSample=split(markers$ENTREZID, markers$cluster) 
xx1 <- compareCluster(gcSample,
                      fun = "enrichKEGG",
                      organism = "hsa", pvalueCutoff = 0.05
)
xx2 <- compareCluster(gcSample,
                      fun = "enrichGO",
                      OrgDb = "org.Hs.eg.db",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05
)

p <- dotplot(xx2)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

#BiocManager::install("GSVA")
library('GSEABase')
library(GSVA)
library(msigdbr)
#
m_df<- msigdbr(species = "human",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

Idents(Thyrocytes)="seurat_clusters"
exp=AverageExpression(Thyrocytes) 
#exp=AverageExpression(scRNA_harmony,add.ident = "orig.ident") 
exp=exp[["RNA"]]

GSVA_hall <- gsva(expr=exp, 
                  gset.idx.list=geneSets, 
                  mx.diff=T, 
                  kcdf="Gaussian", 
                  parallel.sz=10) 
head(GSVA_hall)
#
pheatmap::pheatmap(GSVA_hall,
                   cluster_rows = T,
                   cluster_cols =F,
                   
                   show_colnames=T,
                   scale = "row",
                   color =colorRampPalette(c("#0044BB", "white","#FF7744","red"))(100))
#-------------------------------------------------------------------------------
