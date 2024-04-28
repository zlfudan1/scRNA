BiocManager::install(c("monocle"))
BiocManager::install(c('BiocGenerics', 'DelayedArray', 
                       'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github("cole-trapnell-lab/garnett",ref="monocle3")

library(garnett)
library(monocle3)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(Matrix)
setwd("./garnett")
load(file="./scRNA-umap-celltype.Rdata")
#
DimPlot(scRNA, reduction = "umap", group.by='celltype',label = T,raster=FALSE)
data = GetAssayData(scRNA, assay="RNA", slot = 'counts')
cell_metadata <- scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
# 
cds <- new_cell_data_set(data, 
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 10)
marker_check <- check_markers(cds, "./my_marker_file_all.txt",
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)

sc_seurat_obj_classifier <- train_cell_classifier(cds = cds,
                                                  marker_file = "./my_marker_file_all.txt",
                                                  db = org.Hs.eg.db,
                                                  cds_gene_id_type = "SYMBOL",
                                                  num_unknown = 10,
                                                  marker_file_gene_id_type = "SYMBOL",
                                                  min_observations = 50,
                                                  cores = 64 # linux)
pData(cds)$garnett_cluster <- pData(cds)$seurat_clusters

cds <- classify_cells(cds, sc_seurat_obj_classifier, 
                      db=org.Hs.eg.db, 
                      cluster_extend = TRUE, 
                      cds_gene_id_type = "SYMBOL")

cds.meta <- subset(pData(cds), 
                   select = c("cell_type","cluster_ext_type")) %>% as.data.frame()
sc_seurat_obj <- AddMetaData(scRNA, metadata = cds.meta)

DimPlot(sc_seurat_obj, group.by = "cluster_ext_type",raster=FALSE, 
        label = T, label.size = 3)

uterus <- sc_seurat_obj

df <- uterus@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = uterus@meta.data$celltype)%>%
  cbind(seurat_clusters = uterus@meta.data$seurat_clusters)%>%
  cbind(cluster_ext_type = uterus@meta.data$cluster_ext_type)
colnames(df)

cluster_order <- c('CD4+ T cells', 'CD8+ T cells', 'NK cells','B cells','Myeloid','Endothelial cells','Fibroblasts','Thyrocytes', 'Unknown')
scales::show_col(nature_col)

###
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                          '#D8A326', '#9E732B','#3172A5','#A4C7DA'),
                        c('CD4+ T cells', 'CD8+ T cells', 'NK cells','B cells','Myeloid','Endothelial cells','Fibroblasts','Thyrocytes', 'Unknown'))

#
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(cluster_ext_type, levels = cluster_order))) + 
  geom_point(size = 1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.8,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6))) + 
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+2.5, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+3),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), 
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.8, label = "UMAP1",
           color="black",size = 4) + 
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 1.4, label = "UMAP2",
           color="black",size = 4,angle=90) 