library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(monocle)
library(clustree)
library(ggplot2)
library(ggrepel)
rm(list = ls())
setwd("./samples")
load(file ='ADULT1T.Rdata')
load(file ='ADULT1P.Rdata')
load(file ='ADULT2T.Rdata')
load(file ='ADULT2P.Rdata')
load(file ='ADULT3T.Rdata')
load(file ='ADULT3P.Rdata')
load(file ='ADULT4T.Rdata')
load(file ='ADULT4P.Rdata')
load(file ='ADULT5T.Rdata')
load(file ='ADULT5P.Rdata')
load(file ='ADULT6T.Rdata')
load(file ='ADULT6P.Rdata')
load(file ='CAYA1T.Rdata')
load(file ='CAYA1P.Rdata')
load(file ='CAYA2T.Rdata')
load(file ='CAYA2P.Rdata')
load(file ='CAYA3T.Rdata')
load(file ='CAYA3P.Rdata')
load(file ='CAYA4T.Rdata')
load(file ='CAYA4P.Rdata')
load(file ='CAYA5T.Rdata')
load(file ='CAYA5P.Rdata')

sce.big <- merge(x = ADULT1T, 
                 y = c(ADULT1P,ADULT2T,ADULT2P,ADULT3T,ADULT3P,ADULT4T,ADULT4P,ADULT5T,ADULT5P,ADULT6T,ADULT6P,
                       CAYA1T,CAYA1P,CAYA2T,CAYA2P,CAYA3T,CAYA3P,CAYA4T,CAYA4P,CAYA5T,CAYA5P), 
                 add.cell.ids = c('ADULT1T','ADULT1P','ADULT2T','ADULT2P','ADULT3T','ADULT3P','ADULT4T','ADULT4P','ADULT5T','ADULT5P','ADULT6T','ADULT6P',
                                  'CAYA1T','CAYA1P','CAYA2T','CAYA2P','CAYA3T','CAYA3P','CAYA4T','CAYA4P','CAYA5T','CAYA5P'), 
                 project = "PTC22")
sce.big
sce.list <- SplitObject(sce.big, split.by = "orig.ident")
for (i in 1:length(sce.list)) {
  sce.list[[i]] <- NormalizeData(sce.list[[i]],verbose=FALSE)
  sce.list[[i]] <- FindVariableFeatures(sce.list[[i]], selection.method = "vst",nfeatures = 2000,verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = sce.list)
scRNA.anchors <- FindIntegrationAnchors(object.list = sce.list, anchor.features = features)
scRNA <- IntegrateData(anchorset = scRNA.anchors)
save(scRNA,file = 'scRNA.integratedata.Rdata')

scRNA$HT<-scRNA$orig.ident
scRNA$HT<-as.character(scRNA$HT)
scRNA$HT[scRNA$HT == 'CAYA1P'] <- 'NHT'
scRNA$HT[scRNA$HT == 'CAYA1T'] <- 'NHT'
scRNA$HT[scRNA$HT == 'CAYA2P']  <- "HT"
scRNA$HT[scRNA$HT == 'CAYA2T'] <- 'HT'
scRNA$HT[scRNA$HT == 'CAYA3P'] <- 'NHT'
scRNA$HT[scRNA$HT == 'CAYA3T'] <- 'NHT'
scRNA$HT[scRNA$HT == 'CAYA4P'] <- 'HT'
scRNA$HT[scRNA$HT == 'CAYA4T'] <- 'HT'
scRNA$HT[scRNA$HT == 'CAYA5P'] <- 'HT'
scRNA$HT[scRNA$HT == 'CAYA5T'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT1P'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT1T'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT2T'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT2P'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT3T'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT3P'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT4T'] <- 'NHT'
scRNA$HT[scRNA$HT == 'ADULT4P'] <- 'NHT'
scRNA$HT[scRNA$HT == 'ADULT5T'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT5P'] <- 'HT'
scRNA$HT[scRNA$HT == 'ADULT6T'] <- 'NHT'
scRNA$HT[scRNA$HT == 'ADULT6P'] <- 'NHT'

scRNA$group<-scRNA$orig.ident
scRNA$group<-as.character(scRNA$group)
scRNA$group[scRNA$group == 'CAYA1P'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA1T'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA2P'] <- "CAYA"
scRNA$group[scRNA$group == 'CAYA2T'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA3P'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA3T'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA4P'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA4T'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA5P'] <- 'CAYA'
scRNA$group[scRNA$group == 'CAYA5T'] <- 'CAYA'
scRNA$group[scRNA$group == 'ADULT1P'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT1T'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT2T'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT2P'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT3T'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT3P'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT4T'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT4P'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT5T'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT5P'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT6T'] <- 'ADULT'
scRNA$group[scRNA$group == 'ADULT6P'] <- 'ADULT'

scRNA$metastasis<-scRNA$orig.ident
scRNA$metastasis<-as.character(scRNA$metastasis)
scRNA$metastasis[scRNA$metastasis == 'CAYA1P'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'CAYA1T'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'CAYA2P'] <- "N0"
scRNA$metastasis[scRNA$metastasis == 'CAYA2T'] <- 'N0'
scRNA$metastasis[scRNA$metastasis == 'CAYA3P'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'CAYA3T'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'CAYA4P'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'CAYA4T'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'CAYA5P'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'CAYA5T'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'ADULT1P'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT1T'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT2T'] <- 'N0'
scRNA$metastasis[scRNA$metastasis == 'ADULT2P'] <- 'N0'
scRNA$metastasis[scRNA$metastasis == 'ADULT3T'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT3P'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT4T'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT4P'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT5T'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT5P'] <- 'N1b'
scRNA$metastasis[scRNA$metastasis == 'ADULT6T'] <- 'N1a'
scRNA$metastasis[scRNA$metastasis == 'ADULT6P'] <- 'N1a'

scRNA$spatialdistr<-scRNA$orig.ident
scRNA$spatialdistr<-as.character(scRNA$spatialdistr)
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA1P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA1T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA2P'] <- "Normal"
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA2T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA3P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA3T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA4P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA4T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA5P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'CAYA5T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT1P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT1T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT2T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT2P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT3T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT3P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT4T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT4P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT5T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT5P'] <- 'Normal'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT6T'] <- 'Tumor'
scRNA$spatialdistr[scRNA$spatialdistr == 'ADULT6P'] <- 'Normal'

DefaultAssay(scRNA) <- "RNA"

scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

DefaultAssay(scRNA) <- "integrated"
scRNA<- ScaleData(scRNA, features = VariableFeatures(scRNA))

CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))

g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = VariableFeatures(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = VariableFeatures(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)

scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
p
VlnPlot(scRNAa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)
ggsave("cellcycle_pca.png", p, width = 8, height = 6)
scRNAb <- ScaleData(scRNAa, vars.to.regress = c("S.Score", "G2M.Score"))
scRNA <- scRNAb

scRNA<- ScaleData(scRNA, features = VariableFeatures(scRNA))


scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
scRNA

print(scRNA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA, dims = 1:2, reduction = "pca")

DimHeatmap(scRNA, dims = 1:20, cells = 500)

ElbowPlot(scRNA)
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:15)
JackStrawPlot(scRNA,dims = 1:15)

DimPlot(scRNA, reduction = "pca")+ NoLegend()

scRNA <- FindNeighbors(scRNA, dims = 1:15)
###------------------------------------------------------------------------------

res.used <- seq(0.1,1,by=0.2)
res.used

for(i in res.used){
  scRNA <- FindClusters(object = scRNA, verbose = T, resolution = res.used)
}

library(clustree)
clus.tree.out <- clustree(scRNA) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")

clus.tree.out
##-------------------------------------------------------------------------------
scRNA<- FindClusters(scRNA, resolution = 1) 
pc.num=1:20   
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

scRNA = RunTSNE(scRNA, dims = pc.num)

plot1 = DimPlot(scRNA, reduction = "tsne", label=T) 
plot2 = DimPlot(scRNA, reduction = "tsne", group.by='orig.ident') 
plot3 = DimPlot(scRNA, reduction = "tsne", group.by='HT')
plot4 = DimPlot(scRNA, reduction = "tsne", group.by='group')
plot5 = DimPlot(scRNA, reduction = "tsne", group.by='metastasis') 
plot6 = DimPlot(scRNA, reduction = "tsne", group.by='spatialdistr')
plot7 = DimPlot(scRNA, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
plot8 <- DimPlot(scRNA, reduction = "tsne", group.by='celltype', label = T)

plotc <- plot1+plot2+plot3+plot4+plot5+plot6
plotc

#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
#group_by
plot7 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
plot8 = DimPlot(scRNA, reduction = "umap", group.by='HT')
plot9 = DimPlot(scRNA, reduction = "umap", group.by='age')
plot10 = DimPlot(scRNA, reduction = "umap", group.by='metastasis')
plot11 = DimPlot(scRNA, reduction = "umap", group.by='spatialdistr')
plot12 = DimPlot(scRNA, reduction = "umap", group.by='seurat_clusters',split.by = "age",label = T,raster=FALSE)
plot13 = DimPlot(scRNA, reduction = "umap", group.by='celltype',label = T,raster=FALSE)
plot14 = DimPlot(scRNA, reduction = "umap", group.by='celltype2',label = T)
#combinate
plotd <- plot7+plot8+plot9+plot10+plot12+plot13
plotd
plot13+plot12
#tSNE and UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
plotc
save(scRNA,all.markers,file="scRNA-umap-celltype.Rdata")

####----------------------------------------------------------------------------
DefaultAssay(scRNA) <- "RNA"
rb.genes <- rownames(scRNA)[grep("^RP[SL]",rownames(scRNA))]
C<-GetAssayData(object = scRNA, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
scRNA <- AddMetaData(scRNA, percent.ribo, col.name = "percent.ribo")
scRNAb <- ScaleData(scRNA, vars.to.regress =c("percent.ribo"))
scRNA <- scRNAb

####----------------------------------------------------------------------------

diff.wilcox = FindAllMarkers(scRNA,logfc.threshold = 0.25,only.pos = T,test.use = "MAST")
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)

ClusterMarker_noRibo <- all.markers[!grepl("^RP[SL]", all.markers$gene, ignore.case = F),]
ClusterMarker_noRibo_noMito <- ClusterMarker_noRibo[!grepl("^MT-", ClusterMarker_noRibo$gene, ignore.case = F),]
top = 25   
TopMarkers1 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo_noMito %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarkerTOP25<- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarkerTOP25,'TopMarkersTOP25.csv', row.names=F)

###-----------------------------------------------------------------------------

scRNA1=scRNA
Idents(scRNA1)="seurat_clusters"
scRNA1=RenameIdents(scRNA1,c('0'='B_cells','1'='CD4+T_cells','2'='CD8+T_cells','3'='Thyrocytes','4'='CD4+T_cells','5'='CD8+T_cells','6'='CD4+T_cells','7'='Thyrocytes','8'='CD4+T_cells','9'='Thyrocytes',
                             '10'='Endothelial','11'='Myeloid','12'='Thyrocytes','13'='NK','14'='NK','15'='B_cells','16'='Myeloid','17'='Myeloid','18'='Fibroblasts',
                             '19'='STMN1_cells','20'='Fibroblasts','21'='Fibroblasts','22'='Myeloid','23'='Endothelial','24'='CD4+T_cells','25'='Fibroblasts','26'='B_cells','27'='Myeloid',
                             '28'='Myeloid','29'='Thyrocytes','30'='CD4+T_cells','31'='B_cells'))
scRNA1@meta.data$celltype=scRNA1@active.ident
scRNA1 <- scRNA1
DimPlot(scRNA1, reduction = "umap", group.by='celltype',label = T,raster=FALSE)
DimPlot(scRNA1, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
scRNA <- scRNA1
save(scRNA1,all.markers,file="new_scRNA-umap-celltype.Rdata")
load('./new_scRNA-umap-celltype.Rdata')
##---------------------------------------------------------------------------------------------------------------
##---------------------------------------------------

genes_to_check = c('PTPRC', 'CD3D', 'CD3E', #T cell
                   'CD4','CCR7','IL7R','FOXP3','CTLA4',#CD4
                   'CD8A','CD8B',#CD8
                   'CD19', 'CD79A', 'MS4A1' ,#B cell
                   'IGHG1', 'MZB1', 'SDC1',#plasma cell
                   'KLRB1','NCR1','GNLY','NKG7','KLRD1','CX3CR1','FGFBP2','FCG3RA',# NK
                   'LYZ',#myeloid
                   'CD68', 'CD163', 'CD14', 'C1QA','C1QB', 'CD86','CD206','VSIG4',#TAM
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'LAMP3', 'IDO1','IDO2',# DC3 
                   'CD1E','CD1C', # DC2
                   'IRF7', 'LILRA4','IL3RA',# pDC
                   'RGS5','TAGLN', ## fibo 
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' ,'TG','CLU','CD24')##thyrocyte
p <- DotPlot(scRNA,features =genes_to_check,assay = 'RNA' ,group.by = "seurat_clusters", dot.scale = 3.0)+coord_flip()
p
##---------------------------------------------------
uterus <- scRNA1
uterus <- scRNA1

df <- uterus@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = uterus@meta.data$celltype)%>%
  cbind(seurat_clusters = uterus@meta.data$seurat_clusters)
colnames(df)

###
cluster_order <- c('0','1',"2","3", "4","5","6", "7","8","9","10","11","12",'13',
                   "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", 
                   "24", "25", "26", "27", "28", "29", "30", "31")
###
col_cluster <- setNames(c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A', 
                          '#D8A326', '#9E732B','#656565','#A4C7DA', '#3172A5', 
                          '#ADD28A','#40973E','#EA9393', '#CC161B', '#EFB773',
                          '#E87723','#C4AECF', '#603B8C', '#A3542B','#668DC6',
                          '#7C157E', '#3D5280', '#E35D49','#4BB7BC','#EE56FF', 
                          '#CACA58', '#E130CE','#CE4288','#FF549C', '#60D1FF', 
                          '#7A3C3F','#D0DF5B'),
                        c('0','1',"2","3", "4","5","6", "7","8","9","10","11","12",'13',
                          "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", 
                          "24", "25", "26", "27", "28", "29", "30", "31"))
###
cell <- df %>%group_by(seurat_clusters) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell) <-cell$seurat_clusters
A <- cell[cluster_order,]
a <- c(0:31)
A$ID <- a
A$ID <- as.factor(A$ID)
colnames(A) <- c( "seurat_clusters","UMAP_1", "UMAP_2" ,"ID")
colnames(df)
ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(seurat_clusters, levels = cluster_order))) + 
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
  annotate("text", x = min(df$UMAP_1) + 0.7, y = min(df$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 3, fontface = 'bold', family = 'serif') + 
  annotate("text", x = min(df$UMAP_1) -0.5, y = min(df$UMAP_2) + 0.7, label = "UMAP_2",
           color="black",size = 3, fontface = 'bold', family = 'serif', angle = 90) +
  geom_text_repel(data=A, aes(label=ID),color="black", size=6, point.padding = 0.3)
ggsave("all_umap.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)

#
B <- A
B$x <- 1
B$lab <- c(31:0)
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
ggsave("all_umap_legend.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300)
###--------------------------------------------------------------------------

marker <- c('PTPRC', 'CD3D', 'CD3E', #T cell
            'CD4','CCR7','IL7R','FOXP3','CTLA4',#CD4
            'CD8A','CD8B',#CD8
            'GNLY','NKG7','KLRD1','FGFBP2',# NK
            'CD19', 'CD79A', 'MS4A1', 'CD22',#B cell
            'IGHG1', 'CD38','TNFRSF17','IGHG4',#plasma cell
            'LYZ',#myeloid
            'CD68', 'CD163', 'CD14', 'C1QA','C1QB', 'CD86','VSIG4',#TAM
            'S100A9', 'S100A8', # monocyte
            'LAMP3', 'IDO1','IDO2',# DC3 
            'CD1E','CD1C', # DC2
            'LILRA4','IL3RA',# pDC
            'PECAM1', 'VWF',  ## endo 
            'RGS5','TAGLN', ## fibo 
            'EPCAM' ,'TG','CLU','CD24')
p <- DotPlot(scRNA1, features = marker,
             cols = my36colors, group.by = "seurat_clusters")+coord_flip()
#
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E59CC4', '#AB3282', '#23452F', '#BD956A',  '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69',  '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)

exp$id <- factor(exp$id,levels = c('1',"4","6","8","30",
                                   "2","5","24",
                                   '13',"14",
                                   '0',"15","26","31",
                                   "11","16", "17", "22","27", "28",
                                   "10","23",
                                   "18","20", "21","25",
                                   "3",  "7","9","12","29",
                                   "19"))
exp$features.plot<- factor(exp$features.plot,levels = c('PTPRC', 'CD3D', 'CD3E', #T cell
                                                        'CD4','CCR7','IL7R','FOXP3','CTLA4',#CD4
                                                        'CD8A','CD8B',#CD8
                                                        'GNLY','NKG7','KLRD1','FGFBP2',# NK
                                                        'CD19', 'CD79A', 'MS4A1', 'CD22',#B cell
                                                        'IGHG1', 'CD38','TNFRSF17','IGHG4',#plasma cell
                                                        'LYZ',#myeloid
                                                        'CD68', 'CD163', 'CD14', 'C1QA','C1QB', 'CD86','VSIG4',#TAM
                                                        'S100A9', 'S100A8', # monocyte
                                                        'LAMP3', 'IDO1','IDO2',# DC3 
                                                        'CD1E','CD1C', # DC2
                                                        'LILRA4','IL3RA',# pDC
                                                        'PECAM1', 'VWF',  ## endo 
                                                        'RGS5','TAGLN', ## fibo 
                                                        'EPCAM' ,'TG','CLU','CD24'))
ggplot(exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#6699CC','#FFFF99','#CC3333'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  geom_vline(xintercept=c(5.5,8.5,10.5,14.5,20.5,22.5,26.5,31.5), linetype="dotted",size=1)+
  geom_rect(aes(xmin = 1 - 0.5, xmax = 20 + 0.5, ymin = 1 - 0.5, ymax = 1 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 1 - 0.5, xmax = 5 + 0.5, ymin = 2 - 0.5, ymax = 8 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 5 + 0.5, xmax = 8 + 0.5, ymin = 8 + 0.5, ymax =10 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 8 + 0.5, xmax = 10 + 0.5, ymin = 10 + 0.5, ymax =14 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 10 + 0.5, xmax = 14 + 0.5, ymin = 14 + 0.5, ymax =22 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 14 + 0.5, xmax = 20 + 0.5, ymin = 22 + 0.5, ymax =39 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 20 + 0.5, xmax = 22 + 0.5, ymin = 39 + 0.5, ymax =41 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 22 + 0.5, xmax = 26 + 0.5, ymin = 41 + 0.5, ymax =43 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  geom_rect(aes(xmin = 26 + 0.5, xmax = 31 + 0.5, ymin = 43 + 0.5, ymax =47 + 0.5),
            fill = "transparent", color = "black", size = 0.5)+
  
  theme(legend.direction = "horizontal", legend.position = "bottom")
ggsave("all_bubble_plot.pdf",device = cairo_pdf, width = 11, height = 10, dpi = 300)
###---------------------------------------------------------------------------
#Ucell
library(Seurat)
BiocManager::install("UCell")
library(UCell)
DimPlot(scRNA, reduction = "umap", group.by='celltype',label = T,raster=FALSE)

markers <- list()
markers$CD4<- c("CD4", "CCR7", "IL7R", "CTLA4", "FOXP3")
markers$CD8 <- c("CD8A", "CD8B")
markers$NK <- c("FGFBP2", "KLRD1", "NKG7", "GNLY")
markers$B <- c("MS4A1", "CD79A", "CD19")
markers$Myeloid <- c("LYZ", "CD68", "CD163", "CD14", "C1QA", "C1QB","CD86", "VSIG4", "S100A9", "S100A8", "LAMP3", "IDO1","IDO2", "CD1E", "CD1C", "LILRA4", "IL3RA")
markers$Endo <- c("PECAM1", "VWF")
markers$Fibro <- c("TAGLN","RGS5")
markers$Thyro <- c("CD24", "CLU", "TG", "EPCAM")

#
marker_score <- AddModuleScore_UCell(scRNA,
                                     features=markers)
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")
#
library(stringr)
library(ggplot2)
library(viridis)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
table(marker_score@meta.data[["CD4_UCell"]])
FeaturePlot(marker_score,features = a,order = T, ncol = 4, cols = viridis(256),raster=FALSE)

##---------------------------------------------------
#
table(scRNA@meta.data$seurat_clusters)
cellnumber <- table(scRNA@meta.data$seurat_clusters)
write.csv(cellnumber, "cellnumber.csv", row.names = T)
d<- as.data.frame(cellnumber)
ggplot(data = d, aes(Var1,Freq))+ geom_col()+ labs(x = 'seurat_clusters', y = 'No. of cells')
#
table(scRNA$seurat_clusters)
prop.table(table(scRNA$celltype))
Cellratio <- prop.table(table(scRNA$celltype, scRNA$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#F68080")

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
####----------------------------------------------------------------------------
###
table(scRNA@meta.data$orig.ident)
patientcellnumber <- table(scRNA@meta.data$orig.ident)
write.csv(patientcellnumber, "patientcellnumber.csv", row.names = T)

table(scRNA@meta.data$seurat_clusters)
immunecellnumber <- table(scRNA@meta.data$seurat_clusters)
write.csv(immunecellnumber, "immunecellnumber.csv", row.names = T)

table(scRNA@meta.data$HT)
HTcellnumber <- table(scRNA@meta.data$HT)
write.csv(HTcellnumber, "HTcellnumber.csv", row.names = T)

table(scRNA@meta.data$group)
groupcellnumber <- table(scRNA@meta.data$group)
write.csv(groupcellnumber, "groupcellnumber.csv", row.names = T)

table(scRNA@meta.data$metastasis)
metastasiscellnumber <- table(scRNA@meta.data$metastasis)
write.csv(metastasiscellnumber, "metastasiscellnumber.csv", row.names = T)

table(scRNA@meta.data$spatialdistr)
spatialdistrcellnumber <- table(scRNA@meta.data$spatialdistr)
write.csv(spatialdistrcellnumber, "spatialdistrcellnumber.csv", row.names = T)

setwd("./Figure 1")
data <- read.csv("cell_number_ratio_all.csv", header = T)
ggplot(data, aes(x = var, y= count, fill = flag))+ 
  geom_col()+ labs(x= NULL, y="Number of Count", order = T)+
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
  #scale_fill_brewer(palette = 'Set1') 
  scale_fill_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500"))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 9.5), linetype= "dashed", linewidth = 0.5)
-----------------------------
####
table(scRNA$orig.ident)
prop.table(table(scRNA$celltype))
Cellratio <- prop.table(table(scRNA$celltype, scRNA$orig.ident), margin = 2)#
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#F68080")

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

##
table(scRNA$HT)
prop.table(table(scRNA$celltype))
Cellratio_HT <- prop.table(table(scRNA$celltype, scRNA$HT), margin = 2)
Cellratio_HT <- as.data.frame(Cellratio_HT)
write.csv(Cellratio_HT, "Cellratio_HT.csv", row.names = T)

##
table(scRNA$group)
prop.table(table(scRNA$celltype))
Cellratio_group <- prop.table(table(scRNA$celltype, scRNA$group), margin = 2)
Cellratio_group <- as.data.frame(Cellratio_group)
write.csv(Cellratio_group, "Cellratio_group.csv", row.names = T)

##
table(scRNA$metastasis)
prop.table(table(scRNA$celltype))
Cellratio_metastasis <- prop.table(table(scRNA$celltype, scRNA$metastasis), margin = 2)
Cellratio_metastasis <- as.data.frame(Cellratio_metastasis)
write.csv(Cellratio_metastasis, "Cellratio_metastasis.csv", row.names = T)

##
table(scRNA$spatialdistr)
prop.table(table(scRNA$celltype))
Cellratio_spatialdistr <- prop.table(table(scRNA$celltype, scRNA$spatialdistr), margin = 2)
Cellratio_spatialdistr <- as.data.frame(Cellratio_spatialdistr)
write.csv(Cellratio_spatialdistr, "Cellratio_spatialdistr.csv", row.names = T)

setwd("./Figure 1")
data <- read.csv("cell_ratio_all.csv", header = T)
data$var <- factor(data$var, levels=c('CD4+ T cells', 'CD8+ T cells', 'NK cells','B cells','Myeloid','Endothelial cells','Fibroblasts','Thyrocytes'))
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
  #scale_fill_brewer(palette = 'Set1') 
  scale_fill_manual(values = c('#349573', '#C7591A', '#6F6BA5','#CE2B7D', '#6A9D3A','#D8A326', '#9E732B','#3172A5'))+
  geom_vline(xintercept=c(2.5, 4.5, 6.5, 9.5), linetype= "dashed", linewidth = 0.5)

##------------------------------------------------------------------------------
