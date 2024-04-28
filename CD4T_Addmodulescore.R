library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(monocle3)
library(clustree)
library(ggplot2)
library(harmony)
library(ggpubr)

rm(list = ls())
setwd("./CD4T")

##T细胞
load(file="./final_CD4T.Rdata")
DimPlot(CD4T, reduction = "umap", group.by='celltype',label = T,raster=FALSE)

apoptptic<- getGO("GO:0097190")
apoptptic <- apoptptic[[1]]

MIGRATION<- getGO("GO:0050900")
MIGRATION<- MIGRATION[[1]]

gs=list(
  naive = c( "CCR7", "TCF7", "LEF1", "SELL"), 
  cytotoxicity = c("PRF1", "IFNG", "GNLY"," NKG7"," GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW" , "CST7"), 
  exhaustion= c("LAG3"," TIGIT", "PDCD1", "CTLA4"," HAVCR2","TOX") ,
  metastasis = c('MYH2',	'HAS1',	'ADGRG7',	'WIF1',	'CSRP3',	'SFTPA2',	'EPYC',	'COL11A1',	'KRT6B',	'MMP13',	'MYF6',	
                 'SPRR2A',	'IL36G',	'KRT6C',	'KRT6A',	'AC244102.3',	'LINC01614',	'PAX5',	'KRT75',	'FDCSP',	'CR2',	
                 'HOTAIR',	'PSG3',	'PITX2',	'AF015262.1',	'SLC18A3',	'SPRR2F',	'TCL1A',	'CHAT',	'ROS1'),
  apoptosis=apoptptic,
  migration=MIGRATION
)
gs = lapply(gs, toupper)
sce =  AddModuleScore(object = CD4T,features = gs,name=c("naive","cytotoxicity", "metastasis","exhaustion","apoptosis","migration"))
colnames(sce@meta.data)
colnames(sce@meta.data)[12:17] <-c("naive","cytotoxicity","metastasis2", "exhaustion","apoptosis","migration")

###
library(ggpubr)
#
pnaive <- ggviolin(sce@meta.data, x="celltype", y="naive", width = 0.6, 
                   color = "black",
                   fill="celltype",
                   palette = "npg",
                   add = 'mean_sd',
                   xlab = F, 
                   x.text.angle = 45,
                   bxp.errorbar=T,
                   bxp.errorbar.width=0.5, 小
                   size=1, 
                   outlier.shape=NA, 
                   legend = "right")

pcytotoxicity <- ggviolin(sce@meta.data, x="celltype", y="cytotoxicity", width = 0.6, 
                          color = "black",
                          fill="celltype",
                          palette = "npg",
                          add = 'mean_sd',
                          xlab = F, 
                          x.text.angle = 45,
                          bxp.errorbar=T,
                          bxp.errorbar.width=0.5, 
                          size=1, 
                          outlier.shape=NA, 
                          legend = "right")

pmetastasis2 <- ggviolin(sce@meta.data, x="celltype", y="metastasis2", width = 0.6, 
                         color = "black",
                         fill="celltype",
                         palette = "npg",
                         add = 'mean_sd',
                         xlab = F, 
                         x.text.angle = 45,
                         bxp.errorbar=T,
                         bxp.errorbar.width=0.5, 
                         size=1, 
                         outlier.shape=NA, 
                         legend = "right")

pexhaustion <- ggviolin(sce@meta.data, x="celltype", y="exhaustion", width = 0.6, 
                        color = "black",
                        fill="celltype",
                        palette = "npg",
                        add = 'mean_sd',
                        xlab = F, 
                        x.text.angle = 45,
                        bxp.errorbar=T,
                        bxp.errorbar.width=0.5,
                        size=1, 
                        outlier.shape=NA, 
                        legend = "right")

papoptosis <- ggviolin(sce@meta.data, x="celltype", y="apoptosis", width = 0.6, 
                       color = "black",
                       fill="celltype",
                       palette = "npg",
                       add = 'mean_sd',
                       xlab = F, 
                       x.text.angle = 45,
                       bxp.errorbar=T,
                       bxp.errorbar.width=0.5,
                       size=1, 
                       outlier.shape=NA, 
                       legend = "right")

pmigration <- ggviolin(sce@meta.data, x="celltype", y="migration", width = 0.6, 
                       color = "black",
                       fill="celltype",
                       palette = "npg",
                       add = 'mean_sd',
                       xlab = F,
                       x.text.angle = 45,
                       bxp.errorbar=T,
                       bxp.errorbar.width=0.5, 
                       size=1, 
                       outlier.shape=NA,
                       legend = "right")
p2 <- pnaive + pcytotoxicity + pmetastasis2 + pexhaustion + papoptosis + pmigration
p2 
#
ggsave("eachindex_celltype.pdf",device = cairo_pdf, width = 11, height = 8, dpi = 300, p2)
#
my_comparisons <- list(c('ADULT','CAYA'))
ggboxplot(sce@meta.data, x="celltype", y="cytotoxicity", width = 0.6, 
          color = "black",
          fill="group",
          palette = "colour",
          xlab = F, 
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=1, 
          outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 

#
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272','#fc8452','#9a60b4','#ea7ccc')
ggviolin(sce@meta.data, x = "group", y = "cytotoxicity",
         color = "group",add = 'mean_sd',fill = 'group',
         add.params = list(color = "black")) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
#-------------------------------------------------------------------------------
my_comparisons <- list(c('ADULT','CAYA'))
p1<- ggviolin(sce@meta.data, x="group", y="naive", width = 0.6, 
              color = "black",
              fill="group",
              palette =c("#006699", "#CC3333"),
              add = 'mean_sd',
              xlab = F, 
              x.text.angle = 45,
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=1, 
              outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 


p2<- ggviolin(sce@meta.data, x="group", y="cytotoxicity", width = 0.6, 
              color = "black",
              fill="group",
              palette =c("#006699", "#CC3333"),
              add = 'mean_sd',
              xlab = F, 
              x.text.angle = 45,
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=1, 
              outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 



p3<- ggviolin(sce@meta.data, x="group", y="metastasis2", width = 0.6, 
              color = "black",
              fill="group",
              palette =c("#006699", "#CC3333"),
              add = 'mean_sd',
              xlab = F, 
              x.text.angle = 45,
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=1, 
              outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 



p4<- ggviolin(sce@meta.data, x="group", y="exhaustion", width = 0.6, 
              color = "black",
              fill="group",
              palette =c("#006699", "#CC3333"),
              add = 'mean_sd',
              xlab = F, 
              x.text.angle = 45,
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=1, 
              outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 


p5<- ggviolin(sce@meta.data, x="group", y="apoptosis", width = 0.6, 
              color = "black",
              fill="group",
              palette =c("#006699", "#CC3333"),
              add = 'mean_sd',
              xlab = F, 
              x.text.angle = 45,
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 小
              size=1, 
              outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 


p6<- ggviolin(sce@meta.data, x="group", y="migration", width = 0.6, 
              color = "black",
              fill="group",
              palette =c("#006699", "#CC3333"),
              add = 'mean_sd',
              xlab = F, 
              x.text.angle = 45,
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=1, 
              outlier.shape=NA) +
  theme(axis.text.x.bottom = element_text(angle = 45,vjust = 0.5,hjust = 1)) +
  NoLegend()+stat_compare_means(comparisons = my_comparisons,label = "p.signif") 
library(gridExtra)
p <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
p
ggsave("eachindex_CAYAvsADULT_CD4T.pdf",device = cairo_pdf, width = 10, height = 11, dpi = 300, p)
