library(SCopeLoomR)
library(AUCell)
library(SCENIC)

regulon <- read.csv("auc.csv", header = T, row.names = 1, check.names = F)
#===============================================================================
# 2、calculate CSI
#===============================================================================
CSI_matrix_cal <- function(regulon,
                           CSI_threshold,
                           module_k,
                           module_color=F,
                           Heatmap_col=NULL,
                           legend_parm = c("number","character"),
                           rect_color,
                           label_reg=NULL){
  
  #calculate CSI
  
  Mcor<-cor(regulon)
  n<-nrow(Mcor)
  
  CSI<-matrix(nrow=n,ncol=n)
  
  for (i in 1:n){
    for(j in 1:n){
      
      if(i==j) {
        
        CSI[i,j] <- 1
        
      } else{
        
        nodeA <- names(which(Mcor[i,]>= Mcor[i,j]-0.05))
        nodeB <- names(which(Mcor[,j]>= Mcor[i,j]-0.05))
        CSI[i,j]<- 1-((length(unique(c(nodeA,nodeB))))/n)
        
      }
      
    }
    
  }
  
  rownames(CSI)<-colnames(regulon)
  colnames(CSI)<-colnames(regulon)
  
  CSI_matrix <- as.data.frame(CSI)
  write.csv(CSI_matrix, file = './CSI_matrix.csv')
  
  
  #Heatmap-draw
  require(ComplexHeatmap)
  require(pheatmap)
  
  CSI[CSI <= CSI_threshold]=0
  
  
  if(is.null(Heatmap_col)){
    
    col<-colorRampPalette(c("#FAF9DA","#28245F"))(100)
    
  }else{
    
    col = Heatmap_col
    
  }
  
  
  x=pheatmap::pheatmap(CSI,
                       color=col,
                       clustering_method = "ward.D2",
                       show_rownames=FALSE,
                       show_colnames = FALSE,
                       cutree_rows = module_k,
                       cutree_cols = module_k)
  
  
  annotation_row <- data.frame(Cluster=factor(cutree(x$tree_row, module_k)))
  annotation_col <- data.frame(Cluster=factor(cutree(x$tree_col, module_k)))
  
  row_order <- annotation_row
  row_order$regulon <- rownames(row_order)
  row_order <- row_order[order(row_order$Cluster),]
  write.csv(row_order, file = "./Module.csv")
  
  
  anno_col = annotation_row
  anno_col$TF <- rownames(anno_col)
  
  index <- x$tree_col$order
  TFs <- x$tree_col$labels
  ord_TF <- c()
  for (i in index) {
    
    ord_TF <- append(ord_TF, TFs[i])
    
  }
  
  anno_col <- anno_col[ord_TF,]
  anno_col$Modules <- paste0("Module",anno_col$Cluster)
  
  
  
  if(module_color==F){
    
    calm = c("#C17E73","#688EC1", "#597873")
    
    module_num <- c(1:module_k)  
    cluster_color = setNames(calm[1:module_k],module_num) 
    
    
    
  }else{
    
    
    module_num <- c(1:module_k)  
    cluster_color = setNames(module_color,module_num) 
    
    
  }
  
  
  cluster_color_m <- as.data.frame(cluster_color)
  cluster_color_m$Modules <- paste0("Module",rownames(cluster_color_m))
  rownames(cluster_color_m) <- cluster_color_m$Modules
  
  
  if(legend_parm == "number"){
    
    heatmap_legend_param = list(color_bar = "continuous",
                                legend_direction = "vertical",
                                legend_width = unit(1, "cm"),
                                legend_height = unit(5, "cm"),
                                title = "Connection specificity index (CSI)",
                                title_position="leftcenter-rot",
                                border ="black",
                                at = c(0,0.2,0.4,0.6,0.8,1),
                                labels = c(0,0.2,0.4,0.6,0.8,1),
                                labels_gp = gpar(fontsize = 8,col='black',font = 3))
    
  }
  
  
  if(legend_parm == "character"){
    
    heatmap_legend_param = list(color_bar = "continuous",
                                legend_direction = "vertical",
                                legend_width = unit(1, "cm"),
                                legend_height = unit(5, "cm"),
                                title = "Connection specificity index (CSI)",
                                title_position="leftcenter-rot",
                                border ="black",
                                at = c(0,0.5,1),
                                labels = c("low","mid","high"),
                                labels_gp = gpar(fontsize = 8,col='black',font = 3))
  }
  
  
  
  
  
  hm = ComplexHeatmap::pheatmap(CSI, #ph$tree_row$order
                                annotation_row=annotation_row,
                                annotation_col=annotation_col,
                                clustering_method = "ward.D2",
                                show_rownames=FALSE,
                                show_colnames = FALSE,
                                color=col,
                                name = "ht",
                                treeheight_row = 20,
                                treeheight_col = 20,
                                annotation_names_col = F,
                                annotation_names_row = F,
                                annotation_legend=F,
                                annotation_colors= list(Cluster = cluster_color),
                                heatmap_legend_param = heatmap_legend_param)
  
  
  if(!is.null(label_reg)){
    
    label_reg <- as.data.frame(label_reg)
    hm <- hm+rowAnnotation(link = anno_mark(at = which(rownames(CSI) %in% label_reg$label_reg), 
                                            labels = label_reg$label_reg, labels_gp = gpar(fontsize = 8)))
    
    
  }
  
  
  draw(hm)
  
  ord = anno_col$Cluster
  dup = (which(!duplicated(ord)) - 1)
  fract = dup / nrow(anno_col)
  width =  c(fract[-1], 1) - fract
  
  decorate_heatmap_body("ht", {
    grid.rect(unit(fract, "native"), 
              unit(1-fract, "native"), 
              unit(width, "native"), 
              unit(width, "native"), 
              hjust = 0, 
              vjust = 1, 
              gp = gpar(col = rect_color, lty = 1, lwd = 2, fill=NA))
  })
  
  
  label_m <- unique(anno_col$Modules)
  cluster_color_m <- cluster_color_m[label_m, ]
  
  decorate_heatmap_body("ht", {
    
    
    grid.text(label_m, 
              unit(fract+0.25, "native"), 
              unit(1-fract-0.05, "native"), 
              gp=gpar(fontsize=15, col=cluster_color_m$cluster_color, fontface="bold"))
    
  })
  
  
  return(hm)
}

label <- read.csv("THYlist.csv")

label_reg <- label[,1]

module_TF <- CSI_matrix_cal(regulon = regulon,
                            CSI_threshold = 0.2,
                            module_k = 3,
                            module_color = F,
                            legend_parm = "number",
                            rect_color="red",
                            label_reg = label_reg )

pdf('PyscenicmoduleLable.pdf', width=8, height=6)
module_TF <- CSI_matrix_cal(regulon = regulon,
                            CSI_threshold = 0.2,
                            module_k = 3,
                            module_color = F,
                            legend_parm = "number",
                            rect_color="red",
                            label_reg = label_reg )
draw(module_TF)
dev.off()

#===============================================================================
# 4、tf-tf network
#===============================================================================
CSI_matrix <- read.csv("CSI_matrix.csv", header = T, row.names = 1,check.names = F)
CSI_net <- CSI_matrix
CSI_net$regulon <- rownames(CSI_net)

library(reshape2)
CSI_net<-melt(CSI_net, id.vars = c("regulon"),
              measure.vars = 1:169,
              variable.name = c('tf'),
              value.name = 'CSI')

write.csv(CSI_net, file = "CSI_net.csv")

out<-CSI_net[which(CSI_net[,3]>0.8),]
module <- read.csv("Module.csv", header = T, row.names = 1, check.names = F)

library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)

graph <- as_tbl_graph(out) %>% 
  mutate(interaction = centrality_degree(mode='out'),
         group=module$Cluster)


module_color = c("#C17E73","#688EC1", "#597873")


p2 <- ggraph(graph, layout = 'kk') + 
  geom_edge_fan(color='grey60',
                show.legend=T) + 
  geom_node_point(aes(size = interaction,fill=factor(group)),shape=21)+ 
  geom_node_text(aes(filter= interaction>5,label=name),size=2.5,repel = T)+
  theme_graph()+
  scale_fill_manual(values = module_color)
p2
#===============================================================================
# 5、Averge module activaty score mapped in TSNE
#===============================================================================
module <- read.csv("Module.csv", header = T, row.names = 1, check.names = F)

module_reg <- list()
for (i in 1:3) {
  
  module_r <- subset(module, Cluster==i)
  module_r <- module_r$regulon
  
  module_reg[[i]] <- module_r
}

#计算Module 平均auc
Average_mr <- matrix(nrow=nrow(regulon),ncol=3)
rownames(Average_mr) <- rownames(regulon)

for (i in 1:3) {
  
  Average_M <- regulon[,module_reg[[i]]]
  Average_M$module_A <- rowMeans(Average_M)
  Average_mr[,i] <- Average_M$module_A
  
}

colnames(Average_mr) <- c("M1","M2","M3")

library(Seurat)
sce_test <- Thyrocytes
reduc.df <- as.data.frame(Embeddings(sce_test@reductions$umap))#降维坐标
reduc.df <- cbind(reduc.df, Average_mr)

Mosule_ave <- list()
for (i in 1:3) {
  data = reduc.df[,c(1,2,i+2)]
  colnames(data)[3] <- "Module"
  whitePurple = c('#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b')
  p = ggplot(data, aes(x=UMAP_1, y=UMAP_2))+
    ggrastr::rasterise(geom_point(size=1.5, color="black"), dpi = 200)+
    ggrastr::rasterise(geom_point(size=1, aes(color=Module)), dpi = 200) + 
    theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
    scale_color_gradientn(colours = whitePurple)+
    ggtitle(colnames(Average_mr)[i])
  
  Mosule_ave[[i]] <- p
  
}

p <- Seurat::CombinePlots(Mosule_ave, ncol = 3)
p




