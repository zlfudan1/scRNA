library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)

binary_mtx <- read.csv('binary_mtx.csv', header = T, row.names = 1, encoding = 'uft8', check.names = F)
binary_mtx <- as.data.frame(t(binary_mtx))
head(binary_mtx)
dim(binary_mtx)


sce_data <- Fibroblasts
Idents(sce_data)="orig.ident"
sce_data=subset(sce_data,ident=c("ADULT1P","ADULT1T","ADULT2P","ADULT4P",
                                 "ADULT4T","ADULT5P","ADULT5T","ADULT6P","ADULT6T",
                                 "CAYA1P","CAYA1T","CAYA2P","CAYA2T","CAYA3P","CAYA3T",
                                 "CAYA4P","CAYA4T","CAYA5P","CAYA5T"))
table(sce_data$orig.ident)
table(sce_data$group)

Identify_TF <- function(obj,
                        group,
                        sample,
                        TF_binary_mtx,
                        coverage
                        
){
  
  meta <- obj@meta.data
  meta <- meta[, c(sample,group)]
  colnames(meta) <- c("sample","group")
  meta$cell <- rownames(meta)
  
  identifier <- unique(meta$sample)
  
  
  list_tfkeep <- lapply(identifier, function(x){
    cell <- meta$cell[meta$sample==x]
    bimatrix <- binary_mtx[,cell]
    tfcover <- apply(bimatrix, 1, function(x){sum(x)/length(x)})
    tfkeep <- names(tfcover)[tfcover>=coverage]
  })
  
  names(list_tfkeep) <- identifier
  all_tfkeep <- unique(do.call(c, list_tfkeep))
  
  matrix_tfpos <- sapply(identifier, function(x){
    cell <- meta$cell[meta$sample==x]
    bimatrix <- binary_mtx[all_tfkeep,cell]
    tfpos <- apply(bimatrix, 1, function(x){sum(x)/length(x)})>=coverage
    
  })  
  
  matrix_tfpos <- matrix_tfpos*1
  colnames(matrix_tfpos) <- identifier
  
  return(matrix_tfpos)
  
}


plot_matrix <- Identify_TF(obj = sce_data,
                           group = "group",
                           sample = "orig.ident",
                           TF_binary_mtx = binary_mtx,
                           coverage = 0.5)

anno_col <- data.frame(sample = unique(sce_data$orig.ident),
                       group = c(rep("ADULT",9), rep("CAYA",10)))


rownames(anno_col) <- anno_col$sample
plot_matrix <- plot_matrix[, rownames(anno_col)]

identifier_BM <- colnames(plot_matrix)[anno_col$group=="ADULT"]
incidence_matrix_BM <- plot_matrix[, identifier_BM] 
percentage_BM <- apply(incidence_matrix_BM, 1, function(x){sum(x)/9})


identifier_GM <- colnames(plot_matrix)[anno_col$group=="CAYA"]
incidence_matrix_GM <- plot_matrix[, identifier_GM] 
percentage_GM <- apply(incidence_matrix_GM, 1, function(x){sum(x)/10})

commnon <- rownames(plot_matrix)[apply(plot_matrix, 1, sum)>=15]
specific_BM <- rownames(plot_matrix)[percentage_BM>=1/3 & percentage_GM<=1/3]
specific_GM <- rownames(plot_matrix)[percentage_GM>=1/3 & percentage_BM<=1/3]

anno_row <- data.frame(row.names = 1, 
                       V1 = c(commnon,specific_BM,specific_GM),
                       type = c(rep("commnon", length(commnon)),
                                rep("specific_ADULT", length(specific_BM)),
                                rep("specific_CAYA", length(specific_GM))))

plot_matrix <- plot_matrix[c(commnon,specific_BM,specific_GM),]

groupcolor <- c('#50A5CC','#C0000A') 
names(groupcolor) <- unique(anno_col$group) 
samplecolor <- c("#ebe4d0","#afb79c","#98c39a","#66CC99","#339966","#6b8f78", "#879d95","#9eacb4","#9d8063",
                 "#f0c8b4","#f2c69c","#f2d58c","#f2c54e","#e7ae7c","#c7854d","#bf6a72","#e77a83","#d3a3a0","#e5a8b8")
names(samplecolor) <- anno_col$sample

typecolor <- c("#708090","#9dc2de",'#F3B1A0')
names(typecolor) <- unique(anno_row$type)

ann_colors <- list(group=groupcolor, 
                   sample= samplecolor,
                   type=typecolor) 

'#e0f3db','#a8ddb5','#4eb3d3','#0868ac','#084081'

p <- pheatmap::pheatmap(plot_matrix, 
                        clustering_method = "ward.D2",
                        cluster_cols = F, cluster_rows = F,
                        color = c("#e0f3db","#084081"),
                        annotation_row = anno_row,
                        annotation_col = anno_col,
                        annotation_colors = ann_colors,
                        fontsize = 8)
p


