devtools::load_all("./monocle")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
rm(list = ls())
setwd("./monocle2")
load(file="./harmonythyrocytes-0.1.Rdata")
DimPlot(Thyrocytes, reduction = "umap", group.by='seurat_clusters',label = T,raster=FALSE)
DimPlot(Thyrocytes, reduction = "tsne", group.by='seurat_clusters',label = T,raster=FALSE)
#
data <- as(as.matrix(Thyrocytes@assays$RNA@counts), 'sparseMatrix')
#
pd <- new('AnnotatedDataFrame', data = Thyrocytes@meta.data)
#
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- detectGenes(monocle_cds, min_expr = 1)
monocle_cds <- monocle_cds[fData(monocle_cds)$num_cells_expressed>10,]
print(head(fData(monocle_cds)))

cds <- monocle_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#
disp_table<-dispersionTable(cds)
unsup_clustering_genes<-subset(disp_table,mean_expression>=0.1)
cds<-setOrderingFilter(cds,unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds,return_all=F)#norm_method='log'
#
cds<-reduceDimension(cds,max_components=2,num_dim=15,
                     reduction_method='tSNE',verbose=T)
cds<-clusterCells(cds)#num_clusters=6
plot_cell_clusters(cds,1,2)
table(pData(cds)$Cluster)
colnames(pData(cds))

table(pData(cds)$Cluster,pData(cds)$new)
plot_cell_clusters(cds,1,2)

save(cds,file='input_cds.Rdata')
#-------------------------------------------------------------------------------
rm(list=ls())
options(stringsAsFactors = F)

library(monocle)
library(Seurat)
load(file = 'input_cds.Rdata')

#
colnames(pData(cds))
table(pData(cds)$Cluster)
table(pData(cds)$Cluster,pData(cds)$seurat_clusters)
plot_cell_clusters(cds, 1, 2 )

#-------------------------------------------------------------------------------
pData(cds)$Cluster=pData(cds)$seurat_clusters
table(pData(cds)$Cluster)

Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()
#
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 
#
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )

#
diff_test_res1 <- diff_test_res
diff_test_res<- subset(diff_test_res, qval < 0.01)
ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval,decreasing = F)][1:2000]
save(ordering_genes,cds,file = "ordering_genes.Rdata")
#ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds,root_state=2)

plot_cell_trajectory(cds, color_by = "Cluster") +facet_wrap(~group, nrow = 1)
plot_cell_trajectory(cds, color_by = "Cluster") +facet_wrap(vars(Cluster,group), nrow = 2)
plot_cell_trajectory(cds, color_by = "Pseudotime")+facet_wrap(~group, nrow = 1)
plot_cell_trajectory(cds, color_by = "Cluster")+facet_wrap(~group, nrow = 1)
plot_cell_trajectory(cds, color_by = "Cluster")+facet_wrap(~State,nrow = 1)
plot_cell_trajectory(cds, color_by = "group")
plot_cell_trajectory(cds, color_by = "State")+facet_wrap(~group, nrow = 1)

ggsave('monocle_cell_trajectory_for_seurat.pdf')

length(cg)
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Pseudotime") 
ggsave('monocle_plot_genes_in_pseudotime_for_seurat.pdf')

save( cds,diff_test_res,
      file = '2output_of_cds.Rdata')
load(file = './2output_of_cds.Rdata')

my_cds_subset=cds
head(pData(my_cds_subset))
my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 10)

BEAM_res <- BEAM(cds[ordering_genes,], branch_point = 1, cores = 10)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

my_branched_heatmap = plot_genes_branched_heatmap(
  cds[row.names(subset(BEAM_res, qval<1e-4)),], 
  branch_point=1,
  num_clusters=3, cores=15, 
  use_gene_short_name=T,
  show_rownames=T)

plot_genes_branched_heatmap(cds[row.names(BEAM_res)[1:100]], branch_point = 1, num_clusters = 3, cores=15, use_gene_short_name=TRUE, show_rownames=TRUE)

BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds[BEAM_genes,], branch_point = 1, num_cluster = 3, cores=15, show_rownames = T, return_heatmap = T)

save( cds,diff_test_res,
      file = '2output_of_cds.Rdata')
load(file = './2output_of_cds.Rdata')
#-------------------------------------------------------------------------------
Time_diff <- differentialGeneTest(cds,fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 18)
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]
write.csv(Time_diff, file = "Time_diff.csv",row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds[Time_genes,], num_cluster = 3,  cores=15,show_rownames = T, return_heatmap = T)

ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

p$tree_row
Call:
  hclust(d = d, method = method)
Cluster method : ward.D2 
Number of objects: 2200

clusters <- cutree(p$tree_row, k = 3) #k值可以改变
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering, "Time_clustering_all.csv",row.names = F)