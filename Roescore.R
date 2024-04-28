library(tidyverse)
library(readr)
metainfo <- Fibroblasts@meta.data %>% as.data.frame()
colnames(metainfo)
distribution_Roe = function(
    meta_data=metainfo,
    celltype_column = "majorCluster",
    celltype_level = NULL,
    condition_column = "tissue",
    condition_level = NULL,
    max_threshold = 2,
    min_threshold = 0.01,
    
    add_label = NA, 
    condition_label_angle = 0, #45
    condition_label_hjust = 0.5, #1
    celltype_color = NULL,relative_width = NULL, 
    tile_color = NA,
    tile_fill = c("#f6f8e6","#eb632e")
){
  library(tidyverse)
  
  colnames(meta_data)[which(colnames(meta_data) == celltype_column)] = "celltypE"
  colnames(meta_data)[which(colnames(meta_data) == condition_column)] = "conditioN"
  
  if(is.null(celltype_level)){
    meta_data$celltypE = as.character(meta_data$celltypE)
    meta_data$celltypE = factor(meta_data$celltypE,levels = sort(unique(meta_data$celltypE)))
  } else {
    meta_data$celltypE = factor(meta_data$celltypE,levels = celltype_level)
  }
  
  if(is.null(condition_level)) {
    meta_data$conditioN = as.character(meta_data$conditioN)
    meta_data$conditioN = factor(meta_data$conditioN,levels = sort(unique(meta_data$conditioN)))
  } else {
    meta_data$conditioN = factor(meta_data$conditioN,levels = condition_level)
  }
  
  contengency_table = xtabs(~celltypE+conditioN,data = meta_data)
  chisq <- chisq.test(as.matrix(contengency_table))
  Roe = chisq$observed / chisq$expected
  Roe = as.matrix(Roe) %>% as.data.frame()
  colnames(Roe)[3] = "value"
  Roe$old_value = Roe$value
  Roe$value[Roe$value > max_threshold] = max_threshold
  Roe$value[Roe$value < min_threshold] = 0
  
  Roe$sign = ""
  Roe$sign[Roe$value > 1] = "+++"
  Roe$sign[Roe$value > 0.8 & Roe$value <= 1] = "++"
  Roe$sign[Roe$value >= 0.2 & Roe$value <= 0.8] = "+"
  Roe$sign[Roe$value > 0 & Roe$value < 0.2] = "+/−"
  Roe$sign[Roe$value == 0] = "-"
  
  pb = Roe %>% ggplot(aes(x=conditioN,y=celltypE))+
    geom_tile(aes(fill = value),color = tile_color)+
    scale_y_discrete(expand = c(0,0),position = "right")+
    scale_x_discrete(expand = c(0,0))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y.right = element_text(size = 12,colour = "black"),
      axis.text.x.bottom = element_text(size = 12,colour = "black",angle = condition_label_angle,hjust = condition_label_hjust)
    )
  if (length(tile_fill) == 1) {
    pb=pb+scale_fill_viridis_c("Ro/e",option = tile_fill)
  } else if (length(tile_fill) == 2) {
    pb=pb+scale_fill_gradient("Ro/e",low = tile_fill[1],high = tile_fill[2])
  }
  if(is.na(add_label)) {
    pb=pb
  } else if (add_label == "number") {
    pb=pb+geom_text(aes(label=round(old_value,2)))
  } else if (add_label == "sign") {
    pb=pb+geom_text(aes(label=sign))
  }
  
  
  if (is.null(celltype_color) & is.null(relative_width)) {
    return(pb)
  } else if(!is.null(celltype_color) & !is.null(relative_width)) {
    library(patchwork)
    tmpdf = factor(levels(Roe$celltypE),levels = levels(Roe$celltypE)) %>% as.data.frame()
    colnames(tmpdf) = "celltype"
    
    pa = ggplot(data = tmpdf,aes(x=0,y=celltype))+
      geom_text(aes(label = celltype,color = celltype),hjust = 0,size = 5)+
      scale_color_manual(values = celltype_color)+
      scale_x_continuous(expand = c(0,0),limits = c(0,1))+
      theme_void()+
      theme(
        legend.position = "none"
      )
    
    pb=pb+theme(legend.position = "top",axis.text.y.right = element_blank())
    pc = pb+pa+plot_layout(widths = c(1,relative_width))
    return(pc)
  } else {
    return(print("celltype_color and relative_width do not match!"))
  }
}
celltype_color = c(
  "emCAF_LAMP5" = "#f19138","lpmCAF_CD36" = "#b27da3",
  "iCAF_CFD"="#547eaa","myoCAF_MYH11"="#79b8b3"
)
distribution_Roe(
  meta_data = metainfo,
  celltype_column = "celltype",
  celltype_level = names(celltype_color) %>% rev(),
  condition_column = "group2",
  add_label = "sign",
  celltype_color = celltype_color,relative_width = 0.7,
  tile_color = NA
)


