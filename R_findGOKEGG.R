rm(list = ls())
setwd("./FUJI")

source("getGoTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")
save(GO_DATA, file = "GO_DATA.RData")
###GO
findGO <- function(pattern, method = "key"){
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method == "key"){
    pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
  } else if(method == "gene"){
    pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways) = "pathway"
  
  if(length(pathways) == 0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
}

getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
}

load("GO_DATA.RData") 
findGO("********") # looking for pathway name
findGO("****", method = "gene") # Find the pathway with the specified gene name.

getGO("GO:*******") # Gets the gene set of the specified GO ID.

###KEGG

# 方法一
if(!requireNamespace("BiocManager", quietly = TRUE))      
  install.packages("BiocManager") 
BiocManager::install("KEGGREST", version = "3.10")

library("KEGGREST")
gsInfo = keggGet('hsa*****')[[1]]
names(gsInfo)

geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])
geneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])
names(geneSet) = gsInfo$NAME
geneSet

#
#
library("rjson")
download.file("http://togows.dbcls.jp/entry/pathway/hsa*****/genes.json", "hsa*****.json")
json = fromJSON(file ="hsa*****.json")
geneSet = list(as.character(sapply(json[[1]], function(x) sapply(strsplit(x[1], ";"), function(x) x[1]))))
geneSet

getKEGG <- function(ID){
  
  library("KEGGREST")
  
  gsList = list()
  for(xID in ID){
    
    gsInfo = keggGet(xID)[[1]]
    if(!is.null(gsInfo$GENE)){
      geneSetRaw = sapply(strsplit(gsInfo$GENE, ";"), function(x) x[1])   
      xgeneSet = list(geneSetRaw[seq(2, length(geneSetRaw), 2)])          
      NAME = sapply(strsplit(gsInfo$NAME, " - "), function(x) x[1])
      names(xgeneSet) = NAME
      gsList[NAME] = xgeneSet 
    } else{
      cat(" ", xID, "No corresponding gene set in specific database.\n")
    }
  }
  return(gsList)
}

getKEGG('hsa*****')

