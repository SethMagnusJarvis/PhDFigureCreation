library(tidyverse)
library(EnrichmentBrowser)
library(EnhancedVolcano)
library(patchwork)
library(topGO)


F210IAdultSupport <- read_csv("F210IAdultSupport.csv")

BMOutput <- read_csv("BioMarTMouseWEntrez.csv")
AdultF210IDE <- read_csv("F210IAdultDE.csv") %>% 
  left_join(BMOutput, by=c("gene_id"="ensembl_gene_id"))

KSMat <- function(DEDF){
  tmp <- DEDF %>% drop_na("entrezgene_id") 
  tmp <- filter(tmp, ADJ.PVAL < 0.05)
  geneList <- tmp$PVAL 
  names(geneList) <- tmp$entrezgene_id
  
  selection <- function(allScore){ return(allScore < 0.05)}
  GOdata <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSel=selection,
                annot = annFUN.org, mapping="org.Mm.eg.db" )
  
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
  tabKS <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  
  #par(cex = 0.3)
  print(showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 2, useInfo = "def"))
  return(tabKS)
}




FisherMat <- function(DEDF){
  tmp <- DEDF %>% drop_na("entrezgene_id") 
  geneList <- tmp$PVAL 
  names(geneList) <- tmp$entrezgene_id
  
  selection <- function(allScore){ return(allScore < 0.05)}
  GOdata <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSel=selection,
                annot = annFUN.org, mapping="org.Mm.eg.db" )
  
  
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  tabFisher <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                        numChar = 120)
  
  print(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 2, useInfo = "def"))
  return(tabFisher)
}

svg("GONetworkKS.svg")
CTLVsMutKS <- KSMat(AdultF210IDE)
dev.off()
svg("GONetworkFisher.svg")
CTLVsMutFisher <- FisherMat(AdultF210IDE)
dev.off()

write_csv(CTLVsMutKS, "CTLvsMutGOKS.csv")
write_csv(CTLVsMutFisher, "CTLvsMutGOFisher.csv")
