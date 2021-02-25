library(biomaRt)
library(tidyverse)
library(KEGGREST)
library(pathview)
library(topGO)


BMOutput <- read_csv("BioMarTHumanWEntrez.csv")

CTLVsMutDECell <- read_csv("CtlVsMutDEXCellNeuron.csv") %>% 
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")


CTLVsC9DECell <- read_csv("CtlVsC9DEXCellNeuron.csv") %>% 
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")

CTLVsTAUDECell <- read_csv("CtlVsTAUDEXCellNeuron.csv") %>% 
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")

C9VsTAUDECell <- read_csv("C9VsTAUDEXCellNeuron.csv") %>% 
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")



#TopGo
KFMat <- function(DEDF){
  tmp <- DEDF %>% drop_na("entrezgene_id") 
  tmp <- filter(tmp, padj < 0.05)
  geneList <- tmp$pvalue
  names(geneList) <- tmp$entrezgene_id
  
  selection <- function(allScore){ return(allScore < 0.05)}
  GOdata <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSel=selection,
                annot = annFUN.org, mapping="org.Hs.eg.db" )
  
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
  tabKS <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  
  #par(cex = 0.3)
  print(showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 2, useInfo = "def"))
  return(tabKS)
}




FisherMat <- function(DEDF){
  tmp <- DEDF %>% drop_na("entrezgene_id") 
  geneList <- tmp$pvalue
  names(geneList) <- tmp$entrezgene_id
  
  selection <- function(allScore){ return(allScore < 0.05)}
  GOdata <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSel=selection,
                annot = annFUN.org, mapping="org.Hs.eg.db" )
  
  
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  
  tabFisher <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                        numChar = 120)
  
  print(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 2, useInfo = "def"))
  return(tabFisher)
}



svg("CTLVsMutKSXCell.svg")
CTLVsMutKS <- KFMat(CTLVsMutDECell)
dev.off()
svg("CTLVsMutFisherXCell.svg")
CTLVsMutFisher <- FisherMat(CTLVsMutDECell)
dev.off()

svg("CTLVsC9KSKSXCell.svg")
CTLVsC9KS <- KFMat(CTLVsC9DECell)
dev.off()
svg("CTLVsC9FisherXCell.svg")
CTLVsC9Fisher <- FisherMat(CTLVsC9DECell)
dev.off()

svg("CTLVsTAUKSXCell.svg")
CTLVsTAUKS <- KFMat(CTLVsTAUDECell)
dev.off()
svg("CTLVsTAUFisherXCell.svg")
CTLVsTAUFisher <- FisherMat(CTLVsTAUDECell)
dev.off()

svg("C9VsTAUKSXCell.svg")
C9VsTAUKS <- KFMat(C9VsTAUDECell)
dev.off()
svg("C9VsTAUFisherXCell.svg")
C9VsTAUFisher <- FisherMat(C9VsTAUDECell)
dev.off()

MUTC9KS <- full_join(CTLVsMutKS, CTLVsC9KS, by=c("GO.ID", "Term"), suffix=c(".CTLVsMUT", ".CTLVsC9"))
TAUVsC9xTAUKS <- full_join(CTLVsMutKS, CTLVsC9KS, by=c("GO.ID", "Term"), suffix=c(".CTLVsTAU", ".C9VsTAU"))
AllKS <- full_join(MUTC9KS, TAUVsC9xTAUKS, by=c("GO.ID", "Term"))


MUTC9Fisher <- full_join(CTLVsMutFisher, CTLVsC9Fisher, by=c("GO.ID", "Term"), suffix=c(".CTLVsMUT", ".CTLVsC9"))
TAUVsC9xTAUFisher <- full_join(CTLVsMutFisher, CTLVsC9Fisher, by=c("GO.ID", "Term"), suffix=c(".CTLVsTAU", ".C9VsTAU"))
AllFisher <- full_join(MUTC9Fisher, TAUVsC9xTAUFisher, by=c("GO.ID", "Term"))

write_csv(AllKS, "AllGOTermsKSXCell.csv")
write_csv(AllKS, "AllGOTermsFisherXCell.csv")

write_csv(CTLVsMutKS, "CTLVsMutGOKSXCell.csv")
write_csv(CTLVsMutFisher, "CTLVsMutGOFisherXCell.csv")
write_csv(CTLVsC9KS, "CTLVsC9GOKSXCell.csv")
write_csv(CTLVsC9Fisher, "CTLVsC9GOFisherXCell.csv")
write_csv(CTLVsTAUKS, "CTLVsTAUGOKSXCell.csv")
write_csv(CTLVsTAUFisher, "CTLVsTAUGOFisherXCell.csv")
write_csv(C9VsTAUKS, "C9VsTAUGOKSXCell.csv")
write_csv(C9VsTAUFisher, "C9VsTAUGOFisherXCell.csv")