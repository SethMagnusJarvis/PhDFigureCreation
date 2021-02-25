library(biomaRt)
library(tidyverse)
library(KEGGREST)
library(pathview)
library(topGO)


BMOutput <- read_csv("BioMarTHumanWEntrez.csv")

CTLVsMutDESolo <- read_csv("CtlVsMutDEBothCatAge.csv") %>% 
  rename_at(.vars = vars(ends_with(".Solo")),
            .funs = funs(sub("[.]Solo$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")


CTLVsC9DESolo <- read_csv("CtlVsC9DEBothCatAge.csv") %>% 
  rename_at(.vars = vars(ends_with(".Solo")),
            .funs = funs(sub("[.]Solo$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")

CTLVsTAUDESolo <- read_csv("CtlVsTAUDEBothCatAge.csv") %>% 
  rename_at(.vars = vars(ends_with(".Solo")),
            .funs = funs(sub("[.]Solo$", "", .))) %>% 
  left_join(BMOutput, by="external_gene_name")

C9VsTAUDESolo <- read_csv("C9VsTAUDEBothCatAge.csv") %>% 
  rename_at(.vars = vars(ends_with(".Solo")),
            .funs = funs(sub("[.]Solo$", "", .))) %>% 
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


svg("CTLVsMutKSCatAge.svg")
CTLVsMutKS <- KFMat(CTLVsMutDESolo)
dev.off()
svg("CTLVsMutFisherCatAge.svg")
CTLVsMutFisher <- FisherMat(CTLVsMutDESolo)
dev.off()

svg("CTLVsC9KSKSCatAge.svg")
CTLVsC9KS <- KFMat(CTLVsC9DESolo)
dev.off()
svg("CTLVsC9FisherCatAge.svg")
CTLVsC9Fisher <- FisherMat(CTLVsC9DESolo)
dev.off()

svg("CTLVsTAUKSCatAge.svg")
CTLVsTAUKS <- KFMat(CTLVsTAUDESolo)
dev.off()
svg("CTLVsTAUFisherCatAge.svg")
CTLVsTAUFisher <- FisherMat(CTLVsTAUDESolo)
dev.off()

svg("C9VsTAUKSCatAge.svg")
C9VsTAUKS <- KFMat(C9VsTAUDESolo)
dev.off()
svg("C9VsTAUFisherCatAge.svg")
C9VsTAUFisher <- FisherMat(C9VsTAUDESolo)
dev.off()


MUTC9KS <- full_join(CTLVsMutKS, CTLVsC9KS, by=c("GO.ID", "Term"), suffix=c(".CTLVsMUT", ".CTLVsC9"))
TAUVsC9xTAUKS <- full_join(CTLVsMutKS, CTLVsC9KS, by=c("GO.ID", "Term"), suffix=c(".CTLVsTAU", ".C9VsTAU"))
AllKS <- full_join(MUTC9KS, TAUVsC9xTAUKS, by=c("GO.ID", "Term"))

MUTC9Fisher <- full_join(CTLVsMutFisher, CTLVsC9Fisher, by=c("GO.ID", "Term"), suffix=c(".CTLVsMUT", ".CTLVsC9"))
TAUVsC9xTAUFisher <- full_join(CTLVsMutFisher, CTLVsC9Fisher, by=c("GO.ID", "Term"), suffix=c(".CTLVsTAU", ".C9VsTAU"))
AllFisher <- full_join(MUTC9Fisher, TAUVsC9xTAUFisher, by=c("GO.ID", "Term"))

write_csv(AllKS, "AllGOTermsKS.csv")

write_csv(AllKS, "AllGOTermsFisher.csv")

write_csv(CTLVsMutKS, "CTLVsMutGOKS.csv")
write_csv(CTLVsMutFisher, "CTLVsMutGOFisher.csv")
write_csv(CTLVsC9KS, "CTLVsC9GOKS.csv")
write_csv(CTLVsC9Fisher, "CTLVsC9GOFisher.csv")
write_csv(CTLVsTAUKS, "CTLVsTAUGOKS.csv")
write_csv(CTLVsTAUFisher, "CTLVsTAUGOFisher.csv")
write_csv(C9VsTAUKS, "C9VsTAUGOKS.csv")
write_csv(C9VsTAUFisher, "C9VsTAUGOFisher.csv")