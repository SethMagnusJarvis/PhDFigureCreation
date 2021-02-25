library(tidyverse)
library(topGO)

BMOutput <- read_csv("BioMarTMouseWEntrez.csv")

PR12MDE <- read_csv("DifferentialExpressionPR12MNoType.csv")
GR12MDE <- read_csv("DifferentialExpressionGR12MNoType.csv")

PR12MDEEntrez <- left_join(PR12MDE, BMOutput, by=c("gene_id"="ensembl_gene_id"))
GR12MDEEntrez <- left_join(GR12MDE, BMOutput, by=c("gene_id"="ensembl_gene_id"))

KSMat <- function(DEEntrez){
  tmp <- DEEntrez %>% drop_na("entrezgene_id") 
  tmp <- filter(tmp, ADJ.PVAL < 0.05)
  geneList <- tmp$PVAL
  names(geneList) <- tmp$entrezgene_id
  
  GOdata <- new("topGOdata",
                ontology = "MF",
                allGenes = geneList,
                geneSelectionFun = function(x)x,
                annot = annFUN.org , mapping = "org.Mm.eg.db")
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
  tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  
  par(cex = 0.3)
  showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 4, useInfo = "def")
  return(tab)
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


svg("PR12MGOKS.svg")
GOPR12KS <- KSMat(PR12MDEEntrez)
dev.off()

svg("PR12MGOFisher.svg")
GOPR12Fisher <- FisherMat(PR12MDEEntrez)
dev.off()

svg("GR12MGOKS.svg")
GOGR12KS <- KSMat(GR12MDEEntrez)
dev.off()

svg("GR12MGOFisher.svg")
GOGR12Fisher <- FisherMat(GR12MDEEntrez)
dev.off()

write_csv(GOPR12KS, "PR12MGOKS.csv")
write_csv(GOPR12Fisher, "PR12MGOFisher.csv")
write_csv(GOGR12KS, "GR12MGOKS.csv")
write_csv(GOGR12Fisher, "GR12MGOFisher.csv")