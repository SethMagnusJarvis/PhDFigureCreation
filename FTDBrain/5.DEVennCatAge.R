library(tidyverse)
library(limma)

#Venn Diagrams
AllJoinedDE <- read_csv("CtlVsMutDEBothCatAge.csv")
CTLC9JoinedDE <- read_csv("CtlVsC9DEBothCatAge.csv")
CTLTAUJoinedDE <- read_csv("CtlVsTAUDEBothCatAge.csv")
C9TAUJoinedDE <- read_csv("C9VsTAUDEBothCatAge.csv")

#Function to create a list with a list of all genes significant in either of two datasets, and which datasets they
#are significant in for use by limma's venn function
ListCompare <- function(DEJoin1, DEJoin2){
  DE1SigList <- DEJoin1 %>% dplyr::select(external_gene_name, contains(".Solo")) %>% filter(padj.Solo < 0.05) %>% 
    dplyr::select(external_gene_name)
  DE2SigList <- DEJoin2 %>% dplyr::select(external_gene_name, contains(".Solo")) %>% filter(padj.Solo < 0.05) %>% 
    dplyr::select(external_gene_name)
  ids = sort(unique(c(as.character(unlist(DE1SigList)),  as.character(unlist(DE2SigList) ))))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% unlist(DE1SigList)
    counts[i, 2] = ids[i] %in% unlist(DE2SigList)
  }
  colnames(counts) = c("AllSamples", "Seperate")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

AllVsC9 <- ListCompare(AllJoinedDE, CTLC9JoinedDE)
AllVsTAU <- ListCompare(AllJoinedDE, CTLTAUJoinedDE)
AllVsC9TAU <- ListCompare(AllJoinedDE, C9TAUJoinedDE)






ListCompare2 <- function(DEJoin1, DEJoin2){
  DE1SigList <- DEJoin1 %>% dplyr::select(external_gene_name, contains(".Solo")) %>% filter(padj.Solo < 0.05) %>% 
    dplyr::select(external_gene_name)
  DE2SigList <- DEJoin2 %>% dplyr::select(external_gene_name, contains(".Solo")) %>% filter(padj.Solo < 0.05) %>% 
    dplyr::select(external_gene_name)
  ids = sort(unique(c(as.character(unlist(DE1SigList)),  as.character(unlist(DE2SigList) ))))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% unlist(DE1SigList)
    counts[i, 2] = ids[i] %in% unlist(DE2SigList)
  }
  colnames(counts) = c("C9", "TAU")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

C9VsTAU <- ListCompare2(CTLC9JoinedDE, CTLTAUJoinedDE)

par(mfrow=c(2,2))
vennDiagram(C9VsTAU) + title("DE in C9 vs in TAU")
vennDiagram(AllVsC9) + title("DE in both vs in C9")
vennDiagram(AllVsTAU) + title("DE in both vs in TAU")
vennDiagram(AllVsC9TAU) + title("DE in both vs in C9vTAU")

svg("VennSoloC9VsTau.svg")
vennDiagram(C9VsTAU) + title("DE in C9 vs in TAU")
dev.off()

svg("VennSoloBothVsC9.svg")
vennDiagram(AllVsC9) + title("DE in both vs in C9")
dev.off()

svg("VennSoloBothVsTAU.svg")
vennDiagram(AllVsTAU) + title("DE in both vs in TAU")
dev.off()

svg("VennSoloBothVsC9vTAU.svg")
vennDiagram(AllVsC9TAU) + title("DE in both vs in C9vTAU")
dev.off()