library(tidyverse)
library(limma)
#Venn Diagrams
AllJoinedDE <- read_csv("CtlVsMutDEBothCatAge.csv")
CTLC9JoinedDE <- read_csv("CtlVsC9DEBothCatAge.csv")
CTLTAUJoinedDE <- read_csv("CtlVsTAUDEBothCatAge.csv")
C9TAUJoinedDE <- read_csv("C9VsTAUDEBothCatAge.csv")

ListCompare <- function(DEJoin){
  DE1SigList <- DEJoin %>% dplyr::select(external_gene_name, contains(".Join")) %>% filter(padj.Join < 0.05) %>% 
    dplyr::select(external_gene_name)
  DE2SigList <- DEJoin %>% dplyr::select(external_gene_name, contains(".Solo")) %>% filter(padj.Solo < 0.05) %>% 
    dplyr::select(external_gene_name)
  ids = sort(unique(c(as.character(unlist(DE1SigList)),  as.character(unlist(DE2SigList) ))))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% unlist(DE1SigList)
    counts[i, 2] = ids[i] %in% unlist(DE2SigList)
  }
  colnames(counts) = c("Join", "Solo")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

AllCompared <- ListCompare(AllJoinedDE)
CTLC9Compared <- ListCompare(CTLC9JoinedDE)
CTLTAUCompared <- ListCompare(CTLTAUJoinedDE)
C9TAUCompared <- ListCompare(C9TAUJoinedDE)

par(mfrow=c(2,2))
vennDiagram(AllCompared) + title("Control Vs Mutant Venn")
vennDiagram(CTLC9Compared) + title("Control vs C9 Venn")
vennDiagram(CTLTAUCompared) + title("Control vs TAU Venn")
vennDiagram(C9TAUCompared) + title("C9 vs TAU Venn")
