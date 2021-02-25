library(tidyverse)
library(limma)

AllJoinedDEBIAB <- read_csv("CtlVsMutDEBothCatAge.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
CTLC9JoinedDEBIAB <- read_csv("CtlVsC9DEBothCatAge.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
CTLTAUJoinedDEBIAB <- read_csv("CtlVsTAUDEBothCatAge.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
C9TAUJoinedDEBIAB <- read_csv("C9VsTAUDEBothCatAge.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))

AllJoinedDEXCell <- read_csv("CTLVsMUTDEXCellNeuron.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
CTLC9JoinedDEXCell <- read_csv("CTLVsC9DEXCellNeuron.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
CTLTAUJoinedDEXCell <- read_csv("CTLVsTAUDEXCellNeuron.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
C9TAUJoinedDEXCell <- read_csv("C9VsTAUDEXCellNeuron.csv") %>% select(external_gene_name, EnsemblID, ends_with(".Join")) %>%
  rename_at(.vars = vars(ends_with(".Join")),
            .funs = funs(sub("[.]Join$", "", .)))
            
          
ListCompare <- function(DEBIAB, DEXCell){
  DE1SigList <- DEBIAB %>% filter(padj < 0.05) %>% 
    dplyr::select(external_gene_name)
  DE2SigList <- DEXCell %>% filter(padj < 0.05) %>% 
    dplyr::select(external_gene_name)
  ids = sort(unique(c(as.character(unlist(DE1SigList)),  as.character(unlist(DE2SigList) ))))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% unlist(DE1SigList)
    counts[i, 2] = ids[i] %in% unlist(DE2SigList)
  }
  colnames(counts) = c("BIAB", "XCell")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

AllCompared <- ListCompare(AllJoinedDEBIAB, AllJoinedDEXCell)
CTLC9Compared <- ListCompare(CTLC9JoinedDEBIAB, CTLC9JoinedDEXCell)
CTLTAUCompared <- ListCompare(CTLTAUJoinedDEBIAB, CTLTAUJoinedDEXCell)
C9TAUCompared <- ListCompare(C9TAUJoinedDEBIAB, C9TAUJoinedDEXCell)

par(mfrow=c(2,2))
vennDiagram(AllCompared) + title("Control Vs Mutant Venn")
vennDiagram(CTLC9Compared) + title("Control vs C9 Venn")
vennDiagram(CTLTAUCompared) + title("Control vs TAU Venn")
vennDiagram(C9TAUCompared) + title("C9 vs TAU Venn")
