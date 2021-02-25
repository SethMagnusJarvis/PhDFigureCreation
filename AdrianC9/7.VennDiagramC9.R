library(tidyverse)
library(ggplot2)
library(limma)

PR12MDE <- read_csv("DEPR12M.csv")
GR12MDE <- read_csv("DEGR12M.csv")

DE1 <- PR12MDE
DE2 <- GR12MDE

ListCompare <- function(DE1, DE2){
  DE1SigList <- DE1 %>% filter(ADJ.PVAL < 0.05) %>% select(gene_id)
  DE2SigList <- DE2 %>% filter(ADJ.PVAL < 0.05) %>% select(gene_id)
  ids = sort(unique(c(as.character(unlist(DE1SigList)),  as.character(unlist(DE2SigList) ))))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% unlist(DE1SigList)
    counts[i, 2] = ids[i] %in% unlist(DE2SigList)
  }
  colnames(counts) = c("PR", "GR")
  row.names(counts) = ids
  
  return(as.data.frame(counts))
}


PRGRComp <- ListCompare(PR12MDE, GR12MDE)

vennDiagram(PRGRComp) + title("Significant PR and GR gene Venn Diagram")

PROnly <- PRGRComp %>% rownames_to_column("gene_id") %>% filter(PR == 1 & GR == 0) %>% left_join(PR12MDE) %>%
  select(gene_id, external_gene_name, FC:ADJ.PVAL) %>% arrange(ADJ.PVAL)
GROnly <- PRGRComp %>% rownames_to_column("gene_id") %>% filter(PR == 0 & GR == 1) %>% left_join(PR12MDE) %>%
  select(gene_id, external_gene_name, FC:ADJ.PVAL) %>% arrange(ADJ.PVAL)
Overlap <- PRGRComp %>% rownames_to_column("gene_id") %>% filter(PR == 1 & GR == 1) %>% left_join(PR12MDE) %>%
  select(gene_id, external_gene_name, FC:ADJ.PVAL) %>% arrange(ADJ.PVAL)

write_csv(PROnly, "PROnlySigGenes.csv")
write_csv(GROnly, "GROnlySigGenes.csv")
write_csv(Overlap, "OverlapOnlySigGenes.csv")
