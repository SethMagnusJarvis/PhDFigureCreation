library(tidyverse)
library(ggplot2)

load("D:/GoogleDrive/Work/UCL/F210IAnalysis/sgseq/F210I_adult_sc_res_clean_novel.RData")

SplicingReads <- res.clean[,11:22]
colnames(SplicingReads) <- c("WT1", "WT2", "WT3", "WT4", "WT5", "Het1", "Het2", "Het3", "Het4", "Het5", "Het6", "Het7")

F210IAdultSupport <- read_csv("F210IAdultSupport.csv") %>% mutate(type = case_when(GROUP == 0 ~ "WT",
                                                                                   GROUP == 1 ~ "HET"))
Reads <- SplicingReads
Support <- F210IAdultSupport

MakePCA <- function(Reads, Support){
  SelectedReads <- Reads
  SelectedReads$RowSum <- Reads %>% rowwise() %>% rowSums()
  #Drop all genes with fewer than 10 reads overall
  PCAReads <- SelectedReads %>% filter(!RowSum == 0 ) %>% dplyr::select(-RowSum)
  PCA <- prcomp(t(PCAReads))
  Support <- dplyr::select(Support, sample, type)
  
  PCAData <- as.data.frame(PCA$x) %>% 
    rownames_to_column("sample") %>%
    left_join(Support, by="sample") 
  
  
  ggplot(PCAData, aes(x=PC1,y=PC2), label=sample) + geom_point(aes(color=type), size=4) +
    theme_bw(base_size=32) + 
    theme(legend.position="top") + geom_text(aes(label=sample),hjust=0, vjust=0)
}

svg("SplicingPCA.svg")
MakePCA(SplicingReads, F210IAdultSupport)
dev.off()
