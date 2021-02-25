library(tidyverse)
library(EnrichmentBrowser)
library(EnhancedVolcano)
library(patchwork)

#Import Reads
NormalisedReads <- read_csv("NormalisedReads.csv")

#Create support df
F210IAdultSupport <- read_csv("F210IAdultSupport.csv")


MakePCA <- function(Reads, Support){
  SelectedReads <- dplyr::select(Reads, ensemblID, Support$sample) %>% filter(grepl("ENSMUSG", ensemblID))
  SelectedReads$RowSum <- SelectedReads %>% rowwise() %>% dplyr::select(-ensemblID) %>% rowSums()
  #Drop all genes with fewer than 10 reads overall
  PCAReads <- SelectedReads %>% filter(RowSum < 10) %>% dplyr::select(-RowSum) %>% column_to_rownames("ensemblID") 
  #run principal component analysis 
  #Scaling is on because otherwise one outlier distorts the whole graph
  PCA <- prcomp(t(PCAReads), scale = TRUE)
  PCAData <- as.data.frame(PCA$x) %>% 
    rownames_to_column("sample") %>%
    left_join(Support, by="sample") 
  
  ggplot(PCAData, aes(x=PC1,y=PC2), label=sample) + geom_point(aes(color=GROUP), size=4) +
    theme_bw(base_size=32) + 
    theme(legend.position="top") + geom_text(aes(label=sample),hjust=0, vjust=0)
  
}

svg("NormalisedReadPCA.svg")
MakePCA(NormalisedReads, F210IAdultSupport)
dev.off()
