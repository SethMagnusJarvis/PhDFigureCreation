library(tidyverse)
library(EnrichmentBrowser)
library(EnhancedVolcano)
library(patchwork)

DifferentialExpressionMerge <- function(RawReads, Support){
  #create a summarized experiment with the reads and support file
  data.se <- SummarizedExperiment(list(counts=as.matrix(RawReads)), colData = Support)
  Normalised <- normalize(
    data.se,
    data.type = "rseq",
  )
  return(as.data.frame(assays(Normalised)$counts))
}

#import raw datasets
F210IAdultReads <- read_csv("F210IRawReads.csv")[-c(1:5),] %>%
  column_to_rownames("EnsemblID")


#Create support files for KO and Quant (Summarised experiment needs 0s and 1s)
F210IAdultSupport<-data.frame(SAMPLE = colnames(F210IAdultReads))
F210IAdultSupport$Condition[1:7] <- "HET"
F210IAdultSupport$Condition[8:12] <- "control"
F210IAdultSupport$Condition[1:7] <- "HET"
F210IAdultSupport$GROUP = ifelse(F210IAdultSupport$Condition == "control", 0, 1)
F210IAdultSupport <- select(F210IAdultSupport, SAMPLE, GROUP)

RPKM <- read_csv("rpkm_values.csv") %>% dplyr::select(external_gene_id, ensemblID )

NormalisedReads <- DifferentialExpressionMerge(F210IAdultReads, F210IAdultSupport) %>% rownames_to_column("ensemblID") %>%
  left_join(RPKM, by="ensemblID") %>% select(external_gene_id, everything())

write_csv(NormalisedReads, "NormalisedReads.csv")