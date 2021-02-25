library(tidyverse)
library(DESeq2)
library(BrainInABlender)

MergedReads.csv <- read_csv("MergedReads.csv")

BMOutput <- read_csv("FTDBMOutput.csv")

F210ISupport <- read_csv("FTDSupport.csv")[,1:2]

NormaliseReads <- function (Reads, SupportFile, BioMOutput){
  #Change sample names so they're more human readable and it's immediately obvious which is WT and which is mutant
  Reads <- column_to_rownames(Reads, "EnsemblID")
  SupportFile$condition_A<- ifelse(SupportFile[,2] == "CTL", "CTL","Mut")
  row.names(SupportFile) <- SupportFile$sample
  SupportFile$sample <- NULL
  SupportFile$Condition <- factor(SupportFile$condition_A, levels=c("CTL","Mut"))
  
  #Get Normalised Reads from deseq
  formula0 = as.formula("~ 1") 
  CDS <- DESeqDataSetFromMatrix(countData = Reads, colData = SupportFile, design = as.formula("~ Condition"))
  CDS <- estimateSizeFactors(CDS)
  
  #Extract normalised reads
  normalised_counts <- as.data.frame(counts(CDS, normalized=TRUE)) %>% rownames_to_column("EnsemblID") %>%
    left_join(BMOutput, by=c("EnsemblID" = "ensembl_gene_id")) %>% drop_na("external_gene_name") %>% dplyr::select(external_gene_name, A075:P94)
  return(normalised_counts)
}

NormalisedUnfilteredReads <- NormaliseReads(MergedReads, Support)
FilteredReads <- filter(MergedReads, grepl("ENSG", EnsemblID))
NormalisedFilteredReads <- NormaliseReads(FilteredReads, Support)

BrainType <- Sir_UnMixALot(userInput=NormalisedUnfilteredReads, dataColumns=c(2:ncol(NormalisedUnfilteredReads)), geneColumn=1, species="human")$AveragePrimary_CellTypeIndex
BrainType <- as.data.frame(BrainType) %>% rownames_to_column("CellType")

write_csv(as.data.frame(BrainType), "FTDBrainBIAB.csv")
