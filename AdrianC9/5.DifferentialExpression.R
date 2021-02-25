#Adrian C9 Differential expression
library(tidyverse)
library(EnrichmentBrowser)
library(biomaRt)


Reads <- read_csv("AdrianC9GeneReadsAll.csv")
Support <- read_csv("AdrianC9SupportUpdated.csv")
BMOutput <- read_csv("BioMarTMouseWEntrez.csv")
#Support <- Support %>% filter(Sample != "NM8084_800973", Sample != "NM8093_793422")

SupportPR6M <- dplyr::select(Support, Sample, condition_B, type_3) %>% drop_na()
SupportGR6M <- dplyr::select(Support, Sample, condition_C, type_3) %>% drop_na()
SupportPR12M <- dplyr::select(Support, Sample, condition_D, type_3) %>% drop_na()
SupportGR12M <- dplyr::select(Support, Sample, condition_E, type_3) %>% drop_na()

#create function for getting DE
MakeDE <- function(SupportData, ReadData){
  #Set the condition to be called Group
  colnames(SupportData)[2] <- "GROUP"
  #Drop all files that don't have ensemblIDs and convert to rownames
  FiltReads <- dplyr::select(ReadData, EnsemblID, SupportData$Sample) %>% filter(grepl("ENSMUSG", EnsemblID)) 
  FiltReads <- FiltReads %>% column_to_rownames("EnsemblID")
  
  #use enrichment browser's summarised experiment package to run differential expression after converting 
  #Group column to a binary value
  SEData <- SummarizedExperiment(list(counts=as.matrix(FiltReads)), colData = SupportData)
  SEData$GROUP = ifelse(SupportData$GROUP == "WT", 0, 1)
  SEData <- deAna(SEData, de.method = "DESeq2")
  
  Results = as.data.frame(rowData(SEData)) %>%
    rownames_to_column("gene_id")
  
  return(Results)
}

#Run DE
SEDataPR6M <- MakeDE(SupportPR6M, Reads)
SEDataGR6M <- MakeDE(SupportGR6M, Reads)
SEDataPR12M <- MakeDE(SupportPR12M, Reads)
SEDataGR12M <- MakeDE(SupportGR12M, Reads)

#Merge with the biomart output to get gene name
deseqPR6MJ <- left_join(SEDataPR6M, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)
deseqGR6MJ <- left_join(SEDataGR6M, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)
deseqPR12MJ <- left_join(SEDataPR12M, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)
deseqGR12MJ <- left_join(SEDataGR12M, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)

#Write results
write_csv(deseqPR6MJ, "DifferentialExpressionPR6MNoType.csv")
write_csv(deseqGR6MJ, "DifferentialExpressionGR6MNoType.csv")
write_csv(deseqPR12MJ, "DifferentialExpressionPR12MNoType.csv")
write_csv(deseqGR12MJ, "DifferentialExpressionGR12MNoType.csv")

################################################################################################################
#Run DE on samples of the same age and samples of the same phenotype

Support <- dplyr::select(Support, Sample, condition_A, type_1, type_2, type_3) %>% drop_na()

SixMonthSupport <- filter(Support, type_3 == "6M")
SixMonthSupport$condition_A <- ifelse(SixMonthSupport$condition_A == "WT", "WT", "Mut")
TwelveMonthSupport <- filter(Support, type_3 == "12M")
TwelveMonthSupport$condition_A <- ifelse(TwelveMonthSupport$condition_A == "WT", "WT", "Mut")


GRSupport <- filter(Support, condition_A != "(PR)400")
PRSupport <- filter(Support, condition_A != "(GR)400")


SixMonthDe1 <- MakeDE(SixMonthSupport, Reads)
TwelveMonthDe1 <- MakeDE(TwelveMonthSupport, Reads)

GRDe1 <- MakeDE(GRSupport, Reads)
PRDe1 <- MakeDE(PRSupport, Reads)

deseq6M <- left_join(SixMonthDe1, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)
deseq12M <- left_join(TwelveMonthDe1, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)
deseqGR <- left_join(GRDe1, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)
deseqPR <- left_join(PRDe1, BMOutput, by = c("gene_id"="ensembl_gene_id")) %>% arrange(ADJ.PVAL)

write_csv(deseq6M, "DifferentialExpression6MNoType.csv")
write_csv(deseq12M, "DifferentialExpression12MNoType.csv")
write_csv(deseqGR, "DifferentialExpressionGRNoType.csv")
write_csv(deseqPR, "DifferentialExpressionPRNoType.csv")

