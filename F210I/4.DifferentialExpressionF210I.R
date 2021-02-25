library(tidyverse)
library(EnrichmentBrowser)
library(EnhancedVolcano)
library(patchwork)

MakeSignedZscore <- function(DEMatrix){
  DEMatrix$signedZ <- ifelse(test = DEMatrix$FC > 0,
                             yes = qnorm(1 - (DEMatrix$PVAL / 2)),
                             no = qnorm(DEMatrix$PVAL / 2) )
  # high pass - otherwise very low P is set to Inf
  DEMatrix$signedZ[DEMatrix$signedZ > 6] <- 6
  DEMatrix$signedZ[DEMatrix$signedZ < -6] <- -6
  return(DEMatrix)
}

#A function to make a merged matrix of all 3 methods of differential expression
DifferentialExpressionMerge <- function(RawReads, Support){
  #create a summarized experiment with the reads and support file
  data.se <- SummarizedExperiment(list(counts=as.matrix(RawReads)), colData = Support)
  
  #Run DESeq, then change the rownames to a column to make them easier to manage
  deseq = as.data.frame(rowData(deAna(data.se, de.method = "DESeq2"))) %>%
    rownames_to_column("gene_id")
  #Get ZScore
  deseqZ <- MakeSignedZscore(deseq)
  return(deseqZ)
}

#import raw datasets
F210IAdultReads <- read_csv("F210IRawReads.csv")[-c(1:5),] %>%
  column_to_rownames("EnsemblID")

F210IAdultSupport <- read_csv("F210IAdultSupport.csv")

F210IAdultReadsDE <- DifferentialExpressionMerge(F210IAdultReads, F210IAdultSupport)

write_csv(F210IAdultReadsDE, "F210IAdultDE.csv")


