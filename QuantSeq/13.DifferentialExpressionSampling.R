library(tidyverse)
library(ggplot2)
library(EnrichmentBrowser)

MakeSignedZscore <- function(DEMatrix){
  DEMatrix$signedZ <- ifelse(test = DEMatrix$FC > 0,
                             yes = qnorm(1 - (DEMatrix$PVAL / 2)),
                             no = qnorm(DEMatrix$PVAL / 2) )
  # high pass - otherwise very low P is set to Inf
  DEMatrix$signedZ[DEMatrix$signedZ > 6] <- 6
  DEMatrix$signedZ[DEMatrix$signedZ < -6] <- -6
  return(DEMatrix)
}

DifferentialExpressionMerge <- function(RawReads, Support){
  #create a summarized experiment with the reads and support file
  data.se <- SummarizedExperiment(list(counts=as.matrix(RawReads)), colData = Support)
  
  #Run DESeq, then change the rownames to a column to make them easier to manage
  deseq = as.data.frame(rowData(deAna(data.se, de.method = "DESeq2"))) %>%
    rownames_to_column("gene_id")
  #Get ZScore
  deseqZ <- MakeSignedZscore(deseq)
  #Add a deseq ending to it to make it easier to merge later
  colnames(deseqZ) <- paste(colnames(deseqZ), "DESeq", sep = ".")
  #The same but in edgeR
  edgeR = as.data.frame(rowData(deAna(data.se, de.method = "edgeR")))%>%
    rownames_to_column("gene_id")
  edgeRZ <- MakeSignedZscore(edgeR)
  colnames(edgeRZ) <- paste(colnames(edgeRZ), "edgeR", sep = ".")
  #The same in limma
  limma = as.data.frame(rowData(deAna(data.se, de.method = "limma")))%>%
    rownames_to_column("gene_id")
  limmaZ <- MakeSignedZscore(limma)
  colnames(limmaZ) <- paste(colnames(limmaZ), "limma", sep = ".")
  
  #Merge all of the differential expression dataframes together by their gene ID
  Merge <- full_join(as.data.frame(deseqZ), as.data.frame(edgeRZ), by = c("gene_id.DESeq" = "gene_id.edgeR")) %>%
    full_join(limmaZ, by = c("gene_id.DESeq" = "gene_id.limma"))
  #rename the geneID column to just geneID, then return the whole dataframe
  colnames(Merge)[1] <- "gene_id"
  return(Merge)
}

KOSampledOrigin <- read_csv("KOTotalRawSampledMk2.csv") %>% 
  select(gene_id, k1WT,  k2WT,  k3WT, k4WT, k1KO, k2KO, k3KO, k4KO) %>% column_to_rownames("gene_id") %>% drop_na()
d14SampledOrigin <- read_csv("d14TotalRawSampledMk2.csv") %>% 
  select(gene_id, d1WT, d2WT, d3WT, d4WT, d1HOM, d2HOM, d3HOM, d4HOM) %>% column_to_rownames("gene_id") %>% drop_na()

KOSampled3mil <- read_csv("KOTotalRawSampled3MilMk1.csv") %>% 
  select(gene_id, k1WT,  k2WT,  k3WT, k4WT, k1KO, k2KO, k3KO, k4KO) %>% column_to_rownames("gene_id") %>% drop_na()
d14Sampled3mil <- read_csv("d14TotalRawSampled3MilMk1.csv") %>% 
  select(gene_id, d1WT, d2WT, d3WT, d4WT, d1HOM, d2HOM, d3HOM, d4HOM) %>% column_to_rownames("gene_id") %>% drop_na()

KOSampled10mil <- read_csv("KOTotalRawSampled10MilMK1.csv") %>% 
  select(gene_id, k1WT,  k2WT,  k3WT, k4WT, k1KO, k2KO, k3KO, k4KO) %>% column_to_rownames("gene_id") %>% drop_na()
d14Sampled10mil <- read_csv("d14TotalRawSampled10MilMK1.csv") %>% 
  select(gene_id, d1WT, d2WT, d3WT, d4WT, d1HOM, d2HOM, d3HOM, d4HOM) %>% column_to_rownames("gene_id") %>% drop_na()

KOSampled15mil <- read_csv("KOTotalRawSampled15MilMK1.csv") %>% 
  select(gene_id, k1WT,  k2WT,  k3WT, k4WT, k1KO, k2KO, k3KO, k4KO) %>% column_to_rownames("gene_id") %>% drop_na()
d14Sampled15mil <- read_csv("d14TotalRawSampled15MilMK1.csv") %>% 
  select(gene_id, d1WT, d2WT, d3WT, d4WT, d1HOM, d2HOM, d3HOM, d4HOM) %>% column_to_rownames("gene_id") %>% drop_na()

KOSampled20mil <- read_csv("KOTotalRawSampled20MilMK1.csv") %>% 
  select(gene_id, k1WT,  k2WT,  k3WT, k4WT, k1KO, k2KO, k3KO, k4KO) %>% column_to_rownames("gene_id") %>% drop_na()
d14Sampled20mil <- read_csv("d14TotalRawSampled20MilMK1.csv") %>% 
  select(gene_id, d1WT, d2WT, d3WT, d4WT, d1HOM, d2HOM, d3HOM, d4HOM) %>% column_to_rownames("gene_id") %>% drop_na()

KOQuantRawV1 <- read_csv("KOQuantRawV1.csv") %>% column_to_rownames("gene_id") %>% drop_na()
d14QuantRawV1 <- read_csv("d14QuantRawV1.csv") %>% column_to_rownames("gene_id") %>% drop_na()
KOQuantRawV2 <- read_csv("KOQuantRawV2.csv") %>% column_to_rownames("gene_id") %>% drop_na()
d14QuantRawV2 <- read_csv("d14QuantRawV2.csv") %>% column_to_rownames("gene_id") %>% drop_na()

KOSupport<-data.frame(SAMPLE = colnames(KOSampledOrigin))
KOSupport$Condition[1:4] <- "control"
KOSupport$Condition[5:8] <- "KO"
KOSupport$GROUP = ifelse(KOSupport$Condition == "control", 0, 1)
KOSupport <- select(KOSupport, SAMPLE, GROUP)

d14Support<-data.frame(SAMPLE = colnames(d14SampledOrigin))
d14Support$Condition[1:4] <- "control"
d14Support$Condition[5:8] <- "HOM"
d14Support$GROUP = ifelse(d14Support$Condition == "control", 0, 1)
d14Support <- select(d14Support, SAMPLE, GROUP)

KOQuantV1De <- DifferentialExpressionMerge(KOQuantRawV1, KOSupport)
d14QuantV1De <- DifferentialExpressionMerge(d14QuantRawV1, d14Support)
KOQuantV2De <- DifferentialExpressionMerge(KOQuantRawV2, KOSupport)
d14QuantV2De <- DifferentialExpressionMerge(d14QuantRawV2, d14Support)

KOOriginTotalDe <- DifferentialExpressionMerge(KOSampledOrigin, KOSupport)
d14OriginTotalDe <- DifferentialExpressionMerge(d14SampledOrigin, d14Support)
KO3MilTotalDe <- DifferentialExpressionMerge(KOSampled3mil, KOSupport)
d143MilTotalDe <- DifferentialExpressionMerge(d14Sampled3mil, d14Support)
KO10MilTotalDe <- DifferentialExpressionMerge(KOSampled10mil, KOSupport)
d1410MilTotalDe <- DifferentialExpressionMerge(d14Sampled10mil, d14Support)
KO15MilTotalDe <- DifferentialExpressionMerge(KOSampled15mil, KOSupport)
d1415MilTotalDe <- DifferentialExpressionMerge(d14Sampled15mil, d14Support)
KO20MilTotalDe <- DifferentialExpressionMerge(KOSampled20mil, KOSupport)
d1420MilTotalDe <- DifferentialExpressionMerge(d14Sampled20mil, d14Support)

write_csv(KOQuantV1De, "KOQuantV1De.csv")
write_csv(KOQuantV2De, "KOQuantV2De.csv")
write_csv(d14QuantV1De, "d14QuantV1De.csv")
write_csv(d14QuantV2De, "d14QuantV2De.csv")

write_csv(KOOriginTotalDe, "KOOriginTotalDe.csv")
write_csv(d14OriginTotalDe, "d14OriginTotalDe.csv")
write_csv(KO3MilTotalDe, "KO3MilTotalDe.csv")
write_csv(d143MilTotalDe, "d143MilTotalDe.csv")
write_csv(KO10MilTotalDe, "KO10MilTotalDe.csv")
write_csv(d1410MilTotalDe, "d1410MilTotalDe.csv")
write_csv(KO15MilTotalDe, "KO15MilTotalDe.csv")
write_csv(d1415MilTotalDe, "d1415MilTotalDe.csv")
write_csv(KO20MilTotalDe, "KO20MilTotalDe.csv")
write_csv(d1420MilTotalDe, "d1420MilTotalDe.csv")


DifferentialExpressionMergeNoFilt <- function(RawReads, Support){
  #create a summarized experiment with the reads and support file
  data.se <- SummarizedExperiment(list(counts=as.matrix(RawReads)), colData = Support)
  
  #Run DESeq, then change the rownames to a column to make them easier to manage
  deseq = as.data.frame(rowData(deAna(data.se, de.method = "DESeq2", filter.by.expr = FALSE))) %>%
    rownames_to_column("gene_id")
  #Get ZScore
  deseqZ <- MakeSignedZscore(deseq)
  #Add a deseq ending to it to make it easier to merge later
  colnames(deseqZ) <- paste(colnames(deseqZ), "DESeq", sep = ".")
  #The same but in edgeR
  edgeR = as.data.frame(rowData(deAna(data.se, de.method = "edgeR", filter.by.expr = FALSE)))%>%
    rownames_to_column("gene_id")
  edgeRZ <- MakeSignedZscore(edgeR)
  colnames(edgeRZ) <- paste(colnames(edgeRZ), "edgeR", sep = ".")
  #The same in limma
  limma = as.data.frame(rowData(deAna(data.se, de.method = "limma", filter.by.expr = FALSE)))%>%
    rownames_to_column("gene_id")
  limmaZ <- MakeSignedZscore(limma)
  colnames(limmaZ) <- paste(colnames(limmaZ), "limma", sep = ".")
  
  #Merge all of the differential expression dataframes together by their gene ID
  Merge <- full_join(as.data.frame(deseqZ), as.data.frame(edgeRZ), by = c("gene_id.DESeq" = "gene_id.edgeR")) %>%
    full_join(limmaZ, by = c("gene_id.DESeq" = "gene_id.limma"))
  #rename the geneID column to just geneID, then return the whole dataframe
  colnames(Merge)[1] <- "gene_id"
  return(Merge)
}

KOQuantV1DeNoFilt <- DifferentialExpressionMergeNoFilt(KOQuantRawV1, KOSupport)
d14QuantV1DeNoFilt <- DifferentialExpressionMergeNoFilt(d14QuantRawV1, d14Support)
KOQuantV2DeNoFilt <- DifferentialExpressionMergeNoFilt(KOQuantRawV2, KOSupport)
d14QuantV2DeNoFilt <- DifferentialExpressionMergeNoFilt(d14QuantRawV2, d14Support)

KOOriginTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(KOSampledOrigin, KOSupport)
d14OriginTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(d14SampledOrigin, d14Support)
KO3MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(KOSampled3mil, KOSupport)
d143MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(d14Sampled3mil, d14Support)
KO10MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(KOSampled10mil, KOSupport)
d1410MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(d14Sampled10mil, d14Support)
KO15MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(KOSampled15mil, KOSupport)
d1415MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(d14Sampled15mil, d14Support)
KO20MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(KOSampled20mil, KOSupport)
d1420MilTotalDeNoFilt <- DifferentialExpressionMergeNoFilt(d14Sampled20mil, d14Support)

write_csv(KOQuantV1DeNoFilt, "KOQuantV1DeNoFilt.csv")
write_csv(KOQuantV2DeNoFilt, "KOQuantV2DeNoFilt.csv")
write_csv(d14QuantV1DeNoFilt, "d14QuantV1DeNoFilt.csv")
write_csv(d14QuantV2DeNoFilt, "d14QuantV2DeNoFilt.csv")

write_csv(KOOriginTotalDeNoFilt, "KOOriginTotalDeNoFilt.csv")
write_csv(d14OriginTotalDeNoFilt, "d14OriginTotalDeNoFilt.csv")
write_csv(KO3MilTotalDeNoFilt, "KO3MilTotalDeNoFilt.csv")
write_csv(d143MilTotalDeNoFilt, "d143MilTotalDeNoFilt.csv")
write_csv(KO10MilTotalDeNoFilt, "KO10MilTotalDeNoFilt.csv")
write_csv(d1410MilTotalDeNoFilt, "d1410MilTotalDeNoFilt.csv")
write_csv(KO15MilTotalDeNoFilt, "KO15MilTotalDeNoFilt.csv")
write_csv(d1415MilTotalDeNoFilt, "d1415MilTotalDeNoFilt.csv")
write_csv(KO20MilTotalDeNoFilt, "KO20MilTotalDeNoFilt.csv")
write_csv(d1420MilTotalDeNoFilt, "d1420MilTotalDeNoFilt.csv")

