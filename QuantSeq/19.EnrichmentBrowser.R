library(tidyverse)
library(EnrichmentBrowser)
library(ggplot2)
library(patchwork)

########################################################################################################################################
###################################################Setup functions for later############################################################
########################################################################################################################################

#A function to make a singned Z score whichi is derived from the pvalue
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

########################################################################################################################################
########################################################Import data#####################################################################
########################################################################################################################################

#import raw datasets
d14TotalRaw <- read_csv("d14TotalRawUnsampled.csv") %>%
  column_to_rownames("gene_id")
KOTotalRaw <- read_csv("KOTotalRawUnsampled.csv")%>%
  column_to_rownames("gene_id")

d14QuantRawV1 <- read_csv("d14QuantRawV1.csv")%>%
  column_to_rownames("gene_id")
KOQuantRawV1 <- read_csv("KOQuantRawV1.csv")%>%
  column_to_rownames("gene_id")

d14QuantRawV2 <- read_csv("d14QuantRawV2.csv")%>%
  column_to_rownames("gene_id")
KOQuantRawV2 <- read_csv("KOQuantRawV2.csv")%>%
  column_to_rownames("gene_id")

#Create support files for KO and Quant (Summarised experiment needs 0s and 1s)
KOSupport<-data.frame(SAMPLE = colnames(KOTotalRaw))
KOSupport$Condition[1:4] <- "control"
KOSupport$Condition[5:8] <- "KO"
KOSupport$GROUP = ifelse(KOSupport$Condition == "control", 0, 1)
KOSupport <- select(KOSupport, SAMPLE, GROUP)

d14Support<-data.frame(SAMPLE = colnames(d14TotalRaw))
d14Support$Condition[1:4] <- "control"
d14Support$Condition[5:8] <- "HOM"
d14Support$GROUP = ifelse(d14Support$Condition == "control", 0, 1)
d14Support <- select(d14Support, SAMPLE, GROUP)

########################################################################################################################################
################################################Run Differential expression#############################################################
########################################################################################################################################

#Run differential expression in all datasets
d14QuantV1De <- DifferentialExpressionMerge(d14QuantRawV1, d14Support)
d14QuantV2De <- DifferentialExpressionMerge(d14QuantRawV2, d14Support)
KOQuantV1De <- DifferentialExpressionMerge(KOQuantRawV1, KOSupport)
KOQuantV2De <- DifferentialExpressionMerge(KOQuantRawV2, KOSupport)

d14TotalDe <- DifferentialExpressionMerge(d14TotalRaw, d14Support)
KOTotalDe <- DifferentialExpressionMerge(KOTotalRaw, KOSupport)


########################################################################################################################################
###############################################Only look at highly expressed genes###########################################
########################################################################################################################################

#Import base Means for each dataset
d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")[,1:2]
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")[,1:2]
d14Quant <- read_csv("DESeqD14QuantSeq.csv")[,1:2]
KOQuant <- read_csv("DESeqKOQuantSeq.csv")[,1:2]

#Merge Basemean data with the DE datasets
d14QuantV1DeBM <- left_join(d14QuantV1De, d14Quant, by=c("gene_id" = "EnsemblID"))
KOQuantV1DeBM <- left_join(KOQuantV1De, KOQuant, by=c("gene_id" = "EnsemblID"))
d14TotalDeBM <- left_join(d14TotalDe, d14Unsampled, by=c("gene_id" = "EnsemblID"))
KOTotalDeBM <- left_join(KOTotalDe, KOUnsampled, by=c("gene_id" = "EnsemblID"))

#Merge above into one dataset
d14DMergeV1 <- left_join(d14QuantV1DeBM, d14TotalDeBM, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))
KODMergeV1 <- left_join(KOQuantV1DeBM, KOTotalDeBM, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))

#Sort the datasets by basemean in each data-type
#D14
d14DMergeV1Quant <- arrange(d14DMergeV1,desc(baseMean.QuantSeq))
d14DMergeV1QuantTop <- d14DMergeV1Quant[1:round(nrow(d14DMergeV1Quant)/4),]
d14DMergeV1Total <- arrange(d14DMergeV1,desc(baseMean.TotalSeq))
d14DMergeV1TotalTop <- d14DMergeV1Total[1:round(nrow(d14DMergeV1Total)/4),]

#KO
KODMergeV1Quant <- arrange(KODMergeV1,desc(baseMean.QuantSeq))
KODMergeV1QuantTop <- KODMergeV1Quant[1:round(nrow(KODMergeV1Quant)/4),]
KODMergeV1Total <- arrange(KODMergeV1,desc(baseMean.TotalSeq))
KODMergeV1TotalTop <- KODMergeV1Total[1:round(nrow(KODMergeV1Total)/4),]


d14DMergeV1QuantImage <- ggplot(d14DMergeV1QuantTop, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM d14 ZScore (DEseq) with top Quant genes", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
d14DMergeV1TotalImage <- ggplot(KODMergeV1TotalTop, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM d14 ZScore (DEseq) with top Total genes", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 


KODMergeV1QuantImage <- ggplot(KODMergeV1QuantTop, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM KO ZScore (DEseq) with top Quant genes", x="TotalRNASeq", y="QuantSeq")  + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
KODMergeV1TotalImage <- ggplot(KODMergeV1TotalTop, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + theme(text = element_text(size=20)) +
  theme_bw() + labs(title="FUS ctrl vs HOM KO ZScore (DEseq) with top Total genes", x="TotalRNASeq", y="QuantSeq")  + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 

d14DMergeV1QuantImage+KODMergeV1QuantImage+d14DMergeV1TotalImage+KODMergeV1TotalImage

