library(DESeq2)
library(tidyverse)
library(pheatmap)

MergedReads <- read_csv("MergedReads.csv")
BMOutput <- read_csv("FTDBMOutput.csv")
F210ISupport <- read_csv("FTDSupportAS.csv")
MergedReads <- MergedReads %>% dplyr::select(EnsemblID, F210ISupport$Sample)
colnames(MergedReads)[2:ncol(MergedReads)] <- F210ISupport$UpdatedName


AllSupport <- dplyr::select(F210ISupport, UpdatedName, condition_A) %>% drop_na() %>% dplyr::rename(Condition = condition_A)
CTLC9Support <- dplyr::select(F210ISupport, UpdatedName, condition_CTL_C9) %>% drop_na() %>% dplyr::rename(Condition = condition_CTL_C9)
CTLTAUSupport <- dplyr::select(F210ISupport, UpdatedName, condition_CTL_TAU) %>% drop_na() %>% dplyr::rename(Condition = condition_CTL_TAU)
C9TAUSupport <- dplyr::select(F210ISupport, UpdatedName, condition_C9_TAU) %>% drop_na() %>% dplyr::rename(Condition = condition_C9_TAU)

NormaliseCounts <- function (Reads, SupportFile){
  #select only the reads in the relevant support file
  SelectedReads <- dplyr::select(Reads, EnsemblID, SupportFile$UpdatedName) %>% filter(grepl("EN", EnsemblID)) %>% 
    column_to_rownames("EnsemblID")

    #Get Normalised Reads from deseq
  formula0 = as.formula("~ 1") 
  CDS <- DESeqDataSetFromMatrix(countData = SelectedReads, colData = SupportFile, design = as.formula("~ Condition"))
  CDS <- estimateSizeFactors(CDS)
  
  #Extract normalised reads
  normalized_counts <- as.data.frame(counts(CDS, normalized=TRUE)) %>% rownames_to_column("EnsemblID")
  
  return(normalized_counts)
}

AllReads <- NormaliseCounts(MergedReads, AllSupport)
CTLC9Reads <- NormaliseCounts(MergedReads, CTLC9Support)
CTLTAUReads <- NormaliseCounts(MergedReads, CTLTAUSupport)
C9TAUReads <- NormaliseCounts(MergedReads, C9TAUSupport)

MakePCA <- function(Reads, Support){
  SelectedReads <- dplyr::select(Reads, EnsemblID, Support$UpdatedName)
  SelectedReads$RowSum <- SelectedReads %>% rowwise() %>% dplyr::select(-EnsemblID) %>% rowSums()
  PCAReads <- SelectedReads %>% filter(!RowSum == 0) %>% dplyr::select(-RowSum) %>% column_to_rownames("EnsemblID") 
  PCA <- prcomp(t(PCAReads), scale=TRUE)
  PCAData <- as.data.frame(PCA$x) %>% 
    rownames_to_column("UpdatedName") %>%
    left_join(Support, by="UpdatedName") 
  
  ggplot(PCAData, aes(x=PC1,y=PC2), label=UpdatedName) + geom_point(aes(color=Condition), size=4) +
    theme_bw(base_size=32) + 
    theme(legend.position="top") + geom_text(aes(label=UpdatedName),hjust=0, vjust=0)
  
}

MakePCA(AllReads, AllSupport) + ggtitle("Normalised Reads  PCA in all samples")
MakePCA(CTLC9Reads, CTLC9Support) + ggtitle("Normalised Reads PCA Expression in Control and C9")
MakePCA(CTLTAUReads, CTLTAUSupport) + ggtitle("Normalised Reads PCA Expression in Control and TAU")
MakePCA(C9TAUReads, C9TAUSupport) + ggtitle("Normalised Reads PCA Expression in C9 and TAU")

DESeqHeatmapNonSig <- function (Reads, SupportFile, QuanFilt){

  RowVariable <- Reads %>% mutate(RowVars = rowVars(as.matrix(Reads[,2:ncol(Reads)]))) %>% 
    dplyr::select("EnsemblID", "RowVars") %>% filter(quantile(RowVars, QuanFilt)<RowVars)
  
  SigReads <- filter(Reads,  EnsemblID %in% RowVariable$EnsemblID ) %>% na.omit()%>%
    `row.names<-`(., NULL) %>% column_to_rownames("EnsemblID")
  
  #Make a log2 reads heatmap
  pheatmap(log2(SigReads), show_rownames=F, scale="row", cluster_rows=FALSE, fontsize=15)
}

svg("AllSamplesNormalisedHeatmap.svg")
DESeqHeatmapNonSig(AllReads, AllSupport, 0.99)
dev.off()

svg("CTLvC9NormalisedHearmap.svg")
DESeqHeatmapNonSig(CTLC9Reads, CTLC9Support, 0.99)
dev.off()

svg("CTLvTAUNormalisedHearmap.svg")
DESeqHeatmapNonSig(CTLTAUReads, CTLTAUSupport, 0.99)
dev.off()

svg("C9vTAUNormalisedHearmap.svg")
DESeqHeatmapNonSig(C9TAUReads, C9TAUSupport, 0.99)
dev.off()