library(tidyverse)
library(DESeq2)

MergedReads <- read_csv("MergedReads.csv")
BMOutput <- read_csv("FTDBMOutput.csv")
F210ISupport <- read_csv("FTDSupportAS.csv")
BrainType <- read_csv("FTDBrainBIAB.csv") %>% column_to_rownames("CellType") 
BrainType <- as.data.frame(t(BrainType)) %>% rownames_to_column("Sample")
F210ISupport$Age <- ifelse(F210ISupport$Age < 72, "Under72", "Over71")

SupportWBrain <- left_join(F210ISupport, BrainType, by="Sample")



AllSupport <- dplyr::select(SupportWBrain, Sample, condition_A, Astrocyte:RBC, Age, Sex) %>% drop_na()
CTLC9Support <- dplyr::select(SupportWBrain, Sample, condition_CTL_C9, Astrocyte:RBC, Age, Sex) %>% drop_na()
CTLTAUSupport <- dplyr::select(SupportWBrain, Sample, condition_CTL_TAU, Astrocyte:RBC, Age, Sex) %>% drop_na()
C9TAUSupport <- dplyr::select(SupportWBrain, Sample, condition_C9_TAU, Astrocyte:RBC, Age, Sex) %>% drop_na()


AllReads <- MergedReads
BioM <- BMOutput
Support <- AllSupport
BrainData <- BrainType
CTL <- TRUE

GetFullDEList <- function(AllReads, BioM, Support, BrainData, CTL){
  
  SupportFile <- Support
  
  if(CTL == TRUE){
    print("Converting to Mut")
    #Format support file to just have sample names and Control/Mutant in the other column
    SupportFile[,2] <- ifelse(SupportFile[,2] == "CTL", "CTL","Mut") 
    SupportFile <- SupportFile %>% column_to_rownames("Sample")
    SupportFile$Condition <- factor(SupportFile[,1], levels=c("CTL","Mut"))
  }
  else{
    print("Leaving as is")
    SupportFile <- SupportFile %>% column_to_rownames("Sample")
    SupportFile$Condition <- factor(SupportFile[,1])
  }
  
  CellTypes <- c("Astrocyte", "Endothelial", "Microglia", "Mural", "Neuron_All", "Neuron_Interneuron", "Neuron_Projection", "Oligodendrocyte", "Oligodendrocyte_Immature", "RBC")
  FormulaList <- list()
  
  Reads <- dplyr::select(AllReads, EnsemblID, rownames(SupportFile)) %>% column_to_rownames("EnsemblID")
  
  DE.output <- list()
  
  for (j in 1:length(CellTypes))  { 
    
    cat(paste0("Starting DESeq2 loop ", j, " of ", length(CellTypes)))
    
    
    formula0 = as.formula("~ 1") 
    ResFormula <- as.formula(paste0("~ Age + Sex + ", CellTypes[j], " + Condition"))
    CDS <- DESeqDataSetFromMatrix(countData = Reads, colData = SupportFile, design = ResFormula)
    CDS <- DESeq(CDS, test = "LRT", reduced = formula0, minReplicatesForReplace = 4 ) 
    DESeqRes <- results(CDS)
    DE.output[[j]] <- as.data.frame(DESeqRes) %>% rownames_to_column("EnsemblID")
    
  }
  return(DE.output)
}

AllJoinedDE <- GetFullDEList(MergedReads, BMOutput, AllSupport, BrainType, TRUE)
CTLC9JoinedDE <- GetFullDEList(MergedReads, BMOutput, CTLC9Support, BrainType, TRUE)
CTLTAUJoinedDE <- GetFullDEList(MergedReads, BMOutput, CTLTAUSupport, BrainType, TRUE)
C9TAUJoinedDE <- GetFullDEList(MergedReads, BMOutput, C9TAUSupport, BrainType, FALSE)


save(AllJoinedDE, file = "CTLMUTDECellTypeCatAge.rds")
save(CTLC9JoinedDE, file = "CTLC9DECellTypeCatAge.rds")
save(CTLTAUJoinedDE, file = "CTLTAUDECellTypeCatAge.rds")
save(C9TAUJoinedDE, file = "C9TAUDECellTypeCatAge.rds")


CountSig <- function(SigMatrix){
  CellTypes <- c("Astrocyte", "Endothelial", "Microglia", "Mural", "Neuron_All", "Neuron_Interneuron", "Neuron_Projection", "Oligodendrocyte", "Oligodendrocyte_Immature", "RBC")
  CollapsedSig <- NULL
  CollapsedSig$Names <- unlist(CellTypes)
  for(x in 1:10){
    CollapsedSig$Count[x] <- (nrow(filter(SigMatrix, SigMatrix[,x+1] < 0.05)))
  }
  CollapsedSig <- as.data.frame(CollapsedSig)
  return(CollapsedSig)
}

CellTypes <- c("Astrocyte", "Endothelial", "Microglia", "Mural", "Neuron_All", "Neuron_Interneuron", "Neuron_Projection", "Oligodendrocyte", "Oligodendrocyte_Immature", "RBC")

AllCollapsed <- AllJoinedDE %>% purrr::reduce(inner_join, by = "EnsemblID") %>% dplyr::select(EnsemblID, contains("padj"))
colnames(AllCollapsed)[2:11] <- CellTypes
AllCollapsedSig <- CountSig(AllCollapsed)

CTLC9Collapsed <- CTLC9JoinedDE %>% purrr::reduce(inner_join, by = "EnsemblID") %>% dplyr::select(EnsemblID, contains("padj"))    
colnames(CTLC9Collapsed)[2:11] <- CellTypes
CTLC9CollapsedSig <- CountSig(CTLC9Collapsed)

CTLTAUCollapsed <- CTLTAUJoinedDE %>% purrr::reduce(inner_join, by = "EnsemblID") %>% dplyr::select(EnsemblID, contains("padj"))    
colnames(CTLTAUCollapsed)[2:11] <- CellTypes
CTLTAUCollapsedSig <- CountSig(CTLTAUCollapsed)

C9TAUCollapsed <- C9TAUJoinedDE %>% purrr::reduce(inner_join, by = "EnsemblID") %>% dplyr::select(EnsemblID, contains("padj"))    
colnames(C9TAUCollapsed)[2:11] <- CellTypes
C9TAUCollapsedSig <- CountSig(C9TAUCollapsed)

write_csv(AllCollapsed, "CTLMUTDECellTypeCatAge.csv")
write_csv(AllCollapsedSig, "CTLMUTDECellTypeSigCountCatAge.csv")

write_csv(CTLTAUCollapsed, "CTLC9DECellTypeCatAge.csv")
write_csv(CTLTAUCollapsedSig, "CTLC9DECellTypeSigCountCatAge.csv")

write_csv(CTLTAUCollapsed, "CTLTAUDECellTypeCatAge.csv")
write_csv(CTLTAUCollapsedSig, "CTLTAUDECellTypeSigCountCatAge.csv")

write_csv(C9TAUCollapsed, "C9TAUDECellTypeCatAge.csv")
write_csv(C9TAUCollapsedSig, "C9TAUDECellTypeSigCountCatAge.csv")

MutC9Merged <- full_join(AllCollapsedSig, CTLC9CollapsedSig, by = "Names", suffix=c("WTMut", "WTC9"))
TAUBothMutMerged <- full_join(CTLTAUCollapsedSig, C9TAUCollapsedSig, by = "Names", suffix=c("WTTAU", "C9TAU"))
AllSigMerged <- full_join(MutC9Merged, TAUBothMutMerged)
write_csv(AllSigMerged, "MergedBIABCellTypeSig.csv")
