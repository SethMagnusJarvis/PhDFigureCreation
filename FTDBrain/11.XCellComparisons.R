library(tidyverse)
library(DESeq2)
library(xCell)
library(ggpubr)
library(limma)

############################################################################################################################
############################################################################################################################
                                                          #Run XCell#
############################################################################################################################
############################################################################################################################


MergedReads <- read_csv("MergedReads.csv")

BMOutput <- read_csv("FTDBMOutput.csv")

F210ISupport <- read_csv("FTDSupport.csv")[,1:2]

Reads <- MergedReads
SupportFile <- F210ISupport
BioMOutput <- BMOutput

NormaliseReads <- function (Reads, SupportFile, BioMOutput){
  #Change sample names so they're more human readable and it's immediately obvious which is WT and which is mutant
  Reads <- column_to_rownames(Reads, "EnsemblID")
  SupportFile$condition_A<- ifelse(SupportFile[,2] == "CTL", "CTL","Mut")
  row.names(SupportFile) <- SupportFile$Sample
  SupportFile$Sample <- NULL
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


NormalisedUnfilteredReads <- NormaliseReads(MergedReads, F210ISupport) 
NormalisedUnfilteredReads <- NormalisedUnfilteredReads[!duplicated(NormalisedUnfilteredReads[, "external_gene_name"]), ]
rownames(NormalisedUnfilteredReads) <- NULL
NormalisedUnfilteredReads <- NormalisedUnfilteredReads %>% column_to_rownames("external_gene_name")

XCellType <- as.data.frame(xCellAnalysis(NormalisedUnfilteredReads)) %>% rownames_to_column("CellType")


write_csv(XCellType, "XCellTypesAll.csv")

############################################################################################################################
############################################################################################################################
                                                  #Boxplots#
############################################################################################################################
############################################################################################################################

XCellType <- read_csv("XCellTypesAll.csv") %>% column_to_rownames("CellType") 
#XCellTypeBrain <- XCellType[c("Astrocytes", "Neurons"),]
BrainType <- as.data.frame(t(XCellType)) %>% rownames_to_column("Sample")


BrainMerged <- BrainType %>% right_join(F210ISupport, by="Sample") %>% dplyr::rename(SmoothMuscle ="Smooth muscle") %>% column_to_rownames("Sample")

MyComparisons <- list(c("CTL", "FTD_C9"), c("CTL", "FTD_TAU"), c("FTD_C9", "FTD_TAU"))
AstrocyteBox <- ggboxplot(BrainMerged, x="condition_A", y = "Astrocytes",
                          color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                          add = "jitter", shape = "condition_A")
NeuronsBox <- ggboxplot(BrainMerged, x="condition_A", y = "Neurons",
                        color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                        add = "jitter", shape = "condition_A")
aDCBox <- ggboxplot(BrainMerged, x="condition_A", y = "aDC",
                    color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                    add = "jitter", shape = "condition_A")
FibroblastBox <- ggboxplot(BrainMerged, x="condition_A", y = "Fibroblasts",
                           color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                           add = "jitter", shape = "condition_A")
HSCBox <- ggboxplot(BrainMerged, x="condition_A", y = "HSC",
                    color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                    add = "jitter", shape = "condition_A")
iDCBox <- ggboxplot(BrainMerged, x="condition_A", y = "iDC",
                    color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                    add = "jitter", shape = "condition_A")
MSCBox <- ggboxplot(BrainMerged, x="condition_A", y = "MSC",
                    color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                    add = "jitter", shape = "condition_A")
NKTBox <- ggboxplot(BrainMerged, x="condition_A", y = "NKT",
                    color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                    add = "jitter", shape = "condition_A")
SmoothMuscleBox <- ggboxplot(BrainMerged, x="condition_A", y = "SmoothMuscle",
                             color = "condition_A", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                             add = "jitter", shape = "condition_A")


pdf("XCellBoxplotBrain.pdf")
AstrocyteBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
NeuronsBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
dev.off()

pdf("XCellBoxplotHighMean.pdf")
aDCBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
FibroblastBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
MSCBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
NeuronsBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
NKTBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
SmoothMuscleBox + stat_compare_means(comparisons = MyComparisons) +stat_compare_means(label.y = 0.05)
dev.off()

############################################################################################################################
############################################################################################################################
                                          #Differential Expression#
############################################################################################################################
############################################################################################################################
MergedReads <- read_csv("MergedReads.csv")
BMOutput <- read_csv("FTDBMOutput.csv")
F210ISupport <- read_csv("FTDSupportAS.csv")
XCellType <- read_csv("XCellTypesAll.csv") 

F210ISupport$Age <- ifelse(F210ISupport$Age < 72, "Under72", "Over71")

#Get rowmeans and 
XCellTypeWMean <- XCellType
XCellTypeWMean$RowMean <- XCellTypeWMean %>% rowwise() %>% dplyr::select(-CellType) %>% rowMeans()
XCellTypeFilt <- XCellTypeWMean %>% filter(RowMean > 0.01) %>% dplyr::select(-RowMean) %>% column_to_rownames("CellType")
FiltType <- as.data.frame(t(XCellTypeFilt)) %>% rownames_to_column("Sample")

XCellTypeBrain <- XCellType %>% filter(CellType == "Astrocytes" | CellType == "Neurons") %>% column_to_rownames("CellType")
BrainType <- as.data.frame(t(XCellTypeBrain)) %>% rownames_to_column("Sample")

FiltMerged <- FiltType %>% right_join(F210ISupport, by="Sample") %>% dplyr::rename(SmoothMuscle ="Smooth muscle")
BrainMerged <- BrainType %>% right_join(F210ISupport, by="Sample")

AllSupportMean <- dplyr::select(FiltMerged, Sample, condition_A, aDC:SmoothMuscle, Age, Sex) %>% drop_na()
CTLC9SupportMean <- dplyr::select(FiltMerged, Sample, condition_CTL_C9, aDC:SmoothMuscle, Age, Sex) %>% drop_na()
CTLTAUSupportMean <- dplyr::select(FiltMerged, Sample, condition_CTL_TAU, aDC:SmoothMuscle, Age, Sex) %>% drop_na()
C9TAUSupportMean <- dplyr::select(FiltMerged, Sample, condition_C9_TAU, aDC:SmoothMuscle, Age, Sex) %>% drop_na()

GetFullDEList <- function(AllReads, BioM, Support, BrainData, CTL, Sex, Age){
  
  if(CTL == TRUE){
    print("Converting to Mut")
    SupportFile <- Support
    #Format support file to just have sample names and Control/Mutant in the other column
    SupportFile[,2] <- ifelse(SupportFile[,2] == "CTL", "CTL","Mut") 
    SupportFile <- SupportFile %>% column_to_rownames("Sample")
    SupportFile$Condition <- factor(SupportFile[,1], levels=c("CTL","Mut"))
  }
  else{
    print("Leaving as is")
    SupportFile <- Support
    #Format support file to just have sample names and Control/Mutant in the other column
    SupportFile <- SupportFile %>%  column_to_rownames("Sample")
    SupportFile$Condition <- factor(SupportFile[,1])
  }
  
  ResFormula <- "~ "
  if(Sex == TRUE){ResFormula <- paste(ResFormula, "Sex + ")}
  if(Age == TRUE){ResFormula <- paste(ResFormula, "Age + ")}
  
  
  Reads <- dplyr::select(AllReads, EnsemblID, rownames(SupportFile)) %>% column_to_rownames("EnsemblID")
  formula0 = as.formula("~ 1") 
  
  
  ResFormulaBrain <- as.formula(paste(ResFormula, "Condition:Neurons + Condition"))
  
  CDSBrain <- DESeqDataSetFromMatrix(countData = Reads, colData = SupportFile, design = ResFormulaBrain)
  CDSBrain <- DESeq(CDSBrain, test = "LRT", reduced = formula0, minReplicatesForReplace = 4 ) 
  DESeqResBrain <- results(CDSBrain)
  DESeqDFBrain <- as.data.frame(DESeqResBrain) %>% rownames_to_column("EnsemblID")
  
  #Run DESeq using just the condition as a covariate
  ResFormulaPlain <- as.formula(paste(ResFormula, "Condition"))
  CDSPlain <- DESeqDataSetFromMatrix(countData = Reads, colData = SupportFile, design = ResFormulaPlain)
  CDSPlain <- DESeq(CDSPlain, test = "LRT", reduced = formula0, minReplicatesForReplace = 4 ) 
  DESeqResPlain <- results(CDSPlain)
  DESeqDFPlain <- as.data.frame(DESeqResPlain) %>% rownames_to_column("EnsemblID")
  
  #Join the two differentially expressed datasets together, with a column for which genes are significant 
  #and one for if the pvalue is the same or different in the two datasets (probably not useful)
  #I couldn't get this working and I don't know why
  JoinBoth <- full_join(DESeqDFBrain, DESeqDFPlain, by = "EnsemblID", suffix = c(".Join", ".Solo")) %>%
    left_join(BioM, by = c("EnsemblID" = "ensembl_gene_id")) %>% drop_na("external_gene_name") %>%
    dplyr::select(external_gene_name, everything())
  return(JoinBoth)
}

AllJoinedDEBrain <- GetFullDEList(MergedReads, BMOutput, AllSupportMean, BrainType, TRUE, TRUE, TRUE)
CTLC9JoinedDEBrain <- GetFullDEList(MergedReads, BMOutput, CTLC9SupportMean, BrainType, TRUE, TRUE, TRUE)
CTLTAUJoinedDEBrain <- GetFullDEList(MergedReads, BMOutput, CTLTAUSupportMean, BrainType, TRUE, TRUE, TRUE)
C9TAUJoinedDEBrain <- GetFullDEList(MergedReads, BMOutput, C9TAUSupportMean, BrainType, FALSE, TRUE, TRUE)

write_csv(AllJoinedDEBrain, "CTLVsMUTDEXCellNeuron.csv")
write_csv(CTLC9JoinedDEBrain, "CTLVsC9DEXCellNeuron.csv")
write_csv(CTLTAUJoinedDEBrain, "CTLVsTAUDEXCellNeuron.csv")
write_csv(C9TAUJoinedDEBrain, "C9VsTAUDEXCellNeuron.csv")

############################################################################################################################
############################################################################################################################
                                          #Venn Diagrams Cell Types#
############################################################################################################################
############################################################################################################################
AllJoinedDEBrain <- read_csv("CTLVsMUTDEXCellNeuron.csv")
CTLC9JoinedDEBrain <- read_csv("CTLVsC9DEXCellNeuron.csv")
CTLTAUJoinedDEBrain <- read_csv("CTLVsTAUDEXCellNeuron.csv")
C9TAUJoinedDEBrain <- read_csv("C9VsTAUDEXCellNeuron.csv")

ListCompare <- function(DEJoin){
  DE1SigList <- DEJoin %>% dplyr::select(external_gene_name, contains(".Join")) %>% filter(padj.Join < 0.05) %>% 
    dplyr::select(external_gene_name)
  DE2SigList <- DEJoin %>% dplyr::select(external_gene_name, contains(".Solo")) %>% filter(padj.Solo < 0.05) %>% 
    dplyr::select(external_gene_name)
  ids = sort(unique(c(as.character(unlist(DE1SigList)),  as.character(unlist(DE2SigList) ))))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% unlist(DE1SigList)
    counts[i, 2] = ids[i] %in% unlist(DE2SigList)
  }
  colnames(counts) = c("Join", "Solo")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

AllCompared <- ListCompare(AllJoinedDEBrain)
CTLC9Compared <- ListCompare(CTLC9JoinedDEBrain)
CTLTAUCompared <- ListCompare(CTLTAUJoinedDEBrain)
C9TAUCompared <- ListCompare(C9TAUJoinedDEBrain)

par(mfrow=c(2,2))
vennDiagram(AllCompared) + title("Control Vs Mutant Venn")
vennDiagram(CTLC9Compared) + title("Control vs C9 Venn")
vennDiagram(CTLTAUCompared) + title("Control vs TAU Venn")
vennDiagram(C9TAUCompared) + title("C9 vs TAU Venn")
