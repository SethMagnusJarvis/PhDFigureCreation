library(tidyverse)
library(DESeq2)

############################################################################################################################
#Run DE on 4 possible options

MergedReads <- read_csv("MergedReads.csv")
BMOutput <- read_csv("FTDBMOutput.csv")
F210ISupport <- read_csv("FTDSupportAS.csv")
BrainType <- read_csv("FTDBrainBIAB.csv") %>% column_to_rownames("CellType") 
BrainType <- as.data.frame(t(BrainType)) %>% rownames_to_column("Sample")

#We converted age to a catagorical variable as including it as a continuous variable caused errors in DESeq2,
#either elading to no genes found to be differentially expressed or not working at all
F210ISupport$Age <- ifelse(F210ISupport$Age < 72, "Under72", "Over71")

AllSupport <- dplyr::select(F210ISupport, Sample, condition_A, Age, Sex) %>% drop_na()
CTLC9Support <- dplyr::select(F210ISupport, Sample, condition_CTL_C9, Age, Sex) %>% drop_na()
CTLTAUSupport <- dplyr::select(F210ISupport, Sample, condition_CTL_TAU, Age, Sex) %>% drop_na()
C9TAUSupport <- dplyr::select(F210ISupport, Sample, condition_C9_TAU, Age, Sex) %>% drop_na()


#Function to run differential expression
GetFullDEList <- function(AllReads, BioM, Support, BrainData, CTL, Sex, Age){
  
  #If the samples include control I set all samples in the condition to either control or mut so the fold change is correct
  #If they're not I just convert it to a factor
  if(CTL == TRUE){
    print("Converting to Mut")
    SupportFile <- Support
    #Format support file to just have sample names and Control/Mutant in the other column
    SupportFile[,2] <- ifelse(SupportFile[,2] == "CTL", "CTL","Mut") 
    SupportFile <- SupportFile %>% left_join(BrainData, by = "Sample") %>% column_to_rownames("Sample")
    #
    SupportFile$Condition <- factor(SupportFile[,1], levels=c("CTL","Mut"))
  }
  else{
    print("Leaving as is")
    SupportFile <- Support
    #Format support file to just have sample names and Control/Mutant in the other column
    SupportFile <- SupportFile %>% left_join(BrainData, by = "Sample") %>% column_to_rownames("Sample")
    SupportFile$Condition <- factor(SupportFile[,1])
  }
  
  #If sex is true then sex is included as a catagorical variable, the same is true with age
  ResFormula <- "~ "
  if(Sex == TRUE){ResFormula <- paste(ResFormula, "Sex + ")}
  if(Age == TRUE){ResFormula <- paste(ResFormula, "Age + ")}
  
  #I select the correct rows from reads
  Reads <- dplyr::select(AllReads, EnsemblID, rownames(SupportFile)) %>% column_to_rownames("EnsemblID")
  formula0 = as.formula("~ 1") 
  
  #Run DESeq including the interaction of the proportoins of all brain cell types with the condition, as well 
  #as the condition its self, and age and sex if included
  ResFormulaBrain <- as.formula(paste(ResFormula, "Condition:Astrocyte:Endothelial:Microglia:Mural:Neuron_All:Neuron_Interneuron:Neuron_Projection:Oligodendrocyte:Oligodendrocyte_Immature:RBC + Condition"))
  CDSBrain <- DESeqDataSetFromMatrix(countData = Reads, colData = SupportFile, design = ResFormulaBrain)
  CDSBrain <- DESeq(CDSBrain, test = "LRT", reduced = formula0, minReplicatesForReplace = 4 ) 
  DESeqResBrain <- results(CDSBrain)
  DESeqDFBrain <- as.data.frame(DESeqResBrain) %>% rownames_to_column("EnsemblID")
  
  #Run DESeq using just the condition along with age and sex if required
  ResFormulaPlain <- as.formula(paste(ResFormula, "Condition"))
  CDSPlain <- DESeqDataSetFromMatrix(countData = Reads, colData = SupportFile, design = ResFormulaPlain)
  CDSPlain <- DESeq(CDSPlain, test = "LRT", reduced = formula0, minReplicatesForReplace = 4 ) 
  DESeqResPlain <- results(CDSPlain)
  DESeqDFPlain <- as.data.frame(DESeqResPlain) %>% rownames_to_column("EnsemblID")
  
  #Join the two differentially expressed datasets together, with a column for which genes are significant 
  JoinBoth <- full_join(DESeqDFBrain, DESeqDFPlain, by = "EnsemblID", suffix = c(".Join", ".Solo")) %>%
    left_join(BioM, by = c("EnsemblID" = "ensembl_gene_id")) %>% drop_na("external_gene_name") %>%
    dplyr::select(external_gene_name, everything())
  return(JoinBoth)
}

#Get DE genes in all samples and write to a CSV
AllJoinedDE <- GetFullDEList(MergedReads, BMOutput, AllSupport, BrainType, TRUE, TRUE, TRUE)
CTLC9JoinedDE <- GetFullDEList(MergedReads, BMOutput, CTLC9Support, BrainType, TRUE, TRUE, TRUE)
CTLTAUJoinedDE <- GetFullDEList(MergedReads, BMOutput, CTLTAUSupport, BrainType, TRUE, TRUE, TRUE)
C9TAUJoinedDE <- GetFullDEList(MergedReads, BMOutput, C9TAUSupport, BrainType, FALSE, TRUE, TRUE)

write_csv(AllJoinedDE, "CtlVsMutDEBothCatAge.csv")
write_csv(CTLC9JoinedDE, "CtlVsC9DEBothCatAge.csv")
write_csv(CTLTAUJoinedDE, "CtlVsTAUDEBothCatAge.csv")
write_csv(C9TAUJoinedDE, "C9VsTAUDEBothCatAge.csv")