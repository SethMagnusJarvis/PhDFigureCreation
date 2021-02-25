library(pheatmap)
library(tidyverse)
library(DESeq2)

#Load DF of reads and support file, removing corrupted samples
RawReads <- read_csv("AdrianC9GeneReadsAll.csv")
Support <- read_csv("EasyNameSupport.csv") %>% 
  filter(Sample != "NM8082_800555" & Sample != "NM8084_800973" & Sample != "NM8093_793422") 

#Remove all genes that aren't a mouse ensembl ID and rename columns
RenamedRawReads <- dplyr::select(RawReads, EnsemblID, Support$Sample) %>% filter(grepl("ENSMUSG", EnsemblID)) 
colnames(RenamedRawReads)[2:ncol(RenamedRawReads)] <- Support$UpdatedName

SupportA <- dplyr::select(Support, Sample, UpdatedName, condition_A, type_3) %>% drop_na()

#Function to normalise reads
NormaliseCounts <- function (Reads, SupportFile){
  #select only the reads in the relevant support file
  SelectedReads <- dplyr::select(Reads, EnsemblID, SupportFile$UpdatedName) %>% filter(grepl("ENSMUSG", EnsemblID)) %>% 
    column_to_rownames("EnsemblID")
  
  #Change sample names so they're more human readable and it's immediately obvious which is WT and which is mutant
  coldata <- data.frame(sample = SupportFile$UpdatedName)
  coldata$Condition <- ifelse(SupportFile[,3] == "WT", "WT","Mut")
  row.names(coldata) <- coldata$sample
  coldata$sample <- NULL
  coldata$Condition <- factor(coldata$Condition,levels=c("WT","Mut"))
  
  #Get Normalised Reads from deseq
  formula0 = as.formula("~ 1") 
  CDS <- DESeqDataSetFromMatrix(countData = SelectedReads, colData = coldata, design = as.formula("~ Condition"))
  CDS <- estimateSizeFactors(CDS)
  
  #Extract normalised reads
  normalized_counts <- as.data.frame(counts(CDS, normalized=TRUE)) %>% rownames_to_column("EnsemblID")
  
  return(normalized_counts)
}

#Function to generate PCA plot
MakePCA <- function(Reads, Support){
  #Select only reads with mouse ensembl IDs
  SelectedReads <- dplyr::select(Reads, EnsemblID, Support$UpdatedName) %>% filter(grepl("ENSMUSG", EnsemblID))
  #Sum the rows and remove ones with no expression
  SelectedReads$RowSum <- SelectedReads %>% rowwise() %>% dplyr::select(-EnsemblID) %>% rowSums()
  PCAReads <- SelectedReads %>% filter(!RowSum == 0) %>% dplyr::select(-RowSum) %>% column_to_rownames("EnsemblID") 
  #Calculate PCA data and edit dataframe for plot
  PCA <- prcomp(t(PCAReads), scale=TRUE)
  PCAData <- as.data.frame(PCA$x) %>% 
    rownames_to_column("UpdatedName") %>%
    left_join(Support, by="UpdatedName") %>% 
    unite(SampleType, c("condition_A", "type_3")) 
  
  #Plot PCA
  ggplot(PCAData, aes(x=PC1,y=PC2), label=UpdatedName) + geom_point(aes(color=SampleType), size=4) +
    theme_bw(base_size=32) + 
    theme(legend.position="top") + geom_text(aes(label=UpdatedName),hjust=0, vjust=0)
  
}

#Normalise reads and make combined PCA
NormalisedReads <- NormaliseCounts(RenamedRawReads, Support)
MakePCA(NormalisedReads, SupportA)

#Make 6 and 12 month support files and PCAs
SixMonthSupport <- filter(Support, type_3 == "6M")
TwelveMonthSupport <- filter(Support, type_3 == "12M")

MakePCA(NormalisedReads, SixMonthSupport)
MakePCA(NormalisedReads, TwelveMonthSupport)
