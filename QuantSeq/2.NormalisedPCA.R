library(pheatmap)
library(tidyverse)
library(DESeq2)

d14Reads <- read_csv("RawD14ReadsBoth.csv") %>% dplyr::select(gene_id, !contains("HET"), -X1) %>% 
  dplyr::select(gene_id, contains("WT"), everything())
KOReads <- read_csv("RawKOReadsBoth.csv") %>% dplyr::select(gene_id, !contains("HET"), -X1)  %>% 
  dplyr::select(gene_id, contains("WT"), everything())

d14Support <- NULL
d14Support$sample <- colnames(d14Reads)[-1]
d14Support$condition_HOM <- ifelse(grepl("WT", d14Support$sample), "WT", "HOM")
d14Support <- as.data.frame(d14Support)

KOSupport <- NULL
KOSupport$sample <- colnames(KOReads)[-1]
KOSupport$condition_HOM <- ifelse(grepl("WT", KOSupport$sample), "WT", "HOM")
KOSupport <- as.data.frame(KOSupport)

NormaliseCounts <- function (Reads, SupportFile){
  #select only the reads in the relevant support file
  SelectedReads <- Reads %>% filter(grepl("ENSMUSG", gene_id)) %>% 
    column_to_rownames("gene_id") %>% drop_na()
  
  #Change sample names so they're more human readable and it's immediately obvious which is WT and which is mutant
  coldata <- data.frame(Sample = SupportFile$sample)
  coldata$Condition <- ifelse(SupportFile$condition_HOM == "WT", "WT","Mut")
  row.names(coldata) <- coldata$Sample
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

d14Normalised <- NormaliseCounts(d14Reads, d14Support)
KONormalised <- NormaliseCounts(KOReads, KOSupport)


MakePCA <- function(Reads, Support){
  SelectedReads <- dplyr::select(Reads, EnsemblID, Support$sample) %>% filter(grepl("ENSMUSG", EnsemblID))
  SelectedReads$RowSum <- SelectedReads %>% rowwise() %>% dplyr::select(-EnsemblID) %>% rowSums()
  PCAReads <- SelectedReads %>% filter(!RowSum == 0) %>% dplyr::select(-RowSum) %>% column_to_rownames("EnsemblID") 
  PCA <- prcomp(t(PCAReads), scale=TRUE)
  PCAData <- as.data.frame(PCA$x) %>% 
    rownames_to_column("sample") %>%
    left_join(Support, by="sample") 
  
  ggplot(PCAData, aes(x=PC1,y=PC2), label=sample) + geom_point(aes(color=condition_HOM), size=4) +
    theme_bw(base_size=32) + 
    theme(legend.position="top") + geom_text(aes(label=sample),hjust=0, vjust=0)
  
}

MakePCA(d14Normalised, d14Support)
MakePCA(KONormalised, KOSupport)
