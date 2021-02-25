library(pheatmap)
library(tidyverse)
library(DESeq2)

RawReads <- read_csv("AdrianC9GeneReadsAll.csv")
Support <- read_csv("EasyNameSupport.csv") %>% filter(Sample != "NM8082_800555")

RenamedRawReads <- dplyr::select(RawReads, EnsemblID, Support$Sample) %>% filter(grepl("ENSMUSG", EnsemblID)) 
colnames(RenamedRawReads)[2:ncol(RenamedRawReads)] <- Support$UpdatedName

PR12MDE <- read_csv("DifferentialExpressionPR12MNoType.csv") %>% filter(ADJ.PVAL < 0.05)
GR12MDE <- read_csv("DifferentialExpressionGR12MNoType.csv") %>% filter(ADJ.PVAL < 0.05)

SupportPR6M <- dplyr::select(Support, Sample, UpdatedName, condition_B, type_3) %>% drop_na()
SupportGR6M <- dplyr::select(Support, Sample, UpdatedName, condition_C, type_3) %>% drop_na()
SupportPR12M <- dplyr::select(Support, Sample, UpdatedName, condition_D, type_3) %>% drop_na()
SupportGR12M <- dplyr::select(Support, Sample, UpdatedName, condition_E, type_3) %>% drop_na()

DESeqHeatmapNonSig <- function (Reads, SupportFile, QuanFilt){
  #select only the reads in the relevant support file
  SelectedReads <- dplyr::select(Reads, EnsemblID, SupportFile$UpdatedName) %>% filter(grepl("ENSMUSG", EnsemblID))  
  
  #Filter to the Most Variable significant genes
  RowVariable <- SelectedReads %>% mutate(RowVars = rowVars(as.matrix(SelectedReads[,2:ncol(SelectedReads)]))) %>% 
    dplyr::select("EnsemblID", "RowVars") %>% filter(quantile(RowVars, QuanFilt)<RowVars)
  
  SelectedReads <- SelectedReads %>% column_to_rownames("EnsemblID")
  
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
  
  #Filter to only use the significant reads
  SigReads <- filter(normalized_counts,  EnsemblID %in% RowVariable$EnsemblID ) %>% na.omit()%>%
    `row.names<-`(., NULL) %>% column_to_rownames("EnsemblID")
  
  #Make a log2 reads heatmap
  pheatmap(log2(SigReads), show_rownames=F, scale="row", cluster_rows=FALSE, fontsize = 15)
}

svg("HeatmapPR6M.svg")
DESeqHeatmapNonSig(RenamedRawReads, SupportPR6M, 0.996)
dev.off()

svg("HeatmapGR6M.svg")
DESeqHeatmapNonSig(RenamedRawReads, SupportGR6M, 0.996)
dev.off()

svg("HeatmapPR12M.svg")
DESeqHeatmapNonSig(RenamedRawReads, SupportPR12M, 0.996)
dev.off()

svg("HeatmapGR12M.svg")
DESeqHeatmapNonSig(RenamedRawReads, SupportGR12M, 0.996)
dev.off()

