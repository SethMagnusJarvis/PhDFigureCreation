library(tidyverse)
library(DESeq2)
library(BrainInABlender)
library(pheatmap)

#Import files
MergedReads <- read_csv("AdrianC9GeneReadsAll.csv")

BMOutput <- read_csv("AdrianC9MouseGenes.csv")

Support <- read_csv("AdrianC9SupportUpdated.csv") %>% filter(Sample != "NM8082_800555" & Sample != "NM8084_800973" & Sample != "NM8093_793422") 

SupportFile <- Support %>% dplyr::select(Sample, condition_A)

#Filter the reads
FiltReads <- dplyr::select(MergedReads, EnsemblID, SupportFile$Sample) %>% filter(grepl("ENSMUSG", EnsemblID)) 
FiltReads <- FiltReads %>% column_to_rownames("EnsemblID") 

#Set condition to either wildtype or Mut and format for DESeq usage
SupportFile$condition_A<- ifelse(SupportFile[,2] == "WT", "CTL","Mut")
row.names(SupportFile) <- SupportFile$Sample
SupportFile$Sample <- NULL
SupportFile$Condition <- factor(SupportFile$condition_A, levels=c("CTL","Mut"))

#Get Normalised Reads from deseq
formula0 = as.formula("~ 1") 
CDS <- DESeqDataSetFromMatrix(countData = FiltReads, colData = SupportFile, design = as.formula("~ Condition"))
CDS <- estimateSizeFactors(CDS)

#Extract normalised reads
normalised_counts <- as.data.frame(counts(CDS, normalized=TRUE)) %>% rownames_to_column("EnsemblID") %>%
  left_join(BMOutput, by=c("EnsemblID" = "ensembl_gene_id")) %>% 
  dplyr::select(-one_of("gene_biotype", "EnsemblID")) %>% drop_na("external_gene_name") %>% dplyr::select("external_gene_name", everything())

BrainType <- Sir_UnMixALot(userInput=normalised_counts, dataColumns=c(2:ncol(normalised_counts)), geneColumn=1, species="mouse")$AveragePrimary_CellTypeIndex
BrainType <- as.data.frame(BrainType) %>% rownames_to_column("CellType")

write_csv(as.data.frame(BrainType), "C9BrainBIAB.csv")


###########################################################################################################################
#Cell type heatmaps
#Read in cells
BrainType <- read_csv("C9BrainBIAB.csv")
Support <- read_csv("EasyNameSupport.csv") %>% filter(Sample != "NM8082_800555" & Sample != "NM8084_800973" & Sample != "NM8093_793422") 
SupportBrain <- dplyr::select(Support, Sample, UpdatedName, condition_A, type_3) %>% drop_na()

RenamedBrain <- dplyr::select(BrainType, CellType, Support$Sample)
colnames(RenamedBrain)[2:ncol(RenamedBrain)] <- Support$UpdatedName

#Make heatmap of cell type proportoins in all samples
PheatBrain <- column_to_rownames(RenamedBrain, "CellType")
svg("BIABHeatmap.svg", width=10, height = 8)
pheatmap::pheatmap(PheatBrain, cluster_rows=FALSE,  scale="row", fontsize = 12 )
dev.off()

#make PCA plot of cell type in all samples
PCAType <- prcomp(t(PheatBrain), scale=TRUE)
#
PCARes <- as.data.frame(PCAType$x) %>% 
  rownames_to_column("UpdatedName") %>%
  left_join(SupportBrain, by="UpdatedName")

PCACellType <- ggplot(PCARes, aes(x=PC1,y=PC2)) + geom_point(aes(color=condition_A), size=4) +
  theme_bw(base_size=32) + 
  theme(legend.position="top") + geom_text(aes(size = 8, label=UpdatedName),hjust=0, vjust=0, nudge_x = -1) + ggtitle("Cell Type PCA All Samples") + 
  theme(
    legend.title = element_text(color = "blue", size = 14),
    legend.text = element_text(color = "red", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

svg("PCACellTypeAll.svg", width = 10, height = 7)
PCACellType
dev.off()


#Make 6 and 12 Month old PCAs
SixMonthSupport <- filter(Support, type_3 == "6M")
TwelveMonthSupport <- filter(Support, type_3 == "12M")

PCAReads6M <- dplyr::select(PheatBrain, SixMonthSupport$UpdatedName) 
PCAType6M <- prcomp(t(PCAReads6M), scale=TRUE)
PCARes6Month <- as.data.frame(PCAType6M$x) %>% 
  rownames_to_column("UpdatedName") %>%
  left_join(SixMonthSupport, by="UpdatedName")

PCACellType6M <- ggplot(PCARes6Month, aes(x=PC1,y=PC2)) + geom_point(aes(color=condition_A), size=4) +
  theme_bw(base_size=32) + 
  theme(legend.position="top") + geom_text(aes(size = 10, label=UpdatedName),hjust=0, vjust=0, nudge_x = -0.6) + ggtitle ("Cell Type PCA 6 Month Samples") +
  theme(
    legend.title = element_text(color = "blue", size = 14),
    legend.text = element_text(color = "red", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

svg("PCACellType6M.svg", width = 10, height = 7)
PCACellType6M
dev.off()


PCAReads12M <- dplyr::select(PheatBrain, TwelveMonthSupport$UpdatedName) 
PCAType12M <- prcomp(t(PCAReads12M), scale=TRUE)
PCARes12Month <- as.data.frame(PCAType12M$x) %>% 
  rownames_to_column("UpdatedName") %>%
  left_join(TwelveMonthSupport, by="UpdatedName")

PCACellType12M <- ggplot(PCARes12Month, aes(x=PC1,y=PC2)) + geom_point(aes(color=condition_A), size=4) +
  theme_bw(base_size=32) + 
  theme(legend.position="top") + geom_text(aes(size = 10, label=UpdatedName),hjust=0, vjust=0, nudge_x = -0.75) +  ggtitle("Cell Type PCA 12 Month Samples") +
  theme(
    legend.title = element_text(color = "blue", size = 14),
    legend.text = element_text(color = "red", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

svg("PCACellType12M.svg", width = 10, height = 7)
PCACellType12M
dev.off()

