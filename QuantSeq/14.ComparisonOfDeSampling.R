library(tidyverse)
library(ggplot2)
library(patchwork)
library(limma)

#Load all files without a filter
KOQuantV2De <- read_csv("KOQuantV2De.csv")
d14QuantV2De <- read_csv("d14QuantV2De.csv")
KOOriginTotalDe <- read_csv("KOOriginTotalDe.csv")
d14OriginTotalDe <- read_csv("d14OriginTotalDe.csv")
KO3MilTotalDe <- read_csv("KO3MilTotalDe.csv")
d143MilTotalDe <- read_csv("d143MilTotalDe.csv")
KO10MilTotalDe <- read_csv("KO10MilTotalDe.csv")
d1410MilTotalDe <- read_csv("d1410MilTotalDe.csv")

#Load all files with a filter
KOQuantV1DeNoFilt <- read_csv("KOQuantV1DeNoFilt.csv")
d14QuantV1DeNoFilt <- read_csv("d14QuantV1DeNoFilt.csv")
KOQuantV2DeNoFilt <- read_csv("KOQuantV2DeNoFilt.csv")
d14QuantV2DeNoFilt <- read_csv("d14QuantV2DeNoFilt.csv")
KOOriginTotalDeNoFilt <- read_csv("KOOriginTotalDeNoFilt.csv")
d14OriginTotalDeNoFilt <- read_csv("d14OriginTotalDeNoFilt.csv")
KO3MilTotalDeNoFilt <- read_csv("KO3MilTotalDeNoFilt.csv")
d143MilTotalDeNoFilt <- read_csv("d143MilTotalDeNoFilt.csv")
KO10MilTotalDeNoFilt <- read_csv("KO10MilTotalDeNoFilt.csv")
d1410MilTotalDeNoFilt <- read_csv("d1410MilTotalDeNoFilt.csv")
KO15MilTotalDeNoFilt <- read_csv("KO15MilTotalDeNoFilt.csv")
d1415MilTotalDeNoFilt <- read_csv("d1415MilTotalDeNoFilt.csv")
KO20MilTotalDeNoFilt <- read_csv("KO20MilTotalDeNoFilt.csv")
d1420MilTotalDeNoFilt <- read_csv("d1420MilTotalDeNoFilt.csv")

#Bind Matrixes together
KOOriginMerge <- left_join(KOQuantV2DeNoFilt, KOOriginTotalDeNoFilt, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))
KO3MilMerge <- left_join(KOQuantV2DeNoFilt, KO3MilTotalDeNoFilt, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))
KO10MilMerge <- left_join(KOQuantV2DeNoFilt, KO10MilTotalDeNoFilt, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))
d14OriginMerge <- left_join(d14QuantV2DeNoFilt, d14OriginTotalDeNoFilt, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))
d143MilMerge <- left_join(d14QuantV2DeNoFilt, d143MilTotalDeNoFilt, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))
d1410MilMerge <- left_join(d14QuantV2DeNoFilt, d1410MilTotalDeNoFilt, by = "gene_id", suffix = c(".QuantSeq", ".TotalSeq"))

#Make Z-score plots
d14DESeqOriginImage <- ggplot(d14OriginMerge, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM d14 ZScore (1.5 Million Reads)", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
d14DESeq3MilImage <- ggplot(d143MilMerge, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM d14 ZScore (3 Million Reads)", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
d14DESeq10MilImage <- ggplot(d1410MilMerge, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM d14 ZScore (10 Million Reads)", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
                                                                                        
KODESeqOriginImage <- ggplot(KOOriginMerge, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM KO ZScore (1.5 Million Reads)", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
KODESeq3MilImage <- ggplot(KO3MilMerge, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM KO ZScore (3 Million Reads)", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 
KODESeq10MilImage <- ggplot(KO10MilMerge, aes(x=signedZ.DESeq.TotalSeq, y=signedZ.DESeq.QuantSeq)) + geom_point(alpha=.2) + 
  theme_bw() + labs(title="FUS ctrl vs HOM KO ZScore (10 Million Reads)", x="TotalRNASeq", y="QuantSeq") + theme(text = element_text(size=20)) +
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) 

d14DESeqOriginImage + d14DESeq3MilImage + d14DESeq10MilImage + KODESeqOriginImage + KODESeq3MilImage + KODESeq10MilImage

#Venn Diagrams
SigGenes <- function(DESeqResults){
  QuantDat <- as.data.frame(DESeqResults) %>%
    drop_na
  SigGenes <- filter(QuantDat, ADJ.PVAL.DESeq < 0.05)
  SigList <- SigGenes$gene_id
  return(SigList)
}

ListCompare <- function(Quant, Total){
  QuantSigList <- SigGenes(Quant)
  TotalSigList <- SigGenes(Total)
  ids = sort(unique(c(as.character(QuantSigList),  as.character(TotalSigList) )))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% QuantSigList
    counts[i, 2] = ids[i] %in% TotalSigList
  }
  colnames(counts) = c("QuantSeq", "RNA-Seq")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

d14OriginMerge <- ListCompare(d14QuantV1DeNoFilt, d14OriginTotalDeNoFilt)
d143MilMerge <- ListCompare(d14QuantV1DeNoFilt, d143MilTotalDeNoFilt)
d1410MilMerge <- ListCompare(d14QuantV1DeNoFilt, d1410MilTotalDeNoFilt)
d1415MilMerge <- ListCompare(d14QuantV1DeNoFilt, d1415MilTotalDeNoFilt)
d1420MilMerge <- ListCompare(d14QuantV1DeNoFilt, d1420MilTotalDeNoFilt)
KOOriginMerge <- ListCompare(KOQuantV1DeNoFilt, KOOriginTotalDeNoFilt)
KO3MilMerge <- ListCompare(KOQuantV1DeNoFilt, KO3MilTotalDeNoFilt)
KO10MilMerge <- ListCompare(KOQuantV1DeNoFilt, KO10MilTotalDeNoFilt)
KO15MilMerge <- ListCompare(KOQuantV1DeNoFilt, KO15MilTotalDeNoFilt)
KO20MilMerge <- ListCompare(KOQuantV1DeNoFilt, KO20MilTotalDeNoFilt)


par(mfrow=c(2,3))
vennDiagram(d14OriginMerge) + title("1.5 million sample d14 Venn")
vennDiagram(d143MilMerge) + title("3 million sample d14 Venn")
vennDiagram(d1410MilMerge) + title("10 million sample d14 Venn")
vennDiagram(KOOriginMerge) + title("1.5 million sample KO Venn")
vennDiagram(KO3MilMerge) + title("3 million sample KO Venn")
vennDiagram(KO10MilMerge) + title("10 million sample KO Venn")


SigGenes <- function(DESeqResults){
  QuantDat <- as.data.frame(DESeqResults) %>%
    drop_na
  SigGenes <- filter(QuantDat, ADJ.PVAL.DESeq < 0.05)
  SigList <- SigGenes$gene_id
  return(SigList)
}

ListCompareTotal <- function(Total1Mil, Total3Mil, Total10Mil, Total15Mil, Total20Mil){
  Total1MilSigList <- SigGenes(Total1Mil)
  Total3MilSigList <- SigGenes(Total3Mil)
  Total10MilSigList <- SigGenes(Total10Mil)
  Total15MilSigList <- SigGenes(Total15Mil)
  Total20MilSigList <- SigGenes(Total20Mil)
  ids = sort(unique(c(as.character(Total1MilSigList),  as.character(Total3MilSigList), as.character(Total10MilSigList), 
                      as.character(Total15MilSigList), as.character(Total20MilSigList) )))
  counts = matrix(0, nrow=length(ids), ncol=5)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% Total1MilSigList
    counts[i, 2] = ids[i] %in% Total3MilSigList
    counts[i, 3] = ids[i] %in% Total10MilSigList
    counts[i, 4] = ids[i] %in% Total15MilSigList
    counts[i, 5] = ids[i] %in% Total20MilSigList
  }
  colnames(counts) = c("Total1Mil", "Total3Mil", "Total10Mil", "Total15Mil", "Total20Mil")
  row.names(counts) = ids
  return(as.data.frame(counts))
}

KOMerge <-ListCompareTotal(KOOriginTotalDeNoFilt, KO3MilTotalDeNoFilt, KO10MilTotalDeNoFilt, KO15MilTotalDeNoFilt, KO20MilTotalDeNoFilt)

KOMergeSums <- t(colSums(KOMerge))
Size <- format(c(1500000, 3000000, 10000000, 15000000, 20000000), scientific = FALSE)
KOMergeMat <- rbind(as.numeric(KOMergeSums), as.numeric(Size))
rownames(KOMergeMat) <- colnames(colSums(KOMerge))
KOMergeDF <- as.data.frame(t(KOMergeMat))
colnames(KOMergeDF) <- c("SignificantGenes", "SampleSize")
KOMergeDF <- arrange(KOMergeDF, SampleSize)

ggplot(data=KOMergeDF, aes(y=SignificantGenes, x=SampleSize, group = 1)) +
  geom_line() +
  geom_point()

d14Merge <-ListCompareTotal(d14OriginTotalDeNoFilt, d143MilTotalDeNoFilt, d1410MilTotalDeNoFilt, d1415MilTotalDeNoFilt, d1420MilTotalDeNoFilt)
d14MergeSums <- t(colSums(d14Merge))
Size <- format(c(1500000, 3000000, 10000000, 15000000, 20000000), scientific = FALSE)
d14MergeMat <- rbind(as.numeric(d14MergeSums), as.numeric(Size))
rownames(d14MergeMat) <- colnames(colSums(d14Merge))
d14MergeDF <- as.data.frame(t(d14MergeMat))
colnames(d14MergeDF) <- c("SignificantGenes", "SampleSize")
d14MergeDF <- arrange(d14MergeDF, SampleSize)

ggplot(data=d14MergeDF, aes(y=SignificantGenes, x=SampleSize, group = 1)) +
  geom_line() +
  geom_point()
