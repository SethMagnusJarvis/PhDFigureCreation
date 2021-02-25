library(tidyverse)
library(ggplot2)
library(patchwork)

MakeSignedZscore <- function(Diff){
  signedZ <- Diff %>%
    mutate(ZScore = ifelse(test = Diff$log2FoldChange > 0,
                           yes = qnorm(1 - (Diff$pvalue / 2)),
                           no = qnorm(Diff$pvalue / 2) )) %>%
    na.omit
  # high pass - otherwise very low P is set to Inf
  signedZ$ZScore[signedZ$ZScore > 6] <- 6
  signedZ$ZScore[signedZ$ZScore < -6] <- -6
  signedZ <- select(signedZ, EnsemblID, ZScore)
  return(signedZ)
}

produceZscoreTable <- function(QuantSEQ, TotalSEQ){
  QuantZ <- MakeSignedZscore(QuantSEQ)
  TotalZ <- MakeSignedZscore(TotalSEQ)
  Zscores <- full_join(TotalZ, QuantZ, by = "EnsemblID", suffix = c(".RNASeq", ".Quant")) %>%
    na.omit %>%
    column_to_rownames("EnsemblID")
  return(Zscores)
}

Clusters <- read_tsv("clusters(1).bed", col_names=FALSE)
ClustersTissue <- read_tsv("clusters_withTissueInfo.bed", col_names=FALSE)

colnames(Clusters) <- c("Chr", "ClusterStart", "ClusterEnd", "ClusterID", "NoSupportingProtocols", "Strand", "PolyAInfo", "EnsemblID" ) 
colnames(ClustersTissue) <- c("Chr", "ClusterStart", "ClusterEnd", "ClusterID", "NoSupportingProtocols", "Strand", "PolyAInfo", "EnsemblID", "Tissue" )

ClusterCounts <- Clusters%>%
  group_by(EnsemblID) %>% 
  summarise(polyASignals = n())

ClusterTissueCounts <- ClustersTissue%>%
  group_by(EnsemblID) %>% 
  summarise(polyASignals = n())

ClusterTissWCounts <- left_join(ClustersTissue, ClusterTissueCounts)  %>%
  select(EnsemblID, polyASignals, Tissue)

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv") %>%
  left_join(ClusterCounts) %>%
  filter(polyASignals == 1)
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv") %>%
  left_join(ClusterCounts) %>%
  filter(polyASignals == 1)
d14Quant <- read_csv("DESeqD14QuantSeq.csv") %>%
  left_join(ClusterCounts) %>%
  filter(polyASignals == 1)
KOQuant <- read_csv("DESeqKOQuantSeq.csv") %>%
  left_join(ClusterCounts) %>%
  filter(polyASignals == 1)

d14UnsampledZscore <- produceZscoreTable(d14Quant, d14Unsampled)
KOUnsampledZscore <- produceZscoreTable(KOQuant, KOUnsampled)

d14 <- ggplot(d14UnsampledZscore, aes(x=ZScore.RNASeq, y=ZScore.Quant)) + geom_point(alpha=.2) + theme_bw() + 
  labs(title="FUS ctrl vs HOM d14 ZScore Only 1 polyA", x="TotalRNASeq", y="QuantSeq") + 
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + 
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) +
  theme(text = element_text(size=20)) 

KO <- ggplot(KOUnsampledZscore, aes(x=ZScore.RNASeq, y=ZScore.Quant)) +   geom_point(alpha=.2) + theme_bw() + 
  labs(title="FUS ctrl vs HOM KO ZScore Only 1 polyA", x="TotalRNASeq", y="QuantSeq") + 
  geom_vline(xintercept = c(-2,2)) + geom_hline(yintercept = c(-2,2)) + 
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6,6)) +
  theme(text = element_text(size=20))

d14 + KO
 