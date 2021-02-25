library(tidyverse)
library(ggplot2)
library(rgl)
library(pheatmap)
library(limma)

#import data
HETDE <- read_csv("F210IAdultDE.csv") 
colnames(HETDE) <- c("EnsemblID", "log2FoldChange", "stat", "pvalue", "padj", "signedZ")
HETRPKM <-  read_csv("rpkm_values.csv")[,c(1,3:14)]

HOMDE <- read_tsv("deseq_F210I_HOM_embryo_June_2016_differential_expression.tab") 
HOMRPKM <-  read_csv("rpkm_values_F210I_HOM.csv")[,c(1,3:10)]

#merge tables and run PCA
RPKMJoin <- full_join(HETRPKM, HOMRPKM)
DEJoin <- full_join(HETDE, HOMDE, by="EnsemblID")
ExpressionPCA <- prcomp(t(RPKMJoin[,2:21]))$x
ExpressionPCA <- as.data.frame(ExpressionPCA) %>% rownames_to_column("Sample")

#make PCA plot
Cond <- list()
Cond[1:5] <- "WTHET"
Cond[6:12] <- "HET"
Cond[13:16] <- "WTHOM"
Cond[17:20] <- "HOM"

ExpressionPCA$Type <- unlist(Cond)

ggplot(ExpressionPCA, aes(x=PC1, y=PC2, shape=Type)) + geom_point()  +
  scale_shape_manual(values=c(1,2,3,4)) + geom_text(label = ExpressionPCA$Sample)

##########################################################################################################################

#make venn diagrams
SigGenes <- function(DESeqResults){
  QuantDat <- as.data.frame(DESeqResults) %>%
    drop_na
  SigGenes <- QuantDat[QuantDat$padj < 0.05,]
  SigList <- SigGenes$EnsemblID
  return(SigList)
}

ListCompare <- function(HET, HOM){
  HETSigList <- SigGenes(HET)
  HOMSigList <- SigGenes(HOM)
  ids = sort(unique(c(as.character(HETSigList),  as.character(HOMSigList) )))
  counts = matrix(0, nrow=length(ids), ncol=2)
  for(i in 1:length(ids)){
    counts[i, 1] = ids[i] %in% HETSigList
    counts[i, 2] = ids[i] %in% HOMSigList
  }
  colnames(counts) = c("HET", "HOM")
  row.names(counts) = ids
  
  return(as.data.frame(counts))
}

HetVHOMDE <- ListCompare(HETDE, HOMDE)

svg("HETvHOMVennUpdate.svg")
vennDiagram(HetVHOMDE) + title("Comparison of DE of HET and HOM")
dev.off()

###########################################################################################################################
#Colour Z-score
Diff <- HETDE
MakeSignedZscore <- function(Diff){
  signedZ <- Diff %>%
    mutate(ZScore = ifelse(test = Diff$log2FoldChange > 0,
                           yes = qnorm(1 - (Diff$pvalue / 2)),
                           no = qnorm(Diff$pvalue / 2) )) %>%
    na.omit
  # high pass - otherwise very low P is set to Inf
  signedZ$ZScore[signedZ$ZScore > 20] <- 20
  signedZ$ZScore[signedZ$ZScore < -20] <- -20
  signedZ <- dplyr::select(signedZ, EnsemblID, ZScore)
  return(signedZ)
}

produceZscoreTable <- function(HETSEQ, HOMSEQ){
  HETZ <- MakeSignedZscore(HETSEQ)
  HOMZ <- MakeSignedZscore(HOMSEQ)
  Zscores <- full_join(HOMZ, HETZ, by = "EnsemblID", suffix = c(".HOM", ".HET")) %>%
    na.omit %>%
    column_to_rownames("EnsemblID")
  return(Zscores)
}

ZScore <- produceZscoreTable(HETDE, HOMDE)


#Z-score plot 2
ZScoreColour <- drop_na(ZScore) %>% mutate(
  colour = case_when(
    ZScore.HET > 3 & ZScore.HOM > 3 ~"#1b9e77",
    ZScore.HET < -3 & ZScore.HOM < -3 ~"#1b9e77",
    ZScore.HET < -3 & ZScore.HOM > 3 ~"#e7298a",
    ZScore.HET > 3 & ZScore.HOM < -3 ~"#e7298a",
    ZScore.HET > 3 ~"#7570b3",
    ZScore.HET < -3 ~"#7570b3",
    ZScore.HOM > 3 ~"#d95f02",
    ZScore.HOM < -3 ~"#d95f02",
    ZScore.HET > -3 & ZScore.HET < 3 & ZScore.HOM > -3 & ZScore.HOM < 3 ~"#999999",
  ))


ScorePlotFilt <- ggplot(ZScoreColour, aes(x=ZScore.HOM, y=ZScore.HET, color=colour)) + geom_point(alpha=1) + theme_bw() + 
  labs(title="HET vs HOM F210I Mutant differential expression Zscore plots", x="HOM", y="HET") +  
  scale_colour_identity() +
  geom_vline(xintercept = c(-3,3)) + geom_hline(yintercept = c(-3,3)) +
  theme(text = element_text(size=20)) 

svg("ColouredZscorePlotUpdated.svg")
ScorePlotFilt
dev.off()
