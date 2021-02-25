library(tidyverse)
library(ggplot2)
library(patchwork)

d14Unsampled <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampled <- read_csv("DESeqKORNASeqUnsampled.csv")
d14Quant <- read_csv("DESeqD14QuantSeq.csv")
KOQuant <- read_csv("DESeqKOQuantSeq.csv")


d14Join <- full_join(d14Unsampled, d14Quant, by = "EnsemblID", suffix = c(".RNASeq", ".Quant")) %>%
  na.omit %>% 
  mutate(rank.RNASeq = rank(desc(baseMean.RNASeq)), rank.Quant = rank(desc(baseMean.Quant)))

KOJoin <- full_join(KOUnsampled, KOQuant, by = "EnsemblID", suffix = c(".RNASeq", ".Quant")) %>%
  na.omit %>% 
  mutate(rank.RNASeq = rank(desc(baseMean.RNASeq)), rank.Quant = rank(desc(baseMean.Quant)))

AdjSigPlot <- function(Joined, Filter){
  Joined$col <- "black"
  Joined[ Joined$padj.RNASeq <= 0.05 , "col"] <- "green"
  Joined[ Joined$padj.Quant <= 0.05 , "col"] <- "blue"
  Joined[ Joined$padj.RNASeq < 0.05 & Joined$padj.Quant < 0.05, "col"] <- "red"
  
  if(Filter == TRUE){
    Joined <- filter(Joined, col != "black")
  }
  
  Plot <- ggplot(Joined, aes(x=log(baseMean.RNASeq), y=log(baseMean.Quant), color = col)) + geom_point() + 
    theme_bw() +  scale_color_identity() + 
    coord_cartesian(xlim = c(2, 10), ylim = c(2,10))+
    geom_smooth(method = "lm", se = FALSE)
  return(Plot)
}

UnadjSigPlot <- function(Joined, Filter){
  Joined$col <- "black"
  Joined[ Joined$padj.RNASeq <= 0.05 , "col"] <- "green"
  Joined[ Joined$padj.Quant <= 0.05 , "col"] <- "blue"
  Joined[ Joined$padj.RNASeq < 0.05 & Joined$pvalue.Quant < 0.05, "col" ] <- "red"
  Joined[ Joined$pvalue.RNASeq < 0.05 & Joined$padj.Quant < 0.05, "col" ] <- "red"
  
  if(Filter == TRUE){
    Joined <- filter(Joined, col != "black")
  }
  
  Plot <- ggplot(Joined, aes(x=log(baseMean.RNASeq), y=log(baseMean.Quant), color = col)) + geom_point() + 
    theme_bw() +  scale_color_identity()+ 
    coord_cartesian(xlim = c(2, 10), ylim = c(2,10)) +
    geom_smooth(method = "lm", se = FALSE)
  return(Plot)
}

d14Sig <- AdjSigPlot(d14Join, FALSE)
d14FiltSig <- AdjSigPlot(d14Join, TRUE)
d14UnSig <- UnadjSigPlot(d14Join, FALSE)
d14FiltUnSig <- UnadjSigPlot(d14Join, TRUE)

KOSig <- AdjSigPlot(KOJoin, FALSE)
KOFiltSig <- AdjSigPlot(KOJoin, TRUE)
KOUnSig <- UnadjSigPlot(KOJoin, FALSE)
KOFiltUnSig <- UnadjSigPlot(KOJoin, TRUE)

d14Sig+d14FiltSig+d14UnSig+d14FiltUnSig
KOSig+KOFiltSig+KOUnSig+KOFiltUnSig

#Change rank from descending to ascending
d14Join <- full_join(d14Unsampled, d14Quant, by = "EnsemblID", suffix = c(".RNASeq", ".Quant")) %>%
  na.omit %>% 
  mutate(rank.RNASeq = rank(baseMean.RNASeq), rank.Quant = rank(baseMean.Quant))

KOJoin <- full_join(KOUnsampled, KOQuant, by = "EnsemblID", suffix = c(".RNASeq", ".Quant")) %>%
  na.omit %>% 
  mutate(rank.RNASeq = rank(baseMean.RNASeq), rank.Quant = rank(baseMean.Quant))

d14Sig <- AdjSigPlot(d14Join, FALSE)
d14FiltSig <- AdjSigPlot(d14Join, TRUE)
d14UnSig <- UnadjSigPlot(d14Join, FALSE)
d14FiltUnSig <- UnadjSigPlot(d14Join, TRUE)

KOSig <- AdjSigPlot(KOJoin, FALSE)
KOFiltSig <- AdjSigPlot(KOJoin, TRUE)
KOUnSig <- UnadjSigPlot(KOJoin, FALSE)
KOFiltUnSig <- UnadjSigPlot(KOJoin, TRUE)

d14Sig+d14FiltSig+d14UnSig+d14FiltUnSig
KOSig+KOFiltSig+KOUnSig+KOFiltUnSig
