library(tidyverse)
library(ggplot2)

setwd("D:/GoogleDrive/Work/UCL/Seth/TestingQuantSeq/SigGenes/GoTermLists/Analysis")

Naming <- c("GO biological process", "Mouse #", "#", "expected",  "+/-", "Fold Enrichment", "Raw p-value", "FDR")

d14Quant <- read_tsv("d14QuantSig.txt", skip = 6)
colnames(d14Quant) <- Naming
d14Quant <- d14Quant %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
d14Quant$`GO Term` <- str_remove(d14Quant$`GO Term`, '[)]')
d14QuantSig <- d14Quant %>% filter(FDR < 0.05)

d14Unsampled <- read_tsv("d14UnsampledSig.txt", skip = 6)
colnames(d14Unsampled) <- Naming
d14Unsampled <- d14Unsampled %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
d14Unsampled$`GO Term` <- str_remove(d14Unsampled$`GO Term`, '[)]')
d14UnsampledSig <- d14Unsampled %>% filter(FDR < 0.05)

KOQuant <- read_tsv("KOQuantSig.txt", skip = 6)
colnames(KOQuant) <- Naming
KOQuant <- KOQuant %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
KOQuant$`GO Term` <- str_remove(KOQuant$`GO Term`, '[)]')
KOQuantSig <- KOQuant %>% filter(FDR < 0.05)

KOUnsampled <- read_tsv("KOUnsampledSig.txt", skip = 6)
colnames(KOUnsampled) <- Naming
KOUnsampled <- KOUnsampled %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
KOUnsampled$`GO Term` <- str_remove(KOUnsampled$`GO Term`, '[)]')
KOUnsampledSig <- KOUnsampled %>% filter(FDR < 0.05)

KOJoin <- full_join(KOUnsampled, KOQuant, by="GO Term", suffix=c(".Total", ".Quant")) %>%
  select("Description.Total", "GO Term", "+/-.Total", "Fold Enrichment.Total",
         "FDR.Total", "+/-.Quant", "Fold Enrichment.Quant", "FDR.Quant") %>%
  filter(FDR.Total < 0.05 | FDR.Quant < 0.05) %>%
  arrange(FDR.Total, FDR.Quant)
KOJoinBothSig <- KOJoin %>% filter(FDR.Total < 0.05 & FDR.Quant < 0.05)
d14Join <- full_join(d14Unsampled, d14Quant, by="GO Term", suffix=c(".Total", ".Quant")) %>%
  select("Description.Total", "GO Term", "+/-.Total", "Fold Enrichment.Total", 
         "FDR.Total", "+/-.Quant", "Fold Enrichment.Quant", "FDR.Quant") %>%
  filter(FDR.Total < 0.05 | FDR.Quant < 0.05) %>%
  arrange(FDR.Total, FDR.Quant)
d14JoinBothSig <- d14Join %>% filter(FDR.Total < 0.05 & FDR.Quant < 0.05)

write.csv(d14Join, "d14GOJoin.csv")
write.csv(KOJoin, "KOGOJoin.csv")


''' 
d14QuantAdjSig <- read_tsv("d14QuantAdjSig.txt", skip = 6)
colnames(d14QuantAdjSig) <- Naming
d14QuantAdjSig <- d14QuantAdjSig %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
d14QuantAdjSig$`GO Term` <- str_remove(d14QuantAdjSig$`GO Term`, '[)]')

d14UnsampledAdjSig <- read_tsv("d14UnsampledAdjSig.txt", skip = 6)
colnames(d14UnsampledAdjSig) <- Naming
d14UnsampledAdjSig <- d14UnsampledAdjSig %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
d14UnsampledAdjSig$`GO Term` <- str_remove(d14UnsampledAdjSig$`GO Term`, '[)]')

KOQuantAdjSig <- read_tsv("KOQuantAdjSig.txt", skip = 6)
colnames(KOQuantAdjSig) <- Naming
KOQuantAdjSig <- KOQuantAdjSig %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
KOQuantAdjSig$`GO Term` <- str_remove(KOQuantAdjSig$`GO Term`, '[)]')

KOUnsampledAdjSig <- read_tsv("KOUnsampledAdjSig.txt", skip = 6)
colnames(KOUnsampledAdjSig) <- Naming
KOUnsampledAdjSig <- KOUnsampledAdjSig %>% arrange((FDR)) %>%
  separate("GO biological process", c("Description", "GO Term"), sep = '[(]')
KOUnsampledAdjSig$`GO Term` <- str_remove(KOUnsampledAdjSig$`GO Term`, '[)]')
'''

