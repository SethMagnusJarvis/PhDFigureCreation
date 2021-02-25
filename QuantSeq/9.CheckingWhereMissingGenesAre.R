library(tidyverse)

d14UnsampledDE <- read_csv("DESeqD14RNASeqUnsampled.csv")
KOUnsampledDE <- read_csv("DESeqKORNASeqUnsampled.csv")
d14QuantDE <- read_csv("DESeqD14QuantSeq.csv")
KOQuantDE <- read_csv("DESeqKOQuantSeq.csv")

d14RNAOnlySig <- anti_join(d14UnsampledDE, d14QuantDE, by = "EnsemblID") %>% drop_na() %>% filter(padj < 0.05)
d14QuantOnlySig <- anti_join(d14QuantDE, d14UnsampledDE, by = "EnsemblID") %>% drop_na() %>% filter(padj < 0.05)

KORNAOnlySig <- anti_join(KOUnsampledDE, KOQuantDE, by = "EnsemblID") %>% drop_na() %>% filter(padj < 0.05)
KOQuantOnlySig <- anti_join(KOQuantDE, KOUnsampledDE, by = "EnsemblID") %>% drop_na() %>% filter(padj < 0.05)

d14QuantRaw <- as.data.frame(read.table("jack4.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
d14QuantSeq <- select(d14QuantRaw, "gene_id", "c1.d14_WT1":"t4.d14_HOM4")

KOQuantRaw <- as.data.frame(read.table("jack2.expression_genes.tab", sep="\t",header=TRUE, stringsAsFactors = FALSE ))
KOQuantSeq <- select(KOQuantRaw, "gene_id", "c1.KO_WT1":"t4.KO_HOM4")

d14Join <- left_join(d14RNAOnlySig, d14QuantSeq, by = c("EnsemblID" = "gene_id"))
KOJoin <- inner_join(KORNAOnlySig, KOQuantSeq, by = c("EnsemblID" = "gene_id"))
