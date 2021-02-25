library(tidyverse)
library(biomaRt)

Reads <- read_csv("AdrianC9GeneReadsAll.csv")

mart <- useMart("ensembl")
mart <- useDataset("mmusculus_gene_ensembl", mart)
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "mmusculus_gene_ensembl", 
              host = "www.ensembl.org")


#Set columns you want
attributes <- c("ensembl_gene_id", "gene_biotype", "external_gene_name", 'entrezgene_id')
#Set column to filter by
filters <- "ensembl_gene_id"

BMOutput <- getBM(attributes = attributes, filters = filters, values=Reads$EnsemblID, mart=mart)

write_csv(BMOutput, "BioMarTMouseWEntrez.csv")