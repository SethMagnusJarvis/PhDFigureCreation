library(tidyverse)
library(biomaRt)

#Retrieve all of the exon count files and combine them into a single file
filenames <- list.files("Reads", pattern="*_keep_dups.txt", full.names=TRUE)
ldf <- lapply(filenames, read.table)
Merged <- ldf %>% purrr::reduce(full_join, by = "V1")

names(Merged)[1] <- "EnsemblID"
names(Merged)[2] <- substr(filenames[1], 7, 10)
names(Merged)[3:ncol(Merged)] <- substr(filenames[2:length(filenames)], 7, 9)

#group all of the reads for different exons into a single file
MergedReads <- Merged %>% 
  separate(EnsemblID, c("EnsemblID", NA), sep = ':') %>%
  group_by(EnsemblID) %>%
  summarise_all(funs(sum))

write_csv(MergedReads, "MergedReads.csv")

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "hsapiens_gene_ensembl", 
              host = "www.ensembl.org")

#Set columns you want
attributes <- c("ensembl_gene_id", "external_gene_name")
#Set column to filter by
filters <- "ensembl_gene_id"

BMOutput <- getBM(attributes = attributes, filters = filters, values=MergedReads$EnsemblID, mart=mart)

write_csv(BMOutput, "FTDBMOutput.csv")