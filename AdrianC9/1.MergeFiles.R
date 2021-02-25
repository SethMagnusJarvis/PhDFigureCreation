library(tidyverse)

#Read in all files in the Counts folder that have the keep_dups suffix
filenames <- list.files("Counts", pattern="*_keep_dups.txt", full.names=TRUE)
#Combine all of them into one big list then merge that into a dataframe
ldf <- lapply(filenames, read.table)
Merged <- ldf %>% purrr::reduce(full_join, by = "V1")

#Give the columns more approproate names
names(Merged)[1] <- "EnsemblID"
names(Merged)[2:31] <- substr(filenames, 8, 20)

#Write list of all Exons
write_csv(Merged, "AdrianC9MergedExonsAll.csv")

#Group all exons by gene and sum their expression
MergedGrouped <- Merged %>% 
  separate(EnsemblID, c("EnsemblID", NA), sep = ':') %>%
  group_by(EnsemblID) %>%
  summarise_all(funs(sum))

#Write List of Genes
write_csv(MergedGrouped, "AdrianC9GeneReadsAll.csv")
