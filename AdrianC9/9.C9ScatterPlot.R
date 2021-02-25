#Making Amazon
library(tidyverse)
library(ggplot2)

################################################################################################################
#Load reads and support
Reads <- read_csv("AdrianC9GeneReadsAll.csv")
Support <- read_csv("AdrianC9SupportUpdated.csv") %>% arrange(condition_A)

#Create a 6 and 12 month old support file
SixMonthSupport <- filter(Support, type_3 == "6M") %>% dplyr::select(Sample, condition_A)
TwelveMonthSupport <- filter(Support, type_3 == "12M")  %>% dplyr::select(Sample, condition_A)

#Select the correct samples from the reads and remove all genes other than C9orf72 in both
SixFiltReads <- dplyr::select(Reads, EnsemblID, SixMonthSupport$Sample) %>% filter(EnsemblID == "ENSMUSG00000028300")%>% 
  t() 
TwelveFiltReads <- dplyr::select(Reads, EnsemblID, TwelveMonthSupport$Sample) %>% filter(EnsemblID == "ENSMUSG00000028300") %>% 
  t()

#Format dataframe and rename before making sure it's numeric
SixFiltReadsC9 <- as.data.frame(SixFiltReads) %>% rownames_to_column("Sample") %>% 
  left_join(SixMonthSupport, by = c("Sample"="Sample")) %>% dplyr::select(Sample, V1, condition_A) %>%
  rename("ReadsC9" = "V1")
SixFiltReadsC9 <- SixFiltReadsC9[-1,]
SixFiltReadsC9$ReadsC9 <- as.numeric(levels(SixFiltReadsC9$ReadsC9))[SixFiltReadsC9$ReadsC9]

TwelveFiltReadsC9 <- as.data.frame(TwelveFiltReads) %>% rownames_to_column("Sample") %>% 
  left_join(TwelveMonthSupport, by = c("Sample"="Sample")) %>% dplyr::select(Sample, V1, condition_A) %>%
  rename("ReadsC9" = "V1")
TwelveFiltReadsC9 <- TwelveFiltReadsC9[-1,]
TwelveFiltReadsC9$ReadsC9 <- as.numeric(levels(TwelveFiltReadsC9$ReadsC9))[TwelveFiltReadsC9$ReadsC9]

#Plot both scatter plots
SixFiltScatter<-ggplot(SixFiltReadsC9, aes(x=condition_A, y=ReadsC9)) + 
  geom_dotplot(binaxis='y', stackdir='center')

TwelveFiltScatter<-ggplot(TwelveFiltReadsC9, aes(x=condition_A, y=ReadsC9)) + 
  geom_dotplot(binaxis='y', stackdir='center')

