library(tidyverse)
library(DropletUtils)


d14TotalRaw <- read_csv("d14RawCounts.csv")
KOTotalRaw <- read_csv("KORawCounts.csv")

colnames(d14TotalRaw)[1] <- "gene_id"
colnames(KOTotalRaw)[1]  <- "gene_id"

#Created a function to sample
poplist <- function(Col, SampleSize){
  #Initial setup
  tally <- 1
  ReadCol <- list()
  for(x in 1:nrow(Col)) {
    #Create a list where each geneID is present once per read
    ReadCol[tally:as.integer((unlist(as.integer(Col[x,2]))+tally-1))] <- list(as.character(Col[x,1]))
    tally <- tally+unlist(Col[x,2])
  }
  #Take 1.5 million random rows from the matrix and unlist because of course
  ReadColSample <- sample(ReadCol, SampleSize)
  ReadColSample <- unlist(ReadColSample)
  #Return the number of times each value appeared
  ReadColSample <- table(ReadColSample)
  return(ReadColSample)
}

#A function to allow it to be easily run by just entering a matrix
FullThing <- function(Reads, SampleSize){
  #run sampling on the first column of reads
  listy <- as.data.frame(poplist(Reads[,1:2], SampleSize))
  #Change column names
  colnames(listy) <- colnames(Reads[,1:2])
  #Run sampling on the 2nd column of reads
  listy2 <- as.data.frame(poplist(Reads[,c(1,3)], SampleSize))
  #Change column names
  colnames(listy2) <- colnames(Reads[,c(1,3)])
  
  #Join the two previous datasets together so there's something to join other things to
  Join <- full_join(listy, listy2, by = "gene_id")
  Join2 <- Join
  
  for(x in 4:ncol(Reads)){
    #For each column of reads, sample it then attach it to Join2
    RandCounts <- as.data.frame(poplist(Reads[,c(1,x)], SampleSize))
    colnames(RandCounts) <- colnames(Reads)[c(1,x)]
    Join2 <- full_join(Join2, RandCounts, by = "gene_id")
  }
  return(Join2)
}

set.seed(16)

d14SampledOrigin <- FullThing(d14TotalRaw, 1500000)
KOSampledOrigin <- FullThing(KOTotalRaw,  1500000)
write_csv(KOSampledOrigin, "KOTotalRawSampledMk2.csv")
write_csv(d14SampledOrigin, "d14TotalRawSampledMk2.csv")

d14Sampled3mil <- FullThing(d14TotalRaw, 3000000)
KOSampled3mil <- FullThing(KOTotalRaw,  3000000)
write_csv(KOSampled3mil, "KOTotalRawSampled3MilMk1.csv")
write_csv(d14Sampled3mil, "d14TotalRawSampled3MilMk1.csv")

d14Sampled10mil <- FullThing(d14TotalRaw, 10000000)
KOSampled10mil <- FullThing(KOTotalRaw,  10000000)
write_csv(KOSampled10mil, "KOTotalRawSampled10MilMK1.csv")
write_csv(d14Sampled10mil, "d14TotalRawSampled10MilMK1.csv")

d14Sampled15mil <- FullThing(d14TotalRaw, 15000000)
KOSampled15mil <- FullThing(KOTotalRaw,  15000000)
write_csv(KOSampled15mil, "KOTotalRawSampled15MilMK1.csv")
write_csv(d14Sampled15mil, "d14TotalRawSampled15MilMK1.csv")

d14Sampled20mil <- FullThing(d14TotalRaw, 20000000)
KOSampled20mil <- FullThing(KOTotalRaw,  20000000)
write_csv(KOSampled20mil, "KOTotalRawSampled20MilMK1.csv")
write_csv(d14Sampled20mil, "d14TotalRawSampled20MilMK1.csv")
