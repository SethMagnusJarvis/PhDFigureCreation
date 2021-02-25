library(tidyverse)
library(ggplot2)

setwd("D:/GoogleDrive/Work/UCL/F210IAnalysis/sgseq/")

RESCleanHOM <- read_tsv("F210I_HOM_embryonic_brain_res_clean_novel.tab") %>% dplyr::rename("ensemblName"="geneName")

RESCleanHET <- read_tsv("F210I_HET_adult_brain_res_clean_novel.tab")

merged <- full_join(RESCleanHET,  RESCleanHOM,  by=c("from"="from", "to"="to", "type"="type")) %>%
  drop_na()

############################################
#Splicing Z-score plot

#Make Z score
MakeSignedZscoreGroup <- function(Diff){
  signedZ <- Diff %>%
    mutate(ZScore = ifelse(test = Diff[,10] > 0,
                           yes = qnorm(1 - (Diff$pvalue / 2)),
                           no = qnorm(Diff$pvalue / 2) )) %>%
    na.omit
  # high pass - otherwise very low P is set to Inf
  signedZ$ZScore[signedZ$ZScore > 20] <- 20
  signedZ$ZScore[signedZ$ZScore < -20] <- -20
  signedZ <- dplyr::select(signedZ, ZScore, from, to, type, ensemblName, groupID)
  return(signedZ)
}

produceZscoreTableGroup <- function(HETSEQ, HOMSEQ){
  HETZ <- MakeSignedZscoreGroup(HETSEQ)
  HOMZ <- MakeSignedZscoreGroup(HOMSEQ)
  Zscores <- full_join(HOMZ, HETZ, by = c("from"="from", "to"="to", "type"="type", "ensemblName"="ensemblName"), suffix = c(".HOM", ".HET")) 
  return(Zscores)
}

#Produce HetVHom Z score table
#Convert Z.score columns to vectors
GroupedZScore <- produceZscoreTableGroup(RESCleanHET, RESCleanHOM) 
GroupedZScore$ZScore.HOM <- as.vector(GroupedZScore$ZScore.HOM)
GroupedZScore$ZScore.HET <- as.vector(GroupedZScore$ZScore.HET)
#Group Z-scores and only select one per group
OneGroupZScore <- arrange(GroupedZScore, desc(abs(unlist(ZScore.HET)))) %>% 
  group_by(groupID.HOM) %>% do(head(.,1)) %>%
  ungroup()


#Create column of colours based on Z-scores
ZScoreColourGrouped <- drop_na(OneGroupZScore) %>% mutate(
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

#Select the relevant columns from Rescleanhet and join to the grouped Z-score
ZScoreColourGrouped <- dplyr::select(RESCleanHET, from, to, type, geneName)  %>% 
  right_join(ZScoreColourGrouped, by=c("from"="from", "to"="to", "type"="type"))

#Produce coloured z-score plot
ZScorePlotGrouped <- ggplot(ZScoreColourGrouped, aes(x=ZScore.HOM, y=ZScore.HET, color=colour)) + geom_point(alpha=1) + theme_bw() + 
  labs(title="HET vs HOM F210I Mutant splicing Zscore plots", x="HOM", y="HET") +  
  scale_colour_identity() +
  geom_vline(xintercept = c(-3,3)) + geom_hline(yintercept = c(-3,3)) +
  theme(text = element_text(size=20)) + 
  geom_text(aes(label=ifelse(abs(ZScore.HOM)>3 & abs(ZScore.HET)>3,as.character(geneName),'')),hjust=0.2,vjust=1) +
  coord_cartesian(xlim = c(-20, 20), ylim = c(-20, 20))


svg("SplicingZScoreHETvHOM.svg")
ZScorePlotGrouped
dev.off


#Bar plots
#Make columns that show which fold change and PSI are larger in HET and HOM, and whether they go up or down.
merged <- full_join(RESCleanHET,  RESCleanHOM,  by=c("from"="from", "to"="to", "type"="type"),  suffix = c(".HOM", ".HET")) %>%
  drop_na() %>% 
  group_by(groupID.HET) %>% do(head(.,1)) %>%
  ungroup() %>% rowwise() %>%
  mutate(MeanWTPSIHET = mean(c(F210I_WT_1_psi:F210I_WT_5_psi)), MeanHETPSI = mean(c(F210I_HET_1_psi:F210I_HET_7_psi))) %>% 
  mutate(PSIChangeHET = MeanWTPSIHET-MeanHETPSI) %>% 
  mutate(MeanWTPSIHOM = mean(c(F210I_CTL_1_norm_psi:F210I_CTL_4_norm_psi)), MeanHOMPSI = mean(c(F210I_HOM_1_norm_psi:F210I_HOM_4_norm_psi))) %>% 
  mutate(PSIChangeHOM = MeanWTPSIHOM-MeanHOMPSI) %>%
  mutate(PSIRatio = PSIChangeHOM/PSIChangeHET) %>%
  mutate(FoldChangeRatio = log2fold_HOM_CTL/log2fold_HET_CONTROL) %>%
  ungroup() %>%
  mutate(HOMFoldSign = case_when(log2fold_HOM_CTL<0 ~ "Negative", log2fold_HOM_CTL>0 ~ "Positive")) %>%
  mutate(HETFoldSign = case_when(log2fold_HET_CONTROL<0 ~ "Negative", log2fold_HET_CONTROL>0 ~ "Positive")) %>%
  mutate(HOMPSISign = case_when(PSIChangeHOM<0 ~ "Negative", PSIChangeHOM>0 ~ "Positive")) %>%
  mutate(HETPSISign = case_when(PSIChangeHET<0 ~ "Negative", PSIChangeHET>0 ~ "Positive")) %>%
  mutate(LargerPSI = case_when(abs(PSIChangeHET)<abs(PSIChangeHOM) ~ "HOM", abs(PSIChangeHET)>abs(PSIChangeHOM) ~ "HET")) %>%
  mutate(LargerFold = case_when(abs(log2fold_HET_CONTROL)<abs(log2fold_HOM_CTL) ~ "HOM", abs(log2fold_HET_CONTROL)>abs(log2fold_HOM_CTL) ~ "HET")) %>%
  mutate(FCComp = case_when(
    LargerFold == "HOM" & HOMFoldSign =="Negative" ~ "-HOM",
    LargerFold == "HOM" & HOMFoldSign =="Positive" ~ "+HOM",
    LargerFold == "HET" & HETFoldSign =="Negative" ~ "-HET",
    LargerFold == "HET" & HETFoldSign =="Positive" ~ "+HET",
  )) %>%
  mutate(PSIComp = case_when(
    LargerPSI == "HOM" & HOMPSISign =="Negative" ~ "-HOM",
    LargerPSI == "HOM" & HOMPSISign =="Positive" ~ "+HOM",
    LargerPSI == "HET" & HETPSISign =="Negative" ~ "-HET",
    LargerPSI == "HET" & HETPSISign =="Positive" ~ "+HET",
  )) %>% drop_na()

#Count the number of fold changes which fall into each catagory
FoldCounted<- merged %>% 
  group_by(FCComp) %>%
  summarise(total_genes = n()) 

ggplot(data=FoldCounted, aes(x=FCComp, y=total_genes)) +
  geom_bar(stat="identity") + labs(x = "Comparison of Splicing Fold Change", y="Total Number of Genes") + 
  theme(text = element_text(size = 20))

#Count the number of PSIs that fall into each category and plot
PSICounted <- merged %>% 
  group_by(PSIComp) %>%
  summarise(total_genes = n()) %>% drop_na()

ggplot(data=PSICounted, aes(x=PSIComp , y=total_genes)) +
  geom_bar(stat="identity") + labs(x = "Comparison of PSI ", y="Total Number of Genes") + 
  theme(text = element_text(size = 20))
