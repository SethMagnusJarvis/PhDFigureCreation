library(tidyverse)
library(EnhancedVolcano)

load("D:/GoogleDrive/Work/UCL/F210IAnalysis/sgseq/F210I_adult_sc_res_clean_novel.RData")

DefaultVolcano <- EnhancedVolcano(res.clean,
                                   lab = res.clean$geneName,
                                   x = 'log2fold_HET_CONTROL',
                                   y = 'pvalue',
                                   title = 'F210I Volcano')

LessCutoffVolcano <- EnhancedVolcano(res.clean,
                                   lab = res.clean$geneName,
                                   x = 'log2fold_HET_CONTROL',
                                   y = 'padj',
                                   title = 'F210I Volcano',
                                   pCutoff=0.1,
                                   FCcutoff=0)

#plain
HETSplicing <- res.clean %>% drop_na() %>% rowwise() %>%
  mutate(MeanWTPSI = mean(c(F210I_WT_1_psi:F210I_WT_5_psi)), MeanHETPSI = mean(c(F210I_HET_1_psi:F210I_HET_7_psi))) %>% 
  mutate(PSIChange = MeanWTPSI-MeanHETPSI) %>%
  filter(!geneName =="Sh3bgr") 

DefaultVolcano <- EnhancedVolcano(HETSplicing,
                                  lab = HETSplicing$geneName,
                                  x = 'PSIChange',
                                  y = 'pvalue',
                                  title = 'F210I PSI Volcano',
                                  FCcutoff=0)

LessCutoffVolcano001 <- EnhancedVolcano(HETSplicing,
                                     lab = HETSplicing$geneName,
                                     x = 'PSIChange',
                                     y = 'pvalue',
                                     title = 'F210I PSI Volcano 0.001 Cutoff',
                                     pCutoff=0.001,
                                     FCcutoff=0)

LessCutoffVolcano0001 <- EnhancedVolcano(HETSplicing,
                                     lab = HETSplicing$geneName,
                                     x = 'PSIChange',
                                     y = 'pvalue',
                                     title = 'F210I PSI Volcano 0.0001 Cutoff',
                                     pCutoff=0.0001,
                                     FCcutoff=0)

#One Per group, and Sh3bgr capped
GroupedSplicing <- res.clean %>% drop_na() %>% 
  group_by(groupID) %>% do(head(.,1)) %>%
  ungroup() %>% rowwise() %>%
  mutate(MeanWTPSI = mean(c(F210I_WT_1_psi:F210I_WT_5_psi)), MeanHETPSI = mean(c(F210I_HET_1_psi:F210I_HET_7_psi))) %>% 
  mutate(PSIChange = MeanWTPSI-MeanHETPSI)

GroupedSplicing$pvalue[GroupedSplicing$pvalue < 1e-20] <- 1e-20

DefaultVolcano <- EnhancedVolcano(GroupedSplicing,
                                  lab = GroupedSplicing$geneName,
                                  x = 'PSIChange',
                                  y = 'pvalue',
                                  title = 'F210I PSI Volcano',
                                  FCcutoff=0)

LessCutoffVolcano001 <- EnhancedVolcano(GroupedSplicing,
                                        lab = GroupedSplicing$geneName,
                                        x = 'PSIChange',
                                        y = 'pvalue',
                                        title = 'F210I PSI Volcano 0.001 Cutoff',
                                        pCutoff=0.001,
                                        FCcutoff=0)

LessCutoffVolcano0001 <- EnhancedVolcano(GroupedSplicing,
                                         lab = GroupedSplicing$geneName,
                                         x = 'PSIChange',
                                         y = 'pvalue',
                                         title = 'F210I PSI Volcano 0.0001 Cutoff',
                                         pCutoff=0.0001,
                                         FCcutoff=0)

svg("PSIUncappedVolcano.svg")
DefaultVolcano
dev.off()

svg("PSIUncappedVolcano001.svg")
LessCutoffVolcano001
dev.off()

svg("PSIUncappedVolcano0001.svg")
LessCutoffVolcano0001
dev.off()

#Expression Capped
GroupedSplicing <- res.clean %>% drop_na() %>% 
  group_by(groupID) %>% do(head(.,1)) %>%
  ungroup() %>% rowwise() %>%
  mutate(MeanWTPSI = mean(c(F210I_WT_1_psi:F210I_WT_5_psi)), MeanHETPSI = mean(c(F210I_HET_1_psi:F210I_HET_7_psi))) %>% 
  mutate(PSIChange = MeanWTPSI-MeanHETPSI)

GroupedSplicing$pvalue[GroupedSplicing$pvalue < 1e-20] <- 1e-20

MeanGroupedSplicing <- filter(GroupedSplicing, exonBaseMean > 100)

DefaultVolcano <- EnhancedVolcano(MeanGroupedSplicing,
                                  lab = MeanGroupedSplicing$geneName,
                                  x = 'PSIChange',
                                  y = 'pvalue',
                                  title = 'F210I PSI Volcano',
                                  FCcutoff=0)

LessCutoffVolcano001 <- EnhancedVolcano(MeanGroupedSplicing,
                                        lab = MeanGroupedSplicing$geneName,
                                        x = 'PSIChange',
                                        y = 'pvalue',
                                        title = 'F210I PSI Volcano 0.001 Cutoff',
                                        pCutoff=0.001,
                                        FCcutoff=0)

LessCutoffVolcano0001 <- EnhancedVolcano(MeanGroupedSplicing,
                                         lab = MeanGroupedSplicing$geneName,
                                         x = 'PSIChange',
                                         y = 'pvalue',
                                         title = 'F210I PSI Volcano 0.0001 Cutoff',
                                         pCutoff=0.0001,
                                         FCcutoff=0)

LessCutoffVolcano00001 <- EnhancedVolcano(MeanGroupedSplicing,
                                         lab = MeanGroupedSplicing$geneName,
                                         x = 'PSIChange',
                                         y = 'pvalue',
                                         title = 'F210I PSI Volcano 0.0001 Cutoff',
                                         pCutoff=0.00001,
                                         FCcutoff=0)

svg("PSICappedVolcano.svg")
DefaultVolcano
dev.off()

svg("PSICappedVolcano001.svg")
LessCutoffVolcano001
dev.off()

svg("PSICappedVolcano0001.svg")
LessCutoffVolcano0001
dev.off()

svg("PSICappedVolcano00001.svg")
LessCutoffVolcano00001
dev.off()

ExpressionPSIPlot <- ggplot(GroupedSplicing, aes(x=log(exonBaseMean), y=PSIChange)) + geom_point(alpha=.2) + theme_bw() + 
  labs(title="Reads VS PSI change", x="log Mean Reads", y="PSI Change") +
  theme(text = element_text(size=20)) 


svg("ExpressionVSPSI.svg")
ExpressionPSIPlot
dev.off()
