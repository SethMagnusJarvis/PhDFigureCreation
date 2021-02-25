library(tidyverse)
library(pheatmap)
library(ggplot2)

F210ISupport <- read_csv("FTDSupportAS.csv")
BrainType <- read_csv("FTDBrainBIAB.csv")  %>% column_to_rownames("CellType") %>% 
  dplyr::select(F210ISupport$Sample) 
colnames(BrainType)[1:ncol(BrainType)] <- F210ISupport$UpdatedName

svg("BIABCellTypeHeatmap.svg")
pheatmap::pheatmap(BrainType, scale="row")
dev.off()

PCAType <- prcomp(t(BrainType), scale=TRUE)
PCARes <- as.data.frame(PCAType$x) %>% 
  rownames_to_column("UpdatedName") %>%
  left_join(F210ISupport, by="UpdatedName")

PCACellType <- ggplot(PCARes, aes(x=PC1,y=PC2)) + geom_point(aes(color=condition_A), size=4) +
  theme_bw(base_size=32) + 
  theme(legend.position="top") + geom_text(aes(label=UpdatedName),hjust=0, vjust=0) + theme(
    legend.title = element_text(color = "blue", size = 14),
    legend.text = element_text(color = "red", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

svg("BIABCellTypePCA.svg")
PCACellType
dev.off()