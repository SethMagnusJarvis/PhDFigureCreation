library(tidyverse)
library(EnrichmentBrowser)
library(EnhancedVolcano)
library(patchwork)


############################################################################################################################
#Import data
PipelineDE <- read_tsv("deseq_F210I_adult_differential_expression.tab") 
PipelineNormalised <-  read_csv("NormalisedReads.csv") %>% dplyr::select(-ensemblID)

#Most variant rows
MostVarRows5 <- PipelineNormalised %>% mutate(RowVars = rowVars(as.matrix(PipelineNormalised[,2:ncol(PipelineNormalised)]))) %>% 
  filter(quantile(RowVars, 0.5)<RowVars) %>% dplyr::select(-RowVars) %>% column_to_rownames("external_gene_id")
svg("HeatmapVar5.svg")
pheatmap(MostVarRows5, show_rownames=F, cluster_rows=FALSE, scale="row", fontsize = 15)
dev.off()

MostVarRows9 <- PipelineNormalised %>% mutate(RowVars = rowVars(as.matrix(PipelineNormalised[,2:ncol(PipelineNormalised)]))) %>% 
  filter(quantile(RowVars, 0.9)<RowVars) %>% dplyr::select(-RowVars) %>% column_to_rownames("external_gene_id")
svg("HeatmapVar9.svg")
pheatmap(MostVarRows9, show_rownames=F, cluster_rows=FALSE, scale="row", fontsize = 15)
dev.off()

MostVarRows99 <- PipelineNormalised %>% mutate(RowVars = rowVars(as.matrix(PipelineNormalised[,2:ncol(PipelineNormalised)]))) %>% 
  filter(quantile(RowVars, 0.99)<RowVars) %>% dplyr::select(-RowVars) %>% column_to_rownames("external_gene_id")
svg("HeatmapVar99.svg")
pheatmap(MostVarRows99, show_rownames=F, cluster_rows=FALSE, scale="row", fontsize = 15)
dev.off()

MostVarRows999 <- PipelineNormalised %>% mutate(RowVars = rowVars(as.matrix(PipelineNormalised[,2:ncol(PipelineNormalised)]))) %>% 
  filter(quantile(RowVars, 0.999)<RowVars) %>% dplyr::select(-RowVars) %>% column_to_rownames("external_gene_id")
svg("HeatmapVar999.svg")
pheatmap(MostVarRows999, cluster_rows=FALSE, scale="row", fontsize = 15)
dev.off()