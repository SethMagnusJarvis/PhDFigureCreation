#C9Volcano
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

#Read in files
PR6MDE <- read_csv("DifferentialExpressionPR6MNoType.csv")
GR6MDE <- read_csv("DifferentialExpressionGR6MNoType.csv")
PR12MDE <- read_csv("DifferentialExpressionPR12MNoType.csv")
GR12MDE <- read_csv("DifferentialExpressionGR12MNoType.csv")

#Create volcano plots 
PR6M <- EnhancedVolcano(PR6MDE,
                        lab = PR6MDE$external_gene_name,
                        x = 'FC',
                        y = 'PVAL',
                        title = 'PR 6 Months Volcano')

GR6M <- EnhancedVolcano(GR6MDE,
                        lab = GR6MDE$external_gene_name,
                        x = 'FC',
                        y = 'PVAL',
                        title = 'GR 6 Months Volcano Volcano')

PR12M <- EnhancedVolcano(PR12MDE,
                       lab = PR12MDE$external_gene_name,
                       x = 'FC',
                       y = 'PVAL',
                       title = 'PR 12 Months Volcano')

GR12M <- EnhancedVolcano(GR12MDE,
                       lab = GR12MDE$external_gene_name,
                       x = 'FC',
                       y = 'PVAL',
                       title = 'GR 12 Months Volcano')

#Create a PR6M volcano plot without the gene egf as it distorts the graph
PR6MDEFilt <- PR6MDE %>% filter(external_gene_name != "Egf")
PR6MFilt <- EnhancedVolcano(PR6MDEFilt,
                        lab = PR6MDEFilt$external_gene_name,
                        x = 'FC',
                        y = 'PVAL',
                        title = 'PR 6 Months Volcano No Egf')
