library(tidyverse)
library(ggplot2) 

load("D:/GoogleDrive/Work/UCL/F210IAnalysis/sgseq/F210I_adult_sc_res_clean_novel.RData")

f210i.dir <- "D:/GoogleDrive/Work/UCL/F210IAnalysis/sgseq/"  

makePieChart <- function(sgseqRes, title, FDRlimit, outFolder){
  # filter by FDR < 0.01!
  res.sig <- dplyr::filter(sgseqRes, FDR < FDRlimit) %>% select(one_of(c("groupID", "variantType", "FDR")) ) 
  # for each groupID take one event
  res.sig.by.group <- res.sig %>% group_by(groupID) %>% do(head(.,1)) 
  res.sig <- dplyr::filter(sgseqRes, FDR < FDRlimit) %>% select(one_of(c("groupID", "variantType", "FDR")) ) 
  res.split <- unlist(strsplit(res.sig.by.group$variantType, "+", fixed= TRUE) )
  res.split <- replace(res.split, res.split=="SE:I" | res.split=="SE:S" | res.split=="S2E:I" |res.split=="S2E:S", "Cassette exons")
  res.split <- replace(res.split, res.split=="RI:E" | res.split=="RI:R", "Retained introns")
  res.split <- replace(res.split, res.split=="AS" | res.split=="AFE", "Alternative first exon")
  res.split <- replace(res.split, res.split=="AE" | res.split=="ALE", "Alternative last exon")
  res.split <- replace(res.split, res.split=="A3SS:P" | res.split=="A3SS:D", "Alternative 3' site")
  res.split <- replace(res.split, res.split=="A5SS:P" | res.split=="A5SS:D", "Alternative 5' site")
  res.split <- replace(res.split, res.split=="MXE" , "Mutually exclusive exons")
  
  res.events <- table(res.split) 
  print(res.events)
  print( sum(res.events) )
  
  num.events <- nrow(res.sig.by.group)
  
  
  res.events.plot <- as.data.frame(res.events) %>% column_to_rownames("res.split") 
  names(res.events.plot) <- "varcounts" 
  res.events.plot$variantType <- row.names(res.events.plot)
  res.events.plot <- res.events.plot[ order(res.events.plot$varcounts,decreasing=FALSE),]  
  res.events.plot$variantType <- factor(res.events.plot$variantType, levels = rev(res.events.plot$variantType) )
  res.events.plot <- dplyr::mutate(res.events.plot, pos = cumsum(varcounts) - 0.5*varcounts)
  res.events.plot <-  drop_na(res.events.plot)
  res.events.plot$prop <- signif( (res.events.plot$varcounts / sum(res.events.plot$varcounts) ) * 100, 3) 
  res.events.plot$prop <- paste0(res.events.plot$prop, "%")
  
  pie <- ggplot(res.events.plot, aes(x="", y = varcounts, fill = variantType)) + 
    geom_bar(width = 1, stat="identity")  + 
    geom_text(aes(x=1.6, y = pos, label = prop), size = 3 ) +
    coord_polar(theta="y") + 
    scale_fill_brewer("",palette="Dark2", direction = 1) + 
    theme_void() +
    ggtitle(title)+
    annotate("text", x = 1.8, y = sum(res.events.plot$varcounts) / 2, 
             label = paste0("total number of events at adjusted p < ",FDRlimit,": ", num.events ) )
  
  print(pie)
  ggsave(paste0(outFolder,"/", title, "_sgseq_pie_chart.png") )
  
  write.table( res.events.plot, paste0(outFolder, "/", title, "_sgseq_variant_type_table.tab"), col.names = TRUE, sep = "\t")
} 

makePieChart(res.clean, "F210IHET01", 0.01, f210i.dir)
makePieChart(res.clean, "F210IHET05", 0.05, f210i.dir)
makePieChart(res.clean, "F210IHET1", 0.1, f210i.dir)
makePieChart(res.clean, "F210IHET2", 0.2, f210i.dir)



