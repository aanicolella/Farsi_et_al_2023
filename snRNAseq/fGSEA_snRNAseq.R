# fGSEA code written by Sameer Aryal 

doGoAnalysis <- function(tbl, pathway, tStat){
  library(lazyeval)
  tStat <- as.name(tStat)
  
  bm<-read.csv("ensemblId_mouseGeneName_humanGeneName.csv")
  res <- inner_join(tbl, bm, by=c("gene"="external_gene_name"))
 

  psd.res <- res %>%
    dplyr::select(hsapiens_homolog_associated_gene_name, tStat) %>%
    na.omit() %>%
    distinct() %>%
    rename(tStatistic=tStat) %>%
    group_by(hsapiens_homolog_associated_gene_name) %>%
    slice_max(abs(tStatistic))
 
  ranks <- deframe(psd.res)
  pathways.go <- gmtPathways(pathway)
  fgseaRes.psd <- fgsea(pathways=pathways.go, stats=ranks)
 
  fgseaRes.psd <- fgseaRes.psd %>%
    as_tibble() %>%
    arrange(padj, desc(NES))
 
  return (fgseaRes.psd)
 
}

require(tidyverse)
require(fgsea)
res <- read.csv("/path/to/DE_results.csv")
colnames(res)[1] <- "gene"
pthway <- "/path/to/c5.all.v7.2.symbols.gmt"
statCol <- "insert_logFC_col_name_here"
gseaRes1 <- doGoAnalysis(res, pthway, statCol)
gseaRes3 <- gseaRes1 %>% mutate(leadingEdge = map_chr(leadingEdge, toString))
write_csv(gseaRes3, "/path/to/GSEA_results.csv")