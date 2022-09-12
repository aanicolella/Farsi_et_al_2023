library(dplyr)
library(Seurat)
library(edgeR)
library(fdrtool)

# Load post cell type labeling data object
object <- readRDS("/path/to/object.rds")
Idents(object) <- object@meta.data$cellType
# Identify present cell types
cellTypes <- as.vector(as.data.frame(table(Idents(object)))$Var1)

# Generate pseudobulk counts for each cell type. Adjust replicates as needed.
for(i in cellTypes){
    # Extract cell type counts by replicate
    i_counts <- subset(object, idents=i)

    KO1 <- subset(i_counts, subset= orig.ident == "project_KO1")
    i_counts_KO1 <- data.frame(KO1@assays[["RNA"]]@counts)
    i_total_KO1 <- data.frame(i_KO1=rowSums(i_counts_KO1))
    head(i_total_KO1)

    KO2 <- subset(i_counts, subset= orig.ident == "project_KO2")
    i_counts_KO2 <- data.frame(KO2@assays[["RNA"]]@counts)
    i_total_KO2 <- data.frame(i_KO2=rowSums(i_counts_KO2))
    head(i_total_KO2)

    KO3 <- subset(i_counts, subset= orig.ident == "project_KO3")
    i_counts_KO3 <- data.frame(KO3@assays[["RNA"]]@counts)
    i_total_KO3 <- data.frame(i_KO3=rowSums(i_counts_KO3))
    head(i_total_KO3)
    
    HT1 <- subset(i_counts, subset= orig.ident == "project_HT1")
    i_counts_HT1 <- data.frame(HT1@assays[["RNA"]]@counts)
    i_total_HT1 <- data.frame(i_HT1=rowSums(i_counts_HT1))
    head(i_total_HT1)

    HT2 <- subset(i_counts, subset= orig.ident == "project_HT2")
    i_counts_HT2 <- data.frame(HT2@assays[["RNA"]]@counts)
    i_total_HT2 <- data.frame(i_HT2=rowSums(i_counts_HT2))
    head(i_total_HT2)

    HT3 <- subset(i_counts, subset= orig.ident == "project_HT3")
    i_counts_HT3 <- data.frame(HT3@assays[["RNA"]]@counts)
    i_total_HT3 <- data.frame(i_HT3=rowSums(i_counts_HT3))
    head(i_total_HT3)
    
    WT1 <- subset(i_counts, subset= orig.ident == "project_WT1")
    i_counts_WT1 <- data.frame(WT1@assays[["RNA"]]@counts)
    i_total_WT1 <- data.frame(i_WT1=rowSums(i_counts_WT1))
    head(i_total_WT1)

    WT2 <- subset(i_counts, subset= orig.ident == "project_WT2")
    i_counts_WT2 <- data.frame(WT2@assays[["RNA"]]@counts)
    i_total_WT2 <- data.frame(i_WT2=rowSums(i_counts_WT2))
    head(i_total_WT2)

    WT3 <- subset(i_counts, subset= orig.ident == "project_WT3")
    i_counts_WT3 <- data.frame(WT3@assays[["RNA"]]@counts)
    i_total_WT3 <- data.frame(i_WT3=rowSums(i_counts_WT3))
    head(i_total_WT3)
    
    # Create counts table
        # Save to pseudobulk counts folder ("pseudobulkSheets/" in this example)
        # Replace date to reflect that of your data
    i_totalAll <- cbind(i_total_KO1,i_total_KO2,i_total_KO3,
                        i_total_HT1,i_total_HT2,i_total_HT3,
                        i_total_WT1,i_total_WT2,i_total_WT3)
    colnames(i_totalAll) <- c(paste0(i,"_KO1"),paste0(i,"_KO2"),paste0(i,"_KO3"),
                              paste0(i,"_HT1"),paste0(i,"_HT2"),paste0(i,"_HT3"),
                              paste0(i,"_WT1"),paste0(i,"_WT2"),paste0(i,"_WT3"))
    write.csv(i_totalAll, paste0("pseudobulkSheets/",i,"_totalExpression_date.csv"))
}

##### DE within each cell type ######
for (i in celltypes){
  #load data from pseudobulk sheets
  fp <- paste0("pseudobulkSheets/", i, "_totalExpression_DATE.csv")
  rsem.in <- as.matrix(read.csv(fp, row.names = 1))
  
  #make sure this matches your data
  colData = data.frame(genotype=c(rep(c("KO"), 3), rep(c("HT"), 3), rep(c("WT"), 3)),sample=colnames(rsem.in))
  colData$genotype = as.factor(colData$genotype)
  colData$genotype = relevel(colData$genotype, "WT")
  # Add batch number if data comes from multiple experiments, example below
  #colData$batch = as.factor(c(1,1,1,1,2,1,1,1,1,2,1,1,2,2,2))
  
  dge = DGEList(counts = rsem.in, group = colData$genotype)
  keep = filterByExpr(dge)
  dge = dge[keep, , keep.lib.sizes = FALSE]
  dge = calcNormFactors(dge)
  
  ### GLM (LRT)
  # design matrix for one covariate
  design = model.matrix(~genotype, data = colData)
  # Design matrix for multiple covariates (i.e. genotype and multiple batches)
  #design = model.matrix(~genotype + batch, data = colData)
  dge = estimateDisp(dge, design)
  
  ##LRT
  fit = glmFit(dge, design)
  lrtHT = glmLRT(fit, coef = 2)
  lrtKO = glmLRT(fit, coef = 3)
  
  rem <- as.data.frame(keep[keep==FALSE])
  dropRes <- data.frame(row.names=row.names(rem), matrix(ncol=5, nrow=length(row.names(rem))))
  cols = c("logFC","logCPM","LR","PValue","FDR")
  colnames(dropRes) = cols
  
  resHT = lrtHT$table
  resHT$FDR = p.adjust(resHT$PValue, "fdr")
  # Optional recalculation of nominal and adjusted p-values using fdrtool
  #norm_HT=-sign(resHT$logFC)*qnorm(resHT$PValue/2)
  #pval_HT=fdrtool(norm_HT,plot=F)$pval
  #resHT["pval_fdrtools"]=pval_HT
  #resHT["padj_fdrtools"]=p.adjust(pval_HT,"fdr")
  resHT_merge = rbind(resHT,dropRes)
  dir.create("HT/byCellType",showWarnings=FALSE)
  dir.create(paste0("HT/byCellType/", i),showWarnings=FALSE)
  write.csv(resHT_merge, paste0("HT/byCellType/",i,"/",i,"_HTWT_EdgeR_DATE.csv"))
  
  resKO = lrtKO$table
  resKO$FDR = p.adjust(resKO$PValue, "fdr")
  # Optional recalculation of nominal and adjusted p-values using fdrtool
  #norm_KO=-sign(resKO$logFC)*qnorm(resKO$PValue/2)
  #pval_KO=fdrtool(norm_KO,plot=F)$pval
  #resKO["pval_fdrtools"]=pval_KO
  #resKO["padj_fdrtools"]=p.adjust(pval_KO,"fdr")
  resKO_merge = rbind(resKO,dropRes)
  dir.create("KO/byCellType",showWarnings=FALSE)
  dir.create(paste0("KO/byCellType/", i),showWarnings=FALSE)
  write.csv(resKO_merge, paste0("KO/byCellType/",i,"/",i,"_KOWT_EdgeR_DATE.csv"))
  
}