# Variables to replace
    # EXP: name of experiment, usually the gene of interest
    # REGION: brain region you are studying
    # X, Y, Z: number of replicates for HT, KO, WT

# load in libraries, path to data directory, and list of sample names.
    # i.e.
        # filePath <- "/my/data/dir/"
        # sampleNames <- c("HT1", "HT2", "HT3", "KO1", "KO2", "KO3", "WT1", "WT2", "WT3")
        
library("DESeq2")
library("tximport")
filePath <- "/path/to/data/REGION/"
sampleNames <- c()

# Set up metadata
    # sample = individual sample name
colData = data.frame(genotype=c(rep(c("EXP_HT"), X), rep(c("EXP_KO"), Y), rep(c("EXP_WT"), Z)),sample=sampleNames)
colData$genotype = relevel(colData$genotype, "EXP_WT")
files <- file.path(filePath, paste0(sampleNames, "/", sampleNames, ".genes.results"))
names(files) <- sampleNames
rsem.in <- tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)

# Optional: find genes below expression threshhold. Create 'NA' data frame to re-add lowly expressed 
# genes for easier comparison between analyses.
drop <- rsem.in
drop$abundance <-
    drop$abundance[apply(drop$length,
                             1,
                             function(row) !all(row > 5)),]                              
drop$counts <-
  drop$counts[apply(drop$length,
                            1,
                             function(row) !all(row > 5)),]
drop$length <-
  drop$length[apply(drop$length,
                             1,
                             function(row) !all(row > 5)),]
                    
dropped <- drop$counts
dropRes <- data.frame(row.names=make.unique(row.names(drop$counts)), matrix(ncol=6, nrow=dim(drop$counts)[1]))
cols = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
colnames(dropRes) <- cols
EXPName1 <- data.frame(do.call('rbind', strsplit(as.character(row.names(dropRes)), '_', fixed=TRUE)))
EXPName1 <- EXPName1["X2"]               
dropRes$symbol <- EXPName1$X2

# Filter lowly expressed genes found above out of the main dataset for the DE analysis
rsem.in$abundance <-
  rsem.in$abundance[apply(rsem.in$length,
                             1,
                             function(row) all(row > 5 )),]
rsem.in$counts <-
  rsem.in$counts[apply(rsem.in$length,
                             1,
                             function(row) all(row > 5 )),]
rsem.in$length <-
  rsem.in$length[apply(rsem.in$length,
                             1,
                             function(row) all(row > 5 )),]

# Convert data structure to proper format and run DESeq
dds <- DESeqDataSetFromTximport(rsem.in, colData, ~genotype)
dds <- DESeq(dds)
check <- estimateSizeFactors(dds)
norm <- counts(check, normalized=TRUE)
                       
# Combine normalized count dataframe with the 'NA' dataframe and save
all <- rbind(norm, dropped)                                          
EXPName <- data.frame(do.call('rbind', strsplit(as.character(row.names(all)), '_', fixed=TRUE)))
EXPName <- EXPName["X2"]               
row.names(all) <- EXPName$X2
all <- all[order(row.names(all)),]
head(all)
write.csv(all, "EXP_REGION_normalizedCounts.csv", row.names=TRUE)
                       
# Plot pca using different metadata conditions to see relationships between samples
pdf(file="EXP_REGION_pca.pdf")
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="sample")
dev.off()
pdf(file="EXP_REGION_genotypepca.pdf")
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="genotype")
dev.off()

# Extract DE results from different genotype comparisions
    # Currently set up such that the LFC directionality is in relation to the changes in the first 
    # genotype listed: i.e. + lfc == upregulation in that geno, - lfc == downregulation
resHT <- results(dds, contrast=c("genotype", "EXP_HT","EXP_WT"))
resKO <- results(dds, contrast=c("genotype", "EXP_KO","EXP_WT"))
head(results(dds, tidy=TRUE), n=6)
summary(dds)

# Scale HTvsWT results using lfcShrink(): https://rdrr.io/bioc/DESeq2/man/lfcShrink.html
   # Combine post-shrinkage results with 'NA' dropped dataframe then write to .csv
resHT <- lfcShrink(dds=dds, contrast=c("genotype", "EXP_HT", "EXP_WT"), res=resHT, type="normal")
EXPNameHT <- data.frame(do.call('rbind', strsplit(as.character(row.names(resHT)), '_', fixed=TRUE)))
EXPNameHT <- EXPNameHT["X2"]               
resHT$symbol <- EXPNameHT$X2
resMergeHT <- rbind(resHT, dropRes)
dim(resMergeHT)
resMergeHT <- resMergeHT[order(resMergeHT$pvalue),]
write.csv(resMergeHT, "EXP_REGION_HTWT-DESeq2_allgenes.csv")
                       
# Repeat above for KOvsWT
resKO <- lfcShrink(dds=dds, contrast=c("genotype", "EXP_KO", "EXP_WT"), res=resKO, type="normal")
EXPNameKO <- data.frame(do.call('rbind', strsplit(as.character(row.names(resKO)), '_', fixed=TRUE)))
EXPNameKO <- EXPNameKO["X2"]               
resKO$symbol <- EXPNameKO$X2
resMergeKO <- rbind(resKO, dropRes)
dim(resMergeKO)
resMergeKO <- resMergeKO[order(resMergeKO$pvalue),]
write.csv(resMergeKO, "EXP_REGION_KOWT-DESeq2_allgenes.csv")