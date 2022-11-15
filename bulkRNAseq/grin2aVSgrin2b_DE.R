library(DESeq2)
library(tidyr)
library(dplyr)
library(tximport)

##Compare HT in two different groups--different regions, different genes (Grin2a vs Grin2b), different ages, etc
##dat is count data, meta is meat data, cond is the name of the column with condition, sampType is the name of the group (needs to have 2 levels)
DECompare<-function(dat,meta,cond="condition",sampType="Gene")
{
meta=meta[,c(cond,sampType)] 
colnames(meta)=c("cond","sampType")
dat=dat[,meta$cond %in% c("WT","HT")]
meta=meta[meta$cond %in% c("WT","HT"),]
meta<- meta %>% unite(lab,cond,sampType,sep="_",remove=F)
print(head(meta))
dds <- DESeqDataSetFromMatrix(countData = dat,colData = meta,design= ~0+lab)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
print(resultsNames(dds))
res <- results(dds, contrast=c(1,-1,-1,1))
res=data.frame(res)
res=res[order(res$pvalue),]
res=res[!is.na(res$padj),]
return(res)
}

# load in libraries, path to data directory, and list of sample names.
    # i.e.
        # filePath <- "/my/data/dir/"
        # sampleNames <- c("grin2a-HT1", "grin2a-HT2", "grin2a-WT1", "grin2a-WT2", "grin2a-WT3", ...)
filePath <- c(rep(c("/path/to/grin2a/"),10),
              rep(c("/path/to/grin2b/"),10)
             )
sampleNames <- c(      )
files <- file.path(filePath, paste0(sampleNames, "/", sampleNames, ".genes.results"))
dat_12wPFC <- tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)

dat_12wPFC_counts <- dat_12wPFC$counts          
meta_12wPFC = data.frame(sample=sampleNames,
                         condition=c(rep(c("HT"),5), rep(c("WT"),5), rep(c("HT"),5), rep(c("WT"),5)),
                         Gene=c(rep(c("Grin2a"),10), rep(c("Grin2b"),10))
                         )
runDE <- DECompare(round(dat_12wPFC_counts), meta_12wPFC)
head(runDE)
write.csv(runDE,"grin2aVSgrin2b_DE.csv")