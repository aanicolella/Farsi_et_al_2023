# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(RColorBrewer)

# Downsampling: start from CellRanger aggr output

object <- "projectName"
object.data <- Read10X(data.dir="parentDir/data/aggr/projectName/outs/filtered_feature_bc_matrix/")
object <- CreateSeuratObject(counts=object.data, project = "projectName", min.cells = 3, min.features =200)

# Extract barcodes from cell names for sample mapping

barcodes <- row.names(object@meta.data)
nums <- read.table(text=barcodes, sep="-")
head(nums)

# Map barcode numbers to sample names & genotypes to samples

sampleName <- plyr::mapvalues(
    x = nums$V2, 
    from = c(1,2,3), 
    to = c("samp1","samp2","samp3")
)
object@meta.data$orig.ident <- sampleName

# Add genotype to the metadata
    # put your sample names (called "orig.ident") in the from list and the genotypes of 
    # each sample in the to list
genotype <- plyr::mapvalues(
    x = object@meta.data$orig.ident, 
    from = c("samp1", "samp2", "samp3"), 
    to = c("HT", "KO", "WT")
)

object@meta.data$genotype <- genotype

# Checkpoint 1: save merged raw data

saveRDS(object, "project_aggrObject_DATE.rds")
#object <- readRDS("project_aggrObject_DATE.rds")

# Reassign sample names as replicate name
as.data.frame(table(object@meta.data$orig.ident))
as.data.frame(table(object@meta.data$genotype))
Idents(object) <- object@meta.data$orig.ident


# Visualize QC metrics as violin plot
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol=2, pt.size=0)

# Filter out cells with too few features 
    # Further filtering may be necessary depending on what your data looks like
object <- subset(object, subset = nFeature_RNA > 500)
# Normalizing data
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
# Post filtering, normalization violin plot of QC metrics
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA"), ncol=2, pt.size=0)

# Feature selection (via identification of highly variable features)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

# Scaling data
    # linear transformation (scaling): preprocessing step prior to dimension reduction
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)

# Linear Dimensional Reduction (using PCA)
    # Don't know what PCA is? Read this: https://towardsdatascience.com/a-one-stop-shop-for-principal-component-analysis-5582fb7e0a9c
object <- RunPCA(object, features = VariableFeatures(object=object))
# Look at just top 5 genes for each principle component (PC)
# Look at the PCA visualization colored by genotype and replicate
DimPlot(object, reduction="pca", group.by="genotype")
DimPlot(object, reduction="pca")
ElbowPlot(object)

# Save the PCA plots
    # make sure you have made a "plots" directory
    # replace "project" and "DATE" to reflect your project and the current date
dir.create("plots", showWarnings=F)
pdf("plots/project_pcaBySample_DATE.pdf")
DimPlot(object, reduction="pca")
dev.off()
pdf("plots/project_pcaByGenotype_DATE.pdf")
DimPlot(object, reduction="pca", group.by="genotype")
dev.off()

# Cluster Cells
    # This uses a graph-based clustering approach, built upon the strategies in this paper:
    # https://www.cell.com/fulltext/S0092-8674(15)00549-8
    # Reference the Seurat "Guided Clustering Tutorial" for further explanation of this method

    # Replace "X" with the PC number that is the "elbow" in the ElbowPlot two cells above
        # consult the Seurat tutorial if you do not know what this means
    # Replace "Y" with your desired resolution (consult the Seurat tutorial if you do not know what this means)
object <- FindNeighbors(object, dims=1:X)
object <- FindClusters(object, resolution = Y)

# Run non-linear dimension reduction (UMAP/tSNE)
    # Don't know what UMAP or tSNE are? Read this: https://towardsdatascience.com/how-exactly-umap-works-13e3040e1668
object <- RunUMAP(object, dims=1:X)

# Look at UMAP clusters before running doublet identification
DimPlot(object, reduction="umap", label=TRUE)
DimPlot(object, reduction="umap", group.by="genotype")
DimPlot(object, reduction="umap", group.by="orig.ident")

# Write counts to matrix for double removal

writeMM(object@assays$RNA@counts,"counts.txt")
# Stop here to run RunScrublet.py
    # Use the below command to run this script. Please note, you must install Scrublet 
    # (preferably in a conda environment) before running the script.
       # $ python RunScrublet.py counts.txt output.txt

# Import scrublet results, add to metadata and look at score distribution
object@meta.data["scrublet"] = scan("output.txt")
hist(object@meta.data$scrublet, col="gray")

# Assign doublet label to cells above scrublet score threshhold
    # Adjust THRESHOLD depending on the score distribution 
object@meta.data["doublet"]="singlet"
object@meta.data[object@meta.data[,"scrublet"]>THRESHOLD,"doublet"]="doublet"

# Assess number of doublets and their distribution
    # make sure you create a "dataCleaning" directory to hold the plots below

table(object@meta.data$doublet)
table(object@meta.data$orig.ident)
FeaturePlot(object, features="scrublet")
dir.create("plots/dataCleaning/", showWarnings=F)
pdf("plots/dataCleaning/doublets.pdf")
FeaturePlot(object, features="scrublet")
dev.off()
pdf("plots/dataCleaning/labeledDoublets.pdf")
FeaturePlot(object, features="scrublet", label=TRUE)
dev.off()

# Label small clusters with high % doublets and remove them
    # open the plots created above and see which clusters have a high frequency of 
    # high doublet scores. Put the cluster numbers in the list below 
    # (must be formatted as a list of strings i.e. c("17","18","20"))
doublet_clusters <- c()
object@meta.data[object@active.ident %in% doublet_clusters,"doublet"]="doublet"
# Check that you have the doublet clusters properly identified before removal
DimPlot(object, reduction="umap", group.by="doublet")
pdf("plots/dataCleaning/doubletClassFinal.pdf")
DimPlot(object, reduction="umap", group.by="doublet")
dev.off()

# Remove clusters labeled as doublet clusters above
object=subset(object,doublet=="singlet")

# Check removed clusters to make sure you removed all needed clusters
pdf("plots/dataCleaning/removedDoublets.pdf")
FeaturePlot(object, features="scrublet")
dev.off()
pdf("plots/dataCleaning/labeledremovedDoublets.pdf")
FeaturePlot(object, features="scrublet", label=TRUE)
dev.off()

# Checkpoint 2: Pre-labeled object with doublets removed
    # Replace "project" and "DATE" to reflect your project name and the date
saveRDS(object, "project_preLabel_DATE.rds")
#object <- readRDS("project_preLabel_DATE.rds")

# Evaluate cell type markers and determine cluster cell types
    # Note which cluster numbers express which cell types' marker genes
        # Clusters with high expression of the marker genes for a cell type will be 
        # assigned that cell type's label
    # This template only includes markers for the major cell types found in PFC data. You
        # may add additional tests for your own cell types of interest. 
DimPlot(object, reduction="umap", label=TRUE)
DimPlot(object, reduction="umap", label=TRUE) + NoLegend()

# neurons (Excitory and inhibitory) - 
FeaturePlot(object, features = c("Snap25"), label=TRUE)
FeaturePlot(object, features = c("Rbfox3"), label=TRUE)
# Inhibitory Neurons- 
FeaturePlot(object, features = c("Gad1"), label=TRUE)
FeaturePlot(object, features = c("Gad2"), label=TRUE)
# Excitatory neurons- 
FeaturePlot(object, features = c("Slc17a6"), label=TRUE)
FeaturePlot(object, features = c("Slc17a7"), label=TRUE)
FeaturePlot(object, features = c("Neurod6"), label=TRUE)
# Microglia - 
FeaturePlot(object, features = c("Csf1r"), label=TRUE)
# OPC - 
FeaturePlot(object, features = c("Pdgfra"), label=TRUE)
FeaturePlot(object, features = c("Olig1"), label=TRUE)
# ODC (mODC and nODC) - 
FeaturePlot(object, features = c("Plp1"), label=TRUE)
FeaturePlot(object, features = c("Mbp"), label=TRUE)
FeaturePlot(object, features = c("Mobp"), label=TRUE)
FeaturePlot(object, features = c("Olig1"),label=TRUE)
# Astrocytes -
FeaturePlot(object, features = c("Slc1a3"), label=TRUE)
FeaturePlot(object, features = c("S100b"), label=TRUE)
FeaturePlot(object, features = c("Gfap"), label=TRUE)
# Endothelial
FeaturePlot(object, features = c("Flt1"), label=TRUE)
# Other Vascular Markers (Pericytes, Smooth Muscle, etc) -
FeaturePlot(object, features = c("Bgn"), label=TRUE)
FeaturePlot(object, features = c("Vtn"), label=TRUE)
# Fibroblast-like - 
FeaturePlot(object, features = c("Bnc2"), label=TRUE)

# Assigning Cell Type Identity to clusters
    # Fill the new.cluster.ids2 list with the cell type labels for each cluster in numerical order
    # I usually create an excel sheet with a column for each cell type and put the cluster number
    # of each cluster into its cell type column. Then, I create a separate column with each cluster
    # number in a row and put the cluster number's label next to it. Then, I will simply copy the 
    # column into the new.cluster.ids2 list (make sure that you make this a list of strings by
    # adding quotes around each cell type name)
new.cluster.ids2 <- c()
names(new.cluster.ids2) <- levels(object)
object <- RenameIdents(object, new.cluster.ids2)
DimPlot(object, reduction = "umap", label = TRUE)

# Save commonly needed UMAP plots
    # replace "project" and "DATE" to reflect that of your project
pdf("plots/project_legend_cellTypes_DATE.pdf")
DimPlot(object, reduction = "umap", label = FALSE)
dev.off()
pdf("plots/project_cellTypes_DATE.pdf")
DimPlot(object, reduction = "umap", label = TRUE) + NoLegend()
dev.off()
pdf("plots/project_umapGenotype_DATE.pdf")
DimPlot(object, group.by="genotype")
dev.off()

# Checkpoint 3: Post-labeled object 
    # replace "project" and "DATE" to reflect that of your project
saveRDS(object, "project_postLabel_DATE.rds")
#object <- readRDS("project_postLabel_DATE.rds")