library(Seurat)
library(Signac)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

out.path <- args[3]

#load scRNA data
rna <- readRDS(args[2])
#Load Annotated genome
annotations <- readRDS("hg38.rds")

#Generating paths for atac data
peaks <- paste(args[1], "filtered_peak_bc_matrix.h5", sep="")
frags <- paste(args[1], "fragments.tsv.gz", sep="")
meta <- paste(args[1], "singlecell.csv", sep="")

#initiating atac seurat object
atac.assay <- CreateChromatinAssay(
  counts = peaks,
  sep = c(":", "-"),
  fragments = atac.frags,
  annotation = annotations,
  min.cells = 1)

metadata <- read.csv(
  file = atac.meta,
  header = TRUE,
  row.names = 1,)

atac <- CreateSeuratObject(
  counts = atac.assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata,
  annotation = annotations)

atac <- TSSEnrichment(atac)
atac <- NucleosomeSignal(atac)

atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

#QC filter
atac <- subset(
  x = atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)

#Dimensionallity Reduction of ATAC data
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims=2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# quantify gene activity
gene.activities <- GeneActivity(atac, features = VariableFeatures(rna))
# add gene activities as a new assay
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna,
                                        query = atac,
                                        features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca")
# Transfer Labels
celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = rna$celltypeLabels, 
                                     weight.reduction = atac[["lsi"]], dims = 2:30)

# Add predictions to atac object
atac <- AddMetaData(atac, metadata = celltype.predictions)
#Change cell identities to predicted IDs
Idents(atac) <- "predicted.id"
atac.markers <- FindAllMarkers(object= atac)