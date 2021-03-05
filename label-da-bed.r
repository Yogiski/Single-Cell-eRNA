library(Seurat)
library(Signac)
library(GenomicRanges)
#library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

patient.id <- args[4]
out.path <- args[3]

#load scRNA data
rna <- readRDS(args[2])

#Set Identities to cell type labels and save cell ids in separate csv's for bamslice script
Idents(rna) <- "celltypeLabels"

for (i in levels(rna)) {
    ids <- names(Idents(subset(rna, idents=i)))
    name <- gsub(" ", "-", i)
    path <- paste(out, "/", patient.id, "/cluster-rna-ids/", name, "-ids.csv", sep="")
    f <- file(path)
    writeLines(ids, f)
}

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
#Calc DA for each cluster using wilcox
atac.markers <- FindAllMarkers(object = atac, min.pct = 0.25)

markers.outpath <- paste(out,"/", patient.id, "/unfiltered-da.csv", sep="")
atac.outpath <- paste(out,"/", patient.id, "/atac-labeled.rds", sep="")
write.csv(atac.markers, markers.outpath)
saveRDS(atac.markers, atac.outpath)

bed.to.granges <- function(bed.path){
    
    df <- read.table(bed.path)
    res <- GRanges(seqnames = df$V1,
                   ranges = IRanges(start = df$V2, end = df$V3),
                   names = df$V4,
                   names2 = df$V5,
                   elements = df$V6)
    return(res)
}

distal.path <- "GRCh38-ccREs.dELS.bed"
prox.path <- "GRCh38-ccREs.pELS.bed"

distal <- bed.to.granges(distal.path)
prox <- bed.to.granges(prox.path)

grpd.markers <- atac.markers %>% group_by(cluster) %>% filter(p_val_adj <= 0.05)

for (g in unique(grpd.markers$cluster)) {
    #Make granges object by celltype
    markers.g <- grpd.markers %>% filter(cluster == g)
    g.da.granges <- granges(atac[markers.g$gene,])
    #Find overlaps with annotated enhancers
    g.da.distal <- subsetByOverlaps(g.da.granges, distal, minoverlap=300)
    g.da.prox <- subsetByOverlaps(g.da.granges, proximal, minoverlap=300)
    #Remove duplicates and join proximal and distal into one object
    g.da.prox <- setdiff(g.da.prox, subsetByOverlaps(g.da.distal, g.da.prox))
    g.da.enh <- union(g.da.distal, g.da.prox)
    #Generate filename and save GRanges object as bed file
    name <- gsub(" ", "-", g)
    bed.outpath <- paste(out, "/", patient.id, "/bed-files/", name, "-da-enhancers.bed", sep="")
    export.bed(g.da.enh, bed.outpath)
    }
    
