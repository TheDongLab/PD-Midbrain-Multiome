# =============================================================================
# Step 01 - Build and merge per-sample multiome Seurat objects
# -----------------------------------------------------------------------------
# Purpose : Read Cell Ranger ARC outputs for each sample, build a Seurat object
#           with paired RNA + ATAC (Signac ChromatinAssay) assays, compute
#           per-cell QC metrics, and merge all samples into one object.
# Inputs  : Cell Ranger ARC per-sample outputs under CELLRANGER_DIR
#             - filtered_feature_bc_matrix.h5
#             - atac_fragments.tsv.gz
#             - per_barcode_metrics.csv
# Outputs : <DIR_FINALWNN>/seurat_objects.rda  (per-sample list + fragments)
#           <DIR_FINALWNN>/merged_seurat.rda   (merged object)
# Run     : Rscript Combining-Samples.R
# Notes   : Package installation commands from the original script were removed;
#           install dependencies separately (see README).
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(ggplot2)
  library(hdf5r)
})

# Helper: extract the sample ID (3rd-to-last element of the file path).
last3rd <- function(invect) invect[length(invect) - 2]

# Exclude aggregated batches and the old versions of resequenced samples
# (BN1610SNR_old, BN1862SN_old, BN1959SN_old; *_aggr; *_combine).
exclude_pattern <- "_aggr|_combine|_old"

## ---- Locate per-sample input files ------------------------------------------

fragFiles <- list.files(CELLRANGER_DIR, pattern = "atac_fragments.tsv.gz$",
                        full.names = TRUE, recursive = TRUE)
fragFiles <- fragFiles[grep(exclude_pattern, fragFiles, invert = TRUE)]
names(fragFiles) <- sapply(strsplit(fragFiles, "/"), last3rd)

h5Files <- list.files(CELLRANGER_DIR, pattern = "filtered_feature_bc_matrix.h5$",
                      full.names = TRUE, recursive = TRUE)
h5Files <- h5Files[grep(exclude_pattern, h5Files, invert = TRUE)]
names(h5Files) <- sapply(strsplit(h5Files, "/"), last3rd)

metaFiles <- list.files(CELLRANGER_DIR, pattern = "per_barcode_metrics.csv$",
                        full.names = TRUE, recursive = TRUE)
metaFiles <- metaFiles[grep(exclude_pattern, metaFiles, invert = TRUE)]
names(metaFiles) <- sapply(strsplit(metaFiles, "/"), last3rd)

## ---- Build one Seurat object per sample -------------------------------------

seurat_objects <- list()
listnames <- c()
fragments <- list()

for (sampleID in names(h5Files)[1:101]) {
  print(sampleID)

  # The 10x hdf5 file contains both modalities.
  inputdata.10x <- Read10X_h5(h5Files[sampleID])
  metadata <- read.csv(metaFiles[sampleID], header = TRUE, row.names = 1,
                       stringsAsFactors = FALSE)

  # ---- RNA assay ----
  rna_counts <- inputdata.10x$`Gene Expression`
  pbmc <- CreateSeuratObject(counts = rna_counts, assay = "RNA",
                             project = sampleID, meta.data = metadata)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

  # ---- ATAC assay (standard chromosomes only) ----
  atac_counts <- inputdata.10x$Peaks
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]

  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- GENOME_BUILD

  frag.file <- fragFiles[sampleID]
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = GENOME_BUILD,
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
  )
  pbmc[["ATAC"]] <- chrom_assay

  # ---- Per-cell QC (RNA + ATAC) ----
  p <- VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"),
               ncol = 3, log = TRUE, pt.size = 0) + NoLegend()
  ggsave(paste0("./01-Clusterings3/", sampleID, "_basic_QC.VlnPlot.pdf"), p,
         width = 10, height = 6)

  DefaultAssay(pbmc) <- "ATAC"
  pbmc <- NucleosomeSignal(pbmc)
  pbmc <- TSSEnrichment(pbmc, fast = FALSE)
  pbmc$pct_reads_in_peaks <- pbmc$atac_peak_region_fragments / pbmc$atac_fragments * 100
  pbmc$blacklist_fraction <- FractionCountsInRegion(pbmc, assay = "ATAC",
                                                    regions = blacklist_hg38_unified)

  # Tag cell barcodes with the sample ID.
  obj <- pbmc
  obj$sample_id <- sampleID
  colnames(obj) <- paste(sampleID, colnames(obj), sep = "#")
  print(obj)

  seurat_objects <- c(seurat_objects, obj)
  listnames <- c(listnames, sampleID)
  fragments <- c(fragments, obj@assays$ATAC@fragments)
}

save(seurat_objects, listnames, fragments,
     file = file.path(DIR_FINALWNN, "seurat_objects.rda"))

## ---- Merge all samples ------------------------------------------------------

load(file.path(DIR_FINALWNN, "seurat_objects.rda"))

merged_seurat <- merge(
  x = seurat_objects[[2]],
  y = seurat_objects[3:length(seurat_objects)],
  project = "merged_project"
)
save(merged_seurat, file = file.path(DIR_FINALWNN, "merged_seurat.rda"))
