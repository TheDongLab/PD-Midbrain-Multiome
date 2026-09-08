# =============================================================================
# MACS peak calling, union peak set, and cell-type-specific peaks
# -----------------------------------------------------------------------------
# Purpose : Call ATAC peaks with MACS per cell type and across all cells, remove
#           blacklist / non-standard chromosomes, build an annotated union peak
#           set, quantify a peak x cell count matrix, and derive cell-type-
#           specific peaks.
# Inputs  : SEURAT_OBJECT (labeled, ATAC fragments), MACS_PATH,
#           ENCODE hg38 blacklist (downloaded).
# Outputs : <DIR_FINALWNN>/peaks+Tcells.rda,
#           <DIR_MACS>/union_peaks_with_annotations.narrowPeak,
#           <DIR_MACS>/peaks_macs_count.rda (peaks, union_peaks2, macs_count),
#           cell_type_specific_peaks_list.rds.
# Run     : Rscript Callpeaks-CellTypes.R
# Notes   : Core steps only. Downstream QC/summary figures (genomic-distribution
#           barplot, UpSet, rGREAT GO, accessibility heatmap) were removed.
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
})

load(SEURAT_OBJECT)

## ---- MACS peak calling (per cell type and bulk) -----------------------------

peaks <- CallPeaks(merged_seurat4, assay = "ATAC", macs2.path = MACS_PATH,
                   group.by = "celltype.WNN.after.combine.update2",
                   combine.peaks = FALSE,
                   outdir = DIR_MACS, fragment.tempdir = DIR_MACS, cleanup = FALSE)

peaks_merge <- CallPeaks(merged_seurat4, assay = "ATAC", macs2.path = MACS_PATH,
                         outdir = DIR_MACS, group.by = NULL,
                         fragment.tempdir = DIR_MACS, cleanup = FALSE)

# Remove peaks on non-standard chromosomes / blacklist regions.
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
peaks_merge <- keepStandardChromosomes(peaks_merge, pruning.mode = "coarse")
peaks_merge <- subsetByOverlaps(x = peaks_merge, ranges = blacklist_hg38_unified, invert = TRUE)
save(peaks, peaks_merge, file = file.path(DIR_FINALWNN, "peaks+Tcells.rda"))

## ---- Clean per-cell-type narrowPeak files -----------------------------------

standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
blacklist_gz <- file.path(DIR_MACS, "hg38-blacklist.v2.bed.gz")
download.file("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
              destfile = blacklist_gz)

clean_peaks <- function(narrowPeak_file, blacklist_file, output_file) {
  peaks <- import(narrowPeak_file, format = "narrowPeak")
  peaks_filtered <- peaks[seqnames(peaks) %in% standard_chroms]
  blacklist <- import(blacklist_file, format = "BED")
  peaks_cleaned <- subsetByOverlaps(peaks_filtered, blacklist, invert = TRUE)
  peaks_cleaned_df <- as.data.frame(peaks_cleaned)
  narrowPeak_df <- peaks_cleaned_df[, c("seqnames", "start", "end", "name", "score",
                                        "strand", "signalValue", "pValue", "qValue", "peak")]
  write.table(narrowPeak_df, output_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

peak_files <- file.path(DIR_MACS, c(
  "Astrocytes_peaks.narrowPeak", "Endothelial_cell_peaks.narrowPeak",
  "GABA_neurons_peaks.narrowPeak", "Glu_neurons_peaks.narrowPeak", "DA_neurons_peaks.narrowPeak",
  "Microglia_peaks.narrowPeak", "Monocytes_peaks.narrowPeak",
  "Oligodendrocytes_peaks.narrowPeak", "OPCs_peaks.narrowPeak",
  "Pericytes_peaks.narrowPeak", "merged_project_peaks.narrowPeak", "CD8_T_Cells_peaks.narrowPeak"))
for (file in peak_files) {
  clean_peaks(file, blacklist_gz, gsub("peaks.narrowPeak", "cleaned_peaks.narrowPeak", file))
}

## ---- Build annotated union peak set -----------------------------------------

merged_peaks <- import(file.path(DIR_MACS, "merged_project_cleaned_peaks.narrowPeak"), format = "narrowPeak")
mcols(merged_peaks)$source <- "tissueBulk"

peak_files <- list.files(DIR_MACS, pattern = "*_cleaned_peaks.narrowPeak", full.names = TRUE)
peak_files <- peak_files[!grepl("merged_project_cleaned_peaks.narrowPeak", peak_files)]

cell_type_specific_peaks <- list()
for (file in peak_files) {
  cell_type_peaks <- import(file, format = "narrowPeak")
  cell_type_name <- gsub("_cleaned_peaks.narrowPeak", "", basename(file))
  mcols(cell_type_peaks)$source <- cell_type_name
  cell_type_specific_peaks[[cell_type_name]] <- subsetByOverlaps(cell_type_peaks, merged_peaks, invert = TRUE)
}

union_peaks <- merged_peaks
for (peakss in cell_type_specific_peaks) union_peaks <- c(union_peaks, peakss)

# Per-cell-type membership columns (0/1) on the union set.
union_peaks2 <- union_peaks
for (file in peak_files) {
  cell_type_name <- gsub("_cleaned_peaks.narrowPeak", "", basename(file))
  mcols(union_peaks)[[cell_type_name]] <- 0
}
for (file in peak_files) {
  cell_type_peaks <- import(file, format = "narrowPeak")
  cell_type_name <- gsub("_cleaned_peaks.narrowPeak", "", basename(file))
  overlaps <- findOverlaps(union_peaks, cell_type_peaks)
  mcols(union_peaks)[[cell_type_name]][queryHits(overlaps)] <- 1
}
export(union_peaks, file.path(DIR_MACS, "union_peaks_with_annotations.narrowPeak"), format = "narrowPeak")

## ---- Quantify peak x cell count matrix --------------------------------------

DefaultAssay(merged_seurat4) <- "ATAC"
frags <- Fragments(merged_seurat4)
macs_count <- FeatureMatrix(fragments = frags, features = union_peaks2, cells = colnames(merged_seurat4))
save(peaks, union_peaks2, macs_count, file = MACS_COUNTS)

