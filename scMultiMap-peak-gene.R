# =============================================================================
# Peak-gene enhancer linking with scMultiMap
# -----------------------------------------------------------------------------
# Purpose : For each cell type, link accessible peaks to candidate target genes
#           (cis, +/- 1 Mb) using scMultiMap.
# Inputs  : SEURAT_OBJECT (RNA + ATAC), MACS_COUNTS (peak count matrix),
#           <DIR_MACS>/union_peaks_with_annotations.narrowPeak
# Outputs : <DIR_DIFF>/scMultiMap/
#             scMultiMap_peak_gene_results2.rds, *_results_sig2.rds
# Run     : Rscript scMultiMap-peak-gene.R
# Notes   : Core linking step only (downstream cross-celltype summary and moloc
#           intersection removed). The original script hard-coded a single cell
#           type inside the loop; this version subsets by the loop variable and
#           applies an optional per-cell-type cell cap (max_cells_per_celltype).
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(scMultiMap)
})

load(MACS_COUNTS)
load(SEURAT_OBJECT)

out_dir <- file.path(DIR_DIFF, "scMultiMap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Prepare peak + RNA assays ----------------------------------------------

DefaultAssay(merged_seurat4) <- "ATAC"
frags <- Fragments(merged_seurat4)

macs_count2 <- macs_count[!duplicated(rownames(macs_count)), ]
common_cells <- intersect(colnames(macs_count2), colnames(merged_seurat4))
macs_count2 <- macs_count2[, common_cells]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- GENOME_BUILD

merged_seurat4[["peaks"]] <- CreateChromatinAssay(
  counts = macs_count2, sep = c("-", "-"), genome = GENOME_BUILD,
  fragments = frags, annotation = annotations
)

# Restrict RNA assay to genes with matched annotation.
DefaultAssay(merged_seurat4) <- "RNA"
genes_rna <- rownames(merged_seurat4[["RNA"]])
gene_annot <- annotations[annotations$gene_name %in% genes_rna]
genes_matched <- intersect(genes_rna, gene_annot$gene_name)
gene_annot_matched <- gene_annot[match(genes_matched, gene_annot$gene_name)]
names(gene_annot_matched) <- genes_matched
merged_seurat4 <- subset(merged_seurat4, features = genes_matched)
merged_seurat4[["RNA"]]@meta.features <- as.data.frame(mcols(gene_annot_matched))

Idents(merged_seurat4) <- merged_seurat4[[CELLTYPE_COLUMN]][, 1]

# Cell-type-specific peak membership (one 0/1 column per cell type).
anno_df <- read.delim(file.path(DIR_MACS, "union_peaks_with_annotations.narrowPeak"))
anno_df$peaks <- paste0(anno_df$seqnames, "-", anno_df$start, "-", anno_df$end)
colnames(anno_df)[13:23] <- c("Astrocytes", "T Cells", "DA neurons",
  "Endothelial cell", "GABA neurons", "Glu neurons", "Microglia",
  "Macrophages", "Oligodendrocytes", "OPCs", "Pericytes")

## ---- Peak-gene linking per cell type ----------------------------------------

# Feature-selection / tractability parameters.
min_expr_quantile     <- 0      # quantile threshold on mean expression
min_expr_detected_pct <- 0.03   # min fraction of cells expressing a gene
max_cells_per_celltype <- 30000 # downsample large cell types (Inf to disable)
n_chunks <- 20                  # split pairs to limit peak memory use

cell_types <- c("Glu neurons", "GABA neurons", "Microglia", "OPCs", "T Cells",
                "DA neurons", "Endothelial cell", "Macrophages", "Pericytes",
                "Oligodendrocytes", "Astrocytes")

link_results2 <- list()
link_results_sig2 <- list()

for (ct in cell_types) {
  message("Processing: ", ct)

  colname_match <- colnames(anno_df)[grepl(ct, colnames(anno_df), ignore.case = TRUE)]
  if (length(colname_match) != 1) {
    warning("Cannot uniquely match column for cell type: ", ct)
    next
  }

  # Subset to this cell type (optionally downsampled).
  set.seed(123)
  cells_all <- colnames(merged_seurat4)[
    merged_seurat4[[CELLTYPE_COLUMN]][, 1] == ct]
  if (length(cells_all) < 50) {
    warning("Skipping ", ct, " due to low cell number.")
    next
  }
  if (is.finite(max_cells_per_celltype) && length(cells_all) > max_cells_per_celltype) {
    cells_all <- sample(cells_all, size = max_cells_per_celltype)
  }
  cell_type_subset <- merged_seurat4[, cells_all]

  # ---- Filter genes by expression ----
  expr_counts <- GetAssayData(cell_type_subset, assay = "RNA", slot = "counts")
  avg_expr <- Matrix::rowMeans(expr_counts)
  detected_pct <- Matrix::rowSums(expr_counts > 0) / ncol(expr_counts)
  expr_threshold <- quantile(avg_expr[avg_expr > 0], probs = min_expr_quantile)
  genes_passed <- names(avg_expr)[avg_expr >= expr_threshold &
                                    detected_pct >= min_expr_detected_pct]
  message(" -> ", length(genes_passed), " genes passed expression threshold")
  cell_type_subset[["RNA"]] <- subset(cell_type_subset[["RNA"]], features = genes_passed)

  # ---- Filter peaks to those active in this cell type ----
  peaks_use <- anno_df %>% filter(.data[[colname_match]] == 1) %>% pull(peaks)
  message("Found ", length(peaks_use), " peaks for ", ct)
  cell_type_subset[["peaks"]] <- subset(cell_type_subset[["peaks"]], features = peaks_use)

  # ---- Candidate cis peak-gene pairs ----
  pairs_df <- get_top_peak_gene_pairs(
    cell_type_subset,
    gene_top = nrow(cell_type_subset[["RNA"]]),
    peak_top = nrow(cell_type_subset[["peaks"]]),
    distance = 1e6, gene_assay = "RNA", peak_assay = "peaks"
  )
  message("Found ", nrow(pairs_df), " peak-gene pairs for ", ct)

  # ---- scMultiMap in chunks ----
  pairs_split <- split(pairs_df, cut(seq_len(nrow(pairs_df)), n_chunks, labels = FALSE))
  results_list <- list()
  for (i in seq_along(pairs_split)) {
    message("Running chunk ", i, " of ", length(pairs_split))
    results_list[[i]] <- scMultiMap(cell_type_subset, pairs_split[[i]],
                                    gene_assay = "RNA", peak_assay = "peaks")
  }

  all_results <- do.call(rbind, results_list)
  all_results$BH <- p.adjust(all_results$pval, method = "BH")
  link_results_sig2[[ct]] <- all_results[all_results$BH < 0.1, ]
  link_results2[[ct]] <- all_results
}

saveRDS(link_results2, file.path(out_dir, "scMultiMap_peak_gene_results.rds"))
saveRDS(link_results_sig2, file.path(out_dir, "scMultiMap_peak_gene_results_sig.rds"))
