# =============================================================================
# Step 09 - Single-cell differential expression (NEBULA) across diagnosis
# -----------------------------------------------------------------------------
# Purpose : Per cell type, test differential gene expression between diagnosis
#           groups (PD/ILBD/HC) with NEBULA negative-binomial mixed models,
#           adjusting for technical variation via pseudobulk RUVr factors.
# Method  : scToNeb -> pseudobulk DGEList + RUVr (k=5) -> NEBULA (NBGMM)
#           covariates: Diagnosis + Age + Sex + PMI + BatchLib + W_1..W_5.
# Inputs  : SEURAT_OBJECT (RNA assay, Diagnosis/Age/Sex/PMI/BatchLib metadata)
# Outputs : <DIR_DIFF>/degs/nebula_sc/  (per cell type x comparison)
#             deg_sc_*.{Rds,tsv} and all_deg_results.rds
# Run     : Rscript PD-differential-nebula-deg.R
# Notes   : Core DEG step only. Downstream GO over-representation and GSEA
#           enrichment were removed. Interactive per-iteration overrides
#           (e.g. ct <- "Glu neurons") were removed so the loop iterates over
#           all cell types as intended.
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(nebula)
  library(RUVSeq)
  library(edgeR)
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

load(SEURAT_OBJECT)
DefaultAssay(merged_seurat4) <- "RNA"

out_dir <- file.path(DIR_DIFF, "degs/nebula_sc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comparison_list <- COMPARISONS
cell_types <- c("GABA neurons", "DA neurons", "Microglia", "T Cells",
                "Endothelial cell", "Astrocytes", "Macrophages", "Pericytes",
                "OPCs", "Oligodendrocytes", "Glu neurons")

deg_results <- list()

Idents(merged_seurat4) <- merged_seurat4[[CELLTYPE_COLUMN]][, 1]

for (ct in cell_types) {
  message("Processing cell type: ", ct)
  sub_obj <- subset(merged_seurat4, idents = ct)

  for (comp_name in names(comparison_list)) {
    groups <- comparison_list[[comp_name]]
    message("  Comparison: ", comp_name)

    sub_obj_comp <- subset(sub_obj, subset = Diagnosis %in% groups)
    if (ncol(sub_obj_comp) < 50) {
      warning("Too few cells for ", ct, " in ", comp_name)
      next
    }

    # ---- Convert to NEBULA input ----
    seuratdata <- scToNeb(
      obj = sub_obj_comp,
      assay = "RNA",
      id = "Sample",
      pred = c("Diagnosis", "Age", "Sex", "PMI", "BatchLib", "Sample"),
      offset = "nCount_RNA"
    )

    valid_idx <- complete.cases(seuratdata$pred)
    if (sum(valid_idx) < 100) {
      warning("Too few valid cells after filtering for ", ct, " in ", comp_name)
      next
    }
    seuratdata$count  <- seuratdata$count[, valid_idx, drop = FALSE]
    seuratdata$offset <- seuratdata$offset[valid_idx]
    seuratdata$id     <- seuratdata$id[valid_idx]
    seuratdata$pred   <- seuratdata$pred[valid_idx, , drop = FALSE]

    seuratdata$pred$Diagnosis <- relevel(factor(seuratdata$pred$Diagnosis),
                                         ref = groups[2])

    # ---- Pseudobulk RUVr to estimate unwanted-variation factors ----
    message("Generating pseudobulk for RUV")
    counts <- seuratdata$count
    ids <- seuratdata$id
    samples <- unique(ids)

    pb_mat <- do.call(cbind, lapply(samples, function(s) {
      Matrix::rowSums(counts[, ids == s, drop = FALSE])
    }))
    colnames(pb_mat) <- samples

    dge_pb <- DGEList(counts = pb_mat)
    keep_genes_pb <- rowSums(edgeR::cpm(dge_pb) > 1) >= 5
    dge_pb <- dge_pb[keep_genes_pb, , keep.lib.sizes = FALSE]
    dge_pb <- calcNormFactors(dge_pb)

    pseudo_meta <- seuratdata$pred %>% distinct()
    rownames(pseudo_meta) <- pseudo_meta$Sample
    pseudo_meta <- pseudo_meta[colnames(pb_mat), , drop = FALSE]

    design_pb <- model.matrix(~ Diagnosis + Age + Sex + PMI + BatchLib,
                              data = pseudo_meta)
    dge_pb <- estimateGLMCommonDisp(dge_pb, design_pb)
    dge_pb <- estimateGLMTagwiseDisp(dge_pb, design_pb)
    fit_pb <- glmFit(dge_pb, design_pb)
    residuals_pb <- residuals(fit_pb, type = "deviance")

    ruv_pb <- RUVr(round(dge_pb$counts), rownames(dge_pb), k = 5,
                   res = residuals_pb)

    ruv_df <- as.data.frame(ruv_pb$W)
    colnames(ruv_df) <- paste0("W_", 1:5)
    ruv_df$Sample <- pseudo_meta$Sample
    seuratdata$pred <- left_join(
      seuratdata$pred %>% mutate(Sample = seuratdata$id),
      ruv_df, by = "Sample"
    )

    full_design1 <- model.matrix(
      ~ Diagnosis + Age + Sex + PMI + BatchLib + W_1 + W_2 + W_3 + W_4 + W_5,
      data = seuratdata$pred
    )

    # ---- NEBULA negative-binomial mixed model ----
    neb1 <- nebula(
      count = seuratdata$count,
      id = seuratdata$id,
      pred = full_design1,
      offset = seuratdata$offset,
      model = "NBGMM"
    )
    if (is.null(neb1) || !"summary" %in% names(neb1)) next

    p_col  <- paste0("p_Diagnosis", groups[1])
    fc_col <- paste0("logFC_Diagnosis", groups[1])
    res_name <- paste(ct, comp_name, sep = "_")
    neb_df <- na.omit(neb1$summary)
    if (!(p_col %in% colnames(neb_df))) next

    neb_df$FDR_BH <- p.adjust(neb_df[[p_col]], method = "BH")
    deg_results[[res_name]] <- neb_df

    up_genes <- neb_df$gene[neb_df$FDR_BH <= 0.05 & neb_df[[fc_col]] >= 0.5]
    down_genes <- neb_df$gene[neb_df$FDR_BH <= 0.05 & neb_df[[fc_col]] <= -0.5]
    message("    Up: ", length(up_genes), " | Down: ", length(down_genes))

    saveRDS(neb1, file.path(out_dir, paste0("deg_sc_", gsub(" ", "", res_name), ".Rds")))
    write.table(neb_df, file.path(out_dir, paste0("deg_sc_", gsub(" ", "", res_name), ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# ---- Save master DEG result list ----
saveRDS(deg_results, file.path(out_dir, "all_deg_results.rds"))
