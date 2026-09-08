# =============================================================================
# eRegulon gene-set GSEA core (per cell type)
# -----------------------------------------------------------------------------
# Purpose : Build eRegulon -> target-gene sets from SCENIC+ direct-eRegulon
#           exports (optionally restricted to a heatmap-selected allowlist),
#           then run fgsea (via clusterProfiler::GSEA) on per-comparison NEBULA
#           DEG rankings and write the enrichment table + GSEA objects.
# Inputs  : SCENIC+ r_export (all_peaks_eRegulon_direct.tsv,
#           heatmap_selected_eregulons_gene.tsv) under
#           <DIR_SCENIC>/<CELL_TYPE>/pd_progression_eregulons_activator/r_export/,
#           NEBULA DEGs under <DIR_DIFF>/degs/nebula_sc.
# Outputs : DEG_eRegulon_GSEA_all_terms.tsv and
#           DEG_eRegulon_GSEA_objects_by_group.rds under the cell type's
#           <CELL_TYPE>_eregulon_network_gsea_fgsea_activator output dir.
# Run     : Rscript eRegulon-GSEA-core.R   (default CELL_TYPE = Microglia)
# Notes   : Core GSEA step only. Downstream plotting parameters/helpers and the
#           diagnosis-specific RSS prefilter (used by 13-2..13-4) were removed.
#           Template is Microglia; 13-generate-celltype-pipeline.R generates the
#           equivalent scripts for other cell types.
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(clusterProfiler)
})

# Cell type to process (Microglia template; change to generate for others).
CELL_TYPE <- "Microglia"
cell_slug <- tolower(CELL_TYPE)
focus_celltype <- CELL_TYPE

# PD progression activator arm: r_export + heatmap-selected eRegulon list.
scenic_eregulon_profile <- "pd_progression_eregulons_activator"
r_export_dir <- file.path(DIR_SCENIC, CELL_TYPE, scenic_eregulon_profile, "r_export")
scenic_dir <- file.path(DIR_DIFF, "scenicplus2")
ereg_file <- file.path(r_export_dir, "all_peaks_eRegulon_direct.tsv")
# If present, GSEA term2gene is restricted to eRegulon_name values in this table.
eregulon_gsea_allowlist_file <- file.path(r_export_dir, "heatmap_selected_eregulons_gene.tsv")

out_base <- file.path(DIR_SCENIC, CELL_TYPE, paste0(cell_slug, "_eregulon_network_gsea_fgsea_activator"))
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
gsea_tsv <- file.path(out_base, "DEG_eRegulon_GSEA_all_terms.tsv")
gsea_objects_rds <- file.path(out_base, "DEG_eRegulon_GSEA_objects_by_group.rds")

# Parameters.
min_gssize <- 10L
max_gssize <- 500L
min_eregulon_gene_count <- 1L
deg_fdr_threshold <- 0.1
reproducible_seed <- 20260513L

stopifnot(file.exists(ereg_file))

normalize_ct <- function(x) gsub(" ", "", as.character(x))
normalize_ereg_id <- function(x) gsub("_(\\d+g)_?$", "", as.character(x), perl = TRUE)

## ----- eRegulon triplets -> gene sets ----------------------------------------

ereg_targets <- fread(ereg_file) %>%
  transmute(eRegulon = as.character(eRegulon_name), TF = as.character(TF), Gene = toupper(as.character(Gene))) %>%
  filter(!is.na(Gene), Gene != "", !is.na(TF), TF != "", !is.na(eRegulon), eRegulon != "") %>%
  distinct()

if (!is.na(eregulon_gsea_allowlist_file) && nzchar(eregulon_gsea_allowlist_file)) {
  stopifnot(file.exists(eregulon_gsea_allowlist_file))
  allow_tbl <- fread(eregulon_gsea_allowlist_file)
  idcol <- if ("eRegulon" %in% names(allow_tbl)) "eRegulon" else names(allow_tbl)[1L]
  raw_ids <- trimws(as.character(allow_tbl[[idcol]]))
  raw_ids <- raw_ids[!is.na(raw_ids) & nzchar(raw_ids) & raw_ids != "eRegulon"]
  allow_ids <- unique(normalize_ereg_id(gsub("_\\(\\d+g\\)$", "", raw_ids, perl = TRUE)))
  allow_ids <- allow_ids[nzchar(allow_ids)]
  n_er_all <- dplyr::n_distinct(ereg_targets$eRegulon)
  ereg_targets <- dplyr::filter(ereg_targets, eRegulon %in% allow_ids)
  message("eRegulon GSEA allowlist (", basename(eregulon_gsea_allowlist_file), "): ",
          length(allow_ids), " ID(s) listed -> ", dplyr::n_distinct(ereg_targets$eRegulon), " / ", n_er_all,
          " eRegulons kept from ", basename(ereg_file), ".")
  if (nrow(ereg_targets) == 0L) stop("Allowlist removed all eRegulon rows; check IDs vs eRegulon_name in ", ereg_file)
}

ereg_sets <- ereg_targets %>%
  group_by(eRegulon, TF) %>%
  summarise(Target_genes = list(unique(Gene)), setSize = dplyr::n_distinct(Gene), .groups = "drop")

n_er_before <- dplyr::n_distinct(ereg_targets$eRegulon)
ereg_sets <- dplyr::filter(ereg_sets, .data$setSize >= min_eregulon_gene_count)
if (nrow(ereg_sets) == 0L) {
  stop("No eRegulons left after requiring >=", min_eregulon_gene_count,
       " distinct targets (had ", n_er_before, " eRegulons before filter).")
}
ereg_targets <- dplyr::semi_join(ereg_targets, dplyr::distinct(ereg_sets, eRegulon, TF), by = c("eRegulon", "TF"))
message("eRegulon target filter: kept ", dplyr::n_distinct(ereg_targets$eRegulon), " / ",
        n_er_before, " eRegulons with >=", min_eregulon_gene_count, " genes.")

term2gene_ereg <- ereg_targets %>% transmute(term = eRegulon, gene = Gene) %>% distinct()

## ----- DEG long for the focus cell type --------------------------------------

deg_dir <- file.path(scenic_dir, "..", "degs", "nebula_sc")
deg_rds_candidates <- c(file.path(deg_dir, "all_deg_results+Neurons.rds"), file.path(deg_dir, "all_deg_results.rds"))
deg_rds <- deg_rds_candidates[file.exists(deg_rds_candidates)][1]

load_deg_long_from_rds <- function(path_rds) {
  deg_results <- readRDS(path_rds)
  if (!is.list(deg_results) || length(deg_results) == 0) return(NULL)
  bind_rows(lapply(names(deg_results), function(res_name) {
    deg_df <- deg_results[[res_name]]
    if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
    split_name <- strsplit(res_name, "_")[[1]]
    if (length(split_name) < 3) return(NULL)
    cell_type <- normalize_ct(split_name[1])
    comparison <- paste(split_name[-1], collapse = "_")
    groups <- strsplit(comparison, "_vs_")[[1]]
    if (length(groups) != 2) return(NULL)
    p_col <- paste0("p_Diagnosis", groups[1]); fc_col <- paste0("logFC_Diagnosis", groups[1])
    if (!(p_col %in% names(deg_df)) || !(fc_col %in% names(deg_df))) return(NULL)
    if (!("FDR_BH" %in% names(deg_df))) deg_df$FDR_BH <- p.adjust(deg_df[[p_col]], method = "BH")
    deg_df %>%
      transmute(CellType = cell_type, Comparison = comparison, Gene = toupper(as.character(gene)),
                logFC = as.numeric(.data[[fc_col]]), PValue = as.numeric(.data[[p_col]]), FDR_BH = as.numeric(FDR_BH)) %>%
      filter(!is.na(Gene), Gene != "", is.finite(logFC), is.finite(PValue), PValue > 0)
  })) %>% filter(CellType == normalize_ct(CELL_TYPE))
}

load_deg_long_from_tsv <- function(dir_path) {
  files <- list.files(dir_path, pattern = paste0("^deg_sc_", CELL_TYPE, "_.*\\.tsv$"), full.names = TRUE)
  if (length(files) == 0) return(NULL)
  bind_rows(lapply(files, function(f) {
    comp <- sub("\\.tsv$", "", sub(paste0("^deg_sc_", CELL_TYPE, "_"), "", basename(f)))
    groups <- strsplit(comp, "_vs_")[[1]]
    if (length(groups) != 2) return(NULL)
    grp <- groups[1]
    deg_df <- fread(f)
    p_col <- paste0("p_Diagnosis", grp); fc_col <- paste0("logFC_Diagnosis", grp)
    if (!(p_col %in% names(deg_df)) || !(fc_col %in% names(deg_df))) return(NULL)
    gene_col <- if ("gene" %in% names(deg_df)) "gene" else if ("Gene" %in% names(deg_df)) "Gene" else return(NULL)
    if (!("FDR_BH" %in% names(deg_df))) deg_df[, FDR_BH := stats::p.adjust(.data[[p_col]], method = "BH")]
    deg_df %>%
      transmute(CellType = CELL_TYPE, Comparison = comp, Gene = toupper(as.character(.data[[gene_col]])),
                logFC = as.numeric(.data[[fc_col]]), PValue = as.numeric(.data[[p_col]]), FDR_BH = as.numeric(FDR_BH)) %>%
      filter(!is.na(Gene), Gene != "", is.finite(logFC), is.finite(PValue), PValue > 0)
  }))
}

if (!is.na(deg_rds)) {
  message("Loading ", CELL_TYPE, " DEG from RDS: ", deg_rds)
  deg_long <- load_deg_long_from_rds(deg_rds)
} else {
  message("No all_deg_results*.RDS; building ", CELL_TYPE, " DEG from TSVs in ", deg_dir)
  deg_long <- load_deg_long_from_tsv(deg_dir)
}

if (is.null(deg_long) || nrow(deg_long) == 0) {
  stop("No ", CELL_TYPE, " DEG rows. Provide ", paste(basename(deg_rds_candidates), collapse = " or "),
       " under ", deg_dir, ", or deg_sc_", CELL_TYPE, "_*.tsv files.")
}

deg_ranked <- deg_long %>%
  group_by(CellType, Comparison, Gene) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>% ungroup()

rank_keys <- deg_ranked %>% filter(CellType == normalize_ct(focus_celltype)) %>% distinct(CellType, Comparison)

build_gene_list <- function(df_one) {
  df_one <- df_one %>% mutate(rank_score = logFC * (-log10(PValue)))
  gl <- df_one$rank_score; names(gl) <- df_one$Gene
  gl <- gl[is.finite(gl)]; gl <- sort(gl, decreasing = TRUE)
  dup <- duplicated(names(gl)); if (any(dup)) gl <- gl[!dup]
  gl
}

## ----- GSEA per comparison ---------------------------------------------------

all_tables <- list()
gsea_objects_by_group <- list()

for (k in seq_len(nrow(rank_keys))) {
  ct <- rank_keys$CellType[k]; comp <- rank_keys$Comparison[k]
  key <- paste(ct, comp, sep = "_")
  df_one <- deg_ranked %>% filter(CellType == ct, Comparison == comp) %>% select(Gene, logFC, PValue)
  gl <- build_gene_list(df_one)
  set.seed(reproducible_seed + k)
  gso <- tryCatch(
    clusterProfiler::GSEA(geneList = gl, exponent = 1, minGSSize = min_gssize, maxGSSize = max_gssize,
                          pAdjustMethod = "BH", pvalueCutoff = 1, TERM2GENE = term2gene_ereg, eps = 0,
                          scoreType = "std", by = "fgsea", seed = TRUE, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(gso)) next
  gsea_objects_by_group[[key]] <- gso
  tb <- as.data.frame(gso)
  if (nrow(tb) == 0) next
  tb <- tb %>%
    transmute(ID = as.character(ID), Description = as.character(Description), setSize = as.integer(setSize),
              enrichmentScore = as.numeric(enrichmentScore), NES = as.numeric(NES), pvalue = as.numeric(pvalue),
              p.adjust = as.numeric(p.adjust), qvalue = as.numeric(qvalue), CellType = ct, Comparison = comp,
              Direction = ifelse(NES > 0, "Up", "Down")) %>%
    left_join(ereg_sets %>% transmute(ID = eRegulon, TF = TF, setSize2 = setSize) %>% distinct(ID, TF, setSize2), by = "ID") %>%
    mutate(setSize = ifelse(is.na(setSize), setSize2, setSize)) %>% select(-setSize2)
  all_tables[[key]] <- tb
}

if (length(all_tables) == 0) stop("No GSEA enrichment results for ", CELL_TYPE, " (check eRegulon sets and DEG inputs).")

gsea_all <- bind_rows(all_tables) %>%
  arrange(CellType, Comparison, Direction, p.adjust, desc(abs(NES)), ID)

## ----- Overlap of eRegulon targets with DE up/down genes ---------------------

tf_l <- stats::setNames(ereg_sets$Target_genes, ereg_sets$eRegulon)

up_map <- deg_long %>% filter(FDR_BH <= deg_fdr_threshold, logFC > 0) %>%
  group_by(CellType, Comparison) %>% summarise(up_genes = list(unique(Gene)), .groups = "drop") %>%
  mutate(k = paste(CellType, Comparison, sep = "|"))
up_l <- stats::setNames(up_map$up_genes, up_map$k)

down_map <- deg_long %>% filter(FDR_BH <= deg_fdr_threshold, logFC < 0) %>%
  group_by(CellType, Comparison) %>% summarise(down_genes = list(unique(Gene)), .groups = "drop") %>%
  mutate(k = paste(CellType, Comparison, sep = "|"))
down_l <- stats::setNames(down_map$down_genes, down_map$k)

gsea_all$nOverlap <- as.integer(mapply(function(ct, comp, tid, nes) {
  k <- paste(ct, comp, sep = "|"); tg <- tf_l[[tid]]
  if (is.null(tg)) return(0L)
  if (nes > 0) { ug <- up_l[[k]]; if (is.null(ug)) return(0L); length(intersect(tg, ug)) }
  else { dg <- down_l[[k]]; if (is.null(dg)) return(0L); length(intersect(tg, dg)) }
}, gsea_all$CellType, gsea_all$Comparison, gsea_all$ID, gsea_all$NES, SIMPLIFY = TRUE))

fwrite(as.data.table(gsea_all), gsea_tsv, sep = "\t")
saveRDS(gsea_objects_by_group, gsea_objects_rds)
message("Wrote GSEA enrichment table: ", normalizePath(gsea_tsv, mustWork = FALSE))
message("Wrote GSEA objects RDS: ", normalizePath(gsea_objects_rds, mustWork = FALSE))
message("Done core phase. Wrote GSEA table and objects under out_base.")
