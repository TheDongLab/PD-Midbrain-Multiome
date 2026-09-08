# =============================================================================
# eQTL preprocessing and cis-mapping inputs (FastQTL)
# -----------------------------------------------------------------------------
# Purpose : Build per-cell-type pseudobulk expression matrices per donor,
#           TMM/voom normalize, attach GENCODE gene locations, prepare SNP
#           genotype matrices, estimate PEER hidden factors, assemble FastQTL
#           BED/VCF/covariate inputs, then (shell blocks) run and combine
#           FastQTL cis-eQTL mapping.
# Inputs  : SEURAT_OBJECT, sample covariates
#           GENCODE GFF3 (GENE_GFF3), imputed genotype VCFs (GENOTYPE_DIR).
# Outputs : matrices/locations/covariates and FastQTL inputs under <DIR_EQTL>.
# Run     : Run R sections with Rscript; run shell/conda blocks separately.
#           Example processes the pooled "Neurons" cell type (ct).
# Notes   : Core steps only. QC plots and redundant intermediate writes
#           (per-chromosome matrices, raw pseudobulk, cutoff summary) were removed.
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(edgeR)
  library(rtracklayer)
  library(data.table)
  library(tidyverse)
})

GENOTYPE_DIR <- file.path("./03-caQTL-calling",
                          "genotype_reimputed_allele_aligned_to_GWAS_PD5D_eQTL")
GENE_GFF3 <- file.path(DIR_EQTL, "gencode.v47.annotation.gff3")

ct <- "Neurons"

load(SEURAT_OBJECT)

## ---- Cell-type labels, covariates, gene annotation --------------------------

DefaultAssay(merged_seurat4) <- "RNA"
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- GENOME_BUILD

merged_seurat4$celltype.WNN.after.combine.update3 <- gsub(" ", "", merged_seurat4$celltype.WNN.after.combine.update)
old_ct <- merged_seurat4$celltype.WNN.after.combine.update3
is_neuron <- grepl("(?i)neuron", old_ct)
merged_seurat4$celltype_merged_neuron <- ifelse(is_neuron, "Neurons", old_ct)

Idents(merged_seurat4) <- "celltype_merged_neuron"
merged_seurat4$Sample2 <- substr(merged_seurat4$Sample, 1, 6)

anno.samp2 <- merged_seurat4@meta.data
PCs <- read.delim("./00-data/covariants_100_samples_with_first_three_sva_factors.txt")
anno.samp <- merge(anno.samp2, PCs, by.x = "Sample2", by.y = "library_id")

# Gene locations from GENCODE (keep uniquely named genes).
gff_data <- import(GENE_GFF3, format = "gff3")
gene_data <- subset(gff_data, type == "gene")
gene_annotation <- data.frame(
  gene_name = sapply(mcols(gene_data)$gene_name, function(x) ifelse(is.null(x), NA, x)),
  chr = seqnames(gene_data), start = start(gene_data), end = end(gene_data),
  strand = strand(gene_data))
gene_annotation_unique <- gene_annotation %>% group_by(gene_name) %>%
  dplyr::filter(n() == 1) %>% ungroup()

## ---- Pseudobulk expression, filtering, TMM/voom -----------------------------

print(ct)
seurat_ct <- subset(merged_seurat4, idents = ct)
pseudobulk_counts <- AggregateExpression(seurat_ct, group.by = "Sample2", assays = "RNA",
                                         slot = "counts", return.seurat = FALSE)$RNA
mtx <- pseudobulk_counts
low_5_percent <- quantile(colSums(mtx), 0.05)

filter <- apply(mtx, 2, function(x) sum(x) >= min(1000, low_5_percent))
mtx <- mtx[, filter]
filter <- apply(mtx, 1, function(x) length(x[x >= 2]) >= ncol(mtx) * 0.3)
mtx <- mtx[filter, ]
print(dim(mtx))

y <- DGEList(counts = mtx)
y <- calcNormFactors(y)
y <- voom(y)
y.mtx <- y$E
common_samples <- intersect(anno.samp$Sample2, colnames(y.mtx))
y.mtx <- y.mtx[, common_samples]
write.table(data.frame(gene = rownames(y.mtx), y.mtx),
            file.path(DIR_EQTL, "geneMatrix", paste0(ct, ".TMM.voom.normed.mtx")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# Gene location BED.
gene_annotation_all <- gene_annotation_unique[gene_annotation_unique$gene_name %in% rownames(y.mtx), ]
write.table(gene_annotation_all, file.path(DIR_EQTL, "geneLoc", paste0(ct, ".gene.hg38.bed")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# Sample covariate table (includes RIN; sex coded 0/1; + SVA factors).
sample_metadata <- anno.samp[, c("Sample2", "Age", "Sex", "PMI.x", "RIN.x", "BatchLib",
                                 "Sva_factor1", "Sva_factor2", "Sva_factor3")]
sample_metadata <- unique(sample_metadata)
sample_metadata$Sex <- ifelse(sample_metadata$Sex == "F", 0, 1)
sample_metadata <- t(sample_metadata[match(colnames(y.mtx), sample_metadata$Sample2), ])
sample_metadata <- gsub(" ", "", sample_metadata)
write.table(sample_metadata, file.path(DIR_EQTL, "covariates", paste0(ct, ".sample.metadata.txt")),
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

## ---- SNP genotype matrices (dosage encoding) --------------------------------

inDir <- GENOTYPE_DIR
outDir <- file.path(DIR_EQTL, "snpMatrix")
for (i in 1:22) {
  vcf <- fread(file.path(inDir, paste0("GWAS.chr", i, ".allele_align.all_imputed_20240924.fam.PD5D_include_midBrain.vcf")))
  colnames(vcf)[1:2] <- c("chr", "pos")
  snpLoc <- vcf[, 1:3]
  snpLoc$chr <- paste0("chr", snpLoc$chr)
  snpLoc <- data.frame(snp = paste("snp", snpLoc$chr, snpLoc$pos, snpLoc$ID, sep = "_"), snpLoc)
  write.table(snpLoc[, 1:3], file.path(DIR_EQTL, "snpLoc", paste0("chr", i, "-snp-location-renew.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  rownames(vcf) <- paste("snp", paste0("chr", vcf$chr), vcf$pos, vcf$ID, sep = "_")
  vcf1 <- vcf[, 10:ncol(vcf)]
  vcf3 <- vcf1 %>%
    mutate_all(funs(str_replace(., "1/1", "2"))) %>%
    mutate_all(funs(str_replace(., "0/1", "1"))) %>%
    mutate_all(funs(str_replace(., "1/0", "1"))) %>%
    mutate_all(funs(str_replace(., "0/0", "0"))) %>%
    mutate_all(funs(str_replace(., "./.", "-1")))
  colnames(vcf3) <- substr(colnames(vcf3), 1, 6)
  vcf3 <- data.frame(id = rownames(vcf), vcf3)
  common_samples <- intersect(colnames(vcf3), unique(anno.samp$Sample2))
  vcf3 <- vcf3[, c("id", common_samples)]
  write.table(vcf3, file.path(outDir, paste0("chr", i, "-snp-matrix-renew.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

## ---- PEER hidden factors ----------------------------------------------------
# Run inside the PEER conda environment (see shell block below).

library(peer)

peer_input_file <- file.path(DIR_EQTL, "geneMatrix", paste0(ct, ".TMM.voom.normed.mtx"))
if (!file.exists(peer_input_file)) stop(paste("Transformed matrix file not found for cell type", ct))
peer_mtx <- read.table(peer_input_file, header = TRUE, sep = "\t", row.names = 1)
peer_mtx <- t(as.matrix(peer_mtx))

model <- PEER()
n_hidden <- 20
PEER_setNk(model, n_hidden)
PEER_setPhenoMean(model, peer_mtx)
PEER_setNmax_iterations(model, 1000)

covariates_file <- file.path(DIR_EQTL, "covariates", paste0(ct, ".sample.metadata.txt"))
if (file.exists(covariates_file)) {
  sample_covariates <- read.table(covariates_file, header = TRUE, sep = "\t", row.names = 1)
  sample_covariates <- t(as.matrix(sample_covariates))
  PEER_setCovariates(model, sample_covariates)
}

PEER_update(model)
peer_factors <- PEER_getX(model)
rownames(peer_factors) <- rownames(peer_mtx)
colnames(peer_factors) <- c(colnames(sample_covariates), paste0("InferredCov", 1:20))
peer_residuals <- t(PEER_getResiduals(model))
rownames(peer_residuals) <- colnames(peer_mtx)
colnames(peer_residuals) <- rownames(peer_mtx)

write.table(peer_factors, file.path(DIR_EQTL, "PEER", paste0(ct, ".PEER_covariates.txt")),
            row.names = TRUE, sep = "\t", quote = FALSE)
write.table(peer_residuals, file.path(DIR_EQTL, "PEER", paste0(ct, ".PEER_residuals.txt")),
            row.names = TRUE, sep = "\t", quote = FALSE)

alpha <- PEER_getAlpha(model)
alpha_df <- data.frame(Factor = 1:length(alpha), Alpha = alpha)
write.table(alpha_df, file.path(DIR_EQTL, "PEER", paste0(ct, ".PEER_alpha.txt")),
            row.names = TRUE, sep = "\t", quote = FALSE)

pdf(file.path(DIR_EQTL, "PEER", paste0(ct, "_alpha50_plot.pdf")), height = 8, width = 8)
plot(alpha_df$Factor, alpha_df$Alpha, type = "b", pch = 19, col = "blue",
     main = paste("Alpha Values for PEER Factors -", ct),
     xlab = "Factor", ylab = "Alpha (Inverse Relevance)")
dev.off()

## ---- Assemble FastQTL BED / VCF / covariate inputs --------------------------

peak_loc_dir <- file.path(DIR_EQTL, "geneLoc")
out_dir_bed  <- file.path(DIR_EQTL, "Fastqtl/Genebed")
out_dir_snp  <- file.path(DIR_EQTL, "Fastqtl/SNPvcf")
inDir <- GENOTYPE_DIR

ct_file <- file.path(DIR_EQTL, "PEER", paste0(ct, ".PEER_residuals.txt"))
peak_matrix <- fread(ct_file)
peak_loc <- fread(file.path(peak_loc_dir, paste0(ct, ".gene.hg38.bed")), header = TRUE)
colnames(peak_matrix)[1] <- "gene_name"
bed <- peak_loc %>% left_join(peak_matrix, by = c("gene_name" = "gene_name")) %>%
  dplyr::select(chr, start, end, gene_name, everything()) %>% arrange(chr, start)
colnames(bed)[1] <- "chr"

calculate_maf <- function(genotypes) {
  alleles <- unlist(strsplit(as.character(genotypes), split = "[/]"))
  allele_counts <- table(alleles)
  min(allele_counts / sum(allele_counts))
}

for (i in 1:22) {
  print(i)
  vcf <- fread(file.path(inDir, paste0("GWAS.chr", i, ".allele_align.all_imputed_20240924.fam.PD5D_include_midBrain.vcf")))
  vcf1 <- vcf
  colnames(vcf1)[10:ncol(vcf1)] <- substr(colnames(vcf1)[10:ncol(vcf1)], 1, 6)
  common_samples <- intersect(colnames(bed)[6:ncol(bed)], colnames(vcf1))
  vcf2 <- as.data.frame(vcf1)[, c(colnames(vcf1)[1:9], common_samples)]
  vcf2[, 1] <- paste0("chr", vcf2[, 1])
  genotypes <- vcf2[, 10:ncol(vcf2)]
  vcf2$MAF <- apply(genotypes, 1, calculate_maf)
  filtered_vcf2 <- vcf2[vcf2$MAF > 0.05, ]
  filtered_vcf2 <- filtered_vcf2[, -ncol(filtered_vcf2)]
  write.table(filtered_vcf2, file.path(out_dir_snp, paste0(ct, ".chr", i, "_snp_matrix.vcf")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

bed <- as.data.frame(bed)[, c(colnames(bed)[1:4], common_samples)]
write.table(bed, file.path(out_dir_bed, paste0(ct, "_geneMatrix.hg38.bed")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
for (chr in unique(bed$chr)) {
  bed_chr <- bed %>% dplyr::filter(chr == !!chr)
  write.table(bed_chr, file.path(out_dir_bed, paste0(ct, "_geneMatrix.hg38_", chr, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# FastQTL covariate file = PEER covariates (which already include sample covariates).
peer_file <- file.path(DIR_EQTL, "PEER", paste0(ct, ".PEER_covariates.txt"))
peer_factors <- t(read.table(peer_file, header = TRUE, sep = "\t", row.names = 1))
allcovaraints <- peer_factors[, common_samples]
allcovaraints <- data.frame(ID = rownames(allcovaraints), allcovaraints)
write.table(allcovaraints, file.path(DIR_EQTL, "Fastqtl/Covariates", paste0(ct, ".Covariates.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# =============================================================================
# Shell / conda blocks (run separately, not as R)
# =============================================================================
#
# ---- 1. PEER conda environment (same as caQTL step) ----
# conda create --name peer python=3.5.4 && conda activate peer
# conda config --add channels conda-forge && conda config --add channels r
# conda install zlib=1.2.8 r=3.4.1 r-peer libgfortran=3 conda-forge::r-qtl
#
# ---- 2. Sort, bgzip, tabix VCF and BED files (05-eQTL-calling2/Fastqtl/{SNPvcf,Genebed}) ----
#   grep '^#' $vcf > $sorted; grep -v '^#' $vcf | LC_ALL=C sort -k1,1 -k2,2n >> $sorted
#   bgzip $sorted; tabix -p vcf ${sorted}.gz
#   { head -n1 $bed | sed 's/^/#/'; tail -n +2 $bed | LC_ALL=C sort -k1,1 -k2,2n; } > $sorted_bed
#   bgzip $sorted_bed; tabix -p bed ${sorted_bed}.gz
#
# ---- 3. FastQTL commands (covariates already regressed into PEER residuals) ----
# Permutation : fastQTL.static --vcf VCF --bed BED --out PERM --window 1e6 --permute 1000 10000--region chrN
# Nominal     : fastQTL.static --vcf VCF --bed BED --out NOM  --window 1e6 --region chrN
# Run e.g.    : ParaFly -c ./Scripts/Neurons_fastqtl_commands.txt -CPU 10
#
# ---- 4. Combine per-chromosome results ----
# Permutation header:
#   gene_id num_var beta_shape1 beta_shape2 true_df pval_true_df variant_id
#   tss_distance ma_samples ma_count maf ref_factor pval_nominal slope slope_se
#   pval_perm pval_beta
# =============================================================================
