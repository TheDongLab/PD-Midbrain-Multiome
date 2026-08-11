# =============================================================================
# Step 05 - caQTL preprocessing and cis-mapping inputs (FastQTL)
# -----------------------------------------------------------------------------
# Purpose : Build per-cell-type pseudobulk accessibility matrices per donor,
#           TMM/voom normalize, prepare SNP genotype matrices, estimate PEER
#           hidden factors, assemble FastQTL BED/VCF/covariate inputs, then
#           (in the shell blocks) run and combine FastQTL cis-caQTL mapping.
# Inputs  : SEURAT_OBJECT, MACS_COUNTS, sample covariates
#           (./00-data/covariants_100_samples_with_first_three_sva_factors.txt),
#           imputed genotype VCFs (GENOTYPE_DIR).
# Outputs : matrices/locations/covariates and FastQTL inputs under <DIR_CAQTL>.
# Run     : Run the R sections with Rscript; run the shell/conda blocks at the
#           bottom separately. The example below processes the pooled "Neurons"
#           cell type (ct); loop over cell types as needed.
# Notes   : Core steps only. QC plots and redundant intermediate writes
#           (per-chromosome matrices, raw pseudobulk, cutoff summary, PEER alpha)
#           were removed.
# =============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(edgeR)
  library(data.table)
  library(tidyverse)
})

# Imputed genotype VCFs (GWAS-aligned), one per chromosome.
GENOTYPE_DIR <- file.path("./03-caQTL-calling",
                          "genotype_reimputed_allele_aligned_to_GWAS_PD5D_eQTL")

# Cell type to process (loop over cell types to run all).
ct <- "Neurons"

load(MACS_COUNTS)
load(SEURAT_OBJECT)

## ---- Add MACS peak assay and define a pooled-neuron cell-type label ----------

DefaultAssay(merged_seurat4) <- "ATAC"
frags <- Fragments(merged_seurat4)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- GENOME_BUILD

merged_seurat4[["peaks"]] <- CreateChromatinAssay(counts = macs_count, sep = c("-", "-"),
                                                  genome = GENOME_BUILD, fragments = frags,
                                                  annotation = annotations)
DefaultAssay(merged_seurat4) <- "peaks"

merged_seurat4$celltype.WNN.after.combine.update3 <- gsub(" ", "", merged_seurat4$celltype.WNN.after.combine.update)
old_ct <- merged_seurat4$celltype.WNN.after.combine.update3
is_neuron <- grepl("(?i)neuron", old_ct)         # GABA/Glu/DA neurons -> "Neurons"
merged_seurat4$celltype_merged_neuron <- ifelse(is_neuron, "Neurons", old_ct)

Idents(merged_seurat4) <- "celltype_merged_neuron"
merged_seurat4$Sample2 <- substr(merged_seurat4$Sample, 1, 6)

# Sample covariates + SVA factors.
anno.samp2 <- merged_seurat4@meta.data
PCs <- read.delim("./00-data/covariants_100_samples_with_first_three_sva_factors.txt")
anno.samp <- merge(anno.samp2, PCs, by.x = "Sample2", by.y = "library_id")

## ---- Pseudobulk accessibility matrix, filtering, TMM/voom -------------------

print(ct)
seurat_ct <- subset(merged_seurat4, idents = ct)

pseudobulk_counts <- AggregateExpression(seurat_ct, group.by = "Sample2", assays = "peaks",
                                         slot = "counts", return.seurat = FALSE)$peaks
mtx <- pseudobulk_counts
low_5_percent <- quantile(colSums(mtx), 0.05)

# Sample- and peak-level filtering.
filter <- apply(mtx, 2, function(x) sum(x) >= min(1000, low_5_percent))
mtx <- mtx[, filter]
filter <- apply(mtx, 1, function(x) length(x[x >= 2]) >= ncol(mtx) * 0.3)
mtx <- mtx[filter, ]
print(dim(mtx))
peaks <- rownames(mtx)

# TMM normalization + voom.
y <- DGEList(counts = mtx)
y <- calcNormFactors(y)
y <- voom(y)
y.mtx <- y$E
common_samples <- intersect(anno.samp$Sample2, colnames(y.mtx))
y.mtx <- y.mtx[, common_samples]
write.table(data.frame(peakid = rownames(y.mtx), y.mtx),
            file.path(DIR_CAQTL, "peakMatrix", paste0(ct, ".TMM.voom.normed.mtx")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# Peak location BED.
peak_annotations <- data.frame(peakid = peaks,
                               chr = sapply(strsplit(peaks, "-"), `[`, 1),
                               start = as.integer(sapply(strsplit(peaks, "-"), `[`, 2)),
                               end = as.integer(sapply(strsplit(peaks, "-"), `[`, 3)))
write.table(peak_annotations, file.path(DIR_CAQTL, "peakLoc", paste0(ct, ".peaks.hg38.bed")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# Sample covariate table (sex coded 0/1; + SVA factors).
sample_metadata <- anno.samp[, c("Sample2", "Age", "Sex", "PMI.x", "BatchLib",
                                 "Sva_factor1", "Sva_factor2", "Sva_factor3")]
sample_metadata <- unique(sample_metadata)
sample_metadata$Sex <- ifelse(sample_metadata$Sex == "F", 0, 1)
sample_metadata <- t(sample_metadata[match(colnames(y.mtx), sample_metadata$Sample2), ])
write.table(sample_metadata, file.path(DIR_CAQTL, "covariates", paste0(ct, ".sample.metadata.txt")),
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

## ---- SNP genotype matrices (dosage encoding) --------------------------------

inDir <- GENOTYPE_DIR
outDir <- file.path(DIR_CAQTL, "snpMatrix")
for (i in 1:22) {
  vcf <- fread(file.path(inDir, paste0("GWAS.chr", i, ".allele_align.all_imputed_20240924.fam.PD5D_include_midBrain.vcf")))
  colnames(vcf)[1:2] <- c("chr", "pos")
  snpLoc <- vcf[, 1:3]
  snpLoc$chr <- paste0("chr", snpLoc$chr)
  snpLoc <- data.frame(snp = paste("snp", snpLoc$chr, snpLoc$pos, snpLoc$ID, sep = "_"), snpLoc)
  write.table(snpLoc[, 1:3], file.path(DIR_CAQTL, "snpLoc", paste0("chr", i, "-snp-location-renew.txt")),
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
# Run inside a dedicated PEER conda environment (see shell block below). PEER is
# loaded here; the matrix is samples x peaks (top-variable peaks selected).

library(peer)
library(matrixStats)

peer_input_file <- file.path(DIR_CAQTL, "peakMatrix", paste0(ct, ".TMM.voom.normed.mtx"))
if (!file.exists(peer_input_file)) stop(paste("Transformed matrix file not found for cell type", ct))
peer_mtx <- read.table(peer_input_file, header = TRUE, sep = "\t", row.names = 1)
peer_mtx <- t(as.matrix(peer_mtx))

model <- PEER()
n_samples <- nrow(peer_mtx)
n_hidden <- min(20, max(4, floor(n_samples / 4)))   # K conservative vs sample size
vars <- colVars(peer_mtx)
target_p <- min(50000, ncol(peer_mtx))
sel <- order(vars, decreasing = TRUE)[seq_len(target_p)]
peer_mtx <- peer_mtx[, sel, drop = FALSE]
storage.mode(peer_mtx) <- "double"

PEER_setNk(model, n_hidden)
PEER_setPhenoMean(model, peer_mtx)
PEER_setNmax_iterations(model, 1000)

covariates_file <- file.path(DIR_CAQTL, "covariates", paste0(ct, ".sample.metadata.txt"))
if (file.exists(covariates_file)) {
  sample_covariates <- read.table(covariates_file, header = TRUE, sep = "\t", row.names = 1)
  sample_covariates <- t(as.matrix(sample_covariates))
  PEER_setCovariates(model, sample_covariates)
}

PEER_update(model)
peer_factors <- PEER_getX(model)
rownames(peer_factors) <- rownames(peer_mtx)
colnames(peer_factors) <- paste0("InferredCov", 1:ncol(peer_factors))
peer_residuals <- t(PEER_getResiduals(model))
rownames(peer_residuals) <- colnames(peer_mtx)
colnames(peer_residuals) <- rownames(peer_mtx)

write.table(peer_factors, file.path(DIR_CAQTL, "PEER", paste0(ct, ".PEER_covariates.txt")),
            row.names = TRUE, sep = "\t", quote = FALSE)
write.table(peer_residuals, file.path(DIR_CAQTL, "PEER", paste0(ct, ".PEER_residuals.txt")),
            row.names = TRUE, sep = "\t", quote = FALSE)

## ---- Assemble FastQTL BED / VCF / covariate inputs --------------------------

peak_loc_dir <- file.path(DIR_CAQTL, "peakLoc")
out_dir_bed  <- file.path(DIR_CAQTL, "Fastqtl/Peakbed")
out_dir_snp  <- file.path(DIR_CAQTL, "Fastqtl/SNPvcf")
inDir <- GENOTYPE_DIR

ct_file <- file.path(DIR_CAQTL, "peakMatrix", paste0(ct, ".TMM.voom.normed.mtx"))
peak_matrix <- fread(ct_file)
peak_loc <- fread(file.path(peak_loc_dir, paste0(ct, ".peaks.hg38.bed")), header = TRUE)
bed <- peak_loc %>% left_join(peak_matrix, by = c("peakid" = "peakid")) %>%
  select(chr, start, end, peakid, everything()) %>% arrange(chr, start)
colnames(bed)[1] <- "chr"

calculate_maf <- function(genotypes) {
  alleles <- unlist(strsplit(as.character(genotypes), split = "[/]"))
  allele_counts <- table(alleles)
  freqs <- allele_counts / sum(allele_counts)
  min(freqs)
}

for (i in 1:22) {
  vcf <- fread(file.path(inDir, paste0("GWAS.chr", i, ".allele_align.all_imputed_20240924.fam.PD5D_include_midBrain.vcf")))
  vcf1 <- vcf
  colnames(vcf1)[10:ncol(vcf1)] <- substr(colnames(vcf1)[10:ncol(vcf1)], 1, 6)
  common_samples <- intersect(colnames(bed)[5:ncol(bed)], colnames(vcf1))
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
write.table(bed, file.path(out_dir_bed, paste0(ct, "_peakmatrix.hg38.bed")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
for (chr in unique(bed$chr)) {
  bed_chr <- bed %>% filter(chr == !!chr)
  write.table(bed_chr, file.path(out_dir_bed, paste0(ct, "_peakmatrix.hg38_", chr, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Combine sample covariates + top-20 PEER factors into the FastQTL covariate file.
peer_file <- file.path(DIR_CAQTL, "PEER", paste0(ct, ".PEER_covariates.txt"))
peer_factors <- t(read.table(peer_file, header = TRUE, sep = "\t", row.names = 1))
covariates_file <- file.path(DIR_EQTL, "covariates", paste0(ct, ".sample.metadata.txt"))
sample_covariates <- read.table(covariates_file, header = TRUE, sep = "\t", row.names = 1)
allcovaraints <- rbind(sample_covariates, peer_factors[1:20, ])
allcovaraints <- allcovaraints[, common_samples]
allcovaraints <- data.frame(ID = rownames(allcovaraints), allcovaraints)
write.table(allcovaraints, file.path(DIR_CAQTL, "Fastqtl/Covariates", paste0(ct, ".Covariates.txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)

# =============================================================================
# Shell / conda blocks (run separately, not as R)
# =============================================================================
#
# ---- 1. Create the PEER conda environment ----
# conda create --name peer python=3.5.4
# conda activate peer
# conda config --add channels conda-forge
# conda config --add channels r
# conda install zlib=1.2.8 r=3.4.1 r-peer libgfortran=3 conda-forge::r-qtl
#
# ---- 2. Sort, bgzip, and tabix-index VCF and BED files ----
# cell_types=("Neurons")
# cd 03-caQTL-calling2/Fastqtl/SNPvcf/
# for ct in "${cell_types[@]}"; do
#   for vcf_file in ${ct}*_snp_matrix.vcf; do
#     sorted_vcf=${vcf_file%.vcf}.sort.vcf
#     grep '^#' $vcf_file > $sorted_vcf
#     grep -v '^#' $vcf_file | LC_ALL=C sort -k1,1 -k2,2n >> $sorted_vcf
#     bgzip $sorted_vcf; tabix -p vcf ${sorted_vcf}.gz
#   done
# done
# cd 03-caQTL-calling2/Fastqtl/Peakbed/
# for ct in "${cell_types[@]}"; do
#   for bed_file in ${ct}*.bed; do
#     sorted_bed=${bed_file%.bed}.sort.bed
#     { head -n 1 $bed_file | sed 's/^/#/'; tail -n +2 $bed_file | LC_ALL=C sort -k1,1 -k2,2n; } > $sorted_bed
#     bgzip $sorted_bed; tabix -p bed ${sorted_bed}.gz
#   done
# done
#
# ---- 3. Build FastQTL command files (permutation + nominal, per chromosome) ----
# Permutation : fastQTL.static --vcf VCF --bed BED --cov COV --out PERM --window 1e6 --permute 1000 --region chrN
# Nominal     : fastQTL.static --vcf VCF --bed BED --cov COV --out NOM  --window 1e6 --region chrN
# Then run e.g.: ParaFly -c Neurons_fastqtl_commands.txt -CPU 11
#
# ---- 4. Combine per-chromosome results and count significant caQTLs ----
# Permutation header:
#   gene_id num_var beta_shape1 beta_shape2 true_df pval_true_df variant_id
#   tss_distance ma_samples ma_count maf ref_factor pval_nominal slope slope_se
#   pval_perm pval_beta
# Significant caQTLs: awk '$16 < 0.05' ${ct}_combined_permutations.txt | wc -l
# =============================================================================
