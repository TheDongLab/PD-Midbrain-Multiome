# =============================================================================
# Step 07 - Multi-trait colocalization (moloc) of caQTL + eQTL + PD GWAS
# -----------------------------------------------------------------------------
# Purpose : For each cell type, run moloc three-trait colocalization across PD
#           GWAS, eQTL and caQTL signals around genome-wide-significant lead SNPs.
# Inputs  : GWAS_SUMSTATS, clumped GWAS lead SNPs, SNP locations,
#           combined caQTL/eQTL nominal results (<DIR_CAQTL>/<DIR_EQTL>).
# Outputs : <DIR_MOLOC>/<celltype>_moloc_result-clumped-4.txt
# Run     : Rscript moloc-coloc.R   (example cell type: pooled "Neurons")
# Notes   : Core colocalization only. Downstream LocusZoom inputs, stacked locus
#           plots, and motif analyses (motifmatchr/ggseqlogo/motifbreakR) were
#           removed. Interactive overrides (in-loop cell_type / i=1) were removed.
# =============================================================================

source("config.R")
options(scipen = 999)

suppressPackageStartupMessages({
  library(moloc)
  library(data.table)
  library(GenomicRanges)
  library(dplyr)
})

CAQTL_NOMINALS <- file.path(DIR_CAQTL, "Fastqtl/Output_Combined")
EQTL_NOMINALS  <- file.path(DIR_EQTL, "Fastqtl/Output_Combined")
SNP_LOC        <- file.path("./03-caQTL-calling/snpLoc/all-snp-location-renew.txt")
GWAS_CLUMPED   <- file.path(DIR_GWAS, "PD_GWAS_with_rsIDs_hg19_clumpped2.jma.cojo")

cell_types <- c("Neurons")

## ---- Load and clean PD GWAS summary statistics ------------------------------

gwas_data <- read.table(GWAS_SUMSTATS, header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
colnames(gwas_data) <- c("name", "A1", "A2", "freq", "BETA", "SE", "PVAL",
                         "Ncases", "Ncontrols", "rsID", "CHR", "Start", "POS", "SNP")
gwas_data$other_allele_frequency <- 1 - as.numeric(as.character(gwas_data$freq))
gwas_data$MAF <- as.numeric(apply(data.frame(gwas_data$freq, gwas_data$other_allele_frequency), 1, min))
gwas_data$N <- gwas_data$Ncases + gwas_data$Ncontrols
gwas_data$CHR <- gsub("chr", "", gwas_data$CHR)
gwas_data <- gwas_data[, c("SNP", "CHR", "POS", "BETA", "SE", "PVAL", "A1", "A2", "MAF", "N", "Ncases")]
gwas_data2 <- gwas_data[order(gwas_data$PVAL), ]
gwas_data2 <- gwas_data2[!duplicated(gwas_data2$SNP), ]
gwas_data2 <- na.omit(gwas_data2)
gwas_data2$PVAL <- format(gwas_data2$PVAL, scientific = TRUE, digits = 2)

# Genome-wide significant, LD-clumped lead SNPs.
significant_snps <- read.table(GWAS_CLUMPED, header = FALSE, sep = "\t", skip = 1)
gwas <- merge(significant_snps, gwas_data2, by.x = "V2", by.y = "SNP", all.x = TRUE)
gwas <- gwas[gwas$PVAL < 5e-8, ]

snp_data <- fread(SNP_LOC, header = TRUE)
snp_data[, SNP := sub(".*_(rs[0-9]+)$", "\\1", snp)]

## ---- moloc three-trait colocalization (GWAS x eQTL x caQTL) -----------------

for (cell_type in cell_types) {
  print(cell_type)
  caqtl_data <- fread(file.path(CAQTL_NOMINALS, paste0(cell_type, "_combined_nominals.txt")), header = TRUE, sep = "\t")
  caqtl_data2 <- merge(caqtl_data, snp_data, by.x = "variant_id", by.y = "SNP")
  caqtl_data2 <- caqtl_data2[, c("gene_id", "variant_id", "chr", "pos", "slope", "slope_se", "pval_nominal", "maf")]
  colnames(caqtl_data2) <- c("ProbeID", "SNP", "CHR", "POS", "BETA", "SE", "PVAL", "MAF")
  caqtl_data2$N <- 85
  caqtl_data2$POS <- as.integer(caqtl_data2$POS)
  caqtl_data2$CHR <- as.integer(gsub("chr", "", caqtl_data2$CHR))
  caqtl_data2 <- na.omit(caqtl_data2)

  eqtl_data <- fread(file.path(EQTL_NOMINALS, paste0(cell_type, "_combined_nominals.txt")), header = TRUE, sep = "\t")
  eqtl_data2 <- merge(eqtl_data, snp_data, by.x = "variant_id", by.y = "SNP")
  eqtl_data2 <- eqtl_data2[, c("gene_id", "variant_id", "chr", "pos", "slope", "slope_se", "pval_nominal", "maf")]
  colnames(eqtl_data2) <- c("ProbeID", "SNP", "CHR", "POS", "BETA", "SE", "PVAL", "MAF")
  eqtl_data2$N <- 85
  eqtl_data2$POS <- as.integer(eqtl_data2$POS)
  eqtl_data2$CHR <- as.integer(gsub("chr", "", eqtl_data2$CHR))
  eqtl_data2 <- na.omit(eqtl_data2)

  out_file <- file.path(DIR_MOLOC, paste0(cell_type, "_moloc_result-clumped-4.txt"))
  for (i in 1:nrow(gwas)) {
    snp <- gwas$V2[i]
    win <- 100000
    input1 <- gwas_data2[gwas_data2$CHR == gwas$CHR[i] & gwas_data2$POS <= gwas$POS[i] + win & gwas_data2$POS >= gwas$POS[i] - win, ]
    input2 <- eqtl_data2[eqtl_data2$CHR == gwas$CHR[i] & eqtl_data2$POS <= gwas$POS[i] + win & eqtl_data2$POS >= gwas$POS[i] - win, ]
    input3 <- caqtl_data2[caqtl_data2$CHR == gwas$CHR[i] & caqtl_data2$POS <= gwas$POS[i] + win & caqtl_data2$POS >= gwas$POS[i] - win, ]

    if (nrow(input1) > 0 & nrow(input2) > 0 & nrow(input3) > 0) {
      print(i); print(snp)
      genes <- as.character(unique(input2$ProbeID))
      lncs <- as.character(unique(input3$ProbeID))
      for (gene in genes) {
        for (lnc in lncs) {
          input_gene <- input2[input2$ProbeID == gene, ]
          input_lnc <- input3[input3$ProbeID == lnc, ]
          common <- intersect(input1$SNP, intersect(input_gene$SNP, input_lnc$SNP))
          if (length(common) > 10) {
            l1 <- input1[input1$SNP %in% common, ]; l1 <- l1[!duplicated(l1$SNP), ]
            l2 <- input_gene[input_gene$SNP %in% common, ]
            l3 <- input_lnc[input_lnc$SNP %in% common, ]; l3 <- l3[!duplicated(l3$SNP), ]
            input <- list(l1, l2, l3)
            options(scipen = 1, digits = 2)
            moloc <- moloc_test(input, prior_var = c(0.01, 0.1, 0.5), priors = c(1e-04, 1e-06, 1e-07))
            cat(snp, "\t", gwas$CHR[i], "\t", gwas$POS[i], "\t", gene, "\t", lnc, "\t", input_gene$CHR[1], "\t", nrow(l1),
                "\t", moloc[[1]]["abc", "PPA"], "\t", moloc[[1]]["ab,c", "PPA"], "\t", moloc[[1]]["a,bc", "PPA"], "\t", moloc[[1]]["ac,b", "PPA"],
                "\t", moloc[[3]]["abc", "best.snp.coloc"], "\t", moloc[[3]]["ab", "best.snp.coloc"], "\t", moloc[[3]]["bc", "best.snp.coloc"], "\t", moloc[[3]]["ac", "best.snp.coloc"],
                "\n", file = out_file, append = TRUE)
          }
        }
      }
    }
  }
}
