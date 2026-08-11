# =============================================================================
# config.R - Central configuration for the ASAP Multiome (PD substantia nigra)
#            single-nucleus multiome analysis pipeline.
#
# All cleaned analysis scripts begin with `source("config.R")`. This is the only
# place that contains machine-specific absolute paths and global parameters, so
# the pipeline can be relocated by editing this file alone (or by setting the
# environment variables referenced below).
#
# The scripts in this folder are trimmed to core computational steps only
# (plotting, QC summaries, and downstream interpretation blocks were removed
# from most steps). Step-local run targets (e.g. `ct <- "Neurons"`, `CELL_TYPE`)
# remain at the top of individual scripts.
#
# Usage:
#   - Run scripts from this directory, or set ASAP_MULTIOME_ROOT before sourcing.
#   - Edit the paths / tool locations in the "USER SETTINGS" block as needed.
#   - Shell scripts (08, 11) read ASAP_MULTIOME_ROOT / ASAP_CELLRANGER_DIR /
#     LDSC_HOME directly from the environment (defaults match this file).
# =============================================================================

## ---- USER SETTINGS (edit for your environment) -----------------------------

# Project root: all relative paths below are resolved against this directory.
PROJECT_ROOT <- Sys.getenv(
  "ASAP_MULTIOME_ROOT",
  unset = "/mnt/data/projects/donglab/ming2024/ASAP_Multiome"
)

# Cell Ranger ARC count output (per-sample h5 / fragments / barcode metrics).
CELLRANGER_DIR <- Sys.getenv(
  "ASAP_CELLRANGER_DIR",
  unset = "/mnt/data/projects/donglab/multiome/2023_summer_cellranger-arc_count"
)

# Path to the MACS (>=2) / MACS3 executable used for ATAC peak calling.
MACS_PATH <- Sys.getenv(
  "ASAP_MACS_PATH",
  unset = "/home/ml3277/miniconda3/envs/single_cell_env/bin/macs3"
)

# LDSC installation (used by GWAS-LDSC-QTL.sh).
LDSC_HOME <- Sys.getenv(
  "LDSC_HOME",
  unset = "/home/ml3277/ldsc"
)

# Donor / sample metadata spreadsheet (used by Filtering-Clustering.R).
METADATA_XLSX <- Sys.getenv(
  "ASAP_METADATA_XLSX",
  unset = "/mnt/data/projects/donglab/ming2024/ASAP_scATAC_2024/metadata/2023-07-26-midbrain-scATAC-metadata.xlsx"
)

## ---- DERIVED PATHS (usually no need to edit) --------------------------------

if (!dir.exists(PROJECT_ROOT)) {
  warning("PROJECT_ROOT does not exist: ", PROJECT_ROOT,
          "\nSet the ASAP_MULTIOME_ROOT environment variable or edit config.R.")
}
setwd(PROJECT_ROOT)

# Core single-cell objects / intermediate output directories.
DIR_FINALWNN   <- "./01-FinalWNN"                  # merged + integrated Seurat objects
DIR_MACS       <- "./MACS3_output3+Tcells"         # MACS peaks + counts
DIR_CAQTL      <- "./03-caQTL-calling2"            # caQTL inputs / FastQTL output
DIR_EQTL       <- "./05-eQTL-calling2"             # eQTL inputs / FastQTL output
DIR_GWAS       <- "./04-PDGWAS"                     # PD GWAS, S-LDSC, QTL-LDSC
DIR_MOLOC      <- "./06-moloc-new2"                 # moloc colocalization results
DIR_DIFF       <- "./07-PD-differential"           # DEG / DAR / scMultiMap / SCENIC+
DIR_SCENIC     <- "./08-scenicplus_celltype"       # SCENIC+ eRegulon exports per cell type
DIR_DATA       <- "./00-data"                      # sample covariates, SVA factors, etc.

# Canonical labeled Seurat object and MACS peak count matrix.
SEURAT_OBJECT  <- file.path(DIR_FINALWNN, "merged_seurat4+Tcells.rda")
MACS_COUNTS    <- file.path(DIR_MACS, "peaks_macs_count.rda")

# Sample covariates + SVA factors (steps 05/06).
COVARIATES_FILE <- file.path(DIR_DATA, "covariants_100_samples_with_first_three_sva_factors.txt")

# Imputed genotype VCFs (GWAS-aligned), one per chromosome (steps 05/06).
GENOTYPE_DIR <- file.path(
  "./03-caQTL-calling",
  "genotype_reimputed_allele_aligned_to_GWAS_PD5D_eQTL"
)

# GENCODE annotation for eQTL gene locations (step 06).
GENE_GFF3 <- file.path(DIR_EQTL, "gencode.v47.annotation.gff3")

# SNP location table shared by QTL and moloc steps.
SNP_LOC <- file.path("./03-caQTL-calling/snpLoc/all-snp-location-renew.txt")

# PD GWAS summary statistics and clumped lead SNPs (steps 07/08).
GWAS_SUMSTATS  <- file.path(DIR_GWAS, "S-LDSC/PD_GWAS_with_rsIDs_hg38.tab")
GWAS_CLUMPED   <- file.path(DIR_GWAS, "PD_GWAS_with_rsIDs_hg19_clumpped2.jma.cojo")

# Combined FastQTL nominal QTL results (step 07).
CAQTL_NOMINALS <- file.path(DIR_CAQTL, "Fastqtl/Output_Combined")
EQTL_NOMINALS  <- file.path(DIR_EQTL, "Fastqtl/Output_Combined")

## ---- GLOBAL PARAMETERS / CONSTANTS ------------------------------------------

GENOME_BUILD <- "hg38"

# Metadata column holding the finalized WNN cell-type labels.
CELLTYPE_COLUMN <- "celltype.WNN.after.combine.update"

# Cell types used across QTL / differential / network analyses.
CELL_TYPES <- c(
  "GABA neurons", "Glu neurons", "DA neurons", "Microglia", "Astrocytes",
  "Oligodendrocytes", "OPCs", "T Cells", "Endothelial cell", "Macrophages",
  "Pericytes"
)

# Diagnosis contrasts (case, reference) for differential analysis.
COMPARISONS <- list(
  PD_vs_HC   = c("PD", "HC"),
  ILBD_vs_HC = c("ILBD", "HC"),
  PD_vs_ILBD = c("PD", "ILBD")
)
