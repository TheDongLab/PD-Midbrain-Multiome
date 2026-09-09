# Configuration template for the PD Midbrain Multiome analysis.
# Replace the placeholder paths or set the corresponding environment variables.

## User settings

PROJECT_ROOT <- Sys.getenv(
  "ASAP_MULTIOME_ROOT",
  unset = "/path/to/ASAP_Multiome"
)

CELLRANGER_DIR <- Sys.getenv(
  "ASAP_CELLRANGER_DIR",
  unset = "/path/to/cellranger_arc_outputs"
)

MACS_PATH <- Sys.getenv(
  "ASAP_MACS_PATH",
  unset = "/path/to/bin/macs3"
)

LDSC_HOME <- Sys.getenv(
  "LDSC_HOME",
  unset = "/path/to/ldsc"
)

METADATA_XLSX <- Sys.getenv(
  "ASAP_METADATA_XLSX",
  unset = "/path/to/midbrain_metadata.xlsx"
)

if (!dir.exists(PROJECT_ROOT)) {
  stop(
    "PROJECT_ROOT does not exist: ", PROJECT_ROOT,
    "\nSet ASAP_MULTIOME_ROOT or edit config.R before running the analysis."
  )
}

setwd(PROJECT_ROOT)

## Analysis directories

DIR_FINALWNN <- "./01-FinalWNN"
DIR_MACS     <- "./MACS3_output3+Tcells"
DIR_CAQTL    <- "./03-caQTL-calling2"
DIR_EQTL     <- "./05-eQTL-calling2"
DIR_GWAS     <- "./04-PDGWAS"
DIR_MOLOC    <- "./06-moloc-new2"
DIR_DIFF     <- "./07-PD-differential"
DIR_SCENIC   <- "./08-scenicplus_celltype"
DIR_DATA     <- "./00-data"

## Input and intermediate files

SEURAT_OBJECT <- file.path(DIR_FINALWNN, "merged_seurat4+Tcells.rda")
MACS_COUNTS   <- file.path(DIR_MACS, "peaks_macs_count.rda")

COVARIATES_FILE <- file.path(
  DIR_DATA,
  "covariants_100_samples_with_first_three_sva_factors.txt"
)

GENOTYPE_DIR <- file.path(
  "./03-caQTL-calling",
  "genotype_reimputed_allele_aligned_to_GWAS_PD5D_eQTL"
)

GENE_GFF3 <- file.path(DIR_EQTL, "gencode.v47.annotation.gff3")
SNP_LOC <- "./03-caQTL-calling/snpLoc/all-snp-location-renew.txt"

GWAS_SUMSTATS <- file.path(
  DIR_GWAS,
  "S-LDSC/PD_GWAS_with_rsIDs_hg38.tab"
)
GWAS_CLUMPED <- file.path(
  DIR_GWAS,
  "PD_GWAS_with_rsIDs_hg19_clumpped2.jma.cojo"
)

CAQTL_NOMINALS <- file.path(DIR_CAQTL, "Fastqtl/Output_Combined")
EQTL_NOMINALS  <- file.path(DIR_EQTL, "Fastqtl/Output_Combined")

## Global parameters

GENOME_BUILD <- "hg38"
CELLTYPE_COLUMN <- "celltype.WNN.after.combine.update"

CELL_TYPES <- c(
  "GABA neurons", "Glu neurons", "DA neurons", "Microglia", "Astrocytes",
  "Oligodendrocytes", "OPCs", "T Cells", "Endothelial cell", "Macrophages",
  "Pericytes"
)

COMPARISONS <- list(
  PD_vs_HC   = c("PD", "HC"),
  ILBD_vs_HC = c("ILBD", "HC"),
  PD_vs_ILBD = c("PD", "ILBD")
)
