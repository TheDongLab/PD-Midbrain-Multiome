#!/bin/bash
# =============================================================================
# Step 08 - Partitioned heritability of PD GWAS by cell-type QTL annotations
#           (stratified LD score regression, S-LDSC)
# -----------------------------------------------------------------------------
# Purpose : Build cell-type QTL annotations, compute LD scores against the
#           1000G EUR reference, and run S-LDSC to estimate PD-heritability
#           enrichment for each cell type's caQTL (aQTL) and eQTL SNPs.
# Inputs  : Combined FastQTL significant QTLs (steps 05/06), PD GWAS sumstats,
#           LDSC reference panels (baseline v1.2, weights, frq, 1000G plink).
# Outputs : <BASE>/qtl-LDSC/output/{aqtl,eqtl}.PD_GWAS_partitioned_<ct>.results
# Tools   : LDSC (make_annot.py, ldsc.py, munge_sumstats.py), liftOver.
# Usage   : bash GWAS-LDSC-QTL.sh   (run stages as needed; long-running)
# Notes   : Main steps only. The peak-based annotation sets (all / promoter /
#           non-promoter / aQTL-peaks / cell-type-specific peaks) follow the
#           identical Step 1-3 procedure with different input BEDs.
# =============================================================================
set -euo pipefail

# ---- Settings (edit for your environment) ----
ROOT="${ASAP_MULTIOME_ROOT:-/mnt/data/projects/donglab/ming2024/ASAP_Multiome}"
LDSC_HOME="${LDSC_HOME:-/home/ml3277/ldsc}"
BASE="${ROOT}/04-PDGWAS"
REF="${BASE}/S-LDSC"                 # LDSC reference panels
QTL="${BASE}/qtl-LDSC"               # annotations, LD scores, and outputs
PLINK="${REF}/1000G_EUR_Phase3_plink/1000G.EUR.QC"
FRQ="${REF}/1000G_Phase3_frq/1000G.EUR.QC."
WLD="${REF}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
BASELINE="${REF}/baseline_v1.2/baseline."
SUMSTATS="${QTL}/PD_GWAS.sumstats.gz"

cell_types=("DAneurons" "GABAneurons" "Gluneurons" "Oligodendrocytes" \
            "Microglia" "OPCs" "Astrocytes" "Endothelialcell" "Pericytes" \
            "Macrophages" "TCells")

mkdir -p "${QTL}/output"

# =============================================================================
# Step 0 (one-time): references, munge sumstats, HapMap3 SNPs, liftOver chain
# -----------------------------------------------------------------------------
# Install LDSC once: git clone https://github.com/bulik/ldsc.git
#   cd ldsc && conda env create --file environment.yml && conda activate ldsc
#
# Download reference panels into ${REF} (baseline v1.2 LD scores, HM3 weights,
# 1000G frq, 1000G plink) from the LDSC Zenodo records, e.g.:
#   wget .../1000G_Phase3_baseline_v1.2_ldscores.tgz && tar -xvzf ...
#   wget .../1000G_Phase3_weights_hm3_no_MHC.tgz     && tar -xvzf ...
#   wget .../1000G_Phase3_frq.tgz                    && tar -xvzf ...
#   wget .../1000G_Phase3_plinkfiles.tgz             && tar -xvzf ...
# liftOver chain + binary:
#   wget .../hg38ToHg19.over.chain.gz && gunzip hg38ToHg19.over.chain.gz
#
# Munge the PD GWAS summary statistics (hg19) to LDSC format:
python "${LDSC_HOME}/munge_sumstats.py" \
  --sumstats "${REF}/PD_GWAS_with_rsIDs_hg19.tab" \
  --out "${QTL}/PD_GWAS" \
  --merge-alleles "${REF}/w_hm3.snplist" \
  --chunksize 500000

# HapMap3 SNPs per chromosome (used with --print-snps).
mkdir -p "${REF}/hapmap3_snps"
for chr in {1..22}; do
  awk 'NR==FNR {hm3[$1]; next} ($2 in hm3) {print $2}' \
    "${REF}/w_hm3.snplist" "${PLINK}.${chr}.bim" > "${REF}/hapmap3_snps/hm.${chr}.snp"
done

# -----------------------------------------------------------------------------
# Step 1 (R): build cell-type significant-QTL annotation BEDs (+/- 5 kb windows)
# -----------------------------------------------------------------------------
# Run once in R to produce ${QTL}/<ct>.filtered_significant_{aQTLs,eQTLs}.bed
# from the combined FastQTL significant-QTL tables:
#
#   library(data.table)
#   setwd("${ROOT}")
#   snp <- fread("./03-caQTL-calling/snpLoc/all-snp-location-renew.txt")
#   snp[, SNP := sub(".*_(rs[0-9]+)$", "\\1", snp)]
#   for (ct in cell_types) {
#     for (qtl in c("aQTL", "eQTL")) {
#       path <- if (qtl == "aQTL") "./03-caQTL-calling2/Fastqtl/Output_Combined"
#               else               "./05-eQTL-calling2/Fastqtl/Output_Combined"
#       d <- fread(file.path(path, paste0(ct, ".filtered_significant_", qtl, "s.csv")))
#       d <- merge(d, snp, by.x = "variant_id", by.y = "SNP")
#       d$start <- as.integer(d$pos) - 5000; d$end <- d$start + 5000
#       fwrite(d[, .(chr, start, end, gene_id, variant_id)],
#              paste0("./04-PDGWAS/qtl-LDSC/", ct, ".filtered_significant_", qtl, "s.bed"),
#              sep = "\t", col.names = FALSE)
#     }
#   }

# =============================================================================
# Steps 2-3: LD scores + partitioned heritability per QTL annotation
# -----------------------------------------------------------------------------
# run_sldsc <prefix> <bed_suffix>
#   prefix     : annotation name used in output files (e.g. aqtl, eqtl)
#   bed_suffix : input BED suffix (e.g. filtered_significant_aQTLs)
run_sldsc() {
  local prefix="$1" bed_suffix="$2"

  for ANNOT in "${cell_types[@]}"; do
    echo "[${prefix}] Processing cell type: ${ANNOT}"

    # lift the QTL BED hg38 -> hg19 to match the 1000G reference.
    "${QTL}/liftOver" \
      "${QTL}/${ANNOT}.${bed_suffix}.bed" \
      "${QTL}/hg38ToHg19.over.chain" \
      "${QTL}/${ANNOT}.${bed_suffix}.hg19.bed" \
      "${QTL}/${ANNOT}_${prefix}_unmapped.bed"

    # Step 2: per-chromosome annotation + LD scores.
    for chrom in {1..22}; do
      python "${LDSC_HOME}/make_annot.py" \
        --bed-file "${QTL}/${ANNOT}.${bed_suffix}.hg19.bed" \
        --bimfile "${PLINK}.${chrom}.bim" \
        --annot-file "${QTL}/${prefix}.${ANNOT}.${chrom}.annot.gz"

      python "${LDSC_HOME}/ldsc.py" \
        --l2 \
        --bfile "${PLINK}.${chrom}" \
        --print-snps "${REF}/hapmap3_snps/hm.${chrom}.snp" \
        --ld-wind-cm 1 \
        --annot "${QTL}/${prefix}.${ANNOT}.${chrom}.annot.gz" \
        --thin-annot \
        --out "${QTL}/${prefix}.${ANNOT}.${chrom}"
    done

    # Step 3: partitioned heritability (baseline + this cell type's annotation).
    python "${LDSC_HOME}/ldsc.py" \
      --h2 "${SUMSTATS}" \
      --ref-ld-chr "${BASELINE},${QTL}/${prefix}.${ANNOT}." \
      --frqfile-chr "${FRQ}" \
      --w-ld-chr "${WLD}" \
      --overlap-annot --print-cov --print-coefficients --print-delete-vals \
      --out "${QTL}/output/${prefix}.PD_GWAS_partitioned_${ANNOT}" || {
        echo "LDSC failed for ${prefix} ${ANNOT}; skipping."; continue; }
  done
}

# Chromatin-accessibility QTLs and expression QTLs.
run_sldsc aqtl filtered_significant_aQTLs
run_sldsc eqtl filtered_significant_eQTLs

echo "S-LDSC partitioned heritability complete. Results in ${QTL}/output/"
