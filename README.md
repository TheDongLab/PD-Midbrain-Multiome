# Cleaned analysis scripts (manuscript / reuse)

Reorganized versions of the main step scripts for the ASAP Multiome (PD substantia
nigra) single-nucleus multiome pipeline. Each script keeps the original analysis
logic for its core computational steps, adds a documented header, sources a shared
`config.R`, consolidates `library()` calls, and removes package-install,
interactive (`tmux`/debug), and non-English clutter.

**Scope:** Most scripts were further trimmed to **core steps only** — QC plots,
summary tables, and downstream interpretation blocks (GO/GSEA, LocusZoom, motif
analysis, cross-cell-type summaries, etc.) were removed so the folder is suitable
as a manuscript supplement. `Filtering-Clustering.R` still retains QC and UMAP
figures because clustering and manual annotation depend on them. S-LDSC plotting
remains in the original `../scripts/GWAS-LDSC-QTL-plot+Neurons.R` (not copied
here).

## How to run

1. Edit `config.R` and set `PROJECT_ROOT` (and, if needed, `CELLRANGER_DIR`,
   `MACS_PATH`, `LDSC_HOME`, `METADATA_XLSX`) for your environment, or export the
   matching environment variables (`ASAP_MULTIOME_ROOT`, `ASAP_CELLRANGER_DIR`,
   `ASAP_MACS_PATH`, `LDSC_HOME`, `ASAP_METADATA_XLSX`).
2. Run an R script from this directory, e.g. `Rscript Filtering-Clustering.R`.
   Every R script begins with `source("config.R")`, which sets the working
   directory to `PROJECT_ROOT`; all other paths are relative to it.
3. Shell scripts (`GWAS-LDSC-QTL.sh`, `bw-files-for-ATAC.sh`) read the same
   environment variables and can be run with `bash <script>.sh`.

`config.R` centralizes output directories, canonical object paths (`SEURAT_OBJECT`,
`MACS_COUNTS`), QTL/GWAS inputs (`GENOTYPE_DIR`, `COVARIATES_FILE`, `SNP_LOC`,
`GWAS_SUMSTATS`, `GWAS_CLUMPED`), the cell-type list (`CELL_TYPES`), and diagnosis
contrasts (`COMPARISONS`). Step-local run targets (e.g. `ct <- "Neurons"` in QTL
preprocessing, `CELL_TYPE <- "Microglia"` in eRegulon GSEA) remain at the top of
the relevant script — change them or wrap the script in a loop to process all cell
types.

## Scripts (one per main pipeline step)

| Script | Step | Summary |
|--------|------|---------|
| `config.R` | – | Shared paths, tool locations, and global parameters |
| `Combining-Samples.R` | 01 | Build per-sample multiome Seurat objects and merge |
| `Filtering-Clustering.R` | 02 | QC filtering, SCT/LSI, Harmony, WNN, doublets, annotation *(QC/UMAP plots retained)* |
| `Callpeaks-CellTypes.R` | 03 | MACS peaks, union peak set, counts, cell-type-specific peaks |
| `caQTL-preprocess.R` | 05 | Pseudobulk accessibility, PEER, FastQTL inputs (+ shell blocks) |
| `eQTL-preprocess.R` | 06 | Pseudobulk expression, PEER, FastQTL inputs (+ shell blocks) |
| `moloc-coloc.R` | 07 | moloc caQTL + eQTL + GWAS three-trait colocalization |
| `GWAS-LDSC-QTL.sh` | 08 | S-LDSC partitioned heritability (annotations, LD scores, h²) |
| `PD-differential-nebula-deg.R` | 09 | Single-cell NEBULA DEG across diagnosis contrasts |
| `scMultiMap-peak-gene.R` | 10 | scMultiMap peak–gene cis linking |
| `bw-files-for-ATAC.sh` | 11 | Per-cell-type normalized ATAC bigWig tracks |
| `scenicplus_step1_per_celltype.py` | 12 | SCENIC+ pycisTopic (LDA / topics / DARs) per cell type |
| `scenicplus_step2_per_celltype.py` | 12 | cisTarget motif database construction |
| `scenicplus_step3_per_celltype.py` | 12 | SCENIC+ Snakemake config YAML generation |
| `eRegulon-GSEA-core.R` | 13 | eRegulon gene-set GSEA on NEBULA DEG rankings (per cell type) |

