# Cleaned analysis scripts (manuscript / reuse)

Reorganized versions of the main step scripts for the ASAP Multiome (PD substantia
nigra) single-nucleus multiome pipeline. Each script keeps the original analysis
logic for its core computational steps, adds a documented header, sources a shared
`config.R`, consolidates `library()` calls, and removes package-install,
interactive (`tmux`/debug), and non-English clutter.

**Scope:** Most scripts were further trimmed to **core steps only** — QC plots,
summary tables, and downstream interpretation blocks (GO/GSEA, LocusZoom, motif
analysis, cross-cell-type summaries, etc.) were removed so the folder is suitable
as a manuscript supplement. `02-Filtering-Clustering.R` still retains QC and UMAP
figures because clustering and manual annotation depend on them. S-LDSC plotting
remains in the original `../scripts/08-GWAS-LDSC-QTL-plot+Neurons.R` (not copied
here).

## How to run

1. Edit `config.R` and set `PROJECT_ROOT` (and, if needed, `CELLRANGER_DIR`,
   `MACS_PATH`, `LDSC_HOME`, `METADATA_XLSX`) for your environment, or export the
   matching environment variables (`ASAP_MULTIOME_ROOT`, `ASAP_CELLRANGER_DIR`,
   `ASAP_MACS_PATH`, `LDSC_HOME`, `ASAP_METADATA_XLSX`).
2. Run an R script from this directory, e.g. `Rscript 02-Filtering-Clustering.R`.
   Every R script begins with `source("config.R")`, which sets the working
   directory to `PROJECT_ROOT`; all other paths are relative to it.
3. Shell scripts (`08-GWAS-LDSC-QTL.sh`, `11-bw-files-for-ATAC.sh`) read the same
   environment variables and can be run with `bash <script>.sh`.

`config.R` centralizes output directories, canonical object paths (`SEURAT_OBJECT`,
`MACS_COUNTS`), QTL/GWAS inputs (`GENOTYPE_DIR`, `COVARIATES_FILE`, `SNP_LOC`,
`GWAS_SUMSTATS`, `GWAS_CLUMPED`), the cell-type list (`CELL_TYPES`), and diagnosis
contrasts (`COMPARISONS`). Step-local run targets (e.g. `ct <- "Neurons"` in QTL
preprocessing, `CELL_TYPE <- "Microglia"` in eRegulon GSEA) remain at the top of
the relevant script — change them or wrap the script in a loop to process all cell
types.

## Environment setup

**Verified on this project:** the main analysis used **R 4.3.3 / Python 3.12.3**
(`single_cell_env`); SCENIC+ used a separate **R 4.4.1 / Python 3.11.3**
(`scenicplus` env). Do not mix these up — the submitted README’s “R 4.4.1 /
Python 3.11.3” applies to SCENIC+ only.

| Conda env | R | Python | Used for |
|-----------|---|--------|----------|
| `asap-multiome` (`environment.yml`) | 4.3.3 | 3.12.3 | Steps 01–11, 13 (Seurat, QTL, NEBULA, moloc, scMultiMap, …) |
| `asap-scenicplus` (`environment-scenicplus.yml`) | 4.4.1 | 3.11.3 | Step 12 SCENIC+ only |
| `asap-peer` (`environment-peer.yml`) | 3.4.1 | 3.5.5 | PEER blocks in steps 05/06 only |

```bash
# Main pipeline
conda env create -f environment.yml
conda activate asap-multiome
# GitHub-only R packages (see header comments in environment.yml):
# nebula, moloc, scMultiMap

# SCENIC+ (step 12)
conda env create -f environment-scenicplus.yml
conda activate asap-scenicplus

# PEER (steps 05/06)
conda env create -f environment-peer.yml
conda activate asap-peer
```

| File | Purpose |
|------|---------|
| `environment.yml` | Main `asap-multiome` env (R 4.3.3, Python 3.12.3) |
| `environment-scenicplus.yml` | SCENIC+ `asap-scenicplus` env (Python 3.11.3 + pip stack) |
| `environment-peer.yml` | PEER `asap-peer` env (R 3.4.1, steps 05/06 only) |

Install separately: **Cell Ranger ARC**, **FastQTL**, **LDSC** (`LDSC_HOME` in `config.R`).

## Scripts (one per main pipeline step)

| Script | Step | Summary |
|--------|------|---------|
| `config.R` | – | Shared paths, tool locations, and global parameters |
| `01-Combining-Samples.R` | 01 | Build per-sample multiome Seurat objects and merge |
| `02-Filtering-Clustering.R` | 02 | QC filtering, SCT/LSI, Harmony, WNN, doublets, annotation *(QC/UMAP plots retained)* |
| `03-Callpeaks-CellTypes.R` | 03 | MACS peaks, union peak set, counts, cell-type-specific peaks |
| `05-caQTL-preprocess.R` | 05 | Pseudobulk accessibility, PEER, FastQTL inputs (+ shell blocks) |
| `06-eQTL-preprocess.R` | 06 | Pseudobulk expression, PEER, FastQTL inputs (+ shell blocks) |
| `07-moloc-coloc.R` | 07 | moloc caQTL + eQTL + GWAS three-trait colocalization |
| `08-GWAS-LDSC-QTL.sh` | 08 | S-LDSC partitioned heritability (annotations, LD scores, h²) |
| `09-PD-differential-nebula-deg.R` | 09 | Single-cell NEBULA DEG across diagnosis contrasts |
| `10-scMultiMap-peak-gene.R` | 10 | scMultiMap peak–gene cis linking |
| `11-bw-files-for-ATAC.sh` | 11 | Per-cell-type normalized ATAC bigWig tracks |
| `13-1-eRegulon-GSEA-core.R` | 13 | eRegulon gene-set GSEA on NEBULA DEG rankings (per cell type) |

## Core-only trimming (what was removed)

| Script | Removed (non-core) |
|--------|---------------------|
| `03` | Genomic-distribution barplot, UpSet, rGREAT GO, accessibility heatmap |
| `05` | QC plots, per-chromosome matrix dumps, raw pseudobulk writes, PEER alpha diagnostics |
| `06` | Same as `05` (retains one PEER alpha PDF for factor-count sanity check) |
| `07` | LocusZoom input prep, stacked locus plots, motif analyses |
| `08` | Repetitive shell blocks consolidated into a reusable `run_sldsc` function |
| `09` | GO over-representation and GSEA enrichment on DEGs |
| `10` | Cross-cell-type summary and moloc-intersection sections |
| `13-1` | Plotting helpers and diagnosis-specific RSS prefilter |

Preserved throughout: analysis steps, statistical models, thresholds, and the
primary output files each script writes.

## Dependencies

Key R packages (via `environment.yml`; nebula/moloc/scMultiMap from GitHub — see yml header):

| Step | Packages |
|------|----------|
| 01–03 | Seurat, Signac, GenomicFeatures, EnsDb.Hsapiens.v86, dplyr, edgeR, harmony |
| 05–06 | edgeR, data.table, tidyverse; PEER via `asap-peer` env |
| 07 | moloc, data.table, GenomicRanges |
| 09 | nebula, RUVSeq, edgeR, Matrix |
| 10 | scMultiMap |
| 13-1 | data.table, clusterProfiler |

See also the parent `../README.md` for the full original pipeline dependency list.
