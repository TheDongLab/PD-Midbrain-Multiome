# PD Midbrain Multiome — Analysis Scripts

Analysis scripts for a single-nucleus multiome study of the human substantia
nigra in Parkinson’s disease (PD), incidental Lewy body disease (ILBD), and
healthy controls (HC). Each nucleus has paired snRNA-seq and snATAC-seq
(Cell Ranger ARC; genome build hg38).

The folder contains the core computational steps used in the manuscript: from
multiome object construction and cell-type annotation, through cis-QTL mapping
and GWAS integration, to differential expression, peak–gene linking, and
cell-type gene-regulatory networks (SCENIC+ / eRegulons).

---

## Biological workflow

```
snRNA + snATAC (Cell Ranger ARC)
      │
      ▼
[01] Build & merge multiome objects
      │
      ▼
[02] QC, integration, clustering, cell-type annotation
      │
      ▼
[03] Cell-type ATAC peak calling (MACS)
      │
      ├──────────────────┬──────────────────┐
      ▼                  ▼                  ▼
[05] cis-caQTL          [06] cis-eQTL       [09] Diagnosis DEGs (NEBULA)
      │                  │                  │
      └────────┬─────────┘                  │
               ▼                            │
         [07] Colocalization                │
              (caQTL + eQTL + PD GWAS)      │
               │                            │
               ▼                            │
         [08] Cell-type QTL enrichment      │
              of PD heritability (S-LDSC)   │
                                            │
[10] Peak–gene cis links (scMultiMap)       │
                                            │
[12] Cell-type GRNs (SCENIC+)               │
      │                                     │
      └──────────────┬──────────────────────┘
                     ▼
              [13] eRegulon enrichment
                   in diagnosis DEG ranks
```

**Upstream (01–03)** defines cell types and the accessible peak universe.
**QTL / GWAS modules (05–08)** ask which regulatory variants and cell types
contribute to PD risk. **DEG and linking (09–10)** and **GRN / eRegulon
(12–13)** connect disease-associated expression changes to candidate enhancers
and transcription-factor programs.

---

## Analysis modules

### A. Multiome atlas of the substantia nigra

| Step | Script | Biological aim |
|------|--------|----------------|
| 01 | `Combining-Samples.R` | Load per-donor paired RNA + ATAC, compute QC metrics, merge into one multiome object |
| 02 | `Filtering-Clustering.R` | Filter low-quality nuclei, integrate across donors (SCT / LSI / Harmony / WNN), cluster, and annotate major midbrain cell types |
| 03 | `Callpeaks-CellTypes.R` | Call ATAC peaks per cell type and in bulk; build a union peak set and peak×cell counts for downstream accessibility analyses |

After these steps, analyses use a labeled multiome object and cell-type-aware
peak counts. Major cell types include DA, GABA, and Glu neurons, microglia,
astrocytes, oligodendrocytes, OPCs, endothelial cells, pericytes, macrophages,
and T cells.

### B. Cell-type cis-QTLs and PD GWAS integration

| Step | Script | Biological aim |
|------|--------|----------------|
| 05 | `caQTL-preprocess.R` | Pseudobulk chromatin accessibility by donor and cell type; map cis-caQTLs (FastQTL) |
| 06 | `eQTL-preprocess.R` | Pseudobulk gene expression by donor and cell type; map cis-eQTLs (FastQTL) |
| 07 | `moloc-coloc.R` | Three-trait colocalization of PD GWAS with cell-type caQTL and eQTL signals at lead risk loci |
| 08 | `GWAS-LDSC-QTL.sh` | Stratified LD-score regression: enrichment of PD heritability in cell-type QTL annotations |

Together these steps identify cell-type regulatory variants that may mediate PD
GWAS effects via chromatin accessibility and/or expression.

### C. Disease differential expression and enhancer–gene links

| Step | Script | Biological aim |
|------|--------|----------------|
| 09 | `PD-differential-nebula-deg.R` | Single-cell DEGs between PD, ILBD, and HC within each cell type (NEBULA mixed models), adjusting for age, sex, PMI, batch, and technical factors |
| 10 | `scMultiMap-peak-gene.R` | Link accessible peaks to candidate target genes in *cis* (±1 Mb) per cell type |

Diagnosis contrasts: **PD vs HC**, **ILBD vs HC**, **PD vs ILBD**.

### D. Cell-type gene regulatory networks

| Step | Script | Biological aim |
|------|--------|----------------|
| 12.1 | `scenicplus_step1_per_celltype.py` | Topic modeling of accessibility (pycisTopic), differential accessible regions, and gene activity per cell type |
| 12.2 | `scenicplus_step2_per_celltype.py` | Motif / cisTarget databases for cell-type region sets |
| 12.3 | `scenicplus_step3_per_celltype.py` | Configure SCENIC+ GRN inference (enhancer-driven eRegulons) |
| 13 | `eRegulon-GSEA-core.R` | Test whether SCENIC+ eRegulon target genes are enriched among diagnosis DEG rankings |

This module recovers cell-type transcription-factor programs (eRegulons) and
asks which programs shift along the HC → ILBD → PD axis.

---

---

## Scope

Scripts retain **core analysis steps** for reproducibility. Omitted from most
modules (see original project scripts if needed): package installation,
interactive debugging, many QC/summary figures, and downstream interpretation
panels (e.g. LocusZoom, motif-break analyses, broad GO collections beyond
eRegulon GSEA). 
