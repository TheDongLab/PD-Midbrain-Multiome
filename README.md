# PD Midbrain Multiome

Analysis code and selected results for a single-nucleus multiome atlas of the
human substantia nigra and midbrain across Parkinson's disease progression.

**Repository:** [TheDongLab/PD-Midbrain-Multiome](https://github.com/TheDongLab/PD-Midbrain-Multiome)  
**Reference genome:** GRCh38/hg38  
**Assay:** 10x Genomics Single Cell Multiome ATAC + Gene Expression  
**Study population:** 88 postmortem donors classified as healthy control (HC),
incidental Lewy body disease (ILBD), or Parkinson's disease (PD)  
**README last updated:** 2026-09-09

### Permanent archives

- **Analysis code:** [Zenodo record 21924320](https://zenodo.org/records/21924320)
- **Processed datasets:** [Zenodo record 21937673](https://zenodo.org/records/21937673)
- **Key Resource Table:** [Zenodo record 21938197](https://zenodo.org/records/21938197)

## Project overview

This project maps cell-type-specific gene regulation and genetic risk in PD
using paired single-nucleus RNA sequencing (snRNA-seq) and single-nucleus ATAC
sequencing (snATAC-seq). The analysis includes:

- sample-level quality control, integration, clustering, and cell-type
  annotation;
- cell-type-specific chromatin peak calling;
- chromatin-accessibility and expression quantitative trait locus mapping
  (caQTL and eQTL);
- PD GWAS heritability partitioning and three-trait colocalization;
- diagnosis-associated differential gene expression;
- peak-gene linking; and
- enhancer-driven gene-regulatory network inference and eRegulon enrichment.

The code is intended to document and reproduce the core computational analyses
associated with the study. It is not a general-purpose software package. Most
scripts are designed for a high-performance computing environment and require
controlled-access inputs that are not stored in GitHub.

## Data and outputs

### Data used by the analysis

The primary data are paired snRNA-seq and snATAC-seq profiles generated from
postmortem midbrain tissue. The pipeline also uses donor metadata, imputed
genotypes, PD GWAS summary statistics, GENCODE v47 gene annotations, the ENCODE
hg38 blacklist, and 1000 Genomes European-ancestry LD reference panels.

Raw FASTQ files, Cell Ranger ARC outputs, donor metadata, and individual-level
genotypes are not included in the public GitHub repository. Access to the
human-subject data is controlled. The raw multiome data are available through
the **ASAP CRN Cloud, dataset “Midbrain Multiome,” subject to approval and the
applicable ASAP Access & Use Policy**. Publicly shareable processed datasets are
available in [Zenodo record 21937673](https://zenodo.org/records/21937673).
Large intermediate objects may be requested from the study contact listed
below.

The expected external inputs and their required fields are documented in
[`Multiome_Midbrain_Notebook.txt`](Multiome_Midbrain_Notebook.txt). Reused
resources and identifiers are listed in
[`Key_Resource_Table.csv`](Key_Resource_Table.csv) and archived separately in
[Zenodo record 21938197](https://zenodo.org/records/21938197).

### Results included in the Zenodo deposit

The local Zenodo deposit contains selected analysis-ready result tables:

| Path | Contents | Format |
|---|---|---|
| `results/caQTL/` | Significant cell-type caQTL associations | CSV |
| `results/eQTL/` | Significant cell-type eQTL associations | CSV |
| `results/moloc/` | Three-trait caQTL/eQTL/PD-GWAS colocalization results | Tab-delimited text |
| `results/LDSC/` | Cell-type caQTL/eQTL partitioned-heritability results | Tab-delimited text |
| `results/DEG_nebula/` | NEBULA differential-expression results for PD, ILBD, and HC contrasts | TSV |
| `results/scenicplus/` | Direct enhancer-driven eRegulon region-gene-TF relationships | TSV |

Column definitions for every deposited table are provided in
[`Multiome_Midbrain_Notebook.txt`](Multiome_Midbrain_Notebook.txt). The full
FastQTL nominal association tables (approximately 640 GB) and large Seurat and
QTL objects are omitted because of size. The `scripts.zip` and `results.zip`
files mirror the corresponding directories for convenient Zenodo download.

## Analysis workflow

```text
Cell Ranger ARC outputs
        |
        v
Build and merge multiome objects
        |
        v
QC, integration, clustering, and annotation
        |
        v
Cell-type ATAC peak calling
        |
        +--------------------+--------------------+
        |                    |                    |
        v                    v                    v
      caQTL                eQTL          Differential expression
        |                    |                and peak-gene links
        +----------+---------+                    |
                   v                              v
        PD GWAS colocalization             SCENIC+ eRegulons
        and S-LDSC enrichment                     |
                                                  v
                                      eRegulon enrichment in DEGs
```

The following table identifies the function of each executable script, the
main input, and the main output. Detailed step-level notes are available in
[`scripts/README.md`](scripts/README.md).

| Order | Script | Function | Main output |
|---:|---|---|---|
| 1 | `Combining-Samples.R` | Build per-sample Seurat/Signac objects and merge donors | Merged multiome object |
| 2 | `Filtering-Clustering.R` | QC, SCT/LSI processing, Harmony integration, WNN clustering, doublet removal, and annotation | Labeled Seurat object and QC/UMAP outputs |
| 3 | `Callpeaks-CellTypes.R` | Call MACS peaks by cell type and create a union peak assay | Union peaks and peak-by-cell count matrix |
| 4 | `caQTL-preprocess.R` | Create donor pseudobulk accessibility matrices, PEER factors, and FastQTL inputs | caQTL inputs and association results |
| 5 | `eQTL-preprocess.R` | Create donor pseudobulk expression matrices, PEER factors, and FastQTL inputs | eQTL inputs and association results |
| 6 | `moloc-coloc.R` | Colocalize caQTL, eQTL, and PD GWAS signals | Cell-type moloc tables |
| 7 | `GWAS-LDSC-QTL.sh` | Estimate PD-heritability enrichment in QTL annotations with S-LDSC | Partitioned-heritability tables |
| 8 | `PD-differential-nebula-deg.R` | Fit diagnosis-associated single-cell negative-binomial mixed models | DEG tables for three diagnosis contrasts |
| 9 | `scMultiMap-peak-gene.R` | Link accessible peaks to candidate target genes within 1 Mb | Significant peak-gene links |
| 10a | `scenicplus_step1_per_celltype.py` | Run pycisTopic, topic binarization, DAR detection, and gene activity | Per-cell-type cisTopic objects and regions |
| 10b | `scenicplus_step2_per_celltype.py` | Build cisTarget motif databases | Motif ranking and score databases |
| 10c | `scenicplus_step3_per_celltype.py` | Write and launch SCENIC+ workflow configurations | Enhancer-driven eRegulons |
| 11 | `eRegulon-GSEA-core.R` | Test eRegulon target sets against NEBULA DEG rankings | GSEA tables and R objects |

Steps 1-3 are sequential. After a labeled multiome object and the union peak set
are available, the QTL, differential-expression, peak-gene, and SCENIC+
branches can be run as their required inputs become available.

## System requirements and dependencies

The analysis was performed with **R 4.4.1** and **Python 3.11.3** in a Linux
high-performance computing environment. Other operating systems and software
versions have not been tested. Full reruns are computationally intensive;
memory, CPU, storage, and scheduler requirements depend on the selected module
and input size. The complete dataset requires substantially more storage than
this deposit, including approximately 640 GB for full nominal QTL tables.

### R packages

| Package | Version | Package | Version |
|---|---:|---|---:|
| Seurat | 5.2.1 | Signac | 1.14.0 |
| GenomicFeatures | 1.54.4 | EnsDb.Hsapiens.v86 | 2.99.0 |
| GenomicRanges | 1.54.1 | rtracklayer | 1.62.0 |
| harmony | 1.2.3 | scDblFinder | 1.16.0 |
| BiocParallel | 1.36.0 | hdf5r | 1.3.12 |
| mbkmeans | 1.18.0 | cluster | 2.1.8.1 |
| openxlsx | 4.2.8 | ggplot2 | 3.5.1 |
| dplyr | 1.1.4 | tidyverse | 2.0.0 |
| data.table | 1.17.0 | edgeR | 4.0.16 |
| matrixStats | 1.5.0 | Matrix | 1.6.5 |
| peer | 1.0 | moloc | 0.1.0 |
| nebula | 1.5.3 | RUVSeq | 1.36.0 |
| scMultiMap | 1.0.1 | clusterProfiler | 4.10.1 |
| RColorBrewer | 1.1.3 |  |  |

### Python packages

| Package | Version | Package | Version |
|---|---:|---|---:|
| scanpy | 1.8.2 | pycisTopic | 1.0.1.dev67 |
| pycistarget | 1.1 | scenicplus | 1.0a2 |
| pybiomart | 0.2.0 | pyranges | 0.0.111 |
| numpy | 1.26.4 | pandas | 1.5.0 |
| scipy | 1.12.0 | PyYAML | 6.0.1 |

### Command-line software

| Software | Version |
|---|---:|
| Cell Ranger ARC | 2.0.2 |
| MACS3 | 3.0.3 |
| FastQTL | 2.184_gtex |
| LDSC | 1.0.1 |

The listed versions describe the tested computational environment; dependency
installation is not automated in this repository. Several tools and packages
require their own installation procedures or environment-specific builds.

## Installation

Clone the repository and enter its directory:

```bash
git clone https://github.com/TheDongLab/PD-Midbrain-Multiome.git
cd PD-Midbrain-Multiome
```

Install R, Python, and the dependencies listed above using the package managers
appropriate for your system. SCENIC+, pycisTopic, FastQTL, and LDSC should be
installed in isolated environments following their upstream documentation.
Before a full run, confirm that `Rscript`, `python`, `macs3`, and the required
QTL/GWAS tools resolve from the command line.

## Configuration and data placement

All R scripts source `scripts/config.R`. Either edit its user-settings block or
define the following environment variables:

| Variable | Purpose |
|---|---|
| `ASAP_MULTIOME_ROOT` | Root of the complete analysis project tree |
| `ASAP_CELLRANGER_DIR` | Directory containing per-sample Cell Ranger ARC outputs |
| `ASAP_MACS_PATH` | Path to the MACS3 executable |
| `ASAP_METADATA_XLSX` | Path to the donor/sample metadata workbook |
| `LDSC_HOME` | LDSC installation directory |

Example:

```bash
export ASAP_MULTIOME_ROOT=/path/to/ASAP_Multiome
export ASAP_CELLRANGER_DIR=/path/to/cellranger_arc_count
export ASAP_MACS_PATH=/path/to/bin/macs3
export ASAP_METADATA_XLSX=/path/to/midbrain_metadata.xlsx
export LDSC_HOME=/path/to/ldsc
```

The full project tree must also contain the genotypes, covariates, annotations,
GWAS inputs, reference panels, and intermediate files named in `config.R` and
[`Multiome_Midbrain_Notebook.txt`](Multiome_Midbrain_Notebook.txt). The small
set of public result tables in this deposit is not sufficient to rerun the
pipeline from raw reads.

## Running the code

Run commands from the `scripts/` directory because the R scripts load
`config.R` with a relative path:

```bash
cd scripts

# Sequential atlas construction
Rscript Combining-Samples.R
Rscript Filtering-Clustering.R
Rscript Callpeaks-CellTypes.R

# Examples of downstream modules
Rscript PD-differential-nebula-deg.R
Rscript scMultiMap-peak-gene.R
bash GWAS-LDSC-QTL.sh

# Example SCENIC+ preprocessing for one cell type
python scenicplus_step1_per_celltype.py \
  --proj-dir /path/to/scenicplus_project \
  --celltypes Microglia
```

Some scripts are templates or contain stage-specific blocks:

- `caQTL-preprocess.R`, `eQTL-preprocess.R`, and `moloc-coloc.R` use pooled
  neurons as the documented example; change or loop over the cell-type variable
  for atlas-wide analysis.
- The QTL preprocessing scripts include FastQTL shell blocks that must be run
  separately from their R preprocessing sections.
- `GWAS-LDSC-QTL.sh` is a staged, long-running workflow and should be adapted to
  the local scheduler and reference-panel locations.
- `eRegulon-GSEA-core.R` uses Microglia as its default cell-type template.

Expected outputs and output directories are described in the header of every
script. Successful execution should create non-empty files with the schemas
listed in [`Multiome_Midbrain_Notebook.txt`](Multiome_Midbrain_Notebook.txt).

## Repository structure

```text
.
├── README.md                         # Repository overview and reuse instructions
├── LICENSE                           # MIT License for the analysis code
├── Key_Resource_Table.csv            # Datasets, software, identifiers, and reuse status
├── Multiome_Midbrain_README.txt      # Detailed Zenodo deposit narrative
├── Multiome_Midbrain_metadata.txt    # Study and deposit metadata
├── Multiome_Midbrain_Notebook.txt    # Data dictionary and required external inputs
├── scripts/
│   ├── README.md                     # Analysis workflow details
│   ├── config.R                      # Central paths and global parameters
│   ├── *.R                           # R analysis modules
│   ├── *.py                          # SCENIC+ modules
│   └── *.sh                          # S-LDSC workflow
├── results/
│   ├── README.txt
│   ├── caQTL/
│   ├── eQTL/
│   ├── moloc/
│   ├── LDSC/
│   ├── DEG_nebula/
│   └── scenicplus/
├── scripts.zip                       # Archived copy for Zenodo distribution
└── results.zip                       # Archived copy for Zenodo distribution
```

## Known limitations

- Controlled-access human data and individual-level genotypes cannot be
  redistributed in the public repository.
- The repository does not contain a small public test dataset or an automated
  end-to-end test. Access-controlled inputs are required for a complete run.
- Scripts preserve the core computational analyses, but some exploratory plots,
  interactive debugging blocks, and downstream interpretation panels were
  removed from the cleaned release.
- Paths, scheduler commands, and resource requests must be adapted to the local
  computing environment.
- The selected Zenodo results do not include every intermediate or every output
  produced by the scripts.

## Citation

The code used to analyze the data is permanently archived in
[Zenodo record 21924320](https://zenodo.org/records/21924320). The manuscript
citation was pending when this README was prepared and should be added when it
becomes available.

> **Software citation:** Lu, M., Dong, X., et al. (2026). *PD Midbrain Multiome*
> [Computer software]. Zenodo.
> [https://zenodo.org/records/21924320](https://zenodo.org/records/21924320).

> **Manuscript citation:** **pending**.

When reusing the code, cite both the archived, version-specific Zenodo release
and the associated manuscript. The GitHub URL may be included to direct readers
to subsequent updates, but it should not replace the persistent software DOI.

## License and data-use terms

The analysis code is released under the **MIT License**; see
[`LICENSE`](LICENSE). Data and third-party resources are not relicensed by the
code license. They remain subject to their source licenses, data use agreements,
consent restrictions, and the ASAP Access & Use Policy. Users are responsible
for obtaining all required approvals before accessing controlled human-subject
data.

## Contact

For questions about the study, code, or access to omitted intermediate files:

**Xianjun Dong, PhD**  
Associate Professor, Departments of Neurology and Biomedical Informatics and
Data Science  
Yale School of Medicine, Stephen & Denise Adams Center for Parkinson's Disease
Research  
Email: [xianjun.dong@yale.edu](mailto:xianjun.dong@yale.edu)  
ORCID: [0000-0002-8052-9320](https://orcid.org/0000-0002-8052-9320)

**Code submitter:** Mingming Lu  
Email: [lumm123123123@gmail.com](mailto:lumm123123123@gmail.com)

## Acknowledgments

This project was conducted by the Dong Lab at Yale School of Medicine. The
authors acknowledge Aligning Science Across Parkinson's (ASAP), the study
participants and tissue donors, the contributing brain banks and data
providers, and the developers and maintainers of the open-source software and
reference resources listed in [`Key_Resource_Table.csv`](Key_Resource_Table.csv).
