#!/usr/bin/env python3
"""
SCENIC+ step 1 (per cell type): pycisTopic LDA, topic binarization, DARs, gene activity.

Subsets multiome ATAC data by cell type, optionally filters peaks by prevalence,
runs cisTopic per cell type (whole-genome or per-chromosome), and writes region sets
for downstream cisTarget / Snakemake.

Example:
  python scenicplus_step1_per_celltype.py \\
    --proj-dir ./scenicplus_project \\
    --celltypes Microglia Astrocytes \\
    --no-chrom-split \\
    --group-col Diagnosis
"""

from __future__ import annotations

import argparse
import os
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pybiomart as pbm
import pyranges as pr
import scanpy as sc
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.diff_features import find_diff_features, impute_accessibility
from pycisTopic.gene_activity import get_gene_activity
from pycisTopic.lda_models import evaluate_models, run_cgs_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.topic_qc import compute_topic_metrics
from pycistarget.utils import region_names_to_coordinates
from scipy.sparse import csr_matrix, issparse

os.environ.setdefault("RAY_USAGE_STATS_ENABLED", "0")


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def peak_prevalence_min_cells(n_cells: int) -> int:
    return max(20, int(np.ceil(0.01 * n_cells))) if n_cells > 0 else 1


def subset_peaks_by_cell_prevalence(adata: sc.AnnData, log_label: str) -> sc.AnnData:
    n_cells = adata.shape[0]
    thr = peak_prevalence_min_cells(n_cells)
    X = adata.X
    col_nnz = np.diff(X.tocsc().indptr) if issparse(X) else np.count_nonzero(X, axis=0)
    keep = col_nnz >= thr
    print(f"[{log_label}] Peak filter: keep {keep.sum()}/{len(keep)} (≥{thr} cells)", flush=True)
    return adata[:, keep].copy()


def subsample_cells(adata: sc.AnnData, max_cells: int | None, seed: int, label: str) -> sc.AnnData:
    n = adata.shape[0]
    if max_cells is None or n <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=max_cells, replace=False))
    print(f"[{label}] Subsample {max_cells}/{n} cells", flush=True)
    return adata[idx].copy()


def parse_regions(region_list: list[str]) -> pd.DataFrame:
    rows = []
    for region in region_list:
        parts = region.replace(":", "-").split("-")
        if len(parts) == 3:
            rows.append([parts[0], int(parts[1]), int(parts[2])])
    return pd.DataFrame(rows, columns=["Chromosome", "Start", "End"])


def add_cell_metadata(cistopic_obj, adata_subset, metadata_path: Path | None, ensure_cols: tuple[str, ...] = ()):
    if metadata_path and metadata_path.is_file():
        meta = pd.read_csv(metadata_path, sep="\t", index_col=0)
        inter = cistopic_obj.cell_data.index.intersection(meta.index)
        if len(inter) > 0:
            cistopic_obj.add_cell_data(meta.loc[inter])
            for col in ensure_cols:
                if col not in cistopic_obj.cell_data.columns and col in adata_subset.obs.columns:
                    cistopic_obj.cell_data[col] = adata_subset.obs[col].reindex(cistopic_obj.cell_data.index)
            return
    cistopic_obj.add_cell_data(adata_subset.obs.copy())


def write_region_sets(region_dict: dict, folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for topic, regions in region_dict.items():
        pr.PyRanges(region_names_to_coordinates(regions.index.tolist())).to_csv(
            folder / f"{topic}.bed", sep="\t", header=False
        )


def process_one_peak_group(
    adata_atac_cell: sc.AnnData,
    celltype_label: str,
    group_name: str,
    celltype_out_dir: Path,
    metadata_path: Path | None,
    group_col: str,
    n_topics: list[int],
    lda_n_cpu: int,
    proj_dir: Path,
    peak_mask=None,
    dar_n_cpu: int = 5,
    dar_adjpval_thr: float = 0.05,
    dar_log2fc_thr: float | None = None,
):
    """Run full pycisTopic workflow for one peak group (all_peaks or chrN)."""
    if dar_log2fc_thr is None:
        dar_log2fc_thr = float(np.log2(1.5))

    print(f"[{celltype_label}] {group_name} …", flush=True)
    out_dir = celltype_out_dir / group_name
    out_dir.mkdir(parents=True, exist_ok=True)

    chr_mask = adata_atac_cell.var_names.str.startswith(group_name + ":") if peak_mask is None else peak_mask
    if np.asarray(chr_mask).sum() == 0:
        print(f"[{celltype_label}] Skip {group_name}: no peaks", flush=True)
        return

    adata_chr = adata_atac_cell[:, chr_mask].copy()
    X = adata_chr.X.toarray() if isinstance(adata_chr.X, csr_matrix) else adata_chr.X
    df_matrix = pd.DataFrame(X, index=adata_chr.obs.index, columns=adata_chr.var.index)

    cistopic_obj = create_cistopic_object(df_matrix.T)
    add_cell_metadata(cistopic_obj, adata_chr, metadata_path, ensure_cols=(group_col,))

    bed_df = parse_regions(cistopic_obj.region_names)
    bed_df.to_csv(out_dir / f"{group_name}_cistopic_regions.bed", sep="\t", header=False, index=False)

    models = run_cgs_models(cistopic_obj, n_topics=n_topics, n_cpu=lda_n_cpu, n_iter=500, eta=0.1)
    with open(out_dir / "models.pkl", "wb") as f:
        pickle.dump(models, f)

    model = evaluate_models(models, return_model=True, save=str(out_dir / f"{group_name}_model_selection.pdf"))
    cistopic_obj.add_LDA_model(model)
    with open(out_dir / f"{group_name}_cisTopicObject.pkl", "wb") as f:
        pickle.dump(cistopic_obj, f)

    region_top3k = binarize_topics(cistopic_obj, method="otsu", ntop=3000, plot=True,
                                   save=str(out_dir / f"{group_name}_binarized_topics_top3k.pdf"))
    region_otsu = binarize_topics(cistopic_obj, method="otsu", plot=True, num_columns=5,
                                  save=str(out_dir / f"{group_name}_binarized_topics_otsu.pdf"))
    cell_topic = binarize_topics(cistopic_obj, target="cell", method="li", nbins=60)

    for name, obj in [
        (f"{group_name}_binarized_cell_topic.pkl", cell_topic),
        (f"{group_name}_binarized_topic_region_top_3k.pkl", region_top3k),
        (f"{group_name}_binarized_topic_region_otsu.pkl", region_otsu),
    ]:
        with open(out_dir / name, "wb") as f:
            pickle.dump(obj, f)

    write_region_sets(region_top3k, out_dir / "region_sets" / "Topics_top_3k")
    write_region_sets(region_otsu, out_dir / "region_sets" / "Topics_otsu")

    with open(out_dir / f"{group_name}_Topic_qc_metrics.pkl", "wb") as f:
        pickle.dump(compute_topic_metrics(cistopic_obj), f)

    imputed = impute_accessibility(cistopic_obj, scale_factor=10**6)
    with open(out_dir / f"{group_name}_Imputed_accessibility.pkl", "wb") as f:
        pickle.dump(imputed, f)

    markers_dict = {}
    if group_col in cistopic_obj.cell_data.columns:
        markers_dict = find_diff_features(
            cistopic_obj, imputed, variable=group_col,
            adjpval_thr=dar_adjpval_thr, log2fc_thr=dar_log2fc_thr, n_cpu=dar_n_cpu,
        )
        with open(out_dir / f"{group_name}_DARs.pkl", "wb") as f:
            pickle.dump(markers_dict, f)

    dar_folder = out_dir / "region_sets" / "DARs_group"
    dar_folder.mkdir(parents=True, exist_ok=True)
    for grp, df in markers_dict.items():
        if df.empty:
            continue
        region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
            dar_folder / f"{sanitize_name(grp)}.bed", sep="\t", header=False, index=False
        )

    # Gene activity (Ensembl protein-coding annotation)
    dataset = pbm.Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
    annot = dataset.query(attributes=[
        "chromosome_name", "start_position", "end_position", "strand",
        "external_gene_name", "transcription_start_site", "transcript_biotype",
    ])
    annot.columns = ["Chromosome", "Start", "End", "Strand", "Gene",
                     "Transcription_Start_Site", "Transcript_type"]
    annot = annot[annot["Transcript_type"] == "protein_coding"].dropna()
    annot["Chromosome"] = annot["Chromosome"].astype(str).apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    annot["Strand"] = annot["Strand"].astype(str).replace({"1": "+", "-1": "-"})

    chromsizes = pd.read_csv(proj_dir / "hg38.chrom.sizes", sep="\t", header=None, names=["Chromosome", "End"])
    chromsizes["Start"] = 0
    gene_act, _ = get_gene_activity(imputed, pr.PyRanges(annot), pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]]), scale_factor=1)
    with open(out_dir / f"{group_name}_Gene_activity.pkl", "wb") as f:
        pickle.dump(gene_act, f)

    print(f"[{celltype_label}] Done {group_name}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="SCENIC+ step 1 per cell type (pycisTopic)")
    p.add_argument("--proj-dir", type=Path, default=Path("."))
    p.add_argument("--input-h5ad", default="merged_seurat4.h5ad")
    p.add_argument("--metadata-tsv", default="celldata_mel.tsv")
    p.add_argument("--celltype-col", default="celltype")
    p.add_argument("--group-col", default="Diagnosis", help="Column for DAR contrasts (≥2 levels)")
    p.add_argument("--celltypes", nargs="*", default=None)
    p.add_argument("--chromosomes", nargs="*", default=None)
    p.add_argument("--no-chrom-split", action="store_true", help="One whole-genome model per cell type")
    p.add_argument("--min-cells", type=int, default=200)
    p.add_argument("--max-cells-per-celltype", type=int, default=None)
    p.add_argument("--subsample-seed", type=int, default=666)
    p.add_argument("--lda-n-cpu", type=int, default=2)
    p.add_argument("--n-topics", nargs="*", type=int, default=[10, 15, 20, 25, 30, 35, 40])
    p.add_argument("--out-subdir", default="output_by_celltype")
    p.add_argument("--dar-n-cpu", type=int, default=5)
    p.add_argument("--dar-adjpval", type=float, default=0.05, dest="dar_adjpval_thr")
    p.add_argument("--skip-peak-prevalence-filter", action="store_true")
    p.add_argument(
        "--merge-neuron-subtypes", action="store_true",
        help="Merge DA/GABA/Glu labels into 'Neurons'",
    )
    args = p.parse_args()

    proj_dir = args.proj_dir.resolve()
    h5ad_path = proj_dir / args.input_h5ad
    metadata_path = proj_dir / args.metadata_tsv
    out_dir = proj_dir / args.out_subdir
    chromosomes = args.chromosomes or [f"chr{i}" for i in range(1, 23)]

    adata_atac = sc.read_h5ad(h5ad_path)
    if args.celltype_col not in adata_atac.obs.columns and metadata_path.is_file():
        meta = pd.read_csv(metadata_path, sep="\t", index_col=0)
        adata_atac.obs[args.celltype_col] = meta[args.celltype_col].reindex(adata_atac.obs_names)

    celltypes = adata_atac.obs[args.celltype_col].astype(str)
    if args.merge_neuron_subtypes:
        neuron_labels = {"da neurons", "gaba neurons", "glu neurons"}
        celltypes = celltypes.where(~celltypes.str.lower().isin(neuron_labels), "Neurons")

    selected = sorted(celltypes.unique()) if not args.celltypes else args.celltypes
    print(f"Cell types: {selected}", flush=True)

    for ct in selected:
        mask = celltypes == ct
        n = int(mask.sum())
        if n < args.min_cells:
            print(f"[{ct}] Skip: n={n} < {args.min_cells}", flush=True)
            continue

        ad = adata_atac[mask].copy()
        ad = subsample_cells(ad, args.max_cells_per_celltype, args.subsample_seed, ct)
        if not args.skip_peak_prevalence_filter:
            ad = subset_peaks_by_cell_prevalence(ad, ct)

        ct_dir = out_dir / sanitize_name(ct)
        ct_dir.mkdir(parents=True, exist_ok=True)
        meta = metadata_path if metadata_path.is_file() else None

        if args.no_chrom_split:
            process_one_peak_group(
                ad, ct, "all_peaks", ct_dir, meta, args.group_col,
                args.n_topics, args.lda_n_cpu, proj_dir,
                peak_mask=np.ones(ad.shape[1], dtype=bool),
                dar_n_cpu=args.dar_n_cpu, dar_adjpval_thr=args.dar_adjpval_thr,
            )
        else:
            for chrom in chromosomes:
                process_one_peak_group(
                    ad, ct, chrom, ct_dir, meta, args.group_col,
                    args.n_topics, args.lda_n_cpu, proj_dir,
                    dar_n_cpu=args.dar_n_cpu, dar_adjpval_thr=args.dar_adjpval_thr,
                )

    print("All cell types processed.", flush=True)


if __name__ == "__main__":
    main()
