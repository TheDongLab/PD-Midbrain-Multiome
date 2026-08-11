#!/usr/bin/env python3
"""
SCENIC+ step 3 (per cell type): write Snakemake config YAMLs for motif enrichment + GRN inference.

For each cell type, writes a config pointing to cisTopic pickle, GEX h5ad, region sets,
cisTarget DBs, and distinct outputs under outs_by_celltype/<CellType>/.

Example:
  python scenicplus_step3_per_celltype.py --proj-dir ./scenicplus_project --celltypes Microglia --execute

Manual Snakemake:
  cd scplus_pipeline/Snakemake
  snakemake --cores 20 --configfile ../../configs_by_celltype/config_Microglia_all_peaks.yaml
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path
from typing import Any

import yaml


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def discover_celltype_dirs(base: Path, celltypes: list[str] | None) -> list[Path]:
    if not base.is_dir():
        raise FileNotFoundError(f"Missing step-1 dir: {base}")
    roots = []
    for d in sorted(base.iterdir()):
        if not d.is_dir() or d.name.startswith(("_", ".")):
            continue
        if celltypes:
            want = set(celltypes) | {sanitize_name(ct) for ct in celltypes}
            if d.name not in want:
                continue
        roots.append(d)
    return roots


def base_config_template(*, temp_dir: Path, n_cpu: int, biomart_host: str) -> dict[str, Any]:
    return {
        "params_general": {"temp_dir": str(temp_dir), "n_cpu": int(n_cpu), "seed": 666},
        "params_data_preparation": {
            "bc_transform_func": '"lambda x: x"',
            "is_multiome": True,
            "key_to_group_by": "",
            "nr_cells_per_metacells": 10,
            "direct_annotation": "Direct_annot",
            "extended_annotation": "Orthology_annot",
            "species": "hsapiens",
            "biomart_host": biomart_host,
            "search_space_upstream": "1000 150000",
            "search_space_downstream": "1000 150000",
            "search_space_extend_tss": "10 10",
        },
        "params_motif_enrichment": {
            "species": "homo_sapiens",
            "annotation_version": "v10nr_clust",
            "motif_similarity_fdr": 0.001,
            "orthologous_identity_threshold": 0.0,
            "annotations_to_use": "Direct_annot Orthology_annot",
            "fraction_overlap_w_dem_database": 0.4,
            "dem_max_bg_regions": 500,
            "dem_balance_number_of_promoters": True,
            "dem_promoter_space": 1000,
            "dem_adj_pval_thr": 0.05,
            "dem_log2fc_thr": 1.0,
            "dem_mean_fg_thr": 0.0,
            "dem_motif_hit_thr": 3.0,
            "fraction_overlap_w_ctx_database": 0.4,
            "ctx_auc_threshold": 0.005,
            "ctx_nes_threshold": 3.0,
            "ctx_rank_threshold": 0.05,
        },
        "params_inference": {
            "tf_to_gene_importance_method": "GBM",
            "region_to_gene_importance_method": "GBM",
            "region_to_gene_correlation_method": "SR",
            "order_regions_to_genes_by": "importance",
            "order_TFs_to_genes_by": "importance",
            "gsea_n_perm": 1000,
            "quantile_thresholds_region_to_gene": "0.85 0.90 0.95",
            "top_n_regionTogenes_per_gene": "5 10 15",
            "top_n_regionTogenes_per_region": "",
            "min_regions_per_gene": 0,
            "rho_threshold": 0.05,
            "min_target_genes": 10,
        },
    }


def build_output_data(out_cell_dir: Path, run_tag: str) -> dict[str, str]:
    d = out_cell_dir.resolve()
    def _f(suffix: str) -> str:
        return str(d / f"{run_tag}{suffix}")
    return {
        "combined_GEX_ACC_mudata": _f("_ACC_GEX.h5mu"),
        "dem_result_fname": _f("_dem_results.hdf5"),
        "ctx_result_fname": _f("_ctx_results.hdf5"),
        "output_fname_dem_html": _f("_dem_results.html"),
        "output_fname_ctx_html": _f("_ctx_results.html"),
        "cistromes_direct": _f("_cistromes_direct.h5ad"),
        "cistromes_extended": _f("_cistromes_extended.h5ad"),
        "tf_names": _f("_tf_names.txt"),
        "genome_annotation": _f("_genome_annotation.tsv"),
        "chromsizes": _f("_chromsizes.tsv"),
        "search_space": _f("_search_space.tsv"),
        "tf_to_gene_adjacencies": _f("_tf_to_gene_adj.tsv"),
        "region_to_gene_adjacencies": _f("_region_to_gene_adj.tsv"),
        "eRegulons_direct": _f("_eRegulon_direct.tsv"),
        "eRegulons_extended": _f("_eRegulons_extended.tsv"),
        "AUCell_direct": _f("_AUCell_direct.h5mu"),
        "AUCell_extended": _f("_AUCell_extended.h5mu"),
        "scplus_mdata": _f("_scplusmdata.h5mu"),
    }


def write_one_config(
    *,
    celltype_dir: Path,
    run_name: str,
    db_prefix: str,
    gex_h5ad: Path,
    motif_annot: Path,
    outs_cell_dir: Path,
    config_out: Path,
    temp_dir: Path,
    n_cpu: int,
    biomart_host: str,
) -> Path:
    work = celltype_dir / run_name
    cfg = base_config_template(temp_dir=temp_dir, n_cpu=n_cpu, biomart_host=biomart_host)
    cfg["input_data"] = {
        "cisTopic_obj_fname": str((work / f"{run_name}_cisTopicObject.pkl").resolve()),
        "GEX_anndata_fname": str(gex_h5ad.resolve()),
        "region_set_folder": str((work / "region_sets").resolve()),
        "ctx_db_fname": str((work / f"{db_prefix}.regions_vs_motifs.rankings.feather").resolve()),
        "dem_db_fname": str((work / f"{db_prefix}.regions_vs_motifs.scores.feather").resolve()),
        "path_to_motif_annotations": str(motif_annot.resolve()),
    }
    outs_cell_dir.mkdir(parents=True, exist_ok=True)
    cfg["output_data"] = build_output_data(outs_cell_dir, run_name)
    config_out.parent.mkdir(parents=True, exist_ok=True)
    with open(config_out, "w", encoding="utf-8") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
    print(f"Wrote {config_out}", flush=True)
    return config_out


def main() -> int:
    p = argparse.ArgumentParser(description="SCENIC+ step 3 per cell type: Snakemake configs")
    p.add_argument("--proj-dir", type=Path, default=Path("."))
    p.add_argument("--out-subdir", default="output_by_celltype")
    p.add_argument("--celltypes", nargs="*", default=None)
    p.add_argument("--layout", choices=("all-peaks", "per-chrom"), default="all-peaks")
    p.add_argument("--chromosomes", nargs="*", default=None)
    p.add_argument("--gex-h5ad", type=Path, default=None)
    p.add_argument("--motif-annotations", type=Path, default=None)
    p.add_argument("--configs-dir", type=Path, default=None)
    p.add_argument("--outs-subdir", default="outs_by_celltype")
    p.add_argument("--temp-dir", type=Path, default=None)
    p.add_argument("--n-cpu", type=int, default=20)
    p.add_argument("--biomart-host", default="http://www.ensembl.org")
    p.add_argument("--execute", action="store_true", help="Run snakemake for each config")
    p.add_argument("--snakemake-cwd", type=Path, default=None)
    p.add_argument("--snakemake-cores", type=int, default=20)
    args = p.parse_args()

    proj = args.proj_dir.resolve()
    step1_base = (proj / args.out_subdir).resolve()
    configs_dir = (args.configs_dir or proj / "configs_by_celltype").resolve()
    outs_root = (proj / args.outs_subdir).resolve()
    temp_dir = (args.temp_dir or proj / "tmp").resolve()
    gex = (args.gex_h5ad or proj / "adata_normalized.h5ad").resolve()
    motif_annot = (args.motif_annotations or proj / "motifs" / "motifs-v10-nr.hgnc.tbl").resolve()
    snakemake_cwd = (args.snakemake_cwd or proj / "scplus_pipeline" / "Snakemake").resolve()

    cell_dirs = discover_celltype_dirs(step1_base, args.celltypes)
    if not cell_dirs:
        print(f"No cell-type folders under {step1_base}", flush=True)
        return 1

    chroms = args.chromosomes or ([f"chr{i}" for i in range(1, 23)] if args.layout == "per-chrom" else None)
    written: list[Path] = []

    for cd in cell_dirs:
        out_cell = outs_root / cd.name
        if args.layout == "all-peaks":
            written.append(write_one_config(
                celltype_dir=cd, run_name="all_peaks", db_prefix="1kb_bg_all_peaks",
                gex_h5ad=gex, motif_annot=motif_annot, outs_cell_dir=out_cell,
                config_out=configs_dir / f"config_{cd.name}_all_peaks.yaml",
                temp_dir=temp_dir, n_cpu=args.n_cpu, biomart_host=args.biomart_host,
            ))
        else:
            for c in chroms or []:
                chrom = c if str(c).startswith("chr") else f"chr{c}"
                written.append(write_one_config(
                    celltype_dir=cd, run_name=chrom, db_prefix=f"1kb_bg_{chrom}",
                    gex_h5ad=gex, motif_annot=motif_annot, outs_cell_dir=out_cell,
                    config_out=configs_dir / f"config_{cd.name}_{chrom}.yaml",
                    temp_dir=temp_dir, n_cpu=args.n_cpu, biomart_host=args.biomart_host,
                ))

    if args.execute:
        if not snakemake_cwd.is_dir():
            print(f"Cannot --execute: missing {snakemake_cwd}", flush=True)
            return 1
        for cfg in written:
            cmd = ["snakemake", "--cores", str(args.snakemake_cores),
                   "--configfile", os.path.relpath(cfg, snakemake_cwd)]
            print(f"+ cd {snakemake_cwd} && {' '.join(cmd)}", flush=True)
            subprocess.run(cmd, cwd=str(snakemake_cwd), check=True)

    print(f"\nWrote {len(written)} config(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
