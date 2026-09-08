#!/usr/bin/env python3
"""
SCENIC+ step 2 (per cell type): build cisTarget motif databases from step-1 cisTopic region BEDs.

Expects:
  {proj_dir}/{out_subdir}/{CellType}/all_peaks/all_peaks_cistopic_regions.bed
  or per-chrom: .../chr1/chr1_cistopic_regions.bed

Uses create_cisTarget_databases scripts from the SCENIC+ bundle (FASTA padding + motif scoring).

Example:
  python scenicplus_step2_per_celltype.py --proj-dir ./scenicplus_project --celltypes Microglia
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def _motif_filenames(singletons_dir: Path) -> list[str]:
    names = sorted(p.name for p in singletons_dir.iterdir() if p.is_file() and p.suffix == ".cb")
    if not names:
        raise FileNotFoundError(f"No .cb motif files under {singletons_dir}")
    return names


def resolve_cbust_path(cistarget_scripts_dir: Path, override: Path | None) -> str:
    if override:
        p = override.expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
        raise FileNotFoundError(f"--cbust not executable: {p}")
    bundled = cistarget_scripts_dir / "cbust"
    if bundled.is_file() and os.access(bundled, os.X_OK):
        return str(bundled)
    wh = shutil.which("cbust")
    if wh:
        return wh
    raise FileNotFoundError('Cluster-Buster "cbust" not found on PATH or next to cisTarget scripts.')


def run_cistarget_for_bed(
    *,
    region_bed: Path,
    work_dir: Path,
    database_prefix: str,
    genome_fasta: Path,
    chromsizes: Path,
    cistarget_scripts_dir: Path,
    motif_singletons_dir: Path,
    threads: int,
    bg_padding: int,
    dry_run: bool,
    cbust_path: Path | None,
) -> None:
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta_sh = cistarget_scripts_dir / "create_fasta_with_padded_bg_from_bed.sh"
    db_py = cistarget_scripts_dir / "create_cistarget_motif_databases.py"
    for f in (fasta_sh, db_py, region_bed):
        if not f.is_file():
            raise FileNotFoundError(f"Missing {f}")

    fasta_tag = database_prefix.removeprefix("1kb_bg_")
    fasta_out = work_dir / f"hg38.1kb_bg_padding_{fasta_tag}.fa"
    motif_list = work_dir / f"motifs_{fasta_tag}.txt"
    motif_list.write_text("\n".join(_motif_filenames(motif_singletons_dir)) + "\n", encoding="utf-8")

    cmd_fasta = ["bash", str(fasta_sh), str(genome_fasta), str(chromsizes), str(region_bed), str(fasta_out), str(bg_padding), "yes"]
    cbust = resolve_cbust_path(cistarget_scripts_dir, cbust_path)
    cmd_db = [
        sys.executable, str(db_py), "-f", str(fasta_out), "-M", str(motif_singletons_dir),
        "-m", str(motif_list), "-o", str(work_dir / database_prefix),
        "--bgpadding", str(bg_padding), "-t", str(threads), "-c", cbust,
    ]

    print(f"[step2] {work_dir} / {database_prefix}", flush=True)
    if dry_run:
        print("DRY-RUN:", " ".join(cmd_fasta))
        print("DRY-RUN:", " ".join(cmd_db))
        return
    subprocess.run(cmd_fasta, check=True)
    subprocess.run(cmd_db, check=True)


def discover_run_roots(base: Path, celltypes: list[str] | None) -> list[Path]:
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


def process_run_root(
    run_root: Path, *, layout: str, chromosomes: list[str] | None,
    genome_fasta: Path, chromsizes: Path, cistarget_scripts_dir: Path,
    motif_singletons_dir: Path, threads: int, bg_padding: int,
    dry_run: bool, cbust_path: Path | None,
) -> None:
    print(f"\n=== step2: {run_root} ({layout}) ===", flush=True)
    if layout == "all-peaks":
        run_cistarget_for_bed(
            region_bed=run_root / "all_peaks" / "all_peaks_cistopic_regions.bed",
            work_dir=run_root / "all_peaks",
            database_prefix="1kb_bg_all_peaks",
            genome_fasta=genome_fasta, chromsizes=chromsizes,
            cistarget_scripts_dir=cistarget_scripts_dir,
            motif_singletons_dir=motif_singletons_dir,
            threads=threads, bg_padding=bg_padding, dry_run=dry_run, cbust_path=cbust_path,
        )
        return

    for chrom in (chromosomes or [f"chr{i}" for i in range(1, 23)]):
        c = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
        run_cistarget_for_bed(
            region_bed=run_root / c / f"{c}_cistopic_regions.bed",
            work_dir=run_root / c,
            database_prefix=f"1kb_bg_{c}",
            genome_fasta=genome_fasta, chromsizes=chromsizes,
            cistarget_scripts_dir=cistarget_scripts_dir,
            motif_singletons_dir=motif_singletons_dir,
            threads=threads, bg_padding=bg_padding, dry_run=dry_run, cbust_path=cbust_path,
        )


def main() -> int:
    p = argparse.ArgumentParser(description="SCENIC+ step 2 per cell type: cisTarget DB")
    p.add_argument("--proj-dir", type=Path, default=Path("."))
    p.add_argument("--out-subdir", default="output_by_celltype")
    p.add_argument("--run-root", type=Path, default=None, help="Single cell-type folder")
    p.add_argument("--celltypes", nargs="*", default=None)
    p.add_argument("--layout", choices=("per-chrom", "all-peaks"), default="all-peaks")
    p.add_argument("--chromosomes", nargs="*", default=None)
    p.add_argument("--threads", "-t", type=int, default=20)
    p.add_argument("--bg-padding", type=int, default=1000)
    p.add_argument("--genome-fasta", type=Path, default=None, help="hg38.fa (default: proj-dir/hg38.fa)")
    p.add_argument("--chromsizes", type=Path, default=None)
    p.add_argument("--cistarget-scripts-dir", type=Path, default=None,
                   help="create_cisTarget_databases folder from SCENIC+ bundle")
    p.add_argument("--motif-singletons-dir", type=Path, default=None)
    p.add_argument("--cbust", type=Path, default=None)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    proj = args.proj_dir.resolve()
    gf = (args.genome_fasta or proj / "hg38.fa").resolve()
    cs = (args.chromsizes or proj / "hg38.chrom.sizes").resolve()
    sd = (args.cistarget_scripts_dir or proj / "create_cisTarget_databases").resolve()
    ms = (args.motif_singletons_dir or proj / "motifs" / "singletons").resolve()

    if args.run_root:
        run_roots = [args.run_root.resolve()]
    else:
        run_roots = discover_run_roots((proj / args.out_subdir).resolve(), args.celltypes)
        if not run_roots:
            print("No cell-type folders found.", flush=True)
            return 1

    for rr in run_roots:
        process_run_root(
            rr, layout=args.layout, chromosomes=args.chromosomes,
            genome_fasta=gf, chromsizes=cs, cistarget_scripts_dir=sd,
            motif_singletons_dir=ms, threads=args.threads, bg_padding=args.bg_padding,
            dry_run=args.dry_run, cbust_path=args.cbust,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
