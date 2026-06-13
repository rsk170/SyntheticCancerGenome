#!/usr/bin/env python3
"""Step 4/5: convert manifest MAFs from hg19 to hg38 mutation strings."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

from patient_workflow_utils import (
    DEFAULT_CHAIN,
    DEFAULT_HG38_FASTA,
    FastaReference,
    detect_liftover_tool,
    prepare_clone_conversion,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 4/5: convert clone MAFs listed in a patient manifest from hg19 "
            "to hg38 with VCF-aware liftover, then write SCG mutation strings."
        )
    )
    parser.add_argument("manifest", type=Path, help="patient_manifest.csv from prepare_patient_manifest.py")
    parser.add_argument("--hg19-fasta", type=Path, help="hg19 FASTA; required for reliable conversion")
    parser.add_argument(
        "--drop-indels-without-source-fasta",
        action="store_true",
        help="Explicitly allow a test run without --hg19-fasta; indels that need anchoring are dropped",
    )
    parser.add_argument("--hg38-fasta", type=Path, default=DEFAULT_HG38_FASTA)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument(
        "--liftover-tool",
        default="auto",
        choices=["auto", "picard", "crossmap"],
        help="VCF-aware liftover backend",
    )
    parser.add_argument("--picard-jar", type=Path, help="Path to picard.jar for Picard LiftoverVcf")
    parser.add_argument("--picard-command", help="Picard command wrapper, if available on PATH")
    parser.add_argument("--crossmap-command", help="CrossMap command, if not named CrossMap/CrossMap.py")
    parser.add_argument(
        "--copy-strategy",
        default="copy-1",
        choices=["copy-1", "copy-2", "hash"],
        help="How to assign hg38 mutation-string records to duplicated reference copies",
    )
    return parser.parse_args()


def path_from_manifest(value: str, repo_root: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else repo_root / path


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"Manifest has no active clone rows: {path}")
    required = {"clone_id", "maf_path", "clone_fraction"}
    missing = required - set(rows[0])
    if missing:
        raise SystemExit(f"Manifest missing columns: {', '.join(sorted(missing))}")
    fraction_sum = sum(float(row["clone_fraction"]) for row in rows)
    if abs(fraction_sum - 1.0) > 1e-6:
        raise SystemExit(f"Manifest clone fractions sum to {fraction_sum}, not 1")
    return rows


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd()
    manifest_path = args.manifest.resolve()
    if not manifest_path.exists():
        raise SystemExit(f"Manifest does not exist: {manifest_path}")
    if not args.hg19_fasta and not args.drop_indels_without_source_fasta:
        raise SystemExit(
            "Reliable conversion requires --hg19-fasta. "
            "Use --drop-indels-without-source-fasta only for SNP/MNV-only test runs."
        )
    for required in [args.hg38_fasta, args.chain]:
        if not required.exists():
            raise SystemExit(f"Missing required file: {required}")
    if shutil.which("samtools") is None:
        raise SystemExit("samtools is required")
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools is required")

    out_dir = manifest_path.parent
    ref_cache = out_dir / "ref_cache"
    source_fasta = None
    if args.hg19_fasta:
        if not args.hg19_fasta.exists():
            raise SystemExit(f"Missing hg19 FASTA: {args.hg19_fasta}")
        source_fasta = FastaReference(args.hg19_fasta, ref_cache / "hg19")
    else:
        print("Warning: --hg19-fasta not supplied; indels will be dropped instead of normalized/lifted.")
    target_fasta = FastaReference(args.hg38_fasta, ref_cache / "hg38")
    liftover_tool_name, liftover_tool_command = detect_liftover_tool(
        args.liftover_tool,
        args.picard_jar,
        args.picard_command,
        args.crossmap_command,
    )

    rows = read_manifest(manifest_path)
    print(f"Using VCF-aware liftover tool: {liftover_tool_name}")
    summaries = []
    for row in rows:
        clone_id = row["clone_id"]
        maf = path_from_manifest(row["maf_path"], repo_root)
        if not maf.exists():
            raise SystemExit(f"Missing MAF for {clone_id}: {maf}")
        print(f"Converting {clone_id}: {maf.name}")
        summaries.append(
            prepare_clone_conversion(
                maf=maf,
                clone_id=clone_id,
                out_dir=out_dir,
                source_fasta=source_fasta,
                target_fasta=target_fasta,
                chain=args.chain,
                liftover_tool_name=liftover_tool_name,
                liftover_tool_command=liftover_tool_command,
                copy_strategy=args.copy_strategy,
            )
        )

    summary_path = out_dir / "conversion_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        fieldnames = [
            "clone_id",
            "source_records",
            "hg38_records",
            "scg_mutations",
            "dropped",
            "source_vcf",
            "hg38_vcf",
            "rejected_vcf",
            "scg_mutations_path",
            "dropped_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote conversion summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
