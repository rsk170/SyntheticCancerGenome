#!/usr/bin/env python3
"""Merge per-clone tumor FASTQs into one tumor FASTQ pair."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import BinaryIO, Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 12: merge per-clone tumor FASTQs into one patient-level tumor "
            "FASTQ pair by decompressing each clone stream, concatenating in "
            "manifest order, and recompressing."
        )
    )
    parser.add_argument(
        "manifest",
        type=Path,
        help=(
            "patient_manifest.final_clone_mutations.csv used for tumor clone "
            "generation"
        ),
    )
    parser.add_argument(
        "--clone-fastq-dir",
        type=Path,
        help="Directory from generate_patient_clone_tumor_fastqs.py; defaults to PATIENT_DIR/tumor_clone_fastqs",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory; defaults to PATIENT_DIR/tumor_fastq",
    )
    parser.add_argument("--output-prefix", default="tumor", help="Output FASTQ prefix in --out-dir")
    parser.add_argument("--compresslevel", type=int, default=6, choices=range(1, 10))
    parser.add_argument("--overwrite", action="store_true", help="Replace existing merged FASTQs and metadata")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the clone FASTQs that would be merged without writing outputs",
    )
    return parser.parse_args()


def resolve_path(path: str | Path, repo_root: Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else repo_root / value


def read_manifest(path: Path, repo_root: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"Manifest does not exist: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    if not rows:
        raise SystemExit(f"Manifest has no rows: {path}")

    required = {"patient_id", "sex", "tumor_target_depth", "clone_id", "clone_fraction", "germline_path"}
    missing = required - set(rows[0])
    if missing:
        raise SystemExit(f"Manifest missing columns: {', '.join(sorted(missing))}")

    patient_ids = {row["patient_id"] for row in rows}
    sexes = {row["sex"] for row in rows}
    tumor_depths = {row["tumor_target_depth"] for row in rows}
    germlines = {row["germline_path"] for row in rows}

    if len(patient_ids) != 1:
        raise SystemExit(f"Manifest contains multiple patient IDs: {', '.join(sorted(patient_ids))}")
    if len(sexes) != 1:
        raise SystemExit(f"Manifest contains multiple sex values: {', '.join(sorted(sexes))}")
    if len(tumor_depths) != 1:
        raise SystemExit(f"Manifest contains multiple tumor depths: {', '.join(sorted(tumor_depths))}")
    if len(germlines) != 1:
        raise SystemExit(f"Manifest contains multiple germline paths: {', '.join(sorted(germlines))}")

    germline = resolve_path(next(iter(germlines)), repo_root)
    active_rows = [row for row in rows if float(row["clone_fraction"]) > 0]
    if not active_rows:
        raise SystemExit("Manifest has no active clone rows")

    return {
        "patient_id": next(iter(patient_ids)),
        "sex": next(iter(sexes)),
        "tumor_target_depth": float(next(iter(tumor_depths))),
        "patient_dir": germline.parent,
        "rows": active_rows,
    }


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"Metadata file does not exist: {path}")
    return json.loads(path.read_text())


def existing_path(value: str | None, label: str) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if path.exists() and path.stat().st_size > 0:
        return path
    if path.is_symlink():
        raise SystemExit(f"{label} is a broken symlink: {path}")
    return None


def metadata_fastq_pair(metadata: dict[str, Any]) -> tuple[Path, Path] | None:
    for key in ("output_fastqs", "source_fastqs"):
        fastqs = metadata.get(key) or {}
        read1 = existing_path(fastqs.get("read1"), f"{key}.read1")
        read2 = existing_path(fastqs.get("read2"), f"{key}.read2")
        if read1 and read2:
            return read1, read2
    return None


def discover_clone_fastqs(clone_dir: Path, clone_id: str) -> tuple[Path, Path, Path | None]:
    metadata_matches = sorted(clone_dir.glob("*.tumor_fastq.json"))
    for metadata_path in metadata_matches:
        pair = metadata_fastq_pair(read_json(metadata_path))
        if pair:
            return pair[0], pair[1], metadata_path

    read1_matches = sorted(clone_dir.glob("*_read1.fq.gz")) + sorted(clone_dir.glob("*_read1.fq"))
    read2_matches = sorted(clone_dir.glob("*_read2.fq.gz")) + sorted(clone_dir.glob("*_read2.fq"))
    if len(read1_matches) != 1 or len(read2_matches) != 1:
        raise SystemExit(
            f"Expected exactly one read1/read2 FASTQ for {clone_id} in {clone_dir}; "
            f"found read1={len(read1_matches)}, read2={len(read2_matches)}"
        )
    for fastq in (read1_matches[0], read2_matches[0]):
        if not fastq.exists() or fastq.stat().st_size == 0:
            raise SystemExit(f"Clone FASTQ is missing or empty: {fastq}")
    return read1_matches[0], read2_matches[0], None


def open_input(path: Path) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def merge_streams(sources: list[Path], target: Path, *, compresslevel: int, overwrite: bool) -> dict[str, Any]:
    if target.exists() or target.is_symlink():
        if not overwrite:
            raise SystemExit(f"Output already exists: {target}. Use --overwrite to replace it.")
        target.unlink()

    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + ".tmp")
    if tmp.exists() or tmp.is_symlink():
        if not overwrite:
            raise SystemExit(f"Temporary output already exists: {tmp}. Use --overwrite to replace it.")
        tmp.unlink()

    total_bytes = 0
    total_lines = 0
    source_stats: list[dict[str, Any]] = []
    try:
        with gzip.open(tmp, "wb", compresslevel=compresslevel) as output:
            for source in sources:
                source_bytes = 0
                source_lines = 0
                with open_input(source) as handle:
                    while True:
                        chunk = handle.read(1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)
                        chunk_lines = chunk.count(b"\n")
                        source_bytes += len(chunk)
                        source_lines += chunk_lines
                if source_lines % 4 != 0:
                    raise SystemExit(f"FASTQ line count is not divisible by 4 for {source}: {source_lines}")
                source_stats.append(
                    {
                        "source": str(source),
                        "decompressed_bytes": source_bytes,
                        "decompressed_lines": source_lines,
                        "read_records": source_lines // 4,
                    }
                )
                total_bytes += source_bytes
                total_lines += source_lines
        os.replace(tmp, target)
    except Exception:
        if tmp.exists() or tmp.is_symlink():
            tmp.unlink()
        raise

    return {
        "output": str(target),
        "decompressed_bytes": total_bytes,
        "decompressed_lines": total_lines,
        "read_records": total_lines // 4,
        "sources": source_stats,
    }


def write_metadata(path: Path, metadata: dict[str, Any], overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise SystemExit(f"Metadata already exists: {path}. Use --overwrite to replace it.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2) + "\n")


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd()
    manifest = args.manifest.resolve()
    manifest_info = read_manifest(manifest, repo_root)
    patient_dir = manifest_info["patient_dir"]

    clone_fastq_dir = (args.clone_fastq_dir or patient_dir / "tumor_clone_fastqs").resolve()
    out_dir = (args.out_dir or patient_dir / "tumor_fastq").resolve()
    output_read1 = out_dir / f"{args.output_prefix}_read1.fq.gz"
    output_read2 = out_dir / f"{args.output_prefix}_read2.fq.gz"
    metadata_path = out_dir / f"{args.output_prefix}.tumor_fastq.json"

    clone_entries: list[dict[str, Any]] = []
    read1_sources: list[Path] = []
    read2_sources: list[Path] = []

    for row in manifest_info["rows"]:
        clone_id = row["clone_id"]
        clone_dir = clone_fastq_dir / clone_id
        if not clone_dir.exists():
            raise SystemExit(f"Clone FASTQ directory does not exist for {clone_id}: {clone_dir}")
        read1, read2, clone_metadata = discover_clone_fastqs(clone_dir, clone_id)
        read1_sources.append(read1)
        read2_sources.append(read2)
        clone_entries.append(
            {
                "clone_id": clone_id,
                "clone_fraction": float(row["clone_fraction"]),
                "read1": str(read1),
                "read2": str(read2),
                "metadata": str(clone_metadata) if clone_metadata else None,
            }
        )

    if args.dry_run:
        print("Dry run. Clone FASTQs would be merged in this order:")
        print(json.dumps(clone_entries, indent=2))
        print(f"Planned read1 output: {output_read1}")
        print(f"Planned read2 output: {output_read2}")
        return 0

    print(f"Merging {len(clone_entries)} clone FASTQ pairs for patient {manifest_info['patient_id']}")
    read1_stats = merge_streams(
        read1_sources,
        output_read1,
        compresslevel=args.compresslevel,
        overwrite=args.overwrite,
    )
    read2_stats = merge_streams(
        read2_sources,
        output_read2,
        compresslevel=args.compresslevel,
        overwrite=args.overwrite,
    )

    if read1_stats["read_records"] != read2_stats["read_records"]:
        raise SystemExit(
            "Merged R1/R2 read counts differ: "
            f"read1={read1_stats['read_records']}, read2={read2_stats['read_records']}"
        )

    metadata = {
        "patient_id": manifest_info["patient_id"],
        "sex": manifest_info["sex"],
        "manifest": str(manifest),
        "tumor_target_depth": manifest_info["tumor_target_depth"],
        "clone_fastq_dir": str(clone_fastq_dir),
        "clone_fastqs": clone_entries,
        "output_fastqs": {
            "read1": str(output_read1),
            "read2": str(output_read2),
        },
        "read1_merge": read1_stats,
        "read2_merge": read2_stats,
        "compresslevel": args.compresslevel,
        "created_at": datetime.now().astimezone().isoformat(),
    }
    write_metadata(metadata_path, metadata, args.overwrite)

    print(f"Wrote merged tumor read1: {output_read1}")
    print(f"Wrote merged tumor read2: {output_read2}")
    print(f"Wrote merged tumor metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.exit(1)
