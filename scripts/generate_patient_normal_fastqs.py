#!/usr/bin/env python3
"""Step 8: generate matched normal FASTQs with the repo's normal task."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


RESULT_PREFIX = "STEP8_RESULT_JSON\t"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 8: run SyntheticCancerGenome.normal for one patient using the "
            "Step 6 sex-aware reference and Step 7 sex-aware germline."
        )
    )
    parser.add_argument("manifest", type=Path, help="patient_manifest.csv from Step 2")
    parser.add_argument(
        "--reference-metadata",
        type=Path,
        help="Step 6 reference JSON; defaults to PATIENT_DIR/sex_aware_reference/*.json",
    )
    parser.add_argument(
        "--germline-metadata",
        type=Path,
        help="Step 7 germline JSON; defaults to PATIENT_DIR/sex_aware_germline/*.json",
    )
    parser.add_argument("--out-dir", type=Path, help="Output directory; defaults to PATIENT_DIR/normal_fastq")
    parser.add_argument("--sample-name", help="NEAT sample name; defaults to PATIENT_ID_normal")
    parser.add_argument("--job-name", help="Rbbt job name; defaults to PATIENT_ID_normal")
    parser.add_argument("--depth", type=float, help="Override normal depth from patient_manifest.csv")
    parser.add_argument("--read-length", type=int, help="Override NEAT read length; workflow default is 126")
    parser.add_argument(
        "--neat-cpus",
        type=int,
        help=(
            "Number of chromosome-copy simulations to run in parallel inside "
            "NEATGenReads. This should usually match SLURM_CPUS_PER_TASK."
        ),
    )
    parser.add_argument(
        "--samtools-cpus",
        type=int,
        help="Threads for the final samtools merge step; defaults to the NEAT/Rbbt workflow setting.",
    )
    parser.add_argument(
        "--diploid-reference",
        action="store_true",
        help=(
            "Pass haploid_reference=false. Do not use this for the Step 6 "
            "copy-1/copy-2 sex-aware reference."
        ),
    )
    parser.add_argument("--no-errors", action="store_true", help="Ask NEAT not to simulate sequencing errors")
    parser.add_argument(
        "--no-rename-reads",
        action="store_true",
        help="Disable the repo's default read renaming with position information",
    )
    parser.add_argument(
        "--materialize",
        choices=["symlink", "hardlink", "copy", "none"],
        default="symlink",
        help="How to place FASTQs in --out-dir after the Rbbt job finishes",
    )
    parser.add_argument("--output-prefix", default="normal", help="Output FASTQ prefix in --out-dir")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing output links/files and metadata")
    parser.add_argument("--ruby", default="ruby", help="Ruby executable to use")
    parser.add_argument("--dry-run", action="store_true", help="Print the planned Ruby payload without running NEAT")
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

    required = {"patient_id", "sex", "normal_target_depth", "germline_path"}
    missing = required - set(rows[0])
    if missing:
        raise SystemExit(f"Manifest missing columns: {', '.join(sorted(missing))}")

    patient_ids = {row["patient_id"] for row in rows}
    sexes = {row["sex"] for row in rows}
    normal_depths = {row["normal_target_depth"] for row in rows}
    germlines = {row["germline_path"] for row in rows}

    if len(patient_ids) != 1:
        raise SystemExit(f"Manifest contains multiple patient IDs: {', '.join(sorted(patient_ids))}")
    if len(sexes) != 1:
        raise SystemExit(f"Manifest contains multiple sex values: {', '.join(sorted(sexes))}")
    if len(normal_depths) != 1:
        raise SystemExit(f"Manifest contains multiple normal depths: {', '.join(sorted(normal_depths))}")
    if len(germlines) != 1:
        raise SystemExit(f"Manifest contains multiple germline paths: {', '.join(sorted(germlines))}")

    germline = resolve_path(next(iter(germlines)), repo_root)
    return {
        "patient_id": next(iter(patient_ids)),
        "sex": next(iter(sexes)),
        "normal_target_depth": float(next(iter(normal_depths))),
        "patient_dir": germline.parent,
        "rows": rows,
    }


def default_metadata_path(patient_dir: Path, subdir: str, expected_name: str) -> Path:
    expected = patient_dir / subdir / expected_name
    if expected.exists():
        return expected

    matches = sorted((patient_dir / subdir).glob("*.json"))
    if len(matches) == 1:
        return matches[0]

    raise SystemExit(f"Could not determine {subdir} metadata path. Pass it explicitly.")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"Metadata file does not exist: {path}")
    return json.loads(path.read_text())


def rounded_depth(value: float) -> int:
    rounded = int(math.ceil(value))
    if rounded <= 0:
        raise SystemExit(f"Normal depth must be positive, got: {value}")
    return rounded


def ruby_code() -> str:
    return r"""
require 'json'
require './workflow'

args = JSON.parse(ARGV.fetch(0))
if args["neat_cpus"]
  Scout::Config.set(:cpus, args["neat_cpus"].to_i, :genReads, :NEAT, :gen_reads)
end
if args["samtools_cpus"]
  Scout::Config.set(:cpus, args["samtools_cpus"].to_i, :samtools, :merge)
end
inputs = {
  :germline => args.fetch("germline"),
  :reference => args.fetch("reference"),
  :depth => args.fetch("depth").to_i,
  :sample_name => args.fetch("sample_name"),
  :haploid_reference => args.fetch("haploid_reference"),
  :rename_reads => args.fetch("rename_reads"),
  :no_errors => args.fetch("no_errors"),
  :build => "hg38"
}
inputs[:read_length] = args["read_length"].to_i if args["read_length"]

job = SyntheticCancerGenome.job(:normal, args.fetch("job_name"), inputs)
job.produce
files = job.load
files = files.list if files.respond_to?(:list)
files = Array(files).flatten.compact.collect(&:to_s)

puts "STEP8_RESULT_JSON\t" + JSON.dump({
  "job_path" => job.path.to_s,
  "files" => files
})
"""


def run_normal_job(ruby: str, payload: dict[str, Any]) -> dict[str, Any]:
    command = [ruby, "-Ilib", "-rrbbt-util", "-e", ruby_code(), json.dumps(payload)]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    result_line: str | None = None
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="")
        marker_index = line.find(RESULT_PREFIX)
        if marker_index != -1:
            result_line = line[marker_index + len(RESULT_PREFIX) :].strip()

    returncode = process.wait()
    if returncode != 0:
        raise SystemExit(f"Ruby normal job failed with exit code {returncode}")
    if result_line is None:
        raise SystemExit("Ruby normal job finished but did not report output FASTQs")
    return json.loads(result_line)


def fastq_suffix(path: Path) -> str:
    return ".fq.gz" if path.name.endswith(".fq.gz") else path.suffix


def output_fastq_path(out_dir: Path, prefix: str, read_number: int, source: Path) -> Path:
    return out_dir / f"{prefix}_read{read_number}{fastq_suffix(source)}"


def find_fastq_pair(files: list[str]) -> tuple[Path, Path]:
    read1 = sorted(
        Path(path)
        for path in files
        if path.endswith("_read1.fq.gz") or path.endswith("_read1.fq")
    )
    read2 = sorted(
        Path(path)
        for path in files
        if path.endswith("_read2.fq.gz") or path.endswith("_read2.fq")
    )
    if len(read1) != 1 or len(read2) != 1:
        raise SystemExit(
            "Expected exactly one read1 and one read2 FASTQ from normal job; "
            f"found read1={len(read1)}, read2={len(read2)}"
        )
    for fastq in (read1[0], read2[0]):
        if not fastq.exists() or fastq.stat().st_size == 0:
            raise SystemExit(f"FASTQ output is missing or empty: {fastq}")
    return read1[0], read2[0]


def place_file(source: Path, target: Path, mode: str, overwrite: bool) -> None:
    if mode == "none":
        return
    if target.exists() or target.is_symlink():
        if not overwrite:
            raise SystemExit(f"Output already exists: {target}. Use --overwrite to replace it.")
        target.unlink()

    target.parent.mkdir(parents=True, exist_ok=True)
    if mode == "symlink":
        os.symlink(source, target)
    elif mode == "hardlink":
        os.link(source, target)
    elif mode == "copy":
        shutil.copy2(source, target)
    else:
        raise AssertionError(f"Unknown materialize mode: {mode}")


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

    patient_id = manifest_info["patient_id"]
    sex = manifest_info["sex"]
    patient_dir = manifest_info["patient_dir"]
    sample_name = args.sample_name or f"{patient_id}_normal"
    job_name = args.job_name or sample_name
    requested_depth = args.depth if args.depth is not None else manifest_info["normal_target_depth"]
    depth = rounded_depth(float(requested_depth))

    reference_metadata_path = (
        args.reference_metadata.resolve()
        if args.reference_metadata
        else default_metadata_path(
            patient_dir,
            "sex_aware_reference",
            f"{patient_id}.hg38.{sex}.ploidy.json",
        )
    )
    germline_metadata_path = (
        args.germline_metadata.resolve()
        if args.germline_metadata
        else default_metadata_path(
            patient_dir,
            "sex_aware_germline",
            f"{patient_id}.hg38.{sex}.germline.json",
        )
    )

    reference_metadata = read_json(reference_metadata_path)
    germline_metadata = read_json(germline_metadata_path)
    reference = Path(reference_metadata["reference"]).resolve()
    germline = Path(germline_metadata["output_germline"]).resolve()
    if not reference.exists():
        raise SystemExit(f"Sex-aware reference FASTA does not exist: {reference}")
    if not germline.exists():
        raise SystemExit(f"Sex-aware germline file does not exist: {germline}")

    out_dir = (args.out_dir or patient_dir / "normal_fastq").resolve()
    metadata_path = out_dir / f"{args.output_prefix}.normal_fastq.json"

    payload: dict[str, Any] = {
        "germline": str(germline),
        "reference": str(reference),
        "depth": depth,
        "sample_name": sample_name,
        "job_name": job_name,
        "haploid_reference": not args.diploid_reference,
        "rename_reads": not args.no_rename_reads,
        "no_errors": args.no_errors,
    }
    if args.read_length is not None:
        payload["read_length"] = args.read_length
    if args.neat_cpus is not None:
        if args.neat_cpus <= 0:
            raise SystemExit(f"--neat-cpus must be positive, got: {args.neat_cpus}")
        payload["neat_cpus"] = args.neat_cpus
    if args.samtools_cpus is not None:
        if args.samtools_cpus < 0:
            raise SystemExit(f"--samtools-cpus must be zero or positive, got: {args.samtools_cpus}")
        payload["samtools_cpus"] = args.samtools_cpus

    if args.dry_run:
        print("Dry run. Ruby normal job payload:")
        print(json.dumps(payload, indent=2))
        return 0

    print(f"Generating matched normal FASTQs for patient {patient_id}")
    print(f"Reference: {reference}")
    print(f"Germline: {germline}")
    print(f"Depth: {depth}")
    result = run_normal_job(args.ruby, payload)
    source_read1, source_read2 = find_fastq_pair(result["files"])
    output_read1 = output_fastq_path(out_dir, args.output_prefix, 1, source_read1)
    output_read2 = output_fastq_path(out_dir, args.output_prefix, 2, source_read2)

    place_file(source_read1, output_read1, args.materialize, args.overwrite)
    place_file(source_read2, output_read2, args.materialize, args.overwrite)

    metadata = {
        "patient_id": patient_id,
        "sex": sex,
        "manifest": str(manifest),
        "sample_name": sample_name,
        "job_name": job_name,
        "normal_target_depth": manifest_info["normal_target_depth"],
        "simulated_depth": depth,
        "haploid_reference": not args.diploid_reference,
        "rename_reads": not args.no_rename_reads,
        "no_errors": args.no_errors,
        "read_length": args.read_length,
        "reference_metadata": str(reference_metadata_path),
        "reference": str(reference),
        "germline_metadata": str(germline_metadata_path),
        "germline": str(germline),
        "rbbt_job_path": result["job_path"],
        "source_fastqs": {
            "read1": str(source_read1),
            "read2": str(source_read2),
        },
        "output_fastqs": {
            "read1": str(output_read1) if args.materialize != "none" else None,
            "read2": str(output_read2) if args.materialize != "none" else None,
        },
        "materialize": args.materialize,
        "created_at": datetime.now().astimezone().isoformat(),
    }
    write_metadata(metadata_path, metadata, args.overwrite)

    print(f"Wrote read1: {output_read1 if args.materialize != 'none' else source_read1}")
    print(f"Wrote read2: {output_read2 if args.materialize != 'none' else source_read2}")
    print(f"Wrote metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.exit(1)
