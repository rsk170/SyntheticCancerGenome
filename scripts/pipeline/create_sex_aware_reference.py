#!/usr/bin/env python3
"""Step 6: create one sex-aware hg38 ploidy reference for a patient."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path

from patient_workflow_utils import DEFAULT_HG38_FASTA


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 6: create a patient-specific hg38 reference with duplicated "
            "copy-1_/copy-2_ contigs and XX/XY-aware sex chromosome ploidy."
        )
    )
    parser.add_argument("manifest", type=Path, help="patient_manifest.csv from Step 2")
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_HG38_FASTA)
    parser.add_argument("--out-dir", type=Path, help="Output directory; defaults to PATIENT_DIR/sex_aware_reference")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the existing sex-aware reference if it already exists",
    )
    return parser.parse_args()


def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def read_manifest(path: Path) -> tuple[str, str, Path]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"Manifest has no rows: {path}")

    patient_ids = {row["patient_id"] for row in rows}
    sexes = {row["sex"] for row in rows}
    if len(patient_ids) != 1:
        raise SystemExit(f"Manifest contains multiple patient IDs: {', '.join(sorted(patient_ids))}")
    if len(sexes) != 1:
        raise SystemExit(f"Manifest contains multiple sex values: {', '.join(sorted(sexes))}")

    sex = next(iter(sexes))
    if sex not in {"male", "female"}:
        raise SystemExit(f"Manifest sex must be male or female, got: {sex}")

    germline_path = Path(rows[0]["germline_path"])
    patient_dir = germline_path.parent if not germline_path.is_absolute() else germline_path.parent
    return next(iter(patient_ids)), sex, patient_dir


def canonical_chrom(contig_name: str) -> str:
    name = contig_name.split()[0]
    if name.startswith("chr"):
        name = name[3:]
    name = name.split("_", 1)[0]
    return "M" if name == "MT" else name


def skip_contig(header: str, sex: str) -> bool:
    chrom = canonical_chrom(header[1:])
    return sex == "female" and chrom == "Y"


def write_record(output, header: str, sequence_parts: list[str], copy_number: int, line_width: int = 80) -> None:
    output.write(f">copy-{copy_number}_{header[1:]}")
    if not header.endswith("\n"):
        output.write("\n")

    # NEAT mutation application is case-sensitive, so do not preserve soft-masked
    # lowercase bases from the source reference.
    sequence = "".join(part.strip() for part in sequence_parts).upper()
    for start in range(0, len(sequence), line_width):
        output.write(sequence[start : start + line_width] + "\n")


def second_copy_allowed(header: str, sex: str) -> bool:
    chrom = canonical_chrom(header[1:])
    if sex == "male":
        return chrom not in {"X", "Y", "M"}
    return chrom != "M"


def iter_fasta_records(path: Path):
    header = None
    sequence_parts: list[str] = []
    with open_fasta(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, sequence_parts
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)
        if header is not None:
            yield header, sequence_parts


def build_reference(source_fasta: Path, output_fasta: Path, sex: str) -> dict[str, int]:
    counts = {
        "copy_1_contigs": 0,
        "copy_2_contigs": 0,
        "skipped_y_contigs": 0,
    }
    with output_fasta.open("w") as output:
        for header, sequence_parts in iter_fasta_records(source_fasta):
            if skip_contig(header, sex):
                counts["skipped_y_contigs"] += 1
                continue
            write_record(output, header, sequence_parts, 1)
            counts["copy_1_contigs"] += 1

        for header, sequence_parts in iter_fasta_records(source_fasta):
            if skip_contig(header, sex):
                continue
            if not second_copy_allowed(header, sex):
                continue
            write_record(output, header, sequence_parts, 2)
            counts["copy_2_contigs"] += 1
    return counts


def run_samtools_faidx(path: Path) -> None:
    try:
        subprocess.run(["samtools", "faidx", str(path)], check=True)
    except FileNotFoundError as exc:
        raise SystemExit("samtools is required to index the sex-aware reference") from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"samtools faidx failed for {path}") from exc


def main() -> int:
    args = parse_args()
    manifest = args.manifest.resolve()
    if not manifest.exists():
        raise SystemExit(f"Manifest does not exist: {manifest}")
    if not args.reference_fasta.exists():
        raise SystemExit(f"Reference FASTA does not exist: {args.reference_fasta}")

    patient_id, sex, patient_dir = read_manifest(manifest)
    out_dir = (args.out_dir or patient_dir / "sex_aware_reference").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    output_fasta = out_dir / f"{patient_id}.hg38.{sex}.ploidy.fa"
    metadata_path = out_dir / f"{patient_id}.hg38.{sex}.ploidy.json"
    if output_fasta.exists() and not args.force:
        print(f"Reference already exists: {output_fasta}")
        print("Use --force to rebuild it.")
        return 0

    print(f"Building {sex} sex-aware reference: {output_fasta}")
    counts = build_reference(args.reference_fasta, output_fasta, sex)
    print(f"Indexing reference: {output_fasta}.fai")
    run_samtools_faidx(output_fasta)

    metadata = {
        "patient_id": patient_id,
        "sex": sex,
        "source_reference": str(args.reference_fasta.resolve()),
        "reference": str(output_fasta),
        "fai": str(output_fasta) + ".fai",
        "sequence_case": "uppercase",
        "copy_rule": {
            "male": "autosomes and non-sex/non-M contigs diploid; X, Y, and M haploid",
            "female": "autosomes and X diploid; M haploid; Y contigs absent",
        }[sex],
        **counts,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
