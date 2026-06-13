#!/usr/bin/env python3
"""Step 7: create a sex-aware germline file matching the patient reference."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 7: filter and rename a diploid hg38 germline mutation file so "
            "its contigs match the patient's sex-aware ploidy reference."
        )
    )
    parser.add_argument("manifest", type=Path, help="patient_manifest.csv from Step 2")
    parser.add_argument(
        "--reference-metadata",
        type=Path,
        help="JSON metadata from create_sex_aware_reference.py; defaults to PATIENT_DIR/sex_aware_reference/*.json",
    )
    parser.add_argument("--out-dir", type=Path, help="Output directory; defaults to PATIENT_DIR/sex_aware_germline")
    parser.add_argument("--maternal-copy", default="copy-2", choices=["copy-1", "copy-2"])
    parser.add_argument("--paternal-copy", default="copy-1", choices=["copy-1", "copy-2"])
    parser.add_argument(
        "--strict-diploid",
        action="store_true",
        help="Fail if autosomal germline records do not include both copy-1 and copy-2",
    )
    return parser.parse_args()


def read_manifest(path: Path, repo_root: Path) -> tuple[str, str, Path, Path]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"Manifest has no rows: {path}")

    patient_ids = {row["patient_id"] for row in rows}
    sexes = {row["sex"] for row in rows}
    germlines = {row["germline_path"] for row in rows}
    if len(patient_ids) != 1:
        raise SystemExit(f"Manifest contains multiple patient IDs: {', '.join(sorted(patient_ids))}")
    if len(sexes) != 1:
        raise SystemExit(f"Manifest contains multiple sex values: {', '.join(sorted(sexes))}")
    if len(germlines) != 1:
        raise SystemExit(f"Manifest contains multiple germline paths: {', '.join(sorted(germlines))}")

    sex = next(iter(sexes))
    if sex not in {"male", "female"}:
        raise SystemExit(f"Manifest sex must be male or female, got: {sex}")

    germline = Path(next(iter(germlines)))
    if not germline.is_absolute():
        germline = repo_root / germline
    if not germline.exists():
        raise SystemExit(f"Germline file does not exist: {germline}")
    return next(iter(patient_ids)), sex, germline, germline.parent


def read_reference_metadata(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"Reference metadata does not exist: {path}")
    return json.loads(path.read_text())


def default_reference_metadata(patient_dir: Path, patient_id: str, sex: str) -> Path:
    expected = patient_dir / "sex_aware_reference" / f"{patient_id}.hg38.{sex}.ploidy.json"
    if expected.exists():
        return expected
    matches = sorted((patient_dir / "sex_aware_reference").glob("*.json"))
    if len(matches) == 1:
        return matches[0]
    raise SystemExit(
        "Could not determine reference metadata path. Pass --reference-metadata explicitly."
    )


def read_reference_contigs(metadata: dict) -> set[str]:
    fai = Path(metadata["fai"])
    if not fai.exists():
        raise SystemExit(f"Reference index does not exist: {fai}")
    contigs: set[str] = set()
    with fai.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields:
                contigs.add(fields[0])
    return contigs


def split_copy_contig(contig: str) -> tuple[str | None, str]:
    if contig.startswith("copy-1_"):
        return "copy-1", contig[len("copy-1_") :]
    if contig.startswith("copy-2_"):
        return "copy-2", contig[len("copy-2_") :]
    return None, contig


def canonical_chrom(base_contig: str) -> str:
    name = base_contig.split()[0]
    if name.startswith("chr"):
        name = name[3:]
    name = name.split("_", 1)[0]
    return "M" if name == "MT" else name


def is_autosome(chrom: str) -> bool:
    return chrom.isdigit() and 1 <= int(chrom) <= 22


def transform_contig(
    contig: str,
    sex: str,
    maternal_copy: str,
    paternal_copy: str,
) -> tuple[str | None, str, str]:
    copy, base = split_copy_contig(contig)
    chrom = canonical_chrom(base)

    if copy is None:
        return None, "missing_copy_prefix", contig

    if sex == "female":
        if chrom == "Y":
            return None, "drop_female_y", contig
        if chrom == "M":
            if copy != maternal_copy:
                return None, "drop_nonmaternal_m", contig
            return "copy-1_" + base, "rename_maternal_m_to_haploid", contig
        return contig, "keep", contig

    if chrom == "X":
        if copy != maternal_copy:
            return None, "drop_nonmaternal_x", contig
        return "copy-1_" + base, "rename_maternal_x_to_haploid", contig
    if chrom == "Y":
        if copy != paternal_copy:
            return None, "drop_nonpaternal_y", contig
        return "copy-1_" + base, "keep_paternal_y", contig
    if chrom == "M":
        if copy != maternal_copy:
            return None, "drop_nonmaternal_m", contig
        return "copy-1_" + base, "rename_maternal_m_to_haploid", contig
    return contig, "keep", contig


def check_autosomal_copies(input_copy_by_chrom: dict[str, set[str]]) -> list[str]:
    warnings: list[str] = []
    for chrom in sorted(input_copy_by_chrom, key=lambda c: int(c) if c.isdigit() else 999):
        if not is_autosome(chrom):
            continue
        copies = input_copy_by_chrom[chrom]
        if copies != {"copy-1", "copy-2"}:
            warnings.append(f"autosome chr{chrom} has input copies {sorted(copies)}, expected ['copy-1', 'copy-2']")
    return warnings


def create_germline(
    input_path: Path,
    output_path: Path,
    dropped_path: Path,
    reference_contigs: set[str],
    sex: str,
    maternal_copy: str,
    paternal_copy: str,
) -> tuple[dict, list[str]]:
    stats = Counter()
    action_counts = Counter()
    input_copy_by_chrom: dict[str, set[str]] = {}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with input_path.open() as source, output_path.open("w") as output, dropped_path.open("w") as dropped:
        dropped.write("line_number\tcontig\tpos\talt\treason\tdetail\n")
        for line_number, line in enumerate(source, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            stats["input_records"] += 1
            parts = line.split(":", 2)
            if len(parts) != 3:
                stats["dropped_records"] += 1
                dropped.write(f"{line_number}\t.\t.\t.\tinvalid_record\t{line}\n")
                continue
            contig, pos, alt = parts
            copy, base = split_copy_contig(contig)
            chrom = canonical_chrom(base)
            if copy:
                input_copy_by_chrom.setdefault(chrom, set()).add(copy)

            target_contig, action, detail = transform_contig(contig, sex, maternal_copy, paternal_copy)
            action_counts[action] += 1
            if target_contig is None:
                stats["dropped_records"] += 1
                dropped.write(f"{line_number}\t{contig}\t{pos}\t{alt}\t{action}\t{detail}\n")
                continue
            if target_contig not in reference_contigs:
                stats["dropped_records"] += 1
                action_counts["missing_from_reference"] += 1
                dropped.write(
                    f"{line_number}\t{contig}\t{pos}\t{alt}\tmissing_from_reference\t{target_contig}\n"
                )
                continue
            output.write(f"{target_contig}:{pos}:{alt}\n")
            stats["output_records"] += 1
            if target_contig != contig:
                stats["renamed_records"] += 1

    warnings = check_autosomal_copies(input_copy_by_chrom)
    stats.update({f"action_{key}": value for key, value in action_counts.items()})
    return dict(stats), warnings


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd()
    manifest = args.manifest.resolve()
    if not manifest.exists():
        raise SystemExit(f"Manifest does not exist: {manifest}")

    if args.maternal_copy == args.paternal_copy:
        raise SystemExit("--maternal-copy and --paternal-copy must be different")

    patient_id, sex, germline, patient_dir = read_manifest(manifest, repo_root)
    metadata_path = args.reference_metadata or default_reference_metadata(patient_dir, patient_id, sex)
    reference_metadata = read_reference_metadata(metadata_path)
    reference_contigs = read_reference_contigs(reference_metadata)

    out_dir = (args.out_dir or patient_dir / "sex_aware_germline").resolve()
    output_path = out_dir / f"{patient_id}.hg38.{sex}.germline.scg_mutations"
    dropped_path = out_dir / f"{patient_id}.hg38.{sex}.germline.dropped.tsv"
    metadata_out = out_dir / f"{patient_id}.hg38.{sex}.germline.json"

    print(f"Creating {sex} sex-aware germline: {output_path}")
    stats, warnings = create_germline(
        germline,
        output_path,
        dropped_path,
        reference_contigs,
        sex,
        args.maternal_copy,
        args.paternal_copy,
    )
    if args.strict_diploid and warnings:
        raise SystemExit("Strict diploid check failed:\n" + "\n".join(warnings))

    metadata = {
        "patient_id": patient_id,
        "sex": sex,
        "input_germline": str(germline),
        "reference_metadata": str(metadata_path),
        "reference": reference_metadata.get("reference"),
        "output_germline": str(output_path),
        "dropped_records_path": str(dropped_path),
        "maternal_copy": args.maternal_copy,
        "paternal_copy": args.paternal_copy,
        "copy_assumption": "SyntheticCancerGenome germline task currently writes copy-1 from the father dependency and copy-2 from the mother dependency.",
        "warnings": warnings,
        **stats,
    }
    metadata_out.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote dropped-record report: {dropped_path}")
    print(f"Wrote metadata: {metadata_out}")
    if warnings:
        print("Warnings:")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  - ... {len(warnings) - 10} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
