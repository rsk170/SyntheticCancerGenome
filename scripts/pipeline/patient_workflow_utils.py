#!/usr/bin/env python3
"""Prepare patient-level inputs for SyntheticCancerGenome runs.

This script does two jobs:

1. Build a patient manifest from clone_proportions.csv and clone_*.maf files.
2. Convert clone MAF variants from hg19 coordinates to hg38 VCF records through
   a VCF-aware liftover tool, then write SyntheticCancerGenome mutation strings.

Reliable conversion requires --hg19-fasta so MAF insertions/deletions can be
represented as valid anchored VCF records and normalized before liftover. Runs
without --hg19-fasta stop by default; pass --drop-indels-without-source-fasta
only for explicit SNP/MNV-only test runs where indels may be dropped.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable
from urllib.parse import quote


DEFAULT_CHAIN = Path.home() / ".rbbt/share/lift_over/hg19ToHg38.over.chain.gz"
DEFAULT_HG38_FASTA = Path.home() / ".rbbt/share/organisms/Hsa/hg38/hg38.fa.gz"
SUPPORTED_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
REQUIRED_MAF_COLUMNS = {
    "Chromosome",
    "Start_position",
    "End_position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Variant_Type",
}
KNOWN_INFO_DEFINITIONS = {
    "MAF_ROW": 'Number=1,Type=String,Description="Original MAF row or derived row id"',
    "CLONE_ID": 'Number=1,Type=String,Description="Clone identifier from clone_proportions.csv"',
    "GENE": 'Number=1,Type=String,Description="MAF Hugo_Symbol"',
    "VARIANT_TYPE": 'Number=1,Type=String,Description="MAF Variant_Type"',
    "VARIANT_CLASSIFICATION": 'Number=1,Type=String,Description="MAF Variant_Classification"',
}


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    ident: str
    ref: str
    alt: str
    info: dict[str, str]


@dataclass(frozen=True)
class DroppedRecord:
    clone_id: str
    maf_file: str
    maf_row: str
    reason: str
    detail: str = ""


class VariantDrop(Exception):
    def __init__(self, reason: str, detail: str = ""):
        super().__init__(detail)
        self.reason = reason
        self.detail = detail


def run_command(cmd: list[str], *, stdout=None, stderr=None) -> None:
    try:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
    except FileNotFoundError as exc:
        raise SystemExit(f"Required command not found: {cmd[0]}") from exc
    except subprocess.CalledProcessError as exc:
        pretty = " ".join(cmd)
        raise SystemExit(f"Command failed with exit code {exc.returncode}: {pretty}") from exc


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r", newline="")


def relpath(path: Path, start: Path) -> str:
    try:
        return str(path.resolve().relative_to(start.resolve()))
    except ValueError:
        return str(path.resolve())


def natural_chrom_key(chrom: str) -> tuple[int, str]:
    c = chrom.removeprefix("chr")
    if c.isdigit():
        return (int(c), "")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (order.get(c, 999), c)


def normalize_chrom(raw: str) -> str | None:
    chrom = raw.strip()
    if not chrom:
        return None
    chrom = chrom.removeprefix("chr").removeprefix("CHR")
    chrom = "M" if chrom in {"MT", "Mt", "m", "mt"} else chrom.upper()
    if chrom not in SUPPORTED_CHROMS:
        return None
    if chrom == "MT":
        chrom = "M"
    return f"chr{chrom}"


def clean_allele(value: str) -> str:
    value = (value or "").strip().upper()
    if value in {"", ".", "NA", "N/A", "NULL"}:
        return ""
    return value


def choose_tumor_alt(row: dict[str, str]) -> str:
    ref = clean_allele(row.get("Reference_Allele", ""))
    allele2 = clean_allele(row.get("Tumor_Seq_Allele2", ""))
    allele1 = clean_allele(row.get("Tumor_Seq_Allele1", ""))
    if allele2 and allele2 != ref:
        return allele2
    if allele1 and allele1 != ref:
        return allele1
    return allele2 or allele1


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def valid_dna(seq: str) -> bool:
    return bool(re.fullmatch(r"[ACGTN]+", seq))


def info_value(value: str | int | float | None) -> str:
    if value is None:
        return "."
    text = str(value)
    return quote(text, safe="._:-")


class FastaReference:
    def __init__(self, fasta: Path, cache_dir: Path):
        self.original = fasta.resolve()
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.path = self.cache_dir / self.original.name
        if self.path.exists() or self.path.is_symlink():
            if self.path.resolve() != self.original:
                self.path.unlink()
        if not self.path.exists():
            self.path.symlink_to(self.original)
        self._ensure_index()
        self.lengths = self._read_fai()

    def _ensure_index(self) -> None:
        fai = Path(str(self.path) + ".fai")
        if not fai.exists():
            result = subprocess.run(
                ["samtools", "faidx", str(self.path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0:
                return
            if self.path.suffix == ".gz":
                self.path = self._decompress_gzip_fasta()
                fai = Path(str(self.path) + ".fai")
                if not fai.exists():
                    run_command(["samtools", "faidx", str(self.path)], stdout=subprocess.DEVNULL)
                return
            sys.stderr.write(result.stderr)
            raise SystemExit(f"Could not index FASTA with samtools faidx: {self.path}")

    def _decompress_gzip_fasta(self) -> Path:
        target = self.cache_dir / self.path.name.removesuffix(".gz")
        if target.exists() and target.stat().st_size > 0:
            return target

        tmp = target.with_name(target.name + ".tmp")
        print(f"Decompressing non-BGZF FASTA for indexing: {self.path} -> {target}")
        with gzip.open(self.path, "rb") as source, tmp.open("wb") as sink:
            shutil.copyfileobj(source, sink, length=1024 * 1024)
        tmp.replace(target)
        return target

    def _read_fai(self) -> dict[str, int]:
        lengths: dict[str, int] = {}
        with Path(str(self.path) + ".fai").open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    lengths[fields[0]] = int(fields[1])
        return lengths

    def contig(self, chrom: str) -> str | None:
        candidates = [chrom]
        if chrom.startswith("chr"):
            bare = chrom.removeprefix("chr")
            candidates.extend([bare, "MT" if bare == "M" else bare])
        else:
            candidates.append(f"chr{chrom}")
        for candidate in candidates:
            if candidate in self.lengths:
                return candidate
        return None

    def fetch(self, chrom: str, start: int, end: int) -> str:
        contig = self.contig(chrom)
        if contig is None:
            raise ValueError(f"contig not present in FASTA: {chrom}")
        if start < 1 or end < start:
            raise ValueError(f"invalid FASTA range: {chrom}:{start}-{end}")
        region = f"{contig}:{start}-{end}"
        result = subprocess.run(
            ["samtools", "faidx", str(self.path), region],
            check=True,
            text=True,
            capture_output=True,
        )
        return "".join(
            line.strip() for line in result.stdout.splitlines() if not line.startswith(">")
        ).upper()


def find_clone_mafs(patient_dir: Path) -> dict[str, Path]:
    mafs: dict[str, Path] = {}
    for maf in sorted(patient_dir.glob("clone_*.maf")):
        match = re.match(r"(clone_\d+)(?:_|\.maf)", maf.name)
        if match:
            mafs[match.group(1)] = maf
    return mafs


def timepoint_vaf_columns(fieldnames: Iterable[str]) -> list[str]:
    cols = [name for name in fieldnames if re.fullmatch(r"t\d+_vaf_pct", name)]
    return sorted(cols, key=lambda name: int(name[1:].split("_", 1)[0]))


def choose_timepoint(columns: list[str], requested: str) -> str:
    if not columns:
        raise SystemExit("clone_proportions.csv has no tN_vaf_pct columns")
    if requested == "latest":
        return columns[-1].split("_", 1)[0]
    requested = requested.lower()
    column = f"{requested}_vaf_pct"
    if column not in columns:
        available = ", ".join(col.split("_", 1)[0] for col in columns)
        raise SystemExit(f"Unknown timepoint {requested!r}. Available: {available}")
    return requested


def selected_timepoint_for_patient(patient_dir: Path, requested: str) -> str:
    proportions = patient_dir / "clone_proportions.csv"
    if not proportions.exists():
        raise SystemExit(f"Missing clone proportions file: {proportions}")
    with proportions.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"Empty clone proportions file: {proportions}")
        return choose_timepoint(timepoint_vaf_columns(reader.fieldnames), requested)


def build_manifest(
    *,
    patient_dir: Path,
    repo_root: Path,
    out_dir: Path,
    sex: str,
    normal_depth: float,
    tumor_depth: float,
    timepoint: str,
    include_zero_clones: bool,
) -> tuple[list[dict[str, str]], str]:
    proportions = patient_dir / "clone_proportions.csv"
    germline = patient_dir / "diploid_genotype_germline_hg38.out"
    if not proportions.exists():
        raise SystemExit(f"Missing clone proportions file: {proportions}")
    if not germline.exists():
        raise SystemExit(f"Missing existing germline file: {germline}")

    mafs = find_clone_mafs(patient_dir)
    with proportions.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"Empty clone proportions file: {proportions}")
        vaf_columns = timepoint_vaf_columns(reader.fieldnames)
        selected_timepoint = choose_timepoint(vaf_columns, timepoint)
        selected_column = f"{selected_timepoint}_vaf_pct"
        rows = list(reader)

    total_pct = sum(float(row[selected_column] or 0) for row in rows)
    if total_pct <= 0:
        raise SystemExit(f"Clone fractions at {selected_timepoint} sum to zero")

    manifest_rows: list[dict[str, str]] = []
    for row in rows:
        clone_id = row["clone_id"]
        maf = mafs.get(clone_id)
        if maf is None:
            raise SystemExit(f"Missing MAF for clone_id {clone_id} in {patient_dir}")
        pct = float(row[selected_column] or 0)
        if pct == 0 and not include_zero_clones:
            continue
        clone_fraction = pct / total_pct
        clone_vcf = out_dir / "hg38_vcf" / f"{clone_id}.hg38.normalized.vcf"
        clone_mutations = out_dir / "hg38_mutations" / f"{clone_id}.hg38.scg_mutations"
        manifest_rows.append(
            {
                "patient_id": row["patient_id"],
                "sex": sex,
                "normal_target_depth": f"{normal_depth:g}",
                "tumor_target_depth": f"{tumor_depth:g}",
                "clone_id": clone_id,
                "clone_type": row.get("clone_type", ""),
                "parent_clone_id": row.get("parent_clone_id", ""),
                "timepoint": selected_timepoint,
                "source_vaf_pct": f"{pct:g}",
                "clone_fraction": f"{clone_fraction:.10g}",
                "maf_path": relpath(maf, repo_root),
                "germline_path": relpath(germline, repo_root),
                "source_build": "hg19",
                "simulation_build": "hg38",
                "clone_hg38_vcf_path": relpath(clone_vcf, repo_root),
                "clone_hg38_mutations_path": relpath(clone_mutations, repo_root),
                "simulate_depth": f"{tumor_depth * clone_fraction:.10g}",
            }
        )

    fraction_sum = sum(float(row["clone_fraction"]) for row in manifest_rows)
    if abs(fraction_sum - 1.0) > 1e-6:
        raise SystemExit(f"Manifest clone fractions sum to {fraction_sum}, not 1")
    return manifest_rows, selected_timepoint


def write_manifest(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "patient_id",
        "sex",
        "normal_target_depth",
        "tumor_target_depth",
        "clone_id",
        "clone_type",
        "parent_clone_id",
        "timepoint",
        "source_vaf_pct",
        "clone_fraction",
        "maf_path",
        "germline_path",
        "source_build",
        "simulation_build",
        "clone_hg38_vcf_path",
        "clone_hg38_mutations_path",
        "simulate_depth",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def maf_to_vcf_records(
    maf_path: Path,
    clone_id: str,
    source_fasta: FastaReference | None,
) -> tuple[list[VariantRecord], list[DroppedRecord]]:
    records: list[VariantRecord] = []
    dropped: list[DroppedRecord] = []
    with open_text(maf_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Empty MAF: {maf_path}")
        missing = REQUIRED_MAF_COLUMNS - set(reader.fieldnames)
        if missing:
            raise SystemExit(f"{maf_path} is missing MAF columns: {', '.join(sorted(missing))}")

        for row_number, row in enumerate(reader, start=2):
            if not row or row.get("Chromosome") == "Chromosome":
                continue
            row_id = f"{clone_id}_row{row_number}"
            chrom = normalize_chrom(row["Chromosome"])
            if chrom is None:
                dropped.append(
                    DroppedRecord(
                        clone_id,
                        maf_path.name,
                        str(row_number),
                        "unsupported_contig",
                        row.get("Chromosome", ""),
                    )
                )
                continue

            try:
                start = int(row["Start_position"])
                end = int(row["End_position"])
            except ValueError:
                dropped.append(
                    DroppedRecord(clone_id, maf_path.name, str(row_number), "invalid_coordinates")
                )
                continue

            ref = clean_allele(row["Reference_Allele"])
            alt = choose_tumor_alt(row)
            if not ref or not alt or ref == alt:
                dropped.append(
                    DroppedRecord(clone_id, maf_path.name, str(row_number), "missing_or_same_alt")
                )
                continue

            if "," in alt:
                alts = [clean_allele(part) for part in alt.split(",") if clean_allele(part)]
            else:
                alts = [alt]

            for alt_index, one_alt in enumerate(alts, start=1):
                ident = row_id if len(alts) == 1 else f"{row_id}_alt{alt_index}"
                try:
                    record = maf_row_to_record(row, ident, chrom, start, end, ref, one_alt, source_fasta)
                except VariantDrop as exc:
                    dropped.append(
                        DroppedRecord(clone_id, maf_path.name, str(row_number), exc.reason, exc.detail)
                    )
                    continue
                except ValueError as exc:
                    dropped.append(
                        DroppedRecord(clone_id, maf_path.name, str(row_number), "invalid_variant", str(exc))
                    )
                    continue
                if record is None:
                    dropped.append(
                        DroppedRecord(
                            clone_id,
                            maf_path.name,
                            str(row_number),
                            "requires_hg19_fasta_for_indel",
                            f"REF={ref};ALT={one_alt}",
                        )
                    )
                    continue
                records.append(record)
    return records, dropped


def maf_row_to_record(
    row: dict[str, str],
    ident: str,
    chrom: str,
    start: int,
    end: int,
    ref: str,
    alt: str,
    source_fasta: FastaReference | None,
) -> VariantRecord | None:
    if ref == "-" or alt == "-":
        if source_fasta is None:
            return None
        if ref == "-" and alt != "-":
            if end == start + 1:
                anchor_pos = start
            elif start == end + 1:
                anchor_pos = end
            elif start == end:
                raise VariantDrop(
                    "ambiguous_insertion_coordinates",
                    f"start={start};end={end};variant_type={row.get('Variant_Type', '')}",
                )
            else:
                raise VariantDrop(
                    "ambiguous_insertion_coordinates",
                    f"start={start};end={end};variant_type={row.get('Variant_Type', '')}",
                )
            anchor = source_fasta.fetch(chrom, anchor_pos, anchor_pos)
            pos = anchor_pos
            vcf_ref = anchor
            vcf_alt = anchor + alt
        elif ref != "-" and alt == "-":
            anchor_pos = start - 1
            if anchor_pos < 1:
                raise ValueError("deletion at first base cannot be represented with left anchor")
            anchor = source_fasta.fetch(chrom, anchor_pos, anchor_pos)
            pos = anchor_pos
            vcf_ref = anchor + ref
            vcf_alt = anchor
        else:
            raise ValueError("REF and ALT cannot both be '-'")
    else:
        pos = start
        vcf_ref = ref
        vcf_alt = alt

    if source_fasta is not None:
        source_ref = source_fasta.fetch(chrom, pos, pos + len(vcf_ref) - 1)
        if source_ref != vcf_ref:
            raise ValueError(f"source REF mismatch: VCF={vcf_ref}, FASTA={source_ref}")

    if not valid_dna(vcf_ref) or not valid_dna(vcf_alt):
        raise ValueError(f"unsupported allele characters: REF={vcf_ref};ALT={vcf_alt}")

    info = {
        "MAF_ROW": ident.rsplit("row", 1)[-1],
        "GENE": row.get("Hugo_Symbol", ""),
        "VARIANT_TYPE": row.get("Variant_Type", ""),
        "VARIANT_CLASSIFICATION": row.get("Variant_Classification", ""),
    }
    return VariantRecord(chrom, pos, ident, vcf_ref, vcf_alt, info)


def write_vcf(records: list[VariantRecord], path: Path, lengths: dict[str, int] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    records = sorted(records, key=lambda r: (natural_chrom_key(r.chrom), r.pos, r.ident))
    info_keys = sorted({key for record in records for key in record.info})
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=patient_workflow_utils.py\n")
        for key in info_keys:
            definition = KNOWN_INFO_DEFINITIONS.get(
                key,
                f'Number=1,Type=String,Description="INFO field preserved from liftover: {key}"',
            )
            handle.write(f"##INFO=<ID={key},{definition}>\n")
        if lengths:
            for chrom, length in sorted(lengths.items(), key=lambda item: natural_chrom_key(item[0])):
                if normalize_chrom(chrom) == chrom:
                    handle.write(f"##contig=<ID={chrom},length={length}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for record in records:
            info = ";".join(f"{key}={info_value(value)}" for key, value in record.info.items())
            handle.write(
                f"{record.chrom}\t{record.pos}\t{record.ident}\t{record.ref}\t{record.alt}\t.\tPASS\t{info}\n"
            )


def read_vcf(path: Path) -> list[VariantRecord]:
    records: list[VariantRecord] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, pos, ident, ref, alt, _qual, _filter, info_text = line.rstrip("\n").split("\t")[:8]
            info: dict[str, str] = {}
            if info_text and info_text != ".":
                for item in info_text.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        info[key] = value
            for one_alt in alt.split(","):
                records.append(VariantRecord(chrom, int(pos), ident, ref, one_alt, info))
    return records


def bcftools_norm(input_vcf: Path, fasta: FastaReference, output_vcf: Path) -> None:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "bcftools",
            "norm",
            "-f",
            str(fasta.path),
            "-m",
            "-any",
            "-Ov",
            "-o",
            str(output_vcf),
            str(input_vcf),
        ]
    )


def default_dict_path(fasta: Path) -> Path:
    name = fasta.name
    if name.endswith(".fa.gz"):
        return fasta.with_name(name.removesuffix(".fa.gz") + ".dict")
    if name.endswith(".fasta.gz"):
        return fasta.with_name(name.removesuffix(".fasta.gz") + ".dict")
    if name.endswith(".fa"):
        return fasta.with_suffix(".dict")
    if name.endswith(".fasta"):
        return fasta.with_suffix(".dict")
    return fasta.with_suffix(fasta.suffix + ".dict")


def ensure_sequence_dictionary(fasta: FastaReference) -> None:
    dict_path = default_dict_path(fasta.path)
    if not dict_path.exists():
        with dict_path.open("w") as handle:
            run_command(["samtools", "dict", str(fasta.path)], stdout=handle)

    # HTSJDK versions differ in which dict filename they search for compressed
    # references. Link the common alternatives to the same dictionary.
    aliases = {
        fasta.path.with_name(fasta.path.name + ".dict"),
        fasta.path.with_suffix(".dict"),
        dict_path,
    }
    if fasta.path.name.endswith(".fa.gz"):
        aliases.add(fasta.path.with_name(fasta.path.name.removesuffix(".gz") + ".dict"))
    for alias in aliases:
        if alias == dict_path or alias.exists():
            continue
        try:
            alias.symlink_to(dict_path.name)
        except FileExistsError:
            pass


def detect_liftover_tool(
    requested_tool: str,
    picard_jar: Path | None,
    picard_command: str | None,
    crossmap_command: str | None,
) -> tuple[str, list[str]]:
    if requested_tool in {"auto", "picard"}:
        if picard_jar:
            if not picard_jar.exists():
                raise SystemExit(f"Picard jar does not exist: {picard_jar}")
            return "picard", ["java", "-jar", str(picard_jar)]
        env_picard = os.environ.get("PICARD_JAR")
        if env_picard:
            return "picard", ["java", "-jar", env_picard]
        if picard_command or shutil.which("picard"):
            return "picard", [picard_command or "picard"]

    if requested_tool in {"auto", "crossmap"}:
        command = crossmap_command or shutil.which("CrossMap") or shutil.which("CrossMap.py")
        if command:
            return "crossmap", [command]

    if requested_tool == "picard":
        raise SystemExit("Picard LiftoverVcf was requested, but no Picard command or --picard-jar was found")
    if requested_tool == "crossmap":
        raise SystemExit("CrossMap vcf was requested, but no CrossMap command was found")
    raise SystemExit(
        "No VCF-aware liftover tool found. Install Picard or CrossMap, or pass --picard-jar /path/to/picard.jar."
    )


def run_vcf_liftover(
    *,
    input_vcf: Path,
    lifted_vcf: Path,
    rejected_vcf: Path,
    chain: Path,
    target_fasta: FastaReference,
    tool_name: str,
    tool_command: list[str],
) -> None:
    lifted_vcf.parent.mkdir(parents=True, exist_ok=True)
    rejected_vcf.parent.mkdir(parents=True, exist_ok=True)

    if tool_name == "picard":
        ensure_sequence_dictionary(target_fasta)
        run_command(
            tool_command
            + [
                "LiftoverVcf",
                f"I={input_vcf}",
                f"O={lifted_vcf}",
                f"CHAIN={chain}",
                f"REJECT={rejected_vcf}",
                f"R={target_fasta.path}",
                "CREATE_INDEX=false",
                "VALIDATION_STRINGENCY=LENIENT",
            ]
        )
        return

    if tool_name == "crossmap":
        run_command(tool_command + ["vcf", str(chain), str(input_vcf), str(target_fasta.path), str(lifted_vcf)])
        candidates = [
            Path(str(lifted_vcf) + ".unmap"),
            lifted_vcf.with_suffix(lifted_vcf.suffix + ".unmap"),
            lifted_vcf.with_name(lifted_vcf.name + ".unmap"),
        ]
        for candidate in candidates:
            if candidate.exists():
                if candidate != rejected_vcf:
                    shutil.copyfile(candidate, rejected_vcf)
                return
        write_vcf([], rejected_vcf, target_fasta.lengths)
        return

    raise SystemExit(f"Unsupported liftover tool: {tool_name}")


def dropped_from_rejected_vcf(path: Path, clone_id: str) -> list[DroppedRecord]:
    dropped: list[DroppedRecord] = []
    if not path.exists():
        return dropped
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            ident = fields[2]
            info_text = fields[7]
            maf_row = ident
            detail = fields[6] if len(fields) > 6 else ""
            if info_text and info_text != ".":
                for item in info_text.split(";"):
                    if item.startswith("MAF_ROW="):
                        maf_row = item.split("=", 1)[1]
                    elif item.startswith("FILTER=") or item.startswith("REASON="):
                        detail = item
            dropped.append(DroppedRecord(clone_id, path.name, maf_row, "vcf_liftover_rejected", detail))
    return dropped


def drop_unsupported_target_contigs(
    records: list[VariantRecord],
) -> tuple[list[VariantRecord], list[DroppedRecord]]:
    kept: list[VariantRecord] = []
    dropped: list[DroppedRecord] = []
    for record in records:
        if normalize_chrom(record.chrom) != record.chrom:
            dropped.append(
                DroppedRecord(
                    record.info.get("CLONE_ID", ""),
                    "",
                    record.info.get("MAF_ROW", record.ident),
                    "unsupported_target_contig",
                    record.chrom,
                )
            )
            continue
        kept.append(record)
    return kept, dropped


def verify_target_reference(
    records: list[VariantRecord],
    target_fasta: FastaReference,
) -> tuple[list[VariantRecord], list[DroppedRecord]]:
    kept: list[VariantRecord] = []
    dropped: list[DroppedRecord] = []
    for record in records:
        if not valid_dna(record.ref) or not valid_dna(record.alt):
            dropped.append(
                DroppedRecord(
                    record.info.get("CLONE_ID", ""),
                    "",
                    record.info.get("MAF_ROW", record.ident),
                    "unsupported_final_alleles",
                    f"REF={record.ref};ALT={record.alt}",
                )
            )
            continue
        try:
            target_ref = target_fasta.fetch(record.chrom, record.pos, record.pos + len(record.ref) - 1)
        except (ValueError, subprocess.CalledProcessError) as exc:
            dropped.append(
                DroppedRecord(
                    record.info.get("CLONE_ID", ""),
                    "",
                    record.info.get("MAF_ROW", record.ident),
                    "target_ref_fetch_failed",
                    str(exc),
                )
            )
            continue
        if target_ref != record.ref:
            dropped.append(
                DroppedRecord(
                    record.info.get("CLONE_ID", ""),
                    "",
                    record.info.get("MAF_ROW", record.ident),
                    "target_ref_mismatch",
                    f"VCF={record.ref};hg38={target_ref};{record.chrom}:{record.pos}",
                )
            )
            continue
        kept.append(record)
    return kept, dropped


def record_to_scg_mutation(record: VariantRecord, copy_strategy: str) -> str | None:
    chrom = record.chrom
    if copy_strategy == "hash":
        digest = hashlib.sha1(f"{chrom}:{record.pos}:{record.ref}:{record.alt}".encode()).hexdigest()
        copy = "copy-1_" if int(digest[-1], 16) % 2 == 0 else "copy-2_"
    else:
        copy = f"{copy_strategy}_"
    scg_chrom = copy + chrom

    if len(record.ref) == 1 and len(record.alt) == 1:
        return f"{scg_chrom}:{record.pos}:{record.alt}"
    if record.alt.startswith(record.ref) and len(record.alt) > len(record.ref):
        inserted = record.alt[len(record.ref) :]
        return f"{scg_chrom}:{record.pos}:+{inserted}"
    if record.ref.startswith(record.alt) and len(record.ref) > len(record.alt):
        deleted = record.ref[len(record.alt) :]
        return f"{scg_chrom}:{record.pos + len(record.alt)}:-{deleted}"
    if len(record.ref) == len(record.alt):
        return f"{scg_chrom}:{record.pos}:{record.alt}"
    return None


def write_scg_mutations(records: list[VariantRecord], path: Path, copy_strategy: str) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w") as handle:
        for record in sorted(records, key=lambda r: (natural_chrom_key(r.chrom), r.pos, r.ident)):
            mutation = record_to_scg_mutation(record, copy_strategy)
            if mutation is None:
                continue
            handle.write(mutation + "\n")
            count += 1
    return count


def write_dropped(path: Path, dropped: list[DroppedRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clone_id", "maf_file", "maf_row", "reason", "detail"],
            delimiter="\t",
        )
        writer.writeheader()
        for record in dropped:
            writer.writerow(record.__dict__)


def prepare_clone_conversion(
    *,
    maf: Path,
    clone_id: str,
    out_dir: Path,
    source_fasta: FastaReference | None,
    target_fasta: FastaReference,
    chain: Path,
    liftover_tool_name: str,
    liftover_tool_command: list[str],
    copy_strategy: str,
) -> dict[str, str | int]:
    raw_records, dropped = maf_to_vcf_records(maf, clone_id, source_fasta)
    raw_records = [
        replace(record, info={**record.info, "CLONE_ID": clone_id}) for record in raw_records
    ]

    source_vcf = out_dir / "hg19_vcf" / f"{clone_id}.hg19.source.vcf"
    write_vcf(raw_records, source_vcf, source_fasta.lengths if source_fasta else None)

    if source_fasta is not None:
        normalized_hg19 = out_dir / "hg19_vcf" / f"{clone_id}.hg19.normalized.vcf"
        bcftools_norm(source_vcf, source_fasta, normalized_hg19)
    else:
        normalized_hg19 = source_vcf

    lifted_vcf = out_dir / "hg38_vcf" / f"{clone_id}.hg38.lifted.vcf"
    rejected_vcf = out_dir / "rejected_variants" / f"{clone_id}.rejected.vcf"
    run_vcf_liftover(
        input_vcf=normalized_hg19,
        lifted_vcf=lifted_vcf,
        rejected_vcf=rejected_vcf,
        chain=chain,
        target_fasta=target_fasta,
        tool_name=liftover_tool_name,
        tool_command=liftover_tool_command,
    )
    dropped.extend(dropped_from_rejected_vcf(rejected_vcf, clone_id))

    lifted_records = read_vcf(lifted_vcf) if lifted_vcf.exists() else []
    primary_records, contig_dropped = drop_unsupported_target_contigs(lifted_records)
    dropped.extend(contig_dropped)

    primary_lifted_vcf = out_dir / "hg38_vcf" / f"{clone_id}.hg38.lifted.primary.vcf"
    normalized_hg38 = out_dir / "hg38_vcf" / f"{clone_id}.hg38.normalized.vcf"
    prevalidated_hg38 = out_dir / "hg38_vcf" / f"{clone_id}.hg38.normalized.prevalidate.vcf"
    write_vcf(primary_records, primary_lifted_vcf, target_fasta.lengths)
    if primary_records:
        bcftools_norm(primary_lifted_vcf, target_fasta, prevalidated_hg38)
        normalized_records = read_vcf(prevalidated_hg38)
    else:
        normalized_records = []

    verified_records, verify_dropped = verify_target_reference(normalized_records, target_fasta)
    dropped.extend(verify_dropped)
    write_vcf(verified_records, normalized_hg38, target_fasta.lengths)

    mutation_path = out_dir / "hg38_mutations" / f"{clone_id}.hg38.scg_mutations"
    mutation_count = write_scg_mutations(verified_records, mutation_path, copy_strategy)

    dropped_path = out_dir / "dropped_variants" / f"{clone_id}.dropped.tsv"
    write_dropped(dropped_path, dropped)

    return {
        "clone_id": clone_id,
        "source_records": len(raw_records),
        "hg38_records": len(verified_records),
        "scg_mutations": mutation_count,
        "dropped": len(dropped),
        "source_vcf": str(source_vcf),
        "hg38_vcf": str(normalized_hg38),
        "rejected_vcf": str(rejected_vcf),
        "scg_mutations_path": str(mutation_path),
        "dropped_path": str(dropped_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare patient manifest and hg19-to-hg38 clone MAF conversion outputs."
    )
    parser.add_argument("patient_dir", type=Path, help="Patient directory containing clone_*.maf files")
    parser.add_argument("--out-dir", type=Path, help="Output directory")
    parser.add_argument("--sex", default="unknown", choices=["unknown", "male", "female"])
    parser.add_argument("--normal-depth", type=float, default=30)
    parser.add_argument("--tumor-depth", type=float, default=100)
    parser.add_argument("--timepoint", default="latest", help="Timepoint such as t0/t1/t2, or latest")
    parser.add_argument(
        "--include-zero-clones",
        action="store_true",
        help="Keep zero-fraction clones in the manifest and conversion outputs",
    )
    parser.add_argument(
        "--hg19-fasta",
        type=Path,
        help=(
            "hg19 FASTA for REF validation, indel anchoring, and hg19 normalization; "
            "required for reliable conversion"
        ),
    )
    parser.add_argument(
        "--drop-indels-without-source-fasta",
        action="store_true",
        help=(
            "Explicitly allow a test run without --hg19-fasta; indels that need "
            "source-reference anchoring are dropped"
        ),
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
    parser.add_argument("--manifest-only", action="store_true", help="Only write the manifest")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd()
    patient_dir = args.patient_dir.resolve()
    if not patient_dir.exists():
        raise SystemExit(f"Patient directory does not exist: {patient_dir}")

    selected_timepoint = selected_timepoint_for_patient(patient_dir, args.timepoint)
    out_dir = (args.out_dir or patient_dir / f"prepared_hg38_{selected_timepoint}").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows, selected_timepoint = build_manifest(
        patient_dir=patient_dir,
        repo_root=repo_root,
        out_dir=out_dir,
        sex=args.sex,
        normal_depth=args.normal_depth,
        tumor_depth=args.tumor_depth,
        timepoint=args.timepoint,
        include_zero_clones=args.include_zero_clones,
    )
    manifest_path = out_dir / "patient_manifest.csv"
    write_manifest(manifest_rows, manifest_path)
    print(f"Wrote manifest: {manifest_path}")
    print(f"Selected clone fraction timepoint: {selected_timepoint}")
    print(f"Active clones in manifest: {', '.join(row['clone_id'] for row in manifest_rows)}")

    if args.sex == "unknown":
        print("Warning: manifest sex is 'unknown'; set --sex male or --sex female before sex-aware simulation.")

    if args.manifest_only:
        return 0

    for required in [args.hg38_fasta, args.chain]:
        if not required.exists():
            raise SystemExit(f"Missing required file: {required}")
    if shutil.which("samtools") is None:
        raise SystemExit("samtools is required")
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools is required")
    if not args.hg19_fasta and not args.drop_indels_without_source_fasta:
        raise SystemExit(
            "Reliable conversion requires --hg19-fasta. "
            "Use --drop-indels-without-source-fasta only for SNP/MNV-only test runs."
        )

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
    print(f"Using VCF-aware liftover tool: {liftover_tool_name}")

    mafs = find_clone_mafs(patient_dir)
    active_clone_ids = [row["clone_id"] for row in manifest_rows]
    summaries = []
    for clone_id in active_clone_ids:
        maf = mafs[clone_id]
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
