#!/usr/bin/env python3
"""Count tumor/normal read support for expected truth variants.

This script is intended as a targeted follow-up to somatic variant calling. It
does not call variants. Instead, it checks whether expected truth alleles have
direct read support in existing indexed BAM files.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam


MISSING = object()


@dataclass(frozen=True)
class TruthVariant:
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def pos0(self) -> int:
        return self.pos - 1

    @property
    def key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos, self.ref, self.alt)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count REF/ALT read support for truth variants in tumor and normal "
            "BAMs, and mark which truth variants were recovered as PASS calls."
        )
    )
    parser.add_argument("--truth-vcf", type=Path, required=True)
    parser.add_argument("--called-vcf", type=Path, required=True)
    parser.add_argument("--tumor-bam", type=Path, required=True)
    parser.add_argument("--normal-bam", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--summary-tsv", type=Path, required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=10)
    parser.add_argument(
        "--min-alt-count",
        type=int,
        default=1,
        help="Minimum tumor ALT-supporting reads used for summary support counts.",
    )
    return parser.parse_args()


def check_inputs(paths: list[Path]) -> None:
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            raise SystemExit(f"Required input is missing or empty: {path}")


def variant_type(ref: str, alt: str) -> str:
    if any(token in alt for token in ("<", ">", "*", ",")):
        return "COMPLEX"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    return "COMPLEX"


def load_truth_variants(path: Path) -> list[TruthVariant]:
    variants: list[TruthVariant] = []
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            if not record.alts:
                continue
            for alt in record.alts:
                variants.append(
                    TruthVariant(
                        chrom=record.chrom,
                        pos=record.pos,
                        ref=record.ref.upper(),
                        alt=alt.upper(),
                    )
                )
    return variants


def load_pass_call_keys(path: Path) -> set[tuple[str, int, str, str]]:
    keys: set[tuple[str, int, str, str]] = set()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            filters = set(record.filter.keys())
            if filters and "PASS" not in filters:
                continue
            if not record.alts:
                continue
            for alt in record.alts:
                keys.add((record.chrom, record.pos, record.ref.upper(), alt.upper()))
    return keys


def read_passes(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return (
        not read.is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and not read.is_duplicate
        and read.mapping_quality >= min_mapq
    )


def aligned_pair_maps(read: pysam.AlignedSegment) -> tuple[dict[int, int | None], list[tuple[int | None, int | None]]]:
    pairs = read.get_aligned_pairs(matches_only=False)
    ref_to_qpos: dict[int, int | None] = {}
    for qpos, rpos in pairs:
        if rpos is not None and rpos not in ref_to_qpos:
            ref_to_qpos[rpos] = qpos
    return ref_to_qpos, pairs


def base_quality_ok(read: pysam.AlignedSegment, qpos: int, min_baseq: int) -> bool:
    if read.query_qualities is None:
        return True
    return read.query_qualities[qpos] >= min_baseq


def query_bases_for_ref_span(
    read: pysam.AlignedSegment,
    ref_to_qpos: dict[int, int | None],
    start0: int,
    length: int,
    min_baseq: int,
) -> str | None:
    bases: list[str] = []
    sequence = read.query_sequence
    if sequence is None:
        return None

    for ref_pos in range(start0, start0 + length):
        qpos = ref_to_qpos.get(ref_pos, MISSING)
        if qpos is MISSING or qpos is None:
            return None
        if not base_quality_ok(read, qpos, min_baseq):
            return None
        bases.append(sequence[qpos].upper())
    return "".join(bases)


def inserted_bases_after_anchor(
    read: pysam.AlignedSegment,
    pairs: list[tuple[int | None, int | None]],
    anchor_pos0: int,
    anchor_qpos: int,
    min_baseq: int,
) -> str | None:
    sequence = read.query_sequence
    if sequence is None:
        return None

    anchor_indices = [
        index
        for index, (qpos, rpos) in enumerate(pairs)
        if qpos == anchor_qpos and rpos == anchor_pos0
    ]
    if not anchor_indices:
        return None

    inserted: list[str] = []
    for qpos, rpos in pairs[anchor_indices[-1] + 1 :]:
        if rpos is None and qpos is not None:
            if not base_quality_ok(read, qpos, min_baseq):
                return None
            inserted.append(sequence[qpos].upper())
            continue
        break
    return "".join(inserted)


def classify_read(
    read: pysam.AlignedSegment,
    variant: TruthVariant,
    var_type: str,
    min_baseq: int,
) -> str:
    ref_to_qpos, pairs = aligned_pair_maps(read)
    ref = variant.ref
    alt = variant.alt
    pos0 = variant.pos0
    sequence = read.query_sequence
    if sequence is None:
        return "no_call"

    if var_type in {"SNV", "MNV"}:
        observed = query_bases_for_ref_span(read, ref_to_qpos, pos0, len(ref), min_baseq)
        if observed is None:
            return "no_call"
        if observed == alt:
            return "alt"
        if observed == ref:
            return "ref"
        return "other"

    if var_type == "INS":
        anchor_pos0 = pos0 + len(ref) - 1
        anchor_qpos = ref_to_qpos.get(anchor_pos0, MISSING)
        if anchor_qpos is MISSING or anchor_qpos is None:
            return "no_call"
        if not base_quality_ok(read, anchor_qpos, min_baseq):
            return "no_call"
        if sequence[anchor_qpos].upper() != ref[-1]:
            return "other"

        observed_insert = inserted_bases_after_anchor(read, pairs, anchor_pos0, anchor_qpos, min_baseq)
        if observed_insert is None:
            return "no_call"
        expected_insert = alt[len(ref) :]
        if observed_insert == expected_insert:
            return "alt"
        if observed_insert == "":
            observed_ref = query_bases_for_ref_span(read, ref_to_qpos, pos0, len(ref), min_baseq)
            if observed_ref == ref:
                return "ref"
        return "other"

    if var_type == "DEL":
        anchor_pos0 = pos0 + len(alt) - 1
        anchor_qpos = ref_to_qpos.get(anchor_pos0, MISSING)
        if anchor_qpos is MISSING or anchor_qpos is None:
            return "no_call"
        if not base_quality_ok(read, anchor_qpos, min_baseq):
            return "no_call"

        deleted_positions = range(anchor_pos0 + 1, pos0 + len(ref))
        deleted_qpos = [ref_to_qpos.get(ref_pos, MISSING) for ref_pos in deleted_positions]
        if deleted_qpos and all(qpos is None for qpos in deleted_qpos):
            return "alt"
        if all(qpos is not MISSING and qpos is not None for qpos in deleted_qpos):
            observed_ref = query_bases_for_ref_span(read, ref_to_qpos, pos0, len(ref), min_baseq)
            if observed_ref == ref:
                return "ref"
            if observed_ref is not None:
                return "other"
        return "no_call"

    return "unsupported"


def count_support(
    bam: pysam.AlignmentFile,
    variant: TruthVariant,
    min_mapq: int,
    min_baseq: int,
) -> Counter[str]:
    counts: Counter[str] = Counter({"ref": 0, "alt": 0, "other": 0, "no_call": 0, "unsupported": 0})
    var_type = variant_type(variant.ref, variant.alt)
    if var_type == "COMPLEX":
        counts["unsupported"] = 1
        return counts

    start = max(0, variant.pos0 - 1)
    end = variant.pos0 + max(len(variant.ref), len(variant.alt)) + 1
    try:
        reads = bam.fetch(variant.chrom, start, end)
    except ValueError:
        counts["unsupported"] = 1
        return counts

    seen: set[str] = set()
    for read in reads:
        name = read.query_name
        if name in seen:
            continue
        seen.add(name)
        if not read_passes(read, min_mapq):
            continue
        allele = classify_read(read, variant, var_type, min_baseq)
        counts[allele] += 1
    return counts


def fraction(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.6f}"


def support_label(status: str, tumor_alt: int, normal_alt: int, min_alt_count: int, unsupported: bool) -> str:
    if unsupported:
        return "unsupported_variant_type"
    if status == "called_pass":
        return "called_pass"
    if tumor_alt >= min_alt_count and normal_alt == 0:
        return "truth_only_tumor_alt_support_no_normal_alt"
    if tumor_alt >= min_alt_count and normal_alt > 0:
        return "truth_only_tumor_and_normal_alt_support"
    return "truth_only_no_tumor_alt_support"


def main() -> int:
    args = parse_args()
    check_inputs([args.truth_vcf, args.called_vcf, args.tumor_bam, args.normal_bam])

    truth_variants = load_truth_variants(args.truth_vcf)
    called_keys = load_pass_call_keys(args.called_vcf)

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int]] = []
    summary: Counter[str] = Counter()

    with pysam.AlignmentFile(str(args.tumor_bam), "rb") as tumor_bam, pysam.AlignmentFile(str(args.normal_bam), "rb") as normal_bam:
        for variant in truth_variants:
            var_type = variant_type(variant.ref, variant.alt)
            status = "called_pass" if variant.key in called_keys else "truth_only"
            tumor_counts = count_support(tumor_bam, variant, args.min_mapq, args.min_baseq)
            normal_counts = count_support(normal_bam, variant, args.min_mapq, args.min_baseq)

            unsupported = var_type == "COMPLEX" or tumor_counts["unsupported"] > 0 or normal_counts["unsupported"] > 0
            tumor_informative = tumor_counts["ref"] + tumor_counts["alt"] + tumor_counts["other"]
            normal_informative = normal_counts["ref"] + normal_counts["alt"] + normal_counts["other"]
            label = support_label(status, tumor_counts["alt"], normal_counts["alt"], args.min_alt_count, unsupported)

            summary["truth_total"] += 1
            summary[f"{status}_variants"] += 1
            summary[f"{label}"] += 1
            if tumor_counts["alt"] >= args.min_alt_count:
                summary["truth_variants_with_tumor_alt_support"] += 1
            if normal_counts["alt"] >= args.min_alt_count:
                summary["truth_variants_with_normal_alt_support"] += 1
            if status == "truth_only" and tumor_counts["alt"] >= args.min_alt_count:
                summary["truth_only_with_tumor_alt_support"] += 1
            if status == "truth_only" and tumor_counts["alt"] < args.min_alt_count:
                summary["truth_only_without_tumor_alt_support"] += 1

            rows.append(
                {
                    "chrom": variant.chrom,
                    "pos": variant.pos,
                    "ref": variant.ref,
                    "alt": variant.alt,
                    "variant_type": var_type,
                    "mutect_status": status,
                    "support_label": label,
                    "tumor_ref_count": tumor_counts["ref"],
                    "tumor_alt_count": tumor_counts["alt"],
                    "tumor_other_count": tumor_counts["other"],
                    "tumor_no_call_count": tumor_counts["no_call"],
                    "tumor_informative_depth": tumor_informative,
                    "tumor_alt_fraction": fraction(tumor_counts["alt"], tumor_informative),
                    "normal_ref_count": normal_counts["ref"],
                    "normal_alt_count": normal_counts["alt"],
                    "normal_other_count": normal_counts["other"],
                    "normal_no_call_count": normal_counts["no_call"],
                    "normal_informative_depth": normal_informative,
                    "normal_alt_fraction": fraction(normal_counts["alt"], normal_informative),
                }
            )

    fieldnames = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_type",
        "mutect_status",
        "support_label",
        "tumor_ref_count",
        "tumor_alt_count",
        "tumor_other_count",
        "tumor_no_call_count",
        "tumor_informative_depth",
        "tumor_alt_fraction",
        "normal_ref_count",
        "normal_alt_count",
        "normal_other_count",
        "normal_no_call_count",
        "normal_informative_depth",
        "normal_alt_fraction",
    ]
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary_order = [
        "truth_total",
        "called_pass_variants",
        "truth_only_variants",
        "truth_variants_with_tumor_alt_support",
        "truth_variants_with_normal_alt_support",
        "truth_only_with_tumor_alt_support",
        "truth_only_without_tumor_alt_support",
        "called_pass",
        "truth_only_tumor_alt_support_no_normal_alt",
        "truth_only_tumor_and_normal_alt_support",
        "truth_only_no_tumor_alt_support",
        "unsupported_variant_type",
    ]
    with args.summary_tsv.open("w") as handle:
        handle.write("metric\tcount\n")
        for metric in summary_order:
            handle.write(f"{metric}\t{summary.get(metric, 0)}\n")

    print(f"Wrote read-support table: {args.output_tsv}")
    print(f"Wrote read-support summary: {args.summary_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
