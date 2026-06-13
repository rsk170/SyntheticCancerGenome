#!/usr/bin/env python3
"""Validate paired FASTQ read IDs and report clone/sample ID counts.

This script streams R1/R2 together. It does not modify the FASTQ files and does
not keep all read IDs in memory, so it is suitable for very large gzip FASTQs.
For exact duplicate detection on billion-read FASTQs, use an external-sort based
check in the Slurm wrapper instead.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter
from pathlib import Path
from typing import TextIO


CLONE_RE = re.compile(r"_tumor(?:_t[0-9]+)?_clone_[0-9]+")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stream a paired FASTQ and check that R1/R2 read identifiers match "
            "record-by-record after removing /1 and /2 suffixes."
        )
    )
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument(
        "--max-records",
        type=int,
        help="Only inspect this many read pairs; omit for a full-file check.",
    )
    parser.add_argument(
        "--max-examples",
        type=int,
        default=20,
        help="Maximum mismatch/malformed examples to store in the JSON report.",
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def read_record(handle: TextIO) -> tuple[str, str, str, str] | None:
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        if any(lines[1:]):
            raise ValueError("Truncated FASTQ record at EOF")
        return None
    if any(line == "" for line in lines):
        raise ValueError("Truncated FASTQ record")
    return tuple(line.rstrip("\n") for line in lines)  # type: ignore[return-value]


def normalize_header(header: str, mate_suffix: str) -> str:
    token = header.split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith(mate_suffix):
        token = token[: -len(mate_suffix)]
    return token


def clone_label(base_id: str) -> str:
    match = CLONE_RE.search(base_id)
    if match:
        return match.group(0).lstrip("_")
    return "UNKNOWN"


def main() -> int:
    args = parse_args()
    if not args.read1.exists():
        raise SystemExit(f"R1 does not exist: {args.read1}")
    if not args.read2.exists():
        raise SystemExit(f"R2 does not exist: {args.read2}")

    records = 0
    malformed = 0
    mismatches = 0
    examples: list[dict[str, str | int]] = []
    clone_counts: Counter[str] = Counter()

    with open_text(args.read1) as r1, open_text(args.read2) as r2:
        while True:
            rec1 = read_record(r1)
            rec2 = read_record(r2)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                mismatches += 1
                examples.append(
                    {
                        "record": records + 1,
                        "read1": rec1[0] if rec1 else "EOF",
                        "read2": rec2[0] if rec2 else "EOF",
                    }
                )
                break

            records += 1
            h1, seq1, plus1, qual1 = rec1
            h2, seq2, plus2, qual2 = rec2

            if not h1.startswith("@") or not h2.startswith("@") or not plus1.startswith("+") or not plus2.startswith("+"):
                malformed += 1
                if len(examples) < args.max_examples:
                    examples.append({"record": records, "read1": h1, "read2": h2, "problem": "malformed_header_or_plus"})

            if len(seq1) != len(qual1) or len(seq2) != len(qual2):
                malformed += 1
                if len(examples) < args.max_examples:
                    examples.append({"record": records, "read1": h1, "read2": h2, "problem": "sequence_quality_length_mismatch"})

            base1 = normalize_header(h1, "/1")
            base2 = normalize_header(h2, "/2")
            clone_counts[clone_label(base1)] += 1
            if base1 != base2:
                mismatches += 1
                if len(examples) < args.max_examples:
                    examples.append({"record": records, "read1": h1, "read2": h2, "problem": "mate_id_mismatch"})

            if args.max_records and records >= args.max_records:
                break

    report = {
        "read1": str(args.read1),
        "read2": str(args.read2),
        "records_checked": records,
        "full_file": args.max_records is None,
        "mate_id_mismatches": mismatches,
        "malformed_records": malformed,
        "clone_counts": dict(sorted(clone_counts.items())),
        "examples": examples,
        "status": "PASS" if mismatches == 0 and malformed == 0 else "FAIL",
    }

    text = json.dumps(report, indent=2) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(text)
    print(text, end="")

    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
