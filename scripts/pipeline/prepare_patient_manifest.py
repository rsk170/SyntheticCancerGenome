#!/usr/bin/env python3
"""Step 2: build a patient manifest for one simulation timepoint."""

from __future__ import annotations

import argparse
from pathlib import Path

from patient_workflow_utils import (
    build_manifest,
    selected_timepoint_for_patient,
    write_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Step 2: prepare the patient manifest and validate clone fractions."
    )
    parser.add_argument("patient_dir", type=Path, help="Patient directory containing clone_proportions.csv and clone_*.maf files")
    parser.add_argument("--out-dir", type=Path, help="Output directory; defaults to prepared_hg38_<timepoint>")
    parser.add_argument("--sex", required=True, choices=["male", "female"], help="Patient sex for later sex-aware reference/germline steps")
    parser.add_argument("--normal-depth", type=float, default=30)
    parser.add_argument("--tumor-depth", type=float, default=100)
    parser.add_argument("--timepoint", required=True, help="Timepoint to simulate, such as t0, t1, or t2")
    parser.add_argument(
        "--include-zero-clones",
        action="store_true",
        help="Keep zero-fraction clones in the manifest; by default they are skipped",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path.cwd()
    patient_dir = args.patient_dir.resolve()
    if not patient_dir.exists():
        raise SystemExit(f"Patient directory does not exist: {patient_dir}")

    selected_timepoint = selected_timepoint_for_patient(patient_dir, args.timepoint)
    out_dir = (args.out_dir or patient_dir / f"prepared_hg38_{selected_timepoint}").resolve()

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

    fraction_sum = sum(float(row["clone_fraction"]) for row in manifest_rows)
    print(f"Wrote manifest: {manifest_path}")
    print(f"Selected timepoint: {selected_timepoint}")
    print(f"Active clones: {', '.join(row['clone_id'] for row in manifest_rows)}")
    print(f"Clone fraction sum: {fraction_sum:.10f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
