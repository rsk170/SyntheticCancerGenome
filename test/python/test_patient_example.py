#!/usr/bin/env python3
"""Lightweight checks for the bundled patient-workflow example."""

from __future__ import annotations

import csv
import hashlib
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE_DIR = REPO_ROOT / "scripts" / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

from patient_workflow_utils import REQUIRED_MAF_COLUMNS, build_manifest  # noqa: E402


PATIENT_ID = "79ce1d89-46d2-5513-c704-212aa1ed97d2"
EXAMPLE_DIR = REPO_ROOT / "examples" / "patient_workflow" / PATIENT_ID
EXPECTED_ACTIVE_CLONES = {
    "t0": {"clone_1", "clone_3"},
    "t1": {"clone_1", "clone_2", "clone_3"},
    "t2": {"clone_1", "clone_2", "clone_4"},
}


class PatientExampleTest(unittest.TestCase):
    def test_published_file_checksums(self) -> None:
        checksum_path = EXAMPLE_DIR / "SHA256SUMS"
        entries = []
        for line in checksum_path.read_text().splitlines():
            digest, filename = line.split(maxsplit=1)
            entries.append((digest, filename.lstrip("*")))

        self.assertEqual(len(entries), 5)
        for expected_digest, filename in entries:
            path = EXAMPLE_DIR / filename
            self.assertTrue(path.is_file(), filename)
            observed_digest = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(observed_digest, expected_digest, filename)

    def test_clone_proportions_and_lineage(self) -> None:
        with (EXAMPLE_DIR / "clone_proportions.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual({row["patient_id"] for row in rows}, {PATIENT_ID})
        self.assertEqual({row["clone_id"] for row in rows}, {
            "clone_1", "clone_2", "clone_3", "clone_4"
        })

        by_clone = {row["clone_id"]: row for row in rows}
        self.assertEqual(by_clone["clone_1"]["clone_type"], "founding")
        self.assertEqual(by_clone["clone_1"]["parent_clone_id"], "")
        for clone_id in ("clone_2", "clone_3", "clone_4"):
            self.assertEqual(by_clone[clone_id]["clone_type"], "late")
            self.assertEqual(by_clone[clone_id]["parent_clone_id"], "clone_1")

        for timepoint, expected_active in EXPECTED_ACTIVE_CLONES.items():
            column = f"{timepoint}_vaf_pct"
            self.assertEqual(sum(float(row[column]) for row in rows), 100.0)
            observed_active = {
                row["clone_id"] for row in rows if float(row[column]) > 0
            }
            self.assertEqual(observed_active, expected_active)

    def test_mafs_have_required_columns_and_rows(self) -> None:
        maf_paths = sorted(EXAMPLE_DIR.glob("clone_*.maf"))
        self.assertEqual(len(maf_paths), 4)
        for maf_path in maf_paths:
            with maf_path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertIsNotNone(reader.fieldnames, maf_path.name)
                self.assertTrue(
                    REQUIRED_MAF_COLUMNS.issubset(set(reader.fieldnames or [])),
                    maf_path.name,
                )
                self.assertIsNotNone(next(reader, None), maf_path.name)

    def test_manifest_construction_for_each_timepoint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            patient_dir = Path(tmp_dir_name) / PATIENT_ID
            patient_dir.mkdir()
            for source in EXAMPLE_DIR.glob("clone_*.maf"):
                shutil.copy2(source, patient_dir / source.name)
            shutil.copy2(
                EXAMPLE_DIR / "clone_proportions.csv",
                patient_dir / "clone_proportions.csv",
            )
            (patient_dir / "diploid_genotype_germline_hg38.out").touch()

            for timepoint, expected_active in EXPECTED_ACTIVE_CLONES.items():
                rows, selected_timepoint = build_manifest(
                    patient_dir=patient_dir,
                    repo_root=REPO_ROOT,
                    out_dir=patient_dir / f"prepared_hg38_{timepoint}",
                    sex="female",
                    normal_depth=30,
                    tumor_depth=60,
                    timepoint=timepoint,
                    include_zero_clones=True,
                )

                self.assertEqual(selected_timepoint, timepoint)
                self.assertEqual(len(rows), 4)
                self.assertAlmostEqual(
                    sum(float(row["clone_fraction"]) for row in rows), 1.0
                )
                observed_active = {
                    row["clone_id"]
                    for row in rows
                    if float(row["clone_fraction"]) > 0
                }
                self.assertEqual(observed_active, expected_active)
                for row in rows:
                    self.assertAlmostEqual(
                        float(row["simulate_depth"]),
                        60 * float(row["clone_fraction"]),
                    )


if __name__ == "__main__":
    unittest.main()
