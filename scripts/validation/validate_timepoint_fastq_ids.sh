#!/bin/bash
#SBATCH --job-name=validate_fastq_ids
#SBATCH --output=validate_fastq_ids_%j.out
#SBATCH --error=validate_fastq_ids_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$REPO_ROOT"

TIMEPOINT="${1:?Usage: sbatch scripts/validation/validate_timepoint_fastq_ids.sh t0}"
PATIENT="${PATIENT:?Set PATIENT to the patient directory name before submission}"
PATIENT_DIR="patients/$PATIENT"
MANIFEST="${MANIFEST:-$PATIENT_DIR/prepared_hg38_${TIMEPOINT}/final_clone_mutations/patient_manifest.final_clone_mutations.csv}"
TUMOR_DIR="${TUMOR_DIR:-$PATIENT_DIR/tumor_fastq_${TIMEPOINT}}"
R1="${TUMOR_R1:-$TUMOR_DIR/tumor_read1.fq.gz}"
R2="${TUMOR_R2:-$TUMOR_DIR/tumor_read2.fq.gz}"
OUT_DIR="${OUT_DIR:-$PATIENT_DIR/validation_${TIMEPOINT}/fastq_id_validation}"
SORT_TMP="$OUT_DIR/sort_tmp"
PAIR_JSON="$OUT_DIR/${TIMEPOINT}_pair_order_report.json"
DUP_IDS="$OUT_DIR/${TIMEPOINT}_duplicate_r1_base_ids.txt"
DUP_COUNT="$OUT_DIR/${TIMEPOINT}_duplicate_r1_base_id_count.txt"
SUMMARY="$OUT_DIR/${TIMEPOINT}_fastq_id_validation_summary.txt"
SORT_PARALLEL="${SORT_PARALLEL:-${SLURM_CPUS_PER_TASK:-16}}"
SORT_MEMORY="${SORT_MEMORY:-12G}"
CHECK_DUPLICATES="${CHECK_DUPLICATES:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

if command -v module >/dev/null 2>&1 && [ "${LOAD_MODULES:-1}" = "1" ]; then
  module purge
  module load oneapi hdf5 python/3.12.1
fi

mkdir -p "$OUT_DIR" "$SORT_TMP"

echo "=== Inputs ==="
for f in "$MANIFEST" "$R1" "$R2"; do
  test -s "$f"
  ls -lhL "$f"
done

echo "=== Full R1/R2 pair-order and clone-label check ==="
"$PYTHON_BIN" scripts/validation/validate_fastq_pair_ids.py \
  --read1 "$R1" \
  --read2 "$R2" \
  --output-json "$PAIR_JSON"

"$PYTHON_BIN" - "$PAIR_JSON" "$MANIFEST" "$TIMEPOINT" <<'PY'
import csv
import json
import sys
from pathlib import Path

report_path, manifest_path, timepoint = map(Path, sys.argv[1:])
report = json.loads(report_path.read_text())
with manifest_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

expected = {
    f"tumor_{timepoint}_{row['clone_id']}"
    for row in rows
    if float(row["clone_fraction"]) > 0
}
observed = {label for label, count in report["clone_counts"].items() if count > 0}
if "UNKNOWN" in observed:
    raise SystemExit("FASTQ headers without a recognizable clone label were found")
if observed != expected:
    raise SystemExit(
        "Clone labels do not match the manifest: "
        f"expected={sorted(expected)}, observed={sorted(observed)}"
    )
print("Clone labels match the active manifest clones.")
PY

if [ "$CHECK_DUPLICATES" = "1" ]; then
  echo "=== Exact duplicate read-name check on R1 base IDs ==="
  zcat "$R1" \
    | awk 'NR % 4 == 1 {h=$1; sub(/^@/, "", h); sub(/\/1$/, "", h); print h}' \
    | LC_ALL=C sort --parallel="$SORT_PARALLEL" -T "$SORT_TMP" -S "$SORT_MEMORY" \
    | uniq -d \
    > "$DUP_IDS"
  wc -l < "$DUP_IDS" > "$DUP_COUNT"
else
  echo "Duplicate-ID check disabled with CHECK_DUPLICATES=$CHECK_DUPLICATES"
  : > "$DUP_IDS"
  echo "SKIPPED" > "$DUP_COUNT"
fi

"$PYTHON_BIN" - "$PAIR_JSON" "$DUP_COUNT" "$PATIENT" "$TIMEPOINT" "$TUMOR_DIR" > "$SUMMARY" <<'PY'
import json
import sys
from pathlib import Path

pair_json, duplicate_count, patient, timepoint, tumor_dir = sys.argv[1:]
report = json.loads(Path(pair_json).read_text())
duplicates = Path(duplicate_count).read_text().strip()
print(f"patient={patient}")
print(f"timepoint={timepoint}")
print(f"tumor_dir={tumor_dir}")
print(f"status={report['status']}")
print(f"records_checked={report['records_checked']}")
print(f"mate_id_mismatches={report['mate_id_mismatches']}")
print(f"malformed_records={report['malformed_records']}")
print("clone_counts=" + json.dumps(report["clone_counts"], sort_keys=True))
print(f"duplicate_id_count={duplicates}")
PY

cat "$SUMMARY"
if [ "$(cat "$DUP_COUNT")" != "0" ] && [ "$(cat "$DUP_COUNT")" != "SKIPPED" ]; then
  echo "Duplicate R1 base IDs found. First duplicated IDs:" >&2
  head "$DUP_IDS" >&2
  exit 1
fi

echo "FASTQ ID validation passed."
