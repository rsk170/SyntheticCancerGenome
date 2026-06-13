#!/bin/bash
#SBATCH --job-name=validate_t0_ids
#SBATCH --output=validate_t0_ids_%j.out
#SBATCH --error=validate_t0_ids_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$REPO_ROOT"

if command -v module >/dev/null 2>&1 && [ "${LOAD_MODULES:-1}" = "1" ]; then
  module purge
  module load oneapi hdf5 python/3.12.1
fi

PATIENT="${PATIENT:-79ce1d89-46d2-5513-c704-212aa1ed97d2}"
PATIENT_DIR="patients/$PATIENT"
if [ -d "$PATIENT_DIR/tumor_fastq_t0" ]; then
  TUMOR_DIR="$PATIENT_DIR/tumor_fastq_t0"
else
  TUMOR_DIR="$PATIENT_DIR/tumor_fastq"
fi

R1="$TUMOR_DIR/tumor_read1.fq.gz"
R2="$TUMOR_DIR/tumor_read2.fq.gz"
OUT_DIR="$PATIENT_DIR/validation_t0/fastq_id_validation"
SORT_TMP="$OUT_DIR/sort_tmp"
PAIR_JSON="$OUT_DIR/t0_pair_order_report.json"
DUP_IDS="$OUT_DIR/t0_duplicate_r1_base_ids.txt"
DUP_COUNT="$OUT_DIR/t0_duplicate_r1_base_id_count.txt"
SUMMARY="$OUT_DIR/t0_fastq_id_validation_summary.txt"

mkdir -p "$OUT_DIR" "$SORT_TMP"

echo "=== Inputs ==="
for f in "$R1" "$R2"; do
  test -s "$f"
  ls -lhL "$f"
done

echo "=== Full R1/R2 pair-order and clone-label check ==="
python scripts/validation/validate_fastq_pair_ids.py \
  --read1 "$R1" \
  --read2 "$R2" \
  --output-json "$PAIR_JSON"

echo "=== Exact duplicate read-name check on R1 base IDs ==="
echo "This checks full read-name uniqueness after removing @ and /1 from R1 headers."
echo "It uses external sort and may take a while for billion-read FASTQs."

zcat "$R1" \
  | awk 'NR % 4 == 1 {h=$1; sub(/^@/, "", h); sub(/\/1$/, "", h); print h}' \
  | LC_ALL=C sort --parallel=16 -T "$SORT_TMP" -S 12G \
  | uniq -d \
  > "$DUP_IDS"

wc -l "$DUP_IDS" > "$DUP_COUNT"

{
  echo "patient=$PATIENT"
  echo "tumor_dir=$TUMOR_DIR"
  echo "read1=$R1"
  echo "read2=$R2"
  echo "pair_order_report=$PAIR_JSON"
  echo "duplicate_ids=$DUP_IDS"
  echo "duplicate_id_count=$(awk '{print $1}' "$DUP_COUNT")"
  echo "=== pair-order status ==="
  python - <<PY
import json
from pathlib import Path
data = json.loads(Path("$PAIR_JSON").read_text())
print("status=" + data["status"])
print("records_checked=" + str(data["records_checked"]))
print("mate_id_mismatches=" + str(data["mate_id_mismatches"]))
print("malformed_records=" + str(data["malformed_records"]))
print("clone_counts=" + json.dumps(data["clone_counts"], sort_keys=True))
PY
} > "$SUMMARY"

cat "$SUMMARY"

if [ "$(awk '{print $1}' "$DUP_COUNT")" != "0" ]; then
  echo "Duplicate R1 base IDs found. First duplicated IDs:" >&2
  head "$DUP_IDS" >&2
  exit 1
fi

echo "FASTQ ID validation passed."
