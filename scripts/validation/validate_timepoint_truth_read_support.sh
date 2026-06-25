#!/bin/bash
#SBATCH --job-name=validate_timepoint_support
#SBATCH --output=validate_timepoint_support_%j.out
#SBATCH --error=validate_timepoint_support_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$REPO_ROOT"

TIMEPOINT="${1:?Usage: sbatch scripts/validation/validate_timepoint_truth_read_support.sh t1}"

if command -v module >/dev/null 2>&1 && [ "${LOAD_MODULES:-1}" = "1" ]; then
  module purge
  module load oneapi hdf5 python/3.12.1
  module load samtools/1.19.2
  module load htslib/1.19.1
fi

if [ -n "${PYTHONUSERBASE:-}" ]; then
  export PATH="$PYTHONUSERBASE/bin:$PATH"
fi

PATIENT="${PATIENT:?Set PATIENT to the patient directory name before submission}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
PATIENT_DIR="patients/$PATIENT"
OUT_DIR="${OUT_DIR:-$PATIENT_DIR/validation_${TIMEPOINT}}"
SUPPORT_DIR="$OUT_DIR/read_support"
METRICS_DIR="$PATIENT_DIR/run_metrics"
RUN_ID="validation_${TIMEPOINT}_read_support_${SLURM_JOB_ID:-manual}"
TIME_METRICS="$METRICS_DIR/${RUN_ID}.time.txt"
RUN_SUMMARY="$METRICS_DIR/${RUN_ID}.summary.txt"

TRUTH_VCF="$OUT_DIR/${TIMEPOINT}_truth.vcf.gz"
CALLED_VCF="$OUT_DIR/${TIMEPOINT}.called.pass.norm.vcf.gz"
TUMOR_BAM="$OUT_DIR/tumor_${TIMEPOINT}.bam"
NORMAL_BAM="$OUT_DIR/normal.bam"
SUPPORT_TSV="$SUPPORT_DIR/${TIMEPOINT}_truth_variant_read_support.tsv"
TRUTH_ONLY_TSV="$SUPPORT_DIR/${TIMEPOINT}_truth_only_variant_read_support.tsv"
SUMMARY_TSV="$SUPPORT_DIR/${TIMEPOINT}_truth_variant_read_support_summary.tsv"

if [ -s "scripts/validation/check_truth_variant_read_support.py" ]; then
  SUPPORT_SCRIPT="scripts/validation/check_truth_variant_read_support.py"
elif [ -s "scripts/check_truth_variant_read_support.py" ]; then
  SUPPORT_SCRIPT="scripts/check_truth_variant_read_support.py"
else
  echo "Cannot find check_truth_variant_read_support.py under scripts/ or scripts/validation/" >&2
  exit 2
fi

mkdir -p "$SUPPORT_DIR" "$METRICS_DIR"
START_EPOCH=$(date +%s)

{
  echo "run_id=$RUN_ID"
  echo "started_at=$(date --iso-8601=seconds)"
  echo "patient=$PATIENT"
  echo "timepoint=$TIMEPOINT"
  echo "out_dir=$OUT_DIR"
  echo "support_dir=$SUPPORT_DIR"
  echo "support_script=$SUPPORT_SCRIPT"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_job_name=${SLURM_JOB_NAME:-NA}"
  echo "slurm_ntasks=${SLURM_NTASKS:-NA}"
  echo "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "slurm_job_nodelist=${SLURM_JOB_NODELIST:-NA}"
  if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "=== scontrol show job ==="
    scontrol show job "$SLURM_JOB_ID"
  fi
} > "$RUN_SUMMARY"

echo "=== Tool versions ==="
command -v "$PYTHON_BIN"
"$PYTHON_BIN" --version
"$PYTHON_BIN" - <<'PY'
import pysam
print("pysam", pysam.__version__)
PY
which samtools
samtools --version | head -n 2

echo "=== Input checks ==="
for f in "$TRUTH_VCF" "$CALLED_VCF" "$TUMOR_BAM" "$TUMOR_BAM.bai" "$NORMAL_BAM" "$NORMAL_BAM.bai" "$SUPPORT_SCRIPT"; do
  test -s "$f"
  ls -lh "$f"
done

echo "=== BAM integrity ==="
samtools quickcheck -v "$TUMOR_BAM" "$NORMAL_BAM"

echo "=== Count read support at truth variants ==="
set +e
/usr/bin/time -v -o "$TIME_METRICS" \
"$PYTHON_BIN" "$SUPPORT_SCRIPT" \
  --truth-vcf "$TRUTH_VCF" \
  --called-vcf "$CALLED_VCF" \
  --tumor-bam "$TUMOR_BAM" \
  --normal-bam "$NORMAL_BAM" \
  --output-tsv "$SUPPORT_TSV" \
  --summary-tsv "$SUMMARY_TSV" \
  --min-mapq 20 \
  --min-baseq 10
CMD_EXIT=$?
set -e

if [ "$CMD_EXIT" -eq 0 ]; then
  awk -F'\t' 'NR == 1 || $6 == "truth_only"' "$SUPPORT_TSV" > "$TRUTH_ONLY_TSV"
else
  echo "Read-support counter failed before truth-only table extraction." >&2
fi

END_EPOCH=$(date +%s)
{
  echo "exit_code=$CMD_EXIT"
  echo "finished_at=$(date --iso-8601=seconds)"
  echo "elapsed_seconds=$((END_EPOCH - START_EPOCH))"
  echo "=== output files ==="
  find "$SUPPORT_DIR" -maxdepth 1 -type f -print | sort | xargs -r ls -lh
  echo "=== read-support summary ==="
  cat "$SUMMARY_TSV" 2>/dev/null || true
  echo "=== time metrics file ==="
  cat "$TIME_METRICS" 2>/dev/null || true
} >> "$RUN_SUMMARY"

echo "=== Read-support summary ==="
cat "$SUMMARY_TSV" 2>/dev/null || true
echo "=== Metrics ==="
echo "Summary: $RUN_SUMMARY"
echo "GNU time: $TIME_METRICS"

exit "$CMD_EXIT"
