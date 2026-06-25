#!/bin/bash
#SBATCH --job-name=scg_merge_tumor
#SBATCH --output=merge_tumor_%x_%j.out
#SBATCH --error=merge_tumor_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$REPO_ROOT"

TIMEPOINT="${1:?Usage: sbatch scripts/slurm/merge_timepoint_tumor.sh t2 [--dry-run]}"
MODE="${2:-}"
DRY_RUN_ARGS=()
if [ "$MODE" = "--dry-run" ]; then
  DRY_RUN_ARGS=(--dry-run)
elif [ -n "$MODE" ]; then
  echo "Unknown mode: $MODE" >&2
  echo "Usage: sbatch scripts/slurm/merge_timepoint_tumor.sh t2 [--dry-run]" >&2
  exit 2
fi

if command -v module >/dev/null 2>&1 && [ "${LOAD_MODULES:-1}" = "1" ]; then
  module purge
  module load oneapi hdf5 python/3.12.1
fi

PATIENT="${PATIENT:?Set PATIENT to the patient directory name before submission}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
MANIFEST="${MANIFEST:-patients/$PATIENT/prepared_hg38_${TIMEPOINT}/final_clone_mutations/patient_manifest.final_clone_mutations.csv}"
CLONE_FASTQ_DIR="${CLONE_FASTQ_DIR:-patients/$PATIENT/tumor_clone_fastqs_independent/$TIMEPOINT}"
OUT_DIR="${OUT_DIR:-patients/$PATIENT/tumor_fastq_${TIMEPOINT}}"
METRICS_DIR="patients/$PATIENT/run_metrics"
RUN_ID="merge_${TIMEPOINT}_${SLURM_JOB_ID:-manual}"
TIME_METRICS="$METRICS_DIR/${RUN_ID}.time.txt"
RUN_SUMMARY="$METRICS_DIR/${RUN_ID}.summary.txt"
START_EPOCH=$(date +%s)

mkdir -p "$METRICS_DIR"

{
  echo "run_id=$RUN_ID"
  echo "started_at=$(date --iso-8601=seconds)"
  echo "patient=$PATIENT"
  echo "timepoint=$TIMEPOINT"
  echo "manifest=$MANIFEST"
  echo "clone_fastq_dir=$CLONE_FASTQ_DIR"
  echo "out_dir=$OUT_DIR"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_job_name=${SLURM_JOB_NAME:-NA}"
  echo "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "slurm_job_nodelist=${SLURM_JOB_NODELIST:-NA}"
  echo "mode=${MODE:-run}"
  if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "=== scontrol show job ==="
    scontrol show job "$SLURM_JOB_ID"
  fi
} > "$RUN_SUMMARY"

echo "=== Environment ==="
command -v "$PYTHON_BIN"
"$PYTHON_BIN" --version

echo "=== Input checks ==="
test -s "$MANIFEST"
test -d "$CLONE_FASTQ_DIR"
ls -lh "$MANIFEST"
ls -lh "$CLONE_FASTQ_DIR"

echo "=== Expected active clones from manifest ==="
awk -F, '
  NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
  $h["clone_fraction"]+0 > 0 {
    print $h["timepoint"], $h["clone_id"], "fraction=" $h["clone_fraction"], "simulate_depth=" $h["simulate_depth"]
  }
' "$MANIFEST"

echo "=== Clone FASTQ inputs ==="
awk -F, '
  NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
  $h["clone_fraction"]+0 > 0 {print $h["clone_id"]}
' "$MANIFEST" | while read -r clone_id; do
  clone_dir="$CLONE_FASTQ_DIR/$clone_id"
  echo "--- $clone_id"
  test -d "$clone_dir"
  ls -lhL "$clone_dir"/*_read1.fq.gz "$clone_dir"/*_read2.fq.gz
done

echo "=== Merge dry run ==="
"$PYTHON_BIN" scripts/pipeline/merge_patient_tumor_fastqs.py "$MANIFEST" \
  --clone-fastq-dir "$CLONE_FASTQ_DIR" \
  --out-dir "$OUT_DIR" \
  --output-prefix tumor \
  --compresslevel 1 \
  --dry-run

if [ "$MODE" = "--dry-run" ]; then
  echo "Dry run requested; not writing merged tumour FASTQs."
  exit 0
fi

echo "=== Merge tumour FASTQs ==="
set +e
/usr/bin/time -v -o "$TIME_METRICS" \
"$PYTHON_BIN" scripts/pipeline/merge_patient_tumor_fastqs.py "$MANIFEST" \
  --clone-fastq-dir "$CLONE_FASTQ_DIR" \
  --out-dir "$OUT_DIR" \
  --output-prefix tumor \
  --compresslevel 1 \
  --overwrite
CMD_EXIT=$?
set -e

END_EPOCH=$(date +%s)
{
  echo "exit_code=$CMD_EXIT"
  echo "finished_at=$(date --iso-8601=seconds)"
  echo "elapsed_seconds=$((END_EPOCH - START_EPOCH))"
  echo "=== final outputs ==="
  if [ -e "$OUT_DIR" ]; then
    ls -lh "$OUT_DIR"
  else
    echo "No output directory yet: $OUT_DIR"
  fi
  echo "=== time metrics file ==="
  cat "$TIME_METRICS" 2>/dev/null || true
} >> "$RUN_SUMMARY"

echo "=== Metrics ==="
echo "Summary: $RUN_SUMMARY"
echo "GNU time: $TIME_METRICS"

if [ "$CMD_EXIT" -eq 0 ]; then
  echo "${TIMEPOINT} tumour merge OK"
fi

exit "$CMD_EXIT"
