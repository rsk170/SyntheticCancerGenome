#!/bin/bash
#SBATCH --job-name=scg_normal
#SBATCH --output=serial_%j.out
#SBATCH --error=serial_%j.err
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=64

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
cd "$REPO_ROOT"

PATIENT="${PATIENT:?Set PATIENT to the patient directory name before submission}"
NORMAL_TIMEPOINT="${NORMAL_TIMEPOINT:-t0}"
MANIFEST="${MANIFEST:-patients/$PATIENT/prepared_hg38_${NORMAL_TIMEPOINT}/patient_manifest.csv}"
REFERENCE_METADATA="${REFERENCE_METADATA:-}"
GERMLINE_METADATA="${GERMLINE_METADATA:-}"
OUT_DIR="${OUT_DIR:-}"
NEAT_CPUS="${NEAT_CPUS:-8}"
SAMTOOLS_CPUS="${SAMTOOLS_CPUS:-16}"
MATERIALIZE="${MATERIALIZE:-symlink}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
METRICS_DIR="${METRICS_DIR:-patients/$PATIENT/run_metrics}"
RUN_ID="normal_${SLURM_JOB_ID:-manual}"
TIME_METRICS="$METRICS_DIR/${RUN_ID}.time.txt"
RUN_SUMMARY="$METRICS_DIR/${RUN_ID}.summary.txt"
START_EPOCH=$(date +%s)

mkdir -p "$METRICS_DIR"

{
  echo "run_id=$RUN_ID"
  echo "started_at=$(date --iso-8601=seconds)"
  echo "patient=$PATIENT"
  echo "manifest=$MANIFEST"
  echo "out_dir=${OUT_DIR:-default}"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_job_name=${SLURM_JOB_NAME:-NA}"
  echo "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "neat_cpus=$NEAT_CPUS"
  echo "samtools_cpus=$SAMTOOLS_CPUS"
  echo "materialize=$MATERIALIZE"
} > "$RUN_SUMMARY"

echo "=== Load Python ==="
if command -v module >/dev/null 2>&1 && [ "${LOAD_MODULES:-1}" = "1" ]; then
  module purge
  module load oneapi hdf5 python/3.12.1
  module load samtools/1.19.2
  module load htslib/1.19.1
fi

if [ -n "${PYTHONUSERBASE:-}" ]; then
  export PATH="$PYTHONUSERBASE/bin:$PATH"
fi
RBBT_OPT_BIN="${RBBT_OPT_BIN:-$HOME/.rbbt/software/opt/bin}"
if [ -d "$RBBT_OPT_BIN" ]; then
  export PATH="$RBBT_OPT_BIN:$PATH"
fi

echo "=== Python environment ==="
command -v "$PYTHON_BIN"
"$PYTHON_BIN" --version
"$PYTHON_BIN" -m site --user-site

echo "=== Tool check ==="
which samtools
samtools --version | head -n 2
which bgzip
bgzip --version || true

which gzip
which zcat

which gen_reads.py

echo "=== Python package check ==="
"$PYTHON_BIN" - <<'PY'
import numpy
import Bio
import pandas
import pysam
import matplotlib
import matplotlib_venn
print("Python packages OK")
PY

echo "=== Ruby / Rbbt environment ==="
RBBT_ENV="${RBBT_ENV:-}"
if [ -z "${RUBY_BIN:-}" ]; then
  if [ -n "$RBBT_ENV" ]; then
    RUBY_BIN="$RBBT_ENV/bin/ruby"
  else
    RUBY_BIN="ruby"
  fi
fi
USER_GEM_HOME="${USER_GEM_HOME:-$HOME/.local/share/gem/ruby/$("$RUBY_BIN" -e 'print RbConfig::CONFIG["ruby_version"]')}"

if [ -n "$RBBT_ENV" ]; then
  RUBY_GEM_HOME="${RUBY_GEM_HOME:-$RBBT_ENV/share/rubygems}"
  export GEM_HOME="$RUBY_GEM_HOME"
  if [ -d "$USER_GEM_HOME" ]; then
    export GEM_PATH="$RUBY_GEM_HOME:$USER_GEM_HOME"
  else
    export GEM_PATH="$RUBY_GEM_HOME"
  fi
  export PATH="$RUBY_GEM_HOME/bin:$PATH"
fi

"$RUBY_BIN" --version
"$RUBY_BIN" -e 'puts Gem.ruby; puts Gem.dir; puts Gem.path'
"$RUBY_BIN" -rrbbt-util -e 'puts "rbbt-util OK"'
"$RUBY_BIN" -e 'require "rbbt/sources/organism"; puts "rbbt-sources OK"'
"$RUBY_BIN" -Ilib -rrbbt-util -e 'require "./workflow"; puts SyntheticCancerGenome.tasks.keys.grep(/normal/)'

echo "=== Rbbt jobs path ==="
readlink -f ~/.rbbt/var/jobs

echo "=== RUN ==="
echo "Slurm CPUs reserved: ${SLURM_CPUS_PER_TASK:-64}"
echo "NEAT chromosome-parallel workers: $NEAT_CPUS"
echo "samtools merge threads: $SAMTOOLS_CPUS"
echo "Manifest: $MANIFEST"
echo "Reference metadata: ${REFERENCE_METADATA:-default}"
echo "Germline metadata: ${GERMLINE_METADATA:-default}"
echo "Output directory: ${OUT_DIR:-default}"
echo "Materialize mode: $MATERIALIZE"

EXTRA_ARGS=()
if [ -n "$REFERENCE_METADATA" ]; then
  EXTRA_ARGS+=(--reference-metadata "$REFERENCE_METADATA")
fi
if [ -n "$GERMLINE_METADATA" ]; then
  EXTRA_ARGS+=(--germline-metadata "$GERMLINE_METADATA")
fi
if [ -n "$OUT_DIR" ]; then
  EXTRA_ARGS+=(--out-dir "$OUT_DIR")
fi

set +e
/usr/bin/time -v -o "$TIME_METRICS" \
"$PYTHON_BIN" scripts/pipeline/generate_patient_normal_fastqs.py \
  "$MANIFEST" \
  "${EXTRA_ARGS[@]}" \
  --ruby "$RUBY_BIN" \
  --neat-cpus "$NEAT_CPUS" \
  --samtools-cpus "$SAMTOOLS_CPUS" \
  --no-rename-reads \
  --materialize "$MATERIALIZE" \
  --overwrite
CMD_EXIT=$?
set -e

END_EPOCH=$(date +%s)
{
  echo "exit_code=$CMD_EXIT"
  echo "finished_at=$(date --iso-8601=seconds)"
  echo "elapsed_seconds=$((END_EPOCH - START_EPOCH))"
  echo "=== output size ==="
  if [ -n "$OUT_DIR" ]; then
    du -sh "$OUT_DIR" 2>/dev/null || true
  else
    du -sh "patients/$PATIENT/normal_fastq" 2>/dev/null || true
  fi
  echo "=== time metrics file ==="
  cat "$TIME_METRICS" 2>/dev/null || true
} >> "$RUN_SUMMARY"

echo "=== Metrics ==="
echo "Summary: $RUN_SUMMARY"
echo "GNU time: $TIME_METRICS"
exit "$CMD_EXIT"
