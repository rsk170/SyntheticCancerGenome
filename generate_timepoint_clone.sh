#!/bin/bash
#SBATCH --job-name=scg_timepoint_clone
#SBATCH --chdir=/gpfs/projects/bsc82/bsc720159/SyntheticCancerGenome
#SBATCH --output=timepoint_clone_%j.out
#SBATCH --error=timepoint_clone_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=2-00:00:00
#SBATCH --account=bsc82
#SBATCH --qos=gp_bscls

set -euo pipefail

TIMEPOINT="${1:?Usage: sbatch generate_timepoint_clone.sh t1 clone_1 [--dry-run]}"
CLONE_ID="${2:?Usage: sbatch generate_timepoint_clone.sh t1 clone_1 [--dry-run]}"
MODE="${3:-}"
DRY_RUN_ARGS=()
if [ "$MODE" = "--dry-run" ]; then
  DRY_RUN_ARGS=(--dry-run)
elif [ -n "$MODE" ]; then
  echo "Unknown mode: $MODE" >&2
  echo "Usage: sbatch generate_timepoint_clone.sh t1 clone_1 [--dry-run]" >&2
  exit 2
fi

module purge
module load oneapi hdf5 python/3.12.1
module load samtools/1.19.2
module load htslib/1.19.1

export PYTHONUSERBASE=/gpfs/projects/bsc82/bsc720159/python_packages
export PATH="$PYTHONUSERBASE/bin:$HOME/.rbbt/software/opt/bin:$PATH"

RBBT_ENV=/gpfs/projects/bsc82/bsc720159/conda_envs/rbbt_env
RUBY_BIN="$RBBT_ENV/bin/ruby"
RUBY_GEM_HOME="$RBBT_ENV/share/rubygems"
USER_GEM_HOME="$HOME/.local/share/gem/ruby/3.3.0"

export GEM_HOME="$RUBY_GEM_HOME"
if [ -d "$USER_GEM_HOME" ]; then
  export GEM_PATH="$RUBY_GEM_HOME:$USER_GEM_HOME"
else
  export GEM_PATH="$RUBY_GEM_HOME"
fi
export PATH="$RUBY_GEM_HOME/bin:$PATH"

PATIENT=79ce1d89-46d2-5513-c704-212aa1ed97d2
MANIFEST="patients/$PATIENT/prepared_hg38_${TIMEPOINT}/final_clone_mutations/patient_manifest.final_clone_mutations.csv"
OUT_DIR="patients/$PATIENT/tumor_clone_fastqs_independent/$TIMEPOINT"
METRICS_DIR="patients/$PATIENT/run_metrics"
RUN_ID="${TIMEPOINT}_${CLONE_ID}_${SLURM_JOB_ID:-manual}"
TIME_METRICS="$METRICS_DIR/${RUN_ID}.time.txt"
RUN_SUMMARY="$METRICS_DIR/${RUN_ID}.summary.txt"

NEAT_CPUS=8
SAMTOOLS_CPUS=16
START_EPOCH=$(date +%s)

mkdir -p "$METRICS_DIR"

{
  echo "run_id=$RUN_ID"
  echo "started_at=$(date --iso-8601=seconds)"
  echo "patient=$PATIENT"
  echo "timepoint=$TIMEPOINT"
  echo "clone_id=$CLONE_ID"
  echo "manifest=$MANIFEST"
  echo "out_dir=$OUT_DIR"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_job_name=${SLURM_JOB_NAME:-NA}"
  echo "slurm_ntasks=${SLURM_NTASKS:-NA}"
  echo "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "slurm_job_nodelist=${SLURM_JOB_NODELIST:-NA}"
  echo "neat_cpus=$NEAT_CPUS"
  echo "samtools_cpus=$SAMTOOLS_CPUS"
  echo "mode=${MODE:-run}"
  if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "=== scontrol show job ==="
    scontrol show job "$SLURM_JOB_ID"
  fi
} > "$RUN_SUMMARY"

echo "=== Environment ==="
which python
python --version
which samtools
samtools --version | head -n 2
command -v bgzip
command -v gen_reads.py
"$RUBY_BIN" --version
"$RUBY_BIN" -rrbbt-util -e 'puts "rbbt-util OK"'
"$RUBY_BIN" -e 'require "rbbt/sources/organism"; puts "rbbt-sources OK"'
"$RUBY_BIN" -Ilib -rrbbt-util -e 'require "./workflow"; puts SyntheticCancerGenome.tasks.keys.grep(/tumor|normal/)'

echo "=== Input checks ==="
test -s "$MANIFEST"
ls -lh "$MANIFEST"
awk -F, -v clone="$CLONE_ID" '
  NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
  $h["clone_id"] == clone {
    print "timepoint=" $h["timepoint"]
    print "clone_id=" $h["clone_id"]
    print "clone_fraction=" $h["clone_fraction"]
    print "tumor_target_depth=" $h["tumor_target_depth"]
    print "simulate_depth=" $h["simulate_depth"]
    print "final_mutations=" $h["final_clone_hg38_mutations_path"]
  }
' "$MANIFEST"

echo "=== Rbbt jobs path ==="
readlink -f ~/.rbbt/var/jobs

echo "=== Timepoint clone run ==="
echo "Patient: $PATIENT"
echo "Timepoint: $TIMEPOINT"
echo "Clone: $CLONE_ID"
echo "Manifest: $MANIFEST"
echo "Output directory: $OUT_DIR"
echo "Mode: ${MODE:-run}"
echo "NEAT chromosome-parallel workers: $NEAT_CPUS"
echo "samtools merge threads: $SAMTOOLS_CPUS"

set +e
/usr/bin/time -v -o "$TIME_METRICS" \
python scripts/generate_patient_clone_tumor_fastqs.py \
  "$MANIFEST" \
  --clone-id "$CLONE_ID" \
  --out-dir "$OUT_DIR" \
  --sample-name-template "{patient_id}_tumor_{timepoint}_{clone_id}" \
  --job-name-template "{patient_id}_tumor_{timepoint}_{clone_id}" \
  --output-prefix-template "{clone_id}" \
  --ruby "$RUBY_BIN" \
  --neat-cpus "$NEAT_CPUS" \
  --samtools-cpus "$SAMTOOLS_CPUS" \
  --no-rename-reads \
  --materialize symlink \
  --overwrite \
  "${DRY_RUN_ARGS[@]}"
CMD_EXIT=$?
set -e

END_EPOCH=$(date +%s)
{
  echo "exit_code=$CMD_EXIT"
  echo "finished_at=$(date --iso-8601=seconds)"
  echo "elapsed_seconds=$((END_EPOCH - START_EPOCH))"
  echo "=== output size ==="
  du -sh "$OUT_DIR" 2>/dev/null || true
  echo "=== output files ==="
  if [ -e "$OUT_DIR" ]; then
    find "$OUT_DIR" -maxdepth 2 \( -type f -o -type l \) -print | sort | xargs -r ls -lhL
  else
    echo "No output directory yet: $OUT_DIR"
  fi
  echo "=== time metrics file ==="
  cat "$TIME_METRICS" 2>/dev/null || true
} >> "$RUN_SUMMARY"

echo "=== Metrics ==="
echo "Summary: $RUN_SUMMARY"
echo "GNU time: $TIME_METRICS"
exit "$CMD_EXIT"
