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

PATIENT="${PATIENT:-79ce1d89-46d2-5513-c704-212aa1ed97d2}"
NORMAL_TIMEPOINT="${NORMAL_TIMEPOINT:-t0}"
MANIFEST="${MANIFEST:-patients/$PATIENT/prepared_hg38_${NORMAL_TIMEPOINT}/patient_manifest.csv}"
NEAT_CPUS="${NEAT_CPUS:-8}"
SAMTOOLS_CPUS="${SAMTOOLS_CPUS:-16}"

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
which python
python --version
python -m site --user-site

echo "=== Tool check ==="
which samtools
samtools --version | head -n 2
which bgzip
bgzip --version || true

which gzip
which zcat

which gen_reads.py

echo "=== Python package check ==="
python - <<'PY'
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
USER_GEM_HOME="$HOME/.local/share/gem/ruby/3.3.0"

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

python scripts/pipeline/generate_patient_normal_fastqs.py \
  "$MANIFEST" \
  --ruby "$RUBY_BIN" \
  --neat-cpus "$NEAT_CPUS" \
  --samtools-cpus "$SAMTOOLS_CPUS" \
  --materialize symlink \
  --overwrite
