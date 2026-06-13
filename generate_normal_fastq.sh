#!/bin/bash
#SBATCH --job-name=scg_normal
#SBATCH --chdir=/gpfs/projects/bsc82/bsc720159/SyntheticCancerGenome
#SBATCH --output=serial_%j.out
#SBATCH --error=serial_%j.err
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --account=bsc82
#SBATCH --qos=gp_bscls
#SBATCH --cpus-per-task=64

set -euo pipefail

echo "=== Load Python ==="
module purge
module load oneapi hdf5 python/3.12.1
module load samtools/1.19.2
module load htslib/1.19.1

export PYTHONUSERBASE=/gpfs/projects/bsc82/bsc720159/python_packages
export PATH="$PYTHONUSERBASE/bin:$PATH"
export PATH="$HOME/.rbbt/software/opt/bin:$PATH"

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

"$RUBY_BIN" --version
"$RUBY_BIN" -e 'puts Gem.ruby; puts Gem.dir; puts Gem.path'
"$RUBY_BIN" -rrbbt-util -e 'puts "rbbt-util OK"'
"$RUBY_BIN" -e 'require "rbbt/sources/organism"; puts "rbbt-sources OK"'
"$RUBY_BIN" -Ilib -rrbbt-util -e 'require "./workflow"; puts SyntheticCancerGenome.tasks.keys.grep(/normal/)'

echo "=== Rbbt jobs path ==="
readlink -f ~/.rbbt/var/jobs

echo "=== RUN ==="
NEAT_CPUS=8
SAMTOOLS_CPUS=16
echo "Slurm CPUs reserved: ${SLURM_CPUS_PER_TASK:-64}"
echo "NEAT chromosome-parallel workers: $NEAT_CPUS"
echo "samtools merge threads: $SAMTOOLS_CPUS"

python scripts/generate_patient_normal_fastqs.py \
  patients/79ce1d89-46d2-5513-c704-212aa1ed97d2/prepared_hg38_t0/patient_manifest.csv \
  --ruby "$RUBY_BIN" \
  --neat-cpus "$NEAT_CPUS" \
  --samtools-cpus "$SAMTOOLS_CPUS" \
  --materialize symlink \
  --overwrite
