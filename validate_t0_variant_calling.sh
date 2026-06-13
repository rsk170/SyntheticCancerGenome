#!/bin/bash
#SBATCH --job-name=validate_t0_vc
#SBATCH --chdir=/gpfs/projects/bsc82/bsc720159/SyntheticCancerGenome
#SBATCH --output=validate_t0_vc_%j.out
#SBATCH --error=validate_t0_vc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=2-00:00:00
#SBATCH --account=bsc82
#SBATCH --qos=gp_bscls

set -euo pipefail

module purge
module load oneapi hdf5
module load bwa/0.7.17
module load samtools/1.19.2
module load java-openjdk/17.0.11+9
module load gatk/4.5.0.0
module load bcftools/1.19
module load htslib/1.19.1

PATIENT=79ce1d89-46d2-5513-c704-212aa1ed97d2
PATIENT_DIR="patients/$PATIENT"
OUT_DIR="$PATIENT_DIR/validation_t0"
REF_DIR="$OUT_DIR/reference"
METRICS_DIR="$PATIENT_DIR/run_metrics"
RUN_ID="validation_t0_${SLURM_JOB_ID:-manual}"
TIME_METRICS="$METRICS_DIR/${RUN_ID}.time.txt"
RUN_SUMMARY="$METRICS_DIR/${RUN_ID}.summary.txt"

BWA_THREADS=32
SORT_THREADS=16
INDEX_THREADS=16
SORT_MEM=3G

NORMAL_R1="$PATIENT_DIR/normal_fastq/normal_read1.fq.gz"
NORMAL_R2="$PATIENT_DIR/normal_fastq/normal_read2.fq.gz"
if [ -d "$PATIENT_DIR/tumor_fastq_t0" ]; then
  TUMOR_DIR="$PATIENT_DIR/tumor_fastq_t0"
else
  TUMOR_DIR="$PATIENT_DIR/tumor_fastq"
fi
TUMOR_R1="$TUMOR_DIR/tumor_read1.fq.gz"
TUMOR_R2="$TUMOR_DIR/tumor_read2.fq.gz"

TRUTH_CLONE1="$PATIENT_DIR/prepared_hg38_t0/hg38_vcf/clone_1.hg38.normalized.vcf"
TRUTH_CLONE3="$PATIENT_DIR/prepared_hg38_t0/hg38_vcf/clone_3.hg38.normalized.vcf"

mkdir -p "$OUT_DIR" "$REF_DIR" "$METRICS_DIR"
START_EPOCH=$(date +%s)

{
  echo "run_id=$RUN_ID"
  echo "started_at=$(date --iso-8601=seconds)"
  echo "patient=$PATIENT"
  echo "out_dir=$OUT_DIR"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "slurm_job_name=${SLURM_JOB_NAME:-NA}"
  echo "slurm_ntasks=${SLURM_NTASKS:-NA}"
  echo "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}"
  echo "slurm_job_nodelist=${SLURM_JOB_NODELIST:-NA}"
  echo "bwa_threads=$BWA_THREADS"
  echo "sort_threads=$SORT_THREADS"
  echo "index_threads=$INDEX_THREADS"
  if command -v scontrol >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "=== scontrol show job ==="
    scontrol show job "$SLURM_JOB_ID"
  fi
} > "$RUN_SUMMARY"

echo "=== Tool versions ==="
which bwa
bwa 2>&1 | head -n 3 || true
which samtools
samtools --version | head -n 2
which java
java -version
which gatk
gatk --version
which bcftools
bcftools --version | head -n 2
which bgzip
bgzip --version

echo "=== Input checks ==="
for f in "$NORMAL_R1" "$NORMAL_R2" "$TUMOR_R1" "$TUMOR_R2" "$TRUTH_CLONE1" "$TRUTH_CLONE3"; do
  test -s "$f"
  ls -lhL "$f"
done

find_reference() {
  if [ -n "${HG38_REF:-}" ] && [ -s "$HG38_REF" ]; then
    printf '%s\n' "$HG38_REF"
    return 0
  fi

  local candidates=(
    "$PATIENT_DIR/prepared_hg38_t0/ref_cache/hg38/hg38.fa"
    "$PATIENT_DIR/prepared_hg38_t0/ref_cache/hg38/hg38.fa.gz"
    "$PATIENT_DIR/prepared_hg38_t1/ref_cache/hg38/hg38.fa"
    "$PATIENT_DIR/prepared_hg38_t1/ref_cache/hg38/hg38.fa.gz"
    "$PATIENT_DIR/prepared_hg38_t2/ref_cache/hg38/hg38.fa"
    "$PATIENT_DIR/prepared_hg38_t2/ref_cache/hg38/hg38.fa.gz"
    "$PWD/hg38.fa"
    "$PWD/hg38.fa.gz"
    "/gpfs/projects/bsc82/bsc720159/hg38.fa"
    "/gpfs/projects/bsc82/bsc720159/hg38.fa.gz"
    "/gpfs/projects/bsc82/bsc720159/references/hg38.fa"
    "/gpfs/projects/bsc82/bsc720159/references/hg38.fa.gz"
    "$HOME/.rbbt/share/organisms/Hsa/hg38/hg38.fa"
    "$HOME/.rbbt/share/organisms/Hsa/hg38/hg38.fa.gz"
  )

  local candidate
  for candidate in "${candidates[@]}"; do
    if [ -s "$candidate" ]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

SOURCE_REF="$(find_reference || true)"
if [ -z "$SOURCE_REF" ]; then
  cat >&2 <<'EOF'
Could not find a standard hg38 FASTA.

Set HG38_REF when submitting, for example:
  sbatch --export=ALL,HG38_REF=/path/to/hg38.fa validate_t0_variant_calling.sh

Do not use the patient sex-aware ploidy reference here; validation alignment should use standard hg38 contig names.
EOF
  exit 2
fi

echo "=== Reference ==="
echo "Source reference: $SOURCE_REF"
REF="$REF_DIR/hg38.validation.fa"
if [ ! -s "$REF" ]; then
  rm -f "$REF"
  if [[ "$SOURCE_REF" == *.gz ]]; then
    echo "Decompressing reference copy for alignment: $REF"
    zcat "$SOURCE_REF" > "$REF"
  else
    echo "Linking reference for alignment: $REF"
    ln -s "$(readlink -f "$SOURCE_REF")" "$REF"
  fi
fi

FIRST_CONTIG="$(grep -m 1 '^>' "$REF" | sed 's/^>//; s/[[:space:]].*//')"
if [[ "$FIRST_CONTIG" == copy-* ]]; then
  echo "Reference appears to be a copy-labelled simulation/ploidy FASTA: $FIRST_CONTIG" >&2
  echo "Use a standard hg38 FASTA with contigs like chr1, chr2, ..." >&2
  exit 2
fi
if ! grep -q '^>chr1\([[:space:]]\|$\)' "$REF"; then
  echo "Reference does not appear to use chr-prefixed contigs, but truth VCFs do." >&2
  echo "Use an hg38 FASTA with contigs like chr1, chr2, ..." >&2
  exit 2
fi

if [ ! -s "$REF.fai" ]; then
  samtools faidx "$REF"
fi
if [ ! -s "$REF_DIR/hg38.validation.dict" ]; then
  gatk CreateSequenceDictionary -R "$REF" -O "$REF_DIR/hg38.validation.dict"
fi
if [ ! -s "$REF.bwt" ]; then
  bwa index "$REF"
fi

echo "=== Build t0 truth VCF and target BED ==="
(
  bcftools view -h "$TRUTH_CLONE1"
  bcftools view -H "$TRUTH_CLONE1"
  bcftools view -H "$TRUTH_CLONE3"
) | bcftools sort -Oz -o "$OUT_DIR/t0_truth.raw.vcf.gz"

bcftools norm -d exact -f "$REF" \
  -Oz -o "$OUT_DIR/t0_truth.vcf.gz" \
  "$OUT_DIR/t0_truth.raw.vcf.gz"
bcftools index -f "$OUT_DIR/t0_truth.vcf.gz"

bcftools query -f '%CHROM\t%POS0\t%END\n' "$OUT_DIR/t0_truth.vcf.gz" \
  | awk 'BEGIN{OFS="\t"} {s=$2-20; if(s<0)s=0; print $1,s,$3+20}' \
  | sort -k1,1 -k2,2n \
  > "$OUT_DIR/t0_truth.padded.bed"

echo "Truth variants:"
bcftools view -H "$OUT_DIR/t0_truth.vcf.gz" | wc -l

run_validation() {
  echo "=== Align normal ==="
  if [ ! -s "$OUT_DIR/normal.bam" ]; then
    bwa mem -t "$BWA_THREADS" -R '@RG\tID:normal_t0\tSM:normal\tPL:ILLUMINA' \
      "$REF" "$NORMAL_R1" "$NORMAL_R2" \
      | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -o "$OUT_DIR/normal.bam" -
  fi
  samtools index -@ "$INDEX_THREADS" "$OUT_DIR/normal.bam"

  echo "=== Align tumor t0 ==="
  if [ ! -s "$OUT_DIR/tumor_t0.bam" ]; then
    bwa mem -t "$BWA_THREADS" -R '@RG\tID:tumor_t0\tSM:tumor_t0\tPL:ILLUMINA' \
      "$REF" "$TUMOR_R1" "$TUMOR_R2" \
      | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -o "$OUT_DIR/tumor_t0.bam" -
  fi
  samtools index -@ "$INDEX_THREADS" "$OUT_DIR/tumor_t0.bam"

  echo "=== BAM checks ==="
  samtools quickcheck -v "$OUT_DIR/normal.bam" "$OUT_DIR/tumor_t0.bam"
  samtools flagstat -@ "$INDEX_THREADS" "$OUT_DIR/normal.bam" > "$OUT_DIR/normal.flagstat.txt"
  samtools flagstat -@ "$INDEX_THREADS" "$OUT_DIR/tumor_t0.bam" > "$OUT_DIR/tumor_t0.flagstat.txt"

  echo "=== Mutect2 restricted to t0 truth regions ==="
  gatk Mutect2 \
    -R "$REF" \
    -I "$OUT_DIR/tumor_t0.bam" \
    -tumor tumor_t0 \
    -I "$OUT_DIR/normal.bam" \
    -normal normal \
    -L "$OUT_DIR/t0_truth.padded.bed" \
    -O "$OUT_DIR/t0.mutect2.unfiltered.vcf.gz"

  gatk FilterMutectCalls \
    -R "$REF" \
    -V "$OUT_DIR/t0.mutect2.unfiltered.vcf.gz" \
    -O "$OUT_DIR/t0.mutect2.filtered.vcf.gz"

  echo "=== Normalize PASS calls and compare to truth ==="
  bcftools view -f PASS "$OUT_DIR/t0.mutect2.filtered.vcf.gz" \
    | bcftools norm -f "$REF" -d exact -Oz \
    -o "$OUT_DIR/t0.called.pass.norm.vcf.gz"
  bcftools index -f "$OUT_DIR/t0.called.pass.norm.vcf.gz"

  rm -rf "$OUT_DIR/isec_truth_vs_called"
  bcftools isec \
    -p "$OUT_DIR/isec_truth_vs_called" \
    "$OUT_DIR/t0_truth.vcf.gz" \
    "$OUT_DIR/t0.called.pass.norm.vcf.gz"

  {
    echo -e "category\tcount"
    echo -e "truth_total\t$(bcftools view -H "$OUT_DIR/t0_truth.vcf.gz" | wc -l)"
    echo -e "called_pass_total\t$(bcftools view -H "$OUT_DIR/t0.called.pass.norm.vcf.gz" | wc -l)"
    echo -e "truth_only\t$(grep -vc '^#' "$OUT_DIR/isec_truth_vs_called/0000.vcf" 2>/dev/null || echo 0)"
    echo -e "called_only\t$(grep -vc '^#' "$OUT_DIR/isec_truth_vs_called/0001.vcf" 2>/dev/null || echo 0)"
    echo -e "overlap_truth_side\t$(grep -vc '^#' "$OUT_DIR/isec_truth_vs_called/0002.vcf" 2>/dev/null || echo 0)"
    echo -e "overlap_called_side\t$(grep -vc '^#' "$OUT_DIR/isec_truth_vs_called/0003.vcf" 2>/dev/null || echo 0)"
  } > "$OUT_DIR/t0_truth_vs_called_summary.tsv"
  cat "$OUT_DIR/t0_truth_vs_called_summary.tsv"
}

set +e
TIME_PAYLOAD="$(
  declare -p OUT_DIR REF NORMAL_R1 NORMAL_R2 TUMOR_R1 TUMOR_R2 \
    BWA_THREADS SORT_THREADS SORT_MEM INDEX_THREADS 2>/dev/null
  declare -f run_validation
  echo "run_validation"
)"
/usr/bin/time -v -o "$TIME_METRICS" bash -c "$TIME_PAYLOAD"
CMD_EXIT=$?
set -e

END_EPOCH=$(date +%s)
{
  echo "exit_code=$CMD_EXIT"
  echo "finished_at=$(date --iso-8601=seconds)"
  echo "elapsed_seconds=$((END_EPOCH - START_EPOCH))"
  echo "source_reference=$SOURCE_REF"
  echo "validation_reference=$REF"
  echo "=== output size ==="
  du -sh "$OUT_DIR" 2>/dev/null || true
  echo "=== validation outputs ==="
  find "$OUT_DIR" -maxdepth 2 \( -type f -o -type l \) -print | sort | xargs -r ls -lh
  echo "=== time metrics file ==="
  cat "$TIME_METRICS" 2>/dev/null || true
} >> "$RUN_SUMMARY"

echo "=== Metrics ==="
echo "Summary: $RUN_SUMMARY"
echo "GNU time: $TIME_METRICS"
exit "$CMD_EXIT"
