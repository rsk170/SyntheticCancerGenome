# SyntheticCancerGenome Patient FASTQ Workflow

This repository is a project-specific fork/adaptation of
[Rbbt-Workflows/SyntheticCancerGenome](https://github.com/Rbbt-Workflows/SyntheticCancerGenome).
The upstream repository provides the Ruby/Rbbt workflow used as the core
simulation backend, including germline generation and NEAT-based FASTQ
simulation. This fork adds patient-level Python and Slurm wrappers to prepare
clone-specific inputs, convert clinical MAF files from `hg19` to `hg38`, enforce
sex-aware references and germlines, generate matched normal/tumour FASTQs, and
validate the generated data.

The workflow below documents the exact path used for the patient:

```text
patients/79ce1d89-46d2-5513-c704-212aa1ed97d2
```

The simulation build is `hg38`. Clone MAFs are treated as `hg19` inputs and are
converted before FASTQ generation.

## 🧭 Workflow Overview

The final workflow is:

1. Prepare the per-patient input folder with clone proportions and clone MAFs.
2. Generate or provide one shared diploid germline per patient using the upstream Ruby/Rbbt workflow.
3. Prepare one manifest per patient/timepoint.
4. Convert clone MAF files from `hg19` to normalized `hg38` VCFs.
5. Convert validated `hg38` clone variants to SyntheticCancerGenome mutation strings.
6. Build cumulative final clone mutation files, adding founding mutations to descendant clones.
7. Create one sex-aware patient reference and one sex-aware germline file.
8. Generate one matched normal FASTQ pair.
9. Generate independent tumour clone FASTQ pairs at the required clone depths.
10. Merge active clone FASTQs into one tumour FASTQ pair per timepoint.
11. Validate FASTQ read identifiers, tumour-normal variant calling, and truth-site read support.

Downsampling scripts are present in the repository history/workspace but are not
part of the final recommended workflow. The final strategy is to simulate clone
FASTQs independently at the required timepoint depths, except when an already
generated clone FASTQ exactly matches the required depth.

## 📦 Software Requirements

The workflow was run on a Slurm-managed HPC environment. The exact module names
used on MareNostrum 5 were:

```bash
module load oneapi hdf5 python/3.12.1
module load samtools/1.19.2
module load htslib/1.19.1
module load bcftools/1.19
module load bwa/0.7.17
module load java-openjdk/17.0.11+9
module load gatk/4.5.0.0
```

Additional required software:

```text
Ruby 3.3.6 with Rbbt/SyntheticCancerGenome installed
NEAT gen_reads.py available through the Rbbt setup
Picard LiftoverVcf or a Picard jar
Slurm
GNU time
```

Python packages used by the full workflow:

```text
pysam
numpy
biopython
matplotlib
matplotlib-venn
pandas
```

The repository includes a minimal `requirements.txt` for the Python packages
used by the NEAT/Rbbt simulation environment.

## 🧬 Required Inputs

User-provided/project inputs for each patient are only the clone composition and
clone mutation files:

```text
patients/PATIENT_ID/clone_proportions.csv - supplied from the clinical-genomic framework 
patients/PATIENT_ID/clone_*.maf - supplied from the clinical-genomic framework
```

The shared germline file below is **not** produced by Synthea and is **not** a
clone MAF input. It is generated once per patient by the upstream
SyntheticCancerGenome Ruby/Rbbt task `diploid_genotype_germline_hg38`:

```text
patients/PATIENT_ID/diploid_genotype_germline_hg38.out
```

The current manifest builder records the germline path, so this file must exist
before running `scripts/pipeline/prepare_patient_manifest.py`. In a fresh run, generate it
with the Ruby/Rbbt command in Step 1. If it was generated previously, copy or
symlink the existing result to the path above.

Reference/liftover inputs:

```text
hg19 FASTA
hg38 FASTA
hg19ToHg38.over.chain.gz
picard.jar or a working picard command
```

The patient used here is female and was simulated with:

```text
normal_target_depth = 30x
tumor_target_depth = 100x
```

Clone fractions:

```text
t0: clone_1 45%, clone_3 55%
t1: clone_1 18%, clone_2 33%, clone_3 49%
t2: clone_1 23%, clone_2 51%, clone_4 26%
```

## 🖥️ Where To Run Each Step

Small preprocessing steps can be run locally or on the cluster login node if the
cluster policy allows light file-processing commands. Long-running or I/O-heavy
steps must be run through Slurm.

Run Slurm jobs from the cluster repository root:

```bash
cd /path/to/SyntheticCancerGenome
```

The Slurm wrappers avoid machine-specific absolute paths. If you submit them
from another directory, set `REPO_ROOT=/path/to/SyntheticCancerGenome`. If your
cluster requires an account or QOS, pass those at submission time, for example
`sbatch --account=YOUR_ACCOUNT --qos=YOUR_QOS ...`. On systems without the same
environment modules, set `LOAD_MODULES=0` and provide the required tools on
`PATH`.

The long-running HPC steps are:

```text
normal FASTQ generation
tumour clone FASTQ generation
tumour FASTQ merging
FASTQ ID validation
alignment and Mutect2 validation
truth-site read-support validation
```

## 📁 Step 0: Prepare The Patient Input Folder

Create a patient folder containing only that patient's clone proportions and MAF
files. If the clone proportions file was exported as a multi-patient cohort
table, subset it first so `patients/$PATIENT/clone_proportions.csv` contains
only rows for that patient.

```bash
PATIENT=79ce1d89-46d2-5513-c704-212aa1ed97d2

mkdir -p patients/$PATIENT
```

Expected files before generating the germline:

```text
patients/$PATIENT/clone_proportions.csv
patients/$PATIENT/clone_1*.maf
patients/$PATIENT/clone_2*.maf
...
```

## 🧫 Step 1: Generate Or Provide The Shared Germline

The patient needs one shared diploid `hg38` germline generated by the upstream
Ruby task `diploid_genotype_germline_hg38`. This task samples germline variants
from the 1000 Genomes Phase 3 VCF, lifts selected positions to `hg38`, validates
them against the reference, and writes a diploid mutation list using the
SyntheticCancerGenome copy-labelled chromosome convention:

```text
copy-1_chr... = one haplotype
copy-2_chr... = the other haplotype
```

Run this from the repository root in the Ruby/Rbbt environment:

```bash
export PATIENT=79ce1d89-46d2-5513-c704-212aa1ed97d2
export RUBY_BIN=/path/to/ruby

GERMLINE_JOB_PATH=$("$RUBY_BIN" -Ilib -rrbbt-util -e '
require "./workflow"
patient = ENV.fetch("PATIENT")
job = SyntheticCancerGenome.job(:diploid_genotype_germline_hg38, patient)
job.produce
puts job.path
')

mkdir -p patients/$PATIENT
cp -f "$GERMLINE_JOB_PATH" patients/$PATIENT/diploid_genotype_germline_hg38.out
```

This command assumes the patched Ruby/Rbbt environment is already working,
including access to `liftOver` and the required cached resources under
`~/.rbbt/share`.

On a cluster, use the Ruby executable from the Rbbt environment, for example:

```bash
export RUBY_BIN=/path/to/rbbt_env/bin/ruby
```

The first run may download/cache upstream resources under `~/.rbbt/share`,
including the 1000 Genomes Phase 3 VCF and the `hg19ToHg38` liftOver chain.

Before continuing, confirm the copied germline exists and contains both
haplotypes:

```bash
ls -lh patients/$PATIENT/diploid_genotype_germline_hg38.out
grep -c '^copy-1_' patients/$PATIENT/diploid_genotype_germline_hg38.out
grep -c '^copy-2_' patients/$PATIENT/diploid_genotype_germline_hg38.out
```

Both `grep` counts should be nonzero.

## ⚙️ Step 2: Prepare Patient Manifests

Create one manifest per timepoint. The manifest records patient sex, target
depths, active clones, clone fractions, clone type, parent clone, paths to the
clone MAF files, and the path to the shared diploid germline generated in Step 1.

```bash
PATIENT=79ce1d89-46d2-5513-c704-212aa1ed97d2

python3 scripts/pipeline/prepare_patient_manifest.py patients/$PATIENT \
  --sex female \
  --normal-depth 30 \
  --tumor-depth 100 \
  --timepoint t0

python3 scripts/pipeline/prepare_patient_manifest.py patients/$PATIENT \
  --sex female \
  --normal-depth 30 \
  --tumor-depth 100 \
  --timepoint t1

python3 scripts/pipeline/prepare_patient_manifest.py patients/$PATIENT \
  --sex female \
  --normal-depth 30 \
  --tumor-depth 100 \
  --timepoint t2
```

Outputs:

```text
patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv
patients/$PATIENT/prepared_hg38_t1/patient_manifest.csv
patients/$PATIENT/prepared_hg38_t2/patient_manifest.csv
```

## 🔁 Step 3: Convert Clone MAFs From hg19 To hg38

This step converts active clone MAF variants from `hg19` to `hg38` using a
VCF-aware path: MAF to VCF-like records, normalization against `hg19`, Picard
liftover to `hg38`, normalization against `hg38`, splitting and filtering.

Example:

```bash
python3 scripts/pipeline/convert_patient_mafs_hg19_to_hg38.py \
  patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv \
  --hg19-fasta /path/to/hg19.fa.gz \
  --hg38-fasta /path/to/hg38.fa.gz \
  --chain /path/to/hg19ToHg38.over.chain.gz \
  --picard-jar /path/to/picard.jar
```

Repeat for `t1` and `t2`:

```bash
python3 scripts/pipeline/convert_patient_mafs_hg19_to_hg38.py \
  patients/$PATIENT/prepared_hg38_t1/patient_manifest.csv \
  --hg19-fasta /path/to/hg19.fa.gz \
  --hg38-fasta /path/to/hg38.fa.gz \
  --chain /path/to/hg19ToHg38.over.chain.gz \
  --picard-jar /path/to/picard.jar

python3 scripts/pipeline/convert_patient_mafs_hg19_to_hg38.py \
  patients/$PATIENT/prepared_hg38_t2/patient_manifest.csv \
  --hg19-fasta /path/to/hg19.fa.gz \
  --hg38-fasta /path/to/hg38.fa.gz \
  --chain /path/to/hg19ToHg38.over.chain.gz \
  --picard-jar /path/to/picard.jar
```

Outputs per timepoint:

```text
prepared_hg38_t*/hg19_vcf/
prepared_hg38_t*/hg38_vcf/
prepared_hg38_t*/hg38_mutations/
prepared_hg38_t*/rejected_variants/
prepared_hg38_t*/conversion_summary.tsv
```

## 🧩 Step 4: Build Final Cumulative Clone Mutation Files

The clone MAF files are treated as non-cumulative. The founding clone contains
founding mutations, while descendant clones receive the founding mutations plus
their own clone-specific variants.

```bash
python3 scripts/pipeline/build_final_clone_mutations.py \
  patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv

python3 scripts/pipeline/build_final_clone_mutations.py \
  patients/$PATIENT/prepared_hg38_t1/patient_manifest.csv

python3 scripts/pipeline/build_final_clone_mutations.py \
  patients/$PATIENT/prepared_hg38_t2/patient_manifest.csv
```

Outputs:

```text
patients/$PATIENT/prepared_hg38_t*/final_clone_mutations/
patients/$PATIENT/prepared_hg38_t*/final_clone_mutations/patient_manifest.final_clone_mutations.csv
patients/$PATIENT/prepared_hg38_t*/final_clone_mutations/final_clone_mutation_summary.tsv
```

## 🚻 Step 5: Create Sex-Aware Reference And Germline

Create one sex-aware patient reference:

```bash
python3 scripts/pipeline/create_sex_aware_reference.py \
  patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv \
  --reference-fasta /path/to/hg38.fa.gz
```

Create one sex-aware patient germline:

```bash
python3 scripts/pipeline/create_sex_aware_germline.py \
  patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv \
  --strict-diploid
```

For female patients, autosomes and chromosome `X` are retained as diploid, `Y`
is removed and mitochondrial sequence is retained as haploid. For male patients,
autosomes are diploid and `X`, `Y` and mitochondrial sequence are haploid.

Outputs:

```text
patients/$PATIENT/sex_aware_reference/
patients/$PATIENT/sex_aware_germline/
```

## 🚀 Step 6: Copy Prepared Files To The HPC

If Steps 1 to 5 were run locally, copy the prepared patient folder or changed
subfolders to the cluster before starting Slurm jobs.

Example:

```bash
REMOTE_USER=your_username
REMOTE_HOST=your_transfer_host
REMOTE_REPO=/path/to/SyntheticCancerGenome

rsync -av patients/$PATIENT/prepared_hg38_t0 \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/

rsync -av patients/$PATIENT/prepared_hg38_t1 \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/

rsync -av patients/$PATIENT/prepared_hg38_t2 \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/
```

From this point onward, run jobs from the cluster repository root:

```bash
cd /path/to/SyntheticCancerGenome
```

## 🧪 Step 7: Generate The Matched Normal FASTQs

Submit the normal FASTQ job through Slurm:

```bash
sbatch scripts/slurm/generate_normal_fastq.sh
```

Expected outputs:

```text
patients/$PATIENT/normal_fastq/normal_read1.fq.gz
patients/$PATIENT/normal_fastq/normal_read2.fq.gz
```

## 🧬 Step 8: Generate Tumour Clone FASTQs

The final workflow generates independent clone FASTQ pairs at the depth needed
for each timepoint. The generic Slurm wrapper is:

```text
scripts/slurm/generate_timepoint_clone.sh
```

Usage:

```bash
sbatch scripts/slurm/generate_timepoint_clone.sh TIMEPOINT CLONE_ID
```

Examples:

```bash
sbatch scripts/slurm/generate_timepoint_clone.sh t1 clone_1
sbatch scripts/slurm/generate_timepoint_clone.sh t1 clone_2
sbatch scripts/slurm/generate_timepoint_clone.sh t1 clone_3
sbatch scripts/slurm/generate_timepoint_clone.sh t2 clone_1
```

The wrapper reads:

```text
patients/$PATIENT/prepared_hg38_TIMEPOINT/final_clone_mutations/patient_manifest.final_clone_mutations.csv
```

and writes:

```text
patients/$PATIENT/tumor_clone_fastqs_independent/TIMEPOINT/CLONE_ID/
```

For this patient, existing exact-depth clone FASTQs can be reused:

```text
t0 clone_1 45x
t0 clone_3 55x
t2 clone_2 51x
t2 clone_4 26x
```

For organization, these exact-depth clone folders can be symlinked into the
timepoint-specific independent FASTQ folder. Example for `t2`:

```bash
BASE=patients/$PATIENT
mkdir -p "$BASE/tumor_clone_fastqs_independent/t2"

ln -s ../../tumor_clone_fastqs_max_depth/clone_2 \
  "$BASE/tumor_clone_fastqs_independent/t2/clone_2"

ln -s ../../tumor_clone_fastqs_max_depth/clone_4 \
  "$BASE/tumor_clone_fastqs_independent/t2/clone_4"
```

## 🧵 Step 9: Merge Clone FASTQs Into Timepoint Tumours

After all active clone FASTQs for a timepoint exist, merge them into a single
tumour FASTQ pair.

Submit the merge through Slurm:

```bash
sbatch scripts/slurm/merge_timepoint_tumor.sh t0
sbatch scripts/slurm/merge_timepoint_tumor.sh t1
sbatch scripts/slurm/merge_timepoint_tumor.sh t2
```

The wrapper calls the underlying merge script. A direct command for `t0` is:

```bash
python scripts/pipeline/merge_patient_tumor_fastqs.py \
  patients/$PATIENT/prepared_hg38_t0/final_clone_mutations/patient_manifest.final_clone_mutations.csv \
  --clone-fastq-dir patients/$PATIENT/tumor_clone_fastqs_independent/t0 \
  --out-dir patients/$PATIENT/tumor_fastq_t0 \
  --output-prefix tumor \
  --overwrite
```

Use the corresponding manifest, clone FASTQ directory and output directory for
`t1` and `t2`.

Primary outputs:

```text
patients/$PATIENT/normal_fastq/normal_read1.fq.gz
patients/$PATIENT/normal_fastq/normal_read2.fq.gz
patients/$PATIENT/tumor_fastq_t0/tumor_read1.fq.gz
patients/$PATIENT/tumor_fastq_t0/tumor_read2.fq.gz
```

## ✅ Step 10: Validate FASTQ Read Identifiers

The read-ID validation checks that:

```text
R1 and R2 are in the same order
each read pair has matching identifiers
merged tumour read identifiers are unique
clone labels in read names match the expected clone mixture
```

For `t0`:

```bash
sbatch scripts/validation/validate_t0_fastq_ids.sh
```

Important outputs:

```text
patients/$PATIENT/validation_t0/fastq_id_validation/t0_pair_order_report.json
patients/$PATIENT/validation_t0/fastq_id_validation/t0_fastq_id_validation_summary.txt
patients/$PATIENT/validation_t0/fastq_id_validation/t0_duplicate_r1_base_id_count.txt
```

Expected successful result:

```text
status=PASS
mate_id_mismatches=0
malformed_records=0
duplicate_id_count=0
```

## 🔬 Optional Step 11: Validate Variants With Alignment And Mutect2

This validation aligns the normal and merged tumour FASTQs to standard `hg38`,
runs matched tumour-normal Mutect2 on padded truth-variant intervals, and
compares PASS calls against the expected truth set.

```bash
sbatch --export=ALL,HG38_REF=patients/$PATIENT/prepared_hg38_t1/ref_cache/hg38/hg38.fa.gz \
  scripts/validation/validate_t0_variant_calling.sh
```

Important outputs:

```text
patients/$PATIENT/validation_t0/normal.bam
patients/$PATIENT/validation_t0/tumor_t0.bam
patients/$PATIENT/validation_t0/t0_truth.vcf.gz
patients/$PATIENT/validation_t0/t0.mutect2.filtered.vcf.gz
patients/$PATIENT/validation_t0/t0.called.pass.norm.vcf.gz
patients/$PATIENT/validation_t0/t0_truth_vs_called_summary.tsv
```

For the current `t0` validation:

```text
truth_total = 935
called_pass_total = 439
truth variants recovered as PASS = 422
PASS calls outside truth set = 17
```

## 🧾 Optional Step 12: Check Truth-Site Read Support

This follow-up checks whether expected truth variants have direct tumour
alternate-allele read support in the aligned BAM files.

```bash
sbatch scripts/validation/validate_t0_truth_read_support.sh
```

Important outputs:

```text
patients/$PATIENT/validation_t0/read_support/t0_truth_variant_read_support.tsv
patients/$PATIENT/validation_t0/read_support/t0_truth_only_variant_read_support.tsv
patients/$PATIENT/validation_t0/read_support/t0_truth_variant_read_support_summary.tsv
```

For the current `t0` validation:

```text
truth_total = 935
truth variants with tumour ALT support = 449
Mutect2-missed truth variants with tumour ALT support = 27
Mutect2-missed truth variants without tumour ALT support = 486
```

## 📊 Computational Metrics

The Slurm wrappers write run summaries and GNU time output under:

```text
patients/$PATIENT/run_metrics/
```

For clone-level tumour FASTQ generation, the final workflow used:

```text
1 Slurm task
64 requested CPUs
125 GB requested memory
8 NEAT chromosome-copy workers
16 samtools threads for final merging
```

On MareNostrum 5, Slurm allocated one node for these jobs. In the recorded
runs, that corresponded to 128 allocated CPUs and 125 GB memory.

Example recorded outputs:

```text
t1 clone_1, 18x: 24 h 32 min, 3.3 GB peak RSS, 44 GB compressed FASTQ output
t1 clone_2, 33x: 35 h 43 min, 3.3 GB peak RSS, 80 GB compressed FASTQ output
t0 merged tumour, 100x: 280 GB compressed FASTQ output
```

## 🗂️ Important Scripts

Patient preparation:

```text
scripts/pipeline/prepare_patient_manifest.py
scripts/pipeline/convert_patient_mafs_hg19_to_hg38.py
scripts/pipeline/build_final_clone_mutations.py
scripts/pipeline/create_sex_aware_reference.py
scripts/pipeline/create_sex_aware_germline.py
```

FASTQ generation and merging:

```text
scripts/slurm/generate_normal_fastq.sh
scripts/slurm/generate_timepoint_clone.sh
scripts/slurm/merge_timepoint_tumor.sh
scripts/pipeline/generate_patient_normal_fastqs.py
scripts/pipeline/generate_patient_clone_tumor_fastqs.py
scripts/pipeline/merge_patient_tumor_fastqs.py
```

Validation:

```text
scripts/validation/validate_t0_fastq_ids.sh
scripts/validation/validate_fastq_pair_ids.py
scripts/validation/validate_t0_variant_calling.sh
scripts/validation/validate_t0_truth_read_support.sh
scripts/validation/check_truth_variant_read_support.py
```

## 🧪 Quick Sanity Checks

Check generated clone FASTQs:

```bash
find patients/$PATIENT/tumor_clone_fastqs_independent \
  -maxdepth 3 \( -name "*_read1.fq.gz" -o -name "*_read2.fq.gz" -o -name "*.tumor_fastq.json" \) \
  | sort
```

Check final merged tumour FASTQs:

```bash
ls -lh patients/$PATIENT/tumor_fastq_t0/
zcat patients/$PATIENT/tumor_fastq_t0/tumor_read1.fq.gz | head -4
zcat patients/$PATIENT/tumor_fastq_t0/tumor_read2.fq.gz | head -4
```

Check run status:

```bash
find patients/$PATIENT/run_metrics -name "*.summary.txt" -print | sort
grep -H 'exit_code=' patients/$PATIENT/run_metrics/*.summary.txt
```

## 📝 Notes

- The patient-specific sex-aware reference is used for simulation.
- The standard `hg38` reference is used for downstream validation alignment and
  Mutect2 calling.
- Clone mutation files are cumulative after `build_final_clone_mutations.py`.
- FASTQ files are very large. Avoid running I/O-heavy checks on login nodes.
- Use Slurm for FASTQ generation, merging and validation.
