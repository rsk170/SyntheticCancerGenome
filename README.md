# SyntheticCancerGenome Patient FASTQ Workflow

This repository is a project-specific fork/adaptation of
[Rbbt-Workflows/SyntheticCancerGenome](https://github.com/Rbbt-Workflows/SyntheticCancerGenome).
The upstream repository provides the Ruby/Rbbt workflow used as the core
simulation backend, including germline generation and NEAT-based FASTQ
simulation. This fork adds patient-level Python and Slurm wrappers to prepare
clone-specific inputs, convert clinical MAF files from `hg19` to `hg38`, enforce
sex-aware references and germlines, generate matched normal/tumour FASTQs, and
validate the generated data.

The workflow operates on one patient directory at a time. Set `PATIENT` to the
directory name used under `patients/`. The simulation build is `hg38`. Clone MAFs are treated
as `hg19` inputs and are converted before FASTQ generation.

## Custom Extensions In This Fork

The upstream Ruby/Rbbt simulation backend is retained. The main additions and
modifications maintained in this fork are:

| Area | Location | Extension |
| --- | --- | --- |
| Germline task compatibility | `lib/synthetic_cancer_genome/tasks/germline.rb` | Retains the upstream SyntheticCancerGenome germline-generation workflow, with fork-specific changes to stream variants directly from the 1000 Genomes Phase 3 VCF and perform hg19-to-hg38 liftover in the tested environment. |
| Patient workflow | `scripts/pipeline/` | Adds patient manifests, MAF liftover, cumulative clone mutation lists, sex-aware reference/germline preparation, FASTQ generation wrappers and timepoint FASTQ merging. |
| HPC execution | `scripts/slurm/` | Adds portable Slurm wrappers and computational-resource recording for normal, clone and merged-tumour generation. |
| Validation | `scripts/validation/` | Adds FASTQ pair checks, matched tumour-normal variant calling and truth-site read-support analysis. |
| Example inputs | `examples/patient_workflow/` | Provides a lightweight patient example with clone proportions, clone MAFs and checksums. |
| Compatibility patches | `patches/` | Records the Rbbt HTS reference changes and NEAT Python 3/Biopython compatibility changes required by the tested environment. |
| Environment and tests | `test/python/`, `requirements.txt` | Adds lightweight patient-workflow tests and pinned Python requirements. |

## 🧭 Workflow Overview

The workflow is:

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

The clone FASTQs are simulated independently at the required timepoint depths.

## 📦 Software Requirements

The workflow was tested on Linux in a Slurm-managed HPC environment. The exact
module names used on MareNostrum 5 were:

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

Slurm is required for the supported execution path documented below. The
underlying Python entry points can be invoked directly for development and
small tests, but whole-genome FASTQ generation, merging and validation are
documented and tested as Slurm jobs.

### Rbbt and NEAT dependencies

Rbbt, HTS, NEATGenReads and NEAT are third-party projects and are not vendored
in this repository. Users should follow their upstream installation guidance.
This repository documents the exact integration boundary and tested revisions
needed by this workflow rather than duplicating their full documentation.

The development environment used these versions:

| Component | Version or revision |
| --- | --- |
| Ruby | 3.3.6 |
| `rbbt-util` | 6.0.5 |
| `rbbt-sources` | 3.4.2 |
| [Rbbt HTS](https://github.com/Rbbt-Workflows/HTS) | `5c11b646749076e63204e296ef44674e84234260` |
| [Rbbt NEATGenReads](https://github.com/Rbbt-Workflows/NEATGenReads) | `000b1cf380d668e66beb8fdaa42140a55bddda42` |
| [legacy NEAT simulator](https://github.com/zstephens/neat-genreads) | `a2d7739c9102712f277b99d381308c99d52907f6` |
| Picard | 2.20.6 |

One tested installation layout is:

```bash
gem install rbbt-util -v 6.0.5
gem install rbbt-sources -v 3.4.2

mkdir -p "$HOME/.rbbt/workflows" "$HOME/.rbbt/software/opt/bin"

git clone https://github.com/Rbbt-Workflows/HTS.git \
  "$HOME/.rbbt/workflows/HTS"
git -C "$HOME/.rbbt/workflows/HTS" checkout \
  5c11b646749076e63204e296ef44674e84234260
git -C "$HOME/.rbbt/workflows/HTS" apply \
  "$REPO_ROOT/patches/rbbt-hts-ucsc-reference.patch"

git clone https://github.com/Rbbt-Workflows/NEATGenReads.git \
  "$HOME/.rbbt/workflows/NEAT_gen_reads"
git -C "$HOME/.rbbt/workflows/NEAT_gen_reads" checkout \
  000b1cf380d668e66beb8fdaa42140a55bddda42

git clone https://github.com/zstephens/neat-genreads.git \
  "$HOME/.rbbt/software/opt/NEATGenReads"
```

Use a Ruby environment rather than system-wide gems when required by local
policy. The commands above describe the expected layout. The upstream projects may
have additional platform-specific prerequisites.

The legacy NEAT revision requires small Python 3 and modern-Biopython
compatibility changes. After installing the simulator at the revision above,
apply the included patch:

```bash
export REPO_ROOT=/path/to/SyntheticCancerGenome
export NEAT_DIR=$HOME/.rbbt/software/opt/NEATGenReads

git -C "$NEAT_DIR" checkout a2d7739c9102712f277b99d381308c99d52907f6
git -C "$NEAT_DIR" apply \
  "$REPO_ROOT/patches/neat-genreads-python3-biopython.patch"

ln -sfn ../NEATGenReads/gen_reads.py \
  "$HOME/.rbbt/software/opt/bin/gen_reads.py"
export PATH="$HOME/.rbbt/software/opt/bin:$PATH"
```

The Rbbt environment must then make the `HTS` and `NEATGenReads` workflows
available to `Workflow.require_workflow`, and `gen_reads.py` must be on `PATH`.

The module names above are site-specific. On another cluster, load equivalent
software yourself and submit jobs with `LOAD_MODULES=0`. Install the Python
dependencies into the Python environment used by the jobs:

```bash
python3 -m pip install -r requirements.txt
```

The legacy NEAT interface wrapped here is the `gen_reads.py` interface used by
the Rbbt workflow, not the command-line interface of newer NEAT releases.

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
used by the NEAT/Rbbt simulation environment. It pins a reproducible public
Python environment compatible with the included NEAT patch. It should not be
interpreted as a historical `pip freeze` of every cluster package. It does not
replace a small end-to-end test on the target cluster.

### Reproducibility boundary

The patient manifests, references, mutation lists and validation truth sets are
recorded explicitly. However, this version does not yet expose one seed that
controls both germline sampling and every NEAT/Rbbt random operation. A fresh
uncached simulation can therefore produce different reads. Preserve the input
manifests, generated germline, tool versions and Rbbt job outputs for an exact
reported dataset.

## 🧬 Required Inputs

User-provided/project inputs for each patient are the clone composition and
clone mutation files:

```text
patients/PATIENT_ID/clone_proportions.csv - supplied from the clinical-genomic framework
patients/PATIENT_ID/clone_*.maf - supplied from the clinical-genomic framework
```

These files can be generated by the
[Synthea v3.3 clinical-genomic pipeline](https://github.com/rsk170/synthea-v3.3-genomics-pipeline),
which produces `clone_proportions.csv` and the clone-level MAF files consumed
by this repository.

`clone_proportions.csv` must contain these columns:

```text
patient_id,clone_id,clone_type,parent_clone_id,t0_vaf_pct,t1_vaf_pct,...
```

Timepoint columns must use the pattern `tN_vaf_pct`. `clone_type` should identify
one founding clone, and `parent_clone_id` must describe the clone lineage. Each
clone requires one MAF whose filename begins with its clone ID, such as
`clone_1_founding.maf`. MAF files must be tab-delimited and contain at least:

```text
Chromosome
Start_position
End_position
Reference_Allele
Tumor_Seq_Allele1
Tumor_Seq_Allele2
Variant_Type
```

The MAF coordinates must be `hg19`. The `hg19` and `hg38` FASTAs should use
compatible primary chromosome assemblies. The `hg38` FASTA supplied during
conversion must also be used to build the sex-aware simulation reference and
for validation alignment.

The number of sequencing timepoints and clones is patient-specific. Run only
the timepoints present in `clone_proportions.csv` and submit clone-generation
jobs only for clones with a positive fraction at that timepoint. During manifest
preparation, use `--include-zero-clones` so inactive founding clones remain
available when cumulative descendant mutation lists are constructed. Zero-fraction clones are retained only for lineage bookkeeping. FASTQ generation, tumour FASTQ merging and validation use only clones with positive timepoint fractions.

The shared germline file below is generated once per patient by the upstream
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

The pipeline was developed with the UCSC files below. Use the listed URLs and
checksums because similarly named reference files are not always interchangeable.

| Purpose | UCSC assembly | NCBI assembly accession | Exact archive file | Source | MD5 |
| --- | --- | --- | --- | --- | --- |
| Source coordinates | Human Feb. 2009 (`GRCh37/hg19`) | `GCA_000001405.1` | `hg19.fa.gz` | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz` | `806c02398f5ac5da8ffd6da2d1d5d1a9` |
| Simulation coordinates | Human Dec. 2013 (`GRCh38/hg38`) | `GCA_000001405.15` | `hg38.fa.gz` | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` | `1c9dcaddfa41027f17cd8f7a82c7293b` |
| Coordinate conversion | `GRCh37/hg19` to `GRCh38/hg38` | `GCA_000001405.1` to `GCA_000001405.15` | `hg19ToHg38.over.chain.gz` | `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz` | `35887f73fe5e2231656504d1f6430900` |

Verify the MD5 values after download. MAF liftover used Picard 2.20.6.

The walkthrough uses the bundled patient example under
`examples/patient_workflow/79ce1d89-46d2-5513-c704-212aa1ed97d2/` as a female
patient with:

```text
normal_target_depth = 30x
tumor_target_depth = 60x
```

Clone fractions:

```text
t0: clone_1 45%, clone_3 55%
t1: clone_1 18%, clone_2 33%, clone_3 49%
t2: clone_1 23%, clone_2 51%, clone_4 26%
```

These values are an example, not defaults required by the scripts.

## 🖥️ Where To Run Each Step

Small preprocessing steps can be run locally or on the cluster login node. The supported full
workflow requires Slurm. Long-running or I/O-heavy steps must be submitted as
Slurm jobs.

Run Slurm jobs from the cluster repository root:

```bash
cd /path/to/SyntheticCancerGenome
```

The Slurm wrappers avoid machine-specific absolute paths. Configure the patient
and repository once in the submission shell:

```bash
export PATIENT=PATIENT_ID
export REPO_ROOT=/path/to/SyntheticCancerGenome
export RBBT_ENV=/path/to/rbbt_environment
export RUBY_BIN="$RBBT_ENV/bin/ruby"
export LOAD_MODULES=1
```

`PATIENT` is required - jobs stop immediately if it is missing. On systems without the same
environment modules, provide the required tools
on `PATH`. All exported variables are propagated by the documented `sbatch`
commands.

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
export PATIENT=PATIENT_ID

mkdir -p patients/$PATIENT
```

Expected files before generating the germline:

```text
patients/$PATIENT/clone_proportions.csv
patients/$PATIENT/clone_1*.maf
patients/$PATIENT/clone_2*.maf
...
```

### Bundled lightweight example

The repository includes an example clone proportions and four
clone-level MAF inputs under:

```text
examples/patient_workflow/79ce1d89-46d2-5513-c704-212aa1ed97d2/
```

The clone composition is synthetic, and the MAF rows were derived from the
public ICGC mutation resource used by the
[Synthea v3.3 clinical-genomic pipeline](https://github.com/rsk170/synthea-v3.3-genomics-pipeline). Generated germline files,
references, intermediate files, and FASTQs are intentionally excluded. To use
these inputs as the walkthrough patient:

```bash
export PATIENT=79ce1d89-46d2-5513-c704-212aa1ed97d2
export EXAMPLE_DIR="examples/patient_workflow/$PATIENT"

(cd "$EXAMPLE_DIR" && sha256sum -c SHA256SUMS)

mkdir -p "patients/$PATIENT"
cp "$EXAMPLE_DIR/clone_proportions.csv" \
  "$EXAMPLE_DIR"/clone_*.maf \
  "patients/$PATIENT/"
```

Continue with Step 1 to generate the patient's shared germline. The example
does not bypass any preparation, simulation, or validation step.

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
export PATIENT=PATIENT_ID
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
export PATIENT=PATIENT_ID

python3 scripts/pipeline/prepare_patient_manifest.py patients/$PATIENT \
  --sex female \
  --normal-depth 30 \
  --tumor-depth 60 \
  --timepoint t0 \
  --include-zero-clones

python3 scripts/pipeline/prepare_patient_manifest.py patients/$PATIENT \
  --sex female \
  --normal-depth 30 \
  --tumor-depth 60 \
  --timepoint t1 \
  --include-zero-clones

python3 scripts/pipeline/prepare_patient_manifest.py patients/$PATIENT \
  --sex female \
  --normal-depth 30 \
  --tumor-depth 60 \
  --timepoint t2 \
  --include-zero-clones
```

Create manifests only for timepoints available for that patient. For example,
a patient with only `t0` and `t1` does not require a `t2` manifest or any `t2`
generation, merging or validation jobs.

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

The generated copy-labelled FASTA is written in uppercase before NEAT-based
FASTQ simulation. This prevents soft-masked lowercase bases from the source
`hg38` FASTA from propagating into mutation application.

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

If Steps 1 to 4 were run locally, copy the original patient inputs, generated
diploid germline and prepared timepoint directories to the cluster, then rerun
Step 5 on the cluster. Do not rely on locally generated sex-aware metadata when
the repository path changes between machines.

Example:

```bash
REMOTE_USER=your_username
REMOTE_HOST=your_transfer_host
REMOTE_REPO=/path/to/SyntheticCancerGenome

ssh "$REMOTE_USER@$REMOTE_HOST" \
  "mkdir -p '$REMOTE_REPO/patients/$PATIENT'"

rsync -av patients/$PATIENT/prepared_hg38_t0 \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/

rsync -av patients/$PATIENT/prepared_hg38_t1 \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/

rsync -av patients/$PATIENT/prepared_hg38_t2 \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/

rsync -av patients/$PATIENT/clone_proportions.csv \
  patients/$PATIENT/clone_*.maf \
  patients/$PATIENT/diploid_genotype_germline_hg38.out \
  $REMOTE_USER@$REMOTE_HOST:$REMOTE_REPO/patients/$PATIENT/
```

From this point onward, run commands from the cluster repository root and
recreate the sex-aware files there:

```bash
cd /path/to/SyntheticCancerGenome

python3 scripts/pipeline/create_sex_aware_reference.py \
  patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv \
  --reference-fasta /cluster/path/to/hg38.fa.gz

python3 scripts/pipeline/create_sex_aware_germline.py \
  patients/$PATIENT/prepared_hg38_t0/patient_manifest.csv \
  --strict-diploid
```

If Steps 1 to 5 were all run directly on the cluster, no transfer or metadata
recreation is needed.

## 🧪 Step 7: Generate The Matched Normal FASTQs

Submit the normal FASTQ job through Slurm:

```bash
sbatch scripts/slurm/generate_normal_fastq.sh
```

The wrapper disables the optional Rbbt read-renaming pass. The original NEAT
read identifiers are retained, avoiding an additional whole-FASTQ rewrite.
By default, final FASTQ paths under `patients/` are symbolic links to the Rbbt
job cache. This avoids duplicating hundreds of gigabytes, but the patient folder
is not self-contained: deleting or moving `~/.rbbt/var/jobs` breaks those links.
Set `MATERIALIZE=copy` at submission when independent files are required for
archival, transfer or long-term retention - this needs substantially more
storage and additional copy time. `MATERIALIZE=hardlink` avoids duplicate disk
usage and survives deletion of the cache filename, but works only when both
locations are on the same filesystem.

Expected outputs:

```text
patients/$PATIENT/normal_fastq/normal_read1.fq.gz
patients/$PATIENT/normal_fastq/normal_read2.fq.gz
```

## 🧬 Step 8: Generate Tumour Clone FASTQs

The final workflow generates independent clone FASTQ pairs at the depth needed
for each timepoint. Clone-level depth is calculated from the timepoint clone
fraction and the total tumour target depth recorded in the manifest. The generic
Slurm wrapper is:

```text
scripts/slurm/generate_timepoint_clone.sh
```

Usage:

```bash
sbatch scripts/slurm/generate_timepoint_clone.sh TIMEPOINT CLONE_ID
```

Before a long run, inspect the resolved reference, germline, mutation file and
depth without generating reads:

```bash
sbatch scripts/slurm/generate_timepoint_clone.sh t0 clone_1 --dry-run
```

Examples:

```bash
sbatch scripts/slurm/generate_timepoint_clone.sh t0 clone_1
sbatch scripts/slurm/generate_timepoint_clone.sh t0 clone_3
sbatch scripts/slurm/generate_timepoint_clone.sh t1 clone_1
sbatch scripts/slurm/generate_timepoint_clone.sh t1 clone_2
sbatch scripts/slurm/generate_timepoint_clone.sh t1 clone_3
sbatch scripts/slurm/generate_timepoint_clone.sh t2 clone_1
sbatch scripts/slurm/generate_timepoint_clone.sh t2 clone_2
sbatch scripts/slurm/generate_timepoint_clone.sh t2 clone_4
```

The wrapper reads:

```text
patients/$PATIENT/prepared_hg38_TIMEPOINT/final_clone_mutations/patient_manifest.final_clone_mutations.csv
```

and writes:

```text
patients/$PATIENT/tumor_clone_fastqs_independent/TIMEPOINT/CLONE_ID/
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
python3 scripts/pipeline/merge_patient_tumor_fastqs.py \
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

Run the same validation for every merged timepoint:

```bash
sbatch scripts/validation/validate_timepoint_fastq_ids.sh t0
sbatch scripts/validation/validate_timepoint_fastq_ids.sh t1
sbatch scripts/validation/validate_timepoint_fastq_ids.sh t2
```

Important outputs:

```text
patients/$PATIENT/validation_TIMEPOINT/fastq_id_validation/TIMEPOINT_pair_order_report.json
patients/$PATIENT/validation_TIMEPOINT/fastq_id_validation/TIMEPOINT_fastq_id_validation_summary.txt
patients/$PATIENT/validation_TIMEPOINT/fastq_id_validation/TIMEPOINT_duplicate_r1_base_id_count.txt
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

Set `HG38_REF` to the same standard `hg38` FASTA supplied to the MAF conversion
and sex-aware reference steps. Do not use the copy-labelled sex-aware simulation
reference for alignment.

```bash
export HG38_REF=/path/to/hg38.fa

sbatch --export=ALL,HG38_REF="$HG38_REF" \
  scripts/validation/validate_timepoint_variant_calling.sh t0
```

After `t0` succeeds, reuse its normal BAM for later timepoints instead of
aligning the same shared normal FASTQs again:

```bash
export NORMAL_BAM="patients/$PATIENT/validation_t0/normal.bam"

sbatch --export=ALL,HG38_REF="$HG38_REF",NORMAL_BAM="$NORMAL_BAM" \
  scripts/validation/validate_timepoint_variant_calling.sh t1

sbatch --export=ALL,HG38_REF="$HG38_REF",NORMAL_BAM="$NORMAL_BAM" \
  scripts/validation/validate_timepoint_variant_calling.sh t2
```

Submit later timepoints only after the `t0` job has completed successfully, or
use a Slurm `afterok` dependency. Validation writes BAMs through temporary files
and only promotes them after `samtools quickcheck` succeeds. If a previous
version left a nonempty but corrupt BAM, remove that failed validation directory
or explicitly replace the BAM before restarting.

Important outputs:

```text
patients/$PATIENT/validation_TIMEPOINT/normal.bam
patients/$PATIENT/validation_TIMEPOINT/tumor_TIMEPOINT.bam
patients/$PATIENT/validation_TIMEPOINT/TIMEPOINT_truth.vcf.gz
patients/$PATIENT/validation_TIMEPOINT/TIMEPOINT.mutect2.filtered.vcf.gz
patients/$PATIENT/validation_TIMEPOINT/TIMEPOINT.called.pass.norm.vcf.gz
patients/$PATIENT/validation_TIMEPOINT/TIMEPOINT_truth_vs_called_summary.tsv
```

## 🧾 Optional Step 12: Check Truth-Site Read Support

This follow-up checks whether expected truth variants have direct tumour
alternate-allele read support in the aligned BAM files.

```bash
sbatch scripts/validation/validate_timepoint_truth_read_support.sh t0
sbatch scripts/validation/validate_timepoint_truth_read_support.sh t1
sbatch scripts/validation/validate_timepoint_truth_read_support.sh t2
```

Important outputs:

```text
patients/$PATIENT/validation_TIMEPOINT/read_support/TIMEPOINT_truth_variant_read_support.tsv
patients/$PATIENT/validation_TIMEPOINT/read_support/TIMEPOINT_truth_only_variant_read_support.tsv
patients/$PATIENT/validation_TIMEPOINT/read_support/TIMEPOINT_truth_variant_read_support_summary.tsv
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
8 NEAT chromosome-copy workers
16 samtools threads for final merging
```

On MareNostrum 5, memory was assigned automatically from the CPU request under
the GPP partition policy. Slurm allocated one node for these jobs. In the
recorded runs this was reported as 128 allocated CPUs and approximately 125 GB
memory.

Recorded clone-level tumour FASTQ generation metrics for the bundled example
patient were:

| Time point | Clone | Depth | Wall-clock time | Peak RAM (GNU time) | R1 FASTQ | R2 FASTQ | Total FASTQ size | User CPU time | System CPU time | Avg. CPU use |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| t0 | clone_1 | 27x | 31 h 57 min 42 s | 3.3 GB | 33 GB | 33 GB | 66 GB | 594,680 s | 23,789 s | 537% |
| t0 | clone_3 | 33x | 34 h 32 min 54 s | 3.3 GB | 40 GB | 40 GB | 80 GB | 626,396 s | 20,805 s | 520% |
| t1 | clone_1 | 11x | 20 h 15 min 43 s | 3.3 GB | 15 GB | 15 GB | 30 GB | 404,542 s | 20,201 s | 582% |
| t1 | clone_2 | 20x | 25 h 27 min 48 s | 3.3 GB | 23 GB | 23 GB | 46 GB | 487,439 s | 20,786 s | 554% |
| t1 | clone_3 | 30x | 31 h 53 min 51 s | 3.3 GB | 35 GB | 35 GB | 70 GB | 586,817 s | 20,881 s | 529% |
| t2 | clone_1 | 14x | 21 h 30 min 12 s | 3.3 GB | 16 GB | 16 GB | 32 GB | 425,084 s | 20,326 s | 575% |
| t2 | clone_2 | 31x | 33 h 13 min 38 s | 3.3 GB | 38 GB | 38 GB | 76 GB | 606,367 s | 21,229 s | 524% |
| t2 | clone_4 | 16x | 23 h 05 min 39 s | 3.3 GB | 19 GB | 19 GB | 38 GB | 451,570 s | 21,195 s | 586% |

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
scripts/validation/validate_timepoint_fastq_ids.sh
scripts/validation/validate_fastq_pair_ids.py
scripts/validation/validate_timepoint_variant_calling.sh
scripts/validation/validate_timepoint_truth_read_support.sh
scripts/validation/check_truth_variant_read_support.py
```

## 🧪 Quick Sanity Checks

Run the lightweight tests for the bundled example inputs and manifest logic:

```bash
python3 -m unittest discover -s test/python -p 'test_*.py' -v
```

These tests require only the Python standard library. They do not generate a
germline, download references or run FASTQ simulation.

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
- The sex-aware reference is written in uppercase before simulation.
- The standard `hg38` reference is used for downstream validation alignment and
  Mutect2 calling.
- Clone mutation files are cumulative after `build_final_clone_mutations.py`.
- Use Slurm for FASTQ generation, merging and validation.
- The default FASTQ materialization mode uses links into `~/.rbbt/var/jobs`;
  do not clean that cache unless outputs were materialized with `copy`.
