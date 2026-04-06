# Usage Guide — TTA_HBV-OBI_NGS Pipeline

## Quick Start

```bash
# Activate conda environment
conda activate hbv_ngs

# Run the pipeline
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --outdir results/run_001 \
    -profile conda
```

---

## Samplesheet Format

Create a CSV file with the following columns:

| Column   | Required | Description |
|----------|----------|-------------|
| sample   | ✅        | Unique sample name (alphanumeric, `_`, `-`, `.`) |
| fastq_1  | ✅        | Path to R1 FASTQ (can be gzipped: `.fastq.gz`) |
| fastq_2  | ❌        | Path to R2 FASTQ (leave empty for single-end) |

**Example:**
```csv
sample,fastq_1,fastq_2
OBI_PATIENT_001,/data/fastq/OBI_001_R1.fastq.gz,/data/fastq/OBI_001_R2.fastq.gz
OBI_PATIENT_002,/data/fastq/OBI_002_R1.fastq.gz,/data/fastq/OBI_002_R2.fastq.gz
OBI_PATIENT_003,/data/fastq/OBI_003_SE.fastq.gz,
```

A template is provided at `data/samplesheet_template.csv`.

---

## Pipeline Stages

### 1. Input validation
```bash
nextflow run main.nf --samplesheet samplesheet.csv ...
```
- Validates CSV format, sample names, FASTQ paths
- Fails immediately on errors

### 2. Read QC
- **FastQC** on raw reads
- **fastp** adapter trimming + quality filtering (Q20 by default)
- **MultiQC** cohort summary

To skip trimming:
```bash
nextflow run main.nf --samplesheet ss.csv --skip_trimming ...
```

### 3. Host filtering (optional)
Remove human (or other host) reads before HBV mapping:
```bash
nextflow run main.nf \
    --samplesheet ss.csv \
    --host_filter \
    --host_ref /data/GRCh38/hg38.fa \
    ...
```

### 4. HBV reference mapping
- **bwa-mem2** maps clean reads to HBV reference panel
- **samtools** sorts and indexes BAM files
- Flagstat and depth coverage statistics are computed

### 5. Variant calling
- **iVar** consensus sequence (≥ min_depth coverage, ≥ min_af frequency)
- **iVar** variant TSV with allele frequencies
- Key parameters:
  - `--min_af 0.03` (3% minimum allele frequency)
  - `--min_depth 10` (minimum read depth)

### 6. Genotype inference

**Phylogenetic method (default, `--genotype_method phylo`):**
1. MAFFT aligns all consensus sequences + reference genotype panel
2. IQ-TREE builds ML tree (GTR+G, 1000 ultrafast bootstraps)
3. Genotype assigned by nearest reference in tree

**Identity-based method (`--genotype_method blast`):**
- Pairwise alignment to reference panel
- Faster but less accurate for subgenotypes

```bash
# Use identity method
nextflow run main.nf --genotype_method blast ...
```

### 7. OBI variant annotation
- Each variant mapped to OBI knowledge base (`resources/knowledge_base/obi_mutations.tsv`)
- Confidence categories: A (strong), B (moderate), C (preliminary)
- Flags: MHR/a-determinant mutations, stop codons, pre-S deletions, glycosylation changes

### 8. Reporting
- Per-sample annotated variant tables
- Cohort summary TSV
- Optional HTML report (`--html_report true`)
- MultiQC HTML

---

## All Parameters

```
--samplesheet PATH      Required. Path to CSV samplesheet.
--outdir PATH           Results directory. Default: results
--hbv_ref PATH          HBV reference FASTA panel. Default: resources/hbv_references/hbv_genotype_panel.fasta
--obi_db PATH           OBI mutations knowledge base TSV.
--skip_trimming         Skip fastp trimming.
--fastp_args STRING     Additional fastp arguments.
--host_filter           Enable host-read removal.
--host_ref PATH         Host reference FASTA (required if --host_filter).
--min_af FLOAT          Min allele frequency for minor variants. Default: 0.03
--min_depth INT         Min depth for consensus calling. Default: 10
--genotype_method STR   'phylo' or 'blast'. Default: phylo
--html_report           Generate HTML report. Default: true
--max_memory MEM        Max memory per process. Default: 16.GB
--max_cpus INT          Max CPUs per process. Default: 8
--max_time TIME         Max time per process. Default: 48.h
```

---

## Profiles

| Profile | Description |
|---------|-------------|
| `local` | Run locally (default) |
| `conda` | Use conda environments |
| `docker` | Use Docker containers |
| `singularity` | Use Singularity/Apptainer |
| `slurm` | Run on SLURM cluster |
| `test` | Use built-in test data |

Combine profiles: `-profile conda,slurm`

---

## Output Structure

```
results/
├── fastqc/                  # Per-sample FastQC reports
├── fastp/                   # fastp trimming reports + logs
├── bam/                     # Sorted BAM files + indices
├── qc/
│   ├── flagstat/            # samtools flagstat per sample
│   └── depth/               # samtools depth per sample
├── consensus/               # iVar consensus FASTA per sample
├── variants/                # iVar variant TSV per sample
├── annotation/              # Annotated variants + OBI flags
├── genotyping/
│   ├── alignment/           # MAFFT alignment
│   ├── phylogeny/           # IQ-TREE tree files
│   └── <sample>/            # Per-sample genotype TSV
├── multiqc/                 # MultiQC HTML report
├── report/
│   ├── cohort_summary.tsv   # Cohort summary table
│   ├── per_sample_report/   # Per-sample tables
│   └── hbv_ngs_report.html  # HTML report
└── pipeline_info/           # Nextflow execution reports + traces
```

---

## Example Commands

### Standard paired-end run
```bash
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --hbv_ref resources/hbv_references/hbv_genotype_panel.fasta \
    --outdir results/run_$(date +%Y%m%d) \
    -profile conda \
    -resume
```

### With host filtering (human samples)
```bash
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --host_filter \
    --host_ref /data/references/GRCh38.fasta \
    --outdir results/host_filtered \
    -profile conda
```

### High-sensitivity minor variant detection
```bash
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --min_af 0.01 \
    --min_depth 100 \
    --outdir results/high_sensitivity \
    -profile conda
```

### Dry run (structure check)
```bash
nextflow run main.nf \
    --samplesheet tests/samplesheet_test.csv \
    --hbv_ref tests/data/hbv_ref_small.fasta \
    -profile test \
    -stub-run
```

---

## Limitations

> ⚠️ **Please read before interpreting results:**

1. **Minor variants (ALT_FREQ < 0.05):** These are close to typical sequencing error rates. Treat with caution. Confirm by independent amplicon sequencing or replicate library preparation.

2. **Genotype assignment from partial sequences:** Pre-S1/pre-S2/S regions (~1.2 kb) are generally sufficient for genotype-level assignment (A–I) but subgenotype assignment may be unreliable for borderline cases. Full-genome sequencing is recommended for definitive subgenotyping.

3. **OBI diagnosis:** The presence of OBI-associated mutations does NOT confirm OBI. OBI is a clinical/virological diagnosis (HBsAg-negative, HBV DNA detectable). These mutations are supportive mechanistic evidence only.

4. **Reference panel bias:** Genotyping accuracy depends on the quality and breadth of the reference panel. Ensure it includes locally circulating genotypes.

5. **Recombinant strains:** Recombinant HBV strains (common in some Asian populations) may yield ambiguous genotype assignments.
