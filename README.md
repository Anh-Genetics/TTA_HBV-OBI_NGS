# TTA_HBV-OBI_NGS

> **Production-ready Illumina NGS pipeline for HBV genotyping and OBI-associated mutation analysis**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)

---

## Overview

**TTA_HBV-OBI_NGS** is a reproducible, production-ready Nextflow DSL2 pipeline for:

- Analysing **Illumina paired-end FASTQ** data from HBV amplicons targeting the **pre-S1 / pre-S2 / S** region.
- Determining **HBV genotype** (A–I) with confidence scoring.
- Identifying **OBI-associated variants** — immune/diagnostic escape mutations, stop codons, pre-S deletions, glycosylation changes — using a curated knowledge base.
- Generating comprehensive **QC reports** and **per-sample/cohort summaries**.

**Pipeline stages:**

```
Input FASTQ → FastQC → fastp trimming → [Host filter] → bwa-mem2 mapping
    → iVar consensus + variant calling → Genotype (MAFFT + IQ-TREE)
    → OBI annotation → MultiQC + HTML report
```

---

## Quick Start

### 1. Install prerequisites

```bash
# Java (required for Nextflow)
sudo apt update && sudo apt install -y default-jdk

# Nextflow
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/

# Mamba (recommended over conda for speed)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash && source ~/.bashrc

# Clone and set up environment
git clone https://github.com/Anh-Genetics/TTA_HBV-OBI_NGS.git
cd TTA_HBV-OBI_NGS
mamba env create -f environment.yml
conda activate hbv_ngs
```

### 2. Prepare HBV reference (one-time)

```bash
# Download representative genotype A-I sequences from NCBI
python - <<'EOF'
from Bio import Entrez, SeqIO
Entrez.email = "your@email.com"
accessions = {"A2":"J02203","B1":"AB073858","C1":"AB014381","D1":"X65257",
               "E":"X75657","F1":"X69798","G":"AF160501","H":"AY090454","I":"FJ023669"}
with open("resources/hbv_references/hbv_genotype_panel.fasta","w") as out:
    for subgt, acc in accessions.items():
        h = Entrez.efetch(db="nucleotide",id=acc,rettype="fasta",retmode="text")
        out.write(h.read().replace(f">{acc}",f">HBV_genotype_{subgt}_{acc}",1))
        print(f"  {subgt} ({acc})")
EOF
```

### 3. Run the pipeline

```bash
# Create your samplesheet (see data/samplesheet_template.csv)
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --outdir results/run_001 \
    -profile conda
```

---

## Repository Structure

```
TTA_HBV-OBI_NGS/
├── main.nf                          # Pipeline entry point (Nextflow DSL2)
├── nextflow.config                  # Configuration + profiles
├── environment.yml                  # Conda environment
├── workflow/
│   ├── modules/                     # Nextflow process modules
│   │   ├── fastqc.nf
│   │   ├── fastp.nf
│   │   ├── multiqc.nf
│   │   ├── bwa_mem2.nf
│   │   ├── samtools.nf
│   │   ├── ivar.nf
│   │   ├── mafft.nf
│   │   ├── iqtree.nf
│   │   ├── annotate.nf
│   │   ├── genotype.nf
│   │   ├── report.nf
│   │   ├── host_filter.nf
│   │   └── validate_samplesheet.nf
│   └── bin/                         # Python helper scripts
│       ├── validate_samplesheet.py
│       ├── annotate_variants.py
│       ├── genotype_assign.py
│       └── generate_report.py
├── config/
│   ├── base.config                  # Resource allocation
│   ├── docker.config                # Container assignments
│   ├── slurm.config                 # SLURM settings
│   └── multiqc_config.yml           # MultiQC configuration
├── resources/
│   ├── hbv_references/              # HBV reference FASTA (download separately)
│   └── knowledge_base/
│       └── obi_mutations.tsv        # OBI mutation database (versioned)
├── schema/
│   └── samplesheet_schema.json      # JSON schema for samplesheet validation
├── data/
│   └── samplesheet_template.csv     # Samplesheet template
├── docs/
│   ├── INSTALL.md                   # Detailed installation guide
│   ├── USAGE.md                     # Full usage documentation
│   └── SOP.md                       # Standard Operating Procedure
└── tests/
    ├── smoke_test.sh                # Smoke / structural test
    ├── samplesheet_test.csv         # Test samplesheet
    └── data/                        # Test FASTQ + reference
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | required | CSV samplesheet path |
| `--hbv_ref` | `resources/hbv_references/hbv_genotype_panel.fasta` | HBV reference FASTA |
| `--outdir` | `results` | Results directory |
| `--min_af` | `0.03` | Min allele frequency for variant calling |
| `--min_depth` | `10` | Min depth for consensus |
| `--genotype_method` | `phylo` | `phylo` (MAFFT+IQ-TREE) or `blast` |
| `--host_filter` | `false` | Enable host read removal |
| `--host_ref` | `null` | Host reference FASTA |
| `--skip_trimming` | `false` | Skip fastp trimming |
| `--html_report` | `true` | Generate HTML report |

---

## Profiles

```bash
-profile conda          # Use conda environments (recommended)
-profile docker         # Use Docker containers
-profile singularity    # Use Singularity/Apptainer
-profile slurm          # SLURM cluster (combine with above)
-profile test           # Smoke test with bundled data
```

---

## OBI Mutation Knowledge Base

The file `resources/knowledge_base/obi_mutations.tsv` contains a curated, versioned database of HBV surface antigen variants associated with OBI, immune escape, and diagnostic escape. Each entry includes:

- Genomic position and amino acid change
- Region (MHR, a-determinant, pre-S1, pre-S2)
- Evidence category (A = strong, B = moderate, C = preliminary)
- Effect (diagnostic_escape, immune_escape, secretion_defect, stop_codon, etc.)
- Mechanism and literature references (PMIDs)

---

## Limitations

> ⚠️ **Important clinical and analytical limitations:**

1. **Research use only** — not validated for clinical diagnostic reporting.
2. **Minor variants (< 5% AF)** — may approach sequencing error rates; confirm by independent sequencing.
3. **Partial genome genotyping** — pre-S/S region is generally reliable for genotype A–I, but subgenotype assignment may be unreliable for recombinant strains.
4. **OBI diagnosis** — sequencing evidence alone cannot diagnose OBI; clinical and serological correlation is essential.
5. **Reference panel dependency** — genotyping accuracy depends on the breadth of the reference panel.

---

## Documentation

- [Installation guide](docs/INSTALL.md) — WSL2 Ubuntu setup, Nextflow, conda, references
- [Usage guide](docs/USAGE.md) — All parameters, examples, output structure
- [Standard Operating Procedure](docs/SOP.md) — QC criteria, interpretation guidelines

---

## Testing

```bash
# Run the smoke test (Python + structural validation)
bash tests/smoke_test.sh

# Nextflow stub-run (no tools required)
nextflow run main.nf \
    --samplesheet tests/samplesheet_test.csv \
    --hbv_ref tests/data/hbv_ref_small.fasta \
    -profile test \
    -stub-run
```

---

## Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | ≥ 23.04 | Workflow manager |
| FastQC | 0.12.1 | Read quality control |
| fastp | 0.23.4 | Adapter trimming |
| MultiQC | 1.21 | QC aggregation |
| bwa-mem2 | 2.2.1 | Short-read alignment |
| samtools | 1.19 | BAM processing |
| iVar | 1.4.2 | Consensus + variant calling |
| MAFFT | 7.525 | Multiple sequence alignment |
| IQ-TREE | 2.2.6 | Phylogenetic tree |
| Python | 3.11 | Helper scripts |
| Biopython | 1.83 | Sequence parsing |
| pandas | 2.1 | Data processing |

---

## Citation

If you use this pipeline, please cite the tools it depends on:

- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen S. et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
- **bwa-mem2**: Vasimuddin M. et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. *IPDPS*.
- **iVar**: Grubaugh N.D. et al. (2019). An amplicon-based sequencing framework for accurately measuring intrahost virus diversity. *PLoS Pathog*, 15(1).
- **MAFFT**: Katoh K. & Standley D.M. (2013). MAFFT Multiple Sequence Alignment Software Version 7. *Mol Biol Evol*, 30(4), 772–780.
- **IQ-TREE 2**: Minh B.Q. et al. (2020). IQ-TREE 2: New Models and Methods for Phylogenetic Inference. *Mol Biol Evol*, 37(5), 1530–1534.
- **MultiQC**: Ewels P. et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048.
- **Nextflow**: Di Tommaso P. et al. (2017). Nextflow enables reproducible computational workflows. *Nat Biotechnol*, 35, 316–319.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Author

**Anh-Genetics** — https://github.com/Anh-Genetics
