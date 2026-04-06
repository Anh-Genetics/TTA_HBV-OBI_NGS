# Installation Guide — TTA_HBV-OBI_NGS Pipeline

This guide covers installation on **WSL2 Ubuntu 24.04.1** (recommended) and any Linux/macOS system.

---

## 1. Prerequisites

### System requirements
- Ubuntu 24.04 (WSL2 or native), macOS 12+, or any modern Linux
- Minimum: 8 GB RAM, 4 CPU cores, 50 GB disk
- Recommended: 16 GB RAM, 8 CPU cores, 200 GB disk (for host-filtering with GRCh38)

### WSL2 setup (Windows users)
```powershell
# In PowerShell (as Administrator)
wsl --install -d Ubuntu-24.04
wsl --set-default-version 2
```

---

## 2. Install Nextflow

Nextflow requires Java 11+:

```bash
# Install Java (if not present)
sudo apt update && sudo apt install -y default-jdk

# Verify Java
java -version

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

---

## 3. Install Conda / Mamba (recommended)

```bash
# Install Miniforge (includes mamba)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
source ~/.bashrc

# Verify
mamba --version
```

---

## 4. Clone the repository

```bash
git clone https://github.com/Anh-Genetics/TTA_HBV-OBI_NGS.git
cd TTA_HBV-OBI_NGS
```

---

## 5. Create the conda environment

```bash
mamba env create -f environment.yml
conda activate hbv_ngs
```

This installs: fastqc, multiqc, fastp, bwa-mem2, samtools, ivar, mafft, iqtree, python, biopython, pandas.

---

## 6. Set up HBV reference panel

A reference FASTA is **required** for mapping and genotyping.

```bash
# Download representative sequences from NCBI (requires biopython)
conda activate hbv_ngs
python - <<'EOF'
from Bio import Entrez, SeqIO
import sys

Entrez.email = "your@email.com"  # replace with your email

# Representative genotype A-I accessions
accessions = {
    "A2": "J02203",
    "B1": "AB073858",
    "C1": "AB014381",
    "D1": "X65257",
    "E":  "X75657",
    "F1": "X69798",
    "G":  "AF160501",
    "H":  "AY090454",
    "I":  "FJ023669",
}

with open("resources/hbv_references/hbv_genotype_panel.fasta", "w") as out:
    for subgt, acc in accessions.items():
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
            seq = handle.read()
            # Rename header to include genotype label
            seq = seq.replace(f">{acc}", f">HBV_genotype_{subgt}_{acc}", 1)
            out.write(seq)
            print(f"Downloaded {subgt} ({acc})")
        except Exception as e:
            print(f"ERROR downloading {acc}: {e}", file=sys.stderr)
EOF
```

---

## 7. Verify installation

```bash
# Run the smoke test
bash tests/smoke_test.sh
```

Expected output: all tests PASS (Nextflow stub-run skipped if Nextflow not installed).

---

## 8. Docker (alternative)

If you prefer Docker:

```bash
# Install Docker on Ubuntu
sudo apt install -y docker.io
sudo usermod -aG docker $USER
newgrp docker

# Run pipeline with Docker profile
nextflow run main.nf \
    --samplesheet data/samplesheet.csv \
    --outdir results \
    -profile docker
```

---

## 9. Troubleshooting

| Problem | Solution |
|---------|----------|
| `java.lang.OutOfMemoryError` | Increase Java heap: `export NXF_OPTS="-Xms1g -Xmx4g"` |
| `bwa-mem2: command not found` | Activate conda env: `conda activate hbv_ngs` |
| WSL2 memory limit | Add `%wsl2-memory=12GB` to `%UserProfile%\.wslconfig` |
| Slow I/O on WSL2 | Place data on `/home/username/` not `/mnt/c/` |
