# Standard Operating Procedure (SOP)
## HBV OBI NGS Analysis — TTA_HBV-OBI_NGS Pipeline

**Version:** 1.0.0  
**Date:** 2024-01  
**Repository:** https://github.com/Anh-Genetics/TTA_HBV-OBI_NGS

---

## 1. Purpose and Scope

This SOP describes the procedure for analysing Illumina paired-end (or single-end) FASTQ data from HBV-amplified samples, with the goals of:

1. Determining HBV genotype and subgenotype from pre-S1/pre-S2/S region amplicons.
2. Identifying OBI-associated variants in the HBsAg coding region.
3. Generating quality-controlled, reproducible, documented results.

**Target regions:** pre-S1, pre-S2, S (surface antigen ORF), encompassing MHR and the "a" determinant.

---

## 2. Applicability

- Research applications only (not for clinical reporting without further validation).
- Applicable to samples where OBI is suspected (HBsAg−, anti-HBc±, HBV DNA low/detectable).
- Requires Illumina paired-end sequencing of HBV amplicons (≥ 150 bp reads recommended).

---

## 3. Pre-analytical Requirements

### 3.1 Sample requirements
- DNA extracted from serum/plasma with sufficient quality (A260/A280 ≥ 1.7).
- HBV DNA ≥ 10–100 IU/mL preferred; lower loads may yield incomplete coverage.
- PCR amplification of pre-S1/pre-S2/S using validated primers.

### 3.2 Sequencing requirements
- Illumina paired-end (2 × 150 bp or 2 × 250 bp).
- Minimum 10,000 reads per sample recommended.
- For minor variant detection: ≥ 10,000× average depth recommended.

---

## 4. Pipeline Steps

### Step 1: Prepare samplesheet
```
Column  | Content
--------|------------------------------------------------------
sample  | Unique sample ID (no spaces; letters/numbers/_/-/.)
fastq_1 | Absolute path to R1 FASTQ file (.fastq.gz)
fastq_2 | Absolute path to R2 FASTQ file (empty if single-end)
```

Template: `data/samplesheet_template.csv`

### Step 2: Prepare reference
Download HBV genotype panel as described in `docs/INSTALL.md`, Section 6.

### Step 3: Run pipeline
```bash
conda activate hbv_ngs
nextflow run main.nf \
    --samplesheet /path/to/samplesheet.csv \
    --hbv_ref resources/hbv_references/hbv_genotype_panel.fasta \
    --outdir results/run_YYYYMMDD \
    -profile conda \
    -resume
```

### Step 4: Review QC
- Open `results/multiqc/multiqc_report.html`.
- Check per-sample read counts, Q30 rates, adapter content.
- **Minimum acceptable:** ≥ 1000 reads passing QC, ≥ 80% Q30.

### Step 5: Review mapping
- Check `results/qc/flagstat/*.flagstat`.
- **Minimum acceptable:** ≥ 50% reads mapped to HBV reference.
- Low mapping may indicate failed amplification or sample contamination.

### Step 6: Review coverage
- Check `results/qc/depth/*.depth.tsv`.
- **Minimum for consensus:** 10× average depth.
- **Recommended for minor variants:** 100× average depth.

### Step 7: Review genotype assignment
- Open `results/genotyping/<sample>/<sample>.genotype.tsv`.
- Confidence levels: `high` (≥ 90% identity or branch dist < 0.05), `medium`, `low`.
- Low-confidence results require manual review of phylogenetic tree.

### Step 8: Review OBI-associated variants
- Open `results/annotation/<sample>/<sample>.obi_flags.tsv`.
- Interpret flags using the evidence categories:
  - **Category A:** Strong evidence; commonly reported in OBI/immune-escape literature.
  - **Category B:** Moderate evidence; reported in OBI cohorts, mechanism plausible.
  - **Category C:** Preliminary; single report or limited functional evidence.

### Step 9: Cohort summary
- `results/report/cohort_summary.tsv` provides a per-sample overview.
- HTML report at `results/report/hbv_ngs_report.html`.

---

## 5. Quality Control Acceptance Criteria

| QC metric | Minimum | Recommended |
|-----------|---------|-------------|
| Reads passing QC | ≥ 1,000 | ≥ 10,000 |
| Q30 base rate | ≥ 70% | ≥ 85% |
| HBV mapping rate | ≥ 50% | ≥ 80% |
| Mean coverage depth | ≥ 10× | ≥ 100× |
| Genome fraction covered (>0×) | ≥ 50% | ≥ 90% |

Samples failing minimum criteria should be **flagged as QC-fail** and not used for genotyping or mutation analysis without additional validation.

---

## 6. Interpretation Guidelines

### Genotype assignment
- Report at genotype level (A–I) with confidence.
- Report subgenotype only if confidence is `high` and reference coverage is adequate.
- State: "Based on pre-S/S region (partial genome); full-genome phylogeny recommended for definitive subgenotyping."

### OBI mutation interpretation
- Report variants by category (A/B/C) and mechanism.
- Never state "patient has OBI" based solely on sequencing; OBI is a clinical diagnosis.
- Suggested language:
  - *"This sample carries [variant X] in the [MHR/a-determinant] region, a Category A variant associated with diagnostic escape and immune escape in OBI cases (PMID:XXXXXXX)."*

### Minor variants (ALT_FREQ 0.01–0.10)
- Report with explicit caveat: *"Minor variant detected at [X]% allele frequency; interpretation should be cautious as this may approach sequencing error rates at lower frequencies."*
- Confirm by independent sequencing if clinically significant.

---

## 7. Provenance and Reproducibility

All pipeline runs generate:
- `results/pipeline_info/execution_report_*.html` — Nextflow execution report
- `results/pipeline_info/execution_trace_*.tsv` — per-task trace
- Software versions captured in each module's `versions.yml`
- Full parameters logged in Nextflow log

**Best practices:**
- Record Nextflow command with all parameters in lab notebook.
- Archive `results/pipeline_info/` with results.
- Record pipeline version (git tag/commit).

---

## 8. Known Limitations

1. Pre-S/S partial sequences may misclassify recombinant strains.
2. iVar consensus masks positions below `--min_depth` as `N`.
3. OBI knowledge base (v1.0) is curated from literature up to 2024; new variants may not be included.
4. Host filtering requires a host reference genome (not included; see INSTALL.md).
5. Phylogenetic method (IQ-TREE) may be slow for cohorts > 100 samples.

---

## 9. Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2024-01 | Initial MVP release |
