#!/usr/bin/env nextflow
// ============================================================
// TTA_HBV-OBI_NGS  —  Illumina NGS pipeline for HBV analysis
// Nextflow DSL2
// ============================================================

nextflow.enable.dsl = 2

// ── Parameter defaults ────────────────────────────────────

params.samplesheet       = null          // required: path to CSV samplesheet
params.outdir            = "results"
params.hbv_ref           = "${projectDir}/resources/hbv_references/hbv_genotype_panel.fasta"
params.obi_db            = "${projectDir}/resources/knowledge_base/obi_mutations.tsv"
params.genotype_refs     = "${projectDir}/resources/hbv_references/hbv_genotype_panel.fasta"

// QC / trimming
params.skip_trimming     = false
params.fastp_args        = "--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50"

// Host filtering
params.host_filter       = false
params.host_ref          = null          // path to host reference FASTA (e.g. GRCh38)

// Variant calling
params.min_af            = 0.03          // minimum allele frequency for minor variants
params.min_depth         = 10            // minimum depth for consensus/variant calling
params.ivar_min_freq     = 0.03
params.ivar_min_depth    = 10

// Genotyping
params.genotype_method   = "phylo"       // "phylo" (MAFFT+IQ-TREE) or "blast"

// Reporting
params.multiqc_config    = "${projectDir}/config/multiqc_config.yml"
params.html_report       = true

// Execution
params.max_memory        = '16.GB'
params.max_cpus          = 8
params.max_time          = '48.h'

// ── Help message ──────────────────────────────────────────

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║          TTA HBV-OBI NGS Pipeline  (Illumina MVP)           ║
    ╚══════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf --samplesheet samplesheet.csv [options]

    Required:
        --samplesheet   Path to CSV samplesheet
                        (columns: sample,fastq_1,fastq_2)

    Optional input/reference:
        --hbv_ref       HBV reference panel FASTA
                        [${params.hbv_ref}]
        --obi_db        OBI mutation knowledge base TSV
                        [${params.obi_db}]
        --outdir        Results directory [${params.outdir}]

    QC / trimming:
        --skip_trimming       Skip fastp trimming [${params.skip_trimming}]
        --fastp_args          Arguments passed to fastp

    Host filtering:
        --host_filter         Enable host-read filtering [${params.host_filter}]
        --host_ref            Path to host reference FASTA

    Variant calling:
        --min_af              Min allele frequency [${params.min_af}]
        --min_depth           Min depth for calling [${params.min_depth}]

    Genotyping:
        --genotype_method     "phylo" or "blast" [${params.genotype_method}]

    Reporting:
        --html_report         Generate HTML report [${params.html_report}]

    Profiles:
        -profile conda        Use conda environments
        -profile docker       Use Docker containers
        -profile singularity  Use Singularity/Apptainer
        -profile slurm        Run on SLURM cluster (combine with above)
        -profile test         Run built-in smoke test

    Example:
        nextflow run main.nf \\
            --samplesheet data/samplesheet.csv \\
            --outdir results/run_001 \\
            -profile conda

    Documentation: docs/USAGE.md
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ── Validate core inputs ──────────────────────────────────

if (!params.samplesheet) {
    log.error "ERROR: --samplesheet is required. Run with --help for usage."
    exit 1
}

if (params.host_filter && !params.host_ref) {
    log.error "ERROR: --host_ref must be specified when --host_filter is enabled."
    exit 1
}

// ── Include modules ───────────────────────────────────────

include { VALIDATE_SAMPLESHEET } from './workflow/modules/validate_samplesheet'
include { FASTQC                } from './workflow/modules/fastqc'
include { FASTP                 } from './workflow/modules/fastp'
include { MULTIQC               } from './workflow/modules/multiqc'
include { BWA_MEM2_INDEX        } from './workflow/modules/bwa_mem2'
include { BWA_MEM2_MEM          } from './workflow/modules/bwa_mem2'
include { SAMTOOLS_SORT_INDEX   } from './workflow/modules/samtools'
include { SAMTOOLS_FLAGSTAT     } from './workflow/modules/samtools'
include { SAMTOOLS_DEPTH        } from './workflow/modules/samtools'
include { IVAR_CONSENSUS        } from './workflow/modules/ivar'
include { IVAR_VARIANTS         } from './workflow/modules/ivar'
include { MAFFT                 } from './workflow/modules/mafft'
include { IQTREE                } from './workflow/modules/iqtree'
include { ANNOTATE_VARIANTS     } from './workflow/modules/annotate'
include { GENOTYPE_ASSIGN       } from './workflow/modules/genotype'
include { REPORT                } from './workflow/modules/report'

// Optional host filtering
include { HOST_FILTER           } from './workflow/modules/host_filter'

// ── Log pipeline info ─────────────────────────────────────

log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║          TTA HBV-OBI NGS Pipeline  (Illumina MVP)           ║
    ╚══════════════════════════════════════════════════════════════╝
    samplesheet  : ${params.samplesheet}
    outdir       : ${params.outdir}
    hbv_ref      : ${params.hbv_ref}
    host_filter  : ${params.host_filter}
    genotype_method: ${params.genotype_method}
    min_af       : ${params.min_af}
    min_depth    : ${params.min_depth}
    """.stripIndent()

// ── Main workflow ─────────────────────────────────────────

workflow {

    // 1. Validate samplesheet and create channel
    VALIDATE_SAMPLESHEET(
        file(params.samplesheet, checkIfExists: true)
    )

    // Build reads channel: [meta, [read1, read2]] or [meta, [read]]
    ch_reads = VALIDATE_SAMPLESHEET.out.csv
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [id: row.sample, single_end: (row.fastq_2 == null || row.fastq_2 == '')]
            def reads = meta.single_end
                ? [ file(row.fastq_1, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true) ]
            [ meta, reads ]
        }

    // 2. FastQC on raw reads
    FASTQC(ch_reads)
    ch_fastqc_raw = FASTQC.out.zip

    // 3. Trimming with fastp
    if (!params.skip_trimming) {
        FASTP(ch_reads)
        ch_trimmed   = FASTP.out.reads
        ch_fastp_log = FASTP.out.json
    } else {
        ch_trimmed   = ch_reads
        ch_fastp_log = Channel.empty()
    }

    // 4. (Optional) Host read filtering
    if (params.host_filter) {
        HOST_FILTER(
            ch_trimmed,
            file(params.host_ref, checkIfExists: true)
        )
        ch_clean = HOST_FILTER.out.reads
    } else {
        ch_clean = ch_trimmed
    }

    // 5. Index HBV reference (once, shared)
    ch_hbv_ref = Channel.value(file(params.hbv_ref, checkIfExists: true))
    BWA_MEM2_INDEX(ch_hbv_ref)
    ch_index   = BWA_MEM2_INDEX.out.index

    // 6. Map to HBV reference
    BWA_MEM2_MEM(ch_clean, ch_index)

    // 7. Sort and index BAM
    SAMTOOLS_SORT_INDEX(BWA_MEM2_MEM.out.bam)

    // Collect mapping stats
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT_INDEX.out.bam_bai)
    SAMTOOLS_DEPTH(SAMTOOLS_SORT_INDEX.out.bam_bai)

    // 8. Variant calling with iVar
    IVAR_CONSENSUS(
        SAMTOOLS_SORT_INDEX.out.bam_bai,
        ch_hbv_ref
    )
    IVAR_VARIANTS(
        SAMTOOLS_SORT_INDEX.out.bam_bai,
        ch_hbv_ref
    )

    // 9. Genotype inference
    if (params.genotype_method == "phylo") {
        // Collect all consensus sequences
        ch_all_consensus = IVAR_CONSENSUS.out.fasta.collect { meta, fa -> fa }

        MAFFT(
            ch_all_consensus,
            file(params.genotype_refs, checkIfExists: true)
        )
        IQTREE(MAFFT.out.alignment)
        ch_genotype_result = GENOTYPE_ASSIGN(
            IVAR_CONSENSUS.out.fasta,
            IQTREE.out.tree
        )
    } else {
        // Blast/identity-based assignment
        ch_genotype_result = GENOTYPE_ASSIGN(
            IVAR_CONSENSUS.out.fasta,
            Channel.empty()
        )
    }

    // 10. Annotate OBI-associated variants
    ANNOTATE_VARIANTS(
        IVAR_VARIANTS.out.tsv,
        file(params.obi_db, checkIfExists: true)
    )

    // 11. MultiQC summary
    ch_multiqc_files = ch_fastqc_raw
        .mix(ch_fastp_log)
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat)
        .collect()

    MULTIQC(
        ch_multiqc_files,
        file(params.multiqc_config)
    )

    // 12. Per-sample and cohort report
    REPORT(
        ANNOTATE_VARIANTS.out.annotated.collect(),
        ch_genotype_result.collect(),
        SAMTOOLS_DEPTH.out.depth.collect(),
        MULTIQC.out.report
    )
}

// ── On completion ─────────────────────────────────────────

workflow.onComplete {
    log.info """
    Pipeline complete!
    Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration : ${workflow.duration}
    Results  : ${params.outdir}
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
