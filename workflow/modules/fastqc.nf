// ── fastqc.nf ─────────────────────────────────────────────
// Run FastQC on raw reads (paired or single-end)

process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc/${meta.id}", mode: 'copy', saveAs: { fn -> fn }

    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.zip"),  emit: zip
    tuple val(meta), path("*.html"), emit: html
    path  "versions.yml",            emit: versions

    script:
    def prefix = "${meta.id}"
    def read_files = reads.collect { it }.join(' ')
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${read_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version 2>&1 | sed 's/FastQC //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}_R1_fastqc.zip
    touch ${prefix}_R1_fastqc.html
    touch versions.yml
    """
}
