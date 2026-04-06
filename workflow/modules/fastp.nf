// ── fastp.nf ──────────────────────────────────────────────
// Adapter trimming and quality filtering with fastp

process FASTP {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastp/${meta.id}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.json') || fn.endsWith('.html') ? fn : null }

    conda 'bioconda::fastp=0.23.4'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads
    path  "${meta.id}.fastp.json",                 emit: json
    path  "${meta.id}.fastp.html",                 emit: html
    path  "versions.yml",                          emit: versions

    script:
    def prefix = "${meta.id}"
    def pe_args = meta.single_end ? "" : "--in2 ${reads[1]} --out2 ${prefix}_R2.fastq.gz"
    def in1 = meta.single_end ? reads : reads[0]
    """
    fastp \\
        --in1 ${in1} \\
        --out1 ${prefix}_R1.fastq.gz \\
        ${pe_args} \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread ${task.cpus} \\
        ${params.fastp_args} \\
        2>&1 | tee ${prefix}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | grep -oP '(?<=fastp )[0-9.]+')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    echo "" | gzip > ${prefix}_R1.fastq.gz
    touch ${prefix}.fastp.json
    touch ${prefix}.fastp.html
    touch versions.yml
    """
}
