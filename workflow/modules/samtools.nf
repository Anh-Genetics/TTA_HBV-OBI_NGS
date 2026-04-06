// ── samtools.nf ───────────────────────────────────────────
// Sort/index BAM, flagstat, and depth

process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/bam", mode: 'copy'

    conda 'bioconda::samtools=1.19'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai
    path "versions.yml",                                                                emit: versions

    script:
    def prefix = "${meta.id}"
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.sorted.bam \\
        ${bam}

    samtools index ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch versions.yml
    """
}


process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/qc/flagstat", mode: 'copy'

    conda 'bioconda::samtools=1.19'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.flagstat"), emit: flagstat
    path "versions.yml",                          emit: versions

    script:
    def prefix = "${meta.id}"
    """
    samtools flagstat \\
        -@ ${task.cpus} \\
        ${bam} \\
        > ${prefix}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.flagstat
    touch versions.yml
    """
}


process SAMTOOLS_DEPTH {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/qc/depth", mode: 'copy'

    conda 'bioconda::samtools=1.19'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.depth.tsv"), emit: depth
    path "versions.yml",                           emit: versions

    script:
    def prefix = "${meta.id}"
    """
    samtools depth \\
        -a \\
        -d 0 \\
        ${bam} \\
        > ${prefix}.depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.depth.tsv
    touch versions.yml
    """
}
