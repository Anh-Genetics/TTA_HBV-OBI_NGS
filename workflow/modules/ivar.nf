// ── ivar.nf ───────────────────────────────────────────────
// Consensus sequence and variant calling with iVar

process IVAR_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/consensus", mode: 'copy'

    conda 'bioconda::ivar=1.4.2 bioconda::samtools=1.19'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("${meta.id}.consensus.fa"),    emit: fasta
    tuple val(meta), path("${meta.id}.consensus.qual"),  emit: qual
    path "versions.yml",                                 emit: versions

    script:
    def prefix   = "${meta.id}"
    def min_freq  = params.ivar_min_freq
    def min_depth = params.ivar_min_depth
    """
    samtools mpileup \\
        -aa \\
        -A \\
        -d 0 \\
        -Q 0 \\
        --reference ${fasta} \\
        ${bam} \\
        | ivar consensus \\
            -p ${prefix}.consensus \\
            -m ${min_depth} \\
            -t ${min_freq} \\
            -n N

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | head -1 | sed 's/iVar version //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    echo ">consensus" > ${prefix}.consensus.fa
    echo "ATCGATCGATCG" >> ${prefix}.consensus.fa
    touch ${prefix}.consensus.qual
    touch versions.yml
    """
}


process IVAR_VARIANTS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/variants", mode: 'copy'

    conda 'bioconda::ivar=1.4.2 bioconda::samtools=1.19'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("${meta.id}.ivar.tsv"), emit: tsv
    path "versions.yml",                          emit: versions

    script:
    def prefix    = "${meta.id}"
    def min_freq  = params.ivar_min_freq
    def min_depth = params.ivar_min_depth
    """
    samtools mpileup \\
        -aa \\
        -A \\
        -d 0 \\
        -B \\
        -Q 0 \\
        --reference ${fasta} \\
        ${bam} \\
        | ivar variants \\
            -p ${prefix}.ivar \\
            -m ${min_depth} \\
            -t ${min_freq} \\
            -r ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | head -1 | sed 's/iVar version //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    printf 'REGION\\tPOS\\tREF\\tALT\\tREF_DP\\tREF_RV\\tREF_QUAL\\tALT_DP\\tALT_RV\\tALT_QUAL\\tALT_FREQ\\tTOTAL_DP\\tPVAL\\tPASS\\n' \\
        > ${prefix}.ivar.tsv
    touch versions.yml
    """
}
