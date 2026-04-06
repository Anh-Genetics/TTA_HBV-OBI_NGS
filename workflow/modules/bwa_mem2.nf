// ── bwa_mem2.nf ───────────────────────────────────────────
// Index HBV reference and map reads with bwa-mem2

process BWA_MEM2_INDEX {
    tag "index"
    label 'process_medium'
    publishDir "${params.outdir}/references", mode: 'copy', saveAs: { fn -> fn }

    conda 'bioconda::bwa-mem2=2.2.1'

    input:
    path fasta

    output:
    tuple path(fasta), path("*.{0123,amb,ann,bwt.2bit.64,pac}"), emit: index
    path "versions.yml",                                          emit: versions

    script:
    """
    bwa-mem2 index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.0123 ${fasta}.amb ${fasta}.ann ${fasta}.bwt.2bit.64 ${fasta}.pac
    touch versions.yml
    """
}


process BWA_MEM2_MEM {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/bam", mode: 'copy', saveAs: { fn -> fn.endsWith('.bam') ? fn : null }

    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19'

    input:
    tuple val(meta), path(reads)
    tuple path(fasta), path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "versions.yml",                     emit: versions

    script:
    def prefix    = "${meta.id}"
    def rg        = "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"
    def read_args = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${fasta} \\
        ${read_args} \\
        | samtools view -bS -@ ${task.cpus} - \\
        > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.bam
    touch versions.yml
    """
}
