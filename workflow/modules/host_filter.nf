// ── host_filter.nf ────────────────────────────────────────
// Optional: remove host reads by mapping to host reference

process HOST_FILTER {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/host_filter/${meta.id}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.flagstat') ? fn : null }

    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19'

    input:
    tuple val(meta), path(reads)
    path  host_ref

    output:
    tuple val(meta), path("${meta.id}_filtered*.fastq.gz"), emit: reads
    path  "${meta.id}.host_filter.flagstat",                emit: flagstat
    path  "versions.yml",                                   emit: versions

    script:
    def prefix    = "${meta.id}"
    def rg        = "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA"
    def read_args = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    def r2_out    = meta.single_end ? "" : "-2 ${prefix}_filtered_R2.fastq.gz"
    """
    # Build index if not present
    if [ ! -f "${host_ref}.bwt.2bit.64" ]; then
        bwa-mem2 index ${host_ref}
    fi

    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${host_ref} \\
        ${read_args} \\
        | samtools view -bS -@ ${task.cpus} - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.host.sorted.bam

    samtools index ${prefix}.host.sorted.bam

    # Extract unmapped reads (not mapped to host)
    samtools flagstat ${prefix}.host.sorted.bam > ${prefix}.host_filter.flagstat

    samtools fastq \\
        -f 4 \\
        -@ ${task.cpus} \\
        -1 ${prefix}_filtered_R1.fastq.gz \\
        ${r2_out} \\
        -0 /dev/null \\
        -s /dev/null \\
        ${prefix}.host.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    echo "" | gzip > ${prefix}_filtered_R1.fastq.gz
    touch ${prefix}.host_filter.flagstat
    touch versions.yml
    """
}
