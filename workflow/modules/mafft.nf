// ── mafft.nf ──────────────────────────────────────────────
// Multiple sequence alignment with MAFFT for genotyping

process MAFFT {
    tag "mafft_alignment"
    label 'process_medium'
    publishDir "${params.outdir}/genotyping/alignment", mode: 'copy'

    conda 'bioconda::mafft=7.525'

    input:
    path consensus_files      // collected list of consensus FASTA
    path genotype_refs        // reference sequences for genotypes A-I

    output:
    path "combined_aligned.fasta", emit: alignment
    path "versions.yml",           emit: versions

    script:
    """
    # Concatenate consensus sequences with genotype references
    cat ${genotype_refs} ${consensus_files} > combined_unaligned.fasta

    mafft \\
        --auto \\
        --thread ${task.cpus} \\
        combined_unaligned.fasta \\
        > combined_aligned.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | head -1 | sed 's/MAFFT v//')
    END_VERSIONS
    """

    stub:
    """
    cat ${genotype_refs} > combined_aligned.fasta
    touch versions.yml
    """
}
