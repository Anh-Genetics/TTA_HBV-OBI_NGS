// ── genotype.nf ───────────────────────────────────────────
// Assign HBV genotype from phylogenetic tree or BLAST identity

process GENOTYPE_ASSIGN {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/genotyping/${meta.id}", mode: 'copy'

    conda 'bioconda::python=3.11 conda-forge::biopython=1.83 conda-forge::pandas=2.1'

    input:
    tuple val(meta), path(consensus_fasta)
    path  tree   // optional: IQ-TREE .treefile (may be empty channel)

    output:
    tuple val(meta), path("${meta.id}.genotype.tsv"), emit: genotype
    path "versions.yml",                              emit: versions

    script:
    def prefix    = "${meta.id}"
    def tree_arg  = tree ? "--tree ${tree}" : ""
    def method    = params.genotype_method
    """
    genotype_assign.py \\
        --sample ${prefix} \\
        --consensus ${consensus_fasta} \\
        --refs ${params.genotype_refs} \\
        --method ${method} \\
        ${tree_arg} \\
        --out ${prefix}.genotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        biopython: \$(python -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    printf 'sample\\tgenotype\\tsubgenotype\\tmethod\\tconfidence\\tnotes\\n' \\
        > ${prefix}.genotype.tsv
    echo "${prefix}\\tA\\tA2\\t${params.genotype_method}\\thigh\\tstub" >> ${prefix}.genotype.tsv
    touch versions.yml
    """
}
