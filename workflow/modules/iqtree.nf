// ── iqtree.nf ─────────────────────────────────────────────
// Maximum-likelihood phylogenetic tree with IQ-TREE 2

process IQTREE {
    tag "iqtree"
    label 'process_long'
    publishDir "${params.outdir}/genotyping/phylogeny", mode: 'copy'

    conda 'bioconda::iqtree=2.2.6'

    input:
    path alignment

    output:
    path "*.treefile",    emit: tree
    path "*.iqtree",      emit: log
    path "versions.yml",  emit: versions

    script:
    """
    iqtree2 \\
        -s ${alignment} \\
        -m GTR+G \\
        -B 1000 \\
        -T ${task.cpus} \\
        --prefix iqtree_out \\
        -redo

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 --version 2>&1 | head -1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+')
    END_VERSIONS
    """

    stub:
    """
    touch iqtree_out.treefile
    touch iqtree_out.iqtree
    touch versions.yml
    """
}
