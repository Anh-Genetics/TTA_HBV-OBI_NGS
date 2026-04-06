// ── annotate.nf ───────────────────────────────────────────
// Annotate iVar variant calls against OBI mutation knowledge base

process ANNOTATE_VARIANTS {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/annotation/${meta.id}", mode: 'copy'

    conda 'bioconda::python=3.11 conda-forge::pandas=2.1 conda-forge::biopython=1.83'

    input:
    tuple val(meta), path(variants_tsv)
    path  obi_db

    output:
    tuple val(meta), path("${meta.id}.annotated.tsv"),  emit: annotated
    tuple val(meta), path("${meta.id}.obi_flags.tsv"),  emit: obi_flags
    path "versions.yml",                                emit: versions

    script:
    def prefix = "${meta.id}"
    """
    annotate_variants.py \\
        --sample ${prefix} \\
        --variants ${variants_tsv} \\
        --obi_db ${obi_db} \\
        --out_annotated ${prefix}.annotated.tsv \\
        --out_flags ${prefix}.obi_flags.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    printf 'REGION\\tPOS\\tREF\\tALT\\tALT_FREQ\\taa_change\\tobi_category\\teffect\\tconfidence\\n' \\
        > ${prefix}.annotated.tsv
    touch ${prefix}.obi_flags.tsv
    touch versions.yml
    """
}
