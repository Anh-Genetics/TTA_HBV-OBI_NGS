// ── validate_samplesheet.nf ──────────────────────────────
// Validate samplesheet CSV and emit it as a channel

process VALIDATE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_single'

    conda 'bioconda::python=3.11 conda-forge::pandas=2.1'

    input:
    path samplesheet

    output:
    path "${samplesheet}", emit: csv
    path "versions.yml",   emit: versions

    script:
    """
    validate_samplesheet.py \\
        --samplesheet ${samplesheet} \\
        --schema ${projectDir}/schema/samplesheet_schema.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cp ${samplesheet} ${samplesheet}
    touch versions.yml
    """
}
