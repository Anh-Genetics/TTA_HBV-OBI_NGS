// ── multiqc.nf ────────────────────────────────────────────
// Aggregate QC reports with MultiQC

process MULTIQC {
    tag "cohort"
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    conda 'bioconda::multiqc=1.21'

    input:
    path  qc_files
    path  multiqc_config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data
    path "versions.yml",        emit: versions

    script:
    def config_arg = multiqc_config.name != 'NO_FILE' ? "--config ${multiqc_config}" : ""
    """
    multiqc \\
        --force \\
        ${config_arg} \\
        --filename multiqc_report.html \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    touch versions.yml
    """
}
