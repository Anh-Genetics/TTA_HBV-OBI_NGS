// ── report.nf ─────────────────────────────────────────────
// Generate per-sample and cohort reports

process REPORT {
    tag "cohort_report"
    label 'process_single'
    publishDir "${params.outdir}/report", mode: 'copy'

    conda 'bioconda::python=3.11 conda-forge::pandas=2.1 conda-forge::jinja2=3.1'

    input:
    path annotated_tsvs     // collected per-sample annotated variant TSVs
    path genotype_tsvs      // collected per-sample genotype TSVs
    path depth_tsvs         // collected per-sample depth TSVs
    path multiqc_report     // MultiQC HTML

    output:
    path "cohort_summary.tsv",          emit: cohort_summary
    path "per_sample_report/*.tsv",     emit: per_sample, optional: true
    path "hbv_ngs_report.html",         emit: html, optional: true
    path "pipeline_versions.yml",       emit: versions
    path "versions.yml",                emit: nf_versions

    script:
    def html_flag = params.html_report ? "--html" : ""
    """
    mkdir -p per_sample_report

    generate_report.py \\
        --annotated ${annotated_tsvs} \\
        --genotypes ${genotype_tsvs} \\
        --depths    ${depth_tsvs} \\
        --multiqc   ${multiqc_report} \\
        --outdir    . \\
        ${html_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p per_sample_report
    printf 'sample\\tgenotype\\ttotal_variants\\tobi_flags\\n' > cohort_summary.tsv
    touch hbv_ngs_report.html
    touch pipeline_versions.yml
    touch versions.yml
    """
}
