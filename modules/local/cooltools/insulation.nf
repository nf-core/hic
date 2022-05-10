/*
 * Cooltools - diamond-insulation
 */

process INSULATION {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cooltools=0.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.5.1--py37h37892f8_0' :
        'quay.io/biocontainers/cooltools:0.5.1--py37h37892f8_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    path("*tsv"), emit:results
    path("versions.yml"), emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools insulation ${cool} ${args} > ${prefix}_insulation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
