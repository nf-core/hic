/*
 * Cooltools - diamond-insulation
 */

process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::cooltools=0.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.5.1--py37h37892f8_0' :
        'biocontainers/cooltools:0.5.1--py37h37892f8_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    path("*tsv"), emit:tsv
    path("versions.yml"), emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools insulation ${cool} ${args} > ${prefix}_insulation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version | grep 'cooltools, version ' | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
