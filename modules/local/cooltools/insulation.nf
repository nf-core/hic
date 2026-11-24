/*
 * Cooltools - diamond-insulation
 */

process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_medium'

    // Pin pandas=1.5.3 (last 1.x release, requires Python >=3.8) so conda does not resolve pandas>=2.0 which breaks cooltools 0.5.1.
    // The py37 biocontainer bundles an older pandas 1.x, so minor float/formatting differences in insulation.tsv are ignored in tests/.nftignore.
    conda "bioconda::cooltools=0.5.1 bioconda::cooler=0.8.11 conda-forge::pandas=1.5.3"
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
