/*
 * hicexplorer - hicFindTADs
 */

process HIC_FIND_TADS {
    label 'process_medium'
    tag "$meta.id"

    conda "bioconda::hicexplorer=4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:4.0.0--py312hfffc319_0' :
        'biocontainers/hicexplorer:4.0.0--py312hfffc319_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    path("*hicfindtads*"), emit:results
    path("versions.yml"), emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hicFindTADs --matrix ${cool} \
        --outPrefix ${prefix}_hicfindtads \
        ${args} \
        --numberOfProcessors ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicFindTADs --version 2>&1 | sed 's/hicFindTADs //')
    END_VERSIONS
    """
}
