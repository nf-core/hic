/*
 * hicexplorer - hicFindTADs
 */

process HIC_FIND_TADS {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hicexplorer=3.7.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

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
