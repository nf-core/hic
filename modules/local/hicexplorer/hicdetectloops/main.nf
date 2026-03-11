process HICDETECTLOOPS {
    label 'process_medium'
    tag "$meta.id"

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path('*.bedGraph') , emit: bedgraph
    path("versions.yml")                , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.hicdetectloops"
    """
    hicDetectLoops \\
        ${args} \\
        --threads ${task.cpus} \\
        --matrix ${matrix} \\
        --outFileName ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicDetectLoops --version 2>&1 | sed 's/hicDetectLoops //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.hicdetectloops"
    """
    touch  ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicDetectLoops --version 2>&1 | sed 's/hicDetectLoops //')
    END_VERSIONS
    """
}
