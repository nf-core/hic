process HICPLOTMATRIX {
    label 'process_medium'
    tag "$meta.id"

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(matrix), path(tads), path(loops), path(bigwig)

    output:
    tuple val(meta), path('*.png') , emit: png
    path("versions.yml")           , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.hicplotmatrix"
    def tads_arg = tads ? "--tads ${tads}" : ''
    def loops_arg = loops ? "--loops ${loops}" : ''
    def bigwig_arg = bigwig ? "--bigwig ${bigwig}" : ''
    """
    hicPlotMatrix \\
        ${args} \\
        --matrix ${matrix} \\
        ${tads_arg} \\
        ${loops_arg} \\
        ${bigwig_arg} \\
        --outFileName ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotMatrix --version 2>&1 | sed 's/hicPlotMatrix //')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.hicplotmatrix"
    """
    touch  ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotMatrix --version 2>&1 | sed 's/hicPlotMatrix //')
    END_VERSIONS
    """
}
