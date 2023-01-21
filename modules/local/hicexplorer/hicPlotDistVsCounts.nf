/*
 * hicexplorer - Genomic distance/counts plots
 */

process HIC_PLOT_DIST_VS_COUNTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool)

    output:
    path("*distcount*"), emit:results
    path("versions.yml"), emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hicPlotDistVsCounts --matrices ${cool} \
                        --plotFile ${prefix}_distcount.png \
                        --outFileData ${prefix}_distcount.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotDistVsCounts --version 2>&1 | sed 's/hicPlotDistVsCounts //')
    END_VERSIONS
    """
}
