process COMBINE_MATES {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9  bioconda::pysam=0.19.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' :
        'biocontainers/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*bwt2pairs.bam"), emit:bam
    tuple val(meta), path("*.pairstat"), optional:true, emit:stats
    path("versions.yml"), emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    mergeSAM.py -f ${bam[0]} -r ${bam[1]} -o ${prefix}_bwt2pairs.bam ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
