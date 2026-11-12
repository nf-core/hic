process GET_VALID_INTERACTION_DNASE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9 bioconda::pysam=0.19.0 bioconda::bx-python=0.8.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' :
        'biocontainers/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.validPairs"), emit:valid_pairs
    tuple val(meta), path("*RSstat"), optional: true, emit:stats
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mapped_2hic_dnase.py \\
        -r ${bam} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
