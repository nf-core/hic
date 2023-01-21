process MERGE_VALID_INTERACTION {
    tag "$prefix"
    label 'process_high_memory'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(vpairs)

    output:
    tuple val(meta), path("*.allValidPairs"), emit: valid_pairs
    path("${meta.id}/"), emit:mqc
    path("*mergestat"), emit:stats
    path("versions.yml"), emit: versions

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    """
    hicpro_merge_validpairs.sh ${args} -p ${prefix} ${vpairs}

    ## For MultiQC
    mkdir -p ${prefix}
    cp ${prefix}_allValidPairs.mergestat ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(echo \$(sort --version 2>&1 | head -1 | awk '{print \$NF}' 2>&1))
    END_VERSIONS
    """
}
