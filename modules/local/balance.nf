process COOLER_BALANCE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), val(resolution), path(cool)

    output:
    tuple val(meta), val(resolution), path("*_norm.cool"), emit: cool
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cp ${cool} ${cool.baseName}_norm.cool
    cooler balance \\
        $args \\
        -p ${task.cpus} \\
        ${cool.baseName}_norm.cool

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
