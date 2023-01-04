process HICPRO2PAIRS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pairix=0.3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairix:0.3.7--py36h30a8e3e_3' :
        'quay.io/biocontainers/pairix:0.3.7--py36h30a8e3e_3' }"

    input:
    tuple val(meta), path(vpairs)
    tuple val(meta2), path(chrsize)

    output:
    tuple val(meta), path("*.pairs.gz"), path("*.pairs.gz.px2"), emit: pairs
    path("versions.yml"), emit: versions

    script:
    prefix = "${meta.id}"
    """
    ##columns: readID chr1 pos1 chr2 pos2 strand1 strand2
    awk '{OFS="\t";print \$1,\$2,\$3,\$5,\$6,\$4,\$7}' $vpairs | bgzip -c > ${prefix}_contacts.pairs.gz
    ##sort -k2,2 -k4,4 -k3,3n -k5,5n ${prefix}_contacts.pairs | bgzip -c > ${prefix}_contacts.pairs.gz
    pairix -f ${prefix}_contacts.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairix: \$(echo \$(pairix 2>&1 | grep Version | sed -e 's/Version: //'))
    END_VERSIONS
    """
}
