process SPLIT_COOLER_DUMP {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bedpe)

    output:
    tuple val(meta), path("*.txt"), emit: matrix
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = bedpe.toString() - ~/(\_balanced)?.bedpe$/
    """
    cat ${bedpe} | awk '{OFS="\t"; print \$1,\$2,\$3}' > ${prefix}_raw.txt
    cat ${bedpe} | awk '{OFS="\t"; print \$1,\$2,\$4}' > ${prefix}_balanced.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(awk --version | head -1 | cut -f1 -d, | sed -e 's/GNU Awk //')
    END_VERSIONS
    """
}
