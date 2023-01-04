process TRIM_READS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads)
    val(motif)

    output:
    tuple val(meta), path("*trimmed.fastq.gz"), emit: fastq
    path("versions.yml") , emit: versions

    script:
    """
    zcat ${reads} > tmp.fastq
    cutsite_trimming --fastq tmp.fastq \\
        --cutsite ${motif[0]} \\
        --out ${reads.simpleName}_trimmed.fastq
    gzip ${reads.simpleName}_trimmed.fastq
    /bin/rm -f tmp.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | head -1 | cut -d" " -f2)
    END_VERSIONS
    """
}
