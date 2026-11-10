process TRIM_READS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.10, conda-forge::pyahocorasick, conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pyahocorasick_sed:0125bb306ee03ca4'  :
        'community.wave.seqera.io/library/pyahocorasick_sed:dc5c90f342b65a2a' }"

    input:
    tuple val(meta), path(reads)
    val(motif)

    output:
    tuple val(meta), path("*trimmed.fastq.gz"), emit: fastq

    script:
    """
    cutsite_trimming.py --fastq ${reads} \\
        --cutsite ${motif[0]} \\
        --out ${reads.simpleName}_trimmed.fastq

    gzip ${reads.simpleName}_trimmed.fastq
    """
}
