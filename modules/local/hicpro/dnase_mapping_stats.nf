process MAPPING_STATS_DNASE {
    tag "$sample = $bam"
    label 'process_medium'

    conda "bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"


    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), emit:bam
    tuple val(meta), path("${prefix}.mapstat"), emit:stats

    script:
    prefix = meta.id + "_" + meta.chunk + "_" + meta.mates
    tag = meta.mates
    """
    echo "## ${prefix}" > ${prefix}.mapstat
    echo -n "total_${tag}\t" >> ${prefix}.mapstat
    samtools view -c ${bam} >> ${prefix}.mapstat
    echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
    samtools view -c -F 4 ${bam} >> ${prefix}.mapstat
    echo -n "global_${tag}\t" >> ${prefix}.mapstat
    samtools view -c -F 4 ${bam} >> ${prefix}.mapstat
    echo -n "local_${tag}\t0"  >> ${prefix}.mapstat
    """
}
