process MERGE_BOWTIE2{
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam1), path(bam2)

    output:
    tuple val(meta), path("${prefix}_bwt2merged.bam"), emit: bam
    tuple val(meta), path("${prefix}.mapstat"), emit: stats
    path("versions.yml"), emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    tag = meta.mates
    """
    samtools merge -@ ${task.cpus} \\
        -f ${prefix}_bwt2merged.bam \\
        ${bam1} ${bam2}

    samtools sort -@ ${task.cpus} -m 800M \\
        -n \\
        -o ${prefix}_bwt2merged.sorted.bam \\
        ${prefix}_bwt2merged.bam

    mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam

    echo "## ${prefix}" > ${prefix}.mapstat
    echo -n "total_${tag}\t" >> ${prefix}.mapstat
    samtools view -c ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
    echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
    samtools view -c -F 4 ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
    echo -n "global_${tag}\t" >> ${prefix}.mapstat
    samtools view -c -F 4 ${bam1} >> ${prefix}.mapstat
    echo -n "local_${tag}\t"  >> ${prefix}.mapstat
    samtools view -c -F 4 ${bam2} >> ${prefix}.mapstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
