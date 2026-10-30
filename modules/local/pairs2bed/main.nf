process PAIRS2BED {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(pairs), path(index)

    output:
    tuple val(meta), path("*.reads.bed.gz"), emit: bed
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // One BED6 record per read end (5' position and strand) of every valid pair, the input MACS3 uses
    // to infer the 1D ChIP peaks of a HiChIP library (as in FitHiChIP's PeakInferHiChIP).
    // .pairs columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 ...
    """
    gzip -dc ${pairs} \\
        | awk 'BEGIN{OFS="\\t"} !/^#/ {
            if (\$3 > 0) print \$2, \$3-1, \$3, ".", "0", \$6
            if (\$5 > 0) print \$4, \$5-1, \$5, ".", "0", \$7
        }' \\
        | gzip -c > ${prefix}.reads.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -1 | cut -f1 -d, | sed -e 's/GNU Awk //')
    END_VERSIONS
    """
}
