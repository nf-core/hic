process PAIRS2BED {
    tag "$meta.id"
    label 'process_low'

    // same gawk environment and images as the nf-core gawk module
    conda "conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(pairs), path(index)

    output:
    tuple val(meta), path("*.reads.bed.gz"), emit: bed
    path("versions.yml"), topic: versions, emit: versions

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
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
