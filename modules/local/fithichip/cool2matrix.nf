process FITHICHIP_COOL2MATRIX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::cooler=0.10.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.3--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*_abs.bed"), path("*.matrix"), emit: matrix
    path("versions.yml"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // HiC-Pro style matrix for FitHiChIP: bins (chrom start end id) and raw pixels (id1 id2 count).
    // cooler's pixel ids (bin1_id, bin2_id) are 0-based row numbers of its bin table, so the bin ids written
    // here are NR-1 to match them. FitHiChIP (src/InteractionHicPro.r) joins the two files by id value, not by
    // row position, so only the agreement between the two files matters (HiC-Pro's own files are 1-based on
    // both sides, which works for the same reason).
    """
    cooler dump -t bins ${cool} \\
        | awk 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, NR-1}' \\
        > ${prefix}_abs.bed

    cooler dump -t pixels ${cool} > ${prefix}.matrix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
