/*
 * Filter chromsize
 * Filter out chromosomes smaller than a given threshold
 */

process FILTER_CHROMSIZE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(sizes_file)

    output:
    tuple val(meta), path("filtered.sizes")

    script:
    """
    awk -v MIN=${params.min_size} '\$2 >= MIN' ${sizes_file} > filtered.sizes
    """
}
