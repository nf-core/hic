process SPLIT_COOLER_DUMP {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bedpe)

    output:
    path "*.txt", emit: matrix
    path "versions.yml", emit: versions

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
