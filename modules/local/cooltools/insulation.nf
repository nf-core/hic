/*
 * Cooltools - diamond-insulation
 */

process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_medium'

    // Pin cooler=0.9.3 alongside cooltools=0.7.1 (see open2c/cooltools#522) so cooler>=0.10 does not break _load_attrs / annotate()
    conda "bioconda::cooltools=0.7.1 bioconda::cooler=0.9.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/293465344cd2adc3c57236246b36380a2514a6485bcef57d8d8c5128a2c32b15/data' :
        'community.wave.seqera.io/library/cooltools_cooler:76510a8f37480f17' }"

    input:
    tuple val(meta), path(cool)

    output:
    path("*tsv"), emit:tsv
    path("versions.yml"), emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools insulation ${cool} ${args} > ${prefix}_insulation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version | grep 'cooltools, version ' | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
