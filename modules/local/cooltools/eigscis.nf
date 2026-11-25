/*
 * cooltools - call_compartments
 */

process COOLTOOLS_EIGSCIS {
    tag "${meta.id}"
    label 'process_medium'

    // Pin cooler=0.9.3 alongside cooltools=0.7.1 (see open2c/cooltools#522) so cooler>=0.10 does not break _load_attrs / annotate()
    conda "bioconda::cooltools=0.7.1 bioconda::cooler=0.9.3 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3ffafacb71ffd4570e7da10e4d168adbf54dca242959cdea257caa695d292511/data' :
        'community.wave.seqera.io/library/cooltools_cooler_ucsc-bedgraphtobigwig:940b31e720f1497d' }"

    input:
    tuple val(meta), path(cool), val(resolution)
    path(fasta)
    path(chrsize)

    output:
    path("*compartments*"), emit: results
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools genome binnify --all-names ${chrsize} ${resolution} > genome_bins.txt
    cooltools genome gc genome_bins.txt ${fasta} > genome_gc.txt
    cooltools eigs-cis ${args} --phasing-track genome_gc.txt -o ${prefix}_compartments ${cool}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version | grep 'cooltools, version ' | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
