/*
 * cooltools - call_compartments
 */

process COOLTOOLS_EIGSCIS {
    tag "${meta.id}"
    label 'process_medium'

    // Pin pandas=1.5.3 (last 1.x release, requires Python >=3.8) so conda does not resolve pandas>=2.0 which breaks cooltools 0.5.1.
    // The py37 biocontainer bundles an older pandas 1.x, kept as-is until cooltools is upgraded.
    conda "bioconda::cooltools=0.5.1 bioconda::cooler=0.8.11 conda-forge::pandas=1.5.3 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' :
        'biocontainers/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' }"

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
