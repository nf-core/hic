process GET_RESTRICTION_FRAGMENTS {
    tag "$res_site"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::numpy=1.22.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' :
        'quay.io/biocontainers/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0'}"

    input:
    tuple val(meta), path(fasta)
    val(res_site)

    output:
    tuple val(meta), path("*.bed"), emit: results
    path("versions.yml"), emit: versions

    script:
    """
    digest_genome.py -r ${res_site} -o restriction_fragments.bed ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
