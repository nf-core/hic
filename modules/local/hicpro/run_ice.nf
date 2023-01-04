process ICE_NORMALIZATION{
    tag "$meta.id"
    label 'process_high_memory'

    conda "conda-forge::python=3.9 bioconda::iced=0.5.10 conda-forge::numpy=1.22.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' :
        'quay.io/biocontainers/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0'}"

    input:
    tuple val(meta), val(res), path(rmaps), path(bed)

    output:
    tuple val(meta), val(res), path("*iced.matrix"), path(bed), emit:maps
    path ("*.biases"), emit:bias
    path("versions.yml"), emit: versions

    script:
    prefix = rmaps.toString() - ~/(\.matrix)?$/
    """
    ice --filter_low_counts_perc ${params.ice_filter_low_count_perc} \
        --results_filename ${prefix}_iced.matrix \
        --filter_high_counts_perc ${params.ice_filter_high_count_perc} \
        --max_iter ${params.ice_max_iter} --eps ${params.ice_eps} --remove-all-zeros-loci --output-bias 1 --verbose 1 ${rmaps}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
        iced: \$(python -c "import iced; print(iced.__version__)")
    END_VERSIONS
    """
}
