process HICDETECTLOOPS {
    label 'process_medium'
    tag "$meta.id"

    conda "conda-forge::python=3.9 bioconda::hicexplorer=3.7.2 conda-forge::numpy=1.24.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path('*.bedGraph') , emit: bedgraph
    path("versions.yml"), topic: versions, emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.hicdetectloops"
    """
    # cleanlab >= 2.0 removed LearningWithNoisyLabels, which hicexplorer 3.7.2 imports
    # transitively via hicexplorer.lib.tadClassifier even though hicDetectLoops never uses it.
    python <<'NFCORE_HICDETECTLOOPS'
    import sys
    try:
        import cleanlab.classification
        if not hasattr(cleanlab.classification, "LearningWithNoisyLabels"):
            cleanlab.classification.LearningWithNoisyLabels = getattr(cleanlab.classification, "CleanLearning", object)
    except ImportError:
        pass
    from hicexplorer.hicDetectLoops import main
    sys.argv = ["hicDetectLoops"] + "${args}".split() + [
        "--threads", "${task.cpus}",
        "--matrix", "${matrix}",
        "--outFileName", "${prefix}.bedGraph",
    ]
    main()
    NFCORE_HICDETECTLOOPS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(python -c "from hicexplorer._version import __version__; print(__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.hicdetectloops"
    """
    touch  ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(python -c "from hicexplorer._version import __version__; print(__version__)")
    END_VERSIONS
    """
}
