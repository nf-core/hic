process MERGE_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0' :
        'quay.io/biocontainers/mulled-v2-c6ff206325681cbb9c9ef890bb8de554172c0483:713df51cd897ceb893b9a6e6420f527d83c2ed95-0'}"

    input:
    tuple val(meta), path(fstat)

    output:
    path("${meta.id}/"), emit: mqc
    path("*.{mmapstat,mpairstat,mRSstat}"), emit: stats
    path("versions.yml"), emit:versions

    script:
    if ( (fstat =~ /.mapstat/) ){ ext = "${meta.mates}.mmapstat" }
    if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
    if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }
    """
    mkdir -p ${meta.id}
    merge_statfiles.py -f ${fstat} > ${meta.id}.${ext}
    cp *${ext} ${meta.id}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
