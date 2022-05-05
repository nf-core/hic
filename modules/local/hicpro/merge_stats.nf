process MERGE_STATS {
  label 'process_low'

  conda (params.enable_conda ? "conda-forge::python=3.7.6" : null)

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
