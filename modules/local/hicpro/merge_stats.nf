process MERGE_STATS {
  label 'process_low'

  input:
  tuple val(meta), path(fstat) 

  output:
  path("stats/"), emit:mqc_mstats
  path("*stat"), emit:all_mstats

  script:
  if ( (fstat =~ /.mapstat/) ){ ext = "mmapstat" }
  if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
  if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }
  """
  merge_statfiles.py -f ${fstat} > ${meta.id}.${ext}
  mkdir -p stats/${meta.id}
  cp ${meta.id}.${ext} stats/${meta.id}/
  """
}
