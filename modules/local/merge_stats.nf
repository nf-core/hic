// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process merge_stats {
   tag "$ext"
   label 'process_low'
   publishDir "${params.outdir}/hicpro/", mode: params.publish_dir_mode,
               saveAs: {filename -> if (filename.endsWith("stat")) "stats/$filename"}

   input:
   tuple val(prefix), path(fstat) 

   output:
   path("stats/"), emit:mqc_mstats
   path("*stat"), emit:all_mstats

  script:
  sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)/
  if ( (fstat =~ /.mapstat/) ){ ext = "mmapstat" }
  if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
  if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }
  """
  merge_statfiles.py -f ${fstat} > ${prefix}.${ext}
  mkdir -p stats/${sample}
  cp ${prefix}.${ext} stats/${sample}/
  """
}
