// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process tads_insulation {
  tag "$sample - $res"
  label 'process_medium'
  publishDir "${params.outdir}/tads/insulation", mode: 'copy'

  when:
  !params.skip_tads && params.tads_caller =~ 'insulation'

  input:
  tuple val(sample), val(res), path(cool), val(r) 

  output:
  path("*tsv"), emit:insulation_tads

  script:
  """
  cooltools diamond-insulation --window-pixels ${cool} 15 25 50 > ${sample}_insulation.tsv
  """
}
