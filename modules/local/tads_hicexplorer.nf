// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process tads_hicexplorer {
  tag "$sample - $res"
  label 'process_medium'
  publishDir "${params.outdir}/tads/hicexplorer", mode: 'copy'

  when:
  !params.skip_tads && params.tads_caller =~ 'hicexplorer'

  input:
  tuple val(sample), val(res), path(cool), val(r) 

  output:
  path("*.{bed,bedgraph,gff}"), emit:hicexplorer_tads

  script:
  """
  hicFindTADs --matrix ${cool} \
  	      --outPrefix tad \
	      --correctForMultipleTesting fdr \
	      --numberOfProcessors ${task.cpus}
  """
}
