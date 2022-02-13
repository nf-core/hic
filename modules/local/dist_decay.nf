// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process dist_decay {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/dist_decay", mode: 'copy'

  when:
  !params.skip_dist_decay

  input:
  tuple val(sample), val(res), path(maps), val(r) 
  
  output:
  path("*_distcount.txt")
  path("*.png")


  script:
  """
  hicPlotDistVsCounts --matrices ${maps} \
                      --plotFile ${maps.baseName}_distcount.png \
  		      --outFileData ${maps.baseName}_distcount.txt
  """
}
