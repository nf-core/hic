// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process cooler_balance {
  tag "$sample - ${res}"
  label 'process_medium'

  publishDir "${params.outdir}/contact_maps/", mode: 'copy',
              saveAs: {filename -> filename.endsWith(".cool") ? "norm/cool/$filename" : "norm/txt/$filename"}

  when:
  !params.skip_balancing

  input:
  tuple val(sample), val(res), path(cool) 
  path chrsize 

  output:
  tuple val(sample), val(res), path("${sample}_${res}_norm.cool"), emit:balanced_cool_maps
  path("${sample}_${res}_norm.txt"), emit:norm_txt_maps

  script:
  """
  cp ${cool} ${sample}_${res}_norm.cool
  cooler balance ${sample}_${res}_norm.cool -p ${task.cpus} --force
  cooler dump ${sample}_${res}_norm.cool --balanced --na-rep 0 | awk '{OFS="\t"; print \$1+1,\$2+1,\$4}' > ${sample}_${res}_norm.txt
  """
}
