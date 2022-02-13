// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process cooler_raw {
  tag "$sample - ${res}"
  label 'process_medium'

  publishDir "${params.outdir}/contact_maps/", mode: 'copy',
              saveAs: {filename -> filename.endsWith(".cool") ? "raw/cool/$filename" : "raw/txt/$filename"}

  input:
  tuple val(sample), path(contacts), val(res) 
  path chrsize 

  output:
  tuple val(sample), val(res), path("*cool"), emit:raw_cool_maps
  tuple path("*.bed"), path("${sample}_${res}.txt"), emit:raw_txt_maps

  script:
  """
  cooler makebins ${chrsize} ${res} > ${sample}_${res}.bed
  cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${sample}_${res}.bed ${contacts} ${sample}_${res}.cool
  cooler dump ${sample}_${res}.cool | awk '{OFS="\t"; print \$1+1,\$2+1,\$3}' > ${sample}_${res}.txt
  """
}
