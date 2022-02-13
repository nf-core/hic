// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process build_contact_maps{
   tag "$sample - $mres"
   label 'process_highmem'
   publishDir "${params.outdir}/hicpro/matrix/raw", mode: params.publish_dir_mode

   when:
   !params.skip_maps && params.hicpro_maps

   input:
   tuple val(sample), path(vpairs), val(mres) 
   path chrsize 

   output:
   tuple val(sample), val(mres), path("*.matrix"), path("*.bed"), emit: raw_maps_4cool
   
   script:
   """
   build_matrix --matrix-format upper  --binsize ${mres} --chrsizes ${chrsize} --ipath ${vpairs} --oprefix ${sample}_${mres}
   """
}
