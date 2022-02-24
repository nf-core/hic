// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process cooler_zoomify {
   tag "$sample"
   label 'process_medium'
   publishDir "${params.outdir}/contact_maps/norm/mcool", mode: 'copy'

   when:
   !params.skip_mcool

   input:
   tuple val(sample), path(contacts)  
   path chrsize 

   output:
   path("*mcool"), emit:mcool_maps

   script:
   """
   cooler makebins ${chrsize} ${params.res_zoomify} > bins.bed
   cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 bins.bed ${contacts} ${sample}.cool
   cooler zoomify --nproc ${task.cpus} --balance ${sample}.cool
   """
}
