// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process run_ice{
   tag "$rmaps"
   label 'process_highmem'
   publishDir "${params.outdir}/hicpro/matrix/iced", mode: params.publish_dir_mode

   when:
   !params.skip_maps && !params.skip_balancing && params.hicpro_maps

   input:
   tuple val(sample), val(res), path(rmaps), path(bed) 

   output:
   tuple val(sample), val(res), path("*iced.matrix"), path(bed), emit:hicpro_iced_maps
   path ("*.biases"), emit:hicpro_iced_bias

   script:
   prefix = rmaps.toString() - ~/(\.matrix)?$/
   """
   ice --filter_low_counts_perc ${params.ice_filter_low_count_perc} \
   --results_filename ${prefix}_iced.matrix \
   --filter_high_counts_perc ${params.ice_filter_high_count_perc} \
   --max_iter ${params.ice_max_iter} --eps ${params.ice_eps} --remove-all-zeros-loci --output-bias 1 --verbose 1 ${rmaps}
   """
}
