process ICE_NORMALIZATION{
  tag "$rmaps"
  label 'process_highmem'

  input:
  tuple val(meta), val(res), path(rmaps), path(bed) 

  output:
  tuple val(meta), val(res), path("*iced.matrix"), path(bed), emit:maps
  path ("*.biases"), emit:bias

  script:
  prefix = rmaps.toString() - ~/(\.matrix)?$/
  """
  ice --filter_low_counts_perc ${params.ice_filter_low_count_perc} \
      --results_filename ${prefix}_iced.matrix \
      --filter_high_counts_perc ${params.ice_filter_high_count_perc} \
      --max_iter ${params.ice_max_iter} --eps ${params.ice_eps} --remove-all-zeros-loci --output-bias 1 --verbose 1 ${rmaps}
   """
}
