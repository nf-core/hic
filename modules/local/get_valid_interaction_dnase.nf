// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process get_valid_interaction_dnase{
      tag "$sample"
      label 'process_low'
      publishDir "${params.outdir}/hicpro/valid_pairs", mode: params.publish_dir_mode,
   	      saveAs: {filename -> if (filename.endsWith("RSstat")) "stats/$filename" 
                                   else filename}

      input:
      tuple val(sample), path(pe_bam) 

      output:
      tuple val(sample), path("*.validPairs"), emit:valid_pairs
      tuple val(sample), path("*.validPairs"), emit:valid_pairs_4cool
      tuple val(sample), path("*RSstat"), emit:all_rsstat

      script:
      if (params.split_fastq){
         sample = sample.toString() - ~/(\.[0-9]+)$/
      }

      opts = params.min_cis_dist > 0 ? " -d ${params.min_cis_dist}" : ''
      prefix = pe_bam.toString() - ~/.bam/
      """
      mapped_2hic_dnase.py -r ${pe_bam} ${opts}
      sort -k2,2V -k3,3n -k5,5V -k6,6n -o ${prefix}.validPairs ${prefix}.validPairs
      """
}
