// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process get_valid_interaction{
      tag "$sample"
      label 'process_low'
      publishDir "${params.outdir}/hicpro/valid_pairs", mode: params.publish_dir_mode,
   	      saveAs: {filename -> if (filename.endsWith("RSstat")) "stats/$filename"
                                   else if (filename.endsWith(".validPairs")) filename
                                   else if (params.save_nonvalid_pairs) filename}

      input:
      tuple val(sample), path(pe_bam) 
      path frag_path 

      output:
      tuple val(sample), path("*.validPairs"), emit:valid_pairs
      tuple val(sample), path("*.validPairs"), emit:valid_pairs_4cool
      tuple val(sample), path("*.DEPairs"), emit:de_pairs
      tuple val(sample), path("*.SCPairs"), emit:sc_pairs
      tuple val(sample), path("*.REPairs"), emit:re_pairs
      tuple val(sample), path("*.FiltPairs"), emit:filt_pairs
      tuple val(sample), path("*RSstat"), emit:all_rsstat

      script:
      if (params.split_fastq){
         sample = sample.toString() - ~/(\.[0-9]+)$/
      }

      def opts = ""
      opts += params.min_cis_dist > 0 ? " -d ${params.min_cis_dist}" : ''
      opts += params.min_insert_size > 0 ?  " -s ${params.min_insert_size}" : ''
      opts += params.max_insert_size > 0 ? " -l ${params.max_insert_size}" : ''
      opts += params.min_restriction_fragment_size > 0 ? " -t ${params.min_restriction_fragment_size}" : ''
      opts += params.max_restriction_fragment_size > 0 ? " -m ${params.max_restriction_fragment_size}" : ''
      opts += params.save_interaction_bam ? " --sam" : ''
      prefix = pe_bam.toString() - ~/.bam/
      """
      mapped_2hic_fragments.py -f ${frag_file} -r ${pe_bam} --all ${opts}
      sort -k2,2V -k3,3n -k5,5V -k6,6n -o ${prefix}.validPairs ${prefix}.validPairs
      """
}
