// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process combine_mates{
   tag "$sample = $r1_prefix + $r2_prefix"
   label 'process_low'
   publishDir "${params.outdir}/hicpro/mapping", mode: params.publish_dir_mode,
   	      saveAs: {filename -> filename.endsWith(".pairstat") ? "stats/$filename" : "$filename"}

   input:
   tuple val(sample), path(aligned_bam) 

   output:
   tuple val(oname), path("${sample}_bwt2pairs.bam"), emit:paired_bam
   tuple val(oname), path("*.pairstat"), emit:all_pairstat

   script:
   r1_bam = aligned_bam[0]
   r1_prefix = r1_bam.toString() - ~/_bwt2merged.bam$/
   r2_bam = aligned_bam[1]
   r2_prefix = r2_bam.toString() - ~/_bwt2merged.bam$/
   oname = sample.toString() - ~/(\.[0-9]+)$/

   def opts = "-t"
   if (params.keep_multi) {
     opts="${opts} --multi"
   }else if (params.min_mapq){
     opts="${opts} -q ${params.min_mapq}"
   }
   """
   mergeSAM.py -f ${r1_bam} -r ${r2_bam} -o ${sample}_bwt2pairs.bam ${opts}
   """
}
