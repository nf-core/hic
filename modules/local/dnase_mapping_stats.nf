// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process dnase_mapping_stats{
      tag "$sample = $bam"
      label 'process_medium'
      publishDir "${params.outdir}/hicpro/mapping",  mode: params.publish_dir_mode, 
   	      saveAs: { filename -> if (params.save_aligned_intermediates && filename.endsWith("stat")) "stats/$filename"
	                else if (params.save_aligned_intermediates) filename}

      input:
      tuple val(prefix), path(bam) 

      output:
      tuple val(sample), path(bam), emit:bwt2_merged_bam
      tuple val(oname), path("${prefix}.mapstat"), emit:all_mapstat

      script:
      sample = prefix.toString() - ~/(_R1|_R2)/
      tag = prefix.toString() =~/_R1/ ? "R1" : "R2"
      oname = prefix.toString() - ~/(\.[0-9]+)$/
      """
      echo "## ${prefix}" > ${prefix}.mapstat
      echo -n "total_${tag}\t" >> ${prefix}.mapstat
      samtools view -c ${bam} >> ${prefix}.mapstat
      echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam} >> ${prefix}.mapstat
      echo -n "global_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam} >> ${prefix}.mapstat
      echo -n "local_${tag}\t0"  >> ${prefix}.mapstat
      """
}
