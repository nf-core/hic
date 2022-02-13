// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process bowtie2_merge_mapping_steps{
      tag "$prefix = $bam1 + $bam2"
      label 'process_medium'
      publishDir "${params.outdir}/hicpro/mapping", mode: params.publish_dir_mode,
   	      saveAs: { filename -> if (params.save_aligned_intermediates && filename.endsWith("stat")) "stats/$filename"
			else if (params.save_aligned_intermediates) filename}

      input:
      tuple val(prefix), path(bam1), path(bam2) 

      output:
      tuple val(sample), path("${prefix}_bwt2merged.bam"), emit:bwt2_merged_bam
      tuple val(oname), path("${prefix}.mapstat"), emit:all_mapstat

      script:
      sample = prefix.toString() - ~/(_R1|_R2)/
      tag = prefix.toString() =~/_R1/ ? "R1" : "R2"
      oname = prefix.toString() - ~/(\.[0-9]+)$/
      """
      samtools merge -@ ${task.cpus} \\
    	             -f ${prefix}_bwt2merged.bam \\
                     ${bam1} ${bam2}

      samtools sort -@ ${task.cpus} -m 800M \\
      	            -n  \\
	            -o ${prefix}_bwt2merged.sorted.bam \\
	            ${prefix}_bwt2merged.bam

      mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam

      echo "## ${prefix}" > ${prefix}.mapstat
      echo -n "total_${tag}\t" >> ${prefix}.mapstat
      samtools view -c ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
      echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
      echo -n "global_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam1} >> ${prefix}.mapstat
      echo -n "local_${tag}\t"  >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam2} >> ${prefix}.mapstat
      """
}
