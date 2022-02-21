// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process bowtie2_end_to_end {
   tag "$sample"
   label 'process_medium'
   publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/mapping/bwt2_end2end" : params.outdir },
              saveAs: { filename -> if (params.save_aligned_intermediates) filename }, mode: params.publish_dir_mode

   input:
   tuple val(sample), path(reads)
   path index

   output:
   tuple val(sample), path("${prefix}_unmap.fastq"), emit: unmapped_end_to_end
   tuple val(sample), path("${prefix}.bam"), emit: end_to_end_bam

   script:
   prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
   def bwt2_opts = params.bwt2_opts_end2end
   if (!params.dnase){
   """
   INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
   bowtie2 --rg-id BMG --rg SM:${prefix} \\
	${bwt2_opts} \\
	-p ${task.cpus} \\
	-x \${INDEX} \\
	--un ${prefix}_unmap.fastq \\
 	-U ${reads} | samtools view -F 4 -bS - > ${prefix}.bam
   """
   }else{
   """
   INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
   bowtie2 --rg-id BMG --rg SM:${prefix} \\
	${bwt2_opts} \\
	-p ${task.cpus} \\
	-x \${INDEX} \\
	--un ${prefix}_unmap.fastq \\
 	-U ${reads} > ${prefix}.bam
   """
   }
}
