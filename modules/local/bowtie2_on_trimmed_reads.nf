// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process bowtie2_on_trimmed_reads {
   tag "$sample"
   label 'process_medium'
   publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/mapping/bwt2_trimmed" : params.outdir },
   	      saveAs: { filename -> if (params.save_aligned_intermediates) filename }, mode: params.publish_dir_mode

   when:
   !params.dnase

   input:
   tuple val(sample), path(reads) 
   path index 

   output:
   tuple val(sample), path("${prefix}_trimmed.bam"), emit:trimmed_bam

   script:
   prefix = reads.toString() - ~/(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
   """
   INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
   bowtie2 --rg-id BMG --rg SM:${prefix} \\
           ${params.bwt2_opts_trimmed} \\
           -p ${task.cpus} \\
           -x \${INDEX} \\
           -U ${reads} | samtools view -bS - > ${prefix}_trimmed.bam
   """
}
