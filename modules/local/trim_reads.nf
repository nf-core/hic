// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process trim_reads {
   tag "$sample"
   label 'process_low'
   publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/mapping/bwt2_trimmed" : params.outdir },
              saveAs: { filename -> if (params.save_aligned_intermediates) filename }, mode: params.publish_dir_mode
              
   when:
   !params.dnase

   input:
   tuple val(sample), path(reads) 

   output:
   tuple val(sample), path("${prefix}_trimmed.fastq"), emit:trimmed_reads

   script:
   prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
   """
   cutsite_trimming --fastq $reads \\
                    --cutsite  ${params.ligation_site} \\
                    --out ${prefix}_trimmed.fastq
   """
}
