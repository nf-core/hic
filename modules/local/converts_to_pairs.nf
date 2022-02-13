// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process convert_to_pairs {
   tag "$sample"
   label 'process_medium'

   when:
   !params.skip_maps

   input:
   tuple val(sample), path(vpairs) 
   path chrsize 

   output:
   tuple val(sample), path("*.txt.gz"), emit: cool_build_zoom

   script:
   """
   ## chr/pos/strand/chr/pos/strand
   awk '{OFS="\t";print \$1,\$2,\$3,\$5,\$6,\$4,\$7}' $vpairs > contacts.txt
   gzip contacts.txt
   """
}
