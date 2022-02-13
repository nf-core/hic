// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process remove_duplicates {
   tag "$sample"
   label 'process_highmem'
   publishDir "${params.outdir}/hicpro/valid_pairs", mode: params.publish_dir_mode,
               saveAs: {filename -> if (filename.endsWith("mergestat")) "stats/$filename" 
                                    else if (filename.endsWith("allValidPairs")) "$filename"}
   input:
   tuple val(sample), path(vpairs) 

   output:
   tuple val(sample), path("*.allValidPairs"), emit: ch_vpairs_cool
   path("stats/"), emit:mqc_mergestat
   path("*mergestat"), emit:all_mergestat

   script:
   if ( ! params.keep_dups ){
   """
   mkdir -p stats/${sample}

   ## Sort valid pairs and remove read pairs with same starts (i.e duplicated read pairs)
   sort -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ${vpairs} | \\
   awk -F"\\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=\$2 || c2!=\$5 || s1!=\$3 || s2!=\$6){print;c1=\$2;c2=\$5;s1=\$3;s2=\$6}' > ${sample}.allValidPairs

   echo -n "valid_interaction\t" > ${sample}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> ${sample}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> ${sample}_allValidPairs.mergestat
   cat ${sample}.allValidPairs | wc -l >> ${sample}_allValidPairs.mergestat

   ## Count short range (<20000) vs long range contacts
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> ${sample}_allValidPairs.mergestat
 
   ## For MultiQC
   mkdir -p stats/${sample} 
   cp ${sample}_allValidPairs.mergestat stats/${sample}/
   """
   }else{
   """
   cat ${vpairs} > ${sample}.allValidPairs
   echo -n "valid_interaction\t" > ${sample}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> ${sample}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> ${sample}_allValidPairs.mergestat
   cat ${sample}.allValidPairs | wc -l >> ${sample}_allValidPairs.mergestat

   ## Count short range (<20000) vs long range contacts
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> ${sample}_allValidPairs.mergestat

   ## For MultiQC
   mkdir -p stats/${sample}
   cp ${sample}_allValidPairs.mergestat stats/${sample}/
   """
   }
}
