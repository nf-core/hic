process MERGE_VALID_INTERACTION {
   tag "$prefix"
   label 'process_highmem'

   input:
   tuple val(meta), path(vpairs) 

   output:
   tuple val(meta), path("*.allValidPairs"), emit: valid_pairs
   path("stats/"), emit:mqc
   tuple val(meta), path("*mergestat"), emit:stats

   script:
   prefix = meta.id
   if ( ! params.keep_dups ){
   """
   mkdir -p stats/${prefix}

   ## Sort valid pairs and remove read pairs with same starts (i.e duplicated read pairs)
   sort -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ${vpairs} | \\
   awk -F"\\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=\$2 || c2!=\$5 || s1!=\$3 || s2!=\$6){print;c1=\$2;c2=\$5;s1=\$3;s2=\$6}' > ${prefix}.allValidPairs

   echo -n "valid_interaction\t" > ${prefix}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> ${prefix}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> ${prefix}_allValidPairs.mergestat
   cat ${prefix}.allValidPairs | wc -l >> ${prefix}_allValidPairs.mergestat

   ## Count short range (<20000) vs long range contacts
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${prefix}.allValidPairs >> ${prefix}_allValidPairs.mergestat
 
   ## For MultiQC
   mkdir -p stats/${prefix} 
   cp ${prefix}_allValidPairs.mergestat stats/${prefix}/
   """
   }else{
   """
   cat ${vpairs} > ${prefix}.allValidPairs
   echo -n "valid_interaction\t" > ${prefix}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> ${prefix}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> ${prefix}_allValidPairs.mergestat
   cat ${prefix}.allValidPairs | wc -l >> ${prefix}_allValidPairs.mergestat

   ## Count short range (<20000) vs long range contacts
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${prefix}.allValidPairs >> ${prefix}_allValidPairs.mergestat

   ## For MultiQC
   mkdir -p stats/${prefix}
   cp ${prefix}_allValidPairs.mergestat stats/${prefix}/
   """
   }
}
