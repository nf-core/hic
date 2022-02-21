process GET_RESTRICTION_FRAGMENTS {
  tag "$res_site"
  label 'process_low'
  
   input:
   path fasta 
   val(res_site)

   output:
   path "*.bed", emit: results

   script:
   """
   digest_genome.py -r ${res_site} -o restriction_fragments.bed ${fasta}
   """
}
