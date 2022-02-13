// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process compartment_calling {
  tag "$sample - $res"
  label 'process_medium'
  publishDir "${params.outdir}/compartments", mode: 'copy'

  when:
  !params.skip_compartments

  input:
  tuple val(sample), val(res), path(cool), val(r) 
  path(fasta) 
  path(chrsize) 

  output:
  path("*compartments*") optional true, emit:out_compartments

  script:
  """
  cooltools genome binnify --all-names ${chrsize} ${res} > genome_bins.txt
  cooltools genome gc genome_bins.txt ${fasta} > genome_gc.txt 
  cooltools call-compartments --contact-type cis -o ${sample}_compartments ${cool}
  awk -F"\t" 'NR>1{OFS="\t"; if(\$6==""){\$6=0}; print \$1,\$2,\$3,\$6}' ${sample}_compartments.cis.vecs.tsv | sort -k1,1 -k2,2n > ${sample}_compartments.cis.E1.bedgraph
  """
}
