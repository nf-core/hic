process GET_CHROMSIZE {
  tag "$fasta"
  label 'process_low'

  input:
  path fasta 

  output:
  path "*.size", emit: results

  script:
  """
  samtools faidx ${fasta}
  cut -f1,2 ${fasta}.fai > chrom.size
  """
}
