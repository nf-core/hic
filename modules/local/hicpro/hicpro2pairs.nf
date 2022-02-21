process HICPRO2PAIRS {
  tag "$meta.id"
  label 'process_medium'

  input:
  tuple val(meta), path(vpairs)
  path chrsize 

  output:
  tuple val(meta), path("*.txt.gz"), emit: pairs

  script:
  """
  ## chr/pos/strand/chr/pos/strand
  awk '{OFS="\t";print \$1,\$2,\$3,\$5,\$6,\$4,\$7}' $vpairs > contacts.txt
  gzip contacts.txt
  """
}
