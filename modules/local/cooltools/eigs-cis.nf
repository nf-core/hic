/*
 * cooltools - call_compartments
 */
 
process CALL_COMPARTMENTS {
  label 'process_medium'

  input:
  tuple val(meta), path(cool), val(resolution)
  path(fasta) 
  path(chrsize) 

  output:
  path("*compartments*"), emit: results
  path("versions.yml"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  cooltools genome binnify --all-names ${chrsize} ${resolution} > genome_bins.txt
  cooltools genome gc genome_bins.txt ${fasta} > genome_gc.txt 
  cooltools eigs-cis ${args} -o ${prefix}_compartments ${cool}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      cooltools: \$(cooltools --version 2>&1 | sed 's/cooletools, version //')
  END_VERSIONS
  """
}
