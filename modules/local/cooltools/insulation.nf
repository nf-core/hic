/*
 * Cooltools - diamond-insulation
 */
 
process INSULATION {
  label 'process_medium'

  input:
  tuple val(meta), path(cool)

  output:
  path("*tsv"), emit:results
  path("versions.yml"), emit:versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  cooltools insulation ${cool} ${args} > ${prefix}_insulation.tsv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      cooltools: \$(cooltools --version 2>&1 | sed 's/cooltools, version //')
  END_VERSIONS
  """
}
