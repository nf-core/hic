process GET_VALID_INTERACTION_DNASE {
  tag "$meta.id"
  label 'process_low'
  
  input:
  tuple val(meta), path(bam) 

  output:
  tuple val(meta), path("*.validPairs"), emit:valid_pairs
  tuple val(meta), path("*RSstat"), optional: true, emit:stats
  path("versions.yml"), emit: versions

  script:
  def args = task.ext.args ?: ''
  """
  mapped_2hic_dnase.py \\
    -r ${bam} \\
    ${args}

  sort -k2,2V -k3,3n -k5,5V -k6,6n -o ${bam.baseName}.validPairs ${bam.baseName}.validPairs

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
  END_VERSIONS
  """
}
