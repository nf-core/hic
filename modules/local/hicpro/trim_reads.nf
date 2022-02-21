process TRIM_READS {
  tag "$meta.id"
  label 'process_low'

  input:
  tuple val(meta), path(reads) 
  val(motif)

  output:
  tuple val(meta), path("*trimmed.fastq"), emit: fastq
  path("versions.yml") , emit: versions

  script:
  """
  cutsite_trimming --fastq $reads \\
                   --cutsite ${motif} \\
                   --out ${meta.id}_trimmed.fastq
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      python: \$(echo \$(python --version 2>&1) | sed 's/^Python//; s/ .*\$//')
  END_VERSIONS
  """
}
