process TRIM_READS {
  tag "$meta.id"
  label 'process_low'

  input:
  tuple val(meta), path(reads) 
  val(motif)

  output:
  tuple val(meta), path("*trimmed.fastq.gz"), emit: fastq
  path("versions.yml") , emit: versions

  script:
  """
  zcat ${reads} > tmp.fastq
  cutsite_trimming --fastq tmp.fastq \\
                   --cutsite ${motif[0]} \\
                   --out ${reads.simpleName}_trimmed.fastq
  gzip ${reads.simpleName}_trimmed.fastq
  /bin/rm -f tmp.fastq

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      python: \$(echo \$(python --version 2>&1) | sed 's/^Python//; s/ .*\$//')
  END_VERSIONS
  """
}
