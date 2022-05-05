process COMBINE_MATES {
  tag "$prefix"
  label 'process_low'

  conda (params.enable_conda ? "conda-forge::python=3.7.6  bioconda::pysam=0.15.4" : null)

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*bwt2pairs.bam"), emit:bam
  tuple val(meta), path("*.pairstat"), optional:true, emit:stats
  path("versions.yml"), emit: versions

  script:
  prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  mergeSAM.py -f ${bam[0]} -r ${bam[1]} -o ${prefix}_bwt2pairs.bam ${args}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
  END_VERSIONS
  """
}
