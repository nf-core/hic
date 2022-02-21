process COMBINE_MATES {
  tag "$prefix"
  label 'process_low'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*bwt2pairs.bam"), emit:bam
  tuple val(meta), path("*.pairstat"), optional:true, emit:stats

  script:
  prefix = meta.id
  def args = task.ext.args ?: ''
  """
  mergeSAM.py -f ${bam[0]} -r ${bam[1]} -o ${prefix}_bwt2pairs.bam ${args}
  """
}
