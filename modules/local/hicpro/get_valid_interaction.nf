process GET_VALID_INTERACTION {
  tag "$meta.id"
  label 'process_low'
  
  input:
  tuple val(meta), path(bam) 
  path(resfrag)

  output:
  tuple val(meta), path("*.validPairs"), emit:valid_pairs
  tuple val(meta), path("*.DEPairs"), emit:de_pairs
  tuple val(meta), path("*.SCPairs"), emit:sc_pairs
  tuple val(meta), path("*.REPairs"), emit:re_pairs
  tuple val(meta), path("*.FiltPairs"), emit:filt_pairs
  tuple val(meta), path("*RSstat"), emit:stats

  script:
  if (params.split_fastq){
     sample = sample.toString() - ~/(\.[0-9]+)$/
  }
  def args = task.ext.args ?: ''
  """
  mapped_2hic_fragments.py \\
    -f ${resfrag} \\
    -r ${bam} \\
    --all \\
    ${args}

  sort -k2,2V -k3,3n -k5,5V -k6,6n -o ${bam.baseName}.validPairs ${bam.baseName}.validPairs
  """
}
