process GET_VALID_INTERACTION {
  tag "$meta.id"
  label 'process_low'
  
  input:
  tuple val(meta), path(bam) 
  path(resfrag)

  output:
  tuple val(meta), path("*.validPairs"), emit:valid_pairs
  tuple val(meta), path("*.DEPairs"), optional:true, emit:de_pairs
  tuple val(meta), path("*.SCPairs"), optional: true, emit:sc_pairs
  tuple val(meta), path("*.REPairs"), optional: true, emit:re_pairs
  tuple val(meta), path("*.FiltPairs"), optional: true, emit:filt_pairs
  tuple val(meta), path("*RSstat"), optional: true, emit:stats
  path("versions.yml"), emit: versions

  script:
  def args = task.ext.args ?: ''
  """
  mapped_2hic_fragments.py \\
    -f ${resfrag} \\
    -r ${bam} \\
    --all \\
    ${args}

  sort -k2,2V -k3,3n -k5,5V -k6,6n -o ${bam.baseName}.validPairs ${bam.baseName}.validPairs

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
  END_VERSIONS
  """
}
