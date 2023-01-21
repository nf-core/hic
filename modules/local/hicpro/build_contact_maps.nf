process BUILD_CONTACT_MAPS{
  tag "${meta.id}"
  label 'process_high_memory'

  conda "conda-forge::sed=4.7"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
      'ubuntu:20.04' }"

  input:
  tuple val(meta), path(vpairs), val(resolution) 
  tuple val(meta2), path(chrsize) 

  output:
  tuple val(meta), val(resolution), path("*.matrix"), path("*.bed"), emit: maps
   
  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  build_matrix \\
    --matrix-format upper  \\
    --binsize ${resolution} \\
    --chrsizes ${chrsize} \\
    --ifile ${vpairs} \\
    --oprefix ${prefix}
  """
}
