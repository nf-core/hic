process BUILD_CONTACT_MAPS{
  tag "$meta.id - $res"
  label 'process_highmem'

  conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
      'ubuntu:20.04' }"

  input:
  tuple val(meta), path(vpairs), val(res) 
  path chrsize 

  output:
  tuple val(meta), val(res), path("*.matrix"), path("*.bed"), emit: maps
   
  script:
  """
  build_matrix \\
    --matrix-format upper  \\
    --binsize ${res} \\
    --chrsizes ${chrsize} \\
    --ifile ${vpairs} \\
    --oprefix ${meta.id}_${res}
  """
}
