process BUILD_CONTACT_MAPS{
  tag "$meta.id - $res"
  label 'process_highmem'
  
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
