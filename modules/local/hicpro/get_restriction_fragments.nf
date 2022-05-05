process GET_RESTRICTION_FRAGMENTS {
  tag "$res_site"
  label 'process_low'

  conda (params.enable_conda ? "conda-forge::python=3.7.6  conda-forge::numpy=1.18.1" : null)

  input:
  path fasta 
  val(res_site)

  output:
  path "*.bed", emit: results
  path("versions.yml"), emit: versions

  script:
  """
  digest_genome.py -r ${res_site} -o restriction_fragments.bed ${fasta}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
  END_VERSIONS
  """
}
