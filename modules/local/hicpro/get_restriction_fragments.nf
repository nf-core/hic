process GET_RESTRICTION_FRAGMENTS {
  tag "$res_site"
  label 'process_low'
  
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
