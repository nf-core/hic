include { CALL_COMPARTMENTS } from '../../modules/local/cooltools/eigs-cis'

workflow COMPARTMENTS {

  take:
  cool
  fasta
  chrsize

  main:
  ch_versions = Channel.empty()

  CALL_COMPARTMENTS (
    cool,
    fasta.collect(),
    chrsize.collect()
  )
  ch_versions = ch_versions.mix(CALL_COMPARTMENTS.out.versions)

  emit:
  versions = ch_versions
  compartments = CALL_COMPARTMENTS.out.results
}