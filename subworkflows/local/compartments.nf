include { COOLTOOLS_EIGSCIS } from '../../modules/local/cooltools/eigscis'

workflow COMPARTMENTS {

  take:
  cool
  fasta
  chrsize

  main:
  ch_versions = Channel.empty()

  COOLTOOLS_EIGSCIS(
    cool,
    fasta.map{it -> it[1]}.collect(),
    chrsize.map{it -> it[1]}.collect()
  )
  ch_versions = ch_versions.mix(COOLTOOLS_EIGSCIS.out.versions)

  emit:
  versions = ch_versions
  compartments = COOLTOOLS_EIGSCIS.out.results
}