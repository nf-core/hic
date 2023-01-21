/*
 * COOLER MAIN WORKFLOW
 * INPUT : .pair text file with the list of valid interaction
 * OUTPUT : cooler files
 */

include { COOLER_ZOOMIFY } from '../../modules/nf-core/cooler/zoomify/main'
include { COOLER_DUMP } from '../../modules/nf-core/cooler/dump/main' 
include { COOLER_CLOAD } from '../../modules/nf-core/cooler/cload/main' 
include { COOLER_BALANCE } from '../../modules/nf-core/cooler/balance/main'
include { COOLER_MAKEBINS } from '../../modules/nf-core/cooler/makebins/main'

include { SPLIT_COOLER_DUMP } from '../../modules/local/split_cooler_dump'

// add resolution in meta
def addResolution(row) {
  def meta = [:]
  meta.id = row[0].id
  meta.resolution = row[2]
  return [meta, row[1], row[2]]
}

workflow COOLER {

  take:
  pairs // [meta, pairs, index]
  chromsize // [meta, chromsize]
  cool_bins

  main:
  ch_versions = Channel.empty()

  //*****************************************
  // EXPORT BINS

  COOLER_MAKEBINS(
    chromsize.combine(cool_bins)
  )
  ch_versions = ch_versions.mix(COOLER_MAKEBINS.out.versions)

  //*****************************************
  // BUILD COOL FILE PER RESOLUTION
  // [meta, pairs, resolution]

  COOLER_CLOAD(
    pairs.combine(cool_bins),
    chromsize.map{it -> it[1]}.collect()
  )
  ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

  // Add resolution in meta
  COOLER_CLOAD.out.cool
    .map{ it -> addResolution(it) }
    .set{ ch_cool }

  COOLER_BALANCE(
    ch_cool.map{[it[0], it[1], ""]}
  )
  ch_versions = ch_versions.mix(COOLER_BALANCE.out.versions)

  // Zoomify at minimum bin resolution
  if (!params.res_zoomify){
    ch_res_zoomify = cool_bins.min()
  }else{
    ch_res_zoomify = Channel.from(params.res_zoomify).splitCsv().flatten().unique().toInteger()
  }

  ch_cool
    .combine(ch_res_zoomify)
    .filter{ it[2] == it[3] }
    .map{ it->[it[0], it[1]] }
    .set{ ch_cool_zoomify }

  COOLER_ZOOMIFY(
    ch_cool_zoomify
  )
  ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

  //*****************************************
  // DUMP DATA
  // [meta, cool] / resolution

  COOLER_DUMP(
    COOLER_BALANCE.out.cool.map{[it[0], it[1], ""]}
  )
  ch_versions = ch_versions.mix(COOLER_DUMP.out.versions)

  SPLIT_COOLER_DUMP(
    COOLER_DUMP.out.bedpe
  )
  ch_versions = ch_versions.mix(SPLIT_COOLER_DUMP.out.versions)

  emit:
  versions = ch_versions
  cool = COOLER_BALANCE.out.cool
  mcool = COOLER_ZOOMIFY.out.mcool
}