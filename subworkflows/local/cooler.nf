/*
 * COOLER MAIN WORKFLOW
 * INPUT : pair text file with the list of valid interaction
 * OUTPUT : cooler files
 */

include { COOLER_CLOAD } from '../../modules/nf-core/modules/cooler/cload/main'
include { COOLER_DUMP } from '../../modules/nf-core/modules/cooler/dump/main'
include { COOLER_ZOOMIFY } from '../../modules/nf-core/modules/cooler/zoomify/main'

include { COOLER_BALANCE } from '../../modules/local/cooler/balance'
include { SPLIT_COOLER_DUMP } from '../../modules/local/split_cooler_dump'
include { COOLER_MAKEBINS } from '../../modules/local/cooler/makebins'

workflow COOLER {

  take:
  pairs // [meta, pairs, index]
  chromsize
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
    chromsize.collect()
  )
  ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

  COOLER_BALANCE(
    COOLER_CLOAD.out.cool
  )
  ch_versions = ch_versions.mix(COOLER_BALANCE.out.versions)

  // Zoomify at minimum bin resolution
  if (!params.res_zoomify){
    ch_res_zoomify = cool_bins.min()
  }else{
    ch_res_zoomify = params.res_zoomify
  }
  COOLER_CLOAD.out.cool
    .combine(ch_res_zoomify)
    .filter{ it [1] == it[3] }
    .map{it->[it[0], it[2]]}
    .set{ch_cool_zoomify}

  COOLER_ZOOMIFY(
    ch_cool_zoomify
  )
  ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

  //*****************************************
  // DUMP DATA
  // [meta, cool] / resolution

  COOLER_DUMP(
    COOLER_BALANCE.out.cool.map{[it[0], "", it[2]]}
  )
  ch_versions = ch_versions.mix(COOLER_DUMP.out.versions)

  //COOLER_DUMP(
  //  COOLER_ZOOMIFY.out.mcool.combine(cool_bins).map{it->[it[0], it[2], it[1]]}
  //)

  SPLIT_COOLER_DUMP(
    COOLER_DUMP.out.bedpe
  )
  ch_versions = ch_versions.mix(SPLIT_COOLER_DUMP.out.versions)

  emit:
  versions = ch_versions
  cool = COOLER_BALANCE.out.cool
  mcool = COOLER_ZOOMIFY.out.mcool
}