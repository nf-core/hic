/*
 * COOLER MAIN WORKFLOW
 * INPUT : pair text file with the list of valid interaction
 * OUTPUT : cooler files
 */

include { COOLER_CLOAD } from '../../modules/nf-core/modules/cooler/cload/main'
include { COOLER_DUMP } from '../../modules/nf-core/modules/cooler/dump/main'
include { COOLER_ZOOMIFY } from '../../modules/nf-core/modules/cooler/zoomify/main'

include { COOLER_BALANCE } from '../../modules/local/balance'
include { SPLIT_COOLER_DUMP } from '../../modules/local/split_cooler_dump'
include { COOLER_MAKEBINS } from '../../modules/local/makebins'

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
    
  //*****************************************
  // BUILD COOL FILE PER RESOLUTION
  // [meta, pairs, resolution]

  COOLER_CLOAD(
    pairs.combine(cool_bins),
    chromsize.collect()
  )

  COOLER_BALANCE(
    COOLER_CLOAD.out.cool
  )

  // Zoomify at minimum bin resolution
  COOLER_CLOAD.out.cool
    .combine(cool_bins.min())
    .filter{ it [1] == it[3] }
    .map{it->[it[0], it[2]]}
    .set{ch_cool_zoomify}

  COOLER_ZOOMIFY(
    ch_cool_zoomify
  )

  //*****************************************
  // DUMP DATA
  // [meta, cool] / resolution

  COOLER_DUMP(
    COOLER_BALANCE.out.cool.map{[it[0], "", it[2]]}
  )

  //COOLER_DUMP(
  //  COOLER_ZOOMIFY.out.mcool.combine(cool_bins).map{it->[it[0], it[2], it[1]]}
  //)

  SPLIT_COOLER_DUMP(
    COOLER_DUMP.out.bedpe
  )

  emit:
  versions = ch_versions
  cool = COOLER_BALANCE.out.cool
  mcool = COOLER_ZOOMIFY.out.mcool
}