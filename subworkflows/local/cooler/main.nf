/*
 * COOLER MAIN WORKFLOW
 * INPUT : .pair text file with the list of valid interaction
 * OUTPUT : cooler files
 */

include { COOLER_ZOOMIFY } from '../../../modules/nf-core/cooler/zoomify/main'
include { COOLER_DUMP } from '../../../modules/nf-core/cooler/dump/main'
include { COOLER_CLOAD } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_BALANCE } from '../../../modules/nf-core/cooler/balance/main'
include { COOLER_MAKEBINS } from '../../../modules/nf-core/cooler/makebins/main'

include { FILTER_CHROMSIZE } from '../../../modules/local/filter_chromsize'
include { SPLIT_COOLER_DUMP } from '../../../modules/local/split_cooler_dump'

workflow COOLER {

    take:
    ch_pairs // [meta, pairs, index]
    ch_chromsize // [meta, chromsize]
    ch_cools_bins
    ch_zoom_res

    main:

    //*****************************************
    // FILTER CHROMOSOMES ON SIZE

    if( params.min_size ) {
        ch_chromsize = ch_chromsize | FILTER_CHROMSIZE
    }

    //*****************************************
    // EXPORT BINS

    COOLER_MAKEBINS(
        ch_chromsize.combine(ch_cools_bins)
    )

    //*****************************************
    // BUILD COOL FILE PER RESOLUTION
    COOLER_CLOAD(
        ch_pairs.collect(),
        ch_chromsize.collect(),
        "pairs",
        ch_cools_bins
    )

    // Add resolution in meta
    COOLER_CLOAD.out.cool
        .map { meta, file ->
            def id = (file.baseName =~ /(\d+)(?!.*\d)/)[0][1]
            [meta + [resolution: id.toInteger()], file]
        }
        .set { ch_cool }

    COOLER_BALANCE(
        ch_cool.map{[it[0], it[1], ""]}
    )

    ch_cool
        .combine(ch_zoom_res)
        .filter{ meta, cool, zoom_res -> meta.resolution == zoom_res }
        .map{ it->[it[0], it[1]] }
        .set{ ch_cool_zoomify }

    COOLER_ZOOMIFY(
        ch_cool_zoomify
    )

    //*****************************************
    // DUMP DATA
    // [meta, cool] / resolution

    COOLER_DUMP(
        COOLER_BALANCE.out.cool.map{[it[0], it[1], ""]}
    )

    SPLIT_COOLER_DUMP(
        COOLER_DUMP.out.bedpe
    )

    emit:
    cool = COOLER_BALANCE.out.cool
    mcool = COOLER_ZOOMIFY.out.mcool
}
