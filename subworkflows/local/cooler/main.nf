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
    pairs // [meta, pairs, index]
    chromsize // [meta, chromsize]
    cool_bins

    main:

    //*****************************************
    // FILTER CHROMOSOMES ON SIZE

    if( params.min_size ) {
        chromsize = chromsize | FILTER_CHROMSIZE
    }

    //*****************************************
    // EXPORT BINS

    COOLER_MAKEBINS(
        chromsize.combine(cool_bins)
    )

    //*****************************************
    // BUILD COOL FILE PER RESOLUTION
    pairs_res = pairs.combine(cool_bins)

    cload_inputs = pairs_res.multiMap { meta, pairs_files, index, cool_bin ->
        pairs: [meta, pairs_files, index]
        res: cool_bin
    }

    COOLER_CLOAD(cload_inputs.pairs, chromsize.first(), "pairs", cload_inputs.res)


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

    // Zoomify at minimum bin resolution
    if (!params.res_zoomify){
        ch_res_zoomify = cool_bins.min()
    }else{
        ch_res_zoomify = channel.from(params.res_zoomify).splitCsv().flatten().unique().toInteger()
    }

    COOLER_BALANCE.out.cool
        .combine(ch_res_zoomify)
        .filter{ it -> it[0].resolution == it[2] }
        .map{ it -> [it[0], it[1]] }
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
