include { HICDETECTLOOPS                       } from '../../../modules/local/hicexplorer/hicdetectloops'
include { HICPLOTMATRIX as HICPLOTMATRIX_LOOPS } from '../../../modules/local/hicexplorer/hicplotmatrix'

workflow LOOP_CALLING {

    take:
    ch_cool
    loop_caller

    main:
    ch_versions = channel.empty()

    if (loop_caller =~ 'hicexplorer'){

        HICDETECTLOOPS (
            ch_cool
        )
        ch_loops = HICDETECTLOOPS.out.bedgraph
        ch_versions = ch_versions.mix(HICDETECTLOOPS.out.versions)

        // Create channel: [meta, matrix, tads, loops, bigwig]
        ch_cool
            .map { meta, matrix -> [ meta.id, meta, matrix ] }
            .combine(ch_loops.map{ meta, loops -> [meta.id, loops] }, by: 0)
            .map { id, meta, matrix, loops -> 
                [meta, matrix, [], loops, []]
            }
            .set{ ch_cool_plotmatrix }

        HICPLOTMATRIX_LOOPS (
            ch_cool_plotmatrix
        )
        ch_versions = ch_versions.mix(HICPLOTMATRIX_LOOPS.out.versions)

    }

    emit:
    loops = ch_loops
    plots = HICPLOTMATRIX_LOOPS.out.png
    versions = ch_versions
}
