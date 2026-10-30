/*
 * LOOPS
 * Chromatin loop calling for HiChIP / PLAC-seq data
 * INPUT : valid pairs, cool files at the loop resolution(s), chromosome sizes
 * OUTPUT : 1D anchor peaks (MACS3) and statistically significant interactions (FitHiChIP)
 */

include { PAIRS2BED             } from '../../../modules/local/pairs2bed/main'
include { MACS3_CALLPEAK        } from '../../../modules/nf-core/macs3/callpeak/main'
include { FITHICHIP_COOL2MATRIX } from '../../../modules/local/fithichip/cool2matrix'
include { FITHICHIP             } from '../../../modules/local/fithichip/fithichip'

workflow LOOPS {

    take:
    pairs     // [meta, pairs, index]
    cool      // [meta, cool, resolution]
    chromsize // [meta, chromsize]

    main:

    //
    // 1D anchors: either provided with --loops_peaks, or inferred from the 5' ends of the valid pairs,
    // the approach of FitHiChIP's PeakInferHiChIP (MACS3 on single-end reads, --nomodel --extsize 147)
    //
    if (params.loops_peaks) {
        ch_user_peaks = channel.value(file(params.loops_peaks, checkIfExists: true))
        ch_peaks = cool
            .map { meta, _cool, _res -> meta.id }
            .unique()
            .combine(ch_user_peaks)   // [id, peaks]
        ch_macs_peaks = channel.empty()
    } else {
        PAIRS2BED(
            pairs
        )

        // MACS3 effective genome size: --macs_gsize, else the total length of the chromosome sizes file
        if (params.macs_gsize) {
            ch_gsize = channel.value(params.macs_gsize)
        } else {
            ch_gsize = chromsize
                .map { _meta, sizes -> sizes.readLines().collect { it.tokenize('\t')[1].toLong() }.sum() }
                .first()
        }

        MACS3_CALLPEAK(
            PAIRS2BED.out.bed.map { meta, bed -> [meta, bed, []] },
            ch_gsize
        )
        ch_macs_peaks = MACS3_CALLPEAK.out.peak
        ch_peaks = ch_macs_peaks.map { meta, peak -> [meta.id, peak] }
    }

    //
    // FitHiChIP on the raw contact counts at the loop resolution
    //
    if (params.loops_caller == 'fithichip') {
        FITHICHIP_COOL2MATRIX(
            cool.map { meta, cool_file, _res -> [meta, cool_file] }
        )

        ch_fithichip_input = FITHICHIP_COOL2MATRIX.out.matrix
            .map { meta, bins, matrix -> [meta.id, meta, bins, matrix] }
            .combine(ch_peaks, by: 0)                       // [id, meta, bins, matrix, peaks]
            .combine(chromsize.map { it[1] })               // [id, meta, bins, matrix, peaks, chromsize]
            .map { _id, meta, bins, matrix, peaks, sizes -> [meta, bins, matrix, peaks, sizes] }

        FITHICHIP(
            ch_fithichip_input,
            params.fithichip_lower_dist,
            params.fithichip_upper_dist,
            params.fithichip_qvalue,
            params.fithichip_inttype,
            params.fithichip_bias,
            params.fithichip_p2p_background
        )
        ch_interactions = FITHICHIP.out.interactions
        ch_bedpe = FITHICHIP.out.bedpe
    } else {
        ch_interactions = channel.empty()
        ch_bedpe = channel.empty()
    }

    emit:
    peaks = ch_macs_peaks        // [meta, narrowPeak]
    interactions = ch_interactions // [meta, significant interactions bed]
    bedpe = ch_bedpe             // [meta, IGV bedpe]
}
