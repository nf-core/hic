/*
 * HICPRO
 * MAIN WORKFLOW
 * From the raw sequencing reads to the list of valid interactions
 */

include { HICPRO_MAPPING } from '../hicpro_mapping'
include { GET_VALID_INTERACTION } from '../../../modules/local/hicpro/get_valid_interaction'
include { GET_VALID_INTERACTION_DNASE } from '../../../modules/local/hicpro/get_valid_interaction_dnase'
include { MERGE_VALID_INTERACTION } from '../../../modules/local/hicpro/merge_valid_interaction'
include { MERGE_STATS } from '../../../modules/local/hicpro/merge_stats'
include { HICPRO2PAIRS } from '../../../modules/local/hicpro/hicpro2pairs'
include { BUILD_CONTACT_MAPS } from '../../../modules/local/hicpro/build_contact_maps'
include { ICE_NORMALIZATION } from '../../../modules/local/hicpro/run_ice'

// Remove meta.chunks
def removeChunks(row){
    def meta = row[0].clone()
    meta.remove('chunk')
    return [meta, row[1]]
}

workflow HICPRO {

    take:
    reads // [meta, read1, read2]
    fasta // [meta, fasta]
    index // path
    fragments // path
    chrsize // path
    ligation_site // value
    map_res // values

    main:

    // Fastq to paired-end bam
    HICPRO_MAPPING(
        reads,
        fasta,
        index,
        ligation_site
    )

    //***************************************
    // DIGESTION PROTOCOLS

    if (!params.no_digestion){
        GET_VALID_INTERACTION (
            HICPRO_MAPPING.out.bam,
            fragments.collect()
        )
        ch_valid_pairs = GET_VALID_INTERACTION.out.valid_pairs
        ch_valid_stats = GET_VALID_INTERACTION.out.stats

    }else{

        //****************************************
        // DNASE-LIKE PROTOCOLS

        GET_VALID_INTERACTION_DNASE (
            HICPRO_MAPPING.out.bam
        )
        ch_valid_pairs = GET_VALID_INTERACTION_DNASE.out.valid_pairs
        ch_valid_stats = GET_VALID_INTERACTION_DNASE.out.stats
    }


    //**************************************
    // MERGE AND REMOVE DUPLICATES

    //if (params.split_fastq){
    ch_valid_pairs = ch_valid_pairs
        .map{ meta, valid_pairs ->
            def newMeta = [ id: meta.id, single_end: meta.singleEnd, part:meta.part ]
            [ groupKey(newMeta, meta.part), valid_pairs ]
        }.groupTuple()

    MERGE_VALID_INTERACTION (
        ch_valid_pairs
    )


    ch_hicpro_mappingstats = HICPRO_MAPPING.out.mapstats
        .map{ meta, stats ->
            def newMeta = [ id: meta.id, single_end: meta.singleEnd, part:meta.part ]
            [ groupKey(newMeta, meta.part), stats ]
        }.groupTuple()

    ch_hicpro_pairstats = HICPRO_MAPPING.out.pairstats
        .map{ meta, stats ->
            def newMeta = [ id: meta.id, single_end: meta.singleEnd, part:meta.part ]
            [ groupKey(newMeta, meta.part), stats ]
        }.groupTuple()

    ch_hicpro_validstats = ch_valid_stats
        .map{ meta, stats ->
            def newMeta = [ id: meta.id, single_end: meta.singleEnd, part:meta.part ]
            [ groupKey(newMeta, meta.part), stats ]
        }.groupTuple()

    MERGE_STATS(
        ch_hicpro_mappingstats.concat(ch_hicpro_pairstats, ch_hicpro_validstats)
    )

    //***************************************
    // CONVERTS TO PAIRS

    HICPRO2PAIRS (
        MERGE_VALID_INTERACTION.out.valid_pairs,
        chrsize.collect()
    )

    //***************************************
    // CONTACT MAPS

    if (params.hicpro_maps){

        //build_contact_maps
        BUILD_CONTACT_MAPS(
            MERGE_VALID_INTERACTION.out.valid_pairs.combine(map_res),
            chrsize.collect()
        )
        ch_hicpro_raw_maps = BUILD_CONTACT_MAPS.out.maps

        // run_ice
        ICE_NORMALIZATION(
            BUILD_CONTACT_MAPS.out.maps
        )
        ch_hicpro_iced_maps = ICE_NORMALIZATION.out.maps

    }else{
        ch_hicpro_raw_maps = channel.empty()
        ch_hicpro_iced_maps = channel.empty()
    }

    emit:
    pairs = HICPRO2PAIRS.out.pairs
    mqc = MERGE_VALID_INTERACTION.out.mqc.concat(MERGE_STATS.out.mqc)
    raw_maps = ch_hicpro_raw_maps
    iced_maps = ch_hicpro_iced_maps
}
