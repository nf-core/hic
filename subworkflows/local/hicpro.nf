params.options = [:]

include { BOWTIE2_END_TO_END } from '../../modules/local/bowtie2_end_to_end' addParams( options: params.options )
include { BOWTIE2_ON_TRIMED_READS } from '../../modules/local/bowtie2_on_trimmed_reads' addParams( options: params.options )
include { BOWTIE2_MERGE_MAPPING_STEPS } from '../../modules/local/bowtie2_merge_mapping_steps' addParams( options: params.options )
include { DNASE_MAPPING_STATS } from '../../modules/local/dnase_mapping_stats' addParams( options: params.options )
include { COMBINE_MATES } from '../../modules/local/combine_mates' addParams( options: params.options )
include { GET_VALID_INTERACTION } from '../../modules/local/get_valid_interaction' addParams( options: params.options )
include { GET_VALID_INTERACTION_DNASE } from '../../modules/local/get_valid_interaction_dnase' addParams( options: params.options )
include { REMOVE_DUPLICATES } from '../../modules/local/remove_duplicates' addParams( options: params.options )
include { MERGE_STATS } from '../../modules/local/merge_stats' addParams( options: params.options )
include { BUILD_CONTACT_MAPS } from '../../modules/local/build_contact_maps' addParams( options: params.options )
include { RUN_ICE } from '../../modules/local/run_ice' addParams( options: params.options )
include { CONVERTS_TO_PAIRS } from '../../modules/local/convert_to_pairs' addParams( options: params.options )

workflow HIC_PRO {

    take:


    main:
    

    emit:
    
}
