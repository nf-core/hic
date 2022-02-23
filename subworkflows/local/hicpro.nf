/*
 * HICPRO MAIN WORKFLOW
 * INPUT :  paired-end sequencing data
 * OUTPUT : .pairs file with the list of valid interaction
 */
  
include { HICPRO_MAPPING } from './hicpro_mapping'
include { GET_VALID_INTERACTION } from '../../modules/local/hicpro/get_valid_interaction'
include { MERGE_VALID_INTERACTION } from '../../modules/local/hicpro/merge_valid_interaction'
include { HICPRO2PAIRS } from '../../modules/local/hicpro/hicpro2pairs'
include { BUILD_CONTACT_MAPS } from '../../modules/local/hicpro/build_contact_maps'
include { ICE_NORMALIZATION } from '../../modules/local/hicpro/run_ice'

workflow HICPRO {

  take:
  reads // [meta, read1, read2]
  index
  fragments
  chrsize
  ligation_site
  map_res

  main:
  ch_versions = Channel.empty()

  // fastq to paired-end bam
  HICPRO_MAPPING(
    reads,
    index,
    ligation_site
  )

  // get valid interaction
  GET_VALID_INTERACTION (
    HICPRO_MAPPING.out.bam,
    fragments    
  )

  // merge valid interactions and remove duplicates
  MERGE_VALID_INTERACTION (
    GET_VALID_INTERACTION.out.valid_pairs
  )

  // convert to pairs
  HICPRO2PAIRS (
    MERGE_VALID_INTERACTION.out.valid_pairs,
    chrsize
  )

  //merge stats
  // TODO

  if (params.hicpro_maps){
    
    //build_contact_maps
    BUILD_CONTACT_MAPS(
      MERGE_VALID_INTERACTION.out.valid_pairs.combine(map_res),
      chrsize.collect()
    )
    
    // run_ice
    ICE_NORMALIZATION(
      BUILD_CONTACT_MAPS.out.maps
    )
  }

  emit:
  versions = ch_versions
  pairs = HICPRO2PAIRS.out.pairs
}
