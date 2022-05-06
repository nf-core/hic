/*
 * HICPRO
 * MAIN WORKFLOW
 * From the raw sequencing reads to the list of valid interactions
 */
  
include { HICPRO_MAPPING } from './hicpro_mapping'
include { GET_VALID_INTERACTION } from '../../modules/local/hicpro/get_valid_interaction'
include { GET_VALID_INTERACTION_DNASE } from '../../modules/local/hicpro/get_valid_interaction_dnase'
include { MERGE_VALID_INTERACTION } from '../../modules/local/hicpro/merge_valid_interaction'
include { MERGE_STATS } from '../../modules/local/hicpro/merge_stats'
include { HICPRO2PAIRS } from '../../modules/local/hicpro/hicpro2pairs'
include { BUILD_CONTACT_MAPS } from '../../modules/local/hicpro/build_contact_maps'
include { ICE_NORMALIZATION } from '../../modules/local/hicpro/run_ice'

// Remove meta.chunks
def removeChunks(row){
  meta = row[0].clone()
  meta.remove('chunk')
  return [meta, row[1]]
}

workflow HICPRO {

  take:
  reads // [meta, read1, read2]
  index // path
  fragments // path
  chrsize // path
  ligation_site // value
  map_res // values

  main:
  ch_versions = Channel.empty()

  // Fastq to paired-end bam
  HICPRO_MAPPING(
    reads,
    index,
    ligation_site
  )
  ch_versions = ch_versions.mix(HICPRO_MAPPING.out.versions)

  //***************************************
  // DIGESTION PROTOCOLS

  if (!params.dnase){
    GET_VALID_INTERACTION (
      HICPRO_MAPPING.out.bam,
      fragments.collect()
    )
    ch_versions = ch_versions.mix(GET_VALID_INTERACTION.out.versions)
    ch_valid_pairs = GET_VALID_INTERACTION.out.valid_pairs
    ch_valid_stats = GET_VALID_INTERACTION.out.stats

  }else{

  //****************************************
  // DNASE-LIKE PROTOCOLS

    GET_VALID_INTERACTION_DNASE (
      HICPRO_MAPPING.out.bam
    )
    ch_versions = ch_versions.mix(GET_VALID_INTERACTION_DNASE.out.versions)
    ch_valid_pairs = GET_VALID_INTERACTION_DNASE.out.valid_pairs
    ch_valid_stats = GET_VALID_INTERACTION_DNASE.out.stats
  }
  

  //**************************************
  // MERGE AND REMOVE DUPLICATES
  
  //if (params.split_fastq){
  ch_valid_pairs = ch_valid_pairs.map{ it -> removeChunks(it)}.groupTuple()
  ch_hicpro_stats = HICPRO_MAPPING.out.mapstats.map{it->removeChunks(it)}.groupTuple()
                      .concat(HICPRO_MAPPING.out.pairstats.map{it->removeChunks(it)}.groupTuple(),
		        ch_valid_stats.map{it->removeChunks(it)}.groupTuple())
  //}else{
  //  ch_hicpro_stats = HICPRO_MAPPING.out.mapstats.groupTuple()
  //                      .concat(HICPRO_MAPPING.out.pairstats.groupTuple(),
  //                              ch_valid_stats.groupTuple())
  //}

  MERGE_VALID_INTERACTION (
    ch_valid_pairs
  )
  ch_versions = ch_versions.mix(MERGE_VALID_INTERACTION.out.versions)


  ch_hicpro_stats.view()
  MERGE_STATS(
    ch_hicpro_stats
  )
  ch_versions = ch_versions.mix(MERGE_STATS.out.versions)

  //***************************************
  // CONVERTS TO PAIRS
  HICPRO2PAIRS (
    MERGE_VALID_INTERACTION.out.valid_pairs,
    chrsize.collect()
  )
  ch_versions = ch_versions.mix(HICPRO2PAIRS.out.versions)

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
    ch_versions = ch_versions.mix(ICE_NORMALIZATION.out.versions)

  }else{
    ch_hicpro_raw_maps = Channel.empty()
    ch_hicpro_iced_maps = Channel.empty()
  }

  emit:
  versions = ch_versions
  pairs = HICPRO2PAIRS.out.pairs
  mqc = MERGE_VALID_INTERACTION.out.mqc.concat(MERGE_STATS.out.mqc)
  raw_maps = ch_hicpro_raw_maps
  iced_maps = ch_hicpro_iced_maps
}
