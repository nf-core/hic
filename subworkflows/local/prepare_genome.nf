/*
 * Prepare Annotation Genome for Hi-C data analysis
 */

include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main'
include { GET_RESTRICTION_FRAGMENTS } from '../../modules/local/hicpro/get_restriction_fragments'

workflow PREPARE_GENOME {

  take:
  fasta
  restriction_site

  main:
  ch_versions = Channel.empty()

  //***************************************
  // Bowtie Index
  if(!params.bwt2_index){
    BOWTIE2_BUILD (
      fasta
    )
    ch_index = BOWTIE2_BUILD.out.index
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
  }else{
    Channel.fromPath( params.bwt2_index , checkIfExists: true)
           .map { it -> [[:], it]}
           .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwt2_index}" }
           .set { ch_index }
  }

  //***************************************
  // Chromosome size
  if(!params.chromosome_size){
    CUSTOM_GETCHROMSIZES(
      fasta
    )
    ch_chromsize = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
  }else{
    Channel.fromPath( params.chromosome_size , checkIfExists: true)
           .map { it -> [[:], it]}
           .set {ch_chromsize} 
  }

  //***************************************
  // Restriction fragments
  if(!params.restriction_fragments && !params.dnase){
    GET_RESTRICTION_FRAGMENTS(
      fasta,
      restriction_site
    )
    ch_resfrag = GET_RESTRICTION_FRAGMENTS.out.results
    ch_versions = ch_versions.mix(GET_RESTRICTION_FRAGMENTS.out.versions)
  }else if (!params.dnase){
     Channel.fromPath( params.restriction_fragments, checkIfExists: true )
            .map{ it -> [[:], it] }
            .set {ch_resfrag}
  }else{
    ch_resfrag = Channel.empty()
  }

  emit:
  index = ch_index
  chromosome_size = ch_chromsize
  res_frag = ch_resfrag
  versions = ch_versions
}
