include { BOWTIE2_ALIGN } from '../../modules/nf-core/modules/bowtie2/align/main'
include { TRIM_READS } from '../../modules/local/hicpro/trim_reads'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TRIMMED } from '../../modules/nf-core/modules/bowtie2/align/main'
include { MERGE_BOWTIE2 } from '../../modules/local/hicpro/bowtie2_merge'
include { COMBINE_MATES} from '../../modules/local/hicpro/combine_mates'

//include { BOWTIE2_ON_TRIMED_READS } from '../../modules/local/bowtie2_on_trimmed_reads' addParams( options: params.options )
//include { BOWTIE2_MERGE_MAPPING_STEPS } from '../../modules/local/bowtie2_merge_mapping_steps' addParams( options: params.options )
//include { DNASE_MAPPING_STATS } from '../../modules/local/dnase_mapping_stats' addParams( options: params.options )
//include { COMBINE_MATES } from '../../modules/local/combine_mates' addParams( options: params.options )
//include { GET_VALID_INTERACTION } from '../../modules/local/get_valid_interaction' addParams( options: params.options )
//include { GET_VALID_INTERACTION_DNASE } from '../../modules/local/get_valid_interaction_dnase' addParams( options: params.options )
//include { REMOVE_DUPLICATES } from '../../modules/local/remove_duplicates' addParams( options: params.options )
//include { MERGE_STATS } from '../../modules/local/merge_stats' addParams( options: params.options )
//include { BUILD_CONTACT_MAPS } from '../../modules/local/build_contact_maps' addParams( options: params.options )
//include { RUN_ICE } from '../../modules/local/run_ice' addParams( options: params.options )
//include { CONVERTS_TO_PAIRS } from '../../modules/local/convert_to_pairs' addParams( options: params.options )

// Paired-end to Single-end 
def pairToSingle(row, mates) {
  def meta = [:]
  meta.id = row[0].id
  meta.single_end = true
  meta.mates = mates
  def array = []
  if (mates == "R1") {
    return [meta, [ row[1][0]] ]
  }else if (mates == "R2"){
    return [meta, [ row[1][1]] ]
  }
}

def singleToPair(row){
  def meta = [:]
  meta.id = row[0].id
  meta.single_end = false
  return [ meta, row[1] ]
}

workflow HICPRO_MAPPING {

  take:
  reads // [meta, read1, read2]
  index
  ligation_site

  main:
  ch_versions = Channel.empty()
 
  // Align each mates separetly
  ch_reads_r1 = reads.map{it -> pairToSingle(it,"R1")}
  ch_reads_r2 = reads.map{pairToSingle(it,"R2")}
  ch_reads = ch_reads_r1.concat(ch_reads_r2)

  // bowtie2
  BOWTIE2_ALIGN(
    ch_reads,
    index.collect(),
    Channel.value(true).collect()
  )
  ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

  // trim reads
  TRIM_READS(
    BOWTIE2_ALIGN.out.fastq,
    ligation_site.collect()
  )
  ch_versions = ch_versions.mix(TRIM_READS.out.versions)

  // bowtie2 on trimmed reads
  BOWTIE2_ALIGN_TRIMMED(
    TRIM_READS.out.fastq,
    index.collect(),
    Channel.value(false).collect()
  )
  ch_versions = ch_versions.mix(BOWTIE2_ALIGN_TRIMMED.out.versions)

  // Merge the two mapping steps
  BOWTIE2_ALIGN.out.bam
    .combine(BOWTIE2_ALIGN_TRIMMED.out.bam, by:[0])
    .view()
    .set { ch_bowtie2_align}

  MERGE_BOWTIE2(
    ch_bowtie2_align
  )
  //TODO ch_versions = ch_versions.mix(MERGE_BOWTIE2.out.versions)

  // Combine mates
  MERGE_BOWTIE2.out.bam
    .map { singleToPair(it) }
    .groupTuple()
    .view()
    .set {ch_bams}

  COMBINE_MATES (
    ch_bams
  )
  //TODO ch_versions = ch_versions.mix(COMBINE_MATES.out.versions)

  emit:
  versions = ch_versions
  bam = COMBINE_MATES.out.bam
  mapstats = MERGE_BOWTIE2.out.stats
  pairstats = COMBINE_MATES.out.stats
}
