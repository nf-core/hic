/*
 * PAIRTOOLS
 * MAIN WORKFLOW
 * From the raw sequencing reads to the list of valid interactions
 */
  
//include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { PAIRTOOLS_DEDUP } from '../../modules/nf-core/pairtools/dedup/main'
include { PAIRTOOLS_PARSE } from '../../modules/nf-core/pairtools/parse/main'
include { PAIRTOOLS_RESTRICT } from '../../modules/nf-core/pairtools/restrict/main'
include { PAIRTOOLS_SELECT } from '../../modules/nf-core/pairtools/select/main'
include { PAIRTOOLS_SORT } from '../../modules/nf-core/pairtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { PAIRIX } from '../../modules/nf-core/pairix/main'

include { PAIRTOOLS_MERGE } from '../../modules/local/pairtools/pairtools_merge'
include { PAIRTOOLS_SPLIT } from '../../modules/local/pairtools/pairtools_split'
include { PAIRTOOLS_STATS } from '../../modules/local/pairtools/pairtools_stats'


workflow PAIRTOOLS {

  take:
  reads // [meta, read1, read2]
  index // [meta2, path]
  frag // path
  chrsize // path

  main:
  ch_versions = Channel.empty()

  BWA_MEM(
    reads,
    index,
    Channel.value([])
  )

  PAIRTOOLS_PARSE(
    BWA_MEM.out.bam,
    chrsize
  )

  PAIRTOOLS_RESTRICT(
    PAIRTOOLS_PARSE.out.pairsam,
    frag.map{it->it[1]}
  )

  ch_pairsam = params.dnase ? PAIRTOOLS_PARSE.out.pairsam : PAIRTOOLS_RESTRICT.out.restrict
  PAIRTOOLS_SORT(
    ch_pairsam
  )

  ch_valid_pairs = PAIRTOOLS_SORT.out.sorted
    .view()
    .map{ meta, pairs -> 
      def newMeta = [ id: meta.id, single_end: meta.single_end, part:meta.part ]
      [ groupKey(newMeta, meta.part), pairs ]
    }
    .groupTuple()
    .branch {
      single: it[0].part <=1
      multiple: it[0].part > 1
    }

  PAIRTOOLS_MERGE(
    ch_valid_pairs.multiple
  )

  PAIRTOOLS_DEDUP(
    PAIRTOOLS_MERGE.out.pairs.mix(ch_valid_pairs.single)
  )

  ch_pairsam2split = params.keep_dups ? PAIRTOOLS_MERGE.out.pairs.mix(ch_valid_pairs.single) : PAIRTOOLS_DEDUP.out.pairs
  PAIRTOOLS_SPLIT(
    ch_pairsam2split
  )

  SAMTOOLS_INDEX(
    PAIRTOOLS_SPLIT.out.bam
  )

  PAIRTOOLS_SELECT(
    PAIRTOOLS_SPLIT.out.pairs
  )

  PAIRTOOLS_STATS(
    PAIRTOOLS_SPLIT.out.pairs
  )

  PAIRIX(
    PAIRTOOLS_SELECT.out.selected
  )

  emit:
  versions = ch_versions
  pairs = PAIRIX.out.index
  bam = PAIRTOOLS_SPLIT.out.bam.join(SAMTOOLS_INDEX.out.bai)
  stats = PAIRTOOLS_STATS.out.stats.map{it->it[1]}
}
