/*
 * PAIRTOOLS
 * MAIN WORKFLOW
 * From the raw sequencing reads to the list of valid interactions
 */

//include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { PAIRTOOLS_DEDUP } from '../../modules/nf-core/pairtools/dedup/main'
//include { PAIRTOOLS_PARSE } from '../../modules/nf-core/pairtools/parse/main'
include { PAIRTOOLS_RESTRICT } from '../../modules/nf-core/pairtools/restrict/main'
include { PAIRTOOLS_SELECT } from '../../modules/nf-core/pairtools/select/main'
include { PAIRTOOLS_SORT } from '../../modules/nf-core/pairtools/sort/main'
include { PAIRTOOLS_MERGE } from '../../modules/nf-core/pairtools/merge/main'
include { PAIRTOOLS_STATS } from '../../modules/nf-core/pairtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { PAIRIX } from '../../modules/nf-core/pairix/main'

//include { PAIRTOOLS_MERGE } from '../../modules/local/pairtools/pairtools_merge'
include { PAIRTOOLS_SPLIT } from '../../modules/local/pairtools/pairtools_split'
//include { PAIRTOOLS_STATS } from '../../modules/local/pairtools/pairtools_stats'
include { PAIRTOOLS_PARSE } from '../../modules/local/pairtools/pairtools_parse'

workflow PAIRTOOLS {

    take:
    reads // [meta, read1, read2]
    fasta // [meta, fasta]
    index // [meta2, path]
    frag // path
    chrsize // path

    main:
    ch_versions = Channel.empty()

    BWA_MEM(
        reads,
        index.collect(),
        fasta.collect(),
        Channel.value([])
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    PAIRTOOLS_PARSE(
        BWA_MEM.out.bam,
        chrsize.collect()
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_PARSE.out.versions)

    PAIRTOOLS_RESTRICT(
        PAIRTOOLS_PARSE.out.pairsam,
        frag.map{it->it[1]}.collect()
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_RESTRICT.out.versions)

    ch_pairsam = params.no_digestion ? PAIRTOOLS_PARSE.out.pairsam : PAIRTOOLS_RESTRICT.out.restrict
    PAIRTOOLS_SORT(
        ch_pairsam
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_SORT.out.versions)

    ch_valid_pairs = PAIRTOOLS_SORT.out.sorted
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
    ch_versions = ch_versions.mix(PAIRTOOLS_MERGE.out.versions)

    // Separate pairs/bam files
    PAIRTOOLS_SPLIT(
        PAIRTOOLS_MERGE.out.pairs.mix(ch_valid_pairs.single)
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_SPLIT.out.versions)

    // Manage BAM files
    SAMTOOLS_SORT(
        PAIRTOOLS_SPLIT.out.bam,
        fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    SAMTOOLS_FLAGSTAT(
        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    PAIRTOOLS_DEDUP(
        PAIRTOOLS_SPLIT.out.pairs
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_DEDUP.out.versions)

    ch_pairselect = params.keep_dups ? PAIRTOOLS_SPLIT.out.pairs : PAIRTOOLS_DEDUP.out.pairs
    PAIRTOOLS_SELECT(
        ch_pairselect
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_SELECT.out.versions)

    PAIRTOOLS_STATS(
        PAIRTOOLS_SELECT.out.selected
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_STATS.out.versions)

    PAIRIX(
        PAIRTOOLS_SELECT.out.selected
    )
    ch_versions = ch_versions.mix(PAIRIX.out.versions)

    emit:
    versions = ch_versions
    pairs = PAIRIX.out.index
    bam = PAIRTOOLS_SPLIT.out.bam.join(SAMTOOLS_INDEX.out.bai)
    stats = PAIRTOOLS_STATS.out.stats.map{it->it[1]}
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
}
