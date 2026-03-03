/*
 * HiC-Pro mapping
 * From the raw sequencing reads to a paired-end bam file
 */

include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { TRIM_READS } from '../../modules/local/hicpro/trim_reads'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_TRIMMED } from '../../modules/nf-core/bowtie2/align/main'
include { MERGE_BOWTIE2 } from '../../modules/local/hicpro/bowtie2_merge'
include { COMBINE_MATES} from '../../modules/local/hicpro/combine_mates'
include { MAPPING_STATS_DNASE } from '../../modules/local/hicpro/dnase_mapping_stats'


workflow HICPRO_MAPPING {

    take:
    reads // [meta, read1, read2]
    fasta // [meta, fasta]
    index // [meta, path]
    ligation_site // value

    main:
    ch_versions = Channel.empty()

    // Align each mates separetly and add mates information in [meta]
    ch_reads_r1 = reads
        .map{ meta, fastq ->
            def newMeta = [ id: meta.id, single_end:true, chunk:meta.chunk, part:meta.part, mates:'R1' ]
            [ newMeta, fastq[0] ]
        }

    ch_reads_r2 = reads
        .map{ meta, fastq ->
            def newMeta = [ id: meta.id, single_end:true, chunk:meta.chunk, part:meta.part, mates:'R2' ]
            [ newMeta, fastq[1] ]
        }

    ch_reads = ch_reads_r1.concat(ch_reads_r2)

    // bowtie2 - save_unaligned=true - sort_bam=false
    BOWTIE2_ALIGN(
        ch_reads,
        index.collect(),
        fasta.collect(),
        true,
        false
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    if (!params.no_digestion){

        // trim reads
        TRIM_READS(
            BOWTIE2_ALIGN.out.fastq,
            ligation_site.collect()
        )
        ch_versions = ch_versions.mix(TRIM_READS.out.versions)

        // bowtie2 on trimmed reads - save_unaligned=false - sort_bam=false
        BOWTIE2_ALIGN_TRIMMED(
            TRIM_READS.out.fastq,
            index.collect(),
            fasta.collect(),
            false,
            false
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_TRIMMED.out.versions)

        // Merge the two mapping steps
        ch_bowtie2_align = BOWTIE2_ALIGN.out.bam
            .combine(BOWTIE2_ALIGN_TRIMMED.out.bam, by:[0])

        MERGE_BOWTIE2(
            ch_bowtie2_align
        )
        ch_versions = ch_versions.mix(MERGE_BOWTIE2.out.versions)
        ch_mapping_stats = MERGE_BOWTIE2.out.stats

        ch_bams = MERGE_BOWTIE2.out.bam
            .map{ meta, bam ->
                def newMeta = [ id: meta.id, single_end:false, chunk:meta.chunk, part:meta.part ]
                [ newMeta, bam ]
            }.groupTuple()


    }else{

        MAPPING_STATS_DNASE(
            BOWTIE2_ALIGN.out.aligned
        )
        ch_versions = ch_versions.mix(MAPPING_STATS_DNASE.out.versions)
        ch_mapping_stats = MAPPING_STATS_DNASE.out.stats

        ch_bams = BOWTIE2_ALIGN.out.aligned
            .map{ meta, bam ->
                def newMeta = [ id: meta.id, single_end: false, chunk:meta.chunk, part:meta.part ]
                [ newMeta, bam ]
            }.groupTuple()
    }

    COMBINE_MATES (
        ch_bams
    )
    ch_versions = ch_versions.mix(COMBINE_MATES.out.versions)

    emit:
    versions = ch_versions
    bam = COMBINE_MATES.out.bam
    mapstats = ch_mapping_stats
    pairstats = COMBINE_MATES.out.stats
}
