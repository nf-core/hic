/*
 * Prepare Annotation Genome for Hi-C data analysis
 */

include { BOWTIE2_BUILD } from '../../../modules/nf-core/bowtie2/build/main.nf'
include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { GET_RESTRICTION_FRAGMENTS } from '../../../modules/local/hicpro/get_restriction_fragments'

workflow PREPARE_GENOME {

    take:
    fasta
    bwt2_index
    bwa_index

    main:
    ch_versions = channel.empty()

    //
    // Fasta reference genome
    //
    def genomeName = params.genome ?: fasta.substring(fasta.lastIndexOf(File.separator)+1)
    ch_fasta = channel.fromPath( fasta )
        .ifEmpty { exit 1, "Genome index: Fasta file not found: ${fasta}" }
        .map{it->[[id:genomeName],it]}

    //
    // Bowtie index
    //
    if (params.processing == "hicpro"){
        if(!bwt2_index){
            BOWTIE2_BUILD (
                ch_fasta
            )
            ch_index = BOWTIE2_BUILD.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }else{
            ch_index = channel.fromPath( bwt2_index , checkIfExists: true)
                .map { it -> [[:], it]}
                .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwt2_index}" }
        }
    }

    //
    // Bwa-mem index
    //
    if (params.processing == "pairtools"){
        if(!bwa_index){
            BWA_INDEX (
                ch_fasta
            )
            ch_index = BWA_INDEX.out.index
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }else{
            ch_index = channel.fromPath( bwa_index , checkIfExists: true)
                .map { it -> [[:], it]}
                .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwa_index}" }
        }
    }

    //
    // Chromosome size
    //
    if(!params.chromosome_size){
        SAMTOOLS_FAIDX(
            ch_fasta, 
            ch_index, 
            true
            )
            ch_chromsize = SAMTOOLS_FAIDX.out.sizes
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }else{
        ch_chromsize = channel.fromPath( params.chromosome_size , checkIfExists: true)
            .map { it -> [[:], it]}
    }


    //
    // Digestion parameters
    //
    if (params.digestion){
        restriction_site = params.digestion ? params.digest[ params.digestion ].restriction_site ?: false : false
        ch_restriction_site = channel.value(restriction_site)
        ligation_site = params.digestion ? params.digest[ params.digestion ].ligation_site ?: false : false
        ch_ligation_site = channel.value(ligation_site)
    }else if (params.restriction_site && params.ligation_site){
        ch_restriction_site = channel.value(params.restriction_site)
        ch_ligation_site = channel.value(params.ligation_site)
    }else if (params.no_digestion){
        ch_restriction_site = channel.empty()
        ch_ligation_site = channel.empty()
    }else{
        exit 1, "Ligation motif not found. Please either use the `--digestion` parameters or specify the `--restriction_site` and `--ligation_site`. For DNase/Micro-C Hi-C, please use '--no_digestion' option"
    }


    //
    // Restriction fragments
    //
    if(!params.restriction_fragments && !params.no_digestion){
        GET_RESTRICTION_FRAGMENTS(
            ch_fasta,
            restriction_site
        )
        ch_resfrag = GET_RESTRICTION_FRAGMENTS.out.results
        ch_versions = ch_versions.mix(GET_RESTRICTION_FRAGMENTS.out.versions)
    }else if (!params.no_digestion){
        channel.fromPath( params.restriction_fragments, checkIfExists: true )
            .map { it -> [[:], it] }
            .set {ch_resfrag}
    }else{
        ch_resfrag = channel.empty()
    }

    emit:
    fasta = ch_fasta
    index = ch_index
    chromosome_size = ch_chromsize
    res_frag = ch_resfrag
    restriction_site = ch_restriction_site
    ligation_site = ch_ligation_site
    versions = ch_versions
}
