/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hic_pipeline'

// MODULE: Local to the pipeline
include { HIC_PLOT_DIST_VS_COUNTS } from '../modules/local/hicexplorer/hicPlotDistVsCounts'

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { HICPRO } from '../subworkflows/local/hicpro'
include { PAIRTOOLS } from '../subworkflows/local/pairtools'
include { COOLER } from '../subworkflows/local/cooler'
include { COMPARTMENTS } from '../subworkflows/local/compartments'
include { TADS } from '../subworkflows/local/tads'

//****************************************
// Combine all maps resolution for downstream analysis

ch_map_res = Channel.from( params.bin_size ).splitCsv().flatten().toInteger()

if (params.res_zoomify){
    ch_zoom_res = Channel.from( params.res_zoomify ).splitCsv().flatten().toInteger()
    ch_map_res = ch_map_res.concat(ch_zoom_res)
}

if (params.res_tads && !params.skip_tads){
    ch_tads_res = Channel.from( "${params.res_tads}" ).splitCsv().flatten().toInteger()
    ch_map_res = ch_map_res.concat(ch_tads_res)
}else{
    ch_tads_res=Channel.empty()
    if (!params.skip_tads){
        log.warn "[nf-core/hic] Hi-C resolution for TADs calling not specified. See --res_tads"
    }
}

if (params.res_dist_decay && !params.skip_dist_decay){
    ch_ddecay_res = Channel.from( "${params.res_dist_decay}" ).splitCsv().flatten().toInteger()
    ch_map_res = ch_map_res.concat(ch_ddecay_res)
}else{
    ch_ddecay_res = Channel.empty()
    if (!params.skip_dist_decay){
        log.warn "[nf-core/hic] Hi-C resolution for distance decay not specified. See --res_dist_decay"
    }
}

if (params.res_compartments && !params.skip_compartments){
    ch_comp_res = Channel.from( "${params.res_compartments}" ).splitCsv().flatten().toInteger()
    ch_map_res = ch_map_res.concat(ch_comp_res)
}else{
    ch_comp_res = Channel.empty()
    if (!params.skip_compartments){
        log.warn "[nf-core/hic] Hi-C resolution for compartment calling not specified. See --res_compartments"
    }
}

ch_map_res = ch_map_res.unique()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_index
    ch_chromosome_size
    ch_res_frag
    ch_restriction_site
    ch_ligation_site

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // SUB-WORFLOW: HiC-Pro
    //
    if (params.processing == 'hicpro'){
        HICPRO (
            ch_samplesheet,
            ch_fasta,
            ch_index,
            ch_res_frag,
            ch_chromosome_size,
            ch_ligation_site,
            ch_map_res
        )
        ch_versions = ch_versions.mix(HICPRO.out.versions)
        ch_pairs = HICPRO.out.pairs
        ch_process_mqc = HICPRO.out.mqc
    }else if (params.processing == 'pairtools'){
        PAIRTOOLS(
            ch_samplesheet,
            ch_fasta,
            ch_index,
            ch_res_frag,
            ch_chromosome_size
        )
        ch_versions = ch_versions.mix(PAIRTOOLS.out.versions)
        ch_pairs = PAIRTOOLS.out.pairs
        ch_process_mqc = PAIRTOOLS.out.stats
    }

    //
    // SUB-WORKFLOW: COOLER
    //
    COOLER (
        ch_pairs,
        ch_chromosome_size,
        ch_map_res
    )
    ch_versions = ch_versions.mix(COOLER.out.versions)

    //
    // MODULE: HICEXPLORER/HIC_PLOT_DIST_VS_COUNTS
    //
    if (!params.skip_dist_decay){
        COOLER.out.cool
            .combine(ch_ddecay_res)
            .filter{ it[0].resolution == it[2] }
            .map { it -> [it[0], it[1]]}
            .set{ ch_distdecay }

        HIC_PLOT_DIST_VS_COUNTS(
            ch_distdecay
        )
        ch_versions = ch_versions.mix(HIC_PLOT_DIST_VS_COUNTS.out.versions)
    }

    //
    // SUB-WORKFLOW: COMPARTMENT CALLING
    //
    if (!params.skip_compartments){
        COOLER.out.cool
            .combine(ch_comp_res)
            .filter{ it[0].resolution == it[2] }
            .map { it -> [it[0], it[1], it[2]]}
            .set{ ch_cool_compartments }

        COMPARTMENTS (
            ch_cool_compartments,
            ch_fasta,
            ch_chromosome_size
        )
        ch_versions = ch_versions.mix(COMPARTMENTS.out.versions)
    }

    //
    // SUB-WORKFLOW : TADS CALLING
    //
    if (!params.skip_tads){
        COOLER.out.cool
            .combine(ch_tads_res)
            .filter{ it[0].resolution == it[2] }
            .map { it -> [it[0], it[1]]}
            .set{ ch_cool_tads }

        TADS(
            ch_cool_tads
        )
        ch_versions = ch_versions.mix(TADS.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    if (params.processing == 'hicpro'){
        ch_multiqc_files = ch_multiqc_files.mix(HICPRO.out.mqc)
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        Channel.value([]),
        Channel.value([])
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
