/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowHic.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Digestion parameters
if (params.digestion){
  restriction_site = params.digestion ? params.digest[ params.digestion ].restriction_site ?: false : false
  ch_restriction_site = Channel.value(restriction_site)

  ligation_site = params.digestion ? params.digest[ params.digestion ].ligation_site ?: false : false
  ch_ligation_site = Channel.value(ligation_site)
}else{
  ch_restriction_site = Channel.empty()
  ch_ligation_site = Channel.empty()
}

// Resolutions for contact maps
ch_map_res = Channel.from( params.bin_size ).splitCsv().flatten()
if (params.res_tads && !params.skip_tads){
  Channel.from( "${params.res_tads}" )
    .splitCsv()
    .flatten()
    .set {ch_tads_res}
  ch_map_res = ch_map_res.concat(ch_tads_res)
}else{
  ch_tads_res=Channel.empty()
  if (!params.skip_tads){
    log.warn "[nf-core/hic] Hi-C resolution for TADs calling not specified. See --res_tads" 
  }
}

if (params.res_dist_decay && !params.skip_dist_decay){
  Channel.from( "${params.res_dist_decay}" )
    .splitCsv()
    .flatten()
    .set {ch_ddecay_res}
   ch_map_res = ch_map_res.concat(ch_ddecay_res)
}else{
  ch_ddecay_res = Channel.create()
  if (!params.skip_dist_decay){
    log.warn "[nf-core/hic] Hi-C resolution for distance decay not specified. See --res_dist_decay" 
  }
}

if (params.res_compartments && !params.skip_compartments){
  //Channel.fromPath( params.fasta )
  //  .ifEmpty { exit 1, "Compartments calling: Fasta file not found: ${params.fasta}" }
  //  .set { fasta_for_compartments }
  Channel.from( "${params.res_compartments}" )
    .splitCsv()
    .flatten()
    .set {ch_comp_res}
   ch_map_res = ch_map_res.concat(ch_comp_res)
}else{
  //fasta_for_compartments = Channel.empty()
  ch_comp_res = Channel.create()
  if (!params.skip_compartments){
    log.warn "[nf-core/hic] Hi-C resolution for compartment calling not specified. See --res_compartments" 
  }
}

ch_map_res = ch_map_res.unique()
/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
//def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
//include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'
//include { OUTPUT_DOCUMENTATION } from '../modules/local/output_documentation'
include { HIC_PLOT_DIST_VS_COUNTS } from '../modules/local/hicexplorer/hicPlotDistVsCounts' 
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { HICPRO } from '../subworkflows/local/hicpro'
include { COOLER } from '../subworkflows/local/cooler'
include { COMPARTMENTS } from '../subworkflows/local/compartments'
include { TADS } from '../subworkflows/local/tads'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//def multiqc_options   = modules['multiqc']
//multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'
//include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
  CHANNELS
========================================================================================
*/

Channel.fromPath( params.fasta )
       .ifEmpty { exit 1, "Genome index: Fasta file not found: ${params.fasta}" }
       .set { ch_fasta }

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow HIC {

  ch_software_versions = Channel.empty()

  //
  // SUBWORKFLOW: Read in samplesheet, validate and stage input files
  //
  INPUT_CHECK (
    ch_input
  )

  //
  // SUBWORKFLOW: Prepare genome annotation
  //
  PREPARE_GENOME(
    ch_fasta,
    ch_restriction_site
  )

  //
  // MODULE: Run FastQC
  //
  FASTQC (
    INPUT_CHECK.out.reads
  )

  //
  // SUB-WORFLOW: HiC-Pro
  //
  HICPRO (
    INPUT_CHECK.out.reads,
    PREPARE_GENOME.out.index,
    PREPARE_GENOME.out.res_frag,
    PREPARE_GENOME.out.chromosome_size,
    ch_ligation_site,
    ch_map_res
  )

  //
  // SUB-WORKFLOW: COOLER
  //
  COOLER (
    HICPRO.out.pairs,
    PREPARE_GENOME.out.chromosome_size,
    ch_map_res
  )

  //
  // MODULE: HICEXPLORER/HIC_PLOT_DIST_VS_COUNTS
  //
  if (!params.skip_dist_decay){
    COOLER.out.cool
      .combine(ch_ddecay_res)
      .filter{ it[1] == it[3] }
      .map { it -> [it[0], it[2]]}
      .set{ ch_distdecay }

    HIC_PLOT_DIST_VS_COUNTS(
      ch_distdecay
    )
  }

  //
  // SUB-WORKFLOW: COMPARTMENT CALLING
  //
  if (!params.skip_compartments){
    COOLER.out.cool
      .combine(ch_comp_res)
      .filter{ it[1] == it[3] }
      .map { it -> [it[0], it[1], it[2]]}
      .set{ ch_cool_compartments }

    COMPARTMENTS(
      ch_cool_compartments,
      ch_fasta,
      PREPARE_GENOME.out.chromosome_size
    )
  }

  //
  // SUB-WORKFLOW : TADS CALLING
  //
  if (!params.skip_tads){
    COOLER.out.cool
      .combine(ch_tads_res)
      .filter{ it[1] == it[3] }
      .map { it -> [it[0], it[2]]}
      .set{ ch_cool_tads }
                                                                                                                                                                                                            
    TADS(
      ch_cool_tads
    )
  }

  //
  // MODULE: MultiQC
  //
  workflow_summary    = WorkflowHic.paramsSummaryMultiqc(workflow, summary_params)
  ch_workflow_summary = Channel.value(workflow_summary)

  ch_multiqc_files = Channel.empty()
  ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
  ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
  //ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
  //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

  //MULTIQC (
  //  ch_multiqc_files.collect()
  //)
  //multiqc_report       = MULTIQC.out.report.toList()
  multiqc_report = []
  //ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
  if (params.email || params.email_on_fail) {
      NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
  }
  NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
