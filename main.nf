#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/hic
========================================================================================
 nf-core/hic Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/hic
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/hic --input '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --input [file]                            Path to input data (must be surrounded with quotes)
      --genome [str]                            Name of iGenomes reference
      -profile [str]                            Configuration profile to use. Can use multiple (comma separated)
                                                Available: conda, docker, singularity, awsbatch, test and more.

    References                                  If not specified in the configuration file or you wish to overwrite any of the references.
      --bwt2_index [file]                       Path to Bowtie2 index
      --fasta [file]                            Path to Fasta reference

    Digestion Hi-C                              If not specified in the configuration file or you wish to set up specific digestion protocol
      --digestion [str]                         Digestion Hi-C. Name of restriction enzyme used for digestion pre-configuration. Default: 'hindiii'
      --ligation_site [str]                     Ligation motifs to trim (comma separated) if not available in --digestion. Default: false
      --restriction_site [str]                  Cutting motif(s) of restriction enzyme(s) (comma separated) if not available in --digestion. Default: false
      --chromosome_size [file]                  Path to chromosome size file
      --restriction_fragments [file]            Path to restriction fragment file (bed)
      --save_reference [bool]                   Save reference genome to output folder. Default: False

    DNase Hi-C
      --dnase [bool]                            Run DNase Hi-C mode. All options related to restriction fragments are not considered. Default: False
      --min_cis_dist [int]                      Minimum intra-chromosomal distance to consider. Default: None 

    Alignments
      --bwt2_opts_end2end [str]                 Options for bowtie2 end-to-end mappinf (first mapping step). See hic.config for default.
      --bwt2_opts_trimmed [str]                 Options for bowtie2 mapping after ligation site trimming. See hic.config for default.
      --min_mapq [int]                          Minimum mapping quality values to consider. Default: 10
      --keep_multi [bool]                       Keep multi-mapped reads (--min_mapq is ignored). Default: false
      --keep_dups [bool]                        Keep duplicates. Default: false
      --save_aligned_intermediates [bool]       Save intermediates alignment files. Default: False
      --split_fastq [bool]                      Split fastq files in reads chunks to speed up computation. Default: false
      --fastq_chunks_size [int]                 Size of read chunks if split_fastq is true. Default: 20000000

    Valid Pairs Detection
      --min_restriction_fragment_size [int]     Minimum size of restriction fragments to consider. Default: None
      --max_restriction_fragment_size [int]     Maximum size of restriction fragments to consider. Default: None
      --min_insert_size [int]                   Minimum insert size of mapped reads to consider. Default: None
      --max_insert_size [int]                   Maximum insert size of mapped reads to consider. Default: None
      --save_interaction_bam [bool]             Save BAM file with interaction tags (dangling-end, self-circle, etc.). Default: False

    Contact maps
      --bin_size [str]                          Bin size for contact maps (comma separated). Default: '1000000,500000'
      --ice_max_iter [int]                      Maximum number of iteration for ICE normalization. Default: 100
      --ice_filter_low_count_perc [float]       Percentage of low counts columns/rows to filter before ICE normalization. Default: 0.02
      --ice_filter_high_count_perc [float]      Percentage of high counts columns/rows to filter before ICE normalization. Default: 0
      --ice_eps [float]                         Convergence criteria for ICE normalization. Default: 0.1

    Workflow
      --skip_maps [bool]                        Skip generation of contact maps. Useful for capture-C. Default: False
      --skip_balancing [bool]                   Skip contact maps normalization. Default: False
      --skip_mcool                              Skip mcool file generation. Default: False
      --skip_dist_decay                         Skip distance decay quality control. Default: False
      --skip_tads [bool]                        Skip TADs calling. Default: False
      --skip_multiqc [bool]                     Skip MultiQC. Default: False

    Other options:
      --outdir [file]                           The output directory where the results will be saved
      --publish_dir_mode [str]                  Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. Default: None
      --email_on_fail [email]                   Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]            Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                               Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic. Default: None

    AWSBatch options:
      --awsqueue [str]                          The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]                         The AWS Region for your AWS Batch job to run on
      --awscli [str]                            Path to the AWS CLI tool
    """.stripIndent()
}

/**********************************************************
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (params.digest && params.digestion && !params.digest.containsKey(params.digestion)) {
   exit 1, "Unknown digestion protocol. Currently, the available digestion options are ${params.digest.keySet().join(", ")}. Please set manually the '--restriction_site' and '--ligation_site' parameters."
}

params.restriction_site = params.digestion ? params.digest[ params.digestion ].restriction_site ?: false : false
params.ligation_site = params.digestion ? params.digest[ params.digestion ].ligation_site ?: false : false

// Check Digestion or DNase Hi-C mode
if (!params.dnase && !params.ligation_site) {
   exit 1, "Ligation motif not found. Please either use the `--digestion` parameters or specify the `--restriction_site` and `--ligation_site`. For DNase Hi-C, please use '--dnase' option"
}

// Reference index path configuration
params.bwt2_index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * input read files
 */

if (params.input_paths){

   raw_reads = Channel.create()
   raw_reads_2 = Channel.create()

   Channel
      .from( params.input_paths )
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0] + "_R1", a[1][0]), tuple(a[0] + "_R2", a[1][1])] }

}else{
   raw_reads = Channel.create()
   raw_reads_2 = Channel.create()

   if ( params.split_fastq ){
      Channel
         .fromFilePairs( params.input, flat:true )
         .splitFastq( by: params.fastq_chunks_size, pe:true, file: true, compress:true)
         .separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0] + "_R1", a[1]), tuple(a[0] + "_R2", a[2])] }
   }else{
      Channel
         .fromFilePairs( params.input )
	 .separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0] + "_R1", a[1][0]), tuple(a[0] + "_R2", a[1][1])] }
   }
}

// Update sample name if splitFastq is used
def updateSampleName(x) {
   if ((matcher = x[1] =~ /\s*(\.[\d]+).fastq.gz/)) {
        res = matcher[0][1]
   }
   return [x[0] + res, x[1]]
}

if (params.split_fastq ){
  raw_reads = raw_reads.concat( raw_reads_2 ).map{it -> updateSampleName(it)}.dump(tag:'input')
}else{
  raw_reads = raw_reads.concat( raw_reads_2 ).dump(tag:'input')
}

/*
 * Other input channels
 */

// Reference genome
if ( params.bwt2_index ){
   lastPath = params.bwt2_index.lastIndexOf(File.separator)
   bwt2_dir =  params.bwt2_index.substring(0,lastPath+1)
   bwt2_base = params.bwt2_index.substring(lastPath+1)

   Channel.fromPath( bwt2_dir , checkIfExists: true)
      .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwt2_index}" }
      .into { bwt2_index_end2end; bwt2_index_trim }

}
else if ( params.fasta ) {
   lastPath = params.fasta.lastIndexOf(File.separator)
   fasta_base = params.fasta.substring(lastPath+1)
   bwt2_base = fasta_base.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?(\.fsa)?$/

   Channel.fromPath( params.fasta )
	.ifEmpty { exit 1, "Genome index: Fasta file not found: ${params.fasta}" }
        .set { fasta_for_index }
}
else {
   exit 1, "No reference genome specified!"
}

// Chromosome size
if ( params.chromosome_size ){
   Channel.fromPath( params.chromosome_size , checkIfExists: true)
         .into {chrsize; chrsize_build; chrsize_raw; chrsize_balance; chrsize_zoom}
}
else if ( params.fasta ){
   Channel.fromPath( params.fasta )
	.ifEmpty { exit 1, "Chromosome sizes: Fasta file not found: ${params.fasta}" }
       	.set { fasta_for_chromsize }
}
else {
   exit 1, "No chromosome size specified!"
}

// Restriction fragments
if ( params.restriction_fragments ){
   Channel.fromPath( params.restriction_fragments, checkIfExists: true )
      .set {res_frag_file}
}
else if ( params.fasta && params.restriction_site ){
   Channel.fromPath( params.fasta )
           .ifEmpty { exit 1, "Restriction fragments: Fasta file not found: ${params.fasta}" }
           .set { fasta_for_resfrag }
}
else if (! params.dnase) {
    exit 1, "No restriction fragments file specified!"
}

// Resolutions for contact maps
map_res = Channel.from( params.bin_size ).splitCsv().flatten()
all_res = params.bin_size
if (params.res_tads && !params.skip_tads){
  Channel.from( "${params.res_tads}" )
    .splitCsv()
    .flatten()
    .into {tads_bin; tads_res_hicexplorer; tads_res_insulation}
    map_res = map_res.concat(tads_bin)
    all_res = all_res + ',' + params.res_tads
}else{
  tads_res_hicexplorer=Channel.empty()
  tads_res_insulation=Channel.empty()
  tads_bin=Channel.empty()
  if (!params.skip_tads){
    log.warn "[nf-core/hic] Hi-C resolution for TADs calling not specified. See --res_tads" 
  }
}

if (params.res_dist_decay && !params.skip_dist_decay){
  Channel.from( "${params.res_dist_decay}" )
    .splitCsv()
    .flatten()
    .into {ddecay_res; ddecay_bin }
    map_res = map_res.concat(ddecay_bin)
    all_res = all_res + ',' + params.res_dist_decay
}else{
  ddecay_res = Channel.create()
  ddecay_bin = Channel.create()
  if (!params.skip_dist_decay){
    log.warn "[nf-core/hic] Hi-C resolution for distance decay not specified. See --res_dist_decay" 
  }
}

map_res
  .unique()
  .into { map_res_summary; map_res; map_res_cool }

/**********************************************************
 * SET UP LOGS
 */

// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Input']            = params.input
summary['splitFastq']       = params.split_fastq
if (params.split_fastq)
   summary['Read chunks Size'] = params.fastq_chunks_size
summary['Fasta Ref']        = params.fasta
if (params.restriction_site){
   summary['Digestion']        = params.digestion
   summary['Restriction Motif']= params.restriction_site
   summary['Ligation Motif']   = params.ligation_site
   summary['Min Fragment Size']= params.min_restriction_fragment_size
   summary['Max Fragment Size']= params.max_restriction_fragment_size
   summary['Min Insert Size']  = params.min_insert_size
   summary['Max Insert Size']  = params.max_insert_size
}else{
   summary['DNase Mode']    = params.dnase
   summary['Min CIS dist']  = params.min_cis_dist
}
summary['Min MAPQ']         = params.min_mapq
summary['Keep Duplicates']  = params.keep_dups ? 'Yes' : 'No'
summary['Keep Multihits']   = params.keep_multi ? 'Yes' : 'No'
summary['Maps resolution']  = all_res
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-hic-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/hic Workflow Summary'
    section_href: 'https://github.com/nf-core/hic'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */

process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

   output:
   file 'software_versions_mqc.yaml' into software_versions_yaml
   file "software_versions.csv"

   script:
   """
   echo $workflow.manifest.version > v_pipeline.txt
   echo $workflow.nextflow.version > v_nextflow.txt
   bowtie2 --version > v_bowtie2.txt
   python --version > v_python.txt 2>&1
   samtools --version > v_samtools.txt
   multiqc --version > v_multiqc.txt
   scrape_software_versions.py &> software_versions_mqc.yaml
   """
}

/****************************************************
 * PRE-PROCESSING
 */

if(!params.bwt2_index && params.fasta){
    process makeBowtie2Index {
        tag "$bwt2_base"
        label 'process_highmem'
        publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        file fasta from fasta_for_index

        output:
        file "bowtie2_index" into bwt2_index_end2end
	file "bowtie2_index" into bwt2_index_trim

        script:
        """
        mkdir bowtie2_index
	bowtie2-build ${fasta} bowtie2_index/${bwt2_base}
	"""
      }
 }


if(!params.chromosome_size && params.fasta){
    process makeChromSize {
        tag "$fasta"
	label 'process_low'
        publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        file fasta from fasta_for_chromsize

        output:
        file "*.size" into chrsize, chrsize_build, chrsize_raw, chrsize_balance, chrsize_zoom

        script:
        """
	samtools faidx ${fasta}
	cut -f1,2 ${fasta}.fai > chrom.size
   	"""
      }
 }

if(!params.restriction_fragments && params.fasta && !params.dnase){
    process getRestrictionFragments {
        tag "$fasta ${params.restriction_site}"
	label 'process_low'
        publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        file fasta from fasta_for_resfrag

        output:
        file "*.bed" into res_frag_file

        script:
        """
	digest_genome.py -r ${params.restriction_site} -o restriction_fragments.bed ${fasta}
	"""
      }
 }

/****************************************************
 * MAIN WORKFLOW
 */

/*
 * HiC-pro - Two-steps Reads Mapping
 */

process bowtie2_end_to_end {
   tag "$sample"
   label 'process_medium'
   publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/hicpro/mapping" : params.outdir },
   	      saveAs: { params.save_aligned_intermediates ? it : null }, mode: params.publish_dir_mode

   input:
   set val(sample), file(reads) from raw_reads
   file index from bwt2_index_end2end.collect()

   output:
   set val(sample), file("${prefix}_unmap.fastq") into unmapped_end_to_end
   set val(sample), file("${prefix}.bam") into end_to_end_bam

   script:
   prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
   def bwt2_opts = params.bwt2_opts_end2end
   if (!params.dnase){
   """
   bowtie2 --rg-id BMG --rg SM:${prefix} \\
	${bwt2_opts} \\
	-p ${task.cpus} \\
	-x ${index}/${bwt2_base} \\
	--un ${prefix}_unmap.fastq \\
 	-U ${reads} | samtools view -F 4 -bS - > ${prefix}.bam
   """
   }else{
   """
   bowtie2 --rg-id BMG --rg SM:${prefix} \\
	${bwt2_opts} \\
	-p ${task.cpus} \\
	-x ${index}/${bwt2_base} \\
	--un ${prefix}_unmap.fastq \\
 	-U ${reads} > ${prefix}.bam
   """
   }
}

process trim_reads {
   tag "$sample"
   label 'process_low'
   publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/hicpro/mapping" : params.outdir },
   	      saveAs: { params.save_aligned_intermediates ? it : null }, mode: params.publish_dir_mode

   when:
   !params.dnase

   input:
   set val(sample), file(reads) from unmapped_end_to_end

   output:
   set val(sample), file("${prefix}_trimmed.fastq") into trimmed_reads

   script:
   prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
   """
   cutsite_trimming --fastq $reads \\
                    --cutsite  ${params.ligation_site} \\
                    --out ${prefix}_trimmed.fastq
   """
}

process bowtie2_on_trimmed_reads {
   tag "$sample"
   label 'process_medium'
   publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/hicpro/mapping" : params.outdir },
   	      saveAs: { params.save_aligned_intermediates ? it : null }, mode: params.publish_dir_mode

   when:
   !params.dnase

   input:
   set val(sample), file(reads) from trimmed_reads
   file index from bwt2_index_trim.collect()

   output:
   set val(sample), file("${prefix}_trimmed.bam") into trimmed_bam

   script:
   prefix = reads.toString() - ~/(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
   """
   bowtie2 --rg-id BMG --rg SM:${prefix} \\
           ${params.bwt2_opts_trimmed} \\
           -p ${task.cpus} \\
           -x ${index}/${bwt2_base} \\
           -U ${reads} | samtools view -bS - > ${prefix}_trimmed.bam
   """
}

if (!params.dnase){
   process bowtie2_merge_mapping_steps{
      tag "$prefix = $bam1 + $bam2"
      label 'process_medium'
      publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/hicpro/mapping" : params.outdir },
   	      saveAs: { params.save_aligned_intermediates ? it : null }, mode: params.publish_dir_mode

      input:
      set val(prefix), file(bam1), file(bam2) from end_to_end_bam.join( trimmed_bam ).dump(tag:'merge')

      output:
      set val(sample), file("${prefix}_bwt2merged.bam") into bwt2_merged_bam
      set val(oname), file("${prefix}.mapstat") into all_mapstat

      script:
      //sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)/
      sample = prefix.toString() - ~/(_R1|_R2)/
      //tag = prefix.toString() =~/_R1|_val_1|_1/ ? "R1" : "R2"
      tag = prefix.toString() =~/_R1/ ? "R1" : "R2"
      oname = prefix.toString() - ~/(\.[0-9]+)$/
      """
      samtools merge -@ ${task.cpus} \\
    	             -f ${prefix}_bwt2merged.bam \\
                     ${bam1} ${bam2}

      samtools sort -@ ${task.cpus} -m 800M \\
      	            -n -T /tmp/ \\
	            -o ${prefix}_bwt2merged.sorted.bam \\
	            ${prefix}_bwt2merged.bam

      mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam

      echo "## ${prefix}" > ${prefix}.mapstat
      echo -n "total_${tag}\t" >> ${prefix}.mapstat
      samtools view -c ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
      echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
      echo -n "global_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam1} >> ${prefix}.mapstat
      echo -n "local_${tag}\t"  >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam2} >> ${prefix}.mapstat
      """
   }
}else{
   process dnase_mapping_stats{
      tag "$sample = $bam"
      label 'process_medium'
      publishDir path: { params.save_aligned_intermediates ? "${params.outdir}/hicpro/mapping" : params.outdir },
   	      saveAs: { params.save_aligned_intermediates ? it : null }, mode: params.publish_dir_mode

      input:
      set val(prefix), file(bam) from end_to_end_bam

      output:
      set val(sample), file(bam) into bwt2_merged_bam
      set val(oname), file("${prefix}.mapstat") into all_mapstat

      script:
      //sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)/
      sample = prefix.toString() - ~/(_R1|_R2)/
      //tag = prefix.toString() =~/_R1|_val_1|_1/ ? "R1" : "R2"
      tag = prefix.toString() =~/_R1/ ? "R1" : "R2"
      oname = prefix.toString() - ~/(\.[0-9]+)$/
      """
      echo "## ${prefix}" > ${prefix}.mapstat
      echo -n "total_${tag}\t" >> ${prefix}.mapstat
      samtools view -c ${bam} >> ${prefix}.mapstat
      echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam} >> ${prefix}.mapstat
      echo -n "global_${tag}\t" >> ${prefix}.mapstat
      samtools view -c -F 4 ${bam} >> ${prefix}.mapstat
      echo -n "local_${tag}\t0"  >> ${prefix}.mapstat
      """
   }
}

process combine_mates{
   tag "$sample = $r1_prefix + $r2_prefix"
   label 'process_low'
   publishDir "${params.outdir}/hicpro/mapping", mode: params.publish_dir_mode,
   	      saveAs: {filename -> filename.indexOf(".pairstat") > 0 ? "stats/$filename" : "$filename"}

   input:
   set val(sample), file(aligned_bam) from bwt2_merged_bam.groupTuple().dump(tag:'mates')

   output:
   set val(oname), file("${sample}_bwt2pairs.bam") into paired_bam
   set val(oname), file("*.pairstat") into all_pairstat

   script:
   r1_bam = aligned_bam[0]
   r1_prefix = r1_bam.toString() - ~/_bwt2merged.bam$/
   r2_bam = aligned_bam[1]
   r2_prefix = r2_bam.toString() - ~/_bwt2merged.bam$/
   oname = sample.toString() - ~/(\.[0-9]+)$/

   def opts = "-t"
   if (params.keep_multi) {
     opts="${opts} --multi"
   }else if (params.min_mapq){
     opts="${opts} -q ${params.min_mapq}"
   }
   """
   mergeSAM.py -f ${r1_bam} -r ${r2_bam} -o ${sample}_bwt2pairs.bam ${opts}
   """
}

/*
 * HiC-Pro - detect valid interaction from aligned data
 */

if (!params.dnase){
   process get_valid_interaction{
      tag "$sample"
      label 'process_low'
      publishDir "${params.outdir}/hicpro/valid_pairs", mode: params.publish_dir_mode,
   	      saveAs: {filename -> filename.indexOf(".stat") > 0 ? "stats/$filename" : "$filename"}

      input:
      set val(sample), file(pe_bam) from paired_bam
      file frag_file from res_frag_file.collect()

      output:
      set val(sample), file("*.validPairs") into valid_pairs
      set val(sample), file("*.validPairs") into valid_pairs_4cool
      set val(sample), file("*.DEPairs") into de_pairs
      set val(sample), file("*.SCPairs") into sc_pairs
      set val(sample), file("*.REPairs") into re_pairs
      set val(sample), file("*.FiltPairs") into filt_pairs
      set val(sample), file("*RSstat") into all_rsstat

      script:
      if (params.split_fastq){
         sample = sample.toString() - ~/(\.[0-9]+)$/
      }

      def opts = ""
      opts += params.min_cis_dist > 0 ? " -d ${params.min_cis_dist}" : ''
      opts += params.min_insert_size > 0 ?  " -s ${params.min_insert_size}" : ''
      opts += params.max_insert_size > 0 ? " -l ${params.max_insert_size}" : ''
      opts += params.min_restriction_fragment_size > 0 ? " -t ${params.min_restriction_fragment_size}" : ''
      opts += params.max_restriction_fragment_size > 0 ? " -m ${params.max_restriction_fragment_size}" : ''
      opts += params.save_interaction_bam ? " --sam" : ''
      prefix = pe_bam.toString() - ~/.bam/
      """
      mapped_2hic_fragments.py -f ${frag_file} -r ${pe_bam} --all ${opts}
      sort -T /tmp/ -k2,2V -k3,3n -k5,5V -k6,6n -o ${prefix}.validPairs ${prefix}.validPairs
      """
   }
}
else{
   process get_valid_interaction_dnase{
      tag "$sample"
      label 'process_low'
      publishDir "${params.outdir}/hicpro/valid_pairs", mode: params.publish_dir_mode,
   	      saveAs: {filename -> filename.indexOf(".stat") > 0 ? "stats/$filename" : "$filename"}

      input:
      set val(sample), file(pe_bam) from paired_bam

      output:
      set val(sample), file("*.validPairs") into valid_pairs
      set val(sample), file("*.validPairs") into valid_pairs_4cool
      set val(sample), file("*RSstat") into all_rsstat

      script:
      if (params.split_fastq){
         sample = sample.toString() - ~/(\.[0-9]+)$/
      }

      opts = params.min_cis_dist > 0 ? " -d ${params.min_cis_dist}" : ''
      prefix = pe_bam.toString() - ~/.bam/
      """
      mapped_2hic_dnase.py -r ${pe_bam} ${opts}
      sort -T /tmp/ -k2,2V -k3,3n -k5,5V -k6,6n -o ${prefix}.validPairs ${prefix}.validPairs
      """
   }
}

/*
 * Remove duplicates
 */

process remove_duplicates {
   tag "$sample"
   label 'process_highmem'
   publishDir "${params.outdir}/hicpro/valid_pairs", mode: params.publish_dir_mode,
   	      saveAs: {filename -> filename.indexOf(".stat") > 0 ? "stats/$sample/$filename" : "$filename"}

   input:
   set val(sample), file(vpairs) from valid_pairs.groupTuple().dump(tag:'final')

   output:
   set val(sample), file("*.allValidPairs") into ch_vpairs, ch_vpairs_cool
   file("stats/") into all_mergestat

   script:
   if ( ! params.keep_dups ){
   """
   mkdir -p stats/${sample}

   ## Sort valid pairs and remove read pairs with same starts (i.e duplicated read pairs)
   sort -T /tmp/ -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ${vpairs} | \
   awk -F"\\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=\$2 || c2!=\$5 || s1!=\$3 || s2!=\$6){print;c1=\$2;c2=\$5;s1=\$3;s2=\$6}' > ${sample}.allValidPairs

   echo -n "valid_interaction\t" > stats/${sample}/${sample}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> stats/${sample}/${sample}_allValidPairs.mergestat
   cat ${sample}.allValidPairs | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat

   ## Count short range (<20000) vs long range contacts
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> stats/${sample}/${sample}_allValidPairs.mergestat

   """
   }else{
   """
   mkdir -p stats/${sample}
   cat ${vpairs} > ${sample}.allValidPairs
   echo -n "valid_interaction\t" > stats/${sample}/${sample}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> stats/${sample}/${sample}_allValidPairs.mergestat
   cat ${sample}.allValidPairs | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat

   ## Count short range (<20000) vs long range contacts
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> stats/${sample}/${sample}_allValidPairs.mergestat
   """
   }
}

process merge_stats {
   tag "$ext"
   label 'process_low'
   publishDir "${params.outdir}/hicpro/stats/${sample}", mode: params.publish_dir_mode

   input:
   set val(prefix), file(fstat) from all_mapstat.groupTuple().concat(all_pairstat.groupTuple(), all_rsstat.groupTuple())

   output:
   file("mstats/") into all_mstats

  script:
  sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2|_1|_2)/
  if ( (fstat =~ /.mapstat/) ){ ext = "mmapstat" }
  if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
  if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }
  """
  mkdir -p mstats/${sample}
  merge_statfiles.py -f ${fstat} > mstats/${sample}/${prefix}.${ext}
  """
}

/*
 * HiC-Pro build matrix processes
 * ONGOING VALIDATION - TO REPLACED BY COOLER ?
 */


process build_contact_maps{
   tag "$sample - $mres"
   label 'process_highmem'
   publishDir "${params.outdir}/hicpro/matrix/raw", mode: params.publish_dir_mode

   when:
   !params.skip_maps

   input:
   set val(sample), file(vpairs), val(mres) from ch_vpairs.combine(map_res)
   file chrsize from chrsize.collect()

   output:
   set val(sample), val(mres), file("*.matrix"), file("*.bed") into raw_maps, raw_maps_4cool
   
   script:
   """
   build_matrix --matrix-format upper  --binsize ${mres} --chrsizes ${chrsize} --ifile ${vpairs} --oprefix ${sample}_${mres}
   """
}

process run_ice{
   tag "$rmaps"
   label 'process_highmem'
   publishDir "${params.outdir}/hicpro/matrix/iced", mode: params.publish_dir_mode

   when:
   !params.skip_maps && !params.skip_ice

   input:
   set val(sample), val(res), file(rmaps), file(bed) from raw_maps

   output:
   set val(sample), val(res), file("*iced.matrix"), file(bed) into iced_maps_4h5, iced_maps_4cool
   file ("*.biases") into iced_bias

   script:
   prefix = rmaps.toString() - ~/(\.matrix)?$/
   """
   ice --filter_low_counts_perc ${params.ice_filer_low_count_perc} \
   --results_filename ${prefix}_iced.matrix \
   --filter_high_counts_perc ${params.ice_filer_high_count_perc} \
   --max_iter ${params.ice_max_iter} --eps ${params.ice_eps} --remove-all-zeros-loci --output-bias 1 --verbose 1 ${rmaps}
   """
}


/*
 * Cooler
 */

process convert_to_pairs {
   tag "$sample"
   label 'process_medium'

   when:
   !params.skip_maps

   input:
   set val(sample), file(vpairs) from ch_vpairs_cool
   file chrsize from chrsize_build.collect()

   output:
   set val(sample), file("*.txt.gz") into cool_build, cool_build_zoom

   script:
   """
   ## chr/pos/strand/chr/pos/strand
   awk '{OFS="\t";print \$1,\$2,\$3,\$5,\$6,\$4,\$7}' $vpairs > contacts.txt
   gzip contacts.txt
   """
}


process cooler_raw {
  tag "$sample - ${res}"
  label 'process_medium'

  publishDir "${params.outdir}/contact_maps/", mode: 'copy',
              saveAs: {filename -> filename.indexOf(".cool") > 0 ? "raw/cool/$filename" : "raw/txt/$filename"}

  input:
  set val(sample), file(contacts), val(res) from cool_build.combine(map_res_cool)
  file chrsize from chrsize_raw.collect()

  output:
  set val(sample), val(res), file("*cool") into raw_cool_maps
  set file("*.bed"), file("${sample}_${res}.txt") into raw_txt_maps

  script:
  """
  cooler makebins ${chrsize} ${res} > ${sample}_${res}.bed
  cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${sample}_${res}.bed ${contacts} ${sample}_${res}.cool
  cooler dump ${sample}_${res}.cool | awk '{OFS="\t"; print \$1+1,\$2+1,\$3}' > ${sample}_${res}.txt
  """
}

process cooler_balance {
  tag "$sample - ${res}"
  label 'process_medium'

  publishDir "${params.outdir}/contact_maps/", mode: 'copy',
              saveAs: {filename -> filename.indexOf(".cool") > 0 ? "norm/cool/$filename" : "norm/txt/$filename"}

  when:
  !params.skip_balancing

  input:
  set val(sample), val(res), file(cool) from raw_cool_maps
  file chrsize from chrsize_balance.collect()

  output:
  set val(sample), val(res), file("${sample}_${res}_norm.cool") into norm_cool_maps, norm_cool_maps_h5
  file("${sample}_${res}_norm.txt") into norm_txt_maps

  script:
  """
  cp ${cool} ${sample}_${res}_norm.cool
  cooler balance ${sample}_${res}_norm.cool -p ${task.cpus} --force
  cooler dump ${sample}_${res}_norm.cool --balanced --na-rep 0 | awk '{OFS="\t"; print \$1+1,\$2+1,\$4}' > ${sample}_${res}_norm.txt
  """
}

process cooler_zoomify {
   tag "$sample"
   label 'process_medium'
   publishDir "${params.outdir}/contact_maps/norm/mcool", mode: 'copy'

   when:
   !params.skip_mcool

   input:
   set val(sample), file(contacts)  from cool_build_zoom
   file chrsize from chrsize_zoom.collect()

   output:
   file("*mcool") into mcool_maps

   script:
   """
   cooler makebins ${chrsize} ${params.res_zoomify} > bins.bed
   cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 bins.bed ${contacts} ${sample}.cool
   cooler zoomify --nproc ${task.cpus} --balance ${sample}.cool
   """
}


/*
 * Create h5 file
 */

process convert_to_h5 {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/contact_maps/norm/h5", mode: 'copy'

  input:
  set val(sample), val(res), file(maps)  from norm_cool_maps_h5

  output:
  set val(sample), val(res), file("*.h5") into h5maps_ddecay, h5maps_ccomp, h5maps_tads

  script:
  """
  hicConvertFormat --matrices ${maps} \
  		   --outFileName ${maps.baseName}.h5 \
		   --resolution ${res} \
		   --inputFormat cool \
		   --outputFormat h5 \
  """
}


/****************************************************
 * DOWNSTREAM ANALYSIS
 */

/*
 * Counts vs distance QC
 */

if (!params.skip_dist_decay){
  chddecay = h5maps_ddecay.combine(ddecay_res).filter{ it[1] == it[3] }.dump(tag: "ddecay") 
}else{
  chddecay = Channel.empty()
}

process dist_decay {
  tag "$sample"
  label 'process_medium'
  publishDir "${params.outdir}/dist_decay", mode: 'copy'

  when:
  !params.skip_dist_decay

  input:
  set val(sample), val(res), file(h5mat), val(r) from chddecay
  
  output:
  file("*_distcount.txt")
  file("*.png")


  script:
  """
  hicPlotDistVsCounts --matrices ${h5mat} \
                      --plotFile ${h5mat.baseName}_distcount.png \
  		      --outFileData ${h5mat.baseName}_distcount.txt
  """
}

/*
 * Compartment calling
 */

/*
if(!params.skip_compartments){
  chcomp = iced_maps_comp.combine(comp_res).filter{ it[1] == it[4] }.dump(tag: "comp")
}else{
  chcomp = Channel.empty()
}

process compartment_calling {
  tag "$sample - $res"
  label 'process_medium'
  publishDir "${params.outdir}/compartments", mode: 'copy'

  when:
  !params.skip_compartments

  input:
  set val(sample), val(res), file(mat), file(bed), val(r) from chcomp

  output:
  file("*.bedgraph") optional true into out_compartments

  script:
  """
  call_compartments.r --matrix ${mat} --bed ${bed}
  """
}
*/


/*
 * TADs calling
 */

if (!params.skip_tads){
  chtads = h5maps_tads.combine(tads_res_hicexplorer).filter{ it[1] == it[3] }.dump(tag: "hicexp")
}else{
  chtads = Channel.empty()
}

process tads_hicexplorer {
  tag "$sample - $res"
  label 'process_medium'
  publishDir "${params.outdir}/tads", mode: 'copy'

  when:
  !params.skip_tads && params.tads_caller =~ 'hicexplorer'

  input:
  set val(sample), val(res), file(h5mat), val(r) from chtads

  output:
  file("*.{bed,bedgraph,gff}") into hicexplorer_tads

  script:
  """
  hicFindTADs --matrix ${h5mat} \
  	      --outPrefix tad \
	      --correctForMultipleTesting fdr \
	      --numberOfProcessors ${task.cpus}
  """
}

if (!params.skip_tads){
  chIS = norm_cool_maps.combine(tads_res_insulation).filter{ it[1] == it[3] }.dump(tag : "ins")
}else{
  chIS = Channel.empty()
}

process tads_insulation {
  tag "$sample - $res"
  label 'process_medium'
  publishDir "${params.outdir}/tads", mode: 'copy'

  when:
  !params.skip_tads && params.tads_caller =~ 'insulation'

  input:
  set val(sample), val(res), file(cool), val(r) from chIS

  output:
  file("*tsv") into insulation_tads

  script:
  """
  cooltools diamond-insulation --window-pixels ${cool} 15 25 50 > ${sample}_insulation.tsv
  """
}


/*
 * MultiQC
 */

process multiqc {
   label 'process_low'
   publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

   when:
   !params.skip_multiqc

   input:
   file multiqc_config from ch_multiqc_config
   file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
   file ('input_*/*') from all_mstats.concat(all_mergestat).collect()
   file ('software_versions/*') from software_versions_yaml
   file workflow_summary from ch_workflow_summary.collect()

   output:
   file "*multiqc_report.html" into multiqc_report
   file "*_data"

   script:
   rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
   rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
   """
   multiqc -f $rtitle $rfilename --config $multiqc_config .
   """
}


/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */

workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/hic] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/hic] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/hic] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/hic] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/hic] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/hic] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/hic]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/hic]${c_red} Pipeline completed with errors${c_reset}-"
    }
}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/hic v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
