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
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/hic v${workflow.manifest.version}
    =======================================================

    This pipeline is a Nextflow version of the HiC-Pro pipeline for Hi-C data processing.
    See https://github.com/nservant/HiC-Pro for details.

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/hic --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --readsPath                   Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:


    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/**********************************************************
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
//fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
//if ( params.fasta ){
//    fasta = file(params.fasta)
//    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
//}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")




/**********************************************************
 * SET UP CHANNELS
 */

/*
 * input read files
 */
Channel
        .fromFilePairs( params.readPaths )
	.ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .set { raw_reads_pairs }

raw_reads = Channel.create()
raw_reads_2 = Channel.create()
Channel
        .fromFilePairs( params.readPaths )
        .separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0], a[1][0]), tuple(a[0], a[1][1])] }


// SPlit fastq files
// https://www.nextflow.io/docs/latest/operator.html#splitfastq

/*
 * Other input channels
 */

// Bowtie2 Index
bwt2_file = file("${params.bwt2_index}.1.bt2")
if( !bwt2_file.exists() ) exit 1, "Reference genome Bowtie 2 not found: ${params.bwt2_index}"
bwt2_index = Channel.value( "${params.bwt2_index}" )

// Restriction fragment
res_frag_file = Channel.value( "${params.restriction_fragment_bed}" )

// Chromosome size
chr_size = Channel.value( "${params.chromosome_size}" )



/**********************************************************
 * SET UP LOGS
 */

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/hic v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/hic'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
//summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-hic-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/hic Workflow Summary'
    section_href: 'https://github.com/nf-core/hic'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
     """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    bowtie2 --version > v_bowtie2.txt
    python --version > v_python.txt
    samtools --version > v_samtools.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/****************************************************
 * MAIN WORKFLOW
 */

/*
 * STEP 1 - Two-steps Reads Mapping
*/

raw_reads = raw_reads.concat( raw_reads_2 )

process bowtie2_end_to_end {
   tag "$prefix"
   input:
        set val(sample), file(reads) from raw_reads
        val bt2_index from bwt2_index
 
   output:
	set val(prefix), file("${prefix}_unmap.fastq") into unmapped_end_to_end
     	set val(prefix), file("${prefix}.bam") into end_to_end_bam

   script:
	prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
        def bwt2_opts = params.bwt2_opts_end2end
	"""
        bowtie2 --rg-id BMG --rg SM:${prefix} \\
		${bwt2_opts} \\
		-p ${task.cpus} \\
		-x ${bt2_index} \\
		--un ${prefix}_unmap.fastq \\
	 	-U ${reads} | samtools view -F 4 -bS - > ${prefix}.bam
        """
}

process trim_reads {
   tag "$prefix"
   input:
      set val(prefix), file(reads) from unmapped_end_to_end

   output:
      set val(prefix), file("${prefix}_trimmed.fastq") into trimmed_reads

   script:
      """
      cutsite_trimming --fastq $reads \\
       		       --cutsite  params.ligation_motifs \\
                       --out ${prefix}_trimmed.fastq
      """
}

process bowtie2_on_trimmed_reads {
   tag "$prefix"
   input:
      set val(prefix), file(reads) from trimmed_reads
      val bt2_index from bwt2_index

   output:
      set val(prefix), file("${prefix}_trimmed.bam") into trimmed_bam

   script:
      prefix = reads.toString() - ~/(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
      def bwt2_opts = params.bwt2_opts_trimmed
      """
      bowtie2 --rg-id BMG --rg SM:${prefix} \\
      	      ${bwt2_opts} \\
              -p ${task.cpus} \\
	      -x ${bt2_index} \\
	      -U ${reads} | samtools view -bS - > ${prefix}_trimmed.bam
      """
}

process merge_mapping_steps{
   tag "$bam1 + $bam2"
   input:
      set val(prefix), file(bam1), file(bam2) from end_to_end_bam.join( trimmed_bam )

   output:
      set val(sample), file("${prefix}_bwt2merged.bam") into bwt2_merged_bam

   script:
      sample = prefix.toString() - ~/(_R1)?(_R2)?(_val_1)?(_val_2)?$/
      """
      samtools merge -@ ${task.cpus} \\
       	             -f ${prefix}_bwt2merged.bam \\
	             ${bam1} ${bam2} 

      samtools sort -@ ${task.cpus} -m 800M \\
      	            -n -T /tmp/ \\
	            -o ${prefix}_bwt2merged.sorted.bam \\
	            ${prefix}_bwt2merged.bam
            
      mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam
      """
}


process combine_mapped_files{
   tag "$sample = $r1_prefix + $r2_prefix"
   input:
      set val(sample), file(aligned_bam) from bwt2_merged_bam.groupTuple()

   output:
      set val(sample), file("${sample}_bwt2pairs.bam") into paired_bam

   script:
      r1_bam = aligned_bam[0]
      r1_prefix = r1_bam.toString() - ~/_bwt2merged.bam$/
      r2_bam = aligned_bam[1]
      r2_prefix = r2_bam.toString() - ~/_bwt2merged.bam$/
      """
      mergeSAM.py -f ${r1_bam} -r ${r2_bam} -o ${sample}_bwt2pairs.bam
      """
}

/*
 * STEP2 - DETECT VALID PAIRS
*/

process get_valid_interaction{
   tag "$sample"
   input:
      set val(sample), file(pe_bam) from paired_bam
      val frag_file from res_frag_file

   output:
      set val(sample), file("*.validPairs") into valid_pairs

   script:
      """
      mapped_2hic_fragments.py -f ${frag_file} -r ${pe_bam}
      """
}


/*
 * STEP3 - BUILD MATRIX
*/

process build_contact_maps{
   tag "$sample"
   input:
      set val(sample), file(vpairs) from valid_pairs
      val chrsize from chr_size

   output:
      set val(sample), file("*.matrix") into matrix_file

   script:
   """
   build_matrix --matrix-format upper  --binsize 1000000 --chrsizes ${chrsize} --ifile ${vpairs} --oprefix ${sample}_1000000
   """

}


/*
 // STEP 2 - MultiQC

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


// STEP 3 - Output Description HTML

process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}
*/


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/hic] Successful: $workflow.runName"
    if(!workflow.success){
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
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/hic] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/hic] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/hic] Pipeline Complete"

}
