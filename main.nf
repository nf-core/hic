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

    nextflow run nf-core/hic --reads '*_R{1,2}.fastq.gz' -profile conda

    Mandatory arguments:
      --reads				    Path to input data (must be surrounded with quotes)
      -profile                      	    Configuration profile to use. Can use multiple (comma separated)
                                    	    Available: conda, docker, singularity, awsbatch, test and more.

    References                      	    If not specified in the configuration file or you wish to overwrite any of the references.
      --genome                              Name of iGenomes reference
      --bwt2_index                     	    Path to Bowtie2 index
      --fasta                       	    Path to Fasta reference
      --chromosome_size             	    Path to chromosome size file
      --restriction_fragments    	    Path to restriction fragment file (bed)

    Options:
      --bwt2_opts_end2end		    Options for bowtie2 end-to-end mappinf (first mapping step)
      --bwt2_opts_trimmed	    	    Options for bowtie2 mapping after ligation site trimming
      --min_mapq		    	    Minimum mapping quality values to consider

      --restriction_site	    	    Cutting motif(s) of restriction enzyme(s) (comma separated)
      --ligation_site		    	    Ligation motifs to trim (comma separated)

      --min_restriction_fragment_size	    Minimum size of restriction fragments to consider
      --max_restriction_framgnet_size	    Maximum size of restriction fragmants to consider
      --min_insert_size			    Minimum insert size of mapped reads to consider
      --max_insert_size			    Maximum insert size of mapped reads to consider
      --min_cis_dist			    Minimum intra-chromosomal distance to consider
      --rm_singleton			    Remove singleton reads
      --rm_multi			    Remove multi-mapped reads
      --rm_dup				    Remove duplicates

      --bin_size			    Bin size for contact maps (comma separated)
      --ice_max_iter			    Maximum number of iteration for ICE normalization
      --ice_filter_low_count_perc	    Percentage of low counts columns/rows to filter before ICE normalization
      --ice_filter_high_count_perc	    Percentage of high counts columns/rows to filter before ICE normalization
      --ice_eps				    Convergence criteria for ICE normalization

    Other options:
      --splitFastq			    Size of read chuncks to use to speed up the workflow
      --outdir				    The output directory where the results will be saved
      --email                       	    Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         	    Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    Step options:
      --skip_cool			    Skip generation of cool files
      --skip_multiQC			    Skip MultiQC

    AWSBatch options:
      --awsqueue			    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   	    The AWS Region for your AWS Batch job to run on
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

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Reference index path configuration
params.bwt2_index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false 
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false


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

if (params.readPaths){

   raw_reads = Channel.create()
   raw_reads_2 = Channel.create()

   Channel
      .from( params.readPaths )
      .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
      .separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0], a[1][0]), tuple(a[0], a[1][1])] }
 }else{

   raw_reads = Channel.create()
   raw_reads_2 = Channel.create()

   Channel
      .fromFilePairs( params.reads )
      .separate( raw_reads, raw_reads_2 ) { a -> [tuple(a[0], a[1][0]), tuple(a[0], a[1][1])] }
}

if ( params.splitFastq ){
   raw_reads_full = raw_reads.concat( raw_reads_2 )
   raw_reads = raw_reads_full.splitFastq( by: params.splitFastq , file: true)
 }else{
   raw_reads = raw_reads.concat( raw_reads_2 )
}


// SPlit fastq files
// https://www.nextflow.io/docs/latest/operator.html#splitfastq

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
    bwt2_base = params.fasta.substring(lastPath+1)

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
         .into {chromosome_size; chromosome_size_cool}
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
else {
    exit 1, "No restriction fragments file specified!"
}

// Resolutions for contact maps
map_res = Channel.from( params.bins_size.tokenize(',') )

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

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

summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta


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
 * PRE-PROCESSING
 */

if(!params.bwt2_index && params.fasta){
    process makeBowtie2Index {
        tag "$bwt2_base"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_index

        output:
        file "bowtie2_index" into bwt2_index_end2end
	file "bowtie2_index" into bwt2_index_trim

        script:
        bwt2_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
        """
        mkdir bowtie2_index
	bowtie2-build ${fasta} bowtie2_index/${bwt2_base}
	"""
      }
 }


if(!params.chromosome_size && params.fasta){
    process makeChromSize {
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_chromsize

        output:
        file "*.size" into chromosome_size, chromosome_size_cool 

        script:
        """
	samtools faidx ${fasta}
	cut -f1,2 ${fasta}.fai > chrom.size
   	"""	
      }
 }

if(!params.restriction_fragments && params.fasta){
    process getRestrictionFragments {
        tag "$fasta [${params.restriction_site}]"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

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
 * STEP 1 - Two-steps Reads Mapping
*/

process bowtie2_end_to_end {
   tag "$prefix"
   publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/mapping" : params.outdir },
   	      saveAs: { params.saveAlignedIntermediates ? it : null }, mode: 'copy'

   input:
        set val(sample), file(reads) from raw_reads
        file index from bwt2_index_end2end.collect()
 
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
		-x ${index}/${bwt2_base} \\
		--un ${prefix}_unmap.fastq \\
	 	-U ${reads} | samtools view -F 4 -bS - > ${prefix}.bam
        """
}

process trim_reads {
   tag "$prefix"
   publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/mapping" : params.outdir },
   	      saveAs: { params.saveAlignedIntermediates ? it : null }, mode: 'copy'

   input:
      set val(prefix), file(reads) from unmapped_end_to_end

   output:
      set val(prefix), file("${prefix}_trimmed.fastq") into trimmed_reads

   script:
      """
      cutsite_trimming --fastq $reads \\
       		       --cutsite  ${params.ligation_site} \\
                       --out ${prefix}_trimmed.fastq
      """
}

process bowtie2_on_trimmed_reads {
   tag "$prefix"
   publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/mapping" : params.outdir },
   	      saveAs: { params.saveAlignedIntermediates ? it : null }, mode: 'copy'

   input:
      set val(prefix), file(reads) from trimmed_reads
      file index from bwt2_index_trim.collect()

   output:
      set val(prefix), file("${prefix}_trimmed.bam") into trimmed_bam

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

process merge_mapping_steps{
   tag "$sample = $bam1 + $bam2"
   publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/mapping" : params.outdir },
   	      saveAs: { params.saveAlignedIntermediates ? it : null }, mode: 'copy'

   input:
      set val(prefix), file(bam1), file(bam2) from end_to_end_bam.join( trimmed_bam )

   output:
      set val(sample), file("${prefix}_bwt2merged.bam") into bwt2_merged_bam
      set val(oname), file("${prefix}.mapstat") into all_mapstat

   script:
      sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2)/
      tag = prefix.toString() =~/_R1|_val_1/ ? "R1" : "R2"
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

process combine_mapped_files{
   tag "$sample = $r1_prefix + $r2_prefix"
   publishDir "${params.outdir}/mapping", mode: 'copy',
   	      saveAs: {filename -> filename.indexOf(".pairstat") > 0 ? "stats/$filename" : "$filename"}	      
 
   input:
      set val(sample), file(aligned_bam) from bwt2_merged_bam.groupTuple()

   output:
      set val(sample), file("${sample}_bwt2pairs.bam") into paired_bam
      set val(oname), file("*.pairstat") into all_pairstat

   script:
      r1_bam = aligned_bam[0]
      r1_prefix = r1_bam.toString() - ~/_bwt2merged.bam$/
      r2_bam = aligned_bam[1]
      r2_prefix = r2_bam.toString() - ~/_bwt2merged.bam$/
      oname = sample.toString() - ~/(\.[0-9]+)$/
 
      def opts = "-t"
      opts = params.rm_singleton ? "${opts}" : "--single ${opts}"
      opts = params.rm_multi ? "${opts}" : "--multi ${opts}"
      if ("$params.min_mapq".isInteger()) opts="${opts} -q ${params.min_mapq}"
      """
      mergeSAM.py -f ${r1_bam} -r ${r2_bam} -o ${sample}_bwt2pairs.bam ${opts}
      """
}


/*
 * STEP2 - DETECT VALID PAIRS
*/


process get_valid_interaction{
   tag "$sample"
   publishDir "${params.outdir}/hic_results/data", mode: 'copy',
   	      saveAs: {filename -> filename.indexOf("*stat") > 0 ? "stats/$filename" : "$filename"}	      

   input:
      set val(sample), file(pe_bam) from paired_bam
      file frag_file from res_frag_file.collect()

   output:
      set val(sample), file("*.validPairs") into valid_pairs
      set val(sample), file("*.validPairs") into valid_pairs_4cool
      set val(sample), file("*RSstat") into all_rsstat

   script:
	
      if (params.splitFastq){
      	 sample = sample.toString() - ~/(\.[0-9]+)$/
      }

      def opts = ""
      if ("$params.min_cis_dist".isInteger()) opts="${opts} -d ${params.min_cis_dist}"
      if ("$params.min_insert_size".isInteger()) opts="${opts} -s ${params.min_insert_size}"
      if ("$params.max_insert_size".isInteger()) opts="${opts} -l ${params.max_insert_size}"
      if ("$params.min_restriction_fragment_size".isInteger()) opts="${opts} -t ${params.min_restriction_fragment_size}"
      if ("$params.max_restriction_fragment_size".isInteger()) opts="${opts} -m ${params.max_restriction_fragment_size}"

      """
      mapped_2hic_fragments.py -f ${frag_file} -r ${pe_bam} ${opts}
      """
}

/*
 * STEP3 - BUILD MATRIX
*/

process remove_duplicates {
   tag "$sample"
   publishDir "${params.outdir}/hic_results/data", mode: 'copy',
   	      saveAs: {filename -> filename.indexOf("*stat") > 0 ? "stats/$sample/$filename" : "$filename"}	      

   input:
     set val(sample), file(vpairs) from valid_pairs.groupTuple()

   output:
     set val(sample), file("*.allValidPairs") into all_valid_pairs
     set val(sample), file("*.allValidPairs") into all_valid_pairs_4cool
     file("stats/") into all_mergestat

   script:
   if ( params.rm_dup ){
   """
   mkdir -p stats/${sample}
   sort -T /tmp/ -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ${vpairs} | \
   awk -F"\\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=\$2 || c2!=\$5 || s1!=\$3 || s2!=\$6){print;c1=\$2;c2=\$5;s1=\$3;s2=\$6}' > ${sample}.allValidPairs                   
   echo -n "valid_interaction\t" > stats/${sample}/${sample}_allValidPairs.mergestat
   cat ${vpairs} | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat
   echo -n "valid_interaction_rmdup\t" >> stats/${sample}/${sample}_allValidPairs.mergestat
   cat ${sample}.allValidPairs | wc -l >> stats/${sample}/${sample}_allValidPairs.mergestat
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
   awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} \$2 == \$5{cis=cis+1; d=\$6>\$3?\$6-\$3:\$3-\$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} \$2!=\$5{trans=trans+1}END{print "trans_interaction\\t"trans"\\ncis_interaction\\t"cis"\\ncis_shortRange\\t"sr"\\ncis_longRange\\t"lr}' ${sample}.allValidPairs >> stats/${sample}/${sample}_allValidPairs.mergestat
   """
   }
}

process merge_sample {
   tag "$ext"
   publishDir "${params.outdir}/hic_results/stats/${sample}", mode: 'copy'

   input:
     set val(prefix), file(fstat) from all_mapstat.groupTuple().concat(all_pairstat.groupTuple(), all_rsstat.groupTuple())

  output:
     file("mstats/") into all_mstats

  script:
     sample = prefix.toString() - ~/(_R1|_R2|_val_1|_val_2)/
     if ( (fstat =~ /.mapstat/) ){ ext = "mmapstat" }
     if ( (fstat =~ /.pairstat/) ){ ext = "mpairstat" }
     if ( (fstat =~ /.RSstat/) ){ ext = "mRSstat" }

     """
     mkdir -p mstats/${sample}
     merge_statfiles.py -f ${fstat} > mstats/${sample}/${prefix}.${ext}
     """
}


process build_contact_maps{
   tag "$sample - $mres"
   publishDir "${params.outdir}/hic_results/matrix/raw", mode: 'copy'

   input:
      set val(sample), file(vpairs), val(mres) from all_valid_pairs.combine(map_res)
      file chrsize from chromosome_size.collect()

   output:
      file("*.matrix") into raw_maps
      file "*.bed"

   script:
   """
   build_matrix --matrix-format upper  --binsize ${mres} --chrsizes ${chrsize} --ifile ${vpairs} --oprefix ${sample}_${mres}
   """
}

/*
 * STEP 4 - NORMALIZE MATRIX
*/

process run_ice{
   tag "$rmaps"
   publishDir "${params.outdir}/hic_results/matrix/iced", mode: 'copy'

   input:
      file(rmaps) from raw_maps
      file "*.biases"

   output:
      file("*iced.matrix") into iced_maps

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
 * STEP 5 - COOLER FILE
 */
process generate_cool{
   tag "$sample"
   publishDir "${params.outdir}/export/cool", mode: 'copy'

   when:
      !params.skip_cool

   input:
      set val(sample), file(vpairs) from all_valid_pairs_4cool
      file chrsize from chromosome_size_cool.collect()

   output:
      file("*mcool") into cool_maps

   script:
   """
   hicpro2higlass.sh -i $vpairs -r 5000 -c ${chrsize} -n
   """ 
}


/*
 * STEP 5 - MultiQC
 */ 
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
       !params.skip_multiqc

    input:
       file multiqc_config from ch_multiqc_config
       file ('input_*/*') from all_mstats.concat(all_mergestat).collect()
       file ('software_versions/*') from software_versions_yaml
       file workflow_summary from create_workflow_summary(summary)

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


/****************************************************
 * POST-PROCESSING
 */

/* 
 * Output Description HTML
 */

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
