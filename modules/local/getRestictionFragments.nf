// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process getRestrictionFragments {
        tag "$fasta ${params.restriction_site}"
	label 'process_low'
        publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path fasta 

        output:
        path "*.bed", emit:res_frag_file

        script:
        """
	digest_genome.py -r ${params.restriction_site} -o restriction_fragments.bed ${fasta}
	"""
}
