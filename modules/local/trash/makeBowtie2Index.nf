// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MAKE_BOWTIE2_INDEX {
    tag "$fasta_base"
    label 'process_highmem'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bowtie2=2.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie2:2.3.5--py37he860b03_0"
    } else {
        container "quay.io/biocontainers/bowtie2:2.3.5--py37he860b03_0"
    }

    input:
        path fasta

        output:
        path "bowtie2_index", emit: bwt2_index_end2end
	    path "bowtie2_index", emit: bwt2_index_trim

        script:
        fasta_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?(\.fsa)?$/
        """
        mkdir bowtie2_index
	bowtie2-build ${fasta} bowtie2_index/${fasta_base}
	"""
}
