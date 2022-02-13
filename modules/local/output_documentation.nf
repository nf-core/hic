// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path output_docs 
    path images 

    output:
    path 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}
